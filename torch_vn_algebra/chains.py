"""
Spin-1/2 chains with a conserved number of up spins (U(1) symmetry, total S^z).

The Hilbert space (C^2)^{(x) L} splits into sectors N = 0, ..., L (number of up spins) of
dimension binom(L, N); every operator that conserves S^z_total is an element of

    M = (+)_{N=0}^{L} M_{binom(L, N)}(C),

the algebra obtained by fusing L single-site algebras M_1 (+) M_1 with charges 0, 1
(composite.fused_tensor_product); the basis order of SpinChain coincides with that construction.
Basis states are bit strings: bit i of a label is the spin on site i (1 = up).

    chain = SpinChain(L=10, boundary='periodic')
    H = chain.xxz(J=1.0, Delta=1.0, h=fields)          # fields: (L,) or (batch, L) -> batch of H
    H5 = chain.xxz(..., sector=5)                      # one sector only: no padding
    rho_A = chain.reduced_state(rho, n_left=4)          # state of sites 0..3 (a SpinChain(4) state)

Hamiltonian:  H = sum_<ij> [ J/2 (S+_i S-_j + S-_i S+_j) + J Delta S^z_i S^z_j ] + sum_i h_i S^z_i.
"""
import math
from typing import List, Optional, Sequence, Union

import torch

from .algebra import TypeIAlgebra
from . import cost
from .channels import InterSectorChannel
from .states import DensityMatrix

Operator = TypeIAlgebra.Operator


def _labels(L: int) -> List[List[int]]:
    """Basis labels per sector, in the order of iterated fused products with single sites."""
    labels = [[0], [1]]                                  # one site: N = 0 -> |0>, N = 1 -> |1>
    for site in range(1, L):
        new = []
        for N in range(site + 2):
            # pairs (c, d) with c + d = N in the order (c, d) ascending: (N - 1, 1), then (N, 0)
            sec = []
            if 0 <= N - 1 < len(labels):
                sec += [lab | (1 << site) for lab in labels[N - 1]]
            if N < len(labels):
                sec += labels[N]
            new.append(sec)
        labels = new
    return labels


class SpinChain:
    def __init__(self, L: int, boundary: str = 'open', complex_valued: bool = True,
                 precision: str = 'double', device=None):
        if boundary not in ('open', 'periodic'):
            raise ValueError("boundary must be 'open' or 'periodic'")
        if L > 20:
            raise ValueError("L > 20 is beyond dense exact diagonalisation")
        self.L, self.boundary = L, boundary
        self.complex_valued, self.precision = complex_valued, precision
        self._algebra = None
        self._device_arg = device
        self._sector_algebras = {}
        ref = self.sector_algebra(0)                        # 1 x 1: fixes dtype and device
        self.device = ref.hilbert.device
        self.dtype = ref.hilbert.dtype
        self.real_dtype = ref.hilbert.real_dtype
        self.labels = [torch.tensor(lab, dtype=torch.int64, device=self.device) for lab in _labels(L)]
        index = torch.empty(2 ** L, dtype=torch.int64, device=self.device)
        for lab in self.labels:
            index[lab] = torch.arange(len(lab), device=self.device)
        self.index_of = index                                  # position of a label in its sector

    @property
    def algebra(self) -> TypeIAlgebra:
        """The full algebra (+)_N M_{binom(L, N)} (created on first use; padded to the largest sector)."""
        if self._algebra is None:
            dims = [math.comb(self.L, N) for N in range(self.L + 1)]
            self._algebra = TypeIAlgebra(dims, dims, complex_valued=self.complex_valued, precision=self.precision,
                                         device=self.device, charges=list(range(self.L + 1)))
        return self._algebra
    # ------------------------------------------------------------------
    @property
    def bonds(self) -> List[tuple]:
        b = [(i, i + 1) for i in range(self.L - 1)]
        if self.boundary == 'periodic' and self.L > 2:
            b.append((self.L - 1, 0))
        return b

    def sector_dim(self, N: int) -> int:
        return math.comb(self.L, N)

    def sector_algebra(self, N: int) -> TypeIAlgebra:
        """The single-sector algebra M_{binom(L, N)} (no padding), cached."""
        if N not in self._sector_algebras:
            d = self.sector_dim(N)
            self._sector_algebras[N] = TypeIAlgebra([d], [d], complex_valued=self.complex_valued,
                                                    precision=self.precision,
                                                    device=getattr(self, 'device', self._device_arg), charges=[N])
        return self._sector_algebras[N]

    def bits(self, N: int) -> torch.Tensor:
        """(binom(L, N), L) tensor of 0/1 spins of the basis states of sector N."""
        return (self.labels[N][:, None] >> torch.arange(self.L, device=self.device)) & 1

    # ------------------------------------------------------------------
    def _diag_and_hops(self, N, J, Delta, h, batch):
        lab = self.labels[N]
        k = len(lab)
        sz = self.bits(N).to(self.real_dtype) - 0.5                       # (k, L)
        diag = torch.zeros(batch, k, dtype=self.real_dtype, device=self.device)
        for i, j in self.bonds:
            diag = diag + J * Delta * (sz[:, i] * sz[:, j])
        if h is not None:
            diag = diag + h @ sz.T                                         # (batch, k)
        rows, cols = [], []
        for i, j in self.bonds:
            src = torch.nonzero(((lab >> i) & 1) != ((lab >> j) & 1)).squeeze(-1)
            dst = self.index_of[lab[src] ^ ((1 << i) | (1 << j))]
            rows.append(dst)
            cols.append(src)
        return diag, (torch.cat(rows) if rows else None), (torch.cat(cols) if cols else None)

    def _block(self, N, J, Delta, h, batch):
        k = self.sector_dim(N)
        cost.check_memory(cost.tensor_bytes((batch, k, k), self.dtype), self.device,
                          f"XXZ block N={N} ({k} x {k}, batch {batch})")
        diag, rows, cols = self._diag_and_hops(N, J, Delta, h, batch)
        M = torch.diag_embed(diag.to(self.dtype))
        if rows is not None and len(rows):
            hop = torch.zeros(k, k, dtype=self.dtype, device=self.device)
            hop.index_put_((rows, cols), torch.full_like(rows, 1, dtype=self.dtype) * (J / 2), accumulate=True)
            M = M + hop
        return M

    def _fields(self, h):
        if h is None:
            return None, 1
        h = torch.as_tensor(h, dtype=self.real_dtype, device=self.device)
        if h.dim() == 0:
            h = h.expand(self.L)
        if h.dim() == 1:
            h = h[None]
        return h, h.shape[0]

    def xxz(self, J: float = 1.0, Delta: float = 1.0, h=None, sector: Optional[int] = None) -> Operator:
        """
        XXZ Hamiltonian with fields h (scalar, (L,) or (batch, L): one Hamiltonian per row).
        sector=None: operator in the full algebra (all sectors, padded); sector=N: operator in
        sector_algebra(N) only.
        """
        h, batch = self._fields(h)
        if sector is not None:
            alg = self.sector_algebra(sector)
            return alg.operator(self._block(sector, J, Delta, h, batch)[:, None], is_self_adjoint=True)
        op = self.algebra.from_blocks([self._block(N, J, Delta, h, batch) for N in range(self.L + 1)])
        op._is_self_adjoint = True
        return op

    def heisenberg(self, J: float = 1.0, h=None, sector: Optional[int] = None) -> Operator:
        return self.xxz(J, 1.0, h, sector)

    def random_field_heisenberg(self, W: float, batch_size: int, J: float = 1.0, sector: Optional[int] = None,
                                generator: Optional[torch.Generator] = None) -> Operator:
        """Heisenberg chain with i.i.d. fields h_i ~ U[-W, W]; one disorder realisation per batch element."""
        h = (2 * torch.rand(batch_size, self.L, generator=generator, dtype=self.real_dtype) - 1) * W
        return self.xxz(J, 1.0, h.to(self.device), sector)

    def sz(self, i: int, sector: Optional[int] = None) -> Operator:
        """S^z on site i."""
        return self._diagonal(lambda b: b[:, i] - 0.5, sector)

    def number(self, sector: Optional[int] = None) -> Operator:
        """N = number of up spins (central: N * 1 in sector N)."""
        return self._diagonal(lambda b: b.sum(-1), sector)

    def _diagonal(self, fn, sector):
        if sector is not None:
            alg = self.sector_algebra(sector)
            d = fn(self.bits(sector).to(self.real_dtype)).to(self.dtype)
            return alg.operator(torch.diag_embed(d)[None, None], is_self_adjoint=True)
        op = self.algebra.from_blocks([torch.diag_embed(fn(self.bits(N).to(self.real_dtype)).to(self.dtype))
                                       for N in range(self.L + 1)])
        op._is_self_adjoint = True
        return op

    # ------------------------------------------------------------------
    def basis_state(self, bitstring: Union[str, Sequence[int]]) -> DensityMatrix:
        """|s><s| for a configuration given as '0110...' (site 0 first) or a 0/1 sequence."""
        bits = [int(b) for b in bitstring]
        assert len(bits) == self.L
        label = sum(b << i for i, b in enumerate(bits))
        N = sum(bits)
        mat = torch.zeros(1, self.algebra.C, self.algebra.k_max, self.algebra.k_max, dtype=self.dtype,
                          device=self.device)
        p = self.index_of[label]
        mat[0, N, p, p] = 1.0
        return DensityMatrix(self.algebra, matrix=mat)

    def vector_in_sector(self, bitstring) -> torch.Tensor:
        """Basis vector of a configuration as a column of its sector (shape (binom(L, N),))."""
        bits = [int(b) for b in bitstring]
        label = sum(b << i for i, b in enumerate(bits))
        v = torch.zeros(self.sector_dim(sum(bits)), dtype=self.dtype, device=self.device)
        v[self.index_of[label]] = 1.0
        return v

    # ------------------------------------------------------------------
    def reduced_state(self, rho: Operator, n_left: int, sector: Optional[int] = None) -> DensityMatrix:
        """
        Partial trace over sites n_left, ..., L-1. rho lives in the full algebra (sector=None) or in
        sector_algebra(sector). Returns a state of SpinChain(n_left, ...).algebra; for U(1)
        symmetric states the result is block diagonal in the particle number of the left part.
        """
        left = SpinChain(n_left, 'open', self.complex_valued, self.precision, self.device)
        nB = self.L - n_left
        mask = (1 << n_left) - 1
        mat = rho.matrix
        B = mat.shape[0]
        out = torch.zeros(B, left.algebra.C, left.algebra.k_max, left.algebra.k_max, dtype=mat.dtype,
                          device=mat.device)
        sectors = range(self.L + 1) if sector is None else [sector]
        right_index = SpinChain(nB, 'open', self.complex_valued, self.precision, self.device).index_of if nB else None
        for N in sectors:
            block = mat[:, N if sector is None else 0]
            lab = self.labels[N]
            a, b = lab & mask, lab >> n_left
            ca = torch.tensor([bin(int(x)).count('1') for x in a.tolist()], device=mat.device)
            for c in range(max(0, N - nB), min(N, n_left) + 1):
                sel = torch.nonzero(ca == c).squeeze(-1)
                ia = left.index_of[a[sel]]
                ib = right_index[b[sel]] if nB else torch.zeros_like(ia)
                dA, dB = left.sector_dim(c), math.comb(nB, N - c)
                sub = block[:, sel][:, :, sel]                                   # (B, s, s)
                R = torch.zeros(B, dA, dB, dA, dB, dtype=mat.dtype, device=mat.device)
                R[:, ia[:, None], ib[:, None], ia[None, :], ib[None, :]] = sub
                out[:, c, :dA, :dA] += torch.diagonal(R, dim1=2, dim2=4).sum(-1)
        if isinstance(rho, DensityMatrix):
            return DensityMatrix(left.algebra, matrix=out, validate=False)
        return Operator(left.algebra, matrix=out)

    def entanglement_entropy(self, rho: Operator, n_left: int, sector: Optional[int] = None) -> torch.Tensor:
        """von Neumann entropy of the reduced state of sites 0..n_left-1."""
        from .states import von_neumann_entropy
        return von_neumann_entropy(self.reduced_state(rho, n_left, sector))

    # ------------------------------------------------------------------
    def site_amplitude_damping(self, i: int, gamma: float) -> InterSectorChannel:
        """
        Spin flip down on site i with probability gamma (a map between sectors N -> N-1):
        K_0 = |down><down|_i + sqrt(1 - gamma) |up><up|_i,  K_1 = sqrt(gamma) S-_i.
        """
        alg = self.algebra
        K = torch.zeros(1, 2, alg.C, alg.C, alg.k_max, alg.k_max, dtype=self.dtype, device=self.device)
        for N in range(self.L + 1):
            lab = self.labels[N]
            up = ((lab >> i) & 1).bool()
            idx = torch.arange(len(lab), device=self.device)
            keep = torch.full((len(lab),), 1.0, dtype=self.real_dtype, device=self.device)
            keep[up] = math.sqrt(1 - gamma)
            K[0, 0, N, N, idx, idx] = keep.to(self.dtype)
            if N > 0:
                src = idx[up]
                dst = self.index_of[lab[up] ^ (1 << i)]
                K[0, 1, N - 1, N, dst, src] = gamma ** 0.5
        return InterSectorChannel(alg, alg, K)

    # ------------------------------------------------------------------
    @staticmethod
    def level_spacing_ratio(energies: torch.Tensor) -> torch.Tensor:
        """
        Mean ratio of consecutive level spacings r = <min(s_n, s_{n+1}) / max(s_n, s_{n+1})> per row of
        sorted energies (Oganesyan-Huse); Poisson 2 ln 2 - 1 = 0.386, GOE 0.5307 (Atas et al. 2013).
        """
        e = torch.sort(energies, dim=-1)[0]
        s = torch.diff(e, dim=-1)
        a, b = s[..., :-1], s[..., 1:]
        r = torch.minimum(a, b) / torch.clamp(torch.maximum(a, b), min=1e-300)
        return r.mean(dim=-1)
