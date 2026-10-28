"""
Quantum channels on a Type I algebra M = (+)_c M_{k_c} that preserve the sectors.

A channel is stored by its Kraus operators, one tensor of shape
(batch, r, C, k_max, k_max): r Kraus operators K_i = (+)_c K_{i,c}, block c in the
top-left k_c x k_c corner as for operators. In the Schroedinger picture

    Phi(rho) = sum_i K_i rho K_i^*,

trace preserving iff sum_i K_i^* K_i = 1 and unital iff sum_i K_i K_i^* = 1.
The Heisenberg-picture (dual) channel has Kraus operators K_i^*. Because every block is
mapped to itself, the duality Tr(Phi(rho) A) = Tr(rho Phi^*(A)) holds blockwise and hence for
each of the traces Tr_blunt, Tr_norm and tau_vN.

The Kraus batch dimension is 1 (one channel for all samples) or equal to the batch of the
operators it acts on.
"""
import warnings
from typing import List, Optional, Sequence

import torch

from .algebra import TypeIAlgebra, _same_device, _as_tensor
from . import cost
from .states import DensityMatrix, _trace_weights

Operator = TypeIAlgebra.Operator


def _same_structure(a: TypeIAlgebra, b: TypeIAlgebra) -> bool:
    return a is b or (a.k_factors == b.k_factors and a.n_factors == b.n_factors and a.charges == b.charges
                      and a.dtype == b.dtype and _same_device(a.hilbert.device, b.hilbert.device))


def _check_input(op, alg: TypeIAlgebra, what: str) -> None:
    if not isinstance(op, Operator):
        raise TypeError(f"{what}: expected an Operator or DensityMatrix of {alg}, got {type(op).__name__} "
                        f"(wrap tensors with alg.operator(...))")
    if not _same_structure(op.algebra, alg):
        raise ValueError(f"{what}: the operator belongs to {op.algebra} ({op.algebra.dtype}), the channel to "
                         f"{alg} ({alg.dtype}); convert with op.cast(alg)")


def _kraus_sandwich(K: torch.Tensor, rho: torch.Tensor) -> torch.Tensor:
    """sum_i K_i rho K_i^* for Kraus operators K (B, r, n, k, k) and blocks rho (B, n, k, k)."""
    dtype = torch.promote_types(rho.dtype, K.dtype)
    B = max(rho.shape[0], K.shape[0])
    r, C, k = K.shape[1], K.shape[2], K.shape[-1]
    cost.check_memory(2 * cost.tensor_bytes((B, *K.shape[1:]), dtype), K.device,
                      f"Channel.apply (batch {B}, Kraus rank {r})")
    Kd = K.to(dtype)
    X = Kd @ rho.to(dtype).unsqueeze(1)                                    # (B, r, n, k, k)
    # sum_i X_i K_i^* as one product over the merged index (i, column)
    Xc = X.permute(0, 2, 3, 1, 4).reshape(B, C, k, r * k)
    Kc = Kd.permute(0, 2, 3, 1, 4).reshape(Kd.shape[0], C, k, r * k)
    return Xc @ Kc.conj().transpose(-2, -1)


class Channel:
    """Sector-preserving completely positive map given by Kraus operators."""

    def __init__(self, algebra: TypeIAlgebra, kraus: torch.Tensor, trace_preserving: Optional[bool] = None):
        """trace_preserving: known by construction (not checked); None: checked on first use with a state."""
        if kraus.dim() == 4:                   # (r, C, k, k): a single channel
            kraus = kraus.unsqueeze(0)
        expected = (algebra.C, algebra.k_max, algebra.k_max)
        if kraus.dim() != 5 or tuple(kraus.shape[2:]) != expected:
            raise ValueError(f"Kraus tensor must have shape (batch, r, {algebra.C}, {algebra.k_max}, "
                             f"{algebra.k_max}), got {tuple(kraus.shape)}")
        kraus = algebra._to_field(kraus, "Kraus operators")
        if any(k_c < algebra.k_max for k_c in algebra.k_factors):        # the padding must be zero
            kraus = kraus * algebra._const('block_mask', torch.bool)
        self.algebra = algebra
        self.kraus = kraus
        self._trace_preserving = trace_preserving

    # ------------------------------------------------------------------
    @property
    def batch_size(self) -> int:
        return self.kraus.shape[0]

    @property
    def kraus_rank(self) -> int:
        return self.kraus.shape[1]

    def __repr__(self) -> str:
        return f"Channel(batch={self.batch_size}, kraus_rank={self.kraus_rank}, C={self.algebra.C}, k_max={self.algebra.k_max})"

    # ------------------------------------------------------------------
    def apply(self, op: Operator) -> Operator:
        """
        Schroedinger picture: rho -> sum_i K_i rho K_i^* (lazy).
        A DensityMatrix is mapped to a DensityMatrix if the channel is trace preserving;
        otherwise the result is a plain Operator and a warning is issued.
        """
        _check_input(op, self.algebra, "Channel.apply")
        K = self.kraus
        alg = self.algebra

        def generator():
            if alg._packed:                     # the Kraus blocks of each size applied to the stored blocks
                return [_kraus_sandwich(K[:, :, alg._sel(idx), :k, :k], rho)
                        for (k, idx), rho in zip(alg.size_groups, op.parts)]
            return _kraus_sandwich(K, op.matrix)

        if isinstance(op, DensityMatrix):
            if self._trace_preserving is None:                 # checked once per channel (one synchronisation)
                self._trace_preserving = self.is_trace_preserving()
            if self._trace_preserving:
                return DensityMatrix(self.algebra, generator=generator, validate=False)
            warnings.warn("channel is not trace preserving: the image of a DensityMatrix is "
                          "returned as an unnormalised Operator", stacklevel=2)
        out = Operator(self.algebra, generator=generator)
        if op._is_self_adjoint:
            out._is_self_adjoint = True
        if op._is_positive:
            out._is_positive = True
        return out

    __call__ = apply

    def adjoint(self) -> 'Channel':
        """Heisenberg-picture channel A -> sum_i K_i^* A K_i."""
        return Channel(self.algebra, self.kraus.conj().transpose(-2, -1))

    dual = adjoint

    def compose(self, other: 'Channel') -> 'Channel':
        """self o other: first `other`, then `self` (Kraus operators K_i L_j)."""
        if isinstance(other, InterSectorChannel):
            return self.to_inter_sector().compose(other)
        if not isinstance(other, Channel) or not _same_structure(other.algebra, self.algebra):
            raise ValueError("compose: both channels must act on the same algebra")
        K = self.kraus.unsqueeze(2)            # (B, r1, 1, C, k, k)
        L = other.kraus.unsqueeze(1)           # (B, 1, r2, C, k, k)
        KL = K @ L
        B = KL.shape[0]
        tp = True if (self._trace_preserving and other._trace_preserving) else None
        return Channel(self.algebra, KL.reshape(B, -1, *KL.shape[3:]), trace_preserving=tp)

    __matmul__ = compose

    def to_inter_sector(self) -> 'InterSectorChannel':
        """The same map as an InterSectorChannel with algebra_in = algebra_out."""
        K = self.kraus
        B, r, C, k, _ = K.shape
        G = torch.zeros(B, r, C, C, k, k, dtype=K.dtype, device=K.device)
        idx = torch.arange(C, device=K.device)
        G[:, :, idx, idx] = K
        return InterSectorChannel(self.algebra, self.algebra, G, trace_preserving=self._trace_preserving)

    def mix(self, other: 'Channel', p: float) -> 'Channel':
        """Convex combination (1 - p) self + p other."""
        if not _same_structure(other.algebra, self.algebra):
            raise ValueError("mix: both channels must act on the same algebra")
        if isinstance(p, torch.Tensor) and p.numel() != 1:
            raise ValueError("Channel.mix: p must be a single number (mix channels of a batch one by one)")
        p = float(p)
        if not 0 <= p <= 1:
            raise ValueError("mix: p must lie in [0, 1]")
        B = max(self.batch_size, other.batch_size)
        a = (1 - p) ** 0.5 * self.kraus.expand(B, *self.kraus.shape[1:])
        b = p ** 0.5 * other.kraus.expand(B, *other.kraus.shape[1:]).to(a.dtype)
        tp = True if (self._trace_preserving and other._trace_preserving) else None
        return Channel(self.algebra, torch.cat([a, b], dim=1), trace_preserving=tp)

    # ------------------------------------------------------------------
    def _active_identity(self) -> torch.Tensor:
        return self.algebra._const('eye', self.kraus.dtype)

    def _deviation(self, gram: torch.Tensor) -> torch.Tensor:
        return (gram - self._active_identity()).abs().amax(dim=(-3, -2, -1))

    def trace_preservation_error(self) -> torch.Tensor:
        """max |sum_i K_i^* K_i - 1| per channel in the batch."""
        K = self.kraus
        return self._deviation((K.conj().transpose(-2, -1) @ K).sum(dim=1))

    def unitality_error(self) -> torch.Tensor:
        """max |sum_i K_i K_i^* - 1| per channel in the batch."""
        K = self.kraus
        return self._deviation((K @ K.conj().transpose(-2, -1)).sum(dim=1))

    def _tol(self) -> float:
        return 1e-5 if self.kraus.dtype in (torch.float32, torch.complex64) else 1e-10

    def is_trace_preserving(self, tol: Optional[float] = None) -> bool:
        return bool(torch.all(self.trace_preservation_error() <= (self._tol() if tol is None else tol)))

    def is_unital(self, tol: Optional[float] = None) -> bool:
        return bool(torch.all(self.unitality_error() <= (self._tol() if tol is None else tol)))

    # ------------------------------------------------------------------
    def choi(self) -> List[torch.Tensor]:
        """
        Choi matrices per channel c, J_c = sum_ij |i><j| (x) Phi(|i><j|), shape (batch, k_c^2, k_c^2).
        Positive semidefinite for a completely positive map; Tr_out J = 1 iff trace preserving.
        """
        cost.check_memory(sum(cost.tensor_bytes((self.batch_size, k ** 2, k ** 2), self.kraus.dtype)
                              for k in self.algebra.k_factors), self.kraus.device, "Channel.choi")

        def f(K, k):                                                  # K (B, r, n, k, k), index (a, i)
            v = K.transpose(-2, -1).reshape(*K.shape[:3], k * k)      # index (i, a)
            return torch.einsum('brnx,brny->bnxy', v, v.conj())
        return self._per_size(f)

    def _per_size(self, f) -> List[torch.Tensor]:
        """f(K_group, k) for the Kraus blocks of each block-size group (one call per size), split per sector."""
        alg = self.algebra
        out = [torch.zeros(self.batch_size, 0, 0, dtype=self.kraus.dtype, device=self.kraus.device)] * alg.C
        for k, idx in alg._size_groups():
            res = f(self.kraus[:, :, alg._sel(idx), :k, :k], k)       # (B, n, k^2, k^2)
            for j, c in enumerate(idx):
                out[c] = res[:, j]
        return out

    def superoperator(self) -> List[torch.Tensor]:
        """Matrix of Phi on row-major vec(rho), per channel: S_c = sum_i K_i (x) conj(K_i), (batch, k_c^2, k_c^2)."""
        cost.check_memory(sum(cost.tensor_bytes((self.batch_size, k ** 2, k ** 2), self.kraus.dtype)
                              for k in self.algebra.k_factors), self.kraus.device, "Channel.superoperator")

        def f(K, k):
            return torch.einsum('zrnai,zrnbj->znabij', K, K.conj()).reshape(K.shape[0], K.shape[2], k * k, k * k)
        return self._per_size(f)

    @staticmethod
    def from_superoperator(algebra: TypeIAlgebra, blocks: Sequence[torch.Tensor], tol: Optional[float] = None) -> 'Channel':
        """
        Channel from per-channel superoperators S_c of shape (batch, k_c^2, k_c^2) acting on row-major
        vec(rho). Kraus operators are obtained from the eigendecomposition of the Choi matrix; S_c must
        describe a completely positive map. A Choi matrix J with ||J - J^*|| > 2 tol ||J|| (the map does not
        preserve Hermiticity) or with an eigenvalue below -tol ||J|| raises, where ||.|| is the spectral norm
        (taken at least 1); by default tol = 1e-10 for double and 1e-4 for single-precision input. Small
        negative eigenvalues within the tolerance are clamped to zero.
        """
        if len(blocks) != algebra.C:
            raise ValueError(f"from_superoperator: need {algebra.C} blocks, got {len(blocks)}")
        k_max, C = algebra.k_max, algebra.C
        blocks = [_as_tensor(S) for S in blocks]
        batch = max(S.shape[0] if S.dim() == 3 else 1 for S, k_c in zip(blocks, algebra.k_factors) if k_c)
        r = k_max * k_max
        kraus = torch.zeros(batch, r, C, k_max, k_max, dtype=algebra.hilbert.dtype, device=algebra.hilbert.device)
        checks = []                                    # validated together at the end: one synchronisation
        for k_c, idx in algebra._size_groups():        # one Choi eigendecomposition per block size
            S = torch.stack([(blocks[c] if blocks[c].dim() == 3 else blocks[c].unsqueeze(0)).expand(batch, -1, -1)
                             .to(device=algebra.hilbert.device) for c in idx], dim=1)      # (B, n, k^2, k^2)
            n = len(idx)
            # J[(i,a),(j,b)] = S[(a,b),(i,j)]
            J = S.reshape(batch, n, k_c, k_c, k_c, k_c).permute(0, 1, 4, 2, 5, 3).reshape(batch * n, k_c * k_c,
                                                                                       k_c * k_c)
            default_tol = tol if tol is not None else (1e-10 if J.real.dtype == torch.float64 else 1e-4)
            # per-matrix scale (norm of J, at least 1 so that the zero map passes)
            scale = torch.clamp(torch.linalg.matrix_norm(J, ord=2), min=1.0)
            skew = torch.linalg.matrix_norm(J - J.conj().transpose(-2, -1), ord=2)
            checks.append(("superoperator does not preserve Hermiticity (Choi matrix not Hermitian)",
                           (skew > 2 * default_tol * scale).any()))
            J = 0.5 * (J + J.conj().transpose(-2, -1))
            w, V = cost.batched_call(torch.linalg.eigh, J, f"Choi eigendecomposition (blocks of size {k_c})")
            checks.append(("superoperator is not completely positive",
                           (w.amin(dim=-1) < -default_tol * scale).any()))
            w = torch.clamp(w, min=0.0)
            # K_m[a, i] = sqrt(w_m) V[(i, a), m]
            Km = (V * w.sqrt().unsqueeze(-2).to(V.dtype)).transpose(-2, -1)
            Km = Km.reshape(batch, n, k_c * k_c, k_c, k_c).transpose(-2, -1)
            if torch.is_complex(Km) and not torch.is_complex(kraus):
                checks.append(("complex Kraus operators require a complex-valued algebra",
                               Km.imag.abs().amax() > 1e-5 * torch.clamp(Km.abs().amax(), min=1.0)))
                Km = Km.real
            kraus[:, :k_c * k_c, algebra._sel(idx), :k_c, :k_c] = Km.transpose(1, 2).to(kraus.dtype)
        if checks:
            for (msg, _), bad in zip(checks, torch.stack([c for _, c in checks]).tolist()):
                if bad:
                    raise ValueError(msg)
        return Channel(algebra, kraus)


# ======================================================================
# Constructors
# ======================================================================
def _empty_kraus(alg: TypeIAlgebra, r: int, batch: int = 1, dtype=None) -> torch.Tensor:
    return torch.zeros(batch, r, alg.C, alg.k_max, alg.k_max,
                       dtype=dtype or alg.hilbert.dtype, device=alg.hilbert.device)


def identity_channel(alg: TypeIAlgebra) -> Channel:
    return Channel(alg, alg._const('eye').reshape(1, 1, alg.C, alg.k_max, alg.k_max).clone(), trace_preserving=True)


def unitary_channel(U: Operator) -> Channel:
    """rho -> U rho U^* (U unitary in the algebra, possibly batched)."""
    return Channel(U.algebra, U.matrix.unsqueeze(1))


def kraus_channel(ops: Sequence[Operator]) -> Channel:
    """Channel with the given Kraus operators (all from the same algebra)."""
    alg = ops[0].algebra
    mats = [op.matrix for op in ops]
    B = max(m.shape[0] for m in mats)
    return Channel(alg, torch.stack([m.expand(B, *m.shape[1:]) for m in mats], dim=1))


def dephasing_channel(alg: TypeIAlgebra, p: float) -> Channel:
    """rho_c -> (1 - p) rho_c + p diag(rho_c) in the standard basis of each block (unital)."""
    if not 0 <= p <= 1:
        raise ValueError("dephasing_channel: p must lie in [0, 1]")
    km = alg.k_max
    rows = alg._const('rows').to(alg.dtype)                            # (C, km): active diagonal entries
    K = _empty_kraus(alg, 1 + km)
    K[0, 0] = (1 - p) ** 0.5 * alg._const('eye')
    i = torch.arange(km, device=K.device)
    K[0, 1 + i, :, i, i] = p ** 0.5 * rows.T                          # K_{1+i} = sqrt(p) |i><i| in every block
    return Channel(alg, K, trace_preserving=True)


def depolarizing_channel(alg: TypeIAlgebra, p: float) -> Channel:
    """
    rho_c -> (1 - p) rho_c + p Tr(rho_c) 1/k_c in every block (unital, keeps the weight of each
    sector). p = 1 is the trace-preserving conditional expectation onto the centre of M.
    """
    if not 0 <= p <= 1:
        raise ValueError("depolarizing_channel: p must lie in [0, 1]")
    km = alg.k_max
    K = _empty_kraus(alg, 1 + km * km)
    K[0, 0] = (1 - p) ** 0.5 * alg._const('eye')
    i = torch.arange(km, device=K.device).repeat_interleave(km)
    j = torch.arange(km, device=K.device).repeat(km)
    k = alg._const('k', alg.hilbert.real_dtype)
    amp = torch.where(k > 0, (p / k.clamp(min=1)) ** 0.5, torch.zeros_like(k))            # (C,)
    vals = alg._const('block_mask', torch.bool)[:, i, j].T * amp                          # (km^2, C), zero on padding
    K[0, 1 + i * km + j, :, i, j] = vals.to(K.dtype)                  # K_{ij} = sqrt(p/k_c) |i><j| in every block
    return Channel(alg, K, trace_preserving=True)


def center_expectation(alg: TypeIAlgebra) -> Channel:
    """Conditional expectation onto the centre Z(M): A -> sum_c (Tr A_c / k_c) P_c."""
    return depolarizing_channel(alg, 1.0)


def amplitude_damping_channel(alg: TypeIAlgebra, gamma: float) -> Channel:
    """
    Decay of every level i > 0 to the level 0 of its block with probability gamma:
    K_0 = |0><0| + sqrt(1 - gamma) sum_{i>0} |i><i|,  K_i = sqrt(gamma) |0><i|.
    Trace preserving, not unital for gamma > 0.
    """
    if not 0 <= gamma <= 1:
        raise ValueError("amplitude_damping_channel: gamma must lie in [0, 1]")
    km = alg.k_max
    rows = alg._const('rows').to(alg.dtype)                            # (C, km)
    K = _empty_kraus(alg, km)
    d = torch.full((km,), (1 - gamma) ** 0.5, dtype=alg.dtype, device=K.device)
    d[0] = 1.0
    K[0, 0] = torch.diag_embed(rows * d)
    i = torch.arange(1, km, device=K.device)
    K[0, i, :, 0, i] = gamma ** 0.5 * rows[:, 1:].T                     # K_i = sqrt(gamma) |0><i|
    return Channel(alg, K, trace_preserving=True)


def measurement_channel(projections: Sequence[Operator]) -> Channel:
    """Non-selective Lueders measurement rho -> sum_j P_j rho P_j for a resolution of the identity."""
    return kraus_channel(projections)


def random_channel(alg: TypeIAlgebra, kraus_rank: int, batch_size: int = 1) -> Channel:
    """
    Random trace-preserving channel from a Haar-random Stinespring isometry
    V_c: C^{k_c} -> C^r (x) C^{k_c} in every block (environment first), K_{i,c} = (<i| (x) 1) V_c.
    """
    K = _empty_kraus(alg, kraus_rank, batch_size)
    for k, idx in alg._size_groups():                  # one draw for all sectors of the same size
        V = alg.random_isometry(k * kraus_rank, k, batch_size=batch_size * len(idx))   # (B n, r k, k)
        V = V.reshape(batch_size, len(idx), kraus_rank, k, k).transpose(1, 2)
        K[:, :, alg._sel(idx), :k, :k] = V
    return Channel(alg, K, trace_preserving=True)


def random_mixed_unitary_channel(alg: TypeIAlgebra, n_unitaries: int, batch_size: int = 1,
                                 weights: Optional[torch.Tensor] = None) -> Channel:
    """rho -> sum_j p_j U_j rho U_j^* with Haar unitaries in M and Dirichlet(1) (or given) weights; unital."""
    rdt = alg.hilbert.real_dtype
    if weights is None:                       # in the algebra's precision, so that the channel is TP to round-off
        weights = torch.distributions.Dirichlet(torch.ones(n_unitaries, dtype=rdt)).sample((batch_size,))
    else:
        weights = _as_tensor(weights)
        if weights.shape[-1] != n_unitaries or weights.dim() > 2 or (weights.dim() == 2
                                                                    and weights.shape[0] not in (1, batch_size)):
            raise ValueError(f"weights must have shape ({n_unitaries},) or (batch_size, {n_unitaries})")
        if bool((weights < 0).any()) or bool((weights.sum(dim=-1) <= 0).any()):
            raise ValueError("weights must be non-negative with a positive sum")
    weights = weights.to(device=alg.hilbert.device, dtype=rdt)
    weights = weights.expand(batch_size, n_unitaries) if weights.dim() == 1 or weights.shape[0] == 1 else weights
    weights = weights / weights.sum(dim=-1, keepdim=True)
    # all unitaries in one draw
    Us = alg.random_unitary_operator(batch_size * n_unitaries).matrix.reshape(batch_size, n_unitaries, alg.C,
                                                                              alg.k_max, alg.k_max)
    return Channel(alg, Us * weights.sqrt()[:, :, None, None, None].to(Us.dtype), trace_preserving=True)


# ======================================================================
# Maps between sectors / between different algebras
# ======================================================================
class InterSectorChannel:
    """
    Completely positive map Phi: M_in = (+)_c M_{k_c} -> M_out = (+)_d M_{m_d} between (possibly
    different) Type I algebras, in Schroedinger form

        Phi(rho)_d = sum_c sum_i K_i^{dc} rho_c (K_i^{dc})^*,      K_i^{dc}: C^{k_c} -> C^{m_d}.

    Every normal CP map between such algebras has this form: Phi = sum_{d,c} Phi_{dc} with
    Phi_{dc}: M_{k_c} -> M_{m_d} CP. Kraus operators are stored as one tensor of shape
    (batch, r, D, C, m_max, k_max); blocks of rank < r are zero padded. The classical part of
    Phi is the matrix of sector transitions T[d, c] = Tr Phi_{dc}(rho_c) / Tr rho_c.

    Trace preserving iff sum_{d,i} K_i^{dc*} K_i^{dc} = 1_c for every input sector c.
    """

    def __init__(self, algebra_in: TypeIAlgebra, algebra_out: TypeIAlgebra, kraus: torch.Tensor,
                 trace_preserving: Optional[bool] = None):
        if kraus.dim() == 5:
            kraus = kraus.unsqueeze(0)
        expected = (algebra_out.C, algebra_in.C, algebra_out.k_max, algebra_in.k_max)
        if kraus.dim() != 6 or tuple(kraus.shape[2:]) != expected:
            raise ValueError(f"Kraus tensor must have shape (batch, r, {expected[0]}, {expected[1]}, "
                             f"{expected[2]}, {expected[3]}), got {tuple(kraus.shape)}")
        kraus = algebra_out._to_field(kraus, "Kraus operators")
        if any(k < algebra_out.k_max for k in algebra_out.k_factors) or any(k < algebra_in.k_max
                                                                             for k in algebra_in.k_factors):
            rows_out = algebra_out._const('rows')            # (D, m_max)
            rows_in = algebra_in._const('rows').to(rows_out.device)     # (C, k_max)
            kraus = kraus * (rows_out[:, None, :, None] & rows_in[None, :, None, :])   # zero padding
        self.algebra_in = algebra_in
        self.algebra_out = algebra_out
        self.kraus = kraus
        self._trace_preserving = trace_preserving

    def __repr__(self) -> str:
        return (f"InterSectorChannel(batch={self.kraus.shape[0]}, kraus_rank={self.kraus.shape[1]}, "
                f"C_in={self.algebra_in.C}, C_out={self.algebra_out.C})")

    # ------------------------------------------------------------------
    @staticmethod
    def from_blocks(algebra_in: TypeIAlgebra, algebra_out: TypeIAlgebra, blocks: dict,
                    batch_size: int = 1) -> 'InterSectorChannel':
        """
        Build from {(d, c): [K_1, K_2, ...]} with K_i of shape (m_d, k_c) or (batch, m_d, k_c);
        missing pairs (d, c) are zero maps.
        """
        r = max(len(v) for v in blocks.values())
        dtype = algebra_out.hilbert.dtype
        K = torch.zeros(batch_size, r, algebra_out.C, algebra_in.C, algebra_out.k_max, algebra_in.k_max,
                        dtype=dtype, device=algebra_out.hilbert.device)
        for (d, c), ops in blocks.items():
            m_d, k_c = algebra_out.k_factors[d], algebra_in.k_factors[c]
            for i, op in enumerate(ops):
                op = _as_tensor(op)
                assert op.shape[-2:] == (m_d, k_c), f"block ({d}, {c}) must be {m_d} x {k_c}"
                K[:, i, d, c, :m_d, :k_c] = algebra_out._to_field(op, f"Kraus block ({d}, {c})").to(dtype)
        return InterSectorChannel(algebra_in, algebra_out, K)

    # ------------------------------------------------------------------
    def apply(self, op: Operator) -> Operator:
        """Schroedinger picture (lazy). Trace-preserving maps send a DensityMatrix to a DensityMatrix."""
        _check_input(op, self.algebra_in, "InterSectorChannel.apply")
        K = self.kraus

        def generator():
            rho = op.matrix
            dtype = torch.promote_types(rho.dtype, K.dtype)
            B, r, D, C, m, _ = K.shape
            B = max(B, rho.shape[0])
            cost.check_memory(2 * cost.tensor_bytes((B, r, D, C, m, max(m, rho.shape[-1])), dtype), K.device,
                              f"InterSectorChannel.apply (batch {B}, Kraus rank {r}, {D}x{C} sector pairs)")
            Kd = K.to(dtype)
            X = Kd @ rho.to(dtype)[:, None, None]                          # (B, r, D, C, m, k)
            k = X.shape[-1]
            # sum over (i, c) of X K^* as one product over the merged index (i, c, column)
            Xc = X.permute(0, 2, 4, 1, 3, 5).reshape(B, D, m, r * C * k)
            Kc = Kd.permute(0, 2, 4, 1, 3, 5).reshape(Kd.shape[0], D, m, r * C * k)
            return Xc @ Kc.conj().transpose(-2, -1)                        # (B, D, m, m)

        if isinstance(op, DensityMatrix):
            if self._trace_preserving is None:
                self._trace_preserving = self.is_trace_preserving()
            if self._trace_preserving:
                return DensityMatrix(self.algebra_out, generator=generator, validate=False)
            warnings.warn("map is not trace preserving: the image of a DensityMatrix is returned "
                          "as an unnormalised Operator", stacklevel=2)
        out = Operator(self.algebra_out, generator=generator)
        if op._is_self_adjoint:
            out._is_self_adjoint = True
        if op._is_positive:
            out._is_positive = True
        return out

    __call__ = apply

    def adjoint(self, trace: str = 'blunt') -> 'InterSectorChannel':
        """
        Dual map Phi^*: M_out -> M_in with tr_out(Phi(rho) A) = tr_in(rho Phi^*(A)) for the chosen
        trace ('blunt', 'norm' or 'tau_vN' on both algebras). With weights w_c of the trace,
        Phi^*(A)_c = sum_{d,i} (w^out_d / w^in_c) K_i^{dc*} A_d K_i^{dc}.
        """
        Kd = self.kraus.conj().transpose(-2, -1).transpose(2, 3)          # (B, r, C, D, k, m)
        if trace != 'blunt':
            w_in = _trace_weights(self.algebra_in, trace)
            w_out = _trace_weights(self.algebra_out, trace)
            scale = (w_out[None, :] / w_in[:, None]).sqrt()               # (C, D)
            Kd = Kd * scale[None, None, :, :, None, None].to(device=Kd.device, dtype=Kd.dtype)
        return InterSectorChannel(self.algebra_out, self.algebra_in, Kd)

    dual = adjoint

    def compose(self, other: 'InterSectorChannel') -> 'InterSectorChannel':
        """self o other. Kraus operators K_j^{ed} L_i^{dc}, indexed by (j, i, d)."""
        if isinstance(other, Channel):
            other = other.to_inter_sector()
        if not _same_structure(other.algebra_out, self.algebra_in):
            raise ValueError(f"compose: the output algebra of the first map ({other.algebra_out}) is not the input "
                             f"algebra of the second ({self.algebra_in})")
        K = self.kraus.unsqueeze(2).unsqueeze(5)          # (B, r2, 1, E, D, 1, m_e, m_d)
        L = other.kraus.unsqueeze(1).unsqueeze(3)         # (B, 1, r1, 1, D, C, m_d, k_c)
        KL = K @ L                                        # (B, r2, r1, E, D, C, m_e, k_c)
        B, r2, r1, E, D, C = KL.shape[:6]
        KL = KL.permute(0, 1, 2, 4, 3, 5, 6, 7).reshape(B, r2 * r1 * D, E, C, *KL.shape[-2:])
        tp = True if (self._trace_preserving and other._trace_preserving) else None
        return InterSectorChannel(other.algebra_in, self.algebra_out, KL, trace_preserving=tp)

    __matmul__ = compose

    # ------------------------------------------------------------------
    def _identity_in(self) -> torch.Tensor:
        return self.algebra_in._const('eye', self.kraus.dtype)

    def trace_preservation_error(self) -> torch.Tensor:
        K = self.kraus
        gram = (K.conj().transpose(-2, -1) @ K).sum(dim=(1, 2))           # (B, C, k, k)
        return (gram - self._identity_in()).abs().amax(dim=(-3, -2, -1))

    def is_trace_preserving(self, tol: Optional[float] = None) -> bool:
        tol = (1e-5 if self.kraus.dtype in (torch.float32, torch.complex64) else 1e-10) if tol is None else tol
        return bool(torch.all(self.trace_preservation_error() <= tol))

    def is_unital(self, tol: Optional[float] = None) -> bool:
        """Phi(1) = 1 (Schroedinger picture, blunt trace)."""
        tol = (1e-5 if self.kraus.dtype in (torch.float32, torch.complex64) else 1e-10) if tol is None else tol
        alg = self.algebra_in
        one = alg.operator(self._identity_in().unsqueeze(0).expand(self.kraus.shape[0], -1, -1, -1).clone())
        out = self.apply(one).matrix
        return bool((out - self.algebra_out._const('eye', out.dtype)).abs().max() <= tol)

    def transition_matrix(self, rho: Optional[Operator] = None) -> torch.Tensor:
        """
        Sector transition probabilities T[b, d, c] = Tr Phi_{dc}(rho_c) / Tr rho_c, for the given state
        or, by default, for the maximally mixed state of every input sector. Columns sum to 1 for a
        trace-preserving map.
        """
        alg = self.algebra_in
        if rho is None:
            mat = (self._identity_in() * alg._const('w_norm', self.kraus.real.dtype)[:, None, None]).unsqueeze(0)
        else:
            mat = rho.matrix
            tr = torch.diagonal(mat, dim1=-2, dim2=-1).sum(-1).real
            mat = mat / torch.where(tr > 0, tr, torch.ones_like(tr))[:, :, None, None].to(mat.dtype)
        K = self.kraus
        dtype = torch.promote_types(mat.dtype, K.dtype)
        Kd = K.to(dtype)
        # Tr(K rho K^*) = sum_{m,i,j} K_mi rho_ij conj(K_mj), without forming the images
        return torch.einsum('brdcmi,bcij,brdcmj->bdc', Kd, mat.to(dtype), Kd.conj()).real


def random_inter_sector_channel(algebra_in: TypeIAlgebra, algebra_out: TypeIAlgebra, kraus_rank: int = 1,
                                batch_size: int = 1) -> InterSectorChannel:
    """
    Random trace-preserving map from a Haar-random Stinespring isometry per input sector,
    V_c: C^{k_c} -> C^r (x) ((+)_d C^{m_d}) (environment first), split into the blocks
    K_i^{dc} = (<i| (x) P_d) V_c.
    """
    if algebra_in.hilbert.complex_valued != algebra_out.hilbert.complex_valued:
        raise ValueError("both algebras must be real or both complex")
    m = algebra_out.k_factors
    n_out = sum(m)
    K = torch.zeros(batch_size, kraus_rank, algebra_out.C, algebra_in.C, algebra_out.k_max, algebra_in.k_max,
                    dtype=algebra_out.hilbert.dtype, device=algebra_out.hilbert.device)
    offsets = [sum(m[:d]) for d in range(len(m))]
    for k_c, idx in algebra_in._size_groups():          # one isometry draw per input block size
        if n_out * kraus_rank < k_c:
            raise ValueError("output too small for an isometry: increase kraus_rank")
        n = len(idx)
        V = algebra_out.random_isometry(n_out * kraus_rank, k_c, batch_size=batch_size * n)
        V = V.reshape(batch_size, n, kraus_rank, n_out, k_c).transpose(1, 2)     # (B, r, n, n_out, k_c)
        sel = algebra_in._sel(idx)
        for d, m_d in enumerate(m):
            if m_d:
                K[:, :, d, sel, :m_d, :k_c] = V[:, :, :, offsets[d]:offsets[d] + m_d, :]
    return InterSectorChannel(algebra_in, algebra_out, K, trace_preserving=True)


class SectorJump:
    """
    A sparse jump operator L of M = (+)_c M_{k_c} into itself that maps basis vectors to multiples of basis vectors
    ("monomial"): sector c goes into one sector d(c) (or to zero), and the basis vector a of sector c goes to
    val[c, a] times the basis vector dst[c, a] of d(c), injectively on the support. Spin flips S^+_i, S^-_i, local
    fields S^z_i, particle loss a_i and projections onto basis states are of this form. For a Lindblad term
    g (L rho L^* - 1/2 {L^*L, rho}) the library then never forms L: L rho L^* is a permutation of the entries of rho
    with weights val_a conj(val_b), L^*L is diagonal, and the dual L^* Z L is a gather. Cost O(sum_c k_c^2) per
    application instead of O(C^2 k_max^3) for the same jump as a dense InterSectorChannel.

    target: (C,) target sector d(c) of every sector (-1: L kills the sector); dst: (C, k_max) target positions;
    val: (C, k_max) or (batch, C, k_max) coefficients (zero on the padding and where L kills the basis vector). A
    batch of coefficients gives one jump per sample, e.g. sqrt(rate) folded in for a batch of noise strengths.
    """

    def __init__(self, algebra: TypeIAlgebra, target, dst, val):
        dev = algebra.hilbert.device
        C, k = algebra.C, algebra.k_max
        target = torch.as_tensor(target, dtype=torch.long, device=dev).reshape(C)
        dst = torch.as_tensor(dst, dtype=torch.long, device=dev).reshape(C, k)
        val = torch.as_tensor(val, device=dev)
        if val.dim() == 2:
            val = val.unsqueeze(0)
        if tuple(val.shape[1:]) != (C, k):
            raise ValueError(f"val must have shape (C, k_max) or (batch, C, k_max) = (.., {C}, {k}), got {tuple(val.shape)}")
        val = algebra._to_field(val, "SectorJump coefficients").to(algebra.dtype)
        rows = algebra._const('rows')                                        # (C, k_max)
        live = (val != 0).any(0) & rows & (target >= 0)[:, None]             # support of L, (C, k_max)
        val = val * live
        kk = torch.tensor(algebra.k_factors, device=dev)
        tgt = target.clamp(min=0)
        if bool(((dst >= kk[tgt][:, None]) & live).any()) or bool(((dst < 0) & live).any()):
            raise ValueError("SectorJump: a target position lies outside its target sector")
        flat = tgt[:, None] * k + dst                                        # (C, k_max) row of the target in (C*k_max)
        used = flat[live]
        if used.numel() != torch.unique(used).numel():
            raise ValueError("SectorJump: two basis vectors are mapped to the same basis vector (not injective)")
        self.algebra = algebra
        self.target, self.dst, self.val, self._live = target, dst, val, live
        self._flat_rows = flat
        # per sector c of the support: (c, d(c), support positions a, target positions dst_a); L rho L^* and L^* Z L
        # are gathered and scattered block by block (no index of size sum_c k_c^2 is stored)
        self._sectors = []
        for c in range(C):
            a = torch.nonzero(live[c]).squeeze(-1)
            if a.numel():
                self._sectors.append((c, int(target[c]), a, dst[c, a]))

    def __repr__(self) -> str:
        return f"SectorJump(batch={self.val.shape[0]}, C={self.algebra.C}, nnz={int(self._live.sum())})"

    def scaled(self, coef) -> 'SectorJump':
        """The jump c L with a scalar or a per-sample coefficient c of shape (batch,)."""
        c = torch.as_tensor(coef, device=self.val.device).to(self.val.dtype)
        c = c.reshape(-1, 1, 1) if c.dim() else c
        return SectorJump(self.algebra, self.target, self.dst, self.val * c)

    def _zeros(self, data, B: int):
        alg = self.algebra
        if isinstance(data, (list, tuple)):
            return [x.new_zeros(B, *x.shape[1:]) for x in data]
        return data.new_zeros(B, alg.C, alg.k_max, alg.k_max)

    def sandwich(self, rho):
        """L rho L^* for a block-diagonal rho: a tensor (B, C, k_max, k_max), or in the packed layout the list of
        group tensors (returned as such). O(sum_c (support of L in c)^2) per sample."""
        alg = self.algebra
        if alg._packed and not isinstance(rho, (list, tuple)):       # a padded tensor: through the blocks
            return alg._pad(self.sandwich(alg._unpad(rho)))
        first = rho[0] if isinstance(rho, (list, tuple)) else rho
        B = max(first.shape[0], self.val.shape[0])
        out = self._zeros(rho, B)
        for c, d, a, dd in self._sectors:
            X = alg._block_view(rho, c).index_select(-2, a).index_select(-1, a)     # (b, n, n)
            v = self.val[:, c, a].to(X.dtype)                                        # (bv, n)
            Y = v[:, :, None] * X * v.conj()[:, None, :]
            O = alg._block_view(out, d)
            O[:, dd[:, None], dd[None, :]] += Y                                      # dst is injective
        return out

    def dual_sandwich(self, Z):
        """L^* Z L (Heisenberg picture) for a block-diagonal Z (a padded tensor, or the list of group tensors)."""
        alg = self.algebra
        if alg._packed and not isinstance(Z, (list, tuple)):
            return alg._pad(self.dual_sandwich(alg._unpad(Z)))
        first = Z[0] if isinstance(Z, (list, tuple)) else Z
        B = max(first.shape[0], self.val.shape[0])
        out = self._zeros(Z, B)
        for c, d, a, dd in self._sectors:
            X = alg._block_view(Z, d).index_select(-2, dd).index_select(-1, dd)
            v = self.val[:, c, a].to(X.dtype)
            O = alg._block_view(out, c)
            O[:, a[:, None], a[None, :]] = v.conj()[:, :, None] * X * v[:, None, :]
        return out

    def LdL(self, padded: bool = False):
        """L^*L, diagonal, as a tensor (batch, C, k_max, k_max), or in the packed layout (unless padded=True) as the
        list of group tensors (batch, n_g, k_g, k_g)."""
        d = (self.val.abs() ** 2).to(self.val.dtype)
        if self.algebra._packed and not padded:
            return [torch.diag_embed(d[:, self.algebra._sel(idx), :k]) for k, idx in self.algebra.size_groups]
        return torch.diag_embed(d)

    def to_inter_sector(self) -> InterSectorChannel:
        """The same jump as a dense InterSectorChannel (Kraus rank 1)."""
        alg = self.algebra
        B, C, k = self.val.shape[0], alg.C, alg.k_max
        K = torch.zeros(B, 1, C, C, k, k, dtype=alg.dtype, device=self.val.device)
        c, a = torch.nonzero(self._live, as_tuple=True)
        K[:, 0, self.target[c], c, self.dst[c, a], a] = self.val[:, c, a]
        return InterSectorChannel(alg, alg, K)
