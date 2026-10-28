"""
Krylov-subspace methods for large sectors, where dense matrices no longer fit.

A Hamiltonian is only needed through its action on vectors, `matvec(V)` with V of shape
(batch, k) -> (batch, k) (or (batch, k, m) -> (batch, k, m)); every batch element may be a different operator
(e.g. one disorder realisation each). An operator of one sector is any object with matvec, batch_size, dim, dtype
and device (the protocol of `SectorOperator`); the module provides

    SparseSectorHamiltonian      a batch of diagonals plus one shared sparse matrix (spin and anyon chains)
    MatrixFreeOperator           any function V -> H V
    DenseSectorOperator          a dense block, e.g. one sector of an Operator (from_operator)
    LinearCombination            sum_i c_i H_i with numbers or one coefficient per sample (combine)
    Deflated                     H + shift sum_j |psi_j><psi_j| (deflate; excited states)

and the solvers

    E0, psi, res = krylov.ground_state(H)                # batched Lanczos with restarts
    E, psi = krylov.lowest_eigenvalues(H, levels=4)      # the lowest levels by deflation
    psi_t   = krylov.evolve(H, psi0, times)              # exp(-i t H) psi0, Krylov steps

Lanczos uses full re-orthogonalisation (memory batch * m * k).
"""
import math
import warnings
from typing import Callable, Optional, Tuple

import torch

from . import cost


class SectorOperator:
    """
    Base class (and protocol) of the operators of one sector known through their action on vectors: matvec(V) for
    V of shape (batch, k) or (batch, k, m), and the attributes batch_size, dim, dtype, device. Every solver of this
    module accepts any object with these attributes. Operators of this class can be added, subtracted and scaled
    (by numbers, or by tensors (batch,) with one coefficient per sample): the results are LinearCombinations.
    """
    batch_size: int
    dim: int
    dtype: torch.dtype
    device: torch.device

    def matvec(self, V: torch.Tensor) -> torch.Tensor:
        raise NotImplementedError

    def __call__(self, V: torch.Tensor) -> torch.Tensor:
        return self.matvec(V)

    def __add__(self, other):
        return combine([(1.0, self), (1.0, other)])

    def __sub__(self, other):
        return combine([(1.0, self), (-1.0, other)])

    def __mul__(self, c):
        return combine([(c, self)])

    __rmul__ = __mul__

    def __neg__(self):
        return combine([(-1.0, self)])

    def to_dense(self) -> torch.Tensor:
        """The matrix (batch, k, k), column by column (k matvecs in one call); for checks on small sectors."""
        cost.check_memory(cost.tensor_bytes((self.batch_size, self.dim, self.dim), self.dtype), self.device,
                          f"dense sector operator ({self.dim} x {self.dim}, batch {self.batch_size})")
        eye = torch.eye(self.dim, dtype=self.dtype, device=self.device).expand(self.batch_size, -1, -1)
        return self.matvec(eye)

    def expectation(self, psi: torch.Tensor) -> torch.Tensor:
        """<psi|H|psi> / <psi|psi>, shape (batch,) (complex for a non-Hermitian H)."""
        Hpsi = self.matvec(psi)
        return _dot(psi.to(Hpsi.dtype), Hpsi) / _dot(psi, psi).real


class SparseSectorHamiltonian(SectorOperator):
    """
    H_b = diag(d_b) + A for b = 1..batch, with d of shape (batch, k) (real) and one sparse k x k matrix A
    shared by the batch (hopping terms). Memory O(batch * k + nnz). from_coo builds one from (row, column, value)
    lists, from_sparse from a torch sparse matrix.
    """

    @classmethod
    def from_coo(cls, rows, cols, values, dim: int, diag=None, dtype: torch.dtype = torch.complex128,
                 device=None) -> 'SparseSectorHamiltonian':
        """H = diag(d_b) + A with A[rows, cols] = values (repeated entries add up); diag (k,) or (batch, k), real;
        the diagonal may also be given through rows == cols."""
        dev = torch.device('cpu') if device is None else torch.device(device)
        rows = torch.as_tensor(rows, dtype=torch.long, device=dev).reshape(-1)
        cols = torch.as_tensor(cols, dtype=torch.long, device=dev).reshape(-1)
        vals = torch.as_tensor(values, device=dev).reshape(-1).to(dtype)
        rdt = torch.empty((), dtype=dtype).real.dtype
        d = torch.zeros(1, dim, dtype=rdt, device=dev) if diag is None else \
            torch.as_tensor(diag, dtype=rdt, device=dev).reshape(-1, dim)
        A = torch.sparse_coo_tensor(torch.stack([rows, cols]), vals, (dim, dim)).coalesce() if rows.numel() else None
        return cls(d, A, dtype)

    @classmethod
    def from_sparse(cls, A: torch.Tensor, diag=None) -> 'SparseSectorHamiltonian':
        """H = diag(d_b) + A for a torch sparse (or dense) k x k matrix A shared by the batch."""
        A = A.to_sparse() if not A.is_sparse else A
        A = A.coalesce()
        i, v = A.indices(), A.values()
        return cls.from_coo(i[0], i[1], v, A.shape[0], diag, dtype=v.dtype, device=v.device)

    def __init__(self, diag: torch.Tensor, offdiag: Optional[torch.Tensor], dtype: torch.dtype):
        self.diag = diag
        self.offdiag = offdiag.coalesce() if offdiag is not None else None
        self.dtype = dtype
        self._csr = {}
        self._diag = {}

    @property
    def batch_size(self) -> int:
        return self.diag.shape[0]

    @property
    def dim(self) -> int:
        return self.diag.shape[1]

    @property
    def device(self):
        return self.diag.device

    def _off(self, dtype):
        if self.offdiag is None:
            return None
        if dtype not in self._csr:
            # CSR is ~20x faster than COO for sparse @ dense on CPU (PyTorch 2.x); silence its beta warning
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', UserWarning)
                self._csr[dtype] = self.offdiag.to(dtype).coalesce().to_sparse_csr()
        return self._csr[dtype]

    def matvec(self, V: torch.Tensor) -> torch.Tensor:
        """(batch, k) or (batch, k, m) -> same shape."""
        squeeze = V.dim() == 2
        if squeeze:
            V = V.unsqueeze(-1)
        # never cast H to the type of V: a real V with a complex H would silently drop Im H
        V = V.to(torch.promote_types(V.dtype, self.dtype))
        d = self._diag.get(V.dtype)
        if d is None:                                   # cast once per dtype
            d = self._diag[V.dtype] = self.diag.to(V.dtype).unsqueeze(-1)
        out = d * V
        A = self._off(V.dtype)
        if A is not None:
            B, k, m = V.shape
            AV = (A @ V.permute(1, 0, 2).reshape(k, B * m)).reshape(k, B, m).permute(1, 0, 2)
            out = out + AV
        return out.squeeze(-1) if squeeze else out

    def to_dense(self) -> torch.Tensor:
        cost.check_memory(cost.tensor_bytes((self.batch_size, self.dim, self.dim), self.dtype), self.device,
                          f"dense Hamiltonian ({self.dim} x {self.dim}, batch {self.batch_size})")
        M = torch.diag_embed(self.diag.to(self.dtype))
        if self.offdiag is not None:
            M = M + self.offdiag.to(self.dtype).to_dense()
        return M


class MatrixFreeOperator(SectorOperator):
    """
    A sector operator given by a function: matvec(V) = fn(V) for V of shape (batch, k). With batched_columns=True
    fn takes blocks of vectors (batch, k, m) instead (also for a single vector, as (batch, k, 1)); otherwise blocks
    are applied column by column. hermitian: whether H is self-adjoint (Lanczos and the time evolution assume it; it
    is not checked).
    """

    def __init__(self, fn: Callable[[torch.Tensor], torch.Tensor], dim: int, batch_size: int = 1,
                 dtype: torch.dtype = torch.complex128, device=None, hermitian: bool = True,
                 batched_columns: bool = False):
        self.fn, self.dim, self.batch_size, self.dtype = fn, int(dim), int(batch_size), dtype
        self.device = torch.device('cpu') if device is None else torch.device(device)
        self.hermitian, self.batched_columns = hermitian, batched_columns

    def matvec(self, V: torch.Tensor) -> torch.Tensor:
        if self.batched_columns:
            return self.fn(V.unsqueeze(-1)).squeeze(-1) if V.dim() == 2 else self.fn(V)
        if V.dim() == 2:
            return self.fn(V)
        return torch.stack([self.fn(V[..., j]) for j in range(V.shape[-1])], dim=-1)


class DenseSectorOperator(SectorOperator):
    """A dense block (batch, k, k) as a sector operator (e.g. a sector of an Operator: from_operator)."""

    def __init__(self, block: torch.Tensor):
        self.block = block if block.dim() == 3 else block.unsqueeze(0)
        self.batch_size, self.dim = self.block.shape[0], self.block.shape[-1]
        self.dtype, self.device = self.block.dtype, self.block.device

    def matvec(self, V: torch.Tensor) -> torch.Tensor:
        dtype = torch.promote_types(V.dtype, self.dtype)
        A = self.block.to(dtype)
        return (A @ V.to(dtype).unsqueeze(-1)).squeeze(-1) if V.dim() == 2 else A @ V.to(dtype)

    def to_dense(self) -> torch.Tensor:
        return self.block


def from_operator(op, sector: int) -> DenseSectorOperator:
    """Sector `sector` of an Operator of a TypeIAlgebra as a sector operator (a view of its block)."""
    return DenseSectorOperator(op.block(sector))


class LinearCombination(SectorOperator):
    """sum_i c_i H_i for sector operators of the same dimension; c_i numbers or tensors (batch,) (one coefficient
    per sample, e.g. couplings to be optimised or a disorder average). Built by `combine` and by +, -, *."""

    def __init__(self, terms):
        terms = list(terms)
        if not terms:
            raise ValueError("LinearCombination: need at least one term")
        dims = {H.dim for _, H in terms}
        if len(dims) != 1:
            raise ValueError(f"LinearCombination: all terms must have the same dimension, got {sorted(dims)}")
        self.terms = terms
        self.dim = dims.pop()
        self.device = terms[0][1].device
        bs = [H.batch_size for _, H in terms] + [c.shape[0] for c, _ in terms if isinstance(c, torch.Tensor) and c.dim()]
        self.batch_size = max(bs)
        dtype = terms[0][1].dtype
        for c, H in terms:
            dtype = torch.promote_types(dtype, H.dtype)
            if isinstance(c, complex) or (isinstance(c, torch.Tensor) and torch.is_complex(c)):
                dtype = torch.promote_types(dtype, torch.complex64)
        self.dtype = dtype

    def matvec(self, V: torch.Tensor) -> torch.Tensor:
        out = None
        for c, H in self.terms:
            y = H.matvec(V)
            if isinstance(c, torch.Tensor) and c.dim():
                c = c.to(device=y.device).reshape(-1, *([1] * (y.dim() - 1)))
            y = c * y
            out = y if out is None else out + y
        return out


def combine(terms) -> LinearCombination:
    """sum_i c_i H_i for [(c_i, H_i), ...]; nested combinations are flattened."""
    flat = []
    for c, H in terms:
        if isinstance(H, LinearCombination):
            flat += [(c * ci, Hi) for ci, Hi in H.terms]
        else:
            flat.append((c, H))
    return LinearCombination(flat)


class Deflated(SectorOperator):
    """H + shift * sum_j |psi_j><psi_j| (the vectors normalised, (batch, k) each): moves the given eigenvectors above
    the rest of the spectrum, so that ground_state finds the next level (excited states)."""

    def __init__(self, H, vectors, shift: float):
        self.H, self.vectors, self.shift = H, list(vectors), shift
        self.batch_size, self.dim, self.dtype, self.device = H.batch_size, H.dim, H.dtype, H.device

    def matvec(self, V):
        out = self.H.matvec(V)
        for psi in self.vectors:
            p = psi.to(out.dtype)
            if V.dim() == 2:
                out = out + self.shift * (p.conj() * V).sum(-1, keepdim=True) * p
            else:
                out = out + self.shift * torch.einsum('bk,bkm->bm', p.conj(), V).unsqueeze(1) * p.unsqueeze(-1)
        return out


def deflate(H, vectors, shift: float) -> Deflated:
    return Deflated(H, vectors, shift)


def _spectral_bound(H, m: int = 20, generator=None) -> torch.Tensor:
    """An upper bound of |H| per batch element from a short Lanczos run (largest |Ritz value| plus the residual)."""
    k = H.dim
    v = torch.randn(H.batch_size, k, generator=generator, dtype=torch.float64).to(device=H.device, dtype=H.dtype)
    alpha, beta, Q = lanczos(H.matvec, v, min(m, k))
    theta = torch.linalg.eigvalsh(_ritz_matrix(alpha, beta, Q))
    return theta.abs().amax(-1) + (beta.abs().amax(-1) if beta.shape[-1] else 0)


def lowest_eigenvalues(H, levels: int = 1, m: int = 80, tol: Optional[float] = None, max_restarts: int = 20,
                       generator: Optional[torch.Generator] = None):
    """
    The `levels` lowest eigenvalues and eigenvectors of every operator of the batch, by restarted Lanczos and
    deflation (each level found as the ground state of H + shift sum_{found} |psi><psi|, shift above the spectrum).
    Degenerate levels are returned with their multiplicity. Returns E (batch, levels) ascending and psi
    (batch, levels, k). For sectors that fit in memory a dense eigensolver is exact and often faster.
    """
    levels = min(int(levels), H.dim)
    shift = None
    Es, vecs = [], []
    for _ in range(levels):
        op = H if not vecs else Deflated(H, vecs, shift)
        E, psi, _ = ground_state(op, m=m, tol=tol, max_restarts=max_restarts, generator=generator)
        if shift is None:
            shift = float(2 * _spectral_bound(H, generator=generator).max()) + 1.0
        Es.append(E)
        vecs.append(psi)
    return torch.stack(Es, dim=-1), torch.stack(vecs, dim=1)


def _dot(a, b):
    return (a.conj() * b).sum(-1)


def _norm(x: torch.Tensor) -> torch.Tensor:
    """
    Euclidean norm over the last axis without overflow or underflow of the sum of squares (float32 overflows for
    entries above ~1.8e19): the vector is divided by its largest entry first. 0 for a zero vector.
    """
    scale = cost.abs_max(x, -1)
    safe = torch.where(scale > 0, scale.clamp(min=torch.finfo(scale.dtype).tiny), torch.ones_like(scale))
    return torch.linalg.vector_norm(x * safe.reciprocal().unsqueeze(-1).to(x.dtype), dim=-1) * scale


def lanczos(matvec: Callable, v0: torch.Tensor, m: int) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    m-step Lanczos with full re-orthogonalisation, batched. v0: (batch, k), normalised inside.
    Returns alpha (batch, m), beta (batch, m-1) (real) and the Krylov basis Q (batch, m, k).
    Breakdown (an invariant subspace) is handled by zero beta and zero further vectors.
    """
    B, k = v0.shape
    q = v0 * torch.clamp(_norm(v0), min=1e-300).reciprocal().unsqueeze(-1).to(v0.dtype)   # v0 = 0 stays 0
    w_first = matvec(q)
    dtype = torch.promote_types(v0.dtype, w_first.dtype)     # a complex H makes the basis complex
    q = q.to(dtype)
    cost.check_memory(1.2 * cost.tensor_bytes((B, m, k), dtype), v0.device, f"Lanczos basis (batch {B}, m={m}, k={k})")
    Q = torch.zeros(B, m, k, dtype=dtype, device=v0.device)
    alpha = torch.zeros(B, m, dtype=q.real.dtype, device=v0.device)
    beta = torch.zeros(B, max(m - 1, 0), dtype=q.real.dtype, device=v0.device)
    for j in range(m):
        Q[:, j] = q
        w = w_first.to(dtype) if j == 0 else matvec(q)
        scale = _norm(w)
        alpha[:, j] = _dot(q, w).real
        # full re-orthogonalisation (twice is enough)
        # (the projection conj(Q) w is computed as conj(Q conj(w)): no conjugated copy of the basis)
        Qj = Q[:, :j + 1]
        for _ in range(2):
            coef = torch.einsum('bjk,bk->bj', Qj, w.conj()).conj() if torch.is_complex(w) else \
                torch.einsum('bjk,bk->bj', Qj, w)
            w = w - torch.einsum('bj,bjk->bk', coef, Qj)
        if j == m - 1:
            break
        b = _norm(w)
        # breakdown: the remainder is at the round-off level of |Hq| (an invariant subspace was reached);
        # the coupling is then exactly zero and all further vectors are zero. The threshold is tied to the
        # machine epsilon: a fixed relative threshold would cut off genuine couplings when H carries a large
        # constant shift (|Hq| ~ |shift|, while beta does not depend on it)
        ok = b > 100 * torch.finfo(b.dtype).eps * torch.clamp(scale, min=1e-300)
        beta[:, j] = torch.where(ok, b, torch.zeros_like(b))
        q = torch.where(ok[:, None], w * torch.clamp(b, min=1e-300).reciprocal()[:, None].to(w.dtype), torch.zeros_like(w))
    return alpha, beta, Q


def _ritz_matrix(alpha, beta, Q):
    """
    Tridiagonal matrix of a Lanczos run in which the zero vectors after a breakdown are moved above
    the spectrum (Gershgorin bound), so that they cannot produce spurious low Ritz values.
    """
    alive = Q.abs().amax(dim=-1) > 0                                     # (batch, m)
    bound = (alpha.abs().amax(-1) + 2 * (beta.abs().amax(-1) if beta.shape[-1] else 0) + 1.0)
    alpha = torch.where(alive, alpha, bound[:, None].expand_as(alpha))
    return _tridiag(alpha, beta)


def _tridiag(alpha, beta):
    T = torch.diag_embed(alpha)
    if beta.shape[-1]:
        T = T + torch.diag_embed(beta, 1) + torch.diag_embed(beta, -1)
    return T


def ground_state(H, m: int = 80, tol: Optional[float] = None, max_restarts: int = 20,
                 v0: Optional[torch.Tensor] = None, generator: Optional[torch.Generator] = None):
    """
    Lowest eigenvalue and eigenvector of every Hamiltonian in the batch by restarted Lanczos.
    H: object with matvec, batch_size, dim, dtype, device (e.g. SparseSectorHamiltonian).
    Stops when the residual ||H psi - E psi|| < tol * max(1, |E|) for all batch elements
    (default tol: 1e-10 in double, 1e-5 in single precision).
    Returns E0 (batch,), psi (batch, k), residual (batch,).
    """
    B, k = H.batch_size, H.dim
    m = min(m, k)
    if tol is None:
        tol = 1e-10 if H.dtype in (torch.float64, torch.complex128) else 1e-5
    if v0 is None:
        v = torch.randn(B, k, generator=generator, dtype=torch.float64).to(device=H.device, dtype=H.dtype)
    else:                                                  # never lose Im H to a real start vector
        v = v0.to(device=H.device, dtype=torch.promote_types(v0.dtype, H.dtype))
    for it in range(max_restarts):
        alpha, beta, Q = lanczos(H.matvec, v, m)
        w, Y = torch.linalg.eigh(_ritz_matrix(alpha, beta, Q))
        psi = torch.einsum('bj,bjk->bk', Y[:, :, 0].to(Q.dtype), Q)
        psi = psi * _norm(psi).reciprocal().unsqueeze(-1).to(psi.dtype)
        E = w[:, 0]
        res = _norm(H.matvec(psi) - E[:, None].to(psi.dtype) * psi)
        if torch.all(res < tol * torch.clamp(E.abs(), min=1.0)):
            return E, psi, res
        v = psi
    warnings.warn(f"ground_state: residual {res.max().item():.1e} after {max_restarts} restarts "
                  f"(increase m or max_restarts)", cost.CostWarning, stacklevel=2)
    return E, psi, res


def _tail_over_x(x: torch.Tensor, m: int) -> torch.Tensor:
    """
    S_m(x) / x for x >= 0, with S_m(x) = sum_{j >= m} x^j / j! = e^x P(m, x) (P the regularised lower incomplete
    gamma function); increasing in x. Limit x -> 0: 1 for m = 1, 0 for m >= 2.
    """
    xs = x.to(torch.float64).clamp(min=1e-300)
    P = torch.special.gammainc(torch.full_like(xs, float(m)), xs)
    return torch.exp(xs + torch.log(P) - torch.log(xs))


def krylov_error_bound(beta_m: torch.Tensor, ritz: torch.Tensor, dt: float, m: int) -> torch.Tensor:
    """
    Bound (Theorem 9.18(b) of the lecture notes) for || exp(-i dt H) v - ||v|| Q_m exp(-i dt T_m) e_1 || / ||v|| from
    the data of one m-step Lanczos run: the residual norm beta_m = ||H q_m - Q_m T_m e_m|| (batch,) and the Ritz
    values (batch, m) (eigenvalues of T_m). With w_T the half-width of the Ritz values and x = |dt| w_T / 2 it is
    (4 beta_m / w_T) S_m(x) = 2 beta_m |dt| S_m(x) / x, S_m(x) = sum_{j >= m} x^j / j!. It holds in exact arithmetic
    and needs no estimate of the spectrum of H; it is invariant under H -> H + c.
    """
    wT = ((ritz.amax(-1) - ritz.amin(-1)) / 2).to(torch.float64)
    x = abs(dt) * wT / 2
    return 2 * beta_m.to(torch.float64) * abs(dt) * _tail_over_x(x, m)


def _max_step(beta_m: torch.Tensor, ritz: torch.Tensor, m: int, rate: float) -> torch.Tensor:
    """
    Largest h with krylov_error_bound(beta_m, ritz, h, m) <= rate * h, per batch element (inf if beta_m = 0: the
    Krylov space is invariant and every step is exact). The bound divided by h is increasing in h, so bisection.
    """
    wT = ((ritz.amax(-1) - ritz.amin(-1)) / 2).to(torch.float64)
    b = beta_m.to(torch.float64)
    target = rate / (2 * b.clamp(min=1e-300))            # S_m(x) / x <= rate / (2 beta_m)
    lo = torch.zeros_like(b)
    hi = torch.full_like(b, 4.0 * (m + 1))                # S_m(x)/x >= x^(m-1)/m! is far too large there
    for _ in range(80):
        mid = (lo + hi) / 2
        ok = _tail_over_x(mid, m) <= target
        lo, hi = torch.where(ok, mid, lo), torch.where(ok, hi, mid)
    h = 2 * lo / wT.clamp(min=1e-300)
    exact = (b == 0) | (wT == 0) & (m >= 2)
    return torch.where(exact, torch.full_like(h, math.inf), h)


def _krylov_step(H, psi: torch.Tensor, dt: float, m: int, rate: Optional[float] = None):
    """
    One Krylov step from psi. With rate=None the step is dt; otherwise it is the largest h <= dt whose error bound is
    <= rate * h * ||psi|| for every batch element. Returns (state, error bound (batch,), step taken).
    """
    nrm = _norm(psi)
    alpha, beta, Q = lanczos(H.matvec, psi, m)
    # residual r_m = H q_m - alpha_m q_m - beta_{m-1} q_{m-1} of the Lanczos relation H Q_m = Q_m T_m + r_m e_m^T
    r = H.matvec(Q[:, -1]) - alpha[:, -1, None].to(Q.dtype) * Q[:, -1]
    if m > 1:
        r = r - beta[:, -1, None].to(Q.dtype) * Q[:, -2]
    beta_m = _norm(r)
    # after a breakdown the zero vectors carry alpha = 0, which is not a Ritz value of the live part (whose
    # Krylov space is then invariant and the step exact, beta_m = 0); alpha_0 lies inside the live Ritz range
    alive = Q.abs().amax(dim=-1) > 0
    theta, S = torch.linalg.eigh(_tridiag(torch.where(alive, alpha, alpha[:, :1]), beta))
    if rate is not None:
        h_max = _max_step(beta_m, theta, m, rate).min().item()
        dt = dt if h_max >= dt else h_max
    # exp(-i dt T) e_1 from the eigendecomposition of T, relative to the centre c of the Ritz values, times the phase
    # exp(-i dt c): exact shift invariance (H -> H + c changes only the phase)
    c = (theta.amax(-1) + theta.amin(-1)) / 2
    cdt = torch.promote_types(psi.dtype, torch.complex64)
    ph = torch.exp(-1j * dt * (theta - c[:, None]).to(torch.float64)).to(cdt)
    coef = torch.einsum('bij,bj,bj->bi', S.to(cdt), ph, S[:, 0, :].to(cdt))      # exp(-i dt (T - c)) e_1
    # the zero vectors (after a breakdown) contribute nothing, whatever their coefficients
    out = torch.einsum('bj,bjk->bk', coef, Q.to(cdt)) * (torch.exp(-1j * dt * c.to(torch.float64)).to(cdt)
                                                          * nrm.to(cdt))[:, None]
    err = krylov_error_bound(beta_m, theta, dt, m) * nrm.to(torch.float64)
    return out, err, dt


def expm_multiply_step(H, psi: torch.Tensor, dt: float, m: int = 30) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    exp(-i dt H) psi in an m-dimensional Krylov space, w = ||psi|| Q_m exp(-i dt T_m) e_1. Returns the new state
    and a bound for its error per batch element, ||psi|| (4 beta_m / w_T) S_m(|dt| w_T / 2) (krylov_error_bound;
    Theorem 9.18 of the lecture notes): it holds in exact arithmetic, uses only the Lanczos data (beta_m the norm of
    the Lanczos residual, w_T the half-width of the Ritz values) and needs no estimate of the spectrum of H.
    (Earlier versions returned beta_m |e_m^T exp(-i dt T_m) e_1|, the norm of the defect of the Schroedinger
    equation at time dt: an error rate, not an error.) Rounding errors are not included.
    """
    out, err, _ = _krylov_step(H, psi, dt, m)
    return out, err


def evolve(H, psi0: torch.Tensor, times, m: int = 30, dt: Optional[float] = None, tol: float = 1e-8,
           progress: bool = False) -> torch.Tensor:
    """
    psi(t) = exp(-i t H) psi0 at every time of the increasing grid `times` (starting at times[0]), in Krylov steps.

    Default (dt=None): every step from a state psi is the largest h (up to the next grid time) whose error bound
    (krylov_error_bound) is <= tol * h / (times[-1] - times[0]) * ||psi0||. As every step is an isometry, the local
    errors add up (Corollary 9.19 of the lecture notes): the Krylov error of the state at time t is at most
    tol * ||psi0|| * (t - times[0]) / (times[-1] - times[0]), in exact arithmetic (rounding errors, which grow with the
    number of steps, are not included). The step depends on H only through the Lanczos data, so a constant shift
    H -> H + c does not change it. With a fixed dt the steps are dt (shortened to fit the grid), and a CostWarning
    is issued if the sum of the error bounds exceeds tol * ||psi0||. Returns (T, batch, k), complex.
    """
    if not tol > 0:
        raise ValueError(f"krylov.evolve: tol must be positive, got {tol}")
    ts = [float(t) for t in torch.as_tensor(times, dtype=torch.float64).reshape(-1).tolist()]   # float64: no float32 rounding of t
    psi = psi0.to(torch.promote_types(psi0.dtype, torch.complex64))
    m = min(max(m, 2), H.dim)
    nrm0 = _norm(psi).to(torch.float64)
    span = ts[-1] - ts[0]
    out, err_total = [psi], torch.zeros(psi.shape[0], dtype=torch.float64, device=psi.device)
    timer = None
    if dt is not None:
        n_steps = sum(max(1, math.ceil((b - a) / dt - 1e-9)) for a, b in zip(ts[:-1], ts[1:]))
        timer = cost.StepTimer(n_steps, "krylov.evolve", progress)
    rate = tol / span if span > 0 else math.inf
    for a, b in zip(ts[:-1], ts[1:]):
        if dt is not None:
            n = max(1, math.ceil((b - a) / dt - 1e-9))
            for _ in range(n):
                psi, err, _ = _krylov_step(H, psi, (b - a) / n, m)
                err_total = err_total + err.to(err_total.dtype)
                timer.step(psi.device)
        else:
            t = a
            while b - t > 1e-12 * max(abs(a), abs(b), b - a):
                psi, err, h = _krylov_step(H, psi, b - t, m, rate)
                if not t + h > t:
                    raise ValueError(f"krylov.evolve: the step for tol = {tol:.0e} is {h:.1e} at t = {t:.3e}; "
                                     f"increase tol or m")
                err_total = err_total + err.to(err_total.dtype)
                t = b if h >= b - t else t + h
                if timer is None:                    # the number of steps is known only after the first one
                    est = sum(math.ceil((d - c) / h - 1e-9) for c, d in zip(ts[:-1], ts[1:]))
                    timer = cost.StepTimer(max(est, 1), "krylov.evolve", progress)
                timer.step(psi.device)
        out.append(psi)
    if timer is not None:
        timer.close()
    if torch.any(err_total > tol * nrm0 * (1 + 1e-6)):
        warnings.warn(f"krylov.evolve: error bound {err_total.max().item():.1e} > tol {tol:.0e}; use a "
                      f"smaller dt or a larger m", cost.CostWarning, stacklevel=2)
    return torch.stack(out)
