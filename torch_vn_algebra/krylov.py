"""
Krylov-subspace methods for large sectors, where dense matrices no longer fit.

A Hamiltonian is only needed through its action on vectors, `matvec(V)` with V of shape
(batch, k) -> (batch, k); every batch element may be a different operator (e.g. one disorder
realisation each). `SparseSectorHamiltonian` provides such an operator for spin chains: a
batch of diagonals plus one shared sparse off-diagonal part.

    E0, psi, res = krylov.ground_state(H)                # batched Lanczos with restarts
    psi_t   = krylov.evolve(H, psi0, times)              # exp(-i t H) psi0, Krylov steps

Lanczos uses full re-orthogonalisation (memory batch * m * k).
"""
import math
import warnings
from typing import Callable, Optional, Tuple

import torch

from . import cost


def _version(x):
    try:
        return id(x), x._version
    except RuntimeError:                         # inference tensors have no version counter
        return None


class SparseSectorHamiltonian:
    """
    H_b = diag(d_b) + A for b = 1..batch, with d of shape (batch, k) (real) and one sparse k x k matrix A
    shared by the batch (hopping terms). Memory O(batch * k + nnz).
    """

    def __init__(self, diag: torch.Tensor, offdiag: Optional[torch.Tensor], dtype: torch.dtype):
        self.diag = diag
        self.offdiag = offdiag.coalesce() if offdiag is not None else None
        self.dtype = dtype
        self._csr = {}
        self._centered_parts = {}
        self._abs_csr = {}
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
        version = _version(self.offdiag)
        cacheable = version is not None and not (torch.is_grad_enabled() and self.offdiag.requires_grad)
        cached = self._csr.get(dtype)
        if not cacheable or cached is None or cached[0] != version:
            # CSR is ~20x faster than COO for sparse @ dense on CPU (PyTorch 2.x); silence its beta warning
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', UserWarning)
                out = self.offdiag.to(dtype).coalesce().to_sparse_csr()
            if cacheable:
                self._csr[dtype] = version, out
            return out
        return cached[1]

    def _diagonal(self, dtype):
        version = _version(self.diag)
        cacheable = version is not None and not (torch.is_grad_enabled() and self.diag.requires_grad)
        cached = self._diag.get(dtype)
        if not cacheable or cached is None or cached[0] != version:
            out = self.diag.to(dtype)
            if cacheable:
                self._diag[dtype] = version, out
            return out
        return cached[1]

    def _split_sparse_diagonal(self, dtype):
        """Shared sparse diagonal and strictly off-diagonal CSR, with the same mutation rules as _off."""
        if self.offdiag is None:
            return None, None
        version = _version(self.offdiag)
        cacheable = version is not None and not (torch.is_grad_enabled() and self.offdiag.requires_grad)
        cached = self._centered_parts.get(dtype)
        if not cacheable or cached is None or cached[0] != version:
            indices = self.offdiag.indices()
            values = self.offdiag.values().to(dtype)
            diagonal_mask = indices[0] == indices[1]
            diagonal = torch.zeros(self.dim, dtype=dtype, device=self.device).index_add(
                0, indices[0, diagonal_mask], values[diagonal_mask])
            if diagonal_mask.any():
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', UserWarning)
                    off = torch.sparse_coo_tensor(indices[:, ~diagonal_mask], values[~diagonal_mask],
                                                  self.offdiag.shape).coalesce().to_sparse_csr()
            else:
                off = self._off(dtype)
            if cacheable:
                self._centered_parts[dtype] = version, diagonal, off
            return diagonal, off
        return cached[1], cached[2]

    def _total_diagonal(self, dtype):
        diagonal, _ = self._split_sparse_diagonal(dtype)
        return self._diagonal(dtype) if diagonal is None else self._diagonal(dtype) + diagonal

    def matvec(self, V: torch.Tensor) -> torch.Tensor:
        """(batch, k) or (batch, k, m) -> same shape."""
        return self._matvec(V)

    def matvec_shifted(self, V: torch.Tensor, shift: torch.Tensor) -> torch.Tensor:
        """(H - shift_b I) V, subtracting the shift before multiplying the diagonal."""
        return self._matvec(V, shift)

    def _centered_action_scale(self, V, shift):
        """|| |H - shift I| |V| ||, the sparse action's local roundoff scale, not a global norm."""
        with torch.no_grad():
            v = V.detach().abs()
            diagonal = self._total_diagonal(V.dtype).real - shift[:, None]
            out = diagonal.abs() * v
            if self.offdiag is not None:
                version = _version(self.offdiag)
                cached = self._abs_csr.get(V.dtype)
                if version is None or cached is None or cached[0] != version:
                    _, off = self._split_sparse_diagonal(V.dtype)
                    absolute = off.detach().abs()
                    if version is not None:
                        self._abs_csr[V.dtype] = version, absolute
                else:
                    absolute = cached[1]
                out = out + (absolute @ v.T).T
            return _norm(out)

    def _matvec(self, V, shift=None):
        squeeze = V.dim() == 2
        if squeeze:
            V = V.unsqueeze(-1)
        # never cast H to the type of V: a real V with a complex H would silently drop Im H
        V = V.to(torch.promote_types(V.dtype, self.dtype))
        d = self._diagonal(V.dtype)
        if shift is not None:
            diagonal, A = self._split_sparse_diagonal(V.dtype)
            d = d if diagonal is None else d + diagonal
            d = d - shift.to(V.dtype)[:, None]
        else:
            A = self._off(V.dtype)
        out = d.unsqueeze(-1) * V
        if A is not None:
            B, k, m = V.shape
            AV = (A @ V.permute(1, 0, 2).reshape(k, B * m)).reshape(k, B, m).permute(1, 0, 2)
            out = out + AV
        return out.squeeze(-1) if squeeze else out

    __call__ = matvec

    def to_dense(self) -> torch.Tensor:
        cost.check_memory(cost.tensor_bytes((self.batch_size, self.dim, self.dim), self.dtype), self.device,
                          f"dense Hamiltonian ({self.dim} x {self.dim}, batch {self.batch_size})")
        M = torch.diag_embed(self.diag.to(self.dtype))
        if self.offdiag is not None:
            M = M + self.offdiag.to(self.dtype).to_dense()
        return M


def _dot(a, b):
    return (a.conj() * b).sum(-1)


def _norm(x: torch.Tensor) -> torch.Tensor:
    """
    Euclidean norm over the last axis without overflow or underflow of the sum of squares (float32 overflows for
    entries above ~1.8e19): the vector is divided by its largest entry first. 0 for a zero vector.
    """
    # The scale cancels by homogeneity. Detaching it avoids overflow in reciprocal/abs gradients
    # at extreme magnitudes while preserving the Euclidean norm's derivatives through x.
    scale = cost.abs_max(x, -1).detach()
    safe = torch.where(scale > 0, scale.clamp(min=torch.finfo(scale.dtype).tiny), torch.ones_like(scale))
    return torch.linalg.vector_norm(x * safe.reciprocal().unsqueeze(-1).to(x.dtype), dim=-1) * safe


def lanczos(matvec: Callable, v0: torch.Tensor, m: int) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    m-step Lanczos with full re-orthogonalisation, batched. v0: (batch, k), normalised inside.
    Returns alpha (batch, m), beta (batch, m-1) (real) and the Krylov basis Q (batch, m, k).
    Breakdown (an invariant subspace) is handled by zero beta and zero further vectors.
    Internally a scalar shift is kept separate; the returned alpha still belongs to the original operator.
    """
    alpha, beta, Q, shift = _lanczos_centered(matvec, v0, m)
    alive = Q.abs().amax(dim=-1) > 0
    return torch.where(alive, alpha + shift[:, None], alpha), beta, Q


def _centered_action(matvec, v, shift):
    owner = getattr(matvec, '__self__', None)
    shifted = getattr(owner, 'matvec_shifted', None)
    if shifted is not None:
        return shifted(v, shift)
    return matvec(v) - shift[:, None].to(v.dtype) * v


def _lanczos_centered(matvec, v0, m):
    """Keep the scalar shift separate from Ritz data, including at large spectral offsets."""
    B, k = v0.shape
    owner = getattr(matvec, '__self__', None)
    sparse = isinstance(owner, SparseSectorHamiltonian)
    dtype = torch.promote_types(v0.dtype, owner.dtype) if sparse else v0.dtype
    start = v0.to(dtype)
    norm = _norm(start)
    safe = torch.where(norm > 0, norm, torch.ones_like(norm))
    q = start * safe.reciprocal().unsqueeze(-1).to(dtype)     # v0 = 0 stays 0 in either precision
    if sparse:
        B = torch.broadcast_shapes((B,), (owner.batch_size,))[0]
        q = q.expand(B, k)
        # An occupied diagonal entry avoids estimating a huge scalar shift from rounded Hq, and an
        # unvisited spectator cannot set the shift. Include any diagonal stored in the shared sparse A.
        pivot = q.abs().argmax(dim=-1, keepdim=True)
        diagonal = owner._total_diagonal(dtype).real.expand(B, k)
        shift = diagonal.gather(-1, pivot).squeeze(-1)
        w_first = owner.matvec_shifted(q, shift)
    else:
        w_first = matvec(q)
        dtype = torch.promote_types(dtype, w_first.dtype)    # a complex H makes the basis complex
        q = q.to(dtype)
        if dtype != v0.dtype:
            norm = _norm(q)
            q = q / torch.where(norm > 0, norm, torch.ones_like(norm))[:, None]
            w_first = matvec(q)
        B = max(q.shape[0], w_first.shape[0])
        q = q.expand(B, k)
        norm_sq = _dot(q, q).real
        shift = _dot(q, w_first.to(dtype)).real / torch.where(norm_sq > 0, norm_sq, torch.ones_like(norm_sq))
        w_first = w_first.to(dtype) - shift[:, None].to(dtype) * q
    cost.check_memory(1.2 * cost.tensor_bytes((B, m, k), dtype), v0.device, f"Lanczos basis (batch {B}, m={m}, k={k})")
    Q = torch.zeros(B, m, k, dtype=dtype, device=v0.device)
    alpha = torch.zeros(B, m, dtype=q.real.dtype, device=v0.device)
    beta = torch.zeros(B, max(m - 1, 0), dtype=q.real.dtype, device=v0.device)
    for j in range(m):
        Q[:, j] = q
        w = w_first.to(dtype) if j == 0 else _centered_action(matvec, q, shift)
        alpha[:, j] = _dot(q, w).real
        scale = _norm(w)
        if sparse:
            scale = torch.maximum(scale, owner._centered_action_scale(q, shift))
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
        # Use an absolute-action error scale: it detects cancellation in a numerical eigenvector but does
        # not count an unvisited, decoupled high-energy spectator. Never scale the test by a removable
        # scalar shift. Any discarded defect remains in the evolution error bound.
        ok = b > 100 * torch.finfo(b.dtype).eps * torch.clamp(scale, min=torch.finfo(b.dtype).tiny)
        if j + 1 >= k:                           # the full Hilbert space is already spanned
            ok = torch.zeros_like(ok)
        beta[:, j] = torch.where(ok, b, torch.zeros_like(b))
        safe = torch.where(ok, b, torch.ones_like(b))
        q = torch.where(ok[:, None], w * safe.reciprocal()[:, None].to(w.dtype), torch.zeros_like(w))
    return alpha, beta, Q, shift


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
        alpha, beta, Q, shift = _lanczos_centered(H.matvec, v, m)
        w, Y = torch.linalg.eigh(_ritz_matrix(alpha, beta, Q))
        psi = torch.einsum('bj,bjk->bk', Y[:, :, 0].to(Q.dtype), Q)
        psi = psi * _norm(psi).reciprocal().unsqueeze(-1).to(psi.dtype)
        E = w[:, 0] + shift
        res = _norm(_centered_action(H.matvec, psi, shift) - w[:, :1].to(psi.dtype) * psi)
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
    order = torch.as_tensor(m, dtype=xs.dtype, device=xs.device)
    P = torch.special.gammainc(order, xs)
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
    order = torch.as_tensor(m, dtype=b.dtype, device=b.device)
    target = rate / (2 * b.clamp(min=1e-300))            # S_m(x) / x <= rate / (2 beta_m)
    lo = torch.zeros_like(b)
    hi = torch.ones_like(b) * (4.0 * (order + 1))          # S_m(x)/x >= x^(m-1)/m! is far too large there
    for _ in range(80):
        mid = (lo + hi) / 2
        ok = _tail_over_x(mid, m) <= target
        lo, hi = torch.where(ok, mid, lo), torch.where(ok, hi, mid)
    h = 2 * lo / wT.clamp(min=1e-300)
    exact = (b == 0) | (wT == 0) & (order >= 2)
    return torch.where(exact, torch.full_like(h, math.inf), h)


def _krylov_step(H, psi: torch.Tensor, dt: float, m: int, rate: Optional[float] = None):
    """
    One Krylov step from psi. With rate=None the step is dt; otherwise it is the largest h <= dt whose error bound is
    <= rate * h * ||psi|| for every batch element. Returns (state, error bound (batch,), step taken).
    """
    nrm = _norm(psi)
    alpha, beta, Q, shift = _lanczos_centered(H.matvec, psi, m)
    # Numerical breakdown can leave a small nonzero defect. Measure it at the last LIVE vector; using the
    # trailing zero vector would silently turn that truncation error into an exact-zero certificate.
    alive = Q.abs().amax(dim=-1) > 0
    count = alive.sum(-1)
    last = (count - 1).clamp(min=0)
    batch = torch.arange(Q.shape[0], device=Q.device)
    qlast = Q[batch, last]
    r = _centered_action(H.matvec, qlast, shift) - alpha[batch, last, None].to(Q.dtype) * qlast
    if m > 1:
        prev = (last - 1).clamp(min=0)
        coupling = torch.where(last > 0, beta[batch, prev], torch.zeros_like(beta[:, 0]))
        r = r - coupling[:, None].to(Q.dtype) * Q[batch, prev]
    beta_m = _norm(r)
    # after a breakdown the zero vectors carry alpha = 0, which is not a Ritz value of the live part (whose
    # Krylov space is then invariant and the step exact, beta_m = 0); alpha_0 lies inside the live Ritz range
    theta, S = torch.linalg.eigh(_tridiag(torch.where(alive, alpha, alpha[:, :1]), beta))
    order = count.clamp(min=1)
    if rate is not None:
        h_max = _max_step(beta_m, theta, order, rate).min().item()
        dt = dt if h_max >= dt else h_max
    # exp(-i dt T) e_1 from the eigendecomposition of T, relative to the centre c of the Ritz values, times the phase
    # exp(-i dt c): exact shift invariance (H -> H + c changes only the phase)
    c = (theta.amax(-1) + theta.amin(-1)) / 2
    cdt = torch.promote_types(Q.dtype, torch.complex64)
    ph = torch.exp(-1j * dt * (theta - c[:, None]).to(torch.float64)).to(cdt)
    coef = torch.einsum('bij,bj,bj->bi', S.to(cdt), ph, S[:, 0, :].to(cdt))      # exp(-i dt (T - c)) e_1
    # the zero vectors (after a breakdown) contribute nothing, whatever their coefficients
    out = torch.einsum('bj,bjk->bk', coef, Q.to(cdt)) * (torch.exp(-1j * dt * c.to(torch.float64)).to(cdt)
                                                          * torch.exp(-1j * dt * shift.to(torch.float64)).to(cdt)
                                                          * nrm.to(cdt))[:, None]
    err = krylov_error_bound(beta_m, theta, dt, order) * nrm.to(torch.float64)
    return out, err, dt


def expm_multiply_step(H, psi: torch.Tensor, dt: float, m: int = 30) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    exp(-i dt H) psi in an m-dimensional Krylov space, w = ||psi|| Q_m exp(-i dt T_m) e_1. Returns the new state
    and a bound for its error per batch element, ||psi|| (4 beta_m / w_T) S_m(|dt| w_T / 2) (krylov_error_bound;
    Theorem 9.18 of the lecture notes): it holds in exact arithmetic, uses only the Lanczos data (beta_m the norm of
    the Lanczos residual, w_T the half-width of the Ritz values) and needs no estimate of the spectrum of H.
    (Earlier versions returned beta_m |e_m^T exp(-i dt T_m) e_1|, the norm of the defect of the Schroedinger
    equation at time dt: an error rate, not an error.) Rounding errors are not included.
    After numerical breakdown, the bound uses the live Krylov dimension and its remaining residual.
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
    B = torch.broadcast_shapes((psi.shape[0],), (H.batch_size,))[0]
    psi = psi.expand(B, *psi.shape[1:])
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
            width, elapsed = b - a, 0.0
            while elapsed < width:
                remaining = width - elapsed
                psi, err, h = _krylov_step(H, psi, remaining, m, rate)
                if not elapsed + h > elapsed:
                    raise ValueError(f"krylov.evolve: the step for tol = {tol:.0e} is {h:.1e} at t = {a + elapsed:.3e}; "
                                     f"increase tol or m")
                err_total = err_total + err.to(err_total.dtype)
                elapsed = width if h >= remaining else elapsed + h
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
