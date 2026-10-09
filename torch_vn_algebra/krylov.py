"""
Krylov-subspace methods for large sectors, where dense matrices no longer fit.

A Hamiltonian is only needed through its action on vectors, `matvec(V)` with V of shape
(batch, k) -> (batch, k); every batch element may be a different operator (e.g. one disorder
realisation each). `SparseSectorHamiltonian` provides such an operator for spin chains: a
batch of diagonals plus one shared sparse off-diagonal part.

    E0, psi = krylov.ground_state(H)                     # batched Lanczos with restarts
    psi_t   = krylov.evolve(H, psi0, times)              # exp(-i t H) psi0, Krylov steps

Lanczos uses full re-orthogonalisation (memory batch * m * k).
"""
import math
import warnings
from typing import Callable, Optional, Tuple

import torch

from . import cost


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
        out = self.diag.to(V.dtype).unsqueeze(-1) * V
        A = self._off(V.dtype)
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


def lanczos(matvec: Callable, v0: torch.Tensor, m: int) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    m-step Lanczos with full re-orthogonalisation, batched. v0: (batch, k), normalised inside.
    Returns alpha (batch, m), beta (batch, m-1) (real) and the Krylov basis Q (batch, m, k).
    Breakdown (an invariant subspace) is handled by zero beta and zero further vectors.
    """
    B, k = v0.shape
    q = v0 / torch.linalg.vector_norm(v0, dim=-1, keepdim=True)
    w_first = matvec(q)
    dtype = torch.promote_types(v0.dtype, w_first.dtype)     # a complex H makes the basis complex
    q = q.to(dtype)
    # the basis Q plus the copy made by the re-orthogonalisation
    cost.check_memory(2 * cost.tensor_bytes((B, m, k), dtype), v0.device, f"Lanczos basis (batch {B}, m={m}, k={k})")
    Q = torch.zeros(B, m, k, dtype=dtype, device=v0.device)
    alpha = torch.zeros(B, m, dtype=q.real.dtype, device=v0.device)
    beta = torch.zeros(B, max(m - 1, 0), dtype=q.real.dtype, device=v0.device)
    for j in range(m):
        Q[:, j] = q
        w = w_first.to(dtype) if j == 0 else matvec(q)
        scale = torch.linalg.vector_norm(w, dim=-1)
        alpha[:, j] = _dot(q, w).real
        # full re-orthogonalisation (twice is enough)
        for _ in range(2):
            w = w - torch.einsum('bj,bjk->bk', torch.einsum('bjk,bk->bj', Q[:, :j + 1].conj(), w), Q[:, :j + 1])
        if j == m - 1:
            break
        b = torch.linalg.vector_norm(w, dim=-1)
        # breakdown: the remainder is round-off relative to |Hq| (an invariant subspace was reached);
        # the coupling is then exactly zero and all further vectors are zero
        ok = b > 1e-10 * torch.clamp(scale, min=1e-300)
        beta[:, j] = torch.where(ok, b, torch.zeros_like(b))
        q = torch.where(ok[:, None], w / torch.clamp(b, min=1e-300)[:, None], torch.zeros_like(w))
    return alpha, beta, Q


def _ritz_matrix(alpha, beta, Q):
    """
    Tridiagonal matrix of a Lanczos run in which the zero vectors after a breakdown are moved above
    the spectrum (Gershgorin bound), so that they cannot produce spurious low Ritz values.
    """
    alive = torch.linalg.vector_norm(Q, dim=-1) > 0                       # (batch, m)
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
        psi = psi / torch.linalg.vector_norm(psi, dim=-1, keepdim=True)
        E = w[:, 0]
        res = torch.linalg.vector_norm(H.matvec(psi) - E[:, None].to(psi.dtype) * psi, dim=-1)
        if torch.all(res < tol * torch.clamp(E.abs(), min=1.0)):
            return E, psi, res
        v = psi
    warnings.warn(f"ground_state: residual {res.max().item():.1e} after {max_restarts} restarts "
                  f"(increase m or max_restarts)", cost.CostWarning, stacklevel=2)
    return E, psi, res


def expm_multiply_step(H, psi: torch.Tensor, dt: float, m: int = 30) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    exp(-i dt H) psi in an m-dimensional Krylov space. Returns the new state and an error estimate
    per batch element, beta_m |[exp(-i dt T)]_{m-1, 0}| (the standard a-posteriori bound).
    """
    nrm = torch.linalg.vector_norm(psi, dim=-1)
    alpha, beta, Q = lanczos(H.matvec, psi, m)
    T = _tridiag(alpha, beta).to(torch.promote_types(psi.dtype, torch.complex64))
    E = torch.linalg.matrix_exp(-1j * dt * T)                           # (B, m, m)
    out = torch.einsum('bj,bjk->bk', E[:, :, 0], Q.to(E.dtype)) * nrm[:, None].to(E.dtype)
    # error estimate with the next (unused) Lanczos coefficient: |h_{m+1,m}| * |e_m^T exp(-i dt T) e_1|
    w = H.matvec(Q[:, -1]) - alpha[:, -1, None].to(Q.dtype) * Q[:, -1]
    if m > 1:
        w = w - beta[:, -1, None].to(Q.dtype) * Q[:, -2]
    err = torch.linalg.vector_norm(w, dim=-1) * E[:, -1, 0].abs() * nrm
    return out, err


def evolve(H, psi0: torch.Tensor, times, m: int = 30, dt: Optional[float] = None, tol: float = 1e-8,
           progress: bool = False) -> torch.Tensor:
    """
    psi(t) = exp(-i t H) psi0 at every time of the increasing grid `times` (starting at times[0]),
    in Krylov steps of size dt (default 0.1 * m / ||H||, with ||H|| estimated by a short Lanczos run). Warns if the accumulated error
    estimate exceeds tol. Returns (T, batch, k), complex.
    """
    ts = [float(t) for t in torch.as_tensor(times).reshape(-1).tolist()]
    psi = psi0.to(torch.promote_types(psi0.dtype, torch.complex64))
    if dt is None:
        alpha, beta, _ = lanczos(H.matvec, psi, min(20, H.dim))
        spread = (torch.linalg.eigvalsh(_tridiag(alpha, beta)).abs().max()).item()
        dt = max(1e-3, 0.1 * m / max(spread, 1e-12))
    n_steps = sum(max(1, math.ceil((b - a) / dt - 1e-9)) for a, b in zip(ts[:-1], ts[1:]))
    timer = cost.StepTimer(n_steps, "krylov.evolve", progress)
    out, err_total = [psi], torch.zeros(psi.shape[0], device=psi.device)
    for a, b in zip(ts[:-1], ts[1:]):
        n = max(1, math.ceil((b - a) / dt - 1e-9))
        h = (b - a) / n
        for _ in range(n):
            psi, err = expm_multiply_step(H, psi, h, min(m, H.dim))
            err_total = err_total + err.to(err_total.dtype)
            timer.step(psi.device)
        out.append(psi)
    timer.close()
    if torch.any(err_total > tol):
        warnings.warn(f"krylov.evolve: error estimate {err_total.max().item():.1e} > tol {tol:.0e}; use a "
                      f"smaller dt or a larger m", cost.CostWarning, stacklevel=2)
    return torch.stack(out)
