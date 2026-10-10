"""
Functional calculus on stacks of blocks (n, k, k), at four levels of generality.

1. Self-adjoint X: f(X) = V f(w) V^* from the Hermitian eigensolver; f any function on the real spectrum.
2. Normal X (unitaries, ...): the eigenvectors of the Hermitian matrix H + alpha K, X = H + i K, diagonalise X
   (H and K commute); f on the complex spectrum. Again only the Hermitian eigensolver is used.
3. Singular values of any X from the Hermitian eigensolver of the dilation [[0, X], [X^*, 0]], whose
   eigenvalues are +-s_i: |X|, the pseudo-inverse, the polar decomposition and the norms. Batched GPU SVD
   is limited to blocks of size <= 32; the Hermitian eigensolver is batched for every size.
4. Any X and f holomorphic near the spectrum: f(X) = (2 pi i)^{-1} oint f(z) (z - X)^{-1} dz, by the
   trapezoidal rule on a circle (batched LU solves), with an error estimate from doubling the nodes.

Gradients. Levels 1 and 2 use the Daleckii-Krein formula dF = V (L o V^* dX V) V^* with the Loewner matrix
of divided differences of f (f' for coinciding eigenvalues), which is finite at degenerate spectra; level 2
assumes f holomorphic near the spectrum (the derivative is then the Frechet derivative of the holomorphic
calculus in every direction). Level 3 has explicit backward formulas (L_ij = 1/(s_i + s_j) for |X|, the
constant-rank formula for the pseudo-inverse). Level 4 is differentiated through the linear solves.
"""
import math
import warnings
from typing import Callable, Optional, Tuple

import torch

from . import cost

import contextlib

_CACHE = {'enabled': False}


def cache_enabled() -> bool:
    return _CACHE['enabled']


@contextlib.contextmanager
def caching(enabled: bool = True):
    """
    Within the block, every operator keeps its eigendecomposition after the first use (see
    Operator.cache_spectrum): repeated functions of the same operator then cost matrix products only.
    """
    previous = _CACHE['enabled']
    _CACHE['enabled'] = bool(enabled)
    try:
        yield
    finally:
        _CACHE['enabled'] = previous


# golden-ratio-like irrational mixing coefficients for H + alpha K (level 2)
_ALPHAS = (0.6180339887498949, 0.3247179572447460, 0.7548776662466927, 0.1415926535897932)


# ----------------------------------------------------------------------
# divided differences
# ----------------------------------------------------------------------
def values_and_derivative(f: Callable[[torch.Tensor], torch.Tensor], w: torch.Tensor,
                          create_graph: bool = False) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    f(w) and f'(w) for an elementwise f, with one evaluation of f. For real w, f' is the real derivative;
    for complex w, f' is the complex derivative of a holomorphic f (PyTorch returns conj(f') for the
    gradient of Re f). The values keep their graph to parameters of f. Non-floating results (indicators)
    are converted to the dtype of w; f' is 0 where f does not depend on w.
    """
    with torch.enable_grad():
        wg = w if create_graph else w.detach().requires_grad_(True)
        y = f(wg)
        if not (y.is_floating_point() or y.is_complex()):
            y = y.to(w.dtype)
        if not y.requires_grad:
            return y, torch.zeros_like(w if torch.is_complex(w) or not torch.is_complex(y) else y)

        def d(part):
            g, = torch.autograd.grad(part.sum(), wg, retain_graph=True, create_graph=create_graph, allow_unused=True)
            return torch.zeros_like(wg) if g is None else g
        if torch.is_complex(w):
            dy = d(y.real if torch.is_complex(y) else y).conj()      # holomorphic f: grad Re f = conj(f')
        elif torch.is_complex(y):
            dy = torch.complex(d(y.real), d(y.imag))
        else:
            dy = d(y)
    return y, dy


def divided_differences(w: torch.Tensor, fw: torch.Tensor, dfw: torch.Tensor) -> torch.Tensor:
    """
    Loewner matrix L_ij = (f(w_i) - f(w_j)) / (w_i - w_j), and (f'(w_i) + f'(w_j)) / 2 where w_i and w_j
    coincide up to eps^(1/3) relative to max(|w_i|, |w_j|) (which balances the truncation error of the
    midpoint value against the cancellation error of the quotient). w may be complex. Shape (..., k, k).
    """
    dw = w.unsqueeze(-1) - w.unsqueeze(-2)
    df = fw.unsqueeze(-1) - fw.unsqueeze(-2)
    finfo = torch.finfo(w.real.dtype if torch.is_complex(w) else w.dtype)
    scale = torch.maximum(w.abs().unsqueeze(-1), w.abs().unsqueeze(-2)).clamp(min=finfo.tiny)
    close = dw.abs() <= finfo.eps ** (1 / 3) * scale
    dtype = torch.promote_types(torch.promote_types(df.dtype, dw.dtype), dfw.dtype)
    quotient = df.to(dtype) / torch.where(close, torch.ones_like(dw), dw).to(dtype)
    return torch.where(close, 0.5 * (dfw.unsqueeze(-1) + dfw.unsqueeze(-2)).to(dtype), quotient)


# ----------------------------------------------------------------------
# levels 1 and 2: spectral decompositions with the Hermitian eigensolver
# ----------------------------------------------------------------------
def hermitian_eig(X: torch.Tensor, what: str = "eigh") -> Tuple[torch.Tensor, torch.Tensor]:
    """Eigenvalues (real, ascending) and eigenvectors of the self-adjoint part of a stack X."""
    Xh = 0.5 * (X + X.conj().transpose(-2, -1))
    return cost.batched_call(torch.linalg.eigh, Xh, what)


def normal_eig(X: torch.Tensor, what: str = "eigh", rtol: Optional[float] = None
               ) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Eigenvalues (complex) and a unitary eigenbasis of a stack X of normal matrices, from the Hermitian
    eigensolver: X = H + i K with commuting self-adjoint H, K, and the eigenvectors of H + alpha K diagonalise
    both when alpha separates the eigenvalues (true for all but finitely many alpha; several irrational
    values are tried). Raises ValueError if X is not normal to working precision.
    """
    Xc = X if torch.is_complex(X) else X.to(torch.promote_types(X.dtype, torch.complex64))
    Xs = Xc.conj().transpose(-2, -1)
    H, K = 0.5 * (Xc + Xs), -0.5j * (Xc - Xs)
    eps = torch.finfo(Xc.real.dtype).eps
    tol = (rtol if rtol is not None else 1e3 * eps * max(1, X.shape[-1]))
    scale = Xc.abs().amax(dim=(-2, -1)).clamp(min=torch.finfo(Xc.real.dtype).tiny)
    for alpha in _ALPHAS:
        _, V = cost.batched_call(torch.linalg.eigh, H + alpha * K, what)
        D = V.conj().transpose(-2, -1) @ Xc @ V
        lam = torch.diagonal(D, dim1=-2, dim2=-1)
        off = (D - torch.diag_embed(lam)).abs().amax(dim=(-2, -1))
        if bool((off <= tol * scale).all()):
            return lam, V
    raise ValueError("the operator is not normal to working precision (use holomorphic_function for "
                     "non-normal operators)")


class _SpectralReconstruction(torch.autograd.Function):
    """
    F = V diag(f(w)) V^* as a function of X (diagonalised by V) and of the values f(w). The derivative in X is
    dF = V (L o V^* dX V) V^* with the Loewner matrix L; the values carry the gradient of parameters of f.
    With create_graph the backward pass recomputes the decomposition and L differentiably (second derivatives
    are then correct at non-degenerate spectra).
    """

    @staticmethod
    def forward(ctx, X, fw, V, L, recompute):
        ctx.save_for_backward(X, V, L)
        ctx.recompute = recompute
        ctx.fw_is_real = not torch.is_complex(fw)
        return (V * fw.to(V.dtype).unsqueeze(-2)) @ V.conj().transpose(-2, -1)

    @staticmethod
    def backward(ctx, G):
        X, V, L = ctx.saved_tensors
        if torch.is_grad_enabled():
            V, L = ctx.recompute(X)
        Vh = V.conj().transpose(-2, -1)
        Gt = Vh @ G.to(V.dtype) @ V
        gX = V @ (L.conj().to(Gt.dtype) * Gt) @ Vh
        if not torch.is_complex(X) and torch.is_complex(gX):
            gX = gX.real
        gfw = torch.diagonal(Gt, dim1=-2, dim2=-1)
        if ctx.fw_is_real and torch.is_complex(gfw):
            gfw = gfw.real
        return (gX if ctx.needs_input_grad[0] else None,
                gfw if ctx.needs_input_grad[1] else None, None, None, None)


def spectral_apply(X: torch.Tensor, f: Callable[[torch.Tensor], torch.Tensor], kind: str = 'hermitian',
                   decomposition: Optional[Tuple[torch.Tensor, torch.Tensor]] = None,
                   what: str = "eigh") -> torch.Tensor:
    """
    f(X) for a stack X (n, k, k) of self-adjoint (kind='hermitian') or normal (kind='normal') matrices and an
    elementwise f on the eigenvalues (n, k). `decomposition` = (w, V) may be passed from a cache. For
    self-adjoint X the function is evaluated at the self-adjoint part of X.
    """
    if kind == 'hermitian':
        Xin = 0.5 * (X + X.conj().transpose(-2, -1))
        decompose = lambda Y: torch.linalg.eigh(0.5 * (Y + Y.conj().transpose(-2, -1)))
    elif kind == 'normal':
        Xin = X
        decompose = lambda Y: normal_eig(Y)
    else:
        raise ValueError(f"unknown kind {kind!r}")
    if decomposition is None:
        with torch.no_grad():
            w, V = hermitian_eig(Xin.detach(), what) if kind == 'hermitian' else normal_eig(Xin.detach(), what)
    else:
        w, V = decomposition
    needs_grad = Xin.requires_grad and torch.is_grad_enabled()
    if needs_grad:
        fw, dfw = values_and_derivative(f, w)
    else:
        fw = f(w)
        if not (fw.is_floating_point() or fw.is_complex()):
            fw = fw.to(w.dtype)
    if not needs_grad:
        return (V * fw.to(torch.promote_types(V.dtype, fw.dtype)).unsqueeze(-2)) @ \
            V.conj().transpose(-2, -1).to(torch.promote_types(V.dtype, fw.dtype))
    L = divided_differences(w, fw.detach(), dfw)

    def recompute(Y):                                   # differentiable decomposition for second derivatives
        w2, V2 = decompose(Y)
        fw2, dfw2 = values_and_derivative(f, w2, create_graph=True)
        return V2, divided_differences(w2, fw2, dfw2)
    Vc = V if torch.is_complex(fw) == torch.is_complex(V) or torch.is_complex(V) else V.to(fw.dtype)
    return _SpectralReconstruction.apply(Xin, fw, Vc, L, recompute)


# ----------------------------------------------------------------------
# level 3: singular values from the Hermitian eigensolver of the dilation
# ----------------------------------------------------------------------
def _dilation(X: torch.Tensor) -> torch.Tensor:
    Z = torch.zeros_like(X)
    return torch.cat([torch.cat([Z, X], dim=-1), torch.cat([X.conj().transpose(-2, -1), Z], dim=-1)], dim=-2)


def use_dilation(X: torch.Tensor, method: str = 'auto') -> bool:
    """'auto': the dilation on GPUs for blocks larger than 32 (where batched SVD falls back to a loop)."""
    if method == 'svd':
        return False
    if method == 'eigh':
        return True
    return X.is_cuda and X.shape[-1] > 32


def singular_values(X: torch.Tensor, method: str = 'auto') -> torch.Tensor:
    """Singular values (n, k), descending, of a stack of square blocks."""
    k = X.shape[-1]
    if use_dilation(X, method):
        lam = cost.batched_call(torch.linalg.eigvalsh, _dilation(X), f"eigvalsh (dilation of size {2 * k})")
        return lam[..., k:].flip(-1).clamp(min=0)
    return cost.batched_call(torch.linalg.svdvals, X, f"svdvals (blocks of size {k})", kind='svd')


def svd(X: torch.Tensor, method: str = 'auto') -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    X = U diag(s) V^*, without gradients (callers attach their own backward). With the dilation, the columns
    of U and V belonging to zero singular values are not determined (they never enter the results below).
    """
    k = X.shape[-1]
    with torch.no_grad():
        if use_dilation(X, method):
            lam, Q = cost.batched_call(torch.linalg.eigh, _dilation(X.detach()), f"eigh (dilation of size {2 * k})")
            Qp = Q[..., k:].flip(-1) * math.sqrt(2.0)
            return Qp[..., :k, :], lam[..., k:].flip(-1).clamp(min=0), Qp[..., k:, :]
        U, s, Vh = cost.batched_call(torch.linalg.svd, X.detach(), f"SVD (blocks of size {k})", kind='svd')
        return U, s, Vh.conj().transpose(-2, -1)


class _Abs(torch.autograd.Function):
    """|X| = (X^* X)^{1/2} = V diag(s) V^*; backward through sqrt at X^* X with L_ij = 1/(s_i + s_j)."""

    @staticmethod
    def forward(ctx, X, method):
        U, s, V = svd(X, method)
        ctx.save_for_backward(X, s, V)
        return (V * s.to(V.dtype).unsqueeze(-2)) @ V.conj().transpose(-2, -1)

    @staticmethod
    @torch.autograd.function.once_differentiable
    def backward(ctx, G):
        X, s, V = ctx.saved_tensors
        S = s.unsqueeze(-1) + s.unsqueeze(-2)
        L = torch.where(S > torch.finfo(s.dtype).tiny, 1.0 / torch.where(S > 0, S, torch.ones_like(S)),
                        torch.zeros_like(S))          # kernel of X: |X| is not differentiable there; 0 is used
        Vh = V.conj().transpose(-2, -1)
        gY = V @ (L.to(V.dtype) * (Vh @ G.to(V.dtype) @ V)) @ Vh       # gradient at Y = X^* X
        return X @ (gY + gY.conj().transpose(-2, -1)), None


def abs_(X: torch.Tensor, method: str = 'auto') -> torch.Tensor:
    return _Abs.apply(X, method)


def _pinv_threshold(s: torch.Tensor, k: int, tol: Optional[float]) -> torch.Tensor:
    if tol is not None:
        return torch.as_tensor(tol, dtype=s.dtype, device=s.device)
    return k * torch.finfo(s.dtype).eps * s.amax(dim=-1, keepdim=True)


class _Pinv(torch.autograd.Function):
    """Moore-Penrose pseudo-inverse with a relative cutoff; backward by the constant-rank formula."""

    @staticmethod
    def forward(ctx, X, tol, method):
        U, s, V = svd(X, method)
        thr = _pinv_threshold(s, X.shape[-1], tol)
        big = s > thr
        s_inv = torch.where(big, 1.0 / torch.where(big, s, torch.ones_like(s)), torch.zeros_like(s))
        P = (V * s_inv.to(V.dtype).unsqueeze(-2)) @ U.conj().transpose(-2, -1)
        ctx.save_for_backward(X, P)
        return P

    @staticmethod
    @torch.autograd.function.once_differentiable
    def backward(ctx, G):
        A, P = ctx.saved_tensors
        Ph, Gh = P.conj().transpose(-2, -1), G.conj().transpose(-2, -1)
        I = torch.eye(A.shape[-1], dtype=A.dtype, device=A.device)
        grad = -Ph @ G @ Ph + (I - A @ P) @ Gh @ P @ Ph + Ph @ P @ Gh @ (I - P @ A)
        return grad, None, None


def pinv(X: torch.Tensor, tol: Optional[float] = None, method: str = 'auto') -> torch.Tensor:
    return _Pinv.apply(X, tol, method)


def polar(X: torch.Tensor, method: str = 'auto') -> Tuple[torch.Tensor, torch.Tensor]:
    """
    X = W |X| with W = U V^* (unitary for invertible X; on the kernel of a singular X the factor W is the
    partial isometry from the SVD) and |X| from _Abs (differentiable). W has no gradient.
    """
    U, s, V = svd(X, method)
    thr = _pinv_threshold(s, X.shape[-1], None)
    keep = (s > thr).to(U.dtype).unsqueeze(-2)
    return (U * keep) @ V.conj().transpose(-2, -1), abs_(X, method)


# ----------------------------------------------------------------------
# level 4: holomorphic functional calculus by contour integration
# ----------------------------------------------------------------------
def holomorphic(X: torch.Tensor, f: Callable[[torch.Tensor], torch.Tensor], center=None, radius=None,
                nodes: int = 32, max_nodes: int = 2048, rtol: Optional[float] = None) -> Tuple[torch.Tensor, float]:
    """
    f(X) = (2 pi i)^{-1} oint f(z) (z - X)^{-1} dz on the circle |z - center| = radius, by the trapezoidal
    rule with N nodes, N doubled until two successive results agree to rtol (the nodes of N are reused for
    2N). Default contour: center Tr(X)/k, radius 2 ||X - center||_F, which encloses the spectrum; f must be
    holomorphic on the closed disk and accept complex tensors elementwise. Returns (f(X), error estimate).
    """
    k = X.shape[-1]
    cdt = X.dtype if torch.is_complex(X) else torch.promote_types(X.dtype, torch.complex64)
    Xc = X.to(cdt)
    eye = torch.eye(k, dtype=cdt, device=X.device)
    if center is None:
        center = torch.diagonal(Xc, dim1=-2, dim2=-1).sum(-1) / k
    center = torch.as_tensor(center, dtype=cdt, device=X.device).expand(X.shape[0])
    if radius is None:
        spread = torch.linalg.matrix_norm(Xc - center[:, None, None] * eye).detach()
        radius = 2.0 * spread + torch.finfo(Xc.real.dtype).eps ** 0.5
    radius = torch.as_tensor(radius, dtype=Xc.real.dtype, device=X.device).expand(X.shape[0])
    if rtol is None:
        rtol = 1e-10 if Xc.real.dtype == torch.float64 else 1e-5

    def partial_sum(phases):                              # sum over nodes z = c + r e^{i phi} of f(z)(z-c)(z-X)^{-1}
        e = torch.exp(1j * phases.to(Xc.real.dtype)).to(cdt)                      # (m,)
        z = center[None, :] + radius[None, :] * e[:, None]                         # (m, n)
        fz = f(z)
        if not torch.isfinite(torch.view_as_real(fz.to(cdt))).all():
            raise ValueError("f is not finite on the contour: choose center and radius so that the circle encloses "
                             "the spectrum and f is holomorphic on the closed disk")
        R = cost.batched_call(torch.linalg.inv, (z[..., None, None] * eye - Xc[None]).reshape(-1, k, k),
                              f"resolvents (blocks of size {k})", kind='lu').reshape(len(e), *Xc.shape)
        w = (fz.to(cdt) * (z - center[None, :]))[..., None, None]
        return (w * R).sum(0)

    N = nodes
    total = partial_sum(2 * math.pi * torch.arange(N, device=X.device, dtype=torch.float64) / N)
    F = total / N
    err = float('inf')
    while N < max_nodes:
        total = total + partial_sum(2 * math.pi * (torch.arange(N, device=X.device, dtype=torch.float64) + 0.5) / N)
        N *= 2
        F_new = total / N
        scale = F_new.detach().abs().amax().clamp(min=1.0)
        err = float((F_new - F).detach().abs().amax() / scale)
        F = F_new
        if err <= rtol:
            break
    if err > rtol:
        warnings.warn(f"holomorphic_function: contour integral not converged with {N} nodes (estimated relative "
                      f"error {err:.1e}); adjust the contour (center, radius) or increase max_nodes", RuntimeWarning,
                      stacklevel=3)
    if not torch.is_complex(X):
        if F.detach().imag.abs().amax() > 1e3 * max(rtol, err if err < float('inf') else 0) * \
                F.detach().abs().amax().clamp(min=1.0):
            raise ValueError("f(X) is not real for a real X: f must satisfy f(conj z) = conj f(z) "
                             "(use a complex-valued algebra otherwise)")
        F = F.real
    return F, err
