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
calculus in every direction). Second derivatives (create_graph) recompute the decomposition differentiably
(eigh for self-adjoint, eig for normal X) and are correct at non-degenerate spectra. Level 3 has explicit
backward formulas (L_ij = 1/(s_i + s_j) for |X| on the range, the constant-rank formula for the pseudo-inverse);
level 4 is differentiated through the linear solves.
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
    values = y
    if not create_graph:
        # the values must not carry the graph of the internal leaf wg (they would report requires_grad and keep
        # the graph alive): evaluate f once more on w itself, in the caller's grad mode (elementwise, O(k))
        values = f(w)
        if not (values.is_floating_point() or values.is_complex()):
            values = values.to(w.dtype)
    with torch.enable_grad():

        def d(part):
            g, = torch.autograd.grad(part.sum(), wg, retain_graph=True, create_graph=create_graph, allow_unused=True)
            return torch.zeros_like(wg) if g is None else g
        if torch.is_complex(w):
            dy = d(y.real if torch.is_complex(y) else y).conj()      # holomorphic f: grad Re f = conj(f')
        elif torch.is_complex(y):
            dy = torch.complex(d(y.real), d(y.imag))
        else:
            dy = d(y)
    return values, dy


def divided_differences(w: torch.Tensor, fw: torch.Tensor, dfw: torch.Tensor) -> torch.Tensor:
    """
    Loewner matrix L_ij = (f(w_i) - f(w_j)) / (w_i - w_j), and (f'(w_i) + f'(w_j)) / 2 where w_i and w_j
    coincide up to eps^(1/3) relative to max(|w_i|, |w_j|) (which balances the truncation error of the
    midpoint value against the cancellation error of the quotient). w may be complex. Shape (..., k, k).
    The midpoint value is also used for eigenvalues that are close on the scale of the block (|w_i - w_j| <=
    eps^(1/3) max_l |w_l|) when the quotient is dominated by rounding, eps (|f(w_i)| + |f(w_j)|) >
    eps^(2/3) |w_i - w_j| max(|f'(w_i)|, |f'(w_j)|): e.g. a degenerate eigenvalue 0 computed as +-1e-17 (where
    (exp(1e-17) - exp(-1e-17)) / 2e-17 = 0 instead of 1).
    """
    dw = w.unsqueeze(-1) - w.unsqueeze(-2)
    df = fw.unsqueeze(-1) - fw.unsqueeze(-2)
    finfo = torch.finfo(w.real.dtype if torch.is_complex(w) else w.dtype)
    scale = torch.maximum(w.abs().unsqueeze(-1), w.abs().unsqueeze(-2)).clamp(min=finfo.tiny)
    close = dw.abs() <= finfo.eps ** (1 / 3) * scale
    if w.shape[-1] > 1:
        radius = w.abs().amax(dim=-1, keepdim=True).unsqueeze(-1)              # scale of the block's spectrum
        fa, da = fw.detach().abs(), dfw.detach().abs()
        rounding = finfo.eps * (fa.unsqueeze(-1) + fa.unsqueeze(-2))
        slope = torch.maximum(da.unsqueeze(-1), da.unsqueeze(-2))
        close = close | ((dw.abs() <= finfo.eps ** (1 / 3) * radius)
                         & (rounding > finfo.eps ** (2 / 3) * dw.abs() * slope))
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


def _normal_pass(Xc: torch.Tensor, H: torch.Tensor, K: torch.Tensor, a1: float, a2: float, what: str):
    """Eigenbasis from eigh(H + a1 K), refined inside clusters of nearly equal mu = Re + a1 Im of the eigenvalues
    (where eigh may mix eigenvectors of distinct eigenvalues) by eigh of H_D + a2 K_D restricted to the clusters."""
    mu, V = cost.batched_call(torch.linalg.eigh, H + a1 * K, what)
    D = V.conj().transpose(-2, -1) @ Xc @ V
    eps = torch.finfo(Xc.real.dtype).eps
    scale = torch.diagonal(D, dim1=-2, dim2=-1).abs().amax(dim=-1, keepdim=True).clamp(min=torch.finfo(Xc.real.dtype).tiny)
    brk = torch.cat([torch.ones_like(mu[..., :1], dtype=torch.bool), (mu[..., 1:] - mu[..., :-1]) > eps ** 0.25 * scale],
                    dim=-1)
    if not bool(brk.all()):                 # some cluster has two or more values: refine (one synchronisation)
        cid = torch.cumsum(brk.long(), dim=-1)
        same = cid.unsqueeze(-1) == cid.unsqueeze(-2)             # pairs in the same cluster (sorted mu: contiguous)
        Ds = D.conj().transpose(-2, -1)
        M = 0.5 * (D + Ds) - 0.5j * a2 * (D - Ds)
        _, W = cost.batched_call(torch.linalg.eigh, torch.where(same, M, torch.zeros_like(M)), what)
        V = V @ W
        D = V.conj().transpose(-2, -1) @ Xc @ V
    lam = torch.diagonal(D, dim1=-2, dim2=-1)
    offdiag = D.clone()
    offdiag.diagonal(dim1=-2, dim2=-1).zero_()
    off = offdiag.abs().amax(dim=(-2, -1)) / lam.abs().amax(dim=-1).clamp(min=torch.finfo(Xc.real.dtype).tiny)
    return lam, V, off


def normal_eig(X: torch.Tensor, what: str = "eigh", rtol: Optional[float] = None
               ) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Eigenvalues (complex) and a unitary eigenbasis of a stack X of normal matrices, from the Hermitian
    eigensolver: X = H + i K with commuting self-adjoint H, K, and the eigenvectors of H + alpha K diagonalise
    both when alpha separates the eigenvalues (true for all but finitely many alpha). Eigenvectors of nearly
    tied values of Re + alpha Im are refined with a second coefficient; several irrational coefficients are tried
    and the best result is kept per matrix. Raises ValueError if some X is not normal to working precision
    (relative off-diagonal residual above rtol, default 100 k eps).
    """
    out_dtype = X.dtype if torch.is_complex(X) else torch.promote_types(X.dtype, torch.complex64)
    if out_dtype == torch.complex64 and X.shape[-1] > 64:
        # in single precision the separation of nearly tied Re + alpha Im values fails for blocks of a few hundred;
        # the decomposition is computed in double precision and returned in single precision (smaller blocks are
        # tried in single precision first, see below)
        tol32 = rtol if rtol is not None else 100 * torch.finfo(torch.float32).eps * max(1, X.shape[-1])
        lam, V = normal_eig(X.to(torch.complex128), what, tol32)
        return lam.to(out_dtype), V.to(out_dtype)
    Xc = X.to(out_dtype)
    Xs = Xc.conj().transpose(-2, -1)
    H, K = 0.5 * (Xc + Xs), -0.5j * (Xc - Xs)
    k = X.shape[-1]
    tol = rtol if rtol is not None else 100 * torch.finfo(Xc.real.dtype).eps * max(1, k)
    best = None
    for i, a1 in enumerate(_ALPHAS):
        lam, V, off = _normal_pass(Xc, H, K, a1, _ALPHAS[(i + 1) % len(_ALPHAS)], what)
        if best is None:
            best = [lam, V, off]
        else:
            better = off < best[2]
            best[0] = torch.where(better[:, None], lam, best[0])
            best[1] = torch.where(better[:, None, None], V, best[1])
            best[2] = torch.minimum(off, best[2])
        if bool((best[2] <= tol).all()):
            return best[0], best[1]
    if Xc.dtype == torch.complex64:                    # single precision did not separate: redo in double
        lam, V = normal_eig(X.to(torch.complex128), what, tol)
        return lam.to(torch.complex64), V.to(torch.complex64)
    raise ValueError(f"the operator is not normal to working precision (relative residual "
                     f"{best[2].max().item():.1e} > {tol:.1e}); use holomorphic_function for non-normal operators")


def _warn_if_degenerate(w: torch.Tensor):
    """Second derivatives come from differentiating the eigendecomposition, which is singular at repeated
    eigenvalues: warn if the spectrum is (nearly) degenerate (one synchronisation, only with create_graph)."""
    k = w.shape[-1]
    if k < 2:
        return
    wd = w.detach()
    gap = (wd.unsqueeze(-1) - wd.unsqueeze(-2)).abs()
    gap = gap + torch.eye(k, device=w.device, dtype=gap.dtype) * gap.amax().clamp(min=1)
    scale = wd.abs().amax().clamp(min=torch.finfo(gap.dtype).tiny)
    if bool(gap.amin() <= torch.finfo(gap.dtype).eps ** (1 / 3) * scale):
        warnings.warn("second derivatives of a spectral function at a (nearly) degenerate spectrum are not reliable "
                      "(they differentiate the eigendecomposition); use holomorphic_function for exact second "
                      "derivatives", RuntimeWarning, stacklevel=4)


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
        Vinv = V.conj().transpose(-2, -1)                      # unitary eigenbasis
        if torch.is_grad_enabled():                            # create_graph: differentiable recomputation
            V, Vinv, L = ctx.recompute(X)
        # F = V f(Λ) V^{-1}, dF = V (L o V^{-1} dX V) V^{-1}  =>  grad_X = V^{-*} (conj(L) o V^* G V^{-*}) V^*
        Vh, Vinvh = V.conj().transpose(-2, -1), Vinv.conj().transpose(-2, -1)
        Gt = Vh @ G.to(V.dtype) @ Vinvh
        gX = Vinvh @ (L.conj().to(Gt.dtype) * Gt) @ Vh
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
    correct = kind == 'normal'          # first-order correction for the residual off-diagonal part of V^* X V
    if needs_grad or correct:
        fw, dfw = values_and_derivative(f, w)
    else:
        fw = f(w)
        if not (fw.is_floating_point() or fw.is_complex()):
            fw = fw.to(w.dtype)
    correction = None
    if correct:
        with torch.no_grad():
            Vh = V.conj().transpose(-2, -1)
            E = Vh @ Xin.detach().to(V.dtype) @ V
            E = E - torch.diag_embed(torch.diagonal(E, dim1=-2, dim2=-1))
            L0 = divided_differences(w, fw.detach(), dfw)
            correction = V @ (L0.to(E.dtype) * E) @ Vh      # f(V(Λ+E)V^*) = V(f(Λ) + L o E)V^* + O(|E|^2)
    if not needs_grad:
        dt = torch.promote_types(V.dtype, fw.dtype)
        out = (V.to(dt) * fw.to(dt).unsqueeze(-2)) @ V.conj().transpose(-2, -1).to(dt)
        return out if correction is None else out + correction.to(dt)
    L = divided_differences(w, fw.detach(), dfw)

    def recompute(Y):                                   # differentiable decomposition for second derivatives
        if kind == 'hermitian':
            w2, V2 = decompose(Y)
            V2inv = V2.conj().transpose(-2, -1)
        else:                                           # off the normal matrices: the eigenvectors of Y itself
            w2, V2 = torch.linalg.eig(Y)
            V2inv = torch.linalg.inv(V2)
        _warn_if_degenerate(w2)
        fw2, dfw2 = values_and_derivative(f, w2, create_graph=True)
        return V2, V2inv, divided_differences(w2, fw2, dfw2)
    Vc = V if torch.is_complex(fw) == torch.is_complex(V) or torch.is_complex(V) else V.to(fw.dtype)
    out = _SpectralReconstruction.apply(Xin, fw, Vc, L, recompute)
    return out if correction is None else out + correction.to(out.dtype)


# ----------------------------------------------------------------------
# level 3: singular values from the Hermitian eigensolver of the dilation
# ----------------------------------------------------------------------
def _dilation(X: torch.Tensor) -> torch.Tensor:
    k = X.shape[-1]
    D = X.new_zeros(*X.shape[:-2], 2 * k, 2 * k)
    D[..., :k, k:] = X
    D[..., k:, :k] = X.conj().transpose(-2, -1)
    return D


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
    """
    |X| = (X^* X)^{1/2} = V diag(s) V^*. Backward through the square root at Y = X^* X, using only the right
    singular vectors V_r of the singular values above the relative cutoff of the pseudo-inverse:
    grad_Y = V_r (L o V_r^* G V_r) V_r^* + V_r S_r^{-1} V_r^* G P + P G V_r S_r^{-1} V_r^*, L_ij = 1/(s_i + s_j),
    P = 1 - V_r V_r^* (the derivative of |X| along directions that keep the kernel), and grad_X = X (grad_Y + h.c.).
    """

    @staticmethod
    def forward(ctx, X, method, usv=None):
        U, s, V = svd(X, method) if usv is None else usv
        keep = s > _pinv_threshold(s, X.shape[-1], None)
        ctx.save_for_backward(X, s, V, keep)
        return (V * (s * keep).to(V.dtype).unsqueeze(-2)) @ V.conj().transpose(-2, -1)

    @staticmethod
    @torch.autograd.function.once_differentiable
    def backward(ctx, G):
        X, s, V, keep = ctx.saved_tensors
        Vr = V * keep.to(V.dtype).unsqueeze(-2)                     # kernel columns zeroed (never used)
        s_safe = torch.where(keep, s, torch.ones_like(s))
        S = s_safe.unsqueeze(-1) + s_safe.unsqueeze(-2)
        both = keep.unsqueeze(-1) & keep.unsqueeze(-2)
        L = torch.where(both, 1.0 / S, torch.zeros_like(S))
        Vrh = Vr.conj().transpose(-2, -1)
        G = G.to(V.dtype)
        s_inv = torch.where(keep, 1.0 / s_safe, torch.zeros_like(s)).to(V.dtype)
        # with P = 1 - V_r V_r^*, A = V_r^* G, M = A V_r, Bm = G V_r (no explicit projector):
        # V_r S^-1 V_r^* G P = (V_r S^-1)(A - M V_r^*),  P G V_r S^-1 V_r^* = ((Bm - V_r M) S^-1) V_r^*
        A = Vrh @ G
        M = A @ Vr
        Bm = G @ Vr
        gY = (Vr @ (L.to(V.dtype) * M) @ Vrh + (Vr * s_inv.unsqueeze(-2)) @ (A - M @ Vrh)
              + ((Bm - Vr @ M) * s_inv.unsqueeze(-2)) @ Vrh)
        return X @ (gY + gY.conj().transpose(-2, -1)), None, None


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
    partial isometry from the SVD) and |X| from _Abs (differentiable). W has no gradient. Two Newton-Schulz steps
    W <- (3W - W W^* W)/2 restore orthonormality lost to rounding (the dilation mixes the eigenvectors of +-s
    for singular values small compared with the largest).
    """
    U, s, V = svd(X, method)
    thr = _pinv_threshold(s, X.shape[-1], None)
    keep = (s > thr).to(U.dtype).unsqueeze(-2)
    W = (U * keep) @ V.conj().transpose(-2, -1)
    if use_dilation(X, method):                     # only the dilation loses orthonormality to rounding
        for _ in range(2):
            W = 1.5 * W - 0.5 * W @ W.conj().transpose(-2, -1) @ W
    return W, _Abs.apply(X, method, (U, s, V))      # |X| from the same decomposition


# ----------------------------------------------------------------------
# level 4: holomorphic functional calculus by contour integration
# ----------------------------------------------------------------------
def default_contour(X: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Circle around the spectrum: center c = Tr(X)/k and radius 1.25 max|lambda_i - c| + 0.5 ||X - c||_F / sqrt(k)
    (the second term keeps the circle away from the eigenvalues of non-normal blocks, whose resolvent is large near
    the spectrum, and is positive for X = c + nilpotent). Eigenvalues from torch.linalg.eigvals, without gradients.
    """
    return _contour_and_spectrum(X)[:2]


def _contour_and_spectrum(X: torch.Tensor):
    """(center, radius, eigenvalues) of default_contour, with one eigenvalue computation."""
    k = X.shape[-1]
    with torch.no_grad():
        cdt = X.dtype if torch.is_complex(X) else torch.promote_types(X.dtype, torch.complex64)
        Xc = X.detach().to(cdt)
        lam = torch.linalg.eigvals(Xc)
        c = torch.diagonal(Xc, dim1=-2, dim2=-1).sum(-1) / k
        Y = Xc.clone()
        Y.diagonal(dim1=-2, dim2=-1).sub_(c[..., None])
        rho = (lam - c[..., None]).abs().amax(dim=-1)
        spread = torch.linalg.matrix_norm(Y) / math.sqrt(k)
        r = 1.25 * rho + 0.5 * spread
        # a block that is scalar up to rounding would get a radius at the rounding level: keep the circle at a
        # distance eps^(1/3) relative to the size of X (any circle around c encloses a scalar block's spectrum)
        size = torch.maximum(c.abs(), torch.linalg.matrix_norm(Xc) / math.sqrt(k))
        floor = torch.finfo(size.dtype).eps ** (1 / 3) * size
        r = torch.maximum(r, floor)
        r = torch.where(r > 0, r, torch.ones_like(r))          # X = 0: any circle around 0
    return c, r, lam


def holomorphic(X: torch.Tensor, f: Callable[[torch.Tensor], torch.Tensor], center=None, radius=None,
                nodes: int = 32, max_nodes: int = 2048, rtol: Optional[float] = None,
                memory_bytes: int = 2 ** 28, atol: float = 0.0) -> Tuple[torch.Tensor, float]:
    """
    f(X) = (2 pi i)^{-1} oint f(z) (z - X)^{-1} dz on the circle |z - center| = radius, by the trapezoidal
    rule with N nodes, N doubled until two successive results agree to rtol relative to their size (the nodes of N
    are reused for 2N). Default contour: default_contour(X), a circle just around the spectrum; f must be
    holomorphic on the closed disk and accept complex tensors elementwise. center and radius may be numbers,
    tensors (n,) or callables X -> (n,). Raises ValueError if the estimated error stays above sqrt(rtol) (a
    warning between rtol and sqrt(rtol)). The error is relative to max(|f(X)|, atol): when f(X) may vanish (sin at
    X = 0, z^3 at a nilpotent block) pass atol > 0, a scale for f near the spectrum (e.g. 1.0); the absolute
    accuracy accepted where |f(X)| < atol is then rtol * atol. Returns (f(X), error estimate).
    """
    k, n = X.shape[-1], X.shape[0]
    cdt = X.dtype if torch.is_complex(X) else torch.promote_types(X.dtype, torch.complex64)
    Xc = X.to(cdt)
    c0, r0, lam = _contour_and_spectrum(X)          # one (non-Hermitian) eigenvalue computation
    center = c0 if center is None else center
    radius = r0 if radius is None else radius
    center = center(X) if callable(center) else center
    radius = radius(X) if callable(radius) else radius
    center = torch.as_tensor(center, dtype=cdt, device=X.device).expand(n)
    radius = torch.as_tensor(radius, dtype=Xc.real.dtype, device=X.device).expand(n)
    if rtol is None:
        rtol = 1e-10 if Xc.real.dtype == torch.float64 else 1e-5
    max_nodes = max(max_nodes, 2 * nodes)
    with torch.no_grad():               # |f| on the spectrum: a lower bound for ||f(X)||, used to detect cancellation
        rho_f = f(lam).to(cdt).abs().amax(dim=-1)
    # live per node and matrix: the shifted matrix, LU workspace and identity of inv, the resolvent, |R|
    chunk = max(1, int(memory_bytes // max(1, 5 * n * k * k * Xc.element_size())))
    bad = torch.zeros((), dtype=torch.bool, device=X.device)       # f not finite / singular resolvent (no sync)

    def partial_sum(phases):                    # sum over nodes z = c + r e^{i phi} of f(z)(z-c)(z-X)^{-1}
        nonlocal bad
        acc = torch.zeros_like(Xc)
        mag = torch.zeros(n, dtype=Xc.real.dtype, device=X.device)    # sum of |terms|: rounding error estimate
        for i in range(0, len(phases), chunk):
            e = torch.exp(1j * phases[i:i + chunk]).to(cdt)                         # (m,)
            z = center[None, :] + radius[None, :] * e[:, None]                       # (m, n)
            fz = f(z).to(cdt)
            bad = bad | ~torch.isfinite(torch.view_as_real(fz)).all()
            A = (-Xc).expand(len(e), *Xc.shape).clone()                            # z - X without z * eye
            A.diagonal(dim1=-2, dim2=-1).add_(z[..., None])
            R, info = cost.batched_call(lambda a: torch.linalg.inv_ex(a), A.reshape(-1, k, k),
                                        f"resolvents (blocks of size {k})", kind='lu')
            bad = bad | (info != 0).any()
            R = R.reshape(len(e), *Xc.shape)
            w = fz * (z - center[None, :])                                           # (m, n)
            acc = acc + torch.einsum('mn,mnij->nij', w, R)
            mag = mag + (w.detach().abs() * R.detach().abs().amax(dim=(-2, -1))).sum(0)
        return acc, mag

    N = nodes
    total, mag = partial_sum(2 * math.pi * torch.arange(N, device=X.device, dtype=torch.float64) / N)
    F = total / N
    eps = torch.finfo(Xc.real.dtype).eps
    err = float('inf')
    while N < max_nodes:
        t2, m2 = partial_sum(2 * math.pi * (torch.arange(N, device=X.device, dtype=torch.float64) + 0.5) / N)
        total, mag = total + t2, mag + m2
        N *= 2
        F_new = total / N
        # errors are relative to |f(X)|, or to atol where |f(X)| is smaller (f(X) = 0 at rounding level: sin at
        # X = 0, z^3 at a nilpotent block; there a relative error is meaningless and the rounding noise of the
        # sum cannot be told apart from catastrophic cancellation without an absolute scale from the caller)
        tiny = torch.finfo(Xc.real.dtype).tiny
        size = F_new.detach().abs().amax(dim=(-2, -1)).clamp(min=max(atol, tiny))
        quadrature = (F_new - F).detach().abs().amax(dim=(-2, -1)) / size
        true_size = torch.where(rho_f > 0, torch.minimum(rho_f, size), size).clamp(min=max(atol, tiny))
        rounding = k * eps * (mag / N) / true_size             # cancellation between the terms of the sum
        err_t = torch.maximum(quadrature, rounding).max()
        err_t = torch.where(bad, torch.full_like(err_t, -1.0), err_t)               # -1 flags a failure
        err = float(err_t)                                     # one synchronisation per doubling
        if err < 0:
            raise ValueError("f is not finite on the contour, or the contour passes through the spectrum: choose "
                             "center and radius so that the circle encloses the spectrum and f is holomorphic on "
                             "the closed disk")
        F = F_new
        if err <= rtol:
            break
    if err > rtol ** 0.5:
        raise ValueError(f"holomorphic_function: the contour integral did not converge (estimated relative error "
                         f"{err:.1e} with {N} nodes); choose a circle closer to the spectrum on which f is holomorphic, "
                         f"or use apply_function (self-adjoint or normal operators) or expm; if f(X) may be 0, pass atol")
    if err > rtol:
        warnings.warn(f"holomorphic_function: estimated relative error {err:.1e} with {N} nodes is above rtol "
                      f"{rtol:.0e}", RuntimeWarning, stacklevel=3)
    if not torch.is_complex(X):
        if F.detach().imag.abs().amax() > 1e3 * max(rtol, err) * F.detach().abs().amax().clamp(min=1.0):
            raise ValueError("f(X) is not real for a real X: f must satisfy f(conj z) = conj f(z) "
                             "(use a complex-valued algebra otherwise)")
        F = F.real
    return F, err
