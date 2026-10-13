"""
Quantum Renyi divergences, Petz recovery maps and noncommutative L^p norms on M = (+)_c M_{k_c}.

Conventions
-----------
States are normal positive functionals on M, represented (as everywhere in the library, see states.py) by their
densities with respect to the blunt trace Tr = sum_c Tr_c: rho = (+)_c rho_c >= 0, omega(A) = Tr(rho A). Plain
positive Operators are accepted as well (they need not be normalised); all functions are batched over the leading
dimension (batch 1 broadcasts against batch B) and return tensors of shape (batch,) in the real dtype of the
algebra. Logarithms are natural (divergences in nats). Functions of positive operators are computed by the spectral
calculus of the library, one batched eigendecomposition per distinct block size, with gradients by the
Daleckii-Krein formula (functional.spectral_apply); there are no host synchronisations in the computational paths.

Definitions (rho, sigma >= 0, rho != 0, alpha in (0, 1) u (1, inf); for states Tr rho = 1 and the normalisation
below is void)

    Petz [Petz 1986]:          D_alpha(rho||sigma)  = (alpha-1)^{-1} log( Tr rho^alpha sigma^{1-alpha} / Tr rho )
    sandwiched [MDSFT 2013, WWY 2014]:
                               ~D_alpha(rho||sigma) = (alpha-1)^{-1} log( Tr (sigma^g rho sigma^g)^alpha / Tr rho ),
                               g = (1-alpha)/(2 alpha)
    Umegaki (alpha -> 1 of both families):
                               D(rho||sigma)        = Tr rho (log rho - log sigma) / Tr rho
    max-relative entropy (alpha -> inf of the sandwiched family):
                               D_max(rho||sigma)    = log min{lambda : rho <= lambda sigma}
                                                    = log lambda_max(sigma^{-1/2} rho sigma^{-1/2})
    alpha -> 0 of the Petz family (Datta's "min-relative entropy"; not the sandwiched alpha = 1/2 quantity that
    Tomamichel's book calls D_min):
                               D_0(rho||sigma)      = -log( Tr(Pi_rho sigma) / Tr rho ),  Pi_rho = support of rho
    fidelity: ~D_{1/2}(rho||sigma) = -log F(rho, sigma) with F = ||sqrt(rho) sqrt(sigma)||_1^2 (states.fidelity).

The quantities Q_alpha = Tr rho^alpha sigma^{1-alpha} and ~Q_alpha = Tr (sigma^g rho sigma^g)^alpha are returned per
sector by sector_contributions.

Support conventions
-------------------
Powers (including negative ones) and logarithms of a positive operator are taken on its numerical support: the
eigenvalues above rtol * lambda_max (lambda_max over all sectors of the sample; default rtol = N eps with
N = sum_c k_c, the cut-off of numpy.linalg.pinv) are kept, the others count as exact zeros, and f(0) := 0 for every
power p != 0 as well as for the logarithm (x^0 is the support projection). The spectrum of the sandwich
sigma^g rho sigma^g is computed as the squared singular values of sigma^g rho^{1/2} (accurate down to
(eps ||sigma^g|| ||rho||^{1/2})^2 rather than eps ||sigma^g rho sigma^g||, which matters for small alpha); singular values
below rtol * max(s_max, ||sigma^g|| ||rho||^{1/2}) count as zeros. Accuracy limit: for alpha < 1 the result depends on
all eigenvalues of the sandwich, and eigenvalues below that resolution (e.g. alpha -> 0 with ill-conditioned states,
whose sandwich has condition number up to cond(sigma)^{2g} cond(rho)) cannot be represented in floating point. Hence Q_alpha and ~Q_alpha only involve
the compression of rho to supp sigma. The support condition supp rho <= supp sigma is tested through the weight of
rho outside supp sigma, leak = Tr rho (1 - Pi_sigma); it counts as violated when leak > supp_tol * Tr rho (default
supp_tol = 100 N eps). Then

    alpha > 1 (both families), D (Umegaki) and D_max:  the divergence is +inf;
    alpha < 1:  the divergence is finite unless Q_alpha = 0 (rho and sigma with orthogonal supports): +inf.

These are the standard conventions (Tomamichel 2016, Ch. 4): for alpha < 1 the divergences are continuous in rho on
the boundary, for alpha >= 1 they jump to +inf. The thresholds are relative, so all results are invariant under
rescaling of the inputs.

Trace conventions
-----------------
With a trace tau_w = sum_c w_c Tr_c (Tr_blunt: w = 1, Tr_norm: w_c = 1/k_c, tau_vN: w_c = 1/(C k_c)), the density of
the state omega with respect to tau_w is d = (+)_c rho_c / w_c. Every divergence above, computed with tau_w from the
densities d (powers of d and traces tau_w), equals the one computed with Tr_blunt from rho: the central factors
w_c^{-alpha} w_c^{alpha-1} w_c (Petz), w_c^{-alpha(2g+1)} w_c = 1 (sandwiched) etc. cancel. Divergences are therefore
functions of the states alone. The argument trace='norm' / 'tau_vN' declares plain Operator inputs to be densities
with respect to that trace (d = DensityMatrix.density(trace)). The L^p norms ||A||_{p,w} do depend on the trace; the
state-weighted (Kosaki) norms ||x||_{p,rho} again do not.

Petz recovery map
-----------------
For a completely positive map Phi: M_in -> M_out (Channel or InterSectorChannel) and sigma >= 0 on M_in,

    R_{sigma,Phi}(X) = sigma^{1/2} Phi^*( Phi(sigma)^{-1/2} X Phi(sigma)^{-1/2} ) sigma^{1/2},

with the adjoint Phi^* for Tr_blunt (the trace in which states are densities) and the inverse square root on
supp Phi(sigma). R is completely positive, R(Phi(sigma)) = sigma (on supp sigma) and Tr R(X) = Tr(Pi_{Phi(sigma)} X)
when Phi is trace preserving. The rotated maps R^t = U_{sigma,t} o R o U_{Phi(sigma),-t}, U_{s,t}(X) = s^{it} X s^{-it}
[JRSWW 2018], and the universal recovery map R_univ = int beta_0(t) R^{t/2} dt, beta_0(t) = pi/2 (cosh(pi t) + 1)^{-1},
satisfy D(rho||sigma) - D(Phi rho||Phi sigma) >= -log F(rho, R_univ(Phi(rho))) [JRSWW 2018, Thm. 2.1].

L^p norms
---------
||A||_{p,w} = (sum_c w_c Tr|A_c|^p)^{1/p} for p in [1, inf), ||A||_inf = max_c ||A_c|| (sectors with w_c > 0), and
the state-weighted norms of Kosaki [Kosaki 1984] ||x||_{p,rho,eta} = ||rho^{eta/p} x rho^{(1-eta)/p}||_p (eta = 1/2:
the symmetric embedding used by Beigi 2013, eta = 1 / 0: left / right). For example
~D_alpha(rho||sigma) = alpha/(alpha-1) log ||sigma^{-1/2} rho sigma^{-1/2}||_{alpha,sigma} [Beigi 2013].

Data processing (theorems; tested on random instances in tests/test_divergences.py, explored numerically in
scripts/renyi_threshold.py): for every positive trace-preserving map Phi... more precisely for every CPTP Phi,
~D_alpha(Phi rho||Phi sigma) <= ~D_alpha(rho||sigma) for alpha >= 1/2 [Frank-Lieb 2013, Beigi 2013 (alpha > 1),
MDSFT 2013 and WWY 2014 (1 < alpha <= 2)], and D_alpha(Phi rho||Phi sigma) <= D_alpha(rho||sigma) for alpha in
(0, 2] [Petz 1986]. Both fail in general outside these ranges.

References
----------
D. Petz, Quasi-entropies for finite quantum systems, Rep. Math. Phys. 23, 57 (1986).
M. Mueller-Lennert, F. Dupuis, O. Szehr, S. Fehr, M. Tomamichel, On quantum Renyi entropies: a new generalization
    and some properties, J. Math. Phys. 54, 122203 (2013).                                            [MDSFT 2013]
M. M. Wilde, A. Winter, D. Yang, Strong converse for the classical capacity of entanglement-breaking and Hadamard
    channels via a sandwiched Renyi relative entropy, Commun. Math. Phys. 331, 593 (2014).            [WWY 2014]
R. L. Frank, E. H. Lieb, Monotonicity of a relative Renyi entropy, J. Math. Phys. 54, 122201 (2013).
S. Beigi, Sandwiched Renyi divergence satisfies data processing inequality, J. Math. Phys. 54, 122202 (2013).
M. Junge, R. Renner, D. Sutter, M. M. Wilde, A. Winter, Universal recovery maps and approximate sufficiency of
    quantum relative entropy, Ann. Henri Poincare 19, 2955 (2018).                                  [JRSWW 2018]
H. Kosaki, Applications of the complex interpolation method to a von Neumann algebra: non-commutative L^p-spaces,
    J. Funct. Anal. 56, 29 (1984).
N. Datta, Min- and max-relative entropies and a new entanglement monotone, IEEE Trans. Inf. Theory 55, 2816 (2009).
M. Tomamichel, Quantum Information Processing with Finite Resources, Springer (2016).
"""
import math
from typing import List, NamedTuple, Optional, Tuple, Union

import torch

from .algebra import TypeIAlgebra
from . import cost
from . import functional as _fc
from .states import DensityMatrix, TRACES, _trace_weights
from .channels import Channel, InterSectorChannel, _same_structure

Operator = TypeIAlgebra.Operator
Alpha = Union[float, torch.Tensor]

__all__ = ["petz_renyi_divergence", "sandwiched_renyi_divergence", "renyi_divergence", "umegaki_relative_entropy",
           "max_relative_entropy", "min_relative_entropy", "sector_contributions", "SectorContributions",
           "PetzRecovery", "petz_recovery", "universal_recovery", "universal_recovery_nodes", "lp_norm",
           "kosaki_norm"]


# ======================================================================
# helpers: inputs, batches, tolerances
# ======================================================================
def _as_blunt(op: Operator, trace: str, what: str) -> Operator:
    """The Tr_blunt density of a positive input given as a density with respect to `trace` (positivity is checked
    unless it is known, which reads one value from the device)."""
    if not isinstance(op, Operator):
        raise TypeError(f"{what}: expected an Operator or DensityMatrix, got {type(op).__name__}")
    if trace not in TRACES:
        raise ValueError(f"trace must be one of {TRACES}")
    op._require_positive(f"{what} must be a positive operator")
    if trace == 'blunt':
        return op
    if isinstance(op, DensityMatrix):
        raise ValueError(f"{what}: a DensityMatrix is a density for Tr_blunt; pass trace='blunt', or wrap "
                         f"rho.density('{trace}') with alg.operator(...) to give a density for '{trace}'")
    w = _trace_weights(op.algebra, trace)
    src = op
    out = Operator(op.algebra, generator=lambda: src.matrix * w[None, :, None, None].to(src.matrix.dtype))
    out._is_self_adjoint = out._is_positive = out._is_normal = True   # positive central multiple of a positive op
    out._batch = src._batch
    return out


def _check_pair(rho: Operator, sigma: Operator):
    if not _same_structure(rho.algebra, sigma.algebra):
        raise ValueError(f"rho and sigma must belong to the same algebra ({rho.algebra} with {rho.algebra.dtype} "
                         f"vs {sigma.algebra} with {sigma.algebra.dtype})")


def _batch(*sizes: int) -> int:
    B = max(sizes)
    if any(s not in (1, B) for s in sizes):
        raise ValueError(f"batch sizes {sizes} do not broadcast (each must be 1 or {B})")
    return B


def _alpha(alpha: Alpha, rdt: torch.dtype, dev) -> Tuple[torch.Tensor, Optional[float]]:
    """alpha as a tensor (1,) or (batch,) of the real dtype, and as a number if it was given as one (checked)."""
    if isinstance(alpha, torch.Tensor):
        if alpha.dim() > 1:
            raise ValueError("alpha must be a number or a tensor of shape (batch,)")
        return alpha.to(device=dev, dtype=rdt).reshape(-1), None
    a = float(alpha)
    if not (a > 0 and math.isfinite(a)) or a == 1.0:
        raise ValueError(f"alpha must lie in (0, 1) or (1, inf), got {alpha} (use umegaki_relative_entropy, "
                         f"max_relative_entropy or min_relative_entropy for the limits, or renyi_divergence)")
    return torch.full((1,), a, dtype=rdt, device=dev), a


def _rtol(alg: TypeIAlgebra, rtol: Optional[float]) -> float:
    return alg.total_subspace_dim * torch.finfo(alg.hilbert.real_dtype).eps if rtol is None else float(rtol)


def _supp_tol(alg: TypeIAlgebra, supp_tol: Optional[float]) -> float:
    return 100 * alg.total_subspace_dim * torch.finfo(alg.hilbert.real_dtype).eps if supp_tol is None \
        else float(supp_tol)


def _bcast(t: torch.Tensor, Bo: int, B: int, n: int) -> torch.Tensor:
    """(Bo * n, ...) -> (B * n, ...) for a stack of the n sectors of one size group of Bo samples (Bo in {1, B})."""
    if Bo == B:
        return t
    return t.reshape(Bo, n, *t.shape[1:]).expand(B, n, *t.shape[1:]).reshape(B * n, *t.shape[1:])


def _blocks(op: Operator, B: int) -> List[Tuple[int, list, torch.Tensor]]:
    """[(k, idx, x)] per block-size group with the active blocks x (B n, k, k) (differentiable), broadcast to B."""
    Bo = op.matrix.shape[0]
    return [(k, idx, _bcast(x, Bo, B, len(idx))) for k, idx, x in op._groups()]


def _decomposed(op: Operator, B: int):
    """[(k, idx, x, w, V)] per block-size group: blocks x (B n, k, k) with gradient and their eigendecomposition
    x = V diag(w) V^* (without gradient; from the spectral cache if there is one), broadcast to B."""
    Bo = op.matrix.shape[0]
    out = []
    for k, idx, x, (w, V) in op._decompositions('hermitian'):
        n = len(idx)
        out.append((k, idx, _bcast(x, Bo, B, n), _bcast(w, Bo, B, n), _bcast(V, Bo, B, n)))
    return out


def _global_max(values: List[torch.Tensor], B: int) -> torch.Tensor:
    """Largest value per sample over the groups (each (B n, k)), clamped at 0: shape (B,)."""
    return torch.stack([v.reshape(B, -1).amax(-1) for v in values], dim=-1).amax(-1).clamp(min=0)


def _rep(t: torch.Tensor, n: int) -> torch.Tensor:
    """Per-sample values (B,) as a column (B n, 1) for the n sectors of a size group."""
    return t[:, None].expand(t.shape[0], n).reshape(-1, 1)


def _support_power(p, thr, scale=None):
    """lambda -> (lambda / scale)^p on lambda > thr, 0 elsewhere (p, thr, scale numbers or columns (B n, 1)); p = 0
    gives the support projection."""
    def f(lam):
        big = lam > thr
        safe = torch.where(big, lam, torch.ones_like(lam) if scale is None else scale.expand_as(lam))
        if scale is not None:
            safe = safe / scale
        return torch.where(big, safe ** p, torch.zeros_like(lam))
    return f


def _support_log(thr):
    def f(lam):
        big = lam > thr
        return torch.where(big, torch.log(torch.where(big, lam, torch.ones_like(lam))), torch.zeros_like(lam))
    return f


def _apply(x, w, V, f) -> torch.Tensor:
    """f(x) for a stack of self-adjoint blocks with known decomposition (gradients by Daleckii-Krein)."""
    return _fc.spectral_apply(x, f, 'hermitian', decomposition=(w, V))


def _projection(w, V, thr) -> torch.Tensor:
    """Support projection V 1[w > thr] V^* (no gradient)."""
    with torch.no_grad():
        return (V * (w > thr).to(V.dtype).unsqueeze(-2)) @ V.conj().transpose(-2, -1)


def _tr_prod(a: torch.Tensor, b: torch.Tensor) -> torch.Tensor:
    """Re Tr(a b) for stacks of blocks, without forming a b."""
    dt = torch.promote_types(a.dtype, b.dtype)
    return (a.to(dt) * b.to(dt).transpose(-2, -1)).sum(dim=(-2, -1)).real


def _tr(x: torch.Tensor) -> torch.Tensor:
    return torch.diagonal(x, dim1=-2, dim2=-1).sum(-1).real


def _scatter(alg: TypeIAlgebra, parts, B: int, dtype) -> torch.Tensor:
    """Per-sector values (B, C) from [(idx, values (B n,))] (zero for empty sectors)."""
    out = torch.zeros(B, alg.C, dtype=dtype, device=alg.hilbert.device)
    for idx, v in parts:
        out[:, alg._sel(idx)] = v.reshape(B, len(idx)).to(dtype)
    return out


def _eigvalsh(M: torch.Tensor) -> torch.Tensor:
    return cost.batched_call(torch.linalg.eigvalsh, 0.5 * (M + M.conj().transpose(-2, -1)),
                             f"eigvalsh (blocks of size {M.shape[-1]})")


# ======================================================================
# per-sector quantities
# ======================================================================
class _Terms(NamedTuple):
    terms: torch.Tensor      # (B, C) Renyi: Q_c exp(-log_scale) (Q_c the per-sector Q_alpha or ~Q_alpha);
                             #        Umegaki: Tr rho_c (log rho_c - log sigma_c)
    log_scale: Optional[torch.Tensor]   # (B,) Renyi: log of the common scale of the terms; Umegaki: None
    leak: torch.Tensor       # (B, C) Tr rho_c (1 - Pi_sigma_c) (no gradient)
    p: torch.Tensor          # (B, C) Tr rho_c
    q: torch.Tensor          # (B, C) Tr sigma_c
    alpha: torch.Tensor      # (B,) (ones for Umegaki)


def _scale(values: List[torch.Tensor], thr: torch.Tensor, B: int, smallest: torch.Tensor) -> torch.Tensor:
    """Per sample (B,): the smallest value above thr (where `smallest`) or the largest value, over all groups
    (each (B n, k)); 1 if there is none (zero operator). Without gradient: a constant rescaling."""
    with torch.no_grad():
        big = torch.stack([v.reshape(B, -1).amax(-1) for v in values], dim=-1).amax(-1)
        small = torch.stack([torch.where(v > _rep(thr, v.shape[0] // B), v, torch.full_like(v, math.inf))
                             .reshape(B, -1).amin(-1) for v in values], dim=-1).amin(-1)
        s = torch.where(smallest, small, big)
        return torch.where(torch.isfinite(s) & (s > 0), s, torch.ones_like(s))


def _terms(rho: Operator, sigma: Operator, alpha: Optional[Alpha], kind: str, trace: str,
           rtol: Optional[float]) -> _Terms:
    """
    Per-sector Renyi / Umegaki quantities; one batched decomposition per block size and operator. The Renyi
    quantities are evaluated with rescaled operators, Q = exp(log_scale) sum_c terms_c, so that large alpha (D_max as
    a limit) and single precision do not overflow: Petz with (rho / a)^alpha (sigma / b)^{1-alpha}, sandwiched with
    ((sigma / b)^g rho (sigma / b)^g / m)^alpha, a = lambda_max(rho), m = lambda_max of the sandwich, and b the
    smallest eigenvalue of sigma on its support for alpha > 1, lambda_max(sigma) for alpha < 1 (every factor <= 1).
    """
    rho = _as_blunt(rho, trace, "rho")
    sigma = _as_blunt(sigma, trace, "sigma")
    _check_pair(rho, sigma)
    alg = rho.algebra
    rdt, dev = alg.hilbert.real_dtype, alg.hilbert.device
    if kind == 'umegaki':
        a = torch.ones(1, dtype=rdt, device=dev)
    else:
        a, _ = _alpha(alpha, rdt, dev)
    B = _batch(rho.matrix.shape[0], sigma.matrix.shape[0], a.shape[0])
    a = a.expand(B)
    rt = _rtol(alg, rtol)
    ds = _decomposed(sigma, B)
    ws = [w for *_, w, _ in ds]
    thr_s = rt * _global_max(ws, B)
    parts, leaks, ps, qs = [], [], [], []
    log_scale = None
    if kind != 'umegaki':
        b = _scale(ws, thr_s, B, a > 1)                       # sigma / b: factors (sigma / b)^{1-alpha}, ^g <= 1

    if kind in ('petz', 'umegaki'):
        dr = _decomposed(rho, B)
        wr_all = [w for *_, w, _ in dr]
        thr_r = rt * _global_max(wr_all, B)
        if kind == 'petz':
            ar = _scale(wr_all, thr_r, B, torch.zeros_like(a, dtype=torch.bool))
            log_scale = a * torch.log(ar) + (1 - a) * torch.log(b)
        for (k, idx, xr, wr, Vr), (_, _, xs, ws_, Vs) in zip(dr, ds):
            n = len(idx)
            tr_, ts = _rep(thr_r, n), _rep(thr_s, n)
            if kind == 'petz':
                an, arn, bn = _rep(a, n), _rep(ar, n), _rep(b, n)
                v = _tr_prod(_apply(xr, wr, Vr, _support_power(an, tr_, arn)),
                             _apply(xs, ws_, Vs, _support_power(1 - an, ts, bn)))
            else:                         # Tr rho log rho - Tr rho log sigma, logarithms on the supports
                ent = _tr(_apply(xr, wr, Vr, lambda lam, t=tr_: torch.where(lam > t, torch.special.xlogy(lam, lam),
                                                                            torch.zeros_like(lam))))
                v = ent - _tr_prod(xr, _apply(xs, ws_, Vs, _support_log(ts)))
            parts.append((idx, v))
    elif kind == 'sandwiched':
        # the nonzero spectrum of S rho S (S = (sigma/b)^g) is that of B B^* with B = S rho^{1/2}: its eigenvalues
        # are the squared singular values of B, which the SVD resolves down to (eps ||S|| ||rho^{1/2}||)^2 instead of
        # eps ||S rho S|| for the eigenvalues of the product (they matter: lambda^alpha is not small for small alpha)
        dr = _decomposed(rho, B)
        thr_r = rt * _global_max([w for *_, w, _ in dr], B)
        svs, smax, hmax = [], [], []
        for (k, idx, xr, wr, Vr), (_, _, xs, ws_, Vs) in zip(dr, ds):
            n = len(idx)
            an, bn = _rep(a, n), _rep(b, n)
            f = _support_power((1 - an) / (2 * an), _rep(thr_s, n), bn)
            S = _apply(xs, ws_, Vs, f)
            H = _apply(xr, wr, Vr, _support_power(0.5, _rep(thr_r, n)))
            svs.append(_fc.singular_values(S @ H.to(S.dtype)))
            smax.append(f(ws_).detach())
            hmax.append(wr.clamp(min=0).sqrt())
        m = _global_max([sv.detach() for sv in svs], B)
        m = torch.where(m > 0, m, torch.ones_like(m))
        # singular values below the rounding level of B = S rho^{1/2} count as zeros
        thr_sv = rt * torch.maximum(m, _global_max(smax, B) * _global_max(hmax, B))
        log_scale = (1 - a) * torch.log(b) + 2 * a * torch.log(m)
        for (k, idx, x, w, V), sv in zip(ds, svs):
            n = len(idx)
            parts.append((idx, _support_power(2 * _rep(a, n), _rep(thr_sv, n), _rep(m, n))(sv).sum(-1)))
    else:
        raise ValueError(f"unknown kind {kind!r}")

    for (k, idx, xr), (_, _, xs, ws_, Vs) in zip(_blocks(rho, B), ds):
        n = len(idx)
        xd = xr.detach()
        leaks.append((idx, _tr(xd) - _tr_prod(xd, _projection(ws_, Vs, _rep(thr_s, n)))))
        ps.append((idx, _tr(xr)))
        qs.append((idx, _tr(xs)))
    return _Terms(_scatter(alg, parts, B, rdt), log_scale, _scatter(alg, leaks, B, rdt), _scatter(alg, ps, B, rdt),
                  _scatter(alg, qs, B, rdt), a)


def _from_terms(t: _Terms, supp_tol: float) -> torch.Tensor:
    """Divergence (B,) from the per-sector quantities, with the support conventions of the module docstring."""
    Q = t.terms.sum(-1)
    tr = t.p.sum(-1)
    tr_safe = torch.where(tr > 0, tr, torch.ones_like(tr))
    leak = t.leak.sum(-1) > supp_tol * tr
    inf = torch.full_like(Q, math.inf)
    if t.log_scale is None:                                  # Umegaki
        return torch.where(leak, inf, Q / tr_safe)
    a = t.alpha
    pos = Q > 0
    D = (t.log_scale + torch.log(torch.where(pos, Q, torch.ones_like(Q))) - torch.log(tr_safe)) / (a - 1)
    D = torch.where(pos, D, inf)
    return torch.where(leak & (a > 1), inf, D)


# ======================================================================
# divergences
# ======================================================================
def petz_renyi_divergence(rho: Operator, sigma: Operator, alpha: Alpha, trace: str = 'blunt',
                          rtol: Optional[float] = None, supp_tol: Optional[float] = None) -> torch.Tensor:
    """
    Petz-Renyi divergence D_alpha(rho||sigma) = (alpha-1)^{-1} log(Tr rho^alpha sigma^{1-alpha} / Tr rho), shape
    (batch,).

    rho, sigma: states (DensityMatrix) or positive Operators of the same algebra, batch 1 or B.
    alpha: a number in (0, 1) u (1, inf) (checked), or a tensor (batch,) of such values, one per sample (not checked:
        alpha = 1 gives NaN, see renyi_divergence for the limits).
    trace: the trace for which plain Operator inputs are densities (module docstring; the value does not depend on it).
    rtol, supp_tol: support thresholds (module docstring). +inf for alpha > 1 if supp rho is not in supp sigma, and
        for alpha < 1 if the supports are orthogonal.
    Differentiable in rho, sigma (and alpha) at points where the supports are stable.
    """
    t = _terms(rho, sigma, alpha, 'petz', trace, rtol)
    return _from_terms(t, _supp_tol(rho.algebra, supp_tol))


def sandwiched_renyi_divergence(rho: Operator, sigma: Operator, alpha: Alpha, trace: str = 'blunt',
                                rtol: Optional[float] = None, supp_tol: Optional[float] = None) -> torch.Tensor:
    """
    Sandwiched Renyi divergence ~D_alpha(rho||sigma) = (alpha-1)^{-1} log(Tr (sigma^g rho sigma^g)^alpha / Tr rho),
    g = (1-alpha)/(2 alpha), shape (batch,). Arguments as for petz_renyi_divergence; ~D_{1/2} = -log F (fidelity).
    """
    t = _terms(rho, sigma, alpha, 'sandwiched', trace, rtol)
    return _from_terms(t, _supp_tol(rho.algebra, supp_tol))


def umegaki_relative_entropy(rho: Operator, sigma: Operator, trace: str = 'blunt', rtol: Optional[float] = None,
                             supp_tol: Optional[float] = None) -> torch.Tensor:
    """
    Umegaki relative entropy D(rho||sigma) = Tr rho (log rho - log sigma) / Tr rho, shape (batch,), with the exact
    support convention (+inf if supp rho is not in supp sigma; 0 log 0 = 0). For states with supp rho in supp sigma
    it agrees with states.relative_entropy, which instead clamps the spectrum of sigma (finite values ~ -log eps).
    """
    t = _terms(rho, sigma, None, 'umegaki', trace, rtol)
    return _from_terms(t, _supp_tol(rho.algebra, supp_tol))


def max_relative_entropy(rho: Operator, sigma: Operator, trace: str = 'blunt', rtol: Optional[float] = None,
                         supp_tol: Optional[float] = None) -> torch.Tensor:
    """
    D_max(rho||sigma) = log min{lambda: rho <= lambda sigma} = log lambda_max(sigma^{-1/2} rho sigma^{-1/2}) (the
    largest over all sectors), +inf if supp rho is not in supp sigma; the limit alpha -> inf of the sandwiched
    divergences. Shape (batch,).
    """
    rho = _as_blunt(rho, trace, "rho")
    sigma = _as_blunt(sigma, trace, "sigma")
    _check_pair(rho, sigma)
    alg = rho.algebra
    rdt = alg.hilbert.real_dtype
    B = _batch(rho.matrix.shape[0], sigma.matrix.shape[0])
    rt = _rtol(alg, rtol)
    ds = _decomposed(sigma, B)
    thr_s = rt * _global_max([w for *_, w, _ in ds], B)
    lmax, leaks, ps = [], [], []
    for (k, idx, xr), (_, _, xs, ws, Vs) in zip(_blocks(rho, B), ds):
        n = len(idx)
        ts = _rep(thr_s, n)
        S = _apply(xs, ws, Vs, _support_power(-0.5, ts))
        lmax.append(_eigvalsh(S @ xr.to(S.dtype) @ S)[..., -1].reshape(B, n).amax(-1))
        xd = xr.detach()
        leaks.append((idx, _tr(xd) - _tr_prod(xd, _projection(ws, Vs, ts))))
        ps.append((idx, _tr(xd)))
    lam = torch.stack(lmax, dim=-1).amax(-1)
    tr = _scatter(alg, ps, B, rdt).sum(-1)
    leak = _scatter(alg, leaks, B, rdt).sum(-1) > _supp_tol(alg, supp_tol) * tr
    pos = lam > 0
    D = torch.log(torch.where(pos, lam, torch.ones_like(lam)))
    return torch.where(leak | ~pos, torch.full_like(D, math.inf), D)


def min_relative_entropy(rho: Operator, sigma: Operator, trace: str = 'blunt',
                         rtol: Optional[float] = None) -> torch.Tensor:
    """
    D_0(rho||sigma) = -log(Tr(Pi_rho sigma) / Tr rho), the limit alpha -> 0 of the Petz divergences (Datta's
    min-relative entropy), +inf if the supports are orthogonal. Shape (batch,). Differentiable in sigma only (the
    support projection Pi_rho is piecewise constant in rho).
    """
    rho = _as_blunt(rho, trace, "rho")
    sigma = _as_blunt(sigma, trace, "sigma")
    _check_pair(rho, sigma)
    alg = rho.algebra
    rdt = alg.hilbert.real_dtype
    B = _batch(rho.matrix.shape[0], sigma.matrix.shape[0])
    dr = _decomposed(rho, B)
    thr_r = _rtol(alg, rtol) * _global_max([w for *_, w, _ in dr], B)
    parts, ps = [], []
    for (k, idx, xr, wr, Vr), (_, _, xs) in zip(dr, _blocks(sigma, B)):
        parts.append((idx, _tr_prod(_projection(wr, Vr, _rep(thr_r, len(idx))), xs)))
        ps.append((idx, _tr(xr.detach())))
    Q = _scatter(alg, parts, B, rdt).sum(-1)
    tr = _scatter(alg, ps, B, rdt).sum(-1)
    pos = Q > 0
    D = -torch.log(torch.where(pos, Q, torch.ones_like(Q)) / torch.where(tr > 0, tr, torch.ones_like(tr)))
    return torch.where(pos, D, torch.full_like(D, math.inf))


def renyi_divergence(rho: Operator, sigma: Operator, alpha: Alpha, kind: str = 'sandwiched', **kwargs) -> torch.Tensor:
    """
    D_alpha (kind='petz') or ~D_alpha (kind='sandwiched'), including the limits for a number alpha: alpha = 1 gives
    the Umegaki relative entropy (both kinds), alpha = inf gives D_max (sandwiched), alpha = 0 gives D_0 (Petz).
    Keyword arguments (trace, rtol, supp_tol) as for petz_renyi_divergence.
    """
    if kind not in ('petz', 'sandwiched'):
        raise ValueError("kind must be 'petz' or 'sandwiched'")
    if not isinstance(alpha, torch.Tensor):
        a = float(alpha)
        if a == 1.0:
            return umegaki_relative_entropy(rho, sigma, **kwargs)
        if a == math.inf and kind == 'sandwiched':
            return max_relative_entropy(rho, sigma, **kwargs)
        if a == 0.0 and kind == 'petz':
            kwargs.pop('supp_tol', None)
            return min_relative_entropy(rho, sigma, **kwargs)
    f = petz_renyi_divergence if kind == 'petz' else sandwiched_renyi_divergence
    return f(rho, sigma, alpha, **kwargs)


class SectorContributions(NamedTuple):
    """
    Per-sector decomposition of a divergence of states on M = (+)_c M_{k_c}; all fields (batch, C) except alpha.

    p, q:        sector weights Tr rho_c, Tr sigma_c.
    terms:       Renyi kinds: Q_c = Tr rho_c^alpha sigma_c^{1-alpha} (Petz) or Tr (sigma_c^g rho_c sigma_c^g)^alpha
                 (sandwiched), so that Tr rho^alpha sigma^{1-alpha} = sum_c Q_c (and likewise ~Q); these may
                 overflow for very large alpha, the conditional divergences are computed in logarithms and do not;
                 Umegaki: Tr rho_c (log rho_c - log sigma_c), summing to Tr rho D(rho||sigma).
    conditional: D(rho_c / p_c || sigma_c / q_c), the divergence of the conditional states of sector c: NaN where
                 p_c = 0 (rounding level: p_c <= rtol Tr rho), +inf where q_c = 0 < p_c or where the support
                 conventions of the module give +inf for the sector.
    alpha:       (batch,).
    Identities (sectors with p_c, q_c > 0):  Q_c = p_c^alpha q_c^{1-alpha} exp((alpha-1) conditional_c), hence
    D_alpha(rho||sigma) = (alpha-1)^{-1} log sum_c p_c^alpha q_c^{1-alpha} exp((alpha-1) D_alpha(rho_c/p_c||sigma_c/q_c))
    for states, and for Umegaki terms_c = p_c log(p_c/q_c) + p_c conditional_c, i.e. D(rho||sigma) = D(p||q) +
    sum_c p_c D(rho_c/p_c||sigma_c/q_c).
    """
    p: torch.Tensor
    q: torch.Tensor
    terms: torch.Tensor
    conditional: torch.Tensor
    alpha: torch.Tensor


def sector_contributions(rho: Operator, sigma: Operator, alpha: Optional[Alpha] = None, kind: str = 'petz',
                         trace: str = 'blunt', rtol: Optional[float] = None,
                         supp_tol: Optional[float] = None) -> SectorContributions:
    """Per-sector quantities of D_alpha (kind='petz' / 'sandwiched') or of the Umegaki relative entropy
    (kind='umegaki', alpha ignored); see SectorContributions."""
    if kind not in ('petz', 'sandwiched', 'umegaki'):
        raise ValueError("kind must be 'petz', 'sandwiched' or 'umegaki'")
    t = _terms(rho, sigma, alpha, kind, trace, rtol)
    alg = rho.algebra
    rt, st = _rtol(alg, rtol), _supp_tol(alg, supp_tol)
    p, q, a = t.p, t.q, t.alpha[:, None]
    tr_r, tr_s = p.sum(-1, keepdim=True), q.sum(-1, keepdim=True)
    has_p, has_q = p > rt * tr_r, q > rt * tr_s
    ok = has_p & has_q
    one = torch.ones_like(p)
    ps, qs = torch.where(ok, p, one), torch.where(ok, q, one)
    leak = t.leak > st * torch.where(has_p, p, one)
    inf, nan = torch.full_like(p, math.inf), torch.full_like(p, math.nan)
    if kind == 'umegaki':
        terms = t.terms
        cond = (terms - ps * torch.log(ps / qs)) / ps
        cond = torch.where(leak, inf, cond)
    else:                                         # in logarithms: the scaled terms do not overflow
        ls = t.log_scale[:, None]
        terms = t.terms * torch.exp(ls)
        pos = t.terms > 0
        cond = (ls + torch.log(torch.where(pos, t.terms, one)) - a * torch.log(ps) - (1 - a) * torch.log(qs)) / (a - 1)
        cond = torch.where(pos, cond, inf)
        cond = torch.where(leak & (a > 1), inf, cond)
    cond = torch.where(has_p & ~has_q, inf, cond)
    cond = torch.where(has_p, cond, nan)
    return SectorContributions(p, q, terms, cond, t.alpha)


# ======================================================================
# Petz recovery maps
# ======================================================================
def _spectral_data(op: Operator, rtol: float):
    """(decomposition per size group, support threshold (B,)) of a positive operator, for _spectral_matrix."""
    B = op.matrix.shape[0]
    dec = _decomposed(op, B)
    return dec, rtol * _global_max([w for *_, w, _ in dec], B)


def _spectral_matrix(op: Operator, data, f_of_thr) -> torch.Tensor:
    """Padded matrix (B, C, k_max, k_max) of f(op), f = f_of_thr(threshold column) acting on the numerical support;
    differentiable (Daleckii-Krein)."""
    alg = op.algebra
    dec, thr = data
    B = thr.shape[0]
    blocks = [(k, idx, _apply(x, w, V, f_of_thr(_rep(thr, len(idx))))) for k, idx, x, w, V in dec]
    dt = torch.promote_types(op.matrix.dtype, blocks[0][2].dtype)
    out = torch.zeros(B, alg.C, alg.k_max, alg.k_max, dtype=dt, device=op.matrix.device)
    for k, idx, y in blocks:
        out[:, alg._sel(idx), :k, :k] = y.to(dt).reshape(B, len(idx), k, k)
    return out


def _channel_algebras(channel) -> Tuple[TypeIAlgebra, TypeIAlgebra]:
    if isinstance(channel, Channel):
        return channel.algebra, channel.algebra
    if isinstance(channel, InterSectorChannel):
        return channel.algebra_in, channel.algebra_out
    raise TypeError(f"expected a Channel or InterSectorChannel, got {type(channel).__name__}")


def _complex_power(p: float, t: float):
    """thr -> (lambda -> lambda^(p + i t) on lambda > thr, 0 elsewhere); real for t = 0."""
    def make(thr):
        def f(lam):
            big = lam > thr
            safe = torch.where(big, lam, torch.ones_like(lam))
            val = safe ** p if t == 0 else torch.exp(torch.complex(p * torch.log(safe), t * torch.log(safe)))
            return torch.where(big, val, torch.zeros_like(val))
        return f
    return make


def _flag_like(out: Operator, X: Operator) -> Operator:
    """A map that is completely positive and Hermiticity preserving keeps these flags of its argument."""
    if X._is_self_adjoint:
        out._is_self_adjoint = out._is_normal = True
    if X._is_positive:
        out._is_positive = True
    return out


class PetzRecovery:
    """
    Petz recovery map of sigma and Phi, optionally rotated [JRSWW 2018]:

        R^t(X) = sigma^{1/2 + it} Phi^*( Phi(sigma)^{-1/2 - it} X Phi(sigma)^{-1/2 + it} ) sigma^{1/2 - it},

    i.e. R^t = U_{sigma,t} o R o U_{Phi(sigma),-t} with U_{s,t}(X) = s^{it} X s^{-it}; t = 0 is the Petz map.
    Phi: Channel (M -> M) or InterSectorChannel (M_in -> M_out), batch 1 or B; Phi^* its adjoint for Tr_blunt.
    sigma: positive operator (e.g. a DensityMatrix) on M_in, batch 1 or B. The inverse powers act on the numerical
    support of Phi(sigma) (relative threshold rtol, module docstring). t != 0 needs a complex algebra.

    R^t is completely positive; for trace-preserving Phi, Tr R^t(X) = Tr(Pi_{Phi(sigma)} X) (= Tr X if supp X is in
    supp Phi(sigma)) and R^t(Phi(sigma)) = sigma (on supp sigma). Calling the map on an Operator X of M_out returns
    a lazy Operator of M_in (flagged positive if X is; wrap it as a DensityMatrix when X is a state supported in
    supp Phi(sigma)), differentiable in X, sigma and the Kraus operators. sigma and Phi(sigma) are decomposed once
    per map, also for several rotation parameters (method rotated).
    """

    def __init__(self, sigma: Operator, channel, t: float = 0.0, rtol: Optional[float] = None):
        alg_in, alg_out = _channel_algebras(channel)
        if not isinstance(sigma, Operator) or not _same_structure(sigma.algebra, alg_in):
            raise ValueError("sigma must be a positive operator of the input algebra of the channel")
        sigma._require_positive("sigma must be a positive operator")
        self.alg_in, self.alg_out = alg_in, alg_out
        self.channel, self.sigma = channel, sigma
        self.t = self._check_t(t)
        self._dual = channel.adjoint() if isinstance(channel, Channel) else channel.adjoint('blunt')
        self._rtol_in, self._rtol_out = _rtol(alg_in, rtol), _rtol(alg_out, rtol)
        self._data = None
        self._factors = {}

    def _check_t(self, t: float) -> float:
        t = float(t)
        if t != 0 and not self.alg_in.complex_valued:
            raise ValueError("rotated Petz maps (t != 0) need a complex algebra (alg.like(complex_valued=True))")
        return t

    def _prepare(self):
        if self._data is None:
            plain = Operator(self.alg_in, matrix=self.sigma.matrix)     # no trace-preservation check / warning
            plain._is_self_adjoint = plain._is_positive = plain._is_normal = True
            image = self.channel.apply(plain)
            image._is_self_adjoint = image._is_positive = image._is_normal = True   # Phi is CP
            self._image = image
            self._data = (_spectral_data(self.sigma, self._rtol_in), _spectral_data(image, self._rtol_out))
        return self._data

    def factors(self, t: Optional[float] = None) -> Tuple[torch.Tensor, torch.Tensor]:
        """(sigma^{1/2 + it}, Phi(sigma)^{-1/2 - it}) as padded tensors (B, C, k, k), computed once per t."""
        t = self.t if t is None else self._check_t(t)
        if t not in self._factors:
            ds, di = self._prepare()
            self._factors[t] = (_spectral_matrix(self.sigma, ds, _complex_power(0.5, t)),
                                _spectral_matrix(self._image, di, _complex_power(-0.5, -t)))
        return self._factors[t]

    def _apply_matrix(self, x: torch.Tensor, t: float) -> torch.Tensor:
        A, Bm = self.factors(t)
        dt = torch.promote_types(x.dtype, Bm.dtype)
        Bm = Bm.to(dt)
        Y = Operator(self.alg_out, matrix=Bm @ x.to(dt) @ Bm.conj().transpose(-2, -1))
        Z = self._dual.apply(Y).matrix
        dt = torch.promote_types(Z.dtype, A.dtype)
        A = A.to(dt)
        return (A @ Z.to(dt) @ A.conj().transpose(-2, -1)).to(self.alg_in.dtype)

    def rotated(self, X: Operator, t: float) -> Operator:
        """R^t(X) for another rotation parameter t (same decompositions)."""
        t = self._check_t(t)
        self._check_arg(X)
        return _flag_like(Operator(self.alg_in, generator=lambda: self._apply_matrix(X.matrix, t)), X)

    def _check_arg(self, X):
        if not isinstance(X, Operator) or not _same_structure(X.algebra, self.alg_out):
            raise ValueError("the argument must be an operator of the output algebra of the channel")

    def __call__(self, X: Operator) -> Operator:
        return self.rotated(X, self.t)


def petz_recovery(sigma: Operator, channel, X: Operator, t: float = 0.0, rtol: Optional[float] = None) -> Operator:
    """R^t_{sigma,Phi}(X) (PetzRecovery; t = 0: the Petz recovery map)."""
    return PetzRecovery(sigma, channel, t, rtol)(X)


def universal_recovery_nodes(step: float = 0.125, cutoff: float = 12.0) -> Tuple[List[float], List[float]]:
    """Nodes t_j and weights of the trapezoidal rule for int beta_0(t) g(t) dt on [-cutoff, cutoff] (weights
    renormalised to sum 1; beta_0(t) = pi/2 (cosh(pi t) + 1)^{-1})."""
    n = int(math.floor(cutoff / step))
    ts = [j * step for j in range(-n, n + 1)]
    wts = [0.5 * math.pi / (math.cosh(math.pi * t) + 1) for t in ts]
    total = sum(wts)
    return ts, [w / total for w in wts]


def universal_recovery(sigma: Operator, channel, X: Operator, step: float = 0.125, cutoff: float = 12.0,
                       rtol: Optional[float] = None) -> Operator:
    """
    Universal recovery map of [JRSWW 2018], R(X) = int dt beta_0(t) R^{t/2}(X), beta_0(t) = pi/2 (cosh(pi t) + 1)^{-1}
    (a probability density), with R^t the rotated Petz maps of PetzRecovery. The integral is computed by the
    trapezoidal rule with the given step on [-cutoff, cutoff] (universal_recovery_nodes): the integrand is analytic in
    the strip |Im t| < 1, so the discretisation error decays like exp(-2 pi / step) times a factor that grows with the
    condition of sigma and Phi(sigma), and the truncation error is about 2 pi exp(-pi cutoff) (~1e-15 with the
    defaults). sigma and Phi(sigma) are decomposed once. Complex algebras only. Returns a lazy Operator of M_in
    (positive if X is).
    """
    R = PetzRecovery(sigma, channel, 0.0, rtol)
    if not R.alg_in.complex_valued:
        raise ValueError("the universal recovery map needs a complex algebra (alg.like(complex_valued=True))")
    R._check_arg(X)
    ts, wts = universal_recovery_nodes(step, cutoff)

    def generator():
        x = X.matrix
        acc = None
        for t, w in zip(ts, wts):                       # one term per node (no loop over sectors)
            y = R._apply_matrix(x, t / 2) * w
            acc = y if acc is None else acc + y
            R._factors.clear()                          # the factors of a node are not needed again
        return acc
    return _flag_like(Operator(R.alg_in, generator=generator), X)


# ======================================================================
# noncommutative L^p norms
# ======================================================================
def _weights(alg: TypeIAlgebra, trace: str, weights: Optional[torch.Tensor]) -> torch.Tensor:
    """Sector weights (1 or B, C) of the trace (or the explicit weights, non-negative; not checked)."""
    rdt, dev = alg.hilbert.real_dtype, alg.hilbert.device
    if weights is None:
        return _trace_weights(alg, trace)[None, :]
    w = torch.as_tensor(weights).to(device=dev, dtype=rdt)
    if w.dim() == 1:
        w = w[None, :]
    if w.dim() != 2 or w.shape[-1] != alg.C:
        raise ValueError(f"weights must have shape ({alg.C},) or (batch, {alg.C})")
    return w


def _check_p(p: float) -> float:
    p = float(p)
    if not p >= 1:
        raise ValueError(f"p must lie in [1, inf], got {p}")
    return p


def _schatten(alg: TypeIAlgebra, svals, w: torch.Tensor, p: float, B: int) -> torch.Tensor:
    """(sum_c w_c sum_i s_ci^p)^{1/p} from singular values [(idx, s (B n, k))], scaled by the largest singular value
    against overflow; p = inf: max over sectors with w_c > 0."""
    rdt = alg.hilbert.real_dtype
    smax = _scatter(alg, [(idx, s.amax(-1)) for idx, s in svals], B, rdt)          # (B, C)
    if p == math.inf:
        return torch.where(w > 0, smax, torch.zeros_like(smax)).amax(-1)
    m = smax.amax(-1).clamp(min=torch.finfo(rdt).tiny)                              # (B,)
    sums = _scatter(alg, [(idx, ((s / _rep(m, len(idx))) ** p).sum(-1)) for idx, s in svals], B, rdt)
    return m * (w * sums).sum(-1) ** (1 / p)


def _svals(x: torch.Tensor, self_adjoint: bool) -> torch.Tensor:
    if self_adjoint:
        return _eigvalsh(x).abs()
    return _fc.singular_values(x)


def lp_norm(A: Operator, p: float, trace: str = 'blunt', weights: Optional[torch.Tensor] = None) -> torch.Tensor:
    """
    Noncommutative L^p norm ||A||_{p,w} = (sum_c w_c Tr|A_c|^p)^{1/p} for p in [1, inf), ||A||_inf = max_c ||A_c||_op
    over the sectors with w_c > 0; shape (batch,). The weights are those of `trace` ('blunt': 1, 'norm': 1/k_c,
    'tau_vN': 1/(C k_c)) or explicit non-negative `weights` of shape (C,) or (batch, C). Singular values from the
    Hermitian eigensolver if A is flagged self-adjoint, otherwise from the batched SVD; differentiable (at
    non-degenerate singular values).
    """
    p = _check_p(p)
    alg = A.algebra
    B = A.matrix.shape[0]
    w = _weights(alg, trace, weights)
    _batch(B, w.shape[0])
    sa = A._is_self_adjoint is True
    svals = [(idx, _svals(x, sa)) for k, idx, x in A._groups()]
    return _schatten(alg, svals, w, p, B)


def kosaki_norm(x: Operator, rho: Operator, p: float, eta: float = 0.5, trace: str = 'blunt',
                rtol: Optional[float] = None) -> torch.Tensor:
    """
    State-weighted (Kosaki) L^p norm ||x||_{p,rho,eta} = ||d^{eta/p} x d^{(1-eta)/p}||_{p,trace} of x in M, for the
    density d of the state with respect to `trace` (rho itself for trace='blunt'); eta = 1/2: symmetric embedding
    (Beigi 2013), eta = 1: left, eta = 0: right. Factors with exponent 0 are omitted (identity). The value does not
    depend on the trace: the central weights cancel (module docstring). For a faithful state, ||x||_{p} <= ||x||_inf
    and the norms increase with p. Shape (batch,); differentiable in x and rho.
    """
    p = _check_p(p)
    if not 0 <= eta <= 1:
        raise ValueError("eta must lie in [0, 1]")
    rho = _as_blunt(rho, trace, "rho")
    if not _same_structure(x.algebra, rho.algebra):
        raise ValueError("x and rho must belong to the same algebra")
    alg = rho.algebra
    B = _batch(x.matrix.shape[0], rho.matrix.shape[0])
    if p == math.inf:
        return lp_norm(x, p) if B == x.matrix.shape[0] else lp_norm(x, p).expand(B)
    rt = _rtol(alg, rtol)
    dr = _decomposed(rho, B)
    thr = rt * _global_max([w for *_, w, _ in dr], B)
    svals = []
    for (k, idx, xr, wr, Vr), (_, _, xx) in zip(dr, _blocks(x, B)):
        n = len(idx)
        y = xx
        if eta > 0:
            L = _apply(xr, wr, Vr, _support_power(eta / p, _rep(thr, n)))
            y = L.to(torch.promote_types(L.dtype, y.dtype)) @ y
        if eta < 1:
            R = _apply(xr, wr, Vr, _support_power((1 - eta) / p, _rep(thr, n)))
            y = y @ R.to(y.dtype)
        svals.append((idx, _fc.singular_values(y)))
    return _schatten(alg, svals, torch.ones(1, alg.C, dtype=alg.hilbert.real_dtype, device=alg.hilbert.device), p, B)
