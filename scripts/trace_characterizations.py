#!/usr/bin/env python3
"""
Section "Use case: quantitative failure of trace characterisations" of the paper.

Many inequalities that hold for the trace characterise it: a positive (normal) functional phi on M_n, or on a von
Neumann algebra, satisfies them for all operators if and only if phi is a trace (Gardner 1979; Stolyarov, Tikhonov and
Sherstnev 2002; Bikchentaev and coauthors 2005-2026; catalogue: man/article/usecases/trace_characterizations_catalogue.md).
These theorems are yes/no statements. Here we ask how strongly an inequality LHS <= RHS fails for the faithful
functional phi_W(X) = Tr(W X), W > 0, as a function of the Michelson contrast delta = Delta(W) of its density:

    V(delta) = sup { (LHS - RHS) / N : Delta(W) = delta, all operator variables free },

N a natural normalisation per inequality (so that V does not depend on the scale of W or of the variables). delta = 0
is the normalised trace, where every inequality holds and V(0) = 0.

The density W is diagonal with the spectrum of optimize.SpectrumParam (the exact contrast parametrisation: lambda_max =
1, lambda_min = (1 - delta)/(1 + delta), the other eigenvalues free in between), normalised to Tr W = 1 (to Tr W = n for
the determinant inequality, which needs phi(1) = n); the operator variables are parametrised by optimize.PositiveParam
(free contrast in [0, 1], i.e. every positive operator up to scale, rank-deficient ones included), SelfAdjointParam,
UnitaryParam (polar factors, and projections U diag(1, .., 1, 0, .., 0) U^* of every rank) and positive scalars. One
batch holds the delta grid times the starts; optimize.extremize maximises. Every iterate is an admissible point, so the
computed V is a lower bound for the envelope; it is compared with the closed forms that we prove (paper, Propositions),
with an independent SciPy search (tests), and the best final point of every grid value is recomputed in mpmath, after
re-orthonormalising the unitaries there (certified violation at an exactly admissible point, with the contrast of W
recomputed from its eigenvalues).

Sectors: on M = (+)_c M_{k_c} every trace is sum_c t_c Tr_c, so a density W that is scalar on every block gives a trace
and V must vanish (except for the determinant inequality, which characterises the *canonical* trace among functionals
with phi(1) = n); for the other inequalities the violation is governed by the intra-block contrasts of W. Both are
checked with fixed block densities on M_2 (+) M_3.

Writes generated/trace_characterizations.{json,tex} and figures/trace_characterizations.pdf.

    python scripts/trace_characterizations.py --mode full --device cpu
"""
import math
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
for _p in (str(ROOT), str(ROOT / 'scripts')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from torch_vn_algebra import TypeIAlgebra, cost, optimize as opt   # noqa: E402

DT = torch.complex128


# ----------------------------------------------------------------------------------------------------------------
# parametrised variables
# ----------------------------------------------------------------------------------------------------------------
class ScalarParam(opt.Param):
    """Positive numbers t = exp(s), one per start (relative scale of two operator variables)."""

    def __init__(self, B, dev):
        self.s = (0.5 * torch.randn(B, dtype=torch.float64, device=dev)).requires_grad_(True)

    def parameters(self):
        return [self.s]

    def value(self):
        return torch.exp(self.s)


STRUCTURED = 0.5          # fraction of starts at structured points (optimize.structured_unitary 'pairs', two-level spectra)


class Var:
    """An operator variable: `params` (for extremize), `op()` (an Operator), `export(i)` (numpy data of start i)."""
    kind = None


class PosVar(Var):
    """Positive X = U diag(lam) U^*, lam_max = 1 and lam in [lo, 1] with any contrast in [0, dmax] (lam = 0 allowed
    for dmax = 1): every positive operator up to scale."""
    kind = 'pos'

    def __init__(self, alg, B, dmax=1.0):
        self.alg = alg
        self.p = opt.PositiveParam(alg, B, delta=None, delta_range=(0.0, dmax), init='pairs',
                                   init_fraction=STRUCTURED, spectrum_init='pairs')
        self.params = [self.p]

    def op(self):
        return self.p.operator()

    def lam(self):
        return self.p.spectrum.eigenvalues()

    def fn(self, f):
        """f(X) = U f(lam) U^* without a decomposition (f applied to the parametrised spectrum)."""
        U = self.p.unitary.matrix()
        lam = self.lam()
        mask = self.p.spectrum.mask
        fl = torch.where(mask, f(torch.where(mask, lam, torch.ones_like(lam))), torch.zeros_like(lam))
        return self.alg.operator((U * fl.to(U.dtype).unsqueeze(-2)) @ U.conj().transpose(-2, -1),
                                 is_self_adjoint=True)

    def export(self, i):
        return dict(kind='pos', U=_np(self.p.unitary.matrix()[i]), lam=_np(self.lam()[i]))


class ArbVar(Var):
    """Arbitrary X = U P (polar decomposition, |X| = P exactly): U unitary, P a PosVar."""
    kind = 'arb'

    def __init__(self, alg, B, left=True):
        self.P = PosVar(alg, B)
        self.U = opt.UnitaryParam(alg, B, init='pairs', init_fraction=self.P.p.structured)
        self.left = left                    # X = U P (|X| = P) or X = P U (|X^*| = P)
        self.params = self.P.params + [self.U]

    def op(self):
        u, p = self.U.operator(), self.P.op()
        return u @ p if self.left else p @ u

    def export(self, i):
        d = self.P.export(i)
        return dict(kind='arb', V=_np(self.U.matrix()[i]), U=d['U'], lam=d['lam'], left=self.left)


class SAVar(Var):
    """Self-adjoint X = U diag(s lam) U^* (optimize.SelfAdjointParam, free contrast of |X| in [0, 1])."""
    kind = 'sa'

    def __init__(self, alg, B):
        self.p = opt.SelfAdjointParam(alg, B, delta=None, delta_range=(0.0, 1.0), init='pairs',
                                      init_fraction=STRUCTURED, spectrum_init='pairs')
        self.params = [self.p]

    def op(self):
        return self.p.operator()

    def abs(self):
        return self.p.abs_operator()

    def export(self, i):
        lam = self.p.spectrum.eigenvalues()[i] * self.p.signs[i]
        return dict(kind='sa', U=_np(self.p.unitary.matrix()[i]), lam=_np(lam))


class ProjVar(Var):
    """Projections U D U^*, D = diag(1, .., 1, 0, .., 0) with a random rank r_c in {0, .., k_c} per start and block
    (fixed during the search; the multistart covers the ranks); with `nonzero`, r_0 >= 1."""
    kind = 'proj'

    def __init__(self, alg, B, nonzero=True):
        self.alg = alg
        self.U = opt.UnitaryParam(alg, B, init='pairs', init_fraction=STRUCTURED)
        dev = alg.hilbert.device
        k = torch.tensor(alg.k_factors, device=dev)
        r = torch.floor(torch.rand(B, alg.C, device=dev) * (k + 1)).long().clamp(max=k)
        self.nonzero = nonzero
        if nonzero:
            r[:, 0] = torch.clamp(r[:, 0], min=1)
        self.set_rank(r)
        self.params = [self.U]

    def set_rank(self, r):
        self.rank = r
        self.d = (torch.arange(self.alg.k_max, device=r.device) < r[..., None]).to(torch.float64)   # (B, C, k_max)

    def op(self):
        U = self.U.matrix()
        return self.alg.operator((U * self.d.to(U.dtype).unsqueeze(-2)) @ U.conj().transpose(-2, -1),
                                 is_self_adjoint=True)

    def export(self, i):
        return dict(kind='proj', U=_np(self.U.matrix()[i]), lam=_np(self.d[i]))


def cover_ranks(projs, starts):
    """On a factor M_k: give the starts of every grid value distinct rank patterns of the projections (all patterns
    when `starts` is at least their number, otherwise a random subset without repetition)."""
    import itertools
    alg = projs[0].alg
    if alg.C != 1:
        return
    k = alg.k_factors[0]
    combos = list(itertools.product(*[range(1 if p.nonzero else 0, k + 1) for p in projs]))
    perm = torch.randperm(len(combos)).tolist()
    B = projs[0].rank.shape[0]
    pick = [combos[perm[(b % starts) % len(combos)]] for b in range(B)]
    for j, p in enumerate(projs):
        p.set_rank(torch.tensor([[c[j]] for c in pick], dtype=torch.long, device=p.rank.device))


def _np(t):
    return t.detach().cpu().numpy()


# ----------------------------------------------------------------------------------------------------------------
# the functional phi_W and the density W
# ----------------------------------------------------------------------------------------------------------------
class Weights:
    """Diagonal densities W (eigenvalues w, (B, C, k_max), zero on the padding), with Delta(W) = delta prescribed per
    start (optimize.SpectrumParam; the other eigenvalues are optimised too) or fixed; normalised to Tr W = `total`."""

    def __init__(self, alg, B, delta=None, fixed=None, total=1.0, starts=1):
        self.total = total
        self.starts = starts                        # starts per grid value (consecutive batch elements)
        if fixed is None:
            self.spec = opt.SpectrumParam(alg, B, delta=delta)
            self.params = [self.spec]
        else:
            self.spec = None
            self.fixed = fixed.to(device=alg.hilbert.device, dtype=torch.float64)
            self.params = []

    def w(self):
        lam = self.spec.eigenvalues() if self.spec is not None else self.fixed
        return self.total * lam / lam.sum((-2, -1), keepdim=True)


def phi(op, w):
    """phi_W(X) = Tr(W X) = sum_c sum_i w_ci X_c,ii (complex, (B,))."""
    return (torch.diagonal(op.matrix, dim1=-2, dim2=-1) * w).sum((-2, -1))


def ratio(num, den, tiny=1e-14):
    """num / den, and 0 where den <= tiny (degenerate points such as P = Q = R = 1 in the projection metric, where
    both sides vanish), with finite gradients (no 0/0 in either branch)."""
    ok = den > tiny
    return torch.where(ok, num / torch.where(ok, den, torch.ones_like(den)), torch.zeros_like(num))


def sa(alg, m):
    return alg.operator(m, is_self_adjoint=True)


# ----------------------------------------------------------------------------------------------------------------
# the inequalities: value (LHS - RHS) / N in torch and in mpmath
# ----------------------------------------------------------------------------------------------------------------
def _k1(alg, B, W):
    X = ArbVar(alg, B)

    def f():
        w = W.w()
        p = phi(X.P.op(), w).real
        return (phi(X.op(), w).abs() - p) / p
    return f, [X]


def _k3(alg, B, W):
    A = ArbVar(alg, B, left=False)            # A = P U, |A^*| = P

    def f():
        w = W.w()
        p = phi(A.P.op(), w).real
        return (phi(A.op(), w).real - p) / p      # phi(A + A^*) = 2 Re phi(A) (W self-adjoint)
    return f, [A]


def _k2(alg, B, W):
    X, Y, t = SAVar(alg, B), SAVar(alg, B), ScalarParam(B, alg.hilbert.device)

    def f():
        w, s = W.w(), t.value()
        Z = sa(alg, X.op().matrix + s[:, None, None, None] * Y.op().matrix)
        n = phi(X.abs(), w).real + s * phi(Y.abs(), w).real
        return (phi(Z.abs(), w).real - n) / n
    return f, [X, Y, t]


def _mon(alg, B, W):
    A, C, t = PosVar(alg, B), PosVar(alg, B), ScalarParam(B, alg.hilbert.device)

    def f():
        w, s = W.w(), t.value()
        a = A.op()
        Bm = sa(alg, a.matrix + s[:, None, None, None] * C.op().matrix)
        rhs = phi(Bm @ Bm, w).real
        return (phi(a @ a, w).real - rhs) / rhs
    return f, [A, C, t]


def _met(alg, B, W):
    P, Q, R = ProjVar(alg, B), ProjVar(alg, B, nonzero=False), ProjVar(alg, B)
    cover_ranks([P, Q, R], W.starts)

    def f():
        w = W.w()
        p, q, r = P.op(), Q.op(), R.op()
        rhs = phi(sa(alg, p.matrix - q.matrix).abs(), w).real + phi(sa(alg, q.matrix - r.matrix).abs(), w).real
        return ratio(phi(sa(alg, p.matrix - r.matrix).abs(), w).real - rhs, rhs)
    return f, [P, Q, R]


def _bap(alg, B, W):
    P, Q = ProjVar(alg, B), ProjVar(alg, B)
    cover_ranks([P, Q], W.starts)

    def f():
        w = W.w()
        p, q = P.op(), Q.op()
        lhs = 2 * phi(p @ q, w).real                  # phi(PQ + QP)
        rhs = 2 * phi(q @ p @ q, w).real
        return ratio(lhs - rhs, phi(p, w).real + phi(q, w).real)
    return f, [P, Q]


def _sq(alg, B, W):
    A = ArbVar(alg, B)

    def f():
        w = W.w()
        a, p = A.op(), A.P.op()
        return phi(a @ a, w).abs() / phi(p @ p, w).real - 1
    return f, [A]


def _sgm(alg, B, W):
    A, Bv = PosVar(alg, B, dmax=1 - 1e-8), PosVar(alg, B)        # only A is inverted

    def f():
        w = W.w()
        a, b = A.op(), Bv.op()
        ah, aih = A.fn(torch.sqrt), A.fn(lambda x: x ** -0.5)
        m = sa(alg, (ah @ b @ ah).matrix)
        m._is_positive = True
        g = sa(alg, (aih @ m.power(0.5) @ aih).matrix)          # A^{-1} # B
        g._is_positive = True
        gh = g.power(0.5)
        sgm = gh @ a @ gh                                       # A natural B
        return phi(sgm, w).real / torch.sqrt(phi(a, w).real * phi(b, w).real) - 1
    return f, [A, Bv]


def _det(alg, B, W):
    A = PosVar(alg, B, dmax=0.999)
    n = alg.total_subspace_dim

    def f():
        w = W.w()                                               # Tr W = n (phi(1) = n)
        lam = A.lam()
        mask = A.p.spectrum.mask
        logdet = torch.where(mask, torch.log(torch.where(mask, lam, torch.ones_like(lam))),
                             torch.zeros_like(lam)).sum((-2, -1))
        return 1 - phi(A.op(), w).real / (n * torch.exp(logdet / n))
    return f, [A]


def _alt(alg, B, W):
    A, Bv = PosVar(alg, B), PosVar(alg, B)

    def f():
        w = W.w()
        a, b = A.op(), Bv.op()
        bh = Bv.fn(lambda x: torch.sqrt(x.clamp_min(0)))
        m = bh @ a @ bh
        rhs = phi(b @ a @ a @ b, w).real
        return ratio(phi(m @ m, w).real - rhs, rhs)
    return f, [A, Bv]


def _young(alg, B, W):
    A, Bv, t = ArbVar(alg, B), ArbVar(alg, B), ScalarParam(B, alg.hilbert.device)

    def f():
        w, s = W.w(), t.value()
        x = A.op() @ Bv.op().adjoint()                           # A B^*  (B scaled by s below)
        x._is_self_adjoint = False
        x._is_normal = False
        pa, pb = A.P.op(), Bv.P.op()
        rhs = 0.5 * (phi(pa @ pa, w).real + s ** 2 * phi(pb @ pb, w).real)
        return (s * phi(x.abs(), w).real - rhs) / rhs
    return f, [A, Bv, t]


def _pos_sqrt(alg, m):
    """(m)^{1/2} of a positive matrix tensor (B, C, k, k) by the library's functional calculus."""
    op = alg.operator(m, is_self_adjoint=True, is_positive=True)
    return op.power(0.5)


def _expm_sa(alg, m):
    """exp of a self-adjoint matrix tensor by the spectral theorem."""
    return alg.operator(m, is_self_adjoint=True).apply_function(torch.exp)


def _sub(alg, B, W):
    """Subadditivity (Tikhonov 2005) with f(t) = t^{1/2}: phi((A + B)^{1/2}) <= phi(A^{1/2}) + phi(B^{1/2})."""
    A, Bv, t = PosVar(alg, B, dmax=1 - 1e-10), PosVar(alg, B, dmax=1 - 1e-10), ScalarParam(B, alg.hilbert.device)

    def f():
        w, s = W.w(), t.value()
        lhs = phi(_pos_sqrt(alg, A.op().matrix + s[:, None, None, None] * Bv.op().matrix), w).real
        rhs = phi(A.fn(torch.sqrt), w).real + torch.sqrt(s) * phi(Bv.fn(torch.sqrt), w).real
        return ratio(lhs - rhs, rhs)
    return f, [A, Bv, t]


def _ps(alg, B, W):
    """Generalised Powers-Stormer (Hoa-Osaka-Toan 2013) with f(t) = g(t) = t^{1/2}:
    phi(A + B - |A - B|) <= 2 phi(A^{1/4} B^{1/2} A^{1/4})."""
    A, Bv, t = PosVar(alg, B, dmax=1 - 1e-10), PosVar(alg, B, dmax=1 - 1e-10), ScalarParam(B, alg.hilbert.device)

    def f():
        w, s = W.w(), t.value()
        a, b = A.op().matrix, s[:, None, None, None] * Bv.op().matrix
        tot = phi(alg.operator(a + b, is_self_adjoint=True), w).real
        lhs = tot - phi(sa(alg, a - b).abs(), w).real
        a4 = A.fn(lambda x: x ** 0.25)
        rhs = 2 * torch.sqrt(s) * phi(a4 @ Bv.fn(torch.sqrt) @ a4, w).real
        return ratio(lhs - rhs, tot)               # N = phi(A + B): with N = R the ratio is unbounded
    return f, [A, Bv, t]


HOLDER_P = 3.0


def _hold(alg, B, W):
    """Hoelder restricted to projections, p = 3, q = 3/2: |phi(PQ)| <= phi(P)^{1/p} phi(Q)^{1/q}."""
    P, Q = ProjVar(alg, B), ProjVar(alg, B)
    cover_ranks([P, Q], W.starts)

    def f():
        w = W.w()
        p, q = P.op(), Q.op()
        rhs = phi(p, w).real ** (1 / HOLDER_P) * phi(q, w).real ** (1 - 1 / HOLDER_P)
        return ratio(phi(p @ q, w).abs() - rhs, rhs)
    return f, [P, Q]


E_HALF = math.exp(0.5) - 1.0


def _gt(alg, B, W):
    """Golden-Thompson restricted to projections, symmetric form: phi(e^{P+Q}) <= phi(e^{P/2} e^Q e^{P/2});
    e^{P/2} = 1 + (e^{1/2} - 1) P, e^Q = 1 + (e - 1) Q exactly."""
    P, Q = ProjVar(alg, B), ProjVar(alg, B)
    cover_ranks([P, Q], W.starts)

    def f():
        w = W.w()
        p, q = P.op().matrix, Q.op().matrix
        one = alg._const('eye').to(p.dtype)
        eh = alg.operator(one + E_HALF * p, is_self_adjoint=True)
        eq = alg.operator(one + (math.e - 1) * q, is_self_adjoint=True)
        rhs = phi(eh @ eq @ eh, w).real
        lhs = phi(_expm_sa(alg, p + q), w).real
        return ratio(lhs - rhs, rhs)
    return f, [P, Q]


def _pb(alg, B, W):
    """Peierls-Bogoliubov, symmetric form: phi(e^A) exp(phi(e^{A/2} B e^{A/2}) / phi(e^A)) <= phi(e^{A+B}),
    A, B self-adjoint with free scales."""
    A, Bv = SAVar(alg, B), SAVar(alg, B)
    sA, sB = ScalarParam(B, alg.hilbert.device), ScalarParam(B, alg.hilbert.device)

    def f():
        w = W.w()
        a = sA.value()[:, None, None, None] * A.op().matrix
        b = sB.value()[:, None, None, None] * Bv.op().matrix
        U = A.p.unitary.matrix()
        lam = A.p.spectrum.eigenvalues() * A.p.signs * sA.value()[:, None, None]
        ea = phi(A.p._build(torch.exp(lam) * A.p.spectrum.mask), w).real
        eh = A.p._build(torch.exp(lam / 2) * A.p.spectrum.mask)
        mid = phi(eh @ alg.operator(b, is_self_adjoint=True) @ eh, w).real
        lhs = ea * torch.exp(mid / ea)
        rhs = phi(_expm_sa(alg, a + b), w).real
        return ratio(lhs - rhs, rhs)
    return f, [A, Bv, sA, sB]


# ----------------------------------------------------------------------------------------------------------------
# closed forms (paper, Propositions) and lower-bound constructions; delta = Delta(W), kappa the Kantorovich constant
# ----------------------------------------------------------------------------------------------------------------
def kappa(d):
    """kappa = (a + b)/(2 sqrt(ab)) = 1/sqrt(1 - delta^2) for a/b = (1 + delta)/(1 - delta)."""
    d = np.asarray(d, dtype=float)
    return 1 / np.sqrt(1 - d ** 2)


def kantorovich(d):
    """kappa - 1: Gardner, triangle, A + A^*, |phi(A^2)|, projection metric (sharp, every n), and a lower bound for
    the spectral geometric mean (approached, not attained)."""
    return kappa(d) - 1


def young_exact(d):
    """a/b - 1 = 2 delta/(1 - delta) (sharp, every n)."""
    d = np.asarray(d, dtype=float)
    return 2 * d / (1 - d)


def det_exact(d, n):
    """1 - min over two-level spectra {1 (m times), l (n - m times)} of GM / AM (sharp; n = 2: 1 - sqrt(1 - delta^2))."""
    d = np.atleast_1d(np.asarray(d, dtype=float))
    lo = (1 - d) / (1 + d)
    best = np.zeros_like(d)
    for m in range(1, n):
        best = np.maximum(best, 1 - lo ** ((n - m) / n) / ((m + (n - m) * lo) / n))
    return best


def specht_limit(d):
    """lim_{n -> inf} det_exact = 1 - 1/S(h), S the Specht ratio, h = (1 + delta)/(1 - delta)."""
    d = np.atleast_1d(np.asarray(d, dtype=float))
    out = np.zeros_like(d)
    m = d > 0
    h = (1 + d[m]) / (1 - d[m])
    out[m] = 1 - e_log(h)
    return out


def e_log(h):
    """1/S(h) = e log h / ((h - 1) h^{1/(h-1)})."""
    return math.e * np.log(h) / ((h - 1) * h ** (1 / (h - 1)))


def projections_exact(d):
    """n = 2: s/(1 + s), s = sin(psi/2), delta = sin psi (sharp)."""
    d = np.asarray(d, dtype=float)
    s = np.sqrt((1 - np.sqrt(1 - d ** 2)) / 2)
    return s / (1 + s)


def alt_bound(d):
    """kappa^2 - 1 = delta^2/(1 - delta^2): approached for A = P_u, B -> rank one (every n)."""
    d = np.asarray(d, dtype=float)
    return d ** 2 / (1 - d ** 2)


def monotone_bound(d, grid=20001):
    """max_{0 < y < delta} y^2 (delta - y)^2 / ((1 - 2 delta y + y^2)(1 - y^2)): attained by rank-one A, C
    (every n); = delta^4/16 + O(delta^5)."""
    out = []
    for x in np.atleast_1d(np.asarray(d, dtype=float)):
        if x <= 0:
            out.append(0.0)
            continue
        y = np.linspace(0, x, grid)[1:-1]
        f = y ** 2 * (x - y) ** 2 / ((1 - 2 * x * y + y ** 2) * (1 - y ** 2))
        i = int(np.argmax(f))
        from scipy.optimize import minimize_scalar
        h = y[1] - y[0]
        r = minimize_scalar(lambda t: -(t ** 2 * (x - t) ** 2 / ((1 - 2 * x * t + t ** 2) * (1 - t ** 2))),
                            bounds=(max(y[i] - h, 1e-300), min(y[i] + h, x)), method='bounded',
                            options={'xatol': 1e-14})
        out.append(max(float(f[i]), -float(r.fun)))
    return np.array(out)


# ----------------------------------------------------------------------------------------------------------------
# items 12-16: an independent SciPy search on M_2 (plain NumPy, BFGS from random starts), used as a reference
# ----------------------------------------------------------------------------------------------------------------
def _np_fn(A, f):
    e, V = np.linalg.eigh((A + A.conj().T) / 2)
    return (V * f(e)) @ V.conj().T


def numpy_value(name, w, x):
    """(L - R)/N of items 12-16 on M_2 for W = diag(w) at the real parameter vector x (plain NumPy)."""
    ph = lambda X: float(np.real(np.sum(w * np.diag(X))))                       # noqa: E731
    sq = lambda e: np.sqrt(np.clip(e, 0, None))                                  # noqa: E731
    if name in ('holder', 'gt'):                                                 # rank-one P = uu^*, Q = vv^*
        u = np.array([math.cos(x[0]), np.exp(1j * x[1]) * math.sin(x[0])])
        v = np.array([math.cos(x[2]), np.exp(1j * x[3]) * math.sin(x[2])])
        P, Q = np.outer(u, u.conj()), np.outer(v, v.conj())
        if name == 'holder':
            rhs = ph(P) ** (1 / HOLDER_P) * ph(Q) ** (1 - 1 / HOLDER_P)
            return (abs(np.sum(w * np.diag(P @ Q))) - rhs) / rhs
        I = np.eye(2)
        eh, eq = I + E_HALF * P, I + (math.e - 1) * Q
        rhs = ph(eh @ eq @ eh)
        return (ph(_np_fn(P + Q, np.exp)) - rhs) / rhs
    mats = [(x[8 * i:8 * i + 4] + 1j * x[8 * i + 4:8 * i + 8]).reshape(2, 2) for i in range(2)]
    if name == 'pb':
        A, B = [(m + m.conj().T) / 2 for m in mats]
        ea = ph(_np_fn(A, np.exp))
        eh = _np_fn(A, lambda e: np.exp(e / 2))
        lhs = ea * math.exp(ph(eh @ B @ eh) / ea)
        rhs = ph(_np_fn(A + B, np.exp))
        return (lhs - rhs) / rhs
    A, B = [m @ m.conj().T for m in mats]
    if name == 'subadd':
        rhs = ph(_np_fn(A, sq)) + ph(_np_fn(B, sq))
        return (ph(_np_fn(A + B, sq)) - rhs) / rhs
    if name == 'powers':
        tot = ph(A + B)
        lhs = tot - ph(_np_fn(A - B, np.abs))
        a4 = _np_fn(A, lambda e: np.clip(e, 0, None) ** 0.25)
        return (lhs - 2 * ph(a4 @ _np_fn(B, sq) @ a4)) / tot
    raise KeyError(name)


def scipy_reference(name, delta, starts=40, seed=0):
    """Best value of numpy_value over `starts` BFGS runs from random points (M_2, W = diag(1, l) normalised)."""
    from scipy.optimize import minimize
    rng = np.random.default_rng(seed)
    lo = (1 - delta) / (1 + delta)
    w = np.array([1.0, lo]) / (1 + lo)
    dim = 4 if name in ('holder', 'gt') else 16
    best = -np.inf
    for _ in range(starts):
        r = minimize(lambda x: -numpy_value(name, w, x), rng.normal(size=dim), method='BFGS')
        if np.isfinite(r.fun):
            best = max(best, -float(r.fun))
    return best


# name -> build: the objective; label, norm: LaTeX; ref: bibliography keys; exact(delta, n): sharp value (proved),
# bound(delta, n): proved lower bound (a construction), or None; order: (p, c) with V = c delta^p + o(delta^p);
# kind: 'trace' (every trace sum_c t_c Tr_c of a direct sum satisfies it) or 'canonical' (only the canonical trace)
INEQS = {
    'gardner': dict(build=_k1, label='$|\\varphi(X)|\\le\\varphi(|X|)$', norm='$\\varphi(|X|)$',
                    ref='gardner1979', exact=lambda d, n: kantorovich(d), bound=None, order=(2, 0.5), kind='trace'),
    'triangle': dict(build=_k2, label='$\\varphi(|X+Y|)\\le\\varphi(|X|)+\\varphi(|Y|)$', norm='RHS',
                     ref='gardner1979,stolyarov2002', exact=lambda d, n: kantorovich(d), bound=None, order=(2, 0.5),
                     kind='trace'),
    'realpart': dict(build=_k3, label='$\\varphi(A+A^*)\\le2\\varphi(|A^*|)$', norm='RHS', ref='bikchentaev2019abed',
                     exact=lambda d, n: kantorovich(d), bound=None, order=(2, 0.5), kind='trace'),
    'square': dict(build=_sq, label='$|\\varphi(A^2)|\\le\\varphi(A^*A)$', norm='RHS', ref='petz1988,bikchentaev2026moslehian',
                   exact=lambda d, n: kantorovich(d), bound=None, order=(2, 0.5), kind='trace'),
    'metric': dict(build=_met, label='$\\varphi(|P-R|)\\le\\varphi(|P-Q|)+\\varphi(|Q-R|)$', norm='RHS',
                   ref='bikchentaev2019metric', exact=lambda d, n: kantorovich(d), bound=None, order=(2, 0.5),
                   kind='trace'),
    'sgm': dict(build=_sgm, label='$\\varphi(A\\natural B)\\le\\sqrt{\\varphi(A)\\varphi(B)}$', norm='RHS',
                ref='bikchentaev2026sgm', exact=None, bound=lambda d, n: kantorovich(d), order=(2, 0.5), kind='trace'),
    'alt': dict(build=_alt, label='$\\varphi((B^{1/2}AB^{1/2})^2)\\le\\varphi(BA^2B)$', norm='RHS',
                ref='bikchentaev2011pb', exact=None, bound=lambda d, n: alt_bound(d), order=(2, 1.0), kind='trace'),
    'young': dict(build=_young, label='$\\varphi(|AB^*|)\\le\\frac12\\varphi(|A|^2)+\\frac12\\varphi(|B|^2)$',
                  norm='RHS', ref='bikchentaev2005young', exact=lambda d, n: young_exact(d), bound=None,
                  order=(1, 2.0), kind='trace'),
    'projections': dict(build=_bap, label='$\\varphi(PQ+QP)\\le2\\varphi(QPQ)$',
                        norm='$\\varphi(P)+\\varphi(Q)$', ref='bikchentaev2019abed',
                        exact=lambda d, n: projections_exact(d) if n == 2 else None,
                        bound=lambda d, n: projections_exact(d), order=(1, 0.5), kind='trace'),
    'monotone': dict(build=_mon, label='$0\\le A\\le B\\Rightarrow\\varphi(A^2)\\le\\varphi(B^2)$',
                     norm='$\\varphi(B^2)$', ref='bikchentaev2007mono,sano2006,hoa2010', exact=None,
                     bound=lambda d, n: monotone_bound(d), order=(4, 1 / 16), kind='trace'),
    'det': dict(build=_det, label='$(\\det A)^{1/n}\\le\\frac1n\\varphi(A)$, $\\varphi(1)=n$',
                norm='$(\\det A)^{1/n}$', ref='bikchentaev2020det', exact=lambda d, n: det_exact(d, n), bound=None,
                order=(2, 0.5), kind='canonical'),
    'subadd': dict(build=_sub, label='$\\varphi((A+B)^{1/2})\\le\\varphi(A^{1/2})+\\varphi(B^{1/2})$', norm='RHS',
                   ref='tikhonov2005', exact=None, bound=None, order=None, kind='trace'),
    'powers': dict(build=_ps, label='$\\varphi(A+B-|A-B|)\\le2\\varphi(A^{1/4}B^{1/2}A^{1/4})$',
                   norm='$\\varphi(A+B)$',
                   ref='hoa2013', exact=None, bound=None, order=None, kind='trace'),
    'holder': dict(build=_hold, label='$|\\varphi(PQ)|\\le\\varphi(P)^{1/3}\\varphi(Q)^{2/3}$', norm='RHS',
                   ref='bikchentaev2011mn,bikchentaev2019abed', exact=None, bound=None, order=None, kind='trace'),
    'gt': dict(build=_gt, label='$\\varphi(e^{P+Q})\\le\\varphi(e^{P/2}e^Qe^{P/2})$', norm='RHS',
               ref='bikchentaev2011mn', exact=None, bound=None, order=None, kind='trace'),
    'pb': dict(build=_pb, label='$\\varphi(e^A)e^{\\varphi(e^{A/2}Be^{A/2})/\\varphi(e^A)}\\le\\varphi(e^{A+B})$',
               norm='RHS', ref='bikchentaev2011pb', exact=None, bound=None, order=None, kind='trace'),
}
ORDER = list(INEQS)
INEQS_INDEX = {name: i for i, name in enumerate(ORDER, 1)}


def reference_curve(name, d, n):
    """The proved sharp value if known, otherwise the proved lower bound (None if neither)."""
    e = INEQS[name]['exact'](d, n) if INEQS[name]['exact'] is not None else None
    if e is None and INEQS[name]['bound'] is not None:
        e = INEQS[name]['bound'](d, n)
    return None if e is None else np.asarray(e, dtype=float)


def build(name, alg, B, W):
    f, variables = INEQS[name]['build'](alg, B, W)
    params = [p for v in variables for p in (v.params if isinstance(v, Var) else [v])]
    return f, params + W.params, variables


# ----------------------------------------------------------------------------------------------------------------
# envelopes
# ----------------------------------------------------------------------------------------------------------------
def envelope(name, ks, grid, starts, rounds, steps, seed, dev='cpu', fixed=None, lr=0.05):
    """
    V on the grid of Delta(W) for the algebra (+)_c M_{k_c} (ks), `starts` starts per grid value, all in one batch;
    or, with `fixed` (G, C, k_max) block densities, V for each of these G densities. Returns (V (G,) best over all
    iterates and starts, final values (G, starts), points: exported data of the best final start per grid value).
    """
    torch.manual_seed(seed)
    alg = TypeIAlgebra(ks, ks, precision='double', device=dev)
    G = len(grid) if fixed is None else fixed.shape[0]
    B = G * starts
    total = float(alg.total_subspace_dim) if INEQS[name]['kind'] == 'canonical' else 1.0
    if fixed is None:
        d = torch.as_tensor(np.asarray(grid), dtype=torch.float64, device=dev).repeat_interleave(starts)
        W = Weights(alg, B, delta=d, total=total, starts=starts)
    else:
        W = Weights(alg, B, fixed=fixed.repeat_interleave(starts, 0), total=total, starts=starts)
    f, params, variables = build(name, alg, B, W)
    best = opt.extremize(f, params, maximize=True, rounds=rounds, steps=steps, lr=lr)
    with torch.no_grad():
        final = torch.nan_to_num(f(), nan=-math.inf).reshape(G, starts)
        w = W.w()
    best = torch.nan_to_num(best, nan=-math.inf).reshape(G, starts).max(1).values
    idx = final.argmax(1)
    points = []
    for g in range(G):
        i = g * starts + int(idx[g])
        points.append(dict(w=_np(w[i]), vars=[v.export(i) if isinstance(v, Var) else
                                              dict(kind='scalar', t=float(v.value()[i].detach())) for v in variables]))
    return best.cpu().numpy(), final.cpu().numpy(), points, alg


# ----------------------------------------------------------------------------------------------------------------
# independent recomputation in mpmath (dense matrices on C^n, n = sum k_c)
# ----------------------------------------------------------------------------------------------------------------
def _dense(block, ks):
    """(C, k_max, k_max) blocks -> block-diagonal dense n x n (padding dropped); vectors (C, k_max) -> (n,)."""
    if block.ndim == 2:
        return np.concatenate([block[c, :k] for c, k in enumerate(ks)])
    n = sum(ks)
    out = np.zeros((n, n), dtype=complex)
    o = 0
    for c, k in enumerate(ks):
        out[o:o + k, o:o + k] = block[c, :k, :k]
        o += k
    return out


def certify(name, point, ks, digits=40):
    """(LHS - RHS) / N at an exported point, recomputed in mpmath with `digits` digits; unitaries are re-orthonormalised
    (U (U^* U)^{-1/2}), so the point is exactly admissible. Returns (value, contrast of W, intra-block contrasts)."""
    import mpmath as mp
    with mp.workdps(digits):
        def M(a):
            a = np.asarray(a)
            return mp.matrix([[mp.mpc(complex(v)) for v in row] for row in a])

        def herm(A, f):
            E, Q = mp.eighe((A + A.transpose_conj()) / 2)
            return Q * mp.diag([f(e) for e in E]) * Q.transpose_conj()

        def polar(U):
            return U * herm(U.transpose_conj() * U, lambda e: 1 / mp.sqrt(e))

        def sqrtp(e):
            return mp.sqrt(e) if e > 0 else mp.mpf(0)

        wv = [mp.mpf(float(x)) for x in _dense(point['w'], ks)]
        n = len(wv)

        def ph(A):
            return sum(wv[i] * A[i, i] for i in range(n))

        def rebuild(v):
            if v['kind'] == 'scalar':
                return mp.mpf(v['t'])
            U = polar(M(_dense(v['U'], ks)))
            lam = [mp.mpf(float(x)) for x in _dense(v['lam'], ks)]
            X = U * mp.diag(lam) * U.transpose_conj()
            out = dict(X=X, lam=lam, U=U)
            if v['kind'] == 'arb':
                V = polar(M(_dense(v['V'], ks)))
                out['P'] = X
                out['X'] = V * X if v['left'] else X * V
            return out

        vs = [rebuild(v) for v in point['vars']]
        tiny = mp.mpf(10) ** (-digits // 2)
        absa = lambda A: herm(A, abs)                                       # noqa: E731
        absg = lambda A: herm(A.transpose_conj() * A, sqrtp)                 # noqa: E731
        if name == 'gardner':
            X = vs[0]
            p = mp.re(ph(X['P']))
            val = (abs(ph(X['X'])) - p) / p
        elif name == 'realpart':
            X = vs[0]
            p = mp.re(ph(X['P']))
            val = (mp.re(ph(X['X'])) - p) / p
        elif name == 'triangle':
            X, Y, t = vs
            nrm = mp.re(ph(absa(X['X']))) + t * mp.re(ph(absa(Y['X'])))
            val = (mp.re(ph(absa(X['X'] + t * Y['X']))) - nrm) / nrm
        elif name == 'monotone':
            A, C, t = vs
            Bm = A['X'] + t * C['X']
            rhs = mp.re(ph(Bm * Bm))
            val = (mp.re(ph(A['X'] * A['X'])) - rhs) / rhs
        elif name == 'metric':
            P, Q, R = (v['X'] for v in vs)
            rhs = mp.re(ph(absa(P - Q))) + mp.re(ph(absa(Q - R)))
            val = (mp.re(ph(absa(P - R))) - rhs) / rhs if rhs > mp.mpf(10) ** (-digits // 2) else mp.mpf(0)
        elif name == 'projections':
            P, Q = (v['X'] for v in vs)
            val = (2 * mp.re(ph(P * Q)) - 2 * mp.re(ph(Q * P * Q))) / (mp.re(ph(P)) + mp.re(ph(Q)))
        elif name == 'square':
            A = vs[0]
            val = abs(ph(A['X'] * A['X'])) / mp.re(ph(A['P'] * A['P'])) - 1
        elif name == 'sgm':
            A, B = vs[0]['X'], vs[1]['X']
            ah, aih = herm(A, mp.sqrt), herm(A, lambda e: 1 / mp.sqrt(e))
            g = aih * herm(ah * B * ah, sqrtp) * aih
            gh = herm(g, sqrtp)
            val = mp.re(ph(gh * A * gh)) / mp.sqrt(mp.re(ph(A)) * mp.re(ph(B))) - 1
        elif name == 'det':
            A = vs[0]
            det = mp.re(mp.det(A['X']))
            val = 1 - mp.re(ph(A['X'])) / (n * det ** (mp.mpf(1) / n))
        elif name == 'alt':
            A, B = vs[0]['X'], vs[1]['X']
            bh = herm(B, sqrtp)
            m = bh * A * bh
            rhs = mp.re(ph(B * A * A * B))
            val = mp.re(ph(m * m)) / rhs - 1 if rhs > mp.mpf(10) ** (-digits // 2) else mp.mpf(0)
        elif name == 'young':
            A, B, t = vs
            X = A['X'] * (t * B['X']).transpose_conj()
            rhs = (mp.re(ph(A['P'] * A['P'])) + t ** 2 * mp.re(ph(B['P'] * B['P']))) / 2
            val = (mp.re(ph(absg(X))) - rhs) / rhs
        elif name == 'subadd':
            A, B, t = vs
            rhs = mp.re(ph(herm(A['X'], sqrtp))) + mp.sqrt(t) * mp.re(ph(herm(B['X'], sqrtp)))
            val = (mp.re(ph(herm(A['X'] + t * B['X'], sqrtp))) - rhs) / rhs
        elif name == 'powers':
            A, B, t = vs
            a, b = A['X'], t * B['X']
            a4 = herm(a, lambda e: e ** mp.mpf('0.25') if e > 0 else mp.mpf(0))
            rhs = 2 * mp.re(ph(a4 * herm(b, sqrtp) * a4))
            val = (mp.re(ph(a + b)) - mp.re(ph(absa(a - b))) - rhs) / mp.re(ph(a + b))
        elif name == 'holder':
            P, Q = (v['X'] for v in vs)
            hp = mp.mpf(HOLDER_P)
            rhs = mp.re(ph(P)) ** (1 / hp) * mp.re(ph(Q)) ** (1 - 1 / hp)
            val = (abs(ph(P * Q)) - rhs) / rhs
        elif name == 'gt':
            P, Q = (v['X'] for v in vs)
            I = mp.eye(n)
            eh = I + (mp.exp(mp.mpf(1) / 2) - 1) * P
            eq = I + (mp.e - 1) * Q
            rhs = mp.re(ph(eh * eq * eh))
            val = (mp.re(ph(herm(P + Q, mp.exp))) - rhs) / rhs
        elif name == 'pb':
            A, B, sA, sB = vs
            a, b = sA * A['X'], sB * B['X']
            ea = mp.re(ph(herm(a, mp.exp)))
            eh = herm(a, lambda e: mp.exp(e / 2))
            lhs = ea * mp.exp(mp.re(ph(eh * b * eh)) / ea)
            rhs = mp.re(ph(herm(a + b, mp.exp)))
            val = (lhs - rhs) / rhs
        else:
            raise KeyError(name)
        delta = (max(wv) - min(wv)) / (max(wv) + min(wv))
        intra, o = [], 0
        for k in ks:
            b = wv[o:o + k]
            o += k
            intra.append(float((max(b) - min(b)) / (max(b) + min(b))))
        return float(val), float(delta), intra


# ----------------------------------------------------------------------------------------------------------------
# structure of the extremisers
# ----------------------------------------------------------------------------------------------------------------
def structure(point, ks, tol=1e-4):
    """Rank (eigenvalues above tol of the largest) of every positive / projection variable and of |X|, the relative
    commutator ||[W, X]||_F / (2 ||W||_op ||X||_F) of every operator variable with W, and the eigenvalues of W at the
    ends of its range."""
    w = _dense(point['w'], ks).real
    W = np.diag(w)
    out = dict(ranks=[], comm=[], w_two_level=None)
    for v in point['vars']:
        if v['kind'] == 'scalar':
            continue
        U = _dense(v['U'], ks)
        lam = _dense(v['lam'], ks).real
        X = U @ np.diag(lam) @ U.conj().T
        if v['kind'] == 'arb':
            V = _dense(v['V'], ks)
            X = V @ X if v['left'] else X @ V
        out['ranks'].append(int((np.abs(lam) > tol * max(np.abs(lam).max(), 1e-300)).sum()))
        nx = np.linalg.norm(X)
        out['comm'].append(float(np.linalg.norm(W @ X - X @ W) / (2 * np.abs(w).max() * nx)) if nx > 0 else 0.0)
    span = w.max() - w.min()
    at_end = np.minimum(np.abs(w - w.max()), np.abs(w - w.min())) <= tol * max(span, 1e-300)
    out['w_two_level'] = bool(at_end.all())
    out['w'] = (w / w.max()).tolist()
    return out


# ----------------------------------------------------------------------------------------------------------------
# configurations and the run
# ----------------------------------------------------------------------------------------------------------------
GRID_FULL = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
CONFIGS = {
    # dims: M_n; grid: Delta(W); starts per grid value; rounds x steps of optimize.extremize; sector_starts per fixed
    # block density on M_2 (+) M_3; digits of the mpmath recomputation; workers: processes (CPU, one thread each)
    'quick': dict(dims=[2, 3], grid=[0.0, 0.05, 0.1, 0.2, 0.5, 0.9], starts=8, rounds=5, steps=40, sector_starts=6,
                  digits=30, sector_deltas=[0.5], workers=4, scipy_deltas=[0.5], scipy_starts=6),
    'check': dict(dims=[2, 3, 4], grid=GRID_FULL, starts=4, rounds=2, steps=20, sector_starts=4, digits=30,
                  sector_deltas=[0.3, 0.5, 0.7], workers=4, scipy_deltas=[0.5], scipy_starts=4),
    'full': dict(dims=[2, 3, 4], grid=GRID_FULL, starts=64, rounds=10, steps=80, sector_starts=48, digits=40,
                 sector_deltas=[0.3, 0.5, 0.7], workers=4, scipy_deltas=[0.3, 0.5, 0.9], scipy_starts=40),
}
SECTOR_KS = [2, 3]
N_INDEPENDENT = ('gardner', 'triangle', 'realpart', 'square', 'metric', 'young', 'projections', 'monotone',
                 'subadd', 'holder', 'pb')            # items whose search values do not depend on n (n = 2, 3, 4)
SCIPY_ITEMS = ('subadd', 'powers', 'holder', 'gt', 'pb')     # items without a closed form: SciPy reference on M_2
SECTOR_SCALES = [0.25, 4.0]          # level of the second block relative to lambda_max of the first


def sector_densities(deltas, ks=SECTOR_KS):
    """Fixed block densities on M_2 (+) M_3: (a) central, W = 1 (+) s 1 (a trace for every s, global contrast
    |1 - s|/(1 + s)); (b) diag(1, l) (+) s 1, intra-block contrast delta in the first block only."""
    C, km = len(ks), max(ks)
    rows, meta = [], []
    for s in SECTOR_SCALES:
        w = np.zeros((C, km))
        w[0, :ks[0]] = 1.0
        w[1, :ks[1]] = s
        rows.append(w)
        meta.append(dict(type='central', s=s, delta_intra=0.0, delta=abs(1 - s) / (1 + s)))
    for d in deltas:
        lo = (1 - d) / (1 + d)
        for s in SECTOR_SCALES:
            w = np.zeros((C, km))
            w[0, 0], w[0, 1] = 1.0, lo
            w[1, :ks[1]] = s
            rows.append(w)
            meta.append(dict(type='intra', s=s, delta_intra=d, delta=(max(1, s) - min(lo, s)) / (max(1, s) + min(lo, s))))
    return torch.as_tensor(np.stack(rows)), meta


def det_central(meta, ks=SECTOR_KS):
    """Sharp V of the determinant inequality for a fixed density: 1 - GM/AM of its eigenvalues (with multiplicity)."""
    out = []
    for m in meta:
        lo = (1 - m['delta_intra']) / (1 + m['delta_intra'])
        w = np.array([1.0, lo if m['type'] == 'intra' else 1.0] + [m['s']] * ks[1])
        out.append(1 - np.exp(np.log(w).mean()) / w.mean())
    return out


def _task(job):
    """One envelope (job = ('env', name, n, cfg, seed, dev)) or the sector check (('sec', name, cfg, seed, dev)),
    with the mpmath recomputation of the best final point of every grid value."""
    warnings.simplefilter('ignore', cost.CostWarning)
    warnings.simplefilter('ignore', UserWarning)
    kind, name = job[0], job[1]
    t0 = time.time()
    if kind == 'sci':
        _, _, cfg, seed = job
        warnings.simplefilter('ignore', RuntimeWarning)
        vals = {str(d): scipy_reference(name, d, cfg['scipy_starts'], seed) for d in cfg['scipy_deltas']}
        return dict(job=(kind, name), values=vals, seconds=time.time() - t0)
    if kind == 'env':
        _, _, n, cfg, seed, dev = job
        V, final, points, _ = envelope(name, [n], cfg['grid'], cfg['starts'], cfg['rounds'], cfg['steps'], seed, dev)
        secs = time.time() - t0
        cert = [certify(name, p, [n], cfg['digits']) for p in points]
        ref = reference_curve(name, np.array(cfg['grid']), n)
        return dict(job=(kind, name, n), V_best=V.tolist(), V=final.max(1).tolist(),
                    certified=[c[0] for c in cert], certified_delta=[c[1] for c in cert],
                    reference=None if ref is None else ref.tolist(),
                    structure=[structure(p, [n]) for p in points], points=[_jsonable(p) for p in points],
                    seconds=secs, seconds_certify=time.time() - t0 - secs)
    _, _, cfg, seed, dev = job
    fixed, meta = sector_densities(cfg['sector_deltas'])
    V, final, points, _ = envelope(name, SECTOR_KS, None, cfg['sector_starts'], cfg['rounds'], cfg['steps'], seed,
                                   dev, fixed=fixed)
    cert = [certify(name, p, SECTOR_KS, cfg['digits']) for p in points]
    return dict(job=(kind, name), V_best=V.tolist(), V=final.max(1).tolist(), certified=[c[0] for c in cert],
                intra=[c[2] for c in cert], meta=meta, seconds=time.time() - t0)


def _jsonable(point):
    """An exported point with complex arrays as [real, imag] lists (JSON); from_jsonable inverts it."""
    def conv(x):
        if isinstance(x, np.ndarray):
            return dict(re=x.real.tolist(), im=x.imag.tolist()) if np.iscomplexobj(x) else x.tolist()
        return x
    return dict(w=conv(point['w']), vars=[{k: conv(v) for k, v in var.items()} for var in point['vars']])


def from_jsonable(point):
    def conv(x):
        if isinstance(x, dict) and set(x) == {'re', 'im'}:
            return np.array(x['re']) + 1j * np.array(x['im'])
        return np.array(x) if isinstance(x, list) else x
    return dict(w=conv(point['w']), vars=[{k: conv(v) for k, v in var.items()} for var in point['vars']])


def _init_worker():
    torch.set_num_threads(1)


def run(mode, dev='cpu', seed=0, log=print, names=None, workers=None):
    """All envelopes and sector checks; returns the result dictionary (no files written)."""
    cfg = CONFIGS[mode]
    warnings.simplefilter('ignore', cost.CostWarning)
    names = ORDER if names is None else names
    workers = cfg['workers'] if workers is None else workers
    if not str(dev).startswith('cpu'):
        workers = 1
    jobs = []
    for name in names:
        j = ORDER.index(name)
        for n in cfg['dims']:
            jobs.append(('env', name, n, cfg, seed + 1000 * j + n, dev))
        jobs.append(('sec', name, cfg, seed + 1000 * j + 99, dev))
        if name in SCIPY_ITEMS and cfg.get('scipy_deltas'):
            jobs.append(('sci', name, cfg, seed + 1000 * j + 77))
    jobs.sort(key=lambda jb: -(jb[2] if jb[0] == 'env' else 6 if jb[0] == 'sci' else 5))   # largest first
    t0 = time.time()
    out = []
    if workers > 1:
        import multiprocessing
        from concurrent.futures import ProcessPoolExecutor, as_completed
        with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context('spawn'),
                                 initializer=_init_worker) as ex:
            futs = [ex.submit(_task, jb) for jb in jobs]
            for f in as_completed(futs):
                r = f.result()
                out.append(r)
                _log_task(r, log)
    else:
        for jb in jobs:
            r = _task(jb)
            out.append(r)
            _log_task(r, log)
    res = dict(config={k: v for k, v in cfg.items()}, grid=cfg['grid'], workers=workers, inequalities={})
    for name in names:
        res['inequalities'][name] = dict(
            label=INEQS[name]['label'], norm=INEQS[name]['norm'], ref=INEQS[name]['ref'], kind=INEQS[name]['kind'],
            order=INEQS[name]['order'], sharp=INEQS[name]['exact'] is not None,
            dims={str(r['job'][2]): r for r in out if r['job'][0] == 'env' and r['job'][1] == name},
            sectors=next(r for r in out if r['job'][0] == 'sec' and r['job'][1] == name),
            scipy=next((r for r in out if r['job'][0] == 'sci' and r['job'][1] == name), None))
    res['seconds'] = time.time() - t0
    res['summary'] = summarise(res)
    log(f"total {res['seconds']:.0f} s ({workers} worker{'s' if workers > 1 else ''})")
    return res


def _log_task(r, log):
    if r['job'][0] == 'sci':
        log(f"{r['job'][1]:12s} SciPy reference on M_2: {r['values']} ({r['seconds']:.0f} s)")
        return
    if r['job'][0] == 'env':
        _, name, n = r['job']
        ref = r['reference']
        msg = ' '.join(f"{v:.3g}" for v in r['certified'])
        tail = '' if ref is None else '; ref ' + ' '.join(f"{v:.3g}" for v in ref)
        log(f"{name:12s} n={n}: certified V = [{msg}]{tail} ({r['seconds']:.0f} + {r['seconds_certify']:.0f} s)")
    else:
        cen = [round(v, 10) for v, m in zip(r['certified'], r['meta']) if m['type'] == 'central']
        intra = [round(v, 6) for v, m in zip(r['certified'], r['meta']) if m['type'] == 'intra']
        log(f"{r['job'][1]:12s} M_2+M_3: central W: V = {cen}; intra-block contrast in block 1: V = {intra} "
            f"({r['seconds']:.0f} s)")


def _at(grid, vals, d):
    i = int(np.argmin(np.abs(np.asarray(grid) - d)))
    return None if abs(grid[i] - d) > 1e-12 else vals[i]


def summarise(res):
    """Numbers for the paper: per inequality and overall."""
    grid = np.array(res['grid'])
    dims = sorted(int(n) for n in next(iter(res['inequalities'].values()))['dims'])
    S = dict(per={}, cert_rel_dev=0.0)
    v0 = -np.inf
    sharp_dev, bound_excess, ndev = 0.0, -np.inf, 0.0
    sec_central, sec_excess, sec_det, reached, intra_total = 0.0, -np.inf, 0.0, 0, 0
    sgm_best = dict(ratio=0.0)
    for name, r in res['inequalities'].items():
        p = dict()
        V2 = np.array(r['dims']['2']['certified'])
        for n in dims:
            d = r['dims'][str(n)]
            V, Vc, ref = np.array(d['V']), np.array(d['certified']), d['reference']
            if (grid == 0).any():
                v0 = max(v0, float(Vc[grid == 0].max()))
            big = V > 1e-8
            if big.any():
                S['cert_rel_dev'] = max(S['cert_rel_dev'], float(np.max(np.abs(Vc[big] - V[big]) / V[big])))
            if ref is not None:
                ref = np.array(ref)
                m = (grid >= 0.05) & (ref > 0)
                rel = (Vc[m] - ref[m]) / ref[m]
                p[f'rel_{n}'] = rel.tolist()
                if not m.any():
                    pass
                elif r['sharp'] and not (name == 'projections' and n != 2):
                    sharp_dev = max(sharp_dev, float(np.max(np.abs(rel))))
                elif name == 'sgm' and n >= 3:
                    k = int(np.argmax(rel))
                    if rel[k] + 1 > sgm_best['ratio']:
                        sgm_best = dict(ratio=float(rel[k] + 1), n=n, delta=float(grid[m][k]), V=float(Vc[m][k]),
                                        kappa_minus_1=float(ref[m][k]), point=d['points'][int(np.nonzero(m)[0][k])])
                else:
                    bound_excess = max(bound_excess, float(np.max(rel)))
            if n != 2:
                m = (grid >= 0.1) & (V2 > 1e-4)
                if name in N_INDEPENDENT and m.any():
                    ndev = max(ndev, float(np.max(np.abs(Vc[m] / V2[m] - 1))))
                p[f'ratio_{n}'] = (Vc[m] / V2[m]).tolist()
        v01, v02 = _at(grid, V2, 0.1), _at(grid, V2, 0.2)
        p['slope'] = float(np.log(v02 / v01) / np.log(2)) if v01 and v02 and v01 > 0 and v02 > 0 else None
        ref2 = reference_curve(name, np.array([0.1, 0.5, 0.9]), 2)
        p['V'] = {str(x): _at(grid, V2, x) for x in (0.1, 0.5, 0.9)}
        p['ref'] = None if ref2 is None else ref2.tolist()
        p['structure'] = {str(n): r['dims'][str(n)]['structure'][int(np.argmin(np.abs(grid - 0.5)))] for n in dims}
        sec = r['sectors']
        cen = [v for v, m in zip(sec['certified'], sec['meta']) if m['type'] == 'central']
        if name == 'det':
            exact = det_central(sec['meta'])
            sec_det = max(sec_det, max(abs(v - e) / e for v, e in zip(sec['certified'], exact) if e > 1e-12))
            p['sector_det_values'] = sec['certified']
        else:
            sec_central = max(sec_central, max(cen))
            for v, m in zip(sec['certified'], sec['meta']):
                if m['type'] == 'intra':
                    rc = reference_curve(name, np.array([m['delta_intra']]), 2)
                    ref1 = float(rc[0]) if rc is not None else _at(grid, V2, m['delta_intra'])
                    if not ref1 or ref1 <= 0:
                        continue
                    if rc is not None:
                        sec_excess = max(sec_excess, v / ref1 - 1)
                    if name in ('gardner', 'triangle', 'realpart', 'square', 'metric', 'young', 'projections',
                                'monotone'):
                        intra_total += 1
                        reached += int(v >= (1 - 1e-3) * ref1)
        S['per'][name] = p
    sci_ratios = []
    for name, r in res['inequalities'].items():
        if r.get('scipy'):
            V2 = np.array(r['dims']['2']['certified'])
            for d, v in r['scipy']['values'].items():
                lib = _at(grid, V2, float(d))
                if lib is not None and v > 1e-9:
                    sci_ratios.append((name, float(d), lib / v))
    S['scipy_ratios'] = sci_ratios
    S.update(V0=v0, sharp_dev=sharp_dev, bound_excess=bound_excess, n_dev=ndev, sector_central=sec_central,
             sector_excess=sec_excess, sector_reached=reached, sector_total=intra_total, sector_det=sec_det,
             sgm_best=sgm_best)
    return S


# ----------------------------------------------------------------------------------------------------------------
# output
# ----------------------------------------------------------------------------------------------------------------
def write_outputs(args, res):
    from common import save_json, write_tex, env_macro, sci
    cfg, S = res['config'], res['summary']
    save_json(args, 'trace_characterizations', res)

    def num(x):
        if x is None:
            return '--'
        if x <= 1e-12:
            return 'n.r.'                     # no violation resolved by the search at this contrast
        if abs(x) < 1e-3:
            return f"${sci(x, 2)}$"
        return f"{x:.4f}" if x < 10 else f"{x:.2f}"

    def sc(x, floor=1e-16):
        return sci(max(x, floor), 1)

    sg = S['sgm_best']
    macros = [('Starts', cfg['starts']), ('Steps', cfg['rounds'] * cfg['steps']), ('GridSize', len(cfg['grid'])),
              ('Digits', cfg['digits']), ('Dims', ', '.join(map(str, cfg['dims']))),
              ('SectorStarts', cfg['sector_starts']), ('Workers', res['workers']), ('Seconds', f"{res['seconds']:.0f}"),
              ('Vzero', sc(S['V0'])), ('CertDev', sc(S['cert_rel_dev'])), ('SharpDev', sc(S['sharp_dev'])),
              ('BoundExcess', sc(S['bound_excess'])), ('NDev', sc(S['n_dev'])),
              ('SectorCentral', sc(S['sector_central'])), ('SectorExcess', sc(S['sector_excess'])),
              ('SectorReached', S['sector_reached']), ('SectorTotal', S['sector_total']),
              ('SectorDet', sc(S['sector_det'])),
              ('SciRatioMin', f"{min(x[2] for x in S['scipy_ratios']):.3f}" if S['scipy_ratios'] else '--'),
              ('SciRatioMax', f"{max(x[2] for x in S['scipy_ratios']):.3f}" if S['scipy_ratios'] else '--'),
              ('SciStarts', cfg.get('scipy_starts', '--')),
              ('SGMExcessRatio', f"{sg['ratio']:.2f}" if sg.get('n') else '--'),
              ('SGMExcessV', f"{sg['V']:.3f}" if sg.get('n') else '--'),
              ('SGMExcessKappa', f"{sg['kappa_minus_1']:.3f}" if sg.get('n') else '--'),
              ('SGMExcessDelta', f"{sg['delta']:g}" if sg.get('n') else '--'),
              ('SGMExcessN', sg.get('n', '--'))]
    # the search against the constructions that are only approached (items 6, 7), and against 9, 10
    for name, tag in (('sgm', 'SGM'), ('alt', 'ALT'), ('monotone', 'Mon'), ('projections', 'Proj')):
        if name in S['per']:
            p = S['per'][name]
            for x, xt in ((0.5, 'Half'), (0.9, 'Nine'), (0.1, 'Tenth')):
                v, rf = p['V'][str(x)], p['ref'][[0.1, 0.5, 0.9].index(x)] if p['ref'] else None
                macros.append((f'{tag}Ratio{xt}', '--' if not v or not rf or v <= 0 else f'{v / rf:.3f}'))
    tex = env_macro(args, 'TC')
    for k, v in macros:
        tex += "\\newcommand{\\TC" + k + "}{" + str(v) + "}\n"
    rows = ''
    for i, (name, r) in enumerate(res['inequalities'].items(), 1):
        p = S['per'][name]
        ref = p['ref']
        star = '' if (r['sharp'] and name != 'projections') or not ref else \
            ('$^{a}$' if name == 'projections' else '$^{b}$')
        slope = '--' if p['slope'] is None else f"{p['slope']:.2f}"
        ratios = [abs(v - 1) for n in (3, 4) for v in p.get('ratio_%d' % n, [])]
        ndep = 'yes' if name == 'det' else ('$n\\ge3$' if name == 'sgm' else
                                            (f"${sci(max(ratios), 1)}$" if ratios else '--'))
        cites = ','.join(INEQS[name]['ref'].split(','))
        order = '--' if r['order'] is None else f"$\\delta^{{{r['order'][0]}}}$"
        rows += (f"{i} & {r['label']} & \\cite{{{cites}}} & {num(p['V']['0.1'])} & {num(p['V']['0.5'])} & "
                 f"{num(p['V']['0.9'])} & {num(ref[1]) if ref else '--'}{star} & "
                 f"{order} & {slope} & {ndep} \\\\\n")
    tex += "\\newcommand{\\TCRows}{%\n" + rows + "}\n"
    write_tex(args, 'trace_characterizations', tex)
    make_figure(args, res)


def make_figure(args, res):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    grid = np.array(res['grid'])
    fam = [('$\\kappa-1$: 1--5 (sharp), 6 (lower bound)', ['gardner', 'triangle', 'realpart', 'square', 'metric', 'sgm'],
            lambda d: kantorovich(d), 'C0'),
           ('$\\kappa^2-1$: 7 (lower bound)', ['alt'], alt_bound, 'C1'),
           ('$2\\delta/(1-\\delta)$: 8 (sharp)', ['young'], young_exact, 'C2'),
           ('9 ($n=2$, sharp)', ['projections'], projections_exact, 'C3'),
           ('10 (lower bound)', ['monotone'], monotone_bound, 'C4'),
           ('$1-\\sqrt{1-\\delta^2}$: 11 ($n=2$, sharp)', ['det'], lambda d: det_exact(d, 2), 'C5')]
    markers = dict(gardner='o', triangle='s', realpart='^', square='v', metric='D', sgm='x', alt='o', young='o',
                   projections='o', monotone='o', det='o')
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.4))
    ax = axes[0]
    dd = np.geomspace(0.01, 0.95, 200)
    for label, names, f, col in fam:
        ax.plot(dd, f(dd), color=col, lw=1.1, label=label)
        for name in names:
            if name not in res['inequalities']:
                continue
            V = np.array(res['inequalities'][name]['dims']['2']['certified'])
            m = (grid > 0) & (V > 0)
            ax.plot(grid[m], V[m], markers[name], ms=3.2, mfc='none', color=col, mew=0.8)
    extra = [('subadd', 's', '0.15'), ('powers', 'D', '0.35'), ('holder', '^', '0.5'), ('gt', 'v', '0.65'),
             ('pb', 'P', '0.8')]
    first = INEQS_INDEX
    for name, mk, col in extra:
        if name not in res['inequalities']:
            continue
        V = np.array(res['inequalities'][name]['dims']['2']['certified'])
        m = (grid > 0) & (V > 0)
        ax.plot(grid[m], V[m], ':', marker=mk, ms=3.0, mfc='none', color=col, lw=0.7, mew=0.7,
                label=f"{first[name]} (search)")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(1e-10, 50)
    ax.set_xlabel('$\\delta=\\Delta(W)$')
    ax.set_ylabel('$V(\\delta)$')
    ax.set_title('$M_2$: search (symbols) and closed forms', fontsize=9)
    hs, ls = ax.get_legend_handles_labels()
    nf = len(fam)
    fam_handles = (hs[:nf], ls[:nf])
    if len(hs) > nf:
        ax.legend(hs[nf:], ls[nf:], frameon=False, fontsize=5.6, loc='lower right')
    ax = axes[1]
    dims = sorted(int(n) for n in next(iter(res['inequalities'].values()))['dims'])
    styles = {2: '-', 3: '--', 4: ':'}
    for label, names, f, col in fam:
        for name in names:
            if name not in res['inequalities']:
                continue
            for n in dims:
                d = res['inequalities'][name]['dims'][str(n)]
                ref = reference_curve(name, grid, n)
                V = np.array(d['certified'])
                m = (grid >= 0.05) & (ref > 0) & (V > 0)
                ax.plot(grid[m], V[m] / ref[m], styles[n], color=col, lw=0.9, marker=markers[name], ms=2.5,
                        mfc='none', mew=0.6)
    for n in dims:
        ax.plot([], [], styles[n], color='0.3', lw=0.9, label=f'$n={n}$')
    ax.axhline(1, color='0.6', lw=0.6)
    ax.set_ylim(0, 1.45)
    ax.set_xlabel('$\\delta=\\Delta(W)$')
    ax.set_ylabel('$V/V^\\star$')
    ax.set_title('$M_n$: search / closed form or bound $V^\\star$', fontsize=9)
    leg = ax.legend(frameon=False, fontsize=6.5, loc='upper right', ncol=3)
    ax.add_artist(leg)
    ax.legend(*fam_handles, frameon=False, fontsize=5.6, loc='lower right')
    fig.tight_layout()
    fig.savefig(f"{args.figdir}/trace_characterizations.pdf", metadata={"CreationDate": None})
    plt.close(fig)


def main():
    from common import parse_args
    args = parse_args(__doc__)
    res = run(args.mode, args.device, args.seed, log=lambda s: print(s, flush=True))
    write_outputs(args, res)
    print("written:", args.out, args.figdir)


if __name__ == '__main__':
    main()
