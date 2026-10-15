"""
Numerical analysis behind the library: the Daleckii-Krein formula and the floating-point evaluation of divided
differences (Theorem 1.11, Propositions 1.12 and 1.13 of the lecture notes; functional.divided_differences), and the
error of the Krylov propagator with the step rule built on it (Lemma 9.17, Theorem 9.18, Corollary 9.19;
krylov.expm_multiply_step, krylov.evolve).

Independent references: scipy.linalg.expm_frechet and dense eigendecompositions in NumPy, scipy.special.jv for the
Bessel coefficients, and exact divided differences in high precision (mpmath, skipped if missing) or from expm1.
"""
import math
import warnings

import numpy as np
import pytest
import scipy.linalg as sla
import scipy.special as sps
import torch

from torch_vn_algebra import SpinChain, TypeIAlgebra, cost, krylov
from torch_vn_algebra import functional as fc

f64, c128 = torch.float64, torch.complex128
EPS = float(torch.finfo(f64).eps)
U_ROUND = EPS / 2                                  # unit roundoff of IEEE double


def herm(n, seed, spectrum=None, complex_valued=True):
    """A Hermitian matrix with a prescribed spectrum (or a GUE-like one) and Haar-like eigenvectors (NumPy)."""
    rng = np.random.default_rng(seed)
    G = rng.standard_normal((n, n)) + (1j * rng.standard_normal((n, n)) if complex_valued else 0)
    Q, _ = np.linalg.qr(G)
    lam = rng.standard_normal(n) if spectrum is None else np.asarray(spectrum, dtype=float)
    A = (Q * lam) @ Q.conj().T
    return (A + A.conj().T) / 2


# ==========================================================================================================
# Chapter 1: Daleckii-Krein and divided differences
# ==========================================================================================================
def loewner_derivative(A, E, f, df):
    """U (L o U^* E U) U^* with the library's Loewner matrix (Theorem 1.11)."""
    w, U = np.linalg.eigh(A)
    wt = torch.tensor(w)
    L = fc.divided_differences(wt, f(wt), df(wt)).numpy()
    return U @ (L * (U.conj().T @ E @ U)) @ U.conj().T


SPECTRA = {
    'generic': None,
    'degenerate': [1.0, 1.0, 1.0, -0.5, 2.0],
    'near_degenerate': [1.0, 1.0 + 1e-9, 1.0 - 3e-13, -0.5, 2.0],
    'zero_as_pm_1e-17': [1e-17, -1e-17, 0.4, 1.0, -2.0],
    'shifted_1e3': [1e3 - 0.7, 1e3, 1e3 + 1e-7, 1e3 + 0.4, 1e3 + 1.1],
}


@pytest.mark.parametrize('name', SPECTRA)
def test_theorem_1_11_daleckii_krein_equals_expm_frechet(name):
    """Theorem 1.11: d/dt e^{A+tE}|_0 = U (L o U^* E U) U^* (also at degenerate and shifted spectra) = expm_frechet."""
    spec = SPECTRA[name]
    A = herm(5, 1, spec)
    shift = 0.0 if spec is None else float(np.mean(spec))
    for seed in range(3):
        E = herm(5, 10 + seed)
        # e^{A} = e^{c} e^{A - c}: compare in the scale of e^{A - c} (e^{1000} overflows)
        D = loewner_derivative(A - shift * np.eye(5), E, torch.exp, torch.exp)
        _, ref = sla.expm_frechet(A - shift * np.eye(5), E)
        assert np.abs(D - ref).max() < 1e-12 * max(1.0, np.abs(ref).max())


def test_theorem_1_11_exact_identity_for_finite_t():
    """Proof of Theorem 1.11: U^*(f(A+tE) - f(A)) V = t (L_t o U^* E V) with (L_t)_ij = f[lambda_i, mu_j], exactly."""
    A, E = herm(6, 2, [0.0, 0.0, 1.0, 1.0, 2.5, -1.0]), herm(6, 3)
    t = 0.3
    lam, U = np.linalg.eigh(A)
    mu, V = np.linalg.eigh(A + t * E)
    f = np.sin
    fA, fAt = (U * f(lam)) @ U.conj().T, (V * f(mu)) @ V.conj().T
    Lt = np.where(np.abs(lam[:, None] - mu[None, :]) > 0,
                  (f(lam[:, None]) - f(mu[None, :])) / (lam[:, None] - mu[None, :] + (lam[:, None] == mu[None, :])),
                  np.cos(lam[:, None]))
    lhs = U.conj().T @ (fAt - fA) @ V
    rhs = t * Lt * (U.conj().T @ E @ V)
    assert np.abs(lhs - rhs).max() < 1e-13


@pytest.mark.parametrize('name', SPECTRA)
def test_theorem_1_11_gradient_of_apply_function_matches_expm_frechet(name):
    """Theorem 1.11 in autograd: the gradient of Re Tr(G^* e^A) is the adjoint Frechet derivative (expm_frechet)."""
    spec = SPECTRA[name]
    A = herm(5, 4, spec) - (0.0 if spec is None else float(np.mean(spec))) * np.eye(5)
    alg = TypeIAlgebra([5], [5], precision='double')
    X = torch.tensor(A, dtype=c128)[None, None].requires_grad_(True)
    G = torch.tensor(herm(5, 5) + 1j * herm(5, 6), dtype=c128)
    F = alg.operator(X, is_self_adjoint=True).apply_function(torch.exp).matrix[0, 0]
    (G.conj() * F).sum().real.backward()
    for seed in range(3):
        E = herm(5, 20 + seed)
        _, Lf = sla.expm_frechet(A, E)
        ref = (G.numpy().conj() * Lf).sum().real
        got = (X.grad[0, 0].numpy().conj() * E).sum().real
        assert abs(got - ref) < 1e-11 * max(1.0, abs(ref))


def test_proposition_1_12_midpoint_identity_and_bound():
    """Proposition 1.12: f[a,b] - (f'(a)+f'(b))/2 = -(h^2/2) int_0^1 s(1-s) f'''(a+sh) ds, |.| <= h^2 max|f'''| / 12."""
    mpmath = pytest.importorskip('mpmath')
    mpmath.mp.dps = 40
    cases = [(mpmath.exp, lambda x: mpmath.exp(x)), (mpmath.sin, lambda x: -mpmath.cos(x)),
             (mpmath.log, lambda x: 2 / x ** 3), (lambda x: x ** 3, lambda x: 6 + 0 * x)]
    for f, d3 in cases:
        for a, h in ((0.3, 0.5), (1.0, 1e-2), (2.0, 3.0)):
            a, h = mpmath.mpf(a), mpmath.mpf(h)
            b = a + h
            err = (f(b) - f(a)) / h - (mpmath.diff(f, a) + mpmath.diff(f, b)) / 2
            integral = -h ** 2 / 2 * mpmath.quad(lambda s: s * (1 - s) * d3(a + s * h), [0, 1])
            assert abs(err - integral) < mpmath.mpf(10) ** -30
            M3 = max(abs(d3(a + h * k / 200)) for k in range(201))
            assert abs(err) <= h ** 2 * M3 / 12 * (1 + 1e-12)
    # exact for cubics: -h^2 f''' / 12
    assert abs(float(err) + float(h) ** 2 * 6 / 12) < 1e-12


def test_proposition_1_13a_rounding_bound_of_the_quotient():
    """Proposition 1.13(a): |fl quotient - f[a,b]| <= (delta_a + delta_b)/|a-b| + 5u|q| for exactly rounded f values."""
    mpmath = pytest.importorskip('mpmath')
    mpmath.mp.dps = 50
    rng = np.random.default_rng(7)
    worst = 0.0
    for f in (mpmath.exp, mpmath.log, mpmath.sqrt, mpmath.atan):
        for _ in range(300):
            a = float(rng.uniform(0.5, 3.0) * 10.0 ** rng.integers(-3, 3))
            b = float(a * (1 + 10.0 ** rng.uniform(-15, -1) * rng.choice([-1, 1])))
            if a == b:
                continue
            fa, fb = float(f(a)), float(f(b))                        # exactly rounded values
            da, db = abs(fa - f(a)), abs(fb - f(b))                  # their errors delta_a, delta_b
            assert da <= U_ROUND * abs(fa) * (1 + 1e-15) and db <= U_ROUND * abs(fb) * (1 + 1e-15)
            q = (np.float64(fa) - np.float64(fb)) / (np.float64(a) - np.float64(b))
            x = (f(mpmath.mpf(a)) - f(mpmath.mpf(b))) / (mpmath.mpf(a) - mpmath.mpf(b))
            R = (da + db) / abs(mpmath.mpf(a) - mpmath.mpf(b)) + 5 * U_ROUND * abs(q)
            assert abs(q - x) <= R * (1 + mpmath.mpf(10) ** -30)
            worst = max(worst, float(abs(q - x) / R))
    assert worst > 0.05                                              # the bound is not vacuous


# (name, f (torch), f' (torch), f (mpmath), max |f'''| on [a, b] (mpmath), base point, scales)
def _cases(mpmath):
    def m3_exp(c):
        return lambda a, b: mpmath.exp(max(a, b) - c)
    return [
        ('exp(x-c), c=0', lambda x: torch.exp(x), lambda x: torch.exp(x), mpmath.exp, m3_exp(0), 0.3),
        ('exp(x-c), c=1e4', lambda x: torch.exp(x - 1e4), lambda x: torch.exp(x - 1e4),
         lambda x: mpmath.exp(x - 10000), m3_exp(10000), 1e4 + 0.3),
        ('exp(-7.3ix), x~1e2', lambda x: torch.exp(-7.3j * x), lambda x: -7.3j * torch.exp(-7.3j * x),
         lambda x: mpmath.exp(-7.3j * x), lambda a, b: mpmath.mpf(7.3) ** 3, 100.3),
        ('sin, x~1e3', torch.sin, torch.cos, mpmath.sin, lambda a, b: 1, 1e3 + 0.3),
        ('log, x~1e-8', torch.log, lambda x: 1 / x, mpmath.log, lambda a, b: 2 / min(a, b) ** 3, 1.3e-8),
        ('log, x~1e8', torch.log, lambda x: 1 / x, mpmath.log, lambda a, b: 2 / min(a, b) ** 3, 1.3e8),
        ('sqrt, x~1e-8', torch.sqrt, lambda x: 0.5 / torch.sqrt(x), mpmath.sqrt,
         lambda a, b: mpmath.mpf(3) / 8 * min(a, b) ** -2.5, 1.3e-8),
        ('x^-0.7, x~1e8', lambda x: x ** -0.7, lambda x: -0.7 * x ** -1.7, lambda x: x ** mpmath.mpf(-0.7),
         lambda a, b: abs(mpmath.mpf(-0.7) * -1.7 * -2.7) * min(a, b) ** -3.7, 1.3e8),
        ('x^3, x~1', lambda x: x ** 3, lambda x: 3 * x ** 2, lambda x: x ** 3, lambda a, b: 6, 1.3),
    ]


def test_proposition_1_13_library_choice_meets_the_bound_across_scales():
    """Proposition 1.13(b): the library's Loewner entry has error <= min(T, 3R) (T: Proposition 1.12 plus rounding of
    the midpoint, R: part (a) with the library's evaluation model), for gaps 1e-15..1e-1 relative, scales 1e-8..1e8."""
    mpmath = pytest.importorskip('mpmath')
    mpmath.mp.dps = 50
    for name, f, df, fm, m3, a0 in _cases(mpmath):
        for k in range(-60, -3):
            h = 10.0 ** (k / 4) * (abs(a0) if 'log' in name or 'sqrt' in name or '^' in name else 1.0)
            w = torch.tensor([a0, a0 + h], dtype=f64)
            a, b = (mpmath.mpf(v) for v in w.tolist())
            if a == b:
                continue
            fw, dfw = f(w), df(w)
            # the library's evaluation model must hold for the inputs we feed (Proposition 1.13: delta_a)
            d_model = [EPS * (abs(complex(fw[i])) + abs(w[i].item() * complex(dfw[i]))) for i in range(2)]
            d_true = [abs(complex(fw[i]) - complex(fm(v))) for i, v in enumerate((a, b))]
            assert all(t <= m for t, m in zip(d_true, d_model)), name
            got = complex(fc.divided_differences(w, fw, dfw)[0, 1])
            x = complex((fm(a) - fm(b)) / (a - b))
            q = complex((fw[0] - fw[1]) / (w[0] - w[1]))
            mid = complex((dfw[0] + dfw[1]) / 2)
            dd = [abs(complex(dfw[i]) - complex(mpmath.diff(fm, v))) for i, v in enumerate((a, b))]
            T = float(abs(b - a) ** 2 * m3(a, b) / 12) + (dd[0] + dd[1]) / 2 + 2 * U_ROUND * abs(mid)
            R = (d_model[0] + d_model[1]) / abs(float(b - a)) + 4 * EPS * abs(q)
            assert abs(got - x) <= min(T, 3 * R) * (1 + 1e-9) + 1e-300, (name, h, abs(got - x), T, R)


def test_proposition_1_13_counterexample_to_the_old_threshold():
    """Proposition 1.13 (regression): the former rule |w_i - w_j| <= eps^(1/3) max|w| took the midpoint for
    f(x) = e^{x - c} at c = 1e4 and gap 0.05 (error h^2/12 ~ 2e-4 relative); the quotient is accurate there."""
    for c in (1e2, 1e4, 1e6):
        h = 0.8 * EPS ** (1 / 3) * c
        w = torch.tensor([c, c + h], dtype=f64)
        f = lambda x: torch.exp(x - c)
        exact = math.expm1(float(w[1] - w[0])) / float(w[1] - w[0])          # e^{a-c} = 1
        midpoint_error = abs((1 + math.exp(float(w[1] - w[0]))) / 2 - exact)
        got = fc.divided_differences(w, f(w), f(w))[0, 1].item()
        assert abs(got - exact) < 1e-9 * exact < midpoint_error


def test_divided_differences_unresolved_pairs_and_ties():
    """Proposition 1.13 and the unresolved-pair rule: ties and pairs below k eps max|w| get the midpoint (also for a
    discontinuous f), separated pairs the quotient; in float32 the error stays within the same bound."""
    w = torch.tensor([-1e-17, 1e-17, 0.5, 1.0], dtype=f64)
    ind = lambda x: (x > 0).to(f64)
    L = fc.divided_differences(w, ind(w), torch.zeros_like(w))
    assert L[0, 1] == 0 and L[1, 2] == 0 and abs(L[0, 2].item() - 1 / (0.5 + 1e-17)) < 1e-15
    wt = torch.tensor([2.0, 2.0, 3.0], dtype=f64)
    assert torch.allclose(fc.divided_differences(wt, torch.exp(wt), torch.exp(wt))[0, 1], torch.exp(wt[0]))
    # float32: the quotient error eps/h and the midpoint error h^2/12 balance near h ~ eps^(1/3)
    eps32 = torch.finfo(torch.float32).eps
    worst = 0.0
    for k in range(-28, -2):
        h = 10.0 ** (k / 4)
        w32 = torch.tensor([0.3, 0.3 + h], dtype=torch.float32)
        hh = float(w32[1].double() - w32[0].double())
        if hh == 0:
            continue
        exact = math.exp(float(w32[0])) * math.expm1(hh) / hh
        worst = max(worst, abs(fc.divided_differences(w32, torch.exp(w32), torch.exp(w32))[0, 1].item() - exact))
    assert worst < 3 * (3 * 4 * eps32) ** (2 / 3) * (math.exp(0.31) / 12) ** (1 / 3) + 10 * eps32


# ==========================================================================================================
# Chapter 9: Krylov propagator
# ==========================================================================================================
class DenseOp:
    """A batch of dense Hermitian matrices (batch, k, k) seen only through matvec, as krylov expects."""

    def __init__(self, M):
        self.M = torch.as_tensor(M)
        self.dtype, self.device = self.M.dtype, self.M.device

    @property
    def batch_size(self):
        return self.M.shape[0]

    @property
    def dim(self):
        return self.M.shape[-1]

    def matvec(self, V):
        V = V.to(torch.promote_types(V.dtype, self.dtype))
        return (self.M.to(V.dtype) @ V.unsqueeze(-1)).squeeze(-1)


def adversarial(n=60, seed=0):
    """Spectra that are hard for the half-width estimate: two tight clusters with outliers, a large shift, a huge
    norm, a nearly uniform interval, a single spectral gap."""
    rng = np.random.default_rng(seed)
    return {
        'clusters+outliers': np.concatenate([1e-6 * rng.standard_normal(n // 2 - 1) - 1,
                                             1e-6 * rng.standard_normal(n // 2 - 1) + 1, [-7.0, 9.0]]),
        'shift_1e6': 1e6 + rng.uniform(-2, 2, n),
        'norm_1e8': 1e8 * rng.uniform(-1, 1, n),
        'interval': np.linspace(-3, 3, n),
        'gap': np.concatenate([rng.uniform(-1, -0.9, n // 2), rng.uniform(0.9, 1, n - n // 2)]),
    }


def test_lemma_9_17_chebyshev_coefficients_and_bound():
    """Lemma 9.17: e^{-i rho y} = c_0 + 2 sum c_k T_k(y), c_k = (-i)^k J_k(rho), |c_k| <= (|rho|/2)^k / k!, and the
    truncation after degree m-1 errs by at most 2 S_m(|rho|/2)."""
    y = np.linspace(-1, 1, 2001)
    for rho in (0.1, 1.0, 7.5, 30.0, -12.0):
        c = np.array([(-1j) ** k * sps.jv(k, rho) for k in range(120)])
        series = np.array([(-1j) ** k * sum((-1) ** l * math.exp((k + 2 * l) * math.log(abs(rho) / 2)
                                                                  - math.lgamma(l + 1) - math.lgamma(k + l + 1))
                                            * (np.sign(rho) ** k) for l in range(80)) for k in range(12)])
        assert np.abs(c[:12] - series).max() < 1e-12 * math.exp(abs(rho))
        bound = np.array([math.exp(k * math.log(abs(rho) / 2) - math.lgamma(k + 1)) for k in range(120)])
        assert (np.abs(c) <= bound * (1 + 1e-10) + 1e-300).all()
        for m in (1, 2, 5, 20, 40):
            p = c[0] + 2 * sum(c[k] * np.cos(k * np.arccos(y)) for k in range(1, m))
            tail = 2 * math.exp(abs(rho) / 2) * float(sps.gammainc(m, abs(rho) / 2))   # 2 S_m(|rho|/2)
            assert np.abs(p - np.exp(-1j * rho * y)).max() <= tail * (1 + 1e-9) + 1e-13


def _lanczos_numpy(H, v, m):
    """Plain Lanczos with full re-orthogonalisation in NumPy (reference for the library's Q, T, beta_m)."""
    Q = np.zeros((m, len(v)), dtype=complex)
    a, b = np.zeros(m), np.zeros(m)
    q = v / np.linalg.norm(v)
    for j in range(m):
        Q[j] = q
        r = H @ q
        a[j] = np.vdot(q, r).real
        for _ in range(2):
            r = r - Q[:j + 1].T @ (Q[:j + 1].conj() @ r)
        b[j] = np.linalg.norm(r)
        q = r / b[j]
    return Q, a, b


@pytest.mark.parametrize('name', list(adversarial()))
def test_theorem_9_18_bound_holds_on_adversarial_spectra(name):
    """Theorem 9.18(c): ||e^{-i tau H} v - w_tau|| <= (4 beta_m / w_T) S_m(|tau| w_T / 2) (Ritz half-width w_T), and (d)
    the a priori bound 4 S_m(|tau| w / 2) with the true half-width w is larger."""
    lam = adversarial()[name]
    n = len(lam)
    H = herm(n, 3, lam, complex_valued=False)
    w_true = (lam.max() - lam.min()) / 2
    v = np.random.default_rng(1).standard_normal(n) + 0j
    v /= np.linalg.norm(v)
    Hop = DenseOp(torch.tensor(H)[None])
    worst = 0.0
    for m in (4, 10, 25):
        for rho in (2.0, 0.4 * m, 0.9 * m, 1.5 * m):
            tau = rho / w_true
            out, bound = krylov.expm_multiply_step(Hop, torch.tensor(v)[None], tau, m)
            ref = sla.expm(-1j * tau * (H - lam.mean() * np.eye(n))) @ v * np.exp(-1j * tau * lam.mean())
            err = np.linalg.norm(out[0].numpy() - ref)
            apriori = 4 * math.exp(rho / 2) * float(sps.gammainc(m, rho / 2))
            assert err <= bound.item() * (1 + 1e-6) + 1e-9 * max(1.0, abs(tau * lam.mean()) * 1e-6), (m, rho)
            assert bound.item() <= apriori * (1 + 1e-9)
            worst = max(worst, err / max(bound.item(), 1e-300) if err > 1e-10 else 0.0)
    assert worst > 1e-4                                       # the bound is attained within a few orders


def test_theorem_9_18_b_duhamel_representation_and_defect():
    """Theorem 9.18(b): e^{-i tau H} v - w_tau = -i int_0^tau e^{-i(tau-t)H} r_m phi(t) dt, and beta_m |phi(tau)| is the
    norm of the defect w' + i H w (the former per-step estimate: an error rate, not an error)."""
    n, m, tau = 40, 6, 1.7
    H = herm(n, 8)
    v = np.random.default_rng(2).standard_normal(n) + 0j
    v /= np.linalg.norm(v)
    Q, a, b = _lanczos_numpy(H, v, m)
    T = np.diag(a) + np.diag(b[:m - 1], 1) + np.diag(b[:m - 1], -1)
    r = H @ Q[-1] - a[-1] * Q[-1] - b[m - 2] * Q[-2]
    assert abs(np.linalg.norm(r) - b[-1]) < 1e-12
    wfun = lambda t: Q.T @ sla.expm(-1j * t * T)[:, 0]
    phi = lambda t: sla.expm(-1j * t * T)[-1, 0]
    ts, wts = np.polynomial.legendre.leggauss(60)
    ts, wts = tau * (ts + 1) / 2, wts * tau / 2
    integral = -1j * sum(wt * sla.expm(-1j * (tau - t) * H) @ r * phi(t) for t, wt in zip(ts, wts))
    err = sla.expm(-1j * tau * H) @ v - wfun(tau)
    assert np.linalg.norm(err - integral) < 1e-12
    assert np.linalg.norm(err) <= b[-1] * sum(wt * abs(phi(t)) for t, wt in zip(ts, wts)) + 1e-14
    dt = 1e-6
    defect = (wfun(tau + dt) - wfun(tau - dt)) / (2 * dt) + 1j * H @ wfun(tau)
    assert abs(np.linalg.norm(defect) - b[-1] * abs(phi(tau))) < 1e-7
    # for a short step the endpoint value overestimates the error by about m / (tau * scale)
    _, est_bound = krylov.expm_multiply_step(DenseOp(torch.tensor(H)[None]), torch.tensor(v)[None], 0.05, m)
    small_err = np.linalg.norm(sla.expm(-0.05j * H) @ v - Q.T @ sla.expm(-0.05j * T)[:, 0])
    assert small_err <= est_bound.item() and b[-1] * abs(sla.expm(-0.05j * T)[-1, 0]) > 50 * small_err


def test_theorem_9_18_a_shift_invariance():
    """Theorem 9.18(a): H -> H + c leaves the Lanczos vectors, beta and the error bound unchanged, shifts T_m by c and
    multiplies the Krylov state by e^{-i tau c}."""
    n, m, tau = 60, 12, 0.9
    H = herm(n, 11, complex_valued=False)
    v = torch.tensor(np.random.default_rng(3).standard_normal(n))[None].to(c128)
    a0, b0, Q0 = krylov.lanczos(DenseOp(torch.tensor(H)[None]).matvec, v, m)
    for c in (3.0, -250.0, 1e5):
        a1, b1, Q1 = krylov.lanczos(DenseOp(torch.tensor(H + c * np.eye(n))[None]).matvec, v, m)
        assert torch.allclose(a1, a0 + c, atol=1e-9 * max(1, abs(c))) and torch.allclose(b1, b0, atol=1e-8)
        assert (Q1 - Q0).abs().max() < 1e-8
        w0, e0 = krylov.expm_multiply_step(DenseOp(torch.tensor(H)[None]), v, tau, m)
        w1, e1 = krylov.expm_multiply_step(DenseOp(torch.tensor(H + c * np.eye(n))[None]), v, tau, m)
        assert (w1 - np.exp(-1j * tau * c) * w0).abs().max() < 1e-9
        assert abs(e1.item() - e0.item()) <= 1e-6 * e0.item() + 1e-300


def test_theorem_9_18_d_beta_and_ritz_width_below_half_width():
    """Theorem 9.18(d): beta_m <= w and w_T <= w for the half-width w of the spectrum of H."""
    for name, lam in adversarial(seed=5).items():
        n = len(lam)
        c, w = (lam.max() + lam.min()) / 2, (lam.max() - lam.min()) / 2
        H = herm(n, 6, lam - c, complex_valued=False)
        for m in (3, 8, 20):
            v = np.random.default_rng(m).standard_normal(n) + 0j
            Q, a, b = _lanczos_numpy(H, v, m)
            T = np.diag(a) + np.diag(b[:m - 1], 1) + np.diag(b[:m - 1], -1)
            theta = np.linalg.eigvalsh(T)
            assert b[-1] <= w * (1 + 1e-10) and (theta.max() - theta.min()) / 2 <= w * (1 + 1e-10), name


@pytest.mark.parametrize('name', list(adversarial(n=48)))
def test_corollary_9_19_step_rule_meets_tol(name):
    """Corollary 9.19: with the default step rule the error at every grid time t is <= tol (t - t0) / (t_end - t0)
    ||psi0|| (plus rounding), and no warning is issued; for a batch the common step serves every element."""
    lam = adversarial(n=48)[name]
    n = len(lam)
    Hs = [herm(n, 21, lam, complex_valued=False), herm(n, 22, lam[::-1] * 0.5 + lam.mean() * 0.5,
                                                       complex_valued=False)]
    Hop = DenseOp(torch.tensor(np.stack(Hs)))
    w_true = (lam.max() - lam.min()) / 2
    times = [0.0, 3.0 / w_true, 25.0 / w_true]
    psi0 = torch.randn(2, n, dtype=c128, generator=torch.Generator().manual_seed(4)) * 1.5
    for m, tol in ((8, 1e-4), (20, 1e-8), (30, 1e-11)):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always', cost.CostWarning)
            out = krylov.evolve(Hop, psi0, times, m=m, tol=tol)
        assert not [c for c in caught if 'error bound' in str(c.message)]
        for b in range(2):
            mu, V = np.linalg.eigh(Hs[b] - lam.mean() * np.eye(n))
            for t, o in zip(times, out):
                ref = V @ (np.exp(-1j * t * mu) * (V.T @ psi0[b].numpy())) * np.exp(-1j * t * lam.mean())
                allowed = tol * psi0[b].norm().item() * (t - times[0]) / (times[-1] - times[0])
                # rounding: the matrix-vector products of H = c + (H - c) carry errors eps |c| per unit time
                rounding = (1e-12 * max(1.0, t * w_true) + 20 * EPS * abs(t * lam.mean())) * psi0[b].norm().item()
                assert np.linalg.norm(o[b].numpy() - ref) <= allowed + rounding, (m, tol, t)


def test_step_rule_is_shift_invariant_and_not_wasteful():
    """Theorem 9.18(a) for the step rule: H and H + 1e5 take the same number of Lanczos runs; the steps are long
    (rho = tau w per step well above the former 0.1 m)."""
    ch = SpinChain(10, 'open', complex_valued=False, precision='double')
    psi0 = ch.vector_in_sector('1010101010')[None].to(c128)
    counts = []
    orig = krylov.lanczos
    for h in (0.0, 1e5):
        calls = []
        krylov.lanczos = lambda *a, **k: calls.append(1) or orig(*a, **k)
        try:
            H = ch.xxz_sparse(1.0, 0.5, h, sector=5)
            out = krylov.evolve(H, psi0, [0.0, 20.0], m=30, tol=1e-8)
        finally:
            krylov.lanczos = orig
        w, V = np.linalg.eigh(H.to_dense()[0].numpy())
        ref = V @ (np.exp(-20j * (w - w.mean())) * (V.T @ psi0[0].numpy())) * np.exp(-20j * w.mean())
        assert np.linalg.norm(out[-1, 0].numpy() - ref) < 1e-8 + 1e-15 * 20 * abs(w).max()
        counts.append(len(calls))
        width = (w.max() - w.min()) / 2
    assert counts[0] == counts[1] and counts[0] < 20.0 * width / (0.1 * 30) / 2


def test_evolve_fixed_dt_warns_from_the_bound_and_rejects_bad_tol():
    """evolve with a fixed dt warns when the sum of the bounds of Theorem 9.18(c) exceeds tol; tol <= 0 is refused."""
    ch = SpinChain(8, 'open', complex_valued=False, precision='double')
    H = ch.xxz_sparse(1.0, 0.5, 0.0, sector=4)
    psi0 = ch.vector_in_sector('10101010')[None].to(c128)
    with pytest.warns(cost.CostWarning):
        krylov.evolve(H, psi0, [0.0, 4.0], m=6, dt=2.0, tol=1e-8)
    with warnings.catch_warnings():
        warnings.simplefilter('error', cost.CostWarning)
        krylov.evolve(H, psi0, [0.0, 4.0], m=30, dt=0.5, tol=1e-8)
    with pytest.raises(ValueError):
        krylov.evolve(H, psi0, [0.0, 1.0], tol=0.0)


def test_krylov_error_bound_function_against_direct_sum():
    """krylov_error_bound = 2 beta |dt| S_m(x)/x = (4 beta / w_T) S_m(x), x = |dt| w_T / 2, against a direct sum."""
    ritz = torch.tensor([[-1.5, 0.2, 2.5]], dtype=f64)
    beta = torch.tensor([0.7], dtype=f64)
    for m in (2, 5, 30):
        for dt in (0.01, 0.5, 3.0, -3.0):
            wT, x = 2.0, abs(dt) * 2.0 / 2
            S = sum(math.exp(j * math.log(x) - math.lgamma(j + 1)) for j in range(m, m + 400))
            assert abs(krylov.krylov_error_bound(beta, ritz, dt, m).item() - 4 * 0.7 / wT * S) <= 1e-10 * 4 * 0.7 / wT * S


def test_krylov_scale_safe_norms():
    """Regression: krylov normalises with a scaled norm, so start vectors with entries 1e200 (double) or 1e20
    (single) give finite results equal to 1e200 (1e20) times those of the normalised vector."""
    for prec, big in (('double', 1e200), ('single', 1e20)):
        ch = SpinChain(6, precision=prec)
        H = ch.xxz_sparse(1.0, 0.5, None, sector=3)
        cdt = c128 if prec == 'double' else torch.complex64
        v = torch.randn(1, H.dim, dtype=cdt, generator=torch.Generator().manual_seed(9))
        v = v / v.abs().max()
        out_big = krylov.evolve(H, big * v, [0.0, 0.5])[-1]
        out = krylov.evolve(H, v, [0.0, 0.5])[-1]
        assert torch.isfinite(out_big).all()
        assert ((out_big * (1 / big) - out).abs().max() < (1e-12 if prec == 'double' else 1e-5))   # not / big: older
        # PyTorch divides complex by real without scaling
        E_big, _, _ = krylov.ground_state(H, v0=big * v)
        E, _, _ = krylov.ground_state(H, v0=v)
        assert torch.isfinite(E_big).all() and abs(E_big.item() - E.item()) < (1e-9 if prec == 'double' else 1e-4)
