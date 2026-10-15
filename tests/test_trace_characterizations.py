"""
Quantitative failure of trace characterisations (paper, Section "Use case: quantitative failure of trace
characterisations"; scripts/trace_characterizations.py). The closed forms and lower-bound constructions of the
Propositions, against random sampling, explicit extremisers and brute-force SciPy searches in plain NumPy; the
library objectives against the independent mpmath recomputation; the block (sector) structure. CPU, double
precision, fixed seeds.
"""
import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest
import torch
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location('trace_characterizations', ROOT / 'scripts' / 'trace_characterizations.py')
tc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tc)

pytestmark = [pytest.mark.filterwarnings("ignore::torch_vn_algebra.cost.CostWarning"),
              pytest.mark.filterwarnings("ignore::UserWarning")]


# ----------------------------------------------------------------------------------------------------------------
# plain NumPy helpers (independent of the library)
# ----------------------------------------------------------------------------------------------------------------
def density(delta, n=2, middle=None):
    """diag(1, l, l + t (1 - l) for t in middle) / trace, l = (1 - delta)/(1 + delta): Delta(W) = delta."""
    lo = (1 - delta) / (1 + delta)
    mid = [lo + t * (1 - lo) for t in (middle if middle is not None else [])][:n - 2]   # t in [0, 1]
    w = np.array([1.0, lo] + mid)
    return w / w.sum()


def hf(A, f):
    e, V = np.linalg.eigh((A + A.conj().T) / 2)
    return (V * f(e)) @ V.conj().T


def absm(X):
    return hf(X.conj().T @ X, lambda e: np.sqrt(np.clip(e, 0, None)))


def phi(w, X):
    return np.sum(w * np.diag(X))


def ginibre(rng, n):
    return rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))


def rand_pos(rng, n, rank=None):
    Z = ginibre(rng, n)[:, :rank or n]
    return Z @ Z.conj().T


def rand_proj(rng, n, r):
    Q, _ = np.linalg.qr(ginibre(rng, n))
    return Q[:, :r] @ Q[:, :r].conj().T


def sgm_np(A, B):
    ah, aih = hf(A, np.sqrt), hf(A, lambda e: e ** -0.5)
    g = aih @ hf(ah @ B @ ah, lambda e: np.sqrt(np.clip(e, 0, None))) @ aih
    gh = hf(g, lambda e: np.sqrt(np.clip(e, 0, None)))
    return gh @ A @ gh


def value_np(name, w, *ops):
    """(LHS - RHS)/N of the inequalities of the paper, plain NumPy."""
    if name == 'gardner':
        X, = ops
        return abs(phi(w, X)) / phi(w, absm(X)).real - 1
    if name == 'realpart':
        A, = ops
        return phi(w, A + A.conj().T).real / (2 * phi(w, absm(A.conj().T)).real) - 1
    if name == 'triangle':
        X, Y = ops
        return phi(w, absm(X + Y)).real / (phi(w, absm(X)) + phi(w, absm(Y))).real - 1
    if name == 'square':
        A, = ops
        return abs(phi(w, A @ A)) / phi(w, A.conj().T @ A).real - 1
    if name == 'metric':
        P, Q, R = ops
        return phi(w, absm(P - R)).real / (phi(w, absm(P - Q)) + phi(w, absm(Q - R))).real - 1
    if name == 'young':
        A, B = ops
        rhs = 0.5 * (phi(w, A.conj().T @ A) + phi(w, B.conj().T @ B)).real
        return phi(w, absm(A @ B.conj().T)).real / rhs - 1
    if name == 'projections':
        P, Q = ops
        return 2 * (phi(w, P @ Q).real - phi(w, Q @ P @ Q).real) / (phi(w, P) + phi(w, Q)).real
    if name == 'monotone':
        A, C = ops
        B = A + C
        return phi(w, A @ A).real / phi(w, B @ B).real - 1
    if name == 'alt':
        A, B = ops
        bh = hf(B, lambda e: np.sqrt(np.clip(e, 0, None)))
        M = bh @ A @ bh
        return phi(w, M @ M).real / phi(w, B @ A @ A @ B).real - 1
    if name == 'sgm':
        A, B = ops
        return phi(w, sgm_np(A, B)).real / math.sqrt(phi(w, A).real * phi(w, B).real) - 1
    if name == 'det':
        A, = ops
        n = len(w)
        ww = n * w / w.sum()
        return 1 - phi(ww, A).real / (n * np.linalg.det(A).real ** (1 / n))
    raise KeyError(name)


def random_args(name, rng, n):
    if name in ('gardner', 'realpart', 'square'):
        return (ginibre(rng, n),)
    if name == 'triangle':
        X, Y = ginibre(rng, n), ginibre(rng, n)
        return (X + X.conj().T, Y + Y.conj().T)
    if name == 'metric':
        return tuple(rand_proj(rng, n, int(rng.integers(0, n + 1))) for _ in range(3))
    if name == 'young':
        return (ginibre(rng, n), ginibre(rng, n))
    if name == 'projections':
        return tuple(rand_proj(rng, n, int(rng.integers(1, n + 1))) for _ in range(2))
    if name in ('monotone', 'alt'):
        return (rand_pos(rng, n, int(rng.integers(1, n + 1))), rand_pos(rng, n, int(rng.integers(1, n + 1))))
    if name in ('sgm', 'det'):
        return (rand_pos(rng, n) + 1e-3 * np.eye(n), rand_pos(rng, n) + 1e-3 * np.eye(n))[:1 if name == 'det' else 2]
    raise KeyError(name)


SHARP = ['gardner', 'triangle', 'realpart', 'square', 'metric', 'young', 'det']


# ----------------------------------------------------------------------------------------------------------------
# closed forms: upper bounds never exceeded, explicit extremisers attain them
# ----------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('name', tc.ORDER)
def test_random_points_never_exceed_the_proved_value(name):
    """No random point exceeds V* (the sharp value; for items 6, 7, 10 the lower bound that the search never exceeds
    for n = 2, and a/b - 1 for item 6 with n = 3), for n = 2, 3 and three-level W; V <= 0 for the trace (delta = 0)."""
    rng = np.random.default_rng(len(name))
    for n in (2, 3):
        for d in (0.0, 0.3, 0.8):
            w = density(d, n, middle=[0.6])
            vstar = tc.reference_curve(name, np.array([d]), n)
            vstar = 0.0 if vstar is None else float(vstar[0])
            if name == 'sgm' and n >= 3:
                vstar = float(tc.young_exact(d))              # only a/b - 1 is proved for n >= 3
            for _ in range(120):
                v = value_np(name, w, *random_args(name, rng, n))
                if np.isfinite(v):                            # 0/0 at P = Q = R in the metric
                    assert v <= vstar + 1e-7, (name, n, d, v, vstar)


def _kantorovich_point(d):
    """The extremisers of the Kantorovich class in M_2 for W = diag(a, b), a + b = 2 (Proposition)."""
    a, b = 1 + d, 1 - d
    w = np.array([a, b]) / 2
    v = np.array([math.sqrt(b / (a + b)), math.sqrt(a / (a + b))])          # v^*Wv = harmonic mean of a, b
    Wv = np.diag([a, b]) @ v
    u2 = Wv / np.linalg.norm(Wv)
    U = np.outer(u2, v) + np.outer(np.array([-u2[1], u2[0]]), np.array([-v[1], v[0]]))   # unitary, U v = Wv/|Wv|
    t = math.atan(math.sqrt(b / a))
    ue, ve = np.array([math.sin(t), math.cos(t)]), np.array([-math.sin(t), math.cos(t)])
    return w, U @ np.outer(v, v), np.outer(ue, ue), np.outer(ve, ve)


@pytest.mark.parametrize('d', [0.05, 0.4, 0.9])
def test_kantorovich_class_attained(d):
    """kappa - 1 is attained: Gardner and A + A^* at X = U P_v; the triangle inequality at X = P_u, Y = -P_v; the
    projection metric at P = P_u, Q = 0, R = P_v; |phi(A^2)| at an off-diagonal A with |A_12|/|A_21| = sqrt(a/b)."""
    w, X, Pu, Pv = _kantorovich_point(d)
    k = float(tc.kantorovich(d))
    # (|X| by the square root of X^*X: a zero eigenvalue costs ~1e-8 in double precision)
    assert abs(value_np('gardner', w, X) - k) < 1e-7
    assert abs(value_np('realpart', w, X.conj().T) - k) < 1e-7
    assert abs(value_np('triangle', w, Pu, -Pv) - k) < 1e-7
    assert abs(value_np('metric', w, Pu, np.zeros((2, 2)), Pv) - k) < 1e-7
    a, b = 1 + d, 1 - d
    A = np.array([[0, math.sqrt(a)], [math.sqrt(b), 0]], dtype=complex)
    assert abs(value_np('square', w, A) - k) < 1e-12


@pytest.mark.parametrize('n', [2, 3, 4])
def test_square_bound_in_every_dimension(n):
    """|Tr(W A^2)| <= kappa Tr(W A^* A) in M_n (the pairing proof), sharp at an A supported on the extreme pair."""
    rng = np.random.default_rng(n)
    for _ in range(200):
        w = rng.uniform(0.2, 1.0, n)
        w /= w.sum()
        d = (w.max() - w.min()) / (w.max() + w.min())
        assert value_np('square', w, ginibre(rng, n)) <= tc.kantorovich(d) + 1e-12


@pytest.mark.parametrize('d', [0.1, 0.5, 0.95])
def test_young_exact_and_attained(d):
    """Young (p = q = 2): V = a/b - 1 = 2 delta/(1 - delta), attained at the matrix unit A = B = E_12."""
    w = density(d)
    E = np.array([[0, 1], [0, 0]], dtype=complex)
    assert abs(value_np('young', w, E, E) - tc.young_exact(d)) < 1e-12


@pytest.mark.parametrize('n', [2, 3, 5])
def test_det_two_level_formula(n):
    """Determinant inequality: V = 1 - GM/AM of the eigenvalues of W at A = W^{-1}; the two-level formula is the
    maximum over spectra in [l, 1] (vertices of the cube, against random interior spectra); n = 2: 1 - sqrt(1 - d^2);
    large n: the Specht ratio."""
    rng = np.random.default_rng(n)
    for d in (0.2, 0.7):
        lo = (1 - d) / (1 + d)
        best = 0.0
        for _ in range(2000):
            w = np.concatenate([[1.0, lo], rng.uniform(lo, 1, n - 2)])
            best = max(best, value_np('det', w, np.diag(1 / w)))
        exact = float(tc.det_exact(d, n)[0])
        assert best <= exact + 1e-12
        m = max(range(1, n), key=lambda m: 1 - lo ** ((n - m) / n) / ((m + (n - m) * lo) / n))
        w = np.array([1.0] * m + [lo] * (n - m))
        assert abs(value_np('det', w, np.diag(1 / w)) - exact) < 1e-12
    assert np.allclose(tc.det_exact(np.array([0.3, 0.8]), 2), 1 - np.sqrt(1 - np.array([0.3, 0.8]) ** 2))
    assert abs(float(tc.det_exact(0.6, 4000)[0]) - float(tc.specht_limit(0.6)[0])) < 1e-6


def test_projections_n2_closed_form_by_scipy():
    """phi(PQ + QP) <= 2 phi(QPQ) on M_2: a BFGS search over rank-one P, Q (complex unit vectors) reaches
    s/(1 + s), s = sin(arcsin(delta)/2), and does not exceed it."""
    rng = np.random.default_rng(3)
    for d in (0.3, 0.8):
        w = density(d)

        def f(x):
            u = x[:2] + 1j * x[2:4]
            v = x[4:6] + 1j * x[6:]
            u, v = u / np.linalg.norm(u), v / np.linalg.norm(v)
            return -value_np('projections', w, np.outer(u, u.conj()), np.outer(v, v.conj()))
        best = max(-minimize(f, rng.normal(size=8), method='BFGS', options={'gtol': 1e-12}).fun for _ in range(20))
        assert abs(best - float(tc.projections_exact(d))) < 1e-8


def test_monotone_rank_one_construction():
    """Monotonicity (p = 2): the rank-one construction A = P_u, C = t P_v at angles a, pi - a to the eigenbasis of W
    (cos 2a = -y*) attains the lower bound m(delta) = max_y y^2(d - y)^2/((1 - 2dy + y^2)(1 - y^2)) ~ delta^4/16, and a
    BFGS search over complex A, C >= 0 in M_2 does not exceed it."""
    rng = np.random.default_rng(4)
    for d in (0.4, 0.8):
        w = density(d)
        m = float(tc.monotone_bound(d)[0])
        ys = np.linspace(1e-6, d - 1e-6, 200001)
        f = ys ** 2 * (d - ys) ** 2 / ((1 - 2 * d * ys + ys ** 2) * (1 - ys ** 2))
        y = ys[np.argmax(f)]
        a = 0.5 * math.acos(-y)
        u, v = np.array([math.cos(a), math.sin(a)]), np.array([math.cos(math.pi - a), math.sin(math.pi - a)])
        A, P = np.outer(u, u), np.outer(v, v)
        alpha = gamma = 1 - d * y                         # u^*Wu = v^*Wv for W = diag(1 + d, 1 - d)
        r = y * (d - y)                                   # -Re (u^*v)(v^*Wu)
        t = r / gamma
        assert abs(value_np('monotone', w, A.astype(complex), t * P.astype(complex)) - m) < 1e-7
        assert alpha > 0

        def g(x):
            La = (x[:4] + 1j * x[4:8]).reshape(2, 2)
            Lc = (x[8:12] + 1j * x[12:]).reshape(2, 2)
            return -value_np('monotone', w, La @ La.conj().T, Lc @ Lc.conj().T)
        best = max(-minimize(g, rng.normal(size=16), method='BFGS').fun for _ in range(30))
        assert best <= m + 1e-7
        assert best >= m - 1e-3 * m


def test_alt_limit_construction():
    """ALT (r = 2, q = 1): A = P_u, B = diag(eta, 1) in a basis at 45 degrees to the eigenbasis of W,
    u = f_1 - delta eta f_2: V -> kappa^2 - 1 = delta^2/(1 - delta^2) as eta -> 0 (mpmath, 60 digits)."""
    import mpmath as mp
    with mp.workdps(60):
        for d in (mp.mpf('0.3'), mp.mpf('0.7')):
            W = mp.diag([1 + d, 1 - d])
            R = mp.matrix([[1, -1], [1, 1]]) / mp.sqrt(2)                  # columns f_1, f_2
            Wf = R.T * W * R                                               # [[1, -d], [-d, 1]]
            target = d ** 2 / (1 - d ** 2)
            for eta in ('1e-6', '1e-12'):                   # the error decays like sqrt(eta)
                eta = mp.mpf(eta)
                u = mp.matrix([1, -Wf[0, 1] / Wf[1, 1] * eta])   # u = f_1 + s eta f_2, s = -W_12/W_22
                A = u * u.T
                B = mp.diag([eta, 1])
                bh = mp.diag([mp.sqrt(eta), 1])
                M = bh * A * bh
                tr = lambda X: X[0, 0] * Wf[0, 0] + X[1, 1] * Wf[1, 1] + 2 * X[0, 1] * Wf[0, 1]   # noqa: E731
                v = tr(M * M) / tr(B * A * A * B) - 1
                assert abs(v - target) < 5 * mp.sqrt(eta)


def test_sgm_limit_construction_and_identities():
    """Spectral geometric mean: B = P_v (v balanced between the eigenvectors of W), A = P_x + eps, x ~ W^{-1} v:
    V -> kappa - 1 as eps -> 0 (mpmath); A natural B has the eigenvalues sqrt(eig(AB)) and is symmetric in A, B."""
    import mpmath as mp
    rng = np.random.default_rng(6)
    for _ in range(5):
        A, B = rand_pos(rng, 3) + 0.1 * np.eye(3), rand_pos(rng, 3) + 0.1 * np.eye(3)
        S = sgm_np(A, B)
        assert np.allclose(np.sort(np.linalg.eigvalsh(S)), np.sort(np.sqrt(np.linalg.eigvals(A @ B).real)), atol=1e-9)
        assert np.allclose(S, sgm_np(B, A), atol=1e-9)
    with mp.workdps(60):
        def herm(X, f):
            E, Q = mp.eighe((X + X.transpose_conj()) / 2)
            return Q * mp.diag([f(e) for e in E]) * Q.transpose_conj()
        for d in (mp.mpf('0.2'), mp.mpf('0.6')):
            a, b = 1 + d, 1 - d
            v = mp.matrix([1, 1]) / mp.sqrt(2)
            x = mp.matrix([v[0] / a, v[1] / b])
            x = x / mp.norm(x)
            A, B = x * x.T + mp.mpf('1e-20') * mp.eye(2), v * v.T
            ah, aih = herm(A, mp.sqrt), herm(A, lambda e: 1 / mp.sqrt(e))
            g = aih * herm(ah * B * ah, lambda e: mp.sqrt(e) if e > 0 else 0) * aih
            gh = herm(g, lambda e: mp.sqrt(e) if e > 0 else 0)
            S = gh * A * gh
            ph = lambda X: a * mp.re(X[0, 0]) + b * mp.re(X[1, 1])        # noqa: E731
            val = ph(S) / mp.sqrt(ph(A) * ph(B)) - 1
            assert abs(val - ((a + b) / (2 * mp.sqrt(a * b)) - 1)) < 1e-15


# a point found by the search (n = 3, delta = 0.9; scripts/trace_characterizations.py): W, and A, B as U diag(lam) U^*
SGM3_W = [0.04761904762866577, 0.0476190476185667, 0.9047619047527675]
SGM3_A = ([[-0.9123051700855526, 0.06115997619819226, 0.07071974129893822],
           [-0.2032372541423999, -0.44052350979772115, 0.19072325051940764],
           [-0.012739330304503001, -0.4797602383286061, 0.20981795790812288]],
          [[0.2406429126846855, 0.17974083118539633, -0.2621863351066056],
           [0.2502565389856886, -0.11317876942190475, 0.8079730988111637],
           [0.07546030676294323, -0.7258873246604728, -0.43936042348172444]],
          [1.0, 0.003432720424771716, 5.001713679379445e-09])
SGM3_B = ([[-0.33261224999677996, 0.4130532617450929, 0.2147840805368778],
           [-0.30885624308181064, -0.2663941622080147, -0.13999717399278477],
           [0.15499493107728465, 0.021257394915701895, -0.5579137570724311]],
          [[0.3659874328990947, 0.18178247467317138, 0.7110782121000309],
           [-0.06955196367515866, -0.8463722754219032, 0.30472817144299824],
           [0.7944615830300903, -0.09261954309688579, -0.15655461506106316]],
          [0.002348505422540201, 1.0, 1.7362754367985517e-07])


def test_sgm_exceeds_kantorovich_in_m3():
    """Spectral geometric mean in M_3: at a stored point (Delta(W) = 0.9) the violation, recomputed in mpmath with 60
    digits at the exactly unitary polar factors, exceeds kappa - 1 = 1.294 (V > 2.1): the envelope of item 6 depends
    on n; it stays below the trivial bound a/b - 1 = 18."""
    import mpmath as mp
    with mp.workdps(60):
        def herm(X, f):
            E, Q = mp.eighe((X + X.transpose_conj()) / 2)
            return Q * mp.diag([f(e) for e in E]) * Q.transpose_conj()

        def build(data):
            re, im, lam = data
            U = mp.matrix([[mp.mpc(a, b) for a, b in zip(r, i)] for r, i in zip(re, im)])
            U = U * herm(U.transpose_conj() * U, lambda e: 1 / mp.sqrt(e))
            return U * mp.diag([mp.mpf(x) for x in lam]) * U.transpose_conj()
        A, B = build(SGM3_A), build(SGM3_B)
        w = [mp.mpf(x) for x in SGM3_W]
        ah, aih = herm(A, mp.sqrt), herm(A, lambda e: 1 / mp.sqrt(e))
        g = aih * herm(ah * B * ah, lambda e: mp.sqrt(e) if e > 0 else 0) * aih
        gh = herm(g, lambda e: mp.sqrt(e) if e > 0 else 0)
        S = gh * A * gh
        ph = lambda X: sum(w[i] * mp.re(X[i, i]) for i in range(3))     # noqa: E731
        v = ph(S) / mp.sqrt(ph(A) * ph(B)) - 1
        d = (max(w) - min(w)) / (max(w) + min(w))
        assert abs(d - mp.mpf('0.9')) < 1e-9
        assert v > 2.1 > float(tc.kantorovich(0.9))
        assert v < 18
        assert min(mp.eighe(A)[0]) > 0 and min(mp.eighe(B)[0]) > 0


def test_gardner_scipy_reaches_closed_form():
    """Independent search: BFGS over complex X in M_2 and M_3 reaches kappa - 1 for Gardner's inequality (to 1e-2
    relative: the maximiser has a rank-deficient |X|, where |X| is not smooth) and does not exceed it."""
    rng = np.random.default_rng(7)
    for n in (2, 3):
        d = 0.5
        w = density(d, n, middle=[0.7])

        def f(x):
            return -value_np('gardner', w, (x[:n * n] + 1j * x[n * n:]).reshape(n, n))
        best = max(-minimize(f, rng.normal(size=2 * n * n), method='BFGS').fun for _ in range(30))
        k = float(tc.kantorovich(d))
        assert best <= k + 1e-9 and best > k * (1 - 1e-2)       # BFGS stalls at the rank-deficient |X|


# ----------------------------------------------------------------------------------------------------------------
# the library path
# ----------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('name', tc.ORDER)
def test_library_objective_matches_mpmath_and_numpy(name):
    """At random admissible points of the parametrisation (n = 3, W from SpectrumParam), the library objective, the
    mpmath recomputation (certify) and plain NumPy agree."""
    torch.manual_seed(11)
    alg = tc.TypeIAlgebra([3], [3], precision='double', device='cpu')
    B = 4
    total = 3.0 if name == 'det' else 1.0
    W = tc.Weights(alg, B, delta=torch.full((B,), 0.6, dtype=torch.float64), total=total, starts=B)
    f, params, variables = tc.build(name, alg, B, W)
    with torch.no_grad():
        v = f().numpy()
        w = W.w().numpy()
    for i in range(B):
        point = dict(w=w[i], vars=[x.export(i) if isinstance(x, tc.Var) else dict(kind='scalar', t=float(x.value()[i]))
                                   for x in variables])
        c, delta, _ = tc.certify(name, point, [3], 30)
        assert abs(delta - 0.6) < 1e-12
        if np.isfinite(v[i]):
            assert abs(c - v[i]) < 1e-8 * max(1.0, abs(c)), (name, i, c, v[i])


def test_envelope_reaches_kantorovich_and_vanishes_at_the_trace():
    """optimize-based envelope of Gardner's inequality on M_2: no violation at delta = 0, V(0.5) = kappa - 1."""
    V, final, pts, _ = tc.envelope('gardner', [2], [0.0, 0.5], 8, 6, 50, 0)
    assert final[0].max() < 1e-8
    assert abs(final[1].max() - float(tc.kantorovich(0.5))) < 1e-6
    assert abs(tc.certify('gardner', pts[1], [2], 30)[0] - float(tc.kantorovich(0.5))) < 1e-6


def test_sectors_central_density_is_a_trace():
    """On M_2 (+) M_3 a density that is scalar on each block (a trace) gives V <= 0, also for a large inter-block
    contrast; an intra-block contrast delta in M_2 gives V = kappa - 1; the determinant inequality (canonical trace)
    is violated by central densities, with V = 1 - GM/AM."""
    fixed, meta = tc.sector_densities([0.5])
    V, final, pts, _ = tc.envelope('gardner', tc.SECTOR_KS, None, 6, 4, 40, 1, fixed=fixed)
    for v, m in zip(final.max(1), meta):
        if m['type'] == 'central':
            assert v < 1e-8 and m['delta'] > 0.5
        else:
            assert abs(v - float(tc.kantorovich(m['delta_intra']))) < 1e-5
    Vd, fd, _, _ = tc.envelope('det', tc.SECTOR_KS, None, 8, 6, 60, 2, fixed=fixed)
    assert np.allclose(fd.max(1), tc.det_central(meta), atol=1e-6)


def test_quick_run_and_summary():
    """run() for two inequalities with the quick configuration (one process) produces the summary used by the paper."""
    res = tc.run('quick', 'cpu', 0, log=lambda s: None, names=['gardner', 'det'], workers=1)
    S = res['summary']
    assert S['V0'] < 1e-6
    assert S['sharp_dev'] < 1e-2
    assert S['sector_central'] < 1e-6
    assert S['sector_det'] < 1e-4
