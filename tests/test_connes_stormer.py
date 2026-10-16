"""
Connes-Stormer relative entropy of two MASAs (paper, Section "Use case: the Connes-Stormer entropy of two MASAs";
scripts/connes_stormer.py): the Lemma (reduction to ensembles of states), the Proposition for n = 2 (closed form,
strict inequality over Choda's value, the bound kappa theta), against brute-force partitions in NumPy/SciPy, and the
library path (inclusions, conditional expectations, optimize) against plain NumPy. CPU, double precision, fixed seeds.
"""
import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest
import torch
from scipy.integrate import quad
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location('connes_stormer', ROOT / 'scripts' / 'connes_stormer.py')
cs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cs)

pytestmark = [pytest.mark.filterwarnings("ignore::torch_vn_algebra.cost.CostWarning")]


def haar_np(rng, n):
    A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    Q, R = np.linalg.qr(A)
    return Q * (np.diag(R) / np.abs(np.diag(R)))


def partition_np(Z):
    """x_j = S^{-1/2} Z_j Z_j^* S^{-1/2}, S = sum_j Z_j Z_j^* (Z: (m, n, r))."""
    G = [z @ z.conj().T for z in Z]
    w, V = np.linalg.eigh(sum(G))
    Si = (V * w ** -0.5) @ V.conj().T
    return [Si @ g @ Si for g in G]


def random_partition(rng, m, n, r=None):
    r = n if r is None else r
    return partition_np(rng.normal(size=(m, n, r)) + 1j * rng.normal(size=(m, n, r)))


def brute_force_n2(U, rng, m=3, starts=6):
    """sup over partitions of M_2 with m elements (full rank), by BFGS from random starts (no use of the Lemma)."""
    best = -np.inf
    k = 2 * 2 * m
    for _ in range(starts):
        f = lambda v: -cs.cs_numpy(U, partition_np((v[:k] + 1j * v[k:]).reshape(m, 2, 2)))   # noqa: E731
        r = minimize(f, rng.normal(size=2 * k), method='BFGS', options={'gtol': 1e-10, 'maxiter': 2000})
        best = max(best, -r.fun)
    return best


@pytest.mark.parametrize('c', [0.97, 0.85, 0.7, 0.3])
def test_proposition_n2_against_brute_force_partitions(c):
    """Proposition (n = 2): H = max_alpha [h(cos^2(alpha - theta)) - h(cos^2 alpha)] for complex U with phases; the
    maximum over general partitions (BFGS) reaches the formula and random partitions never exceed it."""
    rng = np.random.default_rng(int(1000 * c))
    U = cs.rotation_with_phases(c, rng)
    H, alpha = cs.formula_n2(c)
    assert abs(brute_force_n2(U, rng) - H) < 1e-7
    for _ in range(300):
        assert cs.cs_numpy(U, random_partition(rng, int(rng.integers(2, 6)), 2)) <= H + 1e-12
    # the optimal partition of the Proposition: the basis rotated by alpha, moved by the left phases of U
    D1 = np.diag(np.exp(1j * np.angle(U[:, 0])))
    v, w = D1 @ np.array([math.cos(alpha), math.sin(alpha)]), D1 @ np.array([-math.sin(alpha), math.cos(alpha)])
    assert abs(cs.cs_numpy(U, [np.outer(v, v.conj()), np.outer(w, w.conj())]) - H) < 1e-13


def test_proposition_n2_phases_do_not_matter():
    """H depends on U only through c = |u_11|^2: the brute-force value is the same for different phases."""
    rng = np.random.default_rng(5)
    vals = [brute_force_n2(cs.rotation_with_phases(0.8, rng), rng, starts=4) for _ in range(2)]
    assert abs(vals[0] - vals[1]) < 1e-7


def test_n2_strictly_above_choda():
    """H > h(c) (Choda's value: the minimal projections of Delta_2) unless c in {0, 1/2, 1}, where H = h(c)."""
    for c in np.linspace(0.02, 0.98, 25):
        if abs(c - 0.5) < 1e-9:
            continue
        H, _ = cs.formula_n2(c)
        assert H > float(cs.h2(c)) + 1e-6
    for c, val in ((0.0, 0.0), (0.5, math.log(2)), (1.0, 0.0)):
        H, _ = cs.formula_n2(c)
        assert abs(H - val) < 1e-12 and abs(float(cs.h2(c)) - val) < 1e-12
    # derivative of phi at alpha = 0: log(tan^2 theta) sin(2 theta), non-zero iff c not in {0, 1/2, 1}
    th, e = math.acos(math.sqrt(0.8)), 1e-6
    d = (cs.phi2(e, th) - cs.phi2(-e, th)) / (2 * e)
    assert abs(d - math.log(math.tan(th) ** 2) * math.sin(2 * th)) < 1e-6


def test_n2_window_integral_and_kappa_bound():
    """phi(alpha, theta) = int_{alpha - theta}^{alpha} psi, psi(b) = sin(2b) log tan^2 b; hence H <= kappa theta,
    H / theta -> kappa as theta -> 0, and H(c) = H(1 - c)."""
    kappa, _ = cs.kappa2()
    assert abs(kappa - 1.3254868) < 1e-6
    for alpha, th in ((-0.3, 0.2), (0.4, 1.1), (-1.0, 0.7)):
        kinks = [j * math.pi / 2 for j in range(-4, 5) if alpha - th < j * math.pi / 2 < alpha]   # beta log beta
        val = quad(lambda b: float(cs.psi2(b)), alpha - th, alpha, points=kinks or None, limit=200,
                   epsabs=1e-13, epsrel=1e-13)[0]
        assert abs(val - float(cs.phi2(alpha, th))) < 1e-10
    for c in (0.999, 0.95, 0.8, 0.6):
        H, _ = cs.formula_n2(c)
        th = math.acos(math.sqrt(c))
        assert H <= kappa * th + 1e-12
        assert abs(H - cs.formula_n2(1 - c)[0]) < 1e-10
    th = math.acos(math.sqrt(1 - 1e-8))
    assert abs(cs.formula_n2(1 - 1e-8)[0] / th - kappa) < 1e-3


def test_bloch_symmetry_and_reduction():
    """n = 2 in Bloch coordinates: x_j = t_j (1 + r_j . sigma) gives F = sum_j t_j f(r_j), f(r) = k(r.q) - k(r.z),
    k(s) = h((1 + s)/2), and f(-r) = f(r); the series of k has non-positive coefficients."""
    rng = np.random.default_rng(3)
    sig = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])]
    U = cs.rotation_with_phases(0.7, rng)
    q = np.real([np.trace(U @ np.diag([1, 0]) @ U.conj().T @ s) for s in sig])
    z = np.array([0.0, 0.0, 1.0])
    k = lambda s: float(cs.h2((1 + s) / 2))                                                   # noqa: E731
    f = lambda r: k(r @ q) - k(r @ z)                                                         # noqa: E731
    for _ in range(50):
        r = rng.normal(size=3)
        r *= rng.random() / np.linalg.norm(r)
        t = rng.random()
        x = t * (np.eye(2) + sum(ri * s for ri, s in zip(r, sig)))
        assert abs(cs.cs_numpy(U, [x]) - t * f(r)) < 1e-13
        assert abs(f(r) - f(-r)) < 1e-14
    s = 0.73
    series = math.log(2) - sum(s ** (2 * m) / (2 * m * (2 * m - 1)) for m in range(1, 400))
    assert abs(series - k(s)) < 1e-12


def test_lemma_ensemble_identity_and_upper_bound():
    """Lemma: with x_j = n t_j rho_j, sum_j F(x_j) = sum_j t_j [S(E_Q rho_j) - S(E_P rho_j)] <= max_j G(rho_j)."""
    rng = np.random.default_rng(7)
    for n in (2, 3, 4):
        U = haar_np(rng, n)

        def G(rho):
            S = lambda d: float(cs.eta_np(d).sum())                                            # noqa: E731
            return S(np.real(np.diag(U.conj().T @ rho @ U))) - S(np.real(np.diag(rho)))

        for _ in range(20):
            m = int(rng.integers(2, 2 * n))
            xs = random_partition(rng, m, n, r=int(rng.integers(-(-n // m), n + 1)))   # m r >= n: S invertible
            assert np.abs(sum(xs) - np.eye(n)).max() < 1e-12 and min(np.linalg.eigvalsh(x).min() for x in xs) > -1e-12
            t = [np.trace(x).real / n for x in xs]
            rhos = [x / np.trace(x).real for x in xs]
            val = cs.cs_numpy(U, xs)
            assert abs(sum(t) - 1) < 1e-12
            assert abs(val - sum(tj * G(r) for tj, r in zip(t, rhos))) < 1e-12
            assert val <= max(G(r) for r in rhos) + 1e-12
            assert val <= math.log(n) + 1e-12


def test_self_and_complementary_masas():
    """H(P|P) = 0 (the functional vanishes identically); for the Fourier matrix, the minimal projections of Delta_n
    give log n, which is the upper bound S(E_Q rho) - S(E_P rho) <= log n of the Lemma (n = 2, 3)."""
    rng = np.random.default_rng(11)
    for n in (2, 3):
        for _ in range(5):
            assert abs(cs.cs_numpy(np.eye(n), random_partition(rng, 4, n))) < 1e-15
        F = cs.fourier(n)
        E = [np.diag(np.eye(n)[j]).astype(complex) for j in range(n)]
        assert abs(cs.cs_numpy(F, E) - math.log(n)) < 1e-13
        assert abs(cs.unistochastic_entropy(F) - math.log(n)) < 1e-13
        for _ in range(50):
            assert cs.cs_numpy(F, random_partition(rng, 5, n)) <= math.log(n) + 1e-12
        # the library: a short ascent from random partitions reaches log n
        best, _, _ = cs.maximise(torch.as_tensor(F)[None], 'partition', m=n * n, starts=4, rounds=4, steps=60)
        assert abs(float(best[0]) - math.log(n)) < 1e-5


def test_library_path_equals_numpy():
    """Pair.terms (Inclusion with unitary, conditional_expectation('norm'), spectra of E(x)) equals the plain NumPy
    functional, for MASAs and for a non-abelian subalgebra (M_2 (+) C c M_3 against its blocks)."""
    rng = np.random.default_rng(13)
    for n in (2, 3, 5):
        U = haar_np(rng, n)
        xs = random_partition(rng, 6, n)
        pair = cs.Pair(torch.as_tensor(U)[None].expand(len(xs), n, n))
        lib = float(pair.terms(torch.as_tensor(np.stack(xs))).sum())
        assert abs(lib - cs.cs_numpy(U, xs)) < 1e-12
    # non-abelian: P = M_2 (+) C in M_3 with the isometries V_c of the embedding: E_P(x)_c = V_c^* x V_c for tau = Tr/3
    # (nu_c = w = 1/3), so tau eta(E_P x) = (1/3) sum_c Tr eta(V_c^* x V_c), and V_c -> U V_c for Q = U P U^*
    U = haar_np(rng, 3)
    xs = random_partition(rng, 4, 3)
    pair = cs.Pair(torch.as_tensor(U)[None].expand(4, 3, 3), sizes=[2, 1], lam=[[1, 1]])
    V0 = pair.P.isometries()[0, 0, 0].numpy()                         # (C, 3, k_max)
    Vs = [V0[0][:, :2], V0[1][:, :1]]
    assert np.allclose(sum(v @ v.conj().T for v in Vs), np.eye(3))

    def tau_eta(x, Uq):
        return sum(float(cs.eta_np(np.linalg.eigvalsh((Uq @ v).conj().T @ x @ (Uq @ v))).sum()) for v in Vs) / 3

    ref = sum(tau_eta(x, U) - tau_eta(x, np.eye(3)) for x in xs)
    lib = float(pair.terms(torch.as_tensor(np.stack(xs))).sum())
    assert abs(lib - ref) < 1e-12


def test_certification_matches_double():
    """certify (mpmath, re-normalised partition and unitary) agrees with the double-precision value."""
    pytest.importorskip("mpmath")
    rng = np.random.default_rng(17)
    U = haar_np(rng, 3)
    Z = rng.normal(size=(4, 3, 1)) + 1j * rng.normal(size=(4, 3, 1))
    val, res, ures = cs.certify(U, Z, digits=30)
    assert abs(val - cs.cs_numpy(U, partition_np(Z))) < 1e-12
    assert res < 1e-25 and ures < 1e-25


def test_quick_comparison_n3_runs():
    """The core of the script (quick configuration) for n = 3 on two Haar unitaries: H >= Choda's value, the bound
    G* >= H of the Lemma, partitions >= bases up to round-off, valid structure of the result."""
    pytest.importorskip("mpmath")
    torch.manual_seed(0)
    cfg = dict(cs.CONFIGS['quick'], digits=20)
    out = cs.compare(3, 2, cfg, 'cpu', 0, log=lambda s: None)
    H = np.maximum(out['H_general'], out['H_bases'])
    assert out['nU'] == 2 and len(out['H_general']) == 2
    assert (H >= np.array(out['choda']) - 1e-9).all()
    assert (np.array(out['G_star']) >= H - 1e-6).all()
    assert (np.array(H) <= math.log(3) + 1e-9).all()
    for c in out['cases']:
        assert c['partition_residual'] < 1e-15


def test_quick_n2_part_runs():
    """part_n2 (quick): the library optimum over general partitions matches the closed form."""
    cfg = dict(cs.CONFIGS['quick'], n2_table=[0.9], n2_starts=4, rounds=4, steps=60)
    out = cs.part_n2(cfg, 'cpu', 0, log=lambda s: None)
    assert out['max_dev_library'] < 1e-6
    assert out['max_dev_numpy'] < 1e-13


def test_paper_listing_runs():
    """Listing lst:cs of the paper (extracted from usecases/connes_stormer.tex, with fewer starts and steps) runs and gives H >= h(U)."""
    tex = (ROOT / 'man' / 'article' / 'usecases' / 'connes_stormer.tex').read_text(encoding='utf-8')
    start = tex.index('label=lst:cs]')
    code = tex[tex.index('\n', start) + 1:tex.index('\\end{lstlisting}', start)]
    assert 'n, B, m = 3, 64, 9' in code and 'rounds=8, steps=150' in code
    code = code.replace('n, B, m = 3, 64, 9', 'n, B, m = 3, 4, 9').replace('rounds=8, steps=150', 'rounds=3, steps=40')
    torch.manual_seed(0)
    ns = {}
    exec(compile(code, 'lst_cs', 'exec'), ns)
    assert float(ns['H']) >= float(ns['h']) - 1e-9
    U = ns['U'].matrix[0, 0].numpy()
    assert float(ns['H']) <= math.log(3) + 1e-12
    assert abs(float(ns['h']) - cs.unistochastic_entropy(U)) < 1e-13
