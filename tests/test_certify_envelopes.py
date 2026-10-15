"""Certified upper envelopes of Exp. 1 and Exp. 3 (scripts/certify_envelopes.py): closed forms, reductions, bounds."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location('certify_envelopes', ROOT / 'scripts' / 'certify_envelopes.py')
ce = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ce)

D = ce.DMAX
B = float(ce.lo(D))


def haar(rng, n, complex_valued=False):
    A = rng.normal(size=(n, n)) + (1j * rng.normal(size=(n, n)) if complex_valued else 0)
    Q, R = np.linalg.qr(A)
    return Q * (np.diag(R) / np.abs(np.diag(R)))


def tn(A):
    return np.linalg.svd(A, compute_uv=False).sum()


def z1_sup_u(X, Y):
    """sup over U of Exp. 1 for given X, Y: ||XY||_1 - Tr(XY) (Theorem 1 (i))."""
    return tn(X @ Y) - np.trace(X @ Y).real


def z3_any(X, Y):
    """Exp. 3 for Hermitian X, Y (real or complex), computed directly."""
    w, V = np.linalg.eigh(X)
    absX = (V * np.abs(w)) @ V.conj().T
    return np.trace(Y @ absX @ Y).real - tn(Y @ X @ Y)


def test_lemma3_double_root():
    """Lemma 3: R(2T - R) - w = 4 a (1 - b^2)^2 (c - c*)^2 with R = s3(a, b), c* = (1 + ab)/((1 + a)(1 + b))."""
    sp = pytest.importorskip('sympy')
    a, b, c = sp.symbols('a b c', positive=True)
    beta = 1 - b ** 2
    T = a * (b ** 2 + beta * c) + (b ** 2 + beta * (1 - c))
    w = 4 * a * beta ** 2 * c * (1 - c)
    R = 2 * a / (1 + a) * (1 - b) ** 2
    cs = (1 + a * b) / ((1 + a) * (1 + b))
    assert sp.simplify(R * (2 * T - R) - w - 4 * a * beta ** 2 * (c - cs) ** 2) == 0


def test_k2_closed_forms_on_dense_grid():
    """Theorem 1 (iii) and Lemma 3: on a 2x2 block the sup over theta is s1(a, b) (Exp. 1) and s3(a, b) (Exp. 3)."""
    th = np.linspace(0, np.pi / 2, 20001)
    for delta, Dy in ((0.3, 0.999), (0.8, 0.5), (0.95, 0.9)):
        a, b = float(ce.lo(delta)), float(ce.lo(Dy))
        q = np.stack([np.cos(th), np.sin(th)], -1)
        Y = b * np.eye(2) + (1 - b) * q[:, :, None] * q[:, None, :]
        # Exp. 1 block: ||XY||_1 - Tr(XY) with X = diag(1, a), from singular values (independent of g(c))
        XY = np.diag([1.0, a]) @ Y
        v1 = np.linalg.svd(XY, compute_uv=False).sum(-1) - np.einsum('nii->n', XY)
        assert v1.max() <= ce.s1(delta, Dy) + 1e-12
        assert v1.max() >= ce.s1(delta, Dy) - 1e-8
        v3 = ce.z3(np.diag([a, -1.0]), Y)
        assert v3.max() <= ce.s3(delta, Dy) + 1e-12
        assert v3.max() >= ce.s3(delta, Dy) - 1e-5                     # square-root peak as b -> 0
        ts = ce.theta_star(delta, Dy)                                  # attained at theta*
        qs = np.array([np.cos(ts), np.sin(ts)])
        assert ce.z3(np.diag([a, -1.0]), b * np.eye(2) + (1 - b) * np.outer(qs, qs)) == \
            pytest.approx(ce.s3(delta, Dy), abs=1e-12)


def test_branch_and_bound_k2_matches_closed_form():
    """Certificate for Z_1 (k = 2): lower <= s3 <= upper <= s3 + tol + EPS."""
    tol = 1e-7
    for delta in (0.0, 0.4, 0.95):
        r = ce.certify_blocks(1, delta, tol=tol)
        ref = float(ce.s3(delta))
        assert r['lower'] <= ref + 1e-12
        assert ref <= r['upper'] <= ref + tol + ce.EPS


def test_branch_and_bound_k4_equals_direct_sum():
    """Certificate for Z_2 (k = 4): the direct sum s3 + s0 is optimal up to the tolerance."""
    tol = 1e-6
    for delta in (0.2, 0.8):
        r = ce.certify_blocks(2, delta, tol=tol)
        ref = float(ce.direct_sum(3, 4, delta))
        assert abs(r['lower'] - ref) < 1e-12
        assert ref <= r['upper'] <= ref + tol + ce.EPS
        c = ce.certified_sup(3, 4, delta, tol=tol)
        assert c['upper'] >= ref


def test_box_bound_dominates_points_in_box():
    """UB(B) >= z(p) for random points p of random boxes B (m = 1, 2, 3), at random contrasts."""
    rng = np.random.default_rng(0)
    for m in (1, 2, 3):
        d = 2 * m - 1
        for _ in range(20):
            a, b = rng.uniform(0.01, 1.0), rng.uniform(0.01, 1.0)
            c = rng.uniform(0, np.pi / 2, d)
            r = rng.uniform(0, 0.4, d)
            lo_, hi_ = (c - r)[None], (c + r)[None]
            zc, ub, _ = ce.box_bounds(lo_, hi_, a, b)
            p = c + r * rng.uniform(-1, 1, (300, d))
            p[:3] = c + r * np.sign(rng.normal(size=(3, d)))        # corners
            X, _, Y, _, _ = ce.block_model(p[:, :m], p[:, m:], a, b)
            assert np.all(ce.z3(X, Y) <= ub[0] + 1e-12)
            assert zc[0] <= ub[0]


def test_second_derivative_bound():
    """|Phi''| along segments in the domain is at most curvature_bound (the certificate's M), by finite differences."""
    rng = np.random.default_rng(1)
    for m in (1, 2, 3):
        d = 2 * m - 1
        for _ in range(30):
            a, b = rng.uniform(0.01, 1.0), rng.uniform(0.01, 1.0)
            p0, p1 = rng.uniform(0, np.pi / 2, (2, d))
            if rng.random() < 0.5:
                p0[m:] = p0[m:] * 0.05                              # near the coordinate singularity phi_1 = 0
            c = (p0 + p1) / 2
            X, absX, Y, _, _ = ce.block_model(c[None, :m], c[None, m:], a, b)
            w, V = np.linalg.eigh(Y[0] @ X[0] @ Y[0])
            T = (V * np.sign(w)) @ V.T
            t = np.linspace(0, 1, 401)
            p = p0 + t[:, None] * (p1 - p0)
            X, absX, Y, _, _ = ce.block_model(p[:, :m], p[:, m:], a, b)
            phi = np.einsum('nii->n', absX @ Y @ Y) - np.einsum('ij,nji->n', T, Y @ X @ Y)
            second = np.abs(np.diff(phi, 2)) / (t[1] - t[0]) ** 2
            h = np.abs(p1 - p0)
            M, share = ce.curvature_bound(h[None, :m], h[None, m:], np.maximum(p0, p1)[None, m:], a, b)
            assert second.max() <= M[0] * (1 + 1e-3) + 1e-6
            assert share.sum() >= M[0] / 2 - 1e-12


def test_lemma2_monotone_in_abs_x_and_convex_in_y2():
    """Lemma 2 (b): z grows with |X| commuting with sign(X); (c): z is convex in the eigenvalues of Y^2."""
    rng = np.random.default_rng(2)
    for n in (3, 4):
        for _ in range(30):
            V, W = haar(rng, n, True), haar(rng, n, True)
            s = np.where(rng.random(n) < 0.5, -1.0, 1.0)
            x = rng.uniform(0.1, 1, n)
            x2 = x + rng.uniform(0, 1, n) * (1 - x)                    # |X| <= |X'|, same eigenbasis and signs
            y = rng.uniform(0.1, 1, n)
            Y = (W * y) @ W.conj().T
            X, X2 = (V * (s * x)) @ V.conj().T, (V * (s * x2)) @ V.conj().T
            assert z3_any(X, Y) <= z3_any(X2, Y) + 1e-12
            t0, t1 = rng.uniform(0.01, 1, n), rng.uniform(0.01, 1, n)  # eigenvalues of Y^2
            Ys = [(W * np.sqrt(t)) @ W.conj().T for t in (t0, t1, (t0 + t1) / 2)]
            assert z3_any(X, Ys[2]) <= (z3_any(X, Ys[0]) + z3_any(X, Ys[1])) / 2 + 1e-12


def test_certified_bounds_dominate_random_complex_configurations():
    """Theorems 1 and 4: z <= certified sup for random complex configurations (two-level and general spectra)."""
    rng = np.random.default_rng(3)
    for k in (2, 3, 4):
        for delta in (0.0, 0.5, 0.9):
            ub1 = ce.certified_sup(1, k, delta)['upper']
            ub3 = ce.certified_sup(3, k, delta, tol=1e-5)['upper']
            a = float(ce.lo(delta))
            for _ in range(60):
                two = rng.random() < 0.5
                x = np.where(rng.random(k) < 0.5, 1.0, a) if two else rng.uniform(a, 1, k)
                x[:2] = 1.0, a
                bmin = rng.uniform(B, 1)
                y = np.where(rng.random(k) < 0.5, 1.0, bmin) if two else rng.uniform(bmin, 1, k)
                y[0] = 1.0
                V, W = haar(rng, k, True), haar(rng, k, True)
                Y = (W * y) @ W.conj().T
                assert z1_sup_u((V * x) @ V.conj().T, Y) <= ub1 + 1e-12
                s = np.where(rng.random(k) < 0.5, -1.0, 1.0)
                assert z3_any((V * (s * x)) @ V.conj().T, Y) <= ub3 + 1e-12


def test_exp3_k3_exceeds_k2_envelope():
    """Theorem 4 (ii), k = 3: X = diag(1, -1, a), Y = b + (1 - b) q q^* on the first two coordinates attains s0 > s3."""
    for delta in (0.3, 0.9):
        a = float(ce.lo(delta))
        q = np.array([1.0, 1.0, 0.0]) / np.sqrt(2)
        Y = B * np.eye(3) + (1 - B) * np.outer(q, q)
        v = z3_any(np.diag([1.0, -1.0, a]), Y)
        assert abs(v - ce.s0()) < 1e-12
        assert v > ce.s3(delta) + 0.1
        assert ce.certified_sup(3, 3, delta)['upper'] == pytest.approx(v, abs=1e-12)


@pytest.mark.filterwarnings('ignore:optimize.extremize')
def test_certified_bounds_dominate_multistart_search():
    """The multistart search of torch_vn_algebra.optimize stays below the certified sup (and reaches it for k = 2)."""
    pytest.importorskip('torch')
    deltas = [0.3, 0.8]
    for exp in (1, 3):
        for k in (2, 3, 4):
            found = ce.search_sup(exp, k, deltas, starts=6, rounds=2, steps=15, seed=k)
            for delta, v in zip(deltas, found):
                r = ce.certified_sup(exp, k, delta, tol=1e-6)
                assert v <= r['upper'] + 1e-9
                if k == 2:
                    assert v >= r['lower'] - 0.05 * max(r['lower'], 0.05)


def test_relaxation_is_valid_but_loose():
    """The concavity (averaging) bound is >= the certified Z_2 and exceeds the direct sum near delta = 1."""
    r = ce.certify_blocks(2, 0.99, tol=1e-6)
    rel = ce.relaxation_bound(2, 0.99)
    assert rel >= r['upper'] - ce.EPS
    assert rel - ce.direct_sum(3, 4, 0.99) > 0.3


def test_several_sectors():
    """Theorems 1 and 4 for C = 2: random block-diagonal complex configurations stay below the certified sup, and the
    direct sum attains it for k = 2."""
    rng = np.random.default_rng(4)
    for k in (2, 3):
        for delta in (0.3, 0.9):
            ub1 = ce.certified_sup(1, k, delta, C=2)['upper']
            ub3 = ce.certified_sup(3, k, delta, C=2)['upper']
            assert ub3 == pytest.approx(float(ce.direct_sum(3, k, delta, C=2)), abs=1e-12)
            a = float(ce.lo(delta))
            for _ in range(40):
                x = np.where(rng.random(2 * k) < 0.5, 1.0, a)
                x[rng.permutation(2 * k)[:2]] = 1.0, a
                y = np.where(rng.random(2 * k) < 0.5, 1.0, B)
                y[0] = 1.0
                V, W = np.zeros((2 * k, 2 * k), complex), np.zeros((2 * k, 2 * k), complex)
                for c in range(2):
                    sl = slice(c * k, (c + 1) * k)
                    V[sl, sl], W[sl, sl] = haar(rng, k, True), haar(rng, k, True)
                Y = (W * y) @ W.conj().T
                s = np.where(rng.random(2 * k) < 0.5, -1.0, 1.0)
                assert z1_sup_u((V * x) @ V.conj().T, Y) <= ub1 + 1e-12
                assert z3_any((V * (s * x)) @ V.conj().T, Y) <= ub3 + 1e-12
