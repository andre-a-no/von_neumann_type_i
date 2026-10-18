"""Checks of scripts/contrast_operations.py (contrasts of composite operators, paper man/trace_envelopes)."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import contrast_operations as co  # noqa: E402


def _pairs(n, S=3000, seed=1):
    rng = np.random.default_rng(seed)
    da, db = rng.uniform(0, 0.999, S), rng.uniform(0, 0.999, S)
    return da, db, co.positive(rng, da, n), co.positive(rng, db, n, scale=True)


@pytest.mark.parametrize('n', [2, 3, 4])
def test_multiplicative_bounds(n):
    da, db, A, B = _pairs(n)
    Ah, Aih = co.fn(A, np.sqrt), co.fn(A, lambda e: 1 / np.sqrt(e))
    for X, w in ((Ah @ B @ Ah, 1.0), (Ah @ co.fn(Aih @ B @ Aih, np.sqrt) @ Ah, 0.5)):
        c, _ = co.contrast(X)
        lo, hi = co.hyp(da, db, w)
        assert np.all(c <= hi + 1e-6) and np.all(c >= lo - 1e-6)
    c, _ = co.contrast(A + B)
    assert np.all(c <= np.maximum(da, db) + 1e-9)


@pytest.mark.parametrize('n', [2, 3, 4])
def test_jordan(n):
    da, db, A, B = _pairs(n)
    c, lam = co.contrast(A @ B + B @ A)
    inside = da ** 2 + db ** 2 <= 1
    assert np.all(lam[inside] > -1e-12)                          # Gustafson's criterion
    ok = lam > 0
    lo, hi = co.jordan_bounds(da, db)
    assert np.all(c[ok] <= hi[ok] + 1e-6) and np.all(c[ok] >= lo[ok] - 1e-6)


def test_jordan_two_by_two():
    # indefinite at c = 1/2 exactly when s^2 + t^2 > 1
    assert np.isnan(co.jordan_sup2(0.8, 0.61)) and not np.isnan(co.jordan_sup2(0.8, 0.59))
    # small contrasts: the commuting value tanh(eta_A + eta_B); larger ones: more on M_n, n >= 3, than on M_2
    assert co.jordan_sup2(0.3, 0.3) == pytest.approx(np.tanh(2 * np.arctanh(0.3)), abs=1e-9)
    assert co.jordan_sup3(0.6, 0.6) > co.jordan_sup2(0.6, 0.6) + 0.01
    assert co.jordan_sup3(0.6, 0.6) < np.sin(2 * np.arcsin(0.6))
    # the construction 1 + (2 x 2 block) attains (2 - mu)/(2 + mu) on M_3
    s = t = 0.6
    a, b = (1 - s) / (1 + s), (1 - t) / (1 + t)
    cs = np.linspace(0, 1, 20001)
    p = (1 - a) * (1 - b)
    T = a + b + p * cs
    c = cs[np.argmin(T - np.sqrt(T * T - 4 * a * b + p * p * cs * (1 - cs)))]
    th = np.arccos(np.sqrt(c))
    v = np.array([np.cos(th), np.sin(th)])
    B2 = b * np.eye(2) + (1 - b) * np.outer(v, v)
    A3, B3 = np.diag([1.0, 1.0, a]), np.eye(3)
    B3[1:, 1:] = B2
    e = np.linalg.eigvalsh(A3 @ B3 + B3 @ A3)
    assert (e[-1] - e[0]) / (e[-1] + e[0]) == pytest.approx(co.jordan_sup3(s, t), abs=1e-6)


@pytest.mark.parametrize('n', [2, 3, 4, 5])
def test_commutator_trace_test_is_sharp(n):
    rng = np.random.default_rng(n)
    dA, dB, dw = 0.6, 0.4, 0.7
    # bound on random samples
    S = 3000
    A, B = co.positive(rng, np.full(S, dA), n, scale=True), co.positive(rng, np.full(S, dB), n, scale=True)
    w = co.spectrum(rng, np.full(S, dw), n)
    w = w / w.sum(1, keepdims=True)
    v = np.abs(np.einsum('si,sii->s', w, A @ B - B @ A)) / (co.mid(A) * co.mid(B))
    assert np.all(v <= 2 * co.trace_test_factor(n, dw) * dA * dB + 1e-9)
    # the extremal configuration: k blocks m(1 + d sigma_z), m(1 + d sigma_x), W aligned with i[A, B]
    k = n // 2
    sz, sx = np.diag([1.0, -1.0]), np.array([[0.0, 1.0], [1.0, 0.0]])
    A, B, U = np.eye(n, dtype=complex), np.eye(n, dtype=complex), np.eye(n, dtype=complex)
    for j in range(k):
        sl = slice(2 * j, 2 * j + 2)
        A[sl, sl] = np.eye(2) + dA * sz
        B[sl, sl] = np.eye(2) + dB * sx
        _, U[sl, sl] = np.linalg.eigh(1j * (A[sl, sl] @ B[sl, sl] - B[sl, sl] @ A[sl, sl]))   # ascending
    h = (1 + dw) / (1 - dw)
    bmin = 1 / (k * h + n - k)
    wd = np.full(n, bmin)
    for j in range(k):
        wd[2 * j + 1] = h * bmin                                    # top eigenvalue of W on the positive eigenvector
    W = U @ np.diag(wd) @ U.conj().T
    assert np.trace(W).real == pytest.approx(1)
    val = abs(np.trace(W @ (A @ B - B @ A)))
    assert val == pytest.approx(2 * co.trace_test_factor(n, dw) * dA * dB, rel=1e-9)
