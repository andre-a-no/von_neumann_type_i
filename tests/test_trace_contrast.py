"""Closed forms of man/trace_envelopes for two prescribed contrasts (scripts/trace_contrast.py), checked with plain
NumPy: Gardner's inequality |phi_W(X)| <= c phi_W(|X|) with Delta(W) = d, Delta(|X|) = d' (Theorem 4.1; random points
never exceed it, the 2x2 maximiser attains it) and the Young lower bound (Proposition 5.2, attained by A = B)."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
tcon = pytest.importorskip('trace_contrast')


def _spec(rng, n, d):
    l = (1 - d) / (1 + d)
    x = rng.uniform(l, 1, n)
    x[0], x[1] = 1.0, l
    return x


def _haar(rng, n):
    q, r = np.linalg.qr(rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)))
    return q * (np.diag(r) / abs(np.diag(r)))


@pytest.mark.parametrize('n', [2, 3, 4])
def test_gardner_two_contrasts_is_an_upper_bound(n):
    rng = np.random.default_rng(n)
    for _ in range(3000):
        d, dp = rng.uniform(0, 0.98, 2)
        w = rng.permutation(_spec(rng, n, d))
        U = _haar(rng, n)
        T = U @ np.diag(_spec(rng, n, dp)) @ U.conj().T
        TW = T @ np.diag(w)
        ratio = np.linalg.svd(TW, compute_uv=False).sum() / np.trace(TW).real    # sup over unitaries / phi(|X|)
        assert ratio - 1 <= tcon.gardner_two(d, dp, n) + 1e-9


@pytest.mark.parametrize('d,dp', [(0.3, 0.2), (0.6, 0.6), (0.9, 0.95)])
def test_gardner_two_contrasts_attained(d, dp):
    a, b = 1 + d, 1 - d
    l = (1 - dp) / (1 + dp)
    S, q = a * l + b, (a + b) * (1 + l)
    g = S / q
    R = np.array([[np.sqrt(g), -np.sqrt(1 - g)], [np.sqrt(1 - g), np.sqrt(g)]])
    TW = R @ np.diag([1, l]) @ R.T @ np.diag([a, b])
    ratio = np.linalg.svd(TW, compute_uv=False).sum() / np.trace(TW)
    assert abs(ratio - 1 - tcon.gardner_two(d, dp, 2)) < 1e-12
    assert tcon.gardner_two(d, dp, 3) < tcon.gardner_two(d, dp, 2)


@pytest.mark.parametrize('d,dp', [(0.3, 0.2), (0.6, 0.6), (0.9, 0.95)])
def test_young_two_contrasts_construction(d, dp):
    a, b = 1 + d, 1 - d
    l = (1 - dp) / (1 + dp)
    W = np.diag([a, b])
    A = np.array([[0, 1], [1, 0]]) @ np.diag([l, 1.0])           # A = U T, T = diag(l, 1), U the swap
    lhs = np.trace(W @ A @ A.T)                                   # phi(|A A^*|) = phi(A A^*)
    rhs = np.trace(W @ A.T @ A)                                   # phi(|A|^2)
    assert abs(lhs / rhs - 1 - tcon.young_two(d, dp)) < 1e-12
