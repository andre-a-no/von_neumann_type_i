"""
Chapter 5, Section 5.4: the merged spectrum of many independent blocks (Lemmas 5.16, 5.17 and Theorem 5.18).

Lemma 5.16: GOE and GUE blocks have a bounded Lipschitz level density and P(delta <= eps) <= K eps^2.
Lemma 5.17: around a fixed point x, the merged spectrum of m independent blocks has
            |P(L, R < inf, r <= t) - 2t/(1+t)| <= A/m, via the exact formula P(R < inf, L > tR) = M - E_1.
Theorem 5.18: pooled over samples, the fraction of the levels in a window I with spacing ratio <= t and their mean
              ratio differ from 2t/(1+t) and 2 ln 2 - 1 (independent levels, Proposition 5.14(a)) by at most A/C.
Monte Carlo checks use z-scores (5 standard errors); the rates are checked as "C * deviation does not grow".
"""
import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import norm

POISSON_MEAN = 2 * math.log(2) - 1


def _goe(rng, n, k):
    """Eigenvalues (ascending) of n independent GOE(k) matrices H = (G + G^T)/2."""
    g = rng.standard_normal((n, k, k))
    return np.linalg.eigvalsh((g + g.transpose(0, 2, 1)) / 2)


def _gue(rng, n, k):
    """Eigenvalues (ascending) of n independent GUE(k) matrices H = (G + G^*)/sqrt(2), G Ginibre."""
    g = (rng.standard_normal((n, k, k)) + 1j * rng.standard_normal((n, k, k))) / math.sqrt(2)
    return np.linalg.eigvalsh((g + g.conj().transpose(0, 2, 1)) / math.sqrt(2))


def _merged(levels, n, C):
    """Merge C consecutive blocks of `levels` (shape (n*C, k)) into n sorted spectra of C*k levels."""
    return np.sort(levels.reshape(n, -1), axis=1)


def _gaps(spec, x):
    """L_m(x), R_m(x): distances from x to the nearest level below / above (inf if none), per sample."""
    R = np.where(spec > x, spec, np.inf).min(1) - x
    L = x - np.where(spec < x, spec, -np.inf).max(1)
    return L, R


def _pooled(spec, lo, hi, ts):
    """Pooled E S(I)/E N(I) and E N_t(I)/E N(I) over the samples (rows), with delta-method standard errors."""
    s = np.diff(spec, axis=1)
    r = np.minimum(s[:, :-1], s[:, 1:]) / np.maximum(s[:, :-1], s[:, 1:])    # r(y) at the inner levels
    y = spec[:, 1:-1]
    N = ((spec >= lo) & (spec <= hi)).sum(1)                                 # N(I), all levels in I
    inner = (y >= lo) & (y <= hi)

    def ratio(S):
        R = S.sum() / N.sum()
        return R, math.sqrt(((S - R * N) ** 2).sum()) / N.sum()

    mean = ratio((r * inner).sum(1))
    cdf = {t: ratio(((r <= t) & inner).sum(1)) for t in ts}
    return mean, cdf


def test_lemma_5_16_goe_gue_level_density_is_lipschitz():
    """Lemma 5.16 (H1): the GOE(2) and GUE(2) level densities of Corollary C.8 are bounded and Lipschitz, and match MC."""
    rng = np.random.default_rng(0)
    for beta, sampler in ((1, _goe), (2, _gue)):
        w = lambda a, b, beta=beta: abs(a - b) ** beta * math.exp(-(a * a + b * b) / 2)
        Z = quad(lambda a: quad(lambda b: w(a, b), a, np.inf)[0], -np.inf, np.inf)[0]   # Z_{2,beta}, ordered
        rho = lambda x, w=w, Z=Z: 2 * quad(lambda y: w(x, y), -np.inf, np.inf)[0] / (2 * Z)   # n int q(x, y) dy
        assert abs(quad(rho, -np.inf, np.inf)[0] - 2) < 1e-6                           # int rho = n
        xs = np.linspace(-4, 4, 161)
        vals = np.array([rho(x) for x in xs])
        assert vals.max() < 1                                                          # bounded
        for h in (1e-2, 1e-3):                                                         # difference quotients
            q = np.array([abs(rho(x + h) - rho(x)) / h for x in xs])
            assert q.max() < 2.0                                                       # bounded, uniformly in h
        lam = sampler(rng, 200_000, 2).ravel()                                         # MC: E N([a, b]) = int rho
        for a, b in ((-0.5, 0.5), (1.0, 2.0)):
            p = ((lam >= a) & (lam <= b)).reshape(-1, 2).sum(1)
            mean, se = p.mean(), p.std() / math.sqrt(p.size)
            assert abs(mean - quad(rho, a, b)[0]) <= 5 * se


def test_lemma_5_16_smallest_spacing_quadratic_bound():
    """Lemma 5.16 (H2) and (3): for GOE(4), P(delta <= eps)/eps^2 stays bounded as eps -> 0, and c^2 E exp(-c delta) <= 2K."""
    rng = np.random.default_rng(1)
    lam = _goe(rng, 400_000, 4)
    delta = np.diff(lam, axis=1).min(1)
    ratios = {}
    for eps in (0.4, 0.1, 0.05, 0.025):
        p = (delta <= eps).mean()
        ratios[eps] = (p / eps ** 2, math.sqrt(p * (1 - p) / delta.size) / eps ** 2)
    K = max(v for v, _ in ratios.values())
    (v1, e1), (v2, e2) = ratios[0.05], ratios[0.025]
    assert abs(v1 - v2) <= 5 * math.hypot(e1, e2)                   # beta = 1: P(delta <= eps)/eps^2 -> constant
    assert v2 <= 1.5 * ratios[0.4][0]                                # ... and stays bounded
    merged = np.sort(lam[:200_000].reshape(100_000, 8), axis=1)      # two blocks merged (Theorem 5.15): no bound
    dm = np.diff(merged, axis=1).min(1)
    assert (dm <= 0.025).mean() / 0.025 ** 2 > 5 * v2                # P(delta <= eps) ~ eps, not eps^2
    eps_grid = np.linspace(0.01, 3, 300)                           # an empirical K for (H2) on a grid of eps
    Kall = max(K, max((delta <= e).mean() / e ** 2 for e in eps_grid))
    for c in (5.0, 20.0, 80.0):                                     # (3): E e^{-c delta} <= 2K/c^2
        v = np.exp(-c * delta)
        assert v.mean() - 5 * v.std() / math.sqrt(v.size) <= 2 * Kall / c ** 2


def _M_gauss(x, m, t):
    """M = m int_0^inf rho(x+b) sigma_b^{m-1} db for one N(0,1) level per block (k = 1), sigma_b = 1 - rho(J_b)."""
    f = lambda b: m * norm.pdf(x + b) * (1 - (norm.cdf(x + b) - norm.cdf(x - t * b))) ** (m - 1)
    s = 1 / (m * (1 + t) * norm.pdf(x))
    cuts = (0, s, 10 * s, 60 * s, 40)
    return sum(quad(f, a, b, limit=200)[0] for a, b in zip(cuts[:-1], cuts[1:]))


def test_lemma_5_17_formula_exact_for_one_level_blocks():
    """Lemma 5.17 (formula): for k = 1 (E_1 = 0), P(R < inf, L > tR) = M exactly; Monte Carlo vs quadrature."""
    rng = np.random.default_rng(2)
    x, n = 0.3, 300_000
    for m, t in ((3, 0.5), (8, 0.2), (8, 0.9)):
        L, R = _gaps(rng.standard_normal((n, m)), x)
        hit = np.isfinite(R) & (L > t * R)
        p, se = hit.mean(), hit.std() / math.sqrt(n)
        assert abs(p - _M_gauss(x, m, t)) <= 5 * se


def test_lemma_5_17_M_rate_deterministic():
    """Lemma 5.17 (M close to 1/(1+t)): ideal integral = 1/(1+t), and m |M - 1/(1+t)| stays bounded (quadrature, k = 1)."""
    for t in (0.0, 0.3, 0.8):
        rho0, m = 0.4, 7                                            # m int_0^{1/a} rho0 (1 - ab)^{m-1} db = 1/(1+t)
        a = (1 + t) * rho0
        assert abs(quad(lambda b: m * rho0 * (1 - a * b) ** (m - 1), 0, 1 / a)[0] - 1 / (1 + t)) < 1e-12
    for t in (0.2, 0.5, 0.9):
        dev = {m: _M_gauss(0.3, m, t) - 1 / (1 + t) for m in (10, 100, 1000, 10_000)}
        scaled = [m * abs(d) for m, d in dev.items()]
        assert max(scaled) < 1.0                                    # |M - 1/(1+t)| <= A/m
        assert scaled[-1] <= 1.05 * scaled[1]                       # m * deviation does not grow
        # both sides: P(L, R < inf, r <= t) = 2 - M(x) - M'(x) up to P(L or R = inf) < 1e-20, M'(x) = M(-x) (rho even)
        p = {m: 2 - _M_gauss(0.3, m, t) - _M_gauss(-0.3, m, t) for m in (100, 1000)}
        for m, v in p.items():
            assert m * abs(v - 2 * t / (1 + t)) < 1.0


def test_lemma_5_17_fixed_point_goe_blocks():
    """Lemma 5.17: GOE(4) blocks, x = 0.5: |P(L, R < inf, r_m(x) <= t) - 2t/(1+t)| <= A/m (z-scores, m up to 32)."""
    rng = np.random.default_rng(3)
    x, n, ts = 0.5, 20_000, (0.25, 0.5, 0.75)
    dev, err = {}, {}
    for m in (2, 4, 8, 16, 32):
        L, R = _gaps(_merged(_goe(rng, n * m, 4), n, m), x)
        fin = np.isfinite(L) & np.isfinite(R)
        r = np.where(fin, np.minimum(L, R) / np.maximum(L, R), 2.0)
        for t in ts:
            hit = r <= t
            dev[m, t] = hit.mean() - 2 * t / (1 + t)
            err[m, t] = hit.std() / math.sqrt(n)
    for t in ts:
        A = 2 * abs(dev[2, t])                                      # the constant found at m = 2 ...
        for m in (4, 8, 16, 32):                                    # ... works for all m (rate 1/m)
            assert m * abs(dev[m, t]) <= 2 * A + 5 * m * err[m, t]
        assert abs(dev[32, t]) <= 5 * err[32, t] + A / 32            # the limit 2t/(1+t)


@pytest.mark.parametrize('ensemble, k, window', [('GOE', 8, 1.5), ('GUE', 4, 1.0)])
def test_theorem_5_18_pooled_ratio_statistics(ensemble, k, window):
    """Theorem 5.18: pooled over samples, E S(I)/E N(I) -> 2 ln 2 - 1 and E N_t(I)/E N(I) -> 2t/(1+t), at rate A/C."""
    rng = np.random.default_rng(4 if ensemble == 'GOE' else 5)
    sampler = _goe if ensemble == 'GOE' else _gue
    ts = (0.1, 0.3, 0.6, 0.9)
    res = {}
    for C in (1, 2, 4, 8, 16, 32):
        n = 40_000 // C
        spec = _merged(sampler(rng, n * C, k), n, C)
        res[C] = _pooled(spec, -window, window, ts)
    (m1, _), _ = res[1]
    (m2, _), _ = res[2]
    assert m1 > m2 > POISSON_MEAN                                   # repulsion is diluted, monotonically at first
    A = 2 * (m2 - POISSON_MEAN)                                     # C * deviation at C = 2
    for C in (4, 8, 16, 32):
        (mean, se), cdf = res[C]
        assert C * abs(mean - POISSON_MEAN) <= A + 5 * C * se       # |E S/E N - (2 ln 2 - 1)| <= A/C
        for t in ts:
            F, seF = cdf[t]
            F2, _ = res[2][1][t]
            A_t = 2 * abs(F2 - 2 * t / (1 + t))
            assert C * abs(F - 2 * t / (1 + t)) <= 2 * A_t + 5 * C * seF
    (mean, se), cdf = res[32]
    assert abs(mean - POISSON_MEAN) <= 5 * se + A / 32


def test_theorem_5_18_expected_number_of_levels():
    """Theorem 5.18, E N(I) = C rho(I): the expected number of merged levels in I is C times that of one block."""
    rng = np.random.default_rng(6)
    one = _goe(rng, 100_000, 8)
    N1 = ((one >= -1.5) & (one <= 1.5)).sum(1)
    for C in (3, 7):
        n = 100_000 // C
        spec = _merged(_goe(rng, n * C, 8), n, C)
        NC = ((spec >= -1.5) & (spec <= 1.5)).sum(1)
        se = math.hypot(NC.std() / math.sqrt(n), C * N1.std() / math.sqrt(N1.size))
        assert abs(NC.mean() - C * N1.mean()) <= 5 * se
