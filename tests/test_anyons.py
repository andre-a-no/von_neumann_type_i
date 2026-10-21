"""Chains of SU(2)_k anyons and Temperley-Lieb algebras (torch_vn_algebra.anyons)."""
import math

import numpy as np
import pytest
import torch

from torch_vn_algebra import AnyonChain


@pytest.mark.parametrize('k,L,periodic', [(2, 6, False), (3, 7, False), (4, 6, True), (3, 6, True)])
def test_temperley_lieb_relations(k, L, periodic):
    ch = AnyonChain(k, L, periodic=periodic)
    d = ch.d
    E = [ch.generator(i).matrix for i in range(1, ch.n_generators + 1)]
    n = len(E)
    nb = (lambda i, j: (i - j) % n in (1, n - 1)) if periodic else (lambda i, j: abs(i - j) == 1)
    for i in range(n):
        assert torch.allclose(E[i] @ E[i], d * E[i], atol=1e-12)
        for j in range(n):
            if i == j:
                continue
            if nb(i, j):
                assert torch.allclose(E[i] @ E[j] @ E[i], E[i], atol=1e-12)
            else:
                assert torch.allclose(E[i] @ E[j], E[j] @ E[i], atol=1e-12)


def test_fibonacci_sizes_and_total_charge():
    ch = AnyonChain(3, 10)                       # Fibonacci anyons: sector sizes are Fibonacci numbers
    assert ch.charges == [0, 2] and ch.algebra.k_factors == [34, 55]
    assert ch.d == pytest.approx((1 + math.sqrt(5)) / 2)
    odd = AnyonChain(3, 9)
    assert odd.charges == [1, 3] and odd.algebra.k_factors == [34, 21]


def test_tower_jones_index_markov_trace_and_basic_construction():
    k, L = 3, 9
    inc = AnyonChain(k, L).inclusion()
    small, big = inc.smaller, AnyonChain(k, L)
    for i in range(1, L - 1):                                        # TL_{L-1} sits in TL_L in standard form
        assert torch.allclose(inc.embed(small.generator(i)).matrix, big.generator(i).matrix, atol=1e-12)
    assert inc.index() == pytest.approx(big.d ** 2)                 # 4 cos^2(pi/5) = phi^2
    t = inc.markov_trace().t.numpy()
    qd = np.array([big.quantum_dimension(b) for b in big.charges])
    assert np.allclose(t / t.sum(), qd / qd.sum())                  # Markov trace = quantum dimensions
    assert inc.basic_construction().M.k_factors == AnyonChain(k, L + 1).algebra.k_factors


def test_lanczos_agrees_with_dense():
    ch = AnyonChain(3, 14)
    dense = ch.low_energies(levels=3)
    lanczos = ch.low_energies(levels=3, dense_max=10)
    for a, b in zip(dense, lanczos):
        assert torch.allclose(a, b, atol=1e-8)


def test_golden_chain_central_charge():
    """Periodic Fibonacci chain: E0(L) = e L - pi v c/(6 L) + O(L^-3) with v = pi sin(pi/5)/(pi/5) gives c = 7/10."""
    Ls = np.array([10, 12, 14, 16], dtype=float)
    E0 = np.array([float(AnyonChain(3, int(L), periodic=True).low_energies()[0][0]) for L in Ls])
    co = np.linalg.lstsq(np.stack([Ls, 1 / Ls, 1 / Ls ** 3], 1), E0, rcond=None)[0]
    c = -co[1] * 6 / (math.pi * AnyonChain.velocity(3))
    assert c == pytest.approx(AnyonChain.central_charge(3), abs=2e-3)
    assert co[0] == pytest.approx(1 - math.sqrt(5), abs=1e-5)       # energy per site of the golden chain


def test_open_chain_gap_ratio_is_a_kac_dimension():
    """Ising (k = 2): the lowest level of the sector b = 1 (spin 1) over the first excitation of b = 0 tends to
    h_{1,3}/2 = 1/4 (the vacuum tower has its first descendant at level 2)."""
    r = []
    for L in (12, 16, 20):
        e = AnyonChain(2, L).low_energies(levels=2)
        r.append(float((e[1][0] - e[0][0]) / (e[0][1] - e[0][0])))
    assert abs(r[-1] - AnyonChain.kac(2, 1, 3) / 2) < 3e-3 and abs(r[-1] - 0.25) < abs(r[0] - 0.25)
