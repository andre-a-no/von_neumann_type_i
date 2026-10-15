"""
Theory tests for Chapters 9-12 of the lecture notes (spin chains, quantum computing, optimisation, equivariant
learning), run through the public API of torch_vn_algebra.

Independent references are built here: dense spin-chain Hamiltonians in the label basis (bit i of a label is
site i), explicit translation matrices, Jordan-Wigner fermions, free-fermion and one-magnon formulas, Bethe
roots, Burnside counts, Kraus-operator algebra in NumPy, and brute-force orbit counts.
"""
import itertools
import math

import numpy as np
import pytest
import scipy.linalg
import torch
from scipy.optimize import brentq

import torch_vn_algebra.decompose as dc
from torch_vn_algebra import (SpinChain, SU2Chain, TypeIAlgebra, DensityMatrix, channels, states, dynamics,
                              inclusions, krylov)
from torch_vn_algebra import optimize as opt

f64, c128 = torch.float64, torch.complex128

# ==========================================================================================================
# dense references (label basis: index b = sum_i b_i 2^i, b_i = 1 for spin up)
# ==========================================================================================================
SP = np.array([[0.0, 0.0], [1.0, 0.0]])          # s^+ : |0> = down -> |1> = up
SM = SP.T.copy()
SZ = np.diag([-0.5, 0.5])
I2 = np.eye(2)


def site_ops(L, ops):
    """Product of one-site operators {site: 2 x 2}; site 0 is the least significant bit."""
    out = np.ones((1, 1))
    for j in reversed(range(L)):
        out = np.kron(out, ops.get(j, I2))
    return out


def bond_xxz(L, i, j, J=1.0, Delta=1.0):
    return J * (0.5 * (site_ops(L, {i: SP, j: SM}) + site_ops(L, {i: SM, j: SP}))
                + Delta * site_ops(L, {i: SZ, j: SZ}))


def xxz_dense(L, bonds, J=1.0, Delta=1.0, h=None):
    H = np.zeros((2 ** L, 2 ** L))
    for i, j in bonds:
        H = H + bond_xxz(L, i, j, J, Delta)
    if h is not None:
        h = np.broadcast_to(np.asarray(h, dtype=float), (L,))
        for i in range(L):
            H = H + h[i] * site_ops(L, {i: SZ})
    return H


def popcount(x):
    x = np.asarray(x)
    return np.array([bin(int(v)).count('1') for v in x.reshape(-1)]).reshape(x.shape)


def rot(x, L, j=1):
    """T^j on labels: T moves site i to site i + 1 (mod L)."""
    j %= L
    mask = (1 << L) - 1
    x = np.asarray(x, dtype=np.int64)
    return ((x << j) | (x >> (L - j))) & mask if j else x


def translation(L):
    T = np.zeros((2 ** L, 2 ** L))
    x = np.arange(2 ** L)
    T[rot(x, L), x] = 1
    return T


def to_full(chain, op, b=0):
    """Dense 2^L matrix of an operator of SpinChain.algebra (block diagonal in the sectors)."""
    L = chain.L
    out = np.zeros((2 ** L, 2 ** L), dtype=complex)
    for N in range(L + 1):
        lab = chain.labels[N].numpy()
        k = len(lab)
        out[np.ix_(lab, lab)] = op.matrix[b, N, :k, :k].numpy()
    return out


def jw_fermions(L):
    """c_j = (prod_{l<j} (-1)^{n_l}) s^-_j as dense 2^L matrices."""
    x = np.arange(2 ** L)
    cs = []
    for j in range(L):
        string = np.diag((-1.0) ** popcount(x & ((1 << j) - 1)))
        cs.append(string @ site_ops(L, {j: SM}))
    return cs


def ground_signs_ok(vec, nB):
    x = (-1.0) ** nB * vec
    x = x * np.sign(x.sum())
    return x.min() > 0


# ==========================================================================================================
# Chapter 9. Spin chains
# ==========================================================================================================
@pytest.mark.parametrize('L,bc', [(2, 'open'), (3, 'open'), (3, 'periodic'), (5, 'periodic'), (6, 'open'),
                                  (6, 'periodic')])
def test_9_1_xxz_blocks_are_sector_restrictions(L, bc):
    """Proposition 9.1: [H, S^z] = 0 and SpinChain.xxz is the full XXZ matrix restricted to the sectors N."""
    rng = np.random.default_rng(L)
    J, Delta, h = 1.3, rng.normal(), rng.normal(size=L)
    chain = SpinChain(L, bc)
    H = xxz_dense(L, chain.bonds, J, Delta, h)
    Sz = sum(site_ops(L, {i: SZ}) for i in range(L))
    assert np.abs(H @ Sz - Sz @ H).max() < 1e-12
    assert np.abs(to_full(chain, chain.xxz(J, Delta, h=torch.tensor(h))) - H).max() < 1e-12
    for N in range(L + 1):
        lab = chain.labels[N].numpy()
        assert sorted(lab.tolist()) == sorted(np.nonzero(popcount(np.arange(2 ** L)) == N)[0].tolist())
        assert np.abs(chain.xxz(J, Delta, h=torch.tensor(h), sector=N).matrix[0, 0].numpy()
                      - H[np.ix_(lab, lab)]).max() < 1e-12


def test_9_2_tower_inclusion_matrix_and_embedding():
    """Proposition 9.2: Lambda_L has lambda_{N,N} = lambda_{N,N-1} = 1, and x -> x (x) 1_2 is the embedding."""
    tower = inclusions.chain_tower(5)
    gen = torch.Generator().manual_seed(0)
    for L, inc in enumerate(tower):
        lam = np.asarray(inc.matrix)
        expect = np.zeros((L + 2, L + 1), dtype=int)
        for N in range(L + 1):
            expect[N, N] = expect[N + 1, N] = 1
        assert lam.shape == (L + 2, L + 1) and (lam == expect).all()
        assert [math.comb(L + 1, N) for N in range(L + 2)] == list(lam @ np.array(inc.N.k_factors))  # Pascal
        # the bipartite graph is a path on 2L + 3 vertices: connected, degrees <= 2, 2L + 2 edges
        B = np.block([[np.zeros((L + 2, L + 2)), lam], [lam.T, np.zeros((L + 1, L + 1))]])
        deg = B.sum(1)
        assert B.sum() / 2 == 2 * L + 2 and deg.max() <= 2 and (deg == 1).sum() == 2 and inc.is_connected
        # the embedding x -> x (x) 1_2 (new site L = top bit), sector by sector in the SpinChain bases
        blocks = [torch.randn(1, k, k, dtype=c128, generator=gen) for k in inc.N.k_factors]
        x = inc.N.from_blocks(blocks)
        if L == 0:
            x_full = blocks[0][0].numpy()
        else:
            x_full = to_full(SpinChain(L), x)
        X = np.kron(I2, x_full)
        Y = inc.embed(x)
        big = SpinChain(L + 1)
        for N in range(L + 2):
            lab = big.labels[N].numpy()
            k = len(lab)
            assert np.abs(Y.matrix[0, N, :k, :k].numpy() - X[np.ix_(lab, lab)]).max() < 1e-12


@pytest.mark.parametrize('n', [1, 2, 3, 7, 12])
def test_9_3_path_spectrum(n):
    """Lemma 9.3: the path A_n has eigenvalues 2 cos(pi k / (n+1)) with eigenvectors sin(pi j k / (n+1))."""
    A = np.eye(n, k=1) + np.eye(n, k=-1)
    j = np.arange(1, n + 1)
    for k in range(1, n + 1):
        v = np.sin(np.pi * j * k / (n + 1))
        assert np.abs(A @ v - 2 * np.cos(np.pi * k / (n + 1)) * v).max() < 1e-12
    assert abs(np.linalg.norm(A, 2) - 2 * np.cos(np.pi / (n + 1))) < 1e-12


@pytest.mark.parametrize('L', range(0, 9))
def test_9_4_index_of_tower(L):
    """Corollary 9.4: [A_{L+1} : A_L] = ||Lambda_L||^2 = 4 cos^2(pi / (2L + 4)), increasing to 4."""
    inc = inclusions.chain_tower(L + 1, L_min=L)[0]
    val = 4 * math.cos(math.pi / (2 * L + 4)) ** 2
    assert abs(inc.index() - val) < 1e-12
    assert abs(np.linalg.norm(np.asarray(inc.matrix, dtype=float), 2) ** 2 - val) < 1e-12
    assert val < 4 * math.cos(math.pi / (2 * L + 6)) ** 2 < 4


@pytest.mark.parametrize('L', [3, 4, 5, 6, 8, 9, 10])
def test_9_5_momentum_dims_burnside_and_orbits(L):
    """Proposition 9.5(b): d_{N,m} = #{orbits: m R_a = 0 mod L}, also = (1/L) sum_j e^{-2 pi i m j/L} Tr(T^j P_N)."""
    ring = SpinChain(L, 'periodic')
    x = np.arange(2 ** L)
    pop = popcount(x)
    for N in range(L + 1):
        fixed = [int(((rot(x, L, j) == x) & (pop == N)).sum()) for j in range(L)]
        burnside = [sum(fixed[j] * np.exp(-2j * np.pi * m * j / L) for j in range(L)) / L for m in range(L)]
        assert max(abs(b - round(b.real)) for b in burnside) < 1e-9
        orbits = {frozenset(rot(a, L, j).item() for j in range(L)) for a in x[pop == N]}
        assert all(L % len(o) == 0 for o in orbits)                                          # (a): R_a | L
        d = [sum((m * len(o)) % L == 0 for o in orbits) for m in range(L)]
        assert d == ring.momentum_dims(N) == [round(b.real) for b in burnside]
        assert sum(d) == math.comb(L, N)


def test_9_5_momentum_blocks_against_translation_eigenspaces():
    """Proposition 9.5: block m of xxz_momentum is H on the e^{2 pi i m/L}-eigenspace of T inside sector N."""
    L, J, Delta, h = 6, 1.1, 0.7, 0.3
    ring = SpinChain(L, 'periodic')
    H = xxz_dense(L, ring.bonds, J, Delta, h)
    T = translation(L)
    assert np.abs(T @ H @ T.T - H).max() < 1e-12 and np.abs(np.linalg.matrix_power(T, L) - np.eye(2 ** L)).max() == 0
    for N in range(L + 1):
        Hk = ring.xxz_momentum(J, Delta, h=h, sector=N)
        lab = ring.labels[N].numpy()
        TN, HN = T[np.ix_(lab, lab)], H[np.ix_(lab, lab)]
        for c, m in enumerate(Hk.algebra.charges):
            Q = sum(np.exp(-2j * np.pi * m * j / L) * np.linalg.matrix_power(TN, j) for j in range(L)) / L
            w, V = np.linalg.eigh(Q)
            B = V[:, w > 0.5]
            assert B.shape[1] == Hk.algebra.k_factors[c]
            ref = np.linalg.eigvalsh(B.conj().T @ HN @ B)
            assert np.abs(Hk.eigenvalues()[c][0].numpy() - ref).max() < 1e-10


def test_9_6_perron_frobenius():
    """Lemma 9.6: a symmetric matrix with nonnegative off-diagonal entries and connected graph has a simple top
    eigenvalue with a positive eigenvector."""
    rng = np.random.default_rng(7)
    for n in (2, 5, 9):
        for _ in range(5):
            M = np.diag(rng.normal(size=n))
            for i in range(1, n):                                      # a random spanning tree
                j = rng.integers(0, i)
                M[i, j] = M[j, i] = rng.uniform(0.1, 1)
            extra = rng.random((n, n)) < 0.2
            W = np.triu(extra * rng.uniform(0, 1, (n, n)), 1)
            M = M + W + W.T
            w, V = np.linalg.eigh(M)
            v = V[:, -1] * np.sign(V[:, -1].sum())
            assert v.min() > 0 and (n == 1 or w[-1] - w[-2] > 1e-10)


@pytest.mark.parametrize('trial', range(3))
def test_9_7_marshall_general_bipartite_graph(trial):
    """Theorem 9.7: for bipartite sign-constrained J_ij on a connected graph, arbitrary K_ij and h_i, every sector
    has a simple ground state with sign (-1)^{N_B(b)}."""
    rng = np.random.default_rng(100 + trial)
    n = 7
    A = set(rng.choice(n, size=int(rng.integers(2, 6)), replace=False).tolist())
    order = rng.permutation(n)
    edges = set()
    for i in range(1, n):
        a, b = int(order[i]), int(order[rng.integers(0, i)])
        edges.add((min(a, b), max(a, b)))
    for a, b in itertools.combinations(range(n), 2):
        if rng.random() < 0.25:
            edges.add((a, b))
    H = np.zeros((2 ** n, 2 ** n))
    for i, j in itertools.combinations(range(n), 2):
        if (i, j) in edges:
            Jij = (-1.0 if ((i in A) == (j in A)) else 1.0) * rng.uniform(0.3, 1.5)
            H += 0.5 * Jij * (site_ops(n, {i: SP, j: SM}) + site_ops(n, {i: SM, j: SP}))
        H += rng.normal() * site_ops(n, {i: SZ, j: SZ})                     # K_ij arbitrary on all pairs
    for i in range(n):
        H += rng.normal() * site_ops(n, {i: SZ})
    maskB = sum(1 << i for i in range(n) if i not in A)
    x = np.arange(2 ** n)
    for N in range(n + 1):
        idx = x[popcount(x) == N]
        w, V = np.linalg.eigh(H[np.ix_(idx, idx)])
        assert ground_signs_ok(V[:, 0], popcount(idx & maskB))
        assert len(w) == 1 or w[1] - w[0] > 1e-9


@pytest.mark.parametrize('L,bc', [(7, 'open'), (8, 'open'), (8, 'periodic'), (10, 'periodic')])
def test_9_7_marshall_xxz_chain_library(L, bc):
    """Theorem 9.7 for the XXZ chain (open, or periodic with L even) with random fields and anisotropy."""
    gen = torch.Generator().manual_seed(L)
    c = SpinChain(L, bc, complex_valued=False)
    for Delta in (-0.8, 0.4, 2.0):
        for N in range(L + 1):
            h = 0.5 * torch.randn(L, dtype=f64, generator=gen)
            w, V = torch.linalg.eigh(c.xxz(1.0, Delta, h=h, sector=N).matrix[0, 0])
            nB = c.bits(N)[:, 1::2].sum(-1).numpy()
            assert ground_signs_ok(V[:, 0].numpy(), nB)
            assert len(w) == 1 or (w[1] - w[0]).item() > 1e-9


@pytest.mark.parametrize('L', [4, 6, 8])
def test_9_8_ground_state_momentum(L):
    """Corollary 9.8: on an even ring with uniform field, T psi_N = (-1)^N psi_N for the ground state of sector N."""
    ring = SpinChain(L, 'periodic')
    for Delta, h in ((0.5, 0.0), (1.0, 0.3), (1.7, -0.2)):
        for N in range(L + 1):
            w, V = torch.linalg.eigh(ring.xxz(1.0, Delta, h=h, sector=N).matrix[0, 0])
            psi = V[:, 0].numpy()
            lab = ring.labels[N].numpy()
            Tpsi = np.zeros_like(psi)
            Tpsi[ring.index_of[torch.as_tensor(rot(lab, L))].numpy()] = psi        # (T psi)(T b) = psi(b)
            assert np.abs(Tpsi - (-1) ** N * psi).max() < 1e-10
            if 0 < N < L:
                Hk = ring.xxz_momentum(1.0, Delta, h=h, sector=N)
                lowest, charges = [e[0, 0].item() for e in Hk.eigenvalues()], Hk.algebra.charges
                assert charges[int(np.argmin(lowest))] == (0 if N % 2 == 0 else L // 2)
                assert abs(min(lowest) - w[0].item()) < 1e-10


def test_9_3_su2_multiplicities():
    """Proposition 4.13 (Section 9.3): m_S = binom(L, L/2-S) - binom(L, L/2-S-1), sum (2S+1) m_S = 2^L."""
    for L in range(1, 11):
        su = SU2Chain(L)
        ms = [math.comb(L, int(L / 2 - S)) - (math.comb(L, int(L / 2 - S - 1)) if L / 2 - S - 1 >= 0 else 0)
              for S in su.spins]
        assert su.multiplicities == ms
        assert sum(int(2 * S + 1) * m for S, m in zip(su.spins, su.multiplicities)) == 2 ** L


@pytest.mark.parametrize('L,bc', [(5, 'open'), (6, 'periodic')])
def test_9_3_su2_blocks_give_full_spectrum(L, bc):
    """Section 9.3: H = (+)_S 1_{2S+1} (x) H_S, also with a uniform field (levels E + h(M - S))."""
    su = SU2Chain(L, bc, complex_valued=False)
    for h in (0.0, 0.37):
        E = su.full_spectrum(su.heisenberg(J=1.2, h=h))[0].numpy()
        ref = np.linalg.eigvalsh(xxz_dense(L, su.bonds, 1.2, 1.0, h))
        assert np.abs(E - ref).max() < 1e-10


@pytest.mark.parametrize('L,bc,S_star', [(7, 'open', 0.5), (8, 'open', 0.0), (8, 'periodic', 0.0),
                                         (9, 'open', 0.5)])
def test_9_9_lieb_mattis_chain(L, bc, S_star):
    """Theorem 9.9: E_0(S_*) < E_0(S_*+1) < ... < E_0(L/2), E_0(S) simple, and the ground state of the sector
    S^z = M has spin max(M, S_*)."""
    su = SU2Chain(L, bc, complex_valued=False)
    H = su.heisenberg(J=1.0)
    E0 = su.lowest_energies(H)[0]
    i0 = su.spins.index(S_star)
    assert (E0[i0 + 1:] > E0[i0:-1] + 1e-9).all() and (E0[:i0] > E0[i0] + 1e-9).all()
    for c, e in enumerate(H.eigenvalues()):
        if su.spins[c] >= S_star and e.shape[1] > 1:
            assert (e[0, 1] - e[0, 0]).item() > 1e-9                                   # simple
    gs = su.ground_state(H)
    assert gs.spin.item() == S_star
    for N in range((L + 1) // 2, L + 1):                                               # M = N - L/2 >= 0
        M = N - L / 2
        lowest = torch.linalg.eigvalsh(su.u1.xxz(1.0, 1.0, sector=N).matrix[0, 0])[0]
        assert abs(lowest - E0[su.spins.index(max(M, S_star))]).item() < 1e-10


def test_9_9_lieb_mattis_bipartite_graph():
    """Theorems 9.7 and 9.9 on a bipartite graph with n_A = 5, n_B = 2 (S_* = 3/2): ground multiplet of spin 3/2."""
    n, A = 7, {0, 1, 2, 3, 4}
    edges = [(0, 5), (1, 5), (2, 6), (3, 6), (4, 5), (4, 6), (0, 1), (2, 3), (5, 6), (1, 3)]
    gen = torch.Generator().manual_seed(5)
    J = torch.tensor([(1.0 if (i in A) != (j in A) else -1.0) * (0.4 + torch.rand(1, generator=gen).item())
                      for i, j in edges], dtype=f64)
    su = SU2Chain(n, complex_valued=False)
    H = su.heisenberg(J=J, bonds=edges)
    E0 = su.lowest_energies(H)[0]
    i0 = su.spins.index(1.5)
    assert (E0[i0 + 1:] > E0[i0:-1] + 1e-9).all() and (E0[:i0] > E0[i0] + 1e-9).all()
    Hd = sum(Jb * bond_xxz(n, i, j) for Jb, (i, j) in zip(J.tolist(), edges))
    w = np.linalg.eigvalsh(Hd)
    assert w[3] - w[0] < 1e-10 and w[4] - w[0] > 1e-9 and abs(w[0] - E0[i0].item()) < 1e-10     # 2 S_* + 1 = 4
    Sp = sum(site_ops(n, {i: SP}) for i in range(n))
    Sz = sum(site_ops(n, {i: SZ}) for i in range(n))
    S2 = Sp.T @ Sp + Sz @ Sz + Sz
    x = np.arange(2 ** n)
    maskB = (1 << 5) | (1 << 6)
    for N in range(4, n + 1):                                                          # M = N - 7/2 >= 1/2
        idx = x[popcount(x) == N]
        ww, U = np.linalg.eigh(Hd[np.ix_(idx, idx)])
        s = max(N - 3.5, 1.5)
        assert ground_signs_ok(U[:, 0], popcount(idx & maskB)) and (len(ww) == 1 or ww[1] - ww[0] > 1e-9)
        assert abs(U[:, 0] @ S2[np.ix_(idx, idx)] @ U[:, 0] - s * (s + 1)) < 1e-9


@pytest.mark.parametrize('L', [4, 6, 8, 10])
def test_9_9_correlations_and_sum_rules(L):
    """Section 9.3 and Exercise 2: <s_0.s_r> = 3<s^z_0 s^z_r>, sign (-1)^r, E_0 = J L <s_0.s_1>, sum_r = 0."""
    su = SU2Chain(L, 'periodic', complex_valued=False)
    H = su.heisenberg(J=1.0)
    gs = su.ground_state(H)
    S0 = gs.spin.item()
    assert S0 == 0
    coef = gs.coefficients[0, :su.multiplicity(0)]
    pairs = [(0, r) for r in range(L)]
    corr = su.expectation(su.spin_correlation(pairs), 0, coef)
    corr_d = su.expectation(su.spin_correlation(pairs, method='direct'), 0, coef)
    assert (corr - corr_d).abs().max() < 1e-10
    psi = su.embed(0, coef).numpy()
    lab = su.u1.labels[L // 2].numpy()
    full = np.zeros(2 ** L)
    full[lab] = psi
    for r in range(1, L):
        zz = full @ site_ops(L, {0: SZ, r: SZ}) @ full
        dot = full @ bond_xxz(L, 0, r) @ full
        assert abs(dot - corr[r].item()) < 1e-10 and abs(dot - 3 * zz) < 1e-10
        assert (-1) ** r * corr[r].item() > 0
    assert abs(corr[0].item() - 0.75) < 1e-12
    assert abs(L * corr[1].item() - gs.energy.item()) < 1e-10 and abs(corr.sum().item()) < 1e-10


@pytest.mark.parametrize('L', [1, 2, 4, 7, 9])
def test_9_10_xx_chain_is_free(L):
    """Theorem 9.10(c): every level of sector N of the open XX chain is a sum of N distinct J cos(pi k/(L+1))."""
    J = 0.8
    xx = SpinChain(L)
    eps = [J * math.cos(math.pi * k / (L + 1)) for k in range(1, L + 1)]
    for N in range(L + 1):
        E = torch.linalg.eigvalsh(xx.xxz(J=J, Delta=0.0, sector=N).matrix[0, 0]).numpy()
        ref = np.sort([sum(K) for K in itertools.combinations(eps, N)])
        assert np.abs(E - ref).max() < 1e-12


@pytest.mark.parametrize('L', [2, 5, 7])
def test_9_10_jordan_wigner(L):
    """Theorem 9.10(a),(b): the CAR, c_j^* c_j = n_j, and H_XX = sum_ij h_ij c_i^* c_j with h = (J/2) A_L."""
    cs = jw_fermions(L)
    one = np.eye(2 ** L)
    for i in range(L):
        for j in range(L):
            assert np.abs(cs[i] @ cs[j].T + cs[j].T @ cs[i] - (i == j) * one).max() < 1e-12
            assert np.abs(cs[i] @ cs[j] + cs[j] @ cs[i]).max() < 1e-12
        assert np.abs(cs[i].T @ cs[i] - site_ops(L, {i: SZ}) - 0.5 * one).max() < 1e-12
    J = 1.4
    h = J / 2 * (np.eye(L, k=1) + np.eye(L, k=-1))
    Hf = sum(h[i, j] * cs[i].T @ cs[j] for i in range(L) for j in range(L))
    chain = SpinChain(L)
    assert np.abs(to_full(chain, chain.xxz(J, 0.0)) - Hf).max() < 1e-12


@pytest.mark.parametrize('L', [6, 8])
def test_9_16_peschel_excited_states(L):
    """Theorems 9.10(c) and 9.16: b_K^*|0> is an eigenvector with energy sum_K eps_k, and its entanglement entropy
    is sum_nu eta(nu) over the eigenvalues of C_A (SVD, Peschel and the library agree)."""
    rng = np.random.default_rng(L)
    J = 1.0
    cs = jw_fermions(L)
    j = np.arange(1, L + 1)
    u = [np.sqrt(2 / (L + 1)) * np.sin(np.pi * j * k / (L + 1)) for k in range(1, L + 1)]   # Lemma 9.3
    eps = [J * math.cos(math.pi * k / (L + 1)) for k in range(1, L + 1)]
    chain = SpinChain(L)
    Hfull = to_full(chain, chain.xxz(J, 0.0)).real
    vac = np.zeros(2 ** L)
    vac[0] = 1.0
    eta = lambda nu: -sum(v * math.log(v) + (1 - v) * math.log(1 - v) for v in nu if 1e-13 < v < 1 - 1e-13)
    for N in (1, L // 2, L - 1):
        K = sorted(rng.choice(L, size=N, replace=False).tolist())
        Phi = vac.copy()
        for k in reversed(K):
            Phi = sum(u[k][i] * cs[i].T for i in range(L)) @ Phi
        assert abs(np.linalg.norm(Phi) - 1) < 1e-12
        assert np.abs(Hfull @ Phi - sum(eps[k] for k in K) * Phi).max() < 1e-12
        C = sum(np.outer(u[k], u[k]) for k in K)
        for i in range(L):                                                        # C_ij = <c_i^* c_j>
            assert abs(Phi @ cs[i].T @ cs[(i + 1) % L] @ Phi - C[i, (i + 1) % L]) < 1e-12
        lab = chain.labels[N].numpy()
        v = torch.tensor(Phi[lab], dtype=c128)
        rho = DensityMatrix(chain.sector_algebra(N), matrix=torch.outer(v, v.conj())[None, None])
        for n in range(1, L):
            s = np.linalg.svd(Phi.reshape(2 ** (L - n), 2 ** n), compute_uv=False) ** 2   # rows: sites >= n
            S_svd = -sum(p * math.log(p) for p in s if p > 1e-15)
            nu = np.linalg.eigvalsh(C[:n, :n])
            assert nu.min() > -1e-12 and nu.max() < 1 + 1e-12
            assert abs(S_svd - eta(nu)) < 1e-9
            assert abs(chain.entanglement_entropy(rho, n, sector=N).item() - S_svd) < 1e-9


@pytest.mark.parametrize('L,J,Delta', [(3, 1.0, 0.5), (5, 1.3, -0.7), (8, 0.6, 1.0)])
def test_9_11_one_magnon(L, J, Delta):
    """Proposition 9.11: the ring sector N = 1 has energies E_F + J(cos k - Delta), k = 2 pi m / L, E_F = J Delta L/4."""
    ring = SpinChain(L, 'periodic')
    EF = J * Delta * L / 4
    assert abs(ring.xxz(J, Delta, sector=0).matrix[0, 0, 0, 0].real.item() - EF) < 1e-12
    E = torch.linalg.eigvalsh(ring.xxz(J, Delta, sector=1).matrix[0, 0]).numpy()
    ks = 2 * np.pi * np.arange(L) / L
    assert np.abs(E - np.sort(EF + J * (np.cos(ks) - Delta))).max() < 1e-12
    Hk = ring.xxz_momentum(J, Delta, sector=1)
    for c, m in enumerate(Hk.algebra.charges):
        assert Hk.algebra.k_factors[c] == 1
        assert abs(Hk.matrix[0, c, 0, 0].real.item() - (EF + J * (math.cos(2 * math.pi * m / L) - Delta))) < 1e-12


def _bethe_vectors(L, Delta):
    """Roots of F(q) on (0, pi) and on q = i kappa, pi + i kappa, with the Bethe vectors and energies (J = 1)."""
    x1, x2 = np.triu_indices(L, 1)

    def roots(g, grid):
        v = np.array([g(t) for t in grid])
        return [brentq(g, grid[i], grid[i + 1], xtol=1e-15) for i in np.nonzero(v[:-1] * v[1:] < 0)[0]]
    out = []
    for n in range(L):
        K, c = 2 * np.pi * n / L, math.cos(np.pi * n / L)
        F = lambda q: (-1) ** n * np.exp(1j * q * L) * (c - Delta * np.exp(-1j * q)) + c - Delta * np.exp(1j * q)
        G = lambda q: (lambda z: z.real + z.imag)(np.exp(-0.5j * q * L) * F(q))
        bound = [s + 1j * t for s in (0, np.pi) for t in roots(lambda t: F(s + 1j * t).real, np.linspace(1e-3, 20, 2000))]
        for q in roots(G, (np.arange(1, 1500) + 0.5 ** 0.5) * np.pi / 1500) + bound:
            k1, k2 = K / 2 + q, K / 2 - q
            A12, A21 = c - Delta * np.exp(1j * q), -(c - Delta * np.exp(-1j * q))
            psi = A12 * np.exp(1j * (k1 * x1 + k2 * x2)) + A21 * np.exp(1j * (k2 * x1 + k1 * x2))
            E = Delta * L / 4 + np.cos(k1) + np.cos(k2) - 2 * Delta
            out.append((K, psi / np.linalg.norm(psi), E))
    return out


@pytest.mark.parametrize('L,Delta', [(9, 0.6), (7, -0.8), (8, 1.5), (6, 0.3)])
def test_9_12_bethe_two_magnons(L, Delta):
    """Proposition 9.12: every root of F gives H Psi = (E_F + J(cos k1 + cos k2 - 2 Delta)) Psi and T Psi = e^{-iK} Psi;
    for even L the vector Phi = sum_x (-1)^x |x, x+1> has energy E_F - J Delta and is orthogonal to all of them."""
    ring = SpinChain(L, 'periodic')
    H2 = ring.xxz(1.0, Delta, sector=2).matrix[0, 0].numpy()
    x1, x2 = np.triu_indices(L, 1)
    rows = ring.index_of[torch.tensor((1 << x1) | (1 << x2))].numpy()
    lab = ring.labels[2].numpy()
    perm = ring.index_of[torch.as_tensor(rot(lab, L))].numpy()                      # T on the sector basis
    vecs = []
    for K, psi, E in _bethe_vectors(L, Delta):
        v = np.zeros(len(lab), dtype=complex)
        v[rows] = psi
        assert abs(E.imag) < 1e-10 and np.abs(H2 @ v - E.real * v).max() < 1e-8
        Tv = np.zeros_like(v)
        Tv[perm] = v
        assert np.abs(Tv - np.exp(-1j * K) * v).max() < 1e-8
        vecs.append(v)
    assert len(vecs) >= L                                                          # roots were found
    if L % 2 == 0:
        Phi = np.zeros(len(lab))
        for x in range(L):
            Phi[ring.index_of[(1 << x) | (1 << ((x + 1) % L))]] = (-1) ** x
        assert np.abs(H2 @ Phi - (Delta * L / 4 - Delta) * Phi).max() < 1e-12
        TPhi = np.zeros_like(Phi)
        TPhi[perm] = Phi
        assert np.abs(TPhi + Phi).max() < 1e-12
        assert max(abs(np.vdot(Phi, v)) for v in vecs) < 1e-8
        assert np.linalg.matrix_rank(np.stack(vecs), tol=1e-6) < math.comb(L, 2)    # never complete
    elif (L, Delta) == (9, 0.6):
        assert np.linalg.matrix_rank(np.stack(vecs), tol=1e-6) == 36              # tutorial's observation


def test_9_13_lanczos_tridiagonal():
    """Proposition 9.13: Q is orthonormal, Q H Q^* = T_m tridiagonal with beta_j > 0, and H Q_m = Q_m T_m + r e_m^T."""
    gen = torch.Generator().manual_seed(3)
    k, m = 50, 20
    A = torch.randn(k, k, dtype=c128, generator=gen)
    H = (A + A.mH) / 2
    v0 = torch.randn(1, k, dtype=c128, generator=gen)
    alpha, beta, Q = krylov.lanczos(lambda V: V @ H.T, v0, m)
    Qm = Q[0]                                                                       # rows q_j
    assert (Qm.conj() @ Qm.T - torch.eye(m, dtype=c128)).abs().max() < 1e-12
    T = Qm.conj() @ H @ Qm.T
    Tri = torch.diag(alpha[0]) + torch.diag(beta[0], 1) + torch.diag(beta[0], -1)
    assert (T - Tri.to(c128)).abs().max() < 1e-12 and (beta > 0).all()
    assert (Qm[0] - v0[0] / v0[0].norm()).abs().max() < 1e-14
    R = H @ Qm.T - Qm.T @ Tri.to(c128)                                              # only the last column
    assert R[:, :-1].abs().max() < 1e-12
    assert (Qm.conj() @ R[:, -1]).abs().max() < 1e-12                               # r_m orthogonal to K_m


def test_9_14_ritz_values_monotone():
    """Theorem 9.14: theta_j^(m) >= lambda_j and theta_j^(m) >= theta_j^(m+1) for all j <= m."""
    L, N = 10, 5
    xx = SpinChain(L, complex_valued=False)
    Hs = xx.xxz_sparse(J=1.0, Delta=0.6, h=0.3 * torch.randn(L, dtype=f64, generator=torch.Generator().manual_seed(1)),
                       sector=N)
    lam = torch.linalg.eigvalsh(Hs.to_dense()[0])
    alpha, beta, Q = krylov.lanczos(Hs.matvec, torch.randn(1, xx.sector_dim(N), dtype=f64,
                                                           generator=torch.Generator().manual_seed(2)), 30)
    T = torch.diag(alpha[0]) + torch.diag(beta[0], 1) + torch.diag(beta[0], -1)
    prev = None
    for m in range(1, 31):
        th = torch.linalg.eigvalsh(T[:m, :m])
        assert (th >= lam[:m] - 1e-10).all()
        assert (th <= lam[-m:] + 1e-10).all()                                        # the mirror statement
        if prev is not None:
            assert (th[:m - 1] <= prev + 1e-10).all()
        prev = th


def test_9_14_invariant_krylov_space_gives_exact_ground_energy():
    """Theorem 9.14 (last claim): if K_{m+1} = K_m and v meets the lambda_1-eigenspace, the lowest Ritz value is
    lambda_1; a start vector in an invariant subspace never leaves it (Exercise 3)."""
    k = 12
    d = torch.linspace(-3, 2, k, dtype=f64)[None]
    Hd = krylov.SparseSectorHamiltonian(d, None, f64)                                # diagonal H
    v = torch.zeros(1, k, dtype=f64)
    v[0, [0, 4, 9]] = torch.tensor([0.3, 1.0, -2.0], dtype=f64)
    E, psi, res = krylov.ground_state(Hd, m=6, v0=v)
    assert abs(E.item() - d[0, 0].item()) < 1e-12 and res.item() < 1e-10
    v[0, 0] = 0.0                                                                    # no lambda_1 component
    E, psi, res = krylov.ground_state(Hd, m=6, v0=v)
    assert abs(E.item() - d[0, 4].item()) < 1e-12 and psi[0, 0].abs().item() < 1e-14


@pytest.mark.parametrize('L', [6, 8])
def test_9_14_exercise_3_symmetric_start(L):
    """Exercise 3 (Chapter 9): from the all-ones vector of sector L/2, ground_state returns J L/4, the top level."""
    J = 1.3
    c = SpinChain(L, 'periodic', complex_valued=False)
    E, _, res = krylov.ground_state(c.xxz_sparse(J, 1.0, sector=L // 2), v0=torch.ones(1, math.comb(L, L // 2), dtype=f64))
    top = torch.linalg.eigvalsh(c.xxz(J, 1.0, sector=L // 2).matrix[0, 0])[-1].item()
    assert abs(E.item() - J * L / 4) < 1e-10 and abs(top - J * L / 4) < 1e-10 and res.item() < 1e-10
    assert abs(torch.linalg.eigvalsh(c.xxz(J, 1.0, sector=L // 2).matrix[0, 0])[0].item() - E.item()) > 1


@pytest.mark.parametrize('J', [1.0, 1.7])
def test_9_14_ring_of_four(J):
    """Section 9.5: the Heisenberg ring with L = 4 has the ground-state energy -2J exactly; the residual certifies it."""
    c = SpinChain(4, 'periodic', complex_valued=False)
    E, psi, res = krylov.ground_state(c.xxz_sparse(J, 1.0, sector=2), generator=torch.Generator().manual_seed(0))
    w = np.linalg.eigvalsh(xxz_dense(4, c.bonds, J, 1.0))
    assert abs(w[0] + 2 * J) < 1e-12 and abs(E.item() + 2 * J) < 1e-10
    assert np.abs(w - E.item()).min() <= res.item() + 1e-14                          # residual bound


def test_9_15_krylov_step_conserves_norm_and_energy():
    """Proposition 9.15: ||w_tau|| = ||v||, <w,Hw> = <v,Hv> for every m and tau; w_tau = e^{-i tau H} v if K_m is invariant."""
    L, N = 8, 4
    c = SpinChain(L, 'periodic')
    Hs = c.xxz_sparse(1.0, 0.8, h=torch.linspace(-0.5, 0.5, L, dtype=f64), sector=N)
    Hd = Hs.to_dense()[0]
    gen = torch.Generator().manual_seed(4)
    v = torch.randn(1, math.comb(L, N), dtype=c128, generator=gen) * 1.7
    for m, tau in ((2, 3.0), (5, 10.0), (12, 0.4)):
        w, _ = krylov.expm_multiply_step(Hs, v, tau, m)
        assert abs(w.norm() - v.norm()).item() < 1e-12
        e_w = (w.conj() * Hs.matvec(w)).sum().real
        e_v = (v.conj() * Hs.matvec(v)).sum().real
        assert abs(e_w - e_v).item() < 1e-10
    w, _ = krylov.expm_multiply_step(Hs, v, 2.5, math.comb(L, N))                   # K_m = whole space
    ref = torch.linalg.matrix_exp(-2.5j * Hd) @ v[0]
    assert (w[0] - ref).abs().max() < 1e-9
    out = krylov.evolve(Hs, v, [0.0, 0.7, 1.5])
    for t, p in zip([0.0, 0.7, 1.5], out):
        assert (p[0] - torch.linalg.matrix_exp(-1j * t * Hd) @ v[0]).abs().max() < 1e-8


def test_9_6_translation_breaking_by_fields():
    """Section 9.6: T H T^* - H = sum_i (h_{i-1} - h_i) s^z_i for the ring with fields h_i."""
    L = 6
    h = np.random.default_rng(2).normal(size=L)
    c = SpinChain(L, 'periodic')
    H = to_full(c, c.xxz(1.0, 1.0, h=torch.tensor(h))).real
    T = translation(L)
    D = sum((h[i - 1] - h[i]) * site_ops(L, {i: SZ}) for i in range(L))
    assert np.abs(T @ H @ T.T - H - D).max() < 1e-12


def test_9_6_poisson_mean_ratio():
    """Proposition 5.14(a), used in Section 9.6: independent exponential spacings give <r> = 2 log 2 - 1."""
    e = torch.distributions.Exponential(torch.tensor(1.0, dtype=f64))
    torch.manual_seed(11)
    levels = torch.cumsum(e.sample((400, 200)), dim=-1)
    r = SpinChain.level_spacing_ratio(levels)
    z = (r.mean().item() - (2 * math.log(2) - 1)) / (r.std().item() / math.sqrt(len(r)))
    assert abs(z) < 5


def test_9_exercise_1_momentum_dims():
    """Exercise 1 (Chapter 9): for L = 4, N = 2 the momentum sectors have dimensions 2, 1, 2, 1."""
    assert SpinChain(4, 'periodic').momentum_dims(2) == [2, 1, 2, 1]


# ==========================================================================================================
# Chapter 10. Quantum computing
# ==========================================================================================================
PX = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
PY = np.array([[0.0, -1j], [1j, 0.0]])
PZ = np.diag([1.0, -1.0]).astype(complex)


def qubit_op(op, j, n):
    """op on qubit j of n (qubit 0 = leftmost tensor factor)."""
    out = np.ones((1, 1), dtype=complex)
    for i in range(n):
        out = np.kron(out, op if i == j else np.eye(2))
    return out


def null_space_dim(Mat, tol=1e-9):
    s = np.linalg.svd(Mat, compute_uv=False)
    return int((s < tol * max(1.0, s.max())).sum()) + max(0, Mat.shape[1] - len(s))


def fixed_point_dim(kraus):
    d = kraus[0].shape[0]
    S = sum(np.kron(E, E.conj()) for E in kraus)                          # row-major vec(E X E^*)
    return null_space_dim(S - np.eye(d * d))


def commutant_dim(ops):
    d = ops[0].shape[0]
    Lm = np.concatenate([np.kron(A, np.eye(d)) - np.kron(np.eye(d), A.T) for A in ops])
    return null_space_dim(Lm)


def apply_kraus(kraus, rho):
    return sum(E @ rho @ E.conj().T for E in kraus)


def dual_kraus(kraus, T):
    return sum(E.conj().T @ T @ E for E in kraus)


@pytest.mark.parametrize('n', [1, 2, 3, 4])
def test_10_3_collective_dephasing(n):
    """Corollary 10.3 and Theorem 10.2: A = span{P_w}, Fix(Phi) = A' = Fix(Phi^*), of dimension sum_w binom(n,w)^2."""
    theta = 0.3
    Ztot = sum(qubit_op(PZ, j, n) for j in range(n))
    U = scipy.linalg.expm(-1j * theta * Ztot)
    kraus = [U / math.sqrt(2), U.conj().T / math.sqrt(2)]
    dim = sum(math.comb(n, w) ** 2 for w in range(n + 1))
    assert dim == math.comb(2 * n, n)                                          # Exercise 2
    assert fixed_point_dim(kraus) == commutant_dim(kraus + [E.conj().T for E in kraus]) == dim
    assert fixed_point_dim([E.conj().T for E in kraus]) == dim
    assert sum(k * k for k in SpinChain(n).algebra.k_factors) == dim
    # Proposition 10.1(b): every H_w is decoherence free
    weight = popcount(np.arange(2 ** n))
    rng = np.random.default_rng(n)
    for w in range(n + 1):
        idx = np.nonzero(weight == w)[0]
        G = np.zeros((2 ** n, 2 ** n), dtype=complex)
        G[np.ix_(idx, idx)] = rng.normal(size=(len(idx), len(idx))) + 1j * rng.normal(size=(len(idx), len(idx)))
        rho = G @ G.conj().T / np.trace(G @ G.conj().T)
        assert np.abs(apply_kraus(kraus, rho) - rho).max() < 1e-12


def test_10_3_library_channel_is_identity_on_sectors():
    """Corollary 10.3: on (+)_w M_binom(n,w) the collective dephasing has central Kraus operators and acts trivially."""
    n, theta = 5, 0.7
    M = SpinChain(n).algebra
    phases = [complex(np.exp(-1j * theta * (n - 2 * w))) for w in range(n + 1)]
    Phi = channels.kraus_channel([M.central([z / math.sqrt(2) for z in phases]),
                                  M.central([z.conjugate() / math.sqrt(2) for z in phases])])
    rho = states.random_density_matrix(M, batch_size=4)
    assert (Phi(rho).matrix - rho.matrix).abs().max() < 1e-12 and Phi.is_trace_preserving() and Phi.is_unital()


def test_10_2_fixed_points_of_random_unital_channels():
    """Theorem 10.2: for a unital channel Fix(Phi) = A' = Fix(Phi^*); without unitality only A' c Fix(Phi^*) holds."""
    rng = np.random.default_rng(8)
    d = 6
    # a random mixed-unitary channel with unitaries in M_2 (x) 1_3: interaction algebra M_2 (x) 1_3, A' = 1_2 (x) M_3
    us = [np.kron(scipy.linalg.expm(1j * (lambda a: a + a.conj().T)(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))),
                  np.eye(3)) for _ in range(3)]
    p = rng.dirichlet(np.ones(3))
    kraus = [math.sqrt(pi) * u for pi, u in zip(p, us)]
    assert np.abs(sum(E @ E.conj().T for E in kraus) - np.eye(d)).max() < 1e-12
    assert fixed_point_dim(kraus) == commutant_dim(kraus + [E.conj().T for E in kraus]) == 9
    assert fixed_point_dim([E.conj().T for E in kraus]) == 9
    T = np.kron(np.eye(2), rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))      # an element of A'
    assert np.abs(apply_kraus(kraus, T) - T).max() < 1e-12 and np.abs(dual_kraus(kraus, T) - T).max() < 1e-12


@pytest.mark.parametrize('g', [0.1, 0.5, 0.9])
def test_10_exercise_1_amplitude_damping(g):
    """Exercise 1 (Chapter 10): for amplitude damping A' = C1, Fix(Phi) = C|0><0|, Fix(Phi^*) = C1."""
    E = [np.array([[1, 0], [0, math.sqrt(1 - g)]], dtype=complex), np.array([[0, math.sqrt(g)], [0, 0]], dtype=complex)]
    assert commutant_dim(E + [e.conj().T for e in E]) == 1
    assert fixed_point_dim(E) == 1 and fixed_point_dim([e.conj().T for e in E]) == 1
    P0 = np.diag([1.0, 0.0]).astype(complex)
    assert np.abs(apply_kraus(E, P0) - P0).max() < 1e-15
    assert np.abs(dual_kraus(E, np.eye(2)) - np.eye(2)).max() < 1e-15
    T = np.array([[0.3, 1.2], [-0.7j, 2.0]])
    assert np.abs(apply_kraus(E, T) - np.array([[0.3 + g * 2.0, math.sqrt(1 - g) * 1.2],
                                                [math.sqrt(1 - g) * -0.7j, (1 - g) * 2.0]])).max() < 1e-14


@pytest.mark.parametrize('n', [2, 3, 4])
def test_10_exercise_2_parity_dephasing(n):
    """Exercise 2 (Chapter 10): for theta = pi/2, A = span{P_even, P_odd} and dim Fix(Phi) = 2 (2^{n-1})^2."""
    U = scipy.linalg.expm(-1j * (math.pi / 2) * sum(qubit_op(PZ, j, n) for j in range(n)))
    kraus = [U / math.sqrt(2), U.conj().T / math.sqrt(2)]
    parity = np.diag((-1.0) ** popcount(np.arange(2 ** n)))
    assert np.abs(U - (-1j) ** n * parity).max() < 1e-12
    assert fixed_point_dim(kraus) == 2 * 4 ** (n - 1)


def test_10_4_collective_rotations():
    """Proposition 10.4 and Section 10.1: for n = 4, dim B = 35, dim B' = 14, dim K = 2; collective channels fix
    the states supported in K and Phi^*(T) = T for T in B' (Proposition 10.1(a))."""
    n = 4
    S = [sum(qubit_op(P, j, n) for j in range(n)) / 2 for P in (PX, PY, PZ)]
    assert commutant_dim(S) == 14
    # dim B: B = (B')' with a basis of B' from the null space
    d = 2 ** n
    Lm = np.concatenate([np.kron(A, np.eye(d)) - np.kron(np.eye(d), A.T) for A in S])
    _, s, Vh = np.linalg.svd(Lm)
    basis = [v.reshape(d, d) for v in Vh[-14:].conj()]
    assert s[-14:].max() < 1e-9 and commutant_dim(basis) == 35
    K = scipy.linalg.null_space(np.concatenate(S))
    assert K.shape[1] == 2
    rng = np.random.default_rng(4)
    axes = rng.normal(size=(3, 3))
    axes /= np.linalg.norm(axes, axis=1, keepdims=True)
    kraus = [scipy.linalg.expm(-1j * t * sum(a * Sa for a, Sa in zip(ax, S))) / math.sqrt(3)
             for t, ax in zip(rng.uniform(0, np.pi, 3), axes)]
    G = K @ (rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))
    rho = G @ G.conj().T / np.trace(G @ G.conj().T)
    assert np.abs(apply_kraus(kraus, rho) - rho).max() < 1e-12
    for T in basis:
        assert np.abs(dual_kraus(kraus, T) - T).max() < 1e-12
    G = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    rho = G @ G.conj().T / np.trace(G @ G.conj().T)
    assert np.abs(apply_kraus(kraus, rho) - rho).max() > 1e-3                    # a generic state moves


def _knill_laflamme_matrix(kraus, P, tol=1e-12):
    lam = np.array([[np.trace(P @ a.conj().T @ b @ P) / np.trace(P) for b in kraus] for a in kraus])
    ok = all(np.abs(P @ a.conj().T @ b @ P - lam[i, j] * P).max() < tol
             for i, a in enumerate(kraus) for j, b in enumerate(kraus))
    return lam if ok else None


def _kl_recovery(kraus, P):
    """The recovery of the proof of Theorem 10.5 (sufficiency), as Kraus operators."""
    lam = _knill_laflamme_matrix(kraus, P)
    d, U = np.linalg.eigh(lam)
    rec, Q = [], np.zeros_like(P)
    for k in range(len(kraus)):
        F = sum(U[i, k] * kraus[i] for i in range(len(kraus)))
        if d[k] > 1e-12:
            V = F @ P / math.sqrt(d[k])
            rec.append(V.conj().T)
            Q = Q + V @ V.conj().T
        else:
            assert np.abs(F @ P).max() < 1e-10                                    # d_k = 0 => F_k P = 0
    rec.append(np.eye(P.shape[0]) - Q)
    return lam, rec


def test_10_5_knill_laflamme_random_correctable_channel():
    """Theorem 10.5: for a code C = V(C^2 (x) e_0) and errors V(1 (x) u_a)V^*, the condition holds with a non-diagonal
    Lambda (Hermitian, positive, trace 1), and the recovery of the proof corrects every code state."""
    rng = np.random.default_rng(12)
    V = np.linalg.qr(rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8)))[0]
    W = np.linalg.qr(rng.normal(size=(12, 4)) + 1j * rng.normal(size=(12, 4)))[0]       # isometry C^4 -> C^3 (x) C^4
    us = [W[4 * a:4 * a + 4] for a in range(3)]                                   # Kraus of a channel on C^4
    kraus = [V @ np.kron(np.eye(2), u) @ V.conj().T for u in us]
    assert np.abs(sum(E.conj().T @ E for E in kraus) - np.eye(8)).max() < 1e-12
    e0 = np.zeros((4, 4))
    e0[0, 0] = 1
    P = V @ np.kron(np.eye(2), e0) @ V.conj().T
    lam, rec = _kl_recovery(kraus, P)
    assert np.abs(lam - lam.conj().T).max() < 1e-12 and np.linalg.eigvalsh(lam).min() > -1e-12
    assert abs(np.trace(lam) - 1) < 1e-12 and np.abs(lam - np.diag(np.diag(lam))).max() > 1e-3
    assert np.abs(sum(R.conj().T @ R for R in rec) - np.eye(8)).max() < 1e-12      # R is a channel
    M8 = TypeIAlgebra([8], [8], precision='double')
    op8 = lambda X: M8.operator(torch.tensor(X, dtype=c128).reshape(1, 1, 8, 8))
    RPhi = channels.kraus_channel([op8(R) for R in rec]) @ channels.kraus_channel([op8(E) for E in kraus])
    for _ in range(3):
        G = P @ (rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))) @ P
        rho = G @ G.conj().T / np.trace(G @ G.conj().T)
        out = RPhi(DensityMatrix(M8, matrix=torch.tensor(rho, dtype=c128).reshape(1, 1, 8, 8)))
        assert np.abs(out.matrix[0, 0].numpy() - rho).max() < 1e-12


def _bitflip():
    ket = lambda s: np.eye(8)[int(s, 2)].astype(complex)
    P = np.outer(ket('000'), ket('000')) + np.outer(ket('111'), ket('111'))
    flips = [np.eye(8, dtype=complex)] + [qubit_op(PX, j, 3) for j in range(3)]
    return ket, P, flips


def _independent(op, prob):
    return [math.sqrt(prob ** sum(s) * (1 - prob) ** (3 - sum(s)))
            * np.linalg.multi_dot([qubit_op(op, j, 3) if s[j] else np.eye(8) for j in range(3)] + [np.eye(8)])
            for s in itertools.product([0, 1], repeat=3)]


@pytest.mark.parametrize('p', [0.02, 0.1, 0.3])
def test_10_5_bit_flip_code(p):
    """Section 10.2: Lambda = diag(p_0..p_3) for at most one flip, R = sum_k P X_k . X_k P corrects it, and for
    independent flips the condition fails and the fidelity of |0_L> is 1 - 3p^2 + 2p^3."""
    ket, P, flips = _bitflip()
    pk = [1 - 3 * p, p, p, p] if p < 1 / 3 else [0.1, 0.3, 0.3, 0.3]
    E_one = [math.sqrt(q) * X for q, X in zip(pk, flips)]
    assert np.abs(_knill_laflamme_matrix(E_one, P) - np.diag(pk)).max() < 1e-12
    R = [P @ X for X in flips]
    assert np.abs(sum(r.conj().T @ r for r in R) - np.eye(8)).max() < 1e-12          # Q = 1
    a, b = 0.6, 0.8j
    psi = a * ket('000') + b * ket('111')
    rho = np.outer(psi, psi.conj())
    assert np.abs(apply_kraus(R, apply_kraus(E_one, rho)) - rho).max() < 1e-12
    E_ind = _independent(PX, p)
    assert _knill_laflamme_matrix(E_ind, P) is None
    z = np.outer(ket('000'), ket('000'))
    F = np.trace(z @ apply_kraus(R, apply_kraus(E_ind, z))).real
    assert abs(F - (1 - 3 * p ** 2 + 2 * p ** 3)) < 1e-12


@pytest.mark.parametrize('phi', [0.2, 1.0, 2.5])
def test_10_exercise_3_linearity_of_correction(phi):
    """Exercise 3 (Chapter 10): R corrects every channel with Kraus operators in span{1, X_1, X_2, X_3}."""
    ket, P, flips = _bitflip()
    R = [P @ X for X in flips]
    rng = np.random.default_rng(int(10 * phi))
    G1 = scipy.linalg.expm(-1j * phi * flips[1])
    assert np.abs(G1 - (math.cos(phi) * flips[0] - 1j * math.sin(phi) * flips[1])).max() < 1e-12
    q = rng.dirichlet(np.ones(3))
    kraus = [math.sqrt(q[0]) * G1, math.sqrt(q[1]) * scipy.linalg.expm(-1j * 0.7 * flips[2]), math.sqrt(q[2]) * flips[3]]
    for _ in range(3):
        v = P @ (rng.normal(size=8) + 1j * rng.normal(size=8))
        rho = np.outer(v, v.conj()) / np.vdot(v, v)
        assert np.abs(apply_kraus(R, apply_kraus(kraus, rho)) - rho).max() < 1e-12


@pytest.mark.parametrize('q', [0.05, 0.2])
def test_10_6_sector_label_correction(q):
    """Proposition 10.6 and its example: condition 2 holds, majority vote reads the bit, P (R o Phi)^*(P_d) P = P_d for
    the measure-and-prepare recovery, and |+_L> has fidelity (1 + (1-2q)^3)/2 after the recovery of Section 10.2."""
    ket, P, flips = _bitflip()
    E_one = [math.sqrt(pk) * X for pk, X in zip((0.7, 0.1, 0.1, 0.1), flips)]
    kraus = [a @ b for a in E_one for b in _independent(PZ, q)]
    assert np.abs(sum(E.conj().T @ E for E in kraus) - np.eye(8)).max() < 1e-12
    P0, P1 = np.outer(ket('000'), ket('000')), np.outer(ket('111'), ket('111'))
    assert _knill_laflamme_matrix(kraus, P) is None
    assert max(np.abs(P0 @ e.conj().T @ f @ P1).max() for e in kraus for f in kraus) == 0
    N0 = np.diag([float(bin(i).count('1') <= 1) for i in range(8)]).astype(complex)
    Ns, Ps = [N0, np.eye(8) - N0], [P0, P1]
    # the measure-and-prepare recovery R(s) = sum_c Tr(N_c s) omega_c, omega_c = P_c (rank one here)
    rec = [np.outer(ket(c * 3), f) @ Nc for c, Nc in zip('01', Ns) for f in np.eye(8)]
    assert np.abs(sum(r.conj().T @ r for r in rec) - np.eye(8)).max() < 1e-12
    for Nd, Pd in zip(Ns, Ps):
        assert np.abs(P @ dual_kraus(kraus, dual_kraus(rec, Pd)) @ P - Pd).max() < 1e-12
        assert np.abs(P @ dual_kraus(kraus, Nd) @ P - Pd).max() < 1e-12
    R = [P @ X for X in flips]
    plus = (ket('000') + ket('111')) / math.sqrt(2)
    rho = np.outer(plus, plus.conj())
    F = (plus.conj() @ apply_kraus(R, apply_kraus(kraus, rho)) @ plus).real
    assert abs(F - (1 + (1 - 2 * q) ** 3) / 2) < 1e-12


def test_10_8_xy_mixer():
    """Proposition 10.8: X(x)X + Y(x)Y = 2(|01><10| + |10><01|), the XY mixer is SpinChain.xxz(J=1, Delta=0), and
    the ansatz stays in its Hamming-weight sector."""
    XX_YY = np.kron(PX, PX) + np.kron(PY, PY)
    ref = np.zeros((4, 4))
    ref[1, 2] = ref[2, 1] = 2
    assert np.abs(XX_YY - ref).max() < 1e-15
    n = 5
    chain = SpinChain(n, 'periodic')
    HM = sum(qubit_op(PX, i, n) @ qubit_op(PX, j, n) + qubit_op(PY, i, n) @ qubit_op(PY, j, n) for i, j in chain.bonds) / 4
    # qubit j = leftmost factor = bit n-1-j of the index; reverse the bits to compare with the label basis
    x = np.arange(2 ** n)
    rev = np.array([int(format(v, f'0{n}b')[::-1], 2) for v in x])
    assert np.abs(to_full(chain, chain.xxz(1.0, 0.0))[np.ix_(rev, rev)] - HM).max() < 1e-12
    H_M = chain.xxz(1.0, 0.0)
    H_C = chain.algebra.zero()
    for i, j in chain.bonds + [(0, 2)]:
        H_C = H_C + (chain.algebra.identity() - chain.sz(i) @ chain.sz(j) * 4.0) * 0.5
    rho = chain.basis_state('10100')
    for gamma, beta in ((0.4, 0.9), (1.3, 0.2)):
        U = dynamics.propagator(H_M, beta) @ dynamics.propagator(H_C, gamma)
        rho = channels.unitary_channel(U)(rho)
    p = rho.sector_probabilities()[0]
    assert abs(p[2].item() - 1) < 1e-12 and abs(p.sum().item() - 1) < 1e-12


@pytest.mark.parametrize('g,p_layers', [(0.05, 2), (0.3, 3)])
def test_10_9_leakage_and_post_selection(g, p_layers):
    """Proposition 10.9: weights stay <= w, the weight w survives with probability (1-g)^{wp}, and the post-selected
    state is the noiseless one."""
    n, w = 5, 2
    chain = SpinChain(n, 'periodic')
    M = chain.algebra
    torch.manual_seed(3)
    layers = [M.random_unitary_operator(1) for _ in range(p_layers)]
    noisy = ideal = chain.basis_state('01001')
    damping = [chain.site_amplitude_damping(i, g) for i in range(n)]
    for U in layers:
        ideal = channels.unitary_channel(U)(ideal)
        noisy = channels.unitary_channel(U)(noisy)
        for D in damping:
            noisy = D(noisy)
    probs = noisy.sector_probabilities()[0]
    assert probs[w + 1:].abs().max() < 1e-14 and abs(probs.sum().item() - 1) < 1e-12
    P_w = M.central([1.0 if v == w else 0.0 for v in range(n + 1)])
    post, accept = noisy.condition_on(P_w)
    assert abs(accept.item() - (1 - g) ** (w * p_layers)) < 1e-12
    assert (post.matrix - ideal.matrix).abs().max() < 1e-12
    # one damping step on a block-diagonal state: P_w D(sigma) P_w = (1-g)^w P_w sigma P_w (proof)
    sigma = states.random_density_matrix(M, batch_size=1)
    out = sigma
    for D in damping:
        out = D(out)
    assert (out.matrix[0, n] - (1 - g) ** n * sigma.matrix[0, n]).abs().max() < 1e-12


def test_10_10_lindblad_post_selection():
    """Proposition 10.10: under loss (rate g1) and dephasing (rate g2), P_w rho(t) P_w = e^{-g1 w t} sigma(t) with sigma
    independent of g1; for g2 = 0, sigma(t) = e^{-iHt} rho(0) e^{iHt}."""
    n, w, T = 4, 2, 0.6
    chain = SpinChain(n, 'periodic')
    H = chain.xxz(1.0, 0.3)
    jumps = [chain.lowering(i) for i in range(n)] + [chain.sz(i) for i in range(n)]
    rho0 = chain.basis_state('1010')
    P_w = chain.algebra.central([1.0 if v == w else 0.0 for v in range(n + 1)])
    ideal = dynamics.von_neumann(rho0, H, [T])[0]
    post = {}
    for g1, g2 in ((0.3, 0.0), (0.0, 0.2), (0.3, 0.2)):
        rho_T = dynamics.lindblad_evolve(rho0, H, jumps, [0.0, T], rates=[g1] * n + [g2] * n, substeps=80)[-1]
        assert rho_T.sector_probabilities()[0, w + 1:].abs().max() < 1e-12
        post[g1, g2], acc = rho_T.condition_on(P_w)
        assert abs(acc.item() - math.exp(-g1 * w * T)) < 1e-9
    assert (post[0.3, 0.0].matrix - ideal.matrix).abs().max() < 1e-8
    assert (post[0.3, 0.2].matrix - post[0.0, 0.2].matrix).abs().max() < 1e-8
    assert abs(post[0.0, 0.2].trace.real.item() - 1) < 1e-10


# ==========================================================================================================
# Chapter 11. Optimisation and trace inequalities
# ==========================================================================================================
ell = lambda delta: (1 - delta) / (1 + delta)


def _blockdiag(blocks):
    return scipy.linalg.block_diag(*blocks)


def _haar(k, rng):
    Z = (rng.normal(size=(k, k)) + 1j * rng.normal(size=(k, k))) / math.sqrt(2)
    Q, R = np.linalg.qr(Z)
    return Q * (np.diag(R) / np.abs(np.diag(R)))


def _psd(spec, rng):
    U = _haar(len(spec), rng)
    return U @ np.diag(spec) @ U.conj().T


def _nuc(A):
    return np.linalg.svd(A, compute_uv=False).sum()


def test_11_ell_is_a_decreasing_involution():
    """Section 11 (notation): ell is a strictly decreasing involution of [0,1], and contrast(lambda_max=1) = ell(lambda_min)."""
    t = np.linspace(0, 1, 101)
    assert np.abs(ell(ell(t)) - t).max() < 1e-15 and (np.diff(ell(t)) < 0).all()
    assert np.abs(1 - ell(t[:-1]) - 2 * t[:-1] / (1 + t[:-1])).max() < 1e-15
    alg = TypeIAlgebra([2, 3], [2, 3], precision='double')
    lo = torch.tensor([0.1, 0.5, 1.0], dtype=f64)
    X = alg.operator_from_eigenvalues(lambda k: torch.cat([torch.ones(3, 1, dtype=f64), lo[:, None].expand(3, k - 1)], 1),
                                      batch_size=3, force_positive=True)
    assert (X.michelson_contrast - ell(lo)).abs().max() < 1e-12


@pytest.mark.parametrize('k', [2, 3, 5])
def test_11_3_rearrangement(k):
    """Lemma 11.3: sum_i p_i q_{k+1-i} <= sum_ij p_i q_j c_ij <= sum_i p_i q_i for doubly stochastic c."""
    rng = np.random.default_rng(k)
    perms = list(itertools.permutations(range(k)))
    for trial in range(50):
        p, q = -np.sort(-rng.normal(size=k)), -np.sort(-rng.normal(size=k))
        if trial % 2:
            c = np.abs(_haar(k, rng)) ** 2                                       # unistochastic
        else:
            wts = rng.dirichlet(np.ones(4))
            c = sum(wt * np.eye(k)[list(perms[rng.integers(len(perms))])] for wt in wts)   # Birkhoff
        assert np.abs(c.sum(0) - 1).max() < 1e-12 and np.abs(c.sum(1) - 1).max() < 1e-12
        val = p @ c @ q
        assert p @ q[::-1] - 1e-12 <= val <= p @ q + 1e-12


@pytest.mark.parametrize('k', [2, 4, 6])
def test_11_4_von_neumann_trace_inequality(k):
    """Theorem 11.4: |Tr(AB)| <= sum_i s_i(A) s_i(B), with equality for A = W S_A V^*, B = V S_B W^*."""
    gen = torch.Generator().manual_seed(k)
    Mk = TypeIAlgebra([k], [k], precision='double')
    A = Mk.operator(torch.randn(500, 1, k, k, dtype=c128, generator=gen))
    B = Mk.operator(torch.randn(500, 1, k, k, dtype=c128, generator=gen))
    sA, sB = torch.linalg.svdvals(A.matrix[:, 0]), torch.linalg.svdvals(B.matrix[:, 0])
    assert ((A @ B).trace.abs() <= (sA * sB).sum(-1) + 1e-12).all()
    rng = np.random.default_rng(k)
    W, V = _haar(k, rng), _haar(k, rng)
    a, b = np.sort(rng.random(k))[::-1], np.sort(rng.random(k))[::-1]
    A1, B1 = W @ np.diag(a) @ V.conj().T, V @ np.diag(b) @ W.conj().T
    assert abs(np.trace(A1 @ B1) - a @ b) < 1e-12


def test_11_5_spectral_trace_bounds_blockwise():
    """Corollary 11.5: sum p^dn q^up <= Tr(PQ) <= sum p^dn q^dn per block, and Tr(PQ) <= ||PQ||_1 <= sum p^dn q^dn for
    P, Q >= 0, on M_1 (+) M_3 (+) M_4."""
    alg = TypeIAlgebra([1, 3, 4], [1, 3, 4], precision='double')
    torch.manual_seed(5)
    P, Q = (alg.operator_from_eigenvalues(lambda k: torch.randn(300, k, dtype=f64), batch_size=300,
                                          force_self_adjoint=True) for _ in range(2))
    t = (P @ Q).trace.real
    lo = sum((p * q.flip(-1)).sum(-1) for p, q in zip(P.eigenvalues(), Q.eigenvalues()))   # increasing order
    hi = sum((p * q).sum(-1) for p, q in zip(P.eigenvalues(), Q.eigenvalues()))
    assert (t >= lo - 1e-12).all() and (t <= hi + 1e-12).all()
    Pp, Qp = P @ P, Q @ Q
    nuc = (Pp @ Qp).trace_norm()
    hi2 = sum((torch.sort(p ** 2)[0] * torch.sort(q ** 2)[0]).sum(-1) for p, q in zip(P.eigenvalues(), Q.eigenvalues()))
    assert ((Pp @ Qp).trace.real <= nuc + 1e-12).all() and (nuc <= hi2 + 1e-10).all()


@pytest.mark.parametrize('ks', [[3], [1, 3, 2], [4, 4]])
def test_11_6_exp2_fixed_spectra(ks):
    """Proposition 11.6: z = ||X_0 V Y_0||_1 - Tr(X_0 Y_0) lies in [-s, s], s = sum x^dn (y^dn - y^up), both attained."""
    rng = np.random.default_rng(len(ks) + sum(ks))
    xs = [np.sort(rng.uniform(0.2, 1, k))[::-1] for k in ks]
    ys = [np.sort(rng.uniform(0.1, 1, k))[::-1] for k in ks]
    s = sum(x @ (y - y[::-1]) for x, y in zip(xs, ys))
    for _ in range(200):
        X0 = _blockdiag([_psd(x, rng) for x in xs])
        Y0 = _blockdiag([_psd(y, rng) for y in ys])
        V = _blockdiag([_haar(k, rng) for k in ks])
        U = _blockdiag([_haar(k, rng) for k in ks])
        X, Y = U @ X0, V @ Y0
        absX = scipy.linalg.sqrtm(X.conj().T @ X)
        assert np.abs(absX - X0).max() < 1e-8                                         # |UX_0| = X_0
        z = _nuc(X @ Y) - np.trace(X0 @ Y0).real
        assert abs(z - (_nuc(X0 @ V @ Y0) - np.trace(X0 @ Y0).real)) < 1e-10
        assert -s - 1e-10 <= z <= s + 1e-10
    J = _blockdiag([np.eye(k)[::-1] for k in ks])
    X0 = _blockdiag([np.diag(x) for x in xs])
    z_max = _nuc(X0 @ J @ _blockdiag([np.diag(y[::-1]) for y in ys])) - sum(x @ y[::-1] for x, y in zip(xs, ys))
    z_min = _nuc(X0 @ J @ _blockdiag([np.diag(y) for y in ys])) - sum(x @ y for x, y in zip(xs, ys))
    assert abs(z_max - s) < 1e-12 and abs(z_min + s) < 1e-12


@pytest.mark.parametrize('ks', [[3, 3], [2, 4, 1], [5], [1, 2]])
def test_11_7_exp2_envelope(ks):
    """Corollary 11.7: sup z = -inf z = sum_c floor(k_c/2)(1 - l_X)(1 - l_Y), attained by two-level spectra; random
    admissible points of the library parametrisation stay inside."""
    dX, dY = 0.5, 0.8
    lX, lY = ell(dX), ell(dY)
    env = sum(k // 2 for k in ks) * (1 - lX) * (1 - lY)

    def spec(k, lo):
        mid = [0.6] if k % 2 else []
        return np.array([1.0] * (k // 2) + mid + [lo] * (k // 2))
    X0 = _blockdiag([np.diag(spec(k, lX)) for k in ks])
    J = _blockdiag([np.eye(k)[::-1] for k in ks])
    Yup = _blockdiag([np.diag(spec(k, lY)[::-1]) for k in ks])
    Ydn = _blockdiag([np.diag(spec(k, lY)) for k in ks])
    assert abs(_nuc(X0 @ J @ Yup) - np.trace(X0 @ Yup).real - env) < 1e-12
    assert abs(_nuc(X0 @ J @ Ydn) - np.trace(X0 @ Ydn).real + env) < 1e-12
    alg = TypeIAlgebra([3 * k for k in ks], ks, precision='double')              # m_c = 3 plays no role
    torch.manual_seed(sum(ks))
    B = 200
    X = opt.PositiveParam(alg, B, delta=torch.full((B,), dX, dtype=f64)).operator()
    Y = opt.PositiveParam(alg, B, delta=torch.full((B,), dY, dtype=f64)).operator()
    V = alg.random_unitary_operator(B)
    M = alg.operator((X @ V @ Y).matrix, is_self_adjoint=False)
    z = M.trace_norm() - (X @ Y).trace.real
    assert (z.abs() <= env + 1e-10).all()


def test_11_8_trace_norm_of_difference():
    """Lemma 11.8: ||A - B||_1 <= Tr A + Tr B for A, B >= 0, with equality iff AB = 0."""
    rng = np.random.default_rng(9)
    n = 6
    for _ in range(50):
        G, H = rng.normal(size=(n, 3)) + 1j * rng.normal(size=(n, 3)), rng.normal(size=(n, 2)) + 1j * rng.normal(size=(n, 2))
        A, B = G @ G.conj().T, H @ H.conj().T
        gap = np.trace(A + B).real - _nuc(A - B)
        assert gap > 1e-8                                                            # AB != 0 generically
        Qb = np.linalg.qr(np.concatenate([G, H], 1))[0]
        A2 = Qb[:, :3] @ np.diag(rng.random(3)) @ Qb[:, :3].conj().T              # orthogonal supports
        B2 = Qb[:, 3:] @ np.diag(rng.random(2)) @ Qb[:, 3:].conj().T
        assert np.abs(A2 @ B2).max() < 1e-12 and abs(np.trace(A2 + B2).real - _nuc(A2 - B2)) < 1e-12


def test_11_9_exp3_inequality_and_equality():
    """Proposition 11.9: Tr|YXY| <= Tr(Y|X|Y), with equality iff Y X_+ Y^2 X_- Y = 0 (commuting pairs; the example
    with real a in (-1, 1))."""
    alg = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')
    torch.manual_seed(6)
    Xs = alg.operator_from_eigenvalues(lambda k: torch.randn(400, k, dtype=f64), batch_size=400, force_self_adjoint=True)
    Yp = alg.operator_from_eigenvalues(lambda k: torch.rand(400, k, dtype=f64), batch_size=400, force_positive=True)
    z = (Yp @ Xs.abs() @ Yp).trace.real - (Yp @ Xs @ Yp).trace_norm()
    assert (z >= -1e-12).all()
    Xp, Xm = (Xs + Xs.abs()) * 0.5, (Xs.abs() - Xs) * 0.5                         # X_+, X_-
    cond = (Yp @ Xp @ Yp @ Yp @ Xm @ Yp).frobenius_norm()
    zero = cond < 1e-10
    assert 0 < zero.sum() < len(zero) and torch.equal(z < 1e-10, zero)            # equality iff the condition
    U = alg.random_unitary_operator(1)                                           # commuting pairs: z = 0
    diag = lambda k: torch.diag_embed(torch.randn(1, k, dtype=f64)).to(c128)
    Xc = U @ alg.from_blocks([diag(k) for k in (1, 2, 3)]) @ U.adjoint()
    Yc = U @ alg.from_blocks([diag(k).abs() for k in (1, 2, 3)]) @ U.adjoint()
    Xc, Yc = alg.operator(Xc.matrix, is_self_adjoint=True), alg.operator(Yc.matrix, is_self_adjoint=True)
    assert abs((Yc @ Xc.abs() @ Yc).trace.real - (Yc @ Xc @ Yc).trace_norm()).item() < 1e-12
    X = np.diag([1.0, -1.0, 0.0])
    for a in (-0.9, 0.3, 0.7):
        Y = np.array([[1, 0, a], [0, 1, 0], [a, 0, 1]])
        assert np.linalg.eigvalsh(Y).min() > 0 and np.abs(X @ Y - Y @ X).max() > 0.1
        Xp, Xm = np.diag([1.0, 0, 0]), np.diag([0, 1.0, 0])
        assert np.abs(Y @ Xp @ Y @ Y @ Xm @ Y).max() < 1e-15
        assert abs(np.trace(Y @ np.abs(X) @ Y) - _nuc(Y @ X @ Y)) < 1e-12
    Y = np.array([[1, 0.3, 0], [0.3, 1, 0], [0, 0, 1]])                         # Y X_+ Y^2 X_- Y != 0: strict
    assert np.trace(Y @ np.abs(X) @ Y) - _nuc(Y @ X @ Y) > 1e-3


def test_11_1_exponential_chart_is_onto():
    """Lemma 11.1: U = U_0 e^{G - G^*} with G = log(U_0^* U)/2 for every unitary U of M_1 (+) M_2 (+) M_3."""
    alg = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')
    torch.manual_seed(7)
    P = opt.UnitaryParam(alg, 2)
    target = alg.random_unitary_operator(2).matrix
    V = P.base.mH @ target
    S = torch.zeros_like(V)
    for b in range(2):
        for c, k in enumerate(alg.k_factors):
            S[b, c, :k, :k] = torch.from_numpy(scipy.linalg.logm(V[b, c, :k, :k].numpy()))
    assert (S + S.mH).abs().max() < 1e-10
    with torch.no_grad():
        P.G.copy_(S / 2)
    assert (P.matrix() - target).abs().max() < 1e-10
    proj = torch.zeros(1, 1, 3, 3, dtype=c128)                                        # e^{2 pi i P} = 1 = e^0
    proj[0, 0, :2, :2] = torch.tensor([[0.5, 0.5], [0.5, 0.5]], dtype=c128)
    M3 = TypeIAlgebra([3], [3], precision='double')
    assert (torch.linalg.matrix_exp(2j * math.pi * proj[0, 0]) - torch.eye(3, dtype=c128)).abs().max() < 1e-12
    assert M3.operator(proj).is_projection


def test_11_1_charts_have_the_same_differential():
    """Section 11.1: exp, Cayley (S -> C(-S/2)) and polar retraction agree with U_0(1 + tS) to second order."""
    rng = np.random.default_rng(3)
    k = 4
    U0 = _haar(k, rng)
    A = rng.normal(size=(k, k)) + 1j * rng.normal(size=(k, k))
    S = A - A.conj().T
    one = np.eye(k)

    def charts(t):
        e = U0 @ scipy.linalg.expm(t * S)
        cay = U0 @ (one + t * S / 2) @ np.linalg.inv(one - t * S / 2)
        Apol = U0 @ (one + t * S)
        w, V = np.linalg.eigh(Apol.conj().T @ Apol)
        pol = Apol @ V @ np.diag(w ** -0.5) @ V.conj().T
        return e, cay, pol
    for chart in range(3):
        errs = [np.abs(charts(t)[chart] - U0 @ (one + t * S)).max() for t in (1e-3, 5e-4)]
        assert 3.5 < errs[0] / errs[1] < 4.5                                     # O(t^2)
        U = charts(0.7)[chart]
        assert np.abs(U.conj().T @ U - one).max() < 1e-12


@pytest.mark.parametrize('ks', [[1, 2], [3], [2, 2, 1]])
def test_11_2_positive_param_constraint_and_slots(ks):
    """Proposition 11.2: X(theta, G) lies in the contrast shell (lambda_max = 1, lambda_min = l(delta), spectrum in
    [l, 1]) with 1 and l in two distinct slots, every ordered pair of slots drawn with probability 1/(D(D-1))."""
    alg = TypeIAlgebra(ks, ks, precision='double')
    D = sum(ks)
    torch.manual_seed(D + len(ks))
    B = 6000
    delta = torch.rand(B, dtype=f64) * 0.95
    P = opt.PositiveParam(alg, B, delta=delta)
    with torch.no_grad():
        P.spectrum.theta.mul_(10.0)                                                  # far from the start
    X = P.operator()
    lo = ell(delta)
    lam = torch.cat([e for e in X.eigenvalues()], dim=1)
    assert (X.lambda_max - 1).abs().max() < 1e-12 and (X.lambda_min - lo).abs().max() < 1e-12
    assert (X.michelson_contrast - delta).abs().max() < 1e-12 and (lam >= lo[:, None] - 1e-12).all()
    slots = [(c, i) for c, k in enumerate(ks) for i in range(k)]
    ip = [slots.index(tuple(s)) for s in P.spectrum.pin_max.tolist()]
    im = [slots.index(tuple(s)) for s in P.spectrum.pin_min.tolist()]
    assert all(a != b for a, b in zip(ip, im))
    counts = np.zeros((D, D))
    np.add.at(counts, (ip, im), 1)
    p = 1 / (D * (D - 1))
    off = counts[~np.eye(D, dtype=bool)]
    z = (off - B * p) / math.sqrt(B * p * (1 - p))
    assert np.abs(z).max() < 4.5


def test_11_2_self_adjoint_param():
    """Section 11.1: X = U diag(sigma lambda) U^* has |X| = U diag(lambda) U^* in the contrast shell; the number of minus
    signs per block is uniform in {0, ..., k_c}."""
    alg = TypeIAlgebra([3, 1], [3, 1], precision='double')
    torch.manual_seed(13)
    B = 8000
    P = opt.SelfAdjointParam(alg, B, delta=torch.full((B,), 0.6, dtype=f64))
    X = P.operator()
    assert (X.abs().matrix - P.abs_operator().matrix).abs().max() < 1e-10
    assert (P.abs_operator().michelson_contrast - 0.6).abs().max() < 1e-12
    for c, k in enumerate(alg.k_factors):
        neg = (P.signs[:, c, :k] < 0).sum(-1).numpy()
        counts = np.bincount(neg, minlength=k + 1)
        pk = 1 / (k + 1)
        z = (counts - B * pk) / math.sqrt(B * pk * (1 - pk))
        assert len(counts) == k + 1 and np.abs(z).max() < 4.5


def test_11_10_11_haar_moments():
    """Lemma 11.10 and Proposition 11.11: for Haar U on M_1 (+) M_2 (+) M_3, E Tr(AU) = 0 and E|Tr(AU)|^2 =
    sum_c ||A_c||_2^2 / k_c; E[U_ij conj(U_lm)] = delta_il delta_jm / k on a block (z-scores)."""
    alg = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')
    torch.manual_seed(17)
    N = 40000
    U = alg.random_unitary_operator(N)
    A = alg.from_blocks([torch.randn(1, k, k, dtype=c128) for k in (1, 2, 3)])
    t = (A @ U).trace
    for part in (t.real, t.imag):
        assert abs(part.mean().item()) < 5 * part.std().item() / math.sqrt(N)
    exact = sum((A.matrix[0, c, :k, :k].abs() ** 2).sum().item() / k for c, k in enumerate((1, 2, 3)))
    t2 = t.abs() ** 2
    assert abs(t2.mean().item() - exact) < 5 * t2.std().item() / math.sqrt(N)
    U3 = U.matrix[:, 2, :3, :3].reshape(N, 9)
    second = torch.einsum('ba,bc->ac', U3, U3.conj()) / N                       # E U_a conj(U_c)
    se = torch.sqrt(torch.einsum('ba,bc->ac', U3.abs() ** 2, U3.abs() ** 2) / N / N)
    assert ((second - torch.eye(9, dtype=c128) / 3).abs() < 5 * se + 1e-12).all()
    V = alg.random_unitary_operator(1)                                            # left invariance: same moments
    t2V = (A @ V @ U).trace.abs() ** 2
    exactV = sum(((A @ V).matrix[0, c, :k, :k].abs() ** 2).sum().item() / k for c, k in enumerate((1, 2, 3)))
    assert abs(exactV - exact) < 1e-10 and abs(t2V.mean().item() - exact) < 5 * t2V.std().item() / math.sqrt(N)


def test_11_exercise_2_real_haar():
    """Exercise 2 (Chapter 11): for Haar O on O(k), E Tr(AO)^2 = ||A||_2^2 / k (real algebra, z-score)."""
    Mr = TypeIAlgebra([4], [4], complex_valued=False, precision='double')
    torch.manual_seed(19)
    N = 30000
    O = Mr.random_unitary_operator(N)
    assert (O.matrix[:, 0] @ O.matrix[:, 0].mT - torch.eye(4, dtype=f64)).abs().max() < 1e-12
    Ar = Mr.operator(torch.randn(1, 1, 4, 4, dtype=f64))
    tr2 = (Ar @ O).trace ** 2
    assert abs(tr2.mean().item() - Ar.frobenius_norm().item() ** 2 / 4) < 5 * tr2.std().item() / math.sqrt(N)
    dets = torch.linalg.det(O.matrix[:, 0])
    assert abs(((dets > 0).double().mean() - 0.5).item()) < 5 * 0.5 / math.sqrt(N)  # both components of O(4)


@pytest.mark.parametrize('k', [3, 6])
def test_11_12_sampling_bounds(k):
    """Corollary 11.12: the rms/max ratio ||YX||_2 / (sqrt(k) ||YX||_1) <= min(1/sqrt k, b^2/(a^2 k)), max_V |Tr(XVY)| =
    ||YX||_1, and P(|Tr(XUY)| >= Tr(XY)) <= b^4/(a^4 k^2)."""
    rng = np.random.default_rng(k)
    a, b = 0.5, 1.0
    X, Y = _psd(rng.uniform(a, b, k), rng), _psd(rng.uniform(a, b, k), rng)
    A = Y @ X
    ratio = np.linalg.norm(A) / (math.sqrt(k) * _nuc(A))
    assert ratio <= min(1 / math.sqrt(k), b ** 2 / (a ** 2 * k)) + 1e-12
    W, s, Vh = np.linalg.svd(A)
    Vbest = (W @ Vh).conj().T                                                    # A = W' |A|, U = W'^*
    assert abs(abs(np.trace(X @ Vbest @ Y)) - _nuc(A)) < 1e-10 and _nuc(A) >= np.trace(X @ Y).real - 1e-12
    Mk = TypeIAlgebra([k], [k], precision='double')
    torch.manual_seed(k)
    U = Mk.random_unitary_operator(20000).matrix[:, 0].numpy()
    vals = np.abs(np.einsum('ij,bjk,ki->b', X, U, Y))
    assert (vals >= np.trace(X @ Y).real).mean() <= b ** 4 / (a ** 4 * k ** 2)
    assert vals.max() <= _nuc(A) + 1e-10


@pytest.mark.parametrize('k', [2, 3, 5])
def test_11_13_zero_trace_unitary(k):
    """Lemma 11.13: U = V P W^* with the cyclic shift P gives Tr(AU) = 0 for A = W Sigma V^*."""
    rng = np.random.default_rng(k)
    A = rng.normal(size=(k, k)) + 1j * rng.normal(size=(k, k))
    W, s, Vh = np.linalg.svd(A)
    P = np.roll(np.eye(k), 1, axis=0)
    assert np.abs(np.diag(P)).max() == 0
    U = Vh.conj().T @ P @ W.conj().T
    assert abs(np.trace(A @ U)) < 1e-12 and np.abs(U.conj().T @ U - np.eye(k)).max() < 1e-12


@pytest.mark.parametrize('ks,delta', [([2, 2], 0.5), ([2, 3], 0.0), ([3, 2, 2], 0.9)])
def test_11_14_lower_envelope_exp1(ks, delta):
    """Proposition 11.14: z >= -(D - 2 delta/(1 + delta)) on the shell, with equality for X = diag(l, 1, ..., 1), Y = 1
    and blockwise zero-trace unitaries."""
    alg = TypeIAlgebra([k * (k + 1) for k in ks], ks, precision='double')         # m_c = k_c + 1 > 1: blunt trace
    D = sum(ks)
    bound = -(D - 2 * delta / (1 + delta))
    torch.manual_seed(D)
    B = 300
    X = opt.PositiveParam(alg, B, delta=torch.full((B,), delta, dtype=f64)).operator()
    Y = opt.PositiveParam(alg, B).operator()
    U = alg.random_unitary_operator(B)
    z = (X @ U @ Y).trace.abs() - (X @ Y).trace.real
    assert (z >= bound - 1e-12).all()
    blocks_X, blocks_U = [], []
    for c, k in enumerate(ks):
        d = np.ones(k)
        if c == 0:
            d[0] = ell(delta)
        blocks_X.append(torch.tensor(np.diag(d), dtype=c128)[None])
        blocks_U.append(torch.tensor(np.roll(np.eye(k), 1, axis=0), dtype=c128)[None])   # Tr(X_c U_c) = 0
    Xe, Ue = alg.from_blocks(blocks_X), alg.from_blocks(blocks_U)
    Ye = alg.identity()
    z_eq = ((Xe @ Ue @ Ye).trace.abs() - (Xe @ Ye).trace.real).item()
    assert abs(z_eq - bound) < 1e-12 and abs(Xe.michelson_contrast.item() - delta) < 1e-12


def test_11_5_extremize_returns_attained_values():
    """Section 11.5: a multistart search returns values attained at admissible points, so the infimum of Exp. 1 found
    is never below the proved minimum of Proposition 11.14."""
    alg = TypeIAlgebra([2, 2], [2, 2], precision='double')
    torch.manual_seed(23)
    S = 6
    X = opt.PositiveParam(alg, S, delta=torch.full((S,), 0.5, dtype=f64))
    Y, U = opt.PositiveParam(alg, S), opt.UnitaryParam(alg, S)

    def z():
        x, y = X.operator(), Y.operator()
        return (x @ U.operator() @ y).trace.abs() - (x @ y).trace.real
    best = opt.extremize(z, [X, Y, U], maximize=False, rounds=3, steps=25)
    exact = -(4 - 2 * 0.5 / 1.5)
    assert (best >= exact - 1e-12).all() and best.min().item() < exact + 0.1
    with torch.no_grad():                                                         # the final point is admissible
        assert (X.operator().michelson_contrast - 0.5).abs().max() < 1e-12
        assert (z() >= best - 1e-12).all()


@pytest.mark.parametrize('k', [1, 3, 5])
def test_11_exercise_1_cayley(k):
    """Exercise 1 (Chapter 11): C(S) = (1 - S)(1 + S)^{-1} is unitary without eigenvalue -1, inverse map S = (1-C)(1+C)^{-1}."""
    rng = np.random.default_rng(k)
    A = rng.normal(size=(k, k)) + 1j * rng.normal(size=(k, k))
    S = 3 * (A - A.conj().T)
    one = np.eye(k)
    C = (one - S) @ np.linalg.inv(one + S)
    assert np.abs(C.conj().T @ C - one).max() < 1e-12 and np.abs(np.linalg.eigvals(C) + 1).min() > 1e-8
    assert np.abs((one - C) @ np.linalg.inv(one + C) - S).max() < 1e-9
    U = _haar(k, rng)
    S2 = (one - U) @ np.linalg.inv(one + U)
    assert np.abs(S2 + S2.conj().T).max() < 1e-10
    assert np.abs((one - S2) @ np.linalg.inv(one + S2) - U).max() < 1e-10


def test_11_exercise_3_direct_sums():
    """Exercise 3 (Chapter 11): z of Exp. 3 is additive over direct sums and Delta(A (+) A') = max(Delta(A), Delta(A'))."""
    torch.manual_seed(29)
    M2 = TypeIAlgebra([2], [2], precision='double')
    M22 = TypeIAlgebra([2, 2], [2, 2], precision='double')
    X1, X2 = (opt.SelfAdjointParam(M2, 5, delta=torch.full((5,), d, dtype=f64)) for d in (0.3, 0.7))
    Y1, Y2 = opt.PositiveParam(M2, 5), opt.PositiveParam(M2, 5)

    def z3(alg, x, xa, y):
        return (y @ xa @ y).trace.real - alg.operator((y @ x @ y).matrix, is_self_adjoint=True).trace_norm()
    zs = z3(M2, X1.operator(), X1.abs_operator(), Y1.operator()) + z3(M2, X2.operator(), X2.abs_operator(), Y2.operator())
    cat = lambda a, b: M22.from_blocks([a.matrix[:, 0], b.matrix[:, 0]])
    zsum = z3(M22, cat(X1.operator(), X2.operator()), cat(X1.abs_operator(), X2.abs_operator()),
              cat(Y1.operator(), Y2.operator()))
    assert (zs - zsum).abs().max() < 1e-12
    dsum = cat(X1.abs_operator(), X2.abs_operator()).michelson_contrast
    assert (dsum - 0.7).abs().max() < 1e-12


# ==========================================================================================================
# Chapter 12. Equivariant machine learning
# ==========================================================================================================
def perm_group(n):
    """S_n as permutations and as permutation matrices pi(sigma) e_x = e_{sigma(x)}."""
    perms = list(itertools.permutations(range(n)))
    mats = np.zeros((len(perms), n, n))
    for a, p in enumerate(perms):
        mats[a, list(p), list(range(n))] = 1
    return perms, mats


def cyclic_group(n):
    perms = [tuple((x + j) % n for x in range(n)) for j in range(n)]
    mats = np.zeros((n, n, n))
    for a, p in enumerate(perms):
        mats[a, list(p), list(range(n))] = 1
    return perms, mats


def tensor_power(mats, p):
    out = np.ones((len(mats), 1, 1))
    for _ in range(p):
        out = np.einsum('gab,gcd->gacbd', out, mats).reshape(len(mats), out.shape[1] * mats.shape[1], -1)
    return out


def dihedral3():
    """The two-dimensional real (non-permutation) representation of D_3 = S_3."""
    c, s = math.cos(2 * math.pi / 3), math.sin(2 * math.pi / 3)
    r = np.array([[c, -s], [s, c]])
    f = np.diag([1.0, -1.0])
    return np.array([np.linalg.matrix_power(r, a) @ np.linalg.matrix_power(f, b) for a in range(3) for b in range(2)])


def twirl_rank(G1, G2):
    """dim Hom_G(pi_1, pi_2) as the rank of the twirl W -> |G|^{-1} sum_g pi_2(g) W pi_1(g)^*."""
    n1, n2 = G1.shape[-1], G2.shape[-1]
    Pm = sum(np.kron(g2, g1.conj()) for g1, g2 in zip(G1, G2)) / len(G1)        # on row-major vec(W)
    assert np.abs(Pm @ Pm - Pm).max() < 1e-10 and np.abs(Pm - Pm.conj().T).max() < 1e-10   # Theorem 12.5
    return int(round(np.trace(Pm).real)), Pm


def orbit_count(perms, shape_points):
    return len({min(tuple(p[i] for i in t) for p in perms) for t in shape_points})


@pytest.mark.parametrize('name', ['S3_vertex_edge', 'S4_vertex', 'C5', 'S3_edge_vertex', 'C4_vertex_edge'])
def test_12_1_hom_dimension_is_number_of_orbits(name):
    """Lemma 12.1: dim Hom_G(pi_1, pi_2) for permutation representations is the number of orbits on Omega_2 x Omega_1,
    and the orbit indicator matrices are equivariant."""
    n = int(name[1])
    perms, P = (cyclic_group if name[0] == 'C' else perm_group)(n)
    p1, p2 = {'S3_vertex_edge': (1, 2), 'S4_vertex': (1, 1), 'C5': (1, 1), 'S3_edge_vertex': (2, 1),
              'C4_vertex_edge': (1, 2)}[name]
    G1, G2 = tensor_power(P, p1), tensor_power(P, p2)
    dim, _ = twirl_rank(G1, G2)
    pts = list(itertools.product(range(n), repeat=p1 + p2))
    assert dim == orbit_count(perms, pts)
    chi1, chi2 = np.trace(G1, axis1=1, axis2=2), np.trace(G2, axis1=1, axis2=2)
    assert abs(np.mean(chi1.conj() * chi2) - dim) < 1e-12                        # Theorem 12.5
    # orbit indicators: T[y, x] = 1 on one orbit of (y, x) (tuples flattened in row-major order)
    rep = {}
    for t in pts:
        rep.setdefault(min(tuple(p[i] for i in t) for p in perms), []).append(t)
    for orb in rep.values():
        T = np.zeros((n ** p2, n ** p1))
        for t in orb:
            y = np.ravel_multi_index(t[:p2], (n,) * p2) if p2 else 0
            x = np.ravel_multi_index(t[p2:], (n,) * p1)
            T[y, x] = 1
        assert max(np.abs(T @ g1 - g2 @ T).max() for g1, g2 in zip(G1, G2)) == 0
    if p1 == p2:
        Gt = torch.tensor(G1, dtype=c128)
        assert dc.commutant_dimension(Gt) == dc.commutant_basis(Gt, group=True).shape[0] == dim


@pytest.mark.parametrize('n', [2, 3, 4, 5])
def test_12_2_deepsets(n):
    """Corollary 12.2(a): pi(S_n)' = span{1, J}; C^n = C1 (+) 1^perp, trivial and standard, inequivalent."""
    _, P = perm_group(n)
    G = torch.tensor(P, dtype=c128)
    basis = dc.commutant_basis(G, group=True)
    assert basis.shape[0] == 2
    J = torch.ones(n, n, dtype=c128)
    span = torch.stack([torch.eye(n, dtype=c128).reshape(-1), J.reshape(-1)], 1)
    for X in basis:                                                               # X in span{1, J}
        coef = torch.linalg.lstsq(span, X.reshape(-1, 1), driver='gelsd').solution
        assert (span @ coef - X.reshape(-1, 1)).abs().max() < 1e-10
    D = dc.from_group(G)
    assert D.k_factors == (1, n - 1) and D.multiplicities == (1, 1)
    assert (D.central_projections()[0] - J / n).abs().max() < 1e-10


def _stirling2(r, b):
    return sum((-1) ** i * math.comb(b, i) * (b - i) ** r for i in range(b + 1)) // math.factorial(b)


@pytest.mark.parametrize('n', [2, 3, 4])
def test_12_2_bell_numbers(n):
    """Corollary 12.2(b): maps from order-p to order-q tensors on n vertices are counted by the partitions of p + q
    elements into at most n blocks (the Bell number B_{p+q} for n >= p + q)."""
    perms, P = perm_group(n)
    bell = {2: 2, 3: 5, 4: 15}
    for r in (2, 3, 4):
        fix = np.trace(P, axis1=1, axis2=2)
        by_char = np.mean(fix ** r)                                                # (1/|G|) sum_g chi_p chi_q
        parts = sum(_stirling2(r, b) for b in range(1, min(n, r) + 1))
        assert abs(by_char - parts) < 1e-12
        assert orbit_count(perms, itertools.product(range(n), repeat=r)) == parts
        if n >= r:
            assert parts == bell[r]
    assert [_stirling2(4, b) for b in (1, 2, 3, 4)] == [1, 7, 6, 1] and [_stirling2(3, b) for b in (1, 2, 3)] == [1, 3, 1]


def test_12_3_4_hom_as_corners_and_block_matrices():
    """Lemma 12.3 and Theorem 12.4: Hom_G(pi_1, pi_2) = Q_2^* pi(G)' Q_1 has dimension sum_c a_c b_c, with a_c, b_c the
    multiplicities (from characters), for S_4 on vertices -> edges (5 = B_3) and S_3 on edges -> vertices."""
    for n, (p1, p2), expect in ((4, (1, 2), 5), (3, (2, 1), 5), (3, (2, 2), 14)):   # 14 = B_4 - 1 (n = 3 < 4)
        _, P = perm_group(n)
        G1, G2 = tensor_power(P, p1), tensor_power(P, p2)
        n1 = G1.shape[-1]
        Gs = torch.tensor(np.array([scipy.linalg.block_diag(a, b) for a, b in zip(G1, G2)]), dtype=c128)
        D = dc.from_group(Gs)
        rho = D.to_algebra(Gs).matrix                                             # blocks rho_c(g)
        chi_rho = [rho[:, c, :k, :k].diagonal(dim1=-2, dim2=-1).sum(-1).numpy() for c, k in enumerate(D.k_factors)]
        chi1, chi2 = np.trace(G1, axis1=1, axis2=2), np.trace(G2, axis1=1, axis2=2)
        a = [np.mean(ch.conj() * chi1) for ch in chi_rho]
        b = [np.mean(ch.conj() * chi2) for ch in chi_rho]
        assert all(abs(x - round(x.real)) < 1e-9 for x in a + b)
        a, b = [round(x.real) for x in a], [round(x.real) for x in b]
        assert [x + y for x, y in zip(a, b)] == list(D.multiplicities)
        P1 = torch.block_diag(torch.eye(n1, dtype=c128), torch.zeros(G2.shape[-1], G2.shape[-1], dtype=c128))
        p = D.to_commutant(P1).matrix[0]
        ranks = [round(torch.linalg.eigvalsh(p[c, :m, :m]).sum().item()) for c, m in enumerate(D.multiplicities)]
        assert ranks == a                                                         # a_c = rank p_c
        dim, _ = twirl_rank(G1, G2)
        corners = D.basis('commutant')[:, n1:, :n1].reshape(D.commutant_dim, -1)
        assert dim == sum(x * y for x, y in zip(a, b)) == expect == torch.linalg.matrix_rank(corners).item()
        for B in D.basis('commutant')[:, n1:, :n1].numpy():                      # corners are equivariant
            assert max(np.abs(B @ g1 - g2 @ B).max() for g1, g2 in zip(G1, G2)) < 1e-10


def test_12_5_twirl_is_orthogonal_projection():
    """Theorem 12.5: the twirl equals the conditional expectation onto pi(G)', is idempotent, self-adjoint, and its
    rectangular version projects onto Hom_G(pi_1, pi_2)."""
    _, P = perm_group(3)
    G = torch.tensor(tensor_power(P, 2), dtype=c128)
    D = dc.from_group(G)
    gen = torch.Generator().manual_seed(1)
    W, S = (torch.randn(9, 9, dtype=c128, generator=gen) for _ in range(2))
    PW = dc.twirl(G, W)
    hs = lambda A, B: (A.conj() * B).sum()
    assert (PW - D.conditional_expectation(W, 'commutant')[0]).abs().max() < 1e-12
    assert (dc.twirl(G, PW) - PW).abs().max() < 1e-12
    assert abs(hs(S, PW) - hs(dc.twirl(G, S), W)) < 1e-10
    assert max((g @ PW - PW @ g).abs().max().item() for g in G) < 1e-12
    for B in D.basis('commutant'):
        assert abs(hs(B, W - PW)) < 1e-10
    _, P4 = perm_group(4)
    G1, G2 = P4, tensor_power(P4, 2)
    dim, Pm = twirl_rank(G1, G2)
    Wr = np.random.default_rng(0).normal(size=(16, 4))
    T = (Pm @ Wr.reshape(-1)).reshape(16, 4)
    assert max(np.abs(T @ g1 - g2 @ T).max() for g1, g2 in zip(G1, G2)) < 1e-12 and dim == 5


def test_12_6_star_isomorphism_and_functional_calculus():
    """Proposition 12.6: Phi(w) = V((+) 1 (x) w_c)V^* is a *-isomorphism onto pi(G)', preserves unitarity and
    positivity, and f(Phi(w)) = Phi(f(w))."""
    _, P = perm_group(3)
    G = torch.tensor(tensor_power(P, 2), dtype=c128)
    D = dc.from_group(G)
    CA = D.commutant_algebra()
    gen = torch.Generator().manual_seed(2)
    a = CA.from_blocks([torch.randn(1, m, m, dtype=c128, generator=gen) for m in D.multiplicities])
    b = CA.from_blocks([torch.randn(1, m, m, dtype=c128, generator=gen) for m in D.multiplicities])
    Phi = lambda x: D.from_commutant(x)[0]
    assert (Phi(a @ b) - Phi(a) @ Phi(b)).abs().max() < 1e-12
    assert (Phi(a.adjoint()) - Phi(a).mH).abs().max() < 1e-12
    assert (Phi(CA.identity()) - torch.eye(9, dtype=c128)).abs().max() < 1e-12
    assert max((g @ Phi(a) - Phi(a) @ g).abs().max().item() for g in G) < 1e-12
    assert (D.to_commutant(Phi(a)).matrix - a.matrix).abs().max() < 1e-12      # injective
    U = Phi(CA.random_unitary_operator(1))
    assert (U.mH @ U - torch.eye(9, dtype=c128)).abs().max() < 1e-12
    pos = Phi(a.adjoint() @ a)
    assert torch.linalg.eigvalsh(pos).min() > -1e-12
    h = CA.operator((a.matrix + a.matrix.mH), is_self_adjoint=True)
    E = Phi(h.apply_function(torch.exp))
    assert (E - torch.linalg.matrix_exp(Phi(h))).abs().max() < 1e-10
    S = Phi(CA.operator(a.matrix - a.matrix.mH))                                   # (c): e^S for skew S
    assert (torch.linalg.matrix_exp(S).mH @ torch.linalg.matrix_exp(S) - torch.eye(9, dtype=c128)).abs().max() < 1e-12
    # a non-positive w gives a non-positive Phi(w): <xi, Phi(w) xi> = <v, w_c v>
    w_neg = CA.from_blocks([-torch.eye(m, dtype=c128)[None] if c == 0 else torch.eye(m, dtype=c128)[None]
                            for c, m in enumerate(D.multiplicities)])
    assert torch.linalg.eigvalsh(Phi(w_neg)).min() < -0.5


@pytest.mark.parametrize('case', ['S3_edges', 'S4_vertices', 'C5', 'D3'])
def test_12_7_invariants(case):
    """Proposition 12.7: P_0 = |G|^{-1} sum_g pi(g) projects onto the invariant vectors (dimension mean chi), the
    invariant Hermitian forms have real dimension sum m_c^2, and the invariant quadratic polynomials number
    |G|^{-1} sum (chi(g)^2 + chi(g^2))/2 (orbits on multisets for permutation representations)."""
    if case == 'D3':
        Gm, perms = dihedral3(), None
    else:
        n = int(case[1])
        perms, P = (cyclic_group if case[0] == 'C' else perm_group)(n)
        p = 2 if case == 'S3_edges' else 1
        Gm = tensor_power(P, p)
    d = Gm.shape[-1]
    chi = np.trace(Gm, axis1=1, axis2=2)
    P0 = Gm.mean(0)
    assert np.abs(P0 @ P0 - P0).max() < 1e-12 and np.abs(P0 - P0.T).max() < 1e-12
    assert abs(np.trace(P0) - chi.mean()) < 1e-12
    G = torch.tensor(Gm, dtype=c128)
    D = dc.from_group(G)
    rho = D.to_algebra(G).matrix
    triv = [c for c, k in enumerate(D.k_factors) if k == 1 and (rho[:, c, 0, 0] - 1).abs().max() < 1e-10]
    if chi.mean() > 0.5:
        assert len(triv) == 1 and (D.central_projections()[triv[0]].numpy() - P0).__abs__().max() < 1e-10
        assert D.multiplicities[triv[0]] == round(chi.mean())
    else:
        assert triv == []
    # invariant Hermitian forms: real solutions of pi(g)^* A pi(g) = A with A = A^*
    basis = []
    for i in range(d):
        for j in range(i, d):
            for ph in ((1.0,) if i == j else (1.0, 1j)):
                A = np.zeros((d, d), dtype=complex)
                A[i, j] = ph
                A[j, i] = np.conj(ph)
                basis.append(A)
    Mlin = np.array([np.concatenate([np.concatenate([(g.conj().T @ A @ g - A).real.ravel(),
                                                     (g.conj().T @ A @ g - A).imag.ravel()]) for g in Gm])
                     for A in basis]).T                                         # real-linear in the d^2 coordinates
    assert d * d - np.linalg.matrix_rank(Mlin, tol=1e-9) == D.commutant_dim
    chi_sq = np.trace(Gm @ Gm, axis1=1, axis2=2)
    q2 = np.mean((chi ** 2 + chi_sq) / 2)
    assert abs(q2 - _invariant_poly_dim(Gm, 2)) < 1e-9
    if perms is not None:
        sym = list(itertools.combinations_with_replacement(range(d), 2))
        act = [g.argmax(0) for g in Gm]
        assert abs(q2 - len({min(tuple(sorted(a[x] for x in m)) for a in act) for m in sym})) < 1e-12
    if case == 'D3':
        assert abs(q2 - 1) < 1e-12                                                # x^2 + y^2


def _invariant_poly_dim(Gm, deg, seed=0):
    """dim C[x]^G_deg by brute force: the trace of the Reynolds operator p -> mean_g p(pi(g) x) on the monomials."""
    d = Gm.shape[-1]
    monos = list(itertools.combinations_with_replacement(range(d), deg))
    rng = np.random.default_rng(seed)
    pts = rng.normal(size=(len(monos) + 10, d))
    ev = lambda X: np.stack([np.prod(X[:, list(m)], axis=1) for m in monos], 1) if deg else np.ones((len(X), 1))
    Ev = ev(pts)
    R = np.linalg.lstsq(Ev, np.mean([ev(pts @ g.T) for g in Gm], 0), rcond=None)[0]
    return np.trace(R).real


def _molien(Gm, dmax):
    lam = np.linalg.eigvals(Gm)
    c = np.zeros((len(Gm), dmax + 1), dtype=complex)
    c[:, 0] = 1
    for i in range(lam.shape[1]):
        for dd in range(1, dmax + 1):
            c[:, dd] = c[:, dd] + lam[:, i] * c[:, dd - 1]
    return c.mean(0)


@pytest.mark.parametrize('case', ['S3', 'C4', 'S3_edges', 'D3'])
def test_12_8_molien(case):
    """Proposition 12.8: the Molien series sum_d dim C[x]^G_d t^d = |G|^{-1} sum_g det(1 - t pi(g))^{-1}; for permutation
    representations the coefficients count orbits of monomials; for D_3 on R^2 the series is 1/((1-t^2)(1-t^3))."""
    if case == 'D3':
        Gm, dmax = dihedral3(), 6
    else:
        n = int(case[1])
        perms, P = (cyclic_group if case[0] == 'C' else perm_group)(n)
        Gm, dmax = tensor_power(P, 2 if case.endswith('edges') else 1), (2 if case.endswith('edges') else 4)
    series = _molien(Gm, dmax)
    t = 0.3
    closed = np.mean([1 / np.linalg.det(np.eye(Gm.shape[-1]) - t * g) for g in Gm])
    assert abs(sum(series[dd] * t ** dd for dd in range(dmax + 1)) - closed) < 3 * t ** (dmax + 1) * Gm.shape[-1] ** dmax
    for dd in range(dmax + 1):
        assert abs(series[dd] - _invariant_poly_dim(Gm, dd)) < 1e-8
        if case != 'D3':
            act = [g.argmax(0) for g in Gm]
            monos = itertools.combinations_with_replacement(range(Gm.shape[-1]), dd)
            assert abs(series[dd] - len({min(tuple(sorted(a[x] for x in m)) for a in act) for m in monos})) < 1e-9
    if case == 'D3':
        assert np.abs(series.real - np.array([1, 0, 1, 1, 1, 1, 2])).max() < 1e-12
    if case == 'S3':
        assert np.abs(series.real - np.array([1, 1, 2, 3, 4])).max() < 1e-12    # partitions into <= 3 parts


@pytest.mark.parametrize('n,d1,d2', [(3, 2, 3), (4, 1, 2)])
def test_12_exercise_1_deepsets_with_channels(n, d1, d2):
    """Exercise 1 (Chapter 12): Hom(P (x) 1_{d1}, P (x) 1_{d2}) = {1 (x) A + J (x) B}, of dimension 2 d1 d2."""
    _, P = perm_group(n)
    G1 = np.array([np.kron(p, np.eye(d1)) for p in P])
    G2 = np.array([np.kron(p, np.eye(d2)) for p in P])
    dim, Pm = twirl_rank(G1, G2)
    assert dim == 2 * d1 * d2
    rng = np.random.default_rng(n)
    A, B = rng.normal(size=(d2, d1)), rng.normal(size=(d2, d1))
    T = np.kron(np.eye(n), A) + np.kron(np.ones((n, n)), B)
    assert max(np.abs(T @ g1 - g2 @ T).max() for g1, g2 in zip(G1, G2)) < 1e-12
    W = rng.normal(size=(n * d2, n * d1))
    T2 = (Pm @ W.reshape(-1)).reshape(n * d2, n * d1).real.reshape(n, d2, n, d1)
    B2 = T2[0, :, 1, :]
    A2 = T2[0, :, 0, :] - B2
    assert all(np.abs(T2[i, :, j, :] - (A2 if i == j else 0) - B2).max() < 1e-12 for i in range(n) for j in range(n))


def test_12_exercise_2_real_layers():
    """Exercise 2 (Chapter 12): for real pi_i the real equivariant matrices have real dimension dim_C Hom_G, and
    T -> Re T maps Hom_G onto them."""
    _, P = perm_group(3)
    G1, G2 = P, tensor_power(P, 2)
    dim, Pm = twirl_rank(G1, G2)
    Mlin = np.concatenate([np.kron(g2, np.eye(3)) - np.kron(np.eye(9), g1.T) for g1, g2 in zip(G1, G2)])  # real system
    assert 27 - np.linalg.matrix_rank(Mlin) == dim
    rng = np.random.default_rng(1)
    T = (Pm @ (rng.normal(size=27) + 1j * rng.normal(size=27))).reshape(9, 3)
    ReT = T.real
    assert max(np.abs(ReT @ g1 - g2 @ ReT).max() for g1, g2 in zip(G1, G2)) < 1e-12
    assert np.abs(Mlin @ ReT.reshape(-1)).max() < 1e-12


def test_12_exercise_3_nonlinearities():
    """Exercise 3 (Chapter 12): coordinatewise f commutes with permutations, ReLU not with the rotation by pi/2, and
    x -> phi(||x||) x with every unitary representation (here D_3 on R^2 and S_3 on edges)."""
    relu = lambda v: np.maximum(v, 0)
    _, P = perm_group(3)
    G = tensor_power(P, 2)
    x = np.random.default_rng(3).normal(size=9)
    assert all(np.array_equal(relu(g @ x), g @ relu(x)) for g in G)
    R = np.array([[0.0, -1.0], [1.0, 0.0]])
    v = np.array([1.0, -1.0])
    assert np.array_equal(relu(R @ v), [1.0, 1.0]) and np.array_equal(R @ relu(v), [0.0, 1.0])
    gate = lambda v: np.tanh(np.linalg.norm(v)) * v
    for Gm in (dihedral3(), G):
        y = np.random.default_rng(4).normal(size=Gm.shape[-1])
        assert max(np.abs(gate(g @ y) - g @ gate(y)).max() for g in Gm) < 1e-12


@pytest.mark.parametrize('n', [3, 4])
def test_12_exercise_4_samples_determine_layers(n):
    """Exercise 4 (Chapter 12): the equivariant T with T x_1 = ... = T x_N = 0 form a space of dimension
    sum_c m_c max(0, m_c - N k_c) for Gaussian samples."""
    _, P = perm_group(n)
    G = torch.tensor(tensor_power(P, 2), dtype=c128)
    D = dc.from_group(G)
    basis = D.basis('commutant')
    gen = torch.Generator().manual_seed(n)
    for N in (1, 2, 3):
        Xs = torch.randn(n * n, N, dtype=c128, generator=gen)
        Mmap = torch.stack([(B @ Xs).reshape(-1) for B in basis], dim=1)
        lost = sum(m * max(0, m - N * k) for k, m in zip(D.k_factors, D.multiplicities))
        assert torch.linalg.matrix_rank(Mmap).item() == D.commutant_dim - lost
    if n == 4:
        assert (D.k_factors, D.multiplicities) == ((1, 2, 3, 3), (2, 1, 1, 3))
        # one sample never suffices: rank S_c <= k_c = 1 < m_c = 2 in the sector (1, 2), for any x_1
        x1 = torch.zeros(16, 1, dtype=c128)
        x1[0] = 1
        Mmap = torch.stack([(B @ x1).reshape(-1) for B in basis], dim=1)
        assert torch.linalg.matrix_rank(Mmap).item() <= D.commutant_dim - 2
