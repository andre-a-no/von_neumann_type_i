"""
Tests for torch_vn_algebra.su2 (SU(2) multiplets of spin-1/2 chains).

Independent references: dense / sparse Kronecker-product Hamiltonians built here with NumPy/SciPy (2^L basis,
site i = bit i), the U(1) route of chains.SpinChain, exact Heisenberg-ring energies tabulated in scripts/chains.py,
the Lieb-Mattis theorem, textbook Clebsch-Gordan values and SymPy.
"""
import ast
import math
import pathlib
from fractions import Fraction

import numpy as np
import pytest
import scipy.sparse as sp
import torch

from torch_vn_algebra import SpinChain
from torch_vn_algebra.su2 import (SU2Chain, clebsch_gordan, clebsch_gordan_exact, multiplicity, two_spins,
                                  wigner_eckart)

ROOT = pathlib.Path(__file__).resolve().parents[1]

# ----------------------------------------------------------------------------------------------------------
# independent references
# ----------------------------------------------------------------------------------------------------------
SP = sp.csr_matrix(np.array([[0.0, 0.0], [1.0, 0.0]]))         # S^+ = |up><down|, basis (|0> = down, |1> = up)
SM = sp.csr_matrix(SP.T)
SZ = sp.csr_matrix(np.diag([-0.5, 0.5]))


def site_op(L, i, A):
    """A on site i of L sites; the 2^L index is the label sum_i bit_i 2^i (site i = bit i)."""
    return sp.kron(sp.kron(sp.identity(2 ** (L - 1 - i)), A), sp.identity(2 ** i)).tocsr()


def heis_bond(L, i, j):
    return (0.5 * (site_op(L, i, SP) @ site_op(L, j, SM) + site_op(L, i, SM) @ site_op(L, j, SP))
            + site_op(L, i, SZ) @ site_op(L, j, SZ)).tocsr()


def nn_bonds(L, boundary):
    b = [(i, i + 1) for i in range(L - 1)]
    if boundary == 'periodic' and L > 2:
        b.append((L - 1, 0))
    return b


def nnn_bonds(L, boundary):
    if boundary == 'open' or L < 3:
        return [(i, i + 2) for i in range(L - 2)]
    pairs = []
    for i in range(L):
        p = frozenset((i, (i + 2) % L))
        if p not in pairs:
            pairs.append(p)
    return [tuple(sorted(p)) for p in pairs]


def sparse_hamiltonian(L, boundary, J=1.0, J2=0.0, h=0.0):
    H = sp.csr_matrix((2 ** L, 2 ** L))
    for i, j in nn_bonds(L, boundary):
        H = H + J * heis_bond(L, i, j)
    if J2:
        for i, j in nnn_bonds(L, boundary):
            H = H + J2 * heis_bond(L, i, j)
    if h:
        for i in range(L):
            H = H + h * site_op(L, i, SZ)
    return H.tocsr()


def casimir_dense(L):
    """S^2 = sum_{i != j} S_i.S_j + 3L/4."""
    S2 = 0.75 * L * np.eye(2 ** L)
    for i in range(L):
        for j in range(L):
            if i != j:
                S2 = S2 + heis_bond(L, i, j).toarray()
    return S2


def popcount(x):
    return np.array([bin(int(v)).count('1') for v in x])


def reference_spectrum(L, boundary, J=1.0, J2=0.0, h=0.0):
    """All 2^L levels: full dense diagonalisation for L <= 10, U(1) sectors of the sparse matrix above."""
    H = sparse_hamiltonian(L, boundary, J, J2, h)
    if L <= 10:
        return np.sort(np.linalg.eigvalsh(H.toarray()))
    pc = popcount(np.arange(2 ** L))
    out = []
    for N in range(L + 1):
        idx = np.nonzero(pc == N)[0]
        out.append(np.linalg.eigvalsh(H[idx][:, idx].toarray()))
    return np.sort(np.concatenate(out))


def u1_spectrum(L, boundary, J=1.0, h=0.0):
    c = SpinChain(L, boundary)
    E = [torch.linalg.eigvalsh(c.xxz(J, 1.0, h, sector=N).matrix[0, 0]) for N in range(L + 1)]
    return torch.sort(torch.cat(E))[0].numpy()


def to_full(chain, N, v):
    """Vector(s) of the U(1) sector N of chain.u1 -> 2^L vectors (index = label)."""
    v = np.asarray(v)
    out = np.zeros(v.shape[:-1] + (2 ** chain.L,), dtype=v.dtype)
    out[..., chain.u1.labels[N].numpy()] = v
    return out


def ring_values():
    """HEISENBERG_RING (exact ground-state energies of the Heisenberg ring, 9 digits) from scripts/chains.py."""
    tree = ast.parse((ROOT / 'scripts' / 'chains.py').read_text())
    for node in tree.body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], 'id', None) == 'HEISENBERG_RING':
            return ast.literal_eval(node.value)
    raise AssertionError("HEISENBERG_RING not found")


# ----------------------------------------------------------------------------------------------------------
# multiplicities and highest-weight bases
# ----------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('L', range(1, 15))
def test_multiplicities(L):
    total = 0
    for t in two_spins(L):
        S = Fraction(t, 2)
        n = L // 2 - S if L % 2 == 0 else Fraction(L, 2) - S
        m = math.comb(L, int(n)) - (math.comb(L, int(n) - 1) if n >= 1 else 0)
        assert multiplicity(L, S) == m == multiplicity(L, float(S))
        total += (t + 1) * m
    assert total == 2 ** L
    assert multiplicity(L, Fraction(L + 2, 2)) == 0 and multiplicity(L, Fraction(L - 1, 2)) == 0
    # m_S = d_{L/2+S} - d_{L/2+S+1}: the S^z sectors
    for t in two_spins(L):
        N = (L + t) // 2
        assert multiplicity(L, t / 2) == math.comb(L, N) - (math.comb(L, N + 1) if N < L else 0)


@pytest.mark.parametrize('L', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
def test_highest_weight_basis_against_dense_null_space(L):
    c = SU2Chain(L)
    Sp = sum(site_op(L, i, SP) for i in range(L)).toarray()
    Sz = sum(site_op(L, i, SZ) for i in range(L)).toarray()
    S2 = casimir_dense(L)
    pc = popcount(np.arange(2 ** L))
    for t, m in zip(c.two_spins, c.multiplicities):
        N = (L + t) // 2
        B = c.hw_basis(t / 2).numpy()
        assert B.shape == (math.comb(L, N), m)
        assert np.abs(B.T @ B - np.eye(m)).max() < 1e-13                      # orthonormal
        F = to_full(c, N, B.T).T                                               # (2^L, m)
        assert np.abs(Sp @ F).max() < 1e-13                                    # highest weight
        assert np.abs(Sz @ F - (t / 2) * F).max() < 1e-13
        assert np.abs(S2 @ F - (t / 2) * (t / 2 + 1) * F).max() < 1e-12
        # the span is the whole null space of S^+ in the sector (SVD reference)
        idx = np.nonzero(pc == N)[0]
        s = np.linalg.svd(Sp[:, idx], compute_uv=False)
        assert np.sum(s > 1e-10) == len(idx) - m
        assert len(c.paths(t / 2)) == m


def test_yamanouchi_paths():
    c = SU2Chain(8)
    for t in c.two_spins:
        P = c.paths(t / 2)
        assert len(set(P)) == len(P)
        for p in P:
            assert p[0] == 1 and p[-1] == t and len(p) == 8
            assert all(abs(a - b) == 1 and b >= 0 for a, b in zip(p[:-1], p[1:]))


@pytest.mark.parametrize('L', [5, 6])
def test_multiplet_basis_and_embedding(L):
    c = SU2Chain(L, 'periodic')
    Sm = sum(site_op(L, i, SM) for i in range(L)).toarray()
    S2 = casimir_dense(L)
    for t in c.two_spins:
        S = Fraction(t, 2)
        F = to_full(c, (L + t) // 2, c.hw_basis(S).numpy().T).T
        for tm in range(t, -t - 1, -2):
            X = c.multiplet_basis(S, Fraction(tm, 2)).numpy()
            G = to_full(c, (L + tm) // 2, X.T).T
            assert np.abs(G - F).max() < 1e-13                               # |S, M> by the dense S^-
            assert np.abs(X.T @ X - np.eye(X.shape[1])).max() < 1e-13
            assert np.abs(S2 @ G - float(S * (S + 1)) * G).max() < 1e-12      # S^2 sanity
            if tm > -t:
                F = Sm @ F / math.sqrt((t + tm) * (t - tm + 2) / 4)
    # embed / multiplet with coefficients
    t = c.two_spins[1]
    coef = torch.randn(c.multiplicity(t / 2), dtype=torch.float64)
    states = c.multiplet(t / 2, coef)
    assert len(states) == t + 1
    for N, v in states:
        assert abs(torch.linalg.vector_norm(v).item() - torch.linalg.vector_norm(coef).item()) < 1e-12
    batch = c.embed(t / 2, torch.stack([coef, 2 * coef]))
    assert torch.allclose(batch[1], 2 * c.embed(t / 2, coef))


# ----------------------------------------------------------------------------------------------------------
# spectra
# ----------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('boundary', ['open', 'periodic'])
@pytest.mark.parametrize('L', range(4, 13))
def test_full_spectrum(L, boundary):
    c = SU2Chain(L, boundary)
    cases = [dict(J=1.0, h=0.0, J2=0.0), dict(J=1.0, h=0.37, J2=0.0), dict(J=0.8, h=-0.21, J2=0.41)]
    for p in cases:
        H = c.heisenberg(**p)
        E = c.full_spectrum(H)[0].numpy()
        assert E.shape == (2 ** L,)
        ref = reference_spectrum(L, boundary, **p)
        assert np.abs(E - ref).max() < 1e-10, p
        if p['J2'] == 0:
            assert np.abs(E - u1_spectrum(L, boundary, p['J'], p['h'])).max() < 1e-10


def test_quantum_numbers_and_degeneracies():
    c = SU2Chain(8, 'periodic')
    E, S, M = c.full_spectrum(c.heisenberg(), return_quantum_numbers=True)
    for t, m in zip(c.two_spins, c.multiplicities):
        sel = S[0] == t / 2
        assert int(sel.sum()) == (t + 1) * m
        for tm in range(-t, t + 1, 2):
            assert int((sel & (M[0] == tm / 2)).sum()) == m
    # with a field: E(S, M) = E(S, S) + h (M - S)
    E2, S2, M2 = c.full_spectrum(c.heisenberg(h=0.3), return_quantum_numbers=True)
    shifted = torch.sort(E2[0] - 0.3 * M2[0])[0]
    assert torch.allclose(shifted, E[0], atol=1e-12)


def test_batched_couplings_and_fields():
    c = SU2Chain(7, 'periodic')
    J = torch.tensor([1.0, -0.5, 2.0], dtype=torch.float64)
    J2 = torch.tensor([0.0, 0.3, 0.5], dtype=torch.float64)
    h = torch.tensor([0.1, 0.0, -0.4], dtype=torch.float64)
    H = c.heisenberg(J=J, J2=J2, h=h)
    E = c.full_spectrum(H)
    for b in range(3):
        Hb = c.heisenberg(J=J[b].item(), J2=J2[b].item(), h=h[b].item())
        assert torch.allclose(H.matrix[b], Hb.matrix[0], atol=1e-13)
        ref = reference_spectrum(7, 'periodic', J[b].item(), J2[b].item(), h[b].item())
        assert np.abs(E[b].numpy() - ref).max() < 1e-10


def test_explicit_bonds_random_couplings():
    L = 8
    c = SU2Chain(L)
    g = torch.Generator().manual_seed(1)
    bonds = [(0, 3), (1, 2), (2, 7), (4, 5), (3, 6), (0, 7), (5, 6), (1, 4)]
    Jr = torch.rand(4, len(bonds), generator=g, dtype=torch.float64) * 2 - 0.5
    H = c.heisenberg(J=Jr, bonds=bonds)
    E = c.full_spectrum(H)
    for b in range(4):
        ref = sum(Jr[b, k].item() * heis_bond(L, i, j) for k, (i, j) in enumerate(bonds))
        assert np.abs(E[b].numpy() - np.linalg.eigvalsh(ref.toarray())).max() < 1e-10
        H1 = c.heisenberg(J=Jr[b], bonds=bonds)                                # (n_bonds,): one Hamiltonian
        assert torch.allclose(H1.matrix[0], H.matrix[b], atol=1e-13)
    # S^2 = 2 sum_{i<j} S_i.S_j + 3L/4: the Casimir
    pairs = [(i, j) for i in range(L) for j in range(i + 1, L)]
    S2 = c.heisenberg(J=2.0, bonds=pairs)
    assert torch.allclose(S2.matrix + 0.75 * L * c.algebra.identity().matrix, c.casimir().matrix, atol=1e-12)
    with pytest.raises(ValueError):
        c.heisenberg(J=torch.ones(3), bonds=bonds)
    with pytest.raises(ValueError):
        c.heisenberg(bonds=[(0, 0)])
    with pytest.raises(ValueError):
        c.heisenberg(J2=0.5, bonds=bonds)


def test_restrict_and_lift_against_u1():
    L = 8
    for boundary in ('open', 'periodic'):
        c = SU2Chain(L, boundary)
        u = c.u1
        H = c.heisenberg(J=1.3, h=0.2)
        R = c.restrict(u.xxz(1.3, 1.0, 0.2))
        assert torch.allclose(R.matrix, H.matrix, atol=1e-13)
        H0 = c.heisenberg(J=1.3)
        lifted = c.lift(H0)
        assert torch.allclose(lifted.matrix, u.xxz(1.3, 1.0).matrix, atol=1e-13)
        assert torch.allclose(c.restrict(lifted).matrix, H0.matrix, atol=1e-13)


def test_single_spin_block():
    c = SU2Chain(10, 'periodic')
    H = c.heisenberg(J2=0.3, h=0.1)
    for k, t in enumerate(c.two_spins):
        Hs = c.heisenberg(J2=0.3, h=0.1, spin=t / 2)
        m = c.multiplicities[k]
        assert Hs.algebra.k_factors == [m] and Hs.algebra.charges == [t]
        assert torch.allclose(Hs.matrix[:, 0], H.matrix[:, k, :m, :m], atol=1e-14)
    with pytest.raises(ValueError):
        c.full_spectrum(c.heisenberg(spin=0))


def test_two_sites_singlet_triplet():
    c = SU2Chain(2)
    E, S, M = c.full_spectrum(c.heisenberg(J=1.7), return_quantum_numbers=True)
    assert torch.allclose(E[0], torch.tensor([-0.75, 0.25, 0.25, 0.25], dtype=torch.float64) * 1.7, atol=1e-15)
    assert S[0].tolist() == [0.0, 1.0, 1.0, 1.0]
    assert c.multiplicities == [1, 1]


@pytest.mark.parametrize('L', [4, 6, 8, 10, 12, 14])
def test_heisenberg_ring_exact_energies(L):
    exact = ring_values()[L]
    c = SU2Chain(L, 'periodic')
    gs = c.ground_state(c.heisenberg())
    assert abs(gs.energy.item() - exact) < 1e-9
    assert gs.spin.item() == 0.0
    # Lanczos on the U(1) sector agrees (independent route)
    if L <= 10:
        from torch_vn_algebra import krylov
        E, _, _ = krylov.ground_state(c.u1.xxz_sparse(sector=L // 2), generator=torch.Generator().manual_seed(0))
        assert abs(E.item() - gs.energy.item()) < 1e-9


@pytest.mark.parametrize('L', range(2, 13))
def test_lieb_mattis(L):
    """Bipartite chains (open, and periodic with even L): the ground state has S = 0 (L even) / 1/2 (L odd) and
    the lowest energy E_0(S) increases strictly with S (Lieb and Mattis 1962)."""
    for boundary in ('open', 'periodic'):
        if boundary == 'periodic' and L % 2:
            continue
        c = SU2Chain(L, boundary)
        H = c.heisenberg()
        gs = c.ground_state(H)
        assert gs.spin.item() == (L % 2) / 2
        e = c.lowest_energies(H)[0]
        assert torch.all(torch.diff(e) > 1e-8), (boundary, e)
    # the ferromagnet: ground multiplet S = L/2 with energy J * (number of bonds) / 4
    c = SU2Chain(L, 'open')
    gs = c.ground_state(c.heisenberg(J=-1.0))
    assert gs.spin.item() == L / 2 and abs(gs.energy.item() + (L - 1) / 4) < 1e-12


def test_ground_state_embedding_and_field():
    L = 10
    c = SU2Chain(L, 'periodic')
    for h in (0.0, 0.3, -0.3, 2.5):                         # 2.5 > saturation field 2J of the ring
        H = c.heisenberg(h=h)
        gs = c.ground_state(H)
        S, sz = gs.spin.item(), gs.sz.item()
        m = c.multiplicity(S)
        psi = c.embed(S, gs.coefficients[0, :m], M=sz)
        N = int(round(L / 2 + sz))
        Hn = c.u1.xxz(1.0, 1.0, h, sector=N).matrix[0, 0]
        res = torch.linalg.vector_norm(Hn @ psi - gs.energy[0] * psi).item()
        assert res < 1e-10 and abs(torch.linalg.vector_norm(psi).item() - 1) < 1e-12
        assert abs(gs.energy.item() - reference_spectrum(L, 'periodic', h=h)[0]) < 1e-10
        if h == 2.5:
            assert S == L / 2 and sz == -L / 2
            assert abs(gs.energy.item() - (L / 4 - h * L / 2)) < 1e-12


def test_majumdar_ghosh_point():
    """J2 = J/2, periodic even L: two exactly degenerate dimer singlets with E = -3 L / 8 J."""
    for L in (8, 10, 12):
        c = SU2Chain(L, 'periodic')
        H = c.heisenberg(J=1.0, J2=0.5)
        e = H.eigenvalues()[0][0]                             # block S = 0
        assert abs(e[0].item() + 3 * L / 8) < 1e-10 and abs(e[1].item() + 3 * L / 8) < 1e-10
        assert c.ground_state(H).spin.item() == 0.0


# ----------------------------------------------------------------------------------------------------------
# Clebsch-Gordan coefficients
# ----------------------------------------------------------------------------------------------------------
def halves(n):
    return [Fraction(k, 2) for k in range(n)]


def ms(j):
    return [Fraction(k, 2) for k in range(-int(2 * j), int(2 * j) + 1, 2)]


def test_clebsch_gordan_textbook_values():
    h = Fraction(1, 2)
    s2 = 1 / math.sqrt(2)
    # 1/2 x 1/2
    assert clebsch_gordan(h, h, h, h, 1, 1) == 1.0
    assert abs(clebsch_gordan(h, h, h, -h, 1, 0) - s2) < 2e-16
    assert abs(clebsch_gordan(h, -h, h, h, 1, 0) - s2) < 2e-16
    assert abs(clebsch_gordan(h, h, h, -h, 0, 0) - s2) < 2e-16
    assert abs(clebsch_gordan(h, -h, h, h, 0, 0) + s2) < 2e-16
    # 1 x 1/2
    assert abs(clebsch_gordan(1, 1, h, -h, Fraction(3, 2), h) - math.sqrt(1 / 3)) < 1e-15
    assert abs(clebsch_gordan(1, 0, h, h, Fraction(3, 2), h) - math.sqrt(2 / 3)) < 1e-15
    assert abs(clebsch_gordan(1, 1, h, -h, h, h) - math.sqrt(2 / 3)) < 1e-15
    assert abs(clebsch_gordan(1, 0, h, h, h, h) + math.sqrt(1 / 3)) < 1e-15
    # 1 x 1
    assert abs(clebsch_gordan(1, 1, 1, -1, 0, 0) - math.sqrt(1 / 3)) < 1e-15
    assert abs(clebsch_gordan(1, 0, 1, 0, 0, 0) + math.sqrt(1 / 3)) < 1e-15
    assert abs(clebsch_gordan(1, 1, 1, -1, 1, 0) - s2) < 1e-15
    assert abs(clebsch_gordan(1, 0, 1, 0, 1, 0)) == 0.0
    assert abs(clebsch_gordan(1, 1, 1, -1, 2, 0) - math.sqrt(1 / 6)) < 1e-15
    assert abs(clebsch_gordan(1, 0, 1, 0, 2, 0) - math.sqrt(2 / 3)) < 1e-15
    assert clebsch_gordan_exact(1, 1, 1, -1, 2, 0) == (1, Fraction(1, 6))
    assert clebsch_gordan_exact(1, 0, 1, 0, 0, 0) == (-1, Fraction(1, 3))
    # selection rules and float / int / Fraction arguments
    assert clebsch_gordan(1, 1, 1, 1, 1, 1) == 0.0                         # M != m1 + m2
    assert clebsch_gordan(1, 0, 1, 0, 3, 0) == 0.0                         # triangle
    assert clebsch_gordan(0.5, 0.5, 0.5, -0.5, 1.0, 0.0) == clebsch_gordan(h, h, h, -h, 1, 0)
    with pytest.raises(ValueError):
        clebsch_gordan(0.3, 0.3, 0.5, 0.5, 1, 1)
    with pytest.raises(ValueError):
        clebsch_gordan(-1, 0, 1, 0, 0, 0)


@pytest.mark.parametrize('j1,j2', [(a, b) for a in halves(6) for b in halves(6)])
def test_clebsch_gordan_orthogonality(j1, j2):
    Js = [abs(j1 - j2) + k for k in range(int(j1 + j2 - abs(j1 - j2)) + 1)]
    rows = [(m1, m2) for m1 in ms(j1) for m2 in ms(j2)]
    cols = [(J, M) for J in Js for M in ms(J)]
    U = np.array([[clebsch_gordan(j1, m1, j2, m2, J, M) for (J, M) in cols] for (m1, m2) in rows])
    assert U.shape[0] == U.shape[1]
    assert np.abs(U.T @ U - np.eye(len(cols))).max() < 1e-14              # sum_{m1 m2}: delta_JJ' delta_MM'
    assert np.abs(U @ U.T - np.eye(len(rows))).max() < 1e-14              # sum_{J M}: delta_m1m1' delta_m2m2'
    # exchange symmetry <j1 m1 j2 m2|J M> = (-1)^{j1+j2-J} <j2 m2 j1 m1|J M>
    for (m1, m2) in rows:
        for (J, M) in cols:
            sign = -1 if int(j1 + j2 - J) % 2 else 1
            assert clebsch_gordan(j1, m1, j2, m2, J, M) == sign * clebsch_gordan(j2, m2, j1, m1, J, M)


def test_clebsch_gordan_against_sympy():
    cg_sympy = pytest.importorskip('sympy.physics.quantum.cg')
    from sympy import Rational
    R = lambda x: Rational(x.numerator, x.denominator)                     # noqa: E731
    for j1 in halves(6):
        for j2 in halves(5):
            for J in [abs(j1 - j2) + k for k in range(int(j1 + j2 - abs(j1 - j2)) + 1)]:
                for m1 in ms(j1):
                    for m2 in ms(j2):
                        if abs(m1 + m2) > J:
                            continue
                        ref = float(cg_sympy.CG(R(j1), R(m1), R(j2), R(m2), R(J), R(m1 + m2)).doit())
                        assert abs(clebsch_gordan(j1, m1, j2, m2, J, m1 + m2) - ref) < 1e-14


# ----------------------------------------------------------------------------------------------------------
# Wigner-Eckart theorem and correlations
# ----------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('L,boundary', [(6, 'open'), (7, 'periodic')])
def test_wigner_eckart_matrix_elements(L, boundary):
    c = SU2Chain(L, boundary)
    R = c.reduced_matrix_elements()
    Sm = sum(site_op(L, i, SM) for i in range(L)).toarray()
    # dense multiplets |S, M, a> by repeated S^- (columns), keyed by (2S, 2M)
    mult = {}
    for t in c.two_spins:
        F = to_full(c, (L + t) // 2, c.hw_basis(t / 2).numpy().T).T
        for tm in range(t, -t - 1, -2):
            mult[(t, tm)] = F
            F = Sm @ F / math.sqrt((t + tm) * (t - tm + 2) / 4) if tm > -t else None
    T = {1: lambda i: -site_op(L, i, SP).toarray() / math.sqrt(2), 0: lambda i: site_op(L, i, SZ).toarray(),
         -1: lambda i: site_op(L, i, SM).toarray() / math.sqrt(2)}
    g = np.random.default_rng(0)
    for i in (0, L // 2):
        for q in (-1, 0, 1):
            Tq = T[q](i)
            for (tp, t), Ri in R.items():
                for tm in range(-t, t + 1, 2):
                    tmp = tm + 2 * q
                    if abs(tmp) > tp:
                        continue
                    direct = mult[(tp, tmp)].T @ Tq @ mult[(t, tm)]
                    we = wigner_eckart(Ri[i].numpy(), tp / 2, tmp / 2, 1, q, t / 2, tm / 2)
                    assert np.abs(direct - we).max() < 1e-12
                    # random multiplets: <psi'|T_q|psi> for random coefficient vectors
                    a, b = g.standard_normal(Ri.shape[1]), g.standard_normal(Ri.shape[2])
                    assert abs(a @ direct @ b - a @ we @ b) < 1e-11
            # selection rule |S' - S| <= 1: other blocks vanish
            for t in c.two_spins:
                for tp in c.two_spins:
                    if abs(tp - t) > 2:
                        assert np.abs(mult[(tp, min(tp, t))].T @ Tq @ mult[(t, min(tp, t))]).max() < 1e-13


def test_reduced_elements_sum_rules():
    """Total spin: sum_i <S'||S_i||S> = delta_SS' sqrt(S(S+1)(2S+1)) 1 (Edmonds 5.4.3); hermiticity
    R_{S,S'} = (-1)^{S'-S} R_{S',S}^T; Edmonds sum rule sum_{M M' q} |<S'M'|T_q|SM>|^2 = |<S'||T||S>|^2."""
    c = SU2Chain(9, 'open')
    R = c.reduced_matrix_elements()
    for (tp, t), Ri in R.items():
        total = Ri.sum(0)
        if tp == t:
            S = t / 2
            assert torch.allclose(total, math.sqrt(S * (S + 1) * (2 * S + 1)) * torch.eye(total.shape[0],
                                  dtype=total.dtype), atol=1e-12)
        else:
            assert total.abs().max() < 1e-12
            assert torch.allclose(R[(t, tp)], -R[(tp, t)].mT, atol=1e-13)    # (-1)^{S'-S} = -1
    # the normalisation behind the sum rule: sum_{M, q} <S M; 1 q | S' M+q>^2 = 2S' + 1
    for (tp, t) in R:
        s = sum(clebsch_gordan(t / 2, tm / 2, 1, q, tp / 2, tm / 2 + q) ** 2 / (tp + 1)
                for tm in range(-t, t + 1, 2) for q in (-1, 0, 1))
        assert abs(s - 1) < 1e-14 or (t == tp == 0)


def test_spin_correlations_wigner_eckart_vs_direct_vs_u1():
    L = 10
    c = SU2Chain(L, 'periodic')
    pairs = [(0, r) for r in range(L)] + [(3, 7), (8, 2)]
    C_we = c.spin_correlation(pairs)
    C_dir = c.spin_correlation(pairs, method='direct')
    assert torch.allclose(C_we.matrix, C_dir.matrix, atol=1e-12)
    # closed form c(S, S'') = (-1)^{S-S''} / (2S+1)
    for t in range(0, 9):
        for tpp in (t - 2, t, t + 2):
            if tpp >= 0 and not (t == tpp == 0):
                assert abs(SU2Chain._scalar_coefficient(t, tpp) - (-1) ** ((t - tpp) // 2) / (t + 1)) < 1e-14
    # ground state of the ring: <S_0.S_1> = E_0 / L, sum_r <S_0.S_r> = <S_0.S_tot> = 0 in a singlet
    H = c.heisenberg()
    gs = c.ground_state(H)
    S, m = gs.spin.item(), c.multiplicity(gs.spin.item())
    corr = c.expectation(C_we, S, gs.coefficients[0, :m].expand(len(pairs), -1))
    assert abs(corr[1].item() - gs.energy.item() / L) < 1e-12
    assert abs(corr[:L].sum().item()) < 1e-12
    assert abs(corr[0].item() - 0.75) < 1e-14
    # the U(1) route: dense <psi|S_i.S_j|psi> on the embedded vector
    psi = to_full(c, L // 2, c.embed(S, gs.coefficients[0, :m]).numpy())
    for k, (i, j) in enumerate(pairs):
        op = heis_bond(L, i, j) if i != j else 0.75 * sp.identity(2 ** L)
        assert abs(psi @ (op @ psi) - corr[k].item()) < 1e-12
    # excited multiplets: S = 1, random coefficient vectors, against the U(1) route
    g = torch.Generator().manual_seed(3)
    a = torch.randn(c.multiplicity(1), generator=g, dtype=torch.float64)
    a = a / torch.linalg.vector_norm(a)
    val = c.expectation(C_we, 1, a.expand(len(pairs), -1))
    for M in (1, 0, -1):
        v = to_full(c, L // 2 + M, c.embed(1, a, M=M).numpy())
        for k, (i, j) in enumerate(pairs[:4]):
            op = heis_bond(L, i, j) if i != j else 0.75 * sp.identity(2 ** L)
            assert abs(v @ (op @ v) - val[k].item()) < 1e-12
    # sum over all pairs: S^2
    allp = [(i, j) for i in range(L) for j in range(L)]
    tot = c.spin_correlation(allp).matrix.sum(0)
    assert torch.allclose(tot, c.casimir().matrix[0], atol=1e-11)


# ----------------------------------------------------------------------------------------------------------
# number fields, precision, properties
# ----------------------------------------------------------------------------------------------------------
def test_real_and_complex_algebras_and_hermiticity():
    for cv in (False, True):
        c = SU2Chain(8, 'periodic', complex_valued=cv)
        H = c.heisenberg(J2=0.3, h=0.2)
        assert H.matrix.dtype == (torch.complex128 if cv else torch.float64)
        assert H.is_self_adjoint
        assert torch.equal(H.matrix, H.matrix.mH)
        if cv:
            assert H.matrix.imag.abs().max() == 0
        assert c.algebra.charges == c.two_spins and c.algebra.k_factors == c.multiplicities
        spectra = spectra + [c.full_spectrum(H)] if cv else [c.full_spectrum(H)]
    assert torch.allclose(spectra[0], spectra[1], atol=1e-12)


def test_single_precision():
    L = 10
    c32 = SU2Chain(L, 'periodic', precision='single')
    c64 = SU2Chain(L, 'periodic', precision='double')
    H32, H64 = c32.heisenberg(J2=0.2, h=0.1), c64.heisenberg(J2=0.2, h=0.1)
    assert H32.matrix.dtype == torch.complex64
    E32, E64 = c32.full_spectrum(H32), c64.full_spectrum(H64)
    assert E32.dtype == torch.float32
    assert (E32.double() - E64).abs().max() < 2e-5
    R = c32.reduced_matrix_elements(sites=[0, 3])
    assert R[(2, 0)].dtype == torch.float32 and R[(2, 0)].shape[0] == 2
    C = c32.spin_correlation([(0, 1)])
    assert (C.matrix.to(torch.complex128) - c64.spin_correlation([(0, 1)]).matrix).abs().max() < 1e-5


def test_invalid_arguments():
    c = SU2Chain(6)
    with pytest.raises(ValueError):
        c.hw_basis(0.5)                                     # half-integer spin for even L
    with pytest.raises(ValueError):
        c.multiplet_basis(1, 2)
    with pytest.raises(ValueError):
        c.sector_index(4)
    with pytest.raises(ValueError):
        SU2Chain(0)
    with pytest.raises(ValueError):
        c.spin_correlation([(0, 1)], method='other')
    assert c.u1_sector(1, -1) == 2 and c.u1_sector(Fraction(2)) == 5


@pytest.mark.skipif(not torch.cuda.is_available(), reason="needs CUDA")
def test_cuda():
    c = SU2Chain(10, 'periodic', device='cuda')
    E = c.full_spectrum(c.heisenberg(J2=0.2)).cpu().numpy()
    assert np.abs(E[0] - reference_spectrum(10, 'periodic', J2=0.2)).max() < 1e-10


def test_no_host_synchronisations():
    """Hamiltonians, spectra, ground states and correlations read no values back to the host (bases cached)."""
    from test_sync_free import count_host_reads
    c = SU2Chain(8, 'periodic')
    c.hw_basis(0)
    c.algebra
    with count_host_reads() as calls:
        H = c.heisenberg(J=torch.tensor([1.0, 0.5], dtype=torch.float64), J2=0.3, h=0.1)
        c.full_spectrum(H)
        c.ground_state(H)
        c.lowest_energies(H)
        c.spin_correlation([(0, 1), (0, 2)])
        c.spin_correlation([(0, 3)], method='direct')
    assert calls == []


def test_pruned_single_spin_bases_match_full_recursion():
    full = SU2Chain(9).build_bases()
    for t in full.two_spins:
        single = SU2Chain(9)
        assert torch.equal(single.hw_basis(t / 2), full.hw_basis(t / 2))
        assert single.paths(t / 2) == full.paths(t / 2)
        assert list(single._bases) == [t]                   # only the requested spin is kept
