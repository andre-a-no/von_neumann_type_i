"""
Review tests of the structure modules (decompose, inclusions, jordan, divergences, su2) against independent
constructions: group representations with known k_c and m_c (characters computed by hand), dense Kronecker
products of spin matrices, scipy.linalg matrix functions, sympy Clebsch-Gordan coefficients, brute-force
Pimsner-Popa constants and Choi matrices. All in double precision on the CPU, with fixed seeds.
"""
import itertools
import math
import warnings
from fractions import Fraction

import numpy as np
import pytest
import scipy.linalg as sla
import torch

import torch_vn_algebra.decompose as dc
import torch_vn_algebra.divergences as dv
import torch_vn_algebra.inclusions as inc
import torch_vn_algebra.jordan as jd
import torch_vn_algebra.su2 as su
from torch_vn_algebra import TypeIAlgebra
from torch_vn_algebra.channels import random_inter_sector_channel, random_channel
from torch_vn_algebra.states import fidelity, random_density_matrix

C128 = torch.complex128
F64 = torch.float64


# ======================================================================================================
# helpers (independent of the modules under test)
# ======================================================================================================
def gen(seed):
    return torch.Generator().manual_seed(seed)


def haar(n, g, dtype=C128):
    Q, R = torch.linalg.qr(torch.randn(n, n, dtype=dtype, generator=g))
    return Q * torch.diagonal(R).sgn()


def herm(n, g, real=False):
    A = torch.randn(n, n, dtype=F64 if real else C128, generator=g)
    return (A + A.mH) / 2


def perm_matrix(p):
    """Permutation matrix e_i -> e_{p[i]}."""
    n = len(p)
    P = torch.zeros(n, n, dtype=F64)
    for i, j in enumerate(p):
        P[j, i] = 1
    return P


def cyclic_regular(n):
    """Regular representation of Z_n."""
    return torch.stack([perm_matrix([(i + s) % n for i in range(n)]) for s in range(n)])


def dihedral_square():
    """D_4 acting on the 4 vertices of a square: 8 permutation matrices."""
    r = [1, 2, 3, 0]
    f = [0, 3, 2, 1]                                   # reflection through vertices 0 and 2
    els, cur = [], list(range(4))
    for _ in range(4):
        els.append(list(cur))
        els.append([cur[f[i]] for i in range(4)])
        cur = [r[cur[i]] for i in range(4)]
    return torch.stack([perm_matrix(p) for p in els])


def quaternion_group():
    """Q_8 acting on R^4 = H by left multiplication (real, quaternionic type)."""
    def L(a, b, c, d):
        return torch.tensor([[a, -b, -c, -d], [b, a, -d, c], [c, d, a, -b], [d, -c, b, a]], dtype=F64)
    units = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
    return torch.stack([L(*[s * x for x in u]) for u in units for s in (1, -1)])


def rotation_group(n_fold, copies=1):
    """Z_n rotating R^2 (complex type) tensored with 1_copies, plus a trivial line."""
    out = []
    for j in range(n_fold):
        c, s = math.cos(2 * math.pi * j / n_fold), math.sin(2 * math.pi * j / n_fold)
        R = torch.kron(torch.tensor([[c, -s], [s, c]], dtype=F64), torch.eye(copies, dtype=F64))
        out.append(torch.block_diag(R, torch.ones(1, 1, dtype=F64)))
    return torch.stack(out)


def character_dim(G):
    """|G|^-1 sum_g |tr g|^2 in numpy."""
    chi = np.trace(G.numpy(), axis1=1, axis2=2)
    return float(np.mean(np.abs(chi) ** 2))


def dense_spin_ops(L, s=0.5):
    """Site spins (L, 3, d^L, d^L) by numpy Kronecker products, site 0 the most significant factor."""
    d = int(round(2 * s)) + 1
    m = np.array([s - a for a in range(d)])
    sp = np.zeros((d, d))
    for a in range(1, d):
        sp[a - 1, a] = math.sqrt(s * (s + 1) - m[a] * (m[a] + 1))
    ops = [(sp + sp.T) / 2, (sp - sp.T) / 2j, np.diag(m)]
    out = np.zeros((L, 3, d ** L, d ** L), dtype=complex)
    for i in range(L):
        for a in range(3):
            out[i, a] = np.kron(np.kron(np.eye(d ** i), ops[a]), np.eye(d ** (L - 1 - i)))
    return out


def label_spin_ops(L):
    """Site spins (L, 3, 2^L, 2^L) in the SpinChain label basis: basis index = label, bit i of the label = site i
    (1 = up), built from Kronecker products with site 0 the least significant factor."""
    sz = np.diag([-0.5, 0.5])                    # index 0 = down, 1 = up
    sp = np.array([[0., 0.], [1., 0.]])          # |up><down|
    sx, sy = (sp + sp.T) / 2, (sp - sp.T) / 2j
    out = np.zeros((L, 3, 2 ** L, 2 ** L), dtype=complex)
    for i in range(L):
        for a, o in enumerate((sx, sy, sz)):
            out[i, a] = np.kron(np.kron(np.eye(2 ** (L - 1 - i)), o), np.eye(2 ** i))
    return out


def op_from(alg, blocks):
    o = alg.from_blocks([torch.as_tensor(b) for b in blocks])
    return o


def pos_op(alg, mats):
    o = alg.from_blocks([torch.as_tensor(m) for m in mats])
    o._is_self_adjoint = o._is_positive = o._is_normal = True
    return o


def blocks_of(op, b=0):
    return [op.matrix[b, c, :k, :k] for c, k in enumerate(op.algebra.k_factors) if k]


def wishart(k, g, cond=None):
    G = torch.randn(k, k, dtype=C128, generator=g)
    W = G @ G.mH
    if cond is not None:
        w, V = torch.linalg.eigh(W)
        w = torch.logspace(0, -math.log10(cond), k, dtype=F64)
        W = (V * w.to(C128)) @ V.mH
    return (W + W.mH) / 2


# ======================================================================================================
# decompose: Wedderburn decomposition
# ======================================================================================================
GROUPS = {
    # name: (elements, expected (k_c, m_c) sorted)        independent derivation in the comments
    'S3_qubits': (lambda: dc.permutation_representation(3, 2), [(1, 4), (2, 2)]),          # Sym^3 C^2 (4), std x2
    'S3_qutrits': (lambda: dc.permutation_representation(3, 3), [(1, 1), (1, 10), (2, 8)]),  # Schur-Weyl GL_3 dims
    'S4_qubits': (lambda: dc.permutation_representation(4, 2), [(1, 5), (2, 1), (3, 3)]),   # (4), (2,2), (3,1)
    'Z5_regular': (lambda: cyclic_regular(5), [(1, 1)] * 5),                                # regular rep, abelian
    'D4_square': (dihedral_square, [(1, 1), (1, 1), (2, 1)]),                               # chi = (4,0,0,2,0)
    'Q8_real': (quaternion_group, [(2, 2)]),                                                # H (x) C = 2 x C^2
}


@pytest.mark.parametrize("name", list(GROUPS))
def test_group_sectors_and_characters(name):
    """Proposition 2.7: k_c = dims of the irreps, m_c = multiplicities; dim pi(G)' = |G|^-1 sum|chi|^2 = sum m_c^2."""
    make, expected = GROUPS[name]
    G = make()
    D = dc.from_group(G)
    assert sorted(zip(D.k_factors, D.multiplicities)) == sorted(expected)
    assert dc.commutant_dimension(G) == round(character_dim(G)) == sum(m * m for _, m in expected)
    assert sum(k * m for k, m in expected) == G.shape[-1]
    assert D.verify()['ok']


@pytest.mark.parametrize("name", ['S3_qubits', 'S4_qubits', 'D4_square'])
def test_group_blocks_are_inequivalent_irreps(name):
    """Proposition 2.7: rho_c is a homomorphism, irreducible (spans M_{k_c}) and inequivalent (Schur orthogonality)."""
    G = GROUPS[name][0]().to(C128)
    D = dc.from_group(G)
    blocks = D.to_algebra(G).matrix                                    # (|G|, C, k, k)
    torch.testing.assert_close(D.from_algebra(D.to_algebra(G)), G, atol=1e-12, rtol=0)
    # homomorphism on all pairs
    prod = torch.einsum('gab,hbc->ghac', G, G).reshape(-1, G.shape[-1], G.shape[-1])
    lhs = D.to_algebra(prod).matrix.reshape(G.shape[0], G.shape[0], D.C, *blocks.shape[-2:])
    rhs = torch.einsum('gcab,hcbd->ghcad', blocks, blocks)
    torch.testing.assert_close(lhs, rhs, atol=1e-12, rtol=0)
    chis = []
    for c, k in enumerate(D.k_factors):
        X = blocks[:, c, :k, :k].reshape(G.shape[0], -1)
        assert torch.linalg.matrix_rank(X) == k * k                    # Burnside: rho_c(G) spans M_k
        chis.append(torch.diagonal(blocks[:, c, :k, :k], dim1=-2, dim2=-1).sum(-1))
    chis = torch.stack(chis)
    gram = chis.conj() @ chis.T / G.shape[0]
    torch.testing.assert_close(gram, torch.eye(D.C, dtype=C128), atol=1e-12, rtol=0)


@pytest.mark.parametrize("L,s,expected", [
    (3, 0.5, [(2, 2), (4, 1)]),
    (4, 0.5, [(1, 2), (3, 3), (5, 1)]),
    (5, 0.5, [(2, 5), (4, 4), (6, 1)]),
    (3, 1.0, [(1, 1), (3, 3), (5, 2), (7, 1)]),       # 1 x 1 x 1 = 0 + 3*1 + 2*2 + 3 (S^z counts 1, 3, 6, 7)
])
def test_total_spin_schur_weyl(L, s, expected):
    """Theorem 2.6 for SU(2) generators: k = 2S + 1 and m_S = #(S^z = S) - #(S^z = S + 1)."""
    D = dc.from_generators(dc.total_spin(L, s))
    assert sorted(zip(D.k_factors, D.multiplicities)) == sorted(expected)
    if s == 0.5:
        for k, m in zip(D.k_factors, D.multiplicities):
            assert m == su.multiplicity(L, Fraction(k - 1, 2))
    # independent: the generators built by numpy Kronecker products give the same algebra
    ops = dense_spin_ops(L, s)
    S_tot = torch.as_tensor(ops.sum(0))
    torch.testing.assert_close(S_tot, dc.total_spin(L, s), atol=1e-12, rtol=0)


def test_spin_matrices_commutation():
    """Spin matrices: [S_x, S_y] = i S_z and S^2 = s(s + 1) 1 for s = 1/2, ..., 2."""
    for s in (0.5, 1, 1.5, 2):
        Sx, Sy, Sz = dc.spin_matrices(s)
        torch.testing.assert_close(Sx @ Sy - Sy @ Sx, 1j * Sz, atol=1e-13, rtol=0)
        d = Sx.shape[0]
        torch.testing.assert_close(Sx @ Sx + Sy @ Sy + Sz @ Sz, s * (s + 1) * torch.eye(d, dtype=C128),
                                   atol=1e-13, rtol=0)


def test_permutation_representation_homomorphism():
    """pi(sigma tau) = pi(sigma) pi(tau) with (sigma tau)(i) = sigma(tau(i)), and the vector in slot i moves to sigma(i)."""
    L, d = 3, 3
    perms = list(itertools.permutations(range(L)))
    P = dc.permutation_representation(L, d)
    idx = {p: i for i, p in enumerate(perms)}
    for a, b in itertools.product(perms, perms):
        ab = tuple(a[b[i]] for i in range(L))
        torch.testing.assert_close(P[idx[a]] @ P[idx[b]], P[idx[ab]], atol=0, rtol=0)
    g = gen(0)
    vs = [torch.randn(d, dtype=C128, generator=g) for _ in range(L)]
    for p in perms:
        x = torch.kron(torch.kron(vs[0], vs[1]), vs[2])
        out = [None] * L
        for i in range(L):
            out[p[i]] = vs[i]
        y = torch.kron(torch.kron(out[0], out[1]), out[2])
        torch.testing.assert_close(P[idx[p]] @ x, y, atol=1e-14, rtol=0)


def test_commutant_side_twirl_and_conditional_expectations():
    """Section 2.6: E onto A' is the twirl; E_A, E_A' are Hilbert-Schmidt projections onto A, A' (bases orthonormal)."""
    G = dc.permutation_representation(3, 2)
    D = dc.from_group(G)
    g = gen(1)
    X = torch.randn(4, 8, 8, dtype=C128, generator=g)
    torch.testing.assert_close(D.conditional_expectation(X, 'commutant'), dc.twirl(G, X), atol=1e-13, rtol=0)
    for which, dim in (('algebra', D.dim), ('commutant', D.commutant_dim)):
        Bs = D.basis(which)
        assert Bs.shape[0] == dim
        gram = torch.einsum('aij,bij->ab', Bs.conj(), Bs)
        torch.testing.assert_close(gram, torch.eye(dim, dtype=C128), atol=1e-13, rtol=0)
        E = D.conditional_expectation(X, which)
        proj = torch.einsum('aij,bij,akl->bkl', Bs.conj(), X, Bs)        # sum_a <B_a, X> B_a
        torch.testing.assert_close(E, proj, atol=1e-13, rtol=0)
    # A and A' commute; the central projections are a partition of unity commuting with both
    P = D.central_projections()
    torch.testing.assert_close(P.sum(0), torch.eye(8, dtype=C128), atol=1e-13, rtol=0)
    torch.testing.assert_close(torch.einsum('cij,djk->cdik', P, P),
                               torch.einsum('cd,cik->cdik', torch.eye(D.C, dtype=C128), P), atol=1e-13, rtol=0)
    assert [round(float(torch.trace(p).real)) for p in P] == list(D.sector_dims)
    for Y in (G.to(C128), D.basis('commutant')):
        torch.testing.assert_close(torch.einsum('cij,gjk->cgik', P, Y), torch.einsum('gij,cjk->cgik', Y, P),
                                   atol=1e-12, rtol=0)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_unitary_invariance_and_seed_independence(seed):
    """Theorem 2.6: (k_c, m_c) are invariants: unchanged under W . W^* and under the random seed."""
    g = gen(10 + seed)
    G = dc.permutation_representation(4, 2).to(C128)
    W = haar(16, g)
    D0 = dc.from_group(G)
    D1 = dc.from_group(W @ G @ W.mH, seed=seed + 5)
    T = dc.permutation_representation(4, 2, 'adjacent').to(C128)     # the transpositions (i, i+1) generate S_4
    D2 = dc.from_generators(W @ T @ W.mH, seed=seed)
    for D in (D1, D2):
        assert (D.k_factors, D.multiplicities) == (D0.k_factors, D0.multiplicities)
        assert D.verify()['ok']


@pytest.mark.parametrize("scale", [1e-6, 1.0, 1e6])
def test_generator_scaling_invariance(scale):
    """The algebra generated by S equals that of c S (c != 0): decomposition invariant under 1e-6 and 1e6."""
    S = dc.total_spin(4)
    D = dc.from_generators(S * scale)
    assert (D.k_factors, D.multiplicities) == ((1, 3, 5), (2, 3, 1))
    assert D.diagnostics['residuals']['algebra'] < 1e-12


def test_single_precision_decomposition():
    """Single precision: same sectors, residuals at the 1e-6 level, complex64 result."""
    D = dc.from_generators(dc.total_spin(4).to(torch.complex64))
    assert (D.k_factors, D.multiplicities) == ((1, 3, 5), (2, 3, 1))
    assert D.V.dtype == torch.complex64
    assert D.diagnostics['residuals']['algebra'] < 1e-5


@pytest.mark.parametrize("delta", [1e-10, 1e-12])
def test_regression_near_identity_generator_respects_tol(delta):
    """Regression (decompose tolerance): diag(1, 1 + delta) with delta far below tol is a multiple of 1 up to tol,
    so its commutant is M_2; the probe clustering used to separate eigenvalues 1e-12 apart (floor 1e3 n eps instead
    of tol) and returned C (+) C."""
    H = torch.diag(torch.tensor([1.0, 1.0 + delta], dtype=C128))
    assert dc.commutant_basis(H).shape[0] == 4
    D = dc.from_generators(H)
    assert (D.k_factors, D.multiplicities) == ((1,), (2,))
    # far above tol: the symmetry is broken, commutant = diagonal
    H2 = torch.diag(torch.tensor([1.0, 1.0 + 1e-5], dtype=C128))
    assert dc.commutant_basis(H2).shape[0] == 2
    assert dc.from_generators(H2).k_factors == (1, 1)


def test_regression_noisy_identity_with_explicit_tol():
    """Regression (decompose tolerance): 1 + 1e-10 N (N Hermitian, generic) with tol = 1e-6 generates C 1."""
    g = gen(3)
    N = herm(5, g)
    S = torch.eye(5, dtype=C128) + 1e-10 * N
    assert dc.commutant_basis(S, tol=1e-6).shape[0] == 25
    D = dc.from_generators(S, tol=1e-6)
    assert (D.k_factors, D.multiplicities) == ((1,), (5,))


def test_non_star_closed_and_edge_cases():
    """The generated algebra is the *-algebra of S and S^*: E_12 generates M_2; n = 1, the identity, zero."""
    E12 = torch.tensor([[0, 1], [0, 0]], dtype=C128)
    D = dc.from_generators(E12)
    assert (D.k_factors, D.multiplicities) == ((2,), (1,))
    assert (dc.from_generators(torch.ones(1, 1)).k_factors, dc.from_group(torch.ones(1, 1, 1)).k_factors) == ((1,), (1,))
    D = dc.from_generators(torch.eye(3))
    assert (D.k_factors, D.multiplicities) == ((1,), (3,))
    D = dc.from_generators(torch.zeros(2, 2))
    assert (D.k_factors, D.multiplicities) == ((1,), (2,))
    # the rotation generator: C over R, C (+) C over C (module docstring)
    D = dc.from_generators(torch.tensor([[0., -1.], [1., 0.]]))
    assert (D.k_factors, D.multiplicities) == ((1, 1), (1, 1))
    assert D.real_structure().types == ('C', 'C')


@pytest.mark.parametrize("ks,ms", [((1,), (1,)), ((2,), (3,)), ((1, 2, 3), (2, 1, 2)), ((2, 2), (1, 3))])
def test_hidden_block_algebra(ks, ms):
    """Theorem 2.6 on random generators of W ((+)_c M_{k_c} (x) 1_{m_c}) W^*: sectors recovered, A' = (+)_c 1 (x) M_{m_c}."""
    g = gen(sum(ks) * 7 + sum(ms))
    n = sum(k * m for k, m in zip(ks, ms))
    W = haar(n, g)
    gens = []
    for _ in range(3):
        X = torch.block_diag(*[torch.kron(torch.randn(k, k, dtype=C128, generator=g), torch.eye(m, dtype=C128))
                               for k, m in zip(ks, ms)])
        gens.append(W @ X @ W.mH)
    D = dc.from_generators(torch.stack(gens))
    assert sorted(zip(D.k_factors, D.multiplicities)) == sorted(zip(ks, ms))
    Y = torch.block_diag(*[torch.kron(torch.eye(k, dtype=C128), torch.randn(m, m, dtype=C128, generator=g))
                           for k, m in zip(ks, ms)])
    Y = W @ Y @ W.mH
    torch.testing.assert_close(D.from_commutant(D.to_commutant(Y))[0], Y, atol=1e-11, rtol=0)
    assert D.commutant_dim == dc.commutant_basis(torch.stack(gens)).shape[0]


# ----- real structure (Section 4.7) -----
@pytest.mark.parametrize("make,types,sizes", [
    (lambda: dc.permutation_representation(3, 2), ('R', 'R'), None),
    (quaternion_group, ('H',), None),
    (lambda: torch.stack([torch.kron(q, torch.eye(2, dtype=F64)) for q in quaternion_group()]), ('H',), None),
    (lambda: rotation_group(3), ('C', 'R', 'C'), None),
    (lambda: rotation_group(4, copies=2), None, None),
])
def test_real_structure_types(make, types, sizes):
    """Section 4.7: types R / C / H from J conj(J) = +-1 and the conjugation partner; real_basis is real orthogonal
    and puts real elements into the forms alpha (x) 1 (R) and [[Re x, Im x], [-Im x, Re x]] (C)."""
    G = make()
    D = dc.from_group(G)
    rs = D.real_structure()
    assert rs.ok
    if types is not None:
        assert rs.types == types
    for c, t in enumerate(rs.types):
        d = rs.partner[c]
        assert rs.partner[d] == c and ((t == 'C') == (d != c))
        if t == 'H':
            assert D.k_factors[c] % 2 == 0 and D.multiplicities[c] % 2 == 0
    O, blocks = D.real_basis()
    assert O.dtype == F64
    torch.testing.assert_close(O.T @ O, torch.eye(O.shape[0], dtype=F64), atol=1e-12, rtol=0)
    Gr = G.real.to(F64) if torch.is_complex(G) else G.to(F64)
    for b in blocks:
        cols = O[:, b['columns']]
        for g in Gr:
            X = cols.T @ g @ cols
            # invariant subspace: the block is all of O^T g O restricted to these columns
            torch.testing.assert_close(g @ cols, cols @ X, atol=1e-12, rtol=0)
            if b['type'] == 'R':
                c = b['sectors'][0]
                k, m = D.k_factors[c], D.multiplicities[c]
                alpha = X.reshape(k, m, k, m)[:, 0, :, 0]
                torch.testing.assert_close(X, torch.kron(alpha, torch.eye(m, dtype=F64)), atol=1e-12, rtol=0)
            if b['type'] == 'C':
                h = X.shape[0] // 2
                torch.testing.assert_close(X[:h, :h], X[h:, h:], atol=1e-12, rtol=0)
                torch.testing.assert_close(X[:h, h:], -X[h:, :h], atol=1e-12, rtol=0)


def test_real_structure_rejects_algebra_not_closed_under_conjugation():
    """Section 4.7 hypothesis conj(A) = A: the projection onto (1, e^{i pi/3})/sqrt2 generates an algebra that is not."""
    v = torch.tensor([1, complex(math.cos(math.pi / 3), math.sin(math.pi / 3))], dtype=C128) / math.sqrt(2)
    D = dc.from_generators(torch.outer(v, v.conj()))
    with pytest.raises(ValueError):
        D.real_structure()


# ======================================================================================================
# inclusions: conditional expectations, index, Pimsner-Popa, basic construction
# ======================================================================================================
INCLUSIONS = [
    ([1, 2], [[1, 1], [2, 1], [0, 2]]),
    ([1, 2], [[3, 0], [2, 1], [0, 3]]),
    ([2], [[1], [2]]),
    ([1], [[3]]),
    ([1, 1, 1], [[1, 1, 1]]),
]


def make_inclusion(ks, lam, seed=0, unitary=True):
    lam = np.asarray(lam)
    N = TypeIAlgebra(list(ks), list(ks), precision='double')
    m = (lam @ np.asarray(ks)).tolist()
    M = TypeIAlgebra(m, m, precision='double')
    U = None
    if unitary:
        torch.manual_seed(seed)
        U = M.random_unitary_operator()
    return inc.Inclusion(N, M, lam, unitary=U)


def rand_op(alg, g, batch=1):
    X = torch.randn(batch, alg.C, alg.k_max, alg.k_max, dtype=C128, generator=g)
    return alg.operator(X)


def tau(alg, w, op):
    """sum_c w_c Tr_c (batch,)."""
    tr = torch.diagonal(op.matrix, dim1=-2, dim2=-1).sum(-1)               # (B, C)
    return (tr * torch.as_tensor(w, dtype=C128)).sum(-1)


@pytest.mark.parametrize("ks,lam", INCLUSIONS)
def test_embedding_is_unital_star_homomorphism(ks, lam):
    """iota = Ad(u) o iota_0 is a unital injective *-homomorphism; m_d = sum_c lambda_dc k_c."""
    I = make_inclusion(ks, lam)
    g = gen(1)
    x, y = rand_op(I.N, g), rand_op(I.N, g)
    torch.testing.assert_close(I.embed(x @ y).matrix, (I.embed(x) @ I.embed(y)).matrix, atol=1e-12, rtol=0)
    torch.testing.assert_close(I.embed(x.adjoint()).matrix, I.embed(x).adjoint().matrix, atol=1e-12, rtol=0)
    torch.testing.assert_close(I.embed(I.N.identity()).matrix, I.M.identity().matrix, atol=1e-12, rtol=0)
    assert float(I.embed(x).frobenius_norm().sum()) > 0


@pytest.mark.parametrize("ks,lam", INCLUSIONS)
def test_conditional_expectation_properties(ks, lam):
    """E is unital, CP (positive), E o iota = id, N-bimodular, tau o E = tau, tau(iota(n) y) = tau|_N(n E(y))."""
    I = make_inclusion(ks, lam, seed=2)
    g = gen(2)
    w = (torch.rand(I.M.C, dtype=F64, generator=g) + 0.2)
    E = I.conditional_expectation(w)
    nu = I.induced_trace(w)
    torch.testing.assert_close(nu, w @ torch.as_tensor(np.asarray(lam), dtype=F64), atol=1e-14, rtol=0)
    y, n1, n2 = rand_op(I.M, g), rand_op(I.N, g), rand_op(I.N, g)
    torch.testing.assert_close(E(I.embed(n1)).matrix, n1.matrix, atol=1e-12, rtol=0)
    torch.testing.assert_close(E(I.M.identity()).matrix, I.N.identity().matrix, atol=1e-12, rtol=0)
    torch.testing.assert_close(E(I.embed(n1) @ y @ I.embed(n2)).matrix, (n1 @ E(y) @ n2).matrix, atol=1e-11, rtol=0)
    torch.testing.assert_close(tau(I.N, nu, E(y)), tau(I.M, w, y), atol=1e-11, rtol=0)
    torch.testing.assert_close(tau(I.M, w, I.embed(n1) @ y), tau(I.N, nu, n1 @ E(y)), atol=1e-11, rtol=0)
    P = y @ y.adjoint()
    for B in blocks_of(E(P)):
        assert float(torch.linalg.eigvalsh(B).min()) > -1e-12
    # E as a channel agrees with E
    ch = E.as_channel()
    torch.testing.assert_close(ch.apply(y).matrix, E(y).matrix, atol=1e-12, rtol=0)


@pytest.mark.parametrize("ks,lam", INCLUSIONS)
def test_watatani_index_and_quasi_basis(ks, lam):
    """Watatani index: (Ind E)_d = (Lambda Lambda^T w)_d / w_d = sum_j u_j u_j^*, y = sum_j u_j iota(E(u_j^* y))."""
    I = make_inclusion(ks, lam, seed=3)
    g = gen(3)
    w = torch.rand(I.M.C, dtype=F64, generator=g) + 0.1
    E = I.conditional_expectation(w)
    L = torch.as_tensor(np.asarray(lam), dtype=F64)
    expected = (L @ L.T @ w) / w
    Ind = E.index()
    torch.testing.assert_close(Ind.matrix, I.M.central(expected.tolist()).matrix, atol=1e-11, rtol=0)
    torch.testing.assert_close(E.index('quasi_basis').matrix, Ind.matrix, atol=1e-11, rtol=0)
    Q = E.quasi_basis()[0]
    y = rand_op(I.M, g)
    acc = torch.zeros_like(y.matrix)
    for j in range(Q.shape[0]):
        u = I.M.operator(Q[j][None])
        acc = acc + (u @ I.embed(E(u.adjoint() @ y))).matrix
    torch.testing.assert_close(acc, y.matrix, atol=1e-11, rtol=0)


def _pp_brute(I, E, g, trials=60):
    """min over pure y = xi xi^* of 1 / xi^* iota(E(y))^+ xi (the best l with iota(E(y)) >= l y)."""
    M = I.M
    best = math.inf
    for d, m in enumerate(M.k_factors):
        for _ in range(trials):
            xi = torch.randn(m, dtype=C128, generator=g)
            y = torch.zeros(1, M.C, M.k_max, M.k_max, dtype=C128)
            y[0, d, :m, :m] = torch.outer(xi, xi.conj())
            A = I.embed(E(M.operator(y))).matrix[0, d, :m, :m]
            val = float((xi.conj() @ torch.linalg.pinv(A, hermitian=True) @ xi).real)
            best = min(best, 1 / val)
    return best


def _cp_brute(I, E):
    """max l with y -> iota(E(y)) - l y completely positive: Choi matrices of the diagonal blocks d -> d."""
    M = I.M
    best = math.inf
    for d, m in enumerate(M.k_factors):
        C = torch.zeros(m * m, m * m, dtype=C128)
        for i, j in itertools.product(range(m), range(m)):
            y = torch.zeros(1, M.C, M.k_max, M.k_max, dtype=C128)
            y[0, d, i, j] = 1
            C[i * m:(i + 1) * m, j * m:(j + 1) * m] = I.embed(E(M.operator(y))).matrix[0, d, :m, :m]
        om = torch.eye(m, dtype=C128).reshape(-1)
        best = min(best, 1 / float((om.conj() @ torch.linalg.pinv(C, hermitian=True) @ om).real))
    return best


@pytest.mark.parametrize("ks,lam", INCLUSIONS[:4])
def test_pimsner_popa_constants_brute_force(ks, lam):
    """Pimsner-Popa: lambda(E)^-1 = max_d sum_c (nu_c/w_d) min(k_c, lambda_dc); complete order: 1/||Ind E|| (Choi)."""
    I = make_inclusion(ks, lam, seed=4)
    g = gen(4)
    w = torch.rand(I.M.C, dtype=F64, generator=g) + 0.3
    E = I.conditional_expectation(w)
    lp = float(E.pimsner_popa_constant()[0])
    lc = float(E.pimsner_popa_constant(complete=True)[0])
    assert abs(_pp_brute(I, E, g) - lp) < 1e-9 * lp + 1e-12          # generic xi attain the minimum
    assert abs(_cp_brute(I, E) - lc) < 1e-9 * lc
    torch.testing.assert_close(torch.tensor(lc, dtype=F64), 1 / E.index().operator_norm().max(), atol=1e-12, rtol=0)


@pytest.mark.parametrize("n", [1, 2, 3, 5])
def test_index_of_factors_and_masa(n):
    """[M_m : M_k] = (m/k)^2, [M_n : C] = n^2, [M_n : C^n] = n; lambda(E) = 1/n and 1/n^2 for C c M_n."""
    M = TypeIAlgebra([n], [n], precision='double')
    sc = inc.scalars(M)
    assert abs(sc.index() - n * n) < 1e-12
    E = sc.conditional_expectation('blunt')
    assert abs(float(E.pimsner_popa_constant()[0]) - 1 / n) < 1e-14
    assert abs(float(E.pimsner_popa_constant(complete=True)[0]) - 1 / n ** 2) < 1e-14
    assert abs(inc.diagonal(n, precision='double').index() - n) < 1e-12
    N = TypeIAlgebra([2], [2], precision='double')
    assert abs(inc.from_tensor_product(N, n).index() - n * n) < 1e-12


@pytest.mark.parametrize("ks,lam", INCLUSIONS)
def test_markov_trace(ks, lam):
    """Markov trace: Lambda Lambda^T t = beta t, s = Lambda^T t, Lambda s = beta t, sum m_d t_d = 1, t > 0; Ind E = beta 1."""
    I = make_inclusion(ks, lam)
    L = np.asarray(lam, dtype=float)
    mt = I.markov_trace()
    t, s = mt.t.numpy(), mt.s.numpy()
    assert np.all(t > 0)
    np.testing.assert_allclose(L @ L.T @ t, mt.beta * t, atol=1e-12)
    np.testing.assert_allclose(L @ s, mt.beta * t, atol=1e-12)
    np.testing.assert_allclose(L.T @ t, s, atol=1e-14)
    assert abs(np.asarray(I.M.k_factors) @ t - 1) < 1e-14
    assert abs(mt.beta - np.linalg.norm(L, 2) ** 2) < 1e-10
    Ind = I.conditional_expectation('markov').index()
    torch.testing.assert_close(Ind.matrix, (I.M.identity() * mt.beta).matrix, atol=1e-11, rtol=0)


@pytest.mark.parametrize("L", [1, 2, 3, 5])
def test_chain_tower_index(L):
    """[A_{L+1} : A_L] = 4 cos^2(pi / (2L + 4)); the embedding x -> x (x) 1_2 is unital."""
    tower = inc.chain_tower(L + 1, L, precision='double')
    assert abs(tower[0].index() - 4 * math.cos(math.pi / (2 * L + 4)) ** 2) < 1e-12


@pytest.mark.parametrize("ks,lam", INCLUSIONS[:3])
def test_jones_projection_and_tower(ks, lam):
    """Basic construction: e_N projection, e j(y) e = j(iota(E(y))) e, [e, j(iota(N))] = 0, Markov property
    tau_1(j(y) e) = tau(y) / beta; the Jones tower has matrices Lambda, Lambda^T, Lambda and constant index."""
    I = make_inclusion(ks, lam, seed=5)
    g = gen(5)
    e = I.jones_projection('markov')
    b = I.basic_construction()
    j = b.embed
    E = I.conditional_expectation('markov')
    y, n = rand_op(I.M, g), rand_op(I.N, g)
    torch.testing.assert_close((e @ e).matrix, e.matrix, atol=1e-12, rtol=0)
    torch.testing.assert_close(e.adjoint().matrix, e.matrix, atol=1e-12, rtol=0)
    torch.testing.assert_close((e @ j(y) @ e).matrix, (j(I.embed(E(y))) @ e).matrix, atol=1e-11, rtol=0)
    jn = j(I.embed(n))
    torch.testing.assert_close((e @ jn).matrix, (jn @ e).matrix, atol=1e-11, rtol=0)
    mt = I.markov_trace()
    lhs = tau(b.M, mt.s / mt.beta, j(y) @ e)
    torch.testing.assert_close(lhs, tau(I.M, mt.t, y) / mt.beta, atol=1e-12, rtol=0)
    tower = I.jones_tower(4)
    L = np.asarray(lam)
    for i, T in enumerate(tower):
        np.testing.assert_array_equal(T.matrix, L if i % 2 == 0 else L.T)
        assert abs(T.index() - I.index()) < 1e-10


@pytest.mark.parametrize("ks,lam", INCLUSIONS)
def test_relative_commutant_is_maximal(ks, lam):
    """N' cap M = (+)_{lambda_dc > 0} M_{lambda_dc}: commutes with iota(N) and has the dimension of the commutant
    computed as a null space, sum lambda_dc^2."""
    I = make_inclusion(ks, lam, seed=6)
    R = I.relative_commutant()
    assert float(I.commutator_error(R).max()) < 1e-12
    L = np.asarray(lam)
    dim = 0
    # independent null space of y -> [iota(x_a), y] over a basis x_a of N, block by block of M
    basis = []
    for c, k in enumerate(I.N.k_factors):
        for a, b2 in itertools.product(range(k), range(k)):
            X = torch.zeros(1, I.N.C, I.N.k_max, I.N.k_max, dtype=C128)
            X[0, c, a, b2] = 1
            basis.append(I.embed(I.N.operator(X)).matrix[0])
    for d, m in enumerate(I.M.k_factors):
        rows = []
        for Xa in basis:
            A = Xa[d, :m, :m]
            rows.append(torch.kron(A.contiguous(), torch.eye(m, dtype=C128))
                        - torch.kron(torch.eye(m, dtype=C128), A.T.contiguous()))
        s = torch.linalg.svdvals(torch.cat(rows))
        dim += int((s < 1e-10).sum())
    assert dim == int((L ** 2).sum()) == sum(k * k for k in R.N.k_factors)


def test_compose_tower():
    """N c M c P: Lambda_PN = Lambda_PM Lambda_MN and iota_PN = iota_PM o iota_MN (with unitaries)."""
    I1 = make_inclusion([1, 2], [[1, 1], [2, 1]], seed=7)
    m = I1.M.k_factors
    lam2 = np.array([[1, 0], [1, 2], [0, 1]])
    p = (lam2 @ np.asarray(m)).tolist()
    P = TypeIAlgebra(p, p, precision='double')
    torch.manual_seed(8)
    I2 = inc.Inclusion(I1.M, P, lam2, unitary=P.random_unitary_operator())
    I3 = I1.compose(I2)
    np.testing.assert_array_equal(I3.matrix, lam2 @ np.asarray([[1, 1], [2, 1]]))
    x = rand_op(I1.N, gen(9))
    torch.testing.assert_close(I3.embed(x).matrix, I2.embed(I1.embed(x)).matrix, atol=1e-12, rtol=0)
    assert abs(I3.index() - np.linalg.norm(I3.matrix, 2) ** 2) < 1e-10


def test_disconnected_inclusion():
    """A disconnected inclusion (the centre of M_2 (+) M_3) has no unique Markov trace: index warns, markov raises."""
    M = TypeIAlgebra([2, 3], [2, 3], precision='double')
    Z = inc.centre(M)
    assert not Z.is_connected
    assert sorted(Z.component_indices()) == [4.0, 9.0]
    with pytest.warns(UserWarning):
        assert abs(Z.index() - 9.0) < 1e-12
    with pytest.raises(ValueError):
        Z.markov_trace()


# ======================================================================================================
# jordan: Jordan products, quadratic representation, generated subalgebras, spin factors
# ======================================================================================================
def jalg():
    return TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')


def herm_op(alg, g):
    X = torch.randn(1, alg.C, alg.k_max, alg.k_max, dtype=C128, generator=g)
    o = alg.operator((X + X.mH) / 2)
    o._is_self_adjoint = True
    return o


def test_jordan_identities_on_operators():
    """Jordan identity, P(y) = 2L(y)^2 - L(y^2), fundamental formula, triple product, Tr(P(y)x) = Tr(y^2 x)."""
    alg = jalg()
    g = gen(11)
    x, y, z = herm_op(alg, g), herm_op(alg, g), herm_op(alg, g)
    J, P = jd.jordan_product, jd.quadratic_representation
    m = lambda o: o.matrix
    x2 = J(x, x)
    torch.testing.assert_close(m(J(J(x2, y), x)), m(J(x2, J(y, x))), atol=1e-11, rtol=0)
    torch.testing.assert_close(m(P(y, x)), m(J(y, J(y, x)) * 2 - J(J(y, y), x)), atol=1e-11, rtol=0)
    torch.testing.assert_close(m(P(P(y, x), z)), m(P(y, P(x, P(y, z)))), atol=1e-10, rtol=0)
    tri = J(J(x, y), z) + J(J(y, z), x) - J(J(x, z), y)
    torch.testing.assert_close(m(jd.jordan_triple(x, y, z)), m(tri), atol=1e-11, rtol=0)
    torch.testing.assert_close(m(jd.jordan_triple(x, y, x)), m(P(x, y)), atol=1e-11, rtol=0)
    torch.testing.assert_close(P(y, x).trace, J(J(y, y), x).trace, atol=1e-11, rtol=0)
    for p in range(0, 6):
        ref = [torch.linalg.matrix_power(B, p) for B in blocks_of(x)]
        torch.testing.assert_close(blocks_of(jd.jordan_power(x, p)), ref, atol=1e-10, rtol=0)


@pytest.mark.parametrize("k,r", [(2, 1), (3, 2), (4, 2), (4, 4)])
def test_quadratic_representation_rank(k, r):
    """P(y) has rank r^2 on Herm_k(C) and r(r+1)/2 on Herm_k(R) for y of rank r; P(y)^-1 = P(y^-1) if r = k."""
    g = gen(12 + k + r)
    for real in (False, True):
        dt = F64 if real else C128
        G = torch.randn(k, r, dtype=dt, generator=g)
        y = G @ G.mH
        basis = []
        for a in range(k):
            for b in range(a, k):
                E = torch.zeros(k, k, dtype=dt)
                E[a, b] = E[b, a] = 1
                basis.append(E)
                if a != b and not real:
                    E = torch.zeros(k, k, dtype=dt)
                    E[a, b], E[b, a] = 1j, -1j
                    basis.append(E)
        img = torch.stack([(y @ E @ y).reshape(-1) for E in basis])
        img = torch.view_as_real(img).reshape(len(basis), -1) if not real else img
        assert int(torch.linalg.matrix_rank(img, rtol=1e-10)) == (r * (r + 1) // 2 if real else r * r)
        if r == k:
            yi = torch.linalg.inv(y)
            for E in basis[:3]:
                torch.testing.assert_close(yi @ (y @ E @ y) @ yi, E, atol=1e-9, rtol=0)


def test_transpose_map_is_jordan_but_not_cp():
    """The transpose preserves Jordan products and spectra, reverses products; its Choi matrix is the swap."""
    alg = jalg()
    g = gen(13)
    x, y = herm_op(alg, g), herm_op(alg, g)
    T = jd.transpose_map
    torch.testing.assert_close(T(jd.jordan_product(x, y)).matrix, jd.jordan_product(T(x), T(y)).matrix, atol=1e-12,
                               rtol=0)
    torch.testing.assert_close(T(x @ y).matrix, (T(y) @ T(x)).matrix, atol=1e-12, rtol=0)
    k = 3
    choi = torch.zeros(k * k, k * k, dtype=C128)
    for i, j in itertools.product(range(k), range(k)):
        E = torch.zeros(k, k, dtype=C128)
        E[i, j] = 1
        choi += torch.kron(E, E.T.contiguous())
    w = torch.linalg.eigvalsh(choi)
    assert abs(float(w.min()) + 1) < 1e-12 and abs(float(w.max()) - 1) < 1e-12


def test_jordan_product_of_positives_not_positive():
    """Module docstring example: det(p o q) = -1/16 for p = diag(1, 0), q = (1/2) [[1, 1], [1, 1]]."""
    alg = TypeIAlgebra([2], [2], precision='double')
    p = alg.operator(torch.tensor([[[[1, 0], [0, 0]]]], dtype=C128))
    q = alg.operator(torch.tensor([[[[0.5, 0.5], [0.5, 0.5]]]], dtype=C128))
    assert abs(complex(torch.linalg.det(jd.jordan_product(p, q).matrix[0, 0])) + 1 / 16) < 1e-15


@pytest.mark.parametrize("scale", [1e-6, 1.0, 1e6])
def test_one_generator_jordan_algebra(scale):
    """One generator: dim J = #distinct eigenvalues (unital), #distinct non-zero eigenvalues (non-unital); scale-free."""
    g = gen(14)
    U = haar(6, g)
    ev = torch.tensor([1., 2., 2., 3., 5., 5.], dtype=F64) * scale
    x = U @ torch.diag(ev).to(C128) @ U.mH
    assert jd.jordan_algebra_generated(x).dim == 4
    assert jd.jordan_algebra_generated(x, unital=False).dim == 4
    x0 = U @ torch.diag(ev - 2 * scale).to(C128) @ U.mH                         # eigenvalue 0 twice
    assert jd.jordan_algebra_generated(x0, unital=False).dim == 3
    assert jd.jordan_algebra_generated(x0).dim == 4


@pytest.mark.parametrize("n", [2, 3, 4])
def test_generic_pairs(n):
    """Two real symmetric matrices: J = Herm_n(R), dim n(n+1)/2, A_sa = Herm_n(C); complex pairs: J = Herm_n(C) for
    n >= 3, Herm_2(R) ~ V_2 for n = 2 (simultaneously real)."""
    g = gen(15 + n)
    a, b = herm(n, g, real=True), herm(n, g, real=True)
    r = jd.compare_with_star_algebra([a, b])
    assert (r['jordan_dim'], r['star_dim']) == (n * (n + 1) // 2, n * n)
    assert r['jordan'].field == 'real' and r['jordan'].is_closed()
    a, b = herm(n, g), herm(n, g)
    r = jd.compare_with_star_algebra([a, b])
    assert r['jordan_dim'] == (3 if n == 2 else n * n) and r['star_dim'] == n * n


def test_block_algebra_two_and_three_generators():
    """A = U (M_2 (x) 1_2 (+) M_3) U^*: two generic elements give V_2 (+) Herm_3 (dim 12), three give A_sa (dim 13)."""
    g = gen(16)
    U = haar(7, g)

    def elem():
        X = torch.block_diag(torch.kron(herm(2, g), torch.eye(2, dtype=C128)), herm(3, g))
        return U @ X @ U.mH
    assert jd.jordan_algebra_generated([elem(), elem()]).dim == 12
    J3 = jd.jordan_algebra_generated([elem(), elem(), elem()])
    assert J3.dim == 13 and J3.is_self_adjoint_part()
    assert jd.star_algebra_self_adjoint_part([elem(), elem()]).dim == 13


@pytest.mark.parametrize("m", [1, 2, 3, 4, 5])
def test_anticommuting_symmetries_spin_factor(m):
    """s_i o s_j = delta_ij 1; J = span{1, s_i} = V_m (dim m + 1); J = A_sa only for m in {1, 3}."""
    s = jd.anticommuting_symmetries(m)
    D = s.shape[-1]
    assert D == (2 if m == 1 else 2 ** (m // 2))
    for i, j in itertools.product(range(m), range(m)):
        torch.testing.assert_close((s[i] @ s[j] + s[j] @ s[i]) / 2, float(i == j) * torch.eye(D, dtype=C128),
                                   atol=1e-15, rtol=0)
    r = jd.compare_with_star_algebra(s)
    assert r['jordan_dim'] == m + 1
    assert r['star_dim'] == (2 if m == 1 else D * D)
    assert r['jordan'].is_self_adjoint_part() == (m in (1, 3))


def test_is_jordan_subalgebra_known_spans():
    """span{1, X, Z} = Herm_2(R) is closed; span{X, Z} is not (X o X = 1); the diagonal matrices are closed."""
    X = torch.tensor([[0., 1.], [1., 0.]], dtype=F64)
    Z = torch.tensor([[1., 0.], [0., -1.]], dtype=F64)
    assert jd.is_jordan_subalgebra([torch.eye(2, dtype=F64), X, Z])
    assert not jd.is_jordan_subalgebra([X, Z])
    assert not jd.is_jordan_subalgebra([torch.diag(torch.tensor([1., 2., 0.], dtype=F64)),
                                        torch.diag(torch.tensor([0., 1., 1.], dtype=F64))])   # squares leave the span
    E = [torch.diag(torch.eye(3, dtype=F64)[i]) for i in range(3)]
    assert jd.is_jordan_subalgebra(E)


def test_generated_algebra_contains_generators_and_is_closed():
    """The closure contains the generators and 1, has an orthonormal basis and is closed under o."""
    g = gen(17)
    gens = [herm(4, g), herm(4, g)]
    J = jd.jordan_algebra_generated(gens)
    gram = torch.einsum('aij,bji->ab', J.basis, J.basis).real
    torch.testing.assert_close(gram, torch.eye(J.dim, dtype=F64), atol=1e-12, rtol=0)
    assert bool(J.contains(torch.stack(gens + [torch.eye(4, dtype=C128)])).all())
    assert J.closure_defect('jordan', exhaustive=True) < 1e-12


class TestSpinFactor:
    m = 4

    def setup_method(self):
        self.V = jd.SpinFactor(self.m)
        self.g = gen(18)

    def test_representation_is_jordan_homomorphism(self):
        """phi(x o y) = phi(x) o phi(y), phi(P(x)y) = phi(x)phi(y)phi(x), Tr phi(x) = (D/2) tr x, from_matrix o phi = id."""
        V = self.V
        x, y = V.random(5, generator=self.g), V.random(5, generator=self.g)
        X, Y = V.to_matrix(x), V.to_matrix(y)
        torch.testing.assert_close(V.to_matrix(V.product(x, y)), (X @ Y + Y @ X) / 2, atol=1e-12, rtol=0)
        torch.testing.assert_close(V.to_matrix(V.quadratic_representation(x, y)), X @ Y @ X, atol=1e-11, rtol=0)
        D = V.representation_dim
        torch.testing.assert_close(torch.diagonal(X, dim1=-2, dim2=-1).sum(-1).real, D / 2 * V.trace(x), atol=1e-12,
                                   rtol=0)
        torch.testing.assert_close(V.from_matrix(X), x, atol=1e-12, rtol=0)
        z = V.random(5, generator=self.g)
        x2 = V.square(x)
        torch.testing.assert_close(V.product(V.product(x2, y), x), V.product(x2, V.product(y, x)), atol=1e-11, rtol=0)
        torch.testing.assert_close(V.triple(x, y, x), V.quadratic_representation(x, y), atol=1e-11, rtol=0)
        torch.testing.assert_close(V.inner(V.product(x, y), z), V.inner(y, V.product(x, z)), atol=1e-11, rtol=0)

    def test_spectral_decomposition_and_functions(self):
        """Eigenvalues a -+ |u| of phi(x) (multiplicity D/2), Jordan frame, f(x) against scipy matrix functions."""
        V = self.V
        x = V.random(4, generator=self.g)
        x[0, 1:] = 0                                                          # u = 0 (degenerate)
        lam, e = V.spectral_decomposition(x)
        X = V.to_matrix(x)
        D = V.representation_dim
        for b in range(4):
            w = torch.linalg.eigvalsh(X[b])
            torch.testing.assert_close(w, lam[b].repeat_interleave(D // 2), atol=1e-12, rtol=0)
        torch.testing.assert_close(V.product(e[:, 0], e[:, 1]), torch.zeros_like(x), atol=1e-14, rtol=0)
        torch.testing.assert_close(V.square(e[:, 1]), e[:, 1], atol=1e-14, rtol=0)
        torch.testing.assert_close(e[:, 0] + e[:, 1], V.unit((4,)), atol=1e-14, rtol=0)
        torch.testing.assert_close((lam[..., None] * e).sum(-2), x, atol=1e-13, rtol=0)
        torch.testing.assert_close(V.determinant(x), lam.prod(-1), atol=1e-12, rtol=0)
        ex = V.exp(x)
        for b in range(4):
            torch.testing.assert_close(V.to_matrix(ex[b]), torch.as_tensor(sla.expm(X[b].numpy())), atol=1e-11, rtol=1e-12)
        c = V.square(x) + V.unit((4,)) * 0.1                                 # interior of the cone
        assert bool(V.in_cone(c, interior=True).all())
        Cm = V.to_matrix(c)
        for b in range(4):
            torch.testing.assert_close(V.to_matrix(V.log(c[b])), torch.as_tensor(sla.logm(Cm[b].numpy())), atol=1e-10,
                                       rtol=0)
            torch.testing.assert_close(V.to_matrix(V.sqrt(c[b])), torch.as_tensor(sla.sqrtm(Cm[b].numpy())), atol=1e-10,
                                       rtol=0)
        torch.testing.assert_close(V.product(V.inverse(c), c), V.unit((4,)), atol=1e-11, rtol=0)

    def test_cone_is_psd_cone(self):
        """Lorentz cone a >= |u| iff phi(x) >= 0 (Koecher-Vinberg example)."""
        V = self.V
        x = V.random(200, generator=self.g)
        psd = torch.linalg.eigvalsh(V.to_matrix(x)).min(-1).values >= 0
        assert bool((V.in_cone(x) == psd).all())

    def test_v3_is_herm2(self):
        """V_3 ~ Herm_2(C): phi is onto the 4-dimensional real space Herm_2(C)."""
        V = jd.SpinFactor(3)
        basis = V.to_matrix(torch.eye(4, dtype=F64))
        assert int(torch.linalg.matrix_rank(torch.view_as_real(basis.reshape(4, -1)).reshape(4, -1))) == 4


# ======================================================================================================
# divergences
# ======================================================================================================
DALG = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')


def rand_state(alg, g, cond=None):
    mats = [wishart(k, g, cond) for k in alg.k_factors]
    tot = sum(float(torch.trace(m).real) for m in mats)
    return pos_op(alg, [m / tot for m in mats])


def scipy_renyi(rho, sig, a, kind):
    tr = sum(float(torch.trace(r).real) for r in blocks_of(rho))
    Q = 0.0
    for r, s in zip(blocks_of(rho), blocks_of(sig)):
        r, s = r.numpy(), s.numpy()
        if kind == 'petz':
            Q += np.trace(sla.fractional_matrix_power(r, a) @ sla.fractional_matrix_power(s, 1 - a)).real
        else:
            sg = sla.fractional_matrix_power(s, (1 - a) / (2 * a))
            M = sg @ r @ sg
            Q += float(np.sum(np.clip(np.linalg.eigvalsh((M + M.conj().T) / 2), 0, None) ** a))
    return math.log(Q / tr) / (a - 1)


@pytest.mark.parametrize("a", [0.3, 0.5, 0.8, 1.3, 2.0, 5.0])
def test_renyi_against_scipy(a):
    """Petz D_alpha = log(Tr rho^a sigma^{1-a})/(a-1), sandwiched ~D_alpha, both against scipy fractional powers."""
    g = gen(20)
    rho, sig = rand_state(DALG, g, cond=30), rand_state(DALG, g, cond=30)
    for kind, f in (('petz', dv.petz_renyi_divergence), ('sandwiched', dv.sandwiched_renyi_divergence)):
        assert abs(float(f(rho, sig, a)) - scipy_renyi(rho, sig, a, kind)) < 1e-10


def test_umegaki_dmax_dmin_against_scipy():
    """D = Tr rho(log rho - log sigma) (scipy logm), D_max = log lambda_max(sigma^-1/2 rho sigma^-1/2), D_0 = -log Tr Pi_rho sigma."""
    g = gen(21)
    rho, sig = rand_state(DALG, g), rand_state(DALG, g)
    ref = sum(np.trace(r.numpy() @ (sla.logm(r.numpy()) - sla.logm(s.numpy()))).real
              for r, s in zip(blocks_of(rho), blocks_of(sig)))
    assert abs(float(dv.umegaki_relative_entropy(rho, sig)) - ref) < 1e-11
    lmax = max(np.linalg.eigvalsh(sla.fractional_matrix_power(s.numpy(), -0.5) @ r.numpy()
                                  @ sla.fractional_matrix_power(s.numpy(), -0.5)).max()
               for r, s in zip(blocks_of(rho), blocks_of(sig)))
    assert abs(float(dv.max_relative_entropy(rho, sig)) - math.log(lmax)) < 1e-11
    assert abs(float(dv.min_relative_entropy(rho, sig))) < 1e-12               # faithful rho: Pi_rho = 1


def test_limits_in_alpha():
    """alpha -> 1: Umegaki (both families); alpha -> inf: D_max (sandwiched); alpha -> 0: D_0 (Petz); ~D_1/2 = -log F."""
    g = gen(22)
    rho, sig = rand_state(DALG, g), rand_state(DALG, g)
    # rank-deficient rho for a non-trivial D_0
    r1 = pos_op(DALG, [torch.ones(1, 1, dtype=C128) * 0.2] +
                [0.4 * (lambda v: torch.outer(v, v.conj()) / v.norm() ** 2)(torch.randn(k, dtype=C128, generator=g))
                 for k in (2, 3)])
    D = dv.umegaki_relative_entropy(rho, sig)
    for a in (1 - 1e-5, 1 + 1e-5):
        for f in (dv.petz_renyi_divergence, dv.sandwiched_renyi_divergence):
            assert abs(float(f(rho, sig, a) - D)) < 1e-4
    assert abs(float(dv.sandwiched_renyi_divergence(rho, sig, 1e5) - dv.max_relative_entropy(rho, sig))) < 1e-3
    assert abs(float(dv.petz_renyi_divergence(r1, sig, 1e-7) - dv.min_relative_entropy(r1, sig))) < 1e-5
    assert abs(float(dv.sandwiched_renyi_divergence(rho, sig, 0.5) + torch.log(fidelity(rho, sig)))) < 1e-12
    assert float(dv.renyi_divergence(r1, sig, 0, 'petz')) == float(dv.min_relative_entropy(r1, sig))


@pytest.mark.parametrize("c,d", [(1e-6, 1e6), (1e6, 1e-6), (1e-6, 1e-6), (1e6, 1e6)])
def test_scaling_covariance(c, d):
    """D(c rho || d sigma) = D(rho || sigma) + log(c/d) for every divergence (relative thresholds)."""
    g = gen(23)
    rho, sig = rand_state(DALG, g), rand_state(DALG, g)
    rc = pos_op(DALG, [b * c for b in blocks_of(rho)])
    sd = pos_op(DALG, [b * d for b in blocks_of(sig)])
    shift = math.log(c / d)
    for a in (0.4, 2.5):
        for f in (dv.petz_renyi_divergence, dv.sandwiched_renyi_divergence):
            assert abs(float(f(rc, sd, a) - f(rho, sig, a)) - shift) < 1e-11
    for f in (dv.umegaki_relative_entropy, dv.max_relative_entropy, dv.min_relative_entropy):
        assert abs(float(f(rc, sd) - f(rho, sig)) - shift) < 1e-11


def test_support_conventions_commuting():
    """Exact commuting cases: supp rho in supp sigma, supp rho not in supp sigma (alpha > 1: inf, alpha < 1: finite)."""
    alg = TypeIAlgebra([3], [3], precision='double')
    U = haar(3, gen(24))
    mk = lambda v: pos_op(alg, [U @ torch.diag(torch.tensor(v, dtype=C128)) @ U.mH])
    rho, sig, sig_def = mk([0.5, 0.5, 0.]), mk([0.2, 0.3, 0.5]), mk([0.4, 0., 0.6])
    Q = math.sqrt(.5 * .2) + math.sqrt(.5 * .3)
    assert abs(float(dv.petz_renyi_divergence(rho, sig, 0.5)) + 2 * math.log(Q)) < 1e-12
    assert abs(float(dv.petz_renyi_divergence(rho, sig_def, 0.5)) + 2 * math.log(math.sqrt(.5 * .4))) < 1e-12
    assert abs(float(dv.min_relative_entropy(rho, sig_def)) + math.log(0.4)) < 1e-12
    for f in (dv.petz_renyi_divergence, dv.sandwiched_renyi_divergence):
        assert math.isinf(float(f(rho, sig_def, 2.0)))
    assert math.isinf(float(dv.umegaki_relative_entropy(rho, sig_def)))
    assert math.isinf(float(dv.max_relative_entropy(rho, sig_def)))


@pytest.mark.parametrize("seed", range(4))
def test_regression_orthogonal_supports_infinite(seed):
    """Regression (divergences support convention): rho, sigma with orthogonal supports in a rotated basis give +inf
    for every alpha (Q_alpha = 0) and for D_0; the Petz pair sum used to return ~ -log(eps^2)/(1-alpha) (e.g. 142 at
    alpha = 1/2) and D_0 ~ 37 from rounding-level overlaps."""
    alg = TypeIAlgebra([2, 4], [2, 4], precision='double')
    g = gen(30 + seed)
    U1, U2 = haar(2, g), haar(4, g)
    rho = pos_op(alg, [0.3 * torch.eye(2, dtype=C128) / 2,
                       U2 @ torch.diag(torch.tensor([0.3, 0.4, 0., 0.], dtype=C128)) @ U2.mH])
    sig = pos_op(alg, [torch.zeros(2, 2, dtype=C128),
                       U2 @ torch.diag(torch.tensor([0., 0., 0.4, 0.6], dtype=C128)) @ U2.mH])
    for a in (0.1, 0.5, 0.9, 2.0):
        for f in (dv.petz_renyi_divergence, dv.sandwiched_renyi_divergence):
            assert math.isinf(float(f(rho, sig, a)))
    assert math.isinf(float(dv.min_relative_entropy(rho, sig)))
    sc = dv.sector_contributions(rho, sig, 0.5, 'petz')
    assert bool(torch.isinf(sc.conditional).all()) and float(sc.terms.abs().max()) == 0.0
    # a genuine (non-rounding) overlap stays finite: rho partly inside supp sigma
    rho2 = pos_op(alg, [0.3 * torch.eye(2, dtype=C128) / 2,
                        U2 @ torch.diag(torch.tensor([0.3, 0.3, 1e-6, 0.], dtype=C128)) @ U2.mH])
    v = float(dv.petz_renyi_divergence(rho2, sig, 0.5))
    ref = -2 * math.log(math.sqrt(1e-6 * 0.4) / (0.9 + 1e-6))
    assert abs(v - ref) < 1e-8
    assert abs(float(dv.min_relative_entropy(rho2, sig)) + math.log(0.4 / (0.9 + 1e-6))) < 1e-10


def test_data_processing_and_order_relations():
    """DPI (sandwiched alpha >= 1/2, Petz alpha in (0, 2]), ALT ~D <= D, monotonicity in alpha, positivity."""
    torch.manual_seed(25)
    A = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')
    Bo = TypeIAlgebra([2, 2], [2, 2], precision='double')
    rho, sig = random_density_matrix(A, 32), random_density_matrix(A, 32)
    ch = random_inter_sector_channel(A, Bo, kraus_rank=2)
    rp, sp = ch.apply(rho), ch.apply(sig)
    for a in (0.5, 0.8, 1.5, 2.0, 3.0):
        assert float((dv.sandwiched_renyi_divergence(rp, sp, a) - dv.sandwiched_renyi_divergence(rho, sig, a)).max()) < 1e-12
    for a in (0.2, 0.5, 0.9, 1.5, 2.0):
        assert float((dv.petz_renyi_divergence(rp, sp, a) - dv.petz_renyi_divergence(rho, sig, a)).max()) < 1e-12
    assert float((dv.umegaki_relative_entropy(rp, sp) - dv.umegaki_relative_entropy(rho, sig)).max()) < 1e-12
    for a in (0.3, 0.7, 1.5, 4.0):
        assert float((dv.sandwiched_renyi_divergence(rho, sig, a) - dv.petz_renyi_divergence(rho, sig, a)).max()) < 1e-12
    al = [0.1, 0.3, 0.5, 0.7, 0.99, 1.0, 1.01, 1.5, 2, 5, 20]
    for kind in ('petz', 'sandwiched'):
        vals = torch.stack([dv.renyi_divergence(rho, sig, a, kind) for a in al])
        assert float((vals[1:] - vals[:-1]).min()) > -1e-12
        assert float(vals.min()) > 0
    assert float(dv.sandwiched_renyi_divergence(rho, rho, 2.0).abs().max()) < 1e-12
    assert float(dv.petz_renyi_divergence(rho, rho, 0.5).abs().max()) < 1e-12


def test_sector_contributions_chain_rule():
    """D_a = (a-1)^-1 log sum_c p_c^a q_c^{1-a} exp((a-1) D_a(rho_c/p_c||sigma_c/q_c)); D = D(p||q) + sum p_c D_c."""
    torch.manual_seed(26)
    rho, sig = random_density_matrix(DALG, 8), random_density_matrix(DALG, 8)
    for kind in ('petz', 'sandwiched'):
        for a in (0.6, 1.7):
            sc = dv.sector_contributions(rho, sig, a, kind)
            Dc = torch.log((sc.p ** a * sc.q ** (1 - a) * torch.exp((a - 1) * sc.conditional)).sum(-1)) / (a - 1)
            torch.testing.assert_close(Dc, dv.renyi_divergence(rho, sig, a, kind), atol=1e-12, rtol=0)
    su_ = dv.sector_contributions(rho, sig, kind='umegaki')
    Du = (su_.p * torch.log(su_.p / su_.q)).sum(-1) + (su_.p * su_.conditional).sum(-1)
    torch.testing.assert_close(Du, dv.umegaki_relative_entropy(rho, sig), atol=1e-12, rtol=0)


def test_trace_independence():
    """Densities w.r.t. tau_w = sum w_c Tr_c (d = rho_c / w_c) give the same divergences ('norm', 'tau_vN')."""
    torch.manual_seed(27)
    rho, sig = random_density_matrix(DALG, 4), random_density_matrix(DALG, 4)
    for trace, w in (('norm', [1, 1 / 2, 1 / 3]), ('tau_vN', [1 / 3, 1 / 6, 1 / 9])):
        wt = torch.tensor(w, dtype=F64)[None, :, None, None]
        dr, ds = DALG.operator(rho.matrix / wt), DALG.operator(sig.matrix / wt)
        for o in (dr, ds):
            o._is_positive = o._is_self_adjoint = True
        for a in (0.5, 2.0):
            torch.testing.assert_close(dv.sandwiched_renyi_divergence(dr, ds, a, trace=trace),
                                       dv.sandwiched_renyi_divergence(rho, sig, a), atol=1e-12, rtol=0)
            torch.testing.assert_close(dv.petz_renyi_divergence(dr, ds, a, trace=trace),
                                       dv.petz_renyi_divergence(rho, sig, a), atol=1e-12, rtol=0)


def test_beigi_kosaki_and_lp_norms():
    """~D_a = a/(a-1) log ||sigma^-1/2 rho sigma^-1/2||_{a,sigma} (Beigi); Kosaki norms <= ||x||_inf, increasing in p;
    ||A||_2 = Frobenius norm, ||A||_1 = trace norm, weights of 'norm'."""
    torch.manual_seed(28)
    rho, sig = random_density_matrix(DALG, 4), random_density_matrix(DALG, 4)
    X = sig.power(-0.5) @ rho @ sig.power(-0.5)
    for a in (1.5, 3.0):
        torch.testing.assert_close(dv.sandwiched_renyi_divergence(rho, sig, a),
                                   a / (a - 1) * torch.log(dv.kosaki_norm(X, sig, a)), atol=1e-12, rtol=0)
    g = gen(28)
    x = rand_op(DALG, g, 4)
    norms = torch.stack([dv.kosaki_norm(x, rho, p) for p in (1, 1.5, 2, 4, 10)])
    assert float((norms[1:] - norms[:-1]).min()) > -1e-12
    assert float((norms[-1] - dv.lp_norm(x, math.inf)).max()) < 1e-12
    torch.testing.assert_close(dv.lp_norm(x, 2), torch.linalg.vector_norm(x.matrix, dim=(-3, -2, -1)), atol=1e-12,
                               rtol=0)
    tn = sum(torch.linalg.svdvals(x.matrix[:, c, :k, :k]).sum(-1) for c, k in enumerate(DALG.k_factors))
    torch.testing.assert_close(dv.lp_norm(x, 1), tn, atol=1e-12, rtol=0)
    tn_norm = sum(torch.linalg.svdvals(x.matrix[:, c, :k, :k]).sum(-1) / k for c, k in enumerate(DALG.k_factors))
    torch.testing.assert_close(dv.lp_norm(x, 1, trace='norm'), tn_norm, atol=1e-12, rtol=0)


def test_petz_recovery_and_universal_recovery_bound():
    """R(Phi(sigma)) = sigma for CPTP Phi; Tr R(X) = Tr(Pi X); rotated maps too; JRSWW:
    D(rho||sigma) - D(Phi rho||Phi sigma) >= -log F(rho, R_univ(Phi rho))."""
    torch.manual_seed(29)
    A = TypeIAlgebra([2, 3], [2, 3], precision='double')
    Bo = TypeIAlgebra([2, 2], [2, 2], precision='double')
    rho, sig = random_density_matrix(A, 6), random_density_matrix(A, 6)
    ch = random_inter_sector_channel(A, Bo, kraus_rank=2)
    rp, sp = ch.apply(rho), ch.apply(sig)
    R = dv.PetzRecovery(sig, ch)
    torch.testing.assert_close(R(sp).matrix, sig.matrix, atol=1e-12, rtol=0)
    torch.testing.assert_close(R.rotated(sp, 0.7).matrix, sig.matrix, atol=1e-12, rtol=0)
    torch.testing.assert_close(R(rp).trace, rp.trace, atol=1e-12, rtol=0)
    lhs = dv.umegaki_relative_entropy(rho, sig) - dv.umegaki_relative_entropy(rp, sp)
    rhs = -torch.log(fidelity(rho, dv.universal_recovery(sig, ch, rp)))
    assert float((rhs - lhs).max()) < 1e-10
    ts, ws = dv.universal_recovery_nodes()
    integral = sum(0.5 * math.pi / (math.cosh(math.pi * t) + 1) for t in ts) * (ts[1] - ts[0])
    assert abs(integral - 1) < 1e-12 and abs(sum(ws) - 1) < 1e-15


# ======================================================================================================
# su2: Clebsch-Gordan, multiplets, total-spin decomposition of spin chains
# ======================================================================================================
def test_clebsch_gordan_against_sympy_all_small():
    """Racah's formula with Condon-Shortley phases equals sympy.physics.quantum.cg.CG for all j1, j2 <= 2."""
    sympy_cg = pytest.importorskip("sympy.physics.quantum.cg")
    from sympy import Rational
    for a, c in itertools.product(range(5), range(5)):
        for e in range(abs(a - c), a + c + 1, 2):
            for b in range(-a, a + 1, 2):
                for d in range(-c, c + 1, 2):
                    f = b + d
                    if abs(f) > e:
                        continue
                    ref = float(sympy_cg.CG(Rational(a, 2), Rational(b, 2), Rational(c, 2), Rational(d, 2),
                                            Rational(e, 2), Rational(f, 2)).doit())
                    val = su.clebsch_gordan(Fraction(a, 2), Fraction(b, 2), Fraction(c, 2), Fraction(d, 2),
                                            Fraction(e, 2), Fraction(f, 2))
                    assert abs(val - ref) < 1e-14


@pytest.mark.parametrize("j1,j2", [(Fraction(1, 2), Fraction(1, 2)), (1, Fraction(3, 2)), (2, 1)])
def test_clebsch_gordan_unitary_exact(j1, j2):
    """The CG matrix (m1, m2) -> (J, M) is real orthogonal: exact rational check of both orthogonality relations."""
    a, c = int(2 * j1), int(2 * j2)
    rows = [(b, d) for b in range(-a, a + 1, 2) for d in range(-c, c + 1, 2)]
    cols = [(e, f) for e in range(abs(a - c), a + c + 1, 2) for f in range(-e, e + 1, 2)]
    assert len(rows) == len(cols)
    M = [[su.clebsch_gordan_exact(Fraction(a, 2), Fraction(b, 2), Fraction(c, 2), Fraction(d, 2), Fraction(e, 2),
                                  Fraction(f, 2)) for (e, f) in cols] for (b, d) in rows]
    # signs times sqrt(squares): products of two entries are rational only up to sqrt; check numerically at 1e-15
    X = np.array([[s * math.sqrt(float(q)) for s, q in r] for r in M])
    np.testing.assert_allclose(X @ X.T, np.eye(len(rows)), atol=1e-14)
    np.testing.assert_allclose(X.T @ X, np.eye(len(rows)), atol=1e-14)
    # squares in each column sum to exactly 1
    for j in range(len(cols)):
        assert sum(M[i][j][1] for i in range(len(rows))) == 1


def test_scalar_coefficient_closed_form():
    """c(S, S'') = sum_q (-1)^q CG CG / sqrt((2S+1)(2S''+1)) = (-1)^{S - S''} / (2S + 1) (S > 0 or S'' > 0)."""
    for t in range(0, 12):
        for tpp in (t - 2, t, t + 2):
            if tpp < 0 or (t == 0 and tpp == 0):
                continue
            c = su.SU2Chain._scalar_coefficient(t, tpp)
            assert abs(c - (-1) ** ((t - tpp) // 2) / (t + 1)) < 1e-14


def hook_count(L, S2):
    """Standard Young tableaux of shape (L/2 + S, L/2 - S) by the hook length formula."""
    r1, r2 = (L + S2) // 2, (L - S2) // 2
    hooks = 1
    for j in range(r1):
        hooks *= (r1 - j) + (1 if j < r2 else 0)
    for j in range(r2):
        hooks *= (r2 - j)
    return math.factorial(L) // hooks


@pytest.mark.parametrize("L", [1, 2, 5, 8, 11])
def test_multiplicities(L):
    """m_S = binom(L, L/2 - S) - binom(L, L/2 - S - 1) = #SYT of shape (L/2 + S, L/2 - S); sum (2S+1) m_S = 2^L."""
    chain = su.SU2Chain(L)
    assert sum((t + 1) * m for t, m in zip(chain.two_spins, chain.multiplicities)) == 2 ** L
    for t, m in zip(chain.two_spins, chain.multiplicities):
        assert m == hook_count(L, t)


@pytest.mark.parametrize("L", [3, 4, 6])
def test_highest_weight_vectors_dense(L):
    """B_S is orthonormal, in ker S^+, has S^z = S and S^2 = S(S+1) (dense Kronecker spin operators)."""
    chain = su.SU2Chain(L)
    ops = label_spin_ops(L).sum(0)                           # total S_x, S_y, S_z in the label basis
    Sp = ops[0] + 1j * ops[1]
    S2 = sum(o @ o for o in ops)
    for t in chain.two_spins:
        N = (L + t) // 2
        lab = chain.u1.labels[N].numpy()
        B = chain.hw_basis(Fraction(t, 2)).numpy()
        X = np.zeros((2 ** L, B.shape[1]), dtype=complex)
        X[lab] = B
        np.testing.assert_allclose(B.T @ B, np.eye(B.shape[1]), atol=1e-13)
        np.testing.assert_allclose(Sp @ X, 0, atol=1e-13)
        np.testing.assert_allclose(ops[2] @ X, t / 2 * X, atol=1e-13)
        np.testing.assert_allclose(S2 @ X, t / 2 * (t / 2 + 1) * X, atol=1e-12)
        for tm in range(-t, t + 1, 2):                          # the multiplets, Condon-Shortley lowering
            Y = chain.multiplet_basis(Fraction(t, 2), Fraction(tm, 2)).numpy()
            Xm = np.zeros((2 ** L, Y.shape[1]), dtype=complex)
            Xm[chain.u1.labels[(L + tm) // 2].numpy()] = Y
            np.testing.assert_allclose(Y.T @ Y, np.eye(Y.shape[1]), atol=1e-12)
            np.testing.assert_allclose(S2 @ Xm, t / 2 * (t / 2 + 1) * Xm, atol=1e-12)
            if tm < t:
                Yu = chain.multiplet_basis(Fraction(t, 2), Fraction(tm + 2, 2)).numpy()
                Xu = np.zeros((2 ** L, Yu.shape[1]), dtype=complex)
                Xu[chain.u1.labels[(L + tm + 2) // 2].numpy()] = Yu
                c = math.sqrt((t / 2 - tm / 2) * (t / 2 + tm / 2 + 1))     # S^+ |S M> = c |S M+1>
                np.testing.assert_allclose(Sp @ Xm, c * Xu, atol=1e-12)


@pytest.mark.parametrize("L,boundary,J2,h", [(4, 'open', 0.0, 0.0), (5, 'periodic', 0.3, 0.0), (6, 'periodic', 0.0, 0.4),
                                             (6, 'open', -0.5, -0.2)])
def test_full_spectrum_against_dense_hamiltonian(L, boundary, J2, h):
    """H = J sum S_i.S_j + J2 sum_nnn S_i.S_j + h S^z: SU(2) block spectra with multiplets E + h(M - S) equal the
    eigenvalues of the dense 2^L Hamiltonian; ground state spin and S^z as in Lieb-Mattis / Zeeman."""
    chain = su.SU2Chain(L, boundary=boundary)
    ops = label_spin_ops(L)
    H = np.zeros((2 ** L, 2 ** L), dtype=complex)
    for (i, j), w in [(b, 1.0) for b in chain.bonds] + [(b, J2) for b in chain.next_nearest_bonds]:
        H += w * sum(ops[i, a] @ ops[j, a] for a in range(3))
    H += h * ops[:, 2].sum(0)
    Hs = chain.heisenberg(J=1.0, J2=J2, h=h)
    E = chain.full_spectrum(Hs)[0].numpy()
    np.testing.assert_allclose(E, np.linalg.eigvalsh(H), atol=1e-11)
    gs = chain.ground_state(Hs)
    assert abs(float(gs.energy[0]) - E[0]) < 1e-11
    psi = chain.embed(gs.spin[0], gs.coefficients[0, :chain.multiplicity(gs.spin[0])], M=gs.sz[0])
    N = int(round(L / 2 + float(gs.sz[0])))
    v = np.zeros(2 ** L, dtype=complex)
    v[chain.u1.labels[N].numpy()] = psi.numpy()
    np.testing.assert_allclose(H @ v, E[0] * v, atol=1e-10)
    if J2 == 0 and h == 0 and L % 2 == 0:
        assert float(gs.spin[0]) == 0.0                           # Lieb-Mattis (bipartite, even L)


@pytest.mark.parametrize("L", [4, 5])
def test_wigner_eckart_reconstructs_all_matrix_elements(L):
    """<S'M'a|T_q|SMb> = <S M; 1 q|S' M'> <S' a||S_i||S b> / sqrt(2S'+1) for all M, M', q, sites (dense check);
    hermiticity R_{S,S'} = (-1)^{S'-S} R_{S',S}^T; Edmonds sum rule."""
    chain = su.SU2Chain(L)
    R = chain.reduced_matrix_elements()
    ops = label_spin_ops(L)
    Tq = {1: -(ops[:, 0] + 1j * ops[:, 1]) / math.sqrt(2), 0: ops[:, 2], -1: (ops[:, 0] - 1j * ops[:, 1]) / math.sqrt(2)}

    def full(t, tm):
        Y = chain.multiplet_basis(Fraction(t, 2), Fraction(tm, 2)).numpy()
        X = np.zeros((2 ** L, Y.shape[1]), dtype=complex)
        X[chain.u1.labels[(L + tm) // 2].numpy()] = Y
        return X
    for (tp, t), red in R.items():
        red = red.numpy()
        assert np.allclose(red.transpose(0, 2, 1) * (-1) ** ((tp - t) // 2), R[(t, tp)].numpy(), atol=1e-12)
        total = 0.0
        for tm in range(-t, t + 1, 2):
            for q in (-1, 0, 1):
                tmp = tm + 2 * q
                if abs(tmp) > tp:
                    continue
                cg = su.clebsch_gordan(Fraction(t, 2), Fraction(tm, 2), 1, q, Fraction(tp, 2), Fraction(tmp, 2))
                A = np.einsum('ai,sab,bj->sij', full(tp, tmp).conj(), Tq[q], full(t, tm))
                np.testing.assert_allclose(A, red * cg / math.sqrt(tp + 1), atol=1e-12)
                total = total + np.abs(A) ** 2
        if t + 1 > 0 and tp + 1 > 0 and np.size(total):
            np.testing.assert_allclose(total, np.abs(red) ** 2, atol=1e-12)


@pytest.mark.parametrize("L", [4, 5])
def test_spin_correlations_and_casimir(L):
    """S_i.S_j by Wigner-Eckart = by direct projection = dense; sum_{i<j} S_i.S_j = (S(S+1) - 3L/4)/2 per block;
    lift o restrict = id on SU(2)'."""
    chain = su.SU2Chain(L)
    pairs = [(i, j) for i in range(L) for j in range(L)]
    Cw = chain.spin_correlation(pairs)
    Cd = chain.spin_correlation(pairs, method='direct')
    torch.testing.assert_close(Cw.matrix, Cd.matrix, atol=1e-12, rtol=0)
    ops = label_spin_ops(L)
    for t, (c, m) in zip(chain.two_spins, enumerate(chain.multiplicities)):
        N = (L + t) // 2
        lab = chain.u1.labels[N].numpy()
        B = chain.hw_basis(Fraction(t, 2)).numpy()
        for p, (i, j) in enumerate(pairs):
            Hd = sum(ops[i, a] @ ops[j, a] for a in range(3))[np.ix_(lab, lab)]
            np.testing.assert_allclose(Cw.matrix[p, c, :m, :m].numpy(), B.T @ Hd @ B, atol=1e-12)
        tot = sum(Cw.matrix[p, c, :m, :m] for p, (i, j) in enumerate(pairs) if i < j)
        S = t / 2
        torch.testing.assert_close(tot, (S * (S + 1) - 0.75 * L) / 2 * torch.eye(m, dtype=tot.dtype), atol=1e-12,
                                   rtol=0)
    H = chain.heisenberg(J=1.0, J2=0.4)
    back = chain.restrict(chain.lift(H))
    torch.testing.assert_close(back.matrix, H.matrix, atol=1e-12, rtol=0)
    torch.testing.assert_close(chain.casimir().matrix,
                               chain.algebra.central([s * (s + 1) for s in chain.spins]).matrix, atol=0, rtol=0)
