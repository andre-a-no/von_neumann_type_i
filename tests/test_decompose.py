"""
Wedderburn decomposition (torch_vn_algebra.decompose) against exact results: the S_3 example of the course,
Schur-Weyl duality for S_L x SU(2) on (C^2)^{(x)L}, the regular representation of Z_n, and random
block-diagonal *-algebras hidden by a Haar unitary, plus robustness of the randomised algorithm.
"""
import contextlib
import itertools
import math
import warnings

import pytest
import torch

import torch_vn_algebra.decompose as dc
from torch_vn_algebra import TypeIAlgebra

C128 = torch.complex128
TOL = 1e-10                      # residuals asserted in double precision


# ----------------------------------------------------------------------------
# helpers (independent of the module under test)
# ----------------------------------------------------------------------------

def haar(n, gen, dtype=C128):
    Z = torch.randn(n, n, dtype=dtype, generator=gen)
    Q, R = torch.linalg.qr(Z)
    return Q * torch.diagonal(R).sgn()


def hidden_algebra(ks, ms, seed, n_gens=3, dtype=C128, kind='random'):
    """
    Random elements of W ((+)_c M_{k_c} (x) 1_{m_c}) W^* with a Haar unitary W. kind='random': Gaussian blocks
    (they generate the algebra); kind='nilpotent': strictly upper triangular blocks plus distinct diagonals
    (a non-*-closed generating set; only with the adjoints does it generate the algebra).
    Returns (generators (n_gens, n, n), W, list of the true sector column ranges).
    """
    g = torch.Generator().manual_seed(seed)
    n = sum(k * m for k, m in zip(ks, ms))
    W = haar(n, g, dtype)
    gens = []
    for _ in range(n_gens):
        blocks = []
        for k, m in zip(ks, ms):
            x = torch.randn(k, k, dtype=dtype, generator=g)
            if kind == 'nilpotent':
                x = torch.triu(x, 1) + torch.diag(torch.randn(k, generator=g).to(dtype))
            blocks.append(torch.kron(x, torch.eye(m, dtype=dtype)))
        gens.append(W @ torch.block_diag(*blocks) @ W.mH)
    ranges, o = [], 0
    for k, m in zip(ks, ms):
        ranges.append((o, o + k * m))
        o += k * m
    return torch.stack(gens), W, ranges


def random_element(ks, ms, W, seed, side='algebra'):
    """A random element of the hidden algebra (side='algebra') or of its commutant."""
    g = torch.Generator().manual_seed(seed)
    blocks = []
    for k, m in zip(ks, ms):
        if side == 'algebra':
            blocks.append(torch.kron(torch.randn(k, k, dtype=W.dtype, generator=g), torch.eye(m, dtype=W.dtype)))
        else:
            blocks.append(torch.kron(torch.eye(k, dtype=W.dtype), torch.randn(m, m, dtype=W.dtype, generator=g)))
    return W @ torch.block_diag(*blocks) @ W.mH


def span_dim(mats, tol=1e-9):
    M = mats.reshape(mats.shape[0], -1)
    s = torch.linalg.svdvals(M)
    return int((s > tol * s[0]).sum())


def in_span(X, basis):
    """Relative distance of the matrices X (B, n, n) from the span of an orthonormal basis (d, n, n)."""
    Xf, Bf = X.reshape(X.shape[0], -1), basis.reshape(basis.shape[0], -1)
    P = (Xf @ Bf.mH) @ Bf
    return float((torch.linalg.vector_norm(Xf - P, dim=-1) / torch.linalg.vector_norm(Xf, dim=-1)).max())


def block_structure_residual(D, a, side='algebra'):
    """
    Independent check of V^* a V = (+)_c a_c (x) 1_{m_c} (side 'algebra') or (+)_c 1_{k_c} (x) x_c (side
    'commutant'): off-diagonal sector blocks vanish and, within a sector, the copies are identical and do not
    mix. Relative Frobenius residual, without using to_algebra.
    """
    b = D.V.mH @ a @ D.V
    expected = torch.zeros_like(b)
    o = 0
    for k, m in zip(D.k_factors, D.multiplicities):
        blk = b[o:o + k * m, o:o + k * m].reshape(k, m, k, m)
        if side == 'algebra':
            core = torch.einsum('ijkj->ik', blk) / m                  # average over the copies
            exp = torch.kron(core, torch.eye(m, dtype=b.dtype))
        else:
            core = torch.einsum('ijil->jl', blk) / k
            exp = torch.kron(torch.eye(k, dtype=b.dtype), core)
        expected[o:o + k * m, o:o + k * m] = exp
        o += k * m
    return float(torch.linalg.matrix_norm(b - expected) / torch.linalg.matrix_norm(a))


def km(D):
    return sorted(zip(D.k_factors, D.multiplicities))


def frob(x):
    return float(torch.linalg.matrix_norm(x))


# ----------------------------------------------------------------------------
# constructors
# ----------------------------------------------------------------------------

def test_permutation_representation_is_a_homomorphism():
    L, d = 3, 3
    perms = list(itertools.permutations(range(L)))
    P = dc.permutation_representation(L, d, perms)
    assert P.shape == (6, 27, 27)
    idx = {p: i for i, p in enumerate(perms)}
    for s, t in itertools.product(perms, perms):
        st = tuple(s[t[i]] for i in range(L))                       # (s t)(i) = s(t(i))
        assert torch.equal(P[idx[s]] @ P[idx[t]], P[idx[st]])
    # the vector in slot i moves to slot sigma(i): sigma = (0 -> 1 -> 2 -> 0)
    e = [torch.eye(d, dtype=C128)[i] for i in range(d)]
    v = torch.kron(torch.kron(e[0], e[1]), e[2])
    s = (1, 2, 0)
    w = torch.kron(torch.kron(e[2], e[0]), e[1])                   # slot 1 <- e0, slot 2 <- e1, slot 0 <- e2
    assert torch.equal(P[idx[s]] @ v, w)
    adj = dc.permutation_representation(4, 2, 'adjacent')
    assert adj.shape == (3, 16, 16)
    with pytest.raises(ValueError):
        dc.permutation_representation(3, 2, [(0, 0, 1)])


def test_total_spin_commutation_relations():
    for spin, L in ((0.5, 3), (1.0, 2)):
        Sx, Sy, Sz = dc.total_spin(L, spin)
        assert frob(Sx @ Sy - Sy @ Sx - 1j * Sz) < 1e-12
        S2 = Sx @ Sx + Sy @ Sy + Sz @ Sz
        jmax = L * spin
        w = torch.linalg.eigvalsh(S2)
        assert abs(float(w.max()) - jmax * (jmax + 1)) < 1e-10
    sx, sy, sz = dc.spin_matrices(1.5)
    assert frob(sx @ sx + sy @ sy + sz @ sz - 1.5 * 2.5 * torch.eye(4, dtype=C128)) < 1e-12


# ----------------------------------------------------------------------------
# S_3 on three qubits (course, Sec. 2.4)
# ----------------------------------------------------------------------------

@pytest.mark.parametrize('route', ['group', 'generators'])
def test_s3_on_three_qubits(route):
    G = dc.permutation_representation(3, 2)
    D = dc.from_group(G) if route == 'group' else dc.from_generators(dc.permutation_representation(3, 2, 'adjacent'))
    assert D.k_factors == (1, 2) and D.multiplicities == (4, 2)    # pi(S3)'' = M_1 (+) M_2, multiplicities 4, 2
    assert sum(k * m for k, m in zip(D.k_factors, D.multiplicities)) == 8
    assert D.dims == (5, 20)
    assert D.algebra.k_factors == [1, 2] and D.commutant_algebra().k_factors == [4, 2]
    # the k = 1 sector is the symmetric subspace: its central projection is the symmetriser
    P = D.central_projections()
    assert frob(P[0] - G.mean(0)) < TOL
    # the blocks rho_c(g) form a representation: rho(g h) = rho(g) rho(h); rho_1 is trivial
    A = D.to_algebra(G).matrix
    perms = list(itertools.permutations(range(3)))
    idx = {p: i for i, p in enumerate(perms)}
    for s, t in itertools.product(perms, perms):
        st = tuple(s[t[i]] for i in range(3))
        assert (A[idx[st]] - A[idx[s]] @ A[idx[t]]).abs().max() < TOL
    assert (A[:, 0, 0, 0] - 1).abs().max() < TOL
    # character of the 2-dim irrep: chi(e) = 2, chi(transposition) = 0, chi(3-cycle) = -1
    chi = torch.diagonal(A[:, 1, :2, :2], dim1=-2, dim2=-1).sum(-1).real
    expect = {(0, 1, 2): 2.0, (1, 0, 2): 0.0, (0, 2, 1): 0.0, (2, 1, 0): 0.0, (1, 2, 0): -1.0, (2, 0, 1): -1.0}
    for p, i in idx.items():
        assert abs(float(chi[i]) - expect[p]) < TOL
    res = D.verify()
    assert res['ok'] and res['algebra'] < TOL and res['commutant'] < TOL and res['commutant_dim'] == (20, 20)


def test_commutant_basis_methods_agree_s3():
    G = dc.permutation_representation(3, 2)
    bases = {m: dc.commutant_basis(G, group=True, method=m) for m in ('twirl', 'nullspace', 'auto')}
    bases['generators'] = dc.commutant_basis(dc.permutation_representation(3, 2, 'adjacent'))
    assert dc.commutant_dimension(G) == 20
    for name, B in bases.items():
        assert B.shape == (20, 8, 8), name
        gram = torch.einsum('aij,bij->ab', B.conj(), B)
        assert (gram - torch.eye(20, dtype=C128)).abs().max() < TOL, name           # orthonormal (HS)
        assert max(frob(g @ X - X @ g) for g in G for X in B) < TOL, name           # commutes with the group
        assert in_span(B, bases['twirl']) < TOL, name                                # same space
    # the twirl is the orthogonal projection onto the commutant
    X = torch.randn(8, 8, dtype=C128, generator=torch.Generator().manual_seed(3))
    PX = dc.twirl(G, X)
    Bf = bases['twirl'].reshape(20, -1)
    coeff = Bf.conj() @ X.reshape(-1)                                  # <B_a, X> = tr(B_a^* X)
    assert frob(PX - (coeff @ Bf).reshape(8, 8)) < TOL
    # double commutant: (pi(S3)')' = pi(S3)'' has dimension 5 and contains the group
    B2 = dc.commutant_basis(bases['twirl'])
    assert B2.shape[0] == 5 and in_span(G, B2) < TOL


# ----------------------------------------------------------------------------
# Schur-Weyl duality: S_L and SU(2) on (C^2)^{(x)L}
# ----------------------------------------------------------------------------

def schur_weyl_pairs(L):
    """(dim of the S_L irrep (L-j, j), dim of the SU(2) irrep of spin L/2 - j), j <= L/2."""
    out = []
    for j in range(L // 2 + 1):
        f = math.comb(L, j) - (math.comb(L, j - 1) if j else 0)        # hook length formula, two rows
        out.append((f, L - 2 * j + 1))
    return sorted(out)


@pytest.mark.parametrize('L', [2, 3, 4, 5, 6])
def test_schur_weyl(L):
    perm = dc.from_generators(dc.permutation_representation(L, 2, 'adjacent'))
    C = dc.commutant_basis(dc.total_spin(L))                           # SU(2)' by the null space
    su2 = dc.from_generators(dc.total_spin(L), commutant=C)
    pairs = schur_weyl_pairs(L)
    assert km(perm) == pairs
    assert km(su2) == sorted((m, k) for k, m in pairs)                 # roles of k and m swapped
    assert perm.dims == su2.dims[::-1]
    for D in (perm, su2):
        r = D.diagnostics['residuals']
        assert r['algebra'] < TOL and r['commutant'] < TOL and r['unitarity'] < TOL
    # the commutant of SU(2) is the algebra generated by the permutations: same dimension, permutations inside
    assert C.shape[0] == perm.dim
    assert in_span(dc.permutation_representation(L, 2), C) < TOL
    # the decompositions share the centre: the same set of central projections
    Pp, Ps = perm.central_projections(), su2.central_projections()
    for P in Pp:
        assert min(frob(P - Q) for Q in Ps) < TOL
    # the permutation operators lie in the commutant of the SU(2) decomposition and vice versa
    assert block_structure_residual(su2, dc.permutation_representation(L, 2, 'adjacent')[0], 'commutant') < TOL
    for S in dc.total_spin(L):
        assert block_structure_residual(perm, S, 'commutant') < TOL


@pytest.mark.parametrize('L', [2, 3, 4, 5])
def test_schur_weyl_group_route(L):
    D = dc.from_group(dc.permutation_representation(L, 2))
    assert km(D) == schur_weyl_pairs(L)
    assert D.diagnostics['commutant_method'] == 'twirl'
    assert D.diagnostics['residuals']['algebra'] < TOL and D.diagnostics['residuals']['commutant'] < TOL


# ----------------------------------------------------------------------------
# Z_n, regular representation
# ----------------------------------------------------------------------------

@pytest.mark.parametrize('n', [2, 5, 8])
def test_cyclic_group_regular_representation(n):
    shift = torch.roll(torch.eye(n, dtype=C128), 1, 0)
    G = torch.stack([torch.linalg.matrix_power(shift, j) for j in range(n)])
    D = dc.from_group(G)
    assert D.k_factors == (1,) * n and D.multiplicities == (1,) * n          # A = A' = C^n
    F = torch.fft.fft(torch.eye(n, dtype=C128)) / math.sqrt(n)               # DFT
    M = (F.mH @ D.V).abs()
    # V = F x (permutation) x (diagonal phases): |F^* V| is a permutation matrix
    assert (M - M.round()).abs().max() < TOL
    assert torch.equal(M.round().sum(0), torch.ones(n, dtype=M.dtype)) and torch.equal(M.round().sum(1), torch.ones(n, dtype=M.dtype))


# ----------------------------------------------------------------------------
# random hidden block structure
# ----------------------------------------------------------------------------

HIDDEN = [
    ((1, 1, 2, 2, 3, 1), (1, 3, 1, 2, 2, 2)),        # k = 1 and m = 1, repeated k with different m
    ((2, 2, 3), (1, 1, 1)),                          # repeated (k, m): inequivalent irreps of equal size
    ((1, 1, 1, 1), (2, 2, 3, 1)),                    # commutative algebra
    ((4,), (3,)),                                    # a factor M_4 (x) 1_3
    ((3, 1, 2), (1, 4, 3)),
]


@pytest.mark.parametrize('ks,ms', HIDDEN)
@pytest.mark.parametrize('seed', [0, 1])
def test_hidden_block_structure(ks, ms, seed):
    S, W, ranges = hidden_algebra(ks, ms, seed)
    D = dc.from_generators(S, seed=seed)
    assert km(D) == sorted(zip(ks, ms))
    assert D.dims == (sum(k * k for k in ks), sum(m * m for m in ms))
    eye = torch.eye(D.n, dtype=C128)
    assert (D.V.mH @ D.V - eye).abs().max() < TOL
    for a in list(S) + [random_element(ks, ms, W, seed + 10)]:
        assert block_structure_residual(D, a, 'algebra') < TOL
    for _ in range(2):
        x = random_element(ks, ms, W, seed + 20 + _, side='commutant')
        assert block_structure_residual(D, x, 'commutant') < TOL
    # the central projections are the true ones
    true_P = [(W[:, s:e] @ W[:, s:e].mH, (k, m)) for (s, e), k, m in zip(ranges, ks, ms)]
    for P, k, m in zip(D.central_projections(), D.k_factors, D.multiplicities):
        assert min(frob(P - Q) for Q, km_ in true_P if km_ == (k, m)) < TOL
    assert D.verify()['ok']


@pytest.mark.parametrize('ks,ms', HIDDEN[:3])
def test_double_commutant(ks, ms):
    S, W, _ = hidden_algebra(ks, ms, 5)
    C1 = dc.commutant_basis(S)
    C2 = dc.commutant_basis(C1)
    assert C1.shape[0] == sum(m * m for m in ms)
    assert C2.shape[0] == sum(k * k for k in ks)                 # A'' = A: dimension ...
    assert in_span(S, C2) < TOL                                    # ... and inclusion
    D = dc.from_generators(S)
    assert in_span(D.basis('commutant'), C1) < TOL and in_span(D.basis('algebra'), C2) < TOL
    for which, dim in (('algebra', D.dim), ('commutant', D.commutant_dim)):
        B = D.basis(which)
        gram = torch.einsum('aij,bij->ab', B.conj(), B)
        assert B.shape[0] == dim and (gram - torch.eye(dim, dtype=C128)).abs().max() < TOL


def test_maps_round_trip_and_products():
    ks, ms = (1, 2, 2, 3), (2, 1, 3, 2)
    S, W, _ = hidden_algebra(ks, ms, 7)
    D = dc.from_generators(S)
    a = torch.stack([random_element(ks, ms, W, s) for s in range(4)])          # batch of 4 elements of A
    b = torch.stack([random_element(ks, ms, W, s + 50) for s in range(4)])
    A, B = D.to_algebra(a), D.to_algebra(b)
    assert isinstance(A, TypeIAlgebra.Operator) and A.algebra is D.algebra
    assert A.matrix.shape == (4, D.C, max(D.k_factors), max(D.k_factors))
    assert (D.from_algebra(A) - a).abs().max() < TOL                            # round trip
    assert (D.to_algebra(a @ b).matrix - (A @ B).matrix).abs().max() < TOL     # homomorphism
    assert (D.to_algebra(a.mH).matrix - A.adjoint().matrix).abs().max() < TOL  # *-homomorphism
    # Tr_{C^n}(a) = sum_c m_c Tr(a_c) (course, after Theorem 2.6)
    tr_blocks = torch.einsum('bcii->bc', A.matrix) @ torch.tensor(D.multiplicities, dtype=C128)
    assert (tr_blocks - torch.einsum('bii->b', a)).abs().max() < TOL
    # single matrix -> batch 1; Operator -> matrices
    assert D.to_algebra(a[0]).matrix.shape[0] == 1 and D.from_algebra(A).shape == (4, D.n, D.n)
    # random operator of the abstract algebra -> an element of A
    U = D.algebra.random_unitary_operator(batch_size=3)
    u = D.from_algebra(U)
    assert max(block_structure_residual(D, x) for x in u) < TOL
    assert (u @ u.mH - torch.eye(D.n, dtype=C128)).abs().max() < TOL
    # the commutant side
    x = torch.stack([random_element(ks, ms, W, s + 80, 'commutant') for s in range(3)])
    y = torch.stack([random_element(ks, ms, W, s + 90, 'commutant') for s in range(3)])
    X, Y = D.to_commutant(x), D.to_commutant(y)
    assert X.algebra.k_factors == list(D.multiplicities)
    assert (D.from_commutant(X) - x).abs().max() < TOL
    assert (D.to_commutant(x @ y).matrix - (X @ Y).matrix).abs().max() < TOL
    a3, x3 = D.from_algebra(A)[:3], D.from_commutant(X)
    assert (a3 @ x3 - x3 @ a3).abs().max() < TOL                                # A and A' commute
    # conditional expectations: trace preserving, idempotent, A-bimodule maps, twirl-like on the commutant
    z = torch.randn(2, D.n, D.n, dtype=C128, generator=torch.Generator().manual_seed(4))
    E = D.conditional_expectation(z)
    assert (torch.einsum('bii->b', E) - torch.einsum('bii->b', z)).abs().max() < TOL
    assert (D.conditional_expectation(E) - E).abs().max() < TOL
    assert (D.conditional_expectation(a[:2] @ z @ b[:2]) - a[:2] @ E @ b[:2]).abs().max() < TOL
    Ec = D.conditional_expectation(z, 'commutant')
    assert (D.conditional_expectation(x[:2] @ z @ y[:2], 'commutant') - x[:2] @ Ec @ y[:2]).abs().max() < TOL
    # central projections: orthogonal, sum to 1, ranks k_c m_c, central in A and A'
    P = D.central_projections()
    assert (P.sum(0) - torch.eye(D.n, dtype=C128)).abs().max() < TOL
    assert (P @ P - P).abs().max() < TOL and (P[0] @ P[1]).abs().max() < TOL
    assert [round(float(torch.einsum('ii->', p).real)) for p in P] == list(D.sector_dims)
    assert (P @ a[0] - a[0] @ P).abs().max() < TOL and (P @ x[0] - x[0] @ P).abs().max() < TOL


def test_maps_are_sync_free():
    """to_algebra / from_algebra / to_commutant / from_commutant read no values back to the host."""
    ks, ms = (1, 2), (2, 2)
    S, W, _ = hidden_algebra(ks, ms, 2)
    D = dc.from_generators(S)
    a = torch.stack([random_element(ks, ms, W, s) for s in range(3)])
    D.to_algebra(a), D.to_commutant(a)                    # build the index maps once
    names = ('item', 'tolist', '__bool__', '__float__', '__int__', '__index__', 'cpu', 'numpy')
    calls = []
    saved = {nm: getattr(torch.Tensor, nm) for nm in names}

    def wrap(nm):
        f = saved[nm]

        def g(self, *args, **kw):
            calls.append(nm)
            return f(self, *args, **kw)
        return g
    try:
        for nm in names:
            setattr(torch.Tensor, nm, wrap(nm))
        A = D.to_algebra(a)
        D.from_algebra(A)
        X = D.to_commutant(a)
        D.from_commutant(X)
        D.conditional_expectation(a)
    finally:
        for nm, f in saved.items():
            setattr(torch.Tensor, nm, f)
    assert calls == []


# ----------------------------------------------------------------------------
# robustness
# ----------------------------------------------------------------------------

def test_input_not_star_closed():
    # a single nilpotent matrix: the *-algebra generated by E_12 and E_21 is M_2 (not the upper triangular algebra)
    E12 = torch.tensor([[0., 1.], [0., 0.]], dtype=torch.float64)
    D = dc.from_generators(E12)
    assert D.k_factors == (2,) and D.multiplicities == (1,)
    # hidden structure given by non-normal (upper triangular) generators only
    ks, ms = (2, 3, 1), (2, 1, 3)
    S, W, _ = hidden_algebra(ks, ms, 3, n_gens=2, kind='nilpotent')
    D = dc.from_generators(S)
    assert km(D) == sorted(zip(ks, ms))
    assert max(block_structure_residual(D, s) for s in S) < TOL
    assert block_structure_residual(D, random_element(ks, ms, W, 1)) < TOL


def test_real_input_is_decomposed_over_C():
    R = torch.tensor([[0., -1.], [1., 0.]], dtype=torch.float64)      # generates C over R, C (+) C over C
    D = dc.from_generators(R)
    assert D.k_factors == (1, 1) and D.multiplicities == (1, 1) and D.V.dtype == C128
    blocks = D.to_algebra(R).matrix[0, :, 0, 0]
    assert sorted(blocks.imag.tolist()) == pytest.approx([-1.0, 1.0], abs=TOL)


def test_reproducible_and_seed_independent():
    ks, ms = (1, 2, 2), (3, 1, 2)
    S, _, _ = hidden_algebra(ks, ms, 11)
    D1, D2, D3 = dc.from_generators(S, seed=4), dc.from_generators(S, seed=4), dc.from_generators(S, seed=9)
    assert torch.equal(D1.V, D2.V)
    assert km(D1) == km(D3)
    for P, Q in zip(D1.central_projections(), D3.central_projections()):    # same sector order (sorted by k, m)
        assert frob(P - Q) < TOL


@contextlib.contextmanager
def patched_sampler(cls, name, replacement_calls):
    """Replace the first len(replacement_calls) calls of cls.name by the given functions of the sampler."""
    original = getattr(cls, name)
    count = [0]

    def method(self):
        i = count[0]
        count[0] += 1
        if i < len(replacement_calls):
            return replacement_calls[i](self)
        return original(self)
    setattr(cls, name, method)
    try:
        yield count
    finally:
        setattr(cls, name, original)


def test_degenerate_random_element_is_refined():
    """A first random element of the commutant that is scalar (maximally degenerate) is split by later ones."""
    ks, ms = (1, 2, 2), (2, 2, 1)
    S, W, _ = hidden_algebra(ks, ms, 12)
    n = S.shape[-1]
    ident = lambda self: torch.eye(n, dtype=C128) / math.sqrt(n)
    with patched_sampler(dc._BasisSampler, 'hermitian', [ident]):
        D = dc.from_generators(S, seed=0)
    assert km(D) == sorted(zip(ks, ms))
    assert D.diagnostics['initial_clusters'] == 1 and D.diagnostics['refinement_rounds'] >= 2
    assert D.diagnostics['residuals']['algebra'] < TOL and D.diagnostics['residuals']['commutant'] < TOL


def test_near_degenerate_random_element():
    """Two copies of inequivalent irreps whose eigenvalues of K differ by 1e-9: merged, then split."""
    ks, ms = (1, 2), (2, 2)
    S, W, _ = hidden_algebra(ks, ms, 13)
    kappa = [torch.diag(torch.tensor([0.3, 0.5], dtype=C128)),
             torch.diag(torch.tensor([0.5 + 1e-9, 0.9], dtype=C128))]
    K = W @ torch.block_diag(*[torch.kron(torch.eye(k, dtype=C128), x) for k, x in zip(ks, kappa)]) @ W.mH
    K = K / torch.linalg.matrix_norm(K)
    with patched_sampler(dc._BasisSampler, 'hermitian', [lambda self: K]):
        D = dc.from_generators(S, seed=0)
    assert km(D) == sorted(zip(ks, ms))
    assert D.diagnostics['initial_clusters'] == 3 and D.diagnostics['refinement_rounds'] >= 2
    for a in S:
        assert block_structure_residual(D, a) < TOL
    assert block_structure_residual(D, random_element(ks, ms, W, 3, 'commutant'), 'commutant') < TOL


def test_failed_draw_is_retried():
    """A degenerate intertwiner sample (zero) makes all copies look inequivalent; verification rejects the
    result (sum m_c^2 != dim A') and a new draw succeeds."""
    ks, ms = (1, 2), (2, 2)
    S, _, _ = hidden_algebra(ks, ms, 14)
    zero = lambda self: torch.zeros(S.shape[-1], S.shape[-1], dtype=C128)
    with patched_sampler(dc._BasisSampler, 'general', [zero]):
        D = dc.from_generators(S, seed=0)
    assert D.diagnostics['attempts'] == 2 and len(D.diagnostics['failures']) == 1
    assert km(D) == sorted(zip(ks, ms))


def test_impossible_raises_clear_error():
    ks, ms = (1, 2), (2, 2)
    S, _, _ = hidden_algebra(ks, ms, 15)
    n = S.shape[-1]
    ident = lambda self: torch.eye(n, dtype=C128)
    with patched_sampler(dc._BasisSampler, 'hermitian', [ident] * 100), \
            patched_sampler(dc._BasisSampler, 'general', [ident] * 100):
        with pytest.raises(dc.DecompositionError, match='no verified decomposition'):
            dc.from_generators(S, max_tries=2)


def test_ill_determined_commutant_warns():
    """A symmetry broken at the 1e-9 level makes the commutant dimension tolerance dependent: warn."""
    G = dc.permutation_representation(3, 2, 'adjacent')
    g = torch.Generator().manual_seed(0)
    H = torch.randn(8, 8, dtype=C128, generator=g)
    with pytest.warns(dc.DecompositionWarning):
        dc.commutant_basis(torch.cat([G, (1e-9 * (H + H.mH))[None] + torch.eye(8, dtype=C128)]))
    with warnings.catch_warnings():
        warnings.simplefilter('error', dc.DecompositionWarning)
        assert dc.commutant_basis(G).shape[0] == 20                   # exact symmetry: no warning


def test_group_input_validation():
    G = dc.permutation_representation(3, 2)
    with pytest.raises(ValueError, match='unitary'):
        dc.from_group(2 * G)
    with pytest.raises(ValueError):                                    # generators only: not a group
        dc.from_group(dc.permutation_representation(3, 2, 'adjacent'))
    with pytest.raises(ValueError):
        dc.commutant_basis(G, method='twirl')                          # twirl needs group=True
    with pytest.raises(ValueError):
        dc.wedderburn(G, method='deterministic')


def test_single_precision():
    ks, ms = (1, 2, 3), (3, 2, 1)
    S, W, _ = hidden_algebra(ks, ms, 16, dtype=torch.complex64)
    D = dc.from_generators(S)
    assert D.precision == 'single' and D.V.dtype == torch.complex64 and D.algebra.precision == 'single'
    assert km(D) == sorted(zip(ks, ms))
    assert D.diagnostics['residuals']['algebra'] < 1e-4 and D.diagnostics['residuals']['commutant'] < 1e-4
    A = D.to_algebra(S)
    assert (D.from_algebra(A) - S).abs().max() < 1e-4
    # in double precision the rounding errors of the input (~1e-7) break the symmetry: with the default
    # tol = sqrt(eps) ~ 1.5e-8 the generators generate M_10 (and the ill-determined rank is reported) ...
    with pytest.warns(dc.DecompositionWarning):
        Dd = dc.from_generators(S, precision='double')
    assert Dd.k_factors == (10,) and Dd.multiplicities == (1,) and Dd.V.dtype == C128
    # ... while a tolerance above the input accuracy recovers the structure
    Dd = dc.from_generators(S, precision='double', tol=1e-4)
    assert km(Dd) == sorted(zip(ks, ms)) and Dd.diagnostics['residuals']['algebra'] < 1e-5


def test_given_commutant_and_verify_arguments():
    ks, ms = (2, 1), (1, 3)
    S, W, _ = hidden_algebra(ks, ms, 17)
    C = dc.commutant_basis(S)
    D = dc.from_generators(S, commutant=C)
    assert D.diagnostics['commutant_method'] == 'given' and km(D) == sorted(zip(ks, ms))
    res = D.verify(generators=S, commutant=C)
    assert res['ok'] and res['algebra'] < TOL and res['commutant'] < TOL
    # a matrix outside A is detected
    bad = D.verify(generators=torch.randn(D.n, D.n, dtype=C128, generator=torch.Generator().manual_seed(0)))
    assert not bad['ok'] and bad['algebra'] > 0.1


# ----------------------------------------------------------------------------
# real structure: types R, C, H of the summands of a real algebra (course, Section 4.7)
# ----------------------------------------------------------------------------

def _quaternion_left(a):
    """Real 4 x 4 matrix of left multiplication by the quaternion unit a in (1, i, j, k) on H = R^4."""
    def mul(p, q):
        a1, b1, c1, d1 = p
        a2, b2, c2, d2 = q
        return [a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2, a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
                a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2, a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2]
    E = torch.eye(4, dtype=torch.float64)
    return torch.tensor([mul(E[a].tolist(), E[i].tolist()) for i in range(4)], dtype=torch.float64).T


def _so3_generators():
    L = torch.zeros(3, 3, 3, dtype=torch.float64)
    for a, (b, c) in enumerate([(1, 2), (2, 0), (0, 1)]):
        L[a, b, c], L[a, c, b] = -1.0, 1.0
    return L


def _realify(X):
    """Complex p x p -> real 2p x 2p [[Re, -Im], [Im, Re]] (the action of X on C^p = R^{2p})."""
    return torch.cat([torch.cat([X.real, -X.imag], -1), torch.cat([X.imag, X.real], -1)], -2)


def _quaternionic(A, B):
    """The quaternion matrix A + B j (A, B complex p x p) as the complex 2p x 2p matrix [[A, -conj B], [B, conj A]]."""
    return torch.cat([torch.cat([A, -B.conj()], -1), torch.cat([B, A.conj()], -1)], -2)


def hidden_real_algebra(spec, seed, n_gens=3):
    """
    Random real generators of O ((+) blocks) O^T with a random real orthogonal O. spec: list of (type, p, copies):
    'R': M_p(R) acting on R^p; 'C': M_p(C) realified on R^{2p}; 'H': M_p(H) realified on R^{4p}.
    Over C these give: R -> (k=p, m=copies); C -> two conjugate sectors (k=p, m=copies); H -> (k=2p, m=2 copies).
    """
    g = torch.Generator().manual_seed(seed)
    gens = []
    for _ in range(n_gens):
        blocks = []
        for t, p, copies in spec:
            if t == 'R':
                X = torch.randn(p, p, dtype=torch.float64, generator=g)
            elif t == 'C':
                X = _realify(torch.randn(p, p, dtype=C128, generator=g))
            else:
                Q = _quaternionic(torch.randn(p, p, dtype=C128, generator=g), torch.randn(p, p, dtype=C128, generator=g))
                X = _realify(Q)
            blocks += [X] * copies
        gens.append(torch.block_diag(*blocks))
    n = gens[0].shape[0]
    O = torch.linalg.qr(torch.randn(n, n, dtype=torch.float64, generator=g))[0]
    return torch.stack([O @ x @ O.T for x in gens])


def _expected_real_types(spec):
    out = []
    for t, p, copies in spec:
        if t == 'R':
            out.append(('R', p, copies))
        elif t == 'C':
            out += [('C', p, copies)] * 2
        else:
            out.append(('H', 2 * p, 2 * copies))
    return sorted(out, key=lambda x: (x[1], x[2], x[0]))


def _found_real_types(D, rs):
    return sorted(zip(rs.types, D.k_factors, D.multiplicities), key=lambda x: (x[1], x[2], x[0]))


def test_real_structure_examples():
    R2 = torch.tensor([[0., -1.], [1., 0.]], dtype=torch.float64)
    D = dc.from_generators(R2)                               # C over R: one complex pair
    rs = D.real_structure()
    assert rs.ok and rs.types == ('C', 'C') and rs.partner == (1, 0)
    D = dc.wedderburn(_so3_generators())                     # spin 1 of SO(3): real type
    rs = D.real_structure()
    assert rs.ok and rs.types == ('R',) and D.k_factors == (3,)
    D = dc.wedderburn(torch.stack([_quaternion_left(a) for a in (1, 2, 3)]))   # SU(2) on R^4 = H: quaternionic
    rs = D.real_structure()
    assert rs.ok and rs.types == ('H',) and (D.k_factors, D.multiplicities) == ((2,), (2,))
    assert torch.allclose(rs.J[0] @ rs.J[0].conj(), -torch.eye(2, dtype=C128), atol=TOL)
    P = torch.roll(torch.eye(3, dtype=torch.float64), 1, 0)  # Z_3 regular: trivial (R) + a conjugate pair
    rs = dc.from_group(torch.stack([torch.eye(3, dtype=torch.float64), P, P @ P])).real_structure()
    assert sorted(rs.types) == ['C', 'C', 'R']
    rs = dc.from_group(dc.permutation_representation(3, 2)).real_structure()   # S_3: all irreps are real
    assert set(rs.types) == {'R'}


def test_real_structure_quaternion_group():
    """Regular representation of Q_8: four real characters and the 2-dimensional quaternionic one, twice."""
    units = []
    for a in range(4):
        for s in (1.0, -1.0):
            units.append(s * _quaternion_left(a))
    G = torch.stack(units)                                   # Q_8 acting on R^4 by left multiplication
    idx = {tuple(g.flatten().tolist()): i for i, g in enumerate(units)}
    reg = torch.zeros(8, 8, 8, dtype=torch.float64)
    for i, g in enumerate(units):
        for j, h in enumerate(units):
            reg[i, idx[tuple((g @ h).flatten().tolist())], j] = 1.0
    D = dc.from_group(reg)
    rs = D.real_structure()
    assert rs.ok and _found_real_types(D, rs) == [('R', 1, 1)] * 4 + [('H', 2, 2)]
    D4 = dc.from_group(G)
    assert D4.real_structure().types == ('H',)


@pytest.mark.parametrize("spec,seed", [
    ([('R', 3, 2), ('C', 2, 1), ('H', 1, 3)], 0),
    ([('R', 1, 1), ('R', 2, 3), ('C', 3, 2)], 1),
    ([('H', 2, 1), ('C', 1, 3), ('R', 4, 1)], 2),
])
def test_hidden_real_algebra_types(spec, seed):
    S = hidden_real_algebra(spec, seed)
    D = dc.wedderburn(S)
    rs = D.real_structure()
    assert rs.ok, rs.residuals
    assert _found_real_types(D, rs) == _expected_real_types(spec)
    for c, d in enumerate(rs.partner):
        assert rs.partner[d] == c and (rs.types[c] == 'C') == (d != c)
        assert torch.allclose(rs.J[c] @ rs.J[c].mH, torch.eye(D.k_factors[c], dtype=C128), atol=1e-8)


@pytest.mark.parametrize("spec,seed", [([('R', 3, 2), ('C', 2, 1), ('H', 1, 3)], 0), ([('C', 3, 2), ('R', 2, 1)], 3)])
def test_real_basis(spec, seed):
    """O is real orthogonal, O^T a O is block diagonal for a in A, real blocks are alpha (x) 1_m with alpha real,
    and a real symmetric element of A' has a Kramers-degenerate spectrum on each quaternionic sector."""
    S = hidden_real_algebra(spec, seed)
    D = dc.wedderburn(S)
    O, blocks = D.real_basis()
    n = D.n
    assert O.dtype == torch.float64 and torch.allclose(O.T @ O, torch.eye(n, dtype=torch.float64), atol=1e-12)
    mask = torch.zeros(n, n, dtype=torch.bool)
    for b in blocks:
        mask[b['columns'], b['columns']] = True
    rs = D.real_structure()
    for s in S:
        G = O.T @ s @ O
        assert G[~mask].abs().max() < 1e-10
        for b in blocks:
            if b['type'] == 'R':
                c = b['sectors'][0]
                k, m = D.k_factors[c], D.multiplicities[c]
                blk = G[b['columns'], b['columns']].reshape(k, m, k, m)
                alpha = blk[:, 0, :, 0]
                assert torch.allclose(blk, torch.einsum('ik,jl->ijkl', alpha, torch.eye(m, dtype=alpha.dtype)), atol=1e-10)
    # a real symmetric element of the commutant
    g = torch.Generator().manual_seed(7)
    X = D.from_commutant(torch.randn(1, D.C, max(D.multiplicities), max(D.multiplicities), dtype=C128, generator=g))[0]
    X = X.real
    X = (X + X.T) / 2
    x = D.to_commutant(X).matrix[0]
    for c, t in enumerate(rs.types):
        m = D.multiplicities[c]
        ev = torch.linalg.eigvalsh(x[c, :m, :m])
        if t == 'H':                                          # Kramers: eigenvalues come in equal pairs
            assert (ev[0::2] - ev[1::2]).abs().max() < 1e-10
        if t == 'C':                                          # conjugate sectors carry the same spectrum
            d = rs.partner[c]
            assert torch.allclose(ev, torch.linalg.eigvalsh(x[d, :m, :m]), atol=1e-10)


def test_real_structure_single_precision():
    S = hidden_real_algebra([('R', 2, 2), ('C', 2, 1), ('H', 1, 2)], 5).float()
    D = dc.wedderburn(S)
    rs = D.real_structure()
    assert rs.ok and _found_real_types(D, rs) == _expected_real_types([('R', 2, 2), ('C', 2, 1), ('H', 1, 2)])


def test_real_structure_requires_conjugation_closed_algebra():
    g = torch.Generator().manual_seed(0)
    U = haar(2, g)
    gens = torch.stack([torch.block_diag(a, U @ a @ U.mH) for a in torch.randn(2, 2, 2, dtype=C128, generator=g)])
    D = dc.wedderburn(gens)
    with pytest.raises(ValueError, match="closed under complex conjugation"):
        D.real_structure()
    # complex input that is closed under conjugation is accepted
    rs = dc.wedderburn(torch.tensor([[1, 0], [0, 1j]], dtype=C128)).real_structure()
    assert rs.types == ('R', 'R')


@pytest.mark.skipif(not torch.cuda.is_available(), reason="needs CUDA")
def test_cuda():
    S = hidden_real_algebra([('R', 3, 2), ('C', 2, 1), ('H', 1, 3)], 0)
    D = dc.wedderburn(S, device='cuda')
    assert D.V.device.type == 'cuda' and D.verify()['ok']
    Dc = dc.wedderburn(S)
    assert (D.k_factors, D.multiplicities) == (Dc.k_factors, Dc.multiplicities)
    rs = D.real_structure()
    assert rs.ok and _found_real_types(D, rs) == _expected_real_types([('R', 3, 2), ('C', 2, 1), ('H', 1, 3)])
    O, _ = D.real_basis()
    assert O.device.type == 'cuda' and torch.allclose(O.T @ O, torch.eye(D.n, dtype=O.dtype, device='cuda'), atol=1e-12)
    G = dc.permutation_representation(4, 2).to('cuda')
    Dg = dc.from_group(G)
    Dh = dc.from_group(G.cpu())
    assert (Dg.k_factors, Dg.multiplicities) == (Dh.k_factors, Dh.multiplicities)
    A = Dg.to_algebra(G)
    assert torch.allclose(Dg.from_algebra(A), G.to(Dg.dtype), atol=1e-10)
