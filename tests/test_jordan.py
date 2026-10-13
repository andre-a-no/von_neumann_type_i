"""
Jordan structure (torch_vn_algebra.jordan), checked against independent references: dense NumPy/SciPy block-diagonal
matrices for the operator-level maps, the identities of Jordan algebras (Jordan identity, P(y) = 2L(y)^2 - L(y^2),
fundamental formula, trace identity, P(y)^{-1} = P(y^{-1}), rank P(y) = (rank y)^2), exact dimension formulas for the
generated Jordan algebras (Herm_n(R) inside Herm_n(C), spin factors, self-adjoint parts of *-algebras, one generator),
the *-algebra dimensions of decompose.wedderburn, and the matrix representation of spin factors (spectra, matrix
functions by scipy/torch, positive semidefiniteness).
"""
import math

import numpy as np
import pytest
import scipy.linalg as sl
import torch

from torch_vn_algebra import TypeIAlgebra, decompose as dc
from torch_vn_algebra import jordan as jd
from test_sync_free import count_host_reads

TOL = {'double': 1e-12, 'single': 2e-5}
KS = (3, 1, 4, 0, 2)                                     # unequal blocks and an empty sector


# ----------------------------------------------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------------------------------------------
def _alg(complex_=True, precision='double', ks=KS):
    return TypeIAlgebra([max(k, 1) for k in ks], list(ks), complex_valued=complex_, precision=precision)


def _rand(alg, B, kind='sa', gen=None):
    M = torch.randn(B, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype, generator=gen)
    if kind == 'general':
        return alg.operator(M)
    if kind == 'sa':
        return alg.operator(M + M.mH, is_self_adjoint=True)
    M = alg.operator(M).matrix                            # padding zero before the product
    return alg.operator(M @ M.mH, is_self_adjoint=True, is_positive=True)


def _blocks(alg, m, b):
    """Active blocks of sample b as NumPy arrays (empty sectors skipped)."""
    return [m[b, c, :k, :k].detach().cpu().numpy() for c, k in enumerate(alg.k_factors) if k]


def _dense(alg, m, b):
    return sl.block_diag(*_blocks(alg, m, b))


def _rel(x, y):
    x = x.matrix if hasattr(x, 'matrix') else x
    y = y.matrix if hasattr(y, 'matrix') else y
    x, y = torch.as_tensor(x), torch.as_tensor(y)
    return float((x - y).abs().max() / max(1.0, float(y.abs().max())))


def _jp(a, b):
    return 0.5 * (a @ b + b @ a)


FIELDS = [(True, 'double'), (False, 'double'), (True, 'single'), (False, 'single')]
FIELD_IDS = ['complex-double', 'real-double', 'complex-single', 'real-single']


# ----------------------------------------------------------------------------------------------------------------
# operator level
# ----------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('complex_,precision', FIELDS, ids=FIELD_IDS)
def test_jordan_product_and_triple_match_dense(complex_, precision):
    torch.manual_seed(0)
    alg, B = _alg(complex_, precision), 4
    A, Bo, Co = _rand(alg, B, 'general'), _rand(alg, B, 'general'), _rand(alg, B, 'general')
    P = jd.jordan_product(A, Bo)
    T = jd.jordan_triple(A, Bo, Co)
    Q = jd.quadratic_representation(A, Bo)
    for b in range(B):
        a, bb, c = (_dense(alg, X.matrix, b) for X in (A, Bo, Co))
        for got, ref in ((P, 0.5 * (a @ bb + bb @ a)), (T, 0.5 * (a @ bb @ c + c @ bb @ a)),
                         (Q, a @ bb @ a.conj().T)):
            g = _dense(alg, got.matrix, b)
            assert np.abs(g - ref).max() <= TOL[precision] * max(1, np.abs(ref).max())
    for X in (P, T, Q):                                   # the padding stays zero
        assert (X.matrix * ~alg._const('block_mask', torch.bool)).abs().max() == 0
        assert X.matrix.dtype == alg.dtype


def test_flags():
    alg = _alg()
    x, y, p, g = _rand(alg, 2), _rand(alg, 2), _rand(alg, 2, 'pos'), _rand(alg, 2, 'general')
    assert jd.jordan_product(x, y)._is_self_adjoint is True
    assert jd.jordan_product(x, y)._is_positive is None
    assert jd.jordan_product(x, x)._is_positive is True                       # x∘x = x^2
    assert jd.jordan_product(x, g)._is_self_adjoint is None
    assert jd.quadratic_representation(g, p)._is_positive is True             # congruence of a positive operator
    assert jd.quadratic_representation(g, x)._is_self_adjoint is True
    assert jd.quadratic_representation(x, g)._is_self_adjoint is None
    assert jd.jordan_triple(x, p, x)._is_positive is True
    assert jd.jordan_triple(x, p, y)._is_positive is None
    assert jd.jordan_power(x, 4)._is_positive is True and jd.jordan_power(x, 3)._is_positive is None
    assert jd.jordan_power(x, 0)._is_projection is True
    # the flags are true: checked by the library's own (synchronising) property tests on fresh operators
    for op in (jd.jordan_product(x, x), jd.quadratic_representation(g, p), jd.jordan_triple(x, p, x),
               jd.jordan_power(x, 4)):
        fresh = alg.operator(op.matrix)
        assert fresh.is_self_adjoint and fresh.is_positive


def test_product_of_positive_operators_is_not_positive():
    """p = diag(1, 0), q = [[1, 1], [1, 1]]/2 are projections; p∘q = [[1, 1/2], [1/2, 0]]/2 has det -1/16 < 0."""
    alg = TypeIAlgebra([2], [2], precision='double')
    p = alg.operator(torch.tensor([[[[1., 0.], [0., 0.]]]]), is_self_adjoint=True, is_positive=True)
    q = alg.operator(torch.tensor([[[[.5, .5], [.5, .5]]]]), is_self_adjoint=True, is_positive=True)
    pq = jd.jordan_product(p, q)
    assert pq._is_positive is None                       # no flag: not guaranteed
    assert pq.is_positive is False
    assert abs(float(pq.lambda_min) - (1 - math.sqrt(2)) / 4) < 1e-15
    # but P(p)q = pqp and P(q)p are positive
    assert jd.quadratic_representation(p, q).is_positive
    assert float(jd.quadratic_representation(q, p).lambda_min) > -1e-15


@pytest.mark.parametrize('complex_,precision', FIELDS, ids=FIELD_IDS)
def test_jordan_identities(complex_, precision):
    """Jordan identity, P(y) = 2L(y)^2 - L(y^2) = {y, ., y}, fundamental formula, triple-product formula, trace."""
    torch.manual_seed(1)
    alg, B = _alg(complex_, precision), 3
    x, y, z = _rand(alg, B), _rand(alg, B), _rand(alg, B)
    tol = TOL[precision]
    jp = jd.jordan_product
    x2 = jp(x, x)
    assert _rel(jp(jp(x2, y), x), jp(x2, jp(y, x))) < tol                                   # Jordan identity
    L = jd.jordan_multiplication(y)
    P = jd.quadratic_representation(y)
    assert _rel(L(L(x)) * 2.0 - jp(jp(y, y), x), P(x)) < tol                                # P = 2L^2 - L(y^2)
    assert _rel(jd.jordan_triple(y, x, y), P(x)) < tol                                      # P(y)x = {y, x, y}
    assert _rel(jd.jordan_triple(x, y, z), jp(jp(x, y), z) + jp(jp(y, z), x) - jp(jp(x, z), y)) < tol
    Px = jd.quadratic_representation(x)
    assert _rel(jd.quadratic_representation(P(x))(z), P(Px(P(z)))) < 10 * tol             # fundamental formula
    tr1, tr2 = P(x).trace, jp(jp(y, y), x).trace                                            # Tr P(y)x = Tr y^2 x
    assert float((tr1 - tr2).abs().max()) < tol * float(tr2.abs().max())
    # commutativity and the unit
    one = alg.identity(B)
    assert _rel(jp(x, y), jp(y, x)) == 0.0 and _rel(jp(one, x), x) < tol


@pytest.mark.parametrize('complex_', [True, False])
def test_quadratic_representation_invertibility(complex_):
    """P(y) is invertible on M_sa iff y is: P(y^{-1}) P(y) = id for invertible y, and on a block M_k a y of rank r
    gives rank P(y) = r^2 on Herm_k(C) (r(r+1)/2 on Herm_k(R)); P(y) kills the kernel projection of y."""
    torch.manual_seed(2)
    alg, B = _alg(complex_), 2
    y = alg.operator_from_eigenvalues(lambda k: torch.where(torch.rand(k) < .5, -1., 1.) * (0.5 + torch.rand(k)),
                                      batch_size=B)
    x = _rand(alg, B)
    yinv = y.inverse()
    assert _rel(jd.quadratic_representation(yinv, jd.quadratic_representation(y, x)), x) < 1e-12
    k = 4
    one = TypeIAlgebra([k], [k], complex_valued=complex_, precision='double')
    field = 'complex' if complex_ else 'real'
    coords = jd._HermCoords(k, field, torch.float64, 'cpu')
    basis = coords.mat(torch.eye(coords.N, dtype=torch.float64))                    # orthonormal basis of Herm_k
    E = one.operator(basis.unsqueeze(1), is_self_adjoint=True)                      # batch N
    for r in (4, 2, 1):
        w = torch.tensor([1.5, -0.7, 2.0, 0.3][:r] + [0.0] * (k - r), dtype=torch.float64)
        U = one.random_unitary(k).to(one.dtype)
        Y = one.operator(((U * w.to(U.dtype)) @ U.mH)[None, None], is_self_adjoint=True)
        img = coords.vec(jd.quadratic_representation(Y, E).matrix[:, 0])           # (N, N): the matrix of P(y)
        s = torch.linalg.svdvals(img)
        rank = int((s > 1e-10 * s[0]).sum())
        assert rank == (r * r if complex_ else r * (r + 1) // 2)
        if r < k:                                                                   # kernel projection
            Pk = U[:, r:] @ U[:, r:].mH
            assert jd.quadratic_representation(Y, one.operator(Pk[None, None])).matrix.abs().max() < 1e-12


@pytest.mark.parametrize('complex_,precision', FIELDS, ids=FIELD_IDS)
def test_jordan_power(complex_, precision):
    torch.manual_seed(3)
    alg, B = _alg(complex_, precision), 3
    x = _rand(alg, B) * 0.5
    tol = TOL[precision]
    for p in range(0, 7):
        ref = torch.linalg.matrix_power(x.matrix, p) * alg._const('block_mask', torch.bool) if p else \
            alg.identity(B).matrix
        assert _rel(jd.jordan_power(x, p), ref) < tol
    xp = lambda p: jd.jordan_power(x, p)
    assert _rel(jd.jordan_product(xp(2), xp(3)), xp(5)) < tol                     # power associativity
    rec = x
    for _ in range(3):
        rec = jd.jordan_product(x, rec)                                             # x∘(x∘(x∘x))
    assert _rel(rec, xp(4)) < tol
    with pytest.raises(ValueError):
        jd.jordan_power(x, -1)


def test_batch_broadcasting():
    torch.manual_seed(4)
    alg = _alg()
    x1, y = _rand(alg, 1), _rand(alg, 5)
    out = jd.jordan_product(x1, y)
    assert out.shape == (5, alg.C, alg.k_max, alg.k_max) and out._batch == 5
    assert _rel(out.matrix[3], _jp(x1.matrix[0], y.matrix[3])) < 1e-13
    assert jd.jordan_triple(x1, y, x1).matrix.shape[0] == 5 and jd.jordan_triple(x1, y, x1)._batch == 5


def test_gradients():
    torch.manual_seed(5)
    alg = _alg(ks=(2, 1, 3))
    mask = alg._const('block_mask', torch.bool)
    X0 = torch.randn(2, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype) * mask
    Y0 = torch.randn(2, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype) * mask

    def f(X, Y):
        x, y = alg.operator(X), alg.operator(Y)
        return (jd.jordan_product(x, y).matrix.abs() ** 2).sum() + \
            jd.quadratic_representation(y, x).trace.real.sum() + jd.jordan_triple(x, y, y).trace.imag.sum()
    assert torch.autograd.gradcheck(f, (X0.requires_grad_(), Y0.requires_grad_()))


def test_no_host_reads():
    torch.manual_seed(6)
    alg = _alg()
    x, y, p, g = _rand(alg, 4), _rand(alg, 4), _rand(alg, 4, 'pos'), _rand(alg, 4, 'general')

    def run():
        return (jd.jordan_product(x, y).matrix, jd.jordan_product(g, g).matrix,
                jd.quadratic_representation(y, p).matrix, jd.jordan_triple(x, y, g).matrix,
                jd.jordan_power(x, 5).matrix, jd.jordan_power(x, 0).matrix,
                jd.transpose_map(g).matrix, jd.jordan_multiplication(x)(y).trace)
    run()
    with count_host_reads() as calls:
        run()
    assert not calls, sorted(set(calls))


# ----------------------------------------------------------------------------------------------------------------
# transpose: a Jordan automorphism that is not a *-automorphism, positive but not completely positive
# ----------------------------------------------------------------------------------------------------------------
def test_transpose_is_jordan_automorphism_reversing_products():
    torch.manual_seed(7)
    alg = _alg()
    x, y, g, h = _rand(alg, 3), _rand(alg, 3), _rand(alg, 3, 'general'), _rand(alg, 3, 'general')
    T = jd.transpose_map
    assert _rel(T(jd.jordan_product(x, y)), jd.jordan_product(T(x), T(y))) < 1e-13
    assert _rel(T(jd.quadratic_representation(y, x)), jd.quadratic_representation(T(y), T(x))) < 1e-13
    assert _rel(T(g @ h), T(h) @ T(g)) < 1e-13                                 # anti-multiplicative
    assert _rel(T(g @ h), T(g) @ T(h)) > 1e-2                                  # not multiplicative
    assert _rel(T(g.adjoint()), T(g).adjoint()) == 0.0
    p = _rand(alg, 3, 'pos')
    tp = T(p)
    assert tp._is_positive is True
    for a, b in zip(p.eigenvalues(), alg.operator(tp.matrix).eigenvalues()):  # same spectrum
        assert a.numel() == 0 or (a - b).abs().max() < 1e-12
    # real algebra: the transpose is the adjoint, the identity on self-adjoint elements
    ralg = _alg(complex_=False)
    xr = _rand(ralg, 2)
    assert _rel(T(xr), xr) == 0.0


def test_transpose_is_not_completely_positive():
    """Choi matrix sum_ij E_ij (x) T(E_ij) of the transpose on every block M_k: the swap, eigenvalues +-1 with
    multiplicities k(k+1)/2 and k(k-1)/2."""
    alg = _alg()
    K = alg.k_max
    E = torch.zeros(K * K, alg.C, K, K, dtype=alg.dtype)
    for i in range(K):
        for j in range(K):
            E[i * K + j, :, i, j] = 1
    TE = jd.transpose_map(alg.operator(E)).matrix                             # padding masked by alg.operator
    for c, k in enumerate(alg.k_factors):
        if not k:
            continue
        choi = sum(np.kron(np.eye(k)[:, [i]] @ np.eye(k)[[j], :], TE[i * K + j, c, :k, :k].numpy())
                   for i in range(k) for j in range(k))
        w = np.linalg.eigvalsh(choi)
        assert np.allclose(np.sort(w), [-1.0] * (k * (k - 1) // 2) + [1.0] * (k * (k + 1) // 2), atol=1e-14)
        if k >= 2:
            assert w.min() < -0.5                                              # not completely positive


# ----------------------------------------------------------------------------------------------------------------
# Jordan subalgebras of Herm_n
# ----------------------------------------------------------------------------------------------------------------
def _herm(n, gen=None, real=False):
    dt = torch.float64 if real else torch.complex128
    a = torch.randn(n, n, dtype=dt, generator=gen)
    return a + a.mH


@pytest.mark.parametrize('field', ['real', 'complex'])
def test_coordinates_are_isometric(field):
    torch.manual_seed(8)
    n = 5
    c = jd._HermCoords(n, field, torch.float64, 'cpu')
    a, b = _herm(n, real=field == 'real'), _herm(n, real=field == 'real')
    assert c.N == (n * (n + 1) // 2 if field == 'real' else n * n)
    assert abs(float(c.vec(a) @ c.vec(b)) - float(torch.trace(a @ b).real)) < 1e-12
    assert (c.mat(c.vec(a)) - a).abs().max() < 1e-14


@pytest.mark.parametrize('n', [2, 3, 5, 6])
def test_real_symmetric_matrices(n):
    """(a) Two generic real symmetric matrices generate J = Herm_n(R), dim n(n+1)/2, while the complex *-algebra
    they generate is M_n(C): dim A_sa = n^2 (= dim_C A from decompose.wedderburn), so J is a proper subspace."""
    gen = torch.Generator().manual_seed(n)
    a, b = _herm(n, gen, real=True), _herm(n, gen, real=True)
    J = jd.jordan_algebra_generated([a, b])
    Jc = jd.jordan_algebra_generated([a, b], field='complex')
    A = jd.star_algebra_self_adjoint_part([a, b])
    assert J.field == 'real' and Jc.field == 'complex' and A.field == 'complex'
    assert J.dim == Jc.dim == n * (n + 1) // 2
    assert A.dim == n * n == dc.wedderburn([a, b]).dim
    assert not J.basis.is_complex()
    assert J.is_closed(exhaustive=True) and Jc.is_closed(exhaustive=True) and A.is_closed()
    assert not J.is_self_adjoint_part() and A.is_self_adjoint_part()
    assert bool(J.contains(_herm(n, gen, real=True))) and not bool(J.contains(1j * (a @ b - b @ a)))
    eye = torch.eye(J.dim, dtype=torch.float64)                                  # orthonormal basis
    assert (torch.einsum('aij,bji->ab', J.basis, J.basis) - eye).abs().max() < 1e-12
    r = jd.compare_with_star_algebra([a, b])
    assert (r['jordan_dim'], r['star_dim'], r['equal']) == (n * (n + 1) // 2, n * n, False)


@pytest.mark.parametrize('m', [1, 2, 3, 4, 5, 6, 7])
def test_spin_factors(m):
    """(b) span{1, s_1..s_m} of anticommuting symmetries is a Jordan algebra of dim m + 1 (s_i∘s_j = δ_ij 1); the
    *-algebra generated is M_D(C), D = 2^floor(m/2) (m >= 2), so J = A_sa exactly for m = 1 and m = 3
    (V_3 = Herm_2(C))."""
    S = jd.anticommuting_symmetries(m)
    D = 2 if m == 1 else 2 ** (m // 2)
    assert S.shape == (m, D, D)
    eye = torch.eye(D, dtype=S.dtype)
    for i in range(m):
        assert (S[i] - S[i].mH).abs().max() == 0
        for j in range(m):
            assert (S[i] @ S[j] + S[j] @ S[i] - 2 * (i == j) * eye).abs().max() == 0
    r = jd.compare_with_star_algebra(S)
    star_dim = 2 if m == 1 else D * D
    assert r['jordan_dim'] == m + 1
    assert r['star_dim'] == star_dim == dc.wedderburn(S).dim
    assert r['equal'] == (m in (1, 3))
    assert r['jordan'].is_self_adjoint_part() == (m in (1, 3))
    assert jd.is_jordan_subalgebra(torch.cat([eye[None], S])) and jd.is_jordan_subalgebra(torch.stack([eye, S[0]]))
    assert not jd.is_jordan_subalgebra(S[:1])                                   # s_1∘s_1 = 1 is missing
    if m >= 2:
        assert jd.is_jordan_subalgebra(torch.stack([eye, S[0], S[1], 1j * S[0] @ S[1]]))     # a copy of V_3
        # span{1, s_1, h}: closed in Herm_2 = V_3 (every subspace of a spin factor containing 1 is), not in Herm_4
        h = _herm(D, torch.Generator().manual_seed(m))
        assert jd.is_jordan_subalgebra(torch.stack([eye, S[0], h])) == (D == 2)


def test_hidden_block_algebra():
    """(c) Three generic Hermitian elements of A = U (⊕_c M_{k_c} ⊗ 1_{m_c}) U^* generate J = A_sa as a Jordan
    algebra.
    Two generate the same *-algebra but a smaller Jordan algebra if some k_c = 2: two Hermitian 2 x 2 matrices are
    simultaneously unitarily similar to real symmetric ones (their Pauli vectors span a plane), so they generate
    Herm_2(R) = V_2 (dim 3) in that block, not Herm_2(C) = V_3; for k_c >= 3 two generic elements suffice. And
    Herm_2(R) ⊕ Herm_3(C), hidden by a unitary, is a Jordan algebra that is not the self-adjoint part of the
    *-algebra M_2(C) ⊕ M_3(C) it generates."""
    gen = torch.Generator().manual_seed(9)
    ks, ms = (1, 2, 3), (2, 1, 2)
    n = sum(k * m for k, m in zip(ks, ms))
    U = torch.linalg.qr(torch.randn(n, n, dtype=torch.complex128, generator=gen))[0]

    def element(real_blocks=()):
        blocks = []
        for c, (k, mult) in enumerate(zip(ks, ms)):
            h = _herm(k, gen, real=c in real_blocks).to(torch.complex128)
            blocks.append(torch.kron(h, torch.eye(mult, dtype=torch.complex128)))
        return U @ torch.block_diag(*blocks) @ U.mH
    G = [element(), element(), element()]
    r = jd.compare_with_star_algebra(G)
    assert r['jordan_dim'] == r['star_dim'] == sum(k * k for k in ks) == dc.wedderburn(G).dim
    assert r['equal'] and r['jordan'].is_self_adjoint_part()
    assert bool(r['jordan'].contains(element()))                                  # a random element of A_sa
    r2 = jd.compare_with_star_algebra(G[:2])
    assert (r2['jordan_dim'], r2['star_dim']) == (1 + 3 + 9, 1 + 4 + 9) and dc.wedderburn(G[:2]).dim == 14
    ks3 = (1, 3, 4)                                                               # no block of size 2
    U3 = torch.linalg.qr(torch.randn(8, 8, dtype=torch.complex128, generator=gen))[0]
    G3 = [U3 @ torch.block_diag(*[_herm(k, gen) for k in ks3]) @ U3.mH for _ in range(2)]
    assert jd.compare_with_star_algebra(G3)['equal'] and jd.jordan_algebra_generated(G3).dim == 1 + 9 + 16

    V = torch.linalg.qr(torch.randn(5, 5, dtype=torch.complex128, generator=gen))[0]
    G2 = [V @ torch.block_diag(_herm(2, gen, real=True).to(torch.complex128), _herm(3, gen)) @ V.mH
          for _ in range(3)]
    r2 = jd.compare_with_star_algebra(G2)
    assert (r2['jordan_dim'], r2['star_dim'], r2['equal']) == (3 + 9, 4 + 9, False)


def test_single_generator():
    """(d) One Hermitian x: J = span of its spectral projections; dim = number of distinct eigenvalues (unital) or of
    distinct non-zero eigenvalues (non-unital: polynomials without constant term, which contain 1 if 0 is not an
    eigenvalue)."""
    gen = torch.Generator().manual_seed(10)
    U = torch.linalg.qr(torch.randn(8, 8, dtype=torch.complex128, generator=gen))[0]
    cases = [([0, 0, 1, 1, 1, 2, -1.5, 2], 4, 3),          # 0 is an eigenvalue: unital adds one dimension
             ([1, 1, 2, 3, 3, 3, -1, -1], 4, 4),           # 0 is not: the same algebra
             ([2] * 8, 1, 1), ([0] * 8, 1, 0),
             ([1, 1 + 1e-13, 2, 2, 2, 2, 2, 2], 2, 2),     # eigenvalues closer than the tolerance are merged
             ([1, 1 + 1e-3, 2, 2, 2, 2, 2, 2], 3, 3)]
    for w, d_unital, d_plain in cases:
        x = (U * torch.tensor(w, dtype=torch.float64).to(U.dtype)) @ U.mH
        assert jd.jordan_algebra_generated(x).dim == d_unital, w
        assert jd.jordan_algebra_generated(x, unital=False).dim == d_plain, w
        J = jd.jordan_algebra_generated(x)
        assert J.is_self_adjoint_part()                    # commutative: J = A_sa (A = polynomials in x)
        assert bool(J.contains(x @ x @ x)) and bool(J.contains(torch.eye(8, dtype=torch.complex128)))
    # eigenvalues 1 and 1 + 1e-7: a rank decision near the threshold. It is flagged, and the inaccurate direction does
    # not cascade to spurious dimensions (before the threshold was tied to the accuracy of the added directions, this
    # input gave dim 64 = all of Herm_8(C))
    U0 = torch.linalg.qr(torch.randn(8, 8, dtype=torch.complex128, generator=torch.Generator().manual_seed(0)))[0]
    x = (U0 * torch.tensor([1, 1 + 1e-7] + [2] * 6, dtype=torch.float64).to(U.dtype)) @ U0.mH
    with pytest.warns(jd.ClosureWarning):
        J = jd.jordan_algebra_generated(x)
    assert J.dim in (2, 3) and J.is_closed()
    # a real symmetric generator: the same dimension in real coordinates
    xr = torch.diag(torch.tensor([0., 0., 1., 3., 3.]))
    assert jd.jordan_algebra_generated(xr).dim == 3 and jd.jordan_algebra_generated(xr, unital=False).dim == 2


def test_generators_conventions():
    gen = torch.Generator().manual_seed(11)
    n = 4
    h1, h2 = _herm(n, gen), _herm(n, gen)
    J1 = jd.jordan_algebra_generated(h1 + 1j * h2)                     # non-Hermitian S = H1 + i H2
    J2 = jd.jordan_algebra_generated([h1, h2])
    assert J1.dim == J2.dim == n * n
    assert jd.jordan_algebra_generated(torch.zeros(n, n)).dim == 1
    assert jd.jordan_algebra_generated(torch.zeros(n, n), unital=False).dim == 0
    p = torch.diag(torch.tensor([1., 1., 0., 0.]))
    assert jd.jordan_algebra_generated(p, unital=False).dim == 1 and jd.jordan_algebra_generated(p).dim == 2
    with pytest.raises(ValueError):
        jd.jordan_algebra_generated(h1, field='real')
    with pytest.raises(ValueError):
        jd.is_jordan_subalgebra(torch.randn(2, n, n))
    a, b = _herm(n, gen, real=True), _herm(n, gen, real=True)
    span = torch.stack([torch.eye(n, dtype=torch.float64), a, b])
    assert not jd.is_jordan_subalgebra(span)
    J = jd.jordan_algebra_generated([a, b])
    assert jd.is_jordan_subalgebra(J.basis) and jd.is_jordan_subalgebra(J.basis, exhaustive=False)
    Jab = jd.JordanAlgebra(span / span.norm(dim=(-2, -1), keepdim=True), 'real', True, 'jordan', J.threshold)
    assert Jab.closure_defect(exhaustive=True) > 1e-3 and Jab.closure_defect(exhaustive=False) > 1e-3
    assert J.closure_defect(exhaustive=True) < 1e-12 and J.diagnostics['closure_defect'] < 1e-12


def test_single_precision_closure():
    gen = torch.Generator().manual_seed(12)
    a, b = _herm(4, gen, real=True).float(), _herm(4, gen, real=True).float()
    J = jd.jordan_algebra_generated([a, b])
    assert J.dim == 10 and J.basis.dtype == torch.float32 and J.diagnostics['precision'] == 'single'
    S = jd.anticommuting_symmetries(4, dtype=torch.complex64)
    r = jd.compare_with_star_algebra(S)
    assert (r['jordan_dim'], r['star_dim']) == (5, 16)


def test_seed_independence():
    gen = torch.Generator().manual_seed(13)
    a, b = _herm(5, gen), _herm(5, gen, real=True)
    dims = {jd.jordan_algebra_generated([a, b], seed=s).dim for s in (0, 1, 2)}
    assert dims == {25}


# ----------------------------------------------------------------------------------------------------------------
# spin factors
# ----------------------------------------------------------------------------------------------------------------
MS = [1, 2, 3, 4, 5]


def _sf_rand(V, B, gen, with_zero=True):
    x = torch.randn(B, V.dim, dtype=V.dtype, generator=gen)
    if with_zero:
        x[0, 1:] = 0                                      # u = 0: a multiple of the unit
    return x


@pytest.mark.parametrize('m', MS)
def test_spin_factor_representation(m):
    """phi(a, u) = a 1 + u·s is an injective Jordan homomorphism into Herm_D; it preserves spectra (λ± with
    multiplicity D/2), traces (Tr phi = D/2 tr) and determinants (det phi = det^(D/2))."""
    gen = torch.Generator().manual_seed(20 + m)
    V = jd.SpinFactor(m)
    x, y = _sf_rand(V, 6, gen), _sf_rand(V, 6, gen)
    X, Y = V.to_matrix(x), V.to_matrix(y)
    D = V.representation_dim
    assert (V.to_matrix(V.product(x, y)) - _jp(X, Y)).abs().max() < 1e-13
    assert (V.to_matrix(V.unit()) - torch.eye(D)).abs().max() == 0
    assert (V.from_matrix(X) - x).abs().max() < 1e-14
    assert (X - X.mH).abs().max() == 0
    lam = V.eigenvalues(x)
    ref = torch.linalg.eigvalsh(X)
    assert (lam.repeat_interleave(D // 2, dim=-1) - ref).abs().max() < 1e-13
    assert (torch.diagonal(X, dim1=-2, dim2=-1).sum(-1).real - D / 2 * V.trace(x)).abs().max() < 1e-13
    det_ref = V.determinant(x) ** (D // 2)
    assert (torch.linalg.det(X).real - det_ref).abs().max() < 1e-11 * (1 + det_ref.abs()).max()
    assert (V.determinant(x) - lam.prod(-1)).abs().max() < 1e-13
    assert (V.trace(x) - lam.sum(-1)).abs().max() < 1e-14
    # quadratic representation: P(x)y = x y x in the representation
    assert (V.to_matrix(V.quadratic_representation(x, y)) - X @ Y @ X).abs().max() < 1e-12
    # V_3 is all of Herm_2(C)
    if m == 3:
        H = _herm(2)
        assert (V.to_matrix(V.from_matrix(H)) - H).abs().max() < 1e-14


@pytest.mark.parametrize('m', MS)
def test_spin_factor_identities(m):
    gen = torch.Generator().manual_seed(30 + m)
    V = jd.SpinFactor(m)
    x, y, z = (_sf_rand(V, 5, gen) for _ in range(3))
    pr, P = V.product, V.quadratic_representation
    x2 = V.square(x)
    assert (pr(pr(x2, y), x) - pr(x2, pr(y, x))).abs().max() < 1e-12                   # Jordan identity
    assert (pr(x, y) - pr(y, x)).abs().max() == 0
    assert (P(P(y, x), z) - P(y, P(x, P(y, z)))).abs().max() < 1e-11 * (1 + P(P(y, x), z).abs().max())
    assert (V.inner(P(y, x), V.unit()) - V.inner(V.square(y), x)).abs().max() < 1e-12  # tr P(y)x = tr y^2 x
    assert (V.triple(y, x, y) - P(y, x)).abs().max() < 1e-12
    assert (V.inner(pr(x, y), z) - V.inner(y, pr(x, z))).abs().max() < 1e-12          # associative trace form
    assert (pr(x, V.inverse(x)) - V.unit()).abs().max() < 1e-10
    assert (P(V.inverse(y), P(y, x)) - x).abs().max() < 1e-9 * (1 + x.abs().max())    # P(y)^{-1} = P(y^{-1})
    assert (V.power(x, 5) - pr(V.power(x, 2), V.power(x, 3))).abs().max() < 1e-11 * (1 + V.power(x, 5).abs().max())
    assert (V.power(x, 0) - V.unit((5,))).abs().max() == 0


@pytest.mark.parametrize('m', MS)
def test_spin_factor_spectral_decomposition(m):
    gen = torch.Generator().manual_seed(40 + m)
    V = jd.SpinFactor(m)
    x = _sf_rand(V, 6, gen)
    lam, e = V.spectral_decomposition(x)
    em, ep = e[..., 0, :], e[..., 1, :]
    one = V.unit((6,))
    assert (V.product(em, em) - em).abs().max() < 1e-15 and (V.product(ep, ep) - ep).abs().max() < 1e-15
    assert V.product(em, ep).abs().max() < 1e-15 and (em + ep - one).abs().max() < 1e-15
    assert (lam[..., :1] * em + lam[..., 1:] * ep - x).abs().max() < 1e-14
    assert (V.trace(em) - 1).abs().max() < 1e-15                                       # primitive idempotents


@pytest.mark.parametrize('m', MS)
def test_spin_factor_functional_calculus(m):
    """f(x) against the matrix function of the representation: expm, the PSD square root, log; at u = 0 exactly
    f(a) 1."""
    gen = torch.Generator().manual_seed(50 + m)
    V = jd.SpinFactor(m)
    x = _sf_rand(V, 6, gen)
    X = V.to_matrix(x)
    E = torch.linalg.matrix_exp(X)
    assert (V.to_matrix(V.exp(x)) - E).abs().max() < 1e-12 * E.abs().max()
    assert (V.exp(x)[0, 1:]).abs().max() == 0 and abs(float(V.exp(x)[0, 0] - torch.exp(x[0, 0]))) < 1e-15
    c = V.random(6, gen, cone=True)                                                    # squares: in the cone
    C = V.to_matrix(c)
    w, Q = torch.linalg.eigh(C)
    sqrtC = (Q * w.clamp(min=0).sqrt().to(Q.dtype).unsqueeze(-2)) @ Q.mH
    assert (V.to_matrix(V.sqrt(c)) - sqrtC).abs().max() < 1e-12 * (1 + sqrtC.abs().max())
    assert (V.square(V.sqrt(c)) - c).abs().max() < 1e-12 * (1 + c.abs().max())
    assert V.in_cone(V.sqrt(c)).all()
    for b in range(3):                                                                  # scipy reference for log
        Lref = sl.logm(C[b].numpy() + 0.5 * np.eye(C.shape[-1]))
        assert np.abs(V.to_matrix(V.log(c[b] + 0.5 * V.unit())).numpy() - Lref).max() < 1e-10
    assert (V.log(V.exp(x)) - x).abs().max() < 1e-12 * (1 + x.abs().max())


@pytest.mark.parametrize('m', MS)
def test_lorentz_cone_is_the_psd_cone(m):
    gen = torch.Generator().manual_seed(60 + m)
    V = jd.SpinFactor(m)
    x = torch.randn(400, V.dim, dtype=V.dtype, generator=gen)
    x[:, 0] = x[:, 0].abs() * 2                                                        # about half in the cone
    edge = torch.cat([torch.ones(1), torch.zeros(m)]) + torch.cat([torch.zeros(1), torch.eye(m)[0]])  # a = |u|
    x = torch.cat([x, edge[None].to(V.dtype), -edge[None].to(V.dtype)])
    psd = torch.linalg.eigvalsh(V.to_matrix(x))[:, 0] >= -1e-12
    inside = V.in_cone(x)
    assert torch.equal(inside, psd)
    assert 0.2 < inside.float().mean() < 0.8 and bool(inside[-2]) and not bool(inside[-1])
    assert not V.in_cone(x[-2:-1], interior=True).any()
    if m >= 2:                                    # the Jordan product of two cone elements need not be in the cone
        p = V.element(1.0, torch.eye(m, dtype=V.dtype)[0])
        q = V.element(1.0, torch.eye(m, dtype=V.dtype)[1])
        assert V.in_cone(p) and V.in_cone(q) and not V.in_cone(V.product(p, q))
        assert V.in_cone(V.quadratic_representation(p, q))                         # P(p) preserves the cone


def test_spin_factor_batch_shapes_and_single_precision():
    gen = torch.Generator().manual_seed(70)
    V = jd.SpinFactor(4, dtype=torch.float32)
    x = torch.randn(2, 3, 5, generator=gen)
    assert V.exp(x).shape == (2, 3, 5) and V.idempotents(x).shape == (2, 3, 2, 5)
    assert V.to_matrix(x).shape == (2, 3, 4, 4)
    X = V.to_matrix(x)
    assert X.dtype == torch.complex64
    E = torch.linalg.matrix_exp(X)
    assert (V.to_matrix(V.exp(x)) - E).abs().max() < 2e-5 * E.abs().max()
    assert (V.eigenvalues(x).repeat_interleave(2, -1) - torch.linalg.eigvalsh(X)).abs().max() < 1e-5
    assert V.element(torch.zeros(3), torch.ones(3, 4)).shape == (3, 5)


@pytest.mark.parametrize('f', ['exp', 'sqrt', 'log', 'sin'])
def test_spin_factor_gradcheck(f):
    """Gradients of the functional calculus (away from u = 0, for sqrt/log inside the cone)."""
    gen = torch.Generator().manual_seed(80)
    V = jd.SpinFactor(3)
    x = torch.randn(4, 4, dtype=torch.float64, generator=gen)
    x[:, 0] = torch.linalg.vector_norm(x[:, 1:], dim=-1) + 0.5 + x[:, 0].abs()    # interior of the cone
    fn = {'exp': V.exp, 'sqrt': V.sqrt, 'log': V.log, 'sin': lambda y: V.apply_function(y, torch.sin)}[f]
    assert torch.autograd.gradcheck(fn, (x.requires_grad_(),))


def test_spin_factor_gradient_at_u_zero():
    """At u = 0 the functional calculus is smooth (divided difference = f'(a)); the gradient agrees with finite
    differences and with the gradient through the matrix representation."""
    V = jd.SpinFactor(3)
    x = torch.tensor([[0.3, 0.0, 0.0, 0.0]], dtype=torch.float64, requires_grad=True)
    assert torch.autograd.gradcheck(V.exp, (x,))
    W = torch.randn(4, dtype=torch.float64)
    g1, = torch.autograd.grad((V.exp(x) * W).sum(), x)
    xm = x.detach().clone().requires_grad_()
    E = torch.linalg.matrix_exp(V.to_matrix(xm))[0]
    g2, = torch.autograd.grad((V.from_matrix(E) * W).sum(), xm)
    assert (g1 - g2).abs().max() < 1e-12


def test_spin_factor_no_host_reads():
    V = jd.SpinFactor(5)
    x, y = torch.randn(8, 6, dtype=torch.float64), torch.randn(8, 6, dtype=torch.float64)

    def run():
        return (V.product(x, y), V.exp(x), V.sqrt(V.square(x)), V.idempotents(x), V.in_cone(x),
                V.quadratic_representation(x, y), V.to_matrix(x), V.inverse(x))
    run()
    with count_host_reads() as calls:
        run()
    assert not calls, sorted(set(calls))
