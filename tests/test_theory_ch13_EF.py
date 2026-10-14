"""
Theory tests for Chapter 13 (Jordan algebras), Appendix E (classification of Euclidean Jordan algebras) and Appendix F
(symmetric cones, order isomorphisms and Jordan maps) of the lecture notes, checked through torch_vn_algebra.jordan
(operator-level Jordan products, quadratic representation, transpose, generated Jordan / *-algebras, SpinFactor) and,
for the non-special and quaternionic algebras that have no matrix model in the library, through structure constants
built here from the Cayley-Dickson product (Section E.6). Everything is deterministic, double precision, CPU.
"""
import functools
import math

import numpy as np
import pytest
import scipy.integrate as si
import scipy.special as ss
import torch
from torch.autograd.functional import hessian, jacobian

from torch_vn_algebra import TypeIAlgebra, SpinFactor, wedderburn
from torch_vn_algebra import jordan as jd

F64, C128 = torch.float64, torch.complex128
LAYOUTS = [([1], [1]), ([3], [3]), ([2, 3], [2, 3]), ([4, 1, 3], [3, 1, 2]), ([2, 2, 3, 1], [2, 2, 3, 1])]


@pytest.fixture(autouse=True)
def _seed():
    torch.manual_seed(0)                                     # library samplers that use the global generator


def gen(seed):
    return torch.Generator().manual_seed(seed)


def herm_t(n, g, dtype=C128):
    A = torch.randn(n, n, dtype=dtype, generator=g)
    return A + A.mH


def mat(A):
    return A.matrix if hasattr(A, 'matrix') else A


def err(A, B):
    return (mat(A) - mat(B)).abs().max().item()


def jm(A, B):
    return 0.5 * (A @ B + B @ A)


def cm(A, B):
    return A @ B - B @ A


def sa_op(alg, B, g, positive=False):
    M = torch.randn(B, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype, generator=g)
    M = alg.operator(M).matrix                                     # zero padding
    if positive:
        return alg.operator(M @ M.mH, is_self_adjoint=True, is_positive=True)
    return alg.operator(M + M.mH, is_self_adjoint=True)


def unit_mat(n, i, j, dtype=C128):
    E = torch.zeros(n, n, dtype=dtype)
    E[i, j] = 1
    return E


def pauli():
    return jd.anticommuting_symmetries(3)                          # X, Y, Z


def bloch(A):
    """Bloch vector (a, u) of a Hermitian 2 x 2 matrix A = a 1 + u.sigma."""
    S = pauli()
    return torch.cat([torch.diagonal(A, dim1=-2, dim2=-1).sum(-1, keepdim=True).real / 2,
                      torch.einsum('kij,...ji->...k', S, A).real / 2], -1)


# ======================================================================================================================
# Composition algebras and Herm_n(D) by structure constants (Sections E.6, E.7)
# ======================================================================================================================
def conj(x):
    if x.shape[-1] == 1:
        return x
    h = x.shape[-1] // 2
    return torch.cat([conj(x[..., :h]), -x[..., h:]], -1)


def cd(x, y):
    """Cayley-Dickson product (a, b)(c, d) = (ac - d*b, da + bc*) on R^(2^k)."""
    if x.shape[-1] == 1:
        return x * y
    h = x.shape[-1] // 2
    a, b, c, d = x[..., :h], x[..., h:], y[..., :h], y[..., h:]
    return torch.cat([cd(a, c) - cd(conj(d), b), cd(d, a) + cd(b, conj(c))], -1)


@functools.lru_cache(None)
def mult_table(d):
    E = torch.eye(d, dtype=F64)
    return cd(E[:, None, :].expand(d, d, d), E[None, :, :].expand(d, d, d))       # T[s, t] = e_s e_t


def dmul(x, y):
    return torch.einsum('...s,...t,stu->...u', x, y, mult_table(x.shape[-1]))


def norm2(x):
    return (x * x).sum(-1)


class HermD:
    """Herm_n(D), D = R, C, H, O (d = 1, 2, 4, 8), coordinates: alpha_1..alpha_n, then a_ij in D for i < j."""

    def __init__(self, n, d):
        self.n, self.d = n, d
        self.pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        self.dim = n + d * len(self.pairs)
        B = torch.zeros(self.dim, n, n, d, dtype=F64)
        for i in range(n):
            B[i, i, i, 0] = 1
        unit = torch.eye(d, dtype=F64)
        for p, (i, j) in enumerate(self.pairs):
            for s in range(d):
                k = n + d * p + s
                B[k, i, j] = unit[s]
                B[k, j, i] = conj(unit[s])
        self.B = B
        X, Y = B[:, None], B[None, :]
        self.C = self.coords(0.5 * (self.matmul(X, Y) + self.matmul(Y, X)))       # C[a, b] = b_a o b_b
        self.e = self.coords(torch.eye(n, dtype=F64)[..., None] * torch.eye(d, dtype=F64)[0])
        self.G = torch.diag(torch.tensor([1.0] * n + [2.0] * (self.dim - n), dtype=F64))   # (X|Y) = Re tr(XY)

    def matmul(self, X, Y):
        return torch.einsum('...ijs,...jkt,stu->...iku', X, Y, mult_table(self.d))

    def coords(self, X):
        return torch.cat([X[..., range(self.n), range(self.n), 0]] + [X[..., i, j, :] for i, j in self.pairs], -1)

    def mat(self, v):
        return torch.einsum('...a,aijs->...ijs', v, self.B)

    def entry(self, a, i, j):
        """a[ij] = a E_ij + a* E_ji."""
        X = torch.zeros(self.n, self.n, self.d, dtype=F64)
        X[i, j], X[j, i] = a, conj(a)
        return self.coords(X)

    def diag(self, i):
        v = torch.zeros(self.dim, dtype=F64)
        v[i] = 1
        return v


class Spin:
    """V_m by structure constants, coordinates (a, u)."""

    def __init__(self, m):
        self.m, self.dim = m, m + 1
        V = SpinFactor(m)
        E = torch.eye(self.dim, dtype=F64)
        self.C = V.product(E[:, None].expand(self.dim, self.dim, self.dim), E[None].expand(self.dim, self.dim, self.dim))
        self.e = E[0]
        self.G = 2 * torch.eye(self.dim, dtype=F64)                 # tr(x o y) = 2(ab + u.v)


def mul(J, x, y):
    return torch.einsum('...a,...b,abc->...c', x, y, J.C)


def Lop(J, x):
    """Matrix of L(x) in the coordinates of J."""
    return torch.einsum('...a,abc->...cb', x, J.C)


def Pop(J, x):
    return 2 * Lop(J, x) @ Lop(J, x) - Lop(J, mul(J, x, x))


def trace_form(J):
    """tau(e_a, e_b) = Tr L(e_a o e_b)."""
    t = torch.einsum('cbb->c', Lop(J, torch.eye(J.dim, dtype=F64)))
    return torch.einsum('abc,c->ab', J.C, t)


def jpow(J, x, k):
    p = J.e.expand_as(x).clone() if k == 0 else x
    for _ in range(k - 1):
        p = mul(J, x, p)
    return p


def spectral(J, x, r):
    """Lemma E.3 / Theorem F.3 from Jordan powers only: the minimal polynomial t^r - sum_k a_k t^k of a generic x (as
    a monic polynomial killing e, x, ..., x^r), its roots lambda_i and the Lagrange idempotents c_i = l_i(x)."""
    pw = torch.stack([jpow(J, x, k) for k in range(r + 1)])            # e, x, ..., x^r
    A = pw[:r].T.numpy()
    coef, *_ = np.linalg.lstsq(A, pw[r].numpy(), rcond=None)
    resid = np.abs(A @ coef - pw[r].numpy()).max()
    lam = np.roots(np.concatenate([[1.0], -coef[::-1]]))
    cs = []
    for i in range(r):
        poly = np.array([1.0])
        for j in range(r):
            if j != i:
                poly = np.polymul(poly, np.array([1.0, -lam[j].real]) / (lam[i].real - lam[j].real))
        cvec = poly[::-1]                                              # coefficients of t^0, t^1, ...
        cs.append(sum(float(cvec[k]) * pw[k] for k in range(len(cvec))))
    return lam, torch.stack(cs), resid


@functools.lru_cache(None)
def albert():
    return HermD(3, 8)


def octonions(k, seed):
    return torch.randn(k, 8, dtype=F64, generator=gen(seed))


# ======================================================================================================================
# Section 13.1: the Jordan product and the quadratic representation
# ======================================================================================================================
@pytest.mark.parametrize('ns,ks', LAYOUTS)
@pytest.mark.parametrize('complex_', [True, False])
def test_lemma_13_1_jordan_identity_and_powers(ns, ks, complex_):
    """Lemma 13.1: the Jordan identity holds in M_sa, and the Jordan powers x∘x^k are the associative powers."""
    alg = TypeIAlgebra(ns, ks, complex_valued=complex_, precision='double')
    g = gen(1)
    x, y = sa_op(alg, 3, g), sa_op(alg, 3, g)
    x2 = jd.jordan_product(x, x)
    assert err(jd.jordan_product(jd.jordan_product(x2, y), x), jd.jordan_product(x2, jd.jordan_product(y, x))) < 1e-11
    assert err(x2, x.matrix @ x.matrix) < 1e-12
    p = x
    for k in range(1, 7):
        assert err(jd.jordan_power(x, k), torch.linalg.matrix_power(x.matrix, k)) < 1e-9 * 10 ** k
        assert err(p, jd.jordan_power(x, k)) < 1e-9 * 10 ** k
        p = jd.jordan_product(x, p)
    assert err(jd.jordan_power(x, 0), alg.identity(3)) == 0
    J = jd.jordan_product(x, y)
    assert err(J, J.matrix.mH) < 1e-12                                                 # closed in M_sa


@pytest.mark.parametrize('ns,ks', LAYOUTS)
def test_prop_13_2_quadratic_representation(ns, ks):
    """Proposition 13.2(1),(2): P(y)x = yxy = 2L(y)^2x - L(y^2)x, the fundamental formula, Tr P(y)x = Tr y^2 x."""
    alg = TypeIAlgebra(ns, ks, precision='double')
    g = gen(2)
    x, y, z = sa_op(alg, 4, g), sa_op(alg, 4, g), sa_op(alg, 4, g)
    jp, P = jd.jordan_product, jd.quadratic_representation
    Pyx = P(y, x)
    assert err(Pyx, y.matrix @ x.matrix @ y.matrix) < 1e-11
    assert err(Pyx, 2 * jp(y, jp(y, x)).matrix - jp(jp(y, y), x).matrix) < 1e-11
    assert err(P(Pyx, z), P(y, P(x, P(y, z)))) < 1e-8
    assert (Pyx.trace - (jp(y, y).matrix @ x.matrix).diagonal(dim1=-2, dim2=-1).sum((-2, -1))).abs().max() < 1e-10
    trip = jd.jordan_triple(x, y, z)                                                   # module docstring
    assert err(trip, jp(jp(x, y), z).matrix + jp(jp(y, z), x).matrix - jp(jp(x, z), y).matrix) < 1e-10
    assert err(jd.jordan_triple(x, y, x), P(x, y)) < 1e-11


@pytest.mark.parametrize('ns,ks', LAYOUTS)
def test_prop_13_2_3_positivity_and_inverse(ns, ks):
    """Proposition 13.2(3): P(y) maps the positive cone into itself for every self-adjoint y (indefinite ones
    included), and for invertible y it is a bijection of the cone with inverse P(y^{-1})."""
    alg = TypeIAlgebra(ns, ks, precision='double')
    g = gen(3)
    y, pos = sa_op(alg, 50, g), sa_op(alg, 50, g, positive=True)
    out = jd.quadratic_representation(y, pos)
    assert out.is_positive
    w = torch.cat([e for e in alg.operator(out.matrix, is_self_adjoint=True).eigenvalues() if e.numel()], -1)
    assert (w > -1e-10 * w.abs().max()).all()
    yinv = y.inverse()
    assert err(yinv.matrix @ y.matrix, alg.identity(1)) < 1e-8
    back = jd.quadratic_representation(yinv, out)
    assert err(back, pos) < 1e-6 * pos.matrix.abs().max().item()
    fwd = jd.quadratic_representation(y, jd.quadratic_representation(yinv, pos))
    assert err(fwd, pos) < 1e-6 * pos.matrix.abs().max().item()


def test_prop_13_2_4_jordan_product_of_positives_not_positive():
    """Proposition 13.2(4): for p = diag(1,0), q = (1/2)[[1,1],[1,1]], p∘q = (1/4)[[2,1],[1,0]] and det(p∘q) = -1/16."""
    M2 = TypeIAlgebra([2], [2], precision='double')
    p = M2.operator(torch.tensor([[[[1.0, 0], [0, 0]]]]), is_self_adjoint=True, is_positive=True)
    q = M2.operator(torch.tensor([[[[0.5, 0.5], [0.5, 0.5]]]]), is_self_adjoint=True, is_positive=True)
    r = jd.jordan_product(p, q)
    assert err(r.matrix[0, 0], torch.tensor([[2.0, 1], [1, 0]], dtype=C128) / 4) < 1e-15
    assert abs(torch.linalg.det(r.matrix[0, 0]).real.item() + 1 / 16) < 1e-15
    w = r.eigenvalues()[0][0]
    assert torch.allclose(w, torch.tensor([(1 - math.sqrt(2)) / 4, (1 + math.sqrt(2)) / 4], dtype=F64))
    assert not r.is_positive


def test_section_13_1_symmetrised_product_is_unique():
    """Section 13.1: among alpha xy + beta yx only alpha = beta = 1/2 is symmetric with x∘x = x^2 (XZ != ZX)."""
    X, Y, Z = pauli()
    assert err(X @ Z, Z @ X) > 1
    for a, b in ((0.3, 0.7), (1.0, 0.0), (0.5, 0.5)):
        sym = err(a * X @ Z + b * Z @ X, a * Z @ X + b * X @ Z) < 1e-15
        assert sym == (a == b)


# ======================================================================================================================
# Section 13.2: spin factors
# ======================================================================================================================
@pytest.mark.parametrize('m', range(1, 9))
def test_lemma_13_3_anticommuting_symmetries(m):
    """Lemma 13.3: s_i = s_i^*, s_i s_j + s_j s_i = 2 delta_ij, D = 2^floor(m/2) (2 for m = 1), Tr(s_i s_j) = D delta_ij,
    the Jordan-Wigner form of the strings, and s_{2j-1} s_{2j} = i Z_j."""
    S = jd.anticommuting_symmetries(m)
    D = 2 if m == 1 else 2 ** (m // 2)
    assert S.shape == (m, D, D)
    eye = torch.eye(D, dtype=C128)
    for i in range(m):
        assert err(S[i], S[i].mH) == 0 and abs(torch.trace(S[i]).item()) < 1e-14
        for j in range(m):
            assert err(S[i] @ S[j] + S[j] @ S[i], 2 * (i == j) * eye) < 1e-14
            assert abs(torch.trace(S[i] @ S[j]).item() - D * (i == j)) < 1e-12
    X, Y, Z = pauli()
    I2 = torch.eye(2, dtype=C128)
    p = m // 2

    def kron_all(fs):
        out = torch.ones(1, 1, dtype=C128)
        for f in fs:
            out = torch.kron(out, f)
        return out
    for j in range(1, p + 1):
        assert err(S[2 * j - 2], kron_all([Z] * (j - 1) + [X] + [I2] * (p - j))) == 0
        assert err(S[2 * j - 1], kron_all([Z] * (j - 1) + [Y] + [I2] * (p - j))) == 0
        assert err(S[2 * j - 2] @ S[2 * j - 1], 1j * kron_all([I2] * (j - 1) + [Z] + [I2] * (p - j))) < 1e-15
    if m % 2:
        assert err(S[-1], kron_all([Z] * p) if m > 1 else Z) == 0


@pytest.mark.parametrize('m', [1, 2, 3, 5, 6])
def test_lemma_13_3_phi_is_injective_jordan_homomorphism(m):
    """Lemma 13.3: phi(a,u) = a1 + sum u_i s_i satisfies phi(x∘y) = phi(x)∘phi(y), is injective, and the s_i generate
    M_D(C) as a *-algebra for m >= 2 (diagonal matrices for m = 1)."""
    V = SpinFactor(m)
    g = gen(4)
    x, y = V.random(100, generator=g), V.random(100, generator=g)
    X, Y = V.to_matrix(x), V.to_matrix(y)
    assert err(V.to_matrix(V.product(x, y)), jm(X, Y)) < 1e-12
    assert err(V.from_matrix(X), x) < 1e-13                                           # left inverse: injective
    D = V.representation_dim
    w = wedderburn(jd.anticommuting_symmetries(m))
    assert w.dim == (2 if m == 1 else D * D)


@pytest.mark.parametrize('p', [1, 2, 3, 4])
def test_section_13_2_majorana_operators(p):
    """Section 13.2: s_{2j-1} = (-1)^{j-1}(c_j + c_j^*), s_{2j} = (-1)^{j-1} i(c_j - c_j^*) for the Jordan-Wigner fermions
    c_j = (prod_{l<j}(-1)^{n_l}) s^-_j of Section 9.4, with (-1)^n = 1 - 2n = -Z and s^- = (X - iY)/2."""
    X, Y, Z = pauli()
    I2 = torch.eye(2, dtype=C128)
    n = 0.5 * Z + 0.5 * I2
    parity = I2 - 2 * n
    assert err(parity, -Z) == 0
    sminus = 0.5 * (X - 1j * Y)
    assert err(sminus, torch.tensor([[0, 0], [1, 0]], dtype=C128)) == 0
    S = jd.anticommuting_symmetries(2 * p)

    def kron_all(fs):
        out = torch.ones(1, 1, dtype=C128)
        for f in fs:
            out = torch.kron(out, f)
        return out
    cs = [kron_all([parity] * (j - 1) + [sminus] + [I2] * (p - j)) for j in range(1, p + 1)]
    for j, c in enumerate(cs, 1):
        sgn = (-1) ** (j - 1)
        assert err(S[2 * j - 2], sgn * (c + c.mH)) < 1e-15 and err(S[2 * j - 1], sgn * 1j * (c - c.mH)) < 1e-15
        for k, d in enumerate(cs, 1):                                                  # CAR (Theorem 9.10(a))
            assert err(c @ d.mH + d.mH @ c, (j == k) * torch.eye(2 ** p, dtype=C128)) < 1e-15


@pytest.mark.parametrize('m', [1, 2, 3, 4, 7])
def test_theorem_13_4_1_idempotents(m):
    """Theorem 13.4(1): x^2 = x iff x = 0, 1 or c_w = (1/2)(1, w) with |w| = 1; c_w∘c_w' = 0 iff w' = -w."""
    V = SpinFactor(m)
    g = gen(5)
    w = torch.nn.functional.normalize(torch.randn(200, m, dtype=F64, generator=g), dim=-1)
    c = V.element(0.5, 0.5 * w)
    assert err(V.square(c), c) < 1e-15
    assert err(V.product(c, V.element(0.5, -0.5 * w)), torch.zeros_like(c)) < 1e-15
    w2 = torch.nn.functional.normalize(torch.randn(200, m, dtype=F64, generator=g), dim=-1)
    prod = V.product(c, V.element(0.5, 0.5 * w2))
    assert torch.allclose(prod, 0.25 * torch.cat([1 + (w * w2).sum(-1, keepdim=True), w + w2], -1))
    if m > 1:
        assert (prod.abs().amax(-1) > 1e-6).all()
    for x in (V.unit() * 0, V.unit()):
        assert err(V.square(x), x) == 0
    a = torch.rand(200, dtype=F64, generator=g)                                         # |u| = 0.3 != 1/2, u != 0
    u = torch.nn.functional.normalize(torch.randn(200, m, dtype=F64, generator=g), dim=-1) * 0.3
    x = V.element(a, u)
    assert (V.square(x) - x).abs().amax(-1).min() > 1e-3                                 # |u| != 1/2: never idempotent


@pytest.mark.parametrize('m', [1, 2, 5, 8])
def test_theorem_13_4_2_3_spectral_theory(m):
    """Theorem 13.4(2),(3): x = lambda+ e+ + lambda- e-, lambda± = a ± |u|, x^k = lambda+^k e+ + lambda-^k e-,
    phi(f(x)) = f(phi(x)), eigenvalues of phi(x) are lambda± with multiplicity D/2, Tr phi(x) = (D/2) tr x,
    tr x = lambda+ + lambda-, det x = lambda+ lambda- = a^2 - |u|^2."""
    V = SpinFactor(m)
    g = gen(6)
    x = V.random(300, generator=g)
    a, u = V.split(x)
    lam, e = V.spectral_decomposition(x)
    assert torch.allclose(lam, torch.stack([a - u.norm(dim=-1), a + u.norm(dim=-1)], -1))
    assert err(lam[:, :1] * e[:, 0] + lam[:, 1:] * e[:, 1], x) < 1e-13
    for k in range(5):
        assert err(V.power(x, k), lam[:, :1] ** k * e[:, 0] + lam[:, 1:] ** k * e[:, 1]) < 1e-10 * 10 ** k
    assert torch.allclose(V.trace(x), lam.sum(-1)) and torch.allclose(V.determinant(x), lam.prod(-1))
    X = V.to_matrix(x)
    D = V.representation_dim
    assert err(torch.linalg.eigvalsh(X), lam.repeat_interleave(D // 2, dim=-1)) < 1e-12
    assert torch.allclose(torch.diagonal(X, dim1=-2, dim2=-1).sum(-1).real, D / 2 * V.trace(x))
    assert err(V.to_matrix(V.exp(x)), torch.linalg.matrix_exp(X)) < 1e-9 * torch.linalg.matrix_exp(X).abs().max()
    wX, UX = torch.linalg.eigh(X)
    assert err(V.to_matrix(V.apply_function(x, torch.abs)), UX @ torch.diag_embed(wX.abs().to(C128)) @ UX.mH) < 1e-11
    assert err(V.to_matrix(V.apply_function(x, torch.tanh)), UX @ torch.diag_embed(torch.tanh(wX).to(C128)) @ UX.mH) < 1e-11
    xz = V.element(torch.randn(5, dtype=F64, generator=g), torch.zeros(5, m, dtype=F64))   # u = 0: f(x) = f(a) 1
    assert err(V.exp(xz), V.element(torch.exp(V.split(xz)[0]), torch.zeros(5, m, dtype=F64))) < 1e-14


@pytest.mark.parametrize('m', [1, 2, 3, 6])
def test_theorem_13_4_4_lorentz_cone_is_self_dual(m):
    """Theorem 13.4(4): the cone of squares is {a >= |u|} = {lambda- >= 0} = {phi(x) >= 0}, and it is self-dual: for
    y not in Q, x = (|v|, -v) (or 1 if v = 0) is a square with tr(x∘y) < 0."""
    V = SpinFactor(m)
    g = gen(7)
    sq = V.random(500, generator=g, cone=True)
    a, u = V.split(sq)
    assert (a >= u.norm(dim=-1) - 1e-14).all() and V.in_cone(sq).all()
    assert (V.eigenvalues(sq)[:, 0] >= -1e-13).all()
    assert (torch.linalg.eigvalsh(V.to_matrix(sq))[:, 0] > -1e-12).all()
    assert err(V.square(V.sqrt(sq)), sq) < 1e-12 * sq.abs().max()                         # every cone point is a square
    y = V.random(2000, generator=g)
    b, v = V.split(y)
    inQ = V.in_cone(y, tol=0.0)
    assert torch.equal(inQ, b >= v.norm(dim=-1))
    clear = V.eigenvalues(y)[:, 0].abs() > 1e-9                                         # away from the boundary
    assert clear.float().mean() > 0.99
    assert torch.equal(inQ[clear], (torch.linalg.eigvalsh(V.to_matrix(y))[:, 0] >= 0)[clear])
    wit = V.element(v.norm(dim=-1), -v)
    assert V.in_cone(wit).all() and (V.inner(wit, y)[~inQ] < 0).all()
    assert (V.inner(sq[:, None], sq[None, :250]) >= -1e-12).all()                      # easy half
    assert V.inner(V.unit(), V.element(-1.0, torch.zeros(m, dtype=F64))) < 0


@pytest.mark.parametrize('m', [1, 3, 4])
def test_theorem_13_4_5_quadratic_representation(m):
    """Theorem 13.4(5): phi(P(y)x) = phi(y) phi(x) phi(y), hence P(y)Q ⊆ Q for every y."""
    V = SpinFactor(m)
    g = gen(8)
    x, y = V.random(300, generator=g), V.random(300, generator=g)
    Y = V.to_matrix(y)
    assert err(V.to_matrix(V.quadratic_representation(y, x)), Y @ V.to_matrix(x) @ Y) < 1e-11
    assert err(V.triple(y, x, y), V.quadratic_representation(y, x)) < 1e-12
    q = V.random(300, generator=g, cone=True)
    assert V.in_cone(V.quadratic_representation(y, q), tol=1e-12).all()


def test_prop_13_5_a_V3_is_Herm2C():
    """Proposition 13.5(a): for m = 3, phi(a, u) = a1 + u.sigma is a Jordan isomorphism V_3 ≅ Herm_2(C)."""
    V = SpinFactor(3)
    assert err(V.symmetries(), pauli()) == 0
    g = gen(9)
    H = torch.stack([herm_t(2, g) for _ in range(50)])
    assert err(V.to_matrix(V.from_matrix(H)), H) < 1e-14                              # onto
    x = V.from_matrix(H)
    assert err(V.to_matrix(V.product(x[:25], x[25:])), jm(H[:25], H[25:])) < 1e-13


@pytest.mark.parametrize('m', [2, 3, 5])
def test_prop_13_5_b_automorphisms_are_orthogonal(m):
    """Proposition 13.5(b): (a, u) -> (a, Ou) is a Jordan automorphism for O in O(m) (det ±1); a map (a, u) -> (a, Au)
    with A not orthogonal is not; and x^2 = 1 iff x = ±1 or x = (0, w), |w| = 1."""
    V = SpinFactor(m)
    g = gen(10)
    x, y = V.random(200, generator=g), V.random(200, generator=g)
    O = torch.linalg.qr(torch.randn(m, m, dtype=F64, generator=g))[0]
    for Om in (O, O @ torch.diag(torch.tensor([-1.0] + [1.0] * (m - 1), dtype=F64))):
        rot = lambda w: torch.cat([w[..., :1], w[..., 1:] @ Om.T], -1)
        assert err(rot(V.product(x, y)), V.product(rot(x), rot(y))) < 1e-12
    A = O + 0.1 * torch.randn(m, m, dtype=F64, generator=g)
    lin = lambda w: torch.cat([w[..., :1], w[..., 1:] @ A.T], -1)
    assert err(lin(V.product(x, y)), V.product(lin(x), lin(y))) > 1e-2
    w = torch.nn.functional.normalize(torch.randn(10, m, dtype=F64, generator=g), dim=-1)
    assert err(V.square(V.element(0.0, w)), V.unit((10,))) < 1e-15
    assert (V.square(V.element(0.5, w)) - V.unit((10,))).abs().amax(-1).min() > 0.1


@pytest.mark.parametrize('m', range(1, 9))
def test_prop_13_5_c_dimensions(m):
    """Proposition 13.5(c): dim phi(V_m) = m + 1, the generated *-algebra has self-adjoint part of dimension
    4^floor(m/2) (m >= 2; 2 for m = 1), and the two agree only for m = 1, 3."""
    r = jd.compare_with_star_algebra(jd.anticommuting_symmetries(m))
    assert r['jordan_dim'] == m + 1
    assert r['star_dim'] == (2 if m == 1 else 4 ** (m // 2))
    assert r['equal'] == (m in (1, 3))
    assert r['jordan'].is_self_adjoint_part() == (m in (1, 3))
    assert jd.is_jordan_subalgebra(r['jordan'].basis)


@pytest.mark.parametrize('m', [2, 3, 6])
def test_section_13_2_states_of_spin_factor(m):
    """Section 13.2 (states of V_m): omega(x) = tr(rho∘x) with rho = (1/2)(1, r) is a state iff |r| <= 1; omega(x) = a + u.r;
    every linear functional is tr(rho∘.) for exactly one rho (Riesz)."""
    V = SpinFactor(m)
    g = gen(11)
    r = torch.randn(400, m, dtype=F64, generator=g) * 0.8
    rho = V.element(0.5, 0.5 * r)
    q = V.random(300, generator=g, cone=True)
    vals = V.inner(rho[:, None], q[None])
    ok = r.norm(dim=-1) <= 1
    assert (vals[ok] >= -1e-12).all()
    assert torch.allclose(V.inner(rho, V.unit()), torch.ones(400, dtype=F64))
    bad = ~ok                                                                          # witness (|r|, -r) in Q
    assert (V.inner(rho[bad], V.element(r[bad].norm(dim=-1), -r[bad])) < 0).all()
    x = V.random(400, generator=g)
    assert torch.allclose(V.inner(rho, x), x[:, 0] + (x[:, 1:] * r).sum(-1))
    phi = torch.randn(m + 1, dtype=F64, generator=g)                                   # Riesz representation
    basis = torch.eye(m + 1, dtype=F64) / math.sqrt(2)                                 # orthonormal for tr(x∘y)
    rep = torch.einsum('k,kd->d', basis @ phi, basis)                                  # rho = sum_k omega(e_k) e_k
    assert torch.allclose(V.inner(rep, x), x @ phi)


# ======================================================================================================================
# Section 13.3: Jordan subalgebras that are not *-subalgebras
# ======================================================================================================================
@pytest.mark.parametrize('n', range(1, 6))
def test_prop_13_6_a_real_symmetric(n):
    """Proposition 13.6(a): two generic real symmetric matrices generate Herm_n(R) (dim n(n+1)/2) as a Jordan algebra and
    M_n(C) (self-adjoint part of dim n^2) as a *-algebra; E_ij = E_ii(E_ij + E_ji)."""
    g = gen(12 + n)
    G = torch.randn(2, n, n, dtype=F64, generator=g)
    r = jd.compare_with_star_algebra(G + G.mT)
    assert (r['jordan_dim'], r['star_dim']) == (n * (n + 1) // 2, n * n)
    assert r['jordan'].field == 'real'
    E = lambda i, j: unit_mat(n, i, j, F64)
    for i in range(n):
        for j in range(n):
            if i != j:
                assert torch.equal(E(i, i) @ (E(i, j) + E(j, i)), E(i, j))


def test_prop_13_6_b_two_hermitian_2x2():
    """Proposition 13.6(b): two Hermitian 2x2 matrices are simultaneously unitarily real symmetric; they generate the
    Jordan algebra span{1, x, y} ≅ V_2 (dim 3) but the *-algebra M_2(C); three generic ones generate Herm_2(C)."""
    g = gen(13)
    for _ in range(5):
        H = torch.stack([herm_t(2, g) for _ in range(3)])
        r2, r3 = jd.compare_with_star_algebra(H[:2]), jd.compare_with_star_algebra(H)
        assert (r2['jordan_dim'], r2['star_dim'], r3['jordan_dim'], r3['star_dim']) == (3, 4, 4, 4)
        span = torch.stack([torch.eye(2, dtype=C128), H[0], H[1]])
        assert jd.is_jordan_subalgebra(span)
        assert bool(r2['jordan'].contains(jm(H[0], H[1]))) and bool(r2['jordan'].contains(H[0] @ H[0]))
        assert not bool(r2['jordan'].contains(0.5j * cm(H[0], H[1])))
        u, v = bloch(H[0])[1:], bloch(H[1])[1:]
        nrm = torch.linalg.cross(u, v)
        nrm = nrm / nrm.norm()
        _, W1 = torch.linalg.eigh(torch.einsum('k,kij->ij', nrm.to(C128), pauli()))
        _, W2 = torch.linalg.eigh(pauli()[1])
        U = W2 @ W1.mH
        assert err(U @ torch.einsum('k,kij->ij', nrm.to(C128), pauli()) @ U.mH, pauli()[1]) < 1e-12
        R = torch.stack([bloch(U @ s @ U.mH)[1:] for s in pauli()], 1)               # R in O(3), Rn = e_2
        assert torch.allclose(R.T @ R, torch.eye(3, dtype=F64))
        assert torch.allclose(R @ nrm, torch.tensor([0., 1, 0], dtype=F64))
        assert all((U @ A @ U.mH).imag.abs().max() < 1e-12 for A in H[:2])
        xy = jm(H[0], H[1])                                                             # x∘y = (ab + u.v)1 + a v.s + b u.s
        a, b = bloch(H[0])[0], bloch(H[1])[0]
        exp = (a * b + u @ v) * torch.eye(2, dtype=C128) + torch.einsum('k,kij->ij', (a * v + b * u).to(C128), pauli())
        assert err(xy, exp) < 1e-12


@pytest.mark.parametrize('ks,ms', [((1, 2, 3), (2, 1, 2)), ((3,), (1,)), ((2, 2), (1, 3)), ((4, 1), (1, 2))])
def test_section_13_3_jordan_vs_star_in_blocks(ks, ms):
    """Section 13.3: in A = U(⊕ M_{k_c} ⊗ 1_{m_c})U^*, three generic self-adjoint elements generate J = A_sa, two
    generate only Herm_2(R) in each block M_2 (so J != A_sa iff some k_c = 2); always J ⊆ A_sa."""
    g = gen(14)
    n = sum(k * m for k, m in zip(ks, ms))
    U = torch.linalg.qr(torch.randn(n, n, dtype=C128, generator=g))[0]

    def element():
        return U @ torch.block_diag(*[torch.kron(herm_t(k, g), torch.eye(m, dtype=C128)) for k, m in zip(ks, ms)]) @ U.mH
    G = [element() for _ in range(3)]
    r3, r2 = jd.compare_with_star_algebra(G), jd.compare_with_star_algebra(G[:2])
    full = sum(k * k for k in ks)
    assert r3['jordan_dim'] == r3['star_dim'] == full == wedderburn(G).dim
    assert r2['star_dim'] == full
    assert r2['jordan_dim'] == full - sum(1 for k in ks if k == 2)
    assert r2['equal'] == (2 not in ks)
    assert bool(r2['star'].contains(r2['jordan'].basis).all())


# ======================================================================================================================
# Section 13.4: Euclidean Jordan algebras, the classification, symmetric cones
# ======================================================================================================================
@pytest.mark.parametrize('ns,ks', LAYOUTS)
def test_section_13_4_M_sa_is_euclidean(ns, ks):
    """Section 13.4: M_sa with Tr(xy) is Euclidean, Tr((x∘y)z) = Tr(y(x∘z)); formally real: Tr(sum x_i^2) = sum ||x_i||^2."""
    alg = TypeIAlgebra(ns, ks, precision='double')
    g = gen(15)
    x, y, z = sa_op(alg, 5, g), sa_op(alg, 5, g), sa_op(alg, 5, g)
    jp = jd.jordan_product
    t1 = (jp(x, y).matrix @ z.matrix).diagonal(dim1=-2, dim2=-1).sum((-2, -1))
    t2 = (y.matrix @ jp(x, z).matrix).diagonal(dim1=-2, dim2=-1).sum((-2, -1))
    assert (t1 - t2).abs().max() < 1e-10 and t1.imag.abs().max() < 1e-10
    s = jp(x, x).trace + jp(y, y).trace
    assert torch.allclose(s.real, x.frobenius_norm() ** 2 + y.frobenius_norm() ** 2)


@pytest.mark.parametrize('m', [1, 2, 4])
def test_section_13_4_spin_factor_is_euclidean(m):
    """Section 13.4: tr((x∘y)∘z) = 2(abc + c u.v + a v.w + b u.w) is symmetric in x, y, z (V_m is Euclidean)."""
    V = SpinFactor(m)
    g = gen(16)
    x, y, z = (V.random(100, generator=g) for _ in range(3))
    t = V.trace(V.product(V.product(x, y), z))
    (a, u), (b, v), (c, w) = V.split(x), V.split(y), V.split(z)
    assert torch.allclose(t, 2 * (a * b * c + c * (u * v).sum(-1) + a * (v * w).sum(-1) + b * (u * w).sum(-1)))
    for perm in ((y, x, z), (z, y, x), (x, z, y)):
        assert torch.allclose(t, V.trace(V.product(V.product(*perm[:2]), perm[2])))


@pytest.mark.parametrize('n,d,dim', [(1, 1, 1), (3, 1, 6), (3, 2, 9), (3, 4, 15), (3, 8, 27), (2, 8, 10), (4, 4, 28)])
def test_theorem_13_7_dimensions_and_rank_two(n, d, dim):
    """Theorem 13.7: dim Herm_r(K) = r + d r(r-1)/2 (Albert algebra: 27), and Herm_2(K) ≅ V_{d+1}: in rank two the
    product of the off-diagonal part is x∘y = <x, y>1 (Theorem E.13)."""
    J = HermD(n, d)
    assert J.dim == dim == n + d * n * (n - 1) // 2
    if n == 2:
        g = gen(17)
        f = J.diag(0) - J.diag(1)
        x, y = torch.randn(2, d, dtype=F64, generator=g)
        X, Y = J.entry(x, 0, 1), J.entry(y, 0, 1)
        assert torch.allclose(mul(J, f, f), J.e) and torch.allclose(mul(J, f, X), torch.zeros(J.dim, dtype=F64))
        assert torch.allclose(mul(J, X, Y), (x @ y) * J.e)


@pytest.mark.parametrize('d', [1, 2, 4, 8])
def test_theorem_E_13_rank_two_isomorphism(d):
    """Theorem E.13 / E.24(3): Herm_2(D) -> V_{d+1}, alpha_1 E_11 + alpha_2 E_22 + a[12] -> ((alpha_1+alpha_2)/2,
    ((alpha_1-alpha_2)/2, a)) is a Jordan isomorphism (f = E_11 - E_22, W = Rf + J_12, w∘w' = B(w,w')1)."""
    J, V = HermD(2, d), SpinFactor(d + 1)
    iso = lambda v: torch.cat([(v[..., :1] + v[..., 1:2]) / 2, (v[..., :1] - v[..., 1:2]) / 2, v[..., 2:]], -1)
    X, Y = torch.randn(2, 100, J.dim, dtype=F64, generator=gen(170 + d))
    assert torch.allclose(iso(mul(J, X, Y)), V.product(iso(X), iso(Y)))
    assert torch.allclose(iso(J.e), V.unit())


def test_section_13_4_cone_of_direct_sum_is_product():
    """Section 13.4: the cone of squares of a direct sum is the product of the cones ((x_1, ..., x_k)^2 = (x_k^2)), and
    for real x >= 0 the positive square root is real (the semidefinite cones of Herm_n(R) are cones of squares)."""
    alg = TypeIAlgebra([2, 3, 1], [2, 3, 1], complex_valued=False, precision='double')
    g = gen(18)
    p = sa_op(alg, 20, g, positive=True)
    r = p.sqrt()
    assert r.matrix.dtype == F64 and err(r.matrix @ r.matrix, p) < 1e-12
    blk = [w for w in r.eigenvalues()]
    assert all((w >= -1e-12).all() for w in blk)
    x = sa_op(alg, 20, g)
    sq = jd.jordan_product(x, x)
    assert all((w >= -1e-10).all() for w in sq.eigenvalues())


# ======================================================================================================================
# Section 13.5: the transpose, Kadison, composite systems, dynamics
# ======================================================================================================================
@pytest.mark.parametrize('ns,ks', LAYOUTS)
def test_prop_13_9_a_transpose_is_jordan_star_antiautomorphism(ns, ks):
    """Proposition 13.9(a): T is a unital Jordan *-automorphism with T(AB) = T(B)T(A), preserves spectra and P."""
    alg = TypeIAlgebra(ns, ks, precision='double')
    g = gen(19)
    x, y = sa_op(alg, 4, g), sa_op(alg, 4, g)
    A = alg.operator(torch.randn(4, alg.C, alg.k_max, alg.k_max, dtype=C128, generator=g))
    T = jd.transpose_map
    assert err(T(alg.identity(1)), alg.identity(1)) == 0
    assert err(T(jd.jordan_product(x, y)), jd.jordan_product(T(x), T(y))) < 1e-12
    assert err(T(jd.quadratic_representation(y, x)), jd.quadratic_representation(T(y), T(x))) < 1e-12
    assert err(T(alg.operator(A.matrix.mH)), T(A).matrix.mH) == 0 and err(T(T(A)), A) == 0
    assert err(T(x @ y), T(y) @ T(x)) < 1e-12
    if max(ks) >= 2:
        assert err(T(x @ y), T(x) @ T(y)) > 1e-3
    w1 = torch.cat([e for e in x.eigenvalues() if e.numel()], -1)
    w2 = torch.cat([e for e in alg.operator(T(x).matrix, is_self_adjoint=True).eigenvalues() if e.numel()], -1)
    assert err(w1, w2) < 1e-12


@pytest.mark.parametrize('k', [1, 2, 3, 4])
def test_prop_13_9_b_choi_matrix_is_swap(k):
    """Proposition 13.9(b): the Choi matrix of T on M_k is the swap F, with eigenvalues +1 (k(k+1)/2 times) and -1
    (k(k-1)/2 times), so T is not completely positive for k >= 2."""
    Mk = TypeIAlgebra([k], [k], precision='double')
    units = torch.eye(k * k, dtype=C128).reshape(k * k, 1, k, k)
    TE = jd.transpose_map(Mk.operator(units)).matrix[:, 0]
    choi = TE.reshape(k, k, k, k).permute(0, 2, 1, 3).reshape(k * k, k * k)
    swap = torch.eye(k * k, dtype=C128).reshape(k, k, k, k).permute(0, 1, 3, 2).reshape(k * k, k * k)
    assert torch.equal(choi, swap)
    w = torch.linalg.eigvalsh(choi)
    assert int((w < 0).sum()) == k * (k - 1) // 2 and int((w > 0).sum()) == k * (k + 1) // 2
    assert torch.allclose(w.abs(), torch.ones(k * k, dtype=F64))


def rot_of(U):
    """R_U in O(3): U (u.sigma) U^* = (R_U u).sigma."""
    return torch.stack([bloch(U @ s @ U.mH)[1:] for s in pauli()], 1)


def test_prop_13_11_jordan_automorphisms_of_Herm2():
    """Proposition 13.11: R_U is in SO(3), R_{UV} = R_U R_V, O_T = diag(1,-1,1); every O in O(3) is R_U (det 1) or
    R_U O_T (det -1), constructed as in the proof (U_1 U_2 with U_2 = diag(e^{-i beta/2}, e^{i beta/2}))."""
    g = gen(20)
    X, Y, Z = pauli()
    S = pauli()
    O_T = torch.stack([bloch(s.T)[1:] for s in S], 1)
    assert torch.equal(O_T, torch.diag(torch.tensor([1.0, -1.0, 1.0], dtype=F64)))
    for _ in range(10):
        U, V = (torch.linalg.qr(torch.randn(2, 2, dtype=C128, generator=g))[0] for _ in range(2))
        assert abs(torch.linalg.det(rot_of(U)).item() - 1) < 1e-12
        assert torch.allclose(rot_of(U @ V), rot_of(U) @ rot_of(V))
        O = torch.linalg.qr(torch.randn(3, 3, dtype=F64, generator=g))[0]
        target = O if torch.linalg.det(O) > 0 else O @ O_T
        n3 = target[:, 2]                                                              # U_1 Z U_1^* = (R e_3).sigma
        _, W1 = torch.linalg.eigh(torch.einsum('k,kij->ij', n3.to(C128), S))
        _, WZ = torch.linalg.eigh(Z)
        U1 = W1 @ WZ.mH
        assert torch.allclose(rot_of(U1)[:, 2], n3)
        Rz = rot_of(U1).T @ target                                                     # rotation about e_3
        beta = math.atan2(Rz[1, 0].item(), Rz[0, 0].item())
        U2 = torch.diag(torch.tensor([complex(math.cos(-beta / 2), math.sin(-beta / 2)),
                                      complex(math.cos(beta / 2), math.sin(beta / 2))], dtype=C128))
        assert err(U2 @ X @ U2.mH, math.cos(beta) * X + math.sin(beta) * Y) < 1e-14
        assert err(U2 @ Y @ U2.mH, -math.sin(beta) * X + math.cos(beta) * Y) < 1e-14
        Uf = U1 @ U2
        assert torch.allclose(rot_of(Uf), target)
        A = herm_t(2, g)
        a = bloch(A)
        alpha = a[0] * torch.eye(2, dtype=C128) + torch.einsum('k,kij->ij', (O @ a[1:]).to(C128), S)
        form = Uf @ A @ Uf.mH if torch.linalg.det(O) > 0 else Uf @ A.T @ Uf.mH
        assert err(alpha, form) < 1e-12


def ptranspose(R, n, m):
    return R.reshape(*R.shape[:-2], n, m, n, m).transpose(-3, -1).reshape(*R.shape[:-2], n * m, n * m)


@pytest.mark.parametrize('n,m', [(2, 2), (2, 3), (3, 2)])
def test_prop_13_12_composite_systems(n, m):
    """Proposition 13.12: (x⊗y)∘(x'⊗y') = (x∘x')⊗(y∘y') + (1/4)[x,x']⊗[y,y']; ∘_theta has the opposite sign, agrees with ∘
    on local observables, is the Jordan product of A⋆B = theta(theta(A)theta(B)) with (x⊗y)⋆(x'⊗y') = xx'⊗y'y, and its
    cone of squares theta(Herm^+) does not contain the Bell state."""
    g = gen(21)
    x, x2, y, y2 = herm_t(n, g), herm_t(n, g), herm_t(m, g), herm_t(m, g)
    kx = torch.kron
    th = lambda R: ptranspose(R, n, m)
    loc = kx(jm(x, x2), jm(y, y2))
    com = 0.25 * kx(cm(x, x2), cm(y, y2))
    assert err(jm(kx(x, y), kx(x2, y2)), loc + com) < 1e-10
    jth = lambda A, B: th(jm(th(A), th(B)))
    assert err(jth(kx(x, y), kx(x2, y2)), loc - com) < 1e-10
    assert err(th(kx(x, y)), kx(x, y.T.contiguous())) == 0 and err(th(th(kx(x, y))), kx(x, y)) == 0
    In, Im = torch.eye(n, dtype=C128), torch.eye(m, dtype=C128)
    for A, B in ((kx(x, Im), kx(In, y)), (kx(x, Im), kx(x2, Im)), (kx(In, y), kx(In, y2))):
        assert err(jth(A, B), jm(A, B)) < 1e-10
    star = lambda A, B: th(th(A) @ th(B))
    assert err(star(kx(x, y), kx(x2, y2)), kx(x @ x2, y2 @ y)) < 1e-10
    W = herm_t(n * m, g)
    W2 = herm_t(n * m, g)
    assert err(jth(W, W2), 0.5 * (star(W, W2) + star(W2, W))) < 1e-10
    sq = jth(W, W)                                                                     # theta(theta(W)^2): in theta(Herm^+)
    assert torch.linalg.eigvalsh(th(sq))[0] > -1e-10
    if (n, m) == (2, 2):
        Om = torch.tensor([1.0, 0, 0, 1.0], dtype=C128)
        omega = 0.5 * torch.outer(Om, Om)
        F = torch.eye(4, dtype=C128).reshape(2, 2, 2, 2).permute(0, 1, 3, 2).reshape(4, 4)
        assert err(th(omega), 0.5 * F) == 0
        assert abs(torch.linalg.eigvalsh(th(omega))[0].item() + 0.5) < 1e-14


def test_prop_13_12_separable_states_are_ppt():
    """Section 13.5 (PPT): a separable state sum p_k a_k ⊗ b_k lies in both cones (theta(rho) >= 0)."""
    g = gen(22)
    rho = torch.zeros(6, 6, dtype=C128)
    for _ in range(5):
        A, B = torch.randn(2, 2, dtype=C128, generator=g), torch.randn(3, 3, dtype=C128, generator=g)
        rho += torch.rand(1, dtype=F64, generator=g).item() * torch.kron(A @ A.mH, B @ B.mH)
    rho /= torch.trace(rho)
    assert torch.linalg.eigvalsh(rho)[0] > -1e-14 and torch.linalg.eigvalsh(ptranspose(rho, 2, 3))[0] > -1e-14


@pytest.mark.parametrize('p', [0.0, 0.2, 1 / 3, 0.5, 1.0])
def test_exercise_13_2_isotropic_states(p):
    """Exercise 13.2: theta(rho_p) = (p/2)F + (1-p)/4 has spectrum (1+p)/4 (x3) and (1-3p)/4: entangled iff p > 1/3."""
    Om = torch.tensor([1.0, 0, 0, 1.0], dtype=C128)
    rho = p * 0.5 * torch.outer(Om, Om) + (1 - p) * torch.eye(4, dtype=C128) / 4
    w = torch.linalg.eigvalsh(ptranspose(rho, 2, 2))
    exp = torch.tensor(sorted([(1 - 3 * p) / 4] + [(1 + p) / 4] * 3), dtype=F64)
    assert err(w, exp) < 1e-14
    assert (w[0] < -1e-15) == (p > 1 / 3 + 1e-12)


def test_prop_13_13_dynamics_from_jordan_product():
    """Proposition 13.13: delta = i[H, .] maps Herm_n to itself and is a derivation of ∘; e^{t delta} z = e^{itH} z e^{-itH}
    is a Jordan automorphism; [L(x), L(y)]z = (1/4)[[x,y],z]; delta = 4 sum_k [L(x_k), L(y_k)] with the x_k, y_k of the
    proof; T(i[H, T(z)]) = -i[H^T, z]."""
    n = 4
    g = gen(23)
    H, z1, z2 = herm_t(n, g), herm_t(n, g), herm_t(n, g)
    delta = lambda A: 1j * cm(H, A)
    assert err(delta(z1), delta(z1).mH) < 1e-12
    assert err(delta(jm(z1, z2)), jm(delta(z1), z2) + jm(z1, delta(z2))) < 1e-10
    t = 0.37
    Ut = torch.linalg.matrix_exp(1j * t * H)
    flow = lambda A: Ut @ A @ Ut.mH
    gen_mat = torch.stack([delta(unit_mat(n, i, j)).reshape(-1) for i in range(n) for j in range(n)], 1)   # delta on M_n
    et = torch.linalg.matrix_exp(t * gen_mat)
    assert err((et @ z1.reshape(-1)).reshape(n, n), flow(z1)) < 1e-11
    assert err(flow(jm(z1, z2)), jm(flow(z1), flow(z2))) < 1e-10
    L = lambda x: (lambda w: jm(x, w))
    assert err(L(z1)(L(z2)(H)) - L(z2)(L(z1)(H)), 0.25 * cm(cm(z1, z2), H)) < 1e-10
    pairs = []
    Eu = lambda i, j: unit_mat(n, i, j)
    for i in range(n):
        for j in range(i + 1, n):
            a, b = Eu(i, j) + Eu(j, i), 1j * (Eu(i, j) - Eu(j, i))
            assert err(cm(a, Eu(i, i)), Eu(j, i) - Eu(i, j)) == 0
            assert err(cm(b, Eu(i, i)), -1j * (Eu(i, j) + Eu(j, i))) == 0
            assert err(cm(a, b), 2j * (Eu(j, j) - Eu(i, i))) == 0
            p_, q_ = H[i, j].real.item(), H[i, j].imag.item()
            pairs += [(-p_ * b, Eu(i, i)), (q_ * a, Eu(i, i))]
    for j in range(1, n):
        a, b = Eu(0, j) + Eu(j, 0), 1j * (Eu(0, j) - Eu(j, 0))
        pairs.append((0.5 * (H[j, j] - torch.trace(H) / n).real.item() * a, b))
    Hp = H - torch.trace(H) / n * torch.eye(n, dtype=C128)
    assert err(sum(cm(x, y) for x, y in pairs), 1j * Hp) < 1e-12
    for z in (z1, z2):
        assert err(4 * sum(L(x)(L(y)(z)) - L(y)(L(x)(z)) for x, y in pairs), delta(z)) < 1e-10
    assert err(1j * cm(H, z1.T).T, -1j * cm(H.T, z1)) < 1e-12


def test_exercise_13_4_real_quantum_mechanics():
    """Exercise 13.4: in Herm_n(R), [L(x),L(y)] = [K, .] with K = [x,y]/4 real antisymmetric, never i[H, .] with H real."""
    g = gen(24)
    n = 4
    x, y, z = (herm_t(n, g, F64) for _ in range(3))
    K = 0.25 * cm(x, y)
    assert err(K.T, -K) < 1e-14
    assert err(jm(x, jm(y, z)) - jm(y, jm(x, z)), cm(K, z)) < 1e-12
    Hr = herm_t(n, g, F64)
    i_comm = 1j * cm(Hr.to(C128), z.to(C128))
    assert i_comm.real.abs().max() == 0 and cm(K, z).abs().max() > 0


# ======================================================================================================================
# Section 13.6: trace inequality
# ======================================================================================================================
@pytest.mark.parametrize('m', [1, 2, 4, 7])
def test_prop_13_14_trace_inequality_in_spin_factor(m):
    """Proposition 13.14: tr|P(y)x| <= tr P(y)|x| in V_m, with equality if x or -x is in Q; via phi this is
    Tr|YXY| <= Tr(Y|X|Y) for Y = Y^* (Proposition 11.9 without Y >= 0)."""
    V = SpinFactor(m)
    g = gen(25)
    x, y = V.random(4000, generator=g), V.random(4000, generator=g)
    absf = lambda w: V.apply_function(w, torch.abs)
    lhs = V.trace(absf(V.quadratic_representation(y, x)))
    rhs = V.trace(V.quadratic_representation(y, absf(x)))
    assert (rhs - lhs >= -1e-12 * rhs.abs().clamp(min=1)).all()
    definite = V.determinant(x) >= 0
    assert definite.any() and (~definite).any()
    assert ((rhs - lhs)[definite].abs() <= 1e-11 * rhs[definite].abs().clamp(min=1)).all()
    if m == 1:                                                                          # V_1 = R + R is associative
        assert ((rhs - lhs).abs() <= 1e-11 * rhs.abs().clamp(min=1)).all()
    else:
        assert ((rhs - lhs)[~definite] > 1e-9).float().mean() > 0.9
    X, Y = V.to_matrix(x[:300]), V.to_matrix(y[:300])
    D = V.representation_dim
    assert torch.allclose(torch.linalg.svdvals(Y @ X @ Y).sum(-1), D / 2 * lhs[:300])
    wx, Ux = torch.linalg.eigh(X)
    absX = Ux @ torch.diag_embed(wx.abs().to(C128)) @ Ux.mH
    assert torch.allclose(V.to_matrix(absf(x[:300])), absX)                             # phi(|x|) = |phi(x)|
    assert torch.allclose((Y @ absX @ Y).diagonal(dim1=-2, dim2=-1).sum(-1).real, D / 2 * rhs[:300])


@pytest.mark.parametrize('ns,ks', LAYOUTS)
def test_section_13_6_prop_11_9_for_self_adjoint_Y(ns, ks):
    """Section 13.6 / Corollary F.7 on M_sa: Tr|YXY| <= Tr(Y|X|Y) for all Y = Y^* (not only Y >= 0), with equality when
    X >= 0, X <= 0, or XY = YX."""
    alg = TypeIAlgebra(ns, ks, precision='double')
    g = gen(26)
    x, y = sa_op(alg, 200, g), sa_op(alg, 200, g)
    P = jd.quadratic_representation
    lhs, rhs = P(y, x).trace_norm(), P(y, x.abs()).trace.real
    assert (rhs - lhs >= -1e-10 * rhs).all()
    if max(ks) >= 2:
        assert (y.lambda_min < 0).any() and (rhs - lhs > 1e-6).any()
    pos = sa_op(alg, 50, g, positive=True)
    y2 = sa_op(alg, 50, g)
    assert torch.allclose(P(y2, pos).trace_norm(), P(y2, pos.abs()).trace.real)
    comm = alg.operator(x.matrix[:50] @ x.matrix[:50] - 0.3 * x.matrix[:50], is_self_adjoint=True)   # commutes with x
    x50 = alg.operator(x.matrix[:50], is_self_adjoint=True)
    assert torch.allclose(P(comm, x50).trace_norm(), P(comm, x50.abs()).trace.real)


# ======================================================================================================================
# Exercises 13.1, 13.3
# ======================================================================================================================
@pytest.mark.parametrize('m', [2, 4])
def test_exercise_13_1_homogeneity_of_lorentz_cone(m):
    """Exercise 13.1: y^{-1} = (a,-u)/det y, P(y)P(y^{-1}) = id, P(y) maps int Q onto itself, P(x^{1/2})1 = x, and
    P(x'^{1/2})P(x^{1/2})^{-1} maps x to x'."""
    V = SpinFactor(m)
    g = gen(27)
    y = V.random(200, generator=g)
    assert err(V.product(y, V.inverse(y)), V.unit((200,))) < 1e-9
    x = V.random(200, generator=g)
    back = V.quadratic_representation(y, V.quadratic_representation(V.inverse(y), x))
    assert err(back, x) < 1e-7 * x.abs().max()
    xi = V.random(200, generator=g, cone=True)
    assert V.in_cone(xi, interior=True).all()
    assert V.in_cone(V.quadratic_representation(y, xi), interior=True).all()
    r = V.sqrt(xi)
    assert err(V.quadratic_representation(r, V.unit((200,))), xi) < 1e-10 * xi.abs().max()
    x2 = V.random(200, generator=g, cone=True)
    mapped = V.quadratic_representation(V.sqrt(x2), V.quadratic_representation(V.inverse(r), xi))
    assert err(mapped, x2) < 1e-8 * x2.abs().max()


@pytest.mark.parametrize('m', [2, 3, 6])
def test_exercise_13_3_products_of_boundary_points(m):
    """Exercise 13.3: (1,w)∘(1,w') = (1 + t, w + w') is in Q iff t = w.w' = ±1."""
    V = SpinFactor(m)
    g = gen(28)
    w, w2 = torch.nn.functional.normalize(torch.randn(2, 500, m, dtype=F64, generator=g), dim=-1)
    one = lambda v: V.element(1.0, v)
    prod = V.product(one(w), one(w2))
    t = (w * w2).sum(-1)
    assert torch.allclose(prod, torch.cat([(1 + t)[:, None], w + w2], -1))
    assert not V.in_cone(prod).any()
    assert V.in_cone(V.product(one(w), one(w))).all()
    assert err(V.product(one(w), one(-w)), torch.zeros_like(prod)) < 1e-15


# ======================================================================================================================
# Appendix E: composition algebras
# ======================================================================================================================
def test_section_E_6_cayley_dickson_tower():
    """Section E.6: CD gives R, C, H, O of dimensions 1, 2, 4, 8; in C, (a,b)(c,d) = (ac - bd, ad + bc); the unit is (1,0)
    and the conjugation (a,b)* = (a*, -b) is 2Re(x)1 - x."""
    g = gen(29)
    z = torch.randn(2, 100, 2, dtype=F64, generator=g)
    zc = torch.complex(z[..., 0], z[..., 1])
    p = cd(z[0], z[1])
    assert torch.allclose(torch.complex(p[..., 0], p[..., 1]), zc[0] * zc[1])
    for d in (1, 2, 4, 8, 16):
        x = torch.randn(50, d, dtype=F64, generator=g)
        one = torch.zeros(d, dtype=F64)
        one[0] = 1
        assert torch.allclose(cd(one.expand_as(x), x), x) and torch.allclose(cd(x, one.expand_as(x)), x)
        assert torch.allclose(conj(x), 2 * x[..., :1] * one - x)


@pytest.mark.parametrize('d', [1, 2, 4, 8])
def test_lemma_E_14_composition_identities(d):
    """Lemma E.14 in R, C, H, O: N(xy) = N(x)N(y), <xy,xz> = N(x)<y,z>, <xy,z> = <y, x*z> = <x, zy*>,
    x*(xy) = N(x)y = (yx)x*,
    x^2 - 2Re(x)x + N(x) = 0, conj(xy) = y*x*, Re(xy) = Re(yx), Re((xy)z) = Re(x(yz)), x(xy) = x^2y, (yx)x = yx^2,
    and the associator is alternating."""
    g = gen(30 + d)
    x, y, z, w = torch.randn(4, 500, d, dtype=F64, generator=g)
    ip = lambda p, q: (p * q).sum(-1)
    N = norm2
    one = torch.zeros(d, dtype=F64)
    one[0] = 1
    assert torch.allclose(N(dmul(x, y)), N(x) * N(y))
    assert torch.allclose(ip(dmul(x, y), dmul(x, z)), N(x) * ip(y, z))
    assert torch.allclose(ip(dmul(y, x), dmul(z, x)), ip(y, z) * N(x))
    assert torch.allclose(ip(dmul(x, y), dmul(w, z)) + ip(dmul(w, y), dmul(x, z)), 2 * ip(x, w) * ip(y, z))
    assert torch.allclose(ip(dmul(x, y), z), ip(y, dmul(conj(x), z)))
    assert torch.allclose(ip(dmul(x, y), z), ip(x, dmul(z, conj(y))))
    assert torch.allclose(dmul(conj(x), dmul(x, y)), N(x)[:, None] * y)
    assert torch.allclose(dmul(dmul(y, x), conj(x)), N(x)[:, None] * y)
    assert torch.allclose(dmul(x, x) - 2 * x[:, :1] * x + N(x)[:, None] * one, torch.zeros_like(x), atol=1e-12)
    assert torch.allclose(dmul(x, y) + dmul(y, x), 2 * x[:, :1] * y + 2 * y[:, :1] * x - 2 * ip(x, y)[:, None] * one)
    assert torch.allclose(dmul(conj(x), dmul(w, y)) + dmul(conj(w), dmul(x, y)), 2 * ip(x, w)[:, None] * y)
    assert torch.allclose(dmul(dmul(y, x), conj(w)) + dmul(dmul(y, w), conj(x)), 2 * ip(x, w)[:, None] * y)
    assert torch.allclose(conj(dmul(x, y)), dmul(conj(y), conj(x)))
    assert torch.allclose(dmul(x, y)[:, 0], dmul(y, x)[:, 0])
    assert torch.allclose(dmul(dmul(x, y), z)[:, 0], dmul(x, dmul(y, z))[:, 0])
    assert torch.allclose(dmul(x, dmul(x, y)), dmul(dmul(x, x), y))
    assert torch.allclose(dmul(dmul(y, x), x), dmul(y, dmul(x, x)))
    A = lambda p, q, r: dmul(dmul(p, q), r) - dmul(p, dmul(q, r))
    assert torch.allclose(A(x, y, z), -A(y, x, z), atol=1e-12) and torch.allclose(A(x, y, z), -A(x, z, y), atol=1e-12)
    assert torch.allclose(A(x, y, z), A(y, z, x), atol=1e-12)


def test_lemma_E_15_quaternions_and_octonions():
    """Lemma E.15: M(a,b) = [[a,-b],[b*,a*]] is an injective homomorphism H -> M_2(C) with M(q*) = M(q)^*; H is not
    commutative ((i,0)(0,1) = (0,i), (0,1)(i,0) = (0,-i)); O is not associative ((pq)s = (0,(0,i)), p(qs) = (0,(0,-i)));
    x*/N(x) is a two-sided inverse; CD(O) is not a composition algebra (Hurwitz, Theorem E.17)."""
    g = gen(31)

    def M(q):
        za, zb = torch.complex(q[..., 0], q[..., 1]), torch.complex(q[..., 2], q[..., 3])
        return torch.stack([torch.stack([za, -zb], -1), torch.stack([zb.conj(), za.conj()], -1)], -2)
    p4, q4 = torch.randn(2, 200, 4, dtype=F64, generator=g)
    assert torch.allclose(M(dmul(p4, q4)), M(p4) @ M(q4)) and torch.allclose(M(conj(p4)), M(p4).mH)
    assert err(M(p4).reshape(200, -1).abs().sum(-1), torch.zeros(200)) > 0
    e4, e8 = torch.eye(4, dtype=F64), torch.eye(8, dtype=F64)
    assert torch.equal(dmul(e4[1], e4[2]), e4[3]) and torch.equal(dmul(e4[2], e4[1]), -e4[3])
    p, q, s = e8[1], e8[2], e8[4]
    assert torch.equal(dmul(dmul(p, q), s), e8[7]) and torch.equal(dmul(p, dmul(q, s)), -e8[7])
    for d in (1, 2, 4, 8):
        x = torch.randn(100, d, dtype=F64, generator=g)
        inv = conj(x) / norm2(x)[:, None]
        one = torch.zeros(d, dtype=F64)
        one[0] = 1
        assert torch.allclose(dmul(x, inv), one.expand_as(x)) and torch.allclose(dmul(inv, x), one.expand_as(x))
    x16, y16 = torch.randn(2, 100, 16, dtype=F64, generator=g)
    assert (norm2(cd(x16, y16)) - norm2(x16) * norm2(y16)).abs().max() > 1e-2


def test_lemma_E_16_doubling_inside_octonions():
    """Lemma E.16: for B = H = span(e0..e3) and i = e4 in O: Bi ⊥ B, N(a+bi) = N(a)+N(b) and
    (a+bi)(c+di) = (ac - d*b) + (da + bc*)i."""
    g = gen(32)
    a, b, c, d = (torch.cat([torch.randn(100, 4, dtype=F64, generator=g), torch.zeros(100, 4, dtype=F64)], -1)
                  for _ in range(4))
    i = torch.eye(8, dtype=F64)[4].expand(100, 8)
    bi = dmul(b, i)
    assert (bi[:, :4].abs().max() < 1e-15) and torch.allclose(norm2(a + bi), norm2(a) + norm2(b))
    lhs = dmul(a + bi, c + dmul(d, i))
    rhs = dmul(a, c) - dmul(conj(d), b) + dmul(dmul(d, a) + dmul(b, conj(c)), i)
    assert torch.allclose(lhs, rhs)


def test_lemma_E_18_artin_and_moufang():
    """Lemma E.18: any two octonions lie in an associative subalgebra span{1, p, q, pq}; Moufang identity
    (ux)(yu) = u((xy)u) = (u(xy))u."""
    g = gen(33)
    p, q, u, x, y = torch.randn(5, 300, 8, dtype=F64, generator=g)
    one = torch.zeros(8, dtype=F64)
    one[0] = 1
    basis = [one.expand_as(p), p, q, dmul(p, q)]
    co = torch.randn(3, 4, 300, 1, dtype=F64, generator=g)
    a, b, c = (sum(co[k, l] * basis[l] for l in range(4)) for k in range(3))
    assert torch.allclose(dmul(dmul(a, b), c), dmul(a, dmul(b, c)), atol=1e-9)
    assert torch.allclose(dmul(dmul(u, x), dmul(y, u)), dmul(u, dmul(dmul(x, y), u)))
    assert torch.allclose(dmul(u, dmul(dmul(x, y), u)), dmul(dmul(u, dmul(x, y)), u))
    assert (dmul(dmul(u, x), y) - dmul(u, dmul(x, y))).abs().max() > 1


def test_theorem_E_19_quaternion_table():
    """Theorem E.19 (proof, step (iv)): in H, I = (i,0), J = (0,1), K = (0,i) satisfy I^2 = J^2 = K^2 = -1, IJ = -JI = K,
    JK = -KJ = I, KI = -IK = J."""
    e = torch.eye(4, dtype=F64)
    one, I, J, K = e
    for X in (I, J, K):
        assert torch.equal(dmul(X, X), -one)
    for X, Y, Z in ((I, J, K), (J, K, I), (K, I, J)):
        assert torch.equal(dmul(X, Y), Z) and torch.equal(dmul(Y, X), -Z)


# ======================================================================================================================
# Appendix E: Jordan algebras by structure constants
# ======================================================================================================================
def algebras():
    return [('Herm3C', HermD(3, 2), 3), ('Herm3H', HermD(3, 4), 3), ('Albert', albert(), 3), ('Herm2O', HermD(2, 8), 2),
            ('V4', Spin(4), 2)]


@pytest.mark.parametrize('idx', range(5))
def test_lemma_E_1_and_F_1_linearised_identities(idx):
    """Lemma E.1(1)-(4) and Lemma F.1(1),(2) as operator identities in Herm_3(C), Herm_3(H), the Albert algebra,
    Herm_2(O) and V_4."""
    name, J, _ = algebras()[idx]
    g = gen(34 + idx)
    x, y, z, w = torch.randn(4, J.dim, dtype=F64, generator=g)
    L = lambda v: Lop(J, v)
    br = lambda A, B: A @ B - B @ A
    x2 = mul(J, x, x)
    tol = dict(atol=1e-9 * float(J.dim))
    assert torch.allclose(br(L(x), L(x2)), torch.zeros(J.dim, J.dim, dtype=F64), **tol)
    assert torch.allclose(br(L(x2), L(z)) + 2 * br(L(mul(J, x, z)), L(x)), torch.zeros(J.dim, J.dim, dtype=F64), **tol)
    assert torch.allclose(L(mul(J, x2, y)), 2 * L(mul(J, x, y)) @ L(x) + L(x2) @ L(y) - 2 * L(x) @ L(y) @ L(x), **tol)
    xw = mul(J, x, w)
    rhs3 = (L(mul(J, w, y)) @ L(x) + L(mul(J, x, y)) @ L(w) + L(xw) @ L(y) - L(w) @ L(y) @ L(x) - L(x) @ L(y) @ L(w))
    assert torch.allclose(L(mul(J, xw, y)), rhs3, **tol)
    d4 = L(mul(J, xw, y)) - L(mul(J, x, mul(J, w, y)))
    assert torch.allclose(d4, br(L(w), br(L(x), L(y))), **tol)
    assert abs(torch.trace(d4).item()) < 1e-9 * J.dim
    a, b, d = x, y, z
    f1 = br(L(a), L(mul(J, b, d))) + br(L(b), L(mul(J, d, a))) + br(L(d), L(mul(J, a, b)))
    assert torch.allclose(f1, torch.zeros(J.dim, J.dim, dtype=F64), **tol)


@pytest.mark.parametrize('idx', range(5))
def test_lemma_E_2_power_associativity(idx):
    """Lemma E.2 / Proposition F.2: x^p∘x^q = x^{p+q} and the L(x^k) commute (also in the non-special Albert algebra)."""
    name, J, _ = algebras()[idx]
    x = torch.randn(J.dim, dtype=F64, generator=gen(40 + idx)) / 2
    pw = [jpow(J, x, k) for k in range(7)]
    for p in range(1, 4):
        for q in range(1, 4):
            assert torch.allclose(mul(J, pw[p], pw[q]), pw[p + q], atol=1e-9 * pw[p + q].abs().max())
    for p in range(1, 4):
        for q in range(1, 4):
            A, B = Lop(J, pw[p]), Lop(J, pw[q])
            assert (A @ B - B @ A).abs().max() < 1e-9 * (A.abs().max() * B.abs().max())
    assert torch.allclose(pw[0], J.e)


@pytest.mark.parametrize('idx', range(5))
def test_lemma_E_3_spectral_lemma(idx):
    """Lemma E.3 / Theorem F.3: the minimal polynomial of a generic x (from Jordan powers only) has distinct real roots,
    the Lagrange polynomials give orthogonal idempotents c_i in R[x] with x = sum lambda_i c_i, x^n = sum lambda_i^n c_i,
    and sum c_i = 1 (the frame size is the rank)."""
    name, J, r = algebras()[idx]
    x = torch.randn(J.dim, dtype=F64, generator=gen(45 + idx))
    lam, cs, resid = spectral(J, x, r)
    assert resid < 1e-9 * np.abs(jpow(J, x, r).numpy()).max()                       # minimal polynomial of degree r
    assert np.abs(lam.imag).max() < 1e-8 and len(set(np.round(lam.real, 6))) == r
    for i in range(r):
        for l in range(r):
            assert torch.allclose(mul(J, cs[i], cs[l]), cs[i] if i == l else torch.zeros(J.dim, dtype=F64), atol=1e-8)
    lam = torch.tensor(lam.real, dtype=F64)
    assert torch.allclose(sum(l * c for l, c in zip(lam, cs)), x, atol=1e-9)
    assert torch.allclose(sum(l ** 3 * c for l, c in zip(lam, cs)), jpow(J, x, 3), atol=1e-8)
    assert torch.allclose(cs.sum(0), J.e, atol=1e-9)
    for c in cs:                                                                       # primitive: trace 1 (Corollary E.28)
        assert abs((r / J.dim) * torch.trace(Lop(J, c)).item() - 1) < 1e-8


@pytest.mark.parametrize('idx', range(5))
def test_prop_E_4_and_lemma_F_4_peirce(idx):
    """Proposition E.4 / Lemma F.4: for an idempotent c, L(c)(2L(c)-1)(L(c)-1) = 0, L(c) >= 0 is symmetric for the trace
    form, the Peirce rules hold, and P(2c-1) is a Jordan automorphism (+1 on J_0+J_1, -1 on J_1/2)."""
    name, J, r = algebras()[idx]
    g = gen(50 + idx)
    if isinstance(J, Spin):
        w = torch.nn.functional.normalize(torch.randn(J.dim - 1, dtype=F64, generator=g), dim=0)
        c = 0.5 * torch.cat([torch.ones(1, dtype=F64), w])
    else:
        _, cs, _ = spectral(J, torch.randn(J.dim, dtype=F64, generator=g), r)
        c = cs[0]
    I = torch.eye(J.dim, dtype=F64)
    Lc = Lop(J, c)
    assert torch.allclose(Lc @ (2 * Lc - I) @ (Lc - I), torch.zeros_like(Lc), atol=1e-9)
    T = trace_form(J)
    assert torch.allclose(T @ Lc, (T @ Lc).T, atol=1e-8)                              # L(c) symmetric for tau
    ev = torch.linalg.eigvals(Lc).real
    assert ((ev.unsqueeze(-1) - torch.tensor([0, 0.5, 1.0], dtype=F64)).abs().min(-1).values < 1e-8).all()
    projs = {}
    for nu in (0.0, 0.5, 1.0):
        others = [m for m in (0.0, 0.5, 1.0) if m != nu]
        Pn = (Lc - others[0] * I) @ (Lc - others[1] * I) / ((nu - others[0]) * (nu - others[1]))
        projs[nu] = Pn
    z, v = torch.randn(2, J.dim, dtype=F64, generator=g)
    comp = {nu: (projs[nu] @ z, projs[nu] @ v) for nu in projs}
    inn = lambda nu, u_: torch.allclose(projs[nu] @ u_, u_, atol=1e-8)
    assert inn(1.0, mul(J, comp[1.0][0], comp[1.0][1])) and inn(0.0, mul(J, comp[0.0][0], comp[0.0][1]))
    assert mul(J, comp[1.0][0], comp[0.0][1]).abs().max() < 1e-8
    assert inn(0.5, mul(J, comp[1.0][0] + comp[0.0][0], comp[0.5][1]))
    h = mul(J, comp[0.5][0], comp[0.5][1])
    assert torch.allclose(projs[0.5] @ h, torch.zeros_like(h), atol=1e-8)
    s = 2 * c - J.e
    Ps = Pop(J, s)
    assert torch.allclose(Ps, projs[0.0] + projs[1.0] - projs[0.5], atol=1e-8)
    assert torch.allclose(Ps @ mul(J, z, v), mul(J, Ps @ z, Ps @ v), atol=1e-8)


@pytest.mark.parametrize('idx', range(5))
def test_theorem_E_5_trace_form(idx):
    """Theorem E.5: tau(x,y) = Tr L(x∘y) is associative and positive definite, Tr L(c) >= 1 for idempotents, the unit e
    satisfies L(e) = id; Euclidean => formally real."""
    name, J, r = algebras()[idx]
    g = gen(55 + idx)
    T = trace_form(J)
    x, y, z = torch.randn(3, J.dim, dtype=F64, generator=g)
    assert torch.allclose(T, T.T) and torch.linalg.eigvalsh(T).min() > 0.5
    assert torch.allclose(mul(J, x, y) @ T @ z, x @ T @ mul(J, y, z))
    assert torch.allclose(Lop(J, J.e), torch.eye(J.dim, dtype=F64))
    xs = torch.randn(4, J.dim, dtype=F64, generator=g)
    s = sum(mul(J, v, v) for v in xs)
    assert abs((s @ J.G @ J.e).item() - sum((v @ J.G @ v).item() for v in xs)) < 1e-9


@pytest.mark.parametrize('n,d', [(2, 1), (3, 1), (4, 1), (3, 2), (4, 2), (3, 4), (4, 4), (2, 8), (3, 8)])
def test_theorem_E_8_prop_E_11_peirce_frames(n, d):
    """Theorem E.8 / Proposition E.11 / Theorem E.24(3): (E_ii) is a Jordan frame of Herm_n(D), the Peirce spaces J_ij =
    {a[ij]} have dimension d = dim D, x^2 = N(x)(c_i + c_j) on J_ij, t_i = Tr L(c_i) = 1 + (n-1)d/2 and dim = n t."""
    J = HermD(n, d)
    I = torch.eye(J.dim, dtype=F64)
    Ls = [Lop(J, J.diag(i)) for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            rank = torch.linalg.matrix_rank(torch.cat([Ls[i] - I / 2, Ls[j] - I / 2])).item()
            assert J.dim - rank == d
    for i in range(n):
        assert abs(torch.trace(Ls[i]).item() - (1 + (n - 1) * d / 2)) < 1e-12
        assert torch.allclose(mul(J, J.diag(i), J.diag(i)), J.diag(i))
    assert J.dim == n * (1 + (n - 1) * d / 2)
    g = gen(60)
    T = trace_form(J)
    for i, j in J.pairs[:3]:
        a, b = torch.randn(2, d, dtype=F64, generator=g)
        X, Y = J.entry(a, i, j), J.entry(b, i, j)
        assert torch.allclose(mul(J, X, X), norm2(a) * (J.diag(i) + J.diag(j)))
        assert torch.allclose(mul(J, X, Y), (a @ b) * (J.diag(i) + J.diag(j)))
        t = 1 + (n - 1) * d / 2
        assert abs((X @ T @ Y).item() - 2 * t * (a @ b).item()) < 1e-10


def test_theorem_E_9_prop_E_11_direct_sum_is_not_simple():
    """Theorem E.9 / Proposition E.11(1): for the direct sum Herm_2(C) ⊕ V_2 with the frame (E_11, E_22, c_w, c_-w), the
    Peirce graph is disconnected (J_ij = 0 across summands), and the summands are orthogonal ideals for tau."""
    A, B = HermD(2, 2), Spin(2)
    dim = A.dim + B.dim
    C = torch.zeros(dim, dim, dim, dtype=F64)
    C[:A.dim, :A.dim, :A.dim] = A.C
    C[A.dim:, A.dim:, A.dim:] = B.C

    class S:
        pass
    J = S()
    J.C, J.dim = C, dim
    J.e = torch.cat([A.e, B.e])
    frame = [torch.cat([A.diag(0), torch.zeros(3, dtype=F64)]), torch.cat([A.diag(1), torch.zeros(3, dtype=F64)]),
             torch.cat([torch.zeros(4, dtype=F64), torch.tensor([0.5, 0.5, 0.0], dtype=F64)]),
             torch.cat([torch.zeros(4, dtype=F64), torch.tensor([0.5, -0.5, 0.0], dtype=F64)])]
    I = torch.eye(dim, dtype=F64)
    Ls = [Lop(J, f) for f in frame]
    dims = {(i, j): dim - torch.linalg.matrix_rank(torch.cat([Ls[i] - I / 2, Ls[j] - I / 2])).item()
            for i in range(4) for j in range(i + 1, 4)}
    assert dims == {(0, 1): 2, (2, 3): 1, (0, 2): 0, (0, 3): 0, (1, 2): 0, (1, 3): 0}
    T = trace_form(J)
    assert T[:A.dim, A.dim:].abs().max() == 0


def test_lemma_E_10_E_12_and_coordinatisation():
    """Lemma E.10 / E.12 and Theorem E.22, step 1, in the Albert algebra: u∘(u∘w) = N(u)w/4, N(m(x,y)) = N(x)N(y),
    m(x, m(x,y)) = N(x)y, <m(x,y),z> = <y,m(x,z)>, and p·q = m(m(p,e_23), m(e_12,q)) on J_13 is the octonion product with
    unit 1' = m(e_12, e_23); in Herm_4(H), (u∘w)∘z = u∘(w∘z) for u, w, z in J_12, J_23, J_34."""
    J = albert()
    g = gen(61)
    a, b, c = octonions(3, 61)
    m = lambda X, Y: 2 * mul(J, X, Y)
    u, w = J.entry(a, 0, 1), J.entry(b, 1, 2)
    Nu = norm2(a)
    assert torch.allclose(mul(J, u, mul(J, u, w)), Nu / 4 * w)
    mxy = m(u, w)
    assert torch.allclose(mxy, J.entry(dmul(a, b), 0, 2))
    assert torch.allclose(norm2(mxy[J.n + 8:J.n + 16]), Nu * norm2(b))
    assert torch.allclose(m(u, m(u, w)), Nu * w)
    z = J.entry(c, 0, 2)
    ip = lambda X, Y: (X @ J.G @ Y) / 2                                               # <a[ij], b[ij]> = <a, b>
    assert abs(ip(m(u, w), z).item() - ip(w, m(u, z)).item()) < 1e-10
    e12, e23 = J.entry(torch.eye(8, dtype=F64)[0], 0, 1), J.entry(torch.eye(8, dtype=F64)[0], 1, 2)
    p, q = J.entry(a, 0, 2), J.entry(b, 0, 2)
    pq = m(m(p, e23), m(e12, q))
    assert torch.allclose(pq, J.entry(dmul(a, b), 0, 2))
    one = m(e12, e23)
    assert torch.allclose(one, J.entry(torch.eye(8, dtype=F64)[0], 0, 2))
    assert torch.allclose(m(m(p, e23), m(e12, one)), p) and torch.allclose(m(m(one, e23), m(e12, q)), q)
    H4 = HermD(4, 4)
    x4, y4, z4 = torch.randn(3, 4, dtype=F64, generator=g)
    U, W, Z = H4.entry(x4, 0, 1), H4.entry(y4, 1, 2), H4.entry(z4, 2, 3)
    assert torch.allclose(mul(H4, mul(H4, U, W), Z), mul(H4, U, mul(H4, W, Z)))
    assert torch.allclose(mul(H4, mul(H4, U, W), Z), H4.entry(dmul(dmul(x4, y4), z4), 0, 3) / 4)


def test_lemma_E_20_table_and_inner_product():
    """Lemma E.20 in Herm_3(O): E_ii∘E_jj = delta_ij E_ii, E_kk∘a[ij] = (delta_ki + delta_kj)a[ij]/2,
    a[ij]∘b[ij] = <a,b>(E_ii+E_jj), a[ij]∘b[jk] = (ab)[ik]/2, and (X|Y) = Re tr(XY) is associative and positive."""
    J = albert()
    g = gen(62)
    a, b = octonions(2, 62)
    for i in range(3):
        for j in range(3):
            assert torch.allclose(mul(J, J.diag(i), J.diag(j)), (i == j) * J.diag(i))
    for (i, j), k in (((0, 1), 2), ((0, 2), 1), ((1, 2), 0)):
        A = J.entry(a, i, j)
        for kk in range(3):
            assert torch.allclose(mul(J, J.diag(kk), A), 0.5 * ((kk == i) + (kk == j)) * A)
        assert torch.allclose(mul(J, A, J.entry(b, i, j)), (a @ b) * (J.diag(i) + J.diag(j)))
    for i, j, k in ((0, 1, 2), (1, 0, 2), (0, 2, 1), (2, 0, 1), (1, 2, 0), (2, 1, 0)):
        assert torch.allclose(mul(J, J.entry(a, i, j), J.entry(b, j, k)), 0.5 * J.entry(dmul(a, b), i, k))
    assert torch.allclose(J.entry(a, 1, 0), J.entry(conj(a), 0, 1))
    X, Y, Z = torch.randn(3, 27, dtype=F64, generator=g)
    G = J.G
    assert torch.allclose(mul(J, X, Y) @ G @ Z, Y @ G @ mul(J, X, Z))
    ReTr = lambda v: v[..., :3].sum(-1)
    assert torch.allclose(ReTr(mul(J, X, Y)), X @ G @ Y)
    assert torch.linalg.eigvalsh(G).min() > 0


def test_lemma_E_21_orientations():
    """Lemma E.21 (ii): in Herm_3(O) with iota_ij(a) = a[ij] and iota_ji(a) = iota_ij(a*), m(iota_ij a, iota_jk b) =
    iota_ik(ab) for all six ordered triples."""
    J = albert()
    a, b = octonions(2, 63)
    iota = lambda x, i, j: J.entry(x, i, j) if i < j else J.entry(conj(x), j, i)
    for i, j, k in ((0, 1, 2), (1, 0, 2), (0, 2, 1), (2, 0, 1), (1, 2, 0), (2, 1, 0)):
        assert torch.allclose(2 * mul(J, iota(a, i, j), iota(b, j, k)), iota(dmul(a, b), i, k))


def moufang_auto(v, u):
    """The automorphism (L_u, R_u, B_u) of Lemma E.23 on Herm_3(O) coordinates (pairs (0,1),(0,2),(1,2))."""
    a12, a13, a23 = v[..., 3:11], v[..., 11:19], v[..., 19:27]
    u = u.expand_as(a12)
    return torch.cat([v[..., :3], dmul(u, a12), dmul(dmul(u, a13), u), dmul(a23, u)], -1)


def test_theorem_E_24_2_albert_jordan_identity_and_lemma_E_23():
    """Theorem E.24(2): the Albert algebra satisfies the Jordan identity; Lemma E.23: (L_u, R_u, B_u) is an automorphism for
    N(u) = 1 (and u(ab)u = (ua)(bu)), and steps (a)-(c) of the proof make the (1,3) and (2,3) entries real."""
    J = albert()
    g = gen(64)
    X, Y = torch.randn(2, 100, 27, dtype=F64, generator=g)
    X2 = mul(J, X, X)
    assert (mul(J, mul(J, X2, Y), X) - mul(J, X2, mul(J, Y, X))).abs().max() < 1e-10
    u = torch.nn.functional.normalize(torch.randn(8, dtype=F64, generator=g), dim=0)
    a, b = octonions(2, 65)
    assert torch.allclose(dmul(dmul(u, dmul(a, b)), u), dmul(dmul(u, a), dmul(b, u)))
    assert torch.allclose(moufang_auto(mul(J, X, Y), u), mul(J, moufang_auto(X, u), moufang_auto(Y, u)))
    e0 = torch.eye(8, dtype=F64)[0]
    for x in torch.randn(10, 27, dtype=F64, generator=g):
        s = x[11:19]
        im = torch.cat([torch.zeros(1, dtype=F64), s[1:]])
        v, th = im / im.norm(), math.atan2(im.norm().item(), s[0].item())
        uu = math.cos(th / 2) * e0 - math.sin(th / 2) * v
        assert torch.allclose(dmul(dmul(uu, s), uu), s.norm() * e0)
        x1 = moufang_auto(x, uu)
        q = x1[19:27]
        gg = torch.cat([torch.zeros(1, dtype=F64), torch.randn(7, dtype=F64, generator=g)])
        gg = gg - (gg @ q) / (q[1:] @ q[1:]) * torch.cat([torch.zeros(1, dtype=F64), q[1:]])
        x2 = moufang_auto(x1, gg / gg.norm())
        assert abs(x2[19].item()) < 1e-12 and x2[12:19].abs().max() < 1e-12
        x3 = moufang_auto(x2, x2[19:27] / x2[19:27].norm())
        assert x3[12:19].abs().max() < 1e-12 and x3[20:27].abs().max() < 1e-12
        assert abs(abs(x3[11].item()) - s.norm().item()) < 1e-12
        assert abs(x3[19].item() + x2[19:27].norm().item()) < 1e-12                    # qu = q^2/|q| = -|q|
        assert torch.allclose(x3[:11], torch.cat([x[:3], x3[3:11]]))


def test_theorem_E_24_1_quaternionic_embedding():
    """Theorem E.24(1): replacing quaternion entries by 2x2 blocks M(q) embeds Herm_n(H) into Herm_2n(C) as a Jordan
    subalgebra (dim n(2n-1)), which is not the self-adjoint part of the *-algebra M_2n(C) it generates."""

    def M(q):
        za, zb = torch.complex(q[..., 0], q[..., 1]), torch.complex(q[..., 2], q[..., 3])
        return torch.stack([torch.stack([za, -zb], -1), torch.stack([zb.conj(), za.conj()], -1)], -2)
    for n in (2, 3):
        J = HermD(n, 4)
        emb = lambda v: M(J.mat(v)).permute(*range(v.dim() - 1), -4, -2, -3, -1).reshape(*v.shape[:-1], 2 * n, 2 * n)
        g = gen(66 + n)
        x, y = torch.randn(2, 50, J.dim, dtype=F64, generator=g)
        X, Y = emb(x), emb(y)
        assert err(X, X.mH) < 1e-14 and err(emb(mul(J, x, y)), jm(X, Y)) < 1e-12
        basis = emb(torch.eye(J.dim, dtype=F64))
        assert jd.is_jordan_subalgebra(basis)
        r = jd.compare_with_star_algebra(basis)
        assert (r['jordan_dim'], r['star_dim']) == (n * (2 * n - 1), 4 * n * n)


def test_corollary_E_25_herm4_octonions_not_jordan():
    """Corollary E.25: in Herm_4(O), (a[12]∘b[23])∘c[34] - a[12]∘(b[23]∘c[34]) = ((ab)c - a(bc))[14]/4 != 0, and the Jordan
    identity fails."""
    J = HermD(4, 8)
    a, b, c = octonions(3, 70)
    A, B, C = J.entry(a, 0, 1), J.entry(b, 1, 2), J.entry(c, 2, 3)
    diff = mul(J, mul(J, A, B), C) - mul(J, A, mul(J, B, C))
    assoc = dmul(dmul(a, b), c) - dmul(a, dmul(b, c))
    assert torch.allclose(diff, J.entry(assoc, 0, 3) / 4) and assoc.abs().max() > 0.1
    X, Y = torch.randn(2, 50, J.dim, dtype=F64, generator=gen(71))
    X2 = mul(J, X, X)
    assert (mul(J, mul(J, X2, Y), X) - mul(J, X2, mul(J, Y, X))).abs().max() > 1


def test_theorem_E_26_mechanism_in_associative_embedding():
    """Theorem E.26 (proof) for Herm_3(H) ⊂ M_6(C): e_i e_j = 0 (i != j), x = sum e_p x e_q with only the (i,j),(j,i)
    blocks for x = phi(a[ij]), alpha_ij(a) alpha_jk(b) = alpha_ik(ab), and alpha_13(ay) = alpha_13(a) pi alpha_13(y)."""
    J = HermD(3, 4)

    def M(q):
        za, zb = torch.complex(q[..., 0], q[..., 1]), torch.complex(q[..., 2], q[..., 3])
        return torch.stack([torch.stack([za, -zb], -1), torch.stack([zb.conj(), za.conj()], -1)], -2)
    phi = lambda v: M(J.mat(v)).permute(-4, -2, -3, -1).reshape(6, 6)
    e = [phi(J.diag(i)) for i in range(3)]
    for i in range(3):
        assert err(e[i] @ e[i], e[i]) == 0
        for j in range(3):
            if i != j:
                assert err(e[i] @ e[j], torch.zeros(6, 6)) == 0
    g = gen(72)
    a, b, y = torch.randn(3, 4, dtype=F64, generator=g)
    iota = lambda x, i, j: J.entry(x, i, j) if i < j else J.entry(conj(x), j, i)
    alpha = lambda x, i, j: e[i] @ phi(iota(x, i, j)) @ e[j]
    X = phi(iota(a, 0, 1))
    assert err(X, alpha(a, 0, 1) + alpha(conj(a), 1, 0)) < 1e-14
    for i, j, k in ((0, 1, 2), (2, 0, 1), (1, 2, 0)):
        assert err(alpha(a, i, j) @ alpha(b, j, k), alpha(dmul(a, b), i, k)) < 1e-12
    pi = alpha(torch.eye(4, dtype=F64)[0], 2, 0)
    assert err(alpha(dmul(a, y), 0, 2), alpha(a, 0, 2) @ pi @ alpha(y, 0, 2)) < 1e-12


def test_theorem_E_27_separating_invariants():
    """Theorem E.27(2),(3): the pairs (s, dim) of R, V_m, Herm_n(R/C/H) (n >= 3) and the Albert algebra are pairwise
    different, and 27 = r + d r(r-1)/2 has no solution with r >= 4, d in {1, 2, 4}."""
    pairs = [(1, 1)] + [(2, m + 1) for m in range(2, 200)] + [(3, 27)]
    for n in range(3, 40):
        pairs += [(n, n * (n + 1) // 2), (n, n * n), (n, n * (2 * n - 1))]
    assert len(pairs) == len(set(pairs))
    sols = [(r, d) for r in range(4, 30) for d in (1, 2, 4) if r + d * r * (r - 1) // 2 == 27]
    assert sols == []
    assert {(r, d) for r in (1, 2, 3) for d in (1, 2, 4, 8) if r + d * r * (r - 1) // 2 == 27} == {(3, 8)}


@pytest.mark.parametrize('idx', range(5))
def test_corollary_E_28_trace(idx):
    """Corollary E.28: tr x = (r/dim J) Tr L(x) gives tr c = 1 on primitive idempotents, tr x = sum of eigenvalues, and
    tr(x∘y) is an associative inner product (for Herm_n(D) it is Re tr, for V_m it is 2a)."""
    name, J, r = algebras()[idx]
    g = gen(73 + idx)
    tr = lambda v: (r / J.dim) * torch.trace(Lop(J, v))
    x, y, z = torch.randn(3, J.dim, dtype=F64, generator=g)
    lam, cs, _ = spectral(J, x, r)
    assert abs(tr(x).item() - lam.real.sum()) < 1e-8
    if isinstance(J, Spin):
        assert abs(tr(x).item() - 2 * x[0].item()) < 1e-12
    else:
        assert abs(tr(x).item() - x[:J.n].sum().item()) < 1e-12
    assert abs(tr(mul(J, mul(J, x, y), z)).item() - tr(mul(J, x, mul(J, y, z))).item()) < 1e-10
    Ib = torch.eye(J.dim, dtype=F64)
    Gm = torch.stack([torch.stack([tr(mul(J, Ib[a], Ib[b])) for b in range(J.dim)]) for a in range(J.dim)])
    assert torch.linalg.eigvalsh(Gm).min() > 0


# ======================================================================================================================
# Appendix F
# ======================================================================================================================
@pytest.mark.parametrize('n', [1, 3, 5])
def test_theorem_F_3_spectral_theorem_in_M_sa(n):
    """Theorem F.3 in Herm_n(C): the Lagrange idempotents l_i(x) are the eigenprojections, c_i∘c_j = delta_ij c_i,
    sum c_i = 1, and p(x) = sum p(lambda_i)c_i; the idempotents of (M_sa, ∘) are the projections."""
    g = gen(80 + n)
    U = torch.linalg.qr(torch.randn(n, n, dtype=C128, generator=g))[0]
    lam = torch.tensor([2.0, -1.0, 0.5, 3.0, -2.5][:n], dtype=F64)
    x = U @ torch.diag(lam.to(C128)) @ U.mH
    distinct = lam.tolist()
    for i, li in enumerate(distinct):
        ell = torch.eye(n, dtype=C128)
        for j, lj in enumerate(distinct):
            if j != i:
                ell = ell @ (x - lj * torch.eye(n, dtype=C128)) / (li - lj)
        proj = U[:, i:i + 1] @ U[:, i:i + 1].mH
        assert err(ell, proj) < 1e-12 and err(jm(ell, ell), ell) < 1e-12
    p3 = x @ x @ x - 2 * x
    assert err(p3, U @ torch.diag((lam ** 3 - 2 * lam).to(C128)) @ U.mH) < 1e-12


@pytest.mark.parametrize('n', [2, 3])
def test_lemma_F_4_4_scalars_on_peirce_spaces(n):
    """Lemma F.4(4) in Herm_n(C): for a = sum a_l c_l, L(a) acts on V_ij as (a_i + a_j)/2 and P(a) as a_i a_j."""
    g = gen(85)
    U = torch.linalg.qr(torch.randn(n, n, dtype=C128, generator=g))[0]
    av = torch.randn(n, dtype=F64, generator=g)
    a = U @ torch.diag(av.to(C128)) @ U.mH
    for i in range(n):
        for j in range(i, n):
            Eij = U[:, i:i + 1] @ U[:, j:j + 1].mH
            z = Eij + Eij.mH                                                            # in V_ij (or V_ii)
            assert err(jm(a, z), 0.5 * (av[i] + av[j]) * z) < 1e-12
            assert err(a @ z @ a, av[i] * av[j] * z) < 1e-12


@pytest.mark.parametrize('kind', ['herm', 'spin'])
def test_theorem_F_5_self_duality(kind):
    """Theorem F.5: y in Q iff all lambda_i >= 0 iff L(y) >= 0 iff <y,x> >= 0 on Q; int Q = {L(y) > 0}; <x,y> = 0 for
    x, y in Q implies x∘y = 0; Q ∩ -Q = 0 (x^2 = 0 => x = 0)."""
    g = gen(86)
    if kind == 'herm':
        n = 3
        coords = jd.jordan_algebra_generated(torch.stack([herm_t(n, g) for _ in range(3)]))
        Bs = coords.basis                                                               # orthonormal basis of Herm_3(C)
        assert coords.dim == 9
        for _ in range(20):
            y = herm_t(n, g)
            Ly = torch.stack([torch.stack([torch.trace(Bs[a] @ jm(y, Bs[b])).real for b in range(9)]) for a in range(9)])
            assert torch.allclose(Ly, Ly.T)
            wL, wy = torch.linalg.eigvalsh(Ly), torch.linalg.eigvalsh(y)
            assert abs(wL[0].item() - wy[0].item()) < 1e-10                              # spec L(y) = (l_i + l_j)/2
            if wy[0] < 0:
                v = torch.linalg.eigh(y)[1][:, :1]
                assert torch.trace(y @ (v @ v.mH)).real < 0                              # witness c_i in Q
        U = torch.linalg.qr(torch.randn(n, n, dtype=C128, generator=g))[0]
        x = U @ torch.diag(torch.tensor([1.0, 2.0, 0], dtype=C128)) @ U.mH
        y = U @ torch.diag(torch.tensor([0, 0, 3.0], dtype=C128)) @ U.mH
        assert abs(torch.trace(x @ y).item()) < 1e-12 and err(jm(x, y), torch.zeros(n, n)) < 1e-12
    else:
        V = SpinFactor(4)
        I5 = torch.eye(5, dtype=F64)
        Lmat = lambda b: V.product(b.expand(5, 5), I5).T
        for _ in range(50):
            y = V.random(generator=g)
            wL = torch.linalg.eigvalsh(Lmat(y))
            assert abs(wL[0].item() - V.eigenvalues(y)[0].item()) < 1e-12               # spec L(y) = {l-, a, l+}
            assert (wL[0] >= 0) == V.in_cone(y, tol=0.0)
        w = torch.nn.functional.normalize(torch.randn(4, dtype=F64, generator=g), dim=0)
        x, y = V.element(2.0, 2.0 * w), V.element(3.0, -3.0 * w)
        assert abs(V.inner(x, y).item()) < 1e-12 and err(V.product(x, y), torch.zeros(5)) < 1e-12
        x0 = V.random(100, generator=g)
        assert (V.inner(V.square(x0), V.unit()) - V.inner(x0, x0)).abs().max() < 1e-12


def test_theorem_F_6_homogeneity_in_M_sa():
    """Theorem F.6 in Herm_n(C): e^{L(b)}1 = exp b (e^{L(b)}x = e^{b/2}xe^{b/2}), P(y) = e^{2L(log y)} for y > 0,
    P(y) = P(s)P(|y|) with s = 2c - 1, and e^{L(b)}, P(y) preserve the positive cone."""
    n = 3
    g = gen(87)
    b, z = herm_t(n, g), herm_t(n, g)
    w, U = torch.linalg.eigh(b)
    eh = U @ torch.diag(torch.exp(w / 2).to(C128)) @ U.mH
    Lmat = torch.einsum('ij,kl->ikjl', 0.5 * b, torch.eye(n, dtype=C128)).reshape(n * n, n * n) + \
        torch.einsum('ij,kl->ikjl', torch.eye(n, dtype=C128), 0.5 * b.T).reshape(n * n, n * n)
    E = torch.linalg.matrix_exp(Lmat)
    assert err((E @ z.reshape(-1)).reshape(n, n), eh @ z @ eh) < 1e-10
    assert err((E @ torch.eye(n, dtype=C128).reshape(-1)).reshape(n, n), torch.linalg.matrix_exp(b)) < 1e-10
    G = torch.randn(n, n, dtype=C128, generator=g)
    y = G @ G.mH + 0.1 * torch.eye(n, dtype=C128)
    wy, Uy = torch.linalg.eigh(y)
    logy = Uy @ torch.diag(torch.log(wy).to(C128)) @ Uy.mH
    L2 = torch.einsum('ij,kl->ikjl', logy, torch.eye(n, dtype=C128)).reshape(n * n, n * n) + \
        torch.einsum('ij,kl->ikjl', torch.eye(n, dtype=C128), logy.T).reshape(n * n, n * n)
    assert err((torch.linalg.matrix_exp(L2) @ z.reshape(-1)).reshape(n, n), y @ z @ y) < 1e-9 * y.abs().max() ** 2
    yi = herm_t(n, g)
    wi, Ui = torch.linalg.eigh(yi)
    s = Ui @ torch.diag(torch.where(wi >= 0, 1.0, -1.0).to(C128)) @ Ui.mH
    ay = Ui @ torch.diag(wi.abs().to(C128)) @ Ui.mH
    assert err(yi @ z @ yi, s @ (ay @ z @ ay) @ s) < 1e-10
    assert err(s @ jm(z, b) @ s, jm(s @ z @ s, s @ b @ s)) < 1e-10                      # P(s) Jordan automorphism


@pytest.mark.parametrize('m', [2, 4])
def test_theorem_F_6_in_spin_factor(m):
    """Theorem F.6 in V_m: e^{L(b)}e = exp b, e^{L(b)}Q ⊆ Q, P(y) = e^{2L(log y)} for y in int Q, P(y) = P(s)P(|y|),
    det P(y) = det(y)^{m+1} and det(P(y)x) = det(y)^2 det(x)."""
    V = SpinFactor(m)
    g = gen(88)
    N = m + 1
    I = torch.eye(N, dtype=F64)
    Lmat = lambda b: V.product(b.expand(N, N), I).T
    Pmat = lambda b: V.quadratic_representation(b.expand(N, N), I).T
    q = V.random(200, generator=g, cone=True)
    for _ in range(10):
        b, yi, yg, x = V.random(generator=g), V.random(generator=g, cone=True), V.random(generator=g), V.random(generator=g)
        E = torch.linalg.matrix_exp(Lmat(b))
        assert torch.allclose(E @ V.unit(), V.exp(b))
        assert V.in_cone((E @ q.T).T, tol=1e-10).all()
        assert torch.allclose(Pmat(yi), torch.linalg.matrix_exp(2 * Lmat(V.log(yi))), atol=1e-9 * Pmat(yi).abs().max())
        lam, idem = V.spectral_decomposition(yg)
        s = (torch.where(lam >= 0, 1.0, -1.0)[:, None] * idem).sum(0)
        assert torch.allclose(Pmat(yg), Pmat(s) @ Pmat(V.apply_function(yg, torch.abs)))
        assert torch.allclose(Pmat(s) @ V.product(x, b), V.product(Pmat(s) @ x, Pmat(s) @ b))
        assert torch.allclose(torch.linalg.det(Pmat(yg)), V.determinant(yg) ** (m + 1))
        assert torch.allclose(V.determinant(V.quadratic_representation(yg, x)), V.determinant(yg) ** 2 * V.determinant(x))
        assert V.in_cone((Pmat(yg) @ q.T).T, tol=1e-10).all()


def test_theorem_F_5_F_6_corollary_F_7_in_albert_algebra():
    """Theorems F.5, F.6(4) and Corollary F.7 in the (non-special) Albert algebra: L(x^2) >= 0, L(y) has a negative
    eigenvalue iff y has one, P(y)x^2 has nonnegative eigenvalues, and tau(|P(y)x|) <= tau(P(y)|x|) with tau = Re tr."""
    J = albert()
    g = gen(89)
    Gh = torch.diag(J.G.diagonal().sqrt())
    sym = lambda A: Gh @ A @ torch.linalg.inv(Gh)                                     # L in orthonormal coordinates
    for _ in range(6):
        x, y = torch.randn(2, 27, dtype=F64, generator=g)
        x2 = mul(J, x, x)
        Lx2 = sym(Lop(J, x2))
        assert torch.allclose(Lx2, Lx2.T, atol=1e-10) and torch.linalg.eigvalsh(Lx2)[0] > -1e-9
        lam, _, _ = spectral(J, y, 3)
        wL = torch.linalg.eigvalsh(sym(Lop(J, y)))
        assert abs(wL[0].item() - lam.real.min()) < 1e-8 and abs(wL[-1].item() - lam.real.max()) < 1e-8
        Px2 = Pop(J, y) @ x2
        lamP, _, _ = spectral(J, Px2, 3)
        assert lamP.real.min() > -1e-8 * np.abs(lamP).max()
        lx, cx, _ = spectral(J, x, 3)
        absx = sum(abs(float(l.real)) * c for l, c in zip(lx, cx))
        z = Pop(J, y) @ x
        lz, cz, _ = spectral(J, z, 3)
        absz = sum(abs(float(l.real)) * c for l, c in zip(lz, cz))
        tau = lambda v: v[:3].sum().item()
        assert abs(tau(absz) - np.abs(lz.real).sum()) < 1e-8 * np.abs(lz).sum()
        assert tau(absz) <= tau(Pop(J, y) @ absx) * (1 + 1e-10)
        xp = x2                                                                          # x in Q: equality
        assert abs(tau(Pop(J, y) @ xp) - np.abs(spectral(J, Pop(J, y) @ xp, 3)[0].real).sum()) < 1e-8 * tau(Pop(J, y) @ xp)


@pytest.mark.parametrize('m', [2, 5])
def test_lemma_F_8_self_dual_cones(m):
    """Lemma F.8 for the Lorentz cone with the standard inner product: <x,y> >= delta ||y|| on the closed cone with
    delta = (a - |u|)/sqrt(2) the distance of x to the boundary; a boundary point x has u in the closed cone, u != 0,
    with <x,u> = 0, and y + tu stays in the cone; g = P(b)O is in G(Omega) and so is g^T."""
    V = SpinFactor(m)
    g = gen(180 + m)
    x, y = V.random(300, generator=g, cone=True), V.random(300, generator=g, cone=True)
    a, u = V.split(x)
    delta = (a - u.norm(dim=-1)) / math.sqrt(2)
    ip = (x[:, None] * y[None]).sum(-1)
    assert (ip >= delta[:, None] * y.norm(dim=-1)[None] * (1 - 1e-12)).all()
    w = torch.nn.functional.normalize(torch.randn(m, dtype=F64, generator=g), dim=0)
    xb, ub = V.element(1.0, w), V.element(1.0, -w)                                     # boundary point, witness
    assert abs((xb * ub).sum().item()) < 1e-15 and V.in_cone(ub)
    assert V.in_cone(y + 3.0 * ub, interior=True).all()
    O = torch.linalg.qr(torch.randn(m, m, dtype=F64, generator=g))[0]
    Om = torch.block_diag(torch.ones(1, 1, dtype=F64), O)
    b = V.random(generator=g)
    Pb = V.quadratic_representation(b.expand(m + 1, m + 1), torch.eye(m + 1, dtype=F64)).T
    for G in (Pb @ Om, (Pb @ Om).T):
        assert V.in_cone((G @ y.T).T, tol=1e-10).all()
        inv = torch.linalg.inv(G)
        assert V.in_cone((inv @ y.T).T, tol=1e-10).all()


def lorentz_phi(a, u_norm):
    """Characteristic function of the Lorentz cone in R^3 by quadrature: (2 pi / c) int_0^inf b e^{-ab} I_1(cb) db."""
    c = u_norm
    f = lambda b: b * ss.ive(1, c * b) * math.exp(-(a - c) * b)
    val, _ = si.quad(f, 0, np.inf, epsabs=0, epsrel=1e-11, limit=200)
    return 2 * math.pi / c * val


def test_lemma_F_9_characteristic_function_of_lorentz_cone():
    """Lemma F.9(2),(5) for the Lorentz cone of V_2 (self-dual for the standard inner product): phi(x) =
    2 pi det(x)^{-3/2} by quadrature, phi(tx) = t^{-3}phi(x), phi(gx) = |det g|^{-1}phi(x) for g = P(y), and phi blows up
    at the boundary."""
    V = SpinFactor(2)
    for a, c in ((1.0, 0.3), (2.0, 1.5), (1.3, 0.05)):
        assert abs(lorentz_phi(a, c) / (2 * math.pi * (a * a - c * c) ** -1.5) - 1) < 1e-8
    x = V.element(1.7, torch.tensor([0.4, -0.9], dtype=F64))
    phi = lambda v: lorentz_phi(v[0].item(), v[1:].norm().item())
    assert abs(phi(2.5 * x) / (2.5 ** -3 * phi(x)) - 1) < 1e-8
    y = V.element(1.1, torch.tensor([0.3, 0.5], dtype=F64))
    Pm = V.quadratic_representation(y.expand(3, 3), torch.eye(3, dtype=F64)).T
    gx = Pm @ x
    assert abs(phi(gx) / (phi(x) / abs(torch.linalg.det(Pm).item())) - 1) < 1e-8
    vals = [phi(V.element(1.0, torch.tensor([1 - eps, 0.0], dtype=F64))) for eps in (1e-1, 1e-2, 1e-3)]
    assert vals[0] < vals[1] < vals[2] and vals[2] > 1e3


def _cone_setup(kind, n):
    """(f, e, N, mat->vec helpers) for f = log phi up to a constant: Herm_n(R) (Example F.13) or the Lorentz cone of V_n."""
    if kind == 'sdp':
        pairs = [(i, j) for i in range(n) for j in range(i, n)]
        Bs = torch.zeros(len(pairs), n, n, dtype=F64)
        for k, (i, j) in enumerate(pairs):
            Bs[k, i, j] = Bs[k, j, i] = 1.0 if i == j else 1 / math.sqrt(2)
        mat_ = lambda v: torch.einsum('k,kij->ij', v, Bs)
        vec_ = lambda A: torch.einsum('kij,ij->k', Bs, A)
        f = lambda v: -(n + 1) / 2 * torch.logdet(mat_(v))
        kappa = math.sqrt((n + 1) / 2)
        e = vec_(kappa * torch.eye(n, dtype=F64))
        prod = lambda u, v: vec_(mat_(u) @ mat_(v) + mat_(v) @ mat_(u)) / (2 * kappa)
        return f, e, len(pairs), prod, mat_, vec_
    V = SpinFactor(n)
    f = lambda v: -(n + 1) / 2 * torch.log(v[0] ** 2 - (v[1:] ** 2).sum())
    kappa = math.sqrt(n + 1)
    e = kappa * V.unit()
    prod = lambda u, v: V.product(u, v) / kappa
    return f, e, n + 1, prod, None, None


@pytest.mark.parametrize('kind,n', [('sdp', 2), ('sdp', 3), ('lorentz', 2), ('lorentz', 4)])
def test_prop_F_11_theorem_F_12_koecher_vinberg(kind, n):
    """Proposition F.11 and Theorem F.12 (Example F.13): e* = e, H(e) = I, x** = x, H(x*) = H(x)^{-1}, H(x)x = x*,
    D^3f(x)[x,u,v] = -2<H(x)u,v>; the product -D^3f(e)/2 has unit e, equals the scaled Jordan product, satisfies the
    Jordan identity, and the inner derivations [L(x),L(y)] are antisymmetric derivations killing e."""
    f, e, N, prod_ref, _, _ = _cone_setup(kind, n)
    star = lambda v: -jacobian(f, v)
    H = lambda v: hessian(f, v)
    assert torch.allclose(star(e), e) and torch.allclose(H(e), torch.eye(N, dtype=F64))
    g = gen(90)
    x = e + 0.3 * torch.randn(N, dtype=F64, generator=g)
    xs = star(x)
    assert torch.allclose(star(xs), x) and torch.allclose(H(xs), torch.linalg.inv(H(x)))
    assert torch.allclose(H(x) @ x, xs) and torch.linalg.eigvalsh(H(x)).min() > 0
    T = jacobian(lambda v: hessian(f, v, create_graph=True), x)
    assert torch.allclose(torch.einsum('abc,a->bc', T, x), -2 * H(x))
    T = jacobian(lambda v: hessian(f, v, create_graph=True), e)
    prod = lambda u, v: torch.einsum('a,b,abc->c', u, v, -0.5 * T)
    u, v, w = torch.randn(3, N, dtype=F64, generator=g)
    assert torch.allclose(prod(e, v), v) and torch.allclose(prod(u, v), prod(v, u))
    assert torch.allclose(prod(u, v), prod_ref(u, v))
    assert abs(prod(u, v) @ w - v @ prod(u, w)) < 1e-10
    u2 = prod(u, u)
    assert torch.allclose(prod(u, prod(u2, v)), prod(u2, prod(u, v)))
    Lm = lambda a: torch.stack([prod(a, b) for b in torch.eye(N, dtype=F64)], 1)
    Dm = Lm(u) @ Lm(v) - Lm(v) @ Lm(u)
    assert torch.allclose(Dm, -Dm.T) and torch.allclose(Dm @ e, torch.zeros(N, dtype=F64), atol=1e-12)
    assert torch.allclose(Dm @ prod(w, x), prod(Dm @ w, x) + prod(w, Dm @ x))


def test_prop_F_11_3_H_inverse_is_quadratic_representation():
    """Proposition F.11(3) / Example F.13 in Herm_3(R): H(x)^{-1} = g g^T with g = P(x^{1/2}/sqrt(kappa)) (g e = x),
    i.e. H(x)^{-1} u = x u x / kappa^2, which lies in G(Omega)."""
    f, e, N, _, mat_, vec_ = _cone_setup('sdp', 3)
    kappa2 = 2.0
    g = gen(91)
    G = torch.randn(3, 3, dtype=F64, generator=g)
    xm = G @ G.T + 0.5 * torch.eye(3, dtype=F64)
    x = vec_(xm)
    Hinv = torch.linalg.inv(hessian(f, x))
    u = torch.randn(N, dtype=F64, generator=g)
    assert torch.allclose(mat_(Hinv @ u), xm @ mat_(u) @ xm / kappa2)
    w, U = torch.linalg.eigh(xm)
    r = U @ torch.diag(w.sqrt()) @ U.T / math.sqrt(math.sqrt(kappa2))
    gm = torch.stack([vec_(r @ mat_(b) @ r) for b in torch.eye(N, dtype=F64)], 1)
    assert torch.allclose(gm @ e, x) and torch.allclose(gm @ gm.T, Hinv)


def test_lemma_F_10_one_parameter_groups():
    """Lemma F.10: c(s)^n with c(s) = I + sZ + o(s) converges to e^{tZ}; e^{sZ}e^{sW} = I + s(Z+W) + O(s^2) and the group
    commutator with a = sqrt(s) is I + s[Z,W] + O(s^{3/2})."""
    g = gen(92)
    Z, W = torch.randn(2, 4, 4, dtype=F64, generator=g)
    t = 0.7
    c = lambda s: torch.eye(4, dtype=F64) + s * Z + s ** 1.5 * W                     # o(s) perturbation
    errs = [err(torch.linalg.matrix_power(c(t / n), n), torch.linalg.matrix_exp(t * Z)) for n in (10, 100, 1000)]
    assert errs[0] > errs[1] > errs[2] and errs[2] < 0.1
    ex = torch.linalg.matrix_exp
    rem1, rem2 = [], []
    for s in (1e-2, 1e-4, 1e-6):
        a = math.sqrt(s)
        gc = ex(a * Z) @ ex(a * W) @ ex(-a * Z) @ ex(-a * W)
        rem1.append((gc - torch.eye(4, dtype=F64) - s * (Z @ W - W @ Z)).abs().max().item() / s ** 1.5)
        rem2.append((ex(s * Z) @ ex(s * W) - torch.eye(4, dtype=F64) - s * (Z + W)).abs().max().item() / s ** 2)
    assert max(rem1) / min(rem1) < 1.1 and max(rem2) / min(rem2) < 1.1                # exactly O(s^{3/2}), O(s^2)
    assert rem1[-1] * 1e-6 ** 0.5 < 0.01 * (Z @ W - W @ Z).abs().max().item()             # o(s): [Z, W] is the derivative


def projections_in(alg, g, B=1):
    """Random projections of every rank in each block (rank drawn per block), as operators."""
    M = torch.randn(B, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype, generator=g)
    M = alg.operator(M + M.mH, is_self_adjoint=True)
    blocks = []
    for c, (w, V) in enumerate(M.eigh()):
        k = alg.k_factors[c]
        r = torch.randint(0, k + 1, (1,), generator=g).item()
        P = V[..., :, :r] @ V[..., :, :r].mH
        blocks.append(P)
    return alg.from_blocks(blocks)


@pytest.mark.parametrize('ns,ks', [([2, 3], [2, 3]), ([3, 1, 2], [3, 1, 2]), ([4], [4])])
def test_lemma_F_14_projections_and_order(ns, ks):
    """Lemma F.14: projections are extreme in [0,1] and a non-projection a is the midpoint of a ± eps q_lambda;
    pq = 0 iff p + q <= 1 iff p∘q = 0; p <= q iff pq = p iff p∘q = p; central projections are sums of block units."""
    alg = TypeIAlgebra(ns, ks, precision='double')
    g = gen(93)
    a = alg.operator_from_eigenvalues(lambda k: torch.rand(k, dtype=F64, generator=g) * 0.8 + 0.1, force_self_adjoint=True)
    wa, Va = a.eigh()[0]
    q = Va[..., :, :1] @ Va[..., :, :1].mH
    lam = wa[..., 0].item()
    eps = min(lam, 1 - lam)
    for sgn in (1, -1):
        b = a.matrix[0, 0, :ks[0], :ks[0]] + sgn * eps * q[0]
        wb = torch.linalg.eigvalsh(b)
        assert wb.min() > -1e-12 and wb.max() < 1 + 1e-12
    I = alg.identity(1).matrix
    H = sa_op(alg, 1, g)
    eigvecs = [V[0] for _, V in H.eigh()]                                              # a common eigenbasis per block

    def proj(subsets):
        out = torch.zeros_like(I)
        for c, (Vc, S) in enumerate(zip(eigvecs, subsets)):
            if S:
                out[0, c, :ks[c], :ks[c]] = Vc[:, S] @ Vc[:, S].mH
        return out
    U = alg.random_unitary_operator(1).matrix
    seen = {'orth': 0, 'le': 0, 'neither': 0}
    for _ in range(40):
        sub_p = [sorted(torch.randperm(k, generator=g)[:torch.randint(0, k + 1, (1,), generator=g).item()].tolist())
                 for k in ks]
        sub_r = [sorted(torch.randperm(k, generator=g)[:torch.randint(0, k + 1, (1,), generator=g).item()].tolist())
                 for k in ks]
        p, r = proj(sub_p), proj(sub_r)
        if torch.rand(1, generator=g).item() < 0.3:                                     # a generic (non-commuting) pair
            r = U @ r @ U.mH
        assert err(p @ p, p) < 1e-12 and err(r @ r, r) < 1e-12
        cond_pq = err(p @ r, torch.zeros_like(p)) < 1e-10
        cond_sum = all((torch.linalg.eigvalsh((I - p - r)[0, c, :k, :k]) > -1e-10).all() for c, k in enumerate(ks))
        cond_j = err(jm(p, r), torch.zeros_like(p)) < 1e-10
        assert cond_pq == cond_sum == cond_j
        cond_le = err(p @ r, p) < 1e-10
        cond_ord = all((torch.linalg.eigvalsh((r - p)[0, c, :k, :k]) > -1e-10).all() for c, k in enumerate(ks))
        assert cond_le == (err(jm(p, r), p) < 1e-10) == cond_ord
        seen['orth' if cond_pq else 'le' if cond_le else 'neither'] += 1
    assert min(seen.values()) > 0                                                       # all cases occurred
    x = sa_op(alg, 1, g).matrix
    for c in range(alg.C):
        Pc = torch.zeros_like(I)
        Pc[0, c] = I[0, c]
        assert err(Pc @ x @ Pc + (I - Pc) @ x @ (I - Pc), x) < 1e-12
    if max(ks) >= 2:
        c = ks.index(max(ks))
        v = torch.zeros(ks[c], dtype=C128)
        v[0] = v[1] = 1 / math.sqrt(2)
        pm = torch.zeros_like(I)
        pm[0, c, :ks[c], :ks[c]] = torch.outer(v, v)
        assert err(pm @ x @ pm + (I - pm) @ x @ (I - pm), x) > 1e-6


def f18_map(Us, perm, transposed):
    """A map of the form of Theorem F.18 on operator matrices (B, C, k, k): block sigma(c) = U_c x_c^(T) U_c^*."""
    def Phi(Xm):
        out = torch.zeros_like(Xm)
        for c, d in enumerate(perm):
            k = Us[c].shape[0]
            xc = Xm[:, c, :k, :k]
            if c in transposed:
                xc = xc.mT
            out[:, d, :k, :k] = Us[c] @ xc @ Us[c].mH
        return out
    return Phi


def test_theorem_F_15_F_18_corollary_F_19_order_isomorphisms():
    """Theorems F.15, F.18 and Corollary F.19 on M_2 ⊕ M_2 ⊕ M_3 ⊕ M_1: a map that exchanges the M_2 blocks and transposes
    some blocks is unital, positive with positive inverse, a Jordan *-isomorphism mapping projections to projections and
    preserving pq = 0 and xyx; it is completely (2-)positive iff no block of size >= 2 is transposed."""
    ks = [2, 2, 3, 1]
    alg = TypeIAlgebra(ks, ks, precision='double')
    g = gen(94)
    Us = [torch.linalg.qr(torch.randn(k, k, dtype=C128, generator=g))[0] for k in ks]
    perm = [1, 0, 2, 3]
    for transposed, cp in (({0, 2}, False), ({3}, True), (set(), True), ({1}, False)):
        Phi = f18_map(Us, perm, transposed)
        inv_perm = [perm.index(d) for d in range(4)]
        Ui = [Us[inv_perm[d]].mT if inv_perm[d] in transposed else Us[inv_perm[d]].mH for d in range(4)]
        Phi_inv = lambda Y: f18_map(Ui, inv_perm, {perm[c] for c in transposed})(Y)
        I = alg.identity(1).matrix
        assert err(Phi(I), I) < 1e-14
        x, y = sa_op(alg, 10, g), sa_op(alg, 10, g)
        assert err(Phi_inv(Phi(x.matrix)), x) < 1e-12
        X, Y = alg.operator(Phi(x.matrix), is_self_adjoint=True), alg.operator(Phi(y.matrix), is_self_adjoint=True)
        assert err(Phi(jd.jordan_product(x, y).matrix), jd.jordan_product(X, Y)) < 1e-11
        assert err(Phi(x.matrix @ y.matrix @ x.matrix), X.matrix @ Y.matrix @ X.matrix) < 1e-10
        mult = err(Phi(x.matrix @ y.matrix), X.matrix @ Y.matrix) < 1e-10
        assert mult == cp
        p = sa_op(alg, 30, g, positive=True)
        for M_ in (Phi(p.matrix), Phi_inv(p.matrix)):
            w = torch.cat([e for e in alg.operator(M_, is_self_adjoint=True).eigenvalues()], -1)
            assert (w > -1e-10).all()
        pr = projections_in(alg, g, 1).matrix
        assert err(Phi(pr) @ Phi(pr), Phi(pr)) < 1e-12
        pr2 = I - pr
        assert err(Phi(pr) @ Phi(pr2), torch.zeros_like(pr)) < 1e-12
        # 2-positivity of each block map a -> U a^(T) U^*: Choi matrix
        for c, k in enumerate(ks):
            units = torch.eye(k * k, dtype=C128).reshape(k * k, k, k)
            ims = torch.stack([Us[c] @ (E.T if c in transposed else E) @ Us[c].mH for E in units])
            choi = ims.reshape(k, k, k, k).permute(0, 2, 1, 3).reshape(k * k, k * k)
            neg = torch.linalg.eigvalsh(choi)[0] < -1e-12
            assert neg == (c in transposed and k >= 2)


@pytest.mark.parametrize('n', [1, 2, 3, 4])
@pytest.mark.parametrize('transpose', [False, True])
def test_theorem_F_17_recover_unitary(n, transpose):
    """Theorem F.17: a Jordan *-automorphism of M_n is x -> UxU^* or x -> Ux^TU^*; the proof recovers U up to a phase from
    the images of the matrix units (epsilon constant, alpha_ij = u_i conj(u_j)); for n >= 2 the forms are distinct."""
    g = gen(95 + n)
    U0 = torch.linalg.qr(torch.randn(n, n, dtype=C128, generator=g))[0]
    Psi = (lambda A: U0 @ A.mT @ U0.mH) if transpose else (lambda A: U0 @ A @ U0.mH)
    E = lambda i, j: unit_mat(n, i, j)
    W = torch.stack([torch.linalg.eigh(Psi(E(i, i)))[1][:, -1] for i in range(n)], 1)
    Psi1 = lambda A: W.mH @ Psi(A) @ W
    for i in range(n):
        assert err(Psi1(E(i, i)), E(i, i)) < 1e-12
    eps = set()
    for i in range(n):
        for j in range(n):
            if i != j:
                im = Psi1(E(i, j))
                assert err(im @ im, torch.zeros(n, n)) < 1e-12                          # alpha beta = 0
                eps.add(im[j, i].abs().item() > 0.5)
    assert len(eps) <= 1
    tr = bool(eps and eps.pop())
    assert tr == (transpose and n >= 2)
    Psi2 = (lambda A: Psi1(A.mT)) if tr else Psi1
    alpha = torch.stack([Psi2(E(i, 0))[i, 0] for i in range(n)])
    assert torch.allclose(alpha.abs(), torch.ones(n, dtype=F64))
    for i in range(n):
        for j in range(n):
            assert err(Psi2(E(i, j)), alpha[i] * alpha[j].conj() * E(i, j)) < 1e-12
    Uf = W @ torch.diag(alpha)
    X = torch.randn(n, n, dtype=C128, generator=g)
    assert err(Psi(X), Uf @ (X.mT if tr else X) @ Uf.mH) < 1e-12
    ph = (U0.mH @ Uf)[0, 0]
    assert err(U0.mH @ Uf, ph * torch.eye(n, dtype=C128)) < 1e-12 or (transpose and n == 1)
    if n >= 2:
        T = lambda A: A.mT
        assert err(T(E(0, 1) @ E(1, 0)), T(E(0, 1)) @ T(E(1, 0))) == 1


def test_lemma_F_16_triple_product_preserved():
    """Lemma F.16: a Jordan *-isomorphism (here the transpose of the library) satisfies Phi(xyx) = Phi(x)Phi(y)Phi(x) and
    Phi(1) = 1 (xyx = 2x∘(x∘y) - x^2∘y)."""
    alg = TypeIAlgebra([3, 2, 1], [3, 2, 1], precision='double')
    g = gen(96)
    x, y = sa_op(alg, 5, g), sa_op(alg, 5, g)
    T = jd.transpose_map
    assert err(T(jd.quadratic_representation(x, y)), T(x).matrix @ T(y).matrix @ T(x).matrix) < 1e-11
    jp = jd.jordan_product
    assert err(jd.quadratic_representation(x, y), 2 * jp(x, jp(x, y)).matrix - jp(jp(x, x), y).matrix) < 1e-11


def test_exercise_F_1_and_F_2():
    """Exercise F.1: R^N with the coordinatewise product is Euclidean with cone R^N_+ and its unital order isomorphisms are
    the coordinate permutations; Exercise F.2: Psi_t is a unital positive bijection whose inverse is not positive and which
    does not map projections to projections."""
    g = gen(97)
    N = 5
    x, y, z = torch.randn(3, N, dtype=F64, generator=g)
    assert torch.allclose((x * y) @ z, y @ (x * z))
    assert ((x * x) >= 0).all()
    perm = torch.randperm(N, generator=g)
    assert torch.allclose((x * y)[perm], x[perm] * y[perm])
    alg = TypeIAlgebra([1] * N, [1] * N, precision='double')                           # C^N = ⊕ M_1
    Us = [torch.ones(1, 1, dtype=C128) * complex(math.cos(0.3 * c), math.sin(0.3 * c)) for c in range(N)]
    Phi = f18_map(Us, perm.tolist(), {0, 2})
    v = alg.from_blocks([torch.tensor([[float(t)]]) for t in x.tolist()])
    out = Phi(v.matrix)[0, :, 0, 0].real
    assert torch.allclose(out[perm], x)
    n, t = 3, 0.25
    Psi = lambda A: (1 - t) * A + t / n * torch.trace(A) * torch.eye(n, dtype=C128)
    Psi_inv = lambda A: (A - t / n * torch.trace(A) * torch.eye(n, dtype=C128)) / (1 - t)
    X = herm_t(n, g)
    assert err(Psi_inv(Psi(X)), X) < 1e-12 and err(Psi(torch.eye(n, dtype=C128)), torch.eye(n)) < 1e-15
    G = torch.randn(n, n, dtype=C128, generator=g)
    assert torch.linalg.eigvalsh(Psi(G @ G.mH))[0] > 0
    p = torch.zeros(n, n, dtype=C128)
    p[0, 0] = 1
    assert abs(torch.linalg.eigvalsh(Psi_inv(p))[0].item() + t / (n * (1 - t))) < 1e-14
    w = torch.linalg.eigvalsh(Psi(p))
    assert abs(w[0].item() - t / n) < 1e-14 and err(Psi(p) @ Psi(p), Psi(p)) > 1e-3
    assert err(jm(Psi(X), Psi(X)), Psi(X @ X)) > 1e-3                                  # not a Jordan map
