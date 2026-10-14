"""
Tests of torch_vn_algebra.inclusions against independent references: dense NumPy block matrices built from the
definition iota_0(x)_d = (+)_c x_c (x) 1_{lambda_dc} (np.kron / scipy block_diag), the axioms of conditional
expectations, exact index formulas (path graphs, factors, masas), Choi matrices, and the spin-chain embedding
x -> x (x) 1_2 on the full Hilbert space (C^2)^{(x) (L+1)}.
"""
import contextlib
import math

import numpy as np
import pytest
import scipy.linalg as sla
import torch

from torch_vn_algebra import TypeIAlgebra, SpinChain, composite as cp
from torch_vn_algebra import inclusions as inc
from torch_vn_algebra.inclusions import Inclusion

B = 3


@pytest.fixture(autouse=True)
def seed():
    torch.manual_seed(0)
    np.random.seed(0)


@pytest.fixture(params=[True, False], ids=['complex', 'real'])
def cplx(request):
    return request.param


@pytest.fixture(params=['double', 'single'])
def prec(request):
    return request.param


def tol(prec):
    return 1e-10 if prec == 'double' else 2e-5


# ----------------------------------------------------------------------
# examples and dense references
# ----------------------------------------------------------------------
LAM = np.array([[1, 2, 0], [0, 1, 1], [2, 0, 1]])
K = [2, 1, 3]


def generic(cplx=True, prec='double', unitary=True, batch=B):
    N = TypeIAlgebra(K, K, complex_valued=cplx, precision=prec)
    m = list(LAM @ np.array(K))
    M = TypeIAlgebra(m, m, complex_valued=cplx, precision=prec)
    U = M.random_unitary_operator(batch) if unitary else None
    return Inclusion(N, M, LAM, unitary=U)


def with_empty_sectors(cplx=True, prec='double'):
    """N and M with an empty sector each (k_1 = 0, m_1 = 0)."""
    N = TypeIAlgebra([2, 1, 1], [2, 0, 1], complex_valued=cplx, precision=prec)
    lam = np.array([[1, 0, 2], [0, 0, 0], [2, 0, 0]])
    M = TypeIAlgebra([4, 1, 4], [4, 0, 4], complex_valued=cplx, precision=prec)
    return Inclusion(N, M, lam, unitary=M.random_unitary_operator(B))


def maxabs(a):
    return np.abs(a).max() if np.size(a) else 0.0


def rand_op(alg, batch=B):
    return alg.operator(torch.randn(batch, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype))


def rand_pos(alg, batch=B):
    G = rand_op(alg, batch).matrix
    return alg.operator(G @ G.mH, is_self_adjoint=True, is_positive=True)


def blocks(op, alg):
    """Active blocks as numpy complex128 arrays: list over sectors of (batch, k, k)."""
    mat = op.matrix.detach().cpu().to(torch.complex128).numpy()
    return [mat[:, c, :k, :k] for c, k in enumerate(alg.k_factors)]


def dense_unitary(I):
    if I.unitary is None:
        return None
    return blocks(I.unitary, I.M)


def dense_embed(I, x_blocks, b):
    """iota(x)_d = u_d blockdiag_c(kron(x_c, 1_lambda_dc)) u_d^* for sample b, from the definition."""
    lam, U = I.matrix, dense_unitary(I)
    out = []
    for d, m_d in enumerate(I.M.k_factors):
        parts = [np.kron(x_blocks[c][b if x_blocks[c].shape[0] > 1 else 0], np.eye(lam[d, c]))
                 for c in range(I.N.C) if lam[d, c] and I.N.k_factors[c]]
        Y = sla.block_diag(*parts) if parts else np.zeros((0, 0))
        if U is not None and m_d:
            u = U[d][b if U[d].shape[0] > 1 else 0]
            Y = u @ Y @ u.conj().T
        out.append(Y)
    return out


def dense_isometries(I, d, c):
    """V_{dc,i}, i < lambda_dc, from the definition: the columns of the summand c of block d, kron(1_k, e_i)."""
    lam, k = I.matrix, I.N.k_factors
    off = sum(lam[d, cc] * k[cc] for cc in range(c))
    Vs = []
    for i in range(lam[d, c]):
        S = np.zeros((I.M.k_factors[d], lam[d, c] * k[c]))
        S[off:off + lam[d, c] * k[c], :] = np.eye(lam[d, c] * k[c])
        Vs.append(S @ np.kron(np.eye(k[c]), np.eye(lam[d, c])[:, [i]]))
    return Vs


def dense_expect(I, w, y_blocks, b):
    """E(y)_c = nu_c^{-1} sum_{d,i} w_d V^* u_d^* y_d u_d V from the definition."""
    lam, U = I.matrix, dense_unitary(I)
    nu = lam.T @ w
    out = []
    for c, k_c in enumerate(I.N.k_factors):
        acc = np.zeros((k_c, k_c), dtype=complex)
        for d in range(I.M.C):
            if not lam[d, c]:
                continue
            y = y_blocks[d][b]
            if U is not None:
                u = U[d][b if U[d].shape[0] > 1 else 0]
                y = u.conj().T @ y @ u
            for V in dense_isometries(I, d, c):
                acc += w[d] * V.T @ y @ V
        out.append(acc / nu[c] if k_c else acc)
    return out


def weights_np(I, trace):
    return I.trace_weights(trace)[0].double().cpu().numpy()


def tau(I, w, op):
    """sum_d w_d Tr(op_d) for weights (D,) or (batch, D)."""
    tr = torch.diagonal(op.matrix, dim1=-2, dim2=-1).sum(-1)
    return (tr * torch.as_tensor(w).to(tr.real.dtype)).sum(-1)


TRACES = ['blunt', 'norm', 'vn', 'markov', 'random']


def trace_spec(I, name):
    if name == 'random':
        return 0.2 + torch.rand(I.M.C, dtype=torch.float64)
    return name


# ----------------------------------------------------------------------
# construction and validation
# ----------------------------------------------------------------------
class TestConstruction:
    def test_rejects_wrong_sizes(self):
        N = TypeIAlgebra([2, 1], [2, 1])
        M = TypeIAlgebra([5], [5])
        with pytest.raises(ValueError, match="unital"):
            Inclusion(N, M, [[1, 2]])                    # 2 + 2 = 4 != 5
        Inclusion(N, M, [[2, 1]])                        # 4 + 1 = 5

    def test_rejects_non_injective_and_bad_entries(self):
        N = TypeIAlgebra([2, 1], [2, 1])
        M = TypeIAlgebra([4], [4])
        with pytest.raises(ValueError, match="injective"):
            Inclusion(N, M, [[2, 0]])
        with pytest.raises(ValueError, match="integer"):
            Inclusion(N, M, [[1.5, 1]])
        with pytest.raises(ValueError, match="shape"):
            Inclusion(N, M, [[2], [0]])

    def test_rejects_field_and_precision_mismatch(self):
        N = TypeIAlgebra([2], [2], complex_valued=False)
        with pytest.raises(ValueError, match="real"):
            Inclusion(N, TypeIAlgebra([4], [4]), [[2]])
        with pytest.raises(ValueError, match="precision"):
            Inclusion(N, TypeIAlgebra([4], [4], complex_valued=False, precision='double'), [[2]])

    def test_rejects_non_unitary(self):
        N = TypeIAlgebra([2], [2], precision='double')
        M = TypeIAlgebra([4], [4], precision='double')
        with pytest.raises(ValueError, match="unitary"):
            Inclusion(N, M, [[2]], unitary=M.identity() * 2.0)

    def test_empty_sectors_need_zero_columns(self):
        N = TypeIAlgebra([2, 1], [2, 0])
        M = TypeIAlgebra([4], [4])
        with pytest.raises(ValueError, match="empty"):
            Inclusion(N, M, [[2, 1]])


# ----------------------------------------------------------------------
# the embedding
# ----------------------------------------------------------------------
class TestEmbedding:
    @pytest.mark.parametrize('make', [generic, with_empty_sectors], ids=['generic', 'empty_sectors'])
    def test_dense(self, cplx, prec, make):
        I = make(cplx, prec)
        x = rand_op(I.N)
        got = blocks(I.embed(x), I.M)
        xb = blocks(x, I.N)
        for b in range(B):
            ref = dense_embed(I, xb, b)
            for d in range(I.M.C):
                assert maxabs(got[d][b] - ref[d]) < tol(prec) * 10
        # padding stays zero
        pad = ~I.M._const('block_mask', torch.bool)
        assert (I.embed(x).matrix[:, pad] == 0).all()

    def test_unital_star_homomorphism(self, cplx, prec):
        I = generic(cplx, prec)
        x, y = rand_op(I.N), rand_op(I.N)
        e = lambda z: I.embed(z).matrix
        t = tol(prec) * 10
        assert (e(x @ y) - e(x) @ e(y)).abs().max() < t
        assert (e(x.adjoint()) - e(x).mH).abs().max() < t
        assert (e(x + y * 0.5) - e(x) - 0.5 * e(y)).abs().max() < t
        assert (e(I.N.identity()) - I.M._const('eye')).abs().max() < t

    def test_flags_and_batches(self):
        I = generic()
        p = rand_pos(I.N, 1)                              # batch 1 against a unitary batch 3
        out = I.embed(p)
        assert out._is_positive and out._is_self_adjoint and out.shape[0] == B
        assert out.lambda_min.min() > -1e-10

    def test_as_channel(self, cplx, prec):
        I = generic(cplx, prec)
        x = rand_op(I.N)
        Phi = I.as_channel()
        assert (Phi(x).matrix - I.embed(x).matrix).abs().max() < tol(prec) * 10
        assert Phi.is_unital(tol=tol(prec) * 10)

    def test_tensor_product_is_kron(self, cplx):
        N = TypeIAlgebra([2, 3], [2, 3], complex_valued=cplx, precision='double')
        I = inc.from_tensor_product(N, 3)
        x = rand_op(N)
        one = I.M.tensor_factors[1].identity(1)
        assert (I.embed(x).matrix - cp.kron(x, one, I.M).matrix).abs().max() < 1e-12
        # a general second factor
        P = TypeIAlgebra([2, 1], [2, 1], complex_valued=cplx, precision='double')
        I2 = inc.from_tensor_product(N, P)
        assert (I2.embed(x).matrix - cp.kron(x, P.identity(1), I2.M).matrix).abs().max() < 1e-12

    def test_blunt_adjoint_is_restriction_of_states(self, cplx):
        N = TypeIAlgebra([2, 3], [2, 3], complex_valued=cplx, precision='double')
        I = inc.from_tensor_product(N, 2)
        y = rand_op(I.M)
        assert (I.as_channel().adjoint('blunt')(y).matrix - cp.partial_trace(y, keep=1).matrix).abs().max() < 1e-12


# ----------------------------------------------------------------------
# conditional expectations
# ----------------------------------------------------------------------
class TestConditionalExpectation:
    @pytest.mark.parametrize('trace', TRACES)
    @pytest.mark.parametrize('make', [generic, with_empty_sectors], ids=['generic', 'empty_sectors'])
    def test_dense(self, cplx, prec, trace, make):
        I = make(cplx, prec)
        if trace == 'markov' and not I.is_connected:
            pytest.skip("disconnected")
        spec = trace_spec(I, trace)
        E = I.conditional_expectation(spec)
        w = weights_np(I, spec)
        y = rand_op(I.M)
        got = blocks(E(y), I.N)
        yb = blocks(y, I.M)
        for b in range(B):
            ref = dense_expect(I, w, yb, b)
            for c in range(I.N.C):
                assert maxabs(got[c][b] - ref[c]) < tol(prec) * 20
        assert (E.as_channel()(y).matrix - E(y).matrix).abs().max() < tol(prec) * 20

    @pytest.mark.parametrize('trace', TRACES)
    def test_axioms(self, cplx, prec, trace):
        I = generic(cplx, prec)
        spec = trace_spec(I, trace)
        E = I.conditional_expectation(spec)
        w = I.trace_weights(spec)
        nu = I.induced_trace(spec)
        t = tol(prec) * 50
        x, n1, n2 = rand_op(I.N), rand_op(I.N), rand_op(I.N)
        y = rand_op(I.M)
        e = I.embed
        assert (E(I.M.identity()).matrix - I.N._const('eye')).abs().max() < t               # E(1) = 1
        assert (E(e(x)).matrix - x.matrix).abs().max() < t                                  # E o iota = id
        assert (E(e(n1) @ y @ e(n2)).matrix - (n1 @ E(y) @ n2).matrix).abs().max() < t      # bimodule
        assert (tau(I, w, e(E(y))) - tau(I, w, y)).abs().max() < t                         # tau o E = tau
        # tau(iota(x) y) = tau|_N(x E(y)) with the restricted trace nu = Lambda^T w
        lhs = tau(I, w, e(x) @ y)
        tr_N = torch.diagonal((x @ E(y)).matrix, dim1=-2, dim2=-1).sum(-1)
        rhs = (tr_N * nu.to(tr_N.real.dtype)).sum(-1)
        assert (lhs - rhs).abs().max() < t
        # idempotent: (iota E)^2 = iota E
        P = lambda z: e(E(z))
        assert (P(P(y)).matrix - P(y).matrix).abs().max() < t
        # positive
        p = E(rand_pos(I.M))
        assert p._is_positive and p.lambda_min.min() > -t

    def test_completely_positive_choi(self, cplx):
        """Choi matrices sum_pq e_pq (x) E(e^d_pq)_c are positive semidefinite for every pair (c, d)."""
        I = generic(cplx, 'double')
        E = I.conditional_expectation(0.2 + torch.rand(I.M.C, dtype=torch.float64))
        mx = I.M.k_max
        for d, m_d in enumerate(I.M.k_factors):
            units = torch.zeros(m_d * m_d, I.M.C, mx, mx, dtype=I.M.dtype)
            p, q = np.divmod(np.arange(m_d * m_d), m_d)
            units[np.arange(m_d * m_d), d, p, q] = 1
            # one inclusion (unitary of sample 0) applied to all matrix units
            I0 = Inclusion(I.N, I.M, I.matrix, unitary=I.M.operator(I.unitary.matrix[:1]))
            img = I0.conditional_expectation(E.weights[0])(I.M.operator(units)).matrix.numpy()
            for c, k_c in enumerate(I.N.k_factors):
                choi = np.zeros((m_d * k_c, m_d * k_c), dtype=complex)
                for j in range(m_d * m_d):
                    Eu = np.zeros((m_d, m_d)); Eu[p[j], q[j]] = 1
                    choi += np.kron(Eu, img[j, c, :k_c, :k_c])
                assert np.linalg.eigvalsh(choi).min() > -1e-12

    @pytest.mark.parametrize('trace', ['blunt', 'norm', 'vn'])
    def test_tensor_product_is_normalised_partial_trace(self, cplx, prec, trace):
        N = TypeIAlgebra([2, 3], [2, 3], complex_valued=cplx, precision=prec)
        m = 3
        I = inc.from_tensor_product(N, m)
        y = rand_op(I.M)
        E = I.conditional_expectation(trace)
        ref = cp.partial_trace(y, keep=1).matrix / m
        assert (E(y).matrix - ref).abs().max() < tol(prec) * 10
        if trace == 'norm':                   # tau|_N is the normalised trace of N: E is iota's 'norm' adjoint
            assert (I.as_channel().adjoint('norm')(y).matrix - ref).abs().max() < tol(prec) * 10

    def test_conditional_expectation_onto_scalars_and_masa(self):
        M = TypeIAlgebra([4], [4], precision='double')
        y = rand_op(M)
        E = inc.scalars(M).conditional_expectation('norm')
        assert (E(y).matrix[:, 0, 0, 0] - torch.diagonal(y.matrix[:, 0], dim1=-2, dim2=-1).sum(-1) / 4).abs().max() < 1e-12
        D = inc.diagonal(4)
        Ed = D.conditional_expectation('norm')
        diag = torch.diagonal(y.matrix[:, 0], dim1=-2, dim2=-1)                 # E = the diagonal part
        assert (Ed(D.M.operator(y.matrix)).matrix[:, :, 0, 0] - diag).abs().max() < 1e-12

    def test_batched_weights_and_unitaries(self):
        I = generic()
        W = 0.2 + torch.rand(B, I.M.C, dtype=torch.float64)
        E = I.conditional_expectation(W)
        y = rand_op(I.M)
        out = E(y).matrix
        for b in range(B):
            Eb = Inclusion(I.N, I.M, I.matrix, unitary=I.M.operator(I.unitary.matrix[b:b + 1]))
            ref = Eb.conditional_expectation(W[b])(I.M.operator(y.matrix[b:b + 1])).matrix
            assert (out[b:b + 1] - ref).abs().max() < 1e-12
        assert I.induced_trace(W).shape == (B, I.N.C)
        assert E.pimsner_popa_constant().shape == (B,)

    def test_rejects_non_faithful_trace(self):
        I = generic()
        with pytest.raises(ValueError, match="positive"):
            I.conditional_expectation(torch.tensor([1.0, 0.0, 1.0], dtype=torch.float64))

    def test_gradcheck(self):
        I = generic(True, 'double', batch=1)
        E = I.conditional_expectation(0.2 + torch.rand(I.M.C, dtype=torch.float64))
        x = torch.randn(1, I.N.C, I.N.k_max, I.N.k_max, dtype=torch.complex128, requires_grad=True)
        y = torch.randn(1, I.M.C, I.M.k_max, I.M.k_max, dtype=torch.complex128, requires_grad=True)
        assert torch.autograd.gradcheck(lambda a: I.embed(I.N.operator(a)).matrix, (x,))
        assert torch.autograd.gradcheck(lambda a: E(I.M.operator(a)).matrix, (y,))


# ----------------------------------------------------------------------
# index
# ----------------------------------------------------------------------
class TestIndex:
    @pytest.mark.parametrize('n,m', [(1, 2), (2, 3), (3, 2)])
    def test_tensor_product(self, n, m):
        N = TypeIAlgebra([n], [n], precision='double')
        assert inc.from_tensor_product(N, m).index() == pytest.approx(m ** 2, rel=1e-14)

    @pytest.mark.parametrize('n', [1, 2, 5])
    def test_scalars_and_masa(self, n):
        M = TypeIAlgebra([n], [n], precision='double')
        assert inc.scalars(M).index() == pytest.approx(n ** 2, rel=1e-14)
        assert inc.diagonal(n).index() == pytest.approx(n, rel=1e-14)

    def test_factor_index_is_jones_index(self):
        N = TypeIAlgebra([3], [3], precision='double')
        M = TypeIAlgebra([12], [12], precision='double')
        assert Inclusion(N, M, [[4]]).index() == pytest.approx((12 / 3) ** 2, rel=1e-14)

    def test_scalars_in_multimatrix_algebra(self):
        """C c (+)_d M_{m_d}: Lambda Lambda^T = m m^T has the eigenvalue |m|^2."""
        M = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')
        assert inc.scalars(M).index() == pytest.approx(14, rel=1e-14)

    @pytest.mark.parametrize('L', range(0, 9))
    def test_chain_tower_path_graph(self, L):
        I = inc.chain_tower(L + 1, L_min=L)[0]
        lam = I.matrix
        # the bipartite graph of Lambda is the path with 2L + 3 vertices
        n = lam.shape[0] + lam.shape[1]
        assert n == 2 * L + 3
        A = np.block([[np.zeros((lam.shape[0],) * 2), lam], [lam.T, np.zeros((lam.shape[1],) * 2)]])
        deg = A.sum(0)
        assert sorted(deg.tolist()) == [1, 1] + [2] * (n - 2) and I.is_connected
        # spectrum of the path graph: 2 cos(pi j / (n + 1)), j = 1..n
        ref = np.sort(2 * np.cos(np.pi * np.arange(1, n + 1) / (n + 1)))
        assert np.abs(np.linalg.eigvalsh(A) - ref).max() < 1e-12
        assert I.index() == pytest.approx(4 * math.cos(math.pi / (2 * L + 4)) ** 2, rel=1e-13, abs=1e-13)

    def test_chain_index_tends_to_four(self):
        idx = [I.index() for I in inc.chain_tower(12)]
        assert all(a < b < 4 for a, b in zip(idx, idx[1:]))

    def test_markov_trace(self):
        I = generic(unitary=False)
        t, s, beta = I.markov_trace()
        lam = I.matrix.astype(float)
        t, s = t.numpy(), s.numpy()
        assert np.all(t > 0) and np.all(s > 0)
        assert np.allclose(lam.T @ t, s, atol=1e-14) and np.allclose(lam @ s, beta * t, atol=1e-13)
        assert np.allclose(lam @ lam.T @ t, beta * t, atol=1e-13)
        assert abs(np.array(I.M.k_factors) @ t - 1) < 1e-14
        assert abs(np.array(I.N.k_factors) @ s - 1) < 1e-14          # tau|_N(1) = 1
        assert beta == pytest.approx(np.linalg.svd(lam, compute_uv=False)[0] ** 2, rel=1e-14)
        assert I.induced_trace('markov').numpy() == pytest.approx(s, rel=1e-14)

    def test_disconnected(self):
        M = TypeIAlgebra([2, 3], [2, 3], precision='double')
        Z = inc.centre(M)
        assert not Z.is_connected and Z.connected_components() == [([0], [0]), ([1], [1])]
        assert Z.component_indices() == pytest.approx([4, 9])
        with pytest.warns(UserWarning, match="not connected"):
            assert Z.index() == pytest.approx(9)
        with pytest.raises(ValueError, match="Markov"):
            Z.markov_trace()
        # the centre is the relative commutant of M in itself and commutes with M
        assert Z.commutator_error(Inclusion(M, M, np.eye(2, dtype=int))).max() < 1e-14


class TestWatatani:
    @pytest.mark.parametrize('trace', TRACES)
    @pytest.mark.parametrize('make', [generic, with_empty_sectors], ids=['generic', 'empty_sectors'])
    def test_quasi_basis(self, cplx, trace, make):
        I = make(cplx, 'double')
        if trace == 'markov' and not I.is_connected:
            pytest.skip("disconnected")
        spec = trace_spec(I, trace)
        E = I.conditional_expectation(spec)
        Q = E.quasi_basis()
        y = rand_op(I.M).matrix
        rec = sum(Q[:, j] @ I.embed(E(I.M.operator(Q[:, j].mH @ y))).matrix for j in range(Q.shape[1]))
        assert (rec - y).abs().max() < 1e-12                              # y = sum_j u_j E(u_j^* y)
        ind_q, ind_f = E.index('quasi_basis').matrix, E.index().matrix
        assert (ind_q - ind_f).abs().max() < 1e-12
        # central: commutes with M, and equal to (Lambda Lambda^T w)_d / w_d
        z = rand_op(I.M).matrix
        assert (ind_f @ z - z @ ind_f).abs().max() < 1e-12
        w = weights_np(I, spec)
        lam = I.matrix
        for d, m_d in enumerate(I.M.k_factors):
            if m_d:
                ref = (lam @ lam.T @ w)[d] / w[d]
                assert np.allclose(ind_f[:, d, :m_d, :m_d].numpy(), ref * np.eye(m_d), atol=1e-12)

    def test_markov_index_is_scalar(self, cplx):
        I = generic(cplx, 'double')
        ind = I.watatani_index('markov').matrix
        assert (ind - I.index() * I.M._const('eye')).abs().max() < 1e-12

    def test_factors_and_known_values(self):
        n, m = 3, 2
        I = inc.from_tensor_product(TypeIAlgebra([n], [n], precision='double'), m)
        assert I.watatani_index('norm').matrix[0, 0, 0, 0].real == pytest.approx(m ** 2, rel=1e-14)
        S = inc.scalars(TypeIAlgebra([4], [4], precision='double'))
        assert S.watatani_index('norm').matrix[0, 0, 0, 0].real == pytest.approx(16, rel=1e-14)
        D = inc.diagonal(5)
        assert (D.watatani_index('norm').matrix - 5 * D.M._const('eye')).abs().max() < 1e-13


class TestPimsnerPopa:
    @staticmethod
    def best_constant(A, P):
        """max{l : A - l P >= 0} for P = xi xi^* by bisection on the smallest eigenvalue (dense)."""
        lo, hi = 0.0, 1.0
        while np.linalg.eigvalsh(A - hi * P).min() >= -1e-13:
            hi *= 2
        for _ in range(60):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if np.linalg.eigvalsh(A - mid * P).min() >= -1e-13 else (lo, mid)
        return lo

    @pytest.mark.parametrize('trace', ['blunt', 'norm', 'markov', 'random'])
    def test_brute_force(self, cplx, trace):
        """lambda(E) = min over rank-one y = xi xi^* of the best constant (attained by generic xi; the
        minimum over all y >= 0 is over the extreme rays)."""
        I = generic(cplx, 'double', unitary=True, batch=1)
        E = I.conditional_expectation(trace_spec(I, trace))
        lam_E = E.pimsner_popa_constant().item()
        mx = I.M.k_max
        per_block = []
        for d, m_d in enumerate(I.M.k_factors):
            vals = []
            for _ in range(4):
                xi = torch.randn(m_d, dtype=I.M.dtype)
                y = torch.zeros(1, I.M.C, mx, mx, dtype=I.M.dtype)
                y[0, d, :m_d, :m_d] = torch.outer(xi, xi.conj())
                img = I.embed(E(I.M.operator(y))).matrix[0, d, :m_d, :m_d].numpy()
                vals.append(self.best_constant(img, y[0, d, :m_d, :m_d].numpy()))
            per_block.append(vals)
        assert min(min(v) for v in per_block) == pytest.approx(lam_E, rel=1e-8)
        # E(y) >= lambda y for random positive y of full rank
        p = rand_pos(I.M, 8)
        diff = I.embed(E(p)).matrix - lam_E * p.matrix
        assert I.M.operator(diff, is_self_adjoint=True).lambda_min.min() > -1e-10

    @pytest.mark.parametrize('trace', ['blunt', 'markov', 'random'])
    def test_complete_constant_from_choi(self, cplx, trace):
        """max{l : iota E - l id completely positive} = 1/||Ind E||, from the Choi matrices of the blocks (d, d)."""
        I = generic(cplx, 'double', unitary=True, batch=1)
        E = I.conditional_expectation(trace_spec(I, trace))
        lc = E.pimsner_popa_constant(complete=True).item()
        assert lc == pytest.approx(1 / I.watatani_index(E).matrix.real.amax().item(), rel=1e-13)
        mx = I.M.k_max
        mins_at, mins_above = [], []
        for d, m_d in enumerate(I.M.k_factors):
            p, q = np.divmod(np.arange(m_d * m_d), m_d)
            units = torch.zeros(m_d * m_d, I.M.C, mx, mx, dtype=I.M.dtype)
            units[np.arange(m_d * m_d), d, p, q] = 1
            img = I.embed(E(I.M.operator(units))).matrix[:, d, :m_d, :m_d].numpy()
            choi = sum(np.kron(np.eye(m_d)[:, [p[j]]] @ np.eye(m_d)[[q[j]], :], img[j]) for j in range(m_d * m_d))
            omega = np.eye(m_d).reshape(-1, 1)
            ident = omega @ omega.T                                     # Choi matrix of the identity map
            mins_at.append(np.linalg.eigvalsh(choi - lc * ident).min())
            mins_above.append(np.linalg.eigvalsh(choi - 1.001 * lc * ident).min())
        assert min(mins_at) > -1e-12 and min(mins_above) < -1e-6

    @pytest.mark.parametrize('n', [2, 3, 5])
    def test_scalars(self, n):
        E = inc.scalars(TypeIAlgebra([n], [n], precision='double')).conditional_expectation('norm')
        assert E.pimsner_popa_constant().item() == pytest.approx(1 / n, rel=1e-14)
        assert E.pimsner_popa_constant(complete=True).item() == pytest.approx(1 / n ** 2, rel=1e-14)

    @pytest.mark.parametrize('n,m', [(3, 2), (2, 3), (2, 2)])
    def test_tensor_product(self, n, m):
        I = inc.from_tensor_product(TypeIAlgebra([n], [n], precision='double'), m)
        E = I.conditional_expectation('norm')
        assert E.pimsner_popa_constant().item() == pytest.approx(1 / (m * min(n, m)), rel=1e-14)
        assert E.pimsner_popa_constant(complete=True).item() == pytest.approx(1 / m ** 2, rel=1e-14)


# ----------------------------------------------------------------------
# relative commutant, basic construction, towers
# ----------------------------------------------------------------------
def commutant_dimension_dense(I):
    """dim {y in M : [y, iota(x)] = 0 for all x in N}, by the null space of the commutator with the images of the
    matrix units of N (dense, sample 0 of the unitary)."""
    gens = []
    for c, k_c in enumerate(I.N.k_factors):
        for a in range(k_c):
            for b in range(k_c):
                x = torch.zeros(1, I.N.C, I.N.k_max, I.N.k_max, dtype=I.N.dtype)
                x[0, c, a, b] = 1
                gens.append(I.embed(I.N.operator(x)).matrix[0].numpy())
    dim = 0
    for d, m_d in enumerate(I.M.k_factors):
        if not m_d:
            continue
        rows = [np.kron(g[d, :m_d, :m_d], np.eye(m_d)) - np.kron(np.eye(m_d), g[d, :m_d, :m_d].T) for g in gens]
        S = np.linalg.svd(np.vstack(rows), compute_uv=False)
        dim += int(np.sum(S < 1e-9)) + max(0, m_d * m_d - len(S))
    return dim


class TestRelativeCommutant:
    @pytest.mark.parametrize('make', [generic, with_empty_sectors], ids=['generic', 'empty_sectors'])
    def test_commutant(self, cplx, make):
        I = make(cplx, 'double')
        R = I.relative_commutant()
        lam = I.matrix
        assert R.N.k_factors == [int(v) for v in lam[lam > 0]]
        assert sum(k * k for k in R.N.k_factors) == int((lam ** 2).sum())
        assert I.commutator_error(R).max() < 1e-13
        I0 = Inclusion(I.N, I.M, lam, unitary=I.M.operator(I.unitary.matrix[:1]))
        assert commutant_dimension_dense(I0) == int((lam ** 2).sum())
        # without the unitary the commutant fails to commute with the conjugated N
        if (lam > 1).any():
            R0 = Inclusion(R.N, R.M, R.matrix, unitary=None)
            assert I.commutator_error(R0).max() > 1e-3

    def test_tensor_product_commutant_is_second_factor(self):
        I = inc.from_tensor_product(TypeIAlgebra([3], [3], precision='double'), 2)
        R = I.relative_commutant()
        y = rand_op(R.N)
        ref = cp.kron(I.M.tensor_factors[0].identity(1), I.M.tensor_factors[1].operator(y.matrix), I.M).matrix
        assert (R.embed(y).matrix - ref).abs().max() < 1e-14


class TestBasicConstruction:
    def test_sizes_index_and_tower(self):
        I = generic(unitary=False)
        tower = I.jones_tower(5)
        lam = I.matrix
        for j, T in enumerate(tower):
            assert np.array_equal(T.matrix, lam if j % 2 == 0 else lam.T)
            assert T.index() == pytest.approx(I.index(), rel=1e-12)
        M1 = tower[1].M
        assert M1.k_factors == list(lam.T @ np.array(I.M.k_factors))
        assert tower[1].N is I.M and tower[2].N is M1
        # Markov trace of M c M_1: t_1 = s / beta
        t, s, beta = I.markov_trace()
        assert torch.allclose(tower[1].markov_trace().t, s / beta, atol=1e-14)

    @pytest.mark.parametrize('trace', ['markov', 'blunt', 'random'])
    def test_jones_projection(self, cplx, trace):
        I = generic(cplx, 'double')
        spec = trace_spec(I, trace)
        j = I.basic_construction()
        e = I.jones_projection(spec).matrix
        E = I.conditional_expectation(spec)
        y, n = rand_op(I.M), rand_op(I.N)
        assert (e @ e - e).abs().max() < 1e-13 and (e - e.mH).abs().max() < 1e-13
        jn = j.embed(I.embed(n)).matrix
        assert (e @ jn - jn @ e).abs().max() < 1e-12                               # e in N'
        lhs = e @ j.embed(y).matrix @ e
        assert (lhs - j.embed(I.embed(E(y))).matrix @ e).abs().max() < 1e-12      # e y e = E(y) e
        # rank of e in sector c of M_1 is k_c (e_N is the projection onto L^2(N), dim sum_c k_c^2)
        ranks = torch.diagonal(e, dim1=-2, dim2=-1).sum(-1).real
        assert torch.allclose(ranks, torch.tensor(I.N.k_factors, dtype=ranks.dtype).expand_as(ranks), atol=1e-12)
        if trace == 'markov':                                                       # tau_1(y e) = tau(y) / beta
            t1 = j.markov_trace().t
            t, _, beta = I.markov_trace()
            lhs = tau(j, t1, j.M.operator(j.embed(y).matrix @ e))
            assert (lhs - tau(I, t, y) / beta).abs().max() < 1e-12

    def test_m1_is_spanned_by_m_e_m(self):
        """M_1 = span j(M) e_N j(M): the dimension equals dim M_1 = sum_c n1_c^2."""
        N = TypeIAlgebra([1, 2], [1, 2], precision='double')
        M = TypeIAlgebra([3, 2], [3, 2], precision='double')
        I = Inclusion(N, M, [[1, 1], [0, 1]], unitary=M.random_unitary_operator(1))
        j = I.basic_construction()
        e = I.jones_projection('markov').matrix[0]
        units = []
        for d, m_d in enumerate(M.k_factors):
            for p in range(m_d):
                for q in range(m_d):
                    x = torch.zeros(1, M.C, M.k_max, M.k_max, dtype=M.dtype)
                    x[0, d, p, q] = 1
                    units.append(j.embed(M.operator(x)).matrix[0])
        span = torch.stack([(a @ e @ b).reshape(-1) for a in units for b in units]).numpy()
        assert np.linalg.matrix_rank(span, tol=1e-9) == sum(n * n for n in j.M.k_factors)


class TestComposition:
    def test_compose_with_unitaries(self, cplx):
        I1 = generic(cplx, 'double')
        I2 = I1.basic_construction()
        P = I2.M
        I2u = Inclusion(I2.N, P, I2.matrix, unitary=P.random_unitary_operator(B))
        comp = I1.compose(I2u)
        assert np.array_equal(comp.matrix, I2.matrix @ I1.matrix)
        x = rand_op(I1.N)
        assert (comp.embed(x).matrix - I2u.embed(I1.embed(x)).matrix).abs().max() < 1e-12
        assert comp.index() == pytest.approx(np.linalg.norm((I2.matrix @ I1.matrix).astype(float), 2) ** 2)

    def test_chain_tower_composition_is_x_tensor_identity(self):
        tower = inc.chain_tower(4, L_min=1)
        comp = tower[0].compose(tower[1]).compose(tower[2])                       # A_1 c A_4
        assert np.array_equal(comp.matrix, tower[2].matrix @ tower[1].matrix @ tower[0].matrix)
        x = rand_op(tower[0].N)
        step = tower[2].embed(tower[1].embed(tower[0].embed(x)))
        assert (comp.embed(x).matrix - step.matrix).abs().max() < 1e-14
        # Lambda of A_1 c A_4: binomial multiplicities binom(3, N - c)
        lam = comp.matrix
        assert all(lam[N, c] == math.comb(3, N - c) if 0 <= N - c <= 3 else lam[N, c] == 0
                   for N in range(5) for c in range(2))

    def test_bratteli_diagram(self):
        tower = inc.chain_tower(3)
        bd = inc.bratteli_diagram(tower)
        assert bd.levels == [[1], [1, 1], [1, 2, 1], [1, 3, 3, 1]]
        assert bd.edges(1) == [(0, 0, 1), (0, 1, 1), (1, 1, 1), (1, 2, 1)]
        assert all(np.array_equal(a, T.matrix) for a, T in zip(bd.matrices, tower))
        with pytest.raises(ValueError, match="tower"):
            inc.bratteli_diagram([tower[0], tower[2]])

    def test_bratteli_draw(self):
        mpl = pytest.importorskip('matplotlib')
        mpl.use('Agg')
        ax = inc.bratteli_diagram(generic(unitary=False).jones_tower(3)).draw()
        assert ax is not None


# ----------------------------------------------------------------------
# spin-chain tower
# ----------------------------------------------------------------------
class TestChainTower:
    def test_dimensions_are_binomials(self):
        tower = inc.chain_tower(8)
        for L, T in enumerate(tower):
            assert T.N.k_factors == [math.comb(L, n) for n in range(L + 1)]
            assert T.M.k_factors == [math.comb(L + 1, n) for n in range(L + 2)]
            assert T.N.charges == list(range(L + 1))
        assert all(a.M is b.N for a, b in zip(tower, tower[1:]))

    @pytest.mark.parametrize('L', [1, 2, 4])
    def test_embedding_is_x_tensor_identity(self, L, cplx):
        """iota(x) = x (x) 1_2 on (C^2)^{(x) (L+1)} (new site L = highest bit), in the basis of SpinChain."""
        T = inc.chain_tower(L + 1, L_min=L, complex_valued=cplx)[0]
        x = rand_op(T.N, 2)
        small, big = SpinChain(L, complex_valued=cplx), SpinChain(L + 1, complex_valued=cplx)
        assert small.algebra.k_factors == T.N.k_factors and big.algebra.k_factors == T.M.k_factors
        H = small.xxz(J=1.0, Delta=0.5)                       # operators of SpinChain(L) embed directly
        assert torch.allclose(T.embed(H).eigenvalues()[0], small.xxz(J=1.0, Delta=0.5).eigenvalues()[0])
        got = T.embed(x).matrix.numpy()
        for b in range(2):
            Xf = np.zeros((2 ** L, 2 ** L), dtype=got.dtype)
            for n, lab in enumerate(small.labels):
                lab = lab.numpy()
                Xf[np.ix_(lab, lab)] = x.matrix[b, n, :len(lab), :len(lab)].numpy()
            Yf = np.kron(np.eye(2), Xf)                           # index s * 2^L + label: x (x) 1 on site L
            for n, lab in enumerate(big.labels):
                lab = lab.numpy()
                assert np.abs(got[b, n, :len(lab), :len(lab)] - Yf[np.ix_(lab, lab)]).max() < 1e-12

    def test_conditional_expectation_is_partial_trace_over_last_site(self):
        """For tau = normalised trace of (C^2)^{(x) (L+1)} restricted to A_{L+1}, E(y) = Tr_{site L}(y) / 2."""
        L = 3
        T = inc.chain_tower(L + 1, L_min=L)[0]
        big, small = SpinChain(L + 1), SpinChain(L)
        w = torch.full((T.M.C,), 2.0 ** -(L + 1), dtype=torch.float64)          # w_d = 1 / 2^{L+1} for all d
        E = T.conditional_expectation(w)
        y = rand_op(T.M, 2)
        got = E(y).matrix.numpy()
        for b in range(2):
            Yf = np.zeros((2 ** (L + 1),) * 2, dtype=complex)
            for n, lab in enumerate(big.labels):
                lab = lab.numpy()
                Yf[np.ix_(lab, lab)] = y.matrix[b, n, :len(lab), :len(lab)].numpy()
            red = np.einsum('sisj->ij', Yf.reshape(2, 2 ** L, 2, 2 ** L)) / 2     # trace over the highest bit
            for n, lab in enumerate(small.labels):
                lab = lab.numpy()
                assert np.abs(got[b, n, :len(lab), :len(lab)] - red[np.ix_(lab, lab)]).max() < 1e-12


# ----------------------------------------------------------------------
# design rules: no host synchronisation, devices
# ----------------------------------------------------------------------
_NAMES = ('item', 'tolist', '__bool__', '__float__', '__int__', '__index__', 'cpu', 'numpy')


@contextlib.contextmanager
def count_host_reads():
    calls, saved = [], {n: getattr(torch.Tensor, n) for n in _NAMES}

    def wrap(name):
        f = saved[name]

        def g(self, *a, **k):
            calls.append(name)
            return f(self, *a, **k)
        return g
    for n in _NAMES:
        setattr(torch.Tensor, n, wrap(n))
    try:
        yield calls
    finally:
        for n, f in saved.items():
            setattr(torch.Tensor, n, f)


def test_no_host_synchronisation():
    I = generic(True, 'double')
    E = I.conditional_expectation('markov')
    x, y = rand_op(I.N), rand_op(I.M)
    I.embed(x).matrix, E(y).matrix                              # index maps are built once, at first use
    with count_host_reads() as calls:
        for _ in range(2):
            I.embed(x).matrix
            E(y).matrix
            E.index().matrix
            E.pimsner_popa_constant()
    assert calls == []


@pytest.mark.skipif(not torch.cuda.is_available(), reason="needs CUDA")
def test_gpu_matches_cpu():
    I = generic(True, 'double', unitary=False)
    Ng = TypeIAlgebra(K, K, precision='double', device='cuda')
    m = list(LAM @ np.array(K))
    Ig = Inclusion(Ng, TypeIAlgebra(m, m, precision='double', device='cuda'), LAM)
    x, y = rand_op(I.N), rand_op(I.M)
    xg, yg = Ng.operator(x.matrix.cuda()), Ig.M.operator(y.matrix.cuda())
    assert (Ig.embed(xg).matrix.cpu() - I.embed(x).matrix).abs().max() < 1e-12
    E, Eg = I.conditional_expectation('markov'), Ig.conditional_expectation('markov')
    assert (Eg(yg).matrix.cpu() - E(y).matrix).abs().max() < 1e-12
