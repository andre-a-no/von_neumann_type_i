"""Tests for maps between sectors and tensor products."""
import pytest
import torch

from torch_vn_algebra import TypeIAlgebra, InterSectorChannel, tensor_product, kron, partial_trace
from torch_vn_algebra import channels as ch
from torch_vn_algebra import composite as cp
from torch_vn_algebra import states as st

B = 5


@pytest.fixture(autouse=True)
def seed():
    torch.manual_seed(0)


@pytest.fixture(params=[True, False], ids=['complex', 'real'])
def cplx(request):
    return request.param


def algs(cplx):
    a = TypeIAlgebra([2, 3, 4], [2, 3, 4], complex_valued=cplx, device='cpu')
    b = TypeIAlgebra([3, 2], [3, 2], complex_valued=cplx, device='cpu')
    return a, b


def rand_sa(alg, batch=B):
    return alg.operator_from_eigenvalues(lambda d: torch.randn(batch, d), batch_size=batch, force_self_adjoint=True)


def close(a, b, atol=1e-4):
    a = torch.as_tensor(a)
    return torch.allclose(a, torch.as_tensor(b).to(a.dtype), atol=atol)


class TestInterSectorChannel:
    def test_random_is_trace_preserving(self, cplx):
        a, b = algs(cplx)
        Phi = ch.random_inter_sector_channel(a, b, 2, B)
        assert Phi.is_trace_preserving()
        rho = st.random_density_matrix(a, B)
        out = Phi(rho)
        assert isinstance(out, st.DensityMatrix) and out.algebra is b
        assert close(out.trace.real, torch.ones(B))
        assert out.lambda_min.min() > -1e-5
        assert close(Phi.transition_matrix().sum(dim=1), torch.ones(B, a.C))
        assert close(Phi.transition_matrix(rho).sum(dim=1), torch.ones(B, a.C))

    @pytest.mark.parametrize('trace', ['Tr_blunt', 'Tr_norm', 'tau_vN'])
    def test_duality(self, cplx, trace):
        a, b = algs(cplx)
        Phi = ch.random_inter_sector_channel(a, b, 2, B)
        X, A = rand_sa(a), rand_sa(b)
        name = {'Tr_blunt': 'blunt', 'Tr_norm': 'norm', 'tau_vN': 'tau_vN'}[trace]
        lhs = getattr(Phi(X) @ A, trace)().real
        rhs = getattr(X @ Phi.adjoint(name)(A), trace)().real
        assert close(lhs, rhs)

    def test_tau_adjoint_of_tp_map_is_not_blunt_adjoint(self, cplx):
        a, b = algs(cplx)
        Phi = ch.random_inter_sector_channel(a, b, 2, 1)
        one = b.identity(1)
        assert Phi.adjoint('blunt')(one).matrix.sub(a.identity(1).matrix).abs().max() < 1e-4  # unital dual
        assert Phi.adjoint('tau_vN')(one).matrix.sub(a.identity(1).matrix).abs().max() > 1e-2

    def test_matches_sector_preserving_channel(self, cplx):
        a, _ = algs(cplx)
        P = ch.random_channel(a, 3, B)
        rho = st.random_density_matrix(a, B)
        assert close(P.to_inter_sector()(rho).matrix, P(rho).matrix)

    def test_composition(self, cplx):
        a, b = algs(cplx)
        P1, P2 = ch.random_inter_sector_channel(a, b, 2, B), ch.random_inter_sector_channel(b, a, 1, B)
        rho = st.random_density_matrix(a, B)
        assert close((P2 @ P1)(rho).matrix, P2(P1(rho)).matrix)
        assert (P2 @ P1).is_trace_preserving()
        P0 = ch.amplitude_damping_channel(a, 0.4)
        assert close((P1 @ P0)(rho).matrix, P1(P0(rho)).matrix)

    def test_from_blocks_particle_loss(self):
        # sectors N = 0, 1, 2 with dims 1, 2, 1; a particle is lost with probability g
        alg = TypeIAlgebra([1, 2, 1], [1, 2, 1], device='cpu')
        g = 0.3
        blocks = {
            (0, 0): [torch.ones(1, 1)],
            (1, 1): [(1 - g) ** 0.5 * torch.eye(2)],
            (0, 1): [g ** 0.5 * torch.tensor([[1.0, 0.0]]), g ** 0.5 * torch.tensor([[0.0, 1.0]])],
            (2, 2): [(1 - g) ** 0.5 * torch.ones(1, 1)],
            (1, 2): [g ** 0.5 * torch.tensor([[1.0], [1.0]]) / 2 ** 0.5],
        }
        Phi = InterSectorChannel.from_blocks(alg, alg, blocks)
        assert Phi.is_trace_preserving()
        T = Phi.transition_matrix()[0]
        expected = torch.tensor([[1.0, g, 0.0], [0.0, 1 - g, g], [0.0, 0.0, 1 - g]])
        assert close(T, expected)

    def test_non_tp_warns(self, cplx):
        a, b = algs(cplx)
        Phi = ch.random_inter_sector_channel(a, b, 1, B)
        half = InterSectorChannel(a, b, 0.5 * Phi.kraus)
        with pytest.warns(UserWarning):
            out = half(st.random_density_matrix(a, B))
        assert not isinstance(out, st.DensityMatrix)


class TestTensorProduct:
    def test_structure(self, cplx):
        a, b = algs(cplx)
        ab = tensor_product(a, b)
        assert ab.C == a.C * b.C
        assert ab.k_factors == [k * m for k in a.k_factors for m in b.k_factors]
        assert cp.sector_index(ab, 2, 1) == 5

    def test_kron_is_multiplicative(self, cplx):
        a, b = algs(cplx)
        ab = tensor_product(a, b)
        A1, A2, B1, B2 = rand_sa(a), rand_sa(a), rand_sa(b), rand_sa(b)
        assert close((kron(A1, B1, ab) @ kron(A2, B2, ab)).matrix, kron(A1 @ A2, B1 @ B2, ab).matrix)
        assert close(kron(A1, B1, ab).trace, A1.trace * B1.trace)

    def test_partial_trace(self, cplx):
        a, b = algs(cplx)
        ab = tensor_product(a, b)
        A, Bo = rand_sa(a), rand_sa(b)
        X = kron(A, Bo, ab)
        assert close(partial_trace(X, keep=1).matrix, A.matrix * Bo.trace.real[:, None, None, None])
        assert close(partial_trace(X, keep=2).matrix, Bo.matrix * A.trace.real[:, None, None, None])

    def test_product_state_marginals(self, cplx):
        a, b = algs(cplx)
        ab = tensor_product(a, b)
        r1, r2 = st.random_density_matrix(a, B), st.random_density_matrix(b, B)
        rho = kron(r1, r2, ab)
        assert isinstance(rho, st.DensityMatrix)
        m1 = partial_trace(rho, keep=1)
        assert isinstance(m1, st.DensityMatrix)
        assert close(m1.matrix, r1.matrix) and close(partial_trace(rho, keep=2).matrix, r2.matrix)

    def test_partial_trace_adjoint_is_embedding(self, cplx):
        a, b = algs(cplx)
        ab = tensor_product(a, b)
        A = rand_sa(a)
        emb = cp.partial_trace_channel(ab, keep=1).adjoint('blunt')(A)
        assert close(emb.matrix, kron(A, b.identity(B), ab).matrix)

    def test_subadditivity_and_mutual_information(self, cplx):
        a, b = algs(cplx)
        ab = tensor_product(a, b)
        rho = st.random_density_matrix(ab, B)
        S = st.von_neumann_entropy
        I = S(partial_trace(rho, 1)) + S(partial_trace(rho, 2)) - S(rho)
        assert torch.all(I >= -1e-4)
