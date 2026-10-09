"""Tests for functional calculus, quantum channels, states and dynamics."""
import pytest
import torch

from torch_vn_algebra import TypeIAlgebra, Channel
from torch_vn_algebra import channels as ch
from torch_vn_algebra import dynamics as dy
from torch_vn_algebra import states as st

B = 6


@pytest.fixture(autouse=True)
def seed():
    torch.manual_seed(0)


@pytest.fixture(params=[True, False], ids=['complex', 'real'])
def alg(request):
    return TypeIAlgebra([3, 5], [3, 4], complex_valued=request.param, device='cpu')


@pytest.fixture
def calg():
    return TypeIAlgebra([3, 5], [3, 4], complex_valued=True, device='cpu')


def random_sa(alg, batch=B):
    return alg.operator_from_eigenvalues(lambda d: torch.randn(batch, d), batch_size=batch,
                                         force_self_adjoint=True)


def random_pos(alg, batch=B, low=0.05):
    return alg.operator_from_eigenvalues(lambda d: low + torch.rand(batch, d), batch_size=batch,
                                         force_self_adjoint=True, force_positive=True)


def close(a, b, atol=1e-4):
    return torch.allclose(torch.as_tensor(a), torch.as_tensor(b).to(torch.as_tensor(a).dtype), atol=atol)


def all_channels(alg):
    return {
        'identity': ch.identity_channel(alg),
        'dephasing': ch.dephasing_channel(alg, 0.3),
        'depolarizing': ch.depolarizing_channel(alg, 0.4),
        'amplitude_damping': ch.amplitude_damping_channel(alg, 0.6),
        'random': ch.random_channel(alg, 3, B),
        'mixed_unitary': ch.random_mixed_unitary_channel(alg, 4, B),
        'center': ch.center_expectation(alg),
    }


# ----------------------------------------------------------------------
# Functional calculus
# ----------------------------------------------------------------------
class TestFunctionalCalculus:
    def test_apply_function_matches_sqrt(self, alg):
        P = random_pos(alg)
        assert close(P.apply_function(torch.sqrt).matrix, P.sqrt().matrix)

    def test_log_exp_roundtrip(self, alg):
        P = random_pos(alg)
        assert close(P.log().expm().matrix, P.matrix)

    def test_padding_stays_zero(self, alg):
        H = random_sa(alg)
        F = H.apply_function(lambda w: w + 1.0).matrix   # f(0) != 0
        assert torch.all(F[:, 0, 3:, :] == 0) and torch.all(F[:, 0, :, 3:] == 0)

    def test_expm_unitary(self, calg):
        H = random_sa(calg)
        U = H.expm(-0.8j).matrix
        assert close(U, H.apply_function(lambda w: torch.exp(-0.8j * w)).matrix)
        UU = U @ U.conj().transpose(-2, -1)
        assert close(UU, calg.identity(B).matrix.to(UU.dtype))

    def test_complex_function_on_real_algebra_raises(self):
        alg = TypeIAlgebra([2], [2], complex_valued=False, device='cpu')
        H = random_sa(alg)
        with pytest.raises(ValueError):
            H.apply_function(lambda w: torch.exp(-1j * w)).matrix

    def test_power(self, alg):
        P = random_pos(alg)
        assert close(P.power(2.0).matrix, (P @ P).matrix)
        assert close((P.power(-1.0) @ P).matrix, alg.identity(B).matrix)


# ----------------------------------------------------------------------
# Channels
# ----------------------------------------------------------------------
class TestChannels:
    def test_trace_preserving_and_cp(self, alg):
        for name, Phi in all_channels(alg).items():
            assert Phi.is_trace_preserving(), name
            for J in Phi.choi():
                assert torch.linalg.eigvalsh(J).min() > -1e-5, name

    def test_unitality(self, alg):
        flags = {name: Phi.is_unital() for name, Phi in all_channels(alg).items()}
        assert flags['identity'] and flags['dephasing'] and flags['depolarizing']
        assert flags['mixed_unitary'] and flags['center']
        assert not flags['amplitude_damping'] and not flags['random']

    def test_maps_states_to_states(self, alg):
        rho = st.random_density_matrix(alg, B)
        for name, Phi in all_channels(alg).items():
            out = Phi(rho)
            assert close(out.trace.real, torch.ones(B)), name
            assert out.lambda_min.min() > -1e-5, name

    @pytest.mark.parametrize('trace', ['Tr_blunt', 'Tr_norm', 'tau_vN'])
    def test_heisenberg_duality(self, alg, trace):
        rho, A = random_pos(alg), random_sa(alg)
        for name, Phi in all_channels(alg).items():
            lhs = getattr(Phi(rho) @ A, trace)().real
            rhs = getattr(rho @ Phi.adjoint()(A), trace)().real
            assert close(lhs, rhs), (name, trace)

    def test_compose(self, alg):
        P1, P2 = ch.random_channel(alg, 2, B), ch.amplitude_damping_channel(alg, 0.3)
        rho = st.random_density_matrix(alg, B)
        assert close((P1 @ P2)(rho).matrix, P1(P2(rho)).matrix)
        assert (P1 @ P2).is_trace_preserving()

    def test_superoperator_roundtrip(self, alg):
        rho = st.random_density_matrix(alg, B)
        for name, Phi in all_channels(alg).items():
            R = Channel.from_superoperator(alg, Phi.superoperator())
            assert close(R(rho).matrix, Phi(rho).matrix), name

    def test_non_cp_superoperator_rejected(self, alg):
        # transpose map: positive but not completely positive
        S = []
        for k in alg.k_factors:
            T = torch.zeros(k * k, k * k)
            for i in range(k):
                for j in range(k):
                    T[j * k + i, i * k + j] = 1.0
            S.append(T)
        with pytest.raises(ValueError):
            Channel.from_superoperator(alg, S)

    def test_center_expectation(self, alg):
        A = random_sa(alg)
        E = ch.center_expectation(alg)(A).matrix
        for c, k in enumerate(alg.k_factors):
            mean = torch.diagonal(A.matrix[:, c, :k, :k], dim1=-2, dim2=-1).sum(-1) / k
            assert close(E[:, c, :k, :k], mean[:, None, None] * torch.eye(k))

    def test_michelson_contrast_decreases_under_unital_channels(self, alg):
        A = random_pos(alg)
        for Phi in (ch.random_mixed_unitary_channel(alg, 3, B), ch.dephasing_channel(alg, 0.7),
                    ch.depolarizing_channel(alg, 0.2)):
            assert torch.all(Phi(A).michelson_contrast <= A.michelson_contrast + 1e-5)

    def test_michelson_contrast_can_increase_without_unitality(self, alg):
        A = alg.identity(1)                      # Delta = 0
        out = ch.amplitude_damping_channel(alg, 0.5)(A)
        assert out.michelson_contrast.item() > 0.1

    def test_lueders_measurement_channel(self, alg):
        P0 = alg.from_blocks([torch.diag(torch.tensor([1.0, 0.0, 0.0])), torch.zeros(4, 4)])
        P1 = alg.identity(1) - P0
        rho = st.random_density_matrix(alg, B)
        probs = st.born_probabilities(rho, [P0, P1])
        assert close(probs.sum(-1), torch.ones(B))
        post, p = st.lueders_update(rho, P0)
        assert close(p, probs[:, 0]) and close(post.trace.real, torch.ones(B))
        assert ch.measurement_channel([P0, P1]).is_trace_preserving()


# ----------------------------------------------------------------------
# States
# ----------------------------------------------------------------------
class TestStates:
    @pytest.mark.parametrize('measure', ['hs', 'bures'])
    def test_random_density(self, alg, measure):
        rho = st.random_density_matrix(alg, B, measure=measure)
        assert close(rho.trace.real, torch.ones(B))
        assert rho.lambda_min.min() > -1e-6
        assert close(st.sector_probabilities(rho).sum(-1), torch.ones(B))

    def test_low_rank(self, alg):
        rho = st.random_density_matrix(alg, B, rank=1)
        assert torch.all(st.purity(rho) <= 1.0 + 1e-5)
        ranks = [torch.linalg.matrix_rank(rho.matrix[:, c]) for c in range(alg.C)]
        assert all(torch.all(r <= 1) for r in ranks)

    def test_gibbs(self, alg):
        H = random_sa(alg)
        g0 = st.gibbs_state(H, 0.0)
        assert close(g0.matrix, st.maximally_mixed_state(alg, B).matrix)
        g = st.gibbs_state(H, 50.0)                 # close to the ground state
        assert close(st.expectation(g, H).real, H.lambda_min, atol=1e-2)

    def test_tracial_state_reproduces_tau(self, alg):
        A = random_sa(alg)
        assert close(st.expectation(st.tracial_state(alg, B), A).real, A.tau_vN().real)

    def test_entropy_and_distances(self, alg):
        rho, sigma = st.random_density_matrix(alg, B), st.random_density_matrix(alg, B)
        assert close(st.von_neumann_entropy(rho), rho.entropy())
        assert close(st.fidelity(rho, rho), torch.ones(B))
        assert torch.all(st.trace_distance(rho, sigma) > 0)
        assert close(st.relative_entropy(rho, rho), torch.zeros(B))
        # Fuchs - van de Graaf
        F, T = st.fidelity(rho, sigma), st.trace_distance(rho, sigma)
        assert torch.all(1 - F.sqrt() <= T + 1e-4) and torch.all(T <= (1 - F).sqrt() + 1e-4)

    def test_data_processing(self, alg):
        rho, sigma = st.random_density_matrix(alg, B), st.random_density_matrix(alg, B)
        Phi = ch.random_channel(alg, 2, B)
        assert torch.all(st.relative_entropy(Phi(rho), Phi(sigma)) <= st.relative_entropy(rho, sigma) + 1e-4)
        assert torch.all(st.trace_distance(Phi(rho), Phi(sigma)) <= st.trace_distance(rho, sigma) + 1e-5)


# ----------------------------------------------------------------------
# Dynamics
# ----------------------------------------------------------------------
def random_ket(alg, batch=B):
    psi = torch.randn(batch, alg.C, alg.k_max, dtype=torch.complex64)
    for c, k in enumerate(alg.k_factors):
        psi[:, c, k:] = 0
    return psi / psi.abs().pow(2).sum((1, 2), keepdim=True).sqrt()


class TestDynamics:
    def test_schrodinger_exact_vs_rk4(self, calg):
        H, psi = random_sa(calg), random_ket(calg)
        ts = [0.0, 0.5, 1.5]
        exact = dy.schrodinger(psi, H, ts)
        assert close(exact, dy.schrodinger_rk4(psi, lambda t: H, ts, substeps=100))
        assert close(exact.abs().pow(2).sum((2, 3)), torch.ones(3, B))
        U = dy.propagator(H, 1.5).matrix
        assert close((U @ psi.unsqueeze(-1)).squeeze(-1), exact[-1])

    def test_sectors_are_conserved(self, calg):
        H, rho = random_sa(calg), st.random_density_matrix(calg, B)
        p0 = st.sector_probabilities(rho)
        for r in dy.von_neumann(rho, H, [0.3, 2.0]):
            assert close(st.sector_probabilities(r), p0)

    def test_energy_conserved(self, calg):
        H, rho = random_sa(calg), st.random_density_matrix(calg, B)
        e0 = st.expectation(rho, H).real
        for r in dy.von_neumann(rho, H, [1.0, 3.0]):
            assert close(st.expectation(r, H).real, e0)

    def test_lindblad_exact_vs_rk4(self, calg):
        H, rho = random_sa(calg), st.random_density_matrix(calg, B)
        jumps, rates = [random_sa(calg), 0.5 * calg.random_unitary_operator(B)], [0.3, 0.7]
        Phi = dy.lindblad_channel(calg, H, jumps, t=1.0, rates=rates)
        assert Phi.is_trace_preserving(1e-4)
        rk = dy.lindblad_evolve(rho, H, jumps, [0.0, 1.0], rates=rates, substeps=100)[-1]
        assert close(Phi(rho).matrix, rk.matrix)

    def test_lindblad_without_jumps_is_unitary(self, calg):
        H, rho = random_sa(calg), st.random_density_matrix(calg, B)
        Phi = dy.lindblad_channel(calg, H, [], t=0.7)
        assert close(Phi(rho).matrix, dy.von_neumann(rho, H, [0.7])[0].matrix)

    def test_pure_dephasing_converges_to_center(self, alg):
        # jump operators = number-like operators diag(0, 1, ..., k-1): kill coherences
        blocks = [torch.diag(torch.arange(k, dtype=torch.float32)) for k in alg.k_factors]
        N = alg.from_blocks(blocks)
        rho = st.random_density_matrix(alg, B)
        final = dy.lindblad_channel(alg, None, [N], t=40.0)(rho).matrix
        diag = torch.diag_embed(torch.diagonal(rho.matrix, dim1=-2, dim2=-1))
        assert close(final, diag, atol=1e-4)

    def test_solve_operator_ode_heisenberg(self, calg):
        H, A = random_sa(calg), random_sa(calg)
        Hm = H.matrix
        ts = [0.0, 1.0]
        At = dy.solve_operator_ode(lambda t, X: 1j * (Hm @ X - X @ Hm), A, ts, substeps=100)[-1]
        U = dy.propagator(H, 1.0).matrix
        assert close(At.matrix, U.conj().transpose(-2, -1) @ A.matrix.to(U.dtype) @ U)

    def test_real_algebra_rejects_hamiltonian_dynamics(self):
        alg = TypeIAlgebra([2], [2], complex_valued=False, device='cpu')
        with pytest.raises(ValueError):
            dy.propagator(random_sa(alg), 1.0)


# ----------------------------------------------------------------------
# DensityMatrix and laziness
# ----------------------------------------------------------------------
class TestDensityMatrix:
    def test_constructors_return_states(self, alg):
        H = random_sa(alg)
        for rho in (st.random_density_matrix(alg, B), st.maximally_mixed_state(alg, B),
                    st.tracial_state(alg, B), st.gibbs_state(H, 1.0)):
            assert isinstance(rho, st.DensityMatrix)
            assert close(rho.trace.real, torch.ones(rho.matrix.shape[0]))

    def test_validation(self, alg):
        bad = 2 * st.maximally_mixed_state(alg, 1).matrix
        with pytest.raises(ValueError):
            st.DensityMatrix(alg, matrix=bad)
        neg = alg.from_blocks([torch.diag(torch.tensor([1.5, -0.5, 0.0])), torch.zeros(4, 4)]).matrix
        with pytest.raises(ValueError):
            st.DensityMatrix(alg, matrix=neg)

    def test_type_propagation(self, alg):
        rho, A = st.random_density_matrix(alg, B), random_sa(alg)
        assert isinstance(ch.random_channel(alg, 2, B)(rho), st.DensityMatrix)
        assert isinstance(rho.mix(st.maximally_mixed_state(alg, B), 0.3), st.DensityMatrix)
        assert not isinstance(rho @ A, st.DensityMatrix)
        assert not isinstance(rho - rho, st.DensityMatrix)
        assert isinstance(rho.condition_on(alg.identity(1))[0], st.DensityMatrix)

    def test_non_trace_preserving_channel_warns(self, alg):
        rho = st.random_density_matrix(alg, B)
        half = Channel(alg, 0.5 * ch.identity_channel(alg).kraus)
        with pytest.warns(UserWarning):
            out = half(rho)
        assert not isinstance(out, st.DensityMatrix)

    def test_dynamics_keep_states(self, calg):
        H, rho = random_sa(calg), st.random_density_matrix(calg, B)
        assert all(isinstance(r, st.DensityMatrix) for r in dy.von_neumann(rho, H, [0.5, 1.0]))
        assert all(isinstance(r, st.DensityMatrix) for r in dy.lindblad_evolve(rho, H, [], [0.0, 0.5]))

    @pytest.mark.parametrize('trace', ['blunt', 'norm', 'tau_vN'])
    def test_density_conventions(self, alg, trace):
        rho, A = st.random_density_matrix(alg, B), random_sa(alg)
        d = alg.operator(rho.density(trace))
        assert close(getattr(d @ A, {'blunt': 'Tr_blunt', 'norm': 'Tr_norm', 'tau_vN': 'tau_vN'}[trace])().real,
                     rho.expectation(A).real)
        back = st.DensityMatrix.from_density(alg, rho.density(trace), trace=trace)
        assert close(back.matrix, rho.matrix)

    def test_trace_of_product(self, alg):
        X, Y = random_sa(alg), random_pos(alg)
        assert close(st.trace_of_product(X, Y), (X @ Y).trace)


def test_materialized_operator_releases_parents():
    import gc
    import weakref
    alg = TypeIAlgebra([4, 4], [4, 4], device='cpu')
    X, Y = random_sa(alg), random_sa(alg)
    mid = X @ Y
    ref = weakref.ref(mid)
    Z = mid @ X
    del mid
    Z.matrix
    gc.collect()
    assert ref() is None


def test_lindblad_channel_rejects_jumps_between_sectors():
    from torch_vn_algebra import SpinChain
    chain = SpinChain(3)
    with pytest.raises(ValueError):
        dy.lindblad_channel(chain.algebra, chain.xxz(1.0, 1.0), [chain.lowering(0)])
