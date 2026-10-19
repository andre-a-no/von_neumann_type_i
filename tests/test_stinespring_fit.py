"""optimize.StinespringParam (every parameter value is a channel) and optimize.fit (learning a channel from
input/output pairs, a small process tomography)."""
import pytest
import torch

from torch_vn_algebra import TypeIAlgebra, optimize, states, channels


@pytest.fixture(scope='module')
def alg():
    return TypeIAlgebra([1, 3, 2], [1, 3, 2], precision='double')


def _blunt_traces(op):
    return torch.diagonal(op.matrix, dim1=-2, dim2=-1).sum((-2, -1)).real


def test_channels_within_one_algebra_are_trace_preserving(alg):
    P = optimize.StinespringParam(alg, 4, kraus_rank=2)
    K = P.kraus()
    assert K.shape == (4, 2, alg.C, alg.k_max, alg.k_max)
    KdK = (K.conj().transpose(-2, -1) @ K).sum(1)
    eye = alg._const('eye').to(KdK.dtype).expand_as(KdK)
    assert torch.allclose(KdK, eye, atol=1e-12)
    rho = states.random_density_matrix(alg, batch_size=4)
    assert torch.allclose(_blunt_traces(P.channel()(rho)), torch.ones(4, dtype=torch.float64), atol=1e-12)


def test_maps_between_algebras_respect_the_allowed_transitions(alg):
    out = TypeIAlgebra([2, 2], [2, 2], precision='double')
    allowed = torch.tensor([[1, 1, 0], [0, 1, 1]], dtype=torch.bool)
    P = optimize.StinespringParam(alg, 3, kraus_rank=1, alg_out=out, transitions=allowed)
    K = P.kraus()
    assert K.shape == (3, 1, out.C, alg.C, out.k_max, alg.k_max)
    assert K[:, :, ~allowed].abs().max() == 0
    rho = states.random_density_matrix(alg, batch_size=3)
    assert torch.allclose(_blunt_traces(P.channel()(rho)), torch.ones(3, dtype=torch.float64), atol=1e-12)
    with pytest.raises(ValueError):
        optimize.StinespringParam(alg, 1, kraus_rank=1, alg_out=out, transitions=[[0, 1, 0], [0, 0, 1]])


def test_gradient_matches_finite_differences(alg):
    torch.manual_seed(0)
    P = optimize.StinespringParam(alg, 1, kraus_rank=2)
    rho = states.random_density_matrix(alg, batch_size=1)
    target = states.random_density_matrix(alg, batch_size=1)

    def value():
        return states.trace_of_product(P.channel()(rho), target).real.sum()
    g, = torch.autograd.grad(value(), P.Z)
    D = torch.randn_like(P.Z) * P._mask
    eps = 1e-6
    with torch.no_grad():
        Z0 = P.Z.clone()
        P.Z.copy_(Z0 + eps * D)
        up = value()
        P.Z.copy_(Z0 - eps * D)
        down = value()
        P.Z.copy_(Z0)
    fd = (up - down) / (2 * eps)
    assert (g.conj() * D).real.sum().item() == pytest.approx(fd.item(), rel=1e-6)


@pytest.mark.filterwarnings('ignore::torch_vn_algebra.cost.CostWarning')
def test_fit_recovers_a_channel_from_input_output_pairs(alg):
    torch.manual_seed(1)
    true = channels.random_channel(alg, kraus_rank=2, batch_size=1)
    inputs = states.random_density_matrix(alg, batch_size=12)
    outputs = true(inputs).matrix
    starts = 4
    P = optimize.StinespringParam(alg, starts, kraus_rank=2)

    def loss():
        K = P.kraus()                                                    # (starts, r, C, k, k)
        X = K[:, None] @ inputs.matrix[None, :, None]                    # (starts, n, r, C, k, k)
        Y = (X @ K[:, None].conj().transpose(-2, -1)).sum(2)             # (starts, n, C, k, k)
        return (Y - outputs[None]).abs().pow(2).sum((1, 2, 3, 4))
    res = optimize.fit(loss, [P], rounds=8, steps=150, lr=0.05)
    assert res.loss[res.best] < 1e-10
    fitted = P.channel()
    rho = states.random_density_matrix(alg, batch_size=1)
    diff = fitted(rho).matrix[res.best] - true(rho).matrix[0]
    assert diff.abs().max() < 1e-5


@pytest.mark.filterwarnings('ignore::torch_vn_algebra.cost.CostWarning')
def test_lindblad_param_learns_a_qubit_generator_from_trajectories():
    from torch_vn_algebra import dynamics
    torch.manual_seed(2)
    q = TypeIAlgebra([2], [2], precision='double')
    sx = torch.tensor([[0, 1], [1, 0]], dtype=torch.complex128)
    sz = torch.tensor([[1, 0], [0, -1]], dtype=torch.complex128)
    sm = torch.tensor([[0, 0], [1, 0]], dtype=torch.complex128)
    op = lambda M: q.operator(M.reshape(1, 1, 2, 2))
    H_true, jumps_true, rates = op(0.5 * sx + 0.3 * sz), [op(sm), op(sz)], [0.2, 0.05]
    rho0 = states.random_density_matrix(q, batch_size=6)
    times = [0.0, 0.7, 1.5]
    data = torch.stack([r.matrix for r in dynamics.lindblad_evolve(rho0, H_true, jumps_true, times, rates=rates,
                                                                     method='exact')[1:]])   # (2, 6, 1, 2, 2)
    starts = 4
    P = optimize.LindbladParam(q, starts, n_jumps=2)
    init = rho0.matrix.repeat(starts, 1, 1, 1)                          # start s, state i -> index s * 6 + i

    def loss():
        H = q.operator(P.hamiltonian().matrix.repeat_interleave(6, 0), is_self_adjoint=True)
        L = [q.operator(J.matrix.repeat_interleave(6, 0)) for J in P.jumps()]
        ys = dynamics.lindblad_evolve(q.operator(init, is_self_adjoint=True), H, L, times, method='exact')[1:]
        diff = torch.stack([y.matrix for y in ys]) - data.repeat(1, starts, 1, 1, 1)
        return diff.abs().pow(2).reshape(len(times) - 1, starts, -1).sum((0, 2))
    res = optimize.fit(loss, [P], rounds=8, steps=100, lr=0.05)
    assert res.loss[res.best] < 1e-8
    test = states.random_density_matrix(q, batch_size=1)
    pred = dynamics.lindblad_evolve(test, q.operator(P.hamiltonian().matrix[res.best:res.best + 1],
                                                     is_self_adjoint=True),
                                    [q.operator(J.matrix[res.best:res.best + 1]) for J in P.jumps()], [0.0, 2.0],
                                    method='exact')[-1].matrix
    true = dynamics.lindblad_evolve(test, H_true, jumps_true, [0.0, 2.0], rates=rates, method='exact')[-1].matrix
    assert (pred - true).abs().max() < 1e-4


def test_lindblad_param_between_sectors_gives_valid_generators():
    from torch_vn_algebra import SpinChain, dynamics
    chain = SpinChain(3, 'open', precision='double')
    C = chain.algebra.C
    lowering = torch.zeros(C, C, dtype=torch.bool)
    for N in range(1, C):
        lowering[N - 1, N] = True                                       # output sector N-1 from input sector N
    P = optimize.LindbladParam(chain.algebra, 2, n_jumps=1, n_inter=2, transitions=lowering)
    K = P.jumps()[1].kraus
    assert K[:, :, ~lowering].abs().max() == 0
    rho = dynamics.lindblad_evolve(chain.basis_state('110'), P.hamiltonian(), P.jumps(), [0.0, 1.0],
                                   method='exact')[-1]
    tr = torch.diagonal(rho.matrix, dim1=-2, dim2=-1).sum((-2, -1)).real
    assert torch.allclose(tr, torch.ones(2, dtype=torch.float64), atol=1e-10)
    assert torch.linalg.eigvalsh(rho.matrix).min() > -1e-10
