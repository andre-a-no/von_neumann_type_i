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
