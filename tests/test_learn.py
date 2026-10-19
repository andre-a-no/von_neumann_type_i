"""learn: channels as trainable networks (every parameter value a channel; composition, residual blocks,
post-selection, feed-forward through the centre; training with torch.optim)."""
import math

import pytest
import torch
from torch import nn

from torch_vn_algebra import TypeIAlgebra, SpinChain, channels, learn, measurement, states


@pytest.fixture(scope='module')
def alg():
    return TypeIAlgebra([1, 3, 2], [1, 3, 2], precision='double')


def _traces(rho):
    return torch.diagonal(rho.matrix, dim1=-2, dim2=-1).sum((-2, -1)).real


def test_every_layer_maps_states_to_states(alg):
    torch.manual_seed(0)
    rho = states.random_density_matrix(alg, batch_size=5)
    out_alg = TypeIAlgebra([2, 2], [2, 2], precision='double')
    layers = [learn.ChannelLayer(alg, 2), learn.UnitaryLayer(alg), learn.LindbladLayer(alg, n_jumps=2),
              learn.Residual(learn.ChannelLayer(alg, 1), p=0.3),
              learn.ChannelLayer(alg, 1, algebra_out=out_alg)]
    for layer in layers:
        out = layer(rho)
        assert torch.allclose(_traces(out), torch.ones(5, dtype=torch.float64), atol=1e-10)
        assert torch.linalg.eigvalsh(out.matrix).min() > -1e-10


def test_a_channel_layer_learns_a_channel_with_torch_optim(alg):
    torch.manual_seed(1)
    true = channels.random_channel(alg, kraus_rank=2, batch_size=1)
    X = states.random_density_matrix(alg, batch_size=16)
    Y = true(X)
    net = learn.ChannelLayer(alg, kraus_rank=2)
    opt = torch.optim.Adam(net.parameters(), lr=0.05)
    for _ in range(600):
        loss = (net(X).matrix - Y.matrix).abs().pow(2).sum()
        opt.zero_grad()
        loss.backward()
        opt.step()
    test = states.random_density_matrix(alg, batch_size=4)
    assert (net(test).matrix - true(test).matrix).abs().max() < 1e-6
    assert learn.count_parameters(net) == 2 * 2 * (1 + 9 + 4)


def test_sequential_network_trains_end_to_end(alg):
    torch.manual_seed(2)
    X = states.random_density_matrix(alg, batch_size=8)
    target = channels.random_channel(alg, kraus_rank=3, batch_size=1)(X)
    net = nn.Sequential(learn.UnitaryLayer(alg), learn.Residual(learn.ChannelLayer(alg, 2), p=0.5, learn_p=True),
                        learn.LindbladLayer(alg, n_jumps=1))
    opt = torch.optim.Adam(net.parameters(), lr=0.05)
    first = None
    for _ in range(40):
        loss = learn.infidelity(net(X), target).mean()
        first = loss.item() if first is None else first
        opt.zero_grad()
        loss.backward()
        assert all(p.grad is not None and torch.isfinite(p.grad).all() for p in net.parameters())
        opt.step()
    assert loss.item() < 0.5 * first


def test_residual_blocks_approach_the_semigroup_of_phi_minus_id(alg):
    torch.manual_seed(3)
    layer = learn.ChannelLayer(alg, kraus_rank=2)
    rho = states.random_density_matrix(alg, batch_size=2)
    t, n = 0.7, 400
    out = rho
    block = learn.Residual(layer, p=t / n)
    with torch.no_grad():
        for _ in range(n):
            out = block(out)
        K = layer.channel().kraus[0]                                   # (r, C, k, k): jumps of Phi - id
        jumps = [alg.operator(K[i:i + 1]) for i in range(K.shape[0])]
        from torch_vn_algebra import dynamics
        exact = dynamics.lindblad_evolve(rho, None, jumps, [0.0, t], method='exact')[-1]
    assert (out.matrix - exact.matrix).abs().max() < 5e-3                # O(t^2 / n)


def test_post_selection_renormalises_and_records_the_acceptance(alg):
    rho = states.random_density_matrix(alg, batch_size=3)
    ps = learn.PostSelect(alg, [0, 1, 1])
    out = ps(rho)
    assert torch.allclose(_traces(out), torch.ones(3, dtype=torch.float64))
    p = states.sector_probabilities(rho)[:, 1:].sum(-1)
    assert torch.allclose(ps.acceptance, p)


def test_feed_forward_through_the_centre_is_a_channel():
    torch.manual_seed(4)
    chain = SpinChain(3, 'open', precision='double')
    alg = chain.algebra
    inst = measurement.Instrument.lueders(measurement.POVM.sectors(alg))
    branches = [learn.ChannelLayer(alg, 1) for _ in range(len(inst))]
    ff = learn.FeedForward(inst, branches)
    rho = states.random_density_matrix(alg, batch_size=2)
    out = ff(rho)
    assert torch.allclose(_traces(out), torch.ones(2, dtype=torch.float64), atol=1e-10)
    # identity branches: the non-selective measurement of a central POVM leaves states of the algebra unchanged
    same = learn.FeedForward(inst, [nn.Identity()] * len(inst))(rho)
    assert torch.allclose(same.matrix, rho.matrix)


def test_negative_log_likelihood_is_minimal_at_the_true_state(alg):
    torch.manual_seed(5)
    rho = states.random_density_matrix(alg, batch_size=1)
    U = alg.random_unitary_operator(3)
    povm = measurement.POVM.combine([measurement.POVM.from_unitary(alg.operator(U.matrix[i:i + 1]))
                                     for i in range(3)])
    probs = povm.probabilities(rho)
    counts = (probs * 1e6).round()
    other = states.random_density_matrix(alg, batch_size=1)
    mix = rho.mix(other, 0.1)
    assert learn.negative_log_likelihood(povm, rho, counts) < learn.negative_log_likelihood(povm, mix, counts)
