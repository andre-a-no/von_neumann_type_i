"""
dynamics.trajectories (Monte Carlo wave functions): the averages over the trajectories agree with the solution of the
Lindblad equation within their statistical errors, for jumps inside the sectors, between sectors (InterSectorChannel,
SectorJump) and for mixed initial states; the result does not depend on the storage layout.
"""
import math

import pytest
import torch

from torch_vn_algebra import TypeIAlgebra, InterSectorChannel, SpinChain, dynamics, states

C128 = torch.complex128


def gen(seed):
    return torch.Generator().manual_seed(seed)


def test_spontaneous_decay():
    alg = TypeIAlgebra([2], [2], precision='double')
    sm = alg.operator(torch.tensor([[[[0, 1.0], [0, 0]]]], dtype=C128))
    Pe = alg.operator(torch.tensor([[[[0, 0], [0, 1.0]]]], dtype=C128), is_self_adjoint=True)
    times = [0.0, 0.5, 1.0, 2.0]
    tr = dynamics.trajectories(torch.tensor([[0, 1.0]], dtype=C128), None, [sm], times, rates=[0.7], n_traj=4000,
                               substeps=100, generator=gen(0))
    m, se = tr.mean(Pe)
    exact = torch.exp(-0.7 * torch.tensor(times, dtype=torch.float64))
    assert ((m - exact).abs() <= 4 * se + 0.01).all() and m[0] == 1
    assert int(tr.n_jumps.max()) == 1                                  # one photon at most
    assert sum(len(i) for _, i, _ in tr.jumps) == int(tr.n_jumps.sum())


def fock_loss():
    """Two bosonic modes with N <= 2 particles, M_1 (+) M_2 (+) M_3, a Bose-Hubbard dimer and loss of either mode."""
    fock = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double', charges=[0, 1, 2])
    r2 = math.sqrt(2)
    a1 = {(0, 1): [[1.0, 0.0]], (1, 2): [[r2, 0.0, 0.0], [0.0, 1.0, 0.0]]}
    a2 = {(0, 1): [[0.0, 1.0]], (1, 2): [[0.0, 1.0, 0.0], [0.0, 0.0, r2]]}
    blocks = lambda a: {dc: torch.tensor(K, dtype=torch.float64) for dc, K in a.items()}
    loss = InterSectorChannel.from_blocks(fock, fock, {dc: [blocks(a1)[dc], blocks(a2)[dc]] for dc in a1})
    H = fock.from_blocks([torch.zeros(1, 1, dtype=torch.float64), torch.tensor([[0.0, 1.0], [1.0, 0.0]]),
                          torch.tensor([[0.8, r2, 0.0], [r2, 0.0, r2], [0.0, r2, 0.8]])])
    return fock, H, loss


def test_jumps_between_sectors_agree_with_the_master_equation():
    fock, H, loss = fock_loss()
    rho0 = states.random_density_matrix(fock, 1)
    times = [0.0, 0.4, 1.2]
    exact = dynamics.lindblad_evolve(rho0, H, [loss], times, rates=[0.6], method='exact')
    tr = dynamics.trajectories(rho0, H, [loss], times, rates=[0.6], n_traj=6000, substeps=60, generator=gen(1))
    for i in range(len(times)):
        rho = tr.density(i)
        assert abs(rho.Tr_blunt().real.item() - 1) < 1e-12
        p, q = rho.sector_probabilities()[0], exact[i].sector_probabilities()[0]
        se = (q * (1 - q) / tr.n_traj).sqrt()
        assert ((p - q).abs() <= 4 * se + 0.01).all(), (i, p, q)
        assert (rho.matrix - exact[i].matrix).abs().max() < 0.05
    assert tr.n_jumps.max() <= 2                                           # at most two particles to lose


def test_sector_jumps_in_a_spin_chain():
    ch = SpinChain(4, 'periodic')
    H = ch.xxz(1.0, 0.5)
    jumps = [ch.sector_jump(i, 'lower') for i in range(4)] + [ch.sz(0)]
    rho0 = ch.basis_state('1101')
    times = [0.0, 0.5, 1.5]
    exact = dynamics.lindblad_evolve(rho0, H, jumps, times, rates=[0.3] * 4 + [0.5], substeps=200)
    tr = dynamics.trajectories(rho0, H, jumps, times, rates=[0.3] * 4 + [0.5], n_traj=5000, substeps=100,
                               generator=gen(2))
    for i in range(len(times)):
        p, q = tr.density(i).sector_probabilities()[0], exact[i].sector_probabilities()[0]
        se = (q * (1 - q) / tr.n_traj).sqrt()
        assert ((p - q).abs() <= 4 * se + 0.01).all(), (i, p, q)
        m, s = tr.mean(ch.sz(1))
        assert abs(m[i] - exact[i].expectation(ch.sz(1)).real[0]) <= 4 * s[i] + 0.01


def test_mixed_initial_state_without_jumps_is_the_unitary_evolution():
    alg = TypeIAlgebra([2, 3], [2, 3], precision='double')
    torch.manual_seed(3)
    rho0 = states.random_density_matrix(alg, 1)
    H = alg.operator_from_eigenvalues(lambda d: torch.randn(1, d, dtype=torch.float64), force_self_adjoint=True)
    tr = dynamics.trajectories(rho0, H, [], [0.0, 1.0], n_traj=20000, substeps=5, generator=gen(4))
    exact = dynamics.von_neumann(rho0, H, [0.0, 1.0])
    for i in range(2):
        assert (tr.density(i).matrix - exact[i].matrix).abs().max() < 0.03
    assert int(tr.n_jumps.sum()) == 0
    norms = torch.linalg.vector_norm(tr.states, dim=(-2, -1))
    assert torch.allclose(norms, torch.ones_like(norms), atol=1e-12)


def test_layouts_give_the_same_trajectories():
    out = []
    for layout in ('padded', 'packed'):
        alg = TypeIAlgebra([1, 4, 6, 4, 1], [1, 4, 6, 4, 1], precision='double', layout=layout)
        torch.manual_seed(5)
        H = alg.operator_from_eigenvalues(lambda d: torch.randn(1, d, dtype=torch.float64), force_self_adjoint=True)
        L = alg.operator(torch.randn(1, 5, 6, 6, dtype=C128) / 3)
        psi0 = torch.zeros(50, 5, 6, dtype=C128)
        psi0[:, 2, 0] = 1
        out.append(dynamics.trajectories(psi0, H, [L], [0.0, 1.0, 2.0], substeps=40, generator=gen(6)))
    assert torch.allclose(out[0].states, out[1].states, atol=1e-10)
    assert torch.equal(out[0].n_jumps, out[1].n_jumps)


def test_errors():
    alg = TypeIAlgebra([2], [2], precision='double')
    H = alg.identity()
    with pytest.raises(ValueError):
        dynamics.trajectories(states.tracial_state(alg), H, [], [0, 1])                 # n_traj missing
    with pytest.raises(ValueError):
        dynamics.trajectories(torch.ones(1, 1, 2), lambda t: H, [], [0, 1])             # time dependent
    with pytest.raises(ValueError):
        dynamics.trajectories(torch.ones(1, 1, 2), H, [alg.identity(3)], [0, 1])       # batched jumps
