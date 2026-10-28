"""SectorJump: sparse jumps between sectors for the Lindblad solvers (agreement with the dense form, batches of rates,
and the exact statement (a) of the noisy-circuit use case of man/article)."""
import math

import pytest
import torch

from torch_vn_algebra import SpinChain, SectorJump, dynamics, states


@pytest.fixture(scope='module')
def chain():
    return SpinChain(5, 'periodic', precision='double')


def _gain(K, rho, alg):
    """The jump part sum K rho K^* of the Lindblad solvers, on a padded tensor."""
    slabs = alg._unpad(rho) if alg._packed else [rho]
    out = dynamics._add_gain([torch.zeros_like(s) for s in slabs], slabs, None, None, [(1.0, K)], alg)
    return alg._pad(out) if alg._packed else out[0]


def _gain_dual(K, Z, alg):
    slabs = alg._unpad(Z) if alg._packed else [Z]
    out = dynamics._gain_dual(slabs, None, None, [(1.0, K)], alg)
    return alg._pad(out) if alg._packed else out[0]


@pytest.mark.parametrize('kind', ['lower', 'raise', 'sz'])
def test_sandwich_dual_and_LdL_agree_with_the_dense_jump(chain, kind):
    rho = states.random_density_matrix(chain.algebra, batch_size=3).matrix
    for i in range(chain.L):
        J = chain.sector_jump(i, kind)
        K = J.to_inter_sector().kraus.to(rho.dtype)
        assert torch.allclose(J.sandwich(rho), _gain(K, rho, chain.algebra), atol=1e-13)
        assert torch.allclose(J.dual_sandwich(rho), _gain_dual(K, rho, chain.algebra), atol=1e-13)
        LdL = (K.conj().transpose(-2, -1) @ K).sum(dim=(1, 2))
        assert torch.allclose(J.LdL(padded=True).to(rho.dtype), LdL, atol=1e-13)


def test_lowering_jump_equals_the_dense_lowering(chain):
    rho = states.random_density_matrix(chain.algebra, batch_size=2).matrix
    for i in range(chain.L):
        dense = chain.lowering(i).kraus.to(rho.dtype)
        assert torch.allclose(chain.sector_jump(i, 'lower').sandwich(rho), _gain(dense, rho, chain.algebra), atol=1e-13)


@pytest.mark.parametrize('method', ['rk4', 'kraus', 'exact'])
def test_lindblad_evolution_agrees_with_dense_jumps(chain, method):
    n = chain.L
    H = chain.xxz(J=1.0, Delta=0.3)
    rho0 = chain.basis_state('11000')
    sparse = [chain.sector_jump(i, 'lower') for i in range(n)] + [chain.sector_jump(i, 'sz') for i in range(n)]
    dense = [chain.lowering(i) for i in range(n)] + [chain.sz(i) for i in range(n)]
    rates = [0.1] * n + [0.05] * n
    a = dynamics.lindblad_evolve(rho0, H, sparse, [0.0, 0.7], rates=rates, substeps=10, method=method)[-1].matrix
    b = dynamics.lindblad_evolve(rho0, H, dense, [0.0, 0.7], rates=rates, substeps=10, method=method)[-1].matrix
    assert torch.allclose(a, b, atol=1e-12)


def test_a_batch_of_coefficients_is_a_batch_of_rates(chain):
    n = chain.L
    H = chain.xxz(J=1.0, Delta=0.0)
    rho0 = chain.basis_state('11000')
    g = torch.tensor([0.0, 0.02, 0.1], dtype=torch.float64)
    jumps = [chain.sector_jump(i, 'lower').scaled(g.sqrt()) for i in range(n)]
    batched = dynamics.lindblad_evolve(rho0, H, jumps, [0.0, 1.0], substeps=10, method='kraus')[-1].matrix
    for b, gb in enumerate(g.tolist()):
        single = dynamics.lindblad_evolve(rho0, H, [chain.sector_jump(i, 'lower') for i in range(n)], [0.0, 1.0],
                                          rates=[gb] * n, substeps=10, method='kraus')[-1].matrix
        assert torch.allclose(batched[b], single[0], atol=1e-13)


def test_post_selection_removes_zero_temperature_damping_exactly():
    n, w, T, g = 6, 3, 2.0, 0.07
    q = SpinChain(n, 'periodic', precision='double')
    H = q.xxz(J=1.0, Delta=0.0)
    rho0 = q.basis_state('111000')
    ideal = dynamics.von_neumann(rho0, H, [T])[0]
    rho = dynamics.lindblad_evolve(rho0, H, [q.sector_jump(i, 'lower') for i in range(n)], [0.0, T],
                                   rates=[g] * n, substeps=40, method='kraus')[-1]
    P = q.algebra.central([1.0 if N == w else 0.0 for N in range(n + 1)])
    post, p = rho.condition_on(P)
    assert states.trace_of_product(post, ideal).real.item() == pytest.approx(1.0, abs=1e-12)
    assert p.item() == pytest.approx(math.exp(-w * g * T), rel=1e-4)


def test_constructor_rejects_non_injective_maps(chain):
    alg = chain.algebra
    C, k = alg.C, alg.k_max
    target = torch.arange(C)
    dst = torch.zeros(C, k, dtype=torch.long)
    val = torch.zeros(C, k)
    val[2, :2] = 1.0                                # two basis vectors of sector 2 onto the same one
    with pytest.raises(ValueError):
        SectorJump(alg, target, dst, val.to(alg.dtype))


@pytest.mark.parametrize('method', ['rk4', 'kraus', 'exact'])
def test_gradients_through_the_evolution_match_finite_differences(method):
    q = SpinChain(4, 'periodic', precision='double')
    H = q.xxz(J=1.0, Delta=0.0)
    rho0 = q.basis_state('1100')
    target = dynamics.von_neumann(rho0, H, [1.0])[0]

    def fidelity(J, g):
        Hm = q.algebra.operator(H.matrix * J, is_self_adjoint=True)
        jumps = [q.sector_jump(i, k).scaled(g.sqrt().reshape(1)) for k in ('lower', 'raise', 'sz') for i in range(4)]
        rho = dynamics.lindblad_evolve(rho0, Hm, jumps, [0.0, 1.0], substeps=20, method=method)[-1]
        return states.trace_of_product(rho, target).real.sum()

    J = torch.tensor(0.7, dtype=torch.float64, requires_grad=True)
    g = torch.tensor(0.05, dtype=torch.float64, requires_grad=True)
    gJ, gg = torch.autograd.grad(fidelity(J, g), (J, g))
    t = lambda x: torch.tensor(x, dtype=torch.float64)
    with torch.no_grad():
        fdJ = (fidelity(t(0.7 + 1e-6), t(0.05)) - fidelity(t(0.7 - 1e-6), t(0.05))) / 2e-6
        fdg = (fidelity(t(0.7), t(0.05 + 1e-6)) - fidelity(t(0.7), t(0.05 - 1e-6))) / 2e-6
    assert torch.isfinite(gJ) and torch.isfinite(gg)
    assert gJ.item() == pytest.approx(fdJ.item(), rel=1e-6)
    assert gg.item() == pytest.approx(fdg.item(), rel=1e-6)
