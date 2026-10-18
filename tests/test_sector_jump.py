"""SectorJump: sparse jumps between sectors for the Lindblad solvers (agreement with the dense form, batches of rates,
and the exact statement (a) of the noisy-circuit use case of man/article)."""
import math

import pytest
import torch

from torch_vn_algebra import SpinChain, SectorJump, dynamics, states


@pytest.fixture(scope='module')
def chain():
    return SpinChain(5, 'periodic', precision='double')


def _gain(K, rho):
    return dynamics._add_gain(torch.zeros_like(rho), rho, None, None, [(1.0, K)])


@pytest.mark.parametrize('kind', ['lower', 'raise', 'sz'])
def test_sandwich_dual_and_LdL_agree_with_the_dense_jump(chain, kind):
    rho = states.random_density_matrix(chain.algebra, batch_size=3).matrix
    for i in range(chain.L):
        J = chain.sector_jump(i, kind)
        K = J.to_inter_sector().kraus.to(rho.dtype)
        assert torch.allclose(J.sandwich(rho), _gain(K, rho), atol=1e-13)
        assert torch.allclose(J.dual_sandwich(rho), dynamics._gain_dual(rho, None, None, [(1.0, K)]), atol=1e-13)
        LdL = (K.conj().transpose(-2, -1) @ K).sum(dim=(1, 2))
        assert torch.allclose(J.LdL().to(rho.dtype), LdL, atol=1e-13)


def test_lowering_jump_equals_the_dense_lowering(chain):
    rho = states.random_density_matrix(chain.algebra, batch_size=2).matrix
    for i in range(chain.L):
        dense = chain.lowering(i).kraus.to(rho.dtype)
        assert torch.allclose(chain.sector_jump(i, 'lower').sandwich(rho), _gain(dense, rho), atol=1e-13)


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
