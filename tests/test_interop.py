"""interop: dense matrices, QuTiP and dynamiqs (conventions of the chain layout, maps between algebras, the
Lindblad solvers against qutip.mesolve and dynamiqs.mesolve)."""
import numpy as np
import pytest
import torch

from torch_vn_algebra import TypeIAlgebra, InterSectorChannel, SpinChain, dynamics, interop, states, channels


@pytest.fixture(scope='module')
def chain():
    return SpinChain(4, 'periodic', precision='double')


def _noisy(chain):
    L = chain.L
    H = chain.xxz(J=1.0, Delta=0.3)
    rho0 = chain.basis_state('1100')
    jumps = [chain.sector_jump(i, 'lower') for i in range(L)] + [chain.sz(i) for i in range(L)]
    rates = [0.1] * L + [0.05] * L
    rho = dynamics.lindblad_evolve(rho0, H, jumps, [0.0, 1.5], rates=rates, substeps=60, method='exact')[-1]
    return H, rho0, jumps, rates, rho


def test_dense_round_trip_with_unequal_blocks():
    alg = TypeIAlgebra([1, 3, 2], [1, 3, 2], precision='double')
    A = states.random_density_matrix(alg, batch_size=3)
    D = interop.to_dense(A)
    assert D.shape == (3, 6, 6)
    assert torch.allclose(torch.diagonal(D, dim1=-2, dim2=-1).sum(-1).real, torch.ones(3, dtype=torch.float64))
    B = interop.from_dense(D, alg, state=True)
    assert torch.allclose(B.matrix, A.matrix)


def test_from_dense_rejects_elements_outside_the_algebra():
    alg = TypeIAlgebra([2, 2], [2, 2], precision='double')
    M = torch.zeros(4, 4, dtype=torch.complex128)
    M[0, 3] = 1.0
    with pytest.raises(ValueError):
        interop.from_dense(M, alg)


def test_kraus_operators_between_different_algebras_act_like_the_channel():
    a_in = TypeIAlgebra([2, 3], [2, 3], precision='double')
    a_out = TypeIAlgebra([1, 2, 2], [1, 2, 2], precision='double')
    ch = channels.random_inter_sector_channel(a_in, a_out, kraus_rank=2, batch_size=2)
    rho = states.random_density_matrix(a_in, batch_size=2)
    K = interop.channel_to_dense(ch)
    r = interop.to_dense(rho)
    dense = (K @ r[:, None] @ K.conj().transpose(-2, -1)).sum(1)
    assert torch.allclose(dense, interop.to_dense(ch(rho)), atol=1e-12)


def test_chain_layout_matches_qutip_conventions(chain):
    qt = pytest.importorskip('qutip')
    L = chain.L

    def site(op, i):
        return qt.tensor([op if j == i else qt.qeye(2) for j in range(L)])
    for i in range(L):
        assert (interop.to_qutip(chain.sz(i), layout=chain) - site(qt.sigmaz() / 2, i)).norm() < 1e-14
    lower, up = interop.jumps_to_qutip([chain.sector_jump(1, 'lower'), chain.sector_jump(2, 'raise')], layout=chain)
    assert (lower - site(qt.sigmam(), 1)).norm() < 1e-14
    assert (up - site(qt.sigmap(), 2)).norm() < 1e-14
    ket = qt.tensor(qt.basis(2, 0), qt.basis(2, 0), qt.basis(2, 1), qt.basis(2, 1))
    assert (interop.to_qutip(chain.basis_state('1100'), layout=chain) - qt.ket2dm(ket)).norm() < 1e-14
    rho = interop.from_qutip(ket, chain.algebra, layout=chain)
    assert torch.allclose(rho.matrix, chain.basis_state('1100').matrix)


def test_lindblad_evolution_agrees_with_qutip_mesolve(chain):
    qt = pytest.importorskip('qutip')
    H, rho0, jumps, rates, rho = _noisy(chain)
    res = qt.mesolve(interop.to_qutip(H, layout=chain), interop.to_qutip(rho0, layout=chain), [0.0, 1.5],
                     c_ops=interop.jumps_to_qutip(jumps, layout=chain, rates=rates),
                     options={'atol': 1e-12, 'rtol': 1e-10})
    assert (interop.to_qutip(rho, layout=chain) - res.states[-1]).norm() < 1e-9
    back = interop.from_qutip(res.states[-1], chain.algebra, layout=chain, state=True)
    assert torch.allclose(back.matrix, rho.matrix, atol=1e-8)


def test_lindblad_evolution_agrees_with_dynamiqs(chain):
    dq = pytest.importorskip('dynamiqs')
    jax = pytest.importorskip('jax')
    jax.config.update('jax_enable_x64', True)
    H, rho0, jumps, rates, rho = _noisy(chain)
    res = dq.mesolve(interop.to_dynamiqs(H, layout=chain)[0], interop.jumps_to_dynamiqs(jumps, layout=chain, rates=rates),
                     interop.to_dynamiqs(rho0, layout=chain)[0], np.array([0.0, 1.5]),
                     method=dq.method.Tsit5(rtol=1e-10, atol=1e-12), progress_meter=False)
    out = np.asarray(res.states[-1].to_jax())
    assert np.abs(out - interop.to_dense(rho, layout=chain)[0].numpy()).max() < 1e-9


def test_inter_sector_jumps_agree_with_qutip_mesolve(chain):
    qt = pytest.importorskip('qutip')
    L = chain.L
    H = chain.xxz(J=1.0, Delta=0.0)
    rho0 = chain.basis_state('1010')
    jumps = [chain.lowering(i) for i in range(L)]                   # dense InterSectorChannel jumps
    rho = dynamics.lindblad_evolve(rho0, H, jumps, [0.0, 1.0], rates=[0.2] * L, substeps=60, method='exact')[-1]
    res = qt.mesolve(interop.to_qutip(H, layout=chain), interop.to_qutip(rho0, layout=chain), [0.0, 1.0],
                     c_ops=interop.jumps_to_qutip(jumps, layout=chain, rates=[0.2] * L),
                     options={'atol': 1e-12, 'rtol': 1e-10})
    assert (interop.to_qutip(rho, layout=chain) - res.states[-1]).norm() < 1e-9


def test_kraus_splitting_converges_to_qutip_mesolve_at_second_order(chain):
    qt = pytest.importorskip('qutip')
    H, rho0, jumps, rates, _ = _noisy(chain)
    ref = qt.mesolve(interop.to_qutip(H, layout=chain), interop.to_qutip(rho0, layout=chain), [0.0, 1.5],
                     c_ops=interop.jumps_to_qutip(jumps, layout=chain, rates=rates),
                     options={'atol': 1e-12, 'rtol': 1e-10}).states[-1]
    err = []
    for substeps in (60, 240):
        rho = dynamics.lindblad_evolve(rho0, H, jumps, [0.0, 1.5], rates=rates, substeps=substeps, method='kraus')[-1]
        err.append((interop.to_qutip(rho, layout=chain) - ref).norm())
    assert err[0] < 2e-5
    assert 10 < err[0] / err[1] < 22                         # 4 times smaller steps: 16 times smaller error


def test_channels_round_trip_through_dense_kraus_operators(chain):
    qt = pytest.importorskip('qutip')
    alg = chain.algebra
    ch = channels.random_channel(alg, kraus_rank=2, batch_size=1)
    rho = states.random_density_matrix(alg, batch_size=1)
    K = interop.channel_to_dense(ch, layout_in=chain, layout_out=chain)
    back = interop.channel_from_qutip([qt.Qobj(k.numpy(), dims=[[2] * chain.L] * 2) for k in K[0]], alg,
                                      layout_in=chain, layout_out=chain)
    assert torch.allclose(back(rho).matrix, ch(rho).matrix, atol=1e-12)
    # a dense channel that does not respect the sectors is restricted and pinched onto them
    X = qt.tensor([qt.sigmax()] + [qt.qeye(2)] * (chain.L - 1))
    flip = interop.channel_from_qutip([X], alg, layout_in=chain, layout_out=chain)
    XrX = torch.as_tensor(X.full()) @ interop.to_dense(rho, layout=chain)[0] @ torch.as_tensor(X.full())
    pinched = interop.from_dense(XrX, alg, layout=chain, atol=float('inf'))      # its blocks: E_M(X rho X)
    assert torch.allclose(flip(rho).matrix, pinched.matrix, atol=1e-12)
    with pytest.raises(TypeError):
        interop.jumps_to_qutip([ch], layout=chain)


def test_states_come_back_from_dynamiqs(chain):
    dq = pytest.importorskip('dynamiqs')
    rho = states.random_density_matrix(chain.algebra, batch_size=2)
    back = interop.from_dynamiqs(interop.to_dynamiqs(rho, layout=chain), chain.algebra, layout=chain, state=True)
    assert torch.allclose(back.matrix, rho.matrix)
