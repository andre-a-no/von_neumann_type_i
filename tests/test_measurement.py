"""measurement: POVMs (checks, Born probabilities, sampling, likelihood), instruments (post-measurement states, the
average channel, the outcome register in the centre) and the classical Fisher information."""
import math

import pytest
import torch

from torch_vn_algebra import TypeIAlgebra, SpinChain, POVM, Instrument, InterSectorChannel, channels, interop, states
from torch_vn_algebra.measurement import fisher_information


@pytest.fixture(scope='module')
def alg():
    return TypeIAlgebra([1, 3, 2], [1, 3, 2], precision='double')


def _qubit():
    return TypeIAlgebra([2], [2], precision='double')


def _bloch(alg, theta, phi=None):
    """rho = (1 + cos(theta) sigma_z + sin(theta) (cos(phi) sigma_x + sin(phi) sigma_y)) / 2, batch = theta's."""
    theta = theta.reshape(-1)
    phi = torch.zeros_like(theta) if phi is None else phi.reshape(-1).expand_as(theta)
    c, s = torch.cos(theta), torch.sin(theta)
    off = 0.5 * s * torch.exp(1j * phi)
    m = torch.stack([torch.stack([0.5 * (1 + c) + 0j, off.conj()], -1),
                     torch.stack([off, 0.5 * (1 - c) + 0j], -1)], -2)
    return states.DensityMatrix(alg, matrix=m.unsqueeze(1), validate=False)


# ----------------------------------------------------------------------
# POVM
# ----------------------------------------------------------------------
def test_check_rejects_invalid_effects(alg):
    E = POVM.basis(alg).effects.clone()
    with pytest.raises(ValueError, match="sum to the identity"):
        POVM(0.9 * E, alg)
    bad = E.clone()
    bad[0, 0] *= -1                                     # negative effect, completed so that the sum stays 1
    bad = torch.cat([bad, -2 * bad[:, :1]], dim=1)
    with pytest.raises(ValueError, match="positive"):
        POVM(bad, alg)
    X = E.clone()
    X[0, 1, 1, 0, 1] = 0.1                             # not self-adjoint
    with pytest.raises(ValueError, match="self-adjoint"):
        POVM(X, alg)
    POVM(E, alg)                                        # the valid one passes


def test_probabilities_sum_to_one_and_match_born(alg):
    rho = states.random_density_matrix(alg, batch_size=4)
    U = alg.random_unitary_operator(4)
    povm = POVM.from_unitary(U)
    POVM(povm.effects, alg)                             # a valid POVM
    p = povm.probabilities(rho)
    assert p.shape == (4, len(povm)) == (4, 6)
    assert torch.allclose(p.sum(-1), torch.ones(4, dtype=torch.float64), atol=1e-12)
    assert torch.allclose(p, states.born_probabilities(rho, [povm[j] for j in range(len(povm))]), atol=1e-12)
    from_ops = POVM([povm[j] for j in range(len(povm))])
    assert torch.allclose(from_ops.probabilities(rho), p, atol=1e-12)


def test_sector_povm_measures_the_charge(alg):
    rho = states.random_density_matrix(alg, batch_size=3)
    povm = POVM.sectors(alg)
    assert povm.outcomes == alg.charges
    assert torch.allclose(povm.probabilities(rho), states.sector_probabilities(rho), atol=1e-12)
    assert povm.is_projective()


def test_basis_povm_on_a_spin_chain_matches_qutip():
    chain = SpinChain(4, precision='double')
    povm = POVM.basis(chain.algebra)
    assert len(povm) == 16
    rho = chain.basis_state('1100')
    p = povm.probabilities(rho)[0]
    label = 0b0011                                      # bit i = spin on site i
    j = povm.outcomes.index((2, int(chain.index_of[label])))
    assert abs(float(p[j]) - 1) < 1e-14 and float(p.sum()) == pytest.approx(1)
    pytest.importorskip('qutip')
    sigma = states.random_density_matrix(chain.algebra, batch_size=1)
    dense = interop.to_qutip(sigma, layout=chain).full()
    pq = torch.as_tensor(dense.diagonal().real.copy())        # computational-basis probabilities
    pos = interop.chain_basis(chain)
    assert torch.allclose(povm.probabilities(sigma)[0], pq[pos], atol=1e-12)
    assert int(pos[j]) == 0b0011                        # |1100> = |0>|0>|1>|1> in QuTiP (|0> = spin up)
    pr = interop.to_qutip(rho, layout=chain).full().diagonal().real
    assert pr[3] == pytest.approx(1) and pr.sum() == pytest.approx(1)


def test_sampling_statistics_and_reproducibility(alg):
    rho = states.random_density_matrix(alg, batch_size=3)
    povm = POVM.basis(alg)
    shots = 200000
    g = torch.Generator().manual_seed(7)
    n = povm.sample(rho, shots, generator=g)
    assert n.dtype == torch.long and n.shape == (3, 6)
    assert bool((n.sum(-1) == shots).all())
    p = povm.probabilities(rho)
    sigma = (p * (1 - p) / shots).sqrt()
    assert bool(((n / shots - p).abs() <= 5 * sigma + 1e-12).all())
    n2 = povm.sample(rho, shots, generator=torch.Generator().manual_seed(7))
    assert torch.equal(n, n2)
    assert povm.sample(rho, 0).sum() == 0


def test_log_likelihood_is_maximal_at_the_true_state(alg):
    rho = states.random_density_matrix(alg, batch_size=1)
    povm = POVM.combine([POVM.basis(alg), POVM.from_unitary(alg.random_unitary_operator())], [0.3, 0.7])
    assert povm.outcomes[0] == (0, (0, 0)) and len(povm) == 12
    p = povm.probabilities(rho)
    counts = 1e5 * p                                    # expected counts: the maximum is exactly at rho
    best = povm.log_likelihood(rho, counts)
    for _ in range(5):
        sigma = rho.mix(states.random_density_matrix(alg), 0.05)
        assert float(povm.log_likelihood(sigma, counts)) < float(best)
    # differentiable: the gradient at the true state vanishes along the state space (only the trace direction)
    m = rho.matrix.clone().requires_grad_(True)
    ll = povm.log_likelihood(states.DensityMatrix(alg, matrix=m, validate=False), counts.reshape(-1))
    ll.sum().backward()
    assert m.grad is not None and torch.isfinite(m.grad).all()


# ----------------------------------------------------------------------
# Instrument
# ----------------------------------------------------------------------
def test_lueders_instrument_matches_lueders_update(alg):
    rho = states.random_density_matrix(alg, batch_size=3)
    povm = POVM.from_unitary(alg.random_unitary_operator())
    inst = Instrument.lueders(povm)
    p = inst.probabilities(rho)
    assert torch.allclose(p.sum(-1), torch.ones(3, dtype=torch.float64), atol=1e-12)
    post, q = inst.post_states(rho)
    assert torch.allclose(q, p, atol=1e-12)
    for j, s in enumerate(post):
        ref, pr = states.lueders_update(rho, povm[j])
        assert torch.allclose(s.matrix, ref.matrix, atol=1e-10)
        assert torch.allclose(pr, p[:, j], atol=1e-12)
        live = pr > 1e-12
        assert torch.allclose(s.trace.real[live], torch.ones_like(pr[live]), atol=1e-10)
    # Instrument.from_povm agrees for projections, and Lueders rejects non-projective effects
    post2, _ = Instrument.from_povm(povm).post_states(rho)
    assert all(torch.allclose(a.matrix, b.matrix, atol=1e-8) for a, b in zip(post, post2))
    with pytest.raises(ValueError):
        Instrument.lueders(POVM.combine([POVM.basis(alg), povm]))


def test_sector_instrument_zero_probability_outcomes(alg):
    rho = states.DensityMatrix.from_density(alg, torch.zeros(1, 3, 3, 3, dtype=torch.complex128).index_put_(
        (torch.tensor([0]), torch.tensor([1]), torch.tensor([0]), torch.tensor([0])),
        torch.tensor([1.0 + 0j], dtype=torch.complex128)))
    post, p = Instrument.lueders(POVM.sectors(alg)).post_states(rho)
    assert torch.allclose(p, torch.tensor([[0.0, 1.0, 0.0]], dtype=torch.float64))
    assert float(post[0].matrix.abs().amax()) == 0 and float(post[2].matrix.abs().amax()) == 0
    assert torch.allclose(post[1].matrix, rho.matrix)


def test_general_instrument_average_channel_and_register(alg):
    # instrument from a random channel: outcome j = Kraus operator j (Phi_j(rho) = K_j rho K_j^*)
    ch = channels.random_channel(alg, kraus_rank=3, batch_size=2)
    inst = Instrument([alg.operator(ch.kraus[:, j]) for j in range(3)])
    assert len(inst) == 3 and inst.batch_size == 2
    rho = states.random_density_matrix(alg, batch_size=2)
    p = inst.probabilities(rho)
    assert torch.allclose(p.sum(-1), torch.ones(2, dtype=torch.float64), atol=1e-12)
    avg = inst.average_channel()
    assert avg.is_trace_preserving()
    assert torch.allclose(avg(rho).matrix, ch(rho).matrix, atol=1e-12)
    post, q = inst.post_states(rho)
    phi = sum(q[:, j, None, None, None] * post[j].matrix for j in range(3))   # sum_j Phi_j(rho)
    assert torch.allclose(phi, ch(rho).matrix, atol=1e-12)
    for s in post:
        assert torch.allclose(s.trace.real, torch.ones(2, dtype=torch.float64), atol=1e-10)
    # one channel into (+)_j M: the register in the centre
    big = inst.as_channel()
    assert isinstance(big, InterSectorChannel) and big.is_trace_preserving()
    reg = inst.register_algebra
    assert reg.C == 3 * alg.C and reg.charges[alg.C + 1] == (1, alg.charges[1])
    out = big(rho)
    assert isinstance(out, states.DensityMatrix)
    sp = out.sector_probabilities().reshape(2, 3, alg.C)
    assert torch.allclose(sp.sum(-1), p, atol=1e-12)
    blocks = out.matrix.reshape(2, 3, alg.C, alg.k_max, alg.k_max)
    for j in range(3):
        K = ch.kraus[:, j]
        assert torch.allclose(blocks[:, j], K @ rho.matrix @ K.mH, atol=1e-12)
    # mixed specifications: a Channel per outcome, Kraus ranks padded
    split = Instrument([channels.Channel(alg, ch.kraus[:, :2]), alg.operator(ch.kraus[:, 2])])
    assert torch.allclose(split.probabilities(rho), torch.stack([p[:, 0] + p[:, 1], p[:, 2]], -1), atol=1e-12)
    with pytest.raises(ValueError, match="trace preserving"):
        Instrument([alg.operator(ch.kraus[:, 0])])


def test_from_povm_is_differentiable(alg):
    U0 = alg.random_unitary_operator().matrix
    t = torch.tensor(0.3, dtype=torch.float64, requires_grad=True)
    H = states.random_density_matrix(alg).matrix                       # a fixed self-adjoint generator
    U = torch.linalg.matrix_exp(-1j * t * H) @ U0
    povm = POVM.combine([POVM.basis(alg), POVM.from_unitary(alg.operator(U))])
    inst = Instrument.from_povm(povm)
    rho = states.random_density_matrix(alg)
    post, p = inst.post_states(rho)
    loss = sum(s.matrix.abs().pow(2).sum() for s in post) + p.pow(2).sum()
    loss.backward()
    assert t.grad is not None and torch.isfinite(t.grad)
    # the post-measurement state is sqrt(E) rho sqrt(E) / p
    E = povm[7]
    sq = E.sqrt()
    ref = (sq @ rho @ sq).matrix.detach() / p[0, 7].detach()
    assert torch.allclose(post[7].matrix.detach(), ref, atol=1e-10)


# ----------------------------------------------------------------------
# Fisher information
# ----------------------------------------------------------------------
def test_fisher_information_of_a_qubit_in_the_z_basis():
    q = _qubit()
    povm = POVM.basis(q)
    theta = torch.tensor([0.3, 1.1, 2.0, 2.9], dtype=torch.float64)
    p = povm.probabilities(_bloch(q, theta))
    assert torch.allclose(p[:, 0], (1 + torch.cos(theta)) / 2, atol=1e-14)
    # per-sample parameters: one theta per sample
    F = fisher_information(lambda t: povm.probabilities(_bloch(q, t)), [theta.clone()], per_sample=True)
    assert F.shape == (4, 1, 1)
    assert torch.allclose(F, torch.ones(4, 1, 1, dtype=torch.float64), atol=1e-10)
    # a shared scalar parameter
    F1 = fisher_information(lambda t: povm.probabilities(_bloch(q, t)), [torch.tensor(0.7, dtype=torch.float64)])
    assert F1.shape == (1, 1, 1) and abs(float(F1) - 1) < 1e-10


def test_fisher_information_is_symmetric_psd_and_handles_complex_parameters():
    q = _qubit()
    U = q.random_unitary_operator()
    povm = POVM.combine([POVM.basis(q), POVM.from_unitary(U)])
    th, ph = torch.tensor(0.8, dtype=torch.float64), torch.tensor(0.4, dtype=torch.float64)
    F = fisher_information(lambda a, b: povm.probabilities(_bloch(q, a, b)), [th, ph])
    assert F.shape == (1, 2, 2)
    assert torch.allclose(F, F.transpose(-2, -1), atol=1e-12)
    assert float(torch.linalg.eigvalsh(F).amin()) > -1e-12
    # z-basis only: phi is invisible, F = diag(1, 0)
    Fz = fisher_information(lambda a, b: POVM.basis(q).probabilities(_bloch(q, a, b)), [th, ph])
    assert torch.allclose(Fz[0], torch.tensor([[1.0, 0.0], [0.0, 0.0]], dtype=torch.float64), atol=1e-10)
    # a complex parameter z = e^{i phi} sin(theta)/2 enters as (Re z, Im z); compare with the real parametrisation
    z = (0.5 * torch.sin(th) * torch.exp(1j * ph)).to(torch.complex128)

    def from_z(zz):
        c = torch.cos(th)
        m = torch.stack([torch.stack([0.5 * (1 + c) + 0j, zz.conj()]), torch.stack([zz, 0.5 * (1 - c) + 0j])])
        return povm.probabilities(states.DensityMatrix(q, matrix=m[None, None], validate=False))
    Fc = fisher_information(from_z, [z])
    assert Fc.shape == (1, 2, 2)
    # chain rule: (Re z, Im z) = sin(th)/2 (cos ph, sin ph) at fixed theta, so F_phi = j^T Fc j
    jac = 0.5 * torch.sin(th) * torch.stack([-torch.sin(ph), torch.cos(ph)])
    assert float(jac @ Fc[0] @ jac) == pytest.approx(float(F[0, 1, 1]), rel=1e-8)


def test_fisher_information_vanishes_for_the_sector_measurement_of_a_unitary_orbit(alg):
    rho = states.random_density_matrix(alg)
    H = states.random_density_matrix(alg).matrix

    def probs(t):
        U = torch.linalg.matrix_exp(-1j * t * H)
        return POVM.sectors(alg).probabilities(states.DensityMatrix(alg, matrix=U @ rho.matrix @ U.mH, validate=False))
    F = fisher_information(probs, [torch.tensor(0.2, dtype=torch.float64)])
    assert float(F.abs().amax()) < 1e-10                # unitaries of M keep the sector weights
