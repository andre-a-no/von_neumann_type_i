"""Regression tests for representation boundaries, checkpoint loading and spectral differentiation."""
import inspect
import math

import pytest
import torch

from torch_vn_algebra import (
    TypeIAlgebra, DensityMatrix, InterSectorChannel, dynamics, divergences, interop, learn, states,
)


@pytest.fixture
def qubit():
    torch.manual_seed(103)
    return TypeIAlgebra([2], [2], precision='double', device='cpu')


@pytest.mark.parametrize('layer', [learn.ChannelLayer, learn.LindbladLayer])
@pytest.mark.parametrize('assign', [False, True])
def test_assigned_checkpoint_controls_the_forward_pass(qubit, layer, assign):
    rho = states.random_density_matrix(qubit, 4)
    model, saved = layer(qubit), layer(qubit)
    expected = saved(rho).matrix.detach()
    supports_assign = 'assign' in inspect.signature(model.load_state_dict).parameters
    if assign and not supports_assign:
        pytest.skip('assign=True requires PyTorch 2.1 or newer')
    model.load_state_dict(saved.state_dict(), **({'assign': assign} if supports_assign else {}))
    assert all(torch.equal(v, saved.state_dict()[k]) for k, v in model.state_dict().items())
    result = model(rho).matrix
    torch.testing.assert_close(result, expected)
    result.abs().square().sum().backward()
    assert all(p.grad is not None for p in model.parameters())


@pytest.mark.parametrize('split_inputs', [False, True])
def test_dense_inter_sector_jumps_preserve_the_lindblad_generator(split_inputs):
    alg = TypeIAlgebra([1, 1], [1, 1], precision='double', device='cpu')
    K = torch.zeros(1, 2 if split_inputs else 1, 2, 2, 1, 1, dtype=alg.dtype)
    K[0, 0, 0, 0, 0, 0] = 1
    K[0, 1 if split_inputs else 0, 0, 1, 0, 0] = 1
    jump = InterSectorChannel(alg, alg, K)
    dense_K = interop.channel_to_dense(jump)
    gram = (dense_K.mH @ dense_K).sum(1)
    torch.testing.assert_close(gram, torch.eye(2, dtype=alg.dtype)[None])
    rho = DensityMatrix(alg, matrix=torch.tensor([.3, .7], dtype=alg.dtype).reshape(1, 2, 1, 1))
    expected = interop.to_dense(alg.operator(dynamics.lindblad_rhs(None, [jump])(0., rho.matrix)))
    dense_rho = interop.to_dense(rho)
    actual = torch.zeros_like(dense_rho)
    for L in interop.channel_to_dense(jump)[0]:
        gram = L.mH @ L
        actual = actual + L @ dense_rho @ L.mH - .5 * (gram @ dense_rho + dense_rho @ gram)
    torch.testing.assert_close(actual, expected)


def test_coupled_spectral_function_has_its_full_diagonal_jacobian(qubit):
    X = torch.diag(torch.tensor([1., 2.], dtype=qubit.dtype)).reshape(1, 1, 2, 2).requires_grad_()
    F = qubit.operator(X, is_self_adjoint=True).apply_function(
        lambda w: w / w.sum(-1, keepdim=True), per_block=True).matrix
    actual, = torch.autograd.grad(F[0, 0, 0, 0].real, X)
    expected = torch.diag(torch.tensor([2 / 9, -1 / 9], dtype=qubit.dtype)).reshape_as(X)
    torch.testing.assert_close(actual, expected)


@pytest.mark.parametrize('kind', ['petz', 'sandwiched'])
@pytest.mark.parametrize('precision', ['single', 'double'])
def test_renyi_alpha_rounding_to_one_uses_the_entropy_limit(kind, precision):
    alg = TypeIAlgebra([1, 1], [1, 1], precision=precision, device='cpu')
    rho = DensityMatrix(alg, matrix=torch.tensor([.9, .1], dtype=alg.dtype).reshape(1, 2, 1, 1))
    sigma = DensityMatrix(alg, matrix=torch.tensor([.5, .5], dtype=alg.dtype).reshape(1, 2, 1, 1))
    # The exact correction to the KL limit at alpha=1+1e-8 is much smaller than single-precision accuracy.
    expected = .9 * math.log(1.8) + .1 * math.log(.2)
    actual = divergences.renyi_divergence(rho, sigma, 1.00000001, kind=kind)
    assert actual.item() == pytest.approx(expected, abs=1e-5)


def test_real_interop_rejects_unrepresentable_imaginary_coherences():
    alg = TypeIAlgebra([2], [2], precision='double', complex_valued=False, device='cpu')
    pure = torch.tensor([[.5, -.5j], [.5j, .5]], dtype=torch.complex128)
    with pytest.raises(ValueError, match='complex|imaginary|real'):
        interop.from_dense(pure, alg, state=True)


def test_pending_lazy_expression_drops_flags_after_parent_mutation(qubit):
    A = qubit.operator(torch.eye(2, dtype=qubit.dtype)[None, None], is_self_adjoint=True, is_positive=True)
    B = A + A
    A.matrix.copy_(torch.tensor([[[[0., 1.], [0., 0.]]]], dtype=qubit.dtype))
    torch.testing.assert_close(B.matrix, 2 * A.matrix)
    assert not B.is_self_adjoint
    torch.testing.assert_close(B.operator_norm(), torch.tensor([2.], dtype=torch.float64))


def test_exported_jumps_match_the_classical_solution_in_qutip():
    qt = pytest.importorskip('qutip')
    alg = TypeIAlgebra([1, 1], [1, 1], precision='double', device='cpu')
    K = torch.zeros(1, 1, 2, 2, 1, 1, dtype=alg.dtype)
    K[0, 0, 0, :, 0, 0] = 1
    jump = InterSectorChannel(alg, alg, K)
    rho = DensityMatrix(alg, matrix=torch.tensor([.3, .7], dtype=alg.dtype).reshape(1, 2, 1, 1))
    initial = interop.to_qutip(rho)
    result = qt.mesolve(0 * initial, initial, [0., .1], c_ops=interop.jumps_to_qutip([jump]),
                        options={'atol': 1e-12, 'rtol': 1e-10}).states[-1]
    expected = torch.diag(torch.tensor([1 - .7 * math.exp(-.1), .7 * math.exp(-.1)], dtype=alg.dtype))
    torch.testing.assert_close(torch.as_tensor(result.full()), expected, atol=1e-9, rtol=1e-9)
    # The exported trajectory must remain in the algebra and be importable again.
    interop.from_qutip(result, alg, state=True)


@pytest.mark.parametrize('complex_valued', [False, True])
@pytest.mark.parametrize('degenerate', [False, True])
def test_coupled_softmax_matches_matrix_exp_and_its_gradients(complex_valued, degenerate):
    alg = TypeIAlgebra([2], [2], precision='double', complex_valued=complex_valued, device='cpu')
    X = (torch.eye(2, dtype=alg.dtype) if degenerate else
         torch.tensor([[2., .3], [.3, 1.]], dtype=alg.dtype))[None, None].requires_grad_()
    t = torch.tensor(.7, dtype=torch.float64, requires_grad=True)
    actual = alg.operator(X, is_self_adjoint=True).apply_function(
        lambda w: torch.softmax(t * w, dim=-1), per_block=True).matrix
    H = .5 * (X + X.mH)
    E = torch.linalg.matrix_exp(t * H)
    expected = E / E.diagonal(dim1=-2, dim2=-1).sum(-1)[..., None, None]
    G = torch.tensor([[1., .2], [.2, -.4]], dtype=alg.dtype)[None, None]
    ga = torch.autograd.grad((actual * G).real.sum(), (X, t), retain_graph=True)
    ge = torch.autograd.grad((expected * G).real.sum(), (X, t))
    torch.testing.assert_close(actual, expected)
    for a, e in zip(ga, ge):
        torch.testing.assert_close(a, e, atol=1e-9, rtol=1e-9)


def test_coupled_spectral_second_derivative(qubit):
    X = torch.tensor([[[[2., .3], [.3, 1.]]]], dtype=qubit.dtype).requires_grad_()
    t = torch.tensor(.7, dtype=torch.float64, requires_grad=True)

    def normalized(m, scale):
        return qubit.operator(m, is_self_adjoint=True).apply_function(
            lambda w: torch.softmax(scale * w, dim=-1), per_block=True).matrix

    assert torch.autograd.gradgradcheck(normalized, (X, t), atol=1e-5)


@pytest.mark.parametrize('kind', ['petz', 'sandwiched'])
def test_renyi_alpha_derivative_at_one_is_half_the_relative_entropy_variance(kind):
    alg = TypeIAlgebra([1, 1], [1, 1], precision='double', device='cpu')
    p = torch.tensor([.9, .1], dtype=torch.float64)
    q = torch.tensor([.5, .5], dtype=torch.float64)
    rho = DensityMatrix(alg, matrix=p.to(alg.dtype).reshape(1, 2, 1, 1))
    sigma = DensityMatrix(alg, matrix=q.to(alg.dtype).reshape(1, 2, 1, 1))
    alpha = torch.tensor(1., dtype=torch.float64, requires_grad=True)
    D = divergences.renyi_divergence(rho, sigma, alpha, kind=kind)
    derivative, = torch.autograd.grad(D.sum(), alpha)
    log_ratio = torch.log(p / q)
    mean = (p * log_ratio).sum()
    variance = (p * (log_ratio - mean).square()).sum()
    torch.testing.assert_close(D[0], mean)
    torch.testing.assert_close(derivative, .5 * variance)


@pytest.mark.parametrize('kind', ['petz', 'sandwiched'])
def test_near_one_sector_contributions_are_consistent(kind):
    alg = TypeIAlgebra([2, 2], [2, 2], precision='single', device='cpu')
    rho, sigma = states.random_density_matrix(alg), states.random_density_matrix(alg)
    alpha = 1.00000001
    parts = divergences.sector_contributions(rho, sigma, alpha, kind=kind)
    reconstructed = (parts.p * (torch.log(parts.p / parts.q) + parts.conditional)).sum(-1)
    torch.testing.assert_close(reconstructed, divergences.renyi_divergence(rho, sigma, alpha, kind=kind),
                               atol=1e-5, rtol=1e-5)


def test_imported_real_channels_check_the_field_and_complex_states_keep_purity(qubit):
    real = qubit.like(complex_valued=False)
    K = torch.diag(torch.tensor([1., 1j], dtype=qubit.dtype))[None]
    with pytest.raises(ValueError, match='complex'):
        interop.channel_from_dense(K, real)
    pure = torch.tensor([[.5, -.5j], [.5j, .5]], dtype=qubit.dtype)
    converted = interop.from_dense(pure, qubit, state=True)
    torch.testing.assert_close(states.purity(converted), torch.ones(1, dtype=torch.float64))


def test_parent_versions_propagate_transitively_and_materialized_results_are_snapshots(qubit):
    A = qubit.operator(torch.eye(2, dtype=qubit.dtype)[None, None], is_self_adjoint=True)
    B = A + A
    C = B + B
    snapshot_op = A + A
    snapshot = snapshot_op.matrix.clone()
    A.matrix.copy_(torch.tensor([[[[0., 1.], [0., 0.]]]], dtype=qubit.dtype))
    assert not C.is_self_adjoint
    torch.testing.assert_close(C.operator_norm(), torch.tensor([4.], dtype=torch.float64))
    torch.testing.assert_close(snapshot, 2 * torch.eye(2, dtype=qubit.dtype)[None, None])
    torch.testing.assert_close(snapshot_op.matrix, snapshot)
