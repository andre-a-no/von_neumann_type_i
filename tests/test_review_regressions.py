"""Regression tests for batched tolerances, physical channels, gradients and numerical edge cases."""
import math

import pytest
import torch

from torch_vn_algebra import TypeIAlgebra, divergences, krylov, learn, optimize, states
from torch_vn_algebra.measurement import fisher_information


@pytest.mark.parametrize('precision', ['single', 'double'])
@pytest.mark.parametrize('layout', ['batch', 'sectors'])
def test_independent_blocks_do_not_share_property_tolerances(precision, layout):
    large = 1e6 if precision == 'single' else 1e12
    sizes = [2] if layout == 'batch' else [2, 2]
    alg = TypeIAlgebra(sizes, sizes, precision=precision, device='cpu')
    small = torch.tensor([[0., 1.], [0., 0.]], dtype=alg.dtype)
    big = large * torch.eye(2, dtype=alg.dtype)
    blocks = [torch.stack([big, small])] if layout == 'batch' else [big, small]
    op = alg.from_blocks(blocks)
    assert not op.is_self_adjoint
    assert not op.is_normal
    assert not op.is_positive
    assert not op.is_projection
    reference = torch.linalg.matrix_norm(op.matrix, ord=2).amax(-1)
    torch.testing.assert_close(op.operator_norm(), reference)
    for flag in ('self_adjoint', 'normal', 'projection', 'invertible'):
        with pytest.raises(ValueError):
            alg.operator(op.matrix, validate=True, **{'is_' + flag: True})


@pytest.mark.parametrize('precision', ['single', 'double'])
def test_small_positive_and_invertible_blocks_keep_their_own_scale(precision):
    alg = TypeIAlgebra([2, 2], [2, 2], precision=precision, device='cpu')
    large = 1e6 if precision == 'single' else 1e12
    big = large * torch.eye(2, dtype=alg.dtype)
    good = alg.from_blocks([big, torch.eye(2, dtype=alg.dtype)])
    assert good.is_invertible
    bad = alg.from_blocks([big, torch.diag(torch.tensor([-1., 1.], dtype=alg.dtype))])
    assert not bad.is_positive
    with pytest.raises(RuntimeError, match='positive'):
        bad.sqrt()
    # Exercise the shared decomposition and entropy validation paths on fresh operators, too.
    with pytest.raises(RuntimeError, match='positive'):
        alg.operator(bad.matrix).sqrt()
    with pytest.raises(RuntimeError, match='positive'):
        alg.operator(bad.matrix).trace_a_log_a()
    with pytest.raises(ValueError):
        alg.operator(bad.matrix, validate=True, is_self_adjoint=True, is_positive=True)


@pytest.mark.parametrize('precision', ['single', 'double'])
@pytest.mark.parametrize('complex_valued', [False, True])
@pytest.mark.parametrize('between_sectors', [False, True])
def test_ill_conditioned_full_rank_stinespring_is_trace_preserving(precision, complex_valued, between_sectors):
    alg = TypeIAlgebra([2], [2], precision=precision, complex_valued=complex_valued, device='cpu')
    out = None
    if between_sectors:
        out = TypeIAlgebra([2, 1], [2, 1], precision=precision, complex_valued=complex_valued, device='cpu')
    net = learn.ChannelLayer(alg, 1, algebra_out=out, transitions=[[True], [False]] if out else None)
    small = 1e-4 if precision == 'single' else 1e-8
    with torch.no_grad():
        net.Z.zero_()
        net.Z[0, 0, 0, 0] = torch.diag(torch.tensor([1., small], dtype=alg.dtype))
    channel = net.channel()
    assert channel.is_trace_preserving()
    rho = states.maximally_mixed_state(alg)
    result = net(rho)
    torch.testing.assert_close(result.Tr_blunt().real, torch.ones(1, dtype=alg.hilbert.real_dtype))
    torch.testing.assert_close(result.matrix[:, 0], rho.matrix[:, 0])
    if out:
        assert torch.count_nonzero(channel.kraus[:, :, 1]) == 0


def test_rank_deficient_stinespring_completion_respects_transitions_and_padding():
    alg = TypeIAlgebra([1, 2, 3], [0, 2, 1], precision='double', device='cpu')
    out = TypeIAlgebra([3, 1], [3, 0], precision='double', device='cpu')
    allowed = [[False, True, True], [False, False, False]]
    net = learn.ChannelLayer(alg, 1, algebra_out=out, transitions=allowed)
    with torch.no_grad():
        net.Z.zero_()
    channel = net.channel()
    assert channel.is_trace_preserving()
    assert torch.count_nonzero(channel.kraus[:, :, 1]) == 0
    assert torch.count_nonzero(channel.kraus[:, :, :, 0]) == 0
    result = net(states.maximally_mixed_state(alg))
    torch.testing.assert_close(result.Tr_blunt().real, torch.ones(1, dtype=torch.float64))
    with pytest.raises(RuntimeError, match='full column rank'):
        channel.kraus.real.sum().backward()


@pytest.mark.parametrize('complex_valued', [False, True])
def test_stinespring_gradient_at_repeated_positive_singular_values(complex_valued):
    alg = TypeIAlgebra([2], [2], precision='double', complex_valued=complex_valued, device='cpu')
    param = optimize.StinespringParam(alg, 1, kraus_rank=2)
    # Every singular value is 1: a native SVD backward would be singular here.
    Z = alg.random_isometry(4, 2).reshape_as(param.Z).detach().requires_grad_(True)

    def kraus(z):
        param.Z = z
        return param.kraus()

    assert torch.autograd.gradcheck(kraus, (Z,), eps=1e-6, atol=1e-5, rtol=1e-4)


@pytest.mark.parametrize('per_sample', [False, True])
def test_fisher_gradient_includes_probability_denominator(per_sample):
    a = torch.tensor(1., dtype=torch.float64, requires_grad=True)
    theta = torch.tensor([0.7, 1.], dtype=torch.float64)

    def probabilities(t):
        p = torch.sigmoid(a * t)
        return torch.stack([p, 1 - p], dim=-1)

    value = fisher_information(probabilities, [theta], per_sample=per_sample, create_graph=True).sum()
    actual, = torch.autograd.grad(value, a)
    # Bernoulli logistic model: exact Fisher information with respect to theta is a^2 p(1-p).
    p = torch.sigmoid(a * theta)
    exact = (a ** 2 * p * (1 - p)).sum()
    expected, = torch.autograd.grad(exact, a)
    torch.testing.assert_close(value, exact)
    torch.testing.assert_close(actual, expected)


def test_in_place_mutation_invalidates_all_operator_caches():
    alg = TypeIAlgebra([2], [2], precision='double', device='cpu')
    op = alg.operator_from_eigenvalues(lambda n: torch.tensor([1., 2.], dtype=torch.float64)).cache_spectrum()
    _ = op.eigenvalues(), op.trace, op.lambda_min, op.lambda_max, op.inv.matrix, op.abs().matrix, op.sqrt().matrix
    op.matrix.mul_(4.)
    assert not op.cache_info()['cached']
    reference = torch.linalg.eigvalsh(op.matrix[:, 0])
    torch.testing.assert_close(op.eigenvalues()[0], reference)
    torch.testing.assert_close(op.lambda_min, reference[:, 0])
    torch.testing.assert_close(op.operator_norm(), reference[:, -1])
    torch.testing.assert_close(op.trace, op.Tr_blunt())
    torch.testing.assert_close(op.inv.matrix, torch.linalg.inv(op.matrix))
    torch.testing.assert_close(op.abs().matrix, op.matrix)
    # Invalidation also drops flags supplied by the user / inferred at construction.
    op.matrix.copy_(torch.tensor([[[[0., 1.], [0., 0.]]]], dtype=alg.dtype))
    assert not op.is_self_adjoint and not op.is_normal and not op.is_positive
    assert not op.is_invertible and not op.is_projection
    assert not (-op).is_self_adjoint and not op.adjoint().is_self_adjoint


@pytest.mark.parametrize('precision', ['single', 'double'])
@pytest.mark.parametrize('p', [1., 2., 10., 100., math.inf])
def test_zero_weight_sectors_do_not_set_lp_scale(precision, p):
    alg = TypeIAlgebra([1, 1], [1, 1], precision=precision, device='cpu')
    op = alg.central([1e10, 1.])
    # Independent weight vectors broadcast over one operator, including the identically zero trace.
    weights = torch.tensor([[0., 1.], [1e-20, 0.], [0., 0.]], dtype=alg.hilbert.real_dtype)
    actual = divergences.lp_norm(op, p, weights=weights)
    second = 1e10 if p == math.inf else 1e10 * 1e-20 ** (1 / p)
    torch.testing.assert_close(actual, torch.tensor([1., second, 0.], dtype=actual.dtype))


@pytest.mark.parametrize('zero', [False, True])
def test_lp_zero_weight_gradients_remain_finite(zero):
    alg = TypeIAlgebra([1, 1], [1, 1], precision='double', device='cpu')
    values = torch.tensor([1e10, 0. if zero else 1.], dtype=torch.float64, requires_grad=True)
    op = alg.operator(values.reshape(1, 2, 1, 1), is_self_adjoint=True)
    value = divergences.lp_norm(op, 100., weights=torch.tensor([0., 1.])).sum()
    grad, = torch.autograd.grad(value, values)
    torch.testing.assert_close(grad, torch.tensor([0., 0. if zero else 1.], dtype=torch.float64))


@pytest.mark.parametrize('dtype', [torch.float32, torch.complex64, torch.float64, torch.complex128])
def test_zero_lanczos_sample_and_breakdown_remain_finite(dtype):
    start = torch.tensor([[0., 0.], [1., 0.]], dtype=dtype)
    diagonal = torch.tensor([[1., 2.], [1., 2.]], dtype=start.real.dtype)
    H = krylov.SparseSectorHamiltonian(diagonal, None, dtype)
    alpha, beta, Q = krylov.lanczos(H.matvec, start, m=2)
    assert torch.isfinite(alpha).all() and torch.isfinite(beta).all() and torch.isfinite(Q).all()
    assert torch.count_nonzero(Q[0]) == 0
    assert torch.count_nonzero(beta) == 0
    times = torch.tensor([0., 0.5, 1.], dtype=torch.float64)
    actual = krylov.evolve(H, start, times)
    reference = torch.exp(-1j * times[:, None, None] * diagonal[None]) * start[None]
    torch.testing.assert_close(actual, reference.to(actual.dtype))
