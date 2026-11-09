"""Regressions for sampled-spectrum ownership and optimisation state contracts."""
import pytest
import torch

from torch_vn_algebra import TypeIAlgebra, optimize


def make_alg(sizes=(2,), complex_valued=False):
    return TypeIAlgebra(list(sizes), list(sizes), complex_valued=complex_valued,
                        precision='double', device='cpu')


@pytest.mark.parametrize('complex_valued', [False, True])
@pytest.mark.parametrize('materialized', [False, True])
def test_sampled_spectrum_is_a_snapshot_of_sampler_data(complex_valued, materialized):
    alg = make_alg(complex_valued=complex_valued)
    sample = torch.tensor([[1., 2.]], dtype=alg.dtype)
    op = alg.operator_from_eigenvalues(lambda k: sample)
    old_matrix = op.matrix.detach().clone() if materialized else None
    with torch.no_grad():
        sample.copy_(torch.tensor([[-1., 5.]], dtype=alg.dtype))

    # A materialised operator stays a snapshot. A pending draw must use the same
    # spectrum on which its already-computed property flags were based.
    dense = op.matrix[:, 0]
    if old_matrix is not None:
        torch.testing.assert_close(op.matrix, old_matrix)
    expected = torch.tensor([[1., 2.]], dtype=torch.float64)
    torch.testing.assert_close(torch.linalg.eigvalsh(dense), expected)
    torch.testing.assert_close(op.eigenvalues()[0], torch.linalg.eigvalsh(dense))
    torch.testing.assert_close(op.operator_norm(), torch.linalg.matrix_norm(dense, ord=2))
    assert op.is_positive
    square_root = op.sqrt().matrix[:, 0]
    torch.testing.assert_close(square_root @ square_root, dense)


@pytest.mark.parametrize('complex_valued', [False, True])
def test_spectrum_sampler_may_reuse_a_buffer_between_sectors(complex_valued):
    alg = make_alg((2, 2), complex_valued)
    buffer = torch.zeros(2, dtype=alg.dtype)
    calls = []

    def sample(k):
        values = torch.tensor([1., 2.], dtype=alg.dtype) + 2*len(calls)
        buffer.copy_(values)
        calls.append(k)
        return buffer

    op = alg.operator_from_eigenvalues(sample)
    assert calls == [2, 2]
    actual = torch.linalg.eigvalsh(op.matrix)
    expected = torch.tensor([[[1., 2.], [3., 4.]]], dtype=torch.float64)
    torch.testing.assert_close(actual, expected)
    for c in range(2):
        torch.testing.assert_close(op.eigenvalues()[c], actual[:, c])


@pytest.mark.parametrize('complex_valued', [False, True])
@pytest.mark.parametrize('path', ['matrix', 'known_spectrum'])
def test_sampled_spectrum_snapshot_keeps_sampler_gradients(complex_valued, path):
    alg = make_alg(complex_valued=complex_valued)
    sample = torch.tensor([[1., 3.]], dtype=alg.dtype, requires_grad=True)
    op = alg.operator_from_eigenvalues(lambda k: sample)
    values = op.matrix if path == 'matrix' else op.eigenvalues()[0]
    loss = values.abs().square().sum()
    gradient, = torch.autograd.grad(loss, sample)
    # Unitary invariance: ||U diag(sample) U*||_F^2 = sum |sample|^2.
    torch.testing.assert_close(gradient, 2*sample.detach())


@pytest.mark.parametrize('complex_valued', [False, True])
def test_unitary_cache_tracks_unversioned_inference_base(complex_valued):
    alg = make_alg(complex_valued=complex_valued)
    with torch.inference_mode():
        base = torch.eye(2, dtype=alg.dtype)[None, None].clone()
    param = optimize.UnitaryParam(alg, 1, base=base)
    with torch.no_grad():
        torch.testing.assert_close(param.matrix(), torch.eye(2, dtype=alg.dtype)[None, None])
    rotation = torch.tensor([[[[0., -1.], [1., 0.]]]], dtype=alg.dtype)
    with torch.inference_mode():
        base.copy_(rotation)
    # G is zero, so the exact result is the new base even without a tensor version counter.
    with torch.no_grad():
        torch.testing.assert_close(param.matrix(), rotation)


@pytest.mark.parametrize('complex_valued', [False, True])
def test_unitary_cache_tracks_in_place_base_changes(complex_valued):
    alg = make_alg(complex_valued=complex_valued)
    base = alg.identity(2).matrix.clone()
    param = optimize.UnitaryParam(alg, 2, base=base)
    with torch.no_grad():
        param.G[0, 0, 0, 1] = .17
        param.G[1, 0, 1, 0] = -.23
    previous = param.matrix().detach().clone()
    with torch.no_grad():
        rotation = torch.tensor([[0., -1.], [1., 0.]], dtype=alg.dtype)
        base.copy_(rotation.expand_as(base))
    skew = param.G-param.G.mH
    reference = base @ torch.linalg.matrix_exp(skew)
    current = param.matrix()
    torch.testing.assert_close(current, reference)
    assert not torch.allclose(current, previous)
    torch.testing.assert_close(current.mH @ current, alg.identity(2).matrix)


@pytest.mark.parametrize('sizes', [(2,), (2, 1)])
@pytest.mark.parametrize('complex_valued', [False, True])
@pytest.mark.parametrize('tensor_scale', [False, True])
@pytest.mark.parametrize('value', [1j, 1+0j, .3])
def test_exponential_scalar_field_contract(sizes, complex_valued, tensor_scale, value):
    alg = make_alg(sizes, complex_valued)
    blocks = [torch.diag(torch.linspace(.5, 1.5, k, dtype=torch.float64)).to(alg.dtype)
              for k in sizes]
    op = alg.from_blocks(blocks)
    scale = torch.tensor(value, dtype=torch.complex128 if isinstance(value, complex) else torch.float64) \
        if tensor_scale else value
    if isinstance(value, complex) and not complex_valued:
        with pytest.raises(ValueError, match='complex'):
            op.expm(scale).matrix
        return
    actual = op.expm(scale).matrix
    assert actual.dtype == alg.dtype
    for c, block in enumerate(blocks):
        k = sizes[c]
        torch.testing.assert_close(actual[:, c, :k, :k], torch.linalg.matrix_exp(scale*block)[None])
        if k < alg.k_max:
            assert torch.count_nonzero(actual[:, c, k:, :]) == 0
            assert torch.count_nonzero(actual[:, c, :, k:]) == 0


def test_complex_tensor_exponential_scale_keeps_its_gradient():
    alg = make_alg(complex_valued=True)
    matrix = torch.tensor([[.4, .1j], [-.1j, 1.2]], dtype=alg.dtype)
    scale = torch.tensor(.2+.3j, dtype=alg.dtype, requires_grad=True)
    actual = alg.from_blocks([matrix]).expm(scale).matrix
    reference = torch.linalg.matrix_exp(scale*matrix)
    actual_grad, = torch.autograd.grad(actual.real.sum(), scale, retain_graph=True)
    reference_grad, = torch.autograd.grad(reference.real.sum(), scale)
    torch.testing.assert_close(actual[0, 0], reference)
    torch.testing.assert_close(actual_grad, reference_grad)


@pytest.mark.parametrize('maximize', [False, True])
@pytest.mark.parametrize('invalid', [float('nan'), float('inf'), -float('inf')])
@pytest.mark.parametrize('finite_first', [False, True])
def test_extremize_ignores_nonfinite_values_when_a_finite_iterate_exists(maximize, invalid, finite_first):
    alg = make_alg()
    param = optimize.UnitaryParam(alg, 2)
    finite = torch.tensor([2., -3.], dtype=torch.float64)
    calls = []

    def objective():
        is_first = not calls
        calls.append(True)
        value = finite if is_first == finite_first else torch.full_like(finite, invalid)
        # Constant objective values with a valid zero gradient isolate best-value
        # bookkeeping from undefined derivatives at an objective singularity.
        return param.G.sum((1, 2, 3))*0+value

    best = optimize.extremize(objective, [param], maximize=maximize, rounds=1, steps=1)
    torch.testing.assert_close(best, finite)
    assert len(calls) == 2


@pytest.mark.parametrize('maximize', [False, True])
@pytest.mark.parametrize('invalid', [float('nan'), float('inf'), -float('inf')])
def test_extremize_reports_nan_when_no_finite_value_was_attained(maximize, invalid):
    alg = make_alg()
    param = optimize.UnitaryParam(alg, 1)
    objective = lambda: param.G.sum((1, 2, 3))*0+invalid
    best = optimize.extremize(objective, [param], maximize=maximize, rounds=1, steps=1)
    assert torch.isnan(best).all()
