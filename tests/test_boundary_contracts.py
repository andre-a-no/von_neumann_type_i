"""Regression tests for lazy domains, inferred properties, module transfers and checkpoints."""
import inspect
import io

import pytest
import torch
from torch import nn

from torch_vn_algebra import TypeIAlgebra, learn, measurement, states

pytestmark = pytest.mark.filterwarnings('ignore:Complex modules are a new feature')


def algebra(precision='double', complex_valued=True):
    return TypeIAlgebra([2], [2], precision=precision, complex_valued=complex_valued, device='cpu')


def basis_state(alg):
    return states.DensityMatrix(alg, matrix=torch.tensor([[[[1., 0.], [0., 0.]]]], dtype=alg.dtype))


@pytest.mark.parametrize('start_precision', ['single', 'double'])
@pytest.mark.parametrize('kind', ['channel', 'unitary', 'postselect'])
def test_module_dtype_conversion_accepts_input_in_the_target_dtype(start_precision, kind):
    torch.manual_seed(6)
    alg = algebra(start_precision)
    constructors = {'channel': learn.ChannelLayer, 'unitary': learn.UnitaryLayer,
                    'postselect': lambda a: learn.PostSelect(a, [1])}
    model = constructors[kind](alg).to(dtype=torch.complex64)
    target = alg.like(precision='single')
    result = model(basis_state(target))
    assert result.algebra.dtype == torch.complex64
    assert result.matrix.dtype == torch.complex64


@pytest.mark.parametrize('different_measurement', [False, True])
def test_feedforward_checkpoint_restores_the_fixed_measurement(different_measurement):
    alg = algebra()
    z_measurement = measurement.POVM.basis(alg)
    H = alg.operator(torch.tensor([[[[1., 1.], [1., -1.]]]], dtype=alg.dtype) / 2 ** .5)
    other = measurement.POVM.from_unitary(H) if different_measurement else z_measurement
    saved = learn.FeedForward(measurement.Instrument.lueders(z_measurement), [nn.Identity(), nn.Identity()])
    restored = learn.FeedForward(measurement.Instrument.lueders(other), [nn.Identity(), nn.Identity()])
    restored.load_state_dict(saved.state_dict())
    rho = basis_state(alg)
    # A nonselective Z measurement leaves |0><0| unchanged.
    torch.testing.assert_close(restored(rho).matrix, rho.matrix)


@pytest.mark.parametrize('coefficient', ['python', 'tensor_scalar', 'tensor_batch'])
@pytest.mark.parametrize('precision', ['single', 'double'])
@pytest.mark.parametrize('complex_valued', [False, True])
def test_zero_multiple_has_the_functional_calculus_of_zero(coefficient, precision, complex_valued):
    alg = algebra(precision, complex_valued)
    A = alg.operator(torch.tensor([[[[0., 1.], [0., 0.]]]], dtype=alg.dtype))
    assert not A.is_self_adjoint and not A.is_normal
    scale = {'python': 0., 'tensor_scalar': torch.tensor(0.), 'tensor_batch': torch.zeros(2)}[coefficient]
    zero = A * scale
    actual = zero.apply_function(torch.exp).matrix
    torch.testing.assert_close(actual, torch.eye(2, dtype=alg.dtype)[None, None].expand_as(actual))


@pytest.mark.parametrize('cutoff', [0., 2.])
@pytest.mark.parametrize('precision', ['single', 'double'])
@pytest.mark.parametrize('complex_valued', [False, True])
def test_pseudoinverse_invertibility_matches_its_actual_matrix(cutoff, precision, complex_valued):
    alg = algebra(precision, complex_valued)
    result = alg.identity().inverse(tol=cutoff)
    expected = torch.eye(2, dtype=alg.dtype) if cutoff < 1 else torch.zeros(2, 2, dtype=alg.dtype)
    torch.testing.assert_close(result.matrix, expected[None, None])
    assert result.is_invertible == (cutoff < 1)


@pytest.mark.parametrize('cutoff', [0., .5])
@pytest.mark.parametrize('precision', ['single', 'double'])
@pytest.mark.parametrize('complex_valued', [False, True])
def test_pseudoinverse_positivity_matches_the_retained_spectrum(cutoff, precision, complex_valued):
    alg = algebra(precision, complex_valued)
    A = alg.operator(torch.diag(torch.tensor([-.1, 2.], dtype=alg.dtype))[None, None])
    assert not A.is_positive
    result = A.inverse(tol=cutoff)
    actual_positive = bool((torch.linalg.eigvalsh(result.matrix) >= 0).all())
    assert result.is_positive == actual_positive
    if actual_positive:
        root = result.sqrt().matrix
        torch.testing.assert_close(root @ root, result.matrix)


@pytest.mark.parametrize('new_first_eigenvalue', [2., -2.])
@pytest.mark.parametrize('precision', ['single', 'double'])
@pytest.mark.parametrize('complex_valued', [False, True])
def test_pending_square_rechecks_its_domain_or_computes_the_square(new_first_eigenvalue, precision, complex_valued):
    alg = algebra(precision, complex_valued)
    A = alg.operator(torch.eye(2, dtype=alg.dtype)[None, None])
    pending = A.power(2)
    with torch.no_grad():
        A.matrix.copy_(torch.diag(torch.tensor([new_first_eigenvalue, 3.], dtype=alg.dtype))[None, None])
    try:
        result = pending.matrix
    except (RuntimeError, ValueError):
        # Rejection is valid if a positive-only operation has lost its precondition.
        if new_first_eigenvalue < 0:
            return
        raise
    torch.testing.assert_close(result, A.matrix @ A.matrix)


@pytest.mark.parametrize('operation', ['sqrt', 'log', 'negative_power'])
@pytest.mark.parametrize('valid', [False, True])
@pytest.mark.parametrize('complex_valued', [False, True])
def test_other_pending_positive_functions_recheck_the_current_domain(operation, valid, complex_valued):
    alg = algebra(complex_valued=complex_valued)
    A = alg.operator(torch.eye(2, dtype=alg.dtype)[None, None])
    pending = A.power(-.5) if operation == 'negative_power' else getattr(A, operation)()
    values = torch.tensor([2. if valid else -2., 3.], dtype=alg.hilbert.real_dtype)
    with torch.no_grad():
        A.matrix.copy_(torch.diag(values).to(alg.dtype)[None, None])
    if not valid:
        with pytest.raises(RuntimeError, match='positive'):
            pending.matrix
    else:
        expected = {'sqrt': lambda: values.sqrt(), 'log': lambda: values.log(),
                    'negative_power': lambda: values.rsqrt()}[operation]()
        torch.testing.assert_close(pending.matrix, torch.diag(expected).to(alg.dtype)[None, None])


def _module_case(kind, alg):
    if kind == 'channel':
        return learn.ChannelLayer(alg, 2)
    if kind == 'intersector':
        out = TypeIAlgebra([2, 3], [2, 3], precision=alg.precision, device='cpu')
        return learn.ChannelLayer(alg, 2, algebra_out=out,
                                  transitions=[[True, True, False], [False, False, True]])
    if kind == 'unitary':
        return learn.UnitaryLayer(alg)
    if kind == 'lindblad':
        return learn.LindbladLayer(alg, time=.1, n_jumps=1, n_inter=1,
                                   transitions=[[True, False, False], [True, True, False], [False, True, True]])
    if kind == 'empty_lindblad':
        return learn.LindbladLayer(alg, hamiltonian=False, n_jumps=0)
    if kind == 'postselect':
        return learn.PostSelect(alg, [0, 1, 1])
    if kind == 'feedforward':
        instrument = measurement.Instrument.lueders(measurement.POVM.sectors(alg))
        return learn.FeedForward(instrument, [learn.ChannelLayer(alg, 1) for _ in range(len(instrument))])
    raise AssertionError(kind)


@pytest.mark.parametrize('kind', ['channel', 'intersector', 'unitary', 'lindblad', 'empty_lindblad',
                                  'postselect', 'feedforward'])
@pytest.mark.parametrize('precision', ['single', 'double'])
@pytest.mark.parametrize('device', ['cpu', pytest.param('cuda', marks=pytest.mark.skipif(
    not torch.cuda.is_available(), reason='CUDA is unavailable'))])
def test_padded_layers_transfer_with_their_metadata_and_gradients(kind, precision, device):
    torch.manual_seed(25)
    source_precision = 'double' if precision == 'single' else 'single'
    alg = TypeIAlgebra([1, 3, 2], [1, 3, 2], precision=source_precision, device='cpu')
    original_dtype = alg.dtype
    alg._const('eye')                    # A populated source cache must remain on the original algebra.
    model = _module_case(kind, alg)
    rho = states.random_density_matrix(alg, batch_size=2)
    with torch.no_grad():
        expected = model(rho).matrix.detach()
    dtype = torch.complex64 if precision == 'single' else torch.complex128
    model.to(device=device, dtype=dtype)
    target = alg.like(precision=precision, device=device)
    result = model(rho.cast(target))
    actual = result.matrix
    assert alg.dtype == original_dtype and str(alg.hilbert.device) == 'cpu'
    assert model.algebra.dtype == dtype and actual.dtype == dtype
    assert actual.device.type == device
    torch.testing.assert_close(actual, expected.to(device=device, dtype=dtype), atol=5e-5, rtol=5e-5)
    trace = actual.diagonal(dim1=-2, dim2=-1).sum((-2, -1)).real
    torch.testing.assert_close(trace, torch.ones_like(trace), atol=5e-5, rtol=5e-5)
    if list(model.parameters()):
        actual.abs().square().sum().backward()
        assert all(p.grad is not None and torch.isfinite(p.grad).all() for p in model.parameters())


@pytest.mark.parametrize('kind', ['channel', 'unitary', 'lindblad'])
@pytest.mark.parametrize('assign', [False, True])
def test_assigned_checkpoint_dtype_updates_algebra_metadata(kind, assign):
    factory = {'channel': learn.ChannelLayer, 'unitary': learn.UnitaryLayer, 'lindblad': learn.LindbladLayer}[kind]
    saved, restored = factory(algebra('single')), factory(algebra('double'))
    supports_assign = 'assign' in inspect.signature(restored.load_state_dict).parameters
    if assign and not supports_assign:
        pytest.skip('assign=True requires PyTorch 2.1 or newer')
    restored.load_state_dict(saved.state_dict(), **({'assign': assign} if supports_assign else {}))
    expected_dtype = torch.complex64 if assign else torch.complex128
    assert restored.algebra.dtype == expected_dtype
    result = restored(basis_state(restored.algebra))
    assert result.matrix.dtype == expected_dtype
    result.matrix.abs().square().sum().backward()
    assert all(p.grad is not None and torch.isfinite(p.grad).all() for p in restored.parameters())


@pytest.mark.parametrize('assign', [False, True])
def test_feedforward_serialized_checkpoint_restores_labels_measurement_and_dtype(assign):
    alg = algebra('single')
    povm = measurement.POVM.basis(alg)
    povm.outcomes = ['zero', 'one']
    saved = learn.FeedForward(measurement.Instrument.lueders(povm), [nn.Identity(), nn.Identity()])
    target = algebra('double')
    H = target.operator(torch.tensor([[[[1., 1.], [1., -1.]]]], dtype=target.dtype) / 2 ** .5)
    restored = learn.FeedForward(measurement.Instrument.lueders(measurement.POVM.from_unitary(H)),
                                 [nn.Identity(), nn.Identity()])
    stream = io.BytesIO()
    torch.save(saved.state_dict(), stream)
    stream.seek(0)
    kwargs = {'weights_only': True} if 'weights_only' in inspect.signature(torch.load).parameters else {}
    checkpoint = torch.load(stream, **kwargs)
    supports_assign = 'assign' in inspect.signature(restored.load_state_dict).parameters
    if assign and not supports_assign:
        pytest.skip('assign=True requires PyTorch 2.1 or newer')
    restored.load_state_dict(checkpoint, **({'assign': assign} if supports_assign else {}))
    assert restored.instrument.outcomes == ['zero', 'one']
    assert restored.algebra.dtype == (torch.complex64 if assign else torch.complex128)
    rho = basis_state(restored.algebra)
    torch.testing.assert_close(restored(rho).matrix, rho.matrix)


def test_legacy_feedforward_checkpoint_requires_explicit_missing_state_handling():
    alg = algebra()
    model = learn.FeedForward(measurement.Instrument.lueders(measurement.POVM.basis(alg)),
                              [nn.Identity(), nn.Identity()])
    with pytest.raises(RuntimeError, match='instrument_kraus'):
        model.load_state_dict({})
    missing = model.load_state_dict({}, strict=False).missing_keys
    assert 'instrument_kraus' in missing and '_extra_state' in missing


@pytest.mark.parametrize('kind', ['channel', 'lindblad', 'postselect'])
def test_real_layers_transfer_precision_without_changing_the_number_field(kind):
    alg = algebra('single', complex_valued=False)
    if kind == 'channel':
        model = learn.ChannelLayer(alg, 2)
    elif kind == 'lindblad':
        model = learn.LindbladLayer(alg, hamiltonian=False)
    else:
        model = learn.PostSelect(alg, [1])
    model.double()
    assert model.algebra.dtype == torch.float64 and not model.algebra.complex_valued
    assert alg.dtype == torch.float32
    result = model(basis_state(model.algebra))
    assert result.matrix.dtype == torch.float64
    torch.testing.assert_close(result.trace, torch.ones_like(result.trace))


@pytest.mark.parametrize('complex_valued', [False, True])
@pytest.mark.parametrize('cutoff', [0., 10.])
def test_general_pseudoinverse_drops_flags_invalidated_by_truncation(complex_valued, cutoff):
    alg = algebra(complex_valued=complex_valued)
    A = alg.operator(torch.tensor([[[[1., .2], [0., 1.]]]], dtype=alg.dtype))
    assert A.is_invertible and not A.is_normal
    result = A.inverse(tol=cutoff)
    expected = torch.tensor([[1., -.2], [0., 1.]], dtype=alg.dtype) if cutoff == 0 else torch.zeros(2, 2, dtype=alg.dtype)
    torch.testing.assert_close(result.matrix, expected[None, None])
    assert result.is_invertible == (cutoff == 0)
    assert result.is_normal == (cutoff != 0)


@pytest.mark.parametrize('coefficient', [1., 2.])
def test_nonzero_tensor_multiples_do_not_make_a_nonnormal_operator_normal(coefficient):
    alg = algebra()
    A = alg.operator(torch.tensor([[[[0., 1.], [0., 0.]]]], dtype=alg.dtype))
    assert not A.is_self_adjoint and not A.is_normal
    result = A * torch.tensor(coefficient)
    assert not result.is_self_adjoint and not result.is_normal


def test_pending_positive_function_can_be_retried_after_restoring_its_domain():
    alg = algebra()
    A = alg.operator(torch.eye(2, dtype=alg.dtype)[None, None])
    pending = A.sqrt()
    with torch.no_grad():
        A.matrix[0, 0, 0, 0] = -1.
    with pytest.raises(RuntimeError, match='positive'):
        pending.matrix
    with torch.no_grad():
        A.matrix[0, 0] = torch.diag(torch.tensor([4., 9.], dtype=alg.dtype))
    torch.testing.assert_close(pending.matrix, torch.diag(torch.tensor([2., 3.], dtype=alg.dtype))[None, None])
