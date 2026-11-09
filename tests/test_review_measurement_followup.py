"""Regressions for differentiable measurements and mutable instrument Kraus data."""
import warnings

import pytest
import torch

from torch_vn_algebra import Instrument, POVM, PVM, TypeIAlgebra, states
from torch_vn_algebra.measurement import fisher_information


@pytest.mark.parametrize('tol', [None, 0.0])
def test_pushforward_preserves_parameter_derivative_at_equal_values(tol):
    alg = TypeIAlgebra([2], [2], precision='double', device='cpu')
    measure = PVM.from_povm(POVM.basis(alg), torch.tensor([-1., 1.], dtype=torch.float64))

    def observable(t):
        return measure.pushforward(lambda values: values.square() + t * values, tol=tol).observable().matrix.real

    t = torch.tensor(0., dtype=torch.float64, requires_grad=True)
    actual = observable(t)
    derivative = torch.autograd.functional.jacobian(observable, t)
    # A = diag(-1, 1), hence f_t(A) = I + t A, even where the two values coincide.
    expected_derivative = torch.diag(torch.tensor([-1., 1.], dtype=torch.float64)).reshape(1, 1, 2, 2)
    torch.testing.assert_close(actual, torch.eye(2, dtype=torch.float64).reshape(1, 1, 2, 2))
    torch.testing.assert_close(derivative, expected_derivative)
    # This step also separates the values under the default numerical grouping tolerance.
    step = 1e-4
    finite_difference = (observable(t.detach() + step) - observable(t.detach() - step)) / (2 * step)
    torch.testing.assert_close(derivative, finite_difference, atol=1e-10, rtol=0)


def test_pushforward_preserves_input_value_derivatives_at_equal_values():
    alg = TypeIAlgebra([2], [2], precision='double', device='cpu')
    basis = POVM.basis(alg)

    def diagonal(values):
        measure = PVM.from_povm(basis, values)
        matrix = measure.pushforward(torch.square, tol=0).observable().matrix[0, 0].real
        return torch.diagonal(matrix)

    values = torch.tensor([-1., 1.], dtype=torch.float64, requires_grad=True)
    derivative = torch.autograd.functional.jacobian(diagonal, values)
    # Squaring a diagonal observable squares each entry independently.
    expected = torch.diag(2 * values.detach())
    torch.testing.assert_close(diagonal(values), torch.ones(2, dtype=torch.float64))
    torch.testing.assert_close(derivative, expected)
    step = 1e-6
    directions = torch.eye(2, dtype=torch.float64)
    finite_difference = torch.stack([
        (diagonal(values.detach() + step * direction) - diagonal(values.detach() - step * direction)) / (2 * step)
        for direction in directions
    ], dim=-1)
    torch.testing.assert_close(derivative, finite_difference, atol=1e-9, rtol=0)


def test_pushforward_still_merges_equal_values_without_gradients():
    alg = TypeIAlgebra([3], [3], precision='double', device='cpu')
    measure = PVM.from_povm(POVM.basis(alg), torch.tensor([-1., 1., 2.], dtype=torch.float64))
    squared = measure.pushforward(torch.square)
    assert len(squared) == 2
    torch.testing.assert_close(squared.values, torch.tensor([[1., 4.]], dtype=torch.float64))
    expected_effects = torch.stack([
        torch.diag(torch.tensor([1., 1., 0.], dtype=alg.dtype)),
        torch.diag(torch.tensor([0., 0., 1.], dtype=alg.dtype)),
    ]).reshape(1, 2, 1, 3, 3)
    torch.testing.assert_close(squared.effects, expected_effects)
    expected_observable = torch.diag(torch.tensor([1., 1., 4.], dtype=alg.dtype)).reshape(1, 1, 3, 3)
    torch.testing.assert_close(squared.observable().matrix, expected_observable)


@pytest.mark.parametrize('per_sample', [False, True])
def test_fisher_information_accepts_conjugate_complex_parameters(per_sample):
    values = [0.2 + 0.3j, -0.4 + 0.1j] if per_sample else 0.2 + 0.3j
    parameter = torch.tensor(values, dtype=torch.complex128).conj()
    assert parameter.is_conj()
    offsets = torch.tensor([-0.5, 0.7], dtype=torch.float64)

    def probabilities(z):
        p = torch.sigmoid(z.real + 2 * z.imag + offsets)
        return torch.stack([p, 1 - p], dim=-1)

    actual = fisher_information(probabilities, [parameter], per_sample=per_sample)
    resolved = fisher_information(probabilities, [parameter.resolve_conj()], per_sample=per_sample)
    # The Bernoulli logit has gradient (1, 2) in the real and imaginary coordinates.
    p = torch.sigmoid(parameter.real + 2 * parameter.imag + offsets)
    direction = torch.tensor([1., 2.], dtype=torch.float64)
    expected = (p * (1 - p))[:, None, None] * torch.outer(direction, direction)
    torch.testing.assert_close(actual, expected)
    torch.testing.assert_close(actual, resolved)


def _instrument_case(known):
    alg = TypeIAlgebra([2], [2], precision='double', device='cpu')
    povm = POVM.basis(alg)
    if known == 'constructed':
        instrument = Instrument.lueders(povm)
    else:
        instrument = Instrument(povm.effects.unsqueeze(2), alg, check=known == 'checked')
    rho = states.DensityMatrix(alg, matrix=torch.tensor([[[[.25, .1], [.1, .75]]]], dtype=alg.dtype))
    return instrument, rho


def _expected_instrument_output(instrument, export):
    if export == 'average_channel':
        return torch.diag(torch.tensor([.25, .75], dtype=instrument.kraus.dtype)).reshape(1, 1, 2, 2)
    # Outcome zero contains |0><0| with weight .25; outcome one contains |1><1| with weight .75.
    return torch.tensor([[[[.25, 0.], [0., 0.]], [[0., 0.], [0., .75]]]], dtype=instrument.kraus.dtype)


@pytest.mark.parametrize('export', ['average_channel', 'as_channel'])
@pytest.mark.parametrize('replace', [False, True])
@pytest.mark.parametrize('known', ['constructed', 'checked', 'unknown'])
def test_instrument_exports_recheck_mutated_kraus(export, replace, known):
    instrument, rho = _instrument_case(known)
    expected = _expected_instrument_output(instrument, export)
    initial = getattr(instrument, export)()(rho)
    assert isinstance(initial, states.DensityMatrix)
    torch.testing.assert_close(initial.matrix, expected)

    if replace:
        instrument.kraus = 2 * instrument.kraus
    else:
        instrument.kraus.mul_(2)
    changed = getattr(instrument, export)()
    assert not changed.is_trace_preserving()
    with pytest.warns(UserWarning, match='not trace preserving'):
        result = changed(rho)
    assert not isinstance(result, states.DensityMatrix)
    # Scaling every Kraus operator by two multiplies every branch by four.
    torch.testing.assert_close(result.matrix, 4 * expected)

    if replace:
        instrument.kraus = instrument.kraus / 2
    else:
        instrument.kraus.div_(2)
    with warnings.catch_warnings(record=True) as caught:
        restored = getattr(instrument, export)()(rho)
        torch.testing.assert_close(restored.matrix, expected)
    assert isinstance(restored, states.DensityMatrix)
    assert not caught


@pytest.mark.parametrize('export', ['average_channel', 'as_channel'])
@pytest.mark.parametrize('replace', [False, True])
def test_instrument_exports_drop_nonpreserving_hint_after_restoration(export, replace):
    instrument, rho = _instrument_case('unknown')
    expected = _expected_instrument_output(instrument, export)
    instrument = Instrument(2 * instrument.kraus, instrument.algebra, check=False)
    # A stored classification may also be False; it describes only the current Kraus data.
    instrument._trace_preserving = False
    with pytest.warns(UserWarning, match='not trace preserving'):
        invalid = getattr(instrument, export)()(rho)
    torch.testing.assert_close(invalid.matrix, 4 * expected)

    if replace:
        instrument.kraus = instrument.kraus / 2
    else:
        instrument.kraus.div_(2)
    with warnings.catch_warnings(record=True) as caught:
        restored = getattr(instrument, export)()(rho)
        torch.testing.assert_close(restored.matrix, expected)
    assert isinstance(restored, states.DensityMatrix)
    assert not caught


@pytest.mark.parametrize('export', ['average_channel', 'as_channel'])
@pytest.mark.parametrize('replace', [False, True])
@pytest.mark.parametrize('known', ['constructed', 'unknown'])
def test_instrument_exports_keep_state_type_after_unitary_kraus_rotation(export, replace, known):
    instrument, rho = _instrument_case(known)
    expected = _expected_instrument_output(instrument, export)
    # Every outcome undergoes a Pauli X after measurement; K_j^* K_j is unchanged.
    swap = torch.tensor([[0., 1.], [1., 0.]], dtype=instrument.kraus.dtype)
    rotated = swap @ instrument.kraus
    if replace:
        instrument.kraus = rotated
    else:
        instrument.kraus.copy_(rotated)
    channel = getattr(instrument, export)()
    assert channel.is_trace_preserving()
    with warnings.catch_warnings(record=True) as caught:
        result = channel(rho)
        torch.testing.assert_close(result.matrix, swap @ expected @ swap.T)
    assert isinstance(result, states.DensityMatrix)
    assert not caught
