"""Mutable Kraus tensors must not leave stale state-preservation classifications."""
import warnings

import pytest
import torch

from torch_vn_algebra import TypeIAlgebra, channels, states


def _case(kind, known=None, complex_valued=True):
    alg = TypeIAlgebra([2], [2], precision='double', complex_valued=complex_valued, device='cpu')
    kraus = torch.eye(2, dtype=alg.dtype)[None, None, None]
    if kind == 'sector':
        phi = channels.Channel(alg, kraus, trace_preserving=known)
    else:
        phi = channels.InterSectorChannel(alg, alg, kraus.unsqueeze(2), trace_preserving=known)
    rho = states.DensityMatrix(alg, matrix=torch.diag(torch.tensor([.2, .8], dtype=alg.dtype))[None, None])
    return alg, phi, rho


@pytest.mark.parametrize('kind', ['sector', 'intersector'])
@pytest.mark.parametrize('known', [None, True])
@pytest.mark.parametrize('complex_valued', [False, True])
@pytest.mark.parametrize('replace', [False, True])
def test_changed_kraus_rechecks_trace_preservation(kind, known, complex_valued, replace):
    _, phi, rho = _case(kind, known, complex_valued)
    torch.testing.assert_close(phi(rho).matrix, rho.matrix)
    with torch.no_grad():
        if replace:
            phi.kraus = 2 * phi.kraus
        else:
            phi.kraus.mul_(2)
    assert not phi.is_trace_preserving()
    with pytest.warns(UserWarning, match='not trace preserving'):
        result = phi(rho)
    assert not isinstance(result, states.DensityMatrix)
    # Multiplying every Kraus operator by two multiplies the CP map by four.
    torch.testing.assert_close(result.matrix, 4 * rho.matrix)


@pytest.mark.parametrize('kind', ['sector', 'intersector'])
def test_nonpreserving_cache_can_become_a_state_map(kind):
    _, phi, rho = _case(kind)
    phi.kraus.mul_(2)
    with pytest.warns(UserWarning, match='not trace preserving'):
        result = phi(rho)
    torch.testing.assert_close(result.matrix, 4 * rho.matrix)
    phi.kraus.div_(2)
    with warnings.catch_warnings(record=True) as caught:
        result = phi(rho)
        torch.testing.assert_close(result.matrix, rho.matrix)
    assert isinstance(result, states.DensityMatrix)
    assert not caught


@pytest.mark.parametrize('kind', ['sector', 'intersector'])
@pytest.mark.parametrize('order', ['left', 'right'])
def test_composition_does_not_propagate_a_stale_preserving_hint(kind, order):
    alg, phi, rho = _case(kind, known=True)
    identity = channels.identity_channel(alg)
    if kind == 'intersector':
        identity = identity.to_inter_sector()
    phi.kraus.mul_(2)
    composed = phi.compose(identity) if order == 'left' else identity.compose(phi)
    with pytest.warns(UserWarning, match='not trace preserving'):
        result = composed(rho)
    assert not isinstance(result, states.DensityMatrix)
    torch.testing.assert_close(result.matrix, 4 * rho.matrix)


def test_sector_conversion_does_not_propagate_a_stale_preserving_hint():
    _, phi, rho = _case('sector', known=True)
    phi.kraus.mul_(2)
    converted = phi.to_inter_sector()
    with pytest.warns(UserWarning, match='not trace preserving'):
        result = converted(rho)
    assert not isinstance(result, states.DensityMatrix)
    torch.testing.assert_close(result.matrix, 4 * rho.matrix)


@pytest.mark.parametrize('order', ['left', 'right'])
def test_mixture_does_not_propagate_a_stale_preserving_hint(order):
    alg, phi, rho = _case('sector', known=True)
    identity = channels.identity_channel(alg)
    phi.kraus.mul_(2)
    mixed = phi.mix(identity, .5) if order == 'left' else identity.mix(phi, .5)
    with pytest.warns(UserWarning, match='not trace preserving'):
        result = mixed(rho)
    assert not isinstance(result, states.DensityMatrix)
    torch.testing.assert_close(result.matrix, 2.5 * rho.matrix)


@pytest.mark.parametrize('kind', ['sector', 'intersector'])
def test_pending_state_rejects_lost_preservation_and_can_be_retried(kind):
    _, phi, rho = _case(kind, known=True)
    pending = phi(rho)
    phi.kraus.mul_(2)
    with pytest.raises(RuntimeError, match='changed.*no longer trace preserving'):
        pending.matrix
    phi.kraus.div_(2)
    torch.testing.assert_close(pending.matrix, rho.matrix)


@pytest.mark.parametrize('kind', ['sector', 'intersector'])
def test_pending_state_follows_a_mutation_that_remains_trace_preserving(kind):
    _, phi, rho = _case(kind, known=True)
    pending = phi(rho)
    swap = torch.tensor([[0., 1.], [1., 0.]], dtype=phi.kraus.dtype)
    phi.kraus.copy_(swap.expand_as(phi.kraus))
    torch.testing.assert_close(pending.matrix, (swap @ rho.matrix @ swap.T))


@pytest.mark.parametrize('kind', ['sector', 'intersector'])
def test_replacing_channel_tensor_does_not_rebind_an_existing_lazy_recipe(kind):
    _, phi, rho = _case(kind, known=True)
    pending = phi(rho)
    phi.kraus = 2 * phi.kraus
    # The recipe captured the original tensor; replacing the attribute does not mutate it.
    torch.testing.assert_close(pending.matrix, rho.matrix)
    with pytest.warns(UserWarning, match='not trace preserving'):
        future = phi(rho)
    torch.testing.assert_close(future.matrix, 4 * rho.matrix)


@pytest.mark.parametrize('kind', ['sector', 'intersector'])
def test_unmodified_versioned_channel_keeps_its_one_time_check(kind, monkeypatch):
    _, phi, rho = _case(kind)
    calls = []
    original = phi.is_trace_preserving

    def counted():
        calls.append(None)
        return original()

    monkeypatch.setattr(phi, 'is_trace_preserving', counted)
    for _ in range(3):
        torch.testing.assert_close(phi(rho).matrix, rho.matrix)
    assert len(calls) == 1
    phi.kraus.mul_(2)
    with pytest.warns(UserWarning, match='not trace preserving'):
        phi(rho)
    assert len(calls) == 2


@pytest.mark.parametrize('kind', ['sector', 'intersector'])
def test_inference_tensors_without_versions_are_rechecked(kind):
    with torch.inference_mode():
        _, phi, rho = _case(kind, known=True)
        torch.testing.assert_close(phi(rho).matrix, rho.matrix)
        phi.kraus.mul_(2)
        with pytest.warns(UserWarning, match='not trace preserving'):
            result = phi(rho)
        assert not isinstance(result, states.DensityMatrix)
        torch.testing.assert_close(result.matrix, 4 * rho.matrix)
