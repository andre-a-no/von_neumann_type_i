"""Checkpointed Kraus buffers must not export obsolete trace-preservation hints."""
import pytest
import torch
from torch import nn

from torch_vn_algebra import Instrument, POVM, TypeIAlgebra, learn, states


@pytest.mark.parametrize('complex_', [False, True])
@pytest.mark.parametrize('export', ['average_channel', 'as_channel'])
def test_feed_forward_instrument_exports_check_the_current_buffer(complex_, export):
    alg = TypeIAlgebra([2], [2], precision='double', complex_valued=complex_)
    instrument = Instrument.lueders(POVM.basis(alg))
    model = learn.FeedForward(instrument, [nn.Identity(), nn.Identity()])
    rho = states.DensityMatrix(alg, matrix=torch.diag(torch.tensor([.2, .8], dtype=alg.dtype))[None, None])
    original = getattr(model.instrument, export)()(rho)
    assert isinstance(original, states.DensityMatrix)
    torch.testing.assert_close(torch.diagonal(original.matrix, dim1=-2, dim2=-1).sum(),
                               torch.tensor(1., dtype=alg.dtype))

    model.instrument_kraus.mul_(2)
    channel = getattr(model.instrument, export)()
    assert not channel.is_trace_preserving()
    with pytest.warns(UserWarning, match='not trace preserving'):
        result = channel(rho)
    assert not isinstance(result, states.DensityMatrix)
    # K -> 2K makes the entire instrument's unnormalised map four times larger.
    torch.testing.assert_close(torch.diagonal(result.matrix, dim1=-2, dim2=-1).sum(),
                               torch.tensor(4., dtype=alg.dtype))


def test_feed_forward_exports_still_recognise_a_trace_preserving_buffer_change():
    alg = TypeIAlgebra([2], [2], precision='double')
    model = learn.FeedForward(Instrument.lueders(POVM.basis(alg)), [nn.Identity(), nn.Identity()])
    model.instrument_kraus.mul_(-1)
    rho = states.maximally_mixed_state(alg)
    result = model.instrument.average_channel()(rho)
    assert isinstance(result, states.DensityMatrix)
    torch.testing.assert_close(result.matrix, rho.matrix)
