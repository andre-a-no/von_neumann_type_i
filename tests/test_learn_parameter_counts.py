"""Real degrees of freedom from block dimensions, independent of padded tensor storage."""
import warnings

import pytest
import torch

from torch_vn_algebra import TypeIAlgebra, learn


def _algebra(sizes, complex_):
    return TypeIAlgebra([max(1, k) for k in sizes], sizes,
                        complex_valued=complex_, precision='double')


@pytest.mark.parametrize('complex_', [False, True])
@pytest.mark.parametrize('sizes', [[2], [1, 3, 0]])
def test_channel_count_uses_the_field_of_each_active_entry(complex_, sizes):
    layer = learn.ChannelLayer(_algebra(sizes, complex_), kraus_rank=2)
    # Each of the two Kraus operators has k_c^2 active entries in sector c.
    expected = (2 if complex_ else 1) * 2 * sum(k * k for k in sizes)
    assert learn.count_parameters(layer) == expected
    assert expected <= learn.count_parameters(layer, effective=False)


@pytest.mark.parametrize('complex_', [False, True])
def test_intersector_channel_count_uses_allowed_rectangular_blocks(complex_):
    source = _algebra([2, 1], complex_)
    target = _algebra([1, 3], complex_)
    transitions = [[True, False], [True, True]]
    layer = learn.ChannelLayer(source, kraus_rank=2, algebra_out=target,
                               transitions=transitions)
    # Allowed rectangular blocks have 1*2, 3*2 and 3*1 entries.
    assert learn.count_parameters(layer) == (2 if complex_ else 1) * 2 * (2 + 6 + 3)


@pytest.mark.parametrize('complex_', [False, True])
@pytest.mark.parametrize('hamiltonian', [False, True])
def test_lindblad_count_respects_the_field_and_self_adjoint_constraint(complex_, hamiltonian):
    sizes = [1, 3, 0]
    layer = learn.LindbladLayer(_algebra(sizes, complex_), n_jumps=2,
                               n_inter=1, transitions=[[1, 0, 0], [1, 1, 0], [0, 0, 0]],
                               hamiltonian=hamiltonian)
    field = 2 if complex_ else 1
    expected = field * (2 * (1 + 9) + (1 + 3 + 9))
    if hamiltonian:
        # Hermitian k-by-k matrices have k^2 real coordinates; real symmetric ones k(k+1)/2.
        expected += sum(k * k if complex_ else k * (k + 1) // 2 for k in sizes)
    assert learn.count_parameters(layer) == expected
    assert expected <= learn.count_parameters(layer, effective=False)


@pytest.mark.parametrize('dtype', [torch.float32, torch.float64, torch.complex64, torch.complex128])
def test_channel_count_tracks_module_dtype_transfers(dtype):
    layer = learn.ChannelLayer(_algebra([2, 1], True), kraus_rank=2)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        layer.to(dtype=dtype)
    expected = (2 if dtype.is_complex else 1) * 2 * (4 + 1)
    assert learn.count_parameters(layer) == expected
    assert expected <= learn.count_parameters(layer, effective=False)
