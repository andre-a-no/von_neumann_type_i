"""GPU code paths (skipped without CUDA). The timing checks are in scripts/gpu_check.py."""
import pytest
import torch

from torch_vn_algebra import TypeIAlgebra, states, channels
from torch_vn_algebra.algebra import _orthonormal_columns, _gram_schmidt_preferred

pytestmark = pytest.mark.skipif(not torch.cuda.is_available(), reason="needs a CUDA device")


@pytest.mark.parametrize("dtype", [torch.float64, torch.complex128])
def test_gram_schmidt_on_gpu_equals_qr_with_phase_fix(dtype):
    torch.manual_seed(0)
    Z = torch.randn(256, 10, 10, dtype=dtype)
    Q, R = torch.linalg.qr(Z)
    d = torch.diagonal(R, dim1=-2, dim2=-1)
    d = d.sgn() if dtype.is_complex else d.sign()
    assert _gram_schmidt_preferred(Z.cuda())
    assert torch.allclose(_orthonormal_columns(Z.cuda()).cpu(), Q * d.unsqueeze(-2), atol=1e-10)


def test_random_objects_on_gpu():
    alg = TypeIAlgebra([3, 5, 5], [3, 5, 5], precision='double', device='cuda')
    U = alg.random_unitary_operator(64)
    assert torch.allclose((U @ U.adjoint()).matrix, alg.identity(64).matrix, atol=1e-12)
    A = alg.operator_from_eigenvalues(lambda k: 1 + torch.rand(64, k, device='cuda'), batch_size=64,
                                      force_positive=True)
    assert all((w >= 1 - 1e-12).all() for w in A.eigenvalues())
    rho = states.random_density_matrix(alg, 64)
    assert torch.allclose(rho.Tr_blunt().real, torch.ones(64, dtype=torch.float64, device='cuda'))
    Phi = channels.random_channel(alg, 2, 8)
    assert Phi.is_trace_preserving()
