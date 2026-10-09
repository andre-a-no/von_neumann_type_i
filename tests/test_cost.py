"""Tests for the time / memory safeguards."""
import warnings

import pytest
import torch

import torch_vn_algebra as tv
from torch_vn_algebra import TypeIAlgebra, cost, dynamics


@pytest.fixture(autouse=True)
def restore_limits():
    old = cost.get_limits()
    yield
    cost.set_limits(**old)


def test_time_warning_for_slow_kernel():
    cost.set_limits(warn_seconds=0.0, probe_flops=0)
    alg = TypeIAlgebra([32], [32], device='cpu')
    H = alg.operator_from_eigenvalues(lambda d: torch.randn(200, d), batch_size=200, force_self_adjoint=True)
    with pytest.warns(cost.CostWarning, match="expected"):
        H.eigh()


def test_probed_result_equals_direct_call():
    cost.set_limits(warn_seconds=1e9, probe_flops=0)
    x = torch.randn(100, 3, 6, 6, dtype=torch.float64)
    x = x + x.transpose(-2, -1)
    w1 = cost.batched_call(torch.linalg.eigvalsh, x, "test")
    assert torch.allclose(w1, torch.linalg.eigvalsh(x))
    U, S, Vh = cost.batched_call(torch.linalg.svd, x, "test", kind='svd')
    assert torch.allclose(U @ torch.diag_embed(S) @ Vh, x)


def test_memory_error_and_warning():
    cost.set_limits(max_memory_fraction=1e-9)
    alg = TypeIAlgebra([64], [64], device='cpu')
    with pytest.raises(cost.InsufficientMemoryError):
        alg.operator_from_eigenvalues(lambda d: torch.rand(10_000, d), batch_size=10_000)
    cost.set_limits(max_memory_fraction=1.0, warn_memory_fraction=1e-9)
    with pytest.warns(cost.CostWarning, match="will use"):
        alg.operator_from_eigenvalues(lambda d: torch.rand(10_000, d), batch_size=10_000)


def test_disabled_context():
    cost.set_limits(warn_seconds=0.0, probe_flops=0, max_memory_fraction=1e-9)
    alg = TypeIAlgebra([32], [32], device='cpu')
    with cost.disabled(), warnings.catch_warnings():
        warnings.simplefilter('error')
        H = alg.operator_from_eigenvalues(lambda d: torch.randn(200, d), batch_size=200, force_self_adjoint=True)
        H.eigh()


def test_rk4_eta_warning_and_progress(capsys):
    cost.set_limits(warn_seconds=0.0)
    y0 = torch.ones(3)
    with pytest.warns(cost.CostWarning, match="rk4"):
        out = dynamics.rk4(lambda t, y: -y, y0, [0.0, 1.0], substeps=50, progress=True)
    assert torch.allclose(out[-1], torch.exp(torch.tensor(-1.0)) * y0, atol=1e-6)


def test_unknown_option():
    with pytest.raises(KeyError):
        tv.cost.set_limits(warn_second=1)
