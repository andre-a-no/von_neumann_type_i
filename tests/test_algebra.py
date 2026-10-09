import math
"""
Tests for TypeIAlgebra and Operator classes, including SU(n) generation and eigenvalue extraction.

Run with: pytest tests/test_algebra.py -v
"""
from torch_vn_algebra import TypeIAlgebra

import torch
import pytest
import numpy as np


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def basic_algebra():
    n_factors = [2, 3]
    k_factors = [2, 3]
    return TypeIAlgebra(
        n_factors=n_factors,
        k_factors=k_factors,
        batch_size=2,
        complex_valued=True,
        device='cpu'
    )


@pytest.fixture
def real_algebra():
    n_factors = [2, 2]
    k_factors = [2, 2]
    return TypeIAlgebra(
        n_factors=n_factors,
        k_factors=k_factors,
        batch_size=1,
        complex_valued=False,
        device='cpu'
    )


@pytest.fixture
def identity_operator(basic_algebra):
    return basic_algebra.identity(batch_size=2)


@pytest.fixture
def zero_operator(basic_algebra):
    return basic_algebra.zero(batch_size=2)


# ============================================================================
# Algebra initialization
# ============================================================================

class TestAlgebraInit:
    def test_valid_init(self):
        n_factors = [2, 3]
        k_factors = [2, 3]
        alg = TypeIAlgebra(n_factors, k_factors)
        assert alg.C == 2
        assert alg.k_max == 3
        assert alg.total_subspace_dim == 5

    def test_invalid_k_factors(self):
        with pytest.raises(ValueError):
            TypeIAlgebra([2, 3], [3, 4])
        with pytest.raises(ValueError):
            TypeIAlgebra([1], [0])                  # no sector with k_c > 0

    def test_mismatched_lengths(self):
        with pytest.raises(AssertionError):
            TypeIAlgebra([2, 3], [2])


# ============================================================================
# Basic operators
# ============================================================================

class TestBasicOperators:
    def test_identity_shape(self, identity_operator):
        assert identity_operator.shape == (2, 2, 3, 3)

    def test_identity_trace(self, identity_operator):
        assert torch.allclose(identity_operator.trace.real, torch.tensor(5.0))

    def test_identity_is_self_adjoint(self, identity_operator):
        assert identity_operator.is_self_adjoint

    def test_identity_is_positive(self, identity_operator):
        assert identity_operator.is_positive

    def test_zero_shape(self, zero_operator):
        assert zero_operator.shape == (2, 2, 3, 3)

    def test_zero_trace(self, zero_operator):
        assert torch.allclose(zero_operator.trace.real, torch.tensor(0.0))

    def test_central_operator(self, basic_algebra):
        central = basic_algebra.central([2.0, 3.0], batch_size=2)
        assert central.shape == (2, 2, 3, 3)
        mat = central.matrix
        assert torch.allclose(mat[0, 0, 0, 0].real, torch.tensor(2.0))
        assert torch.allclose(mat[0, 0, 1, 1].real, torch.tensor(2.0))
        assert torch.allclose(mat[0, 1, 0, 0].real, torch.tensor(3.0))


# ============================================================================
# Random operator generation
# ============================================================================

class TestRandomOperatorGeneration:
    def test_positive_generation(self, basic_algebra):
        def sampler(dim):
            return torch.ones(dim)
        op = basic_algebra.operator_from_eigenvalues(sampler, batch_size=2, force_positive=True)
        assert op.is_positive
        assert torch.allclose(op.trace.real, torch.tensor(5.0))

    def test_self_adjoint_generation(self, basic_algebra):
        def sampler(dim):
            return torch.randn(dim)
        op = basic_algebra.operator_from_eigenvalues(sampler, batch_size=2, force_self_adjoint=True)
        assert op.is_self_adjoint

    def test_projection_generation(self, basic_algebra):
        def sampler(dim):
            if dim == 2:
                return torch.tensor([1.0, 1.0])
            else:
                return torch.tensor([1.0, 1.0, 0.0])
        op = basic_algebra.operator_from_eigenvalues(sampler, batch_size=2, force_projection=True)
        assert op.is_projection

    def test_negative_eigenvalue_error(self, basic_algebra):
        def sampler(dim):
            return torch.tensor([-0.5, 0.5]) if dim == 2 else torch.tensor([0.5, 0.5, 0.5])
        with pytest.raises(ValueError, match="force_positive=True but eigenvalues for channel 0 are negative"):
            basic_algebra.operator_from_eigenvalues(sampler, batch_size=2, force_positive=True)


# ============================================================================
# Operator properties
# ============================================================================

class TestOperatorProperties:
    def test_is_self_adjoint_from_matrix(self, basic_algebra):
        def sampler(dim):
            return torch.randn(dim)
        op = basic_algebra.operator_from_eigenvalues(sampler, force_self_adjoint=True)
        assert op.is_self_adjoint

    def test_is_positive_from_matrix(self, basic_algebra):
        def sampler(dim):
            return torch.abs(torch.randn(dim))
        op = basic_algebra.operator_from_eigenvalues(sampler, force_positive=True)
        assert op.is_positive

    def test_is_normal_self_adjoint(self, basic_algebra):
        def sampler(dim):
            return torch.randn(dim)
        op = basic_algebra.operator_from_eigenvalues(sampler, force_self_adjoint=True)
        assert op.is_normal


# ============================================================================
# Lazy materialization
# ============================================================================

class TestLazyMaterialization:
    def test_lazy_generation(self, basic_algebra):
        def sampler(dim):
            return torch.randn(dim)
        op = basic_algebra.operator_from_eigenvalues(sampler, batch_size=2)
        assert op._matrix is None
        _ = op.matrix
        assert op._matrix is not None


# ============================================================================
# Arithmetic operations
# ============================================================================

class TestArithmetic:
    def test_addition(self, basic_algebra):
        op1 = basic_algebra.identity(batch_size=2)
        op2 = basic_algebra.identity(batch_size=2)
        op_sum = op1 + op2
        assert torch.allclose(op_sum.trace.real, torch.tensor(10.0))

    def test_multiplication(self, basic_algebra):
        op1 = basic_algebra.identity(batch_size=2)
        op2 = basic_algebra.identity(batch_size=2)
        op_prod = op1 @ op2
        assert torch.allclose(op_prod.trace.real, torch.tensor(5.0))

    def test_scalar_multiplication(self, basic_algebra):
        op = basic_algebra.identity(batch_size=2)
        op_scaled = op * 2.0
        assert torch.allclose(op_scaled.trace.real, torch.tensor(10.0))

    def test_adjoint(self, basic_algebra):
        op = basic_algebra.identity(batch_size=2)
        op_adj = op.adjoint()
        assert torch.allclose(op_adj.matrix, op.matrix)


# ============================================================================
# Inverse and absolute value
# ============================================================================

class TestInverseAndAbs:
    def test_inverse_identity(self, basic_algebra):
        op = basic_algebra.identity(batch_size=2)
        inv_op = op.inverse()
        assert torch.allclose(inv_op.matrix, op.matrix)

    def test_abs_positive(self, basic_algebra):
        def sampler(dim):
            return torch.abs(torch.randn(dim))
        op = basic_algebra.operator_from_eigenvalues(sampler, force_positive=True)
        abs_op = op.abs()
        assert torch.allclose(abs_op.matrix, op.matrix, atol=1e-6, rtol=1e-5)

    def test_sqrt_positive(self, basic_algebra):
        def sampler(dim):
            return torch.abs(torch.randn(dim))
        op = basic_algebra.operator_from_eigenvalues(sampler, force_positive=True)
        sqrt_op = op.sqrt()
        assert torch.allclose((sqrt_op @ sqrt_op).matrix, op.matrix, atol=1e-6)


# ============================================================================
# Entropy
# ============================================================================

class TestEntropy:
    def test_entropy_identity(self, basic_algebra):
        op = basic_algebra.identity(batch_size=2)
        entropy = op.entropy()
        assert torch.allclose(entropy, torch.tensor(0.0), atol=1e-6)

    def test_entropy_pure_state(self, basic_algebra):
        def sampler(dim):
            if dim == 2:
                return torch.tensor([1.0, 0.0])
            else:
                return torch.tensor([1.0, 0.0, 0.0])
        op = basic_algebra.operator_from_eigenvalues(sampler, batch_size=2, force_positive=True)
        entropy = op.entropy()
        assert torch.allclose(entropy, torch.tensor(0.0), atol=1e-5)


# ============================================================================
# Michelson contrast
# ============================================================================

class TestMichelsonContrast:
    def test_contrast_identity(self, basic_algebra):
        op = basic_algebra.identity(batch_size=2)
        assert torch.all(op.michelson_contrast.abs() < 1e-6)

    def test_contrast_singular(self, basic_algebra):
        def sampler(dim):
            if dim == 2:
                return torch.tensor([1.0, 0.0])
            else:
                return torch.tensor([1.0, 1.0, 0.0])
        op = basic_algebra.operator_from_eigenvalues(sampler, batch_size=1, force_projection=True)
        assert torch.all((op.michelson_contrast - 1.0).abs() < 1e-6)


# ============================================================================
# Three traces
# ============================================================================

class TestThreeTraces:
    def test_tr_blunt_identity(self, basic_algebra):
        op = basic_algebra.identity(batch_size=2)
        assert torch.allclose(op.Tr_blunt().real, torch.tensor(5.0))

    def test_tr_norm_identity(self, basic_algebra):
        op = basic_algebra.identity(batch_size=2)
        # channel0: trace=2 /2 =1; channel1: trace=3/3=1; sum=2
        expected = torch.tensor(2.0)
        assert torch.allclose(op.Tr_norm().real, expected)

    def test_tau_vN_identity(self, basic_algebra):
        op = basic_algebra.identity(batch_size=2)
        assert torch.allclose(op.tau_vN().real, torch.tensor(1.0))
        assert torch.is_complex(op.tau_vN())   # traces live in the field of the algebra


# ============================================================================
# Error handling
# ============================================================================

class TestErrorHandling:
    def test_sqrt_non_positive(self, basic_algebra):
        def sampler(dim):
            eig = torch.randn(dim)
            eig[0] = -0.5
            return eig
        op = basic_algebra.operator_from_eigenvalues(sampler, force_self_adjoint=True)
        with pytest.raises(RuntimeError, match="sqrt requires positive operator"):
            op.sqrt()

    def test_entropy_non_positive(self, basic_algebra):
        def sampler(dim):
            eig = torch.randn(dim)
            eig[0] = -0.5
            return eig
        op = basic_algebra.operator_from_eigenvalues(sampler, force_self_adjoint=True)
        with pytest.raises(RuntimeError, match="trace_a_log_a requires positive operator"):
            op.entropy()


# ============================================================================
# Unitary measures (including SU(n))
# ============================================================================

class TestUnitaryMeasures:
    def test_haar_unitary(self, basic_algebra):
        U = basic_algebra.random_unitary(4, measure='haar')
        expected = torch.eye(4, dtype=U.dtype, device=U.device)
        assert torch.allclose(U @ U.conj().T, expected, atol=1e-6)

    def test_haar_su_unitary(self, basic_algebra):
        U = basic_algebra.random_unitary(3, measure='haar_su')
        det = torch.linalg.det(U)
        assert torch.allclose(det, torch.tensor(1.0+0.0j, dtype=U.dtype, device=U.device), atol=1e-6)

    def test_haar_su_unitary_real(self, real_algebra):
        U = real_algebra.random_unitary(3, measure='haar_su')
        expected = torch.eye(3, dtype=U.dtype, device=U.device)
        assert torch.allclose(U @ U.T, expected, atol=1e-6)
        det = torch.linalg.det(U)
        assert torch.allclose(torch.abs(det), torch.tensor(1.0, dtype=det.dtype, device=det.device), atol=1e-6)

    def test_coe_symmetric_unitary(self, basic_algebra):
        U = basic_algebra.random_unitary(4, measure='coe')
        expected = torch.eye(4, dtype=U.dtype, device=U.device)
        assert torch.allclose(U @ U.conj().T, expected, atol=1e-5)
        assert torch.allclose(U, U.T, atol=1e-5)

    def test_cse_self_dual_unitary(self, basic_algebra):
        U = basic_algebra.random_unitary(4, measure='cse')
        expected = torch.eye(4, dtype=U.dtype, device=U.device)
        assert torch.allclose(U @ U.conj().T, expected, atol=1e-5)
        I2, O2 = torch.eye(2, dtype=U.dtype), torch.zeros(2, 2, dtype=U.dtype)
        J = torch.cat([torch.cat([O2, I2], 1), torch.cat([-I2, O2], 1)], 0)
        assert torch.allclose(J @ U.T @ J.T, U, atol=1e-5)

    def test_complex_only_measures_reject_real(self, real_algebra):
        for measure in ('coe', 'cse', 'diag'):
            with pytest.raises(ValueError):
                real_algebra.random_unitary(4, measure=measure)

    def test_batched_shape(self, basic_algebra):
        U = basic_algebra.random_unitary(3, measure='haar', batch_size=7)
        assert U.shape == (7, 3, 3)
        expected = torch.eye(3, dtype=U.dtype).expand(7, 3, 3)
        assert torch.allclose(U @ U.conj().transpose(-2, -1), expected, atol=1e-5)

    def test_haar_is_circular(self, basic_algebra):
        # For Haar U(n): E|U_11|^2 = 1/n and E[U_11^2] = 0.
        torch.manual_seed(0)
        n = 4
        U = basic_algebra.random_unitary(n, measure='haar', batch_size=40000)
        u = U[:, 0, 0]
        assert abs((u.abs() ** 2).mean().item() - 1 / n) < 0.01
        assert (u ** 2).mean().abs().item() < 0.01

    def test_haar_su_batched_det_one(self, basic_algebra):
        U = basic_algebra.random_unitary(3, measure='haar_su', batch_size=50)
        det = torch.linalg.det(U)
        assert torch.allclose(det, torch.ones_like(det), atol=1e-5)

    def test_haar_so_real_det_one(self, real_algebra):
        U = real_algebra.random_unitary(3, measure='haar_su', batch_size=50)
        det = torch.linalg.det(U)
        assert torch.allclose(det, torch.ones_like(det), atol=1e-5)

    def test_diag_unitary(self, basic_algebra):
        U = basic_algebra.random_unitary(4, measure='diag')
        expected = torch.eye(4, dtype=U.dtype, device=U.device)
        assert torch.allclose(U @ U.conj().T, expected, atol=1e-6)
        assert torch.allclose(U - torch.diag(torch.diag(U)), torch.zeros(4, 4, dtype=U.dtype), atol=1e-6)


# ============================================================================
# Norms
# ============================================================================

class TestNorms:
    def test_trace_norm_identity(self, basic_algebra):
        op = basic_algebra.identity(batch_size=2)
        expected = torch.tensor([5.0, 5.0])
        assert torch.allclose(op.trace_norm(), expected, atol=1e-6)

    def test_frobenius_norm_identity(self, basic_algebra):
        op = basic_algebra.identity(batch_size=2)
        expected = torch.full((2,), torch.sqrt(torch.tensor(5.0)))
        assert torch.allclose(op.frobenius_norm(), expected, atol=1e-6)

    def test_operator_norm_identity(self, basic_algebra):
        op = basic_algebra.identity(batch_size=2)
        assert torch.allclose(op.operator_norm(), torch.tensor(1.0))


# ============================================================================
# Tests for lambda_max and lambda_min (global extremes)
# ============================================================================

class TestLambdaMaxMin:
    @pytest.fixture
    def single_channel_positive(self):
        n_factors = [5]
        k_factors = [5]
        alg = TypeIAlgebra(n_factors, k_factors, batch_size=3, complex_valued=False, device='cpu')
        return alg

    def test_lambda_max_single_channel(self, single_channel_positive):
        def sampler(dim):
            eig = torch.rand(3, dim)
            max_vals = eig.max(dim=-1, keepdim=True)[0]
            eig = eig / (max_vals + 1e-12)
            return eig
        op = single_channel_positive.operator_from_eigenvalues(sampler, batch_size=3,
                                                                force_positive=True, force_self_adjoint=True)
        _ = op.matrix
        all_eig = torch.cat([eig for eig in op._eigenvalues], dim=1)
        true_max = all_eig.max(dim=1)[0]
        computed_max = op.lambda_max
        assert torch.allclose(computed_max, true_max, atol=1e-4)

    def test_lambda_min_single_channel(self, single_channel_positive):
        def sampler(dim):
            eig = torch.rand(3, dim)
            max_vals = eig.max(dim=-1, keepdim=True)[0]
            eig = eig / (max_vals + 1e-12)
            return eig
        op = single_channel_positive.operator_from_eigenvalues(sampler, batch_size=3,
                                                                force_positive=True, force_self_adjoint=True)
        _ = op.matrix
        all_eig = torch.cat([eig for eig in op._eigenvalues], dim=1)
        true_min = all_eig.min(dim=1)[0]
        computed_min = op.lambda_min
        assert torch.allclose(computed_min, true_min, atol=1e-4)

    def test_lambda_max_multi_channel(self, basic_algebra):
        def sampler(dim):
            eig = torch.rand(2, dim) + 0j
            max_vals = eig.real.max(dim=-1, keepdim=True)[0]
            eig = eig / (max_vals + 1e-12)
            return eig
        op = basic_algebra.operator_from_eigenvalues(sampler, batch_size=2,
                                                      force_positive=True, force_self_adjoint=True)
        _ = op.matrix
        all_eig = torch.cat([eig for eig in op._eigenvalues], dim=1)  # здесь eig уже комплексные
        true_max = all_eig.real.max(dim=1)[0]   # собственные значения вещественные, но хранятся как complex
        computed_max = op.lambda_max
        assert torch.allclose(computed_max, true_max, atol=1e-4)

    def test_lambda_min_multi_channel(self, basic_algebra):
        def sampler(dim):
            eig = torch.rand(2, dim) + 0j
            max_vals = eig.real.max(dim=-1, keepdim=True)[0]
            eig = eig / (max_vals + 1e-12)
            return eig
        op = basic_algebra.operator_from_eigenvalues(sampler, batch_size=2,
                                                      force_positive=True, force_self_adjoint=True)
        _ = op.matrix
        all_eig = torch.cat([eig for eig in op._eigenvalues], dim=1)
        true_min = all_eig.real.min(dim=1)[0]
        computed_min = op.lambda_min
        assert torch.allclose(computed_min, true_min, atol=1e-4)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])

# ============================================================================
# Regression tests
# ============================================================================

class TestRegressions:
    @pytest.mark.parametrize("k", [8, 300])
    def test_lambda_min_negative_spectrum(self, k):
        # k = 300: blocks above the former power-iteration threshold (now exact for every size)
        alg = TypeIAlgebra([k], [k], complex_valued=False, device='cpu')
        op = alg.operator_from_eigenvalues(lambda d: torch.linspace(-2.0, 1.0, d),
                                           force_self_adjoint=True)
        assert abs(op.lambda_max.item() - 1.0) < 1e-3
        assert abs(op.lambda_min.item() + 2.0) < 1e-3

    def test_large_batch_validation_does_not_fail(self):
        # property checks must not accumulate round-off over the whole batch
        alg = TypeIAlgebra([16] * 4, [16] * 4, complex_valued=False, device='cpu')
        op = alg.operator_from_eigenvalues(lambda d: torch.rand(5000, d), batch_size=5000,
                                           force_positive=True, force_self_adjoint=True)
        assert op.matrix.shape == (5000, 4, 16, 16)

    def test_michelson_contrast_global(self, basic_algebra):
        op = basic_algebra.operator_from_eigenvalues(
            lambda d: torch.tensor([1.0, 3.0]) if d == 2 else torch.tensor([2.0, 2.0, 2.0]),
            batch_size=2, force_positive=True)
        expected = torch.full((2,), (3.0 - 1.0) / (3.0 + 1.0))
        assert torch.allclose(op.michelson_contrast, expected, atol=1e-5)

    def test_identity_is_invertible(self, basic_algebra):
        assert basic_algebra.identity().is_invertible

    def test_unitary_is_invertible(self, basic_algebra):
        U = basic_algebra.random_unitary(3, batch_size=2)
        mat = torch.zeros(2, 2, 3, 3, dtype=U.dtype)
        mat[:, 0, :2, :2] = basic_algebra.random_unitary(2, batch_size=2)
        mat[:, 1] = U
        op = TypeIAlgebra.Operator(basic_algebra, matrix=mat)
        assert op.is_invertible


class TestPrecision:
    @pytest.mark.parametrize("cplx", [True, False])
    def test_double_precision_end_to_end(self, cplx):
        alg = TypeIAlgebra([3, 4], [3, 4], complex_valued=cplx, precision='double', device='cpu')
        expected = torch.complex128 if cplx else torch.float64
        assert alg.random_unitary(4).dtype == expected
        op = alg.operator_from_eigenvalues(lambda d: torch.rand(2, d), batch_size=2, force_positive=True)
        assert op.matrix.dtype == expected
        assert op.lambda_max.dtype == torch.float64
        U = alg.random_unitary(6, batch_size=10)
        eye = torch.eye(6, dtype=expected).expand(10, 6, 6)
        assert (U @ U.conj().transpose(-2, -1) - eye).abs().max() < 1e-12

    def test_invalid_precision(self):
        with pytest.raises(ValueError):
            TypeIAlgebra([2], [2], precision='half')


class TestNumerics:
    def test_like_and_cast(self):
        a = TypeIAlgebra([3, 4], [3, 4], complex_valued=False, device='cpu')
        X = a.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(4, d), batch_size=4, force_positive=True)
        b = a.like(complex_valued=True, precision='double')
        assert b.k_factors == a.k_factors and b.dtype == torch.complex128
        Y = X.cast(b)
        assert Y.algebra is b and Y.matrix.dtype == torch.complex128
        assert torch.allclose(Y.matrix.real.float(), X.matrix)
        assert torch.allclose(Y.michelson_contrast.float(), X.michelson_contrast, atol=1e-5)
        assert a.bytes_per_operator(10) * 4 == b.bytes_per_operator(10)

    def test_cast_to_real_requires_real_matrix(self):
        b = TypeIAlgebra([2], [2], complex_valued=True, device='cpu')
        U = b.operator(b.random_unitary(2, batch_size=3)[:, None])
        with pytest.raises(ValueError):
            U.cast(b.like(complex_valued=False))


class TestEigenvalues:
    def test_matches_eigh_and_is_differentiable(self):
        alg = TypeIAlgebra([2, 3, 4], [2, 3, 4], precision='double')
        A = alg.random_unitary_operator(5)
        H = A + A.adjoint()
        for w, (w2, _) in zip(H.eigenvalues(), H.eigh()):
            assert torch.allclose(w, w2, atol=1e-12)
        m = H.matrix.clone().requires_grad_()
        S = sum((w ** 2).sum() for w in alg.operator(m).eigenvalues())
        S.backward()                                       # sum of squared eigenvalues = ||H||_F^2
        assert torch.allclose(m.grad, 2 * H.matrix, atol=1e-10)

    def test_requires_self_adjoint(self):
        alg = TypeIAlgebra([2], [2])
        U = alg.random_unitary_operator(1)
        with pytest.raises(RuntimeError):
            U.eigenvalues()


class TestSecondReview:
    """Regressions found in the second review round."""

    def test_extremes_exact_for_large_and_symmetric_blocks(self):
        alg = TypeIAlgebra([4, 3], [4, 3], precision='double')
        alg.exact_eig_max_dim = 1                      # formerly switched to power iteration
        M = torch.zeros(1, 2, 4, 4, dtype=torch.complex128)
        M[0, 0] = torch.diag(torch.tensor([1.0, -1.0, 0.5, 0.2], dtype=torch.complex128))
        M[0, 1, :3, :3] = torch.diag(torch.tensor([0.3, 0.1, -0.2], dtype=torch.complex128))
        op = alg.operator(M, is_self_adjoint=True)
        assert abs(op.lambda_max.item() - 1) < 1e-12 and abs(op.lambda_min.item() + 1) < 1e-12

    def test_operator_norm_exact(self):
        alg = TypeIAlgebra([5, 3], [5, 3], precision='double')
        A = alg.operator(torch.randn(2, 2, 5, 5, dtype=torch.complex128) * torch.tensor([1.0, 0.0])[None, :, None, None]
                         + alg.zero(2).matrix)
        ref = torch.stack([torch.linalg.matrix_norm(A.matrix[b, 0], ord=2) for b in range(2)])
        assert torch.allclose(A.operator_norm(), ref)

    def test_flags_for_complex_spectra(self):
        alg = TypeIAlgebra([2], [2], precision='double')
        A = alg.operator_from_eigenvalues(lambda d: torch.tensor([1j, 1.0]))
        assert A.is_normal and not A.is_positive and not A.is_projection and not A.is_self_adjoint

    def test_tau_vN_normalised_with_empty_sector(self):
        alg = TypeIAlgebra([2, 1, 3], [2, 0, 3], precision='double')
        assert abs(alg.identity().tau_vN().real.item() - 1) < 1e-12

    def test_entropy_ignores_padding(self):
        from torch_vn_algebra import states
        alg = TypeIAlgebra([1] * 50 + [300], [1] * 50 + [300], precision='double')
        m = torch.zeros(1, 51, 300, 300, dtype=torch.complex128)
        m[0, 50, :2, :2] = torch.diag(torch.tensor([0.3, 0.7], dtype=torch.complex128))
        rho = alg.operator(m, is_self_adjoint=True, is_positive=True)
        exact = -(0.3 * math.log(0.3) + 0.7 * math.log(0.7))
        assert abs(rho.entropy().item() - exact) < 1e-12
        assert abs(states.von_neumann_entropy(rho).item() - exact) < 1e-12

    def test_hs_measure_per_block(self):
        from torch_vn_algebra import states
        alg = TypeIAlgebra([2, 5], [2, 5], precision='double')
        rho = states.random_density_matrix(alg, batch_size=20000)
        b = rho.matrix[:, 0, :2, :2]
        b = b / torch.diagonal(b, dim1=-2, dim2=-1).sum(-1)[:, None, None]
        purity = torch.einsum('bij,bji->b', b, b).real.mean().item()
        assert abs(purity - 0.8) < 0.01                 # HS on M_2: 2k/(k^2+1)

    def test_mixed_unitary_channel_trace_preserving_in_double(self):
        from torch_vn_algebra import channels
        alg = TypeIAlgebra([2, 3], [2, 3], precision='double')
        assert channels.random_mixed_unitary_channel(alg, 5, 4).is_trace_preserving()


class TestThirdReview:
    """Regressions from the third review round (user-level API)."""

    def test_per_sample_scalar(self):
        alg = TypeIAlgebra([2], [2], precision='double')
        X = alg.identity(2) * torch.tensor([1.0, 3.0], dtype=torch.float64)
        assert torch.allclose(X.matrix[0, 0], torch.eye(2, dtype=X.matrix.dtype))
        assert torch.allclose(X.matrix[1, 0], 3 * torch.eye(2, dtype=X.matrix.dtype))
        assert X.matrix.dtype == alg.dtype

    def test_pinv_relative_tolerance_single(self):
        alg = TypeIAlgebra([3], [3])
        U = alg.random_unitary_operator(1)
        P = alg.from_blocks([torch.diag(torch.tensor([1.0, 0.0, 0.0]))])
        A = U @ P @ U.adjoint()
        assert A.inverse().matrix.abs().max().item() < 10

    def test_entropy_rejects_non_positive(self):
        alg = TypeIAlgebra([2], [2], precision='double')
        A = alg.operator(torch.diag(torch.tensor([-0.5, 0.5], dtype=torch.float64))[None, None].to(alg.dtype))
        with pytest.raises(RuntimeError):
            A.entropy()

    def test_real_algebra_refuses_complex_data(self):
        alg = TypeIAlgebra([2], [2], complex_valued=False)
        with pytest.raises(ValueError):
            alg.operator(torch.ones(1, 1, 2, 2, dtype=torch.complex64) * 1j)
        with pytest.raises(ValueError):
            alg.operator_from_eigenvalues(lambda d: torch.tensor([1 + 1j, 2 - 1j]))
        assert alg.operator(torch.ones(1, 1, 2, 2, dtype=torch.float64)).matrix.dtype == alg.dtype

    def test_clear_errors_between_algebras(self):
        a, b = TypeIAlgebra([2], [2]), TypeIAlgebra([3], [3])
        with pytest.raises(ValueError):
            a.identity() + b.identity()
        with pytest.raises(TypeError):
            a.identity() + 1.0
        assert torch.allclose((a.identity() + a.like().identity()).matrix, 2 * a.identity().matrix)

    def test_mix_checks_weights(self):
        from torch_vn_algebra import states
        alg = TypeIAlgebra([2], [2], precision='double')
        rho, sig = states.random_density_matrix(alg, 2), states.random_density_matrix(alg, 2)
        with pytest.raises(ValueError):
            rho.mix(sig, 2.0)
        m = rho.mix(sig, torch.tensor([0.0, 1.0]))
        assert torch.allclose(m.matrix[0], rho.matrix[0]) and torch.allclose(m.matrix[1], sig.matrix[1])

    def test_spacing_ratio_float32_degenerate(self):
        from torch_vn_algebra import SpinChain
        r = SpinChain.level_spacing_ratio(torch.tensor([0.0, 0.0, 0.0, 1.0, 2.0]))
        assert torch.isfinite(r)

    def test_channel_errors_and_mixed_composition(self):
        from torch_vn_algebra import channels, SpinChain
        alg = TypeIAlgebra([2], [2])
        Phi = channels.random_channel(alg, 2)
        with pytest.raises(TypeError):
            Phi(torch.eye(2))
        ch = SpinChain(3)
        L = ch.lowering(0)                                     # InterSectorChannel on the chain algebra
        Id = channels.identity_channel(ch.algebra)
        assert (Id @ L).kraus.shape[-1] == L.kraus.shape[-1]


def test_hilbert_space_inner_product_and_haar():
    from torch_vn_algebra.hilbert_space import HilbertSpace
    H = HilbertSpace(n=3, k=2, device='cpu', precision='double')
    b = H.random_basis()
    assert abs(H.inner_product(b.bra(0), b.ket(0)).item() - 1) < 1e-12
    assert abs(H.inner_product(b.bra(0), b.ket(1)).item()) < 1e-12
    assert abs(H.inner_product(b.ket(0), b.ket(0)).item() - 1) < 1e-12        # kets in: conjugated inside
    Hb = HilbertSpace(n=3, k=3, batch_size=20000, device='cpu', precision='double')
    for V in (Hb.random_basis().V[:, 0], Hb.haar_basis().V[:, 0]):
        assert abs((V[:, 0, 0].abs() ** 4).mean().item() - 2 / 12) < 0.01    # Haar U(3) moment
    assert Hb.random_basis().random_subspace_vector().dtype == torch.complex128


def test_central_with_complex_scalars():
    alg = TypeIAlgebra([2, 3], [2, 3], precision='double')
    Z = alg.central([1j, 2.0])
    assert Z.is_normal and not Z.is_self_adjoint and not Z.is_positive
    assert torch.allclose(Z.matrix[0, 0, :2, :2], 1j * torch.eye(2, dtype=Z.matrix.dtype))
    assert alg.central([0.0, 1.0]).is_projection
    with pytest.raises(ValueError):
        TypeIAlgebra([2], [2], complex_valued=False).central([1j])
    with pytest.raises(ValueError):
        alg.central([1.0])


class TestFourthReview:
    """Regressions from the fourth review round."""

    @pytest.mark.parametrize('complex_valued', [False, True])
    def test_blockwise_functions_gradcheck_with_padding(self, complex_valued):
        alg = TypeIAlgebra([2, 3], [2, 3], complex_valued=complex_valued, precision='double')
        torch.manual_seed(0)
        U = alg.random_unitary_operator(1)
        D = alg.from_blocks([torch.diag(torch.tensor(v, dtype=torch.float64)) for v in ([1.0, 2.0], [0.7, 1.5, 3.0])])
        P = (U @ D @ U.adjoint()).matrix   # positive definite on active blocks

        for f in (lambda X: X.inverse(), lambda X: X.sqrt(), lambda X: X.abs()):
            def g(m):
                m = 0.5 * (m + m.conj().transpose(-2, -1))
                return f(alg.operator(m)).matrix
            x = P.clone().requires_grad_(True)
            assert torch.autograd.gradcheck(g, (x,), eps=1e-6, atol=1e-5)
            out = g(P)
            assert out[0, 0, 2:, :].abs().max() == 0 and out[0, 0, :, 2:].abs().max() == 0   # padding stays zero

    def test_structure_includes_charges(self):
        a = TypeIAlgebra([2, 2], [2, 2], charges=[0, 1])
        b = TypeIAlgebra([2, 2], [2, 2], charges=[1, 2])
        with pytest.raises(ValueError):
            a.identity() + b.identity()

    def test_complex_tensor_scalar_rejected_in_real_algebra(self):
        alg = TypeIAlgebra([2], [2], complex_valued=False)
        with pytest.raises(ValueError):
            alg.identity() * torch.tensor(1j)

    def test_central_keeps_gradients_and_per_sample_values(self):
        alg = TypeIAlgebra([2, 3], [2, 3], precision='double')
        s = torch.tensor([1.0, 2.0], dtype=torch.float64, requires_grad=True)
        Z = alg.central([s, 3.0])
        assert Z.matrix.shape[0] == 2 and Z.is_positive and Z.is_self_adjoint
        assert torch.allclose(Z.matrix[1, 0, :2, :2], 2 * torch.eye(2, dtype=Z.matrix.dtype))
        Z.matrix.real.sum().backward()
        assert torch.allclose(s.grad, torch.tensor([2.0, 2.0], dtype=torch.float64))

    def test_intersector_from_blocks_checks_field(self):
        from torch_vn_algebra.channels import InterSectorChannel
        alg = TypeIAlgebra([1, 2], [1, 2], complex_valued=False)
        with pytest.raises(ValueError):
            InterSectorChannel.from_blocks(alg, alg, {(0, 1): [torch.ones(1, 2) * 1j]})

    def test_inner_product_n1(self):
        from torch_vn_algebra.hilbert_space import HilbertSpace
        H = HilbertSpace(n=1, k=1, device='cpu', precision='double')
        v = torch.full((1, 1, 1, 1), 1j, dtype=torch.complex128)
        with pytest.raises(ValueError):
            H.inner_product(v, v)
        assert abs(H.inner_product(v, v, bra_is_row=False).item() - 1) < 1e-12
        assert abs(H.inner_product(v, v, bra_is_row=True).item() + 1) < 1e-12

    def test_lindblad_broadcasts_intersector_and_callable(self):
        from torch_vn_algebra import SpinChain, dynamics, states
        ch = SpinChain(2)
        alg = ch.algebra
        rho = states.maximally_mixed_state(alg)
        L = ch.lowering(0)
        L4 = type(L)(L.algebra_in, L.algebra_out, L.kraus.expand(4, *L.kraus.shape[1:]).clone())
        out = dynamics.lindblad_evolve(rho, None, [L4], times=[0.0, 0.1])
        assert out[-1].matrix.shape[0] == 4
        U = alg.random_unitary_operator(3)
        H3 = U + U.adjoint()
        out = dynamics.lindblad_evolve(rho, lambda t: H3, [], times=[0.0, 0.1])
        assert out[-1].matrix.shape[0] == 3

    def test_mix_checks_algebra_and_scalar_p(self):
        from torch_vn_algebra import states, channels
        a, b = TypeIAlgebra([2], [2]), TypeIAlgebra([3], [3])
        with pytest.raises(ValueError):
            states.maximally_mixed_state(a).mix(states.maximally_mixed_state(b), 0.5)
        Phi = channels.random_channel(a, 2)
        with pytest.raises(ValueError):
            Phi.mix(Phi, torch.tensor([0.1, 0.2]))
        assert Phi.mix(Phi, torch.tensor(0.3)).kraus.shape[1] == 4
