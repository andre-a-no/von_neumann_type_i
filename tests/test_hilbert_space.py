"""
Tests for Hilbert space module.

Run with: pytest tests/test_hilbert_space.py -v
"""

import torch
import pytest

from torch_vn_algebra.hilbert_space import HilbertSpace, tensor_product_hilbert, direct_sum_hilbert


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def basic_hilbert():
    """Basic Hilbert space with n=5, k=3, complex, CPU."""
    return HilbertSpace(n=5, k=3, batch_size=2, n_channels=2, complex_valued=True, device='cpu')


@pytest.fixture
def real_hilbert():
    """Real Hilbert space, CPU."""
    return HilbertSpace(n=4, k=2, batch_size=1, n_channels=1, complex_valued=False, device='cpu')


# ============================================================================
# Tests for HilbertSpace initialization
# ============================================================================

class TestHilbertSpaceInit:
    def test_valid_initialization(self):
        H = HilbertSpace(n=10, k=5, batch_size=4, n_channels=3, complex_valued=True, device='cpu')
        assert H.n == 10 and H.k == 5 and H.batch_size == 4 and H.n_channels == 3
        assert H.complex_valued is True
        assert H.dtype == torch.complex64

    def test_real_initialization(self):
        H = HilbertSpace(n=10, k=5, complex_valued=False, device='cpu')
        assert H.complex_valued is False and H.dtype == torch.float32

    def test_k_equals_n(self):
        H = HilbertSpace(n=5, k=5, device='cpu')
        assert H.k == H.n

    def test_k_less_than_n(self):
        H = HilbertSpace(n=10, k=3, device='cpu')
        assert H.k < H.n

    def test_invalid_k_greater_than_n(self):
        with pytest.raises(AssertionError):
            HilbertSpace(n=5, k=10, device='cpu')

    def test_invalid_n_zero(self):
        with pytest.raises(AssertionError):
            HilbertSpace(n=0, k=0, device='cpu')

    def test_invalid_k_negative(self):
        with pytest.raises(AssertionError):
            HilbertSpace(n=5, k=-1, device='cpu')

    def test_device_auto_detection(self):
        H = HilbertSpace(n=5, k=3)
        assert H.device is not None

    def test_device_explicit(self):
        device = torch.device('cpu')
        H = HilbertSpace(n=5, k=3, device=device)
        assert H.device == device


# ============================================================================
# Tests for HilbertSpace properties
# ============================================================================

class TestHilbertSpaceProperties:
    def test_dim_property(self, basic_hilbert):
        assert basic_hilbert.dim == 5

    def test_subspace_dim_property(self, basic_hilbert):
        assert basic_hilbert.subspace_dim == 3

    def test_shape_ambient(self, basic_hilbert):
        assert basic_hilbert.shape_ambient == (2, 2, 5, 5)

    def test_shape_ket(self, basic_hilbert):
        assert basic_hilbert.shape_ket == (2, 2, 5, 1)

    def test_shape_bra(self, basic_hilbert):
        assert basic_hilbert.shape_bra == (2, 2, 1, 5)

    def test_shape_ambient_real(self, real_hilbert):
        assert real_hilbert.shape_ambient == (1, 1, 4, 4)


# ============================================================================
# Tests for inner product and norm
# ============================================================================

class TestInnerProductAndNorm:
    def test_inner_product_complex(self, basic_hilbert):
        bra = torch.randn(2, 2, 1, 5, dtype=torch.complex64, device=basic_hilbert.device)
        ket = torch.randn(2, 2, 5, 1, dtype=torch.complex64, device=basic_hilbert.device)
        inner = basic_hilbert.inner_product(bra, ket)
        assert inner.shape == (2, 1, 1, 1)
        assert torch.is_complex(inner)

    def test_inner_product_real(self, real_hilbert):
        bra = torch.randn(1, 1, 1, 4, device=real_hilbert.device)
        ket = torch.randn(1, 1, 4, 1, device=real_hilbert.device)
        inner = real_hilbert.inner_product(bra, ket)
        assert inner.shape == (1, 1, 1, 1)
        assert not torch.is_complex(inner)

    def test_inner_product_orthogonal_vectors(self, basic_hilbert):
        ket1 = torch.zeros(2, 2, 5, 1, dtype=torch.complex64, device=basic_hilbert.device)
        ket2 = torch.zeros(2, 2, 5, 1, dtype=torch.complex64, device=basic_hilbert.device)
        ket1[..., 0, 0] = 1.0
        ket2[..., 1, 0] = 1.0
        bra1 = ket1.conj().transpose(-2, -1)
        inner = basic_hilbert.inner_product(bra1, ket2)
        assert torch.allclose(inner.real, torch.zeros(2, 1, 1, 1, device=basic_hilbert.device))
        assert torch.allclose(inner.imag, torch.zeros(2, 1, 1, 1, device=basic_hilbert.device))

    def test_norm(self, basic_hilbert):
        ket = torch.zeros(2, 2, 5, 1, dtype=torch.complex64, device=basic_hilbert.device)
        # Устанавливаем значение только в первый канал (channel 0)
        ket[:, 0, 0, 0] = 3.0
        norm = basic_hilbert.norm(ket)
        expected = torch.tensor(3.0, device=basic_hilbert.device)
        assert torch.allclose(norm.real, expected)
        assert norm.shape == (2, 1, 1, 1)

    def test_normalize(self, basic_hilbert):
        ket = torch.randn(2, 2, 5, 1, dtype=torch.complex64, device=basic_hilbert.device)
        normalized = basic_hilbert.normalize(ket)
        norm = basic_hilbert.norm(normalized)
        assert torch.all(torch.isfinite(norm))
        expected = torch.ones(2, 1, 1, 1, device=basic_hilbert.device)
        assert torch.allclose(norm.real, expected, atol=1e-6)


# ============================================================================
# Tests for Basis class
# ============================================================================

class TestBasis:
    def test_standard_basis_creation(self, basic_hilbert):
        basis = basic_hilbert.standard_basis()
        assert basis.V.shape == (2, 2, 5, 3)
        assert basis.is_orthonormal()

    def test_random_basis_creation(self, basic_hilbert):
        basis = basic_hilbert.random_basis()
        assert basis.V.shape == (2, 2, 5, 3)
        assert basis.is_orthonormal()

    def test_haar_basis_creation(self, basic_hilbert):
        basis = basic_hilbert.haar_basis()
        assert basis.V.shape == (2, 2, 5, 3)
        assert basis.is_orthonormal()

    def test_basis_with_explicit_V(self, basic_hilbert):
        V = torch.zeros(2, 2, 5, 3, dtype=torch.complex64, device=basic_hilbert.device)
        V[..., 0, 0] = 1.0
        V[..., 1, 1] = 1.0
        V[..., 2, 2] = 1.0
        basis = basic_hilbert.Basis(basic_hilbert, V=V)
        assert basis.is_orthonormal()

    def test_projection_property(self, basic_hilbert):
        basis = basic_hilbert.standard_basis()
        P = basis.projection
        assert P.shape == (2, 2, 5, 5)
        assert torch.allclose(P @ P, P, atol=1e-6)
        P_adj = P.conj().transpose(-2, -1)
        assert torch.allclose(P, P_adj, atol=1e-6)

    def test_restriction_property(self, basic_hilbert):
        basis = basic_hilbert.standard_basis()
        R = basis.restriction
        assert R.shape == (2, 2, 3, 5)

    def test_embed_vector(self, basic_hilbert):
        basis = basic_hilbert.standard_basis()
        v_sub = torch.randn(2, 2, 3, 1, dtype=torch.complex64, device=basic_hilbert.device)
        v_amb = basis.embed_vector(v_sub)
        assert v_amb.shape == (2, 2, 5, 1)
        assert torch.allclose(v_amb[..., :3, :], v_sub, atol=1e-6)
        zeros_last = torch.zeros(2, 2, 2, 1, dtype=torch.complex64, device=basic_hilbert.device)
        assert torch.allclose(v_amb[..., 3:, :], zeros_last, atol=1e-6)

    def test_restrict_vector(self, basic_hilbert):
        basis = basic_hilbert.standard_basis()
        v_amb = torch.randn(2, 2, 5, 1, dtype=torch.complex64, device=basic_hilbert.device)
        v_sub = basis.restrict_vector(v_amb)
        assert v_sub.shape == (2, 2, 3, 1)
        assert torch.allclose(v_sub, v_amb[..., :3, :], atol=1e-6)

    def test_embed_and_restrict_roundtrip(self, basic_hilbert):
        basis = basic_hilbert.standard_basis()
        v_sub = torch.randn(2, 2, 3, 1, dtype=torch.complex64, device=basic_hilbert.device)
        v_amb = basis.embed_vector(v_sub)
        v_sub_round = basis.restrict_vector(v_amb)
        assert torch.allclose(v_sub, v_sub_round, atol=1e-6)

    def test_embed_operator(self, basic_hilbert):
        basis = basic_hilbert.standard_basis()
        A_sub = torch.randn(2, 2, 3, 3, dtype=torch.complex64, device=basic_hilbert.device)
        A_amb = basis.embed_operator(A_sub)
        assert A_amb.shape == (2, 2, 5, 5)
        assert torch.allclose(A_amb[..., :3, :3], A_sub, atol=1e-6)
        zeros_bottom = torch.zeros(2, 2, 2, 5, dtype=torch.complex64, device=basic_hilbert.device)
        zeros_right = torch.zeros(2, 2, 5, 2, dtype=torch.complex64, device=basic_hilbert.device)
        assert torch.allclose(A_amb[..., 3:, :], zeros_bottom, atol=1e-6)
        assert torch.allclose(A_amb[..., :, 3:], zeros_right, atol=1e-6)

    def test_restrict_operator(self, basic_hilbert):
        basis = basic_hilbert.standard_basis()
        A_amb = torch.randn(2, 2, 5, 5, dtype=torch.complex64, device=basic_hilbert.device)
        A_sub = basis.restrict_operator(A_amb)
        assert A_sub.shape == (2, 2, 3, 3)

    def test_bra_ket_methods(self, basic_hilbert):
        basis = basic_hilbert.standard_basis()
        for i in range(3):
            ket_i = basis.ket(i)
            bra_i = basis.bra(i)
            assert ket_i.shape == (2, 2, 5, 1)
            assert bra_i.shape == (2, 2, 1, 5)
            for j in range(3):
                inner = torch.matmul(basis.bra(i), basis.ket(j))
                expected = 1.0 if i == j else 0.0
                assert torch.allclose(inner.real, torch.tensor(expected, device=basic_hilbert.device))
                assert torch.allclose(inner.imag, torch.tensor(0.0, device=basic_hilbert.device))

    def test_outer_product(self, basic_hilbert):
        basis = basic_hilbert.standard_basis()
        outer = basis.outer_product(0, 1)
        assert outer.shape == (2, 2, 5, 5)
        expected = torch.ones(2, 2, device=basic_hilbert.device)
        assert torch.allclose(outer[..., 0, 1].real, expected)
        assert torch.allclose(outer[..., 0, 1].imag, torch.zeros(2, 2, device=basic_hilbert.device))

    def test_gram_matrix(self, basic_hilbert):
        basis = basic_hilbert.standard_basis()
        gram = basis.gram_matrix()
        identity = torch.eye(3, dtype=torch.complex64, device=basic_hilbert.device)
        identity_exp = identity.unsqueeze(0).unsqueeze(0).expand(2, 2, 3, 3)
        assert torch.allclose(gram, identity_exp, atol=1e-6)

    def test_orthonormalize(self, basic_hilbert):
        V = torch.randn(2, 2, 5, 3, dtype=torch.complex64, device=basic_hilbert.device)
        basis = basic_hilbert.Basis(basic_hilbert, V=V)
        assert not basis.is_orthonormal(tol=1e-4)
        basis.orthonormalize()
        assert basis.is_orthonormal()

    def test_random_subspace_vector(self, basic_hilbert):
        basis = basic_hilbert.standard_basis()
        v = basis.random_subspace_vector(normalize=True)
        assert v.shape == (2, 2, 5, 1)
        norm = basic_hilbert.norm(v)
        expected = torch.ones(2, 1, 1, 1, device=basic_hilbert.device)
        assert torch.all(torch.isfinite(norm))
        assert torch.allclose(norm.real, expected, atol=1e-6)
        zeros_last = torch.zeros(2, 2, 2, 1, dtype=torch.complex64, device=basic_hilbert.device)
        assert torch.allclose(v[..., 3:, :], zeros_last, atol=1e-6)


# ============================================================================
# Tests for vector operations
# ============================================================================

class TestVectorOperations:
    def test_embed_vector_with_standard_basis(self, basic_hilbert):
        basis = basic_hilbert.standard_basis()
        v_sub = torch.tensor([[[[1.0], [2.0], [3.0]]]], dtype=torch.complex64, device=basic_hilbert.device)
        v_sub = v_sub.expand(2, 2, 3, 1)
        v_amb = basis.embed_vector(v_sub)
        expected_ones = torch.ones(2, 2, device=basic_hilbert.device)
        expected_two = torch.tensor(2.0, device=basic_hilbert.device)
        expected_three = torch.tensor(3.0, device=basic_hilbert.device)
        assert torch.allclose(v_amb[..., 0, 0].real, expected_ones)
        assert torch.allclose(v_amb[..., 1, 0].real, expected_two)
        assert torch.allclose(v_amb[..., 2, 0].real, expected_three)
        zeros_last = torch.zeros(2, 2, 2, 1, dtype=torch.complex64, device=basic_hilbert.device)
        assert torch.allclose(v_amb[..., 3:, :], zeros_last, atol=1e-6)

    def test_restrict_vector_with_standard_basis(self, basic_hilbert):
        basis = basic_hilbert.standard_basis()
        v_amb = torch.zeros(2, 2, 5, 1, dtype=torch.complex64, device=basic_hilbert.device)
        v_amb[..., 0, 0] = 1.0
        v_amb[..., 1, 0] = 2.0
        v_amb[..., 2, 0] = 3.0
        v_sub = basis.restrict_vector(v_amb)
        expected_ones = torch.ones(2, 2, device=basic_hilbert.device)
        expected_two = torch.tensor(2.0, device=basic_hilbert.device)
        expected_three = torch.tensor(3.0, device=basic_hilbert.device)
        assert torch.allclose(v_sub[..., 0, 0].real, expected_ones)
        assert torch.allclose(v_sub[..., 1, 0].real, expected_two)
        assert torch.allclose(v_sub[..., 2, 0].real, expected_three)


# ============================================================================
# Tests for composite Hilbert spaces
# ============================================================================

class TestCompositeSpaces:
    def test_tensor_product_hilbert(self):
        H1 = HilbertSpace(n=2, k=2, batch_size=1, n_channels=1, device='cpu')
        H2 = HilbertSpace(n=3, k=2, batch_size=1, n_channels=1, device='cpu')
        H_prod = tensor_product_hilbert([H1, H2])
        assert H_prod.n == 6 and H_prod.k == 4

    def test_direct_sum_hilbert(self):
        H1 = HilbertSpace(n=2, k=2, batch_size=1, n_channels=1, device='cpu')
        H2 = HilbertSpace(n=3, k=2, batch_size=1, n_channels=1, device='cpu')
        H_sum = direct_sum_hilbert([H1, H2])
        assert H_sum.n == 5 and H_sum.k == 4


# ============================================================================
# Tests for real Hilbert space
# ============================================================================

class TestRealHilbertSpace:
    def test_real_inner_product(self, real_hilbert):
        bra = torch.ones(1, 1, 1, 4, device=real_hilbert.device)
        ket = torch.ones(1, 1, 4, 1, device=real_hilbert.device)
        inner = real_hilbert.inner_product(bra, ket)
        assert torch.allclose(inner, torch.tensor(4.0, device=real_hilbert.device))

    def test_real_basis(self, real_hilbert):
        basis = real_hilbert.standard_basis()
        assert basis.V.dtype == torch.float32
        assert basis.is_orthonormal()


# ============================================================================
# Tests for edge cases
# ============================================================================

class TestEdgeCases:
    def test_k_zero(self):
        H = HilbertSpace(n=5, k=0, device='cpu')
        assert H.k == 0
        basis = H.standard_basis()
        assert basis.V.shape == (1, 1, 5, 0)

    def test_k_equals_n_standard_basis(self):
        H = HilbertSpace(n=3, k=3, complex_valued=True, device='cpu')
        basis = H.standard_basis()
        identity = torch.eye(3, dtype=torch.complex64, device=H.device)
        assert torch.allclose(basis.projection, identity)

    def test_batch_size_one(self):
        H = HilbertSpace(n=5, k=3, batch_size=1, n_channels=1, device='cpu')
        basis = H.standard_basis()
        assert basis.V.shape == (1, 1, 5, 3)

    def test_channels_three(self):
        H = HilbertSpace(n=5, k=3, batch_size=2, n_channels=3, device='cpu')
        basis = H.standard_basis()
        assert basis.V.shape == (2, 3, 5, 3)


# ============================================================================
# Tests for string representations
# ============================================================================

class TestStringRepresentations:
    def test_repr(self, basic_hilbert):
        assert "HilbertSpace" in repr(basic_hilbert)

    def test_str(self, basic_hilbert):
        assert "Hilbert space" in str(basic_hilbert)

    def test_basis_repr(self, basic_hilbert):
        basis = basic_hilbert.standard_basis()
        assert "Basis" in repr(basis)


# ============================================================================
# Tests for Monte Carlo
# ============================================================================

class TestMonteCarlo:
    def test_batch_vector_generation(self, basic_hilbert):
        basis = basic_hilbert.standard_basis()
        for _ in range(10):
            v = basis.random_subspace_vector(normalize=True)
            assert v.shape == (2, 2, 5, 1)

    def test_batch_trace_consistency(self, basic_hilbert):
        basis = basic_hilbert.standard_basis()
        torch.manual_seed(42)
        v1 = basis.random_subspace_vector(normalize=True)
        torch.manual_seed(42)
        v2 = basis.random_subspace_vector(normalize=True)
        assert torch.allclose(v1, v2)


# ============================================================================
# Run tests if executed directly
# ============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])