"""Functional calculus at four levels: self-adjoint, normal, singular values, holomorphic (contour)."""
import math
import warnings

import pytest
import torch

import torch_vn_algebra as vn
from torch_vn_algebra import TypeIAlgebra, functional as fc


def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol


@pytest.fixture
def alg():
    return TypeIAlgebra([4, 2, 4], [4, 2, 4], precision='double')


# ----------------------------------------------------------------------
# level 3: singular values without SVD (Hermitian dilation)
# ----------------------------------------------------------------------
class TestDilation:
    @pytest.mark.parametrize('k', [5, 40])
    def test_matches_svd(self, k):
        torch.manual_seed(0)
        X = torch.randn(3, k, k, dtype=torch.complex128)
        X[1, :, -2:] = 0                                   # rank deficient
        s_ref = torch.linalg.svdvals(X)
        assert torch.allclose(fc.singular_values(X, 'eigh'), s_ref, atol=1e-10)
        assert torch.allclose(fc.abs_(X, 'eigh'), fc.abs_(X, 'svd'), atol=1e-9)
        assert torch.allclose(fc.pinv(X, None, 'eigh'), torch.linalg.pinv(X), atol=1e-8)
        W, P = fc.polar(X[[0, 2]], 'eigh')
        assert torch.allclose(W @ P, X[[0, 2]], atol=1e-9)
        assert torch.allclose(W.conj().transpose(-2, -1) @ W, torch.eye(k, dtype=X.dtype).expand(2, k, k), atol=1e-9)

    @pytest.mark.parametrize('method', ['eigh', 'svd'])
    def test_abs_gradient(self, method):
        torch.manual_seed(1)
        X = torch.randn(2, 4, 4, dtype=torch.complex128, requires_grad=True)
        assert torch.autograd.gradcheck(lambda Y: fc.abs_(Y, method), (X,), eps=1e-6, atol=1e-5)

    def test_abs_gradient_finite_at_degenerate_singular_values(self):
        torch.manual_seed(2)
        Q1, _ = torch.linalg.qr(torch.randn(3, 3, dtype=torch.float64))
        Q2, _ = torch.linalg.qr(torch.randn(3, 3, dtype=torch.float64))
        A0 = Q1 @ torch.diag(torch.tensor([1.0, 1.0, 2.0], dtype=torch.float64)) @ Q2.T
        A = A0.clone().requires_grad_(True)
        G = torch.randn(3, 3, dtype=torch.float64)
        (fc.abs_(A[None], 'eigh')[0] * G).sum().backward()
        assert torch.isfinite(A.grad).all()
        H, h = torch.randn(3, 3, dtype=torch.float64), 1e-6
        fd = ((fc.abs_((A0 + h * H)[None], 'svd') - fc.abs_((A0 - h * H)[None], 'svd'))[0] / (2 * h) * G).sum()
        assert abs((A.grad * H).sum() - fd) < 1e-5

    def test_pinv_gradient(self):
        torch.manual_seed(3)
        X = torch.randn(2, 4, 4, dtype=torch.float64, requires_grad=True)
        assert torch.autograd.gradcheck(lambda Y: fc.pinv(Y, None, 'eigh'), (X,), eps=1e-6, atol=1e-5)

    def test_operator_methods(self, alg):
        torch.manual_seed(4)
        A = alg.operator(alg.from_blocks([torch.randn(4, 4, dtype=torch.complex128),
                                           torch.randn(2, 2, dtype=torch.complex128),
                                           torch.randn(4, 4, dtype=torch.complex128)]).matrix)
        W, P = A.polar()
        assert close(W @ P, A, 1e-9) and P.is_positive
        assert close(A @ A.inverse(), alg.identity(), 1e-8)
        s = torch.cat([torch.linalg.svdvals(A.matrix[0, c, :k, :k]) for c, k in enumerate(alg.k_factors)])
        assert abs(A.trace_norm().item() - s.sum().item()) < 1e-10
        assert abs(A.operator_norm().item() - s.max().item()) < 1e-10


# ----------------------------------------------------------------------
# level 2: normal operators
# ----------------------------------------------------------------------
class TestNormal:
    def test_unitary_functions(self, alg):
        torch.manual_seed(5)
        U = alg.random_unitary_operator(3)
        L = U.apply_function(torch.log)                    # principal logarithm on the circle
        assert close(L.expm(), U, 1e-9)
        assert close(L + L.adjoint(), 0 * L.matrix, 1e-9)             # log U is i times a self-adjoint operator
        R = U.apply_function(torch.sqrt)
        assert close(R @ R, U, 1e-9)
        assert close(U.apply_function(lambda z: z ** 3), U @ U @ U, 1e-9)

    def test_real_orthogonal(self):
        torch.manual_seed(6)
        alg = TypeIAlgebra([3], [3], complex_valued=False, precision='double')
        O = alg.random_unitary_operator(2)                  # orthogonal
        assert close(O.apply_function(lambda z: z ** 2), O @ O, 1e-9)
        with pytest.raises(ValueError):
            O.apply_function(torch.log).matrix              # complex logarithm in a real algebra

    def test_non_normal_rejected(self, alg):
        A = alg.operator(alg.from_blocks([torch.triu(torch.ones(4, 4, dtype=torch.complex128)),
                                           torch.eye(2, dtype=torch.complex128),
                                           torch.eye(4, dtype=torch.complex128)]).matrix)
        with pytest.raises(RuntimeError):
            A.apply_function(torch.exp).matrix

    def test_gradient_equals_holomorphic_derivative(self):
        """Derivative of f at a unitary in a non-normal direction, against the contour calculus."""
        torch.manual_seed(7)
        alg = TypeIAlgebra([3], [3], precision='double')
        U0 = alg.random_unitary_operator(1).matrix
        f = torch.exp
        X = U0.clone().requires_grad_(True)
        F = alg.operator(X, is_normal=True).apply_function(f).matrix
        G = torch.randn_like(F)
        (F * G.conj()).real.sum().backward()
        H, h = torch.randn_like(U0), 1e-5
        Fp = alg.operator(U0 + h * H).holomorphic_function(f).matrix
        Fm = alg.operator(U0 - h * H).holomorphic_function(f).matrix
        fd = ((Fp - Fm) / (2 * h) * G.conj()).real.sum()
        assert abs((X.grad.conj() * H).real.sum() - fd) < 1e-6


# ----------------------------------------------------------------------
# level 4: holomorphic calculus by contour integration
# ----------------------------------------------------------------------
class TestHolomorphic:
    def test_exp_of_non_normal_and_defective(self, alg):
        torch.manual_seed(8)
        J = torch.diag(torch.ones(3, dtype=torch.complex128), 1) + 0.5 * torch.eye(4, dtype=torch.complex128)
        A = alg.operator(alg.from_blocks([J, torch.randn(2, 2, dtype=torch.complex128),
                                           torch.randn(4, 4, dtype=torch.complex128)]).matrix)
        E = A.holomorphic_function(torch.exp)
        assert close(E, A.expm(), 1e-9)

    def test_log_with_explicit_contour(self):
        torch.manual_seed(9)
        alg = TypeIAlgebra([4], [4], precision='double')
        N = torch.triu(torch.randn(4, 4, dtype=torch.complex128), 1)
        A = alg.operator((3 * torch.eye(4, dtype=torch.complex128) + 0.5 * N)[None, None])
        L = A.holomorphic_function(torch.log, center=3.0, radius=2.0)     # the disk avoids the branch point 0
        assert close(L.expm(), A, 1e-9)

    def test_real_algebra(self):
        torch.manual_seed(10)
        alg = TypeIAlgebra([3], [3], complex_valued=False, precision='double')
        A = alg.operator(torch.randn(1, 1, 3, 3, dtype=torch.float64))
        E = A.holomorphic_function(torch.exp)
        assert not torch.is_complex(E.matrix) and close(E, A.expm(), 1e-9)

    def test_gradient(self):
        torch.manual_seed(11)
        alg = TypeIAlgebra([3], [3], precision='double')
        X = torch.randn(1, 1, 3, 3, dtype=torch.complex128, requires_grad=True)
        g = lambda Y: alg.operator(Y).holomorphic_function(torch.exp).matrix
        assert torch.autograd.gradcheck(g, (X,), eps=1e-6, atol=1e-5)

    def test_diagnostics(self):
        alg = TypeIAlgebra([2], [2], precision='double')
        A = alg.operator(torch.diag(torch.tensor([1.0, 1.5], dtype=torch.complex128))[None, None])
        with pytest.raises(ValueError):                    # the contour passes through the pole z = 3
            A.holomorphic_function(lambda z: 1 / (z - 3), center=1.5, radius=1.5).matrix
        with pytest.raises(ValueError):                    # too few nodes for a large circle: not converged
            A.holomorphic_function(torch.exp, radius=50.0, nodes=4, max_nodes=8).matrix
        B = torch.randn(6, 6, dtype=torch.float64)
        with pytest.raises(ValueError):                    # catastrophic cancellation for exp of a large spread
            TypeIAlgebra([6], [6], precision='double').operator(((B + B.T) * 50)[None, None]).holomorphic_function(
                torch.exp).matrix


# ----------------------------------------------------------------------
# spectral cache
# ----------------------------------------------------------------------
class TestCache:
    def test_one_decomposition_for_many_functions(self, alg, monkeypatch):
        A = alg.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(2, d), batch_size=2, force_positive=True)
        A.matrix
        calls = []
        real = torch.linalg.eigh
        monkeypatch.setattr(torch.linalg, 'eigh', lambda x, *a, **k: (calls.append(1), real(x, *a, **k))[1])
        with vn.caching():
            R, L, P = A.sqrt().matrix, A.log().matrix, A.power(0.3).matrix
            ev = A.eigenvalues()
        assert len(calls) == 2                             # two block sizes, decomposed once
        monkeypatch.setattr(torch.linalg, 'eigh', real)
        assert close(R, A.apply_function(torch.sqrt)) and close(L, A.apply_function(torch.log))
        assert A.cache_info()['cached']

    def test_invalidated_by_in_place_change(self, alg):
        A = alg.operator_from_eigenvalues(lambda d: 1 + torch.rand(1, d), force_positive=True).cache_spectrum()
        before = A.sqrt().matrix
        A.matrix.mul_(4.0)                                 # in-place change of the stored matrix
        assert close(A.sqrt().matrix, 2 * before)

    def test_propagators_reuse_the_decomposition(self, alg, monkeypatch):
        from torch_vn_algebra import dynamics
        H = alg.operator_from_eigenvalues(lambda d: torch.randn(1, d), force_self_adjoint=True).cache_spectrum()
        U1 = dynamics.propagator(H, 0.3).matrix
        calls = []
        real = torch.linalg.eigh
        monkeypatch.setattr(torch.linalg, 'eigh', lambda x, *a, **k: (calls.append(1), real(x, *a, **k))[1])
        Us = [dynamics.propagator(H, t).matrix for t in (0.1, 0.2, 0.3)]
        assert not calls and close(Us[2], U1)


class TestReviewRegressions:
    """Defects found in the review of the functional calculus."""

    @pytest.mark.parametrize('prec,k', [('double', 400), ('single', 400)])
    def test_large_unitaries_are_normal(self, prec, k):
        torch.manual_seed(12)
        alg = TypeIAlgebra([k], [k], precision=prec)
        U = alg.random_unitary_operator(2)
        tol = 1e-12 if prec == 'double' else 1e-5
        assert close(U.apply_function(lambda z: z ** 2), U @ U, tol)

    def test_second_derivatives_normal(self):
        torch.manual_seed(13)
        Q, R = torch.linalg.qr(torch.randn(4, 4, dtype=torch.complex128))
        U = (Q * (R.diagonal() / R.diagonal().abs()))[None]
        G, H = torch.randn_like(U), torch.randn_like(U)

        def hvp(F):
            X = U.clone().requires_grad_(True)
            g, = torch.autograd.grad((F(X) * G.conj()).real.sum(), X, create_graph=True)
            return torch.autograd.grad((g.conj() * H).real.sum(), X)[0]
        assert (hvp(lambda Y: fc.spectral_apply(Y, torch.exp, 'normal')) - hvp(torch.linalg.matrix_exp)).abs().max() < 1e-10

    @pytest.mark.parametrize('method', ['svd', 'eigh'])
    def test_abs_gradient_rank_deficient(self, method):
        torch.manual_seed(14)
        k = 6
        Q1, _ = torch.linalg.qr(torch.randn(k, k, dtype=torch.complex128))
        Q2, _ = torch.linalg.qr(torch.randn(k, k, dtype=torch.complex128))
        s = torch.zeros(k, dtype=torch.complex128); s[:3] = torch.tensor([3.0, 2.0, 1.0])
        X0 = (Q1 @ torch.diag(s) @ Q2.mH)[None]
        M = torch.zeros(k, k, dtype=torch.complex128); M[:3, :3] = torch.randn(3, 3, dtype=torch.complex128)
        E = (Q1 @ M @ Q2.mH)[None]                          # a direction that keeps the kernel
        G = torch.randn_like(X0)
        X = X0.clone().requires_grad_(True)
        (fc.abs_(X, method) * G.conj()).real.sum().backward()
        h = 1e-6
        fd = ((fc.abs_(X0 + h * E, 'svd') - fc.abs_(X0 - h * E, 'svd')) / (2 * h) * G.conj()).real.sum()
        assert abs((X.grad.conj() * E).real.sum() - fd) < 1e-6

    def test_polar_factor_unitary_for_ill_conditioned(self):
        torch.manual_seed(15)
        k = 40
        Q1, _ = torch.linalg.qr(torch.randn(k, k, dtype=torch.complex128))
        Q2, _ = torch.linalg.qr(torch.randn(k, k, dtype=torch.complex128))
        s = torch.ones(k, dtype=torch.complex128); s[-1] = 1e-13
        X = (Q1 @ torch.diag(s) @ Q2.mH)[None]
        W, P = fc.polar(X, 'eigh')
        assert (W.mH @ W - torch.eye(k, dtype=X.dtype)).abs().max() < 1e-12
        assert (W @ P - X).abs().max() < 1e-12

    def test_contour_exact_for_nilpotent_and_tight_default_circle(self):
        N = torch.diag(torch.ones(5, dtype=torch.complex128), 1)
        F, err = fc.holomorphic(N[None], torch.exp)
        assert (F - torch.linalg.matrix_exp(N[None])).abs().max() < 1e-12
        B = torch.randn(6, 6, dtype=torch.float64)
        F, err = fc.holomorphic(((B + B.T) * 2)[None], torch.exp)      # the old radius lost 7 digits here
        ref = torch.linalg.matrix_exp(((B + B.T) * 2)[None])
        assert ((F - ref).abs().max() / ref.abs().max()) < 1e-11
