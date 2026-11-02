"""Numerical guarantees of the paper, checked against analytic and dense references."""
import pytest
import torch

from torch_vn_algebra import TypeIAlgebra, functional, learn


@pytest.mark.parametrize('complex_valued', [False, True])
@pytest.mark.parametrize('precision', ['single', 'double'])
def test_elementwise_mixed_derivative_matches_diagonal_formula(complex_valued, precision):
    alg = TypeIAlgebra([2], [2], complex_valued=complex_valued, precision=precision)
    real_dtype = torch.float32 if precision == 'single' else torch.float64
    values = torch.tensor([.8, 1.7], dtype=real_dtype)
    X = torch.diag(values).to(alg.dtype)[None, None].requires_grad_()
    t = torch.tensor(.6, dtype=real_dtype, requires_grad=True)
    weights = torch.tensor([.7, -.3], dtype=real_dtype)
    Y = alg.operator(X, is_self_adjoint=True).apply_function(lambda w: torch.exp(t * w)).matrix
    loss = (Y.diagonal(dim1=-2, dim2=-1).real * weights).sum()
    dt, = torch.autograd.grad(loss, t, create_graph=True)
    mixed, = torch.autograd.grad(dt, X)
    expected = weights * (1 + t.detach() * values) * torch.exp(t.detach() * values)
    torch.testing.assert_close(mixed.real.diagonal(dim1=-2, dim2=-1).flatten(), expected)


@pytest.mark.parametrize('kind', ['hermitian', 'normal'])
@pytest.mark.parametrize('per_block', [False, True])
@pytest.mark.parametrize('cached', [False, True])
def test_parameter_and_matrix_hessian_matches_matrix_exp(kind, per_block, cached):
    alg = TypeIAlgebra([2], [2], precision='double')
    Q = torch.tensor([[1., 1.j], [1.j, 1.]], dtype=alg.dtype) / 2 ** .5
    spectrum = torch.tensor([.8, 1.7], dtype=alg.dtype)
    if kind == 'normal':
        spectrum = spectrum + torch.tensor([.4j, -.2j], dtype=alg.dtype)
    X = (Q @ torch.diag(spectrum) @ Q.mH)[None, None].requires_grad_()
    t = torch.tensor(.6, dtype=torch.float64, requires_grad=True)
    op = alg.operator(X, is_self_adjoint=kind == 'hermitian', is_normal=True)
    if cached:
        op.cache_spectrum()
    actual = op.apply_function(lambda w: torch.exp(t * w), per_block=per_block).matrix
    H = .5 * (X + X.mH) if kind == 'hermitian' else X
    expected = torch.matrix_exp(t * H)
    G = torch.tensor([[.7, .2 + .1j], [-.3j, -.3]], dtype=alg.dtype)[None, None]
    direction = torch.tensor([[.3, -.2j], [.1j, .4]], dtype=alg.dtype)[None, None]

    def derivatives(Y):
        first = torch.autograd.grad((Y * G.conj()).real.sum(), (X, t),
                                    create_graph=True, retain_graph=True)
        # Both orders of mixed differentiation, as well as the pure second derivatives.
        probes = ((first[0] * direction.conj()).real.sum(), first[1])
        second = [torch.autograd.grad(probe, (X, t), retain_graph=True) for probe in probes]
        return first + tuple(v for pair in second for v in pair)

    for a, e in zip(derivatives(actual), derivatives(expected)):
        torch.testing.assert_close(a, e, atol=1e-9, rtol=1e-9)


@pytest.mark.parametrize('method', ['rk4', 'unsupported'])
def test_lindblad_layer_rejects_methods_without_a_channel_guarantee(method):
    alg = TypeIAlgebra([2], [2], precision='double')
    with pytest.raises(ValueError, match='exact.*kraus'):
        learn.LindbladLayer(alg, method=method)


@pytest.mark.parametrize('complex_valued', [False, True])
def test_elementwise_matrix_and_parameter_gradgradcheck(complex_valued):
    alg = TypeIAlgebra([2], [2], complex_valued=complex_valued, precision='double')
    X = torch.tensor([[[[.8, .15], [.15, 1.7]]]], dtype=alg.dtype).requires_grad_()
    t = torch.tensor(.6, dtype=torch.float64, requires_grad=True)

    def evaluate(matrix, scale):
        return alg.operator(matrix, is_self_adjoint=True).apply_function(lambda w: (scale * w).exp()).matrix

    assert torch.autograd.gradgradcheck(evaluate, (X, t), atol=1e-5)


@pytest.mark.parametrize('coupled', [False, True])
def test_normal_mixed_derivative_with_permuted_cached_eigenvalues(coupled):
    dtype = torch.complex128
    Q = torch.tensor([[1., 1.j], [1.j, 1.]], dtype=dtype) / 2 ** .5
    w = torch.tensor([1.7 - .2j, .8 + .4j], dtype=dtype)[None]
    X = (Q @ torch.diag_embed(w) @ Q.mH).requires_grad_()
    t = torch.tensor(.6, dtype=torch.float64, requires_grad=True)
    actual = functional.spectral_apply(X, lambda z: (t * z).exp(), 'normal',
                                       decomposition=(w, Q[None]), coupled=coupled)
    expected = torch.matrix_exp(t * X)
    G = torch.tensor([[[.7, .2 + .1j], [-.3j, -.3]]], dtype=dtype)
    for Y in (actual, expected):
        dt, = torch.autograd.grad((Y * G.conj()).real.sum(), t, create_graph=True, retain_graph=True)
        mixed, = torch.autograd.grad(dt, X, retain_graph=True)
        if Y is actual:
            result = mixed
        else:
            torch.testing.assert_close(result, mixed, atol=1e-9, rtol=1e-9)


@pytest.mark.parametrize('method', ['exact', 'kraus'])
@pytest.mark.parametrize('precision', ['single', 'double'])
def test_lindblad_layer_preserves_positivity_at_a_large_step(method, precision):
    alg = TypeIAlgebra([2], [2], precision=precision)
    layer = learn.LindbladLayer(alg, hamiltonian=False, n_jumps=1, method=method, time=1., substeps=1)
    with torch.no_grad():
        layer.L.zero_()
        layer.L[..., 0, 1] = 2.
    rho = alg.operator(torch.tensor([[[[0., 0.], [0., 1.]]]], dtype=alg.dtype),
                       is_self_adjoint=True, is_positive=True)
    output = layer(rho)
    tol = 1e-5 if precision == 'single' else 1e-12
    assert torch.linalg.eigvalsh(output.matrix).min() >= -tol
    torch.testing.assert_close(output.trace.real, torch.ones_like(output.trace.real), atol=tol, rtol=tol)
    if method == 'exact':
        excited = torch.exp(torch.tensor(-4., dtype=output.matrix.real.dtype))
        expected = torch.diag(torch.stack([1 - excited, excited])).to(alg.dtype)[None, None]
        torch.testing.assert_close(output.matrix, expected, atol=tol, rtol=tol)


@pytest.mark.parametrize('coupled', [False, True])
def test_degenerate_normal_first_derivative_with_create_graph_remains_finite(coupled):
    X = (1.j * torch.eye(2, dtype=torch.complex128))[None].requires_grad_()
    G = torch.tensor([[[.7, .2 + .1j], [-.3j, -.3]]], dtype=X.dtype)
    Y = functional.spectral_apply(X, torch.exp, 'normal', coupled=coupled)
    with pytest.warns(RuntimeWarning, match='degenerate'):
        actual, = torch.autograd.grad((Y * G.conj()).real.sum(), X, create_graph=True)
    expected, = torch.autograd.grad((torch.matrix_exp(X) * G.conj()).real.sum(), X)
    torch.testing.assert_close(actual, expected)
