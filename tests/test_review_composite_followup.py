"""Regression tests for composite property inference, lazy mutation and precision boundaries."""
import pytest
import torch

from torch_vn_algebra import DensityMatrix, TypeIAlgebra, composite as cp


PRODUCTS = [cp.tensor_product, cp.fused_tensor_product]


def algebra(sizes, precision='double', complex_valued=True, device='cpu'):
    return TypeIAlgebra(sizes, sizes, precision=precision, complex_valued=complex_valued, device=device)


def initial(alg, state):
    if state:
        return DensityMatrix(alg, matrix=alg.identity().matrix / alg.total_subspace_dim)
    return alg.identity()


def replacement(alg, kind):
    # Both matrices have trace one, so mutated state inputs keep their normalisation.
    values = [[-.5, 0.], [0., 1.5]] if kind == 'negative' else [[1., 1.], [0., 0.]]
    return torch.tensor(values, dtype=alg.dtype)


def dense_product(op):
    """Embed stored sectors into the ordinary tensor basis for an independent dense comparison."""
    a, b = op.algebra.tensor_factors
    n = a.total_subspace_dim * b.total_subspace_dim
    out = torch.zeros(n, n, dtype=op.matrix.dtype)
    for s, pairs in enumerate(op.algebra.product_blocks):
        rows = torch.tensor([(sum(a.k_factors[:c]) + i) * b.total_subspace_dim + sum(b.k_factors[:d]) + j
                             for c, d, _ in pairs for i in range(a.k_factors[c]) for j in range(b.k_factors[d])])
        out[rows[:, None], rows[None, :]] = op.matrix[0, s, :len(rows), :len(rows)]
    return out


@pytest.mark.parametrize('product', PRODUCTS)
@pytest.mark.parametrize('keep', [1, 2])
@pytest.mark.parametrize('kind', ['negative', 'nonhermitian'])
@pytest.mark.parametrize('complex_valued', [False, True])
def test_partial_trace_can_restore_properties_cached_false(product, keep, kind, complex_valued):
    one, two = algebra([1], complex_valued=complex_valued), algebra([2], complex_valued=complex_valued)
    ab = product(one, two) if keep == 1 else product(two, one)
    values = [[1., 0.], [0., -1.]] if kind == 'negative' else [[0., 1.], [0., 0.]]
    x = ab.operator(torch.tensor(values, dtype=ab.dtype)[None, None])
    if kind == 'negative':
        assert not x.is_positive
    else:
        assert not x.is_self_adjoint
    reduced = cp.partial_trace(x, keep)
    # The scalar reduction is the ordinary trace, exactly zero for either input.
    zero = torch.zeros(1, 1, 1, 1, dtype=ab.dtype)
    torch.testing.assert_close(reduced.matrix, zero)
    assert reduced.is_self_adjoint and reduced.is_positive
    torch.testing.assert_close(reduced.sqrt().matrix, zero)


@pytest.mark.parametrize('product', PRODUCTS)
@pytest.mark.parametrize('keep', [1, 2])
@pytest.mark.parametrize('kind', ['negative', 'nonhermitian'])
@pytest.mark.parametrize('state', [False, True])
@pytest.mark.parametrize('mutation', ['before_construction', 'pending'])
def test_partial_trace_refreshes_parent_properties_and_keeps_snapshots(product, keep, kind, state, mutation):
    a, b = algebra([2]), algebra([2])
    ab = product(a, b)
    x = initial(ab, state)
    snapshot = cp.partial_trace(x, keep)
    snapshot_matrix = snapshot.matrix.clone()
    target = replacement(a, kind)
    other = torch.eye(2, dtype=ab.dtype) / (2 if state else 1)
    updated = torch.kron(target, other) if keep == 1 else torch.kron(other, target)
    if mutation == 'before_construction':
        x.matrix.copy_(updated[None, None])
    pending = cp.partial_trace(x, keep)
    descendant = pending + pending
    if mutation == 'pending':
        x.matrix.copy_(updated[None, None])
    expected = target * (1 if state else 2)
    assert isinstance(pending, DensityMatrix) == state
    # Query flags before materialisation too: descendants must see parent invalidation transitively.
    assert not pending.is_positive
    assert descendant.is_self_adjoint == (kind == 'negative')
    torch.testing.assert_close(pending.matrix, expected[None, None])
    torch.testing.assert_close(descendant.matrix, (2 * expected)[None, None])
    torch.testing.assert_close(pending.operator_norm(), torch.linalg.matrix_norm(expected, ord=2)[None])
    with pytest.raises(RuntimeError, match='positive'):
        pending.sqrt()
    torch.testing.assert_close(snapshot.matrix, snapshot_matrix)
    assert snapshot.is_positive and snapshot.is_self_adjoint


@pytest.mark.parametrize('product', PRODUCTS)
@pytest.mark.parametrize('factor', [1, 2])
@pytest.mark.parametrize('kind', ['negative', 'nonhermitian'])
@pytest.mark.parametrize('state', [False, True])
@pytest.mark.parametrize('mutation', ['before_construction', 'pending'])
def test_kron_refreshes_both_parents_and_keeps_snapshots(product, factor, kind, state, mutation):
    a, b = algebra([2]), algebra([2])
    ab = product(a, b)
    left, right = initial(a, state), initial(b, state)
    x = left if factor == 1 else right
    snapshot = cp.kron(left, right, ab)
    snapshot_matrix = snapshot.matrix.clone()
    target = replacement(a, kind)
    if mutation == 'before_construction':
        x.matrix.copy_(target[None, None])
    pending = cp.kron(left, right, ab)
    descendant = pending + pending
    if mutation == 'pending':
        x.matrix.copy_(target[None, None])
    other = torch.eye(2, dtype=ab.dtype) / (2 if state else 1)
    expected = torch.kron(target, other) if factor == 1 else torch.kron(other, target)
    assert isinstance(pending, DensityMatrix) == state
    assert not pending.is_positive
    assert descendant.is_self_adjoint == (kind == 'negative')
    torch.testing.assert_close(pending.matrix, expected[None, None])
    torch.testing.assert_close(descendant.matrix, (2 * expected)[None, None])
    torch.testing.assert_close(pending.operator_norm(), torch.linalg.matrix_norm(expected, ord=2)[None])
    with pytest.raises(RuntimeError, match='positive'):
        pending.sqrt()
    torch.testing.assert_close(snapshot.matrix, snapshot_matrix)
    assert snapshot.is_positive and snapshot.is_self_adjoint


@pytest.mark.parametrize('product', PRODUCTS)
@pytest.mark.parametrize('complex_valued', [False, True])
@pytest.mark.parametrize('precision', ['single', 'double'])
def test_precision_mismatch_requires_explicit_conversion(product, complex_valued, precision):
    a = algebra([1, 2], precision, complex_valued)
    b = algebra([2, 1], 'double' if precision == 'single' else 'single', complex_valued)
    with pytest.raises(ValueError, match='same precision'):
        product(a, b)
    compatible = b.like(precision=precision, device=torch.device('cpu'))
    ab = product(a, compatible)
    x = cp.kron(a.identity(), b.identity().cast(compatible), ab)
    expected = torch.kron(torch.eye(3, dtype=a.dtype), torch.eye(3, dtype=a.dtype))
    assert ab.dtype == x.matrix.dtype == a.dtype
    torch.testing.assert_close(dense_product(x), expected)
    torch.testing.assert_close((x @ ab.identity()).matrix, x.matrix)
    assert x.is_positive and x.is_self_adjoint
    full = expected.reshape(3, 3, 3, 3)
    for keep in (1, 2):
        reduced = cp.partial_trace(x, keep)
        dense = torch.einsum('ijlj->il', full) if keep == 1 else torch.einsum('ijil->jl', full)
        assert reduced.algebra.dtype == reduced.matrix.dtype == a.dtype
        actual = torch.block_diag(*[reduced.matrix[0, c, :k, :k] for c, k in enumerate(reduced.algebra.k_factors)])
        torch.testing.assert_close(actual, dense)
        torch.testing.assert_close(cp.partial_trace_channel(ab, keep)(x).matrix, reduced.matrix)


@pytest.mark.parametrize('product', PRODUCTS)
def test_product_requires_both_factor_devices_to_match(product):
    with pytest.raises(ValueError, match='same device'):
        product(algebra([2], device='cpu'), algebra([2], device='meta'))
