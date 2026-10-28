"""
The two storage layouts of TypeIAlgebra (padded (batch, C, k_max, k_max) and packed: one tensor (batch, n_g, k_g, k_g)
per block size) give the same results; the packed layout never forms the padded tensor unless .matrix is requested.
The whole test suite also runs with TORCH_VN_ALGEBRA_LAYOUT=packed (CI job), which makes every algebra packed.
"""
import pytest
import torch

from torch_vn_algebra import TypeIAlgebra, states
from torch_vn_algebra.algebra import _auto_layout

SIZES = [[3, 1, 4, 1, 3], [2, 0, 5, 2], [4, 4], [1, 6, 15, 20, 15, 6, 1]]   # non-contiguous groups, an empty sector


def pair(k, **kw):
    kw.setdefault('precision', 'double')
    n = [max(x, 1) for x in k]
    return (TypeIAlgebra(n, k, layout='padded', **kw), TypeIAlgebra(n, k, layout='packed', **kw))


def same(x, y, tol=1e-10):
    x = x.matrix if hasattr(x, 'matrix') else x
    y = y.matrix if hasattr(y, 'matrix') else y
    return x.shape == y.shape and torch.allclose(x, y, atol=tol, rtol=tol)


def twins(k, kind='herm', B=3, seed=0):
    """The same random operator in both layouts."""
    pad, pak = pair(k)
    torch.manual_seed(seed)
    X = torch.randn(B, pad.C, pad.k_max, pad.k_max, dtype=pad.dtype)
    if kind == 'herm':
        X = X + X.mH
    elif kind == 'pos':
        X = X @ X.mH + 0.1 * torch.eye(pad.k_max, dtype=pad.dtype)
    flags = {'is_self_adjoint': True} if kind in ('herm', 'pos') else {}
    return pad.operator(X, **flags), pak.operator(X, **flags)


@pytest.mark.parametrize('k', SIZES, ids=str)
def test_storage_and_round_trip(k):
    a, b = twins(k, 'any')
    pak = b.algebra
    assert (b._matrix is None or pak._uniform) and len(b.parts) == len(pak.size_groups)
    assert b._bytes_held == pak.bytes_per_operator(3)
    for (kg, idx), p in zip(pak.size_groups, b.parts):
        assert p.shape == (3, len(idx), kg, kg)
        for j, c in enumerate(idx):
            assert torch.equal(p[:, j], a.matrix[:, c, :kg, :kg]) and torch.equal(b.block(c), a.block(c))
    assert same(a, b) and b.shape == a.shape
    assert pak.bytes_per_operator(3) == 3 * (sum(x * x for x in k) if not pak._uniform else pak.C * pak.k_max ** 2) * 16
    c = pak.from_parts(b.parts)
    assert same(c, a)
    d = pak.from_blocks([a.block(i) for i in range(pak.C)])
    assert same(d, a)


@pytest.mark.parametrize('k', SIZES, ids=str)
def test_algebra_operations_agree(k):
    a, b = twins(k, 'any')
    h, g = twins(k, 'herm', seed=1)
    p, q = twins(k, 'pos', seed=2)
    pairs = [(a + h, b + g), (a - h, b - g), (a @ h, b @ g), (-a, -b), (a.adjoint(), b.adjoint()),
             (a * 2.5, b * 2.5), (a * torch.arange(1.0, 4.0, dtype=torch.float64),
                                  b * torch.arange(1.0, 4.0, dtype=torch.float64)),
             (h.apply_function(torch.exp), g.apply_function(torch.exp)), (p.sqrt(), q.sqrt()), (p.log(), q.log()),
             (p.power(-0.5), q.power(-0.5)), (a.expm(0.3), b.expm(0.3)), (a.inverse(), b.inverse()),
             (a.abs(), b.abs()), (a.polar()[0], b.polar()[0]), (a.polar()[1], b.polar()[1]),
             (h.inverse(), g.inverse()), (a.holomorphic_function(torch.exp), b.holomorphic_function(torch.exp))]
    for x, y in pairs:
        assert same(x, y, 1e-8)
    for f in ('Tr_blunt', 'Tr_norm', 'tau_vN', 'trace_norm', 'frobenius_norm', 'operator_norm', 'sector_traces'):
        assert same(getattr(a, f)(), getattr(b, f)()), f
    assert same(p.entropy(), q.entropy()) and same(p.michelson_contrast, q.michelson_contrast)
    assert same(h.lambda_max, g.lambda_max) and same(h.lambda_min, g.lambda_min)
    for x, y in zip(h.eigenvalues(), g.eigenvalues()):
        assert same(x, y)
    for (w1, v1), (w2, v2) in zip(h.eigh(), g.eigh()):
        assert same(w1, w2) and same(v1 @ torch.diag_embed(w1).to(v1.dtype) @ v1.mH,
                                     v2 @ torch.diag_embed(w2).to(v2.dtype) @ v2.mH)
    assert (h.is_self_adjoint, h.is_normal, h.is_positive, a.is_invertible, a.is_projection) == \
           (g.is_self_adjoint, g.is_normal, g.is_positive, b.is_invertible, b.is_projection)


@pytest.mark.parametrize('k', SIZES, ids=str)
def test_packed_never_forms_the_padded_tensor(k):
    _, b = twins(k, 'pos')
    pak = b.algebra
    ops = [b @ b + b * 2.0 - b, b.sqrt(), b.log(), b.apply_function(torch.exp), b.inverse(), b.abs(),
           pak.identity(3), pak.zero(3), pak.central(list(range(pak.C))), pak.random_unitary_operator(3),
           pak.operator_from_eigenvalues(lambda d: torch.rand(3, d, dtype=torch.float64), 3, force_positive=True),
           states.random_density_matrix(pak, 3), states.gibbs_state(b, 0.7), states.tracial_state(pak, 3),
           states.maximally_mixed_state(pak, 3)]
    for op in ops:
        op.Tr_blunt(), op.Tr_norm(), op.frobenius_norm(), op.is_self_adjoint
        assert op._matrix is None or pak._uniform
    rho = states.random_density_matrix(pak, 3)
    for v in (states.von_neumann_entropy(rho), states.purity(rho), states.trace_of_product(rho, b),
              states.sector_probabilities(rho), states.fidelity(rho, rho), states.relative_entropy(rho, rho)):
        assert v.shape[0] == 3
    assert rho._matrix is None or pak._uniform and b._matrix is None or pak._uniform


@pytest.mark.parametrize('k', SIZES, ids=str)
def test_factories_and_states_agree_with_the_same_seed(k):
    pad, pak = pair(k)
    for make in (lambda A: A.identity(2), lambda A: A.zero(2), lambda A: A.central([0.5 * c for c in range(A.C)]),
                 lambda A: states.tracial_state(A, 2), lambda A: states.maximally_mixed_state(A, 2)):
        assert same(make(pad), make(pak))
    torch.manual_seed(5)
    U1 = pad.random_unitary_operator(2).matrix                   # evaluated now (lazy otherwise)
    torch.manual_seed(5)
    U2 = pak.random_unitary_operator(2).matrix                   # evaluated now (lazy otherwise)
    assert same(U1, U2)
    torch.manual_seed(6)
    E1 = pad.operator_from_eigenvalues(lambda d: torch.rand(2, d, dtype=torch.float64), 2, force_positive=True).matrix
    torch.manual_seed(6)
    E2 = pak.operator_from_eigenvalues(lambda d: torch.rand(2, d, dtype=torch.float64), 2, force_positive=True).matrix
    assert same(E1, E2)
    rho = states.random_density_matrix(pak, 4)
    assert same(rho.Tr_blunt(), torch.ones(4, dtype=pak.dtype)) and rho.is_positive
    assert rho.sector_probabilities().shape == (4, pak.C)
    s1, s2 = states.lueders_update(states.tracial_state(pad, 2), pad.identity(2) * 0.0 + pad.identity(2))[0], \
        states.lueders_update(states.tracial_state(pak, 2), pak.identity(2) * 0.0 + pak.identity(2))[0]
    assert same(s1, s2)


def test_gradients_agree():
    k = [3, 1, 4, 1, 3]
    a, b = twins(k, 'herm')
    xa = a.matrix.clone().requires_grad_(True)
    xb = a.matrix.clone().requires_grad_(True)
    la = (a.algebra.operator(xa, is_self_adjoint=True).apply_function(torch.exp).Tr_norm().real.sum())
    lb = (b.algebra.operator(xb, is_self_adjoint=True).apply_function(torch.exp).Tr_norm().real.sum())
    la.backward(), lb.backward()
    assert torch.allclose(xa.grad, xb.grad, atol=1e-10)


def test_cast_between_layouts_and_fields():
    pad, pak = pair([2, 5, 2], complex_valued=False)
    X = pad.operator(torch.randn(2, 3, 5, 5, dtype=torch.float64))
    Y = X.cast(pak.like(complex_valued=True))
    assert Y.algebra.layout == 'packed' and Y._matrix is None and same(Y.matrix.real, X.matrix)
    assert same(Y.cast(pad), X)


def test_in_place_change_of_parts_refreshes_the_padded_view():
    _, b = twins([3, 1, 4], 'any')
    m = b.matrix.clone()
    b.parts[0].mul_(2.0)
    assert not torch.equal(b.matrix, m) and torch.equal(b.matrix[:, 0, :3, :3], 2 * m[:, 0, :3, :3])


def test_layout_choice(monkeypatch):
    from torch_vn_algebra import cost
    monkeypatch.delenv('TORCH_VN_ALGEBRA_LAYOUT', raising=False)      # the CI jobs force a layout
    assert _auto_layout([1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1]) == 'packed'
    assert _auto_layout([4, 4, 4]) == 'padded' and _auto_layout([1, 2, 3]) == 'padded'
    # the cost model: many small block sizes stay padded for small batches and are packed for large ones; more
    # expensive kernel calls (CUDA) need larger batches; much padding is packed whatever the batch
    sizes = list(range(1, 21))
    assert _auto_layout(sizes, typical_batch=1) == 'padded' and _auto_layout(sizes, typical_batch=1000) == 'packed'
    assert _auto_layout(sizes, 'cpu', 100) == 'packed' and _auto_layout(sizes, 'cuda', 100) == 'padded'
    assert _auto_layout([30, 29], typical_batch=16) == 'padded'                     # little padding
    assert _auto_layout([2000, 1, 1], typical_batch=1) == 'packed'                 # 8e6 padding entries per sample
    assert TypeIAlgebra(sizes, sizes, typical_batch=1000).layout == 'packed'
    old = cost.get_limits()
    try:
        cost.set_limits(layout_typical_batch=1)
        assert TypeIAlgebra(sizes, sizes).layout == 'padded'
        cost.set_limits(layout_entries_per_kernel_cpu=1.0)
        assert TypeIAlgebra(sizes, sizes).layout == 'packed'
    finally:
        cost.set_limits(**old)
    assert TypeIAlgebra([2, 3], [2, 3], layout='packed').layout == 'packed'
    monkeypatch.setenv('TORCH_VN_ALGEBRA_LAYOUT', 'packed')
    assert TypeIAlgebra([2, 3], [2, 3]).layout == 'packed'
    assert TypeIAlgebra([2, 3], [2, 3], layout='padded').layout == 'padded'      # an explicit argument wins
    assert TypeIAlgebra([2, 3], [2, 3]).like(precision='double').layout == 'packed'
    monkeypatch.setenv('TORCH_VN_ALGEBRA_LAYOUT', 'bogus')
    with pytest.raises(ValueError):
        TypeIAlgebra([2, 3], [2, 3])
    with pytest.raises(ValueError):
        TypeIAlgebra([2, 3], [2, 3], layout='dense')


def test_modules_agree_across_layouts():
    """Dynamics (all Lindblad methods, with jumps inside and between sectors), channels, tensor products, partial
    traces and spin-chain functions give the same results in both layouts."""
    from torch_vn_algebra import channels, composite, dynamics, SpinChain
    out = {}
    for layout in ('padded', 'packed'):
        pad = TypeIAlgebra([1, 4, 6, 4, 1], [1, 4, 6, 4, 1], precision='double', layout=layout)
        torch.manual_seed(0)
        H = pad.operator(torch.randn(2, 5, 6, 6, dtype=pad.dtype), is_self_adjoint=False)
        H = (H + H.adjoint()) * 0.5
        H._is_self_adjoint = True
        L = pad.operator(torch.randn(1, 5, 6, 6, dtype=pad.dtype) / 3)
        rho = states.DensityMatrix(pad, matrix=states.random_density_matrix(
            TypeIAlgebra([1, 4, 6, 4, 1], [1, 4, 6, 4, 1], precision='double', layout='padded'), 2).matrix)
        jump = channels.random_inter_sector_channel(pad, pad, kraus_rank=2)
        r = {}
        for m in ('rk4', 'kraus', 'exact'):
            r[m] = dynamics.lindblad_evolve(rho, H, [L, jump], [0.0, 0.3, 0.6], rates=[0.5, 0.4], method=m)[-1].matrix
        r['vn'] = dynamics.von_neumann(rho, H, [0.0, 0.7])[-1].matrix
        psi = torch.zeros(2, 5, 6, dtype=pad.dtype)
        psi[:, 2, :6] = 1 / 6 ** 0.5
        r['schrodinger'] = dynamics.schrodinger(psi, H, [0.0, 0.5])
        r['schrodinger_rk4'] = dynamics.schrodinger_rk4(psi, H, [0.0, 0.5])
        Phi = channels.random_channel(pad, kraus_rank=2)
        r['channel'] = Phi(rho).matrix
        r['rhs'] = dynamics.lindblad_rhs(H, [L], rates=[0.5])(0.0, rho.matrix)
        q = TypeIAlgebra([2, 3], [2, 3], precision='double', layout=layout)
        q2 = composite.fused_tensor_product(q, q)
        torch.manual_seed(1)
        a = states.random_density_matrix(TypeIAlgebra([2, 3], [2, 3], precision='double', layout='padded'), 1)
        a = states.DensityMatrix(q, matrix=a.matrix)
        ab = composite.kron(a, a, q2)
        r['kron'] = ab.matrix
        r['ptrace'] = composite.partial_trace(ab, keep=1).matrix
        ch = SpinChain(6, 'periodic', layout=layout)
        assert ch.algebra.layout == layout
        Hc = ch.xxz(1.0, 0.7)
        g = states.gibbs_state(Hc, 0.4)
        r['reduced'] = ch.reduced_state(g, 3).matrix
        r['basis'] = ch.basis_state('110100').matrix
        jumps = [ch.sector_jump(i, 'lower') for i in range(6)]
        r['sector_jump'] = dynamics.lindblad_evolve(g, Hc, jumps, [0.0, 0.2], rates=[0.3] * 6,
                                                    method='kraus')[-1].matrix
        out[layout] = r
    for key in out['padded']:
        assert torch.allclose(out['padded'][key], out['packed'][key], atol=1e-10), key
