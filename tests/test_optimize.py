"""Tests for the constrained optimisation helpers."""
import torch

from torch_vn_algebra import TypeIAlgebra, optimize as opt


def alg(k=3, C=2):
    return TypeIAlgebra([k] * C, [k] * C, complex_valued=False, precision='double', device='cpu')


def test_contrast_constraint_is_exact():
    a = alg()
    d = torch.linspace(0, 0.95, 16, dtype=torch.float64)
    X = opt.PositiveParam(a, 16, delta=d)
    assert torch.allclose(X.operator().michelson_contrast, d, atol=1e-12)
    S = opt.SelfAdjointParam(a, 16, delta=d)
    assert torch.allclose(S.abs_operator().michelson_contrast, d, atol=1e-12)
    assert torch.allclose(S.operator().abs().matrix, S.abs_operator().matrix, atol=1e-10)


def test_free_contrast_stays_in_range():
    a = alg()
    X = opt.PositiveParam(a, 32, delta_range=(0.2, 0.4))
    c = X.operator().michelson_contrast
    assert torch.all((c >= 0.2 - 1e-12) & (c <= 0.4 + 1e-12))


def test_unitary_recenter_keeps_point():
    a = alg()
    U = opt.UnitaryParam(a, 4)
    with torch.no_grad():
        U.G.normal_()
    before = U.operator().matrix.detach()
    U.recenter()
    assert torch.allclose(U.operator().matrix, before, atol=1e-12)
    M = before
    eye = a.identity(4).matrix
    assert torch.allclose(M @ M.transpose(-2, -1), eye, atol=1e-10)


def test_exp1_extremes_at_zero_contrast():
    # Delta = 0: X = Y = 1, z = |Tr U| - kC, sup 0 (U = 1), inf -kC (Tr U = 0)
    torch.manual_seed(0)
    a = alg(k=2, C=2)
    zero = torch.zeros(16, dtype=torch.float64)

    def make():
        X, Y, U = opt.PositiveParam(a, 16, delta=zero), opt.PositiveParam(a, 16, delta=zero), opt.UnitaryParam(a, 16)
        z = lambda: (X.operator() @ U.operator() @ Y.operator()).trace.abs() - (X.operator() @ Y.operator()).trace.real
        return z, [X, Y, U]
    z, p = make()
    assert abs(opt.extremize(z, p, maximize=True, rounds=4, steps=30).max().item()) < 1e-3
    z, p = make()
    assert abs(opt.extremize(z, p, maximize=False, rounds=4, steps=30).min().item() + 4) < 1e-2


def test_ascent_reaches_known_supremum():
    torch.manual_seed(1)
    a = TypeIAlgebra([6], [6], precision='double', device='cpu')
    X = a.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
    Y = a.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
    U = opt.UnitaryParam(a, 16)
    trXY = (X @ Y).trace.real
    z = lambda: (X @ U.operator() @ Y).trace.abs() - trXY
    best = opt.extremize(z, [U], rounds=8, steps=50).max().item()
    sup = ((Y @ X).trace_norm() - trXY).item()
    assert sup > 0 and abs(best - sup) < 1e-3 * abs(trXY.item())
