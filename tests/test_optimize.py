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


def test_pin_sector_skips_one_dimensional_sectors():
    from torch_vn_algebra import SpinChain
    alg = SpinChain(4, complex_valued=False).algebra                 # sector N=0 has dimension 1
    X = opt.PositiveParam(alg, 3, delta=torch.full((3,), 0.4, dtype=torch.float64))
    assert torch.allclose(X.operator().michelson_contrast, torch.full((3,), 0.4, dtype=torch.float64))


def test_extreme_eigenvalues_placed_in_any_sector():
    a = TypeIAlgebra([2, 2], [2, 2], complex_valued=False, precision='double')
    S = opt.SpectrumParam(a, 200, delta=torch.full((200,), 0.5, dtype=torch.float64))
    different = (S.pin_max[:, 0] != S.pin_min[:, 0])
    assert different.any() and (~different).any()                  # both placements occur among the starts
    X = opt.PositiveParam(a, 200, delta=torch.full((200,), 0.5, dtype=torch.float64))
    assert torch.allclose(X.operator().michelson_contrast, torch.full((200,), 0.5, dtype=torch.float64))
    P = opt.SpectrumParam(a, 5, delta=torch.full((5,), 0.5, dtype=torch.float64), pin_sector=1)
    assert (P.pin_max[:, 0] == 1).all() and (P.pin_min[:, 0] == 1).all()


def test_sign_patterns_cover_balanced_sectors():
    # the number of negative eigenvalues is drawn per sector, or one fraction for all sectors: a pattern with one
    # negative eigenvalue in every one of 8 sectors of size 2 must be common (iid signs: probability 2^-8)
    torch.manual_seed(0)
    S = opt.SelfAdjointParam(alg(k=2, C=8), 400)
    balanced = ((S.signs < 0).sum(-1) == 1).all(-1).float().mean().item()
    assert balanced > 0.1


def test_exp3_sign_search_reaches_direct_sum():
    # Exp. 3 at Delta(|X|) = 0: |X| = 1, z = Tr Y^2 - ||Y X Y||_1; a direct sum of C optimal 2x2 blocks gives z -> C
    # (Y -> projection onto (e1 + e2)/sqrt2 per block, X = diag(1, -1)); it needs a mixed sign pattern in every sector
    torch.manual_seed(0)
    a = alg(k=2, C=4)
    d = torch.zeros(32, dtype=torch.float64)
    X, Y = opt.SelfAdjointParam(a, 32, delta=d), opt.PositiveParam(a, 32)

    def z():
        x, y = X.operator(), Y.operator()
        return (y @ X.abs_operator() @ y).trace.real - (y @ x @ y).trace_norm()
    v = opt.extremize(z, [X, Y], maximize=True, rounds=6, steps=40)
    assert v.max().item() > 0.97 * 4
