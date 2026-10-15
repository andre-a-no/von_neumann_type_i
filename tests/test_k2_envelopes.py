"""Exact envelopes of Exp. 1 and Exp. 3 on direct sums of M_2 (lecture notes Proposition 11.15, Corollary 11.16;
article, Proposition 2): closed forms against brute force, against the multistart optimiser, and their maximisers."""
import math

import numpy as np
import pytest
import torch

from torch_vn_algebra import TypeIAlgebra, optimize as opt

D_MAX = 0.999                                                     # cap of the free contrast, as in the article


def ell(d):
    return (1 - d) / (1 + d)


def bound1(x1, x2, y1, y2):
    """Proposition 11.15(a): sup of |Tr(XUY)| - Tr(XY) on M_2 for spectra x1 >= x2, y1 >= y2 >= 0."""
    return (x1 - x2) ** 2 * (y1 - y2) ** 2 / (4 * (x1 + x2) * (y1 + y2))


def bound3(x1, x2, y1, y2):
    """Proposition 11.15(b): sup of Tr(Y|X|Y) - Tr|YXY| on M_2 for |X| with eigenvalues x1 >= x2 >= 0."""
    return 2 * x1 * x2 / (x1 + x2) * (y1 - y2) ** 2


def env1(C, d, dp):
    """Corollary 11.16(a): Delta(X) = d, Delta(Y) = dp on C copies of M_2."""
    return C * d ** 2 * dp ** 2 / ((1 + d) * (1 + dp))


def env3(C, d, dp):
    """Corollary 11.16(b): Delta(|X|) = d, Delta(Y) = dp on C copies of M_2."""
    return (C - d) * (2 * dp / (1 + dp)) ** 2


def angles(x1, x2, y1, y2):
    """cos^2 of the angle between the eigenbases at the maximisers of Exp. 1 and Exp. 3."""
    s, q = x1 * y2 + x2 * y1, (x1 + x2) * (y1 + y2)
    return 0.25 + s / (2 * q), s / q


def rot(c):
    t = math.acos(math.sqrt(c))
    return np.array([[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]])


def nuc(A):
    return np.linalg.svd(A, compute_uv=False).sum(-1)


def z1_sup_u(X, Y):
    """sup_U z of Exp. 1 for fixed X, Y: ||YX||_1 - Tr(XY) (numpy, batched)."""
    return nuc(Y @ X) - np.einsum('...ij,...ji->...', X, Y).real


def z3_np(X, Y):
    w, V = np.linalg.eigh(X)
    aX = (V * np.abs(w)[..., None, :]) @ np.conj(np.swapaxes(V, -1, -2))
    M = Y @ X @ Y
    return np.einsum('...ij,...ji->...', Y @ aX, Y).real - np.abs(np.linalg.eigvalsh(M)).sum(-1)


SPECTRA = [(1.0, 0.3, 1.0, 0.1), (0.7, 0.2, 0.9, 0.0), (1.0, 1e-3, 1.0, 0.5), (1.0, 0.5, 0.8, 0.8),
           (2.0, 0.5, 1.5, 0.25)]


@pytest.mark.parametrize("x1,x2,y1,y2", SPECTRA)
def test_closed_forms_equal_fine_grid_maximum(x1, x2, y1, y2):
    """Prop. 11.15: the bounds equal the maxima over a fine grid of angles (sup over U exact via the trace norm)."""
    th = np.linspace(0, np.pi / 2, 20001)
    R = np.stack([np.stack([np.cos(th), -np.sin(th)], -1), np.stack([np.sin(th), np.cos(th)], -1)], -2)
    Y = R @ np.diag([y1, y2]) @ np.swapaxes(R, -1, -2)
    z1 = z1_sup_u(np.diag([x1, x2]), Y)
    z3 = z3_np(np.diag([x1, -x2]), Y)
    assert z1.max() <= bound1(x1, x2, y1, y2) + 1e-13 and z3.max() <= bound3(x1, x2, y1, y2) + 1e-13
    tol3 = 1e-8 if y2 > 0 else 2e-5                                   # y2 = 0: z3 has a kink at its maximum
    assert abs(z1.max() - bound1(x1, x2, y1, y2)) < 1e-8 and abs(z3.max() - bound3(x1, x2, y1, y2)) < tol3
    c1, c3 = angles(x1, x2, y1, y2)                                   # the grid maximisers sit at the stated angles
    if x1 > x2 and y1 > y2:
        assert abs(np.cos(th[z1.argmax()]) ** 2 - c1) < 1e-3 and abs(np.cos(th[z3.argmax()]) ** 2 - c3) < 1e-3


def test_bounds_hold_for_random_complex_pairs():
    """Prop. 11.15: no random complex X, Y, U (any phases, singular Y included) exceeds the bounds."""
    rng = np.random.default_rng(0)
    n = 20000

    def haar():
        Z = (rng.normal(size=(n, 2, 2)) + 1j * rng.normal(size=(n, 2, 2))) / math.sqrt(2)
        Q, R = map(np.stack, zip(*[np.linalg.qr(z) for z in Z]))      # batched qr needs NumPy >= 1.22
        d = np.diagonal(R, axis1=-2, axis2=-1)
        return Q * (d / np.abs(d))[:, None, :]
    x = -np.sort(-rng.random((n, 2)), -1)
    y = -np.sort(-rng.random((n, 2)), -1)
    y[: n // 5, 1] = 0.0                                              # contrast 1
    V, W, U = haar(), haar(), haar()
    dag = lambda A: np.conj(np.swapaxes(A, -1, -2))
    X = (V * x[:, None, :]) @ dag(V)
    Y = (W * y[:, None, :]) @ dag(W)
    z1 = np.abs(np.einsum('nij,nji->n', X @ U, Y)) - np.einsum('nij,nji->n', X, Y).real
    b1, b3 = bound1(x[:, 0], x[:, 1], y[:, 0], y[:, 1]), bound3(x[:, 0], x[:, 1], y[:, 0], y[:, 1])
    assert (z1 <= z1_sup_u(X, Y) + 1e-12).all() and (z1_sup_u(X, Y) <= b1 + 1e-12).all()
    Xs = (V * (x * np.array([1.0, -1.0]))[:, None, :]) @ dag(V)
    assert (z3_np(Xs, Y) <= b3 + 1e-12).all()
    assert np.abs(z3_np(X, Y)).max() < 1e-12                          # X >= 0: z = 0


@pytest.mark.parametrize("C", [1, 2, 3])
@pytest.mark.parametrize("d,dp", [(0.5, 0.8), (0.9, 0.3), (0.0, 0.7), (0.6, 1.0)])
def test_maximisers_attain_corollary(C, d, dp):
    """Cor. 11.16: the stated real maximisers on C copies of M_2 satisfy the constraints and attain the envelopes."""
    a = TypeIAlgebra([2] * C, [2] * C, complex_valued=False, precision='double')
    lx, ly = ell(d), ell(dp)
    c1, _ = angles(1.0, lx, 1.0, ly)
    X0, Y1 = np.diag([1.0, lx]), rot(c1) @ np.diag([1.0, ly]) @ rot(c1).T
    W, _, Vh = np.linalg.svd(Y1 @ X0)
    U1 = Vh.T @ W.T                                                   # Tr(X U Y) = ||Y X||_1, orthogonal
    X, Y, U = (a.from_blocks([torch.from_numpy(m)] * C) for m in (X0, Y1, U1))
    assert torch.allclose(X.michelson_contrast, torch.tensor([d], dtype=torch.float64), atol=1e-12)
    assert torch.allclose(Y.michelson_contrast, torch.tensor([dp], dtype=torch.float64), atol=1e-12)
    z = ((X @ U @ Y).trace.abs() - (X @ Y).trace.real).item()
    assert abs(z - env1(C, d, dp)) < 1e-12

    # Exp. 3: the first block carries ell(d); the others have |X| = 1 (X = diag(1, -1))
    blocks_x, blocks_y = [], []
    for c in range(C):
        x2 = lx if c == 0 else 1.0
        c3 = angles(1.0, x2, 1.0, ly)[1]
        blocks_x.append(torch.tensor([[1.0, 0.0], [0.0, -x2]], dtype=torch.float64))
        blocks_y.append(torch.from_numpy(rot(c3) @ np.diag([1.0, ly]) @ rot(c3).T))
    X, Y = a.from_blocks(blocks_x), a.from_blocks(blocks_y)
    assert torch.allclose(X.abs().michelson_contrast, torch.tensor([d], dtype=torch.float64), atol=1e-12)
    assert torch.allclose(Y.michelson_contrast, torch.tensor([dp], dtype=torch.float64), atol=1e-12)
    z = ((Y @ X.abs() @ Y).trace.real - (Y @ X @ Y).trace_norm()).item()
    assert abs(z - env3(C, d, dp)) < 1e-12


def test_envelopes_by_brute_force_over_spectra_and_angles():
    """Cor. 11.16 for C = 1, 2: maximum over a grid of spectra (with the prescribed extremes) and angles."""
    d, dp = 0.6, 0.75
    lx, ly = ell(d), ell(dp)
    g = np.linspace(0, 1, 9)
    th = np.linspace(0, np.pi / 2, 1201)
    R = np.stack([np.stack([np.cos(th), -np.sin(th)], -1), np.stack([np.sin(th), np.cos(th)], -1)], -2)

    def best(x1, x2, y1, y2):                                       # best per-block value over the angle
        Y = R @ np.diag([y1, y2]) @ np.swapaxes(R, -1, -2)
        return z1_sup_u(np.diag([x1, x2]), Y).max(), z3_np(np.diag([x1, -x2]), Y).max()
    # C = 1: the spectra are {1, lx}, {1, ly}
    b1, b3 = best(1.0, lx, 1.0, ly)
    assert abs(b1 - env1(1, d, dp)) < 1e-6 and abs(b3 - env3(1, d, dp)) < 1e-6
    # C = 2: one block holds lambda_min, the other lambda_max (possibly the same block), for X and for Y
    cache = {}

    def best_c(x, y):
        if (x, y) not in cache:
            cache[(x, y)] = best(*x, *y)
        return cache[(x, y)]
    xv = [lx + (1 - lx) * t for t in g]
    yv = [ly + (1 - ly) * t for t in g]
    best1 = best3 = -np.inf
    for u0 in xv:
        for u in xv:
            xs = [(u0, lx), (1.0, u)]
            for v0 in yv:
                for v in yv:
                    for ys in ([(v0, ly), (1.0, v)], [(1.0, v), (v0, ly)]):
                        p = [best_c(xs[c], ys[c]) for c in range(2)]
                        best1 = max(best1, p[0][0] + p[1][0])
                        best3 = max(best3, p[0][1] + p[1][1])
    assert best1 <= env1(2, d, dp) + 1e-9 and abs(best1 - env1(2, d, dp)) < 1e-6
    assert best3 <= env3(2, d, dp) + 1e-9 and abs(best3 - env3(2, d, dp)) < 1e-6


def test_lower_envelope_exp1_both_contrasts():
    """Cor. 11.16(a): inf z = -(2C - 1 + ell(d) ell(d')), attained with aligned bases and Tr(X_c U_c Y_c) = 0."""
    for C in (1, 2):
        a = TypeIAlgebra([2] * C, [2] * C, complex_valued=False, precision='double')
        d, dp = 0.5, 0.8
        X = a.from_blocks([torch.diag(torch.tensor([1.0, ell(d) if c == 0 else 1.0], dtype=torch.float64))
                           for c in range(C)])
        Y = a.from_blocks([torch.diag(torch.tensor([1.0, ell(dp) if c == 0 else 1.0], dtype=torch.float64))
                           for c in range(C)])
        U = a.from_blocks([torch.tensor([[0.0, -1.0], [1.0, 0.0]], dtype=torch.float64)] * C)
        z = ((X @ U @ Y).trace.abs() - (X @ Y).trace.real).item()
        assert abs(z + (2 * C - 1 + ell(d) * ell(dp))) < 1e-12


def _search(alg, B, kind, dx, dy=None, dy_range=(0.0, D_MAX), sign=1.0, rounds=4, steps=30):
    if kind == 1:
        X, Y, U = (opt.PositiveParam(alg, B, delta=dx), opt.PositiveParam(alg, B, delta=dy, delta_range=dy_range),
                   opt.UnitaryParam(alg, B))

        def z():
            x, y = X.operator(), Y.operator()
            return sign * ((x @ U.operator() @ y).trace.abs() - (x @ y).trace.real)
        params = [X, Y, U]
    else:
        X, Y = opt.SelfAdjointParam(alg, B, delta=dx), opt.PositiveParam(alg, B, delta=dy, delta_range=dy_range)

        def z():
            y = Y.operator()
            M = alg.operator((y @ X.operator() @ y).matrix, is_self_adjoint=True)
            return sign * ((y @ X.abs_operator() @ y).trace.real - M.trace_norm())
        params = [X, Y]
    return opt.extremize(z, params, maximize=True, rounds=rounds, steps=steps)


@pytest.mark.parametrize("kind", [1, 3])
def test_optimizer_reaches_envelope_with_free_capped_contrast(kind):
    """Cor. 11.16(c): extremize on M_2 (Delta(Y) free in [0, 0.999]) never exceeds and reaches the closed form."""
    torch.manual_seed(3)
    a = TypeIAlgebra([2], [2], complex_valued=False, precision='double')
    grid = torch.tensor([0.2, 0.5, 0.9], dtype=torch.float64)
    S = 6
    best = _search(a, 3 * S, kind, grid.repeat_interleave(S)).reshape(3, S).max(1)[0]
    exact = torch.tensor([(env1 if kind == 1 else env3)(1, d, D_MAX) for d in grid.tolist()], dtype=torch.float64)
    assert (best <= exact + 1e-10).all()
    assert ((exact - best) / exact).max() < 2e-3


@pytest.mark.parametrize("kind", [1, 3])
def test_optimizer_respects_prescribed_second_contrast(kind):
    """Cor. 11.16: with Delta(Y) prescribed the search stays below the two-parameter closed form and reaches it."""
    torch.manual_seed(4)
    a = TypeIAlgebra([2], [2], complex_valued=False, precision='double')
    S = 6
    dx = torch.tensor([0.7, 0.4], dtype=torch.float64).repeat_interleave(S)
    dy = torch.tensor([0.3, 0.8], dtype=torch.float64).repeat_interleave(S)
    best = _search(a, 2 * S, kind, dx, dy).reshape(2, S).max(1)[0]
    exact = torch.tensor([(env1 if kind == 1 else env3)(1, 0.7, 0.3), (env1 if kind == 1 else env3)(1, 0.4, 0.8)],
                         dtype=torch.float64)
    assert (best <= exact + 1e-10).all() and ((exact - best) / exact).max() < 2e-3
    # the capped free contrast gives more: the envelope increases with the cap of Delta(Y)
    assert (exact < torch.tensor([(env1 if kind == 1 else env3)(1, d, D_MAX) for d in (0.7, 0.4)],
                                 dtype=torch.float64)).all()


def test_optimizer_two_sectors_stays_below_and_close():
    """Cor. 11.16 on M_2 (+) M_2: the search never exceeds the proved envelopes and comes within 1 %."""
    torch.manual_seed(5)
    a = TypeIAlgebra([2, 2], [2, 2], complex_valued=False, precision='double')
    S = 12
    dx = torch.full((S,), 0.5, dtype=torch.float64)
    b1 = _search(a, S, 1, dx, rounds=6).max().item()
    b3 = _search(a, S, 3, dx, rounds=6).max().item()
    e1, e3 = env1(2, 0.5, D_MAX), env3(2, 0.5, D_MAX)
    assert b1 <= e1 + 1e-10 and b3 <= e3 + 1e-10
    assert e1 - b1 < 1e-2 * e1 and e3 - b3 < 1e-2 * e3


def test_contrast_one_limit_attained():
    """Cor. 11.16(c), d -> 1: with Y a rank-one projection per block (Delta(Y) = 1) the limits are attained."""
    for C in (1, 2):
        a = TypeIAlgebra([2] * C, [2] * C, complex_valued=False, precision='double')
        d = 0.4
        lx = ell(d)
        c1, c3 = angles(1.0, lx, 1.0, 0.0)
        X0, Y1 = np.diag([1.0, lx]), rot(c1) @ np.diag([1.0, 0.0]) @ rot(c1).T
        W, _, Vh = np.linalg.svd(Y1 @ X0)
        X, Y, U = (a.from_blocks([torch.from_numpy(m)] * C) for m in (X0, Y1, Vh.T @ W.T))
        assert abs(Y.michelson_contrast.item() - 1) < 1e-12
        z = ((X @ U @ Y).trace.abs() - (X @ Y).trace.real).item()
        assert abs(z - C * d ** 2 / (2 * (1 + d))) < 1e-12
        bx = [np.diag([1.0, -(lx if c == 0 else 1.0)]) for c in range(C)]
        by = [rot(angles(1.0, lx if c == 0 else 1.0, 1.0, 0.0)[1]) @ np.diag([1.0, 0.0])
              @ rot(angles(1.0, lx if c == 0 else 1.0, 1.0, 0.0)[1]).T for c in range(C)]
        X, Y = a.from_blocks([torch.from_numpy(m) for m in bx]), a.from_blocks([torch.from_numpy(m) for m in by])
        z = ((Y @ X.abs() @ Y).trace.real - (Y @ X @ Y).trace_norm()).item()
        assert abs(z - (C - d)) < 1e-12
    # and below contrast 1 the value is strictly smaller and increases to the limit
    caps = np.array([0.9, 0.99, 0.999, 0.9999])
    v1, v3 = env1(1, 0.4, caps), env3(1, 0.4, caps)
    assert (np.diff(v1) > 0).all() and (v1 < 0.4 ** 2 / 2.8).all() and (np.diff(v3) > 0).all() and (v3 < 0.6).all()
