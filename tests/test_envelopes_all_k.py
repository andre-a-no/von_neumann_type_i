"""Exact upper envelopes of Exp. 1 and Exp. 3 for every sector size (scripts/certify_envelopes.py, Lemma 6, Theorems 1,
4, 7, Corollary 8; paper, Theorems on the envelopes; lecture notes, Section 11.7): linearisation, certificate, bounds,
maximisers and the multistart search."""
import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest
import torch

from torch_vn_algebra import TypeIAlgebra, optimize as opt

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location('certify_envelopes', ROOT / 'scripts' / 'certify_envelopes.py')
ce = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ce)

D = ce.DMAX


def ell(d):
    return (1 - d) / (1 + d)


def haar(rng, n, complex_valued=False):
    A = rng.normal(size=(n, n)) + (1j * rng.normal(size=(n, n)) if complex_valued else 0)
    Q, R = np.linalg.qr(A)
    return Q * (np.diag(R) / np.abs(np.diag(R)))


def tn(A):
    return np.linalg.svd(A, compute_uv=False).sum(-1)


def z1_sup_u(X, Y):
    """sup over unitaries U of Exp. 1: ||YX||_1 - Tr(XY) (trace-norm duality)."""
    return tn(Y @ X) - np.trace(X @ Y).real


def z3(X, Y):
    """Exp. 3 computed directly: Tr(Y|X|Y) - ||YXY||_1."""
    w, V = np.linalg.eigh(X)
    absX = (V * np.abs(w)) @ V.conj().T
    return np.trace(Y @ absX @ Y).real - tn(Y @ X @ Y)


def blockdiag(blocks):
    n = sum(b.shape[0] for b in blocks)
    out = np.zeros((n, n), dtype=np.result_type(*blocks))
    i = 0
    for b in blocks:
        out[i:i + len(b), i:i + len(b)] = b
        i += len(b)
    return out


def rot(t):
    return np.array([[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]])


def test_linearisation_has_the_spectrum_of_yxy():
    """Lemma 6: Y X Y and the matrix M(A), affine in A = diag(cos^2 theta), have the same eigenvalues (b = 0 too)."""
    rng = np.random.default_rng(0)
    for trial in range(200):
        m = 1 + trial % 4
        b = 0.0 if trial % 5 == 0 else rng.uniform(0, 1)
        W = haar(rng, m, complex_valued=trial % 2 == 1)
        Xp = (W * rng.uniform(0.01, 1, m)) @ W.conj().T
        theta = rng.uniform(0, np.pi / 2, m)
        X, Y = ce.block_model_operators(theta, Xp, b)
        M, _, _ = ce.linearisation(theta, Xp, b)
        assert np.allclose(np.linalg.eigvalsh(Y @ X @ Y), np.linalg.eigvalsh(M), atol=1e-12)


def test_theorem7_certificate_is_a_symmetry_with_constant_value():
    """Theorem 7: T^2 = 1 and Tr(|X| Y^2) - Tr(T M(A)) = (1 - b)^2 sum 2 x_i/(1 + x_i) for every angle."""
    rng = np.random.default_rng(1)
    for trial in range(100):
        m = 1 + trial % 4
        b = rng.uniform(0, 1)
        W = haar(rng, m, complex_valued=True)
        x = rng.uniform(0.01, 1, m)
        Xp = (W * x) @ W.conj().T
        theta = rng.uniform(0, np.pi / 2, m)
        _, T, value = ce.linearisation(theta, Xp, b)
        assert np.allclose(T @ T, np.eye(2 * m), atol=1e-12)
        assert value == pytest.approx(ce.blocks_bound(x, b), abs=1e-12)
        M, _, _ = ce.linearisation(theta, Xp, b)
        _, Y = ce.block_model_operators(theta, Xp, b)                         # Tr(|X| Y^2) computed directly
        assert np.trace(blockdiag([Xp, np.eye(m)]) @ Y @ Y).real == pytest.approx(
            value + np.trace(T @ M).real, abs=1e-12)


def test_theorem7_bound_and_equality():
    """Theorem 7: z <= (1 - b)^2 sum 2 x_i/(1 + x_i) on the m-block model, with equality at the stated angles."""
    rng = np.random.default_rng(2)
    worst = -np.inf
    for trial in range(2000):
        m = 1 + trial % 4
        b = 0.0 if trial % 7 == 0 else rng.uniform(0, 1)
        W = haar(rng, m, complex_valued=trial % 2 == 0)
        x = rng.uniform(0.01, 1, m)
        if trial % 3 == 0:                                    # the case of the envelope: x = (a, 1, ..., 1)
            x[1:] = 1.0
        Xp = (W * x) @ W.conj().T
        theta = rng.uniform(0, np.pi / 2, m)
        X, Y = ce.block_model_operators(theta, Xp, b)
        worst = max(worst, z3(X, Y) - ce.blocks_bound(x, b))
    assert worst <= 1e-12
    for x, b in (([0.2, 1.0, 0.7], 0.3), ([0.05, 1.0, 1.0, 1.0], 1 / 1999), ([0.6], 0.0)):
        X, Y = ce.block_model_operators(ce.blocks_maximiser(x, b), np.diag(x), b)
        assert z3(X, Y) == pytest.approx(ce.blocks_bound(x, b), abs=1e-12)


def random_config(rng, ks, dx, dy, exp, structured):
    """A configuration with Delta(X) (|X| in Exp. 3) = dx and Delta(Y) = dy exactly, on sectors of sizes ks: Haar
    eigenbases and uniform spectra, or two-level spectra with eigenbases close to the maximisers."""
    n = sum(ks)
    a, b = ell(dx), ell(dy)
    if structured:
        x = np.where(rng.random(n) < 0.8, 1.0, a)
        y = np.where(rng.random(n) < 0.5, 1.0, b)
    else:
        x, y = rng.uniform(a, 1, n), rng.uniform(b, 1, n)
    x[rng.permutation(n)[:2]] = 1.0, a
    y[rng.permutation(n)[:2]] = 1.0, b
    V, W = [], []
    for k in ks:
        if structured:                                        # a signed permutation times 2x2 rotations, perturbed
            R = np.eye(k)
            for j in range(k // 2):
                R[2 * j:2 * j + 2, 2 * j:2 * j + 2] = rot(rng.uniform(0, np.pi))
            P = np.eye(k)[rng.permutation(k)] * rng.choice([-1.0, 1.0], k)
            G = rng.normal(size=(k, k)) * 0.05
            V.append(np.eye(k))
            W.append(P @ R @ np.linalg.qr(np.eye(k) + G - G.T)[0])
        else:
            V.append(haar(rng, k, True))
            W.append(haar(rng, k, True))
    V, W = blockdiag(V), blockdiag(W)
    s = np.where(rng.random(n) < 0.5, -1.0, 1.0) if exp == 3 else np.ones(n)
    return (V * (s * x)) @ V.conj().T, (W * y) @ W.conj().T


@pytest.mark.parametrize('ks', [[2], [3], [4], [5], [6], [2, 2], [3, 3], [4, 4], [2, 3], [4, 2], [2, 2, 2]])
def test_no_configuration_exceeds_the_envelopes(ks):
    """Theorem 1 and Corollary 8: random and near-extremal configurations stay below the exact suprema, for both
    contrasts prescribed, and so below the envelopes along either contrast with the other free in [0, D]."""
    rng = np.random.default_rng(sum(ks) * 10 + len(ks))
    for trial in range(150):
        dx, dy = rng.uniform(0, 0.99), rng.uniform(0, D)
        structured = trial % 2 == 0
        X, Y = random_config(rng, ks, dx, dy, 1, structured)
        assert z1_sup_u(X, Y) <= ce.exact_sup(1, ks, dx, other=dy) + 1e-12
        X, Y = random_config(rng, ks, dx, dy, 3, structured)
        v = z3(X, Y)
        assert v <= ce.exact_sup(3, ks, dx, other=dy) + 1e-12
        assert v <= ce.exact_sup(3, ks, dx, D=dy) + 1e-12                     # along Delta(|X|), Delta(Y) <= dy
        assert v <= ce.exact_sup(3, ks, dy, D=dx, along='y') + 1e-12          # along Delta(Y), Delta(|X|) <= dx


def maximisers(ks, dx, dy):
    """Exp. 1 (X, Y, U) and Exp. 3 (X, Y) attaining Theorem 1 and Corollary 8 on sectors of sizes ks (real)."""
    a, b = ell(dx), ell(dy)
    even = all(k % 2 == 0 for k in ks)
    c1 = (1 + a * b + 3 * (a + b)) / (4 * (1 + a) * (1 + b))                  # Theorem 1 (iii)
    X1, Y1, X3, Y3 = [], [], [], []
    a_used = False
    for k in ks:
        x1, y1, x3, y3 = [], [], [], []
        for j in range(k // 2):
            R = rot(math.acos(math.sqrt(c1)))
            x1.append(np.diag([1.0, a]))
            y1.append(R @ np.diag([1.0, b]) @ R.T)
            xa = a if even and not a_used else 1.0
            a_used = a_used or xa == a
            R = rot(float(ce.blocks_maximiser([xa], b)[0]))
            x3.append(np.diag([xa, -1.0]))
            y3.append(R @ np.diag([1.0, b]) @ R.T)
        if k % 2:
            x1.append(np.eye(1))
            y1.append(np.eye(1))
            x3.append(np.eye(1) * (1.0 if a_used else a))
            a_used = True
            y3.append(np.eye(1))
        X1 += x1
        Y1 += y1
        X3 += x3
        Y3 += y3
    X1, Y1, X3, Y3 = map(blockdiag, (X1, Y1, X3, Y3))
    U, _, Vh = np.linalg.svd(Y1 @ X1)                                          # Y X = U S V^T, polar unitary
    return (X1, Y1, Vh.T @ U.T), (X3, Y3)


@pytest.mark.parametrize('ks', [[2], [3], [4], [16], [2, 2], [3, 3], [4, 4, 4], [2, 3], [5, 4]])
def test_maximisers_attain_the_envelopes(ks):
    """Theorem 1 and Corollary 8 are attained by direct sums of 2 x 2 maximisers (with the eigenvalue l(delta) of |X|
    on a 1-dimensional block when some sector is odd); values computed by the library and by NumPy."""
    alg = TypeIAlgebra(ks, ks, complex_valued=False, precision='double')
    starts = np.cumsum([0] + ks)
    for dx, dy in ((0.5, 0.8), (0.9, 0.999), (0.3, 1.0)):
        (X1, Y1, U1), (X3, Y3) = maximisers(ks, dx, dy)
        ops = [alg.from_blocks([torch.from_numpy(np.ascontiguousarray(M[s:e, s:e]))
                                for s, e in zip(starts[:-1], starts[1:])]) for M in (X1, Y1, U1, X3, Y3)]
        x1, y1, u1, x3, y3 = ops
        assert abs(x1.michelson_contrast.item() - dx) < 1e-12 and abs(y1.michelson_contrast.item() - dy) < 1e-12
        assert abs(x3.abs().michelson_contrast.item() - dx) < 1e-12 and abs(y3.michelson_contrast.item() - dy) < 1e-12
        v1 = ((x1 @ u1 @ y1).trace.abs() - (x1 @ y1).trace.real).item()
        v3 = ((y3 @ x3.abs() @ y3).trace.real - (y3 @ x3 @ y3).trace_norm()).item()
        assert v1 == pytest.approx(float(ce.exact_sup(1, ks, dx, other=dy)), abs=1e-12)
        assert v3 == pytest.approx(float(ce.exact_sup(3, ks, dx, other=dy)), abs=1e-12)
        assert z1_sup_u(X1, Y1) == pytest.approx(v1, abs=1e-12) and z3(X3, Y3) == pytest.approx(v3, abs=1e-12)


def test_envelope_formulas():
    """The closed forms of exact_sup: k = 2 agrees with Proposition 2 / Corollary 11.16, odd k does not depend on
    delta, the envelopes along Delta(|X|) decrease and those along Delta(Y) increase."""
    g = np.linspace(0, 0.99, 12)
    for C in (1, 2, 16, 32):
        assert np.allclose(ce.exact_sup(1, 2, g, C=C), C * g ** 2 * D ** 2 / ((1 + g) * (1 + D)))
        assert np.allclose(ce.exact_sup(3, 2, g, C=C), (C - g) * (2 * D / (1 + D)) ** 2)
        assert np.allclose(ce.exact_sup(3, 2, g, C=C, along='y'), C * (2 * g / (1 + g)) ** 2)
        assert np.allclose(ce.exact_sup(3, 16, g, C=C), (8 * C - g) * (2 * D / (1 + D)) ** 2)
        assert np.allclose(ce.exact_sup(3, 3, g, C=C), C * (2 * D / (1 + D)) ** 2)
        assert np.all(np.diff(ce.exact_sup(3, 4, g, C=C)) < 0)
        assert np.all(np.diff(ce.exact_sup(1, 4, g, C=C, along='y')) > 0)
    assert np.allclose(ce.exact_sup(3, 4, g), ce.direct_sum(3, 4, g))


def search(exp, ks, deltas, along, starts, seed, rounds=6, steps=40, decay=0.6):
    """The multistart search of scripts/bounds_search.py (real algebra, the other contrast free in [0, D])."""
    torch.manual_seed(seed)
    alg = TypeIAlgebra(ks, ks, complex_valued=False, precision='double')
    d = torch.tensor(deltas, dtype=torch.float64).repeat_interleave(starts)
    B = len(d)
    dx, dy = (d, None) if along == 'x' else (None, d)
    cls = opt.SelfAdjointParam if exp == 3 else opt.PositiveParam
    X = cls(alg, B, delta=dx, init='diagonal', init_fraction=0.5,            # half of the starts structured
            spectrum_init='high' if exp == 3 else 'pairs')
    Y = opt.PositiveParam(alg, B, delta=dy, init='pairs', init_fraction=X.structured, spectrum_init='pairs')
    if exp == 1:
        U = opt.UnitaryParam(alg, B, init='rotations', init_fraction=X.structured)

        def z():
            x, y = X.operator(), Y.operator()
            return (x @ U.operator() @ y).trace.abs() - (x @ y).trace.real
        params = [X, Y, U]
    else:
        def z():
            x, y = X.operator(), Y.operator()
            M = y @ x @ y
            M._is_self_adjoint = True
            return (y @ X.abs_operator() @ y).trace.real - M.trace_norm()
        params = [X, Y]
    v = opt.extremize(z, params, maximize=True, rounds=rounds, steps=steps, decay=decay).reshape(len(deltas), starts)
    return torch.nan_to_num(v, nan=-math.inf).max(1)[0].numpy()


@pytest.mark.filterwarnings('ignore:optimize.extremize')
@pytest.mark.parametrize('exp,ks,along', [(3, [4], 'x'), (3, [2] * 8, 'x'), (3, [3, 3], 'x'), (3, [4], 'y'),
                                          (1, [4], 'x'), (1, [2] * 8, 'x'), (1, [3], 'y')])
def test_search_reaches_the_envelopes(exp, ks, along):
    """The multistart search (structured initialisation) never exceeds Theorem 1 / Corollary 8 and reaches them to
    3e-3 relative (numerical observation for these seeds; the last digits depend on the BLAS library)."""
    deltas = [0.3, 0.8]
    found = search(exp, ks, deltas, along, starts=8, seed=len(ks) + exp)
    ex = ce.exact_sup(exp, ks, np.array(deltas), along=along)
    assert np.all(found <= ex + 1e-9)
    assert np.all((ex - found) / ex < 3e-3), (found, ex)


@pytest.mark.parametrize('ks,cv', [([2, 3, 4], False), ([3, 3], True), ([5], False)])
def test_structured_unitary(ks, cv):
    """optimize.structured_unitary: unitary blocks (zero padding); 'diagonal' is a signed identity, 'oriented' and
    'rotations' have det +1 per sector in a real algebra, and 'pairs' is a direct sum of 2 x 2 rotations of consecutive
    basis vectors (checked entrywise)."""
    torch.manual_seed(0)
    alg = TypeIAlgebra(ks, ks, complex_valued=cv, precision='double')
    eye = alg.identity(5).matrix
    for kind in opt.INITS:
        U = opt.structured_unitary(alg, 5, kind)
        assert torch.allclose(U @ U.mH, eye.to(U.dtype), atol=1e-12)
    D = opt.structured_unitary(alg, 5, 'diagonal')
    assert torch.equal(D.abs(), eye.to(D.dtype).abs())
    if not cv:
        for kind in ('oriented', 'rotations'):
            O = opt.structured_unitary(alg, 5, kind)
            for c, k in enumerate(ks):
                assert (torch.linalg.det(O[:, c, :k, :k]) > 0).all()
    R = opt.structured_unitary(alg, 5, 'pairs').real.numpy()
    for b in range(5):
        for c, k in enumerate(ks):
            Rc = R[b, c, :k, :k]
            # each column is supported on one pair {2i, 2i + 1} (or on the last slot of an odd sector)
            for col in Rc.T:
                rows = np.nonzero(np.abs(col) > 1e-14)[0]
                assert len(rows) <= 2 and (len(rows) < 2 or (rows[0] // 2 == rows[1] // 2 and rows[1] < k))
    with pytest.raises(ValueError):
        opt.structured_unitary(alg, 2, 'sobol')


def test_default_initialisation_is_unchanged():
    """init='haar' (the default) draws exactly the random numbers of the original parametrisation: same seed, same
    operators as SpectrumParam followed by a Haar UnitaryParam."""
    alg = TypeIAlgebra([3, 2], [3, 2], complex_valued=False, precision='double')
    d = torch.tensor([0.2, 0.7], dtype=torch.float64)
    torch.manual_seed(5)
    X = opt.SelfAdjointParam(alg, 2, delta=d)
    torch.manual_seed(5)
    sp = opt.SpectrumParam(alg, 2, d)
    U = opt.UnitaryParam(alg, 2)
    assert torch.equal(X.spectrum.theta, sp.theta) and torch.equal(X.unitary.base, U.base)
    assert not X.structured.any()


def test_structured_starts_and_spectra():
    """init_fraction as a mask selects the same starts for X and Y; structured self-adjoint starts have one + and one
    - sign in every pair of slots; spectrum_init='pairs' gives one eigenvalue near 1 and one near lambda_min per pair,
    'high' all near 1 (besides the pinned lambda_min)."""
    torch.manual_seed(3)
    alg = TypeIAlgebra([4, 4], [4, 4], complex_valued=False, precision='double')
    d = torch.full((40,), 0.6, dtype=torch.float64)
    X = opt.SelfAdjointParam(alg, 40, delta=d, init='diagonal', init_fraction=0.5, spectrum_init='high')
    Y = opt.PositiveParam(alg, 40, init='pairs', init_fraction=X.structured, spectrum_init='pairs')
    m = X.structured
    assert 0 < m.sum() < 40 and torch.equal(Y.structured, m)
    s = X.signs[m]                                                            # (n, C, 4)
    assert torch.all(s[..., 0::2] * s[..., 1::2] == -1)
    lo = float(ell(0.6))
    lx = X.spectrum.eigenvalues()[m]
    assert ((lx > 0.9) | ((lx - lo).abs() < 1e-12)).all()
    th = Y.spectrum.theta.detach()[m]                                       # unpinned slots: (1 + sin(+-1.5))/2
    assert torch.all(th[..., 0::2] * th[..., 1::2] == -1.5 ** 2)
    assert torch.all(Y.spectrum.eta.detach()[m] == 1.5)                      # free contrast near its upper end
    assert torch.equal(X.unitary.base[m], X.unitary.base[m].diagonal(dim1=-2, dim2=-1).diag_embed())
