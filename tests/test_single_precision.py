"""
Single precision (float32 / complex64) against double precision for the public API.

Every function is run on single-precision data and compared with the double-precision result for the *same*
(rounded) input, so that the difference is the error of the single-precision computation alone. The error is bounded
by a modest constant times its first-order perturbation bound, with u = eps(float32) = 1.19e-7, k the block size and a
backward error ||E|| <= k u ||A|| of the batched LAPACK / cuSOLVER kernels:

    eigenvalues (Weyl)                  |d lambda| <= ||E||
    spectral projectors (Davis-Kahan)   ||dP|| <= ||E|| / gap
    f(A), A self-adjoint (Daleckii-Krein) ||df(A)||_F <= max_ij |f[lambda_i, lambda_j]| ||E||_F  (+ k u max|f|)
    sqrt near 0 (operator monotone)     ||sqrt(A + E) - sqrt(A)|| <= ||E||^{1/2}
    |X| (Araki-Yamagami)                ||d|X|||_F <= sqrt(2) ||E||_F
    X^{-1}                              ||dX^{-1}|| <= ||X^{-1}||^2 ||E||
    scalar functionals F(A)             |dF| <= ||grad F||_1 ||E||_op  (trace norm of the gradient, in double)
    entropy near 0 (Fannes)             |S(rho + E) - S(rho)| <= sum_i eta(||E||), eta(t) = -t log t
    sums / products (Higham)            |fl(sum a_i b_i) - sum a_i b_i| <= n u sum |a_i| |b_i|

The constants are tight enough to catch a silent cast to half precision (u = 1e-3), an unstable formula (an error
growing with a condition number the bound does not contain) or a double-only code path (which fails, returns NaN, or
returns double tensors). Inputs include nearly degenerate spectra, large dynamic ranges, eigenvalues at 0 for
log / sqrt / entropies, and very large and very small norms; single-precision results must never be NaN / Inf where
the double ones are finite and representable in float32. The results must stay in single precision.
"""
import math
import os
import warnings

import numpy as np
import pytest
import torch

import torch_vn_algebra as vn
from torch_vn_algebra import TypeIAlgebra, DensityMatrix
from torch_vn_algebra import channels as ch
from torch_vn_algebra import composite, decompose as dc, divergences as dv, dynamics as dyn, functional as fc
from torch_vn_algebra import inclusions as incl, jordan, krylov, optimize, states as st
from torch_vn_algebra.chains import SpinChain
from torch_vn_algebra.su2 import SU2Chain

pytestmark = [pytest.mark.filterwarnings("ignore:.*not trace preserving"),    # double channels from rounded Kraus
              pytest.mark.filterwarnings("ignore:the inclusion is not connected")]

U = torch.finfo(torch.float32).eps                 # 1.19e-7
F32MAX = torch.finfo(torch.float32).max
SINGLE = (torch.float32, torch.complex64)
REPORT = bool(os.environ.get("VN_PRECISION_REPORT"))
NS = [4, 3, 1]                                     # sectors of different size: padding and a 1 x 1 block
KMAX = max(NS)


# ----------------------------------------------------------------------
# comparison helpers
# ----------------------------------------------------------------------
def d64(t):
    """A tensor (or Operator) as a double / complex double tensor."""
    t = t.matrix if isinstance(t, vn.Operator) else torch.as_tensor(t)
    return t.to(torch.complex128) if torch.is_complex(t) else t.to(torch.float64)


def assert_single(*ts):
    for t in ts:
        t = t.matrix if isinstance(t, vn.Operator) else t
        assert t.dtype in SINGLE, f"result in {t.dtype}, expected single precision"


def assert_finite_where(got, ref):
    """No NaN / Inf in the single result where the double result is finite and representable in float32."""
    got, ref = d64(got), d64(ref)
    ok = torch.isfinite(ref) & (ref.abs() < 1e-3 * F32MAX)
    assert torch.isfinite(got[ok]).all(), "single precision gives NaN/Inf where double precision is finite"


def check(err, bound, C, what=""):
    """err <= C * bound (tensors broadcast; per sample)."""
    err = torch.as_tensor(err).detach().to(torch.float64)
    bound = torch.as_tensor(bound).detach().to(torch.float64).clamp(min=1e-300)
    assert torch.isfinite(err).all(), f"{what}: non-finite error"
    ratio = float((err / bound).max())
    if REPORT:
        print(f"\nRATIO {what}: {ratio:.3g} (C = {C})")
    assert ratio <= C, f"{what}: error / bound = {ratio:.3g} > {C}"


def fro(t):
    """Frobenius norm per sample of (B, ...) tensors."""
    return torch.linalg.vector_norm(t.reshape(t.shape[0], -1), dim=-1)


def compare(got, ref, bound, C, what, norm=fro):
    """Single result `got` against the double result `ref` of the same rounded input: dtype, finiteness, error."""
    assert_single(got)
    assert_finite_where(got, ref)
    check(norm(d64(got) - d64(ref)), bound, C, what)


# ----------------------------------------------------------------------
# inputs
# ----------------------------------------------------------------------
def eye_of(alg, dtype=None):
    """The unit of the algebra as a padded (C, k_max, k_max) tensor (public API: alg.identity())."""
    e = alg.identity().matrix[0]
    return e if dtype is None else e.to(dtype)


def rows_of(alg):
    """(C, k_max) bool: the active rows of every sector."""
    return torch.diagonal(alg.identity().matrix[0], dim1=-2, dim2=-1).real > 0


def mask_of(alg):
    """(C, k_max, k_max) bool: the active blocks (the padding is False)."""
    r = rows_of(alg)
    return r[:, :, None] & r[:, None, :]


def algs(ns=NS, complex_valued=True):
    a64 = TypeIAlgebra(ns, ns, complex_valued=complex_valued, precision='double')
    return a64.like(precision='single'), a64


def haar_np(k, rng, complex_=True):
    Z = rng.standard_normal((k, k)) + (1j * rng.standard_normal((k, k)) if complex_ else 0)
    Q, R = np.linalg.qr(Z)
    d = np.diagonal(R)
    return Q * (d / np.abs(d))


def spectra_matrix(alg, spectra, B=3, seed=0, complex_=True):
    """(B, C, k_max, k_max) double tensor with Hermitian blocks U diag(w) U^*, w = spectra(k, rng)."""
    rng = np.random.default_rng(seed)
    M = np.zeros((B, alg.C, alg.k_max, alg.k_max), dtype=complex if complex_ else float)
    for b in range(B):
        for c, k in enumerate(alg.k_factors):
            Uk = haar_np(k, rng, complex_)
            M[b, c, :k, :k] = (Uk * np.asarray(spectra(k, rng), dtype=float)) @ Uk.conj().T
    M = 0.5 * (M + np.conj(np.swapaxes(M, -1, -2)))
    return torch.as_tensor(M)


def general_matrix(alg, B=3, seed=0, scale=1.0, shift=0.0, complex_=True):
    """Non-normal blocks (Ginibre, scale / sqrt(k)) plus shift * 1, padding zero."""
    gen = torch.Generator().manual_seed(seed)
    dt = torch.complex128 if complex_ else torch.float64
    M = torch.randn(B, alg.C, alg.k_max, alg.k_max, dtype=dt, generator=gen) * scale / math.sqrt(alg.k_max)
    return (M + shift * eye_of(alg, dt)) * mask_of(alg)


def pair(a32, a64, M, **flags):
    """(single operator, double operator of the same rounded data)."""
    m32 = M.to(a32.dtype)
    return a32.operator(m32, **flags), a64.operator(m32.to(a64.dtype), **flags)


def pair_state(a32, a64, M):
    m32 = M.to(a32.dtype)
    return (DensityMatrix(a32, matrix=m32, validate=False),
            DensityMatrix(a64, matrix=m32.to(a64.dtype), validate=False))


def herm_eigs(op64):
    """Eigenvalues (double) of the active blocks per sample, concatenated: (B, sum k_c)."""
    return torch.cat([w for w in op64.eigenvalues() if w.numel()], dim=-1)


def sector_eigs(op64):
    return [w for w in op64.eigenvalues()]


def op_norm(op64):
    """Largest singular value per sample (double)."""
    return op64.operator_norm().to(torch.float64)


def loewner_max(w, f, df):
    """max_ij |f[w_i, w_j]| per row of eigenvalues w (B, k) (f' on (nearly) coinciding values)."""
    a, b = w.unsqueeze(-1), w.unsqueeze(-2)
    close = (a - b).abs() <= 1e-7 * torch.maximum(a.abs(), b.abs()).clamp(min=1e-300)
    q = (f(a) - f(b)) / torch.where(close, torch.ones_like(a), a - b)
    mid = df(0.5 * (a + b))
    return torch.where(close, mid, q).abs().amax(dim=(-2, -1))


def dk_bound(op64, f, df):
    """k u (||A|| max|L| + max|f(lambda)|) per sample, the Daleckii-Krein bound for the backward error k u ||A||."""
    out = torch.zeros(op64.matrix.shape[0], dtype=torch.float64)
    nrm = op_norm(op64)
    for w in sector_eigs(op64):
        if not w.numel():
            continue
        k = w.shape[-1]
        out = torch.maximum(out, KMAX * U * (nrm * loewner_max(w, f, df) + f(w).abs().amax(-1)) * math.sqrt(k))
    return out


# spectral families: (k, rng) -> eigenvalues
def _pairs(g, lo=0.5, hi=2.0):
    def s(k, rng):
        base = rng.uniform(lo, hi, (k + 1) // 2)
        return np.stack([base, base * (1 + g)], -1).reshape(-1)[:k] * rng.choice([-1, 1], 1)
    return s


def _pos_pairs(g):
    def s(k, rng):
        base = rng.uniform(0.5, 2.0, (k + 1) // 2)
        return np.stack([base, base * (1 + g)], -1).reshape(-1)[:k]
    return s


FAMILIES = {
    'generic': lambda k, rng: rng.uniform(-1.0, 2.0, k),
    'gap1e-6': _pairs(1e-6),
    'gap1e-9': _pairs(1e-9),
    'gap1e-12': _pairs(1e-12),
    'large': lambda k, rng: rng.uniform(-1.0, 2.0, k) * 1e15,
    'small': lambda k, rng: rng.uniform(-1.0, 2.0, k) * 1e-15,
    'range': lambda k, rng: rng.permutation(np.logspace(-3, 3, k)) * rng.choice([-1, 1], k),
    'psd_zero': lambda k, rng: np.concatenate([[0.0], np.logspace(-6, 0, k - 1)])[:k],
}
POS_FAMILIES = {
    'pos': lambda k, rng: rng.uniform(0.2, 2.0, k),
    'pos_gap1e-6': _pos_pairs(1e-6),
    'pos_gap1e-12': _pos_pairs(1e-12),
    'pos_range': lambda k, rng: rng.permutation(np.logspace(-3, 1, k)),
    'pos_large': lambda k, rng: rng.uniform(0.2, 2.0, k) * 1e15,
    'pos_small': lambda k, rng: rng.uniform(0.2, 2.0, k) * 1e-15,
}


# ======================================================================
# A. algebra: spectra, functional calculus, norms
# ======================================================================
@pytest.mark.parametrize("complex_", [True, False])
@pytest.mark.parametrize("family", list(FAMILIES))
def test_eigenvalues_weyl(family, complex_):
    """eigenvalues / lambda_max / lambda_min: |d lambda| <= C k u ||A|| (Weyl)."""
    a32, a64 = algs(complex_valued=complex_)
    A32, A64 = pair(a32, a64, spectra_matrix(a64, FAMILIES[family], seed=1, complex_=complex_), is_self_adjoint=True)
    bound = KMAX * U * op_norm(A64)
    for w32, w64 in zip(A32.eigenvalues(), A64.eigenvalues()):
        if w32.numel():
            compare(w32, w64, bound, 2, f"eigenvalues[{family}]", norm=lambda t: t.abs().amax(-1))
    compare(A32.lambda_max, A64.lambda_max, bound, 2, "lambda_max", norm=lambda t: t.abs())
    compare(A32.lambda_min, A64.lambda_min, bound, 2, "lambda_min", norm=lambda t: t.abs())


@pytest.mark.parametrize("family", ['generic', 'gap1e-6', 'gap1e-12', 'range', 'large'])
def test_spectral_projectors_davis_kahan(family):
    """eigh: projectors onto eigenvalue clusters (clusters closer than 1e-4 ||A|| merged) satisfy
    ||dP||_F <= C k u ||A|| / gap (Davis-Kahan)."""
    a32, a64 = algs()
    A32, A64 = pair(a32, a64, spectra_matrix(a64, FAMILIES[family], seed=2), is_self_adjoint=True)
    nrm = op_norm(A64)
    for (w32, V32), (w64, V64) in zip(A32.eigh(), A64.eigh()):
        assert_single(w32, V32)
        k = w64.shape[-1]
        for b in range(w64.shape[0]):
            w = w64[b]
            cuts = [0] + [i + 1 for i in range(k - 1) if w[i + 1] - w[i] > 1e-4 * nrm[b]] + [k]
            for lo, hi in zip(cuts[:-1], cuts[1:]):
                if hi - lo == k:
                    continue
                gap = min(float(w[lo] - w[lo - 1]) if lo > 0 else math.inf,
                          float(w[hi] - w[hi - 1]) if hi < k else math.inf)
                P32 = d64(V32[b, :, lo:hi]) @ d64(V32[b, :, lo:hi]).mH
                P64 = V64[b, :, lo:hi] @ V64[b, :, lo:hi].mH
                check(torch.linalg.matrix_norm(P32 - P64), KMAX * U * nrm[b] / gap, 4, f"projector[{family}]")


SPECTRAL = {
    # name: (method, f, f', families)
    'exp': (lambda A: A.apply_function(torch.exp), torch.exp, torch.exp, ['generic', 'gap1e-6', 'gap1e-12', 'small']),
    'sin': (lambda A: A.apply_function(torch.sin), torch.sin, torch.cos, ['generic', 'gap1e-9', 'range']),
    'abs': (lambda A: A.abs(), torch.abs, torch.sign, ['generic', 'gap1e-9', 'range', 'large', 'small']),
    'inverse': (lambda A: A.inverse(), lambda w: 1 / w, lambda w: -1 / w ** 2, ['gap1e-6', 'gap1e-12', 'range', 'large',
                                                                                  'small']),
    'sqrt': (lambda A: A.sqrt(), torch.sqrt, lambda w: 0.5 / torch.sqrt(w), list(POS_FAMILIES)),
    'log': (lambda A: A.log(), torch.log, lambda w: 1 / w, ['pos', 'pos_gap1e-6', 'pos_gap1e-12', 'pos_range',
                                                            'pos_large', 'pos_small']),
    'power0.3': (lambda A: A.power(0.3), lambda w: w ** 0.3, lambda w: 0.3 * w ** -0.7, list(POS_FAMILIES)),
    'power2.5': (lambda A: A.power(2.5), lambda w: w ** 2.5, lambda w: 2.5 * w ** 1.5, ['pos', 'pos_gap1e-6',
                                                                                      'pos_range']),
    'power-0.5': (lambda A: A.power(-0.5), lambda w: w ** -0.5, lambda w: -0.5 * w ** -1.5, ['pos', 'pos_gap1e-12',
                                                                                           'pos_range', 'pos_large']),
}
SPECTRAL_CASES = [(n, f) for n, v in SPECTRAL.items() for f in v[3]]


@pytest.mark.parametrize("name,family", SPECTRAL_CASES)
def test_spectral_functions_daleckii_krein(name, family):
    """f(A) for self-adjoint A: ||df(A)||_F <= C k u (||A|| max|f[l_i, l_j]| + max|f|) (Daleckii-Krein), and for
    sqrt / A^0.3 also <= C (k u ||A||)^{1/2} (Hoelder, operator monotone), whichever is smaller."""
    method, f, df, _ = SPECTRAL[name]
    fams = {**FAMILIES, **POS_FAMILIES}
    a32, a64 = algs()
    flags = {'is_positive': True} if family.startswith('pos') else {}
    A32, A64 = pair(a32, a64, spectra_matrix(a64, fams[family], seed=3), is_self_adjoint=True, **flags)
    bound = dk_bound(A64, f, df)
    if name in ('sqrt', 'power0.3'):
        p = 0.5 if name == 'sqrt' else 0.3
        bound = torch.minimum(bound, math.sqrt(KMAX) * (KMAX * U * op_norm(A64)) ** p
                              + KMAX * U * op_norm(method(A64)))
    compare(method(A32), method(A64), bound, 4, f"{name}[{family}]")


@pytest.mark.parametrize("name", ['sqrt', 'power0.3', 'log', 'entropy'])
def test_functions_at_zero_eigenvalues(name):
    """Positive semidefinite A with eigenvalues 0 (rounded to +-u ||A|| in single precision): sqrt and A^p obey the
    Hoelder bound C (k u ||A||)^p; log (clamped at 1e-12) is finite; S(rho) obeys the Fannes-type bound
    sum_i eta(C k u), eta(t) = -t log t."""
    a32, a64 = algs()
    M = spectra_matrix(a64, FAMILIES['psd_zero'], seed=4)
    M = M / torch.diagonal(M, dim1=-2, dim2=-1).sum((-2, -1)).real[:, None, None, None]
    A32, A64 = pair(a32, a64, M, is_self_adjoint=True, is_positive=True)
    nrm = op_norm(A64)
    if name in ('sqrt', 'power0.3'):
        p = 0.5 if name == 'sqrt' else 0.3
        F = (lambda A: A.sqrt()) if name == 'sqrt' else (lambda A: A.power(0.3))
        compare(F(A32), F(A64), math.sqrt(KMAX) * (KMAX * U * nrm) ** p, 2, f"{name} at 0",
                norm=lambda t: torch.linalg.matrix_norm(t, ord=2).amax(-1))
    elif name == 'log':
        L32, L64 = A32.log(), A64.log()
        assert_single(L32)
        assert torch.isfinite(L32.matrix).all() and torch.isfinite(L64.matrix).all()
    else:
        b = entropy_bound(herm_eigs(A64), KMAX * U * nrm)
        for F in (lambda r: r.entropy(), st.von_neumann_entropy):
            compare(F(A32), F(A64), b, 2, "entropy at 0", norm=lambda x: x.abs())


@pytest.mark.parametrize("complex_", [True, False])
def test_general_operators(complex_):
    """Non-normal blocks: X^{-1} (||X^{-1}||^2 ||E||), |X| and the polar factor P (sqrt(2) ||E||_F,
    Araki-Yamagami), the polar unitary W (2 ||E|| / (s_min + s_min')), operator / trace / Frobenius norms
    (Lipschitz), expm (||E|| e^{||X||})."""
    a32, a64 = algs(complex_valued=complex_)
    X32, X64 = pair(a32, a64, general_matrix(a64, seed=5, shift=1.5, complex_=complex_))
    nrm = op_norm(X64)
    E = KMAX * U * nrm
    smin = torch.stack([torch.linalg.svdvals(X64.matrix[:, c, :k, :k])[:, -1] for c, k in enumerate(NS)],
                       -1).amin(-1)
    compare(X32.inverse(), X64.inverse(), E / smin ** 2, 2, "inverse(general)")
    compare(X32.abs(), X64.abs(), math.sqrt(2) * math.sqrt(KMAX) * E, 2, "abs(general)")
    W32, P32 = X32.polar()
    W64, P64 = X64.polar()
    compare(P32, P64, math.sqrt(2) * math.sqrt(KMAX) * E, 2, "polar P")
    compare(W32, W64, math.sqrt(KMAX) * E / smin, 2, "polar W")
    compare(X32.operator_norm(), X64.operator_norm(), E, 2, "operator_norm", norm=lambda t: t.abs())
    compare(X32.trace_norm(), X64.trace_norm(), sum(NS) * E, 2, "trace_norm", norm=lambda t: t.abs())
    compare(X32.frobenius_norm(), X64.frobenius_norm(), math.sqrt(sum(NS)) * E, 2, "frobenius", norm=lambda t: t.abs())
    compare(X32.expm(), X64.expm(), math.sqrt(KMAX) * E * torch.exp(nrm) + KMAX * U * op_norm(X64.expm()), 2, "expm")


@pytest.mark.parametrize("scale", [1e-30, 1e-15, 1e15, 1e30])
def test_norms_at_extreme_scales(scale):
    """Norms and |X| of operators of norm 1e-30 .. 1e30 (no overflow / underflow of intermediate squares in float32:
    1e30^2 > float32 max)."""
    a32, a64 = algs()
    X32, X64 = pair(a32, a64, general_matrix(a64, seed=6, scale=scale))
    H32, H64 = pair(a32, a64, spectra_matrix(a64, FAMILIES['generic'], seed=6) * scale, is_self_adjoint=True)
    E = KMAX * U * op_norm(X64)
    for name, F, C, b in [('frobenius', lambda A: A.frobenius_norm(), 2, math.sqrt(sum(NS)) * E),
                          ('operator_norm', lambda A: A.operator_norm(), 2, E),
                          ('trace_norm', lambda A: A.trace_norm(), 2, sum(NS) * E)]:
        compare(F(X32), F(X64), b, C, f"{name}[{scale:g}]", norm=lambda t: t.abs())
    compare(X32.abs(), X64.abs(), math.sqrt(2 * KMAX) * E, 2, f"abs[{scale:g}]")
    EH = KMAX * U * op_norm(H64)
    compare(H32.operator_norm(), H64.operator_norm(), EH, 2, f"op_norm(H)[{scale:g}]", norm=lambda t: t.abs())
    compare(H32.trace_norm(), H64.trace_norm(), sum(NS) * EH, 2, f"trace_norm(H)[{scale:g}]", norm=lambda t: t.abs())
    compare(H32.abs(), H64.abs(), math.sqrt(KMAX) * EH, 4, f"abs(H)[{scale:g}]")


def test_traces_contrast_and_arithmetic():
    """Traces (sums: n u sum|a_ii|), the Michelson contrast (from lambda_max, lambda_min), and +, -, @, scalar
    multiplication, adjoint, central elements (entrywise product bounds |dAB| <= k u |A||B|)."""
    a32, a64 = algs()
    A32, A64 = pair(a32, a64, spectra_matrix(a64, POS_FAMILIES['pos_range'], seed=7), is_self_adjoint=True,
                    is_positive=True)
    X32, X64 = pair(a32, a64, general_matrix(a64, seed=8))
    absdiag = torch.diagonal(X64.matrix, dim1=-2, dim2=-1).abs().sum((-2, -1))
    for name, F in [('trace', lambda A: A.trace), ('Tr_blunt', lambda A: A.Tr_blunt()),
                    ('Tr_norm', lambda A: A.Tr_norm()),
                    ('tau_vN', lambda A: A.tau_vN())]:
        compare(F(X32), F(X64), sum(NS) * U * absdiag, 2, name, norm=lambda t: t.abs())
    lmax, lmin = A64.lambda_max, A64.lambda_min
    cond = 2 * KMAX * U * op_norm(A64) * 2 * torch.maximum(lmax, lmin) / (lmax + lmin) ** 2
    compare(A32.michelson_contrast, A64.michelson_contrast, cond, 2, "contrast", norm=lambda t: t.abs())
    absX = X64.matrix.abs()
    prod_bound = fro(absX @ absX) * KMAX * U * 2
    compare(X32 @ X32, X64 @ X64, prod_bound, 2, "matmul")
    compare(X32 + X32.adjoint(), X64 + X64.adjoint(), fro(absX) * U * 2, 2, "add/adjoint")
    compare((X32 - A32) * 0.3, (X64 - A64) * 0.3, fro(absX + A64.matrix.abs()) * U * 3, 2, "sub/scalar")
    c32 = a32.central([1 / 3, 2.5, -1.25])
    c64 = a64.central([1 / 3, 2.5, -1.25])
    compare(c32, c64, U * 2.5 * math.sqrt(sum(NS)), 1, "central", norm=fro)
    compare(a32.identity(2), a64.identity(2), 1e-300, 0, "identity")


def test_unitary_constructors_single():
    """random_unitary (all measures), random_isometry, random_unitary_operator and operator_from_eigenvalues in
    single precision: ||U^*U - 1|| <= C k u, det = 1 for haar_su, U^T = U for COE, spectra reproduced to k u max|w|."""
    torch.manual_seed(9)
    a32, _ = algs()
    k = 6
    eye = torch.eye(k, dtype=torch.complex128)
    for m in a32.UNITARY_MEASURES:
        Q = a32.random_unitary(k, m, batch_size=5)
        assert_single(Q)
        check(torch.linalg.matrix_norm(d64(Q).mH @ d64(Q) - eye), k * U, 8, f"unitary[{m}]")
        if m == 'haar_su':
            check((torch.linalg.det(d64(Q)) - 1).abs(), k * U, 8, "det haar_su")
        if m == 'coe':
            check(torch.linalg.matrix_norm(d64(Q) - d64(Q).mT), k * U, 8, "coe symmetric")
    V = a32.random_isometry(7, 3, batch_size=4)
    assert_single(V)
    check(torch.linalg.matrix_norm(d64(V).mH @ d64(V) - torch.eye(3, dtype=torch.complex128)), 7 * U, 8, "isometry")
    Uop = a32.random_unitary_operator(4)
    Um = d64(Uop)
    check(fro(Um.mH @ Um - eye_of(a32).to(torch.complex128)), KMAX * U, 8, "unitary operator")
    w = torch.tensor([0.0, 1e-3, 0.5, 1e3])
    A = a32.operator_from_eigenvalues(lambda kk: w[:kk] + 1, batch_size=3)
    for kk, wa in zip(NS, A.eigenvalues()):
        torch.testing.assert_close(d64(wa), (w[:kk] + 1).double().expand(3, kk))     # known spectrum: exact
    Am = a32.like(precision='double').operator(d64(A), is_self_adjoint=True)
    for kk, wa in zip(NS, Am.eigenvalues()):
        check((wa - (w[:kk] + 1).double()).abs().amax(-1), KMAX * U * float(w[:kk].max() + 1), 4,
              "operator_from_eigenvalues")


def test_holomorphic_function_single():
    """holomorphic_function in single precision (default rtol 1e-5) agrees with expm, and with log on a custom contour,
    to C (rtol + k u e^{||X||}) relative."""
    a32, a64 = algs()
    X32, X64 = pair(a32, a64, general_matrix(a64, seed=10, scale=0.5))
    F32 = X32.holomorphic_function(torch.exp)
    E64 = X64.expm()
    compare(F32, E64, (1e-5 + KMAX * U * math.e) * fro(E64.matrix), 4, "holomorphic exp")
    P32, P64 = pair(a32, a64, spectra_matrix(a64, POS_FAMILIES['pos'], seed=10), is_self_adjoint=True,
                    is_positive=True)
    G32 = P32.holomorphic_function(torch.log, center=1.1, radius=1.05)
    compare(G32, P64.log(), (1e-5 + KMAX * U * 5) * fro(P64.log().matrix) + KMAX * U * 5, 8, "holomorphic log")


def test_normal_operator_functions_single():
    """apply_function of a unitary (normal path: eigh of H + alpha K in single precision, refined or redone in
    double when it does not separate), nearly degenerate phases: f(U) = U^2 + U^{-1} against the exact double
    result for U = exp(-iH), to C k u (1 + ||H||) max|f'| (the rounding of U included)."""
    a32, a64 = algs()
    H32, H64 = pair(a32, a64, spectra_matrix(a64, FAMILIES['gap1e-6'], seed=11), is_self_adjoint=True)
    U32 = a32.operator(dyn.propagator(H32, 1.0).matrix, is_normal=True)
    U64 = a64.operator(dyn.propagator(H64, 1.0).matrix, is_normal=True)
    F = lambda A: A.apply_function(lambda z: z ** 2 + 1 / z)
    compare(F(U32), F(U64), 3 * KMAX * U * math.sqrt(KMAX) * (1 + op_norm(H64)), 4, "f(unitary)")


# ======================================================================
# B. functional (levels 1-4 directly, including the GPU paths: dilation for SVD)
# ======================================================================
@pytest.mark.parametrize("method", ['svd', 'eigh'])
def test_singular_value_paths(method):
    """functional.singular_values / svd / abs_ / pinv / polar with method 'svd' and 'eigh' (the dilation used on GPUs):
    singular values to C k u ||X|| (absolute), |X| to C sqrt(2) k u ||X||, X^+ to C k u ||X|| / s_min^2, the polar
    unitary to C k u ||X|| / s_min."""
    a32, a64 = algs([5], )
    X64 = general_matrix(a64, B=4, seed=12, shift=1.0)[:, 0]
    X32 = X64.to(torch.complex64)
    Xr = X32.to(torch.complex128)
    s64 = torch.linalg.svdvals(Xr)
    nrm, smin = s64[:, 0], s64[:, -1]
    E = 5 * U * nrm
    compare(fc.singular_values(X32, method), s64, E, 2, f"singular_values[{method}]",
            norm=lambda t: t.abs().amax(-1))
    Uu, s, V = fc.svd(X32, method)
    assert_single(Uu, s, V)
    rec = d64(Uu) @ torch.diag_embed(d64(s)).to(torch.complex128) @ d64(V).mH
    check(fro(rec - Xr), E, 4, f"svd reconstruction[{method}]")
    compare(fc.abs_(X32, method), fc.abs_(Xr, method), math.sqrt(10) * E, 2, f"abs_[{method}]")
    compare(fc.pinv(X32, None, method), fc.pinv(Xr, None, method), math.sqrt(5) * E / smin ** 2, 2,
            f"pinv[{method}]")
    W32, P32 = fc.polar(X32, method)
    W64, P64 = fc.polar(Xr, method)
    compare(W32, W64, math.sqrt(5) * E / smin, 2, f"polar W[{method}]")
    compare(P32, P64, math.sqrt(10) * E, 2, f"polar P[{method}]")


def test_normal_eig_and_spectral_apply_single():
    """functional.normal_eig on unitary stacks with nearly tied Re + alpha Im values, in single precision: residual
    ||V^* X V - diag(lambda)|| <= C k u and unitary V; spectral_apply of both kinds."""
    torch.manual_seed(13)
    a32, _ = algs([6])
    Q = a32.random_unitary(6, 'haar', batch_size=8)
    phases = torch.exp(1j * torch.tensor([0.1, 0.1 + 1e-6, 1.0, 2.0, -1.0, 3.0])).to(torch.complex64)
    X = (Q * phases) @ Q.mH
    lam, V = fc.normal_eig(X)
    assert_single(lam, V)
    D = d64(V).mH @ d64(X) @ d64(V) - torch.diag_embed(d64(lam))
    check(fro(D), 6 * U, 8, "normal_eig residual")
    check(fro(d64(V).mH @ d64(V) - torch.eye(6, dtype=torch.complex128)), 6 * U, 8, "normal_eig unitary")
    Y = fc.spectral_apply(X, lambda z: z ** 3, 'normal')
    assert_single(Y)
    check(fro(d64(Y) - torch.linalg.matrix_power(d64(X), 3)), 3 * 6 * U * math.sqrt(6), 8, "spectral_apply normal")
    H = 0.5 * (X + X.mH)
    Z = fc.spectral_apply(H, torch.exp, 'hermitian')
    assert_single(Z)
    check(fro(d64(Z) - torch.linalg.matrix_exp(d64(H))), 6 * U * math.e * math.sqrt(6), 8, "spectral_apply herm")


def test_divided_differences_single():
    """functional.divided_differences in single precision: the switch between the quotient and the midpoint value of
    f' at eps^{1/3} balances cancellation (u |f| / |a - b|) and truncation (|a - b|^2 |f'''|), so the Loewner matrix
    of exp is accurate to C u^{2/3} max|f| (relative), for every gap; and values_and_derivative gives f and f'."""
    w64 = torch.tensor([1.0, 1 + 1e-1, 1 + 1e-3, 1 + 3e-5, 1 + 1e-5, 1 + 1e-6, 1 + 1e-7, 0.5, -1.0],
                       dtype=torch.float64)
    w32 = w64.float()
    wr = w32.double()                                  # the representable inputs
    f32, df32 = fc.values_and_derivative(torch.exp, w32)
    assert_single(f32, df32)
    check((d64(df32) - torch.exp(wr)).abs(), U * torch.exp(wr), 4, "values_and_derivative")
    L32 = fc.divided_differences(w32, f32, df32)
    assert_single(L32)
    a, b = wr.unsqueeze(-1), wr.unsqueeze(-2)
    ref = torch.where(a == b, torch.exp(a).expand(-1, len(wr)),
                      (torch.exp(a) - torch.exp(b)) / torch.where(a == b, torch.ones_like(a), a - b))
    check((d64(L32) - ref).abs(), U ** (2 / 3) * math.e, 2, "divided differences")


def test_default_contour_and_holomorphic_level4():
    """functional.default_contour / holomorphic on single-precision stacks: the circle encloses the spectrum, f(X)
    to the default single-precision tolerance."""
    X = (torch.randn(6, 4, 4, dtype=torch.complex64, generator=torch.Generator().manual_seed(14)) * 0.5)
    c, r = fc.default_contour(X)
    assert_single(c, r)
    lam = torch.linalg.eigvals(X.to(torch.complex128))
    assert ((lam - d64(c)[:, None]).abs() < d64(r)[:, None]).all()
    F, err = fc.holomorphic(X, torch.exp)
    assert_single(F)
    ref = torch.linalg.matrix_exp(X.to(torch.complex128))
    check(fro(d64(F) - ref), (1e-5 + 4 * U * math.e) * fro(ref), 4, "holomorphic level 4")


# ======================================================================
# C. states
# ======================================================================
def grad_bound(F, mats, flags=None):
    """First-order bound u (k_max sum_inputs max_c ||A_c||_op sum_c ||grad_{A_c} F||_1 + |F|) per sample for a
    scalar function F of double tensors (gradient by autograd in double, tested in test_gradients.py)."""
    mats = [m.detach().clone().requires_grad_(True) for m in mats]
    val = F(*mats)
    grads = torch.autograd.grad(val.sum(), mats, allow_unused=True)
    out = val.detach().abs().clone()
    for m, g in zip(mats, grads):
        if g is None:
            continue
        nrm = torch.linalg.matrix_norm(m.detach(), ord=2).amax(-1)
        gn = torch.linalg.svdvals(g).sum((-2, -1))
        out = out + KMAX * nrm * gn
    return U * out


def entropy_bound(w, t):
    """First-order / Fannes bound for S = -sum eta(lambda_i) under eigenvalue errors <= t: sum |log lambda + 1| t over
    the eigenvalues above e t, eta(t) for the others (|eta(x + d) - eta(x)| <= eta(|d|) for |d| <= 1/2)."""
    t = torch.as_tensor(t, dtype=torch.float64).reshape(-1, 1)
    big = w > math.e * t
    lin = torch.where(big, (torch.log(w.clamp(min=1e-300)) + 1).abs() * t, -t * torch.log(t))
    return lin.sum(-1)


def states_pair(a32, a64, B=4, seed=0, rank=None, lam_min=0.0):
    """A random state (HS measure, optionally of given rank, optionally mixed with lam_min * 1) in both precisions."""
    torch.manual_seed(seed)
    r = st.random_density_matrix(a64, B, rank=rank).matrix
    if lam_min:
        r = (1 - lam_min * sum(NS)) * r + lam_min * eye_of(a64).to(r.dtype)
    return pair_state(a32, a64, r)


def H_pos(*ms):
    return [vn.TypeIAlgebra(NS, NS, precision='double').operator(m, is_self_adjoint=True, is_positive=True) for m in ms]


def test_random_states_single():
    """random_density_matrix (HS, induced rank, Bures), maximally_mixed_state, tracial_state in single precision:
    Tr rho = 1 to C N u, Hermitian, lambda_min >= -C k u."""
    a32, _ = algs()
    torch.manual_seed(15)
    for rho in (st.random_density_matrix(a32, 6), st.random_density_matrix(a32, 6, rank=2),
                st.random_density_matrix(a32, 6, measure='bures'), st.maximally_mixed_state(a32, 2),
                st.tracial_state(a32, 2)):
        assert_single(rho)
        m = d64(rho)
        check((torch.diagonal(m, dim1=-2, dim2=-1).sum((-2, -1)) - 1).abs(), sum(NS) * U, 4, "trace one")
        check(fro(m - m.mH), U, 1, "hermitian")
        check((-herm_eigs(a32.like(precision='double').operator(m, is_self_adjoint=True))).clamp(min=0).amax(-1),
              KMAX * U, 4, "positive")


@pytest.mark.parametrize("beta", [0.5, 5.0, 50.0])
def test_gibbs_and_partition_function(beta):
    """gibbs_state: ||d rho||_F <= C 2 beta ||E||_F (the Daleckii-Krein constant of exp(-beta x)/Z is <= beta);
    log Z: |d log Z| <= beta k u ||H||; Z itself may overflow float32 (only finite where representable)."""
    a32, a64 = algs()
    H32, H64 = pair(a32, a64, spectra_matrix(a64, FAMILIES['generic'], seed=16) * 2, is_self_adjoint=True)
    E = KMAX * U * op_norm(H64)
    compare(st.gibbs_state(H32, beta), st.gibbs_state(H64, beta), 2 * beta * math.sqrt(KMAX) * E + KMAX * U, 2,
            f"gibbs[{beta}]")
    compare(st.log_partition_function(H32, beta), st.log_partition_function(H64, beta), beta * E + U, 2,
            f"log Z[{beta}]", norm=lambda t: t.abs())
    Z32, Z64 = st.partition_function(H32, beta), st.partition_function(H64, beta)
    assert_single(Z32)
    assert_finite_where(Z32, Z64)
    ok = Z64 < 1e-3 * F32MAX
    check(((d64(Z32) - Z64).abs() / Z64)[ok], (beta * E + 2 * U)[ok], 4, f"Z[{beta}]")


@pytest.mark.parametrize("kind", ['generic', 'near_singular'])
def test_state_functionals(kind):
    """von_neumann_entropy, relative_entropy, fidelity, trace_distance, purity, expectation, Born probabilities,
    sector probabilities, Lueders rule, mix: errors within C times the first-order bound (gradient in double);
    near singular states (rank 2): the entropy obeys the Fannes-type bound, the fidelity the Hoelder bound
    4 k sqrt(k u) (sqrt is not Lipschitz at 0)."""
    a32, a64 = algs()
    if kind == 'generic':
        r32, r64 = states_pair(a32, a64, seed=17, lam_min=0.02)
    else:
        r32, r64 = states_pair(a32, a64, seed=17, rank=2)
    s32, s64 = states_pair(a32, a64, seed=18, lam_min=0.02)
    op = lambda m: a64.operator(m, is_self_adjoint=True, is_positive=True)
    t = KMAX * U
    b_rel = grad_bound(lambda r, s: st.relative_entropy(op(r), op(s)), [r64.matrix, s64.matrix])
    if kind == 'generic':
        b_ent = grad_bound(lambda r: st.von_neumann_entropy(op(r)), [r64.matrix])
        b_fid = grad_bound(lambda r, s: st.fidelity(op(r), op(s)), [r64.matrix, s64.matrix])
    else:
        # entropy: Fannes-type; fidelity: F = ||rho^{1/2} sigma^{1/2}||_1^2 and ||d rho^{1/2}|| <= ||d rho||^{1/2} on
        # each of the n_0 kernel directions of rho, so |dF| <= 2 F^{1/2} n_0 (k u ||rho||)^{1/2} ||sigma||^{1/2}
        w = herm_eigs(r64)
        b_ent = entropy_bound(w, t)
        b_rel = b_rel + b_ent
        n0 = (w < 1e-10).sum(-1)
        F = st.fidelity(r64, s64)
        b_fid = 2 * F.sqrt() * n0 * (t * op_norm(r64)).sqrt() * op_norm(s64).sqrt() + grad_bound(
            lambda r, s: st.fidelity(op(r), op(s)), [s64.matrix, s64.matrix])
    compare(st.von_neumann_entropy(r32), st.von_neumann_entropy(r64), b_ent, 4, f"entropy[{kind}]", norm=torch.abs)
    compare(st.fidelity(r32, s32), st.fidelity(r64, s64), b_fid, 4, f"fidelity[{kind}]", norm=torch.abs)
    compare(st.relative_entropy(r32, s32), st.relative_entropy(r64, s64), b_rel, 4, f"relative_entropy[{kind}]",
            norm=torch.abs)
    compare(st.trace_distance(r32, s32), st.trace_distance(r64, s64), sum(NS) * t, 2, "trace_distance", norm=torch.abs)
    compare(st.purity(r32), st.purity(r64), sum(NS) * U * 2, 2, "purity", norm=torch.abs)
    A32, A64 = pair(a32, a64, spectra_matrix(a64, FAMILIES['generic'], B=4, seed=19), is_self_adjoint=True)
    compare(st.expectation(r32, A32), st.expectation(r64, A64), sum(NS) * U * 2 * op_norm(A64), 2, "expectation",
            norm=torch.abs)
    compare(r32.sector_probabilities(), r64.sector_probabilities(), KMAX * U, 2, "sector_probabilities",
            norm=lambda x: x.abs().amax(-1))
    P = torch.zeros(1, 3, KMAX, KMAX, dtype=torch.complex128)
    P[0, 0, 0, 0] = P[0, 0, 1, 1] = P[0, 1, 0, 0] = 1
    P32, P64 = pair(a32, a64, P, is_self_adjoint=True, is_projection=True)
    compare(st.born_probabilities(r32, [P32, a32.identity() - P32]),
            st.born_probabilities(r64, [P64, a64.identity() - P64]),
            sum(NS) * U * 2, 2, "born", norm=lambda x: x.abs().amax(-1))
    c32, p32 = r32.condition_on(P32)
    c64, p64 = r64.condition_on(P64)
    compare(p32, p64, sum(NS) * U * 2, 2, "lueders probability", norm=torch.abs)
    compare(c32, c64, KMAX * U * 4 / p64, 2, "lueders state")
    compare(r32.mix(s32, 0.3), r64.mix(s64, 0.3), 2 * U, 2, "mix", norm=lambda x: x.abs().amax((-3, -2, -1)))
    for tr in ('norm', 'tau_vN'):
        d32, d64_ = r32.density(tr), r64.density(tr)
        compare(d32, d64_, 3 * KMAX * U, 2, f"density[{tr}]", norm=lambda x: x.abs().amax((-3, -2, -1)))
        back = DensityMatrix.from_density(a32, d32, tr)
        compare(back, r64, 2 * 3 * KMAX * U, 2, f"from_density[{tr}]", norm=lambda x: x.abs().amax((-3, -2, -1)))


# ======================================================================
# D. channels
# ======================================================================
def _abs_apply(K, rho):
    """sum_i |K_i| |rho| |K_i|^T (entrywise absolute values): the Higham bound for K rho K^*."""
    Ka, ra = K.abs(), rho.abs()
    return (Ka @ ra.unsqueeze(1) @ Ka.transpose(-2, -1)).sum(1)


def test_channels_single():
    """Channel constructors and maps in single precision: trace preservation / unitality to C r k u, apply (Higham
    bound (2k + r) u sum |K||rho||K|^T entrywise), adjoint, compose, mix, choi, superoperator, from_superoperator,
    InterSectorChannel (apply, adjoint, compose, transition_matrix, TP / unital checks)."""
    torch.manual_seed(20)
    a32, a64 = algs()
    r32, r64 = states_pair(a32, a64, B=2, seed=21)
    chans = {'random': ch.random_channel(a32, 3, 2), 'dephasing': ch.dephasing_channel(a32, 0.3),
             'depolarizing': ch.depolarizing_channel(a32, 0.4), 'damping': ch.amplitude_damping_channel(a32, 0.2),
             'center': ch.center_expectation(a32), 'identity': ch.identity_channel(a32),
             'mixed_unitary': ch.random_mixed_unitary_channel(a32, 3, 2),
             'unitary': ch.unitary_channel(a32.random_unitary_operator(2))}
    for name, c in chans.items():
        assert_single(c.kraus)
        r = c.kraus_rank
        check(c.trace_preservation_error(), (r + KMAX) * U, 4, f"TP[{name}]")
        assert c.is_trace_preserving()
        c64 = ch.Channel(a64, d64(c.kraus))
        bound = fro(_abs_apply(d64(c.kraus), r64.matrix)) * (2 * KMAX + r) * U
        compare(c.apply(r32), c64.apply(r64), bound, 2, f"apply[{name}]")
        compare(c.adjoint().apply(a32.identity()), c64.adjoint().apply(a64.identity()),
                fro(_abs_apply(d64(c.kraus).conj().transpose(-2, -1), a64.identity().matrix)) * (2 * KMAX + r) * U, 2,
                f"adjoint[{name}]")
    for name in ('dephasing', 'depolarizing', 'mixed_unitary', 'center'):
        check(chans[name].unitality_error(), (chans[name].kraus_rank + KMAX) * U, 4, f"unital[{name}]")
    c1, c2 = chans['random'], chans['damping']
    comp = c1.compose(c2).mix(chans['depolarizing'], 0.25)
    assert_single(comp.kraus)
    check(comp.trace_preservation_error(), comp.kraus_rank * KMAX * U, 4, "TP compose/mix")
    for J32, J64 in zip(c1.choi(), ch.Channel(a64, d64(c1.kraus)).choi()):
        compare(J32, J64, 3 * KMAX * U * fro(J64.abs()), 2, "choi")
    S = c1.superoperator()
    assert_single(*S)
    rebuilt = ch.Channel.from_superoperator(a32, S)
    assert_single(rebuilt.kraus)
    check(rebuilt.trace_preservation_error(), 4 * KMAX ** 2 * U, 4, "from_superoperator TP")
    compare(rebuilt.apply(r32), ch.Channel(a64, d64(c1.kraus)).apply(r64), 8 * KMAX ** 2 * U, 2,
            "from_superoperator apply")
    b32 = TypeIAlgebra([2, 3], [2, 3], precision='single')
    I = ch.random_inter_sector_channel(a32, b32, 2, 2)
    assert_single(I.kraus)
    check(I.trace_preservation_error(), 2 * 5 * U, 4, "inter TP")
    assert I.is_trace_preserving()
    I64 = vn.InterSectorChannel(a64, b32.like(precision='double'), d64(I.kraus))
    compare(I.apply(r32), I64.apply(r64), 4 * 5 * U, 2, "inter apply")
    compare(I.transition_matrix(r32), I64.transition_matrix(r64), 4 * 5 * U, 2, "transition_matrix",
            norm=lambda x: x.abs().amax((-2, -1)))
    compare(I.adjoint('tau_vN').apply(b32.identity()),
            I64.adjoint('tau_vN').apply(b32.like(precision='double').identity()),
            4 * 5 * U * 3, 2, "inter adjoint")
    J = I.compose(chans['random'])
    check(J.trace_preservation_error(), 2 * 3 * 5 * U, 4, "inter compose TP")
    assert vn.InterSectorChannel(a32, a32, c1.to_inter_sector().kraus).is_trace_preserving()


# ======================================================================
# E. dynamics
# ======================================================================
@pytest.mark.parametrize("t", [0.7, 300.0])
def test_unitary_dynamics_single(t):
    """propagator, schrodinger, von_neumann: the phases exp(-i t lambda) carry the absolute error t |d lambda|, so the
    error is C k u (1 + t ||H||) (the first-order bound; the phases themselves are formed in double)."""
    a32, a64 = algs()
    H32, H64 = pair(a32, a64, spectra_matrix(a64, FAMILIES['gap1e-9'], seed=22), is_self_adjoint=True)
    bound = KMAX * U * (1 + t * op_norm(H64)) * math.sqrt(KMAX) * 2
    compare(dyn.propagator(H32, t), dyn.propagator(H64, t), bound, 2, f"propagator[{t}]")
    psi = torch.randn(3, 3, KMAX, dtype=torch.complex128, generator=torch.Generator().manual_seed(23)) \
        * rows_of(a64)
    psi = psi / torch.linalg.vector_norm(psi, dim=(-2, -1), keepdim=True)
    p32 = dyn.schrodinger(psi.to(torch.complex64), H32, [0.0, t])
    p64 = dyn.schrodinger(psi.to(torch.complex64).to(torch.complex128), H64, [0.0, t])
    compare(p32[-1], p64[-1], bound, 2, f"schrodinger[{t}]")
    r32, r64 = states_pair(a32, a64, B=3, seed=24)
    compare(dyn.von_neumann(r32, H32, [0.0, t])[-1], dyn.von_neumann(r64, H64, [0.0, t])[-1], 2 * bound, 2,
            f"von_neumann[{t}]")


def test_ode_solvers_single():
    """rk4 / schrodinger_rk4 / solve_operator_ode / lindblad_evolve (rk4, kraus, exact) / lindblad_rhs in single
    precision against double with the same steps: rounding accumulates to C n_steps k u (step norms <= 1)."""
    a32, a64 = algs()
    H32, H64 = pair(a32, a64, spectra_matrix(a64, FAMILIES['generic'], seed=25), is_self_adjoint=True)
    L32, L64 = pair(a32, a64, general_matrix(a64, seed=26, scale=0.5))
    r32, r64 = states_pair(a32, a64, B=3, seed=27)
    steps = 40
    nH = op_norm(H64)
    bound = 2 * steps * U + KMAX * U * (1 + nH)          # n roundings of y + h k (||y|| <= 1), plus t ||dH||
    psi = (torch.randn(3, 3, KMAX, dtype=torch.complex128, generator=torch.Generator().manual_seed(28))
           * rows_of(a64))
    psi = (psi / torch.linalg.vector_norm(psi, dim=(-2, -1), keepdim=True)).to(torch.complex64)
    compare(dyn.schrodinger_rk4(psi, H32, [0.0, 1.0], substeps=steps)[-1],
            dyn.schrodinger_rk4(psi.to(torch.complex128), H64, [0.0, 1.0], substeps=steps)[-1], bound, 2,
            "schrodinger_rk4")
    f = lambda t, X: -1j * (X @ X - X)
    compare(dyn.solve_operator_ode(lambda t, X: 0.1 * f(t, X), r32, [0.0, 1.0], substeps=steps)[-1],
            dyn.solve_operator_ode(lambda t, X: 0.1 * f(t, X), r64, [0.0, 1.0], substeps=steps)[-1], bound, 2,
            "solve_operator_ode")
    for method in ('rk4', 'kraus', 'exact'):
        o32 = dyn.lindblad_evolve(r32, H32, [L32], times=[0.0, 1.0], substeps=steps, method=method)[-1]
        o64 = dyn.lindblad_evolve(r64, H64, [L64], times=[0.0, 1.0], substeps=steps, method=method)[-1]
        compare(o32, o64, bound + KMAX * U * op_norm(L64) ** 2, 2, f"lindblad_evolve[{method}]")
    rhs32 = dyn.lindblad_rhs(H32, [L32])(0.0, r32.matrix)
    rhs64 = dyn.lindblad_rhs(H64, [L64])(0.0, r64.matrix)
    compare(rhs32, rhs64, 4 * KMAX * U * (nH + op_norm(L64) ** 2), 2, "lindblad_rhs")


def test_lindblad_channel_and_superoperator_single():
    """lindblad_channel of a single-precision algebra (computed in double, Kraus operators returned in single):
    trace preserving to C k^2 u, action within C k^2 u of the double channel; lindblad_superoperator values."""
    a32, a64 = algs()
    H32, H64 = pair(a32, a64, spectra_matrix(a64, FAMILIES['generic'], B=1, seed=29), is_self_adjoint=True)
    L32, L64 = pair(a32, a64, general_matrix(a64, B=1, seed=30, scale=0.5))
    c32 = dyn.lindblad_channel(a32, H32, [L32], t=0.6)
    c64 = dyn.lindblad_channel(a64, H64, [L64], t=0.6)
    assert_single(c32.kraus)
    check(c32.trace_preservation_error(), KMAX ** 2 * U, 4, "lindblad_channel TP")
    r32, r64 = states_pair(a32, a64, B=2, seed=31)
    compare(c32.apply(r32), c64.apply(r64), (2 * KMAX + KMAX ** 2) * U * fro(_abs_apply(d64(c32.kraus), r64.matrix))
            + 2 * KMAX ** 2 * U, 2, "lindblad_channel apply")
    for S32, S64 in zip(dyn.lindblad_superoperator(a32, H32, [L32]), dyn.lindblad_superoperator(a64, H64, [L64])):
        check(fro(d64(S32) - S64), 4 * KMAX * U * (op_norm(H64) + op_norm(L64) ** 2) * KMAX, 2, "superoperator")


# ======================================================================
# F. composite systems
# ======================================================================
def test_composite_single():
    """kron (entrywise products: 2 u |a||b|), partial_trace and partial_trace_channel (sums: m u sum|x|), in plain
    and fused products of single-precision algebras."""
    torch.manual_seed(32)
    b32, b64 = algs([1, 2])
    c32, c64 = algs([2, 1])
    for prod in (composite.tensor_product, composite.fused_tensor_product):
        p32, p64 = prod(b32, c32), prod(b64, c64)
        assert p32.dtype == torch.complex64
        x32, x64 = states_pair(b32, b64, B=2, seed=33)
        y32, y64 = states_pair(c32, c64, B=2, seed=34)
        k32 = composite.kron(x32, y32, p32)
        k64 = composite.kron(x64, y64, p64)
        compare(k32, k64, 3 * U, 2, "kron", norm=lambda x: x.abs().amax((-3, -2, -1)))
        for keep in (1, 2):
            pt32 = composite.partial_trace(k32, keep)
            pt64 = composite.partial_trace(k64, keep)
            compare(pt32, pt64, 4 * 3 * U, 2, f"partial_trace[{keep}]", norm=lambda x: x.abs().amax((-3, -2, -1)))
            ch32 = composite.partial_trace_channel(p32, keep)
            compare(ch32.apply(k32), pt64, 4 * 3 * U, 2, f"partial_trace_channel[{keep}]",
                    norm=lambda x: x.abs().amax((-3, -2, -1)))
        s, off = composite.pair_offset(p32, 1, 0)
        assert s == composite.sector_index(p32, 1, 0)


# ======================================================================
# G. chains, H. krylov, L. su2
# ======================================================================
def test_spin_chain_single():
    """SpinChain(precision='single'): XXZ spectra (full algebra, sectors, momentum blocks) to C k u ||H|| (Weyl),
    level spacing ratio to C k u ||H|| / min spacing, reduced states / entanglement entropy of eigenstates (Fannes),
    sz / number / basis states exact, the damping and lowering maps."""
    L = 6
    c32, c64 = SpinChain(L, 'periodic', precision='single'), SpinChain(L, 'periodic', precision='double')
    h = torch.tensor([[0.3, -0.2, 0.5, 0.1, -0.4, 0.25]], dtype=torch.float64)
    H32, H64 = c32.xxz(1.0, 0.6, h.float()), c64.xxz(1.0, 0.6, h.float().double())
    assert_single(H32)
    nrm = op_norm(H64)
    k = max(c64.algebra.k_factors)
    for w32, w64 in zip(H32.eigenvalues(), H64.eigenvalues()):
        compare(w32, w64, w64.shape[-1] * U * w64.abs().amax(-1), 2, "xxz spectrum", norm=lambda x: x.abs().amax(-1))
    E32, E64 = torch.cat(H32.eigenvalues(), -1), torch.cat(H64.eigenvalues(), -1)
    s64 = torch.diff(torch.sort(E64, -1)[0], dim=-1)
    big = torch.maximum(s64[..., :-1], s64[..., 1:])                 # d r_n <= (|d s_n| + |d s_n+1|) / max(s_n, s_n+1)
    lr = 4 * k * U * nrm[:, None] / big
    compare(SpinChain.level_spacing_ratio(E32), SpinChain.level_spacing_ratio(E64), lr.mean(-1), 2,
            "level spacing ratio", norm=torch.abs)
    for N in range(L + 1):
        M32, M64 = c32.xxz_momentum(1.0, 0.6, 0.2, sector=N), c64.xxz_momentum(1.0, 0.6, 0.2, sector=N)
        w32 = torch.sort(torch.cat([w for w in M32.eigenvalues() if w.numel()], -1))[0]
        w64 = torch.sort(torch.cat([w for w in M64.eigenvalues() if w.numel()], -1))[0]
        compare(w32, w64, k * U * op_norm(M64), 2, "momentum spectrum", norm=lambda x: x.abs().amax(-1))
    for i in range(L):
        torch.testing.assert_close(d64(c32.sz(i)), d64(c64.sz(i)), rtol=0, atol=0)
    torch.testing.assert_close(d64(c32.number()), d64(c64.number()), rtol=0, atol=0)
    torch.testing.assert_close(d64(c32.basis_state('010110')), d64(c64.basis_state('010110')), rtol=0, atol=0)
    # ground state of the half-filled sector as a state, reduced to 3 sites
    w, V = H64.eigh()[3]
    psi = V[:, :, 0]
    rho = torch.zeros(1, L + 1, c64.algebra.k_max, c64.algebra.k_max, dtype=torch.complex128)
    rho[0, 3, :20, :20] = psi[0][:, None] * psi[0].conj()[None, :]
    r32, r64 = pair_state(c32.algebra, c64.algebra, rho)
    red64 = c64.reduced_state(r64, 3)
    compare(c32.reduced_state(r32, 3), red64, 20 * U, 2, "reduced_state", norm=lambda x: x.abs().amax((-3, -2, -1)))
    compare(c32.entanglement_entropy(r32, 3), c64.entanglement_entropy(r64, 3), entropy_bound(herm_eigs(red64), 4 * U),
            2, "entanglement entropy", norm=torch.abs)
    d32 = c32.site_amplitude_damping(2, 0.3)
    assert_single(d32.kraus)
    check(d32.trace_preservation_error(), 4 * U, 4, "site damping TP")
    assert_single(c32.lowering(1).kraus)
    rf = c32.random_field_heisenberg(1.0, 2, generator=torch.Generator().manual_seed(0))
    assert_single(rf)


def test_krylov_single():
    """krylov in single precision: SparseSectorHamiltonian.matvec / to_dense (row sums), lanczos (Ritz values to
    C m k u ||H||), ground_state (E_0 to C (k u ||H|| + tol)), expm_multiply_step and evolve (state error
    C n_steps m u ||H|| dt), without warnings for the default tolerances."""
    L = 10
    c32, c64 = SpinChain(L, precision='single'), SpinChain(L, precision='double')
    h = torch.linspace(-0.5, 0.5, L)
    S32, S64 = c32.xxz_sparse(1.0, 0.7, h, sector=5), c64.xxz_sparse(1.0, 0.7, h.double(), sector=5)
    D32, D64 = S32.to_dense(), S64.to_dense()
    assert_single(D32)
    nrm = torch.linalg.matrix_norm(D64, ord=2)
    check(fro(d64(D32) - D64), U * fro(D64), 1, "to_dense")
    v = torch.randn(2, S64.dim, dtype=torch.float64, generator=torch.Generator().manual_seed(35))
    compare(S32.matvec(v.float()), S64.matvec(v.float().double()), 8 * U * nrm * fro(v), 2, "matvec")
    a32, b32, Q32 = krylov.lanczos(S32.matvec, v.float(), 20)
    a64, b64, Q64 = krylov.lanczos(S64.matvec, v.float().double(), 20)
    assert_single(a32, b32, Q32)
    T32 = torch.linalg.eigvalsh(torch.diag_embed(d64(a32)) + torch.diag_embed(d64(b32), 1)
                                + torch.diag_embed(d64(b32), -1))
    T64 = torch.linalg.eigvalsh(torch.diag_embed(a64) + torch.diag_embed(b64, 1) + torch.diag_embed(b64, -1))
    check((T32[:, 0] - T64[:, 0]).abs(), 8 * U * nrm, 4, "lanczos lowest Ritz value")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        E32, p32, res32 = krylov.ground_state(S32, generator=torch.Generator().manual_seed(0))
    E64, p64, _ = krylov.ground_state(S64, generator=torch.Generator().manual_seed(0))
    assert_single(E32, p32)
    check((d64(E32) - E64).abs(), 8 * U * nrm + d64(res32), 4, "ground_state")
    psi0 = torch.zeros(1, S64.dim, dtype=torch.complex128)
    psi0[0, 0] = 1
    o32, e32 = krylov.expm_multiply_step(S32, psi0.to(torch.complex64), 0.3, m=20)
    o64, _ = krylov.expm_multiply_step(S64, psi0, 0.3, m=20)
    assert_single(o32)
    check(fro(d64(o32) - o64), 20 * U * (1 + 0.3 * nrm), 4, "expm_multiply_step")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ev32 = krylov.evolve(S32, psi0.to(torch.complex64), [0.0, 1.0, 2.0])
    ev64 = krylov.evolve(S64, psi0, [0.0, 1.0, 2.0])
    assert_single(ev32)
    check(fro(d64(ev32[-1]) - ev64[-1]), 8 * U * (1 + 2 * nrm), 4, "evolve")      # t ||dH||, ||dH|| <= 8 u ||H||


def test_su2_chain_single():
    """SU2Chain(precision='single'): Heisenberg J1-J2 spectra (full_spectrum, lowest_energies, ground state energy) to
    C m u ||H||; casimir exact; restrict(lift(A)) = A; reduced matrix elements and spin correlations to C u."""
    L = 6
    s32, s64 = SU2Chain(L, precision='single'), SU2Chain(L, precision='double')
    H32, H64 = s32.heisenberg(1.0, 0.3, 0.4), s64.heisenberg(1.0, 0.3, 0.4)
    assert_single(H32)
    m = max(s64.multiplicities)
    nrm = op_norm(H64)
    compare(s32.full_spectrum(H32), s64.full_spectrum(H64), m * U * nrm + L * U, 2, "full_spectrum",
            norm=lambda x: x.abs().amax(-1))
    compare(s32.lowest_energies(H32), s64.lowest_energies(H64), m * U * nrm + L * U, 2, "lowest_energies",
            norm=lambda x: x.abs().amax(-1))
    g32, g64 = s32.ground_state(H32), s64.ground_state(H64)
    check((d64(g32.energy) - g64.energy).abs(), m * U * nrm + L * U, 2, "ground energy")
    torch.testing.assert_close(d64(s32.casimir()), d64(s64.casimir()), rtol=0, atol=0)
    B = s32.restrict(s32.lift(H32))
    compare(B, H64, 20 * m * U * nrm, 2, "restrict(lift)")
    R32, R64 = s32.reduced_matrix_elements(), s64.reduced_matrix_elements()
    for key in R64:
        check((d64(R32[key]) - R64[key]).abs().amax(), 20 * U * (1 + R64[key].abs().amax()), 2, "reduced elements")
    C32 = s32.spin_correlation([(0, 3)])
    C64 = s64.spin_correlation([(0, 3)])
    compare(C32, C64, 2 * m * U * C64.matrix.abs().amax(), 2, "spin_correlation",
            norm=lambda x: x.abs().amax((-3, -2, -1)))
    hw = s32.hw_basis(1, dtype=torch.float32)
    check(fro((hw.T.double() @ hw.double() - torch.eye(hw.shape[1], dtype=torch.float64))[None]), hw.shape[0] * U, 4,
          "hw basis orthonormal")


# ======================================================================
# I. decompose, K. jordan, M. inclusions, N. optimize
# ======================================================================
def test_wedderburn_single():
    """decompose in single precision (complex64 input): the same sectors and multiplicities as in double, verified
    residuals ~ n u (documented ~1e-6), commutant dimension, twirl to C |G| n u, round trips to C n u, central
    projections idempotent, real structure types as in double."""
    G32 = dc.total_spin(4, dtype=torch.complex64)
    D32 = dc.wedderburn(G32)
    D64 = dc.wedderburn(dc.total_spin(4))
    assert D32.diagnostics['precision'] == 'single'
    assert (D32.k_factors, D32.multiplicities) == (D64.k_factors, D64.multiplicities)
    n = 16
    res = D32.verify()
    assert res['ok']
    for key in ('unitarity', 'algebra', 'commutant'):
        check(res[key], n * U, 8, f"wedderburn {key}")
    A = D32.to_algebra(G32)
    assert_single(A)
    check(fro((d64(D32.from_algebra(A)) - G32.to(torch.complex128))), n * U * 4, 4, "to/from algebra")
    P = D32.central_projections()
    assert_single(P)
    check(fro(d64(P) @ d64(P) - d64(P)), n * U, 8, "central projections")
    perm32 = dc.permutation_representation(3, 2, dtype=torch.complex64)
    perm64 = dc.permutation_representation(3, 2)
    assert dc.commutant_dimension(perm32) == dc.commutant_dimension(perm64)
    X = torch.randn(8, 8, dtype=torch.complex128, generator=torch.Generator().manual_seed(36))
    Xr = X.to(torch.complex64)
    check(fro((d64(dc.twirl(perm32, Xr)) - dc.twirl(perm64, Xr.to(torch.complex128)))[None]),
          2 * 8 * U * fro(X[None]), 2, "twirl")
    B32 = dc.commutant_basis(perm32)
    B64 = dc.commutant_basis(perm64)
    assert B32.shape[0] == B64.shape[0]
    R32, R64 = D32.real_structure(), D64.real_structure()
    assert R32.ok and R32.types == R64.types and R32.partner == R64.partner
    S = dc.spin_matrices(1.0, dtype=torch.complex64)
    comm = d64(S[0]) @ d64(S[1]) - d64(S[1]) @ d64(S[0]) - 1j * d64(S[2])
    check(fro(comm[None]), 3 * U, 4, "spin matrices")


def test_jordan_single():
    """SpinFactor(float32): eigenvalues a +- |u| to C u (|a| + |u|), functional calculus (sqrt, exp, log, inverse,
    power) to C u max|f'| (|a| + |u|) + u |f|, determinant to C u (a^2 + |u|^2); the Jordan maps of the algebra
    (product, quadratic representation, triple, power, transpose) to the Higham bounds; jordan_algebra_generated
    in single has the dimension of the double one; anticommuting symmetries to C u."""
    V32, V64 = vn.SpinFactor(5, torch.float32), vn.SpinFactor(5)
    g = torch.Generator().manual_seed(37)
    u = torch.randn(64, 5, dtype=torch.float64, generator=g)
    a = u.norm(dim=-1) * (1 + torch.rand(64, dtype=torch.float64, generator=g))
    x64 = V64.element(a, u)
    x32 = x64.float()
    xr = x32.double()
    scale = xr[:, 0].abs() + xr[:, 1:].norm(dim=-1)
    lam = V64.eigenvalues(xr)
    compare(V32.eigenvalues(x32), lam, 4 * U * scale, 2, "spin factor eigenvalues", norm=lambda t: t.abs().amax(-1))
    lmin = lam.amin(-1)
    for name, f, b in [('sqrt', 'sqrt', 1 / lmin.sqrt()), ('exp', 'exp', torch.exp(lam.amax(-1))),
                       ('log', 'log', 1 / lmin), ('inverse', 'inverse', 1 / lmin ** 2)]:
        y32, y64 = getattr(V32, f)(x32), getattr(V64, f)(xr)
        compare(y32, y64, 8 * U * (b * scale + y64.abs().amax(-1)), 2, f"spin factor {name}",
                norm=lambda t: t.abs().amax(-1))
    compare(V32.determinant(x32), V64.determinant(xr), 4 * U * scale ** 2, 2, "spin factor determinant", norm=torch.abs)
    assert V32.in_cone(x32).all()
    compare(V32.power(x32, 3), V64.power(xr, 3), 16 * U * scale ** 3, 2, "spin factor power",
            norm=lambda t: t.abs().amax(-1))
    M32 = V32.to_matrix(x32)
    assert_single(M32)
    compare(V32.from_matrix(M32), xr, 8 * U * scale, 2, "to/from matrix", norm=lambda t: t.abs().amax(-1))
    S = jordan.anticommuting_symmetries(5, dtype=torch.complex64)
    for i in range(5):
        for j in range(5):
            ac = (d64(S[i]) @ d64(S[j]) + d64(S[j]) @ d64(S[i])
                  - 2 * (i == j) * torch.eye(S.shape[-1], dtype=torch.complex128))
            check(fro(ac[None]), U, 1, "anticommutation")
    a32, a64 = algs()
    X32, X64 = pair(a32, a64, general_matrix(a64, seed=38))
    Y32, Y64 = pair(a32, a64, general_matrix(a64, seed=39))
    ab = lambda A, B: fro(A.matrix.abs() @ B.matrix.abs())
    compare(jordan.jordan_product(X32, Y32), jordan.jordan_product(X64, Y64), 2 * KMAX * U * ab(X64, Y64), 2,
            "jordan product")
    compare(jordan.quadratic_representation(X32, Y32), jordan.quadratic_representation(X64, Y64),
            8 * KMAX * U * fro(X64.matrix.abs() @ Y64.matrix.abs() @ X64.matrix.abs()), 2, "quadratic repr")
    compare(jordan.jordan_triple(X32, Y32, X32), jordan.jordan_triple(X64, Y64, X64),
            8 * KMAX * U * fro(X64.matrix.abs() @ Y64.matrix.abs() @ X64.matrix.abs()), 2, "jordan triple")
    compare(jordan.jordan_power(X32, 4), jordan.jordan_power(X64, 4),
            3 * KMAX * U * fro(torch.linalg.matrix_power(X64.matrix.abs(), 4)), 2, "jordan power")
    compare(jordan.transpose_map(X32), jordan.transpose_map(X64), 0.0 + 1e-300, 0, "transpose map")
    gens = torch.stack([torch.diag(torch.tensor([1.0, 2.0, 2.0, 3.0])),
                        torch.tensor([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=torch.float32)])
    J32 = jordan.jordan_algebra_generated(gens)
    J64 = jordan.jordan_algebra_generated(gens.double())
    assert J32.dim == J64.dim and J32.basis.dtype in SINGLE


def test_inclusions_single():
    """Inclusions of single-precision algebras: embed (exact), conditional expectation (sums: C m u), Watatani index
    (formula and quasi-basis), Pimsner-Popa constants, Jones projection (e^2 = e to C u), the U(1) chain tower, the
    relative commutant ([iota(N), R] = 0 to C u)."""
    torch.manual_seed(40)
    tower32 = incl.chain_tower(3, precision='single')
    tower64 = incl.chain_tower(3, precision='double')
    for I32, I64 in zip(tower32, tower64):
        X64 = I64.N.operator(general_matrix(I64.N, B=2, seed=42))
        X32 = X64.cast(I32.N)
        Xr = X32.cast(I64.N)
        compare(I32.embed(X32), I64.embed(Xr), 1e-300, 0, "embed")
        for tr in ('blunt', 'markov'):
            E32, E64 = I32.conditional_expectation(tr), I64.conditional_expectation(tr)
            Y32 = I32.embed(X32) @ I32.embed(X32).adjoint()
            Y64 = I64.embed(Xr) @ I64.embed(Xr).adjoint()
            compare(E32(Y32), E64(Y64), 8 * I64.M.k_max * U * fro(Y64.matrix.abs()), 2, f"E[{tr}]")
            compare(E32.index(), E64.index(), 4 * U * E64.index().matrix.abs().amax(), 2, f"index[{tr}]",
                    norm=lambda t: t.abs().amax((-3, -2, -1)))
            compare(E32.index('quasi_basis'), E64.index(), 8 * U * E64.index().matrix.abs().amax() * I64.M.k_max, 2,
                    f"index qb[{tr}]", norm=lambda t: t.abs().amax((-3, -2, -1)))
            compare(E32.pimsner_popa_constant(), E64.pimsner_popa_constant(), 4 * U, 2, f"PP[{tr}]", norm=torch.abs)
        e = I32.jones_projection()
        assert_single(e)
        check(fro(d64(e) @ d64(e) - d64(e)), 4 * U * e.algebra.k_max, 4, "jones projection")
        R = I32.relative_commutant()
        check(R.commutator_error(I32), 8 * U * I64.M.k_max, 4, "relative commutant")


def test_optimize_single():
    """optimize in single precision: UnitaryParam stays unitary to C k u after Adam steps, SpectrumParam keeps
    lambda in [lo, 1] with both ends attained, and extremize reaches the von Neumann trace bound
    max Re Tr(A U) = ||A||_1 to single-precision accuracy."""
    torch.manual_seed(43)
    a32, a64 = algs([3])
    A = a32.operator(general_matrix(a64, B=1, seed=44).to(torch.complex64).expand(4, -1, -1, -1).clone())
    up = optimize.UnitaryParam(a32, 4)
    best = optimize.extremize(lambda: torch.diagonal(A.matrix @ up.matrix(), dim1=-2, dim2=-1).sum((-2, -1)).real,
                              [up], maximize=True, rounds=4, steps=40, lr=0.1)
    assert_single(best, up.matrix())
    Um = d64(up.matrix())
    check(fro(Um.mH @ Um - torch.eye(3, dtype=torch.complex128)), 3 * U * 10, 4, "UnitaryParam unitary")
    target = d64(A.cast(a64).trace_norm())
    check((d64(best) - target).abs(), 1e-4 * target, 4, "extremize")
    P = optimize.PositiveParam(a32, 4, delta=torch.full((4,), 0.5))
    lam = P.spectrum.eigenvalues()
    assert_single(lam)
    lo = (1 - 0.5) / 1.5
    act = lam[P.spectrum.mask.expand_as(lam)]
    assert (act >= lo * (1 - 2 * U)).all() and (act <= 1 + 2 * U).all()
    X = P.operator()
    assert_single(X)
    Xd = a64.operator(d64(X), is_self_adjoint=True)
    check((Xd.lambda_max - 1).abs(), 3 * U * 4, 4, "PositiveParam max")
    check((Xd.lambda_min - lo).abs(), 3 * U * 4, 4, "PositiveParam min")
    S = optimize.SelfAdjointParam(a32, 4)
    B = a32.operator(general_matrix(a64, B=1, seed=45).to(torch.complex64))
    v = optimize.extremize(lambda: torch.diagonal(S.operator().matrix @ B.matrix, dim1=-2, dim2=-1).sum((-2, -1)).real,
                           [S], rounds=2, steps=10)
    assert_single(v, S.operator())
    assert torch.isfinite(v).all()


# ======================================================================
# J. divergences
# ======================================================================
DIV_CASES = [('petz', 0.5), ('petz', 2.0), ('petz', 10.0), ('sandwiched', 0.5), ('sandwiched', 2.0),
             ('sandwiched', 50.0), ('umegaki', None), ('max', None), ('min', None)]


def _div(kind, a):
    if kind == 'petz':
        return lambda r, s: dv.petz_renyi_divergence(r, s, a)
    if kind == 'sandwiched':
        return lambda r, s: dv.sandwiched_renyi_divergence(r, s, a)
    return {'umegaki': dv.umegaki_relative_entropy, 'max': dv.max_relative_entropy,
            'min': dv.min_relative_entropy}[kind]


@pytest.mark.parametrize("family", ['generic', 'near_degenerate'])
@pytest.mark.parametrize("kind,alpha", DIV_CASES)
def test_divergences_single(kind, alpha, family):
    """Renyi (Petz, sandwiched; alpha up to 50), Umegaki, D_max, D_0 in single precision: |dD| <= C times the
    first-order bound u (k ||rho|| ||grad_rho D||_1 + k ||sigma|| ||grad_sigma D||_1 + |D|)."""
    a32, a64 = algs()
    fam = POS_FAMILIES['pos'] if family == 'generic' else POS_FAMILIES['pos_gap1e-6']
    mr = spectra_matrix(a64, fam, B=3, seed=45)
    ms = spectra_matrix(a64, fam, B=3, seed=46)
    mr = mr / torch.diagonal(mr, dim1=-2, dim2=-1).sum((-2, -1)).real[:, None, None, None]
    ms = ms / torch.diagonal(ms, dim1=-2, dim2=-1).sum((-2, -1)).real[:, None, None, None]
    r32, r64 = pair_state(a32, a64, mr)
    s32, s64 = pair_state(a32, a64, ms)
    D = _div(kind, alpha)
    op = lambda m: a64.operator(m, is_self_adjoint=True, is_positive=True)
    b = grad_bound(lambda r, s: D(op(r), op(s)), [r64.matrix, s64.matrix])
    compare(D(r32, s32), D(r64, s64), b, 4, f"{kind}{alpha}[{family}]", norm=torch.abs)


def test_divergence_supports_and_extras_single():
    """Support conventions in single precision (+inf exactly where double gives +inf: supp rho not in supp sigma for
    alpha > 1, Umegaki, D_max; orthogonal supports for alpha < 1), sector_contributions, renyi_divergence limits,
    lp_norm (also at norms 1e-30 .. 1e30), kosaki_norm, Petz and universal recovery maps."""
    a32, a64 = algs()
    rk = lambda k, rng: np.concatenate([[1.0, 0.5], np.zeros(k - 2)])[:k] if k > 2 else np.ones(k)
    ms = spectra_matrix(a64, rk, B=2, seed=47)
    ms = ms / torch.diagonal(ms, dim1=-2, dim2=-1).sum((-2, -1)).real[:, None, None, None]
    s32, s64 = pair_state(a32, a64, ms)                       # singular sigma
    r32, r64 = states_pair(a32, a64, B=2, seed=48, lam_min=0.02)  # faithful rho: supp rho not in supp sigma
    for kind, alpha in [('petz', 2.0), ('sandwiched', 3.0), ('umegaki', None), ('max', None)]:
        v32, v64 = _div(kind, alpha)(r32, s32), _div(kind, alpha)(r64, s64)
        assert torch.isinf(v64).all() and torch.isinf(v32).all(), kind
    for kind, alpha in [('petz', 0.5), ('sandwiched', 0.7)]:
        v32, v64 = _div(kind, alpha)(r32, s32), _div(kind, alpha)(r64, s64)
        assert torch.isfinite(v32).all() and torch.isfinite(v64).all()
    g32, g64 = states_pair(a32, a64, B=2, seed=49, lam_min=0.02)
    op = lambda m: a64.operator(m, is_self_adjoint=True, is_positive=True)
    for kind in ('petz', 'sandwiched', 'umegaki'):
        c32 = dv.sector_contributions(r32, g32, 1.5, kind=kind)
        c64 = dv.sector_contributions(r64, g64, 1.5, kind=kind)
        assert_single(c32.terms, c32.conditional)
        b = grad_bound(lambda r, s: dv.sector_contributions(op(r), op(s), 1.5, kind=kind).conditional.abs().sum(-1),
                       [r64.matrix, g64.matrix])
        check((d64(c32.conditional) - c64.conditional).abs().sum(-1), b, 4, f"sector conditional[{kind}]")
    for alpha, kind in [(1.0, 'petz'), (math.inf, 'sandwiched'), (0.0, 'petz')]:
        compare(dv.renyi_divergence(r32, g32, alpha, kind=kind), dv.renyi_divergence(r64, g64, alpha, kind=kind),
                grad_bound(lambda r, s: dv.renyi_divergence(op(r), op(s), alpha, kind=kind), [r64.matrix, g64.matrix]),
                4, f"renyi limit {alpha}", norm=torch.abs)
    for scale in (1e-30, 1.0, 1e30):
        X32, X64 = pair(a32, a64, general_matrix(a64, seed=50, scale=scale))
        for p in (1.0, 3.0, math.inf):
            compare(dv.lp_norm(X32, p, 'norm'), dv.lp_norm(X64, p, 'norm'), sum(NS) * KMAX * U * op_norm(X64), 2,
                    f"lp_norm[{p}, {scale:g}]", norm=torch.abs)
    X32, X64 = pair(a32, a64, general_matrix(a64, B=2, seed=51))
    for eta in (0.0, 0.5, 1.0):
        compare(dv.kosaki_norm(X32, r32, 2.5, eta), dv.kosaki_norm(X64, r64, 2.5, eta),
                grad_bound(lambda x, r: dv.kosaki_norm(a64.operator(x), op(r), 2.5, eta), [X64.matrix, r64.matrix]), 4,
                f"kosaki[{eta}]", norm=torch.abs)
    torch.manual_seed(52)
    c32 = ch.random_channel(a32, 2)
    c64 = ch.Channel(a64, d64(c32.kraus))
    kappa = lambda X: (lambda w: w.amax(-1) / w.amin(-1))(herm_eigs(X))     # condition numbers of sigma, Phi(sigma)
    R32 = dv.petz_recovery(g32, c32, c32.apply(r32))
    R64 = dv.petz_recovery(g64, c64, c64.apply(r64))
    b = KMAX * U * (kappa(g64) + kappa(c64.apply(g64)) + 4) * fro(R64.matrix)
    compare(R32, R64, b, 2, "petz recovery")
    U32 = dv.universal_recovery(g32, c32, c32.apply(r32), step=0.5, cutoff=3.0)
    U64 = dv.universal_recovery(g64, c64, c64.apply(r64), step=0.5, cutoff=3.0)
    compare(U32, U64, b, 2, "universal recovery")


# ======================================================================
# O. no NaN / Inf in single precision where double precision is finite
# ======================================================================
ADVERSARIAL = {
    'zero': lambda k, rng: np.zeros(k),
    'rank1': lambda k, rng: np.concatenate([[1.0], np.zeros(k - 1)]),
    'tiny': lambda k, rng: rng.uniform(0.2, 2.0, k) * 1e-30,
    'huge': lambda k, rng: rng.uniform(0.2, 2.0, k) * 1e30,
    'range': lambda k, rng: np.logspace(-30, 0, k),
    'degenerate': lambda k, rng: np.ones(k),
    'near_zero': lambda k, rng: np.concatenate([[1e-9, 1e-8], np.ones(max(k - 2, 0))])[:k],
}
NAN_FUNCS = {
    'eigenvalues': lambda A, S: torch.cat(A.eigenvalues(), -1),
    'sqrt': lambda A, S: A.sqrt().matrix,
    'log': lambda A, S: A.log().matrix,
    'power0.5': lambda A, S: A.power(0.5).matrix,
    'power-0.5': lambda A, S: A.power(-0.5).matrix,
    'power3': lambda A, S: A.power(3.0).matrix,
    'inverse': lambda A, S: A.inverse().matrix,
    'abs': lambda A, S: A.abs().matrix,
    'polar': lambda A, S: torch.cat([m.matrix for m in A.polar()], 1),
    'exp': lambda A, S: A.apply_function(lambda w: torch.exp(-w / w.abs().amax().clamp(min=1e-30))).matrix,
    'entropy': lambda A, S: A.entropy(),
    'vn_entropy': lambda A, S: st.von_neumann_entropy(A),
    'trace_norm': lambda A, S: A.trace_norm(),
    'frobenius': lambda A, S: A.frobenius_norm(),
    'operator_norm': lambda A, S: A.operator_norm(),
    'contrast': lambda A, S: A.michelson_contrast,
    'fidelity': lambda A, S: st.fidelity(A, S),
    'trace_distance': lambda A, S: st.trace_distance(A, S),
    'relative_entropy': lambda A, S: st.relative_entropy(S, A),
    'petz0.5': lambda A, S: dv.petz_renyi_divergence(A, S, 0.5),
    'petz2': lambda A, S: dv.petz_renyi_divergence(S, A, 2.0),
    'sandwiched0.5': lambda A, S: dv.sandwiched_renyi_divergence(A, S, 0.5),
    'sandwiched30': lambda A, S: dv.sandwiched_renyi_divergence(S, A, 30.0),
    'umegaki': lambda A, S: dv.umegaki_relative_entropy(S, A),
    'max_rel': lambda A, S: dv.max_relative_entropy(S, A),
    'min_rel': lambda A, S: dv.min_relative_entropy(A, S),
    'lp3': lambda A, S: dv.lp_norm(A, 3.0),
    'kosaki': lambda A, S: dv.kosaki_norm(S, A, 2.0),
    'gibbs': lambda A, S: st.gibbs_state(A, 1.0 / A.operator_norm().clamp(min=1e-30).amax().item()).matrix,
    'log_Z': lambda A, S: st.log_partition_function(A, 1.0),
    'propagator': lambda A, S: dyn.propagator(A, 1.0).matrix,
}


@pytest.mark.parametrize("family", list(ADVERSARIAL))
def test_no_nan_in_single_where_double_finite(family):
    """Adversarial positive inputs (zero, rank one, norms 1e-30 / 1e30, dynamic range 1e-30 .. 1, exact degeneracy,
    eigenvalues 1e-9): every function gives finite single-precision values wherever the double values are finite and
    representable in float32 (and raises in single only if it raises in double)."""
    a32, a64 = algs()
    A32, A64 = pair(a32, a64, spectra_matrix(a64, ADVERSARIAL[family], B=2, seed=53), is_self_adjoint=True,
                    is_positive=True)
    S32, S64 = states_pair(a32, a64, B=2, seed=54, lam_min=0.01)
    for name, F in NAN_FUNCS.items():
        try:
            ref = F(A64, S64)
        except (RuntimeError, ValueError):
            with pytest.raises((RuntimeError, ValueError)):
                F(A32, S32)
            continue
        got = F(A32, S32)
        assert_single(got)
        try:
            assert_finite_where(got, ref)
        except AssertionError:
            raise AssertionError(f"{name}[{family}]: NaN/Inf in single precision where double is finite")


# ======================================================================
# sensitivity of the bounds
# ======================================================================
def test_bounds_detect_half_precision():
    """The bounds used above reject results computed in reduced precision: rounding the single-precision results
    through bfloat16 (as a silent cast would) violates them by orders of magnitude."""
    a32, a64 = algs()
    A32, A64 = pair(a32, a64, spectra_matrix(a64, POS_FAMILIES['pos'], seed=55), is_self_adjoint=True,
                    is_positive=True)
    bound = KMAX * U * op_norm(A64)
    w = torch.cat(A32.eigenvalues(), -1)
    w_half = w.to(torch.bfloat16).to(torch.float32)
    err = (d64(w_half) - torch.cat(A64.eigenvalues(), -1)).abs().amax(-1)
    assert float((err / bound).max()) > 100
    S = A32.sqrt().matrix
    err = fro(d64(S.to(torch.complex64).real.to(torch.bfloat16).to(torch.float32).to(torch.complex64)
                  + 1j * S.imag.to(torch.bfloat16).to(torch.float32)) - d64(A64.sqrt()))
    assert float((err / dk_bound(A64, torch.sqrt, lambda x: 0.5 / torch.sqrt(x))).max()) > 100


# ======================================================================
# further public functions in single precision
# ======================================================================
def test_hilbert_space_single():
    """HilbertSpace(precision='single'): inner products and norms of O(1) vectors to C n u, normalize to norm 1."""
    H32 = vn.HilbertSpace(5, 5, batch_size=3, n_channels=2, precision='single', device='cpu')
    psi = torch.randn(3, 2, 5, 1, dtype=torch.complex128, generator=torch.Generator().manual_seed(60))
    phi = torch.randn(3, 2, 5, 1, dtype=torch.complex128, generator=torch.Generator().manual_seed(61))
    p32, f32 = psi.to(torch.complex64), phi.to(torch.complex64)
    ip = H32.inner_product(p32, f32)
    assert_single(ip)
    ref = (d64(p32).conj() * d64(f32)).sum((-3, -2, -1))
    check((d64(ip).reshape(-1) - ref).abs(), 10 * U * (d64(p32).abs() * d64(f32).abs()).sum((-3, -2, -1)), 2,
          "inner product")
    n = H32.normalize(p32)
    assert_single(n)
    check((d64(H32.norm(n)).reshape(-1) - 1).abs(), 10 * U, 2, "normalize")


@pytest.mark.xfail(strict=True, reason=(
    "HilbertSpace.norm / normalize form sum |psi_i|^2 without scaling: in float32 the squares overflow for entries "
    ">~ 1.8e19 (norm Inf, normalize -> 0) although the norm is representable; frobenius_norm and lp_norm scale"))
def test_hilbert_space_norm_large_entries_single():
    """HilbertSpace.norm and normalize of vectors with entries 1e20 (norm 3e20, representable in float32)."""
    H32 = vn.HilbertSpace(4, 4, batch_size=1, n_channels=2, precision='single', device='cpu')
    psi = torch.full((1, 2, 4, 1), 1e20, dtype=torch.complex64)
    assert torch.isfinite(H32.norm(psi)).all()
    check((d64(H32.norm(H32.normalize(psi))).reshape(-1) - 1).abs(), 10 * U, 2, "normalize huge")


@pytest.mark.xfail(strict=True, reason=(
    "krylov.lanczos / _krylov_step normalise with torch.linalg.vector_norm, which overflows in float32 for entries "
    ">~ 1.8e19: evolve, expm_multiply_step and ground_state return NaN for a start vector that double handles"))
def test_krylov_large_norm_vectors_single():
    """krylov.evolve / ground_state with an unnormalised start vector of entries 1e20 (norm representable)."""
    c32 = SpinChain(6, precision='single')
    S = c32.xxz_sparse(1.0, 0.5, None, sector=3)
    psi = torch.full((1, S.dim), 1e20, dtype=torch.complex64)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = krylov.evolve(S, psi, [0.0, 0.5])
        E, _, _ = krylov.ground_state(S, v0=psi.real)
    assert torch.isfinite(out).all() and torch.isfinite(E).all()


def test_more_channels_and_recovery_single():
    """kraus_channel, measurement_channel, rotated Petz maps (t != 0) and the spectral cache in single precision."""
    a32, a64 = algs()
    r32, r64 = states_pair(a32, a64, B=2, seed=62, lam_min=0.02)
    P = torch.zeros(1, 3, KMAX, KMAX, dtype=torch.complex128)
    P[0, 0, 0, 0] = P[0, 0, 1, 1] = P[0, 1, 0, 0] = 1
    P32, P64 = pair(a32, a64, P, is_self_adjoint=True, is_projection=True)
    m32 = ch.measurement_channel([P32, a32.identity() - P32])
    m64 = ch.measurement_channel([P64, a64.identity() - P64])
    compare(m32.apply(r32), m64.apply(r64), 2 * 2 * KMAX * U, 2, "measurement",
            norm=lambda x: x.abs().amax((-3, -2, -1)))
    X32, X64 = pair(a32, a64, general_matrix(a64, B=2, seed=63, scale=0.3))
    k32 = ch.kraus_channel([X32, X32.adjoint()])
    k64 = ch.kraus_channel([X64, X64.adjoint()])
    compare(k32.apply(r32), k64.apply(r64), fro(_abs_apply(k64.kraus, r64.matrix)) * (2 * KMAX + 2) * U, 2,
            "kraus_channel")
    torch.manual_seed(64)
    c32 = ch.random_channel(a32, 2)
    c64 = ch.Channel(a64, d64(c32.kraus))
    kappa = lambda X: (lambda w: w.amax(-1) / w.amin(-1))(herm_eigs(X))
    R32, R64 = dv.PetzRecovery(r32, c32), dv.PetzRecovery(r64, c64)
    for t in (0.0, 0.7):
        Y32, Y64 = R32.rotated(c32.apply(r32), t), R64.rotated(c64.apply(r64), t)
        b = KMAX * U * (1 + t) * (kappa(r64) + kappa(c64.apply(r64)) + 4) * fro(Y64.matrix)
        compare(Y32, Y64, b, 2, f"rotated petz[{t}]")
    H32, H64 = pair(a32, a64, spectra_matrix(a64, POS_FAMILIES['pos'], seed=65), is_self_adjoint=True, is_positive=True)
    with vn.caching():
        first = [H32.log(), H32.sqrt(), H32.power(0.3)]
        again = [H32.log(), H32.sqrt(), H32.power(0.3)]
    for x, y, z in zip(first, again, [H64.log(), H64.sqrt(), H64.power(0.3)]):
        torch.testing.assert_close(x.matrix, y.matrix, rtol=0, atol=0)
        compare(x, z, 8 * KMAX * U * 10, 2, "cached spectral functions")


def test_standard_inclusions_single():
    """diagonal, scalars, centre, from_tensor_product, compose, basic construction / Jones tower, Bratteli diagram in
    single precision: E o iota = id to C m u, index as in double, Jones projections idempotent to C m u."""
    a32, a64 = algs([2, 3])
    pairs = [(incl.diagonal(4, precision='single'), incl.diagonal(4, precision='double')),
             (incl.scalars(a32), incl.scalars(a64)), (incl.centre(a32), incl.centre(a64)),
             (incl.from_tensor_product(a32, 2), incl.from_tensor_product(a64, 2))]
    I1, I2 = pairs[0][0], pairs[0][0].basic_construction()
    pairs.append((I1.compose(I2), pairs[0][1].compose(pairs[0][1].basic_construction())))
    for I32, I64 in pairs:
        X = general_matrix(I64.N, B=2, seed=66)
        X32 = I32.N.operator(X.to(I32.N.dtype))
        E = I32.conditional_expectation('blunt')
        back = E(I32.embed(X32))
        compare(back, X32.cast(I64.N), 4 * I64.M.k_max * U * X.abs().amax(), 2, "E o iota",
                norm=lambda t: t.abs().amax((-3, -2, -1)))
        assert I32.index() == pytest.approx(I64.index(), rel=1e-12)
        if I32.is_connected:
            e = I32.jones_projection()
            check(fro(d64(e) @ d64(e) - d64(e)), 4 * e.algebra.k_max * U, 4, "jones projection")
    tower = I1.jones_tower(3)
    assert [t.M.dtype for t in tower] == [torch.complex64] * 3
    assert len(incl.bratteli_diagram(tower).levels) == 4


def test_decompose_and_jordan_extras_single():
    """from_group / from_generators / collective_operators / conditional expectations / commutant maps of
    decompose, and star_algebra_self_adjoint_part / closure checks of jordan, with single-precision input."""
    perm32 = dc.permutation_representation(3, 2, dtype=torch.complex64)
    D32, D64 = dc.from_group(perm32), dc.from_group(dc.permutation_representation(3, 2))
    assert (D32.k_factors, D32.multiplicities) == (D64.k_factors, D64.multiplicities)
    a = torch.randn(8, 8, dtype=torch.complex128, generator=torch.Generator().manual_seed(67)).to(torch.complex64)
    for onto in ('algebra', 'commutant'):
        Ea = D32.conditional_expectation(a, onto)
        assert_single(Ea)
        check(fro((d64(D32.conditional_expectation(Ea, onto)) - d64(Ea))[None]), 8 * U * fro(d64(a)[None]), 4,
              f"E idempotent[{onto}]")
    Cm = D32.to_commutant(perm32)
    assert_single(Cm)
    check(fro(d64(D32.from_commutant(Cm)) - d64(D32.conditional_expectation(perm32, 'commutant'))), 8 * 8 * U, 4,
          "commutant round trip")
    G = dc.collective_operators(dc.spin_matrices(0.5, dtype=torch.complex64), 3)
    assert G.dtype == torch.complex64
    D = dc.from_generators(G)
    assert D.k_factors == dc.from_generators(dc.collective_operators(dc.spin_matrices(0.5), 3)).k_factors
    gens = torch.stack([torch.diag(torch.tensor([1.0, 2.0, 2.0, 3.0])),
                        torch.tensor([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=torch.float32)])
    J32 = jordan.jordan_algebra_generated(gens)
    assert J32.is_closed() and bool(J32.contains(gens).all())
    S32 = jordan.star_algebra_self_adjoint_part(gens)
    S64 = jordan.star_algebra_self_adjoint_part(gens.double())
    assert S32.dim == S64.dim


def test_property_checks_single():
    """The tolerance-based property checks (relative 1e-5 in single precision) classify rounded data correctly:
    a rounded projection, a rounded positive semidefinite operator with eigenvalues 0 (computed as -u), a rounded
    unitary (normal), and reject clear violations."""
    a32, a64 = algs()
    proj = lambda k, rng: np.concatenate([np.ones(k // 2 + 1), np.zeros(k)])[:k]
    P32, _ = pair(a32, a64, spectra_matrix(a64, proj, seed=70))
    assert P32.is_self_adjoint and P32.is_projection and P32.is_positive
    Z32, _ = pair(a32, a64, spectra_matrix(a64, FAMILIES['psd_zero'], seed=71))
    assert Z32.is_positive and not Z32.is_invertible
    Ham = a32.operator(spectra_matrix(a64, FAMILIES['generic'], seed=72).to(torch.complex64), is_self_adjoint=True)
    W = a32.operator(dyn.propagator(Ham, 2.0).matrix)
    assert W.is_normal and W.is_invertible and not W.is_self_adjoint
    N32, _ = pair(a32, a64, spectra_matrix(a64, FAMILIES['generic'], seed=73) - 3 * eye_of(a64).to(torch.complex128))
    assert N32.is_self_adjoint and not N32.is_positive
    G32, _ = pair(a32, a64, general_matrix(a64, seed=74))
    assert not G32.is_normal


def test_inter_sector_lindblad_and_su2_expectation_single():
    """lindblad_evolve with jumps between sectors (SpinChain.lowering) in single precision against double, and
    SU2Chain.expectation of H in its ground state = E_0 to C m u ||H||."""
    c32, c64 = SpinChain(3, precision='single'), SpinChain(3, precision='double')
    H32, H64 = c32.xxz(1.0, 0.5, 0.2), c64.xxz(1.0, 0.5, 0.2)
    r32, r64 = pair_state(c32.algebra, c64.algebra, st.random_density_matrix(c64.algebra, 2).matrix)
    o32 = dyn.lindblad_evolve(r32, H32, [c32.lowering(0)], times=[0.0, 1.0], substeps=20)[-1]
    o64 = dyn.lindblad_evolve(r64, H64, [c64.lowering(0)], times=[0.0, 1.0], substeps=20)[-1]
    compare(o32, o64, 2 * 20 * U + 3 * U * (1 + op_norm(H64)), 2, "lindblad inter-sector")
    check((torch.diagonal(d64(o32), dim1=-2, dim2=-1).sum((-2, -1)) - 1).abs(), 2 * 20 * U, 2, "trace kept")
    s32, s64 = SU2Chain(6, precision='single'), SU2Chain(6, precision='double')
    Hs32 = s32.heisenberg(1.0)
    g = s32.ground_state(Hs32)
    S = float(g.spin)
    m = s64.multiplicities[s64.sector_index(S)]
    e = s32.expectation(Hs32, S, g.coefficients[:, :m])
    assert_single(e)
    check((d64(e) - d64(g.energy)).abs(), m * U * op_norm(s64.heisenberg(1.0)), 4, "su2 expectation")
