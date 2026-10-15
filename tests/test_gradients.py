"""
Gradients of the public API in double precision.

* torch.autograd.gradcheck for the differentiable functions of every module, at generic points and at nearly
  degenerate spectra (relative gaps 1e-6 .. 1e-12), where the Daleckii-Krein formula with the midpoint rule for
  close eigenvalues must stay exact.
* Directional derivatives against mpmath (50-80 digits, central differences) where gradcheck (atol 1e-5) is too
  coarse: near-degenerate spectra, logarithms and square roots near 0, the parameter alpha of the Renyi divergences.
* Rank-deficient inputs. The maps |X|, X^+, the polar factor |X|, A^{1/2}, the entropy and the divergences are not
  differentiable at a singular X in every direction. With the SVD X = W diag(S, 0) V^* and E~ = W^* E V split into
  blocks (range, kernel), the true behaviour is:
    - |X| (and the factor |X| of the polar decomposition): |X + tE| = |X| + t Lin(E) + |t| N(E) + o(t) with
      N(E) = V_0 |E~_22| V_0^*, so the one-sided derivatives differ when E~_22 != 0 but the symmetric derivative
      lim (|X + tE| - |X - tE|) / 2t = Lin(E) exists in every direction; the library returns exactly Lin(E).
    - X^+, A^{1/2}, A^p, log A, S(rho), D(rho||sigma): differentiable along rank-preserving curves
      (constant-rank parametrisations X = A B^*, rho = Y Y^* / Tr), where the library must give the exact derivative;
      along directions E~_22 != 0 the derivative does not exist (X^+ jumps, sqrt and entropy have infinite one-sided
      slope), where only finiteness of the returned gradient is required.
Genuine library bugs are marked xfail(strict=True) with the reason.
"""
import math

import numpy as np
import pytest
import torch
from torch.autograd import gradcheck, gradgradcheck

import torch_vn_algebra as vn
from torch_vn_algebra import channels as ch
from torch_vn_algebra import composite, divergences as dv, dynamics as dyn, functional as fc, jordan
from torch_vn_algebra import optimize, states as st
from torch_vn_algebra.chains import SpinChain
from torch_vn_algebra.su2 import SU2Chain

CD = torch.complex128
GAPS = [1e-6, 1e-9, 1e-12]


# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------
def _R(t):
    """Real view of a complex tensor (gradcheck compares real Jacobians)."""
    return torch.view_as_real(t) if torch.is_complex(t) else t


def gc(f, *inputs, eps=1e-6, atol=1e-5, rtol=1e-3, **kw):
    """gradcheck of f(*inputs) with complex outputs viewed as real."""
    return gradcheck(lambda *a: _R(f(*a)), inputs, eps=eps, atol=atol, rtol=rtol, **kw)


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


def herm(x):
    return 0.5 * (x + x.mH)


def haar(k, gen, dtype=CD):
    Z = torch.randn(k, k, dtype=dtype, generator=gen)
    Q, R = torch.linalg.qr(Z)
    d = torch.diagonal(R)
    return Q * (d / d.abs()).unsqueeze(-2)


def spectral_matrix(alg, spectra, seed=0, dtype=CD):
    """(1, C, k_max, k_max) tensor with Hermitian blocks U diag(spectra[c]) U^* (Haar U), padding zero."""
    gen = torch.Generator().manual_seed(seed)
    X = torch.zeros(1, alg.C, alg.k_max, alg.k_max, dtype=dtype)
    for c, (w, k) in enumerate(zip(spectra, alg.k_factors)):
        U = haar(k, gen, dtype)
        X[0, c, :k, :k] = (U * torch.as_tensor(w, dtype=torch.float64).to(dtype)) @ U.mH
    return herm(X)


def cplx(shape, seed, scale=1.0):
    gen = torch.Generator().manual_seed(seed)
    return scale * torch.randn(*shape, dtype=CD, generator=gen)


def hop(alg, x, **flags):
    """Self-adjoint operator from a raw tensor (the Hermitian part: perturbations of gradcheck stay Hermitian)."""
    return alg.operator(herm(x), is_self_adjoint=True, **flags)


def directional(F, A0, E):
    """d/dt F(A0 + t E) at t = 0 by reverse-mode autograd through a real scalar t."""
    t = torch.zeros((), dtype=torch.float64, requires_grad=True)
    out = F(A0 + t * E)
    g, = torch.autograd.grad(out, t)
    return g


@pytest.fixture(scope="module")
def alg():
    """Two sectors of different size (padding) and one of size 1."""
    return vn.TypeIAlgebra([3, 2, 1], [3, 2, 1], precision='double')


@pytest.fixture(scope="module")
def alg_s():
    """A small algebra with padding (fast gradchecks of the expensive functionals)."""
    return vn.TypeIAlgebra([2, 1], [2, 1], precision='double')


@pytest.fixture(scope="module")
def alg4():
    return vn.TypeIAlgebra([4], [4], precision='double')


# ----------------------------------------------------------------------
# mpmath references
# ----------------------------------------------------------------------
def _mp():
    return pytest.importorskip("mpmath")


def mp_matrix(A):
    mp = _mp()
    A = np.asarray(A)
    return mp.matrix([[mp.mpc(complex(z).real, complex(z).imag) for z in row] for row in A])


def mp_hfun(M, f):
    """f(M) for a Hermitian mpmath matrix (f a function of mpf eigenvalues)."""
    mp = _mp()
    Mh = (M + M.transpose_conj()) * mp.mpf('0.5')
    E, Q = mp.eighe(Mh)
    return Q * mp.diag([f(e) for e in E]) * Q.transpose_conj()


def mp_re_inner(G, M):
    """Re sum conj(G_ij) M_ij for a numpy G and an mpmath matrix M."""
    mp = _mp()
    s = mp.mpf(0)
    for i in range(M.rows):
        for j in range(M.cols):
            s += mp.re(mp.conj(mp.mpc(complex(G[i, j]).real, complex(G[i, j]).imag)) * M[i, j])
    return s


def mp_central(F, dps=50, h='1e-20'):
    """(F(h) - F(-h)) / 2h at dps digits for a function F(t) of an mpf t returning an mpf."""
    mp = _mp()
    with mp.workdps(dps):
        hh = mp.mpf(h)
        return float((F(hh) - F(-hh)) / (2 * hh))


# ======================================================================
# 1. functional calculus of the algebra: gradcheck, generic and near-degenerate
# ======================================================================
SPECTRAL = {
    'exp': lambda A: A.apply_function(torch.exp),
    'sin': lambda A: A.apply_function(torch.sin),
    'cube': lambda A: A.apply_function(lambda w: w ** 3),
    'sqrt': lambda A: A.sqrt(),
    'log': lambda A: A.log(),
    'power0.3': lambda A: A.power(0.3),
    'power-0.7': lambda A: A.power(-0.7),
    'power2.5': lambda A: A.power(2.5),
    'inverse': lambda A: A.inverse(),
    'abs_sa': lambda A: A.abs(),
    'expm': lambda A: A.expm(0.5),
}
POSITIVE = {'sqrt', 'log', 'power0.3', 'power-0.7', 'power2.5'}


@pytest.mark.parametrize("name", list(SPECTRAL))
def test_spectral_functions_gradcheck(alg, name):
    """Daleckii-Krein backward of the spectral calculus (and expm) at a generic spectrum (gradcheck)."""
    spectra = [[0.4, 1.3, 2.2], [0.7, 1.9], [1.1]] if name in POSITIVE else [[-1.2, 0.4, 2.0], [-0.6, 1.5], [0.9]]
    x = spectral_matrix(alg, spectra, seed=1).requires_grad_(True)
    flags = {'is_positive': True} if name in POSITIVE else {}
    assert gc(lambda m: SPECTRAL[name](hop(alg, m, **flags)).matrix, x)


@pytest.mark.parametrize("gap", GAPS)
@pytest.mark.parametrize("name", ['exp', 'sin', 'sqrt', 'log', 'power-0.7', 'inverse'])
def test_spectral_functions_gradcheck_near_degenerate(alg, name, gap):
    """Daleckii-Krein backward at nearly degenerate spectra (relative gaps 1e-6 .. 1e-12) (gradcheck)."""
    if name in POSITIVE:
        spectra = [[0.8, 0.8 * (1 + gap), 2.0], [1.5, 1.5 * (1 + gap)], [0.8 * (1 + 2 * gap)]]
    else:
        spectra = [[-0.8, 1.2, 1.2 * (1 + gap)], [-0.8 * (1 + gap), 0.5], [1.2 * (1 - gap)]]
    x = spectral_matrix(alg, spectra, seed=2).requires_grad_(True)
    flags = {'is_positive': True} if name in POSITIVE else {}
    assert gc(lambda m: SPECTRAL[name](hop(alg, m, **flags)).matrix, x)


def test_exactly_degenerate_and_scalar_blocks(alg):
    """f(A) at exactly degenerate spectra (a multiple of the identity, a repeated eigenvalue): gradcheck."""
    x = spectral_matrix(alg, [[1.0, 1.0, 1.0], [2.0, 2.0], [0.5]], seed=3).requires_grad_(True)
    for f in (lambda A: A.apply_function(torch.exp), lambda A: A.log(), lambda A: A.power(-0.5)):
        assert gc(lambda m: f(hop(alg, m, is_positive=True)).matrix, x)
    y = alg.identity().matrix.clone().to(CD).requires_grad_(True)
    assert gc(lambda m: hop(alg, m, is_positive=True).sqrt().matrix, y)


def test_normal_operator_functions_gradcheck():
    """apply_function of a normal (unitary) operator: holomorphic f, derivative along all directions (gradcheck)."""
    alg = vn.TypeIAlgebra([3, 2], [3, 2], precision='double')
    H = spectral_matrix(alg, [[-1.0, 0.3, 1.4], [0.2, 0.9]], seed=4)
    U0 = dyn.propagator(alg.operator(H, is_self_adjoint=True), 1.0).matrix.detach()
    E = cplx(U0.shape, 5) * mask_of(alg)
    t = torch.zeros((), dtype=torch.float64, requires_grad=True)

    def f(t):                                    # unitary curve U0 exp(i t K), K Hermitian: normal for all t
        K = herm(E)
        V = alg.operator(1j * K * t).expm().matrix
        U = alg.operator(U0 @ V, is_normal=True)
        return U.apply_function(lambda z: z ** 2 + 0.5 * z).matrix
    assert gc(f, t)


def test_holomorphic_function_gradcheck(alg_s):
    """holomorphic_function (contour integral) of a non-normal operator: gradcheck against finite differences."""
    alg = alg_s
    x = (cplx((1, 2, 2, 2), 6, 0.4) * mask_of(alg)).requires_grad_(True)
    assert gc(lambda m: alg.operator(m).holomorphic_function(torch.exp).matrix, x, atol=1e-6)
    assert gc(lambda m: alg.operator(m).expm().matrix, x, atol=1e-6)


@pytest.mark.parametrize("name", ['abs', 'inverse', 'polar_P'])
def test_general_blocks_gradcheck(alg, name):
    """|X|, X^{-1} and the factor |X| of the polar decomposition of non-normal invertible blocks (gradcheck)."""
    x = (cplx((1, 3, 3, 3), 7) * mask_of(alg) + 2 * eye_of(alg)).requires_grad_(True)
    f = {'abs': lambda A: A.abs(), 'inverse': lambda A: A.inverse(), 'polar_P': lambda A: A.polar()[1]}[name]
    assert gc(lambda m: f(alg.operator(m)).matrix, x)


def test_polar_unitary_factor_carries_no_gradient(alg):
    """The polar factor W is documented to carry no gradient; |X| does."""
    x = (cplx((1, 3, 3, 3), 8) * mask_of(alg)).requires_grad_(True)
    W, P = alg.operator(x).polar()
    assert not W.matrix.requires_grad and P.matrix.requires_grad


@pytest.mark.parametrize("name", ['eig_sum_exp', 'trace_norm', 'frobenius', 'operator_norm', 'lambda_max', 'lambda_min',
                                  'contrast', 'Tr_norm', 'tau_vN', 'trace', 'entropy', 'trace_a_log_a'])
def test_scalar_functionals_gradcheck(alg, name):
    """Scalar spectral functionals of the algebra (gradcheck at simple extreme eigenvalues)."""
    pos = {'contrast', 'entropy', 'trace_a_log_a'}
    spectra = [[0.1, 0.2, 0.3], [0.15, 0.05], [0.2]] if name in pos else [[-1.3, 0.2, 0.7], [-0.4, 1.6], [0.9]]
    x = spectral_matrix(alg, spectra, seed=9).requires_grad_(True)
    F = {
        'eig_sum_exp': lambda A: sum(torch.exp(w).sum(-1) for w in A.eigenvalues()),
        'trace_norm': lambda A: A.trace_norm(),
        'frobenius': lambda A: A.frobenius_norm(),
        'operator_norm': lambda A: A.operator_norm(),
        'lambda_max': lambda A: A.lambda_max,
        'lambda_min': lambda A: A.lambda_min,
        'contrast': lambda A: A.michelson_contrast,
        'Tr_norm': lambda A: A.Tr_norm().real,
        'tau_vN': lambda A: A.tau_vN().real,
        'trace': lambda A: A.trace.real,
        'entropy': lambda A: A.entropy(),
        'trace_a_log_a': lambda A: A.trace_a_log_a(),
    }[name]
    flags = {'is_positive': True} if name in pos else {}
    assert gc(lambda m: F(hop(alg, m, **flags)), x)


def test_operator_norm_of_non_normal_gradcheck(alg):
    """operator_norm and trace_norm of non-normal operators (singular values, simple largest one)."""
    x = (cplx((1, 3, 3, 3), 10) * mask_of(alg)).requires_grad_(True)
    assert gc(lambda m: alg.operator(m).operator_norm(), x)
    assert gc(lambda m: alg.operator(m).trace_norm(), x)


def test_second_derivatives_non_degenerate(alg_s):
    """Second derivatives (create_graph) of apply_function at a non-degenerate spectrum (gradgradcheck)."""
    alg = alg_s
    x = spectral_matrix(alg, [[0.3, 1.1], [0.6]], seed=11).requires_grad_(True)
    for f in (lambda A: A.apply_function(torch.exp), lambda A: A.log()):
        assert gradgradcheck(lambda m: _R(f(hop(alg, m, is_positive=True)).matrix), (x,), eps=1e-6, atol=1e-5,
                             rtol=1e-3)


def test_gradients_with_spectral_cache(alg):
    """Gradients are unchanged when the spectral cache is switched on (vn.caching())."""
    x = spectral_matrix(alg, [[0.3, 1.1, 2.4], [0.6, 1.7], [1.2]], seed=12).requires_grad_(True)
    G = cplx(x.shape, 13)
    outs = []
    for cache in (False, True):
        with vn.caching(cache):
            A = hop(alg, x, is_positive=True)
            s = (A.log().matrix * G.conj()).real.sum() + (A.sqrt().matrix * G).real.sum()
            outs.append(torch.autograd.grad(s, x)[0])
    torch.testing.assert_close(outs[0], outs[1], rtol=1e-12, atol=1e-12)


def test_cast_is_differentiable():
    """Operator.cast to another field / precision keeps the gradient (gradcheck through complex -> complex)."""
    a = vn.TypeIAlgebra([2, 3], [2, 3], precision='double')
    b = a.like(complex_valued=True)
    x = cplx((1, 2, 3, 3), 14).requires_grad_(True)
    assert gc(lambda m: a.operator(m).cast(b).matrix, x)


# ======================================================================
# 2. near-degenerate spectra, sqrt / log near 0: mpmath directional derivatives
# ======================================================================
MP_FUNS = {
    'exp': (lambda A: A.apply_function(torch.exp), lambda e, mp: mp.exp(e)),
    'log': (lambda A: A.log(), lambda e, mp: mp.log(e)),
    'sqrt': (lambda A: A.sqrt(), lambda e, mp: mp.sqrt(e)),
    'power-0.5': (lambda A: A.power(-0.5), lambda e, mp: e ** mp.mpf(-0.5)),
}


@pytest.mark.parametrize("gap", [1e-4] + GAPS)
@pytest.mark.parametrize("name", list(MP_FUNS))
def test_daleckii_krein_vs_mpmath_near_degenerate(alg4, name, gap):
    """Directional derivative of f(A) at relative gap 1e-4 .. 1e-12 agrees with mpmath to 1e-9 (relative)."""
    mp = _mp()
    lib, mpf = MP_FUNS[name]
    A0 = spectral_matrix(alg4, [[0.5, 0.5 * (1 + gap), 1.7, 1.7 * (1 + gap)]], seed=15)
    E = herm(cplx(A0.shape, 16)) * mask_of(alg4)
    G = cplx(A0.shape, 17)
    d_lib = float(directional(lambda A: (lib(alg4.operator(A, is_self_adjoint=True, is_positive=True)).matrix
                                         * G.conj()).real.sum(), A0, E))
    A_mp, E_mp, Gn = mp_matrix(A0[0, 0].numpy()), mp_matrix(E[0, 0].numpy()), G[0, 0].numpy()
    d_mp = mp_central(lambda t: mp_re_inner(Gn, mp_hfun(A_mp + E_mp * t, lambda e: mpf(e, mp))))
    assert abs(d_lib - d_mp) <= 1e-9 * max(1.0, abs(d_mp)), (d_lib, d_mp)


@pytest.mark.parametrize("lam0", [1e-3, 1e-6, 1e-9])
@pytest.mark.parametrize("name", ['log', 'sqrt', 'power-0.5'])
def test_log_sqrt_near_zero_vs_mpmath(alg4, name, lam0):
    """f(A) with eigenvalues down to 1e-9 (log, sqrt, A^{-1/2}): directional derivative vs mpmath to 1e-7 relative."""
    mp = _mp()
    lib, mpf = MP_FUNS[name]
    A0 = spectral_matrix(alg4, [[lam0, 2 * lam0, 0.3, 1.0]], seed=18)
    E = herm(cplx(A0.shape, 19)) * lam0                    # a direction that keeps A positive
    G = cplx(A0.shape, 20)
    d_lib = float(directional(lambda A: (lib(alg4.operator(A, is_self_adjoint=True, is_positive=True)).matrix
                                         * G.conj()).real.sum(), A0, E))
    A_mp, E_mp, Gn = mp_matrix(A0[0, 0].numpy()), mp_matrix(E[0, 0].numpy()), G[0, 0].numpy()
    d_mp = mp_central(lambda t: mp_re_inner(Gn, mp_hfun(A_mp + E_mp * t, lambda e: mpf(e, mp))), dps=60, h='1e-30')
    assert abs(d_lib - d_mp) <= 1e-7 * max(1.0, abs(d_mp)), (d_lib, d_mp)


@pytest.mark.parametrize("t", [1e2, 1e4])
@pytest.mark.parametrize("gap", [5e-6, 1e-4])
def test_propagator_gradient_long_times_vs_mpmath(alg4, t, gap):
    """d/dH of Re Tr G^* exp(-itH) at long times with nearly degenerate H (t |l_i - l_j| from 5e-3 to 100): the
    Loewner matrix of the oscillating f = exp(-itw) must not use the midpoint rule where its error
    |l_i - l_j|^2 t^3 / 12 is large (mpmath, relative 1e-7 of the scale t)."""
    mp = _mp()
    A0 = spectral_matrix(alg4, [[10.0, 10.0 * (1 + gap), 3.0, -2.0]], seed=21)
    E = herm(cplx(A0.shape, 22))
    G = cplx(A0.shape, 23)
    d_lib = float(directional(lambda A: (dyn.propagator(alg4.operator(A, is_self_adjoint=True), t).matrix
                                         * G.conj()).real.sum(), A0, E))
    A_mp, E_mp, Gn = mp_matrix(A0[0, 0].numpy()), mp_matrix(E[0, 0].numpy()), G[0, 0].numpy()
    with mp.workdps(60):
        tt = mp.mpf(t)

        def val(h):
            M = A_mp + E_mp * h
            ev, Q = mp.eighe((M + M.transpose_conj()) * mp.mpf('0.5'))
            return mp_re_inner(Gn, Q * mp.diag([mp.expj(-tt * e) for e in ev]) * Q.transpose_conj())
        h = mp.mpf('1e-30')
        d_mp = float((val(h) - val(-h)) / (2 * h))
    assert abs(d_lib - d_mp) <= 1e-7 * max(t, abs(d_mp)), (d_lib, d_mp)


def test_divided_differences_accuracy_vs_mpmath():
    """functional.divided_differences: the Loewner matrix of exp is accurate to 1e-10 for gaps 1e-2 down to 1e-14."""
    mp = _mp()
    w = torch.tensor([1.0, 1.0 + 1e-2, 1.0 + 1e-5, 1.0 + 1e-8, 1.0 + 1e-11, 1.0 + 1e-14, -2.0], dtype=torch.float64)
    L = fc.divided_differences(w, torch.exp(w), torch.exp(w))
    with mp.workdps(40):
        for i in range(len(w)):
            for j in range(len(w)):
                a, b = mp.mpf(float(w[i])), mp.mpf(float(w[j]))
                ref = mp.exp(a) if a == b else (mp.exp(a) - mp.exp(b)) / (a - b)
                assert abs(float(L[i, j]) - float(ref)) <= 1e-10 * abs(float(ref)), (i, j)


# ======================================================================
# 3. rank-deficient inputs
# ======================================================================
def _low_rank(seed, k=4, r=2):
    return cplx((k, r), seed), cplx((k, r), seed + 1)


def test_abs_pinv_polar_constant_rank_gradcheck(alg4):
    """|AB^*|, (AB^*)^+ and the polar factor |AB^*| along the constant-rank manifold (A, B of full column rank)."""
    A, B = (t.requires_grad_(True) for t in _low_rank(21))
    op = lambda a, b: alg4.operator((a @ b.mH)[None, None])
    assert gc(lambda a, b: op(a, b).abs().matrix, A, B)
    assert gc(lambda a, b: op(a, b).inverse().matrix, A, B)
    assert gc(lambda a, b: op(a, b).polar()[1].matrix, A, B)


def test_self_adjoint_pinv_and_powers_constant_rank_gradcheck(alg4):
    """A = Y D Y^* of rank 2 (indefinite D) and A = Y Y^*: A^+, |A|, A^{-1/2}, A^{3/2} along the constant-rank
    manifold (gradcheck in Y)."""
    Y = cplx((4, 2), 23).requires_grad_(True)
    D = torch.diag(torch.tensor([1.5, -0.7], dtype=CD))
    op = lambda m, **f: alg4.operator(m[None, None], is_self_adjoint=True, **f)
    assert gc(lambda y: op(y @ D @ y.mH).inverse().matrix, Y)
    assert gc(lambda y: op(y @ D @ y.mH).abs().matrix, Y)
    assert gc(lambda y: op(y @ y.mH, is_positive=True).power(-0.5).matrix, Y)
    assert gc(lambda y: op(y @ y.mH, is_positive=True).power(1.5).matrix, Y)


@pytest.mark.parametrize("kind", ['generic', 'kernel_only'])
def test_abs_symmetric_derivative_at_singular_vs_mpmath(alg4, kind):
    """At singular X, <grad |X|, E> is the symmetric derivative lim (|X+tE| - |X-tE|)/2t (which exists in every
    direction although the one-sided ones differ); for E in the kernel block it is 0 while the one-sided
    derivative is |E~_22| (checked with mpmath)."""
    mp = _mp()
    A, B = _low_rank(25)
    X0 = (A @ B.mH)[None, None]
    U, s, Vh = torch.linalg.svd(X0[0, 0])
    W0, V0 = U[:, 2:], Vh.mH[:, 2:]                      # kernel of X^* and of X
    E = cplx(X0.shape, 26)
    if kind == 'kernel_only':
        E = (W0 @ cplx((2, 2), 27) @ V0.mH)[None, None]
    G = herm(cplx(X0.shape, 28))
    d_lib = float(directional(lambda X: (alg4.operator(X).abs().matrix * G.conj()).real.sum(), X0, E))

    with mp.workdps(60):
        Am, Bm = mp_matrix(A.numpy()), mp_matrix(B.numpy())
        Xm, Em, Gn = Am * Bm.transpose_conj(), mp_matrix(E[0, 0].numpy()), G[0, 0].numpy()

        def absval(t):
            Y = Xm + Em * t
            return mp_re_inner(Gn, mp_hfun(Y.transpose_conj() * Y, lambda e: mp.sqrt(max(e, 0))))
        h = mp.mpf('1e-25')
        sym = float((absval(h) - absval(-h)) / (2 * h))
        right = float((absval(h) - absval(0)) / h)
    assert abs(d_lib - sym) <= 1e-8 * max(1.0, abs(sym)), (d_lib, sym)
    if kind == 'kernel_only':
        assert abs(d_lib) <= 1e-10
        E22 = W0.mH @ E[0, 0] @ V0                       # one-sided derivative: <G, V0 |E22| V0^*>
        ev, Q = torch.linalg.eigh(E22.mH @ E22)
        absE22 = (Q * ev.clamp(min=0).sqrt().to(CD)) @ Q.mH
        one_sided = float((G[0, 0].conj() * (V0 @ absE22 @ V0.mH)).real.sum())
        assert abs(right - one_sided) <= 1e-6 * max(1.0, abs(one_sided))
        assert abs(one_sided) > 1e-3                     # the one-sided derivatives really differ


@pytest.mark.parametrize("name", ['sqrt', 'power0.5', 'entropy_vn', 'entropy_op', 'fidelity'])
def test_rank_preserving_derivatives_vs_mpmath(alg4, name):
    """sqrt, A^{1/2}, entropies and the fidelity at a rank-2 state rho = Y Y^* / Tr: the derivative along the
    rank-preserving direction dY is exact (mpmath; the forward value of a sqrt carries O(eps^{1/2}) noise from the
    rounding-level zero eigenvalues, so gradcheck does not apply)."""
    mp = _mp()
    Y = cplx((4, 2), 30)
    dY = cplx((4, 2), 31)
    sig = st.random_density_matrix(alg4)
    S = sig.matrix[0, 0].numpy()
    G = herm(cplx((4, 4), 32))

    def rho_of(y):
        r = y @ y.mH
        return r / torch.diagonal(r).real.sum()

    def op(m):
        return alg4.operator(m[None, None], is_self_adjoint=True, is_positive=True)
    lib = {
        'sqrt': lambda m: (op(m).sqrt().matrix[0, 0] * G.conj()).real.sum(),
        'power0.5': lambda m: (op(m).power(0.5).matrix[0, 0] * G.conj()).real.sum(),
        'entropy_vn': lambda m: st.von_neumann_entropy(op(m))[0],
        'entropy_op': lambda m: op(m).entropy()[0],
        'fidelity': lambda m: st.fidelity(op(m), sig)[0],
    }[name]
    t = torch.zeros((), dtype=torch.float64, requires_grad=True)
    d_lib = float(torch.autograd.grad(lib(rho_of(Y + t * dY)), t)[0])

    with mp.workdps(80):
        Ym, dYm, Gn = mp_matrix(Y.numpy()), mp_matrix(dY.numpy()), G.numpy()
        Sm = mp_matrix(S)

        def rho_mp(tt):
            y = Ym + dYm * tt
            r = y * y.transpose_conj()
            tr = sum(mp.re(r[i, i]) for i in range(r.rows))
            return r / tr

        def F(tt):
            r = rho_mp(tt)
            if name in ('sqrt', 'power0.5'):
                return mp_re_inner(Gn, mp_hfun(r, lambda e: mp.sqrt(max(e, 0))))
            if name.startswith('entropy'):
                E, _ = mp.eighe((r + r.transpose_conj()) / 2)
                return -sum(e * mp.log(e) for e in E if e > mp.mpf('1e-40'))
            rs = mp_hfun(r, lambda e: mp.sqrt(max(e, 0)))
            M = rs * Sm * rs
            E, _ = mp.eighe((M + M.transpose_conj()) / 2)
            return sum(mp.sqrt(max(e, 0)) for e in E) ** 2
        h = mp.mpf('1e-30')
        d_mp = float((F(h) - F(-h)) / (2 * h))
    assert abs(d_lib - d_mp) <= 1e-6 * max(1.0, abs(d_mp)), (d_lib, d_mp)


@pytest.mark.parametrize("name", ['sqrt', 'log', 'power0.5', 'power-0.5', 'abs', 'inverse', 'polar_P', 'entropy_vn',
                                  'entropy_op', 'fidelity', 'trace_distance', 'relative_entropy', 'petz2',
                                  'sandwiched0.5', 'kosaki', 'max_rel', 'lp1', 'trace_norm', 'contrast'])
@pytest.mark.parametrize("zeros", ['exact', 'rounding'])
def test_gradients_finite_at_singular(alg4, name, zeros):
    """At singular inputs (exact zero eigenvalues of a diagonal state, or rounding-level ones of Y Y^*) every
    gradient is finite (the derivative itself may not exist in the kernel directions)."""
    if zeros == 'exact':
        x0 = torch.diag(torch.tensor([0.6, 0.4, 0.0, 0.0], dtype=CD))[None, None]
    else:
        Y = cplx((4, 2), 33)
        r = Y @ Y.mH
        x0 = (r / torch.diagonal(r).real.sum())[None, None]
    sig = st.random_density_matrix(alg4)
    x = x0.clone().requires_grad_(True)
    A = alg4.operator(x, is_self_adjoint=True, is_positive=True)
    out = {
        'sqrt': lambda: A.sqrt().matrix, 'log': lambda: A.log().matrix, 'power0.5': lambda: A.power(0.5).matrix,
        'power-0.5': lambda: A.power(-0.5).matrix, 'abs': lambda: alg4.operator(x).abs().matrix,
        'inverse': lambda: alg4.operator(x).inverse().matrix, 'polar_P': lambda: alg4.operator(x).polar()[1].matrix,
        'entropy_vn': lambda: st.von_neumann_entropy(A), 'entropy_op': lambda: A.entropy(),
        'fidelity': lambda: st.fidelity(A, sig), 'trace_distance': lambda: st.trace_distance(A, sig),
        'relative_entropy': lambda: st.relative_entropy(A, sig),
        'petz2': lambda: dv.petz_renyi_divergence(A, sig, 2.0),
        'sandwiched0.5': lambda: dv.sandwiched_renyi_divergence(A, sig, 0.5),
        'kosaki': lambda: dv.kosaki_norm(sig, A, 2.0), 'max_rel': lambda: dv.max_relative_entropy(A, sig),
        'lp1': lambda: dv.lp_norm(A, 1.0), 'trace_norm': lambda: A.trace_norm(),
        'contrast': lambda: A.michelson_contrast,
    }[name]()
    G = cplx(out.shape, 34) if torch.is_complex(out) else torch.ones_like(out)
    g, = torch.autograd.grad((out * G.conj()).real.sum(), x)
    assert torch.isfinite(g).all()


@pytest.mark.parametrize("zeros", [
    'rounding',
    pytest.param('exact', marks=pytest.mark.xfail(strict=True, reason=(
        "divergences._terms (Umegaki branch) evaluates where(lam > thr, xlogy(lam, lam), 0) and differentiates it "
        "in functional.values_and_derivative: at an exact zero eigenvalue the backward of xlogy w.r.t. its second "
        "argument is 0 * lam / lam = NaN, which enters the Loewner matrix and makes the whole gradient NaN "
        "(also for alpha = 1 entries of a per-sample alpha tensor)"))),
])
def test_umegaki_gradient_finite_at_singular_rho(alg4, zeros):
    """D(rho||sigma) with faithful sigma and singular rho (pure, diagonal or rank-deficient): finite gradient in rho."""
    if zeros == 'exact':
        x0 = torch.diag(torch.tensor([0.6, 0.4, 0.0, 0.0], dtype=CD))[None, None]
    else:
        Y = cplx((4, 2), 35)
        r = Y @ Y.mH
        x0 = (r / torch.diagonal(r).real.sum())[None, None]
    sig = st.random_density_matrix(alg4)
    x = x0.clone().requires_grad_(True)
    D = dv.umegaki_relative_entropy(alg4.operator(x, is_self_adjoint=True, is_positive=True), sig)
    g, = torch.autograd.grad(D.sum(), x)
    assert torch.isfinite(g).all()


def test_entropy_kernel_direction_conventions(alg4):
    """In the kernel direction of a singular state the entropy has infinite one-sided slope; the two entropy
    functions return finite values by different conventions (von_neumann_entropy: the slope of -w log max(w, 1e-12)
    at the computed eigenvalue, 0 or -log(1e-12) depending on the sign of its rounding; Operator.entropy: 0), and
    agree with the exact derivative along a direction inside the support."""
    x0 = torch.diag(torch.tensor([0.6, 0.4, 0.0, 0.0], dtype=CD))[None, None]
    E_ker = torch.zeros_like(x0)
    E_ker[0, 0, 3, 3] = 1.0
    E_rng = torch.zeros_like(x0)
    E_rng[0, 0, 0, 0], E_rng[0, 0, 1, 1], E_rng[0, 0, 0, 1], E_rng[0, 0, 1, 0] = 1.0, -1.0, 0.3, 0.3
    op = lambda m: alg4.operator(m, is_self_adjoint=True, is_positive=True)
    d_vn = lambda E: float(directional(lambda m: st.von_neumann_entropy(op(m))[0], x0, E))
    d_op = lambda E: float(directional(lambda m: op(m).entropy()[0], x0, E))
    exact = -(math.log(0.6) + 1) + (math.log(0.4) + 1)            # d/dt of -sum x log x along E_rng
    assert abs(d_vn(E_rng) - exact) < 1e-10 and abs(d_op(E_rng) - exact) < 1e-10
    assert math.isfinite(d_vn(E_ker)) and math.isfinite(d_op(E_ker))
    assert abs(d_op(E_ker)) < 1e-12
    # von_neumann_entropy: the slope of -w log max(w, 1e-12) at the computed eigenvalue, 0 for w <= 0 and
    # -log(1e-12) - ... for 0 < w < 1e-12, i.e. it depends on the sign of a rounding-level eigenvalue
    assert d_vn(E_ker) == 0.0 or d_vn(E_ker) == pytest.approx(-math.log(1e-12), rel=1e-6)


# ======================================================================
# 4. states
# ======================================================================
def _state(alg, x):
    r = x @ x.mH
    return alg.operator(r / torch.diagonal(r, dim1=-2, dim2=-1).sum((-2, -1)).real[:, None, None, None],
                        is_self_adjoint=True, is_positive=True)


@pytest.mark.parametrize("name", ['entropy', 'relative_entropy', 'fidelity', 'trace_distance', 'purity', 'expectation',
                                  'born', 'lueders', 'sector_probs', 'mix', 'density'])
def test_state_functionals_gradcheck(alg, name):
    """Functionals of states (entropy, relative entropy, fidelity, ...) w.r.t. the state (gradcheck)."""
    sig = st.random_density_matrix(alg)
    Aobs = hop(alg, cplx((1, 3, 3, 3), 36))
    P = alg.operator(torch.diag_embed(torch.tensor([[1.0, 0.0, 1.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
                                                   dtype=CD))[None], is_self_adjoint=True, is_projection=True)
    x = (cplx((1, 3, 3, 3), 37) * mask_of(alg) + eye_of(alg)).requires_grad_(True)
    F = {
        'entropy': lambda r: st.von_neumann_entropy(r),
        'relative_entropy': lambda r: st.relative_entropy(r, sig) + st.relative_entropy(sig, r),
        'fidelity': lambda r: st.fidelity(r, sig),
        'trace_distance': lambda r: st.trace_distance(r, sig),
        'purity': lambda r: st.purity(r),
        'expectation': lambda r: st.expectation(r, Aobs),
        'born': lambda r: st.born_probabilities(r, [P, alg.identity() - P]),
        'lueders': lambda r: st.lueders_update(r, P)[0].matrix,
        'sector_probs': lambda r: st.sector_probabilities(r),
        'mix': lambda r: vn.DensityMatrix(alg, matrix=r.matrix, validate=False).mix(sig, 0.3).matrix,
        'density': lambda r: vn.DensityMatrix(alg, matrix=r.matrix, validate=False).density('tau_vN'),
    }[name]
    assert gc(lambda m: F(_state(alg, m)), x)


def test_mix_weight_and_from_density_gradcheck(alg):
    """DensityMatrix.mix with a tensor weight p (gradient in p) and DensityMatrix.from_density (normalize=True)."""
    r1, r2 = st.random_density_matrix(alg), st.random_density_matrix(alg)
    p = torch.tensor([0.3], dtype=torch.float64, requires_grad=True)
    assert gc(lambda q: r1.mix(r2, q).matrix, p)
    x = (cplx((1, 3, 3, 3), 38) * mask_of(alg) + 2 * eye_of(alg)).requires_grad_(True)
    assert gc(lambda m: vn.DensityMatrix.from_density(alg, m @ m.mH, 'norm', normalize=True).matrix, x)


def test_partition_function_gradcheck(alg):
    """log Z and Z w.r.t. H (d log Z / dH = rho_beta, gradcheck)."""
    x = spectral_matrix(alg, [[-1.0, 0.2, 0.9], [-0.3, 1.4], [0.5]], seed=39).requires_grad_(True)
    assert gc(lambda m: st.log_partition_function(hop(alg, m), 1.3), x)
    assert gc(lambda m: st.partition_function(hop(alg, m), 0.7), x)
    g, = torch.autograd.grad(st.log_partition_function(hop(alg, x), 1.3).sum(), x)
    rho = st.gibbs_state(hop(alg, x.detach()), 1.3).matrix
    torch.testing.assert_close(g, -1.3 * rho, rtol=1e-10, atol=1e-12)          # Hellmann-Feynman: grad = -beta rho


@pytest.mark.xfail(strict=True, reason=(
    "states.gibbs_state builds exp(-beta H) from H._decompositions('hermitian'), which runs under torch.no_grad: "
    "the state carries no gradient in H (silently detached), unlike log_partition_function and apply_function"))
def test_gibbs_state_gradcheck(alg):
    """rho_beta = exp(-beta H) / Z is differentiable in H (gradcheck)."""
    x = spectral_matrix(alg, [[-1.0, 0.2, 0.9], [-0.3, 1.4], [0.5]], seed=40).requires_grad_(True)
    assert gc(lambda m: st.gibbs_state(hop(alg, m), 1.1).matrix, x)


# ======================================================================
# 5. dynamics
# ======================================================================
@pytest.mark.parametrize("gap", [None] + GAPS)
def test_propagator_gradcheck(alg, gap):
    """U(t) = exp(-itH) w.r.t. H, also at nearly degenerate H (gradcheck)."""
    g = 0.3 if gap is None else gap
    x = spectral_matrix(alg, [[-1.0, 0.5, 0.5 + g], [1.1, 1.1 + g], [0.2]], seed=41).requires_grad_(True)
    assert gc(lambda m: dyn.propagator(hop(alg, m), 0.9).matrix, x)


def test_propagator_time_gradient(alg):
    """dU/dt = -i H U with t a tensor (gradcheck in t)."""
    H = hop(alg, spectral_matrix(alg, [[-1.0, 0.5, 2.0], [1.1, -0.4], [0.2]], seed=42))
    t = torch.tensor(0.7, dtype=torch.float64, requires_grad=True)
    assert gc(lambda s: dyn.propagator(H, s).matrix, t)


def test_schrodinger_and_von_neumann_gradient_in_initial_state(alg):
    """schrodinger and von_neumann are differentiable in psi0 / rho0 (gradcheck)."""
    H = hop(alg, spectral_matrix(alg, [[-1.0, 0.5, 2.0], [1.1, -0.4], [0.2]], seed=43))
    psi = (cplx((1, 3, 3), 44) * rows_of(alg)).requires_grad_(True)
    assert gc(lambda p: dyn.schrodinger(p, H, [0.0, 0.4, 1.0]), psi)
    x = (cplx((1, 3, 3, 3), 45) * mask_of(alg)).requires_grad_(True)
    assert gc(lambda m: dyn.von_neumann(_state(alg, m), H, [0.0, 0.6])[1].matrix, x)


@pytest.mark.xfail(strict=True, reason=(
    "dynamics.schrodinger and dynamics.von_neumann take the eigendecomposition from H._decompositions('hermitian') "
    "(computed under torch.no_grad): the trajectories carry no gradient in H, silently (propagator does)"))
@pytest.mark.parametrize("fn", ['schrodinger', 'von_neumann'])
def test_schrodinger_and_von_neumann_gradient_in_H(alg, fn):
    """psi(t) and rho(t) are differentiable in the Hamiltonian (gradcheck)."""
    x = spectral_matrix(alg, [[-1.0, 0.5, 2.0], [1.1, -0.4], [0.2]], seed=46).requires_grad_(True)
    psi = cplx((1, 3, 3), 47) * rows_of(alg)
    rho0 = st.random_density_matrix(alg)
    if fn == 'schrodinger':
        assert gc(lambda m: dyn.schrodinger(psi, hop(alg, m), [0.0, 0.5]), x)
    else:
        assert gc(lambda m: dyn.von_neumann(rho0, hop(alg, m), [0.0, 0.5])[1].matrix, x)


def test_rk4_solvers_gradcheck(alg):
    """rk4-based solvers (schrodinger_rk4, lindblad_evolve, solve_operator_ode) w.r.t. H, jumps, initial data."""
    x = spectral_matrix(alg, [[-1.0, 0.5, 2.0], [1.1, -0.4], [0.2]], seed=48).requires_grad_(True)
    L = (cplx((1, 3, 3, 3), 49, 0.5) * mask_of(alg)).requires_grad_(True)
    psi = cplx((1, 3, 3), 50) * rows_of(alg)
    rho0 = st.random_density_matrix(alg)
    assert gc(lambda m: dyn.schrodinger_rk4(psi, hop(alg, m), [0.0, 0.3], substeps=3), x)
    assert gc(lambda m, l: dyn.lindblad_evolve(rho0, hop(alg, m), [alg.operator(l)], times=[0.0, 0.3],
                                               substeps=3)[-1].matrix, x, L)
    A = alg.operator(cplx((1, 3, 3, 3), 51) * mask_of(alg))
    assert gc(lambda m: dyn.solve_operator_ode(lambda t, X: A.matrix @ X - X @ A.matrix, alg.operator(m),
                                               [0.0, 0.2], substeps=2)[-1].matrix, x)


def test_lindblad_superoperator_and_channel_gradcheck(alg):
    """lindblad_superoperator (linear) and lindblad_channel = exp(tL) through the Choi eigendecomposition
    (gradcheck w.r.t. H and the jump operator; the channel is compared through its action on a state)."""
    x = spectral_matrix(alg, [[-1.0, 0.5, 2.0], [1.1, -0.4], [0.2]], seed=52).requires_grad_(True)
    L = (cplx((1, 3, 3, 3), 53, 0.5) * mask_of(alg)).requires_grad_(True)
    rho0 = st.random_density_matrix(alg)
    assert gc(lambda m, l: dyn.lindblad_superoperator(alg, hop(alg, m), [alg.operator(l)])[0], x, L)
    assert gc(lambda m, l: dyn.lindblad_channel(alg, hop(alg, m), [alg.operator(l)], t=0.4).apply(rho0).matrix, x, L)


# ======================================================================
# 6. channels, composite systems, inclusions
# ======================================================================
@pytest.mark.filterwarnings("ignore:.*not trace preserving")
def test_channel_maps_gradcheck(alg):
    """Channel.apply / adjoint / compose / choi / superoperator / from_superoperator w.r.t. Kraus operators and
    states; InterSectorChannel.apply, adjoint, transition_matrix (gradcheck)."""
    K = (cplx((1, 2, 3, 3, 3), 54, 0.5) * mask_of(alg)).requires_grad_(True)
    x = (cplx((1, 3, 3, 3), 55) * mask_of(alg)).requires_grad_(True)
    rho = st.random_density_matrix(alg)
    assert gc(lambda k, m: ch.Channel(alg, k).apply(alg.operator(m)).matrix, K, x)
    assert gc(lambda k: ch.Channel(alg, k).adjoint().apply(rho).matrix, K)
    assert gc(lambda k: ch.Channel(alg, k).compose(ch.Channel(alg, k)).apply(rho).matrix, K)
    assert gc(lambda k: ch.Channel(alg, k).choi()[0], K)
    assert gc(lambda k: ch.Channel(alg, k).superoperator()[1], K)
    assert gc(lambda k: ch.Channel.from_superoperator(alg, ch.Channel(alg, k).superoperator()).apply(rho).matrix, K)
    out = vn.TypeIAlgebra([2, 2], [2, 2], precision='double')
    Ki = (cplx((1, 2, 2, 3, 2, 3), 56, 0.5)).requires_grad_(True)
    assert gc(lambda k: vn.InterSectorChannel(alg, out, k).apply(rho).matrix, Ki)
    assert gc(lambda k: vn.InterSectorChannel(alg, out, k).adjoint('norm').apply(out.identity()).matrix, Ki)
    assert gc(lambda k: vn.InterSectorChannel(alg, out, k).transition_matrix(rho), Ki)


def test_composite_kron_partial_trace_gradcheck():
    """kron and partial_trace in plain and fused tensor products (gradcheck)."""
    a1 = vn.TypeIAlgebra([1, 2], [1, 2], precision='double')
    a2 = vn.TypeIAlgebra([2, 1], [2, 1], precision='double')
    for prod in (composite.tensor_product, composite.fused_tensor_product):
        a12 = prod(a1, a2)
        x = cplx((1, 2, 2, 2), 57).requires_grad_(True)
        y = cplx((1, 2, 2, 2), 58).requires_grad_(True)
        assert gc(lambda m, n: composite.kron(a1.operator(m), a2.operator(n), a12).matrix, x, y)
        z = cplx((1, a12.C, a12.k_max, a12.k_max), 59).requires_grad_(True)
        for keep in (1, 2):
            assert gc(lambda m: composite.partial_trace(a12.operator(m), keep).matrix, z)


def test_inclusion_trace_weights_gradcheck():
    """Conditional expectation, Watatani index and Jones projection w.r.t. the trace weights (gradcheck)."""
    N = vn.TypeIAlgebra([1, 2], [1, 2], precision='double')
    M = vn.TypeIAlgebra([3, 4], [3, 4], precision='double')
    I = vn.Inclusion(N, M, [[1, 1], [2, 1]])
    y = M.operator(cplx((1, 2, 4, 4), 60))
    w = torch.tensor([0.7, 1.3], dtype=torch.float64, requires_grad=True)
    assert gc(lambda ww: I.conditional_expectation(ww).apply(y).matrix, w)
    assert gc(lambda ww: I.conditional_expectation(ww).index().matrix, w)
    assert gc(lambda ww: I.conditional_expectation(ww).pimsner_popa_constant(), w)
    assert gc(lambda ww: I.jones_projection(ww).matrix, w)


# ======================================================================
# 7. divergences
# ======================================================================
DIVS = {
    'petz0.5': lambda r, s: dv.petz_renyi_divergence(r, s, 0.5),
    'petz2': lambda r, s: dv.petz_renyi_divergence(r, s, 2.0),
    'sandwiched0.7': lambda r, s: dv.sandwiched_renyi_divergence(r, s, 0.7),
    'sandwiched3': lambda r, s: dv.sandwiched_renyi_divergence(r, s, 3.0),
    'umegaki': lambda r, s: dv.umegaki_relative_entropy(r, s),
    'max': lambda r, s: dv.max_relative_entropy(r, s),
    'min_sigma': lambda r, s: dv.min_relative_entropy(r, s),
    'renyi_petz_tensor_alpha': lambda r, s: dv.renyi_divergence(
        r, s, torch.tensor([0.5, 1.0, 2.0], dtype=torch.float64), kind='petz'),
    'sector_terms': lambda r, s: dv.sector_contributions(r, s, 1.5, kind='sandwiched').terms,
}


@pytest.mark.parametrize("gap", [None, 1e-9])
@pytest.mark.parametrize("name", list(DIVS))
def test_divergences_gradcheck(alg_s, name, gap):
    """Renyi / Umegaki / max / min divergences and sector terms w.r.t. rho and sigma, generic and with nearly
    degenerate spectra of both states (gradcheck)."""
    alg = alg_s
    g = 0.2 if gap is None else gap
    sp_r = [[0.30, 0.30 * (1 + g)], [0.4 - 0.3 * g]]
    sp_s = [[0.25, 0.25 * (1 + g)], [0.1]]            # D_max: the largest ratio (sector 1) is simple
    xr = spectral_matrix(alg, sp_r, seed=61).requires_grad_(True)
    xs = spectral_matrix(alg, sp_s, seed=62).requires_grad_(True)
    if name == 'renyi_petz_tensor_alpha':
        xr = xr.detach().expand(3, -1, -1, -1).clone().requires_grad_(True)
    op = lambda m: alg.operator(herm(m), is_self_adjoint=True, is_positive=True)
    if name == 'min_sigma':                         # D_0 is differentiable in sigma only
        assert gc(lambda s: DIVS[name](op(xr.detach()), op(s)), xs)
    else:
        assert gc(lambda r, s: DIVS[name](op(r), op(s)), xr, xs)


@pytest.mark.parametrize("kind", ['petz', 'sandwiched'])
def test_renyi_alpha_derivative_vs_mpmath(alg4, kind):
    """dD_alpha/dalpha (alpha a tensor with gradient) agrees with an mpmath derivative to 1e-8 (relative)."""
    mp = _mp()
    xr = spectral_matrix(alg4, [[0.1, 0.2, 0.3, 0.4]], seed=63)
    xs = spectral_matrix(alg4, [[0.4, 0.3, 0.2, 0.1]], seed=64)
    rho = alg4.operator(xr, is_self_adjoint=True, is_positive=True)
    sig = alg4.operator(xs, is_self_adjoint=True, is_positive=True)
    a = torch.tensor([1.7], dtype=torch.float64, requires_grad=True)
    D = dv.renyi_divergence(rho, sig, a, kind=kind)
    d_lib = float(torch.autograd.grad(D.sum(), a)[0])
    with mp.workdps(50):
        R, S = mp_matrix(xr[0, 0].numpy()), mp_matrix(xs[0, 0].numpy())

        def Dmp(al):
            if kind == 'petz':
                Q = mp_hfun(R, lambda e: e ** al) * mp_hfun(S, lambda e: e ** (1 - al))
            else:
                Sg = mp_hfun(S, lambda e: e ** ((1 - al) / (2 * al)))
                Q = mp_hfun(Sg * R * Sg, lambda e: e ** al)
            return mp.log(mp.re(sum(Q[i, i] for i in range(4)))) / (al - 1)
        h = mp.mpf('1e-20')
        al0 = mp.mpf(1.7)
        d_mp = float((Dmp(al0 + h) - Dmp(al0 - h)) / (2 * h))
    assert abs(d_lib - d_mp) <= 1e-8 * max(1.0, abs(d_mp)), (d_lib, d_mp)


def test_lp_and_kosaki_norms_gradcheck(alg_s):
    """lp_norm (p = 1, 2.5, inf; flagged self-adjoint and general) and kosaki_norm (eta = 0, 1/2, 1) (gradcheck)."""
    alg = alg_s
    x = (cplx((1, 2, 2, 2), 65) * mask_of(alg)).requires_grad_(True)
    h = spectral_matrix(alg, [[-1.0, 0.4], [0.7]], seed=66).requires_grad_(True)
    rho = spectral_matrix(alg, [[0.3, 0.2], [0.5]], seed=67).requires_grad_(True)
    for p in (1.0, 2.5, math.inf):
        assert gc(lambda m: dv.lp_norm(alg.operator(m), p, 'norm'), x)
        assert gc(lambda m: dv.lp_norm(hop(alg, m), p, 'tau_vN'), h)
    for eta in (0.0, 0.5, 1.0):
        assert gc(lambda m, r: dv.kosaki_norm(alg.operator(m), hop(alg, r, is_positive=True), 3.0, eta), x, rho)


def test_petz_and_universal_recovery_gradcheck(alg_s):
    """Petz recovery (t = 0 and rotated) and the universal recovery map w.r.t. X and sigma (gradcheck)."""
    alg = alg_s
    torch.manual_seed(68)
    chan = ch.random_channel(alg, 2)
    xs = spectral_matrix(alg, [[0.3, 0.2], [0.5]], seed=68).requires_grad_(True)
    X = (cplx((1, 2, 2, 2), 69) * mask_of(alg)).requires_grad_(True)
    op = lambda m: hop(alg, m, is_positive=True)
    assert gc(lambda s, x: dv.petz_recovery(op(s), chan, alg.operator(herm(x))).matrix, xs, X)
    assert gc(lambda s: dv.petz_recovery(op(s), chan, op(X.detach()), t=0.3).matrix, xs)
    assert gc(lambda s: dv.universal_recovery(op(s), chan, op(X.detach()), step=0.5, cutoff=1.5).matrix, xs)


# ======================================================================
# 8. chains, su2, krylov, jordan, optimize
# ======================================================================
def test_chain_fields_hellmann_feynman():
    """d E_0 / d h_i = <psi_0| S^z_i |psi_0> for SpinChain.xxz with tensor fields (eigenvalues() backward), and
    the entanglement entropy of reduced states w.r.t. the state (gradcheck)."""
    chain = SpinChain(4, precision='double')
    h = torch.tensor([0.3, -0.2, 0.5, 0.1], dtype=torch.float64, requires_grad=True)
    H = chain.xxz(1.0, 0.7, h, sector=2)
    E0 = H.eigenvalues()[0][:, 0]
    g, = torch.autograd.grad(E0.sum(), h)
    w, V = torch.linalg.eigh(H.matrix[0, 0].detach())
    psi = V[:, 0]
    ref = torch.stack([(psi.conj() * torch.diagonal(chain.sz(i, sector=2).matrix[0, 0]) * psi).real.sum()
                       for i in range(4)])
    torch.testing.assert_close(g, ref, rtol=1e-10, atol=1e-12)
    small = SpinChain(3, precision='double')
    x = (cplx((1, 4, 3, 3), 70) * mask_of(small.algebra) + eye_of(small.algebra)).requires_grad_(True)
    assert gc(lambda m: small.entanglement_entropy(_state(small.algebra, m), 1), x)
    assert gc(lambda m: small.reduced_state(_state(small.algebra, m), 2).matrix, x)


def test_su2_couplings_gradient():
    """SU2Chain.heisenberg with tensor couplings J, J2, h: lowest energies are differentiable (gradcheck)."""
    chain = SU2Chain(4, precision='double')
    J2 = torch.tensor([0.3], dtype=torch.float64, requires_grad=True)
    h = torch.tensor([0.2], dtype=torch.float64, requires_grad=True)
    assert gc(lambda j2, hh: chain.lowest_energies(chain.heisenberg(1.0, hh, j2)), J2, h)
    Jb = torch.tensor([1.0, 0.8, 1.3], dtype=torch.float64, requires_grad=True)
    assert gc(lambda j: chain.full_spectrum(chain.heisenberg(j, bonds=[(0, 1), (1, 2), (2, 3)]))[:, :3], Jb)


def test_jordan_maps_gradcheck(alg):
    """jordan_product, quadratic_representation, jordan_triple, jordan_power, transpose_map (gradcheck)."""
    x = cplx((1, 3, 3, 3), 72).requires_grad_(True)
    y = cplx((1, 3, 3, 3), 73).requires_grad_(True)
    op = alg.operator
    assert gc(lambda a, b: jordan.jordan_product(op(a), op(b)).matrix, x, y)
    assert gc(lambda a, b: jordan.quadratic_representation(op(a), op(b)).matrix, x, y)
    assert gc(lambda a, b: jordan.jordan_triple(op(a), op(b), op(a)).matrix, x, y)
    assert gc(lambda a: jordan.jordan_power(op(a), 3).matrix, x)
    assert gc(lambda a: jordan.transpose_map(op(a)).matrix, x)


@pytest.mark.parametrize("f", ['exp', 'sqrt', 'log', 'inverse', 'determinant'])
def test_spin_factor_at_degenerate_point_gradcheck(f):
    """SpinFactor functional calculus is smooth also at u = 0 (lambda_+ = lambda_-) (gradcheck there)."""
    V = vn.SpinFactor(3)
    x = torch.tensor([[0.8, 0.0, 0.0, 0.0], [0.9, 0.1, -0.2, 0.05]], dtype=torch.float64, requires_grad=True)
    F = {'exp': V.exp, 'sqrt': V.sqrt, 'log': V.log, 'inverse': V.inverse, 'determinant': V.determinant}[f]
    assert gradcheck(F, (x,), eps=1e-6, atol=1e-6)


def test_optimize_parametrisations_gradcheck():
    """UnitaryParam (U0 exp(G - G^*)), SpectrumParam (with free contrast), PositiveParam, SelfAdjointParam:
    gradcheck of the operators w.r.t. their parameters."""
    torch.manual_seed(74)
    alg = vn.TypeIAlgebra([2, 3], [2, 3], precision='double')
    up = optimize.UnitaryParam(alg, 2)
    G0 = cplx(up.G.shape, 75, 0.3) * mask_of(alg)

    def unitary(g):
        with torch.no_grad():
            up.G.copy_(torch.zeros_like(up.G))
        skew = g - g.mH
        return up.base @ alg.operator(skew).expm().matrix
    assert gc(unitary, G0.clone().requires_grad_(True))
    sp = optimize.SpectrumParam(alg, 2)
    th, eta = (p.detach().clone().requires_grad_(True) for p in sp.parameters())

    def spectrum(t, e):
        sp.theta, sp.eta = t, e
        return sp.eigenvalues()
    assert gc(spectrum, th, eta)
    for cls in (optimize.PositiveParam, optimize.SelfAdjointParam):
        P = cls(alg, 2)
        th2 = P.spectrum.theta.detach().clone().requires_grad_(True)

        def build(t, P=P):
            P.spectrum.theta = t
            return P.operator().matrix
        assert gc(build, th2)


def test_extremize_reaches_known_optimum():
    """extremize over a UnitaryParam maximises Re Tr(A U) to its known value, the trace norm of A (von Neumann)."""
    torch.manual_seed(76)
    alg = vn.TypeIAlgebra([3], [3], precision='double')
    A = alg.operator(cplx((1, 1, 3, 3), 77).expand(4, -1, -1, -1).clone())
    up = optimize.UnitaryParam(alg, 4)
    best = optimize.extremize(lambda: (torch.diagonal(A.matrix @ up.matrix(), dim1=-2, dim2=-1).sum((-2, -1))).real,
                              [up], maximize=True, rounds=6, steps=60, lr=0.1)
    target = A.trace_norm()
    assert torch.all(best <= target + 1e-9) and torch.all(best >= target - 1e-4)
