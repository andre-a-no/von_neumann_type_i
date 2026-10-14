"""Review of the core modules against independent dense numpy references.

Every operator of M = (+)_c M_{k_c} is also realised on the full space H = (+)_c C^{k_c} (x) C^{m_c} as the dense
matrix (+)_c A_c (x) 1_{m_c} (multiplicities m_c > 1 included); the library results are compared with numpy / scipy
computations on these dense matrices. Double precision, CPU, fixed seeds.
"""
import math
import warnings

import numpy as np
import pytest
import scipy.linalg as sla
import torch

import torch_vn_algebra as tv
from torch_vn_algebra import (states as st, channels as ch, dynamics as dy, composite as cp, krylov, optimize,
                              functional as fc, cost)
from torch_vn_algebra.chains import SpinChain
from torch_vn_algebra.hilbert_space import HilbertSpace, tensor_product_hilbert, direct_sum_hilbert

# (k_c, m_c): one sector of size 1, one factor, unequal sizes with multiplicities, equal sizes, very uneven sizes
LAYOUTS = [([1], [1]), ([3], [2]), ([2, 3, 1], [1, 2, 3]), ([3, 3], [2, 1]), ([4, 1, 2], [1, 1, 2])]
IDS = ['k1', 'C1', 'mixed', 'equal', 'uneven']
layouts = pytest.mark.parametrize('ks,ms', LAYOUTS, ids=IDS)
CD = torch.complex128


@pytest.fixture(autouse=True)
def _seed():
    torch.manual_seed(1234)
    np.random.seed(1234)


# ----------------------------------------------------------------------------------------------------------------
# dense references
# ----------------------------------------------------------------------------------------------------------------
def alg_of(ks, complex_=True):
    return tv.TypeIAlgebra(list(ks), list(ks), precision='double', complex_valued=complex_, device='cpu')


def blocks(mat, ks):
    """list over sectors of numpy arrays (B, k_c, k_c)."""
    m = mat.detach().cpu().resolve_conj().resolve_neg().numpy()
    return [m[:, c, :k, :k] for c, k in enumerate(ks)]


def dense(mat, ks, ms):
    """(B, N, N): (+)_c A_c (x) 1_{m_c}."""
    bl = blocks(mat, ks)
    out = []
    for b in range(bl[0].shape[0]):
        out.append(sla.block_diag(*[np.kron(bl[c][b], np.eye(m)) for c, m in enumerate(ms)]))
    return np.stack(out)


def slices(ks, ms):
    off, out = 0, []
    for k, m in zip(ks, ms):
        out.append(slice(off, off + k * m))
        off += k * m
    return out


def undense(D, ks, ms):
    """Inverse of dense: the blocks A_c read off from (B, N, N) (first copy of every multiplicity)."""
    K = max(ks)
    out = np.zeros((D.shape[0], len(ks), K, K), dtype=D.dtype)
    for c, (k, m, sl) in enumerate(zip(ks, ms, slices(ks, ms))):
        out[:, c, :k, :k] = D[:, sl, sl][:, ::m, ::m]
    return out


def sector_sum(vals_per_block, ms):
    """sum_c x_c / m_c for quantities that the dense space counts m_c times."""
    return sum(v / m for v, m in zip(vals_per_block, ms))


def padding_is_zero(mat, ks):
    m = mat.detach()
    for c, k in enumerate(ks):
        if m[:, c, k:, :].abs().max().item() if m[:, c, k:, :].numel() else 0:
            return False
        if m[:, c, :, k:].abs().max().item() if m[:, c, :, k:].numel() else 0:
            return False
    return True


def ginibre(alg, B=2):
    return alg.from_blocks([torch.randn(B, k, k, dtype=alg.dtype) for k in alg.k_factors])


def hermitian(alg, B=2, shift=0.0):
    bl = []
    for k in alg.k_factors:
        x = torch.randn(B, k, k, dtype=alg.dtype)
        bl.append(x + x.mH + shift * torch.eye(k, dtype=alg.dtype))
    return alg.operator(alg.from_blocks(bl).matrix, is_self_adjoint=True)


def positive(alg, B=2, eps=0.1):
    bl = []
    for k in alg.k_factors:
        x = torch.randn(B, k, k, dtype=alg.dtype)
        bl.append(x @ x.mH + eps * torch.eye(k, dtype=alg.dtype))
    return alg.operator(alg.from_blocks(bl).matrix, is_self_adjoint=True, is_positive=True)


def dense_fun(D, f):
    """f(D) for Hermitian D by numpy eigh."""
    w, V = np.linalg.eigh(D)
    return (V * f(w)[..., None, :]) @ V.conj().swapaxes(-1, -2)


def zscore(x, mean):
    x = np.asarray(x, dtype=float)
    return abs(x.mean() - mean) / (x.std(ddof=1) / math.sqrt(len(x)))


# ================================================================================================================
# algebra.py: random unitaries
# ================================================================================================================
@pytest.mark.parametrize('n', [1, 2, 4])
def test_random_unitary_measures_are_unitary(n):
    """random_unitary: every measure gives unitaries; SU(n) has det 1, COE is symmetric, CSE is self-dual."""
    alg = alg_of([2])
    I = torch.eye(n, dtype=CD)
    for m in alg.UNITARY_MEASURES:
        if m == 'cse' and n % 2:
            continue
        U = alg.random_unitary(n, m, batch_size=50)
        assert torch.allclose(U @ U.mH, I.expand_as(U), atol=1e-12), m
        if m == 'haar_su':
            assert torch.allclose(torch.linalg.det(U), torch.ones(50, dtype=CD), atol=1e-12)
        if m == 'coe':
            assert torch.allclose(U, U.transpose(-2, -1), atol=1e-12)
        if m == 'cse':
            h = n // 2
            J = torch.zeros(n, n, dtype=CD)
            J[:h, h:] = torch.eye(h)
            J[h:, :h] = -torch.eye(h)
            assert torch.allclose(J @ U.transpose(-2, -1) @ J.T, U, atol=1e-12)    # S^R = S
    real = alg_of([2], complex_=False)
    O = real.random_unitary(n, 'haar_su', batch_size=50)
    assert O.dtype == torch.float64 and torch.allclose(torch.linalg.det(O), torch.ones(50, dtype=torch.float64))
    with pytest.raises(ValueError):
        real.random_unitary(n, 'coe')


def test_random_unitary_trace_moments():
    """Haar U(3), SU(3), O(3): E|Tr U|^2 = 1; COE(3): 2n/(n+1); CSE(2N): 4N/(2N-1) (z-score < 5)."""
    alg = alg_of([2])
    N = 20000
    cases = [('haar', 3, 1.0), ('haar_su', 3, 1.0), ('coe', 3, 1.5), ('cse', 4, 8 / 3), ('cse', 6, 12 / 5)]
    for m, n, ref in cases:
        t = alg.random_unitary(n, m, batch_size=N).diagonal(dim1=-2, dim2=-1).sum(-1).abs() ** 2
        assert zscore(t, ref) < 5, (m, n)
    cse2 = alg.random_unitary(2, 'cse', batch_size=10).diagonal(dim1=-2, dim2=-1).sum(-1).abs() ** 2
    assert torch.allclose(cse2, torch.full_like(cse2, 4.0))                     # CSE(1) = phase * 1_2
    O = alg_of([2], complex_=False).random_unitary(3, 'haar', batch_size=N)
    assert zscore(O.diagonal(dim1=-2, dim2=-1).sum(-1) ** 2, 1.0) < 5


@pytest.mark.parametrize('complex_', [True, False])
def test_random_isometry(complex_):
    """random_isometry: V^* V = 1 and E|V_11|^2 = 1/n (first column of a Haar unitary)."""
    alg = alg_of([2], complex_)
    V = alg.random_isometry(5, 3, batch_size=20000)
    I = torch.eye(3, dtype=V.dtype)
    assert torch.allclose(V.mH @ V, I.expand(20000, 3, 3), atol=1e-12)
    assert zscore(V[:, 0, 0].abs() ** 2, 1 / 5) < 5
    with pytest.raises(ValueError):
        alg.random_isometry(2, 3)


@layouts
def test_random_unitary_operator(ks, ms):
    """random_unitary_operator: dense (+)_c U_c (x) 1 is unitary, padding zero."""
    alg = alg_of(ks)
    U = alg.random_unitary_operator(3)
    D = dense(U.matrix, ks, ms)
    assert np.allclose(D @ D.conj().swapaxes(-1, -2), np.eye(D.shape[-1]), atol=1e-12)
    assert padding_is_zero(U.matrix, ks) and U.is_normal and U.is_invertible


# ================================================================================================================
# algebra.py: construction and flags
# ================================================================================================================
@layouts
def test_operator_from_eigenvalues_spectrum_and_flags(ks, ms):
    """operator_from_eigenvalues: the dense spectrum is the sampled one with multiplicity m_c; inferred flags."""
    alg = alg_of(ks)
    A = alg.operator_from_eigenvalues(lambda k: torch.linspace(0.5, 2.0, k, dtype=torch.float64), batch_size=2)
    D = dense(A.matrix, ks, ms)
    ref = np.sort(np.concatenate([np.repeat(np.linspace(0.5, 2.0, k), m) for k, m in zip(ks, ms)]))
    assert np.allclose(np.linalg.eigvalsh(D), ref[None], atol=1e-12)
    assert A.is_self_adjoint and A.is_positive and A.is_invertible and not A.is_projection
    lib = torch.cat(A.eigenvalues(), dim=-1).sort(-1).values.numpy()
    assert np.allclose(lib, np.sort(np.concatenate([np.linspace(0.5, 2.0, k) for k in ks]))[None], atol=1e-14)
    P = alg.operator_from_eigenvalues(lambda k: torch.tensor([1.0] + [0.0] * (k - 1), dtype=torch.float64))
    assert P.is_projection and P.is_invertible == (max(ks) == 1)
    Dp = dense(P.matrix, ks, ms)[0]
    assert np.allclose(Dp @ Dp, Dp, atol=1e-12)
    N = alg.operator_from_eigenvalues(lambda k: torch.exp(1j * torch.arange(k, dtype=torch.float64)))
    assert N.is_normal and (N.is_self_adjoint == (max(ks) == 1))


def test_operator_from_eigenvalues_python_list_keeps_double_precision():
    """Regression: a sampler returning a Python list gave float32-rounded eigenvalues in a double algebra."""
    alg = alg_of([2, 3])
    A = alg.operator_from_eigenvalues(lambda k: [1 / 3] * k)
    w = torch.linalg.eigvalsh(A.matrix[0, 1, :3, :3])
    assert (w - 1 / 3).abs().max() < 1e-15
    assert (torch.cat(A.eigenvalues(), -1) - 1 / 3).abs().max() < 1e-15
    B = alg.from_blocks([[[1 / 3, 0.0], [0.0, 1.0]], torch.eye(3, dtype=CD)])
    assert abs(B.matrix[0, 0, 0, 0].real.item() - 1 / 3) < 1e-16


def test_operator_from_eigenvalues_validation():
    """operator_from_eigenvalues: wrong shapes and forced flags that contradict the spectrum raise."""
    alg = alg_of([2, 3])
    with pytest.raises(ValueError):
        alg.operator_from_eigenvalues(lambda k: torch.ones(k + 1))
    with pytest.raises(ValueError):
        alg.operator_from_eigenvalues(lambda k: -torch.ones(k, dtype=torch.float64), force_positive=True)
    with pytest.raises(ValueError):
        alg.operator_from_eigenvalues(lambda k: torch.zeros(k, dtype=torch.float64), force_invertible=True)
    with pytest.raises(ValueError):
        alg_of([2], complex_=False).operator_from_eigenvalues(lambda k: 1j * torch.ones(k, dtype=CD))


@layouts
def test_property_checks_against_dense(ks, ms):
    """is_self_adjoint / is_normal / is_positive / is_invertible / is_projection agree with the dense matrix."""
    alg = alg_of(ks)
    X = ginibre(alg, 1)
    D = dense(X.matrix, ks, ms)[0]
    assert not alg.operator(X.matrix).is_self_adjoint if max(ks) > 1 or len(ks) > 0 else True
    H = alg.operator((X + X.adjoint()).matrix)
    assert H.is_self_adjoint and H.is_normal
    P = alg.operator((X @ X.adjoint()).matrix)
    assert P.is_positive and np.linalg.eigvalsh(D @ D.conj().T).min() > -1e-12
    assert alg.operator(X.matrix).is_invertible == (np.abs(np.linalg.eigvals(D)).min() > 1e-10)
    Q = alg.operator(((-1.0) * P).matrix)
    assert not Q.is_positive
    # a tiny nilpotent is neither self-adjoint nor normal (relative tolerance)
    if max(ks) > 1:
        Nm = torch.zeros(1, len(ks), max(ks), max(ks), dtype=CD)
        c = int(np.argmax(ks))
        Nm[0, c, 0, 1] = 1e-8
        Nop = alg.operator(Nm)
        assert not Nop.is_self_adjoint and not Nop.is_normal and not Nop.is_invertible


def test_validate_flags():
    """alg.operator(..., validate=True) rejects wrong flags."""
    alg = alg_of([2, 3])
    X = ginibre(alg, 1).matrix
    with pytest.raises(ValueError):
        alg.operator(X, validate=True, is_self_adjoint=True)
    with pytest.raises(ValueError):
        alg.operator(-hermitian(alg, 1, shift=20).matrix, validate=True, is_self_adjoint=True, is_positive=True)
    with pytest.raises(ValueError):
        alg.operator(X, validate=True, is_projection=True)
    alg.operator(positive(alg, 1).matrix, validate=True, is_self_adjoint=True, is_positive=True)


@layouts
def test_extreme_eigenvalues_ignore_padding(ks, ms):
    """lambda_max / lambda_min of a negative (resp. positive) definite operator: the zero padding never leaks."""
    alg = alg_of(ks)
    for sign in (1.0, -1.0):
        A = alg.operator((positive(alg, 3, eps=0.5) * sign).matrix, is_self_adjoint=True)
        w = np.linalg.eigvalsh(dense(A.matrix, ks, ms))
        assert np.allclose(A.lambda_max.numpy(), w.max(-1), atol=1e-12)
        assert np.allclose(A.lambda_min.numpy(), w.min(-1), atol=1e-12)
        assert (A.lambda_max * sign > 0).all() and (A.lambda_min * sign > 0).all()


@layouts
def test_traces_against_dense(ks, ms):
    """trace / Tr_blunt = sum_c Tr A_c, Tr_norm = sum_c Tr A_c / k_c, tau_vN = Tr_norm / C, from the dense matrix."""
    alg = alg_of(ks)
    X = ginibre(alg, 3)
    D = dense(X.matrix, ks, ms)
    tr = [np.trace(D[:, sl, sl], axis1=-2, axis2=-1) for sl in slices(ks, ms)]
    blunt = sector_sum(tr, ms)
    norm = sector_sum([t / k for t, k in zip(tr, ks)], ms)
    assert np.allclose(X.trace.numpy(), blunt) and np.allclose(X.Tr_blunt().numpy(), blunt)
    assert np.allclose(X.Tr_norm().numpy(), norm) and np.allclose(X.tau_vN().numpy(), norm / len(ks))
    one = alg.identity(2)
    assert np.allclose(one.tau_vN().numpy(), 1.0) and np.allclose(one.Tr_norm().numpy(), len(ks))


@layouts
def test_norms_against_dense(ks, ms):
    """trace_norm, frobenius_norm, operator_norm of self-adjoint and general operators against the dense matrix."""
    alg = alg_of(ks)
    for X in (ginibre(alg, 3), hermitian(alg, 3)):
        D = dense(X.matrix, ks, ms)
        sv = [np.linalg.svd(D[:, sl, sl], compute_uv=False) for sl in slices(ks, ms)]
        assert np.allclose(X.trace_norm().numpy(), sector_sum([s.sum(-1) for s in sv], ms), atol=1e-12)
        assert np.allclose(X.frobenius_norm().numpy(), np.sqrt(sector_sum([(s ** 2).sum(-1) for s in sv], ms)))
        assert np.allclose(X.operator_norm().numpy(), np.linalg.norm(D, ord=2, axis=(-2, -1)), atol=1e-12)
    # scale invariance of the overflow guard of frobenius_norm
    Y = alg.operator(ginibre(alg, 1).matrix * 1e200)
    assert torch.isfinite(Y.frobenius_norm()).all()


# ================================================================================================================
# algebra.py: inverse and functional calculus
# ================================================================================================================
@layouts
def test_inverse_and_pseudo_inverse(ks, ms):
    """inverse: A^{-1} for invertible A, the Moore-Penrose inverse (pinv of the dense matrix) for singular A."""
    alg = alg_of(ks)
    for X in (ginibre(alg, 2), hermitian(alg, 2)):
        D = dense(X.matrix, ks, ms)
        Xi = X.inverse()
        assert np.allclose(dense(Xi.matrix, ks, ms), np.linalg.inv(D), atol=1e-9)
        assert padding_is_zero(Xi.matrix, ks)
    # rank-deficient blocks (rank k_c - 1, or 0 for k_c = 1)
    bl = []
    for k in ks:
        G = torch.randn(2, k, max(k - 1, 1), dtype=CD) if k > 1 else torch.zeros(2, 1, 1, dtype=CD)
        H = torch.randn(2, max(k - 1, 1), k, dtype=CD) if k > 1 else torch.zeros(2, 1, 1, dtype=CD)
        bl.append(G @ H)
    S = alg.from_blocks(bl)
    D = dense(S.matrix, ks, ms)
    assert np.allclose(dense(S.inverse().matrix, ks, ms), np.linalg.pinv(D, rcond=1e-10), atol=1e-8)
    Sh = alg.operator((S @ S.adjoint()).matrix, is_self_adjoint=True)
    Dh = dense(Sh.matrix, ks, ms)
    assert np.allclose(dense(Sh.inverse().matrix, ks, ms), np.linalg.pinv(Dh, rcond=1e-10, hermitian=True), atol=1e-8)


@layouts
def test_abs_polar_sqrt(ks, ms):
    """abs = (A^*A)^{1/2}, polar A = W|A| with W unitary, sqrt of a positive operator: dense references."""
    alg = alg_of(ks)
    X = ginibre(alg, 2)
    D = dense(X.matrix, ks, ms)
    absD = dense_fun(D.conj().swapaxes(-1, -2) @ D, lambda w: np.sqrt(np.clip(w, 0, None)))
    assert np.allclose(dense(X.abs().matrix, ks, ms), absD, atol=1e-10)
    W, P = X.polar()
    Wd, Pd = dense(W.matrix, ks, ms), dense(P.matrix, ks, ms)
    assert np.allclose(Wd @ Pd, D, atol=1e-10) and np.allclose(Pd, absD, atol=1e-10)
    assert np.allclose(Wd @ Wd.conj().swapaxes(-1, -2), np.eye(D.shape[-1]), atol=1e-10)
    Pos = positive(alg, 2)
    Dp = dense(Pos.matrix, ks, ms)
    S = dense(Pos.sqrt().matrix, ks, ms)
    assert np.allclose(S @ S, Dp, atol=1e-10) and np.allclose(S, dense_fun(Dp, np.sqrt), atol=1e-10)
    H = hermitian(alg, 2)
    assert np.allclose(dense(H.abs().matrix, ks, ms), dense_fun(dense(H.matrix, ks, ms), np.abs), atol=1e-10)
    with pytest.raises(RuntimeError):
        alg.operator((-1.0 * Pos).matrix, is_self_adjoint=True).sqrt()


@layouts
def test_log_power_entropy(ks, ms):
    """log, power (p > 0, p < 0 on the support), trace_a_log_a and entropy against dense eigendecompositions."""
    alg = alg_of(ks)
    A = positive(alg, 2)
    D = dense(A.matrix, ks, ms)
    assert np.allclose(dense(A.log().matrix, ks, ms), dense_fun(D, np.log), atol=1e-10)
    for p in (2.0, 0.5, 0.0, -1.0, -0.5):
        assert np.allclose(dense(A.power(p).matrix, ks, ms), dense_fun(D, lambda w: w ** p), atol=1e-9), p
    tal = [np.sum(np.linalg.eigvalsh(D[:, sl, sl]) * np.log(np.linalg.eigvalsh(D[:, sl, sl])), -1)
           for sl in slices(ks, ms)]
    assert np.allclose(A.trace_a_log_a().numpy(), sector_sum(tal, ms), atol=1e-10)
    assert np.allclose(A.entropy().numpy(), -sector_sum(tal, ms), atol=1e-10)
    # singular positive: negative powers act on the support only (pinv for p = -1)
    S = alg.operator_from_eigenvalues(lambda k: torch.tensor([0.0] + [2.0] * (k - 1), dtype=torch.float64),
                                      batch_size=2)
    Ds = dense(S.matrix, ks, ms)
    assert np.allclose(dense(S.power(-1.0).matrix, ks, ms), np.linalg.pinv(Ds, rcond=1e-10, hermitian=True),
                       atol=1e-9)


@layouts
def test_expm_and_holomorphic_function(ks, ms):
    """expm(scale) and holomorphic_function(exp, sin) of non-normal operators against scipy expm / sinm."""
    alg = alg_of(ks)
    X = ginibre(alg, 2) * 0.5
    D = dense(X.matrix, ks, ms)
    for s in (1.0, -0.7, 0.3j):
        assert np.allclose(dense(X.expm(s).matrix, ks, ms), np.stack([sla.expm(s * d) for d in D]), atol=1e-10)
    E = X.holomorphic_function(torch.exp)
    assert np.allclose(dense(E.matrix, ks, ms), np.stack([sla.expm(d) for d in D]), atol=1e-8)
    assert E.contour_error < 1e-8
    Sn = X.holomorphic_function(torch.sin, atol=1.0)
    assert np.allclose(dense(Sn.matrix, ks, ms), np.stack([sla.sinm(d) for d in D]), atol=1e-8)
    assert padding_is_zero(E.matrix, ks) and padding_is_zero(X.expm().matrix, ks)


@layouts
def test_apply_function_normal_operators(ks, ms):
    """apply_function on unitaries (normal, also with degenerate spectra) against dense functions of D."""
    alg = alg_of(ks)
    U = alg.random_unitary_operator(2)
    D = dense(U.matrix, ks, ms)
    F = U.apply_function(lambda z: z ** 3 + 2 * z)
    assert np.allclose(dense(F.matrix, ks, ms), D @ D @ D + 2 * D, atol=1e-10)
    # degenerate unitary: V diag(1, 1, -1, ...) V^*
    Ud = alg.operator_from_eigenvalues(lambda k: torch.tensor([(-1.0) ** (i // 2) for i in range(k)], dtype=CD)
                                       * torch.exp(0.3j * torch.ones(k, dtype=torch.float64)), batch_size=2)
    Ud = alg.operator(Ud.matrix)            # forget the spectrum: decided numerically as normal
    Dd = dense(Ud.matrix, ks, ms)
    L = Ud.apply_function(torch.log)        # principal log; eigenvalues e^{0.3i}, -e^{0.3i}
    assert np.allclose(dense(L.matrix, ks, ms), np.stack([sla.logm(d) for d in Dd]), atol=1e-9)


def test_apply_function_requires_normal():
    """apply_function of a non-normal operator raises (holomorphic_function is required)."""
    alg = alg_of([3])
    X = ginibre(alg, 1)
    with pytest.raises(RuntimeError):
        X.apply_function(torch.exp).matrix


@layouts
def test_michelson_contrast(ks, ms):
    """michelson_contrast = (l_max - l_min) / (l_max + l_min) over all sectors (dense spectrum); 0 for A = 0."""
    alg = alg_of(ks)
    A = positive(alg, 3)
    w = np.linalg.eigvalsh(dense(A.matrix, ks, ms))
    ref = (w.max(-1) - w.min(-1)) / (w.max(-1) + w.min(-1))
    assert np.allclose(A.michelson_contrast.numpy(), ref)
    assert (alg.zero(2).michelson_contrast == 0).all()
    assert np.allclose(alg.identity(1).michelson_contrast.numpy(), 0.0)


@layouts
def test_central_identity_zero_from_blocks(ks, ms):
    """central(s) = sum_c s_c 1_c (dense: diagonal), identity, zero, from_blocks round trip and flags."""
    alg = alg_of(ks)
    s = [float(c + 1) for c in range(len(ks))]
    Z = alg.central(s)
    ref = np.diag(np.concatenate([np.full(k * m, v) for k, m, v in zip(ks, ms, s)]))
    assert np.allclose(dense(Z.matrix, ks, ms)[0], ref)
    assert Z.is_positive and Z.is_invertible and Z.is_self_adjoint
    Zc = alg.central([1j] * len(ks))
    assert not Zc.is_self_adjoint and Zc.is_normal and Zc.is_invertible
    t = torch.tensor([0.5, 2.0], dtype=torch.float64, requires_grad=True)
    Zt = alg.central([t] * len(ks))
    assert Zt.matrix.shape[0] == 2
    Zt.tau_vN().real.sum().backward()
    assert torch.allclose(t.grad, torch.full((2,), float(len(ks)), dtype=torch.float64) / len(ks))
    assert np.allclose(dense(alg.identity(2).matrix, ks, ms), np.eye(sum(k * m for k, m in zip(ks, ms))))
    assert (alg.zero(2).matrix == 0).all()
    X = ginibre(alg, 2)
    Y = alg.from_blocks([X.matrix[:, c, :k, :k] for c, k in enumerate(ks)])
    assert torch.equal(X.matrix, Y.matrix)
    with pytest.raises(ValueError):
        alg.central([1.0] * (len(ks) + 1))


def test_cast_and_real_algebra():
    """cast: complex -> double / real keeps the matrix, a non-real operator cannot be cast to a real algebra."""
    alg = alg_of([2, 3])
    H = hermitian(alg, 1)
    Hs = H.cast(alg.like(precision='single'))
    assert Hs.matrix.dtype == torch.complex64 and torch.allclose(Hs.matrix.to(CD), H.matrix, atol=1e-6)
    R = alg.operator(H.matrix.real.to(CD))
    Rr = R.cast(alg.like(complex_valued=False))
    assert Rr.matrix.dtype == torch.float64 and torch.equal(Rr.matrix, H.matrix.real)
    with pytest.raises(ValueError):
        (H * 1j).cast(alg.like(complex_valued=False))
    with pytest.raises(ValueError):
        H + Hs                                                     # different precision: explicit cast needed


def test_binary_operations_and_batch_tracking():
    """+, -, @, scalar * and adjoint against dense; regression: the known batch of a lazy result was wrong when one
    operand's batch was unknown (abs of a batch 5 plus a batch-1 identity reported batch 1)."""
    ks, ms = [2, 3], [2, 1]
    alg = alg_of(ks)
    X, Y = ginibre(alg, 5), ginibre(alg, 5)
    Dx, Dy = dense(X.matrix, ks, ms), dense(Y.matrix, ks, ms)
    assert np.allclose(dense((X + Y).matrix, ks, ms), Dx + Dy)
    assert np.allclose(dense((X - Y).matrix, ks, ms), Dx - Dy)
    assert np.allclose(dense((X @ Y).matrix, ks, ms), Dx @ Dy)
    assert np.allclose(dense((X * (2 - 1j)).matrix, ks, ms), (2 - 1j) * Dx)
    assert np.allclose(dense(X.adjoint().matrix, ks, ms), Dx.conj().swapaxes(-1, -2))
    t = torch.arange(5, dtype=torch.float64)
    assert np.allclose(dense((alg.identity(1) * t).matrix, ks, ms), t.numpy()[:, None, None] * np.eye(7))
    Z = X.abs() + alg.identity(1)
    assert Z.shape == Z.matrix.shape == (5, 2, 3, 3)
    W, P = X.polar()
    assert W.shape == (5, 2, 3, 3) and (P @ alg.identity(1)).shape == (5, 2, 3, 3)
    lazy = alg.identity(1) * torch.ones(5, dtype=torch.float64)
    assert lazy.shape == lazy.matrix.shape


# ================================================================================================================
# gradients
# ================================================================================================================
def test_gradcheck_spectral_functions():
    """Daleckii-Krein backward of apply_function (exp, sqrt) and of abs / inverse for general blocks (gradcheck)."""
    alg = alg_of([2, 3])
    H0 = hermitian(alg, 1).matrix.detach()
    P0 = positive(alg, 1).matrix.detach()
    X0 = ginibre(alg, 1).matrix.detach() + 3 * alg.identity(1).matrix

    def f_exp(m):
        return alg.operator(m, is_self_adjoint=True).apply_function(torch.exp).matrix

    def f_sqrt(m):
        return alg.operator(0.5 * (m + m.mH), is_self_adjoint=True, is_positive=True).sqrt().matrix

    def f_abs(m):
        return alg.operator(m, is_self_adjoint=False, is_normal=False).abs().matrix

    def f_inv(m):
        return alg.operator(m, is_self_adjoint=False).inverse().matrix

    for f, x in ((f_exp, H0), (f_sqrt, P0), (f_abs, X0), (f_inv, X0)):
        x = x.clone().requires_grad_(True)
        assert torch.autograd.gradcheck(lambda m: torch.view_as_real(f(m)), (x,), eps=1e-6, atol=1e-5), f.__name__


def test_gradient_of_exp_at_degenerate_zero_eigenvalue():
    """Regression (Daleckii-Krein): eigenvalue 0 of multiplicity 2 computed as +-1e-17; the Loewner entry
    (e^a - e^b)/(a - b) cancelled to 0 instead of 1, so the gradient of Tr(W exp(H)) was wrong."""
    alg = alg_of([3])
    U = alg.random_unitary(3)
    H = U @ torch.diag(torch.tensor([0.0, 0.0, 1.0], dtype=CD)) @ U.mH
    H = 0.5 * (H + H.mH)
    W = torch.randn(3, 3, dtype=CD)
    Hp = H.clone().requires_grad_(True)
    F = alg.operator(Hp[None, None], is_self_adjoint=True).apply_function(torch.exp).matrix[0, 0]
    (W * F).sum().real.backward()
    for _ in range(3):
        E = torch.randn(3, 3, dtype=CD)
        E = E + E.mH
        _, Lf = sla.expm_frechet(H.numpy(), E.numpy())
        ref = (W.numpy() * Lf).sum().real
        assert abs((Hp.grad.conj() * E).sum().real.item() - ref) < 1e-10 * max(1.0, abs(ref))


def test_divided_differences():
    """Loewner matrix: exact quotients for separated points, f' on the diagonal, f' for rounding-level clusters,
    and quotients (not midpoints) for close points of a function with a singularity (log near 0)."""
    w = torch.tensor([-1e-17, 2e-17, 0.5, 1.0], dtype=torch.float64)
    L = fc.divided_differences(w, torch.exp(w), torch.exp(w))
    assert abs(L[0, 1].item() - 1.0) < 1e-12 and abs(L[2, 3].item() - (math.e - math.exp(0.5)) / 0.5) < 1e-14
    assert torch.allclose(L.diagonal(), torch.exp(w))
    v = torch.tensor([1e-8, 2e-8, 1.0], dtype=torch.float64)
    Lg = fc.divided_differences(v, torch.log(v), 1 / v)
    assert abs(Lg[0, 1].item() - math.log(2) / 1e-8) < 1e-6 * math.log(2) / 1e-8
    # exp(-i t w) with t (w_i - w_j) = 2 pi: the quotient (0) and not a midpoint (-i t)
    t = 100.0
    u = torch.tensor([0.0, 2 * math.pi / t], dtype=torch.float64)
    f = lambda x: torch.exp(-1j * t * x)
    Lu = fc.divided_differences(u, f(u), -1j * t * f(u))
    assert abs(Lu[0, 1]) < 1e-12


def test_gradient_of_holomorphic_function_matches_frechet():
    """holomorphic_function is differentiable through the solves: gradient = adjoint Frechet derivative of expm."""
    alg = alg_of([3])
    X0 = 0.5 * ginibre(alg, 1).matrix.detach()
    X = X0.clone().requires_grad_(True)
    W = torch.randn(3, 3, dtype=CD)
    F = alg.operator(X).holomorphic_function(torch.exp).matrix[0, 0]
    (W * F).sum().real.backward()
    E = torch.randn(3, 3, dtype=CD)
    _, Lf = sla.expm_frechet(X0[0, 0].numpy(), E.numpy())
    ref = (W.numpy() * Lf).sum().real
    assert abs((X.grad[0, 0].conj() * E).sum().real.item() - ref) < 1e-8


# ================================================================================================================
# functional.py
# ================================================================================================================
def test_singular_values_svd_and_polar_methods():
    """singular_values / svd / polar: the Hermitian dilation [[0, X], [X^*, 0]] and LAPACK SVD agree with numpy."""
    X = torch.randn(4, 5, 5, dtype=CD)
    ref = np.linalg.svd(X.numpy(), compute_uv=False)
    for method in ('svd', 'eigh'):
        assert np.allclose(fc.singular_values(X, method).numpy(), ref, atol=1e-12)
        U, s, V = fc.svd(X, method)
        assert torch.allclose((U * s.to(CD).unsqueeze(-2)) @ V.mH, X, atol=1e-10)
        W, P = fc.polar(X, method)
        assert torch.allclose(W @ P, X, atol=1e-10)
        assert torch.allclose(W @ W.mH, torch.eye(5, dtype=CD).expand(4, 5, 5), atol=1e-10)
        assert torch.allclose(fc.pinv(X, method=method), torch.linalg.pinv(X), atol=1e-9)


def test_normal_eig_and_rejection_of_non_normal():
    """normal_eig: unitary eigenbasis of a normal stack (also degenerate); a non-normal matrix raises."""
    alg = alg_of([4])
    U = alg.operator_from_eigenvalues(lambda k: torch.tensor([1, 1, 1j, -1], dtype=CD), batch_size=3).matrix[:, 0]
    lam, V = fc.normal_eig(U)
    assert torch.allclose(V.mH @ V, torch.eye(4, dtype=CD).expand(3, 4, 4), atol=1e-12)
    assert torch.allclose((V * lam.unsqueeze(-2)) @ V.mH, U, atol=1e-12)
    with pytest.raises(ValueError):
        fc.normal_eig(torch.tensor([[[1.0, 1.0], [0.0, 1.0]]], dtype=CD))


def test_holomorphic_with_explicit_contour():
    """holomorphic(log) of a non-normal block with spectrum near 1, circle |z - 1| = 0.8: scipy logm."""
    X = torch.eye(3, dtype=CD) + 0.2 * torch.randn(1, 3, 3, dtype=CD)
    F, err = fc.holomorphic(X, torch.log, center=1.0, radius=0.8)
    assert np.allclose(F[0].numpy(), sla.logm(X[0].numpy()), atol=1e-9) and err < 1e-9
    c, r = fc.default_contour(X)
    lam = torch.linalg.eigvals(X)
    assert ((lam - c[:, None]).abs() < r[:, None]).all()
    with pytest.raises(ValueError):
        fc.holomorphic(X, torch.log, center=0.0, radius=2.0)          # branch cut inside the disk


def test_spectral_cache_is_consistent():
    """caching(): functions of a cached operator agree with uncached ones; in-place changes invalidate the cache."""
    alg = alg_of([2, 3])
    H = hermitian(alg, 2)
    a = H.apply_function(torch.exp).matrix
    with tv.caching():
        Hc = alg.operator(H.matrix.clone(), is_self_adjoint=True)
        b = Hc.apply_function(torch.exp).matrix
        c = Hc.apply_function(torch.exp).matrix
        assert Hc.cache_info()['cached']
        Hc.matrix.mul_(2.0)
        d = Hc.apply_function(torch.exp).matrix
    assert torch.allclose(a, b) and torch.allclose(b, c)
    assert torch.allclose(d, (H * 2.0).apply_function(torch.exp).matrix, atol=1e-10)


# ================================================================================================================
# cost.py
# ================================================================================================================
def test_cost_batched_call_and_limits():
    """batched_call returns fn(x) (probed and threaded paths); set_limits validates keys; check_memory raises."""
    x = torch.randn(2000, 6, 6, dtype=torch.float64)                    # large enough for the threaded CPU path
    x = x + x.mT
    old = cost.get_limits()
    try:
        cost.set_limits(probe_flops=0)
        w = cost.batched_call(torch.linalg.eigvalsh, x, "test")
        assert torch.allclose(w, torch.linalg.eigvalsh(x))
        U, s, Vh = cost.batched_call(torch.linalg.svd, x, "test", kind='svd')
        assert torch.allclose((U * s.unsqueeze(-2)) @ Vh, x, atol=1e-12)
        with pytest.raises(KeyError):
            cost.set_limits(nonsense=1)
        with pytest.raises(cost.InsufficientMemoryError):
            cost.check_memory(2.0 ** 70, 'cpu', "test")
        with cost.disabled():
            cost.check_memory(2.0 ** 70, 'cpu', "test")
    finally:
        cost.set_limits(**old)
    assert cost.tensor_bytes((2, 3), torch.complex128) == 96


def test_gram_schmidt_and_device_helpers():
    """GPU code paths run on the CPU: Gram-Schmidt Q (plain and differentiable) equals the QR factor with diag(R) > 0
    (Mezzadri's phase fix), and device comparison treats 'cuda' as the current CUDA device."""
    from torch_vn_algebra.algebra import _orthonormal_columns, _same_device
    Z = torch.randn(7, 5, 3, dtype=CD)
    Q, R = torch.linalg.qr(Z)
    d = torch.diagonal(R, dim1=-2, dim2=-1)
    ref = Q * d.sgn().unsqueeze(-2)
    assert torch.allclose(_orthonormal_columns(Z), ref, atol=1e-12)
    Zg = Z.clone().requires_grad_(True)
    Qg = _orthonormal_columns(Zg)
    assert torch.allclose(Qg, ref, atol=1e-12)
    Qg.abs().sum().backward()
    assert torch.isfinite(Zg.grad).all()
    assert _same_device('cpu', torch.device('cpu')) and not _same_device('cpu', 'cuda')
    assert _same_device('cuda', 'cuda:0') == (not torch.cuda.is_available() or torch.cuda.current_device() == 0)
    assert not _same_device('cuda:1', 'cuda:0')
    U = alg_of([2]).random_unitary(70, batch_size=2).to(torch.complex64)
    lam, V = fc.normal_eig(U)                           # complex64, k > 64: computed in double precision
    assert lam.dtype == torch.complex64 and torch.allclose((V * lam.unsqueeze(-2)) @ V.mH, U, atol=1e-4)


# ================================================================================================================
# hilbert_space.py
# ================================================================================================================
def test_hilbert_space_inner_product_norm_and_bases():
    """inner_product is conjugate-linear in the bra and sums over channels; norm / normalize; bases orthonormal,
    projections of rank k; Haar frames: E|V_11|^2 = 1/n."""
    hs = HilbertSpace(n=4, k=2, batch_size=3, n_channels=2, precision='double', device='cpu')
    x = torch.randn(3, 2, 4, 1, dtype=CD)
    y = torch.randn(3, 2, 4, 1, dtype=CD)
    ip = hs.inner_product(x, y)
    assert ip.shape == (3, 1, 1, 1)
    assert torch.allclose(ip[:, 0, 0, 0], (x.conj() * y).sum((1, 2, 3)))
    assert torch.allclose(hs.inner_product(x.mH, y), ip)                      # a bra row is already conjugated
    assert torch.allclose(hs.norm(x)[:, 0, 0, 0] ** 2, (x.abs() ** 2).sum((1, 2, 3)))
    assert torch.allclose(hs.norm(hs.normalize(x)), torch.ones(3, 1, 1, 1, dtype=torch.float64))
    for basis in (hs.standard_basis(), hs.random_basis(), hs.haar_basis()):
        assert basis.is_orthonormal(1e-12)
        P = basis.projection
        assert torch.allclose(P @ P, P, atol=1e-12)
        assert torch.allclose(P.diagonal(dim1=-2, dim2=-1).sum(-1).real, torch.full((3, 2), 2.0, dtype=torch.float64))
        v = torch.randn(3, 2, 2, 1, dtype=CD)
        assert torch.allclose(basis.restrict_vector(basis.embed_vector(v)), v, atol=1e-12)
        A = torch.randn(3, 2, 2, 2, dtype=CD)
        assert torch.allclose(basis.restrict_operator(basis.embed_operator(A)), A, atol=1e-12)
    one = HilbertSpace(n=1, k=1, precision='double', device='cpu')
    with pytest.raises(ValueError):
        one.inner_product(torch.ones(1, 1, 1, 1, dtype=CD), torch.ones(1, 1, 1, 1, dtype=CD))
    big = HilbertSpace(n=5, k=2, batch_size=20000, precision='double', device='cpu')
    assert zscore(big.random_basis().V[:, 0, 0, 0].abs() ** 2, 1 / 5) < 5
    assert tensor_product_hilbert([hs, big]).n == 20 and direct_sum_hilbert([hs, big]).k == 4
    with pytest.raises(ValueError):
        tensor_product_hilbert([hs, HilbertSpace(2, 1, complex_valued=False, precision='double', device='cpu')])


# ================================================================================================================
# states.py
# ================================================================================================================
def dense_state(rho, ks, ms):
    """Density of the same state on the full space (multiplicities): (+)_c rho_c / m_c (x) 1_{m_c}."""
    scale = torch.tensor([1.0 / m for m in ms], dtype=torch.float64)[None, :, None, None]
    return dense(rho.matrix * scale, ks, ms)


def ent(D):
    w = np.clip(np.linalg.eigvalsh(D), 0, None)
    return -np.sum(np.where(w > 0, w * np.log(np.where(w > 0, w, 1)), 0), -1)


@layouts
def test_random_density_matrix_is_state(ks, ms):
    """random_density_matrix ('hs', 'bures', rank r, real): positive, unit trace, zero padding, rank <= r."""
    for complex_ in (True, False):
        alg = alg_of(ks, complex_)
        for kw in ({}, {'rank': 1}, {'measure': 'bures'} if complex_ else {'rank': 2}):
            rho = st.random_density_matrix(alg, 3, **kw)
            D = dense_state(rho, ks, ms)
            assert np.allclose(np.trace(D, axis1=-2, axis2=-1), 1.0)
            assert np.linalg.eigvalsh(D).min() > -1e-12 and padding_is_zero(rho.matrix, ks)
            if kw.get('rank'):
                for c, k in enumerate(ks):
                    assert np.linalg.matrix_rank(blocks(rho.matrix, ks)[c][0], tol=1e-10) <= kw['rank']


@pytest.mark.parametrize('k', [2, 3])
def test_mean_purity_of_random_states(k):
    """Mean purity: Hilbert-Schmidt E Tr rho^2 = 2k/(k^2+1); Bures (5k^2+1)/(2k(k^2+2)) (z-score < 5)."""
    alg = alg_of([k])
    N = 20000
    assert zscore(st.purity(st.random_density_matrix(alg, N)), 2 * k / (k * k + 1)) < 5
    assert zscore(st.purity(st.random_density_matrix(alg, N, measure='bures')),
                  (5 * k * k + 1) / (2 * k * (k * k + 2))) < 5


@layouts
def test_density_conventions(ks, ms):
    """density(trace) / from_density: omega(A) = trace(d A) for the blunt, normalised and tracial traces."""
    alg = alg_of(ks)
    rho = st.random_density_matrix(alg, 2)
    A = hermitian(alg, 2)
    w = rho.expectation(A)
    for tr in st.TRACES:
        d = rho.density(tr)
        dA = alg.operator(d) @ A
        val = {'blunt': dA.Tr_blunt(), 'norm': dA.Tr_norm(), 'tau_vN': dA.tau_vN()}[tr]
        assert torch.allclose(val, w)
        back = st.DensityMatrix.from_density(alg, d, tr)
        assert torch.allclose(back.matrix, rho.matrix)
    # tracial state: omega = tau_vN; maximally mixed: Tr_blunt / sum k_c
    X = ginibre(alg, 2)
    assert torch.allclose(st.tracial_state(alg).expectation(X), X.tau_vN())
    assert torch.allclose(st.maximally_mixed_state(alg).expectation(X), X.Tr_blunt() / sum(ks))
    with pytest.raises(ValueError):
        st.DensityMatrix(alg, matrix=2 * rho.matrix)
    if sum(ks) > 1:                                    # unit trace, one negative eigenvalue
        with pytest.raises(ValueError):
            st.DensityMatrix(alg, matrix=(2 * st.random_density_matrix(alg, 1, rank=1)
                                          - alg.identity(1) * (1.0 / sum(ks))).matrix)


@layouts
def test_gibbs_state_and_partition_function(ks, ms):
    """gibbs_state = (+)_c exp(-beta H_c) / sum_c Tr exp(-beta H_c) (dense expm); log Z; no overflow at beta = 500."""
    alg = alg_of(ks)
    H = hermitian(alg, 2)
    D = dense(H.matrix, ks, ms)
    for beta in (0.7, -0.3, 500.0):
        G = dense_fun(D, lambda w: np.exp(-beta * (w - w.min(-1, keepdims=True))))
        Z = sector_sum([np.trace(G[:, sl, sl], axis1=-2, axis2=-1) for sl in slices(ks, ms)], ms)
        rho = st.gibbs_state(H, beta)
        assert np.allclose(dense(rho.matrix, ks, ms), G / Z[:, None, None], atol=1e-10)
        assert torch.isfinite(rho.matrix).all()
    w = np.concatenate([np.linalg.eigvalsh(D[:, sl, sl])[:, ::m] for sl, m in zip(slices(ks, ms), ms)], -1)
    assert np.allclose(st.log_partition_function(H, 0.7).numpy(), np.log(np.exp(-0.7 * w).sum(-1)))
    assert np.allclose(st.partition_function(H, 0.7).numpy(), np.exp(-0.7 * w).sum(-1))


@layouts
def test_state_functionals_against_dense(ks, ms):
    """expectation, Born rule, Lueders rule, entropy (+ sum p_c log m_c on the full space), relative entropy,
    fidelity, trace distance and purity against the dense density (+)_c rho_c / m_c (x) 1."""
    alg = alg_of(ks)
    rho, sigma = st.random_density_matrix(alg, 2), st.random_density_matrix(alg, 2)
    Dr, Ds = dense_state(rho, ks, ms), dense_state(sigma, ks, ms)
    A = ginibre(alg, 2)
    Da = dense(A.matrix, ks, ms)
    assert np.allclose(rho.expectation(A).numpy(), np.trace(Dr @ Da, axis1=-2, axis2=-1))
    p = rho.sector_probabilities().numpy()
    assert np.allclose(p, np.stack([np.trace(Dr[:, sl, sl], axis1=-2, axis2=-1).real for sl in slices(ks, ms)], -1))
    S_dense = ent(Dr)
    assert np.allclose(rho.von_neumann_entropy().numpy() + (p * np.log(ms)).sum(-1), S_dense, atol=1e-10)
    logm = lambda D: dense_fun(D, np.log)
    D_ref = np.trace(Dr @ (logm(Dr) - logm(Ds)), axis1=-2, axis2=-1).real
    assert np.allclose(rho.relative_entropy(sigma).numpy(), D_ref, atol=1e-9)
    sq = lambda D: dense_fun(D, lambda w: np.sqrt(np.clip(w, 0, None)))
    F_ref = np.array([np.linalg.svd(sq(Dr)[b] @ sq(Ds)[b], compute_uv=False).sum() ** 2 for b in range(2)])
    assert np.allclose(rho.fidelity(sigma).numpy(), F_ref, atol=1e-10)
    T_ref = 0.5 * np.abs(np.linalg.eigvalsh(Dr - Ds)).sum(-1)
    assert np.allclose(rho.trace_distance(sigma).numpy(), T_ref, atol=1e-12)
    pur = sum(np.trace(Dr[:, sl, sl] @ Dr[:, sl, sl], axis1=-2, axis2=-1).real * m       # Tr rho_c^2 / m_c on H
              for sl, m in zip(slices(ks, ms), ms))
    assert np.allclose(rho.purity().numpy(), pur)
    # projective measurement in the sectors: P = one projection per block (rank 1)
    P = alg.operator_from_eigenvalues(lambda k: torch.tensor([1.0] + [0.0] * (k - 1), dtype=torch.float64))
    Q = alg.identity(1) - P
    probs = st.born_probabilities(rho, [P, Q]).numpy()
    Dp = dense(P.matrix, ks, ms)
    assert np.allclose(probs[:, 0], np.trace(Dr @ Dp, axis1=-2, axis2=-1).real) and np.allclose(probs.sum(-1), 1)
    post, pr = rho.condition_on(P)
    assert np.allclose(pr.numpy(), probs[:, 0])
    ref = Dp @ Dr @ Dp / probs[:, 0, None, None]
    assert np.allclose(dense_state(post, ks, ms), ref, atol=1e-10)
    assert np.allclose(st.fidelity(rho, rho).numpy(), 1.0) and np.allclose(st.trace_distance(rho, rho).numpy(), 0)


def test_mix_states():
    """mix: (1 - p) rho + p sigma with a number or one weight per sample; p outside [0, 1] raises."""
    alg = alg_of([2, 3])
    rho, sigma = st.random_density_matrix(alg, 2), st.random_density_matrix(alg, 2)
    assert torch.allclose(rho.mix(sigma, 0.25).matrix, 0.75 * rho.matrix + 0.25 * sigma.matrix)
    p = [0.1, 0.9]
    m = rho.mix(sigma, p).matrix
    assert torch.allclose(m[1], 0.1 * rho.matrix[1] + 0.9 * sigma.matrix[1], atol=1e-16)
    with pytest.raises(ValueError):
        rho.mix(sigma, 1.5)
    with pytest.raises(ValueError):
        rho.mix(sigma, torch.tensor([0.5, -0.1]))


# ================================================================================================================
# channels.py
# ================================================================================================================
def dense_kraus(chan, ks, ms):
    """Kraus operators of the channel on the full space: (r, B, N, N)."""
    return np.stack([dense(chan.kraus[:, i], ks, ms) for i in range(chan.kraus.shape[1])])


def apply_dense(Ks, D):
    return sum(K @ D @ K.conj().swapaxes(-1, -2) for K in Ks)


@layouts
def test_channel_apply_dual_and_properties(ks, ms):
    """Channel.apply = sum_i K_i rho K_i^* on the full space; duality Tr(Phi(rho) A) = Tr(rho Phi^*(A)); TP / unital."""
    alg = alg_of(ks)
    rho = st.random_density_matrix(alg, 2)
    A = ginibre(alg, 2)
    for chan in (ch.random_channel(alg, 3, batch_size=2), ch.random_mixed_unitary_channel(alg, 3, batch_size=2),
                 ch.dephasing_channel(alg, 0.3), ch.depolarizing_channel(alg, 0.4), ch.amplitude_damping_channel(alg, 0.6),
                 ch.identity_channel(alg)):
        Ks = dense_kraus(chan, ks, ms)
        out = chan(rho)
        assert isinstance(out, st.DensityMatrix)
        assert np.allclose(dense(out.matrix, ks, ms), apply_dense(Ks, dense(rho.matrix, ks, ms)), atol=1e-12)
        lhs = st.trace_of_product(out, A)
        rhs = st.trace_of_product(rho, chan.adjoint()(A))
        assert torch.allclose(lhs, rhs) and chan.is_trace_preserving()
        I = np.eye(Ks.shape[-1])
        assert chan.is_unital() == np.allclose(apply_dense(Ks, I[None]), I, atol=1e-10)


@layouts
def test_named_channels_explicit_formulas(ks, ms):
    """dephasing, depolarizing (= centre expectation at p = 1), amplitude damping: explicit block formulas."""
    alg = alg_of(ks)
    rho = st.random_density_matrix(alg, 1)
    R = blocks(rho.matrix, ks)
    p, g = 0.3, 0.6
    deph = blocks(ch.dephasing_channel(alg, p)(rho).matrix, ks)
    depo = blocks(ch.depolarizing_channel(alg, p)(rho).matrix, ks)
    cent = blocks(ch.center_expectation(alg)(rho).matrix, ks)
    amp = blocks(ch.amplitude_damping_channel(alg, g)(rho).matrix, ks)
    for c, k in enumerate(ks):
        r = R[c][0]
        assert np.allclose(deph[c][0], (1 - p) * r + p * np.diag(np.diag(r)))
        assert np.allclose(depo[c][0], (1 - p) * r + p * np.trace(r) * np.eye(k) / k)
        assert np.allclose(cent[c][0], np.trace(r) * np.eye(k) / k)
        K0 = np.diag([1.0] + [math.sqrt(1 - g)] * (k - 1))
        ref = K0 @ r @ K0
        for i in range(1, k):
            Ki = np.zeros((k, k))
            Ki[0, i] = math.sqrt(g)
            ref = ref + Ki @ r @ Ki.T
        assert np.allclose(amp[c][0], ref)
    for bad in (-0.1, 1.1):
        with pytest.raises(ValueError):
            ch.depolarizing_channel(alg, bad)


@layouts
def test_choi_superoperator_and_reconstruction(ks, ms):
    """choi J = sum_ij |i><j| (x) Phi(|i><j|) (PSD, Tr_out J = 1), superoperator on row-major vec, and
    from_superoperator reproduces the channel; the transpose map (not CP) is rejected."""
    alg = alg_of(ks)
    chan = ch.random_channel(alg, 2, batch_size=2)
    J, S = chan.choi(), chan.superoperator()
    for c, k in enumerate(ks):
        Kc = chan.kraus[:, :, c, :k, :k].numpy()                                  # (B, r, k, k)
        for b in range(2):
            Jref = np.zeros((k * k, k * k), dtype=complex)
            for i in range(k):
                for j in range(k):
                    E = np.zeros((k, k))
                    E[i, j] = 1
                    Jref += np.kron(E, sum(K @ E @ K.conj().T for K in Kc[b]))
            assert np.allclose(J[c][b].numpy(), Jref, atol=1e-12)
            assert np.linalg.eigvalsh(Jref).min() > -1e-12
            assert np.allclose(np.trace(Jref.reshape(k, k, k, k), axis1=1, axis2=3), np.eye(k))
            Sref = sum(np.kron(K, K.conj()) for K in Kc[b])
            assert np.allclose(S[c][b].numpy(), Sref, atol=1e-12)
    rho = st.random_density_matrix(alg, 2)
    back = ch.Channel.from_superoperator(alg, S)
    assert torch.allclose(back(rho).matrix, chan(rho).matrix, atol=1e-12)
    T = [torch.eye(k * k, dtype=CD).reshape(k, k, k, k).permute(0, 1, 3, 2).reshape(k * k, k * k)[None]
         for k in ks]                                                             # rho -> rho^T
    if max(ks) > 1:
        with pytest.raises(ValueError):
            ch.Channel.from_superoperator(alg, T)


def test_channel_compose_mix_kraus_unitary():
    """compose = Phi(Psi(.)), mix = (1-p) Phi + p Psi, kraus_channel / unitary_channel / measurement_channel."""
    ks, ms = [2, 3], [1, 2]
    alg = alg_of(ks)
    rho = st.random_density_matrix(alg, 2)
    P1, P2 = ch.random_channel(alg, 2), ch.random_channel(alg, 3)
    assert torch.allclose((P1 @ P2)(rho).matrix, P1(P2(rho)).matrix, atol=1e-12)
    assert torch.allclose(P1.mix(P2, 0.3)(rho).matrix, 0.7 * P1(rho).matrix + 0.3 * P2(rho).matrix, atol=1e-12)
    U = alg.random_unitary_operator(2)
    assert torch.allclose(ch.unitary_channel(U)(rho).matrix, (U @ rho @ U.adjoint()).matrix, atol=1e-12)
    P = alg.operator_from_eigenvalues(lambda k: torch.tensor([1.0] + [0.0] * (k - 1), dtype=torch.float64))
    Q = alg.identity(1) - P
    M = ch.measurement_channel([P, Q])
    assert M.is_trace_preserving() and M.is_unital()
    assert torch.allclose(M(rho).matrix, (P @ rho @ P + Q @ rho @ Q).matrix, atol=1e-12)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        out = ch.kraus_channel([P])(rho)
    assert not isinstance(out, st.DensityMatrix) and len(w) == 1


def isc_dense(isc):
    """Full Kraus matrices of an InterSectorChannel, (B, r, sum m_d, sum k_c), blocks K^{dc}."""
    ki, ko = isc.algebra_in.k_factors, isc.algebra_out.k_factors
    K = isc.kraus.numpy()
    B, r = K.shape[:2]
    out = np.zeros((B, r, sum(ko), sum(ki)), dtype=K.dtype)
    for d, md in enumerate(ko):
        for c, kc in enumerate(ki):
            out[:, :, sum(ko[:d]):sum(ko[:d + 1]), sum(ki[:c]):sum(ki[:c + 1])] = K[:, :, d, c, :md, :kc]
    return out


def pinch(D, ks):
    out = np.zeros_like(D)
    for sl in slices(ks, [1] * len(ks)):
        out[..., sl, sl] = D[..., sl, sl]
    return out


def test_inter_sector_channel_against_dense():
    """InterSectorChannel: Phi(rho) = pinching of sum_i K_i rho K_i^* (full block Kraus matrices); duality for every
    trace; compose; transition matrix; TP and unitality checks."""
    ki, ko = [2, 3, 1], [3, 2]
    a_in, a_out = alg_of(ki), alg_of(ko)
    isc = ch.random_inter_sector_channel(a_in, a_out, kraus_rank=2, batch_size=2)
    assert isc.is_trace_preserving()
    rho = st.random_density_matrix(a_in, 2)
    K = isc_dense(isc)
    R = dense(rho.matrix, ki, [1] * 3)
    ref = pinch(np.einsum('brij,bjk,brlk->bil', K, R, K.conj()), ko)
    out = isc(rho)
    assert isinstance(out, st.DensityMatrix)
    assert np.allclose(dense(out.matrix, ko, [1, 1]), ref, atol=1e-12)
    A = ginibre(a_out, 2)
    for tr in st.TRACES:
        w_in, w_out = st._trace_weights(a_in, tr), st._trace_weights(a_out, tr)
        lhs = (torch.diagonal((out @ A).matrix, dim1=-2, dim2=-1).sum(-1) * w_out).sum(-1)
        rhs = (torch.diagonal((rho @ isc.adjoint(tr)(A)).matrix, dim1=-2, dim2=-1).sum(-1) * w_in).sum(-1)
        assert torch.allclose(lhs, rhs)
    T = isc.transition_matrix(rho)
    Tref = np.zeros((2, 2, 3))
    for c, sl_c in enumerate(slices(ki, [1] * 3)):
        Rc = np.zeros_like(R)
        Rc[:, sl_c, sl_c] = R[:, sl_c, sl_c]
        img = np.einsum('brij,bjk,brlk->bil', K, Rc, K.conj())
        for d, sl_d in enumerate(slices(ko, [1, 1])):
            Tref[:, d, c] = np.trace(img[:, sl_d, sl_d], axis1=-2, axis2=-1).real / \
                np.trace(Rc, axis1=-2, axis2=-1).real
    assert np.allclose(T.numpy(), Tref) and np.allclose(T.sum(1).numpy(), 1.0)
    isc2 = ch.random_inter_sector_channel(a_out, a_in, kraus_rank=1, batch_size=2)
    comp = isc2 @ isc
    assert torch.allclose(comp(rho).matrix, isc2(isc(rho)).matrix, atol=1e-12) and comp.is_trace_preserving()
    loc = ch.random_channel(a_in, 2)
    assert torch.allclose((isc @ loc)(rho).matrix, isc(loc(rho)).matrix, atol=1e-12)
    blocks_ = {(0, 1): [torch.ones(3, 3, dtype=CD) / 3]}
    fb = ch.InterSectorChannel.from_blocks(a_in, a_out, blocks_)
    assert not fb.is_trace_preserving()
    U = alg_of([2, 2]).random_unitary_operator(1)
    swap = ch.InterSectorChannel.from_blocks(alg_of([2, 2]), alg_of([2, 2]),
                                             {(1, 0): [U.matrix[0, 0]], (0, 1): [U.matrix[0, 1]]})
    assert swap.is_trace_preserving() and swap.is_unital()


# ================================================================================================================
# dynamics.py
# ================================================================================================================
@layouts
def test_propagator_schrodinger_von_neumann(ks, ms):
    """propagator = expm(-i t H), schrodinger psi(t) = U(t) psi0 and von_neumann U rho U^* against dense expm
    (one initial state for a batch of Hamiltonians)."""
    alg = alg_of(ks)
    H = hermitian(alg, 3)
    D = dense(H.matrix, ks, ms)
    t = 0.83
    Ut = np.stack([sla.expm(-1j * t * d) for d in D])
    assert np.allclose(dense(dy.propagator(H, t).matrix, ks, ms), Ut, atol=1e-10)
    psi0 = torch.randn(1, len(ks), max(ks), dtype=CD) * alg._const('rows')
    psi0 = psi0 / psi0.norm()
    out = dy.schrodinger(psi0, H, [0.0, t])
    assert out.shape == (2, 3, len(ks), max(ks))
    flat = lambda v: np.concatenate([v[..., c, :k] for c, k in enumerate(ks)], -1)
    ref = np.einsum('bij,j->bi', np.stack([sla.expm(-1j * t * d) for d in blocks_dense(H, ks)]), flat(psi0.numpy())[0])
    assert np.allclose(flat(out[1].numpy()), ref, atol=1e-10) and np.allclose(flat(out[0].numpy()), flat(psi0.numpy()))
    rho = st.random_density_matrix(alg, 1)
    traj = dy.von_neumann(rho, H, [0.0, t])
    assert isinstance(traj[1], st.DensityMatrix)
    Dr = dense(rho.matrix, ks, ms)
    assert np.allclose(dense(traj[1].matrix, ks, ms), Ut @ Dr @ Ut.conj().swapaxes(-1, -2), atol=1e-10)
    with pytest.raises(ValueError):
        dy.propagator(alg_of(ks, False).operator(torch.zeros(1, len(ks), max(ks), max(ks), dtype=torch.float64)), 1.0)


def blocks_dense(op, ks):
    """(B, sum k, sum k): the operator on (+)_c C^{k_c} (multiplicity 1)."""
    return dense(op.matrix, ks, [1] * len(ks))


def test_rk4_order_and_solvers():
    """rk4 has global error O(h^4); schrodinger_rk4 and solve_operator_ode (Heisenberg equation) match exact solutions."""
    ks = [2, 3]
    alg = alg_of(ks)
    A = torch.tensor([[0.0, 1.0], [-1.0, 0.0]], dtype=torch.float64)
    y0 = torch.tensor([1.0, 0.0], dtype=torch.float64)
    exact = torch.tensor(sla.expm(2.0 * A.numpy()) @ y0.numpy())
    errs = [(dy.rk4(lambda t, y: A @ y, y0, [0.0, 2.0], substeps=n)[-1] - exact).norm().item() for n in (10, 20)]
    assert 12 < errs[0] / errs[1] < 20                                              # 2^4 = 16
    H = hermitian(alg, 2)
    psi0 = torch.randn(2, 2, 3, dtype=CD) * alg._const('rows')
    psi0 = psi0 / psi0.flatten(1).norm(dim=1)[:, None, None]
    a = dy.schrodinger_rk4(psi0, H, [0.0, 0.5, 1.0], substeps=200)
    b = dy.schrodinger(psi0, H, [0.0, 0.5, 1.0])
    assert torch.allclose(a, b, atol=1e-9)
    X0 = ginibre(alg, 2)
    Hm = H.matrix
    sol = dy.solve_operator_ode(lambda t, X: 1j * (Hm @ X - X @ Hm), X0, [0.0, 0.7], substeps=200)
    U = dy.propagator(H, 0.7)
    assert torch.allclose(sol[-1].matrix, (U.adjoint() @ X0 @ U).matrix, atol=1e-9)


def gksl_dense(Hd, Ls, gs):
    """Row-major vec generator on one block."""
    n = Hd.shape[-1]
    I = np.eye(n)
    G = -1j * (np.kron(Hd, I) - np.kron(I, Hd.T))
    for g, L in zip(gs, Ls):
        LdL = L.conj().T @ L
        G = G + g * (np.kron(L, L.conj()) - 0.5 * np.kron(LdL, I) - 0.5 * np.kron(I, LdL.T))
    return G


@layouts
def test_lindblad_superoperator_channel_and_evolve(ks, ms):
    """lindblad_superoperator per block = -i(H(x)1 - 1(x)H^T) + sum g (L(x)conj L - ...); lindblad_channel = exp(tL)
    applied to rho; lindblad_evolve (RK4) agrees; trace preserved."""
    alg = alg_of(ks)
    H = hermitian(alg, 1)
    L1, L2 = ginibre(alg, 1), ginibre(alg, 1)
    rates = [0.3, 0.7]
    G = dy.lindblad_superoperator(alg, H, [L1, L2], rates)
    Hb, A1, A2 = blocks(H.matrix, ks), blocks(L1.matrix, ks), blocks(L2.matrix, ks)
    rho = st.random_density_matrix(alg, 1)
    R = blocks(rho.matrix, ks)
    t = 0.4
    chan = dy.lindblad_channel(alg, H, [L1, L2], t, rates)
    out = blocks(chan(rho).matrix, ks)
    ev = dy.lindblad_evolve(rho, H, [L1, L2], [0.0, t], rates, substeps=100)[-1]
    evb = blocks(ev.matrix, ks)
    for c, k in enumerate(ks):
        Gref = gksl_dense(Hb[c][0], [A1[c][0], A2[c][0]], rates)
        assert np.allclose(G[c][0].numpy(), Gref, atol=1e-12)
        ref = (sla.expm(t * Gref) @ R[c][0].reshape(-1)).reshape(k, k)
        assert np.allclose(out[c][0], ref, atol=1e-10)
        assert np.allclose(evb[c][0], ref, atol=1e-8)
    assert np.allclose(ev.matrix.diagonal(dim1=-2, dim2=-1).sum((-2, -1)).real.numpy(), 1.0)
    assert chan.is_trace_preserving()


def test_lindblad_rhs_with_jumps_between_sectors():
    """lindblad_rhs with an InterSectorChannel jump: drho/dt = E(K rho K^* - 1/2 {K^* K, rho}) with the pinching E
    onto M (full block Kraus matrix); the total trace is conserved by lindblad_evolve."""
    ks = [1, 2, 1]
    alg = alg_of(ks)
    Kb = torch.zeros(1, 1, 3, 3, 2, 2, dtype=CD)
    Kb[0, 0, 0, 1, 0, :2] = torch.tensor([0.6, 0.8j])               # sector 1 -> sector 0
    Kb[0, 0, 2, 1, 0, :2] = torch.tensor([0.3, -0.2])               # sector 1 -> sector 2
    Kb[0, 0, 1, 2, :2, 0] = torch.tensor([1.0, 0.5])                # sector 2 -> sector 1
    jump = ch.InterSectorChannel(alg, alg, Kb)
    H = hermitian(alg, 1)
    rho = st.random_density_matrix(alg, 1)
    f = dy.lindblad_rhs(H, [jump], rates=[0.9])
    out = blocks_dense(alg.operator(f(0.0, rho.matrix)), ks)[0]
    K = isc_dense(jump)[0, 0]
    R, Hd = blocks_dense(rho, ks)[0], blocks_dense(H, ks)[0]
    ref = -1j * (Hd @ R - R @ Hd) + 0.9 * pinch(K @ R @ K.conj().T - 0.5 * (K.conj().T @ K @ R + R @ K.conj().T @ K), ks)
    assert np.allclose(out, ref, atol=1e-12)
    traj = dy.lindblad_evolve(rho, H, [jump], [0.0, 0.5, 1.0], rates=[0.9], substeps=50)
    for r in traj:
        assert abs(r.matrix.diagonal(dim1=-2, dim2=-1).sum().real.item() - 1) < 1e-10
    with pytest.raises(ValueError):
        dy.lindblad_channel(alg, H, [jump])


# ================================================================================================================
# composite.py
# ================================================================================================================
def product_embedding(alg12):
    """For every (sector s, row i) of alg12 the index of e_a (x) e_b in H1 (x) H2 (multiplicity 1)."""
    a1, a2 = alg12.tensor_factors
    off1 = np.cumsum([0] + a1.k_factors)
    off2 = np.cumsum([0] + a2.k_factors)
    n2 = sum(a2.k_factors)
    emb = {}
    for s, pairs in enumerate(alg12.product_blocks):
        for c, d, off in pairs:
            for a in range(a1.k_factors[c]):
                for b in range(a2.k_factors[d]):
                    emb[(s, off + a * a2.k_factors[d] + b)] = (off1[c] + a) * n2 + off2[d] + b
    return emb


def to_full(X, alg12):
    a1, a2 = alg12.tensor_factors
    n = sum(a1.k_factors) * sum(a2.k_factors)
    emb = product_embedding(alg12)
    m = X.matrix.numpy()
    out = np.zeros((m.shape[0], n, n), dtype=m.dtype)
    for (s, i), gi in emb.items():
        for (s2, j), gj in emb.items():
            if s2 == s:
                out[:, gi, gj] = m[:, s, i, j]
    return out


@pytest.mark.parametrize('fused', [False, True])
def test_kron_and_partial_trace_against_full_space(fused):
    """kron = A (x) B on H1 (x) H2; partial_trace = conditional expectation onto alg1 (alg2) of the Hilbert-space
    partial trace of X; partial_trace_channel gives the same map and its blunt adjoint is A -> A (x) 1."""
    a1 = tv.TypeIAlgebra([1, 2, 1], [1, 2, 1], precision='double', device='cpu', charges=[0, 1, 2])
    a2 = tv.TypeIAlgebra([2, 1], [2, 1], precision='double', device='cpu', charges=[0, 1])
    a12 = cp.fused_tensor_product(a1, a2) if fused else cp.tensor_product(a1, a2)
    A, B = ginibre(a1, 2), ginibre(a2, 2)
    K = cp.kron(A, B, a12)
    full = to_full(K, a12)
    Ad, Bd = blocks_dense(A, a1.k_factors), blocks_dense(B, a2.k_factors)
    assert np.allclose(full, np.stack([np.kron(Ad[b], Bd[b]) for b in range(2)]))
    X = ginibre(a12, 2)
    Xf = to_full(X, a12)
    n1, n2 = sum(a1.k_factors), sum(a2.k_factors)
    T = Xf.reshape(2, n1, n2, n1, n2)
    for keep, alg_k, red in ((1, a1, np.einsum('bajcj->bac', T)), (2, a2, np.einsum('bjajc->bac', T))):
        ref = pinch(red, alg_k.k_factors)
        pt = cp.partial_trace(X, keep)
        assert np.allclose(blocks_dense(pt, alg_k.k_factors), ref, atol=1e-12)
        chan = cp.partial_trace_channel(a12, keep)
        assert torch.allclose(chan(X).matrix, pt.matrix, atol=1e-12) and chan.is_trace_preserving()
    emb = cp.partial_trace_channel(a12, 1).adjoint()(A)
    assert torch.allclose(emb.matrix, cp.kron(A, a2.identity(2), a12).matrix, atol=1e-12)
    rho, sigma = st.random_density_matrix(a1, 1), st.random_density_matrix(a2, 1)
    rs = cp.kron(rho, sigma, a12)
    assert isinstance(rs, st.DensityMatrix)
    assert torch.allclose(cp.partial_trace(rs, 1).matrix, rho.matrix, atol=1e-12)
    assert torch.allclose(cp.partial_trace(rs, 2).matrix, sigma.matrix, atol=1e-12)
    if fused:
        assert a12.charges == [0, 1, 2, 3] and a12.k_factors == [2, 5, 4, 1]


# ================================================================================================================
# krylov.py
# ================================================================================================================
def test_krylov_ground_state_and_evolution():
    """SparseSectorHamiltonian.matvec = dense product; ground_state = lowest eigenpair (batch of disorder);
    evolve = expm(-i t H) psi0."""
    chain = SpinChain(8, 'periodic')
    h = torch.rand(3, 8, dtype=torch.float64) * 2 - 1
    Hs = chain.xxz_sparse(1.0, 0.7, h, sector=4)
    Hd = Hs.to_dense()
    V = torch.randn(3, Hs.dim, 2, dtype=CD)
    assert torch.allclose(Hs.matvec(V), Hd @ V, atol=1e-12)
    E, psi, res = krylov.ground_state(Hs, m=40)
    w = torch.linalg.eigvalsh(Hd)
    assert torch.allclose(E, w[:, 0], atol=1e-9) and (res < 1e-8).all()
    assert torch.allclose(torch.linalg.vector_norm(Hd @ psi[..., None] - E[:, None, None] * psi[..., None], dim=(1, 2)),
                          torch.zeros(3, dtype=torch.float64), atol=1e-8)
    psi0 = torch.randn(3, Hs.dim, dtype=CD)
    psi0 = psi0 / psi0.norm(dim=-1, keepdim=True)
    out = krylov.evolve(Hs, psi0, [0.0, 0.5, 1.3], m=20)
    for i, t in enumerate([0.0, 0.5, 1.3]):
        ref = torch.matrix_exp(-1j * t * Hd) @ psi0[..., None]
        assert torch.allclose(out[i], ref[..., 0], atol=1e-8)


def test_krylov_evolve_large_norm_regression():
    """Regression: the default step had a lower bound 1e-3, so for ||H|| >> 100 m every Krylov step was far too
    long (wrong result, only a warning); now dt = 0.1 m / ||H|| as documented."""
    d = torch.linspace(-2e4, 2e4, 60, dtype=torch.float64)[None]
    H = krylov.SparseSectorHamiltonian(d, None, CD)
    psi0 = torch.randn(1, 60, dtype=CD)
    psi0 = psi0 / psi0.norm()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', cost.CostWarning)
        out = krylov.evolve(H, psi0, [0.0, 2e-3], m=20)
    assert torch.allclose(out[-1, 0], torch.exp(-2e-3j * d[0]) * psi0[0], atol=1e-9)


def test_lanczos_breakdown_and_tridiagonal():
    """lanczos: Q orthonormal, Q H Q^* = T; a start vector in a 2-dimensional invariant subspace breaks down with
    beta = 0 and the ground state is still exact."""
    torch.manual_seed(3)
    k = 12
    X = torch.randn(k, k, dtype=torch.float64)
    Hd = X + X.T
    Hm = lambda v: v @ Hd.T
    alpha, beta, Q = krylov.lanczos(Hm, torch.randn(1, k, dtype=torch.float64), 6)
    assert torch.allclose(Q[0] @ Q[0].T, torch.eye(6, dtype=torch.float64), atol=1e-12)
    assert torch.allclose(Q[0] @ Hd @ Q[0].T, krylov._tridiag(alpha, beta)[0], atol=1e-10)
    w, V = torch.linalg.eigh(Hd)
    v0 = (V[:, 0] + V[:, 5])[None]
    alpha, beta, Q = krylov.lanczos(Hm, v0, 5)
    assert beta[0, 1:].abs().max() == 0 and (Q[0, 2:] == 0).all()

    class Dense:
        batch_size, dim, dtype, device = 1, k, torch.float64, torch.device('cpu')
        matvec = staticmethod(Hm)
    E, psi, _ = krylov.ground_state(Dense(), m=5, v0=v0)
    assert abs(E.item() - w[0].item()) < 1e-12


# ================================================================================================================
# optimize.py
# ================================================================================================================
def test_parametrisations_satisfy_constraints():
    """UnitaryParam is unitary; SpectrumParam has max 1, min lo = (1-delta)/(1+delta); PositiveParam has Michelson
    contrast delta; SelfAdjointParam: |X| = abs_operator with that contrast."""
    ks = [2, 3, 1]
    alg = alg_of(ks)
    U = optimize.UnitaryParam(alg, 4)
    with torch.no_grad():
        U.G.normal_()
    Ud = dense(U.operator().matrix, ks, [1, 1, 1])
    assert np.allclose(Ud @ Ud.conj().swapaxes(-1, -2), np.eye(6), atol=1e-12)
    delta = torch.tensor([0.0, 0.3, 0.6, 0.95], dtype=torch.float64)
    for P in (optimize.PositiveParam(alg, 4, delta), optimize.PositiveParam(alg, 4, delta, pin_sector=1)):
        lam = P.spectrum.eigenvalues()
        act = torch.cat([lam[:, c, :k] for c, k in enumerate(ks)], -1)
        assert torch.allclose(act.amax(-1), torch.ones(4, dtype=torch.float64))
        assert torch.allclose(act.amin(-1), (1 - delta) / (1 + delta))
        assert torch.allclose(P.operator().michelson_contrast, delta, atol=1e-12)
    S = optimize.SelfAdjointParam(alg, 4, delta)
    X = S.operator()
    Xd = dense(X.matrix, ks, [1, 1, 1])
    absd = dense_fun(Xd, np.abs)
    assert np.allclose(dense(S.abs_operator().matrix, ks, [1, 1, 1]), absd, atol=1e-12)
    w = np.linalg.eigvalsh(absd)
    assert np.allclose((w.max(-1) - w.min(-1)) / (w.max(-1) + w.min(-1)), delta.numpy())
    free = optimize.SpectrumParam(alg, 4, None, (0.2, 0.7))
    d = free.delta()
    assert ((d >= 0.2) & (d <= 0.7)).all()


def test_extremize_reaches_von_neumann_trace_bound():
    """extremize over unitaries: max_U Re Tr(U A) = ||A||_1 (von Neumann); the result never exceeds it."""
    ks = [2, 3]
    alg = alg_of(ks)
    A = ginibre(alg, 1)
    bound = A.trace_norm().item()
    U = optimize.UnitaryParam(alg, 4)
    best = optimize.extremize(lambda: (U.operator() @ A).Tr_blunt().real, [U], maximize=True, rounds=2, steps=40,
                              lr=0.1)
    assert (best <= bound + 1e-9).all() and best.max().item() > bound - 1e-3


# ================================================================================================================
# chains.py
# ================================================================================================================
def dense_xxz(L, J, Delta, h, bonds):
    """XXZ Hamiltonian on (C^2)^{(x) L} with bit i of the index = spin of site i (1 = up)."""
    n = 2 ** L
    H = np.zeros((n, n))
    for x in range(n):
        s = [((x >> i) & 1) - 0.5 for i in range(L)]
        H[x, x] = sum(J * Delta * s[i] * s[j] for i, j in bonds) + sum(h[i] * s[i] for i in range(L))
        for i, j in bonds:
            if ((x >> i) & 1) != ((x >> j) & 1):
                H[x ^ ((1 << i) | (1 << j)), x] += J / 2
    return H


@pytest.mark.parametrize('L,boundary', [(2, 'periodic'), (4, 'open'), (5, 'periodic')])
def test_xxz_against_dense_spin_chain(L, boundary):
    """SpinChain.xxz (all sectors, one sector, sparse) = the dense XXZ Hamiltonian restricted to the labelled basis."""
    chain = SpinChain(L, boundary)
    h = np.random.uniform(-1, 1, (2, L))
    H = chain.xxz(0.8, 1.3, torch.tensor(h))
    for b in range(2):
        Hd = dense_xxz(L, 0.8, 1.3, h[b], chain.bonds)
        for N in range(L + 1):
            lab = chain.labels[N].numpy()
            k = len(lab)
            ref = Hd[np.ix_(lab, lab)]
            assert np.allclose(H.matrix[b, N, :k, :k].numpy(), ref)
            assert np.allclose(chain.xxz(0.8, 1.3, torch.tensor(h), sector=N).matrix[b, 0].numpy(), ref)
            assert np.allclose(chain.xxz_sparse(0.8, 1.3, torch.tensor(h), sector=N).to_dense()[b].numpy(), ref)
        # off-sector entries of the dense H vanish (U(1) symmetry)
        pop = np.array([bin(x).count('1') for x in range(2 ** L)])
        assert np.abs(Hd[pop[:, None] != pop[None, :]]).max() == 0
    sz = chain.sz(1)
    Nop = chain.number()
    for N in range(L + 1):
        lab = chain.labels[N].numpy()
        assert np.allclose(sz.matrix[0, N, :len(lab), :len(lab)].diagonal().real.numpy(), ((lab >> 1) & 1) - 0.5)
        assert np.allclose(Nop.matrix[0, N, :len(lab), :len(lab)].diagonal().real.numpy(), N)
    if L == 2:
        assert chain.bonds == [(0, 1)]


@pytest.mark.parametrize('L,N', [(6, 3), (6, 2), (5, 2)])
def test_momentum_blocks(L, N):
    """xxz_momentum: Hermitian blocks whose spectra together are the spectrum of the sector N."""
    chain = SpinChain(L, 'periodic')
    Hk = chain.xxz_momentum(1.0, 0.6, 0.3, sector=N)
    assert Hk.algebra.total_subspace_dim == math.comb(L, N)
    assert torch.allclose(Hk.matrix, Hk.matrix.mH, atol=1e-12)
    wk = torch.cat(Hk.eigenvalues(), -1).sort(-1).values
    w = torch.linalg.eigvalsh(chain.xxz(1.0, 0.6, 0.3, sector=N).matrix[0, 0])
    assert torch.allclose(wk[0], w, atol=1e-10)


def test_reduced_state_and_site_channels():
    """reduced_state = numpy partial trace over sites n_left..L-1 of the dense state; site_amplitude_damping and
    lowering act as the dense Kraus maps K_0 = |0><0| + sqrt(1-g)|1><1|, K_1 = sqrt(g) S^-_i on (C^2)^{(x) L}."""
    L, nl = 5, 2
    chain = SpinChain(L, 'open')
    H = chain.xxz(1.0, 0.5, torch.tensor(np.random.uniform(-1, 1, L)))
    rho = st.gibbs_state(H, 0.8)
    n = 2 ** L
    R = np.zeros((n, n), dtype=complex)
    for N in range(L + 1):
        lab = chain.labels[N].numpy()
        R[np.ix_(lab, lab)] = rho.matrix[0, N, :len(lab), :len(lab)].numpy()
    # index x = a + 2^nl * b (a: sites 0..nl-1)
    T = R.reshape(2 ** (L - nl), 2 ** nl, 2 ** (L - nl), 2 ** nl)
    red = np.einsum('babc->ac', T)
    left = SpinChain(nl, 'open')
    out = chain.reduced_state(rho, nl)
    for N in range(nl + 1):
        lab = left.labels[N].numpy()
        assert np.allclose(out.matrix[0, N, :len(lab), :len(lab)].numpy(), red[np.ix_(lab, lab)], atol=1e-12)
    S_ref = ent(red[None])[0]
    assert abs(chain.entanglement_entropy(rho, nl).item() - S_ref) < 1e-10
    i, g = 2, 0.35
    Km = np.zeros((n, n))
    K0 = np.zeros((n, n))
    for x in range(n):
        if (x >> i) & 1:
            Km[x ^ (1 << i), x] = 1.0
            K0[x, x] = math.sqrt(1 - g)
        else:
            K0[x, x] = 1.0
    ref = K0 @ R @ K0.T + g * Km @ R @ Km.T
    amp = chain.site_amplitude_damping(i, g)
    out = amp(rho)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')                 # S^- alone is not trace preserving
        low = chain.lowering(i)(rho)
    for N in range(L + 1):
        lab = chain.labels[N].numpy()
        k = len(lab)
        assert np.allclose(out.matrix[0, N, :k, :k].numpy(), ref[np.ix_(lab, lab)], atol=1e-12)
        assert np.allclose(low.matrix[0, N, :k, :k].numpy(), (Km @ R @ Km.T)[np.ix_(lab, lab)], atol=1e-12)
    assert amp.is_trace_preserving()
    bs = chain.basis_state('01100')
    x = 0b00110
    N = 2
    p = int(chain.index_of[x])
    assert bs.matrix[0, N, p, p] == 1 and bs.matrix.abs().sum() == 1
    assert chain.vector_in_sector('01100')[p] == 1


def test_level_spacing_ratio_exact_values():
    """Proposition 5.14 values: <r> = 2 ln 2 - 1 (Poisson), 4 - 2 sqrt 3 (GOE 3x3), 2 sqrt 3/pi - 1/2 (GUE 3x3)."""
    N = 20000
    pois = torch.cumsum(torch.distributions.Exponential(torch.tensor(1.0, dtype=torch.float64)).sample((N, 3)), -1)
    r = SpinChain.level_spacing_ratio(pois)
    assert zscore(r, 2 * math.log(2) - 1) < 5
    X = torch.randn(N, 3, 3, dtype=torch.float64)
    r = SpinChain.level_spacing_ratio(torch.linalg.eigvalsh(X + X.mT))
    assert zscore(r, 4 - 2 * math.sqrt(3)) < 5
    Z = torch.randn(N, 3, 3, dtype=CD)
    r = SpinChain.level_spacing_ratio(torch.linalg.eigvalsh(Z + Z.mH))
    assert zscore(r, 2 * math.sqrt(3) / math.pi - 0.5) < 5


def test_chain_basis_is_iterated_fused_product():
    """The basis order of SpinChain is that of the iterated fused product of single sites M_1 (+) M_1 (charges 0, 1):
    S^z_i = 1 (x) ... (x) S^z (x) ... (x) 1 built with composite.kron equals chain.sz(i)."""
    L = 4
    site = tv.TypeIAlgebra([1, 1], [1, 1], precision='double', device='cpu', charges=[0, 1])
    sz1, one = site.central([-0.5, 0.5]), site.identity()
    alg, ops = site, [sz1]
    for _ in range(1, L):
        new = cp.fused_tensor_product(alg, site)
        ops = [cp.kron(o, one, new) for o in ops] + [cp.kron(alg.identity(), sz1, new)]
        alg = new
    chain = SpinChain(L)
    assert alg.k_factors == chain.algebra.k_factors
    for i in range(L):
        assert torch.equal(ops[i].matrix, chain.sz(i).matrix)


def test_cost_probe_subtracts_fixed_call_cost():
    """Regression: the time estimate of batched_call extrapolated an 8-item slice linearly, so a fixed cost per call
    (ms for batched LAPACK kernels) gave a ~100x overestimate and a false CostWarning (random_unitary(5, batch=1e5)).
    Simulated kernels: a fixed 20 ms per call and no per-item cost -> no warning; 1e-4 s per item on 2e5 items
    (20 s) -> a warning with an estimate of about 20 s. The full calls are instantaneous."""
    import time
    n = 200000
    x = torch.zeros(n, 1, 1)

    def fixed_cost(y):
        if y.shape[0] < n:
            time.sleep(0.02)
        return y

    def per_item_cost(y):
        if y.shape[0] < n:
            time.sleep(1e-4 * y.shape[0])
        return y

    old = cost.get_limits()
    try:
        cost.set_limits(probe_flops=0, warn_seconds=10.0)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            cost.batched_call(fixed_cost, x, "fixed", kind='matmul')
        assert not [m for m in w if issubclass(m.category, cost.CostWarning)]
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            cost.batched_call(per_item_cost, x, "per item", kind='matmul')
        msg = [str(m.message) for m in w if issubclass(m.category, cost.CostWarning)]
        assert len(msg) == 1 and 'about' in msg[0]
        eta = float(msg[0].split('about ')[1].split(' ')[0])
        assert 15 <= eta <= 60                                   # 20 s; sleep only overshoots on a loaded machine
    finally:
        cost.set_limits(**old)


def test_kronecker_products_accept_non_contiguous_operands():
    """Kronecker products of transposed / conjugated views (torch.kron fails on mixed memory layouts in PyTorch
    2.14): dynamics._kron, the Lindblad superoperator of non-contiguous inputs and composite.kron agree with np.kron."""
    A = torch.randn(2, 3, 3, dtype=CD)
    B = torch.randn(2, 3, 3, dtype=CD)
    for X, Y in ((A, B.mT), (A.mH, B), (A.mT, B.conj()), (A[:, :, :], B.transpose(-2, -1).conj())):
        ref = np.stack([np.kron(X[b].resolve_conj().numpy(), Y[b].resolve_conj().numpy()) for b in range(2)])
        assert np.allclose(dy._kron(X, Y).numpy(), ref)
    ks = [2, 3]
    alg = alg_of(ks)
    H = hermitian(alg, 1)
    L = ginibre(alg, 1)
    Hnc, Lnc = H.matrix.mT.conj(), L.matrix.mH.mH       # the same matrices as non-contiguous / conjugate views
    G1 = dy.lindblad_superoperator(alg, H, [L])
    G2 = dy.lindblad_superoperator(alg, Hnc, [Lnc])
    assert all(torch.allclose(a, b) for a, b in zip(G1, G2))
    a1, a2 = alg_of([2, 1]), alg_of([1, 2])
    X, Y = ginibre(a1, 1), ginibre(a2, 1)
    a12 = cp.tensor_product(a1, a2)
    K1 = cp.kron(a1.operator(X.matrix.mT.contiguous().mT), a2.operator(Y.matrix.mH.mH), a12)
    assert torch.allclose(K1.matrix, cp.kron(X, Y, a12).matrix)


def test_regression_krylov_step_ignores_constant_shift(monkeypatch):
    """krylov.evolve: the default step depends on the half-width of the spectrum, so H + 1e5 takes as many steps as H."""
    from torch_vn_algebra import krylov
    from torch_vn_algebra.chains import SpinChain
    ch = SpinChain(8, 'open', complex_valued=False, precision='double')
    calls = []
    step = krylov.expm_multiply_step
    monkeypatch.setattr(krylov, 'expm_multiply_step', lambda *a, **k: calls.append(1) or step(*a, **k))
    psi0 = torch.ones(1, ch.xxz_sparse(1.0, 0.5, 0.0, sector=2).dim, dtype=torch.float64)
    psi0 = psi0 / psi0.norm()
    for h in (0.0, 1e5):                                         # field h: H -> H + const in the sector N = 2
        H = ch.xxz_sparse(1.0, 0.5, h, sector=2)
        calls.clear()
        out = krylov.evolve(H, psi0, [0.0, 1.0])[-1, 0].numpy()
        w, V = np.linalg.eigh(H.to_dense()[0].numpy())
        assert h == 0 or np.abs(w).min() > 1e4                  # the shift is really there
        ref = V @ (np.exp(-1j * w) * (V.T @ psi0[0].numpy()))
        assert np.abs(out - ref).max() < 1e-9
        assert len(calls) <= 5
