"""
Theory tests for Chapter 6 (quantum channels), Chapter 7 (tensor products), Chapter 8 (dynamics) and Appendix D
(generators of quantum dynamical semigroups, operator algebra error correction) of the lecture notes, checked
through the library (CPU, double precision, fixed seeds). Monte Carlo checks use z-scores against the standard error.
"""
import math

import numpy as np
import pytest
import scipy.linalg as sla
import torch

from torch_vn_algebra import (TypeIAlgebra, DensityMatrix, Channel, InterSectorChannel, channels, dynamics, states,
                              composite, inclusions, tensor_product, fused_tensor_product, kron, partial_trace)

C128, F64 = torch.complex128, torch.float64
LAYOUTS = [[1], [3], [1, 2], [2, 3], [1, 2, 3], [3, 3]]
Z = 5.0                                                   # z-score bound for Monte Carlo checks


# ----------------------------------------------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------------------------------------------
def mk(ks, charges=None, n=None):
    return TypeIAlgebra(list(n or ks), list(ks), precision='double', charges=charges)


def seed(s):
    torch.manual_seed(s)
    return np.random.default_rng(s)


def cplx(rng, *shape):
    return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)


def blocks_of(op, b=0):
    return [op.matrix[b, c, :k, :k].detach().numpy().copy() for c, k in enumerate(op.algebra.k_factors)]


def op_of(alg, blocks):
    return alg.from_blocks([torch.as_tensor(np.asarray(B, dtype=complex)) for B in blocks])


def full(op, b=0):
    return sla.block_diag(*blocks_of(op, b))


def unflat(alg, X):
    """Block-diagonal part of a D x D matrix as an operator of alg."""
    out, o = [], 0
    for k in alg.k_factors:
        out.append(X[o:o + k, o:o + k])
        o += k
    return op_of(alg, out)


def rand_state(rng, k, rank=None):
    G = cplx(rng, k, rank or k)
    R = G @ G.conj().T
    return R / np.trace(R).real


def unit(i, j, k):
    E = np.zeros((k, k), complex)
    E[i, j] = 1
    return E


def choi(f, k):
    """J = sum_ij E_ij (x) f(E_ij) (Theorem 6.2)."""
    return sum(np.kron(unit(i, j, k), f(unit(i, j, k))) for i in range(k) for j in range(k))


def superop(f, k):
    """Row-major matrix of a linear map on k x k matrices."""
    return np.stack([f(unit(i, j, k)).reshape(-1) for i in range(k) for j in range(k)], axis=1)


def choi_from_superop(S, k):
    return S.reshape(k, k, k, k).transpose(2, 0, 3, 1).reshape(k * k, k * k)


def block_map(ch, c):
    """The block Phi_cc of a sector-preserving channel as a numpy function."""
    alg = ch.algebra

    def f(X):
        bl = [np.zeros((k, k), complex) for k in alg.k_factors]
        bl[c] = X
        return blocks_of(ch(op_of(alg, bl)))[c]
    return f


def entropy(w):
    w = np.clip(np.real(w), 0, None)
    w = w[w > 1e-15]
    return float(-(w * np.log(w)).sum())


def zscore(x, exact):
    x = np.asarray(x, dtype=float)
    return abs(x.mean() - exact) / (x.std(ddof=1) / math.sqrt(x.size))


def vec_K(K):
    """(1 (x) K) Omega = sum_i e_i (x) K e_i (Lemma D.1(b))."""
    return K.T.reshape(-1)


def mat_K(v, n):
    return v.reshape(n, n).T


def kron_rows(alg12):
    k1, k2 = alg12.tensor_factors[0].k_factors, alg12.tensor_factors[1].k_factors
    return [np.array([(sum(k1[:c]) + a) * sum(k2) + sum(k2[:d]) + b
                      for c, d, off in pairs for a in range(k1[c]) for b in range(k2[d])])
            for pairs in alg12.product_blocks]


def full12(op, b=0):
    """Operator of a (plain or fused) product algebra as a matrix on H1 (x) H2 (Kronecker order)."""
    f1, f2 = op.algebra.tensor_factors
    N = sum(f1.k_factors) * sum(f2.k_factors)
    out = np.zeros((N, N), complex)
    for s, r in enumerate(kron_rows(op.algebra)):
        out[np.ix_(r, r)] = op.matrix[b, s, :len(r), :len(r)].numpy()
    return out


def from_full12(alg12, X):
    mat = torch.zeros(1, alg12.C, alg12.k_max, alg12.k_max, dtype=C128)
    for s, r in enumerate(kron_rows(alg12)):
        mat[0, s, :len(r), :len(r)] = torch.from_numpy(X[np.ix_(r, r)])
    return alg12.operator(mat)


def random_op(alg, rng):
    return op_of(alg, [cplx(rng, k, k) for k in alg.k_factors])


def herm(rng, k):
    X = cplx(rng, k, k)
    return (X + X.conj().T) / 2


def gksl_np(H, Ls, rho):
    out = -1j * (H @ rho - rho @ H)
    for L in Ls:
        out = out + L @ rho @ L.conj().T - 0.5 * (L.conj().T @ L @ rho + rho @ L.conj().T @ L)
    return out


# ================================================================================================================
# Chapter 6: quantum channels
# ================================================================================================================
@pytest.mark.parametrize('k', [2, 3, 4])
def test_prop_6_1_transpose_positive_not_2_positive(k):
    """Proposition 6.1: the transpose is positive, its Choi matrix is the swap with eigenvalue -1."""
    rng = seed(k)
    for _ in range(5):
        A = rand_state(rng, k, rank=1)
        assert np.linalg.eigvalsh(A.T).min() > -1e-12
    J = choi(lambda X: X.T, k)
    F = np.zeros((k * k, k * k))
    for i in range(k):
        for j in range(k):
            F[j * k + i, i * k + j] = 1
    assert np.abs(J - F).max() == 0
    v = np.zeros(k * k)
    v[0 * k + 1], v[1 * k + 0] = 1, -1                     # e1 (x) e2 - e2 (x) e1
    assert np.allclose(J @ v, -v)
    assert abs(np.linalg.eigvalsh(J).min() + 1) < 1e-12
    alg = mk([k])
    with pytest.raises(ValueError):
        Channel.from_superoperator(alg, [torch.from_numpy(superop(lambda X: X.T, k))[None]])


@pytest.mark.parametrize('ks', LAYOUTS)
@pytest.mark.parametrize('r', [1, 2, 5])
def test_thm_6_2_choi_kraus_and_cor_6_3_rank(ks, r):
    """Theorem 6.2 and Corollary 6.3: the Choi matrix is positive, its rank min(r, k^2) is the Kraus rank, and the
    Kraus operators read off its eigenvectors reproduce the channel."""
    rng = seed(10 * len(ks) + r)
    alg = mk(ks)
    Phi = channels.random_channel(alg, kraus_rank=r)
    Js, Ss = Phi.choi(), Phi.superoperator()
    for c, k in enumerate(ks):
        f = block_map(Phi, c)
        J = Js[c][0].numpy()
        assert np.abs(J - choi(f, k)).max() < 1e-12              # convention of Theorem 6.2
        assert np.abs(Ss[c][0].numpy() - superop(f, k)).max() < 1e-12
        w, V = np.linalg.eigh(J)
        assert w.min() > -1e-12
        rank = int((w > 1e-10).sum())
        assert rank == min(r, k * k)                             # generic Kraus rank of a random channel
        Ks = [math.sqrt(x) * mat_K(V[:, i], k) for i, x in enumerate(w) if x > 1e-10]
        X = cplx(rng, k, k)
        assert np.abs(sum(K @ X @ K.conj().T for K in Ks) - f(X)).max() < 1e-12   # (c) => (d)
        assert sum(np.trace(K.conj().T @ K) for K in Ks).real == pytest.approx(k)  # TP
        assert len(Ks) <= r                                      # Corollary 6.3: r Kraus operators >= rank
    x = random_op(alg, rng)                                      # back from the superoperator (Choi eigenvectors)
    assert (Channel.from_superoperator(alg, Ss)(x).matrix - Phi(x).matrix).abs().max() < 1e-12


@pytest.mark.parametrize('k', [1, 2, 3])
def test_thm_6_2_k_positivity_detects_cp(k):
    """Theorem 6.2 (b) <=> (c): Phi_k applied to |Omega><Omega| is the Choi matrix; the transpose is k-positive only
    for k = 1."""
    Omega = sum(np.kron(np.eye(k)[:, [i]], np.eye(k)[:, [i]]) for i in range(k))
    P = Omega @ Omega.conj().T
    theta_k = np.zeros_like(P)
    for i in range(k):
        for j in range(k):
            theta_k[i * k:(i + 1) * k, j * k:(j + 1) * k] = P[i * k:(i + 1) * k, j * k:(j + 1) * k].T
    assert np.abs(theta_k - choi(lambda X: X.T, k)).max() == 0
    assert (np.linalg.eigvalsh(theta_k).min() >= -1e-12) == (k == 1)


def test_choi_hermitian_iff_hermiticity_preserving():
    """Section 6.1: J(Phi)^* = J(Phi) iff Phi(X^*) = Phi(X)^*; from_superoperator refuses X -> iX."""
    rng = seed(1)
    k = 3
    A, B = cplx(rng, k, k), cplx(rng, k, k)
    hp = lambda X: A @ X @ A.conj().T - B @ X @ B.conj().T     # Hermiticity preserving, not CP
    nhp = lambda X: A @ X @ B.conj().T
    Jh, Jn = choi(hp, k), choi(nhp, k)
    assert np.abs(Jh - Jh.conj().T).max() < 1e-12 and np.abs(Jn - Jn.conj().T).max() > 0.1
    X = cplx(rng, k, k)
    assert np.abs(hp(X.conj().T) - hp(X).conj().T).max() < 1e-12
    assert np.abs(nhp(X.conj().T) - nhp(X).conj().T).max() > 0.1
    alg = mk([k])
    with pytest.raises(ValueError):
        Channel.from_superoperator(alg, [1j * torch.eye(k * k, dtype=C128)[None]])
    with pytest.raises(ValueError):                              # Hermiticity preserving but not CP
        Channel.from_superoperator(alg, [torch.from_numpy(superop(hp, k))[None]])


@pytest.mark.parametrize('ks_in,ks_out', [([2, 3], [1, 2, 2]), ([1], [2, 3]), ([3, 1], [3]), ([2, 2], [1, 3])])
def test_prop_6_4_inter_sector_kraus_and_tp(ks_in, ks_out):
    """Proposition 6.4: Phi(A)_d = sum_{c,l} K^{dc}_l A_c K^{dc*}_l, TP iff sum_{d,l} K^{dc*} K^{dc} = 1 for every c,
    and each Phi_dc is CP."""
    rng = seed(len(ks_in) + 7 * len(ks_out))
    a_in, a_out = mk(ks_in), mk(ks_out)
    Psi = channels.random_inter_sector_channel(a_in, a_out, kraus_rank=2)
    K = Psi.kraus[0].numpy()
    x = random_op(a_in, rng)
    xb = blocks_of(x)
    out = blocks_of(Psi(x))
    for d, kd in enumerate(ks_out):
        ref = sum(K[l, d, c, :kd, :kc] @ xb[c] @ K[l, d, c, :kd, :kc].conj().T
                  for l in range(K.shape[0]) for c, kc in enumerate(ks_in))
        assert np.abs(out[d] - ref).max() < 1e-12
    for c, kc in enumerate(ks_in):
        gram = sum(K[l, d, c, :kd, :kc].conj().T @ K[l, d, c, :kd, :kc]
                   for l in range(K.shape[0]) for d, kd in enumerate(ks_out))
        assert np.abs(gram - np.eye(kc)).max() < 1e-12
        for d, kd in enumerate(ks_out):                          # every Phi_dc is CP
            J = sum(np.kron(unit(i, j, kc), sum(K[l, d, c, :kd, :kc] @ unit(i, j, kc) @ K[l, d, c, :kd, :kc].conj().T
                                                for l in range(K.shape[0])))
                    for i in range(kc) for j in range(kc))
            assert np.linalg.eigvalsh(J).min() > -1e-12
    assert Psi.is_trace_preserving()
    rho = states.random_density_matrix(a_in)
    assert abs(Psi(rho).trace.real.item() - 1) < 1e-12


@pytest.mark.parametrize('ks', [[2], [1, 3], [2, 3]])
def test_thm_6_5_stinespring(ks):
    """Theorem 6.5: V = sum_l e_l (x) K_l gives Phi = Tr_E(V . V^*), Phi^* = V^*(1 (x) .)V, and V is an isometry
    iff Phi is TP."""
    rng = seed(sum(ks))
    alg = mk(ks)
    r = 3
    Phi = channels.random_channel(alg, kraus_rank=r)
    notp = Channel(alg, 0.8 * Phi.kraus)
    for c, k in enumerate(ks):
        K = Phi.kraus[0, :, c, :k, :k].numpy()
        V = K.reshape(r * k, k)                                   # environment index first
        assert np.abs(V.conj().T @ V - np.eye(k)).max() < 1e-12
        A = cplx(rng, k, k)
        TrE = np.einsum('lalb->ab', (V @ A @ V.conj().T).reshape(r, k, r, k))
        assert np.abs(TrE - block_map(Phi, c)(A)).max() < 1e-12
        B = cplx(rng, k, k)
        assert np.abs(V.conj().T @ np.kron(np.eye(r), B) @ V - block_map(Phi.adjoint(), c)(B)).max() < 1e-12
        V2 = notp.kraus[0, :, c, :k, :k].numpy().reshape(r * k, k)
        assert np.abs(V2.conj().T @ V2 - np.eye(k)).max() > 0.1
    assert not notp.is_trace_preserving()


@pytest.mark.parametrize('ks_in,ks_out', [([2, 3], [1, 2, 2]), ([1, 2], [3])])
def test_stinespring_between_sectors(ks_in, ks_out):
    """Section 6.2: V_c = sum_{d,l} e_l (x) iota'_d K^{dc}_l is an isometry for a TP map between sectors, and the
    (d,d) block of sum_c Tr_E(V_c A_c V_c^*) is Phi(A)_d."""
    rng = seed(3)
    a_in, a_out = mk(ks_in), mk(ks_out)
    r = 2
    Psi = channels.random_inter_sector_channel(a_in, a_out, kraus_rank=r)
    K = Psi.kraus[0].numpy()
    Dout = sum(ks_out)
    offs = np.cumsum([0] + ks_out)
    x = random_op(a_in, rng)
    xb = blocks_of(x)
    total = np.zeros((Dout, Dout), complex)
    for c, kc in enumerate(ks_in):
        V = np.zeros((r * Dout, kc), complex)
        for l in range(r):
            for d, kd in enumerate(ks_out):
                V[l * Dout + offs[d]:l * Dout + offs[d] + kd] = K[l, d, c, :kd, :kc]
        assert np.abs(V.conj().T @ V - np.eye(kc)).max() < 1e-12
        total += np.einsum('lalb->ab', (V @ xb[c] @ V.conj().T).reshape(r, Dout, r, Dout))
    out = blocks_of(Psi(x))
    for d, kd in enumerate(ks_out):
        assert np.abs(total[offs[d]:offs[d] + kd, offs[d]:offs[d] + kd] - out[d]).max() < 1e-12


@pytest.mark.parametrize('p,n', [(3, 2), (4, 4), (2, 5)])
def test_lemma_6_6_unitary_from_equal_gram(p, n):
    """Lemma 6.6: XX^* = YY^* implies Y = XW with W unitary (constructed as in the proof)."""
    rng = seed(p * n)
    X = cplx(rng, p, min(p, n) - 1) @ cplx(rng, min(p, n) - 1, n)    # rank deficient on purpose
    W0 = sla.qr(cplx(rng, n, n))[0]
    Y = X @ W0
    assert np.abs(X @ X.conj().T - Y @ Y.conj().T).max() < 1e-10
    # isometry X^* z -> Y^* z on the range of X^*, extended to a unitary U with U X^* = Y^*
    Ux, sx, _ = np.linalg.svd(X.conj().T)
    Uy, sy, _ = np.linalg.svd(Y.conj().T)
    r = int((sx > 1e-10).sum())
    Rx = X.conj().T @ np.linalg.pinv(X.conj().T) @ Ux[:, :r]      # orthonormal basis of ran X^*
    pre = np.linalg.pinv(X.conj().T) @ Rx                         # z with X^* z = basis vectors
    img = Y.conj().T @ pre
    Cx = null(Rx.conj().T)
    Cy = null(img.conj().T)
    U = np.hstack([img, Cy]) @ np.hstack([Rx, Cx]).conj().T
    assert np.abs(U.conj().T @ U - np.eye(n)).max() < 1e-10
    W = U.conj().T
    assert np.abs(Y - X @ W).max() < 1e-10


@pytest.mark.parametrize('ks', [[2], [2, 3], [1, 3]])
def test_thm_6_7_unitary_freedom(ks):
    """Theorem 6.7: L_j = sum_l u_jl K_l gives the same channel; conversely two Kraus families of one channel are
    related by a unitary (after padding)."""
    rng = seed(11 + sum(ks))
    alg = mk(ks)
    Phi = channels.random_channel(alg, kraus_rank=2)
    n = 4
    Kpad = torch.cat([Phi.kraus, torch.zeros_like(Phi.kraus[:, :n - 2])], dim=1)
    u = torch.from_numpy(sla.qr(cplx(rng, n, n))[0])
    PhiL = Channel(alg, torch.einsum('jl,blcxy->bjcxy', u, Kpad))
    rho = states.random_density_matrix(alg)
    assert (PhiL(rho).matrix - Phi(rho).matrix).abs().max() < 1e-12
    for c, k in enumerate(ks):                                    # converse: recover u by least squares
        Kv = np.stack([vec_K(Kpad[0, l, c, :k, :k].numpy()) for l in range(n)], 1)
        Lv = np.stack([vec_K(PhiL.kraus[0, j, c, :k, :k].numpy()) for j in range(n)], 1)
        assert np.abs(Kv @ Kv.conj().T - Lv @ Lv.conj().T).max() < 1e-12     # equal Choi matrices
        if k * k >= 2:
            W = np.linalg.lstsq(Kv[:, :2], Lv, rcond=None)[0]         # w_j = sum_l W_lj v_l
            assert np.abs(Kv[:, :2] @ W - Lv).max() < 1e-12
            assert np.abs(W @ W.conj().T - np.eye(2)).max() < 1e-12      # extends to a unitary
    # non-equivalent families give different channels
    v = sla.qr(cplx(rng, n, n))[0] * np.array([0.5, 1, 1, 1])   # not unitary on the nonzero K_l
    PhiV = Channel(alg, torch.einsum('jl,blcxy->bjcxy', torch.from_numpy(v), Kpad))
    x = alg.operator(rho.matrix)
    assert (PhiV(x).matrix - Phi(x).matrix).abs().max() > 1e-3


TRACE_W = {'blunt': lambda ks: np.ones(len(ks)), 'norm': lambda ks: 1 / np.array(ks, float),
           'tau_vN': lambda ks: 1 / (len(ks) * np.array(ks, float))}


def trw(op, w):
    return sum(wc * np.trace(B) for wc, B in zip(w, blocks_of(op)))


@pytest.mark.parametrize('trace', ['blunt', 'norm', 'tau_vN'])
@pytest.mark.parametrize('ks_in,ks_out', [([2, 3], [1, 2, 2]), ([1, 2], [3, 1]), ([3], [2, 2])])
def test_prop_6_8_dual_formula(trace, ks_in, ks_out):
    """Proposition 6.8: tr_v(Phi(x) y) = tr_w(x Phi#(y)) with Phi#(y) = z_w^{-1} Phi^*(z_v y); Phi# is CP with
    Kraus operators (v_d/w_c)^{1/2} K^{dc*}."""
    rng = seed(len(ks_in) * 13 + len(ks_out))
    a_in, a_out = mk(ks_in), mk(ks_out)
    Psi = channels.random_inter_sector_channel(a_in, a_out, kraus_rank=2)
    w, v = TRACE_W[trace](ks_in), TRACE_W[trace](ks_out)
    dual = Psi.adjoint(trace=trace)
    for _ in range(3):
        x, y = random_op(a_in, rng), random_op(a_out, rng)
        lhs = trw(Psi(x) @ y, v)
        rhs = trw(x @ dual(y), w)
        assert abs(lhs - rhs) < 1e-10
        zw, zv = a_in.central(list(w)), a_out.central(list(v))
        ref = zw.inverse() @ Psi.adjoint()(zv @ y)
        assert (dual(y).matrix - ref.matrix).abs().max() < 1e-12
    K = Psi.kraus[0].numpy()
    yb = blocks_of(y)
    got = blocks_of(dual(y))
    for c, kc in enumerate(ks_in):
        ref = sum(v[d] / w[c] * K[l, d, c, :kd, :kc].conj().T @ yb[d] @ K[l, d, c, :kd, :kc]
                  for l in range(2) for d, kd in enumerate(ks_out))
        assert np.abs(got[c] - ref).max() < 1e-12
    P = rand_state(rng, sum(ks_out))                          # Phi# is positive: maps positive y to positive
    pos = dual(unflat(a_out, P))
    assert min(np.linalg.eigvalsh(B).min() for B in blocks_of(pos)) > -1e-12


@pytest.mark.parametrize('trace', ['blunt', 'norm', 'tau_vN'])
def test_cor_6_9_dual_properties(trace):
    """Corollary 6.9: (a) TP w.r.t. (w,v) iff Phi#(1) = 1; (b) unital iff Phi# TP; (c) Phi## = Phi; (d) sector-
    preserving maps with v = w have Phi# = Phi^* for every trace."""
    rng = seed(5)
    a_in, a_out = mk([2, 3]), mk([1, 2, 2])
    Psi = channels.random_inter_sector_channel(a_in, a_out, kraus_rank=2)
    w, v = TRACE_W[trace](a_in.k_factors), TRACE_W[trace](a_out.k_factors)
    one_dev = (Psi.adjoint(trace=trace)(a_out.identity()).matrix - a_in.identity().matrix).abs().max().item()
    x = random_op(a_in, rng)
    tp_dev = max(abs(trw(Psi(B), v) - trw(B, w)) for B in [x, a_in.identity(), random_op(a_in, rng)])
    assert (one_dev < 1e-12) == (tp_dev < 1e-10)                  # (a)
    assert (one_dev < 1e-12) == (trace == 'blunt')
    dd = Psi.adjoint(trace=trace).adjoint(trace=trace)            # (c)
    assert (dd(x).matrix - Psi(x).matrix).abs().max() < 1e-12
    alg = mk([2, 3])
    for ch, unital in ((channels.random_mixed_unitary_channel(alg, 3), True),
                       (channels.amplitude_damping_channel(alg, 0.4), False)):
        isc = ch.to_inter_sector()
        for t in ('blunt', 'norm', 'tau_vN'):                      # (d)
            y = random_op(alg, rng)
            assert (isc.adjoint(trace=t)(y).matrix - ch.adjoint()(y).matrix).abs().max() < 1e-12
        assert ch.is_unital() == unital
        assert ch.adjoint().is_trace_preserving() == unital         # (b)


def test_particle_loss_dual_depends_on_trace():
    """Section 6.3: particle loss Phi(x) = (x_0 + Tr x_1) (+) 0 is TP for the blunt traces, Phi#(1) = 1 (+) 2 1_2 for
    Tr_norm, and Tr_norm(Phi(x)) = x_0 + Tr x_1 differs from Tr_norm(x)."""
    rng = seed(2)
    a = mk([1, 2])
    loss = InterSectorChannel.from_blocks(a, a, {(0, 0): [torch.ones(1, 1)],
                                                 (0, 1): [torch.tensor([[1.0, 0.0]]), torch.tensor([[0.0, 1.0]])]})
    assert loss.is_trace_preserving()
    assert (loss.adjoint('norm')(a.identity()).matrix - a.central([1.0, 2.0]).matrix).abs().max() < 1e-15
    for _ in range(3):
        x = random_op(a, rng)
        xb = blocks_of(x)
        out = blocks_of(loss(x))
        assert abs(out[0][0, 0] - (xb[0][0, 0] + np.trace(xb[1]))) < 1e-12 and np.abs(out[1]).max() == 0
        assert abs(loss(x).Tr_norm().item() - (xb[0][0, 0] + np.trace(xb[1]))) < 1e-12
        assert abs(x.Tr_norm().item() - (xb[0][0, 0] + 0.5 * np.trace(xb[1]))) < 1e-12
        assert abs(loss(x).Tr_norm().item() - (x @ a.central([1.0, 2.0])).Tr_norm().item()) < 1e-12   # Cor 6.9(a)
    assert torch.allclose(loss.transition_matrix()[0], torch.tensor([[1.0, 1.0], [0.0, 0.0]], dtype=F64))


@pytest.mark.parametrize('trace', ['norm', 'tau_vN'])
def test_densities_in_heisenberg_picture(trace):
    """Section 6.3 (densities in the Heisenberg picture): omega o Psi has the density Psi#(rho), and Psi# is TP when
    Psi is unital."""
    rng = seed(8)
    a_M, a_N = mk([2, 3]), mk([1, 2, 2])
    Phi = channels.random_inter_sector_channel(a_M, a_N, kraus_rank=2)
    Psi = Phi.adjoint()                                           # unital map N -> M (blunt dual of a TP map)
    assert Psi.is_unital()
    w, v = TRACE_W[trace](a_M.k_factors), TRACE_W[trace](a_N.k_factors)
    rho = random_op(a_M, rng)
    Psis = Psi.adjoint(trace=trace)                               # dual w.r.t. (v, w): M -> N
    for _ in range(3):
        y = random_op(a_N, rng)
        assert abs(trw(Psi(y) @ rho, w) - trw(Psis(rho) @ y, v)) < 1e-10
    assert abs(trw(Psis(rho), v) - trw(rho, w)) < 1e-10           # TP w.r.t. (w, v)


@pytest.mark.parametrize('ks', [[1, 2, 3], [3, 3]])
def test_table_6_4_standard_channels(ks):
    """Section 6.4: the Kraus forms and actions of dephasing, depolarizing, centre expectation and amplitude damping;
    all TP, unital except amplitude damping, sum K K^* = 1 + (k-1) gamma E00 - gamma sum_{i>0} E_ii."""
    rng = seed(4)
    alg = mk(ks)
    p, g = 0.37, 0.41
    rho = random_op(alg, rng)
    rb = blocks_of(rho)
    for name, ch, f in (
            ('deph', channels.dephasing_channel(alg, p), lambda R, k: (1 - p) * R + p * np.diag(np.diag(R))),
            ('depol', channels.depolarizing_channel(alg, p), lambda R, k: (1 - p) * R + p * np.trace(R) * np.eye(k) / k),
            ('centre', channels.center_expectation(alg), lambda R, k: np.trace(R) * np.eye(k) / k)):
        out = blocks_of(ch(rho))
        for c, k in enumerate(ks):
            assert np.abs(out[c] - f(rb[c], k)).max() < 1e-12, name
        assert ch.is_trace_preserving() and ch.is_unital()
    ad = channels.amplitude_damping_channel(alg, g)
    out = blocks_of(ad(rho))
    for c, k in enumerate(ks):
        R = rb[c]
        ref = np.empty_like(R)
        for i in range(k):
            for j in range(k):
                if i == 0 and j == 0:
                    ref[i, j] = R[0, 0] + g * sum(R[l, l] for l in range(1, k))
                elif i == 0 or j == 0:
                    ref[i, j] = math.sqrt(1 - g) * R[i, j]
                else:
                    ref[i, j] = (1 - g) * R[i, j]
        assert np.abs(out[c] - ref).max() < 1e-12
        KKs = sum(K @ K.conj().T for K in ad.kraus[0, :, c, :k, :k].numpy())
        ref_KK = np.eye(k) - g * np.diag([0] + [1] * (k - 1)) + (k - 1) * g * unit(0, 0, k)
        assert np.abs(KKs - ref_KK).max() < 1e-12
    assert ad.is_trace_preserving() and (ad.is_unital() == (max(ks) == 1))
    mu = channels.random_mixed_unitary_channel(alg, 3)
    assert mu.is_trace_preserving() and mu.is_unital()
    # the conditional expectation onto the centre is Proposition 2.8, and differs from the sector measurement
    E = channels.center_expectation(alg)(rho)
    assert (E.matrix - alg.central([np.trace(B) / k for B, k in zip(rb, ks)]).matrix).abs().max() < 1e-12


def null(M, tol=1e-9):
    """Orthonormal basis of the null space, with an absolute tolerance on the singular values."""
    _, s, Vh = np.linalg.svd(M)
    s = np.concatenate([s, np.zeros(Vh.shape[0] - s.size)])
    return Vh[s < tol].conj().T


def fixed_space(S, k):
    """Basis (columns, row-major vec) of {X : S vec X = vec X}."""
    return null(S - np.eye(k * k))


def commutant_space(ops, k):
    return null(np.concatenate([np.kron(A, np.eye(k)) - np.kron(np.eye(k), A.T) for A in ops]))


def same_space(A, B):
    if A.shape[1] != B.shape[1]:
        return False
    if A.shape[1] == 0:
        return True
    return np.linalg.matrix_rank(np.hstack([A, B]), tol=1e-8) == A.shape[1]


@pytest.mark.parametrize('ks', [[2, 3], [1, 4]])
def test_thm_6_10_fixed_points_of_unital_channels(ks):
    """Theorem 6.10: for unital channels the fixed points are exactly the commutant of {K_l, K_l^*} (the same space,
    not only the same dimension), and they form an algebra."""
    rng = seed(6)
    alg = mk(ks)
    D = [np.diag(np.exp(1j * rng.uniform(0, 6, k))) for k in ks]   # unitaries sharing an eigenbasis
    Ud = [alg.from_blocks([torch.from_numpy(np.diag(np.exp(1j * rng.uniform(0, 6, k)))) for k in ks]) for _ in range(3)]
    zoo = [channels.dephasing_channel(alg, 0.3), channels.depolarizing_channel(alg, 0.6),
           channels.center_expectation(alg), channels.random_mixed_unitary_channel(alg, 2),
           channels.kraus_channel([U * math.sqrt(1 / 3) for U in Ud])]
    for ch in zoo:
        assert ch.is_unital() and ch.is_trace_preserving()
        for c, k in enumerate(ks):
            S = ch.superoperator()[c][0].numpy()
            Ks = list(ch.kraus[0, :, c, :k, :k].numpy())
            F, Cm = fixed_space(S, k), commutant_space(Ks + [K.conj().T for K in Ks], k)
            assert same_space(F, Cm)
            if F.shape[1] >= 2:                                       # closed under products
                X, Y = F[:, 0].reshape(k, k), F[:, 1].reshape(k, k)
                XY = (X @ Y).reshape(-1)
                assert np.abs(S @ XY - XY).max() < 1e-9
    assert [fixed_space(S[0].numpy(), k).shape[1] for S, k in zip(zoo[0].superoperator(), ks)] == ks   # diagonal
    assert [fixed_space(S[0].numpy(), k).shape[1] for S, k in zip(zoo[1].superoperator(), ks)] == [1] * len(ks)
    assert [fixed_space(S[0].numpy(), k).shape[1] for S, k in zip(zoo[4].superoperator(), ks)] == ks
    del D


@pytest.mark.parametrize('k', [2, 3, 4])
@pytest.mark.parametrize('g', [0.3, 1.0])
def test_amplitude_damping_fixed_points_differ_from_commutant(k, g):
    """Section 6.4: without unitality Theorem 6.10 fails: amplitude damping fixes exactly C E00, its Kraus
    operators commute exactly with C 1."""
    alg = mk([k])
    ad = channels.amplitude_damping_channel(alg, g)
    F = fixed_space(ad.superoperator()[0][0].numpy(), k)
    Ks = [K for K in ad.kraus[0, :, 0, :k, :k].numpy()]
    Cm = commutant_space(Ks + [K.conj().T for K in Ks], k)
    assert F.shape[1] == 1 and Cm.shape[1] == 1
    assert same_space(F, unit(0, 0, k).reshape(-1, 1)) and same_space(Cm, np.eye(k).reshape(-1, 1))


def doubly_stochastic_map(alg_in, alg_out, rng):
    """A unital TP map between sectors with the same D: Kraus blocks P_d U P_c of a unitary U on C^D."""
    D = sum(alg_in.k_factors)
    U = sla.qr(cplx(rng, D, D))[0]
    oi, oo = np.cumsum([0] + alg_in.k_factors), np.cumsum([0] + alg_out.k_factors)
    blocks = {(d, c): [torch.from_numpy(U[oo[d]:oo[d + 1], oi[c]:oi[c + 1]].copy())]
              for d in range(alg_out.C) for c in range(alg_in.C)}
    return InterSectorChannel.from_blocks(alg_in, alg_out, blocks)


@pytest.mark.parametrize('ks_in,ks_out', [([2, 3], [2, 3]), ([1, 1, 3], [2, 3]), ([5], [1, 4]), ([2, 2], [4])])
def test_thm_6_11_unital_channels_increase_entropy(ks_in, ks_out):
    """Theorem 6.11: for a positive, TP, unital map mu = Delta lambda with Delta doubly stochastic, so S increases,
    purity and every sum of a concave function do too."""
    rng = seed(len(ks_in) + 3 * len(ks_out))
    a_in, a_out = mk(ks_in), mk(ks_out)
    Phi = doubly_stochastic_map(a_in, a_out, rng)
    assert Phi.is_trace_preserving() and Phi.is_unital()
    rho = states.random_density_matrix(a_in, rank=2)
    lam_vecs, lam = [], []
    for c, B in enumerate(blocks_of(rho)):
        w, V = np.linalg.eigh(B)
        for i in range(len(w)):
            bl = [np.zeros((k, k), complex) for k in ks_in]
            bl[c] = np.outer(V[:, i], V[:, i].conj())
            lam.append(w[i])
            lam_vecs.append(op_of(a_in, bl))
    out = Phi(rho)
    mu, f_vecs = [], []
    for d, B in enumerate(blocks_of(out)):
        w, V = np.linalg.eigh(B)
        mu += list(w)
        f_vecs += [(d, V[:, j]) for j in range(len(w))]
    Delta = np.array([[np.vdot(fv, blocks_of(Phi(E))[d] @ fv).real for E in lam_vecs] for d, fv in f_vecs])
    assert Delta.min() > -1e-12
    assert np.abs(Delta.sum(0) - 1).max() < 1e-12 and np.abs(Delta.sum(1) - 1).max() < 1e-12
    assert np.abs(Delta @ np.array(lam) - np.array(mu)).max() < 1e-12
    for f in (lambda t: -t * np.log(np.clip(t, 1e-300, None)), lambda t: -t ** 2, np.sqrt, lambda t: np.minimum(t, 0.2)):
        assert f(np.clip(mu, 0, None)).sum() >= f(np.clip(lam, 0, None)).sum() - 1e-12
    assert out.von_neumann_entropy().item() >= rho.von_neumann_entropy().item() - 1e-12
    assert out.purity().item() <= rho.purity().item() + 1e-12


def test_thm_6_11_needs_unitality():
    """Section 6.4: complete amplitude damping lowers the entropy of the maximally mixed state of M_2 (+) M_3 from
    log 5 to H(2/5, 3/5)."""
    alg = mk([2, 3])
    mm = states.maximally_mixed_state(alg)
    after = channels.amplitude_damping_channel(alg, 1.0)(mm).von_neumann_entropy().item()
    assert abs(mm.von_neumann_entropy().item() - math.log(5)) < 1e-12
    assert abs(after - entropy([0.4, 0.6])) < 1e-12


@pytest.mark.parametrize('ks_in,ks_out', [([2, 3], [1, 2, 2]), ([1, 3], [2]), ([2, 2, 1], [3, 1])])
def test_prop_6_12_transition_matrix(ks_in, ks_out):
    """Proposition 6.12: T(rho) is column stochastic, q = T(rho) p, and E_N Phi E_M(rho) = (+)_d (Tp)_d/k'_d 1."""
    a_in, a_out = mk(ks_in), mk(ks_out)
    seed(len(ks_in) + len(ks_out))
    Psi = channels.random_inter_sector_channel(a_in, a_out, kraus_rank=2)
    rho = states.random_density_matrix(a_in)
    T, Tr = Psi.transition_matrix()[0], Psi.transition_matrix(rho)[0]
    assert T.min() >= 0 and torch.allclose(T.sum(0), torch.ones(len(ks_in), dtype=F64))
    assert torch.allclose(Tr.sum(0), torch.ones(len(ks_in), dtype=F64))
    p, q = rho.sector_probabilities()[0], Psi(rho).sector_probabilities()[0]
    assert (q - Tr @ p).abs().max() < 1e-12
    classical = channels.center_expectation(a_out) @ Psi @ channels.center_expectation(a_in)
    Tp = T @ p
    assert (classical(rho).matrix - a_out.central([Tp[d] / k for d, k in enumerate(ks_out)]).matrix).abs().max() < 1e-12
    # T(1) by definition: Tr Phi_dc(1_c)/k_c
    K = Psi.kraus[0].numpy()
    for d, kd in enumerate(ks_out):
        for c, kc in enumerate(ks_in):
            ref = sum(np.trace(K[l, d, c, :kd, :kc] @ K[l, d, c, :kd, :kc].conj().T) for l in range(2)).real / kc
            assert abs(T[d, c].item() - ref) < 1e-12


def test_qubit_coherences_T2():
    """Section 6.4 (qubit coherences): D_p A_gamma = A_gamma D_p, coherence factor (1-p) sqrt(1-gamma),
    1/T2 = 1/(2 T1) + 1/T_phi, and both families compose additively in t."""
    rng = seed(9)
    q2 = mk([2])
    T1, Tphi = 1.7, 2.9
    r = op_of(q2, [rand_state(rng, 2)])
    gam, pd = (lambda t: 1 - math.exp(-t / T1)), (lambda t: 1 - math.exp(-t / Tphi))
    for t, s in ((0.3, 0.9), (1.2, 2.5)):
        D, A = channels.dephasing_channel(q2, pd(t)), channels.amplitude_damping_channel(q2, gam(t))
        assert (D(A(r)).matrix - A(D(r)).matrix).abs().max() < 1e-12
        coh = (D(A(r)).matrix[0, 0, 0, 1] / r.matrix[0, 0, 0, 1]).item()
        assert abs(coh - math.exp(-t * (1 / (2 * T1) + 1 / Tphi))) < 1e-12
        diag = D(A(r)).matrix[0, 0].diagonal().real.numpy()
        assert np.abs(diag - A(r).matrix[0, 0].diagonal().real.numpy()).max() < 1e-12
        for fam in (lambda x: channels.dephasing_channel(q2, pd(x)), lambda x: channels.amplitude_damping_channel(q2, gam(x))):
            assert (fam(t)(fam(s)(r)).matrix - fam(t + s)(r).matrix).abs().max() < 1e-12


@pytest.mark.parametrize('N,k', [(3, 1), (4, 2), (5, 5)])
def test_lemma_6_13_haar_isometry(N, k):
    """Lemma 6.13: random_isometry is an isometry distributed as k columns of a Haar unitary: E V_11 = E V_11^2 = 0,
    E|V_11|^2 = 1/N, E|V_11|^4 = 2/(N(N+1)), and UV has the law of V."""
    seed(N * 10 + k)
    alg = mk([2])
    B = 40_000
    V = alg.random_isometry(N, k, batch_size=B).numpy()
    assert np.abs(np.einsum('bij,bik->bjk', V.conj(), V) - np.eye(k)).max() < 1e-12
    U = sla.qr(np.random.default_rng(1).standard_normal((N, N)) + 1j * np.random.default_rng(2).standard_normal((N, N)))[0]
    for W in (V, np.einsum('ij,bjk->bik', U, V)):
        x = W[:, 0, k - 1]
        for samples, exact in ((x.real, 0), (x.imag, 0), ((x ** 2).real, 0), (np.abs(x) ** 2, 1 / N),
                               (np.abs(x) ** 4, 2 / (N * (N + 1)))):
            if N == 1 and exact == 0:
                continue
            assert zscore(samples, exact) < Z
    phases = np.angle(V[:, 0, 0])                                 # uniform phase (the sign correction of Theorem 5.3)
    assert zscore(np.cos(phases), 0) < Z and zscore(np.sin(phases), 0) < Z


@pytest.mark.parametrize('N', [2, 3, 6])
def test_lemma_6_14_fourth_moments(N):
    """Lemma 6.14: E u_a conj(u_b) u_c conj(u_d) = (d_ab d_cd + d_ad d_cb)/(N(N+1)) for a unitarily invariant unit
    vector, and E u_a conj(u_b) = d_ab/N."""
    seed(N)
    B = 60_000
    u = mk([1]).random_isometry(N, 1, batch_size=B)[:, :, 0].numpy()
    beta = 1 / (N * (N + 1))
    checks = [(np.abs(u[:, 0]) ** 4, 2 * beta), (np.abs(u[:, 0] * u[:, 1]) ** 2, beta),
              ((u[:, 0] * u[:, 1].conj()).real, 0), ((u[:, 0] ** 2 * u[:, 1].conj() ** 2).real, 0),
              ((u[:, 0] * u[:, 0].conj() * u[:, 0] * u[:, 1].conj()).real, 0), (np.abs(u[:, 1]) ** 2, 1 / N)]
    for samples, exact in checks:
        assert zscore(samples, exact) < Z
    assert abs(2 * beta + (N - 1) * beta - 1 / N) < 1e-15           # alpha + (N-1) beta = 1/N, alpha = 2 beta


@pytest.mark.parametrize('ks', [[1], [2, 3], [1, 4]])
@pytest.mark.parametrize('r', [1, 2, 4])
def test_thm_6_15_output_purity_of_random_channel(ks, r):
    """Theorem 6.15: E Tr Phi(psi)^2 = (k+r)/(kr+1) and E Phi(psi) = 1/k for a pure state in a sector of size k."""
    rng = seed(r * 7 + sum(ks))
    alg = mk(ks)
    B = 20_000
    Phi = channels.random_channel(alg, kraus_rank=r, batch_size=B)
    for c, k in enumerate(ks):
        psi = cplx(rng, k, 1)
        psi /= np.linalg.norm(psi)
        bl = [np.zeros((kk, kk), complex) for kk in ks]
        bl[c] = psi @ psi.conj().T
        out = Phi(op_of(alg, bl))
        pur = (out @ out).trace.real.numpy()
        exact = (k + r) / (k * r + 1)
        if r == 1 or k == 1:
            assert np.abs(pur - 1).max() < 1e-12                 # unitary channel (or k = 1)
        else:
            assert zscore(pur, exact) < Z
        mean = out.matrix.mean(0)[c, :k, :k].numpy()
        assert np.abs(mean - np.eye(k) / k).max() < 6 * math.sqrt(1 / B)
        others = [d for d in range(len(ks)) if d != c]
        assert not others or out.matrix[:, others].abs().max() == 0                  # sector preserved
    assert abs((3 + 4) / (3 * 4 + 1) - 7 / 13) < 1e-15            # the qutrit example of the box


@pytest.mark.parametrize('k', [2, 3, 4])
def test_exercise_6_1_werner_holevo(k):
    """Chapter 6, Exercise 1: the Werner-Holevo channel is a unital channel of Kraus rank k(k-1)/2, J = 2P_-/(k-1)."""
    rng = seed(k)
    f = lambda X: (np.trace(X) * np.eye(k) - X.T) / (k - 1)
    J = choi(f, k)
    F = choi(lambda X: X.T, k)
    assert np.abs(J - (np.eye(k * k) - F) / (k - 1)).max() < 1e-12
    w = np.linalg.eigvalsh(J)
    assert w.min() > -1e-12 and int((w > 1e-10).sum()) == k * (k - 1) // 2 and abs(np.trace(F) - k) < 1e-12
    alg = mk([k])
    WH = Channel.from_superoperator(alg, [torch.from_numpy(superop(f, k))[None]])
    assert WH.is_trace_preserving() and WH.is_unital()
    X = cplx(rng, k, k)
    assert np.abs(block_map(WH, 0)(X) - f(X)).max() < 1e-12


@pytest.mark.parametrize('ks', [[2, 3], [1, 4], [3]])
@pytest.mark.parametrize('p', [0.0, 0.25, 1.0])
def test_exercise_6_2_dephasing_kraus_rank(ks, p):
    """Chapter 6, Exercise 2: the dephasing channel has Kraus rank k for 0 < p <= 1 (1 for p = 0) while the library
    stores 1 + k_max operators."""
    alg = mk(ks)
    ch = channels.dephasing_channel(alg, p)
    ranks = [int((np.linalg.eigvalsh(J[0].numpy()) > 1e-10).sum()) for J in ch.choi()]
    assert ranks == (ks if p > 0 else [1] * len(ks))
    assert ch.kraus_rank == 1 + max(ks)
    for J, k in zip(ch.choi(), ks):                              # eigenvalues p and (1-p)k + p on span e_i (x) e_i
        w = np.sort(np.linalg.eigvalsh(J[0].numpy()))[::-1][:k]
        assert abs(w[0] - ((1 - p) * k + p)) < 1e-12 and np.abs(w[1:] - p).max(initial=0) < 1e-12


@pytest.mark.parametrize('k', [2, 3, 5])
def test_exercise_6_3_weyl_twirl(k):
    """Chapter 6, Exercise 3: (1/k^2) sum X^a Z^b rho Z^-b X^-a = Tr(rho) 1/k; the X^a Z^b / sqrt(k) are orthonormal."""
    rng = seed(k)
    Xs = np.roll(np.eye(k), 1, axis=0)
    Zs = np.diag(np.exp(2j * np.pi * np.arange(k) / k))
    W = [np.linalg.matrix_power(Xs, a) @ np.linalg.matrix_power(Zs, b) for a in range(k) for b in range(k)]
    G = np.array([[np.trace(A.conj().T @ B) for B in W] for A in W])
    assert np.abs(G - k * np.eye(k * k)).max() < 1e-10
    alg = mk([k])
    twirl = channels.kraus_channel([alg.from_blocks([torch.from_numpy(Wab / k)]) for Wab in W])
    rho = random_op(alg, rng)
    assert (twirl(rho).matrix - channels.center_expectation(alg)(rho).matrix).abs().max() < 1e-12


@pytest.mark.parametrize('trace', ['blunt', 'norm', 'tau_vN'])
def test_exercise_6_4_dual_of_composition(trace):
    """Chapter 6, Exercise 4: (Psi o Phi)# = Phi# o Psi# with the same trace on the middle algebra."""
    rng = seed(12)
    a, b, c = mk([2, 3]), mk([1, 2, 2]), mk([3, 1])
    Phi = channels.random_inter_sector_channel(a, b, kraus_rank=2)
    Psi = channels.random_inter_sector_channel(b, c, kraus_rank=2)
    z = random_op(c, rng)
    lhs = (Psi @ Phi).adjoint(trace=trace)(z)
    rhs = Phi.adjoint(trace=trace)(Psi.adjoint(trace=trace)(z))
    assert (lhs.matrix - rhs.matrix).abs().max() < 1e-12


# ================================================================================================================
# Chapter 7: tensor products and entanglement
# ================================================================================================================
PAIRS = [([1, 2], [2, 3]), ([2], [3]), ([1], [1, 1]), ([2, 1, 3], [2]), ([3, 2], [1, 2, 2])]


@pytest.mark.parametrize('k1,k2', PAIRS)
def test_thm_7_1_structure_of_product(k1, k2):
    """Theorem 7.1: sectors are the pairs (c,d) of sizes k_c l_d, dim is multiplicative, kron is the Kronecker product
    (multiplicative, *-preserving), and the three traces multiply."""
    rng = seed(len(k1) * 5 + len(k2))
    a1, a2 = mk(k1), mk(k2)
    a12 = tensor_product(a1, a2)
    assert a12.C == len(k1) * len(k2)
    assert a12.k_factors == [kc * ld for kc in k1 for ld in k2]
    assert sum(k * k for k in a12.k_factors) == sum(k * k for k in k1) * sum(l * l for l in k2)
    A, A2, B, B2 = random_op(a1, rng), random_op(a1, rng), random_op(a2, rng), random_op(a2, rng)
    AB = kron(A, B, a12)
    assert np.abs(full12(AB) - np.kron(full(A), full(B))).max() < 1e-12
    assert (AB @ kron(A2, B2, a12) - kron(A @ A2, B @ B2, a12)).matrix.abs().max() < 1e-12
    assert (AB.adjoint() - kron(A.adjoint(), B.adjoint(), a12)).matrix.abs().max() < 1e-12
    assert (kron(a1.identity(), a2.identity(), a12) - a12.identity()).matrix.abs().max() == 0
    for tr in ('Tr_blunt', 'Tr_norm', 'tau_vN'):
        assert abs(getattr(AB, tr)().item() - getattr(A, tr)().item() * getattr(B, tr)().item()) < 1e-10


@pytest.mark.parametrize('k1,k2', PAIRS)
def test_prop_7_2_partial_trace(k1, k2):
    """Proposition 7.2: Tr(Tr_2(rho) A) = Tr(rho (A (x) 1)), the Hilbert-space partial trace, Tr_2(A (x) B) =
    Tr(B) A, the bimodule property, trace preservation, and states go to states (also Tr_1)."""
    rng = seed(len(k1) + 11 * len(k2))
    a1, a2 = mk(k1), mk(k2)
    a12 = tensor_product(a1, a2)
    rho = states.random_density_matrix(a12, batch_size=2)
    D1, D2 = sum(k1), sum(k2)
    for b in range(2):
        R = full12(rho, b).reshape(D1, D2, D1, D2)
        r1, r2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
        assert np.abs(full(r1, b) - np.einsum('ijkj->ik', R)).max() < 1e-12
        assert np.abs(full(r2, b) - np.einsum('ijil->jl', R)).max() < 1e-12
    assert isinstance(partial_trace(rho), DensityMatrix)
    for b in range(2):
        assert min(np.linalg.eigvalsh(B).min() for B in blocks_of(partial_trace(rho), b)) > -1e-12
    A, A2, B, X = random_op(a1, rng), random_op(a1, rng), random_op(a2, rng), random_op(a12, rng)
    one2 = a2.identity()
    lhs = (partial_trace(X) @ A).Tr_blunt()
    rhs = (X @ kron(A, one2, a12)).Tr_blunt()
    assert abs(lhs.item() - rhs.item()) < 1e-10                             # the definition
    trB = sum(np.trace(Bc) for Bc in blocks_of(B))
    assert (partial_trace(kron(A, B, a12)).matrix - A.matrix * trB).abs().max() < 1e-12    # (b)
    left = partial_trace(kron(A, one2, a12) @ X @ kron(A2, one2, a12))
    assert (left.matrix - (A @ partial_trace(X) @ A2).matrix).abs().max() < 1e-10           # (c)
    assert abs(partial_trace(X).Tr_blunt().item() - X.Tr_blunt().item()) < 1e-10           # (d)
    assert (partial_trace(kron(A, B, a12), keep=2).matrix - B.matrix * sum(np.trace(c) for c in blocks_of(A))).abs().max() < 1e-12


@pytest.mark.parametrize('k1,k2', PAIRS)
def test_thm_7_3_partial_trace_is_channel(k1, k2):
    """Theorem 7.3: Tr_2 X = sum_j K_j X K_j^* with sum K_j^* K_j = 1 (a channel) and sum_j K_j^* A K_j = A (x) 1."""
    rng = seed(3 * len(k1) + len(k2))
    a1, a2 = mk(k1), mk(k2)
    a12 = tensor_product(a1, a2)
    for keep, af, other in ((1, a1, a2), (2, a2, a1)):
        P = composite.partial_trace_channel(a12, keep=keep)
        assert P.is_trace_preserving()
        X = random_op(a12, rng)
        assert (P(X).matrix - partial_trace(X, keep=keep).matrix).abs().max() < 1e-12
        A = random_op(af, rng)
        emb = kron(A, other.identity(), a12) if keep == 1 else kron(other.identity(), A, a12)
        assert (P.adjoint('blunt')(A).matrix - emb.matrix).abs().max() < 1e-12
    rho = states.random_density_matrix(a12)
    assert abs(composite.partial_trace_channel(a12)(rho).trace.real.item() - 1) < 1e-12


@pytest.mark.parametrize('k1', [[1, 2], [3], [2, 2, 1]])
@pytest.mark.parametrize('m', [1, 2, 3])
def test_prop_7_4_partial_trace_conditional_expectation(k1, m):
    """Proposition 7.4: E = Tr_2/m is the trace-preserving conditional expectation onto M1 (x) 1 for every trace, has
    the quasi-basis sqrt(m)(1 (x) e_ab), index m^2 and inclusion matrix m 1."""
    rng = seed(len(k1) + 17 * m)
    a1 = mk(k1)
    inc = inclusions.from_tensor_product(a1, m)
    M = inc.M
    assert np.array_equal(inc.matrix, m * np.eye(len(k1), dtype=int))
    assert inc.component_indices() == [float(m * m)] * len(k1)
    X, Y = random_op(M, rng), random_op(M, rng)
    A, A2 = random_op(a1, rng), random_op(a1, rng)
    one = M.tensor_factors[1].identity()
    for trace in ('blunt', 'norm', 'tau_vN'):
        E = inc.conditional_expectation(trace)
        assert (E(X).matrix - partial_trace(X).matrix / m).abs().max() < 1e-12
        tr = {'blunt': 'Tr_blunt', 'norm': 'Tr_norm', 'tau_vN': 'tau_vN'}[trace]
        assert abs(getattr(inc.embed(E(X)), tr)().item() - getattr(X, tr)().item()) < 1e-10
        assert (E.index().matrix - M.identity().matrix * m ** 2).abs().max() < 1e-12
    E = inc.conditional_expectation('blunt')
    bim = E(kron(A, one, M) @ X @ kron(A2, one, M))
    assert (bim.matrix - (A @ E(X) @ A2).matrix).abs().max() < 1e-10               # bimodule
    assert (E(kron(A, one, M)).matrix - A.matrix).abs().max() < 1e-12
    pos = E(unflat(M, rand_state(rng, sum(M.k_factors))))
    assert min(np.linalg.eigvalsh(B).min() for B in blocks_of(pos)) > -1e-12
    u = []
    for a in range(m):
        for b in range(m):
            e = np.zeros((m, m), complex)
            e[a, b] = math.sqrt(m)
            u.append(kron(a1.identity(), M.tensor_factors[1].from_blocks([torch.from_numpy(e)]), M))
    rec = sum((v @ inc.embed(E(v.adjoint() @ Y))).matrix for v in u)
    assert (rec - Y.matrix).abs().max() < 1e-10                                   # (b)
    assert (sum((v @ v.adjoint()).matrix for v in u) - m * m * M.identity().matrix).abs().max() < 1e-12
    Q = E.quasi_basis()[0]                                                       # (c): another quasi-basis, same index
    assert (sum(Q[j] @ Q[j].mH for j in range(Q.shape[0])) - m * m * M.identity().matrix[0]).abs().max() < 1e-12
    if m >= 2:                                   # uniqueness: a non-tracial slice is a conditional expectation that
        w = np.diag(np.linspace(1, 2, m))        # does not preserve the trace
        w = w / np.trace(w)
        R = full12(X)
        D1 = sum(k1)
        F0 = np.einsum('iajb,ba->ij', R.reshape(D1, m, D1, m), w)
        Ftr = np.trace(np.kron(F0, np.eye(m)))
        assert abs(Ftr - X.Tr_blunt().item()) > 1e-3


@pytest.mark.parametrize('L1,L2', [(2, 3), (1, 4), (3, 3), (2, 2)])
def test_prop_7_5_hard_core_bosons(L1, L2):
    """Proposition 7.5 and the hard-core boson example: the fused product of the number-conserving algebras of L1 and
    L2 sites is {N}' of L1+L2 sites (n_q = binom(L1+L2, q), Vandermonde), contains the plain product, and commutes
    with N."""
    rng = seed(L1 * 10 + L2)
    b1 = mk([math.comb(L1, q) for q in range(L1 + 1)], charges=list(range(L1 + 1)))
    b2 = mk([math.comb(L2, q) for q in range(L2 + 1)], charges=list(range(L2 + 1)))
    F, P = fused_tensor_product(b1, b2), tensor_product(b1, b2)
    L = L1 + L2
    assert F.k_factors == [math.comb(L, q) for q in range(L + 1)] and F.charges == list(range(L + 1))
    assert sum(k * k for k in F.k_factors) == math.comb(2 * L, L)
    assert sum(k * k for k in P.k_factors) == (sum(math.comb(L1, q) ** 2 for q in range(L1 + 1))
                                                 * sum(math.comb(L2, q) ** 2 for q in range(L2 + 1)))
    N1 = full(b1.central([float(q) for q in range(L1 + 1)]))
    N2 = full(b2.central([float(q) for q in range(L2 + 1)]))
    N = np.kron(N1, np.eye(2 ** L2)) + np.kron(np.eye(2 ** L1), N2)
    X = full12(random_op(F, rng))
    assert np.abs(X @ N - N @ X).max() < 1e-12
    A, B = random_op(b1, rng), random_op(b2, rng)
    assert np.abs(full12(kron(A, B, F)) - full12(kron(A, B, P))).max() < 1e-12   # (b): plain product inside
    Pi = [np.diag((np.diag(N).real == q).astype(float)) for q in range(L + 1)]
    for q, Pq in enumerate(Pi):                                                   # (c): Pi_q are polynomials in N
        poly = np.eye(N.shape[0])
        for qq in range(L + 1):
            if qq != q:
                poly = poly @ (N - qq * np.eye(N.shape[0])) / (q - qq)
        assert np.abs(poly - Pq).max() < 1e-9


def test_prop_7_5_equality_iff_one_pair_per_charge():
    """Proposition 7.5(b): fused = plain iff every total charge contains one pair."""
    a = mk([1, 2], charges=[0, 5])
    b = mk([2, 3], charges=[0, 1])
    assert fused_tensor_product(a, b).k_factors == tensor_product(a, b).k_factors
    c = mk([2, 3], charges=[0, 5])
    assert len(fused_tensor_product(a, c).k_factors) < len(tensor_product(a, c).k_factors)
    assert sum(k * k for k in fused_tensor_product(a, c).k_factors) > sum(k * k for k in tensor_product(a, c).k_factors)


@pytest.mark.parametrize('rule', ['add', 'parity', 'const'])
def test_thm_7_6_reduced_state_of_fused_product(rule):
    """Theorem 7.6: the restriction of a state of the fused product to M1 has the density E_1(Tr_H2 rho); E_1 removes
    coherences exactly when two charges of M1 fuse identically with one charge of M2."""
    rng = seed({'add': 1, 'parity': 2, 'const': 3}[rule])
    fuse = {'add': lambda q, r: q + r, 'parity': lambda q, r: (q + r) % 2, 'const': lambda q, r: 0}[rule]
    b1 = mk([1, 2, 1], charges=[0, 1, 2])
    b2 = mk([1, 3, 3, 1], charges=[0, 1, 2, 3])
    G = fused_tensor_product(b1, b2, fuse=fuse)
    rho = states.random_density_matrix(G, batch_size=2)
    for b in range(2):
        R = full12(rho, b).reshape(4, 8, 4, 8)
        for keep, XH, af in ((1, np.einsum('ijkj->ik', R), b1), (2, np.einsum('ijil->jl', R), b2)):
            mask = sla.block_diag(*[np.ones((k, k)) for k in af.k_factors])
            EX = XH * mask
            assert np.abs(full(partial_trace(rho, keep=keep), b) - EX).max() < 1e-12
            assert np.abs(full(composite.partial_trace_channel(G, keep=keep)(rho), b) - EX).max() < 1e-12
            lost = np.abs(XH - EX).max()
            assert (lost > 1e-3) == (rule != 'add')
            A = random_op(af, rng)                                # omega(A (x) 1) = Tr(rho_1 A)
            AI = np.kron(full(A), np.eye(8)) if keep == 1 else np.kron(np.eye(4), full(A))
            val = np.trace(full12(rho, b) @ AI)
            assert abs(val - np.trace(EX @ full(A))) < 1e-12
    A = random_op(b1, rng)
    assert (composite.partial_trace_channel(G).adjoint('blunt')(A).matrix
            - kron(A, b2.identity(), G).matrix).abs().max() < 1e-12


def test_thm_7_6_coherence_example():
    """Theorem 7.6 (proof): psi = (x + x') (x) y / sqrt 2 in one H_q has Tr_H2 |psi><psi| = |x+x'><x+x'|/2, whose
    off-diagonal block |x><x'|/2 is removed by E_1."""
    b1 = mk([1, 1], charges=[0, 2])
    b2 = mk([1], charges=[0])
    G = fused_tensor_product(b1, b2, fuse=lambda q, r: (q + r) % 2)
    assert G.k_factors == [2]
    rho = DensityMatrix(G, matrix=torch.full((1, 1, 2, 2), 0.5, dtype=C128))
    r1 = partial_trace(rho)
    assert np.abs(full(r1) - 0.5 * np.eye(2)).max() < 1e-15


@pytest.mark.parametrize('a,b', [(2, 2), (2, 5), (3, 4), (4, 3)])
def test_lemma_7_7_schmidt(a, b):
    """Lemma 7.7: both reduced states of a pure state have the nonzero eigenvalues s_i^2; E(psi) <= log min(a,b) with
    equality for maximally entangled states; product states have E = 0."""
    rng = seed(a * b)
    A12 = tensor_product(mk([a]), mk([b]))
    for kind in ('random', 'product', 'maximal'):
        if kind == 'random':
            psi = cplx(rng, a * b)
        elif kind == 'product':
            psi = np.kron(cplx(rng, a), cplx(rng, b))
        else:
            m = min(a, b)
            U, V = sla.qr(cplx(rng, a, a))[0], sla.qr(cplx(rng, b, b))[0]
            psi = sum(np.kron(U[:, i], V[:, i]) for i in range(m))
        psi = psi / np.linalg.norm(psi)
        rho = DensityMatrix(A12, matrix=torch.from_numpy(np.outer(psi, psi.conj()))[None, None])
        r1, r2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
        s2 = np.linalg.svd(psi.reshape(a, b), compute_uv=False) ** 2
        w1, w2 = np.linalg.eigvalsh(blocks_of(r1)[0]), np.linalg.eigvalsh(blocks_of(r2)[0])
        nz = lambda w: np.sort(w[w > 1e-12])
        assert np.allclose(nz(w1), nz(s2), atol=1e-12) and np.allclose(nz(w2), nz(s2), atol=1e-12)
        E = entropy(s2)
        assert abs(r1.von_neumann_entropy().item() - E) < 1e-10 and abs(r2.von_neumann_entropy().item() - E) < 1e-10
        assert E <= math.log(min(a, b)) + 1e-12
        if kind == 'product':
            assert E < 1e-10
        if kind == 'maximal':
            assert abs(E - math.log(min(a, b))) < 1e-10


@pytest.mark.parametrize('L1,L2,q', [(2, 3, 2), (3, 3, 3), (2, 2, 1)])
def test_prop_7_8_symmetry_resolved_entanglement(L1, L2, q):
    """Proposition 7.8: for a pure state of total charge q, E(psi) = H(p) + sum_c p_c E(psi_c), the reduced state lies
    in M1 with sector weights p_c, and M2 has the weight p_c at d(c) = q - c."""
    rng = seed(L1 + L2 + q)
    b1 = mk([math.comb(L1, x) for x in range(L1 + 1)], charges=list(range(L1 + 1)))
    b2 = mk([math.comb(L2, x) for x in range(L2 + 1)], charges=list(range(L2 + 1)))
    F = fused_tensor_product(b1, b2)
    s = F.charges.index(q)
    n = F.k_factors[s]
    psi = cplx(rng, n)
    psi /= np.linalg.norm(psi)
    mat = torch.zeros(1, F.C, F.k_max, F.k_max, dtype=C128)
    mat[0, s, :n, :n] = torch.from_numpy(np.outer(psi, psi.conj()))
    rho = DensityMatrix(F, matrix=mat)
    r1, r2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
    psi_H = np.zeros(2 ** L1 * 2 ** L2, complex)
    psi_H[kron_rows(F)[s]] = psi
    lam = np.linalg.svd(psi_H.reshape(2 ** L1, 2 ** L2), compute_uv=False) ** 2
    E = entropy(lam)
    assert abs(r1.von_neumann_entropy().item() - E) < 1e-10
    p = r1.sector_probabilities()[0].numpy()
    p2 = r2.sector_probabilities()[0].numpy()
    for c in range(L1 + 1):
        d = q - c
        assert abs(p2[d] - p[c]) < 1e-12 if 0 <= d <= L2 else p[c] < 1e-14
    conf = sum(p[c] * entropy(np.linalg.eigvalsh(B / p[c])) for c, B in enumerate(blocks_of(r1)) if p[c] > 1e-14)
    assert abs(E - entropy(p) - conf) < 1e-10
    # Tr_H2 |psi><psi| has no blocks between sectors of M1 (it commutes with N1)
    R = np.outer(psi_H, psi_H.conj()).reshape(2 ** L1, 2 ** L2, 2 ** L1, 2 ** L2)
    XH = np.einsum('ijkj->ik', R)
    N1 = np.diag(full(b1.central([float(x) for x in range(L1 + 1)]))).real
    assert np.abs(XH * (N1[:, None] != N1[None, :])).max() < 1e-14


@pytest.mark.parametrize('a,b', [(2, 2), (2, 3), (3, 5)])
def test_lemma_7_9_7_10_gaussian_moments(a, b):
    """Lemmas 7.9 and 7.10: E||g||^4 = N(N+1), E Tr(GG^*)^2 = ab(a+b), and ||g|| is uncorrelated with functions of
    g/||g|| (independence)."""
    rng = seed(a * 100 + b)
    B = 80_000
    G = (rng.standard_normal((B, a, b)) + 1j * rng.standard_normal((B, a, b))) / math.sqrt(2)
    n2 = (np.abs(G) ** 2).sum((1, 2))
    W = G @ G.conj().transpose(0, 2, 1)
    trW2 = np.einsum('bij,bji->b', W, W).real
    assert zscore(n2 ** 2, a * b * (a * b + 1)) < Z
    assert zscore(trW2, a * b * (a + b)) < Z
    purity = trW2 / n2 ** 2
    cov = (n2 - n2.mean()) * (purity - purity.mean())
    assert zscore(cov, 0) < Z
    assert zscore(np.abs(G[:, 0, 0]) ** 4, 2) < Z                 # E|g_i|^4 = 2


def page(a, b):
    a, b = min(a, b), max(a, b)
    return sum(1 / k for k in range(b + 1, a * b + 1)) - (a - 1) / (2 * b)


@pytest.mark.parametrize('a,b', [(2, 2), (2, 4), (3, 3), (4, 2)])
def test_thm_7_11_7_13_lubkin_and_page(a, b):
    """Theorems 7.11 and 7.13 (with a > b by symmetry): E Tr rho_A^2 = (a+b)/(ab+1) and E S(rho_A) = Page's formula,
    for Haar-random pure states reduced with partial_trace; Corollary 7.12 bounds the mean."""
    seed(a * 31 + b)
    B = 30_000
    g = torch.randn(B, a * b, dtype=C128)
    psi = g / g.norm(dim=-1, keepdim=True)
    A12 = tensor_product(mk([a]), mk([b]))
    pure = DensityMatrix(A12, matrix=(psi[:, :, None] * psi.conj()[:, None, :])[:, None], validate=False)
    rA = partial_trace(pure, keep=1)
    rB = partial_trace(pure, keep=2)
    pur, S = rA.purity().numpy(), rA.von_neumann_entropy().numpy()
    assert np.abs(S - rB.von_neumann_entropy().numpy()).max() < 1e-9
    assert zscore(pur, (a + b) / (a * b + 1)) < Z
    assert zscore(S, page(a, b)) < Z
    assert math.log((a * b + 1) / (a + b)) <= page(a, b) <= math.log(min(a, b))
    assert (S >= -np.log(pur) - 1e-10).all()                       # S >= -log Tr rho^2 (proof of Cor. 7.12)


def test_page_formula_exact_bounds():
    """Theorem 7.13 / Corollary C.18: log a - (a-1)/(ab+1) - (a-1)/(2b) <= E S <= log a - (a-1)/(2b), the
    1/(2b) bound, two qubits 1/3, Corollary 7.12, and the gap from log min(a,b) is less than one nat."""
    assert abs(page(2, 2) - 1 / 3) < 1e-15
    for a in range(1, 12):
        for b in range(a, 40):
            e = page(a, b)
            assert math.log(a) - (a - 1) / (a * b + 1) - (a - 1) / (2 * b) - 1e-12 <= e <= math.log(a) - (a - 1) / (2 * b) + 1e-12
            assert abs(e - math.log(a) + a / (2 * b)) <= 1 / (2 * b) + 1e-12
            assert math.log((a * b + 1) / (a + b)) - 1e-12 <= e <= math.log(a) + 1e-12
            assert 0 <= math.log(a) - e < 1


@pytest.mark.parametrize('ks,r', [([1, 2, 3], 2), ([2], 3), ([2, 2], 1), ([1, 4], 5)])
def test_prop_7_14_induced_measure_purity(ks, r):
    """Proposition 7.14: E Tr rho^2 = sum_c k_c(k_c+r) / (D(Dr+1)) for random_density_matrix(alg, rank=r)."""
    seed(sum(ks) * 7 + r)
    alg = mk(ks)
    D = sum(ks)
    pur = states.random_density_matrix(alg, 40_000, rank=r).purity().numpy()
    exact = sum(k * (k + r) for k in ks) / (D * (D * r + 1))
    if r == 1 and len(ks) == 1:
        assert np.abs(pur - 1).max() < 1e-12
    else:
        assert zscore(pur, exact) < Z


def test_prop_7_14_is_restriction_of_haar_state():
    """Section 7.5: random_density_matrix(rank=r) has the law of E_1(Tr_2 |psi><psi|) for Haar psi in H1 (x) C^r
    (compared through the mean of p_c and of Tr rho^2)."""
    seed(21)
    ks, r, B = [1, 2], 2, 40_000
    alg = mk(ks)
    G = torch.randn(B, sum(ks), r, dtype=C128)
    G = G / G.flatten(1).norm(dim=-1)[:, None, None]
    rho_h = (G @ G.mH).numpy()
    p_h = np.stack([np.trace(rho_h[:, :1, :1], axis1=1, axis2=2).real,
                    np.trace(rho_h[:, 1:, 1:], axis1=1, axis2=2).real], 1)
    pur_h = (np.abs(rho_h[:, :1, :1]) ** 2).sum((1, 2)) + (np.abs(rho_h[:, 1:, 1:]) ** 2).sum((1, 2))
    lib = states.random_density_matrix(alg, B, rank=r)
    for x, y in ((p_h[:, 0], lib.sector_probabilities()[:, 0].numpy()), (pur_h, lib.purity().numpy())):
        se = math.sqrt(x.var() / B + y.var() / B)
        assert abs(x.mean() - y.mean()) < Z * se
    D = sum(ks)
    assert zscore(p_h[:, 0], 1 / D) < Z                          # E p_c = k_c / D


@pytest.mark.parametrize('trace', ['tau_vN', 'norm', 'blunt'])
def test_exercise_7_1_dual_of_partial_trace(trace):
    """Chapter 7, Exercise 1: the dual of Tr_2 is A (x) Z with Z = (+)_d C_2 l_d 1 for tau_vN, (+)_d l_d 1 for Tr_norm,
    1 for the blunt traces."""
    rng = seed(3)
    a1, a2 = mk([1, 2]), mk([2, 3])
    a12 = tensor_product(a1, a2)
    A = random_op(a1, rng)
    Zs = {'tau_vN': [a2.C * l for l in a2.k_factors], 'norm': [float(l) for l in a2.k_factors], 'blunt': [1.0, 1.0]}
    dual = composite.partial_trace_channel(a12, keep=1).adjoint(trace)
    assert (dual(A).matrix - kron(A, a2.central(Zs[trace]), a12).matrix).abs().max() < 1e-12
    if trace == 'tau_vN':
        assert Zs[trace] == [4, 6]


@pytest.mark.parametrize('k1,k2', [([1, 2], [2, 3]), ([2], [2]), ([3, 1], [2, 1])])
def test_exercise_7_2_subadditivity(k1, k2):
    """Chapter 7, Exercise 2: S(rho_1) + S(rho_2) - S(rho) = S(rho || rho_1 (x) rho_2) >= 0, zero for product states."""
    seed(len(k1) + len(k2))
    a1, a2 = mk(k1), mk(k2)
    a12 = tensor_product(a1, a2)
    rho = states.random_density_matrix(a12, batch_size=4)
    r1, r2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
    mutual = r1.von_neumann_entropy() + r2.von_neumann_entropy() - rho.von_neumann_entropy()
    assert (mutual > 0).all()
    assert torch.allclose(mutual, states.relative_entropy(rho, kron(r1, r2, a12)), atol=1e-10)
    s1, s2 = states.random_density_matrix(a1), states.random_density_matrix(a2)
    prod = kron(s1, s2, a12)
    mi = (partial_trace(prod, 1).von_neumann_entropy() + partial_trace(prod, 2).von_neumann_entropy()
          - prod.von_neumann_entropy())
    assert abs(mi.item()) < 1e-10


def test_exercise_7_3_one_particle_on_two_sites():
    """Chapter 7, Exercise 3: (|01>+|10>)/sqrt2 has E = H(p) = log 2 with zero configurational part, and
    (|00>+|11>)/sqrt2 restricts to the mixture of |00> and |11>."""
    site = mk([1, 1], charges=[0, 1])
    S2 = fused_tensor_product(site, site)
    assert S2.k_factors == [1, 2, 1]
    mat = torch.zeros(1, 3, 2, 2, dtype=C128)
    mat[0, 1] = 0.5
    r1 = partial_trace(DensityMatrix(S2, matrix=mat))
    assert torch.allclose(r1.sector_probabilities()[0], torch.tensor([0.5, 0.5], dtype=F64))
    assert abs(r1.von_neumann_entropy().item() - math.log(2)) < 1e-12
    psi = np.zeros(4)
    psi[[0, 3]] = 1 / math.sqrt(2)                                # (|00> + |11>)/sqrt 2 in Kronecker order
    P = np.outer(psi, psi)
    vals = [np.trace(P @ full12(op)) for op in [S2.from_blocks([torch.eye(1), torch.zeros(2, 2), torch.zeros(1, 1)]),
                                                S2.from_blocks([torch.zeros(1, 1), torch.zeros(2, 2), torch.eye(1)])]]
    assert np.allclose(vals, [0.5, 0.5])
    rng = seed(0)
    for _ in range(3):                                            # no element of the fused algebra sees the coherence
        X = full12(random_op(S2, rng))
        mixture = 0.5 * X[0, 0] + 0.5 * X[3, 3]
        assert abs(psi @ X @ psi - mixture) < 1e-12


# ================================================================================================================
# Chapter 8: dynamics
# ================================================================================================================
def test_lemma_8_1_exponential():
    """Lemma 8.1: d/dt e^{tB} = B e^{tB} = e^{tB} B, e^{(s+t)B} = e^{sB} e^{tB}, ||e^{tB}|| <= e^{|t| ||B||}."""
    rng = seed(1)
    Bm = cplx(rng, 5, 5)
    E = lambda t: sla.expm(t * Bm)
    for t, s in ((0.3, 0.7), (-0.4, 1.1)):
        assert np.abs(E(s + t) - E(s) @ E(t)).max() < 1e-10
        h = 1e-5
        deriv = (E(t + h) - E(t - h)) / (2 * h)
        assert np.abs(deriv - Bm @ E(t)).max() < 1e-6 and np.abs(Bm @ E(t) - E(t) @ Bm).max() < 1e-10
        assert np.linalg.norm(E(t), 2) <= math.exp(abs(t) * np.linalg.norm(Bm, 2)) + 1e-12
    series = sum(np.linalg.matrix_power(0.5 * Bm, n) / math.factorial(n) for n in range(40))
    assert np.abs(series - E(0.5)).max() < 1e-12


@pytest.mark.parametrize('ks', LAYOUTS)
def test_prop_8_2_propagator_and_solutions(ks):
    """Proposition 8.2: U(t) = e^{-itH} is a block-diagonal unitary group with U' = -iHU equal to the exponential
    series; rho(t) = U rho U^* solves the von Neumann equation; the Heisenberg picture."""
    rng = seed(sum(ks) + len(ks))
    alg = mk(ks)
    H = op_of(alg, [herm(rng, k) for k in ks])
    Hf = full(H)
    t, s = 0.6, 0.35
    U = dynamics.propagator(H, t)
    assert np.abs(full(U) - sla.expm(-1j * t * Hf)).max() < 1e-12
    assert np.abs(full(U).conj().T @ full(U) - np.eye(sum(ks))).max() < 1e-12
    assert np.abs(full(dynamics.propagator(H, s + t)) - full(dynamics.propagator(H, s)) @ full(U)).max() < 1e-12
    h = 1e-5
    dU = (full(dynamics.propagator(H, t + h)) - full(dynamics.propagator(H, t - h))) / (2 * h)
    assert np.abs(dU + 1j * Hf @ full(U)).max() < 1e-7
    rho0 = states.random_density_matrix(alg)
    times = [0.0, t, 1.3]
    traj = dynamics.von_neumann(rho0, H, times)
    for tt, r in zip(times, traj):
        Ut = sla.expm(-1j * tt * Hf)
        assert np.abs(full(r) - Ut @ full(rho0) @ Ut.conj().T).max() < 1e-12
        assert min(np.linalg.eigvalsh(B).min() for B in blocks_of(r)) > -1e-12
    rt = full(traj[1])
    dr = (full(dynamics.von_neumann(rho0, H, [t + h])[0]) - full(dynamics.von_neumann(rho0, H, [t - h])[0])) / (2 * h)
    assert np.abs(dr + 1j * (Hf @ rt - rt @ Hf)).max() < 1e-7
    A = op_of(alg, [herm(rng, k) for k in ks])
    At = U.adjoint() @ A @ U
    assert abs(rho0.expectation(At).item() - traj[1].expectation(A).item()) < 1e-12
    dA = (full(dynamics.propagator(H, t + h).adjoint() @ A @ dynamics.propagator(H, t + h))
          - full(dynamics.propagator(H, t - h).adjoint() @ A @ dynamics.propagator(H, t - h))) / (2 * h)
    assert np.abs(dA - 1j * (Hf @ full(At) - full(At) @ Hf)).max() < 1e-7


@pytest.mark.parametrize('ks', [[2, 3, 4], [1, 2], [3]])
def test_thm_8_3_conservation_of_charge(ks):
    """Theorem 8.3: ||P_c psi(t)|| and Tr(Z rho(t)) are constant for time-dependent H in M (exact and RK4); RK4 never
    moves amplitude between sectors; for constant H, Tr(A rho(t)) is constant if [A, H] = 0."""
    rng = seed(len(ks) * 3)
    alg = mk(ks)
    H = op_of(alg, [herm(rng, k) for k in ks])
    Vd = op_of(alg, [herm(rng, k) for k in ks])
    rows = torch.tensor([[1.0] * k + [0.0] * (alg.k_max - k) for k in ks], dtype=C128)
    psi0 = torch.randn(1, alg.C, alg.k_max, dtype=C128) * rows
    psi0 = psi0 / psi0.norm()
    times = [0.0, 0.5, 1.5]
    psi = dynamics.schrodinger(psi0, H, times)
    w = (psi.abs() ** 2).sum(-1)
    assert (w - w[0]).abs().max() < 1e-12
    Hot = lambda s: H.matrix + math.sin(2 * s) * Vd.matrix
    prk = dynamics.schrodinger_rk4(psi0, Hot, times, substeps=60, algebra=alg)
    assert ((prk.abs() ** 2).sum(-1) - w[0]).abs().max() < 1e-7
    if len(ks) > 1:                                               # one sector only: the others stay exactly zero
        one = psi0.clone()
        one[:, 1:] = 0
        one = one / one.norm()
        p1 = dynamics.schrodinger_rk4(one, Hot, times, substeps=10, algebra=alg)
        assert p1[:, :, 1:].abs().max() == 0
    rho0 = states.random_density_matrix(alg)
    Zc = alg.central([float(c + 1) ** 2 for c in range(len(ks))])
    for r in dynamics.von_neumann(rho0, H, times):
        assert abs(r.expectation(Zc).item() - rho0.expectation(Zc).item()) < 1e-12
    A = H.apply_function(lambda x: torch.cos(x) + x ** 3)        # commutes with H
    for r in dynamics.von_neumann(rho0, H, times):
        assert abs(r.expectation(A).item() - rho0.expectation(A).item()) < 1e-11


@pytest.mark.parametrize('x', [0.1, 0.9, 2.0, 2.8, 3.0])
def test_rk4_amplification_factor(x):
    """Section 8.1 (tutorial): one RK4 step multiplies an eigenvector by T(-ih lambda) with
    |T(ix)|^2 = 1 - x^6/72 + x^8/576 < 1 for 0 < |x| < 2 sqrt 2."""
    lam = 1.7
    h = x / lam
    y = dynamics.rk4(lambda t, v: -1j * lam * v, torch.ones(1, dtype=C128), [0.0, h], substeps=1)[-1]
    T = sum((-1j * x) ** n / math.factorial(n) for n in range(5))
    assert abs(y.item() - T) < 1e-14
    val = 1 - x ** 6 / 72 + x ** 8 / 576
    assert abs(abs(T) ** 2 - val) < 1e-12
    assert (val < 1) == (x < 2 * math.sqrt(2))


def test_rk4_does_not_preserve_positivity():
    """Tutorial 8.3: one RK4 step with h gamma = 3 maps diag(0,1) to diag(-3/8, 11/8) for spontaneous decay."""
    sm = torch.tensor([[0.0, 1.0], [0.0, 0.0]], dtype=C128)
    f = dynamics.lindblad_rhs(None, [sm[None, None]], rates=[3.0])
    out = dynamics.rk4(f, torch.diag(torch.tensor([0.0, 1.0], dtype=C128))[None, None], [0.0, 1.0], substeps=1)[-1]
    assert (out[0, 0] - torch.diag(torch.tensor([-3 / 8, 11 / 8], dtype=C128))).abs().max() < 1e-14


def dyson_terms(G0, Phi_s, t, nmax, M=48):
    """S_n(t), n <= nmax, as Taylor coefficients of z -> exp(t (L0 + z Phi)) (discrete Fourier transform)."""
    zs = np.exp(2j * np.pi * np.arange(M) / M)
    vals = np.array([sla.expm(t * (G0 + z * Phi_s)) for z in zs])
    return [np.tensordot(zs ** (-n), vals, axes=(0, 0)) / M for n in range(nmax + 1)]


@pytest.mark.parametrize('k,nj', [(2, 1), (3, 2), (4, 3)])
def test_thm_8_6_gksl_semigroup_and_dyson(k, nj):
    """Theorem 8.6: e^{tL} is CP, TP and a semigroup; S_0 has the Kraus operator e^{-itH_eff}; every Dyson term S_n is
    CP, sum_n S_n = e^{tL}, and ||S_n|| <= e^{t||L0||} (t||Phi||)^n / n!."""
    rng = seed(k * 10 + nj)
    alg = mk([k])
    H = herm(rng, k)
    Ls = [cplx(rng, k, k) / 2 for _ in range(nj)]
    to = lambda X: alg.from_blocks([torch.from_numpy(X)])
    S = lambda t: dynamics.lindblad_channel(alg, to(H), [to(L) for L in Ls], t=t).superoperator()[0][0].numpy()
    G = dynamics.lindblad_superoperator(alg, to(H), [to(L) for L in Ls])[0][0].numpy()
    assert np.abs(G - superop(lambda R: gksl_np(H, Ls, R), k)).max() < 1e-12
    for t in (0.2, 1.0):
        assert np.abs(S(t) - sla.expm(t * G)).max() < 1e-10
        assert np.linalg.eigvalsh(choi_from_superop(S(t), k)).min() > -1e-12
    assert np.abs(S(0.3) @ S(0.5) - S(0.8)).max() < 1e-10
    one = np.eye(k).reshape(-1)
    assert np.abs(one @ G).max() < 1e-12                               # Tr L(rho) = 0
    Heff = H - 0.5j * sum(L.conj().T @ L for L in Ls)
    G0 = superop(lambda R: -1j * (Heff @ R - R @ Heff.conj().T), k)
    Phi_s = superop(lambda R: sum(L @ R @ L.conj().T for L in Ls), k)
    assert np.abs(G0 + Phi_s - G).max() < 1e-12
    t = 0.9
    E0 = sla.expm(-1j * t * Heff)
    assert np.abs(sla.expm(t * G0) - superop(lambda R: E0 @ R @ E0.conj().T, k)).max() < 1e-12
    terms = dyson_terms(G0, Phi_s, t, 8)
    assert np.abs(sum(dyson_terms(G0, Phi_s, t, 30)) - sla.expm(t * G)).max() < 1e-10
    a, b = np.linalg.norm(G0, 2), np.linalg.norm(Phi_s, 2)
    for n, Sn in enumerate(terms):
        assert np.linalg.eigvalsh(choi_from_superop(Sn, k)).min() > -1e-10           # CP
        assert np.linalg.norm(Sn, 2) <= math.exp(t * a) * (t * b) ** n / math.factorial(n) + 1e-10
    jump = sla.expm(t * G) - terms[0]
    assert np.linalg.eigvalsh(choi_from_superop(jump, k)).min() > -1e-10
    # S_1(t) = int_0^t S_0(t-s) Phi S_0(s) ds by Gauss-Legendre quadrature
    xs, ws = np.polynomial.legendre.leggauss(30)
    ss = (xs + 1) * t / 2
    S1 = sum(w * t / 2 * sla.expm((t - s_) * G0) @ Phi_s @ sla.expm(s_ * G0) for w, s_ in zip(ws, ss))
    assert np.abs(S1 - terms[1]).max() < 1e-10


def test_lemma_8_5_generator_of_a_semigroup():
    """Lemma 8.5: L = (T_h - id)/h V_h^{-1} with V_h = (1/h) int_0^h T_s ds, and (T_t - id)/t -> L."""
    rng = seed(4)
    k = 3
    alg = mk([k])
    H, Ls = herm(rng, k), [cplx(rng, k, k) / 2 for _ in range(2)]
    to = lambda X: alg.from_blocks([torch.from_numpy(X)])
    G = dynamics.lindblad_superoperator(alg, to(H), [to(L) for L in Ls])[0][0].numpy()
    T = lambda t: dynamics.lindblad_channel(alg, to(H), [to(L) for L in Ls], t=t).superoperator()[0][0].numpy()
    h = 0.4
    xs, ws = np.polynomial.legendre.leggauss(20)
    Vh = sum(w / 2 * T((x + 1) * h / 2) for w, x in zip(ws, xs))
    assert np.abs((T(h) - np.eye(k * k)) / h @ np.linalg.inv(Vh) - G).max() < 1e-9
    errs = [np.abs((T(t) - np.eye(k * k)) / t - G).max() for t in (1e-2, 1e-3)]
    assert errs[1] < errs[0] / 5


def sector_full_generator(alg, H, inside, between, rates_in, rates_between):
    """GKSL generator on M_D (Proposition 8.8(1)): jumps sqrt(g) L_j and sqrt(g) J^{dc}_l."""
    ks = alg.k_factors
    o = np.cumsum([0] + ks)
    D = sum(ks)
    jumps = [math.sqrt(g) * full(L) for g, L in zip(rates_in, inside)]
    for g, Phi in zip(rates_between, between):
        K = Phi.kraus[0].numpy()
        for l in range(K.shape[0]):
            for d, kd in enumerate(ks):
                for c, kc in enumerate(ks):
                    J = np.zeros((D, D), complex)
                    J[o[d]:o[d] + kd, o[c]:o[c] + kc] = K[l, d, c, :kd, :kc]
                    if np.abs(J).max() > 0:
                        jumps.append(math.sqrt(g) * J)
    Hf = full(H) if H is not None else np.zeros((D, D))
    return Hf, jumps


@pytest.mark.parametrize('ks', [[1, 2, 3], [2, 2], [3]])
def test_prop_8_8_compression_of_jumps_between_sectors(ks):
    """Proposition 8.8: lindblad_rhs is the restriction of the GKSL generator with jumps J^{dc} (1), equals E of the
    generator with summed jumps L_l = sum J^{dc}_l (2), and maps M into M; e^{tL} on M is a channel."""
    rng = seed(sum(ks) * 3)
    alg = mk(ks)
    D = sum(ks)
    H = op_of(alg, [herm(rng, k) for k in ks])
    inside = [random_op(alg, rng)]
    Phi = channels.random_inter_sector_channel(alg, alg, kraus_rank=2)
    gi, gb = [0.7], [1.3]
    rho = states.random_density_matrix(alg)
    on_M = full(alg.operator(dynamics.lindblad_rhs(H, inside + [Phi], rates=gi + gb, algebra=alg)(0.0, rho.matrix)))
    Hf, jumps = sector_full_generator(alg, H, inside, [Phi], gi, gb)
    X = full(rho)
    assert np.abs(gksl_np(Hf, jumps, X) - on_M).max() < 1e-12                          # (1)
    mask = sla.block_diag(*[np.ones((k, k)) for k in ks])
    K = Phi.kraus[0].numpy()
    o = np.cumsum([0] + ks)
    summed = []
    for l in range(K.shape[0]):
        Lsum = np.zeros((D, D), complex)
        for d, kd in enumerate(ks):
            for c, kc in enumerate(ks):
                Lsum[o[d]:o[d] + kd, o[c]:o[c] + kc] = K[l, d, c, :kd, :kc]
        summed.append(math.sqrt(gb[0]) * Lsum)
    tilde = gksl_np(Hf, [math.sqrt(gi[0]) * full(inside[0])] + summed, X)
    assert np.abs(tilde * mask - on_M).max() < 1e-12                                   # (2)
    if len(ks) > 1:
        assert np.abs(tilde * (1 - mask)).max() > 1e-3                                  # L~(rho) not in M
    # (1): exp(t L^) preserves block-diagonal matrices and is a channel on M; compare with RK4 of lindblad_evolve
    Gfull = superop(lambda R: gksl_np(Hf, jumps, R), D)
    E1 = sla.expm(0.7 * Gfull)
    out = (E1 @ X.reshape(-1)).reshape(D, D)
    assert np.abs(out * (1 - mask)).max() < 1e-12
    rk = dynamics.lindblad_evolve(rho, H, inside + [Phi], [0.0, 0.7], rates=gi + gb, substeps=60)[-1]
    assert np.abs(full(rk) - out).max() < 1e-7
    assert np.linalg.eigvalsh(choi_from_superop(E1, D)).min() > -1e-10


def fock_algebra(m, Nmax):
    """Number sectors N = 0..Nmax of m bosonic modes; bases, annihilation operators a_j restricted to N -> N-1."""
    import itertools
    bases = [[occ for occ in itertools.product(range(N + 1), repeat=m) if sum(occ) == N] for N in range(Nmax + 1)]
    alg = mk([len(B) for B in bases], charges=list(range(Nmax + 1)))
    ann = []
    for j in range(m):
        blocks = {}
        for N in range(1, Nmax + 1):
            K = np.zeros((len(bases[N - 1]), len(bases[N])))
            for col, occ in enumerate(bases[N]):
                if occ[j] > 0:
                    new = list(occ)
                    new[j] -= 1
                    K[bases[N - 1].index(tuple(new)), col] = math.sqrt(occ[j])
            blocks[(N - 1, N)] = K
        ann.append(blocks)
    return alg, bases, ann


def test_prop_8_8_3_particle_loss_needs_no_compression():
    """Proposition 8.8(3): with the annihilation operators as jumps the GKSL generator on M_D already maps M into M
    and equals lindblad_rhs."""
    rng = seed(7)
    alg, bases, ann = fock_algebra(2, 2)
    ks = alg.k_factors
    D = sum(ks)
    o = np.cumsum([0] + ks)
    loss = InterSectorChannel.from_blocks(alg, alg, {dc: [torch.from_numpy(a[dc]) for a in ann] for dc in ann[0]})
    rho = states.random_density_matrix(alg)
    on_M = full(alg.operator(dynamics.lindblad_rhs(None, [loss], algebra=alg)(0.0, rho.matrix)))
    Ls = []
    for a in ann:
        L = np.zeros((D, D))
        for (d, c), K in a.items():
            L[o[d]:o[d] + ks[d], o[c]:o[c] + ks[c]] = K
        Ls.append(L)
    assert np.abs(gksl_np(np.zeros((D, D)), Ls, full(rho)) - on_M).max() < 1e-12
    del rng


@pytest.mark.parametrize('m,N0', [(2, 2), (3, 2), (2, 3)])
def test_prop_8_9_binomial_particle_loss(m, N0):
    """Proposition 8.9(2) and Chapter 8, Exercise 3: N0 bosons in m modes, each lost at rate gamma: Phi^*(1) = N,
    Q = R - diag(colsums R), and p_N(t) = binom(N0,N) q^N (1-q)^(N0-N) for any Hamiltonian in M."""
    rng = seed(m * 10 + N0)
    alg, bases, ann = fock_algebra(m, N0)
    ks = alg.k_factors
    loss = InterSectorChannel.from_blocks(alg, alg, {dc: [torch.from_numpy(a[dc]) for a in ann] for dc in ann[0]})
    assert (loss.adjoint()(alg.identity()).matrix - alg.central([float(N) for N in range(N0 + 1)]).matrix).abs().max() < 1e-12
    R = loss.transition_matrix()[0].numpy()
    Q = R - np.diag(R.sum(0))
    Qref = np.zeros((N0 + 1, N0 + 1))
    for N in range(1, N0 + 1):
        Qref[N - 1, N], Qref[N, N] = N, -N
    assert np.abs(Q - Qref).max() < 1e-12
    gamma = 0.45
    H = op_of(alg, [herm(rng, k) for k in ks])
    m0 = torch.zeros(1, alg.C, alg.k_max, alg.k_max, dtype=C128)
    m0[0, N0, :ks[N0], :ks[N0]] = torch.from_numpy(rand_state(rng, ks[N0], rank=1))
    ts = [0.0, 0.6, 1.5]
    traj = dynamics.lindblad_evolve(DensityMatrix(alg, matrix=m0), H, [loss], ts, rates=[gamma], substeps=30)
    for t, r in zip(ts, traj):
        q = math.exp(-gamma * t)
        p = r.sector_probabilities()[0].numpy()
        binom = np.array([math.comb(N0, N) * q ** N * (1 - q) ** (N0 - N) for N in range(N0 + 1)])
        assert np.abs(p - binom).max() < 1e-8
        assert np.abs(p - sla.expm(gamma * t * Q)[:, N0]).max() < 1e-8


def test_prop_8_9_rate_equation_independent_of_states_inside():
    """Proposition 8.9: (1) no jumps between sectors: p constant; (2) if sum_l K^{dc*} K^{dc} = r_dc 1, then
    p(t) = e^{tQ} p(0) for every H, every jump inside and every state with the same p(0)."""
    rng = seed(13)
    ks = [2, 3, 2]
    alg = mk(ks)
    r = {(1, 0): 0.7, (0, 1): 0.3, (2, 1): 0.5, (1, 2): 1.1, (2, 0): 0.2}
    blocks = {}
    for (d, c), rate in r.items():                                 # K^{dc}_l = sqrt(r) blocks of an isometry
        V = sla.qr(cplx(rng, 2 * ks[d], ks[c]), mode='economic')[0].reshape(2, ks[d], ks[c])
        blocks[(d, c)] = [torch.from_numpy(math.sqrt(rate) * V[l]) for l in range(2)]
    blocks[(1, 1)] = [torch.from_numpy(cplx(rng, 3, 3) / 3)]      # a jump inside a sector, also allowed
    Phi = InterSectorChannel.from_blocks(alg, alg, blocks)
    Q = np.zeros((3, 3))
    for (d, c), rate in r.items():
        Q[d, c] += rate
        Q[c, c] -= rate
    p0 = np.array([0.5, 0.2, 0.3])
    ts = [0.0, 0.8]
    results = []
    for trial in range(2):
        H = op_of(alg, [herm(rng, k) for k in ks])
        inside = [random_op(alg, rng)]
        rho = op_of(alg, [p * rand_state(rng, k) for p, k in zip(p0, ks)])
        rho = DensityMatrix(alg, matrix=rho.matrix)
        traj = dynamics.lindblad_evolve(rho, H, inside + [Phi], ts, rates=[0.6, 1.0], substeps=80)
        results.append(traj[-1].sector_probabilities()[0].numpy())
        assert np.abs(results[-1] - sla.expm(0.8 * Q) @ p0).max() < 1e-8
        traj0 = dynamics.lindblad_evolve(rho, H, inside, ts, rates=[0.6], substeps=40)     # (1)
        assert np.abs(traj0[-1].sector_probabilities()[0].numpy() - p0).max() < 1e-12
    R = Phi.transition_matrix()[0].numpy()
    Qlib = R - np.diag(R.sum(0))
    assert np.abs(Qlib - Q).max() < 1e-12


def test_lueders_update():
    """Section 8.4: the Lueders update P rho P / Tr(P rho P) with probability Tr(P rho) = Tr(P rho P)."""
    seed(3)
    alg = mk([2, 3])
    rho = states.random_density_matrix(alg)
    P = alg.central([0.0, 1.0])
    post, prob = states.lueders_update(rho, P)
    assert abs(prob.item() - rho.sector_probabilities()[0, 1].item()) < 1e-12
    assert (post.matrix - (P @ rho @ P).matrix / prob.item()).abs().max() < 1e-12
    assert abs(post.trace.real.item() - 1) < 1e-12


@pytest.mark.parametrize('g1,g2', [(0.6, 0.2), (0.5, 0.5)])
def test_prop_8_10_8_11_no_jump_evolution(g1, g2):
    """Propositions 8.10 and 8.11: for losses from N0 = 2, S_n(t) rho is the sector N0 - n of rho(t); w_0 = p_N0 is
    nonincreasing with w_0' = -Tr(A S_0 rho) and w_0 = ||e^{-itH_eff} psi||^2; the post-selected state is the no-jump
    state (Dyson terms from the counting field on M_D)."""
    rng = seed(int(g1 * 10 + g2 * 100))
    alg, bases, ann = fock_algebra(2, 2)
    ks = alg.k_factors
    D = sum(ks)
    o = np.cumsum([0] + ks)
    J, U = 1.0, 0.8
    r2 = math.sqrt(2)
    Hb = [np.zeros((1, 1)), np.array([[0, J], [J, 0]]), np.array([[U, r2 * J, 0], [r2 * J, 0, r2 * J], [0, r2 * J, U]])]
    H = op_of(alg, Hb)
    loss = [InterSectorChannel.from_blocks(alg, alg, {dc: [torch.from_numpy(K)] for dc, K in a.items()}) for a in ann]
    psi = cplx(rng, 3)
    psi /= np.linalg.norm(psi)
    m0 = torch.zeros(1, 3, 3, 3, dtype=C128)
    m0[0, 2] = torch.from_numpy(np.outer(psi, psi.conj()))
    rho0 = DensityMatrix(alg, matrix=m0)
    Hf, jumps = sector_full_generator(alg, H, [], loss, [], [g1, g2])
    A = sum(L.conj().T @ L for L in jumps)
    Heff = Hf - 0.5j * A
    G0 = superop(lambda R: -1j * (Heff @ R - R @ Heff.conj().T), D)
    Phi_s = superop(lambda R: sum(L @ R @ L.conj().T for L in jumps), D)
    X0 = full(rho0)
    ts = [0.5, 1.0, 2.0]
    traj = dynamics.lindblad_evolve(rho0, H, loss, [0.0] + ts, rates=[g1, g2], substeps=400)[1:]   # RK4 error ~2e-11
    w0s = []
    for t, r in zip(ts, traj):
        terms = dyson_terms(G0, Phi_s, t, 3)
        rt = full(r)
        for n in range(3):
            Sn = (terms[n] @ X0.reshape(-1)).reshape(D, D)
            sec = 2 - n
            Pn = np.zeros((D, D))
            Pn[o[sec]:o[sec] + ks[sec], o[sec]:o[sec] + ks[sec]] = np.eye(ks[sec])
            assert np.abs(Sn - Pn @ rt @ Pn).max() < 1e-9                         # Prop. 8.11
            assert np.trace(Sn).real > -1e-12
        assert np.abs((terms[3] @ X0.reshape(-1))).max() < 1e-10                   # S_n = 0 for n > N0
        w0 = np.linalg.norm(sla.expm(-1j * t * Heff) @ np.concatenate([[0, 0, 0], psi])) ** 2
        S0 = (terms[0] @ X0.reshape(-1)).reshape(D, D)
        assert abs(np.trace(S0).real - w0) < 1e-10
        post, prob = states.lueders_update(r, alg.central([0.0, 0.0, 1.0]))
        assert abs(prob.item() - w0) < 1e-8
        assert np.abs(full(post) - S0 / w0).max() < 1e-7
        h = 1e-5
        dw = (np.trace(sla.expm(-1j * (t + h) * Heff) @ X0 @ sla.expm(-1j * (t + h) * Heff).conj().T)
              - np.trace(sla.expm(-1j * (t - h) * Heff) @ X0 @ sla.expm(-1j * (t - h) * Heff).conj().T)).real / (2 * h)
        assert abs(dw + np.trace(A @ S0).real) < 1e-7 and dw <= 0
        w0s.append(w0)
    assert w0s[0] >= w0s[1] >= w0s[2] > 0
    if g1 == g2:                                                      # A = gamma N is central: w_0 = e^{-2 gamma t}
        assert abs(w0s[-1] - math.exp(-2 * g1 * 2.0)) < 1e-10


def test_prop_8_12_cesaro_mean_and_trace_norm_bound():
    """Proposition 8.12: the Cesaro means sigma_tau are states with ||T_t sigma_tau - sigma_tau||_1 <= 2t/tau, and
    they converge to a stationary state."""
    rng = seed(5)
    k = 3
    alg = mk([k])
    H, Ls = herm(rng, k), [cplx(rng, k, k) / 2]
    to = lambda X: alg.from_blocks([torch.from_numpy(X)])
    G = dynamics.lindblad_superoperator(alg, to(H), [to(L) for L in Ls])[0][0].numpy()
    rho = rand_state(rng, k).reshape(-1)
    n = k * k
    for tau in (1.0, 5.0, 40.0):
        aug = np.zeros((2 * n, 2 * n), complex)
        aug[:n, :n], aug[:n, n:] = G, np.eye(n)
        integ = sla.expm(tau * aug)[:n, n:]                           # int_0^tau e^{sL} ds
        sigma = (integ @ rho / tau).reshape(k, k)
        assert abs(np.trace(sigma) - 1) < 1e-10 and np.linalg.eigvalsh((sigma + sigma.conj().T) / 2).min() > -1e-12
        for t in (0.3, 2.0):
            diff = (sla.expm(t * G) @ sigma.reshape(-1)).reshape(k, k) - sigma
            assert np.abs(np.linalg.eigvalsh((diff + diff.conj().T) / 2)).sum() <= 2 * t / tau + 1e-10
    assert np.abs(G @ sigma.reshape(-1)).max() < 0.1
    ns = null(G)
    assert ns.shape[1] >= 1


@pytest.mark.parametrize('ks', [[2, 3], [1, 2, 3], [3, 3]])
def test_cor_8_13_one_stationary_state_per_sector(ks):
    """Corollary 8.13 and Theorem D.6 sector by sector: without jumps between sectors each block has a stationary
    state; for random data it is faithful, the commutant is trivial and it is unique."""
    rng = seed(sum(ks) + 2)
    alg = mk(ks)
    H = op_of(alg, [herm(rng, k) for k in ks])
    Ls = [random_op(alg, rng) for _ in range(2)]
    gens = dynamics.lindblad_superoperator(alg, H, Ls)
    for c, k in enumerate(ks):
        Gc = gens[c][0].numpy()
        ns = null(Gc)
        assert ns.shape[1] == 1
        sigma = ns[:, 0].reshape(k, k)
        sigma = sigma / np.trace(sigma)
        assert np.abs(sigma - sigma.conj().T).max() < 1e-10 and np.linalg.eigvalsh(sigma).min() > 1e-6
        ops = [blocks_of(H)[c]] + [blocks_of(L)[c] for L in Ls]
        ops += [A.conj().T for A in ops[1:]]
        assert commutant_space(ops, k).shape[1] == 1
        bl = [np.zeros((kk, kk), complex) for kk in ks]
        bl[c] = sigma
        out = dynamics.lindblad_rhs(H, Ls)(0.0, op_of(alg, bl).matrix)
        assert out.abs().max() < 1e-10                                # a stationary state of the whole algebra


@pytest.mark.parametrize('omega,gamma', [(1.3, 0.4), (0.0, 1.0), (2.0, 3.0)])
def test_prop_8_14_spontaneous_emission(omega, gamma):
    """Proposition 8.14: spectrum {0, -gamma, -gamma/2 +- i omega}, e^{tL} = U_t o A_{1-e^{-gamma t}}, the explicit
    solution, the trace distance formula and its bound e^{-gamma t/2}, with e^{gamma t/2} dist -> |rho_ge(0)|."""
    rng = seed(int(10 * omega + gamma))
    q2 = mk([2])
    H = q2.from_blocks([torch.tensor([[0.0, 0.0], [0.0, omega]], dtype=F64)])
    sm = q2.from_blocks([torch.tensor([[0.0, 1.0], [0.0, 0.0]], dtype=F64)])
    G = dynamics.lindblad_superoperator(q2, H, [sm], rates=[gamma])[0][0].numpy()
    ev = np.linalg.eigvals(G)
    ref = np.array([0, -gamma, -gamma / 2 + 1j * omega, -gamma / 2 - 1j * omega])
    assert all(np.abs(ev - x).min() < 1e-12 for x in ref)
    ground = DensityMatrix(q2, matrix=torch.tensor([[[[1.0, 0.0], [0.0, 0.0]]]], dtype=C128))
    ns = null(G)
    assert ns.shape[1] == 1 and np.abs(ns[:, 0] / ns[0, 0] - np.array([1, 0, 0, 0])).max() < 1e-12
    rho0 = rand_state(rng, 2)
    r0 = op_of(q2, [rho0])
    r0 = DensityMatrix(q2, matrix=r0.matrix)
    for t in (0.5, 2.0, 7.0):
        ch = dynamics.lindblad_channel(q2, H, [sm], t=t, rates=[gamma])
        UA = channels.unitary_channel(dynamics.propagator(H, t)).compose(
            channels.amplitude_damping_channel(q2, 1 - math.exp(-gamma * t)))
        assert (ch.superoperator()[0] - UA.superoperator()[0]).abs().max() < 1e-12
        r = blocks_of(ch(r0))[0]
        assert abs(r[1, 1] - math.exp(-gamma * t) * rho0[1, 1]) < 1e-12
        assert abs(r[0, 1] - np.exp((1j * omega - gamma / 2) * t) * rho0[0, 1]) < 1e-12
        exact = math.sqrt(rho0[1, 1].real ** 2 * math.exp(-2 * gamma * t) + abs(rho0[0, 1]) ** 2 * math.exp(-gamma * t))
        dist = states.trace_distance(ch(r0), ground).item()
        assert abs(dist - exact) < 1e-12 and dist <= math.exp(-gamma * t / 2) + 1e-12
    t = 60 / gamma
    ch = dynamics.lindblad_channel(q2, H, [sm], t=t, rates=[gamma])
    dist = states.trace_distance(ch(r0), ground).item()
    assert abs(math.exp(gamma * t / 2) * dist - abs(rho0[0, 1])) < 1e-6


@pytest.mark.parametrize('ks', [[2, 3], [1, 3]])
def test_exercise_8_1_heisenberg_generator(ks):
    """Chapter 8, Exercise 1: Tr(L(rho) A) = Tr(rho L^*(A)) with the stated L^*, L^*(1) = 0, and L^*(P_c) = 0 without
    jumps between sectors (also when a jump given as an InterSectorChannel has only diagonal blocks)."""
    rng = seed(sum(ks) + 9)
    alg = mk(ks)
    H = op_of(alg, [herm(rng, k) for k in ks])
    Ls = [random_op(alg, rng)]
    diag_jump = channels.random_channel(alg, 2).to_inter_sector()
    Phi = channels.random_inter_sector_channel(alg, alg, kraus_rank=2)
    for between, g in (([diag_jump], [0.8]), ([Phi], [0.8])):
        rhs = dynamics.lindblad_rhs(H, Ls + between, rates=[1.3] + g, algebra=alg)
        Phij = between[0]

        def Lstar(A):
            Hf, Af, Lf = full(H), full(A), full(Ls[0])
            out = 1j * (Hf @ Af - Af @ Hf) + 1.3 * (Lf.conj().T @ Af @ Lf
                                                    - 0.5 * (Lf.conj().T @ Lf @ Af + Af @ Lf.conj().T @ Lf))
            P1 = full(Phij.adjoint()(alg.identity()))
            return out + g[0] * (full(Phij.adjoint()(A)) - 0.5 * (P1 @ Af + Af @ P1))
        for _ in range(2):
            rho, A = random_op(alg, rng), random_op(alg, rng)
            lhs = np.trace(full(alg.operator(rhs(0.0, rho.matrix))) @ full(A))
            assert abs(lhs - np.trace(full(rho) @ Lstar(A))) < 1e-10
        assert np.abs(Lstar(alg.identity())).max() < 1e-12
        for c in range(len(ks)):
            Pc = alg.central([1.0 if d == c else 0.0 for d in range(len(ks))])
            dev = np.abs(Lstar(Pc)).max()
            assert (dev < 1e-12) == (Phij is diag_jump)


@pytest.mark.parametrize('J,g1,g2', [(1.0, 0.6, 0.2), (0.05, 0.6, 0.2), (0.1, 0.6, 0.2)])
def test_exercise_8_2_exceptional_point(J, g1, g2):
    """Chapter 8, Exercise 2: H_eff = [[-i g1/2, J], [J, -i g2/2]] has eigenvalues -i(g1+g2)/4 +- sqrt(J^2 - delta^2),
    delta = (g1-g2)/4, and is not diagonalisable at J = |delta|."""
    He = np.array([[-0.5j * g1, J], [J, -0.5j * g2]])
    delta = (g1 - g2) / 4
    lam = np.sort_complex(np.linalg.eigvals(He))
    root = np.sqrt(J ** 2 - delta ** 2 + 0j)
    ref = np.sort_complex(-0.25j * (g1 + g2) + np.array([root, -root]))
    assert np.abs(lam - ref).max() < 1e-7
    Y = He + 0.25j * (g1 + g2) * np.eye(2)
    assert np.abs(Y @ Y - (J ** 2 - delta ** 2) * np.eye(2)).max() < 1e-15
    if abs(J - abs(delta)) < 1e-15:
        assert np.linalg.matrix_rank(Y) == 1                          # one Jordan block of size 2


@pytest.mark.parametrize('omega,gamma', [(1.3, 0.4), (0.0, 0.7)])
def test_exercise_8_4_pure_dephasing(omega, gamma):
    """Chapter 8, Exercise 4: populations constant, rho_ge(t) = e^{(i omega - 2 gamma)t} rho_ge(0), the explicit
    e^{tL}, spectrum {0,0,-2gamma +- i omega} and a two-dimensional space of stationary states (C = 1)."""
    rng = seed(1)
    q2 = mk([2])
    H = q2.from_blocks([torch.tensor([[0.0, 0.0], [0.0, omega]], dtype=F64)])
    sz = q2.from_blocks([torch.diag(torch.tensor([1.0, -1.0], dtype=F64))])
    G = dynamics.lindblad_superoperator(q2, H, [sz], rates=[gamma])[0][0].numpy()
    ev = np.linalg.eigvals(G)
    for x in (-2 * gamma + 1j * omega, -2 * gamma - 1j * omega):
        assert np.abs(ev - x).min() < 1e-12
    assert int((np.abs(ev) < 1e-12).sum()) == 2
    assert null(G).shape[1] == 2
    rho0 = rand_state(rng, 2)
    for t in (0.4, 1.7):
        Ut = np.diag([1, np.exp(-1j * omega * t)])
        Sz = np.diag([1.0, -1.0])
        ref = ((1 + math.exp(-2 * gamma * t)) / 2 * Ut @ rho0 @ Ut.conj().T
               + (1 - math.exp(-2 * gamma * t)) / 2 * Ut @ Sz @ rho0 @ Sz @ Ut.conj().T)
        out = (sla.expm(t * G) @ rho0.reshape(-1)).reshape(2, 2)
        assert np.abs(out - ref).max() < 1e-12
        assert abs(out[0, 1] - np.exp((1j * omega - 2 * gamma) * t) * rho0[0, 1]) < 1e-12
        assert np.abs(np.diag(out) - np.diag(rho0)).max() < 1e-12


# ================================================================================================================
# Appendix D
# ================================================================================================================
@pytest.mark.parametrize('n', [1, 2, 3])
def test_lemma_D_1_choi_building_blocks(n):
    """Lemma D.1: (1 (x) K) Omega <-> K, <(1(x)K)Omega, (1(x)L)Omega> = Tr K^*L, the Choi matrix of A -> K A L^*,
    J(id) = |Omega><Omega|, and Hermiticity preservation <=> J self-adjoint."""
    rng = seed(n)
    K, L = cplx(rng, n, n), cplx(rng, n, n)
    Omega = np.eye(n).reshape(-1)
    assert np.allclose(np.kron(np.eye(n), K) @ Omega, vec_K(K))
    assert abs(np.vdot(vec_K(K), vec_K(L)) - np.trace(K.conj().T @ L)) < 1e-12
    assert abs(np.vdot(Omega, vec_K(K)) - np.trace(K)) < 1e-12
    assert np.abs(choi(lambda A: K @ A @ L.conj().T, n) - np.outer(vec_K(K), vec_K(L).conj())).max() < 1e-12
    assert np.abs(choi(lambda A: A, n) - np.outer(Omega, Omega)).max() == 0
    assert np.allclose(mat_K(vec_K(K), n), K)
    Xs = [cplx(rng, n * n, n * n) for _ in range(2)]                   # a random superoperator from its Choi matrix
    f = lambda A: sum(A[i, j] * Xs[0].reshape(n, n, n, n)[i, :, j, :] for i in range(n) for j in range(n))
    assert np.abs(choi(f, n) - Xs[0]).max() < 1e-12


def test_lemma_D_2_decomposition():
    """Lemma D.2: a Hermiticity preserving Psi with J positive on W^perp decomposes as sum L A L^* + sum(K_c A Q_c^* +
    Q_c A K_c^*) with Tr Q_c^* L_j = 0 and orthogonal L_j."""
    rng = seed(2)
    n = 3
    Qs = [np.diag([1.0, 0, 0]), np.diag([0, 1.0, 1.0]) / math.sqrt(2)]   # orthonormal q_c = (1 (x) Q_c) Omega
    q = np.stack([vec_K(Q) for Q in Qs], 1)
    P = q @ q.conj().T
    Pp = np.eye(n * n) - P
    Xp = cplx(rng, n * n, 4)
    Y = cplx(rng, n * n, n * n)
    Y = Y + Y.conj().T
    J = Pp @ Xp @ Xp.conj().T @ Pp + (P @ Y + Y @ P - P @ Y @ P)            # positive on W^perp, J not positive
    assert np.linalg.eigvalsh(J).min() < -0.1
    psi = lambda A: sum(A[i, j] * J.reshape(n, n, n, n)[i, :, j, :] for i in range(n) for j in range(n))
    w, U = np.linalg.eigh(Pp @ J @ Pp)
    Lj = [math.sqrt(x) * mat_K(U[:, i], n) for i, x in enumerate(w) if x > 1e-10]
    assert len(Lj) == np.linalg.matrix_rank(Pp @ J @ Pp, tol=1e-10)
    xis = [J @ q[:, c] - 0.5 * sum((q[:, d].conj() @ J @ q[:, c]) * q[:, d] for d in range(2)) for c in range(2)]
    Kc = [mat_K(x, n) for x in xis]
    for _ in range(3):
        A = cplx(rng, n, n)
        rec = sum(L @ A @ L.conj().T for L in Lj) + sum(K @ A @ Q.conj().T + Q @ A @ K.conj().T for K, Q in zip(Kc, Qs))
        assert np.abs(rec - psi(A)).max() < 1e-10
    for L in Lj:
        assert all(abs(np.trace(Q.conj().T @ L)) < 1e-10 for Q in Qs)
    for i, Li in enumerate(Lj):
        for j, Lk in enumerate(Lj):
            if i != j:
                assert abs(np.trace(Li.conj().T @ Lk)) < 1e-10


@pytest.mark.parametrize('n', [2, 3])
def test_thm_D_3_conditional_complete_positivity(n):
    """Theorem D.3: GKSL generators are CCP with CP exponentials; (3) => (1) also without trace preservation; the
    transpose minus id is Hermiticity preserving, not CCP, and e^{tL} is positive, TP but not CP."""
    rng = seed(n + 30)
    H, Ls = herm(rng, n), [cplx(rng, n, n) for _ in range(2)]
    S = superop(lambda R: gksl_np(H, Ls, R), n)
    J = choi_from_superop(S, n)
    q = np.eye(n).reshape(-1) / math.sqrt(n)
    Pp = np.eye(n * n) - np.outer(q, q)
    assert np.abs(J - J.conj().T).max() < 1e-12
    assert np.linalg.eigvalsh(Pp @ J @ Pp).min() > -1e-10 and np.linalg.eigvalsh(J).min() < -0.1
    Gm = cplx(rng, n, n)                                              # L = Phi + G . + . G^*, not TP
    S3 = superop(lambda R: sum(L @ R @ L.conj().T for L in Ls) + Gm @ R + R @ Gm.conj().T, n)
    for t in (0.05, 0.5, 2.0):
        assert np.linalg.eigvalsh(choi_from_superop(sla.expm(t * S), n)).min() > -1e-10
        assert np.linalg.eigvalsh(choi_from_superop(sla.expm(t * S3), n)).min() > -1e-9 * np.abs(sla.expm(t * S3)).max()
    Lth = superop(lambda X: X.T, n) - np.eye(n * n)
    v = vec_K(unit(1, 0, n) - unit(0, 1, n))
    assert abs(v @ q) < 1e-15 and abs(v @ choi_from_superop(Lth, n) @ v + 2) < 1e-12
    for t in (0.1, 1.0, 3.0):
        Tt = sla.expm(t * Lth)
        b = (1 - math.exp(-2 * t)) / 2
        assert np.abs(Tt - ((1 - b) * np.eye(n * n) + b * superop(lambda X: X.T, n))).max() < 1e-12
        assert abs(v @ choi_from_superop(Tt, n) @ v + (1 - math.exp(-2 * t))) < 1e-12
        for _ in range(3):                                            # positive on states
            R = rand_state(rng, n)
            assert np.linalg.eigvalsh((Tt @ R.reshape(-1)).reshape(n, n)).min() > -1e-12


def gksl_from_generator(S, n, tol=1e-10):
    J = choi_from_superop(S, n)
    q = np.eye(n).reshape(-1) / math.sqrt(n)
    Pp = np.eye(n * n) - np.outer(q, q)
    w, U = np.linalg.eigh(Pp @ J @ Pp)
    jumps = [math.sqrt(x) * mat_K(U[:, i], n) for i, x in enumerate(w) if x > tol]
    xi = J @ q - 0.5 * (q.conj() @ J @ q) * q
    G = mat_K(xi, n) / math.sqrt(n)
    return 0.5j * (G - G.conj().T), jumps, G


@pytest.mark.parametrize('n,nj', [(2, 7), (3, 2), (4, 1)])
def test_thm_D_4_gksl_recovered_from_generator(n, nj):
    """Theorem D.4: H (traceless) and r <= n^2 - 1 traceless orthogonal jumps recovered from the Choi matrix of a
    generator reproduce it, with Tr G real and G + G^* = -sum L^*L."""
    rng = seed(n * 7 + nj)
    alg = mk([n])
    H, Ls = herm(rng, n), [cplx(rng, n, n) / 2 for _ in range(nj)]
    to = lambda X: alg.from_blocks([torch.from_numpy(X)])
    S = dynamics.lindblad_superoperator(alg, to(H), [to(L) for L in Ls])[0][0].numpy()
    Hr, Lr, G = gksl_from_generator(S, n)
    assert len(Lr) <= n * n - 1 and len(Lr) == min(nj, n * n - 1)
    assert abs(np.trace(Hr)) < 1e-12 and all(abs(np.trace(L)) < 1e-12 for L in Lr)
    assert abs(np.trace(G).imag) < 1e-12
    assert np.abs(G + G.conj().T + sum(L.conj().T @ L for L in Lr)).max() < 1e-10
    for i in range(len(Lr)):
        for j in range(i):
            assert abs(np.trace(Lr[i].conj().T @ Lr[j])) < 1e-10
    S_rec = dynamics.lindblad_superoperator(alg, to(Hr), [to(L) for L in Lr])[0][0].numpy()
    assert np.abs(S_rec - S).max() < 1e-10


def test_thm_D_5_freedom_in_gksl_form():
    """Theorem D.5: (i) H + c1, (ii) unitary mixing of the jumps, (iii) L + a1 with the compensating H all leave L
    unchanged; normalised data with the same generator have the same H and unitarily related jumps."""
    rng = seed(17)
    n = 3
    alg = mk([n])
    to = lambda X: alg.from_blocks([torch.from_numpy(np.asarray(X, complex))])
    gen = lambda H, Ls: dynamics.lindblad_superoperator(alg, to(H), [to(L) for L in Ls])[0][0].numpy()
    H, Ls = herm(rng, n), [cplx(rng, n, n) for _ in range(2)]
    S = gen(H, Ls)
    assert np.abs(gen(H + 0.77 * np.eye(n), Ls) - S).max() < 1e-12                                    # (i)
    u = sla.qr(cplx(rng, 3, 3))[0]
    Lpad = Ls + [np.zeros((n, n))]
    assert np.abs(gen(H, [sum(u[j, l] * Lpad[l] for l in range(3)) for j in range(3)]) - S).max() < 1e-12   # (ii)
    a = cplx(rng, 2)
    H3 = H - 0.5j * sum(np.conj(aj) * L - aj * L.conj().T for aj, L in zip(a, Ls))
    assert np.abs(H3 - H3.conj().T).max() < 1e-12
    assert np.abs(gen(H3, [L + aj * np.eye(n) for aj, L in zip(a, Ls)]) - S).max() < 1e-12             # (iii)
    a0 = [-np.trace(L) / n for L in Ls]                                  # (c): normal form by (iii) then (i)
    Ln = [L + aj * np.eye(n) for L, aj in zip(Ls, a0)]
    Hn = H - 0.5j * sum(np.conj(aj) * L - aj * L.conj().T for aj, L in zip(a0, Ls))
    Hn = Hn - np.trace(Hn) / n * np.eye(n)
    Hr, Lr, _ = gksl_from_generator(S, n)
    assert np.abs(Hr - Hn).max() < 1e-10                                # (b): H unique
    Wr, Wn = np.stack([vec_K(K) for K in Lr], 1), np.stack([vec_K(K) for K in Ln], 1)
    U = np.linalg.lstsq(Wr, Wn, rcond=None)[0]
    assert np.abs(Wr @ U - Wn).max() < 1e-10 and np.abs(U.conj().T @ U - np.eye(2)).max() < 1e-10
    assert np.abs(gen(H + 0.1 * herm(rng, n), Ls) - S).max() > 1e-3       # a different H changes the generator


@pytest.mark.parametrize('example', ['random', 'dephasing', 'absorbing'])
def test_thm_D_6_uniqueness_of_stationary_state(example):
    """Theorem D.6: faithful stationary state and trivial commutant give ker L = C sigma and ker L^* = C 1; pure
    dephasing (commutant too big) and L_j = |e_j><e_3| on M_3 (no faithful stationary state) show both hypotheses
    are needed; the dissipation identity L^*(X^*X) - X^*L^*(X) - L^*(X^*)X = sum [L,X]^*[L,X]."""
    rng = seed(23)
    if example == 'random':
        n = 3
        H, Ls = herm(rng, n), [cplx(rng, n, n) / 2 for _ in range(2)]
    elif example == 'dephasing':
        n = 2
        H, Ls = np.diag([0.0, 1.3]), [np.diag([1.0, -1.0])]
    else:
        n = 3
        H, Ls = np.zeros((3, 3)), [unit(0, 2, 3), unit(1, 2, 3)]
    S = superop(lambda R: gksl_np(H, Ls, R), n)
    Lstar = lambda A: 1j * (H @ A - A @ H) + sum(L.conj().T @ A @ L - 0.5 * (L.conj().T @ L @ A + A @ L.conj().T @ L)
                                               for L in Ls)
    for _ in range(2):
        R, A = cplx(rng, n, n), cplx(rng, n, n)
        assert abs(np.trace((S @ R.reshape(-1)).reshape(n, n) @ A) - np.trace(R @ Lstar(A))) < 1e-10
        X = cplx(rng, n, n)
        lhs = Lstar(X.conj().T @ X) - X.conj().T @ Lstar(X) - Lstar(X.conj().T) @ X
        rhs = sum((L @ X - X @ L).conj().T @ (L @ X - X @ L) for L in Ls)
        assert np.abs(lhs - rhs).max() < 1e-10
    ker = null(S).shape[1]
    kerT = null(S.T).shape[1]
    comm = commutant_space([H] + Ls + [L.conj().T for L in Ls], n).shape[1]
    assert ker == kerT
    assert np.abs(S.T @ np.eye(n).reshape(-1)).max() < 1e-12         # L^*(1) = 0
    if example == 'random':
        sigma = null(S)[:, 0].reshape(n, n)
        sigma = sigma / np.trace(sigma)
        assert np.linalg.eigvalsh((sigma + sigma.conj().T) / 2).min() > 1e-6 and comm == 1 and ker == 1
    elif example == 'dephasing':
        assert comm == 2 and ker == 2                                  # 1/2 is faithful and stationary
        assert np.abs(S @ (np.eye(2) / 2).reshape(-1)).max() < 1e-15
    else:
        assert comm == 1 and ker == 4
        B = null(S)
        for i in range(B.shape[1]):                                    # every stationary matrix lives on span(e1,e2)
            X = B[:, i].reshape(3, 3)
            assert np.abs(X[2, :]).max() < 1e-10 and np.abs(X[:, 2]).max() < 1e-10


@pytest.mark.parametrize('ks', [[1, 2], [2, 2], [1, 1, 2]])
def test_thm_D_7_and_cor_D_8_generators_on_sectors(ks):
    """Theorem D.7 and Corollary D.8: for lindblad_rhs on M with jumps inside and between sectors, J(L o E) is positive
    off the sector vectors, and the five steps give L = Phi - {Phi^*(1), .}/2 - i[H, .] with Phi CP between sectors,
    G + G^* = -Phi^*(1), and the Heisenberg generator Phi^* + G^* . + . G."""
    rng = seed(sum(ks) * 11 + len(ks))
    alg = mk(ks)
    D = sum(ks)
    o = np.cumsum([0] + ks)
    Pc = []
    for c, k in enumerate(ks):
        P = np.zeros((D, D))
        P[o[c]:o[c] + k, o[c]:o[c] + k] = np.eye(k)
        Pc.append(P)
    E = lambda X: sum(P @ X @ P for P in Pc)
    H = op_of(alg, [herm(rng, k) for k in ks])
    inside = random_op(alg, rng)
    between = channels.random_inter_sector_channel(alg, alg, kraus_rank=2)
    rhs = dynamics.lindblad_rhs(H, [inside, between], rates=[0.7, 1.3], algebra=alg)
    LM = lambda X: full(alg.operator(rhs(0.0, unflat(alg, X).matrix)))
    J = choi_from_superop(superop(lambda X: LM(E(X)), D), D)
    q = np.stack([vec_K(P) / math.sqrt(k) for P, k in zip(Pc, ks)], 1)
    Pw = np.eye(D * D) - q @ q.T
    assert np.linalg.eigvalsh(Pw @ J @ Pw).min() > -1e-10                       # Step 2
    w, U = np.linalg.eigh(Pw @ J @ Pw)
    jumps = [math.sqrt(x) * mat_K(U[:, i], D) for i, x in enumerate(w) if x > 1e-10]
    xis = [J @ q[:, c] - 0.5 * sum((q[:, d] @ J @ q[:, c]) * q[:, d] for d in range(len(ks))) for c in range(len(ks))]
    G = E(sum(mat_K(xi, D) @ P / math.sqrt(k) for xi, P, k in zip(xis, Pc, ks)))
    Phi = lambda X: E(sum(K @ X @ K.conj().T for K in jumps))
    Phi1 = E(sum(K.conj().T @ K for K in jumps))
    Hm = 0.5j * (G - G.conj().T)
    assert np.abs(G + G.conj().T + Phi1).max() < 1e-10                         # Step 5
    assert np.abs(Hm - E(Hm)).max() < 1e-14
    for _ in range(3):
        X = E(cplx(rng, D, D))
        gk = Phi(X) - 0.5 * (Phi1 @ X + X @ Phi1) - 1j * (Hm @ X - X @ Hm)
        assert np.abs(LM(X) - gk).max() < 1e-10
        Y = E(cplx(rng, D, D))
        dual = E(sum(K.conj().T @ Y @ K for K in jumps)) + G.conj().T @ Y + Y @ G
        assert abs(np.trace(LM(X) @ Y) - np.trace(X @ dual)) < 1e-10           # Corollary D.8
    one = np.eye(D)
    assert np.abs(E(sum(K.conj().T @ one @ K for K in jumps)) + G.conj().T + G).max() < 1e-10   # Psi(1)+K+K^*=0
    Pom = np.eye(D * D) - np.outer(np.eye(D).reshape(-1), np.eye(D).reshape(-1)) / D
    if len(ks) > 1:
        assert np.linalg.eigvalsh(Pom @ J @ Pom).min() < -1e-3                  # not CCP on Omega^perp


def hybrid_code(rng, dims, D=8, r=3, s=2):
    """A code C (+)_c C^{k_c} (x) C^{m_c} in C^D and noise acting on the gauge factors, mixed by a unitary."""
    dimC = sum(k * m for k, m in dims)
    U8 = sla.qr(cplx(rng, D, D))[0]
    Vc = sla.qr(cplx(rng, D, dimC), mode='economic')[0]
    W, o = [], 0
    for k, m in dims:
        W.append(Vc[:, o:o + k * m])
        o += k * m
    Pcode = sum(Wc @ Wc.conj().T for Wc in W)
    Y, o = [], 0
    for k, m in dims:
        Y.append(U8[:, o:o + k * s])
        o += k * s
    B = [sla.qr(cplx(rng, s * r, m), mode='economic')[0].reshape(r, s, m) for _, m in dims]
    Ecode = [sum(Yc @ np.kron(np.eye(k), Bc[i]) @ Wc.conj().T for Yc, Bc, Wc, (k, m) in zip(Y, B, W, dims))
             for i in range(r)]
    kraus = Ecode + [sla.qr(cplx(rng, D, D))[0] @ (np.eye(D) - Pcode)]
    mix = sla.qr(cplx(rng, r + 1, r + 1))[0]
    kraus = [sum(mix[i, l] * K for l, K in enumerate(kraus)) for i in range(r + 1)]
    return W, Pcode, kraus


@pytest.mark.parametrize('dims', [[(2, 2), (1, 1)], [(1, 2), (1, 1), (1, 1)], [(2, 1)], [(1, 3)]])
def test_thm_D_9_operator_algebra_error_correction(dims):
    """Theorem D.9: P E_i^* E_j P lies in A'; the recovery built in the proof (Steps 1-3) is a channel with
    P (R o Phi)^*(A) P = A for A in A, and R_m E_i P commutes with A (necessity)."""
    rng = seed(len(dims) * 100 + sum(k * m for k, m in dims))
    D = 8
    W, Pcode, kraus = hybrid_code(rng, dims, D=D)
    assert np.abs(sum(K.conj().T @ K for K in kraus) - np.eye(D)).max() < 1e-12
    units = []
    for Wc, (k, m) in zip(W, dims):
        for a1 in range(k):
            for a2 in range(k):
                units.append(Wc @ np.kron(unit(a1, a2, k), np.eye(m)) @ Wc.conj().T)
    comm = lambda X: max(np.abs(X @ A - A @ X).max() for A in units)
    assert max(comm(Pcode @ Ki.conj().T @ Kj @ Pcode) for Ki in kraus for Kj in kraus) < 1e-12
    R_kraus, Qsum = [], np.zeros((D, D), complex)
    for Wc, (k, m) in zip(W, dims):
        F = [Ki @ Wc[:, b::m] for Ki in kraus for b in range(m)]
        Lam = np.array([[(Fa.conj().T @ Fb)[0, 0] for Fb in F] for Fa in F])
        assert max(np.abs(Fa.conj().T @ Fb - Lam[x, y] * np.eye(k)).max()
                   for x, Fa in enumerate(F) for y, Fb in enumerate(F)) < 1e-12               # Step 1
        assert np.abs(Lam - Lam.conj().T).max() < 1e-12 and np.linalg.eigvalsh(Lam).min() > -1e-12   # Step 2
        d, u = np.linalg.eigh(Lam)
        iota = Wc[:, 0::m]
        for kk in range(len(d)):
            if d[kk] > 1e-12:
                Vk = sum(u[x, kk] * Fx for x, Fx in enumerate(F)) / math.sqrt(d[kk])
                assert np.abs(Vk.conj().T @ Vk - np.eye(k)).max() < 1e-10
                R_kraus.append(iota @ Vk.conj().T)
                Qsum += Vk @ Vk.conj().T
    assert np.abs(Qsum @ Qsum - Qsum).max() < 1e-10                                        # Q is a projection
    R_kraus.append(np.eye(D) - Qsum)
    assert np.abs(sum(R.conj().T @ R for R in R_kraus) - np.eye(D)).max() < 1e-10
    dual = lambda X: sum(Ki.conj().T @ R.conj().T @ X @ R @ Ki for R in R_kraus for Ki in kraus)
    for _ in range(3):
        A = sum(complex(z) * Ue for z, Ue in zip(cplx(rng, len(units)), units))
        assert np.abs(Pcode @ dual(A) @ Pcode - A).max() < 1e-10
    assert max(comm(R @ Ki @ Pcode) for R in R_kraus for Ki in kraus) < 1e-10
    dimC = sum(k * m for k, m in dims)
    KL = max(np.abs(Pcode @ Ki.conj().T @ Kj @ Pcode - np.trace(Pcode @ Ki.conj().T @ Kj) / dimC * Pcode).max()
             for Ki in kraus for Kj in kraus)
    assert (KL < 1e-10) == (len(dims) == 1 and dims[0][1] == 1)       # Knill-Laflamme holds only for A = B(C)
    if not (len(dims) == 1 and dims[0][1] == 1):
        Bg = Pcode @ cplx(rng, D, D) @ Pcode
        assert np.abs(Pcode @ dual(Bg) @ Pcode - Bg).max() > 0.05


def test_cor_D_10_bit_flip_code():
    """Corollary D.10(a): for A = B(C) the condition is Knill-Laflamme; the three-qubit code against at most one bit
    flip, Lambda = diag(p), and recovery of every operator on the code."""
    X = np.array([[0, 1], [1, 0]])
    I2 = np.eye(2)
    Xs = [np.eye(8)] + [np.kron(np.kron(*(X if j == i else I2 for j in range(2))), X if i == 2 else I2) for i in range(3)]
    p = np.array([0.7, 0.1, 0.1, 0.1])
    E = [math.sqrt(pi) * Xi for pi, Xi in zip(p, Xs)]
    code = np.zeros((8, 2))
    code[0, 0] = code[7, 1] = 1
    P = code @ code.T
    Lam = np.array([[np.trace(P @ Ei.T @ Ej @ P) / 2 for Ej in E] for Ei in E])
    assert np.abs(Lam - np.diag(p)).max() < 1e-12
    for Ei in E:
        for Ej in E:
            M = P @ Ei.T @ Ej @ P
            assert np.abs(M - np.trace(M) / 2 * P).max() < 1e-12
    R = [P @ Xi for Xi in Xs]
    assert np.abs(sum(Rm.T @ Rm for Rm in R) - np.eye(8)).max() < 1e-12
    rng = seed(0)
    B = P @ cplx(rng, 8, 8) @ P
    dual = sum(Ei.T @ Rm.T @ B @ Rm @ Ei for Rm in R for Ei in E)
    assert np.abs(P @ dual @ P - B).max() < 1e-12
