"""Regression tests for correctness issues found by the independent audit (dense references)."""
import numpy as np
import pytest
import scipy.linalg
import torch

from torch_vn_algebra import TypeIAlgebra, states as st, channels as ch, dynamics as dy


def _real(ks=(2, 3)):
    return TypeIAlgebra(list(ks), list(ks), complex_valued=False, precision='double')


def _cplx(ks=(2, 0, 3)):
    return TypeIAlgebra([max(k, 1) for k in ks], list(ks), precision='double')


def _dense(alg, m):
    """Block-diagonal numpy matrix of one sample."""
    return scipy.linalg.block_diag(*[m[c, :k, :k].numpy() for c, k in enumerate(alg.k_factors) if k])


def test_lindblad_real_algebra_with_jumps():
    """Real jumps on a real algebra: the generator is real; it used to crash on a complex H_eff."""
    torch.manual_seed(0)
    alg = _real()
    rho = st.random_density_matrix(alg, 1)
    L = alg.operator(torch.randn(1, alg.C, alg.k_max, alg.k_max, dtype=torch.float64))
    ys = dy.lindblad_evolve(rho, None, [L], [0.0, 0.7], substeps=200)
    # dense reference: exp(t L) on vec(rho), row-major
    l, r0 = _dense(alg, L.matrix[0]), _dense(alg, rho.matrix[0])
    n = l.shape[0]
    I = np.eye(n)
    gen = np.kron(l, l) - 0.5 * np.kron(l.T @ l, I) - 0.5 * np.kron(I, (l.T @ l).T)
    ref = (scipy.linalg.expm(0.7 * gen) @ r0.reshape(-1)).reshape(n, n)
    assert ys[-1].matrix.dtype == torch.float64
    assert np.allclose(_dense(alg, ys[-1].matrix[0]), ref, atol=1e-8)


def test_dynamics_ignore_padding_garbage():
    torch.manual_seed(1)
    alg = TypeIAlgebra([2, 3], [2, 3], precision='double')
    rho = st.random_density_matrix(alg, 1)
    H = torch.randn(1, alg.C, 3, 3, dtype=alg.dtype); H = H + H.conj().transpose(-2, -1)
    L = torch.randn(1, alg.C, 3, 3, dtype=alg.dtype)
    Hg, Lg = H.clone(), L.clone()
    Hg[:, 0, 2, :] = 5.0; Hg[:, 0, :, 2] = 5.0; Lg[:, 0, 2, 0] = 3.0
    a = dy.lindblad_evolve(rho, alg.operator(H), [alg.operator(L)], [0.0, 0.5])[-1].matrix
    b = dy.lindblad_evolve(rho, Hg, [Lg], [0.0, 0.5])[-1].matrix
    assert torch.allclose(a, b, atol=1e-12)
    psi = torch.randn(1, alg.C, 3, dtype=alg.dtype); psi[:, 0, 2] = 0
    psig = psi.clone(); psig[:, 0, 2] = 0.7
    x = dy.schrodinger_rk4(psi, alg.operator(H), [0.0, 0.5])
    y = dy.schrodinger_rk4(psig, Hg, [0.0, 0.5], algebra=alg)
    assert torch.allclose(x, y, atol=1e-12)


def test_density_matrix_masks_padding():
    alg = _cplx()
    rho = st.random_density_matrix(alg, 2)
    m = 0.5 * rho.matrix.clone()
    m[:, 1, 0, 0] = 0.5                                   # mass in the empty sector: not part of the algebra
    with pytest.raises(ValueError, match="unit trace"):
        st.DensityMatrix(alg, matrix=m)
    d = rho.matrix.clone(); d[:, 0, 2, 2] = 3.0           # garbage in the padding of a 2x2 block
    s = st.DensityMatrix.from_density(alg, d, normalize=True)
    assert torch.allclose(s.matrix, rho.matrix, atol=1e-12)


def test_bures_real_rejected():
    with pytest.raises(ValueError, match="complex"):
        st.random_density_matrix(_real(), 4, measure='bures')


def test_gibbs_negative_beta():
    alg = TypeIAlgebra([2, 3], [2, 3], precision='double')
    H = alg.from_blocks([torch.diag(torch.tensor([0.0, 1000.0])), torch.diag(torch.tensor([0.0, 1.0, 2.0]))])
    H = alg.operator(H.matrix, is_self_adjoint=True)
    g = st.gibbs_state(H, -1.0).matrix
    assert torch.isfinite(g).all()
    assert abs(g[0, 0, 1, 1].real.item() - 1.0) < 1e-12
    g2 = st.gibbs_state(H, 0.3).matrix
    e = np.exp(-0.3 * np.array([0, 1000, 0, 1, 2.0])); e /= e.sum()
    assert np.allclose(np.diag(_dense(alg, g2[0])).real, e)


def test_mixed_unitary_weights_validated():
    alg = _cplx((2, 2))
    with pytest.raises(ValueError):
        ch.random_mixed_unitary_channel(alg, 2, weights=torch.tensor([-0.5, 1.5]))
    with pytest.raises(ValueError):
        ch.random_mixed_unitary_channel(alg, 2, batch_size=1, weights=torch.ones(2, 2))


def test_from_superoperator_double_tolerance():
    """A non-CP map at the 5e-5 level is rejected in double precision (it used to be clamped silently)."""
    alg = TypeIAlgebra([2], [2], precision='double')
    k = 2
    I = np.eye(k * k)
    T = np.zeros((k * k, k * k))                         # transpose map on row-major vec
    for i in range(k):
        for j in range(k):
            T[i * k + j, j * k + i] = 1
    eps = 5e-5
    S = torch.tensor((1 - eps) * I + eps * T, dtype=torch.complex128)
    with pytest.raises(ValueError, match="completely positive"):
        ch.Channel.from_superoperator(alg, [S])
    ch.Channel.from_superoperator(alg, [S], tol=1e-3)     # explicit tolerance still accepts


def test_to_inter_sector_keeps_tp_flag():
    alg = _cplx((2, 3))
    Phi = ch.random_channel(alg, 2)
    assert Phi.to_inter_sector()._trace_preserving is True
