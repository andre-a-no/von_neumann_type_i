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


def test_lanczos_large_shift_not_truncated():
    """A large uniform field is a constant shift in a sector; it must not make Lanczos stop early."""
    from torch_vn_algebra.chains import SpinChain
    from torch_vn_algebra import krylov
    ch_ = SpinChain(10, 'open', complex_valued=False, precision='double')
    H = ch_.xxz_sparse(1.0, 1.0, 1e5, sector=2)
    Hd = H.to_dense()[0].numpy()
    w, V = np.linalg.eigh(Hd)
    rng = np.random.default_rng(0)
    v = V[:, 5] + 1e-6 * rng.standard_normal(V.shape[0]); v /= np.linalg.norm(v)
    psi0 = torch.tensor(v)[None]
    out = krylov.evolve(H, psi0, [0.0, 0.2])[-1, 0].numpy()
    ref = V @ (np.exp(-0.2j * w) * (V.T @ v))
    assert np.abs(out - ref).max() < 1e-9
    E, _, _ = krylov.ground_state(H, v0=torch.tensor(V[:, 5] + 1e-6 * V[:, 0])[None])
    assert abs(E.item() - w[0]) < 1e-6


def test_krylov_zero_vector():
    from torch_vn_algebra.chains import SpinChain
    from torch_vn_algebra import krylov
    H = SpinChain(8, 'open', precision='double').xxz_sparse(1.0, 1.0, None, sector=4)
    psi0 = torch.randn(2, H.dim, dtype=torch.complex128); psi0[1] = 0
    out = krylov.evolve(H, psi0, [0.0, 1.0])
    assert torch.isfinite(out).all() and out[-1, 1].abs().max() == 0


def test_hilbert_norm_flat_vectors():
    from torch_vn_algebra.hilbert_space import HilbertSpace
    hs = HilbertSpace(3, 3, batch_size=2, n_channels=2, precision='double', device='cpu')
    v = torch.randn(2, 2, 3, dtype=torch.complex128)
    n = hs.norm(v).reshape(-1)
    assert torch.allclose(n ** 2, hs.inner_product(v, v).reshape(-1).real)
    assert torch.allclose(hs.norm(hs.normalize(v)).reshape(-1), torch.ones(2, dtype=torch.float64))


def test_krylov_times_in_double():
    """Times given as Python floats must not be rounded to float32 (phase error E * 3e-9 at t = 0.2)."""
    from torch_vn_algebra.chains import SpinChain
    from torch_vn_algebra import krylov
    H = SpinChain(8, 'open', precision='double').xxz_sparse(1.0, 1.0, 1e3, sector=4)
    Hd = H.to_dense()[0].numpy()
    w, V = np.linalg.eigh(Hd)
    v = V[:, 0] + V[:, 3]; v /= np.linalg.norm(v)
    out = krylov.evolve(H, torch.tensor(v, dtype=torch.complex128)[None], [0.0, 0.2])[-1, 0].numpy()
    assert np.abs(out - V @ (np.exp(-0.2j * w) * (V.conj().T @ v))).max() < 1e-11


# ----- algebra / functional calculus -----
def test_inverse_threshold_per_block():
    A = TypeIAlgebra([2, 2], [2, 2], precision='double')
    H = A.from_blocks([torch.diag(torch.tensor([1e10, 1.0], dtype=A.dtype)),
                       torch.diag(torch.tensor([1e-7, 2e-7], dtype=A.dtype))])
    H = A.operator(H.matrix, is_self_adjoint=True)
    d = H.inverse().matrix[0, 1].diagonal().real
    assert torch.allclose(d, torch.tensor([1e7, 5e6], dtype=torch.float64))


def test_power_negative_support_relative():
    torch.manual_seed(0)
    A = TypeIAlgebra([4], [4])                         # single precision
    v = torch.randn(4, 2, dtype=torch.complex64)
    P = A.operator((v @ v.mH)[None, None], is_self_adjoint=True, is_positive=True)
    ref = torch.linalg.pinv(P.matrix[0, 0].to(torch.complex128), hermitian=True, rtol=1e-5)   # support of the float32 matrix
    assert (P.power(-1).matrix[0, 0].to(torch.complex128) - ref).abs().max() < 1e-3 * ref.abs().max()


def test_small_operators_classified_by_relative_tolerance():
    A = TypeIAlgebra([2], [2])
    N = A.from_blocks([1e-6 * torch.tensor([[0.0, 1.0], [0.0, 0.0]])])
    assert not N.is_self_adjoint and not N.is_normal
    assert not A.operator(-1e-6 * A.identity().matrix).is_positive
    assert A.operator(1e-5 * A.identity().matrix).is_invertible


def test_zero_times_operator_flags():
    A = TypeIAlgebra([2], [2], precision='double')
    X = A.from_blocks([torch.tensor([[0.0, 1.0], [0.0, 0.0]], dtype=A.dtype)])
    assert not X.is_self_adjoint
    assert torch.all((0 * X).eigenvalues()[0] == 0)


def test_michelson_scale_invariant():
    A = TypeIAlgebra([2], [2], precision='double')
    P = A.operator(torch.diag(torch.tensor([1.0, 0.25], dtype=A.dtype))[None, None], is_self_adjoint=True,
                   is_positive=True)
    for s in (1.0, 1e-150):
        assert abs((P * s).michelson_contrast.item() - 0.6) < 1e-12


def test_frobenius_no_overflow_single():
    A = TypeIAlgebra([2], [2])
    X = A.operator(1e20 * A.identity().matrix)
    assert abs(X.frobenius_norm().item() / (1e20 * 2 ** 0.5) - 1) < 1e-6


def test_normal_apply_function_no_spurious_grad():
    A = TypeIAlgebra([3], [3], precision='double')
    U = A.random_unitary_operator(2)
    F = U.apply_function(torch.exp).matrix
    assert not F.requires_grad
    t = torch.tensor([0.3, 0.7], dtype=torch.float64, requires_grad=True)
    G = U.apply_function(lambda w: torch.exp(1j * t[:, None] * w)).matrix
    assert G.requires_grad


def test_holomorphic_scalar_block_and_zero_result():
    A = TypeIAlgebra([3], [3], precision='double')
    I = A.operator_from_eigenvalues(lambda k: torch.ones(k, dtype=torch.float64))
    assert torch.allclose(I.holomorphic_function(torch.exp).matrix,
                          np.e * A.identity().matrix, atol=1e-12)
    Z = A.operator(torch.zeros(1, 1, 3, 3, dtype=A.dtype))
    assert Z.holomorphic_function(torch.sin, atol=1.0).matrix.abs().max() < 1e-12


def test_operator_validate_option():
    A = TypeIAlgebra([2], [2], precision='double')
    m = torch.tensor([[[[1.0, 2.0], [0.0, -3.0]]]], dtype=A.dtype)
    with pytest.raises(ValueError):
        A.operator(m, validate=True, is_self_adjoint=True)


def test_central_python_numbers_in_double():
    A = TypeIAlgebra([2, 3], [2, 3], precision='double')
    assert abs(A.central([1 / 3, 2.0]).matrix[0, 0, 0, 0].real.item() - 1 / 3) < 1e-16
