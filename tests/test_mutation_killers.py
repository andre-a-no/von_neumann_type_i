"""
Mutation-killing tests for torch_vn_algebra (from a mutation-testing run over torch_vn_algebra/*.py).

Each test names the surviving mutant(s) it kills in a `# kills:` comment (ids <module>#<n> of the mutation run).
Every test passes on the unmutated library and was checked to fail on its mutant, one mutant at a time. Public API plus numpy/scipy/torch, except three white-box checks that use a private name:
cost._fmt_seconds (the ETA text is only reachable through timing), functional._ALPHAS (spectra built to tie for
every mixing coefficient) and the `.kraus` / `.signs` / `.G` attributes of channels and parameterisations.
Deterministic (fixed seeds), CPU, ~10 s in total.
"""
import math

import numpy as np
import pytest
import torch

import torch_vn_algebra as tv
from torch_vn_algebra import TypeIAlgebra, divergences as dv


def _diag_op(alg, *diags):
    """Operator of `alg` with diagonal blocks (float64 / complex128 lists)."""
    return alg.from_blocks([torch.diag(torch.tensor(d, dtype=torch.float64)) for d in diags])


# kills: divergences#3683
def test_umegaki_rtol_zero_exact_kernel_of_sigma():
    """D(rho||sigma) with rtol=0 and an exactly singular sigma: log is taken on lambda > 0 only, D = log 2."""
    alg = TypeIAlgebra([3], [3], precision='double')
    rho = _diag_op(alg, [1.0, 0.0, 0.0])
    sigma = _diag_op(alg, [0.5, 0.5, 0.0])
    D = dv.umegaki_relative_entropy(rho, sigma, rtol=0.0)
    assert torch.isfinite(D).all()
    assert abs(float(D) - math.log(2)) < 1e-12


def _herm(alg, seed=0, batch=2):
    """Random self-adjoint operator of `alg` (batch samples), from explicit blocks."""
    g = torch.Generator().manual_seed(seed)
    blocks = []
    for k in alg.k_factors:
        x = torch.randn(batch, k, k, 2, generator=g, dtype=torch.float64)
        x = torch.complex(x[..., 0], x[..., 1]) if alg.complex_valued else x[..., 0]
        blocks.append(x + x.mH)
    return alg.from_blocks(blocks)


# kills: algebra#331
def test_shape_of_lazy_operator_does_not_materialize_it():
    """Operator.shape of a lazy operator with known batch does not evaluate it (here: under no_grad, which would
    otherwise drop the gradient of alg.central with tensor values)."""
    alg = TypeIAlgebra([2, 3], [2, 3], precision='double')
    s = torch.tensor(2.0, dtype=torch.float64, requires_grad=True)
    op = alg.central([s, 1.0])
    with torch.no_grad():
        assert tuple(op.shape) == (1, 2, 3, 3)
    op.trace.real.sum().backward()
    assert s.grad is not None and abs(float(s.grad) - 2.0) < 1e-12


# kills: algebra#734
def test_cached_eigh_and_eigenvalues_match_uncached():
    """With the spectral cache on, eigh() and eigenvalues() equal the uncached results (also for an empty sector)."""
    alg = TypeIAlgebra([3, 1, 2, 3], [3, 0, 2, 3], precision='double')
    A = _herm(alg, seed=1, batch=2)
    ref_w = A.eigenvalues()
    ref = A.eigh()
    B = alg.operator(A.matrix.clone(), is_self_adjoint=True).cache_spectrum()
    for _ in range(2):
        got = B.eigh()
        got_w = B.eigenvalues()
        assert len(got) == len(ref) == alg.C
        for c, k in enumerate(alg.k_factors):
            w, V = got[c]
            assert w.shape == ref[c][0].shape and V.shape == ref[c][1].shape
            assert torch.allclose(w, ref[c][0], atol=1e-12)
            assert torch.allclose(got_w[c], ref_w[c], atol=1e-12)
            blk = A.matrix[:, c, :k, :k]
            assert torch.allclose(V @ torch.diag_embed(w).to(V.dtype) @ V.mH, blk, atol=1e-12)




# kills: algebra#270
def test_double_precision_property_tolerance_is_1e_10_relative():
    """Double precision: a relative anti-Hermitian part of 1e-8 is detected by is_self_adjoint and by validate=True."""
    alg = TypeIAlgebra([3], [3], precision='double')
    A = _herm(alg, seed=2, batch=1).matrix
    E = torch.zeros_like(A)
    E[0, 0, 0, 1] = 1e-8 * A.abs().amax()
    M = A + E                                                  # non-Hermitian at relative level 1e-8
    assert alg.operator(M).is_self_adjoint is False
    with pytest.raises(ValueError):
        alg.operator(M, validate=True, is_self_adjoint=True)
    assert alg.operator(A).is_self_adjoint is True


# kills: algebra#696
def test_unary_minus_of_operator():
    """-A has matrix -A.matrix, keeps self-adjointness, and is not flagged positive."""
    alg = TypeIAlgebra([2, 3], [2, 3], precision='double')
    A = _herm(alg, seed=3, batch=2)
    N = -A
    assert torch.equal(N.matrix, -A.matrix)
    assert N.is_self_adjoint
    P = alg.identity()
    assert P.is_positive
    assert torch.allclose((-P).lambda_max, -torch.ones(1, dtype=torch.float64))
    assert not (-P).is_positive




# kills: algebra#718
def test_cast_complex_to_real_tolerance_is_relative_to_the_entries():
    """cast to a real algebra accepts an imaginary part below 1e-6 * max(1, max |entry|) (here 1e-5 for entries ~1e3)."""
    algc = TypeIAlgebra([2], [2], precision='double')
    algr = TypeIAlgebra([2], [2], complex_valued=False, precision='double')
    M = torch.tensor([[1000.0, 2.0], [2.0, -500.0]], dtype=torch.complex128)
    M = M + 1e-5j * torch.tensor([[0.0, 1.0], [-1.0, 0.0]], dtype=torch.complex128)
    X = algc.from_blocks([M])
    Y = X.cast(algr)
    assert not torch.is_complex(Y.matrix)
    assert torch.allclose(Y.matrix[0, 0], M.real)
    with pytest.raises(ValueError):
        algc.from_blocks([M + 1e-2j * torch.eye(2, dtype=torch.complex128)]).cast(algr)


# kills: algebra#390
def test_is_invertible_boundary_singular_value_equal_to_tolerance():
    """is_invertible: a smallest singular value equal to the tolerance 1e-10 max|entry| counts as zero (as in validation)."""
    alg = TypeIAlgebra([2], [2], precision='double')
    A = alg.from_blocks([torch.diag(torch.tensor([1.0, 1e-10], dtype=torch.float64))])
    assert A.is_invertible is False
    B = alg.from_blocks([torch.diag(torch.tensor([1.0, 2e-10], dtype=torch.float64))])
    assert B.is_invertible is True


# kills: algebra#386
def test_is_positive_setter_flag_is_validated_on_materialisation():
    """A user flag set with op.is_positive = True on a lazy operator is checked when the operator is materialised."""
    alg = TypeIAlgebra([2], [2], precision='double')
    M = torch.diag(torch.tensor([1.0, -1.0], dtype=torch.complex128)).reshape(1, 1, 2, 2)
    op = tv.Operator(alg, generator=lambda: M.clone())
    op.is_positive = True
    with pytest.raises(ValueError):
        op.matrix


# kills: algebra#782
def test_clear_cache_drops_the_spectral_cache():
    """clear_cache() empties the cached decomposition (cache_info reports nothing cached, 0 bytes)."""
    alg = TypeIAlgebra([2, 3], [2, 3], precision='double')
    A = _herm(alg, seed=4, batch=1).cache_spectrum()
    A.eigenvalues()
    info = A.cache_info()
    assert info['cached'] and info['bytes'] > 0
    A.clear_cache()
    info = A.cache_info()
    assert info['enabled'] and not info['cached'] and info['bytes'] == 0


# kills: algebra#864
def test_power_result_flags_known_without_evaluation():
    """A.power(p) is flagged self-adjoint by construction: querying the flag does not evaluate the lazy result (here
    under no_grad, which would otherwise drop the gradient)."""
    alg = TypeIAlgebra([2], [2], precision='double')
    M = torch.tensor([[2.0, 0.5], [0.5, 1.0]], dtype=torch.complex128, requires_grad=True)
    A = alg.operator(M.reshape(1, 1, 2, 2), is_self_adjoint=True, is_positive=True)
    P = A.power(0.5)
    with torch.no_grad():
        assert P.is_self_adjoint and P.is_positive
    P.trace.real.sum().backward()
    assert M.grad is not None


# kills: algebra#1207
def test_central_projection_flag():
    """alg.central: is_projection only for real values in {0, 1}; a real central element 2*1 is not a projection."""
    alg = TypeIAlgebra([2, 1], [2, 1], precision='double')
    assert alg.central([1.0, 0.0]).is_projection is True
    assert alg.central([2.0, 1.0]).is_projection is False
    assert alg.central([1j, 0.0]).is_projection is False




# kills: algebra#1097
def test_real_algebra_accepts_complex_input_with_zero_imaginary_part():
    """A real algebra converts complex blocks with negligible imaginary part to their real part."""
    alg = TypeIAlgebra([2], [2], complex_valued=False, precision='double')
    M = torch.tensor([[1.0, 2.0], [2.0, 3.0]], dtype=torch.complex128)
    A = alg.from_blocks([M])
    assert not torch.is_complex(A.matrix)
    assert torch.equal(A.matrix[0, 0], M.real)


# kills: algebra#514
def test_inv_property_is_the_cached_inverse():
    """A.inv equals A.inverse() (A^{-1} A = 1) and is computed once (the same object on repeated access)."""
    alg = TypeIAlgebra([2, 3], [2, 3], precision='double')
    A = _herm(alg, seed=5, batch=1) + alg.identity() * 10.0
    I1 = A.inv
    assert I1 is A.inv
    assert torch.allclose((I1 @ A).matrix, alg.identity().matrix, atol=1e-12)


# kills: cost#2433, cost#2437
def test_cost_eta_formatting_units():
    """ETA strings of CostWarning: seconds below 2 min, minutes below 2 h, hours from 2 h on (private helper)."""
    from torch_vn_algebra.cost import _fmt_seconds
    assert _fmt_seconds(30) == "30 s"
    assert _fmt_seconds(1800) == "30 min"
    assert _fmt_seconds(7199) == "120 min"
    assert _fmt_seconds(7200) == "2.0 h"


# kills: jordan#6688
def test_jordan_star_algebra_method_keeps_unital_flag():
    """J.star_algebra() uses J's `unital` flag: for x = diag(1, 2, 0) without unit, dim A_sa = 2 (no identity)."""
    from torch_vn_algebra import jordan as jd
    x = torch.diag(torch.tensor([1.0, 2.0, 0.0], dtype=torch.float64))
    J = jd.jordan_algebra_generated([x], unital=False)
    A = J.star_algebra()
    assert J.dim == 2 and A.dim == 2 and A.unital is False
    assert jd.jordan_algebra_generated([x]).star_algebra().dim == 3


# kills: su2#8242
def test_su2_full_spectrum_of_derived_real_hamiltonian():
    """full_spectrum of an operator without recorded field (derived by arithmetic) in a real SU2Chain: field 0."""
    chain = tv.SU2Chain(4, complex_valued=False)
    H = chain.heisenberg(J=1.0)
    H2 = H * 2.0
    E = chain.full_spectrum(H)
    E2 = chain.full_spectrum(H2)
    assert E2.shape == (1, 16)
    assert torch.allclose(E2, 2 * E, atol=1e-12)




# kills: functional#5328
def test_holomorphic_error_estimate_tracks_cancellation():
    """functional.holomorphic: the returned error estimate includes the cancellation (rounding) term of the contour
    sum: for exp of a block with spectrum {-40, 0} (terms ~e^40 larger than the result) it is within a factor 2.5 of
    the true error."""
    import warnings
    import scipy.linalg as sl
    from torch_vn_algebra import functional as fc
    X = torch.tensor([[[-40.0, 1.0], [0.0, 0.0]]], dtype=torch.complex128)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        F, err = fc.holomorphic(X, torch.exp)
    ref = sl.expm(X[0].numpy())
    actual = np.abs(F[0].numpy() - ref).max() / np.abs(ref).max()
    assert actual / 2.5 <= err < 1e-5


# kills: divergences#4140
def test_petz_recovery_rejects_operator_of_wrong_algebra():
    """PetzRecovery(...)(X) raises ValueError when X is an operator of another algebra than the channel output."""
    from torch_vn_algebra import channels as ch, states as st
    alg = TypeIAlgebra([2, 3], [2, 3], precision='double')
    other = TypeIAlgebra([3, 2], [3, 2], precision='double')
    torch.manual_seed(0)
    sigma = st.random_density_matrix(alg)
    E = ch.random_channel(alg, 2)
    R = dv.PetzRecovery(sigma, E)
    R(alg.identity())                                  # the right algebra works
    with pytest.raises(ValueError):
        R(other.identity())


# kills: inclusions#6206
def test_chain_tower_needs_two_levels():
    """inclusions.chain_tower(L_max, L_min) requires L_min < L_max (ValueError for L_min == L_max)."""
    from torch_vn_algebra import inclusions as inc
    with pytest.raises(ValueError):
        inc.chain_tower(2, 2)
    assert len(inc.chain_tower(2, 1)) == 1


# kills: composite#2218
def test_partial_trace_of_non_product_algebra_raises_value_error():
    """composite.partial_trace on an algebra not built by tensor_product raises ValueError (not a TypeError)."""
    alg = TypeIAlgebra([2], [2], precision='double')
    with pytest.raises(ValueError):
        tv.partial_trace(alg.identity(), keep=1)


# kills: channels#1746
def test_choi_checks_memory_before_allocating():
    """Channel.choi estimates its memory first: with max_memory_fraction tiny it raises InsufficientMemoryError."""
    from torch_vn_algebra import channels as ch, cost
    alg = TypeIAlgebra([12], [12], precision='double')
    E = ch.random_channel(alg, 1, batch_size=250)          # Choi matrices: 250 * 12^4 * 16 bytes ~ 83 MB
    old = cost.get_limits()
    try:
        cost.set_limits(max_memory_fraction=1e-9, warn_memory_fraction=1e-9)
        with pytest.raises(cost.InsufficientMemoryError):
            E.choi()
    finally:
        cost.set_limits(**old)




# kills: decompose#3042
def test_real_structure_ok_with_zero_tolerance_on_exact_input():
    """RealStructure.ok holds with tol=0 when every residual is exactly 0 (diagonal real generators: ok iff <= tol)."""
    gens = [torch.diag(torch.tensor([1.0, 2.0, 3.0], dtype=torch.float64))]
    rs = tv.wedderburn(gens).real_structure(tol=0.0)
    assert all(v == 0.0 for v in rs.residuals.values())
    assert rs.ok is True
    assert rs.types == ('R', 'R', 'R')


# kills: hilbert_space#5494
def test_basis_does_not_alias_callers_tensor():
    """HilbertSpace.Basis(parent, V) copies V: changing the caller's tensor afterwards leaves the basis unchanged."""
    H = tv.HilbertSpace(n=5, k=3, batch_size=2, n_channels=2, complex_valued=True, device='cpu')
    V = torch.zeros(2, 2, 5, 3, dtype=H.dtype)
    for i in range(3):
        V[..., i, i] = 1.0
    basis = H.Basis(H, V=V)
    V.mul_(2.0)
    assert basis.is_orthonormal()
    assert torch.equal(basis.V[..., 0, 0], torch.ones(2, 2, dtype=H.dtype))


# kills: optimize#7306
def test_positive_param_recenter_moves_reference_point():
    """PositiveParam.recenter(): the operator is unchanged and the unitary generator G is reset to 0."""
    from torch_vn_algebra import optimize as opt
    torch.manual_seed(0)
    alg = TypeIAlgebra([2, 3], [2, 3], precision='double')
    P = opt.PositiveParam(alg, 2)
    with torch.no_grad():
        P.unitary.G.add_(0.1 * torch.randn_like(P.unitary.G))
    X0 = P.operator().matrix.detach().clone()
    P.recenter()
    assert float(P.unitary.G.detach().abs().max()) == 0.0
    assert torch.allclose(P.operator().matrix.detach(), X0, atol=1e-12)




# kills: cost#2482
def test_kernel_flops_cubic_in_block_size():
    """cost.kernel_flops: c k^3 flops per k x k kernel, times 4 for complex (eigh: c = 10)."""
    from torch_vn_algebra import cost
    assert cost.kernel_flops('eigh', 10, False) == 10 * 10 ** 3
    assert cost.kernel_flops('matmul', 8, True) == 2 * 8 ** 3 * 4


# kills: functional#5036, functional#5033
def test_use_dilation_auto_only_on_gpu():
    """functional.use_dilation: 'auto' is False on the CPU for every block size (batched CPU SVD); 'eigh' / 'svd' force it."""
    from torch_vn_algebra import functional as fc
    assert fc.use_dilation(torch.zeros(1, 40, 40)) is False
    assert fc.use_dilation(torch.zeros(1, 4, 4)) is False
    assert fc.use_dilation(torch.zeros(1, 40, 40), 'eigh') is True
    assert fc.use_dilation(torch.zeros(1, 40, 40), 'svd') is False


# kills: chains#1346
def test_xxz_default_couplings_are_one():
    """SpinChain.xxz() defaults to J = 1 and Delta = 1 (the isotropic Heisenberg chain)."""
    chain = tv.SpinChain(4)
    H0 = chain.xxz()
    assert torch.allclose(H0.matrix, chain.xxz(J=1.0, Delta=1.0).matrix)
    assert not torch.allclose(H0.matrix, chain.xxz(J=1.0, Delta=2.0).matrix)
    assert torch.allclose(2 * H0.matrix, chain.xxz(J=2.0).matrix)


# kills: states#7454
def test_from_density_does_not_normalise_by_default():
    """DensityMatrix.from_density: a density of trace != 1 is rejected unless normalize=True (no silent rescaling)."""
    alg = TypeIAlgebra([2, 1], [2, 1], precision='double')
    d = alg.identity().matrix / 3.0
    rho = tv.DensityMatrix.from_density(alg, d)
    assert abs(float(rho.trace.real) - 1.0) < 1e-12
    with pytest.raises(ValueError):
        tv.DensityMatrix.from_density(alg, 2.0 * d)
    rho_n = tv.DensityMatrix.from_density(alg, 2.0 * d, normalize=True)
    assert torch.allclose(rho_n.matrix, rho.matrix)


# kills: channels#1756
def test_choi_and_superoperator_of_empty_sector_are_empty():
    """Channel.choi() / superoperator() give empty (batch, 0, 0) tensors for sectors with k_c = 0."""
    from torch_vn_algebra import channels as ch
    torch.manual_seed(0)
    alg = TypeIAlgebra([2, 1, 3], [2, 0, 3], precision='double')
    E = ch.random_channel(alg, 2, batch_size=3)
    for blocks in (E.choi(), E.superoperator()):
        assert [tuple(b.shape) for b in blocks] == [(3, 4, 4), (3, 0, 0), (3, 9, 9)]


# kills: su2#7638
def test_su2_chain_length_limit():
    """SU2Chain refuses L > 18 (su2.MAX_L) with a ValueError, before building anything."""
    from torch_vn_algebra import su2
    assert su2.MAX_L == 18
    with pytest.raises(ValueError, match="L > 18"):
        tv.SU2Chain(19)


# kills: su2#8291
def test_su2_lowest_energies_requires_operator_of_all_spins():
    """SU2Chain.lowest_energies(H) raises ValueError for an operator of a single-spin algebra."""
    chain = tv.SU2Chain(4)
    Hs = chain.heisenberg(spin=0)
    with pytest.raises(ValueError, match="lowest_energies"):
        chain.lowest_energies(Hs)


# kills: composite#2237
def test_kron_checks_memory_before_allocating():
    """composite.kron estimates the memory of A (x) B first: with max_memory_fraction tiny it raises."""
    from torch_vn_algebra import cost
    a1 = TypeIAlgebra([4], [4], precision='double')
    a12 = tv.tensor_product(a1, a1)
    A = a1.identity(batch_size=6000)
    old = cost.get_limits()
    try:
        cost.set_limits(max_memory_fraction=1e-9, warn_memory_fraction=1e-9)
        with pytest.raises(cost.InsufficientMemoryError):
            tv.kron(A, A, a12).matrix
    finally:
        cost.set_limits(**old)




# kills: channels#2036
def test_inter_sector_channel_zeroes_kraus_padding_when_only_output_is_uneven():
    """InterSectorChannel zeroes the padding of the Kraus tensor whenever either algebra has uneven sectors."""
    from torch_vn_algebra import InterSectorChannel
    a_in = TypeIAlgebra([2], [2], precision='double')                 # no padding
    a_out = TypeIAlgebra([2, 1], [2, 1], precision='double')          # sector 1 padded (1 of 2 rows)
    K = torch.ones(1, 1, 2, 1, 2, 2, dtype=torch.complex128)
    E = InterSectorChannel(a_in, a_out, K)
    assert float(E.kraus[0, 0, 1, 0, 1:, :].abs().max()) == 0.0
    assert float(E.kraus[0, 0, 0, 0].abs().min()) == 1.0
    rho = a_in.identity() * 0.5
    out = E.apply(rho)
    assert float(out.matrix[0, 1, 1:, :].abs().max()) == 0.0 and float(out.matrix[0, 1, :, 1:].abs().max()) == 0.0


# kills: divergences#3688
def test_umegaki_rtol_zero_detects_support_leak_onto_exact_kernel():
    """D(rho||sigma) = inf with rtol=0 when rho has weight on the exact kernel of sigma (support projection on w > 0)."""
    alg = TypeIAlgebra([3], [3], precision='double')
    rho = _diag_op(alg, [0.5, 0.5, 0.0])
    sigma = _diag_op(alg, [1.0, 0.0, 0.0])
    assert math.isinf(float(dv.umegaki_relative_entropy(rho, sigma, rtol=0.0)))
    assert math.isinf(float(dv.umegaki_relative_entropy(rho, sigma)))


# kills: inclusions#5751
def test_inclusion_embed_keeps_known_batch_lazy():
    """Inclusion.embed(x) knows its batch without evaluating (shape under no_grad keeps the gradient of x)."""
    from torch_vn_algebra import inclusions as inc
    N = TypeIAlgebra([2, 1], [2, 1], precision='double')
    I = inc.from_tensor_product(N, 2)
    s = torch.tensor(3.0, dtype=torch.float64, requires_grad=True)
    y = I.embed(N.central([s, 1.0]))
    with torch.no_grad():
        assert tuple(y.shape) == (1, I.M.C, I.M.k_max, I.M.k_max)
    y.trace.real.sum().backward()
    assert s.grad is not None and abs(float(s.grad) - 4.0) < 1e-12




# kills: cost#2507
def test_batched_call_counts_items_of_a_3d_stack():
    """cost.batched_call: a (n, k, k) stack has n items of one k x k kernel each (n * 10 k^3 flops for eigh): below
    probe_flops nothing is probed (no CostWarning even with warn_seconds = 0)."""
    import warnings
    from torch_vn_algebra import cost
    x = torch.eye(8, dtype=torch.float64).expand(64, 8, 8).contiguous()     # 64 * 10 * 8^3 = 327680 flops
    old = cost.get_limits()
    try:
        cost.set_limits(probe_flops=1e6, warn_seconds=0.0)
        with warnings.catch_warnings():
            warnings.simplefilter('error', cost.CostWarning)
            w, V = cost.batched_call(lambda a: torch.linalg.eigh(a), x, "eigh test")
    finally:
        cost.set_limits(**old)
    assert torch.allclose(w, torch.ones(64, 8, dtype=torch.float64))


class _DenseH:
    """Minimal Hamiltonian object for krylov (matvec, batch_size, dim, dtype, device) from a dense (B, k, k)."""

    def __init__(self, M):
        self.M, self.batch_size, self.dim, self.dtype, self.device = M, M.shape[0], M.shape[-1], M.dtype, M.device

    def matvec(self, V):
        return torch.einsum('bij,bj->bi', self.M, V.to(self.M.dtype))


# kills: krylov#7127, krylov#7142
def test_expm_multiply_step_error_estimate_vanishes_on_full_krylov_space():
    """krylov.expm_multiply_step: with m = dim the Krylov space is invariant, the step is exact and the error
    estimate (next Lanczos residual) is ~0."""
    from torch_vn_algebra import krylov
    g = torch.Generator().manual_seed(0)
    A = torch.randn(1, 4, 4, dtype=torch.float64, generator=g)
    M = (A + A.mT).to(torch.complex128)
    psi = torch.randn(1, 4, dtype=torch.complex128, generator=g)
    out, err = krylov.expm_multiply_step(_DenseH(M), psi, 0.3, m=4)
    ref = torch.linalg.matrix_exp(-0.3j * M[0]) @ psi[0]
    assert torch.allclose(out[0], ref, atol=1e-10)
    assert float(err.max()) < 1e-10


# kills: cost#2431, cost#2424
def test_memory_error_message_reports_mebibytes():
    """cost.check_memory: the InsufficientMemoryError message gives the size in binary units (100 * 2^20 B = 100.0 MB)."""
    from torch_vn_algebra import cost
    old = cost.get_limits()
    try:
        cost.set_limits(max_memory_fraction=1e-12)
        with pytest.raises(cost.InsufficientMemoryError, match=r"needs about 100\.0 MB"):
            cost.check_memory(100 * 2 ** 20, 'cpu', "test allocation")
    finally:
        cost.set_limits(**old)


# kills: cost#2647, cost#2648
def test_step_timer_close_closes_progress_bar():
    """cost.StepTimer(progress=True).close() closes its tqdm progress bar."""
    pytest.importorskip("tqdm")
    from torch_vn_algebra import cost
    t = cost.StepTimer(3, "test loop", progress=True)
    for _ in range(3):
        t.step()
    assert t.bar is not None and not t.bar.disable
    t.close()
    assert t.bar.disable


# kills: dynamics#4389
def test_solve_operator_ode_is_silent_by_default(capsys):
    """dynamics.solve_operator_ode shows no progress output unless progress=True."""
    pytest.importorskip("tqdm")
    from torch_vn_algebra import dynamics as dy
    alg = TypeIAlgebra([2], [2], precision='double')
    X0 = alg.identity()
    Xs = dy.solve_operator_ode(lambda t, X: -X, X0, [0.0, 0.1], substeps=2)
    assert len(Xs) == 2
    captured = capsys.readouterr()
    assert captured.err == "" and captured.out == ""


# kills: composite#2206
def test_tensor_product_ambient_dimensions_multiply():
    """composite.tensor_product: n_(c,d) = n_c n_d and k_(c,d) = k_c k_d."""
    a1 = TypeIAlgebra([3, 2], [2, 1], precision='double')
    a2 = TypeIAlgebra([4], [3], precision='double')
    a12 = tv.tensor_product(a1, a2)
    assert a12.n_factors == [12, 8]
    assert a12.k_factors == [6, 3]




# kills: hilbert_space#5536
def test_is_orthonormal_default_tolerance():
    """Basis.is_orthonormal() uses tol = 1e-6 by default: a frame with a 1e-4 overlap is not orthonormal."""
    H = tv.HilbertSpace(n=5, k=3, batch_size=1, n_channels=1, complex_valued=True, device='cpu', precision='double')
    V = torch.zeros(1, 1, 5, 3, dtype=H.dtype)
    for i in range(3):
        V[..., i, i] = 1.0
    assert H.Basis(H, V=V).is_orthonormal()
    V[..., 0, 1] = 1e-4
    B = H.Basis(H, V=V)
    assert not B.is_orthonormal()
    assert B.is_orthonormal(tol=1e-3)


# kills: su2#7836
def test_su2_sector_index_rejects_non_total_spin_with_message():
    """SU2Chain.sector_index(S) raises ValueError naming S when S is not a total spin of the chain."""
    chain = tv.SU2Chain(4)
    assert chain.sector_index(0) == chain.two_spins.index(0)
    with pytest.raises(ValueError, match="is not a total spin"):
        chain.sector_index(3)


# kills: jordan#6846
def test_spin_factor_power_of_unbatched_and_batched_elements():
    """SpinFactor.power(x, p) equals repeated Jordan products, for x of shape (m+1,) and (2, 3, m+1)."""
    from torch_vn_algebra import jordan as jd
    V = jd.SpinFactor(3)
    g = torch.Generator().manual_seed(0)
    for shape in ((4,), (2, 3, 4)):
        x = torch.randn(*shape, dtype=torch.float64, generator=g)
        x3 = V.product(x, V.product(x, x))
        p3 = V.power(x, 3)
        assert p3.shape == x.shape
        assert torch.allclose(p3, x3, atol=1e-12)
        assert torch.allclose(V.power(x, 0), V.unit(shape[:-1]))


# kills: composite#2250
def test_kron_positive_flag_needs_both_factors_positive():
    """kron(A, B) is flagged positive only if both A and B are positive: kron(1, -1) is not positive."""
    a1 = TypeIAlgebra([2], [2], precision='double')
    a12 = tv.tensor_product(a1, a1)
    I = a1.identity()
    assert tv.kron(I, I, a12).is_positive
    assert tv.kron(I, -I, a12).is_positive is False
    assert tv.kron(-I, I, a12).is_positive is False


# kills: su2#8009
def test_su2_lower_accepts_a_single_vector():
    """SU2Chain.lower(x, N) for a vector x equals lower(x[:, None], N)[:, 0]."""
    chain = tv.SU2Chain(4)
    g = torch.Generator().manual_seed(1)
    x = torch.randn(math.comb(4, 2), dtype=torch.float64, generator=g)
    y = chain.lower(x, 2)
    assert y.dim() == 1
    assert torch.allclose(y, chain.lower(x[:, None], 2)[:, 0])


# kills: jordan#6959
def test_spin_factor_symmetries_follow_requested_dtype():
    """SpinFactor.symmetries(dtype) returns the anticommuting symmetries in the requested dtype (cache per dtype)."""
    from torch_vn_algebra import jordan as jd
    V = jd.SpinFactor(3)
    s64 = V.symmetries(torch.complex64)
    s128 = V.symmetries(torch.complex128)
    assert s64.dtype == torch.complex64 and s128.dtype == torch.complex128
    assert torch.allclose(s64.to(torch.complex128), s128, atol=1e-6)
    assert V.symmetries(torch.complex64).dtype == torch.complex64


# kills: composite#2349, composite#2300
def test_partial_trace_channel_with_empty_sector_in_second_factor():
    """partial_trace_channel works when the second factor has an empty sector, and agrees with partial_trace."""
    from torch_vn_algebra import states as st
    torch.manual_seed(0)
    a1 = TypeIAlgebra([2], [2], precision='double')
    a2 = TypeIAlgebra([2, 1], [2, 0], precision='double')
    a12 = tv.tensor_product(a1, a2)
    rho = st.random_density_matrix(a12)
    E = tv.composite.partial_trace_channel(a12, keep=1)
    assert torch.allclose(E.apply(rho).matrix, tv.partial_trace(rho, keep=1).matrix, atol=1e-12)




# kills: chains#1478
def test_random_field_heisenberg_fields_uniform_in_minus_w_w():
    """SpinChain.random_field_heisenberg: fields h_i = (2 u - 1) W with u ~ U[0, 1) from the given generator, i.e.
    h ~ U[-W, W]; equal to xxz(J, 1, h)."""
    chain = tv.SpinChain(4)
    W = 1.5
    H = chain.random_field_heisenberg(W, 3, generator=torch.Generator().manual_seed(7))
    u = torch.rand(3, 4, generator=torch.Generator().manual_seed(7), dtype=torch.float64)
    h = (2 * u - 1) * W
    assert float(h.min()) < 0 < float(h.max())
    assert torch.allclose(H.matrix, chain.xxz(1.0, 1.0, h).matrix, atol=1e-12)


# kills: states#7499
def test_condition_on_default_eps_keeps_small_probabilities():
    """DensityMatrix.condition_on(P): outcomes with probability 1e-10 (> default eps = 1e-12) are conditioned on."""
    alg = TypeIAlgebra([2], [2], precision='double')
    rho = tv.DensityMatrix(alg, matrix=torch.diag(torch.tensor([1 - 1e-10, 1e-10], dtype=torch.complex128))
                           .reshape(1, 1, 2, 2))
    P = _diag_op(alg, [0.0, 1.0])
    post, prob = rho.condition_on(P)
    assert abs(float(prob) - 1e-10) < 1e-20
    assert abs(float(post.trace.real) - 1.0) < 1e-9


# kills: inclusions#5662
def test_inclusion_batch_size_without_unitary_is_one():
    """Inclusion.batch_size is 1 without a unitary; embed then keeps the batch of x."""
    from torch_vn_algebra import inclusions as inc
    N = TypeIAlgebra([2], [2], precision='double')
    I = inc.from_tensor_product(N, 2)
    assert I.batch_size == 1
    assert I.embed(N.identity(batch_size=3)).matrix.shape[0] == 3
    assert I.embed(N.identity()).shape[0] == 1


# kills: divergences#4215
def test_lp_norm_rejects_weights_of_wrong_batch():
    """divergences.lp_norm: weights of shape (batch, C) must broadcast against A's batch (ValueError otherwise)."""
    alg = TypeIAlgebra([2, 3], [2, 3], precision='double')
    A = _herm(alg, seed=6, batch=2)
    ok = dv.lp_norm(A, 2.0, weights=torch.ones(2, 2, dtype=torch.float64))
    assert ok.shape == (2,)
    with pytest.raises(ValueError):
        dv.lp_norm(A, 2.0, weights=torch.ones(3, 2, dtype=torch.float64))


# kills: channels#1761
def test_superoperator_memory_estimate_is_k_squared():
    """Channel.superoperator estimates B k^4 entries: 100 blocks of size 12 (33 MB) stay below the 64 MB floor of
    cost.check_memory, so it runs even with max_memory_fraction tiny."""
    from torch_vn_algebra import channels as ch, cost
    torch.manual_seed(0)
    alg = TypeIAlgebra([12], [12], precision='double')
    E = ch.random_channel(alg, 1, batch_size=100)
    old = cost.get_limits()
    try:
        cost.set_limits(max_memory_fraction=1e-9, warn_memory_fraction=1e-9)
        S = E.superoperator()
    finally:
        cost.set_limits(**old)
    assert S[0].shape == (100, 144, 144)




# kills: algebra#332
def test_is_self_adjoint_flag_known_without_evaluation():
    """Operator.is_self_adjoint returns a known flag without evaluating a lazy operator (under no_grad the gradient
    of alg.central with tensor values would otherwise be lost)."""
    alg = TypeIAlgebra([2, 3], [2, 3], precision='double')
    s = torch.tensor(2.0, dtype=torch.float64, requires_grad=True)
    op = alg.central([s, 1.0])
    with torch.no_grad():
        assert op.is_self_adjoint is True
    op.trace.real.sum().backward()
    assert s.grad is not None and abs(float(s.grad) - 2.0) < 1e-12


# kills: optimize#7368
def test_extremize_ignores_nan_values_when_maximizing():
    """optimize.extremize(maximize=True): a NaN objective value at some iterate is ignored (not taken as +inf)."""
    from torch_vn_algebra import optimize as opt
    torch.manual_seed(0)
    alg = TypeIAlgebra([2], [2], precision='double')
    P = opt.PositiveParam(alg, 2)
    calls = [0]

    def objective():
        calls[0] += 1
        v = P.operator().matrix[:, 0, 0, 0].real
        if calls[0] == 1:                               # one NaN value for start 0 (gradient unaffected)
            v = v + torch.tensor([float('nan'), 0.0], dtype=v.dtype)
        return v
    best = opt.extremize(objective, [P], maximize=True, rounds=1, steps=3)
    assert torch.isfinite(best).all()


# kills: optimize#7385
def test_extremize_discrete_steps_only_between_rounds():
    """optimize.extremize: SelfAdjointParam.discrete_step runs between gradient rounds only (rounds - 1 times):
    2 rounds x 3 steps + 1 discrete step (1 + 8 trial evaluations) + 1 final evaluation = 16 objective calls."""
    from torch_vn_algebra import optimize as opt
    torch.manual_seed(0)
    alg = TypeIAlgebra([2], [2], precision='double')
    P = opt.SelfAdjointParam(alg, 2)
    calls = [0]

    def objective():
        calls[0] += 1
        return P.operator().matrix[:, 0, 0, 0].real
    opt.extremize(objective, [P], maximize=True, rounds=2, steps=3)
    assert calls[0] == 2 * 3 + (1 + 8) + 1


# kills: channels#1694
def test_channel_mix_rejects_batched_p():
    """Channel.mix(other, p) with a tensor p of several entries raises ValueError (asking for a single number)."""
    from torch_vn_algebra import channels as ch
    torch.manual_seed(0)
    alg = TypeIAlgebra([2], [2], precision='double')
    E, F = ch.random_channel(alg, 2), ch.random_channel(alg, 2)
    E.mix(F, torch.tensor(0.3))
    with pytest.raises(ValueError, match="single number"):
        E.mix(F, torch.tensor([0.1, 0.2]))


# kills: channels#1742
def test_channel_single_precision_trace_preservation_tolerance():
    """Channel.is_trace_preserving() in single precision uses tol 1e-5: a deviation of 1e-3 is not trace preserving."""
    from torch_vn_algebra import channels as ch
    alg = TypeIAlgebra([2], [2])                                     # single precision
    E = ch.identity_channel(alg)
    assert E.is_trace_preserving()
    K = E.kraus * (1 + 1e-3)
    F = tv.Channel(alg, K)
    assert float(F.trace_preservation_error()) > 1e-3
    assert not F.is_trace_preserving()


# kills: inclusions#5860
def test_commutator_error_draws_batch_size_samples():
    """Inclusion.commutator_error(other, batch_size) returns batch_size samples for unbatched inclusions (~0 for the
    relative commutant)."""
    from torch_vn_algebra import inclusions as inc
    torch.manual_seed(0)
    N = TypeIAlgebra([2], [2], precision='double')
    I = inc.from_tensor_product(N, 2)
    e = I.commutator_error(I.relative_commutant(), batch_size=5)
    assert e.shape == (5,)
    assert float(e.max()) < 1e-12


# kills: jordan#6813
def test_spin_factor_split_checks_dimension():
    """SpinFactor.split(x) raises ValueError when the last dimension is not m + 1."""
    from torch_vn_algebra import jordan as jd
    V = jd.SpinFactor(3)
    a, u = V.split(torch.arange(4.0, dtype=torch.float64))
    assert float(a) == 0.0 and u.shape == (3,)
    with pytest.raises(ValueError):
        V.split(torch.zeros(5, dtype=torch.float64))




# kills: jordan#6698
def test_star_algebra_self_adjoint_part_is_unital_by_default():
    """jordan.star_algebra_self_adjoint_part adjoins 1 by default: for x = diag(1, 2, 0), dim A_sa = 3 (2 without)."""
    from torch_vn_algebra import jordan as jd
    x = torch.diag(torch.tensor([1.0, 2.0, 0.0], dtype=torch.float64))
    assert jd.star_algebra_self_adjoint_part([x]).dim == 3
    assert jd.star_algebra_self_adjoint_part([x], unital=False).dim == 2


# kills: dynamics#4263
def test_schrodinger_rk4_masks_padding_of_raw_state_with_operator_hamiltonian():
    """dynamics.schrodinger_rk4: with H an Operator of an algebra with uneven sectors, the padding entries of a raw
    psi0 are dropped (they are not part of the state), so the result has zero padding."""
    from torch_vn_algebra import dynamics as dy
    alg = TypeIAlgebra([2, 1], [2, 1], precision='double')
    H = _herm(alg, seed=8, batch=1)
    psi0 = torch.ones(1, 2, 2, dtype=torch.complex128)            # entry [0, 1, 1] is padding
    out = dy.schrodinger_rk4(psi0, H, [0.0, 0.1], substeps=2)
    assert out.shape == (2, 1, 2, 2)
    assert float(out[:, 0, 1, 1].abs().max()) == 0.0


# kills: su2#8022
def test_su2_paths_rejects_non_total_spin():
    """SU2Chain.paths(S) raises ValueError ('is not a total spin') for S not a total spin of the chain."""
    chain = tv.SU2Chain(4)
    assert len(chain.paths(0)) == chain.multiplicity(0)
    with pytest.raises(ValueError, match="is not a total spin"):
        chain.paths(3)


# kills: jordan#6636
def test_jordan_contains_is_inclusive_at_tolerance():
    """JordanAlgebra.contains(x, tol): residual <= tol (inclusive): a basis element with residual exactly 0 is
    contained with tol = 0."""
    from torch_vn_algebra import jordan as jd
    x = torch.diag(torch.tensor([1.0, 2.0, 0.0], dtype=torch.float64))
    J = jd.jordan_algebra_generated([x])
    b = J.basis[0]
    assert float(J.residual(b)) == 0.0
    assert bool(J.contains(b, tol=0.0))




# kills: cost#2451
def test_available_memory_is_in_bytes():
    """cost.available_memory('cpu') converts MemAvailable from kB (1024 bytes) on Linux: a multiple of 1024."""
    import os
    from torch_vn_algebra import cost
    if not os.path.exists('/proc/meminfo'):
        pytest.skip("needs /proc/meminfo")
    avail = cost.available_memory('cpu')
    assert avail is not None and avail > 0
    assert avail % 1024 == 0


# kills: chains#1244
def test_spin_chain_length_limit():
    """SpinChain refuses L > 30 (chains.MAX_L) with a ValueError before allocating the 2^L labels."""
    from torch_vn_algebra import chains
    assert chains.MAX_L == 30
    with pytest.raises(ValueError):
        tv.SpinChain(31)


# kills: optimize#7351
def test_discrete_step_keeps_signs_on_ties_when_minimizing():
    """SelfAdjointParam.discrete_step keeps a sign flip only if it strictly improves the objective (minimize: a
    constant objective leaves the signs unchanged)."""
    from torch_vn_algebra import optimize as opt
    torch.manual_seed(0)
    alg = TypeIAlgebra([2, 3], [2, 3], precision='double')
    P = opt.SelfAdjointParam(alg, 4)
    s0 = P.signs.clone()
    for maximize in (False, True):
        P.discrete_step(lambda: torch.zeros(4, dtype=torch.float64), maximize=maximize)
        assert torch.equal(P.signs, s0)


# kills: states#7459
def test_density_matrix_relative_entropy_default_clamp():
    """DensityMatrix.relative_entropy clamps the eigenvalues of sigma at eps = 1e-12 by default: for a support
    mismatch D = 0.5 log(1e12) - log 2 (large but finite)."""
    alg = TypeIAlgebra([2], [2], precision='double')
    rho = tv.DensityMatrix(alg, matrix=torch.diag(torch.tensor([0.5, 0.5], dtype=torch.complex128)).reshape(1, 1, 2, 2))
    sigma = tv.DensityMatrix(alg, matrix=torch.diag(torch.tensor([1.0, 0.0], dtype=torch.complex128)).reshape(1, 1, 2, 2))
    D = float(rho.relative_entropy(sigma))
    assert abs(D - (0.5 * math.log(1e12) - math.log(2))) < 1e-9


# kills: channels#1745, channels#1744
def test_identity_channel_is_unital_and_trace_preserving_with_zero_tolerance():
    """Channel.is_unital / is_trace_preserving compare error <= tol (inclusive): the identity channel passes tol=0."""
    from torch_vn_algebra import channels as ch
    alg = TypeIAlgebra([2, 3], [2, 3], precision='double')
    E = ch.identity_channel(alg)
    assert float(E.unitality_error()) == 0.0
    assert E.is_unital(tol=0.0) and E.is_trace_preserving(tol=0.0)




# kills: krylov#7097
def test_ground_state_default_max_restarts_is_20():
    """krylov.ground_state stops after max_restarts = 20 restarts by default (tol = 0 never converges): 20 Lanczos
    runs of m matvecs plus one residual matvec each, and a CostWarning."""
    from torch_vn_algebra import krylov, cost
    g = torch.Generator().manual_seed(0)
    A = torch.randn(1, 6, 6, dtype=torch.float64, generator=g)
    H = _DenseH(A + A.mT)
    calls = [0]
    mv = H.matvec

    def counting(V):
        calls[0] += 1
        return mv(V)
    H.matvec = counting
    with pytest.warns(cost.CostWarning):
        krylov.ground_state(H, m=6, tol=0.0, generator=torch.Generator().manual_seed(1))
    assert calls[0] == 20 * (6 + 1)


# kills: functional#5062
def test_abs_gradient_finite_at_zero_matrix():
    """functional.abs_: |X| and its gradient at X = 0 are finite (singular values at the threshold 0 are dropped)."""
    from torch_vn_algebra import functional as fc
    X = torch.zeros(1, 2, 2, dtype=torch.float64, requires_grad=True)
    Y = fc.abs_(X)
    assert float(Y.detach().abs().max()) == 0.0
    Y.sum().backward()
    assert torch.isfinite(X.grad).all()


# kills: channels#1863
def test_identity_channel_kraus_does_not_share_memory_with_the_algebra_cache():
    """identity_channel(alg).kraus is a copy: modifying it in place leaves alg.identity() unchanged."""
    from torch_vn_algebra import channels as ch
    alg = TypeIAlgebra([2, 2], [2, 2], precision='double')
    E = ch.identity_channel(alg)
    E.kraus.mul_(2.0)
    assert torch.equal(alg.identity().matrix[0, 0], torch.eye(2, dtype=alg.dtype))


# kills: algebra#728
def test_eigh_with_spectral_cache_keeps_gradient():
    """With tv.caching() on, eigh() of an operator requiring grad (grad mode on) is differentiable (not from the cache)."""
    alg = TypeIAlgebra([2], [2], precision='double')
    M = torch.tensor([[2.0, 0.5], [0.5, 1.0]], dtype=torch.complex128, requires_grad=True)
    with tv.caching():
        A = alg.operator(M.reshape(1, 1, 2, 2), is_self_adjoint=True)
        w, V = A.eigh()[0]
    assert w.requires_grad
    w.sum().backward()                                   # d(trace)/dM = 1
    assert torch.allclose(M.grad.real, torch.eye(2, dtype=torch.float64), atol=1e-12)


# kills: inclusions#5992
def test_conditional_expectation_rejects_operator_of_wrong_algebra():
    """ConditionalExpectation(y) raises ValueError for an operator that is not in the larger algebra M."""
    from torch_vn_algebra import inclusions as inc
    N = TypeIAlgebra([2], [2], precision='double')
    I = inc.from_tensor_product(N, 2)
    E = I.conditional_expectation()
    assert torch.allclose(E(I.M.identity()).matrix, N.identity().matrix)
    with pytest.raises(ValueError):
        E(N.identity())


# kills: dynamics#4388
def test_solve_operator_ode_default_substeps_is_10():
    """dynamics.solve_operator_ode uses 10 RK4 substeps per interval by default (bitwise equal to substeps=10)."""
    from torch_vn_algebra import dynamics as dy
    alg = TypeIAlgebra([2], [2], precision='double')
    X0 = _herm(alg, seed=9, batch=1)
    f = lambda t, X: -X
    a = dy.solve_operator_ode(f, X0, [0.0, 1.0])
    b = dy.solve_operator_ode(f, X0, [0.0, 1.0], substeps=10)
    assert torch.equal(a[-1].matrix, b[-1].matrix)


# kills: jordan#6273
def test_jordan_triple_self_adjoint_flag_needs_all_three():
    """jordan.jordan_triple(X, Y, Z) is flagged self-adjoint only if X, Y and Z all are: {1, Y, 1} = Y is not."""
    from torch_vn_algebra import jordan as jd
    alg = TypeIAlgebra([2], [2], precision='double')
    Y = alg.operator(torch.tensor([[[[0.0, 1.0], [0.0, 0.0]]]], dtype=torch.complex128))
    I = alg.identity()
    T = jd.jordan_triple(I, Y, I)
    assert torch.allclose(T.matrix, Y.matrix)
    assert T.is_self_adjoint is False


# kills: composite#2228
def test_pair_offset_in_fused_sector_with_repeated_charges():
    """composite.pair_offset(alg12, c, d) finds the pair (c, d) itself (both indices) inside a fused sector."""
    from torch_vn_algebra import composite as cp
    a1 = TypeIAlgebra([1], [1], precision='double', charges=[0])
    a2 = TypeIAlgebra([1, 2], [1, 2], precision='double', charges=[0, 0])
    a12 = tv.fused_tensor_product(a1, a2)
    assert a12.C == 1
    assert cp.pair_offset(a12, 0, 0) == (0, 0)
    assert cp.pair_offset(a12, 0, 1) == (0, 1)


# kills: dynamics#4619
def test_lindblad_superoperator_empty_sector_is_empty():
    """dynamics.lindblad_superoperator gives (batch, 0, 0) generators for sectors with k_c = 0."""
    from torch_vn_algebra import dynamics as dy
    alg = TypeIAlgebra([2, 1], [2, 0], precision='double')
    H = _herm(alg, seed=10, batch=1)
    G = dy.lindblad_superoperator(alg, H, [])
    assert [tuple(g.shape) for g in G] == [(1, 4, 4), (1, 0, 0)]




# kills: functional#4867
def test_normal_eig_clustered_spectrum_reconstructs():
    """functional.normal_eig on normal matrices with three clusters of nearly equal eigenvalues (spread 1e-6, which
    forces the later passes with other mixing coefficients): V diag(w) V^* = X and V unitary."""
    from torch_vn_algebra import functional as fc
    g = torch.Generator().manual_seed(3)
    B, k = 16, 60
    for _ in range(2):
        base = torch.randn(B, 3, dtype=torch.complex128, generator=g)
        idx = torch.randint(0, 3, (B, k), generator=g)
        lam = torch.gather(base, 1, idx) + 1e-6 * torch.randn(B, k, dtype=torch.complex128, generator=g)
        Q, _ = torch.linalg.qr(torch.randn(B, k, k, dtype=torch.complex128, generator=g))
        X = (Q * lam[:, None, :]) @ Q.mH
        w, V = fc.normal_eig(X)
        R = (V * w[:, None, :]) @ V.mH - X
        assert float(R.abs().amax()) < 1e-8 * float(X.abs().amax())
        assert torch.allclose(V.mH @ V, torch.eye(k, dtype=V.dtype).expand(B, k, k), atol=1e-10)


# kills: composite#2247
def test_kron_self_adjoint_flag_needs_both_factors():
    """kron(A, B) is flagged self-adjoint only if both factors are: kron(1, Y) with Y not self-adjoint is not."""
    a1 = TypeIAlgebra([2], [2], precision='double')
    a12 = tv.tensor_product(a1, a1)
    Y = a1.operator(torch.tensor([[[[0.0, 1.0], [0.0, 0.0]]]], dtype=torch.complex128))
    I = a1.identity()
    assert tv.kron(I, I, a12).is_self_adjoint
    assert tv.kron(I, Y, a12).is_self_adjoint is False
    assert tv.kron(Y, I, a12).is_self_adjoint is False


# kills: dynamics#4624
def test_lindblad_channel_with_empty_sector():
    """dynamics.lindblad_channel works for algebras with an empty sector (trace preserving on the others)."""
    from torch_vn_algebra import dynamics as dy
    alg = TypeIAlgebra([2, 1], [2, 0], precision='double')
    H = _herm(alg, seed=11, batch=1)
    L = alg.operator(torch.tensor([[[[0.0, 1.0], [0.0, 0.0]], [[0.0, 0.0], [0.0, 0.0]]]], dtype=torch.complex128))
    E = dy.lindblad_channel(alg, H, [L], 0.5)
    assert float(E.trace_preservation_error().max()) < 1e-10


# kills: divergences#4157
def test_universal_recovery_nodes_follow_beta0_density():
    """divergences.universal_recovery_nodes: weights proportional to beta_0(t) = pi/2 / (cosh(pi t) + 1), sum 1."""
    ts, ws = dv.universal_recovery_nodes(step=0.5, cutoff=3.0)
    assert ts == [0.5 * j for j in range(-6, 7)]
    ref = np.array([0.5 * math.pi / (math.cosh(math.pi * t) + 1) for t in ts])
    ref = ref / ref.sum()
    assert np.allclose(ws, ref, atol=1e-15)
    assert abs(ws[6] / ws[8] - (math.cosh(math.pi) + 1) / 2) < 1e-12


# kills: channels#2143
def test_inter_sector_channel_trace_preserving_with_zero_tolerance():
    """InterSectorChannel.is_trace_preserving(tol=0) holds for an exactly trace-preserving map (inclusive <=)."""
    from torch_vn_algebra import InterSectorChannel
    a_in = TypeIAlgebra([2], [2], precision='double')
    a_out = TypeIAlgebra([2, 1], [2, 1], precision='double')
    E = InterSectorChannel.from_blocks(a_in, a_out, {(0, 0): [torch.eye(2, dtype=torch.complex128)]})
    assert float(E.trace_preservation_error().max()) == 0.0
    assert E.is_trace_preserving(tol=0.0)


# kills: krylov#7011
def test_sparse_hamiltonian_to_dense_checks_memory():
    """SparseSectorHamiltonian.to_dense estimates the dense size first (L = 14, sector N = 7: 3432^2 complex entries,
    ~190 MB): with max_memory_fraction tiny it raises InsufficientMemoryError instead of allocating."""
    from torch_vn_algebra import cost
    chain = tv.SpinChain(14)
    Hs = chain.xxz_sparse(sector=7)
    assert Hs.dim == 3432
    old = cost.get_limits()
    try:
        cost.set_limits(max_memory_fraction=1e-9, warn_memory_fraction=1e-9)
        with pytest.raises(cost.InsufficientMemoryError):
            Hs.to_dense()
    finally:
        cost.set_limits(**old)




def _second_derivative_of_square(A0, E, C):
    """d^2/dt^2 Re sum(F(A0 + t E) o C) at t = 0 for F = spectral_apply(., w -> w^2), via double backward."""
    from torch_vn_algebra import functional as fc
    t = torch.zeros((), dtype=torch.float64, requires_grad=True)
    X = (A0 + t * E).unsqueeze(0)
    F = fc.spectral_apply(X, lambda w: w ** 2)
    L = (F[0] * C).sum().real
    d1, = torch.autograd.grad(L, t, create_graph=True)
    d2, = torch.autograd.grad(d1, t)
    return float(d2)


# kills: functional#5010
def test_spectral_apply_second_derivative_complex_hermitian():
    """functional.spectral_apply: second derivatives (create_graph) for complex Hermitian X are exact at a
    non-degenerate spectrum: for f(w) = w^2, d^2/dt^2 (A + tE)^2 = 2 E^2."""
    A0 = torch.tensor([[1.0, 0.5 - 0.7j, 0.2j], [0.5 + 0.7j, -0.3, 0.4], [-0.2j, 0.4, 2.0]], dtype=torch.complex128)
    E = torch.tensor([[0.3, 0.1 + 0.4j, -0.5], [0.1 - 0.4j, 0.8, 0.2 - 0.1j], [-0.5, 0.2 + 0.1j, -0.6]],
                     dtype=torch.complex128)
    C = torch.tensor([[1.0, 2.0j, -1.0], [0.5, -1.0, 1.0j], [0.3, 0.7, 1.5]], dtype=torch.complex128)
    got = _second_derivative_of_square(A0, E, C)
    ref = float(((2 * E @ E) * C).sum().real)
    assert abs(got - ref) < 1e-9 * max(1.0, abs(ref))


# kills: functional#4884
def test_spectral_apply_warns_for_degenerate_2x2_second_derivatives():
    """functional.spectral_apply: second derivatives at a degenerate 2 x 2 spectrum issue a RuntimeWarning."""
    A0 = torch.eye(2, dtype=torch.complex128)
    E = torch.tensor([[0.0, 1.0], [1.0, 0.0]], dtype=torch.complex128)
    C = torch.ones(2, 2, dtype=torch.complex128)
    with pytest.warns(RuntimeWarning, match="degenerate"):
        _second_derivative_of_square(A0, E, C)


# kills: cost#2633
def test_step_timer_two_steps_no_spurious_warning():
    """cost.StepTimer with 2 steps warns only if the first step took longer than warn_seconds (a fast loop: none)."""
    import warnings
    from torch_vn_algebra import cost
    with warnings.catch_warnings():
        warnings.simplefilter('error', cost.CostWarning)
        t = cost.StepTimer(2, "two fast steps")
        t.step()
        t.step()
        t.close()


# kills: states#7458
def test_density_matrix_entropy_default_eps():
    """DensityMatrix.von_neumann_entropy uses eps = 1e-12 by default: an eigenvalue 1e-10 contributes -p log p."""
    alg = TypeIAlgebra([2], [2], precision='double')
    p = 1e-10
    rho = tv.DensityMatrix(alg, matrix=torch.diag(torch.tensor([1 - p, p], dtype=torch.complex128)).reshape(1, 1, 2, 2))
    ref = -(1 - p) * math.log(1 - p) - p * math.log(p)
    assert abs(float(rho.von_neumann_entropy()) - ref) < 1e-13


# kills: channels#1747
def test_choi_memory_estimate_is_k_squared():
    """Channel.choi estimates B k^4 entries: 100 blocks of size 12 (33 MB) stay below the 64 MB floor of
    cost.check_memory, so it runs even with max_memory_fraction tiny."""
    from torch_vn_algebra import channels as ch, cost
    torch.manual_seed(0)
    alg = TypeIAlgebra([12], [12], precision='double')
    E = ch.random_channel(alg, 1, batch_size=100)
    old = cost.get_limits()
    try:
        cost.set_limits(max_memory_fraction=1e-9, warn_memory_fraction=1e-9)
        J = E.choi()
    finally:
        cost.set_limits(**old)
    assert J[0].shape == (100, 144, 144)




# kills: divergences#4103
def test_petz_recovery_rejects_sigma_of_wrong_algebra():
    """PetzRecovery(sigma, E) raises ValueError when sigma is not an operator of the input algebra of E."""
    from torch_vn_algebra import channels as ch, states as st
    torch.manual_seed(0)
    alg = TypeIAlgebra([2, 3], [2, 3], precision='double')
    other = TypeIAlgebra([3, 2], [3, 2], precision='double')
    E = ch.random_channel(alg, 2)
    dv.PetzRecovery(st.random_density_matrix(alg), E)
    with pytest.raises(ValueError):
        dv.PetzRecovery(st.random_density_matrix(other), E)


# kills: inclusions#5578
def test_from_tensor_product_with_trivial_factor_keeps_ambient_dimensions():
    """inclusions.from_tensor_product(N, 1): N (x) M_1 = N, the new factor has n = k = 1 (ambient sizes unchanged)."""
    from torch_vn_algebra import inclusions as inc
    N = TypeIAlgebra([3, 2], [2, 2], precision='double')
    I = inc.from_tensor_product(N, 1)
    assert I.M.k_factors == N.k_factors
    assert I.M.n_factors == N.n_factors


# kills: divergences#4022
def test_sector_contributions_weight_threshold_relative_to_total_trace():
    """sector_contributions: sigma_c counts as present when q_c > rtol * Tr sigma (the total, not the mean over
    sectors). Here q_2 = 0.75 rtol Tr sigma is above the eigenvalue cut rtol * lambda_max(sigma) = 0.25 rtol (so in
    the support) but at the rounding level of the total: the conditional divergence of sector 2 is +inf."""
    alg = TypeIAlgebra([4, 1], [4, 1], precision='double')
    rt = 5 * torch.finfo(torch.float64).eps                  # default rtol = (sum k_c) eps
    d = 0.75 * rt
    rho = _diag_op(alg, [0.125] * 4, [0.5])
    sigma = _diag_op(alg, [0.25 * (1 - d)] * 4, [d])
    sc = dv.sector_contributions(rho, sigma, kind='umegaki')
    assert torch.isfinite(sc.conditional[0, 0])
    assert math.isinf(float(sc.conditional[0, 1]))


# kills: inclusions#6132
def test_bratteli_draw_labels_only_multiple_edges():
    """BratteliDiagram.draw labels edges only for multiplicities > 1 (show_sizes=False: no other texts)."""
    pytest.importorskip("matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from torch_vn_algebra import inclusions as inc
    D = inc.BratteliDiagram([[1, 1], [2]], [np.array([[1, 1]])])
    fig, ax = plt.subplots()
    D.draw(ax, show_sizes=False)
    assert len(ax.texts) == 0
    D2 = inc.BratteliDiagram([[1], [2]], [np.array([[2]])])
    fig2, ax2 = plt.subplots()
    D2.draw(ax2, show_sizes=False)
    assert [t.get_text() for t in ax2.texts] == ["2"]
    plt.close(fig)
    plt.close(fig2)


# kills: channels#2125
def test_inter_sector_compose_with_non_trace_preserving_map_is_not_flagged_trace_preserving():
    """InterSectorChannel.compose: the composite of a trace-preserving and a non-trace-preserving map is not flagged
    trace preserving (applied to a state it gives an unnormalised Operator with a warning)."""
    from torch_vn_algebra import InterSectorChannel, channels as ch, states as st
    torch.manual_seed(0)
    a = TypeIAlgebra([2], [2], precision='double')
    E = ch.random_inter_sector_channel(a, a)                                  # flagged trace preserving
    F = InterSectorChannel.from_blocks(a, a, {(0, 0): [0.5 * torch.eye(2, dtype=torch.complex128)]})
    G = E.compose(F)
    rho = st.random_density_matrix(a)
    with pytest.warns(UserWarning, match="not trace preserving"):
        out = G.apply(rho)
    assert not isinstance(out, tv.DensityMatrix)
    assert abs(float(out.trace.real) - 0.25) < 1e-12


# kills: channels#2183
def test_random_inter_sector_channel_square_isometry_allowed():
    """random_inter_sector_channel accepts n_out * kraus_rank == k_c (a unitary Stinespring map) and is trace
    preserving; n_out * kraus_rank < k_c raises ValueError."""
    from torch_vn_algebra import channels as ch
    torch.manual_seed(0)
    a_in = TypeIAlgebra([4], [4], precision='double')
    a_out = TypeIAlgebra([2], [2], precision='double')
    E = ch.random_inter_sector_channel(a_in, a_out, kraus_rank=2)
    assert float(E.trace_preservation_error().max()) < 1e-12
    with pytest.raises(ValueError):
        ch.random_inter_sector_channel(a_in, a_out, kraus_rank=1)




# kills: cost#2413, cost#2412
def test_cost_default_limits():
    """cost.get_limits(): the documented defaults (warn at 10 s, warn above 50 % / raise above 90 % of free memory)."""
    from torch_vn_algebra import cost
    lim = cost.get_limits()
    assert lim['enabled'] is True
    assert lim['warn_seconds'] == 10.0
    assert lim['warn_memory_fraction'] == 0.5
    assert lim['max_memory_fraction'] == 0.9




# kills: functional#4670
def test_values_and_derivative_keeps_complex_values_of_real_argument():
    """functional.values_and_derivative: a complex-valued f of a real w keeps complex values (and f' = i e^{iw})."""
    from torch_vn_algebra import functional as fc
    w = torch.tensor([0.0, 0.5, 1.0], dtype=torch.float64)
    v, d = fc.values_and_derivative(lambda x: torch.exp(1j * x), w)
    assert torch.is_complex(v)
    assert torch.allclose(v, torch.exp(1j * w.to(torch.complex128)))
    assert torch.allclose(d, 1j * torch.exp(1j * w.to(torch.complex128)))


# kills: functional#4651, functional#4647
def test_caching_context_restores_previous_state():
    """functional.caching(): the spectral cache is switched on inside the block and restored afterwards."""
    from torch_vn_algebra import functional as fc
    assert not fc.cache_enabled()
    with tv.caching():
        assert fc.cache_enabled()
        with tv.caching(False):
            assert not fc.cache_enabled()
        assert fc.cache_enabled()
    assert not fc.cache_enabled()


# kills: cost#2511
def test_batched_call_does_not_probe_small_batches():
    """cost.batched_call times nothing for fewer than 4 * probe_items items (no CostWarning even with
    warn_seconds = 0 and probe_flops = 0)."""
    import warnings
    from torch_vn_algebra import cost
    x = torch.eye(4, dtype=torch.float64).expand(8, 4, 4).contiguous()
    old = cost.get_limits()
    try:
        cost.set_limits(probe_flops=0.0, warn_seconds=0.0)
        with warnings.catch_warnings():
            warnings.simplefilter('error', cost.CostWarning)
            w = cost.batched_call(lambda a: torch.linalg.eigvalsh(a), x, "eigvalsh small batch")
    finally:
        cost.set_limits(**old)
    assert torch.allclose(w, torch.ones(8, 4, dtype=torch.float64))


# kills: chains#1262
def test_spin_chain_checks_memory_of_labels():
    """SpinChain(L) estimates its 2^L int64 label arrays first (L = 22: ~170 MB): with max_memory_fraction tiny it
    raises InsufficientMemoryError."""
    from torch_vn_algebra import cost
    old = cost.get_limits()
    try:
        cost.set_limits(max_memory_fraction=1e-9, warn_memory_fraction=1e-9)
        with pytest.raises(cost.InsufficientMemoryError):
            tv.SpinChain(22)
    finally:
        cost.set_limits(**old)




# kills: divergences#4143
def test_petz_recovery_default_is_untwirled_t0():
    """divergences.petz_recovery(sigma, E, X) is the Petz map (t = 0 by default), not a rotated one."""
    from torch_vn_algebra import channels as ch, states as st
    torch.manual_seed(1)
    alg = TypeIAlgebra([2, 3], [2, 3], precision='double')
    sigma = st.random_density_matrix(alg)
    E = ch.random_channel(alg, 2)
    X = _herm(alg, seed=12, batch=1)
    R0 = dv.PetzRecovery(sigma, E, t=0.0)(X).matrix
    R1 = dv.PetzRecovery(sigma, E, t=1.0)(X).matrix
    assert not torch.allclose(R0, R1, atol=1e-6)
    assert torch.allclose(dv.petz_recovery(sigma, E, X).matrix, R0, atol=1e-12)




# kills: algebra#866
def test_power_nonnegative_exponent_uses_all_eigenvalues():
    """Operator.power(p) for p >= 0 is w^p on every eigenvalue (no support cut): A^0 = 1 and
    diag(1, 1e-20)^(1/2) = diag(1, 1e-10)."""
    alg = TypeIAlgebra([2], [2], precision='double')
    A = alg.operator(torch.diag(torch.tensor([1.0, 0.0], dtype=torch.complex128)).reshape(1, 1, 2, 2),
                     is_self_adjoint=True, is_positive=True)
    assert torch.allclose(A.power(0.0).matrix, alg.identity().matrix)
    B = alg.operator(torch.diag(torch.tensor([1.0, 1e-20], dtype=torch.complex128)).reshape(1, 1, 2, 2),
                     is_self_adjoint=True, is_positive=True)
    S = B.power(0.5).matrix[0, 0]
    assert abs(float(S[1, 1].real) - 1e-10) < 1e-16


# kills: krylov#7132
def test_expm_multiply_step_with_one_krylov_vector():
    """krylov.expm_multiply_step with m = 1: exp(-i dt <H>) psi, exact for an eigenvector (error estimate 0)."""
    from torch_vn_algebra import krylov
    M = torch.diag(torch.tensor([1.0, 2.0], dtype=torch.float64)).to(torch.complex128)[None]
    psi = torch.tensor([[0.0, 1.0]], dtype=torch.complex128)
    out, err = krylov.expm_multiply_step(_DenseH(M), psi, 0.25, m=1)
    assert torch.allclose(out, torch.exp(torch.tensor(-0.5j)) * psi, atol=1e-14)
    assert float(err.max()) < 1e-14




# kills: cost#2423
def test_cost_disabled_context_restores_checks():
    """cost.disabled(): the checks are off inside the block and on again afterwards."""
    from torch_vn_algebra import cost
    assert cost.get_limits()['enabled']
    with cost.disabled():
        assert not cost.get_limits()['enabled']
    assert cost.get_limits()['enabled']


# kills: decompose#2880
def test_commutant_twirl_rejects_non_group():
    """commutant_basis(..., group=True, method='twirl') raises ValueError for unitaries that are not a whole group
    ({1, diag(1, i)} is not closed: the twirl has the wrong rank)."""
    S = torch.diag(torch.tensor([1.0, 1j], dtype=torch.complex128))
    Z = torch.diag(torch.tensor([1.0, -1.0], dtype=torch.complex128))
    G = torch.stack([torch.eye(2, dtype=torch.complex128), Z])               # a group: Z_2
    B = tv.commutant_basis(G, group=True, method='twirl')
    assert B.shape[0] == 2
    with pytest.raises(ValueError):
        tv.commutant_basis(torch.stack([torch.eye(2, dtype=torch.complex128), S]), group=True, method='twirl')


# kills: functional#4780
def test_normal_eig_with_ties_of_every_mixing_coefficient():
    """functional.normal_eig separates eigenvalues 0 and a - i that tie in Re + a Im for each of the mixing
    coefficients a (white-box: functional._ALPHAS), by refining inside the clusters."""
    from torch_vn_algebra import functional as fc
    lam = torch.tensor([0.0] + [a - 1j for a in fc._ALPHAS], dtype=torch.complex128)
    k = lam.numel()
    g = torch.Generator().manual_seed(5)
    Q, _ = torch.linalg.qr(torch.randn(k, k, dtype=torch.complex128, generator=g))
    X = ((Q * lam) @ Q.mH)[None]
    w, V = fc.normal_eig(X)
    R = (V * w[:, None, :]) @ V.mH - X
    assert float(R.abs().amax()) < 1e-10
    assert torch.allclose(torch.sort(w[0].real).values, torch.sort(lam.real).values, atol=1e-10)
