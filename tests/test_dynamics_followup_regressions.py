"""Physical and analytic references for the dynamics follow-up review."""
import math

import pytest
import torch

from torch_vn_algebra import TypeIAlgebra, dynamics, krylov
from torch_vn_algebra.permutational import PermutationalLindbladian, PermutationalSystem


@pytest.mark.parametrize('dtype', [torch.float32, torch.float64, torch.complex64, torch.complex128])
@pytest.mark.parametrize('magnitude', ['subnormal', 'small-normal', 'large'])
def test_scaled_norm_retains_extreme_values_and_their_gradients(dtype, magnitude):
    real = torch.float32 if dtype in (torch.float32, torch.complex64) else torch.float64
    limits = torch.finfo(real)
    small = {'subnormal': limits.tiny / 16,
             'small-normal': limits.tiny * 2,
             'large': limits.max / 16}[magnitude]
    values = torch.tensor([[small, 2 * small]], dtype=dtype, requires_grad=True)
    if values.detach().real[0, 0].item() == 0:
        pytest.skip('this backend flushes subnormal input values to zero')
    expected = torch.tensor(math.hypot(small, 2 * small), dtype=real)
    actual = krylov._norm(values)
    # Compare relative values: an absolute tolerance would hide this underflow regression.
    torch.testing.assert_close(actual / expected, torch.ones_like(actual))
    derivative, = torch.autograd.grad(actual.sum(), values)
    torch.testing.assert_close(derivative, torch.tensor([[1 / math.sqrt(5), 2 / math.sqrt(5)]], dtype=dtype))


def _pauli_hamiltonian(shift=0.0, scale=1.0, precision='double', batch=1):
    real = torch.float32 if precision == 'single' else torch.float64
    dtype = torch.complex64 if precision == 'single' else torch.complex128
    off = torch.tensor([[0.0, scale], [scale, 0.0]], dtype=real).to_sparse()
    return krylov.SparseSectorHamiltonian(torch.full((batch, 2), shift, dtype=real), off, dtype)


@pytest.mark.parametrize('precision,shift', [('single', 0.0), ('single', 2e5), ('single', 2e12),
                                            ('double', 0.0), ('double', 1e14), ('double', 1e16)])
@pytest.mark.parametrize('scale', [0.125, 1.0, 8.0])
@pytest.mark.parametrize('initial', ['basis', 'superposition'])
def test_lanczos_keeps_rabi_oscillations_under_scalar_shifts(precision, shift, scale, initial):
    H = _pauli_hamiltonian(shift, scale, precision)
    psi = torch.tensor([[1.0, 0.0]] if initial == 'basis' else [[math.sqrt(.3), 1j * math.sqrt(.7)]],
                       dtype=H.dtype)
    out, error = krylov.expm_multiply_step(H, psi, math.pi / (2 * scale), m=2)
    # exp[-i(c I + scale X) pi/(2 scale)] exchanges the two populations, for every c.
    expected = psi.abs().square().flip(-1)
    tol = 2e-5 if precision == 'single' else 2e-12
    torch.testing.assert_close(out.abs().square(), expected, atol=tol, rtol=tol)
    assert torch.isfinite(error).all()


def test_public_lanczos_diagonal_still_represents_original_hamiltonian():
    H = _pauli_hamiltonian(1e14)
    alpha, beta, Q = krylov.lanczos(H.matvec, torch.tensor([[1., 0.]], dtype=H.dtype), 2)
    torch.testing.assert_close(alpha, torch.full_like(alpha, 1e14))
    torch.testing.assert_close(beta, torch.ones_like(beta))
    torch.testing.assert_close(Q @ Q.mH, torch.eye(2, dtype=Q.dtype)[None])


def test_generic_lanczos_matvec_also_keeps_large_shift_coupling():
    dense = torch.tensor([[1e14, 1.0], [1.0, 1e14]], dtype=torch.float64)
    alpha, beta, Q = krylov.lanczos(lambda v: v @ dense.T, torch.tensor([[1., 0.]], dtype=dense.dtype), 2)
    torch.testing.assert_close(beta, torch.ones_like(beta))
    torch.testing.assert_close(Q @ Q.mT, torch.eye(2, dtype=dense.dtype)[None])
    torch.testing.assert_close(alpha, torch.full_like(alpha, 1e14))


@pytest.mark.parametrize('real,dtype,shift', [(torch.float32, torch.complex64, 1e8),
                                          (torch.float64, torch.complex128, 1e16)])
def test_sparse_diagonal_entries_are_centered_before_their_matvec(real, dtype, shift):
    A = torch.tensor([[shift, 1.], [1., shift]], dtype=real).to_sparse()
    H = krylov.SparseSectorHamiltonian(torch.zeros(1, 2, dtype=real), A, dtype)
    psi = torch.tensor([[1., 0.]], dtype=dtype)
    out, error = krylov.expm_multiply_step(H, psi, math.pi / 2, m=2)
    expected = torch.tensor([[0., 1.]], dtype=real)
    tol = 2e-5 if real == torch.float32 else 2e-12
    torch.testing.assert_close(out.abs().square(), expected, atol=tol, rtol=tol)
    assert error.max() < tol


def test_sparse_centering_pivot_is_selected_per_batch_element():
    real, dtype = torch.float32, torch.complex64
    # Two disjoint two-level systems, with different occupied blocks in the two batch elements.
    diagonal = torch.tensor([[0., 0., 1e8, 1e8], [1e8, 1e8, 0., 0.]], dtype=real)
    A = torch.tensor([[0., 1., 0., 0.], [1., 0., 0., 0.],
                      [0., 0., 0., 1.], [0., 0., 1., 0.]], dtype=real).to_sparse()
    H = krylov.SparseSectorHamiltonian(diagonal, A, dtype)
    psi = torch.tensor([[1., 0., 0., 0.], [0., 0., 1., 0.]], dtype=dtype)
    out, error = krylov.expm_multiply_step(H, psi, math.pi / 2, m=2)
    torch.testing.assert_close(out.abs().square(), psi.abs().square().roll(1, dims=-1), atol=2e-5, rtol=2e-5)
    assert error.max() < 2e-5


@pytest.mark.parametrize('precision,spectator', [('single', 1e8), ('double', 1e16)])
@pytest.mark.parametrize('spectator_first', [False, True])
@pytest.mark.parametrize('diagonal_storage', ['diagonal', 'sparse'])
def test_lanczos_rabi_coupling_is_not_scaled_by_an_unvisited_spectator(
        precision, spectator, spectator_first, diagonal_storage):
    real = torch.float32 if precision == 'single' else torch.float64
    dtype = torch.complex64 if precision == 'single' else torch.complex128
    diagonal = torch.tensor([[0., 0., spectator]], dtype=real)
    off = torch.tensor([[0., 1., 0.], [1., 0., 0.], [0., 0., 0.]], dtype=real)
    psi = torch.tensor([[1., 0., 0.]], dtype=dtype)
    expected = torch.tensor([[0., 1., 0.]], dtype=real)
    if spectator_first:
        order = torch.tensor([2, 0, 1])
        diagonal, off = diagonal[:, order], off[order][:, order]
        psi, expected = psi[:, order], expected[:, order]
    if diagonal_storage == 'sparse':
        off = off + torch.diag(diagonal[0])
        diagonal = torch.zeros_like(diagonal)
    H = krylov.SparseSectorHamiltonian(diagonal, off.to_sparse(), dtype)
    out, error = krylov.expm_multiply_step(H, psi, math.pi / 2, m=2)
    tol = 2e-5 if precision == 'single' else 2e-12
    torch.testing.assert_close(out.abs().square(), expected, atol=tol, rtol=tol)
    assert error.max() < tol
    adaptive = krylov.evolve(H, psi, [0., math.pi / 2], m=2, tol=tol)
    torch.testing.assert_close(adaptive[-1].abs().square(), expected, atol=tol, rtol=tol)


def test_early_lanczos_breakdown_does_not_report_zero_truncation_error():
    weak = 1e-7
    off = torch.tensor([[0., 1., 0.], [1., 0., weak], [0., weak, 0.]], dtype=torch.float32).to_sparse()
    H = krylov.SparseSectorHamiltonian(torch.zeros(1, 3), off, torch.complex64)
    psi = torch.tensor([[1., 0., 0.]], dtype=H.dtype)
    out, error = krylov.expm_multiply_step(H, psi, 1., m=3)
    # The weak discarded coupling produces this exact third-coordinate amplitude in a three-site chain.
    omega = math.sqrt(1 + weak ** 2)
    exact_third = weak / (1 + weak ** 2) * (math.cos(omega) - 1)
    assert error[0] > 0
    assert abs(out[0, 2].item() - exact_third) <= error[0].item()


@pytest.mark.parametrize('origin', [0.0, 1e9, 1e12])
@pytest.mark.parametrize('interval', [1.0, 1e-4])
def test_adaptive_krylov_elapsed_time_is_independent_of_grid_origin(origin, interval):
    end = origin + interval
    elapsed = end - origin
    H = _pauli_hamiltonian(scale=1 / elapsed)
    psi = torch.tensor([[1., 0.]], dtype=H.dtype)
    out = krylov.evolve(H, psi, [origin, end], m=2)
    expected = torch.tensor([[math.cos(1.) ** 2, math.sin(1.) ** 2]], dtype=torch.float64)
    torch.testing.assert_close(out[-1].abs().square(), expected, atol=2e-12, rtol=2e-12)


@pytest.mark.parametrize('timedep', [False, True])
def test_schrodinger_rk4_broadcasts_initial_state_over_hamiltonian_batch(timedep):
    alg = TypeIAlgebra([2], [2], precision='double', device='cpu')
    x = torch.tensor([[0., 1.], [1., 0.]], dtype=alg.dtype)
    frequencies = torch.tensor([1., 2.], dtype=torch.float64)

    def H(t):
        return alg.operator((frequencies + (t if timedep else 0))[:, None, None, None] * x,
                            is_self_adjoint=True)

    psi = torch.tensor([[[1., 0.]]], dtype=alg.dtype)
    out = dynamics.schrodinger_rk4(psi, H if timedep else H(0), [0., .1], substeps=64)
    angles = .1 * frequencies + (.005 if timedep else 0)
    expected = torch.stack([angles.cos().to(alg.dtype), -1j * angles.sin()], dim=-1)[:, None]
    assert out.shape == (2, 2, 1, 2)
    torch.testing.assert_close(out[0], psi.expand(2, 1, 2))
    torch.testing.assert_close(out[-1], expected, atol=2e-12, rtol=2e-12)


@pytest.mark.parametrize('N', [1, 2, 5])
def test_coherent_state_pole_gradients_match_spin_moments(N):
    system = PermutationalSystem(N, device='cpu')
    theta = torch.tensor([0., .2], dtype=torch.float64, requires_grad=True)
    phi = .37
    rho = system.coherent_state(theta, phi)
    x = rho.expectation(system.jx()).real.sum()
    y = rho.expectation(system.jy()).real.sum()
    dx, = torch.autograd.grad(x, theta, retain_graph=True)
    dy, = torch.autograd.grad(y, theta)
    torch.testing.assert_close(dx, N / 2 * theta.detach().cos() * math.cos(phi), atol=2e-12, rtol=2e-12)
    torch.testing.assert_close(dy, N / 2 * theta.detach().cos() * math.sin(phi), atol=2e-12, rtol=2e-12)


@pytest.mark.parametrize('kind', ['diagonal', 'offdiagonal'])
def test_sparse_matvec_refreshes_casts_after_input_mutation(kind):
    diagonal = torch.tensor([[1., 2.]], dtype=torch.float32)
    off = torch.tensor([[.2, .3], [.3, -.1]], dtype=torch.float64).to_sparse().coalesce()
    H = krylov.SparseSectorHamiltonian(diagonal, off, torch.complex128)
    v = torch.tensor([[1.+.2j, -.3j]], dtype=H.dtype)
    shift = torch.tensor([.7], dtype=torch.float64)
    H.matvec(v)                     # populate both dtype-conversion caches
    H.matvec_shifted(v, shift)      # also populate the split sparse-diagonal caches
    if kind == 'diagonal':
        diagonal.add_(1)
    else:
        H.offdiag.values().mul_(2)
    expected = (H.to_dense() @ v.unsqueeze(-1)).squeeze(-1)
    torch.testing.assert_close(H.matvec(v), expected)
    torch.testing.assert_close(H.matvec_shifted(v, shift), expected - shift[:, None] * v)


def test_sparse_diagonal_casts_do_not_reuse_a_finished_autograd_graph():
    diagonal = torch.tensor([[1., 2.]], dtype=torch.float64, requires_grad=True)
    H = krylov.SparseSectorHamiltonian(diagonal, None, torch.complex128)
    v = torch.tensor([[1.+1j, 2.-1j]], dtype=H.dtype)
    for _ in range(2):
        loss = H.matvec(v).abs().square().sum()
        grad, = torch.autograd.grad(loss, diagonal)
        torch.testing.assert_close(grad, 2 * diagonal.detach() * v.abs().square())


def test_sparse_diagonal_split_preserves_repeated_gradients():
    values = torch.tensor([.2, .3, .3, -.1], dtype=torch.float64)
    indices = torch.tensor([[0, 0, 1, 1], [0, 1, 0, 1]])
    off = torch.sparse_coo_tensor(indices, values, (2, 2)).coalesce().requires_grad_(True)
    H = krylov.SparseSectorHamiltonian(torch.tensor([[1., 2.]], dtype=torch.float64), off, torch.float64)
    v = torch.tensor([[1.2, -.3]], dtype=H.dtype)
    shift = torch.tensor([.7], dtype=torch.float64)
    for _ in range(2):
        loss = H.matvec_shifted(v, shift).abs().square().sum()
        grad, = torch.autograd.grad(loss, off)
        dense_values = values.detach().clone().requires_grad_(True)
        dense = torch.diag(H.diag[0] - shift) + dense_values.reshape(2, 2)
        reference, = torch.autograd.grad((v @ dense.to(v.dtype).T).abs().square().sum(), dense_values)
        torch.testing.assert_close(grad.to_dense(), reference.reshape(2, 2))


def test_sparse_casts_refresh_unversioned_inference_tensors():
    with torch.inference_mode():
        diagonal = torch.tensor([[1., 2.]], dtype=torch.float32)
        H = krylov.SparseSectorHamiltonian(diagonal, None, torch.complex128)
        v = torch.ones(1, 2, dtype=H.dtype)
        H.matvec(v)
        diagonal.add_(1)
        torch.testing.assert_close(H.matvec(v), torch.tensor([[2., 3.]], dtype=H.dtype))


def test_permutational_rhs_uses_the_current_rate_tensor():
    system = PermutationalSystem(1, device='cpu')
    rate = torch.tensor([.1], dtype=torch.float64)
    generator = PermutationalLindbladian(system, emission=rate)
    rho = system.coherent_state(0).matrix
    first = generator.rhs(0., rho)
    torch.testing.assert_close(first[0, 0].diagonal().real, torch.tensor([-.1, .1], dtype=torch.float64))
    rate.fill_(.5)
    after = generator.rhs(0., rho)
    torch.testing.assert_close(after[0, 0].diagonal().real, torch.tensor([-.5, .5], dtype=torch.float64))
    generator.params['emission'] = torch.tensor([.2], dtype=torch.float64)
    replaced = generator.rhs(0., rho)
    torch.testing.assert_close(replaced[0, 0].diagonal().real, torch.tensor([-.2, .2], dtype=torch.float64))


def test_permutational_rhs_refreshes_unversioned_inference_rates():
    system = PermutationalSystem(1, device='cpu')
    rho = system.coherent_state(0).matrix
    with torch.inference_mode():
        rate = torch.tensor([.1], dtype=torch.float64)
        generator = PermutationalLindbladian(system, emission=rate)
        generator.rhs(0., rho)
        rate.fill_(.5)
        actual = generator.rhs(0., rho)[0, 0].diagonal().real
        torch.testing.assert_close(actual, torch.tensor([-.5, .5], dtype=torch.float64))


@pytest.mark.parametrize('initial_rate', [0., .1])
def test_permutational_rate_derivatives_survive_repeated_forward_backward(initial_rate):
    system = PermutationalSystem(1, device='cpu')
    rate = torch.tensor([initial_rate], dtype=torch.float64, requires_grad=True)
    generator = PermutationalLindbladian(system, emission=rate)
    rho = system.coherent_state(0).matrix
    for _ in range(2):
        excitation_rate = generator.rhs(0., rho)[0, 0, 1, 1].real
        grad, = torch.autograd.grad(excitation_rate, rate)
        torch.testing.assert_close(grad, torch.ones_like(rate))
