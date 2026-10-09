"""Tests for fused (charge-conserving) tensor products and spin chains."""
import itertools
import math

import pytest
import torch

from torch_vn_algebra import (TypeIAlgebra, SpinChain, DensityMatrix, fused_tensor_product, tensor_product,
                              kron, partial_trace, dynamics, states as st)
from torch_vn_algebra import composite as cp


def site():
    return TypeIAlgebra([1, 1], [1, 1], precision='double', charges=[0, 1], device='cpu')


def hopping_hamiltonian(L):
    h = torch.zeros(L, L, dtype=torch.complex128)
    for i in range(L - 1):
        h[i, i + 1] = h[i + 1, i] = 0.5
    return h


class TestFusion:
    def test_dimensions_and_charges(self):
        a = TypeIAlgebra([1, 2, 1], [1, 2, 1], precision='double', charges=[0, 1, 2], device='cpu')
        f = fused_tensor_product(a, a)
        assert f.charges == [0, 1, 2, 3, 4]
        assert f.k_factors == [1, 4, 6, 4, 1]          # binomial: two copies of two sites
        assert tensor_product(a, a).C == 9

    def test_kron_and_partial_trace_in_fused_algebra(self):
        s = site()
        f = fused_tensor_product(s, s)
        r1 = DensityMatrix(s, matrix=torch.tensor([0.3, 0.7], dtype=torch.float64).reshape(1, 2, 1, 1).to(torch.complex128))
        r2 = DensityMatrix(s, matrix=torch.tensor([0.6, 0.4], dtype=torch.float64).reshape(1, 2, 1, 1).to(torch.complex128))
        rho = kron(r1, r2, f)
        assert isinstance(rho, DensityMatrix)
        assert torch.allclose(partial_trace(rho, 1).matrix, r1.matrix)
        assert torch.allclose(partial_trace(rho, 2).matrix, r2.matrix)

    def test_entangled_state_needs_the_fused_algebra(self):
        # (|10> + |01>)/sqrt 2 lives in the merged sector N = 1; its marginals are maximally mixed
        s = site()
        f = fused_tensor_product(s, s)
        psi = torch.tensor([1.0, 1.0], dtype=torch.complex128) / 2 ** 0.5
        mat = torch.zeros(1, 3, 2, 2, dtype=torch.complex128)
        mat[0, 1] = torch.outer(psi, psi)
        rho = DensityMatrix(f, matrix=mat)
        for keep in (1, 2):
            m = partial_trace(rho, keep).matrix
            assert torch.allclose(m[0, :, 0, 0].real, torch.tensor([0.5, 0.5], dtype=torch.float64))
        assert cp.pair_offset(f, 0, 1) == (1, 0) and cp.pair_offset(f, 1, 0) == (1, 1)

    def test_chain_basis_is_iterated_fusion(self):
        c3, c4 = SpinChain(3), SpinChain(4)
        s = site()
        f = fused_tensor_product(c3.algebra, s)
        assert f.k_factors == c4.algebra.k_factors
        nsite = s.from_blocks([torch.full((1, 1), -0.5), torch.full((1, 1), 0.5)])
        for i in range(3):
            assert torch.allclose(kron(c3.sz(i), s.identity(1), f).matrix, c4.sz(i).matrix)
        assert torch.allclose(kron(c3.algebra.identity(1), nsite, f).matrix, c4.sz(3).matrix)


class TestSpinChain:
    def test_xx_chain_is_free_fermions(self):
        L = 7
        ch = SpinChain(L)
        eps = [math.cos(math.pi * k / (L + 1)) for k in range(1, L + 1)]
        for N, (w, _) in enumerate(ch.xxz(J=1.0, Delta=0.0).eigh()):
            exact = sorted(sum(c) for c in itertools.combinations(eps, N)) if N else [0.0]
            assert torch.allclose(torch.sort(w[0])[0], torch.tensor(exact, dtype=torch.float64), atol=1e-10)

    def test_heisenberg_ring(self):
        w = [w for w, _ in SpinChain(4, 'periodic').heisenberg().eigh()]
        assert abs(min(x.min().item() for x in w) + 2.0) < 1e-12
        # SU(2): every level of sector N reappears in sector N + 1 (N < L / 2)
        E = [torch.sort(x[0])[0] for x, _ in SpinChain(8, 'periodic').heisenberg().eigh()]
        for N in range(4):
            assert all((E[N + 1] - e).abs().min() < 1e-9 for e in E[N])

    def test_sector_operator_matches_full(self):
        ch = SpinChain(6, 'periodic')
        h = torch.randn(3, 6, dtype=torch.float64)
        full = ch.xxz(J=1.0, Delta=0.7, h=h)
        sec = ch.xxz(J=1.0, Delta=0.7, h=h, sector=3)
        assert torch.allclose(full.matrix[:, 3, :20, :20], sec.matrix[:, 0])
        assert torch.allclose(full.matrix[:, 3], full.matrix[:, 3].conj().transpose(-2, -1))

    def test_number_operator_is_central(self):
        ch = SpinChain(5)
        Nop = ch.number().matrix
        for N in range(6):
            d = ch.sector_dim(N)
            assert torch.allclose(Nop[0, N, :d, :d], N * torch.eye(d, dtype=Nop.dtype))

    def test_domain_wall_dynamics_vs_free_fermions(self):
        L, N = 8, 4
        ch = SpinChain(L)
        H = ch.xxz(J=1.0, Delta=0.0, sector=N)
        psi = ch.vector_in_sector('1' * N + '0' * (L - N))
        ts = [0.0, 1.3, 3.0]
        out = dynamics.schrodinger(psi[None, None], H, ts)
        bits = ch.bits(N).double()
        C0 = torch.diag(torch.tensor([1.0] * N + [0.0] * (L - N), dtype=torch.complex128))
        for t, p in zip(ts, out):
            U = torch.linalg.matrix_exp(-1j * t * hopping_hamiltonian(L))
            n_exact = (U.conj() @ C0 @ U.T).diagonal().real
            assert torch.allclose((p[0, 0].abs() ** 2) @ bits, n_exact, atol=1e-10)

    def test_entanglement_matches_peschel(self):
        L, N = 8, 4
        ch = SpinChain(L)
        w, V = ch.xxz(J=1.0, Delta=0.0, sector=N).eigh()[0]
        psi = V[0, :, 0]
        rho = DensityMatrix(ch.sector_algebra(N), matrix=torch.outer(psi, psi.conj())[None, None])
        e, U = torch.linalg.eigh(hopping_hamiltonian(L).real)
        Cm = U[:, :N] @ U[:, :N].T
        for nA in (1, 2, 4):
            nu = torch.linalg.eigvalsh(Cm[:nA, :nA]).clamp(1e-14, 1 - 1e-14)
            S = -(nu * torch.log(nu) + (1 - nu) * torch.log(1 - nu)).sum()
            assert abs(ch.entanglement_entropy(rho, nA, sector=N).item() - S.item()) < 1e-10

    def test_reduced_state_of_full_algebra_state(self):
        ch = SpinChain(5)
        rho = st.random_density_matrix(ch.algebra, 3)
        red = ch.reduced_state(rho, 2)
        assert isinstance(red, DensityMatrix)
        assert torch.allclose(red.trace.real, torch.ones(3, dtype=torch.float64))
        assert red.lambda_min.min() > -1e-12

    def test_site_amplitude_damping(self):
        ch = SpinChain(4)
        Phi = ch.site_amplitude_damping(2, 0.3)
        assert Phi.is_trace_preserving()
        T = Phi.transition_matrix()[0]
        # from sector N the flip happens if site 2 is up: probability 0.3 * N / 4
        for N in range(1, 5):
            assert abs(T[N - 1, N].item() - 0.3 * N / 4) < 1e-12
        out = Phi(ch.basis_state('0110'))
        p = out.sector_probabilities()[0]
        assert torch.allclose(p, torch.tensor([0, 0.3, 0.7, 0, 0], dtype=p.dtype))

    def test_level_statistics_limits(self):
        ch = SpinChain(10, 'periodic')
        g = torch.Generator().manual_seed(0)
        r = {}
        for W in (0.5, 10.0):
            w, _ = ch.random_field_heisenberg(W, 20, sector=5, generator=g).eigh()[0]
            k = w.shape[-1]
            r[W] = SpinChain.level_spacing_ratio(w[:, k // 4: 3 * k // 4]).mean().item()
        assert 0.49 < r[0.5] < 0.56          # GOE 0.5307
        assert 0.36 < r[10.0] < 0.41         # Poisson 0.3863


class TestMomentum:
    def test_momentum_blocks_reproduce_sector_spectrum(self):
        for L, N, D, h in [(8, 4, 1.0, 0.0), (9, 4, 0.5, 0.2), (10, 3, 1.3, 0.0)]:
            ch = SpinChain(L, 'periodic')
            Hm = ch.xxz_momentum(1.0, D, h, sector=N)
            assert sum(ch.momentum_dims(N)) == math.comb(L, N)
            ws = torch.sort(torch.cat([w[0] for w, _ in Hm.eigh()]))[0]
            wd = torch.linalg.eigvalsh(ch.xxz(1.0, D, h, sector=N).matrix[0, 0])
            assert torch.allclose(ws, wd, atol=1e-10)

    def test_ground_state_momentum_marshall(self):
        # Heisenberg ring: ground state at k = 0 for L/2 even and k = pi for L/2 odd
        for L in (8, 10, 12):
            ch = SpinChain(L, 'periodic')
            Hm = ch.xxz_momentum(sector=L // 2)
            e0 = [w[0].min().item() for w, _ in Hm.eigh()]
            m = Hm.algebra.charges[min(range(len(e0)), key=lambda i: e0[i])]
            assert m == (0 if (L // 2) % 2 == 0 else L // 2)

    def test_momentum_requires_periodic_complex(self):
        with pytest.raises(ValueError):
            SpinChain(6).xxz_momentum(sector=3)
        with pytest.raises(ValueError):
            SpinChain(6, 'periodic', complex_valued=False).xxz_momentum(sector=3)


def test_lindblad_loss_between_sectors_is_a_death_process():
    ch = SpinChain(5, 'periodic')
    rho0 = ch.basis_state('11011')                    # N = 4
    g = 0.4
    out = dynamics.lindblad_evolve(rho0, ch.xxz(1.0, 0.7), [ch.lowering(i) for i in range(5)],
                                   [0.0, 0.8, 2.0], rates=[g] * 5, substeps=60)
    for t, r in zip([0.0, 0.8, 2.0], out):
        q = math.exp(-g * t)
        binom = torch.tensor([math.comb(4, n) * q ** n * (1 - q) ** (4 - n) for n in range(5)] + [0.0],
                             dtype=torch.float64)
        assert torch.allclose(r.sector_probabilities()[0], binom, atol=1e-8)
        assert abs(r.trace.real.item() - 1) < 1e-10


def test_post_selection_on_hamming_weight_undoes_uniform_t1_noise():
    # N-conserving H plus uniform loss: the no-jump evolution is H - i g N / 2 with N central,
    # so post-selection on the initial Hamming weight returns the ideal state with p = exp(-w g T)
    n, w, T, g = 5, 2, 1.5, 0.05
    q = SpinChain(n, 'periodic')
    H = q.xxz(1.0, 0.0)
    rho0 = q.basis_state('11000')
    ideal = dynamics.von_neumann(rho0, H, [T])[0]
    rho = dynamics.lindblad_evolve(rho0, H, [q.lowering(i) for i in range(n)], [0.0, T], rates=[g] * n,
                                   substeps=80)[-1]
    post, p = rho.condition_on(q.algebra.central([1.0 if N == w else 0.0 for N in range(n + 1)]))
    assert abs(p.item() - math.exp(-w * g * T)) < 1e-8
    assert abs(st.trace_of_product(post, ideal).real.item() - 1.0) < 1e-8


def test_vectorised_labels_match_fusion_order():
    from torch_vn_algebra.chains import _labels, _labels_tensor
    for L in range(1, 11):
        ref = _labels(L)
        new = _labels_tensor(L)
        assert [t.tolist() for t in new] == ref


def test_long_chain_sparse_only():
    from torch_vn_algebra import krylov
    ch = SpinChain(22, 'periodic', complex_valued=False)
    H = ch.xxz_sparse(1.0, 1.0, sector=1)
    E, psi, res = krylov.ground_state(H)
    # one magnon on a ferromagnetic ring: E(k) = J L/4 - J + J cos k, lowest at k = pi: L/4 - 2.
    # The spectrum is positive and degenerate (k, -k), so Lanczos breaks down early: this checks that
    # the zero vectors after the breakdown do not produce a spurious eigenvalue 0.
    assert abs(E.item() - (22 / 4 - 2)) < 1e-9
    assert res.item() < 1e-8
