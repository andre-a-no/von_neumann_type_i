"""
Symmetry-resolved entanglement of XXZ ground states (paper, Section "Use case: symmetry-resolved entanglement";
scripts/symmetry_resolved.py): the sector form of rho_A obtained from the fused tensor product against free fermions
(flux insertion), the full Hilbert space, the library's partial traces and QuSpin; the identities behind
equipartition; the Bethe-ansatz Luttinger parameter. CPU, double precision, fixed seeds.
"""
import importlib.util
import math
import runpy
from pathlib import Path

import numpy as np
import pytest
import torch

from torch_vn_algebra import SpinChain, DensityMatrix, composite

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location('symmetry_resolved', ROOT / 'scripts' / 'symmetry_resolved.py')
sr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sr)

pytestmark = [pytest.mark.filterwarnings("ignore::torch_vn_algebra.cost.CostWarning")]


def dense_ground_states(L, boundary, deltas):
    """Ground states of sector L/2 for every Delta from the library's dense sector Hamiltonians."""
    ch = SpinChain(L, boundary, complex_valued=False, precision='double')
    psis = []
    for d in deltas:
        w, V = ch.xxz(J=1.0, Delta=d, sector=L // 2).eigh()[0]
        assert w[0, 1] - w[0, 0] > 1e-6
        psis.append(V[0, :, 0])
    return ch, torch.stack(psis)


def test_luttinger_parameter():
    """K(Delta) = pi / (2 (pi - arccos Delta)): K = 3/2, 1, 3/4, 1/2 at Delta = -1/2, 0, 1/2, 1."""
    for d, K in ((-0.5, 1.5), (0.0, 1.0), (0.5, 0.75), (1.0, 0.5)):
        assert abs(sr.luttinger_K(d) - K) < 1e-14


def test_batched_delta_hamiltonian():
    """The batch of sector Hamiltonians built from one shared hopping part equals xxz_sparse for each Delta."""
    ch = SpinChain(10, 'periodic', complex_valued=False, precision='double')
    deltas = [-0.7, 0.0, 0.4, 1.0]
    H = sr.xxz_batch(ch, 5, deltas)
    v = torch.randn(len(deltas), ch.sector_dim(5), dtype=torch.float64, generator=torch.Generator().manual_seed(1))
    for i, d in enumerate(deltas):
        ref = ch.xxz_sparse(J=1.0, Delta=d, sector=5).matvec(v[i:i + 1])
        assert torch.allclose(H.matvec(v)[i:i + 1], ref, atol=1e-13)


@pytest.mark.parametrize('boundary', ['periodic', 'open'])
def test_xx_chain_against_flux_insertion(boundary):
    """XX chain: p(q), Tr Pi_q rho log rho, Z_n(q), S(q), S_2(q), S_3(q), b_n, B_1, G_1 equal the free-fermion
    flux-insertion values (charged moments from the correlation matrix, exact Fourier transform) for every block."""
    L = 12
    ch, psi = dense_ground_states(L, boundary, [0.0])
    C, E0 = sr.ff_correlation(L, boundary, L // 2)
    for ell in range(1, L):
        st = sr.sector_statistics(sr.Bipartition(ch, L // 2, ell).reduced_state(psi))
        fl = sr.ff_flux(C, ell, pure_projector_rank=L // 2)
        p = st['p'][0].numpy()
        assert np.abs(p - fl['p']).max() < 1e-13
        assert np.abs(st['T'][0].numpy() - fl['T']).max() < 1e-13
        sel = p > 1e-6
        assert np.abs(st['Sq'][0].numpy()[sel] - fl['Sq'][sel]).max() < 1e-8
        for n in (2, 3):
            Z = st['Z'][n][0].numpy()
            assert np.abs(Z - fl['Z'][n]).max() < 1e-13
            seln = sel & (Z / Z.sum() > 1e-6)
            assert np.abs(st['Snq'][n][0].numpy()[seln] - fl['Snq'][n][seln]).max() < 1e-8
            assert abs(st['b'][n][0].item() - sr.ff_variances(fl['nu'], n)[0]) < 1e-12
        assert abs(st['b'][1][0].item() - sr.ff_variances(fl['nu'], 1)[0]) < 1e-12
        if 2 <= ell <= L // 2:
            B1, G1 = sr.ff_central(fl, ell)
            assert abs(st['B'][1][0].item() - B1) < 1e-10 and abs(st['G'][1][0].item() - G1) < 1e-10


def test_free_fermion_ground_energy_and_lanczos():
    """Jordan-Wigner boundary sign (-1)^(N-1): sum of the N lowest single-particle energies equals the batched
    Lanczos ground-state energy of the periodic XX ring for N = L/2 even and odd."""
    for L in (12, 14):
        ch = SpinChain(L, 'periodic', complex_valued=False, precision='double')
        E, psi, res = sr.ground_states(ch, [0.0], 60, 1, 0, log=lambda s: None)
        assert abs(E[0].item() - sr.ff_correlation(L, 'periodic', L // 2)[1]) < 1e-10


@pytest.mark.parametrize('boundary', ['periodic', 'open'])
def test_full_hilbert_space(boundary):
    """For a batch of Delta, p(q) and S(q) from the sector form agree with Pi_q rho_A Pi_q of the 2^ell x 2^ell
    reduced density matrix of the dense 2^L ground state, and rho_A has no part between different charges."""
    L, deltas = 8, [-0.6, 0.0, 0.3, 1.0]
    ch = SpinChain(L, boundary, complex_valued=False, precision='double')
    E, psi, res = sr.ground_states(ch, deltas, 60, 4, 0, log=lambda s: None)
    for i, d in enumerate(deltas):
        v, Ef = sr.dense_full_space(L, boundary, d)
        assert abs(Ef - E[i].item()) < 1e-10
        for ell in range(1, L):
            st = sr.sector_statistics(sr.Bipartition(ch, L // 2, ell).reduced_state(psi[i:i + 1]))
            p, Sq, off = sr.full_space_sectors(v, L, ell)
            assert off < 1e-12
            assert np.abs(st['p'][0].numpy() - p).max() < 1e-10
            sel = p > 1e-6
            assert np.abs(st['Sq'][0].numpy()[sel] - Sq[sel]).max() < 1e-7


def test_three_library_routes_to_rho_A():
    """rho_A from the pure-state pair blocks, from SpinChain.reduced_state of |psi><psi| and from
    composite.partial_trace in the fused product M_A (x) M_B coincide."""
    L = 10
    ch, psi = dense_ground_states(L, 'open', [0.7])
    N = L // 2
    rhoN = DensityMatrix(ch.sector_algebra(N), matrix=(psi[:, :, None] * psi[:, None, :])[:, None], validate=False)
    for ell in range(1, L):
        bp = sr.Bipartition(ch, N, ell)
        r1 = bp.reduced_state(psi).matrix
        r2 = ch.reduced_state(rhoN, ell, sector=N).matrix
        alg12 = bp.alg12
        k = psi.shape[1]
        mat = torch.zeros(1, alg12.C, alg12.k_max, alg12.k_max, dtype=psi.dtype)
        fv = bp.fused_vector(psi)
        mat[:, bp.sector, :k, :k] = fv[:, :, None] * fv[:, None, :]
        r3 = composite.partial_trace(DensityMatrix(alg12, matrix=mat, validate=False), keep=1).matrix
        assert (r1 - r2).abs().max() < 1e-14 and (r1 - r3).abs().max() < 1e-14


def test_identities_behind_equipartition():
    """S = sum_q p(q) S(q) + H(p); S_n(q) - S_n = [log p_n(q) - n log p(q)]/(1 - n); and
    S(q) - S = log p(q) - d/dn log p_n(q) at n = 1 (central difference in n), which gives G_1 = d/dn (n B_n)."""
    L, ell = 12, 6
    ch, psi = dense_ground_states(L, 'periodic', [-0.4, 0.5])
    rho = sr.Bipartition(ch, L // 2, ell).reduced_state(psi)
    st = sr.sector_statistics(rho)
    p = st['p']
    assert torch.allclose(st['S'], (p * st['Sq']).sum(-1) + st['H'], atol=1e-13)
    pn = st['Z'][2] / st['Z'][2].sum(-1, keepdim=True)
    assert torch.allclose(st['Snq'][2] - st['Sn'][2][:, None], -(torch.log(pn) - 2 * torch.log(p)), atol=1e-12)
    ev = [w.clamp(min=0) for w in rho.eigenvalues()]
    h = 1e-5

    def logpn(n):
        Z = torch.stack([(w ** n).sum(-1) for w in ev], -1)
        return torch.log(Z / Z.sum(-1, keepdim=True))
    deriv = (logpn(1 + h) - logpn(1 - h)) / (2 * h)
    assert torch.allclose(st['Sq'] - st['S'][:, None], torch.log(p) - deriv, atol=1e-8)

    def nB(n):
        return n * sr.central_width(torch.log(torch.stack([(w ** n).sum(-1) for w in ev], -1)), ell)
    assert torch.allclose(st['G'][1], (nB(1 + h) - nB(1 - h)) / (2 * h), atol=1e-7)


def test_gaussian_breaking_exact_at_central_sectors():
    """The integer-sampled Gaussian form reproduces S_n(q) - S_n(q0) exactly at q0 +- 1 (definition of B_n), n = 1, 2."""
    L, ell = 12, 6
    ch, psi = dense_ground_states(L, 'periodic', [0.0, 1.0])
    st = sr.sector_statistics(sr.Bipartition(ch, L // 2, ell).reduced_state(psi))
    q0 = ell // 2
    for n in (1, 2):
        Sq = st['Sq'] if n == 1 else st['Snq'][2]
        g = sr.gaussian_breaking(st, n)
        assert torch.allclose(Sq[:, q0 + 1] - Sq[:, q0], g[:, q0 + 1], atol=1e-12)
        assert torch.allclose(Sq[:, q0 - 1] - Sq[:, q0], g[:, q0 - 1], atol=1e-12)


def test_infinite_xx_chain_constants():
    """Infinite XX chain at ell = 400: Var Q_A = (log 2 ell + 1 + gamma_E)/pi^2 + o(1), G_1 = d/dn (n b_n)|_1 close to
    1/(2 pi^2), and S(q0) - S(q0 + 1) within 10 % of (pi^2/4)/(log 2 ell + 1 + gamma_E)^2 (observations, not proofs)."""
    r = sr.part_ff_inf(dict(ff_inf_ells=[400], ff_inf_flux_max=400), log=lambda s: None)[0]
    assert abs(r['b1_dev']) < 1e-6
    assert abs(r['G1'] - 1 / (2 * math.pi ** 2)) < 2e-5
    assert 0.9 < r['dS1'] / r['dS_asym'] < 1.0


def test_quspin_reduced_density_matrix():
    """QuSpin (external): ground state and reduced density matrix of the XXZ ring, projected on the charges of A,
    give the same p(q) and S(q) as the library."""
    pytest.importorskip('quspin')
    cfg = dict(sr.CONFIGS['quick'], L_quspin=10, delta_quspin=0.6)
    out = sr.part_quspin(cfg, 'cpu', 0, log=lambda s: None)
    for r in out:
        assert r['dE'] < 1e-10 and r['dp'] < 1e-10 and r['dS'] < 1e-8


def test_listing_runs():
    """The listing man/article/listings/symmetry_resolved.py runs and satisfies S = sum_q p(q) S(q) + H(p)."""
    runpy.run_path(str(ROOT / 'man' / 'article' / 'listings' / 'symmetry_resolved.py'), run_name='__main__')


def _bernoulli_sector_reference(nu, ns=(2, 3)):
    """Direct occupation-state enumeration, independent of Fourier flux insertion."""
    p = np.zeros(len(nu)+1)
    T = np.zeros_like(p)
    Z = {n: np.zeros_like(p) for n in ns}
    for occupation in range(1 << len(nu)):
        bits = [(occupation >> j) & 1 for j in range(len(nu))]
        probability = math.prod(v if bit else 1-v for v, bit in zip(nu, bits))
        q = sum(bits)
        p[q] += probability
        if probability > 0:
            T[q] += probability*math.log(probability)
        for n in ns:
            Z[n][q] += probability**n
    return p, T, Z


@pytest.mark.parametrize('rotated', [False, True])
def test_flux_uses_only_rank_forced_pure_modes(rotated):
    # Two occupied orbitals in five sites: one entirely in A and one shared with B.
    # The four-site compression has forced multiplicities (zero:2, one:1).
    occupied = np.zeros((5, 2))
    occupied[0, 0] = 1
    occupied[1, 1] = math.sqrt(.3)
    occupied[4, 1] = math.sqrt(.7)
    C = occupied @ occupied.T
    # Mimic eigensolver/projector-formation noise at known structural endpoints.
    eps = np.finfo(float).eps
    C[2, 2] = 2*eps
    C[3, 3] = 4*eps
    C[0, 0] -= 6*eps
    if rotated:
        unitary = np.eye(5)
        unitary[:4, :4] = np.array([[1, 1, 1, 1], [1, -1, 1, -1],
                                   [1, 1, -1, -1], [1, -1, -1, 1]])/2
        C = unitary @ C @ unitary.T
    result = sr.ff_flux(C, 4, pure_projector_rank=2)
    p, T, Z = _bernoulli_sector_reference([0., 0., .3, 1.])
    assert np.array_equal(result['nu'][:2], np.zeros(2))
    assert result['nu'][-1] == 1
    np.testing.assert_allclose(result['p'], p, rtol=0, atol=2e-15)
    np.testing.assert_allclose(result['T'], T, rtol=0, atol=2e-15)
    for n in (2, 3):
        np.testing.assert_allclose(result['Z'][n], Z[n], rtol=0, atol=2e-15)
    raw = sr.ff_flux(C, 4)
    assert -raw['T'].sum() > -T.sum()+1e-14


@pytest.mark.parametrize('rank', [0, 1, np.int64(2), 3, 5])
def test_full_pure_slater_state_has_no_reference_entropy(rank):
    L = 5
    C = np.diag(np.r_[np.ones(rank), np.zeros(L-rank)])
    result = sr.ff_flux(C, L, pure_projector_rank=rank)
    expected_p = np.zeros(L+1)
    expected_p[rank] = 1
    np.testing.assert_allclose(result['p'], expected_p, rtol=0, atol=2e-15)
    assert np.array_equal(result['T'], np.zeros(L+1))
    for n in (2, 3):
        np.testing.assert_allclose(result['Z'][n], expected_p, rtol=0, atol=2e-15)


def test_genuine_small_entangled_modes_are_not_snapped_by_purity_hint():
    small = 1e-15
    # A genuine rank-two projector with no forced endpoints in its two-site
    # compression. Both tiny probabilities are physical and must be retained.
    occupied = np.zeros((4, 2))
    occupied[0, 0], occupied[2, 0] = math.sqrt(small), math.sqrt(1-small)
    occupied[1, 1], occupied[3, 1] = math.sqrt(1-small), math.sqrt(small)
    C = occupied @ occupied.T
    nu = np.linalg.eigvalsh(C[:2, :2])
    result = sr.ff_flux(C, 2, pure_projector_rank=2)
    assert np.array_equal(result['nu'], nu)
    _, T, _ = _bernoulli_sector_reference(nu)
    np.testing.assert_allclose(result['T'], T, rtol=0, atol=1e-27)
    assert -result['T'].sum() > 1e-14


def test_mixed_reference_keeps_genuine_near_endpoint_occupations():
    nu = np.array([1e-15, 1-1e-15])
    result = sr.ff_flux(np.diag(nu), 2)
    assert np.array_equal(result['nu'], nu)
    _, T, _ = _bernoulli_sector_reference(nu)
    np.testing.assert_allclose(result['T'], T, rtol=0, atol=1e-27)
    assert -result['T'].sum() > 1e-14


def test_flux_rejects_inconsistent_forced_pure_endpoints():
    with pytest.raises(ValueError, match='inconsistent'):
        sr.ff_flux(np.diag([.2, .8]), 2, pure_projector_rank=1)


@pytest.mark.parametrize('rank', [-1, 3, 1.5])
def test_flux_rejects_invalid_pure_rank(rank):
    with pytest.raises(ValueError, match='integer rank'):
        sr.ff_flux(np.eye(2), 2, pure_projector_rank=rank)


def test_wide_pure_schmidt_factor_retains_tiny_mode_and_exact_kernel():
    ch = SpinChain(6, 'open', complex_valued=False, precision='double')
    bp = sr.Bipartition(ch, 3, 4)
    q = 2
    k, m = bp.gather[q].shape
    assert (k, m) == (6, 2)
    basis, _ = torch.linalg.qr(torch.randn(k, m, dtype=torch.float64,
                                          generator=torch.Generator().manual_seed(431)))
    small = 1e-14
    P = basis*torch.sqrt(torch.tensor([1-small, small], dtype=torch.float64))[None]
    psi = torch.zeros(1, ch.sector_dim(3), dtype=torch.float64)
    psi[0, bp.gather[q]] = P
    rho = bp.reduced_state(psi)
    values = rho.eigenvalues()[q][0]
    assert torch.equal(values[:k-m], torch.zeros(k-m, dtype=torch.float64))
    torch.testing.assert_close(values[-2], torch.tensor(small, dtype=torch.float64), rtol=1e-10, atol=0)
    st = sr.sector_statistics(rho)
    reference = torch.tensor([1-small, small], dtype=torch.float64)
    torch.testing.assert_close(st['T'][0, q], torch.special.xlogy(reference, reference).sum(), rtol=0, atol=2e-15)


@pytest.mark.parametrize('complex_valued', [False, True])
def test_schmidt_spectrum_cache_preserves_gradients_and_matrix_invalidation(complex_valued):
    ch = SpinChain(6, 'open', complex_valued=complex_valued, precision='double')
    bp = sr.Bipartition(ch, 3, 4)
    psi = torch.randn(2, ch.sector_dim(3), dtype=ch.algebra.dtype,
                      generator=torch.Generator().manual_seed(583))
    psi = (psi/torch.linalg.vector_norm(psi, dim=-1, keepdim=True)).requires_grad_(True)
    rho = bp.reduced_state(psi)
    assert rho._eigenvalues is not None
    spectral_purity = sum(w.square().sum() for w in rho.eigenvalues())
    dense_purity = rho.matrix.abs().square().sum()
    spectral_grad, = torch.autograd.grad(spectral_purity, psi, retain_graph=True)
    dense_grad, = torch.autograd.grad(dense_purity, psi)
    torch.testing.assert_close(spectral_purity, dense_purity)
    torch.testing.assert_close(spectral_grad, dense_grad)

    q, k = 2, bp.left.algebra.k_factors[2]
    with torch.no_grad():
        rho.matrix[:, q, :k, :k].add_(.1*torch.eye(k, dtype=psi.dtype))
    reference = torch.linalg.eigvalsh(rho.matrix[:, q, :k, :k])
    torch.testing.assert_close(rho.eigenvalues()[q], reference)
    assert rho._eigenvalues is None
