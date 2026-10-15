"""
Cross-checks of torch_vn_algebra against independent, established software on identical inputs.

Reference values in tests/data/crosscheck/*.json were produced by scripts/make_crosscheck_data.py (which does not
import the library) with
    QuTiP   (mesolve, liouvillian, ptrace, entropy_vn, entropy_mutual, fidelity, entropy_relative, tracedist),
    QuSpin  (exact diagonalisation of XXZ chains in magnetisation and momentum sectors, ground-state correlations,
             ent_entropy, hamiltonian.evolve),
    GAP     (character tables, Frobenius-Schur indicators, ScalarProduct of permutation / induced characters),
    SymPy   (sympy.physics.wigner.clebsch_gordan, exact).
The stored-data tests need only numpy and torch. The "live" tests recompute a sample with the external package
(imported by pytest.importorskip: skipped when it is missing).
"""
import json
import math
import pathlib
from fractions import Fraction

import numpy as np
import pytest
import torch

from torch_vn_algebra import (TypeIAlgebra, DensityMatrix, InterSectorChannel, SpinChain, SU2Chain,
                              tensor_product, fused_tensor_product, kron, partial_trace)
from torch_vn_algebra import composite, decompose, divergences, dynamics, krylov, states
from torch_vn_algebra.su2 import clebsch_gordan, clebsch_gordan_exact

DATA = pathlib.Path(__file__).resolve().parent / 'data' / 'crosscheck'


def load(name):
    with open(DATA / name) as f:
        return json.load(f)


def C(d):
    """Complex array from {"re": ..., "im": ...}."""
    return np.asarray(d["re"], dtype=float) + 1j * np.asarray(d["im"], dtype=float)


def offsets(dims):
    return [int(x) for x in np.cumsum([0] + list(dims))]


def blocks_of(X, index_lists):
    """Diagonal blocks X[idx, idx] of a full matrix, one per index list."""
    return [torch.as_tensor(np.ascontiguousarray(X[np.ix_(idx, idx)])) for idx in index_lists]


def full_of(op, index_lists, n, b=0):
    """Full n x n matrix of an operator of an algebra whose sector c sits at index_lists[c]."""
    M = np.zeros((n, n), dtype=complex)
    mat = op.matrix[b].detach().cpu().numpy()
    for c, idx in enumerate(index_lists):
        k = len(idx)
        M[np.ix_(idx, idx)] = mat[c, :k, :k]
    return M


def off_block_weight(X, index_lists):
    inside = np.zeros(X.shape, dtype=bool)
    for idx in index_lists:
        inside[np.ix_(idx, idx)] = True
    return float(np.abs(X[~inside]).max()) if (~inside).any() else 0.0


# ==========================================================================================================
# (1) Lindblad dynamics vs QuTiP mesolve / liouvillian
# ==========================================================================================================
@pytest.fixture(scope='module')
def lb():
    return load('lindblad.json')['cases']


def _block_case(case):
    dims = case['dims']
    o = offsets(dims)
    idx = [list(range(o[c], o[c + 1])) for c in range(len(dims))]
    alg = TypeIAlgebra(dims, dims, precision='double')

    def op(X):
        return alg.from_blocks(blocks_of(X, idx))
    return alg, idx, op


def test_lindblad_evolve_matches_qutip_mesolve_on_direct_sum(lb):
    """lindblad_evolve on M_2 (+) M_3 (+) M_1 equals QuTiP mesolve of the same GKSL equation on C^6."""
    case = lb['blocks']
    alg, idx, op = _block_case(case)
    rho0 = DensityMatrix(alg, matrix=op(C(case['rho0'])).matrix)
    H = op(C(case['H0']))
    jumps = [op(C(L)) for L in case['jumps']]
    out = dynamics.lindblad_evolve(rho0, H, jumps, times=case['times'], rates=case['rates'], substeps=40)
    for rho_t, ref in zip(out, case['rho_t']):
        ref = C(ref)
        assert off_block_weight(ref, idx) < 1e-12                  # QuTiP keeps the block structure
        assert np.abs(full_of(rho_t, idx, ref.shape[0]) - ref).max() < 1e-9


def test_lindblad_evolve_time_dependent_matches_qutip(lb):
    """lindblad_evolve with H(t) = H0 + cos(omega t) H1 (a callable) equals QuTiP mesolve with a coefficient."""
    case = lb['blocks']
    alg, idx, op = _block_case(case)
    rho0 = DensityMatrix(alg, matrix=op(C(case['rho0'])).matrix)
    H0, H1, w = op(C(case['H0'])), op(C(case['H1'])), case['omega']
    jumps = [op(C(L)) for L in case['jumps']]
    out = dynamics.lindblad_evolve(rho0, lambda t: H0 + H1 * math.cos(w * t), jumps, times=case['times'],
                                   rates=case['rates'], substeps=40)
    for rho_t, ref in zip(out, case['rho_t_timedep']):
        ref = C(ref)
        assert np.abs(full_of(rho_t, idx, ref.shape[0]) - ref).max() < 1e-9


def test_lindblad_channel_matches_qutip_liouvillian_expm(lb):
    """lindblad_channel(t) maps every matrix unit E_ij of every block as QuTiP's expm(t * liouvillian)."""
    case = lb['blocks']
    alg, idx, op = _block_case(case)
    H = op(C(case['H0']))
    jumps = [op(C(L)) for L in case['jumps']]
    ch = dynamics.lindblad_channel(alg, H, jumps, t=case['channel_t'], rates=case['rates'])
    n = sum(case['dims'])
    for (c, i, j), ref in zip(case['channel_units'], case['channel_images']):
        ref = C(ref)
        E = np.zeros((n, n), dtype=complex)
        E[idx[c][i], idx[c][j]] = 1.0
        img = full_of(ch.apply(op(E)), idx, n)
        assert off_block_weight(ref, [idx[c]]) < 1e-12             # exp(tL) maps block c into block c
        assert np.abs(img - ref).max() < 1e-11


@pytest.mark.filterwarnings('ignore::torch_vn_algebra.cost.CostWarning')
def test_lindblad_chain_with_losses_matches_qutip_full_space(lb):
    """XXZ chain (L = 4) with single-site and collective losses (jumps N -> N-1), dephasing and a hopping jump:
    lindblad_evolve in (+)_N M_binom(L,N) equals QuTiP mesolve in the full 2^L space."""
    case = lb['chain']
    L = case['L']
    chain = SpinChain(L, case['boundary'], precision='double')
    alg = chain.algebra
    idx = [chain.labels[N].tolist() for N in range(L + 1)]
    n = 2 ** L
    H = chain.xxz(J=case['J'], Delta=case['Delta'], h=case['h'])
    assert np.abs(full_of(H, idx, n) - C(case['H'])).max() < 1e-13         # same Hamiltonian, same basis
    jumps, rates = [], []
    for jp in case['jumps']:
        M = C(jp['matrix'])
        if jp['kind'] == 'lowering':
            K = chain.lowering(jp['sites'])
            for N in range(1, L + 1):                                         # the library's S^-_i is QuTiP's
                blk = M[np.ix_(idx[N - 1], idx[N])]
                assert np.abs(K.kraus[0, 0, N - 1, N, :len(idx[N - 1]), :len(idx[N])].numpy() - blk).max() < 1e-14
        elif jp['kind'] == 'collective_lowering':
            K = InterSectorChannel.from_blocks(alg, alg, {(N - 1, N): [M[np.ix_(idx[N - 1], idx[N])]]
                                                          for N in range(1, L + 1)})
        else:                                                                 # sector preserving: an Operator
            assert off_block_weight(M, idx) == 0.0
            K = alg.from_blocks(blocks_of(M, idx))
        jumps.append(K)
        rates.append(jp['rate'])
    rho0 = DensityMatrix(alg, matrix=alg.from_blocks(blocks_of(C(case['rho0']), idx)).matrix)
    out = dynamics.lindblad_evolve(rho0, H, jumps, times=case['times'], rates=rates, substeps=30)   # RK4 error 3e-10
    assert case['max_off_sector'] < 1e-12               # the full-space evolution never leaves the algebra
    for rho_t, ref in zip(out, case['rho_t']):
        ref = C(ref)
        assert np.abs(full_of(rho_t, idx, n) - ref).max() < 1e-9


def test_lindblad_real_algebra_matches_qutip(lb):
    """On a real algebra M_3(R) (+) M_2(R) with real jumps and H = None, lindblad_evolve and lindblad_channel equal
    QuTiP mesolve / expm(liouvillian), and stay real."""
    case = lb['real_blocks']
    dims = case['dims']
    o = offsets(dims)
    idx = [list(range(o[c], o[c + 1])) for c in range(len(dims))]
    alg = TypeIAlgebra(dims, dims, complex_valued=False, precision='double')
    op = lambda X: alg.from_blocks(blocks_of(np.asarray(X, dtype=float), idx))        # noqa: E731
    rho0 = DensityMatrix(alg, matrix=op(case['rho0']).matrix)
    jumps = [op(L) for L in case['jumps']]
    out = dynamics.lindblad_evolve(rho0, None, jumps, times=case['times'], rates=case['rates'], substeps=40)
    for rho_t, ref in zip(out, case['rho_t']):
        assert not torch.is_complex(rho_t.matrix)
        assert np.abs(full_of(rho_t, idx, 5) - C(ref)).max() < 1e-9
    ch = dynamics.lindblad_channel(alg, None, jumps, t=case['channel_t'], rates=case['rates'])
    assert np.abs(full_of(ch.apply(rho0), idx, 5) - C(case['channel_rho0'])).max() < 1e-12


def test_lindblad_channel_as_jump_matches_qutip_kraus_jumps(lb):
    """A channel passed as a jump (SpinChain.site_amplitude_damping, an InterSectorChannel with Kraus operators
    K_0, K_1) generates Phi(rho) - 1/2 {Phi^*(1), rho}: equal to QuTiP mesolve with c_ops K_0, K_1."""
    case = lb['channel_jump']
    L = case['L']
    chain = SpinChain(L, case['boundary'], precision='double')
    alg = chain.algebra
    idx = [chain.labels[N].tolist() for N in range(L + 1)]
    n = 2 ** L
    Phi = chain.site_amplitude_damping(case['site'], case['gamma'])
    K0, K1 = np.asarray(case['K0']), np.asarray(case['K1'])
    for N in range(L + 1):                              # the library's Kraus blocks are those of the full-space K_i
        k = len(idx[N])
        assert np.abs(Phi.kraus[0, 0, N, N, :k, :k].numpy() - K0[np.ix_(idx[N], idx[N])]).max() < 1e-15
        if N:
            kk = len(idx[N - 1])
            assert np.abs(Phi.kraus[0, 1, N - 1, N, :kk, :k].numpy() - K1[np.ix_(idx[N - 1], idx[N])]).max() < 1e-15
    rho0 = DensityMatrix(alg, matrix=alg.from_blocks(blocks_of(C(case['rho0']), idx)).matrix)
    H = chain.xxz(case['J'], case['Delta'], case['h'])
    out = dynamics.lindblad_evolve(rho0, H, [Phi], times=case['times'], rates=[case['rate']], substeps=80)
    for rho_t, ref in zip(out, case['rho_t']):
        assert np.abs(full_of(rho_t, idx, n) - C(ref)).max() < 1e-9


def test_live_qutip_mesolve_chain_with_loss():
    """Live QuTiP: mesolve of a 3-site XXZ chain with loss on site 1 equals lindblad_evolve (jump between sectors)."""
    qutip = pytest.importorskip('qutip')
    L = 3
    chain = SpinChain(L, 'periodic', precision='double')
    up, dn = qutip.basis(2, 0), qutip.basis(2, 1)
    sp1, sz1, eye = up * dn.dag(), 0.5 * (up * up.dag() - dn * dn.dag()), qutip.qeye(2)
    site = lambda o, i: qutip.tensor([o if j == i else eye for j in range(L)])      # noqa: E731
    Hq = sum(0.5 * (site(sp1, i) * site(sp1, (i + 1) % L).dag() + site(sp1, i).dag() * site(sp1, (i + 1) % L))
             + 0.7 * site(sz1, i) * site(sz1, (i + 1) % L) for i in range(L)) + 0.4 * site(sz1, 0)
    psi = qutip.tensor(up, up, dn)                                        # sites 0, 1 up: label 0b011
    times = [0.0, 0.6, 1.2]
    res = qutip.mesolve(Hq, psi * psi.dag(), times, [math.sqrt(0.5) * site(sp1, 1).dag()],
                        options=dict(atol=1e-12, rtol=1e-11, nsteps=10 ** 6))
    perm = [sum((1 - ((x >> i) & 1)) << (L - 1 - i) for i in range(L)) for x in range(2 ** L)]
    idx = [chain.labels[N].tolist() for N in range(L + 1)]
    h = torch.zeros(L, dtype=torch.float64)
    h[0] = 0.4
    out = dynamics.lindblad_evolve(chain.basis_state('110'), chain.xxz(1.0, 0.7, h), [chain.lowering(1)],
                                   times=times, rates=[0.5], substeps=60)
    for rho_t, ref in zip(out, res.states):
        R = ref.full()[np.ix_(perm, perm)]
        assert np.abs(full_of(rho_t, idx, 2 ** L) - R).max() < 1e-9


# ==========================================================================================================
# (2) tensor products, partial traces and entropies vs QuTiP
# ==========================================================================================================
@pytest.fixture(scope='module')
def cp():
    return load('composite.json')['cases']


def _product(case):
    k1, k2 = case['k1'], case['k2']
    alg1 = TypeIAlgebra(k1, k1, precision='double', charges=case['charges1'])
    alg2 = TypeIAlgebra(k2, k2, precision='double', charges=case['charges2'])
    alg12 = (tensor_product if case['product'] == 'plain' else fused_tensor_product)(alg1, alg2)
    n2 = sum(k2)
    o1, o2 = offsets(k1), offsets(k2)
    idx12 = []
    for pairs in alg12.product_blocks:                      # basis vector off + a m_d + b = e_a (x) e_b
        idx = []
        for c, d, off in pairs:
            for a in range(k1[c]):
                for b in range(k2[d]):
                    idx.append((o1[c] + a) * n2 + o2[d] + b)
        idx12.append(idx)
    idx1 = [list(range(o1[c], o1[c + 1])) for c in range(len(k1))]
    idx2 = [list(range(o2[d], o2[d + 1])) for d in range(len(k2))]
    return alg1, alg2, alg12, idx1, idx2, idx12


@pytest.mark.parametrize('name', ['plain', 'fused', 'fused_repeated_charges'])
def test_partial_trace_and_kron_match_qutip(cp, name):
    """partial_trace (direct and as a channel) and kron equal QuTiP ptrace / tensor on H1 (x) H2; with repeated
    charges the kept factor gets the conditional expectation onto its algebra of QuTiP's partial trace."""
    case = cp[name]
    alg1, alg2, alg12, idx1, idx2, idx12 = _product(case)
    n1, n2 = sum(case['k1']), sum(case['k2'])
    rho_full, X_full = C(case['rho']), C(case['X'])
    assert off_block_weight(rho_full, idx12) == 0.0 and off_block_weight(X_full, idx12) == 0.0
    rho = DensityMatrix(alg12, matrix=alg12.from_blocks(blocks_of(rho_full, idx12)).matrix)
    X = alg12.from_blocks(blocks_of(X_full, idx12))
    for keep, idx, n in ((1, idx1, n1), (2, idx2, n2)):
        for what, op in (('rho', rho), ('X', X)):
            ref = C(case[f'ptrace_{what}_keep{keep}'])
            expected = np.zeros_like(ref)                   # conditional expectation onto the kept algebra
            for i in idx:
                expected[np.ix_(i, i)] = ref[np.ix_(i, i)]
            if name == 'fused_repeated_charges' and keep == 1:
                assert off_block_weight(ref, idx) > 1e-2      # the Hilbert-space partial trace leaves the algebra
            else:
                assert off_block_weight(ref, idx) < 1e-14
            got = full_of(partial_trace(op, keep=keep), idx, n)
            assert np.abs(got - expected).max() < 1e-13
            got_ch = full_of(composite.partial_trace_channel(alg12, keep).apply(op), idx, n)
            assert np.abs(got_ch - expected).max() < 1e-13
    A = alg1.from_blocks(blocks_of(C(case['A']), idx1))
    B = alg2.from_blocks(blocks_of(C(case['B']), idx2))
    AB = C(case['kron_AB'])
    assert off_block_weight(AB, idx12) == 0.0
    assert np.abs(full_of(kron(A, B, alg12), idx12, n1 * n2) - AB).max() < 1e-13


@pytest.mark.parametrize('name', ['plain', 'fused', 'fused_repeated_charges'])
def test_entropies_fidelity_relative_entropy_match_qutip(cp, name):
    """von Neumann entropies (global, reduced, mutual information), Uhlmann fidelity, Umegaki relative entropy,
    trace distance and purity equal QuTiP's entropy_vn, entropy_mutual, fidelity^2, entropy_relative, tracedist."""
    case = cp[name]
    alg1, alg2, alg12, idx1, idx2, idx12 = _product(case)
    rho = DensityMatrix(alg12, matrix=alg12.from_blocks(blocks_of(C(case['rho']), idx12)).matrix)
    sigma = DensityMatrix(alg12, matrix=alg12.from_blocks(blocks_of(C(case['sigma']), idx12)).matrix)
    r1, r2 = partial_trace(rho, 1), partial_trace(rho, 2)
    S, S1, S2 = (float(states.von_neumann_entropy(x)) for x in (rho, r1, r2))
    tol = 1e-12
    assert abs(float(states.von_neumann_entropy(sigma)) - case['entropy_sigma']) < tol
    assert abs(S - case['entropy_rho']) < tol
    assert abs(S2 - case['entropy_rho2']) < tol
    if name == 'fused_repeated_charges':
        # the restriction to alg1 is the compression of the Hilbert-space marginal: higher entropy (Lemma: E_M
        # increases entropy), so QuTiP's value is a lower bound
        assert S1 > case['entropy_rho1'] + 1e-6
    else:
        assert abs(S1 - case['entropy_rho1']) < tol
        assert abs(S1 + S2 - S - case['mutual_information']) < tol
        s1 = partial_trace(sigma, 1)
        assert abs(float(states.relative_entropy(r1, s1)) - case['relative_entropy_reduced1']) < 1e-10
    assert abs(float(states.fidelity(rho, sigma)) - case['fidelity_squared']) < 1e-12
    assert abs(float(states.relative_entropy(rho, sigma)) - case['relative_entropy']) < 1e-10
    assert abs(float(divergences.umegaki_relative_entropy(rho, sigma)) - case['relative_entropy']) < 1e-10
    assert abs(float(states.trace_distance(rho, sigma)) - case['trace_distance']) < 1e-12
    assert abs(float(states.purity(rho)) - case['purity']) < 1e-13


def test_live_qutip_ptrace_and_entropies():
    """Live QuTiP: partial trace, entropy, fidelity and relative entropy of random states on M_2 (x) (M_1 (+) M_2)."""
    qutip = pytest.importorskip('qutip')
    rng = np.random.default_rng(5)
    alg1 = TypeIAlgebra([2], [2], precision='double')
    alg2 = TypeIAlgebra([1, 2], [1, 2], precision='double')
    alg12 = tensor_product(alg1, alg2)
    blocks = []
    for k in alg12.k_factors:
        G = rng.normal(size=(2, k, k)) + 1j * rng.normal(size=(2, k, k))
        blocks.append(G @ G.conj().transpose(0, 2, 1))
    tot = sum(np.trace(b, axis1=1, axis2=2).real for b in blocks)
    mat = alg12.from_blocks([torch.as_tensor(b / tot[:, None, None]) for b in blocks]).matrix
    rho = DensityMatrix(alg12, matrix=mat)
    _, _, _, idx1, idx2, idx12 = _product(dict(k1=[2], k2=[1, 2], charges1=[0], charges2=[0, 1], product='plain'))
    for b in range(2):
        R = qutip.Qobj(full_of(rho, idx12, 6, b), dims=[[2, 3], [2, 3]])
        assert np.abs(full_of(partial_trace(rho, 1), idx1, 2, b) - R.ptrace(0).full()).max() < 1e-13
        assert np.abs(full_of(partial_trace(rho, 2), idx2, 3, b) - R.ptrace(1).full()).max() < 1e-13
        assert abs(float(states.von_neumann_entropy(rho)[b]) - qutip.entropy_vn(R)) < 1e-12
    R0 = qutip.Qobj(full_of(rho, idx12, 6, 0))
    R1 = qutip.Qobj(full_of(rho, idx12, 6, 1))
    rho0 = DensityMatrix(alg12, matrix=rho.matrix[:1])
    rho1 = DensityMatrix(alg12, matrix=rho.matrix[1:])
    assert abs(float(states.fidelity(rho0, rho1)) - qutip.fidelity(R0, R1) ** 2) < 1e-12
    assert abs(float(states.relative_entropy(rho0, rho1)) - qutip.entropy_relative(R0, R1)) < 1e-10


# ==========================================================================================================
# (3) XXZ chains vs QuSpin (magnetisation and momentum sectors), J1-J2 ring vs QuSpin (SU(2) blocks)
# ==========================================================================================================
@pytest.fixture(scope='module')
def xx():
    return load('xxz.json')['cases']


def test_xxz_sector_spectra_match_quspin(xx):
    """Spectra of the XXZ chain (open, L = 10, random fields) in the sectors N = 3, 4, 5 equal QuSpin's."""
    case = xx['open10']
    chain = SpinChain(case['L'], case['boundary'], precision='double')
    h = torch.tensor(case['h'], dtype=torch.float64)
    for N, ref in case['spectra'].items():
        H = chain.xxz(case['J'], case['Delta'], h, sector=int(N))
        ev = H.eigenvalues()[0][0].numpy()
        assert ev.shape == (len(ref),)
        assert np.abs(np.sort(ev) - np.asarray(ref)).max() < 1e-11


def _correlations(chain, psi, N):
    """<S^z_i S^z_j> and <S^+_i S^-_j> of a vector of sector N (library basis: labels of SpinChain)."""
    lab = chain.labels[N]
    L = chain.L
    bits = chain.bits(N).numpy().astype(float)
    p = np.abs(psi) ** 2
    zz = np.einsum('a,ai,aj->ij', p, bits - 0.5, bits - 0.5)
    pm = np.zeros((L, L))
    labn = lab.numpy()
    index = chain.index_of.numpy()
    for i in range(L):
        for j in range(L):
            if i == j:
                pm[i, i] = p @ bits[:, i]
                continue
            sel = (((labn >> j) & 1) == 1) & (((labn >> i) & 1) == 0)       # S^+_i S^-_j |a> = |a ^ bits i, j>
            src = np.nonzero(sel)[0]
            dst = index[labn[src] ^ ((1 << i) | (1 << j))]
            pm[i, j] = float(np.real(np.vdot(psi[dst], psi[src])))
    return zz, pm


def test_xxz_krylov_ground_state_correlations_entanglement_match_quspin(xx):
    """krylov.ground_state on xxz_sparse (N = 5): energy, <S^z_i S^z_j>, <S^+_i S^-_j> and the entanglement
    entropies of all left blocks (SpinChain.entanglement_entropy) equal QuSpin's (different basis convention)."""
    case = xx['open10']
    ref = case['ground_state']
    N = ref['N']
    chain = SpinChain(case['L'], case['boundary'], precision='double')
    h = torch.tensor(case['h'], dtype=torch.float64)
    Hs = chain.xxz_sparse(case['J'], case['Delta'], h, sector=N)
    E0, psi, res = krylov.ground_state(Hs, generator=torch.Generator().manual_seed(0))
    assert abs(float(E0[0]) - ref['E0']) < 1e-10
    psi = psi[0].numpy()
    zz, pm = _correlations(chain, psi, N)
    assert np.abs(zz - np.asarray(ref['SzSz'])).max() < 1e-8
    assert np.abs(pm - np.asarray(ref['SpSm'])).max() < 1e-8
    alg = chain.sector_algebra(N)
    v = torch.as_tensor(psi)
    rho = DensityMatrix(alg, matrix=(v[:, None] * v.conj()[None, :])[None, None])
    for n_left, S_ref in zip(range(1, case['L']), ref['entanglement']):
        assert abs(float(chain.entanglement_entropy(rho, n_left, sector=N)) - S_ref) < 1e-8


def test_xxz_krylov_evolution_matches_quspin(xx):
    """krylov.evolve of the Neel state (N = 5): <S^z_i(t)> and the half-chain entanglement entropy equal those of
    QuSpin's hamiltonian.evolve."""
    case = xx['open10']
    ev = case['evolution']
    N, L = ev['N'], case['L']
    chain = SpinChain(L, case['boundary'], precision='double')
    h = torch.tensor(case['h'], dtype=torch.float64)
    Hs = chain.xxz_sparse(case['J'], case['Delta'], h, sector=N)
    psi0 = chain.vector_in_sector(ev['bits'])[None]
    psit = krylov.evolve(Hs, psi0, ev['times'], m=30, tol=1e-10)
    bits = chain.bits(N).numpy().astype(float) - 0.5
    alg = chain.sector_algebra(N)
    for t, (Sz_ref, S_ref) in enumerate(zip(ev['Sz'], ev['entanglement_half'])):
        v = psit[t, 0]
        p = (v.abs() ** 2).numpy()
        assert np.abs(p @ bits - np.asarray(Sz_ref)).max() < 1e-8
        rho = DensityMatrix(alg, matrix=(v[:, None] * v.conj()[None, :])[None, None])
        assert abs(float(chain.entanglement_entropy(rho, L // 2, sector=N)) - S_ref) < 1e-8


def test_xxz_momentum_sectors_match_quspin(xx):
    """Periodic XXZ ring (L = 12, N = 5, 6): dimensions and spectra of the momentum blocks (xxz_momentum) equal
    QuSpin's kblock bases, and the sector spectrum equals QuSpin's (the spectra of k and -k coincide for a real
    H, so the sign convention of the momentum is not tested)."""
    case = xx['ring12']
    L = case['L']
    chain = SpinChain(L, case['boundary'], precision='double')
    for N, ref in case['sectors'].items():
        N = int(N)
        dims = chain.momentum_dims(N)
        assert dims == [len(ref['momentum_spectra'][str(m)]) for m in range(L)]
        Hk = chain.xxz_momentum(case['J'], case['Delta'], case['h'], sector=N)
        evs = Hk.eigenvalues()
        for c, m in enumerate(Hk.algebra.charges):
            got = np.sort(evs[c][0].numpy())
            assert np.abs(got - np.asarray(ref['momentum_spectra'][str(m)])).max() < 1e-11
        H = chain.xxz(case['J'], case['Delta'], case['h'], sector=N)
        assert np.abs(np.sort(H.eigenvalues()[0][0].numpy()) - np.asarray(ref['spectrum'])).max() < 1e-11


def test_xxz_momentum_blocks_all_small_rings_match_quspin(xx):
    """For every periodic ring L = 2..8 and every sector N, the momentum blocks of xxz_momentum have QuSpin's
    dimensions and spectra (QuSpin kblock bases)."""
    case = xx['momentum_sweep']
    for key, blocks in case['blocks'].items():
        L, N = (int(x) for x in key.split(','))
        chain = SpinChain(L, 'periodic', precision='double')
        assert chain.momentum_dims(N) == [len(blocks[str(m)]) for m in range(L)], key
        Hk = chain.xxz_momentum(case['J'], case['Delta'], case['h'], sector=N)
        evs = Hk.eigenvalues()
        assert sorted(Hk.algebra.charges) == [m for m in range(L) if blocks[str(m)]]
        for c, m in enumerate(Hk.algebra.charges):
            assert np.abs(np.sort(evs[c][0].numpy()) - np.asarray(blocks[str(m)])).max() < 1e-11, (key, m)


def test_su2_heisenberg_chains_all_small_match_quspin(xx):
    """SU2Chain.full_spectrum of heisenberg(J=1, h, J2) equals QuSpin's full spectrum for L = 1..8, open and
    periodic, J2 = 0 and 0.35, with the same bonds (SU2Chain.bonds / next_nearest_bonds)."""
    case = xx['heisenberg_sweep']
    for key, ref in case['chains'].items():
        L, bc, J2 = key.split(',')
        L, J2 = int(L), float(J2)
        chain = SU2Chain(L, bc)
        assert sorted(tuple(sorted(b)) for b in chain.bonds) == sorted(tuple(sorted(b)) for b in ref['nn_bonds'])
        if J2:
            assert sorted(tuple(sorted(b)) for b in chain.next_nearest_bonds) == \
                sorted(tuple(sorted(b)) for b in ref['nnn_bonds']), key
        E = chain.full_spectrum(chain.heisenberg(J=case['J'], h=case['h'], J2=J2))[0].numpy()
        assert E.shape == (2 ** L,) and np.abs(E - np.asarray(ref['spectrum'])).max() < 1e-11, key


def test_su2_j1j2_ring_spectrum_and_correlations_match_quspin(xx):
    """J1-J2 Heisenberg ring (L = 10, J2 = 0.2): SU2Chain.full_spectrum (2^L levels from the multiplicity blocks)
    and the ground-state <S_0 . S_r> (Wigner-Eckart and direct) equal QuSpin's U(1)-sector diagonalisation."""
    case = xx['j1j2_ring10']
    chain = SU2Chain(case['L'], 'periodic')
    H = chain.heisenberg(J=case['J1'], J2=case['J2'])
    E = chain.full_spectrum(H)[0].numpy()
    assert np.abs(E - np.asarray(case['spectrum'])).max() < 1e-11
    gs = chain.ground_state(H)
    assert abs(float(gs.energy[0]) - case['E0']) < 1e-11 and float(gs.spin[0]) == 0.0
    S = float(gs.spin[0])
    coef = gs.coefficients[0, :chain.multiplicity(S)]
    pairs = [(0, r) for r in range(1, case['L'] // 2 + 1)]
    for method in ('wigner_eckart', 'direct'):
        corr = chain.expectation(chain.spin_correlation(pairs, method=method), S, coef).numpy()
        assert np.abs(corr - np.asarray(case['S0_dot_Sr'])).max() < 1e-10


def test_live_quspin_sector_spectrum():
    """Live QuSpin: spectrum of a periodic XXZ ring (L = 8, N = 3, staggered field) equals SpinChain.xxz."""
    pytest.importorskip('quspin')
    from quspin.basis import spin_basis_1d
    from quspin.operators import hamiltonian
    L, N, J, Delta = 8, 3, 0.9, 1.7
    h = [0.25 * (-1) ** i for i in range(L)]
    bonds = [(i, (i + 1) % L) for i in range(L)]
    static = [["+-", [[J / 2, i, j] for i, j in bonds]], ["-+", [[J / 2, i, j] for i, j in bonds]],
              ["zz", [[J * Delta, i, j] for i, j in bonds]], ["z", [[h[i], i] for i in range(L)]]]
    basis = spin_basis_1d(L, Nup=N, pauli=0)
    Hq = hamiltonian(static, [], basis=basis, dtype=np.float64, check_symm=False, check_herm=False,
                     check_pcon=False)
    ref = np.linalg.eigvalsh(Hq.toarray())
    chain = SpinChain(L, 'periodic', precision='double')
    got = np.sort(chain.xxz(J, Delta, torch.tensor(h, dtype=torch.float64), sector=N).eigenvalues()[0][0].numpy())
    assert np.abs(got - ref).max() < 1e-11


# ==========================================================================================================
# (4) Wedderburn decomposition of group representations vs GAP character tables
# ==========================================================================================================
@pytest.fixture(scope='module')
def gp():
    return load('groups.json')['groups']


class _Group:
    def __init__(self, g):
        self.perms = [tuple(p) for p in g['perms']]
        self.index = {p: i for i, p in enumerate(self.perms)}
        self.n = len(self.perms)
        self.classes = g['classes']
        self.irr = np.array([[complex(*v) for v in row] for row in g['irr']])          # (irr, classes)
        self.indicators = g['indicators']
        self.deg = g['deg']

    def mul(self, a, b):
        """GAP's product a * b: apply a, then b."""
        pa, pb = self.perms[a], self.perms[b]
        return self.index[tuple(pb[x] for x in pa)]

    def inv(self, a):
        p = self.perms[a]
        q = [0] * len(p)
        for i, x in enumerate(p):
            q[x] = i
        return self.index[tuple(q)]

    def power(self, a, k):
        r = self.index[tuple(range(self.deg))]
        for _ in range(k):
            r = self.mul(r, a)
        return r

    # ----- representations (homomorphisms g -> unitary matrices) -----
    def natural(self):
        M = np.zeros((self.n, self.deg, self.deg))
        for g, p in enumerate(self.perms):
            M[g, np.arange(self.deg), list(p)] = 1.0          # P(g)[j, g(j)] = 1, i.e. P(g) e_i = e_{g^-1(i)}
        return M.astype(complex)

    def regular(self):
        M = np.zeros((self.n, self.n, self.n))
        for g in range(self.n):
            for x in range(self.n):
                M[g, self.mul(g, x), x] = 1.0                 # L(g) e_x = e_{g x}
        return M.astype(complex)

    def induced(self, h, order, power):
        """Ind_<h>^G lambda with lambda(h^k) = exp(2 pi i k power / order)."""
        Hk = {self.power(h, k): k for k in range(order)}
        assert len(Hk) == order
        reps, seen = [], set()
        for t in range(self.n):
            if t not in seen:
                reps.append(t)
                seen.update(self.mul(t, x) for x in Hk)
        m = len(reps)
        assert m * order == self.n
        inv = [self.inv(t) for t in reps]
        M = np.zeros((self.n, m, m), dtype=complex)
        for g in range(self.n):
            for i in range(m):
                for j in range(m):
                    x = self.mul(self.mul(inv[i], g), reps[j])
                    if x in Hk:
                        M[g, i, j] = np.exp(2j * np.pi * Hk[x] * power / order)
        return M


def _rep(grp, spec):
    if spec['kind'] == 'natural':
        return grp.natural()
    if spec['kind'] == 'natural2':
        P = grp.natural()
        return np.einsum('gij,gkl->gikjl', P, P).reshape(grp.n, grp.deg ** 2, grp.deg ** 2)
    if spec['kind'] == 'regular':
        return grp.regular()
    return grp.induced(spec['h'], spec['order'], spec['power'])


GROUP_REPS = [(g, i) for g, n in [('S4', 5), ('S5', 6), ('A5', 5), ('SL(2,3)', 6), ('D12', 5), ('Q8xC2', 5)]
              for i in range(n)]


@pytest.mark.parametrize('name,i', GROUP_REPS)
def test_wedderburn_of_group_representation_matches_gap_character_table(gp, name, i):
    """from_group: k_c = dimensions of the irreducible characters occurring, m_c = <chi, chi_c> (GAP), every sector
    carries exactly one irreducible character of GAP's table, and for real representations real_structure()
    gives the types R / C / H of the Frobenius-Schur indicators +1 / 0 / -1 with complex-conjugate partners."""
    g = gp[name]
    spec = g['reps'][i]
    grp = _Group(g)
    M = _rep(grp, spec)
    rng = np.random.default_rng(0)
    for a, b in rng.integers(0, grp.n, size=(20, 2)):                           # a homomorphism for GAP's law
        assert np.abs(M[a] @ M[b] - M[grp.mul(a, b)]).max() < 1e-12
    mults = spec['mults']
    dims = [int(round(row[0].real)) for row in grp.irr]
    # sanity of the GAP data against the matrices: <chi_rep, chi_irr> from the traces
    chi = np.trace(M, axis1=1, axis2=2)
    sizes = np.bincount(grp.classes)
    chi_cls = np.array([chi[grp.classes.index(c)] for c in range(len(sizes))])
    assert np.allclose(grp.irr.conj() @ (sizes * chi_cls) / grp.n, mults, atol=1e-10)

    D = decompose.from_group(torch.as_tensor(M))
    expected = sorted((dims[r], mults[r]) for r in range(len(dims)) if mults[r])
    assert sorted(zip(D.k_factors, D.multiplicities)) == expected
    assert D.commutant_dim == sum(m * m for m in mults) == decompose.commutant_dimension(torch.as_tensor(M))
    # character of every sector: tr rho_c(g), constant on classes, one row of GAP's table
    blocks = D.to_algebra(torch.as_tensor(M)).matrix.numpy()
    which = []
    for c, k in enumerate(D.k_factors):
        tr = np.trace(blocks[:, c, :k, :k], axis1=1, axis2=2)
        per_class = np.array([tr[grp.classes.index(cl)] for cl in range(len(sizes))])
        assert np.abs(tr - per_class[grp.classes]).max() < 1e-9
        match = [r for r in range(len(dims)) if np.abs(grp.irr[r] - per_class).max() < 1e-8]
        assert len(match) == 1
        which.append(match[0])
        assert D.multiplicities[c] == mults[match[0]]
    assert sorted(which) == sorted(r for r in range(len(dims)) if mults[r])
    # real structure (the representation is real iff all its matrices are)
    occurring = {r for r in range(len(dims)) if mults[r]}
    conj = {r: next(s for s in range(len(dims)) if np.abs(grp.irr[s] - grp.irr[r].conj()).max() < 1e-9)
            for r in occurring}
    closed = all(conj[r] in occurring for r in occurring)
    if np.abs(M.imag).max() == 0:
        rs = D.real_structure()
        assert rs.ok
        for c, r in enumerate(which):
            assert rs.types[c] == {1: 'R', 0: 'C', -1: 'H'}[g['indicators'][r]]
            assert which[rs.partner[c]] == conj[r]
    elif not closed:
        with pytest.raises(ValueError):
            D.real_structure()
    # the generator path (commutant by a null space) with GAP's SmallGeneratingSet gives the same structure
    if M.shape[1] <= 36:
        Dg = decompose.from_generators(torch.as_tensor(M[g['generators']]))
        assert sorted(zip(Dg.k_factors, Dg.multiplicities)) == expected


def test_live_gap_character_table_s4():
    """Live GAP: the character table of S4 from libgap equals the stored one and from_group of the natural
    representation gives the trivial and the standard (dimension 3) characters."""
    pytest.importorskip('sage.all__sagemath_gap')
    from sage.all__sagemath_gap import libgap
    tbl = libgap.CharacterTable(libgap.SymmetricGroup(4))
    irr = [[int(v) for v in libgap.ValuesOfClassFunction(chi).sage()] for chi in libgap.Irr(tbl)]
    stored = load('groups.json')['groups']['S4']
    assert irr == [[int(round(v[0])) for v in row] for row in stored['irr']]


# ==========================================================================================================
# (5) Clebsch-Gordan coefficients and Wigner 3j symbols vs SymPy
# ==========================================================================================================
def _cg_entries(two_j_max):
    for a in range(two_j_max + 1):
        for b in range(two_j_max + 1):
            for e in range(abs(a - b), a + b + 1, 2):
                for x in range(-a, a + 1, 2):
                    for y in range(-b, b + 1, 2):
                        if abs(x + y) <= e:
                            yield a, b, e, x, y


def test_clebsch_gordan_exact_matches_sympy_table():
    """clebsch_gordan_exact equals sympy.physics.wigner.clebsch_gordan exactly (sign and square) for all
    j1, j2 <= 4 (integer and half-integer), all J and m; clebsch_gordan agrees to rounding."""
    data = load('cg.json')
    table = data['table']
    count = 0
    it = {}
    for a, b, e, x, y in _cg_entries(data['two_j_max']):
        key = f"{a},{b},{e}"
        k = it.get(key, 0)
        it[key] = k + 1
        s = table[key][k]
        ref = Fraction(s.lstrip('-'))
        ref_sign = 0 if ref == 0 else (-1 if s.startswith('-') else 1)
        args = (Fraction(a, 2), Fraction(x, 2), Fraction(b, 2), Fraction(y, 2), Fraction(e, 2), Fraction(x + y, 2))
        sign, sq = clebsch_gordan_exact(*args)
        assert (sign, sq) == (ref_sign, ref), (args, s)
        assert abs(clebsch_gordan(*args) - ref_sign * math.sqrt(ref)) < 1e-15
        count += 1
    assert all(len(v) == it[k] for k, v in table.items())
    assert count == 7809


def test_live_wigner_3j_matches_sympy():
    """Live SymPy: (j1 j2 j3; m1 m2 m3) = (-1)^(j1-j2-m3) <j1 m1 j2 m2 | j3 -m3> / sqrt(2 j3 + 1) built from
    clebsch_gordan_exact equals sympy.physics.wigner.wigner_3j exactly, for all j1, j2 <= 4."""
    wigner = pytest.importorskip('sympy.physics.wigner')
    from sympy import Rational
    for a, b, e, x, y in _cg_entries(8):
        sign, sq = clebsch_gordan_exact(Fraction(a, 2), Fraction(x, 2), Fraction(b, 2), Fraction(y, 2),
                                        Fraction(e, 2), Fraction(x + y, 2))
        phase = -1 if ((a - b + x + y) // 2) % 2 else 1                     # (-1)^(j1 - j2 - m3), m3 = -(m1 + m2)
        ref = wigner.wigner_3j(Rational(a, 2), Rational(b, 2), Rational(e, 2), Rational(x, 2), Rational(y, 2),
                               Rational(-(x + y), 2))
        ref_sign = 0 if ref == 0 else (1 if ref > 0 else -1)
        assert phase * sign == ref_sign, (a, b, e, x, y)
        assert sq / (e + 1) == Fraction(str(ref ** 2)), (a, b, e, x, y)
