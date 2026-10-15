"""
Permutation-invariant qubit ensembles (torch_vn_algebra.permutational; paper, Section "Use case: a boundary time
crystal under local dissipation"; scripts/permutational.py): the Dicke-basis algebra (+)_j M_{2j+1} with the
multiplicities m_j, the local dissipators as maps between sectors, and the banded generator, against the exact
Lindblad generator on (C^2)^{(x) N} built here independently, the library's generic route (InterSectorChannel jumps
in dynamics.lindblad_rhs), closed forms (product states, pure collective decay) and qutip.piqs. CPU, double
precision.
"""
import importlib.util
import math
import runpy
from pathlib import Path

import numpy as np
import pytest
import torch

from torch_vn_algebra import dynamics
from torch_vn_algebra.permutational import PermutationalSystem, PermutationalLindbladian
from torch_vn_algebra.su2 import multiplicity

ROOT = Path(__file__).resolve().parents[1]
pytestmark = [pytest.mark.filterwarnings("ignore::torch_vn_algebra.cost.CostWarning")]


@pytest.fixture(autouse=True, scope='module')
def _single_thread_no_denormals():
    """The tensors here are small and the kernels elementwise: one thread avoids spin-waiting on shared cores;
    the far tails of pure initial states produce denormal numbers, which are flushed to zero."""
    threads = torch.get_num_threads()
    torch.set_num_threads(1)
    torch.set_flush_denormal(True)
    yield
    torch.set_flush_denormal(False)
    torch.set_num_threads(threads)

DT = torch.complex128
ALL_TERMS = dict(omega_x=0.7, omega_y=-0.4, omega_z=0.3, collective_emission=0.2, collective_pumping=0.1,
                 collective_dephasing=0.05, emission=0.3, pumping=0.2, dephasing=0.15)


# ----------------------------------------------------------------------------------------------------------------
# the full Hilbert space (C^2)^{(x) N}, independent of the module under test
# ----------------------------------------------------------------------------------------------------------------
def full_space(N):
    """Site operators (|up> = e_0), collective operators and Dicke isometries V_j of shape (2^N, 2j+1, m_j) from the
    joint eigenvectors of J^2 and J_z (highest weight) and lowering with J_-."""
    sm = torch.tensor([[0, 0], [1, 0]], dtype=DT)
    sz = torch.diag(torch.tensor([1.0, -1.0], dtype=DT))
    I2 = torch.eye(2, dtype=DT)

    def at(o, i):
        out = torch.ones(1, 1, dtype=DT)
        for k in range(N):
            out = torch.kron(out, o if k == i else I2)
        return out
    SM = [at(sm, i) for i in range(N)]
    SP = [x.conj().T.contiguous() for x in SM]
    SZ = [at(sz, i) for i in range(N)]
    Jm, Jp, Jz = sum(SM), sum(SP), sum(SZ) / 2
    J2 = Jp @ Jm + Jz @ Jz - Jz
    V = {}
    for tj in range(N % 2, N + 1, 2):
        j = tj / 2
        idx = torch.nonzero((torch.diagonal(Jz).real - j).abs() < 1e-9).flatten()
        w, U = torch.linalg.eigh(J2[idx][:, idx])
        sel = (w - j * (j + 1)).abs() < 1e-8
        v = torch.zeros(2 ** N, int(sel.sum()), dtype=DT)
        v[idx] = U[:, sel]
        cols = [v]
        for a in range(tj):
            m = j - a
            v = Jm @ v / math.sqrt(j * (j + 1) - m * (m - 1))
            cols.append(v)
        V[tj] = torch.stack(cols, 1)
    return dict(SM=SM, SP=SP, SZ=SZ, Jm=Jm, Jp=Jp, Jz=Jz, V=V)


def embed(ps, sig, F):
    """(+)_j sigma_j (x) 1_{m_j} / m_j as an operator on (C^2)^{(x) N}."""
    rho = torch.zeros(2 ** ps.N, 2 ** ps.N, dtype=DT)
    for c, tj in enumerate(ps.two_j):
        Vj, k = F['V'][tj], tj + 1
        rho = rho + torch.einsum('xia,ij,yja->xy', Vj, sig[c, :k, :k].to(DT), Vj.conj()) / Vj.shape[-1]
    return rho


def project(ps, rho, F):
    """sigma_j = tr_{C^{m_j}} (block j of rho)."""
    out = torch.zeros(ps.C, ps.k_max, ps.k_max, dtype=DT)
    for c, tj in enumerate(ps.two_j):
        Vj, k = F['V'][tj], tj + 1
        out[c, :k, :k] = torch.einsum('xia,xy,yja->ij', Vj.conj(), rho, Vj)
    return out


def full_rhs(F, p, rho):
    H = p['omega_x'] * (F['Jp'] + F['Jm']) / 2 + p['omega_y'] * (F['Jp'] - F['Jm']) / 2j + p['omega_z'] * F['Jz']

    def D(L, r):
        Ld = L.conj().T
        return L @ r @ Ld - 0.5 * (Ld @ L @ r + r @ Ld @ L)
    out = -1j * (H @ rho - rho @ H)
    out = out + p['collective_emission'] * D(F['Jm'], rho) + p['collective_pumping'] * D(F['Jp'], rho) \
        + p['collective_dephasing'] * D(F['Jz'], rho)
    for i in range(len(F['SM'])):
        out = out + p['emission'] * D(F['SM'][i], rho) + p['pumping'] * D(F['SP'][i], rho) \
            + p['dephasing'] * D(F['SZ'][i], rho)
    return out


def random_state(ps, batch, seed):
    g = torch.Generator().manual_seed(seed)
    X = torch.randn(batch, ps.C, ps.k_max, ps.k_max, dtype=DT, generator=g)
    X = X @ X.conj().transpose(-2, -1) * ps.algebra._const('block_mask', torch.bool)
    return X / torch.diagonal(X, dim1=-2, dim2=-1).sum((-2, -1)).real[:, None, None, None]


# ----------------------------------------------------------------------------------------------------------------
def test_algebra_dimensions_and_trace_weights():
    """sum_j (2j+1) m_j = 2^N, Tr_H(rho) = sum_j tr sigma_j for the embedding rho = (+)_j sigma_j (x) 1/m_j."""
    for N in (1, 2, 5, 8, 31):
        ps = PermutationalSystem(N)
        assert sum((t + 1) * m for t, m in zip(ps.two_j, ps.multiplicities)) == 2 ** N
        assert ps.k_max == N + 1 and ps.algebra.charges == ps.two_j
    ps, F = PermutationalSystem(5), full_space(5)
    X = random_state(ps, 1, 0)[0]
    rho = embed(ps, X, F)
    assert abs(torch.trace(rho).real.item() - 1) < 1e-13
    assert (project(ps, rho, F) - X).abs().max() < 1e-13


def test_collective_operators_match_full_space():
    """J_+, J_-, J_x, J_y, J_z of the module are the restrictions of the collective operators to the blocks."""
    N = 5
    ps, F = PermutationalSystem(N), full_space(N)
    pairs = ((ps.jplus(), F['Jp']), (ps.jminus(), F['Jm']), (ps.jz(), F['Jz']),
             (ps.jx(), (F['Jp'] + F['Jm']) / 2), (ps.jy(), (F['Jp'] - F['Jm']) / 2j))
    for op, full in pairs:
        for c, tj in enumerate(ps.two_j):
            Vj = F['V'][tj][:, :, 0]                                  # one copy of V_j
            assert (Vj.conj().T @ full @ Vj - op.matrix[0, c, :tj + 1, :tj + 1]).abs().max() < 1e-12


def test_loss_part_is_sum_of_local_terms():
    """Phi^*(1) = sum K^* K equals sum_i L_i^* L_i: N/2 + J_z, N/2 - J_z and N for sigma_-, sigma_+, sigma_z
    (also with the truncation j_min, where the jumps out of the kept sectors still count)."""
    for N, j_min in ((1, None), (4, None), (7, None), (12, None), (12, 3), (13, 2.5)):
        ps = PermutationalSystem(N, j_min=j_min)
        m = ps.m_values()
        mask = torch.arange(ps.k_max)[None, :] <= torch.tensor(ps.two_j)[:, None]
        for kind, ref in (('-', N / 2 + m), ('+', N / 2 - m), ('z', N + 0 * m)):
            assert ((ps.local_loss_diagonal(kind) - ref) * mask).abs().max() < 1e-12


@pytest.mark.parametrize('N', [2, 3, 4, 5, 6])
def test_banded_generator_equals_full_space_lindbladian(N):
    """The banded generator is the exact GKSL generator on (C^2)^{(x) N} restricted to M (all nine terms)."""
    ps, F = PermutationalSystem(N), full_space(N)
    X = random_state(ps, 2, N)
    ours = PermutationalLindbladian(ps, **ALL_TERMS).rhs(0.0, X)
    for b in range(2):
        full = full_rhs(F, ALL_TERMS, embed(ps, X[b], F))
        assert (embed(ps, project(ps, full, F), F) - full).abs().max() < 1e-12      # the generator maps M into M
        assert (project(ps, full, F) - ours[b]).abs().max() < 1e-12


@pytest.mark.parametrize('kind', ['-', '+', 'z'])
def test_local_jump_channel_is_sum_of_site_conjugations(kind):
    """local_jump(kind) applied to sigma is tr_mult of sum_i L_i rho L_i^* (the inter-sector Kraus form)."""
    N = 5
    ps, F = PermutationalSystem(N), full_space(N)
    L = {'-': F['SM'], '+': F['SP'], 'z': F['SZ']}[kind]
    X = random_state(ps, 1, 7)
    phi = ps.local_jump(kind)
    out = phi(ps.algebra.operator(X)).matrix[0]
    rho = embed(ps, X[0], F)
    full = sum(Li @ rho @ Li.conj().T for Li in L)
    assert (project(ps, full, F) - out).abs().max() < 1e-12


def test_generic_library_route_equals_banded_batched():
    """dynamics.lindblad_rhs with InterSectorChannel jumps (rates folded in, batch of 2 parameter sets) equals rhs."""
    N = 5
    ps = PermutationalSystem(N)
    p = dict(ALL_TERMS, omega_y=[-0.4, 0.9], emission=[0.3, 0.0], dephasing=[0.15, 0.6])
    L = PermutationalLindbladian(ps, **p)
    H, jumps, rates = L.dynamics_arguments()
    f = dynamics.lindblad_rhs(H, jumps, rates, algebra=ps.algebra)
    X = random_state(ps, 2, 3)
    assert (f(0.0, X) - L.rhs(0.0, X)).abs().max() < 1e-12
    traj_lib = dynamics.lindblad_evolve(ps.coherent_state(1.0, 0.3), H, jumps, [0.0, 0.5], rates, substeps=20)
    traj = L.evolve(ps.coherent_state(1.0, 0.3), [0.0, 0.5], substeps=20)
    assert (traj_lib[-1].matrix - traj[-1].matrix).abs().max() < 1e-12


def test_full_space_trajectory():
    """RK4 on (C^2)^{(x) 4} and the banded RK4 give the same blocks sigma_j(t) (same steps, round-off only)."""
    N = 4
    ps, F = PermutationalSystem(N), full_space(N)
    L = PermutationalLindbladian(ps, **ALL_TERMS)
    rho0 = ps.coherent_state(2.0, 0.5)
    ours = L.evolve(rho0, [0.0, 1.0], substeps=50)[-1].matrix[0]
    r, h = embed(ps, rho0.matrix[0], F), 1.0 / 50
    for _ in range(50):
        k1 = full_rhs(F, ALL_TERMS, r)
        k2 = full_rhs(F, ALL_TERMS, r + h / 2 * k1)
        k3 = full_rhs(F, ALL_TERMS, r + h / 2 * k2)
        k4 = full_rhs(F, ALL_TERMS, r + h * k3)
        r = r + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    assert (project(ps, r, F) - ours).abs().max() < 1e-13


def test_real_algebra_with_drive_about_y():
    """With omega_x = omega_z = 0 the generator is real: a real algebra gives the complex result."""
    N = 9
    pr = dict(omega_y=1.3, collective_emission=0.2, emission=[0.1, 0.4], dephasing=0.05, pumping=0.02)
    ps_c, ps_r = PermutationalSystem(N), PermutationalSystem(N, complex_valued=False)
    a = PermutationalLindbladian(ps_c, **pr).evolve(ps_c.coherent_state(math.pi), [0.0, 0.7], substeps=30)[-1]
    b = PermutationalLindbladian(ps_r, **pr).evolve(ps_r.coherent_state(math.pi), [0.0, 0.7], substeps=30)[-1]
    assert b.matrix.dtype == torch.float64
    assert (a.matrix.imag.abs().max() < 1e-13) and (a.matrix.real - b.matrix).abs().max() < 1e-13
    with pytest.raises(ValueError):
        PermutationalLindbladian(ps_r, omega_x=1.0)


def test_coherent_state_moments():
    """<J> of the coherent spin state is N/2 (sin t cos p, sin t sin p, cos t); <J^2> = N/2 (N/2 + 1)."""
    N, th, ph = 11, 1.1, 0.4
    ps = PermutationalSystem(N)
    rho = ps.coherent_state(th, ph)
    ex = [rho.expectation(op).real.item() for op in (ps.jx(), ps.jy(), ps.jz(), ps.j2())]
    ref = [N / 2 * math.sin(th) * math.cos(ph), N / 2 * math.sin(th) * math.sin(ph), N / 2 * math.cos(th),
           N / 2 * (N / 2 + 1)]
    assert np.allclose(ex, ref, atol=1e-12)


def test_product_states_stay_products_under_local_dissipation():
    """Local emission, pumping and dephasing alone (N = 30): <J_z>, <J_x> follow the single-qubit Bloch equations,
    <J^2> = 3N/4 + N(N-1)|s|^2/4, and p_j = m_j sum_m lambda_+^{N/2+m} lambda_-^{N/2-m}."""
    N, ge, gp, gd, th, t = 30, 0.3, 0.1, 0.2, 2 * math.pi / 3, 0.8
    ps = PermutationalSystem(N, complex_valued=False)
    L = PermutationalLindbladian(ps, emission=ge, pumping=gp, dephasing=gd)
    rho = L.evolve(ps.coherent_state(th), [0.0, t], substeps=200)[-1]
    g = ge + gp
    sz = (gp - ge) / g + (math.cos(th) - (gp - ge) / g) * math.exp(-g * t)
    sx = math.sin(th) * math.exp(-(g / 2 + 2 * gd) * t)
    assert abs(rho.expectation(ps.jz()).item() - N / 2 * sz) < 1e-10
    assert abs(rho.expectation(ps.jx()).item() - N / 2 * sx) < 1e-10
    s = math.hypot(sx, sz)
    assert abs(rho.expectation(ps.j2()).item() - (3 * N / 4 + N * (N - 1) * s * s / 4)) < 1e-8
    lp, lm = (1 + s) / 2, (1 - s) / 2
    p_ref = [multiplicity(N, tj / 2) * sum(lp ** ((N + tm) / 2) * lm ** ((N - tm) / 2) for tm in range(-tj, tj + 1, 2))
             for tj in ps.two_j]
    assert np.abs(ps.sector_probabilities(rho)[0].numpy() - np.array(p_ref)).sum() < 1e-8


def test_pure_collective_decay_mean_first_passage_time():
    """Pure collective decay from |N/2, N/2> (N = 12, sector j = N/2 only): int_0^inf P(m > 0) dt = H_N/(N+1)."""
    N = 12
    ps = PermutationalSystem(N, j_min=N / 2)
    assert ps.C == 1
    L = PermutationalLindbladian(ps, collective_emission=1.0)
    T, n = 40.0 / N, 2000
    ts = torch.linspace(0, T, n + 1, dtype=torch.float64)
    traj = dynamics.rk4(L.rhs, ps.dicke_state(N / 2, N / 2).matrix, ts, substeps=2)
    P = torch.diagonal(traj[:, 0, 0], dim1=-2, dim2=-1)[:, :N // 2].sum(-1).real    # m = N/2 .. 1
    integral = float(((P[1:] + P[:-1]) / 2).sum()) * T / n
    exact = sum(1.0 / k for k in range(1, N + 1)) / (N + 1)
    assert abs(integral / exact - 1) < 1e-6


def test_truncation_is_a_compression_with_error_bound():
    """j_min truncation: the kept blocks satisfy 0 <= sigma_trunc <= sigma (as blocks of the full evolution), and
    |<J_z>_trunc - <J_z>| on the kept sectors is at most ||J_z|| (2 (1 - Tr sigma_trunc))."""
    N, j_min = 10, 3
    pr = dict(omega_y=2.0, collective_emission=0.2, emission=0.3, dephasing=0.1)
    full, trunc = PermutationalSystem(N, complex_valued=False), PermutationalSystem(N, complex_valued=False,
                                                                                    j_min=j_min)
    a = PermutationalLindbladian(full, **pr).evolve(full.coherent_state(math.pi), [0.0, 3.0], substeps=300)[-1]
    b = PermutationalLindbladian(trunc, **pr).evolve(trunc.coherent_state(math.pi), [0.0, 3.0], substeps=300)[-1]
    off = full.C - trunc.C
    lost = 1 - b.matrix.diagonal(dim1=-2, dim2=-1).sum().item()
    assert 0 < lost < 1
    for c in range(trunc.C):
        k = trunc.two_j[c] + 1
        diff = a.matrix[0, off + c, :k, :k] - b.matrix[0, c, :k, :k]
        assert torch.linalg.eigvalsh(diff).min() > -1e-12
    jz_full = float((a.matrix[0, off:].diagonal(dim1=-2, dim2=-1) * trunc.m_values()).sum())
    jz_trunc = float((b.matrix[0].diagonal(dim1=-2, dim2=-1) * trunc.m_values()).sum())
    assert abs(jz_full - jz_trunc) <= N / 2 * 2 * lost


def test_piqs_generator_and_steady_state():
    """qutip.piqs (Dicke basis, gamma_D = 4 x dephasing): the same generator (PIQS rates are single precision)
    and the same steady state for a driven ensemble with collective and local decay."""
    qutip = pytest.importorskip('qutip')
    pytest.importorskip('qutip.piqs')
    from qutip import piqs                     # the module with Dicke and jspin (qutip 4 and 5)
    N = 6
    ps = PermutationalSystem(N)
    p = ALL_TERMS
    d = piqs.Dicke(N, emission=p['emission'], pumping=p['pumping'], dephasing=4 * p['dephasing'],
                   collective_emission=p['collective_emission'], collective_pumping=p['collective_pumping'],
                   collective_dephasing=p['collective_dephasing'])
    jx, jy, jz = piqs.jspin(N)
    d.hamiltonian = p['omega_x'] * jx + p['omega_y'] * jy + p['omega_z'] * jz
    Lq = d.liouvillian()
    from scipy.linalg import block_diag

    def to_piqs(X):
        return block_diag(*[X[c, :t + 1, :t + 1].numpy() for c, t in reversed(list(enumerate(ps.two_j)))])

    def from_piqs(M):
        out, o = torch.zeros(ps.C, ps.k_max, ps.k_max, dtype=DT), 0
        for c in reversed(range(ps.C)):
            k = ps.two_j[c] + 1
            out[c, :k, :k] = torch.as_tensor(np.asarray(M)[o:o + k, o:o + k])
            o += k
        return out
    X = random_state(ps, 1, 11)[0]
    rq = qutip.vector_to_operator(Lq * qutip.operator_to_vector(qutip.Qobj(to_piqs(X)))).full()
    ours = PermutationalLindbladian(ps, **p).rhs(0.0, X[None])[0]
    assert (from_piqs(rq) - ours).abs().max() / ours.abs().max() < 1e-5
    # steady state: long evolution here against the null vector of the PIQS Liouvillian on the block entries
    pr = dict(omega_y=2.0, collective_emission=1 / 3, emission=0.5, dephasing=0.1)
    psr = PermutationalSystem(N, complex_valued=False)
    ss = PermutationalLindbladian(psr, **pr).evolve(psr.coherent_state(math.pi), [0.0, 60.0], substeps=6000)[-1]
    d = piqs.Dicke(N, emission=pr['emission'], dephasing=4 * pr['dephasing'],
                   collective_emission=pr['collective_emission'])
    d.hamiltonian = pr['omega_y'] * jy
    A = d.liouvillian().full()
    nds = sum(t + 1 for t in ps.two_j)
    idx, o = [], 0
    for c in reversed(range(ps.C)):
        k = ps.two_j[c] + 1
        idx += [(o + cc) * nds + (o + r) for cc in range(k) for r in range(k)]
        o += k
    _, _, Vh = np.linalg.svd(A[np.ix_(idx, idx)])
    M = np.zeros((nds, nds), dtype=complex)
    M[np.array(idx) % nds, np.array(idx) // nds] = Vh[-1].conj()
    M /= np.trace(M)
    assert (from_piqs(M) - ss.matrix[0].to(DT)).abs().max() < 1e-5


def test_listing_runs():
    """The listing man/article/listings/permutational.py runs: m_z(0) = -1, sum_j p_j = 1, without local emission
    only j = N/2 is populated, and the mean total spin decreases with the local rate."""
    ns = runpy.run_path(str(ROOT / 'man' / 'article' / 'listings' / 'permutational.py'))
    m_z, p_j, ps = ns['m_z'], ns['p_j'], ns['ps']
    assert torch.allclose(m_z[0], -torch.ones(4, dtype=torch.float64))
    assert torch.allclose(p_j.sum(-1), torch.ones(4, dtype=torch.float64), atol=1e-12)
    assert abs(p_j[0, -1].item() - 1) < 1e-12
    j = torch.tensor(ps.two_j, dtype=torch.float64) / 2
    mean_j = (p_j * j).sum(-1)
    assert bool((mean_j[1:] < mean_j[:-1]).all())


def test_lifetime_and_mean_field_helpers_of_the_script():
    """scripts/permutational.py: for -exp(-t/5) cos(2t) the half-swings decay as exp(-t/5) from the first swing
    midpoint t ~ pi/4, so tau ~ 5 + pi/4; the mean-field m_z without local terms stays on its orbit (no decay)."""
    spec = importlib.util.spec_from_file_location('permutational_script', ROOT / 'scripts' / 'permutational.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    t = np.arange(0, 40, 0.01)
    tau = mod.lifetime(t, -np.exp(-t / 5) * np.cos(2 * t))
    assert abs(tau - (5 + math.pi / 4)) < 0.06            # extrema of the damped cosine are shifted slightly
    pytest.importorskip('scipy')
    ts, z = mod.mean_field(0.0, 0.0, 40.0)
    assert math.isinf(mod.lifetime(ts, z)) and abs(z.max() - 1) < 1e-3
