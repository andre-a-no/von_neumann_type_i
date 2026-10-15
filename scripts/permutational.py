#!/usr/bin/env python3
"""
Section "Use case: a boundary time crystal under local dissipation" of the paper.

N two-level systems with collective drive and decay (the boundary time crystal of Iemini et al., PRL 121, 035301
(2018)) and, in addition, local emission sum_i D[sigma_-^i] or local dephasing sum_i D[sigma_z^i]:

    d rho/dt = -i [omega J_y, rho] + (kappa / S) D[J_-] rho
               + gamma sum_i D[sigma_-^i] rho + gamma_phi sum_i D[sigma_z^i] rho,

S = N/2, omega = 2 kappa (the time-crystal regime omega > kappa), initial state all spins down. (The drive about y
is the drive about x of Iemini et al. rotated by exp(-i pi/2 J_z), which commutes with every dissipator; with it
the generator is real and the computation runs in a real algebra.) The state stays permutation invariant, i.e. in
M = (+)_j M_{2j+1} (x) 1_{m_j}; the collective terms act inside the blocks, the local ones couple j with j and
j +- 1 (torch_vn_algebra.permutational). Observable: m_z(t) = <J_z>/S. Lifetime tau: the first time at which the
half-swing |m_z(t_k+1) - m_z(t_k)|/2 between consecutive extrema t_k falls below 1/e of the first half-swing
(log-linear interpolation between swing midpoints; censored at the end of the run).

Parts:
  generator   the banded generator against (i) the exact Lindblad generator on (C^2)^{(x) N}, N <= 8, applied to
              random permutation-invariant states, (ii) the library's generic route (dynamics.lindblad_rhs with
              the local jumps as InterSectorChannel), (iii) a full-space trajectory, and (iv) qutip.piqs (generator,
              trajectories, steady states);
  analytic    (a) pure collective decay from |N/2, N/2>: int_0^infty P(m > 0) dt = H_N / (Gamma (N + 1)) exactly
              (sum of the inverse rates of a pure-death chain), burst delay t_D Gamma N / ln N and peak intensity
              I_max / (Gamma N^2/4) (mean field: 1 and 1); (b) local emission, pumping and dephasing alone keep
              product states product: <J_z>, <J_x>, <J^2> and the full spin distribution
              p_j = m_j sum_m lambda_+^{N/2+m} lambda_-^{N/2-m} in closed form, at large N;
  sweep       tau(N, gamma) for local emission and tau(N, gamma_phi) for local dephasing, one batched run per N
              (all rates at once; finished batch elements are dropped), gamma = 0 in the single sector j = N/2;
              the mean-field limit tau_MF(gamma) (N -> infinity) and the additive rule
              1/tau ~ 1/tau(N, 0) + 1/tau_MF(gamma) as an observation.

Writes generated/permutational.{json,tex} and figures/permutational.pdf.

    OMP_NUM_THREADS=1 python scripts/permutational.py --mode quick  # under a minute on one CPU thread
    OMP_NUM_THREADS=1 python scripts/permutational.py --mode check --out man/article/generated \\
        --figdir man/article/figures                                 # N <= 80, about 20 min on one CPU thread
    python scripts/permutational.py --mode full --device cuda       # N <= 200, designed for a GPU (not run)

PM_FROM_JSON=<path to permutational.json> re-renders the .tex and the figure from a saved run without computing.
On a loaded machine set OMP_NUM_THREADS=1: the elementwise kernels are memory bound and spin-waiting threads
slow them down by orders of magnitude when the cores are shared.
"""
import math
import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
for _p in (str(ROOT), str(ROOT / 'scripts')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from torch_vn_algebra import cost, dynamics                                               # noqa: E402
from torch_vn_algebra.permutational import PermutationalSystem, PermutationalLindbladian  # noqa: E402

OMEGA, KAPPA = 2.0, 1.0
DT_OUT = 0.05                       # sampling of m_z(t)
FIT_NMIN, RULE_NMIN = 20, 30        # fit of tau(N, 0) = a N + b for N >= FIT_NMIN; rule summary for N >= RULE_NMIN
SAFETY = 3.5                        # RK4 step h = SAFETY / (crude upper bound on the spectral radius)

CONFIGS = {
    'quick': dict(Ns=[8, 12, 16, 24], T=30.0, ge=[0.05, 0.1, 0.2, 0.5], gd=[0.0125, 0.025, 0.05, 0.125],
                  full_N=[2, 3, 4, 5, 6], traj_N=5, lib_N=6, piqs_N=[6, 10], piqs_ss_N=[6], prod_N=40,
                  sr_N=[10, 20, 40, 80], T_mf=600.0, conv_N=12),
    'check': dict(Ns=[10, 20, 30, 40, 50, 60, 80], T=45.0, ge=[0.02, 0.05, 0.1, 0.2, 0.5],
                  gd=[0.005, 0.0125, 0.025, 0.05, 0.125], full_N=[2, 3, 4, 5, 6, 7, 8], traj_N=7, lib_N=10,
                  piqs_N=[10, 20, 30], piqs_ss_N=[10, 20], prod_N=100, sr_N=[10, 20, 50, 100, 200], T_mf=600.0,
                  conv_N=30),
    'full': dict(Ns=[10, 20, 30, 40, 60, 80, 100, 120, 160, 200], T=60.0, ge=[0.02, 0.05, 0.1, 0.2, 0.5],
                 gd=[0.005, 0.0125, 0.025, 0.05, 0.125], full_N=[2, 3, 4, 5, 6, 7, 8], traj_N=8, lib_N=12,
                 piqs_N=[10, 20, 30, 40], piqs_ss_N=[10, 20, 30], prod_N=200, sr_N=[10, 20, 50, 100, 200, 400],
                 T_mf=600.0, conv_N=40),
}


# ----------------------------------------------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------------------------------------------
def lifetime(t, z):
    """First time at which the half-swing between consecutive extrema of z falls below 1/e of the first one
    (log-linear interpolation between the swing midpoints); inf if it does not happen before the end of t. The
    end point counts as an extremum only if no extremum occurred for twice the longest gap between extrema."""
    t, z = np.asarray(t), np.asarray(z)
    d = np.diff(z)
    ext = [0] + [i for i in range(1, len(z) - 1) if d[i - 1] * d[i] < 0]
    gap = max(np.diff(ext)) if len(ext) > 1 else 0
    if len(z) - 1 - ext[-1] > 2 * gap:
        ext.append(len(z) - 1)
    if len(ext) < 2:
        return math.inf
    a = np.abs(np.diff(z[ext])) / 2
    tm = (t[ext][1:] + t[ext][:-1]) / 2
    thr = a[0] / math.e
    for i in range(1, len(a)):
        if a[i] < thr:
            l0, l1 = math.log(a[i - 1]), math.log(max(a[i], 1e-300))
            return float(tm[i - 1] + (math.log(thr) - l0) / (l1 - l0) * (tm[i] - tm[i - 1]))
    return math.inf


def mean_field(gamma, gamma_phi, T, dt=DT_OUT, omega=OMEGA, kappa=KAPPA):
    """m_z(t) of the mean-field equations (N -> infinity) for m = <J>/S, drive omega J_y, collective decay
    (kappa/S) D[J_-], local emission gamma and dephasing gamma_phi, from m = (0, 0, -1)."""
    from scipy.integrate import solve_ivp

    def f(t, m):
        x, y, z = m
        tr = gamma / 2 + 2 * gamma_phi
        return [omega * z + kappa * x * z - tr * x, kappa * y * z - tr * y,
                -omega * x - kappa * (x * x + y * y) - gamma * (1 + z)]
    ts = np.arange(0, T + 1e-9, dt)
    s = solve_ivp(f, (0, T), [0.0, 0.0, -1.0], t_eval=ts, rtol=1e-11, atol=1e-13, method='DOP853')
    return ts, s.y[2]


def step_size(L, dt, safety=SAFETY):
    """Substeps per output interval dt for RK4 with h <= safety / (spectral radius bound)."""
    rad = float(L.spectral_radius_bound().max())
    return max(1, int(math.ceil(dt * rad / safety)))


def mz_of(ps, y):
    return ((torch.diagonal(y, dim1=-2, dim2=-1).real * ps.m_values().to(y.device, y.real.dtype)).sum((-2, -1))
            / (ps.N / 2)).cpu().numpy()


# ----------------------------------------------------------------------------------------------------------------
# the sweep
# ----------------------------------------------------------------------------------------------------------------
def btc_batch(N, ge, gd, T, dev, log, j_min=None, early_stop=True, safety=SAFETY):
    """One batched run: element b has local emission ge[b] and local dephasing gd[b]. Returns times, m_z (T, B)
    (nan after an element was dropped), tau (B,), seconds, final spin distributions of the elements still running."""
    ps = PermutationalSystem(N, complex_valued=False, device=dev, j_min=j_min)
    ge, gd = np.asarray(ge, float), np.asarray(gd, float)
    B = len(ge)
    rho = ps.coherent_state(math.pi).matrix.to(dev).expand(B, -1, -1, -1).clone()
    active = np.arange(B)
    nT = int(round(T / DT_OUT)) + 1
    ts = np.arange(nT) * DT_OUT
    mz = np.full((nT, B), np.nan)
    mz[0] = -1.0
    tau = np.full(B, math.inf)
    t0 = time.time()
    L, sub = None, None
    for i in range(1, nT):
        if L is None:
            L = PermutationalLindbladian(ps, omega_y=OMEGA, collective_emission=KAPPA / (N / 2),
                                         emission=torch.tensor(ge[active]), dephasing=torch.tensor(gd[active]))
            sub = step_size(L, DT_OUT, safety)
        rho = dynamics.rk4(L.rhs, rho, [ts[i - 1], ts[i]], substeps=sub)[-1]
        mz[i, active] = mz_of(ps, rho)
        if early_stop and i % 20 == 0:                           # drop the elements whose lifetime is settled
            keep = []
            for b_pos, b in enumerate(active):
                tb = lifetime(ts[:i + 1], mz[:i + 1, b])
                tau[b] = tb
                if not (math.isfinite(tb) and ts[i] > tb + 6.0):
                    keep.append(b_pos)
            if len(keep) < len(active):
                active = active[keep]
                if len(active) == 0:
                    break
                rho = rho[keep].contiguous()
                L = None
    for b in range(B):
        tau[b] = lifetime(ts[~np.isnan(mz[:, b])], mz[~np.isnan(mz[:, b]), b])
    secs = time.time() - t0
    return dict(times=ts, mz=mz, tau=tau, seconds=secs, substeps=sub)


def part_sweep(cfg, dev, log):
    out = dict(N=[], tau0=[], tau_e=[], tau_d=[], seconds=[], seconds0=[], traces={})
    ge, gd = cfg['ge'], cfg['gd']
    for N in cfg['Ns']:
        T0 = 0.6 * N + 15.0                                     # gamma = 0: the single sector j = N/2
        r0 = btc_batch(N, [0.0], [0.0], T0, dev, log, j_min=N / 2)
        ge_all = list(ge) + [0.0] * len(gd)
        gd_all = [0.0] * len(ge) + list(gd)
        r = btc_batch(N, ge_all, gd_all, cfg['T'], dev, log)
        out['N'].append(N)
        out['tau0'].append(float(r0['tau'][0]))
        out['tau_e'].append([float(x) for x in r['tau'][:len(ge)]])
        out['tau_d'].append([float(x) for x in r['tau'][len(ge):]])
        out['seconds'].append(r['seconds'])
        out['seconds0'].append(r0['seconds'])
        out['traces'][N] = dict(t=r['times'][::4].tolist(),
                                mz=[[None if np.isnan(v) else float(v) for v in col] for col in r['mz'][::4].T],
                                t0=r0['times'][::4].tolist(), mz0=r0['mz'][::4, 0].tolist())
        log(f"  N={N:4d}  tau(0)={r0['tau'][0]:7.2f}  emission " +
            ' '.join(f'{x:6.2f}' for x in r['tau'][:len(ge)]) + '  dephasing ' +
            ' '.join(f'{x:6.2f}' for x in r['tau'][len(ge):]) + f"   [{r['seconds']:.0f} s, {r['substeps']} sub]")
    # time-step check: the same batch without early stopping at SAFETY and SAFETY / 2
    Nc = cfg['conv_N']
    ge_all = list(ge) + [0.0] * len(gd)
    gd_all = [0.0] * len(ge) + list(gd)
    Tc = min(cfg['T'], 30.0)
    ra = btc_batch(Nc, ge_all, gd_all, Tc, dev, log, early_stop=False)
    rb = btc_batch(Nc, ge_all, gd_all, Tc, dev, log, early_stop=False, safety=SAFETY / 2)
    fin = np.isfinite(ra['tau']) & np.isfinite(rb['tau'])
    out['conv'] = dict(N=Nc, T=Tc, dmz=float(np.abs(ra['mz'] - rb['mz']).max()),
                       dtau=float(np.abs(ra['tau'][fin] - rb['tau'][fin]).max()) if fin.any() else 0.0,
                       substeps=[ra['substeps'], rb['substeps']])
    log(f"  step check N={Nc}: max |m_z(h) - m_z(h/2)| = {out['conv']['dmz']:.1e}, "
        f"max |tau(h) - tau(h/2)| = {out['conv']['dtau']:.1e}")
    last = out['N'][-1]                                     # keep the traces of the largest N only (figure)
    out['traces'] = {last: {k: [None if v is None else (round(v, 6) if isinstance(v, float) else
                                                         [None if w is None else round(w, 6) for w in v])
                                for v in vals] for k, vals in out['traces'][last].items()}}
    mf_e, mf_d = [], []
    for g in ge:
        t, z = mean_field(g, 0.0, cfg['T_mf'])
        mf_e.append(lifetime(t, z))
    for g in gd:
        t, z = mean_field(0.0, g, cfg['T_mf'])
        mf_d.append(lifetime(t, z))
    out['tau_mf_e'], out['tau_mf_d'] = mf_e, mf_d
    t, z = mean_field(0.0, 0.0, 100.0)
    out['mf_traces'] = {str(g): [round(v, 6) for v in mean_field(g, 0.0, cfg['T'])[1][::4].tolist()] for g in ge}
    log("  mean field: emission " + ' '.join(f'{x:6.2f}' for x in mf_e) + "  dephasing " +
        ' '.join(f'{x:6.2f}' for x in mf_d))
    # linear growth at gamma = 0 and the additive rule
    Ns = np.array(out['N'], float)
    tau0 = np.array(out['tau0'])
    sel = Ns >= FIT_NMIN
    A = np.stack([Ns[sel], np.ones(sel.sum())], 1)
    slope, icpt = np.linalg.lstsq(A, tau0[sel], rcond=None)[0]
    out['slope'], out['intercept'] = float(slope), float(icpt)
    dev_rule = []
    for kind, mf in (('tau_e', mf_e), ('tau_d', mf_d)):
        for iN, N in enumerate(out['N']):
            for ig, tm in enumerate(mf):
                tau = out[kind][iN][ig]
                if math.isfinite(tau) and math.isfinite(tm):
                    rule = 1.0 / (1.0 / out['tau0'][iN] + 1.0 / tm)
                    dev_rule.append(dict(kind=kind, N=N, g=(ge if kind == 'tau_e' else gd)[ig], tau=tau,
                                         rule=rule, rel=tau / rule - 1))
    out['rule'] = dev_rule
    rel = [d['rel'] for d in dev_rule]
    big = [d['rel'] for d in dev_rule if d['N'] >= RULE_NMIN]
    log(f"  tau(N,0) = {slope:.4f} N + {icpt:.3f};  tau/rule - 1 in [{min(rel):+.3f}, {max(rel):+.3f}] "
        f"(N >= {RULE_NMIN}: [{min(big) if big else 0:+.3f}, {max(big) if big else 0:+.3f}])")
    return out


# ----------------------------------------------------------------------------------------------------------------
# validation: full Hilbert space, generic library route, PIQS
# ----------------------------------------------------------------------------------------------------------------
def _full_space(N):
    """Site operators and the Dicke isometries V_j (2^N, 2j+1, m_j) on (C^2)^{(x) N} (|up> = e_0), built from J^2 and
    J_z by eigendecomposition and lowering; independent of torch_vn_algebra.permutational."""
    dt = torch.complex128
    sm = torch.tensor([[0, 0], [1, 0]], dtype=dt)
    sz = torch.diag(torch.tensor([1.0, -1.0], dtype=dt))
    I2 = torch.eye(2, dtype=dt)

    def at(o, i):
        out = torch.ones(1, 1, dtype=dt)
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
        v = torch.zeros(2 ** N, int(sel.sum()), dtype=dt)
        v[idx] = U[:, sel]
        cols = [v]
        for a in range(tj):
            m = j - a
            v = Jm @ v / math.sqrt(j * (j + 1) - m * (m - 1))
            cols.append(v)
        V[tj] = torch.stack(cols, 1)
    return dict(SM=SM, SP=SP, SZ=SZ, Jm=Jm, Jp=Jp, Jz=Jz, V=V)


def _embed(ps, sig, F):
    rho = 0
    for c, tj in enumerate(ps.two_j):
        Vj, k = F['V'][tj], tj + 1
        rho = rho + torch.einsum('xia,ij,yja->xy', Vj, sig[c, :k, :k], Vj.conj()) / Vj.shape[-1]
    return rho


def _project(ps, rho, F):
    out = torch.zeros(ps.C, ps.k_max, ps.k_max, dtype=torch.complex128)
    for c, tj in enumerate(ps.two_j):
        Vj, k = F['V'][tj], tj + 1
        out[c, :k, :k] = torch.einsum('xia,xy,yja->ij', Vj.conj(), rho, Vj)
    return out


def _full_rhs(F, p, rho):
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


ALL_TERMS = dict(omega_x=0.7, omega_y=-0.4, omega_z=0.3, collective_emission=0.2, collective_pumping=0.1,
                 collective_dephasing=0.05, emission=0.3, pumping=0.2, dephasing=0.15)


def part_generator(cfg, log):
    torch.manual_seed(1)
    res = dict(full=[], lib=None, traj=None)
    for N in cfg['full_N']:
        ps = PermutationalSystem(N)
        F = _full_space(N)
        X = torch.randn(4, ps.C, ps.k_max, ps.k_max, dtype=torch.complex128)
        X = X @ X.conj().transpose(-2, -1) * ps.algebra._const('block_mask', torch.bool)
        X = X / torch.diagonal(X, dim1=-2, dim2=-1).sum((-2, -1)).real[:, None, None, None]
        L = PermutationalLindbladian(ps, **ALL_TERMS)
        ours = L.rhs(0.0, X)
        err, err_in = 0.0, 0.0
        for b in range(X.shape[0]):
            full = _full_rhs(F, ALL_TERMS, _embed(ps, X[b], F))
            err = max(err, (_project(ps, full, F) - ours[b]).abs().max().item())
            err_in = max(err_in, (_embed(ps, _project(ps, full, F), F) - full).abs().max().item())
        res['full'].append(dict(N=N, err=err, err_outside=err_in, dim=4 ** N))
        log(f"  N={N}: |banded - full 2^N generator| = {err:.1e}, full generator leaves M by {err_in:.1e}")
    # trajectory in the full space (RK4 there) against the banded RK4
    N = cfg['traj_N']
    ps, F = PermutationalSystem(N), _full_space(N)
    L = PermutationalLindbladian(ps, **ALL_TERMS)
    rho0 = ps.coherent_state(2.0, 0.5)
    times = torch.linspace(0, 2, 11, dtype=torch.float64)
    ours = L.evolve(rho0, times, substeps=40)
    r = _embed(ps, rho0.matrix[0], F)
    h = 2.0 / 10 / 40
    err = 0.0
    for i in range(1, 11):
        for _ in range(40):
            k1 = _full_rhs(F, ALL_TERMS, r)
            k2 = _full_rhs(F, ALL_TERMS, r + h / 2 * k1)
            k3 = _full_rhs(F, ALL_TERMS, r + h / 2 * k2)
            k4 = _full_rhs(F, ALL_TERMS, r + h * k3)
            r = r + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        err = max(err, (_project(ps, r, F) - ours[i].matrix[0]).abs().max().item())
    res['traj'] = dict(N=N, err=err)
    log(f"  N={N}: full-space trajectory vs banded, max block deviation {err:.1e}")
    # the library's generic route: lindblad_evolve with InterSectorChannel jumps (method 'kraus' and 'rk4')
    N = cfg['lib_N']
    ps = PermutationalSystem(N)
    L = PermutationalLindbladian(ps, **ALL_TERMS)
    H, jumps, rates = L.dynamics_arguments()
    rho0 = ps.coherent_state(2.0, 0.5)
    times = torch.linspace(0, 1, 6, dtype=torch.float64)
    t0 = time.time()
    ours = L.evolve(rho0, times, substeps=40)
    t_ours = time.time() - t0
    t0 = time.time()
    lib = dynamics.lindblad_evolve(rho0, H, jumps, times, rates, substeps=40, method='rk4')
    t_lib = time.time() - t0
    lib_k = dynamics.lindblad_evolve(rho0, H, jumps, times, rates, substeps=40, method='kraus')
    e_rk4 = max((a.matrix - b.matrix).abs().max().item() for a, b in zip(ours, lib))
    e_kraus = max((a.matrix - b.matrix).abs().max().item() for a, b in zip(ours, lib_k))
    tr_kraus = max(abs(float(b.matrix.diagonal(dim1=-2, dim2=-1).sum().real) - 1) for b in lib_k)
    res['lib'] = dict(N=N, err_rk4=e_rk4, err_kraus=e_kraus, trace_kraus=tr_kraus, sec_banded=t_ours,
                      sec_generic=t_lib)
    log(f"  N={N}: generic lindblad_evolve (InterSectorChannel jumps) vs banded: rk4 {e_rk4:.1e}, "
        f"kraus {e_kraus:.1e} (second order); {t_lib:.2f} s vs {t_ours:.2f} s")
    return res


def _to_piqs(ps, X):
    from scipy.linalg import block_diag
    return block_diag(*[X[c, :ps.two_j[c] + 1, :ps.two_j[c] + 1].cpu().numpy() for c in reversed(range(ps.C))])


def _from_piqs(ps, M):
    X = torch.zeros(ps.C, ps.k_max, ps.k_max, dtype=torch.complex128)
    o = 0
    for c in reversed(range(ps.C)):
        k = ps.two_j[c] + 1
        X[c, :k, :k] = torch.as_tensor(np.asarray(M)[o:o + k, o:o + k])
        o += k
    return X


def _steady_state(L):
    """Null vector of the (real-linear) generator on the real symmetric matrices of the blocks."""
    ps = L.system
    basis = []
    for c, tj in enumerate(ps.two_j):
        for a in range(tj + 1):
            for b in range(a, tj + 1):
                basis.append((c, a, b))
    n = len(basis)
    E = torch.zeros(n, ps.C, ps.k_max, ps.k_max, dtype=torch.float64)
    for i, (c, a, b) in enumerate(basis):
        E[i, c, a, b] = 1
        E[i, c, b, a] = 1
    R = L.rhs(0.0, E)
    cidx = torch.tensor([x[0] for x in basis])
    aidx = torch.tensor([x[1] for x in basis])
    bidx = torch.tensor([x[2] for x in basis])
    G = R[:, cidx, aidx, bidx].T                                  # column i: generator applied to basis element i
    _, S, Vh = torch.linalg.svd(G)
    v = Vh[-1]
    X = torch.zeros(ps.C, ps.k_max, ps.k_max, dtype=torch.float64)
    X[cidx, aidx, bidx] = v
    X[cidx, bidx, aidx] = v
    X = X / torch.diagonal(X, dim1=-2, dim2=-1).sum()
    return X, float(S[-1] / S[0]), float(S[-2] / S[0])


def _piqs_steady(ps, Lq):
    """Null vector of a PIQS Liouvillian restricted to the entries inside the blocks (the other entries of its
    nds x nds Dicke-basis matrix are not dynamical, so qutip.steadystate sees a singular matrix)."""
    A = Lq.data.as_scipy().tocsr() if hasattr(Lq.data, 'as_scipy') else Lq.data.tocsr()
    nds = sum(t + 1 for t in ps.two_j)
    idx, o = [], 0
    for c in reversed(range(ps.C)):
        k = ps.two_j[c] + 1
        idx += [(o + cc) * nds + (o + r) for cc in range(k) for r in range(k)]      # column stacking
        o += k
    idx = np.array(idx)
    G = A[idx][:, idx].toarray()
    _, _, Vh = np.linalg.svd(G)
    v = Vh[-1].conj()
    M = np.zeros((nds, nds), dtype=complex)
    M[idx % nds, idx // nds] = v
    return M / np.trace(M)


def part_piqs(cfg, log):
    try:
        import qutip
        from qutip import piqs
    except ImportError:
        log("  qutip.piqs not installed: skipped")
        return None
    res = dict(generator=[], traj=[], steady=[], version=qutip.__version__)
    # generator with every term, against Dicke(...).liouvillian(); PIQS writes gamma_D sum_i D[sigma_z^i / 2]
    for N in cfg['piqs_N'][:1]:
        ps = PermutationalSystem(N)
        L = PermutationalLindbladian(ps, **ALL_TERMS)
        p = ALL_TERMS
        d = piqs.Dicke(N, emission=p['emission'], pumping=p['pumping'], dephasing=4 * p['dephasing'],
                       collective_emission=p['collective_emission'], collective_pumping=p['collective_pumping'],
                       collective_dephasing=p['collective_dephasing'])
        jx, jy, jz = piqs.jspin(N)
        d.hamiltonian = p['omega_x'] * jx + p['omega_y'] * jy + p['omega_z'] * jz
        Lq = d.liouvillian()
        torch.manual_seed(2)
        X = torch.randn(ps.C, ps.k_max, ps.k_max, dtype=torch.complex128)
        X = (X + X.conj().transpose(-2, -1)) * ps.algebra._const('block_mask', torch.bool)
        rq = qutip.vector_to_operator(Lq * qutip.operator_to_vector(qutip.Qobj(_to_piqs(ps, X)))).full()
        a, b = _from_piqs(ps, rq), L.rhs(0.0, X[None])[0]
        res['generator'].append(dict(N=N, err=(a - b).abs().max().item() / b.abs().max().item()))
        log(f"  PIQS generator N={N}: relative deviation {res['generator'][-1]['err']:.1e} "
            "(PIQS stores its rates in single precision)")
    # BTC trajectories with local emission and dephasing
    for N in cfg['piqs_N']:
        ge, gd = 0.1, 0.025
        ps = PermutationalSystem(N, complex_valued=False)
        L = PermutationalLindbladian(ps, omega_y=OMEGA, collective_emission=KAPPA / (N / 2), emission=ge,
                                     dephasing=gd)
        times = np.arange(0, 10.0 + 1e-9, 0.5)
        t0 = time.time()
        sub = step_size(L, 0.5, safety=0.3)
        ours = L.evolve(ps.coherent_state(math.pi), torch.tensor(times), substeps=sub)
        t_ours = time.time() - t0
        t0 = time.time()
        d = piqs.Dicke(N, emission=ge, dephasing=4 * gd, collective_emission=KAPPA / (N / 2))
        d.hamiltonian = OMEGA * piqs.jspin(N, 'y')
        Lq = d.liouvillian()
        t_build = time.time() - t0
        rho0 = piqs.dicke(N, N / 2, -N / 2)
        jz = piqs.jspin(N, 'z')
        t0 = time.time()
        out = qutip.mesolve(Lq, rho0, times, e_ops=[jz], options=dict(atol=1e-12, rtol=1e-10, nsteps=10 ** 6,
                                                                       store_final_state=True))
        t_q = time.time() - t0
        mz_q = np.real(np.asarray(out.expect[0])) / (N / 2)
        mz_o = np.array([float((torch.diagonal(r.matrix[0], dim1=-2, dim2=-1) * ps.m_values()).sum())
                         for r in ours]) / (N / 2)
        fin = _from_piqs(ps, out.final_state.full())
        err_blocks = (fin - ours[-1].matrix[0].to(torch.complex128)).abs().max().item()
        res['traj'].append(dict(N=N, err_mz=float(np.abs(mz_q - mz_o).max()), err_blocks=err_blocks,
                                sec_ours=t_ours, sec_piqs_build=t_build, sec_piqs_solve=t_q))
        log(f"  PIQS trajectory N={N}: max |m_z - m_z^PIQS| = {res['traj'][-1]['err_mz']:.1e}, final blocks "
            f"{err_blocks:.1e}; ours {t_ours:.2f} s, PIQS {t_build:.2f} s (Liouvillian) + {t_q:.2f} s (mesolve)")
    for N in cfg['piqs_ss_N']:
        ge, gd = 0.1, 0.025
        ps = PermutationalSystem(N, complex_valued=False)
        L = PermutationalLindbladian(ps, omega_y=OMEGA, collective_emission=KAPPA / (N / 2), emission=ge,
                                     dephasing=gd)
        X, s_min, s_gap = _steady_state(L)
        d = piqs.Dicke(N, emission=ge, dephasing=4 * gd, collective_emission=KAPPA / (N / 2))
        d.hamiltonian = OMEGA * piqs.jspin(N, 'y')
        Y = _from_piqs(ps, _piqs_steady(ps, d.liouvillian()))
        err = (Y - X.to(torch.complex128)).abs().max().item()
        mz = float((torch.diagonal(X, dim1=-2, dim2=-1) * ps.m_values()).sum()) / (N / 2)
        res['steady'].append(dict(N=N, err=err, mz=mz, s_min=s_min, s_gap=s_gap))
        log(f"  PIQS steady state N={N}: max block deviation {err:.1e}, m_z = {mz:+.5f} "
            f"(sigma_min/sigma_max = {s_min:.1e}, next {s_gap:.1e})")
    return res


# ----------------------------------------------------------------------------------------------------------------
# analytic limits
# ----------------------------------------------------------------------------------------------------------------
def part_superradiance(cfg, dev, log):
    """Pure collective decay D[J_-] (Gamma = 1) from |N/2, N/2>, in the single sector j = N/2."""
    rows = []
    for N in cfg['sr_N']:
        ps = PermutationalSystem(N, complex_valued=False, device=dev, j_min=N / 2)
        L = PermutationalLindbladian(ps, collective_emission=1.0)
        T = 40.0 / N
        n_out = 4000
        dt = T / n_out
        sub = step_size(L, dt, safety=0.5)
        rho = ps.dicke_state(N / 2, N / 2).matrix.to(dev)
        m = ps.m_values()[0].to(dev, torch.float64)
        p = ps.raising_coefficients()[0].to(dev, torch.float64)
        rate = torch.zeros_like(p)
        rate[:-1] = p[1:] ** 2                                    # <a|J_+ J_-|a>
        I, Ppos = [], []
        t0 = time.time()
        I.append(float((torch.diagonal(rho[0, 0]) * rate).sum()))
        Ppos.append(1.0)
        for i in range(n_out):
            rho = dynamics.rk4(L.rhs, rho, [i * dt, (i + 1) * dt], substeps=sub)[-1]
            dg = torch.diagonal(rho[0, 0])
            I.append(float((dg * rate).sum()))
            Ppos.append(float(dg[m > 0.25].sum()))
        secs = time.time() - t0
        ts = np.arange(n_out + 1) * dt
        I, Ppos = np.array(I), np.array(Ppos)
        integral = float(np.sum((Ppos[1:] + Ppos[:-1]) / 2) * dt)    # trapezoid
        exact = sum(1.0 / k for k in range(1, N + 1)) / (N + 1)
        ip = int(np.argmax(I))
        if 0 < ip < len(I) - 1:                                  # parabolic refinement of the peak
            y0, y1, y2 = I[ip - 1], I[ip], I[ip + 1]
            off = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
            tp, Imax = ts[ip] + off * dt, y1 - 0.25 * (y0 - y2) * off
        else:
            tp, Imax = ts[ip], I[ip]
        rows.append(dict(N=N, integral=integral, exact=exact, rel=abs(integral / exact - 1), t_peak=float(tp),
                         t_scaled=float(tp * N / math.log(N)), I_scaled=float(Imax / (N * N / 4)), seconds=secs))
        log(f"  N={N:4d}: int P(m>0) = {integral:.8f}, H_N/(N+1) = {exact:.8f} (rel {rows[-1]['rel']:.1e}); "
            f"t_peak N/ln N = {rows[-1]['t_scaled']:.4f}, I_max/(N^2/4) = {rows[-1]['I_scaled']:.4f} [{secs:.1f} s]")
    return rows


def part_product(cfg, dev, log):
    """Local emission, pumping and dephasing alone: the product state stays a product (single-qubit Bloch
    equations); compare <J_z>, <J_x>, <J^2> and the spin distribution p_j with the closed forms."""
    N = cfg['prod_N']
    ge, gp, gd = 0.3, 0.1, 0.2
    theta = 2 * math.pi / 3
    ps = PermutationalSystem(N, complex_valued=False, device=dev)
    L = PermutationalLindbladian(ps, emission=ge, pumping=gp, dephasing=gd)
    times = [0.0, 0.5, 1.0, 2.0]
    sub = step_size(L, 0.5, safety=0.3)
    t0 = time.time()
    rho = ps.coherent_state(theta).matrix.to(dev)
    errs = dict(jz=0.0, jx=0.0, j2=0.0, pj=0.0, trace=0.0)
    jx = ps.jx().matrix.to(dev).real
    lw = torch.tensor([math.log(float(mm)) for mm in ps.multiplicities], dtype=torch.float64)
    for i in range(1, len(times)):
        rho = dynamics.rk4(L.rhs, rho, [times[i - 1], times[i]], substeps=sub * int(round((times[i] - times[i - 1])
                                                                                        / 0.5)))[-1]
        t = times[i]
        g = ge + gp
        sz = (gp - ge) / g + (math.cos(theta) - (gp - ge) / g) * math.exp(-g * t)
        sx = math.sin(theta) * math.exp(-(g / 2 + 2 * gd) * t)
        dg = torch.diagonal(rho[0], dim1=-2, dim2=-1).cpu()
        p_j = dg.sum(-1)
        jz = float((dg * ps.m_values()).sum())
        jxv = float((rho[0] * jx.transpose(-2, -1)).sum())
        jv = torch.tensor([tj / 2 for tj in ps.two_j], dtype=torch.float64)
        j2 = float((p_j * jv * (jv + 1)).sum())
        s2 = sx * sx + sz * sz
        j2_exact = 3 * N / 4 + N * (N - 1) * s2 / 4
        lam_p, lam_m = (1 + math.sqrt(s2)) / 2, (1 - math.sqrt(s2)) / 2
        logp = []
        for c, tj in enumerate(ps.two_j):
            j = tj / 2
            ms = np.arange(-j, j + 1)
            terms = (N / 2 + ms) * math.log(lam_p) + (N / 2 - ms) * math.log(lam_m)
            mx = terms.max()
            logp.append(float(lw[c]) + mx + math.log(np.exp(terms - mx).sum()))
        p_exact = torch.tensor(np.exp(logp))
        errs['jz'] = max(errs['jz'], abs(jz - N / 2 * sz) / (N / 2))
        errs['jx'] = max(errs['jx'], abs(jxv - N / 2 * sx) / (N / 2))
        errs['j2'] = max(errs['j2'], abs(j2 - j2_exact) / j2_exact)
        errs['pj'] = max(errs['pj'], float((p_j - p_exact).abs().sum()))
        errs['trace'] = max(errs['trace'], abs(float(p_j.sum()) - 1))
    secs = time.time() - t0
    log(f"  N={N}: rel. errors <J_z> {errs['jz']:.1e}, <J_x> {errs['jx']:.1e}, <J^2> {errs['j2']:.1e}; "
        f"sum_j |p_j - p_j^exact| = {errs['pj']:.1e}; trace {errs['trace']:.1e} [{secs:.1f} s]")
    return dict(N=N, errs=errs, seconds=secs, rates=[ge, gp, gd], theta=theta)


# ----------------------------------------------------------------------------------------------------------------
def run(mode, dev='cpu', log=print):
    cfg = CONFIGS[mode]
    warnings.simplefilter('ignore', cost.CostWarning)
    torch.set_flush_denormal(True)                       # denormal numbers in the far tails slow the CPU down
    t0 = time.time()
    res = dict(config=cfg, omega=OMEGA, kappa=KAPPA)
    log("generator: full Hilbert space, generic library route")
    res['generator'] = part_generator(cfg, log)
    log("PIQS")
    res['piqs'] = part_piqs(cfg, log)
    log("superradiance (pure collective decay)")
    res['superradiance'] = part_superradiance(cfg, dev, log)
    log("product states under local dissipation")
    res['product'] = part_product(cfg, dev, log)
    log("boundary time crystal: tau(N, gamma)")
    res['sweep'] = part_sweep(cfg, dev, log)
    res['seconds'] = time.time() - t0
    log(f"total {res['seconds']:.0f} s")
    return res


# ----------------------------------------------------------------------------------------------------------------
# output
# ----------------------------------------------------------------------------------------------------------------
def _num(x, d=2):
    return '$>$' if not math.isfinite(x) else f"{x:.{d}f}"


def write_outputs(args, res, rerender=False):
    """Write generated/permutational.{json,tex} and the figure; rerender=True (PM_FROM_JSON) keeps the saved JSON and
    takes the device macro from the environment stored in it."""
    from common import save_json, write_tex, env_macro, sci
    cfg = res['config']
    sw = res['sweep']
    if not rerender:
        save_json(args, 'permutational', res)
    tex = env_macro(args, 'PM')
    env = res.get('environment')
    if rerender and env:
        hw = env.get('gpu') or f"CPU ({env.get('threads')} thread{'' if env.get('threads') == 1 else 's'})"
        tex = tex.split('\n', 1)[1]
        tex = f"\\newcommand{{\\PMDevice}}{{{hw}}}\n" + tex.replace(f"{{\\PMMode}}{{{args.mode}}}",
                                                                    f"{{\\PMMode}}{{{res.get('mode', args.mode)}}}")
    gen = res['generator']
    tex += (f"\\newcommand{{\\PMOmega}}{{{OMEGA:g}}}\n"
            f"\\newcommand{{\\PMFullNmax}}{{{max(r['N'] for r in gen['full'])}}}\n"
            f"\\newcommand{{\\PMFullErr}}{{{sci(max(max(r['err'] for r in gen['full']), 1e-17), 1)}}}\n"
            f"\\newcommand{{\\PMTrajN}}{{{gen['traj']['N']}}}\n"
            f"\\newcommand{{\\PMTrajErr}}{{{sci(max(gen['traj']['err'], 1e-17), 1)}}}\n"
            f"\\newcommand{{\\PMLibN}}{{{gen['lib']['N']}}}\n"
            f"\\newcommand{{\\PMLibErr}}{{{sci(max(gen['lib']['err_rk4'], 1e-17), 1)}}}\n"
            f"\\newcommand{{\\PMLibKrausErr}}{{{sci(max(gen['lib']['err_kraus'], 1e-17), 1)}}}\n"
            f"\\newcommand{{\\PMLibSecGeneric}}{{{gen['lib']['sec_generic']:.2f}}}\n"
            f"\\newcommand{{\\PMLibSecBanded}}{{{gen['lib']['sec_banded']:.2f}}}\n")
    pq = res['piqs']
    if pq:
        tr = pq['traj']
        tex += (f"\\newcommand{{\\PMPiqsVersion}}{{{pq['version']}}}\n"
                f"\\newcommand{{\\PMPiqsGenErr}}{{{sci(max(r['err'] for r in pq['generator']), 1)}}}\n"
                f"\\newcommand{{\\PMPiqsNs}}{{{', '.join(str(r['N']) for r in tr)}}}\n"
                f"\\newcommand{{\\PMPiqsNmax}}{{{tr[-1]['N']}}}\n"
                f"\\newcommand{{\\PMPiqsTrajErr}}{{{sci(max(r['err_mz'] for r in tr), 1)}}}\n"
                f"\\newcommand{{\\PMPiqsBlockErr}}{{{sci(max(r['err_blocks'] for r in tr), 1)}}}\n"
                f"\\newcommand{{\\PMPiqsSecOurs}}{{{tr[-1]['sec_ours']:.1f}}}\n"
                f"\\newcommand{{\\PMPiqsSecBuild}}{{{tr[-1]['sec_piqs_build']:.1f}}}\n"
                f"\\newcommand{{\\PMPiqsSecSolve}}{{{tr[-1]['sec_piqs_solve']:.1f}}}\n"
                f"\\newcommand{{\\PMPiqsSsNs}}{{{', '.join(str(r['N']) for r in pq['steady'])}}}\n"
                f"\\newcommand{{\\PMPiqsSsErr}}{{{sci(max(r['err'] for r in pq['steady']), 1)}}}\n")
    sr = res['superradiance']
    tex += "\\newcommand{\\PMSRRows}{%\n"
    for r in sr:
        tex += (f"{r['N']} & {r['integral']:.6f} & {r['exact']:.6f} & ${sci(max(r['rel'], 1e-16), 1)}$ & "
                f"{r['t_scaled']:.3f} & {r['I_scaled']:.3f} \\\\\n")
    tex += "}\n"
    tex += (f"\\newcommand{{\\PMSRErr}}{{{sci(max(max(r['rel'] for r in sr), 1e-16), 1)}}}\n"
            f"\\newcommand{{\\PMSRNmax}}{{{sr[-1]['N']}}}\n")
    pr = res['product']
    tex += (f"\\newcommand{{\\PMProdN}}{{{pr['N']}}}\n"
            f"\\newcommand{{\\PMProdErr}}{{{sci(max(max(pr['errs'][k] for k in ('jz', 'jx', 'j2')), 1e-16), 1)}}}\n"
            f"\\newcommand{{\\PMProdPjErr}}{{{sci(max(pr['errs']['pj'], 1e-16), 1)}}}\n"
            f"\\newcommand{{\\PMProdSec}}{{{pr['seconds']:.0f}}}\n")
    # the sweep
    ge, gd = cfg['ge'], cfg['gd']
    tex += (f"\\newcommand{{\\PMNs}}{{{', '.join(str(n) for n in sw['N'])}}}\n"
            f"\\newcommand{{\\PMNmax}}{{{sw['N'][-1]}}}\n"
            f"\\newcommand{{\\PMT}}{{{cfg['T']:g}}}\n"
            f"\\newcommand{{\\PMBatch}}{{{len(ge) + len(gd)}}}\n"
            f"\\newcommand{{\\PMGe}}{{{', '.join(f'{g:g}' for g in ge)}}}\n"
            f"\\newcommand{{\\PMGd}}{{{', '.join(f'{g:g}' for g in gd)}}}\n"
            f"\\newcommand{{\\PMConvN}}{{{sw['conv']['N']}}}\n"
            f"\\newcommand{{\\PMConvMz}}{{{sci(max(sw['conv']['dmz'], 1e-17), 1)}}}\n"
            f"\\newcommand{{\\PMConvTau}}{{{sci(max(sw['conv']['dtau'], 1e-17), 1)}}}\n"
            f"\\newcommand{{\\PMFitNmin}}{{{FIT_NMIN}}}\n"
            f"\\newcommand{{\\PMRuleNmin}}{{{RULE_NMIN}}}\n"
            f"\\newcommand{{\\PMSlope}}{{{sw['slope']:.3f}}}\n"
            f"\\newcommand{{\\PMIntercept}}{{{sw['intercept']:.2f}}}\n")
    tex += ("\\newcommand{\\PMTauHeader}{$N$ & $0$ & " + ' & '.join(f'${g:g}$' for g in ge) + ' & ' +
            ' & '.join(f'${g:g}$' for g in gd) + "}\n")
    tex += "\\newcommand{\\PMTauRows}{%\n"
    for iN, N in enumerate(sw['N']):
        tex += (f"{N} & {_num(sw['tau0'][iN])} & " + ' & '.join(_num(x) for x in sw['tau_e'][iN]) + ' & ' +
                ' & '.join(_num(x) for x in sw['tau_d'][iN]) + " \\\\\n")
    tex += ("$\\infty$ (MF) & -- & " + ' & '.join(_num(x) for x in sw['tau_mf_e']) + ' & ' +
            ' & '.join(_num(x) for x in sw['tau_mf_d']) + " \\\\\n}\n")
    rel = [d['rel'] for d in sw['rule']]
    big = [d['rel'] for d in sw['rule'] if d['N'] >= RULE_NMIN] or [math.nan]
    tex += (f"\\newcommand{{\\PMRuleMin}}{{{min(rel):+.3f}}}\n"
            f"\\newcommand{{\\PMRuleMax}}{{{max(rel):+.3f}}}\n"
            f"\\newcommand{{\\PMRuleBigMin}}{{{min(big):+.3f}}}\n"
            f"\\newcommand{{\\PMRuleBigMax}}{{{max(big):+.3f}}}\n"
            f"\\newcommand{{\\PMRuleCount}}{{{len(rel)}}}\n"
            f"\\newcommand{{\\PMSecNmax}}{{{sw['seconds'][-1]:.0f}}}\n"
            f"\\newcommand{{\\PMSecSweep}}{{{sum(sw['seconds']) + sum(sw['seconds0']):.0f}}}\n"
            f"\\newcommand{{\\PMSeconds}}{{{res['seconds']:.0f}}}\n")
    Nmax = sw['N'][-1]
    tex += (f"\\newcommand{{\\PMDimNmax}}{{{sum((t + 1) ** 2 for t in range(Nmax % 2, Nmax + 1, 2))}}}\n"
            f"\\newcommand{{\\PMDimFullNmax}}{{{sci(4.0 ** Nmax, 1)}}}\n")
    write_tex(args, 'permutational', tex)
    make_figure(args, res)


def make_figure(args, res):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter, NullFormatter
    cfg, sw = res['config'], res['sweep']
    ge, gd = cfg['ge'], cfg['gd']
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 4.6))
    ax = axes[0, 0]
    Nmax = sw['N'][-1]
    tr = sw['traces'][Nmax] if Nmax in sw['traces'] else sw['traces'][str(Nmax)]
    ax.plot(tr['t0'], tr['mz0'], color='k', lw=0.8, label='$\\gamma=0$')
    t = np.array(tr['t'])
    for ig, g in enumerate(ge[:4]):
        if ig % 2:
            continue
        z = np.array([np.nan if v is None else v for v in tr['mz'][ig]])
        ok = ~np.isnan(z)
        ax.plot(t[ok], z[ok], color=f'C{ig}', lw=0.9, label=f'$\\gamma={g:g}\\kappa$')
        zm = np.array(sw['mf_traces'][str(g)])
        n = min(len(zm), len(t))
        ax.plot(t[:n], zm[:n], color=f'C{ig}', lw=0.6, ls=':')
    ax.set_xlim(0, min(cfg['T'], 40))
    ax.set_ylim(-1.05, 1.45)
    ax.set_xlabel('$\\kappa t$')
    ax.set_ylabel('$m_z=\\langle J_z\\rangle/S$')
    ax.set_title(f'(a) $N={Nmax}$, local emission (dotted: mean field)', fontsize=9)
    ax.legend(frameon=False, fontsize=6.5, loc='upper center', ncol=3)
    ax = axes[0, 1]
    Ns = np.array(sw['N'])
    ax.plot(Ns, sw['tau0'], 'ko-', ms=3, lw=0.8, label='$\\gamma=0$')
    for ig, g in enumerate(ge):
        tau = np.array([sw['tau_e'][i][ig] for i in range(len(Ns))])
        ok = np.isfinite(tau)
        ax.plot(Ns[ok], tau[ok], 'o-', color=f'C{ig}', ms=3, lw=0.8, label=f'$\\gamma={g:g}\\kappa$')
        tm = sw['tau_mf_e'][ig]
        if math.isfinite(tm):
            ax.axhline(tm, color=f'C{ig}', lw=0.5, ls=':')
            rule = 1 / (1 / np.array(sw['tau0']) + 1 / tm)
            ax.plot(Ns, rule, color=f'C{ig}', lw=0.5, ls='--')
    ax.set_xscale('log')
    ax.set_yscale('log')
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_formatter(ScalarFormatter())
        axis.set_minor_formatter(NullFormatter())
    ax.set_xticks([n for n in (10, 20, 50, 100, 200) if Ns.min() <= n <= Ns.max()] or list(Ns))
    ax.set_yticks([v for v in (2, 5, 10, 20, 50, 100) if v <= 1.5 * max(sw['tau0'])])
    ax.set_xlabel('$N$')
    ax.set_ylabel('$\\kappa\\tau$')
    ax.set_title('(b) lifetime, local emission (dashed: rule, dotted: MF)', fontsize=9)
    ax.legend(frameon=False, fontsize=6, loc='upper left', ncol=2)

    def series(kind):
        mf, gs = (sw['tau_mf_e'], ge) if kind == 'tau_e' else (sw['tau_mf_d'], gd)
        for ig, g in enumerate(gs):
            tm = mf[ig]
            if not math.isfinite(tm):
                continue
            tau = np.array([sw[kind][i][ig] for i in range(len(Ns))])
            ok = np.isfinite(tau)
            t0 = np.array(sw['tau0'])[ok]
            yield ig, g, t0 / tm, tau[ok] / tm, tau[ok] * (1 / t0 + 1 / tm) - 1
    ax = axes[1, 0]
    xs = np.logspace(-2, 1.5, 200)
    ax.plot(xs, xs / (1 + xs), color='0.5', lw=0.8, ls='--', label='$x/(1+x)$')
    for kind, mk in (('tau_e', 'o'), ('tau_d', 's')):
        for ig, g, x, y, _ in series(kind):
            ax.plot(x, y, mk, color=f'C{ig}', ms=3, mfc='none' if kind == 'tau_d' else f'C{ig}')
    ax.plot([], [], 'o', color='0.3', ms=3, label='local emission')
    ax.plot([], [], 's', color='0.3', ms=3, mfc='none', label='local dephasing')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$x=\\tau(N,0)/\\tau_{\\mathrm{MF}}$')
    ax.set_ylabel('$\\tau/\\tau_{\\mathrm{MF}}$')
    ax.set_title('(c) all rates and $N$', fontsize=9)
    ax.legend(frameon=False, fontsize=6.5, loc='lower right')
    ax = axes[1, 1]
    ax.axhline(0, color='0.6', lw=0.6)
    for kind, mk, ls in (('tau_e', 'o', '-'), ('tau_d', 's', '--')):
        for ig, g, x, y, dev in series(kind):
            ax.plot(x, 100 * dev, mk, ls=ls, lw=0.6, color=f'C{ig}', ms=3,
                    mfc='none' if kind == 'tau_d' else f'C{ig}')
    ax.set_xscale('log')
    ax.set_xlabel('$x=\\tau(N,0)/\\tau_{\\mathrm{MF}}$')
    ax.set_ylabel('$\\tau/\\tau_{\\mathrm{rule}}-1$ [%]')
    ax.set_title('(d) deviation from the rule (lines: fixed rate, $N$ grows)', fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{args.figdir}/permutational.pdf", metadata={"CreationDate": None})


def main():
    from common import parse_args
    args = parse_args(__doc__)
    src = os.environ.get('PM_FROM_JSON')                 # re-render tex and figure from a saved run
    if src:
        import json
        with open(src) as f:
            res = json.load(f)
    else:
        res = run(args.mode, args.device, log=lambda s: print(s, flush=True))
    write_outputs(args, res, rerender=bool(src))
    print("written:", args.out, args.figdir)


if __name__ == '__main__':
    main()
