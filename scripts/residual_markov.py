#!/usr/bin/env python3
"""
Residual channel networks express Markovian evolutions, but no channel with negative determinant (man/article,
Section "Residual channel networks and Markovianity", sec:learn-markov).

Algebra M_2 (one qubit). Three target channels, given by their Pauli transfer matrices R_ab = Tr(s_a Phi(s_b))/2
(s = 1, X, Y, Z; R is the superoperator in the orthonormal basis s_a/sqrt2, so ||R - R'||_F is the Hilbert-Schmidt
distance of the superoperators and det R their determinant):

  (i)   'gksl':    exp(L) for a random GKSL generator L (random Hamiltonian, two random jump operators);
  (ii)  'damped':  amplitude damping (gamma = 0.5) after the rotation exp(-i theta n.s/2), theta = 2 pi/3,
                   n = (1, 1, 0)/sqrt2: a product of two Markovian channels (divisible);
  (iii) 'pauli':   rho -> (X rho X + Y rho Y + Z rho Z)/3, R = diag(1, -1/3, -1/3, -1/3), det R = -1/27 < 0: not
                   infinitesimally divisible (Wolf-Cirac, CMP 279 (2008) 147, Prop. 15).

Models (torch_vn_algebra.learn):

  (a) 'channel':   one ChannelLayer of Kraus rank 4 (all qubit channels);
  (b) residual:    nn.Sequential of L blocks Residual(ChannelLayer(alg, 4), p) with p = P_MAX fixed, for
                   L in DEPTHS. A learnable p <= P_MAX adds nothing: the inner layer ranges over all channels,
                   including (1 - s) id + s Psi, and (1 - p) id + p((1 - s) id + s Psi) is the block of weight p s;
  (c) 'lindblad':  one LindbladLayer (time 1, Hamiltonian and three free jump operators: every qubit GKSL
                   generator), a block of continuous depth.

Data: the outputs of the target on the four states |0>, |1>, |+>, |+i> (a tomographically complete set); the loss
is the squared Hilbert-Schmidt distance ||R_net - R_target||_F^2, computed from the outputs by linearity.
Restarts are run in parallel as sectors: the network is built on the algebra (+)_{c<C} M_2 with C = 3 targets x
RESTARTS sectors; every layer of learn is sector preserving with independent parameters in each sector (Kraus
blocks, Hamiltonian and jump blocks), so sector c is an independent network fitted to target c // RESTARTS, and Adam
(elementwise) treats them as independent runs. Every model is trained with Adam (cosine-decayed step) for STEPS
steps; the best restart is reported.

Reported per target and model: the Hilbert-Schmidt error ||R_net - R_t||_F, the worst-case trace distance over pure
inputs sup_psi ||Phi_net(psi) - Phi_t(psi)||_1 / 2 (a lower bound on half the diamond distance; computed on the
Bloch sphere by a grid and projected-gradient refinement), and det R_net.

Whether (i) and (ii) are themselves of the form e^L is decided by the test of Wolf, Eisert, Cubitt and Cirac
(PRL 101 (2008) 150402), run on the branches |k| <= 3 of the logarithm (semigroup_test): (i) passes on the principal
branch, (ii) fails on all of them, so a gap between the Lindblad layer and the Stinespring layer is expected for (ii).

Proved (and checked on every learned map): every block (1 - p) id + p Phi with p < 1/2 has det > 0 (eigenvalues
1 - p + p lambda, |lambda| <= 1, real ones >= 1 - 2p, complex ones in conjugate pairs), and so has every product of
blocks and every exp(L) (det = exp tr L). For the Pauli target, every trace-preserving map with det R >= 0 has
||R - R_t||_F >= 1/3 and worst-case trace distance >= 1/6 (Eckart-Young: the nearest singular 3x3 matrix to -I/3
is at distance 1/3 in the Frobenius and in the operator norm). The Pauli channel diag(1, -1/3, -1/3, 0) attains
both bounds. The trace-distance bound is sharp for infinitesimally divisible channels: Ad_X o (Pauli semigroup
diag(1, e, e, 1/3)) tends to diag(1, 0, 0, -1/3) as e -> 0, at worst-case trace distance 1/6 (Frobenius sqrt2/3).
The semigroup exp(L) = Ad_Z o (Pauli semigroup) with R = diag(1, -a, -a, a^2), a the root of a^3 + 4a/3 - 1/3 = 0,
is the closest among Pauli semigroups composed with a Pauli rotation by pi (|lambda_z| >= |lambda_x lambda_y| for
Pauli semigroups; checked by brute force); it is an upper reference for the floor of the Lindblad layer.

Run time about 6.5 minutes on 4 CPU cores (one process per job, one thread each). Writes
man/article/generated/residual_markov.{json,tex} (macros \\RM...) and man/article/figures/residual_markov.pdf.

    python scripts/residual_markov.py [--steps 2000] [--restarts 6] [--workers 4] [--replot]
"""
import argparse
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import torch
from torch import nn

from common import environment                       # device, versions, matmul precision, commit

from torch_vn_algebra import DensityMatrix, TypeIAlgebra, learn

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'man' / 'article'
DEPTHS = (1, 2, 4, 8, 16, 32)
P_MAX = 0.2                                  # weight of every residual block (< 1/2: det > 0)
TARGETS = ('gksl', 'damped', 'pauli')
LABELS = {'gksl': r'(i) $e^{\mathcal{L}}$, random GKSL', 'damped': r'(ii) damping $\circ$ rotation',
          'pauli': r'(iii) Pauli, $\det=-1/27$'}
INK, GREY = '#0b0b0b', '#7d848c'
BLUE, ORANGE, GREEN, YELLOW = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'

DT = torch.complex128
I2 = torch.eye(2, dtype=DT)
SX = torch.tensor([[0, 1], [1, 0]], dtype=DT)
SY = torch.tensor([[0, -1j], [1j, 0]], dtype=DT)
SZ = torch.tensor([[1, 0], [0, -1]], dtype=DT)
PAULIS = torch.stack([I2, SX, SY, SZ])
_kets = torch.tensor([[1, 0], [0, 1], [1, 1], [1, 1j]], dtype=DT)
_kets = _kets / _kets.norm(dim=1, keepdim=True)
PROBES = torch.einsum('bi,bj->bij', _kets, _kets.conj())          # |0>, |1>, |+>, |+i>


# ----------------------------------------------------------------------
# Pauli transfer matrices
# ----------------------------------------------------------------------
def ptm_from_outputs(Y: torch.Tensor) -> torch.Tensor:
    """Y (4, ..., 2, 2): images of |0>,|1>,|+>,|+i>. Returns R (..., 4, 4), R_ab = Tr(s_a Phi(s_b))/2."""
    phi_i = Y[0] + Y[1]
    phi_z = Y[0] - Y[1]
    phi_x = 2 * Y[2] - phi_i
    phi_y = 2 * Y[3] - phi_i
    imgs = torch.stack([phi_i, phi_x, phi_y, phi_z], dim=-3)       # (..., b, 2, 2)
    return 0.5 * torch.einsum('aij,...bji->...ab', PAULIS, imgs).real


def target_maps(seed: int):
    """The three targets as functions of a batch of 2x2 matrices, plus the data that defines them."""
    g = torch.Generator().manual_seed(seed)

    def cgauss(*shape):
        return torch.complex(torch.randn(*shape, generator=g, dtype=torch.float64),
                             torch.randn(*shape, generator=g, dtype=torch.float64)) / math.sqrt(2)
    A = cgauss(2, 2)
    H = 0.5 * (A + A.conj().T)
    J = [0.5 * cgauss(2, 2) for _ in range(2)]
    # row-major vec: vec(A X B) = (A kron B^T) vec(X)
    Lsup = -1j * (torch.kron(H, I2) - torch.kron(I2, H.T.contiguous()))
    for Jk in J:
        JJ = Jk.conj().T @ Jk
        Lsup = (Lsup + torch.kron(Jk, Jk.conj().resolve_conj()) - 0.5 * torch.kron(JJ, I2)
                - 0.5 * torch.kron(I2, JJ.T.contiguous()))
    E = torch.linalg.matrix_exp(Lsup)

    def gksl(X):
        return (E @ X.reshape(*X.shape[:-2], 4, 1)).reshape(X.shape)

    gamma, theta = 0.5, 2 * math.pi / 3
    n = torch.tensor([1.0, 1.0, 0.0], dtype=torch.float64) / math.sqrt(2)
    U = math.cos(theta / 2) * I2 - 1j * math.sin(theta / 2) * (n[0] * SX + n[1] * SY + n[2] * SZ)
    K0 = torch.tensor([[1, 0], [0, math.sqrt(1 - gamma)]], dtype=DT)
    K1 = torch.tensor([[0, math.sqrt(gamma)], [0, 0]], dtype=DT)

    def damped(X):
        Y = U @ X @ U.conj().T
        return K0 @ Y @ K0.conj().T + K1 @ Y @ K1.conj().T

    def pauli(X):
        return (SX @ X @ SX + SY @ X @ SY + SZ @ X @ SZ) / 3

    info = dict(H_re_im=[H.real.tolist(), H.imag.tolist()], jumps_re_im=[[Jk.real.tolist(), Jk.imag.tolist()] for Jk in J],
                gamma=gamma, theta=theta, axis=n.tolist())
    return dict(gksl=gksl, damped=damped, pauli=pauli), info


def target_ptms(seed: int) -> torch.Tensor:
    maps, _ = target_maps(seed)
    return torch.stack([ptm_from_outputs(maps[t](PROBES)) for t in TARGETS])       # (3, 4, 4)


def semigroup_test(f, kmax: int = 3) -> dict:
    """
    The test of Wolf, Eisert, Cubitt and Cirac (PRL 101 (2008) 150402) for T = e^L on a qubit with non-degenerate
    spectrum: the logarithms of T are V diag(log mu + 2 pi i k_j) V^{-1}, with opposite k on conjugate eigenvalues
    (Hermiticity preserving) and k = 0 on real ones; T is Markovian iff one of them is conditionally completely
    positive, i.e. w J w >= 0 for its Choi matrix J and w the projection onto the complement of the maximally entangled
    vector. Returns, per branch k (|k| <= kmax), the smallest eigenvalue of w J w; >= 0 (up to rounding) means e^L.
    """
    S = torch.zeros(4, 4, dtype=DT)
    for i in range(2):
        for j in range(2):
            E = torch.zeros(1, 2, 2, dtype=DT)
            E[0, i, j] = 1
            S[:, 2 * i + j] = f(E)[0].reshape(4)
    S = S.numpy()
    mu, V = np.linalg.eig(S)
    omega = np.eye(2).reshape(4)
    w = np.eye(4) - np.outer(omega, omega) / 2
    out = {}
    for k in range(-kmax, kmax + 1):
        lg = np.log(mu.astype(complex)) + 2j * np.pi * k * np.sign(np.where(abs(mu.imag) > 1e-9, mu.imag, 0.0))
        L = V @ np.diag(lg) @ np.linalg.inv(V)
        J = L.reshape(2, 2, 2, 2).transpose(2, 0, 3, 1).reshape(4, 4)          # J[(i,a),(j,b)] = L[(a,b),(i,j)]
        herm = float(abs(J - J.conj().T).max())
        out[k] = dict(min_eig=float(np.linalg.eigvalsh(w @ (0.5 * (J + J.conj().T)) @ w).min()), herm_err=herm)
    return dict(branches=out, best=max(v['min_eig'] for v in out.values()),
                degenerate=bool(min(abs(mu[a] - mu[b]) for a in range(4) for b in range(a)) < 1e-6))


def worst_trace_distance(R: np.ndarray, Rt: np.ndarray, grid: int = 4000, iters: int = 300) -> float:
    """sup over pure inputs of ||Phi(psi) - Phi_t(psi)||_1/2 = sup_{|n|=1} |dt + dM n| / 2 (TP maps of a qubit)."""
    D = R - Rt
    a, B = D[1:, 0], D[1:, 1:]
    k = np.arange(grid) + 0.5
    phi, z = np.pi * (1 + 5 ** 0.5) * k, 1 - 2 * k / grid
    n = np.stack([np.sqrt(1 - z ** 2) * np.cos(phi), np.sqrt(1 - z ** 2) * np.sin(phi), z], 1)
    f = np.linalg.norm(a + n @ B.T, axis=1)
    best = n[np.argsort(f)[-16:]]
    for _ in range(iters):                       # projected gradient ascent of |a + B n|^2 on the sphere
        v = a + best @ B.T
        best = best + 0.5 * v @ B
        best /= np.linalg.norm(best, axis=1, keepdims=True)
    return 0.5 * float(max(f.max(), np.linalg.norm(a + best @ B.T, axis=1).max()))


# ----------------------------------------------------------------------
# training (one job = one model, all targets and restarts as sectors)
# ----------------------------------------------------------------------
def build(kind: str, depth: int, alg: TypeIAlgebra, gen: torch.Generator) -> nn.Module:
    if kind == 'channel':
        return learn.ChannelLayer(alg, kraus_rank=4, generator=gen)
    if kind == 'lindblad':
        return learn.LindbladLayer(alg, time=1.0, n_jumps=3, init_scale=0.3, generator=gen)
    return nn.Sequential(*[learn.Residual(learn.ChannelLayer(alg, kraus_rank=4, generator=gen), p=P_MAX)
                           for _ in range(depth)])


def train(job):
    kind, depth, seed, steps, restarts, lr, target_seed = job
    torch.set_num_threads(1)
    t0 = time.time()
    gen = torch.Generator().manual_seed(seed)
    C = len(TARGETS) * restarts
    alg = TypeIAlgebra([2] * C, [2] * C, precision='double')
    X = DensityMatrix(alg, matrix=PROBES[:, None].expand(4, C, 2, 2).contiguous(), validate=False)
    Rt = target_ptms(target_seed).repeat_interleave(restarts, 0)                     # (C, 4, 4)
    net = build(kind, depth, alg, gen)
    opt = torch.optim.Adam(net.parameters(), lr=lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, steps, eta_min=lr * 1e-3)
    history = []
    for step in range(steps):
        R = ptm_from_outputs(net(X).matrix)
        per = ((R - Rt) ** 2).sum((-2, -1))
        opt.zero_grad()
        per.sum().backward()
        opt.step()
        sched.step()
        if step % 50 == 0 or step == steps - 1:
            history.append([step] + per.detach().reshape(len(TARGETS), restarts).amin(1).sqrt().tolist())
    with torch.no_grad():
        R = ptm_from_outputs(net(X).matrix)
    err = ((R - Rt) ** 2).sum((-2, -1)).sqrt().reshape(len(TARGETS), restarts)
    det = torch.linalg.det(R).reshape(len(TARGETS), restarts)
    res = {}
    for i, t in enumerate(TARGETS):
        j = int(err[i].argmin())
        c = i * restarts + j
        td = [worst_trace_distance(R[i * restarts + q].numpy(), Rt[i * restarts + q].numpy()) for q in range(restarts)]
        res[t] = dict(hs=float(err[i, j]), hs_restarts=err[i].tolist(), det=float(det[i, j]),
                      det_min=float(det[i].min()), det_restarts=det[i].tolist(), worst_td=td[j],
                      worst_td_min=min(td), td_restarts=td, ptm=R[c].tolist())
    if kind == 'residual':
        dets = []                                  # determinant of every block of the best Pauli network
        c = TARGETS.index('pauli') * restarts + int(err[TARGETS.index('pauli')].argmin())
        with torch.no_grad():
            for block in net:
                Rb = ptm_from_outputs(block(X).matrix)[c]
                dets.append(float(torch.linalg.det(Rb)))
        res['pauli']['block_det_min'] = min(dets)
    return dict(kind=kind, depth=depth, seconds=time.time() - t0, results=res, history=history)


# ----------------------------------------------------------------------
# bounds for the Pauli target
# ----------------------------------------------------------------------
def pauli_references():
    """Lower bounds over maps with det >= 0, and the closest Markovian Pauli channel diag(1, -a, -a, a^2)."""
    roots = np.roots([1.0, 0.0, 4.0 / 3.0, -1.0 / 3.0])
    a = float(min(r.real for r in roots if abs(r.imag) < 1e-12 and r.real > 0))
    Rm = np.diag([1.0, -a, -a, a * a])
    Rt = np.diag([1.0, -1 / 3, -1 / 3, -1 / 3])
    # brute-force check of the optimum over the Markovian Pauli family (signs of Ad_Z, |l_z| >= |l_x l_y|)
    grid = np.linspace(0.0, 1.0, 401)
    lx, ly = np.meshgrid(grid, grid)
    lz = np.maximum(lx * ly, 0.0)
    brute = np.sqrt((1 / 3 - lx) ** 2 + (1 / 3 - ly) ** 2 + (lz + 1 / 3) ** 2).min()
    # Ad_X o (Pauli semigroup diag(1, e, e, 1/3)), e -> 0: a limit of products of Markovian channels
    Rd = np.diag([1.0, 0.0, 0.0, -1 / 3])
    return dict(hs_lower=1 / 3, td_lower=1 / 6, markov_a=a, markov_hs=float(np.linalg.norm(Rm - Rt)),
                markov_td=worst_trace_distance(Rm, Rt), markov_hs_brute=float(brute),
                dephased_hs=float(np.linalg.norm(Rd - Rt)), dephased_td=worst_trace_distance(Rd, Rt))


# ----------------------------------------------------------------------
# figure and macros
# ----------------------------------------------------------------------
def figure(res, refs, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 4, figsize=(10.0, 2.75))
    depths = np.array(DEPTHS)
    for ax, t in zip(axes[:3], TARGETS):
        hs = [res[('residual', L)]['results'][t]['hs'] for L in DEPTHS]
        ax.plot(depths, hs, '-o', color=BLUE, lw=2, ms=5, mfc='white', mew=1.6, label=f'residual, $p={P_MAX}$, fit error',
                zorder=3)
        ax.axhline(res[('lindblad', 0)]['results'][t]['hs'], color=ORANGE, lw=2, ls='--', label='one Lindblad layer')
        ax.axhline(res[('channel', 0)]['results'][t]['hs'], color=GREEN, lw=2, ls='-.',
                   label='one Stinespring layer')
        if t == 'pauli':                                   # linear scale: the floor is the message
            td = [res[('residual', L)]['results'][t]['worst_td_min'] for L in DEPTHS]
            ax.plot(depths, td, '-s', color=YELLOW, lw=2, ms=4.5, mfc='white', mew=1.6, zorder=3,
                    label='residual, worst-case trace distance')
            for b in (refs['hs_lower'], refs['td_lower']):
                ax.axhline(b, color=GREY, lw=1.4, ls=':', label=r'proved bounds $1/3$, $1/6$ ($\det\geq0$)'
                           if b == refs['hs_lower'] else None)
            for b, txt in ((refs['hs_lower'], '1/3'), (refs['td_lower'], '1/6')):
                ax.text(depths[-1] * 1.12, b, txt, fontsize=7, color=GREY, va='center')
            ax.set_ylim(-0.06, 2.0)
        else:
            ax.set_yscale('log')
        ax.set_xscale('log', base=2)
        ax.set_xticks(depths)
        ax.set_xticklabels([str(L) for L in DEPTHS])
        ax.set_xlabel('depth $L$ (residual blocks)', fontsize=8, color=INK)
        ax.set_title(LABELS[t], fontsize=9, color=INK)
        ax.set_ylabel(r'fit error $\|R_{\rm net}-R\|_F$' if t != 'pauli' else 'distance to the target',
                      fontsize=8, color=INK)
    ax = axes[3]
    dets = [res[('residual', L)]['results']['pauli']['det'] for L in DEPTHS]
    ax.plot(depths, dets, '-o', color=BLUE, lw=2, ms=5, mfc='white', mew=1.6, zorder=3)
    ax.axhline(res[('lindblad', 0)]['results']['pauli']['det'], color=ORANGE, lw=2, ls='--')
    ax.axhline(res[('channel', 0)]['results']['pauli']['det'], color=GREEN, lw=2, ls='-.')
    ax.axhline(0.0, color=GREY, lw=0.8)
    ax.text(1.05, -1 / 27 - 0.015, r'target $-1/27$', fontsize=7, color=GREY, va='top')
    ax.set_ylim(-0.1, 0.42)
    ax.set_xscale('log', base=2)
    ax.set_xticks(depths)
    ax.set_xticklabels([str(L) for L in DEPTHS])
    ax.set_xlabel('depth $L$ (residual blocks)', fontsize=8, color=INK)
    ax.set_ylabel(r'$\det R_{\rm net}$', fontsize=8, color=INK)
    ax.set_title(r'(iii) determinant of the fit', fontsize=9, color=INK)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
        ax.spines[['left', 'bottom']].set_color(GREY)
        ax.tick_params(labelsize=7, colors=GREY, labelcolor=INK)
        ax.tick_params(which='minor', colors=GREY)
        ax.grid(color='#eeeeec', lw=0.6)
    h, lab = axes[2].get_legend_handles_labels()
    fig.legend(h, lab, loc='lower center', ncol=5, frameon=False, fontsize=7.5)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, metadata={'CreationDate': None})
    plt.close(fig)


def sci(x, digits=1):
    if x == 0:
        return '0'
    m, e = f'{x:.{digits}e}'.split('e')
    return f'{m}\\times10^{{{int(e)}}}'


def num(x):
    """Three decimals, or scientific notation below 0.01."""
    return f'{x:.3f}' if x >= 0.01 else sci(x)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--steps', type=int, default=2000)
    ap.add_argument('--restarts', type=int, default=6)
    ap.add_argument('--lr', type=float, default=0.05)
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--seed', type=int, default=0, help='seed of the initial parameters (per job: seed + 100 j)')
    ap.add_argument('--target-seed', type=int, default=7, help='seed of the random GKSL target')
    ap.add_argument('--out', default=str(OUT), help='directory containing generated/ and figures/')
    ap.add_argument('--replot', action='store_true', help='redraw the figure from the JSON of an earlier run')
    args = ap.parse_args()
    out = Path(args.out)
    if args.replot:
        saved = json.loads((out / 'generated' / 'residual_markov.json').read_text())
        figure({(d['kind'], d['depth']): d for d in saved['results']}, saved['pauli_references'],
               out / 'figures' / 'residual_markov.pdf')
        return
    t_all = time.time()
    specs = [('residual', L) for L in sorted(DEPTHS, reverse=True)] + [('lindblad', 0), ('channel', 0)]
    jobs = [(k, L, args.seed + 100 * j, args.steps, args.restarts, args.lr, args.target_seed)
            for j, (k, L) in enumerate(specs)]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        done = list(pool.map(train, jobs))
    res = {(d['kind'], d['depth']): d for d in done}
    seconds = time.time() - t_all
    refs = pauli_references()
    Rt = target_ptms(args.target_seed)
    tdet = {t: float(torch.linalg.det(Rt[i])) for i, t in enumerate(TARGETS)}
    maps, tinfo = target_maps(args.target_seed)
    wecc = {t: semigroup_test(maps[t]) for t in ('gksl', 'damped')}       # Pauli: degenerate spectrum, det < 0
    for (k, L), d in sorted(res.items()):
        print(f"{k:9s} L={L:2d} {d['seconds']:6.0f} s  " + '  '.join(
            f"{t}: hs={r['hs']:.2e} td={r['worst_td']:.2e} det={r['det']:+.3e} (min {r['det_min']:+.3e})"
            for t, r in d['results'].items()))
    print('semigroup test (best branch, min eig of w J w):', {t: v['best'] for t, v in wecc.items()})
    print('targets det', tdet, '\npauli refs', refs, f'\ntotal {seconds:.0f} s')

    # consistency: all residual and Lindblad maps have det > 0, the Pauli fits respect the proved bounds
    for (k, L), d in res.items():
        if k != 'channel':
            assert all(min(r['det_restarts']) > 0 for r in d['results'].values()), (k, L)
            assert d['results']['pauli']['hs'] > refs['hs_lower'] and d['results']['pauli']['worst_td_min'] > refs['td_lower'] - 1e-6

    figure(res, refs, out / 'figures' / 'residual_markov.pdf')
    gen_dir = out / 'generated'
    gen_dir.mkdir(parents=True, exist_ok=True)
    payload = dict(environment=environment('cpu'), args=vars(args), seconds=seconds, depths=DEPTHS, p=P_MAX,
                   targets=dict(info=tinfo, ptm={t: Rt[i].tolist() for i, t in enumerate(TARGETS)}, det=tdet),
                   pauli_references=refs, semigroup_test=wecc,
                   results=[dict(kind=d['kind'], depth=d['depth'], seconds=d['seconds'], results=d['results'],
                                 history=d['history']) for d in done])
    (gen_dir / 'residual_markov.json').write_text(json.dumps(payload, indent=1))

    R = lambda k, L, t, q='hs': res[(k, L)]['results'][t][q]
    names = dict(gksl='Gksl', damped='Damped', pauli='Pauli')
    deep = max(DEPTHS)
    tex = ['% generated by scripts/residual_markov.py -- do not edit by hand',
           f'\\newcommand{{\\RMSteps}}{{{args.steps}}}',
           f'\\newcommand{{\\RMRestarts}}{{{args.restarts}}}',
           f'\\newcommand{{\\RMp}}{{{P_MAX}}}',
           f'\\newcommand{{\\RMDepthMax}}{{{deep}}}',
           f'\\newcommand{{\\RMDepths}}{{{", ".join(map(str, DEPTHS))}}}',
           f'\\newcommand{{\\RMSeconds}}{{{seconds:.0f}}}',
           f'\\newcommand{{\\RMWorkers}}{{{args.workers}}}',
           f'\\newcommand{{\\RMMarkovA}}{{{refs["markov_a"]:.4f}}}',
           f'\\newcommand{{\\RMMarkovHS}}{{{refs["markov_hs"]:.3f}}}',
           f'\\newcommand{{\\RMMarkovTD}}{{{refs["markov_td"]:.3f}}}',
           f'\\newcommand{{\\RMDephasedHS}}{{{refs["dephased_hs"]:.3f}}}',
           f'\\newcommand{{\\RMSemigroupKmax}}{{{max(wecc["damped"]["branches"])}}}',
           f'\\newcommand{{\\RMDampedSemigroupViolation}}{{{-wecc["damped"]["best"]:.3f}}}']
    for t, N in names.items():
        tex += [f'\\newcommand{{\\RM{N}Channel}}{{{sci(R("channel", 0, t))}}}',
                f'\\newcommand{{\\RM{N}Lindblad}}{{{num(R("lindblad", 0, t))}}}',
                f'\\newcommand{{\\RM{N}ResOne}}{{{num(R("residual", 1, t))}}}',
                f'\\newcommand{{\\RM{N}ResDeep}}{{{num(R("residual", deep, t))}}}',
                f'\\newcommand{{\\RM{N}Det}}{{{tdet[t]:.4f}}}']
    pres = [R('residual', L, 'pauli') for L in DEPTHS]
    ptd = [R('residual', L, 'pauli', 'worst_td') for L in DEPTHS]
    tex += [f'\\newcommand{{\\RMPauliResMin}}{{{min(pres):.3f}}}',
            f'\\newcommand{{\\RMPauliResMinTD}}{{{min(ptd):.3f}}}',
            f'\\newcommand{{\\RMPauliResTDAll}}{{{min(R("residual", L, "pauli", "worst_td_min") for L in DEPTHS):.4f}}}',
            f'\\newcommand{{\\RMPauliResTDAllDepth}}{{{min(DEPTHS, key=lambda L: R("residual", L, "pauli", "worst_td_min"))}}}',
            f'\\newcommand{{\\RMPauliLindbladTD}}{{{R("lindblad", 0, "pauli", "worst_td"):.3f}}}',
            f'\\newcommand{{\\RMRestartSpread}}{{{sci(max(max(r["hs_restarts"]) - min(r["hs_restarts"]) for L in DEPTHS for r in res[("residual", L)]["results"].values()) + 1e-16)}}}',
            f'\\newcommand{{\\RMPauliResDetMin}}{{{sci(min(R("residual", L, "pauli", "det_min") for L in DEPTHS))}}}',
            f'\\newcommand{{\\RMPauliResDetMax}}{{{max(R("residual", L, "pauli", "det") for L in DEPTHS):.3f}}}',
            f'\\newcommand{{\\RMPauliLindbladDet}}{{{R("lindblad", 0, "pauli", "det"):.3f}}}',
            f'\\newcommand{{\\RMPauliChannelDet}}{{{R("channel", 0, "pauli", "det"):.4f}}}']
    (gen_dir / 'residual_markov.tex').write_text('\n'.join(tex) + '\n')


if __name__ == '__main__':
    main()
