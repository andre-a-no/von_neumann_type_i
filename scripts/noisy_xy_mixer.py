#!/usr/bin/env python3
"""
Error detection by a conserved charge in a noisy XY mixer (use case of man/article, Section "Noisy circuits with a
conserved charge").

The XY mixer H = (1/2) sum_<ij> (S^+_i S^-_j + S^-_i S^+_j) on a ring of n qubits keeps the Hamming weight, so the
register carries the algebra (+)_N M_binom(n, N). Noise acts by sparse jumps (torch_vn_algebra.SectorJump):
amplitude damping S^-_i (rate g_down, sector N -> N-1), thermal excitation S^+_i (g_up, N -> N+1) and dephasing S^z_i
(g_phi, inside the sectors). Reading out the weight measures a central projection; post-selection on the initial
weight w keeps the accepted part. For a batch of noise models and rates (one batch element each) the script computes in
one Lindblad run per n (CPTP splitting, method='kraus'):

  acceptance p = Tr(P_w rho(T)),  raw fidelity F = <psi|rho(T)|psi>,  post-selected fidelity F_post = F / p,

with psi the noiseless state, and compares with the exact statements of the paper:
  T1 at zero temperature: p = exp(-w g T) and F_post = 1 (the no-jump generator differs from H by a central term);
  thermal T1: 1 - F_post = O(g_down g_up T^2) (only an up and a down jump together return to sector w);
  dephasing: 1 - F_post = g_phi I + O(g_phi^2), I = int_0^T sum_i (1/4 - <S^z_i>_t^2) dt on the ideal trajectory.
Writes man/article/figures/noisy_xy.pdf and man/article/generated/noisy_xy.{json,tex}.

    python scripts/noisy_xy_mixer.py [--dims 4,6,8,10] [--T 2.0] [--substeps 40]
"""
import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
import torch

from torch_vn_algebra import SpinChain, dynamics

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'man' / 'article'
RATES = np.logspace(-3, -1, 7)
MODELS = {                                   # (g_down, g_up, g_phi) per unit rate g
    'T1': (1.0, 0.0, 0.0),
    'thermal': (1.0, 0.2, 0.0),
    'dephasing': (0.0, 0.0, 1.0),
    'mixed': (1.0, 0.1, 1.0),
}
COLORS = {'thermal': '#2a78d6', 'dephasing': '#eb6834', 'mixed': '#1baf7a', 'T1': '#52514e'}
LABELS = {'T1': r'$T_1$, zero temperature', 'thermal': r'$T_1$, $\gamma_\uparrow=0.2\gamma_\downarrow$',
          'dephasing': r'dephasing $\gamma_\varphi$', 'mixed': r'$\gamma_\downarrow=\gamma_\varphi$, $\gamma_\uparrow=0.1\gamma_\downarrow$'}


def run(n, T, substeps):
    w = n // 2
    q = SpinChain(n, 'periodic', precision='double')
    H = q.xxz(J=1.0, Delta=0.0)
    rho0 = q.basis_state('1' * w + '0' * (n - w))
    # ideal trajectory: final state and the dephasing integral I
    ts = np.linspace(0.0, T, 201)
    traj = dynamics.von_neumann(rho0, H, ts.tolist())
    ideal = traj[-1].matrix[0]
    sz = [q.sz(i).matrix[0] for i in range(n)]
    var = [sum(0.25 - float((r.matrix[0] * s.transpose(-2, -1)).sum().real) ** 2 for s in sz) for r in traj]
    var = np.array(var)
    I = float(np.sum(0.5 * (var[1:] + var[:-1]) * np.diff(ts)))         # trapezoidal rule
    # batch: every (model, rate)
    combos = [(m, g) for m in MODELS for g in RATES]
    coef = {k: torch.tensor([math.sqrt(MODELS[m][j] * g) for m, g in combos], dtype=torch.float64)
            for j, k in enumerate(('lower', 'raise', 'sz'))}
    jumps = [q.sector_jump(i, k).scaled(coef[k]) for k in ('lower', 'raise', 'sz') for i in range(n)]
    t0 = time.time()
    rho = dynamics.lindblad_evolve(rho0, H, jumps, [0.0, T], substeps=substeps, method='kraus')[-1].matrix
    seconds = time.time() - t0
    p = rho[:, w].diagonal(dim1=-2, dim2=-1).sum(-1).real                      # weight-w probability
    F = (rho * ideal.transpose(-2, -1)).sum((-3, -2, -1)).real                # Tr(rho ideal)
    Fp = F / p
    res = {m: dict(rate=[], p=[], F=[], Fpost=[]) for m in MODELS}
    for (m, g), a, b, c in zip(combos, p.tolist(), F.tolist(), Fp.tolist()):
        res[m]['rate'].append(g)
        res[m]['p'].append(a)
        res[m]['F'].append(b)
        res[m]['Fpost'].append(c)
    k = q.algebra.k_factors
    info = dict(n=n, w=w, T=T, I=I, seconds=seconds, sector=math.comb(n, w), algebra=sum(x * x for x in k),
                padded=len(k) * max(k) ** 2, dense=4 ** n, substeps=substeps)
    return info, res


def figure(results):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ns = sorted(results)
    fig, axes = plt.subplots(1, len(ns), figsize=(10, 3.4), sharey=True)
    for ax, n in zip(np.atleast_1d(axes), ns):
        info, res = results[n]
        g = np.array(RATES)
        for m in ('thermal', 'dephasing', 'mixed'):
            r = res[m]
            ax.plot(g, 1 - np.array(r['F']), '--', color=COLORS[m], lw=1.2)
            ax.plot(g, np.maximum(1 - np.array(r['Fpost']), 1e-16), '-o', color=COLORS[m], lw=1.5, ms=3.5,
                    mfc='white')
        ax.plot(g, g * info['I'], ':', color='#0b0b0b', lw=1.0)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(1e-8, 1.5)
        ax.set_title(f'$n={n}$, weight {info["w"]}', fontsize=9)
        ax.set_xlabel(r'rate $\gamma$ ($T=' + f'{info["T"]:g}' + '$)', fontsize=8, color='#52514e')
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(labelsize=7, colors='#52514e')
        ax.grid(color='#e6e5e1', lw=0.6)
    np.atleast_1d(axes)[0].set_ylabel('infidelity', fontsize=8, color='#52514e')
    handles = [plt.Line2D([], [], color=COLORS[m], lw=1.5, label=LABELS[m]) for m in ('thermal', 'dephasing', 'mixed')]
    handles += [plt.Line2D([], [], color='#52514e', ls='--', label='without post-selection'),
                plt.Line2D([], [], color='#52514e', ls='-', marker='o', mfc='white', label='post-selected on the weight'),
                plt.Line2D([], [], color='#0b0b0b', ls=':', label=r'$\gamma_\varphi I$ (first order)')]
    fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.16, 1, 1))
    (OUT / 'figures').mkdir(exist_ok=True)
    fig.savefig(OUT / 'figures' / 'noisy_xy.pdf', metadata={'CreationDate': None})
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dims', default='4,6,8,10')
    ap.add_argument('--T', type=float, default=2.0)
    ap.add_argument('--substeps', type=int, default=40)
    args = ap.parse_args()
    torch.set_num_threads(max(1, torch.get_num_threads()))
    results = {}
    rows, summary = [], {}
    for n in [int(x) for x in args.dims.split(',')]:
        info, res = run(n, args.T, args.substeps)
        results[n] = (info, res)
        w, T = info['w'], info['T']
        g = np.array(RATES)
        t1 = res['T1']
        acc_err = float(np.max(np.abs(np.array(t1['p']) - np.exp(-w * g * T)) / np.exp(-w * g * T)))
        t1_post = float(np.max(np.abs(1 - np.array(t1['Fpost']))))
        th = 1 - np.array(res['thermal']['Fpost'])
        slope = float(np.polyfit(np.log(g[:3]), np.log(th[:3]), 1)[0])
        c2 = float(th[0] / g[0] ** 2)
        de = 1 - np.array(res['dephasing']['Fpost'])
        ratio = float(de[0] / (g[0] * info['I']))
        mixed_raw, mixed_post = 1 - res['mixed']['F'][4], 1 - res['mixed']['Fpost'][4]   # g = 0.01
        summary[n] = dict(info, acc_err=acc_err, t1_post=t1_post, thermal_slope=slope, thermal_c2=c2,
                          deph_ratio=ratio, mixed_raw=mixed_raw, mixed_post=mixed_post)
        print(f"n={n}: {info['seconds']:.1f} s; T1: max rel. err. of p {acc_err:.1e}, |1-F_post| {t1_post:.1e}; "
              f"thermal slope {slope:.3f}; dephasing (1-F_post)/(g I) {ratio:.4f}; mixed g=0.01 raw {mixed_raw:.4f} "
              f"post {mixed_post:.4f}")
        rows.append(f"{n} & {info['sector']} & {info['algebra']} & {info['dense']} & {info['seconds']:.0f} & "
                    f"{acc_err:.0e} & {t1_post:.0e} & {slope:.2f} & {info['I']:.3f} & {ratio:.3f} & "
                    f"{mixed_raw:.4f} & {mixed_post:.4f}\\\\".replace('e-0', 'e-'))
    figure(results)
    gen = OUT / 'generated'
    gen.mkdir(exist_ok=True)
    (gen / 'noisy_xy.json').write_text(json.dumps(dict(rates=RATES.tolist(), models=MODELS, summary=summary,
                                                       results={n: r for n, (i, r) in results.items()}), indent=1))
    worst_acc = max(s['acc_err'] for s in summary.values())
    worst_t1 = max(s['t1_post'] for s in summary.values())
    tex = ['% generated by scripts/noisy_xy_mixer.py -- do not edit by hand',
           '\\newcommand{\\NXYRows}{%', *rows, '}',
           f'\\newcommand{{\\NXYT}}{{{args.T:g}}}',
           f'\\newcommand{{\\NXYSubsteps}}{{{args.substeps}}}',
           f'\\newcommand{{\\NXYAccErr}}{{{worst_acc:.0e}}}'.replace('e-0', 'e-'),
           f'\\newcommand{{\\NXYTOnePost}}{{{worst_t1:.0e}}}'.replace('e-0', 'e-'),
           f'\\newcommand{{\\NXYBatch}}{{{len(MODELS) * len(RATES)}}}']
    (gen / 'noisy_xy.tex').write_text('\n'.join(tex) + '\n')
    print('wrote figures/noisy_xy.pdf and generated/noisy_xy.{json,tex}')


if __name__ == '__main__':
    main()
