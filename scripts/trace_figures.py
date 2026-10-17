#!/usr/bin/env python3
"""
Figures of man/trace_envelopes, grouped five inequalities at a time, one panel per dimension.

Reads the certified envelopes (man/article/generated/trace_characterizations.json, scripts/trace_characterizations.py),
the Monte Carlo summaries (man/trace_envelopes/generated/trace_montecarlo.json, scripts/trace_montecarlo.py) and the
two-contrast envelopes (man/trace_envelopes/generated/trace_contrast.json, scripts/trace_contrast.py), and writes to
man/trace_envelopes/figures:

  envelope_<group>.pdf   top row: V(delta) on log-log axes for M_2, M_3, M_4 (proved value: line, certified search:
                         markers, grey guides delta^k); bottom row: the best Monte Carlo sample in each bin of delta
                         divided by the envelope, for random operators (solid) and scaled projections (dashed);
  operator_<group>.pdf   V(delta, delta') against the contrast delta' of the operators, M_2 and M_3, delta = 0.9 (solid)
                         and 0.6 (dashed).

    python scripts/trace_figures.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
OUT = ROOT / 'man' / 'trace_envelopes'

COLORS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4']       # categorical slots 1-5 (validated palette)
MARKERS = ['o', 's', '^', 'D', 'v']
STYLES = ['-', '--', '-.', ':', (0, (5, 1, 1, 1))]
INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e6e5e1'

GROUPS = {
    'linear': ('Envelopes of order $\\delta$', ['young', 'projections', 'meet', 'join', 'lp']),
    'gardner': ('Gardner\'s inequality and its relatives (all five envelopes equal $\\kappa-1$)', ['gardner', 'triangle', 'realpart', 'square', 'metric']),
    'projection': ('Inequalities for two projections', ['pqp', 'absdiff', 'expqpq', 'repos', 'reabs']),
    'operator': ('Inequalities for positive operators', ['sgm', 'alt', 'det', 'subadd', 'powers']),
    'higher': ('Higher order and search-only items', ['monotone', 'holder', 'gt', 'pb', 'wmon']),
}
CONTRAST_GROUPS = {
    'smooth': ('Smooth growth in both contrasts', ['gardner', 'triangle', 'realpart', 'young']),
    'threshold': ('A threshold in the contrast of the operators', ['sgm', 'monotone', 'subadd', 'powers', 'wmon']),
    'other': ('A peak, or no dependence on the operators', ['square', 'alt', 'det', 'pb', 'lp']),
}
NAMES = dict(gardner='1 Gardner', triangle='2 triangle', realpart='3 $\\varphi(A+A^*)$', square='4 $|\\varphi(A^2)|$',
             metric='5 metric', sgm='6 geom. mean', alt='7 Araki-Lieb-Thirring', young='8 Young (ext.)',
             projections='9 Bikchentaev-Abed', monotone='10 monotonicity', det='11 AM-GM', subadd='12 subadditivity',
             powers='13 Powers-St\\o rmer', holder='14 Hoelder', gt='15 Golden-Thompson', pb='16 Peierls-Bogoliubov',
             pqp='17 $\\varphi(PQP)\\leq\\varphi(Q)$', absdiff='18 $\\varphi(|P-Q|)$', meet='19 $P\\wedge Q$',
             expqpq='20 $e^{QPQ}$', repos='21 Re$\\,PQ\\geq0$', reabs='22 Re$\\,PQ\\leq|PQ|$', join='23 $P\\vee Q$',
             lp='24 power mean', wmon='25 weighted monot.')


def style(ax):
    ax.spines[['top', 'right']].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=7)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def envelope_figures():
    import trace_characterizations as tc
    env = json.loads((ROOT / 'man' / 'article' / 'generated' / 'trace_characterizations.json').read_text())
    mc = json.loads((OUT / 'generated' / 'trace_montecarlo.json').read_text())
    grid = np.array(env['grid'], dtype=float)
    edges = np.array(mc['delta_bins'])
    mid = np.sqrt(np.maximum(edges[:-1], 1e-3) * edges[1:])
    upper = edges[1:]                                  # the envelope at the upper bin edge bounds every sample in the bin
    dims = [2, 3, 4]
    xs = np.logspace(-2, np.log10(0.95), 200)
    for key, (title, items) in GROUPS.items():
        fig, axes = plt.subplots(2, 3, figsize=(11, 6.4), gridspec_kw=dict(height_ratios=[1.35, 1]))
        lo, hi = np.inf, 0
        for j, n in enumerate(dims):
            ax, bx = axes[0, j], axes[1, j]
            for i, name in enumerate(items):
                c = COLORS[i]
                V = np.array(env['inequalities'][name]['dims'][str(n)]['V'], dtype=float)
                ok = (grid > 0) & (V > 1e-9)
                ax.plot(grid[ok], V[ok], MARKERS[i], color=c, ms=4.5, mfc='white', mew=1.2, zorder=3)
                ref, dashed = None, False
                if tc.INEQS[name]['exact'] is not None:
                    ref = tc.INEQS[name]['exact'](xs, n)
                if ref is None and tc.INEQS[name]['bound'] is not None:
                    ref, dashed = tc.INEQS[name]['bound'](xs, n), True
                if ref is not None:
                    ref = np.asarray(ref, dtype=float)
                    m = ref > 1e-14
                    ax.plot(xs[m], ref[m], ls=':' if dashed else '-', color=c, lw=1.4, zorder=2)
                if np.any(ok):
                    lo, hi = min(lo, V[ok].min()), max(hi, V[ok].max())
                ax.plot([], [], ls=':' if dashed else '-', marker=MARKERS[i], color=c, mfc='white', label=NAMES[name])
                # Monte Carlo: best sample per bin / envelope at the bin centre
                Vm = np.interp(upper, grid, V)
                b = mc['binned'][str(n)][name]
                for mode, ls in (('random', '-'), ('projection', '--')):
                    if mode not in b:
                        continue
                    mx = np.array([np.nan if x is None else x for x in b[mode]['delta']['max']], dtype=float)
                    r = np.where((Vm > 1e-9) & (upper <= grid[-1]), mx / np.maximum(Vm, 1e-300), np.nan)
                    r = np.where(r < 0, 0.0, r)
                    bx.plot(mid, r, ls=ls, color=c, lw=1.4, marker=MARKERS[i], ms=3)
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_title(f'$M_{n}$', fontsize=10, color=INK)
            bx.set_xscale('log')
            bx.set_ylim(-0.03, 1.08)
            bx.axhline(1, color=MUTED, lw=0.8)
            bx.set_xlabel(r'$\delta=\Delta(W)$', fontsize=8, color=MUTED)
            style(ax)
            style(bx)
        for ax in axes[0]:
            ax.set_ylim(max(lo / 3, 1e-12), hi * 3)
            for k in (1, 2, 4, 8):
                g = hi * (xs / xs[-1]) ** k              # guides through the upper right corner
                m = g > max(lo / 3, 1e-12)
                ax.plot(xs[m], g[m], color='#cfcec8', lw=0.7, zorder=1)
                if m.sum() > 5:
                    ax.text(xs[m][0], g[m][0], f'$\\delta^{k}$' if k > 1 else r'$\delta$', fontsize=7,
                            color=MUTED, ha='right', va='center')
        axes[0, 0].set_ylabel(r'envelope $V(\delta)$', fontsize=8, color=MUTED)
        axes[1, 0].set_ylabel('best Monte Carlo sample / $V$', fontsize=8, color=MUTED)
        h, l = axes[0, 0].get_legend_handles_labels()
        h += [plt.Line2D([], [], color=MUTED, ls='-'), plt.Line2D([], [], color=MUTED, ls='--')]
        l += ['random operators (bottom)', 'scaled projections (bottom)']
        fig.legend(h, l, loc='lower center', ncol=4, frameon=False, fontsize=8)
        fig.suptitle(title, fontsize=11, color=INK)
        fig.tight_layout(rect=(0, 0.08, 1, 0.97))
        fig.savefig(OUT / 'figures' / f'envelope_{key}.pdf', metadata={'CreationDate': None})
        plt.close(fig)


def operator_figures():
    d = json.loads((OUT / 'generated' / 'trace_contrast.json').read_text())
    cfg = d['config']
    tab = {}
    for r in d['B']:
        if r['kind'] == 'free':
            tab.setdefault((r['name'], r['n']), {})[r['vdelta']] = np.array(r['cert'], dtype=float)
    vg = np.array(cfg['vgrid'])
    gb = cfg['grid_b']
    for key, (title, items) in CONTRAST_GROUPS.items():
        fig, axes = plt.subplots(1, len(cfg['dims']), figsize=(10, 4.2), sharey=True)
        for ax, n in zip(axes, cfg['dims']):
            for i, name in enumerate(items):
                for g, ls in ((0.9, '-'), (0.6, '--')):
                    y = np.array([tab[(name, n)][x][gb.index(g)] for x in vg])
                    y = np.where(y > 1e-6, y, np.nan)
                    ax.plot(vg, y, ls=ls, color=COLORS[i], marker=MARKERS[i], ms=3.5, mfc='white', lw=1.4,
                            label=NAMES[name] if g == 0.9 else None)
            ax.set_yscale('log')
            ax.set_xlim(-0.02, 1.02)
            ax.set_title(f'$M_{n}$', fontsize=10, color=INK)
            ax.set_xlabel(r"contrast $\delta'$ of the operators", fontsize=8, color=MUTED)
            style(ax)
        axes[0].set_ylabel(r"$V(\delta,\delta')$", fontsize=8, color=MUTED)
        h, l = axes[0].get_legend_handles_labels()
        h += [plt.Line2D([], [], color=MUTED, ls='-'), plt.Line2D([], [], color=MUTED, ls='--')]
        l += [r'$\delta=0.9$', r'$\delta=0.6$']
        fig.legend(h, l, loc='lower center', ncol=4, frameon=False, fontsize=8)
        fig.suptitle(title, fontsize=11, color=INK)
        fig.tight_layout(rect=(0, 0.13, 1, 0.95))
        fig.savefig(OUT / 'figures' / f'operator_{key}.pdf', metadata={'CreationDate': None})
        plt.close(fig)


if __name__ == '__main__':
    (OUT / 'figures').mkdir(exist_ok=True)
    envelope_figures()
    operator_figures()
    print('wrote', ', '.join(sorted(p.name for p in (OUT / 'figures').glob('envelope_*.pdf'))), 'and',
          ', '.join(sorted(p.name for p in (OUT / 'figures').glob('operator_*.pdf'))))
