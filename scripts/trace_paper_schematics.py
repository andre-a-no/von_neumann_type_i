"""Schematic and summary figures for man/trace_envelopes/trace_inequalities_contrast.tex.

Writes three PDFs to man/trace_envelopes/figures/:

  ladder_of_orders.pdf   computed envelopes V(delta) on M_2 of all 25 items on log-log axes, grouped by the
                         order of vanishing at the trace, with guide slopes 1, 2, 4, 8
                         (data: man/article/generated/trace_characterizations.json);
  operator_angle.pdf     the operator angle (a ball of radius Delta(A)|x| about x seen from the origin), the
                         Jordan product (angles add) and the Gardner quotient as the secant of arg(1 + r e^{i alpha}),
                         r = delta delta';
  contrast_square.pdf    the contrast square (delta, delta') with the level curves of the Gardner envelope, the
                         proved no-violation band of item 12, the exact violation thresholds of item 24 and the
                         numerically observed thresholds (data: man/trace_envelopes/generated/trace_contrast.json);
                         and the plane (Delta(A), Delta(B)) with the Jordan disc and the larger item-12 region.

Plain matplotlib/NumPy; runs in a few seconds:  python scripts/trace_paper_schematics.py
"""
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc, Circle, Wedge

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'man' / 'trace_envelopes' / 'figures'
DATA_ENV = ROOT / 'man' / 'article' / 'generated' / 'trace_characterizations.json'
DATA_TC = ROOT / 'man' / 'trace_envelopes' / 'generated' / 'trace_contrast.json'

BLUE, ORANGE, GREEN, YELLOW = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
GREY, INK, MUTED, GRID = '#7d848c', '#0b0b0b', '#52514e', '#e6e5e1'

plt.rcParams.update({
    'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 9, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'axes.edgecolor': MUTED, 'axes.labelcolor': MUTED, 'xtick.color': MUTED, 'ytick.color': MUTED,
    'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False, 'legend.fontsize': 7,
    'font.family': 'serif', 'font.serif': ['cmr10', 'DejaVu Serif'], 'axes.formatter.use_mathtext': True,
    'axes.unicode_minus': False, 'mathtext.fontset': 'cm', 'pdf.fonttype': 42, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02,
})

# order of vanishing of the envelope at the trace (Table 'classes' of the paper; items 12-16 by the search)
ORDER = {1: [8, 9, 19, 23, 24], 2: [1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22],
         4: [10], 8: [25]}
ORDER_COLOR = {1: BLUE, 2: ORANGE, 4: GREEN, 8: YELLOW}


def ladder_of_orders():
    d = json.loads(DATA_ENV.read_text())
    grid = np.array(d['grid'])
    names = list(d['inequalities'])                     # items 1..25 in table order
    fig, ax = plt.subplots(figsize=(4.6, 3.3))
    for k, items in ORDER.items():
        for i in items:
            v = np.array([np.nan if x is None else x for x in d['inequalities'][names[i - 1]]['dims']['2']['V']])
            keep = (grid > 0) & (v > 1e-8)
            if i == 25:                                 # the search resolves item 25 only from delta = 0.4 on
                keep &= grid >= 0.4
            ax.plot(grid[keep], v[keep], '-', color=ORDER_COLOR[k], lw=0.9, alpha=0.75, marker='o', ms=2.2)
    # guide slopes, anchored a factor below a representative item of each order
    anchors = {1: (19, 0.4, (0.01, 0.6)), 2: (17, 0.25, (0.01, 0.6)), 4: (10, 0.25, (0.05, 0.8)),
               8: (25, 0.25, (0.4, 0.95))}
    for k, (item, fac, (lo, hi)) in anchors.items():
        v = np.array(d['inequalities'][names[item - 1]]['dims']['2']['V'])
        j = int(np.argmin(abs(grid - (0.5 if k < 8 else 0.8))))
        c = fac * v[j] / grid[j] ** k
        xs = np.geomspace(lo, hi, 50)
        ax.plot(xs, c * xs ** k, '--', color=GREY, lw=0.8)
        ax.text(xs[0] * 0.93, c * xs[0] ** k, rf'$\delta^{{{k}}}$' if k > 1 else r'$\delta$', color=GREY,
                ha='right', va='center', fontsize=8)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlim(0.006, 1.05); ax.set_ylim(1e-9, 300)
    ax.set_xlabel(r'contrast $\delta=\Delta(W)$ of the density')
    ax.set_ylabel(r'computed envelope $V_2(\delta)$')
    ax.grid(True, which='major', color=GRID, lw=0.5)
    handles = [plt.Line2D([], [], color=ORDER_COLOR[k], marker='o', ms=2.5, lw=0.9,
                          label={1: r'order $\delta$: items 8, 9, 19, 23, 24',
                                 2: r'order $\delta^2$: 18 items',
                                 4: r'order $\delta^4$: item 10 (monotonicity)',
                                 8: r'order $\delta^8$: item 25 (weighted monotonicity)'}[k]) for k in ORDER]
    ax.legend(handles=handles, loc='upper left', handlelength=1.6)
    fig.savefig(OUT / 'ladder_of_orders.pdf')
    plt.close(fig)


def _arrow(ax, p, q, color, lw=1.2):
    ax.annotate('', xy=q, xytext=p, arrowprops=dict(arrowstyle='-|>', color=color, lw=lw, shrinkA=0, shrinkB=0,
                                                   mutation_scale=8))


def _label(ax, x, y, text, color, **kw):
    ax.text(x, y, text, color=color, fontsize=kw.pop('fontsize', 8), **kw)


def operator_angle():
    fig, axes = plt.subplots(1, 3, figsize=(6.6, 2.45))
    for ax in axes:
        ax.set_aspect('equal'); ax.axis('off')

    # (a) the angle of a positive operator: A x / lambda-bar(A) lies in the ball of radius Delta(A)|x| about x
    ax = axes[0]
    dA = 0.5
    psi = math.asin(dA)
    ax.add_patch(Circle((1, 0), dA, facecolor=BLUE, alpha=0.10, edgecolor=BLUE, lw=0.8))
    for s in (1, -1):
        ax.plot([0, 1.25 * math.cos(psi)], [0, s * 1.25 * math.sin(psi)], color=GREY, lw=0.7, ls='--')
    _arrow(ax, (0, 0), (1, 0), INK)
    _label(ax, 1.03, -0.13, r'$x$', INK)
    q = (0.86, 0.40)
    _arrow(ax, (0, 0), q, BLUE)
    _label(ax, 0.20, 0.52, r'$Ax/\bar\lambda(A)$', BLUE)
    ax.add_patch(Arc((0, 0), 0.8, 0.8, theta1=0, theta2=math.degrees(psi), color=MUTED, lw=0.8))
    _label(ax, 0.43, 0.06, r'$\psi_A$', MUTED, fontsize=7)
    ax.plot([1, 1], [0, -dA], color=BLUE, lw=0.6)
    _label(ax, 1.05, -0.36, r'$\Delta(A)\,|x|$', BLUE, fontsize=7)
    ax.set_xlim(-0.1, 1.65); ax.set_ylim(-0.8, 0.85)

    # (b) Jordan product: B turns x by at most psi_B, A turns Bx by at most psi_A
    ax = axes[1]
    pB, pA = math.radians(26), math.radians(32)
    tot = pA + pB
    ax.add_patch(Wedge((0, 0), 1.2, -math.degrees(tot), math.degrees(tot), facecolor=GREEN, alpha=0.10,
                       edgecolor='none'))
    for s in (1, -1):
        ax.plot([0, 1.2 * math.cos(tot)], [0, s * 1.2 * math.sin(tot)], color=GREEN, lw=0.7, ls='--')
    ax.plot([0, 0], [-1.05, 1.05], color=GREY, lw=0.6, ls=':')
    _arrow(ax, (0, 0), (1, 0), INK)
    _label(ax, 1.03, -0.13, r'$x$', INK)
    bx = (0.95 * math.cos(pB), 0.95 * math.sin(pB))
    _arrow(ax, (0, 0), bx, BLUE)
    _label(ax, bx[0] + 0.04, bx[1] - 0.02, r'$Bx$', BLUE)
    abx = (0.92 * math.cos(tot), 0.92 * math.sin(tot))
    _arrow(ax, (0, 0), abx, ORANGE)
    _label(ax, abx[0] - 0.12, abx[1] + 0.06, r'$ABx$', ORANGE)
    bax = (0.85 * math.cos(-0.55 * tot), 0.85 * math.sin(-0.55 * tot))
    _arrow(ax, (0, 0), bax, ORANGE)
    _label(ax, bax[0] + 0.03, bax[1] - 0.12, r'$BAx$', ORANGE)
    ax.add_patch(Arc((0, 0), 0.55, 0.55, theta1=0, theta2=math.degrees(pB), color=BLUE, lw=0.8))
    ax.add_patch(Arc((0, 0), 0.85, 0.85, theta1=math.degrees(pB), theta2=math.degrees(tot), color=ORANGE, lw=0.8))
    _label(ax, 0.30, 0.03, r'$\psi_B$', BLUE, fontsize=7)
    _label(ax, 0.40, 0.33, r'$\psi_A$', ORANGE, fontsize=7)
    _label(ax, 0.62, -1.08, r'cone of half-angle $\psi_A+\psi_B<\pi/2$', GREEN, fontsize=7, ha='center')
    ax.set_xlim(-0.25, 1.45); ax.set_ylim(-1.15, 1.1)

    # (c) Gardner's quotient on a two-dimensional block: |z| / Re z for z = 1 + delta delta' e^{i alpha}
    ax = axes[2]
    r = 0.6
    th = math.asin(r)
    ax.plot([-0.1, 1.75], [0, 0], color=GREY, lw=0.6); ax.plot([0, 0], [-0.85, 0.85], color=GREY, lw=0.6)
    ax.add_patch(Circle((1, 0), r, facecolor=ORANGE, alpha=0.06, edgecolor=ORANGE, lw=1.0))
    for s in (1, -1):
        ax.plot([0, 1.15 * math.cos(th)], [0, s * 1.15 * math.sin(th)], color=GREY, lw=0.7, ls='--')
    alpha = 2.35
    z = (1 + r * math.cos(alpha), r * math.sin(alpha))
    _arrow(ax, (0, 0), z, ORANGE)
    ax.plot([1, z[0]], [0, z[1]], color=ORANGE, lw=0.6)
    _label(ax, 0.05, 0.72, r'$z=1+\delta\delta^\prime e^{i\alpha}$', ORANGE, fontsize=7.5)
    _label(ax, 1.12, 0.14, r'$\delta\delta^\prime$', ORANGE, fontsize=7)
    ax.plot([1], [0], 'o', color=INK, ms=2)
    _label(ax, 1.0, -0.06, r'$1$', INK, fontsize=7, ha='center', va='top')
    ax.add_patch(Arc((0, 0), 0.62, 0.62, theta1=0, theta2=math.degrees(th), color=MUTED, lw=0.8))
    _label(ax, 0.33, 0.035, r'$\theta$', MUTED, fontsize=7)
    _label(ax, 0.82, -0.86, r'$\sin\theta=\delta\delta^\prime$,  $|z|/\mathrm{Re}\,z\leq\sec\theta$', INK,
           fontsize=7, ha='center')
    ax.set_xlim(-0.15, 1.75); ax.set_ylim(-0.95, 0.95)
    fig.subplots_adjust(wspace=0.06, top=0.86)
    for ax, title in zip(axes, (r'(a) $\sin\psi_A=\Delta(A)$', r'(b) $\langle (AB+BA)x,x\rangle>0$',
                                r'(c) Gardner: $|z|/\mathrm{Re}\,z$')):
        box = ax.get_position()
        fig.text(0.5 * (box.x0 + box.x1), 0.95, title, ha='center', va='top', fontsize=9, color=INK)
    fig.savefig(OUT / 'operator_angle.pdf')
    plt.close(fig)


def item24_quotient(n, d, dp):
    """Exact sup of phi(A^2)^{1/2}/phi(A) - 1 for phi(1)=n, Delta(W)=d, Delta(A)=dp (the proposition of Sect. 9)."""
    l, lp = (1 - d) / (1 + d), (1 - dp) / (1 + dp)
    xlo = n * lp + (1 - lp) * n * l / (n - 1 + l)
    x = max(xlo, 2 * n * lp / (1 + lp))
    return math.sqrt((1 + lp) * x - n * lp) / x - 1


def item24_threshold(n, d):
    lo, hi = 0.0, 1 - 1e-12
    if item24_quotient(n, d, hi) <= 0:
        return 1.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if item24_quotient(n, d, mid) <= 0 else (lo, mid)
    return 0.5 * (lo + hi)


def contrast_square():
    tc = json.loads(DATA_TC.read_text())
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 3.1), gridspec_kw={'wspace': 0.30})
    box = dict(facecolor='white', edgecolor='none', alpha=0.85, pad=0.6)

    # (a) the contrast square
    ax = axes[0]
    g = np.linspace(0, 1, 401)
    D, Dp = np.meshgrid(g, g)
    prod = np.clip(D * Dp, 0, 0.999)
    cs = ax.contour(D, Dp, 1 / np.sqrt(1 - prod ** 2) - 1, levels=[0.02, 0.1, 0.3, 1.0], colors=BLUE,
                    linewidths=0.7)
    ax.clabel(cs, fmt=lambda v: f'{v:g}', fontsize=6, inline_spacing=2,
              manual=[(0.30, 0.66), (0.50, 0.85), (0.83, 0.97), (0.985, 0.995)])
    c12 = 2 * math.sqrt(2) / 3
    ax.fill_between([0, 1], 0, c12, color=GREEN, alpha=0.09, lw=0)
    ax.axhline(c12, color=GREEN, lw=1.0)
    _label(ax, 0.02, 0.10, 'item 12: no violation below\n' + r'$\delta^\prime=2\sqrt{2}/3$ (proved)', GREEN,
           fontsize=6.5, bbox=box)
    ds = np.linspace(0.001, 0.999, 300)
    for n, ls in ((2, '-'), (3, '--')):
        ax.plot(ds, [item24_threshold(n, d) for d in ds], color=ORANGE, lw=1.0, ls=ls)
    _label(ax, 0.55, 0.655, r'item 24, $n=2$', ORANGE, fontsize=6.5, bbox=box)
    _label(ax, 0.55, 0.765, r'item 24, $n=3$', ORANGE, fontsize=6.5, bbox=box)
    # numerically observed onset of item-12 violations on M_2 (last grid value without, first with a violation)
    rows = [r for r in tc['B'] if r['name'] == 'subadd' and r['n'] == 2 and r['kind'] == 'free']
    for j, d in enumerate(rows[0]['grid']):
        vals = sorted((r['vdelta'], r['cert'][j]) for r in rows)
        first = next(vp for vp, c in vals if c > 1e-4)
        last = max(vp for vp, c in vals if vp < first)
        ax.plot([d, d], [last, first], color=INK, lw=2.0, solid_capstyle='butt', clip_on=False)
    ax.plot([], [], color=INK, lw=2.0, label=r'item 12: observed onset on $M_2$')
    ax.plot([], [], color=BLUE, lw=0.7, label=r'Gardner: $(1-\delta^2\delta^{\prime 2})^{-1/2}-1$')
    ax.plot([], [], color=ORANGE, lw=1.0, label='item 24: exact threshold')
    ax.legend(loc='center left', bbox_to_anchor=(0.0, 0.40), fontsize=6.3, facecolor='white', framealpha=0.85,
              frameon=True, edgecolor='none')
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.0)
    ax.set_xlabel(r'contrast $\delta$ of the density')
    ax.set_ylabel(r"contrast $\delta^\prime$ of the operators")
    ax.set_title('(a) the contrast square', color=INK, pad=9)

    # (b) the plane of the two operator contrasts
    ax = axes[1]
    t = np.linspace(0, math.pi / 2, 300)
    u = np.linspace(0, 1, 400)                          # u = Delta(A^{1/2}), v = Delta(B^{1/2}) = sqrt(1 - u^2)
    v = np.sqrt(1 - u ** 2)
    up, vp = 2 * u / (1 + u ** 2), 2 * v / (1 + v ** 2)  # Delta(A) = tanh(2 artanh u)
    ax.fill_between(up, 0, vp, color=GREEN, alpha=0.09, lw=0)
    ax.plot(up, vp, color=GREEN, lw=1.0)
    ax.fill_between(np.cos(t), 0, np.sin(t), color=BLUE, alpha=0.10, lw=0)
    ax.plot(np.cos(t), np.sin(t), color=BLUE, lw=1.0)
    ax.plot([0, c12, c12], [c12, c12, 0], color=GREEN, lw=0.6, ls=':')
    ax.plot([c12], [c12], 'o', color=GREEN, ms=3)
    _label(ax, 0.14, 0.34, r'$AB+BA\geq0$' + '\nfor all positions', BLUE, fontsize=6.5)
    _label(ax, 0.40, 0.84, r'$(A+B)^{1/2}\leq A^{1/2}+B^{1/2}$', GREEN, fontsize=6.5, ha='center', bbox=box)
    ax.set_xticks([0, 0.2, 0.4, 0.6, c12]); ax.set_xticklabels(['0', '0.2', '0.4', '0.6', r'$\frac{2\sqrt{2}}{3}$'])
    ax.set_yticks([0, 0.2, 0.4, 0.6, c12]); ax.set_yticklabels(['0', '0.2', '0.4', '0.6', r'$\frac{2\sqrt{2}}{3}$'])
    ax.set_xlim(0, 1.02); ax.set_ylim(0, 1.02); ax.set_aspect('equal')
    ax.set_xlabel(r'$\Delta(A)$'); ax.set_ylabel(r'$\Delta(B)$')
    ax.set_title('(b) the two operator contrasts', color=INK)
    fig.savefig(OUT / 'contrast_square.pdf')
    plt.close(fig)


if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    ladder_of_orders()
    operator_angle()
    contrast_square()
    print('wrote', ', '.join(str(OUT / f) for f in ('ladder_of_orders.pdf', 'operator_angle.pdf',
                                                     'contrast_square.pdf')))
