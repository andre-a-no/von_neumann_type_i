#!/usr/bin/env python3
"""
Schematic figures of the library paper (man/article/vn_algebra_pytorch.tex) and the summary figure of the noisy
XY mixer, drawn from data that already exist (no computation of the library is run):

  figures/architecture.pdf   user code -> TypeIAlgebra / Operator (generator, flags, caches) -> size-group
                             dispatcher -> torch.linalg, sparse and gather kernels on CPU/GPU
  figures/layout.pdf         the padded (B, C, k_max, k_max) tensor and the size groups, contrasted with one
                             block-diagonal D x D matrix
  figures/ladder.pdf         the sector ladder of the noisy register: N = 0..n, S^- (N -> N-1), S^+ (N -> N+1),
                             S^z inside the sectors, the post-selection projection P_w
  figures/noisy_summary.pdf  the integral I and the gain factor of post-selection against n, and the infidelities
                             for the largest register, from generated/noisy_xy.json (scripts/noisy_xy_mixer.py)

    python scripts/library_paper_schematics.py
"""
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                     # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle   # noqa: E402
import numpy as np                                                  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / 'man' / 'article'
FIG = ART / 'figures'

BLUE, ORANGE, GREEN, YELLOW = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
GREY, INK = '#7d848c', '#0b0b0b'
LIGHT = '#eef0f2'                      # light fill for boxes and padding
META = {'CreationDate': None}

plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'mathtext.fontset': 'cm', 'axes.edgecolor': GREY,
                     'axes.labelcolor': INK, 'xtick.color': GREY, 'ytick.color': GREY, 'text.color': INK,
                     'axes.linewidth': 0.6, 'pdf.fonttype': 42})


def box(ax, x, y, w, h, title, body=None, edge=GREY, fill='white', title_color=INK, fs=8, lw=1.0):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=0.12', lw=lw, ec=edge, fc=fill))
    if body is None:
        ax.text(x + w / 2, y + h / 2, title, ha='center', va='center', fontsize=fs, color=title_color)
    else:
        ax.text(x + w / 2, y + h - 0.16, title, ha='center', va='top', fontsize=fs, color=title_color,
                fontweight='bold')
        ax.text(x + w / 2, y + h - 0.48, body, ha='center', va='top', fontsize=fs - 1.2, color=INK,
                linespacing=1.35)


def arrow(ax, p, q, color=GREY, lw=1.0, style='-|>', rad=0.0, ms=8, ls='-'):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle=style, mutation_scale=ms, color=color, lw=lw,
                                 connectionstyle=f'arc3,rad={rad}', linestyle=ls, shrinkA=0, shrinkB=0))


# ------------------------------------------------------------------------------------------------ architecture
def architecture():
    fig, ax = plt.subplots(figsize=(6.8, 3.7))
    ax.set_xlim(0, 13.6)
    ax.set_ylim(0, 7.4)
    ax.axis('off')
    c1, w1, c2, w2, c3, w3 = 1.45, 4.55, 6.3, 3.45, 10.05, 3.5
    for y, lab in ((6.75, 'user code'), (5.05, 'data types'), (3.0, 'dispatch'), (0.95, 'kernels')):
        ax.text(0.0, y, lab, ha='left', va='center', fontsize=7.5, color=GREY, style='italic')
    box(ax, c1, 6.35, w1, 0.8, 'user scripts and examples', fs=7.5)
    box(ax, c2, 6.35, w2 + 0.3 + w3, 0.8, 'models: SpinChain, SU2Chain, Inclusion, ...', fs=7.5)
    box(ax, c1, 4.05, w1, 1.95, 'TypeIAlgebra',
        'sector sizes $k_c$, field, precision,\ndevice; cached constants:\nmasks, trace weights $w_c$, size groups',
        edge=BLUE, title_color=BLUE)
    box(ax, c2, 4.05, w2 + 0.3 + w3, 1.95, 'Operator and its subclasses',
        'matrix $(B,C,k_{\\max},k_{\\max})$, or a generator (lazy recipe)\n'
        'flags: self-adjoint, normal, positive, invertible, ...\n'
        'caches: $\\lambda_{\\min}$, $\\lambda_{\\max}$, trace, inverse, $|A|$, spectrum',
        edge=BLUE, title_color=BLUE)
    box(ax, c1, 2.25, w1, 1.45, 'size groups',
        'stack the sectors of equal size;\none batched call per size', edge=ORANGE, title_color=ORANGE)
    box(ax, c2, 2.25, w2, 1.45, 'kernel by flags',
        'eigvalsh, eigh, dilation,\ncontour integral', edge=ORANGE, title_color=ORANGE)
    box(ax, c3, 2.25, w3, 1.45, 'cost module',
        'time and memory\nestimates, warnings', edge=ORANGE, title_color=ORANGE)
    box(ax, c1, 0.25, w1, 1.45, 'torch.linalg',
        'LAPACK (CPU);\ncuSOLVER, cuBLAS, MAGMA (GPU)', edge=GREEN, title_color=GREEN)
    box(ax, c2, 0.25, w2, 1.45, 'sparse and gathers',
        'CSR Lanczos, Krylov;\nindex maps (SectorJump)', edge=GREEN, title_color=GREEN)
    box(ax, c3, 0.25, w3, 1.45, 'autograd',
        'Daleckii-Krein\nbackward passes', edge=GREEN, title_color=GREEN)
    arrow(ax, (c1 + w1 / 2, 6.35), (c1 + w1 / 2, 6.0))
    arrow(ax, (c2 + 3.4, 6.35), (c2 + 3.4, 6.0))
    arrow(ax, (c2 + 0.6, 6.35), (c1 + w1 - 0.4, 6.0))
    arrow(ax, (c1 + w1, 5.0), (c2, 5.0), style='<|-|>')
    for x in (c1 + w1 / 2, c2 + w2 / 2, c3 + w3 / 2):
        arrow(ax, (x, 4.05), (x, 3.7))
    arrow(ax, (c1 + w1, 2.97), (c2, 2.97), style='<|-|>')
    for x in (c1 + w1 / 2, c2 + w2 / 2):
        arrow(ax, (x, 2.25), (x, 1.7))
    arrow(ax, (c2 + w2 * 0.8, 2.25), (c3 + 0.3, 1.7), color=GREY)
    ax.text(c2 + w2 / 2 + 0.15, 3.87, 'on first read', ha='left', va='center', fontsize=6.5, color=GREY,
            style='italic')
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / 'architecture.pdf', bbox_inches='tight', metadata=META)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------ data layout
def layout():
    ks = [2, 3, 3, 4]                       # block sizes k_c
    C, kmax, D = len(ks), max(ks), sum(ks)
    colors = [BLUE, ORANGE, ORANGE, GREEN]  # equal sizes share a colour (one size group)
    fig, axes = plt.subplots(1, 3, figsize=(6.8, 2.5), gridspec_kw=dict(width_ratios=[1.0, 1.55, 1.2]))
    for ax in axes:
        ax.set_aspect('equal')
        ax.axis('off')
    # (a) one block-diagonal D x D matrix
    ax = axes[0]
    ax.add_patch(Rectangle((0, 0), D, D, fc='white', ec=GREY, lw=0.8))
    off = 0
    for k, col in zip(ks, colors):
        ax.add_patch(Rectangle((off, D - off - k), k, k, fc=col, ec='white', lw=1.0, alpha=0.85))
        off += k
    for i in range(1, D):
        ax.plot([i, i], [0, D], color=LIGHT, lw=0.4, zorder=0)
        ax.plot([0, D], [i, i], color=LIGHT, lw=0.4, zorder=0)
    ax.set_xlim(-0.6, D + 0.6)
    ax.set_ylim(-3.2, D + 2.4)
    ax.text(D / 2, D + 0.6, '(a) one matrix in $M_D$, $D=%d$' % D, ha='center', fontsize=7.5)
    ax.text(D / 2, -0.7, '%d of %d entries are\nstructural zeros' % (D * D - sum(k * k for k in ks), D * D),
            ha='center', va='top', fontsize=7, color=GREY)
    # (b) the padded tensor (B, C, k_max, k_max): C slices, a few batch layers behind
    ax = axes[1]
    gap, depth = 0.9, 3
    for b in reversed(range(depth)):
        d = 0.35 * b
        for c, (k, col) in enumerate(zip(ks, colors)):
            x0 = c * (kmax + gap) + d
            ax.add_patch(Rectangle((x0, d), kmax, kmax, fc='white', ec=GREY, lw=0.6, zorder=-b))
            if b == 0:
                ax.add_patch(Rectangle((x0, 0), kmax, kmax, fc='none', ec=GREY, lw=0.0, hatch='////',
                                       zorder=0.1, alpha=0.35))
                ax.add_patch(Rectangle((x0, kmax - k), k, k, fc=col, ec='white', lw=1.0, alpha=0.85, zorder=1))
                ax.text(x0 + kmax / 2, -0.45, '$c=%d$\n$k_c=%d$' % (c + 1, k), ha='center', va='top', fontsize=6.8)
    W = C * (kmax + gap) - gap
    ax.set_xlim(-0.5, W + 1.6)
    ax.set_ylim(-5.0, kmax + 3.3)
    ax.text(W / 2, kmax + 1.7, r'(b) library: tensor $(B,\,C,\,k_{\max},\,k_{\max})$', ha='center', fontsize=7.5)
    ax.text(W + 1.25, kmax + 0.95, '$B$', fontsize=7.5, color=GREY, ha='center')
    arrow(ax, (W + 0.1, kmax + 0.1), (W + 0.95, kmax + 0.95), color=GREY, ms=6, lw=0.8)
    ax.text(W / 2, -3.3, 'hatched: padding (zeros)', ha='center', va='top', fontsize=7, color=GREY)
    # (c) size groups: one batched call per distinct block size
    ax = axes[2]
    groups = [(2, [1], BLUE), (3, [2, 3], ORANGE), (4, [4], GREEN)]
    y = 0.0
    for k, cs, col in reversed(groups):
        for j in range(3):                           # a few stacked matrices of the batch
            ax.add_patch(Rectangle((0.3 * j, y + 0.3 * j), k * 0.55, k * 0.55, fc=col, ec='white', lw=0.8,
                                   alpha=0.85, zorder=-j))
        lab = '$k=%d$: sector%s %s' % (k, 's' if len(cs) > 1 else '', ', '.join(str(c) for c in cs))
        mult = '' if len(cs) == 1 else str(len(cs))
        ax.text(3.3, y + 0.55 * k / 2 + 0.25, lab + '\n' + r'eigh on $(%sB,%d,%d)$' % (mult, k, k), va='center',
                fontsize=6.8)
        y += 0.55 * k + 1.2
    ax.set_xlim(-0.2, 9.0)
    ax.set_ylim(-0.6, y + 0.9)
    ax.text(4.3, y + 0.2, '(c) one batched call per size', ha='center', fontsize=7.5)
    fig.tight_layout(w_pad=0.3)
    fig.savefig(FIG / 'layout.pdf', bbox_inches='tight', metadata=META)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------ sector ladder
def ladder(n=6, w=3):
    fig, ax = plt.subplots(figsize=(6.6, 2.4))
    ax.axis('off')
    step = 2.0
    hmax = 1.5
    dims = [math.comb(n, N) for N in range(n + 1)]
    ys = {}
    for N, d in enumerate(dims):
        x = N * step
        h = 0.35 + hmax * math.sqrt(d / max(dims))       # height grows with the block size
        fill = '#fdf3dc' if N == w else LIGHT
        edge = YELLOW if N == w else GREY
        ax.add_patch(FancyBboxPatch((x - 0.55, -h / 2), 1.1, h, boxstyle='round,pad=0,rounding_size=0.08',
                                    fc=fill, ec=edge, lw=1.6 if N == w else 0.8))
        ax.text(x, 0.0, '$N=%d$\n$M_{%d}$' % (N, d), ha='center', va='center', fontsize=7.5)
        ys[N] = h / 2
    for N in range(n):                                   # S^+ above (N -> N+1), S^- below (N+1 -> N)
        x0, x1 = N * step, (N + 1) * step
        arrow(ax, (x0 + 0.3, ys[N] + 0.02), (x1 - 0.3, ys[N + 1] + 0.02), color=BLUE, rad=-0.45, lw=1.1)
        arrow(ax, (x1 - 0.3, -ys[N + 1] - 0.02), (x0 + 0.3, -ys[N] - 0.02), color=ORANGE, rad=-0.45, lw=1.1)
    for N in range(n + 1):                               # S^z and the mixer act inside every sector
        x = N * step
        ax.add_patch(FancyArrowPatch((x + 0.56, 0.14), (x + 0.56, -0.14), arrowstyle='-|>', mutation_scale=6,
                                     color=GREEN, lw=1.0, connectionstyle='arc3,rad=-1.6'))
    top = max(ys.values())
    yl = top + 0.85
    ax.plot([w * step, w * step], [yl + 0.1, top + 0.03], color=YELLOW, lw=1.4)
    ax.text(w * step + 0.15, yl + 0.1, 'post-selection $P_w$ on weight $w=%d$' % w, ha='left', va='center',
            fontsize=7.5)
    ax.plot([-0.55, -0.1], [yl + 0.1] * 2, color=BLUE, lw=1.6)
    ax.text(0.0, yl + 0.1, r'$S^+_i$: $N\to N+1$, rate $\gamma_\uparrow(n-N)$', fontsize=7.5, va='center')
    ax.plot([-0.55, -0.1], [-yl - 0.1] * 2, color=ORANGE, lw=1.6)
    ax.text(0.0, -yl - 0.1, r'$S^-_i$: $N\to N-1$, rate $\gamma_\downarrow N$', fontsize=7.5, va='center')
    ax.plot([w * step + 0.15, w * step + 0.6], [-yl - 0.1] * 2, color=GREEN, lw=1.6)
    ax.text(w * step + 0.7, -yl - 0.1, r'$S^z_i$ (dephasing) and the mixer $H$: inside each sector', fontsize=7.5,
            va='center')
    ax.set_xlim(-0.8, n * step + 0.9)
    ax.set_ylim(-yl - 0.5, yl + 0.5)
    ax.set_aspect('equal')
    fig.savefig(FIG / 'ladder.pdf', bbox_inches='tight', metadata=META)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------ noisy summary
def noisy_summary():
    data = json.loads((ART / 'generated' / 'noisy_xy.json').read_text())
    summ = data['summary']
    ns = sorted(int(n) for n in summ)
    rates = np.array(data['rates'])
    I = np.array([summ[str(n)]['I'] for n in ns])
    T = summ[str(ns[0])]['T']
    gd, gu, gp = data['models']['mixed']
    gain = []                                 # infidelity without / with post-selection at the smallest rate
    for n in ns:
        r = data['results'][str(n)]['mixed']
        gain.append((1 - r['F'][0]) / (1 - r['Fpost'][0]))
    gain = np.array(gain)
    w = np.array([n // 2 for n in ns])
    gain_first = 1 + (gd * w + gu * (np.array(ns) - w)) * T / (gp * I)     # first-order prediction

    fig, axes = plt.subplots(1, 3, figsize=(6.8, 2.75), gridspec_kw=dict(width_ratios=[1, 1, 1.35]))

    def style(ax):
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(labelsize=7)
        ax.grid(color='#e6e8ea', lw=0.5)
        ax.set_axisbelow(True)

    ax = axes[0]
    ax.plot(ns, I, '-o', color=ORANGE, lw=1.6, ms=4.5, mfc='white', mew=1.4)
    ax.set_xticks(ns)
    ax.set_ylim(0, 1.7)
    ax.set_xlabel('qubits $n$', fontsize=8)
    ax.set_title('(a) dephasing integral $I$', fontsize=8)
    style(ax)

    ax = axes[1]
    ax.plot(ns, gain_first, ':', color=INK, lw=1.0)
    ax.plot(ns, gain, '-o', color=GREEN, lw=1.6, ms=4.5, mfc='white', mew=1.4)
    ax.text(ns[1], 8.6, 'computed (circles)\nfirst order (dotted)', fontsize=6.6, color=INK, ha='center', va='center')
    ax.set_xticks(ns)
    ax.set_ylim(0, 10)
    ax.set_xlabel('qubits $n$', fontsize=8)
    ax.set_title(r'(b) gain $(1-F)/(1-F_{\rm post})$', fontsize=8)
    style(ax)

    ax = axes[2]
    n = ns[-1]
    res = data['results'][str(n)]
    ax.plot(rates, 1 - np.array(res['mixed']['F']), '--', color=GREEN, lw=1.2)
    for m, col in (('thermal', BLUE), ('dephasing', ORANGE), ('mixed', GREEN)):
        ax.plot(rates, np.maximum(1 - np.array(res[m]['Fpost']), 1e-16), '-o', color=col, lw=1.6, ms=3.5,
                mfc='white', mew=1.2)
    ax.plot(rates, rates * summ[str(n)]['I'], ':', color=INK, lw=1.0)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(1e-7, 1.5)
    ax.set_xlim(rates[0] / 1.4, rates[-1] * 1.4)
    ax.set_xlabel(r'rate $\gamma$', fontsize=8)
    ax.set_title('(c) infidelity, $n=%d$' % n, fontsize=8)
    handles = [plt.Line2D([], [], color=col, lw=1.6, label=lab) for col, lab in
               ((BLUE, r'thermal $T_1$, $\gamma_\uparrow=0.2\gamma_\downarrow$'), (ORANGE, r'dephasing $\gamma_\varphi$'),
                (GREEN, r'mixed, $\gamma_\downarrow=\gamma_\varphi$, $\gamma_\uparrow=0.1\gamma_\downarrow$'))]
    handles += [plt.Line2D([], [], color=GREY, ls='--', lw=1.2, label='without post-selection'),
                plt.Line2D([], [], color=GREY, ls='-', marker='o', mfc='white', ms=3.5, lw=1.6,
                           label='post-selected on the weight'),
                plt.Line2D([], [], color=INK, ls=':', lw=1.0, label=r'first order $\gamma_\varphi I$')]
    fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, fontsize=6.8, handlelength=2.4,
               columnspacing=1.6)
    style(ax)
    fig.tight_layout(w_pad=1.0, rect=(0, 0.17, 1, 1))
    fig.savefig(FIG / 'noisy_summary.pdf', bbox_inches='tight', metadata=META)
    plt.close(fig)
    for n_, g_, f_ in zip(ns, gain, gain_first):
        print(f'n={n_}: gain at gamma={rates[0]:g}: {g_:.2f} (first order {f_:.2f})')


def main():
    architecture()
    layout()
    ladder()
    noisy_summary()
    print('wrote', ', '.join(str(FIG / f) for f in ('architecture.pdf', 'layout.pdf', 'ladder.pdf',
                                                     'noisy_summary.pdf')))


if __name__ == '__main__':
    main()
