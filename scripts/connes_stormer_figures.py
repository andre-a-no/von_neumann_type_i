"""Figures for man/connes_stormer/connes_stormer_masas.tex.

Writes man/connes_stormer/figures/cs_bloch.pdf, cs_window.pdf, cs_asymptotics.pdf and cs_numerics.pdf.

The n = 2 figures are computed with NumPy/SciPy directly from the formulas of the paper (the lemma on psi, the
variational formula, the proposition on the optimal angle and the asymptotics of Theorem B), independently of the
library: psi(beta) = sin(2 beta) log tan^2 beta,
phi(alpha) = int_{alpha-theta}^{alpha} psi, the optimal angle alpha* is the root of psi(alpha) = psi(alpha - theta)
in the bracket given there, and H = phi(alpha*). The n >= 3 figure only plots the stored results of
scripts/connes_stormer.py (man/article/generated/connes_stormer.json); it does not rerun that computation.

Usage: python scripts/connes_stormer_figures.py
"""
import json
import os
import warnings

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'man', 'connes_stormer', 'figures')
DATA = os.path.join(ROOT, 'man', 'article', 'generated', 'connes_stormer.json')

BLUE, ORANGE, GREEN, YELLOW = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
GREY, INK, LIGHT = '#7d848c', '#0b0b0b', '#e6e8ea'

plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['cmr10', 'DejaVu Serif'], 'mathtext.fontset': 'cm',
    'axes.formatter.use_mathtext': True, 'font.size': 8.5, 'axes.titlesize': 9, 'axes.labelsize': 8.5,
    'legend.fontsize': 7.5, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
    'axes.edgecolor': GREY, 'axes.labelcolor': INK, 'xtick.color': GREY, 'ytick.color': GREY,
    'xtick.labelcolor': INK, 'ytick.labelcolor': INK, 'axes.linewidth': 0.6,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 3, 'ytick.major.size': 3,
    'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
    'lines.linewidth': 1.4, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02,
})
WIDTH = 6.3


# ---------------------------------------------------------------- formulas of the paper (n = 2)

def eta(t):
    t = np.asarray(t, dtype=float)
    return np.where(t > 0, -t * np.log(np.where(t > 0, t, 1.0)), 0.0)


def h2(p):
    """Binary entropy h_2(p) = eta(p) + eta(1 - p)."""
    return eta(p) + eta(1.0 - np.asarray(p, dtype=float))


def psi(b):
    """psi(beta) = sin(2 beta) log tan^2 beta, extended by 0 on (pi/2) Z."""
    b = np.asarray(b, dtype=float)
    s, c = np.sin(b), np.cos(b)
    ok = (np.abs(s) > 1e-300) & (np.abs(c) > 1e-300)
    with np.errstate(divide='ignore', invalid='ignore'):
        v = np.sin(2 * b) * np.log((s / np.where(ok, c, 1.0)) ** 2)
    return np.where(ok, v, 0.0)


def phi(a, th):
    """phi(alpha) = h_2(cos^2(alpha - theta)) - h_2(cos^2 alpha) = int_{alpha-theta}^alpha psi."""
    return h2(np.cos(a - th) ** 2) - h2(np.cos(a) ** 2)


Y = brentq(lambda y: y * np.arctanh(y) - 1.0, 0.1, 1 - 1e-12)
KAPPA = 2 * np.sqrt(1 - Y * Y) / Y            # max psi
BETA_K = -np.arccos(Y) / 2                     # the maximiser of psi in (-pi/4, 0)


def alpha_star(th):
    """The optimal angle alpha* of the paper, in (-pi/4, pi/4]."""
    if th <= 0 or th >= np.pi / 2:
        return np.nan
    if abs(th - np.pi / 4) < 1e-14:
        return 0.0
    if th > np.pi / 4:
        return -alpha_star(np.pi / 2 - th)
    lo, hi = max(BETA_K, th - np.pi / 4), min(BETA_K + th, 0.0)
    g = lambda a: psi(a) - psi(a - th)
    eps = 1e-14 * max(1.0, hi - lo)
    return brentq(g, lo + eps, hi - eps, xtol=1e-15)


def H2(th):
    """H(Delta_2 | U Delta_2 U^*) = phi(alpha*) for the angle theta, by quadrature of psi over the optimal window."""
    if th <= 0 or th >= np.pi / 2:
        return 0.0
    a = alpha_star(th)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')        # the log-singularity of psi' at 0 only limits the attainable rel. error
        return quad(lambda b: float(psi(b)), a - th, a, epsabs=1e-16, epsrel=1e-12, limit=200)[0]


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, name), metadata={'CreationDate': None})
    plt.close(fig)


def diverging():
    return LinearSegmentedColormap.from_list('bo', [BLUE, '#f1f2f3', ORANGE])


# ---------------------------------------------------------------- figure 1: Bloch circle and optimal angle

def fig_bloch(c=0.9):
    th = np.arccos(np.sqrt(c))
    a = alpha_star(th)
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 2.85), gridspec_kw=dict(width_ratios=[1.0, 1.25]))

    ax = axes[0]
    w = np.linspace(-np.pi, np.pi, 1441)
    val = phi(w / 2, th)                        # value of the basis with Bloch vector at angle omega from z
    pts = np.column_stack([np.sin(w), np.cos(w)])
    seg = np.stack([pts[:-1], pts[1:]], axis=1)
    m = np.abs(val).max()
    lc = LineCollection(seg, cmap=diverging(), norm=TwoSlopeNorm(0.0, -m, m), linewidths=7, capstyle='butt')
    lc.set_array(0.5 * (val[:-1] + val[1:]))
    ax.add_collection(lc)
    ax.plot([0], [0], 'o', ms=2, color=GREY)
    ax.plot([-1.25, 1.25], [0, 0], color=LIGHT, lw=0.6, zorder=0)
    ax.plot([0, 0], [-1.25, 1.25], color=LIGHT, lw=0.6, zorder=0)
    ax.text(1.27, 0, '$x$', va='center', color=GREY)
    ax.text(0, 1.29, '$z$', ha='center', color=GREY)

    def basis(omega, color, lw, ls='-'):
        x, z = np.sin(omega), np.cos(omega)
        ax.plot([-x, x], [-z, z], color=color, lw=lw, ls=ls, zorder=3)
        ax.plot([x, -x], [z, -z], 'o', ms=4.5, color=color, mec='white', mew=0.8, zorder=4)

    basis(0.0, INK, 1.0)
    basis(2 * th, BLUE, 1.0)
    basis(2 * a, ORANGE, 1.4)

    def arc(w0, w1, r, color):
        t = np.linspace(w0, w1, 60)
        ax.plot(r * np.sin(t), r * np.cos(t), color=color, lw=0.8)
    arc(0, 2 * th, 0.42, BLUE)
    arc(2 * a, 0, 0.30, ORANGE)
    ax.text(0.47 * np.sin(th), 0.47 * np.cos(th), r'$2\theta$', color=BLUE, fontsize=8, ha='left', va='bottom')
    ax.text(-0.13, 0.36, r'$2\alpha^*$', color=ORANGE, fontsize=8, ha='right', va='bottom')
    ax.text(0.05, 1.10, r'$\pm z$: $\Delta_2$, Choda', color=INK, fontsize=7.5, ha='left')
    ax.text(1.12 * np.sin(2 * th), 1.06 * np.cos(2 * th) - 0.08, r'$\pm q$: $Q$', color=BLUE, fontsize=7.5,
            ha='left')
    ax.text(1.10 * np.sin(2 * a), 1.10 * np.cos(2 * a) + 0.02, r'$\pm r^*$', color=ORANGE, fontsize=7.5,
            ha='right')
    ax.set_xlim(-1.35, 1.45)
    ax.set_ylim(-1.3, 1.4)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(rf'(a) the $x$-$z$ circle, $c={c}$', fontsize=8.5)
    cb = fig.colorbar(lc, ax=ax, fraction=0.05, pad=0.0, shrink=0.75)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=6.5, length=2)
    cb.set_label('value of the basis', fontsize=7)

    ax = axes[1]
    t = np.linspace(1e-4, np.pi / 2 - 1e-4, 1200)
    al = np.array([alpha_star(x) for x in t])
    left = t < np.pi / 4
    tl = t[left]
    lo = np.maximum(BETA_K, tl - np.pi / 4)
    hi = np.minimum(BETA_K + tl, 0.0)
    ax.fill_between(tl, lo, hi, color=ORANGE, alpha=0.13, lw=0, label='bracket of Prop. 4.6')
    ax.fill_between(np.pi / 2 - tl, -hi, -lo, color=ORANGE, alpha=0.13, lw=0)
    ax.plot(t[left], al[left], color=ORANGE, label=r'optimal angle $\alpha^*$')
    ax.plot(t[~left], al[~left], color=ORANGE)
    ax.plot(t, np.where(np.abs(t - np.pi / 4) < 2e-3, np.nan, np.where(t < np.pi / 4, t, t - np.pi / 2)),
            color=BLUE, lw=1.0, ls='--',
            label=r'basis of $Q$ ($\alpha=\theta$ mod $\frac{\pi}{2}$)')
    ax.axhline(0, color=INK, lw=0.8, label=r'basis of $\Delta_2$ (Choda)')
    ax.axhline(BETA_K, color=GREY, lw=0.6, ls=':')
    ax.axhline(-BETA_K, color=GREY, lw=0.6, ls=':')
    ax.text(1.0, BETA_K - 0.065, r'$\beta_\kappa$', color=GREY, ha='center', fontsize=7.5)
    ax.text(np.pi / 2 - 0.5, -BETA_K + 0.02, r'$-\beta_\kappa$', color=GREY, ha='center', fontsize=7.5)
    ax.axvline(th, color=GREY, lw=0.5, ls='-', alpha=0.5)
    ax.set_xlim(0, np.pi / 2)
    ax.set_ylim(-np.pi / 4, np.pi / 4)
    ax.set_xticks([0, np.pi / 8, np.pi / 4, 3 * np.pi / 8, np.pi / 2])
    ax.set_xticklabels(['0', r'$\pi/8$', r'$\pi/4$', r'$3\pi/8$', r'$\pi/2$'])
    ax.set_yticks([-np.pi / 4, -np.pi / 8, 0, np.pi / 8, np.pi / 4])
    ax.set_yticklabels([r'$-\pi/4$', r'$-\pi/8$', '0', r'$\pi/8$', r'$\pi/4$'])
    ax.set_xlabel(r'$\theta$  ($c=\cos^2\theta$)')
    ax.set_ylabel(r'angle $\alpha$ (radians)')
    ax.legend(loc='upper left', fontsize=6.5, handlelength=1.8, bbox_to_anchor=(0.0, 1.02))
    ax.set_title(r'(b) the optimal basis against $\theta$', fontsize=8.5)
    fig.tight_layout(w_pad=1.0)
    save(fig, 'cs_bloch.pdf')
    return th, a


# ---------------------------------------------------------------- figure 2: psi and the sliding window

def fig_window(c=0.9):
    th = np.arccos(np.sqrt(c))
    a = alpha_star(th)
    H = H2(th)
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 2.45))

    ax = axes[0]
    b = np.linspace(-np.pi / 2, np.pi / 2, 2001)
    ax.axhline(0, color=GREY, lw=0.5)
    w = np.linspace(a - th, a, 300)
    ax.fill_between(w, 0, psi(w), color=ORANGE, alpha=0.35, lw=0)
    w0 = np.linspace(-th, 0, 300)
    ax.fill_between(w0, 0, psi(w0), facecolor='none', edgecolor=BLUE, hatch='//////', lw=0)
    ax.plot(b, psi(b), color=INK, lw=1.2)
    ax.plot([BETA_K], [KAPPA], 'o', ms=3.5, color=INK)
    ax.annotate(r'$[\alpha^*-\theta,\alpha^*]$: area $H$', xy=(a - th / 2, 0.55), xytext=(-1.45, 1.9),
                fontsize=7.5, color=ORANGE, arrowprops=dict(arrowstyle='-', color=ORANGE, lw=0.6))
    ax.annotate(r'$[-\theta,0]$: area $h(U)$', xy=(-0.03, 0.35), xytext=(0.05, 1.95),
                fontsize=7.5, color=BLUE, arrowprops=dict(arrowstyle='-', color=BLUE, lw=0.6))
    ax.set_xlim(-np.pi / 2, np.pi / 2)
    ax.set_ylim(-1.9, 2.3)
    ax.set_xticks([-np.pi / 2, -np.pi / 4, 0, np.pi / 4, np.pi / 2])
    ax.set_xticklabels([r'$-\pi/2$', r'$-\pi/4$', '0', r'$\pi/4$', r'$\pi/2$'])
    ax.set_xlabel(r'$\beta$')
    ax.set_ylabel(r'$\psi(\beta)=\sin 2\beta\,\log\tan^2\beta$')
    ax.set_title(rf'(a) windows of length $\theta$, $c={c}$', fontsize=8.5)

    ax = axes[1]
    al = np.linspace(-np.pi / 2, np.pi / 2, 2001)
    ax.axhline(0, color=GREY, lw=0.5)
    ax.axhline(KAPPA * th, color=GREY, lw=0.7, ls='--')
    ax.text(-0.95, KAPPA * th + 0.015, r'$\kappa\theta$', color=GREY, ha='center', va='bottom', fontsize=7.5)
    ax.plot(al, phi(al, th), color=INK, lw=1.2)
    hc = float(h2(c))
    ax.plot([0], [hc], 'o', ms=4.5, color=BLUE, mec='white', mew=0.8, zorder=4)
    ax.plot([th], [-hc], 'o', ms=4.5, color=BLUE, mec='white', mew=0.8, zorder=4)
    ax.plot([a, a + np.pi / 2, a - np.pi / 2], [H] * 3, 'o', ms=5, color=ORANGE, mec='white', mew=0.8, zorder=4)
    ax.annotate(r'$\alpha=0$: $h(U)$', xy=(0.02, hc - 0.01), xytext=(0.35, 0.25), color=BLUE, fontsize=7.5,
                arrowprops=dict(arrowstyle='-', color=BLUE, lw=0.6))
    ax.annotate(r'$\alpha=\theta$: $-h(U)$', xy=(th + 0.02, -hc - 0.01), xytext=(0.75, -0.36), color=BLUE,
                fontsize=7.5, arrowprops=dict(arrowstyle='-', color=BLUE, lw=0.6))
    ax.text(a, H + 0.04, r'$\alpha^*$: $H$', color=ORANGE, fontsize=7.5, ha='center', va='bottom')
    ax.set_xlim(-np.pi / 2, np.pi / 2)
    ax.set_ylim(-0.5, 0.62)
    ax.set_xticks([-np.pi / 2, -np.pi / 4, 0, np.pi / 4, np.pi / 2])
    ax.set_xticklabels([r'$-\pi/2$', r'$-\pi/4$', '0', r'$\pi/4$', r'$\pi/2$'])
    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'$\varphi(\alpha)$ (nats)')
    ax.set_title(r'(b) the value $\varphi(\alpha)$ of the basis at angle $\alpha$', fontsize=8.5)
    fig.tight_layout(w_pad=1.5)
    save(fig, 'cs_window.pdf')
    return H, hc


# ---------------------------------------------------------------- figure 3: H against Choda's value

def fig_asymptotics():
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 2.3), gridspec_kw=dict(width_ratios=[1.15, 1.0, 1.0]))
    t = np.linspace(1e-4, np.pi / 2 - 1e-4, 600)
    H = np.array([H2(x) for x in t])
    hc = h2(np.cos(t) ** 2)
    ax = axes[0]
    ax.fill_between(t, hc, H, color=ORANGE, alpha=0.18, lw=0)
    ax.plot(t, np.minimum(KAPPA * np.minimum(t, np.pi / 2 - t), np.log(2) + 0.6), color=GREY, lw=0.8, ls='--')
    ax.axhline(np.log(2), color=GREY, lw=0.6, ls=':')
    ax.plot(t, H, color=ORANGE, label=r'$H(\Delta_2|U\Delta_2U^*)$')
    ax.plot(t, hc, color=BLUE, label=r'$h(U)$ (Choda)')
    ax.text(0.10, 0.42, r'$\kappa\theta$', color=GREY, fontsize=7.5)
    ax.text(np.pi / 4, np.log(2) + 0.03, r'$\log 2$', color=GREY, fontsize=7.5, ha='center')
    ax.set_xlim(0, np.pi / 2)
    ax.set_ylim(0, 0.85)
    ax.set_xticks([0, np.pi / 8, np.pi / 4, 3 * np.pi / 8, np.pi / 2])
    ax.set_xticklabels(['0', r'$\pi/8$', r'$\pi/4$', r'$3\pi/8$', r'$\pi/2$'])
    ax.set_xlabel(r'$\theta$')
    ax.set_ylabel('nats')
    ax.legend(loc='lower center', fontsize=6.5)
    ax.set_title('(a) $H$ and $h(U)$', fontsize=8.5)

    ax = axes[1]
    s = np.geomspace(1e-3, np.pi / 4 * 0.98, 160)
    Hs = np.array([H2(x) for x in s])
    hs = h2(np.cos(s) ** 2)
    ax.loglog(s, KAPPA * s, color=GREY, lw=0.8, ls='--')
    ax.loglog(s, 2 * s ** 2 * np.log(1 / s) + s ** 2, color=GREY, lw=0.8, ls=':')
    ax.loglog(s, Hs, color=ORANGE, label='$H$')
    ax.loglog(s, hs, color=BLUE, label='$h(U)$')
    ax.text(8e-3, 3e-2, r'$\kappa\theta$', color=GREY, fontsize=7.5)
    ax.text(1.3e-2, 2e-5, r'$2\theta^2\log\frac{1}{\theta}+\theta^2$', color=GREY, fontsize=7.5)
    ax.set_xlabel(r'$\theta$')
    ax.legend(loc='upper left', fontsize=6.5)
    ax.set_title('(b) small angles', fontsize=8.5)

    ax = axes[2]
    s = np.geomspace(5e-3, np.pi / 4 * 0.98, 120)
    Hs = np.array([H2(x) for x in s])
    lim = (KAPPA + 4 / KAPPA) / 6
    ax.axhline(lim, color=GREY, lw=0.8, ls='--')
    ax.semilogx(s, (KAPPA * s - Hs) / s ** 3, color=ORANGE)
    ax.text(6e-3, lim + 0.02, r'$(\kappa+4/\kappa)/6$', color=GREY, fontsize=7.5, va='bottom')
    ax.set_ylim(0.66, 0.76)
    ax.set_xlabel(r'$\theta$')
    ax.set_ylabel(r'$(\kappa\theta-H)/\theta^3$')
    ax.set_title('(c) the cubic term', fontsize=8.5)
    fig.tight_layout(w_pad=0.8)
    save(fig, 'cs_asymptotics.pdf')


# ---------------------------------------------------------------- figure 4: stored results for n = 3, 4, 5

def fig_numerics():
    with open(DATA) as f:
        res = json.load(f)
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 2.45))
    colors = {3: BLUE, 4: ORANGE, 5: GREEN}
    markers = {3: 'o', 4: 's', 5: 'D'}
    ax = axes[0]
    for n in (3, 4, 5):
        r = res[f'n{n}']
        H = np.maximum(r['H_general'], r['H_bases'])
        ch = np.array(r['choda'])
        ax.plot(ch / np.log(n), H - ch, markers[n], ms=2.4, color=colors[n], alpha=0.75, mew=0,
                label=f'$n={n}$')
    ax.set_xlabel(r'$h(U)/\log n$')
    ax.set_ylabel(r'$H-h(U)$ (nats)')
    ax.set_ylim(bottom=0)
    ax.legend(loc='upper right', markerscale=2, fontsize=7)
    ax.set_title(r"(a) excess over Choda's value, Haar-random $U$", fontsize=8.5)

    ax = axes[1]
    rng = np.random.default_rng(0)
    for k, n in enumerate((3, 4, 5)):
        gaps = np.array([cs['certified_gap'] for cs in res[f'n{n}']['cases']])
        x = k + rng.uniform(-0.12, 0.12, size=gaps.size)
        ax.semilogy(x, gaps, markers[n], ms=4.2, color=colors[n], mec='white', mew=0.5)
        ax.text(k, 2.2e-2, f"{gaps.size}/{res[f'n{n}']['nU']}", ha='center', fontsize=7, color=INK)
    ax.axhline(res['n3']['flag'], color=GREY, lw=0.7, ls='--')
    ax.text(2.45, res['n3']['flag'] * 1.6, 'flag threshold', color=GREY, fontsize=7, ha='right')
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['$n=3$', '$n=4$', '$n=5$'])
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim(3e-9, 6e-2)
    ax.set_ylabel('partition $-$ best basis (nats)')
    ax.set_title('(b) unitaries where a partition beats the bases', fontsize=8.5)
    fig.tight_layout(w_pad=1.5)
    save(fig, 'cs_numerics.pdf')


def main():
    th, a = fig_bloch()
    H, hc = fig_window()
    fig_asymptotics()
    fig_numerics()
    print(f'kappa = {KAPPA:.10f}, beta_kappa = {BETA_K:.10f}; c = 0.9: theta = {th:.6f}, alpha* = {a:.6f}, '
          f'H = {H:.6f}, h = {hc:.6f}')
    for c in (0.99, 0.95, 0.9, 0.8, 0.7, 0.6):
        t = np.arccos(np.sqrt(c))
        print(f'c = {c}: alpha* = {alpha_star(t):+.4f}, H = {H2(t):.6f}')
    print('written:', OUT)


if __name__ == '__main__':
    main()
