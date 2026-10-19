"""Figure for Theorem C of man/trace_envelopes (Gardner's inequality with two prescribed contrasts).

Plain NumPy, independent of the library: random samples of z = |Tr(XUY)| - Tr(XY) in M_k(R) with ||X|| = ||Y|| = 1
(Haar-orthogonal eigenbases and U, spectra with a uniformly drawn contrast and the other eigenvalues uniform in
between), the sharp bound nu d^2 d'^2 / ((1+d)(1+d')) with the other contrast free in [0, 1] (nu = floor(k/2), the
paper's notation), and the value of the explicit maximisers of the proof, evaluated numerically; it prints the largest
difference between the two (quoted in the caption of the figure). Writes
man/trace_envelopes/figures/gardner_contrasts.pdf.

    python scripts/trace_envelopes_figures.py [--samples 100000] [--seed 0]
"""
import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--samples', type=int, default=100000)
p.add_argument('--seed', type=int, default=0)
p.add_argument('--out', default=str(ROOT / 'man' / 'trace_envelopes' / 'figures' / 'gardner_contrasts.pdf'))
args = p.parse_args()
rng = np.random.default_rng(args.seed)


def haar_orthogonal(m, k):
    q, r = np.linalg.qr(rng.standard_normal((m, k, k)))
    return q * np.sign(np.diagonal(r, axis1=1, axis2=2))[:, None, :]


def spectra(m, k):
    """lambda_max = 1, contrast d uniform in [0, 1), the other eigenvalues uniform in [ell, 1]."""
    d = rng.uniform(0, 1, m)
    ell = (1 - d) / (1 + d)
    lam = ell[:, None] + (1 - ell[:, None]) * rng.uniform(0, 1, (m, k))
    lam[:, 0], lam[:, 1] = 1.0, ell
    return lam, d


def samples(k, m):
    lx, dx = spectra(m, k)
    ly, dy = spectra(m, k)
    ox, oy, u = haar_orthogonal(m, k), haar_orthogonal(m, k), haar_orthogonal(m, k)
    X = ox @ (lx[:, :, None] * ox.transpose(0, 2, 1))
    Y = oy @ (ly[:, :, None] * oy.transpose(0, 2, 1))
    z = np.abs(np.einsum('nij,njk,nki->n', X, u, Y)) - np.einsum('nij,nji->n', X, Y)
    return dx, dy, z


def bound(N, d, dp):
    return N * d**2 * dp**2 / ((1 + d) * (1 + dp))


def maximiser_value(N, d, dp):
    """z at the maximiser of the proof: N blocks X = diag(1, l), Y with spectrum {1, l'} at the angle of the
    2x2 proposition, U the polar unitary of YX (adjoint)."""
    l, lp = (1 - d) / (1 + d), (1 - dp) / (1 + dp)
    S, q = l + lp, (1 + l) * (1 + lp)
    g = 0.25 + S / (2 * q)
    c, s = np.sqrt(g), np.sqrt(1 - g)
    R = np.array([[c, -s], [s, c]])
    X = np.diag([1.0, l])
    Y = R @ np.diag([1.0, lp]) @ R.T
    w, _, vh = np.linalg.svd(Y @ X)
    U = (w @ vh).T                                           # Tr(Y X U) = ||YX||_1
    return N * (abs(np.trace(X @ U @ Y)) - np.trace(X @ Y))


grid = np.linspace(0, 0.99, 100)
marks = np.linspace(0.1, 0.9, 9)
fig, axes = plt.subplots(2, 2, figsize=(9, 6.2), sharex=True)
for row, k in enumerate((2, 16)):
    N = k // 2
    dx, dy, z = samples(k, args.samples)
    for col, (dd, lab) in enumerate(((dx, r'$\Delta(X)$'), (dy, r'$\Delta(Y)$'))):
        ax = axes[row, col]
        ax.scatter(dd, z, s=0.8, color='#a9b4bf', alpha=0.15, lw=0, rasterized=True)
        ax.plot([], [], 'o', ms=3, color='#a9b4bf', label=f'random samples ({len(z)})')
        if col == 0:
            ax.plot(grid, bound(N, grid, 1.0), color='#1b1b1b', lw=1.6, label=r'bound, $\Delta(Y)$ free (Thm. C)')
            vals = [maximiser_value(N, d, 1.0) for d in marks]
            ax.plot(marks, vals, 'o', ms=5, mfc='#c8402f', mec='white', mew=0.8, label='explicit maximiser')
        else:
            ax.plot(grid, bound(N, 1.0, grid), color='#1b1b1b', lw=1.6, ls='--',
                    label=r'supremum, $\Delta(X)\to1$ (not attained)')
        ax.axhline(0, color='#555555', lw=0.6)
        ax.set_title(f'$M_{{{k}}}(\\mathbb{{R}})$, $\\nu={N}$', fontsize=9)
        if row == 1:
            ax.set_xlabel(lab)
        if col == 0:
            ax.set_ylabel(r'$z=|\mathrm{Tr}(XUY)|-\mathrm{Tr}(XY)$')
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='y', color='#e5e5e5', lw=0.5)
        ax.legend(frameon=False, fontsize=7, loc='lower right')
    print(f'k={k}: largest sample z = {z.max():+.4f}, bound at Delta(X)->1: {bound(N, 1, 1):.4f}; '
          f'max |maximiser - bound| = {max(abs(maximiser_value(N, d, 1.0) - bound(N, d, 1.0)) for d in marks):.1e}')
fig.tight_layout()
Path(args.out).parent.mkdir(parents=True, exist_ok=True)
fig.savefig(args.out, metadata={'CreationDate': None})
print('wrote', args.out)
