#!/usr/bin/env python3
"""
Sharp bounds instead of Monte Carlo: for the three trace functionals of the original experiments

    Exp 1  z = |Tr(X U Y)| - Tr(X Y)          X, Y >= 0, U unitary in M
    Exp 2  z = Tr|X Y| - Tr(|X| |Y|),  X = U X0, Y = V Y0   (z = ||X0 V Y0||_1 - Tr(X0 Y0))
    Exp 3  z = Tr(Y |X| Y) - Tr|Y X Y|         X = X*, Y >= 0

we estimate the upper envelope  sup { z : Delta(X) = delta }  and the lower envelope
inf { z : Delta(X) = delta } as functions of delta (the contrast of the other operator is
optimised as well), and likewise as functions of Delta(Y). Operators are normalised as in the
Monte Carlo experiments (lambda_max = 1, real algebra, blunt trace). Every grid point is
optimised from several random starts in one batch (torch_vn_algebra.optimize).

If the Monte Carlo samples of scripts/experiment.py are available in results/experiments, they are
overlaid and the script checks that all of them lie between the two envelopes.

Writes generated/bounds.tex, generated/bounds.json, figures/bounds_*.pdf.

    python scripts/paper/bounds_search.py --mode full --device cuda
"""
import math
import time
from pathlib import Path

import numpy as np
import torch

from common import ROOT, parse_args, save_json, write_tex, env_macro, sync

from torch_vn_algebra import TypeIAlgebra, optimize as opt

args = parse_args(__doc__)
dev = args.device
FULL = args.mode == 'full'
CONFIGS = ([(2, 1), (2, 2), (2, 16), (2, 32), (16, 1), (16, 2), (16, 16), (16, 32)] if FULL
           else [(2, 1), (2, 2), (4, 1)])
GRID = np.linspace(0.0, 0.99, 23 if FULL else 9)
GRID2 = np.linspace(0.0, 0.95, 11 if FULL else 5)        # joint (Delta(X), Delta(Y)) grid
CONFIGS_2D = [(2, 1), (2, 2), (4, 1)] if FULL else [(2, 1)]
STARTS_2D = 32 if FULL else 8
STARTS = 64 if FULL else 12
ROUNDS, STEPS = (10, 50) if FULL else (6, 40)
MC_DIR = ROOT / 'results' / 'experiments'


def trace_norm(op):
    return op.trace_norm()


def build(exp, alg, batch, delta_x, delta_y):
    """Parametrised operators and the objective z for one experiment."""
    if exp == 3:
        X = opt.SelfAdjointParam(alg, batch, delta=delta_x)
    else:
        X = opt.PositiveParam(alg, batch, delta=delta_x)
    Y = opt.PositiveParam(alg, batch, delta=delta_y)
    U = opt.UnitaryParam(alg, batch)
    if exp == 1:
        def z():
            x, y = X.operator(), Y.operator()
            return (x @ U.operator() @ y).trace.abs() - (x @ y).trace.real
        params = [X, Y, U]
    elif exp == 2:
        def z():
            x0, y0 = X.operator(), Y.operator()
            return trace_norm(x0 @ U.operator() @ y0) - (x0 @ y0).trace.real
        params = [X, Y, U]
    else:
        def z():
            x, y = X.operator(), Y.operator()
            absx = X.abs_operator()                       # |X| = U diag(lambda) U^*, known exactly
            return (y @ absx @ y).trace.real - trace_norm(y @ x @ y)
        params = [X, Y]
    return z, params


def envelope(exp, k, C, along):
    """sup and inf of z on the grid of Delta(X) (along='x') or Delta(Y) (along='y')."""
    alg = TypeIAlgebra([k] * C, [k] * C, complex_valued=False, precision='double', device=dev)
    G = len(GRID)
    batch = G * STARTS
    grid = torch.tensor(GRID, dtype=torch.float64, device=dev).repeat_interleave(STARTS)
    out = {}
    for name, maximize in (('sup', True), ('inf', False)):
        torch.manual_seed(args.seed + 17 * exp + (0 if maximize else 1))
        dx, dy = (grid, None) if along == 'x' else (None, grid)
        z, params = build(exp, alg, batch, dx, dy)
        vals = opt.extremize(z, params, maximize=maximize, rounds=ROUNDS, steps=STEPS)
        vals = vals.reshape(G, STARTS)
        vals = torch.nan_to_num(vals, nan=-math.inf if maximize else math.inf)
        best = vals.max(dim=1)[0] if maximize else vals.min(dim=1)[0]
        out[name] = best.cpu().numpy()
    return out


def envelope_2d(exp, k, C):
    """sup and inf of z with both contrasts prescribed: Delta(X) = a, Delta(Y) = b on GRID2 x GRID2."""
    alg = TypeIAlgebra([k] * C, [k] * C, complex_valued=False, precision='double', device=dev)
    G = len(GRID2)
    g = torch.tensor(GRID2, dtype=torch.float64, device=dev)
    dx = g.repeat_interleave(G * STARTS_2D)
    dy = g.repeat_interleave(STARTS_2D).repeat(G)
    out = {}
    for name, maximize in (('sup', True), ('inf', False)):
        torch.manual_seed(args.seed + 31 * exp + (0 if maximize else 1))
        z, params = build(exp, alg, len(dx), dx, dy)
        vals = opt.extremize(z, params, maximize=maximize, rounds=ROUNDS, steps=STEPS).reshape(G, G, STARTS_2D)
        vals = torch.nan_to_num(vals, nan=-math.inf if maximize else math.inf)
        out[name] = (vals.max(dim=-1)[0] if maximize else vals.min(dim=-1)[0]).cpu().numpy()
    return out


def load_mc(exp, k, C):
    path = MC_DIR / f'exp{exp}' / 'tables' / f'exp{exp}_dim{k}_ch{C}.csv.gz'
    if not path.exists():
        return None
    import pandas as pd
    return pd.read_csv(path)


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LABELS = {1: ('$\\Delta(X)$', '$\\Delta(Y)$'), 2: ('$\\Delta(|X|)$', '$\\Delta(|Y|)$'),
          3: ('$\\Delta(|X|)$', '$\\Delta(Y)$')}
results, rows = [], []
t_start = time.time()
for exp in (1, 2, 3):
    for k, C in CONFIGS:
        t0 = time.time()
        env = {along: envelope(exp, k, C, along) for along in ('x', 'y')}
        sync(dev)
        mc = load_mc(exp, k, C)
        outside = None
        if mc is not None:
            tol = 1e-6 + 1e-3 * np.abs(mc['z']).max()
            mc = mc[mc['deltaX'] <= GRID[-1]]                    # only inside the grid
            ub = np.interp(mc['deltaX'], GRID, env['x']['sup'])   # linear interpolation between
            lb = np.interp(mc['deltaX'], GRID, env['x']['inf'])   # grid points (can cut corners)
            bad = mc[(mc['z'] > ub + tol) | (mc['z'] < lb - tol)]
            interp_outside = len(bad)
            # linear interpolation between grid points can cut corners of the envelope: re-optimise at
            # the exact contrast of every such sample and count those that are still outside
            outside = 0
            for _, row in bad.head(50).iterrows():
                above = row['z'] > np.interp(row['deltaX'], GRID, env['x']['sup'])
                alg1 = TypeIAlgebra([k] * C, [k] * C, complex_valued=False, precision='double', device=dev)
                z1, p1 = build(exp, alg1, STARTS, torch.full((STARTS,), float(row['deltaX']), dtype=torch.float64,
                                                             device=dev), None)
                v = opt.extremize(z1, p1, maximize=above, rounds=ROUNDS, steps=STEPS)
                bound = v.max().item() if above else v.min().item()
                outside += (row['z'] > bound + tol) if above else (row['z'] < bound - tol)
            outside += max(0, len(bad) - 50)
        r = dict(experiment=exp, k=k, C=C, grid=GRID.tolist(),
                 sup_vs_x=env['x']['sup'].tolist(), inf_vs_x=env['x']['inf'].tolist(),
                 sup_vs_y=env['y']['sup'].tolist(), inf_vs_y=env['y']['inf'].tolist(),
                 mc_samples=None if mc is None else len(mc), mc_outside=outside,
                 mc_outside_interpolated=None if mc is None else interp_outside,
                 mc_max=None if mc is None else float(mc['z'].max()),
                 mc_min=None if mc is None else float(mc['z'].min()),
                 seconds=time.time() - t0)
        results.append(r)
        print(f"exp{exp} k={k:2d} C={C:2d}: sup z = {max(r['sup_vs_x']):+.4g}, inf z = {min(r['inf_vs_x']):+.4g}"
              + ("" if mc is None else f"; Monte Carlo range [{r['mc_min']:+.4g}, {r['mc_max']:+.4g}], "
                 f"{interp_outside} of {len(mc)} samples outside the interpolated envelopes, "
                 f"{outside} after re-optimising at their exact contrast")
              + f"  ({r['seconds']:.0f} s)")

        fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
        for ax, along, col, label in ((axes[0], 'x', 'deltaX', LABELS[exp][0]),
                                      (axes[1], 'y', 'deltaY', LABELS[exp][1])):
            if mc is not None:
                ax.scatter(mc[col], mc['z'], s=1, alpha=0.25, color='#9ab', rasterized=True,
                           label=f'Monte Carlo ({len(mc)})')
            ax.plot(GRID, env[along]['sup'], 'r-', lw=1.5, label='sup (optimisation)')
            ax.plot(GRID, env[along]['inf'], 'b-', lw=1.5, label='inf (optimisation)')
            ax.axhline(0, color='k', lw=0.5)
            ax.set_xlabel(label)
        axes[0].set_ylabel('$z$')
        axes[0].legend(frameon=False, fontsize=7, markerscale=5)
        fig.suptitle(f'Experiment {exp}: $k={k}$, $C={C}$', fontsize=10)
        fig.tight_layout()
        fig.savefig(f"{args.figdir}/bounds_exp{exp}_k{k}_C{C}.pdf")
        plt.close(fig)

# joint constraints: both contrasts prescribed
results_2d = []
for exp in (1, 2, 3):
    for k, C in CONFIGS_2D:
        t0 = time.time()
        e2 = envelope_2d(exp, k, C)
        results_2d.append(dict(experiment=exp, k=k, C=C, grid=GRID2.tolist(),
                               sup=e2['sup'].tolist(), inf=e2['inf'].tolist(), seconds=time.time() - t0))
        print(f"exp{exp} k={k} C={C}, Delta(X), Delta(Y) prescribed: sup on grid in "
              f"[{e2['sup'].min():+.3g}, {e2['sup'].max():+.3g}], inf in [{e2['inf'].min():+.3g}, "
              f"{e2['inf'].max():+.3g}]  ({time.time() - t0:.0f} s)")
        fig, axes = plt.subplots(1, 2, figsize=(8, 3.4))
        for ax, name in zip(axes, ('sup', 'inf')):
            im = ax.imshow(e2[name], origin='lower', extent=[GRID2[0], GRID2[-1], GRID2[0], GRID2[-1]],
                           aspect='auto', cmap='coolwarm' if name == 'inf' else 'viridis')
            ax.set_xlabel(LABELS[exp][1])
            ax.set_ylabel(LABELS[exp][0])
            ax.set_title(f'{name} z')
            fig.colorbar(im, ax=ax)
        fig.suptitle(f'Experiment {exp}, $k={k}$, $C={C}$: both contrasts prescribed', fontsize=10)
        fig.tight_layout()
        fig.savefig(f"{args.figdir}/bounds2d_exp{exp}_k{k}_C{C}.pdf")
        plt.close(fig)

# analytic check: for Exp 1, inf over U and Y of z at Delta(X) = d equals -(kC - 2d/(1+d))
# (Tr(XY) <= Tr X <= kC - 1 + lambda_min, attained for Y = 1; |Tr(UYX)| can be made 0)
exact_dev = 0.0
for r in results:
    if r['experiment'] == 1:
        g = np.array(r['grid'])
        exact = -(r['k'] * r['C'] - 2 * g / (1 + g))
        r['inf_exact'] = exact.tolist()
        exact_dev = max(exact_dev, float(np.max(np.abs(np.array(r['inf_vs_x']) - exact) / np.abs(exact))))
print(f"Exp 1 lower envelope vs exact -(kC - 2d/(1+d)): max relative deviation {exact_dev:.1e}")

save_json(args, 'bounds', dict(results=results, results_2d=results_2d, starts=STARTS, steps=ROUNDS * STEPS,
                               seconds=time.time() - t_start))
tex = env_macro(args, 'Bounds')
tex += f"\\newcommand{{\\BoundsStarts}}{{{STARTS}}}\n\\newcommand{{\\BoundsSteps}}{{{ROUNDS * STEPS}}}\n"
tex += f"\\newcommand{{\\BoundsGrid}}{{{len(GRID)}}}\n"
tex += f"\\newcommand{{\\BoundsExactDev}}{{{exact_dev:.1e}}}\n"
tex += "\\newcommand{\\BoundsRows}{%\n"
for r in results:
    mc = ('--' if r['mc_samples'] is None else
          f"$[{r['mc_min']:+.3g},\\,{r['mc_max']:+.3g}]$ & {r['mc_outside']}")
    if r['mc_samples'] is None:
        mc = '-- & --'
    tex += (f"{r['experiment']} & {r['k']} & {r['C']} & ${min(r['inf_vs_x']):+.4g}$ & ${max(r['sup_vs_x']):+.4g}$ "
            f"& {mc} \\\\\n")
tex += "}\n"
write_tex(args, 'bounds', tex)
print("written:", args.out, args.figdir, f"total {time.time() - t_start:.0f} s")
