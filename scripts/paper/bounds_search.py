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

from common import ROOT, parse_args, save_json, write_tex, env_macro, sync, sci, load_partial, save_partial, finish_partial

from torch_vn_algebra import TypeIAlgebra, optimize as opt

args = parse_args(__doc__)
dev = args.device
FULL = args.mode in ('full', 'check')            # check: the sizes of full, minimal repetitions
CHECK = args.mode == 'check'
CONFIGS = ([(2, 1), (2, 2), (2, 16), (2, 32), (16, 1), (16, 2), (16, 16), (16, 32)] if FULL
           else [(2, 1), (2, 2), (4, 1)])
GRID = np.append(np.linspace(0.0, 0.99, 23 if FULL else 9), 0.999)   # last point: the samples near contrast 1
GRID2 = np.linspace(0.0, 0.95, 11 if FULL else 5)        # joint (Delta(X), Delta(Y)) grid
CONFIGS_2D = [(2, 1), (2, 2), (4, 1)] if FULL else [(2, 1)]
STARTS_2D = 32 if FULL else 8
STARTS = 64 if FULL else 12
ROUNDS, STEPS = (1, 2) if CHECK else (10, 50) if FULL else (6, 40)
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
        mc = load_mc(exp, k, C)
        mc_dropped = None
        if mc is not None:                                       # only inside the grid, along both axes
            inside = (mc['deltaX'] <= GRID[-1]) & (mc['deltaY'] <= GRID[-1])
            mc_dropped = int((~inside).sum())
            mc = mc[inside]
        r = load_partial(args, f'exp{exp}_k{k}_C{C}')
        if r is not None:                                        # finished in an interrupted earlier run
            env = {along: {n: np.array(r[f'{n}_vs_{along}']) for n in ('sup', 'inf')} for along in ('x', 'y')}
            outside, interp_outside = r['mc_outside'], r['mc_outside_interpolated']
        else:
            torch.manual_seed(args.seed + 100 * exp + 10 * k + C)    # per configuration: resumable runs agree
            env = {along: envelope(exp, k, C, along) for along in ('x', 'y')}
            sync(dev)
            outside = interp_outside = None
        if mc is not None and r is None:
            tol = 1e-6 + 1e-3 * np.abs(mc['z']).max()
            # every sample is checked against both envelopes: along Delta(X) and along Delta(Y)
            alg1 = TypeIAlgebra([k] * C, [k] * C, complex_valued=False, precision='double', device=dev)
            still_out, interp_out = set(), set()
            for along, col in (('x', 'deltaX'), ('y', 'deltaY')):
                ub = np.interp(mc[col], GRID, env[along]['sup'])   # linear interpolation between
                lb = np.interp(mc[col], GRID, env[along]['inf'])   # grid points (can cut corners)
                bad = mc[(mc['z'] > ub + tol) | (mc['z'] < lb - tol)]
                interp_out |= set(bad.index)
                # interpolation can cut corners of the envelope: re-optimise at the exact contrast of every
                # such sample (at most 50 per axis) and count those that are still outside
                for n_row, (idx, row) in enumerate(bad.iterrows()):
                    if n_row >= 50:
                        still_out.add(idx)
                        continue
                    above = row['z'] > np.interp(row[col], GRID, env[along]['sup'])
                    d = torch.full((STARTS,), float(row[col]), dtype=torch.float64, device=dev)
                    z1, p1 = build(exp, alg1, STARTS, *((d, None) if along == 'x' else (None, d)))
                    v = opt.extremize(z1, p1, maximize=above, rounds=ROUNDS, steps=STEPS)
                    v = torch.nan_to_num(v, nan=-math.inf if above else math.inf)  # a failed start must not hide a sample
                    bound = v.max().item() if above else v.min().item()
                    if (row['z'] > bound + tol) if above else (row['z'] < bound - tol):
                        still_out.add(idx)
            interp_outside, outside = len(interp_out), len(still_out)
        if r is None:
            r = dict(experiment=exp, k=k, C=C, grid=GRID.tolist(),
                     sup_vs_x=env['x']['sup'].tolist(), inf_vs_x=env['x']['inf'].tolist(),
                     sup_vs_y=env['y']['sup'].tolist(), inf_vs_y=env['y']['inf'].tolist(),
                     mc_samples=None if mc is None else len(mc), mc_dropped=mc_dropped, mc_outside=outside,
                     mc_outside_interpolated=interp_outside,
                     mc_max=None if mc is None else float(mc['z'].max()),
                     mc_min=None if mc is None else float(mc['z'].min()),
                     seconds=time.time() - t0)
            save_partial(args, f'exp{exp}_k{k}_C{C}', r)
        results.append(r)
        print(f"exp{exp} k={k:2d} C={C:2d}: sup z = {max(r['sup_vs_x'] + r['sup_vs_y']):+.4g}, "
              f"inf z = {min(r['inf_vs_x'] + r['inf_vs_y']):+.4g}"
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
        fig.savefig(f"{args.figdir}/bounds_exp{exp}_k{k}_C{C}.pdf", metadata={"CreationDate": None})
        plt.close(fig)

# joint constraints: both contrasts prescribed
h2 = 0.5 * (GRID2[1] - GRID2[0])                         # half a grid cell: imshow cells centred on the grid
results_2d = []
for exp in (1, 2, 3):
    for k, C in CONFIGS_2D:
        t0 = time.time()
        r2 = load_partial(args, f'2d_exp{exp}_k{k}_C{C}')
        if r2 is None:
            torch.manual_seed(args.seed + 1000 + 100 * exp + 10 * k + C)
            e2 = envelope_2d(exp, k, C)
            r2 = dict(experiment=exp, k=k, C=C, grid=GRID2.tolist(),
                      sup=e2['sup'].tolist(), inf=e2['inf'].tolist(), seconds=time.time() - t0)
            save_partial(args, f'2d_exp{exp}_k{k}_C{C}', r2)
        e2 = {n: np.array(r2[n]) for n in ('sup', 'inf')}
        results_2d.append(r2)
        print(f"exp{exp} k={k} C={C}, Delta(X), Delta(Y) prescribed: sup on grid in "
              f"[{e2['sup'].min():+.3g}, {e2['sup'].max():+.3g}], inf in [{e2['inf'].min():+.3g}, "
              f"{e2['inf'].max():+.3g}]  ({time.time() - t0:.0f} s)")
        fig, axes = plt.subplots(1, 2, figsize=(8, 3.4))
        for ax, name in zip(axes, ('sup', 'inf')):
            im = ax.imshow(e2[name], origin='lower', extent=[GRID2[0] - h2, GRID2[-1] + h2, GRID2[0] - h2, GRID2[-1] + h2],
                           aspect='auto', cmap='coolwarm' if name == 'inf' else 'viridis')
            ax.set_xlabel(LABELS[exp][1])
            ax.set_ylabel(LABELS[exp][0])
            ax.set_title(f'{name} z')
            fig.colorbar(im, ax=ax)
        fig.suptitle(f'Experiment {exp}, $k={k}$, $C={C}$: both contrasts prescribed', fontsize=10)
        fig.tight_layout()
        fig.savefig(f"{args.figdir}/bounds2d_exp{exp}_k{k}_C{C}.pdf", metadata={"CreationDate": None})
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
# analytic check: for Exp 2, |X| = X_0 and |Y| = Y_0, so z = ||X_0 V Y_0||_1 - Tr(X_0 Y_0); von Neumann's trace
# inequality, its lower counterpart and rearrangement give, sector by sector,
# sup z = sum_c sum_i x_ci^desc (y_ci^desc - y_ci^asc) = -inf z for given spectra. For k = 2, C = 1 this is
# (1 - lambda_min(X)) (1 - lambda_min(Y)); with Delta(|X|) = d and the other contrast free in [0, DMAX] the sup is
# 2d/(1+d) * (1 - lo(DMAX)), which tends to 2d/(1+d) as DMAX -> 1.
DMAX = 0.999                                                  # upper end of the free contrast (optimize default)
lo_max = (1 - DMAX) / (1 + DMAX)
exact_dev2 = None
for r in results:
    if r['experiment'] == 2 and r['k'] == 2 and r['C'] == 1:
        g = np.array(r['grid'])
        exact = 2 * g / (1 + g) * (1 - lo_max)
        mask = exact > 1e-3
        exact_dev2 = float(max(np.max(np.abs(np.array(r['sup_vs_x'])[mask] - exact[mask]) / exact[mask]),
                               np.max(np.abs(np.array(r['inf_vs_x'])[mask] + exact[mask]) / exact[mask])))
        print(f"Exp 2 envelopes (k=2, C=1) vs exact +-2d/(1+d)(1-lo): max relative deviation {exact_dev2:.1e}")
# conjectured closed forms for k = 2, C = 1 along Delta(X) = d (no proof; compared numerically):
# Exp 1: sup z = d^2 / (2 (1 + d)),  Exp 3: sup z = 1 - d
conj_dev = {}
for r in results:
    if r['k'] == 2 and r['C'] == 1 and r['experiment'] in (1, 3):
        g = np.array(r['grid'])
        conj = g ** 2 / (2 * (1 + g)) if r['experiment'] == 1 else 1 - g
        mask = conj > 1e-3
        conj_dev[r['experiment']] = float(np.max(np.abs(np.array(r['sup_vs_x'])[mask] - conj[mask]) / conj[mask]))
        print(f"Exp {r['experiment']} upper envelope (k=2, C=1) vs conjectured closed form: "
              f"max relative deviation {conj_dev[r['experiment']]:.1e}")

save_json(args, 'bounds', dict(results=results, results_2d=results_2d, starts=STARTS, steps=ROUNDS * STEPS,
                               seconds=time.time() - t_start))
tex = env_macro(args, 'Bounds')
tex += f"\\newcommand{{\\BoundsStarts}}{{{STARTS}}}\n\\newcommand{{\\BoundsSteps}}{{{ROUNDS * STEPS}}}\n"
tex += f"\\newcommand{{\\BoundsGrid}}{{{len(GRID)}}}\n"
tex += f"\\newcommand{{\\BoundsExactDev}}{{{exact_dev:.1e}}}\n"
tex += f"\\newcommand{{\\BoundsExactDevTwo}}{{{'--' if exact_dev2 is None else f'{exact_dev2:.1e}'}}}\n"
for e, name in ((1, 'One'), (3, 'Three')):
    tex += f"\\newcommand{{\\BoundsConjDev{name}}}{{{f'{conj_dev[e]:.1e}' if e in conj_dev else '--'}}}\n"
exp3_inf = max((abs(v) for r in results if r['experiment'] == 3 for v in r['inf_vs_x'] + r['inf_vs_y']), default=0.0)
mc_checked = sum(r['mc_samples'] or 0 for r in results)
mc_dropped = sum(r.get('mc_dropped') or 0 for r in results)
mc_out = sum(r['mc_outside'] or 0 for r in results)
tex += f"\\newcommand{{\\BoundsExpThreeInf}}{{{sci(exp3_inf, 1)}}}\n"          # max |inf z| of Exp. 3 (proven 0)
tex += f"\\newcommand{{\\BoundsMcChecked}}{{{mc_checked}}}\n\\newcommand{{\\BoundsMcDropped}}{{{mc_dropped}}}\n"
tex += f"\\newcommand{{\\BoundsMcOutside}}{{{mc_out}}}\n"
def num(x, d=4):
    return '0' if abs(x) < 1e-6 else f"{x:+.{d}g}"


tex += "\\newcommand{\\BoundsRows}{%\n"
for r in results:
    mc = ('--' if r['mc_samples'] is None else
          f"$[{num(r['mc_min'], 3)},\\,{num(r['mc_max'], 3)}]$ & {r['mc_outside']}")
    if r['mc_samples'] is None:
        mc = '-- & --'
    lo, hi = min(r['inf_vs_x'] + r['inf_vs_y']), max(r['sup_vs_x'] + r['sup_vs_y'])   # over both envelopes
    tex += (f"{r['experiment']} & {r['k']} & {r['C']} & ${num(lo)}$ & ${num(hi)}$ "
            f"& {mc} \\\\\n")
tex += "}\n"
write_tex(args, 'bounds', tex)
finish_partial(args)
print("written:", args.out, args.figdir, f"total {time.time() - t_start:.0f} s")
