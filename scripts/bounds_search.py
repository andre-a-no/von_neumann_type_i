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

    python scripts/bounds_search.py --mode full --device cuda
"""
import math
import time
from pathlib import Path

import numpy as np
import torch

from common import ROOT, parse_args, save_json, write_tex, env_macro, sync, sci, load_partial, save_partial, finish_partial

from torch_vn_algebra import TypeIAlgebra, optimize as opt
from bounds_checks import DMAX, structured_samples, direct_sum_bound, exact_envelope

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
STRUCTURED = 400 if CHECK else 4000 if FULL else 1000     # structured samples per configuration (k C <= 64)
REOPT = 200                                              # outside samples re-optimised per axis and direction
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
            M = x0 @ U.operator() @ y0
            M._is_self_adjoint = False                    # not self-adjoint: singular values, no check per step
            return trace_norm(M) - (x0 @ y0).trace.real
        params = [X, Y, U]
    else:
        def z():
            x, y = X.operator(), Y.operator()
            absx = X.abs_operator()                       # |X| = U diag(lambda) U^*, known exactly
            M = y @ x @ y
            M._is_self_adjoint = True                     # Y X Y is self-adjoint: eigenvalues, no check per step
            return (y @ absx @ y).trace.real - trace_norm(M)
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
        # independent starting points for every configuration, axis and direction (reproducible)
        torch.manual_seed(args.seed + 100000 * exp + 1000 * k + 10 * C + (0 if along == 'x' else 5)
                          + (0 if maximize else 1))
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
        # the envelope is the multistart search; the direct sums of the k = 2, C = 1 solution (computed first) are
        # admissible configurations and serve as an independent check (and are drawn in the figure)
        base = next((b for b in results if b['experiment'] == exp and b['k'] == 2 and b['C'] == 1), None)
        dsum = None if base is None or (k, C) == (2, 1) or k % 2 else \
            {along: {n: direct_sum_bound(base, along, n, k * C // 2) for n in ('sup', 'inf')} for along in ('x', 'y')}

        def outside_count(samples):
            """Samples outside the search envelopes (along either contrast), before and after re-optimising at
            their exact contrast (all outside samples of an axis and direction in one batched search, at most
            REOPT of them; further ones count as outside)."""
            tol = 1e-6 + 1e-3 * np.abs(samples['z']).max()
            alg1 = TypeIAlgebra([k] * C, [k] * C, complex_valued=False, precision='double', device=dev)
            still_out, interp_out = set(), set()
            for along, col in (('x', 'deltaX'), ('y', 'deltaY')):
                ub = np.interp(samples[col], GRID, env[along]['sup'])   # linear interpolation between
                lb = np.interp(samples[col], GRID, env[along]['inf'])   # grid points (can cut corners)
                for above, bad in ((True, samples[samples['z'] > ub + tol]), (False, samples[samples['z'] < lb - tol])):
                    interp_out |= set(bad.index)
                    if bad.empty:
                        continue
                    still_out |= set(bad.index[REOPT:])
                    bad = bad.iloc[:REOPT]
                    # batched: as many samples per search as fit the batch of one envelope (G * STARTS starts)
                    per = max(1, (len(GRID) * STARTS) // STARTS)
                    for i0 in range(0, len(bad), per):
                        part = bad.iloc[i0:i0 + per]
                        n = len(part)
                        d = torch.tensor(part[col].to_numpy(), dtype=torch.float64, device=dev).repeat_interleave(STARTS)
                        torch.manual_seed(args.seed + 7 + i0)
                        z1, p1 = build(exp, alg1, n * STARTS, *((d, None) if along == 'x' else (None, d)))
                        v = opt.extremize(z1, p1, maximize=above, rounds=ROUNDS, steps=STEPS).reshape(n, STARTS)
                        bound = (v.max(dim=1)[0] if above else v.min(dim=1)[0]).cpu().numpy()
                        zs = part['z'].to_numpy()
                        out = zs > bound + tol if above else zs < bound - tol
                        still_out |= set(part.index[out])
            return len(interp_out), len(still_out)

        st = structured_samples(exp, k, C, STRUCTURED if k * C <= 64 else STRUCTURED // 4, args.seed + 7 * exp + k + C,
                                GRID[-1])
        if r is None:
            if mc is not None:
                interp_outside, outside = outside_count(mc)
            st_interp, st_outside = outside_count(st)
        else:
            st_interp, st_outside = r['structured_outside_interpolated'], r['structured_outside']
        if r is None:
            r = dict(experiment=exp, k=k, C=C, grid=GRID.tolist(),
                     sup_vs_x=env['x']['sup'].tolist(), inf_vs_x=env['x']['inf'].tolist(),
                     sup_vs_y=env['y']['sup'].tolist(), inf_vs_y=env['y']['inf'].tolist(),
                     mc_samples=None if mc is None else len(mc), mc_dropped=mc_dropped, mc_outside=outside,
                     mc_outside_interpolated=interp_outside,
                     mc_max=None if mc is None else float(mc['z'].max()),
                     mc_min=None if mc is None else float(mc['z'].min()),
                     structured_samples=len(st), structured_outside=st_outside,
                     structured_outside_interpolated=st_interp,
                     structured_max=float(st['z'].max()), structured_min=float(st['z'].min()),
                     seconds=time.time() - t0)
            save_partial(args, f'exp{exp}_k{k}_C{C}', r)
        results.append(r)
        print(f"exp{exp} k={k:2d} C={C:2d}: sup z = {max(r['sup_vs_x'] + r['sup_vs_y']):+.4g}, "
              f"inf z = {min(r['inf_vs_x'] + r['inf_vs_y']):+.4g}"
              + ("" if mc is None else f"; Monte Carlo range [{r['mc_min']:+.4g}, {r['mc_max']:+.4g}], "
                 f"{interp_outside} of {len(mc)} samples outside the interpolated envelopes, "
                 f"{outside} after re-optimising at their exact contrast")
              + f"; structured samples [{r['structured_min']:+.4g}, {r['structured_max']:+.4g}], "
              f"{st_outside} of {len(st)} outside"
              + f"  ({r['seconds']:.0f} s)")
        if dsum is not None:
            # how far the computed envelope falls short of the direct-sum configurations (relative to the range of
            # z); positive: the optimiser missed a configuration that exists
            short = []
            scale = max(np.abs(env[a][n]).max() for n in ('sup', 'inf') for a in ('x', 'y'))
            for along in ('x', 'y'):
                for n, sgn in (('sup', 1), ('inf', -1)):
                    ref, got = dsum[along][n], env[along][n]
                    short.append(float(np.max(sgn * (ref - got)) / max(np.abs(ref).max(), 1e-3 * scale)))
            r['direct_sum_shortfall'] = max(short)
            print(f"    search vs direct sums of the k=2, C=1 solutions: largest shortfall "
                  f"{100 * r['direct_sum_shortfall']:.2f}% of the largest |value| of that curve")

        fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
        for ax, along, col, label in ((axes[0], 'x', 'deltaX', LABELS[exp][0]),
                                      (axes[1], 'y', 'deltaY', LABELS[exp][1])):
            if mc is not None:
                ax.scatter(mc[col], mc['z'], s=1, alpha=0.25, color='#9ab', rasterized=True,
                           label=f'Monte Carlo, Haar ({len(mc)})')
            ax.scatter(st[col], st['z'], s=1, alpha=0.35, color='#e8a33d', rasterized=True,
                       label=f'two-level spectra, structured bases ({len(st)})')
            ax.plot(GRID, env[along]['sup'], 'r-', lw=1.5, label='sup (optimisation)')
            ax.plot(GRID, env[along]['inf'], 'b-', lw=1.5, label='inf (optimisation)')
            for curve in (exact_envelope(exp, k, C, GRID) or {}).values():
                ax.plot(GRID, curve, 'k--', lw=0.9, label='exact')
            if dsum is not None and exp in (1, 3):
                ax.plot(GRID, dsum[along]['sup'], 'k:', lw=1.2, label='direct sum of $k=2$ solutions')
            ax.axhline(0, color='k', lw=0.5)
            ax.set_xlabel(label)
        axes[0].set_ylabel('$z$')
        handles, labels = axes[0].get_legend_handles_labels()
        seen = {}
        for h, lab in zip(handles, labels):
            seen.setdefault(lab, h)
        axes[0].legend(seen.values(), seen.keys(), frameon=False, fontsize=6, markerscale=5)
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
        # along Delta(Y) the same (Tr(X U Y) = Tr(Y U^T X) in a real algebra: X and Y exchange roles)
        for along in ('x', 'y'):
            exact_dev = max(exact_dev, float(np.max(np.abs(np.array(r[f'inf_vs_{along}']) - exact) / np.abs(exact))))
print(f"Exp 1 lower envelopes (both contrasts) vs exact -(kC - 2d/(1+d)): max relative deviation {exact_dev:.1e}")
# analytic check: for Exp 2, |X| = X_0 and |Y| = Y_0, so z = ||X_0 V Y_0||_1 - Tr(X_0 Y_0); von Neumann's trace
# inequality, its lower counterpart and rearrangement give, sector by sector,
# sup z = sum_c sum_i x_ci^desc (y_ci^desc - y_ci^asc) = -inf z for given spectra, and
# sum_i x_i^desc (y_i^desc - y_i^asc) = sum_{i <= k/2} (x_i^desc - x_{k+1-i}^desc)(y_i^desc - y_{k+1-i}^desc)
# <= floor(k/2) (1 - lo_X)(1 - lo_Y), attained for two-level spectra. With Delta(|X|) = d prescribed and the other
# contrast free in [0, DMAX]: sup z = -inf z = C floor(k/2) 2d/(1+d) (1 - lo(DMAX)).
exact_dev2 = 0.0
for r in results:
    if r['experiment'] == 2:
        for along in ('x', 'y'):
            ex = exact_envelope(2, r['k'], r['C'], GRID)['sup']
            mask = ex > 1e-3
            for n, sgn in (('sup', 1), ('inf', -1)):
                exact_dev2 = max(exact_dev2, float(np.max(np.abs(sgn * np.array(r[f'{n}_vs_{along}'])[mask] - ex[mask])
                                                          / ex[mask])))
print(f"Exp 2 envelopes (all k, C) vs exact +-C floor(k/2) 2d/(1+d)(1-lo): max relative deviation {exact_dev2:.1e}")
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
tex += f"\\newcommand{{\\BoundsExactDevTwo}}{{{exact_dev2:.1e}}}\n"
for e, name in ((1, 'One'), (3, 'Three')):
    tex += f"\\newcommand{{\\BoundsConjDev{name}}}{{{f'{conj_dev[e]:.1e}' if e in conj_dev else '--'}}}\n"
exp3_inf = max((abs(v) for r in results if r['experiment'] == 3 for v in r['inf_vs_x'] + r['inf_vs_y']), default=0.0)
mc_checked = sum(r['mc_samples'] or 0 for r in results)
mc_dropped = sum(r.get('mc_dropped') or 0 for r in results)
mc_out = sum(r['mc_outside'] or 0 for r in results)
tex += f"\\newcommand{{\\BoundsExpThreeInf}}{{{sci(exp3_inf, 1)}}}\n"          # max |inf z| of Exp. 3 (proven 0)
tex += f"\\newcommand{{\\BoundsMcChecked}}{{{mc_checked}}}\n\\newcommand{{\\BoundsMcDropped}}{{{mc_dropped}}}\n"
tex += f"\\newcommand{{\\BoundsMcOutside}}{{{mc_out}}}\n"
st_checked = sum(r['structured_samples'] for r in results)
st_out = sum(r['structured_outside'] for r in results)
shortfall = max((r.get('direct_sum_shortfall') or 0.0 for r in results), default=0.0)
tex += f"\\newcommand{{\\BoundsStChecked}}{{{st_checked}}}\n\\newcommand{{\\BoundsStOutside}}{{{st_out}}}\n"
tex += f"\\newcommand{{\\BoundsShortfall}}{{{100 * max(shortfall, 0.0):.1f}}}\n"         # percent of max|z|
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
