"""
Independent checks of the envelopes of scripts/bounds_search.py (no optimisation involved):

  exact_envelope       envelopes known in closed form
  direct_sum_bound     values attained by direct sums of copies of the k = 2, C = 1 solution
  structured_samples   NumPy samples with two-level spectra and structured bases, which reach the
                       neighbourhood of the extremes in every dimension (unlike Haar samples)

    python scripts/bounds_checks.py            # checks paper/generated/bounds.json, a few minutes on a CPU
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

DMAX = 0.999                                                  # upper end of the free contrast (optimize default)


def exact_envelope(exp, k, C, grid, along='x'):
    """Envelopes known in closed form (proofs in the paper), on the contrast grid: Exp. 1 inf along Delta(X),
    Exp. 2 sup and inf, Exp. 3 inf; None where nothing is known."""
    g = np.asarray(grid)
    if exp == 2:                          # sup = -inf = C floor(k/2) (1 - lo(delta)) (1 - lo(DMAX))
        v = C * (k // 2) * 2 * g / (1 + g) * (1 - (1 - DMAX) / (1 + DMAX))
        return {'sup': v, 'inf': -v}
    if exp == 1:                          # along either contrast (X and Y exchange roles in a real algebra)
        return {'inf': -(k * C - 2 * g / (1 + g))}
    if exp == 3:
        return {'inf': np.zeros_like(g)}
    return None


def structured_samples(exp, k, C, n, seed, dmax=DMAX):
    """
    Samples concentrated near the extremal configurations, computed with plain NumPy (independent of the
    library): every spectrum takes only the two values 1 and lambda_min (both present, so the contrast is
    exact), and the eigenbases and U are random signed permutations times direct sums of 2x2 rotations (one
    angle for all pairs, or independent angles, or no rotation). Unlike Haar samples these reach the
    neighbourhood of the extremes in every dimension, so they test the optimiser: none may lie outside the
    envelopes.
    """
    rng = np.random.default_rng(seed)

    def basis(n):
        th = rng.uniform(0, np.pi, (n, C, k // 2))
        common = rng.random(n) < 0.5
        th[common] = th[common][:, :1, :1]
        R = np.zeros((n, C, k, k))
        c, s, i = np.cos(th), np.sin(th), np.arange(k // 2)
        R[..., 2 * i, 2 * i] = c
        R[..., 2 * i + 1, 2 * i + 1] = c
        R[..., 2 * i, 2 * i + 1] = -s
        R[..., 2 * i + 1, 2 * i] = s
        if k % 2:
            R[..., -1, -1] = 1
        R[rng.random(n) < 0.25] = np.eye(k)
        P = np.zeros((n, C, k, k))
        np.put_along_axis(P, np.argsort(rng.random((n, C, k)), axis=-1)[..., None], 1.0, axis=-1)
        return (P * rng.choice([-1.0, 1.0], (n, C, k, 1))) @ R

    def two_level(delta):
        lam = (1 - delta) / (1 + delta)
        e = np.where(rng.random((n, C * k)) < rng.random((n, 1)), 1.0, lam[:, None])
        e[:, 0], e[:, 1] = 1.0, lam
        return np.take_along_axis(e, np.argsort(rng.random((n, C * k)), axis=1), axis=1).reshape(n, C, k)

    def op(e, V):
        return (V * e[..., None, :]) @ np.swapaxes(V, -1, -2)

    def tn(A):
        return np.linalg.svd(A, compute_uv=False).sum(axis=(-1, -2))

    def tr(A):
        return np.trace(A, axis1=-2, axis2=-1).sum(-1)

    dx, dy = rng.uniform(0, dmax, n), rng.uniform(0, dmax, n)
    dx[rng.random(n) < 0.2] = dmax
    dy[rng.random(n) < 0.2] = dmax
    ex, ey = two_level(dx), two_level(dy)
    Y = op(ey, basis(n))
    if exp == 3:
        V = basis(n)
        sg = np.where(rng.random((n, C, k)) < rng.random((n, 1, 1)), -1.0, 1.0)
        z = tr(Y @ op(ex, V) @ Y) - tn(Y @ op(ex * sg, V) @ Y)
    else:
        X, U = op(ex, basis(n)), basis(n)
        z = np.abs(tr(X @ U @ Y)) - tr(X @ Y) if exp == 1 else tn(X @ U @ Y) - tr(X @ Y)
    return pd.DataFrame(dict(deltaX=dx, deltaY=dy, z=z))


def direct_sum_bound(base, along, name, m):
    """
    Envelope of a direct sum of m copies of the k = 2, C = 1 problem (z is additive over orthogonal blocks; the
    contrast of the sum is the largest block contrast): with s the k = 2, C = 1 envelope along one contrast,
    sup z(delta) >= s(delta) + (m - 1) max_{delta' <= delta} s(delta')  (and the same with min for inf).
    A configuration with this value exists, so the computed envelope must not lie below it (above, for inf).
    """
    s = np.array(base[f'{name}_vs_{along}'])
    best = np.maximum.accumulate(s) if name == 'sup' else np.minimum.accumulate(s)
    return s + (m - 1) * best




def check(results, n_samples=4000, seed=0):
    """Compare every envelope in `results` (the list in bounds.json) with the three checks; returns rows."""
    rows = []
    by_key = {(r['experiment'], r['k'], r['C']): r for r in results}
    for r in results:
        exp, k, C, g = r['experiment'], r['k'], r['C'], np.array(r['grid'])
        env = {a: {n: np.array(r[f'{n}_vs_{a}']) for n in ('sup', 'inf')} for a in ('x', 'y')}
        scale = max(np.abs(c).max() for e in env.values() for c in e.values())
        row = dict(experiment=exp, k=k, C=C)
        ex = exact_envelope(exp, k, C, g)
        # shortfalls relative to the largest |value| of the reference curve
        if ex:
            row['exact_dev'] = max(float(np.max(sgn * (ex[n] - env['x'][n])) / max(np.abs(ex[n]).max(), 1e-3 * scale))
                                   for n, sgn in (('sup', 1), ('inf', -1)) if n in ex)
        base = by_key.get((exp, 2, 1))
        if base is not None and (k, C) != (2, 1) and k % 2 == 0:
            row['direct_sum_shortfall'] = max(
                float(np.max(sgn * (ref - env[a][n])) / max(np.abs(ref).max(), 1e-3 * scale))
                for a in ('x', 'y') for n, sgn in (('sup', 1), ('inf', -1))
                for ref in [direct_sum_bound(base, a, n, k * C // 2)])
        st = structured_samples(exp, k, C, n_samples if k * C <= 64 else n_samples // 4, seed + 7 * exp + k + C,
                                g[-1])
        tol = 1e-6 + 1e-3 * scale
        out = 0
        for a, col in (('x', 'deltaX'), ('y', 'deltaY')):
            out += int(((st['z'] > np.interp(st[col], g, env[a]['sup']) + tol)
                        | (st['z'] < np.interp(st[col], g, env[a]['inf']) - tol)).sum())
        row.update(structured=len(st), structured_outside=out, structured_max=float(st['z'].max()),
                   structured_min=float(st['z'].min()), envelope_max=float(max(env['x']['sup'].max(), env['y']['sup'].max())),
                   envelope_min=float(min(env['x']['inf'].min(), env['y']['inf'].min())))
        rows.append(row)
    return rows


if __name__ == '__main__':
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1] / 'paper' / 'generated' / 'bounds.json'
    data = json.load(open(path))
    print(f"{path} ({data.get('mode')} run on {data.get('environment', {}).get('gpu') or 'CPU'})")
    print("shortfall: how far an envelope stays inside a value that is attained, in % of the largest |value| of that\n"
          "curve (0 = fine); structured outside: samples beyond the envelopes (should be 0)")
    print(f"{'exp':>3s} {'k':>3s} {'C':>3s} {'exact':>8s} {'dir.sum':>8s} {'structured outside':>19s} "
          f"{'envelope range':>22s} {'structured range':>22s}")
    for w in check(data['results']):
        f = lambda v: '--' if v is None else f"{100 * max(v, 0):.2f}%"
        print(f"{w['experiment']:3d} {w['k']:3d} {w['C']:3d} {f(w.get('exact_dev')):>8s} "
              f"{f(w.get('direct_sum_shortfall')):>8s} {w['structured_outside']:>9d} of {w['structured']:<6d} "
              f"[{w['envelope_min']:+9.4g}, {w['envelope_max']:+9.4g}] [{w['structured_min']:+9.4g}, {w['structured_max']:+9.4g}]")
