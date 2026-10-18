#!/usr/bin/env python3
"""
How operations change the Michelson contrast (man/trace_envelopes, Section "Contrasts of composite operators").

Plain NumPy/SciPy, independent of torch_vn_algebra. For positive A, B on C^n with contrasts Delta(A), Delta(B) the
script samples random relative positions (Haar unitaries, the other eigenvalues uniform between the extreme ones) and
compares the contrast of

  sum         A + tB (t log-uniform)            0 <= Delta <= max(Delta(A), Delta(B))
  congruence  A^{1/2} B A^{1/2}                 tanh|eta_A - eta_B| <= Delta <= tanh(eta_A + eta_B),  eta = artanh Delta
  geom. mean  A # B                             tanh|eta_A - eta_B|/2 <= Delta <= tanh((eta_A + eta_B)/2)
  Jordan      AB + BA                           tanh|eta_A - eta_B| <= Delta <= sin(psi_A + psi_B),  psi = arcsin Delta;
                                                positive for every relative position iff Delta(A)^2 + Delta(B)^2 <= 1
  commutator  |phi_W(AB - BA)| / (m(A) m(B))    <= 2 Delta(W) Delta(A) Delta(B) (n even), m = (lambda_max + lambda_min)/2

with the proved bounds. The largest contrast of the Jordan product is computed exactly on M_2 (one angle) and by
optimisation (scipy, many starts) on M_3 and M_4. Writes man/trace_envelopes/figures/ops_*.pdf and
man/trace_envelopes/generated/contrast_operations.{json,tex}.

    python scripts/contrast_operations.py [--samples 100000] [--starts 24] [--seed 0]
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'man' / 'trace_envelopes'
DIMS = [2, 3, 4]
DB = 0.6                                         # contrast of B in the figures
GREY, ORANGE, INK, MUTED, GRID = '#7d848c', '#eb6834', '#0b0b0b', '#52514e', '#e6e5e1'


# ---------------------------------------------------------------------------------------------------- sampling
def adj(X):
    return np.conj(np.swapaxes(X, -1, -2))


def haar(rng, S, n):
    Z = rng.standard_normal((S, n, n)) + 1j * rng.standard_normal((S, n, n))
    try:
        Q, R = np.linalg.qr(Z)
    except np.linalg.LinAlgError:                    # numpy < 1.22: no stacked QR
        Q, R = (np.array(x) for x in zip(*(np.linalg.qr(z) for z in Z)))
    d = np.diagonal(R, axis1=1, axis2=2)
    return Q * (d / np.abs(d))[:, None, :]


def spectrum(rng, d, n):
    """lambda_max = 1, lambda_min = (1 - d)/(1 + d), the others uniform in between."""
    ell = (1 - d) / (1 + d)
    lam = ell[:, None] + (1 - ell[:, None]) * rng.uniform(0, 1, (len(d), n))
    lam[:, 0], lam[:, 1] = 1.0, ell
    return lam


def positive(rng, d, n, scale=False):
    U = haar(rng, len(d), n)
    lam = spectrum(rng, d, n)
    if scale:
        lam = lam * np.exp(rng.uniform(-1.5, 1.5, len(d)))[:, None]
    return (U * lam[:, None, :]) @ adj(U)


def fn(X, f):
    e, V = np.linalg.eigh((X + adj(X)) / 2)
    return (V * f(e)[..., None, :]) @ adj(V)


def extremes(X):
    e = np.linalg.eigvalsh((X + adj(X)) / 2)
    return e[..., 0], e[..., -1]


def contrast(X):
    lo, hi = extremes(X)
    return (hi - lo) / (hi + lo), lo


def mid(X):
    lo, hi = extremes(X)
    return (hi + lo) / 2


# ---------------------------------------------------------------------------------------------------- bounds
H = np.arctanh


def hyp(da, db, w=1.0):
    """(tanh(w|eta_A - eta_B|), tanh(w(eta_A + eta_B))): congruence (w = 1) and geometric mean (w = 1/2)."""
    return np.tanh(w * np.abs(H(da) - H(db))), np.tanh(w * (H(da) + H(db)))


def jordan_bounds(da, db):
    up = np.sin(np.minimum(np.arcsin(da) + np.arcsin(db), np.pi / 2))
    return np.tanh(np.abs(H(da) - H(db))), up


def jordan_sup2(s, t, grid=200001):
    """Largest contrast of AB + BA on M_2: A = diag(1, a), B with eigenvalues 1, b at the angle theta, c = cos^2 theta;
    det(AB + BA) = 4ab - p^2 c(1 - c), tr = 2(a + b + p c), p = (1 - a)(1 - b)."""
    a, b = (1 - s) / (1 + s), (1 - t) / (1 + t)
    p = (1 - a) * (1 - b)
    c = np.linspace(0, 1, grid)
    det, tr = 4 * a * b - p * p * c * (1 - c), a + b + p * c
    if np.any(det <= 0):
        return float('nan')
    return float(np.sqrt(np.clip(1 - det / tr ** 2, 0, None)).max())


def jordan_mu(s, t, grid=200001):
    """Smallest eigenvalue of AB + BA over 2 x 2 positive A, B with eigenvalues {1, a}, {1, b} (one angle)."""
    a, b = (1 - s) / (1 + s), (1 - t) / (1 + t)
    p = (1 - a) * (1 - b)
    c = np.linspace(0, 1, grid)
    T = a + b + p * c
    return float(np.min(T - np.sqrt(np.clip(T * T - 4 * a * b + p * p * c * (1 - c), 0, None))))


def jordan_sup3(s, t):
    """Largest contrast of AB + BA on M_n, n >= 3: A = 1 + A', B = 1 + B' (direct sums) with the 2 x 2 block of
    smallest eigenvalue mu; the eigenvalue 2 comes from the common eigenvector."""
    mu = jordan_mu(s, t)
    return (2 - mu) / (2 + mu)


def trace_test_factor(n, dw):
    """sup over states phi_W with Delta(W) = dw of sum_{i<=n/2} (w_i - w_{n+1-i}): dw for even n, k(kappa-1)/(k(kappa+1)+1)
    for n = 2k + 1."""
    if n % 2 == 0:
        return dw
    k, kap = n // 2, (1 + dw) / (1 - dw)
    return k * (kap - 1) / (k * (kap + 1) + 1)


# ---------------------------------------------------------------------------------------------------- Jordan optimisation
def jordan_sup_opt(n, s, t, starts, rng):
    """Largest contrast of AB + BA on M_n with Delta(A) = s, Delta(B) = t (s^2 + t^2 < 1, so AB + BA > 0): A diagonal,
    B = V diag V^*, V = exp(K) with K anti-Hermitian; the inner eigenvalues free in [ell, 1]."""
    a, b = (1 - s) / (1 + s), (1 - t) / (1 + t)
    m = n - 2

    def unpack(x):
        z = 1 / (1 + np.exp(-x[:2 * m]))
        la = np.concatenate([[1.0, a], a + (1 - a) * z[:m]])
        lb = np.concatenate([[1.0, b], b + (1 - b) * z[m:]])
        R = x[2 * m:].reshape(n, n)
        K = (R - R.T) / 2 + 1j * (R + R.T) / 2
        V = expm(K)
        return np.diag(la), (V * lb) @ V.conj().T

    def obj(x):
        A, B = unpack(x)
        e = np.linalg.eigvalsh(A @ B + B @ A)
        return -(e[-1] - e[0]) / (e[-1] + e[0])

    best = 0.0
    for _ in range(starts):
        x0 = rng.standard_normal(2 * m + n * n)
        r = minimize(obj, x0, method='L-BFGS-B', options=dict(maxiter=400))
        best = max(best, -r.fun)
    return best


# ---------------------------------------------------------------------------------------------------- figures
def style(ax):
    ax.spines[['top', 'right']].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=7)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def panel_row(name, title_y, xlabel, data, curves, legend, ylim=(-0.02, 1.02), xlim=(0, 1)):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(DIMS), figsize=(10, 3.3))
    for ax, n in zip(axes, DIMS):
        x, y = data[n]
        ax.scatter(x, y, s=0.8, color=GREY, alpha=0.12, lw=0, rasterized=True)
        for c in curves(n):
            ax.plot(*c['xy'], c.get('ls', '-'), color=c.get('color', INK), lw=1.4, marker=c.get('marker'),
                    ms=4, mfc='white', zorder=3)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_title(f'$M_{n}$', fontsize=9, color=INK)
        ax.set_xlabel(xlabel, fontsize=8, color=MUTED)
        style(ax)
    axes[0].set_ylabel(title_y, fontsize=8, color=MUTED)
    fig.legend(handles=legend(plt), loc='lower center', ncol=len(legend(plt)), frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(OUT / 'figures' / f'ops_{name}.pdf', metadata={'CreationDate': None}, dpi=200)
    plt.close(fig)


def positivity_figure(rng, S):
    """Random pairs (Delta(A), Delta(B)) in random relative position: is AB + BA positive? One panel per dimension."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(DIMS), figsize=(10, 3.6))
    t = np.linspace(0, np.pi / 2, 200)
    for ax, n in zip(axes, DIMS):
        sa, sb = rng.uniform(0, 1, S), rng.uniform(0, 1, S)
        P, Q = positive(rng, sa, n), positive(rng, sb, n)
        neg = np.linalg.eigvalsh(P @ Q + Q @ P)[:, 0] < 0
        ax.scatter(sa[~neg], sb[~neg], s=0.8, color=GREY, alpha=0.12, lw=0, rasterized=True)
        ax.scatter(sa[neg], sb[neg], s=1.0, color=ORANGE, alpha=0.35, lw=0, rasterized=True)
        ax.plot(np.cos(t), np.sin(t), color=INK, lw=1.6)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect('equal')
        ax.set_title(f'$M_{n}$', fontsize=9, color=INK)
        ax.set_xlabel(r'$\Delta(A)$', fontsize=8, color=MUTED)
        style(ax)
    axes[0].set_ylabel(r'$\Delta(B)$', fontsize=8, color=MUTED)
    handles = [plt.Line2D([], [], ls='', marker='o', ms=4, color=GREY, label=r'$AB+BA\geq0$'),
               plt.Line2D([], [], ls='', marker='o', ms=4, color=ORANGE, label=r'$AB+BA$ indefinite'),
               plt.Line2D([], [], color=INK, label=r'$\Delta(A)^2+\Delta(B)^2=1$')]
    fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(OUT / 'figures' / 'ops_jordan_positivity.pdf', metadata={'CreationDate': None}, dpi=200)
    plt.close(fig)


def sci(x):
    m, e = f'{x:.0e}'.split('e')
    return f'{m}\\times10^{{{int(e)}}}'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--samples', type=int, default=100000)
    ap.add_argument('--starts', type=int, default=24)
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    S = args.samples
    xs = np.linspace(0, 0.999, 400)
    res = dict(samples=S, db=DB, excess={}, jordan={}, positivity={})
    data = {k: {} for k in ('sum', 'congruence', 'geomean', 'jordan', 'commutator')}
    for n in DIMS:
        da = rng.uniform(0, 1, S)
        db = np.full(S, DB)
        A, B = positive(rng, da, n), positive(rng, db, n, scale=True)
        # sum
        c, _ = contrast(A + B)
        data['sum'][n] = (da, c)
        ex = {'sum': float(np.max(c - np.maximum(da, db)))}
        # congruence and geometric mean
        Ah, Aih = fn(A, np.sqrt), fn(A, lambda e: 1 / np.sqrt(e))
        c, _ = contrast(Ah @ B @ Ah)
        lo, hi = hyp(da, db)
        data['congruence'][n] = (da, c)
        ex['congruence'] = float(max(np.max(c - hi), np.max(lo - c)))
        G = Ah @ fn(Aih @ B @ Aih, np.sqrt) @ Ah
        c, _ = contrast(G)
        lo, hi = hyp(da, db, 0.5)
        data['geomean'][n] = (da, c)
        ex['geomean'] = float(max(np.max(c - hi), np.max(lo - c)))
        # Jordan product: contrast where positive; positivity over the whole square of contrasts
        J = A @ B + B @ A
        c, lam = contrast(J)
        ok = lam > 0
        lo, hi = jordan_bounds(da, db)
        data['jordan'][n] = (da[ok], c[ok])
        ex['jordan'] = float(max(np.max((c - hi)[ok]), np.max((lo - c)[ok])))
        sa, sb = rng.uniform(0, 1, S), rng.uniform(0, 1, S)
        P, Q = positive(rng, sa, n), positive(rng, sb, n)
        lmin = np.linalg.eigvalsh(P @ Q + Q @ P)[:, 0]
        inside = sa ** 2 + sb ** 2 <= 1
        res['positivity'][str(n)] = dict(indefinite_inside=int(np.sum((lmin < 0) & inside)),
                                         indefinite_outside=float(np.mean(lmin[~inside] < 0)))
        # commutator: trace test on states, Delta(A) = Delta(B) = DB, Delta(W) uniform
        dw = rng.uniform(0, 0.999, S)
        w = spectrum(rng, dw, n)
        w = w / w.sum(1, keepdims=True)
        A2, B2 = positive(rng, np.full(S, DB), n, scale=True), positive(rng, np.full(S, DB), n, scale=True)
        v = np.abs(np.einsum('si,sii->s', w, A2 @ B2 - B2 @ A2)) / (mid(A2) * mid(B2))
        data['commutator'][n] = (dw, v)
        bound = 2 * DB * DB * np.array([trace_test_factor(n, x) for x in dw])
        ex['commutator'] = float(np.max(v - bound))
        res['excess'][str(n)] = ex
        print(f'n={n}: largest excess over the bounds', {k: f'{x:.1e}' for k, x in ex.items()})

    # Jordan: largest contrast, M_2 exact, M_3 and M_4 by optimisation, on a grid of Delta(A) (Delta(B) = DB)
    jgrid = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.79]
    res['jordan']['grid'] = jgrid
    res['jordan']['2'] = [jordan_sup2(s, DB) for s in jgrid]
    for n in DIMS[1:]:
        res['jordan'][str(n)] = [jordan_sup_opt(n, s, DB, args.starts, rng) for s in jgrid]
        print(f'Jordan sup n={n}:', np.round(res['jordan'][str(n)], 4))
    print('Jordan sup n=2:', np.round(res['jordan']['2'], 4))
    res['jordan']['formula3'] = [jordan_sup3(s, DB) for s in jgrid]
    res['jordan']['dev3'] = float(max(abs(x - y) for n in DIMS[1:]
                                      for x, y in zip(res['jordan'][str(n)], res['jordan']['formula3'])))
    print('Jordan sup n>=3 formula:', np.round(res['jordan']['formula3'], 4), 'max deviation', res['jordan']['dev3'])
    res['jordan']['symmetric'] = {str(n): (jordan_sup2(DB, DB) if n == 2 else jordan_sup_opt(n, DB, DB, args.starts, rng))
                                  for n in DIMS}

    # figures
    (OUT / 'figures').mkdir(exist_ok=True)
    xa = r'$\Delta(A)$, with $\Delta(B)=' + f'{DB}' + '$'

    def leg(*items):
        return lambda plt: [plt.Line2D([], [], ls='', marker='o', ms=4, color=GREY, label='random positions')] + [
            plt.Line2D([], [], color=c, ls=ls, marker=m, ms=4, mfc='white', label=lab) for lab, c, ls, m in items]

    panel_row('sum', r'$\Delta(A+tB)$', xa, data['sum'],
              lambda n: [dict(xy=(xs, np.maximum(xs, DB)))],
              leg((r'$\max(\Delta(A),\Delta(B))$', INK, '-', None)))
    lo, hi = hyp(xs, DB)
    panel_row('congruence', r'$\Delta(A^{1/2}BA^{1/2})$', xa, data['congruence'],
              lambda n: [dict(xy=(xs, hi)), dict(xy=(xs, lo), ls='--')],
              leg((r'$\tanh(\eta_A+\eta_B)$', INK, '-', None), (r'$\tanh|\eta_A-\eta_B|$', INK, '--', None)))
    lo, hi = hyp(xs, DB, 0.5)
    panel_row('geomean', r'$\Delta(A\,\sharp\,B)$', xa, data['geomean'],
              lambda n: [dict(xy=(xs, hi)), dict(xy=(xs, lo), ls='--')],
              leg((r'$\tanh\frac{\eta_A+\eta_B}{2}$', INK, '-', None), (r'$\tanh\frac{|\eta_A-\eta_B|}{2}$', INK, '--', None)))
    xj = np.linspace(0, np.sqrt(1 - DB ** 2), 300)
    lo, hi = jordan_bounds(xj, DB)
    sup2 = np.array([jordan_sup2(s, DB, 20001) for s in xj])
    sup3 = np.array([jordan_sup3(s, DB) for s in xj])
    panel_row('jordan', r'$\Delta(AB+BA)$', xa, data['jordan'],
              lambda n: [dict(xy=(xj, hi), ls=':'), dict(xy=(xj, lo), ls='--'),
                         dict(xy=(xj, sup2 if n == 2 else sup3), color=ORANGE)]
              + ([dict(xy=(jgrid, res['jordan'][str(n)]), ls='', marker='o')] if n > 2 else []),
              leg((r'$\sin(\psi_A+\psi_B)$', INK, ':', None), (r'$\tanh|\eta_A-\eta_B|$', INK, '--', None),
                  ('largest value (Prop.)', ORANGE, '-', None), ('optimisation', INK, '', 'o')))
    panel_row('commutator', r'$|\varphi_W(AB-BA)|/(m(A)m(B))$', r'$\Delta(W)$, with $\Delta(A)=\Delta(B)=' + f'{DB}$',
              data['commutator'],
              lambda n: [dict(xy=(xs, 2 * DB * DB * np.array([trace_test_factor(n, x) for x in xs])))],
              leg(('sharp bound', INK, '-', None)), ylim=(-0.02, 0.8))

    positivity_figure(rng, S)

    gen = OUT / 'generated'
    gen.mkdir(exist_ok=True)
    (gen / 'contrast_operations.json').write_text(json.dumps(res, indent=1))
    worst = max(x for e in res['excess'].values() for x in e.values())
    sym = res['jordan']['symmetric']
    tex = ['% generated by scripts/contrast_operations.py -- do not edit by hand',
           f'\\newcommand{{\\COSamples}}{{{S}}}',
           f'\\newcommand{{\\COWorst}}{{{sci(worst)}}}',
           f"\\newcommand{{\\COIndefInside}}{{{sum(p['indefinite_inside'] for p in res['positivity'].values())}}}",
           f"\\newcommand{{\\COJordanTwo}}{{{sym['2']:.4f}}}",
           f"\\newcommand{{\\COJordanThree}}{{{sym['3']:.4f}}}",
           f"\\newcommand{{\\COJordanFour}}{{{sym['4']:.4f}}}",
           f"\\newcommand{{\\COJordanDev}}{{{sci(res['jordan']['dev3'])}}}",
           f"\\newcommand{{\\COJordanSin}}{{{np.sin(2 * np.arcsin(DB)):.4f}}}",
           f"\\newcommand{{\\COJordanTanh}}{{{np.tanh(2 * np.arctanh(DB)):.4f}}}"]
    (gen / 'contrast_operations.tex').write_text('\n'.join(tex) + '\n')
    print('wrote figures/ops_*.pdf and generated/contrast_operations.{json,tex}; symmetric Jordan sup:', sym)


if __name__ == '__main__':
    main()
