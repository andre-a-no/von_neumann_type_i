#!/usr/bin/env python3
"""
Monte Carlo study for man/trace_envelopes (plain NumPy, independent of torch_vn_algebra).

For each of the 25 trace-characterising inequalities L <= R of scripts/trace_characterizations.py, random functionals
phi_W = Tr(W .) and random operators are drawn and the normalised violation (L - R)/N is recorded together with the
Michelson contrast of the density, Delta(W), and of the operators, Delta(.):

  density   W = diag(w), lambda_max = 1, lambda_min = (1 - d)/(1 + d), the others uniform in between, d ~ U[0, 1);
  positive  U diag(lam) U^*, U Haar, lam_max = 1, lam_min = (1 - d')/(1 + d'), d' ~ U[0, 1], times a log-uniform scale;
  arbitrary V |X| with V Haar and |X| positive as above;  self-adjoint: U diag(+-lam) U^*;
  projection U diag(1, .., 1, 0, .., 0) U^* of rank 1 .. n-1.

A second sample replaces every positive variable (and every modulus) by a scaled projection of random rank
('projection sample'). For each item the script reports the fraction of samples that violate the inequality, the largest
sampled violation relative to the envelope.

For every dimension the samples are summarised in bins of Delta(W) and of the contrast of the operators (largest
violation, 99th percentile, fraction violating). The largest ratio of a sample to the computed envelope uses the
envelope interpolated linearly in Delta(W) at the sample's own contrast. The figures (one per item, function figures()
below) show the samples with the computed envelopes, i.e. the certified (mpmath) values of
man/article/generated/trace_characterizations.json (top rows) and man/trace_envelopes/generated/trace_contrast.json
(bottom rows); run scripts/trace_characterizations.py and scripts/trace_contrast.py first.

Writes man/trace_envelopes/generated/trace_montecarlo.{json,tex}, generated/trace_montecarlo_figures.tex and
figures/mc_*.pdf.

    python scripts/trace_montecarlo.py [--samples 200000] [--dims 2,3,4] [--seed 0]
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
OUT = ROOT / 'man' / 'trace_envelopes'

# ---------------------------------------------------------------------------------------------------- linear algebra
def herm(A):
    return (A + np.conj(np.swapaxes(A, -1, -2))) / 2


def adj(A):
    return np.conj(np.swapaxes(A, -1, -2))


def fn(A, f):
    """f(A) for a batch of (nearly) self-adjoint matrices."""
    e, V = np.linalg.eigh(herm(A))
    return (V * f(e)[..., None, :]) @ adj(V)


def absm(A):
    """|A| = (A^* A)^{1/2}."""
    return fn(adj(A) @ A, lambda e: np.sqrt(np.clip(e, 0, None)))


def sqrtp(A):
    return fn(A, lambda e: np.sqrt(np.clip(e, 0, None)))


def proj_where(A, pred):
    """Spectral projection of a self-adjoint A onto the eigenvalues e with pred(e)."""
    e, V = np.linalg.eigh(herm(A))
    return (V * pred(e).astype(float)[..., None, :]) @ adj(V)


def phi(w, X):
    return np.real(np.einsum('si,sii->s', w, X))


# ---------------------------------------------------------------------------------------------------- sampling
class Sampler:
    def __init__(self, rng, S, n, proj_mode, dW=None, dop=None):
        self.rng, self.S, self.n, self.proj = rng, S, n, proj_mode
        self.dW, self.dop = dW, dop                 # optional: fixed Delta(W); one contrast per sample for every variable
        self.contrasts = []                         # contrasts of the operator variables of the current item

    def haar(self):
        Z = self.rng.standard_normal((self.S, self.n, self.n)) + 1j * self.rng.standard_normal((self.S, self.n, self.n))
        Q, R = np.linalg.qr(Z)
        d = np.diagonal(R, axis1=1, axis2=2)
        return Q * (d / np.abs(d))[:, None, :]

    def spectrum(self, d):
        l = (1 - d) / (1 + d)
        lam = l[:, None] + (1 - l[:, None]) * self.rng.uniform(0, 1, (self.S, self.n))
        lam[:, 0], lam[:, 1] = 1.0, l
        return self.rng.permuted(lam, axis=1)

    def density(self, total=1.0):
        d = self.rng.uniform(0, 1, self.S) * 0.999 if self.dW is None else np.full(self.S, self.dW)
        w = self.spectrum(d)
        return total * w / w.sum(1, keepdims=True), d

    def _projector_spectrum(self):
        r = self.rng.integers(1, self.n, self.S)
        return (np.arange(self.n)[None, :] < r[:, None]).astype(float)

    def pos(self, scale=True, invertible=False):
        if self.proj and not invertible:
            lam, d = self._projector_spectrum(), np.ones(self.S)
        else:
            d = self.rng.uniform(0, 1, self.S) if self.dop is None else self.dop
            d = np.minimum(d, 0.999) if invertible else d
            lam = self.spectrum(d)
        self.contrasts.append(d)
        U = self.haar()
        A = (U * lam[:, None, :]) @ adj(U)
        if scale:
            A = A * np.exp(self.rng.uniform(-1.5, 1.5, self.S))[:, None, None]
        return A

    def arb(self, left=True):
        P, V = self.pos(), self.haar()
        return (V @ P, P) if left else (P @ V, P)

    def sa(self, scale=1.0):
        if self.proj:
            lam, d = self._projector_spectrum(), np.ones(self.S)
        else:
            d = self.rng.uniform(0, 1, self.S) if self.dop is None else self.dop
            lam = self.spectrum(d)
        self.contrasts.append(d)
        lam = lam * np.where(self.rng.uniform(0, 1, (self.S, self.n)) < 0.5, -1.0, 1.0)
        U = self.haar()
        return (U * lam[:, None, :]) @ adj(U) * (scale * np.exp(self.rng.uniform(-1.0, 1.0, self.S)))[:, None, None]

    def projection(self):
        lam = self._projector_spectrum()
        U = self.haar()
        return (U * lam[:, None, :]) @ adj(U)


# ---------------------------------------------------------------------------------------------------- the 25 items
E_HALF = math.exp(0.5) - 1.0
I2 = None


def item_value(name, smp):
    """Returns ((L - R)/N per sample, Delta(W) per sample, operator contrast per sample or None)."""
    n, S = smp.n, smp.S
    I = np.broadcast_to(np.eye(n), (S, n, n))
    smp.contrasts = []
    total = float(n) if name in ('det', 'lp') else 1.0
    w, dW = smp.density(total)
    ph = lambda X: phi(w, X)                                       # noqa: E731
    if name == 'gardner':
        X, P = smp.arb()
        v = (np.abs(np.einsum('si,sii->s', w, X)) - ph(P)) / ph(P)
    elif name == 'triangle':
        X, Y = smp.sa(), smp.sa()
        r = ph(absm(X)) + ph(absm(Y))
        v = (ph(absm(X + Y)) - r) / r
    elif name == 'realpart':
        A, P = smp.arb(left=False)
        v = (ph(A + adj(A)) - 2 * ph(P)) / (2 * ph(P))
    elif name == 'square':
        A, P = smp.arb()
        v = np.abs(np.einsum('si,sii->s', w, A @ A)) / ph(P @ P) - 1
    elif name == 'metric':
        P, Q, R = smp.projection(), smp.projection(), smp.projection()
        r = ph(absm(P - Q)) + ph(absm(Q - R))
        v = (ph(absm(P - R)) - r) / r
    elif name == 'sgm':
        A, B = smp.pos(invertible=True), smp.pos()
        ah, aih = sqrtp(A), fn(A, lambda e: np.clip(e, 1e-300, None) ** -0.5)
        G = aih @ sqrtp(ah @ B @ ah) @ aih
        Gh = sqrtp(G)
        r = np.sqrt(ph(A) * ph(B))
        v = ph(Gh @ A @ Gh) / r - 1
    elif name == 'alt':
        A, B = smp.pos(), smp.pos()
        bh = sqrtp(B)
        M = bh @ A @ bh
        r = ph(B @ A @ A @ B)
        v = (ph(M @ M) - r) / r
    elif name == 'young':
        (A, PA), (B, PB) = smp.arb(), smp.arb()
        r = 0.5 * (ph(PA @ PA) + ph(PB @ PB))
        v = (ph(absm(A @ adj(B))) - r) / r
    elif name == 'projections':
        P, Q = smp.projection(), smp.projection()
        v = (2 * ph(P @ Q) - 2 * ph(Q @ P @ Q)) / (ph(P) + ph(Q))
    elif name == 'monotone':
        A, C = smp.pos(), smp.pos()
        B = A + C
        r = ph(B @ B)
        v = (ph(A @ A) - r) / r
    elif name == 'det':
        A = smp.pos(invertible=True)
        g = np.exp(np.mean(np.log(np.clip(np.linalg.eigvalsh(herm(A)), 1e-300, None)), axis=1))
        v = 1 - ph(A) / (n * g)
    elif name == 'subadd':
        A, B = smp.pos(), smp.pos()
        r = ph(sqrtp(A)) + ph(sqrtp(B))
        v = (ph(sqrtp(A + B)) - r) / r
    elif name == 'powers':
        A, B = smp.pos(), smp.pos()
        tot = ph(A + B)
        a4 = fn(A, lambda e: np.clip(e, 0, None) ** 0.25)
        v = (tot - ph(absm(A - B)) - 2 * ph(a4 @ sqrtp(B) @ a4)) / tot
    elif name == 'holder':
        P, Q = smp.projection(), smp.projection()
        r = ph(P) ** (1 / 3) * ph(Q) ** (2 / 3)
        v = (np.abs(np.einsum('si,sii->s', w, P @ Q)) - r) / r
    elif name == 'gt':
        P, Q = smp.projection(), smp.projection()
        eh, eq = I + E_HALF * P, I + (math.e - 1) * Q
        r = ph(eh @ eq @ eh)
        v = (ph(fn(P + Q, np.exp)) - r) / r
    elif name == 'pb':
        A, B = smp.pos(scale=False), smp.pos(scale=False)            # positive H, K as in Bikchentaev's form
        s = np.exp(smp.rng.uniform(-1.0, 1.5, (2, S)))
        A, B = A * s[0][:, None, None], B * s[1][:, None, None]
        ea = ph(fn(A, np.exp))
        eh = fn(A, lambda e: np.exp(e / 2))
        r = ph(fn(A + B, np.exp))
        v = (ea * np.exp(ph(eh @ B @ eh) / ea) - r) / r
    elif name == 'pqp':
        P, Q = smp.projection(), smp.projection()
        v = ph(P @ Q @ P) / ph(Q) - 1
    elif name == 'absdiff':
        P, Q = smp.projection(), smp.projection()
        v = ph(absm(P - Q)) / ph(P + Q) - 1
    elif name == 'meet':
        P, Q = smp.projection(), smp.projection()
        M = proj_where(P + Q, lambda e: e > 2 - 1e-9)
        r = ph(I + M)
        v = (ph(P + Q) - r) / r
    elif name == 'expqpq':
        P, Q = smp.projection(), smp.projection()
        r = ph(I + (math.e - 1) * P)
        v = (ph(fn(Q @ P @ Q, np.exp)) - r) / r
    elif name == 'repos':
        P, Q = smp.projection(), smp.projection()
        v = -ph(herm(P @ Q)) / (ph(P) + ph(Q))
    elif name == 'reabs':
        P, Q = smp.projection(), smp.projection()
        r = ph(absm(P @ Q))
        v = ph(herm(P @ Q)) / r - 1
    elif name == 'join':
        P, Q = smp.projection(), smp.projection()
        J = proj_where(P + Q, lambda e: e > 1e-9)
        r = ph(P) + ph(Q)
        v = (ph(J) - r) / r
    elif name == 'lp':
        A = smp.pos()
        v = np.sqrt(ph(A @ A)) / ph(A) - 1
    elif name == 'wmon':
        A, C = smp.pos(), smp.pos()
        B = A + C
        A2 = A @ A
        r = ph(A @ B @ B @ A)
        v = (ph(A2 @ A2) - r) / r
    else:
        raise KeyError(name)
    dop = np.min(np.stack(smp.contrasts), axis=0) if smp.contrasts else None
    return v, dW, dop


ITEMS = ['gardner', 'triangle', 'realpart', 'square', 'metric', 'sgm', 'alt', 'young', 'projections', 'monotone',
         'det', 'subadd', 'powers', 'holder', 'gt', 'pb', 'pqp', 'absdiff', 'meet', 'expqpq', 'repos', 'reabs', 'join',
         'lp', 'wmon']
OPERATOR = ['gardner', 'triangle', 'realpart', 'square', 'sgm', 'alt', 'young', 'monotone', 'det', 'subadd', 'powers',
            'pb', 'lp', 'wmon']
SHORT = dict(gardner='Gardner', triangle='triangle', realpart=r'$\varphi(A{+}A^*)$', square=r'$|\varphi(A^2)|$',
             metric='metric', sgm='geom. mean', alt='ALT', young='Young (ext.)', projections='Bikchentaev-Abed',
             monotone='monotonicity', det='AM-GM', subadd='subadditivity', powers='Powers-Stormer', holder='Hoelder',
             gt='Golden-Thompson', pb='Peierls-Bogoliubov', pqp=r'$\varphi(PQP)$', absdiff=r'$\varphi(|P-Q|)$',
             meet=r'$P\wedge Q$', expqpq=r'$e^{QPQ}$', repos=r'Re$PQ\geq0$', reabs=r'Re$PQ\leq|PQ|$',
             join=r'$P\vee Q$', lp='power mean', wmon='weighted mon.')


def envelope_curve(name, n):
    """(grid, certified search values, proved value or bound or None) from the main computation."""
    import trace_characterizations as tc
    d = json.loads((ROOT / 'man' / 'article' / 'generated' / 'trace_characterizations.json').read_text())
    grid = np.array(d['grid'], dtype=float)
    V = np.array(d['inequalities'][name]['dims'][str(n)]['certified'], dtype=float)
    xs = np.linspace(0, 0.95, 120)
    try:
        ref = tc.reference_curve(name, xs, n)
    except Exception:                                       # noqa: BLE001
        ref = None
    return grid, V, xs, ref


DELTA_BINS = np.concatenate([[0.0], np.logspace(-2, np.log10(0.95), 10), [1.0]])
OP_BINS = np.linspace(0, 1, 11)


def binned(x, v, edges):
    """Per bin of x: the largest v, the 99th percentile of v and the fraction of violating samples (v > 0)."""
    idx = np.digitize(x, edges) - 1
    out = dict(max=[], q99=[], fraction=[], count=[])
    for b in range(len(edges) - 1):
        m = idx == b
        if m.sum() < 20:
            for k in ('max', 'q99', 'fraction'):
                out[k].append(None)
            out['count'].append(int(m.sum()))
            continue
        out['max'].append(float(v[m].max()))
        out['q99'].append(float(np.quantile(v[m], 0.99)))
        out['fraction'].append(float(np.mean(v[m] > 1e-12)))
        out['count'].append(int(m.sum()))
    return out


def sample_all(n, samples, seed):
    res = {}
    for k, name in enumerate(ITEMS):
        out = {}
        for mode in ('random', 'projection'):
            if mode == 'projection' and (name not in OPERATOR or name in ('det', 'sgm')):
                continue                                  # projections by definition, or an invertible operator needed
            rng = np.random.default_rng(seed + 1000 * k + 100 * n + (7 if mode == 'projection' else 0))
            smp = Sampler(rng, samples, n, mode == 'projection')
            with np.errstate(all='ignore'):
                v, dW, dop = item_value(name, smp)
            ok = np.isfinite(v)
            out[mode] = dict(v=v[ok], dW=dW[ok], dop=None if dop is None else dop[ok])
        res[name] = out
    return res



# ---------------------------------------------------------------------------------------------------- figures
FIXED_DW = 0.9                                      # Delta(W) of the second sample (contrast of the operators varies)
GREY, ORANGE, INK, MUTED, GRID = '#7d848c', '#eb6834', '#0b0b0b', '#52514e', '#e6e5e1'
YLABEL = r'$(\mathcal{L}-\mathcal{R})/\mathcal{N}$'           # the paper's normalised violation


def sample_fixed(n, samples, seed):
    """Operator items: Delta(W) = FIXED_DW and every operator variable of one contrast delta' ~ U[0, 1)."""
    res = {}
    for k, name in enumerate(OPERATOR):
        rng = np.random.default_rng(seed + 50000 + 1000 * k + 100 * n)
        smp = Sampler(rng, samples, n, False, dW=FIXED_DW, dop=rng.uniform(0, 1, samples))
        with np.errstate(all='ignore'):
            v, _, dop = item_value(name, smp)
        ok = np.isfinite(v)
        res[name] = dict(v=v[ok], dop=dop[ok])
    return res


def _style(ax):
    ax.spines[['top', 'right']].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=7)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.axhline(0, color=MUTED, lw=0.7)


def _limits(curve, samples):
    top = max(np.max(curve), 0.0)
    if top <= 1e-12:                                 # no violation: scale to the samples
        top = max(float(np.max(samples)) if len(samples) else 0.0, 1e-6)
    return -0.2 * top, 1.12 * top


def figures(res, fixed, dims, fixed_dims):
    """One figure per item: Monte Carlo samples and the computed envelope (certified values, a lower bound for the
    supremum), one panel per dimension."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import trace_characterizations as tc
    num = {nm: i for i, nm in enumerate(tc.INEQS, 1)}
    env = json.loads((ROOT / 'man' / 'article' / 'generated' / 'trace_characterizations.json').read_text())
    grid = np.array(env['grid'], dtype=float)
    con = json.loads((OUT / 'generated' / 'trace_contrast.json').read_text())
    vg, gb = con['config']['vgrid'], con['config']['grid_b']
    two = {}
    for r in con['B']:
        if r['kind'] == 'free':
            two.setdefault((r['name'], r['n']), {})[r['vdelta']] = r['cert'][gb.index(FIXED_DW)]
    figdir = OUT / 'figures'
    figdir.mkdir(exist_ok=True)
    for old in figdir.glob('mc_*.pdf'):
        old.unlink()
    tex = []
    for name in ITEMS:
        op = name in OPERATOR
        fig, axes = plt.subplots(2 if op else 1, len(dims), figsize=(10, 5.1 if op else 2.55), squeeze=False)
        for j, n in enumerate(dims):
            ax = axes[0, j]
            V = np.array(env['inequalities'][name]['dims'][str(n)]['certified'], dtype=float)
            allv = []
            for mode, c, a in (('projection', ORANGE, 0.18), ('random', GREY, 0.12)):   # random on top
                if mode not in res[n][name]:
                    continue
                r = res[n][name][mode]
                m = r['dW'] <= grid[-1]
                ax.scatter(r['dW'][m], r['v'][m], s=0.8, color=c, alpha=a, lw=0, rasterized=True)
                allv.append(r['v'][m])
            ax.plot(grid, V, '-o', color=INK, lw=1.4, ms=3.5, mfc='white', mew=1.0, zorder=3)
            ax.set_xlim(0, grid[-1] + 0.01)
            ax.set_ylim(*_limits(V, np.concatenate(allv)))
            ax.set_title(f'$M_{n}$', fontsize=9, color=INK)
            ax.set_xlabel(r'$\Delta(W)$', fontsize=8, color=MUTED)
            _style(ax)
        axes[0, 0].set_ylabel(YLABEL, fontsize=8, color=MUTED)
        handles = [plt.Line2D([], [], ls='', marker='o', ms=4, color=GREY, label='random operators'),
                   plt.Line2D([], [], ls='', marker='o', ms=4, color=ORANGE, label='scaled projections'),
                   plt.Line2D([], [], color=INK, marker='o', ms=3.5, mfc='white', label='computed envelope (lower bound)')]
        if not any('projection' in res[n][name] for n in dims):
            handles.pop(1)
        if op:
            for j, n in enumerate(fixed_dims):
                ax = axes[1, j]
                r = fixed[n][name]
                ax.scatter(r['dop'], r['v'], s=0.8, color=GREY, alpha=0.12, lw=0, rasterized=True)
                y = np.array([two[(name, n)][x] for x in vg], dtype=float)
                ax.plot(vg, y, '-o', color=INK, lw=1.4, ms=3.5, mfc='white', mew=1.0, zorder=3)
                ax.set_xlim(0, 1.01)
                ax.set_ylim(*_limits(y, r['v']))
                ax.set_title(f'$M_{n}$, $\\Delta(W)={FIXED_DW}$', fontsize=9, color=INK)
                ax.set_xlabel(r"contrast $\delta'$ of the operators", fontsize=8, color=MUTED)
                _style(ax)
            axes[1, 0].set_ylabel(YLABEL, fontsize=8, color=MUTED)
            for ax in axes[1, len(fixed_dims):]:
                ax.axis('off')
            axes[1, -1].legend(handles=handles, loc='center', frameon=False, fontsize=8)
            fig.tight_layout()
        else:
            fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, fontsize=8)
            fig.tight_layout(rect=(0, 0.1, 1, 1))
        fname = f'mc_{num[name]:02d}_{name}.pdf'
        fig.savefig(figdir / fname, metadata={'CreationDate': None}, dpi=200)
        plt.close(fig)
        cap = (f'Item {num[name]}, {tc.INEQS[name]["label"]}: random samples and the computed envelope (a certified lower bound)'
               + (r'; top: against $\Delta(W)$, bottom: against the contrast $\delta^\prime$ of the operators at '
                  f'$\\Delta(W)={FIXED_DW}$.' if op else r' against $\Delta(W)$.'))
        tex += [r'\begin{figure}[!htbp]', r'\centering', f'\\includegraphics[width=\\textwidth]{{{fname}}}',
                f'\\caption[]{{{cap}}}', f'\\label{{fig:mc-{name}}}', r'\end{figure}', '']
    (OUT / 'generated' / 'trace_montecarlo_figures.tex').write_text(
        '% generated by scripts/trace_montecarlo.py -- do not edit by hand\n' + '\n'.join(tex))

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--samples', type=int, default=200000)
    ap.add_argument('--dims', default='2,3,4')
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()
    dims = [int(x) for x in args.dims.split(',')]
    data = dict(samples=args.samples, dims=dims, seed=args.seed, delta_bins=DELTA_BINS.tolist(),
                op_bins=OP_BINS.tolist(), binned={}, summary={})
    allres, fixed = {}, {}
    fixed_dims = [n for n in dims if n in (2, 3)]          # the two-contrast envelopes are computed on M_2 and M_3
    for n in dims:
        res = allres[n] = sample_all(n, args.samples, args.seed)
        if n in fixed_dims:
            fixed[n] = sample_fixed(n, args.samples // 2, args.seed)
        data['binned'][str(n)] = {}
        data['summary'][str(n)] = {}
        for name in ITEMS:
            grid, V, _, _ = envelope_curve(name, n)
            data['binned'][str(n)][name] = {}
            data['summary'][str(n)][name] = {}
            for mode, r in res[name].items():
                env = np.interp(r['dW'], grid, V)
                inside = r['dW'] <= grid[-1]               # the envelope is computed up to its last grid value
                rel = np.where((env > 1e-6) & inside, r['v'] / np.maximum(env, 1e-300), np.nan)
                data['summary'][str(n)][name][mode] = dict(
                    fraction=float(np.mean(r['v'] > 1e-12)), max=float(r['v'].max()),
                    max_rel=float(np.nanmax(rel)) if np.isfinite(rel).any() else float('nan'), n=int(len(r['v'])))
                data['binned'][str(n)][name][mode] = dict(
                    delta=binned(r['dW'], r['v'], DELTA_BINS),
                    operator=None if r['dop'] is None else binned(r['dop'], r['v'], OP_BINS))
            s = data['summary'][str(n)][name]
            print(f"n={n} {name:12s} " + '; '.join(f"{m}: {100 * x['fraction']:.1f}% violating, max/V {x['max_rel']:.2f}"
                                                   for m, x in s.items()))
    gen = OUT / 'generated'
    gen.mkdir(parents=True, exist_ok=True)
    (gen / 'trace_montecarlo.json').write_text(json.dumps(data, indent=1))

    import trace_characterizations as tc
    num = {nm: i for i, nm in enumerate(tc.INEQS, 1)}

    def pct(x):
        return '--' if x is None or not np.isfinite(x) else (f'{100 * x:.1f}' if x >= 0.001 else
                                                              ('0' if x == 0 else '$<0.1$'))

    def rel(x):
        return '--' if x is None or not np.isfinite(x) else (f'{x:.2f}' if x > 0 else '$\\le0$')

    rows = []
    for name in ITEMS:
        cells = [str(num[name]), SHORT[name]]
        for n in dims:
            s = data['summary'][str(n)][name]
            p = s.get('projection')
            cells += [pct(s['random']['fraction']), pct(p['fraction']) if p else '--']
        rows.append(' & '.join(cells) + ' \\\\')
    worst = max(x['max_rel'] for n in dims for s in data['summary'][str(n)].values() for x in s.values()
                if np.isfinite(x['max_rel']))
    tex = ['% generated by scripts/trace_montecarlo.py -- do not edit by hand',
           f"\\newcommand{{\\TMCSamples}}{{{args.samples}}}",
           f"\\newcommand{{\\TMCDims}}{{{', '.join(map(str, dims))}}}",
           f"\\newcommand{{\\TMCWorstRel}}{{{worst:.3f}}}",
           '\\newcommand{\\TMCRows}{%', *rows, '}']
    (gen / 'trace_montecarlo.tex').write_text('\n'.join(tex) + '\n')
    figures(allres, fixed, dims, fixed_dims)
    print('wrote', gen / 'trace_montecarlo.json', 'and figures/mc_*.pdf')


if __name__ == '__main__':
    main()
