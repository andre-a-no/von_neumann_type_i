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
sampled violation relative to the envelope, and writes three figures:

  montecarlo_delta.pdf     (L - R)/N against Delta(W) for all 25 items, random and projection samples, with the envelope;
  montecarlo_operator.pdf  (L - R)/N against the contrast of the operators for the items with operator variables;
  montecarlo_rate.pdf      fraction of violating samples against Delta(W), random against projection samples.

Writes man/trace_envelopes/generated/trace_montecarlo.{json,tex} and the figures to man/trace_envelopes/figures.

    python scripts/trace_montecarlo.py [--samples 40000] [--n 2] [--seed 0]
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
    def __init__(self, rng, S, n, proj_mode):
        self.rng, self.S, self.n, self.proj = rng, S, n, proj_mode
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
        d = self.rng.uniform(0, 1, self.S) * 0.999
        w = self.spectrum(d)
        return total * w / w.sum(1, keepdims=True), d

    def _projector_spectrum(self):
        r = self.rng.integers(1, self.n, self.S)
        return (np.arange(self.n)[None, :] < r[:, None]).astype(float)

    def pos(self, scale=True, invertible=False):
        if self.proj and not invertible:
            lam, d = self._projector_spectrum(), np.ones(self.S)
        else:
            d = self.rng.uniform(0, 1, self.S) * (0.999 if invertible else 1.0)
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
            d = self.rng.uniform(0, 1, self.S)
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
    V = np.array(d['inequalities'][name]['dims'][str(n)]['V'], dtype=float)
    xs = np.linspace(0, 0.95, 120)
    try:
        ref = tc.reference_curve(name, xs, n)
    except Exception:                                       # noqa: BLE001
        ref = None
    return grid, V, xs, ref


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--samples', type=int, default=40000)
    ap.add_argument('--n', type=int, default=2)
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    res = {}
    for k, name in enumerate(ITEMS):
        out = {}
        for mode in ('random', 'projection'):
            if mode == 'projection' and name not in OPERATOR:
                continue
            if mode == 'projection' and name in ('det', 'sgm'):
                continue                                  # need an invertible operator
            rng = np.random.default_rng(args.seed + 1000 * k + (7 if mode == 'projection' else 0))
            smp = Sampler(rng, args.samples, args.n, mode == 'projection')
            with np.errstate(all='ignore'):
                v, dW, dop = item_value(name, smp)
            ok = np.isfinite(v)
            out[mode] = dict(v=v[ok], dW=dW[ok], dop=None if dop is None else dop[ok])
        res[name] = out
        print(f"{name:12s} random: violating {np.mean(out['random']['v'] > 1e-12):.3f}, max {out['random']['v'].max():+.4f}"
              + (f"; projections: violating {np.mean(out['projection']['v'] > 1e-12):.3f}, "
                 f"max {out['projection']['v'].max():+.4f}" if 'projection' in out else ''))

    # summary: fraction violating, largest sample relative to the certified envelope at its contrast
    summary = {}
    for name in ITEMS:
        grid, V, _, _ = envelope_curve(name, args.n)
        row = {}
        for mode, r in res[name].items():
            env = np.interp(r['dW'], grid, V)
            inside = r['dW'] <= grid[-1]                   # the envelope is computed up to its last grid value
            rel = np.where((env > 1e-6) & inside, r['v'] / np.maximum(env, 1e-300), np.nan)
            row[mode] = dict(fraction=float(np.mean(r['v'] > 1e-12)), max=float(r['v'].max()),
                             max_rel=float(np.nanmax(rel)) if np.isfinite(rel).any() else float('nan'),
                             n=int(len(r['v'])))
        summary[name] = row
    gen = OUT / 'generated'
    gen.mkdir(parents=True, exist_ok=True)
    (gen / 'trace_montecarlo.json').write_text(json.dumps(dict(samples=args.samples, n=args.n, seed=args.seed,
                                                               summary=summary), indent=1))
    import trace_characterizations as tc
    num = {nm: i for i, nm in enumerate(tc.INEQS, 1)}

    def pct(x):
        return '--' if x is None or not np.isfinite(x) else (f'{100 * x:.1f}' if x >= 0.001 else
                                                              ('0' if x == 0 else '$<0.1$'))

    def rel(x):
        return '--' if x is None or not np.isfinite(x) else (f'{x:.2f}' if x > 0 else '$\\le0$')

    rows = []
    for name in ITEMS:
        s = summary[name]
        p = s.get('projection')
        rows.append(f"{num[name]} & {SHORT[name]} & {pct(s['random']['fraction'])} & {rel(s['random']['max_rel'])} & "
                    f"{pct(p['fraction']) if p else '--'} & {rel(p['max_rel']) if p else '--'} \\\\")
    tex = ['% generated by scripts/trace_montecarlo.py -- do not edit by hand',
           f"\\newcommand{{\\TMCSamples}}{{{args.samples}}}", f"\\newcommand{{\\TMCn}}{{{args.n}}}",
           '\\newcommand{\\TMCRows}{%', *rows, '}']
    (gen / 'trace_montecarlo.tex').write_text('\n'.join(tex) + '\n')

    figdir = OUT / 'figures'
    figdir.mkdir(exist_ok=True)
    grey, orange = '#9aa7b4', '#e8a33d'

    # figure 1: all items against Delta(W)
    fig, axes = plt.subplots(5, 5, figsize=(12, 12.5))
    for ax, name in zip(axes.flat, ITEMS):
        r = res[name]
        m = min(4000, len(r['random']['v']))
        ax.scatter(r['random']['dW'][:m], r['random']['v'][:m], s=1, color=grey, alpha=0.4, rasterized=True)
        if 'projection' in r:
            ax.scatter(r['projection']['dW'][:m], r['projection']['v'][:m], s=1, color=orange, alpha=0.4,
                       rasterized=True)
        grid, V, xs, ref = envelope_curve(name, args.n)
        ax.plot(grid, V, 'k.', ms=4)
        if ref is not None:
            ax.plot(xs, ref, 'k-', lw=1)
        ax.axhline(0, color='#555', lw=0.5)
        top = max(np.nanmax(V), 1e-3)
        ax.set_ylim(-1.1 * top, 1.25 * top)
        ax.set_title(f"{num[name]}: {SHORT[name]}", fontsize=8)
        ax.tick_params(labelsize=6)
        ax.spines[['top', 'right']].set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel(r'$\Delta(W)$', fontsize=8)
    for ax in axes[:, 0]:
        ax.set_ylabel(r'$(L-R)/N$', fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / 'montecarlo_delta.pdf', metadata={'CreationDate': None}, dpi=200)
    plt.close(fig)

    # figure 2: operator items against the contrast of the operators, coloured by Delta(W)
    fig, axes = plt.subplots(4, 4, figsize=(12, 10.5))
    for ax, name in zip(axes.flat, OPERATOR):
        r = res[name]['random']
        m = min(6000, len(r['v']))
        sc = ax.scatter(r['dop'][:m], r['v'][:m], s=1.5, c=r['dW'][:m], cmap='viridis', vmin=0, vmax=1,
                        rasterized=True)
        grid, V, _, _ = envelope_curve(name, args.n)
        top = max(np.nanmax(V), 1e-3)
        ax.set_ylim(-1.1 * top, 1.25 * top)
        ax.axhline(0, color='#555', lw=0.5)
        ax.set_title(f"{num[name]}: {SHORT[name]}", fontsize=8)
        ax.tick_params(labelsize=6)
        ax.spines[['top', 'right']].set_visible(False)
    for ax in axes.flat[len(OPERATOR):]:
        ax.axis('off')
    for ax in axes[-1]:
        ax.set_xlabel("contrast of the operators", fontsize=8)
    cax = fig.add_axes([0.80, 0.08, 0.15, 0.015])
    fig.colorbar(sc, cax=cax, orientation='horizontal').set_label(r'$\Delta(W)$', fontsize=8)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(figdir / 'montecarlo_operator.pdf', metadata={'CreationDate': None}, dpi=200)
    plt.close(fig)

    # figure 3: fraction of violating samples against Delta(W), random against projection samples
    bins = np.linspace(0, 1, 11)
    mid = 0.5 * (bins[1:] + bins[:-1])
    fig, axes = plt.subplots(4, 4, figsize=(12, 10))
    for ax, name in zip(axes.flat, OPERATOR):
        for mode, col, lab in (('random', grey, 'random operators'), ('projection', orange, 'scaled projections')):
            if mode not in res[name]:
                continue
            r = res[name][mode]
            idx = np.digitize(r['dW'], bins) - 1
            frac = [np.mean(r['v'][idx == b] > 1e-12) if np.any(idx == b) else np.nan for b in range(len(mid))]
            ax.plot(mid, frac, 'o-', color=col, ms=3, lw=1.2, label=lab)
        ax.set_ylim(-0.02, 1.02)
        ax.set_title(f"{num[name]}: {SHORT[name]}", fontsize=8)
        ax.tick_params(labelsize=6)
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(color='#e8e8e8', lw=0.5)
    axes.flat[0].legend(frameon=False, fontsize=7)
    for ax in axes.flat[len(OPERATOR):]:
        ax.axis('off')
    for ax in axes[-1]:
        ax.set_xlabel(r'$\Delta(W)$', fontsize=8)
    for ax in axes[:, 0]:
        ax.set_ylabel('fraction violating', fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / 'montecarlo_rate.pdf', metadata={'CreationDate': None})
    plt.close(fig)
    print('wrote', gen / 'trace_montecarlo.json', 'and three figures')


if __name__ == '__main__':
    main()
