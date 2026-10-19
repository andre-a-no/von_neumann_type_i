#!/usr/bin/env python3
"""
Noise-aware, post-selected state preparation in a noisy XY register (use case of man/article, Section "Noise-aware
control with post-selection", sec:control; continues scripts/noisy_xy_mixer.py).

A ring of n qubits starts in |1^w 0^{n-w}>, w = n/2, and is steered towards the Dicke state |D^n_w> (the uniform
superposition of the weight-w sector, a standard initial state of XY-QAOA) by S piecewise-constant segments of
duration T/S of

  H_s = sum_<ij> J_{ij,s} (S^+_i S^-_j + S^-_i S^+_j)/2 + sum_i h_{i,s} S^z_i + V_s sum_<ij> S^z_i S^z_j,

|J|, |h|, |V| <= 1 and 0 < T <= T_max = 2n (the bond terms are SpinChain.bond). Every control keeps the weight, so the
register carries (+)_N M_binom(n, N) and the noise of noisy_xy_mixer.py (S^-_i, S^+_i, S^z_i as SectorJumps) acts on
it. The rescaled evolution (unit time per segment, H_s -> (T/S) H_s, jump coefficients sqrt(gamma T/S)) makes the
final state differentiable in all controls and in T. Three designs, each by batched gradient ascent
(optimize.extremize, one start per batch element):

  (i)   noiseless:  maximise F_0 - LAMBDA T / T_max without noise, F_0 = <D|rho(T)|D> (multistart): the fastest
                    noiseless solution, the most favourable noise-unaware choice;
  (ii)  raw:        maximise the mean over a batch of noise settings of F = <D|rho(T)|D>;
  (iii) post:       maximise the mean of the post-selected fidelity F_post = F / p, p = Tr(P_w rho(T)) the acceptance
                    of the weight measurement (a central projection, DensityMatrix.condition_on).

(ii) and (iii) start from the R best noiseless solutions (and (iii) also from the solution of (ii)), so the search is a
local, noise-aware refinement. The noise batch is 4 models x 3 rates (one Lindblad run for all settings and starts). Every reported number comes
from an independent forward evaluation (physical time, rates passed as `rates`, EVAL_SUBSTEPS steps per segment),
compared with the optimiser's values; p is compared with the birth-death process of Proposition (a) of the paper
(it depends on T and the rates only), one point with lindblad_evolve(method='exact'), and the gradient once with
central differences. Writes man/article/figures/noisy_control.pdf and man/article/generated/noisy_control.{json,tex}.

    python scripts/noisy_xy_control.py [--dims 4,6] [--segments 4] [--substeps 2]
"""
import argparse
import json
import math
import time
import warnings
from pathlib import Path

import numpy as np
import torch

from torch_vn_algebra import DensityMatrix, SpinChain, dynamics, optimize, states
from torch_vn_algebra.cost import CostWarning

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'man' / 'article'
MODELS = {                                   # (g_down, g_up, g_phi) per unit rate g, as in noisy_xy_mixer.py
    'T1': (1.0, 0.0, 0.0),
    'thermal': (1.0, 0.2, 0.0),
    'dephasing': (0.0, 0.0, 1.0),
    'mixed': (1.0, 0.1, 1.0),
}
RATES = (0.003, 0.01, 0.03)                  # training batch: every model at every rate
GRID = np.logspace(-3, -1, 9)                # rates of the figure (evaluation only)
KINDS = ('lower', 'raise', 'sz')
DESIGNS = ('noiseless', 'raw', 'post')
F0_MIN = 0.99                                # a noiseless solution
LAMBDA = 0.05                                # time penalty of the noiseless design: F_0 - LAMBDA T / T_max
EVAL_SUBSTEPS = 20
COLORS = {'noiseless': '#7d848c', 'raw': '#eb6834', 'post': '#2a78d6'}
LABELS = {'noiseless': '(i) noiseless design', 'raw': '(ii) noise-aware, raw $F$',
          'post': r'(iii) noise-aware, $F_{\mathrm{post}}$'}
INK, GREY = '#0b0b0b', '#7d848c'
NAMES = {4: 'Four', 5: 'Five', 6: 'Six', 7: 'Seven', 8: 'Eight'}


class Register:
    """The chain, its control Hamiltonians, the jumps, |1^w 0^{n-w}>, the Dicke state and the projection P_w."""

    def __init__(self, n: int, segments: int):
        self.n, self.w, self.S, self.T_max = n, n // 2, segments, 2.0 * n
        q = SpinChain(n, 'periodic', precision='double')
        self.q, self.alg = q, q.algebra
        self.hop = torch.stack([q.bond(i, j, 1.0, 0.0).matrix[0] for i, j in q.bonds])    # (bonds, C, k, k)
        self.field = torch.stack([q.sz(i).matrix[0] for i in range(n)])                  # (n, C, k, k)
        self.zz = sum(q.bond(i, j, 0.0, 1.0).matrix[0] for i, j in q.bonds)              # (C, k, k)
        self.jumps = {kind: [q.sector_jump(i, kind) for i in range(n)] for kind in KINDS}
        self.rho0 = q.basis_state('1' * self.w + '0' * (n - self.w))
        self.P = self.alg.central([float(N == self.w) for N in range(n + 1)])
        k = math.comb(n, self.w)
        D = torch.zeros(1, self.alg.C, self.alg.k_max, self.alg.k_max, dtype=self.alg.dtype)
        D[0, self.w, :k, :k] = 1.0 / k                                                    # |D><D|, D = sum |s>/sqrt(k)
        self.dicke = DensityMatrix(self.alg, matrix=D)

    def hamiltonians(self, J, h, V) -> torch.Tensor:
        """H_s for controls J (R, S, bonds), h (R, S, n), V (R, S): shape (R, S, C, k, k)."""
        dt = self.hop.dtype
        return (torch.einsum('rsb,bcij->rscij', J.to(dt), self.hop) + torch.einsum('rsl,lcij->rscij', h.to(dt), self.field)
                + V.to(dt)[..., None, None, None] * self.zz)

    def measure(self, rho):
        """(F, F_post, p) per batch element."""
        post, p = rho.condition_on(self.P)
        return (states.trace_of_product(rho, self.dicke).real, states.trace_of_product(post, self.dicke).real, p)


class Controls(optimize.Param):
    """Bounded controls of R starts: J = tanh(x), h = tanh(z), V = tanh(y), T = T_max sigmoid(t)."""

    def __init__(self, raw):
        self.raw = {k: v.detach().clone().requires_grad_() for k, v in raw.items()}

    @classmethod
    def random(cls, R, S, n, gen):
        r = lambda *shape: 0.5 * torch.randn(*shape, generator=gen, dtype=torch.float64)
        return cls(dict(x=r(R, S, n), z=r(R, S, n), y=r(R, S), t=r(R)))

    def parameters(self):
        return list(self.raw.values())

    def values(self, T_max):
        r = self.raw
        return torch.tanh(r['x']), torch.tanh(r['z']), torch.tanh(r['y']), T_max * torch.sigmoid(r['t'])

    def select(self, idx):
        return Controls({k: v[idx] for k, v in self.raw.items()})


def noise_batch(rates=RATES, models=MODELS):
    """(K, 3) rates (g_down, g_up, g_phi) and their (model, rate) labels."""
    combos = [(m, g) for m in models for g in rates]
    return torch.tensor([[c * g for c in MODELS[m]] for m, g in combos], dtype=torch.float64), combos


def evolve(reg, ctrl, noise, substeps):
    """rho(T) for every start r and noise setting k (batch index r K + k), each segment rescaled to unit time:
    T enters through H and the jump coefficients sqrt(gamma T / S), so rho(T) is differentiable in T."""
    J, h, V, T = ctrl.values(reg.T_max)
    K = noise.shape[0]
    tau = T / reg.S
    H = (reg.hamiltonians(J, h, V) * tau.to(reg.hop.dtype)[:, None, None, None, None]).repeat_interleave(K, 0)
    jumps = []
    for j, kind in enumerate(KINDS):
        if bool((noise[:, j] > 0).any()):
            coef = (noise[None, :, j].sqrt() * tau[:, None].sqrt()).reshape(-1)          # (R K,)
            jumps += [L.scaled(coef) for L in reg.jumps[kind]]
    rho = reg.rho0
    for s in range(reg.S):
        Hs = reg.alg.operator(H[:, s], is_self_adjoint=True)
        rho = dynamics.lindblad_evolve(rho, Hs, jumps, [0.0, 1.0], substeps=substeps, method='kraus')[-1]
    return rho


def objective_values(reg, ctrl, noise, design, substeps):
    """The objective of a design per start (R,), and (F, F_post, p) of shape (R, K). The noiseless design (one
    noise-free setting) maximises F_0 - LAMBDA T / T_max: the fastest noiseless solution."""
    F, Fp, p = (x.reshape(-1, noise.shape[0]) for x in reg.measure(evolve(reg, ctrl, noise, substeps)))
    if design == 'noiseless':
        return F[:, 0] - LAMBDA * ctrl.values(reg.T_max)[3] / reg.T_max, (F, Fp, p)
    return (Fp if design == 'post' else F).mean(-1), (F, Fp, p)


def optimise(reg, ctrl, noise, design, substeps, rounds, steps, lr=0.05, decay=0.5):
    optimize.extremize(lambda: objective_values(reg, ctrl, noise, design, substeps)[0], [ctrl], maximize=True,
                       rounds=rounds, steps=steps, lr=lr, decay=decay)
    with torch.no_grad():
        return objective_values(reg, ctrl, noise, design, substeps)


def evaluate(reg, ctrl, g, substeps, method='kraus'):
    """Independent forward evaluation of one design (R = 1) under rates g = (g_down, g_up, g_phi): physical time
    [0, T/S] per segment, unit jumps with the rates passed as `rates`. Returns (F, F_post, p) as floats."""
    with torch.no_grad():
        J, h, V, T = ctrl.values(reg.T_max)
        H = reg.hamiltonians(J, h, V)[0]
        jumps, rates = [], []
        for kind, gk in zip(KINDS, g):
            if gk > 0:
                jumps += reg.jumps[kind]
                rates += [float(gk)] * reg.n
        rho = reg.rho0
        for s in range(reg.S):
            rho = dynamics.lindblad_evolve(rho, reg.alg.operator(H[s:s + 1], is_self_adjoint=True), jumps,
                                           [0.0, float(T.detach()) / reg.S], rates=rates or None, substeps=substeps,
                                           method=method)[-1]
        return tuple(float(x[0]) for x in reg.measure(rho))


def birth_death_acceptance(n, w, g, T):
    """p_w(T) of the birth-death process of Proposition (a): death rate g_down N, birth rate g_up (n - N)."""
    gd, gu = g[0], g[1]
    Q = np.zeros((n + 1, n + 1))
    for N in range(n + 1):
        if N > 0:
            Q[N - 1, N] += gd * N
            Q[N, N] -= gd * N
        if N < n:
            Q[N + 1, N] += gu * (n - N)
            Q[N, N] -= gu * (n - N)
    p0 = np.zeros(n + 1)
    p0[w] = 1.0
    return float((torch.linalg.matrix_exp(torch.tensor(Q * T)) @ torch.tensor(p0))[w])


def gradient_check(reg, ctrl, noise, substeps, gen, eps=1e-6):
    """Relative difference of the autograd directional derivative of the (iii) objective and central differences."""
    v = {k: torch.randn(t.shape, generator=gen, dtype=torch.float64) for k, t in ctrl.raw.items()}

    def f(shift):
        c = Controls({k: t.detach() + shift * v[k] for k, t in ctrl.raw.items()})
        return objective_values(reg, c, noise, 'post', substeps)[0].sum()
    c = Controls({k: t.detach() for k, t in ctrl.raw.items()})
    val = objective_values(reg, c, noise, 'post', substeps)[0].sum()
    grads = torch.autograd.grad(val, c.parameters())
    auto = sum((g_ * v[k]).sum() for g_, k in zip(grads, c.raw)).item()
    with torch.no_grad():
        fd = ((f(eps) - f(-eps)) / (2 * eps)).item()
    return abs(auto - fd) / abs(fd)


def run(n, args, gen):
    reg = Register(n, args.segments)
    noise, combos = noise_batch()
    starts, warm = (args.starts, args.warm) if n <= 4 else (args.starts, args.warm_large)
    timing = {}
    t0 = time.time()
    # (i) noiseless multistart (no jumps: closed evolution through the same solver)
    ctrl0 = Controls.random(starts, args.segments, n, gen)
    obj0, (F0, _, _) = optimise(reg, ctrl0, torch.zeros(1, 3, dtype=torch.float64), 'noiseless', args.substeps,
                                rounds=4, steps=args.steps0)
    F0 = F0[:, 0]
    T0 = ctrl0.values(reg.T_max)[3].detach()
    order = torch.argsort(obj0, descending=True)
    chosen = order[:warm]
    designs = {'noiseless': ctrl0.select(order[:1])}
    timing['noiseless'] = time.time() - t0
    # (ii), (iii): noise-aware refinement of the `warm` fastest noiseless solutions, all noise settings in one batch
    values = {}
    for d in ('raw', 'post'):
        t1 = time.time()
        ctrl = ctrl0.select(chosen)
        if d == 'post':                                        # the solution of (ii) is one more start of (iii)
            ctrl = Controls({k: torch.cat([v.detach(), designs['raw'].raw[k].detach()]) for k, v in ctrl.raw.items()})
        obj, _ = optimise(reg, ctrl, noise, d, args.substeps, rounds=args.rounds, steps=args.steps)
        best = int(torch.argmax(obj))
        designs[d] = ctrl.select([best])
        values[d] = obj.tolist()
        timing[d] = time.time() - t1
    # independent evaluation of every reported number
    t1 = time.time()
    res, dev, bd = {}, 0.0, 0.0
    for d, c in designs.items():
        T = float(c.values(reg.T_max)[3].detach())
        with torch.no_grad():
            _, (Fo, Fpo, po) = objective_values(reg, c, noise, 'post', args.substeps)    # the optimiser's values
        rows = [evaluate(reg, c, noise[i].tolist(), EVAL_SUBSTEPS) for i in range(len(combos))]
        F, Fp, p = (np.array(x) for x in zip(*rows))
        dev = max(dev, float(np.max(np.abs(np.stack([F, Fp, p]) - torch.stack([Fo[0], Fpo[0], po[0]]).numpy()))))
        bd = max(bd, max(abs(p[i] - birth_death_acceptance(n, reg.w, noise[i].tolist(), T)) for i in range(len(combos))))
        F_0 = evaluate(reg, c, (0.0, 0.0, 0.0), EVAL_SUBSTEPS)[0]
        J, h, V, _ = (x.detach()[0] for x in c.values(reg.T_max))
        res[d] = dict(T=T, F0=F_0, F=F.tolist(), Fpost=Fp.tolist(), p=p.tolist(), meanF=float(F.mean()),
                      meanFpost=float(Fp.mean()), meanp=float(p.mean()), minF=float(F.min()), minFpost=float(Fp.min()),
                      J=J.tolist(), h=h.tolist(), V=V.tolist())
    # the figure: every model on a finer grid of rates (batched, evaluation only)
    gnoise, gcombos = noise_batch(GRID.tolist())
    for d, c in designs.items():
        with torch.no_grad():
            _, (F, Fp, p) = objective_values(reg, c, gnoise, 'post', EVAL_SUBSTEPS // 2)
        res[d]['grid'] = {m: dict(F=[F[0, i].item() for i, (mm, _) in enumerate(gcombos) if mm == m],
                                  Fpost=[Fp[0, i].item() for i, (mm, _) in enumerate(gcombos) if mm == m],
                                  p=[p[0, i].item() for i, (mm, _) in enumerate(gcombos) if mm == m]) for m in MODELS}
    timing['evaluation'] = time.time() - t1
    # one point with the exact propagator exp(h L), and one gradient check
    t1 = time.time()
    i_top = len(combos) - 1                                    # mixed model, largest rate
    exact = evaluate(reg, designs['post'], noise[i_top].tolist(), 1, method='exact')
    kraus = (res['post']['F'][i_top], res['post']['Fpost'][i_top], res['post']['p'][i_top])
    exact_diff = max(abs(a - b) for a, b in zip(exact, kraus))
    timing['exact'] = time.time() - t1
    grad_err = gradient_check(reg, designs['post'], noise, args.substeps, gen) if n == min(args.dims) else None
    info = dict(n=n, w=reg.w, sector=math.comb(n, reg.w), algebra=sum(k * k for k in reg.alg.k_factors),
                starts=starts, warm=int(len(chosen)), noiseless_solutions=int((F0 >= F0_MIN).sum()),
                F0_starts=F0.tolist(), T_starts=T0.tolist(), objective_values=values, recheck=dev, birth_death=bd,
                exact=list(exact), exact_point=list(combos[i_top]), exact_diff=exact_diff, grad_rel_err=grad_err,
                seconds=timing, T_max=reg.T_max, batch=len(combos) * int(len(chosen)))
    return info, res


def figure(results):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    n = max(results)
    info, res = results[n]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.3))
    ax = axes[0]
    for d in DESIGNS:
        g = res[d]['grid']['mixed']
        ax.plot(GRID, 1 - np.array(g['Fpost']), '-o', color=COLORS[d], lw=1.6, ms=3.5, mfc='white')
        ax.plot(GRID, 1 - np.array(g['F']), '--', color=COLORS[d], lw=1.1)
    for r in RATES:
        ax.axvline(r, color='#e6e5e1', lw=3, zorder=0)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'rate $\gamma$ (mixed noise: $\gamma_\downarrow=\gamma_\varphi=\gamma$, $\gamma_\uparrow=0.1\gamma$)',
                  fontsize=8, color=INK)
    ax.set_ylabel('infidelity with $|D^{%d}_{%d}\\rangle$' % (n, info['w']), fontsize=8, color=INK)
    ax.set_title(f'(a) $n={n}$: solid after post-selection, dashed without', fontsize=9, color=INK)
    ax = axes[1]
    ref = res['noiseless']['grid']['mixed']
    ax.axhline(1.0, color=COLORS['noiseless'], lw=1.6)
    for d in ('raw', 'post'):
        g = res[d]['grid']['mixed']
        ax.plot(GRID, (1 - np.array(g['Fpost'])) / (1 - np.array(ref['Fpost'])), '-o', color=COLORS[d], lw=1.6,
                ms=3.5, mfc='white')
        ax.plot(GRID, (1 - np.array(g['F'])) / (1 - np.array(ref['F'])), '--', color=COLORS[d], lw=1.1)
    for r in RATES:
        ax.axvline(r, color='#e6e5e1', lw=3, zorder=0)
    ax.set_xscale('log')
    ax.set_xlabel(r'rate $\gamma$ (mixed noise)', fontsize=8, color=INK)
    ax.set_ylabel('infidelity relative to design (i)', fontsize=8, color=INK)
    ax.set_title('(b) gain of the noise-aware designs', fontsize=9, color=INK)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
        ax.spines[['left', 'bottom']].set_color(GREY)
        ax.tick_params(labelsize=7, colors=GREY, labelcolor=INK)
        ax.grid(color='#eeeeec', lw=0.6)
    handles = [plt.Line2D([], [], color=COLORS[d], lw=1.6,
                          label=LABELS[d] + f', $T={res[d]["T"]:.2f}$, ' + r'$\bar p=' + f'{res[d]["meanp"]:.2f}$')
               for d in DESIGNS]
    fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, fontsize=7.5)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    (OUT / 'figures').mkdir(exist_ok=True)
    fig.savefig(OUT / 'figures' / 'noisy_control.pdf', metadata={'CreationDate': None})
    plt.close(fig)


def sci(x):
    m, e = f'{x:.0e}'.split('e')
    return f'{m}\\times10^{{{int(e)}}}'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dims', default='4,6')
    ap.add_argument('--segments', type=int, default=4)
    ap.add_argument('--substeps', type=int, default=2, help='Kraus steps per segment during the optimisation')
    ap.add_argument('--starts', type=int, default=32, help='noiseless multistart')
    ap.add_argument('--warm', type=int, default=8, help='noise-aware starts for n <= 4')
    ap.add_argument('--warm-large', type=int, default=4, help='noise-aware starts for n > 4')
    ap.add_argument('--steps0', type=int, default=100, help='steps per round, noiseless (4 rounds)')
    ap.add_argument('--rounds', type=int, default=3)
    ap.add_argument('--steps', type=int, default=50)
    ap.add_argument('--seed', type=int, default=1)
    args = ap.parse_args()
    args.dims = [int(x) for x in args.dims.split(',')]
    warnings.simplefilter('ignore', CostWarning)
    torch.set_num_threads(max(1, torch.get_num_threads()))
    gen = torch.Generator().manual_seed(args.seed)
    t_all = time.time()
    results, rows, tex = {}, [], []
    for n in args.dims:
        info, res = run(n, args, gen)
        results[n] = (info, res)
        s = info['seconds']
        print(f"n={n}: {sum(s.values()):.0f} s ({', '.join(f'{k} {v:.0f}' for k, v in s.items())}); "
              f"recheck {info['recheck']:.1e}, birth-death {info['birth_death']:.1e}, exact {info['exact_diff']:.1e}, "
              f"grad {info['grad_rel_err']}")
        for d in DESIGNS:
            r = res[d]
            print(f"  {d:9s} T={r['T']:.3f} F0={r['F0']:.4f}  mean F={r['meanF']:.4f}  mean Fpost={r['meanFpost']:.4f} "
                  f"mean p={r['meanp']:.4f}  min Fpost={r['minFpost']:.4f}  min F={r['minF']:.4f}")
            rows.append(f"{n} & {LABELS[d].split(' ')[0]} & {r['T']:.2f} & {r['F0']:.4f} & {r['meanp']:.3f} & "
                        f"{r['meanF']:.3f} & {r['meanFpost']:.3f} & {r['minF']:.3f} & {r['minFpost']:.3f}\\\\")
        N = NAMES[n]
        ri, rii, riii = res['noiseless'], res['raw'], res['post']
        ratio = lambda i: ((1 - riii['grid']['mixed']['Fpost'][i]) / (1 - ri['grid']['mixed']['Fpost'][i]))
        tex += [f'\\newcommand{{\\NC{N}MinFPost{k}}}{{{res[d]["minFpost"]:.3f}}}' for k, d in zip(('One', 'Two', 'Three'), DESIGNS)]
        tex += [f'\\newcommand{{\\NC{N}T{k}}}{{{res[d]["T"]:.2f}}}' for k, d in zip(('One', 'Two', 'Three'), DESIGNS)]
        tex += [f'\\newcommand{{\\NC{N}FPost{k}}}{{{res[d]["meanFpost"]:.3f}}}' for k, d in zip(('One', 'Two', 'Three'), DESIGNS)]
        tex += [f'\\newcommand{{\\NC{N}F{k}}}{{{res[d]["meanF"]:.3f}}}' for k, d in zip(('One', 'Two', 'Three'), DESIGNS)]
        tex += [f'\\newcommand{{\\NC{N}P{k}}}{{{res[d]["meanp"]:.3f}}}' for k, d in zip(('One', 'Two', 'Three'), DESIGNS)]
        tex += [f'\\newcommand{{\\NC{N}FZero{k}}}{{{res[d]["F0"]:.4f}}}' for k, d in zip(('One', 'Two', 'Three'), DESIGNS)]
        tex += [f'\\newcommand{{\\NC{N}GainPost}}{{{(1 - ri["meanFpost"]) / (1 - riii["meanFpost"]):.2f}}}',
                f'\\newcommand{{\\NC{N}GainRaw}}{{{(1 - ri["meanF"]) / (1 - rii["meanF"]):.2f}}}',
                f'\\newcommand{{\\NC{N}GainTotal}}{{{(1 - rii["meanF"]) / (1 - riii["meanFpost"]):.2f}}}',
                f'\\newcommand{{\\NC{N}PostFactor}}{{{(1 - ri["meanF"]) / (1 - ri["meanFpost"]):.1f}}}',
                f'\\newcommand{{\\NC{N}RatioLow}}{{{ratio(0):.2f}}}',
                f'\\newcommand{{\\NC{N}RatioHigh}}{{{ratio(-1):.2f}}}',
                f'\\newcommand{{\\NC{N}Solutions}}{{{info["noiseless_solutions"]}}}',
                f'\\newcommand{{\\NC{N}Seconds}}{{{sum(s.values()):.0f}}}']
    seconds = time.time() - t_all
    figure(results)
    gen_dir = OUT / 'generated'
    gen_dir.mkdir(exist_ok=True)
    (gen_dir / 'noisy_control.json').write_text(json.dumps(dict(
        models=MODELS, rates=RATES, grid=GRID.tolist(), args=vars(args), seconds=seconds, threads=torch.get_num_threads(),
        results={n: dict(info=i, designs=r) for n, (i, r) in results.items()}), indent=1))
    infos = [i for i, _ in results.values()]
    tex = ['% generated by scripts/noisy_xy_control.py -- do not edit by hand',
           '\\newcommand{\\NCRows}{%', *rows, '}', *tex,
           f'\\newcommand{{\\NCSegments}}{{{args.segments}}}',
           f'\\newcommand{{\\NCSubsteps}}{{{args.substeps}}}',
           f'\\newcommand{{\\NCEvalSubsteps}}{{{EVAL_SUBSTEPS}}}',
           f'\\newcommand{{\\NCStarts}}{{{args.starts}}}',
           f'\\newcommand{{\\NCWarm}}{{{args.warm}}}',
           f'\\newcommand{{\\NCWarmLarge}}{{{args.warm_large}}}',
           f'\\newcommand{{\\NCSettings}}{{{len(MODELS) * len(RATES)}}}',
           f'\\newcommand{{\\NCSteps}}{{{args.rounds * args.steps}}}',
           f'\\newcommand{{\\NCStepsZero}}{{{4 * args.steps0}}}',
           f'\\newcommand{{\\NCFZeroMin}}{{{F0_MIN}}}',
           f'\\newcommand{{\\NCLambda}}{{{LAMBDA:g}}}',
           f'\\newcommand{{\\NCRecheck}}{{{sci(max(i["recheck"] for i in infos))}}}',
           f'\\newcommand{{\\NCBirthDeath}}{{{sci(max(i["birth_death"] for i in infos))}}}',
           f'\\newcommand{{\\NCExactDiff}}{{{sci(max(i["exact_diff"] for i in infos))}}}',
           f'\\newcommand{{\\NCGradErr}}{{{sci(min(i["grad_rel_err"] for i in infos if i["grad_rel_err"] is not None))}}}',
           f'\\newcommand{{\\NCSeconds}}{{{seconds:.0f}}}',
           f'\\newcommand{{\\NCMinutes}}{{{seconds / 60:.0f}}}',
           f'\\newcommand{{\\NCThreads}}{{{torch.get_num_threads()}}}',
           f'\\newcommand{{\\NCnMax}}{{{max(results)}}}']
    (gen_dir / 'noisy_control.tex').write_text('\n'.join(tex) + '\n')
    print(f'total {seconds:.0f} s; wrote figures/noisy_control.pdf and generated/noisy_control.{{json,tex}}')


if __name__ == '__main__':
    main()
