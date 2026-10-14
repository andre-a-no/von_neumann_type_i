#!/usr/bin/env python3
"""
Section "Rediscovering the data-processing thresholds of the Renyi divergences" of paper.

The sandwiched Renyi divergences satisfy the data-processing inequality (DPI) ~D_alpha(Phi rho||Phi sigma) <=
~D_alpha(rho||sigma) for every channel Phi exactly when alpha >= 1/2 (Frank-Lieb 2013, Beigi 2013, Mueller-Lennert
et al. 2013, Wilde-Winter-Yang 2014), the Petz divergences D_alpha for alpha in (0, 2] (Petz 1986), and both fail in
general outside these ranges. Without using this knowledge we maximise the violation

    v(alpha) = D(Phi rho||Phi sigma) - D(rho||sigma)

on a grid of alpha in [0.2, 3] by batched multistart gradient ascent (optimize.extremize, one start per batch element,
all alpha of the grid in one batch through per-sample alpha), over states rho, sigma on M_d and over two channel
families: the complete dephasing (the conditional expectation onto the diagonal maximal abelian subalgebra; fixed)
and all channels with a d-dimensional environment (Stinespring isometries V = Z (Z^* Z)^{-1/2}, Z a free complex
dr x d matrix). The states are parametrised as rho = (1 - eps) G G^* / Tr(G G^*) + eps 1/d (G free complex d x d),
i.e. over the compact set {rho >= eps/d}: the supremum is finite and the sandwich sigma^g rho sigma^g,
g = (1 - alpha)/(2 alpha), stays resolved in double precision for alpha >= 0.2 (its condition number is at most
(d/eps)^(2|g|+1); for smaller alpha eigenvalues at the rounding level would contribute lambda^alpha and the computed
divergence would be meaningless). The brackets [last violating alpha, first non-violating alpha] of both thresholds
are then refined by further searches at interior points: random starts as before plus starts at small perturbations
of the best point of the violating end of the bracket (continuation in alpha; the violation vanishes continuously at
the threshold, so random starts alone often miss the small violations next to it). Every reported value v is attained
at a known point (the final iterate of the best start) and is recomputed independently with mpmath at 50 significant
digits (dense matrices, eigendecompositions in mpmath; the Kraus operators re-orthonormalised in mpmath, so that the
certified point is an exactly trace-preserving channel); a violation counts only if the certified value exceeds TOL,
which is far above the error of the certification: a certified violation is a proof that the DPI fails at that
alpha. The other end of a bracket (no violation above TOL found) is evidence only: a local search can miss
violations, and violations below TOL are not counted.

Writes generated/renyi_threshold.tex, generated/renyi_threshold.json and figures/renyi_threshold.pdf.
"""
import math
import time
import warnings

import mpmath as mp
import numpy as np
import torch

from common import parse_args, save_json, write_tex, env_macro, sci, load_partial, save_partial, finish_partial

from torch_vn_algebra import TypeIAlgebra, Channel, DensityMatrix, channels, cost, divergences as dv, optimize as opt

args = parse_args(__doc__)
dev = args.device
FULL = args.mode in ('full', 'check')            # check: the sizes of full, minimal repetitions
CHECK = args.mode == 'check'
DIMS = (2, 3) if FULL else (2,)
STARTS = 64 if FULL else 16                      # independent starts per value of alpha
ROUNDS, STEPS, LR, DECAY = (1, 4, 0.1, 0.6) if CHECK else (8, 75, 0.1, 0.6) if FULL else (5, 60, 0.1, 0.6)
REFINE_ROUNDS, REFINE_POINTS = (1, 3) if CHECK else (3, 7)
INIT_NOISE, WARM_LR = 0.02, 0.2                  # continuation starts: relative perturbation, learning rate factor
EPS = 0.01                                       # states rho >= EPS / d
TOL = 1e-12                                      # a certified v above TOL is a violation (certification error ~1e-45)
DIGITS = 50                                      # mpmath precision of the certification
ALPHAS = [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.48, 0.49, 0.5, 0.51, 0.52, 0.55, 0.6, 0.7, 0.8, 0.9, 1.1, 1.25, 1.5,
          1.75, 1.9, 1.95, 1.98, 2.0, 2.02, 2.05, 2.1, 2.25, 2.5, 2.75, 3.0]
FAMILIES = ('dephasing', 'stinespring')
KINDS = ('sandwiched', 'petz')
KNOWN = {'sandwiched': 0.5, 'petz': 2.0}         # DPI holds for alpha >= 1/2 (sandwiched), alpha <= 2 (Petz)
warnings.simplefilter('ignore', cost.CostWarning)     # run-time estimates of extremize (expected, not an issue)


# ----------------------------------------------------------------------------------------------------------------
# parametrisations (cheap on the CPU: no matrix exponentials; optimize.Param interface for optimize.extremize)
# ----------------------------------------------------------------------------------------------------------------
class StateParam(opt.Param):
    """rho = (1 - EPS) G G^* / Tr(G G^*) + EPS 1/d, one state per start (G complex Ginibre at the start)."""

    def __init__(self, alg, batch):
        self.alg = alg
        self.G = torch.randn(batch, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype, device=alg.hilbert.device,
                             requires_grad=True)

    def parameters(self):
        return [self.G]

    def state(self) -> DensityMatrix:
        P = self.G @ self.G.conj().transpose(-2, -1)
        P = P / torch.diagonal(P, dim1=-2, dim2=-1).sum((-2, -1)).real[:, None, None, None].to(P.dtype)
        d = self.alg.total_subspace_dim
        return DensityMatrix(self.alg, matrix=(1 - EPS) * P + (EPS / d) * self.alg._const('eye'), validate=False)


class StinespringParam(opt.Param):
    """Channels of M_d with an r-dimensional environment: isometry V = Z (Z^* Z)^{-1/2}: C^d -> C^d (x) C^r and Kraus
    operators K_j = (1 (x) <j|) V."""

    def __init__(self, alg, batch, r):
        self.alg, self.r = alg, r
        d = alg.k_max
        self.Z = torch.randn(batch, d * r, d, dtype=alg.dtype, device=alg.hilbert.device, requires_grad=True)

    def parameters(self):
        return [self.Z]

    def kraus(self) -> torch.Tensor:
        Z, d, r, B = self.Z, self.alg.k_max, self.r, self.Z.shape[0]
        gram = self.alg.operator((Z.conj().transpose(-2, -1) @ Z).unsqueeze(1), is_self_adjoint=True,
                                 is_positive=True)
        V = Z @ gram.power(-0.5).matrix[:, 0]                                    # (B, d r, d), V^* V = 1
        return V.reshape(B, d, r, d).permute(0, 2, 1, 3).unsqueeze(2)          # (B, r, 1, d, d)

    def channel(self) -> Channel:
        return Channel(self.alg, self.kraus(), trace_preserving=True)          # trace preserving by construction


DIVERGENCE = {'sandwiched': dv.sandwiched_renyi_divergence, 'petz': dv.petz_renyi_divergence}


def search(d, family, kind, alphas, seed, starts=STARTS, lr=LR, init=None):
    """Maximise v over states (and channels) for every alpha in `alphas` (`starts` starts each, one batch): random
    starts, or with init (the parameter tensors of one point) the first start at init and the others at perturbations
    of it. Returns per alpha: the best value over all iterates, the value at the final iterate of the best start, and
    that point (matrices and parameters)."""
    torch.manual_seed(seed)
    alg = TypeIAlgebra([d], [d], precision='double', device=dev)
    A = torch.tensor(alphas, dtype=torch.float64, device=dev).repeat_interleave(starts)
    B = A.shape[0]
    X, Y = StateParam(alg, B), StateParam(alg, B)
    params = [X, Y]
    if family == 'dephasing':
        deph = channels.dephasing_channel(alg, 1.0)
        get_channel = lambda: deph
    else:
        W = StinespringParam(alg, B, r=d)
        params.append(W)
        get_channel = W.channel
    tensors = [t for p in params for t in p.parameters()]
    if init is not None:                             # continuation: perturbations of a known point
        with torch.no_grad():
            for t, t0 in zip(tensors, init):
                t0 = torch.as_tensor(t0, device=t.device)
                scale = INIT_NOISE * t0.abs().pow(2).mean().sqrt()
                noise = torch.randn_like(t) * scale
                noise[::starts] = 0                  # the first start of every alpha exactly at the point
                t.copy_(t0.expand_as(t) + noise)
    D = DIVERGENCE[kind]

    def objective():
        rho, sigma, Phi = X.state(), Y.state(), get_channel()
        return D(Phi(rho), Phi(sigma), A) - D(rho, sigma, A)

    best = opt.extremize(objective, params, maximize=True, rounds=ROUNDS, steps=STEPS, lr=lr, decay=DECAY)
    with torch.no_grad():
        final = objective()
        rho, sigma = X.state().matrix[:, 0], Y.state().matrix[:, 0]
        K = get_channel().kraus.expand(B, -1, -1, -1, -1)[:, :, 0]
    n = len(alphas)
    best, final = best.reshape(n, starts), final.reshape(n, starts)
    i = torch.nan_to_num(final, nan=-math.inf).argmax(dim=1) + starts * torch.arange(n, device=final.device)
    out = []
    for j, a in enumerate(alphas):
        out.append(dict(alpha=a, best_iterate=float(best[j].max()), value=float(final.reshape(-1)[i[j]]),
                        rho=rho[i[j]].cpu().numpy(), sigma=sigma[i[j]].cpu().numpy(), kraus=K[i[j]].cpu().numpy(),
                        params=[t[i[j]].detach().clone() for t in tensors]))
    return out


# ----------------------------------------------------------------------------------------------------------------
# independent certification in mpmath
# ----------------------------------------------------------------------------------------------------------------
def _mp(m):
    return mp.matrix([[mp.mpc(complex(x)) for x in row] for row in m])


def _mp_power(M, p):
    """M^p of a Hermitian positive definite mpmath matrix (the states here satisfy rho >= EPS/d, also their images)."""
    E, Q = mp.eighe(M)
    return Q * mp.diag([e ** p for e in E]) * Q.transpose_conj()


def _mp_divergence(kind, a, rho, sigma):
    a = mp.mpf(a)
    tr = sum(rho[i, i] for i in range(rho.rows)).real
    if kind == 'petz':
        Q = sum((_mp_power(rho, a) * _mp_power(sigma, 1 - a))[i, i] for i in range(rho.rows)).real
    else:
        S = _mp_power(sigma, (1 - a) / (2 * a))
        M = S * rho * S
        M = (M + M.transpose_conj()) / 2
        Q = sum(e ** a for e in mp.eighe(M)[0])
    return mp.log(Q / tr) / (a - 1)


def certify(kind, a, rho, sigma, kraus):
    """v = D(Phi rho||Phi sigma) - D(rho||sigma) at the given point, with mpmath at DIGITS significant digits."""
    with mp.workdps(DIGITS):
        r, s = _mp(rho), _mp(sigma)
        r, s = (r + r.transpose_conj()) / 2, (s + s.transpose_conj()) / 2
        Ks = [_mp(K) for K in kraus]
        # re-orthonormalise in mpmath, K_j -> K_j S^{-1/2} with S = sum_j K_j^* K_j (= 1 up to double rounding): the
        # certified point is an exactly trace-preserving channel (a CP map with S != 1 could violate the DPI)
        Sm = sum((K.transpose_conj() * K for K in Ks), mp.zeros(Ks[0].cols, Ks[0].cols))
        Si = _mp_power((Sm + Sm.transpose_conj()) / 2, -0.5)
        Ks = [K * Si for K in Ks]
        phi = lambda x: sum((K * x * K.transpose_conj() for K in Ks), mp.zeros(x.rows, x.cols))
        v = _mp_divergence(kind, a, phi(r), phi(s)) - _mp_divergence(kind, a, r, s)
        return float(v)


def run(d, family, kind, alphas, seed, init=None):
    """search (random starts; with init also a continuation search from that point, keeping per alpha the better
    final point), then certification of the final points."""
    t0 = time.time()
    res = search(d, family, kind, alphas, seed)
    if init is not None:
        warm = search(d, family, kind, alphas, seed + 7919, starts=max(2, STARTS // 2), lr=LR * WARM_LR, init=init)
        res = [dict(w if w['value'] > r['value'] else r, best_iterate=max(r['best_iterate'], w['best_iterate']))
               for r, w in zip(res, warm)]
    for r in res:
        r['certified'] = certify(kind, r['alpha'], r['rho'], r['sigma'], r['kraus'])
        r['rho_min_eig'] = float(np.linalg.eigvalsh(r['rho'])[0])
        r['sigma_min_eig'] = float(np.linalg.eigvalsh(r['sigma'])[0])
        for key in ('rho', 'sigma', 'kraus'):          # (the parameters are kept for continuation searches)
            del r[key]
    return res, time.time() - t0


def bracket(points, kind):
    """(lo, hi): the threshold between the last grid point on one side and the first on the other, from the
    certified values; None if the search found no violation (or no non-violating point)."""
    pts = sorted(points, key=lambda r: r['alpha'])
    viol = [r['alpha'] for r in pts if r['certified'] > TOL]
    if not viol:
        return None
    if kind == 'sandwiched':                         # violations expected below the threshold
        lo = max(a for a in viol if a < 1.0) if any(a < 1.0 for a in viol) else None
        if lo is None:
            return None
        above = [r['alpha'] for r in pts if r['alpha'] > lo and r['certified'] <= TOL]
        return (lo, min(above)) if above else None
    hi = min(viol)                                   # Petz: violations expected above the threshold
    below = [r['alpha'] for r in pts if r['alpha'] < hi and r['certified'] <= TOL]
    return (max(below), hi) if below else None


# ----------------------------------------------------------------------------------------------------------------
results = {}
t_start = time.time()
for d in DIMS:
    for fi, family in enumerate(FAMILIES):
        for ki, kind in enumerate(KINDS):
            key = f'd{d}_{family}_{kind}'
            cached = load_partial(args, key)
            if cached is not None:
                results[key] = cached
                continue
            seed = args.seed + 1000 * d + 100 * fi + 10 * ki
            pts, secs = run(d, family, kind, ALPHAS, seed)
            refinements = []
            for rr in range(REFINE_ROUNDS):
                br = bracket(pts, kind)
                if br is None:
                    break
                lo, hi = br
                inner = [lo + (hi - lo) * (j + 1) / (REFINE_POINTS + 1) for j in range(REFINE_POINTS)]
                end = lo if kind == 'sandwiched' else hi           # the violating end of the bracket
                start = max((r for r in pts if r['alpha'] == end), key=lambda r: r['certified'])
                new, s2 = run(d, family, kind, inner, seed + rr + 1, init=start['params'])
                secs += s2
                pts += new
                refinements.append(dict(bracket=[lo, hi], points=inner))
            for r in pts:
                del r['params']
            pts.sort(key=lambda r: r['alpha'])
            br = bracket(pts, kind)
            dpi_side = [r for r in pts if (r['alpha'] >= 0.5 if kind == 'sandwiched' else r['alpha'] <= 2.0)]
            out = dict(d=d, family=family, kind=kind, points=pts, bracket=br, refinements=refinements,
                       seconds=secs, n_violations=sum(r['certified'] > TOL for r in pts),
                       max_violation=max(r['certified'] for r in pts),
                       max_v_on_dpi_side=max(r['certified'] for r in dpi_side),
                       max_iterate_on_dpi_side=max(r['best_iterate'] for r in dpi_side),
                       max_cert_discrepancy=max(abs(r['value'] - r['certified']) for r in pts))
            results[key] = out
            save_partial(args, key, out)
            print(f"d={d} {family:11s} {kind:10s}: {out['n_violations']:2d}/{len(pts)} alphas violate "
                  f"(certified > {TOL:g}); threshold bracket {br}; max v {out['max_violation']:.3e}; "
                  f"max v on the DPI side of {KNOWN[kind]}: {out['max_v_on_dpi_side']:.1e} (best iterate "
                  f"{out['max_iterate_on_dpi_side']:.1e}); |double - mpmath| <= {out['max_cert_discrepancy']:.1e}; "
                  f"{secs:.0f} s", flush=True)
            print("    " + "  ".join(f"{r['alpha']:.4g}:{r['certified']:+.2e}" for r in pts), flush=True)
total = time.time() - t_start

# ----------------------------------------------------------------------------------------------------------------
# consistency with the theory (reported, not assumed)
sand = [results[f'd{d}_{f}_sandwiched'] for d in DIMS for f in FAMILIES]
petz = [results[f'd{d}_{f}_petz'] for d in DIMS for f in FAMILIES]
# combined estimates (the threshold is the same for every channel): a certified violation at alpha proves that the
# DPI fails there, so the sandwiched threshold is >= the largest certified violating alpha below 1 and the Petz
# threshold <= the smallest certified violating alpha above 1, over all families; the other end of the bracket of that
# family is the first point beyond it at which it found no violation above TOL (evidence, not proof)
brs = [r['bracket'] for r in sand if r['bracket']]
sand_lo, sand_hi = max(brs, key=lambda b: (b[0], -b[1])) if brs else (float('nan'), float('nan'))
brp = [r['bracket'] for r in petz if r['bracket']]
petz_lo, petz_hi = min(brp, key=lambda b: (b[1], -b[0])) if brp else (float('nan'), float('nan'))


def v_at(rows, a):
    """Largest certified v at alpha = a over the families (NaN if a was not searched)."""
    vs = [p['certified'] for r in rows for p in r['points'] if p['alpha'] == a]
    return max(vs) if vs else float('nan')


sand_lo_v, petz_hi_v = v_at(sand, sand_lo), v_at(petz, petz_hi)
max_dpi_side = max(r['max_v_on_dpi_side'] for r in sand + petz)
print(f"sandwiched: DPI violated (certified) at alpha = {sand_lo:.6f} (v = {sand_lo_v:.2e}), no violation > {TOL:g} "
      f"found at {sand_hi:.6f} -> threshold in [{sand_lo:.6f}, {sand_hi:.6f}] (theory 1/2); Petz: no violation found at "
      f"{petz_lo:.6f}, violated (certified) at alpha = {petz_hi:.6f} (v = {petz_hi_v:.2e}) -> threshold in "
      f"[{petz_lo:.6f}, {petz_hi:.6f}] (theory 2); largest certified v where the theorems guarantee the DPI: "
      f"{max_dpi_side:.1e}; total {total:.0f} s", flush=True)

summary = dict(sandwiched_bracket=[sand_lo, sand_hi], sandwiched_v_at_lower_end=sand_lo_v,
               petz_bracket=[petz_lo, petz_hi], petz_v_at_upper_end=petz_hi_v,
               bracket_semantics='sandwiched: certified violation at the lower end, no violation > tol found at the '
                                 'upper end; Petz: no violation > tol found at the lower end, certified violation at '
                                 'the upper end',
               max_v_where_dpi_holds=max_dpi_side, seconds=total)
save_json(args, 'renyi_threshold', dict(results=results, summary=summary, alphas=ALPHAS, dims=DIMS, starts=STARTS,
                                        steps=ROUNDS * STEPS, eps=EPS, tol=TOL, digits=DIGITS,
                                        refine_rounds=REFINE_ROUNDS, refine_points=REFINE_POINTS,
                                        init_noise=INIT_NOISE, warm_lr=WARM_LR))


def dec(x, n=6):                                     # bracket ends: the refinement resolves ~1e-5 in alpha
    return f"{x:.{n}f}"


def sci_or_nan(x):
    return sci(x) if math.isfinite(x) else 'NaN'


tex = env_macro(args, 'Renyi')
tex += (f"\\newcommand{{\\RenyiStarts}}{{{STARTS}}}\n\\newcommand{{\\RenyiSteps}}{{{ROUNDS * STEPS}}}\n"
        f"\\newcommand{{\\RenyiEps}}{{{EPS:g}}}\n\\newcommand{{\\RenyiTol}}{{{sci(TOL, 1)}}}\n"
        f"\\newcommand{{\\RenyiDigits}}{{{DIGITS}}}\n\\newcommand{{\\RenyiNumAlphas}}{{{len(ALPHAS)}}}\n"
        f"\\newcommand{{\\RenyiDims}}{{{', '.join(map(str, DIMS))}}}\n"
        f"\\newcommand{{\\RenyiSandLow}}{{{dec(sand_lo)}}}\n\\newcommand{{\\RenyiSandHigh}}{{{dec(sand_hi)}}}\n"
        f"\\newcommand{{\\RenyiPetzLow}}{{{dec(petz_lo)}}}\n\\newcommand{{\\RenyiPetzHigh}}{{{dec(petz_hi)}}}\n"
        f"% certified violation at the violating end of the brackets (RenyiSandLow, RenyiPetzHigh)\n"
        f"\\newcommand{{\\RenyiSandLowV}}{{{sci_or_nan(sand_lo_v)}}}\n\\newcommand{{\\RenyiPetzHighV}}{{{sci_or_nan(petz_hi_v)}}}\n"
        f"\\newcommand{{\\RenyiMaxDPISide}}{{{sci(max(max_dpi_side, 1e-300), 1) if max_dpi_side > 0 else '0'}}}\n"
        f"\\newcommand{{\\RenyiCertErr}}{{{sci(max(r['max_cert_discrepancy'] for r in sand + petz), 1)}}}\n"
        f"\\newcommand{{\\RenyiSeconds}}{{{total:.0f}}}\n")
d0 = DIMS[0]
tex += "\\newcommand{\\RenyiRows}{%\n"               # alpha & v sandwiched (deph, Stinespring) & v Petz (deph, Stine)
cols = [results[f'd{d0}_{f}_{k}'] for k in KINDS for f in FAMILIES]
for a in ALPHAS:
    vals = []
    for c in cols:
        v = next(r['certified'] for r in c['points'] if r['alpha'] == a)
        vals.append(f"${sci(v)}$" if v > TOL else "--")
    tex += f"{a:g} & " + " & ".join(vals) + " \\\\\n"
tex += "}\n"
write_tex(args, 'renyi_threshold', tex)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
colors = {'dephasing': 'C0', 'stinespring': 'C1'}
markers = {'dephasing': ('s', 4.5), 'stinespring': ('o', 2.5)}   # dephasing drawn larger: the curves can coincide
labels = {'dephasing': 'complete dephasing', 'stinespring': 'channels, environment $d$'}
FLOOR = TOL / 30                     # certified v <= TOL (no violation found) is drawn on this line
for ax, kind in zip(axes, KINDS):
    t = KNOWN[kind]
    lo, hi = (t, 3.05) if kind == 'sandwiched' else (0.15, t)
    ax.axvspan(lo, hi, color='0.93', lw=0)
    ax.axvline(t, color='k', ls='--', lw=0.8)
    ax.axhline(TOL, color='0.6', ls=':', lw=0.8)
    for d in DIMS:
        for f in FAMILIES:
            pts = results[f'd{d}_{f}_{kind}']['points']
            a = np.array([r['alpha'] for r in pts])
            v = np.array([r['certified'] for r in pts])
            viol = v > TOL
            ls = '-' if d == DIMS[0] else ':'
            mk, ms = markers[f]
            ax.plot(a[viol], v[viol], mk + ls, color=colors[f], ms=ms, lw=0.8, label=f"{labels[f]}, $d={d}$")
            ax.plot(a[~viol], np.full((~viol).sum(), FLOOR), mk, color=colors[f], ms=ms, mfc='none', lw=0)
    ax.set_yscale('log')
    ax.set_ylim(FLOOR / 3, 30)
    ax.set_xlim(0.15, 3.05)
    ax.set_xlabel('$\\alpha$')
    ax.set_title('sandwiched $\\widetilde D_\\alpha$' if kind == 'sandwiched' else 'Petz $D_\\alpha$', fontsize=9)
    ax.text(0.98 if kind == 'sandwiched' else 0.02, 0.4, 'DPI theorem', transform=ax.transAxes, fontsize=7,
            ha='right' if kind == 'sandwiched' else 'left', color='0.4')
axes[0].set_ylabel('largest certified violation $v(\\alpha)$')
axes[0].text(0.02, 0.06, f'open symbols: no violation ($v \\leq {sci(TOL, 1)}$)',
             transform=axes[0].transAxes, fontsize=7, color='0.3')
axes[1].legend(frameon=False, fontsize=7, loc='upper left')
fig.tight_layout()
fig.savefig(f"{args.figdir}/renyi_threshold.pdf", metadata={"CreationDate": None})
finish_partial(args)
print("written:", args.out, args.figdir)
