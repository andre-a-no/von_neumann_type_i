#!/usr/bin/env python3
"""
Section "Searching for violations" of paper/v2: random sampling versus gradient ascent.

For random positive X, Y in M_k(C) the inequality |Tr(XUY)| <= Tr(XY) fails for suitable
unitaries U; the exact supremum of z(U) = |Tr(XUY)| - Tr(XY) is ||YX||_1 - Tr(XY) > 0.
We compare, for growing k, the best z found by Haar sampling and by batched gradient ascent
(multi-start, exponential chart re-centred every 50 steps) with that supremum.

Writes generated/search.tex, generated/search.json and figures/search.pdf.
"""
import math
import time

import torch

from common import parse_args, save_json, write_tex, env_macro, sync, sci

from torch_vn_algebra import TypeIAlgebra

args = parse_args(__doc__)
dev = args.device
FULL = args.mode == 'full'
DIMS = (2, 4, 8, 16, 32, 64) if FULL else (2, 4, 8, 16)
PAIRS = 16 if FULL else 4
N_SAMPLES = 1_000_000 if FULL else 20_000
STARTS = 128 if FULL else 32
ROUNDS, STEPS, LR, DECAY = 10, 50, 0.05, 0.6   # learning rate LR * DECAY**round

results = []
for k in DIMS:
    alg = TypeIAlgebra([k], [k], precision='double', device=dev)
    pos = lambda d: 0.2 + torch.rand(1, d)
    rel_sup, rel_mc, rel_opt, err_opt, found_mc, found_opt, t_mc, t_opt = [], [], [], [], 0, 0, 0.0, 0.0
    for _ in range(PAIRS):
        X = alg.operator_from_eigenvalues(pos, force_positive=True, force_self_adjoint=True)
        Y = alg.operator_from_eigenvalues(pos, force_positive=True, force_self_adjoint=True)
        trXY = (X @ Y).trace.real
        sup = ((Y @ X).trace_norm() - trXY).item()

        def z(U):
            return (X @ U @ Y).trace.abs() - trXY

        sync(dev)
        t0 = time.time()
        best = -math.inf
        left = N_SAMPLES
        while left > 0:
            b = min(left, 50_000 if k <= 16 else 5_000)
            best = max(best, z(alg.random_unitary_operator(b)).max().item())
            left -= b
        sync(dev)
        t_mc += time.time() - t0

        t0 = time.time()
        U0 = alg.random_unitary_operator(STARTS).matrix
        for rnd in range(ROUNDS):
            G = torch.zeros(STARTS, 1, k, k, dtype=alg.hilbert.dtype, device=dev, requires_grad=True)
            opt = torch.optim.Adam([G], lr=LR * DECAY ** rnd)
            for _ in range(STEPS):
                U = alg.operator(U0) @ alg.operator(G - G.conj().transpose(-2, -1)).expm()
                loss = -z(U).sum()
                opt.zero_grad()
                loss.backward()
                opt.step()
            U0 = U.matrix.detach()
        best_opt = z(alg.operator(U0)).max().item()
        sync(dev)
        t_opt += time.time() - t0

        t = trXY.item()
        rel_sup.append(sup / t)
        rel_mc.append(best / t)
        rel_opt.append(best_opt / t)
        err_opt.append((sup - best_opt) / sup)
        found_mc += best > 0
        found_opt += best_opt > 0

    def median(x):
        return float(torch.tensor(x).median())
    r = dict(k=k, pairs=PAIRS, found_by_sampling=found_mc, found_by_ascent=found_opt,
             sup_over_trXY=median(rel_sup), sampling_over_trXY=median(rel_mc), ascent_over_trXY=median(rel_opt),
             ascent_rel_error_median=median(err_opt), ascent_rel_error_max=max(err_opt),
             seconds_sampling=t_mc / PAIRS, seconds_ascent=t_opt / PAIRS)
    results.append(r)
    print(f"k={k:3d}: violation found by sampling {found_mc}/{PAIRS}, by ascent {found_opt}/{PAIRS}; "
          f"median z/Tr(XY): sup {r['sup_over_trXY']:.2e}, sampling {r['sampling_over_trXY']:+.3f}, "
          f"ascent {r['ascent_over_trXY']:.2e} (rel. error median {r['ascent_rel_error_median']:.1e}, "
          f"max {r['ascent_rel_error_max']:.1e}); {r['seconds_sampling']:.1f}s / {r['seconds_ascent']:.1f}s")

save_json(args, 'search', dict(results=results, samples=N_SAMPLES, starts=STARTS,
                               steps=ROUNDS * STEPS, pairs=PAIRS))
tex = env_macro(args, 'Search')
tex += (f"\\newcommand{{\\SearchSamples}}{{{N_SAMPLES:,}}}\n".replace(',', '\\,')
        + f"\\newcommand{{\\SearchStarts}}{{{STARTS}}}\n\\newcommand{{\\SearchSteps}}{{{ROUNDS * STEPS}}}\n"
        + f"\\newcommand{{\\SearchPairs}}{{{PAIRS}}}\n")
tex += "\\newcommand{\\SearchRows}{%\n"
for r in results:
    tex += (f"{r['k']} & {r['found_by_sampling']}/{r['pairs']} & {r['found_by_ascent']}/{r['pairs']} "
            f"& ${sci(r['sup_over_trXY'])}$ & ${r['sampling_over_trXY']:+.2f}$ "
            f"& ${sci(r['ascent_rel_error_median'])}$ & ${sci(r['ascent_rel_error_max'])}$ "
            f"& {r['seconds_sampling']:.2g} & {r['seconds_ascent']:.2g} \\\\\n")
tex += "}\n"
write_tex(args, 'search', tex)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(4.8, 3.2))
ks = [r['k'] for r in results]
ax.axhline(0, color='k', lw=0.6)
ax.plot(ks, [r['sup_over_trXY'] for r in results], 'k-', label='exact $\\sup_U z$')
ax.plot(ks, [r['ascent_over_trXY'] for r in results], 'o', label=f'gradient ascent ({STARTS} starts)')
ax.plot(ks, [r['sampling_over_trXY'] for r in results], 's-', label=f'best of {N_SAMPLES:.0e} Haar samples')
ax.set_xscale('log', base=2)
ax.set_yscale('symlog', linthresh=1e-3)
ax.set_xlabel('$k$')
ax.set_ylabel('$z / \\mathrm{Tr}(XY)$ (median over pairs)')
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(f"{args.figdir}/search.pdf")
print("written:", args.out, args.figdir)
