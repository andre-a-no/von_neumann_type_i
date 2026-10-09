#!/usr/bin/env python3
"""
Section "Performance" of paper: wall-clock time of library operations (not raw kernels),
per Monte Carlo sample, on the selected device and on the CPU.

Operations (each on a batch of B operators in M = (+)_{c<C} M_k(C), single precision):
  generate     operator_from_eigenvalues(...) incl. Haar unitaries, materialised
  f(A) eigh    apply_function(sqrt) (spectral theorem, torch.linalg.eigh)
  |A| SVD      abs() (batched SVD)
  contrast     michelson_contrast (eigvalsh)
  channel      random Stinespring channel of Kraus rank 4 applied to a state
  exp(tL)      Lindblad channel (k^2 x k^2 matrix exponential), k <= 16
  opt. step    one forward/backward step of the Exp. 1 objective in optimize.extremize

The batch size B is chosen so that B * C * k^2 is about 2^22 (2^18 in quick mode).
Writes generated/benchmark.tex, generated/benchmark.json and figures/benchmark.pdf.

    python scripts/benchmark_library.py --mode full --device cuda
"""
import math
import statistics
import time

import torch

from common import parse_args, save_json, write_tex, env_macro, sync, num3

from torch_vn_algebra import TypeIAlgebra, channels, cost, dynamics, optimize as opt, states

args = parse_args(__doc__)
FULL = args.mode in ('full', 'check')            # check: the sizes of full, minimal repetitions
CHECK = args.mode == 'check'
DEVICES = [args.device] + (['cpu'] if args.device != 'cpu' else [])
DIMS = (4, 16, 64) if FULL else (4, 16)
CHANNELS = (1, 16, 256) if FULL else (1, 16)
TARGET = 2 ** 22 if FULL else 2 ** 18
REPEATS = 1 if CHECK else 5 if FULL else 3


def timeit(fn, device):
    fn()                               # warm-up
    sync(device)
    times = []
    for _ in range(REPEATS):
        t0 = time.perf_counter()
        fn()
        sync(device)
        times.append(time.perf_counter() - t0)
    return statistics.median(times)


def operations(alg, B, dev, k, C):
    def gen():
        return alg.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(B, d, device=dev), batch_size=B,
                                             force_positive=True, force_self_adjoint=True)
    X = gen()
    X.matrix
    rho = states.random_density_matrix(alg, B)
    Phi = channels.random_channel(alg, 4, B)
    ops = {
        'generate': lambda: gen().matrix,
        'f(A) eigh': lambda: alg.operator(X.matrix, is_self_adjoint=True).apply_function(torch.sqrt).matrix,
        '|A| SVD': lambda: alg.operator(X.matrix).abs().matrix,
        'contrast': lambda: alg.operator(X.matrix, is_self_adjoint=True, is_positive=True).michelson_contrast,
        'channel': lambda: Phi(rho).matrix,
    }
    if k <= 16:
        # k^2 x k^2 superoperators: use a smaller batch (B_L * C * k^4 ~ TARGET elements)
        BL = max(1, min(B, TARGET // (C * k ** 4)))
        H = alg.operator(X.matrix[:BL], is_self_adjoint=True)
        L = alg.operator(rho.matrix[:BL])
        ops['exp(tL)'] = (lambda: dynamics.lindblad_channel(alg, H, [L], t=1.0).kraus, BL)

    Xp, Yp, Up = opt.PositiveParam(alg, B, delta=torch.full((B,), 0.5, device=dev)), opt.PositiveParam(alg, B), \
        opt.UnitaryParam(alg, B)

    def step():
        x, y = Xp.operator(), Yp.operator()
        z = (x @ Up.operator() @ y).trace.abs() - (x @ y).trace.real
        z.sum().backward()
    ops['opt. step'] = step
    return ops


rows = []
cost.set_limits(warn_seconds=float('inf'), probe_flops=float('inf'))   # no timing probes; memory checks stay on
if True:
    for k in DIMS:
        for C in CHANNELS:
            B = max(1, min(4096, TARGET // (C * k * k)))
            per_dev, batch_of = {}, {}
            for dev in DEVICES:
                torch.manual_seed(args.seed)
                alg = TypeIAlgebra([k] * C, [k] * C, complex_valued=True, device=dev)
                for name, fn in operations(alg, B, dev, k, C).items():
                    fn, b = fn if isinstance(fn, tuple) else (fn, B)
                    per_dev.setdefault(name, {})[dev] = timeit(fn, dev) / b
                    batch_of[name] = b                    # exp(tL) runs with a smaller batch
            for name, t in per_dev.items():
                row = dict(op=name, k=k, C=C, batch=batch_of[name], **{f't_{d}': t[d] for d in DEVICES})
                if len(DEVICES) == 2:
                    row['speedup'] = t['cpu'] / t[args.device]
                rows.append(row)
                print(f"k={k:3d} C={C:4d} B={B:5d} {name:10s} " +
                      "  ".join(f"{d}: {1e6 * t[d]:9.2f} us/sample" for d in DEVICES) +
                      (f"  speedup {row['speedup']:.1f}x" if 'speedup' in row else ''))

save_json(args, 'benchmark', dict(rows=rows, devices=DEVICES, target_elements=TARGET))
tex = env_macro(args, 'Bench')
tex += f"\\newcommand{{\\BenchCpuThreads}}{{{torch.get_num_threads()}}}\n"
tex += f"\\newcommand{{\\BenchTargetLog}}{{{int(round(math.log2(TARGET)))}}}\n"     # batch holds ~2^this entries
two = len(DEVICES) == 2
tex += f"\\newcommand{{\\BenchHasGPU}}{{{'1' if two else '0'}}}\n"
tex += "\\newcommand{\\BenchmarkRows}{%\n"
for r in rows:
    main = 1e6 * r[f't_{args.device}']
    tex += f"{r['op']} & {r['k']} & {r['C']} & {r['batch']} & {num3(main)}"
    if two:
        tex += f" & {num3(1e6 * r['t_cpu'])} & {r['speedup']:.1f}"
    tex += " \\\\\n"
tex += "}\n"
write_tex(args, 'benchmark', tex)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ops = list(dict.fromkeys(r['op'] for r in rows))
configs = list(dict.fromkeys((r['k'], r['C']) for r in rows))
fig, ax = plt.subplots(figsize=(7.5, 3.4))
width = 0.8 / len(configs)
for j, (k, C) in enumerate(configs):
    vals = []
    for op in ops:
        r = next((r for r in rows if r['op'] == op and r['k'] == k and r['C'] == C), None)
        vals.append(float('nan') if r is None else (r['speedup'] if two else 1e6 * r[f't_{args.device}']))
    ax.bar([i + j * width for i in range(len(ops))], vals, width, label=f'k={k}, C={C}')
ax.set_xticks([i + 0.4 - width / 2 for i in range(len(ops))])
ax.set_xticklabels(ops, rotation=20, fontsize=8)
ax.set_yscale('log')
ax.set_ylabel('speedup CPU / GPU' if two else 'time per sample [us]')
if two:
    ax.axhline(1, color='k', lw=0.6)
ax.legend(frameon=False, fontsize=7, ncol=2)
fig.tight_layout()
fig.savefig(f"{args.figdir}/benchmark.pdf", metadata={"CreationDate": None})
print("written:", args.out, args.figdir)
