"""
Padded versus packed storage (TypeIAlgebra(layout=...)) on the full algebra (+)_N M_binom(L,N) of a spin chain: memory
of one operator, peak memory of the process and run times of typical operations, each configuration in a fresh
process (so that the peak resident memory belongs to it alone).

    python scripts/layout_bench.py [--L 10 12 13] [--device cpu] [--out .quick_output/layout_bench.json]
    python scripts/layout_bench.py --calibrate [--device cuda]

--calibrate measures, on the given device, where packing starts to pay off for algebras with several block sizes
and prints the constant entries_per_kernel of the layout choice (TypeIAlgebra(layout='auto'), cost.set_limits):
for every algebra the batch B* at which packed becomes faster gives B* W / (G - 1), W the padding entries per sample
and G the number of block sizes; the median over the algebras is suggested.
"""
import argparse
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run_case(L: int, layout: str, device: str) -> dict:
    os.environ['TORCH_VN_ALGEBRA_LAYOUT'] = layout
    sys.path.insert(0, str(ROOT))
    import torch
    from torch_vn_algebra import SpinChain, dynamics, states

    def timed(fn, sync=None):
        t0 = time.perf_counter()
        out = fn()
        if device.startswith('cuda'):
            torch.cuda.synchronize()
        return out, time.perf_counter() - t0

    ch = SpinChain(L, 'periodic', device=device)
    alg = ch.algebra
    res = {'L': L, 'layout': alg.layout, 'C': alg.C, 'k_max': alg.k_max,
           'bytes_per_operator': alg.bytes_per_operator(1)}
    H, res['t_build'] = timed(lambda: ch.xxz(1.0, 0.5) * 1.0 + ch.sz(0) * 0.1)
    _, res['t_build'] = timed(lambda: H._canonical())
    H2, res['t_product'] = timed(lambda: (H @ H)._canonical())
    rho, res['t_gibbs'] = timed(lambda: states.gibbs_state(H, 0.5)._canonical())
    rho = states.gibbs_state(H, 0.5)
    _, res['t_entropy'] = timed(lambda: states.von_neumann_entropy(rho))
    jumps = [ch.sector_jump(i, 'lower') for i in range(L)]
    _, res['t_lindblad_kraus_step'] = timed(
        lambda: dynamics.lindblad_evolve(rho, H, jumps, [0.0, 0.05], rates=[0.1] * L, substeps=1,
                                         method='kraus')[-1]._canonical())
    res['peak_rss_MB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    return res


def calibrate(device: str) -> float:
    import math
    import statistics
    sys.path.insert(0, str(ROOT))
    import torch
    from torch_vn_algebra import TypeIAlgebra

    def bench(fn, n=30):
        for _ in range(3):
            fn()
        ts = []
        for _ in range(n):
            t0 = time.perf_counter()
            fn()
            if device.startswith('cuda'):
                torch.cuda.synchronize()
            ts.append(time.perf_counter() - t0)
        return statistics.median(ts)

    configs = {'1..20': list(range(1, 21)), '2..30 step 4': list(range(2, 31, 4)), '[30,5,10,15,20]': [30, 5, 10, 15, 20],
               'chain L=6': [math.comb(6, n) for n in range(7)], 'chain L=7': [math.comb(7, n) for n in range(8)],
               'chain L=8': [math.comb(8, n) for n in range(9)]}
    batches = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
    kappas = []
    for name, ks in configs.items():
        C, km = len(ks), max(ks)
        W, G = C * km * km - sum(k * k for k in ks), len(set(ks))
        prev, found = None, None
        for B in batches:
            r = {}
            for layout in ('padded', 'packed'):
                alg = TypeIAlgebra(ks, ks, precision='double', layout=layout, device=device)
                A = alg.operator(torch.randn(B, C, km, km, dtype=alg.dtype, device=device))
                U = alg.operator(torch.randn(B, C, km, km, dtype=alg.dtype, device=device))
                A._canonical(), U._canonical()
                r[layout] = bench(lambda: ((A + U * 2.0) @ U).Tr_norm())
            if r['padded'] > r['packed']:
                found = B
                break
            prev = B
        if found is None:
            print(f"{name:16s} padded faster up to batch {batches[-1]}")
            continue
        lo = (prev or found / 2) * W / (G - 1)
        k = math.sqrt(lo * found * W / (G - 1))
        kappas.append(k)
        print(f"{name:16s} W={W:7d} G={G:2d}: packed faster from batch {found}  -> entries_per_kernel ~ {k:.2e}")
    kappa = statistics.median(kappas) if kappas else float('nan')
    key = 'layout_entries_per_kernel_' + ('cuda' if device.startswith('cuda') else 'cpu')
    print(f"\nsuggested: torch_vn_algebra.cost.set_limits({key}={kappa:.2g})")
    return kappa


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, nargs='+', default=[10, 12, 13])
    ap.add_argument('--device', default='cpu')
    ap.add_argument('--out', default=str(ROOT / '.quick_output' / 'layout_bench.json'))
    ap.add_argument('--case', nargs=2, help=argparse.SUPPRESS)
    ap.add_argument('--calibrate', action='store_true')
    a = ap.parse_args()
    if a.calibrate:
        calibrate(a.device)
        return
    if a.case:
        print(json.dumps(run_case(int(a.case[0]), a.case[1], a.device)))
        return
    rows = []
    for L in a.L:
        for layout in ('padded', 'packed'):
            out = subprocess.run([sys.executable, __file__, '--case', str(L), layout, '--device', a.device],
                                 capture_output=True, text=True)
            if out.returncode:
                print(out.stderr, file=sys.stderr)
                continue
            rows.append(json.loads(out.stdout.strip().splitlines()[-1]))
            r = rows[-1]
            print(f"L={L:2d} {layout:6s}  operator {r['bytes_per_operator'] / 2**20:8.1f} MB  peak RSS "
                  f"{r['peak_rss_MB']:8.0f} MB  build {r['t_build']:.3f}s  H@H {r['t_product']:.3f}s  "
                  f"gibbs {r['t_gibbs']:.3f}s  entropy {r['t_entropy']:.3f}s  "
                  f"lindblad step {r['t_lindblad_kraus_step']:.3f}s", flush=True)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(rows, indent=1))


if __name__ == '__main__':
    main()
