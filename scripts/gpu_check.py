#!/usr/bin/env python3
"""
Short GPU check (a few minutes): is the library correct on the GPU, is it free of per-sector work and host
synchronisations, and is it faster on the GPU than on all CPU cores for large batches?

    python scripts/gpu_check.py                  # cuda:0; writes gpu_check.json next to this script's output
    python scripts/gpu_check.py --device cuda:1
    python scripts/gpu_check.py --dry-run-cpu    # tests the script itself without a GPU (CPU against CPU)

Three parts, each with a verdict (exit code 1 if part 1 or 2 fails):

1. Correctness. Inputs are drawn on the CPU and copied to the GPU; every GPU result is compared with the CPU
   result for the same input (f(A), |A| by SVD and by the Hermitian dilation, polar factors, inverse, contrast,
   channels, Lindblad channel, holomorphic calculus, normal operators, gradients), in single and double
   precision. Haar sampling on the GPU (batched Gram-Schmidt) is compared with QR on the same Ginibre matrices
   and its moments are tested (z-scores).
2. Structure (no timing noise). For every operation the host synchronisations (torch.cuda.set_sync_debug_mode)
   and the kernel launches (profiler) are counted for C = 1, 16 and 256 sectors with the same total size. They
   must not grow with C: per-sector Python work and per-sector synchronisations, the cause of the slow random
   generation in the first full run, show up here directly.
3. Timing. Median and quartiles of the time per sample on the GPU and on the CPU (all threads, interleaved runs),
   for the operations and sizes of the paper's benchmark. Hypothesis: the GPU is faster for large batches; it is
   counted as supported when the GPU's upper quartile is below the CPU's lower quartile. The batched random
   generation is also compared with the per-sector algorithm of the first full run (a reference implementation in
   this script).
"""
import argparse
import json
import platform
import statistics
import sys
import time
import warnings
from pathlib import Path

import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from torch_vn_algebra import TypeIAlgebra, states, channels, dynamics, optimize as opt, cost   # noqa: E402
from torch_vn_algebra.algebra import _orthonormal_columns                                      # noqa: E402

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--device', default='cuda')
p.add_argument('--repeats', type=int, default=5, help='timed runs per operation and device')
p.add_argument('--threads', type=int, default=None, help='CPU threads (default: all)')
p.add_argument('--out', default=str(ROOT / 'gpu_check.json'))
p.add_argument('--dry-run-cpu', action='store_true', help='CPU against CPU, to test this script without a GPU')
args = p.parse_args()
if args.dry_run_cpu:
    args.device = 'cpu'
elif not torch.cuda.is_available():
    sys.exit("no CUDA device: this check needs a GPU (or use --dry-run-cpu)")
if args.threads:
    torch.set_num_threads(args.threads)
GPU = torch.device(args.device)
CPU = torch.device('cpu')
IS_CUDA = GPU.type == 'cuda'
cost.set_limits(warn_seconds=float('inf'), probe_flops=float('inf'))      # no timing probes inside timed code
warnings.simplefilter('ignore', cost.CostWarning)
report = {'environment': {}, 'correctness': [], 'structure': [], 'timing': [], 'reference': []}


def cpu_model():
    try:
        for line in open('/proc/cpuinfo'):
            if line.startswith('model name'):
                return line.split(':', 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or platform.machine()


env = dict(torch=torch.__version__, cuda=torch.version.cuda, cpu=cpu_model(), threads=torch.get_num_threads(),
           gpu=torch.cuda.get_device_name(GPU) if IS_CUDA else 'none (dry run)')
report['environment'] = env
print(f"torch {env['torch']}, CUDA {env['cuda']}, GPU {env['gpu']}, CPU {env['cpu']} ({env['threads']} threads)\n")
failures = []


def sync(dev):
    if dev.type == 'cuda':
        torch.cuda.synchronize(dev)


# ---------------------------------------------------------------------------------------------------------------
# 1. correctness
# ---------------------------------------------------------------------------------------------------------------
print("1. correctness: GPU result against CPU result for the same input")


def compare(name, f, precision, tol):
    """f(device, generator_seed) -> tensor; the inputs are drawn on the CPU inside f with a fixed seed."""
    a = f(CPU).detach().to(torch.complex128)
    b = f(GPU).detach().cpu().to(torch.complex128)
    err = ((a - b).abs().max() / a.abs().max().clamp(min=1e-300)).item()
    ok = err <= tol
    report['correctness'].append(dict(name=name, precision=precision, rel_error=err, tol=tol, ok=ok))
    print(f"  [{'ok' if ok else 'FAIL'}] {name:38s} {precision:6s} relative difference {err:.1e} (tol {tol:.0e})")
    if not ok:
        failures.append(f"correctness: {name} ({precision})")


def make_inputs(dev, precision, k, C, B=16, seed=0, uneven=False):
    ks = [k] * C if not uneven else [max(1, k - (c % 3)) for c in range(C)]
    alg = TypeIAlgebra(ks, ks, precision=precision, device=dev)
    g = torch.Generator().manual_seed(seed)
    dt = alg.dtype
    rdt = alg.hilbert.real_dtype
    Z = torch.complex(torch.randn(B, C, k, k, generator=g, dtype=rdt), torch.randn(B, C, k, k, generator=g, dtype=rdt))
    mask = alg._const('block_mask', torch.bool).cpu()
    Z = (Z * mask).to(dt)
    return alg, Z.to(dev)


for precision, tol_f, tol_s in (('single', 2e-4, 2e-4), ('double', 1e-10, 1e-9)):
    for k in (16, 64):                          # k = 64: |A| through the Hermitian dilation on the GPU
        def herm(dev, k=k):
            alg, Z = make_inputs(dev, precision, k, 4, uneven=True)
            return alg, alg.operator(Z + Z.mH, is_self_adjoint=True)

        compare(f"f(A) = exp(-A/4), k={k}", lambda d: herm(d)[1].apply_function(lambda w: torch.exp(-0.25 * w)).matrix,
                precision, tol_f)
        compare(f"|A| (general A), k={k}", lambda d: make_inputs(d, precision, k, 4, uneven=True)[0].operator(
            make_inputs(d, precision, k, 4, uneven=True)[1]).abs().matrix, precision, tol_s)
        compare(f"trace norm, k={k}", lambda d: make_inputs(d, precision, k, 4)[0].operator(
            make_inputs(d, precision, k, 4)[1]).trace_norm(), precision, tol_s)
        compare(f"polar factor |A|, k={k}", lambda d: make_inputs(d, precision, k, 4)[0].operator(
            make_inputs(d, precision, k, 4)[1]).polar()[1].matrix, precision, tol_s)
    compare("contrast of A A^* + 1/10", lambda d: (lambda alg, Z: alg.operator(Z @ Z.mH + 0.1 * alg._const('eye'),
                                                                               is_self_adjoint=True, is_positive=True)
                                                    .michelson_contrast)(*make_inputs(d, precision, 16, 8)),
            precision, tol_f)
    compare("inverse of A A^* + 1", lambda d: (lambda alg, Z: alg.operator(Z @ Z.mH + alg._const('eye'),
                                                                           is_self_adjoint=True).inverse().matrix)(
        *make_inputs(d, precision, 16, 8)), precision, tol_f)
    compare("expm(A / 10)", lambda d: (lambda alg, Z: alg.operator(Z).expm(0.1).matrix)(
        *make_inputs(d, precision, 16, 8)), precision, tol_f)
    compare("holomorphic exp(A / 10)", lambda d: (lambda alg, Z: alg.operator(Z * 0.1).holomorphic_function(
        torch.exp).matrix)(*make_inputs(d, precision, 8, 4)), precision, 50 * tol_f)

    def unitary_power(d):
        alg, Z = make_inputs(d, precision, 12, 4)
        H = alg.operator(Z + Z.mH, is_self_adjoint=True)
        U = H.apply_function(lambda w: torch.exp(-1j * w))
        return U.apply_function(lambda z: z ** 3).matrix
    compare("normal operator: U^3", unitary_power, precision, 10 * tol_f)

    def channel_apply(d):
        alg, Z = make_inputs(d, precision, 8, 6)
        alg2, Kz = make_inputs(d, precision, 8, 6, B=3 * 16, seed=1)
        K = Kz.reshape(16, 3, 6, 8, 8)
        rho = alg.operator(Z @ Z.mH, is_self_adjoint=True)       # an operator: K is not trace preserving
        return channels.Channel(alg, K)(rho).matrix
    compare("channel (Kraus rank 3)", channel_apply, precision, tol_f)

    def lindblad(d):
        alg, Z = make_inputs(d, precision, 4, 3)
        H = alg.operator(Z + Z.mH, is_self_adjoint=True)
        L = alg.operator(0.3 * Z)
        rho = states.maximally_mixed_state(alg, 16)
        return dynamics.lindblad_channel(alg, H, [L], t=0.5)(rho).matrix
    compare("Lindblad channel exp(tL)", lindblad, precision, 1e-9 if precision == 'double' else 1e-4)

    def gradient(d):
        alg, Z = make_inputs(d, precision, 16, 4)
        X = (Z + Z.mH).requires_grad_(True)
        A = alg.operator(X, is_self_adjoint=True)
        A.apply_function(torch.sin).trace.real.sum().backward()
        return X.grad
    compare("gradient of Tr sin(A)", gradient, precision, 10 * tol_f)

# Haar sampling on the GPU path: Gram-Schmidt = QR with phase fix on the same matrices; moments
for dtype, tol in ((torch.complex64, 2e-5), (torch.complex128, 1e-12), (torch.float32, 2e-5)):
    for n in (8, 64, 256):
        Z = torch.randn(max(2 * n, 256), n, n, dtype=dtype, generator=torch.Generator().manual_seed(n))
        Q, R = torch.linalg.qr(Z.to(torch.complex128 if dtype.is_complex else torch.float64))
        dg = torch.diagonal(R, dim1=-2, dim2=-1)
        ref = Q * (dg.sgn() if dtype.is_complex else dg.sign()).unsqueeze(-2)
        got = _orthonormal_columns(Z.to(GPU)).cpu().to(ref.dtype)
        err = (got - ref).abs().max().item()
        ok = err <= tol * n ** 0.5
        report['correctness'].append(dict(name=f"Gram-Schmidt = QR, n={n}", precision=str(dtype), rel_error=err,
                                          tol=tol * n ** 0.5, ok=ok))
        print(f"  [{'ok' if ok else 'FAIL'}] Gram-Schmidt = QR with phase fix, n={n:3d} {str(dtype):15s} "
              f"max difference {err:.1e}")
        if not ok:
            failures.append(f"correctness: Gram-Schmidt n={n} {dtype}")
for cplx, n in ((True, 6), (False, 6), (True, 64)):
    alg = TypeIAlgebra([n], [n], complex_valued=cplx, precision='double', device=GPU)
    N = 200000 if n < 32 else 20000
    u2 = (alg.random_unitary(n, batch_size=N)[:, 0, 0].abs() ** 2).cpu()
    m4_th = 2 / (n * (n + 1)) if cplx else 3 / (n * (n + 2))       # E|U11|^4 for U(n) and O(n)
    zs = [((u2.mean() - 1 / n) / (u2.std() / N ** 0.5)).item(),
          (((u2 ** 2).mean() - m4_th) / ((u2 ** 2).std() / N ** 0.5)).item()]
    ok = all(abs(z) < 5 for z in zs)
    report['correctness'].append(dict(name=f"Haar moments {'U' if cplx else 'O'}({n})", z=zs, ok=ok))
    print(f"  [{'ok' if ok else 'FAIL'}] Haar moments on the GPU path, {'U' if cplx else 'O'}({n}): "
          f"z = {zs[0]:+.2f} (E|U11|^2), {zs[1]:+.2f} (E|U11|^4)")
    if not ok:
        failures.append(f"correctness: Haar moments {'U' if cplx else 'O'}({n})")


# ---------------------------------------------------------------------------------------------------------------
# operations of the benchmark
# ---------------------------------------------------------------------------------------------------------------
def operations(alg, B, dev):
    k = alg.k_max
    X = alg.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(B, d, device=dev), batch_size=B,
                                      force_positive=True).matrix
    rho = states.random_density_matrix(alg, B)
    Phi = channels.random_channel(alg, 4, B)
    Xp, Yp, Up = opt.PositiveParam(alg, B, delta=torch.full((B,), 0.5, device=dev)), opt.PositiveParam(alg, B), \
        opt.UnitaryParam(alg, B)

    def step():
        x, y = Xp.operator(), Yp.operator()
        z = (x @ Up.operator() @ y).trace.abs() - (x @ y).trace.real
        z.sum().backward()

    ops = {
        'generate': lambda: alg.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(B, d, device=dev),
                                                         batch_size=B, force_positive=True).matrix,
        'Haar unitary': lambda: alg.random_unitary_operator(B).matrix,
        'random state': lambda: states.random_density_matrix(alg, B).matrix,
        'random channel': lambda: channels.random_channel(alg, 4, B).kraus,
        'f(A) eigh': lambda: alg.operator(X, is_self_adjoint=True).apply_function(torch.sqrt).matrix,
        '|A|': lambda: alg.operator(X).abs().matrix,
        'contrast': lambda: alg.operator(X, is_self_adjoint=True, is_positive=True).michelson_contrast,
        'channel': lambda: Phi(rho).matrix,
        'opt. step': step,
    }
    if k <= 16:
        BL = max(1, min(B, 2 ** 22 // (alg.C * k ** 4)))
        H = alg.operator(X[:BL], is_self_adjoint=True)
        L = alg.operator(rho.matrix[:BL])
        ops['exp(tL)'] = lambda: dynamics.lindblad_channel(alg, H, [L], t=1.0).kraus
    return ops


# ---------------------------------------------------------------------------------------------------------------
# 2. structure: synchronisations and kernel launches do not grow with the number of sectors
# ---------------------------------------------------------------------------------------------------------------
print("\n2. structure: host synchronisations and kernel launches for C = 1, 16, 256 sectors (same total size)")
if IS_CUDA:
    from torch.profiler import profile, ProfilerActivity

    def count(fn):
        fn()
        sync(GPU)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            torch.cuda.set_sync_debug_mode('warn')
            try:
                fn()
            finally:
                torch.cuda.set_sync_debug_mode(0)
        n_sync = sum('synchroniz' in str(x.message) for x in w)
        with profile(activities=[ProfilerActivity.CPU]) as prof:
            fn()
            sync(GPU)
        n_launch = sum(e.count for e in prof.key_averages() if e.key.startswith('cudaLaunchKernel'))
        return n_sync, n_launch

    k = 8
    configs = [(1, 4096), (16, 256), (256, 16)]                    # B * C fixed
    counts = {}
    for C, B in configs:
        torch.manual_seed(0)
        alg = TypeIAlgebra([k] * C, [k] * C, device=GPU)
        for name, fn in operations(alg, B, GPU).items():
            counts.setdefault(name, {})[C] = count(fn)
    for name, byC in counts.items():
        s1, l1 = byC[1]
        s256, l256 = byC[256]
        ok = s256 <= s1 + 1 and l256 <= 3 * max(l1, 1)
        report['structure'].append(dict(name=name, syncs={C: v[0] for C, v in byC.items()},
                                        launches={C: v[1] for C, v in byC.items()}, ok=ok))
        print(f"  [{'ok' if ok else 'FAIL'}] {name:15s} synchronisations C=1/16/256: "
              f"{'/'.join(str(byC[C][0]) for C in (1, 16, 256)):10s} kernel launches: "
              f"{'/'.join(str(byC[C][1]) for C in (1, 16, 256))}")
        if not ok:
            failures.append(f"structure: {name} grows with the number of sectors")
else:
    print("  (needs CUDA; skipped in the dry run)")


# ---------------------------------------------------------------------------------------------------------------
# 3. timing
# ---------------------------------------------------------------------------------------------------------------
print(f"\n3. timing: time per sample [us], median (quartiles) of {args.repeats} interleaved runs; "
      f"speedup = CPU median / GPU median")


def time_pair(fn_gpu, fn_cpu):
    fn_gpu(); sync(GPU); fn_cpu()
    tg, tc = [], []
    for _ in range(args.repeats):
        sync(GPU)
        t0 = time.perf_counter(); fn_gpu(); sync(GPU); tg.append(time.perf_counter() - t0)
        t0 = time.perf_counter(); fn_cpu(); tc.append(time.perf_counter() - t0)
    q = lambda t: statistics.quantiles(t, n=4) if len(t) >= 2 else [t[0]] * 3
    return q(tg), q(tc)


SIZES = [(4, 256, 1024), (16, 16, 1024), (16, 256, 64), (64, 16, 64)]   # (k, C, batch) as in the benchmark
print(f"  {'operation':15s} {'k':>3s} {'C':>4s} {'batch':>5s} {'GPU':>22s} {'CPU':>22s} {'speedup':>8s}")
for k, C, B in SIZES:
    torch.manual_seed(0)
    og = operations(TypeIAlgebra([k] * C, [k] * C, device=GPU), B, GPU)
    oc = operations(TypeIAlgebra([k] * C, [k] * C, device=CPU), B, CPU)
    for name in og:
        (g1, g2, g3), (c1, c2, c3) = time_pair(og[name], oc[name])
        sp = c2 / g2
        supported = g3 < c1
        report['timing'].append(dict(name=name, k=k, C=C, batch=B, gpu_us=[1e6 * x / B for x in (g1, g2, g3)],
                                     cpu_us=[1e6 * x / B for x in (c1, c2, c3)], speedup=sp, gpu_faster=supported))
        print(f"  {name:15s} {k:3d} {C:4d} {B:5d} {1e6 * g2 / B:9.2f} ({1e6 * g1 / B:.2f}-{1e6 * g3 / B:.2f}) "
              f"{1e6 * c2 / B:9.2f} ({1e6 * c1 / B:.2f}-{1e6 * c3 / B:.2f}) {sp:7.1f}x{'' if supported else '  (not faster)'}")


# per-sector reference: the generation algorithm of the first full run (QR per sector, host checks per sector)
def generate_per_sector(alg, B, dev):
    mat = torch.zeros(B, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype, device=dev)
    for c, k_c in enumerate(alg.k_factors):
        eig = (0.1 + torch.rand(B, k_c, device=dev)).to(alg.dtype)
        bool(torch.any(eig.real < 0))                                  # the per-sector flag checks (host sync)
        Zg = torch.randn(B, k_c, k_c, dtype=alg.dtype, device=dev)
        Q, R = torch.linalg.qr(Zg)
        d = torch.diagonal(R, dim1=-2, dim2=-1)
        U = Q * d.sgn().unsqueeze(-2)
        mat[:, c, :k_c, :k_c] = (U * eig.unsqueeze(-2)) @ U.mH
    return mat


print("\n   batched generation against the per-sector algorithm of the first full run, on the GPU:")
for k, C, B in [(4, 256, 1024), (16, 16, 1024), (64, 16, 64)]:
    alg = TypeIAlgebra([k] * C, [k] * C, device=GPU)
    (n1, n2, n3), (o1, o2, o3) = time_pair(
        lambda: alg.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(B, d, device=GPU), batch_size=B,
                                              force_positive=True).matrix,
        lambda: generate_per_sector(alg, B, GPU))
    report['reference'].append(dict(k=k, C=C, batch=B, new_us=1e6 * n2 / B, old_us=1e6 * o2 / B, speedup=o2 / n2))
    print(f"  k={k:3d} C={C:4d} batch={B:5d}: batched {1e6 * n2 / B:9.2f} us/sample, per-sector "
          f"{1e6 * o2 / B:9.2f} us/sample, {o2 / n2:6.1f}x faster")

# verdict
big = [r for r in report['timing'] if r['name'] not in ('random channel',)]
n_sup = sum(r['gpu_faster'] for r in big)
print(f"\nhypothesis 'GPU faster than all CPU cores for large batches': supported in {n_sup} of {len(big)} "
      f"(operation, size) pairs (GPU upper quartile below CPU lower quartile)")
# clearly slower on the GPU: its lower quartile above the CPU's upper quartile
slow = [f"{r['name']} k={r['k']} C={r['C']}" for r in big if r['gpu_us'][0] > r['cpu_us'][2]]
if slow:
    print("  clearly slower on the GPU: " + "; ".join(slow))
report['verdict'] = dict(supported=n_sup, total=len(big), slower_on_gpu=slow, failures=failures)
with open(args.out, 'w') as f:
    json.dump(report, f, indent=1)
print(f"\nwritten: {args.out}")
if failures:
    print("FAILED:\n  " + "\n  ".join(failures))
    sys.exit(1)
print("parts 1 and 2 passed")
