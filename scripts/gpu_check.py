#!/usr/bin/env python3
"""
Quick check of a GPU (about a minute): correctness of the GPU code paths and CPU/GPU timings of the
operations of the paper's benchmark, at a few sizes.

    python scripts/gpu_check.py               # cuda:0
    python scripts/gpu_check.py --device cuda:1

Checks (exit code 1 if one fails):
  * Haar sampling on the GPU (batched Gram-Schmidt) gives the same unitaries as QR with the phase
    correction, and the results are unitary;
  * random operators, states and channels have the right properties on the GPU;
  * the GPU is not slower than the CPU where it should be faster: batched generation, eigh, |A|,
    channels and the optimisation step for a large batch (a GPU that is slower there points to a
    per-matrix loop somewhere, as in the generation code before it was batched).
Then a table of times per sample is printed. Nothing is written to the paper.
"""
import argparse
import sys
import time
from pathlib import Path

import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import warnings                                                                      # noqa: E402
from torch_vn_algebra import TypeIAlgebra, states, channels, optimize as opt        # noqa: E402
from torch_vn_algebra.algebra import _orthonormal_columns                           # noqa: E402
from torch_vn_algebra.cost import CostWarning                                       # noqa: E402

warnings.simplefilter('ignore', CostWarning)

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--device', default='cuda')
p.add_argument('--repeats', type=int, default=3)
p.add_argument('--dry-run-cpu', action='store_true', help=argparse.SUPPRESS)    # test this script without a GPU
args = p.parse_args()
if args.dry_run_cpu:
    args.device = 'cpu'
elif not torch.cuda.is_available():
    sys.exit("no CUDA device: this check needs a GPU")
gpu = torch.device(args.device)
print(f"torch {torch.__version__}, CUDA {torch.version.cuda}, "
      f"{torch.cuda.get_device_name(gpu) if gpu.type == 'cuda' else 'no GPU (dry run)'}, "
      f"CPU threads {torch.get_num_threads()}")
failures = []


def check(name, ok, detail=''):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}{'' if not detail else ': ' + detail}")
    if not ok:
        failures.append(name)


# ---- correctness of the GPU paths ----------------------------------------------------------------
print("correctness")
torch.manual_seed(0)
for dtype in (torch.float64, torch.complex128):
    Z = torch.randn(512, 12, 12, dtype=dtype)
    Q, R = torch.linalg.qr(Z)
    d = torch.diagonal(R, dim1=-2, dim2=-1)
    d = d.sgn() if dtype.is_complex else d.sign()
    err = (_orthonormal_columns(Z.to(gpu)).cpu() - Q * d.unsqueeze(-2)).abs().max().item()
    check(f"Gram-Schmidt on GPU = QR with phase fix ({dtype})", err < 1e-10, f"max difference {err:.1e}")

for precision, tol in (('single', 1e-4), ('double', 1e-12)):
    alg = TypeIAlgebra([8] * 16, [8] * 16, precision=precision, device=gpu)
    U = alg.random_unitary_operator(256).matrix
    eye = torch.eye(8, dtype=U.dtype, device=gpu)
    err = (U @ U.mH - eye).abs().max().item()
    check(f"random unitaries are unitary ({precision})", err < tol, f"{err:.1e}")
    A = alg.operator_from_eigenvalues(lambda k: 0.1 + torch.rand(256, k, device=gpu), batch_size=256,
                                      force_positive=True)
    w = A.eigenvalues()
    check(f"generated operators are positive with the requested spectrum ({precision})",
          bool((torch.cat([x.reshape(-1) for x in w]) > 0.1 - 10 * tol).all()))
    rho = states.random_density_matrix(alg, 256)
    tr = rho.Tr_blunt().real
    check(f"random states have trace 1 ({precision})", (tr - 1).abs().max().item() < 10 * tol)
    Phi = channels.random_channel(alg, 3, 64)
    tp = (Phi.kraus.mH @ Phi.kraus).sum(1) - eye
    check(f"random channels are trace preserving ({precision})", tp.abs().max().item() < 10 * tol)

# Haar moments on the GPU path: E|U_11|^2 = 1/n, E|U_11|^4 = 2/(n(n+1)) for U(n)
alg = TypeIAlgebra([6], [6], precision='double', device=gpu)
u = alg.random_unitary(6, batch_size=200000)[:, 0, 0].abs()
m2, m4 = (u ** 2).mean().item(), (u ** 4).mean().item()
se4 = (u ** 4).std().item() / 200000 ** 0.5
check("Haar moments on the GPU path", abs(m2 - 1 / 6) < 5e-3 and abs(m4 - 2 / 42) < 5 * se4,
      f"E|U11|^2 = {m2:.4f} (1/6 = {1 / 6:.4f}), E|U11|^4 = {m4:.5f} ({2 / 42:.5f})")


# ---- timings -------------------------------------------------------------------------------------
def timed(fn, device):
    fn()
    if device.type == 'cuda':
        torch.cuda.synchronize(device)
    best = float('inf')
    for _ in range(args.repeats):
        t0 = time.perf_counter()
        fn()
        if device.type == 'cuda':
            torch.cuda.synchronize(device)
        best = min(best, time.perf_counter() - t0)
    return best


def operations(k, C, B, device):
    alg = TypeIAlgebra([k] * C, [k] * C, device=device)
    X = alg.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(B, d, device=device), batch_size=B,
                                      force_positive=True).matrix
    rho = states.random_density_matrix(alg, B)
    Phi = channels.random_channel(alg, 4, B)
    Xp, Yp, Up = opt.PositiveParam(alg, B), opt.PositiveParam(alg, B), opt.UnitaryParam(alg, B)

    def step():
        x, y = Xp.operator(), Yp.operator()
        z = (x @ Up.operator() @ y).trace.abs() - (x @ y).trace.real
        z.sum().backward()

    return {
        'generate': lambda: alg.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(B, d, device=device),
                                                         batch_size=B, force_positive=True).matrix,
        'Haar unitary': lambda: alg.random_unitary_operator(B).matrix,
        'random state': lambda: states.random_density_matrix(alg, B).matrix,
        'f(A) eigh': lambda: alg.operator(X, is_self_adjoint=True).apply_function(torch.sqrt).matrix,
        '|A|': lambda: alg.operator(X).abs().matrix,
        'channel': lambda: Phi(rho).matrix,
        'opt. step': step,
    }


SHOULD_GAIN = {'generate', 'Haar unitary', 'f(A) eigh', '|A|', 'opt. step'}
print("\ntime per sample [us]; speedup = CPU / GPU")
print(f"{'operation':14s} {'k':>3s} {'C':>4s} {'batch':>6s} {'GPU':>10s} {'CPU':>10s} {'speedup':>8s}")
for k, C, B in ((4, 16, 4096), (16, 16, 1024), (64, 4, 256)):
    torch.manual_seed(0)
    ops_gpu, ops_cpu = operations(k, C, B, gpu), operations(k, C, B, torch.device('cpu'))
    for name in ops_gpu:
        tg, tc = timed(ops_gpu[name], gpu) / B, timed(ops_cpu[name], torch.device('cpu')) / B
        print(f"{name:14s} {k:3d} {C:4d} {B:6d} {1e6 * tg:10.2f} {1e6 * tc:10.2f} {tc / tg:8.1f}")
        if name in SHOULD_GAIN and k * k * C * B >= 2 ** 18:
            if tc / tg < 1:
                failures.append(f"{name} k={k} C={C}: GPU slower than CPU ({tc / tg:.2f}x)")

print()
if failures:
    print("FAILED:")
    for f in failures:
        print("  " + f)
    sys.exit(1)
print("all checks passed")
