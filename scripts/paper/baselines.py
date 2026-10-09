#!/usr/bin/env python3
"""
Section "Comparison with NumPy/SciPy" of paper/v2: the same computations written directly with
NumPy/SciPy and with the library.

  A. Monte Carlo over M = (+)_{c<C} M_k(R): B pairs X, Y >= 0 with prescribed Michelson contrast and
     Haar-orthogonal eigenvectors per block, a Haar-orthogonal U in M, and z = |Tr(XUY)| - Tr(XY).
       (a1) NumPy/SciPy, loop over samples and sectors (scipy.stats.ortho_group)
       (a2) NumPy, hand-vectorised (stacked QR and matmul)
       (a3) the library, one batched call (CPU, and GPU if selected)
  B. Ground states of sector Hamiltonians (Heisenberg ring and random-field chain, sector N = L/2):
       (b1) scipy.sparse.linalg.eigsh, one call per Hamiltonian
       (b2) krylov.ground_state, all disorder realisations in one batch

Writes generated/baselines.tex and generated/baselines.json.
"""
import math
import time
import warnings

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import torch
from scipy.stats import ortho_group

from common import parse_args, save_json, write_tex, env_macro, sync

from torch_vn_algebra import TypeIAlgebra, SpinChain, krylov, cost

warnings.simplefilter('ignore', cost.CostWarning)
args = parse_args(__doc__)
dev = args.device
FULL = args.mode in ('full', 'check')            # check: the sizes of full, minimal repetitions
CHECK = args.mode == 'check'
rng = np.random.default_rng(args.seed)
rec, rows_a, rows_b = {}, [], []


def spectra_np(B, C, k):
    delta = rng.random((B, C, 1))
    lo = (1 - delta) / (1 + delta)
    lam = lo + (1 - lo) * rng.random((B, C, k))
    lam[..., 0], lam[..., 1] = 1.0, lo[..., 0]
    return lam


# A. Monte Carlo ----------------------------------------------------------------------------------
def numpy_loop(B, C, k):
    z = np.empty(B)
    for b in range(B):
        lx, ly = spectra_np(1, C, k)[0], spectra_np(1, C, k)[0]
        t1 = t2 = 0.0
        for c in range(C):
            Vx, Vy, U = ortho_group.rvs(k, random_state=rng), ortho_group.rvs(k, random_state=rng), \
                ortho_group.rvs(k, random_state=rng)
            X = (Vx * lx[c]) @ Vx.T
            Y = (Vy * ly[c]) @ Vy.T
            t1 += np.trace(X @ U @ Y)
            t2 += np.trace(X @ Y)
        z[b] = abs(t1) - t2
    return z


def haar_np(shape, k):
    Z = rng.standard_normal((*shape, k, k))
    Q, R = np.linalg.qr(Z)
    return Q * np.sign(np.diagonal(R, axis1=-2, axis2=-1))[..., None, :]


def numpy_vectorised(B, C, k):
    lx, ly = spectra_np(B, C, k), spectra_np(B, C, k)
    Vx, Vy, U = haar_np((B, C), k), haar_np((B, C), k), haar_np((B, C), k)
    X = (Vx * lx[..., None, :]) @ np.swapaxes(Vx, -1, -2)
    Y = (Vy * ly[..., None, :]) @ np.swapaxes(Vy, -1, -2)
    t1 = np.trace(X @ U @ Y, axis1=-2, axis2=-1).sum(-1)
    t2 = np.trace(X @ Y, axis1=-2, axis2=-1).sum(-1)
    return np.abs(t1) - t2


def library(B, C, k, device):
    alg = TypeIAlgebra([k] * C, [k] * C, complex_valued=False, precision='double', device=device)

    def sampler(d):
        delta = torch.rand(B, 1, device=device)
        lo = (1 - delta) / (1 + delta)
        lam = lo + (1 - lo) * torch.rand(B, d, device=device)
        lam[:, 0], lam[:, 1] = 1.0, lo[:, 0]
        return lam
    X = alg.operator_from_eigenvalues(sampler, batch_size=B, force_positive=True, force_self_adjoint=True)
    Y = alg.operator_from_eigenvalues(sampler, batch_size=B, force_positive=True, force_self_adjoint=True)
    U = alg.random_unitary_operator(B)
    return ((X @ U @ Y).trace.abs() - (X @ Y).trace.real).cpu().numpy()


def timed(fn, *a):
    fn(*a)                                                  # warm-up
    t0 = time.perf_counter()
    out = fn(*a)
    if 'device' in fn.__code__.co_varnames:
        sync(a[-1])
    return time.perf_counter() - t0, out


B_LOOP = 20 if CHECK else 2000 if FULL else 300
B_MAX = 100_000 if FULL else 10_000
ENTRIES = 2 ** 24 if FULL else 2 ** 21     # entries per array; the vectorised versions hold ~10 such arrays
for C, k in ((1, 2), (16, 16), (32, 16)) if FULL else ((1, 2), (8, 16)):
    B_VEC = max(1000, min(B_MAX, ENTRIES // (C * k * k)))
    t_loop, z_loop = timed(numpy_loop, B_LOOP, C, k)
    t_vec, z_vec = timed(numpy_vectorised, B_VEC, C, k)
    t_cpu, z_cpu = timed(library, B_VEC, C, k, 'cpu')
    r = dict(C=C, k=k, batch=B_VEC, us_loop=1e6 * t_loop / B_LOOP, us_numpy_vec=1e6 * t_vec / B_VEC,
             us_lib_cpu=1e6 * t_cpu / B_VEC, mean_z=dict(loop=float(z_loop.mean()), vec=float(z_vec.mean()),
                                                         lib=float(z_cpu.mean())))
    if dev != 'cpu':
        t_gpu, _ = timed(library, B_VEC, C, k, dev)
        r['us_lib_gpu'] = 1e6 * t_gpu / B_VEC
    rows_a.append(r)
    print(f"A: C={C:2d} k={k:2d}: NumPy loop {r['us_loop']:.1f} us/sample, NumPy vectorised "
          f"{r['us_numpy_vec']:.2f}, library CPU {r['us_lib_cpu']:.2f}"
          + (f", library GPU {r['us_lib_gpu']:.3f}" if 'us_lib_gpu' in r else '')
          + f";  mean z: {r['mean_z']['loop']:.3f} / {r['mean_z']['vec']:.3f} / {r['mean_z']['lib']:.3f}")


# B. ground states ---------------------------------------------------------------------------------
def scipy_matrix(H, b):
    A = H.offdiag.coalesce()
    idx = A.indices().cpu().numpy()
    off = sp.csr_matrix((A.values().real.cpu().numpy(), (idx[0], idx[1])), shape=(H.dim, H.dim))
    return off + sp.diags(H.diag[b].cpu().numpy())


for L, R in (((16, 1), (16, 2 if CHECK else 20), (18, 1), (20, 1)) if FULL else ((14, 1), (14, 10), (16, 1))):
    ch = SpinChain(L, 'periodic', complex_valued=False, device=dev)
    W = 0.0 if R == 1 else 3.0
    h = (2 * torch.rand(R, L, dtype=torch.float64) - 1) * W
    H = ch.xxz_sparse(1.0, 1.0, h, sector=L // 2)
    mats = [scipy_matrix(H, b) for b in range(R)]
    t0 = time.perf_counter()
    e_sp = [spla.eigsh(M, k=1, which='SA', tol=1e-12)[0][0] for M in mats]
    t_sp = time.perf_counter() - t0
    krylov.ground_state(H)                                   # warm-up
    sync(dev)
    t0 = time.perf_counter()
    E, _, _ = krylov.ground_state(H)
    sync(dev)
    t_kr = time.perf_counter() - t0
    diff = float(np.max(np.abs(E.cpu().numpy() - np.array(e_sp))))
    rows_b.append(dict(L=L, realisations=R, dim=H.dim, s_scipy=t_sp, s_library=t_kr, max_diff=diff))
    print(f"B: L={L} ({R} Hamiltonian{'s' if R > 1 else ''}, dim {H.dim}): eigsh {t_sp:.2f} s, "
          f"library {t_kr:.2f} s ({dev}), |dE| <= {diff:.1e}")

save_json(args, 'baselines', dict(monte_carlo=rows_a, ground_states=rows_b, b_loop=B_LOOP))
tex = env_macro(args, 'Base')
gpu = dev != 'cpu'
tex += f"\\newcommand{{\\BaseHasGPU}}{{{'1' if gpu else '0'}}}\n"
tex += "\\newcommand{\\BaseMCRows}{%\n"
for r in rows_a:
    tex += (f"{r['C']} & {r['k']} & {r['us_loop']:.3g} & {r['us_numpy_vec']:.3g} & {r['us_lib_cpu']:.3g}"
            + (f" & {r['us_lib_gpu']:.3g}" if gpu else " & --") + " \\\\\n")
tex += "}\n\\newcommand{\\BaseGSRows}{%\n"
for r in rows_b:
    tex += (f"{r['L']} & {r['dim']} & {r['realisations']} & {r['s_scipy']:.3g} & {r['s_library']:.3g} "
            f"& {r['max_diff']:.1e} \\\\\n")
tex += "}\n"
write_tex(args, 'baselines', tex)
print("written:", args.out)
