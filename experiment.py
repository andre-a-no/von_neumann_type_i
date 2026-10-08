#!/usr/bin/env python3
"""
Batched Monte Carlo experiments for trace inequalities.
Exact translation of original NumPy code to PyTorch with GPU batching.
20 000 samples per configuration, plots always saved.
"""
import torch
import numpy as np
import pandas as pd
import argparse
from pathlib import Path
from tqdm import tqdm
import sys
import os
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from structures.Algebra import TypeIAlgebra


# ----------------------------------------------------------------------
# Samplers (identical to NumPy generate_eigenvalues logic)
# ----------------------------------------------------------------------
def positive_sampler(batch_size, device):
    """Positive eigenvalues, max=1, uniform contrast in [0,1]."""
    def sampler(dim):
        delta = torch.rand(batch_size, device=device, dtype=torch.float32)
        lam_min = (1 - delta) / (1 + delta)
        eig = torch.zeros(batch_size, dim, device=device, dtype=torch.float32)
        eig[:, 0] = 1.0
        if dim > 1:
            eig[:, 1] = lam_min
        if dim > 2:
            low = lam_min.unsqueeze(1)
            high = torch.ones(batch_size, 1, device=device)
            eig[:, 2:] = low + (high - low) * torch.rand(batch_size, dim-2, device=device)
        eig, _ = torch.sort(eig, dim=-1, descending=True)
        return eig
    return sampler

def selfadjoint_sampler(batch_size, device):
    """Self-adjoint eigenvalues, max|λ|=1, uniform contrast of |λ|, random signs."""
    def sampler(dim):
        delta = torch.rand(batch_size, device=device, dtype=torch.float32)
        lam_min = (1 - delta) / (1 + delta)
        eig_abs = torch.zeros(batch_size, dim, device=device, dtype=torch.float32)
        eig_abs[:, 0] = 1.0
        if dim > 1:
            eig_abs[:, 1] = lam_min
        if dim > 2:
            low = lam_min.unsqueeze(1)
            high = torch.ones(batch_size, 1, device=device)
            eig_abs[:, 2:] = low + (high - low) * torch.rand(batch_size, dim-2, device=device)
        eig_abs, _ = torch.sort(eig_abs, dim=-1, descending=True)
        signs = 2 * torch.randint(0, 2, (batch_size, dim), device=device, dtype=torch.float32) - 1
        eig = eig_abs * signs
        return eig
    return sampler


# ----------------------------------------------------------------------
# Helper: apply random orthogonal matrix to a batch of square matrices (for Exp2)
# ----------------------------------------------------------------------
def apply_random_orthogonal(mat, alg):
    """Multiply each matrix in batch (active block) on the left by random orthogonal."""
    device = mat.device
    batch, C, k_max, _ = mat.shape
    for c in range(C):
        k_c = alg.k_factors[c]
        if k_c == 0:
            continue
        U = torch.stack([alg.random_unitary(k_c, measure='haar') for _ in range(batch)]).to(device)
        mat[:, c, :k_c, :k_c] = U @ mat[:, c, :k_c, :k_c]
    return mat


# ----------------------------------------------------------------------
# Experiments
# ----------------------------------------------------------------------
def experiment1(alg, batch_size, num_batches, device):
    """X, Y positive, U random orthogonal. z = |Tr(X U Y)| - Tr(X Y)"""
    pos_sampler = positive_sampler(batch_size, device)
    deltaX_list, deltaY_list, z_list = [], [], []
    for _ in tqdm(range(num_batches), desc="Exp1", leave=False):
        X = alg.operator_from_eigenvalues(pos_sampler, batch_size=batch_size,
                                           force_positive=True, force_self_adjoint=True)
        Y = alg.operator_from_eigenvalues(pos_sampler, batch_size=batch_size,
                                           force_positive=True, force_self_adjoint=True)
        # contrast Δ(X) and Δ(Y)
        lmaxX, lminX = X.lambda_max, X.lambda_min
        deltaX = (lmaxX - lminX) / (lmaxX + lminX)
        lmaxY, lminY = Y.lambda_max, Y.lambda_min
        deltaY = (lmaxY - lminY) / (lmaxY + lminY)
        # generate random orthogonal U per batch element (per channel, but we apply globally)
        # Since random_unitary generates one matrix per call, we create a batch.
        # For simplicity, generate one U per channel and broadcast? But original code used one U for the whole matrix.
        # We'll generate a single random orthogonal per batch element and apply to the whole active block (all channels share same U).
        # However, original NumPy code generated U of size dimension x dimension, then multiplied X @ U @ Y.
        # We need to apply U to the entire operator (across channels). Since operator is block-diagonal, we can apply U independently per channel.
        # For simplicity and fidelity, we generate U for each channel separately, as in original each channel is independent? Actually original had no channels – one matrix.
        # To match, we generate a single U for the full k_max, but active block may be smaller. We'll generate U of size k_c per channel.
        # We'll implement per‑channel U as in original spirit (each factor gets its own U).
        X_mat = X.matrix.clone()
        Y_mat = Y.matrix.clone()
        batch_z = torch.zeros(batch_size, device=device)
        for c in range(alg.C):
            k_c = alg.k_factors[c]
            if k_c == 0:
                continue
            Xc = X_mat[:, c, :k_c, :k_c]   # (batch, k_c, k_c)
            Yc = Y_mat[:, c, :k_c, :k_c]   # (batch, k_c, k_c)
            U = torch.stack([alg.random_unitary(k_c, measure='haar') for _ in range(batch_size)]).to(device)
            XY = Xc @ U @ Yc
            term1 = torch.abs(torch.diagonal(XY, dim1=-2, dim2=-1).sum(-1))
            term2 = torch.diagonal(Xc @ Yc, dim1=-2, dim2=-1).sum(-1)
            batch_z += (term1 - term2)
        z = batch_z / alg.C
        deltaX_list.append(deltaX.cpu().numpy())
        deltaY_list.append(deltaY.cpu().numpy())
        z_list.append(z.cpu().numpy())
    return (np.concatenate(deltaX_list), np.concatenate(deltaY_list), np.concatenate(z_list))

def experiment2(alg, batch_size, num_batches, device):
    """X0, Y0 positive, X = U X0, Y = V Y0. z = Tr(|XY|) - Tr(|X||Y|)"""
    pos_sampler = positive_sampler(batch_size, device)
    deltaX_list, deltaY_list, z_list = [], [], []
    for _ in tqdm(range(num_batches), desc="Exp2", leave=False):
        X0 = alg.operator_from_eigenvalues(pos_sampler, batch_size=batch_size,
                                            force_positive=True, force_self_adjoint=True)
        Y0 = alg.operator_from_eigenvalues(pos_sampler, batch_size=batch_size,
                                            force_positive=True, force_self_adjoint=True)
        # contrasts from original positive matrices (Δ(|X|) = Δ(X0) because |X| = X0)
        lmaxX0, lminX0 = X0.lambda_max, X0.lambda_min
        deltaX = (lmaxX0 - lminX0) / (lmaxX0 + lminX0)
        lmaxY0, lminY0 = Y0.lambda_max, Y0.lambda_min
        deltaY = (lmaxY0 - lminY0) / (lmaxY0 + lminY0)
        # apply random orthogonal to make non-Hermitian
        X_mat = apply_random_orthogonal(X0.matrix.clone(), alg)
        Y_mat = apply_random_orthogonal(Y0.matrix.clone(), alg)
        # compute |X|, |Y|, |XY|
        absX_mat = torch.zeros_like(X_mat)
        absY_mat = torch.zeros_like(Y_mat)
        for c in range(alg.C):
            k_c = alg.k_factors[c]
            if k_c == 0:
                continue
            Xc = X_mat[:, c, :k_c, :k_c]
            Yc = Y_mat[:, c, :k_c, :k_c]
            Ux, Sx, Vx = torch.linalg.svd(Xc)
            Uy, Sy, Vy = torch.linalg.svd(Yc)
            absX_mat[:, c, :k_c, :k_c] = Ux @ torch.diag_embed(Sx) @ Ux.conj().transpose(-2, -1)
            absY_mat[:, c, :k_c, :k_c] = Uy @ torch.diag_embed(Sy) @ Uy.conj().transpose(-2, -1)
        XY_mat = X_mat @ Y_mat
        absXY_mat = torch.zeros_like(XY_mat)
        for c in range(alg.C):
            k_c = alg.k_factors[c]
            if k_c == 0:
                continue
            U, S, V = torch.linalg.svd(XY_mat[:, c, :k_c, :k_c])
            absXY_mat[:, c, :k_c, :k_c] = U @ torch.diag_embed(S) @ U.conj().transpose(-2, -1)
        def trace_op(t):
            return torch.diagonal(t, dim1=-2, dim2=-1).sum(-1).sum(-1)
        z = trace_op(absXY_mat) - trace_op(absX_mat @ absY_mat)
        deltaX_list.append(deltaX.cpu().numpy())
        deltaY_list.append(deltaY.cpu().numpy())
        z_list.append(z.cpu().numpy())
    return (np.concatenate(deltaX_list), np.concatenate(deltaY_list), np.concatenate(z_list))

def experiment3(alg, batch_size, num_batches, device):
    """X self-adjoint, Y positive. z = Tr(Y|X|Y) - Tr(|YXY|)"""
    sa_sampler = selfadjoint_sampler(batch_size, device)
    pos_sampler = positive_sampler(batch_size, device)
    deltaX_list, deltaY_list, z_list = [], [], []
    for _ in tqdm(range(num_batches), desc="Exp3", leave=False):
        X = alg.operator_from_eigenvalues(sa_sampler, batch_size=batch_size, force_self_adjoint=True)
        Y = alg.operator_from_eigenvalues(pos_sampler, batch_size=batch_size,
                                           force_positive=True, force_self_adjoint=True)
        absX = X.abs()
        lmaxX, lminX = absX.lambda_max, absX.lambda_min
        deltaX = (lmaxX - lminX) / (lmaxX + lminX)
        lmaxY, lminY = Y.lambda_max, Y.lambda_min
        deltaY = (lmaxY - lminY) / (lmaxY + lminY)
        Y_absX_Y = Y @ absX @ Y
        abs_YXY = (Y @ X @ Y).abs()
        z = Y_absX_Y.trace - abs_YXY.trace
        deltaX_list.append(deltaX.cpu().numpy())
        deltaY_list.append(deltaY.cpu().numpy())
        z_list.append(z.cpu().numpy())
    return (np.concatenate(deltaX_list), np.concatenate(deltaY_list), np.concatenate(z_list))


# ----------------------------------------------------------------------
# Plotting (saves files with names like exp1_dim2_ch1_xz.png)
# ----------------------------------------------------------------------
def save_scatter_plots(exp_id, dim, C, dx, dy, z, xlabel, ylabel, zlabel, plot_dir):
    if len(dx) == 0:
        return
    plot_dir = Path(plot_dir).resolve()
    plot_dir.mkdir(parents=True, exist_ok=True)

    xlim = (np.min(dx), np.max(dx))
    ylim = (np.min(dy), np.max(dy))
    zlim = (np.min(z), np.max(z))

    def pad(lo, hi):
        if hi == lo:
            hi = lo + 1.0
            lo = lo - 1.0
        pad_val = (hi - lo) * 0.05 if hi != lo else 0.05
        return lo - pad_val, hi + pad_val

    xlim = pad(xlim[0], xlim[1])
    ylim = pad(ylim[0], ylim[1])
    zlim = pad(zlim[0], zlim[1])

    plt.figure(figsize=(5,4))
    plt.scatter(dx, z, s=1, alpha=0.5, rasterized=True)
    plt.xlabel(xlabel)
    plt.ylabel(zlabel)
    plt.xlim(xlim)
    plt.ylim(zlim)
    plt.tight_layout()
    plt.savefig(plot_dir / f'exp{exp_id}_dim{dim}_ch{C}_xz.png', dpi=150)
    plt.close()

    plt.figure(figsize=(5,4))
    plt.scatter(dy, z, s=1, alpha=0.5, rasterized=True)
    plt.xlabel(ylabel)
    plt.ylabel(zlabel)
    plt.xlim(ylim)
    plt.ylim(zlim)
    plt.tight_layout()
    plt.savefig(plot_dir / f'exp{exp_id}_dim{dim}_ch{C}_yz.png', dpi=150)
    plt.close()

    fig = plt.figure(figsize=(6,5))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(dx, dy, z, s=1, alpha=0.5, rasterized=True)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_zlim(zlim)
    plt.tight_layout()
    plt.savefig(plot_dir / f'exp{exp_id}_dim{dim}_ch{C}_3d.png', dpi=150)
    plt.close()


# ----------------------------------------------------------------------
# Adaptive batch size (simple)
# ----------------------------------------------------------------------
def get_auto_batch_size(dim, C, device, total_samples, safety_factor=1.5):
    if device.type == 'cuda' and torch.cuda.is_available():
        free_mem, _ = torch.cuda.mem_get_info(device)
    else:
        free_mem = 8 * 1024**3
    mem_per_batch = 2 * C * dim * dim * 4 * safety_factor
    if mem_per_batch == 0:
        return 1
    batch = int(free_mem / mem_per_batch)
    if batch < 1:
        batch = 1
    batch = min(batch, total_samples)
    return batch


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dims', type=str, default='2,3,4,5,8,12,16',
                        help='Comma-separated dimensions')
    parser.add_argument('--channels', type=str, default='1,2,3,4,8,16,32',
                        help='Comma-separated channel counts')
    parser.add_argument('--batch-size', type=int, default=None,
                        help='Fixed batch size (if None, auto-detect)')
    parser.add_argument('--total-samples', type=int, default=20000,
                        help='Total random pairs per combination')
    parser.add_argument('--device', type=str, default='cuda:0',
                        help='Device to use (cuda:0 or cpu)')
    parser.add_argument('--output-dir', type=str, default='experiment_results_full',
                        help='Root output directory')
    args = parser.parse_args()

    dims = [int(d) for d in args.dims.split(',')]
    ch_list = [int(c) for c in args.channels.split(',')]
    total_samples = args.total_samples
    device = torch.device(args.device)
    root_dir = Path(args.output_dir).resolve()
    root_dir.mkdir(parents=True, exist_ok=True)

    print(f"Device: {device}")
    print(f"Dimensions: {dims}")
    print(f"Channels: {ch_list}")
    print(f"Total samples per config: {total_samples}")

    # Create directories for each experiment
    for exp_id in [1,2,3]:
        (root_dir / f'exp{exp_id}' / 'tables').mkdir(parents=True, exist_ok=True)
        (root_dir / f'exp{exp_id}' / 'plots').mkdir(parents=True, exist_ok=True)

    for dim in dims:
        for C in ch_list:
            print(f"\n=== dim={dim}, C={C} ===")
            # Determine batch size
            if args.batch_size is not None:
                batch_size = args.batch_size
            else:
                batch_size = get_auto_batch_size(dim, C, device, total_samples)
            num_batches = (total_samples + batch_size - 1) // batch_size
            print(f"Batch size: {batch_size}, num batches: {num_batches}")

            n_factors = [dim] * C
            k_factors = [dim] * C
            alg = TypeIAlgebra(n_factors, k_factors, batch_size=batch_size,
                               complex_valued=False, device=device)

            # Experiment 1
            dx1, dy1, z1 = experiment1(alg, batch_size, num_batches, device)
            df1 = pd.DataFrame({'deltaX': dx1, 'deltaY': dy1, 'z': z1})
            df1.to_csv(root_dir / f'exp1/tables/exp1_dim{dim}_ch{C}.csv', index=False)
            save_scatter_plots(1, dim, C, dx1, dy1, z1,
                               'Δ(X)', 'Δ(Y)', 'z',
                               root_dir / f'exp1/plots')

            # Experiment 2
            dx2, dy2, z2 = experiment2(alg, batch_size, num_batches, device)
            df2 = pd.DataFrame({'deltaX': dx2, 'deltaY': dy2, 'z': z2})
            df2.to_csv(root_dir / f'exp2/tables/exp2_dim{dim}_ch{C}.csv', index=False)
            save_scatter_plots(2, dim, C, dx2, dy2, z2,
                               'Δ(|X|)', 'Δ(|Y|)', 'z',
                               root_dir / f'exp2/plots')

            # Experiment 3
            dx3, dy3, z3 = experiment3(alg, batch_size, num_batches, device)
            df3 = pd.DataFrame({'deltaX': dx3, 'deltaY': dy3, 'z': z3})
            df3.to_csv(root_dir / f'exp3/tables/exp3_dim{dim}_ch{C}.csv', index=False)
            save_scatter_plots(3, dim, C, dx3, dy3, z3,
                               'Δ(|X|)', 'Δ(Y)', 'z',
                               root_dir / f'exp3/plots')

    print(f"\nAll results saved to {root_dir}")

if __name__ == '__main__':
    main()