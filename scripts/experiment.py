#!/usr/bin/env python3
"""
Batched Monte Carlo experiments on trace inequalities in M = (+)_c M_k(R).

For every pair (k, C) the script draws `--total-samples` random operator pairs whose
eigenvalues give a uniformly distributed Michelson contrast in [0, 1), conjugated by
independent Haar-random orthogonal matrices, and records

    Exp 1  (X, Y >= 0, U unitary in M):   z = |Tr(X U Y)| - Tr(X Y)
    Exp 2  (X = U X0, Y = V Y0):           z = Tr|X Y| - Tr(|X| |Y|)
    Exp 3  (X = X*, Y >= 0):               z = Tr(Y |X| Y) - Tr|Y X Y|

together with the contrasts Delta of the relevant positive operators.
Tr is the blunt trace (sum over channels), |A| = (A* A)^{1/2}.

Outputs (under --output-dir):
    exp{1,2,3}/tables/exp{i}_dim{k}_ch{C}.csv.gz  raw samples (deltaX, deltaY, z)
    exp{1,2,3}/plots/...png                     scatter plots
    summary.csv                                 share of z > 0, extremes of z, per config

Usage:
    python scripts/experiment.py --dims 2,16 --channels 1,2,16,32 --device cuda
"""
import argparse
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from tqdm import tqdm

warnings.filterwarnings("ignore", category=UserWarning)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from torch_vn_algebra import TypeIAlgebra


# ----------------------------------------------------------------------
# Eigenvalue samplers
# ----------------------------------------------------------------------
def _abs_spectrum(batch_size, dim, device):
    """max = 1, min = (1 - delta)/(1 + delta) with delta ~ U[0, 1), rest uniform in between,
    so that the Michelson contrast of the spectrum is exactly delta."""
    delta = torch.rand(batch_size, device=device)
    lam_min = (1 - delta) / (1 + delta)
    eig = torch.empty(batch_size, dim, device=device)
    eig[:, 0] = 1.0
    if dim > 1:
        eig[:, 1] = lam_min
    if dim > 2:
        low = lam_min.unsqueeze(1)
        eig[:, 2:] = low + (1 - low) * torch.rand(batch_size, dim - 2, device=device)
    return torch.sort(eig, dim=-1, descending=True)[0]


def positive_sampler(batch_size, device):
    """Positive eigenvalues with uniform Michelson contrast in [0, 1) (per channel)."""
    return lambda dim: _abs_spectrum(batch_size, dim, device)


def selfadjoint_sampler(batch_size, device):
    """Same absolute values as positive_sampler, with independent random signs."""
    def sampler(dim):
        signs = 2 * torch.randint(0, 2, (batch_size, dim), device=device).float() - 1
        return _abs_spectrum(batch_size, dim, device) * signs
    return sampler


def contrast(op):
    return op.michelson_contrast


def blunt_trace(op):
    return op.trace.real


# ----------------------------------------------------------------------
# Experiments
# ----------------------------------------------------------------------
def experiment1(alg, batch_size, num_batches, device):
    """X, Y positive, U Haar unitary in the algebra.  z = |Tr(X U Y)| - Tr(X Y)."""
    pos = positive_sampler(batch_size, device)
    out = []
    for _ in tqdm(range(num_batches), desc="Exp1", leave=False):
        X = alg.operator_from_eigenvalues(pos, batch_size=batch_size,
                                          force_positive=True, force_self_adjoint=True)
        Y = alg.operator_from_eigenvalues(pos, batch_size=batch_size,
                                          force_positive=True, force_self_adjoint=True)
        U = alg.random_unitary_operator(batch_size)
        z = (X @ U @ Y).trace.abs() - blunt_trace(X @ Y)
        out.append((contrast(X), contrast(Y), z))
    return _collect(out)


def experiment2(alg, batch_size, num_batches, device):
    """X = U X0, Y = V Y0 (X0, Y0 positive).  z = Tr|XY| - Tr(|X||Y|)."""
    pos = positive_sampler(batch_size, device)
    out = []
    for _ in tqdm(range(num_batches), desc="Exp2", leave=False):
        X0 = alg.operator_from_eigenvalues(pos, batch_size=batch_size,
                                           force_positive=True, force_self_adjoint=True)
        Y0 = alg.operator_from_eigenvalues(pos, batch_size=batch_size,
                                           force_positive=True, force_self_adjoint=True)
        X = alg.random_unitary_operator(batch_size) @ X0
        Y = alg.random_unitary_operator(batch_size) @ Y0
        absX, absY = X.abs(), Y.abs()          # |X| = (X*X)^{1/2} = X0
        z = blunt_trace((X @ Y).abs()) - blunt_trace(absX @ absY)
        out.append((contrast(X0), contrast(Y0), z))
    return _collect(out)


def experiment3(alg, batch_size, num_batches, device):
    """X self-adjoint, Y positive.  z = Tr(Y|X|Y) - Tr|YXY|."""
    sa = selfadjoint_sampler(batch_size, device)
    pos = positive_sampler(batch_size, device)
    out = []
    for _ in tqdm(range(num_batches), desc="Exp3", leave=False):
        X = alg.operator_from_eigenvalues(sa, batch_size=batch_size, force_self_adjoint=True)
        Y = alg.operator_from_eigenvalues(pos, batch_size=batch_size,
                                          force_positive=True, force_self_adjoint=True)
        absX = X.abs()
        z = blunt_trace(Y @ absX @ Y) - blunt_trace((Y @ X @ Y).abs())
        out.append((contrast(absX), contrast(Y), z))
    return _collect(out)


def _collect(out):
    return tuple(np.concatenate([o[i].detach().cpu().numpy() for o in out]) for i in range(3))


EXPERIMENTS = {
    1: (experiment1, 'Δ(X)', 'Δ(Y)'),
    2: (experiment2, 'Δ(|X|)', 'Δ(|Y|)'),
    3: (experiment3, 'Δ(|X|)', 'Δ(Y)'),
}


# ----------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------
def save_scatter_plots(exp_id, dim, C, dx, dy, z, xlabel, ylabel, plot_dir):
    plot_dir = Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)
    stem = f'exp{exp_id}_dim{dim}_ch{C}'
    for suffix, x, label in (('xz', dx, xlabel), ('yz', dy, ylabel)):
        plt.figure(figsize=(5, 4))
        plt.scatter(x, z, s=1, alpha=0.4, rasterized=True)
        plt.axhline(0.0, color='k', lw=0.6)
        plt.xlabel(label)
        plt.ylabel('z')
        plt.title(f'Exp {exp_id}: k={dim}, C={C}')
        plt.tight_layout()
        plt.savefig(plot_dir / f'{stem}_{suffix}.png', dpi=120)
        plt.close()
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(dx, dy, z, s=1, alpha=0.4, rasterized=True)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel('z')
    plt.tight_layout()
    plt.savefig(plot_dir / f'{stem}_3d.png', dpi=120)
    plt.close()


def summarize(exp_id, dim, C, dx, dy, z, tol):
    small = np.maximum(dx, dy) < 0.1
    return {
        'experiment': exp_id, 'dim': dim, 'channels': C, 'samples': len(z),
        'frac_z_positive': float(np.mean(z > tol)),
        'frac_z_negative': float(np.mean(z < -tol)),
        'z_min': float(z.min()), 'z_median': float(np.median(z)), 'z_max': float(z.max()),
        # behaviour close to scalar operators (both contrasts below 0.1)
        'z_max_contrast_lt_0.1': float(z[small].max()) if small.any() else float('nan'),
    }


def get_auto_batch_size(dim, C, device, total_samples, safety_factor=8.0):
    if device.type == 'cuda':
        free_mem, _ = torch.cuda.mem_get_info(device)
    else:
        free_mem = 2 * 1024 ** 3
    mem_per_sample = C * dim * dim * 4 * safety_factor * 4   # a handful of live operators
    return max(1, min(int(free_mem / mem_per_sample), total_samples))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--dims', type=str, default='2,16')
    parser.add_argument('--channels', type=str, default='1,2,16,32')
    parser.add_argument('--experiments', type=str, default='1,2,3')
    parser.add_argument('--batch-size', type=int, default=None, help='fixed batch size (default: auto)')
    parser.add_argument('--total-samples', type=int, default=20000)
    parser.add_argument('--device', type=str, default='cuda' if torch.cuda.is_available() else 'cpu')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--tol', type=float, default=1e-4, help='|z| below tol counts as zero')
    parser.add_argument('--no-plots', action='store_true')
    parser.add_argument('--output-dir', type=str, default='experiment_results')
    args = parser.parse_args()

    torch.manual_seed(args.seed)
    dims = [int(d) for d in args.dims.split(',')]
    channels = [int(c) for c in args.channels.split(',')]
    exps = [int(e) for e in args.experiments.split(',')]
    device = torch.device(args.device)
    root = Path(args.output_dir)
    print(f"Device: {device}; dims {dims}; channels {channels}; {args.total_samples} samples per config")

    rows = []
    for dim in dims:
        for C in channels:
            batch_size = args.batch_size or get_auto_batch_size(dim, C, device, args.total_samples)
            num_batches = -(-args.total_samples // batch_size)
            print(f"=== dim={dim}, C={C}: batch {batch_size} x {num_batches}")
            alg = TypeIAlgebra([dim] * C, [dim] * C, complex_valued=False, device=device)
            for exp_id in exps:
                func, xlabel, ylabel = EXPERIMENTS[exp_id]
                dx, dy, z = func(alg, batch_size, num_batches, device)
                dx, dy, z = dx[:args.total_samples], dy[:args.total_samples], z[:args.total_samples]
                tables = root / f'exp{exp_id}' / 'tables'
                tables.mkdir(parents=True, exist_ok=True)
                pd.DataFrame({'deltaX': dx, 'deltaY': dy, 'z': z}).to_csv(
                    tables / f'exp{exp_id}_dim{dim}_ch{C}.csv.gz', index=False, float_format='%.7g')
                if not args.no_plots:
                    save_scatter_plots(exp_id, dim, C, dx, dy, z, xlabel, ylabel, root / f'exp{exp_id}' / 'plots')
                row = summarize(exp_id, dim, C, dx, dy, z, args.tol)
                rows.append(row)
                print(f"  exp{exp_id}: z>0 in {row['frac_z_positive']:.1%}, z<0 in {row['frac_z_negative']:.1%}, "
                      f"z in [{row['z_min']:.3g}, {row['z_max']:.3g}]")

    pd.DataFrame(rows).to_csv(root / 'summary.csv', index=False)
    print(f"Results saved to {root.resolve()}")


if __name__ == '__main__':
    main()
