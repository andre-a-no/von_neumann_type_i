"""
CPU vs GPU timing of the batched kernels behind the library operations.

Each test tensor has shape (channels, dim, dim) and holds random *positive definite*
matrices A = G G^T / dim + 0.1 I, so that every metric below is well defined
(inverse, entropy, Michelson contrast). Speedup = t_CPU / t_GPU, including the
.item() synchronisation but excluding host-to-device transfer.

Results in archive/v1/benchmark_results/ were produced by an earlier version of this script
that used Gaussian (non-symmetric) matrices; see archive/v1/benchmark_results/README.md.

Usage:  python scripts/benchmark.py --dims 2,4,8,16,32,64 --channels 1,2,4,...,1024 --cpu-threads 1
"""
import torch
import numpy as np
import pandas as pd
import timeit
from itertools import product
import argparse
import os
import warnings
warnings.filterwarnings('ignore')

# ---------- METRIC DEFINITIONS (WITH BATCH/MULTI-CHANNEL HANDLING) ----------
def trace_batch(A):
    """Matrix trace (averaged over channels)"""
    return torch.diagonal(A, dim1=-2, dim2=-1).sum(-1).mean().item()

def lambda_max_batch(A):
    """Largest eigenvalue (exact batched diagonalisation, torch.linalg.eigvalsh)"""
    eigvals = torch.linalg.eigvalsh(A)
    return eigvals.max(dim=-1)[0].mean().item()

def svd_abs_batch(A):
    """Mean absolute value of singular values"""
    s = torch.linalg.svdvals(A)
    return s.abs().mean().item()

def inverse_batch(A):
    """Matrix inverse (only for square). Return mean absolute element value."""
    inv = torch.linalg.inv(A)
    return inv.abs().mean().item()

def frobenius_norm_batch(A):
    """Frobenius norm"""
    return torch.linalg.norm(A, dim=(-2, -1)).mean().item()

def trace_norm_batch(A):
    """Nuclear norm = sum of singular values"""
    s = torch.linalg.svdvals(A)
    return s.sum(dim=-1).mean().item()

def addition_batch(A):
    """Element-wise addition A + A (return mean absolute value)"""
    return (A + A).abs().mean().item()

def multiplication_batch(A):
    """Matrix multiplication A @ A"""
    return (A @ A).abs().mean().item()

def von_neumann_entropy(A):
    """von Neumann entropy -Tr(rho log rho) of rho = A / Tr A"""
    w = torch.linalg.eigvalsh(A).clamp_min(1e-30)
    p = w / w.sum(dim=-1, keepdim=True)
    return (-(p * torch.log(p)).sum(dim=-1)).mean().item()

def michelson_contrast(A):
    """Michelson contrast (lambda_max - lambda_min) / (lambda_max + lambda_min) per matrix"""
    w = torch.linalg.eigvalsh(A)
    lmin, lmax = w[..., 0], w[..., -1]
    return ((lmax - lmin) / (lmax + lmin)).mean().item()

METRICS = {
    'trace': trace_batch,
    'lambda_max': lambda_max_batch,
    'svd_abs': svd_abs_batch,
    'inverse': inverse_batch,
    'frobenius_norm': frobenius_norm_batch,
    'trace_norm': trace_norm_batch,
    'addition': addition_batch,
    'multiplication': multiplication_batch,
    'entropy': von_neumann_entropy,
    'michelson': michelson_contrast
}

# ---------- TIME MEASUREMENT WITH DATA TRANSFER ----------
def measure_time(func, tensor, device='cpu', repeats=10):
    if device == 'cuda' and not torch.cuda.is_available():
        raise RuntimeError("CUDA not available")
    x = tensor.to(device)
    def wrapper():
        return func(x)
    for _ in range(3):
        wrapper()
    times = timeit.repeat(wrapper, number=1, repeat=repeats)
    mean_t = np.mean(times)
    std_t = np.std(times)
    val = wrapper()
    return mean_t, std_t, val

# ---------- GENERATE TENSOR WITH FIXED SEED ----------
def generate_tensor(dim, channels, seed):
    torch.manual_seed(seed)
    G = torch.randn(channels, dim, dim)
    return G @ G.transpose(-2, -1) / dim + 0.1 * torch.eye(dim)

# ---------- MAIN BENCHMARK LOOP ----------
def run_benchmark(dim_list, channels_list, seeds, repeats=10, device_cpu='cpu', device_gpu='cuda'):
    results = []
    for dim, ch, seed in product(dim_list, channels_list, seeds):
        print(f"Running dim={dim}, channels={ch}, seed={seed}")
        A = generate_tensor(dim, ch, seed)
        
        row_cpu = {'dim': dim, 'channels': ch, 'seed': seed, 'device': device_cpu}
        for name, func in METRICS.items():
            t_mean, t_std, val = measure_time(func, A, device=device_cpu, repeats=repeats)
            row_cpu[f'{name}_time_mean'] = t_mean
            row_cpu[f'{name}_time_std'] = t_std
            row_cpu[f'{name}_value'] = val
        results.append(row_cpu)
        
        row_gpu = {'dim': dim, 'channels': ch, 'seed': seed, 'device': device_gpu}
        for name, func in METRICS.items():
            t_mean, t_std, val = measure_time(func, A, device=device_gpu, repeats=repeats)
            row_gpu[f'{name}_time_mean'] = t_mean
            row_gpu[f'{name}_time_std'] = t_std
            row_gpu[f'{name}_value'] = val
        results.append(row_gpu)
    
    return pd.DataFrame(results)

# ---------- POST-PROCESSING ----------
def pairwise_comparison(df, metric_names):
    cpu = df[df['device'] == 'cpu'].copy()
    gpu = df[df['device'] == 'cuda'].copy()
    merged = cpu.merge(gpu, on=['dim', 'channels', 'seed'], suffixes=('_cpu', '_cuda'))
    for m in metric_names:
        merged[f'{m}_speedup'] = merged[f'{m}_time_mean_cpu'] / merged[f'{m}_time_mean_cuda']
        merged[f'{m}_value_diff_abs'] = np.abs(merged[f'{m}_value_cpu'] - merged[f'{m}_value_cuda'])
        merged[f'{m}_value_diff_rel'] = merged[f'{m}_value_diff_abs'] / (np.abs(merged[f'{m}_value_cpu']) + 1e-12)
    return merged

def summary_speedups(pair_df, metric_names):
    rows = []
    for (dim, ch), group in pair_df.groupby(['dim', 'channels']):
        for m in metric_names:
            sp = group[f'{m}_speedup']
            rows.append({
                'dim': dim,
                'channels': ch,
                'operation': m,
                'mean_speedup': sp.mean(),
                'median_speedup': sp.median(),
                'min_speedup': sp.min(),
                'max_speedup': sp.max(),
                'std_speedup': sp.std(),
                'gpu_always_faster': (sp > 1).all(),
                'gpu_never_faster': (sp <= 1).all()
            })
    return pd.DataFrame(rows)

# ---------- ESTIMATE MAXIMUM CHANNELS (BASED ON GPU MEMORY) ----------
def estimate_max_channels(dim, safety_factor=0.5, max_test_channels=1024):
    if not torch.cuda.is_available():
        return 1
    free_mem, _ = torch.cuda.mem_get_info()
    available_mem = free_mem * safety_factor
    bytes_per_element = 4
    overhead = 3
    max_elements = available_mem / (bytes_per_element * overhead)
    max_channels = int(max_elements // (dim * dim))
    return min(max_channels, max_test_channels)

# ---------- HEATMAP GENERATION ----------
def plot_speedup_heatmap(summary_df, metric='inverse', save_path=None):
    """
    Creates a heatmap of GPU/CPU speedup for a given metric over dim and channels.
    summary_df: output of summary_speedups().
    metric: operation name (e.g., 'inverse', 'trace', etc.).
    save_path: if None, auto-generates filename.
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError as e:
        print(f"Warning: cannot generate heatmap because {e.name} is not installed.")
        return
    
    sub = summary_df[summary_df['operation'] == metric]
    if sub.empty:
        print(f"Warning: metric '{metric}' not found in summary. Available: {summary_df['operation'].unique()}")
        return
    
    pivot = sub.pivot(index='dim', columns='channels', values='mean_speedup')
    pivot.sort_index(ascending=False, inplace=True)
    pivot = pivot.reindex(sorted(pivot.columns), axis=1)
    
    if save_path is None:
        save_path = f'speedup_heatmap_{metric}.png'
    
    plt.figure(figsize=(12, 6))
    ax = sns.heatmap(pivot, annot=True, fmt='.2f', cmap='viridis',
                     cbar_kws={'label': f'GPU/CPU speedup ({metric})'})
    ax.set_title(f'Heatmap of {metric} speedup over dimension and channels')
    ax.set_xlabel('Number of channels')
    ax.set_ylabel('Dimension (dim)')
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.close()
    print(f"Heatmap saved as {save_path}")

# ---------- MAIN FUNCTION ----------
def main():
    parser = argparse.ArgumentParser(description='Benchmark matrix operations on CPU vs GPU with heatmaps for all metrics.')
    parser.add_argument('--dims', type=str, default='2,4,8,16,32,64',
                        help='Comma-separated dimensions')
    parser.add_argument('--channels', type=str, default=None,
                        help='Comma-separated channel counts (if None, auto-detected from GPU memory)')
    parser.add_argument('--seeds', type=str, default='42,123,456',
                        help='Comma-separated random seeds')
    parser.add_argument('--repeats', type=int, default=10,
                        help='Number of time measurement repetitions per operation')
    parser.add_argument('--no-heatmap', action='store_true',
                        help='Disable heatmap generation')
    parser.add_argument('--cpu-threads', type=int, default=None,
                        help='torch.set_num_threads for the CPU runs (default: PyTorch default)')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Directory for CSV files and heatmaps')
    args = parser.parse_args()
    if args.cpu_threads is not None:
        torch.set_num_threads(args.cpu_threads)
    if not torch.cuda.is_available():
        raise SystemExit("benchmark.py compares CPU with GPU and needs a CUDA device")
    print(f"CPU threads: {torch.get_num_threads()}")
    out = args.output_dir
    os.makedirs(out, exist_ok=True)
    
    dims = [int(d) for d in args.dims.split(',')]
    seeds = [int(s) for s in args.seeds.split(',')]
    repeats = args.repeats
    
    # Auto-detect channels if not provided
    if args.channels is None:
        max_ch = estimate_max_channels(max(dims))
        if max_ch < 1:
            max_ch = 1
        channels_list = [1]
        while channels_list[-1] * 2 <= max_ch:
            channels_list.append(channels_list[-1] * 2)
    else:
        channels_list = [int(c) for c in args.channels.split(',')]
        max_ch = max(channels_list)
    
    if torch.cuda.is_available():
        print(f"Available GPU memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
    print(f"Maximum number of channels for dim={max(dims)}: {max_ch}")
    print(f"Channels to be used: {channels_list}")
    
    print("\nStarting benchmark...")
    df_full = run_benchmark(dims, channels_list, seeds, repeats)
    df_full.to_csv(os.path.join(out, 'full_benchmark.csv'), index=False)
    print("Saved full_benchmark.csv")
    
    metric_names = list(METRICS.keys())
    pair_df = pairwise_comparison(df_full, metric_names)
    pair_df.to_csv(os.path.join(out, 'pairwise_comparisons.csv'), index=False)
    print("Saved pairwise_comparisons.csv")
    
    summary_df = summary_speedups(pair_df, metric_names)
    summary_df.to_csv(os.path.join(out, 'summary_speedups.csv'), index=False)
    print("Saved summary_speedups.csv")
    
    print("\n=== AVERAGE GPU/CPU SPEEDUP PER OPERATION ===")
    for m in metric_names:
        avg_speed = pair_df[f'{m}_speedup'].mean()
        print(f"{m:20s}: {avg_speed:.2f}x")
    
    print("\n=== CHECKING RESULT AGREEMENT (average relative difference) ===")
    for m in metric_names:
        avg_rel_diff = pair_df[f'{m}_value_diff_rel'].mean()
        print(f"{m:20s}: {avg_rel_diff:.2e}")
    
    # Generate heatmaps for ALL metrics (unless disabled)
    if not args.no_heatmap:
        for metric in metric_names:
            plot_speedup_heatmap(summary_df, metric=metric,
                                 save_path=os.path.join(out, f'speedup_heatmap_{metric}.png'))

if __name__ == '__main__':
    main()