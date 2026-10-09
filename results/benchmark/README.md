# GPU vs CPU benchmark results

Produced by an earlier version of `scripts/benchmark.py` (before the metric fixes of v0.2.0).
Configuration: `dim ∈ {2, 4, 8, 16, 32, 64}`, `channels ∈ {1, 2, 4, …, 1024}`, seeds 42, 123, 456,
10 timing repeats; NVIDIA Tesla P100 (16 GB) vs Intel Xeon E5-2650 v4. The CPU thread count was
not fixed by that script (PyTorch default).

Speedup = mean CPU time / mean GPU time per operation, measured on raw batched PyTorch kernels
applied to a tensor of shape `(channels, dim, dim)` filled with i.i.d. Gaussian entries.
With that input some metric names do not match their mathematical meaning, which the current
script fixes:

| metric | measured in this data | current script |
|---|---|---|
| `lambda_max` | `eigvalsh(A Aᵀ)` (exact diagonalisation, not power iteration) | `eigvalsh(A)`, A positive definite |
| `entropy` | entropy of normalised singular values | von Neumann entropy of A / Tr A |
| `michelson` | (max − min)/(max + min) over matrix **entries** | over eigenvalues |
| `inverse` | `torch.linalg.inv` | same, on positive definite A |

The timings of `trace`, `addition`, `multiplication`, `frobenius_norm`, `inverse`, `svd_abs`,
`trace_norm` do not depend on these definitions and remain representative.

Files:

- `full_benchmark.csv`, `pairwise_comparisons.csv`, `summary_speedups.csv`, `speedup_heatmap_*.png` – run A
  (the heatmaps were drawn from this run).
- `run_b/` – a second run with the same configuration that had been committed on top of run A as
  an unresolved git merge conflict. Its median speedup ratio to run A is 1.12; individual
  configurations differ by more than an order of magnitude, so single cells of the heatmaps
  should not be over-interpreted.
