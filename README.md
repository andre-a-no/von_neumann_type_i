# torch_vn_algebra – Type I von Neumann algebras in PyTorch

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyTorch](https://img.shields.io/badge/PyTorch-1.9+-red.svg)](https://pytorch.org)

Batched, GPU-friendly Monte Carlo for finite-dimensional Type I von Neumann algebras

$$\mathcal M = \bigoplus_{c=1}^{C} M_{n_c}(\mathbb C)\quad\text{acting on}\quad \mathcal H=\bigoplus_{c=1}^{C}\mathbb C^{k_c}.$$

Operators are stored as one tensor of shape `(batch, C, k_max, k_max)`: the batch axis holds
Monte Carlo samples, the channel axis holds the direct summands, and block `c` occupies the
top-left `k_c × k_c` corner (the rest is zero padding).

Companion paper: I. Nikolaeva, A. Novikov, *Finite-Dimensional Type I von Neumann Algebras in
PyTorch: A GPU-Accelerated Framework for Random Block-Diagonal Operators*, arXiv:2606.15882.

## Features

- **Random operators with any spectrum** – `operator_from_eigenvalues(sampler, ...)` builds
  `A_c = U_c diag(λ) U_c*` from a user callable `sampler(dim) -> (dim,) or (batch, dim)`.
- **Unitary ensembles** – Haar on U(n)/O(n), SU(n)/SO(n), COE, CSE, random diagonal phases;
  batched (`random_unitary(n, measure, batch_size)`), and block-diagonal unitaries in the algebra
  (`random_unitary_operator`).
- **Lazy evaluation** – `X @ Y`, `X + Y`, `X.abs()`, `X.sqrt()`, `X.inv` build a recipe; nothing is
  computed until `.matrix` (or a scalar functional) is requested.
- **Functional calculus** – `abs`, `sqrt`, `inverse` (pseudo-inverse), `entropy`, `trace_a_log_a`
  via batched SVD; `lambda_max`/`lambda_min` by exact diagonalisation for blocks up to
  `alg.exact_eig_max_dim = 256`, shifted power iteration above that.
- **Three trace functionals** – `Tr_blunt` (Σ_c Tr A_c), `Tr_norm` (Σ_c Tr A_c / k_c) and the
  tracial state `tau_vN` ((1/C) Σ_c Tr A_c / k_c); norms `trace_norm`, `frobenius_norm`,
  `operator_norm`; `michelson_contrast`.
- **Hilbert space utilities** – `HilbertSpace` with inner products, orthonormal bases
  (standard / random / Haar), embedding and restriction of vectors and operators.

## Installation

```bash
git clone https://github.com/andre-a-no/von_neumann_type_i.git
cd von_neumann_type_i
pip install -e .              # library only (torch, numpy)
pip install -e ".[scripts]"   # + pandas, tqdm, matplotlib, seaborn for scripts/
pip install -e ".[test]"      # + pytest
```

## Quick start

```python
import torch
from torch_vn_algebra import TypeIAlgebra

device = 'cuda' if torch.cuda.is_available() else 'cpu'

# M = M_4(C) (+) M_6(C), acting on C^4 (+) C^6
alg = TypeIAlgebra(n_factors=[4, 6], k_factors=[4, 6], device=device)

# 100 random positive operators: eigenvalues ~ U[0, 1], Haar-random eigenvectors per block
X = alg.operator_from_eigenvalues(lambda dim: torch.rand(100, dim, device=device),
                                  batch_size=100, force_positive=True, force_self_adjoint=True)

print(X.Tr_blunt().real.mean())            # sum of block traces
print(X.Tr_norm().mean())                  # sum of normalised block traces
print(X.tau_vN().mean())                   # tracial state (1/C) sum_c Tr(A_c)/k_c
print(X.michelson_contrast.mean())         # (lmax - lmin)/(lmax + lmin), shape (100,)

# operations stay lazy until .matrix or a scalar functional is requested
U = alg.random_unitary_operator(batch_size=100)        # block-diagonal Haar unitary
Z = X.sqrt() @ U @ X.sqrt()
print(Z.trace_norm().mean(), X.entropy().mean())
```

`complex_valued=False` gives real algebras (orthogonal instead of unitary groups). Each
`force_*` flag (`self_adjoint`, `positive`, `normal`, `invertible`, `projection`) both tags the
operator and checks the property when the matrix is materialised.

More in [`examples/`](examples): random Hamiltonian with a parity symmetry, free additive
convolution, a trace inequality, unitary ensembles, a Zipf density matrix.

## Reproducing the paper

| Script | What it does |
|---|---|
| `scripts/validation.py` | Haar moments, power iteration vs. spectral gap, SVD square root accuracy (Sec. 4.1) |
| `scripts/experiment.py` | Monte Carlo experiments on trace inequalities (Sec. 5); writes raw samples, plots and `summary.csv` |
| `scripts/benchmark.py` | CPU vs GPU timings and speedup heatmaps (Sec. 4.2; needs CUDA) |

```bash
python scripts/validation.py
python scripts/experiment.py --dims 2,16 --channels 1,2,16,32 --output-dir results/experiments
python scripts/benchmark.py --cpu-threads 1 --output-dir results/benchmark_new
```

Stored outputs live in [`results/`](results): `results/experiments/` (current code) and
`results/benchmark/` (GPU timings, see its README). The corrected paper text is in
[`paper/`](paper).

## Repository layout

```
torch_vn_algebra/    library: algebra.py (TypeIAlgebra, Operator), hilbert_space.py
tests/               pytest suite (CPU, runs in CI)
examples/            short usage examples
scripts/             validation, experiments and benchmarks from the paper
results/             stored experiment and benchmark outputs
paper/               paper source and list of corrections
```

## Tests

```bash
pytest -q
```

## Limitations

- SVD dominates the cost for `k_max ≳ 200` with large batches.
- Power iteration (blocks above `exact_eig_max_dim`) converges linearly in the spectral gap and
  can fail when the dominant eigenvalues are ±λ; its stopping rule is on the change of the
  estimate, not on the error.
- No automatic differentiation guarantees, finite dimensions and Type I only.

## Citation

```bibtex
@misc{nikolaeva_novikov_2026_torch_vn_algebra,
  title  = {Finite-Dimensional Type I von Neumann Algebras in PyTorch:
            A GPU-Accelerated Framework for Random Block-Diagonal Operators},
  author = {Nikolaeva, Irina and Novikov, Andrej},
  year   = {2026},
  eprint = {2606.15882},
  archivePrefix = {arXiv},
  primaryClass  = {cs.MS}
}
```

## License

MIT.
