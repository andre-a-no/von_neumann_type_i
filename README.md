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

### Channels, states and dynamics

```python
import torch
from torch_vn_algebra import TypeIAlgebra, channels, dynamics, states

alg = TypeIAlgebra([4, 4], [4, 4], complex_valued=True)          # two superselection sectors
rho = states.random_density_matrix(alg, batch_size=100)            # HS-random states on M
H = alg.operator_from_eigenvalues(lambda d: torch.randn(100, d), batch_size=100,
                                  force_self_adjoint=True)

Phi = channels.random_channel(alg, kraus_rank=2, batch_size=100)   # Haar Stinespring isometry
print(Phi.is_trace_preserving(), Phi.is_unital())
print(states.relative_entropy(Phi(rho), Phi(states.gibbs_state(H, 1.0))))

psi0 = torch.zeros(100, 2, 4, dtype=torch.complex64)               # state vectors, block layout
psi0[:, 0, 0] = 1.0
psi_t = dynamics.schrodinger(psi0, H, times=[0.0, 1.0, 2.0])       # exact, shape (3, 100, 2, 4)

N = alg.from_blocks([torch.diag(torch.arange(4.0))] * 2)           # dephasing jump operator
rho_t = dynamics.lindblad_evolve(rho, H, jumps=[N], times=[0, 1, 2])  # RK4, list of states
Lt = dynamics.lindblad_channel(alg, H, jumps=[N], t=2.0)            # exact exp(tL) as a Channel
print((Lt(rho) - rho_t[-1]).frobenius_norm().max())
```

| Module | Contents |
|---|---|
| `Operator` | `apply_function(f)` (spectral theorem), `expm`, `log`, `power`, `eigh`; `alg.operator(tensor)`, `alg.from_blocks([...])` |
| `channels` | `Channel` (Kraus form: apply, `adjoint` = Heisenberg picture, composition `@`, `mix`, `choi`, `superoperator`, `from_superoperator`, TP / unitality checks); identity, unitary, dephasing, depolarizing, amplitude damping, Lüders measurement, conditional expectation onto the centre, random (Stinespring) and random mixed-unitary channels |
| `states` | `DensityMatrix` (an `Operator` subclass: positive, unit `Tr_blunt`; `expectation`, `mix`, `condition_on`, `density(trace)` / `from_density(..., trace)` for the `Tr_norm` and `tau_vN` conventions; preserved by trace-preserving channels and by the dynamics), random density matrices (Hilbert–Schmidt, Bures, fixed rank), Gibbs states and partition functions, tracial state, sector probabilities, Born probabilities, Lüders update, entropy, relative entropy, fidelity, trace distance, purity |
| `channels.InterSectorChannel` | CP maps between sectors and between different algebras, Φ(ρ)_d = Σ_c Σ_i K_i^{dc} ρ_c K_i^{dc*}: duals for each of the three traces, composition, sector transition matrix, `from_blocks`, `random_inter_sector_channel` |
| `composite` | `tensor_product(alg1, alg2)` (sectors = pairs (c, d)), `fused_tensor_product(alg1, alg2, fuse=add)` (pairs with equal total charge merged into one sector), `kron(A, B, alg12)`, `partial_trace` and `partial_trace_channel` (an `InterSectorChannel`; its dual is the embedding A ↦ A ⊗ 1) |
| `optimize` | constrained, batched multistart optimisation: `UnitaryParam`, `PositiveParam` / `SelfAdjointParam` with a prescribed or bounded Michelson contrast (hard constraint), `extremize` |
| `cost` | safeguards: ETA warnings for large batched decompositions and integrators, memory checks with `InsufficientMemoryError`, `set_limits`, `disabled()` |
| `chains` | `SpinChain`: spin-1/2 chains with conserved S^z, sectors N = 0..L (basis = iterated fused products), XXZ Hamiltonians with batched disorder, full algebra or single sectors, reduced states and entanglement, site amplitude damping between sectors, level-spacing ratio |
| `dynamics` | exact `propagator`, `schrodinger`, `von_neumann`; RK4 `schrodinger_rk4` (time-dependent H), `solve_operator_ode` (any dX/dt = f(t, X) in M), Lindblad `lindblad_evolve`, `lindblad_superoperator`, `lindblad_channel` |

Lazy operators are evaluated once and cached; after that they drop their recipe, so the
intermediate results of an expression such as `X @ U @ Y` can be garbage-collected. Laziness
defers and skips work but does not fuse operations; where it matters the library uses direct
formulas instead (e.g. `Tr(ρA)` in O(k²) via `states.trace_of_product`).

Channels map each sector to itself, so the Heisenberg dual is the same for `Tr_blunt`, `Tr_norm`
and `tau_vN`, and Hamiltonian / Lindblad dynamics with generators in M conserve the sector
probabilities `Tr rho_c`. Unitary dynamics requires `complex_valued=True`.

**Number field and precision** are two independent switches of every algebra:
`complex_valued` (real / complex) and `precision` (`'single'` / `'double'`), i.e. float32, float64,
complex64 or complex128; `alg.like(complex_valued=..., precision=...)` and `op.cast(alg2)` move between
them, `alg.bytes_per_operator(B)` reports the memory. Rules of thumb (details and measurements in the
paper, `scripts/paper/numerics.py`):

| use | when |
|---|---|
| real | the problem is real: real symmetric spectra, ground states without complex couplings, GOE statistics, orthogonal conjugation. Real is *exact* there and halves memory. Note that it changes the model otherwise (O(n) instead of U(n), GOE instead of GUE) |
| complex | unitary dynamics e^{-iHt}, Lindblad with a Hamiltonian, Haar U(n) / CUE / CSE / SU(n), random pure states, broken time reversal |
| single | Monte Carlo statistics whose statistical error is far above 1e-6, exploratory optimisation |
| double | effects that are small compared with the operators: inequalities near equality, small gaps, entropies of nearly pure states, contrast near 1, exact-result checks, long RK4 runs. On consumer GPUs FP64 is 32-64x slower than FP32 |

`complex_valued=False` gives real algebras (orthogonal instead of unitary groups). Each
`force_*` flag (`self_adjoint`, `positive`, `normal`, `invertible`, `projection`) both tags the
operator and checks the property when the matrix is materialised.

More in [`examples/`](examples): random Hamiltonian with a parity symmetry, free additive
convolution, a trace inequality, unitary ensembles, a Zipf density matrix, decoherence inside
superselection sectors, the Michelson contrast under unital and non-unital channels, and the
spectral form factor of GUE vs. Poisson Hamiltonians.

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

**Paper v2** ([`paper/v2/main.tex`](paper/v2/main.tex), compiled `main.pdf`): a usage-centred description
with runnable listings ([`examples/paper`](examples/paper), run by the tests), validation against exact
results, sampling vs. constrained optimisation and library-level benchmarks. All numbers and figures come
from `scripts/paper/`:

```bash
bash scripts/paper/run_all.sh              # quick mode on CPU, a few minutes
bash scripts/paper/run_all.sh full cuda    # sizes of the paper, on a GPU
cd paper/v2 && pdflatex main.tex && pdflatex main.tex
```

Stored outputs live in [`results/`](results): `results/experiments/` (current code) and
`results/benchmark/` (GPU timings, see its README). The corrected paper text is in
[`paper/`](paper).

## Repository layout

```
torch_vn_algebra/    library: algebra.py (TypeIAlgebra, Operator), hilbert_space.py, channels.py,
                     states.py, dynamics.py, composite.py, chains.py, optimize.py, cost.py
tests/               pytest suite (CPU, runs in CI)
examples/            short usage examples
scripts/             validation, experiments and benchmarks from the paper
results/             stored experiment and benchmark outputs
paper/               corrected text of arXiv:2606.15882v1 + list of corrections; paper/v2: new text
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
- `Channel`, Hamiltonians and Lindblad generators preserve the sectors; maps between sectors are
  `InterSectorChannel`s. Only abelian charges are fused; non-abelian symmetries are not used to reduce
  blocks further, and operators on all sectors of a chain are stored padded (work in single sectors for
  large chains).
- `lindblad_channel`, `choi` and `superoperator` work with k_c² × k_c² matrices per block, which
  limits them to blocks of a few tens; use `lindblad_evolve` for larger blocks.
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
