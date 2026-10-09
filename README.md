# torch_vn_algebra – Type I von Neumann algebras in PyTorch

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.0+-red.svg)](https://pytorch.org)

Batched operators, states, channels and constrained optimisation in finite-dimensional Type I von Neumann
algebras, on CPUs and GPUs:

$$\mathcal M = \bigoplus_{c=1}^{C} M_{k_c}(\mathbb C)\quad\text{acting on}\quad \mathcal H=\bigoplus_{c=1}^{C}\mathbb C^{k_c}.$$

(Nominal factor sizes `n_factors` $n_c\ge k_c$ describe computations in the corner $p\mathcal M_0p$ of a larger
algebra $\mathcal M_0=\bigoplus_c M_{n_c}$ with rank-$k_c$ projections $p_c$; everything refers to the corner.)

Operators are stored as one tensor of shape `(batch, C, k_max, k_max)`: the batch axis holds
samples (or optimisation starts), the sector axis holds the direct summands, and block `c` occupies the
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
pip install -e ".[scripts]"   # + scipy, pandas, tqdm, matplotlib, seaborn for scripts/
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
print(X.Tr_norm().real.mean())             # sum of normalised block traces
print(X.tau_vN().real.mean())              # tracial state (1/C) sum_c Tr(A_c)/k_c (complex in general)
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
| `Operator` | `apply_function(f)` (spectral theorem), `expm`, `log`, `power`, `eigh`, `eigenvalues` (no eigenvectors, faster); `alg.operator(tensor)`, `alg.from_blocks([...])` |
| `channels` | `Channel` (Kraus form: apply, `adjoint` = Heisenberg picture, composition `@`, `mix`, `choi`, `superoperator`, `from_superoperator`, TP / unitality checks); identity, unitary, dephasing, depolarizing, amplitude damping, Lüders measurement, conditional expectation onto the centre, random (Stinespring) and random mixed-unitary channels |
| `states` | `DensityMatrix` (an `Operator` subclass: positive, unit `Tr_blunt`; `expectation`, `mix`, `condition_on`, `density(trace)` / `from_density(..., trace)` for the `Tr_norm` and `tau_vN` conventions; preserved by trace-preserving channels and by the dynamics), random density matrices (Hilbert–Schmidt, Bures, fixed rank), Gibbs states and partition functions, tracial state, sector probabilities, Born probabilities, Lüders update, entropy, relative entropy, fidelity, trace distance, purity |
| `channels.InterSectorChannel` | CP maps between sectors and between different algebras, Φ(ρ)_d = Σ_c Σ_i K_i^{dc} ρ_c K_i^{dc*}: duals for each of the three traces, composition, sector transition matrix, `from_blocks`, `random_inter_sector_channel` |
| `composite` | `tensor_product(alg1, alg2)` (sectors = pairs (c, d)), `fused_tensor_product(alg1, alg2, fuse=add)` (pairs with equal total charge merged into one sector), `kron(A, B, alg12)`, `partial_trace` and `partial_trace_channel` (an `InterSectorChannel`; its dual is the embedding A ↦ A ⊗ 1) |
| `optimize` | constrained, batched multistart optimisation: `UnitaryParam`, `PositiveParam` / `SelfAdjointParam` with a prescribed or bounded Michelson contrast (hard constraint), `extremize` |
| `cost` | safeguards: ETA warnings for large batched decompositions and integrators, memory checks with `InsufficientMemoryError`, `set_limits`, `disabled()` |
| `chains` | `SpinChain`: spin-1/2 chains with conserved S^z, sectors N = 0..L (basis = iterated fused products), XXZ Hamiltonians with batched disorder, full algebra or single sectors, sparse form (`xxz_sparse`), momentum sectors of rings (`xxz_momentum`), reduced states and entanglement, site amplitude damping and jump operators `lowering(i)` between sectors, level-spacing ratio |
| `krylov` | `SparseSectorHamiltonian` (batch of diagonals + shared sparse hopping), batched Lanczos `ground_state`, Krylov `evolve` (exp(-itH) psi with error estimate); with `SpinChain.xxz_sparse` for sectors beyond dense matrices |
| `dynamics` | (Lindblad jumps may be `InterSectorChannel`s, e.g. particle loss) exact `propagator`, `schrodinger`, `von_neumann`; RK4 `schrodinger_rk4` (time-dependent H), `solve_operator_ode` (any dX/dt = f(t, X) in M), Lindblad `lindblad_evolve`, `lindblad_superoperator`, `lindblad_channel` |

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

## Physics background

[Physics for mathematicians](docs/physics_for_mathematicians.md) explains the physical models used in the
library (superselection, open systems, spin chains, quantum circuits, random matrices) in the language of
the algebra ⊕_c M_{k_c}(C), with pointers to the corresponding functions.
Translations: [Español](docs/physics_for_mathematicians.es.md) · [Français](docs/physics_for_mathematicians.fr.md) ·
[Deutsch](docs/physics_for_mathematicians.de.md) · [中文](docs/physics_for_mathematicians.zh.md) ·
[日本語](docs/physics_for_mathematicians.ja.md) · [Русский](docs/physics_for_mathematicians.ru.md).
The English version is authoritative.

## Reproducing the paper

**Paper v2** ([`paper/v2/main.tex`](paper/v2/main.tex), compiled `main.pdf`): every number, table and figure
comes from the scripts in [`scripts/paper`](scripts/paper), run by one driver:

```bash
bash scripts/paper/run_all.sh                    # quick mode, about 4 minutes on a 4-thread CPU
bash scripts/paper/run_all.sh full cuda          # sizes of the paper, on a GPU
cd paper/v2 && pdflatex main.tex && pdflatex main.tex
```

### Running the full experiments on a GPU server

```bash
git clone https://github.com/andre-a-no/von_neumann_type_i && cd von_neumann_type_i
git checkout claude/funny-ptolemy-confo7          # until it is merged into main
python -m venv .venv && source .venv/bin/activate
pip install torch                                  # the CUDA build matching the server's driver, see pytorch.org
pip install -e ".[scripts,test]"
python -m pytest -q                                # ~30 s, must pass
bash scripts/paper/run_all.sh check cuda           # full sizes, minimal repetitions: tests every code path and the GPU memory
nohup bash scripts/paper/run_all.sh full cuda > full_run.log 2>&1 &   # the full run (or inside tmux)
```

- The driver first prints the GPU, its memory, CUDA and cuDNN versions (also saved to
  `paper/v2/generated/logs/environment.txt`) and stops if CUDA is not available.
- Every script writes its own log to `paper/v2/generated/logs/<script>.log`, including run-time estimates
  from the `cost` module. A failing script does not stop the others; the run ends with a summary.
- `ONLY="chains baselines" bash scripts/paper/run_all.sh full cuda` reruns a subset.
- Duration. Measured on a 4-thread CPU, the `check` run takes about 35 minutes, and the full run is estimated at
  35-45 hours, most of it in `bounds_search` (~25 h) and `inequality_search` (~8 h). A GPU is several times
  faster; to estimate the full run on your GPU, multiply these CPU estimates by the ratio of the `check` times on
  the GPU and on the CPU (CPU `check` times: validate 107 s, numerics 90 s, baselines 51 s, benchmark 205 s,
  chains 580 s, inequality_search 82 s, bounds_search 1016 s; every log ends with a `[resources]` line).
- Interruptions. `inequality_search` and `bounds_search` save every finished part to `generated/partial/`; after a
  crash or reboot, rerunning the same command continues where it stopped and gives the same numbers.
- Several GPUs. The scripts are independent, e.g.
  `CUDA_VISIBLE_DEVICES=0 ONLY=bounds_search bash scripts/paper/run_all.sh full cuda` and
  `CUDA_VISIBLE_DEVICES=1 ONLY="validate_known_results numerics baselines benchmark_library chains inequality_search" bash scripts/paper/run_all.sh full cuda`
  in two terminals (the examples run only when `ONLY` is not set).
- `... full cuda tf32` repeats the run with TF32 tensor-core matrix products (Ampere or newer), to measure
  their effect; the default is true FP32.
- Results land in `paper/v2/generated/` (tables as `.tex` and `.json`) and `paper/v2/figures/`. Commit those
  two directories and the logs, then rebuild the PDF; the tables state the device and mode they came from.
- Double precision is used throughout the validation and chain scripts; a data-centre GPU (A100, H100) is much
  faster there than a consumer GPU (FP64 at 1/2 versus 1/32-1/64 of the FP32 rate).

### First version

`scripts/experiment.py` (Monte Carlo samples of the trace inequalities, stored in `results/experiments/` and
used for comparison in paper v2) and `scripts/benchmark.py` (kernel timings of the first version, stored in
`results/benchmark/`) reproduce the data of arXiv:2606.15882v1; the corrected v1 text is in [`paper/`](paper).

## Repository layout

```
torch_vn_algebra/    library: algebra.py (TypeIAlgebra, Operator), hilbert_space.py, channels.py,
                     states.py, dynamics.py, composite.py, chains.py, krylov.py, optimize.py, cost.py
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
