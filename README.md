# torch_vn_algebra – Type I von Neumann Algebras with PyTorch

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyTorch](https://img.shields.io/badge/PyTorch-1.9+-red.svg)](https://pytorch.org)

GPU‑accelerated Monte Carlo simulations for finite‑dimensional Type I von Neumann algebras.

This library provides a flexible, batched, and GPU‑friendly framework for working with direct sums of matrix algebras (Type I factors). It targets:

- **Operator algebra theorists** – numerical tests of trace inequalities, approximations of hyperfinite factors, and explorations of the Connes embedding problem.
- **Quantum physicists** – simulations of systems with superselection rules (charge, parity, angular momentum), decoherence, and random Hamiltonians with arbitrary eigenvalue distributions.
- **Computational scientists** – examples of batched linear algebra, lazy evaluation, and diagonalisation‑free functional calculus on GPUs.

## Key features

- **Direct sum (block‑diagonal) structure** – native support for channels with independent dimensions.
- **Arbitrary eigenvalue distributions** – supply any callable that generates eigenvalues.
- **Three trace functionals** – blunt trace, normalised subspace trace, and the von Neumann tracial state.
- **Functional calculus without full diagonalisation** – power iteration for extreme eigenvalues, SVD for |A|, sqrt(A), inv(A), and Tr(A log A).
- **Random unitary matrices** – Haar measure on U(n)/O(n), SU(n), COE, CSE, diagonal random phases.
- **Lazy evaluation** – postpone matrix construction until needed, save memory.
- **GPU batching** – process thousands of random operator pairs in parallel on a single GPU.
- **Modular design** – easy to extend with new samplers, unitary ensembles, or functional calculus methods.

## Installation

```bash
# Clone the repository
git clone https://gitlab.com/a.hobukov/von_neumann_type_i.git
cd von_neumann_type_i

# Install dependencies
pip install torch numpy pandas tqdm matplotlib seaborn

# (Optional) Install in editable mode
pip install -e .


# Quick start

Here is a minimal example that creates a random positive operator and computes its trace and Michelson contrast:


import torch
from torch_vn_algebra import Operator, HilbertSpace

# Set up the algebra: two channels of sizes 4 and 6
alg = HilbertSpace(k_factors=[4, 6], device='cuda' if torch.cuda.is_available() else 'cpu')

# Define a sampler for eigenvalues (uniform in [0,1])
def uniform_positive(dim):
    return torch.rand(dim, device=alg.device)

# Generate a batch of 100 random positive operators
batch_size = 100
X = Operator.from_eigenvalues(
    alg, 
    eigenvalue_sampler=uniform_positive, 
    batch_size=batch_size,
    force_positive=True,
    force_self_adjoint=True
)

# Compute traces
blunt_trace = X.trace_blunt()          # sum of traces over channels
norm_trace = X.trace_norm()            # normalised by subspace dimensions
vN_trace = X.trace_vN()                # von Neumann tracial state (average over channels)

# Compute Michelson contrast
lmax, lmin = X.lambda_max, X.lambda_min
contrast = (lmax - lmin) / (lmax + lmin)

print(f"Blunt trace: {blunt_trace.mean().item():.4f}")
print(f"von Neumann trace: {vN_trace.mean().item():.4f}")
print(f"Mean contrast: {contrast.mean().item():.4f}")