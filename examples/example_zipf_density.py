#!/usr/bin/env python3
"""Example: Random density matrix with power‑law (Zipf) eigenvalues."""
import torch
from structures.Algebra import TypeIAlgebra

def zipf(dim):
    """Eigenvalues proportional to 1/k, normalised to sum = 1."""
    eig = 1.0 / torch.arange(1, dim + 1, dtype=torch.float32)
    return eig / eig.sum()

def main():
    # Single channel, dimension 100.
    alg = TypeIAlgebra([100], [100], batch_size=1,
                       complex_valued=False, device='cuda')

    rho = alg.operator_from_eigenvalues(zipf, batch_size=1,
                                        force_positive=True)
    entropy = rho.entropy()
    print(f"von Neumann entropy of the random density matrix: {entropy:.6f}")
    # Also show the first few eigenvalues of the density operator
    # (the eigenvalues are the same as the input, up to unitary rotation)
    # For a single sample, we can inspect the eigenvalues stored internally:
    eig = rho._eigenvalues[0][0]   # first channel, first batch
    print("First 5 eigenvalues:")
    print(eig[:5])

if __name__ == '__main__':
    main()