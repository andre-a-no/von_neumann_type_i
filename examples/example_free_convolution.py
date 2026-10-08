#!/usr/bin/env python3
"""Example: Free additive convolution of a semicircle and a Bernoulli distribution."""
import torch
from structures.Algebra import TypeIAlgebra

def semicircle(dim):
    eig = torch.randn(dim)
    return eig / eig.abs().max()

def bernoulli(dim):
    """Eigenvalues ±1 with equal probability."""
    return 2 * torch.bernoulli(0.5 * torch.ones(dim)) - 1

def main():
    # Single channel, dimension 500 for good resolution.
    alg = TypeIAlgebra([500], [500], batch_size=1,
                       complex_valued=False, device='cuda')

    A = alg.operator_from_eigenvalues(semicircle, batch_size=1,
                                      force_self_adjoint=True)
    B = alg.operator_from_eigenvalues(bernoulli, batch_size=1,
                                      force_self_adjoint=True)
    C = A + B
    eigvals = torch.linalg.eigvalsh(C.matrix[0, 0])
    print(f"Number of eigenvalues: {len(eigvals)}")
    print(f"Mean of eigenvalues: {eigvals.mean().item():.4f}")
    print(f"Standard deviation:  {eigvals.std().item():.4f}")
    print("First 10 eigenvalues:")
    print(eigvals[:10])

if __name__ == '__main__':
    main()