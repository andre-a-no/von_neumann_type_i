#!/usr/bin/env python3
"""Example: Random Hamiltonian with parity symmetry (two blocks of size 50)."""
import torch
from structures.Algebra import TypeIAlgebra

def semicircle(dim):
    """Eigenvalues from the semicircle law (Wigner) scaled to max 1."""
    eig = torch.randn(dim)
    return eig / eig.abs().max()

def main():
    # Two channels, each acting on a 50‑dimensional subspace.
    n_factors = [50, 50]
    k_factors = [50, 50]
    alg = TypeIAlgebra(n_factors, k_factors, batch_size=1,
                       complex_valued=False, device='cuda')

    H = alg.operator_from_eigenvalues(semicircle, batch_size=1,
                                      force_self_adjoint=True)
    # Materialise the operator (first batch, first channel)
    mat = H.matrix[0, 0]   # shape (50, 50)
    eigvals = torch.linalg.eigvalsh(mat)
    print(f"First 5 eigenvalues of the first block:\n{eigvals[:5]}")
    print(f"Smallest eigenvalue: {eigvals.min().item():.4f}")
    print(f"Largest eigenvalue:  {eigvals.max().item():.4f}")

if __name__ == '__main__':
    main()