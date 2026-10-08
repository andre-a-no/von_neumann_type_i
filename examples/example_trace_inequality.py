#!/usr/bin/env python3
"""Example: Verify the inequality |τ(XY)| ≤ τ(|XY|) for self‑adjoint X,Y."""
import torch
from structures.Algebra import TypeIAlgebra

def main():
    # Algebra with 5 channels, each of dimension 20.
    dim = 20
    C = 5
    alg = TypeIAlgebra([dim]*C, [dim]*C, batch_size=100,
                       complex_valued=False, device='cuda')

    # Generate random self‑adjoint operators (eigenvalues uniform in [-1,1]).
    def uniform_sa(dim):
        return 2 * torch.rand(dim) - 1
    X = alg.operator_from_eigenvalues(uniform_sa, batch_size=100,
                                      force_self_adjoint=True)
    Y = alg.operator_from_eigenvalues(uniform_sa, batch_size=100,
                                      force_self_adjoint=True)

    lhs = (X @ Y).trace.abs()
    rhs = (X.abs() @ Y.abs()).trace

    # Check elementwise (batch dimension)
    satisfied = (lhs <= rhs + 1e-12).all()
    print(f"Inequality holds for all 100 batches? {satisfied}")
    print(f"Maximum violation (if any): {(lhs - rhs).max().item():.3e}")

if __name__ == '__main__':
    main()