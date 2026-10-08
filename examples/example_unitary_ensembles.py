#!/usr/bin/env python3
"""Example: Generate Haar‑random unitary matrices from different ensembles."""
import torch
from structures.Algebra import TypeIAlgebra

def main():
    # We need an algebra instance for the random_unitary method.
    # Use a dummy algebra with 1 channel of size 1 (unused for the unitary generator).
    alg = TypeIAlgebra([1], [1], batch_size=1,
                       complex_valued=True, device='cuda')

    # U(4) – Haar measure
    U_haar = alg.random_unitary(4, measure='haar')
    print("U(4) Haar – determinant magnitude:", torch.linalg.det(U_haar).abs().item())

    # SU(3) – determinant forced to 1
    U_su = alg.random_unitary(3, measure='haar_su')
    print("SU(3) Haar – determinant:", torch.linalg.det(U_su).item())

    # Circular Orthogonal Ensemble (real orthogonal)
    if not alg.complex_valued:
        U_coe = alg.random_unitary(4, measure='coe')
        print("COE (4×4) – real orthogonal, determinant:", torch.linalg.det(U_coe).item())

    # Circular Symplectic Ensemble (even dimension)
    U_cse = alg.random_unitary(4, measure='cse')
    print("CSE (4×4) – quaternion-like, determinant magnitude:", torch.linalg.det(U_cse).abs().item())

    # Diagonal random phases
    U_diag = alg.random_unitary(5, measure='diag')
    print("Diagonal random phases (5×5):")
    print(torch.diag(U_diag).real)   # phases on the diagonal

if __name__ == '__main__':
    main()