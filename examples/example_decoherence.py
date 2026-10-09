#!/usr/bin/env python3
"""
Example: decoherence inside superselection sectors.

A system with a parity-like symmetry (two sectors of dimension 4) evolves under a random
Hamiltonian H in M plus pure dephasing L = N (number operator in each sector). We follow the
purity, the von Neumann entropy and the Michelson contrast of the state, and check that the
sector probabilities p_c = Tr rho_c are conserved (superselection).
"""
import torch
from torch_vn_algebra import TypeIAlgebra, dynamics, states

DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'


def main():
    torch.manual_seed(1)
    batch = 200
    alg = TypeIAlgebra([4, 4], [4, 4], complex_valued=True, device=DEVICE)

    H = alg.operator_from_eigenvalues(lambda d: torch.randn(batch, d, device=DEVICE),
                                      batch_size=batch, force_self_adjoint=True)
    N = alg.from_blocks([torch.diag(torch.arange(4.0)), torch.diag(torch.arange(4.0))])
    rho0 = states.random_density_matrix(alg, batch, rank=1)          # pure within each sector

    times = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0]
    print(" t     purity   entropy  contrast   max|dp_c|")
    p0 = states.sector_probabilities(rho0)
    for t in times:
        Phi = dynamics.lindblad_channel(alg, H, [N], t=t, rates=[0.5])
        rho = Phi(rho0)
        dp = (states.sector_probabilities(rho) - p0).abs().max().item()
        print(f"{t:4.1f}   {states.purity(rho).mean().item():.4f}   "
              f"{states.von_neumann_entropy(rho).mean().item():.4f}   "
              f"{rho.michelson_contrast.mean().item():.4f}    {dp:.1e}")


if __name__ == '__main__':
    main()
