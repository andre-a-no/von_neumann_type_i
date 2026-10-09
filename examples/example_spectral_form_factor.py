#!/usr/bin/env python3
"""
Example: spectral form factor K(t) = |Tr exp(-iHt)|^2 / N^2 of random Hamiltonians.

GUE-distributed H (random Hermitian Ginibre matrices) show the characteristic
dip - ramp - plateau, while Hamiltonians with independent (Poisson) eigenvalues have no ramp.
With C sectors (superselection), Tr exp(-iHt) is the sum of the sector traces.
"""
import torch
from torch_vn_algebra import TypeIAlgebra, dynamics

DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'


def form_factor(H, times):
    N = H.algebra.total_subspace_dim
    return torch.stack([(dynamics.propagator(H, t).trace.abs() ** 2).mean() / N ** 2 for t in times])


def main():
    torch.manual_seed(0)
    batch, k = 400, 64
    alg = TypeIAlgebra([k], [k], complex_valued=True, device=DEVICE)

    G = torch.randn(batch, 1, k, k, dtype=torch.complex64, device=DEVICE)
    gue = alg.operator((G + G.conj().transpose(-2, -1)) / (2 * k) ** 0.5, is_self_adjoint=True)
    poisson = alg.operator_from_eigenvalues(lambda d: 4 * torch.rand(batch, d, device=DEVICE) - 2,
                                            batch_size=batch, force_self_adjoint=True)

    times = [0.5, 2, 5, 10, 20, 40, 80, 160, 320]
    Kg, Kp = form_factor(gue, times), form_factor(poisson, times)
    print("    t      K_GUE    K_Poisson   (plateau 1/N = %.4f)" % (1 / k))
    for t, a, b in zip(times, Kg.tolist(), Kp.tolist()):
        print(f"{t:6.1f}   {a:.5f}    {b:.5f}")


if __name__ == '__main__':
    main()
