#!/usr/bin/env python3
"""
Example: Michelson contrast under quantum channels.

For a unital positive map Phi, lambda_min 1 <= A <= lambda_max 1 implies the same bounds for
Phi(A), hence Delta(Phi(A)) <= Delta(A). Non-unital channels (amplitude damping, generic random
channels) can increase the contrast. We estimate how often and by how much, by Monte Carlo.
"""
import torch
from torch_vn_algebra import TypeIAlgebra, channels

DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'


def main():
    torch.manual_seed(0)
    batch = 5000
    alg = TypeIAlgebra([3, 3, 3], [3, 3, 3], complex_valued=True, device=DEVICE)
    A = alg.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(batch, d, device=DEVICE),
                                      batch_size=batch, force_positive=True, force_self_adjoint=True)
    d0 = A.michelson_contrast

    tests = {
        'mixed unitary (unital)': channels.random_mixed_unitary_channel(alg, 3, batch),
        'dephasing p=0.5 (unital)': channels.dephasing_channel(alg, 0.5),
        'random Stinespring, rank 2': channels.random_channel(alg, 2, batch),
        'amplitude damping g=0.5': channels.amplitude_damping_channel(alg, 0.5),
    }
    print(f"{'channel':30s} {'unital':>7s} {'P[Delta grows]':>15s} {'max increase':>13s}")
    for name, Phi in tests.items():
        out = Phi(A)
        diff = out.michelson_contrast - d0
        print(f"{name:30s} {str(Phi.is_unital()):>7s} {(diff > 1e-6).float().mean().item():15.3f} "
              f"{diff.max().item():13.3e}")


if __name__ == '__main__':
    main()
