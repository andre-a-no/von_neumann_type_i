#!/usr/bin/env python3
"""
Example: SU(2) multiplets of the periodic Heisenberg chain H = J sum_i S_i . S_{i+1}, L = 12.

The U(1) symmetry (total S^z) splits the 2^12 = 4096 states into sectors N = 0..12 of size binom(12, N)
(at most 924). SU(2) splits them further: the algebra of invariant operators is (+)_S M_{m_S}, one block of
size m_S per total spin S (at most 297), and every level is (2S+1)-fold degenerate. The example compares the
block sizes and timings of the two routes, finds the ground state (energy and spin, Lieb-Mattis: S = 0) and the
singlet-triplet gap, and evaluates the spin correlations <S_0 . S_r> of the ground state from reduced matrix
elements (Wigner-Eckart).
"""
import time

import torch
from torch_vn_algebra import SpinChain
from torch_vn_algebra.su2 import SU2Chain

DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'


def main():
    L = 12
    t0 = time.perf_counter()
    su2 = SU2Chain(L, 'periodic', complex_valued=False, device=DEVICE)
    H = su2.heisenberg(J=1.0)
    E_su2 = su2.full_spectrum(H)[0]
    t_su2 = time.perf_counter() - t0

    t0 = time.perf_counter()
    u1 = SpinChain(L, 'periodic', complex_valued=False, device=DEVICE)
    E_u1 = torch.sort(torch.cat([torch.linalg.eigvalsh(u1.heisenberg(sector=N).matrix[0, 0])
                                 for N in range(L + 1)]))[0]
    t_u1 = time.perf_counter() - t0

    print(f"L = {L}, periodic Heisenberg chain, 2^L = {2 ** L} states")
    print(f"U(1)  blocks binom(L, N), N = 0..{L}: {[u1.sector_dim(N) for N in range(L + 1)]}")
    print(f"SU(2) blocks m_S, S = {su2.spins}: {su2.multiplicities}")
    print(f"largest block {max(u1.sector_dim(N) for N in range(L + 1))} -> {max(su2.multiplicities)}; "
          f"sum of cubes (dense eigensolver cost) reduced "
          f"{sum(u1.sector_dim(N) ** 3 for N in range(L + 1)) / sum(m ** 3 for m in su2.multiplicities):.0f}x")
    print(f"full spectrum: SU(2) {t_su2:.2f} s (bases, projection and eigenvalues), U(1) {t_u1:.2f} s; "
          f"max |difference| {(E_su2 - E_u1).abs().max().item():.1e}")

    gs = su2.ground_state(H)
    lowest = su2.lowest_energies(H)[0]
    print(f"ground state: E0 = {gs.energy.item():.9f} (exact -5.387390917), E0/L = {gs.energy.item() / L:.6f}, "
          f"S = {gs.spin.item():g}")
    print("lowest energy per total spin S (Lieb-Mattis: increasing in S):")
    for S, e in zip(su2.spins, lowest.tolist()):
        print(f"   S = {S:g}: {e:+.9f}")
    print(f"singlet-triplet gap E0(S=1) - E0(S=0) = {(lowest[1] - lowest[0]).item():.9f}")

    m = su2.multiplicity(gs.spin.item())
    pairs = [(0, r) for r in range(L // 2 + 1)]
    C = su2.spin_correlation(pairs)                       # S_0 . S_r in (+)_S M_{m_S}, batch = r
    corr = su2.expectation(C, gs.spin.item(), gs.coefficients[0, :m].expand(len(pairs), -1))
    print("ground-state correlations <S_0 . S_r> (Wigner-Eckart):")
    print("   " + "  ".join(f"r={r}: {c:+.6f}" for (_, r), c in zip(pairs, corr.tolist())))
    print(f"   check: <S_0 . S_1> = E0 / L: {abs(corr[1].item() - gs.energy.item() / L):.1e}")


if __name__ == '__main__':
    main()
