"""Schur-Weyl duality from the symmetries alone, and the SU(2) blocks of a Heisenberg ring (sec:symmetries)."""
import torch_vn_algebra.decompose as dc
from torch_vn_algebra import SU2Chain

S4 = dc.from_group(dc.permutation_representation(4, 2))     # all 24 permutations of 4 qubits
SU2 = dc.from_generators(dc.total_spin(4))                   # S_x, S_y, S_z
print(S4.k_factors, S4.multiplicities)                       # (1, 2, 3), (5, 1, 3)
print(SU2.k_factors, SU2.multiplicities)                     # (1, 3, 5), (2, 3, 1): swapped
chain = SU2Chain(14, 'periodic', complex_valued=False)
print(chain.multiplicities)                                  # m_S for S = 0..7, largest 1001
H = chain.heisenberg(J=1.0)                                  # (+)_S H_S in (+)_S M_{m_S}
gs = chain.ground_state(H)
print(gs.energy, gs.spin, chain.lowest_energies(H))          # S = 0; E_0(S) increases with S
