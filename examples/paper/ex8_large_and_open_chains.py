"""Beyond dense matrices: Lanczos, Krylov evolution, momentum sectors and particle loss."""
import math
import torch
from torch_vn_algebra import SpinChain, krylov, dynamics

dev = 'cuda' if torch.cuda.is_available() else 'cpu'

# Lanczos in a sector of dimension binom(16, 8) = 12870, sparse Hamiltonian
ring = SpinChain(16, 'periodic', complex_valued=False, device=dev)
E0, psi, res = krylov.ground_state(ring.xxz_sparse(sector=8))
print(f"Heisenberg ring L=16: E0 = {E0.item():.9f} (exact -7.142296361), residual {res.item():.1e}")

# Krylov evolution of a domain wall in the XX chain, L = 14 (sector dimension 3432)
chain = SpinChain(14, device=dev)
psi0 = chain.vector_in_sector('1' * 7 + '0' * 7)
psi_t = krylov.evolve(chain.xxz_sparse(J=1.0, Delta=0.0, sector=7), psi0[None], [0.0, 3.0])
mag = (psi_t[-1, 0].abs() ** 2) @ chain.bits(7).to(torch.float64) - 0.5
print("S^z_i(t=3):", [round(x, 3) for x in mag.tolist()])

# translation symmetry: blocks of fixed momentum k = 2 pi m / L inside the sector N = 5 of a ring of 10
ring10 = SpinChain(10, 'periodic', device=dev)
Hk = ring10.xxz_momentum(J=1.0, Delta=1.0, sector=5)
print("momentum block sizes:", ring10.momentum_dims(5))
e0 = {m: w[0].min().item() for m, (w, _) in zip(Hk.algebra.charges, Hk.eigh())}
print("ground state momentum m =", min(e0, key=e0.get), "(pi for L/2 odd)")

# particle loss S-_i on every site with rate g: the sector distribution stays binomial
small = SpinChain(5, 'periodic', device=dev)
rho = dynamics.lindblad_evolve(small.basis_state('11011'), small.xxz(1.0, 0.7),
                               [small.lowering(i) for i in range(5)], [0.0, 1.0], rates=[0.4] * 5)[-1]
q = math.exp(-0.4)
print("p_N(t=1):", [round(x, 6) for x in rho.sector_probabilities()[0].tolist()])
print("binomial:", [round(math.comb(4, n) * q ** n * (1 - q) ** (4 - n), 6) for n in range(5)] + [0.0])
