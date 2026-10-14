"""Spin chains with conserved S^z: sectors, disorder averages and entanglement."""
import torch
from torch_vn_algebra import SpinChain, DensityMatrix

torch.manual_seed(0)                               # reproducible output
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
chain = SpinChain(L=10, boundary='periodic', device=dev)        # sectors N = 0..10, dims binom(10, N)
print(chain.algebra.k_factors)

# Heisenberg ring in a uniform field: the ground state changes sector as h grows
for h in (0.0, 1.0, 2.5):
    w = [blk[0] for blk in chain.heisenberg(h=h).eigenvalues()]
    print(f"h={h}: ground state in sector N={min(range(11), key=lambda N: w[N].min())}")

# random fields h_i ~ U[-W, W]: 100 disorder realisations in one batch, sector N = 5 only
for W in (0.5, 8.0):
    H = chain.random_field_heisenberg(W, batch_size=100, sector=5)
    E = H.eigenvalues()[0]                                      # (100, 252)
    r = SpinChain.level_spacing_ratio(E[:, 63:189]).mean()
    print(f"W={W}: <r> = {r:.3f}   (GOE 0.531, Poisson 0.386)")

# entanglement entropy of the ground state of the clean chain, left block of n sites
E, V = chain.heisenberg(sector=5).eigh()[0]
psi = V[0, :, 0]
rho = DensityMatrix(chain.sector_algebra(5), matrix=torch.outer(psi, psi.conj())[None, None])
print([round(chain.entanglement_entropy(rho, n, sector=5).item(), 4) for n in range(1, 6)])
