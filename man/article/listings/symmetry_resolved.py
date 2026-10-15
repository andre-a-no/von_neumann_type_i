"""Symmetry-resolved entanglement: XXZ ground states, rho_A in sector form, S(q) and H(p)."""
import warnings
import torch
from torch_vn_algebra import SpinChain, DensityMatrix, composite, krylov, cost
warnings.simplefilter("ignore", cost.CostWarning)   # AB is used only for its layout

L, ell, N, f64 = 16, 8, 8, torch.float64                  # ring, S^z = 0, A = sites 0..7
D = torch.tensor([-0.5, 0.0, 0.5, 1.0], dtype=f64)
chain = SpinChain(L, 'periodic', complex_valued=False)
H1 = chain.xxz_sparse(Delta=1.0, sector=N)               # Delta only scales the diagonal
H = krylov.SparseSectorHamiltonian(D[:, None] * H1.diag, H1.offdiag, f64)   # batch of 4
E0, psi, res = krylov.ground_state(H, generator=torch.Generator().manual_seed(0))

A, B = (SpinChain(n, 'open', complex_valued=False) for n in (ell, L - ell))
AB = composite.fused_tensor_product(A.algebra, B.algebra)   # sector N: pairs (q, N - q)
a, b = chain.labels[N] & (2 ** ell - 1), chain.labels[N] >> ell
qa = sum((a >> i) & 1 for i in range(ell))               # charge of A
rho = torch.zeros(4, A.algebra.C, A.algebra.k_max, A.algebra.k_max, dtype=f64)
for q, d, off in AB.product_blocks[AB.charges.index(N)]:
    m = qa == q
    Psi = torch.zeros(4, A.sector_dim(q), B.sector_dim(d), dtype=f64)
    Psi[:, A.index_of[a[m]], B.index_of[b[m]]] = psi[:, m]   # pair block of |psi>
    rho[:, q, :Psi.shape[1], :Psi.shape[1]] = Psi @ Psi.mT   # rho_A on charge q
rho_A = DensityMatrix(A.algebra, matrix=rho)

lam = [w.clamp(min=0) for w in rho_A.eigenvalues()]      # spectra of the charge sectors
p = torch.stack([w.sum(-1) for w in lam], -1)            # p(q)
T = torch.stack([torch.special.xlogy(w, w).sum(-1) for w in lam], -1)
S_q, H_p = -T / p + torch.log(p), -torch.special.xlogy(p, p).sum(-1)
S = rho_A.von_neumann_entropy()
assert torch.allclose(S, (p * S_q).sum(-1) + H_p)        # S = sum_q p(q) S(q) + H(p)
for i in range(4):
    print(f"Delta={D[i]:+.1f}: S={S[i]:.4f} H(p)={H_p[i]:.4f} S(q), q=2..6:",
          ' '.join(f"{x:.4f}" for x in S_q[i, 2:7].tolist()))
