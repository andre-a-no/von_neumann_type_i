"""States on M, sector probabilities, Lueders conditioning and trace conventions."""
import torch
from torch_vn_algebra import TypeIAlgebra, DensityMatrix, states

dev = 'cuda' if torch.cuda.is_available() else 'cpu'
alg = TypeIAlgebra([2, 3], [2, 3], device=dev)        # two superselection sectors

rho = states.random_density_matrix(alg, batch_size=4)  # Hilbert-Schmidt random states
print(rho.sector_probabilities())                      # p_c = Tr rho_c, rows sum to 1

# measure the projection onto the first basis vector of sector 0
P = alg.from_blocks([torch.diag(torch.tensor([1.0, 0.0])), torch.zeros(3, 3)])
post, prob = rho.condition_on(P)                       # rho|P = P rho P / Tr(rho P)
print(prob, post.sector_probabilities())               # all weight now in sector 0

# the same state as a density w.r.t. tau_vN: omega(A) = tau_vN(d A) = Tr(rho A)
A = alg.operator_from_eigenvalues(lambda d: torch.randn(4, d), batch_size=4,
                                  force_self_adjoint=True)
d_tau = rho.density('tau_vN')
print(rho.expectation(A).real, (alg.operator(d_tau) @ A).tau_vN().real)
back = DensityMatrix.from_density(alg, d_tau, trace='tau_vN')
print((back.matrix - rho.matrix).abs().max())
