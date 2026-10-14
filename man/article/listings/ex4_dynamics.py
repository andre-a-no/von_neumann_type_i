"""Spontaneous emission of a two-level atom: Lindblad evolution against the exact solution."""
import math
import torch
from torch_vn_algebra import TypeIAlgebra, DensityMatrix, dynamics

torch.manual_seed(0)                               # reproducible output
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
alg = TypeIAlgebra([2], [2], device=dev, precision='double')
omega, gamma = 1.0, 0.3
H = alg.from_blocks([torch.diag(torch.tensor([0.5 * omega, -0.5 * omega]))])   # basis |e>, |g>
L = alg.from_blocks([torch.tensor([[0.0, 0.0], [1.0, 0.0]])])                  # sigma_- = |g><e|
psi = torch.tensor([0.6, 0.8], dtype=torch.complex128)
rho0 = DensityMatrix(alg, matrix=torch.outer(psi, psi.conj())[None, None].to(dev))

times = [0.0, 1.0, 2.0, 5.0]
rk4 = dynamics.lindblad_evolve(rho0, H, [L], times, rates=[gamma], substeps=100)
for t, rho in zip(times, rk4):
    exact = dynamics.lindblad_channel(alg, H, [L], t=t, rates=[gamma])(rho0)   # exp(tL)
    p_e = 0.36 * math.exp(-gamma * t)                       # analytic excited population
    coh = 0.48 * math.exp(-gamma * t / 2)                   # analytic |rho_eg|
    m = rho.matrix[0, 0]
    print(f"t={t}: P_e {m[0, 0].real:.6f} (exact {p_e:.6f}), |rho_eg| {m[0, 1].abs():.6f} "
          f"(exact {coh:.6f}), RK4 vs exp(tL): {(rho.matrix - exact.matrix).abs().max():.1e}")
