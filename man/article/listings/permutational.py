"""Boundary time crystal with local emission: one batched run over the local rate."""
import math
import torch
from torch_vn_algebra import dynamics
from torch_vn_algebra.permutational import PermutationalSystem, PermutationalLindbladian

N = 30
ps = PermutationalSystem(N, complex_valued=False)   # blocks j = 0..15 of size 2j+1, real
gamma = torch.tensor([0.0, 0.05, 0.1, 0.2], dtype=torch.float64)   # one rate per batch element
L = PermutationalLindbladian(ps, omega_y=2.0, collective_emission=1 / (N / 2), emission=gamma)
rho0 = ps.coherent_state(math.pi).matrix.expand(len(gamma), -1, -1, -1)   # all spins down
times = torch.arange(0, 10.001, 0.05, dtype=torch.float64)
traj = dynamics.rk4(L.rhs, rho0, times, substeps=2)  # sigma_j(t): (T, 4, C, k_max, k_max)
diag = traj.diagonal(dim1=-2, dim2=-1)               # populations of |j, m>
m_z = (diag * ps.m_values()).sum((-2, -1)) / (N / 2) # <J_z>/S, shape (T, 4)
p_j = diag[-1].sum(-1)                               # distribution of j at t = 10, (4, C)
