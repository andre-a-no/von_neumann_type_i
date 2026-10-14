"""
Quantum-computing example: an XY mixer (QAOA for constrained problems) on 6 qubits keeps the
Hamming weight, i.e. it acts inside one sector of M = (+)_N M_binom(6,N). Amplitude damping (T1)
moves weight to lower sectors, pure dephasing does not; measuring the Hamming weight - a central
projection of M - detects the first kind of error and post-selection restores part of the fidelity.
"""
import torch
from torch_vn_algebra import SpinChain, dynamics, states

torch.manual_seed(0)                               # reproducible output
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
n, w, T = 6, 3, 2.0                                         # qubits, Hamming weight, evolution time
qubits = SpinChain(n, 'periodic', device=dev)               # ring of qubits, sectors = Hamming weights
H = qubits.xxz(J=1.0, Delta=0.0)                            # XY mixer sum (X_i X_j + Y_i Y_j) / 4
rho0 = qubits.basis_state('111000')                          # feasible initial state, weight 3
ideal = dynamics.von_neumann(rho0, H, [T])[0]                # noiseless reference (pure)
P = qubits.algebra.central([1.0 if N == w else 0.0 for N in range(n + 1)])   # Hamming-weight projector

T1_jumps = [qubits.lowering(i) for i in range(n)]            # |1> -> |0>: between sectors
dephasing_jumps = [qubits.sz(i) for i in range(n)]           # pure dephasing: inside sectors
print(" gamma1  gamma2   P(weight 3)   fidelity   fidelity after post-selection")
for g1, g2 in ((0.0, 0.0), (0.02, 0.0), (0.0, 0.02), (0.02, 0.02), (0.1, 0.1)):
    rho = dynamics.lindblad_evolve(rho0, H, T1_jumps + dephasing_jumps, [0.0, T],
                                   rates=[g1] * n + [g2] * n, substeps=80)[-1]
    post, p_ok = rho.condition_on(P)                         # Lueders rule for the central projection
    F = states.trace_of_product(rho, ideal).real.item()      # <psi_ideal| rho |psi_ideal>
    F_post = states.trace_of_product(post, ideal).real.item()
    print(f" {g1:5.2f}   {g2:5.2f}   {p_ok.item():10.4f}   {F:8.4f}   {F_post:8.4f}")
