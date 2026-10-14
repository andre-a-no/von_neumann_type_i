"""Searching for violations of |Tr(XUY)| <= Tr(XY): random sampling versus gradient ascent."""
import torch
from torch_vn_algebra import TypeIAlgebra

dev = 'cuda' if torch.cuda.is_available() else 'cpu'
torch.manual_seed(0)
k = 16
alg = TypeIAlgebra([k], [k], device=dev, precision='double')
pos = lambda d: 0.2 + torch.rand(1, d)
X = alg.operator_from_eigenvalues(pos, force_positive=True, force_self_adjoint=True)
Y = alg.operator_from_eigenvalues(pos, force_positive=True, force_self_adjoint=True)
trXY = (X @ Y).trace.real


def z(U):                                   # z(U) = |Tr(XUY)| - Tr(XY), batched over U
    return (X @ U @ Y).trace.abs() - trXY


# 1) Monte Carlo over Haar unitaries
print("sampling, 10^5 Haar U:", z(alg.random_unitary_operator(100_000)).max().item())

# 2) gradient ascent over U = U0 exp(G - G^*), 64 random starts optimised in parallel;
#    the exponential chart is re-centred at the current point every 50 steps
U0 = alg.random_unitary_operator(64).matrix
for _ in range(10):
    G = torch.zeros(64, 1, k, k, dtype=alg.hilbert.dtype, device=dev, requires_grad=True)
    opt = torch.optim.Adam([G], lr=0.05)
    for step in range(50):
        U = alg.operator(U0) @ alg.operator(G - G.conj().transpose(-2, -1)).expm()
        loss = -z(U).sum()
        opt.zero_grad()
        loss.backward()
        opt.step()
    U0 = U.matrix.detach()
print("gradient ascent:     ", z(alg.operator(U0)).max().item())

# exact supremum: sup_U |Tr(U YX)| = ||YX||_1
print("exact sup:           ", ((Y @ X).trace_norm() - trXY).item())
