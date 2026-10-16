"""H(Delta_3|U Delta_3 U^*) for a Haar-random U by batched ascent over partitions of unity (sec:cs)."""
import torch
from torch_vn_algebra import TypeIAlgebra, inclusions as inc, optimize as opt

n, B, m = 3, 64, 9                                      # M_3; B starts; m elements per partition
M = TypeIAlgebra([n], [n], precision='double')
N = TypeIAlgebra([1] * n, [1] * n, precision='double')  # C^n, embedded as the diagonal
U = M.random_unitary_operator(1)
EP = inc.Inclusion(N, M, [[1] * n]).conditional_expectation('norm')             # tau = Tr/n
EQ = inc.Inclusion(N, M, [[1] * n], unitary=U).conditional_expectation('norm')  # Q = U P U^*

def tau_eta(E, x):                                      # sum_c nu_c Tr eta(E(x)_c)
    return sum(E.nu[:, c] * torch.special.entr(w.clamp_min(1e-300)).sum(-1)
               for c, w in enumerate(E(x).eigenvalues()))

Z = torch.randn(B, m, n, n, dtype=M.dtype, requires_grad=True)
class Partition(opt.Param):
    def parameters(self): return [Z]

def objective():                                        # x_j = S^{-1/2} Z_j Z_j^* S^{-1/2}
    G = Z @ Z.mH
    S = M.operator(G.sum(1, keepdim=True), is_self_adjoint=True, is_positive=True).power(-0.5).matrix
    x = M.operator((S @ G @ S).reshape(B * m, 1, n, n), is_self_adjoint=True, is_positive=True)
    return (tau_eta(EQ, x) - tau_eta(EP, x)).reshape(B, m).sum(-1)

H = opt.extremize(objective, [Partition()], rounds=8, steps=150).max()
h = torch.special.entr(U.matrix[0, 0].abs() ** 2).sum() / n                     # Choda's h(U) <= H
print(f"H = {float(H):.6f} >= h(U) = {float(h):.6f}")
