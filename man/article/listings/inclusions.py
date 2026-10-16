"""Inclusions, conditional expectations and indices (sec:symmetries)."""
import math, torch
from torch_vn_algebra import TypeIAlgebra, inclusions as inc

torch.manual_seed(0)                                         # reproducible output

N = TypeIAlgebra([2, 1], [2, 1], precision='double')         # M_2 (+) C
M = TypeIAlgebra([4, 3], [4, 3], precision='double')         # M_4 (+) M_3
I = inc.Inclusion(N, M, [[1, 2], [1, 1]], unitary=M.random_unitary_operator(3))
print(I.index())                                             # ||Lambda||^2 = (7 + sqrt 45)/2
E = I.conditional_expectation('markov')                      # tau-preserving E: M -> N
y = M.operator(torch.randn(3, M.C, M.k_max, M.k_max, dtype=M.dtype))
x = N.operator(torch.randn(3, N.C, N.k_max, N.k_max, dtype=N.dtype))
print((E(I.embed(x) @ y).matrix - (x @ E(y)).matrix).abs().max())   # bimodule property
print(E.index().matrix[0, :, 0, 0].real)                     # Watatani index: [M:N] 1
print(E.pimsner_popa_constant(), E.pimsner_popa_constant(complete=True))
J, e = I.basic_construction(), I.jones_projection('markov')  # M c M_1 = <M, e_N>
eye = e.matrix @ J.embed(y).matrix @ e.matrix                # e_N y e_N = E(y) e_N
print((eye - J.embed(I.embed(E(y))).matrix @ e.matrix).abs().max())
for L, T in enumerate(inc.chain_tower(6)):                   # A_L c A_{L+1}, x -> x (x) 1_2
    print(L, T.index(), 4 * math.cos(math.pi / (2 * L + 4)) ** 2)
