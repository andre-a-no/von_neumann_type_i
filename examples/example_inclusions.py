#!/usr/bin/env python3
"""
Example: inclusions N c M of finite-dimensional von Neumann algebras.

1. The index [M : N] = ||Lambda||^2 of a few standard inclusions.
2. The U(1) spin-chain tower A_L c A_{L+1} (x -> x (x) 1_2): index 4 cos^2(pi / (2L + 4)) -> 4.
3. The trace-preserving conditional expectation E: M -> N, its defining properties, the Watatani index
   (= [M : N] 1 for the Markov trace) and the Pimsner-Popa constants.
4. The basic construction M c M_1 = <M, e_N>: inclusion matrix Lambda^T, same index, e_N y e_N = E(y) e_N.
Runs in a few seconds on a CPU.
"""
import math

import numpy as np
import torch

from torch_vn_algebra import TypeIAlgebra
from torch_vn_algebra import inclusions as inc


def err(a, b) -> float:
    return (a - b).abs().max().item()


def main():
    torch.manual_seed(0)
    opts = dict(complex_valued=True, precision='double')

    print("1. Index [M : N] = ||Lambda||^2")
    examples = {
        "C c M_4 (scalars)": inc.scalars(TypeIAlgebra([4], [4], **opts)),
        "C^4 c M_4 (diagonal masa)": inc.diagonal(4),
        "M_2 (x) 1 c M_2 (x) M_3": inc.from_tensor_product(TypeIAlgebra([2], [2], **opts), 3),
        "M_2 + M_1 c M_4 + M_3 (Lambda = [[1, 2], [1, 1]])":
            inc.Inclusion(TypeIAlgebra([2, 1], [2, 1], **opts), TypeIAlgebra([4, 3], [4, 3], **opts), [[1, 2], [1, 1]]),
    }
    for name, I in examples.items():
        print(f"   {name:50s} {I.index():9.5f}")

    print("\n2. Spin-chain tower A_L c A_{L+1}: the Bratteli diagram is Pascal's triangle")
    tower = inc.chain_tower(8)
    print(f"   {'L':>3s} {'[A_L+1 : A_L]':>15s} {'4 cos^2(pi/(2L+4))':>20s}")
    for L, I in enumerate(tower):
        print(f"   {L:3d} {I.index():15.12f} {4 * math.cos(math.pi / (2 * L + 4)) ** 2:20.12f}")
    print("   levels:", inc.bratteli_diagram(tower[:4]).levels)

    print("\n3. Conditional expectation for an inclusion conjugated by a random unitary (batch of 4)")
    N = TypeIAlgebra([2, 1, 3], [2, 1, 3], **opts)
    lam = np.array([[1, 2, 0], [0, 1, 1], [2, 0, 1]])
    m = [int(v) for v in lam @ np.array(N.k_factors)]
    M = TypeIAlgebra(m, m, **opts)
    I = inc.Inclusion(N, M, lam, unitary=M.random_unitary_operator(4))
    rand = lambda alg: alg.operator(torch.randn(4, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype))
    x, n1, n2, y = rand(N), rand(N), rand(N), rand(M)
    t, s, beta = I.markov_trace()
    print(f"   Lambda = {lam.tolist()}, block sizes of M = {m}, [M : N] = {beta:.6f}")
    print(f"   Markov trace: t = {np.round(t.numpy(), 5).tolist()}, s = Lambda^T t = {np.round(s.numpy(), 5).tolist()}")
    for trace in ('blunt', 'norm', 'vn', 'markov'):
        E = I.conditional_expectation(trace)
        w, nu = I.trace_weights(trace), I.induced_trace(trace)
        tau = lambda op: (torch.diagonal(op.matrix, dim1=-2, dim2=-1).sum(-1) * w).sum(-1)
        tau_N = lambda op: (torch.diagonal(op.matrix, dim1=-2, dim2=-1).sum(-1) * nu).sum(-1)
        checks = [err(E(I.embed(x)).matrix, x.matrix),                                         # E o iota = id
                  err(E(I.embed(n1) @ y @ I.embed(n2)).matrix, (n1 @ E(y) @ n2).matrix),       # bimodule
                  err(tau(I.embed(E(y))), tau(y)),                                             # tau o E = tau
                  err(tau(I.embed(x) @ y), tau_N(x @ E(y)))]                                   # adjoint of iota
        ind = E.index().matrix[0, :, 0, 0].real
        print(f"   {trace:6s}: max error of (E o iota = id, bimodule, tau o E = tau, adjointness) = {max(checks):.1e}; "
              f"Ind E = {np.round(ind.numpy(), 4).tolist()}; lambda(E) = {E.pimsner_popa_constant().item():.4f}, "
              f"complete: {E.pimsner_popa_constant(complete=True).item():.4f}")

    print("\n4. Basic construction and Jones tower")
    J = I.jones_tower(4)
    print("   inclusion matrices:", [T.matrix.tolist() for T in J])
    print("   indices:", [round(T.index(), 10) for T in J])
    e = I.jones_projection('markov')
    j, E = I.basic_construction(), I.conditional_expectation('markov')
    eye_ = e.matrix @ j.embed(y).matrix @ e.matrix
    print(f"   |e y e - E(y) e| = {err(eye_, j.embed(I.embed(E(y))).matrix @ e.matrix):.1e}")
    R = I.relative_commutant()
    print(f"   relative commutant N' c M: blocks {R.N.k_factors}, dim {sum(k * k for k in R.N.k_factors)}, "
          f"|[N, N' c M]| = {I.commutator_error(R).max().item():.1e}")


if __name__ == '__main__':
    main()
