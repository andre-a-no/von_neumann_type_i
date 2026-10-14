#!/usr/bin/env python3
"""
Example: Jordan algebras of Hermitian matrices and spin factors (torch_vn_algebra.jordan).

The self-adjoint matrices form a Jordan algebra under x∘y = (xy + yx)/2. A set of Hermitian matrices generates a
Jordan algebra J and a complex *-algebra A; always J ⊆ A_sa (the self-adjoint part of A), and the two can differ:

  1. real symmetric matrices generate J = Herm_n(R), dim n(n+1)/2, but A = M_n(C), dim A_sa = n^2;
  2. the Pauli matrices span the spin factor V_3 = Herm_2(C) (J = A_sa), while five anticommuting symmetries on C^4
     span V_5, a 6-dimensional Jordan subalgebra of the 16-dimensional Herm_4(C);
  3. two Hermitian 2 x 2 matrices are simultaneously unitarily similar to real symmetric ones, so they generate only
     V_2 = Herm_2(R) as a Jordan algebra, although they generate M_2(C) as a *-algebra;
  4. the positive cone of V_m is the Lorentz cone a >= |u| = the PSD cone of its matrix representation; it is
     preserved by the quadratic representation P(y)x = yxy, but not by the Jordan product;
  5. the transpose preserves ∘ but reverses products: a Jordan automorphism, positive, not completely positive.
"""
import time

import torch

import torch_vn_algebra.decompose as dc
import torch_vn_algebra.jordan as jd
from torch_vn_algebra import TypeIAlgebra


def herm(n, gen, real=False):
    a = torch.randn(n, n, dtype=torch.float64 if real else torch.complex128, generator=gen)
    return a + a.mH


def main():
    t0 = time.time()
    gen = torch.Generator().manual_seed(0)

    print("1. Real symmetric generators (two generic ones)")
    print(f"   {'n':>3s} {'dim J':>6s} {'n(n+1)/2':>9s} {'dim A_sa':>9s} {'wedderburn dim A':>17s}")
    for n in (2, 3, 4, 6, 8):
        G = [herm(n, gen, real=True), herm(n, gen, real=True)]
        r = jd.compare_with_star_algebra(G)
        print(f"   {n:3d} {r['jordan_dim']:6d} {n * (n + 1) // 2:9d} {r['star_dim']:9d} {dc.wedderburn(G).dim:17d}")

    print("\n2. Spin factors from anticommuting symmetries s_i∘s_j = δ_ij 1")
    print(f"   {'m':>3s} {'D':>3s} {'dim J':>6s} {'dim A_sa':>9s}  J = A_sa")
    for m in range(1, 8):
        S = jd.anticommuting_symmetries(m)
        r = jd.compare_with_star_algebra(S)
        print(f"   {m:3d} {S.shape[-1]:3d} {r['jordan_dim']:6d} {r['star_dim']:9d}  {r['equal']}")
    V5 = jd.jordan_algebra_generated(jd.anticommuting_symmetries(5))
    print(f"   V_5 on C^4: {V5}, closed: {V5.is_closed()}, self-adjoint part of a *-algebra: "
          f"{V5.is_self_adjoint_part()}")

    print("\n3. The exceptional block size 2")
    for J in (2, 3):
        dims = [jd.jordan_algebra_generated([herm(k, gen) for _ in range(J)]).dim for k in (2, 3, 4)]
        print(f"   {J} generic Hermitian generators: dim J = {dims} for k = 2, 3, 4 (dim Herm_k(C) = [4, 9, 16])")

    print("\n4. The Lorentz cone of V_5 is the PSD cone of its representation on C^4")
    V = jd.SpinFactor(5)
    x = torch.randn(100000, 6, dtype=torch.float64, generator=gen)
    x[:, 0] *= 3
    X = V.to_matrix(x)
    w = torch.linalg.eigvalsh(X)
    lam = V.eigenvalues(x)
    print(f"   spectra: max |eig phi(x) - (λ-, λ-, λ+, λ+)| = {(w - lam.repeat_interleave(2, -1)).abs().max():.1e}")
    agree = (V.in_cone(x) == (w[:, 0] >= -1e-12)).float().mean()
    print(f"   cone membership a >= |u| agrees with phi(x) >= 0 on {100 * agree:.1f}% of {len(x)} samples "
          f"({100 * V.in_cone(x).float().mean():.1f}% in the cone)")
    p, q = V.element(1.0, torch.eye(5, dtype=torch.float64)[0]), V.element(1.0, torch.eye(5, dtype=torch.float64)[1])
    print(f"   p = (1, e_1), q = (1, e_2) in the cone; p∘q = {V.product(p, q).tolist()} in the cone: "
          f"{bool(V.in_cone(V.product(p, q)))}; P(p)q in the cone: {bool(V.in_cone(V.quadratic_representation(p, q)))}")
    y = V.random(4, gen, cone=True)
    err = (V.to_matrix(V.sqrt(y)) @ V.to_matrix(V.sqrt(y)) - V.to_matrix(y)).abs().max()
    err_exp = (V.to_matrix(V.exp(x[:4])) - torch.linalg.matrix_exp(X[:4])).abs().max()
    print(f"   functional calculus: |phi(sqrt y)^2 - phi(y)| = {err:.1e}; |phi(exp x) - expm(phi(x))| = {err_exp:.1e}")

    print("\n5. The transpose on M_3 (+) M_2: a Jordan automorphism that reverses products")
    alg = TypeIAlgebra([3, 2], [3, 2], precision='double')
    a, b = (alg.operator(torch.randn(1, 2, 3, 3, dtype=alg.dtype, generator=gen)) for _ in range(2))
    T = jd.transpose_map
    diff = lambda X, Y: float((X.matrix - Y.matrix).abs().max())
    print(f"   |T(a∘b) - T(a)∘T(b)| = {diff(T(jd.jordan_product(a, b)), jd.jordan_product(T(a), T(b))):.1e}, "
          f"|T(ab) - T(b)T(a)| = {diff(T(a @ b), T(b) @ T(a)):.1e}, "
          f"|T(ab) - T(a)T(b)| = {diff(T(a @ b), T(a) @ T(b)):.2f}")
    k = 3
    E = torch.zeros(k * k, 2, 3, 3, dtype=alg.dtype)
    for i in range(k):
        for j in range(k):
            E[i * k + j, 0, i, j] = 1
    TE = T(alg.operator(E)).matrix[:, 0].contiguous()
    choi = sum(torch.kron(E[i * k + j, 0], TE[i * k + j]) for i in range(k) for j in range(k))
    eig = sorted(set(torch.linalg.eigvalsh(choi).round(decimals=12).tolist()))
    print(f"   Choi matrix of the transpose on M_3: eigenvalues {eig} -> positive, not completely positive")
    print(f"\ntime {time.time() - t0:.2f} s")


if __name__ == '__main__':
    main()
