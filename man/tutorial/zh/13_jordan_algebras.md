[English](../13_jordan_algebras.md) · [Español](../es/13_jordan_algebras.md) · [Français](../fr/13_jordan_algebras.md) · [Deutsch](../de/13_jordan_algebras.md) · **中文** · [日本語](../ja/13_jordan_algebras.md) · [Русский](../ru/13_jordan_algebras.md)

# 教程 13. 若尔当代数：量子力学为何使用 *-代数

本教程配合讲义的[第 13 章](../../lecture_notes/zh/13_jordan_algebras.md)。该章研究若尔当积
$`x\circ y=\frac12(xy+yx)`$ 与二次表示 $`P(y)x=yxy`$、自旋因子 $`V_m`$ 与洛伦兹锥，
并说明若尔当结构只在相差符号的意义下决定交换子，而复合系统能察觉这一符号
（[命题 13.12](../../lecture_notes/zh/13_jordan_algebras.md#r-13-12)）。这里我们用模块 `torch_vn_algebra.jordan` 与类 `SpinFactor`
数值地检验这些结果。

<a id="s-13-1"></a>

## 13.1 可观测量与若尔当积

*讲义：* [13.1 节](../../lecture_notes/zh/13_jordan_algebras.md#s-13-1)：[引理 13.1](../../lecture_notes/zh/13_jordan_algebras.md#r-13-1)（若尔当恒等式）与[命题 13.2](../../lecture_notes/zh/13_jordan_algebras.md#r-13-2)（二次表示）。

模块 `torch_vn_algebra.jordan` 在 `TypeIAlgebra` 的算子上实现若尔当积 $`x\circ y`$、若尔当幂与二次表示
$`P(y)=2L(y)^2-L(y^2)`$，惰性求值，并在所有扇区上批量进行。
代码在 $`M_2\oplus M_3`$ 的随机自伴元上检验[引理 13.1](../../lecture_notes/zh/13_jordan_algebras.md#r-13-1)与[命题 13.2](../../lecture_notes/zh/13_jordan_algebras.md#r-13-2)：

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, SpinFactor
from torch_vn_algebra import jordan as jd

torch.manual_seed(0)
dt = torch.complex128
jp, P = jd.jordan_product, jd.quadratic_representation

def close(A, B, tol=1e-10):
    A = A.matrix if hasattr(A, 'matrix') else A
    B = B.matrix if hasattr(B, 'matrix') else B
    return (A - B).abs().max().item() < tol

alg = TypeIAlgebra([2, 3], [2, 3], precision='double')           # M_2 (+) M_3
def sa(B=4):
    return alg.operator_from_eigenvalues(lambda d: torch.randn(B, d), batch_size=B, force_self_adjoint=True)
x, y, z = sa(), sa(), sa()
x2 = jp(x, x)
assert close(jp(jp(x2, y), x), jp(x2, jp(y, x)))                  # Jordan identity (Lemma 13.1)
assert close(jd.jordan_power(x, 5), x @ x @ x @ x @ x)            # power associativity
assert close(P(y, x), 2 * jp(y, jp(y, x)) - jp(jp(y, y), x))      # P(y) = 2L(y)^2 - L(y^2)
assert close(P(P(y, x), z), P(y, P(x, P(y, z))))                  # fundamental formula
assert close(P(y, x).trace, (jp(y, y) @ x).trace)                 # Tr P(y)x = Tr y^2 x
pos = alg.operator_from_eigenvalues(lambda d: torch.rand(4, d), batch_size=4, force_positive=True)
assert all((e > -1e-12).all() for e in P(y, pos).eigenvalues())  # Proposition 13.2(3)

M2 = TypeIAlgebra([2], [2], precision='double')
p = M2.operator(torch.tensor([[[[1.0, 0], [0, 0]]]]), is_self_adjoint=True, is_positive=True)
q = M2.operator(torch.tensor([[[[0.5, 0.5], [0.5, 0.5]]]]), is_self_adjoint=True, is_positive=True)
print(jp(p, q).eigenvalues()[0])                                  # (1 -+ sqrt 2)/4: not positive
assert abs(torch.linalg.det(jp(p, q).matrix[0, 0]).real.item() + 1 / 16) < 1e-14
```

<a id="s-13-2"></a>

## 13.2 自旋因子与洛伦兹锥

*讲义：* [13.2 节](../../lecture_notes/zh/13_jordan_algebras.md#s-13-2)：[引理 13.3](../../lecture_notes/zh/13_jordan_algebras.md#r-13-3)（反对易的对称元）、[定理 13.4](../../lecture_notes/zh/13_jordan_algebras.md#r-13-4)（$`V_m`$ 的谱理论）与[命题 13.5](../../lecture_notes/zh/13_jordan_algebras.md#r-13-5)。

类 `SpinFactor(m)` 把 $`V_m`$ 的元素 $`x=(a,u)`$ 存储为形状为 `(batch, m+1)` 的张量，并提供乘积、
谱分解以及[定理 13.4](../../lecture_notes/zh/13_jordan_algebras.md#r-13-4)的函数演算；`to_matrix` 是[引理 13.3](../../lecture_notes/zh/13_jordan_algebras.md#r-13-3)的映射 $`\phi`$，
由 `jordan.anticommuting_symmetries(m)` 给出的 $`s_i`$ 构造。代码在 $`V_5`$ 中检验若尔当恒等式
与定理 13.4，对随机正交矩阵 $`O`$ 检验[命题 13.5](../../lecture_notes/zh/13_jordan_algebras.md#r-13-5)(b)，对 $`m\le6`$ 检验 (c)，并对 $`p=3`$ 检验
9.4 节的 Majorana 算子（[讲义，13.2 节](../../lecture_notes/zh/13_jordan_algebras.md#s-13-2)）：

```python
V = SpinFactor(5)
xs, ys = V.random(1000), V.random(1000)                          # (1000, 6): batches of (a, u)
x2 = V.square(xs)
assert (V.product(V.product(x2, ys), xs) - V.product(x2, V.product(ys, xs))).abs().max() < 1e-12
lam, e = V.spectral_decomposition(xs)                            # (lambda-, lambda+), (e-, e+)
assert (lam[:, 1, None] * e[:, 1] + lam[:, 0, None] * e[:, 0] - xs).abs().max() < 1e-12
Xm, Ym, Dm = V.to_matrix(xs), V.to_matrix(ys), V.representation_dim  # phi into Herm_4(C)
assert (V.to_matrix(V.product(xs, ys)) - 0.5 * (Xm @ Ym + Ym @ Xm)).abs().max() < 1e-12
assert (torch.linalg.eigvalsh(Xm) - lam.repeat_interleave(Dm // 2, dim=-1)).abs().max() < 1e-12
assert (V.to_matrix(V.exp(xs)) - torch.linalg.matrix_exp(Xm)).abs().max() < 1e-9       # phi(f(x)) = f(phi(x))
cone = V.random(1000, cone=True)                                 # squares
assert V.in_cone(cone).all() and (torch.linalg.eigvalsh(V.to_matrix(cone)) > -1e-12).all()
assert (V.square(V.sqrt(cone)) - cone).abs().max() < 1e-10
O, _ = torch.linalg.qr(torch.randn(5, 5, dtype=torch.float64))  # Proposition 13.5(b)
rot = lambda w: torch.cat([w[..., :1], w[..., 1:] @ O.T], dim=-1)
assert (rot(V.product(xs, ys)) - V.product(rot(xs), rot(ys))).abs().max() < 1e-12

for m in range(1, 7):                                            # Proposition 13.5(c)
    r = jd.compare_with_star_algebra(jd.anticommuting_symmetries(m))
    print(f"m = {m}: dim V_m = {r['jordan_dim']}, dim of the self-adjoint part of the *-algebra = {r['star_dim']}")
    assert r['jordan_dim'] == m + 1 and r['equal'] == (m in (1, 3))

pq = 3                                                           # Majorana operators of Section 9.4, p = 3 sites
Xp, Yp, Zp = jd.anticommuting_symmetries(3)
sminus, I2 = 0.5 * (Xp - 1j * Yp), torch.eye(2, dtype=dt)
parity = I2 - 2 * (0.5 * Zp + 0.5 * I2)                          # (-1)^n with n = s^z + 1/2
def kron_all(fs):
    out = torch.ones(1, 1, dtype=dt)
    for f_ in fs:
        out = torch.kron(out, f_)
    return out
s6 = jd.anticommuting_symmetries(2 * pq)
for j in range(1, pq + 1):
    c = kron_all([parity] * (j - 1) + [sminus] + [I2] * (pq - j))
    sign = (-1) ** (j - 1)
    assert close(s6[2 * j - 2], sign * (c + c.mH)) and close(s6[2 * j - 1], sign * 1j * (c - c.mH))
```

对 $`m=1`$，$`s_1=Z`$ 生成对角矩阵，维数为 $`2=m+1`$，而 $`V_1\cong\mathbb R\oplus\mathbb R`$ 不是
单的。打印出的维数是：$`m=2`$ 时为 $`3,4`$；$`m=4`$ 时为 $`5,16`$；$`m=6`$ 时为 $`7,64`$。

<a id="s-13-3"></a>

## 13.3 不是 *-子代数的若尔当子代数

*讲义：* [13.3 节](../../lecture_notes/zh/13_jordan_algebras.md#s-13-3)：[命题 13.6](../../lecture_notes/zh/13_jordan_algebras.md#r-13-6)。

`jordan.compare_with_star_algebra` 计算由给定厄米矩阵与 $`1`$ 生成的若尔当代数的维数，
以及它们生成的 $`*`$-代数的自伴部分的维数。代码对一对随机实对称矩阵（$`n=2,\dots,5`$）检验
[命题 13.6](../../lecture_notes/zh/13_jordan_algebras.md#r-13-6)(a)，并检验 (b)：对两个随机 $`2\times2`$ 厄米矩阵，维数分别为 $`3`$ 与 $`4`$，而证明中的酉矩阵 $`U`$ 使两者都变为实矩阵。

```python
for n in range(2, 6):                                            # Proposition 13.6(a), a random real pair
    G = torch.randn(2, n, n, dtype=torch.float64)
    r = jd.compare_with_star_algebra(G + G.mT)
    assert (r['jordan_dim'], r['star_dim']) == (n * (n + 1) // 2, n * n)

H = torch.randn(3, 2, 2, dtype=dt)
H = H + H.mH                                                     # three random Hermitian 2 x 2 matrices
r2, r3 = jd.compare_with_star_algebra(H[:2]), jd.compare_with_star_algebra(H)
print("two generators:", r2['jordan_dim'], r2['star_dim'], "  three:", r3['jordan_dim'], r3['star_dim'])
assert (r2['jordan_dim'], r2['star_dim'], r3['jordan_dim']) == (3, 4, 4)

V3 = SpinFactor(3)
u, v = V3.from_matrix(H[0])[1:], V3.from_matrix(H[1])[1:]        # Bloch vectors
nrm = torch.linalg.cross(u, v)
nrm = nrm / nrm.norm()
_, W1 = torch.linalg.eigh(V3.to_matrix(torch.cat([torch.zeros(1, dtype=torch.float64), nrm])))
_, W2 = torch.linalg.eigh(jd.anticommuting_symmetries(3)[1])     # Y
Um = W2 @ W1.mH                                                  # U (n . sigma) U* = Y
assert all((Um @ A @ Um.mH).imag.abs().max() < 1e-12 for A in H[:2])
```

三个一般的 $`2\times2`$ 厄米矩阵生成整个 $`\mathrm{Herm}_2(\mathbb C)`$：两列中的维数都是 $`4`$。
因此在 $`\mathcal M`$ 的一个 $`M_2`$ 块中，两个可观测量永远不会揭示复结构；第三个则
可以。

<a id="s-13-5"></a>

## 13.5 若尔当映射与 *-映射

*讲义：* [13.5 节](../../lecture_notes/zh/13_jordan_algebras.md#s-13-5)：[命题 13.9](../../lecture_notes/zh/13_jordan_algebras.md#r-13-9)（转置）、[命题 13.11](../../lecture_notes/zh/13_jordan_algebras.md#r-13-11)、[命题 13.12](../../lecture_notes/zh/13_jordan_algebras.md#r-13-12)（复合系统）与[命题 13.13](../../lecture_notes/zh/13_jordan_algebras.md#r-13-13)（动力学）。

`jordan.transpose_map` 是逐块作用的转置 $`T`$。代码检验[命题 13.9](../../lecture_notes/zh/13_jordan_algebras.md#r-13-9)（用 $`M_3`$ 上 $`T`$ 的
Choi 矩阵）、[命题 13.11](../../lecture_notes/zh/13_jordan_algebras.md#r-13-11)的证明中转置在 Bloch 向量上的矩阵 $`O_T=\mathrm{diag}(1,-1,1)`$、
[命题 13.12](../../lecture_notes/zh/13_jordan_algebras.md#r-13-12)的乘积公式与 $`\circ_\theta`$（连同 Bell 态的部分转置的谱），
以及[命题 13.13](../../lecture_notes/zh/13_jordan_algebras.md#r-13-13)(a)、(b)：

```python
x, y = sa(), sa()
T = jd.transpose_map
assert close(T(jp(x, y)), jp(T(x), T(y))) and close(T(x @ y), T(y) @ T(x))   # Proposition 13.9(a)
assert not close(T(x @ y), T(x) @ T(y))

k = 3
Mk = TypeIAlgebra([k], [k], precision='double')
units = torch.eye(k * k, dtype=dt).reshape(k * k, 1, k, k)                     # E_ij, batch index i*k + j
TE = T(Mk.operator(units)).matrix[:, 0]
choi = TE.reshape(k, k, k, k).permute(0, 2, 1, 3).reshape(k * k, k * k)       # block (i, j) = T(E_ij)
swap = torch.eye(k * k, dtype=dt).reshape(k, k, k, k).permute(0, 1, 3, 2).reshape(k * k, k * k)
w = torch.linalg.eigvalsh(choi)
assert torch.equal(choi, swap) and (w < 0).sum() == k * (k - 1) // 2          # Proposition 13.9(b)

sig = jd.anticommuting_symmetries(3)                                           # X, Y, Z
O_T = torch.stack([V3.from_matrix(s.T)[1:] for s in sig], dim=1)               # transpose on Bloch vectors
assert torch.allclose(O_T, torch.diag(torch.tensor([1.0, -1.0, 1.0], dtype=torch.float64)))

def ptranspose(R, n, m):                                                       # theta = id (x) T
    return R.reshape(n, m, n, m).permute(0, 3, 2, 1).reshape(n * m, n * m)
def herm(n):
    A = torch.randn(n, n, dtype=dt)
    return A + A.mH
jm = lambda A, B: 0.5 * (A @ B + B @ A)
cm = lambda A, B: A @ B - B @ A
a, a2, b, b2 = herm(2), herm(2), herm(3), herm(3)
local = torch.kron(jm(a, a2), jm(b, b2))
assert close(jm(torch.kron(a, b), torch.kron(a2, b2)), local + 0.25 * torch.kron(cm(a, a2), cm(b, b2)))
theta = lambda R: ptranspose(R, 2, 3)
assert close(theta(jm(theta(torch.kron(a, b)), theta(torch.kron(a2, b2)))),
             local - 0.25 * torch.kron(cm(a, a2), cm(b, b2)))                  # Proposition 13.12
Om = torch.tensor([1.0, 0, 0, 1.0], dtype=dt) / math.sqrt(2)
omega = torch.outer(Om, Om.conj())                                             # Bell state
print("spectrum of the partial transpose of the Bell state:", torch.linalg.eigvalsh(ptranspose(omega, 2, 2)))
assert abs(torch.linalg.eigvalsh(ptranspose(omega, 2, 2))[0] + 0.5) < 1e-12

Hm, z1, z2 = herm(3), herm(3), herm(3)
delta = lambda A: 1j * cm(Hm, A)
assert close(delta(jm(z1, z2)), jm(delta(z1), z2) + jm(z1, delta(z2)))         # Proposition 13.13(a)
assert close(jm(z1, jm(z2, Hm)) - jm(z2, jm(z1, Hm)), 0.25 * cm(cm(z1, z2), Hm))   # (b)
```

<a id="s-13-6"></a>

## 13.6 回到迹不等式

*讲义：* [13.6 节](../../lecture_notes/zh/13_jordan_algebras.md#s-13-6)：[命题 13.14](../../lecture_notes/zh/13_jordan_algebras.md#r-13-14)。

由[讲义，命题 13.14](../../lecture_notes/zh/13_jordan_algebras.md#r-13-14)，在 $`V_m`$ 中 $`\mathrm{tr}|P(y)x|\le\mathrm{tr}\,P(y)|x|`$，且当 $`x\in Q`$ 或
$`-x\in Q`$ 时等号成立。经由满足 $`\phi(P(y)x)=\phi(y)\phi(x)\phi(y)`$ 与 $`\phi(|x|)=|\phi(x)|`$ 的 $`\phi`$，它也就是
$`Y=Y^*`$ 时命题 11.9 的不等式 $`\mathrm{Tr}|YXY|\le\mathrm{Tr}(Y|X|Y)`$（[讲义，13.6 节](../../lecture_notes/zh/13_jordan_algebras.md#s-13-6)）。
代码对 $`m=2,4,7`$ 检验这两者：

```python
for m in (2, 4, 7):
    Vm = SpinFactor(m)
    xs, ys = Vm.random(20_000), Vm.random(20_000)
    lhs = Vm.trace(Vm.apply_function(Vm.quadratic_representation(ys, xs), torch.abs))
    rhs = Vm.trace(Vm.quadratic_representation(ys, Vm.apply_function(xs, torch.abs)))
    definite = Vm.determinant(xs) >= 0                                 # x or -x in the cone
    print(f"m = {m}: min gap {(rhs - lhs).min().item():.1e}, max gap {(rhs - lhs).max().item():.1f}, "
          f"definite x: {definite.double().mean().item():.3f}")
    assert (rhs - lhs >= -1e-12 * rhs).all()                            # Proposition 13.14
    assert ((rhs - lhs)[definite].abs() <= 1e-12 * rhs[definite]).all()
    Xs, Ys = Vm.to_matrix(xs[:200]), Vm.to_matrix(ys[:200])
    tn = torch.linalg.svdvals(Ys @ Xs @ Ys).sum(-1)                     # Tr |phi(y) phi(x) phi(y)|
    assert (tn - Vm.representation_dim / 2 * lhs[:200]).abs().max() < 1e-10 * tn.max()
    wx, Ux = torch.linalg.eigh(Xs)
    absX = Ux @ torch.diag_embed(wx.abs().to(dt)) @ Ux.mH                  # |phi(x)|
    tr_rhs = (Ys @ absX @ Ys).diagonal(dim1=-2, dim2=-1).sum(-1).real    # Tr phi(y) |phi(x)| phi(y)
    assert (torch.linalg.eigvalsh(Ys)[:, 0] < 0).any()                  # some phi(y) are indefinite
    assert (tr_rhs - tn >= -1e-10 * tr_rhs).all()                       # Proposition 11.9 for Y = Y*
```

差值 $`\mathrm{tr}\,P(y)|x|-\mathrm{tr}|P(y)x|`$ 对定号的 $`x`$ 为零（命题 13.14）。*数值观察*
（不是定理）：对这些 Gauss 样本，打印出的最大差值的量级为 $`10^1`$–$`10^2`$，所以对不定的 $`x`$，该不等式远非等式。

<a id="s-13-7"></a>

## 13.7 可与之比较的工具

据我们所知，目前没有主流的若尔当代数数值库；本章中的运算
通常是用矩阵手工写出的。计算机代数涵盖了代数方面：SageMath 的 `JordanAlgebra` 可以构造
特殊若尔当代数、对称双线性形式的若尔当代数（其中包括自旋因子），以及由一个
八元数代数构造的 27 维例外若尔当代数；Magma 也有若尔当代数的构造函数，
包括例外的那一个。这些是精确的代数工具，而不是批量化的数值工具。洛伦兹锥在
优化中是标准的：SOCP 求解器（ECOS、SCS、Clarabel、MOSEK，可通过 cvxpy 的 `SOC` 约束调用）使用 $`V_m`$ 的谱
分解（用于向锥投影）或它的二次表示（用于内点法的 Nesterov–Todd 缩放），
半定规划求解器对 $`\mathrm{Herm}_n`$ 也是如此。本库补充的是：若尔当积、二次
表示与转置，作为 $`\bigoplus_cM_{k_c}`$ 上批量化、可微的算子；带有
函数演算与矩阵表示的自旋因子；以及对给定矩阵所生成的若尔当代数的随机化计算，
并与它们所生成的 $`*`$-代数进行比较。

<a id="x-13"></a>

## 习题

讲义[习题](../../lecture_notes/zh/13_jordan_algebras.md#x-13) 1–3 的数值检验：

- *习题 1*（洛伦兹锥的齐性）：在 $`V_4`$ 中 $`P(x^{1/2})1=x`$、$`P(y^{-1})P(y)=\mathrm{id}`$ 以及 $`P(y)Q\subseteq Q`$。
- *习题 2*（各向同性态）：$`\theta(\rho_p)`$ 的最小特征值为 $`\frac{1-3p}4`$。
- *习题 3*（边界点的乘积）：对随机单位向量 $`w,w'`$ 有 $`(1,w)\circ(1,w')\notin Q`$，而
  $`(1,w)\circ(1,w)\in Q`$。

<details>
<summary>解答</summary>

*习题 1.*

```python
V4 = SpinFactor(4)
xi = V4.random(100, cone=True)                                   # interior points (almost surely)
r = V4.sqrt(xi)
assert (V4.quadratic_representation(r, V4.unit((100,))) - xi).abs().max() < 1e-10
yv = V4.random(100)
back = V4.quadratic_representation(V4.inverse(yv), V4.quadratic_representation(yv, xi))
assert (back - xi).abs().max() < 1e-8 * xi.abs().max()
assert V4.in_cone(V4.quadratic_representation(yv, xi)).all()
```

*习题 2.*

```python
for pp in (0.2, 1 / 3, 0.5):
    rho_p = pp * omega + (1 - pp) * torch.eye(4, dtype=dt) / 4
    lam_min = torch.linalg.eigvalsh(ptranspose(rho_p, 2, 2))[0].item()
    assert abs(lam_min - (1 - 3 * pp) / 4) < 1e-12
```

*习题 3.*

```python
V3c = SpinFactor(3)
wv = torch.nn.functional.normalize(torch.randn(2, 1000, 3, dtype=torch.float64), dim=-1)
prod = V3c.product(V3c.element(1.0, wv[0]), V3c.element(1.0, wv[1]))
assert not V3c.in_cone(prod).any()
assert V3c.in_cone(V3c.product(V3c.element(1.0, wv[0]), V3c.element(1.0, wv[0]))).all()
```

</details>

**下一章：** [教程，附录 A. 分析与概率中的若干事实](A_analysis.md)
