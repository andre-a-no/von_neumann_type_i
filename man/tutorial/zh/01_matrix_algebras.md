[English](../01_matrix_algebras.md) · [Español](../es/01_matrix_algebras.md) · [Français](../fr/01_matrix_algebras.md) · [Deutsch](../de/01_matrix_algebras.md) · **中文** · [日本語](../ja/01_matrix_algebras.md) · [Русский](../ru/01_matrix_algebras.md)

# 教程 1. 矩阵代数与函数演算

本教程配合讲义的[第 1 章](../../lecture_notes/zh/01_matrix_algebras.md)，该章证明了 $`M_n(\mathbb C)`$ 中的谱定理、函数演算、极分解
以及迹范数的对偶性。这里用库逐一检验这些结果。

**你将学到**如何构造具有指定谱的算子，如何把一个函数作用于算子，算子的批量是如何存储的，
以及如何用精确结果检验一个计算。

<a id="s-1-5"></a>

## 1.5 库中的算子

*讲义：* [1.2 节](../../lecture_notes/zh/01_matrix_algebras.md#s-1-2)和[1.3 节](../../lecture_notes/zh/01_matrix_algebras.md#s-1-3)：[定理 1.3](../../lecture_notes/zh/01_matrix_algebras.md#r-1-3)、[推论 1.4](../../lecture_notes/zh/01_matrix_algebras.md#r-1-4)、[定理 1.5](../../lecture_notes/zh/01_matrix_algebras.md#r-1-5)、[定理 1.6](../../lecture_notes/zh/01_matrix_algebras.md#r-1-6)、[定理 1.11](../../lecture_notes/zh/01_matrix_algebras.md#r-1-11)、[命题 1.12](../../lecture_notes/zh/01_matrix_algebras.md#r-1-12)、[命题 1.13](../../lecture_notes/zh/01_matrix_algebras.md#r-1-13)。

本库处理代数 $`\bigoplus_cM_{k_c}(\mathbb C)`$（第 2 章）。单个块是 $`C=1`$ 的情形：
`TypeIAlgebra([n], [n])` 就是 $`M_n(\mathbb C)`$。（两个列表分别是名义块大小和实际块大小；
目前二者相同。）

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra

torch.manual_seed(0)
M3 = TypeIAlgebra([3], [3], precision='double')       # the algebra M_3(C), complex128

# a positive operator with eigenvalues 1, 2, 4 and Haar-random eigenvectors
A = M3.operator_from_eigenvalues(lambda d: torch.tensor([[1.0, 2.0, 4.0]]),
                                 force_positive=True)
print(A.matrix.shape)          # (batch, C, k_max, k_max) = (1, 1, 3, 3)
print(A.eigenvalues()[0])      # one tensor per block, eigenvalues in ascending order
```

每个算子存储为一个形状为 `(batch, C, k_max, k_max)` 的张量：一批算子，每个算子是 `C` 个块
组成的列表。这里 `batch = C = 1`。这一构造恰好是把定理 1.3 反过来读：$`A=U\Lambda U^*`$，
其中 $`U`$ 服从 Haar 分布（库通过定理 5.3 的 QR 构造抽取 $`U`$）。标志 `force_positive=True` 会根据
指定的特征值加以检查（若其中某个小于 $`-10^{-12}`$，库会报错）；由推论 1.4(3)，所得结果就是正的
（当指定的特征值都 $`\ge0`$ 时严格如此，如此处；一般情形下在该容差范围内成立），库会记录这一点，
而不再次检验。出于同样的原因，上面的 `A.eigenvalues()` 返回排好序的指定谱，不做任何计算；
它并不是检验。下面对 $`\mathrm{spec}\,f(A)=f(\mathrm{spec}\,A)`$ 的检验则确实要计算：
`apply_function` 对角化 $`A`$ 的矩阵，所得矩阵 $`e^A`$ 的特征值是从头计算的。

自伴算子的函数演算由 `apply_function` 实现（它用 `eigh` 对角化每个块）。函数 $`f`$ 必须逐元素作用：
它把同一大小 $`k`$ 的所有扇区的特征值作为一个形状为 `(batch, n·k)` 的张量接收（$`n`$ 为这类扇区的个数），
并返回同样形状的张量。常用函数有各自专门的方法。我们来检验定理 1.5 和定理 1.6：

```python
def close(X, Y, tol=1e-10):
    """Maximal entry of |X - Y| for two operators (or tensors) is below tol."""
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

R = A.sqrt()
assert close(R @ R, A)                                  # Theorem 1.6
assert close(A.log().expm(), A)                         # exp(log A) = A, Theorem 1.5(e)
assert close(A.apply_function(lambda w: w ** 3), A @ A @ A)   # Theorem 1.5(b)
assert close(A.power(-1), A.inverse())                  # f(t) = 1/t
spec = A.apply_function(torch.exp).eigenvalues()[0]     # spec f(A) = f(spec A)
assert torch.allclose(spec, torch.exp(torch.tensor([[1.0, 2.0, 4.0]], dtype=spec.dtype)))
```

`R @ R` 此时还不进行任何计算：乘积、和以及函数只构建一个计算方案，在请求 `.matrix`
或诸如迹之类的数值时才进行求值。

**导数。** 由[定理 1.11](../../lecture_notes/zh/01_matrix_algebras.md#r-1-11)（Daleckii–Krein），$`f(A)`$ 沿自伴方向 $`E`$
的导数是 $`U(L\circ U^*EU)U^*`$，其中 $`L_{ij}=f^{[1]}(\lambda_i,\lambda_j)`$ 是差商构成的 Loewner 矩阵（本征值相等时取
$`f'(\lambda_i)`$）。`functional.divided_differences` 由本征值以及 $`f`$、$`f'`$ 的值计算 $`L`$，autograd 用它对 `apply_function`
求导。我们在简并谱处把公式与稠密指数函数的中心差商比较，后者的误差为 $`O(s^2)`$ 加上 $`\varepsilon/s`$ 量级的舍入：

```python
from torch_vn_algebra import functional

lam = torch.tensor([1.0, 1.0, 2.0], dtype=torch.float64)          # a degenerate spectrum
Q = torch.linalg.qr(torch.randn(3, 3, dtype=torch.complex128))[0]  # a unitary
A0 = Q @ torch.diag(lam).to(torch.complex128) @ Q.mH
E = torch.randn(3, 3, dtype=torch.complex128)
E = (E + E.mH) / 2                                                 # a self-adjoint direction
L = functional.divided_differences(lam, torch.exp(lam), torch.exp(lam))
D = Q @ (L * (Q.mH @ E @ Q)) @ Q.mH                                # Theorem 1.11 for f = exp
s = 1e-6
fd = (torch.linalg.matrix_exp(A0 + s * E) - torch.linalg.matrix_exp(A0 - s * E)) / (2 * s)
print(f"{(D - fd).abs().max().item():.1e}")
assert (D - fd).abs().max() < 1e-8
```

对距离为 $`h`$ 的相近本征值 $`a,b`$，商 $`\hat q`$ 因相消而损失精度，而中点值 $`\hat m=(f'(a)+f'(b))/2`$ 的误差至多为
$`T=h^2\max|f'''|/12`$ 加上舍入（[命题 1.12](../../lecture_notes/zh/01_matrix_algebras.md#r-1-12)）。`divided_differences` 当且仅当
$`|\hat q-\hat m|\le2R`$ 时取 $`\hat m`$，其中 $`R`$ 是[命题 1.13](../../lecture_notes/zh/01_matrix_algebras.md#r-1-13) 的舍入界，求值误差取
$`\varepsilon(|\hat f_a|+|a\hat f'_a|)`$；于是其误差至多为 $`\min(T,3R)`$。我们对 $`c=10^4`$ 附近的 $`f(x)=e^{x-c}`$ 检验这一点，
此时 $`f^{[1]}(a,b)=e^{a-c}(e^h-1)/h`$ 可用 `expm1` 精确计算。最后一个距离 $`h=0.06`$ 低于早期版本的阈值
$`\varepsilon^{1/3}\max(|a|,|b|)`$，早期版本在那里取的是中点值：

```python
eps = torch.finfo(torch.float64).eps
c = 1e4
f = lambda x: torch.exp(x - c)
for h in (1e-12, 1e-9, 1e-6, 1e-4, 1e-2, 0.06):
    w = torch.tensor([c, c + h], dtype=torch.float64)
    h = (w[1] - w[0]).item()                                       # the distance actually stored
    fw = f(w)                                                      # f' = f
    exact = math.expm1(h) / h                                      # f[a, b], as e^(a - c) = 1
    q = ((fw[1] - fw[0]) / h).item()                               # quotient
    m = (fw.sum() / 2).item()                                      # midpoint value
    y = functional.divided_differences(w, fw, fw)[0, 1].item()     # the library's choice
    R = eps * (fw.abs() + (w * fw).abs()).sum().item() / h + 4 * eps * abs(q)
    T = h ** 2 * math.exp(h) / 12 + 2 * eps                        # Proposition 1.12, plus rounding
    print(f"h = {h:.0e}: quotient {abs(q - exact):.0e}, midpoint {abs(m - exact):.0e}, library {abs(y - exact):.0e}")
    assert abs(y - exact) <= min(T, 3 * R) + 2 * eps              # Proposition 1.13(b)
```

*数值观察*：对 $`h\le10^{-6}`$，商的误差可达 $`5\cdot10^{-10}`$，在 $`h=0.06`$ 处为 $`10^{-15}`$；中点值对 $`h\le10^{-9}`$
在舍入意义下是精确的，在 $`h=0.06`$ 处误差为 $`3\cdot10^{-4}`$。除 $`h=10^{-4}`$ 外，库在每一行都取得两者中较小的误差；
在 $`h=10^{-4}`$ 处它取中点值（$`8\cdot10^{-10}`$），尽管商更好（$`10^{-12}`$）：$`R`$ 允许自变量有
$`\varepsilon|a|\approx2\cdot10^{-12}`$ 的扰动，这对 `torch.exp` 偏于保守，但对 $`e^{-i\tau x}`$ 这类要对 $`\tau x`$ 舍入的函数是必要的。
保证 $`\min(T,3R)`$ 在每一行都成立。

<a id="s-1-6"></a>

## 1.6 范数与迹范数对偶

*讲义：* [1.4 节](../../lecture_notes/zh/01_matrix_algebras.md#s-1-4)：[定理 1.8](../../lecture_notes/zh/01_matrix_algebras.md#r-1-8)、三种范数、[定理 1.9](../../lecture_notes/zh/01_matrix_algebras.md#r-1-9)。

一般的（非正规）矩阵通过 `alg.operator(tensor)` 引入。我们用 PyTorch 直接计算的奇异值来
检验三种范数，并检验极分解（定理 1.8）和对偶性（定理 1.9）：

```python
G = torch.randn(1, 1, 3, 3, dtype=torch.complex128)    # Ginibre matrix: no symmetry at all
X = M3.operator(G)
s = torch.linalg.svdvals(G[0, 0])                       # singular values, descending

assert abs(X.operator_norm().item() - s[0].item()) < 1e-12
assert abs(X.frobenius_norm().item() - s.norm().item()) < 1e-12
assert abs(X.trace_norm().item() - s.sum().item()) < 1e-12

absX = X.abs()
W = X @ absX.inverse()                                  # |X| is invertible almost surely, see below
assert close(W.adjoint() @ W, M3.identity())            # W is unitary
assert abs((X @ W.adjoint()).trace.item() - s.sum().item()) < 1e-10    # Tr(X W*) = ||X||_1
assert abs(X.trace.item()) <= X.trace_norm().item()     # |Tr X| <= ||X||_1
```

带复数 dtype 的 `torch.randn` 抽取相互独立的元素，其实部与虚部相互独立且
$`\sim N(0,\tfrac12)`$，所以 `G` 是 Ginibre 矩阵，它有密度（命题 A.16(b)）。
$`|X|`$ 几乎必然可逆：$`(\det|X|)^2=\det(|X|^2)=\det(X^*X)=|\det X|^2`$，而 $`|\det X|^2`$ 是各元素实部与虚部的
多项式，在 $`X=1`$ 处等于 $`1`$，所以它为零的概率为零（推论 A.7(b)；这也是定理 5.3 证明的第一步）。

<a id="s-1-7"></a>

## 1.7 批量：同时处理许多算子

*讲义：* [1.7 节](../../lecture_notes/zh/01_matrix_algebras.md#s-1-7)、[命题 1.10](../../lecture_notes/zh/01_matrix_algebras.md#r-1-10)。

批量轴正是本库存在的理由：随机矩阵理论、蒙特卡罗方法和多起点优化需要成千上万个算子，
而一次批量调用就取代了对它们的 Python 循环。作为例子，我们检验**高斯酉系综**（GUE）的一个精确矩。

由[讲义，命题 1.10](../../lecture_notes/zh/01_matrix_algebras.md#r-1-10)，若 $`G`$ 的各元素相互独立，其实部与虚部相互独立且
$`\sim N(0,\tfrac12)`$，并且 $`X=(G+G^*)/\sqrt{2n}`$，则 $`\mathbb E\,\tfrac1n\mathrm{Tr}\,X^4=2+\tfrac1{n^2}`$。

```python
n, B = 4, 200_000
G = torch.complex(torch.randn(B, 1, n, n), torch.randn(B, 1, n, n)).to(torch.complex128) / math.sqrt(2)
M4 = TypeIAlgebra([n], [n], precision='double')
X = M4.operator((G + G.conj().transpose(-2, -1)) / math.sqrt(2 * n), is_self_adjoint=True)

X2 = X @ X
m4 = (X2 @ X2).trace.real / n                           # shape (B,): one value per sample
mean, err = m4.mean().item(), m4.std().item() / math.sqrt(B)
print(f"E tr X^4 / n = {mean:.4f} +- {err:.4f}   (exact {2 + 1 / n**2:.4f})")
assert abs(mean - (2 + 1 / n**2)) < 5 * err
```

蒙特卡罗检验是一种统计检验：估计值必须在几个标准误差的范围内与精确值吻合。变量
$`\tfrac1n\mathrm{Tr}\,X^4`$ 是高斯变量的多项式，所以它有有限的方差 $`\sigma^2`$（引理 A.23：所有高斯矩都有限；
乘积用 $`|y_1\cdots y_r|\le\sum_i|y_i|^r`$ 处理），由 Chebyshev 不等式，$`B`$ 个独立样本的均值偏离精确值
超过 $`5\sigma/\sqrt B`$ 的概率至多为 $`1/25`$。代码用样本标准差代替 $`\sigma`$，固定的随机种子使检验可以复现。
这是本教程中所有与已证明的结果对照的检验所遵循的模式：一个精确的命题、一个数值估计和一个明确的容差。

<a id="s-1-8"></a>

## 1.8 可与之比较的工具

对单个矩阵而言，`numpy.linalg` 和 `scipy.linalg`（`sqrtm`、`expm`、`logm`、
`polar`）是参考实现，上面的检验都可以用它们重复。本库增加了三项在后文中很重要的功能：
批量处理（1.7 节）、带有块结构的算子（第 2 章），以及贯穿其各项运算的自动微分（PyTorch 的 autograd）。

<a id="x-1"></a>

## 习题

讲义的[习题 1 和习题 3](../../lecture_notes/zh/01_matrix_algebras.md#x-1)，数值检验：
- 对随机自伴矩阵 $`A`$ 以及 $`f=\sin`$、$`g=\exp`$，检验 $`(fg)(A)=f(A)g(A)`$（习题 1）。
- $`A,B>0`$ 的**几何平均** $`A\#B=A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}`$ 是方程 $`XA^{-1}X=B`$ 的唯一正解
  $`X`$（习题 3）。对这个方程和 $`A\#B=B\#A`$ 进行数值检验。

<details>
<summary>解答</summary>

```python
Hs = M3.operator_from_eigenvalues(lambda d: torch.randn(1, d), force_self_adjoint=True)
assert close(Hs.apply_function(lambda w: torch.sin(w) * torch.exp(w)),
             Hs.apply_function(torch.sin) @ Hs.apply_function(torch.exp))
```

```python
def geometric_mean(A, B):
    a, ai = A.sqrt(), A.power(-0.5)
    return a @ (ai @ B @ ai).sqrt() @ a

P = M3.operator_from_eigenvalues(lambda d: 0.5 + torch.rand(1, d), force_positive=True)
Q = M3.operator_from_eigenvalues(lambda d: 0.5 + torch.rand(1, d), force_positive=True)
GM = geometric_mean(P, Q)
assert close(GM, geometric_mean(Q, P), 1e-9)
assert close(GM @ P.inverse() @ GM, Q, 1e-9)
```

</details>

**下一章：** [教程 2](02_von_neumann_algebras.md)：库中的交换子、中心与代数 $`\bigoplus_cM_{k_c}(\mathbb C)`$。
