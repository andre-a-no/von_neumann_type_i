[English](../01_matrix_algebras.md) · [Español](../es/01_matrix_algebras.md) · [Français](../fr/01_matrix_algebras.md) · [Deutsch](../de/01_matrix_algebras.md) · **中文** · [日本語](../ja/01_matrix_algebras.md) · [Русский](../ru/01_matrix_algebras.md)

# 1. 矩阵代数与函数演算

本章固定记号，并证明后续各章都要用到的关于 $M_n(\mathbb C)$ 的事实：谱定理、函数演算、
极分解以及由奇异值定义的范数。每个结果之后都给出相应的库代码和数值检验。

**你将学到**如何构造具有指定谱的算子，如何把一个函数作用于算子，算子的批量是如何存储的，
以及如何用精确结果检验一个计算。

**约定。** 内积 $\langle x,y\rangle=\sum_i\bar x_iy_i$ 关于第二个变元是线性的。$A^*$ 是共轭
转置，因此 $\langle x,Ay\rangle=\langle A^*x,y\rangle$。$e_1,\dots,e_n$ 是标准基，$1$ 是单位
矩阵，$\operatorname{Tr}$ 是通常的迹。

## 1.1 代数 $M_n(\mathbb C)$

$M_n(\mathbb C)$ 是复 $n\times n$ 矩阵构成的代数，带有对合
$A\mapsto A^*$ 和算子范数 $\|A\|=\sup_{\|x\|=1}\|Ax\|$。称元素 $A$ 为

| 名称 | 定义 |
|---|---|
| 正规的 | $A^*A=AA^*$ |
| 自伴的 | $A=A^*$ |
| 正的，$A\ge0$ | 对所有 $x$，$\langle x,Ax\rangle\ge0$ |
| 酉的 | $A^*A=AA^*=1$ |
| 投影 | $A=A^*=A^2$ |

若 $B-A\ge0$，记 $A\le B$；若 $A\ge0$ 且 $A$ 可逆，记 $A>0$。

**引理 1.1（C\*-恒等式）。** $\|A^*\|=\|A\|$ 且 $\|A^*A\|=\|A\|^2$。

*证明。* $\|Ax\|=\sup_{\|y\|=1}|\langle y,Ax\rangle|=\sup_{\|y\|=1}|\langle A^*y,x\rangle|$，
所以 $\|A\|=\sup_{\|x\|=\|y\|=1}|\langle A^*y,x\rangle|=\|A^*\|$。其次，当 $\|x\|=1$ 时
$\|Ax\|^2=\langle x,A^*Ax\rangle\le\|A^*A\|$，因此 $\|A\|^2\le\|A^*A\|\le\|A^*\|\,\|A\|=\|A\|^2$。$\square$

> **物理意义**（物理解释，而非本课程的定理）。在量子力学中，具有 $n$ 个可区分构型的系统由 $\mathbb C^n$ 描述。自伴矩阵是
> *可观测量*（能量、自旋、格点上的位置）；它们的特征值是测量的可能结果。投影对应于是/否
> 问题。酉矩阵是对称变换和时间演化，对哈密顿量 $H$ 而言即 $e^{-itH}$。迹为一的正矩阵是
> *态*（密度矩阵，第 3 章）。

## 1.2 谱定理

**引理 1.2（Schur）。** 对每个 $A\in M_n(\mathbb C)$，存在酉矩阵 $U$ 使得
$U^*AU$ 为上三角矩阵。

*证明。* 对 $n$ 归纳；$n=1$ 是平凡的。$A$ 的特征多项式有一个根 $\lambda$（代数基本定理），
所以对某个单位向量 $v$ 有 $Av=\lambda v$。将 $v$ 扩充为规范正交基 $v,w_2,\dots,w_n$，并令
$V$ 为以这些向量为列的酉矩阵。则
$V^*AV=\begin{pmatrix}\lambda&*\\0&A'\end{pmatrix}$，其中 $A'\in M_{n-1}(\mathbb C)$。
由归纳假设，存在酉矩阵 $U'$ 使 $U'^*A'U'$ 为上三角矩阵，于是
$U=V\,(1\oplus U')$ 满足要求。$\square$

**定理 1.3（谱定理）。** $A$ 是正规的当且仅当对某个酉矩阵 $U$ 有
$A=U\,\mathrm{diag}(\lambda_1,\dots,\lambda_n)\,U^*$。这些
$\lambda_i$ 是 $A$ 的计重数的特征值。

*证明。* 若 $A=U\Lambda U^*$，$\Lambda$ 为对角矩阵，则
$A^*A=U\bar\Lambda\Lambda U^*=U\Lambda\bar\Lambda U^*=AA^*$。反之，设 $A$
正规，$T=U^*AU$ 为上三角矩阵（引理 1.2）；$T$ 是正规的。比较
$(1,1)$ 元：$(T^*T)_{11}=|t_{11}|^2$，$(TT^*)_{11}=\sum_j|t_{1j}|^2$，所以当 $j>1$ 时
$t_{1j}=0$。于是 $T=t_{11}\oplus T'$，其中 $T'$ 正规且为上三角矩阵，
由归纳法知 $T$ 是对角矩阵。$U^*AU$ 的对角线由
$\det(A-\lambda)=\det(T-\lambda)$ 的根组成。$\square$

**推论 1.4.** 设 $A$ 正规，特征值为 $\lambda_i$。则

1. $\|A\|=\max_i|\lambda_i|$；
2. $A$ 自伴当且仅当所有 $\lambda_i\in\mathbb R$；$A$ 为酉矩阵当且仅当所有 $|\lambda_i|=1$；
   $A$ 为投影当且仅当所有 $\lambda_i\in\{0,1\}$；
3. $A\ge0$ 当且仅当 $A$ 自伴且所有 $\lambda_i\ge0$。此外，每个正矩阵都是自伴的，因而是正规的，所以 (3)
   在*所有*矩阵中刻画了正性，并且下面的函数演算适用于每个 $A\ge0$。

*证明。* 令 $y=U^*x$，则 $\|Ax\|^2=\sum_i|\lambda_i|^2|y_i|^2$ 且 $\|y\|=\|x\|$，
由此得 (1)，当 $x=Ue_i$ 且 $|\lambda_i|$ 最大时等号成立。在 (2) 中，$A^*=U\bar\Lambda U^*$，
$A^*A=U|\Lambda|^2U^*$，$A^2=U\Lambda^2U^*$，这些条件成立当且仅当它们对 $\Lambda$
逐元成立。对于 (3)，首先，每个正矩阵 $A$ 都是自伴的，因为对所有 $x$，
$\langle x,(A-A^*)x\rangle=\langle x,Ax\rangle-\overline{\langle x,Ax\rangle}=0$，
而对所有 $x$ 满足 $\langle x,Dx\rangle=0$ 的矩阵 $D$ 必为零（极化恒等式：
$4\langle x,Dy\rangle=\sum_{k=0}^3 i^{-k}\langle x+i^ky,D(x+i^ky)\rangle$）。其次，
$\langle x,Ax\rangle=\sum_i\lambda_i|y_i|^2$，它对所有 $x$ 都 $\ge0$ 当且仅当所有 $\lambda_i\ge0$。$\square$

## 1.3 函数演算

设 $A=U\Lambda U^*$ 正规，$f$ 是 $\operatorname{spec}A$ 上的复值函数。定义
$$
f(A)=U\,\mathrm{diag}\big(f(\lambda_1),\dots,f(\lambda_n)\big)\,U^*.
$$

**定理 1.5.** (a) $f(A)$ 与 $U$ 的选取无关。(b) 对多项式 $p$，$p(A)$ 就是通常意义下
$A$ 的多项式。(c) $f\mapsto f(A)$ 是线性且可乘的，$\bar f(A)=f(A)^*$，且 $\operatorname{spec}f(A)=f(\operatorname{spec}A)$。
(d) 若 $B$ 与 $A$ 可交换，则 $B$ 与 $f(A)$ 可交换。(e) 若 $g$ 定义在
$f(\operatorname{spec}A)$ 上，则 $g(f(A))=(g\circ f)(A)$。

*证明。* (b) $A^k=U\Lambda^kU^*$，由线性得 $p(A)=Up(\Lambda)U^*$。(a) 谱是有限集，
所以由 Lagrange 插值可得多项式 $p$，使得在 $\operatorname{spec}A$ 上 $p=f$：若 $\mu_1,\dots,\mu_r$ 是互不相同的特征值，取
$p(t)=\sum_if(\mu_i)\prod_{j\neq i}(t-\mu_j)/(\mu_i-\mu_j)$。于是对*每个*
对角化 $A$ 的酉矩阵 $U$ 都有 $U\,f(\Lambda)\,U^*=U\,p(\Lambda)\,U^*=p(A)$，而 $p(A)$ 的定义
不依赖于 $U$。(c) 对对角矩阵 $f(\Lambda)$ 成立，且在用 $U$ 共轭时保持；$Uf(\Lambda)U^*$
的谱就是 $f(\Lambda)$ 的对角线。(d) $B$ 与 $A$ 的每个多项式可交换，特别地与
$p(A)=f(A)$ 可交换。(e) $f(A)=Uf(\Lambda)U^*$ 是正规的，且由同一个 $U$ 对角化，所以由 (a)
得 $g(f(A))=U\,g(f(\Lambda))\,U^*=(g\circ f)(A)$。$\square$

**定理 1.6（平方根）。** 对 $A\ge0$，恰有一个 $B\ge0$ 满足
$B^2=A$，即由函数演算定义的 $B=\sqrt A$。

*证明。* 由定理 1.5(c) 和推论 1.4(3)，$\sqrt A\ge0$ 且 $(\sqrt A)^2=A$。
设 $B\ge0$，$B^2=A$，$B=V\,\mathrm{diag}(\mu_j)\,V^*$，$\mu_j\ge0$。则
$A=V\,\mathrm{diag}(\mu_j^2)\,V^*$，所以 $\operatorname{spec}A=\{\mu_j^2\}$。设 $p$ 是
在 $\operatorname{spec}A$ 上满足 $p(t)=\sqrt t$ 的多项式。则
$p(A)=\sqrt A$，且 $p(A)=p(B^2)=V\,\mathrm{diag}(p(\mu_j^2))\,V^*=V\,\mathrm{diag}(\mu_j)\,V^*=B$。$\square$

**推论 1.7.** $A\ge0$ 当且仅当对某个 $T$ 有 $A=T^*T$。（若 $A\ge0$，取 $T=\sqrt A$；
反之 $\langle x,T^*Tx\rangle=\|Tx\|^2\ge0$。）

## 1.4 极分解与奇异值

对任意 $A$，矩阵 $A^*A$ 是正的（推论 1.7）；它的平方根
$|A|=(A^*A)^{1/2}$ 称为 $A$ 的**绝对值**。

**定理 1.8（极分解）。** 每个 $A\in M_n(\mathbb C)$ 都可以写成
$A=W|A|$，其中 $W$ 为酉矩阵。

*证明。* 对每个 $x$，$\||A|x\|^2=\langle x,|A|^2x\rangle=\langle x,A^*Ax\rangle=\|Ax\|^2$。
因此 $W_0:|A|x\mapsto Ax$ 是从 $|A|$ 的值域 $R_1$ 到 $A$ 的值域 $R_2$ 上的良定义的线性
等距映射：若 $|A|x=|A|x'$，则
$\|A(x-x')\|=\||A|(x-x')\|=0$。特别地 $\dim R_1=\dim R_2$，所以二者的正交补也有相同
维数；令 $W_1:R_1^\perp\to R_2^\perp$ 为它们之间的任一等距映射（把一组规范正交基映为一组
规范正交基）。令
$W=W_0\oplus W_1$，即对 $u_1\in R_1$、$u_2\in R_1^\perp$ 有 $W(u_1+u_2)=W_0u_1+W_1u_2$。由于
$W_0u_1\in R_2\perp R_2^\perp\ni W_1u_2$，由勾股定理，对所有 $u$ 有 $\|Wu\|^2=\|u_1\|^2+\|u_2\|^2=\|u\|^2$。
因此对所有 $u$ 有 $\langle u,(W^*W-1)u\rangle=0$，由极化恒等式（推论 1.4(3) 的证明）得 $W^*W=1$；于是
$W$ 可逆且 $W^{-1}=W^*$，所以也有 $WW^*=1$：$W$ 是酉的，且 $W|A|x=W_0|A|x=Ax$。$\square$

$|A|$ 的特征值 $s_1\ge\dots\ge s_n\ge0$ 称为 $A$ 的**奇异值**。全书使用以下三种范数：
$$
\|A\|=s_1,\qquad \|A\|_2=\Big(\sum_is_i^2\Big)^{1/2}=(\operatorname{Tr}A^*A)^{1/2},
\qquad\|A\|_1=\sum_is_i=\operatorname{Tr}|A|.
$$
第一个等式成立是因为对每个 $x$ 有 $\|Ax\|=\||A|x\|$，所以由推论 1.4(1)，
$\|A\|=\||A|\|=s_1$；第二个等式成立是因为
$\operatorname{Tr}A^*A=\operatorname{Tr}|A|^2=\sum_is_i^2$。

**定理 1.9（迹范数的对偶性）。**
$\|A\|_1=\max_{U}|\operatorname{Tr}(AU)|$（最大值取遍酉矩阵 $U$），且最大值在
$U=W^*$ 处取得，其中 $A=W|A|$。特别地 $|\operatorname{Tr}A|\le\|A\|_1$。

*证明。* $\operatorname{Tr}(AW^*)=\operatorname{Tr}(W|A|W^*)=\operatorname{Tr}|A|$。对任意酉矩阵 $U$，
由迹的循环性以及关于 $\langle X,Y\rangle_{HS}=\operatorname{Tr}(X^*Y)$ 的 Cauchy–Schwarz 不等式，
$$
|\operatorname{Tr}(AU)|=|\operatorname{Tr}(|A|^{1/2}\,UW\,|A|^{1/2})|
=|\langle|A|^{1/2},\,UW|A|^{1/2}\rangle_{HS}|
\le\||A|^{1/2}\|_2\,\|UW|A|^{1/2}\|_2=\operatorname{Tr}|A|,
$$
这是因为对酉矩阵 $UW$ 有 $\|UW X\|_2=\|X\|_2$，且 $\||A|^{1/2}\|_2^2=\operatorname{Tr}|A|$。$\square$

特别地，上述三个量都是范数。$\|\cdot\|$ 是算子范数；
$\|\cdot\|_2$ 是内积 $\langle X,Y\rangle_{HS}$ 诱导的范数；对于
$\|\cdot\|_1$，由定理 1.9 得
$\|A+B\|_1=\max_U|\operatorname{Tr}(AU)+\operatorname{Tr}(BU)|\le\|A\|_1+\|B\|_1$ 且
$\|cA\|_1=\max_U|c\operatorname{Tr}(AU)|=|c|\,\|A\|_1$，而
$\|A\|_1=0$ 迫使 $|A|=0$，从而 $A=W|A|=0$。对 $P\ge0$，由定理 1.6 有 $|P|=(P^2)^{1/2}=P$，
所以 $\|P\|_1=\operatorname{Tr}P$。

## 1.5 库中的算子

本库处理代数 $\bigoplus_cM_{k_c}(\mathbb C)$（第 2 章）。单个块是 $C=1$ 的情形：
`TypeIAlgebra([n], [n])` 就是 $M_n(\mathbb C)$。（两个列表分别是名义块大小和实际块大小；
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
组成的列表。这里 `batch = C = 1`。这一构造恰好是把定理 1.3 反过来读：$A=U\Lambda U^*$，
其中 $U$ 服从 Haar 分布（库通过定理 5.3 的 QR 构造抽取 $U$）。标志 `force_positive=True` 会根据
指定的特征值加以检查（若其中某个小于 $-10^{-12}$，库会报错）；由推论 1.4(3)，所得结果就是正的
（当指定的特征值都 $\ge0$ 时严格如此，如此处；一般情形下在该容差范围内成立），库会记录这一点，
而不再次检验。出于同样的原因，上面的 `A.eigenvalues()` 返回排好序的指定谱，不做任何计算；
它并不是检验。下面对 $\operatorname{spec}f(A)=f(\operatorname{spec}A)$ 的检验则确实要计算：
`apply_function` 对角化 $A$ 的矩阵，所得矩阵 $e^A$ 的特征值是从头计算的。

自伴算子的函数演算由 `apply_function` 实现（它用 `eigh` 对角化每个块）。函数 $f$ 必须逐元素作用：
它把同一大小 $k$ 的所有扇区的特征值作为一个形状为 `(batch, n·k)` 的张量接收（$n$ 为这类扇区的个数），
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

## 1.6 范数与迹范数对偶

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
$\sim N(0,\tfrac12)$，所以 `G` 是 Ginibre 矩阵，它有密度（命题 A.16(b)）。
$|X|$ 几乎必然可逆：$(\det|X|)^2=\det(|X|^2)=\det(X^*X)=|\det X|^2$，而 $|\det X|^2$ 是各元素实部与虚部的
多项式，在 $X=1$ 处等于 $1$，所以它为零的概率为零（推论 A.7(b)；这也是定理 5.3 证明的第一步）。

## 1.7 批量：同时处理许多算子

批量轴正是本库存在的理由：随机矩阵理论、蒙特卡罗方法和多起点优化需要成千上万个算子，
而一次批量调用就取代了对它们的 Python 循环。作为例子，我们检验**高斯酉系综**（GUE）的一个精确矩。

设 $G$ 的各元素相互独立，且实部与虚部相互独立、服从
$N(0,\tfrac12)$，从而 $\mathbb E|G_{ij}|^2=1$，并令 $H=(G+G^*)/\sqrt2$。

**命题 1.10.** $\mathbb E\,\operatorname{Tr}H^4=2n^3+n$。因此对
$X=H/\sqrt n$，$\ \mathbb E\,\tfrac1n\operatorname{Tr}X^4=2+\tfrac1{n^2}$。

*证明。* *Isserlis（Wick）公式*（定理 A.24，$m=4$ 的情形）。$G$ 的各元素的实部与虚部是相互独立的
中心化高斯变量，所以它们是同一个中心化高斯向量的坐标，而 $H$ 与 $\bar H$ 的元素是它们的复线性组合。
于是定理 A.24 对任意四个这样的变量 $Z_1,\dots,Z_4$（允许重复）给出
$$
\mathbb E[Z_1Z_2Z_3Z_4]=\mathbb E[Z_1Z_2]\,\mathbb E[Z_3Z_4]+\mathbb E[Z_1Z_3]\,\mathbb E[Z_2Z_4]+\mathbb E[Z_1Z_4]\,\mathbb E[Z_2Z_3],
$$
即对 $\{1,2,3,4\}$ 的三种配对求和，每种配对贡献各对期望的乘积（不取复共轭）。附录 A 中的证明
把矩母函数的幂级数与高斯矩相比较，不需要在期望号下求导。

*矩的计算。* 配对期望为
$\mathbb E[H_{ij}H_{kl}]=\tfrac12\mathbb E[(G_{ij}+\bar G_{ji})(G_{kl}+\bar G_{lk})]=\delta_{il}\delta_{jk}$，
因为 $\mathbb E[G_{ij}G_{kl}]=0$ 且 $\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}$（命题 A.16(a)；
当 $(i,j)\neq(k,l)$ 时这两个元素相互独立且中心化，而对 $G_{ij}=x+iy$ 有
$\mathbb E[G_{ij}^2]=\mathbb E[x^2-y^2]+2i\,\mathbb E[xy]=0$ 和 $\mathbb E|G_{ij}|^2=\tfrac12+\tfrac12=1$）。
在 $\operatorname{Tr}H^4=\sum_{ijkl}H_{ij}H_{jk}H_{kl}H_{li}$ 中，各配对给出：
$\{12\}\{34\}$：$\delta_{ik}\delta_{ki}$，求和得 $n^3$（自由指标为 $i=k$、$j$、$l$）；
$\{14\}\{23\}$：$\delta_{jl}\delta_{jl}$，求和得 $n^3$；
$\{13\}\{24\}$：$\delta_{il}\delta_{jk}\delta_{ji}\delta_{kl}$，迫使 $i=j=k=l$，求和得 $n$。
除以 $n\cdot n^2$ 即得第二个结论。$\square$

两种*非交叉*配对 $\{12\}\{34\}$ 和 $\{14\}\{23\}$（这些对可以画成点 $1,2,3,4$ 上方互不相交的弧）
对每个 $n$ 都给出首项 $2$；交叉配对 $\{13\}\{24\}$ 给出 $1/n^2$ 修正项。

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
$\tfrac1n\operatorname{Tr}X^4$ 是高斯变量的多项式，所以它有有限的方差 $\sigma^2$（引理 A.23：所有高斯矩都有限；
乘积用 $|y_1\cdots y_r|\le\sum_i|y_i|^r$ 处理），由 Chebyshev 不等式，$B$ 个独立样本的均值偏离精确值
超过 $5\sigma/\sqrt B$ 的概率至多为 $1/25$。代码用样本标准差代替 $\sigma$，固定的随机种子使检验可以复现。
这是本教程中所有与已证明的结果对照的检验所遵循的模式：一个精确的命题、一个数值估计和一个明确的容差。

## 1.8 可与之比较的工具

对单个矩阵而言，`numpy.linalg` 和 `scipy.linalg`（`sqrtm`、`expm`、`logm`、
`polar`）是参考实现，上面的检验都可以用它们重复。本库增加了三项在后文中很重要的功能：
批量处理（1.7 节）、带有块结构的算子（第 2 章），以及贯穿其各项运算的自动微分（PyTorch 的 autograd）。

## 习题

1. 对随机自伴矩阵 $A$ 以及 $f=\sin$、$g=\exp$，检验 $(fg)(A)=f(A)g(A)$。
2. 利用定理 1.9 证明 $\|A\|_1=\|A^*\|_1$ 且 $\|AB\|_1\le\|A\|\,\|B\|_1$。
3. $A,B>0$ 的**几何平均**定义为 $A\#B=A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}$。
   证明 $A\#B$ 是方程 $XA^{-1}X=B$ 的唯一正解 $X$，由此推出
   $A\#B=B\#A$，并对两者进行数值检验。

<details>
<summary>解答</summary>

1. 由定理 1.5(c) 即得。

```python
Hs = M3.operator_from_eigenvalues(lambda d: torch.randn(1, d), force_self_adjoint=True)
assert close(Hs.apply_function(lambda w: torch.sin(w) * torch.exp(w)),
             Hs.apply_function(torch.sin) @ Hs.apply_function(torch.exp))
```

2. $|\operatorname{Tr}(A^*U)|=|\overline{\operatorname{Tr}(U^*A)}|=|\operatorname{Tr}(AU^*)|$，而
   $U\mapsto U^*$ 是酉群到自身的双射，所以定理 1.9 中的两个最大值相等。对第二个结论，设
   $B=V|B|$。对酉矩阵 $U$，
   $|\operatorname{Tr}(ABU)|=|\operatorname{Tr}(|B|^{1/2}UAV|B|^{1/2})|
   \le\||B|^{1/2}\|_2\,\|UAV|B|^{1/2}\|_2\le\|A\|\operatorname{Tr}|B|$，这里用到了
   $\|YX\|_2\le\|Y\|\,\|X\|_2$（对 $X$ 的每一列应用 $\|Yx\|\le\|Y\|\,\|x\|$）以及对酉矩阵
   $U,V$ 有 $\|UAV\|=\|A\|$（酉矩阵保持向量的长度）。

3. 记 $T=A^{-1/2}XA^{-1/2}$。则 $XA^{-1}X=B$ 当且仅当 $T^2=A^{-1/2}BA^{-1/2}$，并且
   $X\ge0$ 当且仅当 $T\ge0$，因为 $\langle x,Xx\rangle=\langle A^{1/2}x,TA^{1/2}x\rangle$ 且
   $A^{1/2}$ 可逆（这里 $A^{\pm1/2}$ 由 $\operatorname{spec}A\subset(0,\infty)$ 上的函数 $t^{\pm1/2}$ 通过函数演算
   给出，所以由定理 1.5(c) 有 $A^{-1/2}A^{1/2}=1$）。此外
   $A^{-1/2}BA^{-1/2}\ge0$，因为 $\langle x,A^{-1/2}BA^{-1/2}x\rangle=\langle A^{-1/2}x,BA^{-1/2}x\rangle\ge0$。
   由定理 1.6，$T=(A^{-1/2}BA^{-1/2})^{1/2}$ 是唯一的正解，即
   $X=A\#B$。方程 $XA^{-1}X=B$ 的每个解 $X$ 都可逆，
   因为 $\det(X)^2\det(A)^{-1}=\det B\neq0$。所以 $XA^{-1}X=B$ 等价于
   $XB^{-1}X=A$（由 $A^{-1}=X^{-1}BX^{-1}$ 两边取逆），再由唯一性得
   $A\#B=B\#A$。

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

**下一章：** [第 2 章](02_von_neumann_algebras.md)把 $M_n$ 换成矩阵代数的直和，并证明每个
有限维冯·诺依曼代数都具有这种形式。
