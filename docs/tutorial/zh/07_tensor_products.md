[English](../07_tensor_products.md) · [Español](../es/07_tensor_products.md) · [Français](../fr/07_tensor_products.md) · [Deutsch](../de/07_tensor_products.md) · **中文** · [日本語](../ja/07_tensor_products.md) · [Русский](../ru/07_tensor_products.md)

# 7. 张量积与纠缠

两个量子系统合在一起，由它们代数的张量积来描述。有扇区时存在两种积：普通的积，其中每个部分
各自保持自己的荷；以及*融合*积，其中只有总荷守恒。本章证明这两种积的结构，把偏迹定义为一个
嵌入的对偶，并研究纠缠。

**你将学到**扇区如何相乘以及本库如何为它们排序，为什么偏迹是一个信道、归一化后又是一个有限指标的
条件期望，荷守恒态的约化态是什么，纠缠如何分成数的部分与构型的部分，以及随机态的平均纯度与平均熵。

在本章中，$\mathcal M_1=\bigoplus_{c=1}^{C_1}M_{k_c}(\mathbb C)$ 作用在
$H_1=\bigoplus_c\mathbb C^{k_c}$ 上，极小中心投影为 $P_c$；
$\mathcal M_2=\bigoplus_{d=1}^{C_2}M_{l_d}(\mathbb C)$ 作用在 $H_2=\bigoplus_d\mathbb C^{l_d}$ 上，
极小中心投影为 $Q_d$。$\operatorname{Tr}$ 是算子所作用的空间上的迹；在 $\mathcal M_i$ 上它就是第 3 章的
blunt 迹。对 $\mathbb C^k\otimes\mathbb C^l$ 上的 $X$，记
$X_{(a,b),(a',b')}=\langle e_a\otimes e_b,X\,e_{a'}\otimes e_{b'}\rangle$，
并通过 $e_a\otimes e_b\mapsto e_{(a-1)l+b}$ 把 $\mathbb C^k\otimes\mathbb C^l$ 与 $\mathbb C^{kl}$ 等同
（Kronecker 约定），从而 $(A\otimes B)_{(a,b),(a',b')}=A_{aa'}B_{bb'}$。

## 7.1 两个代数的张量积

对 $H_1$ 上的 $*$-代数 $\mathcal M_1$ 和 $H_2$ 上的 $*$-代数 $\mathcal M_2$，**张量积**
$\mathcal M_1\otimes\mathcal M_2$ 是 $H_1\otimes H_2$ 上的算子 $A\otimes B$（$A\in\mathcal M_1$，$B\in\mathcal M_2$）
的线性张成。由于 $(A\otimes B)(A'\otimes B')=AA'\otimes BB'$，$(A\otimes B)^*=A^*\otimes B^*$ 且 $1\otimes1=1$，
它是一个 $*$-代数，因而是一个冯·诺依曼代数（定理 2.2）。

**定理 7.1.** (a) $H_1\otimes H_2=\bigoplus_{c,d}\mathbb C^{k_c}\otimes\mathbb C^{l_d}$，这是投影
$P_c\otimes Q_d$ 的值域的正交和，并且
$$
\mathcal M_1\otimes\mathcal M_2=\bigoplus_{c,d}M_{k_c}\otimes M_{l_d}\cong\bigoplus_{c,d}M_{k_cl_d}(\mathbb C):
$$
它由所有满足 $X=\sum_{c,d}(P_c\otimes Q_d)X(P_c\otimes Q_d)$ 的算子组成，在各个
$\mathbb C^{k_c}\otimes\mathbb C^{l_d}$ 上的块是任意的。(b) 它的极小中心投影是 $C_1C_2$ 个投影
$P_c\otimes Q_d$，并且 $\dim(\mathcal M_1\otimes\mathcal M_2)=\dim\mathcal M_1\dim\mathcal M_2$。
(c) 对 $\mathcal M_1$ 上的迹 $\operatorname{tr}_w$ 和 $\mathcal M_2$ 上的迹 $\operatorname{tr}_v$（定理 3.1），
$\operatorname{tr}_w(A)\operatorname{tr}_v(B)=\operatorname{tr}_{w\otimes v}(A\otimes B)$，其中
$(w\otimes v)_{cd}=w_cv_d$。特别地，积代数上的 `Tr_blunt`、`Tr_norm` 和 `tau_vN` 分别是各因子上同名迹的乘积。

*证明。* (a) 诸 $P_c\otimes Q_d$ 是两两正交的投影，其和为 $1\otimes1$，且
$P_c\otimes Q_d$ 的值域是 $\mathbb C^{k_c}\otimes\mathbb C^{l_d}$。由于 $P_c$、$Q_d$ 是中心元，
$(P_c\otimes Q_d)(A\otimes B)=AP_c\otimes BQ_d=(A\otimes B)(P_c\otimes Q_d)$，所以
$\mathcal M_1\otimes\mathcal M_2$ 的每个元素都是块对角的。反之，设 $E^c_{aa'}\in\mathcal M_1$ 和
$F^d_{bb'}\in\mathcal M_2$ 分别是块 $c$ 与块 $d$ 的矩阵单位，则在
$\mathbb C^{k_c}\otimes\mathbb C^{l_d}$ 内 $E^c_{aa'}\otimes F^d_{bb'}=|e_a\otimes e_b\rangle\langle e_{a'}\otimes e_{b'}|$；
这些正是 $M_{k_cl_d}(\mathbb C)$ 的矩阵单位，所以每个支撑在块 $(c,d)$ 中的算子都属于
$\mathcal M_1\otimes\mathcal M_2$。(b) 全矩阵代数直和的中心由到各直和项上的投影张成（第 2 章习题 1），且
$\sum_{c,d}k_c^2l_d^2=\sum_ck_c^2\sum_dl_d^2$。(c)
$\operatorname{Tr}(X\otimes Y)=\sum_{a,b}X_{aa}Y_{bb}=\operatorname{Tr}X\operatorname{Tr}Y$，所以
$\operatorname{tr}_{w\otimes v}(A\otimes B)=\sum_{c,d}w_cv_d\operatorname{Tr}A_c\operatorname{Tr}B_d=\operatorname{tr}_w(A)\operatorname{tr}_v(B)$。
三种迹的权重 $1$、$1/k_c$、$1/(C_1k_c)$（3.1 节）相乘后得到 $1$、
$1/(k_cl_d)$、$1/(C_1C_2k_cl_d)$，正是积代数上同名迹的权重。$\square$

`tensor_product(alg1, alg2)` 构造这个代数。扇区 $(c,d)$ 的编号为 $s=cC_2+d$（代码中编号
从 $0$ 开始），荷为 $(q_c,r_d)$，在扇区内部 $e_a\otimes e_b$ 的编号为 $a\,l_d+b$。
`kron(A, B, alg12)` 计算 $A\otimes B$。扇区顺序与 $H_1\otimes H_2$ 的 Kronecker 顺序相差一个基的置换，
下面的 `kron_rows` 显式地构造这个置换；借助它，我们与全矩阵的 Kronecker 积进行比较。

```python
import math
import numpy as np
import torch
from torch_vn_algebra import (TypeIAlgebra, DensityMatrix, states, composite,
                              tensor_product, fused_tensor_product, kron, partial_trace)

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def full(X):
    """An operator of (+)_c M_{k_c} as a block-diagonal matrix on (+)_c C^{k_c}, shape (batch, D, D)."""
    return torch.stack([torch.block_diag(*[x[c, :k, :k] for c, k in enumerate(X.algebra.k_factors)])
                        for x in X.matrix])

def kron_rows(alg12):
    """For every sector of a product algebra, the positions of its basis vectors in H1 (x) H2 (Kronecker order)."""
    k1, k2 = alg12.tensor_factors[0].k_factors, alg12.tensor_factors[1].k_factors
    return [torch.tensor([(sum(k1[:c]) + a) * sum(k2) + sum(k2[:d]) + b
                          for c, d, off in pairs for a in range(k1[c]) for b in range(k2[d])])
            for pairs in alg12.product_blocks]              # pair blocks (c, d, offset), by offset

def full12(X):
    """An operator of a (plain or fused) product algebra as a matrix on H1 (x) H2 in Kronecker order."""
    f1, f2 = X.algebra.tensor_factors
    N = sum(f1.k_factors) * sum(f2.k_factors)
    out = torch.zeros(X.matrix.shape[0], N, N, dtype=X.matrix.dtype)
    for s, r in enumerate(kron_rows(X.algebra)):
        out[:, r[:, None], r[None, :]] = X.matrix[:, s, :len(r), :len(r)]
    return out

def random_operator(alg):
    return alg.from_blocks([torch.randn(k, k, dtype=dt) for k in alg.k_factors])

alg1 = TypeIAlgebra([1, 2], [1, 2], precision='double')            # C (+) M_2
alg2 = TypeIAlgebra([2, 3], [2, 3], precision='double')            # M_2 (+) M_3
alg12 = tensor_product(alg1, alg2)
print(alg12.k_factors, alg12.charges, sum(k * k for k in alg12.k_factors))     # 4 sectors, dimension 65
assert sum(k * k for k in alg12.k_factors) == 5 * 13                           # Theorem 7.1(b)

A, A2, B, B2 = random_operator(alg1), random_operator(alg1), random_operator(alg2), random_operator(alg2)
AB = kron(A, B, alg12)
assert close(full12(AB)[0], torch.kron(full(A)[0], full(B)[0]))               # Theorem 7.1(a)
assert close(AB @ kron(A2, B2, alg12), kron(A @ A2, B @ B2, alg12))
for tr in ('Tr_blunt', 'Tr_norm', 'tau_vN'):                                   # Theorem 7.1(c)
    assert (getattr(AB, tr)() - getattr(A, tr)() * getattr(B, tr)()).abs().max() < 1e-10
```

> **物理意义。** 如果两个系统各自服从自己的超选择定则——两个不能交换粒子的盒子——那么联合系统
> 对每一*对*荷都有一个扇区，并且没有可观测量连接不同的荷对。一个联合扇区的权重是两个权重之积。

## 7.2 偏迹

**定义。** 设 $\rho\in\mathcal M_1\otimes\mathcal M_2$。**偏迹**
$\operatorname{Tr}_2\rho\in\mathcal M_1$ 是满足下式的元素：
$$
\operatorname{Tr}\big((\operatorname{Tr}_2\rho)\,A\big)=\operatorname{Tr}\big(\rho\,(A\otimes1)\big)\qquad\text{for all }A\in\mathcal M_1 .
$$
由定理 3.2（应用于线性泛函 $A\mapsto\operatorname{Tr}(\rho(A\otimes1))$ 和 $\mathcal M_1$ 的忠实迹
$\operatorname{Tr}$），它存在且唯一：偏迹是嵌入 $A\mapsto A\otimes1$ 关于 blunt 迹的对偶。
$\operatorname{Tr}_1$ 的定义与此相同。对全代数 $\mathcal M_i=B(H_i)$，这就是
**Hilbert 空间偏迹** $\operatorname{Tr}_{H_2}:B(H_1\otimes H_2)\to B(H_1)$。

**命题 7.2.** (a) 对 $\rho=\bigoplus_{c,d}\rho_{cd}$，
$$
(\operatorname{Tr}_2\rho)_c=\sum_d\operatorname{tr}_{l_d}\rho_{cd},\qquad
(\operatorname{tr}_lX)_{aa'}=\sum_{b=1}^lX_{(a,b),(a',b)}\quad(X\text{ on }\mathbb C^k\otimes\mathbb C^l),
$$
并且 $\operatorname{Tr}_2\rho=\operatorname{Tr}_{H_2}\rho$，其中 $(\operatorname{Tr}_{H_2}X)_{ii'}=\sum_jX_{(i,j),(i',j)}$。
(b) $\operatorname{Tr}_2(A\otimes B)=(\operatorname{Tr}B)\,A$。
(c) 对 $A,A'\in\mathcal M_1$，$\operatorname{Tr}_2\big((A\otimes1)X(A'\otimes1)\big)=A\,(\operatorname{Tr}_2X)\,A'$。
(d) $\operatorname{Tr}(\operatorname{Tr}_2X)=\operatorname{Tr}X$，并且 $\operatorname{Tr}_2$ 把态映为态。

*证明。* (a) 取 $A=\bigoplus_cA_c$，则 $A\otimes1=\bigoplus_{c,d}A_c\otimes1_{l_d}$，且
$$
\operatorname{Tr}\big(\rho(A\otimes1)\big)=\sum_{c,d}\sum_{a,a',b}(\rho_{cd})_{(a,b),(a',b)}(A_c)_{a'a}
=\sum_c\operatorname{Tr}\Big(\Big(\sum_d\operatorname{tr}_{l_d}\rho_{cd}\Big)A_c\Big),
$$
所以右边的元素满足定义，即为 $\operatorname{Tr}_2\rho$。对 $B(H_i)$ 作同样的计算即得
$\operatorname{Tr}_{H_2}$ 的公式；对块对角的 $\rho$，$i,i'$ 属于不同块的矩阵元
$\rho_{(i,j),(i',j)}$ 为零，所以 $\operatorname{Tr}_{H_2}\rho$ 具有 (a) 中的块。(b) $\operatorname{tr}_l(X\otimes Y)=(\operatorname{Tr}Y)X$。
(c) 对 $B\in\mathcal M_1$，
$\operatorname{Tr}\big(\operatorname{Tr}_2((A\otimes1)X(A'\otimes1))B\big)=\operatorname{Tr}\big(X(A'BA\otimes1)\big)
=\operatorname{Tr}\big(A(\operatorname{Tr}_2X)A'B\big)$；再用唯一性。(d) 取 $A=1$。若 $\rho$ 是态，则
$A\mapsto\operatorname{Tr}(\rho(A\otimes1))$ 是正的，因为 $A^*A\otimes1=(A\otimes1)^*(A\otimes1)$，并且它在
$A=1$ 处取值 $1$，所以它的密度 $\operatorname{Tr}_2\rho$ 是一个态（定理 3.2）。$\square$

信道是第 6 章的主题；这里只需要如下定义：线性映射 $\Phi$ 称为**完全正**的，如果对每个 $n$，
$\mathrm{id}_n\otimes\Phi$ 在 $M_n(\mathbb C)\otimes(\cdot)$ 上是正的；如果它还保迹，则称为**信道**。

**定理 7.3.** 设 $j$ 取遍 $H_2$ 的标准基，$K_j=1\otimes\langle j|:H_1\otimes H_2\to H_1$，
$K_j(x\otimes y)=\langle j,y\rangle x$。则
$$
\operatorname{Tr}_{H_2}X=\sum_jK_jXK_j^*,\qquad\sum_jK_j^*K_j=1,\qquad\sum_jK_j^*AK_j=A\otimes1 .
$$
因此 $\operatorname{Tr}_2$ 是从 $\mathcal M_1\otimes\mathcal M_2$ 到 $\mathcal M_1$ 的信道，并且它关于
blunt 迹的对偶是 $A\mapsto A\otimes1$。

*证明。* $K_j^*x=x\otimes e_j$，所以 $\langle e_a,K_jXK_j^*e_{a'}\rangle=X_{(a,j),(a',j)}$；对 $j$ 求和
即得命题 7.2(a)。此外 $\sum_jK_j^*K_j(x\otimes y)=\sum_j\langle j,y\rangle x\otimes e_j=x\otimes y$，
且 $\sum_jK_j^*AK_j=\sum_jA\otimes|j\rangle\langle j|=A\otimes1$。对 $Y\ge0$，
$(\mathrm{id}_n\otimes\operatorname{Tr}_{H_2})(Y)=\sum_j(1_n\otimes K_j)Y(1_n\otimes K_j)^*\ge0$，因为
$\langle x,TYT^*x\rangle=\langle T^*x,YT^*x\rangle\ge0$；这一性质在子代数
$\mathcal M_1\otimes\mathcal M_2$ 上依然成立。最后 $\operatorname{Tr}\sum_jK_jXK_j^*A=\operatorname{Tr}X\sum_jK_j^*AK_j$，
当 $A=1$ 时这就是保迹性，一般情形下则是对偶关系。$\square$

`partial_trace(X, keep=1)` 把每个块 $\rho_{cd}$ 重排为形状 $(k_c,l_d,k_c,l_d)$ 再求迹；
`composite.partial_trace_channel(alg12, keep=1)` 是定理 7.3 中的信道，其 Kraus 算子为从扇区 $(c,d)$
到扇区 $c$ 的 $1_{k_c}\otimes\langle j|$。我们把两者与定义、一个独立的 `numpy.einsum` 以及对偶关系进行核对：

```python
rho = states.random_density_matrix(alg12, batch_size=3)
rho1, rho2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
assert isinstance(rho1, DensityMatrix)                                      # states give states
R = full12(rho).numpy().reshape(3, 3, 5, 3, 5)                              # indices (batch, i, j, i', j')
assert np.abs(full(rho1).numpy() - np.einsum('xijkj->xik', R)).max() < 1e-12
assert np.abs(full(rho2).numpy() - np.einsum('xijil->xjl', R)).max() < 1e-12
one2 = alg2.identity()
assert (states.expectation(rho, kron(A, one2, alg12)) - states.expectation(rho1, A)).abs().max() < 1e-12
assert close(partial_trace(AB, keep=1), A * B.trace.item())                 # Proposition 7.2(b)

Ptr = composite.partial_trace_channel(alg12, keep=1)
print(Ptr)
assert Ptr.is_trace_preserving()
assert close(Ptr(rho), rho1)
assert close(Ptr.adjoint('blunt')(A), kron(A, one2, alg12))                 # Theorem 7.3: the dual
```

> **物理意义。** $\operatorname{Tr}_2\rho$ 是**约化态**：一个只能测量第一个系统的观察者所能得知的
> 全部信息——定义表明，它重现了每个可观测量 $A\otimes1$ 的期望值。约化态相同的联合态在关联上有所不同。

**作为条件期望的偏迹。** 设 $\mathcal M_2=M_m(\mathbb C)$，于是 $\mathcal M_1\otimes M_m$ 的扇区
就是 $\mathcal M_1$ 的扇区（定理 7.1）；令 $\iota(A)=A\otimes1_m$，
$\mathcal E=\frac1m\operatorname{Tr}_2:\mathcal M_1\otimes M_m\to\mathcal M_1$。设 $e_{ab}$ 为 $M_m$ 的矩阵单位。

**命题 7.4.** (a) $\iota\circ\mathcal E$ 是 $\mathcal M_1\otimes M_m$ 到 $\mathcal M_1\otimes1$ 上的条件期望
（2.6 节），并且对 $\mathcal M_1\otimes M_m$ 的每个忠实迹 $\operatorname{tr}_w$，它是唯一保持
$\operatorname{tr}_w$ 的条件期望。(b) $m^2$ 个元素 $u_{ab}=\sqrt m\,(1\otimes e_{ab})$ 对所有 $X$ 满足
$X=\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*X))$，且 $\sum_{a,b}u_{ab}u_{ab}^*=m^2\,1$。(c) 每个对所有 $X$
满足 $X=\sum_kv_k\,\iota(\mathcal E(v_k^*X))$ 的族 $(v_k)$ 都有 $\sum_kv_kv_k^*=m^2\,1$。

*证明。* (a) $\mathcal E$ 是正的（定理 7.3），$\iota$ 也是正的
（$T^*T\otimes1=(T\otimes1)^*(T\otimes1)$）；由命题 7.2(b)，$\mathcal E(A\otimes1)=A$，而双模性质
就是命题 7.2(c)。对 $\mathcal M_1$ 上具有相同权重 $w_c$ 的迹也记为 $\operatorname{tr}_w$；则
$\operatorname{tr}_w(Y\otimes1)=m\operatorname{tr}_w(Y)$。把 $\operatorname{Tr}_2$ 的定义应用于 $AP_c$，乘以 $w_c$ 并对 $c$ 求和，得到
$\operatorname{tr}_w(X(A\otimes1))=\operatorname{tr}_w((\operatorname{Tr}_2X)A)$，取 $A=1$ 得
$\operatorname{tr}_w(\iota\mathcal E(X))=\operatorname{tr}_w(X)$。若 $F=\iota\circ F_0$ 是到 $\mathcal M_1\otimes1$ 上另一个保持
$\operatorname{tr}_w$ 的条件期望，则
$\operatorname{tr}_w(X(A\otimes1))=\operatorname{tr}_w(F(X)(A\otimes1))=m\operatorname{tr}_w(F_0(X)A)$，所以对所有
$A\in\mathcal M_1$ 有 $\operatorname{tr}_w((\operatorname{Tr}_2X-mF_0(X))A)=0$，由忠实性得 $F_0=\mathcal E$
（取 $A=(\operatorname{Tr}_2X-mF_0(X))^*$）。(b) 写 $X=\sum_{a,b}X_{ab}\otimes e_{ab}$，其中
$X_{ab}\in\mathcal M_1$（定理 7.1(a)）。由命题 7.2(b)，
$\mathcal E(u_{ab}^*X)=\frac{\sqrt m}m\sum_{b'}\operatorname{Tr}_2(X_{ab'}\otimes e_{bb'})=X_{ab}/\sqrt m$，所以
$u_{ab}\,\iota(\mathcal E(u_{ab}^*X))=X_{ab}\otimes e_{ab}$；并且
$\sum_{a,b}u_{ab}u_{ab}^*=m\sum_{a,b}1\otimes e_{aa}=m^21$。(c) 由命题 7.2(a)，$\operatorname{Tr}_2(Y^*)=(\operatorname{Tr}_2Y)^*$，
所以 $Y=u_{ab}$ 的展开式取伴随后为
$u_{ab}^*=\sum_k\iota(\mathcal E(u_{ab}^*v_k))\,v_k^*$，再由 (b)，
$\sum_{a,b}u_{ab}u_{ab}^*=\sum_k\big(\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*v_k))\big)v_k^*=\sum_kv_kv_k^*$。$\square$

像 (b) 中那样的族称为 $\mathcal E$ 的一个*拟基*，而 $\sum_kv_kv_k^*$ 由 (c) 不依赖于拟基的选取，它就是
**Watatani 指标** $\operatorname{Ind}\mathcal E=m^2$（Y. Watatani, Mem. Amer. Math.
Soc. 83, no. 424, 1990）。对因子而言，$M_k\otimes1\subset M_k\otimes M_m$ 以重数 $m$ 在 $\mathbb C^{km}$ 中
包含 $\mathbb C^k$：其包含矩阵为 $\Lambda=(m)$，而
$m^2=\|\Lambda\|^2=\dim(M_k\otimes M_m)/\dim M_k$ 是它的 Jones 指标 $[\mathcal M:\mathcal N]$（V. F. R. Jones,
Invent. Math. 72, 1, 1983；F. M. Goodman, P. de la Harpe and V. F. R. Jones, *Coxeter graphs and
towers of algebras*, Springer, 1989）。对一般的 $\mathcal M_1$，这个包含是 $C_1$ 个包含
$M_{k_c}\otimes1\subset M_{k_c}\otimes M_m$ 的直和，$\Lambda=m1_{C_1}$。模块
`torch_vn_algebra.inclusions` 实现了包含、它们的条件期望和指标；
`from_tensor_product(alg1, m)` 构造上面这个包含：

```python
from torch_vn_algebra import inclusions
m = 3
inc = inclusions.from_tensor_product(alg1, m)                # alg1 (x) 1 c alg1 (x) M_3
X = random_operator(inc.M)
for trace in ('Tr_blunt', 'Tr_norm', 'tau_vN'):
    E = inc.conditional_expectation(trace.replace('Tr_', ''))
    assert close(E(X), partial_trace(X, keep=1) * (1 / m))   # Proposition 7.4(a): the same E for every trace
    assert abs(getattr(inc.embed(E(X)), trace)() - getattr(X, trace)()).max() < 1e-10
    assert close(E.index(), inc.M.identity() * m**2)         # Proposition 7.4(c)

def unit(a, b):
    e = torch.zeros(m, m, dtype=dt)
    e[a, b] = 1
    return math.sqrt(m) * kron(alg1.identity(), inc.M.tensor_factors[1].from_blocks([e]), inc.M)

u = [unit(a, b) for a in range(m) for b in range(m)]         # the quasi-basis of Proposition 7.4(b)
assert close(sum((v @ inc.embed(E(v.adjoint() @ X))).matrix for v in u), X)
assert close(sum((v @ v.adjoint()).matrix for v in u), inc.M.identity() * m**2)
print(inc.matrix.tolist(), inc.component_indices())           # Lambda = 3 * 1_2, index 9 per component
```

## 7.3 荷的融合

现在设 $\mathcal M_1$ 的扇区 $c$ 带有荷 $q_c$，$\mathcal M_2$ 的扇区 $d$ 带有荷 $r_d$，并且一个
**融合规则** $f$ 给出总荷 $f(q_c,r_d)$——对粒子数是加法，对费米子宇称是模 $2$ 加法，也可以是任何
其他映射。令
$$
H_q=\bigoplus_{f(q_c,r_d)=q}\mathbb C^{k_c}\otimes\mathbb C^{l_d},\qquad
\Pi_q=\sum_{f(q_c,r_d)=q}P_c\otimes Q_d,\qquad n_q=\sum_{f(q_c,r_d)=q}k_cl_d .
$$
**融合积** $\mathcal M_1\otimes_f\mathcal M_2=\{\Pi_q\}'$ 由 $H_1\otimes H_2$ 上保持总荷守恒的算子组成。

**命题 7.5.** (a) $\mathcal M_1\otimes_f\mathcal M_2=\bigoplus_qB(H_q)\cong\bigoplus_qM_{n_q}(\mathbb C)$，
极小中心投影为 $\Pi_q$。(b) $\mathcal M_1\otimes\mathcal M_2\subseteq\mathcal M_1\otimes_f\mathcal M_2$，
等号成立当且仅当每个 $H_q$ 只包含一个对 $(c,d)$。(c) 若荷是实数且 $f$ 是加法，则
$\mathcal M_1\otimes_f\mathcal M_2=\{N\}'$，其中
$N=N_1\otimes1+1\otimes N_2$，$N_1=\sum_cq_cP_c$，$N_2=\sum_dr_dQ_d$。

*证明。* (a) 每个对恰好属于一个 $H_q$，所以诸 $\Pi_q$ 是两两正交的投影，其和为 $1$，而 $X$ 与它们
全都可交换当且仅当 $X=\sum_q\Pi_qX\Pi_q$，其中 $B(H_q)$ 中的块是任意的；中心同定理 7.1(b)。(b)
$\mathcal M_1\otimes\mathcal M_2$ 的每个元素都与每个 $P_c\otimes Q_d$ 可交换，因而与它们的和 $\Pi_q$ 可交换。
若每个 $H_q$ 都只是一个 $\mathbb C^{k_c}\otimes\mathbb C^{l_d}$，则两个代数都是
$\bigoplus_{c,d}B(\mathbb C^{k_c}\otimes\mathbb C^{l_d})$；若 $H_q$ 包含两个对，则 $B(H_q)$ 包含把一个对的块
映到另一个对的块的算子，它们不是块对角的。(c)
$N=\sum_{c,d}(q_c+r_d)P_c\otimes Q_d=\sum_qq\,\Pi_q$，所以非零的 $\Pi_q$ 是 $N$ 的谱投影。它们是 $N$ 的
多项式（定理 1.5(a) 的证明），而 $N$ 是它们的线性组合，所以 $X$ 与 $N$ 可交换当且仅当它与每个
$\Pi_q$ 可交换。$\square$

`fused_tensor_product(alg1, alg2, fuse)`（默认的 `fuse` 为加法）按融合后的荷对扇区排序，扇区内部的
对 $(c,d)$ 按字典序排列，每个对都采用基 $e_a\otimes e_b$。
例子：硬核玻色子。两个格点上保持粒子数的算子构成
$M_1\oplus M_2\oplus M_1$（荷 $0,1,2$），三个格点上则构成 $M_1\oplus M_3\oplus M_3\oplus M_1$；它们的融合积
在基的顺序意义下就是五个格点的代数：$n_q=\sum_c\binom2c\binom3{q-c}=\binom5q$。

```python
bos1 = TypeIAlgebra([1, 2, 1], [1, 2, 1], precision='double', charges=[0, 1, 2])
bos2 = TypeIAlgebra([1, 3, 3, 1], [1, 3, 3, 1], precision='double', charges=[0, 1, 2, 3])
F = fused_tensor_product(bos1, bos2)                         # U(1): total particle number
P = tensor_product(bos1, bos2)                               # particle numbers conserved separately
print(F.charges, F.k_factors, sum(k * k for k in F.k_factors), sum(k * k for k in P.k_factors))
assert F.k_factors == [math.comb(5, q) for q in range(6)]

N1, N2 = full(bos1.central([0.0, 1.0, 2.0]))[0], full(bos2.central([0.0, 1.0, 2.0, 3.0]))[0]
N = torch.kron(N1, torch.eye(8, dtype=dt)) + torch.kron(torch.eye(4, dtype=dt), N2)
X = full12(F.operator(torch.randn(1, F.C, F.k_max, F.k_max, dtype=dt)))[0]
assert close(X @ N, N @ X)                                   # Proposition 7.5(c)
A, B = random_operator(bos1), random_operator(bos2)
assert close(full12(kron(A, B, F)), full12(kron(A, B, P)))   # Proposition 7.5(b)
```

融合代数的维数是 $252$，普通积的维数是 $120$；两者之差由在两部分之间移动粒子的算子构成，例如跃迁项。

$\mathcal M_1\otimes_f\mathcal M_2$ 上的态可以在不同的对块之间具有相干性，因此命题 7.2 不适用。
在 $B(H_1)$ 上令 $E_1(X)=\sum_cP_cXP_c$（2.6 节"物理意义"中的条件期望）；它映到 $\mathcal M_1$ 上，
并且对 $A\in\mathcal M_1$ 有 $\operatorname{Tr}(E_1(X)A)=\operatorname{Tr}(XA)$，因为 $A=\sum_cP_cAP_c$。

**定理 7.6.** 设 $\rho$ 是 $\mathcal M_1\otimes_f\mathcal M_2$ 上的态 $\omega$ 的密度。它到 $\mathcal M_1$ 的
限制 $A\mapsto\omega(A\otimes1)$ 的密度为
$$
\rho_1=E_1\big(\operatorname{Tr}_{H_2}\rho\big).
$$
$\operatorname{Tr}_{H_2}\rho$ 的块 $(c,c')$ 等于 $\sum_d\operatorname{tr}_{l_d}\rho_{(c,d),(c',d)}$，求和遍及满足
$f(q_c,r_d)=f(q_{c'},r_d)$ 的 $d$，其中 $\rho_{(c,d),(c',d)}$ 是 $\rho$ 从对 $(c',d)$ 到对 $(c,d)$ 的块。
因此，$E_1$ 对所有 $\operatorname{Tr}_{H_2}\rho$ 都平凡地作用，当且仅当不存在 $c\neq c'$ 和 $d$ 使得
$f(q_c,r_d)=f(q_{c'},r_d)$；例如当 $f$ 是加法且诸 $q_c$ 两两不同时即是如此。

*证明。* $A\otimes1\in\mathcal M_1\otimes_f\mathcal M_2$（命题 7.5(b)），且
$\omega(A\otimes1)=\operatorname{Tr}(\rho(A\otimes1))=\operatorname{Tr}\big((\operatorname{Tr}_{H_2}\rho)A\big)
=\operatorname{Tr}\big(E_1(\operatorname{Tr}_{H_2}\rho)A\big)$。由于 $E_1(\operatorname{Tr}_{H_2}\rho)\in\mathcal M_1$，
由定理 3.2 中的唯一性，它就是该密度。块公式是 $\operatorname{Tr}_{H_2}$ 的命题 7.2(a)，其中 $i$ 在块 $c$ 中，
$i'$ 在块 $c'$ 中，$j$ 在块 $d$ 中；除非两个对位于同一个 $H_q$ 中，否则 $\rho$ 在 $(c',d)$ 与 $(c,d)$ 之间的块为零。
若存在 $c\neq c'$ 和 $d$ 使得 $f(q_c,r_d)=f(q_{c'},r_d)=q$，取单位向量 $x\in\mathbb C^{k_c}$、
$x'\in\mathbb C^{k_{c'}}$、$y\in\mathbb C^{l_d}$，并令 $\psi=(x+x')\otimes y/\sqrt2\in H_q$。则
$|\psi\rangle\langle\psi|\in B(H_q)$，而 $\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|=\frac12|x+x'\rangle\langle x+x'|$
具有非零块 $\frac12|x\rangle\langle x'|$，它被 $E_1$ 去掉。$\square$

本库的两种实现对任意融合规则都计算这个映射：`partial_trace` 只对对角的对块 $\rho_{(c,d),(c,d)}$ 求迹，
而 `partial_trace_channel` 对同一扇区中 $c$ 相同的对使用各自独立的 Kraus 指标；它的对偶仍然是
$A\mapsto A\otimes1$。在我们的例子中，取加法时 $E_1$ 什么也不做；取宇称规则 $f(q,r)=q+r\bmod2$ 时，
`bos1` 的荷 $0$ 和 $2$ 与每个 $r$ 的融合结果都相同，于是 $E_1$ 去掉了真正的相干性：

```python
G = fused_tensor_product(bos1, bos2, fuse=lambda q, r: (q + r) % 2)       # only the parity is conserved
print(G.k_factors)

def block_mask(alg):
    """Ones in the diagonal blocks of (+)_c C^{k_c}: X * block_mask(alg) = sum_c P_c X P_c."""
    return full(alg.operator(torch.ones(1, alg.C, alg.k_max, alg.k_max, dtype=dt)))

for alg, coherent in ((F, False), (G, True)):
    rho = states.random_density_matrix(alg, batch_size=2)
    R = full12(rho).reshape(2, 4, 8, 4, 8)
    for keep, X_H, alg_k in ((1, torch.einsum('xijkj->xik', R), bos1), (2, torch.einsum('xijil->xjl', R), bos2)):
        E_X = X_H * block_mask(alg_k)                                     # E(Tr_H rho), Theorem 7.6
        assert close(full(partial_trace(rho, keep=keep)), E_X)
        assert close(full(composite.partial_trace_channel(alg, keep=keep)(rho)), E_X)
        lost = (X_H - E_X).abs().max().item()                             # coherences removed by E
        print(f"keep={keep}: largest coherence removed: {lost:.3e}")
        assert (lost > 1e-3) if coherent else (lost == 0)
    assert close(composite.partial_trace_channel(alg, keep=1).adjoint('blunt')(A), kron(A, bos2.identity(), alg))
```

> **物理意义。** 在超导体的 BCS 描述中，只有电子数的宇称守恒，因此整体态可以是一个区域的粒子数
> $0$ 与 $2$ 的叠加。一个其可观测量保持该区域粒子数守恒的观察者看不到这种相干性：对他而言，
> 约化态是定理 7.6 中经过块对角化（pinching）的那个态，即不同粒子数的混合。

## 7.4 纯态的纠缠

**引理 7.7（Schmidt 分解）。** 每个单位向量 $\psi\in\mathbb C^a\otimes\mathbb C^b$ 都可以写成
$\psi=\sum_{i=1}^rs_i\,u_i\otimes v_i$，其中 $s_i>0$，$\sum_is_i^2=1$，且 $u_i\in\mathbb C^a$、$v_i\in\mathbb C^b$
分别规范正交。于是 $\operatorname{Tr}_2|\psi\rangle\langle\psi|=\sum_is_i^2|u_i\rangle\langle u_i|$，
$\operatorname{Tr}_1|\psi\rangle\langle\psi|=\sum_is_i^2|v_i\rangle\langle v_i|$：两个约化态具有相同的非零特征值 $s_i^2$。

*证明。* 设 $\Psi$ 为 $a\times b$ 矩阵 $\Psi_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,\psi\rangle$。由谱定理，
$\Psi\Psi^*=\sum_{i=1}^a\lambda_i|u_i\rangle\langle u_i|$，其中 $u_i$ 是一组规范正交基，
$\lambda_i=\|\Psi^*u_i\|^2\ge0$，并排序使得 $\lambda_i>0$ 恰好对 $i\le r$ 成立。对 $i\le r$ 令
$s_i=\sqrt{\lambda_i}$，$w_i=\Psi^*u_i/s_i$；则
$\langle w_i,w_j\rangle=\langle u_i,\Psi\Psi^*u_j\rangle/(s_is_j)=\delta_{ij}$，且当 $i>r$ 时 $\Psi^*u_i=0$。
因此 $\Psi=\sum_i|u_i\rangle\langle u_i|\Psi=\sum_{i\le r}|u_i\rangle\langle\Psi^*u_i|=\sum_{i\le r}s_i|u_i\rangle\langle w_i|$，
即 $\Psi_{\alpha\beta}=\sum_is_i(u_i)_\alpha\overline{(w_i)_\beta}$：这就是 $v_i=\bar w_i$ 时的结论，并且
$\sum_is_i^2=\operatorname{Tr}\Psi\Psi^*=\|\psi\|^2=1$。最后
$|\psi\rangle\langle\psi|=\sum_{i,j}s_is_j|u_i\rangle\langle u_j|\otimes|v_i\rangle\langle v_j|$，
由命题 7.2(b) 即得约化态。$\square$

$\psi$ 的**纠缠熵**为
$E(\psi)=S(\operatorname{Tr}_2|\psi\rangle\langle\psi|)=S(\operatorname{Tr}_1|\psi\rangle\langle\psi|)=-\sum_is_i^2\log s_i^2$。
它为零当且仅当 $r=1$，即当且仅当 $\psi=u\otimes v$；并且 $E(\psi)\le\log\min(a,b)$，等号成立当且仅当
较小的那个约化态是最大混合态（推论 3.6）。

有荷时，纯态位于融合积的某一个扇区 $H_q$ 中；记
$\psi=\sum_{(c,d)}\psi_{cd}$，其中 $\psi_{cd}\in\mathbb C^{k_c}\otimes\mathbb C^{l_d}$。

**命题 7.8（对称性分辨纠缠）。** 设 $f$ 是加法，诸 $q_c$ 两两不同，诸 $r_d$ 也两两不同，从而对每个 $c$
至多有一个 $d=d(c)$ 满足 $q_c+r_d=q$。令 $\psi_c=\psi_{c\,d(c)}$（或 $0$），$p_c=\|\psi_c\|^2$，并在 $p_c>0$ 时令
$\hat\psi_c=\psi_c/\sqrt{p_c}$。则 $\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|$ 属于 $\mathcal M_1$（它与 $N_1$
可交换），并且等于 $\bigoplus_cp_c\operatorname{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|$；$\mathcal M_2$ 上的约化态在
$d(c)$ 处具有扇区概率 $p_c$；并且
$$
E(\psi)=H(p)+\sum_cp_c\,E(\hat\psi_c),\qquad H(p)=-\sum_cp_c\log p_c .
$$

*证明。* 由定理 7.6 的块公式，$\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|$ 的块 $(c,c')$
是 $\operatorname{tr}_{l_d}|\psi_{cd}\rangle\langle\psi_{c'd}|$ 对满足 $q_c+r_d=q=q_{c'}+r_d$ 的 $d$ 之和，
这迫使 $c=c'$ 且 $d=d(c)$。所以约化态为
$\bigoplus_c\operatorname{Tr}_2|\psi_c\rangle\langle\psi_c|$，其块的迹为 $p_c$；交换两个因子的角色即得关于
$\mathcal M_2$ 的结论。由引理 7.7，$E(\psi)$ 是这个态的熵，再对
$\sigma_c=\operatorname{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|$、$S(\sigma_c)=E(\hat\psi_c)$ 应用定理 3.3 即得公式。$\square$

$H(p)$ 称为**数纠缠**，$\sum_cp_cE(\hat\psi_c)$ 称为**构型纠缠**。对五个格点上 $2$ 个粒子的随机态，
我们分别由两侧的偏迹以及 $4\times8$ 矩阵 $\Psi$ 的奇异值计算 $E(\psi)$：

```python
s = F.charges.index(2)                                      # the sector of two particles, n_q = 10
n = F.k_factors[s]
psi = torch.randn(n, dtype=dt)
psi = psi / psi.norm()
mat = torch.zeros(1, F.C, F.k_max, F.k_max, dtype=dt)
mat[0, s, :n, :n] = torch.outer(psi, psi.conj())
rho = DensityMatrix(F, matrix=mat)
r1, r2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)

psi_H = torch.zeros(4 * 8, dtype=dt)
psi_H[kron_rows(F)[s]] = psi                                # psi as a vector of H1 (x) H2
lam = torch.linalg.svdvals(psi_H.reshape(4, 8)) ** 2        # squared Schmidt coefficients
lam = lam[lam > 1e-14]
E_psi = -(lam * lam.log()).sum().item()
assert abs(r1.von_neumann_entropy().item() - E_psi) < 1e-10
assert abs(r2.von_neumann_entropy().item() - E_psi) < 1e-10  # Lemma 7.7

def entropy(sigma):
    w = torch.linalg.eigvalsh(sigma).clamp(min=1e-300)
    return -(w * w.log()).sum().item()

p = r1.sector_probabilities()[0]                            # c = 0, 1, 2 particles in the first part
assert torch.allclose(r2.sector_probabilities()[0, [2, 1, 0]], p)   # d(c) = 2 - c
H_num = -(p * p.log()).sum().item()
S_conf = sum(p[c].item() * entropy(r1.matrix[0, c, :k, :k] / p[c]) for c, k in enumerate(bos1.k_factors))
print(f"E = {E_psi:.4f} = H(p) {H_num:.4f} + configurational {S_conf:.4f}")
assert abs(E_psi - H_num - S_conf) < 1e-10                  # Proposition 7.8
```

> **物理意义。** 受粒子数超选择定则约束的各方只能把构型部分用作资源：这就是 H. M. Wiseman 和 J. A.
> Vaccaro 的*粒子纠缠*（Phys. Rev. Lett. 91, 097902, 2003）。$H(p)$ 在冷原子实验中作为*数熵*
> 被测量（A. Lukin et al., Science 364, 256, 2019）。

## 7.5 随机态与 Page 公式

设 $g\in\mathbb C^N$ 的各坐标相互独立，其实部与虚部是独立的 $N(0,\tfrac12)$ 随机变量，于是
$\mathbb E|g_i|^2=1$。**Haar 随机纯态**是 $\psi=g/\|g\|$；它的分布是酉不变的，因为密度
$\pi^{-N}e^{-\|g\|^2}$ 与 Lebesgue 测度都是酉不变的（它是球面上唯一的不变分布，这一事实我们不会用到）。

**引理 7.9.** $\|g\|$ 与 $\psi=g/\|g\|$ 相互独立。

*证明。* 设 $h\ge0$ 与 $\phi\ge0$ 分别是 $(0,\infty)$ 和单位球面 $S^{2N-1}$ 上的可测函数。
在 $\mathbb R^{2N}$ 中用极坐标积分（$\sigma$ 为面积测度）得
$$
\mathbb E\big[h(\|g\|)\phi(\psi)\big]=\pi^{-N}\int_0^\infty h(t)e^{-t^2}t^{2N-1}dt\int_{S^{2N-1}}\phi\,d\sigma=:\alpha(h)\beta(\phi).
$$
先取 $\phi=1$，再取 $h=1$，再取两者都为 $1$，得到 $\mathbb Eh(\|g\|)=\alpha(h)\beta(1)$、
$\mathbb E\phi(\psi)=\alpha(1)\beta(\phi)$ 以及 $\alpha(1)\beta(1)=1$，因此
$\mathbb E[h(\|g\|)\phi(\psi)]=\mathbb Eh(\|g\|)\,\mathbb E\phi(\psi)$。$\square$

**引理 7.10.** $\mathbb E\|g\|^4=N(N+1)$。若 $G$ 是一个 $a\times b$ 矩阵，其各元素独立且与 $g_i$ 同分布
（*Ginibre* 矩阵），则 $\mathbb E\operatorname{Tr}(GG^*)^2=ab(a+b)$。

*证明。* 记 $g_i=x+iy$，则 $\mathbb E|g_i|^4=\mathbb Ex^4+2\mathbb Ex^2\mathbb Ey^2+\mathbb Ey^4=\tfrac34+\tfrac12+\tfrac34=2$，
所以 $\mathbb E\|g\|^4=\sum_i\mathbb E|g_i|^4+\sum_{i\ne j}\mathbb E|g_i|^2\mathbb E|g_j|^2=2N+N(N-1)$。其次，
$\operatorname{Tr}(GG^*)^2=\sum_{i,i',j,j'}G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}$。诸元素及其共轭是独立实高斯变量的
复线性组合，且 $\mathbb E[G_{ij}G_{kl}]=0$，$\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}$，所以 Isserlis 公式
（命题 1.10 的证明）给出 $\mathbb E[G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}]=\delta_{ii'}+\delta_{jj'}$。
求和得 $ab^2+a^2b$。$\square$

**定理 7.11（E. Lubkin, J. Math. Phys. 19, 1028, 1978）。** 对 Haar 随机的
$\psi\in\mathbb C^a\otimes\mathbb C^b$ 和 $\rho_A=\operatorname{Tr}_2|\psi\rangle\langle\psi|$，
$$
\mathbb E\operatorname{Tr}\rho_A^2=\frac{a+b}{ab+1}.
$$

*证明。* 令 $G_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,g\rangle$，这是一个 Ginibre 矩阵；$\psi$ 的矩阵为
$\Psi=G/\|g\|$，由引理 7.7 的证明，$\rho_A=\Psi\Psi^*=GG^*/\|g\|^2$。因此
$\operatorname{Tr}(GG^*)^2=\|g\|^4\operatorname{Tr}\rho_A^2$，其中 $\operatorname{Tr}\rho_A^2\in[0,1]$ 是 $\psi$ 的函数。
由引理 7.9，$\mathbb E\operatorname{Tr}(GG^*)^2=\mathbb E\|g\|^4\,\mathbb E\operatorname{Tr}\rho_A^2$，再由引理 7.10，
$\mathbb E\operatorname{Tr}\rho_A^2=ab(a+b)/(ab(ab+1))$。$\square$

第 6 章由 Haar 向量的四阶矩得到同一个均值，作为随机信道的输出纯度。

**推论 7.12.** $\log\frac{ab+1}{a+b}\le\mathbb E\,S(\rho_A)\le\log\min(a,b)$。

*证明。* 由推论 3.6，$S(\rho_A)\le\log a$；由引理 7.7，
$S(\rho_A)=S(\operatorname{Tr}_1|\psi\rangle\langle\psi|)\le\log b$。对凹函数 $\log$ 以
$\rho_A$ 的非零特征值 $\lambda_i>0$ 为权重应用 Jensen 不等式，得
$\sum_i\lambda_i\log\lambda_i\le\log\sum_i\lambda_i^2$，即 $S(\rho_A)\ge-\log\operatorname{Tr}\rho_A^2$；对凸函数
$-\log$ 应用 Jensen 不等式并结合定理 7.11，得
$\mathbb E[-\log\operatorname{Tr}\rho_A^2]\ge-\log\mathbb E\operatorname{Tr}\rho_A^2$。$\square$

当 $a$ 固定而 $b\to\infty$ 时，上下界趋于一致：随机纯态的一小部分几乎是最大混合的。精确的均值由
Page 公式给出，我们只陈述而不加证明。

**定理 7.13（Page 公式）。** 对 $a\le b$，
$$
\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}.
$$
它由 D. N. Page 提出猜想（Phys. Rev. Lett. 71, 1291, 1993），并由 S. K. Foong 和
S. Kanno（Phys. Rev. Lett. 72, 1148, 1994）、J. Sánchez-Ruiz（Phys. Rev. E 52, 5653, 1995）以及
S. Sen（Phys. Rev. Lett. 77, 1, 1996）证明。当 $1\ll a\le b$ 时它约为 $\log a-\frac a{2b}$；对两个量子比特，
它等于 $\frac13$。

$\rho_A$ 的分布是**诱导测度**（K. Życzkowski and H.-J. Sommers, J. Phys. A 34,
7111, 2001）。`states.random_density_matrix(alg, rank=r)` 用 $k_c\times r$ 的 Ginibre 块 $G_c$ 构造
$G_cG_c^*$，再除以 $\sum_c\operatorname{Tr}G_cG_c^*$。对单个块，这就是上面 $b=r$ 时 $\rho_A$ 的分布。
对多个块，它是 $H_1\otimes\mathbb C^r$ 上一个 Haar 随机纯态到 $\mathcal M_1$ 的限制：取 $\mathcal M_2=M_r$ 和一个
常值融合规则应用定理 7.6，使得 $\mathcal M_1\otimes_f\mathcal M_2=B(H_1\otimes\mathbb C^r)$，其中 $G_c$ 是 $G$
位于块 $c$ 中的那些行。

**命题 7.14.** 对由 `random_density_matrix(alg, rank=r)` 在 $\bigoplus_cM_{k_c}$ 上抽取的 $\rho$，
$D=\sum_ck_c$，
$$
\mathbb E\operatorname{Tr}\rho^2=\frac{\sum_ck_c(k_c+r)}{D(Dr+1)} .
$$

*证明。* $\operatorname{Tr}\rho^2=\sum_c\operatorname{Tr}(G_cG_c^*)^2/\|g\|^4$，其中 $g\in\mathbb C^{Dr}$ 是全部矩阵元
构成的向量；与定理 7.11 一样，引理 7.9 和 7.10 给出 $\sum_ck_cr(k_c+r)/(Dr(Dr+1))$。$\square$

我们在 $\mathbb C^3\otimes\mathbb C^6$ 的 $10^4$ 个经 `partial_trace` 约化的纯态上以及诱导测度上检验定理 7.11
和 7.13，并对 $M_1\oplus M_2\oplus M_3$、$r=2$ 检验命题 7.14：

```python
a, b, n_samples = 3, 6, 10_000
g = torch.randn(n_samples, a * b, dtype=dt)                 # complex Gaussian: E|g_i|^2 = 1
psi = g / g.norm(dim=-1, keepdim=True)
AB_alg = tensor_product(TypeIAlgebra([a], [a], precision='double'), TypeIAlgebra([b], [b], precision='double'))
pure = DensityMatrix(AB_alg, matrix=(psi[:, :, None] * psi.conj()[:, None, :])[:, None],
                     validate=False)                        # pure by construction
rho_A = partial_trace(pure, keep=1)
purity_exact = (a + b) / (a * b + 1)
page = sum(1 / k for k in range(b + 1, a * b + 1)) - (a - 1) / (2 * b)

def mean_se(x):
    return x.mean().item(), x.std().item() / math.sqrt(x.numel())

for name, rho_s in (('reduced pure states', rho_A),
                    ('induced measure', states.random_density_matrix(rho_A.algebra, 40_000, rank=b))):
    (m_p, e_p), (m_S, e_S) = mean_se(rho_s.purity()), mean_se(rho_s.von_neumann_entropy())
    print(f"{name}: E Tr rho^2 = {m_p:.4f} +- {e_p:.4f} (exact {purity_exact:.4f}), "
          f"E S = {m_S:.4f} +- {e_S:.4f} (Page {page:.4f})")
    assert abs(m_p - purity_exact) < 5 * e_p and abs(m_S - page) < 5 * e_S     # Theorems 7.11, 7.13
    assert math.log((a * b + 1) / (a + b)) < m_S < math.log(a)                # Corollary 7.12

alg3 = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')
m_p, e_p = mean_se(states.random_density_matrix(alg3, 40_000, rank=2).purity())
print(f"E Tr rho^2 = {m_p:.4f} +- {e_p:.4f} (exact {26 / 78:.4f})")
assert abs(m_p - 26 / 78) < 5 * e_p                                          # Proposition 7.14
```

> **物理意义。** 典型的纯态在系统中任意一个不超过一半的部分与其余部分之间几乎是最大纠缠的。
> 随着该部分增大，纠缠像 $\log a$ 一样增长，越过中点后又下降：这就是 *Page 曲线*，它描述了信息
> 如何随 Hawking 辐射离开黑洞，也解释了为什么封闭系统的各个部分看起来是热的。

## 7.6 可与之比较的工具

在单个 Hilbert 空间上，QuTiP（`Qobj.ptrace`、`entropy_vn`、`rand_ket`、带 Ginibre 分布和秩参数的
`rand_dm`）和 `numpy`（如上所示，一次 reshape 加一次 `einsum`）可以计算偏迹、熵和诱导态，并且在
普通积上与本库结果一致。带守恒荷的张量网络库——ITensor（量子数块稀疏张量）、TeNPy（带荷的
腿、按荷分解的纠缠谱），以及针对非阿贝尔对称性的 TensorKit.jl（融合树）——对大系统以同样的方式
融合扇区，并给出对称性分辨纠缠。本库增加的是代数层面：积是带有扇区和迹的冯·诺依曼代数，
融合规则可以是标签的任意映射，偏迹是一个信道，对每种迹都有其对偶，并且对融合积给出正确的
限制（定理 7.3、7.6），而且这一切都是批量化的。本库没有张量网络，并且 `fused_tensor_product`
接受的规则对每一对荷只给出一个结果，正如阿贝尔对称性那样；像 $SU(2)$ 的 $j_1\otimes j_2=\bigoplus_jj$
这样带有 Clebsch–Gordan 系数的非阿贝尔规则，在这里不是融合积（本库对自旋 $\frac12$ 链上的 $SU(2)$
另行处理，见第 4 章）。

## 习题

1. 计算 $\operatorname{Tr}_2$ 关于 $\mathcal M_1\otimes\mathcal M_2$ 和 $\mathcal M_1$ 上的
   $\tau_{\mathrm{vN}}$ 的对偶，并用 `partial_trace_channel(alg12).adjoint('tau_vN')` 加以验证。
2. （次可加性。）设 $\rho$ 是 $\mathcal M_1\otimes\mathcal M_2$ 上的态，其约化态 $\rho_1,\rho_2$ 可逆。
   证明 $S(\rho)\le S(\rho_1)+S(\rho_2)$，等号成立当且仅当 $\rho=\rho_1\otimes\rho_2$。
3. 设 `site` 为荷为 $0,1$ 的 $\mathbb C\oplus\mathbb C$（空格点或占据格点），并按加法与自身融合。
   对 $\psi=(|01\rangle+|10\rangle)/\sqrt2$（两个格点共享一个粒子）计算 $E(\psi)$ 及其两部分，并描述
   $(|00\rangle+|11\rangle)/\sqrt2$ 的向量态在融合代数上的限制。

<details>
<summary>解答</summary>

1. $\mathcal M_1$ 上的权重为 $1/(C_1k_c)$，积代数上的权重为 $1/(C_1C_2k_cl_d)$。对
   $A\in\mathcal M_1$，由命题 7.2(a)，
   $\tau_{\mathrm{vN}}((\operatorname{Tr}_2X)A)=\sum_c\frac1{C_1k_c}\sum_d\operatorname{Tr}(X_{cd}(A_c\otimes1))
   =\sum_{c,d}\frac1{C_1C_2k_cl_d}\operatorname{Tr}\big(X_{cd}\,C_2l_d(A_c\otimes1)\big)$，所以对偶为
   $A\mapsto\bigoplus_{c,d}C_2l_d\,A_c\otimes1_{l_d}=A\otimes Z$，其中 $Z=\bigoplus_dC_2l_d1_{l_d}$ 是
   $\mathcal M_2$ 的中心元。它不保单位：在本库的三种迹中，只有 blunt 迹使 $\operatorname{Tr}_2$ 成为
   $A\mapsto A\otimes1$ 的对偶（一般地，由同样的计算，$\mathcal M_1$ 上的 $\operatorname{tr}_w$ 与积代数上的
   权重 $w_{cd}=w_c$ 也能做到这一点）。

```python
dual_tau = composite.partial_trace_channel(alg12, keep=1).adjoint('tau_vN')
A = random_operator(alg1)
Z = alg2.central([alg2.C * l for l in alg2.k_factors])
assert close(dual_tau(A), kron(A, Z, alg12))
```

2. $\rho_1\otimes\rho_2$ 是一个态（定理 7.1(c)），可逆，特征向量为 $u_i\otimes v_j$，
   特征值为 $\lambda_i\mu_j$，所以 $\log(\rho_1\otimes\rho_2)=\log\rho_1\otimes1+1\otimes\log\rho_2$
   （定理 1.5(a)）。由偏迹的定义，
   $\operatorname{Tr}\rho(\log\rho_1\otimes1)=\operatorname{Tr}\rho_1\log\rho_1$，对第二个因子同理。
   因此 $S(\rho\|\rho_1\otimes\rho_2)=S(\rho_1)+S(\rho_2)-S(\rho)$，它 $\ge0$，等号成立当且仅当
   $\rho=\rho_1\otimes\rho_2$（推论 3.6）。这个差就是互信息。

```python
rho = states.random_density_matrix(alg12, batch_size=4)
rho1, rho2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
mutual = rho1.von_neumann_entropy() + rho2.von_neumann_entropy() - rho.von_neumann_entropy()
assert (mutual > 0).all()
assert torch.allclose(mutual, states.relative_entropy(rho, kron(rho1, rho2, alg12)))
```

3. 扇区为 $N=0,1,2$，大小分别为 $1,2,1$。$\psi$ 位于 $N=1$ 中，$\psi_0=|01\rangle/\sqrt2$，
   $\psi_1=|10\rangle/\sqrt2$，所以 $p=(\frac12,\frac12)$；由于各块都是一维的，
   $E(\hat\psi_c)=0$：由命题 7.8，$E(\psi)=H(p)=\log2$ 全部是数纠缠，可用的构型部分为 $0$。
   向量 $(|00\rangle+|11\rangle)/\sqrt2$ 是 $N=0$ 与 $N=2$ 的叠加；在块对角的融合代数上，它的向量态是
   乘积态的混合 $\frac12|00\rangle\langle00|+\frac12|11\rangle\langle11|$（第 3 章习题 3）。

```python
site = TypeIAlgebra([1, 1], [1, 1], precision='double', charges=[0, 1])
S2 = fused_tensor_product(site, site)
assert S2.k_factors == [1, 2, 1]
mat = torch.zeros(1, 3, 2, 2, dtype=dt)
mat[0, 1] = 0.5                                       # |psi><psi| with psi = (|01> + |10>)/sqrt 2
r1 = partial_trace(DensityMatrix(S2, matrix=mat), keep=1)
assert torch.allclose(r1.sector_probabilities()[0], torch.tensor([0.5, 0.5], dtype=torch.float64))
assert abs(r1.von_neumann_entropy().item() - math.log(2)) < 1e-12
```

</details>

**下一章：** 第 8 章（编写中）讨论动力学：带扇区的代数中的 Schrödinger 方程与 Lindblad 方程、
扇区之间的跃迁以及后选择。
