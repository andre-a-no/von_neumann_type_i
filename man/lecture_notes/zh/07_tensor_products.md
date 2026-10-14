[English](../07_tensor_products.md) · [Español](../es/07_tensor_products.md) · [Français](../fr/07_tensor_products.md) · [Deutsch](../de/07_tensor_products.md) · **中文** · [日本語](../ja/07_tensor_products.md) · [Русский](../ru/07_tensor_products.md)

# 7. 张量积与纠缠

本章研究两个带扇区的代数的张量积（在物理上，它是两个系统合在一起的代数；见*物理意义*方框）。
有扇区时存在两种积：普通的积，其中每个部分各自保持自己的荷；以及*融合*积，其中只有总荷守恒。
我们证明这两种积的结构，把偏迹定义为一个嵌入的对偶，并研究纠缠。

**你将学到**扇区如何相乘，为什么偏迹是一个信道、归一化后又是一个有限指标的
条件期望，荷守恒态的约化态是什么，纠缠如何分成数的部分与构型的部分，以及随机态的平均纯度与平均熵。
本章的数值检验见[教程 7](../../tutorial/zh/07_tensor_products.md)。

在本章中，$`\mathcal M_1=\bigoplus_{c=1}^{C_1}M_{k_c}(\mathbb C)`$ 作用在
$`H_1=\bigoplus_c\mathbb C^{k_c}`$ 上，极小中心投影为 $`P_c`$；
$`\mathcal M_2=\bigoplus_{d=1}^{C_2}M_{l_d}(\mathbb C)`$ 作用在 $`H_2=\bigoplus_d\mathbb C^{l_d}`$ 上，
极小中心投影为 $`Q_d`$。$`\mathrm{Tr}`$ 是算子所作用的空间上的迹；在 $`\mathcal M_i`$ 上它就是第 3 章的
blunt 迹。对 $`\mathbb C^k\otimes\mathbb C^l`$ 上的 $`X`$，记
$`X_{(a,b),(a',b')}=\langle e_a\otimes e_b,X\,e_{a'}\otimes e_{b'}\rangle`$，
并通过 $`e_a\otimes e_b\mapsto e_{(a-1)l+b}`$ 把 $`\mathbb C^k\otimes\mathbb C^l`$ 与 $`\mathbb C^{kl}`$ 等同
（Kronecker 约定），从而 $`(A\otimes B)_{(a,b),(a',b')}=A_{aa'}B_{bb'}`$。

<a id="s-7-1"></a>

## 7.1 两个代数的张量积

对 $`H_1`$ 上的 $`*`$-代数 $`\mathcal M_1`$ 和 $`H_2`$ 上的 $`*`$-代数 $`\mathcal M_2`$，**张量积**
$`\mathcal M_1\otimes\mathcal M_2`$ 是 $`H_1\otimes H_2`$ 上的算子 $`A\otimes B`$（$`A\in\mathcal M_1`$，$`B\in\mathcal M_2`$）
的线性张成。由于 $`(A\otimes B)(A'\otimes B')=AA'\otimes BB'`$，$`(A\otimes B)^*=A^*\otimes B^*`$ 且 $`1\otimes1=1`$，
它是一个 $`*`$-代数，因而是一个冯·诺依曼代数（定理 2.2）。

<a id="r-7-1"></a>

**定理 7.1.** (a) $`H_1\otimes H_2=\bigoplus_{c,d}\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$，这是投影
$`P_c\otimes Q_d`$ 的值域的正交和，并且
```math
\mathcal M_1\otimes\mathcal M_2=\bigoplus_{c,d}M_{k_c}\otimes M_{l_d}\cong\bigoplus_{c,d}M_{k_cl_d}(\mathbb C):
```
它由所有满足 $`X=\sum_{c,d}(P_c\otimes Q_d)X(P_c\otimes Q_d)`$ 的算子组成，在各个
$`\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$ 上的块是任意的。(b) 它的极小中心投影是 $`C_1C_2`$ 个投影
$`P_c\otimes Q_d`$，并且 $`\dim(\mathcal M_1\otimes\mathcal M_2)=\dim\mathcal M_1\dim\mathcal M_2`$。
(c) 对 $`\mathcal M_1`$ 上的迹 $`\mathrm{tr}_w`$ 和 $`\mathcal M_2`$ 上的迹 $`\mathrm{tr}_v`$（定理 3.1），
$`\mathrm{tr}_w(A)\mathrm{tr}_v(B)=\mathrm{tr}_{w\otimes v}(A\otimes B)`$，其中
$`(w\otimes v)_{cd}=w_cv_d`$。特别地，积代数上的 `Tr_blunt`、`Tr_norm` 和 `tau_vN` 分别是各因子上同名迹的乘积。

*证明。* (a) 诸 $`P_c\otimes Q_d`$ 是两两正交的投影，其和为 $`1\otimes1`$，且
$`P_c\otimes Q_d`$ 的值域是 $`\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$。由于 $`P_c`$、$`Q_d`$ 是中心元，
$`(P_c\otimes Q_d)(A\otimes B)=AP_c\otimes BQ_d=(A\otimes B)(P_c\otimes Q_d)`$，所以
$`\mathcal M_1\otimes\mathcal M_2`$ 的每个元素都是块对角的。反之，设 $`E^c_{aa'}\in\mathcal M_1`$ 和
$`F^d_{bb'}\in\mathcal M_2`$ 分别是块 $`c`$ 与块 $`d`$ 的矩阵单位，则在
$`\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$ 内 $`E^c_{aa'}\otimes F^d_{bb'}=|e_a\otimes e_b\rangle\langle e_{a'}\otimes e_{b'}|`$；
这些正是 $`M_{k_cl_d}(\mathbb C)`$ 的矩阵单位，所以每个支撑在块 $`(c,d)`$ 中的算子都属于
$`\mathcal M_1\otimes\mathcal M_2`$。(b) 全矩阵代数直和的中心由到各直和项上的投影张成（第 2 章习题 1），且
$`\sum_{c,d}k_c^2l_d^2=\sum_ck_c^2\sum_dl_d^2`$。(c)
$`\mathrm{Tr}(X\otimes Y)=\sum_{a,b}X_{aa}Y_{bb}=\mathrm{Tr}\,X\mathrm{Tr}\,Y`$，所以
$`\mathrm{tr}_{w\otimes v}(A\otimes B)=\sum_{c,d}w_cv_d\mathrm{Tr}\,A_c\mathrm{Tr}\,B_d=\mathrm{tr}_w(A)\mathrm{tr}_v(B)`$。
三种迹的权重 $`1`$、$`1/k_c`$、$`1/(C_1k_c)`$（3.1 节）相乘后得到 $`1`$、
$`1/(k_cl_d)`$、$`1/(C_1C_2k_cl_d)`$，正是积代数上同名迹的权重。$`\square`$

*库中：* 积是 `tensor_product`，$`A\otimes B`$ 是 `kron`；在 $`(\mathbb C\oplus M_2)\otimes(M_2\oplus M_3)`$ 上检验定理 7.1——[教程，7.1 节](../../tutorial/zh/07_tensor_products.md#s-7-1)。

> **物理意义**（物理解释，而非本课程的定理）。如果两个系统各自服从自己的超选择定则——两个不能交换粒子的盒子——那么联合系统
> 对每一*对*荷都有一个扇区，并且没有可观测量连接不同的荷对。一个联合扇区的权重是两个权重之积
> （定理 7.1(c)）。

<a id="s-7-2"></a>

## 7.2 偏迹

**定义。** 设 $`\rho\in\mathcal M_1\otimes\mathcal M_2`$。**偏迹**
$`\mathrm{Tr}_2\,\rho\in\mathcal M_1`$ 是满足下式的元素：
```math
\mathrm{Tr}\big((\mathrm{Tr}_2\,\rho)\,A\big)=\mathrm{Tr}\big(\rho\,(A\otimes1)\big)\qquad\text{for all }A\in\mathcal M_1 .
```
由定理 3.2（应用于线性泛函 $`A\mapsto\mathrm{Tr}(\rho(A\otimes1))`$ 和 $`\mathcal M_1`$ 的忠实迹
$`\mathrm{Tr}`$），它存在且唯一：偏迹是嵌入 $`A\mapsto A\otimes1`$ 关于 blunt 迹的对偶。
$`\mathrm{Tr}_1`$ 的定义与此相同；对态 $`\rho`$，$`\mathrm{Tr}_2\,\rho`$ 与
$`\mathrm{Tr}_1\,\rho`$ 称为**约化态**。对全代数 $`\mathcal M_i=B(H_i)`$，这就是
**Hilbert 空间偏迹** $`\mathrm{Tr}_{H_2}:B(H_1\otimes H_2)\to B(H_1)`$。

<a id="r-7-2"></a>

**命题 7.2.** (a) 对 $`\rho=\bigoplus_{c,d}\rho_{cd}`$，
```math
(\mathrm{Tr}_2\,\rho)_c=\sum_d\mathrm{tr}_{l_d}\,\rho_{cd},\qquad
(\mathrm{tr}_l\,X)_{aa'}=\sum_{b=1}^lX_{(a,b),(a',b)}\quad(X\text{ on }\mathbb C^k\otimes\mathbb C^l),
```
并且 $`\mathrm{Tr}_2\,\rho=\mathrm{Tr}_{H_2}\,\rho`$，其中 $`(\mathrm{Tr}_{H_2}\,X)_{ii'}=\sum_jX_{(i,j),(i',j)}`$。
(b) $`\mathrm{Tr}_2(A\otimes B)=(\mathrm{Tr}\,B)\,A`$。
(c) 对 $`A,A'\in\mathcal M_1`$，$`\mathrm{Tr}_2\big((A\otimes1)X(A'\otimes1)\big)=A\,(\mathrm{Tr}_2\,X)\,A'`$。
(d) $`\mathrm{Tr}(\mathrm{Tr}_2\,X)=\mathrm{Tr}\,X`$，并且 $`\mathrm{Tr}_2`$ 把态映为态。

*证明。* (a) 取 $`A=\bigoplus_cA_c`$，则 $`A\otimes1=\bigoplus_{c,d}A_c\otimes1_{l_d}`$，且
```math
\mathrm{Tr}\big(\rho(A\otimes1)\big)=\sum_{c,d}\sum_{a,a',b}(\rho_{cd})_{(a,b),(a',b)}(A_c)_{a'a}
=\sum_c\mathrm{Tr}\Big(\Big(\sum_d\mathrm{tr}_{l_d}\,\rho_{cd}\Big)A_c\Big),
```
所以右边的元素满足定义，即为 $`\mathrm{Tr}_2\,\rho`$。对 $`B(H_i)`$ 作同样的计算即得
$`\mathrm{Tr}_{H_2}`$ 的公式；对块对角的 $`\rho`$，$`i,i'`$ 属于不同块的矩阵元
$`\rho_{(i,j),(i',j)}`$ 为零，所以 $`\mathrm{Tr}_{H_2}\,\rho`$ 具有 (a) 中的块。(b) $`\mathrm{tr}_l(X\otimes Y)=(\mathrm{Tr}\,Y)X`$。
(c) 对 $`B\in\mathcal M_1`$，
$`\mathrm{Tr}\big(\mathrm{Tr}_2((A\otimes1)X(A'\otimes1))B\big)=\mathrm{Tr}\big(X(A'BA\otimes1)\big) =\mathrm{Tr}\big(A(\mathrm{Tr}_2\,X)A'B\big)`$；再用唯一性。(d) 取 $`A=1`$。若 $`\rho`$ 是态，则
$`A\mapsto\mathrm{Tr}(\rho(A\otimes1))`$ 是正的，因为 $`A^*A\otimes1=(A\otimes1)^*(A\otimes1)`$，并且它在
$`A=1`$ 处取值 $`1`$，所以它的密度 $`\mathrm{Tr}_2\,\rho`$ 是一个态（定理 3.2）。$`\square`$

*库中：* $`\mathrm{Tr}_2`$ 是 `partial_trace`；对照 `numpy.einsum` 检验命题 7.2(a)、(b)——[教程，7.2 节](../../tutorial/zh/07_tensor_products.md#s-7-2)。

信道是第 6 章的主题；这里只需要如下定义：线性映射 $`\Phi`$ 称为**完全正**的，如果对每个 $`n`$，
$`\mathrm{id}_n\otimes\Phi`$ 在 $`M_n(\mathbb C)\otimes(\cdot)`$ 上是正的；如果它还保迹，则称为**信道**。

<a id="r-7-3"></a>

**定理 7.3.** 设 $`j`$ 取遍 $`H_2`$ 的标准基，$`K_j=1\otimes\langle j|:H_1\otimes H_2\to H_1`$，
$`K_j(x\otimes y)=\langle j,y\rangle x`$。则
```math
\mathrm{Tr}_{H_2}\,X=\sum_jK_jXK_j^*,\qquad\sum_jK_j^*K_j=1,\qquad\sum_jK_j^*AK_j=A\otimes1 .
```
因此 $`\mathrm{Tr}_2`$ 是从 $`\mathcal M_1\otimes\mathcal M_2`$ 到 $`\mathcal M_1`$ 的信道，并且它关于
blunt 迹的对偶是 $`A\mapsto A\otimes1`$。

*证明。* $`K_j^*x=x\otimes e_j`$，所以 $`\langle e_a,K_jXK_j^*e_{a'}\rangle=X_{(a,j),(a',j)}`$；对 $`j`$ 求和
即得命题 7.2(a)。此外 $`\sum_jK_j^*K_j(x\otimes y)=\sum_j\langle j,y\rangle x\otimes e_j=x\otimes y`$，
且 $`\sum_jK_j^*AK_j=\sum_jA\otimes|j\rangle\langle j|=A\otimes1`$。对 $`Y\ge0`$，
$`(\mathrm{id}_n\otimes\mathrm{Tr}_{H_2})(Y)=\sum_j(1_n\otimes K_j)Y(1_n\otimes K_j)^*\ge0`$，因为
$`\langle x,TYT^*x\rangle=\langle T^*x,YT^*x\rangle\ge0`$；这一性质在子代数
$`\mathcal M_1\otimes\mathcal M_2`$ 上依然成立。最后 $`\mathrm{Tr}\,\sum_jK_jXK_j^*A=\mathrm{Tr}\,X\sum_jK_j^*AK_j`$，
当 $`A=1`$ 时这就是保迹性，一般情形下则是对偶关系。$`\square`$

*库中：* 定理 7.3 的信道是 `composite.partial_trace_channel`，连同其对偶——[教程，7.2 节](../../tutorial/zh/07_tensor_products.md#s-7-2)。

> **物理意义**（物理解释，而非本课程的定理）。约化态 $`\mathrm{Tr}_2\,\rho`$ 是一个只能测量第一个系统的观察者所能得知的
> 全部信息——定义表明，它重现了每个可观测量 $`A\otimes1`$ 的期望值。约化态相同的联合态在关联上有所不同。

**作为条件期望的偏迹。** 设 $`\mathcal M_2=M_m(\mathbb C)`$，于是 $`\mathcal M_1\otimes M_m`$ 的扇区
就是 $`\mathcal M_1`$ 的扇区（定理 7.1）；令 $`\iota(A)=A\otimes1_m`$，
$`\mathcal E=\frac1m\mathrm{Tr}_2:\mathcal M_1\otimes M_m\to\mathcal M_1`$。设 $`e_{ab}`$ 为 $`M_m`$ 的矩阵单位。

<a id="r-7-4"></a>

**命题 7.4.** (a) $`\iota\circ\mathcal E`$ 是 $`\mathcal M_1\otimes M_m`$ 到 $`\mathcal M_1\otimes1`$ 上的条件期望
（2.6 节），并且对 $`\mathcal M_1\otimes M_m`$ 的每个忠实迹 $`\mathrm{tr}_w`$，它是唯一保持
$`\mathrm{tr}_w`$ 的条件期望。(b) $`m^2`$ 个元素 $`u_{ab}=\sqrt m\,(1\otimes e_{ab})`$ 对所有 $`X`$ 满足
$`X=\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*X))`$，且 $`\sum_{a,b}u_{ab}u_{ab}^*=m^2\,1`$。(c) 每个对所有 $`X`$
满足 $`X=\sum_kv_k\,\iota(\mathcal E(v_k^*X))`$ 的族 $`(v_k)`$ 都有 $`\sum_kv_kv_k^*=m^2\,1`$。

*证明。* (a) $`\mathcal E`$ 是正的（定理 7.3），$`\iota`$ 也是正的
（$`T^*T\otimes1=(T\otimes1)^*(T\otimes1)`$）；由命题 7.2(b)，$`\mathcal E(A\otimes1)=A`$，而双模性质
就是命题 7.2(c)。对 $`\mathcal M_1`$ 上具有相同权重 $`w_c`$ 的迹也记为 $`\mathrm{tr}_w`$；则
$`\mathrm{tr}_w(Y\otimes1)=m\mathrm{tr}_w(Y)`$。把 $`\mathrm{Tr}_2`$ 的定义应用于 $`AP_c`$，乘以 $`w_c`$ 并对 $`c`$ 求和，得到
$`\mathrm{tr}_w(X(A\otimes1))=\mathrm{tr}_w((\mathrm{Tr}_2\,X)A)`$，取 $`A=1`$ 得
$`\mathrm{tr}_w(\iota\mathcal E(X))=\mathrm{tr}_w(X)`$。若 $`F=\iota\circ F_0`$ 是到 $`\mathcal M_1\otimes1`$ 上另一个保持
$`\mathrm{tr}_w`$ 的条件期望，则由于 $`F(X(A\otimes1))=F(X)(A\otimes1)`$，
$`\mathrm{tr}_w(X(A\otimes1))=\mathrm{tr}_w(F(X)(A\otimes1))=m\mathrm{tr}_w(F_0(X)A)`$，所以对所有
$`A\in\mathcal M_1`$ 有 $`\mathrm{tr}_w((\mathrm{Tr}_2\,X-mF_0(X))A)=0`$，由忠实性得 $`F_0=\mathcal E`$
（取 $`A=(\mathrm{Tr}_2\,X-mF_0(X))^*`$）。(b) 写 $`X=\sum_{a,b}X_{ab}\otimes e_{ab}`$，其中
$`X_{ab}\in\mathcal M_1`$（定理 7.1(a)）。由命题 7.2(b)，
$`\mathcal E(u_{ab}^*X)=\frac{\sqrt m}m\sum_{b'}\mathrm{Tr}_2(X_{ab'}\otimes e_{bb'})=X_{ab}/\sqrt m`$，所以
$`u_{ab}\,\iota(\mathcal E(u_{ab}^*X))=X_{ab}\otimes e_{ab}`$；并且
$`\sum_{a,b}u_{ab}u_{ab}^*=m\sum_{a,b}1\otimes e_{aa}=m^21`$。(c) 由命题 7.2(a)，$`\mathrm{Tr}_2(Y^*)=(\mathrm{Tr}_2\,Y)^*`$，
所以 $`Y=u_{ab}`$ 的展开式取伴随后为
$`u_{ab}^*=\sum_k\iota(\mathcal E(u_{ab}^*v_k))\,v_k^*`$，再由 (b)，
$`\sum_{a,b}u_{ab}u_{ab}^*=\sum_k\big(\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*v_k))\big)v_k^*=\sum_kv_kv_k^*`$。$`\square`$

像 (c) 中那样的族称为 $`\mathcal E`$ 的一个*拟基*，而 $`\sum_kv_kv_k^*`$ 由 (c) 不依赖于拟基的选取、由 (b)
存在，它就是 $`\mathcal E`$ 的**指标** $`\mathrm{Ind}\,\mathcal E=m^2`$（这是一个定义；名称沿用 Y. Watatani, Mem. Amer. Math.
Soc. 83, no. 424, 1990）。把 $`\bigoplus_cM_{k_c}`$ 包含到 $`\bigoplus_sM_{n_s}`$ 中的一个包含，其**包含矩阵**的元素
$`\Lambda_{sc}`$ 是子代数的块 $`c`$ 在较大代数的块 $`s`$ 中出现的副本个数（第 9 章把它用于代数塔）。对因子而言，
$`M_k\otimes1\subset M_k\otimes M_m`$：把 $`\mathbb C^{km}`$ 的基 $`e_a\otimes e_b`$ 先按 $`b`$ 排序，则
$`A\otimes1_m=\mathrm{diag}(A,\dots,A)`$，含 $`m`$ 个副本，所以 $`\Lambda=(m)`$，并且
$`m^2=\|\Lambda\|^2=\dim(M_k\otimes M_m)/\dim M_k`$。（延伸阅读，此处不使用：V. F. R. Jones,
Invent. Math. 72, 1, 1983，关于子因子的指标；F. M. Goodman, P. de la Harpe and V. F. R. Jones, *Coxeter graphs and
towers of algebras*, Springer, 1989，关于包含矩阵。）对一般的 $`\mathcal M_1`$，这个包含是 $`C_1`$ 个包含
$`M_{k_c}\otimes1\subset M_{k_c}\otimes M_m`$ 的直和，$`\Lambda=m1_{C_1}`$。

*库中：* `inclusions.from_tensor_product`，连同命题 7.4 的条件期望、拟基和指标——[教程，7.2 节](../../tutorial/zh/07_tensor_products.md#s-7-2)。

<a id="s-7-3"></a>

## 7.3 荷的融合

现在设 $`\mathcal M_1`$ 的扇区 $`c`$ 带有荷 $`q_c`$，$`\mathcal M_2`$ 的扇区 $`d`$ 带有荷 $`r_d`$，并且一个
**融合规则** $`f`$ 给出总荷 $`f(q_c,r_d)`$——对粒子数是加法，对费米子宇称是模 $`2`$ 加法，也可以是任何
其他映射。令
```math
H_q=\bigoplus_{f(q_c,r_d)=q}\mathbb C^{k_c}\otimes\mathbb C^{l_d},\qquad
\Pi_q=\sum_{f(q_c,r_d)=q}P_c\otimes Q_d,\qquad n_q=\sum_{f(q_c,r_d)=q}k_cl_d .
```
**融合积** $`\mathcal M_1\otimes_f\mathcal M_2=\{\Pi_q\}'`$ 由 $`H_1\otimes H_2`$ 上保持总荷守恒的算子组成。

<a id="r-7-5"></a>

**命题 7.5.** (a) $`\mathcal M_1\otimes_f\mathcal M_2=\bigoplus_qB(H_q)\cong\bigoplus_qM_{n_q}(\mathbb C)`$，
极小中心投影为 $`\Pi_q`$。(b) $`\mathcal M_1\otimes\mathcal M_2\subseteq\mathcal M_1\otimes_f\mathcal M_2`$，
等号成立当且仅当每个 $`H_q`$ 只包含一个对 $`(c,d)`$。(c) 若荷是实数且 $`f`$ 是加法，则
$`\mathcal M_1\otimes_f\mathcal M_2=\{N\}'`$，其中
$`N=N_1\otimes1+1\otimes N_2`$，$`N_1=\sum_cq_cP_c`$，$`N_2=\sum_dr_dQ_d`$。

*证明。* (a) 每个对恰好属于一个 $`H_q`$，所以诸 $`\Pi_q`$ 是两两正交的投影，其和为 $`1`$，而 $`X`$ 与它们
全都可交换当且仅当 $`X=\sum_q\Pi_qX\Pi_q`$，其中 $`B(H_q)`$ 中的块是任意的；中心同定理 7.1(b)。(b)
$`\mathcal M_1\otimes\mathcal M_2`$ 的每个元素都与每个 $`P_c\otimes Q_d`$ 可交换，因而与它们的和 $`\Pi_q`$ 可交换。
若每个 $`H_q`$ 都只是一个 $`\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$，则两个代数都是
$`\bigoplus_{c,d}B(\mathbb C^{k_c}\otimes\mathbb C^{l_d})`$；若 $`H_q`$ 包含两个对，则 $`B(H_q)`$ 包含把一个对的块
映到另一个对的块的算子，它们不是块对角的。(c)
$`N=\sum_{c,d}(q_c+r_d)P_c\otimes Q_d=\sum_qq\,\Pi_q`$，所以非零的 $`\Pi_q`$ 是 $`N`$ 的谱投影。它们是 $`N`$ 的
多项式（定理 1.5(a) 的证明），而 $`N`$ 是它们的线性组合，所以 $`X`$ 与 $`N`$ 可交换当且仅当它与每个
$`\Pi_q`$ 可交换。$`\square`$

**例：硬核玻色子。** 在 $`L`$ 个格点上，$`(\mathbb C^2)^{\otimes L}`$ 带有占据数
$`N=\sum_in^{(i)}`$，$`n=\mathrm{diag}(0,1)`$；本征空间 $`N=q`$ 由恰有 $`q`$ 个格点被占据的 $`\binom Lq`$ 个乘积
基向量张成，所以保持粒子数的算子 $`\{N\}'`$ 构成
$`\bigoplus_qM_{\binom Lq}`$（命题 4.10(2)）：两个格点时为 $`M_1\oplus M_2\oplus M_1`$（荷 $`0,1,2`$），
三个格点时为 $`M_1\oplus M_3\oplus M_3\oplus M_1`$。由命题 7.5(c)，在基的顺序意义下，融合积就是五个格点的总占据数
所对应的 $`\{N\}'`$，并且确实有
$`n_q=\sum_c\binom2c\binom3{q-c}=\binom5q`$（Vandermonde 恒等式：五个格点的一个 $`q`$ 元子集由前两个格点中的
$`c`$ 个和后三个格点中的 $`q-c`$ 个组成）。

*库中：* 融合积是 `fused_tensor_product`；硬核玻色子与命题 7.5——[教程，7.3 节](../../tutorial/zh/07_tensor_products.md#s-7-3)。

融合代数的维数是 $`\sum_q\binom5q^2=\binom{10}5=252`$（又是 Vandermonde 恒等式），普通积的维数是
$`(1+4+1)(1+9+9+1)=120`$（定理 7.1(b)）。由定理 7.1(a) 和命题 7.5(b)，普通积由融合积中关于对 $`(c,d)`$
块对角的元素组成；其余 $`132`$ 维由这样的算子张成：它们把一个对块 $`(c,d)`$ 映到同一个 $`H_q`$ 中另一个对块
$`(c',d')`$，$`c'\ne c`$（因而 $`d'\ne d`$）。它们在保持总荷守恒的同时改变每个部分的荷，即在这里把粒子在两部分之间移动。

$`\mathcal M_1\otimes_f\mathcal M_2`$ 上的态可以在不同的对块之间具有相干性，因此命题 7.2 不适用。
在 $`B(H_1)`$ 上令 $`E_1(X)=\sum_cP_cXP_c`$（2.6 节的条件期望）；它映到 $`\mathcal M_1`$ 上，
并且对 $`A\in\mathcal M_1`$ 有 $`\mathrm{Tr}(E_1(X)A)=\mathrm{Tr}(XA)`$，因为 $`A=\sum_cP_cAP_c`$。

<a id="r-7-6"></a>

**定理 7.6.** 设 $`\rho`$ 是 $`\mathcal M_1\otimes_f\mathcal M_2`$ 上的态 $`\omega`$ 的密度。它到 $`\mathcal M_1`$ 的
限制 $`A\mapsto\omega(A\otimes1)`$ 的密度为
```math
\rho_1=E_1\big(\mathrm{Tr}_{H_2}\,\rho\big).
```
$`\mathrm{Tr}_{H_2}\,\rho`$ 的块 $`(c,c')`$ 等于 $`\sum_d\mathrm{tr}_{l_d}\,\rho_{(c,d),(c',d)}`$，求和遍及满足
$`f(q_c,r_d)=f(q_{c'},r_d)`$ 的 $`d`$，其中 $`\rho_{(c,d),(c',d)}`$ 是 $`\rho`$ 从对 $`(c',d)`$ 到对 $`(c,d)`$ 的块。
因此，$`E_1`$ 对所有 $`\mathrm{Tr}_{H_2}\,\rho`$ 都平凡地作用，当且仅当不存在 $`c\neq c'`$ 和 $`d`$ 使得
$`f(q_c,r_d)=f(q_{c'},r_d)`$；例如当 $`f`$ 是加法且诸 $`q_c`$ 两两不同时即是如此。

*证明。* $`A\otimes1\in\mathcal M_1\otimes_f\mathcal M_2`$（命题 7.5(b)），且
$`\omega(A\otimes1)=\mathrm{Tr}(\rho(A\otimes1))=\mathrm{Tr}\big((\mathrm{Tr}_{H_2}\,\rho)A\big) =\mathrm{Tr}\big(E_1(\mathrm{Tr}_{H_2}\,\rho)A\big)`$。由于 $`E_1(\mathrm{Tr}_{H_2}\,\rho)\in\mathcal M_1`$，
由定理 3.2 中的唯一性，它就是该密度。块公式是 $`\mathrm{Tr}_{H_2}`$ 的命题 7.2(a)，其中 $`i`$ 在块 $`c`$ 中，
$`i'`$ 在块 $`c'`$ 中，$`j`$ 在块 $`d`$ 中；除非两个对位于同一个 $`H_q`$ 中，否则 $`\rho`$ 在 $`(c',d)`$ 与 $`(c,d)`$ 之间的块为零。
若存在 $`c\neq c'`$ 和 $`d`$ 使得 $`f(q_c,r_d)=f(q_{c'},r_d)=q`$，取单位向量 $`x\in\mathbb C^{k_c}`$、
$`x'\in\mathbb C^{k_{c'}}`$、$`y\in\mathbb C^{l_d}`$，并令 $`\psi=(x+x')\otimes y/\sqrt2\in H_q`$。则
$`|\psi\rangle\langle\psi|\in B(H_q)`$，而 $`\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|=\frac12|x+x'\rangle\langle x+x'|`$
具有非零块 $`\frac12|x\rangle\langle x'|`$，它被 $`E_1`$ 去掉。$`\square`$

由定理 7.6 的块公式，$`(\rho_1)_c=\sum_d\mathrm{tr}_{l_d}\,\rho_{(c,d),(c,d)}`$。

*库中：* `partial_trace` 和 `composite.partial_trace_channel` 对任意融合规则计算 $`\rho_1`$；对加法和宇称规则检验定理 7.6——[教程，7.3 节](../../tutorial/zh/07_tensor_products.md#s-7-3)。

> **物理意义**（物理解释，而非本课程的定理）。在超导体的 BCS 描述中，只有电子数的宇称守恒，因此整体态可以是一个区域的粒子数
> $`0`$ 与 $`2`$ 的叠加。一个其可观测量保持该区域粒子数守恒的观察者看不到这种相干性：对他而言，
> 约化态是定理 7.6 中经过块对角化（pinching）的那个态（这是本方框的数学内容），即不同粒子数的混合。

<a id="s-7-4"></a>

## 7.4 纯态的纠缠

<a id="r-7-7"></a>

**引理 7.7（Schmidt 分解）。** 每个单位向量 $`\psi\in\mathbb C^a\otimes\mathbb C^b`$ 都可以写成
$`\psi=\sum_{i=1}^rs_i\,u_i\otimes v_i`$，其中 $`s_i>0`$，$`\sum_is_i^2=1`$，且 $`u_i\in\mathbb C^a`$、$`v_i\in\mathbb C^b`$
分别规范正交。于是 $`\mathrm{Tr}_2|\psi\rangle\langle\psi|=\sum_is_i^2|u_i\rangle\langle u_i|`$，
$`\mathrm{Tr}_1|\psi\rangle\langle\psi|=\sum_is_i^2|v_i\rangle\langle v_i|`$：两个约化态具有相同的非零特征值 $`s_i^2`$。

*证明。* 设 $`\Psi`$ 为 $`a\times b`$ 矩阵 $`\Psi_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,\psi\rangle`$。由谱定理，
$`\Psi\Psi^*=\sum_{i=1}^a\lambda_i|u_i\rangle\langle u_i|`$，其中 $`u_i`$ 是一组规范正交基，
$`\lambda_i=\|\Psi^*u_i\|^2\ge0`$，并排序使得 $`\lambda_i>0`$ 恰好对 $`i\le r`$ 成立。对 $`i\le r`$ 令
$`s_i=\sqrt{\lambda_i}`$，$`w_i=\Psi^*u_i/s_i`$；则
$`\langle w_i,w_j\rangle=\langle u_i,\Psi\Psi^*u_j\rangle/(s_is_j)=\delta_{ij}`$，且当 $`i>r`$ 时 $`\Psi^*u_i=0`$。
因此 $`\Psi=\sum_i|u_i\rangle\langle u_i|\Psi=\sum_{i\le r}|u_i\rangle\langle\Psi^*u_i|=\sum_{i\le r}s_i|u_i\rangle\langle w_i|`$，
即 $`\Psi_{\alpha\beta}=\sum_is_i(u_i)_\alpha\overline{(w_i)_\beta}`$：这就是 $`v_i=\bar w_i`$ 时的结论，并且
$`\sum_is_i^2=\mathrm{Tr}\,\Psi\Psi^*=\|\psi\|^2=1`$。最后
$`|\psi\rangle\langle\psi|=\sum_{i,j}s_is_j|u_i\rangle\langle u_j|\otimes|v_i\rangle\langle v_j|`$，
由命题 7.2(b) 即得约化态。$`\square`$

$`\psi`$ 的**纠缠熵**为
$`E(\psi)=S(\mathrm{Tr}_2|\psi\rangle\langle\psi|)=S(\mathrm{Tr}_1|\psi\rangle\langle\psi|)=-\sum_is_i^2\log s_i^2`$。
它为零当且仅当 $`r=1`$（每一项 $`-s_i^2\log s_i^2`$ 都 $`\ge0`$，且仅当 $`s_i^2=1`$ 时为 $`0`$），即当且仅当
$`\psi=u\otimes v`$（此时 $`\Psi=uv^T`$ 的秩为 $`1`$，而 $`r=\mathrm{rank}\,\Psi\Psi^*=\mathrm{rank}\,\Psi`$）；并且
$`E(\psi)\le\log\min(a,b)`$，等号成立当且仅当较小的那个约化态是最大混合态
（推论 3.6，应用于两个约化态）。

有荷时，我们取具有确定总荷的单位向量 $`\psi`$，即位于融合积的某一个扇区 $`H_q`$ 中（不同 $`H_q`$ 的叠加在
$`\mathcal M_1\otimes_f\mathcal M_2`$ 上给出与混合相同的态，第 3 章习题 3）；记 $`\psi=\sum_{(c,d)}\psi_{cd}`$，其中
$`\psi_{cd}\in\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$。

<a id="r-7-8"></a>

**命题 7.8（对称性分辨纠缠）。** 设 $`f`$ 是加法，诸 $`q_c`$ 两两不同，诸 $`r_d`$ 也两两不同，从而对每个 $`c`$
至多有一个 $`d=d(c)`$ 满足 $`q_c+r_d=q`$。令 $`\psi_c=\psi_{c\,d(c)}`$（或 $`0`$），$`p_c=\|\psi_c\|^2`$，并在 $`p_c>0`$ 时令
$`\hat\psi_c=\psi_c/\sqrt{p_c}`$。则 $`\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|`$ 属于 $`\mathcal M_1`$（它与 $`N_1`$
可交换），并且等于 $`\bigoplus_cp_c\mathrm{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|`$；$`\mathcal M_2`$ 上的约化态在
$`d(c)`$ 处具有扇区概率 $`p_c`$；并且
```math
E(\psi)=H(p)+\sum_cp_c\,E(\hat\psi_c),\qquad H(p)=-\sum_cp_c\log p_c .
```

*证明。* 由定理 7.6 的块公式，$`\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|`$ 的块 $`(c,c')`$
是 $`\mathrm{tr}_{l_d}|\psi_{cd}\rangle\langle\psi_{c'd}|`$ 对满足 $`q_c+r_d=q=q_{c'}+r_d`$ 的 $`d`$ 之和，
这迫使 $`c=c'`$ 且 $`d=d(c)`$。所以约化态为
$`\bigoplus_c\mathrm{Tr}_2|\psi_c\rangle\langle\psi_c|`$，其块的迹为 $`p_c`$；交换两个因子的角色即得关于
$`\mathcal M_2`$ 的结论。由引理 7.7，$`E(\psi)`$ 是这个态的熵，再对
$`\sigma_c=\mathrm{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|`$、$`S(\sigma_c)=E(\hat\psi_c)`$ 应用定理 3.3 即得公式。$`\square`$

$`H(p)`$ 称为**数纠缠**，$`\sum_cp_cE(\hat\psi_c)`$ 称为**构型纠缠**。

*计算：* 对五个格点上两个粒子的随机态检验引理 7.7 和命题 7.8——[教程，7.4 节](../../tutorial/zh/07_tensor_products.md#s-7-4)。

> **物理意义**（物理解释，而非本课程的定理）。受粒子数超选择定则约束的各方只能把构型部分用作资源：这就是 H. M. Wiseman 和 J. A.
> Vaccaro 的*粒子纠缠*（Phys. Rev. Lett. 91, 097902, 2003）。$`H(p)`$ 在冷原子实验中作为*数熵*
> 被测量（A. Lukin et al., Science 364, 256, 2019）。

<a id="s-7-5"></a>

## 7.5 随机态与 Page 公式

设 $`g\in\mathbb C^N`$ 的各坐标相互独立，其实部与虚部是独立的 $`N(0,\tfrac12)`$ 随机变量（标准复高斯向量，
命题 A.16(a)），于是 $`\mathbb E|g_i|^2=1`$；$`g\neq0`$ 几乎必然成立（定理 A.20）。**Haar 随机纯态**是
$`\psi=g/\|g\|`$。它的分布是酉不变的：对每个酉矩阵 $`U`$，$`Ug`$ 与 $`g`$ 同分布（命题 A.16(a)），并且
$`Ug/\|Ug\|=U\psi`$。由命题 A.22，它是单位球面上唯一的酉不变分布；并且如该命题的证明所示，它是
Haar 酉矩阵（定理 5.3）第一列的分布。

<a id="r-7-9"></a>

**引理 7.9.** $`\|g\|`$ 与 $`\psi=g/\|g\|`$ 相互独立。

*证明。*（这是定理 A.20(b) 的一部分。）设 $`h\ge0`$ 与 $`\phi\ge0`$ 分别是 $`(0,\infty)`$ 和
$`\mathbb C^N\cong\mathbb R^{2N}`$ 的单位球面 $`S^{2N-1}`$ 上的 Borel 函数。由于 $`g`$ 有密度
$`\pi^{-N}e^{-|x|^2}`$（命题 A.16(a)），极坐标（取 $`n=2N`$ 的定理 A.18，以及其中定义的测度 $`\mu`$）给出
```math
\mathbb E\big[h(\|g\|)\phi(\psi)\big]=\pi^{-N}\!\int_{\mathbb R^{2N}}\!h(|x|)e^{-|x|^2}\phi\Big(\frac x{|x|}\Big)dx
=2N\pi^{-N}\int_0^\infty h(t)e^{-t^2}t^{2N-1}dt\int_{S^{2N-1}}\phi\,d\mu=:\alpha(h)\beta(\phi).
```
先取 $`\phi=1`$，再取 $`h=1`$，再取两者都为 $`1`$，得到 $`\mathbb Eh(\|g\|)=\alpha(h)\beta(1)`$、
$`\mathbb E\phi(\psi)=\alpha(1)\beta(\phi)`$ 以及 $`\alpha(1)\beta(1)=1`$，因此
$`\mathbb E[h(\|g\|)\phi(\psi)]=\mathbb Eh(\|g\|)\,\mathbb E\phi(\psi)`$。$`\square`$

<a id="r-7-10"></a>

**引理 7.10.** $`\mathbb E\|g\|^4=N(N+1)`$。若 $`G`$ 是一个 $`a\times b`$ 矩阵，其各元素独立且与 $`g_i`$ 同分布
（*Ginibre* 矩阵），则 $`\mathbb E\mathrm{Tr}(GG^*)^2=ab(a+b)`$。

*证明。* 记 $`g_i=x+iy`$，其中 $`x,y`$ 是独立的 $`N(0,\frac12)`$ 随机变量，则 $`\mathbb Ex^4=3(\mathbb Ex^2)^2=\frac34`$
（定理 A.24：三个配对），所以
$`\mathbb E|g_i|^4=\mathbb Ex^4+2\mathbb Ex^2\mathbb Ey^2+\mathbb Ey^4=\tfrac34+\tfrac12+\tfrac34=2`$，并且
$`\mathbb E\|g\|^4=\sum_i\mathbb E|g_i|^4+\sum_{i\ne j}\mathbb E|g_i|^2\mathbb E|g_j|^2=2N+N(N-1)`$。其次，
$`\mathrm{Tr}(GG^*)^2=\sum_{i,i',j,j'}G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}`$。诸元素及其共轭是同一个中心化
高斯向量（所有元素的实部与虚部）的坐标的复线性组合，且
$`\mathbb E[G_{ij}G_{kl}]=\mathbb E[\bar G_{ij}\bar G_{kl}]=0`$，$`\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}`$
（命题 A.16(a)）。由 Isserlis–Wick 公式（定理 A.24 的复形式，取四个因子的三个配对），
$`\mathbb E[G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}]=\delta_{ii'}\delta_{i'i}+0+\delta_{jj'}\delta_{jj'}=\delta_{ii'}+\delta_{jj'}`$。
对 $`i,i'\le a`$ 和 $`j,j'\le b`$ 求和得 $`ab^2+a^2b`$。$`\square`$

<a id="r-7-11"></a>

**定理 7.11（E. Lubkin, J. Math. Phys. 19, 1028, 1978）。** 对 Haar 随机的
$`\psi\in\mathbb C^a\otimes\mathbb C^b`$ 和 $`\rho_A=\mathrm{Tr}_2|\psi\rangle\langle\psi|`$，
```math
\mathbb E\mathrm{Tr}\,\rho_A^2=\frac{a+b}{ab+1}.
```

*证明。* 令 $`G_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,g\rangle`$，这是一个 Ginibre 矩阵；$`\psi`$ 的矩阵为
$`\Psi=G/\|g\|`$，由命题 7.2(a)，$`(\rho_A)_{\alpha\alpha'}=\sum_\beta\Psi_{\alpha\beta}\bar\Psi_{\alpha'\beta}`$，即
$`\rho_A=\Psi\Psi^*=GG^*/\|g\|^2`$。因此
$`\mathrm{Tr}(GG^*)^2=\|g\|^4\mathrm{Tr}\,\rho_A^2`$，其中 $`\mathrm{Tr}\,\rho_A^2\in[0,1]`$ 是 $`\psi`$ 的函数。
由引理 7.9，$`\mathbb E\mathrm{Tr}(GG^*)^2=\mathbb E\|g\|^4\,\mathbb E\mathrm{Tr}\,\rho_A^2`$，再由引理 7.10，
$`\mathbb E\mathrm{Tr}\,\rho_A^2=ab(a+b)/(ab(ab+1))`$。$`\square`$

第 6 章由酉不变单位向量的四阶矩（引理 6.14）得到同一个均值，作为随机信道的输出纯度
（定理 6.15 及其后的注记）。

<a id="r-7-12"></a>

**推论 7.12.** $`\log\frac{ab+1}{a+b}\le\mathbb E\,S(\rho_A)\le\log\min(a,b)`$。

*证明。* 由推论 3.6，$`S(\rho_A)\le\log a`$；由引理 7.7，
$`S(\rho_A)=S(\mathrm{Tr}_1|\psi\rangle\langle\psi|)\le\log b`$。对凹函数 $`\log`$ 以
$`\rho_A`$ 的非零特征值 $`\lambda_i>0`$ 为权重应用有限形式的 Jensen 不等式（第 3 章习题 2 的解答），得
$`\sum_i\lambda_i\log\lambda_i\le\log\sum_i\lambda_i^2`$，即 $`S(\rho_A)\ge-\log X`$，其中
$`X=\mathrm{Tr}\,\rho_A^2\in(0,1]`$。令 $`m=\mathbb EX=\frac{a+b}{ab+1}`$（定理 7.11）。对 $`(0,\infty)`$ 上的凸函数
$`-\log`$ 在点 $`m`$ 处应用引理 3.4，得 $`-\log X\ge-\log m-(X-m)/m`$；两边都有界（由 Cauchy–Schwarz 不等式，
$`X\ge1/a`$），取期望即得 $`\mathbb E[-\log X]\ge-\log m`$。$`\square`$

当 $`a`$ 固定而 $`b\to\infty`$ 时，上下界趋于一致：随机纯态的一小部分的平均熵趋于其最大值 $`\log a`$。
精确的均值由 Page 公式给出，它在附录 C 中证明。

<a id="r-7-13"></a>

**定理 7.13（Page 公式）。** 对 $`a\le b`$，
```math
\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}.
```
此外 $`\log a-\frac{a-1}{ab+1}-\frac{a-1}{2b}\le\mathbb E\,S(\rho_A)\le\log a-\frac{a-1}{2b}`$，所以
$`\big|\mathbb E\,S(\rho_A)-\log a+\frac a{2b}\big|\le\frac1{2b}`$；对两个量子比特，$`\mathbb E\,S(\rho_A)=\frac13`$。

*证明。* 这就是定理 C.17 连同推论 C.18。附录 C 中的证明从
$`\rho_A=GG^*/\mathrm{Tr}\,GG^*`$（定理 7.11 的证明）出发，计算 Wishart 矩阵 $`GG^*`$ 的特征值的联合密度
（定理 C.13）以及 $`\rho_A`$ 的特征值的联合密度（定理 C.14），再把 $`x\log x`$ 对 Laguerre 能级密度积分
（引理 C.15 和 C.16）。$`\square`$

*计算：* 在 $`\mathbb C^3\otimes\mathbb C^6`$ 的 $`10^4`$ 个随机纯态上检验定理 7.11、7.13 和推论 7.12——[教程，7.5 节](../../tutorial/zh/07_tensor_products.md#s-7-5)。

该公式由 D. N. Page 提出猜想（Phys. Rev. Lett. 71, 1291, 1993），并由 S. K. Foong 和
S. Kanno（Phys. Rev. Lett. 72, 1148, 1994）、J. Sánchez-Ruiz（Phys. Rev. E 52, 5653, 1995）以及
S. Sen（Phys. Rev. Lett. 77, 1, 1996）证明。

当 $`a>b`$ 时，$`a`$ 与 $`b`$ 的角色互换。事实上，从 $`\mathbb C^a\otimes\mathbb C^b`$ 到 $`\mathbb C^b\otimes\mathbb C^a`$ 上的酉映射
$`e_\alpha\otimes e_\beta\mapsto e_\beta\otimes e_\alpha`$ 置换 $`g`$ 的坐标，所以它把 Haar 随机纯态 $`\psi`$ 映为
$`\mathbb C^b\otimes\mathbb C^a`$ 的 Haar 随机纯态 $`\psi'`$；后者的矩阵是 $`\Psi^T`$，所以
$`\mathrm{Tr}_2|\psi'\rangle\langle\psi'|=\Psi^T\bar\Psi=\mathrm{Tr}_1|\psi\rangle\langle\psi|`$（命题 7.2(a)），
由引理 7.7，它的熵是 $`S(\rho_A)`$。因此，记 $`a'=\min(a,b)`$，$`b'=\max(a,b)`$，定理 7.13 表明平均纠缠熵比其最大值
$`\log a'`$ 至少小 $`\frac{a'-1}{2b'}`$，至多小 $`\frac{a'-1}{2b'}+\frac{a'-1}{a'b'+1}`$，而后者小于一个 nat：当 $`b'=1`$ 时它为 $`0`$，
否则 $`\frac{a'-1}{2b'}\le\frac12`$ 且 $`\frac{a'-1}{a'b'+1}<\frac1{b'}\le\frac12`$。所以当 $`ab`$ 固定时，作为 $`a`$ 的函数，
平均值与 $`\log\min(a,b)=\min(\log a,\log(ab)-\log a)`$ 相差不到一个 nat。

$`\rho_A`$ 的分布是**诱导测度**（这一名称来自 K. Życzkowski 和 H.-J. Sommers,
J. Phys. A 34, 7111, 2001；定理 C.14 给出其特征值的密度）。
`states.random_density_matrix(alg, rank=r)` 用相互独立的 $`k_c\times r`$ Ginibre 块 $`G_c`$ 构造
$`G_cG_c^*`$，再除以 $`\sum_c\mathrm{Tr}\,G_cG_c^*`$。对单个块，这就是上面 $`b=r`$ 时 $`\rho_A`$ 的分布
（定理 7.11 的证明）。对多个块，它是 $`H_1\otimes\mathbb C^r`$ 的 Haar 随机纯态 $`\psi=g/\|g\|`$ 的向量态
到 $`\mathcal M_1`$ 的限制的密度的分布。事实上，令 $`G_{i\beta}=\langle e_i\otimes e_\beta,g\rangle`$，这是一个
$`D\times r`$ Ginibre 矩阵（$`D=\sum_ck_c`$），并令 $`G_c`$ 为它位于块 $`c`$ 中的那些行：它们是相互独立的
$`k_c\times r`$ Ginibre 矩阵。取 $`\mathcal M_2=M_r`$（一个扇区）和一个常值融合规则；则只有一个总荷，
$`\Pi=1`$，并且 $`\mathcal M_1\otimes_f\mathcal M_2=B(H_1\otimes\mathbb C^r)`$。由定理 7.6，该限制的密度为
$`E_1(\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|)`$，而 $`\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|=GG^*/\|g\|^2`$
（命题 7.2(a)，与定理 7.11 的证明相同），它的对角块是 $`G_cG_c^*/\|g\|^2`$，其中
$`\|g\|^2=\sum_c\mathrm{Tr}\,G_cG_c^*`$。

*库中：* `states.random_density_matrix(alg, rank=r)`——[教程，7.5 节](../../tutorial/zh/07_tensor_products.md#s-7-5)。

<a id="r-7-14"></a>

**命题 7.14.** 对由 `random_density_matrix(alg, rank=r)` 在 $`\bigoplus_cM_{k_c}`$ 上抽取的 $`\rho`$，
$`D=\sum_ck_c`$，
```math
\mathbb E\mathrm{Tr}\,\rho^2=\frac{\sum_ck_c(k_c+r)}{D(Dr+1)} .
```

*证明。* $`\mathrm{Tr}\,\rho^2=\sum_c\mathrm{Tr}(G_cG_c^*)^2/\|g\|^4`$，其中 $`g\in\mathbb C^{Dr}`$ 是全部矩阵元
构成的向量；它是 $`g/\|g\|`$ 的函数，取值于 $`[0,1]`$。与定理 7.11 一样，引理 7.9 给出
$`\sum_c\mathbb E\mathrm{Tr}(G_cG_c^*)^2=\mathbb E\|g\|^4\,\mathbb E\mathrm{Tr}\,\rho^2`$，而引理 7.10 给出
$`\mathbb E\mathrm{Tr}(G_cG_c^*)^2=k_cr(k_c+r)`$ 和 $`\mathbb E\|g\|^4=Dr(Dr+1)`$，所以
$`\mathbb E\mathrm{Tr}\,\rho^2=\sum_ck_cr(k_c+r)/(Dr(Dr+1))`$。$`\square`$

*计算：* 对 $`M_1\oplus M_2\oplus M_3`$、$`r=2`$ 检验命题 7.14——[教程，7.5 节](../../tutorial/zh/07_tensor_products.md#s-7-5)。

> **物理意义**（物理解释，而非本课程的定理）。随机纯态的平均纠缠熵在总维数固定时作为该部分大小的函数
> 所画出的曲线（定理 7.13 之后描述的曲线）称为 *Page 曲线*。在物理中，它被用来描述信息如何随 Hawking
> 辐射离开黑洞，以及为什么封闭系统的各个部分看起来是热的；这些是解释，而不是本课程的结果。

<a id="s-7-6"></a>

## 7.6 可与之比较的工具

*库中：* 与其他软件的比较——[教程，7.6 节](../../tutorial/zh/07_tensor_products.md#s-7-6)。

<a id="x-7"></a>

## 习题

1. 计算 $`\mathrm{Tr}_2`$ 关于 $`\mathcal M_1\otimes\mathcal M_2`$ 和 $`\mathcal M_1`$ 上的
   $`\tau_{\mathrm{vN}}`$ 的对偶。
2. （次可加性。）设 $`\rho`$ 是 $`\mathcal M_1\otimes\mathcal M_2`$ 上的态，其约化态 $`\rho_1,\rho_2`$ 可逆。
   证明 $`S(\rho)\le S(\rho_1)+S(\rho_2)`$，等号成立当且仅当 $`\rho=\rho_1\otimes\rho_2`$。
3. 设 `site` 为荷为 $`0,1`$ 的 $`\mathbb C\oplus\mathbb C`$（空格点或占据格点），并按加法与自身融合。
   对 $`\psi=(|01\rangle+|10\rangle)/\sqrt2`$（两个格点共享一个粒子）计算 $`E(\psi)`$ 及其两部分，并描述
   $`(|00\rangle+|11\rangle)/\sqrt2`$ 的向量态在融合代数上的限制。

*计算：* 三道习题的数值检验——[教程，习题](../../tutorial/zh/07_tensor_products.md#x-7)。

<details>
<summary>解答</summary>

1. $`\mathcal M_1`$ 上的权重为 $`1/(C_1k_c)`$，积代数上的权重为 $`1/(C_1C_2k_cl_d)`$。对
   $`A\in\mathcal M_1`$，由命题 7.2(a)，
   $`\tau_{\mathrm{vN}}((\mathrm{Tr}_2\,X)A)=\sum_c\frac1{C_1k_c}\sum_d\mathrm{Tr}(X_{cd}(A_c\otimes1))    =\sum_{c,d}\frac1{C_1C_2k_cl_d}\mathrm{Tr}\big(X_{cd}\,C_2l_d(A_c\otimes1)\big)`$，所以对偶为
   $`A\mapsto\bigoplus_{c,d}C_2l_d\,A_c\otimes1_{l_d}=A\otimes Z`$，其中 $`Z=\bigoplus_dC_2l_d1_{l_d}`$ 是
   $`\mathcal M_2`$ 的中心元。它保单位当且仅当 $`Z=1`$，即 $`\mathcal M_2=\mathbb C`$；在本例中不保单位
   （$`Z=4\cdot1_2\oplus6\cdot1_3`$）。对 `Tr_norm`，同样的计算给出 $`A\otimes\bigoplus_dl_d1_{l_d}`$，它保单位当且仅当
   所有 $`l_d=1`$。blunt 迹总是使 $`\mathrm{Tr}_2`$ 成为 $`A\mapsto A\otimes1`$ 的对偶；一般地，由同样的计算，
   $`\mathcal M_1`$ 上的 $`\mathrm{tr}_w`$ 与积代数上的权重 $`w_{cd}=w_c`$ 也能做到这一点。

2. $`\rho_1\otimes\rho_2`$ 是一个态（定理 7.1(c)），可逆，特征向量为 $`u_i\otimes v_j`$，
   特征值为 $`\lambda_i\mu_j`$，所以 $`\log(\rho_1\otimes\rho_2)=\log\rho_1\otimes1+1\otimes\log\rho_2`$
   （定理 1.5(a)）。由偏迹的定义，
   $`\mathrm{Tr}\,\rho(\log\rho_1\otimes1)=\mathrm{Tr}\,\rho_1\log\rho_1`$，对第二个因子同理。
   因此 $`S(\rho\|\rho_1\otimes\rho_2)=S(\rho_1)+S(\rho_2)-S(\rho)`$，它 $`\ge0`$，等号成立当且仅当
   $`\rho=\rho_1\otimes\rho_2`$（推论 3.6）。这个差就是互信息。

3. 扇区为 $`N=0,1,2`$，大小分别为 $`1,2,1`$。$`\psi`$ 位于 $`N=1`$ 中，$`\psi_0=|01\rangle/\sqrt2`$，
   $`\psi_1=|10\rangle/\sqrt2`$，所以 $`p=(\frac12,\frac12)`$；由于各块都是一维的，
   $`E(\hat\psi_c)=0`$：由命题 7.8，$`E(\psi)=H(p)=\log2`$ 全部是数纠缠，构型部分为 $`0`$。
   向量 $`(|00\rangle+|11\rangle)/\sqrt2`$ 是 $`N=0`$ 与 $`N=2`$ 的叠加；在块对角的融合代数上，它的向量态是
   乘积态的混合 $`\frac12|00\rangle\langle00|+\frac12|11\rangle\langle11|`$（第 3 章习题 3）。

</details>

**下一章：** [第 8 章](08_dynamics.md)讨论动力学：带扇区的代数中的 Schrödinger 方程与 Lindblad 方程、
扇区之间的跃迁以及后选择。
