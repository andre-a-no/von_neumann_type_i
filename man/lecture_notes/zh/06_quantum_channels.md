[English](../06_quantum_channels.md) · [Español](../es/06_quantum_channels.md) · [Français](../fr/06_quantum_channels.md) · [Deutsch](../de/06_quantum_channels.md) · **中文** · [日本語](../ja/06_quantum_channels.md) · [Русский](../ru/06_quantum_channels.md)

# 6. 量子信道

量子信道是一个完全正且保迹的线性映射。酉矩阵的共轭作用、到中心上的条件期望、由等距算子 $V$
给出的映射 $\rho\mapsto\mathrm{Tr}_E(V\rho V^*)$，以及把权重从一个扇区移到另一个扇区的映射都是
例子；它们的物理解读（演化、结果被遗忘的测量、环境、粒子的丢失）只在*物理意义*框中讨论。本章对有限维冯·诺依曼代数之间的
这类映射证明其标准描述——Choi 矩阵、Kraus 算子、Stinespring 扩张，然后研究直和情形的特殊之处：
扇区之间的映射、关于第 3 章中各种迹的对偶，以及信道的经典部分。第 2 章和第 3 章已经用到了本库的
信道；这里给出它们的理论。

**你将学到**为什么正性是不够的，如何在 Choi 矩阵、Kraus 算子、超算子与等距算子之间相互转换，
为什么一旦映射改变扇区，"保迹"就依赖于迹的选取，保单位信道使哪些算子保持不变，以及随机信道
平均而言对纯态有什么作用。

本章的数值检验见[教程 6](../../tutorial/zh/06_quantum_channels.md)。

在本章中，$\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$，
$\mathcal N=\bigoplus_{d=1}^{C'}M_{k'_d}(\mathbb C)$，其极小中心投影分别为 $P_c$ 和
$Q_d$；$A=\bigoplus_cA_c$；$M_{k'\times k}$ 是复 $k'\times k$ 矩阵构成的空间；
$E_{ij}=|i\rangle\langle j|$ 为矩阵单位；$\mathrm{Tr}$ 是 blunt 迹
$\sum_c\mathrm{Tr}\,A_c$，而 $\mathrm{tr}_w=\sum_cw_c\mathrm{Tr}\,A_c$ 是一个忠实
迹（定理 3.1）。分块矩阵作用在 $\mathbb C^n\otimes\mathbb C^k$ 上，块指标在前。

<a id="s-6-1"></a>

## 6.1 完全正映射与 Choi 定理

若 $A\ge0$ 蕴含 $\Phi(A)\ge0$，则称线性映射 $\Phi:\mathcal M\to\mathcal N$ 是**正的**。
元素属于 $\mathcal M$ 的 $n\times n$ 矩阵 $[A_{ij}]$ 再次构成矩阵代数的直和，
经由 $[A_{ij}]\mapsto([(A_{ij})_c]_{ij})_c$ 有 $M_n(\mathcal M)\cong\bigoplus_cM_n(M_{k_c})$，
所以 $[A_{ij}]\ge0$ 意味着每个分块矩阵 $[(A_{ij})_c]\in M_{nk_c}$ 都是正的。若
$\Phi_n:[A_{ij}]\mapsto[\Phi(A_{ij})]$ 是正的，则称 $\Phi$ 是 **$n$-正的**；若它对每个 $n$ 都是
$n$-正的，则称它是**完全正的**（CP）。若 $\mathrm{Tr}\,\Phi(A)=\mathrm{Tr}\,A$，则称它是
**保迹的**（TP）；若 $\Phi(1)=1$，则称它是**保单位的**；若它既是 CP 又是 TP，则称它为**信道**。
我们会反复用到：对每个长方矩阵 $T$，$X\ge0$ 蕴含 $TXT^*\ge0$，因为
$\langle y,TXT^*y\rangle=\langle T^*y,XT^*y\rangle$。

<a id="r-6-1"></a>

**命题 6.1.** 对 $k\ge2$，$M_k$ 上的转置 $\theta(A)=A^T$ 是正的，但不是
$2$-正的。

*证明。* 对 $A\ge0$，$\langle x,A^Tx\rangle=\sum_{ij}\bar x_iA_{ji}x_j=\langle\bar x,A\bar x\rangle\ge0$。
取 $\Omega=e_1\otimes e_1+e_2\otimes e_2\in\mathbb C^2\otimes\mathbb C^k$，则
$|\Omega\rangle\langle\Omega|=[E_{ij}]_{i,j=1}^2\ge0$，而 $\theta_2$ 把它映为
$F=[E_{ji}]_{i,j=1}^2$，后者把 $e_1\otimes e_2-e_2\otimes e_1$ 映为它的相反数。$\square$

*计算：* 转置的 Choi 矩阵有特征值 $-1$——[教程，6.1 节](../../tutorial/zh/06_quantum_channels.md#s-6-1)。

<a id="r-6-2"></a>

**定理 6.2（Choi）。** 对线性映射 $\Phi:M_k\to M_{k'}$，以下各条等价：
(a) $\Phi$ 是 CP 的；(b) $\Phi$ 是 $k$-正的；(c) **Choi 矩阵**
$$
J(\Phi)=\sum_{i,j=1}^kE_{ij}\otimes\Phi(E_{ij})=[\Phi(E_{ij})]_{i,j=1}^k\in M_k(M_{k'})
$$
是正的；(d) 存在 $K_1,\dots,K_r\in M_{k'\times k}$（**Kraus 算子**），使得
$\Phi(A)=\sum_lK_lAK_l^*$，并且可以取 $r=\mathrm{rank}\,J(\Phi)\le kk'$。

*证明。* (a)$\Rightarrow$(b) 是平凡的。(b)$\Rightarrow$(c)：对
$\Omega=\sum_ie_i\otimes e_i\in\mathbb C^k\otimes\mathbb C^k$，
$|\Omega\rangle\langle\Omega|=[E_{ij}]\ge0$，且 $J(\Phi)=\Phi_k([E_{ij}])$。(c)$\Rightarrow$(d)：由
定理 1.3，$J(\Phi)=\sum_{l=1}^rv_lv_l^*$，其中 $r=\mathrm{rank}\,J(\Phi)$（$v_l$ 为非零特征值的
特征向量，乘以特征值的平方根）。每个 $v\in\mathbb C^k\otimes\mathbb C^{k'}$
都可唯一地写成 $v=\sum_ie_i\otimes Ke_i$，其中 $K\in M_{k'\times k}$（$Ke_i$ 是 $v$ 的第 $i$ 块），
而 $vv^*$ 的 $(i,j)$ 块为 $Ke_ie_j^*K^*=KE_{ij}K^*$。把 $v_l$ 以这种方式写成对应于
$K_l$ 的形式并比较各块，得到 $\Phi(E_{ij})=\sum_lK_lE_{ij}K_l^*$，由线性即得 (d)。
(d)$\Rightarrow$(a)：$\Phi_n([A_{ij}])=\sum_l(1_n\otimes K_l)[A_{ij}](1_n\otimes K_l)^*\ge0$。$\square$

<a id="r-6-3"></a>

**推论 6.3.** 若 $\Phi(A)=\sum_{l=1}^sL_lAL_l^*$，则 $s\ge\mathrm{rank}\,J(\Phi)$。因此
Kraus 算子的最少个数，即 **Kraus 秩**，等于 $\mathrm{rank}\,J(\Phi)$。

*证明。* 把上面的计算倒过来读，$J(\Phi)=\sum_{l=1}^sw_lw_l^*$，其中
$w_l=\sum_ie_i\otimes L_le_i$，这是 $s$ 个秩 $\le1$ 的矩阵之和。$\square$

*计算：* 随机信道的 Choi 矩阵与 Kraus 秩——[教程，6.1 节](../../tutorial/zh/06_quantum_channels.md#s-6-1)。

直和之间的映射是块与块之间映射构成的矩阵。设 $\iota_c:M_{k_c}\to\mathcal M$
把一个矩阵放入块 $c$，$\pi_d:\mathcal N\to M_{k'_d}$ 取出块 $d$，并令
$\Phi_{dc}=\pi_d\circ\Phi\circ\iota_c:M_{k_c}\to M_{k'_d}$。由于 $A=\sum_c\iota_c(A_c)$，
$\Phi(A)_d=\sum_c\Phi_{dc}(A_c)$。

<a id="r-6-4"></a>

**命题 6.4.** $\Phi:\mathcal M\to\mathcal N$ 是 $n$-正的当且仅当每个 $\Phi_{dc}$ 都是 $n$-正的。因此
$\Phi$ 是 CP 的当且仅当存在 $K^{dc}_l\in M_{k'_d\times k_c}$，$l=1,\dots,r$，使得
$$
\Phi(A)_d=\sum_c\sum_{l=1}^rK^{dc}_l\,A_c\,K^{dc*}_l ,
$$
并且取 $r=\max_{d,c}\mathrm{rank}\,J(\Phi_{dc})$ 就足够了。$\Phi$ 是 TP 的当且仅当对每个 $c$ 有
$\sum_{d,l}K^{dc*}_lK^{dc}_l=1_{k_c}$。

*证明。* 若所有 $\Phi_{dc}$ 都是 $n$-正的且 $[A_{ij}]\ge0$，则每个 $X^c=[(A_{ij})_c]\ge0$，而
$\Phi_n([A_{ij}])$ 的块 $d$ 为 $\sum_c(\Phi_{dc})_n(X^c)\ge0$。反之，对 $M_n(M_{k_c})$ 中的
$[X_{ij}]\ge0$，在 $M_n(\mathcal M)$ 中有 $[\iota_c(X_{ij})]\ge0$（它的块 $c$ 是 $[X_{ij}]$，其余块
为零），而它的像的块 $d$ 是 $[\Phi_{dc}(X_{ij})]$。Kraus 形式由对每一对 $(d,c)$ 应用定理 6.2
得到，再用零算子补齐到共同的 $r$。最后，由迹的循环性，
$\mathrm{Tr}\,\Phi(A)=\sum_c\mathrm{Tr}\big(A_c\sum_{d,l}K^{dc*}_lK^{dc}_l\big)$，它对所有 $A$
都等于 $\sum_c\mathrm{Tr}\,A_c$ 当且仅当每个和都是 $1_{k_c}$（取矩阵单位：
$\mathrm{Tr}(E_{ij}X)=X_{ji}$）。$\square$

线性映射 $\Phi:M_k\to M_{k'}$ 对所有 $X$ 满足 $\Phi(X^*)=\Phi(X)^*$ 当且仅当
$J(\Phi)^*=J(\Phi)$：$J^*$ 的 $(i,j)$ 块是 $\Phi(E_{ji})^*$，由线性性，它对所有 $i,j$ 都等于
$\Phi(E_{ij})$ 当且仅当 $\Phi$ 保持伴随。

*库中：* 保扇区的映射是 `Channel`，命题 6.4 中的一般映射是 `InterSectorChannel`；`Channel.choi()`、`Channel.superoperator()` 和 `Channel.from_superoperator` 在各种表示之间转换——[教程，6.1 节](../../tutorial/zh/06_quantum_channels.md#s-6-1)。

> **物理意义**（物理解释，而非本课程的定理）。系统从来不是孤立的：它可能与一个未被触及的、具有 $n$ 个能级的辅助系统
> 相关联。物理过程在这一对系统上作用为 $\Phi_n=\mathrm{id}_n\otimes\Phi$，并且必须把态映为态，
> 所以 $\Phi$ 必须是完全正的。转置是正的，但它把纠缠态 $|\Omega\rangle\langle\Omega|/2$ 映为
> 一个具有特征值 $-1/2$ 的算子。Choi 矩阵是最大纠缠态的像：单单一个态就决定了完全正性。

<a id="s-6-2"></a>

## 6.2 Stinespring 扩张与 Kraus 算子的自由度

对从 $M_k$ 到 $M_{k'}$ 的 $\Phi(A)=\sum_{l=1}^rK_lAK_l^*$，令 $\Phi^*(B)=\sum_lK_l^*BK_l$；由
循环性 $\mathrm{Tr}(\Phi(A)B)=\mathrm{Tr}(A\Phi^*(B))$，所以 $\Phi^*$ 是 $\Phi$
关于迹的对偶（6.3 节），并且不依赖于 $K_l$ 的选取（若对所有 $A$ 有
$\mathrm{Tr}(AX)=\mathrm{Tr}(AX')$，则 $X=X'$：取 $A=E_{ji}$）。在
$\mathbb C^r\otimes\mathbb C^{k'}$（**环境** $\mathbb C^r$ 在前）上令
$\mathrm{Tr}_E\,X=\sum_l(e_l^*\otimes1)X(e_l\otimes1)$。

<a id="r-6-5"></a>

**定理 6.5（Stinespring，有限维情形）。** 令
$V=\sum_le_l\otimes K_l:\mathbb C^k\to\mathbb C^r\otimes\mathbb C^{k'}$。则
$$
\Phi(A)=\mathrm{Tr}_E(VAV^*),\qquad\Phi^*(B)=V^*(1_r\otimes B)V,
$$
并且 $\Phi$ 是 TP 的当且仅当 $V$ 是等距算子。反之，每个 $V$ 都通过第一个公式定义一个 CP 映射。

*证明。* $(e_l^*\otimes1)V=K_l$，所以 $\mathrm{Tr}_E(VAV^*)=\sum_lK_lAK_l^*$；反之，
$\mathrm{Tr}_E(VAV^*)$ 是以 $K_l=(e_l^*\otimes1)V$ 为算子的 Kraus 形式。其次，
$V^*(1\otimes B)V=\sum_{l,l'}\langle e_l,e_{l'}\rangle K_l^*BK_{l'}=\Phi^*(B)$，并且
$V^*V=\sum_lK_l^*K_l$，它等于 $1$ 当且仅当 $\Phi$ 是 TP 的（命题 6.4）。$\square$

因此，使得对某个 $V:\mathbb C^k\to\mathbb C^r\otimes\mathbb C^{k'}$ 有
$\Phi(A)=\mathrm{Tr}_E(VAV^*)$ 的最小 $r$ 就是 Kraus 秩 $\mathrm{rank}\,J(\Phi)$：每个这样的
$V$ 都给出一个含 $r$ 个算子的 Kraus 形式，所以 $r\ge\mathrm{rank}\,J(\Phi)$（推论 6.3），而由
定理 6.2(d)，$r=\mathrm{rank}\,J(\Phi)$ 是可以达到的。

*计算：* 对一个随机信道验证定理 6.5——[教程，6.2 节](../../tutorial/zh/06_quantum_channels.md#s-6-2)。

带表示 $\pi(B)=1_r\otimes B$ 的公式 $\Phi^*(B)=V^*\pi(B)V$，正是 W. F. Stinespring（Proc. AMS 6 (1955) 211）
陈述该定理时所用的形式，他的设定比这里证明的有限维情形要一般得多（历史归属；延伸阅读：V. Paulsen,
*Completely Bounded Maps and Operator Algebras*, Cambridge 2002, Ch. 4）。对于直和之间的映射，设
$\iota'_d$ 为 $\bigoplus_d\mathbb C^{k'_d}$ 的第 $d$ 个直和项的嵌入，并令
$V_c=\sum_{d,l}e_l\otimes\iota'_dK^{dc}_l:\mathbb C^{k_c}\to\mathbb C^r\otimes\bigoplus_d\mathbb C^{k'_d}$。
则 $(e_l^*\otimes1)V_c=\sum_d\iota'_dK^{dc}_l$，所以
$\sum_c\mathrm{Tr}_E(V_cA_cV_c^*)=\sum_{c,l}\big(\sum_d\iota'_dK^{dc}_l\big)A_c\big(\sum_{d'}\iota'_{d'}K^{d'c}_l\big)^*$
的 $(d,d)$ 块是 $\sum_{c,l}K^{dc}_lA_cK^{dc*}_l=\Phi(A)_d$；又因为 $\iota'^*_d\iota'_{d'}=\delta_{dd'}1$，
$V_c^*V_c=\sum_{d,l}K^{dc*}_lK^{dc}_l$，所以 $\Phi$ 是 TP 的当且仅当每个 $V_c$ 都是等距算子
（命题 6.4）。`random_channel` 和 `random_inter_sector_channel` 正是这样抽取信道的
（6.5 节）。

*库中：* `channels.random_channel` 和 `channels.random_inter_sector_channel`——[教程，6.5 节](../../tutorial/zh/06_quantum_channels.md#s-6-5)。

<a id="r-6-6"></a>

**引理 6.6.** 若 $X,Y\in M_{p\times n}$ 且 $XX^*=YY^*$，则对某个酉矩阵 $W\in M_n$ 有 $Y=XW$。

*证明。* 对所有 $z$，$\|X^*z\|^2=\langle z,XX^*z\rangle=\|Y^*z\|^2$，所以 $X^*z\mapsto Y^*z$ 是从
$X^*$ 的值域到 $Y^*$ 的值域上的一个良定义的线性等距映射。与定理 1.8 的证明一样，它可以延拓为
满足 $UX^*=Y^*$ 的酉矩阵 $U$；取 $W=U^*$。$\square$

<a id="r-6-7"></a>

**定理 6.7（酉自由度）。** 设 $K_1,\dots,K_n$ 与 $L_1,\dots,L_n$ 属于 $M_{k'\times k}$
（较短的一族用零补齐）。则对所有 $A$ 有 $\sum_lK_lAK_l^*=\sum_jL_jAL_j^*$ 当且仅当对某个酉矩阵
$u\in M_n$ 有 $L_j=\sum_lu_{jl}K_l$。

*证明。* 若 $L_j=\sum_lu_{jl}K_l$，则
$\sum_jL_jAL_j^*=\sum_{l,l'}(u^*u)_{l'l}K_lAK_{l'}^*=\sum_lK_lAK_l^*$。反之，相等的映射具有
相等的 Choi 矩阵，即 $\sum_lv_lv_l^*=\sum_jw_jw_j^*$，其中 $v_l=\sum_ie_i\otimes K_le_i$，
$w_j=\sum_ie_i\otimes L_je_i$。对以 $v_l$ 和 $w_j$ 为列的矩阵应用引理 6.6，得到
$w_j=\sum_lW_{lj}v_l$；由于 $K\mapsto\sum_ie_i\otimes Ke_i$ 是线性单射，
$L_j=\sum_lW_{lj}K_l$，而 $u=W^T$ 是酉的。$\square$

*计算：* $n=3$ 时的定理 6.7——[教程，6.2 节](../../tutorial/zh/06_quantum_channels.md#s-6-2)。

> **物理意义**（物理解释，而非本课程的定理）。每个信道都是系统与处于固定态的环境之间的
> 一次等距相互作用，随后忽略环境。
> 在基 $e_l$ 下测量环境，以概率 $\mathrm{Tr}\,K_l\rho K_l^*$ 得到结果 $l$，并留下态
> $K_l\rho K_l^*/\mathrm{Tr}(K_l\rho K_l^*)$：Kraus 算子就是可能发生的"量子跳跃"。
> 换一组基会得到别的跳跃，但信道相同（定理 6.7），而最小的环境的维数
> 等于 Kraus 秩，这已在定理 6.5 之后证明。

<a id="s-6-3"></a>

## 6.3 关于迹的对偶

**定义。** 设 $\mathrm{tr}_w$、$\mathrm{tr}_v$ 分别是 $\mathcal M$、$\mathcal N$ 上的
忠实迹。线性映射 $\Phi:\mathcal M\to\mathcal N$ 关于 $(w,v)$ 的**对偶**是满足
$$
\mathrm{tr}_v\big(\Phi(x)\,y\big)=\mathrm{tr}_w\big(x\,\Phi^\sharp(y)\big)\qquad(x\in\mathcal M,\ y\in\mathcal N).
$$
的线性映射 $\Phi^\sharp:\mathcal N\to\mathcal M$。它存在且唯一：对固定的 $y$，
$x\mapsto\mathrm{tr}_v(\Phi(x)y)$ 是一个线性泛函，所以由定理 3.2，它对唯一的
$\rho=:\Phi^\sharp(y)$ 等于 $\mathrm{tr}_w(\rho x)=\mathrm{tr}_w(x\rho)$，而唯一性使得
$\Phi^\sharp$ 是线性的。对 blunt 迹，我们记作 $\Phi^*$。
令 $z_w=\sum_cw_cP_c$，$z_v=\sum_dv_dQ_d$（它们是中心元且可逆）。

<a id="r-6-8"></a>

**命题 6.8.** 若 $\Phi$ 具有 Kraus 算子 $K^{dc}_l$（命题 6.4），则
$$
\Phi^\sharp(y)_c=\sum_d\frac{v_d}{w_c}\sum_lK^{dc*}_l\,y_d\,K^{dc}_l,\qquad\text{i.e.}\qquad
\Phi^\sharp(y)=z_w^{-1}\,\Phi^*(z_v\,y).
$$
特别地，$\Phi^\sharp$ 是 CP 的，其 Kraus 算子为 $(v_d/w_c)^{1/2}K^{dc*}_l$。

*证明。* 由每个块中迹的循环性，
$$
\mathrm{tr}_v(\Phi(x)y)=\sum_dv_d\sum_{c,l}\mathrm{Tr}\big(K^{dc}_lx_cK^{dc*}_ly_d\big)
=\sum_cw_c\mathrm{Tr}\Big(x_c\sum_d\frac{v_d}{w_c}\sum_lK^{dc*}_ly_dK^{dc}_l\Big),
$$
由唯一性即得该公式。对 $v=w=1$，它就是 $\Phi^*(y)_c=\sum_{d,l}K^{dc*}_ly_dK^{dc}_l$；
因子 $v_d$ 和 $1/w_c$ 就是乘以 $z_v$ 和 $z_w^{-1}$。$\square$

若对所有 $x$ 有 $\mathrm{tr}_v\,\Phi(x)=\mathrm{tr}_w(x)$，则称 $\Phi$ **关于 $(w,v)$ 保迹**。

<a id="r-6-9"></a>

**推论 6.9.** (a) $\Phi$ 关于 $(w,v)$ 保迹当且仅当 $\Phi^\sharp(1)=1$，即对每个 $c$ 有
$\sum_{d,l}(v_d/w_c)K^{dc*}_lK^{dc}_l=1_{k_c}$。(b) $\Phi$ 是保单位的当且仅当 $\Phi^\sharp$
关于 $(v,w)$ 保迹。(c) $\Phi^\sharp$ 关于 $(v,w)$ 的对偶是 $\Phi$。(d) 若
$\Phi:\mathcal M\to\mathcal M$ 是保扇区的且 $v=w$，则对每个 $w$ 有 $\Phi^\sharp=\Phi^*$，
并且 $\Phi$ 对一个迹保迹当且仅当它对所有迹都保迹。

*证明。* (a) $\mathrm{tr}_v\,\Phi(x)=\mathrm{tr}_w(x\Phi^\sharp(1))$，而
$\mathrm{tr}_w(x(\Phi^\sharp(1)-1))=0$ 对所有 $x$ 成立当且仅当 $\Phi^\sharp(1)=1$（取
$x=(\Phi^\sharp(1)-1)^*$；$\mathrm{tr}_w$ 是忠实的）；然后应用命题 6.8。(b)
$\mathrm{tr}_w\,\Phi^\sharp(y)=\mathrm{tr}_v(\Phi(1)y)$ 对所有 $y$ 都等于 $\mathrm{tr}_v(y)$
当且仅当 $\Phi(1)=1$。(c) 由迹的性质，$\mathrm{tr}_w(\Phi^\sharp(y)x)=\mathrm{tr}_v(y\Phi(x))$；
应用唯一性。(d) 只出现 $d=c$ 的对，此时 $v_d/w_c=1$。$\square$

**Heisenberg 绘景中的密度。** 设 $\Psi:\mathcal N\to\mathcal M$ 是线性映射，$\omega$ 是
$\mathcal M$ 上关于 $\mathrm{tr}_w$ 具有密度 $\rho$ 的线性泛函（定理 3.2），
$\Psi^\sharp:\mathcal M\to\mathcal N$ 是 $\Psi$ 关于 $(v,w)$ 的对偶。则
$\omega(\Psi(y))=\mathrm{tr}_w(\Psi(y)\rho)=\mathrm{tr}_v(y\Psi^\sharp(\rho))=\mathrm{tr}_v(\Psi^\sharp(\rho)y)$，
所以 $\omega\circ\Psi$ 关于 $\mathrm{tr}_v$ 的密度为 $\Psi^\sharp(\rho)$；若 $\Psi$ 是保单位的，则由
推论 6.9(b)（交换 $\mathcal M$ 与 $\mathcal N$ 的角色），$\Psi^\sharp$ 关于 $(w,v)$ 保迹。

对于改变扇区的映射，推论 6.9(d) 不成立。设 $\mathcal M=\mathcal N=M_1\oplus M_2$（扇区 $0$ 和 $1$），
并设 $\Phi$ 是称为**粒子丢失**的映射（见下面的框），
$\Phi(x)=(x_0+\mathrm{Tr}\,x_1)\oplus0$，其 Kraus 算子为 $K^{00}=1$，$K^{01}_l=\langle l|$。
它对 blunt 迹是保迹的。对 $\mathrm{tr}_w$，$w=v=(1,\tfrac12)$，有 $\Phi^\sharp(1)=1\oplus2\cdot1_2$，
且 $\mathrm{tr}_v\,\Phi(x)=x_0+\mathrm{Tr}\,x_1\ne x_0+\tfrac12\mathrm{Tr}\,x_1=\mathrm{tr}_w(x)$。

*库中：* `Channel.adjoint()` 是 $\Phi^*$，`InterSectorChannel.adjoint(trace=...)` 是 $\Phi^\sharp$；命题 6.8、推论 6.9(a) 以及粒子丢失的检验——[教程，6.3 节](../../tutorial/zh/06_quantum_channels.md#s-6-3)。

> **物理意义**（物理解释，而非本课程的定理）。与迹无关的对象是作用在可观测量上的 Heisenberg 映射
> $\Psi:\mathcal N\to\mathcal M$：过程之前的态 $\omega$ 在过程之后变为
> $\omega\circ\Psi$，而 $\Psi$ 是 CP 且保单位的。密度依赖于迹
> （定理 3.2），而密度按 $\Psi$ 关于所选迹的对偶变换（推论 6.9 之后的段落）。如果代数以一定
> 重数作用，物理的迹的权重为 $w_c=m_c$（2.3 节）。在上面的例子中，扇区 $0$ 被解读为"没有粒子"，
> 扇区 $1$ 被解读为"两个格点上的一个粒子"；于是该映射移走了这个粒子。它对数态个数的 blunt 迹
> 保持概率守恒，而对平等地看待各扇区的 `Tr_norm` 则不然。

<a id="s-6-4"></a>

## 6.4 退相干

本库的标准信道在每个块内部作用（$p,\gamma\in[0,1]$）：

| 信道 | $\Phi(\rho)_c$ | 块 $c$ 中的 Kraus 算子 |
|---|---|---|
| `dephasing_channel` | $(1-p)\rho_c+p\mathrm{diag}\,\rho_c$ | $\sqrt{1-p}\,1$, $\sqrt p\,E_{ii}$ |
| `depolarizing_channel` | $(1-p)\rho_c+p\mathrm{Tr}(\rho_c)\,1/k_c$ | $\sqrt{1-p}\,1$, $\sqrt{p/k_c}\,E_{ij}$ |
| `center_expectation` | $\mathrm{Tr}(\rho_c)\,1/k_c$ | $\sqrt{1/k_c}\,E_{ij}$ |
| `amplitude_damping_channel` | 能级 $i>0$ 衰减到 $0$ | $E_{00}+\sqrt{1-\gamma}\sum_{i>0}E_{ii}$, $\sqrt\gamma\,E_{0i}$ |
| `random_mixed_unitary_channel` | $\sum_jp_jU_j\rho U_j^*$ | $\sqrt{p_j}\,U_j$ |

这些 Kraus 形式可以直接计算验证，例如
$\sum_{ij}\frac1kE_{ij}\rho E_{ji}=\frac1k\sum_{ij}\rho_{jj}E_{ii}$。这五个信道都是 TP 的，并且除
$\gamma>0$ 时的振幅阻尼外都是保单位的（命题 6.4：计算 $\sum_lK_l^*K_l$ 和
$\sum_lK_lK_l^*$；对大小为 $k$ 的块中的振幅阻尼，它们分别是 $1$ 和
$1+(k-1)\gamma E_{00}-\gamma\sum_{i>0}E_{ii}$，而对混合酉信道有 $\sum_jp_j=1$）。
具体地，振幅阻尼的作用为 $\rho_{00}\mapsto\rho_{00}+\gamma\sum_{i>0}\rho_{ii}$，
$\rho_{ii}\mapsto(1-\gamma)\rho_{ii}$，对两个都为正的 $i\ne j$ 有 $\rho_{ij}\mapsto(1-\gamma)\rho_{ij}$，
对 $j>0$ 有 $\rho_{0j}\mapsto\sqrt{1-\gamma}\,\rho_{0j}$，$\rho_{j0}\mapsto\sqrt{1-\gamma}\,\rho_{j0}$。到中心上的条件期望 $E$
（命题 2.8）就是 $p=1$ 的退极化。与扇区的测量
$X\mapsto\sum_cP_cXP_c$（它在 $\mathcal M$ 上是恒等映射，见 3.3 节）不同，$E$ 只保留 $p_c$。

*库中：* 表中的五个信道都是 `channels` 中的函数——[教程，6.4 节](../../tutorial/zh/06_quantum_channels.md#s-6-4)。

<a id="r-6-10"></a>

**定理 6.10（保单位信道的不动点）。** 设 $\Phi(X)=\sum_lK_lXK_l^*$，$K_l\in\mathcal M$，
是 TP 且保单位的。则 $\Phi(X)=X$ 当且仅当 $X$ 与每个 $K_l$ 和 $K_l^*$ 可交换。特别地，
不动点构成 $\mathcal M$ 的一个冯·诺依曼子代数。

*证明。* 若 $X$ 与所有 $K_l$ 可交换，则 $\Phi(X)=X\sum_lK_lK_l^*=X$。反之，$\Phi(X^*)=\Phi(X)^*$，
所以不动点的自伴部分 $(X+X^*)/2$ 和 $(X-X^*)/2i$ 也是不动点；不妨设
$X=X^*$。令 $Y_l=K_lX-XK_l$，则
$$
\sum_lY_lY_l^*=\sum_l\big(K_lX^2K_l^*-K_lXK_l^*X-XK_lXK_l^*+XK_lK_l^*X\big)
=\Phi(X^2)-\Phi(X)X-X\Phi(X)+X\Phi(1)X=\Phi(X^2)-X^2,
$$
这里用到了 $Y_l^*=XK_l^*-K_l^*X$、$\Phi(X)=X$ 和 $\Phi(1)=1$。由于 $\Phi$ 是 TP 的，它的迹为 $0$。因为
$\mathrm{Tr}\,Y_lY_l^*=\|Y_l\|_2^2\ge0$，每个 $Y_l=0$；取伴随可知 $X$ 也与 $K_l^*$ 可交换。这些 $X$ 构成一个自伴集合的交换子与 $\mathcal M$
的交：这是一个包含 $1$ 的 $*$-代数，因而是冯·诺依曼代数
（定理 2.2）。$\square$

因此，退相位（$p>0$）恰好保留每个块的对角部分，而退极化（$p>0$）和 $E$
恰好保留中心。没有保单位性，该定理就不成立：振幅阻尼固定 $E_{00}$，
而 $E_{00}$ 与跳跃算子 $E_{01}$ 不可交换。对 $0<\gamma\le1$ 和 $k\ge2$，它在一个块中的不动点空间
与其 Kraus 算子的交换子维数都是 $1$，但它们分别是 $\mathbb CE_{00}$ 和
$\mathbb C1$：由上面的显式形式，不动点满足对 $i\ne j$ 有 $\rho_{ij}=0$，对 $i>0$ 有
$\rho_{ii}=0$（这些元素被乘以 $\sqrt{1-\gamma}\ne1$ 或 $1-\gamma\ne1$），而 $E_{00}$ 的每个倍数都是
不动点；另一方面，与所有 $E_{0i}$、$E_{i0}$（$i>0$）可交换的元素也与它们的乘积 $E_{ij}=E_{i0}E_{0j}$
和 $E_{00}=E_{01}E_{10}$ 可交换，从而与 $M_k$ 可交换，所以它是 $1$ 的倍数
（第 2 章习题 1）。

*计算：* 五个信道的不动点与交换子——[教程，6.4 节](../../tutorial/zh/06_quantum_channels.md#s-6-4)。

<a id="r-6-11"></a>

**定理 6.11（保单位信道使熵增加）。** 设 $\Phi:\mathcal M\to\mathcal N$ 是正的、
TP 且保单位的，$\rho$ 和 $\Phi(\rho)$ 是态，其特征值分别为 $\lambda_i$ 和 $\mu_j$。
则对某个双随机矩阵 $\Delta$ 有 $\mu=\Delta\lambda$，从而对 $[0,1]$ 上的每个凹函数 $f$ 有
$\sum_jf(\mu_j)\ge\sum_if(\lambda_i)$。特别地，$S(\Phi(\rho))\ge S(\rho)$，且
$\mathrm{Tr}\,\Phi(\rho)^2\le\mathrm{Tr}\,\rho^2$。

*证明。* 对每个块对角化：$\rho=\sum_i\lambda_i|e_i\rangle\langle e_i|$，
$\Phi(\rho)=\sum_j\mu_j|f_j\rangle\langle f_j|$，其中 $\bigoplus_c\mathbb C^{k_c}$
与 $\bigoplus_d\mathbb C^{k'_d}$ 的规范正交基中每个向量都只位于一个直和项内，从而
$|e_i\rangle\langle e_i|\in\mathcal M$。于是 $\mu_j=\langle f_j,\Phi(\rho)f_j\rangle=\sum_i\Delta_{ji}\lambda_i$，
其中 $\Delta_{ji}=\langle f_j,\Phi(|e_i\rangle\langle e_i|)f_j\rangle\ge0$，
$\sum_j\Delta_{ji}=\mathrm{Tr}\,\Phi(|e_i\rangle\langle e_i|)=1$，且
$\sum_i\Delta_{ji}=\langle f_j,\Phi(1)f_j\rangle=1$（所以 $\Delta$ 是方阵：对所有元素求和即知）。由
有限 Jensen 不等式 $f(\sum_i\Delta_{ji}\lambda_i)\ge\sum_i\Delta_{ji}f(\lambda_i)$（它可由凹性的定义
对项数归纳得到，同第 3 章习题 2 的解答），
$\sum_jf(\mu_j)\ge\sum_j\sum_i\Delta_{ji}f(\lambda_i)=\sum_if(\lambda_i)$。取 $f(t)=-t\log t$
（$f(0)=0$）和 $f(t)=-t^2$；由引理 3.4，它们在 $[0,1]$ 上是凹的：$t\log t$ 和 $t^2$ 在 $[0,1]$ 上
连续，且在 $(0,1)$ 上的二阶导数分别为 $1/t>0$ 和 $2>0$。$\square$

*计算：* 保单位信道下的熵与纯度，以及振幅阻尼下熵的下降——[教程，6.4 节](../../tutorial/zh/06_quantum_channels.md#s-6-4)。

扇区之间映射的经典部分是它的**转移矩阵**
$T_{dc}(\rho)=\mathrm{Tr}\,\Phi_{dc}(\rho_c)/\mathrm{Tr}\,\rho_c$（$\rho_c\ne0$）；本库的
`transition_matrix()` 是 $T=T(1)$，对应于每个输入扇区的最大混合态。

<a id="r-6-12"></a>

**命题 6.12.** 设 $\Phi:\mathcal M\to\mathcal N$ 是正的且 TP 的，$\rho$ 是满足
$p_c=\mathrm{Tr}\,\rho_c>0$ 的态，$q_d=\mathrm{Tr}\,\Phi(\rho)_d$。则 $T(\rho)$ 是列随机矩阵，
且 $q=T(\rho)p$。记 $E_{\mathcal M}$、$E_{\mathcal N}$ 为到各自中心上的条件期望（命题 2.8），则
$$
E_{\mathcal N}\Phi E_{\mathcal M}(\rho)=\bigoplus_d\frac{(Tp)_d}{k'_d}\,1_{k'_d},\qquad T=T(1).
$$

*证明。* 由于 $\Phi_{dc}$ 是正的，$T_{dc}\ge0$，并且
$\sum_dT_{dc}=\mathrm{Tr}\,\Phi(\iota_c\rho_c)/\mathrm{Tr}\,\rho_c=1$。其次，
$q_d=\sum_c\mathrm{Tr}\,\Phi_{dc}(\rho_c)=\sum_cT_{dc}(\rho)p_c$。最后，
$E_{\mathcal M}\rho=\bigoplus_cp_c1_{k_c}/k_c$，所以
$\Phi(E_{\mathcal M}\rho)_d=\sum_cp_c\Phi_{dc}(1_{k_c}/k_c)$ 的迹为 $\sum_cT_{dc}p_c$，而 $E_{\mathcal N}$
把块 $d$ 替换为它的迹乘以 $1/k'_d$。$\square$

*库中：* 转移矩阵是 `InterSectorChannel.transition_matrix()`；命题 6.12 的检验——[教程，6.4 节](../../tutorial/zh/06_quantum_channels.md#s-6-4)。

**量子比特的相干项。** 在 $M_2$ 上，设 $D_p$ 为退相位，$A_\gamma$ 为振幅阻尼。由上表和上面的
显式形式，$D_p$ 把 $\rho_{01}$ 和 $\rho_{10}$ 乘以 $1-p$ 并保持对角部分不变，而 $A_\gamma$ 把
$\rho_{01}$ 和 $\rho_{10}$ 乘以 $\sqrt{1-\gamma}$，并把对角部分映为
$(\rho_{00}+\gamma\rho_{11},(1-\gamma)\rho_{11})$，与相干项无关。因此
$D_pA_\gamma=A_\gamma D_p$：二者都把相干项乘以 $(1-p)\sqrt{1-\gamma}$，并在对角部分上
像 $A_\gamma$ 一样作用。特别地，对 $T_1,T_\varphi>0$，$t\ge0$，$1-p=e^{-t/T_\varphi}$ 和 $1-\gamma=e^{-t/T_1}$，
该复合映射把相干项乘以 $e^{-t/T_2}$，其中
$$
\frac1{T_2}=\frac1{2T_1}+\frac1{T_\varphi},
$$
并且这两族映射关于 $t$ 可加地复合：$D_{p(t)}D_{p(s)}=D_{p(t+s)}$，
$A_{\gamma(t)}A_{\gamma(s)}=A_{\gamma(t+s)}$，因为因子 $e^{-t/T}$ 相乘，且迹保持不变。

*计算：* 相干因子与 $1/T_2=1/(2T_1)+1/T_\varphi$——[教程，6.4 节](../../tutorial/zh/06_quantum_channels.md#s-6-4)。

> **物理意义**（物理解释，而非本课程的定理）。对一个量子比特，
> 退相位模拟时间常数为 $T_\varphi$ 的纯退相位，振幅阻尼模拟激发能级以时间常数 $T_1$ 的衰减；
> 于是上面的计算给出同时受这两者作用的量子比特的相干时间 $T_2$。定理 6.10 是无退相干信息的
> 代数形式：与所有跳跃算子可交换的东西得以幸存（第 10 章）。
> 命题 6.12 表明，透过中心来看，扇区之间的映射是荷上的一个经典 Markov
> 矩阵。

<a id="s-6-5"></a>

## 6.5 随机信道

与本库（`random_isometry`）一样，**Haar 随机等距** $V:\mathbb C^k\to\mathbb C^N$ 定义为
$G=QR$ 中的因子 $Q$，其中 $G$ 是 $N\times k$ 的 Ginibre 矩阵（各元素独立，密度为
$\pi^{-1}e^{-|z|^2}$），$Q^*Q=1$，$R$ 是对角元为正的上三角矩阵（在 CPU 上，本库先取 LAPACK
分解 $G=Q'R'$，再把 $Q'$ 的第 $i$ 列乘以 $R'_{ii}/|R'_{ii}|$；与定理 5.3 中一样，只要所有
$R'_{ii}\ne0$，即只要 $G$ 的秩为 $k$，这就给出上述分解，而由引理 6.13 的证明，这几乎必然成立；
对 GPU 上的大批量，本库使用 Gram–Schmidt 正交化，它直接给出该分解）。
`random_channel(alg, r)` 在每个块中独立地抽取这样的 $V_c:\mathbb C^{k_c}\to\mathbb C^r\otimes\mathbb C^{k_c}$，
并使用定理 6.5 中的 Kraus 算子。

*库中：* `random_isometry` 和 `channels.random_channel`——[教程，6.5 节](../../tutorial/zh/06_quantum_channels.md#s-6-5)。

<a id="r-6-13"></a>

**引理 6.13.** $V$ 与 $U(N)$ 的 Haar 酉矩阵的前 $k$ 列同分布。
特别地，对每个酉矩阵 $U\in M_N$，$UV$ 与 $V$ 同分布。

*证明。* 设 $Z$ 是一个 $N\times N$ 的 Ginibre 矩阵，其前 $k$ 列为 $G$，并设 $E\in M_{N\times k}$
为 $1_N$ 的前 $k$ 列。由定理 5.3，$Z$ 几乎必然可逆，且 $Z=Q_ZR_Z$，其中
$Q_Z$ 服从 Haar 分布。由于 $R_Z$ 是上三角的，$G=ZE=(Q_ZE)(E^*R_ZE)$ 是所要求类型的
分解，而且是唯一的：若 $G=QR=Q'R'$，则 $R^*R=G^*G=R'^*R'$，所以 $R'R^{-1}$
是酉的、上三角的且对角元为正，因而等于 $1$（引理 5.2 的证明）。所以 $V=Q_ZE$，
并且 $UV=(UQ_Z)E$ 与 $V$ 同分布，因为 $UQ_Z$ 与 $Q_Z$ 同分布。$\square$

<a id="r-6-14"></a>

**引理 6.14.** 设 $u$ 是 $\mathbb C^N$ 中的随机单位向量，其分布在所有
酉变换下不变。则 $\mathbb E[u_a\bar u_b]=\delta_{ab}/N$，且
$$
\mathbb E\big[u_a\bar u_b\,u_c\bar u_d\big]=\frac{\delta_{ab}\delta_{cd}+\delta_{ad}\delta_{cb}}{N(N+1)} .
$$

*证明。* 在 $\mathrm{diag}(e^{i\theta_1},\dots,e^{i\theta_N})$ 下的不变性把这两个矩分别乘以
$e^{i(\theta_a-\theta_b)}$ 和 $e^{i(\theta_a-\theta_b+\theta_c-\theta_d)}$，所以除非 $a=b$，
或相应地 $\{a,c\}=\{b,d\}$（作为多重集），它们都等于零。由置换矩阵下的不变性，
$\mathbb E|u_a|^2=1/N$，并且只剩下 $\alpha=\mathbb E|u_1|^4$ 和 $\beta=\mathbb E|u_1|^2|u_2|^2$。
首先，$\sum_c|u_c|^2=1$ 给出 $\alpha+(N-1)\beta=\mathbb E|u_1|^2=1/N$。其次，$(u_1+u_2)/\sqrt2$ 是
某个酉矩阵 $U$ 下 $Uu$ 的第一个坐标，所以 $\mathbb E|u_1+u_2|^4=4\alpha$。令
$s=|u_1|^2+|u_2|^2$，$z=u_1\bar u_2$，则 $|u_1+u_2|^4=s^2+2s(z+\bar z)+(z+\bar z)^2$；
不变性 $u_1\mapsto-u_1$ 和 $u_1\mapsto iu_1$ 给出 $\mathbb E[s(z+\bar z)]=0$ 和
$\mathbb Ez^2=0$，所以 $4\alpha=(2\alpha+2\beta)+2\beta$，即 $\alpha=2\beta$。因此
$\beta=1/(N(N+1))$，$\alpha=2\beta$，这在所有情形下都给出该公式。（这里用到了 $N\ge2$；对
$N=1$，$|u_1|=1$，两边都等于 $1$。）$\square$

<a id="r-6-15"></a>

**定理 6.15（随机信道的输出纯度）。** 设 $\psi$ 是大小为 $k$ 的扇区中的纯态，
$\Phi$ 由 `random_channel` 以 Kraus 秩 $r$ 抽取。则
$$
\mathbb E\mathrm{Tr}\,\Phi(\psi)^2=\frac{k+r}{kr+1},\qquad\mathbb E\,\Phi(\psi)=\frac1k1_k .
$$

*证明。* 在 $\psi$ 所在的扇区中，$\Phi(\psi)=\mathrm{Tr}_E|u\rangle\langle u|$，其中 $u=V\psi$
是单位向量，由引理 6.13，其分布是酉不变的。用指标对 $(e,x)$ 标记
$\mathbb C^r\otimes\mathbb C^k$。则 $\Phi(\psi)_{xy}=\sum_eu_{ex}\bar u_{ey}$，其均值为 $r\delta_{xy}/(kr)$，且
$\mathrm{Tr}\,\Phi(\psi)^2=\sum_{x,y,e,f}u_{ex}\bar u_{ey}u_{fy}\bar u_{fx}$。由引理 6.14（取
$N=kr$），第一项迫使 $x=y$（$kr^2$ 个指标组），第二项迫使 $e=f$（$k^2r$ 个指标组），所以
均值为 $(kr^2+k^2r)/(kr(kr+1))$。$\square$

*计算：* 对 $r=1,2,4$ 用蒙特卡罗方法检验定理 6.15——[教程，6.5 节](../../tutorial/zh/06_quantum_channels.md#s-6-5)。

对 $r=1$，信道是酉的，纯度为 $1$；当 $r\to\infty$ 时，纯度趋于 $1/k$，即
最大混合态的纯度。证明只用到了 $u$ 的分布的酉不变性，
因此它也给出了 Lubkin 关于随机纯态 $u=g/\|g\|$（$g$ 为高斯向量，由定理 A.20(b)，其分布是酉不变的）
的约化态的平均纯度，第 7 章定理 7.11 用高斯积分证明了这一结果；引理 6.14 中的数
$\alpha,\beta$ 就是推论 5.11 中关于 Haar 酉矩阵元素的那些数，而推论 A.21 对 $u=g/\|g\|$
重新计算了它们。

> **物理意义**（物理解释，而非本课程的定理）。随机信道模拟来自维数为 $r$ 的环境的一般的、
> 无结构的噪声。仅 $r=4$ 就使一个纯的 qutrit 态的纯度平均降到 $7/13\approx0.54$（定理 6.15），
> 而无噪声时为 $1$，最小值为 $1/3$（由 Cauchy–Schwarz 不等式，
> $\mathrm{Tr}\,\rho^2\ge(\mathrm{Tr}\,\rho)^2/k$）。在超选择系统中，噪声在每个扇区中
> 分别作用，相关的维数是扇区的维数，而不是整个空间的维数。

<a id="s-6-6"></a>

## 6.6 可与之比较的工具

*库中：* 与其他软件的比较——[教程，6.6 节](../../tutorial/zh/06_quantum_channels.md#s-6-6)。

<a id="x-6"></a>

## 习题

1. （Werner–Holevo）证明 $M_k$（$k\ge2$）上的 $\Phi(\rho)=\big(\mathrm{Tr}(\rho)1-\rho^T\big)/(k-1)$
   是一个 Kraus 秩为 $k(k-1)/2$ 的保单位信道，尽管它是由
   转置构造的。
2. 证明对 $0<p\le1$，$M_k$ 上的退相位信道的 Kraus 秩为 $k$，比本库存储的
   Kraus 算子个数少一。
3. 设 $X$ 是 $\mathbb C^k$ 上的循环移位，$Z=\mathrm{diag}(1,\omega,\dots,\omega^{k-1})$，
   $\omega=e^{2\pi i/k}$。证明 $\frac1{k^2}\sum_{a,b=0}^{k-1}X^aZ^b\rho Z^{-b}X^{-a}=\mathrm{Tr}(\rho)1/k$：
   到 $M_k$ 的中心上的条件期望是一个混合酉信道。
4. 证明当中间的代数上使用同一个迹时，$(\Psi\circ\Phi)^\sharp=\Phi^\sharp\circ\Psi^\sharp$。

*计算：* 四道习题的数值检验——[教程，习题](../../tutorial/zh/06_quantum_channels.md#x-6)。

<details>
<summary>解答</summary>

1. 记 $F$ 为交换算子，则 $J(\theta)=\sum_{ij}E_{ij}\otimes E_{ji}=F$，$J(\rho\mapsto\mathrm{Tr}(\rho)1)=1\otimes1$，
   所以 $J(\Phi)=(1-F)/(k-1)=2P_-/(k-1)$，其中 $P_-=(1-F)/2$ 是到
   反对称向量上的投影（$F=F^*$，$F^2=1$），其秩为 $\mathrm{Tr}\,P_-=(k^2-\mathrm{Tr}\,F)/2$，
   其中 $\mathrm{Tr}\,F=\sum_{ij}\langle e_i\otimes e_j,e_j\otimes e_i\rangle=k$。所以 $J\ge0$，秩为 $k(k-1)/2$（定理 6.2、
   推论 6.3），并且可直接验证 $\mathrm{Tr}\,\Phi(\rho)=\mathrm{Tr}\,\rho$，$\Phi(1)=1$。

2. $J=(1-p)|\Omega\rangle\langle\Omega|+p\sum_iE_{ii}\otimes E_{ii}$ 在
   $e_i\otimes e_i$ 张成的空间之外为零，而在这个 $k$ 维空间上它是矩阵 $(1-p)\mathbf 1\mathbf 1^T+p1$
   （$\mathbf 1$ 为全一向量），其特征值为 $p$ 和 $(1-p)k+p$，都是正的。所以
   $\mathrm{rank}\,J=k$，而本库存储 $1+k_{\max}$ 个算子。

3. $\mathrm{Tr}((X^aZ^b)^*X^{a'}Z^{b'})=k\,\delta_{aa'}\delta_{bb'}$：除非 $a=a'$，
   $Z^{-b}X^{a'-a}Z^{b'}$ 的对角线为零，而 $\mathrm{Tr}\,Z^{b'-b}=\sum_j\omega^{j(b'-b)}=k\delta_{bb'}$。所以
   $B_{ab}=X^aZ^b/\sqrt k$ 构成 $M_k$ 关于 $\mathrm{Tr}(A^*B)$ 的一组规范正交基，而任意
   两组规范正交基之间由一个酉矩阵联系；由定理 6.7，$\sum_mB_m\rho B_m^*$ 对所有规范正交基
   都相同。对基 $E_{ab}$，它等于
   $\sum_{ab}E_{ab}\rho E_{ba}=\mathrm{Tr}(\rho)1$。再除以 $k$ 即可。

4. 对所有 $x$ 和 $z$，$\mathrm{tr}(\Psi(\Phi(x))z)=\mathrm{tr}(\Phi(x)\Psi^\sharp(z))=\mathrm{tr}(x\,\Phi^\sharp(\Psi^\sharp(z)))$；
   由对偶的唯一性即得结论。

</details>

**下一章：** [第 7 章](07_tensor_products.md)讨论张量积与纠缠：偏迹——
一个不同代数之间的信道，其对偶是嵌入 $A\mapsto A\otimes1$——Page
公式，以及荷的融合。
