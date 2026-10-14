[English](../12_equivariant_learning.md) · [Español](../es/12_equivariant_learning.md) · [Français](../fr/12_equivariant_learning.md) · [Deutsch](../de/12_equivariant_learning.md) · **中文** · [日本語](../ja/12_equivariant_learning.md) · [Русский](../ru/12_equivariant_learning.md)

# 12. 等变机器学习

神经网络的一层称为**等变的**，如果先用一个对称变换作用于输入再应用该层，与先应用该层再变换输出
得到相同的结果：对图的顶点重新编号，会以同样的方式对其特征重新编号。对线性层而言这是一个交换
关系，因此等变层构成一个交换子，由第 2 章，在相差一个基变换的意义下它就是代数 $\bigoplus_cM_{m_c}(\mathbb C)$。
本章计算它们的参数个数，用 `torch_vn_algebra.decompose` 构造并训练它们，并利用代数结构得到按构造
即为等变的酉层、正层与谱层。

**你将学到**为什么两个表示之间的等变线性映射是分块矩阵，其块的大小由重数给出；如何通过轨道或
特征标来计数它们（DeepSets 以及图网络的 Bell 数）；为什么在群上取平均就是到等变映射上的正交投影；
函数演算与酉参数化如何在交换子上逐块作用；以及如何计数不变向量、不变二次型以及各次数的不变多项式
（Molien 公式）。

这些层在[教程 12](../../tutorial/zh/12_equivariant_learning.md)中用本库构造和训练，该教程也数值地检验本章的结果。

在本章中，$G$ 是有限群，$G$ 在 $\mathbb C^n$ 上的**表示**是取酉值的同态
$\pi:G\to M_n(\mathbb C)$（这不是限制：可逆的 $\pi(g)$ 对内积
$\sum_g\langle\pi(g)x,\pi(g)y\rangle$ 是酉的，见 4.1 节）。它的**特征标**是 $\chi_\pi(g)=\mathrm{Tr}\,\pi(g)$。不可约
表示、等价表示，以及不可约表示 $\rho$ 在 $\pi$ 中的**重数** $m_\rho(\pi)$（若 $\rho$ 不出现则为 $0$）
均如 2.4 节所述，该节证明了 $m_\rho(\pi)$ 不依赖于分解，并证明了 Schur 引理；[第 4 章](04_symmetries.md)
对二者有详细讨论（4.1–4.2 节）。在长方矩阵上，$\langle S,T\rangle=\mathrm{Tr}(S^*T)$ 是
**Hilbert–Schmidt** 内积。

<a id="s-12-1"></a>

## 12.1 等变线性映射

**定义。** 设 $\pi_1,\pi_2$ 是 $G$ 在 $\mathbb C^{n_1}$、$\mathbb C^{n_2}$ 上的表示。若对所有
$g\in G$ 有 $T\pi_1(g)=\pi_2(g)T$，则称矩阵 $T\in M_{n_2\times n_1}(\mathbb C)$ 是**等变的**（或称
**缠结算子**）。它们构成向量空间 $\mathrm{Hom}_G(\pi_1,\pi_2)$；当 $\pi_1=\pi_2=\pi$ 时，它就是交换子
$\mathrm{Hom}_G(\pi,\pi)=\pi(G)'$。

机器学习中最常见的表示是置换坐标的表示。设 $G$ 作用在有限集 $\Omega$ 上；$\mathbb C^\Omega$ 上的
**置换表示**为 $\pi(g)e_x=e_{gx}$。例子：$S_n$ 置换一个集合的 $n$ 个点（$\Omega=[n]=\{1,\dots,n\}$），
以及 $S_n$ 对图的顶点重新编号，它作用在 $\Omega=[n]^2$ 上，即通过 $X\mapsto P_\sigma XP_\sigma^T$ 作用在
邻接矩阵上。

<a id="r-12-1"></a>

**引理 12.1（轨道）。** 设 $\pi_i$ 是 $G$ 在 $\mathbb C^{\Omega_i}$ 上的置换表示。则
$T\in\mathrm{Hom}_G(\pi_1,\pi_2)$ 当且仅当对所有 $g$、$y\in\Omega_2$、$x\in\Omega_1$ 有 $T_{gy,gx}=T_{yx}$。
$G$ 在 $\Omega_2\times\Omega_1$ 上各轨道的示性矩阵构成一组基，并且
$\dim\mathrm{Hom}_G(\pi_1,\pi_2)$ 等于这些轨道的个数。

*证明。* $\pi_2(g)T\pi_1(g)^*e_x=\pi_2(g)Te_{g^{-1}x}=\sum_yT_{y,g^{-1}x}e_{gy}$，所以
$\pi_2(g)T\pi_1(g)^*$ 的 $(gy,x)$ 元为 $T_{y,g^{-1}x}$。因此 $T$ 等变，即对所有 $g$ 有
$\pi_2(g)T\pi_1(g)^*=T$，当且仅当对所有 $g,y,x$ 有 $T_{gy,x}=T_{y,g^{-1}x}$，而这（把 $x$ 换成 $gx$）
正是所述条件：$T$ 在轨道上为常数。各轨道的示性矩阵支撑互不相交，因而线性无关，并且每个在轨道上
为常数的函数都是它们的线性组合。$\square$

*计算：* 对作用在顶点特征与边特征上的 $S_4$，通过计数轨道验证引理 12.1——[教程，12.1 节](../../tutorial/zh/12_equivariant_learning.md#s-12-1)。

<a id="r-12-2"></a>

**推论 12.2.** (a)（*DeepSets。*）对作用在 $\mathbb C^n$ 上的 $S_n$，$n\ge2$，交换子为
$\mathrm{span}\{1,J\}$，其中 $J$ 是全 1 矩阵，并且 $\mathbb C^n=\mathbb C\mathbf 1\oplus\mathbf 1^\perp$
是两个不等价不可约表示之和：平凡表示，以及维数为 $n-1$ 的**标准**表示。(b)（*图网络。*）从
$\mathbb C^{[n]^p}$ 到 $\mathbb C^{[n]^q}$（$n$ 个顶点上的 $p$ 阶与 $q$ 阶张量）的等变映射空间的维数，
等于把一个 $p+q$ 元集合划分为至多 $n$ 块的划分个数；当 $n\ge p+q$ 时，这个数是 **Bell 数** $B_{p+q}$，
即一个 $p+q$ 元集合的全部划分的个数（$B_2=2$，$B_3=5$，$B_4=15$）。

*证明。* (b) 两个元组 $x,y\in[n]^{p+q}$ 属于同一轨道，当且仅当它们有相同的**相等模式**，即位置集
$\{1,\dots,p+q\}$ 按 $i\sim j\Leftrightarrow x_i=x_j$ 所得的划分。模式显然是不变的。反之，若 $x,y$
有相同的模式，则 $x_i\mapsto y_i$ 是从 $\{x_i\}$ 到 $[n]$ 的一个定义良好的单射，它可以延拓为满足
$\sigma x=y$ 的置换 $\sigma$。每个至多有 $n$ 块的划分都是某个模式（给各块赋以不同的值）。由引理 12.1
即得结论。当 $n\ge p+q$ 时，对块数的限制不起作用。按块数计数：$B_2=1+1$；$B_3=1+3+1$（$3$ 元集合的两块
划分由其单元素块决定）；$B_4=1+7+6+1$（两块：大小为 $3+1$ 的有 $4$ 个，大小为 $2+2$ 的有 $3$ 个；三块：
那个二元块有 $\binom42=6$ 种选法）。(a) 对 $p=q=1$，$[n]^2$ 上有两个轨道：对角线及其补集，示性矩阵为 $1$ 与 $J-1$。所以
$\dim\pi(S_n)'=2=\sum_cm_c^2$（定理 2.6，命题 2.7），这迫使有两个 $m_c=1$ 的扇区，即两个不等价的
不可约直和项，并且 $\pi(S_n)'$ 由它们的投影 $P_1,P_2$ 张成。到不变直线 $\mathbb C\mathbf 1$ 上的投影
属于 $\pi(S_n)'$（引理 2.1），所以它是 $P_1$ 或 $P_2$：承载平凡表示的这条直线是一个扇区，
$\mathbf 1^\perp$ 是另一个。$\square$

DeepSets 层 $x\mapsto\lambda x+\gamma\,(\mathbf 1^Tx)\mathbf 1$（Zaheer et al. 2017）就是 (a)；
Maron et al.（2019）关于图网络的基就是 (b)。

*计算：* $n=4$ 时的推论 12.2（$B_2=2$，$B_3=5$，$B_4=15$），以及用 `dc.from_group` 对 $S_4$ 的顶点特征与边特征所作的分解（边特征：维数为 $k=(1,2,3,3)$、重数为 $m=(2,1,1,3)$ 的扇区）——[教程，12.1 节](../../tutorial/zh/12_equivariant_learning.md#s-12-1)。

对于两个表示，一切都归结为一个交换子：

<a id="r-12-3"></a>

**引理 12.3（直和）。** 设 $\pi=\pi_1\oplus\pi_2$ 作用在 $\mathbb C^{n_1}\oplus\mathbb C^{n_2}$ 上，等距
嵌入为 $Q_i$，投影为 $P_i=Q_iQ_i^*$。则 $P_1,P_2\in\pi(G)'$，
$\mathrm{Hom}_G(\pi_1,\pi_2)=Q_2^*\,\pi(G)'\,Q_1$，并且 $T\mapsto Q_2TQ_1^*$ 是从
$\mathrm{Hom}_G(\pi_1,\pi_2)$ 到 $P_2\,\pi(G)'P_1$ 上的线性双射。

*证明。* 由于 $\pi(g)$ 是分块对角的，$\pi(g)Q_i=Q_i\pi_i(g)$ 且 $Q_i^*\pi(g)=\pi_i(g)Q_i^*$；因此
$P_i\pi(g)=Q_i\pi_i(g)Q_i^*=\pi(g)P_i$。若 $X\in\pi(G)'$，则
$Q_2^*XQ_1\pi_1(g)=Q_2^*X\pi(g)Q_1=Q_2^*\pi(g)XQ_1=\pi_2(g)Q_2^*XQ_1$。若 $T$ 等变，则
$\pi(g)Q_2TQ_1^*=Q_2\pi_2(g)TQ_1^*=Q_2T\pi_1(g)Q_1^*=Q_2TQ_1^*\pi(g)$，所以 $Q_2TQ_1^*\in\pi(G)'$，它等于
$P_2(Q_2TQ_1^*)P_1$，且 $Q_2^*(Q_2TQ_1^*)Q_1=T$。反之，每个 $Y=P_2XP_1$（$X\in\pi(G)'$）都属于
$\pi(G)'$，且 $Y=Q_2(Q_2^*YQ_1)Q_1^*$。$\square$

*计算：* 顶点特征与边特征的角块 $Q_2^*\,\pi(G)'\,Q_1$——[教程，12.2 节](../../tutorial/zh/12_equivariant_learning.md#s-12-2)。

<a id="r-12-4"></a>

**定理 12.4（等变映射是分块矩阵）。** 设 $\rho_1,\dots,\rho_C$ 是在 $\pi_1\oplus\pi_2$ 中出现的两两
不等价的不可约表示，$k_c=\dim\rho_c$，$a_c=m_{\rho_c}(\pi_1)$，$b_c=m_{\rho_c}(\pi_2)$。存在酉算子
$V_1:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{a_c}\to\mathbb C^{n_1}$ 与
$V_2:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{b_c}\to\mathbb C^{n_2}$，满足
$V_i^*\pi_i(g)V_i=\bigoplus_c\rho_c(g)\otimes1$，使得
$$
\mathrm{Hom}_G(\pi_1,\pi_2)=\Big\{V_2\Big(\bigoplus_{c=1}^C1_{k_c}\otimes y_c\Big)V_1^*:\ y_c\in M_{b_c\times a_c}(\mathbb C)\Big\},
$$
并且 $(y_c)_c\mapsto T$ 是单射。因此 $\mathrm{Hom}_G(\pi_1,\pi_2)\cong\bigoplus_cM_{b_c\times a_c}(\mathbb C)$，
$$
\dim\mathrm{Hom}_G(\pi_1,\pi_2)=\sum_\rho m_\rho(\pi_1)\,m_\rho(\pi_2),
$$
其中求和遍及等价意义下的所有不可约表示；当 $\pi_1=\pi_2$ 时，这就是
$\pi(G)'\cong\bigoplus_cM_{m_c}$（命题 2.7），并且 $\mathrm{Hom}_G(\pi_1,\pi_2)=0$ 当且仅当没有
不可约表示同时出现在两者之中。

*证明。* 对 $\pi=\pi_1\oplus\pi_2$ 应用定理 2.6 与命题 2.7：存在酉算子 $V$，使得
$V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$ 且 $V^*\pi(G)'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$。由重数的
唯一性（2.4 节），在 $\pi$ 中出现的每个不可约表示恰与一个 $\rho_c$ 等价，所以这些 $\rho_c$ 是定理陈述中
各等价类的代表。任何其他代表都形如 $\rho_c'=w_c\rho_cw_c^*$，其中 $w_c$ 为酉矩阵（引理 4.2(c)），而把 $V$
换成 $V\bigoplus_c(w_c^*\otimes1)$ 会把第一个标准形变为 $\bigoplus_c\rho'_c(g)\otimes1_{m_c}$，并且不改变第二个，
因为 $w_c^*\otimes1$ 与 $1\otimes M_{m_c}$ 可交换；所以我们可以使用命题 2.7 中的 $\rho_c$。由于
$X\mapsto V^*XV$ 是 $*$-同构且 $P_1\in\pi(G)'$（引理 12.3），有 $V^*P_1V=\bigoplus_c1\otimes p_c$，其中
$p_c$ 为投影。选取酉矩阵 $u_c$ 使 $u_c^*p_cu_c=1_{a_c}\oplus0_{b_c}$（定理 1.3），这里暂时令
$a_c=\mathrm{rank}\,p_c$，$b_c=m_c-a_c$，并把 $V$ 换成 $V\bigoplus_c(1\otimes u_c)$；这不改变两个
标准形，因为 $1\otimes u_c$ 与 $\rho_c(g)\otimes1$ 可交换。设 $\iota_c:\mathbb C^{a_c}\to\mathbb C^{m_c}$ 与
$\kappa_c:\mathbb C^{b_c}\to\mathbb C^{m_c}$ 分别是前 $a_c$ 个坐标与后 $b_c$ 个坐标的嵌入，并令
$$
V_1=Q_1^*V\bigoplus_c(1\otimes\iota_c),\qquad V_2=Q_2^*V\bigoplus_c(1\otimes\kappa_c).
$$
$W_1=V\bigoplus_c(1\otimes\iota_c)$ 是值域为 $\mathrm{range}\,P_1=\mathrm{range}\,Q_1$ 的等距算子，所以
$V_1^*V_1=W_1^*P_1W_1=1$，又因为 $\sum_ck_ca_c=\mathrm{rank}\,P_1=n_1$，$V_1$ 是酉的。此外
$\pi_1(g)V_1=Q_1^*\pi(g)W_1=Q_1^*V\big(\bigoplus_c\rho_c(g)\otimes1_{m_c}\big)\bigoplus_c(1\otimes\iota_c)
=V_1\bigoplus_c\rho_c(g)\otimes1_{a_c}$。所以 $\pi_1\cong\bigoplus_c\rho_c\otimes1_{a_c}$，其中各 $\rho_c$
为不等价的不可约表示，由重数的唯一性（2.4 节）得 $a_c=m_{\rho_c}(\pi_1)$；对 $V_2$ 与 $b_c$ 同理。
最后，对 $X=V\big(\bigoplus_c1\otimes x_c\big)V^*\in\pi(G)'$，有
$P_2XP_1=V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*$，其中 $y_c=\kappa_c^*x_c\iota_c$ 是 $x_c$ 左下角的
$b_c\times a_c$ 子块，它是任意矩阵。由引理 12.3，
$\mathrm{Hom}_G(\pi_1,\pi_2)=Q_2^*P_2\,\pi(G)'P_1Q_1$ 由矩阵
$Q_2^*V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*Q_1=V_2\big(\bigoplus_c1\otimes y_c\big)V_1^*$ 构成，并且由于
$V_1,V_2$ 是酉的，$(y_c)\mapsto T$ 是单射。不在 $\pi_1\oplus\pi_2$ 中出现的不可约表示 $\rho$ 满足
$m_\rho(\pi_1)m_\rho(\pi_2)=0$。$\square$

*库中：* 作为算法的定理 12.4 的证明，建立在 `dc.from_group`、`D.commutant_algebra()` 与 `D.from_commutant` 之上的层 `EquivariantLinear` 与 `EquivariantMap`；对 $S_4$ 从顶点特征到边特征，$\dim\mathrm{Hom}_G=1\cdot2+1\cdot3=5=B_3$——[教程，12.2 节](../../tutorial/zh/12_equivariant_learning.md#s-12-2)。

> **物理意义**（物理解释，而非本课程的定理）。定理 12.4 是**选择定则**的线性代数核心，也是 Wigner–Eckart 定理的标量情形：与
> 对称性可交换的算子只能连接对称类型相同的态，并且在同一类型内部它以 $1_{k_c}\otimes y_c$ 的方式作用，
> 即通过一个对多重态内部的位置不敏感的“约化矩阵元” $y_c$ 作用。（本身按非平凡表示变换的算子，例如
> 矢量算子，还需要 Clebsch–Gordan 系数。）对称的哈密顿量是 $\pi_1=\pi_2$ 的情形；它的简并度 $k_c$
> 就是各不可约表示的维数。

<a id="s-12-2"></a>

## 12.2 用库构造等变层

*库中：* 以定理 12.4 的块为参数的等变层——[教程，12.2 节](../../tutorial/zh/12_equivariant_learning.md#s-12-2)。

使层等变的第二种方法是在群上对它取平均。**twirl** 定义为
$$
\mathcal P(W)=\frac1{|G|}\sum_{g\in G}\pi_2(g)\,W\,\pi_1(g)^*,\qquad W\in M_{n_2\times n_1}(\mathbb C).
$$

<a id="r-12-5"></a>

**定理 12.5（twirl 是正交投影）。** $\mathcal P$ 是 $M_{n_2\times n_1}(\mathbb C)$ 到
$\mathrm{Hom}_G(\pi_1,\pi_2)$ 上关于 Hilbert–Schmidt 内积的正交投影，并且
$$
\dim\mathrm{Hom}_G(\pi_1,\pi_2)=\frac1{|G|}\sum_{g\in G}\overline{\chi_{\pi_1}(g)}\,\chi_{\pi_2}(g).
$$

*证明。* 对 $h\in G$，$\pi_2(h)\mathcal P(W)\pi_1(h)^*=\frac1{|G|}\sum_g\pi_2(hg)W\pi_1(hg)^*=\mathcal P(W)$，因为
$g\mapsto hg$ 是 $G$ 的双射；所以 $\mathcal P(W)$ 是等变的。若 $W$ 等变，则每一项都等于 $W$，所以
$\mathcal P(W)=W$。因此 $\mathcal P^2=\mathcal P$，其值域为 $\mathrm{Hom}_G(\pi_1,\pi_2)$。它是自伴的：
$\langle S,\mathcal P(W)\rangle=\frac1{|G|}\sum_g\mathrm{Tr}\big((\pi_2(g)^*S\pi_1(g))^*W\big)$
$=\langle\mathcal P(S),W\rangle$，这里作了代换 $g\mapsto g^{-1}$ 并用到 $\pi_i(g^{-1})=\pi_i(g)^*$。自伴
幂等算子就是到其值域上的正交投影。它的迹等于它的秩。$M_{n_2\times n_1}$ 上的线性映射
$W\mapsto AWB$ 的迹为 $\mathrm{Tr}\,A\mathrm{Tr}\,B$（矩阵单位 $E_{ij}$ 在 $AE_{ij}B$ 中的系数为
$A_{ii}B_{jj}$）；取 $A=\pi_2(g)$，$B=\pi_1(g)^*$ 即得该公式。$\square$

结合定理 12.4，$\sum_\rho m_\rho(\pi_1)m_\rho(\pi_2)=|G|^{-1}\sum_g\overline{\chi_{\pi_1}(g)}\chi_{\pi_2}(g)$；
当 $\pi_1=\pi_2$ 时，这就是 4.2 节中 $\dim\pi(G)'$ 的特征标公式。
参数化 $W\mapsto\mathcal P(W)$ 能取到每一个等变映射，$\mathcal P(W)$ 是离 $W$ 最近的等变映射，并且由于
$\mathcal P$ 是自伴的，关于 $W$ 的梯度就是关于 $T=\mathcal P(W)$ 的梯度在 $\mathcal P$ 下的像。但它的参数中有
$n_1n_2-\dim\mathrm{Hom}_G$ 个浪费在核上。

*计算：* 在边特征上验证定理 12.5，以及从顶点特征到边特征的维数公式——[教程，12.2 节](../../tutorial/zh/12_equivariant_learning.md#s-12-2)。*库中：* twirl `dc.twirl`，由分解计算的同一个投影 `D.conditional_expectation(W, 'commutant')`，以及特征标公式 `dc.commutant_dimension`——[教程，12.2 节](../../tutorial/zh/12_equivariant_learning.md#s-12-2)。

上面的层都是线性的。网络让它们与非线性交替出现，而这些非线性也必须是等变的。对置换表示，逐坐标
作用的函数是等变的，因为 $(\pi(g)x)_{gy}=x_y$；对其他表示则不然，此时人们使用作用于范数的非线性，
或者作用于一组被 $G$ 置换的基的系数上的非线性（习题 3）。

> **物理意义**（物理解释，而非本课程的定理）。twirl 是到交换子上的条件期望（第 2 章习题 3）：把一个耦合在所有对称操作上取平均，
> 恰好保留它的对称部分。通过对对称群求和来“对称化”哈密顿量或相互作用，就是同一个投影。

<a id="s-12-3"></a>

## 12.3 带约束的层与谱层

等变层 $x\mapsto Tx$（$T\in\pi(G)'$）的权重 $w=(w_c)_c$ 位于 $*$-代数 $\bigoplus_cM_{m_c}$ 中（定理 12.4），
所以第 1 章与第 11 章的每个构造都适用于它们。对 $w\in\bigoplus_cM_{m_c}$ 记
$\Phi(w)=V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*$。

<a id="r-12-6"></a>

**命题 12.6.** (a) $\Phi$ 是从 $\bigoplus_cM_{m_c}$ 到 $\pi(G)'$ 上的 $*$-同构，并且 $\Phi(w)$ 是自伴的、
酉的或正的，当且仅当每个 $w_c$ 都是。(b) 对正规的 $w$ 以及
$\bigcup_c\mathrm{spec}\,w_c$ 上的函数 $f$，有 $f(\Phi(w))=\Phi\big((f(w_c))_c\big)$；特别地，对每个正规的
$X\in\pi(G)'$ 有 $f(X)\in\pi(G)'$。(c) 对任意固定的酉元 $u_0\in\bigoplus_cM_{m_c}$，等变酉算子恰为
$\Phi(u_0e^{s-s^*})$，其中 $s$ 取遍 $\bigoplus_cM_{m_c}$；特别地，它们恰为反自伴的 $S\in\pi(G)'$ 的指数 $e^S$。

*证明。* (a) 由于 $(1\otimes x)(1\otimes y)=1\otimes xy$ 且 $(1\otimes x)^*=1\otimes x^*$，$\Phi$ 是线性、
保单位、可乘且保持 $*$ 的；它是单射，并由定理 2.6 映满 $\pi(G)'$。因此 $\Phi(w)^*=\Phi(w)$ 当且仅当
$w^*=w$，且 $\Phi(w)^*\Phi(w)=1=\Phi(w)\Phi(w)^*$ 当且仅当 $w^*w=1=ww^*$。若所有 $w_c=t_c^*t_c$，则
$\Phi(w)=\Phi(t)^*\Phi(t)\ge0$（推论 1.7）。反之，若 $\Phi(w)\ge0$，则对 $v\in\mathbb C^{m_c}$ 以及扇区 $c$
中的向量 $\xi=V(e_1\otimes v)$，有 $0\le\langle\xi,\Phi(w)\xi\rangle=\langle v,w_cv\rangle$。
(b) 设 $w_c=u_c\Lambda_cu_c^*$（定理 1.3）。则 $\Phi(w)=U\big(\bigoplus_c1\otimes\Lambda_c\big)U^*$，其中酉算子
$U=V\bigoplus_c(1\otimes u_c)$；这是一个对角化，$w_c$ 的特征值各重复 $k_c$ 次；由定义与定理 1.5(a)，
$f(\Phi(w))=U\big(\bigoplus_c1\otimes f(\Lambda_c)\big)U^*=\Phi\big((f(w_c))_c\big)$。
(c) 由 (a)，等变酉算子是 $\Phi(u)$，其中 $u$ 是 $\bigoplus_cM_{m_c}$ 中的酉元，而由引理 11.1，这些
$u$ 恰为 $u_0e^{s-s^*}$。取 $u_0=1$，(b) 给出 $\Phi(e^{s-s^*})=e^{\Phi(s-s^*)}$，而 $\Phi(s-s^*)$ 取遍
$\pi(G)'$ 的反自伴元素。$\square$

*库中：* 用 `D.commutant_algebra()` 上的 `optimize.UnitaryParam` 与 `PositiveParam` 构造酉的和正的等变层，用 `apply_function` 构造谱层；训练一个等变酉层，以及对 $e^{-tL}$ 及其关于 $t$ 的导数验证 (b)——[教程，12.3 节](../../tutorial/zh/12_equivariant_learning.md#s-12-3)。

> **物理意义**（物理解释，而非本课程的定理）。对称的哈密顿量 $H\in\pi(G)'$ 生成对称的时间演化
> $e^{-itH}\in\pi(G)'$（第 8 章），它从不混合扇区，只作用在重数空间上。等变酉层就是这样一种演化，
> 只不过它是学习得到的，而不是给定的。

<a id="s-12-4"></a>

## 12.4 不变量

若对所有 $g$ 有 $\pi(g)x=x$，则称向量 $x$ 是**不变的**；若 $q(\pi(g)x)=q(x)$，则称 $\mathbb C^n$ 上的
函数 $q$ 是不变的。不变特征就是不变网络（例如图分类器）所能输出的东西。

<a id="r-12-7"></a>

**命题 12.7.** (a) 不变向量构成正交投影 $P_0=|G|^{-1}\sum_g\pi(g)$ 的值域；其维数为
$m_{\mathrm{triv}}(\pi)=|G|^{-1}\sum_g\chi_\pi(g)$，并且若平凡表示作为 $\rho_{c_0}$ 出现，则 $P_0$ 是
$\pi(G)''$ 的极小中心投影 $P_{c_0}$。
(b) Hermite 型 $q_A(x)=\langle x,Ax\rangle$（$A=A^*$）是不变的，当且仅当 $A\in\pi(G)'$；不变 Hermite 型
构成维数为 $\sum_cm_c^2$ 的实向量空间。
(c) 关于坐标 $x_1,\dots,x_n$ 的 $2$ 次齐次不变多项式构成维数为
$$
\frac1{|G|}\sum_{g\in G}\frac{\chi_\pi(g)^2+\chi_\pi(g^2)}2 .
$$
的空间。对置换表示，它等于 $G$ 在由 $\Omega$ 的两个元素组成的多重集 $\{x,y\}$ 上的轨道个数。

*证明。* (a) 不变向量就是 $\mathrm{Hom}_G(\mathrm{triv},\pi)$，其中 $\mathrm{triv}$ 是
$\mathbb C=M_{1\times1}$ 上的平凡表示；对它而言 twirl 为 $\mathcal P(x)=P_0x$，所以定理 12.5 给出投影与
维数，定理 12.4 给出 $m_{\mathrm{triv}}(\pi)$。在定理 2.6 的基下，把不变向量写成 $\xi=\sum_c\xi_c$，其中
$\xi_c=\sum_jv_{c,j}\otimes e_j\in\mathbb C^{k_c}\otimes\mathbb C^{m_c}$。由于 $\pi(g)$ 保持各扇区，对所有 $g$ 有
$\rho_c(g)v_{c,j}=v_{c,j}$。非零的 $v_{c,j}$ 张成一条承载平凡表示的不变直线，于是 $\rho_c$ 作为不可约
表示是平凡的，且 $c=c_0$。因此不变向量就是扇区 $c_0$ 中的向量，$\pi(g)$ 在那里作用为 $1$：它们构成
$P_{c_0}$ 的值域，因此 $P_{c_0}$ 就是 $P_0$。
(b) $q_A(\pi(g)x)=\langle x,\pi(g)^*A\pi(g)x\rangle$，而对所有 $x$ 满足 $\langle x,Dx\rangle=0$ 的矩阵 $D$ 必为零
（极化，见推论 1.4 的证明）。所以 $q_A$ 不变当且仅当对所有 $g$ 有 $\pi(g)^*A\pi(g)=A$，即 $A\in\pi(G)'$。
由命题 12.6(a)，$\pi(G)'$ 的自伴元素是各 $w_c$ 为 Hermite 矩阵的 $\Phi(w)$，它们构成维数为
$\sum_cm_c^2$ 的实空间。
(c) 二次多项式 $p(x)=x^TBx$ 由其对称矩阵 $B$ 决定。令 $(g\cdot p)(x)=p(\pi(g)^{-1}x)$；则
$g\cdot(h\cdot p)=(gh)\cdot p$，且 $g\cdot p$ 的矩阵为 $\sigma(g)(B)=A^TBA$，$A=\pi(g)^{-1}$。平均
$M=|G|^{-1}\sum_g\sigma(g)$ 满足 $\sigma(h)M=M$，并且对不变的 $B$ 有 $MB=B$，所以（与定理 12.5 的证明
一样）它是幂等的，其值域为不变多项式，而不变多项式空间的维数为
$\mathrm{Tr}\,M=|G|^{-1}\sum_g\mathrm{Tr}\,\sigma(g)$（幂等映射在与值域 $\oplus$ 核相适应的基下是
$1\oplus0$，所以它的迹等于它的秩）。$A$ 是酉的，因此 $A=R\Lambda R^{-1}$，其中
$\Lambda=\mathrm{diag}(\lambda)$（定理 1.3）。对称矩阵
$B_{ij}=(R^{-1})^T(e_ie_j^T+e_je_i^T)R^{-1}$（$i\le j$）构成一组基，并且
$A^TB_{ij}A=(R^{-1})^T\Lambda(e_ie_j^T+e_je_i^T)\Lambda R^{-1}=\lambda_i\lambda_jB_{ij}$。所以
$\mathrm{Tr}\,\sigma(g)=\sum_{i\le j}\lambda_i\lambda_j=\frac12\big((\mathrm{Tr}\,A)^2+\mathrm{Tr}\,A^2\big)
=\frac12\big(\chi_\pi(g^{-1})^2+\chi_\pi(g^{-2})\big)$，再把对 $g$ 的求和换成对 $g^{-1}$ 的求和即得该公式。
对置换表示，$G$ 置换单项式 $x_ax_b$，所以多项式不变当且仅当其系数在单项式的轨道上为常数。$\square$

*计算：* 对 $S_4$ 的边特征验证 (a)、(b) 以及 (c) 的两种计数——[教程，12.4 节](../../tutorial/zh/12_equivariant_learning.md#s-12-4)。

对于图，不变边特征由对角线及其补集的示性向量张成（$B_2=2$），不变的线性读出是对角元之和与非对角元
之和。在顶点特征 $\mathbb C^4$ 上，(c) 计数的是 $S_4$ 在由两个顶点组成的多重集上的轨道：共有两个。

在 $\mathbb C^4$ 上，它们是 $\sum_ix_i^2$ 与 $\sum_{i<j}x_ix_j$。(c) 的证明对每个次数都适用：

<a id="r-12-8"></a>

**命题 12.8（Molien 公式）。** 设 $\mathbb C[x]_d$ 是关于 $x_1,\dots,x_n$ 的 $d$ 次齐次多项式构成的空间，
$\mathbb C[x]^G_d$ 是其中不变多项式（对所有 $g$ 有 $p(\pi(g)x)=p(x)$）构成的子空间。对 $|t|<1$，
$$
\sum_{d\ge0}\dim\mathbb C[x]^G_d\,t^d=\frac1{|G|}\sum_{g\in G}\frac1{\det\big(1-t\pi(g)\big)} .
$$
对置换表示，$\dim\mathbb C[x]^G_d$ 等于 $G$ 在 $d$ 次单项式（即由 $\Omega$ 的 $d$ 个元素组成的多重集）上的
轨道个数。$t^2$ 的系数就是命题 12.7(c)。

*证明。* 对可逆矩阵 $M$，令 $\phi_M(p)(x)=p(Mx)$（代入）。$\phi_M$ 是线性的，把 $\mathbb C[x]_d$ 映到自身，
并且 $\phi_N\circ\phi_M=\phi_{MN}$，所以 $\phi_M$ 是双射，其逆为 $\phi_{M^{-1}}$。因此
$\sigma(g)=\phi_{\pi(g)^{-1}}$ 满足 $\sigma(g)\sigma(h)=\sigma(gh)$，并且 $p$ 不变当且仅当对所有 $g$ 有
$\sigma(g)p=p$。与命题 12.7(c) 的证明一样，$\mathbb C[x]_d$ 上的 $M_d=|G|^{-1}\sum_g\sigma(g)$ 是幂等的，
其值域为 $\mathbb C[x]^G_d$，所以 $\dim\mathbb C[x]^G_d=\mathrm{Tr}\,M_d=|G|^{-1}\sum_g\mathrm{Tr}\,\sigma_d(g)$，
其中 $\sigma_d(g)$ 是 $\sigma(g)$ 在 $\mathbb C[x]_d$ 上的限制。固定 $g$，把酉矩阵 $A=\pi(g)^{-1}$ 写成
$A=R\Lambda R^{-1}$，$\Lambda=\mathrm{diag}(\lambda_1,\dots,\lambda_n)$，$|\lambda_i|=1$（定理 1.3）。借助线性型
$\ell_i(x)=(R^{-1}x)_i$，乘积 $\ell^\alpha=\prod_i\ell_i^{\alpha_i}=\phi_{R^{-1}}(x^\alpha)$（$|\alpha|=d$）构成
$\mathbb C[x]_d$ 的一组基（单项式基在双射 $\phi_{R^{-1}}$ 下的像），并且由于 $R^{-1}A=\Lambda R^{-1}$，
$$
(\sigma(g)\ell^\alpha)(x)=\prod_i\big((R^{-1}Ax)_i\big)^{\alpha_i}=\prod_i\big(\lambda_i\ell_i(x)\big)^{\alpha_i}=\lambda^\alpha\ell^\alpha(x).
$$
所以 $\mathrm{Tr}\,\sigma_d(g)=\sum_{|\alpha|=d}\lambda^\alpha$。对 $|t|<1$，$n$ 个绝对收敛级数
$(1-t\lambda_i)^{-1}=\sum_{a\ge0}t^a\lambda_i^a$ 的乘积可以展开并按次数重新分组：
$$
\frac1{\det(1-tA)}=\prod_{i=1}^n\frac1{1-t\lambda_i}=\sum_\alpha t^{|\alpha|}\lambda^\alpha=\sum_{d\ge0}t^d\mathrm{Tr}\,\sigma_d(g).
$$
对 $g$ 取平均并把 $g$ 换成 $g^{-1}$（使得 $A=\pi(g)$）即得该公式；两边都收敛，因为
$\dim\mathbb C[x]^G_d\le\dim\mathbb C[x]_d=\binom{n+d-1}d$。对置换表示，$(\pi(g)^{-1}x)_y=x_{gy}$，所以 $\sigma(g)$
把单项式 $\prod_jx_{y_j}$ 映为单项式 $\prod_jx_{gy_j}$：$G$ 置换 $d$ 次单项式，多项式不变当且仅当其系数在
各轨道上为常数，而各轨道上单项式之和构成 $\mathbb C[x]^G_d$ 的一组基。$\square$

（这是 1897 年的 Molien 定理；延伸阅读参见 B. Sturmfels, *Algorithms in Invariant Theory*, Springer 2008，
2.2 节。）对作用在 $\mathbb C^4$ 上的 $S_4$，$d\le4$ 时的计数 $1,1,2,3,5$ 就是把 $d$ 分拆为至多四个部分的分拆个数
（多重集的轨道由其排序后的重数决定）。

*计算：* 对顶点特征与边特征，把展开到 $t^4$ 的 Molien 级数与单项式的轨道作比较——[教程，12.4 节](../../tutorial/zh/12_equivariant_learning.md#s-12-4)。

> **物理意义**（物理解释，而非本课程的定理）。在 Landau 相变理论中，自由能是实序参量 $x$ 的不变多项式，
> $G$ 通过实正交矩阵作用在 $x$ 上。它的二次部分是 $x^TBx$，其中 $B$ 属于交换子（因为
> $\pi(g)^T=\pi(g)^{-1}$），所以由命题 12.6(a)，它在各扇区上是分块对角的。Landau 理论进而假定：对称相在
> $B$ 的某个特征值首先改变符号的那个扇区中变得不稳定，并且在一般情形下，该扇区的不可约表示决定了有序相
> 的对称性。

<a id="s-12-5"></a>

## 12.5 可与之比较的工具

*库中：* 与其他软件（escnn、e3nn、EMLP）的比较，以及本库补充了什么——[教程，12.5 节](../../tutorial/zh/12_equivariant_learning.md#s-12-5)。

<a id="x-12"></a>

## 习题

1. （**带通道的 DeepSets。**）$S_n$（$n\ge2$）通过 $P_\sigma\otimes1_d$ 作用在 $\mathbb C^n\otimes\mathbb C^d$ 上。
   证明从 $\mathbb C^n\otimes\mathbb C^{d_1}$ 到 $\mathbb C^n\otimes\mathbb C^{d_2}$ 的等变映射恰为
   $1\otimes A+J\otimes B$，其中 $A,B\in M_{d_2\times d_1}$。
2. （**实的层。**）设 $\pi_1,\pi_2$ 的矩阵是实的（例如置换表示）。证明实的等变矩阵构成维数为
   $\dim_{\mathbb C}\mathrm{Hom}_G(\pi_1,\pi_2)$ 的实向量空间，并且 $T\mapsto\mathrm{Re}\,T$ 把
   $\mathrm{Hom}_G(\pi_1,\pi_2)$ 映满这个空间。
3. （**非线性。**）(a) 证明对每个坐标应用函数 $f:\mathbb C\to\mathbb C$ 与每个置换表示都可交换。
   (b) 证明对 $f(t)=\max(t,0)$，它与 $\mathbb R^2$ 上转角为 $\pi/2$ 的旋转不可交换。(c) 证明
   $x\mapsto\varphi(\|x\|)x$ 与每个酉表示都可交换。
4. （**多少个样本能确定一个层？**）设 $T=\Phi(w)\in\pi(G)'$，$x_1,\dots,x_N\in\mathbb C^n$。把 $V^*x_i$
   在扇区 $c$ 中的部分写成 $\sum_{a,j}(X_{i,c})_{aj}e_a\otimes e_j$，其中 $X_{i,c}\in M_{k_c\times m_c}$，并设
   $S_c\in M_{Nk_c\times m_c}$ 是把矩阵 $X_{1,c},\dots,X_{N,c}$ 上下堆叠所得的矩阵。证明满足
   $Tx_1=\dots=Tx_N=0$ 的等变 $T$ 构成维数为 $\sum_cm_c(m_c-\mathrm{rank}\,S_c)$ 的空间；对相互独立的
   标准复高斯向量 $x_i$，它几乎必然等于 $\sum_cm_c\max(0,m_c-Nk_c)$。由此推出：等变的 $T$ 由
   $Tx_1,\dots,Tx_N$（在等变映射之中）确定，当且仅当这个维数为 $0$；并且对边表示，两个高斯样本几乎必然
   确定一个等变层，一个样本则永远不能；而任意的 $16\times16$ 矩阵永远不能由少于十六个样本确定。

*计算：* 习题 1、3 和 4 的数值检验——[教程，习题](../../tutorial/zh/12_equivariant_learning.md#x-12)。

<details>
<summary>解答</summary>

1. $P_\sigma J=J=JP_\sigma$，所以这些映射是等变的。它们线性无关：若
   $1\otimes A+J\otimes B=0$，则非对角块给出 $B=0$，进而 $A=0$。由推论 12.2(a) 与
   定理 2.6，$P_\sigma\otimes1_d\cong(\mathrm{triv}\otimes1_d)\oplus(\mathrm{std}\otimes1_d)$，所以两个不可约
   表示的重数都是 $d$，而定理 12.4 给出 $\dim\mathrm{Hom}_G=2d_1d_2$，这正是由形如
   $1\otimes A+J\otimes B$ 的映射构成的空间的维数。

2. 若 $T$ 等变，则 $\bar T$ 也等变，因为 $\overline{\pi_i(g)}=\pi_i(g)$。因此
   $\mathrm{Re}\,T=(T+\bar T)/2$ 与 $\mathrm{Im}\,T=(T-\bar T)/2i$ 是实的且等变，所以映射
   $T\mapsto\mathrm{Re}\,T$ 映满实的等变矩阵（在它们上面它是恒等映射），并且
   $\mathrm{Hom}_G(\pi_1,\pi_2)=H_{\mathbb R}\oplus iH_{\mathbb R}$，其中 $H_{\mathbb R}$ 是实的等变
   矩阵全体。于是 $H_{\mathbb R}$ 的一组实基就是 $\mathrm{Hom}_G$ 的一组复基：它张成整个空间，而对
   实的 $T_j$，复线性关系 $\sum_j(\alpha_j+i\beta_j)T_j=0$ 可以拆成实部与虚部。在引理 12.1 中，轨道的
   示性矩阵就是这样一组实基。（即使数据是实的，定理 12.4 的块 $y_c$ 也是复的，因为库计算出的 $V$ 不一定
   是实的；用复的 $y_c$ 参数化 $\mathrm{Re}\,T$ 所用的实参数是所需的两倍。）

3. (a) $(\pi(g)f(x))_{gy}=f(x)_y=f(x_y)=f((\pi(g)x)_{gy})$。(b) 对 $R(x_1,x_2)=(-x_2,x_1)$ 与 $x=(1,-1)$，
   $f(Rx)=f(1,1)=(1,1)$，但 $Rf(x)=R(1,0)=(0,1)$。(c) $\|\pi(g)x\|=\|x\|$，所以
   $\varphi(\|\pi(g)x\|)\pi(g)x=\pi(g)\varphi(\|x\|)x$。

4. $(1\otimes w_c)(e_a\otimes e_j)=\sum_le_a\otimes(w_c)_{lj}e_l$，所以 $V^*Tx_i$ 在扇区 $c$ 中的系数矩阵为
   $X_{i,c}w_c^T$。因此对所有 $i$ 有 $Tx_i=0$，当且仅当对所有 $c$ 有 $S_cw_c^T=0$，即当且仅当 $w_c^T$ 的
   每一列都属于 $\ker S_c$，这是一个维数为 $m_c-\mathrm{rank}\,S_c$ 的空间：具有这一性质的 $w_c$ 构成
   维数为 $m_c(m_c-\mathrm{rank}\,S_c)$ 的空间，而 $\Phi$ 是单射。对高斯的 $x_i$，$V^*x_i$ 的坐标就是
   $x_i$ 在由 $V$ 的各列构成的标准正交基下的坐标，它们仍是相互独立的标准复高斯变量（命题 A.16(a)），并且
   对不同的 $i$ 相互独立；所以 $S_c$ 是一个 $Nk_c\times m_c$ 的 Ginibre 矩阵，它有密度（命题 A.16(b)）。它的
   秩至多为 $\min(Nk_c,m_c)$，并且除非该尺寸的左上角方子矩阵奇异，否则等号成立；该子矩阵的 $|\det|^2$ 是
   各元素实部与虚部的实多项式，且不恒为零（当该子矩阵为单位矩阵时它等于 $1$），所以由推论 A.7(b)，秩
   几乎必然为 $\min(Nk_c,m_c)$。对边表示，$(k,m)=(1,2),(2,1),(3,1),(3,3)$（计算出的分解），只有扇区
   $(1,2)$ 有贡献：$N=1$ 时为 $2\cdot(2-1)=2$，$N=2$ 时为 $0$。任意的 $16\times16$ 矩阵是平凡群的情形（一个
   扇区，$k=1$，$m=16$）：对任意选取的 $N<16$ 个样本，维数为 $16(16-\mathrm{rank}\,S)\ge16(16-N)>0$，因为
   $\mathrm{rank}\,S\le N$；而对 $N=16$ 个高斯样本，它几乎必然为 $0$。由于 $T\mapsto(Tx_1,\dots,Tx_N)$ 是
   线性的，两个等变映射在样本上取值相同当且仅当它们的差属于刚才计算的空间，这就给出了判别准则。当 $N=1$
   时，对*每个* $x_1$，该空间的维数都至少为 $2$，因为在扇区 $(1,2)$ 中 $\mathrm{rank}\,S_c\le Nk_c=1<m_c=2$。

</details>

**下一步。** 最后一章，[第 13 章](13_jordan_algebras.md)，讨论量子力学为何使用 *-代数而不是若尔当代数。[第 4 章](04_symmetries.md)发展了这里用到的表示论：Schur 引理、
特征标，以及 `dc.from_group` 背后的随机化算法（4.1–4.3 节）。[第 9 章](09_spin_chains.md)把同样的扇区应用于
自旋链，而[第 5 章](05_random_matrices.md)中 $\prod_cU(m_c)$ 上的 Haar 测度给出了等变酉层的自然随机
初始化。数学方面参见 J.-P. Serre, *Linear Representations of Finite Groups* (Springer 1977) 与
W. Fulton and J. Harris, *Representation Theory* (Springer 1991)；等变学习方面参见 M. Bronstein, J. Bruna,
T. Cohen and P. Veličković, *Geometric Deep Learning* (2021) 与 M. Weiler et al., *Equivariant and Coordinate
Independent Convolutional Networks* (2023)；分解算法方面参见 `torch_vn_algebra/decompose.py` 中的参考文献。
