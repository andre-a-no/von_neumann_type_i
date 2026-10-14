[English](../02_von_neumann_algebras.md) · [Español](../es/02_von_neumann_algebras.md) · [Français](../fr/02_von_neumann_algebras.md) · [Deutsch](../de/02_von_neumann_algebras.md) · **中文** · [日本語](../ja/02_von_neumann_algebras.md) · [Русский](../ru/02_von_neumann_algebras.md)

# 2. 有限维冯·诺依曼代数

冯·诺依曼代数是在"取两次交换子"运算下封闭的算子代数。在有限维情形，这是一个纯代数的
条件，而本章的主要定理断言：在相差一个酉基变换的意义下，每个这样的代数都是若干全矩阵代数
的直和，每个全矩阵代数以一定的重数重复出现。这正是本库所依据的结构。

*库中：* 代数 $\bigoplus_cM_{k_c}(\mathbb C)$ 就是类 `TypeIAlgebra`——[教程，2.5 节](../../tutorial/zh/02_von_neumann_algebras.md#s-2-5)。

**你将学到**双交换子定理、带完整证明的结构定理、对称群如何产生这样的代数、什么是中心，
以及什么是条件期望。本章的数值检验见
[教程 2](../../tutorial/zh/02_von_neumann_algebras.md)。

在本章中，*$\mathbb C^n$ 上的 $*$-代数*指的是一个复线性子空间
$\mathcal M\subseteq M_n(\mathbb C)$，它在乘积和伴随运算下封闭，并且包含
单位元 $1$。

<a id="s-2-1"></a>

## 2.1 交换子与双交换子定理

对集合 $S\subseteq M_n(\mathbb C)$，其**交换子**为
$S'=\{X\in M_n(\mathbb C):XA=AX\ \forall A\in S\}$（即与 $S$ 中所有元素可交换的矩阵之集）。
它总是一个含单位元的代数；若 $S$ 在伴随运算下封闭（$S^*=S$），则 $S'$ 是一个 $*$-代数，因为
对所有 $A\in S$ 有 $XA=AX$ 蕴含对所有 $A\in S$ 有 $X^*A^*=A^*X^*$。$S$ 的每个元素都与 $S'$
的每个元素可交换，所以 $S\subseteq S''$；把这一点应用于 $S'$ 得到
$S'\subseteq S'''$，而对 $S\subseteq S''$ 取交换子（取交换子使包含关系反向）
得到 $S'''\subseteq S'$。因此 $S'=S'''$。

*计算：* 交换子作为零空间的数值计算——[教程，2.4 节](../../tutorial/zh/02_von_neumann_algebras.md#s-2-4)。*库中：* `decompose.commutant_basis`、`decompose.commutant_dimension`——[教程，4.2 节](../../tutorial/zh/04_symmetries.md#s-4-2)。

当一个代数作用在 $\mathbb C^n$ 的子空间 $E$ 上时，我们借助一组规范正交基把 $E$ 等同于
$\mathbb C^{\dim E}$，从而在 $\mathbb C^n$ 上陈述的所有结果都适用于它。

**定义。** $\mathbb C^n$ 上的**冯·诺依曼代数**是满足 $\mathcal M=\mathcal M''$ 的
$*$-代数 $\mathcal M$。它的**中心**是
$Z(\mathcal M)=\mathcal M\cap\mathcal M'$。**因子**是满足
$Z(\mathcal M)=\mathbb C1$ 的冯·诺依曼代数。

这一定义是为无穷维希尔伯特空间上的算子而设的，而本课程不讨论无穷维情形。在有限维情形，
条件 $\mathcal M=\mathcal M''$ 是自动成立的：

<a id="r-2-1"></a>

**引理 2.1.** 设 $\mathcal M$ 是 $\mathbb C^n$ 上的 $*$-代数，$K\subseteq\mathbb C^n$
是满足 $\mathcal MK\subseteq K$ 的子空间。则到 $K$ 上的正交投影 $P$ 属于
$\mathcal M'$。

*证明。* 对 $A\in\mathcal M$，有 $AK\subseteq K$ 且 $A^*K\subseteq K$（因为 $A^*\in\mathcal M$）。
第二个包含关系给出 $AK^\perp\subseteq K^\perp$：对 $y\in K^\perp$，$x\in K$，
$\langle x,Ay\rangle=\langle A^*x,y\rangle=0$。因此 $A$ 把 $K$ 映入 $K$，把 $K^\perp$
映入 $K^\perp$，这意味着 $AP=PAP=PA$。$\square$

<a id="r-2-2"></a>

**定理 2.2（双交换子定理）。** $\mathbb C^n$ 上的每个 $*$-代数 $\mathcal M$ 都满足
$\mathcal M=\mathcal M''$。

*证明。* 已知 $\mathcal M\subseteq\mathcal M''$；设 $T\in\mathcal M''$。

*第 1 步：一个向量。* 对 $x\in\mathbb C^n$，子空间 $K=\mathcal Mx$ 在 $\mathcal M$ 下不变，
且包含 $x=1x$。由引理 2.1，它的投影 $P$ 属于 $\mathcal M'$，
所以 $TP=PT$，且 $Tx=TPx=PTx\in K$：存在 $A\in\mathcal M$ 使得 $Tx=Ax$。

*第 2 步：同时处理所有向量。* 在 $\mathbb C^n\otimes\mathbb C^n=(\mathbb C^n)^{\oplus n}$ 上
考虑 $*$-代数 $\tilde{\mathcal M}=\{A\oplus\dots\oplus A:A\in\mathcal M\}$。分块矩阵
$[X_{ij}]_{i,j=1}^n$ 与每个 $A\oplus\dots\oplus A$ 可交换当且仅当对所有 $i,j$ 有 $X_{ij}A=AX_{ij}$，
所以 $\tilde{\mathcal M}'=\{[X_{ij}]:X_{ij}\in\mathcal M'\}$。算子
$\tilde T=T\oplus\dots\oplus T$ 与每个这样的 $[X_{ij}]$ 可交换，因为
$TX_{ij}=X_{ij}T$。因此 $\tilde T\in\tilde{\mathcal M}''$。对 $\tilde{\mathcal M}$、
$\tilde T$ 和向量 $\xi=(e_1,\dots,e_n)$ 应用第 1 步：存在 $A\in\mathcal M$ 使得
$(Te_1,\dots,Te_n)=(Ae_1,\dots,Ae_n)$，即 $T=A\in\mathcal M$。$\square$

*计算：* 定理 2.2 的一个例子——[教程，2.4 节](../../tutorial/zh/02_von_neumann_algebras.md#s-2-4)。

所以在有限维情形，"冯·诺依曼代数"与"含单位元的矩阵 $*$-代数"是同一回事。尽管如此，
双交换子定理仍然提供了正确的视角：在物理上，可观测量代数常常*以交换子的形式给出*——
即与某个对称性可交换的算子全体（2.4 节）。

<a id="s-2-2"></a>

## 2.2 中心

<a id="r-2-3"></a>

**引理 2.3.** 设 $\mathcal Z$ 是 $\mathbb C^n$ 上的交换 $*$-代数。则存在
非零、两两正交且满足 $\sum_cP_c=1$ 的投影 $P_1,\dots,P_C\in\mathcal Z$，
使得 $\mathcal Z=\mathrm{span}\{P_1,\dots,P_C\}$。它们是 $\mathcal Z$ 的
极小非零投影，并且在不计顺序的意义下唯一。

*证明。* 每个 $z\in\mathcal Z$ 都是正规的（$z^*\in\mathcal Z$ 与 $z$ 可交换）。设
$z_1,\dots,z_r$ 张成 $\mathcal Z$。我们把 $\mathbb C^n$ 分解为公共特征子空间。$z_1$ 的
特征子空间两两正交并张成 $\mathbb C^n$（定理 1.3）；它们在 $z_2,\dots,z_r$ 及其伴随（与
$z_1$ 可交换）下不变，所以每个 $z_j$ 限制在 $z_1$ 的特征子空间上仍是正规的，于是可以把
$z_1$ 的每个特征子空间再分解为 $z_2$ 的特征子空间，依此类推。结果得到一个正交分解
$\mathbb C^n=E_1\oplus\dots\oplus E_C$，其中各子空间非零，每个
$z\in\mathcal Z$ 在其上作用为一个数乘 $\chi_c(z)$；这里我们已把对应于同一函数
$\chi_c:\mathcal Z\to\mathbb C$ 的子空间合并，使得当 $c\neq d$ 时 $\chi_c\neq\chi_d$。令
$P_c$ 为到 $E_c$ 上的投影。则对每个 $z\in\mathcal Z$ 有 $z=\sum_c\chi_c(z)P_c$。

还需证明 $P_c\in\mathcal Z$。对 $d\neq c$，选取 $z_{cd}\in\mathcal Z$ 使得
$\chi_c(z_{cd})\neq\chi_d(z_{cd})$。乘积
$\prod_{d\neq c}\big(z_{cd}-\chi_d(z_{cd})1\big)/\big(\chi_c(z_{cd})-\chi_d(z_{cd})\big)$
属于 $\mathcal Z$，在 $E_c$ 上作用为 $1$，在每个 $E_d$（$d\neq c$）上作用为 $0$；因此它就是
$P_c$（当 $C=1$ 时 $P_1=1$：空积等于 $1$）。极小性与唯一性：投影 $\sum_c a_cP_c$ 满足 $a_c\in\{0,1\}$，所以
$\mathcal Z$ 的非零投影恰好是若干个（非空集合的）$P_c$ 之和。$\square$

将此引理应用于 $\mathcal Z=Z(\mathcal M)$，投影 $P_c$ 称为 $\mathcal M$ 的**极小中心
投影**。每个 $A\in\mathcal M$ 可分解为
$A=\sum_cAP_c$，其中 $AP_c=P_cAP_c$，所以
$\mathcal M=\bigoplus_c\mathcal MP_c$（作用在
$E_c=P_c\mathbb C^n$ 上的代数的正交直和）。

*库中：* 中心元素 $\sum_cz_cP_c$（其中包括极小中心投影）由 `TypeIAlgebra.central` 构造——[教程，2.5 节](../../tutorial/zh/02_von_neumann_algebras.md#s-2-5)。

<a id="r-2-4"></a>

**引理 2.4.** 每个 $\mathcal M_c=\mathcal MP_c$，作为 $E_c$ 上以 $P_c$ 为单位元的代数，都是
因子。

*证明。* 设 $y\in\mathcal M_c$ 与 $\mathcal M_c$ 可交换。对 $A\in\mathcal M$，
$Ay=AP_cy=(AP_c)y=y(AP_c)=yP_cA=yA$，这里用到 $P_c$ 是中心元以及 $y=yP_c=P_cy$。
所以 $y\in Z(\mathcal M)$，$y=\sum_d\chi_d(y)P_d$，从而 $y=yP_c=\chi_c(y)P_c$。$\square$

<a id="s-2-3"></a>

## 2.3 结构定理

<a id="r-2-5"></a>

**引理 2.5（因子）。** 设 $\mathcal N$ 是空间 $E$ 上的因子。则存在整数
$k,m\ge1$ 和酉映射 $V:\mathbb C^k\otimes\mathbb C^m\to E$，使得
$V^*\mathcal NV=M_k(\mathbb C)\otimes1_m$。此外 $V^*\mathcal N'V=1_k\otimes M_m(\mathbb C)$。

*证明。* *第 1 步：极小投影。* $\mathcal N$ 含有非零投影（例如 $1$）。
设 $e\in\mathcal N$ 是秩最小的非零投影。我们断言
$e\mathcal Ne=\mathbb Ce$。设 $x=x^*\in e\mathcal Ne$。$x$ 对应于非零特征值 $\lambda$ 的每个
谱投影 $q$ 都可写成 $q=p(x)$，其中 $p$ 是满足
$p(\lambda)=1$ 且对 $\mu\in(\mathrm{spec}\,x\cup\{0\})\setminus\{\lambda\}$ 有 $p(\mu)=0$
的 Lagrange 多项式（定理 1.5(a),(b)：$p$ 在 $\mathrm{spec}\,x$ 上与 $\{\lambda\}$ 的示性函数一致）。
由于 $p(0)=0$，对某个多项式 $r$ 有 $p(t)=t\,r(t)$，所以 $q=x\,r(x)=r(x)\,x\in\mathcal N$，而 $x=ex=xe$ 给出
$q=eq=qe$；所以 $q$ 的值域包含在 $e$ 的值域中。若 $q\neq0$，由极小性得
$\mathrm{rank}\,q\ge\mathrm{rank}\,e$，所以两个值域重合，$q=e$。由于对应于不同特征值的谱投影
相互正交，$x$ 至多有一个非零特征值 $\lambda$，其谱投影为 $e$，所以 $x=\lambda e$（若没有非零特征值，
则 $x=0=0\cdot e$）。$e\mathcal Ne$ 的一般元素为 $x_1+ix_2$，其中
$x_{1,2}=x_{1,2}^*\in e\mathcal Ne$。

*第 2 步：$e$ 生成 $E$。* 令 $F=\mathrm{span}\,\mathcal Ne E$。它在 $\mathcal N$ 下
不变，所以它的投影 $p_F$ 属于 $\mathcal N'$（引理 2.1）。它在 $*$-代数 $\mathcal N'$ 下也
不变：对 $b\in\mathcal N'$，$bAe\xi=Aeb\xi\in F$；所以对 $\mathcal N'$ 应用引理 2.1，
得 $p_F\in\mathcal N''=\mathcal N$（定理 2.2）。因此
$p_F\in Z(\mathcal N)=\mathbb C1$，又 $p_F\neq0$（因为 $eE\subseteq F$），所以 $F=E$。

*第 3 步：矩阵单位。* 称满足 $v^*v=e$ 的 $v\in\mathcal N$ 为 *$e$-等距*；此时 $vv^*$
是 $\mathcal N$ 中的投影。（首先 $v=ve$：对每个 $x$，
$\|v(1-e)x\|^2=\langle(1-e)x,v^*v(1-e)x\rangle=\langle(1-e)x,e(1-e)x\rangle=0$。因此
$vv^*v=ve=v$，且 $(vv^*)^2=vev^*=vv^*$。对下面的 $v_i$，这给出 $e_iv_i=v_i$ 和
$v_i^*e_i=v_i^*$。）选取
值域 $e_i=v_iv_i^*$ 两两正交的 $e$-等距 $v_1=e,v_2,\dots,v_k$，使 $k$
最大（这些值域是维数都为 $\mathrm{rank}\,e$ 的两两正交子空间，
所以 $k\le\dim E/\mathrm{rank}\,e$）。我们断言 $\sum_ie_i=1$。否则 $g=1-\sum_ie_i\neq0$ 是
$\mathcal N$ 中的投影。由第 2 步，对某个 $A\in\mathcal N$ 有 $gAe\neq0$（若对所有 $A$ 都有
$gAe=0$，则对所有 $A,\xi$ 有 $gAe\xi=0$，于是 $g$ 在 $F=E$ 上为零）。令 $y=gAe$。则 $y^*y\in e\mathcal Ne$ 是正的且非零，
所以 $y^*y=\lambda e$，$\lambda>0$，而 $v=y/\sqrt\lambda$ 是一个值域含于 $gE$、
与所有 $e_i$ 正交的 $e$-等距——这与 $k$ 的最大性矛盾。

对 $i\neq j$，$v_i^*v_j=v_i^*e_ie_jv_j=0$。令 $E_{ij}=v_iv_j^*$。则
$E_{ij}E_{kl}=v_iv_j^*v_kv_l^*=\delta_{jk}v_iev_l^*=\delta_{jk}E_{il}$，
$E_{ij}^*=E_{ji}$，且 $\sum_iE_{ii}=1$。对 $A\in\mathcal N$，
$A=\sum_{ij}e_iAe_j=\sum_{ij}v_i(v_i^*Av_j)v_j^*$，而 $v_i^*Av_j\in e\mathcal Ne=\mathbb Ce$，
所以 $A=\sum_{ij}a_{ij}E_{ij}$，$a_{ij}\in\mathbb C$。于是 $\mathcal N=\mathrm{span}\{E_{ij}\}$。

*第 4 步：酉映射。* 令 $m=\mathrm{rank}\,e$，$f_1,\dots,f_m$ 为 $eE$ 的一组规范正交基。
向量 $v_if_r$ 是规范正交的：
$\langle v_if_r,v_jf_s\rangle=\langle f_r,v_i^*v_jf_s\rangle=\delta_{ij}\langle f_r,ef_s\rangle=\delta_{ij}\delta_{rs}$。
它们张成 $E$，因为 $E=\bigoplus_ie_iE$ 且 $e_iE=v_ieE$。定义 $V(|i\rangle\otimes f_r)=v_if_r$，其中 $|1\rangle,\dots,|k\rangle$ 是
$\mathbb C^k$ 的标准基，并借助 $f_r$ 把 $\mathbb C^m$ 等同于 $eE$。
则 $E_{ij}v_lf_r=\delta_{jl}v_if_r$，即 $V^*E_{ij}V=|i\rangle\langle j|\otimes1_m$，从而
$V^*\mathcal NV=M_k(\mathbb C)\otimes1_m$。

*第 5 步：交换子。* 把 $X\in M_k\otimes M_m$ 写成分块矩阵
$X=\sum_{ij}|i\rangle\langle j|\otimes X_{ij}$。对固定的 $a,b$，
$(|a\rangle\langle b|\otimes1)X=\sum_j|a\rangle\langle j|\otimes X_{bj}$，
$X(|a\rangle\langle b|\otimes1)=\sum_i|i\rangle\langle b|\otimes X_{ia}$。若二者相等，
比较位置 $(a,j)$（$j\neq b$）上的块得 $X_{bj}=0$，比较位置
$(a,b)$ 上的块得 $X_{bb}=X_{aa}$。由于 $a,b$ 任意，$X=1_k\otimes Y$，其中 $Y=X_{11}$；
反之，每个 $1_k\otimes Y$ 都与 $M_k\otimes1_m$ 可交换。$\square$

<a id="r-2-6"></a>

**定理 2.6（有限维冯·诺依曼代数的结构）。** 设
$\mathcal M$ 是 $\mathbb C^n$ 上的冯·诺依曼代数。则存在 $C\ge1$、满足 $\sum_ck_cm_c=n$ 的整数
$k_c,m_c\ge1$，以及酉映射
$V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n$，使得
$$
V^*\mathcal MV=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},\qquad
V^*\mathcal M'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),
$$
并且 $V^*Z(\mathcal M)V=\mathrm{span}\{V^*P_1V,\dots,V^*P_CV\}$，其中 $P_1,\dots,P_C$ 是
$\mathcal M$ 的极小中心投影（引理 2.3），而 $V^*P_cV$ 是到第 $c$ 个直和项上的投影。数对 $(k_c,m_c)$ 在不计顺序的意义下唯一：对*任何*具有第一条性质的酉映射 $V$，
$V(\cdot)V^*$ 都把到各直和项上的投影映为诸 $P_c$（见证明），
所以 $k_c^2=\dim\mathcal MP_c$ 与 $m_c=\mathrm{rank}\,P_c/k_c$ 由 $\mathcal M$ 决定。特别地，作为 $*$-代数有 $\mathcal M\cong\bigoplus_cM_{k_c}(\mathbb C)$，并且
$\dim\mathcal M=\sum_ck_c^2$，$\dim\mathcal M'=\sum_cm_c^2$。

*证明。* 由引理 2.3 及其后的注记，$\mathcal M=\bigoplus_c\mathcal M_c$，其中
$\mathcal M_c$ 作用在 $E_c$ 上，且每个 $\mathcal M_c$ 都是因子（引理 2.4）。引理 2.5
给出酉映射 $V_c$，使得 $V_c^*\mathcal M_cV_c=M_{k_c}\otimes1_{m_c}$；令
$V=\bigoplus_cV_c$。关于交换子：$X\in\mathcal M'$ 与每个 $P_c$ 可交换，所以
$X=\sum_cP_cXP_c$，且 $P_cXP_c$ 与 $\mathcal M_c$ 可交换；反之，任何这样的和都与
$\mathcal M$ 可交换。引理 2.5 确定了 $E_c$ 上的 $\mathcal M_c'$，再用 $V$ 共轭即得
$V^*\mathcal M'V$ 的公式。

*唯一性。* 设 $V$ 是任一满足 $V^*\mathcal MV=\mathcal B:=\bigoplus_{c=1}^{C'}M_{k_c}\otimes1_{m_c}$ 的酉映射，
$Q_c$ 是到第 $c$ 个直和项上的投影。$\mathcal B$ 的元素 $\bigoplus_cX_c$（$X_c\in M_{k_c}\otimes1_{m_c}$）
是中心元当且仅当每个 $X_c$ 都与 $M_{k_c}\otimes1_{m_c}$ 可交换，即（引理 2.5，第 5 步）$X_c\in(M_{k_c}\otimes1_{m_c})\cap(1_{k_c}\otimes M_{m_c})=\mathbb C1$
（元素 $a\otimes1=1\otimes b$ 的所有对角块都等于 $b$，所有非对角块都为 $0$，所以 $a$ 是数量）。
因此 $Z(\mathcal B)=\mathrm{span}\{Q_1,\dots,Q_{C'}\}$，而 $Z(\mathcal M)=VZ(\mathcal B)V^*$ 由非零、
两两正交且和为 $1$ 的投影 $VQ_cV^*$ 张成。它们每一个都是 $Z(\mathcal M)$ 的非零投影，因此（引理 2.3）
是诸 $P_d$ 中某个非空集合 $S_c$ 之和；诸 $S_c$ 两两不交（正交性）且覆盖所有 $d$（和为 $1$），又由于
$C'=\dim Z(\mathcal M)=C$，每个 $S_c$ 都是单元素集。所以在重新编号后 $VQ_cV^*=P_c$，从而
$\dim\mathcal MP_c=\dim(M_{k_c}\otimes1_{m_c})=k_c^2$，$\mathrm{rank}\,P_c=k_cm_c$。$\square$

*库中：* 代数 $\bigoplus_cM_{k_c}(\mathbb C)$ 是一个 `TypeIAlgebra`——[教程，2.5 节](../../tutorial/zh/02_von_neumann_algebras.md#s-2-5)；酉映射 $V$ 与数对 $(k_c,m_c)$ 由 `decompose.wedderburn` 计算（前端接口为 `from_group`、`from_generators`）——[教程，4.3 节](../../tutorial/zh/04_symmetries.md#s-4-3)。

数 $k_c$ 描述代数本身。**重数** $m_c$ 描述代数如何嵌入 $M_n(\mathbb C)$：它们在代数结构中
不可见，但决定了 $\mathbb C^n$ 的迹在 $\mathcal M$ 上的限制，
$$
\mathrm{Tr}_{\mathbb C^n}(A)=\sum_cm_c\mathrm{Tr}(A_c),\qquad A=V\Big(\bigoplus_cA_c\otimes1_{m_c}\Big)V^*.
$$
这就是本库使用一族迹的原因（第 3 章）。同样，对自伴的 $A$，
$A$ 在 $V(\mathbb C^{k_c}\otimes\mathbb C^{m_c})$ 上的限制具有块 $A_c$ 的特征值，每个特征值的重数是它在
$A_c$ 中重数的 $m_c$ 倍，因为 $A_c\otimes1_{m_c}$ 就是 $A_c$ 的 $m_c$ 个拷贝。

*库中：* 迹 `Tr_blunt`、`Tr_norm`、`tau_vN`——[教程，3.1 节](../../tutorial/zh/03_traces_and_states.md#s-3-1)。

> **物理意义**（物理解释，而非本课程的定理）。$P_c$ 是**超选择扇区**或**对称扇区**：$\mathcal M$ 中没有可观测量连接不同的
> 扇区，因此"处于哪个扇区"的测量与其他任何测量都是相容的。中心是系统的*经典*部分。
> 重数 $m_c$ 是 $\mathcal M$ 中每个可观测量的一种简并度。

<a id="s-2-4"></a>

## 2.4 例子：对称群及其交换子

设有限群 $G$ 通过酉表示 $\pi$ 作用在 $\mathbb C^n$ 上。与该对称性可交换的算子构成
$\pi(G)'$，这是一个冯·诺依曼代数；它的交换子 $\pi(G)''$ 是诸 $\pi(g)$ 的线性张成（对该
张成应用定理 2.2，它是一个 $*$-代数：包含 $\pi(e)=1$，对乘积 $\pi(g)\pi(h)=\pi(gh)$ 和伴随
$\pi(g)^*=\pi(g^{-1})$ 封闭）。

<a id="r-2-7"></a>

**命题 2.7.** 对 $\mathcal M=\pi(G)''$ 应用定理 2.6。则
$V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$，其中每个 $\rho_c$ 是 $G$ 在 $\mathbb C^{k_c}$
上的不可约酉表示，且当 $c\neq d$ 时 $\rho_c,\rho_d$ 不等价。因此
$$
\pi\cong\bigoplus_c\rho_c\otimes1_{m_c},\qquad
\pi(G)''\cong\bigoplus_cM_{\dim\rho_c},\qquad
\pi(G)'\cong\bigoplus_cM_{m_c}:
$$
对称算子 $\pi(G)'$ 的块大小等于各不可约表示的重数，并且在一个代数与其交换子之间，
$k_c$ 与 $m_c$ 的角色互换。

*证明。* $\pi(g)\in\mathcal M$，所以对某些
$\rho_c(g)\in M_{k_c}$ 有 $V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$，并且由于 $\pi$ 是酉表示，
$g\mapsto\rho_c(g)$ 也是酉表示。由于诸 $\pi(g)$ 张成 $\mathcal M$，诸 $\rho_c(g)$ 张成 $M_{k_c}$；
更强地，映射
$\sum_ga_g\pi(g)\mapsto\big(\sum_ga_g\rho_c(g)\big)_c$ 是到 $\bigoplus_cM_{k_c}$ 上的满射。
*不可约性：* $\mathbb C^{k_c}$ 中在所有 $\rho_c(g)$ 下不变的子空间在它们的张成 $M_{k_c}$ 下
也不变，因此为 $0$ 或 $\mathbb C^{k_c}$。*不等价性：* 若对所有 $g$ 和某个可逆的 $T$，
$c\neq d$，有 $\rho_d(g)=T\rho_c(g)T^{-1}$，则像中的每个元素
$(x_{c'})_{c'}$ 都满足 $x_d=Tx_cT^{-1}$，这就排除了满足
$x_c=1$、$x_d=0$ 的元素——与满射性矛盾。最后两个同构即
定理 2.6。$\square$

*库中：* 有限群表示的分解，`decompose.from_group`——[教程，4.3 节](../../tutorial/zh/04_symmetries.md#s-4-3)。

表示分解为不可约表示的方式在等价意义下唯一，从而 $m_c$ 是任何分解中*唯一确定的*重数。
*证明。* 设 $\rho,\rho'$ 不可约，$T$ 是一个缠结算子，即对所有 $g$ 有 $T\rho(g)=\rho'(g)T$。$T$ 的核与
像都是不变子空间，所以 $T=0$ 或 $T$ 可逆（Schur 引理）；若 $\rho'=\rho$，对 $T$ 的一个特征值
$\lambda$，同样的结论适用于 $T-\lambda1$，故 $T=\lambda1$。因此从 $\rho$ 到 $\rho'$ 的缠结算子构成的
空间，当 $\rho\cong\rho'$ 时维数为 $1$（对一个可逆缠结算子 $S$，$T\mapsto S^{-1}T$
把它们与从 $\rho$ 到自身的缠结算子等同起来，而后者是 $\mathbb C1$），否则为 $0$（非零的缠结算子
会是可逆的，即一个等价）。若 $\pi\cong\bigoplus_j\sigma_j$，其中 $\sigma_j$
不可约（允许重复），则从不可约表示 $\rho$ 到 $\bigoplus_j\sigma_j$ 的缠结算子是由从 $\rho$ 到 $\sigma_j$ 的
缠结算子组成的组 $(T_j)_j$，而与该等价复合不改变维数；所以从 $\rho$ 到 $\pi$ 的缠结算子构成维数为 $\#\{j:\sigma_j\cong\rho\}$
的空间，这个数只依赖于 $\pi$ 和 $\rho$。对于命题 2.7 的分解，其中各 $\rho_c$ 两两不等价，这个数在
$\rho\cong\rho_c$ 时为 $m_c$，在 $\rho$ 不与任何 $\rho_c$ 等价时为 $0$。$\square$

我们取 $G=S_3$，它置换 $(\mathbb C^2)^{\otimes3}$ 的三个张量因子，即三个
量子比特或三个自旋 $\tfrac12$：对 $x\in\{0,1\}^3$，$\pi(\sigma)|x\rangle=|\sigma\cdot x\rangle$，其中
$(\sigma\cdot x)_{\sigma(i)}=x_i$。令 $F=\{v:\pi(g)v=v\ \forall g\}$ 为**对称子空间**。其
分解为：$F$（维数 $4$）承载重数为 $4$ 的平凡表示；它的
正交补（维数 $4$）承载重数为 $2$ 的二维**标准表示**（$S_3$ 在 $\mathbb C^3$ 中坐标和为 $0$ 的向量上的
置换表示）；符号表示不出现。所以
$\pi(S_3)''\cong M_1\oplus M_2$（维数 $5$），$\pi(S_3)'\cong M_4\oplus M_2$
（维数 $20$）。

*证明。* 写 $V^*\pi(S_3)''V=\bigoplus_cM_{k_c}\otimes1_{m_c}$（定理 2.6），于是
$\sum_ck_cm_c=8$，$\sum_cm_c^2=\dim\pi(S_3)'$。

*第 1 步：$\dim\pi(S_3)'=20$。* 由于 $\pi(g)^*|g\cdot x\rangle=|x\rangle$，矩阵
$\pi(g)X\pi(g)^*$ 在位置 $(g\cdot x,g\cdot y)$ 上的元素为 $X_{xy}$。因此 $X\in\pi(S_3)'$，即对所有 $g$ 有
$\pi(g)X\pi(g)^*=X$，当且仅当 $(x,y)\mapsto X_{xy}$ 在 $S_3$ 通过 $g\cdot(x,y)=(g\cdot x,g\cdot y)$ 作用于
数对 $(x,y)$ 的轨道上为常数；各轨道的示性矩阵构成
$\pi(S_3)'$ 的一组基，所以它的维数等于轨道的个数。数对 $(x,y)$ 就是四元集合 $\{0,1\}^2$ 中元素的三元组
$\big((x_1,y_1),(x_2,y_2),(x_3,y_3)\big)$，而 $S_3$ 置换这三个分量，所以一个轨道就是从四个元素中取
大小为 $3$ 的多重集：形如 $\{a,a,a\}$ 的有 $4$ 个，形如 $\{a,a,b\}$（$b\neq a$）的有 $4\cdot3=12$ 个，
形如 $\{a,b,c\}$（元素互不相同）的有 $4$ 个，共 $20$ 个。

*第 2 步：$k=(1,2)$，$m=(4,2)$。* 由于 $k_c\ge1$，$\sum_cm_c\le\sum_ck_cm_c=8$。令 $M=\max_cm_c$。
则 $20=\sum_cm_c^2\le M\sum_cm_c\le8M$ 且 $M^2\le20$，所以 $M\in\{3,4\}$。若 $M=4$，其余
$m_c$ 的平方和为 $4$，所以它们是 $(2)$ 或 $(1,1,1,1)$。若 $M=3$，其余 $m_c$ 的平方和为 $11$，
和至多为 $5$；若它们都 $\le2$，则其平方和至多为 $2\cdot5=10$，所以其中有一个等于 $3$，而剩下的
平方和为 $2$：它们是 $(1,1)$。
因此在不计顺序的意义下，$m$ 是 $(4,2)$、$(4,1,1,1,1)$ 或 $(3,3,1,1)$。在后两种情形，
$\sum_cm_c=8=\sum_ck_cm_c$ 迫使所有 $k_c=1$，于是 $\pi(S_3)''\cong\bigoplus_cM_1$ 将是
交换的，而对换 $(12)$ 与 $(23)$ 并不交换：
$\pi((12))\pi((23))|001\rangle=|100\rangle$，而 $\pi((23))\pi((12))|001\rangle=|010\rangle$。所以
$m=(4,2)$（把直和项编号使得 $m_1=4$），而由 $4k_1+2k_2=8$ 及 $k_1,k_2\ge1$ 得 $k=(1,2)$。特别地，
$\dim\pi(S_3)''=\sum_ck_c^2=5$。

*第 3 步：平凡表示与符号表示。* $\rho_2$ 不固定任何非零向量 $u$：否则直线
$\mathbb Cu$ 将是 $\mathbb C^2$ 的不变子空间，与不可约性矛盾
（命题 2.7）。写 $V^*v=v_1\oplus v_2$，其中 $v_1\in\mathbb C^1\otimes\mathbb C^4$，
$v_2=\sum_{r=1}^2u_r\otimes f_r\in\mathbb C^2\otimes\mathbb C^2$（$f_1,f_2$ 为标准基）。则
$v\in F$ 当且仅当对所有 $g$ 和 $r$ 有 $\rho_1(g)v_1=v_1$ 且 $\rho_2(g)u_r=u_r$，所以每个 $v\in F$ 都有
$v_2=0$。非零向量 $|000\rangle$ 属于 $F$，所以数量 $\rho_1(g)$ 固定一个非零向量：
$\rho_1$ 是平凡表示，$F=V(\mathbb C^1\otimes\mathbb C^4)$ 的维数为 $4$，而它的正交
补是 $V(\mathbb C^2\otimes\mathbb C^2)$。符号表示既不与
$\rho_1$ 等价（因为 $\mathrm{sgn}(12)=-1$），也不与 $\rho_2$ 等价（维数不同），所以由上面证明的
重数唯一性，它出现 $0$ 次。

*第 4 步：$\rho_2$ 是标准表示。* 向量 $e_1=|100\rangle$、$e_2=|010\rangle$、
$e_3=|001\rangle$ 在诸 $\pi(\sigma)$ 下按 $\sigma$ 置换位置的方式被置换，所以
$W=\{\sum_ix_ie_i:\sum_ix_i=0\}$ 是不变的，且 $\pi$ 在其上的作用就是标准表示。
$W\perp F$：对 $v\in F$ 和 $w=\sum_ix_ie_i\in W$，对每个 $g$ 有
$\langle v,w\rangle=\langle\pi(g)v,w\rangle=\langle v,\pi(g)^{-1}w\rangle$，
所以 $\langle v,w\rangle=\langle v,\frac16\sum_g\pi(g)w\rangle$，而
$\frac16\sum_g\pi(g)w=\frac13\big(\sum_ix_i\big)(e_1+e_2+e_3)=0$（六个置换中恰有两个把每个 $e_i$ 映到
每个 $e_j$）。因此 $V^*W\subseteq\mathbb C^2\otimes\mathbb C^2$。映射
$T_r:\mathbb C^2\otimes\mathbb C^2\to\mathbb C^2$，$u\otimes f\mapsto\langle f_r,f\rangle u$，满足
$T_r(\rho_2(g)\otimes1)=\rho_2(g)T_r$。选取 $r$ 使 $T_rV^*|_W\neq0$（若两者都为零，则 $V^*W$
将为 $0$）。于是 $T_rV^*|_W$ 是从 $\pi|_W$ 到 $\rho_2$ 的非零缠结算子；它的像是 $\mathbb C^2$ 的非零
不变子空间，因而是 $\mathbb C^2$，又由于 $\dim W=2$，它是双射。所以 $\rho_2$
等价于标准表示（顺便可知标准表示是不可约的）。$\square$

*计算：* 维数 $20$ 与 $5$，以及定理 2.2 在本例中的验证——[教程，2.4 节](../../tutorial/zh/02_von_neumann_algebras.md#s-2-4)。*库中：* `decompose.permutation_representation`、`decompose.commutant_dimension`——[教程，4.2 节](../../tutorial/zh/04_symmetries.md#s-4-2)。

重数在谱中是可见的。每个自伴的 $K\in\pi(S_3)'$ 都具有形式
$V^*KV=(1_1\otimes\kappa_1)\oplus(1_2\otimes\kappa_2)$，其中 $\kappa_1\in M_4$、$\kappa_2\in M_2$ 自伴
（定理 2.6）。在 $F$ 上，它的特征值就是 $\kappa_1$ 的特征值；在 $F$ 的补空间上，它们是 $\kappa_2$ 的
特征值，每个重复两次，因为 $1_2\otimes\kappa_2$ 就是 $\kappa_2$ 的两个拷贝。所以 $K$ 在 $F$ 的补空间上的
*每个*特征值都具有偶数重数。对随机的 $K$，结论更强。对群取平均，
$E(X)=\frac1{|G|}\sum_g\pi(g)X\pi(g)^*$，把 $M_8(\mathbb C)$ 映满交换子，并且在交换子上是
恒等映射（习题 3），所以对 Ginibre 矩阵 $Z$（元素为相互独立的标准复高斯变量），
$K=E(Z+Z^*)$ 在自伴元素构成的实向量空间 $\pi(S_3)'_h$ 上有密度
（定理 A.15 及其后的注记）。由命题 4.7(a)（这是向后文的引用；其证明用到
定理 2.6 和附录 A，而不用本段），以概率一，$\kappa_1$ 与 $\kappa_2$ 的六个
特征值两两不同：此时 $K$ 在 $F$ 上有 $4$ 个单重特征值，在其补空间上有
$2$ 个二重简并特征值。最后，$P_F=\frac1{|G|}\sum_g\pi(g)$ 是到 $F$ 上的正交投影：对所有 $h$ 有 $\pi(h)P_F=P_F$，所以它的
值域含于 $F$；对 $v\in F$ 有 $P_Fv=v$；并且通过代换 $g\to g^{-1}$ 得 $P_F^*=P_F$。

*计算：* 随机 $K$ 的二重简并特征值——[教程，2.4 节](../../tutorial/zh/02_von_neumann_algebras.md#s-2-4)。

<a id="s-2-5"></a>

## 2.5 库中的代数

*库中：* `TypeIAlgebra` 如何存储 $\bigoplus_cM_{k_c}(\mathbb C)$、它的算子及其中心——[教程，2.5 节](../../tutorial/zh/02_von_neumann_algebras.md#s-2-5)。

<a id="s-2-6"></a>

## 2.6 条件期望

从 $\mathcal M$ 到冯·诺依曼子代数
$\mathcal N\subseteq\mathcal M$（包含于 $\mathcal M$ 且具有相同单位元的 $*$-代数）上的**条件期望**是一个线性映射 $E:\mathcal M\to\mathcal N$，
它是正的（$A\ge0\Rightarrow E(A)\ge0$），对 $B\in\mathcal N$ 满足 $E(B)=B$，并且是
$\mathcal N$-双模映射：对 $B_1,B_2\in\mathcal N$ 有 $E(B_1AB_2)=B_1E(A)B_2$。

<a id="r-2-8"></a>

**命题 2.8.** 对 $\bigoplus_c\mathbb C^{k_c}$ 上的 $\mathcal M=\bigoplus_cM_{k_c}$（$P_c$ 为到第 $c$ 个块上的投影）
和 $\mathcal N=Z(\mathcal M)=\mathrm{span}\{P_1,\dots,P_C\}$（习题 1），
$$
E(A)=\sum_c\frac{\mathrm{Tr}\,A_c}{k_c}\,P_c
$$
是一个条件期望。它保持每个形如
$\sum_cw_c\mathrm{Tr}\,A_c$ 的迹。

*证明。* 线性是显然的；$A\ge0$ 蕴含 $\mathrm{Tr}\,A_c=\sum_i\langle e_i,Ae_i\rangle\ge0$（对第 $c$ 个块的
标准基向量 $e_i$ 求和），所以 $E(A)\ge0$。对
$B=\sum_cb_cP_c$，$E(B)=\sum_cb_ck_cP_c/k_c=B$。对 $B_i=\sum_cb^{(i)}_cP_c$，
$\mathrm{Tr}(B_1AB_2)_c=b^{(1)}_cb^{(2)}_c\mathrm{Tr}\,A_c$，由此得双模
性质。最后，对每个 $c$ 有 $\mathrm{Tr}(E(A)_c)=\mathrm{Tr}\,A_c$。$\square$

*库中：* $E$ 就是 `channels.center_expectation`——[教程，2.6 节](../../tutorial/zh/02_von_neumann_algebras.md#s-2-6)。

2.4 节中的群平均是 $M_n(\mathbb C)$ 中到交换子
$\pi(G)'$ 上的条件期望，并且保持 $\mathrm{Tr}$（习题 3）。类似地，在完整的
矩阵代数 $M_D(\mathbb C)$, $D=\sum_ck_c$ 上（它以其分块对角矩阵的形式包含 $\mathcal M=\bigoplus_cM_{k_c}$），
映射 $E_{\mathcal M}(X)=\sum_cP_cXP_c$ 是到 $\mathcal M$ 上保持 $\mathrm{Tr}$ 的条件期望。
*证明。* $P_cXP_c$ 在第 $c$ 个对角块之外为零，所以 $E_{\mathcal M}(X)\in\mathcal M$。对 $X\ge0$，
$\langle v,P_cXP_cv\rangle=\langle P_cv,XP_cv\rangle\ge0$。元素
$B\in\mathcal M$ 是分块对角的，所以 $BP_c=P_cB$ 且 $B=\sum_cP_cBP_c$；因此 $E_{\mathcal M}(B)=B$，并且对
$B_1,B_2\in\mathcal M$ 有 $E_{\mathcal M}(B_1XB_2)=\sum_cB_1P_cXP_cB_2=B_1E_{\mathcal M}(X)B_2$。最后
$\sum_c\mathrm{Tr}(P_cXP_c)=\sum_c\mathrm{Tr}(P_cX)=\mathrm{Tr}\,X$。$\square$

*库中：* 群平均是 `decompose.twirl`——[教程，4.3 节](../../tutorial/zh/04_symmetries.md#s-4-3)。

对 $\mathcal M$ 上的**态** $\omega$（满足对所有 $A$ 有 $\omega(A^*A)\ge0$ 且 $\omega(1)=1$ 的线性泛函；
态在第 3 章中研究），由命题 2.8 和线性性得
$\omega(E(A))=\sum_c\omega(P_c)\mathrm{Tr}(A_c)/k_c$，并且在中心上 $\omega\circ E=\omega$，因为在那里
$E(z)=z$。（条件期望一般不唯一：到 $\mathbb C1$ 上，每个态
$\varphi$ 都给出一个条件期望 $A\mapsto\varphi(A)1$；在 $M_2(\mathbb C)$ 上，态 $A\mapsto A_{11}$ 与 $A\mapsto A_{22}$
给出两个不同的条件期望。）它们都是同一思想的实例：
投影到"看不见"某种结构的算子上。

> **物理意义**（物理解释，而非本课程的定理）。映射 $E_{\mathcal M}$ 抹去超选择扇区之间的相干性，
> 这正是测量荷（扇区）而不记录结果时发生的事情。到中心上的条件期望 $E$ 走得更远：它只保留
> 经典信息，即每个扇区的权重 $\omega(P_c)$，并把每个扇区内部的态换成最大混合态（归一化的迹
> $\mathrm{Tr}(A_c)/k_c$）。

<a id="s-2-7"></a>

## 2.7 可与之比较的工具

*计算：* 与其他软件（GAP、SageMath、QuSpin、QuTiP）的比较——[教程，2.7 节](../../tutorial/zh/02_von_neumann_algebras.md#s-2-7)。

<a id="x-2"></a>

## 习题

1. 证明 $M_k(\mathbb C)$ 的中心是 $\mathbb C1$（利用矩阵单位
   $|i\rangle\langle j|$），并证明 $Z(\bigoplus_cM_{k_c})$ 的维数为 $C$。
2. 设 $S$ 是 $\mathbb C^3$ 上的循环移位，$Se_i=e_{i+1\bmod3}$。确定
   $\{S\}'$ 与 $\{S\}''$。
3. 对有限群的酉表示 $\pi$，证明
   $E(X)=\frac1{|G|}\sum_g\pi(g)X\pi(g)^*$ 是从 $M_n(\mathbb C)$ 到 $\pi(G)'$ 上的
   条件期望，并且保持 $\mathrm{Tr}$。

<details>
<summary>解答</summary>

1. 若 $X$ 与所有 $|i\rangle\langle j|$ 可交换，则
   $X|i\rangle\langle j|=|i\rangle\langle j|X$；把两边作用于 $e_j$ 并比较，
   得对所有 $i,j$ 有 $Xe_i=X_{jj}e_i$，所以 $X=X_{11}1$。在直和中，元素
   $\bigoplus X_c$ 是中心元当且仅当每个 $X_c$ 都是，所以中心为
   $\{\bigoplus_cz_c1_{k_c}\}$，维数为 $C$。

2. $S$ 是酉矩阵，特征值为 $1,\omega,\omega^2$，$\omega=e^{2\pi i/3}$，均为单重：
   向量 $u_j=\sum_i\omega^{-ij}e_i$（$j=0,1,2$）满足 $Su_j=\omega^ju_j$，且相互正交。
   与 $S$ 可交换的算子保持 $S$ 的特征子空间，所以在特征基下它是
   对角的，而每个对角矩阵都与 $S$ 可交换：$\{S\}'\cong\mathbb C^3$。它包含
   $1,S,S^2$，在特征基下它们是 $\mathrm{diag}(\lambda^0)$、$\mathrm{diag}(\lambda)$、
   $\mathrm{diag}(\lambda^2)$，$\lambda=(1,\omega,\omega^2)$，它们线性无关（Vandermonde：特征值
   互不相同）；所以 $\{S\}'=\mathrm{span}\{1,S,S^2\}$，即循环矩阵全体。于是 $\{S\}''$ 也是
   $\{S\}'$，因为在 $S$ 的特征基下，与所有对角矩阵可交换的矩阵是对角的（对对角元
   互不相同的 $D$，比较 $XD$ 与 $DX$ 的元素）。

3. $E$ 是线性且正的（它是映射 $X\mapsto UXU^*$ 之和）。对 $h\in G$，
   $\pi(h)E(X)\pi(h)^*=\frac1{|G|}\sum_g\pi(hg)X\pi(hg)^*=E(X)$，所以 $E(X)\in\pi(G)'$。若
   $X\in\pi(G)'$，则每一项都等于 $X$，所以 $E(X)=X$。对 $B_1,B_2\in\pi(G)'$，
   $\pi(g)B_1XB_2\pi(g)^*=B_1\pi(g)X\pi(g)^*B_2$，由此得双模性质。
   最后，对每个 $g$ 有 $\mathrm{Tr}(\pi(g)X\pi(g)^*)=\mathrm{Tr}\,X$，所以
   $\mathrm{Tr}\,E(X)=\mathrm{Tr}\,X$。

*计算：* 习题 2 的数值检验——[教程，习题](../../tutorial/zh/02_von_neumann_algebras.md#x-2)。

</details>

**下一章：** [第 3 章](03_traces_and_states.md)对
$\bigoplus_cM_{k_c}$ 上的迹进行分类，引入态与密度，并证明后文要用到的熵不等式与
迹不等式。
