[English](../02_von_neumann_algebras.md) · [Español](../es/02_von_neumann_algebras.md) · [Français](../fr/02_von_neumann_algebras.md) · [Deutsch](../de/02_von_neumann_algebras.md) · **中文** · [日本語](../ja/02_von_neumann_algebras.md) · [Русский](../ru/02_von_neumann_algebras.md)

# 2. 有限维冯·诺依曼代数

冯·诺依曼代数是在"取两次交换子"运算下封闭的算子代数。在有限维情形，这是一个纯代数的
条件，而本章的主要定理断言：在相差一个酉基变换的意义下，每个这样的代数都是若干全矩阵代数
的直和，每个全矩阵代数以一定的重数重复出现。这正是本库所依据的结构。

**你将学到**双交换子定理、带完整证明的结构定理、对称群如何产生这样的代数、什么是中心，
以及本库如何表示这一切。

在本章中，*$\mathbb C^n$ 上的 $*$-代数*指的是一个复线性子空间
$\mathcal M\subseteq M_n(\mathbb C)$，它在乘积和伴随运算下封闭，并且包含
单位元 $1$。

## 2.1 交换子与双交换子定理

对集合 $S\subseteq M_n(\mathbb C)$，其**交换子**为
$S'=\{X\in M_n(\mathbb C):XA=AX\text{ for all }A\in S\}$（即与 $S$ 中所有元素可交换的矩阵之集）。
它总是一个含单位元的代数；若 $S$ 在伴随运算下封闭（$S^*=S$），则 $S'$ 是一个 $*$-代数，因为
对所有 $A\in S$ 有 $XA=AX$ 蕴含对所有 $A\in S$ 有 $X^*A^*=A^*X^*$。$S$ 的每个元素都与 $S'$
的每个元素可交换，所以 $S\subseteq S''$；把这一点应用于 $S'$ 得到
$S'\subseteq S'''$，而对 $S\subseteq S''$ 取交换子（取交换子使包含关系反向）
得到 $S'''\subseteq S'$。因此 $S'=S'''$。

当一个代数作用在 $\mathbb C^n$ 的子空间 $E$ 上时，我们借助一组规范正交基把 $E$ 等同于
$\mathbb C^{\dim E}$，从而在 $\mathbb C^n$ 上陈述的所有结果都适用于它。

**定义。** $\mathbb C^n$ 上的**冯·诺依曼代数**是满足 $\mathcal M=\mathcal M''$ 的
$*$-代数 $\mathcal M$。它的**中心**是
$Z(\mathcal M)=\mathcal M\cap\mathcal M'$。**因子**是满足
$Z(\mathcal M)=\mathbb C1$ 的冯·诺依曼代数。

在无穷维情形，$\mathcal M=\mathcal M''$ 是一个拓扑条件（在弱算子拓扑下闭）。在有限维情形
它是自动成立的：

**引理 2.1.** 设 $\mathcal M$ 是 $\mathbb C^n$ 上的 $*$-代数，$K\subseteq\mathbb C^n$
是满足 $\mathcal MK\subseteq K$ 的子空间。则到 $K$ 上的正交投影 $P$ 属于
$\mathcal M'$。

*证明。* 对 $A\in\mathcal M$，有 $AK\subseteq K$ 且 $A^*K\subseteq K$（因为 $A^*\in\mathcal M$）。
第二个包含关系给出 $AK^\perp\subseteq K^\perp$：对 $y\in K^\perp$，$x\in K$，
$\langle x,Ay\rangle=\langle A^*x,y\rangle=0$。因此 $A$ 把 $K$ 映入 $K$，把 $K^\perp$
映入 $K^\perp$，这意味着 $AP=PAP=PA$。$\square$

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

所以在有限维情形，"冯·诺依曼代数"与"含单位元的矩阵 $*$-代数"是同一回事。尽管如此，
双交换子定理仍然提供了正确的视角：在物理上，可观测量代数常常*以交换子的形式给出*——
即与某个对称性可交换的算子全体（2.4 节）。

## 2.2 中心

**引理 2.3.** 设 $\mathcal Z$ 是 $\mathbb C^n$ 上的交换 $*$-代数。则存在
非零、两两正交且满足 $\sum_cP_c=1$ 的投影 $P_1,\dots,P_C\in\mathcal Z$，
使得 $\mathcal Z=\operatorname{span}\{P_1,\dots,P_C\}$。它们是 $\mathcal Z$ 的
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
$P_c$。极小性与唯一性：投影 $\sum_c a_cP_c$ 满足 $a_c\in\{0,1\}$，所以
$\mathcal Z$ 的非零投影恰好是若干个（非空集合的）$P_c$ 之和。$\square$

将此引理应用于 $\mathcal Z=Z(\mathcal M)$，投影 $P_c$ 称为 $\mathcal M$ 的**极小中心
投影**。每个 $A\in\mathcal M$ 可分解为
$A=\sum_cAP_c$，其中 $AP_c=P_cAP_c$，所以
$\mathcal M=\bigoplus_c\mathcal MP_c$（作用在
$E_c=P_c\mathbb C^n$ 上的代数的正交直和）。

**引理 2.4.** 每个 $\mathcal M_c=\mathcal MP_c$，作为 $E_c$ 上以 $P_c$ 为单位元的代数，都是
因子。

*证明。* 设 $y\in\mathcal M_c$ 与 $\mathcal M_c$ 可交换。对 $A\in\mathcal M$，
$Ay=AP_cy=(AP_c)y=y(AP_c)=yP_cA=yA$，这里用到 $P_c$ 是中心元以及 $y=yP_c=P_cy$。
所以 $y\in Z(\mathcal M)$，$y=\sum_d\chi_d(y)P_d$，从而 $y=yP_c=\chi_c(y)P_c$。$\square$

## 2.3 结构定理

**引理 2.5（因子）。** 设 $\mathcal N$ 是空间 $E$ 上的因子。则存在整数
$k,m\ge1$ 和酉映射 $V:\mathbb C^k\otimes\mathbb C^m\to E$，使得
$V^*\mathcal NV=M_k(\mathbb C)\otimes1_m$。此外 $V^*\mathcal N'V=1_k\otimes M_m(\mathbb C)$。

*证明。* *第 1 步：极小投影。* $\mathcal N$ 含有非零投影（例如 $1$）。
设 $e\in\mathcal N$ 是秩最小的非零投影。我们断言
$e\mathcal Ne=\mathbb Ce$。设 $x=x^*\in e\mathcal Ne$。$x$ 对应于非零特征值 $\lambda$ 的每个
谱投影 $q$ 都可写成 $q=p(x)$，其中 $p$ 是满足
$p(\lambda)=1$ 且对 $\mu\in(\operatorname{spec}x\cup\{0\})\setminus\{\lambda\}$ 有 $p(\mu)=0$
的 Lagrange 多项式（见定理 1.5(a) 的证明）。由于 $p(0)=0$ 且 $x=ex=xe$，有 $q\in\mathcal N$
且 $q=eqe$，所以 $q$ 的值域包含在 $e$ 的值域中。由极小性，$q=e$ 或 $q=0$；因此 $x$
至多有一个非零特征值，且其特征子空间就是 $e$ 的值域，所以
$x=\lambda e$。$e\mathcal Ne$ 的一般元素为 $x_1+ix_2$，其中
$x_{1,2}=x_{1,2}^*\in e\mathcal Ne$。

*第 2 步：$e$ 生成 $E$。* 令 $F=\operatorname{span}\mathcal Ne E$。它在 $\mathcal N$ 下
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
最大（这些值域是维数都为 $\operatorname{rank}e$ 的两两正交子空间，
所以 $k\le\dim E/\operatorname{rank}e$）。我们断言 $\sum_ie_i=1$。否则 $g=1-\sum_ie_i\neq0$ 是
$\mathcal N$ 中的投影。由第 2 步，对某个 $A\in\mathcal N$ 有 $gAe\neq0$（若对所有 $A$ 都有
$gAe=0$，则对所有 $A,\xi$ 有 $gAe\xi=0$，于是 $g$ 在 $F=E$ 上为零）。令 $y=gAe$。则 $y^*y\in e\mathcal Ne$ 是正的且非零，
所以 $y^*y=\lambda e$，$\lambda>0$，而 $v=y/\sqrt\lambda$ 是一个值域含于 $gE$、
与所有 $e_i$ 正交的 $e$-等距——这与 $k$ 的最大性矛盾。

对 $i\neq j$，$v_i^*v_j=v_i^*e_ie_jv_j=0$。令 $E_{ij}=v_iv_j^*$。则
$E_{ij}E_{kl}=v_iv_j^*v_kv_l^*=\delta_{jk}v_iev_l^*=\delta_{jk}E_{il}$，
$E_{ij}^*=E_{ji}$，且 $\sum_iE_{ii}=1$。对 $A\in\mathcal N$，
$A=\sum_{ij}e_iAe_j=\sum_{ij}v_i(v_i^*Av_j)v_j^*$，而 $v_i^*Av_j\in e\mathcal Ne=\mathbb Ce$，
所以 $A=\sum_{ij}a_{ij}E_{ij}$，$a_{ij}\in\mathbb C$。于是 $\mathcal N=\operatorname{span}\{E_{ij}\}$。

*第 4 步：酉映射。* 令 $m=\operatorname{rank}e$，$f_1,\dots,f_m$ 为 $eE$ 的一组规范正交基。
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

**定理 2.6（有限维冯·诺依曼代数的结构）。** 设
$\mathcal M$ 是 $\mathbb C^n$ 上的冯·诺依曼代数。则存在 $C\ge1$、满足 $\sum_ck_cm_c=n$ 的整数
$k_c,m_c\ge1$，以及酉映射
$V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n$，使得
$$
V^*\mathcal MV=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},\qquad
V^*\mathcal M'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),
$$
并且 $V^*Z(\mathcal M)V=\operatorname{span}\{P_1,\dots,P_C\}$，其中 $P_c$ 是到第 $c$ 个
直和项上的投影。数对 $(k_c,m_c)$ 在不计顺序的意义下唯一：$P_c$ 是唯一的（引理 2.3），
$k_c^2=\dim\mathcal MP_c$，且 $m_c=\operatorname{rank}P_c/k_c$。特别地，作为 $*$-代数有 $\mathcal M\cong\bigoplus_cM_{k_c}(\mathbb C)$，并且
$\dim\mathcal M=\sum_ck_c^2$，$\dim\mathcal M'=\sum_cm_c^2$。

*证明。* 由引理 2.3 及其后的注记，$\mathcal M=\bigoplus_c\mathcal M_c$，其中
$\mathcal M_c$ 作用在 $E_c$ 上，且每个 $\mathcal M_c$ 都是因子（引理 2.4）。引理 2.5
给出酉映射 $V_c$，使得 $V_c^*\mathcal M_cV_c=M_{k_c}\otimes1_{m_c}$；令
$V=\bigoplus_cV_c$。关于交换子：$X\in\mathcal M'$ 与每个 $P_c$ 可交换，所以
$X=\sum_cP_cXP_c$，且 $P_cXP_c$ 与 $\mathcal M_c$ 可交换；反之，任何这样的和都与
$\mathcal M$ 可交换。引理 2.5 确定了 $E_c$ 上的 $\mathcal M_c'$。$\square$

数 $k_c$ 描述代数本身。**重数** $m_c$ 描述代数如何嵌入 $M_n(\mathbb C)$：它们在代数结构中
不可见，但决定了 $\mathbb C^n$ 的迹在 $\mathcal M$ 上的限制，
$$
\operatorname{Tr}_{\mathbb C^n}(A)=\sum_cm_c\operatorname{Tr}(A_c),\qquad A=V\Big(\bigoplus_cA_c\otimes1_{m_c}\Big)V^*.
$$
这就是本库使用一族迹的原因（第 3 章）。

> **物理意义。** $P_c$ 是**超选择扇区**或**对称扇区**：$\mathcal M$ 中没有可观测量连接不同的
> 扇区，因此"处于哪个扇区"的测量与其他任何测量都是相容的。中心是系统的*经典*部分。
> 重数 $m_c$ 是一种简并度：块 $A_c$ 的每个特征值在 $\mathbb C^n$ 中出现的重数是它在 $A_c$
> 中重数的 $m_c$ 倍。

## 2.4 例子：对称群及其交换子

设有限群 $G$ 通过酉表示 $\pi$ 作用在 $\mathbb C^n$ 上。与该对称性可交换的算子构成
$\pi(G)'$，这是一个冯·诺依曼代数；它的交换子 $\pi(G)''$ 是诸 $\pi(g)$ 的线性张成（对该
张成应用定理 2.2，它是一个 $*$-代数：包含 $\pi(e)=1$，对乘积 $\pi(g)\pi(h)=\pi(gh)$ 和伴随
$\pi(g)^*=\pi(g^{-1})$ 封闭）。

**命题 2.8.** 对 $\mathcal M=\pi(G)''$ 应用定理 2.6。则
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

表示分解为不可约表示的方式在等价意义下唯一——从而 $m_c$ 是任何分解中*唯一确定的*
重数——这一点将在第 4 章证明。

我们取 $G=S_3$，它置换 $(\mathbb C^2)^{\otimes3}$ 的三个张量因子，即三个
量子比特或三个自旋 $\tfrac12$。表示论（一个特征标计算，第 4 章）
预言：对称子空间
（总自旋 $\tfrac32$，维数 $4$）承载重数为 $4$ 的平凡表示；它的补空间（维数 $4$）承载重数为
$2$ 的二维标准表示；符号表示不出现。所以
$\pi(S_3)''\cong M_1\oplus M_2$（维数 $5$），$\pi(S_3)'\cong M_4\oplus M_2$
（维数 $20$）。我们用线性代数检验这些维数：集合 $S$ 的交换子是
线性映射 $X\mapsto(AX-XA)_{A\in S}$ 的零空间。

```python
import itertools
import torch

torch.manual_seed(0)
dt = torch.complex128

def permutation_operator(perm, n_qubits=3):
    """Unitary on (C^2)^{(x)3} that permutes the tensor factors."""
    dim = 2 ** n_qubits
    P = torch.zeros(dim, dim, dtype=dt)
    for bits in itertools.product([0, 1], repeat=n_qubits):
        src = int(''.join(map(str, bits)), 2)
        dst = int(''.join(str(bits[perm[i]]) for i in range(n_qubits)), 2)
        P[dst, src] = 1
    return P

group = [permutation_operator(p) for p in itertools.permutations(range(3))]

def commutant_basis(ops, tol=1e-10):
    """Orthonormal basis (Hilbert-Schmidt) of {X : X A = A X for all A in ops}."""
    d = ops[0].shape[0]
    eye = torch.eye(d, dtype=dt)
    # vec(AX - XA) = (A (x) 1 - 1 (x) A^T) vec(X) for row-major vec
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    _, s, Vh = torch.linalg.svd(L)
    null = Vh[(s > tol).sum():].conj().resolve_conj()   # rows span the null space
    return [v.reshape(d, d) for v in null]

def span_dim(ops, tol=1e-10):
    return int((torch.linalg.svdvals(torch.stack([A.reshape(-1) for A in ops])) > tol).sum())

comm = commutant_basis(group)
print(len(comm), span_dim(group), len(commutant_basis(comm)))
assert len(comm) == 20                                   # pi(S3)'  = M_4 (+) M_2
assert span_dim(group) == 5                              # pi(S3)'' = M_1 (+) M_2
assert len(commutant_basis(comm)) == 5                   # double commutant theorem
```

重数在谱中是可见的。对群取平均，
$X\mapsto\frac1{|G|}\sum_g\pi(g)X\pi(g)^*$，把任意算子映入交换子
（习题 3）。因此 $\pi(S_3)'$ 的一个随机自伴元素（几乎必然）在对称子空间上有 $4$ 个单重
特征值，在其补空间上有 $2$ 个二重简并特征值；
代码检验了这些简并对：
它们来自 $M_2\otimes1_2$ 中的因子 $1_2$。

```python
def group_average(X):
    return sum(g @ X @ g.conj().T for g in group) / len(group)

H0 = torch.randn(8, 8, dtype=dt)
Hsym = group_average(H0 + H0.conj().T)                  # a random symmetric Hamiltonian
assert all(torch.allclose(g @ Hsym, Hsym @ g) for g in group)

P_sym = sum(group) / len(group)                         # projection onto the symmetric subspace
w, Q = torch.linalg.eigh(P_sym)
Q_rest = Q[:, w.real < 0.5]                             # orthonormal basis of the complement (dim 4)
e_rest = torch.linalg.eigvalsh(Q_rest.conj().T @ Hsym @ Q_rest)
print(e_rest)
assert Q_rest.shape[1] == 4
assert torch.allclose(e_rest[0], e_rest[1]) and torch.allclose(e_rest[2], e_rest[3])
```

## 2.5 库中的代数

本库在 $\bigoplus_c\mathbb C^{k_c}$ 上存储定理 2.6 中的抽象代数 $\bigoplus_cM_{k_c}(\mathbb C)$，
不带重数：一个算子就是它的块 $A_c$ 组成的列表，各块填充到共同的大小。上例中的交换子就是
`TypeIAlgebra([4, 2], [4, 2])`。

```python
from torch_vn_algebra import TypeIAlgebra, channels

alg = TypeIAlgebra([4, 2], [4, 2], precision='double')    # M_4 (+) M_2
print(alg.k_factors, alg.C)

A = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
print(A.matrix.shape)              # (1, 2, 4, 4): block 1 occupies the top-left 2 x 2 corner
assert A.matrix[0, 1, 2:, :].abs().max() == 0 and A.matrix[0, 1, :, 2:].abs().max() == 0
```

中心由元素 $\sum_cz_cP_c$ 组成；`alg.central` 构造这些元素，而它们与一切元素可交换：

```python
def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

P0, P1 = alg.central([1.0, 0.0]), alg.central([0.0, 1.0])   # minimal central projections
Z = alg.central([2.0, -1.0])
assert P0.is_projection and P1.is_projection
assert close(P0 + P1, alg.identity())
assert close(Z @ A, A @ Z)
assert close(A @ P0 + A @ P1, A)                         # A = sum_c A P_c
```

## 2.6 条件期望

从 $\mathcal M$ 到冯·诺依曼子代数
$\mathcal N\subseteq\mathcal M$ 上的**条件期望**是一个线性映射 $E:\mathcal M\to\mathcal N$，
它是正的（$A\ge0\Rightarrow E(A)\ge0$），对 $B\in\mathcal N$ 满足 $E(B)=B$，并且是
$\mathcal N$-双模映射：对 $B_1,B_2\in\mathcal N$ 有 $E(B_1AB_2)=B_1E(A)B_2$。

**命题 2.7.** 对 $\mathcal M=\bigoplus_cM_{k_c}$ 和 $\mathcal N=Z(\mathcal M)$，
$$
E(A)=\sum_c\frac{\operatorname{Tr}A_c}{k_c}\,P_c
$$
是一个条件期望。它保持每个形如
$\sum_cw_c\operatorname{Tr}A_c$ 的迹。

*证明。* 线性是显然的；$A\ge0$ 蕴含 $\operatorname{Tr}A_c\ge0$，所以 $E(A)\ge0$。对
$B=\sum_cb_cP_c$，$E(B)=\sum_cb_ck_cP_c/k_c=B$。对 $B_i=\sum_cb^{(i)}_cP_c$，
$\operatorname{Tr}(B_1AB_2)_c=b^{(1)}_cb^{(2)}_c\operatorname{Tr}A_c$，由此得双模
性质。最后，对每个 $c$ 有 $\operatorname{Tr}(E(A)_c)=\operatorname{Tr}A_c$。$\square$

在本库中，这个映射是信道 `channels.center_expectation(alg)`。把它作用于算子，得到其
经典（按扇区平均的）部分：

```python
E = channels.center_expectation(alg)
EA = E(A)
assert close(E(EA), EA)                                  # idempotent
assert close(E(Z @ A @ P1), Z @ EA @ P1)                 # bimodule property
A2 = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
assert close(EA @ A2, A2 @ EA)                           # E(A) commutes with a fresh element
for c in range(alg.C):                                   # Tr E(A)_c = Tr A_c
    assert abs(EA.matrix[0, c].diagonal().sum() - A.matrix[0, c].diagonal().sum()) < 1e-10
```

2.4 节中的群平均是 $M_n(\mathbb C)$ 中到交换子
$\pi(G)'$ 上的条件期望，并且保持 $\operatorname{Tr}$（习题 3）。
（条件期望一般不唯一：到 $\mathbb C1$ 上，每个态
$\varphi$ 都给出一个条件期望 $A\mapsto\varphi(A)1$。）两者都是同一思想的实例：
投影到"看不见"某种结构的算子上。

## 2.7 可与之比较的工具

2.4 节的例子可以用任何线性代数软件包检验；
表示论软件（GAP、SageMath）能以符号方式给出分解
$\pi\cong\bigoplus\rho_\lambda\otimes1_{m_\lambda}$。物理中的严格对角化代码（QuSpin、QuTiP
的对称性模块）使用同样的块结构，逐个扇区地对角化哈密顿量。本库以这一结构为
出发点：一旦块已知，后续所有计算——函数、迹、态、信道、动力学——都停留在代数内部。

## 习题

1. 证明 $M_k(\mathbb C)$ 的中心是 $\mathbb C1$（利用矩阵单位
   $|i\rangle\langle j|$），并证明 $Z(\bigoplus_cM_{k_c})$ 的维数为 $C$。
2. 设 $S$ 是 $\mathbb C^3$ 上的循环移位，$Se_i=e_{i+1\bmod3}$。确定
   $\{S\}'$ 与 $\{S\}''$，并用数值方法检验其维数。
3. 对有限群的酉表示 $\pi$，证明
   $E(X)=\frac1{|G|}\sum_g\pi(g)X\pi(g)^*$ 是从 $M_n(\mathbb C)$ 到 $\pi(G)'$ 上的
   条件期望，并且保持 $\operatorname{Tr}$。

<details>
<summary>解答</summary>

1. 若 $X$ 与所有 $|i\rangle\langle j|$ 可交换，则
   $X|i\rangle\langle j|=|i\rangle\langle j|X$；把两边作用于 $e_j$ 并比较，
   得对所有 $i,j$ 有 $Xe_i=X_{jj}e_i$，所以 $X=X_{11}1$。在直和中，元素
   $\bigoplus X_c$ 是中心元当且仅当每个 $X_c$ 都是，所以中心为
   $\{\bigoplus_cz_c1_{k_c}\}$，维数为 $C$。

2. $S$ 是酉矩阵，特征值为 $1,\omega,\omega^2$，$\omega=e^{2\pi i/3}$，均为单重。
   与 $S$ 可交换的算子保持 $S$ 的特征子空间，所以在特征基下它是
   对角的：$\{S\}'\cong\mathbb C^3$，即循环矩阵全体。于是 $\{S\}''$ 也
   $\cong\mathbb C^3$，因为在 $S$ 的特征基下，与所有对角矩阵可交换的矩阵是对角的（对对角元
   互不相同的 $D$，比较 $XD$ 与 $DX$ 的元素）；它由 $1,S,S^2$ 张成（Vandermonde：特征值互不相同）。

```python
S = torch.roll(torch.eye(3, dtype=dt), 1, dims=0)
assert len(commutant_basis([S])) == 3
assert len(commutant_basis(commutant_basis([S]))) == 3
assert span_dim([torch.eye(3, dtype=dt), S, S @ S]) == 3
```

3. $E$ 是线性且正的（它是映射 $X\mapsto UXU^*$ 之和）。对 $h\in G$，
   $\pi(h)E(X)\pi(h)^*=\frac1{|G|}\sum_g\pi(hg)X\pi(hg)^*=E(X)$，所以 $E(X)\in\pi(G)'$。若
   $X\in\pi(G)'$，则每一项都等于 $X$，所以 $E(X)=X$。对 $B_1,B_2\in\pi(G)'$，
   $\pi(g)B_1XB_2\pi(g)^*=B_1\pi(g)X\pi(g)^*B_2$，由此得双模性质。
   最后，对每个 $g$ 有 $\operatorname{Tr}(\pi(g)X\pi(g)^*)=\operatorname{Tr}X$，所以
   $\operatorname{Tr}E(X)=\operatorname{Tr}X$。

</details>

**下一章：** [第 3 章](03_traces_and_states.md)对
$\bigoplus_cM_{k_c}$ 上的迹进行分类，引入态与密度，并证明后文要用到的熵不等式与
迹不等式。
