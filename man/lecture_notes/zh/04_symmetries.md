[English](../04_symmetries.md) · [Español](../es/04_symmetries.md) · [Français](../fr/04_symmetries.md) · [Deutsch](../de/04_symmetries.md) · **中文** · [日本語](../ja/04_symmetries.md) · [Русский](../ru/04_symmetries.md)

# 4. 对称性与表示

第 2 章表明，一个对称群产生一个冯·诺依曼代数，即对称算子的线性张成，以及它的交换子，即对称的
可观测量（命题 2.7），并手工分解了一个例子。本章给出其背后的表示论，把它变成本库可以作用于任意
矩阵列表的算法，并把它应用于自旋系统的两种对称性：$U(1)$（粒子数或磁化强度）和 $SU(2)$（自旋的转动）。

**你将学到**完全可约性与 Schur 引理，分解为不可约表示的方式如何决定交换子与生成的代数，特征标如何
计算重数，为什么交换子的一个随机元素能数值地找到分解，什么是荷扇区，不可约自旋表示以及 $L$ 个
自旋 $\frac12$ 如何分解，为什么 $SU(2)$ 不变哈密顿量的块比 $U(1)$ 不变哈密顿量的块小得多，以及复共轭如何
把由实矩阵生成的代数的扇区分为实型、复型与四元数型。本章的数值检验见
[教程 4](../../tutorial/zh/04_symmetries.md)。

在本章中，$S'$ 是 $S\subseteq M_n(\mathbb C)$ 的交换子（第 2 章）。若 $S^*=S$，则 $S''$ 是**由** $S$
**生成的冯·诺依曼代数**：它是包含 $S$ 的 $*$-代数，并且每个 $*$-代数 $\mathcal B\supseteq S$ 都满足
$\mathcal B=\mathcal B''\supseteq S''$（定理 2.2）。$\mathcal A'_h=\{K\in\mathcal A':K=K^*\}$ 是一个实
向量空间。$\mathbb C^k\otimes\mathbb C^m$ 采用第 7 章的 Kronecker 顺序（$e_a\otimes f_j$ 是第 $am+j$ 个基向量，
从 $0$ 开始编号），因此 $1_k\otimes x$ 是含 $k$ 个 $x$ 拷贝的分块对角矩阵。$|G|$ 是有限群的阶，并取 $\hbar=1$。

<a id="s-4-1"></a>

## 4.1 表示与 Schur 引理

**定义。** 设 $X$ 是带有对合 $x\mapsto x^*$（$x^{**}=x$）的集合。$X$ 在 $\mathbb C^k$ 上的一个
**表示**是满足 $\rho(x^*)=\rho(x)^*$ 的映射 $\rho:X\to M_k(\mathbb C)$。若对所有 $x$ 有
$\rho(x)K\subseteq K$，则称子空间 $K$ 是**不变的**；若 $k\ge1$ 且 $0,\mathbb C^k$ 是仅有的不变子空间，
则称 $\rho$ 是**不可约的**。从 $\rho$ 到 $\mathbb C^{k'}$ 上的表示 $\rho'$ 的**缠结算子**是对所有 $x$ 满足
$T\rho(x)=\rho'(x)T$ 的 $k'\times k$ 矩阵 $T$；它们构成一个向量空间 $\mathrm{Hom}(\rho,\rho')$，并且
$\mathrm{Hom}(\rho,\rho)=\rho(X)'$。若存在可逆的缠结算子，则称 $\rho$ 与 $\rho'$ **等价**，记作 $\rho\cong\rho'$。

会出现三种情形：群 $G$，取 $g^*=g^{-1}$，此时要求 $\rho$ 是同态 $G\to U(k)$（**酉表示**；于是
$\rho(g^{-1})=\rho(g)^{-1}=\rho(g)^*$）；带有其伴随运算的 $*$-代数，此时 $\rho$ 是保单位的 $*$-同态；以及
满足 $x^*=x$ 的生成元标签，例如 $\{x,y,z\}$，用厄米矩阵组来表示，比如自旋矩阵（4.5 节）。若 $K$ 不变，
且 $Q\in M_{k\times d}$ 的列构成 $K$ 的一组规范正交基，则 $\rho(x)Q=QQ^*\rho(x)Q$，因此
$\rho_K(x)=Q^*\rho(x)Q$ 仍是一个表示（若 $\rho$ 是同态，则它也是同态），称为 $\rho$ 在 $K$ 上的**限制**。

<a id="r-4-1"></a>

**引理 4.1（完全可约性）。** 设 $\rho$ 是 $X$ 在 $\mathbb C^n$ 上的表示。
(a) 若 $K$ 不变，则 $K^\perp$ 也不变。(b) $\mathbb C^n$ 是若干不变子空间的正交直和，$\rho$ 在每个
这样的子空间上都不可约。

*证明。* (a) 对 $y\in K^\perp$ 和 $z\in K$，$\langle z,\rho(x)y\rangle=\langle\rho(x^*)z,y\rangle=0$，因为
$\rho(x^*)z\in K$。(b) 对 $n$ 归纳。维数最小的非零不变子空间 $K_1$ 是不可约的；由 (a)，$K_1^\perp$
是不变的，而限制到 $K_1^\perp$ 上的表示的不变子空间都是 $\rho$ 的不变子空间，所以归纳假设适用于它。$\square$

对有限群，酉性不需要任何代价：若诸 $\rho(g)$ 仅仅可逆，则 $P=\sum_g\rho(g)^*\rho(g)$ 满足
$P\ge\rho(e)^*\rho(e)=1$，所以 $P>0$，并且 $\rho(h)^*P\rho(h)=\sum_g\rho(gh)^*\rho(gh)=P$（作替换 $g\to gh^{-1}$）；
因此 $(P^{1/2}\rho(h)P^{-1/2})^*(P^{1/2}\rho(h)P^{-1/2})=P^{-1/2}\rho(h)^*P\rho(h)P^{-1/2}=1$，从而
$h\mapsto P^{1/2}\rho(h)P^{-1/2}$ 是与 $\rho$ 等价的酉表示（Weyl 技巧）。

<a id="r-4-2"></a>

**引理 4.2（Schur 引理）。** 设 $\rho,\rho'$ 是 $X$ 在 $\mathbb C^k,\mathbb C^{k'}$ 上的不可约表示，
$T\in\mathrm{Hom}(\rho,\rho')$。
(a) $T=0$ 或 $T$ 可逆。(b) $\mathrm{Hom}(\rho,\rho)=\mathbb C1$。
(c) $T^*\in\mathrm{Hom}(\rho',\rho)$，并且若 $T\neq0$，则 $T=\lambda U$，其中 $\lambda>0$，$U$ 为酉矩阵。
因此，当 $\rho\cong\rho'$ 时 $\dim\mathrm{Hom}(\rho,\rho')$ 为 $1$，否则为 $0$；并且等价的不可约表示是
酉等价的。

*证明。* (a) 和 (b) 已在命题 2.7 之后对群证明过，而证明只用到 $\ker T$ 在 $\rho$ 下不变
（若 $Tv=0$，则 $T\rho(x)v=\rho'(x)Tv=0$）以及 $T$ 的像在 $\rho'$ 下不变：所以
$\ker T\in\{0,\mathbb C^k\}$ 且 $\mathrm{ran}\,T\in\{0,\mathbb C^{k'}\}$；对 (b)，取 $\lambda$ 为 $T$ 的一个
特征值，$T-\lambda1$ 是不可逆的缠结算子，因此为 $0$。(c) 对 $T\rho(x^*)=\rho'(x^*)T$ 取伴随得
$\rho(x)T^*=T^*\rho'(x)$。于是 $T^*T\in\mathrm{Hom}(\rho,\rho)$，所以若 $T\neq0$，则 $T^*T=\lambda^21$，
$\lambda>0$（$T^*T\ge0$ 非零），而由 (a)，$U=T/\lambda$ 是方阵且 $U^*U=1$。最后，若
$T_1,T_2\in\mathrm{Hom}(\rho,\rho')$ 非零，则 $T_1^{-1}T_2\in\mathrm{Hom}(\rho,\rho)=\mathbb C1$。$\square$

若所有 $\rho(x)$ 彼此可交换（阿贝尔群），则不可约的 $\rho$ 是一维的：每个 $\rho(x)$ 都属于
$\mathrm{Hom}(\rho,\rho)=\mathbb C1$，所以每个子空间都是不变的。

<a id="s-4-2"></a>

## 4.2 同型分解与交换子

<a id="r-4-3"></a>

**定理 4.3（同型分解）。** 设 $\rho$ 是 $X$ 在 $\mathbb C^n$ 上的表示。存在 $X$ 在 $\mathbb C^{k_c}$ 上
两两不等价的不可约表示 $\rho_1,\dots,\rho_C$、满足 $\sum_ck_cm_c=n$ 的整数 $m_c\ge1$，以及酉算子
$V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n$，使得对所有 $x$ 有
$V^*\rho(x)V=\bigoplus_c\rho_c(x)\otimes1_{m_c}$。对每个这样的分解，
$$
V^*\rho(X)'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),\qquad
V^*\rho(X)''V=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},
$$
并且 $m_c=\dim\mathrm{Hom}(\rho_c,\rho)$；$X$ 的不可约表示 $\sigma$ 满足
$\mathrm{Hom}(\sigma,\rho)\neq0$ 当且仅当恰有一个 $c$ 使 $\sigma\cong\rho_c$。

*证明。* *存在性。* 引理 4.1 给出两两正交、和为 $\mathbb C^n$ 的不变子空间 $W_1,\dots,W_M$，其限制
$\sigma_i(x)=Q_i^*\rho(x)Q_i$ 不可约。把诸 $W_i$ 按限制的等价性分成若干类，并在每一类 $c$ 中选取一个参照
表示 $\rho_c=\sigma_{i_0}$。对该类中其他每个 $i$，引理 4.2(c) 给出酉矩阵 $U_i$ 使
$\sigma_i(x)=U_i\rho_c(x)U_i^*$；把 $Q_i$ 换成 $Q_iU_i$ 就有 $Q_i^*\rho(x)Q_i=\rho_c(x)$。把类 $c$ 的成员编号为
$j=0,\dots,m_c-1$，把它们的基记作 $Q_{c,j}$，并在第 $c$ 个直和项上令 $V(e_a\otimes f_j)=Q_{c,j}e_a$。由于
诸 $W_i$ 不变且相互正交，$Q_{c,j}^*\rho(x)Q_{d,l}=\delta_{cd}\delta_{jl}\rho_c(x)$，这正是所断言的形式。

*交换子。* 现在设 $V$ 是任意一个如定理所述的分解，并令 $Q_{c,j}$ 为 $n\times k_c$ 等距算子
$e_a\mapsto V(e_a\otimes f_j)$（在第 $c$ 个直和项中）；对刚才构造的 $V$，它们就是上面的 $Q_{c,j}$。恒等式
$V^*\rho(x)V=\bigoplus_c\rho_c(x)\otimes1_{m_c}$ 说的是 $\rho(x)Q_{d,l}=Q_{d,l}\rho_d(x)$，而把它的伴随应用于 $x^*$ 得
$Q_{c,j}^*\rho(x)=\rho_c(x)Q_{c,j}^*$。由于 $V$ 是酉的，$Y\in M_n$ 与所有 $\rho(x)$ 可交换，当且仅当每个块
$Y_{cj,dl}=Q_{c,j}^*YQ_{d,l}$ 满足 $Y_{cj,dl}\rho_d(x)=\rho_c(x)Y_{cj,dl}$，即属于
$\mathrm{Hom}(\rho_d,\rho_c)$，而它在 $c\neq d$ 时为 $0$，在 $c=d$ 时为 $\mathbb C1$（引理 4.2）。所以
$Y\in\rho(X)'$ 当且仅当 $Y_{cj,dl}=\delta_{cd}(y_c)_{jl}1_{k_c}$，即 $V^*YV=\bigoplus_c1_{k_c}\otimes y_c$。

*双交换子。* 若 $Z$ 与 $\rho(X)'$ 可交换，则它与到各直和项上的投影可交换（取 $y_c=1$，并对
$d\neq c$ 取 $y_d=0$），所以 $V^*ZV=\bigoplus_cZ_c$，且 $Z_c$ 与 $1_{k_c}\otimes M_{m_c}$ 可交换；由引理 2.5
第 5 步（交换两个张量因子），$Z_c\in M_{k_c}\otimes1_{m_c}$。反之，这些算子都与 $\rho(X)'$ 可交换。

*重数。* 对 $T\in\mathrm{Hom}(\sigma,\rho)$，其中 $\sigma$ 不可约，同样的论证表明
$Q_{c,j}^*T\in\mathrm{Hom}(\sigma,\rho_c)$，当 $\sigma\cong\rho_c$ 时它等于 $\mathbb Ct_{c}$（$t_c$ 可逆），否则为
$0$；并且 $T=\sum_{c,j}Q_{c,j}Q_{c,j}^*T$。由于诸 $\rho_c$ 两两不等价，$\sigma$ 至多与其中一个等价，设为
$\rho_c$（否则 $T=0$），于是 $T=\sum_ja_jQ_{c,j}t_c$；反之，每个 $Q_{c,j}t_c$ 都属于
$\mathrm{Hom}(\sigma,\rho)$，而这 $m_c$ 个映射的像相互正交，所以 $\dim\mathrm{Hom}(\sigma,\rho)=m_c$。$\square$

*库中：* 如定理 4.3 所述的分解由 `decompose.from_group` 和 `decompose.from_generators` 计算——[教程，4.3 节](../../tutorial/zh/04_symmetries.md#s-4-3)。

对群而言，$\rho(G)''$ 是 $\rho(G)$ 的线性张成（2.4 节），而定理 4.3 就是反过来读的命题 2.7：在那里由
结构定理得出不可约表示，在这里由不可约表示得出结构；两种情形下都有 $\dim\rho(X)'=\sum_cm_c^2$ 和
$\dim\rho(X)''=\sum_ck_c^2$。$\rho(X)''$ 的极小中心投影 $P_c=\sum_jQ_{c,j}Q_{c,j}^*$ 的像是 $\rho_c$ 的
**同型分量**。
由定理 4.3，每个厄米的 $H\in\rho(X)'$ 都具有形式 $V(\bigoplus_c1_{k_c}\otimes h_c)V^*$，其中 $h_c\in M_{m_c}$
厄米，而 $1_{k_c}\otimes h_c$ 是含 $k_c$ 个 $h_c$ 拷贝的分块对角矩阵：$H$ 的谱是诸 $m_c\times m_c$ 块 $h_c$ 的谱的
并，$h_c$ 的每个特征值重复 $k_c$ 次。

> **物理意义**（物理解释，而非本课程的定理）。若 $\rho(X)$ 是一个对称性，$H\in\rho(X)'$ 是哈密顿量，则刚才
> 描述的 $k_c$ 重简并能级就是该对称性的**多重态**，而需要对角化的块的大小是重数 $m_c$，而不是多重态的大小。

对有限群，重数不需要分解任何东西就能算出。

<a id="r-4-4"></a>

**定理 4.4（Schur 正交关系）。** 设 $\rho,\rho'$ 是有限群 $G$ 在 $\mathbb C^k,\mathbb C^{k'}$ 上的不可约酉
表示。则
$$
\frac1{|G|}\sum_{g\in G}\rho(g)_{ij}\,\overline{\rho'(g)_{i'j'}}=
\begin{cases}\delta_{ii'}\delta_{jj'}/k&\text{if }\rho=\rho',\\0&\text{if }\rho\not\cong\rho'.\end{cases}
$$

*证明。* 对 $k\times k'$ 矩阵 $Y$，令 $T=\frac1{|G|}\sum_g\rho(g)Y\rho'(g)^*$。对 $h\in G$，作替换 $g\to hg$
并利用 $\rho'(hg)^*=\rho'(g)^*\rho'(h)^{-1}$，得 $\rho(h)T=T\rho'(h)$，所以 $T\in\mathrm{Hom}(\rho',\rho)$。
若 $\rho\not\cong\rho'$，则 $T=0$；若 $\rho=\rho'$，则 $T=\lambda1$，且 $k\lambda=\mathrm{Tr}\,T=\mathrm{Tr}\,Y$
（引理 4.2）。对 $Y=|j\rangle\langle j'|$，矩阵元 $T_{ii'}$ 就是左边，而 $\mathrm{Tr}\,Y=\delta_{jj'}$。$\square$

表示的**特征标**是 $\chi_\rho(g)=\mathrm{Tr}\,\rho(g)$；等价的表示有相同的特征标，因为
$\mathrm{Tr}(T\rho(g)T^{-1})=\mathrm{Tr}\,\rho(g)$。记
$\langle\varphi,\psi\rangle_G=\frac1{|G|}\sum_g\overline{\varphi(g)}\psi(g)$。

<a id="r-4-5"></a>

**推论 4.5（特征标）。** 设 $G$ 有限。
(a) 对不可约酉表示 $\rho,\rho'$：当 $\rho\cong\rho'$ 时 $\langle\chi_\rho,\chi_{\rho'}\rangle_G$ 为 $1$，否则为 $0$。
(b) 对按定理 4.3 分解的酉表示 $\pi$ 和不可约表示 $\rho$，$\rho$ 在 $\pi$ 中的重数为
$$
m_\rho=\frac1{|G|}\sum_{g\in G}\overline{\chi_\rho(g)}\,\chi_\pi(g)
$$
（若 $\rho\cong\rho_c$ 则 $m_\rho=m_c$；若 $\rho$ 不与任何 $\rho_c$ 等价则 $m_\rho=0$）。
(c) $\dim\pi(G)'=\sum_cm_c^2=\frac1{|G|}\sum_{g\in G}|\chi_\pi(g)|^2$。特别地，$\pi$ 不可约当且仅当这个数为 $1$。

*证明。* (a) 把 $\rho'$ 换成一个等价表示（特征标不变）。对定理 4.4 在 $i=j$、$i'=j'$ 上求和并取共轭：
若 $\rho=\rho'$，结果为 $\sum_{i,i'}\delta_{ii'}/k=1$；若 $\rho\not\cong\rho'$，结果为 $0$。(b) 迹是酉不变的，且
$\mathrm{Tr}(a\otimes1_m)=m\mathrm{Tr}\,a$，所以 $\chi_\pi=\sum_cm_c\chi_{\rho_c}$；应用 (a)。(c) 由同样的展开和 (a)，
$\langle\chi_\pi,\chi_\pi\rangle_G=\sum_cm_c^2$，由定理 4.3 它等于 $\dim\pi(G)'$。$\square$

*计算：* 对置换三个量子比特的 $S_3$ 验证推论 4.5——[教程，4.2 节](../../tutorial/zh/04_symmetries.md#s-4-2)；对 $S_3$ 的标准表示验证定理 4.4——[教程，4.3 节](../../tutorial/zh/04_symmetries.md#s-4-3)。

(c) 的第二个证明正是本库对交换子进行采样的方式。**twirl**
$P(Y)=\frac1{|G|}\sum_g\pi(g)Y\pi(g)^*$ 是幂等的，其像为 $\pi(G)'$（第 2 章习题 3），并且关于
$\langle Y,Z\rangle_{HS}=\mathrm{Tr}(Y^*Z)$ 是自伴的（作替换 $g\to g^{-1}$），因此它是到 $\pi(G)'$ 上的正交
投影，而 $\dim\pi(G)'$ 就是它的迹。映射 $Y\mapsto AYB$ 的迹为 $\mathrm{Tr}\,A\mathrm{Tr}\,B$（其对角元为
$A_{ii}B_{jj}$），由此得到 $\frac1{|G|}\sum_g|\chi_\pi(g)|^2$。

*库中：* twirl 即 `decompose.twirl`——[教程，4.3 节](../../tutorial/zh/04_symmetries.md#s-4-3)；$\dim\pi(G)'$ 即 `decompose.commutant_dimension`——[教程，4.2 节](../../tutorial/zh/04_symmetries.md#s-4-2)。

我们用特征标重做 2.4 节的例子。置换矩阵的迹等于它固定的基向量的个数，而 $\pi(\sigma)$ 固定
$|b_1b_2b_3\rangle$ 当且仅当 $b$ 在 $\sigma$ 的每个轮换上为常数（引理 B.1(c)），所以
$\chi_\pi(\sigma)=2^{\#\text{cycles}(\sigma)}$。手头有 $S_3$ 的三个酉表示。设 $P_\sigma e_i=e_{\sigma(i)}$ 是 $\mathbb C^3$ 上的
置换表示（到置换矩阵中的同态）。**平凡**表示是 $\sigma\mapsto1$。**符号**表示是
$\mathrm{sgn}\,\sigma=\det P_\sigma$，它是到 $\{\pm1\}\subseteq U(1)$ 中的同态，因为行列式是可乘的，而置换矩阵
的行列式为 $\pm1$。直线 $\mathbb C(1,1,1)$ 是不变的，因此 $E=(1,1,1)^\perp$ 也是不变的（引理 4.1(a)）；**标准**
表示是 $P$ 在 $E$ 上的限制。在适配于 $\mathbb C(1,1,1)\oplus E$ 的规范正交基下，$P_\sigma=1\oplus\mathrm{std}(\sigma)$，
所以 $\chi_{\rm std}(\sigma)=\mathrm{Tr}\,P_\sigma-1$，即不动点个数减一。$S_3$ 的六个元素是恒等元、三个对换和
两个 $3$-轮换，它们分别有 $3,2,1$ 个轮换和 $3,1,0$ 个不动点；对换交换 $1_3$ 的两列，所以其符号为 $-1$，而
$3$-轮换是两个对换的乘积（$(1\,2\,3)=(1\,3)(1\,2)$），所以其符号为 $+1$。由此得到下表

| | 恒等元（$1$） | 对换（$3$） | $3$-轮换（$2$） |
|---|---|---|---|
| $\chi_{\rm triv}$ | $1$ | $1$ | $1$ |
| $\mathrm{sgn}$ | $1$ | $-1$ | $1$ |
| $\chi_{\rm std}$ | $2$ | $0$ | $-1$ |
| $\chi_\pi$ | $8$ | $4$ | $2$ |

（所以 $\mathrm{sgn}\,\sigma=(-1)^{3-\#\text{cycles}(\sigma)}$；所有特征标都是实的）。逐个共轭类计算
$\langle\varphi,\psi\rangle_G$：
$\langle\chi_{\rm triv},\chi_{\rm triv}\rangle_G=\frac{1+3+2}6=1$，$\langle\mathrm{sgn},\mathrm{sgn}\,\rangle_G=\frac{1+3+2}6=1$，
$\langle\chi_{\rm std},\chi_{\rm std}\rangle_G=\frac{4+0+2}6=1$，所以由推论 4.5(c)，这三个表示都是不可约的，并且
由于它们的特征标不同，它们两两不等价。推论 4.5(b) 给出它们在 $\pi$ 中的重数：
$$
m_{\rm triv}=\tfrac{8+3\cdot4+2\cdot2}6=4,\qquad m_{\rm sgn}=\tfrac{8-3\cdot4+2\cdot2}6=0,\qquad
m_{\rm std}=\tfrac{2\cdot8+0-2\cdot2}6=2 .
$$
所以平凡表示和标准表示分别等价于 $\pi$ 的成分 $\rho_c$（定理 4.3），其 $m_c=4$ 和 $m_c=2$；它们占据
$4\cdot1+2\cdot2=8=\sum_ck_cm_c$ 维，所以没有其他成分：$k=(1,2)$，$m=(4,2)$，并且
$\dim\pi(S_3)'=4^2+2^2=\frac16(64+3\cdot16+2\cdot4)=20$（推论 4.5(c)）。

这就是 2.4 节的分解，现在由每个表示的一张六个数的表得到。

*计算：* 由 $\pi$ 重新计算这张表——[教程，4.2 节](../../tutorial/zh/04_symmetries.md#s-4-2)。

<a id="s-4-3"></a>

## 4.3 从生成元到代数：一个随机化算法

定理 4.3 描述了分解，但并没有计算它：它的证明需要不可约的不变子空间。Murota、Kanno、Kojima 和
Kojima（*Japan J. Indust. Appl. Math.* 27 (2010) 125–160；另见 Maehara and Murota, *ibid.* 263–293）
观察到，*交换子的一个随机元素*的特征子空间能找到它们。设 $\mathcal A$ 是由给定矩阵 $S_1,\dots,S_J$ 及其
伴随生成的冯·诺依曼代数，其分解为 $V^*\mathcal AV=\bigoplus_cM_{k_c}\otimes1_{m_c}$，
$V^*\mathcal A'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$（定理 2.6）；只有分析中用到 $V$。

<a id="r-4-6"></a>

**引理 4.6.** 设 $K\in\mathcal A'_h$，$V^*KV=\bigoplus_c1_{k_c}\otimes\kappa_c$。$K$ 的每个特征子空间都是
$\mathcal A$-不变的。若所有 $\kappa_c$ 合在一起的 $N=\sum_cm_c$ 个特征值两两不同，则 $K$ 的特征子空间恰为
$N$ 个子空间 $V(\mathbb C^{k_c}\otimes e)$（在第 $c$ 个直和项中，$e$ 为 $\kappa_c$ 的特征向量），并且
$\mathcal A$ 在其中每一个上都不可约地作用，如同 $M_{k_c}$。

*证明。* 若 $Kv=\lambda v$ 且 $a\in\mathcal A$，则 $Kav=aKv=\lambda av$。写
$\kappa_c=\sum_j\mu_{cj}|e_{cj}\rangle\langle e_{cj}|$，其中 $e_{cj}$ 规范正交（定理 1.3）。则
$V^*KV=\sum_{c,j}\mu_{cj}R_{cj}$，其中 $R_{cj}$ 是到 $\mathbb C^{k_c}\otimes e_{cj}$ 上的投影；诸 $R_{cj}$ 相互
正交、和为 $1$，且诸 $\mu_{cj}$ 互不相同，所以它们就是谱投影。在 $\mathbb C^{k_c}\otimes e_{cj}$ 上，元素
$V^*aV=\bigoplus_da_d\otimes1$ 作用为 $a_c$，而诸 $a_c$ 取遍 $M_{k_c}$，它除 $0$ 与 $\mathbb C^{k_c}$ 外没有
不变子空间。$\square$

<a id="r-4-7"></a>

**命题 4.7（一般性）。** 设 $K$ 和 $Y$ 分别是 $\mathcal A'_h$ 和 $\mathcal A'$ 中相互独立的随机元素，其分布在这些
实向量空间上关于 Lebesgue 测度有密度（引理 A.3）。则以概率一：(a) 引理 4.6 的假设成立；(b) 对 $K$ 的任意两个
像位于同一直和项中的不同谱投影 $\Pi_i,\Pi_j$，有 $\Pi_iY\Pi_j\neq0$。

*证明。* (a) 厄米矩阵 $\kappa=\bigoplus_c\kappa_c\in M_N$ 实线性地依赖于 $K$。$1,\kappa,\dots,\kappa^{N-1}$ 关于
$\langle\cdot,\cdot\rangle_{HS}$ 的 Gram 行列式 $\delta(K)=\det[\mathrm{Tr}\,\kappa^{i+j}]_{i,j=0}^{N-1}$ 是 $K$ 的坐标的实
多项式，它非零当且仅当这些幂线性无关。由于对 $\kappa=U\mathrm{diag}(\lambda)U^*$ 有
$p(\kappa)=U\mathrm{diag}(p(\lambda_l))U^*$，一个次数 $<N$ 的非零多项式零化 $\kappa$ 当且仅当它在
$\mathrm{spec}\,\kappa$ 上为零，而这当且仅当 $\kappa$ 的不同特征值少于 $N$ 个时才可能（对这些特征值取
$\prod_\mu(t-\mu)$）。所以 $\delta(K)\neq0$ 就是引理 4.6 的假设，并且 $\delta\not\equiv0$：它在
$V(\bigoplus_c1\otimes\mathrm{diag}(\text{互不相同的实数}))V^*\in\mathcal A'_h$ 处非零。因此由推论 A.7(b)，
$\mathbb P(\delta(K)=0)=0$。

(b) 集合 $O=\{\delta\neq0\}\subseteq\mathcal A'_h$ 是开集（$\delta$ 连续）。对 $K\in O$，$K$ 恰有 $N$ 个互不相同的特征值
$\mu_1(K)<\dots<\mu_N(K)$（即 $\kappa$ 的特征值，重数为诸 $k_c$），所以它按特征值递增编号的谱投影
$\Pi_1(K),\dots,\Pi_N(K)$ 在 $O$ 上连续（命题 A.13）。固定 $l\neq l'$，令 $\mathcal N_{ll'}$ 为满足下列条件的
$(K,Y)\in O\times\mathcal A'$ 的集合：对某个 $c$ 有 $P_c\Pi_l(K)=\Pi_l(K)$ 且 $P_c\Pi_{l'}(K)=\Pi_{l'}(K)$（$P_c$ 为中心投影），
并且 $\Pi_l(K)Y\Pi_{l'}(K)=0$；它由 Borel 集 $O\times\mathcal A'$ 上的连续方程定义，因此是 Borel 集。固定一个使两个像
都位于直和项 $c$ 中的 $K\in O$。由引理 4.6，$V^*\Pi_lV$ 和 $V^*\Pi_{l'}V$ 分别投影到 $\mathbb C^{k_c}\otimes e$ 和
$\mathbb C^{k_c}\otimes e'$ 上，其中 $e,e'$ 是 $\kappa_c$ 的规范正交特征向量，而对 $V^*YV=\bigoplus_d1\otimes y_d$，
$\Pi_lY\Pi_{l'}=0$ 当且仅当 $\langle e,y_ce'\rangle=0$。映射 $Y\mapsto\langle e,y_ce'\rangle$ 是实线性映射
$\mathcal A'\to\mathbb C$ 且非零（它在 $V(0\oplus\dots\oplus1\otimes|e\rangle\langle e'|\oplus\dots\oplus0)V^*$ 处的值为 $1$），
所以它的核是真子空间，是零测集（推论 A.7(c)）。所以每个截面 $\{Y:(K,Y)\in\mathcal N_{ll'}\}$ 都是零测集（若两个像
位于不同的直和项中，它是空集），而由于 $Y$ 与 $K$ 独立，由推论 A.7(d) 得 $\mathbb P((K,Y)\in\mathcal N_{ll'})=0$。
(b) 中的事件只在 $\{\delta(K)=0\}\cup\bigcup_{l\neq l'}\{(K,Y)\in\mathcal N_{ll'}\}$ 上不成立，而这是有限个概率为 $0$ 的
事件之并。$\square$

<a id="r-4-8"></a>

**引理 4.8（压缩；对齐各拷贝）。** (a) 设 $W$ 是一个 $\mathcal A$-不变子空间，其规范正交基为
$Q\in M_{n\times d}$，$\sigma(a)=Q^*aQ$ 为限制。则 $\sigma(\mathcal A)$ 是 $*$-代数，
$\{Q^*YQ:Y\in\mathcal A'\}=\sigma(\mathcal A)'$，且 $\{Q^*KQ:K\in\mathcal A'_h\}=\sigma(\mathcal A)'_h$。$\mathcal A$ 在 $W$ 上
不可约地作用，当且仅当对每个 $Y\in\mathcal A'$ 有 $Q^*YQ\in\mathbb C1$。若不然，且 $X$ 是 $\mathcal A'$ 或
$\mathcal A'_h$ 中具有密度的随机元素，则以概率一 $Q^*XQ\notin\mathbb C1$。

(b) 设 $W_i,W_j$ 是 $\mathcal A$ 在其上不可约地作用的不变子空间，其规范正交基为
$Q_i,Q_j$，限制为 $\sigma_i(a)=Q_i^*aQ_i$，并设 $Y\in\mathcal A'$。则 $T=Q_i^*YQ_j\in\mathrm{Hom}(\sigma_j,\sigma_i)$。
若 $T\neq0$，则 $T=\lambda U$，其中 $\lambda>0$，$U$ 为酉矩阵，并且 $Q_j'=Q_jU^*$ 对所有 $a\in\mathcal A$ 满足
$Q_j'^*aQ_j'=\sigma_i(a)$。若 $W_i,W_j$ 位于不同的直和项中（$W_i\subseteq P_c\mathbb C^n$，
$W_j\subseteq P_d\mathbb C^n$，$c\neq d$），则 $T=0$。

*证明。* (a) 与定理 4.3 的证明一样（$W$ 在 $a$ 和 $a^*$ 下不变），对 $a\in\mathcal A$ 有 $aQ=Q\sigma(a)$ 且
$Q^*a=\sigma(a)Q^*$。因此 $\sigma(ab)=Q^*aQ\sigma(b)=\sigma(a)\sigma(b)$，$\sigma(a^*)=\sigma(a)^*$，
$\sigma(1)=Q^*Q=1_d$，所以 $\sigma$ 是保单位的 $*$-同态，$\sigma(\mathcal A)$ 是 $*$-代数；并且对 $Y\in\mathcal A'$，
$Q^*YQ\sigma(a)=Q^*YaQ=Q^*aYQ=\sigma(a)Q^*YQ$。反之，对 $T\in\sigma(\mathcal A)'$，令 $Y=QTQ^*$；则
$aY=Q\sigma(a)TQ^*=QT\sigma(a)Q^*=Ya$，所以 $Y\in\mathcal A'$，$Q^*YQ=T$，并且若 $T$ 厄米则 $Y$ 也厄米。若 $\sigma$
不可约，则 $\sigma(\mathcal A)'=\mathrm{Hom}(\sigma,\sigma)=\mathbb C1$（引理 4.2(b)）。若 $E\subseteq\mathbb C^d$ 是不同于
$0,\mathbb C^d$ 的 $\sigma$-不变子空间，则 $E^\perp$ 也是（引理 4.1(a)），而到 $E$ 上的正交投影是
$\sigma(\mathcal A)'_h$ 中的一个非数乘元素。在这种情形下，$\mathbb C1\cap\sigma(\mathcal A)'$ 和
$\mathbb R1\cap\sigma(\mathcal A)'_h$ 是真实子空间，它们在满射实线性映射 $X\mapsto Q^*XQ$ 下的原像是 $\mathcal A'$ 和
$\mathcal A'_h$ 的真子空间，而这些都是零测集（推论 A.7(c)）。

(b) 与 (a) 一样，$aQ_j=Q_j\sigma_j(a)$ 且 $Q_i^*a=\sigma_i(a)Q_i^*$，所以
$T\sigma_j(a)=Q_i^*YaQ_j=Q_i^*aYQ_j=\sigma_i(a)T$。引理 4.2 给出 $T=\lambda U$ 和
$U\sigma_j(a)U^*=\sigma_i(a)$，后者就是 $Q_j'^*aQ_j'$。最后，$P_c\in\mathcal A\cap\mathcal A'$，$Q_i=P_cQ_i$，
$P_cQ_j=0$，所以 $T=Q_i^*P_cYQ_j=Q_i^*YP_cQ_j=0$。$\square$

*计算：* 对 $\pi(S_3)'$ 的一个随机元素验证引理 4.6 和引理 4.8(b)——[教程，4.3 节](../../tutorial/zh/04_symmetries.md#s-4-3)。

于是算法（`decompose.wedderburn`，前端为 `from_group` 和 `from_generators`）如下：

1. *对交换子采样。* 对有限群的全部元素，取 $Y=P(Z)$ 和 $K=P(Z'+Z'^*)$，其中 $P$ 为 twirl，$Z,Z'$ 为相互
   独立的 Ginibre 矩阵。$P$ 是 $M_n$ 到 $\mathcal A'$ 上的线性投影，且 $P(Z)^*=P(Z^*)$，所以 $Z\mapsto P(Z)$ 和
   $Z'\mapsto P(Z'+Z'^*)$ 是到 $\mathcal A'$ 和到 $\mathcal A'_h$ 上的实线性满射（对 $H\in\mathcal A'_h$，
   $P(\frac H2+\frac H2)=H$）；Ginibre 矩阵是 $M_n$ 中的非退化中心化高斯向量（命题 A.16），所以 $Y$ 和 $K$ 是
   $\mathcal A'$ 和 $\mathcal A'_h$ 中的非退化中心化高斯向量，从而有密度（定理 A.15）。对生成元，
   $\mathcal A'=\{X:[S_j,X]=[S_j^*,X]=0\}$ 是 $M_n$ 上一个线性映射的零空间（本库使用厄米部分
   $\frac12(S_j+S_j^*)$、$\frac1{2i}(S_j-S_j^*)$，它们与 $S_j,S_j^*$ 张成同一个复空间，因而有相同的交换子）；取它的一组基
   $(B_r)_{r\le D}$，令 $Y=\sum_rc_rB_r$，$K=Y'+Y'^*$，其中 $Y'$ 是 $Y$ 的独立拷贝，诸 $c_r$ 是相互独立的标准复
   高斯变量：它们是 $\mathbb C^D$ 中非退化高斯向量在到 $\mathcal A'$ 上的实线性双射 $c\mapsto\sum_rc_rB_r$ 和到
   $\mathcal A'_h$ 上的实线性满射 $c\mapsto X+X^*$（$X=\sum_rc_rB_r$）下的像，因此同样有密度。每个随机元素都是
   新的抽样，与之前的抽样独立；本库把它们重新缩放为单位范数，这不改变任何特征子空间，也不改变任何块是否为零。
2. *分裂。* $K$ 的特征子空间是不可约的（引理 4.6，命题 4.7(a)）；它们的维数就是诸 $k_c$。
3. *分组与对齐。* 维数相同的特征子空间之间的块 $Q_i^*YQ_j$ 恰好在同一直和项内非零（命题 4.7(b)，
   引理 4.8(b)）；关系"由一个非零块相连"的各等价类给出扇区，类的大小给出 $m_c$，而极分解因子把各拷贝对齐到
   每类的一个参照拷贝（引理 4.8(b)；对齐到一个已对齐拷贝的拷贝也就对齐到参照拷贝）。
4. *组装* Kronecker 顺序下的 $V$，使得 $V^*aV=\bigoplus_ca_c\otimes1_{m_c}$。

所以，在精确算术中，若 $K$ 的特征子空间被精确地分开，则算法以概率一返回一个如定理 4.3 所述的分解：由引理
4.6，直和项 $c$ 包含 $m_c$ 个维数为 $k_c$ 的特征子空间，而对齐之后，这一构造正是定理 4.3 证明中的存在性部分。
在浮点运算中，本库把相差小于某个容差的特征值视为相等，所以计算出的"特征子空间" $W$ 可能是 $K$ 的若干特征
子空间之和；它仍是 $\mathcal A$-不变的（引理 4.6），本库把一个新的随机元素 $X\in\mathcal A'_h$ 压缩到它上面。由
引理 4.8(a)，只要 $\mathcal A$ 在 $W$ 上不可约地作用，$Q^*XQ$ 就是 $1$ 的倍数，否则以概率一不是；此时它的特征子
空间在 $\sigma(\mathcal A)$ 下不变（把引理 4.6 应用于 $\sigma(\mathcal A)$），把 $W$ 分裂为更小的不变子空间，本库重复
这一步，直到每个压缩都是数乘。无论浮点运算中发生什么，结果最后都会得到验证：

<a id="r-4-9"></a>

**命题 4.9（证书）。** 设 $\mathcal A$ 由 $S_1,\dots,S_J$ 及其伴随生成，
$\mathcal B=V(\bigoplus_cM_{k_c}\otimes1_{m_c})V^*$，其中 $V$ 为酉矩阵。若所有 $S_j\in\mathcal B$ 且
$\dim\mathcal A'=\sum_cm_c^2$，则 $\mathcal A=\mathcal B$。

*证明。* $\mathcal B$ 是包含所有 $S_j,S_j^*$ 的 $*$-代数，所以 $\mathcal A\subseteq\mathcal B$，
$\mathcal B'\subseteq\mathcal A'$。把定理 4.3 应用于 $\mathcal B$ 在 $\mathbb C^n$ 上的表示（分块映射
$b\mapsto b_c$ 是不可约的，并且互不等价，因为第 $c$ 个中心投影被其中一个映为 $1$，被其余的映为 $0$），得
$\dim\mathcal B'=\sum_cm_c^2=\dim\mathcal A'$。所以 $\mathcal B'=\mathcal A'$，
$\mathcal B=\mathcal B''=\mathcal A''=\mathcal A$。$\square$

*库中：* 该算法 `decompose.wedderburn`（前端为 `from_group`、`from_generators`）、命题 4.9 的判据 `Decomposition.verify()` 及其代价——[教程，4.3 节](../../tutorial/zh/04_symmetries.md#s-4-3)。

<a id="s-4-4"></a>

## 4.4 $U(1)$ 与荷扇区

<a id="r-4-10"></a>

**命题 4.10.** 设 $Q$ 是 $\mathbb C^n$ 上具有整数特征值 $q$（**荷**）的厄米矩阵，其谱投影为 $P_q$，
$d_q=\mathrm{rank}\,P_q$，并令 $U_\theta=e^{i\theta Q}$；记 $U(1)'=\{U_\theta:\theta\in\mathbb R\}'$，并用
$U(1)''$ 表示它的交换子。

1. $\theta\mapsto U_\theta$ 是 $U(1)=\mathbb R/2\pi\mathbb Z$ 的酉表示，并且 $U(1)$ 在 $\mathbb C^n$ 上的每个
   连续酉表示都具有这种形式。
2. 对 $X\in M_n$，下列条件等价：(i) 对所有 $\theta$ 有 $XU_\theta=U_\theta X$；(ii) $[X,Q]=0$；
   (iii) $X=\sum_qP_qXP_q$。因此 $U(1)'\cong\bigoplus_qM_{d_q}(\mathbb C)$（分块对角算子），且
   $U(1)''=\mathrm{span}\{P_q\}$：所有 $k_c=1$，而重数就是诸 $d_q$。
3. （*守恒荷。*）对厄米的 $H$，$e^{itH}Qe^{-itH}=Q$ 对所有 $t$ 成立当且仅当 $[H,Q]=0$，即当且仅当
   $H\in U(1)'$。此时 $\rho(t)=e^{-itH}\rho e^{itH}$ 的每个扇区权重 $\mathrm{Tr}(P_q\rho(t))$ 都是常数。

*证明。* (1) $U_\theta=\sum_qe^{iq\theta}P_q$ 是酉的，$U_\theta U_{\theta'}=U_{\theta+\theta'}$（定理 1.5(c)），且
$U_{2\pi}=1$。反之，对连续表示，$(U_\theta)_{\theta\ge0}$ 是矩阵的连续半群，所以由引理 8.5（在第 8 章证明；
证明只用到有限维）对 $\theta\ge0$ 有 $U_\theta=e^{\theta B}$，而对 $\theta<0$ 有
$U_\theta=U_{-\theta}^{-1}=e^{\theta B}$。在 $0$ 处对 $U_\theta^*U_\theta=1$ 求导（引理 8.1）得 $B^*+B=0$，所以
$B=iQ$，其中 $Q$ 厄米。指数级数 $e^{i\theta Q}$ 等于函数演算（命题 8.2(1)，取 $H=-Q$），所以
$U_{2\pi}=\sum_qe^{2\pi iq}P_q=1$ 迫使 $Q$ 的每个特征值 $q$ 满足 $e^{2\pi iq}=1$，$q\in\mathbb Z$。(2) (i)$\Rightarrow$(ii)：在
$\theta=0$ 处求导。(ii)$\Rightarrow$(iii)：$X$ 与 $P_q=f(Q)$ 可交换（定理 1.5(d)），所以当 $q\neq q'$ 时
$P_qXP_{q'}=XP_qP_{q'}=0$。(iii)$\Rightarrow$(i)：$XU_\theta=\sum_qe^{iq\theta}P_qXP_q=U_\theta X$。分块对角算子构成
$\bigoplus_qM_{d_q}\otimes1_1$（在适配于诸 $P_q$ 的像的规范正交基下），其交换子为
$\bigoplus_q1_{d_q}\otimes M_1=\mathrm{span}\{P_q\}$：与它们可交换的算子与每个 $P_q$ 可交换，所以是分块对角的，
而每个块与 $M_{d_q}$ 可交换，因而是数乘（引理 2.5 第 5 步，取 $m=1$）。
(3) $\frac d{dt}e^{itH}Qe^{-itH}=ie^{itH}[H,Q]e^{-itH}$（命题 8.2）对所有 $t$ 为零当且仅当 $[H,Q]=0$。此时
$e^{itH}$ 与 $P_q$ 可交换（定理 1.5(d)，用两次），并且
$\mathrm{Tr}(P_q\rho(t))=\mathrm{Tr}(e^{itH}P_qe^{-itH}\rho)=\mathrm{Tr}(P_q\rho)$。$\square$

因此 $U(1)$ 的不可约表示就是特征标 $\theta\mapsto e^{iq\theta}$，$q\in\mathbb Z$（与每个阿贝尔群一样，它们是
一维的），而荷 $q$ 的重数为 $d_q$。荷在张量积下相加：若 $u_k$ 和 $w_l$ 分别是 $Q_1$ 和 $Q_2$ 的规范正交
特征基，特征值为 $\alpha_k$ 和 $\beta_l$，则规范正交基 $u_k\otimes w_l$ 由 $Q_1\otimes1+1\otimes Q_2$ 的特征向量（特征值
$\alpha_k+\beta_l$）组成，同时也由 $e^{i\theta Q_1}\otimes e^{i\theta Q_2}$ 的特征向量（特征值 $e^{i\theta\alpha_k}e^{i\theta\beta_l}$）
组成，所以由函数演算的定义，
$$
e^{i\theta(Q_1\otimes1+1\otimes Q_2)}=e^{i\theta Q_1}\otimes e^{i\theta Q_2},
$$
它在 $P_{q_1}\otimes P_{q_2}$ 的像上以荷 $q_1+q_2$ 作用。这就是第 7 章融合张量积的融合规则
$(q_1,q_2)\mapsto q_1+q_2$。

对 $n$ 个量子比特，荷是**数算子** $\hat N=\sum_jn_j$，其中 $n_j=|1\rangle\langle1|$ 作用在量子比特 $j$ 上：
$\hat N|b\rangle=|b|\,|b\rangle$，$|b|$ 为汉明重量。它的扇区维数为 $\binom nw$，因此保持汉明重量守恒的算子构成
$$
U(1)'\cong\bigoplus_{w=0}^nM_{\binom nw}(\mathbb C),
$$
这就是 $n$ 个荷为 $0,1$ 的 $M_1\oplus M_1$ 拷贝的融合积（比特 $1$ = 自旋向上，
$\hat N=S^z+\frac n2$）。每个 XXZ 哈密顿量
$\sum_{\langle ij\rangle}J(s^x_is^x_j+s^y_is^y_j+\Delta s^z_is^z_j)$ 都保持 $\hat N$ 守恒：$s^z_is^z_j$ 在基 $|b\rangle$ 下是对角的，
而跃迁项 $s^x_is^x_j+s^y_is^y_j=\frac12(s^+_is^-_j+s^-_is^+_j)$（展开 $s^\pm=s^x\pm is^y$；交叉项相消）把 $|b\rangle$ 映为
$0$，或映为一个向上的自旋从格点 $i,j$ 中的一个移到另一个所得的构型，所以它把 $\hat N$ 的每个特征子空间映入自身。

*计算：* 对四个自旋验证命题 4.10：仅由 $\hat N$ 找出荷扇区，并逐个扇区对角化一个 XXZ 环——[教程，4.4 节](../../tutorial/zh/04_symmetries.md#s-4-4)。*库中：* $n$ 个量子比特的 $U(1)'$ 是 `SpinChain(n).algebra`，XXZ 哈密顿量是 `SpinChain.xxz`——[教程，4.4 节](../../tutorial/zh/04_symmetries.md#s-4-4)。

> **物理意义**（物理解释，而非本课程的定理）。$U(1)$ 对称性就是粒子数、荷或磁化强度的守恒。命题 4.10(3) 是这一情形下的 Noether
> 定理：荷守恒当且仅当哈密顿量是对称的，此时动力学永远不会离开一个扇区（定理 8.3）。若此外没有任何东西
> 能测量扇区之间的相干性，这些扇区就是超选择扇区（第 3 章）。

<a id="s-4-5"></a>

## 4.5 $SU(2)$：自旋

**定义。** $\mathbb C^d$ 上的一个**自旋表示**是满足
$$
[S_x,S_y]=iS_z,\qquad[S_y,S_z]=iS_x,\qquad[S_z,S_x]=iS_y,
$$
的厄米矩阵三元组 $S_x,S_y,S_z$，即标签 $\{x,y,z\}$ 满足这些关系的表示（4.1 节）。基本例子是
$\mathbb C^2$ 上的 $s_a=\frac12\sigma_a$（Pauli 矩阵）。若 $(S_a)$ 与 $(T_a)$ 是自旋表示，则
$S_a\otimes1+1\otimes T_a$ 也是（交叉项可交换）；因此 $L$ 个自旋 $\frac12$ 的**总自旋**
$S_a=\sum_is_a^{(i)}$（`decompose.total_spin(L)`）是一个自旋表示。转动 $e^{-i\theta\,\mathbf n\cdot\mathbf S}$
（$\theta\in\mathbb R$，$\mathbf n\in\mathbb R^3$ 为单位向量，函数演算）对所有 $\theta,\mathbf n$ 都与 $X$ 可交换，当且仅当
$X$ 与三个 $S_a$ 可交换（在 $\theta=0$ 处求导，命题 8.2(1)；反过来用定理 1.5(d)），所以"$SU(2)$ 作用"的交换子是
$\{S_x,S_y,S_z\}'$，并且由定理 2.2，转动与诸 $S_a$ 生成同一个代数（两个集合都在取伴随下封闭）。对 $L$ 个自旋
$\frac12$，转动就是 $V^{\otimes L}$，$V\in SU(2)$。事实上，$\mathbf n\cdot\mathbf S=\sum_ih^{(i)}$，其中
$h=\mathbf n\cdot\mathbf s$，而与 4.4 节中的荷一样（对 $L$ 归纳），$e^{-i\theta\,\mathbf n\cdot\mathbf S}=(e^{-i\theta h})^{\otimes L}$，
且 $e^{-i\theta h}\in SU(2)$（它的特征值为 $e^{\mp i\theta/2}$，因为 $h$ 无迹且 $h^2=\frac14$）。反之，$V\in SU(2)$ 是正规的，
其特征值的模为 $1$、乘积为 $1$，所以 $V=W\mathrm{diag}(e^{i\varphi},e^{-i\varphi})W^*$，其中 $W$ 为酉矩阵，并且对无迹
厄米矩阵 $A=W\mathrm{diag}(\varphi,-\varphi)W^*$ 有 $V=e^{iA}$。无迹厄米 $2\times2$ 矩阵是 Pauli 矩阵的实线性组合
$A=2\,\mathbf c\cdot\mathbf s$，所以 $V=e^{-i\theta\,\mathbf n\cdot\mathbf s}$，其中 $\theta=2|\mathbf c|$，
$\mathbf n=-\mathbf c/|\mathbf c|$（若 $\mathbf c=0$ 则 $\theta=0$）。定理 B.6(a),(c) 与推论 B.7（取 $d=2$）直接表明，
诸 $V^{\otimes L}$（$V\in SU(2)$）张成由诸 $S_a$ 生成的代数。

令 $S_\pm=S_x\pm iS_y$，于是 $S_-=S_+^*$，并令 $\mathbf S^2=S_x^2+S_y^2+S_z^2$。由这些关系得
$$
[S_z,S_\pm]=\pm S_\pm,\qquad[S_+,S_-]=2S_z,\qquad\mathbf S^2=S_-S_++S_z^2+S_z,
$$
并且 $\mathbf S^2$ 与每个 $S_a$ 可交换（例如
$[\mathbf S^2,S_z]=S_x[S_x,S_z]+[S_x,S_z]S_x+S_y[S_y,S_z]+[S_y,S_z]S_y=-i\{S_x,S_y\}+i\{S_y,S_x\}=0$）。

<a id="r-4-11"></a>

**定理 4.11（不可约自旋表示）。** 对每个 $s\in\{0,\frac12,1,\frac32,\dots\}$，`decompose.spin_matrices(s)` 的
矩阵——作用在以 $|m\rangle$（$m=s,s-1,\dots,-s$）为基的 $\mathbb C^{2s+1}$ 上，
$$
S_z|m\rangle=m|m\rangle,\qquad S_+|m\rangle=\sqrt{s(s+1)-m(m+1)}\;|m+1\rangle,\qquad S_-=S_+^*,
$$
构成一个不可约自旋表示（**自旋 $s$**），且 $\mathbf S^2=s(s+1)1$。每个不可约自旋表示都恰与其中一个酉等价。

*证明。* *自旋 $s$ 是自旋表示。* $[S_z,S_+]=S_+$ 是显然的。$S_+S_-|m\rangle=(s(s+1)-(m-1)m)|m\rangle$，
$S_-S_+|m\rangle=(s(s+1)-m(m+1))|m\rangle$，所以 $[S_+,S_-]=2S_z$；这两个关系及其伴随等价于三个定义关系。
$\mathbf S^2|m\rangle=(s(s+1)-m(m+1)+m^2+m)|m\rangle=s(s+1)|m\rangle$。*不可约性：* 不变子空间 $W\neq0$ 在
$S_z$ 下不变，而 $S_z$ 的特征值都是单重的，所以 $W$ 包含某个 $|m\rangle$；$S_\pm$ 的系数只在 $m=\pm s$ 处为零，
所以反复作用 $S_\pm$ 得到所有 $|m'\rangle$。

*分类。* 设 $(S_a)$ 在 $\mathbb C^d$ 上不可约，$\lambda$ 是 $S_z$ 的最大特征值，$v_0$ 是相应的单位特征向量。
$S_zS_+v_0=S_+(S_z+1)v_0=(\lambda+1)S_+v_0$，所以由最大性 $S_+v_0=0$。令 $v_j=S_-^jv_0$；则
$S_zv_j=(\lambda-j)v_j$，并且利用 $S_+S_-=S_-S_++2S_z$ 归纳可得
$$
S_+v_j=j(2\lambda-j+1)\,v_{j-1}
$$
（$S_+v_j=S_-S_+v_{j-1}+2S_zv_{j-1}=[(j-1)(2\lambda-j+2)+2(\lambda-j+1)]v_{j-1}$）。非零的 $v_j$ 是 $S_z$ 属于
不同特征值的特征向量，所以存在最后一个非零的 $v_r$。$v_0,\dots,v_r$ 的张成在 $S_z,S_\pm$ 下不变，因此等于
$\mathbb C^d$，且 $d=r+1$。由 $0=S_+v_{r+1}=(r+1)(2\lambda-r)v_r$ 得 $\lambda=r/2=:s$。最后，
$\|v_j\|^2=\langle v_{j-1},S_+v_j\rangle=j(2s-j+1)\|v_{j-1}\|^2$，所以单位向量 $|s-j\rangle=v_j/\|v_j\|$ 满足
$S_+|m\rangle=\sqrt{(s-m)(s+m+1)}\,|m+1\rangle$，其中 $m=s-j$，这正是上面的公式。不同的 $s$ 有不同的维数。$\square$

*库中：* 自旋 $s$ 即 `decompose.spin_matrices(s)`，在[教程，4.6 节](../../tutorial/zh/04_symmetries.md#s-4-6)中对照定理 4.11 检验；$L$ 个自旋 $\frac12$ 的总自旋即 `decompose.total_spin(L)`——[教程，4.5 节](../../tutorial/zh/04_symmetries.md#s-4-5)。

<a id="r-4-12"></a>

**推论 4.12.** 每个 $\mathbb C^n$ 上的自旋表示都酉等价于 $\bigoplus_s(\text{spin }s)\otimes1_{m_s}$，并且记
$d_M=\dim\ker(S_z-M)$，有
$$
m_s=d_s-d_{s+1}.
$$

*证明。* 定理 4.3 与定理 4.11。若 $|M|\le s'$ 且 $s'-M\in\mathbb Z$，则自旋 $s'$ 的每个拷贝恰含 $S_z$ 属于特征值
$M$ 的一个特征向量，否则不含；所以 $d_M=\sum_{s'\ge|M|,\,s'-M\in\mathbb Z}m_{s'}$，而两个相邻值之差就是 $m_s$。$\square$

<a id="r-4-13"></a>

**命题 4.13（$L$ 个自旋 $\frac12$）。** $(\mathbb C^2)^{\otimes L}$ 上的总自旋包含自旋 $S$ 当且仅当
$S\in\{\frac L2,\frac L2-1,\dots\}$，$S\ge0$，其重数为
$$
m_S=\binom L{\frac L2-S}-\binom L{\frac L2-S-1}\qquad\Big(\tbinom L{-1}=0\Big),
$$
所以由 $S_x,S_y,S_z$ 生成的代数是 $\bigoplus_SM_{2S+1}\otimes1_{m_S}$，其交换子是
$\bigoplus_S1_{2S+1}\otimes M_{m_S}$。

*证明。* $S_z=\sum_is_z^{(i)}$ 在乘积基下是对角的，在 $u$ 个自旋向上、$L-u$ 个自旋向下的乘积态上特征值为
$u-\frac L2$，所以若 $\frac L2+M\in\{0,\dots,L\}$，则
$d_M=\binom L{L/2+M}=\binom L{L/2-M}$，否则 $d_M=0$。应用推论 4.12 与定理 4.3。这些值是正的，因为当
$j<\frac{L-1}2$ 时 $\binom L{j+1}/\binom Lj=\frac{L-j}{j+1}>1$。$\square$

对 $L=2$：$\frac12\otimes\frac12=1\oplus0$（**Clebsch–Gordan**）：三重态 $S_-^j|{\uparrow\uparrow}\rangle$，$j=0,1,2$，
在两个自旋的交换 $F$ 下对称；单态 $\frac1{\sqrt2}(|{\uparrow\downarrow}\rangle-|{\downarrow\uparrow}\rangle)$ 反对称，
并被 $S_+$ 和 $S_z$ 零化。对 $L=4$：$m_2=1$，$m_1=4-1=3$，$m_0=6-4=2$，所以
$$
\{S_x,S_y,S_z\}''\cong M_5\oplus M_3\oplus M_1\ \text{ with multiplicities }1,3,2,\qquad
\{S_x,S_y,S_z\}'\cong M_1\oplus M_3\oplus M_2,
$$
维数分别为 $35$ 和 $14$——正是第 10 章中得到的数。

*计算：* 仅由三个生成元求 $L=4$ 时的分解，并验证命题 4.13 以及自旋 $S$ 上 $\mathbf S^2=S(S+1)$——[教程，4.5 节](../../tutorial/zh/04_symmetries.md#s-4-5)。

哪些算子与所有转动可交换？格点置换就是：用 $\pi(\sigma)$ 共轭把 $s^{(i)}_a$ 映为 $s^{(\sigma(i))}_a$（引理 B.1(b)），所以它
保持 $\sum_is^{(i)}_a$ 不变。反之，除此之外再没有别的：

<a id="r-4-14"></a>

**定理 4.14（Schur–Weyl 对偶）。** 设 $\pi$ 是 $S_L$ 置换 $(\mathbb C^2)^{\otimes L}$ 各张量因子的表示。
(a) $\mathrm{span}\,\pi(S_L)\subseteq\{S_x,S_y,S_z\}'$。(b) 等号成立；等价地（定理 2.2），由总自旋生成的
代数是 $\pi(S_L)'$。

*证明。* (a) 部分刚刚已证。(b) 部分就是推论 B.7，即定理 B.6 在 $d=2$ 时的情形，而定理 B.6 对所有 $d$ 和 $L$
证明了 $(\mathbb C^d)^{\otimes L}$ 与 $U(d)$ 的对角作用的对偶。$\square$

（延伸阅读：W. Fulton and J. Harris, *Representation Theory*, Springer 1991, Lecture 6；R. Goodman and
N. Wallach, *Symmetry, Representations, and Invariants*, Springer 2009, Section 9.1。）

*计算：* 对 $L\le5$ 验证定理 4.14，$\{S_a\}'$ 的维数 $2,5,14,42$（推论 B.10），以及 $\pi(S_4)''$——其数对 $(k_c,m_c)$ 正是总自旋的数对互换角色（推论 B.9）——[教程，4.5 节](../../tutorial/zh/04_symmetries.md#s-4-5)。

> **物理意义**（物理解释，而非本课程的定理）。对每个转动不变的哈密顿量，一个自旋 $S$ 多重态由 $2S+1$ 个能量相同的态组成，它们由
> $S^z=M$ 区分。两个自旋 $\frac12$ 构成一个三重态和一个单态；四个自旋构成一个五重态、三个三重态和两个
> 单态。两个单态的重数空间和三个三重态的重数空间就是第 10 章中无退相干的量子比特与无噪声的三能级系统（qutrit）。

<a id="s-4-6"></a>

## 4.6 自旋系统

$L$ 个自旋 $\frac12$ 上的 **Heisenberg 模型**是 $H=\sum_{i<j}J_{ij}\,\mathbf s^{(i)}\cdot\mathbf s^{(j)}$，其中耦合
$J_{ij}$ 为实数。

<a id="r-4-15"></a>

**命题 4.15.** (a)（*Dirac 恒等式*）$\mathbf s^{(i)}\cdot\mathbf s^{(j)}=\frac12F_{ij}-\frac14$，其中 $F_{ij}$
交换自旋 $i$ 与 $j$。因此每个 Heisenberg 哈密顿量都属于 $\mathrm{span}\,\pi(S_L)\subseteq\{S_a\}'$，并且在命题
4.13 的分解中 $V^*HV=\bigoplus_S1_{2S+1}\otimes H_S$，其中块 $H_S\in M_{m_S}$：$H_S$ 的每个特征值都是 $H$ 的一个
$(2S+1)$ 重简并能级。
(b) **最高权向量**空间 $\ker S_+\cap\ker(S_z-S)$ 的维数为 $m_S$，它在 $H$ 下不变，且 $H$ 在其上作用为
$H_S$（相差一个酉基变换）。

*证明。* (a) 在两个自旋上，$\mathbf s^{(1)}\cdot\mathbf s^{(2)}=\frac12(\mathbf S^2-\mathbf s^{(1)2}-\mathbf s^{(2)2})=\frac12\mathbf S^2-\frac34$
（定理 4.11，$s=\frac12$）。三重态向量 $S_-^j|{\uparrow\uparrow}\rangle$ 是 $|{\uparrow\uparrow}\rangle$、
$|{\downarrow\uparrow}\rangle+|{\uparrow\downarrow}\rangle$、$2|{\downarrow\downarrow}\rangle$，它们与单态一起构成
$\mathbb C^2\otimes\mathbb C^2$ 的一组正交基。在三重态上 $\mathbf S^2=2$ 且 $F=1$：$\mathbf S^2=S_-S_++S_z^2+S_z$ 在
$|{\uparrow\uparrow}\rangle$ 上为 $0+1+1=2$，而 $\mathbf S^2$ 和 $F$ 与 $S_-$ 可交换，$F$ 固定 $|{\uparrow\uparrow}\rangle$。在单态上
$\mathbf S^2=0$（它被 $S_+$ 和 $S_z$ 零化）且 $F=-1$。两边在这四个向量上一致。对 $L$ 个格点中的格点 $i,j$，与单位
算子作张量积即可。其余部分是定理 4.14(a) 与定理 4.3。(b) 我们证明更多的结论：对按推论 4.12 分解的任意自旋表示
和任意 $\mu\in\frac12\mathbb Z$，若 $\mu\ge0$，则 $\ker S_+\cap\ker(S_z-\mu)$ 是 $|\mu,\mu\rangle\otimes\mathbb C^{m_\mu}$（在推论
4.12 的酉算子下的像），若 $\mu<0$，则为 $0$。事实上，在该基下，满足 $S_z=\mu$ 的向量是对满足 $S'-\mu\in\mathbb Z$ 的
$S'\ge|\mu|$ 求和的 $\sum_{S'}|S',\mu\rangle\otimes w_{S'}$，而由定理 4.11，$S_+$ 把它映为
$\sum_{S'}\sqrt{(S'-\mu)(S'+\mu+1)}\,|S',\mu+1\rangle\otimes w_{S'}$，其中 $S'+\mu+1\ge1$；所以它为零当且仅当对所有
$S'\neq\mu$ 有 $w_{S'}=0$。对 Heisenberg 哈密顿量 $H$，$1\otimes H_S$ 在 $|S,S\rangle\otimes\mathbb C^{m_S}$ 上作用为
$H_S$。$\square$

*计算：* 在四个自旋的环和八个自旋的环上验证命题 4.15——[教程，4.6 节](../../tutorial/zh/04_symmetries.md#s-4-6)。

(b) 部分就是模块 `su2`（自旋链，第 9 章）不经数值分解而构造 $H_S$ 的方式：它逐个格点地添加，用
$j\otimes\frac12$ 的 Clebsch–Gordan 系数构造最高权向量的一组规范正交基。我们从最高权出发通过下降来推导这些
系数。由定理 4.11 和 $S_-=S_+^*$，自旋 $j$ 满足
$$
S_\pm|j,m\rangle=\sqrt{(j\mp m)(j\pm m+1)}\;|j,m\pm1\rangle ,
$$
并且对 $|m|>j$ 我们令 $|j,m\rangle=0$（这与该公式相容，其系数在 $m=\pm j$ 处为零）。

<a id="r-4-16"></a>

**定理 4.16（$j\otimes\frac12$ 的 Clebsch–Gordan 系数）。** 设 $j\in\{\frac12,1,\frac32,\dots\}$，令
$S_a=S^{(j)}_a\otimes1+1\otimes s_a$ 作用在 $\mathbb C^{2j+1}\otimes\mathbb C^2$ 上（定理 4.11 的自旋 $j$，基为 $|j,m\rangle$，
与一个自旋 $\frac12$ 耦合，$|{\uparrow}\rangle=|\frac12,\frac12\rangle$，$|{\downarrow}\rangle=|\frac12,-\frac12\rangle$），并令
$\alpha_M=\sqrt{(j+M+\frac12)/(2j+1)}$，$\beta_M=\sqrt{(j-M+\frac12)/(2j+1)}$。对 $J=j\pm\frac12$ 和 $M=J,J-1,\dots,-J$，令
$$
\begin{aligned}
|j+\tfrac12,M\rangle\!\rangle&=\alpha_M\,|j,M-\tfrac12\rangle\otimes|{\uparrow}\rangle+\beta_M\,|j,M+\tfrac12\rangle\otimes|{\downarrow}\rangle,\\
|j-\tfrac12,M\rangle\!\rangle&=-\beta_M\,|j,M-\tfrac12\rangle\otimes|{\uparrow}\rangle+\alpha_M\,|j,M+\tfrac12\rangle\otimes|{\downarrow}\rangle .
\end{aligned}
$$
(a) 这 $2(2j+1)$ 个向量构成 $\mathbb C^{2j+1}\otimes\mathbb C^2$ 的一组规范正交基。
(b) $|j+\frac12,j+\frac12\rangle\!\rangle=|j,j\rangle\otimes|{\uparrow}\rangle$，而在相差一个模为一的因子的意义下，
$|j-\frac12,j-\frac12\rangle\!\rangle$ 是满足 $S_zx=(j-\frac12)x$ 和 $S_+x=0$ 的唯一单位向量 $x$；该因子由
$|j,j\rangle\otimes|{\downarrow}\rangle$ 的系数为正来确定（Condon–Shortley 约定）。
(c) $S_z|J,M\rangle\!\rangle=M|J,M\rangle\!\rangle$，$S_-|J,M\rangle\!\rangle=\sqrt{(J+M)(J-M+1)}\,|J,M-1\rangle\!\rangle$：态
$|J,M\rangle\!\rangle$ 由 $|J,J\rangle\!\rangle$ 经下降并归一化得到。
(d) 因此对 $a=x,y,z$ 有 $S_a|J,M\rangle\!\rangle=\sum_{M'}(S^{(J)}_a)_{M'M}|J,M'\rangle\!\rangle$：酉算子
$|J,M\rangle\mapsto|J,M\rangle\!\rangle$ 把自旋 $(j+\frac12)\oplus$ 自旋 $(j-\frac12)$ 映为 $S_a$，所以
$j\otimes\frac12\cong(j+\frac12)\oplus(j-\frac12)$，而 Clebsch–Gordan 系数 $\langle j,m;\frac12,\sigma|J,M\rangle$
（$|J,M\rangle\!\rangle$ 在乘积基下的坐标）就是上面显示的 $\pm\alpha_M,\beta_M$。

*证明。* (a) 若上式中某个基向量为 $0$，则其系数也为 $0$：$|j,M-\frac12\rangle=0$ 仅当 $M=-j-\frac12$，此时
$\alpha_M=0$；$|j,M+\frac12\rangle=0$ 仅当 $M=j+\frac12$，此时 $\beta_M=0$。所以每个向量的范数平方为
$\alpha_M^2+\beta_M^2=1$。$M$ 不同的向量涉及不同的乘积基向量，而对相同的 $M$，内积为
$-\alpha_M\beta_M+\beta_M\alpha_M=0$。它们共有 $(2j+2)+2j=2(2j+1)$ 个。
(b) 第一个结论就是 $\alpha_{j+1/2}=1$，$\beta_{j+1/2}=0$。特征子空间 $S_z=j-\frac12$ 由
$|j,j-1\rangle\otimes|{\uparrow}\rangle$ 和 $|j,j\rangle\otimes|{\downarrow}\rangle$ 张成，而由 $s_+|{\uparrow}\rangle=0$，$s_+|{\downarrow}\rangle=|{\uparrow}\rangle$，
$$
S_+\big(p\,|j,j-1\rangle\otimes|{\uparrow}\rangle+q\,|j,j\rangle\otimes|{\downarrow}\rangle\big)=(p\sqrt{2j}+q)\,|j,j\rangle\otimes|{\uparrow}\rangle ,
$$
它为零当且仅当 $q=-\sqrt{2j}\,p$；结合 $|p|^2+|q|^2=1$ 和 $q>0$（因而 $p<0$），得 $p=-\beta_{j-1/2}=-1/\sqrt{2j+1}$，
$q=\alpha_{j-1/2}=\sqrt{2j/(2j+1)}$。(c) $S_z$：两项都有 $S_z=(M\mp\frac12)\pm\frac12=M$。对 $S_-$，令
$x_M=p\,|j,M-\frac12\rangle\otimes|{\uparrow}\rangle+q\,|j,M+\frac12\rangle\otimes|{\downarrow}\rangle$；由于 $s_-|{\uparrow}\rangle=|{\downarrow}\rangle$，$s_-|{\downarrow}\rangle=0$，
$$
S_-x_M=p\sqrt{(j+M-\tfrac12)(j-M+\tfrac32)}\;|j,M-\tfrac32\rangle\otimes|{\uparrow}\rangle
+\Big(p+q\sqrt{(j+M+\tfrac12)(j-M+\tfrac12)}\Big)|j,M-\tfrac12\rangle\otimes|{\downarrow}\rangle .
$$
我们利用 $\beta_M\sqrt{(j+M+\frac12)(j-M+\frac12)}=(j-M+\frac12)\alpha_M$ 和 $\alpha_M\sqrt{(j+M+\frac12)(j-M+\frac12)}=(j+M+\frac12)\beta_M$。
对 $J=j+\frac12$，$(p,q)=(\alpha_M,\beta_M)$，$c=\sqrt{(J+M)(J-M+1)}=\sqrt{(j+M+\frac12)(j-M+\frac32)}$：第一个系数为
$\sqrt{(j+M+\frac12)(j+M-\frac12)(j-M+\frac32)/(2j+1)}=c\,\alpha_{M-1}$，第二个为 $\alpha_M(j-M+\frac32)=c\,\beta_{M-1}$。
对 $J=j-\frac12$，$(p,q)=(-\beta_M,\alpha_M)$，$c=\sqrt{(j+M-\frac12)(j-M+\frac12)}$：第一个系数为
$-\sqrt{(j-M+\frac12)(j+M-\frac12)(j-M+\frac32)/(2j+1)}=-c\,\beta_{M-1}$，第二个为 $\beta_M(j+M-\frac12)=c\,\alpha_{M-1}$。
在两种情形下，对 $M>-J$ 都有 $S_-|J,M\rangle\!\rangle=c\,|J,M-1\rangle\!\rangle$。对 $M=-J$，两个系数都为零（若 $J=j+\frac12$，
则 $p=\alpha_M=0$；若 $J=j-\frac12$，则 $j+M-\frac12=0$），所以 $S_-|J,-J\rangle\!\rangle=0$。(d) 由 (a) 和 (c)，在规范正交基
$|J,M\rangle\!\rangle$ 下，$S_z$ 和 $S_-$ 的矩阵就是自旋 $(j+\frac12)\oplus$ 自旋 $(j-\frac12)$ 的矩阵（定理 4.11 与上面的
公式）；$S_+=S_-^*$ 的矩阵也是如此，因而 $S_x=\frac12(S_++S_-)$ 和 $S_y=\frac1{2i}(S_+-S_-)$ 的矩阵也是如此。$\square$

*计算：* 对 $j\le2$ 用稠密矩阵并对照 `su2.clebsch_gordan` 验证定理 4.16——[教程，4.6 节](../../tutorial/zh/04_symmetries.md#s-4-6)。

习题 2 把这一分解（但不包括系数）推广到 $j_1\otimes j_2$。对本库来说，只有定理 4.16 的两个最高权向量是重要的，
而它们在任意自旋表示中都适用，无论是否不可约：

<a id="r-4-17"></a>

**推论 4.17（逐个格点构造最高权向量）。** 设 $(T_a)$ 是 $E$ 上的自旋表示，令
$\mathcal H_\mu(E)=\ker T_+\cap\ker(T_z-\mu)$，其规范正交基为 $B_\mu$（若 $\mathcal H_\mu(E)=0$，则 $B_\mu=\emptyset$），并令
$S_a=T_a\otimes1+1\otimes s_a$ 作用在 $E\otimes\mathbb C^2$ 上。对每个 $S\in\{0,\frac12,1,\dots\}$，向量
$$
v\otimes|{\uparrow}\rangle\quad(v\in B_{S-1/2}),\qquad
\sqrt{\tfrac{2S+1}{2S+2}}\;v\otimes|{\downarrow}\rangle-\tfrac1{\sqrt{(2S+1)(2S+2)}}\;T_-v\otimes|{\uparrow}\rangle\quad(v\in B_{S+1/2})
$$
构成 $\ker S_+\cap\ker(S_z-S)$ 的一组规范正交基（对 $S=0$，第一族为空）。从一个自旋出发
（$E=\mathbb C^2$，$B_{1/2}=\{|{\uparrow}\rangle\}$，其余所有 $B_\mu$ 为空）并添加 $L-1$ 个自旋，这就给出 $L$ 个自旋 $\frac12$ 的
最高权向量空间的规范正交基，其维数为 $m_S$（命题 4.15(b)）。

*证明。* 记 $a=\sqrt{(2S+1)/(2S+2)}$，$b=-1/\sqrt{(2S+1)(2S+2)}$，$u(v)=a\,v\otimes|{\downarrow}\rangle+b\,T_-v\otimes|{\uparrow}\rangle$。对
$v\in\mathcal H_{S+1/2}(E)$，$T_+T_-v=(T_-T_++2T_z)v=(2S+1)v$，且 $a+(2S+1)b=0$。*最高权：* $S_z$ 在两族上都等于 $S$
（在 $v$ 上 $T_z=S\mp\frac12$，在 $T_-v$ 上 $T_z=S-\frac12$）；$S_+(v\otimes|{\uparrow}\rangle)=T_+v\otimes|{\uparrow}\rangle=0$，并且
$S_+u(v)=a\,T_+v\otimes|{\downarrow}\rangle+a\,v\otimes|{\uparrow}\rangle+b\,T_+T_-v\otimes|{\uparrow}\rangle=(a+(2S+1)b)\,v\otimes|{\uparrow}\rangle=0$。
*规范正交：* 对 $v\in B_{S-1/2}$，$\langle v\otimes|{\uparrow}\rangle,u(v')\rangle=b\langle v,T_-v'\rangle=b\langle T_+v,v'\rangle=0$，并且
$\langle u(v),u(v')\rangle=a^2\langle v,v'\rangle+b^2\langle v,T_+T_-v'\rangle=(a^2+(2S+1)b^2)\langle v,v'\rangle=\langle v,v'\rangle$。
*张成：* 设 $S_zw=Sw$，$S_+w=0$，写 $w=x\otimes|{\uparrow}\rangle+y\otimes|{\downarrow}\rangle$。则 $T_zx=(S-\frac12)x$，
$T_zy=(S+\frac12)y$，且 $0=S_+w=(T_+x+y)\otimes|{\uparrow}\rangle+T_+y\otimes|{\downarrow}\rangle$，所以 $T_+y=0$，即
$y\in\mathcal H_{S+1/2}(E)$，并且 $T_+x=-y$。于是 $x'=x+\frac1{2S+1}T_-y$ 满足 $T_+x'=-y+y=0$，所以
$x'\in\mathcal H_{S-1/2}(E)$，而对 $S=0$ 它为 $0$（命题 4.15(b) 的证明）。由于 $b/a=-1/(2S+1)$，
$w=x'\otimes|{\uparrow}\rangle+a^{-1}u(y)$（$u$ 线性延拓）属于两族的张成。对 $L$ 个自旋，对 $n=2,\dots,L$ 把这一结论应用于
$E=(\mathbb C^2)^{\otimes(n-1)}$ 及其总自旋 $T_a$；在一个自旋上，$\mathcal H_{1/2}=\mathbb C|{\uparrow}\rangle$，其余所有
$\mathcal H_\mu=0$（唯一的自旋是 $\frac12$）。$\square$

第二族由定理 4.16 中 $j=S+\frac12$ 时的向量 $|S,S\rangle\!\rangle$ 组成，它们写在由诸 $T_-^kv$ 张成的自旋 $j$ 拷贝之中
（定理 4.11 的证明，$|j,j-1\rangle=T_-v/\sqrt{2j}$），而第一族是 $j=S-\frac12$ 时的 $|j+\frac12,j+\frac12\rangle\!\rangle$。用中间
自旋序列 $S_1=\frac12,S_2,\dots,S_L=S$（$|S_n-S_{n-1}|=\frac12$）标记每个基向量，就得到 Yamanouchi 基。

*库中：* `su2.SU2Chain` 用这些系数构造 Yamanouchi 基，它位于 `SpinChain` 的 $U(1)$ 扇区 $\hat N=\frac L2+S$ 之中；推论 4.17 在八个自旋上得到检验，同时给出 Heisenberg 环的 $U(1)$ 块与 $SU(2)$ 块及其代价——[教程，4.6 节](../../tutorial/zh/04_symmetries.md#s-4-6)。

> **物理意义**（物理解释，而非本课程的定理）。Heisenberg 模型描述磁性绝缘体中局域电子自旋的交换相互作用；
> $J_{ij}>0$ 倾向于使相邻自旋反平行。它的谱由多重态组成，而磁场 $h\,S^z$ 把每个多重态分裂为 $2S+1$ 个等间距
> 能级（习题 4）。定理 9.9（二分耦合；对偶数长度的反铁磁链和环，基态是单态）可追溯到 E. Lieb 和 D. Mattis
> （*J. Math. Phys.* 3 (1962) 749）。

<a id="s-4-7"></a>

## 4.7 实代数：实型、复型与四元数型

对称性常常由实矩阵给出：$\mathbb R^3$ 的转动、置换、在某组基下为实的哈密顿量。它们生成的复代数由定理 2.6
分解，但其中的*实*矩阵构成一个实代数，其单的部分有三种：$M_k(\mathbb R)$、$M_k(\mathbb C)$ 与 $M_k(\mathbb H)$。
我们由复分解再加上一个数据把它们推导出来：复共轭对扇区做了什么。

对 $a\in M_n(\mathbb C)$，令 $\bar a$ 为逐元素的复共轭（对向量同样处理，对集合则令 $\bar S=\{\bar s:s\in S\}$）。
共轭是反线性的双射，$\overline{ab}=\bar a\bar b$，$\overline{a^*}=\bar a^*$，$\bar{\bar a}=a$；$\bar a=a$ 当且仅当 $a$
是实的，此时 $a^*=a^T$。对矩阵（或向量）的集合 $\mathcal B$，令 $\mathcal B_{\mathbb R}=\mathcal B\cap M_n(\mathbb R)$（或
$\mathcal B\cap\mathbb R^n$）。内积为 $\langle u,v\rangle=\sum_i\bar u_iv_i$，所以
$\langle\bar u,\bar v\rangle=\overline{\langle u,v\rangle}=\langle v,u\rangle$。

<a id="r-4-18"></a>

**引理 4.18（实部）。** (a) 设 $\mathcal B$ 是 $M_n(\mathbb C)$ 或 $\mathbb C^n$ 的复子空间，满足
$\bar{\mathcal B}=\mathcal B$。则 $\mathcal B=\mathcal B_{\mathbb R}\oplus i\mathcal B_{\mathbb R}$，且
$\dim_{\mathbb R}\mathcal B_{\mathbb R}=\dim_{\mathbb C}\mathcal B$。
(b) 设 $\mathcal A$ 是 $\mathbb C^n$ 上满足 $\bar{\mathcal A}=\mathcal A$ 的 $*$-代数。则 $\mathcal A_{\mathbb R}$ 是带单位元、
在转置下封闭的实代数；$\overline{\mathcal A'}=\mathcal A'$；并且 $(\mathcal A')_{\mathbb R}$ 是对所有 $a\in\mathcal A_{\mathbb R}$
满足 $Xa=aX$ 的实矩阵 $X$ 的集合。
(c) 若 $S\subseteq M_n(\mathbb C)$ 满足 $S^*=S$ 与 $\bar S=S$——例如 $S=\{S_1,\dots,S_J,S_1^T,\dots,S_J^T\}$，
其中 $S_j$ 为实矩阵，或者一个由实正交矩阵构成的群——则 $\overline{S''}=S''$。

*证明。* (a) 对 $b\in\mathcal B$，$\mathrm{Re}\,b=\frac12(b+\bar b)$ 与 $\mathrm{Im}\,b=\frac1{2i}(b-\bar b)$ 是
$\mathcal B$ 的实元素，且 $b=\mathrm{Re}\,b+i\mathrm{Im}\,b$；而 $\mathcal B_{\mathbb R}\cap i\mathcal B_{\mathbb R}$ 中的元素
既是实的又是纯虚的，因而为 $0$。所以实空间 $\mathcal B_{\mathbb R}$ 的一组基 $b_1,\dots,b_r$ 在 $\mathbb C$ 上张成
$\mathcal B$，并且它在 $\mathbb C$ 上线性无关：若 $\sum_j(\alpha_j+i\beta_j)b_j=0$，其中 $\alpha_j,\beta_j$ 为实数，取实部与虚部
即得 $\sum_j\alpha_jb_j=\sum_j\beta_jb_j=0$。(b) $\mathcal A_{\mathbb R}$ 含有 $1$，在乘积下封闭，也在
$a\mapsto a^T=a^*$ 下封闭。$Y$ 与所有 $t\in T$ 可交换当且仅当 $\bar Y$ 与所有 $\bar t$ 可交换，所以对任意集合 $T$ 有
$\overline{T'}=\bar T'$，从而 $\overline{\mathcal A'}=\bar{\mathcal A}'=\mathcal A'$。与 $\mathcal A_{\mathbb R}$ 可交换的实矩阵 $X$
与 $\mathcal A_{\mathbb R}+i\mathcal A_{\mathbb R}=\mathcal A$ 可交换（由 (a)）。(c) 两次应用
$\overline{T'}=\bar T'$：$\overline{S''}=\bar S''=S''$。$\square$

**通用记号。** 在本节余下部分中，$\mathcal A$ 是 $\mathbb C^n$ 上满足 $\bar{\mathcal A}=\mathcal A$ 的 $*$-代数，$V$ 是如
定理 2.6 中的分解。$V_c$ 是 $V$ 中属于第 $c$ 个直和项的 $n\times k_cm_c$ 列块，所以 $V_c^*V_d=\delta_{cd}1$，
$V_cV_c^*=P_c$ 为极小中心投影，并且
$$
aV_c=V_c\big(\pi_c(a)\otimes1_{m_c}\big)\quad(a\in\mathcal A),\qquad XV_c=V_c\big(1_{k_c}\otimes x_c\big)\quad(X\in\mathcal A'),
$$
其中 $\pi_c(a)\in M_{k_c}$，$x_c\in M_{m_c}$。每个 $\pi_c$ 都是 $\mathcal A$ 到 $M_{k_c}$ 上的保单位的 $*$-同态（满射）：
对 $x\in M_{k_c}$，$\mathcal A$ 的元素 $V_c(x\otimes1)V_c^*$ 满足 $\pi_c=x$，且对 $e\neq c$ 有 $\pi_e=0$。
最后令 $\Omega_2=\begin{pmatrix}0&-1\\1&0\end{pmatrix}$，$\Omega_{2r}=1_r\otimes\Omega_2$（含 $r$ 个 $\Omega_2$ 拷贝的分块
对角矩阵）。

子空间 $E\subseteq\mathbb C^k$ 上的映射 $\theta:E\to E$ 称为**反酉的**，若
$\theta(\alpha u+\beta v)=\bar\alpha\theta u+\bar\beta\theta v$ 且 $\langle\theta u,\theta v\rangle=\langle v,u\rangle$。对
酉矩阵 $J\in U(k)$，$\theta_Jv=J\bar v$ 是 $\mathbb C^k$ 上的反酉映射
（$\langle J\bar u,J\bar v\rangle=\langle\bar u,\bar v\rangle=\langle v,u\rangle$），且 $\theta_J^2v=J\bar Jv$。

<a id="r-4-19"></a>

**引理 4.19（反酉映射）。** 设 $\theta$ 是 $E\neq0$ 上的反酉映射，$\theta^2=\varepsilon\,\mathrm{id}$，
$\varepsilon\in\{1,-1\}$。
(a) 若 $\varepsilon=1$，则 $E$ 有一组由满足 $\theta w=w$ 的向量 $w$ 构成的规范正交基。
(b) 若 $\varepsilon=-1$，则对所有 $u\in E$ 有 $\langle u,\theta u\rangle=0$；对每个单位向量 $u$，$\mathrm{span}\{u,\theta u\}$
是 $\theta$-不变平面；并且 $E$ 有规范正交基 $w_1,\theta w_1,\dots,w_r,\theta w_r$。特别地，
$\dim E=2r$ 是偶数。
(c) 设 $J\in U(k)$ 满足 $J\bar J=\varepsilon1$。若 $\varepsilon=1$，则存在酉矩阵 $W$ 使 $J\bar W=W$，即
$J=WW^T$。若 $\varepsilon=-1$，则 $k$ 为偶数，并存在酉矩阵 $W$ 使 $J\bar W=W\Omega_k$，即 $J=W\Omega_kW^T$。

*证明。* (a) 令 $F=\{v\in E:\theta v=v\}$，这是一个实子空间。每个 $v\in E$ 都可写成 $v_1+iv_2$，其中
$v_1=\frac12(v+\theta v)$，$v_2=\frac1{2i}(v-\theta v)$，两者都属于 $F$
（$\theta v_2=\frac1{-2i}(\theta v-v)=v_2$）。对 $u,v\in F$，$\langle u,v\rangle=\langle\theta u,\theta v\rangle=\langle v,u\rangle$ 是
实数，所以 $\langle\cdot,\cdot\rangle$ 是 $F$ 上的实内积；设 $w_1,\dots,w_s$ 是 $F$ 的一组规范正交基
（Gram–Schmidt）。诸 $w_j$ 在 $E$ 中规范正交，因而在 $\mathbb C$ 上线性无关，并且由于 $E=F+iF$，它们在 $\mathbb C$
上张成 $E$。(b) 把 $\langle\theta x,\theta y\rangle=\langle y,x\rangle$ 应用于 $x=\theta u$，$y=u$：
$\langle u,\theta u\rangle=\langle\theta^2u,\theta u\rangle=-\langle u,\theta u\rangle$。所以对单位向量 $u$，$\theta u$ 是与 $u$
正交的单位向量，且 $E_1=\mathrm{span}\{u,\theta u\}$ 是不变的：
$\theta(\alpha u+\beta\theta u)=\bar\alpha\theta u-\bar\beta u$。$E_1$ 在 $E$ 中的正交补也是不变的：对 $x\perp E_1$ 与
$y\in E_1$，$y=-\theta(\theta y)$，且
$\langle\theta x,y\rangle=-\langle\theta x,\theta(\theta y)\rangle=-\langle\theta y,x\rangle=0$，因为 $\theta y\in E_1$。对
$\dim E$ 归纳即得这组基。(c) $\theta_J$ 是 $\mathbb C^k$ 上的反酉映射，$\theta_J^2=J\bar J=\varepsilon$。设 $W$ 是以 (a)
的基、或者 (b) 的基 $w_1,\theta_Jw_1,w_2,\theta_Jw_2,\dots$ 为列的酉矩阵。$J\bar W$ 的列是 $W$ 的列在
$\theta_J$ 下的像。在情形 (a) 中这就是 $J\bar W=W$。在情形 (b) 中，像为 $\theta_Jw_j$ 与 $\theta_J^2w_j=-w_j$，它们正是
$W\Omega_k$ 的第 $2j-1$ 列与第 $2j$ 列（因为 $\Omega_2e_1=e_2$，$\Omega_2e_2=-e_1$）。由于 $\bar W^{-1}=\overline{W^*}=W^T$，
得 $J=WW^T$，或相应地 $J=W\Omega_kW^T$。$\square$

<a id="r-4-20"></a>

**定理 4.20（扇区的共轭）。** (a) 存在扇区的一个满足 $\sigma\circ\sigma=\mathrm{id}$ 的置换 $\sigma$，以及酉矩阵
$J_c\in U(k_c)$，$M_c\in U(m_c)$，使得 $\bar P_c=P_{\sigma(c)}$，$k_{\sigma(c)}=k_c$，$m_{\sigma(c)}=m_c$，并且对每个 $c$ 与
所有 $a\in\mathcal A$，
$$
\overline{V_c}=V_{\sigma(c)}\,(J_c\otimes M_c),\qquad\pi_{\sigma(c)}(\bar a)=J_c\,\overline{\pi_c(a)}\,J_c^*.
$$
$J_c\otimes M_c$ 是唯一的，而 $J_c,M_c$ 在 $(J_c,M_c)\mapsto(\lambda J_c,\bar\lambda M_c)$（$|\lambda|=1$）的意义下唯一。
(b) 若 $\sigma(c)=c$，则存在 $\varepsilon_c\in\{1,-1\}$，使 $J_c\bar J_c=\varepsilon_c1_{k_c}$ 且
$M_c\bar M_c=\varepsilon_c1_{m_c}$。若 $\varepsilon_c=-1$，则 $k_c$ 与 $m_c$ 都是偶数。
(c) $\sigma$ 与诸 $\varepsilon_c$ 不依赖于 $V$ 的选取。
(d) 对 $y\in M_{m_c}$，$\overline{V_c(1\otimes y)V_c^*}=V_{\sigma(c)}\big(1\otimes M_c\bar yM_c^*\big)V_{\sigma(c)}^*$。因此，若
$X\in\mathcal A'$ 的块为 $x_c$，则 $\bar X\in\mathcal A'$ 的块为 $(\bar X)_{\sigma(c)}=M_c\bar x_cM_c^*$：在交换子上，
$M_c$ 扮演 $J_c$ 在 $\mathcal A$ 上所扮演的角色。

*证明。* (a) 共轭把 $\mathcal A$ 映满 $\mathcal A$，并保持和、乘积与伴随，所以它把中心 $Z(\mathcal A)$ 映满自身，并把
$Z(\mathcal A)$ 的非零投影的集合双射地映到自身，且保持序 $p\le q$（$pq=p$）。因此它置换极小投影，而极小投影就是
诸 $P_c$（引理 2.3）：$\bar P_c=P_{\sigma(c)}$，并且由于 $\bar{\bar P}_c=P_c$，$\sigma\circ\sigma=\mathrm{id}$。令 $d=\sigma(c)$。
$\overline{V_c}$ 的像是 $P_c$ 的像的共轭，即 $\bar P_c=P_d$ 的像；所以
$\overline{V_c}=P_d\overline{V_c}=V_dT$，其中 $T=V_d^*\overline{V_c}$，且
$T^*T=\overline{V_c}^*P_d\overline{V_c}=\overline{V_c^*V_c}=1$。由于 $k_cm_c=\mathrm{rank}\,P_c=\mathrm{rank}\,\bar P_c=k_dm_d$，
等距算子 $T$ 是方阵，因而是酉的。对 $aV_c=V_c(\pi_c(a)\otimes1)$ 取共轭，并利用
$\bar aV_d=V_d(\pi_d(\bar a)\otimes1)$，得 $V_d(\pi_d(\bar a)\otimes1)T=V_dT(\overline{\pi_c(a)}\otimes1)$；左乘
$V_d^*$，
$$
\big(\pi_d(\bar a)\otimes1_{m_d}\big)\,T=T\,\big(\overline{\pi_c(a)}\otimes1_{m_c}\big)\qquad(a\in\mathcal A).
$$
把 $T$ 写成 $T=\sum_{j,l}T_{jl}\otimes|j\rangle\langle l|$，其中 $T_{jl}$ 是 $k_d\times k_c$ 块。上式说的是
$\pi_d(\bar a)T_{jl}=T_{jl}\overline{\pi_c(a)}$：每个 $T_{jl}$ 都是从 $\rho(a)=\overline{\pi_c(a)}$ 到
$\rho'(a)=\pi_d(\bar a)$ 的缠结算子。它们是带有伴随运算的 $\mathcal A$ 在 4.1 节意义下的表示
（$\rho(a^*)=\rho(a)^*$，$\rho'(a^*)=\rho'(a)^*$），并且是不可约的，因为它们的像分别是整个 $M_{k_c}$ 与
$M_{k_d}$。由于 $T\neq0$，某个 $T_{jl}\neq0$，所以 $k_c=k_d$，且对某个酉矩阵 $J_c$ 有 $\mathrm{Hom}(\rho,\rho')=\mathbb CJ_c$
（引理 4.2）。因此 $T_{jl}=\mu_{jl}J_c$，即 $T=J_c\otimes M_c$，其中 $M_c=(\mu_{jl})$；$1=T^*T=1\otimes M_c^*M_c$，
且 $m_c=k_cm_c/k_c=m_d$，所以 $M_c$ 是酉的，而上式就变成 $\pi_d(\bar a)$ 的公式。$T$ 由 $V$ 决定；
若又有 $T=J'\otimes M'$，其中 $J'$ 为酉矩阵，则块 $T_{jl}=M'_{jl}J'$ 表明 $J'\in\mathrm{Hom}(\rho,\rho')$，所以
$J'=\lambda J_c$，$|\lambda|=1$，进而 $M'=\bar\lambda M_c$。

(b) 设 $\sigma(c)=c$，并令 $J=J_c$，$M=M_c$。对 $\overline{V_c}=V_c(J\otimes M)$ 取共轭得
$V_c=\overline{V_c}(\bar J\otimes\bar M)=V_c(J\bar J\otimes M\bar M)$，所以 $J\bar J\otimes M\bar M=1$。两次应用 (a) 的公式，并利用
$\overline{J^*}=J^T$，
$$
\pi_c(a)=\pi_c(\bar{\bar a})=J\,\overline{J\,\overline{\pi_c(a)}\,J^*}\,J^*=(J\bar J)\,\pi_c(a)\,(J\bar J)^*,
$$
所以酉矩阵 $J\bar J$ 与 $\pi_c(\mathcal A)=M_{k_c}$ 可交换，因而是标量 $\lambda1$（引理 2.5 第 5 步，取
$m=1$），$|\lambda|=1$。于是 $\bar JJ=\overline{J\bar J}=\bar\lambda1$，且 $\lambda J=(J\bar J)J=J(\bar JJ)=\bar\lambda J$，所以
$\lambda=\bar\lambda\in\{1,-1\}$；记之为 $\varepsilon_c$。它在 $J\mapsto\lambda J$ 下不变，而
$\varepsilon_c1\otimes M\bar M=1$ 给出 $M\bar M=\varepsilon_c1$。若 $\varepsilon_c=-1$，则
$(-1)^{k_c}=\det(J\bar J)=\det J\cdot\overline{\det J}=|\det J|^2>0$，所以 $k_c$ 是偶数，$m_c$ 同理。

(c) 诸 $P_c$ 不依赖于 $V$（定理 2.6），所以 $\sigma$ 也不依赖于 $V$。设 $V'$ 是另一个分解，其编号使得
$V'_cV_c'^*=P_c$（定理 2.6），相应的同态为 $\pi'_c$，数据为 $J'_c,M'_c$，并设 $\sigma(c)=c$。则
$V'_c=P_cV'_c=V_cY$，其中 $Y=V_c^*V'_c$ 是方阵形式的等距算子，因而是酉的，并且对 $a\in\mathcal A$ 有
$(\pi_c(a)\otimes1)Y=Y(\pi'_c(a)\otimes1)$。把 (a) 的论证用于不可约表示 $\pi'_c,\pi_c$（代替 $\rho,\rho'$），得
$Y=u\otimes w$，其中 $u,w$ 为酉矩阵。因此
$\overline{V'_c}=\overline{V_c}(\bar u\otimes\bar w)=V_c(J_c\bar u\otimes M_c\bar w)=V'_c(u^*J_c\bar u\otimes w^*M_c\bar w)$，所以由 (a) 中的
唯一性，$J'_c=\lambda u^*J_c\bar u$，$|\lambda|=1$，并且由于 $\bar uu^T=\overline{uu^*}=1$，
$J'_c\bar J'_c=u^*J_c\bar uu^T\bar J_cu=u^*J_c\bar J_cu=\varepsilon_c1$。

(d) 由 (a)，$\overline{V_c(1\otimes y)V_c^*}=V_{\sigma(c)}(J_c\otimes M_c)(1\otimes\bar y)(J_c\otimes M_c)^*V_{\sigma(c)}^*$，并且
$(J_c\otimes M_c)(1\otimes\bar y)(J_c\otimes M_c)^*=J_cJ_c^*\otimes M_c\bar yM_c^*=1\otimes M_c\bar yM_c^*$。对
$X=\sum_cV_c(1\otimes x_c)V_c^*$ 关于 $c$ 求和即可。$\square$

若 $\sigma(c)=c$ 且 $\varepsilon_c=1$，则称扇区 $c$ 为**实型**的；若 $\sigma(c)=c$ 且 $\varepsilon_c=-1$，则称为
**四元数型**的；若 $\sigma(c)\neq c$，则称为**复型**的。由定理 4.20(c)，这是 $\mathcal A$ 与 $P_c$ 的性质。由定理
4.20(b)，满足 $\sigma(c)=c$ 且 $k_c$ 为奇数或 $m_c$ 为奇数的扇区是实型的。

*库中：* $\sigma$、$J_c$、$M_c$、$\varepsilon_c$ 以及各扇区的类型由 `Decomposition.real_structure()` 计算——[教程，4.7 节](../../tutorial/zh/04_symmetries.md#s-4-7)。

<a id="r-4-21"></a>

**定理 4.21（实代数）。** 设 $\Gamma$ 包含每个满足 $\sigma(c)=c$ 的扇区 $c$，以及每一对满足 $\sigma(c)\neq c$ 的
$\{c,\sigma(c)\}$ 中的一个扇区。对 $\sigma(c)=c$，令 $R_c=\{x\in M_{k_c}(\mathbb C):xJ_c=J_c\bar x\}$。则
$a\mapsto(\pi_c(a))_{c\in\Gamma}$ 是实代数的同构，它把转置映为伴随，
$$
\mathcal A_{\mathbb R}\;\cong\bigoplus_{c\in\Gamma,\ \sigma(c)=c}R_c\;\oplus\bigoplus_{c\in\Gamma,\ \sigma(c)\neq c}M_{k_c}(\mathbb C),
$$
其中最后的那些直和项视为实代数（**复型**）。此外：
(a)（*实型*）若 $\varepsilon_c=1$，写 $J_c=W_JW_J^T$，$M_c=W_MW_M^T$，其中 $J_c\bar W_J=W_J$，$M_c\bar W_M=W_M$
（引理 4.19(c)）。则 $x\mapsto W_J^*xW_J$ 把 $R_c$ 同构地映到 $M_{k_c}(\mathbb R)$ 上。
$O_c=V_c(W_J\otimes W_M)$ 的列是实的且规范正交，它们在 $\mathbb C$ 上张成 $P_c$ 的像，并且对 $a\in\mathcal A_{\mathbb R}$，
$aO_c=O_c(\alpha_c(a)\otimes1_{m_c})$，其中 $\alpha_c(a)=W_J^*\pi_c(a)W_J$ 是实矩阵。
(b)（*四元数型*）若 $\varepsilon_c=-1$，写 $J_c=W_J\Omega_{k_c}W_J^T$，其中 $J_c\bar W_J=W_J\Omega_{k_c}$
（引理 4.19(c)）。则 $x\mapsto\mathbb M^{-1}(W_J^*xW_J)$ 是实代数的同构 $R_c\cong M_{k_c/2}(\mathbb H)$，它把伴随
映为四元数伴随 $(q^*)_{ij}=\bar q_{ji}$；这里 $\mathbb M:M_r(\mathbb H)\to M_{2r}(\mathbb C)$ 把每个矩阵元 $q$
换成引理 E.15(2) 的 $2\times2$ 块 $M(q)$，其像为
$\{y\in M_{2r}(\mathbb C):y\Omega_{2r}=\Omega_{2r}\bar y\}$。
特别地，$\dim_{\mathbb R}\mathcal A_{\mathbb R}=\sum_ck_c^2=\dim_{\mathbb C}\mathcal A$。

*证明。* 对实的 $a$ 有 $\bar a=a$，于是定理 4.20(a) 成为 $\pi_{\sigma(c)}(a)=J_c\overline{\pi_c(a)}J_c^*$。对
$\sigma(c)=c$，这说明 $\pi_c(a)\in R_c$。所以 $\Phi(a)=(\pi_c(a))_{c\in\Gamma}$ 把 $\mathcal A_{\mathbb R}$ 映入右边；它是
实代数的同态，且 $\Phi(a^T)=\Phi(a^*)=\Phi(a)^*$。若 $\Phi(a)=0$，则同一公式对 $c\in\Gamma$ 给出
$\pi_{\sigma(c)}(a)=0$，所以所有 $\pi_e(a)=0$，从而 $a=0$。*满射：* 设 $c\in\Gamma$，当 $\sigma(c)=c$ 时取 $x\in R_c$，
否则取 $x\in M_{k_c}$，并令 $a=V_c(x\otimes1)V_c^*\in\mathcal A$。由对扇区 $\sigma(e)$ 的定理 4.20(a)，
$\pi_e(\bar a)=J_{\sigma(e)}\overline{\pi_{\sigma(e)}(a)}J_{\sigma(e)}^*$，除非 $\sigma(e)=c$，它都为 $0$。若
$\sigma(c)=c$，则 $\pi_c(\bar a)=J_c\bar xJ_c^*=x$（$R_c$ 的定义），所以实矩阵 $\mathrm{Re}\,a=\frac12(a+\bar a)\in\mathcal A_{\mathbb R}$
满足：$\Phi(\mathrm{Re}\,a)$ 在位置 $c$ 上为 $x$，在其他位置上为 $0$。若 $\sigma(c)\neq c$，则 $\pi_c(\bar a)=0$，而 $2\mathrm{Re}\,a$
满足 $\pi_c=x$，且对 $e\notin\{c,\sigma(c)\}$ 有 $\pi_e=0$；由于 $\sigma(c)\notin\Gamma$，$\Phi(2\mathrm{Re}\,a)$ 在位置
$c$ 上为 $x$，在其他位置上为 $0$。这些元素之和给出整个右边。

(a) 用 $W_J^*$ 左乘、用 $\bar W_J$ 右乘 $xW_JW_J^T=W_JW_J^T\bar x$；由于
$W_J^T\bar W_J=\overline{W_J^*W_J}=1$，它变为 $W_J^*xW_J=W_J^T\bar x\bar W_J=\overline{W_J^*xW_J}$，而这些步骤可以
逆推。所以 $R_c=W_JM_{k_c}(\mathbb R)W_J^*$，并且 $x\mapsto W_J^*xW_J$ 是与伴随相容的实代数同构。其次
$\overline{O_c}=\overline{V_c}(\bar W_J\otimes\bar W_M)=V_c(J_c\bar W_J\otimes M_c\bar W_M)=V_c(W_J\otimes W_M)=O_c$，$O_c^*O_c=1$，并且
$aO_c=V_c(\pi_c(a)W_J\otimes W_M)=O_c(W_J^*\pi_c(a)W_J\otimes1)$，其中对 $a\in\mathcal A_{\mathbb R}$，$W_J^*\pi_c(a)W_J$ 是实的。

(b) 对 $J_c=W_J\Omega W_J^T$（$\Omega=\Omega_{k_c}$）作同样的计算表明，$xJ_c=J_c\bar x$ 当且仅当 $y=W_J^*xW_J$ 满足
$y\Omega=\Omega\bar y$。把 $y$ 写成 $2\times2$ 块 $y_{ij}$；由于 $\Omega$ 是分块对角的，$y\Omega=\Omega\bar y$ 当且仅当
对所有 $i,j$ 有 $y_{ij}\Omega_2=\Omega_2\bar y_{ij}$，而对 $y_{ij}=\begin{pmatrix}\alpha&\beta\\\gamma&\delta\end{pmatrix}$，这就是
$\begin{pmatrix}\beta&-\alpha\\\delta&-\gamma\end{pmatrix}=\begin{pmatrix}-\bar\gamma&-\bar\delta\\\bar\alpha&\bar\beta\end{pmatrix}$，即
$\delta=\bar\alpha$ 且 $\beta=-\bar\gamma$：$y_{ij}=\begin{pmatrix}\alpha&-\bar\gamma\\\gamma&\bar\alpha\end{pmatrix}=M(\alpha,\bar\gamma)$，其中
$M(a,b)=\begin{pmatrix}a&-b\\\bar b&\bar a\end{pmatrix}$。所以 $\mathbb M$ 的像如所述。由引理 E.15(2)，$M$ 是满足
$M(\bar q)=M(q)^*$ 的实线性单同态 $\mathbb H\to M_2(\mathbb C)$；由于矩阵乘积只用到矩阵元的和与积，
$\mathbb M$ 是满足 $\mathbb M(q^*)=\mathbb M(q)^*$ 的实代数单同态（同定理 E.24(1) 的证明）。因此 $x\mapsto\mathbb M^{-1}(W_J^*xW_J)$
就是所断言的同构。*维数：* 两种情形下都有 $\dim_{\mathbb R}R_c=k_c^2$（$\dim_{\mathbb R}M_r(\mathbb H)=4r^2$，
$2r=k_c$），而一对扇区贡献 $\dim_{\mathbb R}M_{k_c}(\mathbb C)=2k_c^2=k_c^2+k_{\sigma(c)}^2$。$\square$

交换子也是如此，只是两个因子互换：$(\mathcal A')_{\mathbb R}$ 是在同样的扇区上、以同样的类型取的
$M_{m_c}(\mathbb R)$、$M_{m_c/2}(\mathbb H)$ 与 $M_{m_c}(\mathbb C)$ 的直和。事实上，$\mathcal A'$ 是满足
$\overline{\mathcal A'}=\mathcal A'$ 的 $*$-代数（引理 4.18(b)），而上面的证明只用到三个事实，它们在把 $\pi_c$ 换成
$X\mapsto x_c$、把 $J_c$ 换成 $M_c$ 后对 $\mathcal A'$ 成立：$X\mapsto x_c$ 是到 $M_{m_c}$ 上的保单位的 $*$-同态
（$V_c(1\otimes y)V_c^*\in\mathcal A'$ 在扇区 $c$ 中的块为 $y$，在其他扇区中为 $0$），$(\bar X)_{\sigma(c)}=M_c\bar x_cM_c^*$
（定理 4.20(d)），以及对 $\sigma(c)=c$ 有 $M_c\bar M_c=\varepsilon_c1$（定理 4.20(b)）。

*库中：* 适配于定理 4.21 的实正交基，`Decomposition.real_basis()`——[教程，4.7 节](../../tutorial/zh/04_symmetries.md#s-4-7)。

三个例子。(i) $\mathbb R^3$ 的转动的生成元 $(L_a)_{bc}=-\epsilon_{abc}$，以
$L_a\otimes1_2$ 的方式作用于 $\mathbb R^3\otimes\mathbb R^2$：只有一个扇区，$k=3$，$m=2$（$iL_x,iL_y,iL_z$ 是 $\mathbb C^3$ 上的一个自旋表示，其中 $iL_z$ 有单重特征值
$1,0,-1$，所以由推论 4.12 它是自旋 $1$）；由于 $k$ 是奇数，它是
实型的，且 $\mathcal A_{\mathbb R}\cong M_3(\mathbb R)$。(ii) $\mathbb R^2$ 的转动生成元 $\Omega_2$，以
$\Omega_2\otimes1_2$ 的方式作用于 $\mathbb R^2\otimes\mathbb R^2$：由于 $\Omega_2^2=-1$ 且
$\Omega_2^*=-\Omega_2$，该代数为 $\mathrm{span}\{1,\Omega_2\otimes1\}$；$\Omega_2$ 对应于 $\pm i$ 的特征向量
$(1,\mp i)/\sqrt2$ 互为共轭，它们给出两个 $k=1$，$m=2$ 的扇区，在共轭下互换：这是一对复型扇区，且
$\mathcal A_{\mathbb R}=\{\alpha+\beta\,\Omega_2\otimes1\}\cong\mathbb C$。
(iii) 四元数 $q=(a,b)\in\mathbb H$（E.6 节）在 $\mathbb H=\mathbb R^4$ 上的左乘 $L_q$，坐标为
$(\mathrm{Re}\,a,\mathrm{Im}\,a,\mathrm{Re}\,b,\mathrm{Im}\,b)$。$q\mapsto L_q$ 是实代数的单同态
（$\mathbb H$ 是结合的，引理 E.15(2)），且 $L_{\bar q}=L_q^T$，因为对范数 $N(a,b)=|a|^2+|b|^2$（它就是
$\mathbb R^4$ 的欧几里得范数）有 $\langle qx,y\rangle=\langle x,\bar qy\rangle$（引理 E.14(2)）。所以由
$L_I,L_J,L_K$（$I=(i,0)$，$J=(0,1)$，$K=(0,i)$）及其转置生成的实代数是
$\{L_q:q\in\mathbb H\}\cong\mathbb H$；复代数 $\mathcal A$ 是生成元乘积的复线性张成，即
$\{L_q\}+i\{L_q\}$，所以 $\mathcal A_{\mathbb R}=\{L_q\}\cong\mathbb H$（引理 4.18(a)）。这已经确定了分解：
$\sum_ck_c^2=\dim_{\mathbb C}\mathcal A=\dim_{\mathbb R}\mathcal A_{\mathbb R}=4$（定理 4.21），并且并非所有 $k_c=1$，否则
$\mathcal A$，从而 $\mathcal A_{\mathbb R}\cong\mathbb H$，将是交换的；所以只有一个扇区，$k=2$ 且
$\sigma(c)=c$，并且它不是实型的，否则定理 4.21(a) 将给出 $\mathcal A_{\mathbb R}\cong M_2(\mathbb R)$，而它有零因子，
$\mathbb H$ 却没有零因子（引理 E.15(3)）：它是四元数型的，在 $\mathbb R^4$ 上 $m=4/2=2$。（单位四元数在这里作为 $SU(2)$
作用：$M(q)$ 是行列式为 $|a|^2+|b|^2=1$ 的酉矩阵，并且每个 $U=\begin{pmatrix}\alpha&\beta\\\gamma&\delta\end{pmatrix}\in SU(2)$
都具有这种形式，因为 $U^*=U^{-1}=\begin{pmatrix}\delta&-\beta\\-\gamma&\alpha\end{pmatrix}$ 给出 $\delta=\bar\alpha$，$\beta=-\bar\gamma$，
所以 $U=M(\alpha,\bar\gamma)$。）在 $\mathbb H^2=\mathbb R^8$ 上，$q$ 以 $1_2\otimes L_q$ 的方式作用，同样的论证给出一个
$k=2$，$m=4$ 的四元数型扇区，$\mathcal A_{\mathbb R}\cong M_1(\mathbb H)=\mathbb H$。

*计算：* 这三个例子，使用 `Decomposition.real_structure()` 和 `Decomposition.real_basis()`——[教程，4.7 节](../../tutorial/zh/04_symmetries.md#s-4-7)。

在由实矩阵给出的群表示论中，下面是 Schur 引理的实形式；它表明这三种类型可以从实表示本身读出。

<a id="r-4-22"></a>

**定理 4.22（实 Schur 引理）。** 设 $U\subseteq\mathbb R^n$ 是非零子空间，对所有 $a\in\mathcal A_{\mathbb R}$ 满足
$aU\subseteq U$，并且是**不可约的**：$0$ 与 $U$ 是它仅有的具有这一性质的子空间。令 $D_U$ 为满足
$Tau=aTu$（对所有 $a\in\mathcal A_{\mathbb R}$，$u\in U$）的实线性映射 $T:U\to U$ 的集合。
(a) $D_U$ 是带单位元、没有零因子的结合实代数，因而同构于 $\mathbb R$、$\mathbb C$ 或 $\mathbb H$（定理 E.19）。
(b) 对 $\sigma$ 的恰好一个轨道 $\mathcal O=\{c,\sigma(c)\}$，$U$ 落在 $P_{\mathcal O}=\sum_{e\in\mathcal O}P_e$ 的像中，并且
按 $c$ 为实型、复型或四元数型，分别有 $\dim_{\mathbb R}D_U=1,2,4$ 与 $\dim_{\mathbb R}U=k_c,2k_c,2k_c$。所以
分别有 $D_U\cong\mathbb R,\mathbb C,\mathbb H$。
(c) 对每个轨道 $\mathcal O$，在 $P_{\mathcal O}$ 的像中存在这样的 $U$。$\mathbb R^n$ 本身不可约，当且仅当 $\sigma$ 只有
一个轨道 $\{c,\sigma(c)\}$，且 $m_c=1$（实型与复型），或相应地 $m_c=2$（四元数型）。

*证明。* (a) 对 $T\in D_U$，$\ker T$ 与 $TU$ 是不变子空间，所以 $T=0$ 或 $T$ 是双射，此时
$T^{-1}\in D_U$。$D_U$ 在复合下封闭，复合是结合的，$D_U$ 含有 $\mathrm{id}_U\neq0$，并且由 $ST=0$ 与
$S\neq0$ 推出 $T=S^{-1}ST=0$。

(b) *不变子空间与实投影。* 若 $U'\subseteq\mathbb R^n$ 是不变的，则它在 $\mathbb R^n$ 中的正交补也是不变的（对
$y\perp U'$，$u\in U'$ 与 $a\in\mathcal A_{\mathbb R}$，$u\cdot ay=a^Tu\cdot y=0$，因为
$a^T\in\mathcal A_{\mathbb R}$），所以到 $U'$ 上的正交投影 $e$（一个实对称幂等矩阵）与 $\mathcal A_{\mathbb R}$ 可交换；
由引理 4.18(b)，$e\in(\mathcal A')_{\mathbb R}$。反之，对每个投影 $e\in(\mathcal A')_{\mathbb R}$，$e\mathbb R^n$ 是不变的。
所以 $U'\mapsto e$ 是从不变子空间到 $(\mathcal A')_{\mathbb R}$ 中投影的集合 $\mathcal P$ 上的保持包含关系的双射，而
$U$ 不可约当且仅当其投影 $e$ 是 $\mathcal P$ 的极小非零元素。把引理 4.18(a) 用于像 $e\mathbb C^n$（其实向量为
$e\mathbb R^n$），得 $\dim_{\mathbb R}U=\mathrm{rank}\,e$。

*$D_U\cong e(\mathcal A')_{\mathbb R}e$。* 对 $T\in D_U$，令 $\hat T$ 为实 $n\times n$ 矩阵 $x\mapsto T(ex)$。对
$a\in\mathcal A_{\mathbb R}$，$\hat Tax=T(eax)=T(aex)=aT(ex)=a\hat Tx$，所以 $\hat T\in(\mathcal A')_{\mathbb R}$（引理 4.18(b)），且
$\hat T=e\hat Te$。反之，对 $X\in(\mathcal A')_{\mathbb R}$，$eXe$ 限制在 $U$ 上属于 $D_U$。这两个映射互逆且保持
乘法。此外，$e\mathcal A'e$ 在共轭下封闭，其实元素是 $e(\mathcal A')_{\mathbb R}e$ 中的元素（若 $eXe$ 是实的，则
$eXe=e\mathrm{Re}(X)e$，其中 $\mathrm{Re}\,X\in(\mathcal A')_{\mathbb R}$），所以由引理 4.18(a)，
$\dim_{\mathbb R}D_U=\dim_{\mathbb C}e\mathcal A'e$。

*分解中的投影。* $\mathcal A'$ 中的投影 $e$ 形如 $e=\sum_cV_c(1\otimes f_c)V_c^*$，其中 $f_c\in M_{m_c}$ 为投影
（定理 2.6）。于是 $\mathrm{rank}\,e=\sum_ck_c\mathrm{rank}\,f_c$，
$e\mathcal A'e=\sum_cV_c(1\otimes f_cM_{m_c}f_c)V_c^*$ 的复维数为 $\sum_c(\mathrm{rank}\,f_c)^2$，并且由
定理 4.20(d)，$\bar e=e$ 当且仅当对所有 $c$ 有 $f_{\sigma(c)}=M_c\bar f_cM_c^*$。

*极小元。* 设 $e$ 在 $\mathcal P$ 中极小。对 $\sigma$ 的每个轨道 $\mathcal O$，$P_{\mathcal O}$ 是
$Z(\mathcal A)\subseteq\mathcal A'$ 中的实中心投影（由定理 4.20(a)，$\bar P_{\mathcal O}=P_{\mathcal O}$），它与 $e$ 可交换，所以
$eP_{\mathcal O}\in\mathcal P$ 且 $eP_{\mathcal O}\le e$；由于 $e=\sum_{\mathcal O}eP_{\mathcal O}$，恰有一个轨道 $\mathcal O$
满足 $eP_{\mathcal O}=e$，即 $U=e\mathbb R^n$ 落在 $P_{\mathcal O}$ 的像中，且对 $c\notin\mathcal O$ 有 $f_c=0$。若
$\mathcal O=\{c,d\}$，$d\neq c$，则 $f_d=M_c\bar f_cM_c^*$ 与 $f_c$ 秩相同；对秩一投影 $g\le f_c$，$g'=Y+\bar Y$（其中
$Y=V_c(1\otimes g)V_c^*$）是 $\mathcal A'$ 的两个正交投影之和（分别在扇区 $c$ 与 $d$ 中），它是实的，并且
$g'\le e$（因为 $Y\le e$ 且 $\bar Y\le\bar e=e$）；所以 $g'=e$，$\mathrm{rank}\,f_c=\mathrm{rank}\,f_d=1$，$\dim_{\mathbb C}e\mathcal A'e=2$，$\mathrm{rank}\,e=2k_c$。
若 $\mathcal O=\{c\}$，令 $f=f_c$，$M=M_c$，$\theta v=M\bar v$，这是 $\mathbb C^{m_c}$ 上的反酉映射，
$\theta^2=M\bar M=\varepsilon_c$（定理 4.20(b)）。由 $fM=M\bar f$，$\theta(fv)=M\bar f\bar v=fM\bar v=f\theta v$，所以 $f$ 的
像 $F$ 是 $\theta$-不变的。若 $F_1\subseteq F$ 是 $\theta$-不变的，其投影为 $g$，则 $M\bar gM^*$ 是到
$M\bar F_1=\theta F_1=F_1$ 上的投影（$\theta$ 是单射），所以 $g=M\bar gM^*$，且 $V_c(1\otimes g)V_c^*\in\mathcal P$
位于 $e$ 之下。由极小性，$F$ 除 $0$ 与 $F$ 外没有 $\theta$-不变子空间。若 $\varepsilon_c=1$，$\theta$ 在 $F$ 中的一个
不动单位向量（引理 4.19(a)）张成一个这样的子空间，所以 $\mathrm{rank}\,f=1$，$\dim_{\mathbb C}e\mathcal A'e=1$，
$\mathrm{rank}\,e=k_c$。若 $\varepsilon_c=-1$，单位向量 $u\in F$ 的平面 $\mathrm{span}\{u,\theta u\}$ 就是这样的子空间
（引理 4.19(b)），所以 $\mathrm{rank}\,f=2$，$\dim_{\mathbb C}e\mathcal A'e=4$，$\mathrm{rank}\,e=2k_c$。由于 $\mathbb R$、
$\mathbb C$、$\mathbb H$ 的维数为 $1,2,4$，由 (a) 分别得 $D_U\cong\mathbb R,\mathbb C,\mathbb H$。

(c) $P_{\mathcal O}\in\mathcal P$ 非零，而在它之下、秩最小的 $\mathcal P\setminus\{0\}$ 的元素是极小的。$\mathbb R^n$
不可约当且仅当 $e=1$ 是极小的；它的块为 $f_c=1_{m_c}$，所以由 (b) 的证明，这要求只有一个轨道，并且在三种情形下
分别有 $m_c=1,1,2$。反之，在这些条件下，非零的 $e'\in\mathcal P$ 在该轨道的扇区上有
$f'_c\neq0$（$f'_{\sigma(c)}=M_c\bar f'_cM_c^*$），所以若 $m_c=1$ 则 $f'_c=1$；在四元数型情形下，$f'_c$ 的像是
$\theta$-不变的，因而维数为偶数（引理 4.19(b)），所以 $f'_c=1_2$。于是 $e'=1$。$\square$

对以实正交矩阵作用的群 $G$，$\mathcal A=\mathrm{span}\,\pi(G)$ 在共轭下封闭（引理 4.18(c)），
$\mathcal A_{\mathbb R}$ 是 $\pi(G)$ 的实线性张成（引理 4.18(a)：$\mathcal A_{\mathbb R}$ 由 $\mathcal A$ 中元素的实部
构成），而 $D_U$ 是 $G$ 的不可约实表示的交换子：它是 $\mathbb R$、$\mathbb C$ 或 $\mathbb H$。由 (b) 的证明，复化
$U+iU=e\mathbb C^n=V_c(\mathbb C^{k_c}\otimes F)$（对一对扇区再加上它的共轭）在三种情形下分别承载
$\pi_c$、$\pi_c\oplus\pi_{\sigma(c)}$ 或 $\pi_c\oplus\pi_c$：它是不可约的、两个不等价不可约表示之和，或两个等价
不可约表示之和。

*计算：* 在三个例子的不可约实模 $\mathbb R^3$、$\mathbb R^2$、$\mathbb R^4$ 上验证定理 4.22——[教程，4.7 节](../../tutorial/zh/04_symmetries.md#s-4-7)。

<a id="r-4-23"></a>

**命题 4.23（交换子的实元素）。** 设 $X\in\mathcal A'$ 是实的，其块为 $x_c$。则对每个 $c$ 有
$x_{\sigma(c)}=M_c\bar x_cM_c^*$。若此外 $X=X^T$：
(a) 对实型的 $c$，$y_c=W_M^*x_cW_M$（$W_M$ 如定理 4.21(a)）是实对称矩阵，并且
$XO_c=O_c(1_{k_c}\otimes y_c)$，其中 $O_c$ 是定理 4.21(a) 的实规范正交列；
(b) 对四元数型的 $c$，$x_c$ 的每个特征值都有偶数重数，所以 $X$ 在 $P_c$ 的像上的每个特征值的重数都能被
$2k_c$ 整除（**Kramers 简并**）；
(c) 对复型的 $c$，$x_c$ 与 $x_{\sigma(c)}$ 有相同的特征值及相同的重数。
不假设 $X=X^T$ 时，(a) 中的 $y_c$ 仍然是实的。

*证明。* $X=\bar X$，所以定理 4.20(d) 给出 $x_{\sigma(c)}=M_c\bar x_cM_c^*$。若 $X=X^T$，则 $X^*=X$，且
每个 $x_c$ 都是厄米的。(a) 由 $M_c=W_MW_M^T$，$M_c^*=\bar W_MW_M^*$ 与 $W_M^T\bar W_M=1$：
$\bar y_c=W_M^T\bar x_c\bar W_M=W_M^TM_c^*x_cM_c\bar W_M=W_M^*x_cW_M=y_c$，而 $y_c$ 是厄米的，因而是实对称的。此外
$XO_c=XV_c(W_J\otimes W_M)=V_c(W_J\otimes x_cW_M)=O_c(1\otimes y_c)$。(b) 令 $\theta v=M_c\bar v$，它是反酉的，
$\theta^2=-1$（定理 4.20(b)）。由 $x_cM_c=M_c\bar x_c$，$x_c\theta v=M_c\bar x_c\bar v=\theta(x_cv)$，所以对 $x_c$ 的特征值
$\lambda\in\mathbb R$，特征子空间 $E_\lambda$ 满足 $\theta E_\lambda\subseteq E_\lambda$（$x_c\theta v=\theta(\lambda v)=\lambda\theta v$），
由引理 4.19(b)，$\dim E_\lambda$ 是偶数。在 $P_c$ 的像上，$X$ 作用为 $V_c(1_{k_c}\otimes x_c)V_c^*$，其特征值就是
$x_c$ 的特征值，每个重复 $k_c$ 次。(c) $\bar x_c$ 具有 $x_c$ 的特征值 $\bar\lambda=\lambda$，特征向量为原特征向量的
共轭，而 $x_{\sigma(c)}$ 与 $\bar x_c$ 酉等价。$\square$

在四元数型情形下，$y_c=W_M^*x_cW_M$（其中 $M_c=W_M\Omega_{m_c}W_M^T$ 如引理 4.19(c)）满足
$y_c\Omega=\Omega\bar y_c$（即把定理 4.21(b) 中的计算里的 $J_c$ 换成 $M_c$），所以 $y_c=\mathbb M(h)$，其中 $h$ 是一个
四元数矩阵；若 $X=X^T$，则 $h^*=h$，因为 $\mathbb M(h^*)=y_c^*=y_c$ 且 $\mathbb M$ 是单射。

*计算：* 在三个例子中各取交换子的一个随机实对称元素，验证命题 4.23——[教程，4.7 节](../../tutorial/zh/04_symmetries.md#s-4-7)。

> **物理意义**（物理解释，而非本课程的定理）。若哈密顿量与一个对称性可交换，并且在使该对称性为实的基下是实的
> ——这是时间反演对称性的效果，它在那里作用为复共轭——那么命题 4.23 说明了它的块是什么样子：在实型扇区中为实对称
> 矩阵，在复型扇区对中为带有同谱伙伴的复厄米矩阵，在四元数型扇区中为能级二重简并的四元数厄米矩阵，例如奇数个
> 自旋 $\frac12$ 的时间反演，其平方为 $-1$（Kramers 定理）。Dyson 的*三重分类*（threefold way；F. J. Dyson，
> *J. Math. Phys.* 3 (1962) 1199）用相应类别的随机矩阵来模拟每一类一般哈密顿量的块——即 5.3 节的 GOE 和
> GUE，以及本课程不予定义的高斯辛系综（GSE）——并据此预言每个块的能级统计。

<a id="s-4-8"></a>

## 4.8 可与之比较的工具

*库中：* 与其他软件的比较——[教程，4.8 节](../../tutorial/zh/04_symmetries.md#s-4-8)。

<a id="x-4"></a>

## 习题

1. 证明对 $S_L$ 在 $(\mathbb C^d)^{\otimes L}$ 上的表示，$\chi_\pi(\sigma)=d^{\#\text{cycles}(\sigma)}$ 且
   $\dim\pi(S_L)'=\binom{d^2+L-1}L$。验证 $20$（2.4 节）和 $35$（$=$ 四个自旋的 $\dim\{S_a\}''$）。
2. （*Clebsch–Gordan。*）对自旋表示 $S^{(j_1)}_a\otimes1+1\otimes S^{(j_2)}_a$ 证明
   $j_1\otimes j_2\cong\bigoplus_jj$，$j=|j_1-j_2|,|j_1-j_2|+1,\dots,j_1+j_2$——这就是第 7 章提到的 $SU(2)$ 的
   非阿贝尔融合规则。
3. 设 $\theta/2\pi$ 为无理数。证明对单个 $U_\theta=e^{i\theta Q}$ 有 $\{U_\theta\}'=U(1)'$，但对 $n$ 个量子比特有
   $\{e^{i\pi\hat N}\}'\cong M_{2^{n-1}}\oplus M_{2^{n-1}}$。
4. （*Zeeman 分裂。*）设 $H$ 是 $L\ge1$ 个自旋上的 Heisenberg 哈密顿量，$H_h=H+hS_z$。证明 $H_h$ 与 $S_z$
   和 $\mathbf S^2$ 可交换，但当 $h\neq0$ 时与 $S_x$ 不可交换，并且它的特征值为 $E+hM$，$M=-S,\dots,S$，其中
   $E$ 取遍 $H_S$ 的特征值。
5. （*四元数群。*）设 $Q_8=\{\pm1,\pm I,\pm J,\pm K\}\subseteq\mathbb H$，其中 $I=(i,0)$，$J=(0,1)$，$K=(0,i)$
   （E.6 节），并设 $\pi$ 是它在 $\mathbb R^8$ 上的正则表示，$\pi(g)e_h=e_{gh}$。证明
   $\mathcal A=\mathrm{span}\,\pi(Q_8)$ 有四个 $k=m=1$ 的实型扇区和一个 $k=m=2$ 的四元数型扇区，因而
   $\mathcal A_{\mathbb R}\cong\mathbb R^4\oplus\mathbb H$。

<details>
<summary>解答</summary>

1. 由引理 B.1(c)，$\pi(\sigma)$ 是置换矩阵，它固定 $|b_1\dots b_L\rangle$ 当且仅当 $b$ 在 $\sigma$ 的各轮换上为常数，
   而这样的 $b$ 共有 $d^{\#\text{cycles}(\sigma)}$ 个。由推论 4.5(c)，
   $\dim\pi(S_L)'=\frac1{L!}\sum_\sigma(d^2)^{\#\text{cycles}(\sigma)}=\frac{x(x+1)\cdots(x+L-1)}{L!}$，其中 $x=d^2$
   （恒等式 $\sum_\sigma x^{\#\text{cycles}(\sigma)}=x(x+1)\cdots(x+L-1)$ 在第 5 章习题 2 的解答中证明），这就是
   $\binom{d^2+L-1}L$（命题 B.3 通过对多重集计数证明了同一维数）：$\binom63=20$，$\binom74=35$。$35$ 也等于
   $\dim\{S_a\}''$，这就是 $L=4$ 时的定理 4.14(b)：$\{S_a\}''=\pi(S_4)'$。

2. $S_a=S^{(j_1)}_a\otimes1+1\otimes S^{(j_2)}_a$ 是自旋表示（4.5 节），而酉算子 $x\otimes y\mapsto y\otimes x$ 把它与
   $(j_2,j_1)$ 的自旋表示缠结起来；所以可设 $j_1\ge j_2$。$S_z$ 在 $|m_1\rangle\otimes|m_2\rangle$ 上的特征值为 $m_1+m_2$，
   所以它的特征值属于 $j_1+j_2+\mathbb Z$，而对这样的 $M$，$d_M$ 是满足 $|M-m_2|\le j_1$ 的 $m_2\in\{-j_2,\dots,j_2\}$
   的个数（此时自动有 $m_1=M-m_2\in j_1+\mathbb Z$）；对 $M\notin j_1+j_2+\mathbb Z$，$d_M=0$。对 $M\ge0$，条件为
   $M-j_1\le m_2\le M+j_1$，且 $M+j_1\ge j_2$，所以 $m_2$ 从 $\max(-j_2,M-j_1)$ 取到 $j_2$：当 $0\le M\le j_1-j_2$ 时
   $d_M=2j_2+1$，当 $j_1-j_2\le M\le j_1+j_2$ 时 $d_M=j_1+j_2-M+1$，当 $M>j_1+j_2$ 时 $d_M=0$。由推论 4.12，自旋
   $j\ge0$ 的重数为 $m_j=d_j-d_{j+1}$。若 $j\notin j_1+j_2+\mathbb Z$ 或 $j>j_1+j_2$，它为 $0$；$m_{j_1+j_2}=1-0=1$；对
   $j_1-j_2\le j<j_1+j_2$（其中 $j\in j_1+j_2+\mathbb Z$，所以 $j+1\le j_1+j_2$），
   $m_j=(j_1+j_2-j+1)-(j_1+j_2-j)=1$；而对 $0\le j<j_1-j_2$（所以 $j+1\le j_1-j_2$），$m_j=(2j_2+1)-(2j_2+1)=0$。
   因此 $j_1\otimes j_2\cong\bigoplus_{j=j_1-j_2}^{j_1+j_2}j$，每个自旋恰好出现一次。对 $j_2=\frac12$，这就是
   定理 4.16(d)，后者还给出了系数。

3. $U_\theta=\sum_qe^{iq\theta}P_q$，而 $e^{iq\theta}=e^{iq'\theta}$ 当且仅当 $(q-q')\theta\in2\pi\mathbb Z$，对无理的
   $\theta/2\pi$ 这迫使 $q=q'$。所以 $U_\theta$ 的特征子空间就是诸 $P_q$ 的像，而一个算子与正规矩阵 $U_\theta$ 可交换
   当且仅当它保持其特征子空间（定理 1.5(d) 以及如 2.2 节中的 Lagrange 插值），即当且仅当 $X=\sum_qP_qXP_q$；
   由命题 4.10，这就是 $U(1)'$。论证只用到对出现的荷 $q$ 而言诸 $e^{iq\theta}$ 互不相同。
   对 $\theta=\pi$，特征值 $(-1)^w$ 只看汉明重量的奇偶性：有两个维数为 $2^{n-1}$ 的特征子空间（$n\ge1$ 时，
   即偶权重的构型与奇权重的构型；翻转第一个比特是它们之间的一个双射），交换子为 $M_{2^{n-1}}\oplus M_{2^{n-1}}$。

4. $S_z$ 和 $\mathbf S^2$ 与 $H$ 可交换（命题 4.15(a)），也与 $S_z$ 可交换，因而与 $H_h$ 可交换；并且
   $[H_h,S_x]=h[S_z,S_x]=ihS_y\neq0$。在命题 4.13 的分解中，
   $V^*S_zV=\bigoplus_SS^{(S)}_z\otimes1_{m_S}$，$V^*HV=\bigoplus_S1_{2S+1}\otimes H_S$（命题 4.15），所以在自旋 $S$
   的直和项上 $V^*H_hV=hS^{(S)}_z\otimes1+1\otimes H_S$，这是可交换算子之和，由特征向量的乘积对角化，特征值为
   $hM+E$。

5. 由定理 E.19 的证明，$I^2=J^2=K^2=-1$，且 $IJ=-JI=K$，$JK=-KJ=I$，$KI=-IK=J$，所以 $Q_8$ 是一个群，并且
   $h\mapsto gh$ 是 $Q_8$ 的双射：诸 $\pi(g)$ 是置换矩阵，构成一个由实矩阵给出的酉表示，而
   $\mathcal A=\pi(Q_8)''$ 在共轭下封闭（引理 4.18(c)）。$\chi_\pi(g)$ 是满足 $gh=h$ 的 $h$ 的个数，即 $g=1$ 时为 $8$，
   否则为 $0$，所以由推论 4.5(b)，不可约酉表示 $\rho$ 的重数为 $m_\rho=\frac18\overline{\chi_\rho(1)}\cdot8=\dim\rho$。
   *一维表示：* 忽略符号时，$I,J,K$ 的乘法与 $\mathbb Z_2\times\mathbb Z_2$ 的非单位元相同，所以对 $s,t\in\{\pm1\}$，
   $\chi(\pm1)=1$，$\chi(\pm I)=s$，$\chi(\pm J)=t$，$\chi(\pm K)=st$ 是一个同态；这是四个不同的特征标，因而给出四个
   互不等价的不可约表示，每个的 $m=1$。*二维表示：* $\rho(g)=M(g)$（引理 E.15(2)）是酉表示
   （$M(g)^*M(g)=M(\bar gg)=1$），$\chi_\rho(\pm1)=\pm2$，在其他元素上 $\chi_\rho=0$（$\mathrm{Tr}\,M(a,b)=2\mathrm{Re}\,a$），
   所以 $\langle\chi_\rho,\chi_\rho\rangle_{Q_8}=\frac{4+4}8=1$，$\rho$ 不可约（推论 4.5(c)），且 $m=2$。对如定理 2.6 中
   的分解 $V$，诸 $\pi_c\circ\pi$ 不可约且两两不等价（命题 2.7），所以 $V$ 是如定理 4.3 中的分解，而上面五个表示
   中的每一个（$m_\rho\ge1$）都恰与一个 $\pi_c\circ\pi$ 等价，不同的表示对应不同的 $c$（推论 4.5(b)）。由于
   $4\cdot1\cdot1+2\cdot2=8=\sum_ck_cm_c$，不存在其他扇区。由定理 4.20(c)，我们可以使用一个对这些表示满足
   $\pi_c(\pi(g))=\rho_c(g)$ 的分解（引理 4.2(c) 给出酉矩阵 $u$ 使 $\rho_c=u\,\pi_c(\pi(\cdot))\,u^*$；
   把 $V_c$ 换成 $V_c(u^*\otimes1)$）。由于 $\pi(g)$ 是实的，定理 4.20(a) 给出
   $\pi_{\sigma(c)}(\pi(g))=J_c\overline{\rho_c(g)}J_c^*$。对特征标 $\chi$，这就是 $\chi(g)$，所以 $\pi_{\sigma(c)}\circ\pi\cong\pi_c\circ\pi$，
   $\sigma(c)=c$（命题 2.7：不同的扇区承载不等价的表示），且 $J_c\bar J_c=|J_c|^2=1$：实型。对
   $\rho=M$，$\sigma(c)=c$（这是唯一 $k=2$ 的扇区），且 $\overline{M(g)}=\Omega_2^*M(g)\Omega_2$（由定理 4.21(b) 的证明，$M(q)\Omega_2=\Omega_2\overline{M(q)}$），所以
   $M(g)=(J_c\Omega_2^*)M(g)(J_c\Omega_2^*)^*$。矩阵 $M(1)=1$，$M(I)=\mathrm{diag}(i,-i)$，$M(J)=\Omega_2$，
   $M(K)=\begin{pmatrix}0&-i\\-i&0\end{pmatrix}$ 线性无关，所以它们张成 $M_2(\mathbb C)$；因此
   $J_c\Omega_2^*=\lambda1$，$J_c=\lambda\Omega_2$，且 $J_c\bar J_c=|\lambda|^2\Omega_2^2=-1$：四元数型。由
   定理 4.21，$\mathcal A_{\mathbb R}\cong\mathbb R^4\oplus M_1(\mathbb H)$，其维数为 $8=|Q_8|$。

</details>

*计算：* 习题 1–5 的数值检验——[教程，习题](../../tutorial/zh/04_symmetries.md#x-4)。

**下一章：** [第 5 章](05_random_matrices.md)在代数中抽取随机矩阵：$\prod_cU(k_c)$ 的 Haar 酉矩阵、Weingarten
演算——其二阶情形依赖于 $\{V\otimes V\}$ 的交换子，即 Schur–Weyl 对偶最小的实例——以及逐个扇区的能级统计。
