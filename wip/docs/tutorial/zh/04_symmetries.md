[English](../04_symmetries.md) · [Español](../es/04_symmetries.md) · [Français](../fr/04_symmetries.md) · [Deutsch](../de/04_symmetries.md) · **中文** · [日本語](../ja/04_symmetries.md) · [Русский](../ru/04_symmetries.md)

# 4. 对称性与表示

第 2 章表明，一个对称群产生一个冯·诺依曼代数，即对称算子的线性张成，以及它的交换子，即对称的
可观测量（命题 2.7），并手工分解了一个例子。本章给出其背后的表示论，把它变成本库可以作用于任意
矩阵列表的算法，并把它应用于自旋系统的两种对称性：$U(1)$（粒子数或磁化强度）和 $SU(2)$（自旋的转动）。

**你将学到**完全可约性与 Schur 引理，分解为不可约表示的方式如何决定交换子与生成的代数，特征标如何
计算重数，为什么交换子的一个随机元素能数值地找到分解，什么是荷扇区，不可约自旋表示以及 $L$ 个
自旋 $\frac12$ 如何分解，以及为什么 $SU(2)$ 不变哈密顿量的块比 $U(1)$ 不变哈密顿量的块小得多。

在本章中，$S'$ 是 $S\subseteq M_n(\mathbb C)$ 的交换子（第 2 章）。若 $S^*=S$，则 $S''$ 是**由** $S$
**生成的冯·诺依曼代数**：它是包含 $S$ 的 $*$-代数，并且每个 $*$-代数 $\mathcal B\supseteq S$ 都满足
$\mathcal B=\mathcal B''\supseteq S''$（定理 2.2）。$\mathcal A'_h=\{K\in\mathcal A':K=K^*\}$ 是一个实
向量空间。$\mathbb C^k\otimes\mathbb C^m$ 采用第 7 章的 Kronecker 顺序（$e_a\otimes f_j$ 是第 $am+j$ 个基向量，
从 $0$ 开始编号），因此 $1_k\otimes x$ 是含 $k$ 个 $x$ 拷贝的分块对角矩阵。$|G|$ 是有限群的阶，并取 $\hbar=1$。

## 4.1 表示与 Schur 引理

**定义。** 设 $X$ 是带有对合 $x\mapsto x^*$（$x^{**}=x$）的集合。$X$ 在 $\mathbb C^k$ 上的一个
**表示**是满足 $\rho(x^*)=\rho(x)^*$ 的映射 $\rho:X\to M_k(\mathbb C)$。若对所有 $x$ 有
$\rho(x)K\subseteq K$，则称子空间 $K$ 是**不变的**；若 $k\ge1$ 且 $0,\mathbb C^k$ 是仅有的不变子空间，
则称 $\rho$ 是**不可约的**。从 $\rho$ 到 $\mathbb C^{k'}$ 上的表示 $\rho'$ 的**缠结算子**是对所有 $x$ 满足
$T\rho(x)=\rho'(x)T$ 的 $k'\times k$ 矩阵 $T$；它们构成一个向量空间 $\operatorname{Hom}(\rho,\rho')$，并且
$\operatorname{Hom}(\rho,\rho)=\rho(X)'$。若存在可逆的缠结算子，则称 $\rho$ 与 $\rho'$ **等价**，记作 $\rho\cong\rho'$。

会出现三种情形：群 $G$，取 $g^*=g^{-1}$，此时要求 $\rho$ 是同态 $G\to U(k)$（**酉表示**；于是
$\rho(g^{-1})=\rho(g)^{-1}=\rho(g)^*$）；带有其伴随运算的 $*$-代数，此时 $\rho$ 是保单位的 $*$-同态；以及
满足 $x^*=x$ 的生成元标签，例如 $\{x,y,z\}$，用厄米矩阵组来表示，比如自旋矩阵（4.5 节）。若 $K$ 不变，
且 $Q\in M_{k\times d}$ 的列构成 $K$ 的一组规范正交基，则 $\rho(x)Q=QQ^*\rho(x)Q$，因此
$\rho_K(x)=Q^*\rho(x)Q$ 仍是一个表示（若 $\rho$ 是同态，则它也是同态），称为 $\rho$ 在 $K$ 上的**限制**。

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

**引理 4.2（Schur 引理）。** 设 $\rho,\rho'$ 是 $X$ 在 $\mathbb C^k,\mathbb C^{k'}$ 上的不可约表示，
$T\in\operatorname{Hom}(\rho,\rho')$。
(a) $T=0$ 或 $T$ 可逆。(b) $\operatorname{Hom}(\rho,\rho)=\mathbb C1$。
(c) $T^*\in\operatorname{Hom}(\rho',\rho)$，并且若 $T\neq0$，则 $T=\lambda U$，其中 $\lambda>0$，$U$ 为酉矩阵。
因此，当 $\rho\cong\rho'$ 时 $\dim\operatorname{Hom}(\rho,\rho')$ 为 $1$，否则为 $0$；并且等价的不可约表示是
酉等价的。

*证明。* (a) 和 (b) 已在命题 2.7 之后对群证明过，而证明只用到 $\ker T$ 在 $\rho$ 下不变
（若 $Tv=0$，则 $T\rho(x)v=\rho'(x)Tv=0$）以及 $T$ 的像在 $\rho'$ 下不变：所以
$\ker T\in\{0,\mathbb C^k\}$ 且 $\operatorname{ran}T\in\{0,\mathbb C^{k'}\}$；对 (b)，取 $\lambda$ 为 $T$ 的一个
特征值，$T-\lambda1$ 是不可逆的缠结算子，因此为 $0$。(c) 对 $T\rho(x^*)=\rho'(x^*)T$ 取伴随得
$\rho(x)T^*=T^*\rho'(x)$。于是 $T^*T\in\operatorname{Hom}(\rho,\rho)$，所以若 $T\neq0$，则 $T^*T=\lambda^21$，
$\lambda>0$（$T^*T\ge0$ 非零），而由 (a)，$U=T/\lambda$ 是方阵且 $U^*U=1$。最后，若
$T_1,T_2\in\operatorname{Hom}(\rho,\rho')$ 非零，则 $T_1^{-1}T_2\in\operatorname{Hom}(\rho,\rho)=\mathbb C1$。$\square$

若所有 $\rho(x)$ 彼此可交换（阿贝尔群），则不可约的 $\rho$ 是一维的：每个 $\rho(x)$ 都属于
$\operatorname{Hom}(\rho,\rho)=\mathbb C1$，所以每个子空间都是不变的。

## 4.2 同型分解与交换子

**定理 4.3（同型分解）。** 设 $\rho$ 是 $X$ 在 $\mathbb C^n$ 上的表示。存在 $X$ 在 $\mathbb C^{k_c}$ 上
两两不等价的不可约表示 $\rho_1,\dots,\rho_C$、满足 $\sum_ck_cm_c=n$ 的整数 $m_c\ge1$，以及酉算子
$V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n$，使得对所有 $x$ 有
$V^*\rho(x)V=\bigoplus_c\rho_c(x)\otimes1_{m_c}$。对每个这样的分解，
$$
V^*\rho(X)'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),\qquad
V^*\rho(X)''V=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},
$$
并且 $m_c=\dim\operatorname{Hom}(\rho_c,\rho)$；$X$ 的不可约表示 $\sigma$ 满足
$\operatorname{Hom}(\sigma,\rho)\neq0$ 当且仅当恰有一个 $c$ 使 $\sigma\cong\rho_c$。

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
$\operatorname{Hom}(\rho_d,\rho_c)$，而它在 $c\neq d$ 时为 $0$，在 $c=d$ 时为 $\mathbb C1$（引理 4.2）。所以
$Y\in\rho(X)'$ 当且仅当 $Y_{cj,dl}=\delta_{cd}(y_c)_{jl}1_{k_c}$，即 $V^*YV=\bigoplus_c1_{k_c}\otimes y_c$。

*双交换子。* 若 $Z$ 与 $\rho(X)'$ 可交换，则它与到各直和项上的投影可交换（取 $y_c=1$，并对
$d\neq c$ 取 $y_d=0$），所以 $V^*ZV=\bigoplus_cZ_c$，且 $Z_c$ 与 $1_{k_c}\otimes M_{m_c}$ 可交换；由引理 2.5
第 5 步（交换两个张量因子），$Z_c\in M_{k_c}\otimes1_{m_c}$。反之，这些算子都与 $\rho(X)'$ 可交换。

*重数。* 对 $T\in\operatorname{Hom}(\sigma,\rho)$，其中 $\sigma$ 不可约，同样的论证表明
$Q_{c,j}^*T\in\operatorname{Hom}(\sigma,\rho_c)$，当 $\sigma\cong\rho_c$ 时它等于 $\mathbb Ct_{c}$（$t_c$ 可逆），否则为
$0$；并且 $T=\sum_{c,j}Q_{c,j}Q_{c,j}^*T$。由于诸 $\rho_c$ 两两不等价，$\sigma$ 至多与其中一个等价，设为
$\rho_c$（否则 $T=0$），于是 $T=\sum_ja_jQ_{c,j}t_c$；反之，每个 $Q_{c,j}t_c$ 都属于
$\operatorname{Hom}(\sigma,\rho)$，而这 $m_c$ 个映射的像相互正交，所以 $\dim\operatorname{Hom}(\sigma,\rho)=m_c$。$\square$

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

**定理 4.4（Schur 正交关系）。** 设 $\rho,\rho'$ 是有限群 $G$ 在 $\mathbb C^k,\mathbb C^{k'}$ 上的不可约酉
表示。则
$$
\frac1{|G|}\sum_{g\in G}\rho(g)_{ij}\,\overline{\rho'(g)_{i'j'}}=
\begin{cases}\delta_{ii'}\delta_{jj'}/k&\text{if }\rho=\rho',\\0&\text{if }\rho\not\cong\rho'.\end{cases}
$$

*证明。* 对 $k\times k'$ 矩阵 $Y$，令 $T=\frac1{|G|}\sum_g\rho(g)Y\rho'(g)^*$。对 $h\in G$，作替换 $g\to hg$
并利用 $\rho'(hg)^*=\rho'(g)^*\rho'(h)^{-1}$，得 $\rho(h)T=T\rho'(h)$，所以 $T\in\operatorname{Hom}(\rho',\rho)$。
若 $\rho\not\cong\rho'$，则 $T=0$；若 $\rho=\rho'$，则 $T=\lambda1$，且 $k\lambda=\operatorname{Tr}T=\operatorname{Tr}Y$
（引理 4.2）。对 $Y=|j\rangle\langle j'|$，矩阵元 $T_{ii'}$ 就是左边，而 $\operatorname{Tr}Y=\delta_{jj'}$。$\square$

表示的**特征标**是 $\chi_\rho(g)=\operatorname{Tr}\rho(g)$；等价的表示有相同的特征标，因为
$\operatorname{Tr}(T\rho(g)T^{-1})=\operatorname{Tr}\rho(g)$。记
$\langle\varphi,\psi\rangle_G=\frac1{|G|}\sum_g\overline{\varphi(g)}\psi(g)$。

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
$\operatorname{Tr}(a\otimes1_m)=m\operatorname{Tr}a$，所以 $\chi_\pi=\sum_cm_c\chi_{\rho_c}$；应用 (a)。(c) 由同样的展开和 (a)，
$\langle\chi_\pi,\chi_\pi\rangle_G=\sum_cm_c^2$，由定理 4.3 它等于 $\dim\pi(G)'$。$\square$

(c) 的第二个证明正是本库对交换子进行采样的方式。**twirl**
$P(Y)=\frac1{|G|}\sum_g\pi(g)Y\pi(g)^*$ 是幂等的，其像为 $\pi(G)'$（第 2 章习题 3），并且关于
$\langle Y,Z\rangle_{HS}=\operatorname{Tr}(Y^*Z)$ 是自伴的（作替换 $g\to g^{-1}$），因此它是到 $\pi(G)'$ 上的正交
投影，而 $\dim\pi(G)'$ 就是它的迹。映射 $Y\mapsto AYB$ 的迹为 $\operatorname{Tr}A\operatorname{Tr}B$（其对角元为
$A_{ii}B_{jj}$），由此得到 $\frac1{|G|}\sum_g|\chi_\pi(g)|^2$。

我们用特征标重做 2.4 节的例子。置换矩阵的迹等于它固定的基向量的个数，而 $\pi(\sigma)$ 固定
$|b_1b_2b_3\rangle$ 当且仅当 $b$ 在 $\sigma$ 的每个轮换上为常数（引理 B.1(c)），所以
$\chi_\pi(\sigma)=2^{\#\text{cycles}(\sigma)}$。手头有 $S_3$ 的三个酉表示。设 $P_\sigma e_i=e_{\sigma(i)}$ 是 $\mathbb C^3$ 上的
置换表示（到置换矩阵中的同态）。**平凡**表示是 $\sigma\mapsto1$。**符号**表示是
$\operatorname{sgn}\sigma=\det P_\sigma$，它是到 $\{\pm1\}\subseteq U(1)$ 中的同态，因为行列式是可乘的，而置换矩阵
的行列式为 $\pm1$。直线 $\mathbb C(1,1,1)$ 是不变的，因此 $E=(1,1,1)^\perp$ 也是不变的（引理 4.1(a)）；**标准**
表示是 $P$ 在 $E$ 上的限制。在适配于 $\mathbb C(1,1,1)\oplus E$ 的规范正交基下，$P_\sigma=1\oplus\operatorname{std}(\sigma)$，
所以 $\chi_{\rm std}(\sigma)=\operatorname{Tr}P_\sigma-1$，即不动点个数减一。$S_3$ 的六个元素是恒等元、三个对换和
两个 $3$-轮换，它们分别有 $3,2,1$ 个轮换和 $3,1,0$ 个不动点；对换交换 $1_3$ 的两列，所以其符号为 $-1$，而
$3$-轮换是两个对换的乘积（$(1\,2\,3)=(1\,3)(1\,2)$），所以其符号为 $+1$。由此得到下表

| | 恒等元（$1$） | 对换（$3$） | $3$-轮换（$2$） |
|---|---|---|---|
| $\chi_{\rm triv}$ | $1$ | $1$ | $1$ |
| $\operatorname{sgn}$ | $1$ | $-1$ | $1$ |
| $\chi_{\rm std}$ | $2$ | $0$ | $-1$ |
| $\chi_\pi$ | $8$ | $4$ | $2$ |

（所以 $\operatorname{sgn}\sigma=(-1)^{3-\#\text{cycles}(\sigma)}$；所有特征标都是实的）。逐个共轭类计算
$\langle\varphi,\psi\rangle_G$：
$\langle\chi_{\rm triv},\chi_{\rm triv}\rangle_G=\frac{1+3+2}6=1$，$\langle\operatorname{sgn},\operatorname{sgn}\rangle_G=\frac{1+3+2}6=1$，
$\langle\chi_{\rm std},\chi_{\rm std}\rangle_G=\frac{4+0+2}6=1$，所以由推论 4.5(c)，这三个表示都是不可约的，并且
由于它们的特征标不同，它们两两不等价。推论 4.5(b) 给出它们在 $\pi$ 中的重数：
$$
m_{\rm triv}=\tfrac{8+3\cdot4+2\cdot2}6=4,\qquad m_{\rm sgn}=\tfrac{8-3\cdot4+2\cdot2}6=0,\qquad
m_{\rm std}=\tfrac{2\cdot8+0-2\cdot2}6=2 .
$$
所以平凡表示和标准表示分别等价于 $\pi$ 的成分 $\rho_c$（定理 4.3），其 $m_c=4$ 和 $m_c=2$；它们占据
$4\cdot1+2\cdot2=8=\sum_ck_cm_c$ 维，所以没有其他成分：$k=(1,2)$，$m=(4,2)$，并且
$\dim\pi(S_3)'=4^2+2^2=\frac16(64+3\cdot16+2\cdot4)=20$（推论 4.5(c)）。代码重新计算这张表：它把 $P_\sigma$ 求为
$\pi(\sigma)$ 在 $|100\rangle,|010\rangle,|001\rangle$ 张成的不变子空间上的限制（$\pi(\sigma)$ 在其上置换唯一的那个 $1$
的位置，引理 B.1(a)），并且 $\#\text{cycles}(\sigma)=\dim\ker(P_\sigma-1)$（在各轮换上为常数的向量）。

```python
import math
import torch
import torch_vn_algebra.decompose as dc

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    return (torch.as_tensor(X) - torch.as_tensor(Y)).abs().max().item() < tol

G = dc.permutation_representation(3, 2)                  # (6, 8, 8): pi(sigma) permutes three qubits
chi_pi = G.diagonal(dim1=-2, dim2=-1).sum(-1).real       # characters Tr pi(sigma)
P3 = G[:, [4, 2, 1]][:, :, [4, 2, 1]].real              # on |100>, |010>, |001>: S_3 permuting C^3
cycles = 3 - torch.linalg.matrix_rank(P3 - torch.eye(3, dtype=torch.float64)).double()   # dim ker(P - 1)
assert close(chi_pi, 2.0 ** cycles)
irreps = {'trivial': torch.ones(6, dtype=torch.float64), 'sign': (-1.0) ** (3 - cycles),
          'standard': P3.diagonal(dim1=-2, dim2=-1).sum(-1) - 1}                 # fixed points - 1
chis = torch.stack(list(irreps.values()))
assert close(chis @ chis.T / 6, torch.eye(3, dtype=torch.float64))    # Corollary 4.5(a): all irreducible
mult = chis @ chi_pi / 6                                               # Corollary 4.5(b)
print(dict(zip(irreps, mult.tolist())), (chi_pi ** 2).mean().item())
assert close(mult, torch.tensor([4.0, 0.0, 2.0], dtype=torch.float64))
assert dc.commutant_dimension(G) == 20                                # Corollary 4.5(c): M_4 (+) M_2
```

这就是 2.4 节的分解，现在由每个表示的一张六个数的表得到。

## 4.3 从生成元到代数：一个随机化算法

定理 4.3 描述了分解，但并没有计算它：它的证明需要不可约的不变子空间。Murota、Kanno、Kojima 和
Kojima（*Japan J. Indust. Appl. Math.* 27 (2010) 125–160；另见 Maehara and Murota, *ibid.* 263–293）
观察到，*交换子的一个随机元素*的特征子空间能找到它们。设 $\mathcal A$ 是由给定矩阵 $S_1,\dots,S_J$ 及其
伴随生成的冯·诺依曼代数，其分解为 $V^*\mathcal AV=\bigoplus_cM_{k_c}\otimes1_{m_c}$，
$V^*\mathcal A'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$（定理 2.6）；只有分析中用到 $V$。

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

**命题 4.7（一般性）。** 设 $K$ 和 $Y$ 分别是 $\mathcal A'_h$ 和 $\mathcal A'$ 中相互独立的随机元素，其分布在这些
实向量空间上关于 Lebesgue 测度有密度（引理 A.3）。则以概率一：(a) 引理 4.6 的假设成立；(b) 对 $K$ 的任意两个
像位于同一直和项中的不同谱投影 $\Pi_i,\Pi_j$，有 $\Pi_iY\Pi_j\neq0$。

*证明。* (a) 厄米矩阵 $\kappa=\bigoplus_c\kappa_c\in M_N$ 实线性地依赖于 $K$。$1,\kappa,\dots,\kappa^{N-1}$ 关于
$\langle\cdot,\cdot\rangle_{HS}$ 的 Gram 行列式 $\delta(K)=\det[\operatorname{Tr}\kappa^{i+j}]_{i,j=0}^{N-1}$ 是 $K$ 的坐标的实
多项式，它非零当且仅当这些幂线性无关。由于对 $\kappa=U\operatorname{diag}(\lambda)U^*$ 有
$p(\kappa)=U\operatorname{diag}(p(\lambda_l))U^*$，一个次数 $<N$ 的非零多项式零化 $\kappa$ 当且仅当它在
$\operatorname{spec}\kappa$ 上为零，而这当且仅当 $\kappa$ 的不同特征值少于 $N$ 个时才可能（对这些特征值取
$\prod_\mu(t-\mu)$）。所以 $\delta(K)\neq0$ 就是引理 4.6 的假设，并且 $\delta\not\equiv0$：它在
$V(\bigoplus_c1\otimes\operatorname{diag}(\text{互不相同的实数}))V^*\in\mathcal A'_h$ 处非零。因此由推论 A.7(b)，
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

**引理 4.8（压缩；对齐各拷贝）。** (a) 设 $W$ 是一个 $\mathcal A$-不变子空间，其规范正交基为
$Q\in M_{n\times d}$，$\sigma(a)=Q^*aQ$ 为限制。则 $\sigma(\mathcal A)$ 是 $*$-代数，
$\{Q^*YQ:Y\in\mathcal A'\}=\sigma(\mathcal A)'$，且 $\{Q^*KQ:K\in\mathcal A'_h\}=\sigma(\mathcal A)'_h$。$\mathcal A$ 在 $W$ 上
不可约地作用，当且仅当对每个 $Y\in\mathcal A'$ 有 $Q^*YQ\in\mathbb C1$。若不然，且 $X$ 是 $\mathcal A'$ 或
$\mathcal A'_h$ 中具有密度的随机元素，则以概率一 $Q^*XQ\notin\mathbb C1$。

(b) 设 $W_i,W_j$ 是 $\mathcal A$ 在其上不可约地作用的不变子空间，其规范正交基为
$Q_i,Q_j$，限制为 $\sigma_i(a)=Q_i^*aQ_i$，并设 $Y\in\mathcal A'$。则 $T=Q_i^*YQ_j\in\operatorname{Hom}(\sigma_j,\sigma_i)$。
若 $T\neq0$，则 $T=\lambda U$，其中 $\lambda>0$，$U$ 为酉矩阵，并且 $Q_j'=Q_jU^*$ 对所有 $a\in\mathcal A$ 满足
$Q_j'^*aQ_j'=\sigma_i(a)$。若 $W_i,W_j$ 位于不同的直和项中（$W_i\subseteq P_c\mathbb C^n$，
$W_j\subseteq P_d\mathbb C^n$，$c\neq d$），则 $T=0$。

*证明。* (a) 与定理 4.3 的证明一样（$W$ 在 $a$ 和 $a^*$ 下不变），对 $a\in\mathcal A$ 有 $aQ=Q\sigma(a)$ 且
$Q^*a=\sigma(a)Q^*$。因此 $\sigma(ab)=Q^*aQ\sigma(b)=\sigma(a)\sigma(b)$，$\sigma(a^*)=\sigma(a)^*$，
$\sigma(1)=Q^*Q=1_d$，所以 $\sigma$ 是保单位的 $*$-同态，$\sigma(\mathcal A)$ 是 $*$-代数；并且对 $Y\in\mathcal A'$，
$Q^*YQ\sigma(a)=Q^*YaQ=Q^*aYQ=\sigma(a)Q^*YQ$。反之，对 $T\in\sigma(\mathcal A)'$，令 $Y=QTQ^*$；则
$aY=Q\sigma(a)TQ^*=QT\sigma(a)Q^*=Ya$，所以 $Y\in\mathcal A'$，$Q^*YQ=T$，并且若 $T$ 厄米则 $Y$ 也厄米。若 $\sigma$
不可约，则 $\sigma(\mathcal A)'=\operatorname{Hom}(\sigma,\sigma)=\mathbb C1$（引理 4.2(b)）。若 $E\subseteq\mathbb C^d$ 是不同于
$0,\mathbb C^d$ 的 $\sigma$-不变子空间，则 $E^\perp$ 也是（引理 4.1(a)），而到 $E$ 上的正交投影是
$\sigma(\mathcal A)'_h$ 中的一个非数乘元素。在这种情形下，$\mathbb C1\cap\sigma(\mathcal A)'$ 和
$\mathbb R1\cap\sigma(\mathcal A)'_h$ 是真实子空间，它们在满射实线性映射 $X\mapsto Q^*XQ$ 下的原像是 $\mathcal A'$ 和
$\mathcal A'_h$ 的真子空间，而这些都是零测集（推论 A.7(c)）。

(b) 与 (a) 一样，$aQ_j=Q_j\sigma_j(a)$ 且 $Q_i^*a=\sigma_i(a)Q_i^*$，所以
$T\sigma_j(a)=Q_i^*YaQ_j=Q_i^*aYQ_j=\sigma_i(a)T$。引理 4.2 给出 $T=\lambda U$ 和
$U\sigma_j(a)U^*=\sigma_i(a)$，后者就是 $Q_j'^*aQ_j'$。最后，$P_c\in\mathcal A\cap\mathcal A'$，$Q_i=P_cQ_i$，
$P_cQ_j=0$，所以 $T=Q_i^*P_cYQ_j=Q_i^*YP_cQ_j=0$。$\square$

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

**命题 4.9（证书）。** 设 $\mathcal A$ 由 $S_1,\dots,S_J$ 及其伴随生成，
$\mathcal B=V(\bigoplus_cM_{k_c}\otimes1_{m_c})V^*$，其中 $V$ 为酉矩阵。若所有 $S_j\in\mathcal B$ 且
$\dim\mathcal A'=\sum_cm_c^2$，则 $\mathcal A=\mathcal B$。

*证明。* $\mathcal B$ 是包含所有 $S_j,S_j^*$ 的 $*$-代数，所以 $\mathcal A\subseteq\mathcal B$，
$\mathcal B'\subseteq\mathcal A'$。把定理 4.3 应用于 $\mathcal B$ 在 $\mathbb C^n$ 上的表示（分块映射
$b\mapsto b_c$ 是不可约的，并且互不等价，因为第 $c$ 个中心投影被其中一个映为 $1$，被其余的映为 $0$），得
$\dim\mathcal B'=\sum_cm_c^2=\dim\mathcal A'$。所以 $\mathcal B'=\mathcal A'$，
$\mathcal B=\mathcal B''=\mathcal A''=\mathcal A$。$\square$

`Decomposition.verify()` 报告生成元（以及交换子的新元素）到 $\mathcal B$（相应地 $\mathcal B'$）的相对距离，
并将 $\sum_cm_c^2$ 与 $\dim\mathcal A'$ 比较（对群由特征标得到，推论 4.5(c)；对生成元由零空间得到）。我们
对 $S_3$ 与三个量子比特运行它，然后看看内部：一个随机 $K$ 的特征子空间和一个随机 $Y$ 的缠结算子。

```python
D = dc.from_group(G)                                     # the twirl samples pi(S_3)'
print(D)
assert D.k_factors == (1, 2) and D.multiplicities == (4, 2)            # Section 2.4: k = (1, 2), m = (4, 2)
assert D.verify()['ok']
rho = D.to_algebra(G).matrix                             # (6, 2, 2, 2): the blocks rho_c(sigma), padded
assert close(rho[:, 0, 0, 0], torch.ones(6))                           # sector 1: trivial representation
assert close(rho[:, 1].diagonal(dim1=-2, dim2=-1).sum(-1), chis[2])    # sector 2: the standard one
schur = torch.einsum('gij,gkl->ijkl', rho[:, 1], rho[:, 1].conj()) / 6
assert close(schur.reshape(4, 4), torch.eye(4) / 2)                    # Theorem 4.4 with k = 2
assert close(D.from_algebra(D.to_algebra(G)), G)

Z = torch.randn(8, 8, dtype=dt)
K = dc.twirl(G, Z + Z.mH)                                # a random element of pi(S_3)'_h
w, Q = torch.linalg.eigh(K)
cuts = [0] + [i + 1 for i in range(7) if w[i + 1] - w[i] > 1e-8] + [8]
spaces = [Q[:, a:b] for a, b in zip(cuts[:-1], cuts[1:])]
assert sorted(W.shape[1] for W in spaces) == [1, 1, 1, 1, 2, 2]      # Lemma 4.6: m_c pieces of size k_c
assert all(close(G @ W, W @ (W.mH @ G @ W)) for W in spaces)           # invariant subspaces
Y = dc.twirl(G, torch.randn(8, 8, dtype=dt))
Wa, Wb = [W for W in spaces if W.shape[1] == 2]
s = torch.linalg.svdvals(Wa.mH @ Y @ Wb)                               # Lemma 4.8: a multiple of a unitary
assert s[1] > 1e-3 and abs(s[0] - s[1]) < 1e-10
assert all(close(W.mH @ Y @ Wa, 0) for W in spaces if W.shape[1] == 1)   # different summands: zero
```

主要开销在交换子：一个 twirl 样本需要 $|G|$ 次 $n\times n$ 矩阵乘积。对生成元，本库首先对角化 $\mathcal A$ 的
一个厄米元素 $R=U\operatorname{diag}(\lambda)U^*$（生成元厄米部分的一个随机实线性组合，加上两个这样的组合的对称化
乘积的一个倍数）。每个 $X\in\mathcal A'$ 都与 $R$ 可交换，因而与它的谱投影可交换（定理 1.5(d)），所以 $U^*XU$
关于 $R$ 的特征子空间（维数为 $d_\lambda$）是分块对角的；这一点被预先施加，于是零空间问题只有
$N_R=\sum_\lambda d_\lambda^2$ 个未知数，而不是 $n^2$ 个。把若干特征值归并在一起只会使块变大，所以约化后的问题有
相同的解集 $\mathcal A'$。在分解中，$V^*RV=\bigoplus_cr_c\otimes1_{m_c}$，其中 $r_c\in M_{k_c}$ 厄米；若 $\lambda$ 在
$r_c$ 中的重数为 $\mu_c(\lambda)$，则 $d_\lambda=\sum_c\mu_c(\lambda)m_c$，并且由于对整数 $\mu\ge0$ 有 $\mu^2\ge\mu$，
$$
N_R=\sum_\lambda\Big(\sum_c\mu_c(\lambda)m_c\Big)^2\ge\sum_{\lambda,c}\mu_c(\lambda)^2m_c^2\ge\sum_{\lambda,c}\mu_c(\lambda)m_c^2=\sum_ck_cm_c^2,
$$
等号成立当且仅当所有 $r_c$ 的特征值都是单重的且两两不同。对 $L=3,4,5$ 个自旋的总自旋，下面的代码发现等号
成立（*数值观察*）。即使是最小值也很大：对八个自旋 $\frac12$（$n=256$），$\sum_ck_cm_c^2=4900$，而约化后的
方程组是一个 $(n^2+N_R)\times N_R$ 矩阵，有 $3.5\cdot10^8$ 个复元素。因此对自旋链，4.6 节使用模块 `su2`，它事先
就知道答案。

```python
PLACEHOLDER
```

## 4.4 $U(1)$ 与荷扇区

**命题 4.10.** 设 $Q$ 是 $\mathbb C^n$ 上具有整数特征值 $q$（**荷**）的厄米矩阵，其谱投影为 $P_q$，
$d_q=\operatorname{rank}P_q$，并令 $U_\theta=e^{i\theta Q}$；记 $U(1)'=\{U_\theta:\theta\in\mathbb R\}'$，并用
$U(1)''$ 表示它的交换子。

1. $\theta\mapsto U_\theta$ 是 $U(1)=\mathbb R/2\pi\mathbb Z$ 的酉表示，并且 $U(1)$ 在 $\mathbb C^n$ 上的每个
   连续酉表示都具有这种形式。
2. 对 $X\in M_n$，下列条件等价：(i) 对所有 $\theta$ 有 $XU_\theta=U_\theta X$；(ii) $[X,Q]=0$；
   (iii) $X=\sum_qP_qXP_q$。因此 $U(1)'\cong\bigoplus_qM_{d_q}(\mathbb C)$（分块对角算子），且
   $U(1)''=\operatorname{span}\{P_q\}$：所有 $k_c=1$，而重数就是诸 $d_q$。
3. （*守恒荷。*）对厄米的 $H$，$e^{itH}Qe^{-itH}=Q$ 对所有 $t$ 成立当且仅当 $[H,Q]=0$，即当且仅当
   $H\in U(1)'$。此时 $\rho(t)=e^{-itH}\rho e^{itH}$ 的每个扇区权重 $\operatorname{Tr}(P_q\rho(t))$ 都是常数。

*证明。* (1) $U_\theta=\sum_qe^{iq\theta}P_q$ 是酉的，$U_\theta U_{\theta'}=U_{\theta+\theta'}$（定理 1.5(c)），且
$U_{2\pi}=1$。反之，对连续表示，$(U_\theta)_{\theta\ge0}$ 是矩阵的连续半群，所以由引理 8.5（在第 8 章证明；
证明只用到有限维）对 $\theta\ge0$ 有 $U_\theta=e^{\theta B}$，而对 $\theta<0$ 有
$U_\theta=U_{-\theta}^{-1}=e^{\theta B}$。在 $0$ 处对 $U_\theta^*U_\theta=1$ 求导（引理 8.1）得 $B^*+B=0$，所以
$B=iQ$，其中 $Q$ 厄米。指数级数 $e^{i\theta Q}$ 等于函数演算（命题 8.2(1)，取 $H=-Q$），所以
$U_{2\pi}=\sum_qe^{2\pi iq}P_q=1$ 迫使 $Q$ 的每个特征值 $q$ 满足 $e^{2\pi iq}=1$，$q\in\mathbb Z$。(2) (i)$\Rightarrow$(ii)：在
$\theta=0$ 处求导。(ii)$\Rightarrow$(iii)：$X$ 与 $P_q=f(Q)$ 可交换（定理 1.5(d)），所以当 $q\neq q'$ 时
$P_qXP_{q'}=XP_qP_{q'}=0$。(iii)$\Rightarrow$(i)：$XU_\theta=\sum_qe^{iq\theta}P_qXP_q=U_\theta X$。分块对角算子构成
$\bigoplus_qM_{d_q}\otimes1_1$（在适配于诸 $P_q$ 的像的规范正交基下），其交换子为
$\bigoplus_q1_{d_q}\otimes M_1=\operatorname{span}\{P_q\}$：与它们可交换的算子与每个 $P_q$ 可交换，所以是分块对角的，
而每个块与 $M_{d_q}$ 可交换，因而是数乘（引理 2.5 第 5 步，取 $m=1$）。
(3) $\frac d{dt}e^{itH}Qe^{-itH}=ie^{itH}[H,Q]e^{-itH}$（命题 8.2）对所有 $t$ 为零当且仅当 $[H,Q]=0$。此时
$e^{itH}$ 与 $P_q$ 可交换（定理 1.5(d)，用两次），并且
$\operatorname{Tr}(P_q\rho(t))=\operatorname{Tr}(e^{itH}P_qe^{-itH}\rho)=\operatorname{Tr}(P_q\rho)$。$\square$

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
这就是 `chains.SpinChain(n).algebra`，即 $n$ 个荷为 $0,1$ 的 $M_1\oplus M_1$ 拷贝的融合积（比特 $1$ = 自旋向上，
$\hat N=S^z+\frac n2$；而在 `decompose.total_spin` 中，自旋向上是 $\mathbb C^2$ 的*第一个*基向量，这只是给基
重新编号，并不改变扇区）。每个 XXZ 哈密顿量
$\sum_{\langle ij\rangle}J(s^x_is^x_j+s^y_is^y_j+\Delta s^z_is^z_j)$ 都保持 $\hat N$ 守恒：$s^z_is^z_j$ 在基 $|b\rangle$ 下是对角的，
而跃迁项 $s^x_is^x_j+s^y_is^y_j=\frac12(s^+_is^-_j+s^-_is^+_j)$（展开 $s^\pm=s^x\pm is^y$；交叉项相消）把 $|b\rangle$ 映为
$0$，或映为一个向上的自旋从格点 $i,j$ 中的一个移到另一个所得的构型，所以它把 $\hat N$ 的每个特征子空间映入自身。

```python
import functools
from torch_vn_algebra import SpinChain

L = 4
Sx, Sy, Sz = dc.total_spin(L)                              # (16, 16) each, site 0 most significant
N_op = Sz + L / 2 * torch.eye(2 ** L, dtype=dt)            # number of up spins
DU = dc.from_generators(N_op)                              # the algebra generated by U(1)
print(DU)
chain = SpinChain(L, boundary='periodic')                  # U(1)' = (+)_w M_binom(L, w)
assert set(DU.k_factors) == {1} and sorted(DU.multiplicities) == sorted(chain.algebra.k_factors)
assert sorted(chain.algebra.k_factors) == sorted(math.comb(L, w) for w in range(L + 1))

sx, sy, sz = dc.spin_matrices(0.5)
def two_site(a, b, i, j, n=L):                             # a on site i, b on site j of n sites
    ops = [a if k == i else b if k == j else torch.eye(2, dtype=dt) for k in range(n)]
    return functools.reduce(torch.kron, ops)

H = sum(two_site(sx, sx, i, j) + two_site(sy, sy, i, j) + 0.5 * two_site(sz, sz, i, j)
        for i, j in chain.bonds)                           # XXZ ring, Delta = 0.5
assert close(H @ N_op, N_op @ H)                           # [H, N] = 0 ...
assert close(DU.from_commutant(DU.to_commutant(H)), H)     # ... so H is block diagonal (Proposition 4.10)
Hc = chain.xxz(J=1.0, Delta=0.5)                           # the same H, built sector by sector
E_sec = torch.cat([torch.linalg.eigvalsh(Hc.matrix[0, N, :k, :k]) for N, k in enumerate(chain.algebra.k_factors)])
assert close(torch.sort(E_sec).values, torch.linalg.eigvalsh(H))
```

> **物理意义**（物理解释，而非本课程的定理）。$U(1)$ 对称性就是粒子数、荷或磁化强度的守恒。命题 4.10(3) 是这一情形下的 Noether
> 定理：荷守恒当且仅当哈密顿量是对称的，此时动力学永远不会离开一个扇区（定理 8.3）。若此外没有任何东西
> 能测量扇区之间的相干性，这些扇区就是超选择扇区（第 3 章）。

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
其特征值的模为 $1$、乘积为 $1$，所以 $V=W\operatorname{diag}(e^{i\varphi},e^{-i\varphi})W^*$，其中 $W$ 为酉矩阵，并且对无迹
厄米矩阵 $A=W\operatorname{diag}(\varphi,-\varphi)W^*$ 有 $V=e^{iA}$。无迹厄米 $2\times2$ 矩阵是 Pauli 矩阵的实线性组合
$A=2\,\mathbf c\cdot\mathbf s$，所以 $V=e^{-i\theta\,\mathbf n\cdot\mathbf s}$，其中 $\theta=2|\mathbf c|$，
$\mathbf n=-\mathbf c/|\mathbf c|$（若 $\mathbf c=0$ 则 $\theta=0$）。定理 B.6(a),(c) 与推论 B.7（取 $d=2$）直接表明，
诸 $V^{\otimes L}$（$V\in SU(2)$）张成由诸 $S_a$ 生成的代数。

令 $S_\pm=S_x\pm iS_y$，于是 $S_-=S_+^*$，并令 $\mathbf S^2=S_x^2+S_y^2+S_z^2$。由这些关系得
$$
[S_z,S_\pm]=\pm S_\pm,\qquad[S_+,S_-]=2S_z,\qquad\mathbf S^2=S_-S_++S_z^2+S_z,
$$
并且 $\mathbf S^2$ 与每个 $S_a$ 可交换（例如
$[\mathbf S^2,S_z]=S_x[S_x,S_z]+[S_x,S_z]S_x+S_y[S_y,S_z]+[S_y,S_z]S_y=-i\{S_x,S_y\}+i\{S_y,S_x\}=0$）。

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

**推论 4.12.** 每个 $\mathbb C^n$ 上的自旋表示都酉等价于 $\bigoplus_s(\text{spin }s)\otimes1_{m_s}$，并且记
$d_M=\dim\ker(S_z-M)$，有
$$
m_s=d_s-d_{s+1}.
$$

*证明。* 定理 4.3 与定理 4.11。若 $|M|\le s'$ 且 $s'-M\in\mathbb Z$，则自旋 $s'$ 的每个拷贝恰含 $S_z$ 属于特征值
$M$ 的一个特征向量，否则不含；所以 $d_M=\sum_{s'\ge|M|,\,s'-M\in\mathbb Z}m_{s'}$，而两个相邻值之差就是 $m_s$。$\square$

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
维数分别为 $35$ 和 $14$——正是第 10 章中得到的数。本库仅由三个生成元就得到同样的结果（扇区按
$(k_c,m_c)$ 排序）：

```python
S4 = dc.total_spin(4)
assert close(S4[0] @ S4[1] - S4[1] @ S4[0], 1j * S4[2])                # [S_x, S_y] = i S_z
D4 = dc.from_generators(S4)
print(D4)
assert D4.k_factors == (1, 3, 5) and D4.multiplicities == (2, 3, 1)     # M_1, M_3, M_5 with m = 2, 3, 1
assert D4.verify()['ok'] and D4.dims == (35, 14)
m_S = {S: math.comb(4, 2 - S) - (math.comb(4, 1 - S) if S < 2 else 0) for S in (0, 1, 2)}
assert [m_S[(k - 1) // 2] for k in D4.k_factors] == list(D4.multiplicities)    # Proposition 4.13
casimir = D4.to_algebra((S4 @ S4).sum(0))                                  # S^2 = S(S+1) on spin S
for c, k in enumerate(D4.k_factors):
    assert close(casimir.matrix[0, c, :k, :k], (k - 1) / 2 * (k + 1) / 2 * torch.eye(k))
```

哪些算子与所有转动可交换？格点置换就是：用 $\pi(\sigma)$ 共轭把 $s^{(i)}_a$ 映为 $s^{(\sigma(i))}_a$（引理 B.1(b)），所以它
保持 $\sum_is^{(i)}_a$ 不变。反之，除此之外再没有别的：

**定理 4.14（Schur–Weyl 对偶）。** 设 $\pi$ 是 $S_L$ 置换 $(\mathbb C^2)^{\otimes L}$ 各张量因子的表示。
(a) $\operatorname{span}\pi(S_L)\subseteq\{S_x,S_y,S_z\}'$。(b) 等号成立；等价地（定理 2.2），由总自旋生成的
代数是 $\pi(S_L)'$。

*证明。* (a) 部分刚刚已证。(b) 部分就是推论 B.7，即定理 B.6 在 $d=2$ 时的情形，而定理 B.6 对所有 $d$ 和 $L$
证明了 $(\mathbb C^d)^{\otimes L}$ 与 $U(d)$ 的对角作用的对偶。$\square$

（延伸阅读：W. Fulton and J. Harris, *Representation Theory*, Springer 1991, Lecture 6；R. Goodman and
N. Wallach, *Symmetry, Representations, and Invariants*, Springer 2009, Section 9.1。）对 $L\le5$，代码数值地比较
这两个空间：由 (a)，只需它们维数相同。

```python
for n in range(2, 6):
    comm = dc.commutant_basis(dc.total_spin(n)).reshape(-1, 4 ** n)          # a basis of {S_a}'
    perm = dc.permutation_representation(n, 2).reshape(-1, 4 ** n)          # all n! permutations
    dims = [torch.linalg.matrix_rank(A).item() for A in (comm, perm, torch.cat([comm, perm]))]
    print(n, dims)
    assert dims[0] == dims[1] == dims[2]                                     # span pi(S_n) = {S_a}'
Dp = dc.from_group(dc.permutation_representation(4, 2))                   # pi(S_4)'' = {S_a}'
print(Dp)
assert sorted(zip(Dp.multiplicities, Dp.k_factors)) == sorted(zip(D4.k_factors, D4.multiplicities))
```

打印出的维数是 $2,5,14,42$，即 Catalan 数 $\binom{2L}L/(L+1)$；推论 B.10 对每个 $L$ 证明了
$\dim\{S_a\}'=\sum_Sm_S^2=\binom{2L}L/(L+1)$。最后一行展示了对偶的作用：$\pi(S_4)''$
的块为 $k=(1,2,3)$，重数为 $(5,1,3)$——正是总自旋的数对 $(k_c,m_c)$ 互换角色。一般地（推论 B.9），自旋 $S$ 的
重数空间承载 $S_L$ 的一个维数为 $m_S$ 的不可约表示，并且这些表示两两不等价；这里对 $S=2,1,0$ 其维数分别为 $1,3,2$。

> **物理意义**（物理解释，而非本课程的定理）。对每个转动不变的哈密顿量，一个自旋 $S$ 多重态由 $2S+1$ 个能量相同的态组成，它们由
> $S^z=M$ 区分。两个自旋 $\frac12$ 构成一个三重态和一个单态；四个自旋构成一个五重态、三个三重态和两个
> 单态。两个单态的重数空间和三个三重态的重数空间就是第 10 章中无退相干的量子比特与无噪声的三能级系统（qutrit）。

## 4.6 自旋系统

$L$ 个自旋 $\frac12$ 上的 **Heisenberg 模型**是 $H=\sum_{i<j}J_{ij}\,\mathbf s^{(i)}\cdot\mathbf s^{(j)}$，其中耦合
$J_{ij}$ 为实数。

**命题 4.15.** (a)（*Dirac 恒等式*）$\mathbf s^{(i)}\cdot\mathbf s^{(j)}=\frac12F_{ij}-\frac14$，其中 $F_{ij}$
交换自旋 $i$ 与 $j$。因此每个 Heisenberg 哈密顿量都属于 $\operatorname{span}\pi(S_L)\subseteq\{S_a\}'$，并且在命题
4.13 的分解中 $V^*HV=\bigoplus_S1_{2S+1}\otimes H_S$，其中块 $H_S\in M_{m_S}$：$H_S$ 的每个特征值都是 $H$ 的一个
$(2S+1)$ 重简并能级。
(b) **最高权向量**空间 $\ker S_+\cap\ker(S_z-S)$ 的维数为 $m_S$，它在 $H$ 下不变，且 $H$ 在其上作用为
$H_S$（相差一个酉基变换）。

*证明。* (a) 在两个自旋上，$\mathbf s^{(1)}\cdot\mathbf s^{(2)}=\frac12(\mathbf S^2-\mathbf s^{(1)2}-\mathbf s^{(2)2})=\frac12\mathbf S^2-\frac34$
（定理 4.11，$s=\frac12$）。在三重态上 $\mathbf S^2=2$ 且 $F=1$（三重态由 $S_-^j|{\uparrow\uparrow}\rangle$ 张成，而 $F$
与 $S_-$ 可交换并固定 $|{\uparrow\uparrow}\rangle$）；在单态上 $\mathbf S^2=0$ 且 $F=-1$。两边在这四个向量上一致。对
$L$ 个格点中的格点 $i,j$，与单位算子作张量积即可。其余部分是定理 4.14(a) 与定理 4.3。(b) 在推论 4.12 的分解中，
满足 $S_z=S$ 的向量为 $\sum_{S'\ge S}|S',S\rangle\otimes w_{S'}$，而由定理 4.11，当 $S'>S$ 时 $S_+|S',S\rangle\neq0$；所以
$\ker S_+\cap\ker(S_z-S)=|S,S\rangle\otimes\mathbb C^{m_S}$，$1\otimes H_S$ 在其上作用为 $H_S$。$\square$

(b) 部分就是模块 `su2`（自旋链，第 9 章）不经数值分解而构造 $H_S$ 的方式：在 `SpinChain` 的 $U(1)$ 扇区
$\hat N=\frac L2+S$ 内部，通过逐个格点的 Clebsch–Gordan 耦合构造最高权向量的一组规范正交基。我们用 `D4`
在四个自旋的环上检验 (a)，然后比较八个自旋时的块大小。

```python
def swap(i, j, n):
    p = list(range(n))
    p[i], p[j] = j, i
    return dc.permutation_representation(n, 2, [p])[0]

ring = [(i, (i + 1) % 4) for i in range(4)]
H4 = sum(swap(i, j, 4) / 2 - torch.eye(16, dtype=dt) / 4 for i, j in ring)        # Dirac's identity
assert close(H4, sum(two_site(a, a, i, j) for a in (sx, sy, sz) for i, j in ring))
h = D4.to_commutant(H4)                                   # blocks H_S of sizes m_S = 2, 3, 1
assert close(D4.from_commutant(h), H4)                    # H4 lies in {S_a}'
levels = torch.cat([torch.linalg.eigvalsh(h.matrix[0, c, :m, :m]).repeat_interleave(k)
                    for c, (k, m) in enumerate(zip(D4.k_factors, D4.multiplicities))])
assert close(torch.sort(levels).values, torch.linalg.eigvalsh(H4))   # level of H_S repeated 2S+1 times

from torch_vn_algebra import su2
u1, sym = SpinChain(8, boundary='periodic'), su2.SU2Chain(8, boundary='periodic')
H_u1, H_su2 = u1.xxz(J=1.0, Delta=1.0), sym.heisenberg(J=1.0)
E_u1 = torch.sort(torch.cat([torch.linalg.eigvalsh(H_u1.matrix[0, N, :k, :k])
                             for N, k in enumerate(u1.algebra.k_factors)])).values
E_su2, S_lev, _ = sym.full_spectrum(H_su2, return_quantum_numbers=True)   # level of block S, 2S+1 times
E_su2, S_lev = E_su2[0], S_lev[0]
assert close(E_u1, E_su2)
assert S_lev[0] == 0 and (S_lev[1:4] == 1).all() and E_su2[4] - E_su2[3] > 1e-3     # singlet, then triplet
assert close(sym.restrict(H_u1).matrix, H_su2.matrix)      # Proposition 4.15(b)
flops = lambda ks: sum(k ** 3 for k in ks)
print(u1.algebra.k_factors, sym.multiplicities, flops(u1.algebra.k_factors), flops(sym.multiplicities))
print(f"ground state energy {E_su2[0].item():.6f}, first excited {E_su2[1].item():.6f} (x3)")
```

对八个自旋的环，$U(1)$ 块的大小为 $1,8,28,56,70,56,28,8,1$，$SU(2)$ 块的大小为 $m_S=14,28,20,7,1$
（$S=0,\dots,4$）。对角化的开销 $\sum k^3$ 从 $739162$ 降到 $33040$，减小为约 $1/22$；代数的维数 $\sum k^2$
从 $\binom{16}8=12870$ 降到 $1430$。基态是能量为 $-3.651093$ 的单态，第一激发能级是位于 $-3.128419$ 的三重态。

> **物理意义。** Heisenberg 模型描述磁性绝缘体中局域电子自旋的交换相互作用；$J_{ij}>0$ 倾向于使相邻
> 自旋反平行。它的谱由多重态组成，而磁场 $h\,S^z$ 把每个多重态分裂为 $2S+1$ 个等间距能级（习题 4）。
> 若二分格子的两个子格之间的耦合 $J_{ij}>0$，子格内部的耦合 $J_{ij}\le0$，且两个子格大小相等，则基态是
> 单态（Lieb–Mattis 定理，E. Lieb and D. Mattis, *J. Math. Phys.* 3 (1962) 749；此处不加证明地陈述）。
> 上面八个自旋的环就是一个例子。

## 4.7 可与之比较的工具

对有限群，GAP 和 SageMath 精确计算特征标表，然后由推论 4.5 得到重数；RepLAB（MATLAB；D. Rosset,
F. Montealegre-Mora, J.-D. Bancal）数值地分解有限群和紧群的表示，同样使用交换子的随机元素；
SymbolicWedderburn.jl（Julia）为半定规划的对称性约化计算 Wedderburn 分解。在物理中，QuSpin 和其他精确
对角化程序对 $U(1)$ 和晶格对称性进行块对角化，ITensor 和 TeNPy 在张量网络中携带阿贝尔荷，TensorKit.jl
用融合树处理 $SU(2)$。本库的 `decompose` 是面向小矩阵（$n$ 至多几百）的稠密方法；它接受任意矩阵列表
（群元素或生成元），验证其结果（命题 4.9），并把代数及其交换子作为带有 `to_algebra`、`from_algebra`、
`to_commutant`、`from_commutant` 的 `TypeIAlgebra` 对象返回，使第 3 章、第 5–8 章的内容可以直接应用。对
自旋 $\frac12$ 链，`su2.SU2Chain` 直接给出 $SU(2)$ 块，对 $L$ 至多 $16$ 的所有自旋值都适用。

## 习题

1. 证明对 $S_L$ 在 $(\mathbb C^d)^{\otimes L}$ 上的表示，$\chi_\pi(\sigma)=d^{\#\text{cycles}(\sigma)}$ 且
   $\dim\pi(S_L)'=\binom{d^2+L-1}L$。验证 $20$（2.4 节）和 $35$（$=$ 四个自旋的 $\dim\{S_a\}''$）。
2. （*Clebsch–Gordan。*）对自旋表示 $S^{(j_1)}_a\otimes1+1\otimes S^{(j_2)}_a$ 证明
   $j_1\otimes j_2\cong\bigoplus_jj$，$j=|j_1-j_2|,|j_1-j_2|+1,\dots,j_1+j_2$——这就是第 7 章提到的 $SU(2)$ 的
   非阿贝尔融合规则——并用 `from_generators` 验证 $1\otimes\frac12$ 和 $1\otimes1$。
3. 设 $\theta/2\pi$ 为无理数。证明对单个 $U_\theta=e^{i\theta Q}$ 有 $\{U_\theta\}'=U(1)'$，但对 $n$ 个量子比特有
   $\{e^{i\pi\hat N}\}'\cong M_{2^{n-1}}\oplus M_{2^{n-1}}$。
4. （*Zeeman 分裂。*）设 $H$ 是 $L\ge1$ 个自旋上的 Heisenberg 哈密顿量，$H_h=H+hS_z$。证明 $H_h$ 与 $S_z$
   和 $\mathbf S^2$ 可交换，但当 $h\neq0$ 时与 $S_x$ 不可交换，并且它的特征值为 $E+hM$，$M=-S,\dots,S$，其中
   $E$ 取遍 $H_S$ 的特征值。对八个自旋的环验证这一点。

<details>
<summary>解答</summary>

1. 与 $d=2$ 时一样：$\pi(\sigma)$ 是置换矩阵，它固定 $|b_1\dots b_L\rangle$ 当且仅当 $b$ 在各轮换上为常数，而这样的
   $b$ 共有 $d^{\#\text{cycles}}$ 个。由推论 4.5(c)，
   $\dim\pi(S_L)'=\frac1{L!}\sum_\sigma(d^2)^{\#\text{cycles}(\sigma)}=\frac{x(x+1)\cdots(x+L-1)}{L!}$，其中 $x=d^2$
   （第 5 章习题 2），这就是 $\binom{d^2+L-1}L$：$\binom63=20$，$\binom74=35$。$35$ 也等于
   $\dim\{S_a\}''$，这就是 $L=4$ 时的定理 4.14。

```python
for n, value in ((3, 20), (4, 35)):
    assert dc.commutant_dimension(dc.permutation_representation(n, 2)) == math.comb(n + 3, n) == value
```

2. 设 $j_1\ge j_2$。$S_z$ 在 $|m_1\rangle\otimes|m_2\rangle$ 上的特征值为 $m_1+m_2$，所以对
   $M\in j_1+j_2+\mathbb Z$，$d_M$ 是满足 $|M-m_2|\le j_1$ 的 $m_2\in\{-j_2,\dots,j_2\}$ 的个数：若 $|M|\le j_1-j_2$，
   则全部 $2j_2+1$ 个都满足；若 $j_1-j_2\le|M|\le j_1+j_2$，则为 $j_1+j_2-|M|+1$ 个（对 $M\ge0$，条件为
   $m_2\ge M-j_1$）；超出此范围则为 $0$。由推论 4.12，$m_j=d_j-d_{j+1}$ 在 $j_1-j_2\le j\le j_1+j_2$ 时为 $1$，
   在 $j<j_1-j_2$ 时为 $(2j_2+1)-(2j_2+1)=0$。

```python
def coupled(a, b):                                        # S^(a) (x) 1 + 1 (x) S^(b)
    Sa, Sb = dc.spin_matrices(a), dc.spin_matrices(b)
    ea, eb = torch.eye(Sa.shape[-1], dtype=dt), torch.eye(Sb.shape[-1], dtype=dt)
    return torch.stack([torch.kron(x, eb) + torch.kron(ea, y) for x, y in zip(Sa, Sb)])
for a, b, ks in ((1.0, 0.5, (2, 4)), (1.0, 1.0, (1, 3, 5))):  # 1 (x) 1/2 = 1/2 (+) 3/2, 1 (x) 1 = 0 (+) 1 (+) 2
    D_cg = dc.from_generators(coupled(a, b))
    assert D_cg.k_factors == ks and set(D_cg.multiplicities) == {1}
```

3. $U_\theta=\sum_qe^{iq\theta}P_q$，而 $e^{iq\theta}=e^{iq'\theta}$ 当且仅当 $(q-q')\theta\in2\pi\mathbb Z$，对无理的
   $\theta/2\pi$ 这迫使 $q=q'$。所以 $U_\theta$ 的特征子空间就是诸 $P_q$ 的像，而一个算子与正规矩阵 $U_\theta$ 可交换
   当且仅当它保持其特征子空间（定理 1.5(d) 以及如 2.2 节中的 Lagrange 插值），即当且仅当 $X=\sum_qP_qXP_q$；
   由命题 4.10，这就是 $U(1)'$。
   对 $\theta=\pi$，特征值 $(-1)^w$ 只看汉明重量的奇偶性：有两个维数为 $2^{n-1}$ 的特征子空间（$n\ge1$ 时，
   即偶子集与奇子集），交换子为 $M_{2^{n-1}}\oplus M_{2^{n-1}}$。

```python
for theta, m in ((1.0, (1, 1, 4, 4, 6)), (math.pi, (8, 8))):
    D_theta = dc.from_generators(torch.linalg.matrix_exp(1j * theta * N_op))
    assert D_theta.multiplicities == m and set(D_theta.k_factors) == {1}
```

4. $S_z$ 和 $\mathbf S^2$ 与 $H$ 可交换（命题 4.15(a)），也与 $S_z$ 可交换，因而与 $H_h$ 可交换；并且
   $[H_h,S_x]=h[S_z,S_x]=ihS_y\neq0$。在命题 4.13 的分解中，
   $V^*S_zV=\bigoplus_SS^{(S)}_z\otimes1_{m_S}$，$V^*HV=\bigoplus_S1_{2S+1}\otimes H_S$（命题 4.15），所以在自旋 $S$
   的直和项上 $V^*H_hV=hS^{(S)}_z\otimes1+1\otimes H_S$，这是可交换算子之和，由特征向量的乘积对角化，特征值为
   $hM+E$。

```python
assert not close(H4 @ Sx - Sx @ H4 + 0.3 * (Sz @ Sx - Sx @ Sz), 0)          # [H4 + 0.3 S_z, S_x] != 0
E_field = sym.full_spectrum(sym.heisenberg(J=1.0, h=0.3))[0]                # E + h M on every multiplet
H_field = u1.xxz(J=1.0, Delta=1.0, h=0.3)
E_u1_field = torch.cat([torch.linalg.eigvalsh(H_field.matrix[0, N, :k, :k])
                        for N, k in enumerate(u1.algebra.k_factors)])
assert close(torch.sort(E_u1_field).values, E_field)
```

</details>

**下一章：** [第 5 章](05_random_matrices.md)在代数中抽取随机矩阵：$\prod_cU(k_c)$ 的 Haar 酉矩阵、Weingarten
演算——其二阶情形依赖于 $\{V\otimes V\}$ 的交换子，即 Schur–Weyl 对偶最小的实例——以及逐个扇区的能级统计。
