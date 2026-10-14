[English](../08_dynamics.md) · [Español](../es/08_dynamics.md) · [Français](../fr/08_dynamics.md) · [Deutsch](../de/08_dynamics.md) · **中文** · [日本語](../ja/08_dynamics.md) · [Русский](../ru/08_dynamics.md)

# 8. 动力学

在量子力学中，封闭系统用一个单参数酉群来建模；开放系统——与一个不被观测的环境耦合的系统——则用一个
信道半群来建模。在有扇区的代数中，这两种演化都具有一个在单个矩阵代数上看不到的特点：
$\mathcal M$ 中的哈密顿量永远不能改变扇区，所以荷是守恒的；而开放系统可以在扇区之间跳跃，
例如由于丢失一个粒子。本章证明关于这两种演化的基本事实，写出带有扇区间跳跃的 Lindblad 方程，
并说明"没有粒子丢失"是对荷的一次测量，其后选择的态按照一个非厄米哈密顿量演化。

**你将学到**为什么扇区权重是运动常数，为什么 Gorini–Kossakowski–Sudarshan–Lindblad（GKSL）
形式的生成元产生信道，扇区间的跳跃如何进入 Lindblad 方程，以及扇区权重何时服从经典的速率方程，
无跳跃演化与后选择意味着什么，以及如何求稳态与弛豫速率。本章的数值检验见
[教程 8](../../tutorial/zh/08_dynamics.md)。

在本章中，$\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ 作用在
$\mathcal H=\bigoplus_c\mathbb C^{k_c}$ 上，$D=\sum_ck_c$，$P_c$ 为极小中心投影，
态由其 blunt 密度 $\rho$ 给出，扇区权重为 $p_c=\mathrm{Tr}\,\rho_c$
（第 3 章），$[A,B]=AB-BA$，$\{A,B\}=AB+BA$，并取 $\hbar=1$。**超算子**是 $\mathcal M$ 或
$M_D(\mathbb C)$ 上的线性映射；超算子的范数一律取算子范数，因而是次乘性的。对有限维空间上的
线性映射 $B$，$e^{tB}=\sum_{n\ge0}t^nB^n/n!$。

<a id="s-8-1"></a>

## 8.1 封闭系统

**哈密顿量**是一个自伴元 $H\in\mathcal M$。向量按照 **Schrödinger
方程** $i\psi'(t)=H\psi(t)$ 演化，密度按照**冯·诺依曼方程** $\rho'(t)=-i[H,\rho(t)]$ 演化。
二者都是有限维空间上的常系数线性方程，下面的引理就是我们所需要的全部微分方程理论。

<a id="r-8-1"></a>

**引理 8.1.** 设 $B$ 是有限维赋范空间 $X$ 上的线性映射。则
$\frac d{dt}e^{tB}=Be^{tB}=e^{tB}B$，$e^{(s+t)B}=e^{sB}e^{tB}$，$\|e^{tB}\|\le e^{|t|\,\|B\|}$，并且对
每个 $x_0\in X$，映射 $x(t)=e^{tB}x_0$ 是区间
$I\ni0$ 上满足 $x'=Bx$ 与 $x(0)=x_0$ 的唯一可微函数 $x$。

*证明。* 该级数绝对收敛，并满足上述范数界；作为 $t$ 的幂级数，它可以
逐项求导，这就给出第一个公式。于是
$\frac d{dt}(e^{tB}e^{-tB})=Be^{tB}e^{-tB}-e^{tB}Be^{-tB}=0$，所以 $e^{tB}e^{-tB}=1$。若 $x'=Bx$，则
$\frac d{dt}(e^{-tB}x(t))=e^{-tB}(x'-Bx)=0$，所以 $e^{-tB}x(t)=x_0$，从而 $x(t)=e^{tB}x_0$。最后，
$t\mapsto e^{(s+t)B}x_0$ 与 $t\mapsto e^{tB}e^{sB}x_0$ 都是 $y'=By$ 的解，且初值相同。$\square$

<a id="r-8-2"></a>

**命题 8.2.** 设 $H=H^*\in\mathcal M$，$U(t)=e^{-itH}$，它由 1.3 节的函数演算
以 $f_t(\lambda)=e^{-it\lambda}$ 定义。

1. $U(t)=\bigoplus_ce^{-itH_c}$ 是 $\mathcal M$ 中的酉元，$U(s+t)=U(s)U(t)$，且 $U'(t)=-iHU(t)$；
   $U(t)$ 等于 $-itH$ 的指数级数。
2. $\psi(t)=U(t)\psi_0$ 是 $i\psi'=H\psi$ 满足 $\psi(0)=\psi_0$ 的唯一解，而
   $\rho(t)=U(t)\rho_0U(t)^*$ 是 $\rho'=-i[H,\rho]$ 在 $\mathcal M$ 中满足
   $\rho(0)=\rho_0$ 的唯一解。若 $\rho_0$ 是态，则每个 $\rho(t)$ 也是态。
3. （*Heisenberg 绘景*）$\mathrm{Tr}(\rho(t)A)=\mathrm{Tr}(\rho_0A(t))$，其中
   $A(t)=U(t)^*AU(t)$，且 $A'=i[H,A]$。

*证明。* (1) 对每个块对角化，$H_c=V_c\Lambda_cV_c^*$；则 $V=\bigoplus_cV_c$ 使
$H$ 对角化，且 $U(t)=Ve^{-it\Lambda}V^*$ 是块对角的（由定理 1.5(a)，结果不依赖于
所用的对角化酉矩阵）。由于 $|f_t|^2=1$ 且 $f_sf_t=f_{s+t}$，定理 1.5(c) 给出酉性与
群律。对 $e^{-it\Lambda}$ 逐元素求导得 $U'=V(-i\Lambda)e^{-it\Lambda}V^*=-iHU$，所以由引理 8.1
（$X=M_D$，$B$ 为左乘 $-iH$），$U(t)$ 就是级数 $e^{-itH}$。
(2) $(U\psi_0)'=-iHU\psi_0$，且 $(U\rho_0U^*)'=-iHU\rho_0U^*+U\rho_0U^*(iH)$，因为 $(U')^*=iU^*H$。
唯一性即引理 8.1，其中分别取 $\mathcal H$ 上的 $B=-iH$ 与 $\mathcal M$ 上的 $B=-i[H,\cdot\,]$。正性：
$\langle x,U\rho_0U^*x\rangle=\langle U^*x,\rho_0U^*x\rangle\ge0$；迹由循环性得到。(3) 迹的
循环性，以及 $(U^*AU)'=iU^*HAU-iU^*AHU=i[H,A(t)]$。$\square$

*库中：* 逐块计算的 $U(t)=e^{-itH}$ 是 `dynamics.propagator`；命题 8.2(2) 中的解是 `dynamics.schrodinger` 和 `dynamics.von_neumann`——[教程，8.1 节](../../tutorial/zh/08_dynamics.md#s-8-1)。

<a id="r-8-3"></a>

**定理 8.3（荷守恒）。** 设 $t\mapsto H(t)$ 在区间 $I$ 上取值于 $\mathcal M$ 的
自伴元。

1. 若 $\psi:I\to\mathcal H$ 可微且 $i\psi'=H(t)\psi$，则每个 $\|P_c\psi(t)\|$ 都是常数。
2. 若 $\rho:I\to\mathcal M$ 可微且 $\rho'=-i[H(t),\rho]$，则对每个中心元 $Z$，$\mathrm{Tr}(Z\rho(t))$
   都是常数；特别地，每个扇区权重 $p_c(t)=\mathrm{Tr}(P_c\rho(t))$ 都是常数。
3. 若 $H$ 是常数且 $A\in\mathcal M$ 与 $H$ 可交换，则 $\mathrm{Tr}(A\rho(t))$ 是常数。

*证明。* (1) 内积关于第一个变元是共轭线性的，所以
$\frac d{dt}\langle\psi,P_c\psi\rangle=\langle-iH\psi,P_c\psi\rangle+\langle\psi,-iP_cH\psi\rangle
=i\langle\psi,HP_c\psi\rangle-i\langle\psi,P_cH\psi\rangle=0$，因为 $P_c$ 是中心元。(2)
由循环性，$\frac d{dt}\mathrm{Tr}(Z\rho)=-i\mathrm{Tr}(Z[H,\rho])=-i\mathrm{Tr}([Z,H]\rho)=0$。
(3) 对 $A$ 作同样的计算。$\square$

> **物理意义**（物理解释，而非本课程的定理）。$\mathcal M$ 中的哈密顿量与荷（粒子数、
> 磁化、对称性量子数）可交换，所以每个荷值的概率都是运动常数，
> 即使对受驱动的、含时的哈密顿量也是如此。

在 $\mathcal M$ 上，单位向量 $\psi=\bigoplus_c\psi_c$ 定义的是混合态
$\bigoplus_c|\psi_c\rangle\langle\psi_c|$ 的态：每个 $A\in\mathcal M$ 都是块对角的，所以
$\langle\psi,A\psi\rangle=\sum_c\langle\psi_c,A_c\psi_c\rangle=\mathrm{Tr}\big(A\bigoplus_c|\psi_c\rangle\langle\psi_c|\big)$
（与第 3 章习题 3 相同）。若 $\psi(t)$ 是 $i\psi'=H\psi$ 的解，则由于 $U(t)$ 是块对角的，
$\bigoplus_c|\psi_c(t)\rangle\langle\psi_c(t)|=U(t)\big(\bigoplus_c|\psi_c(0)\rangle\langle\psi_c(0)|\big)U(t)^*$。

*计算：* 定理 8.3，以及 $\psi$ 与 $\bigoplus_c|\psi_c\rangle\langle\psi_c|$ 的演化相一致，分别用精确求解器和用 Runge–Kutta 求解器 `dynamics.schrodinger_rk4`（对受驱动的哈密顿量）检验——[教程，8.1 节](../../tutorial/zh/08_dynamics.md#s-8-1)。

<a id="s-8-2"></a>

## 8.2 量子动力学半群

CP 映射、TP 映射以及信道在 6.1 节（[第 6 章](06_quantum_channels.md)）中定义。除
定理 6.2 和命题 6.4（映射 $A\mapsto\sum_lK_lAK_l^*$——包括扇区之间的映射——是 CP 的）之外，我们还需要
三条封闭性质。

<a id="r-8-4"></a>

**引理 8.4.** (a) CP 映射的和、非负倍数与复合都是 CP 的。设 $\Phi$ 是
$\mathcal A=M_D$ 或 $\mathcal A=\mathcal M$ 上的 CP 映射，并设 $\mathcal N\subseteq\mathcal A$ 是 $\mathcal M$（当
$\mathcal A=M_D$ 时）或者某个角 $P_c\mathcal MP_c\cong M_{k_c}$（以 $P_c$ 为单位）。若 $\Phi$ 把
$\mathcal N$ 映入自身，则 $\Phi|_{\mathcal N}$ 是 CP 的。
(b) 若 CP 映射 $\Phi_m$ 收敛于 $\Phi$（即对所有 $A$ 有 $\Phi_m(A)\to\Phi(A)$），则 $\Phi$ 是 CP 的。
(c) 若 $s\mapsto\Psi(s)$ 是从 $[a,b]$ 到 CP 映射的连续映射，则 $\int_a^b\Psi(s)\,ds$ 是 CP 的。

*证明。* (a) $\Phi\mapsto\Phi_n$ 是线性的，且 $(\Phi\Psi)_n=\Phi_n\Psi_n$；正矩阵之和是
正的。关于限制：$M_n(\mathcal N)\subseteq M_n(\mathcal A)$，而正性在两者中是同一个
条件。事实上，由 6.1 节，$M_n(\mathcal M)$ 的元素 $[A_{ij}]$ 是正的当且仅当每个
块矩阵 $[(A_{ij})_c]$ 是正的；对 $\mathcal A=M_D$，这就是 $[A_{ij}]\in M_{nD}$ 的正性，因为经过
基的置换，$[A_{ij}]$ 是这些块的直和（所有 $A_{ij}$ 都是块对角的）；
而 $M_n(P_c\mathcal MP_c)$ 的元素只有唯一一个非零块 $[(A_{ij})_c]$，它的正性就是
在 $M_n(M_{k_c})$ 中的正性。(b) 正矩阵构成闭集（对每个 $x$，$\langle x,Xx\rangle\ge0$
是闭条件），且 $(\Phi_m)_n(X)\to\Phi_n(X)$。(c) 由 (a)，Riemann 和是 CP 的，并且
它们收敛于该积分（引理 A.4(a)）；应用 (b)。$\square$

**定义。** $\mathcal M$ 上的**量子动力学半群**是一族信道
$\mathcal M\to\mathcal M$，$(T_t)_{t\ge0}$，满足 $T_0=\mathrm{id}$、$T_{s+t}=T_sT_t$，且 $t\mapsto T_t$ 连续。

<a id="r-8-5"></a>

**引理 8.5.** 每个连续的超算子半群 $(T_t)_{t\ge0}$（$T_0=\mathrm{id}$，$T_{s+t}=T_sT_t$）都可以写成
$T_t=e^{t\mathcal L}$，其中超算子 $\mathcal L$ 是唯一确定的，称为它的**生成元**，并且
$\mathcal L=\lim_{t\to0}(T_t-\mathrm{id})/t$。对任意有限维空间上（例如矩阵上）的线性映射构成的连续半群，
同样的结论成立，证明也相同。

*证明。*（取值为超算子的连续函数的积分就是引理 A.4 中的积分。）
对 $V_h=\frac1h\int_0^hT_s\,ds$ 有 $\|V_h-\mathrm{id}\|\le\max_{0\le s\le h}\|T_s-\mathrm{id}\|\to0$（引理 A.4(c)），所以
由行列式的连续性，对某个固定的小 $h$ 有 $\det V_h\neq0$。于是
$T_tV_h=\frac1h\int_0^hT_{t+s}\,ds=\frac1h\int_t^{t+h}T_s\,ds$（引理 A.4(b)、(e)）关于 $t$ 可微，导数为
$\frac1h(T_{t+h}-T_t)=T_t\frac{T_h-\mathrm{id}}h$。因此 $T_t$ 可微，且 $T_t'=T_t\mathcal L$，
其中 $\mathcal L=\frac{T_h-\mathrm{id}}hV_h^{-1}$，再由引理 8.1（取 $B$ 为
右复合 $\mathcal L$）得 $T_t=e^{t\mathcal L}$。唯一性：$\mathcal L=T'(0)$。$\square$

<a id="r-8-6"></a>

**定理 8.6（GKSL 生成元）。** 设 $H=H^*$ 与 $L_1,\dots,L_J$ 属于 $M_D$，并且
$$
\mathcal L(\rho)=-i[H,\rho]+\sum_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
$$
则 $(e^{t\mathcal L})_{t\ge0}$ 是 $M_D$ 上的量子动力学半群。更确切地，令
$H_{\mathrm{eff}}=H-\frac i2\sum_jL_j^*L_j$，$\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)$，
$\Phi(\rho)=\sum_jL_j\rho L_j^*$，于是 $\mathcal L=\mathcal L_0+\Phi$，并定义
$$
S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},\qquad
S_n(t)=\int_0^tS_0(t-s)\,\Phi\,S_{n-1}(s)\,ds\quad(n\ge1).
$$
则每个 $S_n(t)$ 都是 CP 的，且 $e^{t\mathcal L}=\sum_{n\ge0}S_n(t)$（**Dyson 展开**），并有
$\|S_n(t)\|\le e^{t\|\mathcal L_0\|}(t\|\Phi\|)^n/n!$。

*证明。* *第 1 步。* $X(t)=e^{-itH_{\mathrm{eff}}}\rho e^{itH_{\mathrm{eff}}^*}$ 满足 $X'=\mathcal L_0(X)$
与 $X(0)=\rho$，所以由引理 8.1，$S_0(t)=e^{t\mathcal L_0}$；它是 CP 的，只有一个 Kraus 算子
$e^{-itH_{\mathrm{eff}}}$。*第 2 步。* 记 $a=\|\mathcal L_0\|$，$b=\|\Phi\|$。由归纳法，各 $S_n$
连续，且
$\|S_n(t)\|\le\int_0^te^{(t-s)a}\,b\,e^{sa}\frac{(bs)^{n-1}}{(n-1)!}ds=e^{ta}\frac{(bt)^n}{n!}$，所以
级数在每个 $[0,T]$ 上一致收敛。*第 3 步。* 由于
$S_n(t)=S_0(t)\int_0^te^{-s\mathcal L_0}\Phi S_{n-1}(s)\,ds$，乘积法则与微积分基本
定理（引理 A.4(e)）给出 $S_n'=\mathcal L_0S_n+\Phi S_{n-1}$（$n\ge1$）以及 $S_0'=\mathcal L_0S_0$。由这些
界，导数的级数在 $[0,T]$ 上一致收敛，所以 $Y(t)=\sum_nS_n(t)$
可微，$Y'=\mathcal L_0Y+\Phi Y=\mathcal LY$，$Y(0)=\mathrm{id}$，再由引理 8.1 得 $Y(t)=e^{t\mathcal L}$。
*第 4 步。* $\Phi$ 是 CP 的（定理 6.2），所以由归纳法，被积函数
$s\mapsto S_0(t-s)\Phi S_{n-1}(s)$ 是一族连续的 CP 映射，从而 $S_n(t)$ 是 CP 的（引理 8.4）；
部分和及其极限 $e^{t\mathcal L}$ 也是如此。*第 5 步。* 由循环性，
$\mathrm{Tr}\,\mathcal L(\rho)=\sum_j\big(\mathrm{Tr}(L_j\rho L_j^*)-\mathrm{Tr}(L_j^*L_j\rho)\big)=0$，
因此 $\frac d{dt}\mathrm{Tr}\,e^{t\mathcal L}\rho=\mathrm{Tr}\,\mathcal L(e^{t\mathcal L}\rho)=0$：每个
$e^{t\mathcal L}$ 都是 TP 的。半群律与连续性来自引理 8.1。$\square$

*计算：* 对 $M_3$ 上任意的跳跃检验定理 8.6，包括跳跃部分 $\sum_{n\ge1}S_n(t)$ 的完全正性——[教程，8.2 节](../../tutorial/zh/08_dynamics.md#s-8-2)。*库中：* $e^{t\mathcal L}$ 是 `dynamics.lindblad_channel`，其矩阵 $\mathcal L$ 是 `dynamics.lindblad_superoperator`——[教程，8.2 节](../../tutorial/zh/08_dynamics.md#s-8-2)。

<a id="r-8-7"></a>

**定理 8.7（Gorini–Kossakowski–Sudarshan；Lindblad）。** (1) 反之，$M_D$ 上每个量子动力学半群的
生成元都具有定理 8.6 的形式，且跳跃个数 $J\le D^2-1$。可以选取
$\mathrm{Tr}\,H=0$、$\mathrm{Tr}\,L_j=0$，并且对 $j\neq k$ 有 $\mathrm{Tr}(L_j^*L_k)=0$；在这些
归一化之下，$H$ 与 $\Phi=\sum_jL_j\cdot L_j^*$ 由 $\mathcal L$ 唯一确定，而各 $L_j$ 在相差一个
$L'_j=\sum_ku_{jk}L_k$（$u$ 为酉矩阵；较短的一族先用零补齐）的意义下唯一确定。

(2) $\mathcal M=\bigoplus_cM_{k_c}$ 上每个量子动力学半群的生成元都是
$$
\mathcal L(\rho)=\Phi(\rho)-\tfrac12\{\Phi^*(1),\rho\}-i[H,\rho]
$$
其中 $\Phi:\mathcal M\to\mathcal M$ 是一个 CP 映射（扇区之间的 Kraus 算子，命题 6.4），且
$H=H^*\in\mathcal M$：这就是 8.3 节的形式，其中只有一个扇区之间的跳跃 $\Phi_1=\Phi$，速率为 $1$。

*证明。* (1) 是定理 D.4 连同定理 D.5(b)、(c) 中的唯一性陈述，(2) 是定理 D.7，
二者都在[附录 D](D_gksl.md)中证明。$\square$

关键的一步（定理 D.3）是：CP 映射半群的生成元的 Choi 矩阵在最大纠缠向量的正交补上是正的；
仅有半群性质而各 $T_t$ 只是正映射，是不够的（见定理 D.3 之后的例子）。对于 (2)，映射 $T_t\circ E$
（其中 $E(X)=\sum_cP_cXP_c$ 是条件期望）是 $M_D$ 上的 CP 映射；它们在 $t=0$ 处的导数
$\mathcal L\circ E$ 的 Choi 矩阵在各向量 $(1\otimes P_c)\Omega$ 的正交补上是正的，这给出
$\mathcal L\circ E$ 的一个分解（引理 D.2），再用 $E$ 把它压缩回 $\mathcal M$。
在 Heisenberg 绘景中（推论 D.8），$\mathcal M$ 上由保单位 CP 映射构成的连续半群的生成元为 $A\mapsto\Psi(A)+K^*A+AK$，其中 $\Psi$ 是 CP 的，$K\in\mathcal M$。历史上，(1)
归功于 V. Gorini、A. Kossakowski、E. C. G. Sudarshan [*J. Math. Phys.* 17, 821 (1976)] 和 G. Lindblad
[*Commun. Math. Phys.* 48, 119 (1976)]，而 Heisenberg 绘景中的形式（对一般冯·诺依曼代数上的范数连续半群）
归功于 E. Christensen 和 D. E. Evans [*J. London Math. Soc.* 20, 358 (1979)]。

> **物理意义**（物理解释，而非本课程的定理）。一个与大环境弱相互作用、而环境的关联
> 衰减很快（Born–Markov 近似）的系统，预期会丧失对其过去的记忆：这时它在时刻 $s+t$ 的态
> 由时刻 $s$ 的态经过一个只依赖于 $t$ 的信道得到。*如果*情况确实如此，定理 8.7 表明，
> 该动力学是一个哈密顿部分加上"跳跃" $L_j$（光子发射、粒子丢失、
> 退相位），其速率由 $L_j$ 的归一化确定。

<a id="s-8-3"></a>

## 8.3 有扇区的代数中的 Lindblad 方程

在 $\mathcal M$ 上有两种跳跃。**扇区内部**的跳跃是一个算子
$L_j\in\mathcal M$。**扇区之间**的跳跃是一个 CP 映射 $\Phi:\mathcal M\to\mathcal M$，其 Kraus
算子为 $K^{dc}_l:\mathbb C^{k_c}\to\mathbb C^{k_d}$（命题 6.4，即一个 `InterSectorChannel`）：
$\Phi(\rho)_d=\sum_{c,l}K^{dc}_l\rho_cK^{dc*}_l$。由命题 6.8，它的 blunt 对偶满足
$\Phi^*(1)_c=\sum_{d,l}K^{dc*}_lK^{dc}_l$，这是 $\mathcal M$ 的一个块对角元。取速率
$g_j,g_m>0$，`lindblad_rhs` 的生成元为
$$
\mathcal L(\rho)=-i[H,\rho]+\sum_jg_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big)
+\sum_mg_m\Big(\Phi_m(\rho)-\tfrac12\{\Phi_m^*(1),\rho\}\Big).
$$
由下面的命题 8.8(1)，它生成 $\mathcal M$ 上的一个量子动力学半群；反之，由
定理 8.7(2)（定理 D.7），$\mathcal M$ 上每个量子动力学半群的生成元都具有这种形式。

*库中：* 扇区之间的跳跃是一个 `InterSectorChannel`，而 $\mathcal L(\rho)$ 是 `dynamics.lindblad_rhs`——[教程，8.3 节](../../tutorial/zh/08_dynamics.md#s-8-3)。

<a id="r-8-8"></a>

**命题 8.8.** 设 $J^{dc}_l\in M_D$ 是在位置 $(d,c)$ 处为块 $K^{dc}_l$、
其余处为零的矩阵。

1. $\mathcal L$ 把 $\mathcal M$ 映入 $\mathcal M$，并且是 $M_D$ 上以 $\sqrt{g_j}L_j$ 和 $\sqrt{g_m}J^{dc}_{m,l}$（所有 $m,l,d,c$）为跳跃的 GKSL 生成元
   $\hat{\mathcal L}$ 在 $\mathcal M$ 上的限制。因此
   $(e^{t\mathcal L})_{t\ge0}$ 是 $\mathcal M$ 上的量子动力学半群。
2. 另设 $\tilde{\mathcal L}$ 是 $M_D$ 上的 GKSL 生成元，其中同一个 Kraus
   指标的各块被加成一个跳跃，$L_{m,l}=\sum_{d,c}J^{dc}_{m,l}$（哈密顿量相同，跳跃 $\sqrt{g_j}L_j$
   也相同，另有跳跃 $\sqrt{g_m}L_{m,l}$）。则对 $\rho\in\mathcal M$ 有
   $\mathcal L(\rho)=E\big(\tilde{\mathcal L}(\rho)\big)$，其中 $E(X)=\sum_dP_dXP_d$ 是
   $M_D$ 到 $\mathcal M$ 上的条件期望；一般而言 $\tilde{\mathcal L}(\rho)\notin\mathcal M$。
3. 若每个 $L_{m,l}$ 在每个块行和每个块列中至多有一个非零块（例如
   它只把扇区 $c$ 映入扇区 $c-1$，就像下面例子中的湮灭算子那样），
   则对 $\rho\in\mathcal M$ 有 $\tilde{\mathcal L}(\rho)=\mathcal L(\rho)$。

*证明。* 由于 $\rho=\sum_cP_c\rho P_c$ 且 $J^{dc}P_{c'}=\delta_{cc'}J^{dc}$，
$J^{dc}\rho J^{dc*}$ 就是位于位置 $(d,d)$ 的块 $K^{dc}\rho_cK^{dc*}$，而 $J^{dc*}J^{dc}$ 是位于位置 $(c,c)$ 的
$K^{dc*}K^{dc}$。对 $d,l$ 求和得到 $\Phi(\rho)$ 与 $\Phi^*(1)$，这就证明了 (1) 中的
恒等式，并且所有项都属于 $\mathcal M$。所以 $\hat{\mathcal L}$ 及其所有幂都把 $\mathcal M$
映入 $\mathcal M$，$e^{t\hat{\mathcal L}}|_{\mathcal M}=e^{t\mathcal L}$，而由
定理 8.6 和引理 8.4(a)，它是 $\mathcal M$ 的一个信道。(2) 对 $L=\sum_{d,c}J^{dc}$ 有
$L\rho L^*=\sum_{d,d',c}J^{dc}\rho J^{d'c*}$，而
$E$ 保留 $d=d'$ 的项。对 $B\in\mathcal M$，$P_dXBP_d=P_dXP_dB$，所以 $E(XB)=E(X)B$，
$E(BX)=BE(X)$；于是 $E\{L^*L,\rho\}=\{E(L^*L),\rho\}$，并且 $E(L^*L)=\sum_cP_cL^*LP_c=\sum_{c,d}J^{dc*}J^{dc}$，
因为当 $d\ne d'$ 时 $J^{dc*}J^{d'c}=0$。最后 $E(-i[H,\rho])=-i[H,\rho]$，并且跳跃
$\sqrt{g_j}L_j\in\mathcal M$ 的各项都属于 $\mathcal M$。(3) 设 $L=\sum_{(d,c)\in S}J^{dc}$，其中指标对的集合 $S$
对给定的 $d$ 至多含一个对，对给定的 $c$ 也至多含一个对。对 $\rho\in\mathcal M$，
$J^{dc}\rho J^{d'c'*}=J^{dc}P_c\rho P_{c'}J^{d'c'*}$ 除非 $c=c'$ 否则为零，而 $c=c'$ 时又有 $d=d'$；同样
$J^{dc*}J^{d'c'}$ 除非 $d=d'$ 否则为零，而 $d=d'$ 时又有 $c=c'$。所以 $L\rho L^*$ 和 $L^*L$ 都是块对角的，
$\tilde{\mathcal L}(\rho)$ 的所有项都属于 $\mathcal M$，再由 (2) 得 $\tilde{\mathcal L}(\rho)=E(\tilde{\mathcal L}(\rho))=\mathcal L(\rho)$。$\square$

> **物理意义**（物理解释，而非本课程的定理）。第 (2) 部分表明，$\mathcal M$ 上的 Lindblad 方程是
> $M_D$ 上的方程在跳跃所产生的扇区间相干性被立即抹去之后剩下的部分——
> 这正是超选择定则的内容。对于粒子丢失，没有任何东西被抹去（第 (3) 部分）：以湮灭算子为跳跃的
> 通常的 Lindblad 方程本身已经是 $\mathcal M$ 上的方程。

*计算：* 对扇区之间的一个随机映射检验命题 8.8(2)，对粒子丢失检验命题 8.8(3)，与 $M_6$ 上的 GKSL 生成元进行比较——[教程，8.3 节](../../tutorial/zh/08_dynamics.md#s-8-3)。

<a id="r-8-9"></a>

**命题 8.9（扇区权重）。** 设 $\rho(t)=e^{t\mathcal L}\rho_0$，$p(t)$ 为其扇区权重。

1. 若没有扇区之间的跳跃（对所有 $m,l$ 以及 $d\ne c$ 有 $K^{dc}_{m,l}=0$），则每个 $p_c(t)$ 都是常数。
2. 若对每个 $m$ 以及每一对 $d\ne c$，有 $\sum_lK^{dc*}_{m,l}K^{dc}_{m,l}=r_{m,dc}1_{k_c}$，其中
   $r_{m,dc}$ 为数，则 $p'=Qp$，其中对 $d\ne c$ 有 $Q_{dc}=\sum_mg_mr_{m,dc}$，且 $Q_{cc}=-\sum_{d\ne c}Q_{dc}$。
   因此 $p(t)=e^{tQ}p(0)$ 与 $H$、各 $L_j$ 以及扇区内部的态都无关。

*证明。* 由引理 8.1，$p_d'=\mathrm{Tr}(P_d\mathcal L(\rho))$。哈密顿项如定理 8.3 中那样为零。
对 $L_j\in\mathcal M$，$\mathrm{Tr}(P_dL_j\rho L_j^*)-\mathrm{Tr}(P_dL_j^*L_j\rho)
=\mathrm{Tr}(L_j^*P_dL_j\rho)-\mathrm{Tr}(L_j^*L_jP_d\rho)=0$，因为 $P_d$ 是中心元。对扇区之间的
跳跃 $\Phi$，$\mathrm{Tr}(P_d\Phi(\rho))=\sum_c\mathrm{Tr}\big(\sum_lK^{dc*}_lK^{dc}_l\rho_c\big)$，
且 $\mathrm{Tr}(P_d\tfrac12\{\Phi^*(1),\rho\})=\sum_{d'}\mathrm{Tr}\big(\sum_lK^{d'd*}_lK^{d'd}_l\rho_d\big)$；
$c=d$ 与 $d'=d$ 的项相互抵消。这证明了 (1)；在 (2) 的假设下，剩下的项
为 $\sum_{c\ne d}r_{dc}p_c-\sum_{d'\ne d}r_{d'd}p_d$。由引理 8.1 得 $p(t)=e^{tQ}p(0)$。$\square$

在 (2) 的假设下，对 $d\ne c$，$r_{m,dc}$ 对每个满足 $p_c>0$ 的态 $\rho$ 都等于 $\Phi_m$ 的转移
矩阵（在命题 6.12 之前定义）的元素 $T_{dc}(\rho)$，因为
$\mathrm{Tr}\,\Phi_{m,dc}(\rho_c)=\mathrm{Tr}\big(\sum_lK^{dc*}_{m,l}K^{dc}_{m,l}\rho_c\big)=r_{m,dc}p_c$。所以，对速率为
$1$ 的单个跳跃，$Q=R-\mathrm{diag}(\sum_dR_{dc})_c$，其中 $R$ = `transition_matrix()`（$R$ 的对角元相互抵消）。

*库中：* $R$ 是 `InterSectorChannel.transition_matrix`——[教程，8.3 节](../../tutorial/zh/08_dynamics.md#s-8-3)。

**例：粒子丢失。** 两个玻色模式，至多两个粒子；扇区为
$N=0,1,2$，基分别为 $|00\rangle$；$|10\rangle,|01\rangle$；$|20\rangle,|11\rangle,|02\rangle$，所以
$\mathcal M=M_1\oplus M_2\oplus M_3$。Bose–Hubbard 二聚体
$H=J(a_1^*a_2+a_2^*a_1)+\frac U2\sum_jn_j(n_j-1)$ 与 $N$ 可交换，并且每个模式以
速率 $\gamma$ 丢失粒子：这是扇区之间的一个跳跃，其 Kraus 算子为 $a_1,a_2$，限制在 $N\to N-1$ 上。
由于 $a_1^*a_1+a_2^*a_2=N$，命题 8.9(2) 适用，$Q_{N-1,N}=\gamma N$，并且从 $N=2$ 出发
$$
p_2(t)=q^2,\qquad p_1(t)=2q(1-q),\qquad p_0(t)=(1-q)^2,\qquad q=e^{-\gamma t},
$$
这是参数为 $2$ 和 $q$ 的二项分布（若每个粒子都以概率 $q$ 独立地存活，这就是存活粒子数的分布），
无论隧穿与相互作用如何：两边都是 $p'=Qp$（求导即可验证）满足 $p(0)=(0,0,1)$ 的解，所以由引理 8.1
它们相等。

*计算：* 这个例子，检验 $\Phi^*(1)=0\oplus1\oplus2$，并用两个公式检验 $p(t)$——[教程，8.3 节](../../tutorial/zh/08_dynamics.md#s-8-3)。

<a id="s-8-4"></a>

## 8.4 后选择与无跳跃演化

把 8.3 节生成元的跳跃汇集为 $\Phi_{\mathrm{tot}}(\rho)=\sum_jg_jL_j\rho L_j^*+\sum_mg_m\Phi_m(\rho)$，
令 $\mathcal M$ 中的 $A=\Phi_{\mathrm{tot}}^*(1)=\sum_jg_jL_j^*L_j+\sum_mg_m\Phi_m^*(1)\ge0$，并令
$H_{\mathrm{eff}}=H-\frac i2A$。把定理 8.6 应用于命题 8.8(1) 中 $M_D$ 上的生成元 $\hat{\mathcal L}$：
由命题 8.8 的证明，它的跳跃部分把 $\mathcal M$ 映入 $\mathcal M$，并在 $\mathcal M$ 上等于 $\Phi_{\mathrm{tot}}$，
而它对所有跳跃的和 $\sum\hat L^*\hat L$ 等于 $A$，所以它的 $H_{\mathrm{eff}}$ 就是 $H-\frac i2A\in\mathcal M$。因此它的 $S_n(t)$
把 $\mathcal M$ 映入 $\mathcal M$（对 $n$ 归纳），并且对 $\rho\in\mathcal M$ 有
$$
\rho(t)=\sum_{n\ge0}S_n(t)\rho,\qquad S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},
\qquad S_n(t)=\int_0^tS_0(t-s)\,\Phi_{\mathrm{tot}}\,S_{n-1}(s)\,ds .
$$

<a id="r-8-10"></a>

**命题 8.10.** 对态 $\rho$，令 $w_n(t)=\mathrm{Tr}\,S_n(t)\rho$。则 $w_n(t)\ge0$，
$\sum_nw_n(t)=1$，且 $w_0$ 不增，$w_0'(t)=-\mathrm{Tr}\big(A\,S_0(t)\rho\big)$。对
纯态 $\rho=|\psi\rangle\langle\psi|$，$w_0(t)=\|e^{-itH_{\mathrm{eff}}}\psi\|^2$。

*证明。* $S_n(t)$ 是 CP 的，因而是正的，且 $\sum_nw_n=\mathrm{Tr}\,\rho(t)=1$。令
$X=S_0(t)\rho\ge0$，则 $w_0'=\mathrm{Tr}(-iH_{\mathrm{eff}}X+iXH_{\mathrm{eff}}^*)=-\mathrm{Tr}(AX)$，
且 $\mathrm{Tr}(AX)=\mathrm{Tr}(X^{1/2}AX^{1/2})\ge0$。最后一个结论就是
$\mathrm{Tr}|\phi\rangle\langle\phi|=\|\phi\|^2$。$\square$

> **物理意义**（物理解释，而非本课程的定理）。设一个探测器记录每一次跳跃——每个发射的光子、每个丢失的原子。
> 则 $w_n(t)$ 是到时刻 $t$ 为止恰好有 $n$ 次计数的概率，$S_n(t)\rho/w_n(t)$ 是
> 给定 $n$ 次计数时的态，而 $S_0(t)\rho/w_0(t)$ 是给定**没有计数**时的态：零结果同样会
> 改变态，这种改变通过非厄米的 $H_{\mathrm{eff}}$ 实现。这种计数过程诠释 [E. B. Davies,
> *Quantum Theory of Open Systems*, 1976; M. D. Srinivas and E. B. Davies, *Optica Acta* 28, 981 (1981)]
> 是量子轨迹（Monte Carlo 波函数）方法的基础 [J. Dalibard, Y. Castin,
> K. Mølmer, *Phys. Rev. Lett.* 68, 580 (1992)]；这里我们只用到它的数学部分。

对于损失，Dyson 项就是 $\rho(t)$ 的各扇区分量（用上面方框的语言说：计数
就等同于测量荷）。态 $\rho$ 关于满足 $\mathrm{Tr}(P\rho)>0$ 的投影 $P$ 的 **Lüders 更新**是态
$P\rho P/\mathrm{Tr}(P\rho P)$（注意 $\mathrm{Tr}(P\rho P)=\mathrm{Tr}(P\rho)$）。

*库中：* Lüders 更新是 `states.lueders_update`——[教程，8.4 节](../../tutorial/zh/08_dynamics.md#s-8-4)。

<a id="r-8-11"></a>

**命题 8.11.** 设扇区标记为 $N=0,\dots,C-1$，扇区内部没有跳跃，
每个跳跃都把荷降低一（除非 $d=c-1$，否则 $K^{dc}_{m,l}=0$），并且
$\rho=P_{N_0}\rho P_{N_0}$。则对 $0\le n\le N_0$ 有 $S_n(t)\rho=P_{N_0-n}\,\rho(t)\,P_{N_0-n}$，
对 $n>N_0$ 有 $S_n(t)\rho=0$。特别地，$p_{N_0}(t)=w_0(t)>0$，并且 $\rho(t)$ 关于 $P_{N_0}$ 的
Lüders 更新（"没有粒子丢失"）就是无跳跃态：
$$
\frac{P_{N_0}\rho(t)P_{N_0}}{p_{N_0}(t)}=\frac{e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*}}{w_0(t)},
\qquad p_{N_0}(t)=w_0(t).
$$

*证明。* 由于 $H_{\mathrm{eff}}\in\mathcal M$，$S_0(s)$ 把支撑在扇区 $N$ 中的算子映为支撑在扇区 $N$ 中的算子；
$\Phi_{\mathrm{tot}}$ 把扇区 $N\ge1$ 映到扇区 $N-1$，并把扇区 $0$ 映为零。
由归纳法，对所有 $s$，$S_n(s)\rho$ 支撑在扇区 $N_0-n$ 中，而当 $n>N_0$ 时为 $0$。用 $P_{N_0-n}$ 压缩
$\rho(t)=\sum_nS_n(t)\rho$ 恰好保留第 $n$ 项。最后 $w_0(t)>0$，因为
$e^{-itH_{\mathrm{eff}}}$ 可逆（引理 8.1），所以 $S_0(t)\rho$ 是一个非零的正矩阵。$\square$

*计算：* 对以不同速率泄漏的两个模式，检验 $n=0$ 与 $n=1$ 时的命题 8.11——[教程，8.4 节](../../tutorial/zh/08_dynamics.md#s-8-4)。

<a id="s-8-5"></a>

## 8.5 稳态与弛豫

若对所有 $t\ge0$ 有 $e^{t\mathcal L}\sigma=\sigma$，则称态 $\sigma$ 是**稳态**，等价地
$\mathcal L(\sigma)=0$（在 $t=0$ 处求导；反之利用级数）。

<a id="r-8-12"></a>

**命题 8.12.** $\mathcal M$ 上的每个量子动力学半群都有稳态。

*证明。* 固定一个态 $\rho$，令 $\sigma_\tau=\frac1\tau\int_0^\tau T_s\rho\,ds$。由引理 A.4(d)，它是一个态
（态构成闭凸集）。由引理 A.4(b)、(e)，
$T_t\sigma_\tau-\sigma_\tau=\frac1\tau\big(\int_\tau^{\tau+t}-\int_0^t\big)T_s\rho\,ds$；由于 $\|T_s\rho\|_1=1$，
由引理 A.4(c)，它的迹范数至多为 $2t/\tau$。态在有限维
空间中构成有界闭集，所以某个序列 $\sigma_{\tau_m}$，$\tau_m\to\infty$，收敛于一个态 $\sigma$
（Bolzano–Weierstrass），并且 $T_t\sigma-\sigma=\lim_m(T_t\sigma_{\tau_m}-\sigma_{\tau_m})=0$。$\square$

<a id="r-8-13"></a>

**推论 8.13.** 若没有扇区之间的跳跃（如命题 8.9(1) 中那样），则稳态张成的空间
维数至少为 $C$：每个扇区都有它自己的稳态。

*证明。* $\mathcal L$ 把 $P_c\mathcal MP_c\cong M_{k_c}$ 映入自身（$H$ 与各 $L_j$ 属于 $\mathcal M$，并且
$\Phi_m(X)=\sum_lK^{cc}_{m,l}X_cK^{cc*}_{m,l}$ 与 $\Phi_m^*(1)$ 都是块对角的），因此每个 $e^{t\mathcal L}$ 也是如此。
它的限制是 CP 的（引理 8.4(a)）且是 TP 的，所以是 $M_{k_c}$ 上的量子动力学半群，
由命题 8.12，它有稳态 $\sigma_c$；作为 $\mathcal M$ 的元素，$\sigma_c$ 是 $e^{t\mathcal L}$ 的稳态。
各 $\sigma_c$ 位于不同的块中，因而线性无关。$\square$

定理 D.6（A. Frigerio 的一个判据的有限维情形）给出唯一性的一个充分条件：若 $M_k$ 上以
$(H;L_1,\dots,L_r)$ 为数据的 GKSL 生成元有一个忠实的稳态 $\sigma>0$，并且与 $H$ 以及所有 $L_j$、$L_j^*$
可交换的矩阵只有 $1$ 的倍数，则 $\sigma$ 是唯一的稳态。没有扇区之间的跳跃时，$\mathcal L$ 在
$P_c\mathcal MP_c\cong M_{k_c}$ 上的限制是以 $(H_c;\sqrt{g_1}L_{1,c},\sqrt{g_2}L_{2,c},\dots)$ 为数据的 GKSL 生成元，
所以该定理可以逐扇区地应用。习题 4 表明，
即使在因子中，唯一性也可能不成立：那里第二个假设不成立。

*计算：* 对 $M_2\oplus M_3$ 扇区内部的随机跳跃检验推论 8.13 与定理 D.6 的假设——[教程，8.5 节](../../tutorial/zh/08_dynamics.md#s-8-5)。

<a id="r-8-14"></a>

**命题 8.14（自发辐射）。** 在以 $|g\rangle,|e\rangle$ 为基的 $M_2$ 上，令
$H=\omega|e\rangle\langle e|$，单个跳跃 $L=|g\rangle\langle e|$，速率 $\gamma>0$。则
$$
\rho_{ee}(t)=e^{-\gamma t}\rho_{ee}(0),\qquad\rho_{ge}(t)=e^{(i\omega-\gamma/2)t}\rho_{ge}(0),
\qquad e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}},
$$
其中 $\mathcal U_t(\rho)=e^{-itH}\rho e^{itH}$，$\mathcal A_p$ 是参数为
$p$ 的振幅阻尼信道（6.4 节）。$\mathcal L$ 的谱为 $\{0,-\gamma,-\gamma/2\pm i\omega\}$，
$|g\rangle\langle g|$ 是唯一的稳态，并且
$\tfrac12\|\rho(t)-|g\rangle\langle g|\|_1=\big(\rho_{ee}(0)^2e^{-2\gamma t}+|\rho_{ge}(0)|^2e^{-\gamma t}\big)^{1/2}$。

*证明。* 这里 $L^*L=|e\rangle\langle e|$，$L\rho L^*=\rho_{ee}|g\rangle\langle g|$，所以逐元素地有
$\rho_{ee}'=-\gamma\rho_{ee}$，$\rho_{gg}'=\gamma\rho_{ee}$，$\rho_{ge}'=(i\omega-\gamma/2)\rho_{ge}$，这就
给出了解。在基 $(\rho_{gg},\rho_{ge},\rho_{eg},\rho_{ee})$ 下，$\mathcal L$ 的矩阵是上
三角的，对角线为 $(0,-\gamma/2+i\omega,-\gamma/2-i\omega,-\gamma)$。$\mathcal L(\sigma)=0$ 迫使
$\sigma_{ee}=\sigma_{ge}=0$。振幅阻尼信道把 $\rho_{ee}$ 乘以 $1-p$，把 $\rho_{ge}$ 乘以
$\sqrt{1-p}$（其 Kraus 算子为 $|g\rangle\langle g|+\sqrt{1-p}\,|e\rangle\langle e|$ 和
$\sqrt p\,|g\rangle\langle e|$），而 $\mathcal U_t$ 把 $\rho_{ge}$ 乘以 $e^{i\omega t}$；二者都是 TP 的，这就
确定了 $\rho_{gg}$。最后，$\rho(t)-|g\rangle\langle g|$ 是无迹的自伴算子，对角线为
$(-\rho_{ee}(t),\rho_{ee}(t))$，所以它的特征值为 $\pm(\rho_{ee}(t)^2+|\rho_{ge}(t)|^2)^{1/2}$。$\square$

相干项以布居数一半的速率衰减：取 $T_1=1/\gamma$ 和 $T_2=2/\gamma$，这就是 6.4 节中去掉退相位项的关系
$\frac1{T_2}=\frac1{2T_1}+\frac1{T_\varphi}$。最慢的非零衰减速率，即
**谱隙** $\min\{-\mathrm{Re}\,\mu:\mu\ne0\text{ 为 }\mathcal L\text{ 的特征值}\}=\gamma/2$，在这里就是
弛豫的速率：由上面最后一个公式，到 $|g\rangle\langle g|$ 的距离至多为 $e^{-\gamma t/2}$（由于正性给出
$|\rho_{ge}|^2\le\rho_{gg}\rho_{ee}$，有 $\rho_{ee}^2+|\rho_{ge}|^2\le\rho_{ee}\le1$），并且当 $t\to\infty$ 时，它乘以 $e^{\gamma t/2}$ 趋于 $|\rho_{ge}(0)|$。

*计算：* 用 `dynamics.lindblad_superoperator` 和 `dynamics.lindblad_channel` 检验 $\mathcal L$ 的谱、分解 $e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}}$ 以及命题 8.14 的迹距离——[教程，8.5 节](../../tutorial/zh/08_dynamics.md#s-8-5)。

<a id="s-8-6"></a>

## 8.6 可与之比较的工具

*计算：* 与其他软件（QuTiP、dynamiqs、QuantumOptics.jl、SciPy）的比较——[教程，8.6 节](../../tutorial/zh/08_dynamics.md#s-8-6)。

<a id="x-8"></a>

## 习题

1. （Heisenberg 绘景）对 8.3 节的生成元证明 $\mathrm{Tr}(\mathcal L(\rho)A)=\mathrm{Tr}(\rho\,\mathcal L^*(A))$，其中
   $\mathcal L^*(A)=i[H,A]+\sum_jg_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})+\sum_mg_m(\Phi_m^*(A)-\frac12\{\Phi_m^*(1),A\})$，
   证明 $\mathcal L^*(1)=0$，并且若没有扇区之间的跳跃，则 $\mathcal L^*(P_c)=0$。
2. （例外点）在 8.4 节的扇区 $N=1$ 中，
   $H_{\mathrm{eff}}=\begin{pmatrix}-i\gamma_1/2&J\\J&-i\gamma_2/2\end{pmatrix}$。求它的特征值，并证明
   当 $J=|\gamma_1-\gamma_2|/4>0$ 时它不可对角化。
3. 设 $m$ 个模式中的 $N_0$ 个玻色子各自以速率 $\gamma$ 丢失（跳跃 $a_1,\dots,a_m$ 限制在
   $N\to N-1$ 上）。证明对 $\mathcal M$ 中的每个哈密顿量，$p_N(t)=\binom{N_0}Nq^N(1-q)^{N_0-N}$，其中 $q=e^{-\gamma t}$。
4. （纯退相位）在 $M_2$ 上令 $H=\omega|e\rangle\langle e|$，$L=\sigma_z=\mathrm{diag}(1,-1)$，速率为
   $\gamma$。证明布居数不变，$\rho_{ge}(t)=e^{(i\omega-2\gamma)t}\rho_{ge}(0)$，并且尽管 $C=1$，
   稳态仍张成由对角矩阵构成的二维空间。

<details>
<summary>解答</summary>

1. 由循环性，$\mathrm{Tr}(L\rho L^*A)=\mathrm{Tr}(\rho L^*AL)$，
   $\mathrm{Tr}(\{L^*L,\rho\}A)=\mathrm{Tr}(\rho\{L^*L,A\})$，且
   $\mathrm{Tr}(-i[H,\rho]A)=\mathrm{Tr}(\rho\,i[H,A])$；对扇区之间的跳跃，
   $\mathrm{Tr}(\Phi(\rho)A)=\mathrm{Tr}(\rho\Phi^*(A))$ 就是 $\Phi^*$ 的定义
   （6.3 节）。于是 $\mathcal L^*(1)=\sum_jg_j(L_j^*L_j-L_j^*L_j)+\sum_mg_m(\Phi_m^*(1)-\Phi_m^*(1))=0$，
   这又一次是保迹性。没有扇区之间的跳跃时，$P_c$ 与 $H$ 及
   所有 $L_j$ 可交换，所以 $L_j^*P_cL_j=L_j^*L_jP_c=\frac12\{L_j^*L_j,P_c\}$，从而 $\mathcal L^*(P_c)=0$；于是
   $p_c'=\mathrm{Tr}(\rho\mathcal L^*(P_c))=0$，这是命题 8.9(1) 的第二个证明。

2. $H_{\mathrm{eff}}=-\frac i4(\gamma_1+\gamma_2)1+\begin{pmatrix}-i\delta&J\\J&i\delta\end{pmatrix}$，其中
   $\delta=(\gamma_1-\gamma_2)/4$；第二个矩阵的迹为 $0$，行列式为 $\delta^2-J^2$，所以
   特征值为 $-\frac i4(\gamma_1+\gamma_2)\pm\sqrt{J^2-\delta^2}$：当 $J<|\delta|$ 时为两个衰减速率，
   当 $J>|\delta|$ 时为一个速率和两个频率。当 $J=|\delta|>0$ 时，第二个矩阵 $Y$ 满足
   $Y^2=(J^2-\delta^2)1=0$ 但 $Y\ne0$，所以 $H_{\mathrm{eff}}$ 只有一个特征值
   $-\frac i4(\gamma_1+\gamma_2)$，却不是 $1$ 的倍数：它不可对角化。

3. 在扇区 $N$ 上 $\sum_ja_j^*a_j=N1$，并且唯一的扇区对是 $(N-1,N)$，所以命题
   8.9(2) 给出 $p_N'=-\gamma Np_N+\gamma(N+1)p_{N+1}$。由 $q'=-\gamma q$，二项式
   项的导数为 $-\gamma Np_N+\gamma(N_0-N)\binom{N_0}Nq^{N+1}(1-q)^{N_0-N-1}$，而
   $(N+1)\binom{N_0}{N+1}=(N_0-N)\binom{N_0}N$ 表明第二项为 $\gamma(N+1)p_{N+1}$。在 $t=0$
   时它是集中在 $N_0$ 处的点质量；唯一性由引理 8.1 得到。

4. $L\rho L^*-\frac12\{L^*L,\rho\}=\sigma_z\rho\sigma_z-\rho$ 把非对角元乘以 $-2$ 并
   消去对角元，且 $-i[H,\rho]_{ge}=i\omega\rho_{ge}$。所以 $\mathcal L(\sigma)=0$ 当且仅当 $\sigma_{ge}=0$（对
   $\gamma>0$），谱为 $\{0,0,-2\gamma\pm i\omega\}$。记 $U_t=e^{-itH}$，则
   $e^{t\mathcal L}(\rho)=\frac{1+e^{-2\gamma t}}2U_t\rho U_t^*+\frac{1-e^{-2\gamma t}}2U_t\sigma_z\rho\sigma_zU_t^*$：由于
   $\sigma_z$ 与 $U_t$ 可交换且 $\sigma_z^2=1$，求导表明右边 $R(t)$ 满足
   $R'=-i[H,R]+\gamma(\sigma_zR\sigma_z-R)=\mathcal L(R)$ 以及 $R(0)=\rho$，所以引理 8.1 适用。对 $t>0$，这是一个保单位信道，其不动
   点是与它的 Kraus 算子可交换的矩阵（定理 6.10），即对角矩阵。

</details>

*计算：* 习题 2 和 4 的数值检验——[教程，习题](../../tutorial/zh/08_dynamics.md#x-8)。

**下一章：** [第 9 章](09_spin_chains.md)讨论自旋链：磁化扇区与动量扇区，
大扇区的 Lanczos 与 Krylov 时间演化，以及无序。
