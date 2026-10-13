[English](../08_dynamics.md) · [Español](../es/08_dynamics.md) · [Français](../fr/08_dynamics.md) · [Deutsch](../de/08_dynamics.md) · **中文** · [日本語](../ja/08_dynamics.md) · [Русский](../ru/08_dynamics.md)

# 8. 动力学

封闭量子系统按照一个单参数酉群演化；开放系统——与一个不被观测的环境耦合的系统——则按照一个
信道半群演化。在有扇区的代数中，这两种演化都具有一个在单个矩阵代数上看不到的特点：
$\mathcal M$ 中的哈密顿量永远不能改变扇区，所以荷是守恒的；而开放系统可以在扇区之间跳跃，
例如由于丢失一个粒子。本章证明关于这两种演化的基本事实，写出带有扇区间跳跃的 Lindblad 方程，
并说明"没有粒子丢失"是对荷的一次测量，其后选择的态按照一个非厄米哈密顿量演化。

**你将学到**为什么扇区权重是运动常数，为什么 Gorini–Kossakowski–Sudarshan–Lindblad（GKSL）
形式的生成元产生信道，扇区间的跳跃如何进入 Lindblad 方程，以及扇区权重何时服从经典的速率方程，
无跳跃演化与后选择意味着什么，以及如何求稳态与弛豫速率。

在本章中，$\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ 作用在
$\mathcal H=\bigoplus_c\mathbb C^{k_c}$ 上，$D=\sum_ck_c$，$P_c$ 为极小中心投影，
态由其 blunt 密度 $\rho$ 给出，扇区权重为 $p_c=\operatorname{Tr}\rho_c$
（第 3 章），$[A,B]=AB-BA$，$\{A,B\}=AB+BA$，并取 $\hbar=1$。**超算子**是 $\mathcal M$ 或
$M_D(\mathbb C)$ 上的线性映射；超算子的范数一律取算子范数，因而是次乘性的。对有限维空间上的
线性映射 $B$，$e^{tB}=\sum_{n\ge0}t^nB^n/n!$。

## 8.1 封闭系统

**哈密顿量**是一个自伴元 $H\in\mathcal M$。向量按照 **Schrödinger
方程** $i\psi'(t)=H\psi(t)$ 演化，密度按照**冯·诺依曼方程** $\rho'(t)=-i[H,\rho(t)]$ 演化。
二者都是有限维空间上的常系数线性方程，下面的引理就是我们所需要的全部微分方程理论。

**引理 8.1.** 设 $B$ 是有限维赋范空间 $X$ 上的线性映射。则
$\frac d{dt}e^{tB}=Be^{tB}=e^{tB}B$，$e^{(s+t)B}=e^{sB}e^{tB}$，$\|e^{tB}\|\le e^{|t|\,\|B\|}$，并且对
每个 $x_0\in X$，映射 $x(t)=e^{tB}x_0$ 是区间
$I\ni0$ 上满足 $x'=Bx$ 与 $x(0)=x_0$ 的唯一可微函数 $x$。

*证明。* 该级数绝对收敛，并满足上述范数界；作为 $t$ 的幂级数，它可以
逐项求导，这就给出第一个公式。于是
$\frac d{dt}(e^{tB}e^{-tB})=Be^{tB}e^{-tB}-e^{tB}Be^{-tB}=0$，所以 $e^{tB}e^{-tB}=1$。若 $x'=Bx$，则
$\frac d{dt}(e^{-tB}x(t))=e^{-tB}(x'-Bx)=0$，所以 $e^{-tB}x(t)=x_0$，从而 $x(t)=e^{tB}x_0$。最后，
$t\mapsto e^{(s+t)B}x_0$ 与 $t\mapsto e^{tB}e^{sB}x_0$ 都是 $y'=By$ 的解，且初值相同。$\square$

**命题 8.2.** 设 $H=H^*\in\mathcal M$，$U(t)=e^{-itH}$，它由 1.3 节的函数演算
以 $f_t(\lambda)=e^{-it\lambda}$ 定义。

1. $U(t)=\bigoplus_ce^{-itH_c}$ 是 $\mathcal M$ 中的酉元，$U(s+t)=U(s)U(t)$，且 $U'(t)=-iHU(t)$；
   $U(t)$ 等于 $-itH$ 的指数级数。
2. $\psi(t)=U(t)\psi_0$ 是 $i\psi'=H\psi$ 满足 $\psi(0)=\psi_0$ 的唯一解，而
   $\rho(t)=U(t)\rho_0U(t)^*$ 是 $\rho'=-i[H,\rho]$ 在 $\mathcal M$ 中满足
   $\rho(0)=\rho_0$ 的唯一解。若 $\rho_0$ 是态，则每个 $\rho(t)$ 也是态。
3. （*Heisenberg 绘景*）$\operatorname{Tr}(\rho(t)A)=\operatorname{Tr}(\rho_0A(t))$，其中
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

**定理 8.3（荷守恒）。** 设 $t\mapsto H(t)$ 在区间 $I$ 上取值于 $\mathcal M$ 的
自伴元。

1. 若 $\psi:I\to\mathcal H$ 可微且 $i\psi'=H(t)\psi$，则每个 $\|P_c\psi(t)\|$ 都是常数。
2. 若 $\rho:I\to\mathcal M$ 可微且 $\rho'=-i[H(t),\rho]$，则对每个中心元 $Z$，$\operatorname{Tr}(Z\rho(t))$
   都是常数；特别地，每个扇区权重 $p_c(t)=\operatorname{Tr}(P_c\rho(t))$ 都是常数。
3. 若 $H$ 是常数且 $A\in\mathcal M$ 与 $H$ 可交换，则 $\operatorname{Tr}(A\rho(t))$ 是常数。

*证明。* (1) 内积关于第一个变元是共轭线性的，所以
$\frac d{dt}\langle\psi,P_c\psi\rangle=\langle-iH\psi,P_c\psi\rangle+\langle\psi,-iP_cH\psi\rangle
=i\langle\psi,HP_c\psi\rangle-i\langle\psi,P_cH\psi\rangle=0$，因为 $P_c$ 是中心元。(2)
由循环性，$\frac d{dt}\operatorname{Tr}(Z\rho)=-i\operatorname{Tr}(Z[H,\rho])=-i\operatorname{Tr}([Z,H]\rho)=0$。
(3) 对 $A$ 作同样的计算。$\square$

> **物理意义。** $\mathcal M$ 中的哈密顿量与荷（粒子数、
> 磁化、对称性量子数）可交换，所以每个荷值的概率都是运动常数，
> 即使对受驱动的、含时的哈密顿量也是如此。在计算上，演化分解为
> $C$ 个独立的问题：一次对角化的代价约为 $\sum_ck_c^3$ 次运算，而不是 $D^3$。

本库对每种块大小做一次特征分解来计算 $U(t)$（`dynamics.propagator`），而
`dynamics.von_neumann` 在本征基中一次性把密度演化到许多时刻，在该基中
$\rho(t)_{ij}=e^{-it(\lambda_i-\lambda_j)}\rho_{ij}$。我们用 `scipy.linalg.expm` 加以检验：

```python
import math
import numpy as np
import scipy.integrate
import scipy.linalg as sla
import torch
from torch_vn_algebra import TypeIAlgebra, DensityMatrix, InterSectorChannel, channels, dynamics, states

torch.manual_seed(0)
f64, c128 = torch.float64, torch.complex128
alg = TypeIAlgebra([2, 3, 4], [2, 3, 4], precision='double')     # M_2 (+) M_3 (+) M_4, D = 9
H = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
t = 0.7
U = dynamics.propagator(H, t)                                      # exp(-itH), block by block
for c, k in enumerate(alg.k_factors):                              # independent reference
    ref = sla.expm(-1j * t * H.matrix[0, c, :k, :k].numpy())
    assert np.abs(U.matrix[0, c, :k, :k].numpy() - ref).max() < 1e-12
assert (U.adjoint() @ U - alg.identity()).matrix.abs().max() < 1e-12

rho0 = states.random_density_matrix(alg)
times = [0.0, 0.5, 1.0, 2.0]
traj = dynamics.von_neumann(rho0, H, times)                        # exact, all times at once
for s, rho_s in zip(times, traj):
    U_s = dynamics.propagator(H, s)
    assert (rho_s.matrix - (U_s @ rho0 @ U_s.adjoint()).matrix).abs().max() < 1e-12
    assert (rho_s.sector_probabilities() - rho0.sector_probabilities()).abs().max() < 1e-12   # Theorem 8.3
A = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
A_t = U.adjoint() @ A @ U                                          # Heisenberg picture
assert abs(rho0.expectation(A_t) - dynamics.von_neumann(rho0, H, [t])[0].expectation(A)) < 1e-12
```

在 $\mathcal M$ 上，向量 $\psi=\bigoplus_c\psi_c$ 定义的是混合态
$\bigoplus_c|\psi_c\rangle\langle\psi_c|$ 的态（第 3 章习题 3），所以 `dynamics.schrodinger`（向量采用
块布局 `(batch, C, k_max)`）与 `von_neumann` 在它上面必须一致。对含时的 $H(t)$，
本库使用固定步长的经典 Runge–Kutta 方法（`dynamics.rk4`、`schrodinger_rk4`）：每一步把
$\psi$ 乘以 $H(t)\in\mathcal M$ 若干取值的一个多项式，所以它从不在扇区之间转移振幅，
但它不是严格酉的，扇区权重会以离散化误差的量级漂移。
（对于一个大的扇区，`krylov.evolve` 只需要矩阵–向量乘积；见第 9 章。）

```python
rows = torch.tensor([[1, 1, 0, 0], [1, 1, 1, 0], [1, 1, 1, 1]], dtype=c128)   # entries of C^2 (+) C^3 (+) C^4
psi0 = torch.randn(1, alg.C, alg.k_max, dtype=c128) * rows
psi0 = psi0 / psi0.norm()
psi = dynamics.schrodinger(psi0, H, times)                         # shape (T, batch, C, k_max)
p_vec = (psi.abs() ** 2).sum(-1)                                   # ||P_c psi(t)||^2
assert (p_vec - p_vec[0]).abs().max() < 1e-12
mix = DensityMatrix(alg, matrix=psi0[..., :, None] * psi0[..., None, :].conj())   # the state of psi0 on M
outer = psi[..., :, None] * psi[..., None, :].conj()
assert all((o - r.matrix).abs().max() < 1e-12 for o, r in zip(outer, dynamics.von_neumann(mix, H, times)))

V = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
H_of_t = lambda s: H.matrix + math.cos(3 * s) * V.matrix          # driven, in M at every time
psi_rk = dynamics.schrodinger_rk4(psi0, H_of_t, times, substeps=50, algebra=alg)
drift = ((psi_rk.abs() ** 2).sum(-1) - p_vec[0]).abs().max().item()
print(f"RK4, time-dependent H: largest change of a sector weight {drift:.1e}")
assert drift < 1e-8 and (dynamics.schrodinger_rk4(psi0, H, times, substeps=50) - psi).abs().max() < 1e-7

```

在 RK4 下扇区权重的变化为 $1.1\cdot10^{-9}$，而在精确求解器下除舍入误差外
没有任何变化。

## 8.2 量子动力学半群

CP 映射、TP 映射以及信道在 6.1 节（[第 6 章](06_quantum_channels.md)）中定义。除
定理 6.2 和命题 6.4（映射 $A\mapsto\sum_lK_lAK_l^*$——包括扇区之间的映射——是 CP 的）之外，我们还需要
三条封闭性质。

**引理 8.4.** (a) CP 映射的和、非负倍数与复合都是 CP 的。若
$\Phi$ 是 $M_D$ 上的 CP 映射，并把一个 $*$-子代数 $\mathcal N\subseteq M_D$（例如 $\mathcal M$，或
以 $P_c$ 为单位的 $P_c\mathcal MP_c\cong M_{k_c}$）映入自身，则 $\Phi|_{\mathcal N}$ 是 CP 的。
(b) 若 CP 映射 $\Phi_m$ 收敛于 $\Phi$（即对所有 $A$ 有 $\Phi_m(A)\to\Phi(A)$），则 $\Phi$ 是 CP 的。
(c) 若 $s\mapsto\Psi(s)$ 是从 $[a,b]$ 到 CP 映射的连续映射，则 $\int_a^b\Psi(s)\,ds$ 是 CP 的。

*证明。* (a) $\Phi\mapsto\Phi_n$ 是线性的，且 $(\Phi\Psi)_n=\Phi_n\Psi_n$；正矩阵之和是
正的。关于限制：$M_n(\mathcal N)\subseteq M_n(M_D)$，而正性在两者中是同一个
条件（经过基的置换，$M_n(\mathcal N)$ 的元素是
块 $[(A_{ij})_c]$ 与零的直和）。(b) 正矩阵构成闭集（对每个 $x$，$\langle x,Xx\rangle\ge0$
是闭条件），且 $(\Phi_m)_n(X)\to\Phi_n(X)$。(c) 由 (a)，Riemann 和是 CP 的，并且由于
$\Psi$ 一致连续，它们收敛于该积分；应用 (b)。$\square$

**定义。** $\mathcal M$ 上的**量子动力学半群**是一族信道
$\mathcal M\to\mathcal M$，$(T_t)_{t\ge0}$，满足 $T_0=\mathrm{id}$、$T_{s+t}=T_sT_t$，且 $t\mapsto T_t$ 连续。

**引理 8.5.** 每个连续的超算子半群（$T_0=\mathrm{id}$，$T_{s+t}=T_sT_t$）都可以写成 $T_t=e^{t\mathcal L}$，
其中超算子 $\mathcal L$ 是唯一确定的，称为它的**生成元**，并且 $\mathcal L=\lim_{t\to0}(T_t-\mathrm{id})/t$。

*证明。* 当 $h\to0$ 时 $V_h=\frac1h\int_0^hT_s\,ds\to\mathrm{id}$，所以对某个固定的小 $h$，
$V_h$ 可逆。于是 $T_tV_h=\frac1h\int_t^{t+h}T_s\,ds$ 关于 $t$ 可微，导数为
$\frac1h(T_{t+h}-T_t)=T_t\frac{T_h-\mathrm{id}}h$。因此 $T_t$ 可微，且 $T_t'=T_t\mathcal L$，
其中 $\mathcal L=\frac{T_h-\mathrm{id}}hV_h^{-1}$，再由引理 8.1（取 $B$ 为
右复合 $\mathcal L$）得 $T_t=e^{t\mathcal L}$。唯一性：$\mathcal L=T'(0)$。$\square$

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
定理给出 $S_n'=\mathcal L_0S_n+\Phi S_{n-1}$（$n\ge1$）以及 $S_0'=\mathcal L_0S_0$。由这些
界，导数的级数在 $[0,T]$ 上一致收敛，所以 $Y(t)=\sum_nS_n(t)$
可微，$Y'=\mathcal L_0Y+\Phi Y=\mathcal LY$，$Y(0)=\mathrm{id}$，再由引理 8.1 得 $Y(t)=e^{t\mathcal L}$。
*第 4 步。* $\Phi$ 是 CP 的（定理 6.2），所以由归纳法，被积函数
$s\mapsto S_0(t-s)\Phi S_{n-1}(s)$ 是一族连续的 CP 映射，从而 $S_n(t)$ 是 CP 的（引理 8.4）；
部分和及其极限 $e^{t\mathcal L}$ 也是如此。*第 5 步。* 由循环性，
$\operatorname{Tr}\mathcal L(\rho)=\sum_j\big(\operatorname{Tr}(L_j\rho L_j^*)-\operatorname{Tr}(L_j^*L_j\rho)\big)=0$，
因此 $\frac d{dt}\operatorname{Tr}e^{t\mathcal L}\rho=\operatorname{Tr}\mathcal L(e^{t\mathcal L}\rho)=0$：每个
$e^{t\mathcal L}$ 都是 TP 的。半群律与连续性来自引理 8.1。$\square$

**定理 8.7（Gorini–Kossakowski–Sudarshan；Lindblad；不加证明地陈述）。** 反之，$M_D$ 上每个量子
动力学半群的生成元都具有定理 8.6 的形式。[V. Gorini, A. Kossakowski,
E. C. G. Sudarshan, *J. Math. Phys.* 17, 821 (1976); G. Lindblad, *Commun. Math. Phys.* 48, 119
(1976).] 对一般冯·诺依曼代数 $\mathcal M$ 上由保单位 CP 映射构成的范数连续半群，
Christensen 与 Evans [*J. London Math. Soc.* 20, 358 (1979)] 证明了其生成元
为 $A\mapsto\Psi(A)+K^*A+AK$，其中 $\Psi$ 是 CP 的，$K\in\mathcal M$。

在 $\mathcal M=\bigoplus_cM_{k_c}$ 上，把这一结果应用于 Heisenberg 绘景中的半群 $(T_t^*)$：它的
映射是 CP 且保单位的（命题 6.8，推论 6.9(b)），其生成元为 $\mathcal L^*$。
再取 blunt 对偶回来，得到 $\mathcal L(\rho)=\Phi(\rho)+G\rho+\rho G^*$，其中 $\Phi=\Psi^*$ 是
$\mathcal M$ 上的 CP 映射（扇区之间的 Kraus 算子，命题 6.4），$G=K$。保迹性，即对所有 $\rho\in\mathcal M$ 有 $\operatorname{Tr}\big((\Phi^*(1)+G+G^*)\rho\big)=0$，迫使
$G+G^*=-\Phi^*(1)$，即 $G=-iH-\frac12\Phi^*(1)$，其中 $H=H^*\in\mathcal M$：这恰好是
下一节的生成元。

> **物理意义。** 一个与大环境弱相互作用、而环境的关联
> 衰减很快（Born–Markov 近似）的系统会丧失对其过去的记忆：它在时刻 $s+t$ 的态
> 由时刻 $s$ 的态经过一个只依赖于 $t$ 的信道得到。定理 8.7 表明，
> 这样的动力学总是一个哈密顿部分加上"跳跃" $L_j$（光子发射、粒子丢失、
> 退相位），其速率由 $L_j$ 的归一化确定。

本库有四个 Lindblad 例程：`lindblad_rhs` 以函数形式返回右端项，供
`rk4` 使用；`lindblad_evolve` 用 RK4 积分；`lindblad_superoperator` 对每个扇区返回 $\mathcal L$ 作用在
按行展开的 $\operatorname{vec}\rho_c$ 上的 $k_c^2\times k_c^2$ 矩阵；
`lindblad_channel` 逐扇区精确计算 $e^{t\mathcal L}$，结果为一个 `Channel`。我们对
$M_3$ 上任意的跳跃检验定理 8.6：$e^{t\mathcal L}$ 的 Choi 矩阵是正的（定理
6.2），该映射是 TP 的，半群律成立，RK4 与指数一致，并且跳跃部分
$e^{t\mathcal L}-S_0(t)=\sum_{n\ge1}S_n(t)$ 本身是 CP 的：

```python
M3 = TypeIAlgebra([3], [3], precision='double')
H3 = M3.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
jumps = [M3.operator(torch.randn(1, 1, 3, 3, dtype=c128) / 2) for _ in range(2)]     # arbitrary L_j
T = lambda s: dynamics.lindblad_channel(M3, H3, jumps, t=s)        # exp(sL), exact
S = lambda s: T(s).superoperator()[0][0]                           # its 9 x 9 matrix on row-major vec
choi = lambda X: X.reshape(3, 3, 3, 3).permute(2, 0, 3, 1).reshape(9, 9)   # J_{(i,a),(j,b)} = S_{(a,b),(i,j)}
assert (choi(S(0.8)) - T(0.8).choi()[0][0]).abs().max() < 1e-12
assert torch.linalg.eigvalsh(choi(S(0.8))).min() > -1e-12         # completely positive
assert T(0.8).trace_preservation_error().max() < 1e-12             # trace preserving
assert (S(0.3) @ S(0.5) - S(0.8)).abs().max() < 1e-12              # semigroup

rho3 = states.random_density_matrix(M3)
rk = dynamics.lindblad_evolve(rho3, H3, jumps, [0.0, 0.8], substeps=40)[-1]
assert (rk.matrix - T(0.8)(rho3).matrix).abs().max() < 1e-7        # RK4 against exp(tL)
assert abs(dynamics.lindblad_rhs(H3, jumps)(0.0, rho3.matrix)[0, 0].trace()) < 1e-12   # Tr L(rho) = 0

H_eff = H3.matrix - 0.5j * sum(L.adjoint().matrix @ L.matrix for L in jumps)
S0 = channels.Channel(M3, torch.linalg.matrix_exp(-0.8j * H_eff)[:, None]).superoperator()[0][0]
jump_part = torch.linalg.eigvalsh(choi(S(0.8) - S0))
print(f"smallest eigenvalue of the Choi matrix of the jump part: {jump_part.min().item():.1e}")
assert jump_part.min() > -1e-12                                    # sum_{n >= 1} S_n(t) is CP
```

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

**命题 8.8.** 设 $J^{dc}_l\in M_D$ 是在位置 $(d,c)$ 处为块 $K^{dc}_l$、
其余处为零的矩阵。

1. $\mathcal L$ 把 $\mathcal M$ 映入 $\mathcal M$，并且是 $M_D$ 上以 $\sqrt{g_j}L_j$ 和 $\sqrt{g_m}J^{dc}_{m,l}$（所有 $m,l,d,c$）为跳跃的 GKSL 生成元
   $\hat{\mathcal L}$ 在 $\mathcal M$ 上的限制。因此
   $(e^{t\mathcal L})_{t\ge0}$ 是 $\mathcal M$ 上的量子动力学半群。
2. 另设 $\tilde{\mathcal L}$ 是 $M_D$ 上的 GKSL 生成元，其中同一个 Kraus
   指标的各块被加成一个跳跃，$L_{m,l}=\sum_{d,c}J^{dc}_{m,l}$。则对 $\rho\in\mathcal M$ 有
   $\mathcal L(\rho)=E\big(\tilde{\mathcal L}(\rho)\big)$，其中 $E(X)=\sum_dP_dXP_d$ 是
   $M_D$ 到 $\mathcal M$ 上的条件期望；一般而言 $\tilde{\mathcal L}(\rho)\notin\mathcal M$。

*证明。* 由于 $\rho=\sum_cP_c\rho P_c$ 且 $J^{dc}P_{c'}=\delta_{cc'}J^{dc}$，
$J^{dc}\rho J^{dc*}$ 就是位于位置 $(d,d)$ 的块 $K^{dc}\rho_cK^{dc*}$，而 $J^{dc*}J^{dc}$ 是位于位置 $(c,c)$ 的
$K^{dc*}K^{dc}$。对 $d,l$ 求和得到 $\Phi(\rho)$ 与 $\Phi^*(1)$，这就证明了 (1) 中的
恒等式，并且所有项都属于 $\mathcal M$。所以 $\hat{\mathcal L}$ 及其所有幂都把 $\mathcal M$
映入 $\mathcal M$，$e^{t\hat{\mathcal L}}|_{\mathcal M}=e^{t\mathcal L}$，而由
定理 8.6 和引理 8.4(a)，它是 $\mathcal M$ 的一个信道。(2) 对 $L=\sum_{d,c}J^{dc}$ 有
$L\rho L^*=\sum_{d,d',c}J^{dc}\rho J^{d'c*}$，而
$E$ 保留 $d=d'$ 的项。对 $B\in\mathcal M$，$P_dXBP_d=P_dXP_dB$，所以 $E(XB)=E(X)B$，
$E(BX)=BE(X)$；于是 $E\{L^*L,\rho\}=\{E(L^*L),\rho\}$，并且 $E(L^*L)=\sum_cP_cL^*LP_c=\sum_{c,d}J^{dc*}J^{dc}$，
因为当 $d\ne d'$ 时 $J^{dc*}J^{d'c}=0$。最后 $E(-i[H,\rho])=-i[H,\rho]$。$\square$

> **物理意义。** 第 (2) 部分表明，$\mathcal M$ 上的 Lindblad 方程是
> $M_D$ 上的方程在跳跃所产生的扇区间相干性被立即抹去之后剩下的部分——
> 这正是超选择定则的内容。对于粒子丢失，没有任何东西被抹去：湮灭
> 算子 $a_j$ 恰好把粒子数降低一，所以 $a_j\rho a_j^*$ 和 $a_j^*a_j=n_j$ 都是
> 块对角的，$\mathcal M$ 上的方程就是通常的方程。

**命题 8.9（扇区权重）。** 设 $\rho(t)=e^{t\mathcal L}\rho_0$，$p(t)$ 为其扇区权重。

1. 若没有扇区之间的跳跃，则每个 $p_c(t)$ 都是常数。
2. 若对每个 $m$ 以及每一对 $d\ne c$，有 $\sum_lK^{dc*}_{m,l}K^{dc}_{m,l}=r_{m,dc}1_{k_c}$，其中
   $r_{m,dc}$ 为数，则 $p'=Qp$，其中对 $d\ne c$ 有 $Q_{dc}=\sum_mg_mr_{m,dc}$，且 $Q_{cc}=-\sum_{d\ne c}Q_{dc}$。
   因此 $p(t)=e^{tQ}p(0)$ 与 $H$、各 $L_j$ 以及扇区内部的态都无关。

*证明。* $p_d'=\operatorname{Tr}(P_d\mathcal L(\rho))$。哈密顿项如定理 8.3 中那样为零。
对 $L_j\in\mathcal M$，$\operatorname{Tr}(P_dL_j\rho L_j^*)-\operatorname{Tr}(P_dL_j^*L_j\rho)
=\operatorname{Tr}(L_j^*P_dL_j\rho)-\operatorname{Tr}(L_j^*L_jP_d\rho)=0$，因为 $P_d$ 是中心元。对扇区之间的
跳跃 $\Phi$，$\operatorname{Tr}(P_d\Phi(\rho))=\sum_c\operatorname{Tr}\big(\sum_lK^{dc*}_lK^{dc}_l\rho_c\big)$，
且 $\operatorname{Tr}(P_d\tfrac12\{\Phi^*(1),\rho\})=\sum_{d'}\operatorname{Tr}\big(\sum_lK^{d'd*}_lK^{d'd}_l\rho_d\big)$；
$c=d$ 与 $d'=d$ 的项相互抵消。这证明了 (1)；在 (2) 的假设下，剩下的项
为 $\sum_{c\ne d}r_{dc}p_c-\sum_{d'\ne d}r_{d'd}p_d$。由引理 8.1 得 $p(t)=e^{tQ}p(0)$。$\square$

在 (2) 的假设下，对 $d\ne c$，$r_{m,dc}$ 对每个态 $\rho$ 都等于 $\Phi_m$ 的转移
矩阵（在命题 6.12 之前定义）的元素 $T_{dc}(\rho)$。所以，对速率为
$1$ 的单个跳跃，$Q=R-\operatorname{diag}(\sum_dR_{dc})_c$，其中 $R$ = `transition_matrix()`（$R$ 的对角元相互抵消）。

**例：粒子丢失。** 两个玻色模式，至多两个粒子；扇区为
$N=0,1,2$，基分别为 $|00\rangle$；$|10\rangle,|01\rangle$；$|20\rangle,|11\rangle,|02\rangle$，所以
$\mathcal M=M_1\oplus M_2\oplus M_3$。Bose–Hubbard 二聚体
$H=J(a_1^*a_2+a_2^*a_1)+\frac U2\sum_jn_j(n_j-1)$ 与 $N$ 可交换，并且每个模式以
速率 $\gamma$ 丢失粒子：这是扇区之间的一个跳跃，其 Kraus 算子为 $a_1,a_2$，限制在 $N\to N-1$ 上。
由于 $a_1^*a_1+a_2^*a_2=N$，命题 8.9(2) 适用，$Q_{N-1,N}=\gamma N$，并且从 $N=2$ 出发
$$
p_2(t)=q^2,\qquad p_1(t)=2q(1-q),\qquad p_0(t)=(1-q)^2,\qquad q=e^{-\gamma t},
$$
这是二项分布：无论隧穿与相互作用如何，每个粒子都以概率 $e^{-\gamma t}$ 独立地存活
（通过求导验证 $p'=Qp$）。

```python
fock = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double', charges=[0, 1, 2])   # N = 0, 1, 2
r2 = math.sqrt(2)
a1 = {(0, 1): [[1.0, 0.0]], (1, 2): [[r2, 0.0, 0.0], [0.0, 1.0, 0.0]]}           # a_1 : N -> N - 1
a2 = {(0, 1): [[0.0, 1.0]], (1, 2): [[0.0, 1.0, 0.0], [0.0, 0.0, r2]]}           # a_2 : N -> N - 1
blocks = lambda a: {dc: torch.tensor(K, dtype=f64) for dc, K in a.items()}
loss = InterSectorChannel.from_blocks(fock, fock, {dc: [blocks(a1)[dc], blocks(a2)[dc]] for dc in a1})
J, U_int = 1.0, 0.8                                                # hopping and interaction
H_bh = fock.from_blocks([torch.zeros(1, 1, dtype=f64), torch.tensor([[0.0, J], [J, 0.0]], dtype=f64),
                         torch.tensor([[U_int, r2 * J, 0.0], [r2 * J, 0.0, r2 * J], [0.0, r2 * J, U_int]], dtype=f64)])

assert (loss.adjoint()(fock.identity()).matrix - fock.central([0.0, 1.0, 2.0]).matrix).abs().max() < 1e-12
R = loss.transition_matrix()[0]                                    # r_dc (columns: from sector c)
Q = R - torch.diag(R.sum(0))                                       # Proposition 8.9(2), for rate 1
print(Q)

gamma = 0.5
psi2 = torch.randn(3, dtype=c128)
psi2 = psi2 / psi2.norm()
m0 = torch.zeros(1, 3, 3, 3, dtype=c128)
m0[0, 2] = torch.outer(psi2, psi2.conj())                          # a pure state with N = 2
ts = [0.0, 0.5, 1.0, 2.0, 4.0]
traj = dynamics.lindblad_evolve(DensityMatrix(fock, matrix=m0), H_bh, [loss], ts, rates=[gamma], substeps=100)
for s, rho_s in zip(ts, traj):
    q = math.exp(-gamma * s)
    p = rho_s.sector_probabilities()[0]
    assert (p - torch.tensor([(1 - q) ** 2, 2 * q * (1 - q), q ** 2], dtype=f64)).abs().max() < 1e-8
    assert (p - torch.from_numpy(sla.expm(gamma * s * Q.numpy())[:, 2])).abs().max() < 1e-8
neg = min(w.min().item() for rho_s in traj for w in rho_s.eigenvalues())
print(f"smallest eigenvalue of rho(t) under RK4: {neg:.1e}")
assert neg > -1e-7                                                 # positive up to the RK4 error
```

RK4 不是严格完全正的：$\rho(t)$ 的特征值可以在
离散化误差的量级上为负（这里为 $\approx-2\cdot10^{-8}$，步长至多 $0.02$）；`lindblad_channel` 是精确的，
但只处理扇区内部的跳跃。最后，我们对扇区之间的一个随机映射检验命题 8.8(2)，
并与用 numpy 构造的 $M_6$ 上的 GKSL 生成元进行比较：

```python
Phi = channels.random_inter_sector_channel(fock, fock, kraus_rank=2)
rho = states.random_density_matrix(fock)
on_M = dynamics.lindblad_rhs(None, [Phi], algebra=fock)(0.0, rho.matrix)[0].numpy()
ks, K = fock.k_factors, Phi.kraus[0].numpy()                       # K[l, d, c] = K_l^{dc}, zero padded
embed = lambda Y: sla.block_diag(*[Y[c, :k, :k] for c, k in enumerate(ks)])    # (+)_c Y_c in M_6
Ls = [np.block([[K[l, d, c, :kd, :kc] for c, kc in enumerate(ks)] for d, kd in enumerate(ks)])
      for l in range(2)]                                           # L_l = sum_{d,c} J_l^{dc}
X = embed(rho.matrix[0].numpy())
full = sum(L @ X @ L.conj().T - 0.5 * (L.conj().T @ L @ X + X @ L.conj().T @ L) for L in Ls)   # GKSL on M_6
E_full = full * embed(np.ones((3, 3, 3)))                          # E(X) = sum_d P_d X P_d
assert np.abs(E_full - embed(on_M)).max() < 1e-12                  # Proposition 8.8(2)
print(f"coherences between sectors created on M_6: {np.abs(full - E_full).max():.2f}")
```

## 8.4 后选择与无跳跃演化

把 8.3 节生成元的跳跃汇集为 $\Phi_{\mathrm{tot}}(\rho)=\sum_jg_jL_j\rho L_j^*+\sum_mg_m\Phi_m(\rho)$，
令 $\mathcal M$ 中的 $A=\Phi_{\mathrm{tot}}^*(1)=\sum_jg_jL_j^*L_j+\sum_mg_m\Phi_m^*(1)\ge0$，并令
$H_{\mathrm{eff}}=H-\frac i2A$。由命题 8.8(1) 和定理 8.6（其 $S_n$ 把 $\mathcal M$ 映入
$\mathcal M$，因为 $H_{\mathrm{eff}}\in\mathcal M$），
$$
\rho(t)=\sum_{n\ge0}S_n(t)\rho,\qquad S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},
\qquad S_n(t)=\int_0^tS_0(t-s)\,\Phi_{\mathrm{tot}}\,S_{n-1}(s)\,ds .
$$

**命题 8.10.** 对态 $\rho$，令 $w_n(t)=\operatorname{Tr}S_n(t)\rho$。则 $w_n(t)\ge0$，
$\sum_nw_n(t)=1$，且 $w_0$ 不增，$w_0'(t)=-\operatorname{Tr}\big(A\,S_0(t)\rho\big)$。对
纯态 $\rho=|\psi\rangle\langle\psi|$，$w_0(t)=\|e^{-itH_{\mathrm{eff}}}\psi\|^2$。

*证明。* $S_n(t)$ 是 CP 的，因而是正的，且 $\sum_nw_n=\operatorname{Tr}\rho(t)=1$。令
$X=S_0(t)\rho\ge0$，则 $w_0'=\operatorname{Tr}(-iH_{\mathrm{eff}}X+iXH_{\mathrm{eff}}^*)=-\operatorname{Tr}(AX)$，
且 $\operatorname{Tr}(AX)=\operatorname{Tr}(X^{1/2}AX^{1/2})\ge0$。最后一个结论就是
$\operatorname{Tr}|\phi\rangle\langle\phi|=\|\phi\|^2$。$\square$

> **物理意义。** 设一个探测器记录每一次跳跃——每个发射的光子、每个丢失的原子。
> 则 $w_n(t)$ 是到时刻 $t$ 为止恰好有 $n$ 次计数的概率，$S_n(t)\rho/w_n(t)$ 是
> 给定 $n$ 次计数时的态，而 $S_0(t)\rho/w_0(t)$ 是给定**没有计数**时的态：零结果同样会
> 改变态，这种改变通过非厄米的 $H_{\mathrm{eff}}$ 实现。这种计数过程诠释 [E. B. Davies,
> *Quantum Theory of Open Systems*, 1976; M. D. Srinivas and E. B. Davies, *Optica Acta* 28, 981 (1981)]
> 是量子轨迹（Monte Carlo 波函数）方法的基础 [J. Dalibard, Y. Castin,
> K. Mølmer, *Phys. Rev. Lett.* 68, 580 (1992)]；这里我们只用到它的数学部分。

对于损失，计数就等同于测量荷：

**命题 8.11.** 设扇区标记为 $N=0,\dots,C-1$，扇区内部没有跳跃，
每个跳跃都把荷降低一（除非 $d=c-1$，否则 $K^{dc}_{m,l}=0$），并且
$\rho=P_{N_0}\rho P_{N_0}$。则对 $0\le n\le N_0$ 有 $S_n(t)\rho=P_{N_0-n}\,\rho(t)\,P_{N_0-n}$，
对 $n>N_0$ 有 $S_n(t)\rho=0$。特别地，$\rho(t)$ 关于 $P_{N_0}$ 的 Lüders 更新（"没有粒子丢失"）
就是无跳跃态：
$$
\frac{P_{N_0}\rho(t)P_{N_0}}{p_{N_0}(t)}=\frac{e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*}}{w_0(t)},
\qquad p_{N_0}(t)=w_0(t).
$$

*证明。* 由于 $H_{\mathrm{eff}}\in\mathcal M$，$S_0(s)$ 把支撑在扇区 $N$ 中的算子映为支撑在扇区 $N$ 中的算子；
$\Phi_{\mathrm{tot}}$ 把扇区 $N\ge1$ 映到扇区 $N-1$，并把扇区 $0$ 映为零。
由归纳法，对所有 $s$，$S_n(s)\rho$ 支撑在扇区 $N_0-n$ 中，而当 $n>N_0$ 时为 $0$。用 $P_{N_0-n}$ 压缩
$\rho(t)=\sum_nS_n(t)\rho$ 恰好保留第 $n$ 项。$\square$

本库的 `states.lueders_update(rho, P)` 返回 $P\rho P/\operatorname{Tr}(P\rho P)$ 和
$\operatorname{Tr}(\rho P)$。现在让两个模式以不同的速率 $\gamma_1>\gamma_2$ 泄漏；则
$A=\gamma_1n_1+\gamma_2n_2$ 不是中心元，命题 8.9(2) 不适用，而对
$N=2$ 的后选择给出真正非厄米的演化。我们对 $n=0$ 检验命题 8.11，并通过用 Simpson 公式计算
$S_1$ 的积分，对 $n=1$ 检验它：

```python
g1, g2 = 0.6, 0.2                                                  # mode 1 leaks faster
loss1 = InterSectorChannel.from_blocks(fock, fock, {dc: [K] for dc, K in blocks(a1).items()})
loss2 = InterSectorChannel.from_blocks(fock, fock, {dc: [K] for dc, K in blocks(a2).items()})
m1 = torch.zeros(1, 3, 3, 3, dtype=c128)
m1[0, 2, 0, 0] = 1.0                                               # |20>: both particles in mode 1
ts = [0.0, 0.5, 1.0, 2.0]
traj = dynamics.lindblad_evolve(DensityMatrix(fock, matrix=m1), H_bh, [loss1, loss2], ts, rates=[g1, g2],
                                substeps=100)

n_1 = [np.diag(x) for x in ([0.0], [1.0, 0.0], [2.0, 1.0, 0.0])]   # n_1 in each sector; n_2 = N - n_1
H_eff = [H_bh.matrix[0, N, :N + 1, :N + 1].numpy() - 0.5j * (g1 * n + g2 * (N * np.eye(N + 1) - n))
         for N, n in enumerate(n_1)]

def S0(s, N, X):                                  # e^{-isH_eff} X e^{isH_eff^*} in sector N, s a time or an array
    E = torch.linalg.matrix_exp(torch.from_numpy(-1j * np.asarray(s)[..., None, None] * H_eff[N])).numpy()
    return E @ X @ E.conj().swapaxes(-1, -2)

rho_20 = m1[0, 2].numpy()
for s, rho_s in zip(ts, traj):
    w0 = np.trace(S0(s, 2, rho_20)).real                           # probability of no loss
    post, prob = states.lueders_update(rho_s, fock.central([0.0, 0.0, 1.0]))   # post-select N = 2
    assert abs(prob.item() - w0) < 1e-8                            # Proposition 8.11, n = 0
    assert np.abs(post.matrix[0, 2].numpy() - S0(s, 2, rho_20) / w0).max() < 1e-7
closed = dynamics.von_neumann(DensityMatrix(fock, matrix=m1), H_bh, [2.0])[0]   # no losses at all
dist = 0.5 * np.abs(np.linalg.eigvalsh(S0(2.0, 2, rho_20) / w0 - closed.matrix[0, 2].numpy())).sum()
print(f"t = 2: P(no loss) = {w0:.3f} (without hopping: {math.exp(-2 * g1 * 2):.3f}); "
      f"trace distance post-selected/closed system: {dist:.3f}")

K1, K2 = blocks(a1)[(1, 2)].numpy(), blocks(a2)[(1, 2)].numpy()   # S_1(2) rho = int_0^2 S_0(2-s) Phi S_0(s) rho ds
jump = lambda X: g1 * K1 @ X @ K1.T + g2 * K2 @ X @ K2.T
s_grid = np.linspace(0.0, 2.0, 201)
S1 = scipy.integrate.simpson(S0(2.0 - s_grid, 1, jump(S0(s_grid, 2, rho_20))), x=s_grid, axis=0)
assert np.abs(S1 - traj[-1].matrix[0, 1, :2, :2].numpy()).max() < 1e-8     # Proposition 8.11, n = 1
```

在 $t=2$ 时，没有损失的概率为 $0.230$；若没有隧穿，它将是 $e^{-2\gamma_1t}=0.091$，所以
这里扇区权重确实依赖于哈密顿量。后选择的态与封闭系统的态之间的迹距离为
$0.135$：没有丢失粒子本身就是信息。

## 8.5 稳态与弛豫

若对所有 $t\ge0$ 有 $e^{t\mathcal L}\sigma=\sigma$，则称态 $\sigma$ 是**稳态**，等价地
$\mathcal L(\sigma)=0$（在 $t=0$ 处求导；反之利用级数）。

**命题 8.12.** $\mathcal M$ 上的每个量子动力学半群都有稳态。

*证明。* 固定一个态 $\rho$，令 $\sigma_\tau=\frac1\tau\int_0^\tau T_s\rho\,ds$。它是一个态（态的
凸组合的极限；态构成闭集）。由线性，
$T_t\sigma_\tau-\sigma_\tau=\frac1\tau\big(\int_\tau^{\tau+t}-\int_0^t\big)T_s\rho\,ds$，由于 $\|T_s\rho\|_1=1$，
其迹范数至多为 $2t/\tau$。态在有限维
空间中构成有界闭集，所以某个序列 $\sigma_{\tau_m}$，$\tau_m\to\infty$，收敛于一个态 $\sigma$
（Bolzano–Weierstrass），并且 $T_t\sigma-\sigma=\lim_m(T_t\sigma_{\tau_m}-\sigma_{\tau_m})=0$。$\square$

**推论 8.13.** 若没有扇区之间的跳跃，则稳态张成的空间
维数至少为 $C$：每个扇区都有它自己的稳态。

*证明。* $\mathcal L$ 把 $P_c\mathcal MP_c\cong M_{k_c}$ 映入自身（$H$、$L_j$ 全都是块对角的），
所以 $e^{t\mathcal L}$ 的限制是 $M_{k_c}$ 上的量子动力学半群（引理 8.4(a)），
由命题 8.12，它有稳态 $\sigma_c$。各 $\sigma_c$ 的支撑互不相交。$\square$

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

相干项以布居数一半的速率衰减，这就是 6.4 节中没有
纯退相位时的关系 $T_2=2T_1$；
最慢的非零速率，即**谱隙** $\gamma/2$，控制着弛豫。我们用
`lindblad_superoperator` 检验谱，用 `lindblad_channel` 检验信道与迹距离，并对
$M_2\oplus M_3$ 扇区内部的随机跳跃检验推论 8.13：

```python
q2 = TypeIAlgebra([2], [2], precision='double')                    # basis |g>, |e>
omega, gamma = 1.3, 0.4
H_at = q2.from_blocks([torch.tensor([[0.0, 0.0], [0.0, omega]], dtype=f64)])
sigma_minus = q2.from_blocks([torch.tensor([[0.0, 1.0], [0.0, 0.0]], dtype=f64)])   # |g><e|
gen = dynamics.lindblad_superoperator(q2, H_at, [sigma_minus], rates=[gamma])[0][0]
spec = sorted(torch.linalg.eigvals(gen).tolist(), key=lambda z: (z.real, z.imag))
print([f"{z.real:+.2f}{z.imag:+.2f}i" for z in spec])
assert np.allclose(spec, [-gamma, -gamma / 2 - 1j * omega, -gamma / 2 + 1j * omega, 0.0], atol=1e-12)

psi_at = torch.tensor([0.6, 0.8], dtype=c128)                      # rho_ee(0) = 0.64, rho_ge(0) = 0.48
rho_at = DensityMatrix(q2, matrix=torch.outer(psi_at, psi_at.conj())[None, None])
ground = DensityMatrix(q2, matrix=torch.tensor([[[[1.0, 0.0], [0.0, 0.0]]]], dtype=c128))
for s in (0.5, 2.0, 10.0):
    ch = dynamics.lindblad_channel(q2, H_at, [sigma_minus], t=s, rates=[gamma])
    UA = channels.unitary_channel(dynamics.propagator(H_at, s)).compose(
        channels.amplitude_damping_channel(q2, 1 - math.exp(-gamma * s)))
    assert (ch.superoperator()[0] - UA.superoperator()[0]).abs().max() < 1e-12
    r = ch(rho_at).matrix[0, 0]
    assert abs(r[1, 1].real - 0.64 * math.exp(-gamma * s)) < 1e-12
    assert abs(r[0, 1] - 0.48 * np.exp((1j * omega - gamma / 2) * s)) < 1e-12
    exact = math.sqrt(0.64 ** 2 * math.exp(-2 * gamma * s) + 0.48 ** 2 * math.exp(-gamma * s))
    assert abs(states.trace_distance(ch(rho_at), ground).item() - exact) < 1e-12

two = TypeIAlgebra([2, 3], [2, 3], precision='double')
H_two = two.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
L_two = [two.operator(torch.randn(1, 2, 3, 3, dtype=c128)) for _ in range(2)]
kernel = [int((torch.linalg.eigvals(G_c[0]).abs() < 1e-10).sum())
          for G_c in dynamics.lindblad_superoperator(two, H_two, L_two)]
print("stationary states per sector:", kernel)
assert all(dim >= 1 for dim in kernel)                             # Corollary 8.13
```

对这些随机跳跃，代码在每个扇区中恰好找到一个稳态。唯一性的
判据——例如，一个忠实的稳态再加上 $\{H,L_j,L_j^*\}'=\mathbb C1$
[A. Frigerio, *Commun. Math. Phys.* 63, 269 (1978)]——这里不予证明；习题 4 表明，
即使在因子中，唯一性也可能不成立。

## 8.6 可与之比较的工具

QuTiP（`sesolve`；`mesolve` 以自适应 ODE 求解器求解 Lindblad 方程；`mcsolve` 用于
量子轨迹；`steadystate`；`liouvillian`）、dynamiqs（JAX；在 GPU 上成批且可微的
Schrödinger 与 Lindblad 求解器）、QuantumOptics.jl 以及 SciPy（`scipy.linalg.expm`，以及对稀疏的
向量化生成元计算 $e^{t\mathcal L}\rho$ 的 `scipy.sparse.linalg.expm_multiply`）在单个 Hilbert 空间上
求解同样的方程，对于上面的问题它们必须给出相同的数值。
在 $\mathcal M$ 上，需要把 $D\times D$ 的块对角矩阵传给它们；这时向量化的生成元
位于一个 $D^2$ 维空间上（稀疏，但各扇区在一次求解中耦合在一起），除非用户
手动把问题拆分到各扇区。本库内建了块结构：逐扇区的精确传播子
与 Lindblad 信道（$k_c^2\times k_c^2$ 超算子），作为
`InterSectorChannel` 的扇区间跳跃以及命题 8.8 的压缩，成批的哈密顿量与跳跃
算子（用于无序平均），以及自动微分（后两者 dynamiqs 也有）。它
没有自适应步长控制（RK4 是固定步长的，并且只是近似保正），没有轨迹
求解器，并且 `lindblad_channel` 不接受扇区之间的跳跃；对于这些情形，以 $M_D$ 上的 QuTiP 或 dynamiqs
为参照。

## 习题

1. （Heisenberg 绘景）对 8.3 节的生成元证明 $\operatorname{Tr}(\mathcal L(\rho)A)=\operatorname{Tr}(\rho\,\mathcal L^*(A))$，其中
   $\mathcal L^*(A)=i[H,A]+\sum_jg_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})+\sum_mg_m(\Phi_m^*(A)-\frac12\{\Phi_m^*(1),A\})$，
   证明 $\mathcal L^*(1)=0$，并且若没有扇区之间的跳跃，则 $\mathcal L^*(P_c)=0$。
2. （例外点）在 8.4 节的扇区 $N=1$ 中，
   $H_{\mathrm{eff}}=\begin{pmatrix}-i\gamma_1/2&J\\J&-i\gamma_2/2\end{pmatrix}$。求它的特征值，并证明
   当 $J=|\gamma_1-\gamma_2|/4>0$ 时它不可对角化。
3. 设 $m$ 个模式中的 $N_0$ 个玻色子各自以速率 $\gamma$ 丢失（跳跃 $a_1,\dots,a_m$ 限制在
   $N\to N-1$ 上）。证明对 $\mathcal M$ 中的每个哈密顿量，$p_N(t)=\binom{N_0}Nq^N(1-q)^{N_0-N}$，其中 $q=e^{-\gamma t}$。
4. （纯退相位）在 $M_2$ 上令 $H=\omega|e\rangle\langle e|$，$L=\sigma_z=\operatorname{diag}(1,-1)$，速率为
   $\gamma$。证明布居数不变，$\rho_{ge}(t)=e^{(i\omega-2\gamma)t}\rho_{ge}(0)$，并且尽管 $C=1$，
   稳态仍张成由对角矩阵构成的二维空间。

<details>
<summary>解答</summary>

1. 由循环性，$\operatorname{Tr}(L\rho L^*A)=\operatorname{Tr}(\rho L^*AL)$，
   $\operatorname{Tr}(\{L^*L,\rho\}A)=\operatorname{Tr}(\rho\{L^*L,A\})$，且
   $\operatorname{Tr}(-i[H,\rho]A)=\operatorname{Tr}(\rho\,i[H,A])$；对扇区之间的跳跃，
   $\operatorname{Tr}(\Phi(\rho)A)=\operatorname{Tr}(\rho\Phi^*(A))$ 就是 $\Phi^*$ 的定义
   （6.3 节）。于是 $\mathcal L^*(1)=\sum_jg_j(L_j^*L_j-L_j^*L_j)+\sum_mg_m(\Phi_m^*(1)-\Phi_m^*(1))=0$，
   这又一次是保迹性。没有扇区之间的跳跃时，$P_c$ 与 $H$ 及
   所有 $L_j$ 可交换，所以 $L_j^*P_cL_j=L_j^*L_jP_c=\frac12\{L_j^*L_j,P_c\}$，从而 $\mathcal L^*(P_c)=0$；于是
   $p_c'=\operatorname{Tr}(\rho\mathcal L^*(P_c))=0$，这是命题 8.9(1) 的第二个证明。

2. $H_{\mathrm{eff}}=-\frac i4(\gamma_1+\gamma_2)1+\begin{pmatrix}-i\delta&J\\J&i\delta\end{pmatrix}$，其中
   $\delta=(\gamma_1-\gamma_2)/4$；第二个矩阵的迹为 $0$，行列式为 $\delta^2-J^2$，所以
   特征值为 $-\frac i4(\gamma_1+\gamma_2)\pm\sqrt{J^2-\delta^2}$：当 $J<|\delta|$ 时为两个衰减速率，
   当 $J>|\delta|$ 时为一个速率和两个频率。当 $J=|\delta|>0$ 时，第二个矩阵 $Y$ 满足
   $Y^2=(J^2-\delta^2)1=0$ 但 $Y\ne0$，所以 $H_{\mathrm{eff}}$ 只有一个特征值
   $-\frac i4(\gamma_1+\gamma_2)$，却不是 $1$ 的倍数：它不可对角化。

```python
lam = np.linalg.eigvals(H_eff[1])
delta = (g1 - g2) / 4
root = np.sqrt(J ** 2 - delta ** 2 + 0j)
assert np.allclose(sorted(lam, key=np.real), sorted(-0.25j * (g1 + g2) + np.array([root, -root]), key=np.real))
H_ep = np.array([[-0.5j * g1, delta], [delta, -0.5j * g2]])         # J = delta: the exceptional point
Y = H_ep + 0.25j * (g1 + g2) * np.eye(2)                           # H_ep minus its only eigenvalue
assert np.abs(Y @ Y).max() < 1e-15 and np.linalg.matrix_rank(Y) == 1   # one Jordan block of size 2
```

3. 在扇区 $N$ 上 $\sum_ja_j^*a_j=N1$，并且唯一的扇区对是 $(N-1,N)$，所以命题
   8.9(2) 给出 $p_N'=-\gamma Np_N+\gamma(N+1)p_{N+1}$。由 $q'=-\gamma q$，二项式
   项的导数为 $-\gamma Np_N+\gamma(N_0-N)\binom{N_0}Nq^{N+1}(1-q)^{N_0-N-1}$，而
   $(N+1)\binom{N_0}{N+1}=(N_0-N)\binom{N_0}N$ 表明第二项为 $\gamma(N+1)p_{N+1}$。在 $t=0$
   时它是集中在 $N_0$ 处的点质量；唯一性由引理 8.1 得到。

4. $L\rho L^*-\frac12\{L^*L,\rho\}=\sigma_z\rho\sigma_z-\rho$ 把非对角元乘以 $-2$ 并
   消去对角元，且 $-i[H,\rho]_{ge}=i\omega\rho_{ge}$。所以 $\mathcal L(\sigma)=0$ 当且仅当 $\sigma_{ge}=0$（对
   $\gamma>0$），谱为 $\{0,0,-2\gamma\pm i\omega\}$。由于 $H$ 与 $\sigma_z$ 可交换，且
   $X\mapsto\sigma_zX\sigma_z$ 是对合，$e^{t\mathcal L}(\rho)=\frac{1+e^{-2\gamma t}}2U_t\rho U_t^*
   +\frac{1-e^{-2\gamma t}}2U_t\sigma_z\rho\sigma_zU_t^*$，其中 $U_t=e^{-itH}$：这是一个保单位信道，其不动
   点是与它的 Kraus 算子可交换的矩阵（定理 6.10），即对角矩阵。

```python
deph = dynamics.lindblad_superoperator(q2, H_at, [q2.from_blocks([torch.diag(torch.tensor([1.0, -1.0], dtype=f64))])],
                                       rates=[gamma])[0][0]
ev = torch.linalg.eigvals(deph)
assert int((ev.abs() < 1e-12).sum()) == 2
assert np.allclose(sorted(ev.tolist(), key=lambda z: (z.real, z.imag))[:2],
                   [-2 * gamma - 1j * omega, -2 * gamma + 1j * omega], atol=1e-12)
```

</details>

**下一章：** 第 9 章（编写中）讨论自旋链：磁化扇区与动量扇区，
大扇区的 Lanczos 与 Krylov 时间演化，以及无序。
