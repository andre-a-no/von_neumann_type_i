[English](../10_quantum_computing.md) · [Español](../es/10_quantum_computing.md) · [Français](../fr/10_quantum_computing.md) · [Deutsch](../de/10_quantum_computing.md) · **中文** · [日本語](../ja/10_quantum_computing.md) · [Русский](../ru/10_quantum_computing.md)

# 10. 量子计算

量子计算机上的噪声是一个信道，而信道生成一个代数：由其 Kraus 算子生成的 $*$-代数。这个代数的
结构（第 2 章）告诉我们哪些信息能毫发无损地经受住噪声；一个关于 Kraus 算子乘积、涉及另一个代数的
条件告诉我们恢复操作能够还原哪些信息。最后一部分讨论守恒某个荷的算法：它们的态处在荷扇区的代数中，
而振幅阻尼会使态在扇区之间移动。

**你将学到**为什么保单位信道与其对偶具有相同的不动点（噪声的交换子），带完整证明的
Knill–Laflamme 定理，如何只保护扇区的标记，以及量子交替算子拟设的 XY 混合算子、振幅阻尼和后选择
在汉明重量扇区的代数 $\bigoplus_wM_{\binom nw}$ 中是什么样子。本章的数值检验见
[教程 10](../../tutorial/zh/10_quantum_computing.md)。

在本章中，$M_D(\mathbb C)$ 上的**信道**是映射 $\Phi(\rho)=\sum_{i=1}^rE_i\rho E_i^*$，其
**Kraus 算子**满足 $\sum_iE_i^*E_i=1$（保迹）；若此外还有
$\sum_iE_iE_i^*=1$，则称它是**保单位的**。它的**对偶** $\Phi^*(T)=\sum_iE_i^*TE_i$ 满足
$\mathrm{Tr}(\Phi(\rho)T)=\mathrm{Tr}(\rho\,\Phi^*(T))$（6.2 节）。（由定理 6.2，每个完全
正映射都具有这种形式；这里信道是*以* Kraus 算子*给出*的。）$\Phi$ 的
**相互作用代数** $\mathcal A$ 是由 $1$ 和各 $E_i$ 生成的 $*$-代数，即 $E_i$ 与 $E_i^*$ 的所有乘积的
线性张成；一个矩阵与 $\mathcal A$ 可交换当且仅当它与所有
$E_i$ 和 $E_i^*$ 可交换。若 $\rho=P_K\rho P_K$，则称态 $\rho$ **支撑在**投影为 $P_K$ 的子空间 $K$ 中。
在 $n$ 个量子比特上，$D=2^n$；$X,Y,Z$ 是 Pauli 矩阵，$Z=\mathrm{diag}(1,-1)$，$X_j$ 是作用在量子比特 $j$ 上的 $X$（编号为 $1,\dots,n$；
10.4 节沿用本库的约定，格点编号为 $0,\dots,n-1$），$|b\rangle$（$b\in\{0,1\}^n$）是计算基，
$|b|$ 是**汉明重量**，$P_w$ 是到 $H_w=\mathrm{span}\{|b\rangle:|b|=w\}$ 上的投影。

<a id="s-10-1"></a>

## 10.1 噪声代数、无噪声子系统与不动点

由定理 2.6，存在酉矩阵 $V$，使得 $V^*\mathcal AV=\bigoplus_cM_{k_c}\otimes1_{m_c}$ 且
$V^*\mathcal A'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$：每个 Kraus 算子都只作用在每个直和项的第一个
张量因子上。

<a id="r-10-1"></a>

**命题 10.1（无噪声子系统；无退相干子空间）。** (a) 记
$V^*E_iV=\bigoplus_cE_{i,c}\otimes1_{m_c}$，并令 $\rho=V(\sigma\otimes\tau)V^*$ 位于直和项
$c$ 中，其中 $\sigma$ 是 $\mathbb C^{k_c}$ 上的态，$\tau$ 是 $\mathbb C^{m_c}$ 上的态。则
$\Phi(\rho)=V(\Phi_c(\sigma)\otimes\tau)V^*$，其中 $\Phi_c(\sigma)=\sum_iE_{i,c}\sigma E_{i,c}^*$ 是
$M_{k_c}$ 上的信道：态 $\tau$ 丝毫不受影响。此外，对每个 $T\in\mathcal A'$ 有 $\Phi^*(T)=T$。
(b) 若子空间 $K$ 对所有 $i$ 满足 $E_iP_K=e_iP_K$，$e_i\in\mathbb C$，则对每个支撑在 $K$ 中的
态 $\rho$ 有 $\Phi(\rho)=\rho$。

*证明。* (a) $\Phi(\rho)=V\sum_i(E_{i,c}\sigma E_{i,c}^*\otimes\tau)V^*$，且 $\Phi_c$ 是保迹的，
因为 $V^*(\sum_iE_i^*E_i)V=1$ 的第 $c$ 个直和项是 $(\sum_iE_{i,c}^*E_{i,c})\otimes1_{m_c}$。对
$T\in\mathcal A'$，$\Phi^*(T)=T\sum_iE_i^*E_i=T$。(b) $P_KE_i^*=(E_iP_K)^*=\bar e_iP_K$，所以
$E_i\rho E_i^*=E_iP_K\rho P_KE_i^*=|e_i|^2\rho$，且 $\sum_i|e_i|^2P_K=P_K\sum_iE_i^*E_iP_K=P_K$。$\square$

$k_c=1$ 的直和项就是情形 (b)：**无退相干子空间**（Zanardi–Rasetti 1997，
Lidar–Chuang–Whaley 1998）。当 $k_c>1$ 时，因子 $\mathbb C^{m_c}$ 是一个**无噪声子系统**
（Knill–Laflamme–Viola 2000）：噪声确实在起作用，但只作用在另一个因子上。命题 10.1 对每个信道
都成立；对保单位信道，被保持的可观测量*恰好*是 $\mathcal A'$ 中的那些：

<a id="r-10-2"></a>

**定理 10.2（保单位信道的不动点）。** 设 $\Phi$ 是保单位的。则
$$
\mathrm{Fix}(\Phi)=\{T\in M_D(\mathbb C):\Phi(T)=T\}=\mathcal A'=\mathrm{Fix}(\Phi^*).
$$
特别地，$\mathrm{Fix}(\Phi)$ 是一个冯·诺依曼代数，并且 $*$-同构于 $\bigoplus_cM_{m_c}$。

*证明。* 第一个等式就是定理 6.10 应用于代数 $M_D(\mathbb C)$（只有一个块）的情形：$T$ 是不动点
当且仅当它与所有 $E_i$ 和 $E_i^*$ 可交换，即当且仅当 $T\in\mathcal A'$。对偶 $\Phi^*$ 也是一个
保单位信道：它的 Kraus 算子 $E_i^*$ 满足 $\sum_iE_iE_i^*=1$（保迹）和
$\sum_iE_i^*E_i=1$（保单位），并且生成同一个 $*$-代数 $\mathcal A$。因此对 $\Phi^*$ 应用定理 6.10
得到 $\mathrm{Fix}(\Phi^*)=\mathcal A'$。最后一个结论就是定理 2.6。$\square$

*计算：* 集体退相位的 $\mathrm{Fix}(\Phi)=\mathcal A'$，其中交换子与不动点均作为零空间计算——[教程，10.1 节](../../tutorial/zh/10_quantum_computing.md#s-10-1)。

对每个信道都有 $\mathcal A'\subseteq\mathrm{Fix}(\Phi^*)$（命题 10.1(a)），但没有
保单位性，该定理就不成立：对振幅阻尼，$\mathrm{Fix}(\Phi)=\mathbb C|0\rangle\langle0|$，而
$\mathcal A'=\mathrm{Fix}(\Phi^*)=\mathbb C1$（6.4 节与习题 1）。（历史上：Kribs 2003，
Arias–Gheondea–Gudder 2002。）

<a id="r-10-3"></a>

**推论 10.3（集体退相位）。** 设 $U=e^{-i\theta Z_{\mathrm{tot}}}$，
$Z_{\mathrm{tot}}=\sum_jZ_j$，其中对 $d=1,\dots,n$ 有 $e^{2i\theta d}\neq1$，并设
$\Phi(\rho)=\tfrac12(U\rho U^*+U^*\rho U)$。则 $\mathcal A=\mathrm{span}\{P_0,\dots,P_n\}$，每个 $H_w$
都是无退相干子空间，并且 $\mathrm{Fix}(\Phi)=\bigoplus_wB(H_w)\cong\bigoplus_wM_{\binom nw}$，
其维数为 $\sum_w\binom nw^2$。

*证明。* $Z_{\mathrm{tot}}|b\rangle=(n-2|b|)|b\rangle$，所以 $U=\sum_we^{-i\theta(n-2w)}P_w$，其
特征值互不相同（它们的比为 $e^{2i\theta(w-w')}$，$0<|w-w'|\le n$）。定理 1.5(a) 证明中的
Lagrange 多项式给出 $P_w=p_w(U)\in\mathcal A$；反之，$U,U^*$ 属于代数
$\mathrm{span}\{P_w\}$。所以 $\mathcal A=\mathrm{span}\{P_w\}$（$k_c=1$，$m_c=\binom nw$），
命题 10.1(b) 适用于每个 $H_w$。若 $T$ 与所有 $P_w$ 可交换，则当
$v\neq v'$ 时 $P_vTP_{v'}=TP_vP_{v'}=0$，即 $T$ 是块对角的，反之亦然；其余部分由定理 10.2 给出。$\square$

*计算：* $n=4$ 时的维数 $5$、$70$、$70$——[教程，10.1 节](../../tutorial/zh/10_quantum_computing.md#s-10-1)。*库中：* $\bigoplus_wM_{\binom nw}$ 就是 `SpinChain(n).algebra`，该信道在其上是恒等映射——[教程，10.1 节](../../tutorial/zh/10_quantum_computing.md#s-10-1)。

> **物理意义**（物理解释，而非本课程的定理）。一个在所有量子比特上都相同的涨落场使 $|b\rangle$ 的相位按 $n-2|b|$ 的比例
> 移动。同一汉明重量的态获得共同的相位，因此 $|01\rangle,|10\rangle$
> 编码了一个对这种噪声免疫的量子比特（离子阱：Kielpinski et al., *Science* 291, 2001），而
> 不同重量之间的相干性会衰减：噪声施加了一个超选择定则（3.3 节）。

**集体退相干。** 若一个信道的 Kraus 算子都是转动
$e^{-i\theta\,\mathbf n\cdot\mathbf S}$（$\theta\in\mathbb R$，$\mathbf n\in\mathbb R^3$ 为单位向量）的倍数，其中
$S_\alpha=\tfrac12\sum_j\sigma_{\alpha,j}$，$\sigma\in\{X,Y,Z\}$，我们就称该信道是**集体的**；它模拟方向也在涨落的场
（这是模型，而不是定理）。由于 $\mathbf n\cdot\mathbf S$ 是自伴的，$e^{-i\theta\,\mathbf n\cdot\mathbf S}$ 是
$\mathbf n\cdot\mathbf S$ 的多项式（定理 1.5(a) 的证明），所以这些 Kraus 算子属于由
$1,S_x,S_y,S_z$ 生成的 $*$-代数 $\mathcal B$。

<a id="r-10-4"></a>

**命题 10.4（单态子空间）。** 设 $K=\ker S_x\cap\ker S_y\cap\ker S_z$。每个 Kraus 算子属于
$\mathcal B$ 的信道都固定每个支撑在 $K$ 中的态。

*证明。* $\mathcal B$ 由 $S_x,S_y,S_z$ 的字（以及 $1$）张成；长度 $\ge1$ 的字把 $K$ 映为零，
所以 $\mathcal B$ 的每个元素在 $K$ 上都作用为一个数乘。应用命题 10.1(b)。$\square$

对 $n=4$，命题 4.13（[第 4 章](04_symmetries.md)，$L=4$ 的情形）给出
$\mathcal B\cong M_5\oplus M_3\oplus M_1$，重数为 $1,3,2$（总自旋 $2,1,0$），所以 $\dim\mathcal B=35$，而
$\mathcal B'\cong M_1\oplus M_3\oplus M_2$ 的维数为 $14$。此外 $\dim K=m_0=2$：各 $S_a$ 逐个直和项地作用，
所以 $K$ 是它与各直和项之交的和；在自旋 $S$ 的直和项上，它们作用为 $S^{(S)}_a\otimes1$（推论 4.12），而
$\sum_bu_b\otimes f_b$（$f_b$ 是 $\mathbb C^{m_S}$ 的一组基）被它们全部映为零，当且仅当每个 $u_b$ 都被
所有 $S^{(S)}_a$ 映为零。这样的 $u$ 构成不可约自旋 $S$ 表示的一个不变子空间（定理 4.11），当 $S>0$ 时它不是整个空间
（因为 $S^{(S)}_z\neq0$），因此为 $0$，而当 $S=0$ 时它是 $\mathbb C^1$。命题 10.1(a) 及其证明
对任何包含 Kraus 算子、形如 $V(\bigoplus_cM_{k_c}\otimes1_{m_c})V^*$ 的 $*$-代数逐字成立，这里即
$\mathcal B$：对每个 Kraus 算子属于 $\mathcal B$ 的信道，四个量子比特承载一个无退相干的量子比特
（$k=1$，$m=2$）和一个无噪声的三能级系统（qutrit，$k=3$，$m=3$）。

*计算：* $\dim\mathcal B'=14$、$\dim\mathcal B=35$、$\dim K=2$，以及针对随机集体转动的命题 10.4——[教程，10.1 节](../../tutorial/zh/10_quantum_computing.md#s-10-1)。

<a id="s-10-2"></a>

## 10.2 纠错：Knill–Laflamme 定理

**码**是一个投影为 $P$ 的非零子空间 $\mathcal C\subseteq\mathbb C^D$。若存在信道 $\mathcal R$
（**恢复**），使得对每个支撑在 $\mathcal C$ 中的态 $\rho$ 有 $\mathcal R(\Phi(\rho))=\rho$，则称
$\mathcal C$ 能**纠正** $\Phi$。

<a id="r-10-5"></a>

**定理 10.5（Knill–Laflamme）。** $\mathcal C$ 能纠正 $\Phi$ 当且仅当存在数 $\lambda_{ij}$，
使得
$$
PE_i^*E_jP=\lambda_{ij}P\qquad\text{for all }i,j.
$$

*证明。* *必要性。* 设 $\mathcal R$ 的 Kraus 算子为 $R_m$；则 $\mathcal R\circ\Phi$ 的
Kraus 算子为 $R_mE_i$。对单位向量 $\psi\in\mathcal C$，
$\sum_{m,i}R_mE_i|\psi\rangle\langle\psi|E_i^*R_m^*=|\psi\rangle\langle\psi|$，而对 $\phi\perp\psi$，
取在 $\phi$ 中的期望值得到 $\sum_{m,i}|\langle\phi,R_mE_i\psi\rangle|^2=0$，所以 $R_mE_i\psi\in\mathbb C\psi$。
因此 $R_mE_iP=c_{mi}P$：对所有 $x$ 都满足 $Tx\in\mathbb Cx$ 的线性映射 $T$ 是 $1$ 的倍数（若
$Tx=ax$，$Ty=by$，$T(x+y)=c(x+y)$，且 $x,y$ 线性无关，则 $a=c=b$）。由于 $\sum_mR_m^*R_m=1$，
$PE_i^*E_jP=\sum_m(R_mE_iP)^*(R_mE_jP)=\sum_m\bar c_{mi}c_{mj}P$。

*充分性。* $\Lambda=(\lambda_{ij})$ 是 Hermite 的，因为 $\lambda_{ji}P=(PE_i^*E_jP)^*=\bar\lambda_{ij}P$；
它也是正的，因为对 $F=\sum_jc_jE_j$ 有 $\sum_{ij}\bar c_ic_j\lambda_{ij}P=PF^*FP\ge0$。写 $\Lambda=UDU^*$，
其中 $U$ 为酉矩阵，$D=\mathrm{diag}(d_k)$，$d_k\ge0$（定理 1.3），并令 $F_k=\sum_iu_{ik}E_i$。则
$\sum_kF_k\rho F_k^*=\sum_{ij}(UU^*)_{ij}E_i\rho E_j^*=\Phi(\rho)$，且 $PF_k^*F_lP=(U^*\Lambda U)_{kl}P=d_k\delta_{kl}P$。
若 $d_k=0$，则对所有 $x$ 有 $\|F_kPx\|^2=0$，所以 $F_kP=0$。若 $d_k>0$，令 $V_k=F_kP/\sqrt{d_k}$；则
$V_k^*V_l=\delta_{kl}P$，所以 $V_kP=V_k$，各 $Q_k=V_kV_k^*$ 是两两正交的投影，且
$Q=\sum_kQ_k$ 是投影。恢复
$$
\mathcal R(\sigma)=\sum_{k:d_k>0}V_k^*\sigma V_k+(1-Q)\sigma(1-Q)
$$
是一个信道，因为 $\sum_kV_kV_k^*+(1-Q)^2=1$。对 $\rho=P\rho P$，$\Phi(\rho)=\sum_{d_k>0}d_kV_k\rho V_k^*$；
由于 $V_l^*V_k=\delta_{kl}P$ 且 $(1-Q)V_k=(1-Q)Q_kV_k=0$，有 $\mathcal R(\Phi(\rho))=\sum_kd_k\rho$，而
$\sum_kd_k=\mathrm{Tr}\,\Lambda=1$，因为 $\sum_i\lambda_{ii}P=P\sum_iE_i^*E_iP=P$。$\square$

（Knill–Laflamme 1997；Nielsen–Chuang，定理 10.1。）$\mathcal R$ 测量**综合征**——即两两正交的
子空间 $Q_k\mathbb C^D$ 中哪一个包含该态——然后把态转回码中。

**例：比特翻转码。** $\mathcal C=\mathrm{span}\{|000\rangle,|111\rangle\}$，噪声为"至多一个比特
翻转"，$E_0=\sqrt{p_0}\,1$，$E_j=\sqrt{p_j}\,X_j$。当 $i\neq j$ 时，$X_iX_j$（其中 $X_0=1$）翻转一个或两个
比特，并把 $\mathcal C$ 映入 $\mathcal C^\perp$；所以 $\Lambda=\mathrm{diag}(p_0,\dots,p_3)$，$V_k=X_kP$，
$Q=1$，且 $\mathcal R(\sigma)=\sum_{k=0}^3PX_k\sigma X_kP$。若每个量子比特以概率 $p$ 独立地
翻转，则该条件不成立（Kraus 算子中包含 $X_1X_2$ 和 $X_3$ 的倍数，而
$PX_1X_2X_3P=|000\rangle\langle111|+|111\rangle\langle000|$ 不是 $P$ 的倍数）；当发生 $\ge2$ 次翻转时，$\mathcal R$ 把 $|0_L\rangle=|000\rangle$ 映为 $|1_L\rangle$，
否则保持不变；因此恢复后 $|0_L\rangle$ 的保真度就是
至多发生一次翻转的概率，$p=0.1$ 时 $(1-p)^3+3p(1-p)^2=1-3p^2+2p^3=0.972$。

*计算：* 在这个例子上检验定理 10.5：矩阵 $\Lambda$、恢复 $\mathcal R$ 与保真度 $0.972$——[教程，10.2 节](../../tutorial/zh/10_quantum_computing.md#s-10-2)。*库中：* 由 Kraus 算子给出的信道 `channels.kraus_channel` 及其复合 `R @ Phi`——[教程，10.2 节](../../tutorial/zh/10_quantum_computing.md#s-10-2)。

<a id="s-10-3"></a>

## 10.3 算子代数纠错

Heisenberg 绘景使我们可以只保护比整个码更少的东西。设 $\mathfrak A$ 是 $\mathcal C=P\mathbb C^D$ 上的
冯·诺依曼代数（$\mathcal C$ 上包含 $P$ 的算子 $*$-代数），即需要保护的可观测量。若存在信道
$\mathcal R$，使得
$$
P\,(\mathcal R\circ\Phi)^*(A)\,P=A\qquad\text{for all }A\in\mathfrak A,
$$
即对所有支撑在 $\mathcal C$ 中的态 $\rho$
和所有 $A\in\mathfrak A$ 有 $\mathrm{Tr}(\mathcal R(\Phi(\rho))A)=\mathrm{Tr}(\rho A)$（左边等于
$\mathrm{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(A)P)$，而对所有这些 $\rho$ 都满足 $\mathrm{Tr}(\rho Y)=0$ 的算子
$Y\in B(\mathcal C)$ 必为零：对所有 $x$ 有 $\langle x,Yx\rangle=\langle Px,YPx\rangle=0$，再像推论 1.4 的证明那样
使用极化恒等式），则称 $\mathfrak A$ 对 $\Phi$ 是**可纠正的**。对 $\mathfrak A=B(\mathcal C)$，这就是 10.2 节意义下的纠错：若
$\sigma=\mathcal R(\Phi(\rho))$ 满足 $P\sigma P=\rho$（取所有 $A\in B(\mathcal C)$），则由于 $\mathcal R\circ\Phi$ 保迹，
$\mathrm{Tr}(1-P)\sigma(1-P)=\mathrm{Tr}\,\sigma-\mathrm{Tr}\,P\sigma P=1-1=0$，所以 $\sigma^{1/2}(1-P)=0$，从而
$\sigma=\rho$；反之，由 $\mathcal R(\Phi(\rho))=\rho$ 可得迹的等式。另一个极端是交换代数：只有一个经典标记受到保护。

<a id="r-10-6"></a>

**命题 10.6（纠正扇区标记）。** 设 $P=\sum_{c=1}^CP_c$，其中 $P_c\neq0$ 是两两正交的
投影，并设 $\mathfrak A=\mathrm{span}\{P_1,\dots,P_C\}$。以下各条等价：
1. 存在 $N_1,\dots,N_C\ge0$，$\sum_cN_c=1$，使得对每个 $c$ 和
   每个支撑在 $P_c\mathbb C^D$ 中的态 $\rho$ 有 $\mathrm{Tr}(N_c\Phi(\rho))=1$；
2. 对所有 $i,j$ 及 $c\neq d$ 有 $P_cE_i^*E_jP_d=0$，即 $PE_i^*E_jP\in\mathfrak A'$，这里 $\mathfrak A'$ 是
   $\mathfrak A$ 在 $B(\mathcal C)$ 中的交换子；
3. $\mathfrak A$ 对 $\Phi$ 是可纠正的。

*证明。* $T\in B(\mathcal C)$ 与所有 $P_c$ 可交换当且仅当对 $c\neq d$ 有 $P_cTP_d=0$（同推论 10.3）。

(2)$\Rightarrow$(1)。子空间 $S_c=\mathrm{span}\{E_iP_cx\}$ 两两正交，因为
$\langle E_iP_cx,E_jP_dy\rangle=\langle x,P_cE_i^*E_jP_dy\rangle=0$。对 $c<C$，令 $N_c$ 为到 $S_c$ 上的投影，
并令 $N_C=1-\sum_{c<C}N_c$，它也是投影。对支撑在 $P_c\mathbb C^D$ 中的 $\rho$，
$\Phi(\rho)=\sum_iE_iP_c\rho P_cE_i^*$ 的值域含于 $S_c$，而 $N_c$ 在 $S_c$ 上是恒等映射（对 $c=C$，
这是因为 $S_C\perp S_{c'}$，$c'<C$），所以 $\mathrm{Tr}(N_c\Phi(\rho))=\mathrm{Tr}\,\Phi(\rho)=1$。
(1)$\Rightarrow$(2)。对单位向量 $x\in P_c\mathbb C^D$ 和 $d\neq c$，诸数
$\mathrm{Tr}(N_{d'}\Phi(|x\rangle\langle x|))\ge0$ 之和为 $1$，且 $d'=c$ 对应的那个数等于 $1$，所以
$0=\mathrm{Tr}(N_d\Phi(|x\rangle\langle x|))=\sum_i\|N_d^{1/2}E_ix\|^2$；因此 $N_dE_iP_c=0$，再对
$d$ 求和得 $N_cE_iP_c=E_iP_c$。从而当 $c\neq d$ 时 $P_cE_i^*E_jP_d=P_cE_i^*N_cE_jP_d=0$。
(1),(2)$\Rightarrow$(3)。令 $\omega_c=P_c/\mathrm{Tr}\,P_c$，$\mathcal R(\sigma)=\sum_c\mathrm{Tr}(N_c\sigma)\,\omega_c$
（"先测量，再制备"），这是一个信道，其 Kraus 算子为 $\sqrt\mu\,|u\rangle\langle f|N_c^{1/2}$，其中
$(\mu,u)$ 取遍 $\omega_c$ 的特征对，$f$ 取遍一组规范正交基中的向量。由于
$\mathrm{Tr}(\mathcal R(\sigma)P_d)=\mathrm{Tr}(N_d\sigma)$，有 $\mathcal R^*(P_d)=N_d$。利用 $N_dE_iP_c=\delta_{cd}E_iP_c$
（由上一步得到），再利用 $\sum_iE_i^*E_i=1$，
$$
P(\mathcal R\circ\Phi)^*(P_d)P=P\,\Phi^*(N_d)\,P=\sum_i\sum_{c,c'}P_{c'}E_i^*N_dE_iP_c=\sum_i\sum_{c'}P_{c'}E_i^*E_iP_d=P_d\Big(\sum_iE_i^*E_i\Big)P_d=P_d .
$$
(3)$\Rightarrow$(1)。对 $c<C$ 令 $N_c=\mathcal R^*(P_c)$，并令 $N_C=\mathcal R^*(1-\sum_{c<C}P_c)$；它们是
正的，且其和为 $\mathcal R^*(1)=1$。对支撑在 $P_c\mathbb C^D$ 中的 $\rho$ 和 $d<C$，
$\mathrm{Tr}(N_d\Phi(\rho))=\mathrm{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(P_d)P)=\mathrm{Tr}(\rho P_d)=\delta_{cd}$，因此也有
$\mathrm{Tr}(N_C\Phi(\rho))=\delta_{cC}$。$\square$

<a id="r-10-7"></a>

**定理 10.7（Bény–Kempf–Kribs）。** $\mathcal C$ 上的冯·诺依曼代数 $\mathfrak A$
对 $\Phi$ 是可纠正的，当且仅当对所有 $i,j$ 有 $PE_i^*E_jP\in\mathfrak A'$，其中
$\mathfrak A'$ 是 $\mathfrak A$ 在 $B(\mathcal C)$ 中的交换子。

*证明。* 这就是定理 D.9，在[附录 D](D_gksl.md)中证明；其证明还构造了一个恢复信道：它测量
$\mathfrak A\cong\bigoplus_cM_{k_c}$ 的扇区和一个综合征，并像定理 10.5 的证明那样撤销每个因子
$\mathbb C^{k_c}$ 上的错误。$\square$

对 $\mathfrak A=B(\mathcal C)$，$\mathfrak A'=\mathbb CP$，这就是定理 10.5；命题 10.6 是
交换代数 $\mathfrak A=\mathrm{span}\{P_c\}$ 的情形（推论 D.10(a),(b)）。对
$\mathcal C=\mathbb C^k\otimes\mathbb C^m$ 上的 $\mathfrak A=M_k\otimes1_m$，得到算子（子系统）
码；10.1 节的无噪声子系统就是这样的码，由 $\mathcal R=\mathrm{id}$ 纠正
（推论 D.10(c)，交换两个张量因子）。由定理 2.6，每个 $\mathfrak A$ 都
$*$-同构于某个 $\bigoplus_cM_{k_c}$；当 $C\ge2$ 且某个 $k_c\ge2$ 时，可纠正的 $\mathfrak A$ 既保护扇区标记 $c$，
又保护每个 $M_{k_c}$ 的量子信息——这就是**混合**码（见推论 D.10 之后的说明，那里有一个
Knill–Laflamme 条件不成立而定理 10.7 适用的例子）。
历史上：C. Bény, A. Kempf, D. W. Kribs, *Phys. Rev. Lett.* 98, 100502 (2007) 以及 *Phys. Rev. A* 76,
042303 (2007)；算子码：Kribs–Laflamme–Poulin, *Phys. Rev. Lett.* 94, 180501 (2005)。

**例。** 重复码先受到独立的相位翻转（每个量子比特的概率为 $q$），再受到
至多一个比特翻转。Kraus 算子是 $X_kZ^s$（$s\in\{0,1\}^3$）的倍数，每个
$E_i^*E_j$ 至多翻转两个比特，所以条件 2 对 $P_0=|000\rangle\langle000|$、
$P_1=|111\rangle\langle111|$ 成立：经典比特得以保存，可以通过多数表决读出，
$N_0=\sum_{|b|\le1}|b\rangle\langle b|$，$N_1=1-N_0$。量子比特则不然：$PZ_1P$ 不是 $P$ 的倍数，
并且在 10.2 节的恢复 $\mathcal R$ 之后，态 $|{+_L}\rangle=(|000\rangle+|111\rangle)/\sqrt2$ 恰好在发生奇数次
相位翻转时遭受逻辑相位翻转 $Z_1Z_2Z_3$；其保真度为
$(1+(1-2q)^3)/2=0.8645$（$q=0.05$ 时）。（该信道是 $\Phi_{\mathrm{one}}\circ\mathcal Z$，其中 $\Phi_{\mathrm{one}}$ 是
比特翻转例子中"至多一个比特翻转"的信道，$\mathcal Z$ 是独立的相位翻转；
在 $\mathcal C$ 上每个 $Z^s$ 都作用为 $Z_1^{|s|}$，所以 $\mathcal Z(\rho)=\sum_sq^{|s|}(1-q)^{3-|s|}Z_1^{|s|}\rho Z_1^{|s|}$
仍支撑在 $\mathcal C$ 中，而由定理 10.5，$\mathcal R\circ\Phi_{\mathrm{one}}$ 在这样的态上是恒等映射。
$Z_1$ 把 $|{+_L}\rangle$ 映为 $|{-_L}\rangle\perp|{+_L}\rangle$，而 $|s|$ 为偶数的概率是
$\frac12\big((1-q+q)^3+(1-q-q)^3\big)$。）

*计算：* 在这个例子上检验命题 10.6：条件 2、多数表决 $N_0,N_1$ 以及 $|{+_L}\rangle$ 的保真度 $0.8645$——[教程，10.3 节](../../tutorial/zh/10_quantum_computing.md#s-10-3)。

> **物理意义**（物理解释，而非本课程的定理）。保护一个交换代数意味着保护经典信息（一个
> 比特、一个超选择扇区的标记），这比保护一个量子比特便宜得多：与扇区投影
> 可交换的噪声是无害的（此时对 $c\neq d$ 有 $P_cE_i^*E_jP_d=E_i^*E_jP_cP_d=0$，即命题 10.6 的条件 2）。

<a id="s-10-4"></a>

## 10.4 带约束的混合算子

量子交替算子拟设（Hadfield et al., *Algorithms* 12, 34, 2019）制备
$|\beta,\gamma\rangle=e^{-i\beta_pH_M}e^{-i\gamma_pH_C}\cdots e^{-i\beta_1H_M}e^{-i\gamma_1H_C}|b_0\rangle$，其中 $H_C$ 是对角的
代价哈密顿量，$H_M$ 是**混合算子**。对于约束 $|b|=w$（图二分、由 $w$ 种资产构成的投资
组合），取 $|b_0|=w$，并在环上取 **XY 混合算子** $H_M=\tfrac14\sum_{(i,j)}(X_iX_j+Y_iY_j)$。

<a id="r-10-8"></a>

**命题 10.8.** $H_M$、$H_C$、$e^{-i\beta H_M}$ 和 $e^{-i\gamma H_C}$ 都属于
$\mathcal M=\{P_0,\dots,P_n\}'=\bigoplus_wB(H_w)\cong\bigoplus_wM_{\binom nw}$，因此 $|\beta,\gamma\rangle\in H_w$：
拟设永远不会离开可行子空间。

*证明。* 在量子比特 $i,j$ 上，$X\otimes X+Y\otimes Y=2(|01\rangle\langle10|+|10\rangle\langle01|)$（利用 $Y|0\rangle=i|1\rangle$，
$Y|1\rangle=-i|0\rangle$）：当比特 $i,j$ 不同时，$X_iX_j+Y_iY_j$ 交换它们，否则把 $|b\rangle$
映为零，因此它把每个 $H_w$ 映入自身。对角的 $H_C$ 与对角的 $P_w$ 可交换，
指数函数也是如此（定理 1.5(d)）。$\square$

*库中：* $\mathcal M$ 就是 `SpinChain(n).algebra`，XY 混合算子就是取 $\Delta=0$ 的 `SpinChain.xxz`——[教程，10.4 节](../../tutorial/zh/10_quantum_computing.md#s-10-4)。*计算：* 六个量子比特上图二分问题的命题 10.8，分别用本库和稠密矩阵计算——[教程，10.4 节](../../tutorial/zh/10_quantum_computing.md#s-10-4)。

量子比特 $i$ 的**振幅阻尼**（以概率 $g$ 发生 $|1\rangle\to|0\rangle$）具有 Kraus 算子
$K_0^{(i)}=|0\rangle\langle0|_i+\sqrt{1-g}\,|1\rangle\langle1|_i$ 和 $K_1^{(i)}=\sqrt g\,|0\rangle\langle1|_i$，后者降低
重量：这是扇区之间的映射。

*库中：* `SpinChain.site_amplitude_damping`，一个 `InterSectorChannel`——[教程，10.4 节](../../tutorial/zh/10_quantum_computing.md#s-10-4)。

<a id="r-10-9"></a>

**命题 10.9（泄漏与后选择）。** 设 $\rho_0$ 支撑在 $H_w$ 中，$U_l$ 是
$\mathcal M$ 中的酉算子，$\mathcal D$ 是对每个量子比特的阻尼，并令 $\rho_l=\mathcal D(U_l\rho_{l-1}U_l^*)$。
则每个 $\rho_l$ 都是块对角的，其重量 $\le w$，并且
$$
P_w\rho_pP_w=(1-g)^{wp}\;U\rho_0U^*,\qquad U=U_p\cdots U_1 :
$$
即以概率 $(1-g)^{wp}$ 测得重量 $w$，而后选择得到的态就是无噪声的态。

*证明。* $\mathcal D$ 的 Kraus 算子为 $K_s=\prod_iK^{(i)}_{s_i}$，$s\in\{0,1\}^n$；$K_s$ 把
$H_v$ 映入 $H_{v-|s|}$，即 $P_vK_s=K_sP_{v+|s|}$（当 $v>n$ 时 $P_v=0$）。若 $\sigma$ 是重量
$\le w$ 的块对角算子，则 $P_vK_s\sigma K_s^*P_{v'}=K_sP_{v+|s|}\sigma P_{v'+|s|}K_s^*$ 在 $v\neq v'$ 时以及
$v>w-|s|$ 时为零：$\mathcal D(\sigma)$ 是重量 $\le w$ 的块对角算子，且其重量为 $w$ 的块
只来自 $s=0$。由于 $K_0|b\rangle=(1-g)^{|b|/2}|b\rangle$，$P_w\mathcal D(\sigma)P_w=(1-g)^wP_w\sigma P_w$。各 $U_l$
与 $P_v$ 可交换；对 $l$ 归纳即可。$\square$

*计算：* 在二分线路上检验命题 10.9，并用 `DensityMatrix.condition_on` 进行后选择——[教程，10.4 节](../../tutorial/zh/10_quantum_computing.md#s-10-4)。

**连续时间。** 同样的结论对 Lindblad 方程（定理 8.6，[第 8 章](08_dynamics.md)）也成立：
$$
\dot\rho=-i[H,\rho]+\gamma_1\sum_i\Big(S^-_i\rho S^+_i-\tfrac12\{S^+_iS^-_i,\rho\}\Big)+\gamma_2\sum_k\Big(L_k\rho L_k^*-\tfrac12\{L_k^*L_k,\rho\}\Big)
$$
其中 $H\in\mathcal M$，扇区内部的跳跃 $L_k\in\mathcal M$，例如纯退相位 $L_k=S^z_k$。
这里 $\sum_iS^+_iS^-_i=\sum_i|1\rangle\langle1|_i=N$ 是**数算子**，$NP_v=vP_v$：它是中心元。

<a id="r-10-10"></a>

**命题 10.10.** 若 $\rho(0)$ 支撑在 $H_w$ 中，则 $P_w\rho(t)P_w=e^{-\gamma_1wt}\sigma(t)$，
其中 $\sigma$ 是同一方程在 $\gamma_1=0$、$\sigma(0)=\rho(0)$ 时的解。因此接受概率为
$e^{-\gamma_1wt}$，而后选择的态 $\sigma(t)$ 不依赖于 $\gamma_1$。

*证明。* 令 $B_{vv'}=P_v\rho P_{v'}$。由于 $H,L_k,N$ 与 $P_v$ 可交换且 $P_vS^-_i=S^-_iP_{v+1}$，
$\dot B_{vv'}=-i[H,B_{vv'}]+\gamma_1\sum_iS^-_iB_{v+1,v'+1}S^+_i-\tfrac{\gamma_1}2(v+v')B_{vv'}+\gamma_2\mathcal L_2(B_{vv'})$，
其中 $\mathcal L_2$ 是退相位部分。$\max(v,v')>w$ 的块满足一个初值为零的封闭线性方程组，
所以它们为零（引理 8.1 中的唯一性）。于是 $\dot B_{ww}=-i[H,B_{ww}]-\gamma_1wB_{ww}+\gamma_2\mathcal L_2(B_{ww})$，所以
$e^{\gamma_1wt}B_{ww}$ 是 $\sigma$ 所满足的方程的解，且初值相同，从而等于 $\sigma$
（引理 8.1）。该方程仍是 GKSL 形式的，所以 $\mathrm{Tr}\,\sigma(t)=1$（定理 8.6）。$\square$

请与命题 8.11 比较：它允许非中心的 $\sum_jL_j^*L_j$，但不允许扇区内部的跳跃；这里
损失项 $\gamma_1N$ 是中心元，因此无跳跃因子是数 $e^{-\gamma_1wt}$，并且允许
扇区内部的退相位。

当 $\gamma_2=0$ 时，$\sigma(t)=e^{-iHt}\rho(0)e^{iHt}$：后选择是完美的。

*计算：* 六个量子比特上 XY 混合算子的命题 10.10，使用 `dynamics.lindblad_evolve`——[教程，10.4 节](../../tutorial/zh/10_quantum_computing.md#s-10-4)。

> **物理意义**（物理解释，而非本课程的定理）。能量弛豫（$T_1$）是超导量子比特的主要误差之一，它是
> 跃迁到较低扇区的跳跃。测量汉明重量——一个中心投影，在扇区内部是无害的
> （3.3 节）——可以探测到它，而后选择以接受率为代价将其去除。
> 退相位与扇区投影可交换，因而探测不到。

<a id="s-10-5"></a>

## 10.5 可与之比较的工具

*库中：* 与其他软件（Qiskit、Cirq、PennyLane、QuTiP、stim）的比较——[教程，10.5 节](../../tutorial/zh/10_quantum_computing.md#s-10-5)。

<a id="x-10"></a>

## 习题

1. （定理 10.2 需要保单位性；比较 6.4 节。）对单个量子比特上的振幅阻尼，
   $E_0=|0\rangle\langle0|+\sqrt{1-g}\,|1\rangle\langle1|$，$E_1=\sqrt g\,|0\rangle\langle1|$，$0<g<1$，证明
   $\mathcal A'=\mathbb C1$，$\mathrm{Fix}(\Phi)=\mathbb C|0\rangle\langle0|$ 且 $\mathrm{Fix}(\Phi^*)=\mathbb C1$。
2. 证明 $\sum_w\binom nw^2=\binom{2n}n$。对推论 10.3 中取 $\theta=\pi/2$ 的信道（它违反了该推论的
   假设），$\mathcal A$ 和
   $\mathrm{Fix}(\Phi)$ 是什么？对 $n=4$ 检验 $\dim\mathrm{Fix}(\Phi)$。
3. （纠错的线性。）设对有限多个算子 $E_a$（不一定
   是某个信道的 Kraus 算子）有 $PE_a^*E_bP=\lambda_{ab}P$，并按定理 10.5 的证明由 $\Lambda$ 构造
   $\mathcal R$。证明 $\mathcal R$ 能纠正每个 Kraus 算子属于
   $\mathrm{span}\{E_a\}$ 的信道。由此推出比特翻转码能纠正相干误差 $e^{-i\phi X_1}$。

<details>
<summary>解答</summary>

1. $E_1$、$E_1^*$ 及其乘积给出所有矩阵单位，所以 $\mathcal A=M_2$，$\mathcal A'=\mathbb C1$。
   对 $T=(x_{kl})$，$\Phi(T)=\begin{pmatrix}x_{00}+gx_{11}&\sqrt{1-g}\,x_{01}\\\sqrt{1-g}\,x_{10}&(1-g)x_{11}\end{pmatrix}=T$
   当且仅当 $x_{11}=x_{01}=x_{10}=0$，而
   $\Phi^*(T)=\begin{pmatrix}x_{00}&\sqrt{1-g}\,x_{01}\\\sqrt{1-g}\,x_{10}&(1-g)x_{11}+gx_{00}\end{pmatrix}$
   等于 $T$ 当且仅当 $T=x_{00}1$。（$\Phi$ 不是保单位的，而 $\Phi^*$ 是保单位的但不保迹。）

2. 按前一半中元素的个数 $w$ 对 $\{1,\dots,2n\}$ 的 $n$ 元子集计数：
   $\binom{2n}n=\sum_w\binom nw\binom n{n-w}=\sum_w\binom nw^2$。对 $\theta=\pi/2$，
   $U=(-i)^n\sum_w(-1)^wP_w$ 只看 $w$ 的奇偶性，所以 $\mathcal A=\mathrm{span}\{P_{\mathrm{even}},P_{\mathrm{odd}}\}$，
   其中 $P_{\mathrm{even}}=\sum_{w\text{ 为偶数}}P_w$，并且由定理 10.2，$\mathrm{Fix}(\Phi)$ 由关于奇偶性
   块对角的算子组成：对 $n=4$，维数为 $2\cdot(2^{n-1})^2=128$。

3. 令 $F_k=\sum_au_{ak}E_a$，则 $\mathrm{span}\{E_a\}=\mathrm{span}\{F_k\}$，所以信道的每个 Kraus 算子都形如
   $G_m=\sum_kh_{mk}F_k$，且 $G_mP=\sum_{d_k>0}h_{mk}\sqrt{d_k}V_k$。由于 $V_l^*V_k=\delta_{lk}P$ 且 $(1-Q)V_k=0$，
   对 $\rho=P\rho P$ 有 $\mathcal R(G_m\rho G_m^*)=\sum_l|h_{ml}|^2d_l\,\rho$，所以 $\mathcal R(\Phi(\rho))=c\rho$；取迹得
   $c=1$。对 $E_a\in\{1,X_1,X_2,X_3\}$，$\Lambda=1$，$V_k=X_kP$，而 $\mathcal R$ 就是
   10.2 节的恢复；$e^{-i\phi X_1}=\cos\phi\,1-i\sin\phi\,X_1$ 属于该张成空间。

</details>

*计算：* 习题 1–3 的数值验证——[教程，习题](../../tutorial/zh/10_quantum_computing.md#x-10)。

**下一章：** [第 11 章](11_optimisation.md)：代数的酉群上的优化，以及迹不等式。
