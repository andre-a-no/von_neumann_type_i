[English](../05_random_matrices.md) · [Español](../es/05_random_matrices.md) · [Français](../fr/05_random_matrices.md) · [Deutsch](../de/05_random_matrices.md) · **中文** · [日本語](../ja/05_random_matrices.md) · [Русский](../ru/05_random_matrices.md)

# 5. 代数中的随机矩阵

随机矩阵理论关心"典型"算子是什么样子。在代数 $\bigoplus_cM_{k_c}$ 中，答案必须尊重扇区：
自然的随机酉算子是该代数的酉元，而不是大矩阵代数 $M_D$ 的酉矩阵；自然的随机哈密顿量
具有相互独立的块。本章构造代数酉群上的 Haar 测度，精确计算它的一阶矩与二阶矩，然后说明
忽略扇区时谱统计会出现什么问题。

**你将学到**为什么高斯矩阵的 QR 分解给出 Haar 酉矩阵（以及为什么必须固定相位），在代数的
酉群上求平均就是到中心上的条件期望，二阶 Weingarten 公式，Wigner 猜测，间距比及其对 $3\times3$
高斯矩阵的精确值，以及为什么在混合了扇区的谱上计算的能级统计看起来像独立能级的统计。

*计算：* 本章的数值检验——[教程 5](../../tutorial/zh/05_random_matrices.md)。

在本章中，$\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$，$D=\sum_ck_c$，$P_c$ 为
极小中心投影，$E$ 是命题 2.8 中到中心上的条件期望。$U(k)$ 是 $k\times k$ 酉矩阵构成的群，
$U(\mathcal M)=\{U\in\mathcal M:U^*U=1\}$ 是 $\mathcal M$ 的酉群。随机变量定义在一个
概率空间上，期望记为 $\mathbb E$；随机矩阵的期望按元素逐个取。*标准复高斯*变量的实部与虚部
相互独立且 $\sim N(0,\tfrac12)$，从而 $\mathbb E|z|^2=1$；*Ginibre 矩阵*的元素是相互独立的标准复高斯变量。

<a id="s-5-1"></a>

## 5.1 Haar 测度

设 $G$ 是一个紧矩阵群（某个 $U(n)$ 的闭子群）。$G$ 上的 Borel 概率测度
$\mu$ 称为**左不变的**，如果对所有 $g\in G$ 和 Borel 集 $B$ 有 $\mu(gB)=\mu(B)$，
等价地，对每个有界 Borel 函数 $f$ 有 $\int f(gh)\,d\mu(h)=\int f\,d\mu$；
**右不变**用 $Bg$ 和 $f(hg)$ 类似地定义。

<a id="r-5-1"></a>

**引理 5.1（唯一性）。** 若 $\mu$ 是 $G$ 上左不变的 Borel 概率测度，$\nu$ 是右不变的 Borel
概率测度，则 $\mu=\nu$。因此 $G$ 上至多有一个左不变的 Borel
概率测度；若它存在，则它也是右不变的，并且在 $g\mapsto g^{-1}$ 下不变。

*证明。* 设 $f:G\to\mathbb R$ 连续（从而有界）。由 Fubini 定理，
$\int\!\big(\int f(gh)\,d\mu(h)\big)d\nu(g)=\int\!\big(\int f(gh)\,d\nu(g)\big)d\mu(h)$。左边的内层
积分由左不变性等于 $\int f\,d\mu$，右边的内层积分由右不变性等于 $\int f\,d\nu$；
所以 $\int f\,d\mu=\int f\,d\nu$。度量空间上两个 Borel 概率测度若对连续函数的积分都相同，
则二者相等（推论 A.2：连续函数
$\max(0,1-n\mathrm{dist}(\cdot,K))$ 递减地收敛到闭集 $K$ 的示性函数，而由 Dynkin 引理（引理 A.1），
闭集决定测度）。若 $\mu_1,\mu_2$ 都是左不变的，则
$\nu(B)=\mu_2(B^{-1})$ 是右不变的，因为 $(Bg)^{-1}=g^{-1}B^{-1}$；因此 $\mu_1=\nu=\mu_2$，
所以 $\mu_1$ 是右不变且在取逆下不变的。$\square$

**定义。** $G$ 的 **Haar 测度**是它的左不变 Borel 概率测度。以此为分布的
随机元素称为**服从 Haar 分布**。

引理 5.1 并没有说明存在性。本课程只需要 $U(k)$、$U(1)$、$U(\mathcal M)$、$SU(k)$、$O(k)$
以及引理 C.2 中的对角群 $T$ 的 Haar 测度，而它们每一个都被显式地构造出来：$U(k)$ 在定理 5.3 中，
$U(\mathcal M)$ 在推论 5.4 中，$SU(k)$ 紧接其后，$O(k)$ 在定理 B.17 中，$T$ 在引理 C.2 中（见 B.7 节）。

<a id="r-5-2"></a>

**引理 5.2（QR 分解）。** 每个可逆的 $Z\in M_k(\mathbb C)$ 都可以写成 $Z=QR$，
其中 $Q$ 是酉矩阵，$R$ 是对角元为正的上三角矩阵。$Q$ 与 $R$ 是唯一的，并且
连续地依赖于 $Z$。

*证明。* 对 $Z$ 的线性无关的列 $z_1,\dots,z_k$ 施行 Gram–Schmidt 正交化：
$q_j=w_j/\|w_j\|$，其中 $w_j=z_j-\sum_{i<j}\langle q_i,z_j\rangle q_i$；这里 $w_j\neq0$，因为
$z_j\notin\mathrm{span}\{z_1,\dots,z_{j-1}\}=\mathrm{span}\{q_1,\dots,q_{j-1}\}$。于是
$z_j=\sum_{i<j}\langle q_i,z_j\rangle q_i+\|w_j\|q_j$，即 $Z=QR$，其中 $R_{jj}=\|w_j\|>0$，并且
每一步在可逆矩阵上都是连续的。若 $Q_1R_1=Q_2R_2$，则
$T=Q_2^*Q_1=R_2R_1^{-1}$ 是对角元为正的酉上三角矩阵。它的逆
$T^*$ 是下三角的，同时也是上三角的（上三角矩阵的逆是上三角的），
所以 $T$ 是对角的；对角元为正的对角酉矩阵就是 $1$。$\square$

<a id="r-5-3"></a>

**定理 5.3（由高斯矩阵得到 Haar 酉矩阵）。** 设 $Z$ 是 $k\times k$
Ginibre 矩阵。则 $Z$ 以概率一可逆，并且引理 5.2 中的因子 $Q$ 在 $U(k)$ 上
服从 Haar 分布。此外，若 $Z=Q'R'$ 是任意一个 QR 分解，其中 $Q'$ 为酉矩阵、$R'$
为上三角矩阵（例如 `torch.linalg.qr` 返回的分解），则
$Q=Q'\Lambda$（其中 $\Lambda=\mathrm{diag}(R'_{ii}/|R'_{ii}|)$）服从 Haar 分布。

*证明。* *可逆性。* $|\det Z|^2$ 是 $Z$ 的 $2k^2$ 个实坐标的实多项式，不恒为零
（它在 $Z=1$ 处等于 $1$），而 $Z$ 有密度 $\pi^{-k^2}e^{-\|Z\|_2^2}$
（命题 A.16(b)）；所以由推论 A.7(b)，$P(\det Z=0)=0$。

*高斯分布的不变性。* 对 $V\in U(k)$，映射 $Z\mapsto VZ$ 是实线性的且保持 $\|Z\|_2$
（1.4 节），所以它是 $M_k(\mathbb C)\cong\mathbb R^{2k^2}$ 的正交映射，它保持 Lebesgue 测度
（引理 A.3(b)）并保持密度：$VZ$ 与 $Z$ 同分布（命题 A.16(b)）。

*等变性。* 若 $Z=QR$，则 $VZ=(VQ)R$ 是 $VZ$ 的引理 5.2 意义下的分解，所以
$Q(VZ)=VQ(Z)$。由于 $Q$ 在可逆矩阵上连续（引理 5.2），$Q(Z)$ 是 $U(k)$ 中的随机元素，
并且对每个 Borel 集 $B\subseteq U(k)$，
$P(VQ(Z)\in B)=P(Q(VZ)\in B)=P(Q(Z)\in B)$：$Q(Z)$ 的分布是左不变的，由
引理 5.1 它就是 Haar 测度（这也证明了 $U(k)$ 的 Haar 测度存在）。
最后 $Q'R'=(Q'\Lambda)(\Lambda^*R')$，而 $\Lambda^*R'$ 是对角元为
$|R'_{ii}|>0$ 的上三角矩阵，所以由唯一性得 $Q'\Lambda=Q(Z)$。$\square$

对实高斯矩阵作同样的证明，就得到正交群 $O(k)$ 的 Haar 测度，
修正因子为 $\mathrm{diag}(\mathrm{sgn}\,R'_{ii})$；这在定理 B.17 中详细写出。
若不作修正 $\Lambda$，$Q'$ 的分布依赖于 QR 例程的符号约定，未必是 Haar 测度（带 $\Lambda$ 的做法出自 F. Mezzadri, *How to
generate random matrices from the classical compact groups*, Notices AMS 54 (2007), 592–604）。

*计算：* 本库的采样器 `TypeIAlgebra.random_unitary` 就是定理 5.3 中经过修正的 QR，而 `torch.linalg.qr` 给出的未经修正的 $Q'$ 不服从 Haar 分布——[教程，5.1 节](../../tutorial/zh/05_random_matrices.md#s-5-1)。

<a id="r-5-4"></a>

**推论 5.4（代数的 Haar 测度）。** $U(\mathcal M)=\{\bigoplus_cU_c:U_c\in U(k_c)\}\cong\prod_cU(k_c)$，
并且 $U=\bigoplus_cU_c$ 在 $U(\mathcal M)$ 上服从 Haar 分布当且仅当各块 $U_c$ 相互独立，
且每个 $U_c$ 在 $U(k_c)$ 上服从 Haar 分布。

*证明。* $U^*U=\bigoplus_cU_c^*U_c$，所以 $U$ 是酉的当且仅当每个块都是酉的。设各 $U_c$
相互独立且服从 Haar 分布。对 $V=\bigoplus_cV_c\in U(\mathcal M)$，$VU$ 的各块 $V_cU_c$
相互独立（独立变量的函数），且 $V_cU_c$ 与 $U_c$ 同分布，所以
$VU$ 与 $U$ 同分布。于是 $U$ 的分布是左不变的，由引理 5.1 它就是紧群
$U(\mathcal M)$ 的 Haar 测度。反之，由唯一性，每个服从 Haar 分布的 $U$ 都具有
这一分布。$\square$

*库中：* $U(\mathcal M)$ 上各块相互独立的 Haar 酉元由 `TypeIAlgebra.random_unitary_operator` 抽取——[教程，5.1 节](../../tutorial/zh/05_random_matrices.md#s-5-1)。

群 $U(\mathcal M)$ 比 $U(D)$ 小得多：当 $C\ge2$ 时，$U(D)$ 的 Haar 酉矩阵 $U$ 以概率零
落在 $\mathcal M$ 中。事实上，由定理 5.3 和引理 5.1，$U$ 与 $D\times D$ Ginibre 矩阵 $Z$ 的因子 $Q$
同分布，而由引理 5.2 中的 Gram–Schmidt 构造，$Q$ 的第一列是 $z/\|z\|$，其中 $z$ 是 $Z$ 的第一列，
即 $\mathbb C^D$ 中具有密度的标准复高斯向量（命题 A.16(a)）。对每个 $i$，$\{z_i=0\}$ 是
$\mathbb C^D$ 的真子空间，所以 $P(z_i=0)=0$（推论 A.7(c)），从而第一列的所有元素几乎必然
非零。而 $\mathcal M$ 的元素在第一列中、第一个块的行之外的元素都为零，当 $C\ge2$ 时这样的行是存在的。

$SU(k)$ 的 Haar 测度可以用同样的思路得到：若 $U$ 在 $U(k)$ 上服从 Haar 分布且 $\theta=\arg\det U\in(-\pi,\pi]$（$U$ 的 Borel 函数），则 $V=e^{-i\theta/k}U\in SU(k)$，并且对 $W\in SU(k)$，
$\det(WU)=\det U$ 给出 $V(WU)=WV(U)$；由于 $WU$ 与 $U$ 同分布，$V$ 的分布在 $SU(k)$ 上是
左不变的，因此由引理 5.1 它就是 $SU(k)$ 的 Haar 测度。

<a id="r-5-5"></a>

**推论 5.5（元素的一阶矩）。** 若 $U$ 在 $U(k)$ 上服从 Haar 分布，则
$\mathbb EU_{ij}=0$ 且 $\mathbb E|U_{ij}|^2=1/k$。

*证明。* $-1\in U(k)$，所以 $-U$ 与 $U$ 同分布，$\mathbb EU=-\mathbb EU$。对置换
矩阵 $P$，$PU$ 与 $U$ 同分布，所以 $\mathbb E|U_{ij}|^2$ 与 $i$ 无关；由于
$\sum_i|U_{ij}|^2=1$，它等于 $1/k$。$\square$

*计算：* 用蒙特卡罗方法检验推论 5.4 和 5.5，以及 $SU(k)$ 的 Haar 测度（`measure='haar_su'`）——[教程，5.1 节](../../tutorial/zh/05_random_matrices.md#s-5-1)。

<a id="s-5-2"></a>

## 5.2 矩与 Weingarten 演算

<a id="r-5-6"></a>

**引理 5.6.** 若 $T\in M_k(\mathbb C)$ 与所有对角酉矩阵和所有置换
矩阵可交换，则 $T\in\mathbb C1$。特别地，$U(k)$ 的交换子是 $\mathbb C1$。

*证明。* 对 $\Delta=\mathrm{diag}(e^{i\theta_1},\dots,e^{i\theta_k})$，其中各
$e^{i\theta_j}$ 互不相同，$\Delta T=T\Delta$ 给出 $(e^{i\theta_i}-e^{i\theta_j})T_{ij}=0$，所以 $T$ 是
对角的。对置换矩阵 $P_\sigma e_j=e_{\sigma(j)}$，$P_\sigma TP_\sigma^{-1}=T$ 给出
$T_{\sigma(j)\sigma(j)}=T_{jj}$，所以所有对角元都相等。$\square$

<a id="r-5-7"></a>

**命题 5.7（Haar twirl 就是到中心上的条件期望）。** 设 $U$ 在
$U(\mathcal M)$ 上服从 Haar 分布。则对每个 $X\in\mathcal M$，
$$
\mathbb E\,[UXU^*]=E(X)=\bigoplus_c\frac{\mathrm{Tr}\,X_c}{k_c}\,1_{k_c}.
$$

*证明。* $(UXU^*)_c=U_cX_cU_c^*$，且 $U_c$ 在 $U(k_c)$ 上服从 Haar 分布（推论 5.4）。令
$T_c=\mathbb E[U_cX_cU_c^*]$（各元素有界）。对 $V\in U(k_c)$，$VU_c$ 与
$U_c$ 同分布，所以 $VT_cV^*=\mathbb E[(VU_c)X_c(VU_c)^*]=T_c$。由引理 5.6，$T_c=\lambda1$，并且
$k_c\lambda=\mathrm{Tr}\,T_c=\mathbb E\mathrm{Tr}(U_cX_cU_c^*)=\mathrm{Tr}\,X_c$。$\square$

这里的归一化与命题 2.8 以及 `channels.center_expectation` 相同。若
$X\in\mathcal M$ 与 $U(\mathcal M)$ 可交换，则 $UXU^*=X$，所以 $X=E(X)$ 是中心元；Haar
平均是第 2 章习题 3 中群平均的连续类比。若 $U$ 改为在 $U(D)$ 上服从 Haar 分布，
同样的证明（取 $C=1$）给出 $\mathbb E[UXU^*]=\frac{\mathrm{Tr}\,X}{D}1$：
扇区权重丢失了，而且 $UXU^*$ 甚至不一定落在 $\mathcal M$ 中（对 $C\ge2$，取 $X=P_1$，$U$ 为
$e_1$ 与 $e_{k_1+1}$ 张成的平面内转角为 $\pi/4$ 的旋转，则 $UP_1U^*$ 在位置 $(1,k_1+1)$ 处的元素为 $\frac12$）。

*计算：* 用蒙特卡罗方法检验命题 5.7，并与 `channels.center_expectation` 对照——[教程，5.2 节](../../tutorial/zh/05_random_matrices.md#s-5-2)。

> **物理意义**（物理解释，而非本课程的定理）。代数的 Haar 随机酉元是一个遵守超选择定则的随机演化。
> 对它求平均会把每个扇区内部的一切打乱，而恰好保留经典部分：每个荷的
> 权重。这就是量子信息中的"twirl"，它正是第 2 章中到中心上的条件期望。

对于二阶矩，我们需要 $\{V\otimes V\}$ 在 $\mathbb C^k\otimes\mathbb C^k$ 上的交换子。
设 $F$ 为**交换算子**（swap），$F(x\otimes y)=y\otimes x$；它是自伴且酉的，并且
$\mathrm{Tr}\,1=k^2$，$\mathrm{Tr}\,F=\sum_{i,j}\langle e_i\otimes e_j,e_j\otimes e_i\rangle=k$，$F^2=1$。

<a id="r-5-8"></a>

**引理 5.8.** 设 $k\ge2$。$\mathbb C^k\otimes\mathbb C^k$ 上的矩阵 $X$ 与所有 $V\in U(k)$ 的 $V\otimes V$
可交换，当且仅当 $X=a1+bF$，其中 $a,b\in\mathbb C$；并且 $1,F$ 线性无关。

*证明。* $F(V\otimes V)(x\otimes y)=Vy\otimes Vx=(V\otimes V)F(x\otimes y)$，所以 $1$ 和 $F$ 与
所有 $V\otimes V$ 可交换。反之，设 $X$ 与它们可交换，并记
$X_{ij,pq}=\langle e_i\otimes e_j,X(e_p\otimes e_q)\rangle$。
*对角酉矩阵：* $\Delta\otimes\Delta$ 把 $e_i\otimes e_j$ 乘以 $e^{i(\theta_i+\theta_j)}$，
所以 $(e^{i(\theta_i+\theta_j)}-e^{i(\theta_p+\theta_q)})X_{ij,pq}=0$。若多重集 $\{i,j\}$ 与
$\{p,q\}$ 不同，则某个指标 $m$ 在前者中比在后者中多出现 $n\in\{\pm1,\pm2\}$ 次；取
$\theta_m=\pi/2$、其余 $\theta=0$，括号等于 $e^{i(\theta_p+\theta_q)}(e^{in\pi/2}-1)\neq0$，
所以 $X_{ij,pq}=0$。因此 $X(e_i\otimes e_i)=a_ie_i\otimes e_i$，并且对 $i\neq j$，
$X(e_i\otimes e_j)=b_{ij}\,e_i\otimes e_j+c_{ij}\,e_j\otimes e_i$。
*置换：* 把 $X(P_\sigma\otimes P_\sigma)=(P_\sigma\otimes P_\sigma)X$ 作用于 $e_i\otimes e_j$，得到
$a_{\sigma(i)}=a_i$，$b_{\sigma(i)\sigma(j)}=b_{ij}$，$c_{\sigma(i)\sigma(j)}=c_{ij}$。置换群在由不同指标
构成的有序对上传递地作用，所以 $a_i=a'$，$b_{ij}=b'$，$c_{ij}=c'$；于是在 $i\ne j$ 的
$e_i\otimes e_j$ 张成的子空间上 $X=b'1+c'F$。
*一个旋转：* 设 $R$ 是满足 $Re_1=(e_1+e_2)/\sqrt2$，$Re_2=(-e_1+e_2)/\sqrt2$，以及当 $m\ge3$ 时 $Re_m=e_m$ 的酉矩阵。
则 $(R\otimes R)(e_1\otimes e_1)=\frac12(e_1\otimes e_1+e_1\otimes e_2+e_2\otimes e_1+e_2\otimes e_2)$，而
$X(R\otimes R)(e_1\otimes e_1)=\frac12\big(a'e_1\otimes e_1+(b'+c')(e_1\otimes e_2+e_2\otimes e_1)+a'e_2\otimes e_2\big)$
必须等于 $(R\otimes R)X(e_1\otimes e_1)=a'(R\otimes R)(e_1\otimes e_1)$。比较
$e_1\otimes e_2$ 的系数得 $a'=b'+c'$，这正是 $b'1+c'F$ 在 $e_i\otimes e_i$ 上的值。所以 $X=b'1+c'F$。
最后，$1$ 和 $F$ 把 $e_1\otimes e_2$ 分别映为线性无关的向量 $e_1\otimes e_2$ 和 $e_2\otimes e_1$。$\square$

*计算：* 用线性代数检验引理 5.8，计算三个抽样得到的 $V\otimes V$ 的交换子——[教程，5.2 节](../../tutorial/zh/05_random_matrices.md#s-5-2)。

<a id="r-5-9"></a>

**推论 5.9（二阶 Schur–Weyl 对偶）。** 设 $k\ge2$，$P_\pm=\frac12(1\pm F)$。
它们是到对称子空间与反对称子空间
$\mathrm{Sym}=\{x:Fx=x\}$ 和 $\mathrm{Alt}=\{x:Fx=-x\}$ 上的投影，维数分别为 $k(k+1)/2$ 和 $k(k-1)/2$。
$U(k)$ 的表示 $V\mapsto V\otimes V$ 在 $\mathrm{Sym}$ 上和在 $\mathrm{Alt}$ 上都是不可约的，
并且这两个表示不等价。

*证明。* $F=F^*$ 和 $F^2=1$ 使 $P_\pm$ 成为投影，满足 $P_++P_-=1$，值域分别为 $\mathrm{Sym}$、
$\mathrm{Alt}$；它们的维数是 $\mathrm{Tr}\,P_\pm=(k^2\pm k)/2$。
$\{V\otimes V\}$ 的线性张成 $\mathcal A$ 是一个 $*$-代数（$1\otimes1$、乘积 $VW\otimes VW$、伴随 $V^*\otimes V^*$），
由引理 5.8，其交换子为 $\mathcal A'=\mathrm{span}\{1,F\}=\mathrm{span}\{P_+,P_-\}$。设
$K\subseteq\mathrm{Sym}$ 在所有 $V\otimes V$ 下不变。由引理 2.1，它的投影属于 $\mathcal A'$，所以它
等于 $xP_++yP_-$，其中 $x,y\in\{0,1\}$（因为它是投影）；它在 $\mathrm{Alt}$ 上为零，所以 $y=0$，
从而 $K=0$ 或 $K=\mathrm{Sym}$。同样的论证适用于 $\mathrm{Alt}$。若
$T:\mathrm{Sym}\to\mathrm{Alt}$ 是这两个表示之间的缠结算子，则 $TP_+$（看作
$\mathbb C^k\otimes\mathbb C^k$ 上的算子）与 $\mathcal A$ 可交换，所以 $TP_+=xP_++yP_-$；它
在 $\mathrm{Alt}$ 上为零，所以 $y=0$，而它的值域含于 $\mathrm{Alt}\cap\mathrm{Sym}=0$，所以 $T=0$。$\square$

<a id="r-5-10"></a>

**定理 5.10（二阶 Weingarten 公式）。** 设 $k\ge2$，$U$ 在 $U(k)$ 上服从 Haar 分布，
$X$ 是 $\mathbb C^k\otimes\mathbb C^k$ 上的矩阵。则
$$
\mathbb E\,[(U\otimes U)X(U\otimes U)^*]=\alpha\,1+\beta F,\qquad
\alpha=\frac{\mathrm{Tr}\,X-\mathrm{Tr}(FX)/k}{k^2-1},\quad
\beta=\frac{\mathrm{Tr}(FX)-\mathrm{Tr}\,X/k}{k^2-1}.
$$
等价地，在群 $S_2=\{e,\tau\}$ 上令 $\mathrm{Wg}(e)=\frac1{k^2-1}$，$\mathrm{Wg}(\tau)=-\frac1{k(k^2-1)}$，则
$$
\mathbb E\big[U_{i_1j_1}U_{i_2j_2}\overline{U_{i'_1j'_1}U_{i'_2j'_2}}\big]
=\sum_{\sigma,\pi\in S_2}\delta_{i_1i'_{\sigma(1)}}\delta_{i_2i'_{\sigma(2)}}\,
\delta_{j_1j'_{\pi(1)}}\delta_{j_2j'_{\pi(2)}}\,\mathrm{Wg}(\sigma\pi^{-1}).
$$

*证明。* 令 $T=\mathbb E[(U\otimes U)X(U\otimes U)^*]$。与命题 5.7 一样，$T$ 与
所有 $V\otimes V$ 可交换，所以由引理 5.8，$T=\alpha1+\beta F$。取迹得
$\mathrm{Tr}\,T=\mathrm{Tr}\,X$；又因为 $F$ 与 $U\otimes U$ 可交换，
$\mathrm{Tr}(FT)=\mathbb E\mathrm{Tr}((U\otimes U)^*F(U\otimes U)X)=\mathrm{Tr}(FX)$。
这给出 $k^2\alpha+k\beta=\mathrm{Tr}\,X$ 和 $k\alpha+k^2\beta=\mathrm{Tr}(FX)$，这是一个
行列式为 $k^2(k^2-1)\neq0$ 的线性方程组，其解如上所述。对于
矩阵元，取 $X=|e_{j_1}\otimes e_{j_2}\rangle\langle e_{j'_1}\otimes e_{j'_2}|$：则
$\langle e_{i_1}\otimes e_{i_2},T\,e_{i'_1}\otimes e_{i'_2}\rangle$ 就是左边，$\mathrm{Tr}\,X=\delta_{j_1j'_1}\delta_{j_2j'_2}$
（$\pi=e$ 的项），$\mathrm{Tr}(FX)=\delta_{j_1j'_2}\delta_{j_2j'_1}$（$\pi=\tau$ 的项），而
$1$、$F$ 分别贡献 $\delta_{i_1i'_1}\delta_{i_2i'_2}$（$\sigma=e$ 的项）和
$\delta_{i_1i'_2}\delta_{i_2i'_1}$（$\sigma=\tau$）。合并各项得
$\alpha=\sum_\pi\mathrm{Wg}(e\pi^{-1})[\pi]$，$\beta=\sum_\pi\mathrm{Wg}(\tau\pi^{-1})[\pi]$，
其中 $[\pi]$ 表示 $\pi$ 对应的 $j$-δ 乘积。$\square$

矩阵 $\big(\mathrm{Wg}(\sigma^{-1}\pi)\big)_{\sigma,\pi}$ 是 $1$ 与 $F$ 关于
$\langle A,B\rangle=\mathrm{Tr}(A^*B)$ 的 Gram 矩阵
$\begin{pmatrix}k^2&k\\k&k^2\end{pmatrix}$ 的逆；这正是可以推广的形式。

*计算：* 对 $k=3$ 用蒙特卡罗方法检验定理 5.10——[教程，5.2 节](../../tutorial/zh/05_random_matrices.md#s-5-2)。

<a id="r-5-11"></a>

**推论 5.11.** 对 $k\ge2$：
$$
\mathbb E|U_{11}|^4=\frac2{k(k+1)},\qquad
\mathbb E|U_{11}|^2|U_{12}|^2=\frac1{k(k+1)},\qquad
\mathbb E\big[U_{11}U_{22}\overline{U_{12}U_{21}}\big]=-\frac1{(k-1)k(k+1)}.
$$

*证明。* 在定理 5.10 中，(i) 所有指标都等于 $1$：每个 δ 乘积都是 $1$，所以该
矩等于 $2(\mathrm{Wg}(e)+\mathrm{Wg}(\tau))=2\frac{k-1}{k(k^2-1)}$。(ii) $i=(1,1)$，$j=(1,2)$，
$i'=(1,1)$，$j'=(1,2)$：只有 $\pi=e$ 的项留下，得到 $\mathrm{Wg}(e)+\mathrm{Wg}(\tau)$。
(iii) $i=(1,2)$，$j=(1,2)$，$i'=(1,2)$，$j'=(2,1)$：只有 $\sigma=e$，$\pi=\tau$ 的项留下，
得到 $\mathrm{Wg}(\tau)$。$\square$

*计算：* 对 $k=3$ 用蒙特卡罗方法检验推论 5.11 的三个矩——[教程，5.2 节](../../tutorial/zh/05_random_matrices.md#s-5-2)。

最后一个矩是负的：Haar 酉矩阵的元素不是相互独立的，其
相关性的相对量级为 $1/k$。对于代数，不同块的元素是
相互独立的（推论 5.4），所以跨扇区的混合矩分解为单个块的矩的乘积。

对任意阶的矩，只要知道
$\{V^{\otimes n}\}$ 的交换子，同样的推理就适用；它由 $n$ 个张量因子的置换张成
（Schur–Weyl 对偶，定理 B.6(c)）。附录 B 完成了这一推导。

<a id="r-5-12"></a>

**定理 5.12（任意阶的 Weingarten 公式）。** 设 $n\ge1$，$U$ 在 $U(k)$ 上服从 Haar 分布。则
$$
\mathbb E\big[U_{i_1j_1}\cdots U_{i_nj_n}\overline{U_{i'_1j'_1}\cdots U_{i'_nj'_n}}\big]
=\sum_{\sigma,\pi\in S_n}\prod_{m=1}^n\delta_{i_mi'_{\sigma(m)}}\delta_{j_mj'_{\pi(m)}}\;\mathrm{Wg}(\sigma\pi^{-1},k),
$$
其中 $\mathrm{Wg}(\cdot,k)$ 是 $S_n$ 上只依赖于循环型的函数。当 $n\le k$ 时，
矩阵 $\big(\mathrm{Wg}(\sigma^{-1}\pi,k)\big)_{\sigma,\pi}$ 是 Gram
矩阵 $\big(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)}\big)_{\sigma,\pi}$ 的逆；当 $n>k$ 时该矩阵是奇异的，
此时取其伪逆。$U$ 与 $\overline U$ 个数不同的矩为零。

*证明。* 这就是 $d=k$ 时的定理 B.13：(a) 部分给出 $W_{\sigma\pi}=\mathrm{Wg}(\sigma^{-1}\pi)$，其中 $\mathrm{Wg}$
是循环型的函数；(b) 部分给出 Gram 矩阵可逆当且仅当 $n\le k$；(d) 部分给出矩的
公式；(e) 部分给出不平衡矩为零。$\square$

（该公式出自 D. Weingarten, J. Math. Phys. 19 (1978), 999–1001，以及 B. Collins, Int. Math. Res.
Not. 2003, 953–982；另见 B. Collins, P. Śniady, Comm. Math. Phys. 264 (2006), 773–795。）

*计算：* $k=3$ 时的 $3$ 阶 Weingarten 矩阵与 $\mathbb E|U_{11}|^6=\frac1{10}$（习题 2）——[教程，习题](../../tutorial/zh/05_random_matrices.md#x-5)。

<a id="s-5-3"></a>

## 5.3 能级统计

设 $H$ 是一个随机自伴 $n\times n$ 矩阵，其特征值为
$\lambda_1\le\dots\le\lambda_n$。它的**能级间距**为 $s_i=\lambda_{i+1}-\lambda_i$。有两个标准的
系综：**GOE**，$H=(G+G^T)/2$，其中 $G$ 是元素相互独立且服从 $N(0,1)$ 的实矩阵；以及
**GUE**，$H=(G+G^*)/\sqrt2$，其中 $G$ 是 Ginibre 矩阵（1.7 节）。在这两种情形中，对角元都服从
$N(0,1)$，对角线以上的元素在 GOE 中服从 $N(0,\tfrac12)$，在 GUE 中其实部与虚部相互独立且都服从 $N(0,\tfrac12)$，
元素 $h_{ij}$（$i\le j$）相互独立，而
$H$ 的密度正比于 $e^{-\mathrm{Tr}\,H^2/2}$（命题 A.17；正比关系来自
$\mathrm{Tr}\,H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2$）。

<a id="r-5-13"></a>

**命题 5.13（Wigner 猜测）。** 设 $n=2$，$s=\lambda_2-\lambda_1$，$S=s/\mathbb Es$。
对 GOE（$\beta=1$）与 GUE（$\beta=2$），$S$ 的密度分别为
$$
p_1(S)=\frac\pi2\,S\,e^{-\pi S^2/4},\qquad p_2(S)=\frac{32}{\pi^2}\,S^2\,e^{-4S^2/\pi},
$$
分布函数分别为 $P(S\le t)=1-e^{-\pi t^2/4}$ 和
$P(S\le t)=\mathrm{erf}(2t/\sqrt\pi)-\frac{4t}\pi e^{-4t^2/\pi}$。特别地，当 $\varepsilon\to0$ 时，
$P(S\le\varepsilon)$ 的量级为 $\varepsilon^{\beta+1}$。

*证明。* 对 $H=\begin{pmatrix}a&b\\\bar b&c\end{pmatrix}$，特征值为
$\frac{a+c}2\pm\frac12\sqrt{(a-c)^2+4|b|^2}$，所以 $s=\sqrt{(a-c)^2+4|b|^2}$。由命题 A.17，$a$、$c$ 与
$b$ 的实坐标相互独立，$a,c\sim N(0,1)$，且 $b$ 的每个实坐标服从 $N(0,\frac12)$；所以
$a-c\sim N(0,2)$（命题 A.14(c)），而 $2\mathrm{Re}\,b$、$2\mathrm{Im}\,b$ 服从 $N(0,2)$，它们全都相互独立。
*GOE：* $x=a-c$ 与 $y=2b$ 是相互独立的 $N(0,2)$ 变量，在极坐标下（推论 A.19，取
$n=2$，$s_2=2\pi$）
$P(s\le t)=\int_{x^2+y^2\le t^2}\frac{e^{-(x^2+y^2)/4}}{4\pi}\,dx\,dy=\int_0^t\frac r2e^{-r^2/4}dr=1-e^{-t^2/4}$。
于是 $\mathbb Es=\int_0^\infty P(s>t)\,dt=\int_0^\infty e^{-t^2/4}dt=\sqrt\pi$，且
$P(S\le t)=P(s\le\sqrt\pi t)=1-e^{-\pi t^2/4}$；再求导即可。*GUE：* $x=a-c$，$y=2\mathrm{Re}\,b$，
$y'=2\mathrm{Im}\,b$ 是相互独立的 $N(0,2)$ 变量，由定理 A.20（$n=3$，$s_3=4\pi$），
$s=\sqrt{x^2+y^2+y'^2}$ 的密度为 $4\pi r^2(4\pi)^{-3/2}e^{-r^2/4}=\frac{r^2e^{-r^2/4}}{2\sqrt\pi}$。其均值为
$\frac1{2\sqrt\pi}\int_0^\infty r^3e^{-r^2/4}dr=\frac{8}{2\sqrt\pi}=\frac4{\sqrt\pi}$。代入
$r=4S/\sqrt\pi$ 得到 $p_2$，再对 $p_2$ 分部积分
（$\int x^2e^{-ax^2}dx=-\frac{x}{2a}e^{-ax^2}+\frac1{2a}\int e^{-ax^2}dx$，其中 $a=4/\pi$）即得
分布函数。小 $\varepsilon$ 时的行为由 $p_\beta(S)\sim cS^\beta$ 得出。$\square$

*计算：* 对 $2\times2$ 的 GOE 与 GUE 用蒙特卡罗方法检验命题 5.13——[教程，5.3 节](../../tutorial/zh/05_random_matrices.md#s-5-3)。

指数是在数维数：实对称 $2\times2$ 矩阵简并当且仅当两个实
条件成立（$a=c$，$b=0$），复 Hermite 矩阵则当且仅当三个条件成立。对每个 $n$，GOE（$\beta=1$）
与 GUE（$\beta=2$）的有序特征值 $\lambda_1<\dots<\lambda_n$ 具有联合密度
$$
Z_{n,\beta}^{-1}\prod_{i<j}(\lambda_j-\lambda_i)^\beta\,e^{-\sum_i\lambda_i^2/2}
$$
（定理 C.7，在附录 C 中借助 Weyl 积分公式证明），并且最小间距 $\delta$
对 $0<\varepsilon\le1$ 满足 $c\,\varepsilon^{\beta+1}\le P(\delta<\varepsilon)\le C\varepsilon^{\beta+1}$，其中常数
$0<c\le C$ 依赖于 $n$ 和 $\beta$（命题 C.9）。本课程不讨论 $n\to\infty$ 的极限
（见 M. L. Mehta, *Random Matrices*, 3rd ed., Elsevier, 2004）。

> **物理意义**（物理解释，而非本课程的定理）。Wigner 提出，
> 复杂量子系统（重原子核）的能级在统计上类似于具有恰当对称性的随机哈密顿量的特征值：
> 若存在时间反演对称性则为实对称矩阵（$\beta=1$），否则为复 Hermite 矩阵
> （$\beta=2$）。于是能级相互*排斥*：同一对称性的两个能级几乎从不靠近
> （von Neumann 与 Wigner，1929）。能级排斥是量子混沌的典型特征，而独立
> 能级（无排斥）则是可积系统或局域化系统的特征；这就是物理学中的 Bohigas–Giannoni–Schmit
> 猜想与 Berry–Tabor 猜想，本课程既不证明也不使用它们。

间距依赖于局部能级密度，而后者沿谱变化；要将间距与 $p_\beta$ 比较，需要进行
**展开**（unfolding），即借助计数函数的光滑近似 $\bar N$ 作变量替换 $\lambda\mapsto\bar N(\lambda)$，
使平均间距为 $1$。结果依赖于 $\bar N$ 的选取。相邻间距之比则避免了这一选择。

**定义。** **间距比**为 $r_i=\min(s_i,s_{i+1})/\max(s_i,s_{i+1})\in[0,1]$
（V. Oganesyan, D. A. Huse, Phys. Rev. B 75 (2007), 155111）。它在
$\lambda\mapsto a\lambda+b$（$a>0$）下不变，所以它不依赖于局部平均间距；这就是它
无需展开即可使用的原因。

*库中：* 每一行能级的平均间距比由 `SpinChain.level_spacing_ratio` 给出——[教程，5.3 节](../../tutorial/zh/05_random_matrices.md#s-5-3)。

一个比值需要三个能级，所以具有间距比的最小高斯矩阵是 $3\times3$ 的；对它们，
平均间距比可以精确计算，正如对 $2\times2$ 矩阵可以精确计算 Wigner 猜测一样。

<a id="r-5-14"></a>

**命题 5.14（平均间距比）。** (a) *独立能级。* 设 $s,s'$ 相互独立且服从同一指数
分布（如同由独立的指数变量 $e_j$ 构造的能级 $e_0+\dots+e_{k-1}$ 的两个相邻间距，A.7 节）。则
$r=\min(s,s')/\max(s,s')$ 在 $[0,1]$ 上满足 $P(r\le t)=\frac{2t}{1+t}$，密度为 $\frac2{(1+t)^2}$，并且
$$
\mathbb Er=2\ln2-1\approx0.3863 .
$$
(b) *$3\times3$ GOE 与 GUE。* 设 $H$ 是 $3\times3$ GOE（$\beta=1$）或 GUE（$\beta=2$）矩阵，其特征值为
$\lambda_1<\lambda_2<\lambda_3$，间距为 $s_1=\lambda_2-\lambda_1$，$s_2=\lambda_3-\lambda_2$，并令 $r=\min(s_1,s_2)/\max(s_1,s_2)$。
则 $x=s_2/s_1$ 的密度正比于
$$
f_\beta(x)=\frac{(x+x^2)^\beta}{(1+x+x^2)^{1+3\beta/2}}\qquad\text{on }(0,\infty),
$$
$r$ 在 $(0,1)$ 上的密度为 $f_\beta/I_0$，其中 $I_0=\int_0^1f_\beta$，并且
$$
\mathbb Er=4-2\sqrt3\approx0.5359\quad(\text{GOE}),\qquad
\mathbb Er=\frac{2\sqrt3}\pi-\frac12\approx0.6027\quad(\text{GUE}).
$$

*证明。* (a) 由尺度变换，可设速率为 $1$。对 $x=s/s'$，由独立性和 Tonelli 定理，
$P(x\le t)=\int_0^\infty e^{-y}(1-e^{-ty})\,dy=1-\frac1{1+t}=\frac t{1+t}$。
对 $t<1$，事件 $\{x\le t\}$ 与 $\{1/x\le t\}$ 不相交，$x$ 与 $1/x$ 同分布
（$s,s'$ 可交换），且 $\{r\le t\}$ 是它们的并，所以 $P(r\le t)=\frac{2t}{1+t}$。于是
$\mathbb Er=\int_0^1\frac{2t}{(1+t)^2}dt=2\Big[\ln(1+t)+\frac1{1+t}\Big]_0^1=2\ln2-1$。

(b) *$x$ 的分布。* 由定理 C.7，向量 $\lambda=(\lambda_1,\lambda_2,\lambda_3)$ 在
$C=\{\lambda_1<\lambda_2<\lambda_3\}$ 上具有密度 $Z^{-1}\Delta(\lambda)^\beta e^{-|\lambda|^2/2}$，其中 $\Delta(\lambda)=\prod_{i<j}(\lambda_j-\lambda_i)$。线性映射
$$
(m,s_1,s_2)\mapsto\Big(m-\tfrac{2s_1+s_2}3,\ m+\tfrac{s_1-s_2}3,\ m+\tfrac{s_1+2s_2}3\Big)
$$
的行列式为 $1$，并把 $\mathbb R\times(0,\infty)^2$ 映满 $C$；它的逆为 $m=\frac13\sum_i\lambda_i$，
$s_1=\lambda_2-\lambda_1$，$s_2=\lambda_3-\lambda_2$。在这些变量下 $\Delta=s_1s_2(s_1+s_2)$，且
$|\lambda|^2=3m^2+\frac23(s_1^2+s_1s_2+s_2^2)$（展开即得）。由引理 A.3(c) 和 Fubini 定理（积掉 $m$），
$(s_1,s_2)$ 在 $(0,\infty)^2$ 上具有密度 $c\,\big(s_1s_2(s_1+s_2)\big)^\beta e^{-(s_1^2+s_1s_2+s_2^2)/3}$，其中 $c>0$
为常数。令 $q(x)=1+x+x^2$。对 Borel 集 $A\subseteq(0,\infty)$，在内层积分中作替换 $s_2=xs_1$ 并
利用 Tonelli 定理，
$$
P(x\in A)=c\int_0^\infty\!\!\int_A s_1^{3\beta+1}(x+x^2)^\beta e^{-s_1^2q(x)/3}\,dx\,ds_1
=c\,\gamma_\beta\int_A(x+x^2)^\beta\Big(\frac3{q(x)}\Big)^{1+3\beta/2}dx ,
$$
其中因子来自替换 $s_1=(3/q(x))^{1/2}u$，而 $\gamma_\beta=\int_0^\infty u^{3\beta+1}e^{-u^2}du\in(0,\infty)$。
所以 $x$ 的密度为 $f_\beta/\int_0^\infty f_\beta$。

*化归到 $[0,1]$。* $f_\beta(1/x)=x^2f_\beta(x)$（分子分母同乘 $x$ 的适当幂次）。令
$x=1/y$，对 Borel 集 $A\subseteq(0,1)$ 有 $\int_{1/A}f_\beta(x)\,dx=\int_Af_\beta(1/y)y^{-2}dy=\int_Af_\beta$，其中
$1/A=\{1/y:y\in A\}$；特别地，$\int_0^\infty f_\beta=2I_0$，这里 $I_j=\int_0^1x^jf_\beta(x)\,dx$。由于 $r=\min(x,1/x)$ 且
$P(x=1)=0$，有 $\{r\in A\}=\{x\in A\}\cup\{x\in1/A\}$（不相交），所以
$$
P(r\in A)=\frac{2\int_Af_\beta}{2I_0}=\int_A\frac{f_\beta}{I_0},\qquad \mathbb Er=\frac{I_1}{I_0} .
$$
由于 $x+x^2=q-1$，有 $f_\beta=g_\beta(q)$，其中 $g_1(u)=u^{-3/2}-u^{-5/2}$，$g_2(u)=(u-1)^2u^{-4}$。把
$x$ 写成 $x=\frac12(q'(x)-1)$（其中 $q'=2x+1$），并作替换 $u=q(x)$（从 $q(0)=1$ 到 $q(1)=3$）：
$$
I_1=\frac12\int_1^3g_\beta(u)\,du-\frac12I_0,\qquad \mathbb Er=\frac{\int_1^3g_\beta}{2I_0}-\frac12 .
$$
这里 $\int_1^3g_1=\big[-2u^{-1/2}+\frac23u^{-3/2}\big]_1^3=\frac43-\frac{16}{9\sqrt3}$，
$\int_1^3g_2=\big[-u^{-1}+u^{-2}-\frac13u^{-3}\big]_1^3=\frac8{81}$。

*积分 $I_0$。* 令 $L_a=\int_0^1q^{-a}dx$。由 $(2x+1)^2=4q-3$，
$\frac d{dx}\big[(2x+1)q^{-a}\big]=2q^{-a}-a(2x+1)^2q^{-a-1}=(2-4a)q^{-a}+3aq^{-a-1}$，在 $[0,1]$ 上积分得
$$
3a\,L_{a+1}=3^{1-a}-1+(4a-2)L_a .
$$
GOE：$f_1=q^{-3/2}-q^{-5/2}$。取 $a=\frac12$ 得 $L_{3/2}=\frac23(\sqrt3-1)$；取 $a=\frac32$ 得
$L_{5/2}=\frac29\big(\frac1{\sqrt3}-1+4L_{3/2}\big)=\frac{18\sqrt3-22}{27}$。所以 $I_0=L_{3/2}-L_{5/2}=\frac4{27}$，且
$\mathbb Er=\frac{27}8\big(\frac43-\frac{16}{9\sqrt3}\big)-\frac12=\frac92-2\sqrt3-\frac12=4-2\sqrt3$。
GUE：$f_2=q^{-2}-2q^{-3}+q^{-4}$。由于 $q=(x+\frac12)^2+\frac34$，
$L_1=\frac2{\sqrt3}\big[\arctan\frac{2x+1}{\sqrt3}\big]_0^1=\frac2{\sqrt3}\big(\frac\pi3-\frac\pi6\big)=\frac\pi{3\sqrt3}$，而
取 $a=1,2,3$ 的递推关系给出 $L_2=\frac23L_1=\frac{2\pi}{9\sqrt3}$，$L_3=L_2-\frac19$，$L_4=\frac19\big(10L_3-\frac89\big)$。所以
$I_0=L_2-2L_3+L_4=\frac{2\pi}{81\sqrt3}$，且 $\mathbb Er=\frac{81\sqrt3}{4\pi}\cdot\frac8{81}-\frac12=\frac{2\sqrt3}\pi-\frac12$。$\square$

(b) 中的数值由 Y. Y. Atas、E. Bogomolny、O. Giraud 和 G. Roux 得到（Phys. Rev. Lett. 110
(2013), 084101）。对更大的 $n$，本课程不计算平均间距比。

(a) 部分对区间上 $n$ 个独立均匀分布的能级精确成立：它们的 $n+1$ 个间距（包括两端的
两个空隙）与 $(e_0,\dots,e_n)/\sum_je_j$ 同分布，其中 $e_j$ 是独立的指数变量，因此相邻间距的每个
比值 $s_i/s_{i+1}$ 与 $e_i/e_{i+1}$ 同分布，而每个间距比 $r_i$ 与 (a) 中的 $r$ 同分布
（定理 A.25 与推论 A.26(b)、(d)）。

*计算：* 用蒙特卡罗方法和对 $f_\beta$ 的数值积分检验命题 5.14，以及数值观察到的 $n=64$ 的 GOE 的平均间距比——[教程，5.3 节](../../tutorial/zh/05_random_matrices.md#s-5-3)。

<a id="s-5-4"></a>

## 5.4 为何不能混合扇区

设 $H=\bigoplus_cH_c\in\mathcal M$ 是一个具有**相互独立的块**的随机自伴元素，
例如 GOE 块。作为 $\mathbb C^D$ 上的算子，它的谱是各块谱的
并。在一个块内部能级相互排斥；我们将证明不同块的能级之间并不排斥。

假设每个块都有有界连续的**能级密度** $\rho_c$：
对有界 Borel 函数 $\varphi$ 有 $\mathbb E\sum_i\varphi(\lambda_i(H_c))=\int\varphi\,\rho_c$（从而
$\int\rho_c=k_c$）。对 GOE 块和 GUE 块，这由推论 C.8 保证，它从定理 C.7 的联合
密度导出 $\rho_c$。

<a id="r-5-15"></a>

**定理 5.15（扇区之间没有排斥）。** 对 $c\neq d$ 以及每个有界 Borel 函数
$\psi:\mathbb R^2\to\mathbb R$，
$$
\mathbb E\sum_{i,j}\psi\big(\lambda_i(H_c),\lambda_j(H_d)\big)=\iint\psi(x,y)\,\rho_c(x)\rho_d(y)\,dx\,dy .
$$
因此，若 $N_\varepsilon$ 是来自不同块、距离
$<\varepsilon$ 的能级对的个数，则当 $\varepsilon\to0$ 时 $\mathbb EN_\varepsilon/(2\varepsilon)\to\sum_{c<d}\int\rho_c\rho_d$，
并且合并谱的最小间距 $\delta$ 满足
$$
P(\delta<\varepsilon)\ \ge\ \frac{\mathbb EN_\varepsilon}{\sum_{c<d}k_ck_d}
=\frac{2\varepsilon\sum_{c<d}\int\rho_c\rho_d}{\sum_{c<d}k_ck_d}+o(\varepsilon).
$$

*证明。* 排好序的特征值是矩阵的连续函数（推论 A.11），所以
$\Lambda_c=(\lambda_i(H_c))_i$ 与 $\Lambda_d$ 是相互独立的随机向量。由 Fubini 定理（它们的联合分布是
各自分布的乘积），$\mathbb E\sum_{i,j}\psi(\Lambda_{c,i},\Lambda_{d,j})=\mathbb E\sum_j g(\Lambda_{d,j})$，其中
$g(y)=\mathbb E\sum_i\psi(\Lambda_{c,i},y)=\int\psi(x,y)\rho_c(x)\,dx$ 是有界 Borel 函数；对 $g$ 应用
$\rho_d$ 的定义即可。取 $\psi(x,y)=1\{|x-y|<\varepsilon\}$，
$$
\frac{\mathbb EN_\varepsilon}{2\varepsilon}=\sum_{c<d}\int\rho_c(x)\Big(\frac1{2\varepsilon}\int_{x-\varepsilon}^{x+\varepsilon}\rho_d(y)\,dy\Big)dx ;
$$
括号内的量由连续性趋于 $\rho_d(x)$，且以 $\sup\rho_d$ 为界，而 $\rho_c$ 可积，
所以由控制收敛定理得到该极限。最后 $N_\varepsilon\le\sum_{c<d}k_ck_d$，所以
$\mathbb EN_\varepsilon\le\big(\sum_{c<d}k_ck_d\big)P(N_\varepsilon\ge1)$，而 $N_\varepsilon\ge1$ 蕴含 $\delta<\varepsilon$。$\square$

因此，对两个 GOE 或 GUE 块，$P(\delta<\varepsilon)$ 至少是 $\varepsilon$ 的量级：它们的能级密度
处处为正（在推论 C.8 中，除去 $(x,y)$ 的某两个坐标相等的零测集之外，$q(x,\cdot)>0$），
所以 $\int\rho_c\rho_d>0$。而在一个 $2\times2$ GOE 块内部它等于 $1-e^{-\varepsilon^2/4}\le\varepsilon^2/4$
（命题 5.13，未归一化的间距），在一个 $n\times n$ GOE 块内部它至多是常数乘以
$\varepsilon^2$（命题 C.9）。尽管每个块都是完美的 GOE，
合并谱却通不过能级排斥检验。

*计算：* 对 $C$ 个独立的 GOE 块，随着 $C$ 增大，合并谱的平均间距比下降并趋向独立能级的值 $2\ln2-1$；带随机场的 Heisenberg 链（第 9 章）的各块并不独立，却显示出同样的效应；二者都是数值观察——[教程，5.4 节](../../tutorial/zh/05_random_matrices.md#s-5-4)。

> **物理意义**（物理解释，而非本课程的定理）。具有不同守恒量子数的能级之间没有相互作用，
> 所以它们可以自由交叉；只有同一对称性的能级才相互排斥。因此，混合了对称扇区的谱
> 即使在每个扇区内的动力学完全混沌时，也会呈现类 Poisson 统计，
> 而一个未被分辨的对称性可能被误认为可积性或局域化。业内的经验
> 规则——在计算能级统计之前先分辨所有对称性——说的正是：应当
> 在代数 $\bigoplus_cM_{k_c}$ 中计算，而不是在 $M_D$ 中计算。

<a id="s-5-5"></a>

## 5.5 可与之比较的工具

*库中：* 与其他软件的比较——[教程，5.5 节](../../tutorial/zh/05_random_matrices.md#s-5-5)。

<a id="x-5"></a>

## 习题

1. 设 $U$ 在 $U(\mathcal M)$ 上服从 Haar 分布，$\mathrm{Tr}\,U=\sum_c\mathrm{Tr}\,U_c$。证明
   $\mathbb E|\mathrm{Tr}\,U|^2=C$ 且 $\mathbb E|\mathrm{Tr}\,U|^4=\sum_cm(k_c)+2C(C-1)$，其中
   $m(1)=1$，当 $k\ge2$ 时 $m(k)=2$；并与 $U(D)$（$D\ge2$）上 Haar 酉矩阵的 $1$ 和 $2$ 比较，再对
   $\mathcal M=M_1\oplus M_2\oplus M_3$ 求出其值。
2. 利用定理 5.12，证明当 $n\le k$ 时 $\mathbb E|U_{11}|^{2n}=\frac{n!}{k(k+1)\cdots(k+n-1)}$
   （推论 B.14 把它推广到 $n>k$）。
3. **圆正交系综**是 $S=W^TW$，其中 $W$ 在 $U(k)$ 上服从 Haar 分布。证明
   $S$ 是酉且对称的，并且它的分布对每个 $V\in U(k)$ 在 $S\mapsto V^TSV$ 下不变。设
   $\lambda_1,\dots,\lambda_k$ 为 $S$ 的特征值（计入重数），对单位圆 $\mathbb T$ 的 Borel 集 $B$ 令
   $N_B(S)=\#\{j:\lambda_j\in B\}$。证明 $N_B(S)$ 是随机变量，$\nu(B)=\frac1k\mathbb EN_B(S)$ 是归一化的
   弧长（本征相位的平均密度是均匀的），并且对每个整数 $m\neq0$ 有 $\mathbb E\mathrm{Tr}\,S^m=0$。

<details>
<summary>解答</summary>

1. 记 $t_c=\mathrm{Tr}\,U_c$；各 $t_c$ 相互独立（推论 5.4）。把 $U_c$ 乘以
   $e^{i\theta}$ 不改变其分布，并把 $t_c$ 乘以 $e^{i\theta}$；所以矩
   $\mathbb E[t_at_b\bar t_c\bar t_d]$ 为零，除非每个扇区在 $(a,b)$ 中出现的次数与在
   $(c,d)$ 中相同，即作为多重集 $\{a,b\}=\{c,d\}$；同理，当 $a\neq b$ 时 $\mathbb E[t_a\bar t_b]=0$。由
   命题 5.7，$\mathbb E[(U_c)_{ii}\overline{(U_c)_{jj}}]=\langle e_i,\mathbb E[U_ce_ie_j^*U_c^*]e_j\rangle=\delta_{ij}/k_c$，
   所以 $\mathbb E|t_c|^2=1$，$\mathbb E|\mathrm{Tr}\,U|^2=C$。在 $\mathbb E|\sum_at_a|^4$ 中，
   $a=b$ 的项给出 $\sum_a\mathbb E|t_a|^4$，$a\neq b$ 的项给出
   $2\sum_{a\neq b}\mathbb E|t_a|^2\,\mathbb E|t_b|^2=2C(C-1)$（$(c,d)$ 有两种选法）。对 $k=1$，$|t|=1$。
   对 $k\ge2$，$\mathbb E|t|^4=\sum_{i,j,i',j'}\mathbb E[U_{ii}U_{jj}\overline{U_{i'i'}U_{j'j'}}]$；在
   定理 5.10 中，$\sigma=\pi$ 的两项以同一种方式迫使 $\{i',j'\}=\{i,j\}$，求和后各给出
   $k^2\mathrm{Wg}(e)$；$\sigma\neq\pi$ 的两项迫使 $i=j=i'=j'$，各给出
   $k\,\mathrm{Wg}(\tau)$：总计 $2\big(\frac{k^2}{k^2-1}-\frac1{k^2-1}\big)=2$。$C=1$，
   $k=D$ 的情形对 $U(D)$ 给出 $1$ 和 $2$。对 $M_1\oplus M_2\oplus M_3$：
   $\mathbb E|\mathrm{Tr}\,U|^4=1+2+2+12=17$，而 $U(6)$ 上的 Haar 酉矩阵为 $2$。

2. 当所有指标都等于 $1$ 时，每个 δ 都是 $1$，所以
   $\mathbb E|U_{11}|^{2n}=\sum_{\sigma,\pi}\mathrm{Wg}(\sigma\pi^{-1})=\mathbf 1^TW\mathbf 1$，其中
   $W=\Gamma^{-1}=(\mathrm{Wg}(\sigma^{-1}\pi))$，$\Gamma=(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)})$。这里用到
   $\mathrm{Wg}$ 只依赖于循环型，以及 $\sigma\pi^{-1}=\sigma(\pi^{-1}\sigma)\sigma^{-1}$
   共轭于 $\pi^{-1}\sigma=(\sigma^{-1}\pi)^{-1}$，后者与 $\sigma^{-1}\pi$ 具有相同的循环型。$\Gamma$ 的每个
   行和都是 $c_n(k)=\sum_{\rho\in S_n}k^{\#\mathrm{cycles}(\rho)}$，所以 $\Gamma\mathbf 1=c_n(k)\mathbf 1$，
   $W\mathbf 1=\mathbf 1/c_n(k)$，$\mathbf 1^TW\mathbf 1=n!/c_n(k)$。最后，由归纳法
   $c_n(x)=x(x+1)\cdots(x+n-1)$：$\{1,\dots,n\}$ 的一个置换可由
   $\{1,\dots,n-1\}$ 的一个置换得到，或者把 $n$ 作为不动点加入（多一个循环，因子 $x$），或者把 $n$
   在其循环中插入到 $n-1$ 个元素之一的后面（$n-1$ 种方式，循环数不变）；这是一个双射，所以
   $c_n(x)=(x+n-1)c_{n-1}(x)$。

3. $S^*S=W^*\bar WW^TW=1$，因为 $W^T$ 是酉的；并且 $S^T=W^TW=S$。对 $V\in U(k)$，$WV$ 与 $W$
   同分布（右不变性，引理 5.1），而 $(WV)^T(WV)=V^TSV$。取 $V=e^{i\alpha/2}1$：对每个实数 $\alpha$，
   $e^{i\alpha}S$ 与 $S$ 同分布。

   *可测性。* 由定理 1.3，$S=Y\mathrm{diag}(\lambda_1,\dots,\lambda_k)Y^*$，其中 $Y$ 是酉矩阵。若 $v$ 是
   满足 $Sv=\mu v$ 的单位向量，则把 $v$ 写成 $v=\sum_jc_jYe_j$ 得 $c_j(\lambda_j-\mu)=0$，所以只要 $\mu\in K$，$v$ 就属于
   $E_K=\mathrm{span}\{Ye_j:\lambda_j\in K\}$；并且 $\dim E_K=N_K(S)$。因此，对 $K\subseteq\mathbb T$
   和 $m\ge1$：$N_K(S)\ge m$ 当且仅当 $S$ 有 $m$ 个规范正交的特征向量，其特征值属于 $K$。若 $K$
   是闭集，则具有这一性质的酉矩阵构成的集合 $F_{K,m}$ 在 $U(k)$ 中是闭的：若 $S_l\to S$，且 $S_l$ 有
   规范正交特征向量 $v_i^{(l)}$，其特征值 $\mu_i^{(l)}\in K$（$i\le m$），则存在子列使得对所有 $i\le m$ 有
   $v_i^{(l)}\to v_i$ 且 $\mu_i^{(l)}\to\mu_i\in K$（在紧集 $(\text{单位球面}\times K)^m$ 中应用
   Bolzano–Weierstrass 定理），而极限是规范正交的，且 $Sv_i=\mu_iv_i$。
   设 $0\le\theta_1(S)\le\dots\le\theta_k(S)<2\pi$ 为本征相位，于是特征值就是诸
   $e^{i\theta_j(S)}$（推论 1.4：$|\lambda_j|=1$）。对 $0\le t<2\pi$，
   弧 $K_t=\{e^{i\varphi}:0\le\varphi\le t\}$ 是闭的，且 $\{\theta_j\le t\}=\{N_{K_t}\ge j\}=F_{K_t,j}$，所以每个
   $\theta_j$ 都是 $U(k)$ 上的 Borel 函数。因此 $N_B(S)=\sum_j1_B(e^{i\theta_j(S)})$ 是随机变量，而
   $\nu(B)=\frac1k\sum_jP(e^{i\theta_j(S)}\in B)$ 是 $\mathbb T$ 上的 Borel 概率测度，即 $k$ 个特征值的
   分布的平均。

   *均匀性。* $e^{i\alpha}S$ 的特征值是诸 $e^{i\alpha}\lambda_j$，所以 $N_B(e^{i\alpha}S)=N_{e^{-i\alpha}B}(S)$；
   又由于 $e^{i\alpha}S$ 与 $S$ 同分布，$\nu(e^{-i\alpha}B)=\frac1k\mathbb EN_B(e^{i\alpha}S)=\nu(B)$。所以 $\nu$
   在旋转下不变，由命题 A.22（取 $N=1$）它就是归一化的弧长：
   对每段长度为 $|B|$ 的弧，$\mathbb EN_B(S)=k\,|B|/2\pi$。同样，$\mathrm{Tr}\,S^m=\sum_j\lambda_j^m$（对
   $m<0$ 也成立，因为 $S^{-1}=S^*$）是 $S$ 的有界连续函数，并且对所有 $\alpha$ 有
   $\mathbb E\mathrm{Tr}\,S^m=\mathbb E\mathrm{Tr}(e^{i\alpha}S)^m=e^{im\alpha}\,\mathbb E\mathrm{Tr}\,S^m$，所以当
   $m\neq0$ 时它为零：这些量是 $\nu$ 的 Fourier 系数的 $k$ 倍。特别地，平均
   本征相位密度是常数，绕圆周的间距不需要展开。

</details>

*计算：* 习题 1–3 的数值部分——[教程，习题](../../tutorial/zh/05_random_matrices.md#x-5)。

**下一章：** [第 6 章](06_quantum_channels.md)讨论量子信道：Kraus 表示与 Stinespring
表示、关于迹的对偶，以及扇区之间的映射。本章中作为 Haar 平均出现的条件
期望 $E$ 是第一个例子。
