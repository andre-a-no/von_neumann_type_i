[English](../05_random_matrices.md) · [Español](../es/05_random_matrices.md) · [Français](../fr/05_random_matrices.md) · [Deutsch](../de/05_random_matrices.md) · **中文** · [日本語](../ja/05_random_matrices.md) · [Русский](../ru/05_random_matrices.md)

# 5. 代数中的随机矩阵

随机矩阵理论关心"典型"算子是什么样子。在代数 $\bigoplus_cM_{k_c}$ 中，答案必须尊重扇区：
自然的随机酉算子是该代数的酉元，而不是大矩阵代数 $M_D$ 的酉矩阵；自然的随机哈密顿量
具有相互独立的块。本章构造代数酉群上的 Haar 测度，精确计算它的一阶矩与二阶矩，然后说明
忽略扇区时谱统计会出现什么问题。

**你将学到**为什么高斯矩阵的 QR 分解给出 Haar 酉矩阵（以及为什么必须固定相位），在代数的
酉群上求平均就是到中心上的条件期望，二阶 Weingarten 公式，Wigner 猜测与间距比，以及为什么
在混合了扇区的谱上计算的能级统计看起来像独立能级的统计。

在本章中，$\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$，$D=\sum_ck_c$，$P_c$ 为
极小中心投影，$E$ 是命题 2.8 中到中心上的条件期望。$U(k)$ 是 $k\times k$ 酉矩阵构成的群，
$U(\mathcal M)=\{U\in\mathcal M:U^*U=1\}$ 是 $\mathcal M$ 的酉群。随机变量定义在一个
概率空间上，期望记为 $\mathbb E$；随机矩阵的期望按元素逐个取。*标准复高斯*变量的实部与虚部
相互独立且 $\sim N(0,\tfrac12)$，从而 $\mathbb E|z|^2=1$（这正是 `torch.randn` 对复数
dtype 抽取的分布）；*Ginibre 矩阵*的元素是相互独立的标准复高斯变量。

## 5.1 Haar 测度

设 $G$ 是一个紧矩阵群（某个 $U(n)$ 的闭子群）。$G$ 上的 Borel 概率测度
$\mu$ 称为**左不变的**，如果对所有 $g\in G$ 和 Borel 集 $B$ 有 $\mu(gB)=\mu(B)$，
等价地，对每个有界 Borel 函数 $f$ 有 $\int f(gh)\,d\mu(h)=\int f\,d\mu$；
**右不变**用 $Bg$ 和 $f(hg)$ 类似地定义。

**引理 5.1（唯一性）。** 若 $\mu$ 是 $G$ 上左不变的 Borel 概率测度，$\nu$ 是右不变的 Borel
概率测度，则 $\mu=\nu$。因此 $G$ 上至多有一个左不变的 Borel
概率测度；若它存在，则它也是右不变的，并且在 $g\mapsto g^{-1}$ 下不变。

*证明。* 设 $f:G\to\mathbb R$ 连续（从而有界）。由 Fubini 定理，
$\int\!\big(\int f(gh)\,d\mu(h)\big)d\nu(g)=\int\!\big(\int f(gh)\,d\nu(g)\big)d\mu(h)$。左边的内层
积分由左不变性等于 $\int f\,d\mu$，右边的内层积分由右不变性等于 $\int f\,d\nu$；
所以 $\int f\,d\mu=\int f\,d\nu$。度量空间上两个 Borel 概率测度若对连续函数的积分都相同，
则二者相等：连续函数
$\max(0,1-n\operatorname{dist}(\cdot,K))$ 递减地收敛到闭集 $K$ 的示性函数，所以两测度
在闭集上一致（控制收敛定理），而闭集生成 Borel 集且对交封闭
（Dynkin 的 $\pi$–$\lambda$ 定理）。若 $\mu_1,\mu_2$ 都是左不变的，则
$\nu(B)=\mu_2(B^{-1})$ 是右不变的，因为 $(Bg)^{-1}=g^{-1}B^{-1}$；因此 $\mu_1=\nu=\mu_2$，
所以 $\mu_1$ 是右不变且在取逆下不变的。$\square$

**定义。** $G$ 的 **Haar 测度**是它的左不变 Borel 概率测度。以此为分布的
随机元素称为**服从 Haar 分布**。

每个紧群都有 Haar 测度（见 G. B. Folland, *A Course in Abstract Harmonic
Analysis*, 2nd ed., CRC Press, 2016, Section 2.2；此处不加证明地引用）。对于 $U(k)$，我们
不需要这个一般定理：定理 5.3 直接构造出该测度。

**引理 5.2（QR 分解）。** 每个可逆的 $Z\in M_k(\mathbb C)$ 都可以写成 $Z=QR$，
其中 $Q$ 是酉矩阵，$R$ 是对角元为正的上三角矩阵。$Q$ 与 $R$ 是唯一的，并且
连续地依赖于 $Z$。

*证明。* 对 $Z$ 的线性无关的列 $z_1,\dots,z_k$ 施行 Gram–Schmidt 正交化：
$q_j=w_j/\|w_j\|$，其中 $w_j=z_j-\sum_{i<j}\langle q_i,z_j\rangle q_i$；这里 $w_j\neq0$，因为
$z_j\notin\operatorname{span}\{z_1,\dots,z_{j-1}\}=\operatorname{span}\{q_1,\dots,q_{j-1}\}$。于是
$z_j=\sum_{i<j}\langle q_i,z_j\rangle q_i+\|w_j\|q_j$，即 $Z=QR$，其中 $R_{jj}=\|w_j\|>0$，并且
每一步在可逆矩阵上都是连续的。若 $Q_1R_1=Q_2R_2$，则
$T=Q_2^*Q_1=R_2R_1^{-1}$ 是对角元为正的酉上三角矩阵。它的逆
$T^*$ 是下三角的，同时也是上三角的（上三角矩阵的逆是上三角的），
所以 $T$ 是对角的；对角元为正的对角酉矩阵就是 $1$。$\square$

**定理 5.3（由高斯矩阵得到 Haar 酉矩阵）。** 设 $Z$ 是 $k\times k$
Ginibre 矩阵。则 $Z$ 以概率一可逆，并且引理 5.2 中的因子 $Q$ 在 $U(k)$ 上
服从 Haar 分布。此外，若 $Z=Q'R'$ 是任意一个 QR 分解，其中 $Q'$ 为酉矩阵、$R'$
为上三角矩阵（例如 LAPACK 返回的分解，其对角元为实数，符号可正可负），则
$Q=Q'\Lambda$（其中 $\Lambda=\operatorname{diag}(R'_{ii}/|R'_{ii}|)$）服从 Haar 分布。

*证明。* *可逆性。* $\det Z$ 是 $Z$ 的 $2k^2$ 个实坐标的非零多项式，
而非零实多项式的零点集的 Lebesgue 测度为零（对变量个数归纳：
写 $p=\sum_jp_j(x')x_n^j$，其中某个 $p_j\neq0$；对 $x'$ 位于 $p_j$ 的零测零点集之外的情形，
$x_n$ 只有有限多个根；再由 Fubini 定理得出结论）。$Z$ 有密度，所以
$P(\det Z=0)=0$。

*高斯分布的不变性。* $Z$ 关于
$M_k(\mathbb C)\cong\mathbb R^{2k^2}$ 上 Lebesgue 测度的密度为 $\pi^{-k^2}e^{-\|Z\|_2^2}$。对 $V\in U(k)$，映射
$Z\mapsto VZ$ 是实线性的且保持 $\|Z\|_2$（1.4 节），所以它是
$\mathbb R^{2k^2}$ 的正交映射，Jacobi 行列式为 $1$，并且保持密度：$VZ$ 与 $Z$ 同分布。

*等变性。* 若 $Z=QR$，则 $VZ=(VQ)R$ 是 $VZ$ 的引理 5.2 意义下的分解，所以
$Q(VZ)=VQ(Z)$。因此对每个 Borel 集 $B\subseteq U(k)$，
$P(VQ(Z)\in B)=P(Q(VZ)\in B)=P(Q(Z)\in B)$：$Q(Z)$ 的分布是左不变的，由
引理 5.1 它就是 Haar 测度（这也证明了 $U(k)$ 的 Haar 测度存在）。
最后 $Q'R'=(Q'\Lambda)(\Lambda^*R')$，而 $\Lambda^*R'$ 是对角元为
$|R'_{ii}|>0$ 的上三角矩阵，所以由唯一性得 $Q'\Lambda=Q(Z)$。$\square$

对实高斯矩阵作同样的证明，就得到正交群 $O(k)$ 的 Haar 测度。
若不作修正 $\Lambda$，结果就*不*服从 Haar 分布：LAPACK 对 $R'_{ii}$ 符号的约定
会引入偏差，下面的代码展示了这一点
（F. Mezzadri, *How to generate random matrices from the classical compact groups*, Notices AMS 54
(2007), 592–604）。

**推论 5.4（代数的 Haar 测度）。** $U(\mathcal M)=\{\bigoplus_cU_c:U_c\in U(k_c)\}\cong\prod_cU(k_c)$，
并且 $U=\bigoplus_cU_c$ 在 $U(\mathcal M)$ 上服从 Haar 分布当且仅当各块 $U_c$ 相互独立，
且每个 $U_c$ 在 $U(k_c)$ 上服从 Haar 分布。

*证明。* $U^*U=\bigoplus_cU_c^*U_c$，所以 $U$ 是酉的当且仅当每个块都是酉的。设各 $U_c$
相互独立且服从 Haar 分布。对 $V=\bigoplus_cV_c\in U(\mathcal M)$，$VU$ 的各块 $V_cU_c$
相互独立（独立变量的函数），且 $V_cU_c$ 与 $U_c$ 同分布，所以
$VU$ 与 $U$ 同分布。于是 $U$ 的分布是左不变的，由引理 5.1 它就是紧群
$U(\mathcal M)$ 的 Haar 测度。反之，由唯一性，每个服从 Haar 分布的 $U$ 都具有
这一分布。$\square$

群 $U(\mathcal M)$ 比 $U(D)$ 小得多：当 $C\ge2$ 时，$U(D)$ 的 Haar 酉矩阵以概率零
落在 $\mathcal M$ 中（由定理 5.3 以及引理 5.2 中的 Gram–Schmidt 构造，
它的第一列是 $z/\|z\|$，其中 $z$ 是高斯向量，所以它的所有元素几乎必然非零）。
$SU(k)$ 的 Haar 测度可以用同样的思路得到：若 $U$ 在
$U(k)$ 上服从 Haar 分布且 $\theta=\arg\det U$，则 $V=e^{-i\theta/k}U\in SU(k)$，并且对 $W\in SU(k)$，
$\det(WU)=\det U$ 给出 $V(WU)=WV(U)$；由于 $WU$ 与 $U$ 同分布，$V$ 的分布在 $SU(k)$ 上是
左不变的，因此由引理 5.1 它就是 $SU(k)$ 的 Haar 测度。

**推论 5.5（元素的一阶矩）。** 若 $U$ 在 $U(k)$ 上服从 Haar 分布，则
$\mathbb EU_{ij}=0$ 且 $\mathbb E|U_{ij}|^2=1/k$。

*证明。* $-1\in U(k)$，所以 $-U$ 与 $U$ 同分布，$\mathbb EU=-\mathbb EU$。对置换
矩阵 $P$，$PU$ 与 $U$ 同分布，所以 $\mathbb E|U_{ij}|^2$ 与 $i$ 无关；由于
$\sum_i|U_{ij}|^2=1$，它等于 $1/k$。$\square$

本库用 `random_unitary_operator` 抽取 $U(\mathcal M)$ 上的 Haar 酉元：对每种块大小，
一次批量 QR 分解同时处理所有扇区和所有样本（在 GPU 上则用 Gram–Schmidt，由引理 5.2 中的
唯一性，它给出相同的 $Q$）。测度 `'haar'` 就是在每个块中应用定理 5.3；对复代数，
`'haar_su'` 如上乘以 $e^{-i\theta/k}$；`'coe'` 和
`'cse'` 是 Dyson 的圆系综（前者见习题 3）。我们检验推论 5.4 和 5.5，
并检验本库的采样器恰好就是定理 5.3 中经过修正的 QR：

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, SpinChain, channels

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def mc(samples, exact, name):
    """Monte Carlo mean of real samples with its standard error; must agree within 5 standard errors."""
    mean, err = samples.mean().item(), samples.std().item() / math.sqrt(samples.numel())
    print(f"{name}: {mean:.4f} +- {err:.4f}   (theory {exact:.4f})")
    assert abs(mean - exact) <= 5 * err + 1e-12
    return mean

alg = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')        # M_1 (+) M_2 (+) M_3, D = 6
B = 100_000
U = alg.random_unitary_operator(B)                                   # Haar on U(1) x U(2) x U(3)
assert close(U.adjoint() @ U, alg.identity())                        # unitary elements of the algebra
u = [U.matrix[:, c, :k, :k] for c, k in enumerate(alg.k_factors)]   # the blocks U_c
for c, k in enumerate(alg.k_factors):                                # Corollary 5.5
    mc(u[c][:, 0, 0].real, 0.0, f"E Re (U_{c})_11")
    mc(u[c][:, 0, 0].abs() ** 2, 1 / k, f"E |(U_{c})_11|^2")
mc(u[1][:, 0, 0].abs() ** 2 * u[2][:, 0, 0].abs() ** 2, 1 / 6, "E |(U_1)_11|^2 |(U_2)_11|^2")   # independent blocks

SU = alg.random_unitary_operator(4, measure='haar_su').matrix
assert torch.allclose(torch.linalg.det(SU[:, 2, :3, :3]), torch.ones(4, dtype=dt))

torch.manual_seed(1)
Z = torch.randn(B, 2, 2, dtype=dt)                                   # Ginibre matrices
torch.manual_seed(1)
V = alg.random_unitary(2, batch_size=B)                              # the library, same random numbers
Q, R = torch.linalg.qr(Z)                                            # LAPACK: diag(R) real, either sign
d = R.diagonal(dim1=-2, dim2=-1)
assert close(V, Q * (d / d.abs()).unsqueeze(-2))                     # Theorem 5.3: V = Q Lambda
raw = Q.diagonal(dim1=-2, dim2=-1).sum(-1).real                      # Re Tr Q without the correction
print(f"E Re Tr Q without the phase correction: {raw.mean():.4f} +- {raw.std() / math.sqrt(B):.4f}")
assert abs(raw.mean()) > 20 * raw.std() / math.sqrt(B)               # far from 0: not Haar
```

若不作修正，对 $k=2$ 有 $\mathbb E\operatorname{Re}\operatorname{Tr}Q\approx-0.84$，而
Haar 酉矩阵满足 $\mathbb EU=0$。

## 5.2 矩与 Weingarten 演算

**引理 5.6.** 若 $T\in M_k(\mathbb C)$ 与所有对角酉矩阵和所有置换
矩阵可交换，则 $T\in\mathbb C1$。特别地，$U(k)$ 的交换子是 $\mathbb C1$。

*证明。* 对 $\Delta=\operatorname{diag}(e^{i\theta_1},\dots,e^{i\theta_k})$，其中各
$e^{i\theta_j}$ 互不相同，$\Delta T=T\Delta$ 给出 $(e^{i\theta_i}-e^{i\theta_j})T_{ij}=0$，所以 $T$ 是
对角的。对置换矩阵 $P_\sigma e_j=e_{\sigma(j)}$，$P_\sigma TP_\sigma^{-1}=T$ 给出
$T_{\sigma(j)\sigma(j)}=T_{jj}$，所以所有对角元都相等。$\square$

**命题 5.7（Haar twirl 就是到中心上的条件期望）。** 设 $U$ 在
$U(\mathcal M)$ 上服从 Haar 分布。则对每个 $X\in\mathcal M$，
$$
\mathbb E\,[UXU^*]=E(X)=\bigoplus_c\frac{\operatorname{Tr}X_c}{k_c}\,1_{k_c}.
$$

*证明。* $(UXU^*)_c=U_cX_cU_c^*$，且 $U_c$ 在 $U(k_c)$ 上服从 Haar 分布（推论 5.4）。令
$T_c=\mathbb E[U_cX_cU_c^*]$（各元素有界）。对 $V\in U(k_c)$，$VU_c$ 与
$U_c$ 同分布，所以 $VT_cV^*=\mathbb E[(VU_c)X_c(VU_c)^*]=T_c$。由引理 5.6，$T_c=\lambda1$，并且
$k_c\lambda=\operatorname{Tr}T_c=\mathbb E\operatorname{Tr}(U_cX_cU_c^*)=\operatorname{Tr}X_c$。$\square$

这里的归一化与命题 2.8 以及 `channels.center_expectation` 相同。若
$X\in\mathcal M$ 与 $U(\mathcal M)$ 可交换，则 $UXU^*=X$，所以 $X=E(X)$ 是中心元；Haar
平均是第 2 章习题 3 中群平均的连续类比。若 $U$ 改为在 $U(D)$ 上服从 Haar 分布，
同样的证明（取 $C=1$）给出 $\mathbb E[UXU^*]=\frac{\operatorname{Tr}X}{D}1$：
扇区权重丢失了，而且 $UXU^*$ 甚至不在 $\mathcal M$ 中。

> **物理意义。** 代数的 Haar 随机酉元是一个遵守超选择定则的随机演化。
> 对它求平均会把每个扇区内部的一切打乱，而恰好保留经典部分：每个荷的
> 权重。这就是量子信息中的"twirl"，它正是第 2 章中到中心上的条件期望。

对于二阶矩，我们需要 $\{V\otimes V\}$ 在 $\mathbb C^k\otimes\mathbb C^k$ 上的交换子。
设 $F$ 为**交换算子**（swap），$F(x\otimes y)=y\otimes x$；它是自伴且酉的，并且
$\operatorname{Tr}1=k^2$，$\operatorname{Tr}F=\sum_{i,j}\langle e_i\otimes e_j,e_j\otimes e_i\rangle=k$，$F^2=1$。

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

**推论 5.9（二阶 Schur–Weyl 对偶）。** 设 $k\ge2$，$P_\pm=\frac12(1\pm F)$。
它们是到对称子空间与反对称子空间
$\mathrm{Sym}=\{x:Fx=x\}$ 和 $\mathrm{Alt}=\{x:Fx=-x\}$ 上的投影，维数分别为 $k(k+1)/2$ 和 $k(k-1)/2$。
$U(k)$ 的表示 $V\mapsto V\otimes V$ 在 $\mathrm{Sym}$ 上和在 $\mathrm{Alt}$ 上都是不可约的，
并且这两个表示不等价。

*证明。* $F=F^*$ 和 $F^2=1$ 使 $P_\pm$ 成为投影，满足 $P_++P_-=1$，值域分别为 $\mathrm{Sym}$、
$\mathrm{Alt}$；它们的维数是 $\operatorname{Tr}P_\pm=(k^2\pm k)/2$。
$\{V\otimes V\}$ 的线性张成 $\mathcal A$ 是一个 $*$-代数（$1\otimes1$、乘积 $VW\otimes VW$、伴随 $V^*\otimes V^*$），
由引理 5.8，其交换子为 $\mathcal A'=\operatorname{span}\{1,F\}=\operatorname{span}\{P_+,P_-\}$。设
$K\subseteq\mathrm{Sym}$ 在所有 $V\otimes V$ 下不变。由引理 2.1，它的投影属于 $\mathcal A'$，所以它
等于 $xP_++yP_-$，其中 $x,y\in\{0,1\}$（因为它是投影）；它在 $\mathrm{Alt}$ 上为零，所以 $y=0$，
从而 $K=0$ 或 $K=\mathrm{Sym}$。同样的论证适用于 $\mathrm{Alt}$。若
$T:\mathrm{Sym}\to\mathrm{Alt}$ 是这两个表示之间的缠结算子，则 $TP_+$（看作
$\mathbb C^k\otimes\mathbb C^k$ 上的算子）与 $\mathcal A$ 可交换，所以 $TP_+=xP_++yP_-$；它
在 $\mathrm{Alt}$ 上为零，所以 $y=0$，而它的值域含于 $\mathrm{Alt}\cap\mathrm{Sym}=0$，所以 $T=0$。$\square$

**定理 5.10（二阶 Weingarten 公式）。** 设 $k\ge2$，$U$ 在 $U(k)$ 上服从 Haar 分布，
$X$ 是 $\mathbb C^k\otimes\mathbb C^k$ 上的矩阵。则
$$
\mathbb E\,[(U\otimes U)X(U\otimes U)^*]=\alpha\,1+\beta F,\qquad
\alpha=\frac{\operatorname{Tr}X-\operatorname{Tr}(FX)/k}{k^2-1},\quad
\beta=\frac{\operatorname{Tr}(FX)-\operatorname{Tr}X/k}{k^2-1}.
$$
等价地，在群 $S_2=\{e,\tau\}$ 上令 $\mathrm{Wg}(e)=\frac1{k^2-1}$，$\mathrm{Wg}(\tau)=-\frac1{k(k^2-1)}$，则
$$
\mathbb E\big[U_{i_1j_1}U_{i_2j_2}\overline{U_{i'_1j'_1}U_{i'_2j'_2}}\big]
=\sum_{\sigma,\pi\in S_2}\delta_{i_1i'_{\sigma(1)}}\delta_{i_2i'_{\sigma(2)}}\,
\delta_{j_1j'_{\pi(1)}}\delta_{j_2j'_{\pi(2)}}\,\mathrm{Wg}(\sigma\pi^{-1}).
$$

*证明。* 令 $T=\mathbb E[(U\otimes U)X(U\otimes U)^*]$。与命题 5.7 一样，$T$ 与
所有 $V\otimes V$ 可交换，所以由引理 5.8，$T=\alpha1+\beta F$。取迹得
$\operatorname{Tr}T=\operatorname{Tr}X$；又因为 $F$ 与 $U\otimes U$ 可交换，
$\operatorname{Tr}(FT)=\mathbb E\operatorname{Tr}((U\otimes U)^*F(U\otimes U)X)=\operatorname{Tr}(FX)$。
这给出 $k^2\alpha+k\beta=\operatorname{Tr}X$ 和 $k\alpha+k^2\beta=\operatorname{Tr}(FX)$，这是一个
行列式为 $k^2(k^2-1)\neq0$ 的线性方程组，其解如上所述。对于
矩阵元，取 $X=|e_{j_1}\otimes e_{j_2}\rangle\langle e_{j'_1}\otimes e_{j'_2}|$：则
$\langle e_{i_1}\otimes e_{i_2},T\,e_{i'_1}\otimes e_{i'_2}\rangle$ 就是左边，$\operatorname{Tr}X=\delta_{j_1j'_1}\delta_{j_2j'_2}$
（$\pi=e$ 的项），$\operatorname{Tr}(FX)=\delta_{j_1j'_2}\delta_{j_2j'_1}$（$\pi=\tau$ 的项），而
$1$、$F$ 分别贡献 $\delta_{i_1i'_1}\delta_{i_2i'_2}$（$\sigma=e$ 的项）和
$\delta_{i_1i'_2}\delta_{i_2i'_1}$（$\sigma=\tau$）。合并各项得
$\alpha=\sum_\pi\mathrm{Wg}(e\pi^{-1})[\pi]$，$\beta=\sum_\pi\mathrm{Wg}(\tau\pi^{-1})[\pi]$，
其中 $[\pi]$ 表示 $\pi$ 对应的 $j$-δ 乘积。$\square$

矩阵 $\big(\mathrm{Wg}(\sigma^{-1}\pi)\big)_{\sigma,\pi}$ 是 $1$ 与 $F$ 关于
$\langle A,B\rangle=\operatorname{Tr}(A^*B)$ 的 Gram 矩阵
$\begin{pmatrix}k^2&k\\k&k^2\end{pmatrix}$ 的逆；这正是可以推广的形式。

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

最后一个矩是负的：Haar 酉矩阵的元素不是相互独立的，其
相关性的相对量级为 $1/k$。对于代数，不同块的元素是
相互独立的（推论 5.4），所以跨扇区的混合矩分解为单个块的矩的乘积。

对任意阶的矩，只要知道
$\{V^{\otimes n}\}$ 的交换子，同样的推理就适用；它由 $n$ 个张量因子的置换张成
（Schur–Weyl 对偶，见 R. Goodman, N. R. Wallach, *Symmetry, Representations, and Invariants*,
Springer, 2009, Chapter 9）。我们不加证明地叙述这一结果。

**定理 5.12（Weingarten 公式；不加证明）。** 设 $n\ge1$，$U$ 在 $U(k)$ 上服从 Haar 分布。则
$$
\mathbb E\big[U_{i_1j_1}\cdots U_{i_nj_n}\overline{U_{i'_1j'_1}\cdots U_{i'_nj'_n}}\big]
=\sum_{\sigma,\pi\in S_n}\prod_{m=1}^n\delta_{i_mi'_{\sigma(m)}}\delta_{j_mj'_{\pi(m)}}\;\mathrm{Wg}(\sigma\pi^{-1},k),
$$
其中 $\mathrm{Wg}(\cdot,k)$ 是 $S_n$ 上只依赖于循环型的函数。当 $n\le k$ 时，
矩阵 $\big(\mathrm{Wg}(\sigma^{-1}\pi,k)\big)_{\sigma,\pi}$ 是 Gram
矩阵 $\big(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)}\big)_{\sigma,\pi}$ 的逆；当 $n>k$ 时该矩阵是奇异的，
此时取其伪逆。$U$ 与 $\overline U$ 个数不同的矩为零。（D. Weingarten,
J. Math. Phys. 19 (1978), 999–1001; B. Collins, Int. Math. Res. Not. 2003, 953–982; B. Collins,
P. Śniady, Comm. Math. Phys. 264 (2006), 773–795.）

我们将命题 5.7 与 `channels.center_expectation` 对照检验，用精确的线性代数检验引理 5.8 中的
交换子，并用上面样本中的 $U(3)$ 块通过蒙特卡罗方法检验定理 5.10 和推论 5.11：

```python
X = alg.from_blocks([torch.randn(k, k, dtype=dt) for k in alg.k_factors])      # a fixed element
T = (U @ X @ U.adjoint()).matrix                                    # B samples of U X U*
mean, err = T.mean(0), T.std(0) / math.sqrt(B)
EX = channels.center_expectation(alg)(X).matrix[0]
random = err > 1e-12                                                 # entries that are not constant
z = ((mean - EX).abs()[random] / err[random]).max().item()
print(f"twirl over U(M) vs E(X): largest deviation {z:.2f} standard errors")
assert ((mean - EX).abs() <= 5 * err + 1e-12).all()                 # Proposition 5.7

def commutant_dim(ops, tol=1e-9):
    """Dimension of {X : XA = AX for all A in ops}, as the null space of X -> (AX - XA)_A."""
    d = ops[0].shape[0]
    eye = torch.eye(d, dtype=dt)
    ops = [A.contiguous() for A in ops]                              # torch.kron needs contiguous input
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    return int((torch.linalg.svdvals(L) < tol).sum())

k = 3
F = torch.zeros(k * k, k * k, dtype=dt)
for i in range(k):
    for j in range(k):
        F[j * k + i, i * k + j] = 1                                  # F (e_i (x) e_j) = e_j (x) e_i
VV = [torch.kron(W, W) for W in alg.random_unitary(k, batch_size=3).contiguous()]
assert all(close(F @ W, W @ F) for W in VV)
assert commutant_dim(VV) == 2                                        # Lemma 5.8: span{1, F}

u3 = u[2]                                                            # Haar on U(3)
mc(u3[:, 0, 0].abs() ** 4, 2 / (k * (k + 1)), "E|U_11|^4")
mc(u3[:, 0, 0].abs() ** 2 * u3[:, 0, 1].abs() ** 2, 1 / (k * (k + 1)), "E|U_11|^2 |U_12|^2")
mc((u3[:, 0, 0] * u3[:, 1, 1] * (u3[:, 0, 1] * u3[:, 1, 0]).conj()).real,
   -1 / ((k - 1) * k * (k + 1)), "E U_11 U_22 conj(U_12 U_21)")

n_s = 20_000                                                         # Theorem 5.10 for a random X
X2 = torch.randn(k * k, k * k, dtype=dt)
trX, trFX = X2.diagonal().sum(), (F @ X2).diagonal().sum()
alpha, beta = (trX - trFX / k) / (k ** 2 - 1), (trFX - trX / k) / (k ** 2 - 1)
UU = torch.einsum('bij,bkl->bikjl', u3[:n_s], u3[:n_s]).reshape(n_s, k * k, k * k)   # U (x) U
T2 = UU @ X2 @ UU.mH
mean, err = T2.mean(0), T2.std(0) / math.sqrt(n_s)
assert ((mean - alpha * torch.eye(k * k, dtype=dt) - beta * F).abs() <= 5 * err).all()
```

## 5.3 能级统计

设 $H$ 是一个随机自伴 $n\times n$ 矩阵，其特征值为
$\lambda_1\le\dots\le\lambda_n$。它的**能级间距**为 $s_i=\lambda_{i+1}-\lambda_i$。有两个标准的
系综：**GOE**，$H=(G+G^T)/2$，其中 $G$ 是元素相互独立且服从 $N(0,1)$ 的实矩阵；以及
**GUE**，$H=(G+G^*)/\sqrt2$，其中 $G$ 是 Ginibre 矩阵（1.7 节）。在这两种情形中，对角元都服从
$N(0,1)$，对角线以上元素的实部与虚部服从 $N(0,\tfrac12)$，而
$H$ 的密度正比于 $e^{-\operatorname{Tr}H^2/2}$。

**命题 5.13（Wigner 猜测）。** 设 $n=2$，$s=\lambda_2-\lambda_1$，$S=s/\mathbb Es$。
对 GOE（$\beta=1$）与 GUE（$\beta=2$），$S$ 的密度分别为
$$
p_1(S)=\frac\pi2\,S\,e^{-\pi S^2/4},\qquad p_2(S)=\frac{32}{\pi^2}\,S^2\,e^{-4S^2/\pi},
$$
分布函数分别为 $P(S\le t)=1-e^{-\pi t^2/4}$ 和
$P(S\le t)=\operatorname{erf}(2t/\sqrt\pi)-\frac{4t}\pi e^{-4t^2/\pi}$。特别地，当 $\varepsilon\to0$ 时，
$P(S\le\varepsilon)$ 的量级为 $\varepsilon^{\beta+1}$。

*证明。* 对 $H=\begin{pmatrix}a&b\\\bar b&c\end{pmatrix}$，特征值为
$\frac{a+c}2\pm\frac12\sqrt{(a-c)^2+4|b|^2}$，所以 $s=\sqrt{(a-c)^2+4|b|^2}$。*GOE：* $x=a-c$ 与
$y=2b$ 是相互独立的 $N(0,2)$ 变量，在极坐标下
$P(s\le t)=\int_{x^2+y^2\le t^2}\frac{e^{-(x^2+y^2)/4}}{4\pi}\,dx\,dy=\int_0^t\frac r2e^{-r^2/4}dr=1-e^{-t^2/4}$。
于是 $\mathbb Es=\int_0^\infty P(s>t)\,dt=\int_0^\infty e^{-t^2/4}dt=\sqrt\pi$，且
$P(S\le t)=P(s\le\sqrt\pi t)=1-e^{-\pi t^2/4}$；再求导即可。*GUE：* $x=a-c$，$y=2\operatorname{Re}b$，
$y'=2\operatorname{Im}b$ 是相互独立的 $N(0,2)$ 变量，在球坐标下
$s=\sqrt{x^2+y^2+y'^2}$ 的密度为 $4\pi r^2(4\pi)^{-3/2}e^{-r^2/4}=\frac{r^2e^{-r^2/4}}{2\sqrt\pi}$。其均值为
$\frac1{2\sqrt\pi}\int_0^\infty r^3e^{-r^2/4}dr=\frac{8}{2\sqrt\pi}=\frac4{\sqrt\pi}$。代入
$r=4S/\sqrt\pi$ 得到 $p_2$，再对 $p_2$ 分部积分
（$\int x^2e^{-ax^2}dx=-\frac{x}{2a}e^{-ax^2}+\frac1{2a}\int e^{-ax^2}dx$，其中 $a=4/\pi$）即得
分布函数。小 $\varepsilon$ 时的行为由 $p_\beta(S)\sim cS^\beta$ 得出。$\square$

指数是在数维数：实对称 $2\times2$ 矩阵简并当且仅当两个实
条件成立（$a=c$，$b=0$），复 Hermite 矩阵则当且仅当三个条件成立。对每个 $n$，特征值的联合
密度正比于 $\prod_{i<j}|\lambda_i-\lambda_j|^\beta e^{-\sum_i\lambda_i^2/2}$，并且
当 $n\to\infty$ 时，谱主体中的间距分布收敛到一个极限（Gaudin–Mehta），它与
Wigner 猜测接近但不相等（M. L. Mehta, *Random Matrices*, 3rd ed., Elsevier, 2004；不加
证明地引用）。

> **物理意义。** Wigner 提出，复杂量子系统（重原子核）的能级在统计上
> 类似于具有恰当对称性的随机哈密顿量的特征值：
> 若存在时间反演对称性则为实对称矩阵（$\beta=1$），否则为复 Hermite 矩阵
> （$\beta=2$）。于是能级相互*排斥*：同一对称性的两个能级几乎从不靠近
> （von Neumann 与 Wigner，1929）。能级排斥是量子混沌的典型特征，而独立
> 能级（无排斥）则是可积系统或局域化系统的特征（Bohigas–Giannoni–Schmit 猜想与
> Berry–Tabor 猜想）。

间距依赖于局部能级密度，而后者沿谱变化；要将间距与 $p_\beta$ 比较，需要进行
**展开**（unfolding），即借助计数函数的光滑近似 $\bar N$ 作变量替换 $\lambda\mapsto\bar N(\lambda)$，
使平均间距为 $1$。$\bar N$ 的选取是偏差的一个来源。相邻间距之比则避免了这一问题。

**定义。** **间距比**为 $r_i=\min(s_i,s_{i+1})/\max(s_i,s_{i+1})\in[0,1]$
（V. Oganesyan, D. A. Huse, Phys. Rev. B 75 (2007), 155111）。它在
$\lambda\mapsto a\lambda+b$（$a>0$）下不变，所以在密度在两个间距的范围内大致不变的地方，
它不需要展开。

**命题 5.14（独立能级）。** 设 $s,s'$ 相互独立且服从同一指数
分布（如同 Poisson 过程的两个相邻间距）。则 $r=\min(s,s')/\max(s,s')$ 在 $[0,1]$ 上满足
$P(r\le t)=\frac{2t}{1+t}$，密度为 $\frac2{(1+t)^2}$，并且
$$
\mathbb Er=2\ln2-1\approx0.3863 .
$$

*证明。* 由尺度变换，可设速率为 $1$。对 $x=s/s'$，$P(x\le t)=\mathbb E\,P(s\le ts'\mid s')=\mathbb E(1-e^{-ts'})=1-\frac1{1+t}=\frac t{1+t}$。
对 $t<1$，事件 $\{x\le t\}$ 与 $\{1/x\le t\}$ 不相交，$x$ 与 $1/x$ 同分布
（$s,s'$ 可交换），且 $\{r\le t\}$ 是它们的并，所以 $P(r\le t)=\frac{2t}{1+t}$。于是
$\mathbb Er=\int_0^1\frac{2t}{(1+t)^2}dt=2\Big[\ln(1+t)+\frac1{1+t}\Big]_0^1=2\ln2-1$。$\square$

对区间上 $n$ 个独立均匀分布的能级，这一结论同样精确成立，因为它们的间距
与 $(e_1,\dots,e_{n+1})/\sum_je_j$ 同分布，其中 $e_j$ 是独立的指数变量，因此每个
比值 $s_i/s_{i+1}=e_i/e_{i+1}$（W. Feller, *An Introduction to Probability Theory and Its
Applications*, Vol. II, Wiley, 1971, Chapter I；不加证明地引用）。对 GOE 与 GUE，在
$n\to\infty$ 的极限下，$\langle r\rangle=0.5307$ 和 $0.5996$（数值结果；Wigner 猜测的 $3\times3$ 类比
给出 $4-2\sqrt3\approx0.5359$ 和 $2\sqrt3/\pi-\frac12\approx0.6027$；Y. Y. Atas,
E. Bogomolny, O. Giraud, G. Roux, Phys. Rev. Lett. 110 (2013), 084101；不加证明地引用）。

函数 `SpinChain.level_spacing_ratio`（位于 `chains` 中）返回每一行能级的平均间距比。
我们检验命题 5.13 和 5.14 以及 GOE 的数值：

```python
ratio = SpinChain.level_spacing_ratio
n2 = 200_000
alg_r = TypeIAlgebra([2], [2], complex_valued=False, precision='double')     # real symmetric 2 x 2
alg_c = TypeIAlgebra([2], [2], precision='double')
G = torch.randn(n2, 1, 2, 2, dtype=torch.float64)
Gc = torch.randn(n2, 1, 2, 2, dtype=dt)
goe2 = alg_r.operator((G + G.mT) / 2, is_self_adjoint=True)
gue2 = alg_c.operator((Gc + Gc.mH) / math.sqrt(2), is_self_adjoint=True)
cdf = {1: lambda t: 1 - math.exp(-math.pi * t ** 2 / 4),
       2: lambda t: math.erf(2 * t / math.sqrt(math.pi)) - 4 * t / math.pi * math.exp(-4 * t ** 2 / math.pi)}
for beta, H, mean_s in ((1, goe2, math.sqrt(math.pi)), (2, gue2, 4 / math.sqrt(math.pi))):
    lam = H.eigenvalues()[0]                                         # (n2, 2), ascending
    S = (lam[:, 1] - lam[:, 0]) / mean_s
    mc(S, 1.0, f"beta={beta}: E S")
    for t in (0.25, 1.0, 2.0):
        mc((S <= t).double(), cdf[beta](t), f"beta={beta}: P(S <= {t})")

levels = torch.rand(4000, 50, dtype=torch.float64)                   # independent levels
mc(ratio(levels), 2 * math.log(2) - 1, "<r>, independent uniform levels")

n = 64
alg64 = TypeIAlgebra([n], [n], complex_valued=False, precision='double')
G = torch.randn(800, 1, n, n, dtype=torch.float64)
lam = alg64.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()[0]
mc(ratio(lam[:, n // 4: 3 * n // 4]), 0.5307, "<r>, GOE(64), middle half of the spectrum")
```

## 5.4 为何不能混合扇区

设 $H=\bigoplus_cH_c\in\mathcal M$ 是一个具有**相互独立的块**的随机自伴元素，
例如 GOE 块。作为 $\mathbb C^D$ 上的算子，它的谱是各块谱的
并。在一个块内部能级相互排斥；我们将证明不同块的能级之间并不排斥。

假设每个块都有有界连续的**能级密度** $\rho_c$：
对有界 Borel 函数 $\varphi$ 有 $\mathbb E\sum_i\varphi(\lambda_i(H_c))=\int\varphi\,\rho_c$（从而
$\int\rho_c=k_c$）。对 GOE 和 GUE，这由上面引用的联合密度得出。

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

*证明。* 排好序的特征值是矩阵的连续函数，所以
$\Lambda_c=(\lambda_i(H_c))_i$ 与 $\Lambda_d$ 是相互独立的随机向量。由 Fubini 定理（对
$\Lambda_d$ 取条件），$\mathbb E\sum_{i,j}\psi(\Lambda_{c,i},\Lambda_{d,j})=\mathbb E\sum_j g(\Lambda_{d,j})$，其中
$g(y)=\mathbb E\sum_i\psi(\Lambda_{c,i},y)=\int\psi(x,y)\rho_c(x)\,dx$ 是有界 Borel 函数；对 $g$ 应用
$\rho_d$ 的定义即可。取 $\psi(x,y)=1\{|x-y|<\varepsilon\}$，
$$
\frac{\mathbb EN_\varepsilon}{2\varepsilon}=\sum_{c<d}\int\rho_c(x)\Big(\frac1{2\varepsilon}\int_{x-\varepsilon}^{x+\varepsilon}\rho_d(y)\,dy\Big)dx ;
$$
括号内的量由连续性趋于 $\rho_d(x)$，且以 $\sup\rho_d$ 为界，而 $\rho_c$ 可积，
所以由控制收敛定理得到该极限。最后 $N_\varepsilon\le\sum_{c<d}k_ck_d$，所以
$\mathbb EN_\varepsilon\le\big(\sum_{c<d}k_ck_d\big)P(N_\varepsilon\ge1)$，而 $N_\varepsilon\ge1$ 蕴含 $\delta<\varepsilon$。$\square$

因此，对于谱相互重叠的两个块，$P(\delta<\varepsilon)$ 至少是 $\varepsilon$ 的量级，
而在一个 $2\times2$ GOE 块内部它等于 $1-e^{-\varepsilon^2/4}\le\varepsilon^2/4$
（命题 5.13，未归一化的间距），对 $n\times n$ GOE 块则是
$\varepsilon^2$ 的量级（由联合密度得出；此处不证）。尽管每个块都是完美的 GOE，
合并谱却通不过能级排斥检验。

*多个扇区（一个推导，而非证明）。* 对平稳谱，在长度为 $s$ 的区间内
找不到任何能级的概率 $\mathcal G(s)$ 决定了间距密度：在
平均间距为 1 时 $p(s)=\mathcal G''(s)$（Mehta, *op. cit.*）。独立谱叠加时，这些空隙
概率相乘。若 $C$ 个扇区各自贡献能级总数的 $1/C$，且在各自展开后的单位下都有空隙
概率 $\mathcal G_1$，则 $\mathcal G(s)=\mathcal G_1(s/C)^C$，其中
$\mathcal G_1(x)=1-x+\int_0^x(x-t)p_1(t)\,dt$。在 $0$ 处求两次导数，
$$
p(0)=1-\frac1C+\frac{p_1(0)}C=1-\frac1C ,
$$
这是定理 5.15 对展开后的谱的类比；并且 $\mathcal G_1(s/C)^C\to e^{-s}$，即 Poisson 过程的
空隙概率：当 $C\to\infty$ 时，合并谱看起来像独立能级
（N. Rosenzweig, C. E. Porter, Phys. Rev. 120 (1960), 1698–1714; M. V. Berry, M. Robnik, J. Phys. A 17
(1984), 2413–2421；Poisson 极限的严格版本是关于稀疏点过程叠加的极限定理，
见 B. Grigelionis, Theory Probab. Appl. 8 (1963)）。

在本库中，各块从一开始就是分开的：`eigenvalues()` 对每个扇区返回一个谱，
混合它们需要显式地调用 `torch.cat`。我们把 $C$ 个大小为 $64$ 的独立 GOE 块
放进一个代数，并比较每个扇区的平均间距比与合并谱的平均间距比：

```python
k, merged_r = 64, {}
for C in (1, 2, 4, 8, 16):
    alg_C = TypeIAlgebra([k] * C, [k] * C, complex_valued=False, precision='double')
    G = torch.randn(800 // C, C, k, k, dtype=torch.float64)
    lam = alg_C.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()    # one spectrum per sector
    per_sector = torch.stack([ratio(l[:, k // 4: 3 * k // 4]) for l in lam])
    merged = torch.cat(lam, dim=-1).sort(dim=-1).values                       # the spectrum on C^D
    r_m = ratio(merged[:, k * C // 4: 3 * k * C // 4])
    merged_r[C] = r_m.mean().item()
    print(f"C = {C:2d}: <r> per sector {per_sector.mean():.4f}, merged {merged_r[C]:.4f}"
          f" +- {r_m.std() / math.sqrt(len(r_m)):.4f}")
    assert abs(per_sector.mean() - 0.5307) < 0.01                    # every sector is GOE
assert merged_r[1] > merged_r[2] > merged_r[4] > merged_r[16]
assert abs(merged_r[16] - (2 * math.log(2) - 1)) < 0.01              # close to independent levels
```

每个扇区都有 $\langle r\rangle\approx0.53$；合并谱在两个扇区时降到 $0.42$，在 $8$ 和 $16$ 个
扇区时降到 $0.39$，即独立能级的数值。

在物理模型中也会发生同样的事情。带随机场的 Heisenberg 链（第 9 章）守恒
磁化强度，所以它的哈密顿量属于 $\bigoplus_{N=0}^LM_{\binom LN}$；对弱无序，每个扇区内部的
能级具有 GOE 统计（数值结果；A. Pal, D. A. Huse, Phys. Rev. B 82 (2010),
174411）。它的各块*不是*相互独立的——它们都依赖于同一组随机场——所以
定理 5.15 不能照原样应用；实验显示出同样的效应。

```python
chain = SpinChain(10, complex_valued=False)                         # sectors N = 0, ..., 10
lam = chain.random_field_heisenberg(W=0.5, batch_size=20).eigenvalues()
k5 = chain.sector_dim(5)                                            # 252
r_sector = ratio(lam[5][:, k5 // 4: 3 * k5 // 4]).mean().item()
spectrum = torch.cat(lam, dim=-1).sort(dim=-1).values               # all 1024 levels mixed
r_all = ratio(spectrum[:, 256:768]).mean().item()
print(f"Heisenberg chain, L = 10: <r> in sector N = 5: {r_sector:.4f}, all sectors mixed: {r_all:.4f}")
assert r_sector > 0.49 and r_all < 0.42
```

在扇区 $N=5$ 内部，间距比为 $0.52$（一个有限尺寸值，略低于 $0.5307$）；
所有扇区混合后的谱给出 $0.39$。

> **物理意义。** 具有不同守恒量子数的能级之间没有相互作用，所以它们可以
> 自由交叉；只有同一对称性的能级才相互排斥。因此，混合了对称扇区的谱
> 即使在每个扇区内的动力学完全混沌时，也会呈现类 Poisson 统计，
> 而一个未被分辨的对称性可能被误认为可积性或局域化。业内的经验
> 规则——在计算能级统计之前先分辨所有对称性——说的正是：应当
> 在代数 $\bigoplus_cM_{k_c}$ 中计算，而不是在 $M_D$ 中计算。

## 5.5 可与之比较的工具

Haar 酉矩阵：`scipy.stats.unitary_group`（即定理 5.3 中经过修正的 QR）和 `ortho_group`，
QuTiP 的 Haar 随机酉矩阵，以及任何 `numpy`/`torch` 的 QR 再加上相位修正。
Weingarten 演算：符号计算软件包 IntU（Mathematica；Z. Puchała, J. A. Miszczak）、RTNI
（Mathematica 与 Python；M. Fukuda, R. König, I. Nechita, J. Phys. A 52 (2019), 425303）以及
Python 软件包 Haarpy 能符号地计算任意阶的矩；本章证明的公式是它们的
二阶情形。能级统计通常用几行 `numpy` 代码计算，
而 QuSpin 等精确对角化程序逐个扇区地提供谱。本库
增加的是把代数的酉群作为一个对象：在批量调用中抽取所有大小的独立 Haar 块，
把 twirl 等同于条件期望
`channels.center_expectation`，并按扇区返回谱，使得混合扇区永远不是
默认行为。

## 习题

1. 设 $U$ 在 $U(\mathcal M)$ 上服从 Haar 分布，$\operatorname{Tr}U=\sum_c\operatorname{Tr}U_c$。证明
   $\mathbb E|\operatorname{Tr}U|^2=C$ 且 $\mathbb E|\operatorname{Tr}U|^4=\sum_cm(k_c)+2C(C-1)$，其中
   $m(1)=1$，当 $k\ge2$ 时 $m(k)=2$；并与 $U(D)$（$D\ge2$）上 Haar 酉矩阵的 $1$ 和 $2$ 比较。对
   $\mathcal M=M_1\oplus M_2\oplus M_3$ 加以检验。
2. 利用定理 5.12，证明当 $n\le k$ 时 $\mathbb E|U_{11}|^{2n}=\frac{n!}{k(k+1)\cdots(k+n-1)}$。
   对 $k=3$ 数值计算 $3$ 阶 Weingarten 矩阵，并检验 $\mathbb E|U_{11}|^6=\frac1{10}$。
3. **圆正交系综**（`measure='coe'`）是 $S=W^TW$，其中 $W$ 在 $U(k)$ 上服从 Haar 分布。证明
   $S$ 是酉且对称的，它的分布对每个
   $V\in U(k)$ 在 $S\mapsto V^TSV$ 下不变，并且其本征相位的平均密度是均匀的。数值检验它的
   本征相位满足 $\langle r\rangle\approx0.53$，与 GOE 相同。

<details>
<summary>解答</summary>

1. 记 $t_c=\operatorname{Tr}U_c$；各 $t_c$ 相互独立（推论 5.4）。把 $U_c$ 乘以
   $e^{i\theta}$ 不改变其分布，并把 $t_c$ 乘以 $e^{i\theta}$；所以矩
   $\mathbb E[t_at_b\bar t_c\bar t_d]$ 为零，除非每个扇区在 $(a,b)$ 中出现的次数与在
   $(c,d)$ 中相同，即作为多重集 $\{a,b\}=\{c,d\}$；同理，当 $a\neq b$ 时 $\mathbb E[t_a\bar t_b]=0$。由
   命题 5.7，$\mathbb E[(U_c)_{ii}\overline{(U_c)_{jj}}]=\langle e_i,\mathbb E[U_ce_ie_j^*U_c^*]e_j\rangle=\delta_{ij}/k_c$，
   所以 $\mathbb E|t_c|^2=1$，$\mathbb E|\operatorname{Tr}U|^2=C$。在 $\mathbb E|\sum_at_a|^4$ 中，
   $a=b$ 的项给出 $\sum_a\mathbb E|t_a|^4$，$a\neq b$ 的项给出
   $2\sum_{a\neq b}\mathbb E|t_a|^2\,\mathbb E|t_b|^2=2C(C-1)$（$(c,d)$ 有两种选法）。对 $k=1$，$|t|=1$。
   对 $k\ge2$，$\mathbb E|t|^4=\sum_{i,j,i',j'}\mathbb E[U_{ii}U_{jj}\overline{U_{i'i'}U_{j'j'}}]$；在
   定理 5.10 中，$\sigma=\pi$ 的两项以同一种方式迫使 $\{i',j'\}=\{i,j\}$，求和后各给出
   $k^2\mathrm{Wg}(e)$；$\sigma\neq\pi$ 的两项迫使 $i=j=i'=j'$，各给出
   $k\,\mathrm{Wg}(\tau)$：总计 $2\big(\frac{k^2}{k^2-1}-\frac1{k^2-1}\big)=2$。$C=1$，
   $k=D$ 的情形对 $U(D)$ 给出 $1$ 和 $2$。对 $M_1\oplus M_2\oplus M_3$：
   $\mathbb E|\operatorname{Tr}U|^4=1+2+2+12=17$，而 $U(6)$ 上的 Haar 酉矩阵为 $2$。

```python
tr = U.trace                                                        # Tr U = sum_c Tr U_c, shape (B,)
mc(tr.abs() ** 2, 3.0, "E|Tr U|^2")
mc(tr.abs() ** 4, 17.0, "E|Tr U|^4")
```

2. 当所有指标都等于 $1$ 时，每个 δ 都是 $1$，所以
   $\mathbb E|U_{11}|^{2n}=\sum_{\sigma,\pi}\mathrm{Wg}(\sigma\pi^{-1})=\mathbf 1^TW\mathbf 1$，其中
   $W=\Gamma^{-1}=(\mathrm{Wg}(\sigma^{-1}\pi))$，$\Gamma=(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)})$。这里用到
   $\mathrm{Wg}$ 只依赖于循环型，以及 $\sigma\pi^{-1}=\sigma(\pi^{-1}\sigma)\sigma^{-1}$
   共轭于 $\pi^{-1}\sigma=(\sigma^{-1}\pi)^{-1}$，后者与 $\sigma^{-1}\pi$ 具有相同的循环型。$\Gamma$ 的每个
   行和都是 $c_n(k)=\sum_{\rho\in S_n}k^{\#\mathrm{cycles}(\rho)}$，所以 $\Gamma\mathbf 1=c_n(k)\mathbf 1$，
   $W\mathbf 1=\mathbf 1/c_n(k)$，$\mathbf 1^TW\mathbf 1=n!/c_n(k)$。最后，由归纳法
   $c_n(x)=x(x+1)\cdots(x+n-1)$：$\{1,\dots,n\}$ 的一个置换可由
   $\{1,\dots,n-1\}$ 的一个置换得到，或者把 $n$ 作为不动点加入（多一个循环，因子 $x$），或者把 $n$
   插入到某个元素在其循环中的后面（$n-1$ 种方式，循环数不变）；这是一个双射，所以
   $c_n(x)=(x+n-1)c_{n-1}(x)$。

```python
import itertools

def n_cycles(p):                                                    # number of cycles of a permutation
    seen, count = set(), 0
    for i in range(len(p)):
        count += i not in seen
        while i not in seen:
            seen.add(i)
            i = p[i]
    return count

perms = list(itertools.permutations(range(3)))                      # S_3; s.index(j) = s^{-1}(j)
Gram = torch.tensor([[3.0 ** n_cycles([s.index(p[i]) for i in range(3)]) for p in perms] for s in perms],
                    dtype=torch.float64)                            # k^{#cycles(s^{-1} p)}, k = 3
Wg3 = torch.linalg.inv(Gram)                                        # the Weingarten matrix, n = k = 3
assert abs(Wg3.sum().item() - 6 / (3 * 4 * 5)) < 1e-12
mc(u3[:, 0, 0].abs() ** 6, 0.1, "E|U_11|^6")
```

3. $S^*S=W^*\bar WW^TW=1$，因为 $W^T$ 是酉的；并且 $S^T=W^TW=S$。对 $V\in U(k)$，$WV$ 与 $W$
   同分布（右不变性，引理 5.1），而 $(WV)^T(WV)=V^TSV$。取 $V=e^{i\theta/2}1$，$e^{i\theta}S$
   与 $S$ 同分布，所以本征相位集合的分布在旋转下不变，平均密度
   也是如此，因此它是均匀的。特别地，间距不需要展开。COE 的局部
   统计在 $k$ 很大时与 GOE 的一致（F. J. Dyson, J. Math. Phys. 3
   (1962), 140–156; Mehta, *op. cit.*；不加证明地引用）。

```python
alg_s = TypeIAlgebra([32], [32], precision='double')
S = alg_s.random_unitary_operator(400, measure='coe').matrix[:, 0]
assert close(S, S.mT) and close(S @ S.mH, torch.eye(32, dtype=dt))
phases = torch.linalg.eigvals(S).angle().sort(dim=-1).values
s = torch.diff(torch.cat([phases, phases[:, :1] + 2 * math.pi], dim=-1), dim=-1)   # around the circle
r = torch.minimum(s, s.roll(-1, -1)) / torch.maximum(s, s.roll(-1, -1))
mc(r.mean(-1), 0.5307, "<r>, COE(32) eigenphases")
```

</details>

**下一章：** [第 6 章](06_quantum_channels.md)讨论量子信道：Kraus 表示与 Stinespring
表示、关于迹的对偶，以及扇区之间的映射。本章中作为 Haar 平均出现的条件
期望 $E$ 是第一个例子。
