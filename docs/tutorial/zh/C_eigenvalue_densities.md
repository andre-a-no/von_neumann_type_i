[English](../C_eigenvalue_densities.md) · [Español](../es/C_eigenvalue_densities.md) · [Français](../fr/C_eigenvalue_densities.md) · [Deutsch](../de/C_eigenvalue_densities.md) · **中文** · [日本語](../ja/C_eigenvalue_densities.md) · [Русский](../ru/C_eigenvalue_densities.md)

# 附录 C. 特征值密度与 Page 公式

第 5 章用到高斯系综特征值的联合密度，第 7 章用到随机纯态平均纠缠熵的 Page 公式。本附录证明这两者。
路线是：Weyl 积分公式（"对角化"的 Jacobi 行列式）、高斯系综、借助正交多项式求行列式型密度的单点函数、
Wishart–Laguerre 系综与约化态上的诱导测度，最后是用 Laguerre 多项式进行的精确计算。

**你将学到**为什么不变随机矩阵的特征值按
$\prod_{i<j}|\lambda_i-\lambda_j|^\beta$ 相互排斥，为什么其特征向量此时服从 Haar 分布并且与特征值相互独立，
正交归一多项式构成的行列式如何被积掉，以及
$\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}$ 如何由少数几个积分得出。

除了本课程的预备知识（测度论；对导数可逆的单射 $C^1$ 映射成立的 $\mathbb R^N$ 中的变量替换公式；反函数定理与
隐函数定理）之外，我们还用到：第 5 章的引理 5.1（Haar 测度的唯一性与右不变性）、
引理 5.2 和定理 5.3（借助高斯 QR 分解得到 $U(n)$ 上 Haar 测度的存在性），附录 B 的
定理 B.17（$O(n)$ 的相应结论），附录 A 的命题 A.16(a)（复高斯向量的酉不变性）
和定理 A.20(b)（$|z|^2$ 的分布），以及第 7 章定理 7.11 的证明中的表示
$\rho_A=GG^*/\operatorname{Tr}GG^*$。

**记号。** $\beta\in\{1,2\}$，相应地 $\mathbb F=\mathbb R$ 或 $\mathbb C$。$\mathcal H$ 是
$n\times n$ 实对称矩阵（$\beta=1$）或复 Hermite 矩阵（$\beta=2$）构成的空间，其上的
Lebesgue 测度 $dH$ 取在坐标 $h_{ii}$、$\operatorname{Re}h_{ij}$ 以及（当 $\beta=2$ 时）
$\operatorname{Im}h_{ij}$（$i<j$）中；它的维数是 $n+\beta n(n-1)/2$。$K=O(n)$ 或 $U(n)$，$T$ 是
$K$ 的对角元构成的群（符号或相位），$\mathfrak k=\{X\in M_n(\mathbb F):X^*=-X\}$，
$\mathfrak t$ 是它的对角部分，$\mathfrak m$ 是它的非对角部分（对角线为零），坐标为
$\operatorname{Re}X_{ij}$、$\operatorname{Im}X_{ij}$（$i<j$）。对 $\lambda\in\mathbb R^n$，
$\Lambda=\operatorname{diag}(\lambda)$，$\Delta(\lambda)=\prod_{i<j}(\lambda_j-\lambda_i)$；$C=\{\lambda:\lambda_1<\dots<\lambda_n\}$，
在其上 $\Delta>0$。对 $H\in\mathcal H$，$\lambda(H)$ 是按递增顺序排列的特征值构成的向量。
$\|X\|_2^2=\operatorname{Tr}X^*X$。对 $V\in K$，映射 $H\mapsto VHV^*$ 是 $\mathcal H$ 上的线性映射，
并保持 $\|\cdot\|_2$；线性等距满足 $|\det|=1$，因此它**保持 $dH$**。

## C.1 单谱与特征标架

**引理 C.1.** (a) 具有 $n$ 个互不相同特征值的矩阵构成的集合 $\Omega\subseteq\mathcal H$ 是
开集，并且 $\mathcal H\setminus\Omega$ 是 Lebesgue 零测集。(b) 在 $\Omega$ 上映射 $H\mapsto\lambda(H)$ 是
$C^\infty$ 的，并且谱投影 $P_j(H)=\prod_{k\neq j}\frac{H-\lambda_k(H)}{\lambda_j(H)-\lambda_k(H)}$ 是连续的。

*证明。* (a) 令 $A_{ij}=\lambda_j^{i-1}$，其中 $\lambda=\lambda(H)$。则 $\det A=\Delta(\lambda)$（Vandermonde），并且
$(AA^T)_{ij}=\sum_k\lambda_k^{i+j-2}=\operatorname{Tr}H^{i+j-2}$，所以
$\operatorname{disc}(H)=\det\big(\operatorname{Tr}H^{i+j-2}\big)_{i,j=1}^n=\Delta(\lambda(H))^2$ 是 $H$ 的坐标的多项式，
它在 $\operatorname{diag}(1,\dots,n)$ 处不为零，并且 $\Omega=\{\operatorname{disc}\neq0\}$ 是开集。
$\mathbb R^N$ 上非零多项式 $p$ 的零点集是零测集，这对 $N$ 归纳可得：$N=1$ 时它是
有限集；否则写 $p=\sum_jp_j(x')x_N^j$，其中某个 $p_{j_0}\neq0$；对 $p_{j_0}$ 的零点集（由归纳假设为零测集）之外的 $x'$，
$x_N\mapsto p(x',x_N)$ 只有有限多个零点，对零点集的示性函数应用 Tonelli 定理
即得测度为 $0$。
(b) $\chi(H,t)=\det(t-H)=\prod_i(t-\lambda_i)$ 是 $(H,t)$ 的实多项式。在 $H_0\in\Omega$ 处，
$\partial_t\chi(H_0,\lambda_j(H_0))=\prod_{k\neq j}(\lambda_j-\lambda_k)\neq0$，所以隐函数定理给出 $H_0$ 附近的
$C^\infty$ 函数 $r_j$，满足 $\chi(H,r_j(H))=0$ 和 $r_j(H_0)=\lambda_j(H_0)$。由连续性，在 $H_0$ 附近
$r_1<\dots<r_n$：它们是一个 $n$ 次多项式的 $n$ 个不同的根，因此按顺序就是它的全部根，
从而 $r_j=\lambda_j$。若 $H=\sum_k\lambda_kQ_k$ 是谱分解，则 Lagrange 插值给出
$P_j(H)=Q_j$；连续性由 $\lambda$ 的连续性得出。$\square$

对 $H\in\Omega$，**特征标架**是满足 $H=V\Lambda V^*$（$\Lambda=\operatorname{diag}\lambda(H)$）的 $V\in K$。若 $V,V'$
都是特征标架，则 $M=V'^*V$ 与 $\Lambda$ 交换，所以 $(\lambda_i-\lambda_j)M_{ij}=0$ 且 $M\in T$；反之，对每个 $D\in T$，$VD$
都是特征标架。取
$P_j(H)e_m/\|P_j(H)e_m\|$ 作为第 $j$ 列，其中 $m$ 是使 $P_j(H)e_m\neq0$ 的最小指标，就得到一个 **Borel 特征标架** $F(H)$：它在这些指标固定的
各个 Borel 集上连续，其各列是正交归一的特征向量，并且当 $H$ 为实矩阵时它也是实的。

**引理 C.2（Haar 测度）。** $U(n)$、$O(n)$ 和 $T$ 都有 Haar 测度，并且对所有 $W\in K$ 和有界 Borel 函数 $f$，
$\int_Kf(VW)\,dV=\int_Kf\,dV$。

*证明。* $U(n)$：定理 5.3。$O(n)$：定理 B.17（实高斯 QR 构造：实高斯矩阵几乎必然
可逆，并且其 QR 分解的正交因子服从 Haar 分布）。$T$：$\{\pm1\}^n$ 上的均匀测度，或各圆周上
归一化弧长的乘积，在乘以固定元素下不变。右不变性即引理 5.1。$\square$

**引理 C.3（不变系综的特征标架）。** 设 $H$ 是 $\mathcal H$ 中的随机元素，
$P(H\in\Omega)=1$，并且其分布对每个 $V\in K$ 在 $H\mapsto VHV^*$ 下不变。设 $D$ 在 $T$ 上均匀分布且
与 $H$ 独立，$\tilde U=F(H)D$。则 $\tilde U$ 是 $H$ 的特征标架，它在 $K$ 上服从 Haar 分布，
并且与 $\lambda(H)$ 独立。

*证明。* 固定 $V\in K$。$F(VHV^*)$ 与 $VF(H)$ 都是 $VHV^*$ 的特征标架，所以
$F(VHV^*)=VF(H)E(H)$，其中 $E(H)=(VF(H))^*F(VHV^*)\in T$ 是 $H$ 的 Borel 函数。对 $\mathbb R^n$ 上的有界 Borel 函数
$\varphi$ 和 $K$ 上的有界 Borel 函数 $\psi$，由于 $VHV^*$ 与 $H$ 同分布且 $\lambda(VHV^*)=\lambda(H)$，
$$
\mathbb E\big[\varphi(\lambda(H))\psi(\tilde U)\big]=\mathbb E\big[\varphi(\lambda(H))\psi(VF(H)E(H)D)\big]
=\mathbb E\big[\varphi(\lambda(H))\psi(V\tilde U)\big],
$$
其中第二步先对 $D$ 积分（Fubini）：对固定的 $H$，$E(H)D$ 在 $T$ 上均匀分布。现在设 $V$
服从 Haar 分布且与其余一切独立；由 Fubini 定理和引理 C.2，
$\mathbb E[\varphi(\lambda(H))\psi(\tilde U)]=\mathbb E\big[\varphi(\lambda(H))\int_K\psi(V\tilde U)\,dV\big]
=\mathbb E[\varphi(\lambda(H))]\int_K\psi\,dV$。取 $\varphi=1$ 知 $\tilde U$ 服从 Haar 分布；乘积形式就是
独立性。$\square$

## C.2 Weyl 积分公式

对 $\varepsilon>0$，令 $\mathfrak m_\varepsilon=\{X\in\mathfrak m:\|X\|_2<\varepsilon\}$ 以及
$$
\Psi(X,\lambda)=e^X\Lambda e^{-X},\qquad (X,\lambda)\in\mathfrak m_\varepsilon\times C .
$$
对 $X\in\mathfrak k$，$(e^X)^*=e^{-X}=(e^X)^{-1}$，所以 $e^X\in K$。

**引理 C.4（一个坐标卡）。** 存在 $\varepsilon>0$ 和连续函数 $j:\overline{\mathfrak m_\varepsilon}\to(0,\infty)$，
$j(0)=1$，使得：(a) $|\det D\Psi(X,\lambda)|=\Delta(\lambda)^\beta j(X)$；(b) $\Psi$ 是单射；(c) 对
$H\in\Omega$，$H\in\Psi(\mathfrak m_\varepsilon\times C)$ 当且仅当 $F(H)\in S=\{e^XD:X\in\mathfrak m_\varepsilon,D\in T\}$，并且
$S$ 是满足 $ST=S$ 的 Borel 集。

*证明。* (a) $\Psi$ 是 $C^\infty$ 的（幂级数）。对 $Y\in\mathfrak m$，令 $A_X(Y)=e^{-X}\frac{d}{dt}\big|_0e^{X+tY}$。
对 $W(t)=e^{X+tY}\in K$ 微分 $W(t)^*W(t)=1$ 得 $A_X(Y)\in\mathfrak k$，再利用
$\frac{d}{dt}W^{-1}=-W^{-1}W'W^{-1}$，对 $g=e^X$ 有
$$
D\Psi(X,\lambda)(Y,\mu)=g\big([A_X(Y),\Lambda]+\operatorname{diag}(\mu)\big)g^{-1}.
$$
对角矩阵与 $\Lambda$ 交换，所以 $[A_X(Y),\Lambda]=[B_XY,\Lambda]$，其中 $B_X=\pi\circ A_X:\mathfrak m\to\mathfrak m$，
$\pi$ 是 $\mathfrak k$ 沿 $\mathfrak t$ 到 $\mathfrak m$ 上的投影。于是 $D\Psi=\operatorname{Ad}_g\circ L_\lambda\circ(B_X\times1)$，
其中 $L_\lambda(Z,\mu)=[Z,\Lambda]+\operatorname{diag}(\mu)$。由于 $([Z,\Lambda])_{ij}=Z_{ij}(\lambda_j-\lambda_i)$，在坐标中
$L_\lambda$ 把 $\operatorname{Re}Z_{ij}$（以及 $\operatorname{Im}Z_{ij}$）乘以 $\lambda_j-\lambda_i$，并把 $\mu$ 映到对角线上：
$|\det L_\lambda|=\Delta(\lambda)^\beta$。$\operatorname{Ad}_g$ 保持 $dH$。令 $j(X)=|\det B_X|$：它是连续的，并且
$B_0=1$，因为 $\frac d{dt}e^{tY}|_0=Y$；选取 $\varepsilon$ 使得在 $\overline{\mathfrak m_\varepsilon}$ 上 $j>0$。
(b) 若 $\Psi(X,\lambda)=\Psi(X',\lambda')=H$，则 $\lambda=\lambda'=\lambda(H)$，并且 $e^X,e^{X'}$ 是 $H\in\Omega$ 的特征标架，
所以 $e^X=e^{X'}D$，其中 $D\in T$。在 $\mathfrak m\times\mathfrak t$ 上令 $G(X,t)=e^Xe^t$（对 $\beta=1$，$\mathfrak t=\{0\}$）；$DG(0,0)(Y,s)=Y+s$ 是单射，
并且由于 $Y$ 是非对角的而 $s$ 是对角的，$\|Y+s\|_2^2=\|Y\|_2^2+\|s\|_2^2$，所以对乘积范数有 $\|DG(0,0)v\|=\|v\|$；令 $m=1$。选取 $\delta>0$ 使得当 $\|z\|<\delta$ 时 $\|DG(z)-DG(0)\|\le m/2$；对此球中的
$u,v$，由 $G(u)-G(v)=\int_0^1DG(v+s(u-v))(u-v)\,ds$ 得 $\|G(u)-G(v)\|\ge\frac m2\|u-v\|$，所以 $G$ 在那里是单射。
（$\delta$ 只依赖于 $n,\beta$；$\varepsilon$ 在它之后选取。）
由幂级数和
$\|AB\|_2\le\|A\|_2\|B\|_2$，有 $\|D-1\|_2=\|e^{X'}(D-1)\|_2=\|e^X-e^{X'}\|_2\le\|e^X-1\|_2+\|e^{X'}-1\|_2\le2(e^\varepsilon-1)$。对 $\beta=1$，$D=\operatorname{diag}(\pm1)$，若有一个元素为 $-1$ 就会有
$\|D-1\|_2\ge2$；当 $\varepsilon<\log2$ 时，$D=1=e^t$，$t=0$。对 $\beta=2$，$D=e^t$，其中 $t=i\operatorname{diag}(\theta)$，
$\theta_j\in(-\pi,\pi]$，并且 $|\theta_j|\le\frac\pi2|e^{i\theta_j}-1|$（因为 $|e^{i\theta}-1|=2|\sin\frac\theta2|$，且在 $[0,\pi/2]$ 上 $\sin x\ge2x/\pi$），所以
$\|t\|_2\le\frac\pi2\|D-1\|_2\le\pi(e^\varepsilon-1)$。当 $\varepsilon$ 足够小（$\varepsilon^2+\pi^2(e^\varepsilon-1)^2<\delta^2$）时，$(X,0)$ 与 $(X',t)$ 位于 $\delta$-球中，并且
$G(X,0)=e^X=e^{X'}D=G(X',t)$，所以 $X=X'$。
(c) 若 $H=\Psi(X,\lambda')$，则 $e^X$ 是特征标架，所以 $F(H)\in e^XT\subseteq S$。若 $F(H)=e^XD$，则
$H=F(H)\Lambda F(H)^*=\Psi(X,\lambda(H))$。$S=G(\mathfrak m_\varepsilon\times\mathfrak t)$（或 $e^{\mathfrak m_\varepsilon}T$，$T$ 有限）是
$\sigma$-紧集的连续像，因而是 $\sigma$-紧的，从而是 Borel 集；$ST=S$ 是显然的。$\square$

**定理 C.5（Weyl 积分公式）。** 存在常数 $c_{n,\beta}>0$，使得对每个 Borel 函数
$f:\mathcal H\to[0,\infty]$，
$$
\int_{\mathcal H}f(H)\,dH=c_{n,\beta}\int_C\Delta(\lambda)^\beta\Big(\int_Kf(V\Lambda V^*)\,dV\Big)d\lambda .
$$
等价地，$\Phi(V,\lambda)=V\Lambda V^*$ 把 $K\times C$ 上的测度 $c_{n,\beta}\,dV\otimes\Delta(\lambda)^\beta d\lambda$ 映为
$\mathcal H$ 上的 Lebesgue 测度。因子 $\Delta(\lambda)^\beta$ 是坐标卡 $\Psi$ 在 $X=0$ 处的 Jacobi 行列式
（引理 C.4(a)）；常数 $c_{n,\beta}$ 反映 Haar 测度的归一化，下文从不需要它的显式
值，因为由该公式得到的每个密度都用其总质量归一化。

*证明。* *参考系综。* 设 $H_0$ 具有密度 $p_0(H)=Z_0^{-1}e^{-\operatorname{Tr}H^2/2}$；由于
$\operatorname{Tr}H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2$，$Z_0$ 是有限个高斯积分的乘积。$p_0$ 和 $dH$
在 $H\mapsto VHV^*$ 下不变，因此 $H_0$ 的分布也不变；由引理 C.1，$P(H_0\in\Omega)=1$。由引理 C.3，
$\tilde U$ 服从 Haar 分布，与 $\lambda=\lambda(H_0)$ 独立，并且 $H_0=\tilde U\Lambda\tilde U^*$。
*$\lambda$ 的分布。* 设 $A\subseteq C$ 是 Borel 集。由引理 C.4(c) 和 $ST=S$，
除一个零概率事件外 $\{H_0\in\Psi(\mathfrak m_\varepsilon\times A)\}=\{\lambda\in A,\ \tilde U\in S\}$，其概率为
$P(\lambda\in A)\,h$，其中 $h=\int_K1_S\,dV$。另一方面，$\Psi$ 在开集 $\mathfrak m_\varepsilon\times C$ 上是单射、$C^1$ 的且 Jacobi 行列式非零，
并且 $p_0(\Psi(X,\lambda))=p_0(\Lambda)$；变量替换公式和引理 C.4(a) 给出
$$
P\big(H_0\in\Psi(\mathfrak m_\varepsilon\times A)\big)=\int_{\mathfrak m_\varepsilon}\!\int_Ap_0(\Lambda)\Delta(\lambda)^\beta j(X)\,d\lambda\,dX
=\kappa\int_Ap_0(\Lambda)\Delta(\lambda)^\beta\,d\lambda,\qquad \kappa=\int_{\mathfrak m_\varepsilon}j>0 .
$$
取 $A=C$ 得 $h=\kappa\int_Cp_0\Delta^\beta>0$。因此 $\lambda$ 在 $C$ 上具有密度 $c\,p_0(\Lambda)\Delta(\lambda)^\beta$，$c=\kappa/h$。
*一般的 $f$。* 由于 $p_0>0$，$\int f\,dH=\mathbb E\big[(f/p_0)(\tilde U\Lambda\tilde U^*)\big]$。$(\tilde U,\lambda)$ 的联合分布是
$dV\otimes c\,p_0(\Lambda)\Delta^\beta d\lambda$，并且 $p_0(V\Lambda V^*)=p_0(\Lambda)$；Tonelli 定理给出 $c_{n,\beta}=c$ 时的公式。
常数由公式本身确定（取 $f=p_0$）。$\square$

**推论 C.6（不变密度）。** 设 $H$ 在 $\mathcal H$ 上具有密度 $p$，且对所有
$V\in K$ 有 $p(VHV^*)=p(H)$。则 $\lambda(H)$ 在 $C$ 上具有密度 $c_{n,\beta}\,p(\Lambda)\Delta(\lambda)^\beta$；引理 C.3 的特征标架 $\tilde U$
服从 Haar 分布并与 $\lambda(H)$ 独立；并且对 $\mathbb R^n$ 上对称的 Borel 函数 $g\ge0$，
$$
\mathbb E\,g(\lambda(H))=\frac{c_{n,\beta}}{n!}\int_{\mathbb R^n}g(x)\,p(\operatorname{diag}x)\,|\Delta(x)|^\beta\,dx .
$$

*证明。* 由引理 C.1，$P(H\in\Omega)=1$。对 Borel 集 $A\subseteq C$，把定理 C.5 应用于 $f=1_A(\lambda(\cdot))\,p$，对它有
$f(V\Lambda V^*)=1_A(\lambda)p(\Lambda)$。$H$ 的分布是不变的，所以引理 C.3 适用。除去两个
坐标相等的零测集外，$\mathbb R^n$ 是 $C$ 在坐标置换下的 $n!$ 个像的不交并，并且
$g$、$|\Delta|$ 和 $x\mapsto p(\operatorname{diag}x)$ 都是置换不变的（置换矩阵属于 $K$）。$\square$

下面的检验对 $n=2,3$ 和两种 $\beta$ 数值地对 $\Psi$ 求导：在 $X=0$ 处 Jacobi 行列式是 $\Delta(\lambda)^\beta$，而
在固定的 $X\neq0$ 处它与 $\Delta(\lambda)^\beta$ 之比不依赖于 $\lambda$（引理 C.4(a)）。

```python
import math
import torch
from fractions import Fraction
from torch_vn_algebra import TypeIAlgebra, states

torch.manual_seed(0)
dt, f64 = torch.complex128, torch.float64

def mc(samples, exact, name):
    """Monte Carlo mean of real samples with its standard error; must agree within 5 standard errors."""
    mean, err = samples.mean().item(), samples.std().item() / math.sqrt(samples.numel())
    print(f"{name}: {mean:.4f} +- {err:.4f}   (theory {exact:.4f})")
    assert abs(mean - exact) <= 5 * err + 1e-12
    return mean

def vandermonde(lam):
    """Delta(lam) = prod_{i<j} (lam_j - lam_i) along the last axis."""
    n = lam.shape[-1]
    return math.prod(lam[..., j] - lam[..., i] for i in range(n) for j in range(i + 1, n))

def coords(H, beta):
    """Coordinates of dH: h_ii, Re h_ij (i < j) and, for beta = 2, Im h_ij."""
    iu = torch.triu_indices(H.shape[-1], H.shape[-1], 1)
    off = H[..., iu[0], iu[1]]
    return torch.cat([H.diagonal(dim1=-2, dim2=-1).real, off.real] + ([off.imag] if beta == 2 else []), dim=-1)

def from_m(c, n, beta):
    """The element X of m (X* = -X, zero diagonal) with coordinates Re X_ij (, Im X_ij), i < j."""
    iu = torch.triu_indices(n, n, 1)
    p = iu.shape[1]
    X = torch.zeros(n, n, dtype=dt)
    X[iu[0], iu[1]] = c[:p] + 1j * c[p:] if beta == 2 else c[:p].to(dt)
    return X - X.mH

def Psi(p, n, beta):                                       # (X, lambda) -> e^X Lambda e^{-X}
    d = beta * n * (n - 1) // 2
    E = torch.linalg.matrix_exp(from_m(p[:d], n, beta))
    return coords(E @ torch.diag(p[d:]).to(dt) @ E.mH, beta)

def jacobian(p, n, beta, h=1e-6):
    J = torch.stack([(Psi(p + h * e, n, beta) - Psi(p - h * e, n, beta)) / (2 * h)
                     for e in torch.eye(len(p), dtype=f64)], dim=1)
    return torch.linalg.det(J).abs().item()

for n in (2, 3):
    for beta in (1, 2):
        d = beta * n * (n - 1) // 2
        X = 0.4 * torch.randn(d, dtype=f64)
        ratios = []
        for trial in range(3):
            lam = torch.randn(n, dtype=f64).sort().values
            D = vandermonde(lam).item() ** beta
            assert abs(jacobian(torch.cat([torch.zeros(d, dtype=f64), lam]), n, beta) / D - 1) < 1e-6
            ratios.append(jacobian(torch.cat([X, lam]), n, beta) / D)
        print(f"n = {n}, beta = {beta}: j(X) = {ratios[0]:.6f} for a random X")
        assert max(ratios) - min(ratios) < 1e-6              # j(X) does not depend on lambda
```

## C.3 高斯系综

**定理 C.7（GOE 与 GUE）。** 设 $H=(G+G^T)/2$，其中 $G$ 是元素相互独立且服从 $N(0,1)$ 的实矩阵（$\beta=1$，GOE），或
$H=(G+G^*)/\sqrt2$，其中 $G$ 是 Ginibre 矩阵（$\beta=2$，GUE），如 5.3 节所述。则 $H$ 在 $\mathcal H$ 上具有密度
$(2\pi)^{-n/2}\pi^{-\beta n(n-1)/4}e^{-\operatorname{Tr}H^2/2}$，其特征值具有密度
$$
Z_{n,\beta}^{-1}\,\Delta(\lambda)^\beta\,e^{-\sum_i\lambda_i^2/2}\quad\text{on }C,\qquad Z_{n,\beta}=\int_C\Delta^\beta e^{-|\lambda|^2/2}\,d\lambda\in(0,\infty),
$$
并且特征标架 $\tilde U$ 在 $K$ 上服从 Haar 分布，与特征值独立。

*证明。* GOE：$h_{ii}=G_{ii}\sim N(0,1)$，对 $i<j$ 有 $h_{ij}=(G_{ij}+G_{ji})/2\sim N(0,\frac12)$。GUE：$h_{ii}=\sqrt2\operatorname{Re}G_{ii}\sim N(0,1)$，
$\operatorname{Re}h_{ij}=(\operatorname{Re}G_{ij}+\operatorname{Re}G_{ji})/\sqrt2$ 和 $\operatorname{Im}h_{ij}=(\operatorname{Im}G_{ij}-\operatorname{Im}G_{ji})/\sqrt2$ 服从 $N(0,\frac12)$。
不同的坐标是互不相交的独立变量组的函数，所以它们相互独立，
密度为 $\prod_i(2\pi)^{-1/2}e^{-h_{ii}^2/2}\prod_{i<j}\pi^{-\beta/2}e^{-|h_{ij}|^2}$，它就是所述的密度，因为
$\operatorname{Tr}H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2$。它是不变的，推论 C.6 适用；常数由
归一化确定。$Z_{n,\beta}<\infty$，因为 $|\lambda_j-\lambda_i|\le(1+|\lambda_i|)(1+|\lambda_j|)$，所以
$|\Delta(\lambda)|\le\prod_i(1+|\lambda_i|)^{n-1}$。$\square$

**推论 C.8（能级密度）。** 对大小为 $n$ 的 GOE 或 GUE，在 $\mathbb R^n$ 上令 $q(x)=(n!Z_{n,\beta})^{-1}|\Delta(x)|^\beta e^{-|x|^2/2}$，
并令 $\rho(x)=n\int_{\mathbb R^{n-1}}q(x,y)\,dy$。则对有界 Borel 函数
$\varphi$ 有 $\mathbb E\sum_i\varphi(\lambda_i(H))=\int\varphi\rho$，$\int\rho=n$，并且 $\rho$ 连续，满足 $0\le\rho(x)\le C_n(1+|x|)^{\beta(n-1)}e^{-x^2/2}$。特别地，
对相互独立的 GOE 或 GUE 块，定理 5.15 的假设成立。

*证明。* 由推论 C.6，取 $g(x)=\sum_i\varphi(x_i)$（对 $\varphi\ge0$；一般情形则应用于 $\varphi^\pm=\max(\pm\varphi,0)$），并利用 $q$ 的对称性，
$\mathbb E\sum_i\varphi(\lambda_i)=n\int\varphi(x_1)q(x)\,dx=\int\varphi\rho$（Fubini）。由 $|\Delta|$ 的界，
$q(x,y)\le(n!Z_{n,\beta})^{-1}(1+|x|)^{\beta(n-1)}e^{-x^2/2}\prod_{i}(1+|y_i|)^{\beta(n-1)}e^{-y_i^2/2}$，对 $y$ 积分后即得所述的界；
又由于 $\sup_x(1+|x|)^{\beta(n-1)}e^{-x^2/2}<\infty$，它还给出一个与 $x$ 无关的可积控制函数：连续性
由控制收敛定理得出。$\square$

**命题 C.9（最小间距）。** 对大小为 $n\ge2$ 的 GOE 或 GUE，令 $\delta=\min_i(\lambda_{i+1}-\lambda_i)$。存在
常数 $0<c\le C$（依赖于 $n,\beta$），使得对 $0<\varepsilon\le1$ 有 $c\,\varepsilon^{\beta+1}\le P(\delta<\varepsilon)\le C\varepsilon^{\beta+1}$。

*证明。* *上界。* $P(\delta<\varepsilon)\le\sum_iP(\lambda_{i+1}-\lambda_i<\varepsilon)$。在 $\{0<u=x_{i+1}-x_i<\varepsilon\}$ 上，把数对
$(i,i+1)$ 的因子用 $u^\beta$ 估计，其余因子如上估计，并利用 $1+|x_{i+1}|\le2(1+|x_i|)$ 和 $e^{-x_{i+1}^2/2}\le1$；代入
$x_{i+1}=x_i+u$，
$$
P(\lambda_{i+1}-\lambda_i<\varepsilon)\le\frac{2^{\beta(n-1)}}{Z_{n,\beta}}\int_0^\varepsilon u^\beta du\int(1+|x|)^{2\beta(n-1)}e^{-\frac{x^2}2}dx
\Big(\int(1+|x|)^{\beta(n-1)}e^{-\frac{x^2}2}dx\Big)^{n-2}.
$$
*下界。* 令 $R_\varepsilon=\{x:0<x_1<1,\ 0<x_2-x_1<\varepsilon,\ 3k<x_k<3k+1\ (3\le k\le n)\}\subseteq C$。在 $R_\varepsilon$ 上，
$\delta\le x_2-x_1$，除 $(1,2)$ 外的每一对的距离都 $\ge1$，并且 $|x|^2\le M$，其中 $M$ 与 $\varepsilon\le1$ 无关；所以
$P(\delta<\varepsilon)\ge P(\lambda\in R_\varepsilon)\ge Z_{n,\beta}^{-1}e^{-M/2}\int_0^\varepsilon u^\beta du$。$\square$

对 GOE，这就是 5.4 节末尾引用的 $\varepsilon^2$ 量级。我们对 $n=3$ 检验定理 C.7：GOE 的
能级密度由联合密度在网格上对另外两个特征值积分算出，并与
每个区间内特征值的平均个数比较；我们还检验特征向量与特征值不相关
（$|\tilde U_{11}|^2$ 不依赖于相位 $D$，所以可以从
`eigh` 返回的任何特征标架读出）。对服从 Haar 分布的 $U\in O(3)$，$\mathbb E|U_{11}|^2=\frac13$：左乘一个置换矩阵会置换
第一列的元素并保持 Haar 测度，所以三个 $\mathbb E|U_{i1}|^2$ 相等，而它们的和为 $1$。

```python
B, n = 200_000, 3
G = torch.randn(B, n, n, dtype=f64)
goe_lam, goe_vec = torch.linalg.eigh((G + G.mT) / 2)                    # Section 5.3 normalisation
edges = torch.arange(-4.0, 4.01, 1.0, dtype=f64)

h = 0.04                                                                # grid for the other two eigenvalues
xs = torch.arange(-6, 6 + h / 2, h, dtype=f64)
Y1, Y2 = torch.meshgrid(xs, xs, indexing='ij')
joint = lambda x1, x2, x3: (vandermonde(torch.stack([x1, x2, x3], -1)).abs()
                            * torch.exp(-(x1 ** 2 + x2 ** 2 + x3 ** 2) / 2))
marg = torch.stack([joint(x + 0 * Y1, Y1, Y2).sum() * h * h for x in xs])   # int |Delta| e^{-|x|^2/2} dy dz
rho_goe = 3 * marg / (marg.sum() * h)                                   # Corollary C.8, normalised by n! Z
for lo, hi in zip(edges[:-1], edges[1:]):
    sel = (xs >= lo - 1e-9) & (xs <= hi + 1e-9)
    count = ((goe_lam >= lo) & (goe_lam < hi)).sum(-1).double()
    mc(count, torch.trapezoid(rho_goe[sel], xs[sel]).item(), f"GOE(3): E #eigenvalues in [{lo:.0f}, {hi:.0f})")

u11 = goe_vec[:, 0, 0] ** 2                                             # |U_11|^2, Haar: mean 1/3
mc(u11, 1 / 3, "GOE(3): E |U_11|^2")
mc((u11 - 1 / 3) * goe_lam[:, 0], 0.0, "GOE(3): E (|U_11|^2 - 1/3) lambda_1")   # independence (Theorem C.7)
```

## C.4 正交多项式与单点函数

**引理 C.10（单点函数）。** 设 $I\subseteq\mathbb R$ 是区间，$w\ge0$ 是 $I$ 上的 Borel 函数，在一个
正测度集上为正，并且对所有 $m$ 有 $\int_I|x|^mw<\infty$。设 $p_0,\dots,p_{n-1}$ 是实多项式，$\deg p_k=k$，
首项系数 $\kappa_k\neq0$，在 $L^2(w)$ 中正交归一（它们存在：对 $1,x,\dots,x^{n-1}$ 应用 Gram–Schmidt 过程，这些函数
在 $L^2(w)$ 中线性无关，因为非零多项式只有有限多个零点；每个 $p_k$ 的符号可以任意选取）。在 $I^n$ 上令 $q=Z_n^{-1}\Delta(x)^2\prod_iw(x_i)$，
其中 $Z_n=\int_{I^n}\Delta^2\prod_iw(x_i)\,dx$。则 $Z_n=n!\prod_k\kappa_k^{-2}$，并且
$$
n\int_{I^{n-1}}q(x,y)\,dy=K_n(x)\,w(x),\qquad K_n(x)=\sum_{k=0}^{n-1}p_k(x)^2=v(x)^TM^{-1}v(x),
$$
其中 $v(x)=(1,x,\dots,x^{n-1})$，$M_{ij}=\int_Ix^{i+j}w$，$0\le i,j<n$。

*证明。* $\Delta(x)=\det(x_j^{i-1})_{i,j}$。对每一行加上其上方各行的倍数，可把 $x_j^{i-1}$ 换成
$p_{i-1}(x_j)/\kappa_{i-1}$，所以 $\Delta(x)=(\prod_k\kappa_k)^{-1}\det(p_{i-1}(x_j))$，并且
$$
\Delta(x)^2=\Big(\prod_k\kappa_k\Big)^{-2}\sum_{\sigma,\tau\in S_n}\operatorname{sgn}\sigma\operatorname{sgn}\tau\prod_{j=1}^np_{\sigma(j)-1}(x_j)\,p_{\tau(j)-1}(x_j).
$$
对 $x_2,\dots,x_n$ 关于 $\prod_{j\ge2}w(x_j)$ 积分，因子 $j$ 给出 $\delta_{\sigma(j)\tau(j)}$，所以只有 $\sigma=\tau$
的项保留下来，结果为 $(\prod_k\kappa_k)^{-2}\sum_\sigma p_{\sigma(1)-1}(x_1)^2=(\prod_k\kappa_k)^{-2}(n-1)!\,K_n(x_1)$。
再关于 $w(x_1)$ 积分得 $Z_n=(\prod_k\kappa_k)^{-2}n!$，相除即得公式。最后，若
$p_k(x)=\sum_iC_{ki}x^i$，则正交归一性写成 $CMC^T=1$，其中 $C$ 可逆，所以 $M^{-1}=C^TC$，并且
$v^TM^{-1}v=|Cv|^2=\sum_kp_k^2$。$\square$

由推论 C.6，GUE 的无序特征值具有引理 C.10 中的密度 $q$，其中 $w(x)=e^{-x^2/2}$，定义在
$\mathbb R$ 上，所以 GUE 的能级密度是 $K_n(x)e^{-x^2/2}$。我们对 $n=3$ 用矩
$\int x^{2m}e^{-x^2/2}dx=\sqrt{2\pi}(2m-1)(2m-3)\cdots1$ 来检验这一点：

```python
def one_point(x, n, moment, w):
    """K_n(x) w(x) = v(x)^T M^{-1} v(x) w(x), M the matrix of moments (Lemma C.10)."""
    M = torch.tensor([[moment(i + j) for j in range(n)] for i in range(n)], dtype=f64)
    v = torch.stack([x ** i for i in range(n)], dim=-1)
    return (v @ torch.linalg.inv(M) * v).sum(-1) * w(x)

gauss_moment = lambda m: math.sqrt(2 * math.pi) * math.prod(range(m - 1, 0, -2)) if m % 2 == 0 else 0.0
Gc = torch.randn(B, n, n, dtype=dt)
gue_lam = torch.linalg.eigvalsh((Gc + Gc.mH) / math.sqrt(2))
for lo, hi in zip(edges[:-1], edges[1:]):
    x = torch.linspace(lo, hi, 401, dtype=f64)
    exact = torch.trapezoid(one_point(x, n, gauss_moment, lambda t: torch.exp(-t ** 2 / 2)), x).item()
    mc(((gue_lam >= lo) & (gue_lam < hi)).sum(-1).double(), exact, f"GUE(3): E #eigenvalues in [{lo:.0f}, {hi:.0f})")
```

## C.5 Wishart 矩阵与诱导测度

设 $a\le b$，$G$ 是 $a\times b$ 的 Ginibre 矩阵（元素相互独立，其实部与虚部相互独立且服从 $N(0,\frac12)$，
密度为 $\pi^{-ab}e^{-\operatorname{Tr}GG^*}$），$W=GG^*$。设 $T_+$ 是对角元为正的 $a\times a$ 下三角矩阵
构成的集合，它在坐标 $L_{ii}$、$\operatorname{Re}L_{ij}$、$\operatorname{Im}L_{ij}$（$i>j$）中是 $\mathbb R^{a^2}$ 的开子集；
$P_a$ 是正定矩阵构成的集合，它在 $a\times a$ Hermite 矩阵中是开集。

**引理 C.11（Bartlett 分解）。** 几乎必然地，$G$ 的各行 $r_1,\dots,r_a$ 线性无关，并且
Gram–Schmidt 过程给出 $G=LV$，其中 $L\in T_+$，$V$ 的各行正交归一，所以 $W=LL^*$。元素 $L_{ij}$（$i\ge j$）
相互独立，$L_{ij}$（$i>j$）是标准复高斯变量，而 $L_{ii}^2$ 在 $(0,\infty)$ 上具有密度 $x^{b-i}e^{-x}/(b-i)!$。
因此 $L$ 在 $T_+$ 上具有密度 $c_L\prod_iL_{ii}^{2(b-i)+1}e^{-\operatorname{Tr}LL^*}$。

*证明。* 令 $\mathcal F_{i-1}=\sigma(r_1,\dots,r_{i-1})$，并假设（归纳）$r_1,\dots,r_{i-1}$ 几乎必然线性无关，其
Gram–Schmidt 向量为 $v_1,\dots,v_{i-1}$。设 $u_1,\dots,u_b$ 是对
$v_1,\dots,v_{i-1},e_1,\dots,e_b$ 应用 Gram–Schmidt 过程（舍去线性相关的向量）得到的正交归一基；对 $k<i$ 有 $u_k=v_k$，并且各 $u_k$ 是
$r_1,\dots,r_{i-1}$ 的 Borel 函数。令 $c_k=\langle u_k,r_i\rangle$。由于 $r_i$ 与 $\mathcal F_{i-1}$ 独立，且其密度 $\pi^{-b}e^{-\|r\|^2}$
在酉变换下不变（命题 A.16(a)），在给定 $\mathcal F_{i-1}$ 的条件下，向量 $c$ 是标准复高斯向量（Fubini）。此时
对 $j<i$ 有 $L_{ij}=c_j$，并且 $r_i-\sum_{j<i}c_jv_j=\sum_{k\ge i}c_ku_k$，所以几乎必然 $L_{ii}^2=\sum_{k=i}^b|c_k|^2>0$。因此，给定 $\mathcal F_{i-1}$ 时
$L$ 的第 $i$ 行的条件分布是一个固定的分布，所以第 $i$ 行与前面各行独立，并且其坐标
$c_1,\dots,c_{i-1}$ 与 $L_{ii}$ 相互独立；由归纳，$L$ 的所有行相互独立。由定理 A.20(b)（取 $N=m=b-i+1$），
$L_{ii}^2=\|(c_i,\dots,c_b)\|^2$ 具有密度 $s^{m-1}e^{-s}/(m-1)!$。于是 $L_{ii}=\sqrt s$ 的密度是
$2t^{2(b-i)+1}e^{-t^2}/(b-i)!$，并且 $\operatorname{Tr}LL^*=\sum_{i\ge j}|L_{ij}|^2$。$\square$

**引理 C.12（Cholesky 映射）。** $\chi(L)=LL^*$ 是从 $T_+$ 到 $P_a$ 上的 $C^\infty$ 双射，并且
$|\det D\chi(L)|=2^a\prod_{i=1}^aL_{ii}^{2(a-i)+1}$。

*证明。* *满射：* 对 $W\in P_a$，把引理 5.2 应用于 $W^{1/2}=QR$；则 $W=R^*R$，且 $L=R^*\in T_+$。*单射：* 若
$LL^*=MM^*$，则 $M^{-1}L=M^*(L^*)^{-1}$ 既是下三角的又是上三角的，因而是对角元为正的对角矩阵，并且是
酉的，因而等于 $1$。*Jacobi 行列式：* 对 $i\ge j$，$W_{ij}=L_{ij}L_{jj}+\sum_{k<j}L_{ik}\overline{L_{jk}}$。把 $L$ 的坐标和
$W$ 的下三角部分的坐标逐列排列，并在第 $j$ 列内把对角元排在最前。于是 $W_{ij}$ 只依赖于
$L_{ij}$ 和排在前面的坐标，所以 $D\chi$ 是三角矩阵，其对角元为 $2L_{jj}$（$\partial W_{jj}/\partial L_{jj}$），以及
对每个 $i>j$ 两次出现的 $L_{jj}$（$W_{ij}$ 的实部与虚部）。乘积为 $\prod_j2L_{jj}\cdot L_{jj}^{2(a-j)}$。$\square$

**定理 C.13（Wishart 系综与 Laguerre 系综）。** $W=GG^*$ 在 $P_a$ 上具有密度 $c_W\det(W)^{b-a}e^{-\operatorname{Tr}W}$
（在其外为零），其特征值具有密度
$$
Z_{a,b}^{-1}\,\Delta(\lambda)^2\prod_{i=1}^a\lambda_i^{\,b-a}e^{-\lambda_i}\quad\text{on }C_+=\{0<\lambda_1<\dots<\lambda_a\}.
$$

*证明。* 对 Borel 函数 $g\ge0$，引理 C.11 和变量替换 $W=\chi(L)$（引理 C.12）给出
$\mathbb E\,g(W)=\int_{P_a}g(W)\,c_L\prod_iL_{ii}^{2(b-i)+1}e^{-\operatorname{Tr}W}\big/\big(2^a\prod_iL_{ii}^{2(a-i)+1}\big)dW$，其中
$L=\chi^{-1}(W)$，并且 $\prod_iL_{ii}^{2(b-a)}=\det(W)^{b-a}$。这个密度在 $W\mapsto VWV^*$（$V\in U(a)$）下不变，并且
$\Lambda\in P_a$ 当且仅当 $\lambda_1>0$；以 $\beta=2$ 应用推论 C.6 即可。$\square$

**定理 C.14（诱导测度）。** 设 $\psi$ 是 $\mathbb C^a\otimes\mathbb C^b$ 的 Haar 随机纯态（7.5 节），
$\rho_A=\operatorname{Tr}_2|\psi\rangle\langle\psi|=W/s$，其中 $s=\operatorname{Tr}W$（定理 7.11 的证明），
$\mu_1<\dots<\mu_a$ 是 $\rho_A$ 的特征值。则 $s$ 与 $\mu$ 相互独立，$s$ 具有密度
$s^{ab-1}e^{-s}/(ab-1)!$，并且 $(\mu_1,\dots,\mu_{a-1})$ 具有密度
$$
Z'^{-1}_{a,b}\,\Delta(\mu)^2\prod_{i=1}^a\mu_i^{\,b-a}\quad\text{on }D_a=\{0<\mu_1<\dots<\mu_{a-1}<\mu_a\},\qquad\mu_a=1-\sum_{i<a}\mu_i ,
$$
这就是物理文献中的密度"$\propto\delta(1-\sum_i\mu_i)\,\Delta(\mu)^2\prod_i\mu_i^{b-a}$"。

*证明。* $\Theta(s,\mu_1,\dots,\mu_{a-1})=s\,(\mu_1,\dots,\mu_{a-1},1-\sum_{i<a}\mu_i)$ 是从 $(0,\infty)\times D_a$ 到
$C_+$ 上的 $C^\infty$ 双射，其逆为 $\lambda\mapsto(\sum\lambda_i,\lambda_1/\sum\lambda_i,\dots)$。把 $D\Theta$ 的所有行加到最后一行上，最后一行变为
$(1,0,\dots,0)$，所以 $|\det D\Theta|=s^{a-1}$。由于 $\Delta$ 是 $a(a-1)/2$ 次齐次的，$(s,\mu)$ 的密度为
$$
Z_{a,b}^{-1}\,s^{a(a-1)+a(b-a)+a-1}e^{-s}\cdot\Delta(\mu)^2\prod_i\mu_i^{b-a}=Z_{a,b}^{-1}\,s^{ab-1}e^{-s}\cdot\Delta(\mu)^2\prod_i\mu_i^{b-a},
$$
这是一个乘积，从而证明了独立性，并在归一化后确定了两个因子。$\square$

我们在 `torch.linalg.cholesky` 计算出的 Cholesky 因子上检验引理 C.11（由引理 C.12，它就是唯一的 $L\in T_+$），
并对 $a=2$ 检验定理 C.14，此时 $\mu_1\in(0,\frac12)$ 的密度 $\propto(1-2t)^2(t(1-t))^{b-2}$：

```python
a, b, B = 3, 5, 100_000
G = torch.randn(B, a, b, dtype=dt)                        # E|G_ij|^2 = 1
W = G @ G.mH
L = torch.linalg.cholesky(W)                              # W = L L*, L in T_+
for i in range(a):                                        # 0-based i: L_ii^2 ~ Gamma(b - i) (Lemma C.11, index i + 1)
    mc(L[:, i, i].real ** 2, b - i, f"E L_{i + 1}{i + 1}^2")
    mc(L[:, i, i].real ** 4, (b - i) * (b - i + 1), f"E L_{i + 1}{i + 1}^4")
mc(L[:, 2, 0].abs() ** 2 * L[:, 1, 1].real ** 2, 1.0 * (b - 1), "E |L_31|^2 L_22^2")   # independence

for b2 in (2, 4):
    rho = states.random_density_matrix(TypeIAlgebra([2], [2], precision='double'), 100_000, rank=b2)
    mu1 = rho.eigenvalues()[0][:, 0]                      # smallest eigenvalue of rho_A
    t = torch.linspace(0, 0.5, 20_001, dtype=f64)
    dens = (1 - 2 * t) ** 2 * (t * (1 - t)) ** (b2 - 2)
    dens = dens / torch.trapezoid(dens, t)
    for lo, hi in ((0.0, 0.1), (0.1, 0.25), (0.25, 0.4)):
        sel = (t >= lo) & (t <= hi)
        mc(((mu1 >= lo) & (mu1 < hi)).double(), torch.trapezoid(dens[sel], t[sel]).item(),
           f"a = 2, b = {b2}: P({lo} <= mu_1 < {hi})")
```

## C.6 Laguerre 多项式与 Page 公式

**引理 C.15（Laguerre 多项式）。** 对整数 $\alpha,k\ge0$，令
$L_k^{(\alpha)}(x)=\sum_{j=0}^k(-1)^j\binom{k+\alpha}{k-j}\frac{x^j}{j!}$，并令 $L_{-1}^{(\alpha)}=0$。则：

(a) $L_k^{(\alpha)}(x)\,x^\alpha e^{-x}=\frac1{k!}\frac{d^k}{dx^k}\big(x^{k+\alpha}e^{-x}\big)$；$L_k^{(\alpha)}$ 的次数为 $k$，首项系数为 $(-1)^k/k!$。

(b) $\frac{d}{dx}L_k^{(\alpha)}=-L_{k-1}^{(\alpha+1)}$，因此对 $i\le k$ 有 $\frac{d^i}{dx^i}L_k^{(\alpha)}=(-1)^iL_{k-i}^{(\alpha+i)}$；并且 $L_k^{(\alpha)}=L_k^{(\alpha+1)}-L_{k-1}^{(\alpha+1)}$。

(c) 对 $m<k$ 有 $\int_0^\infty x^mL_k^{(\alpha)}x^\alpha e^{-x}dx=0$，并且 $\int_0^\infty L_k^{(\alpha)}L_l^{(\alpha)}x^\alpha e^{-x}dx=\delta_{kl}\frac{(k+\alpha)!}{k!}$。

因此，由引理 C.10 和定理 C.13，$W$ 的特征值具有能级密度
（对 Borel 函数 $\varphi\ge0$，以及对每个满足 $\int|\varphi|\rho_{a,b}<\infty$ 的 $\varphi$，有 $\mathbb E\sum_i\varphi(\lambda_i(W))=\int_0^\infty\varphi\,\rho_{a,b}$）
$$
\rho_{a,b}(x)=\sum_{k=0}^{a-1}\frac{k!}{(k+\alpha)!}L_k^{(\alpha)}(x)^2\,x^\alpha e^{-x},\qquad\alpha=b-a .
$$

*证明。* (a) 由 Leibniz 公式，$\frac{d^k}{dx^k}(x^{k+\alpha}e^{-x})=\sum_j\binom kj(-1)^j\frac{(k+\alpha)!}{(\alpha+j)!}x^{\alpha+j}e^{-x}$（其中 $j$ 次导数
作用在 $e^{-x}$ 上），并且 $\frac1{k!}\binom kj\frac{(k+\alpha)!}{(\alpha+j)!}=\binom{k+\alpha}{k-j}\frac1{j!}$。(b) 逐项求导，
$-\sum_{i=0}^{k-1}(-1)^i\binom{k+\alpha}{k-1-i}\frac{x^i}{i!}=-L_{k-1}^{(\alpha+1)}$；第二个恒等式是 Pascal 法则
$\binom{k+\alpha+1}{k-j}-\binom{k+\alpha}{k-1-j}=\binom{k+\alpha}{k-j}$。(c) 把 (a) 与 $x^m$ 相乘后分部积分 $k$ 次：
边界项为零，因为 $\frac{d^r}{dx^r}(x^{k+\alpha}e^{-x})$（$r<k$）含有因子 $x^{k+\alpha-r}$，并且在 $\infty$ 处衰减。
结果是 $\frac{(-1)^k}{k!}\int\frac{d^kx^m}{dx^k}x^{k+\alpha}e^{-x}dx$，它在 $m<k$ 时为 $0$，在 $m=k$ 时为 $(-1)^k(k+\alpha)!$。对 $l<k$，
$L_l^{(\alpha)}$ 是 $x^m$（$m<k$）的线性组合；对 $l=k$，只有首项系数 $(-1)^k/k!$ 起作用。所以
$p_k=\big(k!/(k+\alpha)!\big)^{1/2}L_k^{(\alpha)}$ 关于 $(0,\infty)$ 上的 $w(x)=x^\alpha e^{-x}$ 正交归一，并且由推论 C.6，$W$ 的无序
特征值具有引理 C.10（取 $n=a$）中的密度 $q$，于是与推论 C.8 的证明一样，
对 Borel 函数 $\varphi\ge0$ 有 $\mathbb E\sum_i\varphi(\lambda_i)=a\int\varphi(x_1)q(x)\,dx=\int\varphi\,K_aw=\int\varphi\,\rho_{a,b}$，再由线性性推广到
满足 $\int|\varphi|\rho_{a,b}<\infty$ 的 $\varphi$。$\square$

**引理 C.16（对数积分）。** 令 $I_m=\int_0^\infty x^m\log x\,e^{-x}dx$，$H_m=\sum_{i=1}^m\frac1i$。则
$I_m=m!\,(I_0+H_m)$。对整数 $\alpha\ge0$ 和 $0\le k\le m$，
$$
A^{(\alpha)}_{k,m}:=\int_0^\infty\log x\,L_k^{(\alpha)}(x)L_m^{(\alpha)}(x)\,x^\alpha e^{-x}dx=
\begin{cases}I_{k+\alpha}/k!, & m=k,\\ -\dfrac{(k+\alpha)!}{k!\,(m-k)}, & m>k.\end{cases}
$$
（$I_0$ 在定理 C.17 中相互抵消，所以不需要它的值。）

*证明。* 对 $m\ge1$，利用 $(-e^{-x})'=e^{-x}$ 分部积分得 $I_m=mI_{m-1}+(m-1)!$，即
$I_m/m!=I_{m-1}/(m-1)!+1/m$。对 $A_{k,m}$，记 $L_m^{(\alpha)}x^\alpha e^{-x}=\frac1{m!}g^{(m)}$，其中 $g=x^{m+\alpha}e^{-x}$
（引理 C.15(a)），并令 $f=\log x\,L_k^{(\alpha)}$，然后分部积分 $m$ 次。边界项
$f^{(i)}g^{(m-1-i)}$ 在 $\infty$ 处和在 $0$ 处都为零，在原点附近 $f^{(i)}=O(|\log x|+x^{-i})$ 且 $g^{(m-1-i)}=O(x^{\alpha+1+i})$。所以
$A_{k,m}=\frac{(-1)^m}{m!}\int_0^\infty f^{(m)}g\,dx$（被积函数在 $0$ 处为 $O(x^\alpha(1+|\log x|))$），其中
$f^{(m)}=\sum_{i\le k}\binom mi(\log x)^{(m-i)}\frac{d^i}{dx^i}L_k^{(\alpha)}$。当 $m-i\ge1$ 时，$(\log x)^{(m-i)}=(-1)^{m-i-1}(m-i-1)!\,x^{i-m}$，
并且由引理 C.15(b,c)，
$$
\int_0^\infty x^{i-m}\frac{d^iL_k^{(\alpha)}}{dx^i}\,x^{m+\alpha}e^{-x}dx=(-1)^i\int_0^\infty L_{k-i}^{(\alpha+i)}\,x^{\alpha+i}e^{-x}dx=(-1)^i(\alpha+i)!\,\delta_{ik}.
$$
若 $m>k$，只有 $i=k$ 有贡献：$A_{k,m}=\frac{(-1)^m}{m!}\binom mk(-1)^{m-k-1}(m-k-1)!(-1)^k(\alpha+k)!=-\frac{(k+\alpha)!}{k!(m-k)}$。
若 $m=k$，$i<k$ 的项为零，而 $i=k$ 的项是 $\log x\cdot(-1)^k$，从而 $A_{k,k}=\frac{(-1)^k}{k!}(-1)^kI_{k+\alpha}$。$\square$

**定理 C.17（Page 公式）。** 对 $\mathbb C^a\otimes\mathbb C^b$（$a\le b$）的 Haar 随机纯态，
$$
\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}.
$$

*证明。* *化归。* 采用定理 C.14 的记号，$\lambda_i=s\mu_i$ 且 $\sum_i\mu_i=1$，所以
$\sum_i\lambda_i\log\lambda_i=s\log s-s\,S(\rho_A)$。所有项都可积（$|x\log x|\le1+x^2$，并且 $\mathbb Es^2<\infty$，
$0\le S(\rho_A)\le a/e$，因为在 $[0,1]$ 上 $0\le-x\log x\le1/e$），由引理 C.15 有 $\mathbb E\sum_i\lambda_i\log\lambda_i=\int x\log x\,\rho_{a,b}(x)\,dx$，并且 $s$ 与
$S(\rho_A)$ 独立，所以 $\mathbb E\sum_i\lambda_i\log\lambda_i=\mathbb E[s\log s]-\mathbb E[s]\,\mathbb E\,S(\rho_A)$。由 $s$ 的密度，
$\mathbb Es=ab$，并且 $\mathbb E[s\log s]=I_{ab}/(ab-1)!=ab\,(I_0+H_{ab})$（引理 C.16）。因此
$$
\mathbb E\,S(\rho_A)=I_0+H_{ab}-\frac1{ab}\int_0^\infty x\log x\,\rho_{a,b}(x)\,dx .
$$
*积分。* 由引理 C.15(b)，$L_k^{(\alpha)}=L_k^{(\alpha+1)}-L_{k-1}^{(\alpha+1)}$，所以利用 $x\cdot x^\alpha=x^{\alpha+1}$ 以及对
$\alpha+1$ 应用的引理 C.16，
$$
\frac{k!}{(k+\alpha)!}\int_0^\infty x\log x\,L_k^{(\alpha)}(x)^2x^\alpha e^{-x}dx=\frac{k!}{(k+\alpha)!}\Big(A^{(\alpha+1)}_{k,k}-2A^{(\alpha+1)}_{k-1,k}+A^{(\alpha+1)}_{k-1,k-1}\Big)
$$
$$
=(k+\alpha+1)(I_0+H_{k+\alpha+1})+2k+k\,(I_0+H_{k+\alpha}),
$$
这对 $k=0$ 也成立。由于 $\sum_{k=0}^{a-1}(2k+\alpha+1)=a(a-1)+a(b-a+1)=ab$，
$$
\int_0^\infty x\log x\,\rho_{a,b}\,dx=ab\,I_0+\Sigma,\qquad \Sigma=\sum_{k=0}^{a-1}\big[(k+\alpha+1)H_{k+\alpha+1}+kH_{k+\alpha}+2k\big].
$$
*求和。* 令 $m=\alpha+1=b-a+1$，$j=k+m\in\{m,\dots,b\}$。由于 $(j-m)H_{j-1}=(j-m)H_j-1+\frac mj$，
$\Sigma=\sum_{j=m}^b(2j-m)H_j-a+m(H_b-H_{m-1})+a(a-1)$。令 $F(j)=j(j+1-m)$，则 $F(j)-F(j-1)=2j-m$，$F(m-1)=0$，
$F(b)=ab$，分部求和给出
$$
\sum_{j=m}^b(2j-m)H_j=abH_b-\sum_{j=m}^{b-1}\frac{F(j)}{j+1}=abH_b-\sum_{j=m}^{b-1}\Big(j-m+\frac m{j+1}\Big)=abH_b-\frac{(a-1)(a-2)}2-m(H_b-H_m).
$$
由于 $m(H_m-H_{m-1})=1$，$\Sigma=abH_b-\frac{(a-1)(a-2)}2+1-a+a(a-1)=abH_b+\frac{a(a-1)}2$。因此
$\mathbb E\,S(\rho_A)=H_{ab}-H_b-\frac{a-1}{2b}$。$\square$

**推论 C.18.** 对 $a\le b$，
$\log a-\frac{a-1}{ab+1}-\frac{a-1}{2b}\le\mathbb E\,S(\rho_A)\le\log a-\frac{a-1}{2b}$；特别地，
$\big|\mathbb E\,S(\rho_A)-\log a+\frac a{2b}\big|\le\frac1{2b}$。对两个量子比特，$\mathbb E\,S(\rho_A)=\frac13$。

*证明。* 由于 $1/x$ 递减，$\int_k^{k+1}\frac{dx}x\le\frac1k\le\int_{k-1}^k\frac{dx}x$；对 $b<k\le ab$ 求和，
$\log\frac{ab+1}{b+1}\le H_{ab}-H_b\le\log a$，并且 $\log\frac{ab+1}{b+1}=\log a-\log\big(1+\frac{a-1}{ab+1}\big)\ge\log a-\frac{a-1}{ab+1}$。
最后 $\frac{a-1}{ab+1}\le\frac1b$，而 $a=b=2$ 给出 $\frac13+\frac14-\frac14$。$\square$

该公式由 D. N. Page 提出猜想（Phys. Rev. Lett. 71, 1291, 1993），并由 S. K. Foong 和 S. Kanno
（Phys. Rev. Lett. 72, 1148, 1994）、J. Sánchez-Ruiz（Phys. Rev. E 52, 5653, 1995）以及 S. Sen（Phys. Rev. Lett. 77, 1, 1996）证明；
上面的计算沿用 Sen 化归到 Laguerre 系综的做法，但逐项对单点密度求和，
而不使用 Christoffel–Darboux 公式。诱导测度来自 K. Życzkowski 和 H.-J. Sommers
（J. Phys. A 34, 7111, 2001）；引理 C.11 的分解归功于 M. S. Bartlett（1933），定理 C.5 归功于 H. Weyl。

我们检验：Laguerre 能级密度与直方图的比较；引理 C.16 的值 $A^{(\alpha)}_{k,m}$ 以及中间
结果 $\int x\log x\,\rho_{a,b}=ab\,I_0+\Sigma$（用数值积分，其中 $I_0$ 本身也用数值积分计算）；$\Sigma$ 的恒等式
（精确算术）；Page 公式，与 $a=2$ 时诱导密度的直接积分（精确有理数值）以及蒙特卡罗结果比较；
最后是推论 C.18：

```python
def laguerre(k, al, x):
    return sum((-1) ** j * math.comb(k + al, k - j) * x ** j / math.factorial(j) for j in range(k + 1))

def rho_laguerre(x, a, b):                                # Lemma C.15
    al = b - a
    return sum(math.factorial(k) / math.factorial(k + al) * laguerre(k, al, x) ** 2 for k in range(a)) * x ** al * torch.exp(-x)

lam = torch.linalg.eigvalsh(W)                            # a = 3, b = 5 from above
x = torch.linspace(0, 80, 160_001, dtype=f64)
r = rho_laguerre(x, a, b)
assert abs(torch.trapezoid(r, x).item() - a) < 1e-8 and abs(torch.trapezoid(x * r, x).item() - a * b) < 1e-8
for lo, hi in ((0, 1), (1, 2), (2, 4), (4, 8), (8, 12), (12, 20)):
    sel = (x >= lo) & (x <= hi)
    mc(((lam >= lo) & (lam < hi)).sum(-1).double(), torch.trapezoid(r[sel], x[sel]).item(),
       f"Laguerre(3, 5): E #eigenvalues in [{lo}, {hi})")

def Hn(m):
    return sum(Fraction(1, i) for i in range(1, m + 1))

def page(a, b):
    return Hn(a * b) - Hn(b) - Fraction(a - 1, 2 * b)

vv = torch.linspace(0, 12, 1_200_001, dtype=f64)[1:]      # x = v^2 smooths log x at 0
def log_integral(g):                                      # int_0^oo log(x) g(x) dx
    xv = vv ** 2
    return torch.trapezoid(xv.log() * g(xv) * 2 * vv, vv).item()
I0 = log_integral(lambda t: torch.exp(-t))
print(f"I_0 = {I0:.10f} (quadrature)")
for al in range(3):                                       # Lemma C.16
    for k in range(4):
        for m in range(k, 5):
            A = log_integral(lambda t: laguerre(k, al, t) * laguerre(m, al, t) * t ** al * torch.exp(-t))
            Im = math.factorial(k + al) * (I0 + float(Hn(k + al)))
            assert abs(A - (Im / math.factorial(k) if m == k else -math.factorial(k + al) / (math.factorial(k) * (m - k)))) < 1e-6
for a_, b_ in ((2, 2), (2, 3), (3, 3), (3, 5), (4, 6)):
    al = b_ - a_
    Sigma = sum((k + al + 1) * Hn(k + al + 1) + k * Hn(k + al) + 2 * k for k in range(a_))
    assert Sigma == a_ * b_ * Hn(b_) + Fraction(a_ * (a_ - 1), 2)            # proof of Theorem C.17, exactly
    integral = log_integral(lambda t: t * rho_laguerre(t, a_, b_))
    assert abs(integral - (a_ * b_ * I0 + float(Sigma))) < 1e-6

print("a = 2:", ", ".join(f"b = {b_}: {page(2, b_)}" for b_ in range(2, 7)))
assert page(2, 2) == Fraction(1, 3)
for b_ in range(2, 7):                                    # a = 2: integrate the density of Theorem C.14
    v = torch.linspace(0, math.sqrt(0.5), 200_001, dtype=f64)[1:]
    t = v ** 2                                            # mu_1 = v^2 smooths the logarithm at 0
    dens = (1 - 2 * t) ** 2 * (t * (1 - t)) ** (b_ - 2) * 2 * v
    S = -(t * t.log() + (1 - t) * (1 - t).log())
    assert abs((torch.trapezoid(dens * S, v) / torch.trapezoid(dens, v)).item() - float(page(2, b_))) < 1e-9

for a_, b_ in ((2, 2), (2, 3), (3, 3), (3, 5)):
    rho = states.random_density_matrix(TypeIAlgebra([a_], [a_], precision='double'), 40_000, rank=b_)
    mc(rho.von_neumann_entropy(), float(page(a_, b_)), f"E S(rho_A), a = {a_}, b = {b_}")   # Theorem C.17

for a_ in range(1, 12):
    for b_ in range(a_, 40):
        p = float(page(a_, b_))
        assert math.log(a_) - (a_ - 1) / (a_ * b_ + 1) - (a_ - 1) / (2 * b_) - 1e-12 <= p <= math.log(a_) - (a_ - 1) / (2 * b_) + 1e-12
        assert abs(p - math.log(a_) + a_ / (2 * b_)) <= 1 / (2 * b_) + 1e-12  # Corollary C.18
```

**返回：** [第 5 章](05_random_matrices.md)（定理 C.7、推论 C.8 和命题 C.9），
[第 7 章](07_tensor_products.md)（定理 C.14、C.17 和推论 C.18）。
