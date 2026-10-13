[English](../A_analysis.md) · [Español](../es/A_analysis.md) · [Français](../fr/A_analysis.md) · [Deutsch](../de/A_analysis.md) · **中文** · [日本語](../ja/A_analysis.md) · [Русский](../ru/A_analysis.md)

# 附录 A. 分析与概率中的若干事实

各章用到了少数几条分析与概率中的事实，它们超出了本课程的预备知识，但本身又算不上独立的
定理：多项式的零点集是零测集；Hermite 矩阵按大小排列的特征值连续地依赖于该矩阵；高斯向量的
非退化线性像具有密度；高斯向量的方向在球面上均匀分布并且与其长度独立；高斯矩是对配对的求和；
均匀分布点的间距是归一化的指数变量。本附录证明这些事实，以便各章按编号引用它们。

**所用的预备知识。** 线性代数；实分析（极限、连续性、紧性以及 $\mathbb R^n$ 中的
Bolzano–Weierstrass 定理、单变量幂级数、连续函数的 Riemann 积分、变量替换公式）；测度与积分
（Lebesgue 测度、乘积测度与 Fubini–Tonelli 定理、单调收敛定理与控制收敛定理、前推测度、
密度、独立性）；初等概率论（分布 $N(0,s^2)$，其密度为
$(2\pi s^2)^{-1/2}e^{-y^2/2s^2}$，矩母函数为 $\mathbb Ee^{aY}=e^{a^2s^2/2}$）。对 $\mathbb R^m$ 的
可逆线性映射 $L$ 和 Borel 集 $E$，有 $\lambda_m(LE)=|\det L|\,\lambda_m(E)$，其中
$\lambda_m$ 是 Lebesgue 测度；这是变量替换公式的线性情形。

全文中，"可测"指 Borel 可测，$|x|$ 是 $x\in\mathbb R^n$ 或 $\mathbb C^n$ 的 Euclid 范数，
$\|A\|=\sup_{|x|=1}|Ax|$ 是算子范数，$\|A\|_2=(\operatorname{Tr}A^*A)^{1/2}$ 是
Frobenius 范数。

## A.1 测度与积分

**引理 A.1（Dynkin）。** 设 $\mu,\nu$ 是可测空间 $(\Omega,\mathcal F)$ 上的有限测度，
$\mu(\Omega)=\nu(\Omega)$，又设 $\mathcal P\subseteq\mathcal F$ 对有限交封闭，并且作为
$\sigma$-代数生成 $\mathcal F$。若对所有 $P\in\mathcal P$ 有 $\mu(P)=\nu(P)$，则 $\mu=\nu$。

*证明。* 称 $\mathcal D\subseteq\mathcal F$ 为 *Dynkin 系*，如果 $\Omega\in\mathcal D$；只要
$A\subseteq B$ 都属于 $\mathcal D$，就有 $B\setminus A\in\mathcal D$；并且对每个递增序列
$A_k\in\mathcal D$ 有 $\bigcup_kA_k\in\mathcal D$。族 $\{A\in\mathcal F:\mu(A)=\nu(A)\}$ 是 Dynkin 系，因为对有限测度有
$\mu(B\setminus A)=\mu(B)-\mu(A)$，并且测度沿递增并是连续的。设 $d$ 是所有包含 $\mathcal P$ 的 Dynkin 系的交，
它本身也是 Dynkin 系。族 $\mathcal D_1=\{A\in d:\text{对所有 }P\in\mathcal P\text{ 有 }A\cap P\in d\}$ 是 Dynkin 系（对
$A\subseteq A'$，有 $(A'\setminus A)\cap P=(A'\cap P)\setminus(A\cap P)$）并且包含 $\mathcal P$，所以
$\mathcal D_1=d$。于是 $\mathcal D_2=\{A\in d:\text{对所有 }B\in d\text{ 有 }A\cap B\in d\}$ 是包含
$\mathcal P$ 的 Dynkin 系（因为 $\mathcal D_1=d$），所以 $\mathcal D_2=d$：$d$ 对有限交封闭。
因此 $d$ 是 $\sigma$-代数：$\Omega\setminus A\in d$，$A\cup B=\Omega\setminus\big((\Omega\setminus A)\cap(\Omega\setminus B)\big)\in d$，
而可数并是有限并的递增并。所以
$\mathcal F=\sigma(\mathcal P)\subseteq d\subseteq\{\mu=\nu\}$。$\square$

**推论 A.2.** 度量空间 $M$ 上的两个有限 Borel 测度 $\mu,\nu$，若对每个有界连续函数
$f:M\to\mathbb R$ 都有 $\int f\,d\mu=\int f\,d\nu$，则二者相等。

*证明。* 取 $f=1$ 得 $\mu(M)=\nu(M)$。对非空闭集 $K$（$K=\emptyset$ 时无需证明），连续函数
$f_k(x)=\max(0,1-k\operatorname{dist}(x,K))$ 取值于 $[0,1]$，并递减地收敛到 $K$ 的示性函数
（对 $x\notin K$ 有 $\operatorname{dist}(x,K)>0$）。由控制收敛定理得 $\mu(K)=\nu(K)$。闭集族对有限交
封闭并且生成 Borel 集；应用引理 A.1 即可。$\square$

引理 5.1 用的正是这一论证。

**引理 A.3（向量空间上的 Lebesgue 测度）。** 设 $W$ 是 $m$ 维实向量空间。对线性双射
$\varphi:\mathbb R^m\to W$，对 Borel 集 $E\subseteq W$ 令 $\lambda_\varphi(E)=\lambda_m(\varphi^{-1}E)$。

(a) 对另一个线性双射 $\psi$，有 $\lambda_\psi=|\det(\psi^{-1}\varphi)|\,\lambda_\varphi$。因此 $W$ 的零测集，以及
$W$ 中的随机向量是否具有密度，都不依赖于 $\varphi$ 的选取。

(b) 对每个线性双射 $L:W\to W$，有 $\lambda_\varphi(LE)=|\det L|\,\lambda_\varphi(E)$。若 $L$ 对 $W$ 上的某个内积
是正交的，则 $|\det L|=1$，从而 $L$ 保持 $\lambda_\varphi$。

(c) 若 $W$ 中的随机向量 $X$ 关于 $\lambda_\varphi$ 有密度 $f$，则 $LX$ 有密度
$y\mapsto f(L^{-1}y)/|\det L|$。

*证明。* (a) $\lambda_\psi(E)=\lambda_m\big((\psi^{-1}\varphi)\varphi^{-1}E\big)=|\det(\psi^{-1}\varphi)|\,\lambda_m(\varphi^{-1}E)$。
(b) $\lambda_\varphi(LE)=\lambda_m\big((\varphi^{-1}L\varphi)\varphi^{-1}E\big)=|\det(\varphi^{-1}L\varphi)|\,\lambda_\varphi(E)$，并且
$\det(\varphi^{-1}L\varphi)=\det L$ 与基的选取无关。在标准正交基下，正交的 $L$ 的矩阵 $O$ 满足
$O^TO=1$，所以 $(\det L)^2=1$。(c) 由 (b)，$\lambda_\varphi$ 在 $L$ 下的前推测度是
$|\det L|^{-1}\lambda_\varphi$，所以对 Borel 集 $E$，
$\mathbb P(LX\in E)=\int1_E(Lx)f(x)\,d\lambda_\varphi(x)=|\det L|^{-1}\int1_E(y)f(L^{-1}y)\,d\lambda_\varphi(y)$。$\square$

*约定。* 实向量空间上的"Lebesgue 测度"指任意一个 $\lambda_\varphi$；由 (a)，对零测集和密度而言，这一选取
从不影响结果。在 $\mathbb C^n$ 和 $M_{a\times b}(\mathbb C)$ 上，我们以各元素的实部与虚部为坐标；它们对实内积
$\operatorname{Re}\operatorname{Tr}(X^*Y)$ 是标准正交的，而对这个内积，每个映射 $X\mapsto VXW$（$V,W$ 为酉矩阵）都是正交的
（它保持 $\operatorname{Tr}(X^*Y)$）。在 Hermite 矩阵空间 $\operatorname{Herm}_n$ 上，我们使用坐标
$h_{ii}$ 和 $\operatorname{Re}h_{ij},\operatorname{Im}h_{ij}$（$i<j$）；它们对实内积 $\operatorname{Tr}(XY)$ 不是标准正交的，
但由 (b)，对该内积正交的映射 $H\mapsto UHU^*$（$U$ 为酉矩阵）仍然保持 Lebesgue 测度。同样，在实对称矩阵空间
$\operatorname{Sym}_n(\mathbb R)$ 上，我们使用坐标 $h_{ii}$ 和 $h_{ij}$（$i<j$），而对 $\operatorname{Tr}(XY)$ 正交的映射
$H\mapsto OHO^T$（$O\in O(n)$）保持 Lebesgue 测度。

**引理 A.4（向量值函数的积分）。** 设 $V$ 是有限维实或复向量空间，$F:[a,b]\to V$ 连续。在一组基下
逐坐标地（实部与虚部分开）定义 $\int_a^bF(t)\,dt$。则：

(a) 对分划 $a=t_0<\dots<t_N=b$ 及标记点 $\tau_i\in[t_{i-1},t_i]$，当 $\max_i(t_i-t_{i-1})\to0$ 时，Riemann 和
$\sum_iF(\tau_i)(t_i-t_{i-1})$ 收敛到 $\int_a^bF$；特别地，积分不依赖于基的选取；

(b) 对每个映到有限维空间 $V'$ 的（实）线性映射 $L:V\to V'$，有 $L\int_a^bF=\int_a^bLF$；

(c) 对 $V$ 上的每个范数，有 $\|\int_a^bF\|\le\int_a^b\|F(t)\|\,dt$；

(d) 若所有 $F(t)$ 都属于闭凸集 $K$，则 $\frac1{b-a}\int_a^bF\in K$；若 $K$ 还是锥，则
$\int_a^bF\in K$；

(e) $t\mapsto\int_a^tF$ 可微且导数为 $F(t)$；$\int_a^bF=\int_a^cF+\int_c^bF$，并且
$\int_a^bF(s+c)\,ds=\int_{a+c}^{b+c}F(s)\,ds$。

*证明。* $F$ 的每个坐标都是连续实函数，所以 (a) 与 (e) 逐坐标成立；(b) 由 (a) 推出，因为 $L$ 把 Riemann
和映为 Riemann 和并且是连续的。范数是连续的
（$|\|x\|-\|y\||\le\|x-y\|\le\sum_k|x_k-y_k|\,\|e_k\|$），所以 $t\mapsto\|F(t)\|$ 连续，而 (c) 由 (a) 和 Riemann 和的
三角不等式推出。对 (d)，取 $N$ 等分的分划：Riemann 和除以 $b-a$ 是凸组合
$\frac1N\sum_iF(\tau_i)\in K$，而 $K$ 是闭的。$\square$

这为引理 8.4(c)（CP 映射的 Riemann 和）以及命题 8.12 的 Cesàro 平均
$\sigma_\tau=\frac1\tau\int_0^\tau T_s\rho\,ds$ 提供了依据：由 (d)，$\sigma_\tau$ 是一个态；由 (b)、(c) 和 (e)，
$\|T_t\sigma_\tau-\sigma_\tau\|_1\le2t/\tau$。

## A.2 解析函数的零点集

对多重指标 $\alpha\in\mathbb N^n$，记 $|\alpha|=\alpha_1+\dots+\alpha_n$，$y^\alpha=y_1^{\alpha_1}\cdots y_n^{\alpha_n}$。
对 $a\in\mathbb R^n$ 和 $0<\rho\le\infty$，令 $Q_\rho(a)=\{x\in\mathbb R^n:\max_i|x_i-a_i|<\rho\}$ 为开立方体
（$Q_\infty(a)=\mathbb R^n$）。

**定义。** 函数 $f:Q_\rho(a)\to\mathbb R$ 称为 **$Q_\rho(a)$ 上的幂级数**，如果
$f(x)=\sum_\alpha c_\alpha(x-a)^\alpha$，且对每个 $0<r<\rho$ 有 $\sum_\alpha|c_\alpha|r^{|\alpha|}<\infty$。

多项式是 $\mathbb R^n$ 上的幂级数（只有有限多个 $c_\alpha\neq0$）；在多变量情形，本课程只对多项式应用
定理 A.6。（对于在连通开集上仅局部为幂级数的函数，还需要多变量的恒等定理，而本课程不使用它。）
在半径 $r<\rho$ 的闭立方体上，各项以 $|c_\alpha|r^{|\alpha|}$ 为界，所以级数在那里绝对且一致收敛
（Weierstrass M 判别法）：其和与各项的次序无关，并且 $f$ 连续。

**引理 A.5（单变量）。** 设 $g(t)=\sum_kc_k(t-a)^k$ 是 $I=(a-\rho,a+\rho)$ 上的幂级数。

(a) 对 $b\in I$ 和 $|t-b|<\rho_b=\rho-|b-a|$，有 $g(t)=\sum_jd_j(t-b)^j$，其中
$d_j=\sum_{k\ge j}\binom kjc_k(b-a)^{k-j}$，这是 $(b-\rho_b,b+\rho_b)$ 上的幂级数。

(b) 若某个 $d_j\neq0$，则 $g$ 在 $b$ 的某个去心邻域内没有零点。

(c) 若 $g$ 的零点在 $I$ 的某一点处聚集，则 $g$ 在 $I$ 上恒为 $0$。若 $g$ 在 $I$ 上不恒为零，则其零点集
可数，从而是 Lebesgue 零测集。若 $g$ 在 $I$ 上恒为零，则所有 $c_k=0$。

*证明。* (a) 对 $0\le s<\rho_b$，
$\sum_k\sum_{j\le k}|c_k|\binom kjs^j|b-a|^{k-j}=\sum_k|c_k|(s+|b-a|)^k<\infty$，所以把
$(t-a)^k=\big((t-b)+(b-a)\big)^k$ 按二项式展开得到的二重级数在 $|t-b|\le s$ 时绝对收敛，可以按任一次序
求和；先对 $k$ 求和即得 $d_j$，且 $\sum_j|d_j|s^j<\infty$。
(b) 设 $j_0$ 是使 $d_j\ne0$ 的最小的 $j$。则 $g(t)=(t-b)^{j_0}h(t)$，其中
$h(t)=\sum_{j\ge j_0}d_j(t-b)^{j-j_0}$ 在 $b$ 附近连续（M 判别法）且 $h(b)=d_{j_0}\neq0$；所以
$h$ 在 $b$ 附近不为零。(c) 设 $A$ 是由这样的 $b\in I$ 组成的集合：$g$ 在 $b$ 附近恒为零；$A$ 是开集。若
零点在 $b\in I$ 处聚集，则由 (b) 所有 $d_j$ 都为零，再由 (a) 得 $b\in A$。若 $b\in I$ 是 $A$ 中点的极限，
则零点在 $b$ 处聚集，所以 $b\in A$：$A$ 在 $I$ 中是闭的。由于 $I$ 连通，$A$ 要么为空，要么等于 $I$。
若 $g\not\equiv0$，则 $A=\emptyset$，所以 $I$ 中没有点是零点的聚点；于是每个紧区间 $J\subset I$ 只含有限多个
零点（由 Bolzano–Weierstrass 定理，$J$ 的无限子集在 $J$ 中有聚点），而 $I$ 是可数多个这样的 $J$ 的并。
若 $g\equiv0$，则在 $b=a$ 处（此时 $d_j=c_j$）应用 (b) 可知所有 $c_k=0$。$\square$

**定理 A.6（零点集是零测集）。** 设 $f$ 是立方体 $Q\subseteq\mathbb R^n$ 上的幂级数。若 $f$ 在 $Q$ 的一个
Lebesgue 测度为正的子集上为零，则 $f$ 在 $Q$ 上恒为 $0$。等价地：若 $f$ 在 $Q$ 上不恒为零，则
$\{x\in Q:f(x)=0\}$ 是 Lebesgue 零测集。

*证明。* 对 $n$ 归纳；$n=1$ 的情形就是引理 A.5(c)。设 $n\ge2$，写 $x=(x',t)$，其中
$x'\in\mathbb R^{n-1}$，$Q=Q'\times I$，$a=(a',a_n)$，并按 $t-a_n$ 的幂次将各项分组：
$$
f(x',t)=\sum_kh_k(x')(t-a_n)^k,\qquad h_k(x')=\sum_{\alpha'}c_{(\alpha',k)}(x'-a')^{\alpha'} .
$$
这是允许的：对 $x'\in Q'$ 及 $\max_i|x'_i-a_i|\le r<\rho$，
$\sum_k\sum_{\alpha'}|c_{(\alpha',k)}||x'-a'|^{\alpha'}r^k\le\sum_\alpha|c_\alpha|r^{|\alpha|}<\infty$（其中
$|y|^{\alpha'}=\prod_i|y_i|^{\alpha'_i}$）。同一个界说明每个 $h_k$ 都是 $Q'$ 上的幂级数，并且对固定的
$x'\in Q'$，$g_{x'}(t)=f(x',t)$ 是 $I$ 上以 $h_k(x')$ 为系数的幂级数。

令 $Z=\{x\in Q:f(x)=0\}$，由于 $f$ 连续，它在 $Q$ 中是闭的；假设 $\lambda_n(Z)>0$。由 Tonelli 定理，
$\lambda_n(Z)=\int_{Q'}\lambda_1(Z_{x'})\,dx'$，其中 $Z_{x'}=\{t\in I:f(x',t)=0\}$，所以集合
$F=\{x'\in Q':\lambda_1(Z_{x'})>0\}$ 的 $(n-1)$ 维测度为正。对 $x'\in F$，引理 A.5(c) 给出 $I$ 上
$g_{x'}\equiv0$，进而对所有 $k$ 有 $h_k(x')=0$。所以每个 $h_k$ 在 $F$ 上为零，由归纳假设，$h_k$ 在
$Q'$ 上恒为 $0$。因此在 $Q$ 上 $f(x',t)=\sum_kh_k(x')(t-a_n)^k=0$。$\square$

**推论 A.7.** (a) $\mathbb R^n$ 上不恒为零的实多项式 $p$ 只在一个 Lebesgue 零测集上为零。系数不全为零的
多项式不恒为零。

(b) 设 $W$ 是有限维实向量空间，$p:W\to\mathbb R$ 是关于某组基的坐标的多项式（从而关于每组基都是多项式），
且不恒为零。若 $W$ 中的随机向量 $X$ 有密度，则 $\mathbb P(p(X)=0)=0$。

(c) 真线性子空间 $U\subsetneq W$ 是零测集，具有密度的随机向量几乎必然落在 $U$ 之外。特别地，若
$E\neq0$ 是 $\mathbb C^n$ 的子空间，而 $v$ 在 $\mathbb C^n$ 上有密度（例如高斯初始向量），则 $v$ 几乎必然
在 $E$ 中有非零分量。

(d) 设 $X$ 与 $Y$ 是有限维实空间 $W_1$ 与 $W_2$ 中相互独立的随机向量，$Y$ 有密度，
$N\subseteq W_1\times W_2$ 是 Borel 集，且对所有 $x$，其截面 $N_x=\{y:(x,y)\in N\}$ 都是零测集。则
$\mathbb P((X,Y)\in N)=0$。

*证明。* (a) 第一个结论是 $Q=\mathbb R^n$ 时的定理 A.6。对第二个结论，对 $n$ 归纳：写
$p(x',t)=\sum_kp_k(x')t^k$；某个 $p_k$ 有非零系数，所以由归纳假设，对某个 $x'_0$ 有 $p_k(x'_0)\neq0$
（$n=1$ 时 $p_k$ 是常数，没有 $x'$）；于是 $t\mapsto p(x'_0,t)$ 是非零的单变量多项式，只有有限多个根。
(b) 基变换是线性的，所以 $p$ 仍是多项式，而由 (a) 和引理 A.3(a)，$\{p=0\}$ 是零测集；
于是 $\mathbb P(p(X)=0)=\int_{\{p=0\}}f\,d\lambda=0$。(c) 取一个在 $U$ 上为零的非零线性泛函 $\ell$，
对 $p=\ell$ 应用 (b)。对于例子，$v\mapsto P_Ev$（$P_E$ 为正交投影）是非零的实线性映射，而
$\{P_Ev=0\}$ 是它的核，一个真子空间。(d) $(X,Y)$ 的分布是两个分布的乘积，所以由
Tonelli 定理，$\mathbb P((X,Y)\in N)=\int\mathbb P(Y\in N_x)\,d\mathbb P_X(x)=0$。$\square$

本课程把推论 A.7 应用于 Ginibre 矩阵 $Z$ 的 $\det Z$（$|\det Z|^2$ 是各元素实部与虚部的多项式，在
$Z=1$ 处等于 $1$；定理 5.3，1.6 节），命题 4.7(a) 的 Gram 行列式，高斯矩阵的秩（第 12 章习题 4），
以及 Lanczos 方法的初始向量（第 9 章）。(d) 部分就是命题 4.7(b) 中"由于 $Y$ 与 $K$ 独立，Tonelli 定理
给出概率 $0$"这一步。

```python
import math
import numpy as np

rng = np.random.default_rng(0)


def mc(samples, exact, label, kind="exact"):
    """Monte Carlo mean with standard error; the reference value must lie within 5 standard errors."""
    samples = np.asarray(samples, dtype=float)
    mean, err = samples.mean(), samples.std() / math.sqrt(samples.size)
    print(f"{label}: {mean:.4f} +- {err:.4f}   ({kind} {exact:.4f})")
    assert abs(mean - exact) <= 5 * err + 1e-12, label


# Corollary A.7: det Z = 0 has probability 0 for a real Gaussian 3 x 3 matrix, so by continuity of the measure
# P(|det Z| < eps) -> 0 as eps -> 0. Numerical observation (not a theorem of this appendix): the fractions
# below shrink by about a factor 10 per factor 10 in eps.
d = np.abs(np.linalg.det(rng.standard_normal((200_000, 3, 3))))
fractions = [(d < eps).mean() for eps in (1e-1, 1e-2, 1e-3)]
print("P(|det Z| < eps) for eps = 0.1, 0.01, 0.001:", fractions)
assert fractions[1] < 0.2 * fractions[0] and fractions[2] < 0.2 * fractions[1]
```

## A.3 Hermite 矩阵的特征值

对 $A\in\operatorname{Herm}_n$，令 $\lambda_1(A)\le\dots\le\lambda_n(A)$ 为其按重数计的特征值，
$v_1,\dots,v_n$ 为由特征向量组成的标准正交基（定理 1.3，推论 1.4）。对单位向量 $x$，
$\langle x,Ax\rangle=\sum_i\lambda_i|\langle v_i,x\rangle|^2$ 且 $\sum_i|\langle v_i,x\rangle|^2=1$，所以
$\lambda_1(A)\le\langle x,Ax\rangle\le\lambda_n(A)$；并且 $\|A\|=\max_i|\lambda_i(A)|$（推论 1.4）。

**定理 A.8（Courant–Fischer）。** 对 $A\in\operatorname{Herm}_n$ 和 $1\le j\le n$，
$$
\lambda_j(A)=\min_{\dim U=j}\ \max_{x\in U,\,|x|=1}\langle x,Ax\rangle=\max_{\dim U=n-j+1}\ \min_{x\in U,\,|x|=1}\langle x,Ax\rangle ,
$$
其中 $U$ 取遍 $\mathbb C^n$ 中给定维数的子空间。（内层的最大值与最小值存在，因为 $U$ 的单位球面是紧的；
外层的最大值与最小值可以达到。）

*证明。* 设 $\dim U=j$，$W=\operatorname{span}(v_j,\dots,v_n)$。由于 $\dim U+\dim W=n+1>n$，存在单位向量
$x\in U\cap W$，并且 $\langle x,Ax\rangle=\sum_{i\ge j}\lambda_i|\langle v_i,x\rangle|^2\ge\lambda_j$。所以在 $U$ 上的最大值
$\ge\lambda_j$。对 $U=\operatorname{span}(v_1,\dots,v_j)$，有 $\langle x,Ax\rangle=\sum_{i\le j}\lambda_i|\langle v_i,x\rangle|^2\le\lambda_j$，
且在 $x=v_j$ 处取等号；所以对 $U$ 的最小值可以达到并且等于 $\lambda_j$。对第二个公式，
把第一个公式应用于 $-A$，其特征值为 $\lambda_j(-A)=-\lambda_{n+1-j}(A)$：
$-\lambda_j(A)=\lambda_{n+1-j}(-A)=\min_{\dim U=n+1-j}\max_x(-\langle x,Ax\rangle)=-\max_{\dim U=n+1-j}\min_x\langle x,Ax\rangle$。$\square$

**推论 A.9（压缩）。** 设 $Q\in M_{n\times m}(\mathbb C)$ 满足 $Q^*Q=1_m$，$A\in\operatorname{Herm}_n$。则
$\lambda_j(Q^*AQ)=\min\max\langle x,Ax\rangle$，其中 $U$ 取遍 $\operatorname{ran}Q$ 中的 $j$ 维子空间，$x$ 取遍 $U$ 中的单位向量，并且
$$
\lambda_j(A)\le\lambda_j(Q^*AQ)\le\lambda_{j+n-m}(A)\qquad(1\le j\le m).
$$

*证明。* $U'\mapsto QU'$ 是从 $\mathbb C^m$ 的 $j$ 维子空间到 $\operatorname{ran}Q$ 的 $j$ 维子空间的双射，
$|Qy|=|y|$，且 $\langle y,Q^*AQy\rangle=\langle Qy,AQy\rangle$；所以对 $Q^*AQ$ 应用定理 A.8 即得该公式。
在更大的族——所有 $j$ 维子空间 $U\subseteq\mathbb C^n$——上取最小值得到 $\lambda_j(A)$，因此它更小。
同样，$\lambda_j(Q^*AQ)$ 的 max–min 公式取遍 $(m-j+1)$ 维的 $U\subseteq\operatorname{ran}Q$，而在 $\mathbb C^n$ 中
所有这一维数的 $U$ 上，它给出更大的数 $\lambda_{n-(m-j+1)+1}(A)=\lambda_{j+n-m}(A)$。$\square$

**定理 A.10（Weyl）。** 对 $A,E\in\operatorname{Herm}_n$ 和所有 $j$，
$$
\lambda_j(A)+\lambda_1(E)\le\lambda_j(A+E)\le\lambda_j(A)+\lambda_n(E),
\qquad\text{hence}\qquad
|\lambda_j(A+E)-\lambda_j(A)|\le\|E\|\le\|E\|_2 .
$$

*证明。* 对单位向量 $x$，$\langle x,(A+E)x\rangle\le\langle x,Ax\rangle+\lambda_n(E)$。先对 $j$ 维子空间 $U$ 的单位向量
取最大值，再对 $U$ 取最小值（定理 A.8）：$\lambda_j(A+E)\le\lambda_j(A)+\lambda_n(E)$。将此应用于
$A+E$ 与 $-E$，得到 $\lambda_j(A)\le\lambda_j(A+E)+\lambda_n(-E)=\lambda_j(A+E)-\lambda_1(E)$。因此
$|\lambda_j(A+E)-\lambda_j(A)|\le\max(|\lambda_1(E)|,|\lambda_n(E)|)=\|E\|$，并且 $\|E\|^2=\max_i\lambda_i(E)^2\le\sum_i\lambda_i(E)^2=\|E\|_2^2$。$\square$

**推论 A.11.** 每个 $\lambda_j:\operatorname{Herm}_n\to\mathbb R$ 关于 $\|\cdot\|$ 是 $1$-Lipschitz 的，从而连续且
Borel 可测。对随机 Hermite 矩阵 $H$，$(\lambda_1(H),\dots,\lambda_n(H))$ 是随机向量，其间距与间距比（在其分母
不为零的 Borel 集上）也是随机向量；相互独立的随机 Hermite 矩阵的特征值向量相互独立。具有
$n$ 个互不相同特征值的矩阵构成开集。

*证明。* 取 $E=B-A$ 应用定理 A.10。相互独立的随机变量的函数相互独立。若对所有 $j$ 有
$\lambda_j(A)<\lambda_{j+1}(A)$，则当 $\|B-A\|<\frac12\min_j(\lambda_{j+1}(A)-\lambda_j(A))$ 时这些严格不等式保持成立。$\square$

**引理 A.12（Lagrange 公式）。** 设 $A\in\operatorname{Herm}_n$ 有互不相同的特征值 $\mu_1<\dots<\mu_r$ 与
谱投影 $P_1,\dots,P_r$（到各特征空间上的正交投影，从而 $A=\sum_i\mu_iP_i$，
$P_iP_l=\delta_{il}P_i$，$\sum_iP_i=1$）。则对每个多项式 $p$ 有 $p(A)=\sum_ip(\mu_i)P_i$，并且
$$
P_i=\prod_{l\neq i}\frac{A-\mu_l}{\mu_i-\mu_l} .
$$

*证明。* 对 $k$ 归纳（利用 $P_iP_l=\delta_{il}P_i$）得 $A^k=\sum_i\mu_i^kP_i$，再由线性得 $p(A)=\sum_ip(\mu_i)P_i$。
多项式 $\ell_i(t)=\prod_{l\neq i}(t-\mu_l)/(\mu_i-\mu_l)$ 满足 $\ell_i(\mu_l)=\delta_{il}$，所以 $\ell_i(A)=P_i$。$\square$

**命题 A.13（特征值与谱投影的连续性）。** 设 $S$ 是有限维实向量空间的子集，$A:S\to\operatorname{Herm}_n$
连续，并假设对每个 $x\in S$，$A(x)$ 恰有 $r$ 个互不相同的特征值 $\mu_1(x)<\dots<\mu_r(x)$。则它们的重数在
$S$ 上局部常值，而 $\mu_i$ 以及 $A(x)$ 的谱投影 $P_i(x)$ 在 $S$ 上连续。特别地，在具有单谱的 Hermite 矩阵
构成的开集上，特征值与秩一谱投影是连续的。

*证明。* 若 $r=1$，则 $A(x)=\mu_1(x)1$，其中 $\mu_1=\operatorname{Tr}A/n$，$P_1=1$。设 $r\ge2$，固定 $x_0\in S$，令
$d=\min_i(\mu_{i+1}(x_0)-\mu_i(x_0))>0$。令 $J_i$ 为满足 $\lambda_j(A(x_0))=\mu_i(x_0)$ 的指标 $j$ 的集合，它是一段
连续整数。由 $A$ 的连续性，存在 $x_0$ 在 $S$ 中的邻域 $N$，使得对 $x\in N$ 有 $\|A(x)-A(x_0)\|<d/2$；
由定理 A.10，对 $j\in J_i$，$\lambda_j(A(x))$ 落在区间 $I_i=(\mu_i(x_0)-d/2,\mu_i(x_0)+d/2)$ 中，而这
$r$ 个区间两两不交。所以 $A(x)$ 在每个 $I_i$ 中至少有一个特征值，而它恰有 $r$ 个互不相同的特征值，
所以在每个 $I_i$ 中恰有一个：$\lambda_j(A(x))$ 对所有 $j\in J_i$ 都相同。因此在 $N$ 上重数为常数 $|J_i|$，并且
$\mu_i(x)=\lambda_{\max J_i}(A(x))$ 连续（推论 A.11）。由引理 A.12，
$P_i(x)=\prod_{l\neq i}(A(x)-\mu_l(x))/(\mu_i(x)-\mu_l(x))$ 在 $N$ 上连续，因为分母始终远离 $0$。
最后一个结论是 $r=n$、$S$ 为推论 A.11 中那个开集的情形。$\square$

在命题 4.7 中，$S$ 是 Gram 行列式 $\delta$ 不为零的集合，在其上 $K$ 恰有 $N$ 个互不相同的特征值；所以
$K\mapsto\Pi_i(K)$ 在那里连续，从而命题 4.7(b) 中的事件是 Borel 集。定理 5.15
用到推论 A.11，定理 9.14 用到推论 A.9。

```python
def herm(n, batch):
    G = rng.standard_normal((batch, n, n)) + 1j * rng.standard_normal((batch, n, n))
    return (G + G.conj().transpose(0, 2, 1)) / 2


n = 6
A, E = herm(n, 2000), 0.1 * herm(n, 2000)
lam_A, lam_AE, lam_E = np.linalg.eigvalsh(A), np.linalg.eigvalsh(A + E), np.linalg.eigvalsh(E)   # ascending
normE = np.linalg.norm(E, ord=2, axis=(1, 2))                                    # operator norm
assert (np.abs(lam_AE - lam_A).max(axis=1) <= normE + 1e-12).all()               # Theorem A.10
assert (lam_AE >= lam_A + lam_E[:, :1] - 1e-12).all() and (lam_AE <= lam_A + lam_E[:, -1:] + 1e-12).all()

a0, j = A[0], 3                          # Theorem A.8: every 3-dim subspace has max <x, Ax> >= lambda_3
for _ in range(200):
    Q, _ = np.linalg.qr(rng.standard_normal((n, j)) + 1j * rng.standard_normal((n, j)))
    assert np.linalg.eigvalsh(Q.conj().T @ a0 @ Q)[-1] >= lam_A[0, j - 1] - 1e-12
m = 4                                    # Corollary A.9: interlacing for a compression to 4 dimensions
Q, _ = np.linalg.qr(rng.standard_normal((n, m)) + 1j * rng.standard_normal((n, m)))
theta = np.linalg.eigvalsh(Q.conj().T @ a0 @ Q)
assert (theta >= lam_A[0, :m] - 1e-12).all() and (theta <= lam_A[0, n - m:] + 1e-12).all()

U, _ = np.linalg.qr(rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)))
mu, mult = np.array([-1.0, 0.5, 2.0]), [2, 3, 1]          # Lemma A.12 with multiplicities 2, 3, 1
spec = np.repeat(mu, mult)
A1 = U @ np.diag(spec) @ U.conj().T
for i, m_i in enumerate(mult):
    P = np.eye(n, dtype=complex)
    for l in range(len(mu)):
        if l != i:
            P = P @ (A1 - mu[l] * np.eye(n)) / (mu[i] - mu[l])
    Vi = U[:, spec == mu[i]]
    assert np.abs(P - Vi @ Vi.conj().T).max() < 1e-10 and abs(np.trace(P).real - m_i) < 1e-10
```

## A.4 高斯向量

$\mathbb R^m$ 中的**标准高斯向量**是坐标相互独立且服从 $N(0,1)$ 的 $g=(g_1,\dots,g_m)$。它的分布是各坐标
分布的乘积，所以由 Tonelli 定理，它有密度
$\varphi_m(x)=(2\pi)^{-m/2}e^{-|x|^2/2}$。

**命题 A.14（正交不变性）。** (a) 对每个正交矩阵 $O\in O(m)$，$Og$ 与 $g$ 同分布。

(b) 设 $W$ 是 $m$ 维实内积空间。对标准正交基 $(w_i)$，$\sum_ig_iw_i$ 的分布不依赖于基的选取；称这样的向量为
**$W$ 的标准高斯向量**。它在每组标准正交基下的坐标都是相互独立的 $N(0,1)$，它的分布在 $W$ 的正交映射下不变，
并且若 $W=W_1\oplus W_2$ 是正交直和，则它的两个分量是 $W_1$ 与 $W_2$ 的相互独立的标准高斯向量。

(c) 对 $v\in\mathbb R^m$，$\langle v,g\rangle$ 服从 $N(0,|v|^2)$；特别地，$\sum_ia_ig_i\sim N(0,\sum_ia_i^2)$。

*证明。* (a) 由引理 A.3(c)，$Og$ 有密度 $\varphi_m(O^{-1}y)/|\det O|=\varphi_m(y)$。(b) 若 $w'_i=\sum_jO_{ji}w_j$
是另一组标准正交基，则 $O$ 是正交的，且 $\sum_ig_iw'_i=\sum_j(Og)_jw_j$，由 (a)，它与
$\sum_jg_jw_j$ 同分布。$X=\sum_ig_iw_i$ 在另一组标准正交基 $(w'_i)$ 下的坐标是 $O^Tg$，由 (a) 仍是
标准的。正交映射 $L$ 把 $(w_i)$ 映为标准正交基，并且 $LX=\sum_ig_iLw_i$。对最后一个结论，取与
$W_1\oplus W_2$ 相适应的 $(w_i)$。(c) 当 $v=0$ 时两边都是 $0$ 处的点质量。当 $v\neq0$ 时，
把 $v/|v|$ 扩充为标准正交基；由 (b)，坐标 $\langle v/|v|,g\rangle$ 服从 $N(0,1)$，所以
$\langle v,g\rangle=|v|\,\langle v/|v|,g\rangle$ 有密度 $(2\pi|v|^2)^{-1/2}e^{-y^2/2|v|^2}$（$m=1$ 时的引理 A.3(c)）。$\square$

**定义。** 有限维实向量空间 $W$ 中的**中心化高斯向量**是 $X=Ag$，其中 $A:\mathbb R^M\to W$ 是线性映射，
$g$ 是 $\mathbb R^M$ 的标准高斯向量。在坐标下（$W=\mathbb R^m$），它的**协方差**是 $C=\mathbb E[XX^T]=AA^T$。
若 $A$ 是满射，等价地 $C$ 正定（$\langle x,Cx\rangle=|A^Tx|^2$），则称它是**非退化的**。

本课程中每个由相互独立的高斯元素线性构造出的随机向量——复高斯向量、Ginibre 矩阵、GOE 与 GUE 矩阵，
以及它们在 twirl 下的像——都是中心化高斯向量。

**定理 A.15（线性像与密度）。** 设 $X=Ag$ 是 $W$ 中的中心化高斯向量。

(a) 对每个线性映射 $L:W\to W'$，$LX$ 是中心化高斯向量；若 $X$ 非退化且 $L$ 是满射，则 $LX$ 非退化。

(b) 若 $X$ 非退化，则它关于 $W$ 上的 Lebesgue 测度有密度；在坐标
$W=\mathbb R^m$ 下，
$$
f_X(x)=(2\pi)^{-m/2}(\det C)^{-1/2}\,e^{-\langle x,C^{-1}x\rangle/2}.
$$
若 $X$ 退化，则它落在真子空间 $\operatorname{ran}A$ 中，这是零测集（推论 A.7(c)），因而 $X$ 没有密度。

(c) $X$ 的分布只依赖于它的协方差。

*证明。* (a) $LX=(LA)g$，并且若 $A$ 与 $L$ 都是满射，则 $LA$ 是满射。(b) 令 $K=\ker A$，于是
$\mathbb R^M=K^\perp\oplus K$。由命题 A.14(b)，$g=g_1+g_2$，其中 $g_1$ 是 $K^\perp$ 的标准高斯向量，并且
$X=Ag_1=Bg_1$，其中 $B=A|_{K^\perp}$ 是单射，又因 $A$ 是满射，$B$ 是双射
$K^\perp\to W$；所以 $\dim K^\perp=m$。在由 $K^\perp$ 的一组标准正交基再补上 $K$ 的一组标准正交基所得的基下，
$A$ 的矩阵是 $(B\ 0)$，所以 $C=AA^T=BB^T$。$g_1$ 的坐标有密度 $\varphi_m$，由引理 A.3(c)，$X$ 有密度
$\varphi_m(B^{-1}x)/|\det B|$，其中 $|B^{-1}x|^2=\langle x,(BB^T)^{-1}x\rangle$，$|\det B|=(\det C)^{1/2}$。若
$X$ 退化却有密度 $f$，则 $1=\mathbb P(X\in\operatorname{ran}A)=\int_{\operatorname{ran}A}f\,d\lambda=0$。
(c) $\ker C=\ker A^T$（若 $AA^Tx=0$，则 $|A^Tx|^2=0$），所以 $\operatorname{ran}C=(\ker C)^\perp=(\ker A^T)^\perp=\operatorname{ran}A$
包含 $X$。令 $d=\operatorname{rank}C$（若 $d=0$，则 $X=0$），取 $\operatorname{ran}C$ 的一组标准正交基
$p_1,\dots,p_d$（它只依赖于 $C$），并令 $P$ 为以 $p_k^T$ 为行的 $d\times m$ 矩阵，于是 $PP^T=1_d$，
而 $P^TP$ 是到 $\operatorname{ran}C$ 上的正交投影。于是 $PX=(PA)g$，其中 $PA$ 映满 $\mathbb R^d$，这是协方差为
$PCP^T$ 的非退化中心化高斯向量；由 (b)，它的分布只依赖于 $C$，从而
$X=P^T(PX)$ 的分布也只依赖于 $C$。$\square$

例如，在第 4 章中，twirl $P$ 是 $M_n(\mathbb C)$ 到 $\mathcal A'$ 上的线性投影，所以
$Z\mapsto P(Z)$ 与 $Z\mapsto P(Z+Z^*)$ 分别是到 $\mathcal A'$ 和 $\mathcal A'_h$ 上的实线性满射（对
$H\in\mathcal A'_h$，$P(\frac H2+\frac H2)=H$）；对 Ginibre 矩阵 $Z$（由命题 A.16，它是非退化的），由定理 A.15，
$Y=P(Z)$ 与 $K=P(Z+Z^*)$ 在这些空间上有密度。$\mathcal A'$ 的一组基 $(B_k)$ 的高斯组合
$\sum_k\xi_kB_k$ 是标准高斯向量 $\xi$ 在一个线性双射下的像。

**命题 A.16（复高斯向量，Ginibre 矩阵）。** (a) 设 $z\in\mathbb C^n$ 的坐标相互独立，且每个坐标的实部与虚部
相互独立并服从 $N(0,\frac12)$（**标准复高斯向量**）。则 $\sqrt2z$ 是带内积
$\operatorname{Re}\langle u,w\rangle$ 的 $\mathbb C^n\cong\mathbb R^{2n}$ 的标准高斯向量；$z$ 有密度 $\pi^{-n}e^{-|z|^2}$；对每个
$U\in U(n)$，$Uz$ 与 $z$ 同分布，所以在任意标准正交基下的坐标 $\langle w_i,z\rangle$ 仍是相互独立的标准复高斯变量；
并且 $\mathbb E[z_i\bar z_j]=\delta_{ij}$，$\mathbb E[z_iz_j]=0$。

(b) $a\times b$ Ginibre 矩阵 $G$（元素为相互独立的标准复高斯变量）在 $M_{a\times b}(\mathbb C)$ 上有密度
$\pi^{-ab}e^{-\|G\|_2^2}$，并且对所有 $V\in U(a)$，$W\in U(b)$，$VGW$ 与 $G$ 同分布。

(c) 元素相互独立且服从 $N(0,1)$ 的 $a\times b$ 实矩阵 $G$ 有密度 $(2\pi)^{-ab/2}e^{-\|G\|_2^2/2}$，并且对所有
$O\in O(a)$，$O'\in O(b)$，$OGO'$ 与 $G$ 同分布。

*证明。* (a) $2n$ 个数 $\sqrt2\operatorname{Re}z_i,\sqrt2\operatorname{Im}z_i$ 相互独立且服从 $N(0,1)$，它们是
$\sqrt2z$ 在标准正交基 $(e_i,ie_i)$ 下的坐标。由 $L=2^{-1/2}$ 时的引理 A.3(c)，$z$ 有密度
$(2\pi)^{-n}e^{-|\sqrt2z|^2/2}2^{n}=\pi^{-n}e^{-|z|^2}$。酉矩阵 $U$ 是实线性的且保持
$\operatorname{Re}\langle u,w\rangle$，因而是正交的；应用命题 A.14(b)。在基 $(w_i)$ 下的坐标是
$W^*z$，其中 $W$ 是以 $w_i$ 为列的酉矩阵。最后，记 $z_j=x_j+iy_j$，则
$\mathbb E[z_i\bar z_j]=\mathbb E[x_ix_j+y_iy_j]+i\,\mathbb E[y_ix_j-x_iy_j]=\delta_{ij}$，$\mathbb E[z_iz_j]=\mathbb E[x_ix_j-y_iy_j]+i\,\mathbb E[x_iy_j+y_ix_j]=0$。
(b) $M_{a\times b}(\mathbb C)\cong\mathbb C^{ab}$，且 $|G|=\|G\|_2$，而 $X\mapsto VXW$ 对
$\operatorname{Tr}(X^*Y)$ 是酉的（引理 A.3 之后的约定）；应用 (a)。(c) 同理，用 $\mathbb R^{ab}$ 和命题 A.14。$\square$

**命题 A.17（GOE 与 GUE）。** 设 $H=(G+G^T)/2$，其中 $G$ 是元素相互独立且服从 $N(0,1)$ 的 $n\times n$ 实矩阵
（GOE），或 $H=(G+G^*)/\sqrt2$，其中 $G$ 是 Ginibre 矩阵（GUE）。则元素 $h_{ij}$（$i\le j$）相互独立，
$h_{ii}\sim N(0,1)$，而对 $i<j$，$h_{ij}\sim N(0,\frac12)$（GOE），或者 $\operatorname{Re}h_{ij},\operatorname{Im}h_{ij}$
相互独立且服从 $N(0,\frac12)$（GUE）。$H$ 在引理 A.3 之后所定的坐标下的密度为
$$
(2\pi)^{-n/2}\pi^{-n(n-1)/4}e^{-\operatorname{Tr}H^2/2}\ \ \text{(GOE)},\qquad
(2\pi)^{-n/2}\pi^{-n(n-1)/2}e^{-\operatorname{Tr}H^2/2}\ \ \text{(GUE)},
$$
并且 $OHO^T$（$O\in O(n)$）或 $UHU^*$（$U\in U(n)$）与 $H$ 同分布。

*证明。* GOE：$G$ 是带内积 $\operatorname{Tr}(X^TY)$ 的 $M_n(\mathbb R)$ 的标准高斯向量，矩阵
$e_{ii}$，$(e_{ij}\pm e_{ji})/\sqrt2$（$i<j$）构成标准正交基。由命题 A.14(b)，坐标
$G_{ii}$，$(G_{ij}+G_{ji})/\sqrt2$，$(G_{ij}-G_{ji})/\sqrt2$ 相互独立且服从 $N(0,1)$。而 $h_{ii}=G_{ii}$，
$h_{ij}=2^{-1/2}\cdot(G_{ij}+G_{ji})/\sqrt2$。GUE：$\sqrt2G$ 是带内积
$\operatorname{Re}\operatorname{Tr}(X^*Y)$ 的 $M_n(\mathbb C)$ 的标准高斯向量（命题 A.16），标准正交基为 $e_{ii}$，$ie_{ii}$，
$(e_{ij}\pm e_{ji})/\sqrt2$，$i(e_{ij}\pm e_{ji})/\sqrt2$；$\sqrt2G$ 的坐标 $\sqrt2\operatorname{Re}G_{ii}$，
$\sqrt2\operatorname{Im}G_{ii}$，$\operatorname{Re}(G_{ij}\pm G_{ji})$，$\operatorname{Im}(G_{ij}\pm G_{ji})$ 相互独立且服从 $N(0,1)$。而
$h_{ii}=\sqrt2\operatorname{Re}G_{ii}$，$\operatorname{Re}h_{ij}=2^{-1/2}\operatorname{Re}(G_{ij}+G_{ji})$，
$\operatorname{Im}h_{ij}=2^{-1/2}\operatorname{Im}(G_{ij}-G_{ji})$。在两种情形下，密度都是一维密度的乘积：
$(2\pi)^{-1/2}e^{-h_{ii}^2/2}$，以及对每个方差为 $\frac12$ 的实坐标 $t$ 的 $\pi^{-1/2}e^{-t^2}$；并且
$\operatorname{Tr}H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2$。映射 $H\mapsto UHU^*$ 是线性、可逆的，并且对
$\operatorname{Tr}(XY)$ 正交，所以由引理 A.3(b)、(c)，$UHU^*$ 有密度 $f(U^*HU)=f(H)$，因为 $f$ 只依赖于
$\operatorname{Tr}H^2$；对 $O$ 同理。$\square$

```python
Amat = rng.standard_normal((2, 4))                       # Theorem A.15: X = A g in R^2, A of rank 2
C = Amat @ Amat.T
X = rng.standard_normal((400_000, 4)) @ Amat.T
h = 1 / 400                                              # midpoint rule for the density on [0, 1]^2
grid = (np.arange(400) + 0.5) * h
xx = np.stack(np.meshgrid(grid, grid, indexing='ij'), -1).reshape(-1, 2)
dens = np.exp(-0.5 * np.einsum('ni,ij,nj->n', xx, np.linalg.inv(C), xx)) / (2 * math.pi * math.sqrt(np.linalg.det(C)))
mc(((X >= 0) & (X <= 1)).all(axis=1), dens.sum() * h * h, "P(X in [0,1]^2) against the integral of f_X",
   kind="midpoint rule")

nn, B = 4, 100_000                                       # Proposition A.17
G = rng.standard_normal((B, nn, nn))
goe = (G + G.transpose(0, 2, 1)) / 2
Gc = (rng.standard_normal((B, nn, nn)) + 1j * rng.standard_normal((B, nn, nn))) / math.sqrt(2)
gue = (Gc + Gc.conj().transpose(0, 2, 1)) / math.sqrt(2)
mc(np.einsum('bij,bji->b', goe, goe), nn * (nn + 1) / 2, "E Tr H^2, GOE")      # n + 2 * n(n-1)/2 * 1/2
mc(np.einsum('bij,bji->b', gue, gue).real, nn ** 2, "E Tr H^2, GUE")           # n + 2 * n(n-1)/2 * 1
O, _ = np.linalg.qr(rng.standard_normal((nn, nn)))
rot = O @ goe @ O.T                                                             # has the law of the GOE
mc(rot[:, 0, 0] ** 4, 3.0, "E (O H O^T)_11^4")
mc(rot[:, 0, 1] ** 2, 0.5, "E (O H O^T)_12^2")
```

## A.5 高斯向量的方向

令 $S^{n-1}=\{x\in\mathbb R^n:|x|=1\}$，$V_n=\lambda_n(\{|x|\le1\})$，并对 Borel 集 $B\subseteq S^{n-1}$ 令
$C(B)=\{x\in\mathbb R^n:0<|x|\le1,\ x/|x|\in B\}$，这是 Borel 集（$x\mapsto(|x|,x/|x|)$ 在 $\mathbb R^n\setminus0$ 上连续）。

**定理 A.18（极坐标）。** 对 Borel 集 $A\subseteq(0,\infty)$ 和 $B\subseteq S^{n-1}$，
$$
\lambda_n\big(\{x\neq0:|x|\in A,\ x/|x|\in B\}\big)=n\,\lambda_n(C(B))\int_Ar^{n-1}\,dr .
$$
因此，对 Borel 函数 $h:(0,\infty)\to[0,\infty]$ 和 $\phi:S^{n-1}\to[0,\infty]$，记 $\mu(B)=\lambda_n(C(B))$，有
$$
\int_{\mathbb R^n}h(|x|)\,\phi(x/|x|)\,dx=n\int_0^\infty h(r)r^{n-1}\,dr\int_{S^{n-1}}\phi\,d\mu .
$$

*证明。* 固定 $B$，令 $\nu_1(A)$ 和 $\nu_2(A)$ 为等式两边，它们是 $(0,\infty)$ 上的 Borel 测度。伸缩
$x\mapsto rx$ 的行列式为 $r^n$，并且把 $C(B)$ 映满 $\{0<|x|\le r,\ x/|x|\in B\}$，因为 $rx/|rx|=x/|x|$；所以
$\nu_1((0,r])=r^n\lambda_n(C(B))=\nu_2((0,r])$。对 $R>0$，$\nu_1,\nu_2$ 在 $(0,R]$ 上的限制是有限的
（至多为 $\lambda_n(\{|x|\le R\})$），并且在区间 $(0,r]$（$r\le R$）上相等；这些区间对交封闭，
包含 $(0,R]$，并且生成 $(0,R]$ 的 Borel 集（它们给出 $(s,r]=(0,r]\setminus(0,s]$，而 $(0,R]$ 的每个相对开
子集都是可数多个这样的区间的并）。由引理 A.1，$\nu_1=\nu_2$ 在 $(0,R]$ 上成立；令 $R\to\infty$。
积分公式对 $h=1_A$，$\phi=1_B$ 成立（即刚证明的恒等式；$\{0\}$ 是零测集），由线性性对非负简单函数
$h,\phi$ 成立，再沿简单函数 $h_k\uparrow h$，$\phi_k\uparrow\phi$ 用单调收敛定理，对一般的 $h,\phi$ 成立。$\square$

**推论 A.19（径向积分）。** 对 $(0,\infty)$ 上的 Borel 函数 $h\ge0$，
$\int_{\mathbb R^n}h(|x|)\,dx=s_n\int_0^\infty h(r)r^{n-1}\,dr$，其中 $s_n=nV_n=(2\pi)^{n/2}/I_n$，
$I_n=\int_0^\infty r^{n-1}e^{-r^2/2}dr$。这里 $I_1=\sqrt{\pi/2}$，$I_2=1$，$I_{n+2}=nI_n$，所以 $s_1=2$，$s_2=2\pi$，$s_3=4\pi$，
$s_4=2\pi^2$；一般地，$s_n=2\pi^{n/2}/\Gamma(n/2)$。

*证明。* 该公式是 $\phi=1$ 时的定理 A.18，$\mu(S^{n-1})=V_n$。对 $h(r)=e^{-r^2/2}$，左边是
$\big(\int_{\mathbb R}e^{-t^2/2}dt\big)^n=(2\pi)^{n/2}$（Fubini），由此得到 $s_n$。$I_1=\frac12\sqrt{2\pi}$，$I_2=1$，
分部积分给出 $I_{n+2}=\big[-r^ne^{-r^2/2}\big]_0^\infty+n\int_0^\infty r^{n-1}e^{-r^2/2}dr=nI_n$；代换
$r=\sqrt{2u}$ 给出 $I_n=2^{n/2-1}\Gamma(n/2)$。$\square$

当 $n=2$ 和 $n=3$ 时，这就是命题 5.13 中的极坐标与球坐标。

**定义。** $S^{n-1}$ 上的**均匀分布**是 $\sigma_{n-1}(B)=\lambda_n(C(B))/V_n$。在
$\mathbb C^N\cong\mathbb R^{2N}$ 的单位球面上，它是 $\sigma_{2N-1}$。

**定理 A.20（径向分布的极分解）。** 设 $X$ 是 $\mathbb R^n$ 中的随机向量，其密度形如 $f(|x|)$。则
$X\neq0$ 几乎必然成立，$|X|$ 在 $(0,\infty)$ 上有密度 $s_nr^{n-1}f(r)$，$X/|X|$ 服从
$\sigma_{n-1}$，并且 $|X|$ 与 $X/|X|$ 相互独立。此外，对每个 $O\in O(n)$ 有 $\sigma_{n-1}(OB)=\sigma_{n-1}(B)$。
特别地：

(a) 对 $\mathbb R^n$ 的标准高斯向量 $g$，$g/|g|\sim\sigma_{n-1}$，且与 $|g|$ 独立；

(b) 对标准复高斯向量 $z\in\mathbb C^N$，$\psi=z/|z|\sim\sigma_{2N-1}$，这一分布在 $U(N)$ 下不变，且与
$|z|$ 独立；$|z|^2$ 在 $(0,\infty)$ 上有密度 $t^{N-1}e^{-t}/(N-1)!$，所以 $\mathbb E|z|^4=N(N+1)$。

*证明。* 由 $h=f\cdot1_A$，$\phi=1_B$ 时的定理 A.18：
$\mathbb P(|X|\in A,\ X/|X|\in B)=n\,\lambda_n(C(B))\int_Af(r)r^{n-1}dr$。取 $A=(0,\infty)$ 和 $B=S^{n-1}$，这就是
$1=nV_n\int_0^\infty f(r)r^{n-1}dr$（并且 $\mathbb P(X=0)=0$）。取 $A=(0,\infty)$ 得 $\mathbb P(X/|X|\in B)=\sigma_{n-1}(B)$，
取 $B=S^{n-1}$ 得 $\mathbb P(|X|\in A)=s_n\int_Af(r)r^{n-1}dr$，而一般情形是二者的乘积。
对 $O\in O(n)$，$OC(B)=C(OB)$ 且 $\lambda_n(OE)=\lambda_n(E)$（引理 A.3(b)）。(a) $\varphi_n(x)=(2\pi)^{-n/2}e^{-|x|^2/2}$ 是径向的。
(b) $z=g/\sqrt2$，其中 $g$ 是 $\mathbb R^{2N}$ 中的标准高斯向量（命题 A.16），所以 $z/|z|=g/|g|$，$|z|=|g|/\sqrt2$；
酉矩阵在 $\mathbb R^{2N}$ 上是正交的。$|g|$ 的密度是 $s_{2N}(2\pi)^{-N}r^{2N-1}e^{-r^2/2}$，其中 $s_{2N}=2\pi^N/(N-1)!$，
代换 $t=r^2/2$ 给出 $|z|^2$ 的密度；于是 $\mathbb E|z|^4=\int_0^\infty t^{N+1}e^{-t}dt/(N-1)!=(N+1)N$。$\square$

引理 7.9 就是 (b) 部分。当 $N=1$ 时，$\sigma_1$ 是归一化的弧长，即 $[0,2\pi)$ 上的
$\frac1{2\pi}\lambda_1$ 在 $\theta\mapsto e^{i\theta}$ 下的前推测度：对 Borel 集 $B\subseteq S^1$，令
$\Theta(B)=\{\theta\in(0,2\pi):e^{i\theta}\in B\}$。映射 $(r,\theta)\mapsto re^{i\theta}$ 是从开集
$(0,1)\times(0,2\pi)$ 到去掉线段 $[0,1)$ 的开单位圆盘上的 $C^1$ 双射，Jacobi 行列式为 $r$，而 $C(B)$ 与
$(0,1)\times\Theta(B)$ 的像只相差该线段与单位圆周的一个子集，这是零测集。由变量替换公式，
$\lambda_2(C(B))=\int_{\Theta(B)}\int_0^1r\,dr\,d\theta=\frac12\lambda_1(\Theta(B))$，再由 $V_2=\pi$ 得
$\sigma_1(B)=\lambda_1(\Theta(B))/2\pi$（点 $\theta=0$ 是零测集）。

**推论 A.21（均匀单位向量的矩）。** 若 $u\sim\sigma_{2N-1}$ 在 $\mathbb C^N$ 的单位球面上，则
$\mathbb E[u_a\bar u_b]=\delta_{ab}/N$，$\mathbb E|u_1|^4=\frac2{N(N+1)}$，$\mathbb E|u_1|^2|u_2|^2=\frac1{N(N+1)}$（$N\ge2$）。若
$u\sim\sigma_{n-1}$ 在 $S^{n-1}\subset\mathbb R^n$ 上，则 $\mathbb Eu_1^2=\frac1n$，$\mathbb Eu_1^4=\frac3{n(n+2)}$，
$\mathbb Eu_1^2u_2^2=\frac1{n(n+2)}$。

*证明。* 如定理 A.20(b) 中那样取 $u=z/|z|$。由于 $z_a\bar z_b=|z|^2u_a\bar u_b$，且 $|z|$ 与 $u$ 相互独立，
$\delta_{ab}=\mathbb E[z_a\bar z_b]=\mathbb E|z|^2\,\mathbb E[u_a\bar u_b]=N\,\mathbb E[u_a\bar u_b]$（由命题 A.16(a)，
$\mathbb E|z|^2=N$）。同样 $\mathbb E|z_1|^4=N(N+1)\mathbb E|u_1|^4$，
$\mathbb E|z_1|^2|z_2|^2=N(N+1)\mathbb E|u_1|^2|u_2|^2$，其中由独立性 $\mathbb E|z_1|^2|z_2|^2=1$，而对相互独立的
$x,y\sim N(0,\frac12)$，$\mathbb E|z_1|^4=\mathbb E(x^2+y^2)^2=\frac34+2\cdot\frac14+\frac34=2$。（对
$Y\sim N(0,s^2)$，分部积分给出 $\mathbb EY^4=3s^2\,\mathbb EY^2=3s^4$，因为
$y^3\cdot ye^{-y^2/2s^2}=-s^2y^3\frac d{dy}e^{-y^2/2s^2}$。）在实情形，$\mathbb E|g|^2=n$，
$\mathbb E|g|^4=s_n(2\pi)^{-n/2}I_{n+4}=I_{n+4}/I_n=n(n+2)$（推论 A.19），$\mathbb Eg_1^4=3$，$\mathbb Eg_1^2g_2^2=1$，并且
$|g|^2$，$|g|^4$ 与 $u=g/|g|$ 独立（定理 A.20(a)）。$\square$

这些与引理 6.14 一致，后者仅由酉不变性推出它们。反过来，不变性决定了分布：

**命题 A.22（不变分布的唯一性）。** 设 $u$ 是 $\mathbb C^N$ 中的随机单位向量，其分布在所有酉矩阵下不变。
则 $u\sim\sigma_{2N-1}$。对 $N=1$：圆周上每个旋转不变的 Borel 概率测度都是归一化的弧长。

*证明。* 设 $Q$ 在 $U(N)$ 上服从 Haar 分布（定理 5.3），且与 $u$ 独立。对固定的单位向量 $v$，取
$W\in U(N)$ 使 $We_1=v$。由引理 5.1，$QW$ 服从 Haar 分布，所以 $Qv=(QW)e_1$ 与 $Qe_1$ 同分布。在
定理 5.3 中，$Q$ 是 Ginibre 矩阵的 Gram–Schmidt 因子（引理 5.2），所以对其第一列 $z$ 有 $Qe_1=z/|z|$，
由定理 A.20(b)，它服从 $\sigma_{2N-1}$。对 $u$ 取条件（Fubini，独立性），对每个 Borel 集 $B$ 有
$\mathbb P(Qu\in B)=\sigma_{2N-1}(B)$。若改为对 $Q$ 取条件，并利用对每个固定的酉矩阵 $U$，$Uu$ 与 $u$ 同分布，
则得 $\mathbb P(Qu\in B)=\mathbb P(u\in B)$。圆周上的概率测度就是恒等映射的分布，而
$U(1)$ 以旋转的方式作用。$\square$

```python
N, B = 3, 400_000                                   # Theorem A.20 (b), Corollary A.21, in C^3
z = (rng.standard_normal((B, N)) + 1j * rng.standard_normal((B, N))) / math.sqrt(2)
R2 = (np.abs(z) ** 2).sum(axis=1)
u = z / np.sqrt(R2)[:, None]
t = 2.0
p_gamma = 1 - math.exp(-t) * sum(t ** k / math.factorial(k) for k in range(N))   # P(|z|^2 <= t) for Gamma(N)
mc(R2 <= t, p_gamma, "P(|z|^2 <= 2)")
mc(R2 ** 2, N * (N + 1), "E|z|^4")
mc(np.abs(u[:, 0]) ** 4, 2 / (N * (N + 1)), "E|u_1|^4")
mc(np.abs(u[:, 0] * u[:, 1]) ** 2, 1 / (N * (N + 1)), "E|u_1|^2 |u_2|^2")
mc((R2 <= t) * np.abs(u[:, 0]) ** 4, p_gamma * 2 / (N * (N + 1)), "E[1{|z|^2 <= 2} |u_1|^4], independence")
mc((u[:, 0] * u[:, 1].conj()).real, 0.0, "E Re u_1 conj(u_2)")

x3 = rng.uniform(-1, 1, (B, 3))                     # Corollary A.19: V_3 = s_3 / 3 = 4 pi / 3
mc((x3 ** 2).sum(axis=1) <= 1, math.pi / 6, "P(uniform point of [-1,1]^3 in the unit ball)")
g3 = rng.standard_normal((B, 3))
u3 = g3 / np.linalg.norm(g3, axis=1, keepdims=True)
mc(u3[:, 0] ** 4, 3 / 15, "E u_1^4 on S^2")
phase = np.angle(z[:, 0])                           # sigma_1 is normalised arc length
mc((phase >= 0) & (phase < 1), 1 / (2 * math.pi), "P(arg z_1 in [0, 1))")
```

## A.6 高斯矩：Isserlis–Wick 公式

**引理 A.23.** 设 $Y$ 是 $\mathbb R^n$ 中协方差为 $C=(c_{ij})$ 的中心化高斯向量。则
$c_{ij}=\mathbb E[Y_iY_j]$，并且对所有 $t\in\mathbb R^n$ 有 $\mathbb Ee^{\langle t,Y\rangle}=e^{\langle t,Ct\rangle/2}$ 和
$\mathbb Ee^{|\langle t,Y\rangle|}<\infty$，所有矩 $\mathbb E|Y_i|^k$ 都有限。

*证明。* 写 $Y=Ag$。则 $\mathbb E[Y_iY_j]=\sum_{k,l}A_{ik}A_{jl}\mathbb E[g_kg_l]=(AA^T)_{ij}$。令 $a=A^Tt$，则
$\langle t,Y\rangle=\sum_ka_kg_k$，由独立性
$\mathbb Ee^{\langle t,Y\rangle}=\prod_k\mathbb Ee^{a_kg_k}=\prod_ke^{a_k^2/2}=e^{|A^Tt|^2/2}=e^{\langle t,Ct\rangle/2}$。最后
$e^{|s|}\le e^s+e^{-s}$，且 $|Y_i|^k\le k!\,e^{|Y_i|}$。$\square$

$\{1,\dots,2p\}$ 的一个**配对**是把它划分为 $p$ 个二元块；这样的配对共有 $(2p-1)!!=(2p)!/(2^pp!)$ 个
（$1$ 的配偶有 $2p-1$ 种选法，然后对其余 $2p-2$ 个元素配对；对 $p$ 归纳）。

**定理 A.24（Isserlis–Wick）。** 设 $Y_1,\dots,Y_m$ 是一个中心化高斯向量的坐标（允许重复）。
若 $m$ 为奇数，则 $\mathbb E[Y_1\cdots Y_m]=0$；若 $m=2p$，则
$$
\mathbb E[Y_1\cdots Y_{2p}]=\sum_\pi\prod_{\{a,b\}\in\pi}\mathbb E[Y_aY_b],
$$
其中求和取遍 $\{1,\dots,2p\}$ 的所有配对 $\pi$。对于复随机变量 $Z_1,\dots,Z_m$，若它们是同一个中心化高斯向量的
坐标的复线性组合（例如 Ginibre 矩阵或 GUE 矩阵的元素及其复共轭），则同一公式成立，配对因子取
$\mathbb E[Z_aZ_b]$（不取共轭）。

*证明。* 向量 $(Y_1,\dots,Y_m)$ 是给定高斯向量在一个线性映射下的像，因而是中心化高斯向量
（定理 A.15(a)）；设 $C$ 为其协方差。固定 $t\in\mathbb R^m$，令 $S=\langle t,Y\rangle$。对 $s\in\mathbb R$，
$\sum_k(sS)^k/k!$ 的部分和以 $e^{|sS|}$ 为界，由引理 A.23 它可积，所以由控制收敛定理和引理 A.23，
$$
\sum_{k\ge0}\frac{s^k}{k!}\mathbb E[S^k]=\mathbb Ee^{sS}=e^{s^2q/2}=\sum_{p\ge0}\frac{s^{2p}q^p}{2^pp!},\qquad q=\langle t,Ct\rangle .
$$
两边都是 A.2 节意义下 $\mathbb R$ 上的幂级数：右边显然如此，而对左边，对每个 $\rho>0$ 有
$\sum_k\frac{\rho^k}{k!}|\mathbb E[S^k]|\le\mathbb E\sum_k\frac{\rho^k|S|^k}{k!}=\mathbb Ee^{\rho|S|}<\infty$
（单调收敛定理和引理 A.23）。所以它们的系数相同（对二者之差应用引理 A.5(c)）。因此对每个
$t\in\mathbb R^m$，$\mathbb E\langle t,Y\rangle^m/m!$ 当 $m$ 为奇数时等于 $0$，当 $m=2p$ 时等于 $\langle t,Ct\rangle^p/(2^pp!)$。
两边都是 $t$ 的多项式，所以它们的系数相同（对二者之差应用推论 A.7(a)）。
在 $\mathbb E\langle t,Y\rangle^m=\sum_{i_1,\dots,i_m}t_{i_1}\cdots t_{i_m}\mathbb E[Y_{i_1}\cdots Y_{i_m}]$ 中，单项式 $t_1t_2\cdots t_m$
来自那 $m!$ 个是 $(1,\dots,m)$ 的排列的序列 $(i_1,\dots,i_m)$，所以它在
$\mathbb E\langle t,Y\rangle^m/m!$ 中的系数是 $\mathbb E[Y_1\cdots Y_m]$。因此当 $m$ 为奇数时它为 $0$。当 $m=2p$ 时，
$q^p=\sum c_{i_1j_1}\cdots c_{i_pj_p}\,t_{i_1}t_{j_1}\cdots t_{i_p}t_{j_p}$，求和取遍由 $p$ 个有序对组成的所有序列；单项式
$t_1\cdots t_{2p}$ 恰在无序对 $\{i_k,j_k\}$ 构成一个配对 $\pi$ 时出现，而每个配对来自
$p!\,2^p$ 个序列（各对的次序，以及每对内部的次序），由于 $C$ 对称，每个序列贡献 $\prod_{\{a,b\}\in\pi}c_{ab}$。
除以 $2^pp!$ 即得该公式，因为 $c_{ab}=\mathbb E[Y_aY_b]$。
对复的情形，写 $Z_a=\sum_i\alpha_{ai}Y_i$。两边关于 $(Z_1,\dots,Z_m)$ 都是多重线性的——左边由 $\mathbb E$ 的
线性性，右边是因为每个 $Z_a$ 在每一项中恰好出现在一个因子里——并且当每个 $Z_a$ 都是某个 $Y_i$ 时两边相等；
展开即可。$\square$

例如，对标准复高斯变量 $z$，$\mathbb E|z|^4=\mathbb E[z\bar zz\bar z]=1\cdot1+0+1\cdot1=2$。命题 1.10
和引理 7.10 用到 $m=4$ 的情形；定理 A.24 无需在期望号下求导就证明了它。

```python
def pairings(idx):
    """All pairings of the list idx, as lists of pairs."""
    if not idx:
        yield []
        return
    for i in range(1, len(idx)):
        for rest in pairings(idx[1:i] + idx[i + 1:]):
            yield [(idx[0], idx[i])] + rest


def wick(cov, idx):                                       # sum over pairings of the positions of idx
    return sum(math.prod(cov[idx[a], idx[b]] for a, b in p) for p in pairings(list(range(len(idx)))))


assert [len(list(pairings(list(range(2 * p))))) for p in range(1, 5)] == [1, 3, 15, 105]    # (2p-1)!!
Amat = rng.standard_normal((3, 3)) / math.sqrt(3)
cov = Amat @ Amat.T
Y = rng.standard_normal((1_000_000, 3)) @ Amat.T
for idx in ((0, 1, 2, 0), (0, 0, 0, 0), (0, 1, 1, 2, 2, 0)):
    mc(np.prod(Y[:, list(idx)], axis=1), wick(cov, idx), f"E Y{idx}")
mc(np.prod(Y[:, [0, 1, 2]], axis=1), 0.0, "E Y(0, 1, 2), odd")
M = np.array([[1, 1j], [1, -1j]])                  # (z, conj z) = M (x, y) with x, y independent N(0, 1/2)
cz = M @ M.T / 2                                     # E[Z_a Z_b] = [[0, 1], [1, 0]]
assert abs(wick(cz, (0, 1, 0, 1)) - 2) < 1e-12                           # E|z|^4 = 2
```

## A.7 指数变量与均匀间距

若随机变量 $e$ 在 $(0,\infty)$ 上有密度 $e^{-y}$，则称它是**指数变量**（$e\sim\operatorname{Exp}(1)$）；
$e/\beta$ 是速率为 $\beta$ 的指数变量。令 $\Delta_n=\{x\in\mathbb R^n:0<x_1<\dots<x_n<1\}$。

**定理 A.25.** 设 $e_0,\dots,e_n$（$n\ge1$）是相互独立的 $\operatorname{Exp}(1)$ 变量，$T=e_0+\dots+e_n$，并对
$k=1,\dots,n$ 令 $X_k=(e_0+\dots+e_{k-1})/T$。则 $(X_1,\dots,X_n)$ 在 $\Delta_n$ 上有密度 $n!$，$T$ 在 $(0,\infty)$ 上有密度
$t^ne^{-t}/n!$，并且 $(X_1,\dots,X_n)$ 与 $T$ 相互独立。特别地，
$(e_0,\dots,e_n)/T$ 与 $T$ 独立。

*证明。* 映射 $\Phi(x,t)=t\,(x_1,\ x_2-x_1,\ \dots,\ x_n-x_{n-1},\ 1-x_n)$ 是从开集
$\Delta_n\times(0,\infty)$ 到 $(0,\infty)^{n+1}$ 上的双射，其逆为 $y\mapsto\big(((y_0+\dots+y_{k-1})/\sum_iy_i)_{k=1}^n,\ \sum_iy_i\big)$；
二者都是 $C^1$ 的。在 $\Phi$ 的 Jacobi 矩阵中，$\partial/\partial t$ 那一列的元素之和为 $1$，而每一列
$\partial/\partial x_j$（元素为 $t$ 和 $-t$）的元素之和为 $0$。把最后一行换成所有行之和（这
不改变行列式），得到行 $(0,\dots,0,1)$，而剩下的 $n\times n$ 块是对角元为 $t$ 的下三角矩阵；
所以 $|\det D\Phi|=t^n$。向量 $e=(e_0,\dots,e_n)$ 有密度 $e^{-\sum_iy_i}$，并且
$\sum_i\Phi(x,t)_i=t$。由变量替换公式，对 Borel 集 $E\subseteq\Delta_n\times(0,\infty)$，
$$
\mathbb P\big((X,T)\in E\big)=\mathbb P\big(e\in\Phi(E)\big)=\int_{\Phi(E)}e^{-\sum_iy_i}dy=\int_Ee^{-t}t^n\,dx\,dt .
$$
密度 $e^{-t}t^n=n!\cdot\frac{t^ne^{-t}}{n!}$ 是乘积形式，且 $\int_0^\infty t^ne^{-t}dt=n!$，所以
$\lambda_n(\Delta_n)=1/n!$；由 Fubini 定理，$\mathbb P(X\in A,T\in A')=\mathbb P(X\in A)\mathbb P(T\in A')$。最后，$(e_0,\dots,e_n)/T$ 等于
$(X_1,X_2-X_1,\dots,1-X_n)$，是 $X$ 的函数。$\square$

**推论 A.26（均匀间距）。** 设 $U_1,\dots,U_n$ 相互独立且在 $(0,1)$ 上均匀分布，其顺序统计量为
$U_{(1)}<\dots<U_{(n)}$，间距为 $s_0=U_{(1)}$，$s_k=U_{(k+1)}-U_{(k)}$（$1\le k\le n-1$），$s_n=1-U_{(n)}$。则：

(a) 出现相等值的概率为 $0$，并且 $(U_{(1)},\dots,U_{(n)})$ 在 $\Delta_n$ 上有密度 $n!$；

(b) $(s_0,\dots,s_n)$ 与 $(e_0,\dots,e_n)/T$ 同分布，其中 $e_i$，$T$ 如定理 A.25；

(c) $\mathbb Es_i=\frac1{n+1}$，$\mathbb Es_i^2=\frac2{(n+1)(n+2)}$，并且对 $i\neq j$ 有 $\mathbb E[s_is_j]=\frac1{(n+1)(n+2)}$；

(d) 对 $i\neq j$，$s_i/s_j$ 与 $e_i/e_j$ 同分布，且 $\mathbb P(e_i/e_j\le t)=\frac t{1+t}$（$t>0$）；所以对
$0\le i\le n-1$，间距比 $\min(s_i,s_{i+1})/\max(s_i,s_{i+1})$ 与命题 5.14(a) 中求得的 $\min(e_i,e_{i+1})/\max(e_i,e_{i+1})$
同分布：在 $[0,1]$ 上 $\mathbb P(r\le t)=\frac{2t}{1+t}$，均值为 $2\ln2-1$。

对区间 $[a,b]$ 上 $n$ 个相互独立的均匀分布点，仿射映射 $u\mapsto a+(b-a)u$ 把所有间距乘以
$b-a$：间距除以 $b-a$ 后满足 (b) 与 (c)，而 (d) 中的比值不变。

*证明。* (a) 两个坐标相等的集合是有限多个超平面的并，是零测集（推论 A.7(c)）。在其之外，
恰有一个置换 $\pi$ 使 $U_{\pi(1)}<\dots<U_{\pi(n)}$，所以对 Borel 集 $E\subseteq\Delta_n$，事件
$\{(U_{(1)},\dots,U_{(n)})\in E\}$ 在相差一个零测集的意义下，是对 $\pi$ 取的不交并
$\{(U_{\pi(1)},\dots,U_{\pi(n)})\in E\}$。每一个的概率都是 $\lambda_n(E)$，因为 $(U_1,\dots,U_n)$ 有密度
$1_{(0,1)^n}$，而坐标置换同时保持 $\lambda_n$ 和 $(0,1)^n$。所以该概率为 $n!\,\lambda_n(E)$。
(b) 由 (a) 和定理 A.25，向量 $(U_{(k)})_k$ 与 $(X_k)_k$ 同分布，而间距是二者的同一个线性
函数。(c) 令 $w=(e_0,\dots,e_n)/T$，它与 $T$ 独立（定理 A.25），且 $\mathbb ET=n+1$，
$\mathbb ET^2=\int_0^\infty t^{n+2}e^{-t}dt/n!=(n+1)(n+2)$。由 $e_i=Tw_i$：$1=\mathbb Ee_i=\mathbb ET\,\mathbb Ew_i$，
$2=\mathbb Ee_i^2=\mathbb ET^2\,\mathbb Ew_i^2$，$1=\mathbb E[e_ie_j]=\mathbb ET^2\,\mathbb E[w_iw_j]$。(d) $s_i/s_j=w_i/w_j=e_i/e_j$，并由
Tonelli 定理，$\mathbb P(e_i\le te_j)=\int_0^\infty e^{-y}(1-e^{-ty})\,dy=1-\frac1{1+t}$。间距比是
$s_i/s_{i+1}$ 的函数（即在 $x=s_i/s_{i+1}$ 处取值的 $\min(x,1/x)$，由于所有 $e_k>0$，它几乎必然有定义），所以它与
$e_i/e_{i+1}$ 的同一函数同分布，而后者正是命题 5.14(a) 中对相互独立的 $\operatorname{Exp}(1)$ 变量
$e_i,e_{i+1}$ 的比值 $r$。$\square$

特别地，考虑"Poisson 能级" $\lambda_k=e_0+\dots+e_{k-1}$，$k=1,\dots,n+1$（相互独立的指数
间距 $e_0,\dots,e_n$），以及 $n$ 个排好序的均匀分布点及其 $n+1$ 个间距 $s_0,\dots,s_n$。相邻间距比的向量
$r_i=\min(d_i,d_{i+1})/\max(d_i,d_{i+1})$，$0\le i<n$，对 $d=(e_0,\dots,e_n)$ 和对 $d=(s_0,\dots,s_n)$ 有相同的联合分布，
因为把所有 $d_i$ 都除以 $T$ 时这些比值不变，而
$(e_0,\dots,e_n)/T$ 与 $(s_0,\dots,s_n)$ 同分布（(b) 部分）；并且由定理 A.25，前 $n$ 个 Poisson 能级除以第 $n+1$ 个
所得的结果，恰好与 $n$ 个排好序的均匀分布点同分布。

```python
n, B = 5, 200_000                                       # Corollary A.26
levels = np.sort(rng.uniform(0, 1, (B, n)), axis=1)
s = np.diff(np.concatenate([np.zeros((B, 1)), levels, np.ones((B, 1))], axis=1), axis=1)    # s_0, ..., s_n
e = rng.exponential(1.0, (B, n + 1))
w = e / e.sum(axis=1, keepdims=True)                    # normalised exponentials
for name, sp in (("uniform spacings", s), ("normalised exponentials", w)):
    mc(sp[:, 2], 1 / (n + 1), f"{name}: E s_2")
    mc(sp[:, 2] ** 2, 2 / ((n + 1) * (n + 2)), f"{name}: E s_2^2")
    mc(sp[:, 0] * sp[:, n], 1 / ((n + 1) * (n + 2)), f"{name}: E s_0 s_n")
    for tt in (0.5, 2.0):
        mc(sp[:, 1] / sp[:, 2] <= tt, tt / (1 + tt), f"{name}: P(s_1 / s_2 <= {tt})")
r = np.minimum(s[:, 1:-2], s[:, 2:-1]) / np.maximum(s[:, 1:-2], s[:, 2:-1])    # ratios of interior spacings
mc(r.mean(axis=1), 2 * math.log(2) - 1, "<r>, uniform levels")
T = e.sum(axis=1)
mc(w[:, 0] * T ** 2, n + 2, "E[w_0 T^2] = E w_0 E T^2 (independence)")
```

**返回**[目录](README.md)。
