[English](../03_traces_and_states.md) · [Español](../es/03_traces_and_states.md) · [Français](../fr/03_traces_and_states.md) · [Deutsch](../de/03_traces_and_states.md) · **中文** · [日本語](../ja/03_traces_and_states.md) · [Русский](../ru/03_traces_and_states.md)

# 3. 迹、态与熵

在 $M_n(\mathbb C)$ 上本质上只有一个迹。在 $\bigoplus_cM_{k_c}$ 上则有一个
$C$-参数族，而且迹的选取很重要：密度、熵以及映射的对偶都依赖于它。本章对迹进行分类，
用密度表示态，证明基本的熵不等式，最后给出一个看似合理却不成立的迹不等式。

**你将学到**本库中的三种迹以及它们之间的换算，态的经典部分，冯·诺依曼熵在各扇区上的
分解，Klein 不等式及其推论，以及 Michelson 对比度。

在本章中，$\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$，所有 $k_c\ge1$，
$D=\sum_ck_c$，$P_c$ 为极小中心投影，$A=\bigoplus_cA_c$。

## 3.1 迹

$\mathcal M$ 上的**迹**是对所有 $A,B$ 满足 $\tau(AB)=\tau(BA)$ 的线性泛函 $\tau$。若对所有
$A$ 有 $\tau(A^*A)\ge0$，则称它是**正的**；若此外
$\tau(A^*A)=0$ 蕴含 $A=0$，则称它是**忠实的**。

**定理 3.1.** $\mathcal M$ 上的每个迹都具有形式
$\operatorname{tr}_w(A)=\sum_cw_c\operatorname{Tr}A_c$，其中 $w\in\mathbb C^C$。它是正的
当且仅当所有 $w_c\ge0$，是忠实的当且仅当所有 $w_c>0$。

*证明。* 设 $E^c_{ij}$ 为第 $c$ 个块的矩阵单位。对迹 $\tau$，
$\tau(E^c_{ij})=\tau(E^c_{i1}E^c_{1j})=\tau(E^c_{1j}E^c_{i1})=\delta_{ij}\tau(E^c_{11})$。
令 $w_c=\tau(E^c_{11})$；由线性得 $\tau(A)=\sum_c\sum_i(A_c)_{ii}w_c$。反之，
每个 $\operatorname{tr}_w$ 都是迹。正性：$E^c_{11}=(E^c_{11})^*E^c_{11}$，所以
$w_c=\tau\big((E^c_{11})^*E^c_{11}\big)\ge0$ 是必要的；又因 $E^c_{11}\neq0$，$w_c>0$ 是忠实性的必要条件。
反之，$\operatorname{tr}_w(A^*A)=\sum_cw_c\|A_c\|_2^2$：若所有 $w_c\ge0$，它 $\ge0$；若所有 $w_c>0$，它仅当
$A=0$ 时为零。$\square$

本库实现了三种忠实迹：

| 名称 | 权重 $w_c$ | 在 $1$ 上的值 | 含义 |
|---|---|---|---|
| `Tr_blunt` | $1$ | $D$ | $\bigoplus_c\mathbb C^{k_c}$ 的迹（重数 $m_c=1$） |
| `Tr_norm` | $1/k_c$ | $C$ | 每个块都归一化 |
| `tau_vN` | $1/(Ck_c)$ | $1$ | 各扇区权重均为 $1/C$ 的迹态 |

由 2.3 节，$\mathcal M$ 以重数 $m_c$ 作用的空间的迹就是
$w_c=m_c$ 的 $\operatorname{tr}_w$。因此，选取一个忠实迹就是选定"每个扇区占多大分量"。

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, DensityMatrix, states, channels

torch.manual_seed(0)
alg = TypeIAlgebra([2, 3, 5], [2, 3, 5], precision='double')    # M_2 (+) M_3 (+) M_5, D = 10, C = 3
U = alg.random_unitary_operator(4)
A = alg.operator_from_eigenvalues(lambda d: torch.randn(4, d), batch_size=4, force_self_adjoint=True)
B = U @ A                                                       # a batch of four non-normal operators
one = alg.identity()

for f, value in ((lambda X: X.Tr_blunt(), 10), (lambda X: X.Tr_norm(), 3), (lambda X: X.tau_vN(), 1)):
    assert abs(f(one).item() - value) < 1e-12                   # normalisation
    assert (f(A @ B) - f(B @ A)).abs().max() < 1e-10            # trace property
assert (A.Tr_norm() - 3 * A.tau_vN()).abs().max() < 1e-12       # Tr_norm = C tau_vN
```

## 3.2 态与密度

**态**是 $\mathcal M$ 上一个正的
（$\omega(A^*A)\ge0$）且归一化的（$\omega(1)=1$）线性泛函 $\omega$。

**定理 3.2.** 设 $\operatorname{tr}_w$ 是忠实迹。对每个线性泛函
$\omega$，存在唯一的 $\rho_w\in\mathcal M$，使得对所有
$A$ 有 $\omega(A)=\operatorname{tr}_w(\rho_wA)$。泛函 $\omega$ 是正的当且仅当 $\rho_w\ge0$；
它是态当且仅当此外还有
$\operatorname{tr}_w(\rho_w)=1$。对两个忠实迹，相应的密度按块满足
$$
w_c\,(\rho_w)_c=w'_c\,(\rho_{w'})_c .
$$

*证明。* $\langle X,Y\rangle_w=\operatorname{tr}_w(X^*Y)$ 是 $\mathcal M$ 上的内积
（定理 3.1）。由有限维空间 $\mathcal M$ 中的 Riesz 表示定理，
存在唯一的 $R$ 使得 $\omega(A)=\langle R,A\rangle_w=\operatorname{tr}_w(R^*A)$；令
$\rho_w=R^*$。若 $\rho_w\ge0$，则
$\omega(A^*A)=\sum_cw_c\operatorname{Tr}(A_c\,(\rho_w)_c\,A_c^*)\ge0$，因为
$A_c(\rho_w)_cA_c^*\ge0$（取 $T=(\rho_w)_c^{1/2}A_c^*$ 应用推论 1.7），正矩阵的迹
$\ge0$（其特征值 $\ge0$，推论 1.4(3)），并且 $w_c>0$。反之，若
$\omega$ 是正的，取 $x\in\mathbb C^{k_c}$，并令 $A$ 为块 $c$ 中的秩一算子
$|x\rangle\langle x|$（$x\neq0$）；则 $A=A^*A/\|x\|^2$，且
$0\le\omega(A)=w_c\langle x,(\rho_w)_cx\rangle$，所以 $(\rho_w)_c\ge0$。归一化条件就是
$\omega(1)=\operatorname{tr}_w(\rho_w)$。最后，$\sum_cw_c\operatorname{Tr}((\rho_w)_cA_c)$ 与
$\sum_cw'_c\operatorname{Tr}((\rho_{w'})_cA_c)$ 对所有 $A$ 都等于 $\omega(A)$，所以 $\bigoplus_cw_c(\rho_w)_c$ 与
$\bigoplus_cw'_c(\rho_{w'})_c$ 都是 $\omega$ 关于 blunt 迹（$w_c=1$）的密度，由唯一性得
$w_c(\rho_w)_c=w'_c(\rho_{w'})_c$。$\square$

本库通过态关于 `Tr_blunt` 的密度 $\rho=\rho_{\mathrm{blunt}}$ 来存储态
（类 `DensityMatrix`）；`density(trace)` 和 `DensityMatrix.from_density`
负责换算。态的取值与约定无关：

```python
rho = states.random_density_matrix(alg, batch_size=4)          # Hilbert-Schmidt random states
for trace, tr in (('blunt', lambda X: X.Tr_blunt()), ('norm', lambda X: X.Tr_norm()),
                  ('tau_vN', lambda X: X.tau_vN())):
    d = alg.operator(rho.density(trace))
    assert (tr(d @ A) - (rho @ A).Tr_blunt()).abs().max() < 1e-10       # omega(A) for all traces
    assert (tr(d) - 1).abs().max() < 1e-10                               # tr_w(rho_w) = 1
    back = DensityMatrix.from_density(alg, rho.density(trace), trace=trace)
    assert (back.matrix - rho.matrix).abs().max() < 1e-12
```

> **物理意义**（物理解释，而非本课程的定理）。$\omega(A)$ 是可观测量 $A$ 的期望值。
> blunt 迹的密度 $\rho$ 就是量子力学中 $\bigoplus_c\mathbb C^{k_c}$ 上通常的密度矩阵。
> 当同一个代数以一定重数作用时（第 2 章），就会出现关于其他迹的密度，例如在自旋系统的
> 对称部分与非对称部分上。

## 3.3 态的经典部分

$\omega$ 在中心上的限制是各扇区上的概率分布，
$p_c=\omega(P_c)=\operatorname{Tr}\rho_c$。若 $p_c>0$，扇区 $c$ 中的**条件态**
在 $\mathbb C^{k_c}$ 上的密度为 $\sigma_c=\rho_c/p_c$，并且
$$
\rho=\bigoplus_cp_c\,\sigma_c
$$
（若 $p_c=0$，则 $\rho_c=0$，因为 $\rho_c\ge0$ 的特征值 $\ge0$ 且其和为
$\operatorname{Tr}\rho_c=0$（推论 1.4(3)），从而全为 $0$，而 $\rho_c$ 可酉对角化（定理 1.3）；此时略去该项）。
因此 $\mathcal M$ 上的每个态都是各扇区的态的凸组合。映射
$\rho\mapsto\sum_cP_c\rho P_c$ 在 $\mathcal M$ 上是恒等映射，因为
$\mathcal M$ 的每个元素都是块对角的。

> **物理意义**（物理解释，而非本课程的定理）。$\mathcal M$ 上的态是各扇区中量子态的经典混合。
> 映射 $\rho\mapsto\sum_cP_c\rho P_c$ 是对"处于哪个扇区"的测量，它不会扰动 $\mathcal M$ 的
> 任何态。这正是超选择定则的操作含义。

```python
p = rho.sector_probabilities()                                   # shape (batch, C)
assert torch.allclose(p.sum(-1), torch.ones(4, dtype=p.dtype))
P = [alg.central([1.0 if c == d else 0.0 for d in range(alg.C)]) for c in range(alg.C)]
measure_sector = channels.measurement_channel(P)
assert (measure_sector(rho).matrix - rho.matrix).abs().max() < 1e-12   # no disturbance
```

## 3.4 熵

具有（blunt）密度 $\rho$ 的态的**冯·诺依曼熵**为
$S(\rho)=-\operatorname{Tr}\rho\log\rho=-\sum_i\lambda_i\log\lambda_i$，求和遍及
$\rho$ 的特征值，约定 $0\log0=0$；也就是说，对 $f(t)=t\log t$，$f(0)=0$，令 $\rho\log\rho:=f(\rho)$，
并且由 $f(\rho)=U\,\mathrm{diag}(f(\lambda_i))\,U^*$ 的定义（1.3 节）有 $\operatorname{Tr}f(\rho)=\sum_if(\lambda_i)$。

**定理 3.3（扇区混合的熵）。** 设 $p_c$ 与 $\sigma_c$ 如上
（略去 $p_c=0$ 的扇区），则
$$
S(\rho)=H(p)+\sum_cp_c\,S(\sigma_c),\qquad H(p)=-\sum_cp_c\log p_c .
$$

*证明。* $\rho$ 的非零特征值（计重数）是数 $p_c\mu$，其中 $\mu$ 取遍各 $\sigma_c$（$p_c>0$ 的扇区）的
非零特征值；零特征值对两边的贡献都是 $0$。因此
$-\sum_c\sum_{\mu>0}p_c\mu\log(p_c\mu)=-\sum_cp_c\log p_c\sum_{\mu>0}\mu-\sum_cp_c\sum_{\mu>0}\mu\log\mu$，
而 $\sum_{\mu>0}\mu=\operatorname{Tr}\sigma_c=1$。$\square$

> **物理意义**（物理解释，而非本课程的定理）。这个公式把关于系统的不确定性分成两部分：$H(p)$，即关于荷（系统处于哪个扇区）
> 的经典不确定性；以及 $\sum_cp_cS(\sigma_c)$，即扇区内部的平均量子不确定性。在超选择定则下，
> 第一部分是一个经典随机变量——荷——的熵。

相对于另一个迹的熵为 $S_w(\rho)=-\operatorname{tr}_w(\rho_w\log\rho_w)$，采用同样的约定
$\rho_w\log\rho_w:=f(\rho_w)$。由
定理 3.2，$(\rho_w)_c=\rho_c/w_c$，所以
$$
S_w(\rho)=-\sum_cw_c\operatorname{Tr}\Big(\frac{\rho_c}{w_c}\log\frac{\rho_c}{w_c}\Big)
=S(\rho)+\sum_cp_c\log w_c .
$$
熵是态*与*迹共同的性质：同一个态关于 `Tr_blunt` 和 `tau_vN` 具有不同的熵。

```python
S = rho.von_neumann_entropy()
sigma = [rho.matrix[:, c, :k, :k] / p[:, c, None, None] for c, k in enumerate(alg.k_factors)]
S_sec = [-(w * torch.log(w)).sum(-1) for w in (torch.linalg.eigvalsh(s) for s in sigma)]
H = -(p * torch.log(p)).sum(-1)
assert torch.allclose(S, H + sum(p[:, c] * S_sec[c] for c in range(alg.C)))    # Theorem 3.3

w_tau = torch.tensor([1 / (alg.C * k) for k in alg.k_factors], dtype=torch.float64)
d_tau = alg.operator(rho.density('tau_vN'), is_positive=True)
S_tau = -(d_tau @ d_tau.log()).tau_vN().real
assert torch.allclose(S_tau, S + (p * torch.log(w_tau)).sum(-1))
```

## 3.5 相对熵与 Klein 不等式

对密度为 $\rho,\sigma$ 且 $\sigma>0$ 的态，**相对熵**为
$S(\rho\|\sigma)=\operatorname{Tr}\rho\log\rho-\operatorname{Tr}\rho\log\sigma$，其中
$\operatorname{Tr}\rho\log\rho:=\sum_i\lambda_i\log\lambda_i$，求和遍及 $\rho$ 的特征值，约定
$0\log0=0$（在每个块中也同样如此）。它与迹的选取无关：由定理 3.2，
$(\rho_w)_c=\rho_c/w_c$，$(\sigma_w)_c=\sigma_c/w_c$，所以 $(\rho_w)_c$ 的特征值是 $\rho_c$ 的特征值除以 $w_c$，
$\log(\sigma_c/w_c)=\log\sigma_c-(\log w_c)1$（若 $\sigma_c=U\,\mathrm{diag}(s_i)\,U^*$，则
$\sigma_c/w_c=U\,\mathrm{diag}(s_i/w_c)\,U^*$，再用定理 1.5(a)），从而 $\operatorname{tr}_w(\rho_w\log\rho_w)$ 与 $\operatorname{tr}_w(\rho_w\log\sigma_w)$ 分别与
$\operatorname{Tr}\rho\log\rho$ 和 $\operatorname{Tr}\rho\log\sigma$ 相差同一个量 $-\sum_cp_c\log w_c$
（按约定 $0\log0=0$，$\rho$ 的对数只在其支撑上需要）。二者之差不变。

**引理 3.4（切线不等式）。** 设 $f:I\to\mathbb R$ 在区间 $I$ 上是凸的，且在
$b\in I$ 处可微（若 $b$ 是 $I$ 的端点，则指单侧可微）。则对所有 $a\in I$ 有 $f(a)\ge f(b)+(a-b)f'(b)$；
若 $f$ 严格凸，即对 $I$ 中的 $x\neq y$ 及 $t\in(0,1)$ 有
$f((1-t)x+ty)<(1-t)f(x)+tf(y)$，则当 $a\neq b$ 时不等式严格成立。此外，若 $f$ 在 $I$ 上连续，
且在 $I$ 的内部二阶可微并满足 $f''>0$，则 $f$ 严格凸。

*证明。* 对 $t\in(0,1]$，凸性给出 $f(b+t(a-b))\le(1-t)f(b)+tf(a)$，即
$\big(f(b+t(a-b))-f(b)\big)/t\le f(a)-f(b)$；令 $t\to0$ 即得该不等式。若 $f$
严格凸且 $a\neq b$，令 $m=(a+b)/2$；则
$f(a)-f(b)>2\big(f(m)-f(b)\big)\ge2(m-b)f'(b)=(a-b)f'(b)$，第一步用到在 $m$ 处的严格
凸性，第二步用到刚刚证明的不等式（应用于 $m$）。对于最后一个结论，设 $I$ 中 $x<y$，
$t\in(0,1)$，$z=(1-t)x+ty$。由中值定理，存在 $\xi_1\in(x,z)$ 与 $\xi_2\in(z,y)$，使得
$f(z)-f(x)=t(y-x)f'(\xi_1)$ 且 $f(y)-f(z)=(1-t)(y-x)f'(\xi_2)$；又因在 $[\xi_1,\xi_2]$ 上 $f''>0$，
对 $f'$ 在 $[\xi_1,\xi_2]$ 上应用中值定理得 $f'(\xi_1)<f'(\xi_2)$。乘以 $t(1-t)(y-x)>0$ 得
$(1-t)\big(f(z)-f(x)\big)<t\big(f(y)-f(z)\big)$，即 $f(z)<(1-t)f(x)+tf(y)$。$\square$

**定理 3.5（Klein 不等式）。** 设 $f:I\to\mathbb R$ 在区间 $I$ 上是凸的，
$X,Y$ 是谱含于 $I$ 的自伴算子，且 $f$ 在 $Y$ 的每个特征值处可微
（在 $I$ 的端点处指单侧可微）。
则 $\operatorname{Tr}[f(X)-f(Y)-(X-Y)f'(Y)]\ge0$。若 $f$ 严格凸，则等号
仅当 $X=Y$ 时成立。

*证明。* 设 $X=\sum_ia_i|u_i\rangle\langle u_i|$，$Y=\sum_jb_j|v_j\rangle\langle v_j|$，其中
两组特征基都是规范正交的，并令 $c_{ij}=|\langle u_i,v_j\rangle|^2$，于是 $\sum_ic_{ij}=\sum_jc_{ij}=1$。
则 $\operatorname{Tr}f(X)=\sum_{ij}c_{ij}f(a_i)$，$\operatorname{Tr}f(Y)=\sum_{ij}c_{ij}f(b_j)$，
$\operatorname{Tr}Xf'(Y)=\sum_{ij}c_{ij}a_if'(b_j)$，$\operatorname{Tr}Yf'(Y)=\sum_{ij}c_{ij}b_jf'(b_j)$。
因此左边等于 $\sum_{ij}c_{ij}[f(a_i)-f(b_j)-(a_i-b_j)f'(b_j)]$，由引理 3.4，这是
非负项之和。若 $f$ 严格凸且该和为零，则
只要 $c_{ij}>0$ 就有 $a_i=b_j$，所以对每个 $j$ 有 $Xv_j=\sum_ia_i\langle u_i,v_j\rangle u_i=b_jv_j$，
从而 $X=Y$。$\square$

**推论 3.6.** $S(\rho\|\sigma)\ge0$，等号成立当且仅当 $\rho=\sigma$。特别地
$S(\rho)\le\log D$，等号成立当且仅当 $\rho=1/D$。

*证明。* 对 $f(t)=t\log t$，$f(0)=0$ 应用定理 3.5；它在
$I=[0,\infty)$ 上严格凸（由引理 3.4：$f$ 在 $0$ 处连续，因为当 $t\to0^+$ 时 $t\log t\to0$，且
在 $(0,\infty)$ 上 $f''(t)=1/t>0$），并在
$(0,\infty)\supseteq\operatorname{spec}\sigma$ 上可微，$f'(t)=\log t+1$；取
$X=\rho$，$Y=\sigma$：左边为
$\operatorname{Tr}\rho\log\rho-\operatorname{Tr}\sigma\log\sigma-\operatorname{Tr}(\rho-\sigma)(\log\sigma+1)
=S(\rho\|\sigma)$，因为 $\operatorname{Tr}\rho=\operatorname{Tr}\sigma=1$。对于第二个结论，取
$\sigma=1/D$：$0\le S(\rho\|1/D)=-S(\rho)+\log D$。$\square$

**定理 3.7（链式法则；在到中心的条件期望下的单调性）。** 设 $\rho,\sigma$ 是态，
$\sigma>0$，$p_c=\operatorname{Tr}\rho_c$，$q_c=\operatorname{Tr}\sigma_c>0$，$\hat\sigma_c=\sigma_c/q_c$，并且
当 $p_c>0$ 时令 $\hat\rho_c=\rho_c/p_c$。则在对满足 $p_c>0$ 的 $c$ 求和的意义下（其余
扇区由 $0\log0=0$ 贡献为 $0$），
$$
S(\rho\|\sigma)=\sum_cp_c\log\frac{p_c}{q_c}+\sum_cp_c\,S(\hat\rho_c\|\hat\sigma_c)
\ \ge\ \sum_cp_c\log\frac{p_c}{q_c}=S(E\rho\|E\sigma),
$$
其中 $E$ 是到中心上的条件期望（命题 2.8）。

*证明。* 对 $p_c>0$，$\operatorname{Tr}\rho_c\log\rho_c=p_c\log p_c+p_c\operatorname{Tr}\hat\rho_c\log\hat\rho_c$
（特征值为 $p_c\mu$，同定理 3.3），且 $\log\sigma_c=(\log q_c)1+\log\hat\sigma_c$（将
$\hat\sigma_c=U\,\mathrm{diag}(s_i)\,U^*$ 对角化，则 $\sigma_c=U\,\mathrm{diag}(q_cs_i)\,U^*$，再用定理 1.5(a)），所以
$\operatorname{Tr}\rho_c\log\sigma_c=p_c\log q_c+p_c\operatorname{Tr}\hat\rho_c\log\hat\sigma_c$，从而
$\operatorname{Tr}\rho_c\log\rho_c-\operatorname{Tr}\rho_c\log\sigma_c=p_c\log(p_c/q_c)+p_cS(\hat\rho_c\|\hat\sigma_c)$。
若 $p_c=0$，则 $\rho_c=0$，该扇区贡献为 $0$。
对 $c$ 求和；由推论 3.6（在 $M_{k_c}$ 中应用于态 $\hat\rho_c$ 与 $\hat\sigma_c>0$），第二个和 $\ge0$。最后
$E\rho=\bigoplus_c(p_c/k_c)1_{k_c}$ 与 $E\sigma=\bigoplus_c(q_c/k_c)1_{k_c}>0$ 都是态；$E\rho$ 是态 $\omega\circ E$
（其中 $\omega(A)=\operatorname{Tr}\rho A$）的密度，因为
$\operatorname{Tr}(\rho E(A))=\sum_cp_c\operatorname{Tr}(A_c)/k_c=\operatorname{Tr}(E(\rho)A)$。并且
$S(E\rho\|E\sigma)=\sum_ck_c\frac{p_c}{k_c}\big(\log\frac{p_c}{k_c}-\log\frac{q_c}{k_c}\big)=\sum_cp_c\log\frac{p_c}{q_c}$，这里
$E\rho$ 在块 $c$ 中的特征值为 $p_c/k_c$（$k_c$ 重），且 $p_c=0$ 的扇区贡献为 $0$。$\square$

因此，把两个态都替换为它们在 $E$ 下的像只会减小相对熵；剩下的是各扇区分布的
Kullback–Leibler 散度 $\sum_cp_c\log(p_c/q_c)$。

> **物理意义**（物理解释，而非本课程的定理）。相对熵度量两个态能被区分的程度。忘掉扇区内部的
> 量子信息只会使两个态更难区分；剩下的是荷分布的经典可区分性。

```python
sigma_state = states.random_density_matrix(alg, batch_size=4)
q = sigma_state.sector_probabilities()
D_rel = states.relative_entropy(rho, sigma_state)
kl = (p * torch.log(p / q)).sum(-1)
E = channels.center_expectation(alg)
assert (D_rel >= kl - 1e-12).all()                                       # Theorem 3.7
assert torch.allclose(states.relative_entropy(E(rho), E(sigma_state)), kl)
assert (rho.von_neumann_entropy() <= math.log(10) + 1e-12).all()         # Corollary 3.6
assert abs(states.maximally_mixed_state(alg).von_neumann_entropy().item() - math.log(10)) < 1e-12
```

## 3.6 Michelson 对比度

对 $A\ge0$，$A\neq0$，**Michelson 对比度**为
$\Delta(A)=(\lambda_{\max}-\lambda_{\min})/(\lambda_{\max}+\lambda_{\min})\in[0,1]$，其中
极端特征值取遍所有扇区（并约定 $\Delta(0)=0$）。它度量 $A$ 偏离单位算子的倍数的程度：
$\Delta(A)=0$ 当且仅当 $A=\lambda1$。

**命题 3.8.** 若 $\Phi:\mathcal M\to\mathcal M$ 是线性、正且保单位的
（$\Phi(1)=1$），则对所有 $A\ge0$ 有 $\Delta(\Phi(A))\le\Delta(A)$。

*证明。* $\lambda_{\min}1\le A\le\lambda_{\max}1$（对 $A-\lambda_{\min}1$ 与 $\lambda_{\max}1-A$ 应用推论 1.4(3)），
由正性、线性和保单位性得
$\lambda_{\min}1\le\Phi(A)\le\lambda_{\max}1$。因此 $\Phi(A)$ 是自伴的（正矩阵都是自伴的，推论 1.4(3)），
$\Phi(A)\ge0$（因为 $\lambda_{\min}\ge0$），并且由推论 1.4(3) 得 $\lambda_{\min}(\Phi(A))\ge\lambda_{\min}$ 且
$\lambda_{\max}(\Phi(A))\le\lambda_{\max}$。若 $\Phi(A)=0$，结论是平凡的。否则 $A\neq0$，
$\lambda_{\max}$ 与 $\lambda_{\max}(\Phi(A))$ 都 $>0$，令 $r=\lambda_{\min}/\lambda_{\max}\in[0,1]$，则
$\Delta=(1-r)/(1+r)$ 关于 $r$ 递减，而 $r(\Phi(A))\ge r(A)$。$\square$

```python
A_pos = alg.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(4, d), batch_size=4, force_positive=True)
Phi = channels.random_mixed_unitary_channel(alg, 3, batch_size=4)       # unital
assert (Phi(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
assert (E(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
```

## 3.7 一个看似合理却不成立的不等式

对 $M_k$ 中的正矩阵 $X,Y$ 和酉矩阵 $U$，是否有 $|\operatorname{Tr}(XUY)|\le\operatorname{Tr}(XY)$？
它对 $U=1$ 成立，看起来像一个 Cauchy–Schwarz 型的命题。但只要 $X$
与 $Y$ 不可交换，它就不成立。

**命题 3.9.** 设 $X,Y\ge0$ 属于 $M_k$。则 $\max_U|\operatorname{Tr}(XUY)|=\|YX\|_1\ge\operatorname{Tr}(XY)$（最大值取遍酉矩阵
$U\in M_k$），等号成立当且仅当 $XY=YX$。

*证明。* $\operatorname{Tr}(XUY)=\operatorname{Tr}(YX\,U)$，所以由定理 1.9 最大值为 $\|YX\|_1$，
并且 $\|YX\|_1\ge|\operatorname{Tr}(YX)|=\operatorname{Tr}(XY)$，因为
$\operatorname{Tr}(XY)=\operatorname{Tr}(X^{1/2}YX^{1/2})\ge0$，这是由于
$X^{1/2}YX^{1/2}=(Y^{1/2}X^{1/2})^*(Y^{1/2}X^{1/2})\ge0$（推论 1.7）。若 $XY=YX$，则
$X^{1/2}$ 与 $Y$ 可交换（定理 1.5(d)），所以 $YX=X^{1/2}YX^{1/2}\ge0$，且
$\|YX\|_1=\operatorname{Tr}(YX)$（对 $P\ge0$，$\|P\|_1=\operatorname{Tr}P$，见 1.4 节）。反之，
设 $\|YX\|_1=\operatorname{Tr}(YX)$。若 $|YX|=0$，则 $YX=0=(YX)^*=XY$。否则写
$YX=W|YX|$（定理 1.8），并令 $R=|YX|^{1/2}\neq0$。由迹的循环性，
$\operatorname{Tr}(YX)=\operatorname{Tr}(RWR)=\langle R,WR\rangle_{HS}$，而
$\|R\|_2^2=\|WR\|_2^2=\operatorname{Tr}|YX|=\|YX\|_1$。所以 $\langle R,WR\rangle_{HS}=\|R\|_2\|WR\|_2$：
Cauchy–Schwarz 不等式取等号，因此对某个标量 $\lambda$ 有 $WR=\lambda R$，而
$\lambda\|R\|_2^2=\langle R,WR\rangle_{HS}=\|R\|_2^2$ 给出 $\lambda=1$。所以
$W|YX|^{1/2}=|YX|^{1/2}$，从而 $YX=W|YX|=|YX|\ge0$。正矩阵是自伴的，
所以 $YX=(YX)^*=XY$。$\square$

本库借助极分解求出最大值点：由定理 1.9，对 $YX=W|YX|$，酉矩阵 $U=W^*$ 达到 $\|YX\|_1$。
随机搜索平均而言更差：对按 Haar 分布的 $U$（定理 5.3；本库正是这样对 $U$ 采样的），取 $A=YX$ 应用命题 11.11 得
$\big(\mathbb E|\operatorname{Tr}(XUY)|^2\big)^{1/2}=\|YX\|_2/\sqrt k\le\|YX\|_1/\sqrt k=\max_V|\operatorname{Tr}(XVY)|/\sqrt k$
（因为 $\sum_is_i^2\le(\sum_is_i)^2$），而推论 11.12 表明，对谱含于 $[a,b]$（$a>0$）的 $X,Y$，单个样本达到
$\operatorname{Tr}(XY)$ 的概率至多为 $b^4/(a^4k^2)$（仅当 $k>b^2/a^2$ 时才有信息量；在下面的运行中 $b/a=6$ 且
$k=6$，所以它什么也没说）。代码对一对 $k=6$ 的不可交换矩阵计算精确的上确界以及 $20000$ 个 Haar
样本中的最佳值，二者都减去 $\operatorname{Tr}(XY)$。由命题 3.9，没有样本能超过上确界；最佳样本低于上确界
多少，是本次运行的*数值观察*，而不是定理。

```python
Mk = TypeIAlgebra([6], [6], precision='double')
X = Mk.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
Y = Mk.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
assert (X @ Y - Y @ X).matrix.abs().max().item() > 1e-6             # X and Y do not commute
YX = Y @ X
W = YX @ YX.abs().inverse()                                          # polar part of YX
gap = (X @ W.adjoint() @ Y).trace.abs() - (X @ Y).trace.real
assert abs(gap.item() - (YX.trace_norm() - (X @ Y).trace.real).item()) < 1e-10   # Theorem 1.9
assert gap.item() > 0                                  # Proposition 3.9: the inequality fails
samples = Mk.random_unitary_operator(20_000)
best = (X @ samples @ Y).trace.abs().max() - (X @ Y).trace.real
print(f"sup - Tr(XY) = {gap.item():.4f}; numerical observation: best of 20000 Haar samples "
      f"- Tr(XY) = {best.item():.4f}")
assert best.item() <= gap.item() + 1e-10               # Proposition 3.9: no sample beats the sup
```

## 3.8 可与之比较的工具

QuTiP（`entropy_vn`、`entropy_relative`）以及其他量子信息软件包（例如
`toqito`）用单个 Hilbert 空间的迹计算该空间上密度矩阵的熵与相对熵；在
$\bigoplus_c\mathbb C^{k_c}$ 上，这就是这里的 blunt 迹（至多相差对数的底，不同软件包的取法不同）。它们不提供的是定理 3.1
中的迹族以及定理 3.2 中的换算，而一旦代数以一定重数作用，或者各扇区必须赋予不同权重，
这些就变得必不可少。

## 习题

1. 证明迹态 $\tau_{\mathrm{vN}}$（各扇区权重相等）关于
   $\tau_{\mathrm{vN}}$ 的密度为 $1$，并计算它的 blunt 密度及其熵 $S$。
2. 设 $\operatorname{tr}_w$ 忠实。证明 $S_w(\rho)\le\log\operatorname{tr}_w(1)$，等号成立当且仅当
   $\rho_w=1/\operatorname{tr}_w(1)$。
3. 证明 $M_D$（而不是
   $\mathcal M$！）上的测量信道 $\rho\mapsto\sum_cP_c\rho P_c$ 确实会扰动态，并对两个扇区的
   叠加态计算测量前后的 $S(\rho)$。

<details>
<summary>解答</summary>

1. $\tau_{\mathrm{vN}}(A)=\tau_{\mathrm{vN}}(1\cdot A)$，由唯一性得 $\rho_{\tau}=1$。由
   定理 3.2，blunt 密度为 $\bigoplus_c1_{k_c}/(Ck_c)$，其中 $p_c=1/C$，
   $\sigma_c=1/k_c$，所以由定理 3.3 得 $S=\log C+\frac1C\sum_c\log k_c$。

```python
tau = states.tracial_state(alg)
assert (tau.density('tau_vN') - one.matrix).abs().max() < 1e-12
assert abs(tau.von_neumann_entropy().item() - (math.log(3) + sum(math.log(k) for k in (2, 3, 5)) / 3)) < 1e-12
```

2. 求和遍及满足 $p_c>0$ 的 $c$。则
   $S_w(\rho)=S(\rho)+\sum_cp_c\log w_c=-\sum_cp_c\log p_c+\sum_cp_cS(\sigma_c)+\sum_cp_c\log w_c
   \le\sum_cp_c\log\frac{k_cw_c}{p_c}\le\log\sum_{p_c>0}k_cw_c\le\log\operatorname{tr}_w(1)$，这里用到
   $S(\sigma_c)\le\log k_c$（在 $M_{k_c}$ 中应用推论 3.6）以及有限 Jensen 不等式
   $\sum_{c=1}^mp_c\log x_c\le\log\sum_{c=1}^mp_cx_c$（$p_c>0$，$\sum_cp_c=1$，$x_c>0$），等号仅当所有 $x_c$
   相等时成立，这里应用于 $x_c=k_cw_c/p_c$。对 $m$ 归纳证明：$m=1$ 时无需证明。$m\ge2$ 时令
   $q=1-p_m\in(0,1)$，$\bar x=\sum_{c<m}(p_c/q)x_c$；由归纳假设
   $\sum_{c<m}(p_c/q)\log x_c\le\log\bar x$，等号仅当 $x_1=\dots=x_{m-1}$ 时成立；由 $\log$ 的严格凹性
   （对 $-\log$ 应用引理 3.4，其二阶导数为 $1/t^2>0$）
   $q\log\bar x+p_m\log x_m\le\log(q\bar x+p_mx_m)$，
   等号仅当 $\bar x=x_m$ 时成立。将第一个不等式乘以 $q$ 后相加，得
   $\sum_cp_c\log x_c\le\log\sum_cp_cx_c$，而等号迫使 $x_1=\dots=x_{m-1}=\bar x=x_m$。
   各处都取等号
   迫使每个 $p_c>0$（因为所有 $k_cw_c>0$）、$\sigma_c=1/k_c$ 且 $k_cw_c/p_c$ 为常数，即
   $p_c=k_cw_c/\operatorname{tr}_w(1)$，亦即 $(\rho_w)_c=\rho_c/w_c=1_{k_c}/\operatorname{tr}_w(1)$。
   反之，对 $\rho_w=1/\operatorname{tr}_w(1)$，
   $S_w(\rho)=-\operatorname{tr}_w(1)\cdot\frac1{\operatorname{tr}_w(1)}\log\frac1{\operatorname{tr}_w(1)}=\log\operatorname{tr}_w(1)$。

3. 在 $M_D$ 中，由两个扇区中的单位向量构成的纯叠加态 $\psi=(\psi_1\oplus\psi_2)/\sqrt2$
   满足 $S=0$；测量之后的态为
   $\frac12|\psi_1\rangle\langle\psi_1|\oplus\frac12|\psi_2\rangle\langle\psi_2|$，$S=\log2$。然而，限制到
   $\mathcal M$ 上时，$\psi$ 的向量态与该混合态一致：对所有 $A\in\mathcal M$，
   $\langle\psi,A\psi\rangle=\tfrac12\langle\psi_1,A_1\psi_1\rangle+\tfrac12\langle\psi_2,A_2\psi_2\rangle$。
   $\mathcal M$ 中没有元素能探测到这种相干性（在物理上，这正是
   超选择定则的内容）。

</details>

**下一章：** [第 4 章](04_symmetries.md)讨论对称性与表示：物理系统上的群作用如何产生扇区，以及本库
如何为自旋系统和量子比特寄存器构造相应的代数。
