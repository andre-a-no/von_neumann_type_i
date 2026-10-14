[English](../D_gksl.md) · [Español](../es/D_gksl.md) · [Français](../fr/D_gksl.md) · [Deutsch](../de/D_gksl.md) · **中文** · [日本語](../ja/D_gksl.md) · [Русский](../ru/D_gksl.md)

# 附录 D. 量子动力学半群的生成元

第 8 章证明了 Gorini–Kossakowski–Sudarshan–Lindblad（GKSL）形式的生成元产生量子动力学半群
（定理 8.6），并使用了其逆命题（定理 8.7）以及它在 $\mathcal M=\bigoplus_cM_{k_c}$ 上的版本。本附录在
有限维情形下证明这两者，同时给出 GKSL 形式中的自由度以及稳态唯一性的一个判据。所用的工具是
第 6 章的 Choi 矩阵：生成元一般不是完全正的，但它的 Choi 矩阵在最大纠缠向量的正交补上是正的。
同一类论证——先证明一族算子的正性，再作分解——在最后一节证明了第 10 章所陈述的算子代数纠错定理
（定理 10.7），并给出一个显式的恢复信道。

**你将学到**为什么完全正半群的生成元是条件完全正的，如何从生成元的 Choi 矩阵读出 $H$ 和跳跃算子，
它们在什么意义下是唯一的，为什么有扇区的代数上的每个生成元都具有 8.3 节所用的形式，稳态何时唯一，
以及恢复信道能够保护哪些可观测量代数。

全文中，$M_n=M_n(\mathbb C)$，$E_{ij}=|i\rangle\langle j|$ 是它的矩阵单位，**超算子**是 $M_n$（或
$\mathcal M$）上的线性映射，而
$$
\Omega=\sum_{i=1}^ne_i\otimes e_i\in\mathbb C^n\otimes\mathbb C^n,\qquad
J(\Psi)=\sum_{i,j=1}^nE_{ij}\otimes\Psi(E_{ij})=[\Psi(E_{ij})]_{i,j=1}^n
$$
分别是最大纠缠向量和超算子 $\Psi$ 的 **Choi 矩阵**（定理 6.2；这里 $\Psi$ 不必是正的）。对向量
$v,w$，$|v\rangle\langle w|$ 表示算子 $x\mapsto\langle w,x\rangle v$。若对所有 $A$ 有
$\Psi(A^*)=\Psi(A)^*$，则称超算子 $\Psi$ 是**保厄米的**。量子动力学半群及其生成元是 8.2 节
（引理 8.5）中的那些：$T_t=e^{t\mathcal L}$，其中 $T_t$ 是信道。

## D.1 基本构件的 Choi 矩阵

**引理 D.1.** 设 $\Psi$ 是 $M_n$ 上的超算子。

(a) $\Psi\mapsto J(\Psi)$ 是从 $M_n$ 上的超算子到 $M_n(M_n)=M_{n^2}$ 上的线性双射。
(b) 每个 $v\in\mathbb C^n\otimes\mathbb C^n$ 都可以对恰好一个 $K\in M_n$ 写成 $v=(1\otimes K)\Omega$，并且
$\langle(1\otimes K)\Omega,(1\otimes L)\Omega\rangle=\operatorname{Tr}(K^*L)$。特别地，
$\langle\Omega,(1\otimes K)\Omega\rangle=\operatorname{Tr}K$。
(c) 对 $K,L\in M_n$，超算子 $A\mapsto KAL^*$ 的 Choi 矩阵是
$|(1\otimes K)\Omega\rangle\langle(1\otimes L)\Omega|$。特别地，$J(\mathrm{id})=|\Omega\rangle\langle\Omega|$。
(d) $\Psi$ 是保厄米的当且仅当 $J(\Psi)$ 是自伴的。每个 CP 映射都是保厄米的。

*证明。* (a) $J(\Psi)$ 的 $(i,j)$ 块是 $\Psi(E_{ij})$，而各 $E_{ij}$ 构成 $M_n$ 的一组基，所以
$J$ 是线性单射，并且每个块矩阵 $[X_{ij}]$ 都是把 $E_{ij}$ 映为 $X_{ij}$ 的超算子的 Choi 矩阵。
(b) $(1\otimes K)\Omega=\sum_ie_i\otimes Ke_i$ 的第 $i$ 块是 $Ke_i$，即 $K$ 的第 $i$ 列，它们决定了 $K$，
并且每个 $v$ 都以这种方式出现。于是
$\langle\sum_ie_i\otimes Ke_i,\sum_je_j\otimes Le_j\rangle=\sum_i\langle Ke_i,Le_i\rangle=\operatorname{Tr}(K^*L)$；
对最后一个断言，把 $(K,L)$ 换成 $(1,K)$。(c) 右端的 $(i,j)$ 块是
$Ke_i(Le_j)^*=KE_{ij}L^*$。(d) $J(\Psi)^*$ 的 $(i,j)$ 块是 $\Psi(E_{ji})^*$。若 $\Psi$ 是保厄米的，则
$\Psi(E_{ji})^*=\Psi(E_{ji}^*)=\Psi(E_{ij})$。反之，若对所有 $i,j$ 有 $\Psi(E_{ji})^*=\Psi(E_{ij})$，且
$A=\sum a_{ij}E_{ij}$，则
$\Psi(A)^*=\sum\bar a_{ij}\Psi(E_{ij})^*=\sum\bar a_{ij}\Psi(E_{ji})=\Psi(A^*)$。CP 映射具有 Kraus 形式
（定理 6.2），而 $(\sum_lK_lAK_l^*)^*=\sum_lK_lA^*K_l^*$。$\square$

**引理 D.2（分解）。** 设 $q_1,\dots,q_C$ 是 $\mathbb C^n\otimes\mathbb C^n$ 中的规范正交向量，写成
$q_c=(1\otimes Q_c)\Omega$，令 $W=\operatorname{span}\{q_c\}$，并设 $\Psi$ 是保厄米的超算子，且对所有
$v\perp W$ 有 $\langle v,J(\Psi)v\rangle\ge0$。令 $J=J(\Psi)$，$\xi_c=Jq_c-\frac12\sum_d\langle q_d,Jq_c\rangle q_d$，
并由 $(1\otimes K_c)\Omega=\xi_c$ 定义 $K_c\in M_n$。则存在 $L_1,\dots,L_r\in M_n$，
$r=\operatorname{rank}(P_W^\perp JP_W^\perp)$，其中 $P_W^\perp$ 是到 $W^\perp$ 上的投影，使得
$$
\Psi(A)=\sum_{j=1}^rL_jAL_j^*+\sum_{c=1}^C\big(K_cAQ_c^*+Q_cAK_c^*\big)\qquad(A\in M_n),
$$
对所有 $c,j$ 有 $\operatorname{Tr}(Q_c^*L_j)=0$，并且当 $j\neq k$ 时 $\operatorname{Tr}(L_j^*L_k)=0$。

*证明。* 由引理 D.1(d)，$J$ 是自伴的。令 $P=\sum_c|q_c\rangle\langle q_c|$ 为到 $W$ 上的投影，
$P^\perp=P_W^\perp=1-P$。展开 $J=(P+P^\perp)J(P+P^\perp)$，并利用 $PJP^\perp=PJ-PJP$ 与
$P^\perp JP=JP-PJP$，得到
$$
J=P^\perp JP^\perp+\big(PJ+JP-PJP\big).
$$
第一项是自伴的，并且由于 $P^\perp v\perp W$，有
$\langle v,P^\perp JP^\perp v\rangle=\langle P^\perp v,JP^\perp v\rangle\ge0$；所以 $P^\perp JP^\perp\ge0$。由定理 1.3，
$P^\perp JP^\perp=\sum_{j=1}^r|l_j\rangle\langle l_j|$，其中 $l_j$ 是非零特征值 $\mu_j$ 的特征向量乘以
$\sqrt{\mu_j}$：它们两两正交，并且位于 $P^\perp JP^\perp$ 的值域中，而该值域与 $W$ 正交。写
$l_j=(1\otimes L_j)\Omega$；由引理 D.1(b)，$\operatorname{Tr}(Q_c^*L_j)=\langle q_c,l_j\rangle=0$，并且当 $j\ne k$ 时
$\operatorname{Tr}(L_j^*L_k)=\langle l_j,l_k\rangle=0$；由引理 D.1(c)，$P^\perp JP^\perp$ 是
$A\mapsto\sum_jL_jAL_j^*$ 的 Choi 矩阵。对第二项，
$JP=\sum_c|Jq_c\rangle\langle q_c|$，$PJP=\sum_{c,d}\langle q_d,Jq_c\rangle|q_d\rangle\langle q_c|$，所以
$$
\sum_c|\xi_c\rangle\langle q_c|=JP-\tfrac12PJP,\qquad\sum_c|q_c\rangle\langle\xi_c|=PJ-\tfrac12PJP
$$
（第二式是第一式的伴随；$PJP$ 是自伴的），二者之和为 $PJ+JP-PJP$。由引理 D.1(c)，
$|\xi_c\rangle\langle q_c|$ 与 $|q_c\rangle\langle\xi_c|$ 分别是 $A\mapsto K_cAQ_c^*$ 与 $A\mapsto Q_cAK_c^*$ 的
Choi 矩阵。所以所断言的恒等式两边具有相同的 Choi 矩阵 $J$，再由引理 D.1(a) 即得该恒等式。$\square$

## D.2 条件完全正性

若对每个 $v\perp\Omega$ 有 $\langle v,J(\Psi)v\rangle\ge0$，则称 $M_n$ 上的超算子 $\Psi$ 是
**条件完全正的**（CCP）。由引理 D.1(b)，这些 $v$ 恰好是满足 $\operatorname{Tr}K=0$ 的向量
$(1\otimes K)\Omega$。

**定理 D.3.** 对 $M_n$ 上的超算子 $\mathcal L$，以下各条等价：

1. 对每个 $t\ge0$，$e^{t\mathcal L}$ 是 CP 的；
2. $\mathcal L$ 是保厄米且条件完全正的；
3. 存在 $M_n$ 上的 CP 映射 $\Phi$ 和 $G\in M_n$，使得对所有 $A$ 有 $\mathcal L(A)=\Phi(A)+GA+AG^*$。

在 3 中可以取 $\Phi(A)=\sum_{j=1}^rL_jAL_j^*$，其中 $r\le n^2-1$，$\operatorname{Tr}L_j=0$，且当 $j\neq k$ 时
$\operatorname{Tr}(L_j^*L_k)=0$，并且可以取 $\operatorname{Tr}G\in\mathbb R$ 的 $G$。

*证明。* (1)$\Rightarrow$(2)。由引理 8.1，当 $t\to0$ 时 $(e^{t\mathcal L}-\mathrm{id})/t\to\mathcal L$，而 $J$ 是
有限维空间上的线性映射，因而连续，所以
$J(e^{t\mathcal L})=J(\mathrm{id})+t\,J\big((e^{t\mathcal L}-\mathrm{id})/t\big)$，其中最后一个 Choi 矩阵趋于
$J(\mathcal L)$。对 $v\perp\Omega$，$\langle v,J(\mathrm{id})v\rangle=|\langle\Omega,v\rangle|^2=0$（引理 D.1(c)），所以
对 $t>0$，由于 $J(e^{t\mathcal L})\ge0$（定理 6.2），
$$
0\le\frac1t\langle v,J(e^{t\mathcal L})v\rangle=\big\langle v,J\big((e^{t\mathcal L}-\mathrm{id})/t\big)v\big\rangle
\;\xrightarrow[t\to0]{}\;\langle v,J(\mathcal L)v\rangle .
$$
每个 $e^{t\mathcal L}$ 都是 CP 的，从而是保厄米的（引理 D.1(d)），所以
$\mathcal L(A^*)=\lim_t(e^{t\mathcal L}(A)^*-A^*)/t=\mathcal L(A)^*$。

(2)$\Rightarrow$(3)。对 $C=1$ 和 $q_1=\Omega/\sqrt n$（即 $Q_1=1/\sqrt n$）应用引理 D.2：
$\mathcal L(A)=\sum_jL_jAL_j^*+GA+AG^*$，其中 $G=K_1/\sqrt n$，而 $\Phi=\sum_jL_j\cdot L_j^*$ 是 CP 的
（定理 6.2）。此外 $\operatorname{Tr}L_j=\sqrt n\operatorname{Tr}(Q_1^*L_j)=0$，$r\le\operatorname{rank}P_W^\perp=n^2-1$，
并且由引理 D.1(b)，
$\operatorname{Tr}G=\langle\Omega,\xi_1\rangle/\sqrt n=\langle q_1,\xi_1\rangle=\frac12\langle q_1,J(\mathcal L)q_1\rangle\in\mathbb R$，
因为 $J(\mathcal L)$ 是自伴的。

(3)$\Rightarrow$(1)。定理 8.6 证明的第 1–4 步只用到生成元为
$\mathcal L_0+\Phi$，其中对某个矩阵 $H_{\mathrm{eff}}$ 有 $\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)$，
且 $\Phi$ 是 CP 的。取 $H_{\mathrm{eff}}=iG$，则 $\mathcal L_0(\rho)=G\rho+\rho G^*$，而这几步表明
$e^{t\mathcal L}=\sum_{m\ge0}S_m(t)$ 是 CP 的。$\square$

仅有正性是不够的。设 $n\ge2$，$\theta$ 是 $M_n$ 上的转置，$\mathcal L=\theta-\mathrm{id}$。由于 $\theta^2=\mathrm{id}$，
$$
e^{t\mathcal L}=e^{-t}\big(\cosh t\,\mathrm{id}+\sinh t\,\theta\big)=\tfrac{1+e^{-2t}}2\,\mathrm{id}+\tfrac{1-e^{-2t}}2\,\theta ,
$$
这是正的保迹映射的凸组合（命题 6.1）。但 $J(\theta)=[E_{ji}]=\sum_{i,j}E_{ij}\otimes E_{ji}$ 是交换算子
$e_k\otimes e_l\mapsto e_l\otimes e_k$，所以它把与 $\Omega$ 正交的
$v=e_1\otimes e_2-e_2\otimes e_1$ 映为 $-v$（同命题 6.1 的证明）；由
$\|v\|^2=2$ 与 $\langle v,J(\mathrm{id})v\rangle=|\langle\Omega,v\rangle|^2=0$ 得
$\langle v,J(\mathcal L)v\rangle=-2$：$\mathcal L$ 不是 CCP 的，而且确实对每个 $t>0$ 有
$\langle v,J(e^{t\mathcal L})v\rangle=-(1-e^{-2t})<0$。由正的保迹映射构成的半群
不一定由信道构成。

我们在 $M_3$ 上一个带任意跳跃的 GKSL 生成元上检验定理 D.3，它由 `lindblad_superoperator` 构造
（即 $\mathcal L$ 作用在按行展开的 $\operatorname{vec}\rho$ 上的矩阵 $S$，因此与 6.1 节一样
$J_{(i,a),(j,b)}=S_{(a,b),(i,j)}$），并且也在 $\theta-\mathrm{id}$ 上检验：

```python
import math
import numpy as np
import scipy.linalg as sla
import torch
from torch_vn_algebra import TypeIAlgebra, channels, dynamics

torch.manual_seed(0)
rng = np.random.default_rng(0)
f64, c128 = torch.float64, torch.complex128

def choi_of(S, n):
    """J_{(i,a),(j,b)} = S_{(a,b),(i,j)}: the Choi matrix of the map with row-major matrix S (Section 6.1)."""
    return S.reshape(n, n, n, n).transpose(2, 0, 3, 1).reshape(n * n, n * n)

def superop(f, n):
    """Row-major matrix of a linear map f on n x n matrices."""
    E = np.eye(n * n).reshape(n * n, n, n)
    return np.stack([f(X).reshape(-1) for X in E], axis=1)

vec_of = lambda K: K.T.reshape(-1)                        # (1 (x) K) Omega = sum_i e_i (x) K e_i
op_of = lambda v, n: v.reshape(n, n).T                    # its inverse (Lemma D.1(b))
cplx = lambda *s: rng.standard_normal(s) + 1j * rng.standard_normal(s)

n = 3
M3 = TypeIAlgebra([n], [n], precision='double')
H = cplx(n, n)
H = (H + H.conj().T) / 2
Ls = [cplx(n, n) / 2 for _ in range(2)]                   # arbitrary jumps, not traceless
to_op = lambda X: M3.operator(torch.from_numpy(X).reshape(1, 1, n, n))
S = dynamics.lindblad_superoperator(M3, to_op(H), [to_op(L) for L in Ls])[0][0].numpy()
J = choi_of(S, n)
q = vec_of(np.eye(n)) / math.sqrt(n)                      # the unit vector Omega / sqrt(n)
Pperp = np.eye(n * n) - np.outer(q, q)
assert np.abs(J - J.conj().T).max() < 1e-12               # Hermiticity preserving (Lemma D.1(d))
print(f"min eigenvalue of J(L): {np.linalg.eigvalsh(J).min():.2f}; "
      f"on the complement of Omega: {np.linalg.eigvalsh(Pperp @ J @ Pperp).min():.1e}")
assert np.linalg.eigvalsh(J).min() < -0.1                 # L itself is not CP ...
assert np.linalg.eigvalsh(Pperp @ J @ Pperp).min() > -1e-12        # ... but CCP (Theorem D.3)
for t in (0.01, 0.1, 1.0):
    assert np.linalg.eigvalsh(choi_of(sla.expm(t * S), n)).min() > -1e-12

L_theta = superop(lambda X: X.T, n) - np.eye(n * n)     # L = transpose - id: Hermiticity preserving, not CCP
anti = vec_of(np.eye(n)[:, [1]] @ np.eye(n)[[0], :] - np.eye(n)[:, [0]] @ np.eye(n)[[1], :])   # e1(x)e2 - e2(x)e1
assert abs(anti @ q) < 1e-15 and abs(anti @ choi_of(L_theta, n) @ anti + 2) < 1e-12
X = cplx(n, n)
for t in (0.01, 0.1, 1.0):
    T_t = sla.expm(t * L_theta)
    b = (1 - math.exp(-2 * t)) / 2
    assert abs(anti @ choi_of(T_t, n) @ anti + (1 - math.exp(-2 * t))) < 1e-12        # <v, J(e^{tL}) v>
    assert np.abs(T_t - ((1 - b) * np.eye(n * n) + b * superop(lambda Y: Y.T, n))).max() < 1e-12
    assert abs(np.trace((T_t @ X.reshape(-1)).reshape(n, n)) - np.trace(X)) < 1e-12   # trace preserving
    assert np.linalg.eigvalsh(choi_of(T_t, n)).min() < -0.9 * b                     # but not CP
```

## D.3 GKSL 定理

**定理 D.4（Gorini–Kossakowski–Sudarshan；Lindblad；定理 8.7）。** $M_n$ 上的超算子 $\mathcal L$ 是
$M_n$ 上某个量子动力学半群的生成元，当且仅当存在 $M_n$ 中的 $H=H^*$ 与 $L_1,\dots,L_r$，使得
$$
\mathcal L(\rho)=-i[H,\rho]+\sum_{j=1}^r\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
$$
可以选取 $r\le n^2-1$，$\operatorname{Tr}H=0$，$\operatorname{Tr}L_j=0$，并且当 $j\neq k$ 时 $\operatorname{Tr}(L_j^*L_k)=0$。

*证明。* "当"的方向就是定理 8.6。反之，设 $(T_t)$ 是以 $\mathcal L$ 为生成元的量子动力学半群，
于是 $T_t=e^{t\mathcal L}$（引理 8.5）。由于每个 $T_t$ 都是 CP 的，定理 D.3 给出
$\mathcal L(\rho)=\sum_jL_j\rho L_j^*+G\rho+\rho G^*$，其中 $r\le n^2-1$ 个 $L_j$ 无迹且两两正交，并且
$\operatorname{Tr}G\in\mathbb R$。令 $H=\frac i2(G-G^*)$，$R=\frac12(G+G^*)$，二者都是自伴的；于是
$G=-iH+R$，$G^*=iH+R$，且 $G\rho+\rho G^*=-i[H,\rho]+\{R,\rho\}$。此外
$\operatorname{Tr}H=\frac i2(\operatorname{Tr}G-\overline{\operatorname{Tr}G})=0$。由于每个 $T_t$ 都是 TP 的，
$\operatorname{Tr}\mathcal L(\rho)=\frac d{dt}\operatorname{Tr}T_t(\rho)\big|_{t=0}=0$（引理 8.1），再由循环性，
$$
0=\operatorname{Tr}\mathcal L(\rho)=\operatorname{Tr}\Big(\rho\,\Big(\sum_jL_j^*L_j+2R\Big)\Big)\qquad\text{for all }\rho .
$$
取 $\rho=X=\sum_jL_j^*L_j+2R$，它是自伴的，得 $\operatorname{Tr}(X^*X)=0$，所以 $X=0$，
$R=-\frac12\sum_jL_j^*L_j$，这就是 GKSL 形式。$\square$

**定理 D.5（GKSL 形式中的自由度）。** 设 $(H;L_1,\dots,L_r)$ 与 $(H';L'_1,\dots,L'_s)$ 是 GKSL
数据（$H,H'$ 自伴），按定理 D.4 给出生成元 $\mathcal L$、$\mathcal L'$。

(a) 在以下每种情形下都有 $\mathcal L'=\mathcal L$：
(i) $L'_j=L_j$，且 $H'=H+c1$，$c\in\mathbb R$；
(ii) $H'=H$，并且在用零补齐较短的那一族之后，$L'_j=\sum_ku_{jk}L_k$，其中 $u$ 是酉矩阵；
(iii) $L'_j=L_j+a_j1$，$a_j\in\mathbb C$，且 $H'=H-\frac i2\sum_j(\bar a_jL_j-a_jL_j^*)$，它是自伴的。

(b) 若 $\operatorname{Tr}H=\operatorname{Tr}H'=0$，所有 $L_j,L'_j$ 都无迹，且 $\mathcal L=\mathcal L'$，则 $H=H'$，并且
在补齐之后 $L'_j=\sum_ku_{jk}L_k$，其中 $u$ 是酉矩阵。

(c) 因此，两组 GKSL 数据具有相同的生成元，当且仅当其中一组可以由另一组经过有限多步 (i)–(iii)
得到。特别地，定理 D.4 的归一化形式中的算子 $H$ 及其耗散部分 $\Phi(\rho)=\sum_jL_j\rho L_j^*$ 由
$\mathcal L$ 唯一确定。

*证明。* (a) (i) $[c1,\rho]=0$。(ii) 由定理 6.7，$\Phi'=\sum_jL'_j\cdot L_j'^*$ 等于
$\Phi=\sum_jL_j\cdot L_j^*$；对偶不依赖于 Kraus 算子的选取（6.2 节），所以
$\sum_jL_j'^*L'_j=\Phi^*(1)=\sum_jL_j^*L_j$。(iii) 对单个跳跃 $L$ 和 $a\in\mathbb C$，展开得到
$$
(L+a)\rho(L+a)^*-\tfrac12\{(L+a)^*(L+a),\rho\}-\Big(L\rho L^*-\tfrac12\{L^*L,\rho\}\Big)=\tfrac12[B,\rho],
\qquad B=\bar aL-aL^* ,
$$
因为第一个乘积中的项 $\bar aL\rho+a\rho L^*$ 与
$-\frac12(aL^*\rho+a\rho L^*+\bar aL\rho+\bar a\rho L)$ 合并，而含 $|a|^2$ 的项相互抵消。由于 $B^*=-B$，
矩阵 $\frac i2B$ 是自伴的，所以 $H'=H-\frac i2\sum_jB_j$ 是自伴的，并且
$-i[H'-H,\rho]=-\frac12\sum_j[B_j,\rho]$ 抵消了耗散项的变化。

(b) 令 $\Phi=\sum_jL_j\cdot L_j^*$，$G=-iH-\frac12\Phi^*(1)$，于是 $\mathcal L(\rho)=\Phi(\rho)+G\rho+\rho G^*$。
由引理 D.1(c)，
$J(\mathcal L)=\sum_j|l_j\rangle\langle l_j|+|g\rangle\langle\Omega|+|\Omega\rangle\langle g|$，其中 $l_j=(1\otimes L_j)\Omega$，
$g=(1\otimes G)\Omega$，并且由于 $\operatorname{Tr}L_j=0$ 有 $l_j\perp\Omega$（引理 D.1(b)）。记 $P^\perp$ 为到
$\Omega^\perp$ 上的投影，于是 $P^\perp J(\mathcal L)P^\perp=\sum_j|l_j\rangle\langle l_j|=J(\Phi)$，同样
$P^\perp J(\mathcal L')P^\perp=J(\Phi')$。所以 $J(\Phi)=J(\Phi')$，$\Phi=\Phi'$（引理 D.1(a)），再由定理 6.7 得到
酉矩阵 $u$。于是 $\Phi^*(1)=\Phi'^*(1)$，从而对所有 $\rho$ 有 $[H-H',\rho]=0$。与所有 $E_{ij}$ 可交换的矩阵
$X$ 是 $1$ 的倍数（比较 $XE_{ij}$ 与 $E_{ij}X$ 的矩阵元），所以 $H-H'=\lambda1$，取迹即得
$\lambda=0$。

(c) 由 (a)，步骤 (i)–(iii) 不改变生成元，并且每一步的逆仍是同一种步骤（对 (iii) 取 $-a_j$，因为
$\bar a_jL'_j-a_jL_j'^*=\bar a_jL_j-a_jL_j^*$）。任何数据都可以先经 $a_j=-\operatorname{Tr}L_j/n$ 的 (iii)，
再经 $c=-\operatorname{Tr}H''/n$ 的 (i) 化为归一化形式，其中 $H''$ 是 (iii) 之后的哈密顿量；由 (b)，
具有相同生成元的两组归一化数据由 (ii) 联系起来。$\square$

定理 D.4 的证明就是一个算法：由 $J(\mathcal L)$ 在 $\Omega^\perp$ 上的压缩得到 $\Phi$，再由 $\xi_1$ 得到 $G$，
然后得到 $H$。我们对上面的生成元运行这个算法，用本库重建生成元，并与按定理 D.5(iii) 和 (i) 计算出的
原始数据的归一化形式进行比较：

```python
def gksl_from_generator(S, n, tol=1e-10):
    """(H, [L_j]) with Tr H = Tr L_j = 0 from the matrix of a GKSL generator: the proof of Theorem D.4."""
    J = choi_of(S, n)
    q = vec_of(np.eye(n)) / math.sqrt(n)
    Pperp = np.eye(n * n) - np.outer(q, q)
    w, U = np.linalg.eigh(Pperp @ J @ Pperp)              # = J(Phi) >= 0
    jumps = [math.sqrt(x) * op_of(U[:, k], n) for k, x in enumerate(w) if x > tol]
    xi = J @ q - 0.5 * (q.conj() @ J @ q) * q
    G = op_of(xi, n) / math.sqrt(n)
    return 0.5j * (G - G.conj().T), jumps

H_rec, L_rec = gksl_from_generator(S, n)
S_rec = dynamics.lindblad_superoperator(M3, to_op(H_rec), [to_op(L) for L in L_rec])[0][0].numpy()
assert np.abs(S_rec - S).max() < 1e-12                    # the same generator
assert abs(np.trace(H_rec)) < 1e-12 and all(abs(np.trace(L)) < 1e-12 for L in L_rec)

a = [-np.trace(L) / n for L in Ls]                        # Theorem D.5(iii): make the jumps traceless
L_norm = [L + aj * np.eye(n) for L, aj in zip(Ls, a)]
H_norm = H - 0.5j * sum(np.conj(aj) * L - aj * L.conj().T for L, aj in zip(Ls, a))
H_norm = H_norm - np.trace(H_norm) / n * np.eye(n)        # Theorem D.5(i)
assert np.abs(H_norm - H_norm.conj().T).max() < 1e-12
assert np.abs(H_rec - H_norm).max() < 1e-12               # H is unique (Theorem D.5(b))
jump_choi = lambda Ks: sum(np.outer(vec_of(K), vec_of(K).conj()) for K in Ks)
assert np.abs(jump_choi(L_rec) - jump_choi(L_norm)).max() < 1e-12      # so is sum_j L_j . L_j^*
Wr, Wn = np.stack([vec_of(K) for K in L_rec], 1), np.stack([vec_of(K) for K in L_norm], 1)
u = np.linalg.lstsq(Wr, Wn, rcond=None)[0]                # L_norm_j = sum_k u_kj L_rec_k
assert np.abs(Wr @ u - Wn).max() < 1e-12 and np.abs(u.conj().T @ u - np.eye(2)).max() < 1e-12   # unitary
print("eigenvalues of the normalised H:", np.round(np.linalg.eigvalsh(H_rec), 4), " jumps:", len(L_rec))
```

**定理 D.6（稳态的唯一性）。** 设 $\mathcal L$ 是定理 D.4 中以 $(H;L_1,\dots,L_r)$ 为数据的 GKSL
生成元。若 $\mathcal L$ 有一个忠实的稳态 $\sigma$（$\sigma>0$），并且与 $H$ 以及所有 $L_j$、$L_j^*$ 可交换的
矩阵只有 $1$ 的倍数，则 $\ker\mathcal L=\mathbb C\sigma$：$\sigma$ 是唯一的稳态。此外 $\ker\mathcal L^*=\mathbb C1$。

*证明。* 关于迹配对 $\operatorname{Tr}(\mathcal L(\rho)A)=\operatorname{Tr}(\rho\,\mathcal L^*(A))$ 的对偶是
$\mathcal L^*(A)=i[H,A]+\sum_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})$：由循环性，
$\operatorname{Tr}(-i[H,\rho]A)=\operatorname{Tr}(\rho\,i[H,A])$，$\operatorname{Tr}(L\rho L^*A)=\operatorname{Tr}(\rho L^*AL)$，
$\operatorname{Tr}(\{L^*L,\rho\}A)=\operatorname{Tr}(\rho\{L^*L,A\})$（这是第 8 章的习题 1）。
直接展开（其中由于 $A\mapsto i[H,A]$ 是导子，哈密顿量的项相互抵消）对每个 $X$ 给出
$$
\mathcal L^*(X^*X)-X^*\mathcal L^*(X)-\mathcal L^*(X^*)X=\sum_j[L_j,X]^*[L_j,X]:
$$
对单个跳跃，六个反对易子项中有四个两两抵消，另外两个
$\frac12X^*L^*LX+\frac12X^*L^*LX$ 相加，所以左边是 $L^*X^*XL-X^*L^*XL-L^*X^*LX+X^*L^*LX$，即
$(XL-LX)^*(XL-LX)$。设 $X\in\ker\mathcal L^*$。由于 $\mathcal L^*$ 是保厄米的，$\mathcal L^*(X^*)=\mathcal L^*(X)^*=0$，
所以 $\mathcal L^*(X^*X)=\sum_j[L_j,X]^*[L_j,X]=:Y\ge0$，并且
$\operatorname{Tr}(\sigma Y)=\operatorname{Tr}(\mathcal L(\sigma)X^*X)=0$。由于 $\sigma>0$，
$\operatorname{Tr}(\sigma Y)=\operatorname{Tr}(\sigma^{1/2}Y\sigma^{1/2})$，而 $\sigma^{1/2}Y\sigma^{1/2}\ge0$，它只有在
$\sigma^{1/2}Y\sigma^{1/2}=0$，即 $Y=0$ 时才为零；取迹得 $\sum_j\operatorname{Tr}([L_j,X]^*[L_j,X])=0$，各项非负，
所以对所有 $j$ 有 $[L_j,X]=0$。对 $X^*\in\ker\mathcal L^*$ 用同样的论证
得到 $[L_j,X^*]=0$，即 $[L_j^*,X]=0$。于是 $X$ 与 $L_j^*L_j$ 可交换，$\mathcal L^*(X)$ 的耗散部分
为零，从而 $0=\mathcal L^*(X)=i[H,X]$。由假设 $X\in\mathbb C1$；反之 $\mathcal L^*(1)=0$。所以
$\ker\mathcal L^*=\mathbb C1$。在关于配对 $(\rho,A)\mapsto\operatorname{Tr}(\rho A)$ 互为对偶的 $M_n$ 的基 $(E_{ij})$ 与
$(E_{ji})$ 下，$\mathcal L^*$ 的矩阵是 $\mathcal L$ 的矩阵的转置；所以二者的秩相同，
$\dim\ker\mathcal L=\dim\ker\mathcal L^*=1$。$\mathbb C\sigma$ 中迹为 $1$ 的态就是 $\sigma$。$\square$

这是 Frigerio 的一个判据 [*Commun. Math. Phys.* 63, 269 (1978)] 的有限维情形，
8.5 节用到了它。两个假设都不可缺少。在第 8 章的习题 4 中，态 $\tfrac12\,1$ 是忠实的
稳态，但 $H$ 和 $\sigma_z$ 与所有对角矩阵可交换，而稳态张成二维的对角矩阵空间。
反过来，在 $M_3$ 上令 $H=0$，$L_1=|e_1\rangle\langle e_3|$，$L_2=|e_2\rangle\langle e_3|$。由 $L_1,L_2$ 生成的
$*$-代数包含 $L_1L_2^*=E_{12}$ 与 $L_1^*L_1=E_{33}$，从而包含所有矩阵单位，所以与 $L_j,L_j^*$ 可交换的
只有 $1$ 的倍数；但满足 $\rho=P\rho P$（$P=E_{11}+E_{22}$）的每个 $\rho$ 都是稳态，因为
$L_j\rho=0=\rho L_j^*$ 且 $\{E_{33},\rho\}=0$，所以 $\dim\ker\mathcal L\ge4$；由定理 D.6，没有一个稳态
是忠实的。对上面的随机生成元，代码检验了假设与结论（$\mathcal L^*$ 在 $A^T$ 的按行展开向量上的矩阵是
$S^T$: $\operatorname{vec}(\mathcal L^*(A)^T)=S^T\operatorname{vec}(A^T)$），并检验了第二个例子：

```python
def kernel_dim(S, tol=1e-9):
    return int((np.linalg.svd(S, compute_uv=False) < tol).sum())

def commutant_dim(ops, n):                    # dim {X : XA = AX for all A in ops}
    return kernel_dim(np.concatenate([np.kron(A, np.eye(n)) - np.kron(np.eye(n), A.T) for A in ops]))

w, V = np.linalg.eig(S)
sigma = V[:, np.argmin(abs(w))].reshape(n, n)             # row-major: vec(X)_{(a,b)} = X_ab
sigma = sigma / np.trace(sigma)
assert np.abs(S @ sigma.reshape(-1)).max() < 1e-12        # stationary
assert np.abs(sigma - sigma.conj().T).max() < 1e-10 and np.linalg.eigvalsh(sigma).min() > 1e-3   # faithful
assert commutant_dim([H] + Ls + [L.conj().T for L in Ls], n) == 1
assert kernel_dim(S) == 1 and kernel_dim(S.T) == 1        # Theorem D.6
assert np.abs(S.T @ np.eye(n).reshape(-1)).max() < 1e-12  # ker L^* = C 1

e3 = np.eye(3)
L_abs = [np.outer(e3[0], e3[2]), np.outer(e3[1], e3[2])]       # L_1 = |e_1><e_3|, L_2 = |e_2><e_3|
S_abs = superop(lambda R: sum(L @ R @ L.T - 0.5 * (L.T @ L @ R + R @ L.T @ L) for L in L_abs), 3)
assert commutant_dim(L_abs + [L.T for L in L_abs], 3) == 1
print("second example: dim ker L =", kernel_dim(S_abs))
assert kernel_dim(S_abs) == 4                             # no faithful stationary state, no uniqueness
```

## D.4 有扇区的代数

设 $\mathcal M=\bigoplus_{c=1}^CM_{k_c}\subseteq M_D$ 为块对角矩阵，$D=\sum_ck_c$，其极小中心投影为
$P_c$，并设 $E(X)=\sum_cP_cXP_c$ 是 $M_D$ 到 $\mathcal M$ 上的条件期望。与命题 8.8 的证明一样，对
$B\in\mathcal M$ 有 $E(B)=B$，$E(XB)=E(X)B$，$E(BX)=BE(X)$，并且 $E(X^*)=E(X)^*$。对 CP 映射
$\Phi:\mathcal M\to\mathcal M$，$\Phi^*$ 是它的 blunt 对偶（命题 6.8）。

**定理 D.7.** 设 $\mathcal L$ 是 $\mathcal M$ 上的量子动力学半群 $(T_t)$ 的生成元。则存在 CP 映射
$\Phi:\mathcal M\to\mathcal M$ 和 $H=H^*\in\mathcal M$，使得
$$
\mathcal L(\rho)=\Phi(\rho)-\tfrac12\{\Phi^*(1),\rho\}-i[H,\rho]\qquad(\rho\in\mathcal M).
$$
由命题 6.4，$\Phi(\rho)_d=\sum_{c,l}K^{dc}_l\rho_cK^{dc*}_l$，其 Kraus 算子在扇区之间作用；
反之，每个这样的 $\mathcal L$ 都生成 $\mathcal M$ 上的量子动力学半群（命题 8.8(1)）。
所以 8.3 节的生成元（其中 $\Phi=\sum_jg_jL_j\cdot L_j^*+\sum_mg_m\Phi_m$）恰好就是 $\mathcal M$ 上
量子动力学半群的生成元。

*证明。* *第 1 步：$M_D$ 上的一个 CP 映射。* 由命题 6.4，
$T_t(\rho)_d=\sum_{c,l}K^{dc}_l(t)\rho_cK^{dc}_l(t)^*$。设 $N^{dc}_l(t)\in M_D$ 在位置 $(d,c)$ 处为块
$K^{dc}_l(t)$，其余处为零。对 $X\in M_D$，$N^{dc}_lXN^{dc*}_l$ 是位于位置 $(d,d)$ 的块
$K^{dc}_lX_{cc}K^{dc*}_l$，其中 $X_{cc}$ 是第 $c$ 个对角块，即 $E(X)$ 的第 $c$ 个块。因此
$\hat T_t=T_t\circ E=\sum_{d,c,l}N^{dc}_l\cdot N^{dc*}_l$ 是 $M_D$ 上的 CP 映射（定理 6.2），且
$\hat T_0=E$。由引理 8.5，$(\hat T_t-E)/t=\big((T_t-\mathrm{id})/t\big)\circ E\to\mathcal L\circ E$。

*第 2 步：在扇区向量之外的正性。* 令 $p_c=(1\otimes P_c)\Omega\in\mathbb C^D\otimes\mathbb C^D$，
$W=\operatorname{span}\{p_c\}$。由引理 D.1(c)，$J(E)=\sum_c|p_c\rangle\langle p_c|$，所以对 $v\perp W$ 和 $t>0$，
$0\le\frac1t\langle v,J(\hat T_t)v\rangle=\langle v,J((\hat T_t-E)/t)v\rangle\to\langle v,J(\mathcal L\circ E)v\rangle$。
每个 $T_t$ 都是保厄米的（由其 Kraus 形式），因此 $\mathcal L$ 也是（同定理 D.3），$E$ 与
$\mathcal L\circ E$ 也是。

*第 3 步：分解。* 由引理 D.1(b)，$\langle p_c,p_d\rangle=\operatorname{Tr}(P_cP_d)=\delta_{cd}k_c$，所以
$q_c=p_c/\sqrt{k_c}=(1\otimes Q_c)\Omega$（$Q_c=P_c/\sqrt{k_c}$）是规范正交的。对
$\Psi=\mathcal L\circ E$ 应用引理 D.2，得到 $M_D$ 上的 CP 映射 $\hat\Phi=\sum_j\hat L_j\cdot\hat L_j^*$ 和 $K_c\in M_D$，使得
$$
\mathcal L(E(A))=\hat\Phi(A)+\sum_c\big(K_cAP_c+P_cAK_c^*\big)/\sqrt{k_c}\qquad(A\in M_D).
$$
对 $A\in\mathcal M$，$E(A)=A$ 且 $AP_c=P_cA$，所以这个和等于 $\tilde GA+A\tilde G^*$，其中
$\tilde G=\sum_cK_cP_c/\sqrt{k_c}$。

*第 4 步：压缩。* 对 $A\in\mathcal M$，$\mathcal L(A)\in\mathcal M$，所以
$\mathcal L(A)=E(\mathcal L(A))=E(\hat\Phi(A))+GA+AG^*$，其中 $G=E(\tilde G)\in\mathcal M$。$\mathcal M$ 上的映射
$\Phi=E\circ\hat\Phi$ 为
$\Phi(A)_d=\sum_{c,j}(P_d\hat L_jP_c)A_c(P_d\hat L_jP_c)^*$：由命题 6.4 它是 CP 的，其 Kraus 算子是
各 $\hat L_j$ 的块。

*第 5 步：保迹性。* 与定理 D.4 中一样，对 $A\in\mathcal M$ 有 $\operatorname{Tr}\mathcal L(A)=0$，并且
$\operatorname{Tr}\Phi(A)=\operatorname{Tr}(A\Phi^*(1))$，所以对所有 $A\in\mathcal M$ 有 $\operatorname{Tr}(AX)=0$，其中
$X=\Phi^*(1)+G+G^*\in\mathcal M$ 是自伴的。取 $A=X$ 得 $X=0$（blunt 迹是忠实的，定理 3.1）。
令 $H=\frac i2(G-G^*)\in\mathcal M$，则 $G=-iH-\frac12\Phi^*(1)$，且 $G\rho+\rho G^*=-i[H,\rho]-\frac12\{\Phi^*(1),\rho\}$。$\square$

向量 $p_c$ 取代了 $\Omega=\sum_cp_c$，因为 $\hat T_0=E$ 而不是 $\mathrm{id}$：$J(\mathcal L\circ E)$ 在 $\Omega^\perp$
上不一定是正的，在下面的例子中它就不是正的。

**推论 D.8（Heisenberg 绘景）。** 设 $(S_t)_{t\ge0}$ 是 $\mathcal M$ 上一族连续的保单位 CP 映射，满足
$S_0=\mathrm{id}$ 与 $S_{s+t}=S_sS_t$。则它的生成元为 $A\mapsto\Psi(A)+K^*A+AK$，其中
$\Psi$ 是 $\mathcal M$ 上的 CP 映射，$K\in\mathcal M$，且 $\Psi(1)+K+K^*=0$。

*证明。* 令 $T_t=S_t^*$（blunt 对偶）。它们是 CP 的（命题 6.8）和 TP 的（推论 6.9(b)），
$T_0=\mathrm{id}$，并且 $T_{s+t}=(S_tS_s)^*=S_s^*S_t^*=T_sT_t$，因为
$\operatorname{Tr}(S_tS_s(x)y)=\operatorname{Tr}(xS_s^*S_t^*(y))$ 且 $S_tS_s=S_{s+t}$。由于 $\Psi\mapsto\Psi^*$ 是
有限维空间上的线性映射，$t\mapsto T_t$ 是连续的，且 $(T_t-\mathrm{id})/t\to\mathcal L_S^*$，其中 $\mathcal L_S$ 是
$(S_t)$ 的生成元。定理 D.7 给出 $\mathcal L_S^*(\rho)=\Phi(\rho)+G\rho+\rho G^*$，其中 $G+G^*=-\Phi^*(1)$；而由
循环性，$\operatorname{Tr}((G\rho+\rho G^*)A)=\operatorname{Tr}(\rho(G^*A+AG))$，所以由推论 6.9(c)，
$\mathcal L_S(A)=\Phi^*(A)+G^*A+AG$。取 $\Psi=\Phi^*$（由命题 6.8 它是 CP 的）和 $K=G$。$\square$

这是定理 8.7 之后所引用的 Christensen 与 Evans 定理的有限维情形。我们在 $M_1\oplus M_2\subseteq M_3$
上，对 `lindblad_rhs` 的一个生成元检验定理 D.7，它有一个扇区内部的跳跃和一个随机的扇区之间的跳跃，
检验按照证明的五个步骤进行：

```python
alg = TypeIAlgebra([1, 2], [1, 2], precision='double')    # M = M_1 (+) M_2 in M_3
D, ks = 3, alg.k_factors
P = [np.diag([1.0, 0, 0]), np.diag([0, 1.0, 1.0])]
E = lambda X: sum(Pc @ X @ Pc for Pc in P)               # the conditional expectation onto M
def blocks(X):                                            # dense block-diagonal 3 x 3 -> library layout
    out = torch.zeros(1, 2, 2, 2, dtype=c128)
    out[0, 0, :1, :1] = torch.from_numpy(X[:1, :1])
    out[0, 1] = torch.from_numpy(X[1:, 1:])
    return out
dense = lambda Y: sla.block_diag(Y[0, 0, :1, :1].numpy(), Y[0, 1].numpy())
Hs = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
inside = alg.operator(torch.randn(1, 2, 2, 2, dtype=c128))       # a jump inside the sectors
between = channels.random_inter_sector_channel(alg, alg, kraus_rank=2)
rhs = dynamics.lindblad_rhs(Hs, [inside, between], rates=[0.7, 1.3], algebra=alg)
L_M = lambda X: dense(rhs(0.0, blocks(X)))               # the generator on M (X block diagonal)
J = choi_of(superop(lambda X: L_M(E(X)), D), D)          # J(L o E) on M_3
p = [vec_of(Pc) for Pc in P]
Q = np.stack([pc / np.linalg.norm(pc) for pc in p], 1)    # orthonormal q_c = p_c / sqrt(k_c)
Pw = np.eye(D * D) - Q @ Q.T                              # projection onto the complement of W
assert np.linalg.eigvalsh(Pw @ J @ Pw).min() > -1e-12     # Step 2
w, U = np.linalg.eigh(Pw @ J @ Pw)                        # Step 3 (Lemma D.2)
jumps = [math.sqrt(x) * op_of(U[:, k], D) for k, x in enumerate(w) if x > 1e-10]
xis = [J @ Q[:, c] - 0.5 * sum((Q[:, d] @ J @ Q[:, c]) * Q[:, d] for d in range(2)) for c in range(2)]
G = E(sum(op_of(xi, D) @ Pc / math.sqrt(k) for xi, Pc, k in zip(xis, P, ks)))   # Step 4
Phi = lambda X: E(sum(K @ X @ K.conj().T for K in jumps))
Phi_dual_one = E(sum(K.conj().T @ K for K in jumps))      # Phi^*(1)
H_M = 0.5j * (G - G.conj().T)
assert np.abs(G + G.conj().T + Phi_dual_one).max() < 1e-12      # Step 5
for _ in range(3):
    X = E(cplx(D, D))
    gksl = Phi(X) - 0.5 * (Phi_dual_one @ X + X @ Phi_dual_one) - 1j * (H_M @ X - X @ H_M)
    assert np.abs(L_M(X) - gksl).max() < 1e-12            # Theorem D.7
    dual = lambda A: E(sum(K.conj().T @ A @ K for K in jumps)) + G.conj().T @ A + A @ G   # Corollary D.8
    assert abs(np.trace(L_M(X) @ (Y := E(cplx(D, D)))) - np.trace(X @ dual(Y))) < 1e-12
Pomega = np.eye(D * D) - np.outer(sum(p), sum(p)) / D     # projection onto the complement of Omega
print(f"numerical observation: min eigenvalue of J(L o E) on the complement of Omega: "
      f"{np.linalg.eigvalsh(Pomega @ J @ Pomega).min():.2f}")
assert np.linalg.eigvalsh(Pomega @ J @ Pomega).min() < -0.1
```

## D.5 算子代数纠错

我们沿用 10.3 节的设定。$\Phi(\rho)=\sum_{i=1}^rE_i\rho E_i^*$ 是 $M_D$ 上的信道，$\mathcal C=P\mathbb C^D$
是一个码，$B(\mathcal C)=\{X\in M_D:X=PXP\}$ 是 $\mathcal C$ 上的算子，$\mathfrak A\subseteq B(\mathcal C)$ 是以
$P$ 为单位的 $*$-代数（由定理 2.2，它是 $\mathcal C$ 上的冯·诺依曼代数），$\mathfrak A'$ 是它在
$B(\mathcal C)$ 中的交换子。若存在 $M_D$ 上的信道 $\mathcal R$，使得对所有 $A\in\mathfrak A$ 有
$P(\mathcal R\circ\Phi)^*(A)P=A$，则称 $\mathfrak A$ 对 $\Phi$ 是**可纠正的**。若 $\mathcal R$ 的 Kraus 算子为 $R_m$，则
$\mathcal R\circ\Phi$ 的 Kraus 算子为 $R_mE_i$，并且 $\sum_mR_m^*R_m=1$。

**定理 D.9（Bény–Kempf–Kribs；定理 10.7）。** $\mathfrak A$ 对 $\Phi$ 是可纠正的，当且仅当对所有 $i,j$ 有
$PE_i^*E_jP\in\mathfrak A'$。在这种情况下，证明中构造的信道 $\mathcal R$ 纠正 $\mathfrak A$。

*证明。* *必要性。* 设 $\mathcal R$ 纠正 $\mathfrak A$，令 $W_{mi}=R_mE_iP$，于是对所有 $X\in M_D$ 有
$\sum_{m,i}W_{mi}^*XW_{mi}=P(\mathcal R\circ\Phi)^*(X)P$。对 $X\in\mathfrak A$ 它等于 $X$；对 $X=1$ 它等于 $P$，
因为 $(\mathcal R\circ\Phi)^*$ 是保单位的（推论 6.9(a)）。设 $p\in\mathfrak A$ 是投影；则 $p\le P$，
并且 $\sum W_{mi}^*(1-p)W_{mi}=P-p$。对 $x=px$，
$\sum_{m,i}\|(1-p)W_{mi}x\|^2=\langle x,(P-p)x\rangle=0$；对 $x=(P-p)x$，
$\sum_{m,i}\|pW_{mi}x\|^2=\langle x,px\rangle=0$。因此 $(1-p)W_{mi}p=0$ 且 $pW_{mi}(P-p)=0$，又由于
$W_{mi}=W_{mi}P$，
$$
W_{mi}p=pW_{mi}p=pW_{mi}P=pW_{mi}.
$$
$\mathfrak A$ 中每个自伴的 $A$ 都是 $\mathfrak A$ 中投影的线性组合：在 $\mathcal C$ 上应用定理 1.3，
$A=\sum_k\lambda_kp_k$，其中 $\lambda_k$ 互不相同，$p_k$ 是正交投影，$\sum_kp_k=P$，并且
$p_k=\prod_{l\neq k}(A-\lambda_lP)/(\lambda_k-\lambda_l)\in\mathfrak A$（把两边作用到 $\mathcal C$ 中的特征向量上以及作用到 $\mathcal C^\perp$ 上；空乘积为 $P$）。$\mathfrak A$ 的
每个元素都是 $A_1+iA_2$，其中 $A_1,A_2\in\mathfrak A$ 自伴。所以对所有 $A\in\mathfrak A$ 有 $W_{mi}A=AW_{mi}$，
取伴随得 $W_{mi}^*A=AW_{mi}^*$。由于 $\sum_mW_{mi}^*W_{mj}=PE_i^*(\sum_mR_m^*R_m)E_jP=PE_i^*E_jP$，
$$
A\,PE_i^*E_jP=\sum_mW_{mi}^*AW_{mj}=PE_i^*E_jP\,A\qquad(A\in\mathfrak A),
$$
即 $PE_i^*E_jP\in\mathfrak A'$。

*充分性。* 对 $\mathcal C$ 上的 $\mathfrak A$ 应用定理 2.6，存在 $\mathcal C$ 的一组规范正交基 $|c;a,b\rangle$
（$c=1,\dots,C$，$a=1,\dots,k_c$，$b=1,\dots,m_c$）——即第 $c$ 个直和项中的 $e_a\otimes f_b$ 在酉映射
$V$ 下的像——在这组基下 $\mathfrak A$ 由算子 $\sum_{c,a,a',b}(A_c)_{aa'}|c;a,b\rangle\langle c;a',b|$
（$A_c\in M_{k_c}$）组成，而 $\mathfrak A'$ 由算子 $\sum_{c,a,b,b'}(Y_c)_{bb'}|c;a,b\rangle\langle c;a,b'|$（$Y_c\in M_{m_c}$）组成。
因此假设给出矩阵 $\Lambda^c_{ij}\in M_{m_c}$，使得
$$
\langle c;a,b|E_i^*E_j|d;a',b'\rangle=\delta_{cd}\,\delta_{aa'}\,(\Lambda^c_{ij})_{bb'} .
$$
*第 1 步。* 对每个 $c$ 和每一对 $\alpha=(i,b)$，令 $F^c_\alpha:\mathbb C^{k_c}\to\mathbb C^D$，
$F^c_\alpha e_a=E_i|c;a,b\rangle$，并对 $\beta=(j,b')$ 令 $\Lambda^c_{\alpha\beta}=(\Lambda^c_{ij})_{bb'}$。上面的显示式表明
$F^{c*}_\alpha F^d_\beta=\delta_{cd}\Lambda^c_{\alpha\beta}1_{k_c}$：这是关于码 $\mathbb C^{k_c}$ 上的"错误" $F^c_\alpha$ 的
Knill–Laflamme 条件，此外还有扇区之间的正交性。

*第 2 步。* 与定理 10.5 的证明一样，$\Lambda^c$ 是自伴的，因为
$\Lambda^c_{\beta\alpha}1=F^{c*}_\beta F^c_\alpha=(F^{c*}_\alpha F^c_\beta)^*=\bar\Lambda^c_{\alpha\beta}1$；它也是正的，因为
$\sum\bar z_\alpha z_\beta\Lambda^c_{\alpha\beta}1=(\sum z_\alpha F^c_\alpha)^*(\sum z_\beta F^c_\beta)\ge0$。写
$\Lambda^c=U^c\operatorname{diag}(d^c_\kappa)U^{c*}$，其中 $U^c$ 是酉矩阵，$d^c_\kappa\ge0$，并令
$G^c_\kappa=\sum_\alpha u^c_{\alpha\kappa}F^c_\alpha$。则 $G^{c*}_\kappa G^d_\lambda=\delta_{cd}\delta_{\kappa\lambda}d^c_\kappa1$，所以
若 $d^c_\kappa=0$ 则 $G^c_\kappa=0$，否则 $V^c_\kappa=G^c_\kappa/\sqrt{d^c_\kappa}$ 是等距，满足
$V^{c*}_\kappa V^d_\lambda=\delta_{cd}\delta_{\kappa\lambda}1$。反解得 $F^c_\alpha=\sum_\kappa\bar u^c_{\alpha\kappa}G^c_\kappa
=\sum_{\kappa:d^c_\kappa>0}\bar u^c_{\alpha\kappa}\sqrt{d^c_\kappa}\,V^c_\kappa$，因此
$$
V^{c'*}_\kappa F^d_\beta=\delta_{c'd}\,\bar u^d_{\beta\kappa}\sqrt{d^d_\kappa}\,1,\qquad
F^{c*}_\alpha V^{c'}_\kappa=\delta_{cc'}\,u^c_{\alpha\kappa}\sqrt{d^c_\kappa}\,1 .
$$
*第 3 步：恢复。* 各 $V^c_\kappa V^{c*}_\kappa$（$d^c_\kappa>0$）是两两正交的投影；令 $Q$ 为它们的
和。令 $\iota_c:\mathbb C^{k_c}\to\mathbb C^D$，$\iota_ce_a=|c;a,1\rangle$，这是一个等距，并令
$$
\mathcal R(\sigma)=\sum_{c}\sum_{\kappa:d^c_\kappa>0}\iota_cV^{c*}_\kappa\,\sigma\,V^c_\kappa\iota_c^*+(1-Q)\sigma(1-Q).
$$
它是一个信道：$\sum_{c,\kappa}V^c_\kappa\iota_c^*\iota_cV^{c*}_\kappa+(1-Q)^2=Q+1-Q=1$。它测量扇区 $c$ 和
综合征 $\kappa$，撤销 $\mathbb C^{k_c}$ 上的错误，并把结果写入副本 $b=1$。

*第 4 步：正确性。* $\mathcal R^*(X)=\sum_{c,\kappa}V^c_\kappa\iota_c^*X\iota_cV^{c*}_\kappa+(1-Q)X(1-Q)$，并且对
块为 $A_c$ 的 $A\in\mathfrak A$ 有 $\iota_c^*A\iota_c=A_c$。由于每个 $F^d_\beta$ 都是各 $V^d_\kappa$ 的线性组合，
$(1-Q)F^d_\beta=0$。结合第 2 步，对 $\alpha=(i,b)$，$\beta=(j,b')$，
$$
\langle c;a,b|E_i^*\mathcal R^*(A)E_j|d;a',b'\rangle=\sum_{c',\kappa}\langle e_a,F^{c*}_\alpha V^{c'}_\kappa A_{c'}V^{c'*}_\kappa F^d_\beta e_{a'}\rangle
=\delta_{cd}\sum_\kappa u^c_{\alpha\kappa}d^c_\kappa\bar u^c_{\beta\kappa}(A_c)_{aa'}=\delta_{cd}(\Lambda^c_{ij})_{bb'}(A_c)_{aa'} .
$$
取 $j=i$ 并对 $i$ 求和：由上面带 $\sum_iE_i^*E_i=1$ 的显示式，$\sum_i\Lambda^c_{ii}=1_{m_c}$，所以
$$
\langle c;a,b|(\mathcal R\circ\Phi)^*(A)|d;a',b'\rangle=\delta_{cd}\,\delta_{bb'}(A_c)_{aa'}=\langle c;a,b|A|d;a',b'\rangle .
$$
$P(\mathcal R\circ\Phi)^*(A)P$ 与 $A$ 都属于 $B(\mathcal C)$，并且在基 $|c;a,b\rangle$ 下有相同的矩阵，所以二者
相等。$\square$

**推论 D.10.** (a) 对 $\mathfrak A=B(\mathcal C)$，$\mathfrak A'=\mathbb CP$，定理 D.9 就是 Knill–Laflamme
定理 10.5（$C=1$，$m_1=1$；第 3 步就是定理 10.5 的恢复）。(b) 对 $\mathfrak A=\operatorname{span}\{P_1,\dots,P_C\}$，其中
$P_c\ne0$ 是两两正交的投影，$\sum_cP_c=P$，它就是命题 10.6 中 2 与 3 的等价性。(c)（*无噪声
子系统。*）沿用命题 10.1(a) 的记号，设 $\mathcal C=V(\mathbb C^{k_c}\otimes\mathbb C^{m_c})$ 是 $\Phi$ 的相互作用代数
$\mathcal A$ 的第 $c$ 个直和项，而 $\mathfrak A=V(1_{k_c}\otimes M_{m_c})V^*$ 作用在 $\mathcal C$ 上，交换两个因子之后，它是
形如 $M_{m}\otimes1_{k}$ 的代数。则 $\mathfrak A$ 被 $\mathcal R=\mathrm{id}$ 纠正，并且
$PE_i^*E_jP=V(E_{i,c}^*E_{j,c}\otimes1_{m_c})V^*\in\mathfrak A'$ ($P=VV^*$)。

*证明。* (a) $\mathfrak A'=\mathbb CP$ 是上面对 $\mathfrak A'$ 的描述中 $C=1$，$m_1=1$ 的情形，而
$PE_i^*E_jP\in\mathbb CP$ 就是 Knill–Laflamme 条件。$B(\mathcal C)$ 的可纠正性就是 10.2 节意义下的纠错：
若对所有支撑在 $\mathcal C$ 中的态 $\rho$ 有 $\mathcal R(\Phi(\rho))=\rho$，则对所有这样的 $\rho$ 有
$\operatorname{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(A)P)=\operatorname{Tr}(\rho A)$，而这些 $\rho$ 张成 $B(\mathcal C)$；
逆命题在 10.3 节中证明。当 $C=1$，$m_1=1$ 时，第 3 步的 Kraus 算子 $\iota_1V^{1*}_\kappa$ 就是
定理 10.5 中的那些，因为 $\iota_1\iota_1^*=P$。(b) $T\in B(\mathcal C)$ 与所有 $P_c$ 可交换当且仅当对 $c\ne d$ 有
$P_cTP_d=0$，这就是命题 10.6 的条件 2。(c) 每个 $A\in\mathfrak A$ 在其他直和项上延拓为 $0$ 后都属于
$\mathcal A'=V(\bigoplus_c1_{k_c}\otimes M_{m_c})V^*$（定理 2.6），所以由命题 10.1(a) 有 $\Phi^*(A)=A$，并且
$PAP=A$。$PE_i^*E_jP$ 的公式来自 $V^*E_iV=\bigoplus_cE_{i,c}\otimes1_{m_c}$，而
$E_{i,c}^*E_{j,c}\otimes1_{m_c}$ 与 $1_{k_c}\otimes M_{m_c}$ 可交换。$\square$

由定理 2.6，每个 $\mathfrak A$ 都 $*$-同构于某个 $\bigoplus_cM_{k_c}$。当 $C\ge2$ 且某个 $k_c\ge2$ 时，
可纠正的 $\mathfrak A$ 既保护扇区标记 $c$（经典信息），又保护每个 $M_{k_c}$ 的量子信息；
这样的码称为**混合**码。作为例子，我们在 $\mathbb C^8$ 中一个随机的五维码上取
$\mathfrak A=(M_2\otimes1_2)\oplus M_1$，噪声通过 $\sum_iB_{c,i}^*B_{c,i}=1_{m_c}$ 作用在规范因子 $\mathbb C^{m_c}$ 上，
把各扇区映到相互正交的子空间，打乱 $\mathcal C^\perp$，然后再被一个随机酉矩阵混合（定理 6.7）。
Knill–Laflamme 条件不成立，定理 D.9 的条件成立，而证明中的恢复还原了每个 $A\in\mathfrak A$，
却不能还原 $\mathcal C$ 上的一般算子：

```python
D, dims = 8, [(2, 2), (1, 1)]                             # A = (M_2 (x) 1_2) (+) (M_1 (x) 1_1), dim C = 5
U8 = np.linalg.qr(cplx(D, D))[0]
Vc = np.linalg.qr(cplx(D, 5))[0]                          # orthonormal basis |c; a, b> of the code
W = [Vc[:, :4], Vc[:, 4:]]                                # sector c: columns indexed by (a, b), row-major
Pcode = sum(Wc @ Wc.conj().T for Wc in W)
Y = [U8[:, :4], U8[:, 4:6]]                               # orthogonal images C^{k_c} (x) C^{s_c}, s_c = 2
r = 3
B = [np.linalg.qr(cplx(2 * r, m))[0].reshape(r, 2, m) for _, m in dims]   # sum_i B_i^* B_i = 1_{m_c}
E_code = [sum(Yc @ np.kron(np.eye(k), Bc[i]) @ Wc.conj().T for Yc, Bc, Wc, (k, m) in zip(Y, B, W, dims))
          for i in range(r)]
kraus = E_code + [np.linalg.qr(cplx(D, D))[0] @ (np.eye(D) - Pcode)]
mix = np.linalg.qr(cplx(r + 1, r + 1))[0]                 # Theorem 6.7: an equivalent family
kraus = [sum(mix[i, l] * K for l, K in enumerate(kraus)) for i in range(r + 1)]
assert np.abs(sum(K.conj().T @ K for K in kraus) - np.eye(D)).max() < 1e-12      # a channel

units = []                                                 # matrix units of A
for Wc, (k, m) in zip(W, dims):
    for a1 in range(k):
        for a2 in range(k):
            units.append(Wc @ np.kron(np.eye(k)[:, [a1]] @ np.eye(k)[[a2], :], np.eye(m)) @ Wc.conj().T)
comm = lambda X: max(np.abs(X @ A - A @ X).max() for A in units)
assert max(comm(Pcode @ Ki.conj().T @ Kj @ Pcode) for Ki in kraus for Kj in kraus) < 1e-12   # Theorem D.9
KL_dev = max(np.abs(Pcode @ Ki.conj().T @ Kj @ Pcode - np.trace(Pcode @ Ki.conj().T @ Kj) / 5 * Pcode).max()
             for Ki in kraus for Kj in kraus)
print(f"Knill-Laflamme condition violated by {KL_dev:.2f}")
assert KL_dev > 0.05

R_kraus, Qsum = [], np.zeros((D, D), complex)             # Steps 1-3 of the proof
for Wc, (k, m) in zip(W, dims):
    F = [Ki @ Wc[:, b::m] for Ki in kraus for b in range(m)]   # F_(i,b) e_a = E_i |c; a, b>
    Lam = np.array([[(Fa.conj().T @ Fb)[0, 0] for Fb in F] for Fa in F])
    assert max(np.abs(Fa.conj().T @ Fb - Lam[x, y] * np.eye(k)).max()
               for x, Fa in enumerate(F) for y, Fb in enumerate(F)) < 1e-12      # Step 1
    d, u = np.linalg.eigh(Lam)
    iota = Wc[:, 0::m]                                    # e_a -> |c; a, 1>
    for kk in range(len(d)):
        if d[kk] > 1e-12:
            Vk = sum(u[x, kk] * Fx for x, Fx in enumerate(F)) / math.sqrt(d[kk])
            R_kraus.append(iota @ Vk.conj().T)
            Qsum += Vk @ Vk.conj().T
R_kraus.append(np.eye(D) - Qsum)
assert np.abs(sum(R.conj().T @ R for R in R_kraus) - np.eye(D)).max() < 1e-10    # R is a channel
dual = lambda X: sum(Ki.conj().T @ R.conj().T @ X @ R @ Ki for R in R_kraus for Ki in kraus)
for _ in range(3):
    A = sum(complex(z) * Ue for z, Ue in zip(cplx(len(units)), units))
    assert np.abs(Pcode @ dual(A) @ Pcode - A).max() < 1e-10                       # A is corrected
B_gen = Pcode @ cplx(D, D) @ Pcode
print(f"a generic operator on the code: error {np.abs(Pcode @ dual(B_gen) @ Pcode - B_gen).max():.2f}")
assert np.abs(Pcode @ dual(B_gen) @ Pcode - B_gen).max() > 0.1
assert max(comm(R @ Ki @ Pcode) for R in R_kraus for Ki in kraus) < 1e-10   # proof of necessity: W_mi in A'
```

本附录的结果归功于 V. Gorini、A. Kossakowski 和 E. C. G. Sudarshan
[*J. Math. Phys.* 17, 821 (1976)] 以及 G. Lindblad [*Commun. Math. Phys.* 48, 119 (1976)]（定理 D.3–D.5），
E. Christensen 和 D. E. Evans [*J. London Math. Soc.* 20, 358 (1979)]（一般冯·诺依曼代数中的推论 D.8），
A. Frigerio（定理 D.6），以及 C. Bény、A. Kempf 和 D. W. Kribs [*Phys. Rev. Lett.* 98, 100502
(2007)；*Phys. Rev. A* 76, 042303 (2007)]（定理 D.9）；上面的证明对于本课程所用的有限维陈述是完整的。
