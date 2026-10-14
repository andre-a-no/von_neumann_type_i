[English](../D_gksl.md) · [Español](../es/D_gksl.md) · [Français](../fr/D_gksl.md) · [Deutsch](../de/D_gksl.md) · [中文](../zh/D_gksl.md) · **日本語** · [Русский](../ru/D_gksl.md)

# 付録 D. 量子動力学的半群の生成子

第8章では、ゴリーニ–コサコフスキー–スダルシャン–リンドブラッド（GKSL）形の生成子が量子動力学的半群を
生み出すこと（定理 8.6）を証明し、その逆（定理 8.7）と、$\mathcal M=\bigoplus_cM_{k_c}$ に対するその版を用いた。
本付録では、GKSL 形における自由度と、定常状態の一意性の判定条件とともに、有限次元においてこれら両方を
証明する。道具は第6章のチョイ行列である：生成子は一般には完全正値ではないが、そのチョイ行列は
最大エンタングルベクトルの直交補空間上で正である。同じ種類の議論——作用素の族の正値性、続いて分解——により、
最後の節では第10章で述べた作用素環による誤り訂正の定理（定理 10.7）を、明示的な回復チャネルとともに証明する。

**本付録で学ぶこと**：完全正値な半群の生成子が条件付き完全正値である理由、生成子のチョイ行列から $H$ と
ジャンプ作用素を読み取る方法、それらがどのような意味で一意であるか、セクターをもつ環上のすべての生成子が
8.3節で用いた形をもつ理由、定常状態が一意となるのはいつか、そして回復チャネルがどのオブザーバブルの環を
保護できるか。

全体を通じて、$M_n=M_n(\mathbb C)$、$E_{ij}=|i\rangle\langle j|$ はその行列単位、**超演算子** とは
$M_n$（または $\mathcal M$）上の線形写像のことであり、
$$
\Omega=\sum_{i=1}^ne_i\otimes e_i\in\mathbb C^n\otimes\mathbb C^n,\qquad
J(\Psi)=\sum_{i,j=1}^nE_{ij}\otimes\Psi(E_{ij})=[\Psi(E_{ij})]_{i,j=1}^n
$$
はそれぞれ最大エンタングルベクトルと超演算子 $\Psi$ の **チョイ行列** である（定理 6.2。ここでは
$\Psi$ は正値である必要はない）。ベクトル $v,w$ に対し、$|v\rangle\langle w|$ は作用素
$x\mapsto\langle w,x\rangle v$ を表す。すべての $A$ に対し $\Psi(A^*)=\Psi(A)^*$ であるとき、超演算子 $\Psi$ は
**エルミート性を保つ** という。量子動力学的半群とその生成子は 8.2節（補題 8.5）のものである：
$T_t=e^{t\mathcal L}$ であり、$T_t$ はチャネルである。

## D.1 構成要素のチョイ行列

**補題 D.1.** $\Psi$ を $M_n$ 上の超演算子とする。

(a) $\Psi\mapsto J(\Psi)$ は、$M_n$ 上の超演算子全体から $M_n(M_n)=M_{n^2}$ の上への線形全単射である。
(b) すべての $v\in\mathbb C^n\otimes\mathbb C^n$ は、ちょうど一つの $K\in M_n$ により $v=(1\otimes K)\Omega$ と書け、
$\langle(1\otimes K)\Omega,(1\otimes L)\Omega\rangle=\operatorname{Tr}(K^*L)$ である。特に
$\langle\Omega,(1\otimes K)\Omega\rangle=\operatorname{Tr}K$ である。
(c) $K,L\in M_n$ に対し、超演算子 $A\mapsto KAL^*$ のチョイ行列は
$|(1\otimes K)\Omega\rangle\langle(1\otimes L)\Omega|$ である。特に $J(\mathrm{id})=|\Omega\rangle\langle\Omega|$ である。
(d) $\Psi$ がエルミート性を保つための必要十分条件は、$J(\Psi)$ が自己共役であることである。任意の CP 写像はエルミート性を保つ。

*証明.* (a) $J(\Psi)$ の $(i,j)$ ブロックは $\Psi(E_{ij})$ であり、$E_{ij}$ たちは $M_n$ の基底をなすので、
$J$ は線形かつ単射であり、任意のブロック行列 $[X_{ij}]$ は $E_{ij}\mapsto X_{ij}$ で定まる超演算子のチョイ行列である。
(b) $(1\otimes K)\Omega=\sum_ie_i\otimes Ke_i$ の第 $i$ ブロックは $Ke_i$、すなわち $K$ の第 $i$ 列であり、
これらが $K$ を決定し、すべての $v$ はこのようにして得られる。すると
$\langle\sum_ie_i\otimes Ke_i,\sum_je_j\otimes Le_j\rangle=\sum_i\langle Ke_i,Le_i\rangle=\operatorname{Tr}(K^*L)$ である。
最後の主張には $(K,L)$ を $(1,K)$ で置き換えればよい。(c) 右辺の $(i,j)$ ブロックは
$Ke_i(Le_j)^*=KE_{ij}L^*$ である。(d) $J(\Psi)^*$ の $(i,j)$ ブロックは $\Psi(E_{ji})^*$ である。$\Psi$ が
エルミート性を保つならば $\Psi(E_{ji})^*=\Psi(E_{ji}^*)=\Psi(E_{ij})$ である。逆に、すべての $i,j$ に対し
$\Psi(E_{ji})^*=\Psi(E_{ij})$ であり $A=\sum a_{ij}E_{ij}$ ならば、
$\Psi(A)^*=\sum\bar a_{ij}\Psi(E_{ij})^*=\sum\bar a_{ij}\Psi(E_{ji})=\Psi(A^*)$ である。CP 写像はクラウス表示をもち
（定理 6.2）、$(\sum_lK_lAK_l^*)^*=\sum_lK_lA^*K_l^*$ である。$\square$

**補題 D.2（分解）.** $q_1,\dots,q_C$ を $\mathbb C^n\otimes\mathbb C^n$ の正規直交ベクトルとし、
$q_c=(1\otimes Q_c)\Omega$ と書き、$W=\operatorname{span}\{q_c\}$ とする。また $\Psi$ を、すべての $v\perp W$ に対し
$\langle v,J(\Psi)v\rangle\ge0$ をみたす、エルミート性を保つ超演算子とする。
$J=J(\Psi)$、$\xi_c=Jq_c-\frac12\sum_d\langle q_d,Jq_c\rangle q_d$ とおき、$K_c\in M_n$ を
$(1\otimes K_c)\Omega=\xi_c$ により定める。このとき、$P_W^\perp$ を $W^\perp$ の上への射影として
$r=\operatorname{rank}(P_W^\perp JP_W^\perp)$ 個の $L_1,\dots,L_r\in M_n$ が存在して
$$
\Psi(A)=\sum_{j=1}^rL_jAL_j^*+\sum_{c=1}^C\big(K_cAQ_c^*+Q_cAK_c^*\big)\qquad(A\in M_n),
$$
かつすべての $c,j$ に対し $\operatorname{Tr}(Q_c^*L_j)=0$、$j\neq k$ に対し $\operatorname{Tr}(L_j^*L_k)=0$ となる。

*証明.* 補題 D.1(d) により $J$ は自己共役である。$P=\sum_c|q_c\rangle\langle q_c|$ を $W$ の上への射影とし、
$P^\perp=P_W^\perp=1-P$ とする。$J=(P+P^\perp)J(P+P^\perp)$ を展開し、$PJP^\perp=PJ-PJP$ と
$P^\perp JP=JP-PJP$ を用いると
$$
J=P^\perp JP^\perp+\big(PJ+JP-PJP\big).
$$
第一項は自己共役であり、$P^\perp v\perp W$ であるから
$\langle v,P^\perp JP^\perp v\rangle=\langle P^\perp v,JP^\perp v\rangle\ge0$ である。よって $P^\perp JP^\perp\ge0$ である。
定理 1.3 により $P^\perp JP^\perp=\sum_{j=1}^r|l_j\rangle\langle l_j|$ である。ここで $l_j$ は零でない固有値 $\mu_j$ の
固有ベクトルを $\sqrt{\mu_j}$ 倍したものである：それらは互いに直交し、$P^\perp JP^\perp$ の値域に属し、
この値域は $W$ と直交する。$l_j=(1\otimes L_j)\Omega$ と書く。補題 D.1(b) により
$\operatorname{Tr}(Q_c^*L_j)=\langle q_c,l_j\rangle=0$ であり、$j\ne k$ に対し $\operatorname{Tr}(L_j^*L_k)=\langle l_j,l_k\rangle=0$ である。
また補題 D.1(c) により、$P^\perp JP^\perp$ は $A\mapsto\sum_jL_jAL_j^*$ のチョイ行列である。第二項については、
$JP=\sum_c|Jq_c\rangle\langle q_c|$、$PJP=\sum_{c,d}\langle q_d,Jq_c\rangle|q_d\rangle\langle q_c|$ であるから
$$
\sum_c|\xi_c\rangle\langle q_c|=JP-\tfrac12PJP,\qquad\sum_c|q_c\rangle\langle\xi_c|=PJ-\tfrac12PJP
$$
（第二式は第一式の随伴である。$PJP$ は自己共役である）であり、それらの和は $PJ+JP-PJP$ である。
補題 D.1(c) により、$|\xi_c\rangle\langle q_c|$ と $|q_c\rangle\langle\xi_c|$ はそれぞれ
$A\mapsto K_cAQ_c^*$ と $A\mapsto Q_cAK_c^*$ のチョイ行列である。したがって主張する等式の両辺はチョイ行列 $J$ をもち、
補題 D.1(a) により等式が成り立つ。$\square$

## D.2 条件付き完全正値性

$M_n$ 上の超演算子 $\Psi$ は、すべての $v\perp\Omega$ に対し $\langle v,J(\Psi)v\rangle\ge0$ であるとき
**条件付き完全正値**（CCP）であるという。補題 D.1(b) により、これらの $v$ は $\operatorname{Tr}K=0$ をみたす
ベクトル $(1\otimes K)\Omega$ にほかならない。

**定理 D.3.** $M_n$ 上の超演算子 $\mathcal L$ について、次は同値である：

1. すべての $t\ge0$ に対し $e^{t\mathcal L}$ は CP である；
2. $\mathcal L$ はエルミート性を保ち、条件付き完全正値である；
3. $M_n$ 上の CP 写像 $\Phi$ と $G\in M_n$ が存在して、すべての $A$ に対し $\mathcal L(A)=\Phi(A)+GA+AG^*$ となる。

3 においては、$r\le n^2-1$、$\operatorname{Tr}L_j=0$、$j\neq k$ に対し $\operatorname{Tr}(L_j^*L_k)=0$ をみたす
$\Phi(A)=\sum_{j=1}^rL_jAL_j^*$ と、$\operatorname{Tr}G\in\mathbb R$ をみたす $G$ をとることができる。

*証明.* (1)$\Rightarrow$(2). 補題 8.1 により $t\to0$ のとき $(e^{t\mathcal L}-\mathrm{id})/t\to\mathcal L$ であり、$J$ は
有限次元空間上の線形写像なので連続である。よって
$J(e^{t\mathcal L})=J(\mathrm{id})+t\,J\big((e^{t\mathcal L}-\mathrm{id})/t\big)$ であり、最後のチョイ行列は
$J(\mathcal L)$ に収束する。$v\perp\Omega$ に対し $\langle v,J(\mathrm{id})v\rangle=|\langle\Omega,v\rangle|^2=0$ であるから
（補題 D.1(c)）、$t>0$ に対し、$J(e^{t\mathcal L})\ge0$（定理 6.2）より
$$
0\le\frac1t\langle v,J(e^{t\mathcal L})v\rangle=\big\langle v,J\big((e^{t\mathcal L}-\mathrm{id})/t\big)v\big\rangle
\;\xrightarrow[t\to0]{}\;\langle v,J(\mathcal L)v\rangle .
$$
各 $e^{t\mathcal L}$ は CP であり、したがってエルミート性を保つ（補題 D.1(d)）。よって
$\mathcal L(A^*)=\lim_t(e^{t\mathcal L}(A)^*-A^*)/t=\mathcal L(A)^*$ である。

(2)$\Rightarrow$(3). 補題 D.2 を $C=1$、$q_1=\Omega/\sqrt n$、すなわち $Q_1=1/\sqrt n$ として適用すると、
$G=K_1/\sqrt n$ として $\mathcal L(A)=\sum_jL_jAL_j^*+GA+AG^*$ となり、$\Phi=\sum_jL_j\cdot L_j^*$ は CP である
（定理 6.2）。さらに $\operatorname{Tr}L_j=\sqrt n\operatorname{Tr}(Q_1^*L_j)=0$、$r\le\operatorname{rank}P_W^\perp=n^2-1$ であり、
補題 D.1(b) により
$\operatorname{Tr}G=\langle\Omega,\xi_1\rangle/\sqrt n=\langle q_1,\xi_1\rangle=\frac12\langle q_1,J(\mathcal L)q_1\rangle\in\mathbb R$
である。$J(\mathcal L)$ が自己共役だからである。

(3)$\Rightarrow$(1). 定理 8.6 の証明のステップ 1–4 は、ある行列 $H_{\mathrm{eff}}$ により
$\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)$ として生成子が $\mathcal L_0+\Phi$ であり、
$\Phi$ が CP であることしか用いない。$H_{\mathrm{eff}}=iG$ とすれば $\mathcal L_0(\rho)=G\rho+\rho G^*$ であり、
これらのステップにより $e^{t\mathcal L}=\sum_{m\ge0}S_m(t)$ は CP である。$\square$

正値性だけでは不十分である。$n\ge2$ とし、$\theta$ を $M_n$ 上の転置、$\mathcal L=\theta-\mathrm{id}$ とする。$\theta^2=\mathrm{id}$ であるから
$$
e^{t\mathcal L}=e^{-t}\big(\cosh t\,\mathrm{id}+\sinh t\,\theta\big)=\tfrac{1+e^{-2t}}2\,\mathrm{id}+\tfrac{1-e^{-2t}}2\,\theta ,
$$
これは正値かつトレース保存な写像の凸結合である（命題 6.1）。しかし $J(\theta)=[E_{ji}]=\sum_{i,j}E_{ij}\otimes E_{ji}$ は
入れ替え $e_k\otimes e_l\mapsto e_l\otimes e_k$ であるから、$\Omega$ と直交する
$v=e_1\otimes e_2-e_2\otimes e_1$ を $-v$ に写す（命題 6.1 の証明と同様）。
$\|v\|^2=2$ と $\langle v,J(\mathrm{id})v\rangle=|\langle\Omega,v\rangle|^2=0$ より
$\langle v,J(\mathcal L)v\rangle=-2$ である：$\mathcal L$ は CCP ではなく、実際すべての $t>0$ に対し
$\langle v,J(e^{t\mathcal L})v\rangle=-(1-e^{-2t})<0$ である。正値かつトレース保存な写像からなる半群は
チャネルからなるとは限らない。

定理 D.3 を、`lindblad_superoperator` で構成した $M_3$ 上の任意のジャンプをもつ GKSL 生成子
（行優先の $\operatorname{vec}\rho$ に作用する $\mathcal L$ の行列 $S$ であり、6.1節と同様に
$J_{(i,a),(j,b)}=S_{(a,b),(i,j)}$ となる）と、$\theta-\mathrm{id}$ について確かめる：

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

## D.3 GKSL の定理

**定理 D.4（ゴリーニ–コサコフスキー–スダルシャン；リンドブラッド；定理 8.7）.** $M_n$ 上の超演算子 $\mathcal L$ が
$M_n$ 上の量子動力学的半群の生成子であるための必要十分条件は、$M_n$ の元 $H=H^*$ と
$L_1,\dots,L_r$ が存在して
$$
\mathcal L(\rho)=-i[H,\rho]+\sum_{j=1}^r\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
$$
となることである。$r\le n^2-1$、$\operatorname{Tr}H=0$、$\operatorname{Tr}L_j=0$、$j\neq k$ に対し $\operatorname{Tr}(L_j^*L_k)=0$ となるように選ぶことができる。

*証明.* 「ならば」の向きは定理 8.6 である。逆に、$(T_t)$ を生成子 $\mathcal L$ をもつ量子動力学的半群とすると、
$T_t=e^{t\mathcal L}$ である（補題 8.5）。各 $T_t$ は CP であるから、定理 D.3 により、トレースが零で互いに直交する
$r\le n^2-1$ 個の $L_j$ と $\operatorname{Tr}G\in\mathbb R$ をみたす $G$ により
$\mathcal L(\rho)=\sum_jL_j\rho L_j^*+G\rho+\rho G^*$ となる。$H=\frac i2(G-G^*)$、$R=\frac12(G+G^*)$ とおくと、
どちらも自己共役であり、$G=-iH+R$、$G^*=iH+R$、$G\rho+\rho G^*=-i[H,\rho]+\{R,\rho\}$ である。さらに
$\operatorname{Tr}H=\frac i2(\operatorname{Tr}G-\overline{\operatorname{Tr}G})=0$ である。各 $T_t$ は TP であるから
$\operatorname{Tr}\mathcal L(\rho)=\frac d{dt}\operatorname{Tr}T_t(\rho)\big|_{t=0}=0$ であり（補題 8.1）、巡回性により
$$
0=\operatorname{Tr}\mathcal L(\rho)=\operatorname{Tr}\Big(\rho\,\Big(\sum_jL_j^*L_j+2R\Big)\Big)\qquad\text{for all }\rho .
$$
自己共役な $\rho=X=\sum_jL_j^*L_j+2R$ をとると $\operatorname{Tr}(X^*X)=0$ となるので、$X=0$ であり
$R=-\frac12\sum_jL_j^*L_j$ である。これが GKSL 形である。$\square$

**定理 D.5（GKSL 形における自由度）.** $(H;L_1,\dots,L_r)$ と $(H';L'_1,\dots,L'_s)$ を GKSL
データ（$H,H'$ は自己共役）とし、定理 D.4 のように生成子 $\mathcal L$、$\mathcal L'$ を定める。

(a) 次のいずれの場合にも $\mathcal L'=\mathcal L$ である：
(i) $L'_j=L_j$ かつ $H'=H+c1$、$c\in\mathbb R$；
(ii) $H'=H$ であり、短い方の族を零で補ったうえで、ユニタリ $u$ により $L'_j=\sum_ku_{jk}L_k$；
(iii) $a_j\in\mathbb C$ により $L'_j=L_j+a_j1$ であり、$H'=H-\frac i2\sum_j(\bar a_jL_j-a_jL_j^*)$（これは自己共役である）。

(b) $\operatorname{Tr}H=\operatorname{Tr}H'=0$ で、すべての $L_j,L'_j$ のトレースが零であり、$\mathcal L=\mathcal L'$ ならば、$H=H'$ であり、
零で補ったうえで、ユニタリ $u$ により $L'_j=\sum_ku_{jk}L_k$ となる。

(c) したがって、二つの GKSL データが同じ生成子をもつための必要十分条件は、一方が他方から有限回の
ステップ (i)–(iii) で得られることである。特に、定理 D.4 の正規化された形における作用素 $H$ と、その
散逸部分 $\Phi(\rho)=\sum_jL_j\rho L_j^*$ は $\mathcal L$ により決まる。

*証明.* (a) (i) $[c1,\rho]=0$ である。(ii) 定理 6.7 により $\Phi'=\sum_jL'_j\cdot L_j'^*$ は
$\Phi=\sum_jL_j\cdot L_j^*$ に等しい。双対はクラウス作用素の選び方によらない（6.2節）ので、
$\sum_jL_j'^*L'_j=\Phi^*(1)=\sum_jL_j^*L_j$ である。(iii) 一つのジャンプ $L$ と $a\in\mathbb C$ に対し、展開すると
$$
(L+a)\rho(L+a)^*-\tfrac12\{(L+a)^*(L+a),\rho\}-\Big(L\rho L^*-\tfrac12\{L^*L,\rho\}\Big)=\tfrac12[B,\rho],
\qquad B=\bar aL-aL^* ,
$$
となる。最初の積から出る項 $\bar aL\rho+a\rho L^*$ が
$-\frac12(aL^*\rho+a\rho L^*+\bar aL\rho+\bar a\rho L)$ と組み合わさり、$|a|^2$ を含む項は打ち消し合うからである。
$B^*=-B$ であるから行列 $\frac i2B$ は自己共役であり、したがって $H'=H-\frac i2\sum_jB_j$ は自己共役であって、
$-i[H'-H,\rho]=-\frac12\sum_j[B_j,\rho]$ が散逸項の変化を打ち消す。

(b) $\Phi=\sum_jL_j\cdot L_j^*$、$G=-iH-\frac12\Phi^*(1)$ とおくと、$\mathcal L(\rho)=\Phi(\rho)+G\rho+\rho G^*$ である。
補題 D.1(c) により、$l_j=(1\otimes L_j)\Omega$、$g=(1\otimes G)\Omega$ として
$J(\mathcal L)=\sum_j|l_j\rangle\langle l_j|+|g\rangle\langle\Omega|+|\Omega\rangle\langle g|$ であり、
$\operatorname{Tr}L_j=0$ であるから $l_j\perp\Omega$ である（補題 D.1(b)）。$P^\perp$ を $\Omega^\perp$ の上への射影とすると、
$P^\perp J(\mathcal L)P^\perp=\sum_j|l_j\rangle\langle l_j|=J(\Phi)$ であり、同様に
$P^\perp J(\mathcal L')P^\perp=J(\Phi')$ である。よって $J(\Phi)=J(\Phi')$、$\Phi=\Phi'$ であり（補題 D.1(a)）、定理 6.7 が
ユニタリ $u$ を与える。すると $\Phi^*(1)=\Phi'^*(1)$ であり、したがってすべての $\rho$ に対し $[H-H',\rho]=0$ である。
すべての $E_{ij}$ と可換な行列 $X$ は $1$ の定数倍である（$XE_{ij}$ と $E_{ij}X$ の成分を比較せよ）から、
$H-H'=\lambda1$ であり、トレースをとれば $\lambda=0$ である。

(c) (a) により、ステップ (i)–(iii) は生成子を変えず、各ステップの逆は同じ種類のステップである
（(iii) については $-a_j$ をとる。$\bar a_jL'_j-a_jL_j'^*=\bar a_jL_j-a_jL_j^*$ だからである）。任意のデータは、
$a_j=-\operatorname{Tr}L_j/n$ とした (iii) に続けて、(iii) の後のハミルトニアンを $H''$ として $c=-\operatorname{Tr}H''/n$ とした (i)
を施すことで正規化された形になる。同じ生成子をもつ二つの正規化されたデータは、(b) により (ii) で結ばれる。$\square$

定理 D.4 の証明はアルゴリズムである：$J(\mathcal L)$ の $\Omega^\perp$ への圧縮から $\Phi$ を、次に $\xi_1$ から
$G$ を、そして $H$ を得る。これを上の生成子に対して実行し、ライブラリで生成子を再構成し、定理 D.5(iii) と (i)
で計算した元のデータの正規化された形と比較する：

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

**定理 D.6（定常状態の一意性）.** $\mathcal L$ をデータ $(H;L_1,\dots,L_r)$ をもつ定理 D.4 の GKSL 生成子とする。
$\mathcal L$ が忠実な定常状態 $\sigma$（$\sigma>0$）をもち、$H$ およびすべての $L_j$、$L_j^*$ と可換な行列が
$1$ の定数倍だけであるならば、$\ker\mathcal L=\mathbb C\sigma$ である：$\sigma$ は唯一の定常状態である。
さらに $\ker\mathcal L^*=\mathbb C1$ である。

*証明.* トレースによるペアリング $\operatorname{Tr}(\mathcal L(\rho)A)=\operatorname{Tr}(\rho\,\mathcal L^*(A))$ に関する双対は
$\mathcal L^*(A)=i[H,A]+\sum_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})$ である：巡回性により
$\operatorname{Tr}(-i[H,\rho]A)=\operatorname{Tr}(\rho\,i[H,A])$、$\operatorname{Tr}(L\rho L^*A)=\operatorname{Tr}(\rho L^*AL)$、
$\operatorname{Tr}(\{L^*L,\rho\}A)=\operatorname{Tr}(\rho\{L^*L,A\})$ だからである（これは第8章の演習問題 1 である）。
直接展開すると（$A\mapsto i[H,A]$ は導分であるから、ハミルトニアンの項は打ち消し合う）、すべての $X$ に対し
$$
\mathcal L^*(X^*X)-X^*\mathcal L^*(X)-\mathcal L^*(X^*)X=\sum_j[L_j,X]^*[L_j,X]:
$$
一つのジャンプについては、六つの反交換子の項のうち四つが対ごとに打ち消し合い、残りの二つ
$\frac12X^*L^*LX+\frac12X^*L^*LX$ は足し合わされる。よって左辺は $L^*X^*XL-X^*L^*XL-L^*X^*LX+X^*L^*LX$ となり、これは
$(XL-LX)^*(XL-LX)$ である。$X\in\ker\mathcal L^*$ とする。$\mathcal L^*$ はエルミート性を保つので
$\mathcal L^*(X^*)=\mathcal L^*(X)^*=0$ であり、したがって $\mathcal L^*(X^*X)=\sum_j[L_j,X]^*[L_j,X]=:Y\ge0$ かつ
$\operatorname{Tr}(\sigma Y)=\operatorname{Tr}(\mathcal L(\sigma)X^*X)=0$ である。$\sigma>0$ であるから、
$\sigma^{1/2}Y\sigma^{1/2}\ge0$ に対する $\operatorname{Tr}(\sigma Y)=\operatorname{Tr}(\sigma^{1/2}Y\sigma^{1/2})$ が零となるのは
$\sigma^{1/2}Y\sigma^{1/2}=0$、すなわち $Y=0$ のときに限る。トレースをとると、非負の項からなる和について
$\sum_j\operatorname{Tr}([L_j,X]^*[L_j,X])=0$ となるので、すべての $j$ に対し $[L_j,X]=0$ である。$X^*\in\ker\mathcal L^*$ に
同じ議論を適用すると $[L_j,X^*]=0$、すなわち $[L_j^*,X]=0$ を得る。すると $X$ は $L_j^*L_j$ と可換であり、
$\mathcal L^*(X)$ の散逸部分は消え、$0=\mathcal L^*(X)=i[H,X]$ となる。仮定により $X\in\mathbb C1$ である。逆に
$\mathcal L^*(1)=0$ である。よって $\ker\mathcal L^*=\mathbb C1$ である。ペアリング
$(\rho,A)\mapsto\operatorname{Tr}(\rho A)$ に関して互いに双対な $M_n$ の基底 $(E_{ij})$ と $(E_{ji})$ において、$\mathcal L^*$ の
行列は $\mathcal L$ の行列の転置である。したがって両者の階数は等しく、$\dim\ker\mathcal L=\dim\ker\mathcal L^*=1$ である。
$\mathbb C\sigma$ に属するトレース $1$ の状態は $\sigma$ である。$\square$

これはフリジェリオの判定条件 [*Commun. Math. Phys.* 63, 269 (1978)] の有限次元の場合であり、
8.5節で用いられる。どちらの仮定も必要である。第8章の演習問題 4 では、状態 $\tfrac12\,1$ は忠実かつ
定常であるが、$H$ と $\sigma_z$ はすべての対角行列と可換であり、定常状態は対角行列からなる
二次元の空間を張る。逆に、$M_3$ 上で $H=0$、$L_1=|e_1\rangle\langle e_3|$、
$L_2=|e_2\rangle\langle e_3|$ とする。$L_1,L_2$ が生成する $*$-環は $L_1L_2^*=E_{12}$ と $L_1^*L_1=E_{33}$ を含み、
したがってすべての行列単位を含むので、$L_j,L_j^*$ と可換なのは $1$ の定数倍だけである。しかし
$P=E_{11}+E_{22}$ に対し $\rho=P\rho P$ をみたすすべての $\rho$ は定常である。$L_j\rho=0=\rho L_j^*$ かつ
$\{E_{33},\rho\}=0$ だからである。よって $\dim\ker\mathcal L\ge4$ であり、定理 D.6 により、どの定常状態も
忠実ではない。上のランダムな生成子について、コードは仮定と結論を確かめ（行優先に並べた $A^T$ のベクトルに
作用する $\mathcal L^*$ の行列は $S^T$: $\operatorname{vec}(\mathcal L^*(A)^T)=S^T\operatorname{vec}(A^T)$ である）、さらに第二の例を確かめる：

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

## D.4 セクターをもつ環

$\mathcal M=\bigoplus_{c=1}^CM_{k_c}\subseteq M_D$ をブロック対角行列のなす環、$D=\sum_ck_c$、その極小中心射影を
$P_c$ とし、$E(X)=\sum_cP_cXP_c$ を $M_D$ から $\mathcal M$ の上への条件付き期待値とする。命題 8.8 の証明と同様に、
$B\in\mathcal M$ に対し $E(B)=B$、$E(XB)=E(X)B$、$E(BX)=BE(X)$ であり、$E(X^*)=E(X)^*$ である。
CP 写像 $\Phi:\mathcal M\to\mathcal M$ に対し、$\Phi^*$ はその blunt 双対である（命題 6.8）。

**定理 D.7.** $\mathcal L$ を $\mathcal M$ 上の量子動力学的半群 $(T_t)$ の生成子とする。このとき
CP 写像 $\Phi:\mathcal M\to\mathcal M$ と $H=H^*\in\mathcal M$ が存在して
$$
\mathcal L(\rho)=\Phi(\rho)-\tfrac12\{\Phi^*(1),\rho\}-i[H,\rho]\qquad(\rho\in\mathcal M).
$$
となる。命題 6.4 により、セクター間のクラウス作用素を用いて $\Phi(\rho)_d=\sum_{c,l}K^{dc}_l\rho_cK^{dc*}_l$ と書ける。
逆に、このような $\mathcal L$ はすべて $\mathcal M$ 上の量子動力学的半群を生成する（命題 8.8(1)）。
したがって 8.3節の生成子（そこでは $\Phi=\sum_jg_jL_j\cdot L_j^*+\sum_mg_m\Phi_m$）は、ちょうど $\mathcal M$ 上の
量子動力学的半群の生成子である。

*証明.* *ステップ1：$M_D$ 上の CP 写像.* 命題 6.4 により
$T_t(\rho)_d=\sum_{c,l}K^{dc}_l(t)\rho_cK^{dc}_l(t)^*$ である。$N^{dc}_l(t)\in M_D$ を、位置 $(d,c)$ にブロック
$K^{dc}_l(t)$ をもち、それ以外は零である行列とする。$X\in M_D$ に対し、$N^{dc}_lXN^{dc*}_l$ は位置 $(d,d)$ の
ブロック $K^{dc}_lX_{cc}K^{dc*}_l$ である。ここで $X_{cc}$ は第 $c$ 対角ブロック、すなわち $E(X)$ の第 $c$ ブロックである。
ゆえに $\hat T_t=T_t\circ E=\sum_{d,c,l}N^{dc}_l\cdot N^{dc*}_l$ は $M_D$ 上の CP 写像であり（定理 6.2）、
$\hat T_0=E$ である。補題 8.5 により $(\hat T_t-E)/t=\big((T_t-\mathrm{id})/t\big)\circ E\to\mathcal L\circ E$ である。

*ステップ2：セクターベクトルの外での正値性.* $p_c=(1\otimes P_c)\Omega\in\mathbb C^D\otimes\mathbb C^D$、
$W=\operatorname{span}\{p_c\}$ とする。補題 D.1(c) により $J(E)=\sum_c|p_c\rangle\langle p_c|$ であるから、$v\perp W$ と $t>0$ に対し
$0\le\frac1t\langle v,J(\hat T_t)v\rangle=\langle v,J((\hat T_t-E)/t)v\rangle\to\langle v,J(\mathcal L\circ E)v\rangle$ である。
各 $T_t$ はエルミート性を保つ（そのクラウス表示による）ので、$\mathcal L$ もそうであり（定理 D.3 と同様）、$E$ と
$\mathcal L\circ E$ もそうである。

*ステップ3：分解.* 補題 D.1(b) により $\langle p_c,p_d\rangle=\operatorname{Tr}(P_cP_d)=\delta_{cd}k_c$ であるから、
$Q_c=P_c/\sqrt{k_c}$ とした $q_c=p_c/\sqrt{k_c}=(1\otimes Q_c)\Omega$ は正規直交である。
$\Psi=\mathcal L\circ E$ に補題 D.2 を適用すると、$M_D$ 上の CP 写像 $\hat\Phi=\sum_j\hat L_j\cdot\hat L_j^*$ と $K_c\in M_D$ で
$$
\mathcal L(E(A))=\hat\Phi(A)+\sum_c\big(K_cAP_c+P_cAK_c^*\big)/\sqrt{k_c}\qquad(A\in M_D).
$$
をみたすものが得られる。$A\in\mathcal M$ に対しては $E(A)=A$ かつ $AP_c=P_cA$ であるから、和は
$\tilde G=\sum_cK_cP_c/\sqrt{k_c}$ として $\tilde GA+A\tilde G^*$ となる。

*ステップ4：圧縮.* $A\in\mathcal M$ に対し $\mathcal L(A)\in\mathcal M$ であるから、$G=E(\tilde G)\in\mathcal M$ として
$\mathcal L(A)=E(\mathcal L(A))=E(\hat\Phi(A))+GA+AG^*$ である。$\mathcal M$ 上の写像 $\Phi=E\circ\hat\Phi$ は
$\Phi(A)_d=\sum_{c,j}(P_d\hat L_jP_c)A_c(P_d\hat L_jP_c)^*$ である：これは命題 6.4 により CP であり、そのクラウス作用素は
$\hat L_j$ たちのブロックである。

*ステップ5：トレース保存.* 定理 D.4 と同様に $A\in\mathcal M$ に対し $\operatorname{Tr}\mathcal L(A)=0$ であり、
$\operatorname{Tr}\Phi(A)=\operatorname{Tr}(A\Phi^*(1))$ であるから、自己共役な $X=\Phi^*(1)+G+G^*\in\mathcal M$ について
すべての $A\in\mathcal M$ に対し $\operatorname{Tr}(AX)=0$ である。$A=X$ とすれば $X=0$ である（blunt トレースは忠実である。定理 3.1）。
$H=\frac i2(G-G^*)\in\mathcal M$ とおくと、$G=-iH-\frac12\Phi^*(1)$ かつ $G\rho+\rho G^*=-i[H,\rho]-\frac12\{\Phi^*(1),\rho\}$ である。$\square$

$\hat T_0$ が $\mathrm{id}$ ではなく $E$ であるため、ベクトル $p_c$ が $\Omega=\sum_cp_c$ に取って代わる：$J(\mathcal L\circ E)$ は
$\Omega^\perp$ 上で正であるとは限らず、下の例では実際に正ではない。

**系 D.8（ハイゼンベルク描像）.** $(S_t)_{t\ge0}$ を $\mathcal M$ 上の単位的 CP 写像の連続な族で、
$S_0=\mathrm{id}$ かつ $S_{s+t}=S_sS_t$ をみたすものとする。その生成子は、$\mathcal M$ 上の CP 写像 $\Psi$ と
$K\in\mathcal M$ で $\Psi(1)+K+K^*=0$ をみたすものにより $A\mapsto\Psi(A)+K^*A+AK$ と表される。

*証明.* $T_t=S_t^*$（blunt 双対）とおく。これらは CP（命題 6.8）かつ TP（系 6.9(b)）であり、
$T_0=\mathrm{id}$ であって、$T_{s+t}=(S_tS_s)^*=S_s^*S_t^*=T_sT_t$ である。
$\operatorname{Tr}(S_tS_s(x)y)=\operatorname{Tr}(xS_s^*S_t^*(y))$ かつ $S_tS_s=S_{s+t}$ だからである。$\Psi\mapsto\Psi^*$ は
有限次元空間上の線形写像であるから、$t\mapsto T_t$ は連続であり、$\mathcal L_S$ を $(S_t)$ の生成子として
$(T_t-\mathrm{id})/t\to\mathcal L_S^*$ である。定理 D.7 により $G+G^*=-\Phi^*(1)$ をみたす
$\mathcal L_S^*(\rho)=\Phi(\rho)+G\rho+\rho G^*$ が得られ、巡回性により
$\operatorname{Tr}((G\rho+\rho G^*)A)=\operatorname{Tr}(\rho(G^*A+AG))$ であるから、系 6.9(c) により
$\mathcal L_S(A)=\Phi^*(A)+G^*A+AG$ である。$\Psi=\Phi^*$（命題 6.8 により CP）と $K=G$ をとればよい。$\square$

これは定理 8.7 の後で引用したクリステンセンとエヴァンスの定理の有限次元の場合である。セクター内の
ジャンプとセクター間のランダムなジャンプをもつ `lindblad_rhs` の生成子について、$M_1\oplus M_2\subseteq M_3$ 上で
証明の五つのステップに沿って定理 D.7 を確かめる：

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

## D.5 作用素環による誤り訂正

10.3節の設定を用いる。$\Phi(\rho)=\sum_{i=1}^rE_i\rho E_i^*$ は $M_D$ 上のチャネル、$\mathcal C=P\mathbb C^D$ は
符号、$B(\mathcal C)=\{X\in M_D:X=PXP\}$ は $\mathcal C$ 上の作用素、$\mathfrak A\subseteq B(\mathcal C)$ は単位 $P$ を
もつ $*$-環（定理 2.2 により $\mathcal C$ 上のフォン・ノイマン環）、$\mathfrak A'$ は $B(\mathcal C)$ におけるその
可換子環である。すべての $A\in\mathfrak A$ に対し $P(\mathcal R\circ\Phi)^*(A)P=A$ となる $M_D$ 上のチャネル
$\mathcal R$ が存在するとき、$\mathfrak A$ は $\Phi$ に対して **訂正可能** であるという。$\mathcal R$ がクラウス作用素
$R_m$ をもつならば、$\mathcal R\circ\Phi$ はクラウス作用素 $R_mE_i$ をもち、$\sum_mR_m^*R_m=1$ である。

**定理 D.9（Bény–Kempf–Kribs；定理 10.7）.** $\mathfrak A$ が $\Phi$ に対して訂正可能であるための必要十分条件は、
すべての $i,j$ に対し $PE_i^*E_jP\in\mathfrak A'$ となることである。このとき、証明で構成するチャネル $\mathcal R$ が $\mathfrak A$ を訂正する。

*証明.* *必要性.* $\mathcal R$ が $\mathfrak A$ を訂正するとし、$W_{mi}=R_mE_iP$ とおくと、すべての $X\in M_D$ に対し
$\sum_{m,i}W_{mi}^*XW_{mi}=P(\mathcal R\circ\Phi)^*(X)P$ である。これは $X\in\mathfrak A$ に対しては $X$ に等しく、
$X=1$ に対しては $P$ に等しい。$(\mathcal R\circ\Phi)^*$ は単位的だからである（系 6.9(a)）。$p\in\mathfrak A$ を射影とする。
すると $p\le P$ であり、$\sum W_{mi}^*(1-p)W_{mi}=P-p$ である。$x=px$ に対しては
$\sum_{m,i}\|(1-p)W_{mi}x\|^2=\langle x,(P-p)x\rangle=0$ であり、$x=(P-p)x$ に対しては
$\sum_{m,i}\|pW_{mi}x\|^2=\langle x,px\rangle=0$ である。ゆえに $(1-p)W_{mi}p=0$ かつ $pW_{mi}(P-p)=0$ であり、
$W_{mi}=W_{mi}P$ であるから
$$
W_{mi}p=pW_{mi}p=pW_{mi}P=pW_{mi}.
$$
$\mathfrak A$ のすべての自己共役元 $A$ は $\mathfrak A$ の射影の線形結合である：$\mathcal C$ 上で定理 1.3 を用いると、
相異なる $\lambda_k$ と直交射影 $p_k$（$\sum_kp_k=P$）により $A=\sum_k\lambda_kp_k$ となり、
$p_k=\prod_{l\neq k}(A-\lambda_lP)/(\lambda_k-\lambda_l)\in\mathfrak A$ である（両辺を $\mathcal C$ の固有ベクトルと $\mathcal C^\perp$ に作用させよ。空積は $P$ とする）。
$\mathfrak A$ のすべての元は、自己共役な $A_1,A_2\in\mathfrak A$ により $A_1+iA_2$ と書ける。よってすべての
$A\in\mathfrak A$ に対し $W_{mi}A=AW_{mi}$ であり、随伴をとって $W_{mi}^*A=AW_{mi}^*$ である。
$\sum_mW_{mi}^*W_{mj}=PE_i^*(\sum_mR_m^*R_m)E_jP=PE_i^*E_jP$ であるから
$$
A\,PE_i^*E_jP=\sum_mW_{mi}^*AW_{mj}=PE_i^*E_jP\,A\qquad(A\in\mathfrak A),
$$
すなわち $PE_i^*E_jP\in\mathfrak A'$ である。

*十分性.* $\mathcal C$ 上の $\mathfrak A$ に定理 2.6 を適用すると、$\mathcal C$ の正規直交基底 $|c;a,b\rangle$
（$c=1,\dots,C$、$a=1,\dots,k_c$、$b=1,\dots,m_c$）——第 $c$ 直和成分における $e_a\otimes f_b$ のユニタリ $V$ による像——
が存在し、この基底において $\mathfrak A$ は作用素 $\sum_{c,a,a',b}(A_c)_{aa'}|c;a,b\rangle\langle c;a',b|$
（$A_c\in M_{k_c}$）からなり、$\mathfrak A'$ は作用素 $\sum_{c,a,b,b'}(Y_c)_{bb'}|c;a,b\rangle\langle c;a,b'|$（$Y_c\in M_{m_c}$）からなる。
したがって仮定により、次をみたす行列 $\Lambda^c_{ij}\in M_{m_c}$ が得られる：
$$
\langle c;a,b|E_i^*E_j|d;a',b'\rangle=\delta_{cd}\,\delta_{aa'}\,(\Lambda^c_{ij})_{bb'} .
$$
*ステップ1.* 各 $c$ と各組 $\alpha=(i,b)$ に対し $F^c_\alpha:\mathbb C^{k_c}\to\mathbb C^D$ を
$F^c_\alpha e_a=E_i|c;a,b\rangle$ で定め、$\beta=(j,b')$ に対し $\Lambda^c_{\alpha\beta}=(\Lambda^c_{ij})_{bb'}$ とおく。上の式は
$F^{c*}_\alpha F^d_\beta=\delta_{cd}\Lambda^c_{\alpha\beta}1_{k_c}$ を意味する：これは符号 $\mathbb C^{k_c}$ 上の「誤り」$F^c_\alpha$ に対する
クニル–ラフラム条件であり、加えてセクター間の直交性も含んでいる。

*ステップ2.* 定理 10.5 の証明と同様に、$\Lambda^c$ は自己共役である。
$\Lambda^c_{\beta\alpha}1=F^{c*}_\beta F^c_\alpha=(F^{c*}_\alpha F^c_\beta)^*=\bar\Lambda^c_{\alpha\beta}1$ だからである。また正である。
$\sum\bar z_\alpha z_\beta\Lambda^c_{\alpha\beta}1=(\sum z_\alpha F^c_\alpha)^*(\sum z_\beta F^c_\beta)\ge0$ だからである。
ユニタリ $U^c$ と $d^c_\kappa\ge0$ により $\Lambda^c=U^c\operatorname{diag}(d^c_\kappa)U^{c*}$ と書き、
$G^c_\kappa=\sum_\alpha u^c_{\alpha\kappa}F^c_\alpha$ とおく。すると $G^{c*}_\kappa G^d_\lambda=\delta_{cd}\delta_{\kappa\lambda}d^c_\kappa1$ であるから、
$d^c_\kappa=0$ ならば $G^c_\kappa=0$ であり、そうでなければ $V^c_\kappa=G^c_\kappa/\sqrt{d^c_\kappa}$ は
$V^{c*}_\kappa V^d_\lambda=\delta_{cd}\delta_{\kappa\lambda}1$ をみたす等長作用素である。逆に解くと $F^c_\alpha=\sum_\kappa\bar u^c_{\alpha\kappa}G^c_\kappa
=\sum_{\kappa:d^c_\kappa>0}\bar u^c_{\alpha\kappa}\sqrt{d^c_\kappa}\,V^c_\kappa$ であり、ゆえに
$$
V^{c'*}_\kappa F^d_\beta=\delta_{c'd}\,\bar u^d_{\beta\kappa}\sqrt{d^d_\kappa}\,1,\qquad
F^{c*}_\alpha V^{c'}_\kappa=\delta_{cc'}\,u^c_{\alpha\kappa}\sqrt{d^c_\kappa}\,1 .
$$
*ステップ3：回復.* $V^c_\kappa V^{c*}_\kappa$（$d^c_\kappa>0$）は互いに直交する射影である。$Q$ をそれらの
和とする。等長作用素 $\iota_c:\mathbb C^{k_c}\to\mathbb C^D$ を $\iota_ce_a=|c;a,1\rangle$ で定め、
$$
\mathcal R(\sigma)=\sum_{c}\sum_{\kappa:d^c_\kappa>0}\iota_cV^{c*}_\kappa\,\sigma\,V^c_\kappa\iota_c^*+(1-Q)\sigma(1-Q).
$$
とおく。これはチャネルである：$\sum_{c,\kappa}V^c_\kappa\iota_c^*\iota_cV^{c*}_\kappa+(1-Q)^2=Q+1-Q=1$。これはセクター $c$ と
シンドローム $\kappa$ を測定し、$\mathbb C^{k_c}$ 上の誤りを元に戻し、結果をコピー $b=1$ に書き込む。

*ステップ4：正しさ.* $\mathcal R^*(X)=\sum_{c,\kappa}V^c_\kappa\iota_c^*X\iota_cV^{c*}_\kappa+(1-Q)X(1-Q)$ であり、
ブロック $A_c$ をもつ $A\in\mathfrak A$ に対し $\iota_c^*A\iota_c=A_c$ である。各 $F^d_\beta$ は $V^d_\kappa$ たちの線形結合であるから
$(1-Q)F^d_\beta=0$ である。ステップ2とあわせて、$\alpha=(i,b)$、$\beta=(j,b')$ に対し
$$
\langle c;a,b|E_i^*\mathcal R^*(A)E_j|d;a',b'\rangle=\sum_{c',\kappa}\langle e_a,F^{c*}_\alpha V^{c'}_\kappa A_{c'}V^{c'*}_\kappa F^d_\beta e_{a'}\rangle
=\delta_{cd}\sum_\kappa u^c_{\alpha\kappa}d^c_\kappa\bar u^c_{\beta\kappa}(A_c)_{aa'}=\delta_{cd}(\Lambda^c_{ij})_{bb'}(A_c)_{aa'} .
$$
$j=i$ とおいて $i$ について和をとる：$\sum_iE_i^*E_i=1$ と上の式により $\sum_i\Lambda^c_{ii}=1_{m_c}$ であるから
$$
\langle c;a,b|(\mathcal R\circ\Phi)^*(A)|d;a',b'\rangle=\delta_{cd}\,\delta_{bb'}(A_c)_{aa'}=\langle c;a,b|A|d;a',b'\rangle .
$$
$P(\mathcal R\circ\Phi)^*(A)P$ と $A$ はともに $B(\mathcal C)$ に属し、基底 $|c;a,b\rangle$ に関して同じ行列をもつので、
等しい。$\square$

**系 D.10.** (a) $\mathfrak A=B(\mathcal C)$ に対しては $\mathfrak A'=\mathbb CP$ であり、定理 D.9 はクニル–ラフラムの
定理 10.5 である（$C=1$、$m_1=1$。ステップ3は定理 10.5 の回復である）。(b) 互いに直交する射影 $P_c\ne0$、
$\sum_cP_c=P$ による $\mathfrak A=\operatorname{span}\{P_1,\dots,P_C\}$ に対しては、これは命題 10.6 における 2 と 3 の同値性である。
(c)（*ノイズレス部分系.*）命題 10.1(a) の記号で、$\mathcal C=V(\mathbb C^{k_c}\otimes\mathbb C^{m_c})$ を $\Phi$ の
相互作用環 $\mathcal A$ に関する第 $c$ 直和成分とし、$\mathcal C$ 上で $\mathfrak A=V(1_{k_c}\otimes M_{m_c})V^*$ とする。これは
二つの因子を入れ替えれば $M_{m}\otimes1_{k}$ の形の環である。このとき $\mathfrak A$ は $\mathcal R=\mathrm{id}$ により訂正され、
$PE_i^*E_jP=V(E_{i,c}^*E_{j,c}\otimes1_{m_c})V^*\in\mathfrak A'$ ($P=VV^*$) である。

*証明.* (a) $\mathfrak A'=\mathbb CP$ は上で述べた $\mathfrak A'$ の記述の $C=1$、$m_1=1$ の場合であり、
$PE_i^*E_jP\in\mathbb CP$ はクニル–ラフラム条件である。$B(\mathcal C)$ の訂正可能性は 10.2節の意味での訂正である：
$\mathcal C$ に台をもつすべての状態 $\rho$ に対し $\mathcal R(\Phi(\rho))=\rho$ ならば、それらすべてに対し
$\operatorname{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(A)P)=\operatorname{Tr}(\rho A)$ であり、これらの $\rho$ は $B(\mathcal C)$ を張る。
逆は 10.3節で示されている。$C=1$、$m_1=1$ のとき、ステップ3のクラウス作用素 $\iota_1V^{1*}_\kappa$ は
定理 10.5 のものである。$\iota_1\iota_1^*=P$ だからである。(b) $T\in B(\mathcal C)$ がすべての $P_c$ と可換であるための
必要十分条件は、$c\ne d$ に対し $P_cTP_d=0$ となることであり、これは命題 10.6 の条件 2 である。
(c) すべての $A\in\mathfrak A$ は、他の直和成分上で $0$ と拡張すれば
$\mathcal A'=V(\bigoplus_c1_{k_c}\otimes M_{m_c})V^*$ に属する（定理 2.6）。よって命題 10.1(a) により $\Phi^*(A)=A$ であり、
$PAP=A$ である。$PE_i^*E_jP$ の式は $V^*E_iV=\bigoplus_cE_{i,c}\otimes1_{m_c}$ から従い、
$E_{i,c}^*E_{j,c}\otimes1_{m_c}$ は $1_{k_c}\otimes M_{m_c}$ と可換である。$\square$

定理 2.6 により、すべての $\mathfrak A$ はある $\bigoplus_cM_{k_c}$ と $*$-同型である。$C\ge2$ で、ある $k_c\ge2$ であるとき、
訂正可能な $\mathfrak A$ はセクターのラベル $c$（古典情報）を、各 $M_{k_c}$ の量子情報とともに保護する。
このような符号は **ハイブリッド** 符号と呼ばれる。例として、$\mathbb C^8$ 内のランダムな五次元の符号上の
$\mathfrak A=(M_2\otimes1_2)\oplus M_1$ と、次のようなノイズをとる：ゲージ因子 $\mathbb C^{m_c}$ に
$\sum_iB_{c,i}^*B_{c,i}=1_{m_c}$ を通じて作用し、セクターを互いに直交する部分空間に写し、$\mathcal C^\perp$ をかき混ぜ、
その後ランダムなユニタリによって混ぜ合わされる（定理 6.7）。クニル–ラフラム条件は成り立たず、定理 D.9 の条件は
成り立ち、証明の回復はすべての $A\in\mathfrak A$ を復元するが、$\mathcal C$ 上の一般の作用素は復元しない：

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

本付録の結果は、V. Gorini, A. Kossakowski, E. C. G. Sudarshan
[*J. Math. Phys.* 17, 821 (1976)] と G. Lindblad [*Commun. Math. Phys.* 48, 119 (1976)]（定理 D.3–D.5）、
E. Christensen と D. E. Evans [*J. London Math. Soc.* 20, 358 (1979)]（一般のフォン・ノイマン環における系 D.8）、
A. Frigerio（定理 D.6）、そして C. Bény, A. Kempf, D. W. Kribs [*Phys. Rev. Lett.* 98, 100502
(2007)；*Phys. Rev. A* 76, 042303 (2007)]（定理 D.9）による。上の証明は、本講義で用いる有限次元の主張に対しては完全である。
