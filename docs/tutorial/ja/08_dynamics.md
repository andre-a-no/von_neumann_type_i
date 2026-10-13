[English](../08_dynamics.md) · [Español](../es/08_dynamics.md) · [Français](../fr/08_dynamics.md) · [Deutsch](../de/08_dynamics.md) · [中文](../zh/08_dynamics.md) · **日本語** · [Русский](../ru/08_dynamics.md)

# 8. ダイナミクス

閉じた量子系はユニタリの一径数群によって時間発展し、開いた量子系——観測されない環境と
結合した系——はチャネルの半群によって時間発展する。セクターをもつ環では、どちらの発展にも
単一の行列環の上では見えない特徴がある：$\mathcal M$ に属するハミルトニアンは決してセクターを変えられず、
したがって電荷は保存されるが、開いた系は、たとえば粒子を失うことによって、セクター間をジャンプしうる。
本章では両種の時間発展についての基本的事実を証明し、セクター間のジャンプを含むリンドブラッド方程式を書き下し、
「粒子が一つも失われなかった」ことが電荷の測定であって、その事後選択された状態が非エルミートな
ハミルトニアンによって発展することを示す。

**本章で学ぶこと**：セクターの重みが運動の定数である理由、
ゴリーニ–コサコフスキー–スダルシャン–リンドブラッド（GKSL）形の生成子がチャネルを生み出す理由、セクター間のジャンプが
リンドブラッド方程式にどのように入り、セクターの重みがいつ古典的なレート方程式に従うか、ジャンプなしの
時間発展と事後選択が何を意味するか、そして定常状態と緩和レートの求め方。

本章を通じて、$\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ は
$\mathcal H=\bigoplus_c\mathbb C^{k_c}$ に作用し、$D=\sum_ck_c$、$P_c$ は極小中心射影とする。
状態はその blunt 密度 $\rho$ で与え、セクターの重みを $p_c=\operatorname{Tr}\rho_c$ とする
（第3章）。また $[A,B]=AB-BA$、$\{A,B\}=AB+BA$、$\hbar=1$ とする。**超演算子** とは $\mathcal M$ または
$M_D(\mathbb C)$ 上の線形写像のことである。超演算子のノルムはすべて作用素ノルムとし、したがって
劣乗法的である。有限次元空間上の線形写像 $B$ に対し $e^{tB}=\sum_{n\ge0}t^nB^n/n!$ とする。

## 8.1 閉じた系

**ハミルトニアン** とは自己共役な $H\in\mathcal M$ のことである。ベクトルは **シュレーディンガー
方程式** $i\psi'(t)=H\psi(t)$ に従って、密度は **フォン・ノイマン方程式** $\rho'(t)=-i[H,\rho(t)]$ に従って発展する。
どちらも有限次元空間上の定数係数線形方程式であり、
次の補題が我々に必要な微分方程式論のすべてである。

**補題 8.1.** $B$ を有限次元ノルム空間 $X$ 上の線形写像とする。このとき
$\frac d{dt}e^{tB}=Be^{tB}=e^{tB}B$、$e^{(s+t)B}=e^{sB}e^{tB}$、$\|e^{tB}\|\le e^{|t|\,\|B\|}$ であり、
任意の $x_0\in X$ に対し、写像 $x(t)=e^{tB}x_0$ は、区間
$I\ni0$ 上で $x'=Bx$ かつ $x(0)=x_0$ を満たす微分可能な $x$ としてただ一つのものである。

*証明.* 級数は上のノルム評価とともに絶対収束し、$t$ のべき級数として
項別微分できるので、第一の公式が得られる。すると
$\frac d{dt}(e^{tB}e^{-tB})=Be^{tB}e^{-tB}-e^{tB}Be^{-tB}=0$ であるから $e^{tB}e^{-tB}=1$ である。$x'=Bx$ ならば
$\frac d{dt}(e^{-tB}x(t))=e^{-tB}(x'-Bx)=0$ であるから、$e^{-tB}x(t)=x_0$、すなわち $x(t)=e^{tB}x_0$ である。最後に
$t\mapsto e^{(s+t)B}x_0$ と $t\mapsto e^{tB}e^{sB}x_0$ はともに同じ初期値をもつ $y'=By$ の解である。$\square$

**命題 8.2.** $H=H^*\in\mathcal M$ とし、$U(t)=e^{-itH}$ を 1.3節の関数解析により
$f_t(\lambda)=e^{-it\lambda}$ として定める。

1. $U(t)=\bigoplus_ce^{-itH_c}$ は $\mathcal M$ のユニタリであり、$U(s+t)=U(s)U(t)$ かつ $U'(t)=-iHU(t)$ である。
   $U(t)$ は $-itH$ の指数級数に等しい。
2. $\psi(t)=U(t)\psi_0$ は $\psi(0)=\psi_0$ を満たす $i\psi'=H\psi$ のただ一つの解であり、
   $\rho(t)=U(t)\rho_0U(t)^*$ は
   $\rho(0)=\rho_0$ を満たす $\rho'=-i[H,\rho]$ の $\mathcal M$ におけるただ一つの解である。$\rho_0$ が状態ならば、すべての $\rho(t)$ も状態である。
3. （*ハイゼンベルク描像*）
   $A(t)=U(t)^*AU(t)$ として $\operatorname{Tr}(\rho(t)A)=\operatorname{Tr}(\rho_0A(t))$ であり、$A'=i[H,A]$ である。

*証明.* (1) 各ブロックを対角化する：$H_c=V_c\Lambda_cV_c^*$。すると $V=\bigoplus_cV_c$ は
$H$ を対角化し、$U(t)=Ve^{-it\Lambda}V^*$ はブロック対角である（定理 1.5(a) により、結果は
対角化に用いるユニタリによらない）。$|f_t|^2=1$ かつ $f_sf_t=f_{s+t}$ であるから、定理 1.5(c) によりユニタリ性と
群法則が得られる。$e^{-it\Lambda}$ を成分ごとに微分すると $U'=V(-i\Lambda)e^{-it\Lambda}V^*=-iHU$ であるから、
補題 8.1（$X=M_D$、$B$ は $-iH$ による左からの乗法）により $U(t)$ は級数 $e^{-itH}$ である。
(2) $(U\psi_0)'=-iHU\psi_0$ であり、$(U')^*=iU^*H$ であるから $(U\rho_0U^*)'=-iHU\rho_0U^*+U\rho_0U^*(iH)$ である。
一意性は補題 8.1 を、$\mathcal H$ 上の $B=-iH$、あるいは $\mathcal M$ 上の $B=-i[H,\cdot\,]$ に適用したものである。正値性：
$\langle x,U\rho_0U^*x\rangle=\langle U^*x,\rho_0U^*x\rangle\ge0$。トレースは巡回性による。(3) トレースの
巡回性と、$(U^*AU)'=iU^*HAU-iU^*AHU=i[H,A(t)]$ による。$\square$

**定理 8.3（電荷の保存）.** $t\mapsto H(t)$ を区間 $I$ 上で $\mathcal M$ の自己共役元に値をとる
写像とする。

1. $\psi:I\to\mathcal H$ が微分可能で $i\psi'=H(t)\psi$ を満たすならば、すべての $\|P_c\psi(t)\|$ は一定である。
2. $\rho:I\to\mathcal M$ が微分可能で $\rho'=-i[H(t),\rho]$ を満たすならば、すべての中心元 $Z$ に対し $\operatorname{Tr}(Z\rho(t))$ は
   一定である。特に、すべてのセクターの重み $p_c(t)=\operatorname{Tr}(P_c\rho(t))$ は一定である。
3. $H$ が定数で $A\in\mathcal M$ が $H$ と可換ならば、$\operatorname{Tr}(A\rho(t))$ は一定である。

*証明.* (1) 内積は第一引数について反線形であるから、
$\frac d{dt}\langle\psi,P_c\psi\rangle=\langle-iH\psi,P_c\psi\rangle+\langle\psi,-iP_cH\psi\rangle
=i\langle\psi,HP_c\psi\rangle-i\langle\psi,P_cH\psi\rangle=0$ である。$P_c$ は中心に属するからである。(2)
巡回性により $\frac d{dt}\operatorname{Tr}(Z\rho)=-i\operatorname{Tr}(Z[H,\rho])=-i\operatorname{Tr}([Z,H]\rho)=0$ である。
(3) $A$ について同じ計算をすればよい。$\square$

> **物理的意味.** $\mathcal M$ に属するハミルトニアンは電荷（粒子数、
> 磁化、対称性の量子数）と可換であるから、各電荷の確率は運動の
> 定数である。これは駆動された時間依存のハミルトニアンに対しても成り立つ。計算の上では、時間発展は
> $C$ 個の独立な問題に分かれる：一回の対角化のコストは $D^3$ ではなく $\sum_ck_c^3$ 程度の演算数である。

ライブラリはブロックサイズごとに一回の固有値分解で $U(t)$ を計算し（`dynamics.propagator`）、
`dynamics.von_neumann` は固有基底において密度を多数の時刻へ一度に発展させる。固有基底では
$\rho(t)_{ij}=e^{-it(\lambda_i-\lambda_j)}\rho_{ij}$ である。`scipy.linalg.expm` と照合する：

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

$\mathcal M$ 上では、ベクトル $\psi=\bigoplus_c\psi_c$ は混合
$\bigoplus_c|\psi_c\rangle\langle\psi_c|$ の状態を定める（第3章の演習問題 3）。したがって `dynamics.schrodinger`（ブロック配置
`(batch, C, k_max)` のベクトル）と `von_neumann` はその上で一致しなければならない。時間依存の $H(t)$ に対しては、
ライブラリは固定刻みの古典的ルンゲ–クッタ法（`dynamics.rk4`、`schrodinger_rk4`）を用いる：一つのステップは
$\psi$ に $H(t)\in\mathcal M$ のいくつかの値の多項式を掛けるので、セクター間で振幅を移すことは決してないが、
厳密にはユニタリでなく、セクターの重みは離散化誤差の程度でドリフトする。
（一つの大きなセクターに対しては、`krylov.evolve` は行列–ベクトル積しか必要としない。第9章を参照。）

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

セクターの重みは RK4 のもとでは $1.1\cdot10^{-9}$ だけ変化し、厳密なソルバーのもとでは
丸め誤差を超える変化はない。

## 8.2 量子動力学的半群

CP 写像、TP 写像およびチャネルは 6.1節（[第6章](06_quantum_channels.md)）で定義した。
定理 6.2 と命題 6.4（写像 $A\mapsto\sum_lK_lAK_l^*$ は、セクター間のものも含めて CP である）のほかに、
三つの閉性が必要である。

**補題 8.4.** (a) CP 写像の和、非負の定数倍、合成は CP である。
$\Phi$ が $M_D$ 上の CP 写像で、$*$-部分環 $\mathcal N\subseteq M_D$（たとえば $\mathcal M$、あるいは
単位元 $P_c$ をもつ $P_c\mathcal MP_c\cong M_{k_c}$）をそれ自身に写すならば、$\Phi|_{\mathcal N}$ は CP である。
(b) CP 写像 $\Phi_m$ が $\Phi$ に収束する（すなわちすべての $A$ に対し $\Phi_m(A)\to\Phi(A)$）ならば、$\Phi$ は CP である。
(c) $s\mapsto\Psi(s)$ が $[a,b]$ から CP 写像への連続写像ならば、$\int_a^b\Psi(s)\,ds$ は CP である。

*証明.* (a) $\Phi\mapsto\Phi_n$ は線形であり、$(\Phi\Psi)_n=\Phi_n\Psi_n$ である。正の行列の和は
正である。制限について：$M_n(\mathcal N)\subseteq M_n(M_D)$ であり、正値性は両者で同じ
条件である（基底を並べ替えると、$M_n(\mathcal N)$ の元は
ブロック $[(A_{ij})_c]$ と零の直和になる）。(b) 正の行列全体は閉集合をなし（各
$x$ に対する $\langle x,Xx\rangle\ge0$ は閉条件である）、$(\Phi_m)_n(X)\to\Phi_n(X)$ である。(c) リーマン和は (a) により CP であり、
$\Psi$ は一様連続であるから積分に収束する。(b) を適用すればよい。$\square$

**定義.** $\mathcal M$ 上の **量子動力学的半群** とは、
$T_0=\mathrm{id}$、$T_{s+t}=T_sT_t$ を満たし、$t\mapsto T_t$ が連続であるようなチャネル $\mathcal M\to\mathcal M$ の族 $(T_t)_{t\ge0}$ のことである。

**補題 8.5.** 超演算子の任意の連続な半群（$T_0=\mathrm{id}$、$T_{s+t}=T_sT_t$）は、ただ一つの超演算子 $\mathcal L$ により
$T_t=e^{t\mathcal L}$ と書ける。$\mathcal L$ をその **生成子** といい、$\mathcal L=\lim_{t\to0}(T_t-\mathrm{id})/t$ である。

*証明.* $h\to0$ のとき $V_h=\frac1h\int_0^hT_s\,ds\to\mathrm{id}$ であるから、ある固定した
小さな $h$ に対し $V_h$ は可逆である。ここで $T_tV_h=\frac1h\int_t^{t+h}T_s\,ds$ は $t$ について微分可能で、導関数は
$\frac1h(T_{t+h}-T_t)=T_t\frac{T_h-\mathrm{id}}h$ である。したがって $T_t$ は微分可能で、
$\mathcal L=\frac{T_h-\mathrm{id}}hV_h^{-1}$ として $T_t'=T_t\mathcal L$ であり、補題 8.1（$B$ は
$\mathcal L$ を右から合成する写像）により $T_t=e^{t\mathcal L}$ である。一意性：$\mathcal L=T'(0)$。$\square$

**定理 8.6（GKSL 生成子）.** $H=H^*$ および $L_1,\dots,L_J$ を $M_D$ の元とし、
$$
\mathcal L(\rho)=-i[H,\rho]+\sum_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
$$
とする。このとき $(e^{t\mathcal L})_{t\ge0}$ は $M_D$ 上の量子動力学的半群である。より詳しく、
$H_{\mathrm{eff}}=H-\frac i2\sum_jL_j^*L_j$、$\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)$、
$\Phi(\rho)=\sum_jL_j\rho L_j^*$ とおいて $\mathcal L=\mathcal L_0+\Phi$ とし、
$$
S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},\qquad
S_n(t)=\int_0^tS_0(t-s)\,\Phi\,S_{n-1}(s)\,ds\quad(n\ge1).
$$
と定める。このときすべての $S_n(t)$ は CP であり、$e^{t\mathcal L}=\sum_{n\ge0}S_n(t)$（**ダイソン展開**）であって、
$\|S_n(t)\|\le e^{t\|\mathcal L_0\|}(t\|\Phi\|)^n/n!$ が成り立つ。

*証明.* *ステップ 1.* $X(t)=e^{-itH_{\mathrm{eff}}}\rho e^{itH_{\mathrm{eff}}^*}$ は $X'=\mathcal L_0(X)$
と $X(0)=\rho$ を満たすので、補題 8.1 により $S_0(t)=e^{t\mathcal L_0}$ である。これはただ一つのクラウス作用素
$e^{-itH_{\mathrm{eff}}}$ をもつ CP 写像である。*ステップ 2.* $a=\|\mathcal L_0\|$、$b=\|\Phi\|$ とおく。帰納法により $S_n$
は連続で、
$\|S_n(t)\|\le\int_0^te^{(t-s)a}\,b\,e^{sa}\frac{(bs)^{n-1}}{(n-1)!}ds=e^{ta}\frac{(bt)^n}{n!}$ であるから、
級数は任意の $[0,T]$ 上で一様収束する。*ステップ 3.*
$S_n(t)=S_0(t)\int_0^te^{-s\mathcal L_0}\Phi S_{n-1}(s)\,ds$ であるから、積の微分法則と微積分学の基本
定理により $S_n'=\mathcal L_0S_n+\Phi S_{n-1}$（$n\ge1$）および $S_0'=\mathcal L_0S_0$ を得る。上の
評価により導関数の級数は $[0,T]$ 上で一様収束するので、$Y(t)=\sum_nS_n(t)$ は
微分可能で $Y'=\mathcal L_0Y+\Phi Y=\mathcal LY$、$Y(0)=\mathrm{id}$ であり、補題 8.1 により $Y(t)=e^{t\mathcal L}$ である。
*ステップ 4.* $\Phi$ は CP であるから（定理 6.2）、帰納法により被積分関数
$s\mapsto S_0(t-s)\Phi S_{n-1}(s)$ は CP 写像の連続な族であり、$S_n(t)$ は CP である（補題 8.4）。
部分和とその極限 $e^{t\mathcal L}$ も同様である。*ステップ 5.* 巡回性により
$\operatorname{Tr}\mathcal L(\rho)=\sum_j\big(\operatorname{Tr}(L_j\rho L_j^*)-\operatorname{Tr}(L_j^*L_j\rho)\big)=0$ であり、
したがって $\frac d{dt}\operatorname{Tr}e^{t\mathcal L}\rho=\operatorname{Tr}\mathcal L(e^{t\mathcal L}\rho)=0$ である：すべての
$e^{t\mathcal L}$ は TP である。半群法則と連続性は補題 8.1 による。$\square$

**定理 8.7（ゴリーニ–コサコフスキー–スダルシャン；リンドブラッド；証明なしで述べる）.** 逆に、$M_D$ 上の任意の量子
動力学的半群の生成子は定理 8.6 の形をもつ。[V. Gorini, A. Kossakowski,
E. C. G. Sudarshan, *J. Math. Phys.* 17, 821 (1976); G. Lindblad, *Commun. Math. Phys.* 48, 119
(1976).] 一般のフォン・ノイマン環 $\mathcal M$ 上の単位的 CP 写像からなるノルム連続な半群については、
クリステンセンとエヴァンス [*J. London Math. Soc.* 20, 358 (1979)] が、生成子が
$\Psi$ を CP、$K\in\mathcal M$ として $A\mapsto\Psi(A)+K^*A+AK$ の形であることを示した。

$\mathcal M=\bigoplus_cM_{k_c}$ 上では、これをハイゼンベルク描像の半群 $(T_t^*)$ に適用する。その
写像は CP かつ単位的であり（命題 6.8、系 6.9(b)）、生成子は $\mathcal L^*$ である。
blunt 双対をとって戻ると、$\mathcal L(\rho)=\Phi(\rho)+G\rho+\rho G^*$ となる。ここで $\Phi=\Psi^*$ は
$\mathcal M$ 上の CP 写像（セクター間のクラウス作用素、命題 6.4）であり、$G=K$ である。トレース保存、すなわちすべての $\rho\in\mathcal M$ に対する $\operatorname{Tr}\big((\Phi^*(1)+G+G^*)\rho\big)=0$ から
$G+G^*=-\Phi^*(1)$、すなわち $H=H^*\in\mathcal M$ により $G=-iH-\frac12\Phi^*(1)$ となる：これはちょうど
次節の生成子である。

> **物理的意味.** 相関が速く減衰する大きな環境と弱く相互作用する系
> （ボルン–マルコフ近似）は、過去の記憶を失う：時刻 $s+t$ における状態は、
> 時刻 $s$ における状態から $t$ だけに依存するチャネルによって得られる。定理 8.7 は、
> そのようなダイナミクスが常にハミルトニアン部分と「ジャンプ」$L_j$（光子放出、粒子損失、
> 位相緩和）の和であり、そのレートは $L_j$ の規格化によって決まることを述べている。

ライブラリには四つのリンドブラッド・ルーチンがある：`lindblad_rhs` は右辺を
`rk4` 用の関数として返す。`lindblad_evolve` は RK4 で積分する。`lindblad_superoperator` は各セクターについて、
行優先の $\operatorname{vec}\rho_c$ に作用する $\mathcal L$ の $k_c^2\times k_c^2$ 行列を返す。そして
`lindblad_channel` は $e^{t\mathcal L}$ をセクターごとに厳密に計算し、`Channel` として返す。
$M_3$ 上の任意のジャンプについて定理 8.6 を確かめる：$e^{t\mathcal L}$ のチョイ行列は正であり（定理
6.2）、写像は TP であり、半群法則が成り立ち、RK4 は指数関数と一致し、ジャンプ部分
$e^{t\mathcal L}-S_0(t)=\sum_{n\ge1}S_n(t)$ それ自体が CP である：

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

## 8.3 セクターをもつ環におけるリンドブラッド方程式

$\mathcal M$ 上には二種類のジャンプがある。**セクター内** のジャンプは作用素
$L_j\in\mathcal M$ である。**セクター間** のジャンプは、クラウス
作用素 $K^{dc}_l:\mathbb C^{k_c}\to\mathbb C^{k_d}$ をもつ CP 写像 $\Phi:\mathcal M\to\mathcal M$ である（命題 6.4、`InterSectorChannel`）：
$\Phi(\rho)_d=\sum_{c,l}K^{dc}_l\rho_cK^{dc*}_l$。命題 6.8 により、その blunt 双対は
$\Phi^*(1)_c=\sum_{d,l}K^{dc*}_lK^{dc}_l$ を満たし、これは $\mathcal M$ のブロック対角な元である。レート
$g_j,g_m>0$ をもつ `lindblad_rhs` の生成子は
$$
\mathcal L(\rho)=-i[H,\rho]+\sum_jg_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big)
+\sum_mg_m\Big(\Phi_m(\rho)-\tfrac12\{\Phi_m^*(1),\rho\}\Big).
$$
である。

**命題 8.8.** $J^{dc}_l\in M_D$ を、位置 $(d,c)$ にブロック $K^{dc}_l$ をもち、
それ以外は零である行列とする。

1. $\mathcal L$ は $\mathcal M$ を $\mathcal M$ に写し、ジャンプ $\sqrt{g_j}L_j$ と $\sqrt{g_m}J^{dc}_{m,l}$（すべての $m,l,d,c$）をもつ $M_D$ 上の GKSL 生成子
   $\hat{\mathcal L}$ の $\mathcal M$ への制限である。したがって
   $(e^{t\mathcal L})_{t\ge0}$ は $\mathcal M$ 上の量子動力学的半群である。
2. 代わりに $\tilde{\mathcal L}$ を、一つのクラウス
   添字のブロックを足し合わせて一つのジャンプ $L_{m,l}=\sum_{d,c}J^{dc}_{m,l}$ とした $M_D$ 上の GKSL 生成子とする。このとき
   $\rho\in\mathcal M$ に対し $\mathcal L(\rho)=E\big(\tilde{\mathcal L}(\rho)\big)$ である。ここで $E(X)=\sum_dP_dXP_d$ は
   $M_D$ から $\mathcal M$ の上への条件付き期待値である。一般には $\tilde{\mathcal L}(\rho)\notin\mathcal M$ である。

*証明.* $\rho=\sum_cP_c\rho P_c$ かつ $J^{dc}P_{c'}=\delta_{cc'}J^{dc}$ であるから、
$J^{dc}\rho J^{dc*}$ は位置 $(d,d)$ にあるブロック $K^{dc}\rho_cK^{dc*}$ であり、$J^{dc*}J^{dc}$ は位置 $(c,c)$ にある
$K^{dc*}K^{dc}$ である。$d,l$ について和をとると $\Phi(\rho)$ と $\Phi^*(1)$ が得られ、これで (1) の
等式が示され、すべての項は $\mathcal M$ に属する。したがって $\hat{\mathcal L}$ とそのすべてのべきは $\mathcal M$
を $\mathcal M$ に写し、$e^{t\hat{\mathcal L}}|_{\mathcal M}=e^{t\mathcal L}$ であり、これは
定理 8.6 と補題 8.4(a) により $\mathcal M$ のチャネルである。(2) $L=\sum_{d,c}J^{dc}$ に対し $L\rho L^*=\sum_{d,d',c}J^{dc}\rho J^{d'c*}$ であり、
$E$ は $d=d'$ の項を残す。$B\in\mathcal M$ に対し $P_dXBP_d=P_dXP_dB$ であるから $E(XB)=E(X)B$ かつ
$E(BX)=BE(X)$ である。よって $E\{L^*L,\rho\}=\{E(L^*L),\rho\}$ であり、$d\ne d'$ に対し $J^{dc*}J^{d'c}=0$ であるから
$E(L^*L)=\sum_cP_cL^*LP_c=\sum_{c,d}J^{dc*}J^{dc}$ である。最後に $E(-i[H,\rho])=-i[H,\rho]$ である。$\square$

> **物理的意味.** (2) は、$\mathcal M$ 上のリンドブラッド方程式が、ジャンプによって生じるセクター間の
> コヒーレンスを直ちに消去したときに $M_D$ 上の方程式から残るものであることを述べている——
> これが超選択則の内容である。粒子損失の場合には何も消去されない：消滅
> 作用素 $a_j$ は粒子数をちょうど一つ下げるので、$a_j\rho a_j^*$ と $a_j^*a_j=n_j$ は
> ブロック対角であり、$\mathcal M$ 上の方程式は通常のものである。

**命題 8.9（セクターの重み）.** $\rho(t)=e^{t\mathcal L}\rho_0$ とし、$p(t)$ をそのセクターの重みとする。

1. セクター間のジャンプがなければ、すべての $p_c(t)$ は一定である。
2. すべての $m$ とすべての組 $d\ne c$ に対し、数
   $r_{m,dc}$ により $\sum_lK^{dc*}_{m,l}K^{dc}_{m,l}=r_{m,dc}1_{k_c}$ となるならば、$d\ne c$ に対し $Q_{dc}=\sum_mg_mr_{m,dc}$、$Q_{cc}=-\sum_{d\ne c}Q_{dc}$ として $p'=Qp$ である。
   したがって $p(t)=e^{tQ}p(0)$ は $H$ にも、$L_j$ にも、セクター内の状態にもよらない。

*証明.* $p_d'=\operatorname{Tr}(P_d\mathcal L(\rho))$ である。ハミルトニアンの項は定理 8.3 と同様に消える。
$L_j\in\mathcal M$ に対しては、$P_d$ が中心に属するので $\operatorname{Tr}(P_dL_j\rho L_j^*)-\operatorname{Tr}(P_dL_j^*L_j\rho)
=\operatorname{Tr}(L_j^*P_dL_j\rho)-\operatorname{Tr}(L_j^*L_jP_d\rho)=0$ である。セクター間の
ジャンプ $\Phi$ に対しては $\operatorname{Tr}(P_d\Phi(\rho))=\sum_c\operatorname{Tr}\big(\sum_lK^{dc*}_lK^{dc}_l\rho_c\big)$
かつ $\operatorname{Tr}(P_d\tfrac12\{\Phi^*(1),\rho\})=\sum_{d'}\operatorname{Tr}\big(\sum_lK^{d'd*}_lK^{d'd}_l\rho_d\big)$ であり、
$c=d$ の項と $d'=d$ の項は打ち消し合う。これで (1) が示され、(2) の仮定のもとでは残りの項は
$\sum_{c\ne d}r_{dc}p_c-\sum_{d'\ne d}r_{d'd}p_d$ である。補題 8.1 により $p(t)=e^{tQ}p(0)$ である。$\square$

(2) の仮定のもとでは、$d\ne c$ に対し、$r_{m,dc}$ は任意の状態 $\rho$ について $\Phi_m$ の遷移
行列（命題 6.12 の前で定義した）の成分 $T_{dc}(\rho)$ に等しい。したがって、レート
$1$ の一つのジャンプに対しては、$R$ = `transition_matrix()` として $Q=R-\operatorname{diag}(\sum_dR_{dc})_c$ である（$R$ の対角成分は打ち消し合う）。

**例：粒子損失.** 粒子数が高々二の二つのボース・モードを考える。セクターは
$N=0,1,2$ で、基底はそれぞれ $|00\rangle$；$|10\rangle,|01\rangle$；$|20\rangle,|11\rangle,|02\rangle$ であるから、
$\mathcal M=M_1\oplus M_2\oplus M_3$ である。ボース–ハバード二量体
$H=J(a_1^*a_2+a_2^*a_1)+\frac U2\sum_jn_j(n_j-1)$ は $N$ と可換であり、各モードは
レート $\gamma$ で粒子を失う：これはクラウス作用素 $a_1,a_2$ を $N\to N-1$ に制限した、セクター間の一つのジャンプである。
$a_1^*a_1+a_2^*a_2=N$ であるから、命題 8.9(2) が $Q_{N-1,N}=\gamma N$ として適用でき、$N=2$ から出発すると
$$
p_2(t)=q^2,\qquad p_1(t)=2q(1-q),\qquad p_0(t)=(1-q)^2,\qquad q=e^{-\gamma t},
$$
となる。これは二項分布である：ホッピングと相互作用が何をしようとも、各粒子は独立に確率 $e^{-\gamma t}$ で生き残る
（$p'=Qp$ は微分して確かめよ）。

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

RK4 は厳密には完全正値でない：$\rho(t)$ の固有値は離散化誤差の程度で
負になりうる（ここでは刻み幅が高々 $0.02$ で $\approx-2\cdot10^{-8}$）。`lindblad_channel` は厳密であるが、
セクター内のジャンプしか扱わない。最後に、セクター間のランダムな写像について命題 8.8(2) を、
numpy で組み立てた $M_6$ 上の GKSL 生成子と照合して確かめる：

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

## 8.4 事後選択とジャンプなしの時間発展

8.3節の生成子のジャンプをまとめて $\Phi_{\mathrm{tot}}(\rho)=\sum_jg_jL_j\rho L_j^*+\sum_mg_m\Phi_m(\rho)$ とし、
$\mathcal M$ の元 $A=\Phi_{\mathrm{tot}}^*(1)=\sum_jg_jL_j^*L_j+\sum_mg_m\Phi_m^*(1)\ge0$ と
$H_{\mathrm{eff}}=H-\frac i2A$ をとる。命題 8.8(1) と定理 8.6（その $S_n$ は、$H_{\mathrm{eff}}\in\mathcal M$ であるから
$\mathcal M$ を $\mathcal M$ に写す）により、
$$
\rho(t)=\sum_{n\ge0}S_n(t)\rho,\qquad S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},
\qquad S_n(t)=\int_0^tS_0(t-s)\,\Phi_{\mathrm{tot}}\,S_{n-1}(s)\,ds .
$$
である。

**命題 8.10.** 状態 $\rho$ に対し $w_n(t)=\operatorname{Tr}S_n(t)\rho$ とおく。このとき $w_n(t)\ge0$、
$\sum_nw_n(t)=1$ であり、$w_0$ は非増加で $w_0'(t)=-\operatorname{Tr}\big(A\,S_0(t)\rho\big)$ である。
純粋状態 $\rho=|\psi\rangle\langle\psi|$ に対しては $w_0(t)=\|e^{-itH_{\mathrm{eff}}}\psi\|^2$ である。

*証明.* $S_n(t)$ は CP であるから正値であり、$\sum_nw_n=\operatorname{Tr}\rho(t)=1$ である。
$X=S_0(t)\rho\ge0$ とおくと $w_0'=\operatorname{Tr}(-iH_{\mathrm{eff}}X+iXH_{\mathrm{eff}}^*)=-\operatorname{Tr}(AX)$
であり、$\operatorname{Tr}(AX)=\operatorname{Tr}(X^{1/2}AX^{1/2})\ge0$ である。最後の主張は
$\operatorname{Tr}|\phi\rangle\langle\phi|=\|\phi\|^2$ にほかならない。$\square$

> **物理的意味.** 検出器があらゆるジャンプ——放出されたすべての光子、失われたすべての原子——を記録するとしよう。
> このとき $w_n(t)$ は時刻 $t$ までにちょうど $n$ 回のクリックが起こる確率であり、$S_n(t)\rho/w_n(t)$ は
> $n$ 回のクリックが起きたという条件のもとでの状態、$S_0(t)\rho/w_0(t)$ は **クリックがない** という条件のもとでの状態である：
> 何も起こらなかったという結果も、非エルミートな $H_{\mathrm{eff}}$ を通じて状態を変える。この計数過程による解釈 [E. B. Davies,
> *Quantum Theory of Open Systems*, 1976; M. D. Srinivas and E. B. Davies, *Optica Acta* 28, 981 (1981)]
> は量子軌跡法（モンテカルロ波動関数法）の基礎である [J. Dalibard, Y. Castin,
> K. Mølmer, *Phys. Rev. Lett.* 68, 580 (1992)]。ここではその数学的な部分だけを用いる。

損失の場合、クリックを数えることは電荷を測定することと同じである：

**命題 8.11.** セクターに $N=0,\dots,C-1$ とラベルをつけ、セクター内の
ジャンプはなく、すべてのジャンプが電荷を一つ下げ（$d=c-1$ でない限り $K^{dc}_{m,l}=0$）、
$\rho=P_{N_0}\rho P_{N_0}$ であるとする。このとき $0\le n\le N_0$ に対し $S_n(t)\rho=P_{N_0-n}\,\rho(t)\,P_{N_0-n}$ であり、
$n>N_0$ に対し $S_n(t)\rho=0$ である。特に、$\rho(t)$ の $P_{N_0}$ によるリューダース更新（「粒子が一つも失われなかった」）
はジャンプなしの状態である：
$$
\frac{P_{N_0}\rho(t)P_{N_0}}{p_{N_0}(t)}=\frac{e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*}}{w_0(t)},
\qquad p_{N_0}(t)=w_0(t).
$$

*証明.* $H_{\mathrm{eff}}\in\mathcal M$ であるから、$S_0(s)$ はセクター $N$ に台をもつ作用素をセクター $N$ に台をもつ作用素に写す。
$\Phi_{\mathrm{tot}}$ はセクター $N\ge1$ をセクター $N-1$ に写し、セクター $0$ を零にする。
帰納法により、すべての $s$ に対し $S_n(s)\rho$ はセクター $N_0-n$ に台をもち、$n>N_0$ ならば $0$ である。
$\rho(t)=\sum_nS_n(t)\rho$ を $P_{N_0-n}$ で圧縮すると、ちょうど第 $n$ 項が残る。$\square$

ライブラリの `states.lueders_update(rho, P)` は $P\rho P/\operatorname{Tr}(P\rho P)$ と
$\operatorname{Tr}(\rho P)$ を返す。今度は二つのモードが異なるレート $\gamma_1>\gamma_2$ で漏れ出すとする。すると
$A=\gamma_1n_1+\gamma_2n_2$ は中心に属さず、命題 8.9(2) は適用できず、
$N=2$ への事後選択は真に非エルミートな時間発展を与える。命題 8.11 を $n=0$ について確かめ、さらに
$S_1$ の積分をシンプソン則で計算して $n=1$ についても確かめる：

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

$t=2$ において損失のない確率は $0.230$ である。ホッピングがなければ $e^{-2\gamma_1t}=0.091$ となるはずであるから、
ここではセクターの重みがハミルトニアンに依存している。事後選択された状態は、閉じた系の状態から
トレース距離 $0.135$ だけ離れている：粒子を失わなかったことは情報である。

## 8.5 定常状態と緩和

状態 $\sigma$ は、すべての $t\ge0$ に対し $e^{t\mathcal L}\sigma=\sigma$ となるとき **定常** であるという。これは
$\mathcal L(\sigma)=0$ と同値である（$t=0$ で微分すればよく、逆は級数を用いる）。

**命題 8.12.** $\mathcal M$ 上の任意の量子動力学的半群は定常状態をもつ。

*証明.* 状態 $\rho$ を固定し、$\sigma_\tau=\frac1\tau\int_0^\tau T_s\rho\,ds$ とおく。これは状態である（状態の
凸結合の極限であり、状態全体は閉集合をなす）。線形性により
$T_t\sigma_\tau-\sigma_\tau=\frac1\tau\big(\int_\tau^{\tau+t}-\int_0^t\big)T_s\rho\,ds$ であり、$\|T_s\rho\|_1=1$ であるから
そのトレースノルムは高々 $2t/\tau$ である。状態全体は有限次元空間の有界閉集合であるから、
ある列 $\sigma_{\tau_m}$（$\tau_m\to\infty$）は状態 $\sigma$ に収束し
（ボルツァーノ–ワイエルシュトラス）、$T_t\sigma-\sigma=\lim_m(T_t\sigma_{\tau_m}-\sigma_{\tau_m})=0$ である。$\square$

**系 8.13.** セクター間のジャンプがなければ、定常状態の張る空間の
次元は少なくとも $C$ である：各セクターがそれ自身の定常状態をもつ。

*証明.* $\mathcal L$ は $P_c\mathcal MP_c\cong M_{k_c}$ をそれ自身に写す（$H$、$L_j$ はすべてブロック対角である）から、
$e^{t\mathcal L}$ の制限は $M_{k_c}$ 上の量子動力学的半群であり（補題 8.4(a)）、命題 8.12 により
定常状態 $\sigma_c$ をもつ。$\sigma_c$ たちの台は互いに交わらない。$\square$

**命題 8.14（自然放出）.** 基底 $|g\rangle,|e\rangle$ をもつ $M_2$ 上で、
$H=\omega|e\rangle\langle e|$、一つのジャンプ $L=|g\rangle\langle e|$、レート $\gamma>0$ とする。このとき
$$
\rho_{ee}(t)=e^{-\gamma t}\rho_{ee}(0),\qquad\rho_{ge}(t)=e^{(i\omega-\gamma/2)t}\rho_{ge}(0),
\qquad e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}},
$$
である。ここで $\mathcal U_t(\rho)=e^{-itH}\rho e^{itH}$ であり、$\mathcal A_p$ はパラメータ
$p$ の振幅減衰チャネル（6.4節）である。$\mathcal L$ のスペクトルは $\{0,-\gamma,-\gamma/2\pm i\omega\}$ であり、
$|g\rangle\langle g|$ がただ一つの定常状態であって、
$\tfrac12\|\rho(t)-|g\rangle\langle g|\|_1=\big(\rho_{ee}(0)^2e^{-2\gamma t}+|\rho_{ge}(0)|^2e^{-\gamma t}\big)^{1/2}$ である。

*証明.* ここで $L^*L=|e\rangle\langle e|$、$L\rho L^*=\rho_{ee}|g\rangle\langle g|$ であるから、成分ごとに
$\rho_{ee}'=-\gamma\rho_{ee}$、$\rho_{gg}'=\gamma\rho_{ee}$、$\rho_{ge}'=(i\omega-\gamma/2)\rho_{ge}$ であり、
これから解が得られる。基底 $(\rho_{gg},\rho_{ge},\rho_{eg},\rho_{ee})$ において $\mathcal L$ の行列は
上三角で、対角成分は $(0,-\gamma/2+i\omega,-\gamma/2-i\omega,-\gamma)$ である。$\mathcal L(\sigma)=0$ からは
$\sigma_{ee}=\sigma_{ge}=0$ が従う。振幅減衰チャネルは $\rho_{ee}$ に $1-p$ を、$\rho_{ge}$ に
$\sqrt{1-p}$ を掛け（そのクラウス作用素は $|g\rangle\langle g|+\sqrt{1-p}\,|e\rangle\langle e|$ と
$\sqrt p\,|g\rangle\langle e|$）、$\mathcal U_t$ は $\rho_{ge}$ に $e^{i\omega t}$ を掛ける。どちらも TP であり、これで
$\rho_{gg}$ も決まる。最後に、$\rho(t)-|g\rangle\langle g|$ はトレース零の自己共役作用素で、対角成分は
$(-\rho_{ee}(t),\rho_{ee}(t))$ であるから、その固有値は $\pm(\rho_{ee}(t)^2+|\rho_{ge}(t)|^2)^{1/2}$ である。$\square$

コヒーレンスは占有数の半分のレートで減衰する。これは純粋位相緩和がない場合の 6.4節の関係
$T_2=2T_1$ である。
最も遅い非零のレート、すなわち **スペクトルギャップ** $\gamma/2$ が緩和を支配する。スペクトルを
`lindblad_superoperator` で、チャネルとトレース距離を `lindblad_channel` で確かめ、さらに
$M_2\oplus M_3$ のセクター内のランダムなジャンプについて系 8.13 を確かめる：

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

これらのランダムなジャンプについて、コードはセクターごとにちょうど一つの定常状態を見つける。
一意性の判定条件——たとえば、忠実な定常状態が存在し、かつ $\{H,L_j,L_j^*\}'=\mathbb C1$ であること
[A. Frigerio, *Commun. Math. Phys.* 63, 269 (1978)]——はここでは証明しない。演習問題 4 は、
因子環においてさえ一意性が成り立たないことがあることを示す。

## 8.6 比較の対象

QuTiP（`sesolve`、適応的 ODE ソルバーによるリンドブラッド方程式用の `mesolve`、量子軌跡用の `mcsolve`、
`steadystate`、`liouvillian`）、dynamiqs（JAX。GPU 上でバッチ化され微分可能な
シュレーディンガー・ソルバーとリンドブラッド・ソルバー）、QuantumOptics.jl、および SciPy（`scipy.linalg.expm`、および
疎なベクトル化された生成子に対して $e^{t\mathcal L}\rho$ を計算する `scipy.sparse.linalg.expm_multiply`）は、
同じ方程式を一つのヒルベルト空間上で解く。上の問題に対しては、それらは同じ数値を与えなければならない。
$\mathcal M$ 上では、それらに $D\times D$ のブロック対角行列を渡すことになる。このときベクトル化された生成子は
$D^2$ 次元の空間に作用し（疎ではあるが、一回の求解の中でセクターどうしが結合される）、ユーザーが
問題を手作業でセクターに分割しない限りそうなる。ライブラリはブロック構造を組み込んでいる：セクターごとの
厳密なプロパゲーターとリンドブラッド・チャネル（$k_c^2\times k_c^2$ の超演算子）、命題 8.8 の圧縮を伴う
`InterSectorChannel` としてのセクター間のジャンプ、ハミルトニアンとジャンプ
作用素のバッチ（不規則性の平均のため）、そして自動微分（最後の二つは dynamiqs にもある）。
適応的な刻み幅制御はなく（RK4 は固定刻みで、近似的にしか正値性を保たない）、軌跡
ソルバーもなく、`lindblad_channel` はセクター間のジャンプを受け付けない。それらについては、$M_D$ 上の QuTiP や dynamiqs
が基準となる。

## 演習問題

1. （ハイゼンベルク描像）8.3節の生成子について、
   $\mathcal L^*(A)=i[H,A]+\sum_jg_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})+\sum_mg_m(\Phi_m^*(A)-\frac12\{\Phi_m^*(1),A\})$
   として $\operatorname{Tr}(\mathcal L(\rho)A)=\operatorname{Tr}(\rho\,\mathcal L^*(A))$ であること、
   $\mathcal L^*(1)=0$ であること、およびセクター間のジャンプがなければ $\mathcal L^*(P_c)=0$ であることを示せ。
2. （例外点）8.4節のセクター $N=1$ では
   $H_{\mathrm{eff}}=\begin{pmatrix}-i\gamma_1/2&J\\J&-i\gamma_2/2\end{pmatrix}$ である。その固有値を求め、
   $J=|\gamma_1-\gamma_2|/4>0$ のとき対角化可能でないことを示せ。
3. $m$ 個のモードにある $N_0$ 個のボース粒子がそれぞれレート $\gamma$ で失われるとする（ジャンプ $a_1,\dots,a_m$ を
   $N\to N-1$ に制限したもの）。$\mathcal M$ に属する任意のハミルトニアンに対し、$q=e^{-\gamma t}$ として
   $p_N(t)=\binom{N_0}Nq^N(1-q)^{N_0-N}$ であることを示せ。
4. （純粋位相緩和）$M_2$ 上で $H=\omega|e\rangle\langle e|$、$L=\sigma_z=\operatorname{diag}(1,-1)$、レート
   $\gamma$ とする。占有数は一定であり、$\rho_{ge}(t)=e^{(i\omega-2\gamma)t}\rho_{ge}(0)$ であること、そして
   $C=1$ であるにもかかわらず定常状態が対角行列からなる二次元の空間を張ることを示せ。

<details>
<summary>解答</summary>

1. 巡回性により $\operatorname{Tr}(L\rho L^*A)=\operatorname{Tr}(\rho L^*AL)$、
   $\operatorname{Tr}(\{L^*L,\rho\}A)=\operatorname{Tr}(\rho\{L^*L,A\})$、
   $\operatorname{Tr}(-i[H,\rho]A)=\operatorname{Tr}(\rho\,i[H,A])$ である。セクター間のジャンプについては、
   $\operatorname{Tr}(\Phi(\rho)A)=\operatorname{Tr}(\rho\Phi^*(A))$ が $\Phi^*$ の定義である
   （6.3節）。すると $\mathcal L^*(1)=\sum_jg_j(L_j^*L_j-L_j^*L_j)+\sum_mg_m(\Phi_m^*(1)-\Phi_m^*(1))=0$ であり、
   これは再びトレース保存にほかならない。セクター間のジャンプがなければ、$P_c$ は $H$ および
   すべての $L_j$ と可換であるから、$L_j^*P_cL_j=L_j^*L_jP_c=\frac12\{L_j^*L_j,P_c\}$ となり $\mathcal L^*(P_c)=0$ である。したがって
   $p_c'=\operatorname{Tr}(\rho\mathcal L^*(P_c))=0$ であり、これは命題 8.9(1) の第二の証明である。

2. $\delta=(\gamma_1-\gamma_2)/4$ として
   $H_{\mathrm{eff}}=-\frac i4(\gamma_1+\gamma_2)1+\begin{pmatrix}-i\delta&J\\J&i\delta\end{pmatrix}$ である。第二の行列はトレース $0$、
   行列式 $\delta^2-J^2$ をもつので、固有値は $-\frac i4(\gamma_1+\gamma_2)\pm\sqrt{J^2-\delta^2}$ である：
   $J<|\delta|$ ならば二つの減衰レート、$J>|\delta|$ ならば一つのレートと二つの振動数である。$J=|\delta|>0$ のとき、第二の行列 $Y$ は
   $Y^2=(J^2-\delta^2)1=0$ を満たすが $Y\ne0$ であるから、$H_{\mathrm{eff}}$ はただ一つの固有値
   $-\frac i4(\gamma_1+\gamma_2)$ をもつが $1$ の定数倍ではない：したがって対角化可能でない。

```python
lam = np.linalg.eigvals(H_eff[1])
delta = (g1 - g2) / 4
root = np.sqrt(J ** 2 - delta ** 2 + 0j)
assert np.allclose(sorted(lam, key=np.real), sorted(-0.25j * (g1 + g2) + np.array([root, -root]), key=np.real))
H_ep = np.array([[-0.5j * g1, delta], [delta, -0.5j * g2]])         # J = delta: the exceptional point
Y = H_ep + 0.25j * (g1 + g2) * np.eye(2)                           # H_ep minus its only eigenvalue
assert np.abs(Y @ Y).max() < 1e-15 and np.linalg.matrix_rank(Y) == 1   # one Jordan block of size 2
```

3. セクター $N$ 上で $\sum_ja_j^*a_j=N1$ であり、セクターの組は $(N-1,N)$ だけであるから、命題
   8.9(2) により $p_N'=-\gamma Np_N+\gamma(N+1)p_{N+1}$ である。$q'=-\gamma q$ を用いると、二項分布の
   項の導関数は $-\gamma Np_N+\gamma(N_0-N)\binom{N_0}Nq^{N+1}(1-q)^{N_0-N-1}$ であり、
   $(N+1)\binom{N_0}{N+1}=(N_0-N)\binom{N_0}N$ により第二項は $\gamma(N+1)p_{N+1}$ である。$t=0$ では
   これは $N_0$ に集中した点質量であり、一意性は補題 8.1 による。

4. $L\rho L^*-\frac12\{L^*L,\rho\}=\sigma_z\rho\sigma_z-\rho$ は非対角成分に $-2$ を掛け、
   対角成分を消す。また $-i[H,\rho]_{ge}=i\omega\rho_{ge}$ である。したがって
   （$\gamma>0$ のとき）$\mathcal L(\sigma)=0$ であるための必要十分条件は $\sigma_{ge}=0$ であり、スペクトルは $\{0,0,-2\gamma\pm i\omega\}$ である。$H$ と $\sigma_z$ は可換で、
   $X\mapsto\sigma_zX\sigma_z$ は対合であるから、$U_t=e^{-itH}$ として $e^{t\mathcal L}(\rho)=\frac{1+e^{-2\gamma t}}2U_t\rho U_t^*
   +\frac{1-e^{-2\gamma t}}2U_t\sigma_z\rho\sigma_zU_t^*$ である：これは単位的チャネルであり、その不動点は
   クラウス作用素と可換な行列（定理 6.10）、すなわち対角行列である。

```python
deph = dynamics.lindblad_superoperator(q2, H_at, [q2.from_blocks([torch.diag(torch.tensor([1.0, -1.0], dtype=f64))])],
                                       rates=[gamma])[0][0]
ev = torch.linalg.eigvals(deph)
assert int((ev.abs() < 1e-12).sum()) == 2
assert np.allclose(sorted(ev.tolist(), key=lambda z: (z.real, z.imag))[:2],
                   [-2 * gamma - 1j * omega, -2 * gamma + 1j * omega], atol=1e-12)
```

</details>

**次へ：**第9章（準備中）ではスピン鎖を扱う：磁化と運動量のセクター、
大きなセクターに対するランチョス法とクリロフ法による時間発展、そして不規則性。
