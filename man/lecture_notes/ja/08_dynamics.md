[English](../08_dynamics.md) · [Español](../es/08_dynamics.md) · [Français](../fr/08_dynamics.md) · [Deutsch](../de/08_dynamics.md) · [中文](../zh/08_dynamics.md) · **日本語** · [Русский](../ru/08_dynamics.md)

# 8. ダイナミクス

量子力学では、閉じた系はユニタリの一径数群によって、開いた系——観測されない環境と
結合した系——はチャネルの半群によってモデル化される。セクターをもつ環では、どちらの発展にも
単一の行列環の上では見えない特徴がある：$`\mathcal M`$ に属するハミルトニアンは決してセクターを変えられず、
したがって電荷は保存されるが、開いた系は、たとえば粒子を失うことによって、セクター間をジャンプしうる。
本章では両種の時間発展についての基本的事実を証明し、セクター間のジャンプを含むリンドブラッド方程式を書き下し、
「粒子が一つも失われなかった」ことが電荷の測定であって、その事後選択された状態が非エルミートな
ハミルトニアンによって発展することを示す。

**本章で学ぶこと**：セクターの重みが運動の定数である理由、
ゴリーニ–コサコフスキー–スダルシャン–リンドブラッド（GKSL）形の生成子がチャネルを生み出す理由、セクター間のジャンプが
リンドブラッド方程式にどのように入り、セクターの重みがいつ古典的なレート方程式に従うか、ジャンプなしの
時間発展と事後選択が何を意味するか、そして定常状態と緩和レートの求め方。本章の数値的検証は
[チュートリアル 8](../../tutorial/ja/08_dynamics.md) にある。

本章を通じて、$`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$ は
$`\mathcal H=\bigoplus_c\mathbb C^{k_c}`$ に作用し、$`D=\sum_ck_c`$、$`P_c`$ は極小中心射影とする。
状態はその blunt 密度 $`\rho`$ で与え、セクターの重みを $`p_c=\mathrm{Tr}\,\rho_c`$ とする
（第3章）。また $`[A,B]=AB-BA`$、$`\{A,B\}=AB+BA`$、$`\hbar=1`$ とする。**超演算子** とは $`\mathcal M`$ または
$`M_D(\mathbb C)`$ 上の線形写像のことである。超演算子のノルムはすべて作用素ノルムとし、したがって
劣乗法的である。有限次元空間上の線形写像 $`B`$ に対し $`e^{tB}=\sum_{n\ge0}t^nB^n/n!`$ とする。

<a id="s-8-1"></a>

## 8.1 閉じた系

**ハミルトニアン** とは自己共役な $`H\in\mathcal M`$ のことである。ベクトルは **シュレーディンガー
方程式** $`i\psi'(t)=H\psi(t)`$ に従って、密度は **フォン・ノイマン方程式** $`\rho'(t)=-i[H,\rho(t)]`$ に従って発展する。
どちらも有限次元空間上の定数係数線形方程式であり、
次の補題が我々に必要な微分方程式論のすべてである。

<a id="r-8-1"></a>

**補題 8.1.** $`B`$ を有限次元ノルム空間 $`X`$ 上の線形写像とする。このとき
$`\frac d{dt}e^{tB}=Be^{tB}=e^{tB}B`$、$`e^{(s+t)B}=e^{sB}e^{tB}`$、$`\|e^{tB}\|\le e^{|t|\,\|B\|}`$ であり、
任意の $`x_0\in X`$ に対し、写像 $`x(t)=e^{tB}x_0`$ は、区間
$`I\ni0`$ 上で $`x'=Bx`$ かつ $`x(0)=x_0`$ を満たす微分可能な $`x`$ としてただ一つのものである。

*証明.* 級数は上のノルム評価とともに絶対収束し、$`t`$ のべき級数として
項別微分できるので、第一の公式が得られる。すると
$`\frac d{dt}(e^{tB}e^{-tB})=Be^{tB}e^{-tB}-e^{tB}Be^{-tB}=0`$ であるから $`e^{tB}e^{-tB}=1`$ である。$`x'=Bx`$ ならば
$`\frac d{dt}(e^{-tB}x(t))=e^{-tB}(x'-Bx)=0`$ であるから、$`e^{-tB}x(t)=x_0`$、すなわち $`x(t)=e^{tB}x_0`$ である。最後に
$`t\mapsto e^{(s+t)B}x_0`$ と $`t\mapsto e^{tB}e^{sB}x_0`$ はともに同じ初期値をもつ $`y'=By`$ の解である。$`\square`$

<a id="r-8-2"></a>

**命題 8.2.** $`H=H^*\in\mathcal M`$ とし、$`U(t)=e^{-itH}`$ を 1.3節の関数解析により
$`f_t(\lambda)=e^{-it\lambda}`$ として定める。

1. $`U(t)=\bigoplus_ce^{-itH_c}`$ は $`\mathcal M`$ のユニタリであり、$`U(s+t)=U(s)U(t)`$ かつ $`U'(t)=-iHU(t)`$ である。
   $`U(t)`$ は $`-itH`$ の指数級数に等しい。
2. $`\psi(t)=U(t)\psi_0`$ は $`\psi(0)=\psi_0`$ を満たす $`i\psi'=H\psi`$ のただ一つの解であり、
   $`\rho(t)=U(t)\rho_0U(t)^*`$ は
   $`\rho(0)=\rho_0`$ を満たす $`\rho'=-i[H,\rho]`$ の $`\mathcal M`$ におけるただ一つの解である。$`\rho_0`$ が状態ならば、すべての $`\rho(t)`$ も状態である。
3. （*ハイゼンベルク描像*）
   $`A(t)=U(t)^*AU(t)`$ として $`\mathrm{Tr}(\rho(t)A)=\mathrm{Tr}(\rho_0A(t))`$ であり、$`A'=i[H,A]`$ である。

*証明.* (1) 各ブロックを対角化する：$`H_c=V_c\Lambda_cV_c^*`$。すると $`V=\bigoplus_cV_c`$ は
$`H`$ を対角化し、$`U(t)=Ve^{-it\Lambda}V^*`$ はブロック対角である（定理 1.5(a) により、結果は
対角化に用いるユニタリによらない）。$`|f_t|^2=1`$ かつ $`f_sf_t=f_{s+t}`$ であるから、定理 1.5(c) によりユニタリ性と
群法則が得られる。$`e^{-it\Lambda}`$ を成分ごとに微分すると $`U'=V(-i\Lambda)e^{-it\Lambda}V^*=-iHU`$ であるから、
補題 8.1（$`X=M_D`$、$`B`$ は $`-iH`$ による左からの乗法）により $`U(t)`$ は級数 $`e^{-itH}`$ である。
(2) $`(U\psi_0)'=-iHU\psi_0`$ であり、$`(U')^*=iU^*H`$ であるから $`(U\rho_0U^*)'=-iHU\rho_0U^*+U\rho_0U^*(iH)`$ である。
一意性は補題 8.1 を、$`\mathcal H`$ 上の $`B=-iH`$、あるいは $`\mathcal M`$ 上の $`B=-i[H,\cdot\,]`$ に適用したものである。正値性：
$`\langle x,U\rho_0U^*x\rangle=\langle U^*x,\rho_0U^*x\rangle\ge0`$。トレースは巡回性による。(3) トレースの
巡回性と、$`(U^*AU)'=iU^*HAU-iU^*AHU=i[H,A(t)]`$ による。$`\square`$

*ライブラリでは：* ブロックごとの $`U(t)=e^{-itH}`$ は `dynamics.propagator` であり、命題 8.2(2) の解は `dynamics.schrodinger` と `dynamics.von_neumann` である — [チュートリアル、8.1節](../../tutorial/ja/08_dynamics.md#s-8-1)。

<a id="r-8-3"></a>

**定理 8.3（電荷の保存）.** $`t\mapsto H(t)`$ を区間 $`I`$ 上で $`\mathcal M`$ の自己共役元に値をとる
写像とする。

1. $`\psi:I\to\mathcal H`$ が微分可能で $`i\psi'=H(t)\psi`$ を満たすならば、すべての $`\|P_c\psi(t)\|`$ は一定である。
2. $`\rho:I\to\mathcal M`$ が微分可能で $`\rho'=-i[H(t),\rho]`$ を満たすならば、すべての中心元 $`Z`$ に対し $`\mathrm{Tr}(Z\rho(t))`$ は
   一定である。特に、すべてのセクターの重み $`p_c(t)=\mathrm{Tr}(P_c\rho(t))`$ は一定である。
3. $`H`$ が定数で $`A\in\mathcal M`$ が $`H`$ と可換ならば、$`\mathrm{Tr}(A\rho(t))`$ は一定である。

*証明.* (1) 内積は第一引数について反線形であるから、
$`\frac d{dt}\langle\psi,P_c\psi\rangle=\langle-iH\psi,P_c\psi\rangle+\langle\psi,-iP_cH\psi\rangle =i\langle\psi,HP_c\psi\rangle-i\langle\psi,P_cH\psi\rangle=0`$ である。$`P_c`$ は中心に属するからである。(2)
巡回性により $`\frac d{dt}\mathrm{Tr}(Z\rho)=-i\mathrm{Tr}(Z[H,\rho])=-i\mathrm{Tr}([Z,H]\rho)=0`$ である。
(3) $`A`$ について同じ計算をすればよい。$`\square`$

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. $`\mathcal M`$ に属するハミルトニアンは電荷（粒子数、
> 磁化、対称性の量子数）と可換であるから、各電荷の確率は運動の
> 定数である。これは駆動された時間依存のハミルトニアンに対しても成り立つ。

$`\mathcal M`$ 上では、単位ベクトル $`\psi=\bigoplus_c\psi_c`$ は混合
$`\bigoplus_c|\psi_c\rangle\langle\psi_c|`$ の状態を定める：すべての $`A\in\mathcal M`$ はブロック対角であるから
$`\langle\psi,A\psi\rangle=\sum_c\langle\psi_c,A_c\psi_c\rangle=\mathrm{Tr}\big(A\bigoplus_c|\psi_c\rangle\langle\psi_c|\big)`$
である（第3章の演習問題 3 と同様）。$`\psi(t)`$ が $`i\psi'=H\psi`$ の解ならば、$`U(t)`$ はブロック対角であるから
$`\bigoplus_c|\psi_c(t)\rangle\langle\psi_c(t)|=U(t)\big(\bigoplus_c|\psi_c(0)\rangle\langle\psi_c(0)|\big)U(t)^*`$ である。

*計算：* 定理 8.3、および $`\psi`$ と $`\bigoplus_c|\psi_c\rangle\langle\psi_c|`$ の時間発展の一致を、厳密なソルバーと、駆動されたハミルトニアンに対するルンゲ–クッタ・ソルバー `dynamics.schrodinger_rk4` で確かめる — [チュートリアル、8.1節](../../tutorial/ja/08_dynamics.md#s-8-1)。

<a id="s-8-2"></a>

## 8.2 量子動力学的半群

CP 写像、TP 写像およびチャネルは 6.1節（[第6章](06_quantum_channels.md)）で定義した。
定理 6.2 と命題 6.4（写像 $`A\mapsto\sum_lK_lAK_l^*`$ は、セクター間のものも含めて CP である）のほかに、
三つの閉性が必要である。

<a id="r-8-4"></a>

**補題 8.4.** (a) CP 写像の和、非負の定数倍、合成は CP である。$`\Phi`$ を
$`\mathcal A=M_D`$ または $`\mathcal A=\mathcal M`$ 上の CP 写像とし、$`\mathcal N\subseteq\mathcal A`$ を（$`\mathcal A=M_D`$
のときは）$`\mathcal M`$、あるいはコーナー $`P_c\mathcal MP_c\cong M_{k_c}`$（単位元 $`P_c`$ をもつ）の一つとする。$`\Phi`$ が
$`\mathcal N`$ をそれ自身に写すならば、$`\Phi|_{\mathcal N}`$ は CP である。
(b) CP 写像 $`\Phi_m`$ が $`\Phi`$ に収束する（すなわちすべての $`A`$ に対し $`\Phi_m(A)\to\Phi(A)`$）ならば、$`\Phi`$ は CP である。
(c) $`s\mapsto\Psi(s)`$ が $`[a,b]`$ から CP 写像への連続写像ならば、$`\int_a^b\Psi(s)\,ds`$ は CP である。

*証明.* (a) $`\Phi\mapsto\Phi_n`$ は線形であり、$`(\Phi\Psi)_n=\Phi_n\Psi_n`$ である。正の行列の和は
正である。制限について：$`M_n(\mathcal N)\subseteq M_n(\mathcal A)`$ であり、正値性は両者で同じ
条件である。実際、6.1節により $`M_n(\mathcal M)`$ の元 $`[A_{ij}]`$ が正であるのは、すべての
ブロック行列 $`[(A_{ij})_c]`$ が正であるときかつそのときに限る。$`\mathcal A=M_D`$ に対しては、これは $`[A_{ij}]\in M_{nD}`$ の
正値性である。基底を並べ替えると $`[A_{ij}]`$ はこれらのブロックの直和になるからである（すべての $`A_{ij}`$ はブロック対角である）。
また $`M_n(P_c\mathcal MP_c)`$ の元はただ一つの零でないブロック $`[(A_{ij})_c]`$ をもち、その正値性は
$`M_n(M_{k_c})`$ における正値性である。(b) 正の行列全体は閉集合をなし（各
$`x`$ に対する $`\langle x,Xx\rangle\ge0`$ は閉条件である）、$`(\Phi_m)_n(X)\to\Phi_n(X)`$ である。(c) リーマン和は (a) により CP であり、
積分に収束する（補題 A.4(a)）。(b) を適用すればよい。$`\square`$

**定義.** $`\mathcal M`$ 上の **量子動力学的半群** とは、
$`T_0=\mathrm{id}`$、$`T_{s+t}=T_sT_t`$ を満たし、$`t\mapsto T_t`$ が連続であるようなチャネル $`\mathcal M\to\mathcal M`$ の族 $`(T_t)_{t\ge0}`$ のことである。

<a id="r-8-5"></a>

**補題 8.5.** 超演算子の任意の連続な半群 $`(T_t)_{t\ge0}`$（$`T_0=\mathrm{id}`$、$`T_{s+t}=T_sT_t`$）は、ただ一つの超演算子 $`\mathcal L`$ により
$`T_t=e^{t\mathcal L}`$ と書ける。$`\mathcal L`$ をその **生成子** といい、$`\mathcal L=\lim_{t\to0}(T_t-\mathrm{id})/t`$ である。
任意の有限次元空間上の線形写像（たとえば行列）からなる連続な半群についても、同じ証明で同じことが成り立つ。

*証明.*（超演算子値の連続関数の積分は補題 A.4 のものである。）
$`V_h=\frac1h\int_0^hT_s\,ds`$ に対し $`\|V_h-\mathrm{id}\|\le\max_{0\le s\le h}\|T_s-\mathrm{id}\|\to0`$ であるから（補題 A.4(c)）、
行列式の連続性により、ある固定した小さな $`h`$ に対し $`\det V_h\neq0`$ である。ここで
$`T_tV_h=\frac1h\int_0^hT_{t+s}\,ds=\frac1h\int_t^{t+h}T_s\,ds`$（補題 A.4(b),(e)）は $`t`$ について微分可能で、導関数は
$`\frac1h(T_{t+h}-T_t)=T_t\frac{T_h-\mathrm{id}}h`$ である。したがって $`T_t`$ は微分可能で、
$`\mathcal L=\frac{T_h-\mathrm{id}}hV_h^{-1}`$ として $`T_t'=T_t\mathcal L`$ であり、補題 8.1（$`B`$ は
$`\mathcal L`$ を右から合成する写像）により $`T_t=e^{t\mathcal L}`$ である。一意性：$`\mathcal L=T'(0)`$。$`\square`$

<a id="r-8-6"></a>

**定理 8.6（GKSL 生成子）.** $`H=H^*`$ および $`L_1,\dots,L_J`$ を $`M_D`$ の元とし、
```math
\mathcal L(\rho)=-i[H,\rho]+\sum_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
```
とする。このとき $`(e^{t\mathcal L})_{t\ge0}`$ は $`M_D`$ 上の量子動力学的半群である。より詳しく、
$`H_{\mathrm{eff}}=H-\frac i2\sum_jL_j^*L_j`$、$`\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)`$、
$`\Phi(\rho)=\sum_jL_j\rho L_j^*`$ とおいて $`\mathcal L=\mathcal L_0+\Phi`$ とし、
```math
S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},\qquad
S_n(t)=\int_0^tS_0(t-s)\,\Phi\,S_{n-1}(s)\,ds\quad(n\ge1).
```
と定める。このときすべての $`S_n(t)`$ は CP であり、$`e^{t\mathcal L}=\sum_{n\ge0}S_n(t)`$（**ダイソン展開**）であって、
$`\|S_n(t)\|\le e^{t\|\mathcal L_0\|}(t\|\Phi\|)^n/n!`$ が成り立つ。

*証明.* *ステップ 1.* $`X(t)=e^{-itH_{\mathrm{eff}}}\rho e^{itH_{\mathrm{eff}}^*}`$ は $`X'=\mathcal L_0(X)`$
と $`X(0)=\rho`$ を満たすので、補題 8.1 により $`S_0(t)=e^{t\mathcal L_0}`$ である。これはただ一つのクラウス作用素
$`e^{-itH_{\mathrm{eff}}}`$ をもつ CP 写像である。*ステップ 2.* $`a=\|\mathcal L_0\|`$、$`b=\|\Phi\|`$ とおく。帰納法により $`S_n`$
は連続で、
$`\|S_n(t)\|\le\int_0^te^{(t-s)a}\,b\,e^{sa}\frac{(bs)^{n-1}}{(n-1)!}ds=e^{ta}\frac{(bt)^n}{n!}`$ であるから、
級数は任意の $`[0,T]`$ 上で一様収束する。*ステップ 3.*
$`S_n(t)=S_0(t)\int_0^te^{-s\mathcal L_0}\Phi S_{n-1}(s)\,ds`$ であるから、積の微分法則と微積分学の基本
定理（補題 A.4(e)）により $`S_n'=\mathcal L_0S_n+\Phi S_{n-1}`$（$`n\ge1`$）および $`S_0'=\mathcal L_0S_0`$ を得る。上の
評価により導関数の級数は $`[0,T]`$ 上で一様収束するので、$`Y(t)=\sum_nS_n(t)`$ は
微分可能で $`Y'=\mathcal L_0Y+\Phi Y=\mathcal LY`$、$`Y(0)=\mathrm{id}`$ であり、補題 8.1 により $`Y(t)=e^{t\mathcal L}`$ である。
*ステップ 4.* $`\Phi`$ は CP であるから（定理 6.2）、帰納法により被積分関数
$`s\mapsto S_0(t-s)\Phi S_{n-1}(s)`$ は CP 写像の連続な族であり、$`S_n(t)`$ は CP である（補題 8.4）。
部分和とその極限 $`e^{t\mathcal L}`$ も同様である。*ステップ 5.* 巡回性により
$`\mathrm{Tr}\,\mathcal L(\rho)=\sum_j\big(\mathrm{Tr}(L_j\rho L_j^*)-\mathrm{Tr}(L_j^*L_j\rho)\big)=0`$ であり、
したがって $`\frac d{dt}\mathrm{Tr}\,e^{t\mathcal L}\rho=\mathrm{Tr}\,\mathcal L(e^{t\mathcal L}\rho)=0`$ である：すべての
$`e^{t\mathcal L}`$ は TP である。半群法則と連続性は補題 8.1 による。$`\square`$

*計算：* $`M_3`$ 上の任意のジャンプについての定理 8.6。ジャンプ部分 $`\sum_{n\ge1}S_n(t)`$ の完全正値性を含む — [チュートリアル、8.2節](../../tutorial/ja/08_dynamics.md#s-8-2)。*ライブラリでは：* $`e^{t\mathcal L}`$ は `dynamics.lindblad_channel`、その行列 $`\mathcal L`$ は `dynamics.lindblad_superoperator` である — [チュートリアル、8.2節](../../tutorial/ja/08_dynamics.md#s-8-2)。

<a id="r-8-7"></a>

**定理 8.7（ゴリーニ–コサコフスキー–スダルシャン；リンドブラッド）.** (1) 逆に、$`M_D`$ 上の任意の量子
動力学的半群の生成子は、$`J\le D^2-1`$ 個のジャンプをもつ定理 8.6 の形をもつ。
$`\mathrm{Tr}\,H=0`$、$`\mathrm{Tr}\,L_j=0`$、および $`j\neq k`$ に対し $`\mathrm{Tr}(L_j^*L_k)=0`$ となるように選ぶことができる。
これらの規格化のもとで、$`H`$ と $`\Phi=\sum_jL_j\cdot L_j^*`$ は $`\mathcal L`$ によって定まり、$`L_j`$ はユニタリ $`u`$ による
$`L'_j=\sum_ku_{jk}L_k`$ を除いて定まる（短い方の族を零で補った後で）。

(2) $`\mathcal M=\bigoplus_cM_{k_c}`$ 上の任意の量子動力学的半群の生成子は
```math
\mathcal L(\rho)=\Phi(\rho)-\tfrac12\{\Phi^*(1),\rho\}-i[H,\rho]
```
である。ここで $`\Phi:\mathcal M\to\mathcal M`$ は CP 写像（セクター間のクラウス作用素、命題 6.4）、
$`H=H^*\in\mathcal M`$ である：これは 8.3節の形であり、セクター間のジャンプは一つ $`\Phi_1=\Phi`$、レートは $`1`$ である。

*証明.* (1) は定理 D.4 に定理 D.5(b),(c) の一意性の主張を合わせたものであり、(2) は定理 D.7 である。
どちらも[付録 D](D_gksl.md)で証明される。$`\square`$

鍵となるステップ（定理 D.3）は、CP 写像の半群の生成子のチョイ行列が最大エンタングルベクトルの直交補空間上で
正であることである。正値写像 $`T_t`$ についての半群の性質だけでは十分ではない（定理 D.3 の後の例）。
(2) については、条件付き期待値 $`E(X)=\sum_cP_cXP_c`$ との合成 $`T_t\circ E`$ が $`M_D`$ 上の CP 写像である。$`t=0`$ における
その導関数 $`\mathcal L\circ E`$ のチョイ行列はベクトル $`(1\otimes P_c)\Omega`$ たちの直交補空間上で正であり、これが
$`\mathcal L\circ E`$ の分解（補題 D.2）を与え、それを $`E`$ によって $`\mathcal M`$ に圧縮して戻す。
ハイゼンベルク描像では（系 D.8）、$`\mathcal M`$ 上の単位的 CP 写像の連続な半群の生成子は、$`\Psi`$ を CP、$`K\in\mathcal M`$ として $`A\mapsto\Psi(A)+K^*A+AK`$ である。歴史的には、(1) は
V. Gorini, A. Kossakowski, E. C. G. Sudarshan [*J. Math. Phys.* 17, 821 (1976)] と G. Lindblad
[*Commun. Math. Phys.* 48, 119 (1976)] により、ハイゼンベルク描像の形は、一般のフォン・ノイマン環上のノルム連続な
半群について、E. Christensen と D. E. Evans [*J. London Math. Soc.* 20, 358 (1979)] による。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. 相関が速く減衰する大きな環境と弱く相互作用する系
> （ボルン–マルコフ近似）は、過去の記憶を失うと期待される：そのとき時刻 $`s+t`$ における状態は、
> 時刻 $`s`$ における状態から $`t`$ だけに依存するチャネルによって得られる。*もし* そうであれば、定理 8.7 は、
> そのダイナミクスがハミルトニアン部分と「ジャンプ」$`L_j`$（光子放出、粒子損失、
> 位相緩和）の和であり、そのレートは $`L_j`$ の規格化によって決まることを述べている。

<a id="s-8-3"></a>

## 8.3 セクターをもつ環におけるリンドブラッド方程式

$`\mathcal M`$ 上には二種類のジャンプがある。**セクター内** のジャンプは作用素
$`L_j\in\mathcal M`$ である。**セクター間** のジャンプは、クラウス
作用素 $`K^{dc}_l:\mathbb C^{k_c}\to\mathbb C^{k_d}`$ をもつ CP 写像 $`\Phi:\mathcal M\to\mathcal M`$ である（命題 6.4、`InterSectorChannel`）：
$`\Phi(\rho)_d=\sum_{c,l}K^{dc}_l\rho_cK^{dc*}_l`$。命題 6.8 により、その blunt 双対は
$`\Phi^*(1)_c=\sum_{d,l}K^{dc*}_lK^{dc}_l`$ を満たし、これは $`\mathcal M`$ のブロック対角な元である。レート
$`g_j,g_m>0`$ をもつ `lindblad_rhs` の生成子は
```math
\mathcal L(\rho)=-i[H,\rho]+\sum_jg_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big)
+\sum_mg_m\Big(\Phi_m(\rho)-\tfrac12\{\Phi_m^*(1),\rho\}\Big).
```
である。下の命題 8.8(1) により、これは $`\mathcal M`$ 上の量子動力学的半群を生成する。逆に、
定理 8.7(2)（定理 D.7）により、$`\mathcal M`$ 上の任意の量子動力学的半群はこの形の生成子をもつ。

*ライブラリでは：* セクター間のジャンプは `InterSectorChannel` であり、$`\mathcal L(\rho)`$ は `dynamics.lindblad_rhs` である — [チュートリアル、8.3節](../../tutorial/ja/08_dynamics.md#s-8-3)。

<a id="r-8-8"></a>

**命題 8.8.** $`J^{dc}_l\in M_D`$ を、位置 $`(d,c)`$ にブロック $`K^{dc}_l`$ をもち、
それ以外は零である行列とする。

1. $`\mathcal L`$ は $`\mathcal M`$ を $`\mathcal M`$ に写し、ジャンプ $`\sqrt{g_j}L_j`$ と $`\sqrt{g_m}J^{dc}_{m,l}`$（すべての $`m,l,d,c`$）をもつ $`M_D`$ 上の GKSL 生成子
   $`\hat{\mathcal L}`$ の $`\mathcal M`$ への制限である。したがって
   $`(e^{t\mathcal L})_{t\ge0}`$ は $`\mathcal M`$ 上の量子動力学的半群である。
2. 代わりに $`\tilde{\mathcal L}`$ を、一つのクラウス
   添字のブロックを足し合わせて一つのジャンプ $`L_{m,l}=\sum_{d,c}J^{dc}_{m,l}`$ とした $`M_D`$ 上の GKSL 生成子とする
   （同じハミルトニアン、同じジャンプ $`\sqrt{g_j}L_j`$、およびジャンプ $`\sqrt{g_m}L_{m,l}`$ をもつ）。このとき
   $`\rho\in\mathcal M`$ に対し $`\mathcal L(\rho)=E\big(\tilde{\mathcal L}(\rho)\big)`$ である。ここで $`E(X)=\sum_dP_dXP_d`$ は
   $`M_D`$ から $`\mathcal M`$ の上への条件付き期待値である。一般には $`\tilde{\mathcal L}(\rho)\notin\mathcal M`$ である。
3. すべての $`L_{m,l}`$ が各ブロック行および各ブロック列に高々一つの零でないブロックしかもたないならば（たとえば、
   下の例の消滅作用素のように、セクター $`c`$ をセクター $`c-1`$ にのみ写す場合）、
   $`\rho\in\mathcal M`$ に対し $`\tilde{\mathcal L}(\rho)=\mathcal L(\rho)`$ である。

*証明.* $`\rho=\sum_cP_c\rho P_c`$ かつ $`J^{dc}P_{c'}=\delta_{cc'}J^{dc}`$ であるから、
$`J^{dc}\rho J^{dc*}`$ は位置 $`(d,d)`$ にあるブロック $`K^{dc}\rho_cK^{dc*}`$ であり、$`J^{dc*}J^{dc}`$ は位置 $`(c,c)`$ にある
$`K^{dc*}K^{dc}`$ である。$`d,l`$ について和をとると $`\Phi(\rho)`$ と $`\Phi^*(1)`$ が得られ、これで (1) の
等式が示され、すべての項は $`\mathcal M`$ に属する。したがって $`\hat{\mathcal L}`$ とそのすべてのべきは $`\mathcal M`$
を $`\mathcal M`$ に写し、$`e^{t\hat{\mathcal L}}|_{\mathcal M}=e^{t\mathcal L}`$ であり、これは
定理 8.6 と補題 8.4(a) により $`\mathcal M`$ のチャネルである。(2) $`L=\sum_{d,c}J^{dc}`$ に対し $`L\rho L^*=\sum_{d,d',c}J^{dc}\rho J^{d'c*}`$ であり、
$`E`$ は $`d=d'`$ の項を残す。$`B\in\mathcal M`$ に対し $`P_dXBP_d=P_dXP_dB`$ であるから $`E(XB)=E(X)B`$ かつ
$`E(BX)=BE(X)`$ である。よって $`E\{L^*L,\rho\}=\{E(L^*L),\rho\}`$ であり、$`d\ne d'`$ に対し $`J^{dc*}J^{d'c}=0`$ であるから
$`E(L^*L)=\sum_cP_cL^*LP_c=\sum_{c,d}J^{dc*}J^{dc}`$ である。最後に $`E(-i[H,\rho])=-i[H,\rho]`$ であり、ジャンプ
$`\sqrt{g_j}L_j\in\mathcal M`$ の項は $`\mathcal M`$ に属する。(3) $`L=\sum_{(d,c)\in S}J^{dc}`$ とする。ここで組の集合 $`S`$ は、
与えられた $`d`$ をもつ組を高々一つ、与えられた $`c`$ をもつ組を高々一つしか含まない。$`\rho\in\mathcal M`$ に対し、
$`J^{dc}\rho J^{d'c'*}=J^{dc}P_c\rho P_{c'}J^{d'c'*}`$ は $`c=c'`$ でない限り消え、$`c=c'`$ ならば $`d=d'`$ である。また
$`J^{dc*}J^{d'c'}`$ は $`d=d'`$ でない限り消え、$`d=d'`$ ならば $`c=c'`$ である。したがって $`L\rho L^*`$ と $`L^*L`$ はブロック対角であり、
$`\tilde{\mathcal L}(\rho)`$ のすべての項は $`\mathcal M`$ に属し、(2) により $`\tilde{\mathcal L}(\rho)=E(\tilde{\mathcal L}(\rho))=\mathcal L(\rho)`$
である。$`\square`$

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. (2) は、$`\mathcal M`$ 上のリンドブラッド方程式が、ジャンプによって生じるセクター間の
> コヒーレンスを直ちに消去したときに $`M_D`$ 上の方程式から残るものであることを述べている——
> これが超選択則の内容である。粒子損失の場合には何も消去されない（(3)）：消滅作用素をジャンプとする
> 通常のリンドブラッド方程式は、すでに $`\mathcal M`$ 上の方程式である。

*計算：* セクター間のランダムな写像についての命題 8.8(2) と粒子損失についての命題 8.8(3) を、$`M_6`$ 上の GKSL 生成子と照合する — [チュートリアル、8.3節](../../tutorial/ja/08_dynamics.md#s-8-3)。

<a id="r-8-9"></a>

**命題 8.9（セクターの重み）.** $`\rho(t)=e^{t\mathcal L}\rho_0`$ とし、$`p(t)`$ をそのセクターの重みとする。

1. セクター間のジャンプがなければ（すべての $`m,l`$ と $`d\ne c`$ に対し $`K^{dc}_{m,l}=0`$）、すべての $`p_c(t)`$ は一定である。
2. すべての $`m`$ とすべての組 $`d\ne c`$ に対し、数
   $`r_{m,dc}`$ により $`\sum_lK^{dc*}_{m,l}K^{dc}_{m,l}=r_{m,dc}1_{k_c}`$ となるならば、$`d\ne c`$ に対し $`Q_{dc}=\sum_mg_mr_{m,dc}`$、$`Q_{cc}=-\sum_{d\ne c}Q_{dc}`$ として $`p'=Qp`$ である。
   したがって $`p(t)=e^{tQ}p(0)`$ は $`H`$ にも、$`L_j`$ にも、セクター内の状態にもよらない。

*証明.* 補題 8.1 により $`p_d'=\mathrm{Tr}(P_d\mathcal L(\rho))`$ である。ハミルトニアンの項は定理 8.3 と同様に消える。
$`L_j\in\mathcal M`$ に対しては、$`P_d`$ が中心に属するので $`\mathrm{Tr}(P_dL_j\rho L_j^*)-\mathrm{Tr}(P_dL_j^*L_j\rho) =\mathrm{Tr}(L_j^*P_dL_j\rho)-\mathrm{Tr}(L_j^*L_jP_d\rho)=0`$ である。セクター間の
ジャンプ $`\Phi`$ に対しては $`\mathrm{Tr}(P_d\Phi(\rho))=\sum_c\mathrm{Tr}\big(\sum_lK^{dc*}_lK^{dc}_l\rho_c\big)`$
かつ $`\mathrm{Tr}(P_d\tfrac12\{\Phi^*(1),\rho\})=\sum_{d'}\mathrm{Tr}\big(\sum_lK^{d'd*}_lK^{d'd}_l\rho_d\big)`$ であり、
$`c=d`$ の項と $`d'=d`$ の項は打ち消し合う。これで (1) が示され、(2) の仮定のもとでは残りの項は
$`\sum_{c\ne d}r_{dc}p_c-\sum_{d'\ne d}r_{d'd}p_d`$ である。補題 8.1 により $`p(t)=e^{tQ}p(0)`$ である。$`\square`$

(2) の仮定のもとでは、$`d\ne c`$ に対し、$`r_{m,dc}`$ は $`p_c>0`$ を満たす任意の状態 $`\rho`$ について $`\Phi_m`$ の遷移
行列（命題 6.12 の前で定義した）の成分 $`T_{dc}(\rho)`$ に等しい。
$`\mathrm{Tr}\,\Phi_{m,dc}(\rho_c)=\mathrm{Tr}\big(\sum_lK^{dc*}_{m,l}K^{dc}_{m,l}\rho_c\big)=r_{m,dc}p_c`$ だからである。したがって、レート
$`1`$ の一つのジャンプに対しては、$`R`$ = `transition_matrix()` として $`Q=R-\mathrm{diag}(\sum_dR_{dc})_c`$ である（$`R`$ の対角成分は打ち消し合う）。

*ライブラリでは：* $`R`$ は `InterSectorChannel.transition_matrix` である — [チュートリアル、8.3節](../../tutorial/ja/08_dynamics.md#s-8-3)。

**例：粒子損失.** 粒子数が高々二の二つのボース・モードを考える。セクターは
$`N=0,1,2`$ で、基底はそれぞれ $`|00\rangle`$；$`|10\rangle,|01\rangle`$；$`|20\rangle,|11\rangle,|02\rangle`$ であるから、
$`\mathcal M=M_1\oplus M_2\oplus M_3`$ である。ボース–ハバード二量体
$`H=J(a_1^*a_2+a_2^*a_1)+\frac U2\sum_jn_j(n_j-1)`$ は $`N`$ と可換であり、各モードは
レート $`\gamma`$ で粒子を失う：これはクラウス作用素 $`a_1,a_2`$ を $`N\to N-1`$ に制限した、セクター間の一つのジャンプである。
$`a_1^*a_1+a_2^*a_2=N`$ であるから、命題 8.9(2) が $`Q_{N-1,N}=\gamma N`$ として適用でき、$`N=2`$ から出発すると
```math
p_2(t)=q^2,\qquad p_1(t)=2q(1-q),\qquad p_0(t)=(1-q)^2,\qquad q=e^{-\gamma t},
```
となる。これはパラメータ $`2`$ と $`q`$ の二項分布（各粒子が独立に確率 $`q`$ で生き残るときの生存粒子数の分布）であり、
ホッピングと相互作用が何をしようとも成り立つ：両辺はともに $`p(0)=(0,0,1)`$ を満たす $`p'=Qp`$ の解である
（微分せよ）から、補題 8.1 により一致する。

*計算：* この例。$`\Phi^*(1)=0\oplus1\oplus2`$ を確かめ、$`p(t)`$ を両方の公式と照合する — [チュートリアル、8.3節](../../tutorial/ja/08_dynamics.md#s-8-3)。

<a id="s-8-4"></a>

## 8.4 事後選択とジャンプなしの時間発展

8.3節の生成子のジャンプをまとめて $`\Phi_{\mathrm{tot}}(\rho)=\sum_jg_jL_j\rho L_j^*+\sum_mg_m\Phi_m(\rho)`$ とし、
$`\mathcal M`$ の元 $`A=\Phi_{\mathrm{tot}}^*(1)=\sum_jg_jL_j^*L_j+\sum_mg_m\Phi_m^*(1)\ge0`$ と
$`H_{\mathrm{eff}}=H-\frac i2A`$ をとる。定理 8.6 を命題 8.8(1) の $`M_D`$ 上の生成子 $`\hat{\mathcal L}`$ に適用する：
命題 8.8 の証明により、そのジャンプ部分は $`\mathcal M`$ を $`\mathcal M`$ に写し、そこで $`\Phi_{\mathrm{tot}}`$ に等しく、
ジャンプにわたる和 $`\sum\hat L^*\hat L`$ は $`A`$ であるから、その $`H_{\mathrm{eff}}`$ は $`H-\frac i2A\in\mathcal M`$ である。
したがってその $`S_n(t)`$ は $`\mathcal M`$ を $`\mathcal M`$ に写し（$`n`$ についての帰納法）、$`\rho\in\mathcal M`$ に対し
```math
\rho(t)=\sum_{n\ge0}S_n(t)\rho,\qquad S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},
\qquad S_n(t)=\int_0^tS_0(t-s)\,\Phi_{\mathrm{tot}}\,S_{n-1}(s)\,ds .
```
である。

<a id="r-8-10"></a>

**命題 8.10.** 状態 $`\rho`$ に対し $`w_n(t)=\mathrm{Tr}\,S_n(t)\rho`$ とおく。このとき $`w_n(t)\ge0`$、
$`\sum_nw_n(t)=1`$ であり、$`w_0`$ は非増加で $`w_0'(t)=-\mathrm{Tr}\big(A\,S_0(t)\rho\big)`$ である。
純粋状態 $`\rho=|\psi\rangle\langle\psi|`$ に対しては $`w_0(t)=\|e^{-itH_{\mathrm{eff}}}\psi\|^2`$ である。

*証明.* $`S_n(t)`$ は CP であるから正値であり、$`\sum_nw_n=\mathrm{Tr}\,\rho(t)=1`$ である。
$`X=S_0(t)\rho\ge0`$ とおくと $`w_0'=\mathrm{Tr}(-iH_{\mathrm{eff}}X+iXH_{\mathrm{eff}}^*)=-\mathrm{Tr}(AX)`$
であり、$`\mathrm{Tr}(AX)=\mathrm{Tr}(X^{1/2}AX^{1/2})\ge0`$ である。最後の主張は
$`\mathrm{Tr}|\phi\rangle\langle\phi|=\|\phi\|^2`$ にほかならない。$`\square`$

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. 検出器があらゆるジャンプ——放出されたすべての光子、失われたすべての原子——を記録するとしよう。
> このとき $`w_n(t)`$ は時刻 $`t`$ までにちょうど $`n`$ 回のクリックが起こる確率であり、$`S_n(t)\rho/w_n(t)`$ は
> $`n`$ 回のクリックが起きたという条件のもとでの状態、$`S_0(t)\rho/w_0(t)`$ は **クリックがない** という条件のもとでの状態である：
> 何も起こらなかったという結果も、非エルミートな $`H_{\mathrm{eff}}`$ を通じて状態を変える。この計数過程による解釈 [E. B. Davies,
> *Quantum Theory of Open Systems*, 1976; M. D. Srinivas and E. B. Davies, *Optica Acta* 28, 981 (1981)]
> は量子軌跡法（モンテカルロ波動関数法）の基礎である [J. Dalibard, Y. Castin,
> K. Mølmer, *Phys. Rev. Lett.* 68, 580 (1992)]。ここではその数学的な部分だけを用いる。

損失の場合、ダイソン項は $`\rho(t)`$ のセクター成分である（囲みの言葉でいえば、クリックを数えることは
電荷を測定することと同じである）。$`\mathrm{Tr}(P\rho)>0`$ を満たす射影 $`P`$ による状態 $`\rho`$ の
**リューダース更新** とは、状態 $`P\rho P/\mathrm{Tr}(P\rho P)`$ のことである（$`\mathrm{Tr}(P\rho P)=\mathrm{Tr}(P\rho)`$ に注意）。

*ライブラリでは：* リューダース更新は `states.lueders_update` である — [チュートリアル、8.4節](../../tutorial/ja/08_dynamics.md#s-8-4)。

<a id="r-8-11"></a>

**命題 8.11.** セクターに $`N=0,\dots,C-1`$ とラベルをつけ、セクター内の
ジャンプはなく、すべてのジャンプが電荷を一つ下げ（$`d=c-1`$ でない限り $`K^{dc}_{m,l}=0`$）、
$`\rho=P_{N_0}\rho P_{N_0}`$ であるとする。このとき $`0\le n\le N_0`$ に対し $`S_n(t)\rho=P_{N_0-n}\,\rho(t)\,P_{N_0-n}`$ であり、
$`n>N_0`$ に対し $`S_n(t)\rho=0`$ である。特に $`p_{N_0}(t)=w_0(t)>0`$ であり、$`\rho(t)`$ の $`P_{N_0}`$ によるリューダース更新（「粒子が一つも失われなかった」）
はジャンプなしの状態である：
```math
\frac{P_{N_0}\rho(t)P_{N_0}}{p_{N_0}(t)}=\frac{e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*}}{w_0(t)},
\qquad p_{N_0}(t)=w_0(t).
```

*証明.* $`H_{\mathrm{eff}}\in\mathcal M`$ であるから、$`S_0(s)`$ はセクター $`N`$ に台をもつ作用素をセクター $`N`$ に台をもつ作用素に写す。
$`\Phi_{\mathrm{tot}}`$ はセクター $`N\ge1`$ をセクター $`N-1`$ に写し、セクター $`0`$ を零にする。
帰納法により、すべての $`s`$ に対し $`S_n(s)\rho`$ はセクター $`N_0-n`$ に台をもち、$`n>N_0`$ ならば $`0`$ である。
$`\rho(t)=\sum_nS_n(t)\rho`$ を $`P_{N_0-n}`$ で圧縮すると、ちょうど第 $`n`$ 項が残る。最後に、
$`e^{-itH_{\mathrm{eff}}}`$ は可逆であるから（補題 8.1）$`S_0(t)\rho`$ は零でない正の行列であり、$`w_0(t)>0`$ である。$`\square`$

*計算：* 異なるレートで漏れ出す二つのモードについて、$`n=0`$ と $`n=1`$ に対する命題 8.11 — [チュートリアル、8.4節](../../tutorial/ja/08_dynamics.md#s-8-4)。

<a id="s-8-5"></a>

## 8.5 定常状態と緩和

状態 $`\sigma`$ は、すべての $`t\ge0`$ に対し $`e^{t\mathcal L}\sigma=\sigma`$ となるとき **定常** であるという。これは
$`\mathcal L(\sigma)=0`$ と同値である（$`t=0`$ で微分すればよく、逆は級数を用いる）。

<a id="r-8-12"></a>

**命題 8.12.** $`\mathcal M`$ 上の任意の量子動力学的半群は定常状態をもつ。

*証明.* 状態 $`\rho`$ を固定し、$`\sigma_\tau=\frac1\tau\int_0^\tau T_s\rho\,ds`$ とおく。補題 A.4(d) により、これは状態である
（状態全体は閉凸集合をなす）。補題 A.4(b),(e) により
$`T_t\sigma_\tau-\sigma_\tau=\frac1\tau\big(\int_\tau^{\tau+t}-\int_0^t\big)T_s\rho\,ds`$ であり、$`\|T_s\rho\|_1=1`$ であるから、補題 A.4(c) により
そのトレースノルムは高々 $`2t/\tau`$ である。状態全体は有限次元空間の有界閉集合であるから、
ある列 $`\sigma_{\tau_m}`$（$`\tau_m\to\infty`$）は状態 $`\sigma`$ に収束し
（ボルツァーノ–ワイエルシュトラス）、$`T_t\sigma-\sigma=\lim_m(T_t\sigma_{\tau_m}-\sigma_{\tau_m})=0`$ である。$`\square`$

<a id="r-8-13"></a>

**系 8.13.** セクター間のジャンプがなければ（命題 8.9(1) のように）、定常状態の張る空間の
次元は少なくとも $`C`$ である：各セクターがそれ自身の定常状態をもつ。

*証明.* $`\mathcal L`$ は $`P_c\mathcal MP_c\cong M_{k_c}`$ をそれ自身に写す（$`H`$ と $`L_j`$ は $`\mathcal M`$ に属し、
$`\Phi_m(X)=\sum_lK^{cc}_{m,l}X_cK^{cc*}_{m,l}`$ と $`\Phi_m^*(1)`$ はブロック対角である）から、すべての $`e^{t\mathcal L}`$ も同様である。
その制限は CP（補題 8.4(a)）かつ TP であるから、$`M_{k_c}`$ 上の量子動力学的半群であり、命題 8.12 により
定常状態 $`\sigma_c`$ をもつ。$`\mathcal M`$ の元として、$`\sigma_c`$ は $`e^{t\mathcal L}`$ の定常状態である。
$`\sigma_c`$ たちは異なるブロックにあるから、一次独立である。$`\square`$

定理 D.6（A. Frigerio の判定条件の有限次元の場合）は一意性のための十分条件を与える：データ $`(H;L_1,\dots,L_r)`$ をもつ
$`M_k`$ 上の GKSL 生成子が忠実な定常状態 $`\sigma>0`$ をもち、$`H`$ およびすべての $`L_j`$、$`L_j^*`$ と可換な行列が
$`1`$ の定数倍だけならば、$`\sigma`$ はただ一つの定常状態である。セクター間のジャンプがなければ、$`\mathcal L`$ の
$`P_c\mathcal MP_c\cong M_{k_c}`$ への制限はデータ $`(H_c;\sqrt{g_1}L_{1,c},\sqrt{g_2}L_{2,c},\dots)`$ をもつ GKSL 生成子であるから、
定理はセクターごとに適用できる。演習問題 4 は、
因子環においてさえ一意性が成り立たないことがあることを示す：そこでは第二の仮定が成り立たない。

*計算：* $`M_2\oplus M_3`$ のセクター内のランダムなジャンプについての系 8.13 と定理 D.6 の仮定 — [チュートリアル、8.5節](../../tutorial/ja/08_dynamics.md#s-8-5)。

<a id="r-8-14"></a>

**命題 8.14（自然放出）.** 基底 $`|g\rangle,|e\rangle`$ をもつ $`M_2`$ 上で、
$`H=\omega|e\rangle\langle e|`$、一つのジャンプ $`L=|g\rangle\langle e|`$、レート $`\gamma>0`$ とする。このとき
```math
\rho_{ee}(t)=e^{-\gamma t}\rho_{ee}(0),\qquad\rho_{ge}(t)=e^{(i\omega-\gamma/2)t}\rho_{ge}(0),
\qquad e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}},
```
である。ここで $`\mathcal U_t(\rho)=e^{-itH}\rho e^{itH}`$ であり、$`\mathcal A_p`$ はパラメータ
$`p`$ の振幅減衰チャネル（6.4節）である。$`\mathcal L`$ のスペクトルは $`\{0,-\gamma,-\gamma/2\pm i\omega\}`$ であり、
$`|g\rangle\langle g|`$ がただ一つの定常状態であって、
$`\tfrac12\|\rho(t)-|g\rangle\langle g|\|_1=\big(\rho_{ee}(0)^2e^{-2\gamma t}+|\rho_{ge}(0)|^2e^{-\gamma t}\big)^{1/2}`$ である。

*証明.* ここで $`L^*L=|e\rangle\langle e|`$、$`L\rho L^*=\rho_{ee}|g\rangle\langle g|`$ であるから、成分ごとに
$`\rho_{ee}'=-\gamma\rho_{ee}`$、$`\rho_{gg}'=\gamma\rho_{ee}`$、$`\rho_{ge}'=(i\omega-\gamma/2)\rho_{ge}`$ であり、
これから解が得られる。基底 $`(\rho_{gg},\rho_{ge},\rho_{eg},\rho_{ee})`$ において $`\mathcal L`$ の行列は
上三角で、対角成分は $`(0,-\gamma/2+i\omega,-\gamma/2-i\omega,-\gamma)`$ である。$`\mathcal L(\sigma)=0`$ からは
$`\sigma_{ee}=\sigma_{ge}=0`$ が従う。振幅減衰チャネルは $`\rho_{ee}`$ に $`1-p`$ を、$`\rho_{ge}`$ に
$`\sqrt{1-p}`$ を掛け（そのクラウス作用素は $`|g\rangle\langle g|+\sqrt{1-p}\,|e\rangle\langle e|`$ と
$`\sqrt p\,|g\rangle\langle e|`$）、$`\mathcal U_t`$ は $`\rho_{ge}`$ に $`e^{i\omega t}`$ を掛ける。どちらも TP であり、これで
$`\rho_{gg}`$ も決まる。最後に、$`\rho(t)-|g\rangle\langle g|`$ はトレース零の自己共役作用素で、対角成分は
$`(-\rho_{ee}(t),\rho_{ee}(t))`$ であるから、その固有値は $`\pm(\rho_{ee}(t)^2+|\rho_{ge}(t)|^2)^{1/2}`$ である。$`\square`$

コヒーレンスは占有数の半分のレートで減衰する：$`T_1=1/\gamma`$、$`T_2=2/\gamma`$ とすると、これは 6.4節の関係
$`\frac1{T_2}=\frac1{2T_1}+\frac1{T_\varphi}`$ から位相緩和の項を除いたものである。最も遅い非零の減衰レート、
すなわち **スペクトルギャップ** $`\min\{-\mathrm{Re}\,\mu:\mu\ne0\text{ は }\mathcal L\text{ の固有値}\}=\gamma/2`$ が、
ここでの緩和のレートである：最後の公式により、$`|g\rangle\langle g|`$ までの距離は高々 $`e^{-\gamma t/2}`$ であり
（正値性により $`|\rho_{ge}|^2\le\rho_{gg}\rho_{ee}`$ であるから $`\rho_{ee}^2+|\rho_{ge}|^2\le\rho_{ee}\le1`$）、その $`e^{\gamma t/2}`$ 倍は
$`t\to\infty`$ のとき $`|\rho_{ge}(0)|`$ に収束する。

*計算：* `dynamics.lindblad_superoperator` と `dynamics.lindblad_channel` による、$`\mathcal L`$ のスペクトル、分解 $`e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}}`$、および命題 8.14 のトレース距離 — [チュートリアル、8.5節](../../tutorial/ja/08_dynamics.md#s-8-5)。

<a id="s-8-6"></a>

## 8.6 比較の対象

*計算：* 他のソフトウェア（QuTiP、dynamiqs、QuantumOptics.jl、SciPy）との比較 — [チュートリアル、8.6節](../../tutorial/ja/08_dynamics.md#s-8-6)。

<a id="s-8-7"></a>

## 8.7 チャネルである時間ステップ

数値ソルバーは $`e^{t\mathcal L}`$ を写像 $`R_h`$ の $`n`$ ステップで置き換える（$`h=t/n`$）。各 $`e^{h\mathcal L}`$ は
$`\mathcal M`$ のチャネルである（命題 8.8(1)）。各ステップもチャネルであれば、$`h`$ が何であっても、計算された密度は
状態である。$`\rho'=\mathcal L(\rho)`$ に対する古典的ルンゲ–クッタ法（RK4）は、$`k_1=\mathcal L(\rho)`$、
$`k_2=\mathcal L(\rho+\frac h2k_1)`$、$`k_3=\mathcal L(\rho+\frac h2k_2)`$、$`k_4=\mathcal L(\rho+hk_3)`$ によるステップ
$`\rho\mapsto\rho+\frac h6(k_1+2k_2+2k_3+k_4)`$ である。各段を次の段に代入すると
$`k_2=(\mathcal L+\frac h2\mathcal L^2)\rho`$、$`k_3=(\mathcal L+\frac h2\mathcal L^2+\frac{h^2}4\mathcal L^3)\rho`$、
$`k_4=(\mathcal L+h\mathcal L^2+\frac{h^2}2\mathcal L^3+\frac{h^3}4\mathcal L^4)\rho`$ となるので
```math
R_h=T_4(h\mathcal L),\qquad T_4(z)=1+z+\frac{z^2}2+\frac{z^3}6+\frac{z^4}{24}.
```
すべての $`X`$ に対して $`\mathrm{Tr}\,\mathcal L(X)=0`$ であるから（定理 8.6 の証明のステップ 5 と命題 8.8(1)）、$`R_h`$ は
TP であり、$`\mathcal L(X^*)=\mathcal L(X)^*`$ であるから自己共役元を自己共役元に写す。$`R_h`$ は $`e^{h\mathcal L}`$ と指数級数の
$`\ge5`$ 次の項だけ異なるが、正値であるとは限らない。

<a id="r-8-15"></a>

**命題 8.15（RK4 は正値でない）.** 基底 $`|g\rangle,|e\rangle`$ をもつ $`M_2`$ 上で、$`\mathcal L`$ を命題 8.14 の生成子
（$`H=\omega|e\rangle\langle e|`$、$`\omega\in\mathbb R`$、一つのジャンプ $`|g\rangle\langle e|`$、レート $`\gamma>0`$）とし、
$`\lambda=i\omega-\gamma/2`$、$`h>0`$、$`\tau=T_4(-h\gamma)`$、$`a=T_4(h\lambda)`$ とする。

1. $`R_h(\rho)=\begin{pmatrix}\rho_{gg}+(1-\tau)\rho_{ee}&a\rho_{ge}\\\bar a\rho_{eg}&\tau\rho_{ee}\end{pmatrix}`$ であり、次は
   同値である：$`R_h`$ は正値；$`R_h`$ は CP；$`\tau\le1`$ かつ $`|a|^2\le\tau`$。
2. $`\tau\le1`$ であるのは $`h\gamma\le x_*`$ のとき、かつそのときに限る。ここで $`x_*=2.7853\ldots`$ は
   $`x^3-4x^2+12x-24`$ のただ一つの実根である。$`h\gamma>x_*`$ のとき、状態 $`|e\rangle\langle e|`$ は
   $`\mathrm{diag}(1-\tau,\tau)`$ に写され、その第一成分は負である。
3. $`u=\gamma^2+4\omega^2`$ とおくと
   ```math
   \tau-|a|^2=h^5\Big(-\frac{\gamma u(4\omega^2-3\gamma^2)}{384}+\frac{u^2(4\omega^2-5\gamma^2)}{4608}\,h
   +\frac{\gamma u^3}{9216}\,h^2-\frac{u^4}{147456}\,h^3\Big).
   ```
   したがって、$`2\omega>\sqrt3\,\gamma`$ ならば、$`0<h<h_0`$ で $`R_h`$ が正値でないような $`h_0>0`$ が存在する。
   $`\omega=\gamma`$ ならば、どの $`h>0`$ に対しても $`R_h`$ は正値でない。$`\omega=0`$ ならば、$`R_h`$ が CP であるのは
   $`h\gamma\le x_*`$ のとき、かつそのときに限る。

*証明.* (1) 命題 8.14 の証明により、$`\mathcal L(|g\rangle\langle g|)=0`$、
$`\mathcal L(|e\rangle\langle e|)=\gamma(|g\rangle\langle g|-|e\rangle\langle e|)`$、$`\mathcal L(|g\rangle\langle e|)=\lambda|g\rangle\langle e|`$、
$`\mathcal L(|e\rangle\langle g|)=\bar\lambda|e\rangle\langle g|`$ である。帰納法により $`n\ge1`$ に対して
$`\mathcal L^n(|e\rangle\langle e|)=(-\gamma)^n(|e\rangle\langle e|-|g\rangle\langle g|)`$ であるから
$`R_h(|e\rangle\langle e|)=|e\rangle\langle e|+(\tau-1)(|e\rangle\langle e|-|g\rangle\langle g|)`$ となり、一方 $`R_h`$ は
$`|g\rangle\langle g|`$ を固定し、$`|g\rangle\langle e|`$ に $`T_4(h\lambda)=a`$ を、$`|e\rangle\langle g|`$ に
$`T_4(h\bar\lambda)=\bar a`$ を掛ける（$`T_4`$ の係数は実数）。これが公式である。$`\tau\le1`$ かつ $`|a|^2\le\tau`$ ならば、
クラウス作用素 $`|g\rangle\langle g|+\bar a|e\rangle\langle e|`$、$`\sqrt{1-\tau}\,|g\rangle\langle e|`$、
$`\sqrt{\tau-|a|^2}\,|e\rangle\langle e|`$ が $`R_h`$ を与える（四つの成分を比較せよ）ので、$`R_h`$ は CP（定理 6.2）、特に
正値である。逆に $`R_h`$ が正値であるとする。$`R_h(|e\rangle\langle e|)\ge0`$ から $`1-\tau\ge0`$ を得る。$`c,s>0`$、
$`c^2+s^2=1`$ である $`\psi=c|g\rangle+s|e\rangle`$ に対し、行列 $`R_h(|\psi\rangle\langle\psi|)\ge0`$ の行列式は
$`s^2\big(c^2(\tau-|a|^2)+\tau(1-\tau)s^2\big)\ge0`$ である。$`s\to0`$ とすると $`\tau-|a|^2\ge0`$ を得る。

(2) $`x=h\gamma`$ とおくと $`\tau-1=\frac x{24}q(x)`$、$`q(x)=x^3-4x^2+12x-24`$ である。$`q'(x)=3x^2-8x+12>0`$（判別式は
$`64-144<0`$）であるから $`q`$ は狭義単調増加で、$`q(2)=-8<0<3=q(3)`$ である。よって $`q`$ はちょうど一つの実根
$`x_*\in(2,3)`$ をもち、$`\tau\le1`$ であるのは $`x\le x_*`$ のとき、かつそのときに限る。

(3) $`\tau-|a|^2=T_4(-h\gamma)-T_4(h\lambda)T_4(h\bar\lambda)`$ は $`h`$ の $`8`$ 次多項式である。
$`T_4(h\lambda)T_4(h\bar\lambda)`$ における $`h^n`$ の係数は、$`j+l=n`$、$`0\le j,l\le4`$ にわたる和
$`\sum\lambda^j\bar\lambda^l/(j!\,l!)`$ である。$`n\le4`$ ではそのような組がすべて現れ、二項定理により和は
$`(\lambda+\bar\lambda)^n/n!=(-\gamma)^n/n!`$、すなわち $`T_4(-h\gamma)`$ における $`h^n`$ の係数である。よって
$`\tau-|a|^2`$ の $`h^0,\dots,h^4`$ の係数は消える。$`n=5,\dots,8`$ に対しては、$`\sigma=\lambda+\bar\lambda=-\gamma`$ と
$`\pi=\lambda\bar\lambda=u/4`$ を用いると、和はそれぞれ
$`\frac{\pi(\lambda^3+\bar\lambda^3)}{24}+\frac{\pi^2\sigma}{12}=\frac{\pi\sigma^3-\pi^2\sigma}{24}`$、
$`\frac{\pi^2(\lambda^2+\bar\lambda^2)}{48}+\frac{\pi^3}{36}=\frac{\pi^2\sigma^2}{48}-\frac{\pi^3}{72}`$、$`\frac{\pi^3\sigma}{144}`$、
$`\frac{\pi^4}{576}`$ である（$`\lambda^2+\bar\lambda^2=\sigma^2-2\pi`$、$`\lambda^3+\bar\lambda^3=\sigma^3-3\pi\sigma`$ を用いる）。
$`T_4(-h\gamma)`$ にはそのような項はない。$`\sigma`$ と $`\pi`$ を代入すると公式を得る。$`2\omega>\sqrt3\gamma`$ ならば
$`h^5`$ の係数は負であるから、小さいすべての $`h>0`$ に対して $`\tau-|a|^2<0`$ であり、(1) により $`R_h`$ は正値でない。
$`\omega=\gamma`$ ならば $`u=5\gamma^2`$ であり、$`y=5h\gamma`$ とおくと公式は
$`\tau-|a|^2=-\frac{5(h\gamma)^5}{147456}\,p(y)`$、$`p(y)=y^3-16y^2+32y+384`$ となる。ここで $`y>0`$ に対し $`p(y)>0`$ である：
$`y\le12`$ では $`p(y)=y(y-8)^2+(384-32y)`$ は同時には消えない二つの非負の項の和である。$`12\le y\le16`$ では
$`y(y-8)^2\ge192`$（$`y\ge8`$ で増加）かつ $`384-32y\ge-128`$ である。$`y\ge16`$ では $`p(y)=y^2(y-16)+32y+384>0`$ である。
よってすべての $`h>0`$ に対して $`\tau<|a|^2`$ である。$`\omega=0`$ ならば $`u=\gamma^2`$ であり、$`x=h\gamma`$ とおくと
$`0<x\le3`$ に対して
$`\tau-|a|^2=x^5\big(\frac1{128}-\frac{5x}{4608}+\frac{x^2}{9216}-\frac{x^3}{147456}\big)\ge
x^5\big(\frac1{128}-\frac{15}{4608}-\frac{27}{147456}\big)>0`$ である。(1) と (2) により、$`R_h`$ が CP であるのは
$`x\le x_*`$ のとき、かつそのときに限る。$`\square`$

したがって、弱く減衰する二準位系（$`2\omega>\sqrt3\,\gamma`$）では、ステップをどれほど小さくしても RK4 のステップは正値で
ない。一方、減衰だけ（$`\omega=0`$、(3) による）でも回転だけでも、小さいステップではこの欠陥は生じない：$`\gamma=0`$ では
(1) の証明が $`\tau=1`$ で適用でき、$`x^2\le8`$ ならば $`|T_4(ix)|^2=1-\frac{x^6}{72}+\frac{x^8}{576}\le1`$ である（展開せよ）。
破れは小さく $`\tau-|a|^2=O(h^5)`$ であるが、丸め誤差の効果ではなく、$`R_h`$ のトレース保存性によっても解消されない。RK4 の
安定性理論（$`|T_4(z)|\le1`$ となる領域）もこれを捉えない：ここでは $`h`$ が小さいとき $`h\lambda`$ と $`-h\gamma`$ はその領域に
ある。

*計算：* `lindblad_evolve(..., method='rk4')` による命題 8.15：$`T_4(h\mathcal L)`$ としてのステップ、負の固有値をもつ行列に
写される純粋状態、および 8.3 節の粒子損失に対する負の固有値をもつチョイ行列 — [チュートリアル、8.7節](../../tutorial/ja/08_dynamics.md#s-8-7)。

次のスキームは、どの刻み幅に対してもチャネルである。これは 8.4 節の分解 $`\mathcal L=\mathcal L_0+\Phi_{\mathrm{tot}}`$、
$`\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)`$（$`H_{\mathrm{eff}}=H-\frac i2A`$、$`A=\Phi_{\mathrm{tot}}^*(1)`$）と、
$`\mathcal M`$ 上の超演算子 $`\Theta`$ のトレースに関する双対 $`\Theta^*`$、$`\mathrm{Tr}(\Theta(X)Z)=\mathrm{Tr}(X\Theta^*(Z))`$ を
用いる。$`\mathcal M`$ 上の超演算子 $`\Theta`$ に対して、1.4 節のトレースノルムを用いて
$`\|\Theta\|_1=\max\{\|\Theta(X)\|_1:X\in\mathcal M,\ \|X\|_1\le1\}`$ と書く。このノルムは劣乗法的である。

<a id="r-8-16"></a>

**定理 8.16（各ステップがチャネルである分割法）.** $`h\ge0`$ に対して $`E_h=e^{-ihH_{\mathrm{eff}}/2}`$、
$`S_h(\rho)=E_h\rho E_h^*`$ とし、
```math
\Psi_h=S_h\circ\Big(\mathrm{id}+h\Phi_{\mathrm{tot}}+\frac{h^2}2\Phi_{\mathrm{tot}}^2\Big)\circ S_h,\qquad M_h=\Psi_h^*(1),\qquad
\tilde\Psi_h(\rho)=\Psi_h\big(M_h^{-1/2}\rho M_h^{-1/2}\big)
```
とおく。

1. $`M_h\in\mathcal M`$ かつ $`M_h\ge e^{ihH_{\mathrm{eff}}^*}e^{-ihH_{\mathrm{eff}}}`$ であり、右辺は正定値である。$`\tilde\Psi_h`$ は
   すべての $`h\ge0`$ に対して $`\mathcal M`$ のチャネルである。
2. （*局所誤差.*）$`\mathcal L`$ に依存する定数 $`C`$ が存在して、$`0\le h\le1`$ に対して
   $`\|\tilde\Psi_h-e^{h\mathcal L}\|_1\le Ch^3`$ が成り立つ。
3. （*大域誤差.*）すべての状態 $`\rho`$、すべての $`t>0`$、すべての整数 $`n\ge t`$ に対して、$`h=t/n`$ とおくと
   $`\|\tilde\Psi_h^n(\rho)-e^{t\mathcal L}\rho\|_1\le Cth^2`$ が成り立つ。

*証明.* *ステップ 0.* $`\Theta`$ が $`\mathcal M`$ 上の正値かつ TP な写像で $`Y=Y^*\in\mathcal M`$ ならば
$`\|\Theta(Y)\|_1\le\|Y\|_1`$ である。実際、$`Y=Y_+-Y_-`$ と書く。ここで $`Y_\pm\ge0`$ は、スペクトル分解のうち正の固有値に
対応する部分と負の固有値に対応する部分である。$`Y`$ の特異値はその固有値の絶対値である（$`|Y|=(Y^2)^{1/2}`$）から
$`\|Y\|_1=\mathrm{Tr}\,Y_++\mathrm{Tr}\,Y_-`$ であり、三角不等式と $`P\ge0`$ に対する $`\|P\|_1=\mathrm{Tr}\,P`$（1.4節）により
$`\|\Theta(Y)\|_1\le\mathrm{Tr}\,\Theta(Y_+)+\mathrm{Tr}\,\Theta(Y_-)=\|Y\|_1`$ である。任意の $`X\in\mathcal M`$ は、トレースノルムが
$`\le\|X\|_1`$ である自己共役な $`Y_1=\frac12(X+X^*)`$、$`Y_2=\frac1{2i}(X-X^*)`$ により $`Y_1+iY_2`$ と書ける
（$`\|X^*\|_1=\|X\|_1`$、第1章、演習問題 2）。よって $`\|\Theta\|_1\le2`$ である。

*ステップ 1：(1).* $`H_{\mathrm{eff}}\in\mathcal M`$（8.4節）であるから $`E_h\in\mathcal M`$ である（その指数級数の項はすべて部分空間
$`\mathcal M`$ にある）。$`\Phi_{\mathrm{tot}}`$ は $`\mathcal M`$ を $`\mathcal M`$ に写し、クラウス作用素 $`\sqrt{g_j}L_j`$ と
$`\sqrt{g_m}K^{dc}_{m,l}`$ をもつ CP 写像である（命題 6.4）。よって $`\mathrm{id}+h\Phi_{\mathrm{tot}}+\frac{h^2}2\Phi_{\mathrm{tot}}^2`$
は CP（補題 8.4(a)）、$`S_h`$ はクラウス作用素 $`E_h`$ をもつ CP 写像であり、$`\Psi_h`$ は CP である。トレースの巡回性により
$`S_h^*(Z)=E_h^*ZE_h`$ であり、
$`\Psi_h^*=S_h^*\circ(\mathrm{id}+h\Phi_{\mathrm{tot}}^*+\frac{h^2}2(\Phi_{\mathrm{tot}}^*)^2)\circ S_h^*`$ である。ここで
$`\Phi_{\mathrm{tot}}^*(Z)=\sum_jg_jL_j^*ZL_j+\sum_mg_m\Phi_m^*(Z)`$（演習問題 1）は正値である（命題 6.8：$`\Phi_m^*`$ のクラウス
作用素は $`K^{dc*}_{m,l}`$ である）。$`W=S_h^*(1)=E_h^*E_h\ge0`$ とすると、中央の因子は $`W`$ を $`W`$ と正の元の和に写し、
$`S_h^*`$ は順序を保つから、$`M_h\ge S_h^*(W)=(E_h^2)^*E_h^2=e^{ihH_{\mathrm{eff}}^*}e^{-ihH_{\mathrm{eff}}}`$ である
（補題 8.1 により $`E_h^2=e^{-ihH_{\mathrm{eff}}}`$）。この行列は正定値である：$`E_h^2`$ は逆 $`e^{ihH_{\mathrm{eff}}}`$ をもつ
（補題 8.1）から、$`x\ne0`$ に対し $`\langle x,(E_h^2)^*E_h^2x\rangle=\|E_h^2x\|^2>0`$ である。よって自己共役な
$`M_h\in\mathcal M`$ の固有値はすべて正であり、$`N_h=M_h^{-1/2}\in\mathcal M`$（関数計算、1.3節）は $`N_hM_hN_h=1`$ を満たす。写像
$`\rho\mapsto N_h\rho N_h`$ は CP であるから $`\tilde\Psi_h`$ も CP であり（補題 8.4(a)）、
$`\mathrm{Tr}\,\tilde\Psi_h(\rho)=\mathrm{Tr}(N_h\rho N_hM_h)=\mathrm{Tr}(\rho N_hM_hN_h)=\mathrm{Tr}\,\rho`$ である。

*ステップ 2：テイラー係数.* $`\mathcal M`$ 上の超演算子 $`\alpha=\frac12\mathcal L_0`$、$`\beta=\Phi_{\mathrm{tot}}`$ をとる。定理 8.6
の証明のステップ 1 と同様に、$`X(s)=E_s\rho E_s^*`$ は $`X'=\alpha(X)`$ を満たす（補題 8.1 により
$`E_s=e^{-isH_{\mathrm{eff}}/2}`$ を微分する）ので、補題 8.1 により $`S_h=e^{h\alpha}`$ であり、
$`\Psi_h=e^{h\alpha}(1+h\beta+\frac{h^2}2\beta^2)e^{h\alpha}`$、一方 $`e^{h\mathcal L}=e^{h(2\alpha+\beta)}`$ である。絶対収束する
級数を掛け合わせると $`\Psi_h=\sum_{n\ge0}h^nP_n`$、$`e^{h\mathcal L}=\sum_{n\ge0}h^nQ_n`$ であり、
$`\sum_n\|P_n\|_1\le e^{2\|\alpha\|_1}(1+\|\beta\|_1+\frac12\|\beta\|_1^2)`$、$`\sum_n\|Q_n\|_1\le e^{\|2\alpha+\beta\|_1}`$ である。
ここで $`P_0=Q_0=\mathrm{id}`$、$`P_1=2\alpha+\beta=Q_1`$、
$`P_2=\frac{\alpha^2}2+\alpha\beta+\alpha^2+\beta\alpha+\frac{\beta^2}2+\frac{\alpha^2}2=\frac12(2\alpha+\beta)^2=Q_2`$ である。よって
$`0\le h\le1`$ に対し、$`C_1=\sum_{n\ge3}(\|P_n\|_1+\|Q_n\|_1)`$ として $`\|\Psi_h-e^{h\mathcal L}\|_1\le C_1h^3`$ であり、
$`\|\Psi_h\|_1\le K=\sum_n\|P_n\|_1`$ である。

*ステップ 3：(2).* $`x`$ を $`M_h`$ の固有値 $`\mu`$ に対する単位固有ベクトルとする。$`M_h`$ はブロック対角であるから $`x`$ は
一つのセクターの中にとれ、$`X=|x\rangle\langle x|\in\mathcal M`$、$`\|X\|_1=1`$ である。$`e^{h\mathcal L}`$ は TP であるから
$`\mu-1=\mathrm{Tr}(XM_h)-\mathrm{Tr}\,X=\mathrm{Tr}\,\Psi_h(X)-\mathrm{Tr}\,e^{h\mathcal L}(X)`$ であり、$`|\mathrm{Tr}\,Y|\le\|Y\|_1`$
（1.4節）から $`|\mu-1|\le C_1h^3`$ を得る。$`h_1=\min\{1,(2C_1)^{-1/3}\}`$ とおく。$`h\le h_1`$ ならばすべての
$`\mu\ge\frac12`$ であり、$`N_h`$ の固有値 $`\mu^{-1/2}\le\sqrt2`$ は
$`|\mu^{-1/2}-1|=\frac{|1-\mu|}{\sqrt\mu(1+\sqrt\mu)}\le|\mu-1|`$ を満たす（$`\sqrt\mu(1+\sqrt\mu)\ge2^{-1/2}(1+2^{-1/2})>1`$ であるから）。
自己共役行列の作用素ノルムは固有値の絶対値の最大値であるから、$`\|N_h-1\|\le C_1h^3`$、$`\|N_h\|\le\sqrt2`$ である。ここで
```math
\tilde\Psi_h(X)-e^{h\mathcal L}(X)=\Psi_h\big((N_h-1)XN_h+X(N_h-1)\big)+\big(\Psi_h-e^{h\mathcal L}\big)(X)
```
であり、$`\|AB\|_1\le\|A\|\,\|B\|_1`$、$`\|BA\|_1\le\|A\|\,\|B\|_1`$（第1章、演習問題 2。後者は随伴をとって得る）により、
$`h\le h_1`$ に対して $`\|\tilde\Psi_h-e^{h\mathcal L}\|_1\le\big((1+\sqrt2)K+1\big)C_1h^3`$ である。$`h_1<h\le1`$ では両方の写像が
正値かつ TP であるから、ステップ 0 により $`\|\tilde\Psi_h-e^{h\mathcal L}\|_1\le4\le4h^3/h_1^3`$ である。
$`C=\max\{((1+\sqrt2)K+1)C_1,\,4h_1^{-3}\}`$ とすればよい。

*ステップ 4：(3).* 和
$`\tilde\Psi_h^n-e^{nh\mathcal L}=\sum_{k=0}^{n-1}\tilde\Psi_h^{\,n-1-k}\big(\tilde\Psi_h-e^{h\mathcal L}\big)e^{kh\mathcal L}`$
は望遠鏡和である。$`\rho`$ に適用すると、第 $`k`$ 項は $`\tilde\Psi_h^{\,n-1-k}(Y_k)`$ であり、
$`Y_k=(\tilde\Psi_h-e^{h\mathcal L})(e^{kh\mathcal L}\rho)`$ は自己共役で（クラウス形により CP 写像は自己共役性を保つ）、状態
$`e^{kh\mathcal L}\rho`$ のトレースノルムは $`1`$ であるから (2) により $`\|Y_k\|_1\le Ch^3`$ である。(1) とステップ 0 により
$`\|\tilde\Psi_h^{\,n-1-k}(Y_k)\|_1\le Ch^3`$ であり、$`n`$ 個の項の和は高々 $`nCh^3=Cth^2`$ である。$`\square`$

ステップがチャネルであることが効いてくるのはステップ 4 である：各ステップで生じた誤差は、$`n`$ が何であっても決して
増幅されず、$`h`$ に対する安定性条件は必要ない。定義を展開すると、$`\tilde\Psi_h`$ のクラウス作用素は $`E_h^2N_h`$、
$`\sqrt h\,E_hJE_hN_h`$、$`\frac h{\sqrt2}E_hJ'JE_hN_h`$ である。ここで $`J,J'`$ は $`\Phi_{\mathrm{tot}}`$ のクラウス作用素を
動く：8.4 節のジャンプなしの時間発展による半ステップ、高々二回のジャンプ、残りの半ステップ、そして確率の補正 $`N_h`$
（ステップ 3 により $`\|N_h-1\|\le C_1h^3`$）である。

*ライブラリでは：* $`\tilde\Psi_h`$ は `lindblad_evolve(..., method='kraus')` である。`method='exact'` は $`\mathcal M`$ 上の
$`\mathcal L`$ の行列（サイズ $`\sum_ck_c^2`$、セクター間のジャンプを含む）から $`e^{h\mathcal L}`$ そのものを計算し、そのステップは
命題 8.8(1) によりチャネルである。既定の `method='rk4'` は $`R_h`$ である — [チュートリアル、8.7節](../../tutorial/ja/08_dynamics.md#s-8-7)。

<a id="x-8"></a>

## 演習問題

1. （ハイゼンベルク描像）8.3節の生成子について、
   $`\mathcal L^*(A)=i[H,A]+\sum_jg_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})+\sum_mg_m(\Phi_m^*(A)-\frac12\{\Phi_m^*(1),A\})`$
   として $`\mathrm{Tr}(\mathcal L(\rho)A)=\mathrm{Tr}(\rho\,\mathcal L^*(A))`$ であること、
   $`\mathcal L^*(1)=0`$ であること、およびセクター間のジャンプがなければ $`\mathcal L^*(P_c)=0`$ であることを示せ。
2. （例外点）8.4節のセクター $`N=1`$ では
   $`H_{\mathrm{eff}}=\begin{pmatrix}-i\gamma_1/2&J\\J&-i\gamma_2/2\end{pmatrix}`$ である。その固有値を求め、
   $`J=|\gamma_1-\gamma_2|/4>0`$ のとき対角化可能でないことを示せ。
3. $`m`$ 個のモードにある $`N_0`$ 個のボース粒子がそれぞれレート $`\gamma`$ で失われるとする（ジャンプ $`a_1,\dots,a_m`$ を
   $`N\to N-1`$ に制限したもの）。$`\mathcal M`$ に属する任意のハミルトニアンに対し、$`q=e^{-\gamma t}`$ として
   $`p_N(t)=\binom{N_0}Nq^N(1-q)^{N_0-N}`$ であることを示せ。
4. （純粋位相緩和）$`M_2`$ 上で $`H=\omega|e\rangle\langle e|`$、$`L=\sigma_z=\mathrm{diag}(1,-1)`$、レート
   $`\gamma`$ とする。占有数は一定であり、$`\rho_{ge}(t)=e^{(i\omega-2\gamma)t}\rho_{ge}(0)`$ であること、そして
   $`C=1`$ であるにもかかわらず定常状態が対角行列からなる二次元の空間を張ることを示せ。

<details>
<summary>解答</summary>

1. 巡回性により $`\mathrm{Tr}(L\rho L^*A)=\mathrm{Tr}(\rho L^*AL)`$、
   $`\mathrm{Tr}(\{L^*L,\rho\}A)=\mathrm{Tr}(\rho\{L^*L,A\})`$、
   $`\mathrm{Tr}(-i[H,\rho]A)=\mathrm{Tr}(\rho\,i[H,A])`$ である。セクター間のジャンプについては、
   $`\mathrm{Tr}(\Phi(\rho)A)=\mathrm{Tr}(\rho\Phi^*(A))`$ が $`\Phi^*`$ の定義である
   （6.3節）。すると $`\mathcal L^*(1)=\sum_jg_j(L_j^*L_j-L_j^*L_j)+\sum_mg_m(\Phi_m^*(1)-\Phi_m^*(1))=0`$ であり、
   これは再びトレース保存にほかならない。セクター間のジャンプがなければ、$`P_c`$ は $`H`$ および
   すべての $`L_j`$ と可換であるから、$`L_j^*P_cL_j=L_j^*L_jP_c=\frac12\{L_j^*L_j,P_c\}`$ となり、各 $`\Phi_m`$ はブロック $`K^{cc}_{m,l}`$ しか持たないから
   $`\Phi_m^*(P_c)=\Phi_m^*(1)P_c=\frac12\{\Phi_m^*(1),P_c\}`$ となって、$`\mathcal L^*(P_c)=0`$ である。したがって
   $`p_c'=\mathrm{Tr}(\rho\mathcal L^*(P_c))=0`$ であり、これは命題 8.9(1) の第二の証明である。

2. $`\delta=(\gamma_1-\gamma_2)/4`$ として
   $`H_{\mathrm{eff}}=-\frac i4(\gamma_1+\gamma_2)1+\begin{pmatrix}-i\delta&J\\J&i\delta\end{pmatrix}`$ である。第二の行列はトレース $`0`$、
   行列式 $`\delta^2-J^2`$ をもつので、固有値は $`-\frac i4(\gamma_1+\gamma_2)\pm\sqrt{J^2-\delta^2}`$ である：
   $`J<|\delta|`$ ならば二つの減衰レート、$`J>|\delta|`$ ならば一つのレートと二つの振動数である。$`J=|\delta|>0`$ のとき、第二の行列 $`Y`$ は
   $`Y^2=(J^2-\delta^2)1=0`$ を満たすが $`Y\ne0`$ であるから、$`H_{\mathrm{eff}}`$ はただ一つの固有値
   $`-\frac i4(\gamma_1+\gamma_2)`$ をもつが $`1`$ の定数倍ではない：したがって対角化可能でない。

3. セクター $`N`$ 上で $`\sum_ja_j^*a_j=N1`$ であり、セクターの組は $`(N-1,N)`$ だけであるから、命題
   8.9(2) により $`p_N'=-\gamma Np_N+\gamma(N+1)p_{N+1}`$ である。$`q'=-\gamma q`$ を用いると、二項分布の
   項の導関数は $`-\gamma Np_N+\gamma(N_0-N)\binom{N_0}Nq^{N+1}(1-q)^{N_0-N-1}`$ であり、
   $`(N+1)\binom{N_0}{N+1}=(N_0-N)\binom{N_0}N`$ により第二項は $`\gamma(N+1)p_{N+1}`$ である。$`t=0`$ では
   これは $`N_0`$ に集中した点質量であり、一意性は補題 8.1 による。

4. $`L\rho L^*-\frac12\{L^*L,\rho\}=\sigma_z\rho\sigma_z-\rho`$ は非対角成分に $`-2`$ を掛け、
   対角成分を消す。また $`-i[H,\rho]_{ge}=i\omega\rho_{ge}`$ である。したがって
   （$`\gamma>0`$ のとき）$`\mathcal L(\sigma)=0`$ であるための必要十分条件は $`\sigma_{ge}=0`$ であり、スペクトルは $`\{0,0,-2\gamma\pm i\omega\}`$ である。$`U_t=e^{-itH}`$ とすると
   $`e^{t\mathcal L}(\rho)=\frac{1+e^{-2\gamma t}}2U_t\rho U_t^*+\frac{1-e^{-2\gamma t}}2U_t\sigma_z\rho\sigma_zU_t^*`$ である：
   $`\sigma_z`$ は $`U_t`$ と可換で $`\sigma_z^2=1`$ であるから、微分すると右辺 $`R(t)`$ が
   $`R'=-i[H,R]+\gamma(\sigma_zR\sigma_z-R)=\mathcal L(R)`$ と $`R(0)=\rho`$ を満たすことがわかり、補題 8.1 が適用できる。$`t>0`$ に対し、これは単位的チャネルであり、その不動点は
   クラウス作用素と可換な行列（定理 6.10）、すなわち対角行列である。

</details>

*計算：* 演習問題 2 と 4 の数値的検証 — [チュートリアル、演習問題](../../tutorial/ja/08_dynamics.md#x-8)。

**次へ：**[第9章](09_spin_chains.md)ではスピン鎖を扱う：磁化と運動量のセクター、
大きなセクターに対するランチョス法とクリロフ法による時間発展、そして不規則性。
