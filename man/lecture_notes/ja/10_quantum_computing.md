[English](../10_quantum_computing.md) · [Español](../es/10_quantum_computing.md) · [Français](../fr/10_quantum_computing.md) · [Deutsch](../de/10_quantum_computing.md) · [中文](../zh/10_quantum_computing.md) · **日本語** · [Русский](../ru/10_quantum_computing.md)

# 10. 量子計算

量子計算機上のノイズはチャネルであり、チャネルは一つの環を生成する：そのクラウス作用素が生成する $`*`$-環である。
その構造（第2章）は、どの情報がノイズを受けても無傷で残るかを教えてくれる。また、別の環が関わる、クラウス作用素の積に
ついての条件は、回復操作によってどの情報を復元できるかを教えてくれる。最後の部分では電荷を保存するアルゴリズムを扱う：
その状態は電荷セクターの環の中にあり、振幅減衰は状態をセクター間で移動させる。

**本章で学ぶこと**：単位的チャネルとその双対が同じ不動点（ノイズの可換子環）をもつ理由、完全な証明付きの
クニル–ラフラム（Knill–Laflamme）の定理、セクターのラベルだけを保護する方法、そして量子交互作用素アンザッツの
XY ミキサー、振幅減衰、事後選択がハミング重みのセクターの環 $`\bigoplus_wM_{\binom nw}`$ の中でどのように見えるか。本章の数値的な検算は
[チュートリアル 10](../../tutorial/ja/10_quantum_computing.md)にある。

本章を通じて、$`M_D(\mathbb C)`$ 上の **チャネル** とは写像 $`\Phi(\rho)=\sum_{i=1}^rE_i\rho E_i^*`$ であって、その
**クラウス作用素** が $`\sum_iE_i^*E_i=1`$（トレース保存）を満たすものをいう。さらに
$`\sum_iE_iE_i^*=1`$ であるとき **単位的** であるという。その **双対** $`\Phi^*(T)=\sum_iE_i^*TE_i`$ は
$`\mathrm{Tr}(\Phi(\rho)T)=\mathrm{Tr}(\rho\,\Phi^*(T))`$ を満たす（6.2節）。（定理 6.2 により、任意の完全
正値写像はこの形をもつ。ここではチャネルはクラウス作用素によって *与えられる* ものとする。）$`\Phi`$ の
**相互作用環** $`\mathcal A`$ とは、$`1`$ と $`E_i`$ たちが生成する $`*`$-環、すなわち $`E_i`$ と $`E_i^*`$ のすべての積の
線形包である。行列が $`\mathcal A`$ と可換であるための必要十分条件は、それがすべての
$`E_i`$ および $`E_i^*`$ と可換であることである。射影 $`P_K`$ をもつ部分空間 $`K`$ について $`\rho=P_K\rho P_K`$ であるとき、状態 $`\rho`$ は $`K`$ に **台をもつ** という。
$`n`$ 量子ビット上では $`D=2^n`$ であり、$`X,Y,Z`$ はパウリ行列、$`Z=\mathrm{diag}(1,-1)`$、$`X_j`$ は量子ビット $`j`$ に作用する $`X`$ である（番号は $`1,\dots,n`$。
10.4節はライブラリに従い、サイトを $`0,\dots,n-1`$ とする）。$`|b\rangle`$（$`b\in\{0,1\}^n`$）は計算
基底、$`|b|`$ は **ハミング重み**、$`P_w`$ は $`H_w=\mathrm{span}\{|b\rangle:|b|=w\}`$ の上への射影である。

<a id="s-10-1"></a>

## 10.1 ノイズ環、ノイズレス部分系、不動点

定理 2.6 により、$`V^*\mathcal AV=\bigoplus_cM_{k_c}\otimes1_{m_c}`$ かつ
$`V^*\mathcal A'V=\bigoplus_c1_{k_c}\otimes M_{m_c}`$ となるユニタリ $`V`$ が存在する：各クラウス作用素は、各直和成分の第一テンソル
因子にのみ作用する。

<a id="r-10-1"></a>

**命題 10.1（ノイズレス部分系；デコヒーレンスフリー部分空間）.** (a)
$`V^*E_iV=\bigoplus_cE_{i,c}\otimes1_{m_c}`$ と書き、$`\mathbb C^{k_c}`$ 上の状態 $`\sigma`$ と $`\mathbb C^{m_c}`$ 上の状態 $`\tau`$ に対し、直和成分
$`c`$ に置かれた $`\rho=V(\sigma\otimes\tau)V^*`$ を考える。このとき
$`M_{k_c}`$ 上のチャネル $`\Phi_c(\sigma)=\sum_iE_{i,c}\sigma E_{i,c}^*`$ により
$`\Phi(\rho)=V(\Phi_c(\sigma)\otimes\tau)V^*`$ となる：状態 $`\tau`$ はまったく影響を受けない。さらに、すべての $`T\in\mathcal A'`$ に対し $`\Phi^*(T)=T`$ である。
(b) 部分空間 $`K`$ がすべての $`i`$ に対し $`E_iP_K=e_iP_K`$（$`e_i\in\mathbb C`$）を満たすならば、$`K`$ に台をもつ
すべての状態 $`\rho`$ に対し $`\Phi(\rho)=\rho`$ である。

*証明.* (a) $`\Phi(\rho)=V\sum_i(E_{i,c}\sigma E_{i,c}^*\otimes\tau)V^*`$ であり、$`\Phi_c`$ はトレース保存である。
なぜなら $`V^*(\sum_iE_i^*E_i)V=1`$ の直和成分 $`c`$ は $`(\sum_iE_{i,c}^*E_{i,c})\otimes1_{m_c}`$ だからである。
$`T\in\mathcal A'`$ に対しては $`\Phi^*(T)=T\sum_iE_i^*E_i=T`$ である。(b) $`P_KE_i^*=(E_iP_K)^*=\bar e_iP_K`$ であるから
$`E_i\rho E_i^*=E_iP_K\rho P_KE_i^*=|e_i|^2\rho`$ であり、$`\sum_i|e_i|^2P_K=P_K\sum_iE_i^*E_iP_K=P_K`$ である。$`\square`$

$`k_c=1`$ の直和成分は場合 (b) であり、**デコヒーレンスフリー部分空間** と呼ばれる（Zanardi–Rasetti 1997、
Lidar–Chuang–Whaley 1998）。$`k_c>1`$ のとき、因子 $`\mathbb C^{m_c}`$ は **ノイズレス部分系** である
（Knill–Laflamme–Viola 2000）：ノイズは作用するが、もう一方の因子にのみ作用する。命題 10.1 は
任意のチャネルについて成り立つ。単位的チャネルについては、保存されるオブザーバブルは *ちょうど* $`\mathcal A'`$ のものである：

<a id="r-10-2"></a>

**定理 10.2（単位的チャネルの不動点）.** $`\Phi`$ を単位的とする。このとき
```math
\mathrm{Fix}(\Phi)=\{T\in M_D(\mathbb C):\Phi(T)=T\}=\mathcal A'=\mathrm{Fix}(\Phi^*).
```
特に $`\mathrm{Fix}(\Phi)`$ はフォン・ノイマン環であり、$`\bigoplus_cM_{m_c}`$ と $`*`$-同型である。

*証明.* 最初の等式は、環 $`M_D(\mathbb C)`$（ブロックは一つ）に対する定理 6.10 である：$`T`$ が不動
点であるための必要十分条件は、$`T`$ がすべての $`E_i`$ および $`E_i^*`$ と可換であること、すなわち $`T\in\mathcal A'`$ である。双対 $`\Phi^*`$ もまた
単位的チャネルである：そのクラウス作用素 $`E_i^*`$ は $`\sum_iE_iE_i^*=1`$（トレース保存）および
$`\sum_iE_i^*E_i=1`$（単位的）を満たし、同じ $`*`$-環 $`\mathcal A`$ を生成する。したがって定理 6.10 を
$`\Phi^*`$ に適用すると $`\mathrm{Fix}(\Phi^*)=\mathcal A'`$ を得る。最後の主張は定理 2.6 である。$`\square`$

*計算：* 集団的位相緩和に対する $`\mathrm{Fix}(\Phi)=\mathcal A'`$。可換子環と不動点は零空間として計算する — [チュートリアル、10.1節](../../tutorial/ja/10_quantum_computing.md#s-10-1)。

任意のチャネルについて $`\mathcal A'\subseteq\mathrm{Fix}(\Phi^*)`$ である（命題 10.1(a)）が、
単位性がなければ定理は成り立たない：振幅減衰については $`\mathrm{Fix}(\Phi)=\mathbb C|0\rangle\langle0|`$ である一方、
$`\mathcal A'=\mathrm{Fix}(\Phi^*)=\mathbb C1`$ である（6.4節と演習問題 1）。（歴史的には：Kribs 2003、
Arias–Gheondea–Gudder 2002。）

<a id="r-10-3"></a>

**系 10.3（集団的位相緩和）.** $`U=e^{-i\theta Z_{\mathrm{tot}}}`$、
$`Z_{\mathrm{tot}}=\sum_jZ_j`$ とし、$`d=1,\dots,n`$ に対し $`e^{2i\theta d}\neq1`$ であるとし、
$`\Phi(\rho)=\tfrac12(U\rho U^*+U^*\rho U)`$ とする。このとき $`\mathcal A=\mathrm{span}\{P_0,\dots,P_n\}`$ であり、すべての $`H_w`$
はデコヒーレンスフリー部分空間であり、$`\mathrm{Fix}(\Phi)=\bigoplus_wB(H_w)\cong\bigoplus_wM_{\binom nw}`$ で、
その次元は $`\sum_w\binom nw^2`$ である。

*証明.* $`Z_{\mathrm{tot}}|b\rangle=(n-2|b|)|b\rangle`$ であるから $`U=\sum_we^{-i\theta(n-2w)}P_w`$ であり、その
固有値は相異なる（それらの比は $`e^{2i\theta(w-w')}`$、$`0<|w-w'|\le n`$ である）。定理 1.5(a) の証明における
ラグランジュ多項式により $`P_w=p_w(U)\in\mathcal A`$ である。逆に $`U,U^*`$ は環
$`\mathrm{span}\{P_w\}`$ に属する。したがって $`\mathcal A=\mathrm{span}\{P_w\}`$（$`k_c=1`$、$`m_c=\binom nw`$）であり、命題
10.1(b) が各 $`H_w`$ に適用できる。$`T`$ がすべての $`P_w`$ と可換ならば、
$`v\neq v'`$ に対し $`P_vTP_{v'}=TP_vP_{v'}=0`$、すなわち $`T`$ はブロック対角であり、逆も成り立つ。残りは定理 10.2 による。$`\square`$

*計算：* $`n=4`$ に対する次元 $`5`$、$`70`$、$`70`$ — [チュートリアル、10.1節](../../tutorial/ja/10_quantum_computing.md#s-10-1)。*ライブラリでは：* $`\bigoplus_wM_{\binom nw}`$ は `SpinChain(n).algebra` であり、このチャネルはその上で恒等写像である — [チュートリアル、10.1節](../../tutorial/ja/10_quantum_computing.md#s-10-1)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. すべての量子ビットで同じである揺らぐ場は、$`|b\rangle`$ の位相を
> $`n-2|b|`$ に比例してずらす。同じハミング重みの状態は共通の位相を獲得するので、$`|01\rangle,|10\rangle`$
> はこのノイズの影響を受けない量子ビットを符号化する（イオントラップ：Kielpinski et al., *Science* 291, 2001）。一方、
> 異なる重みの間のコヒーレンスは減衰する：ノイズが超選択則を課すのである（3.3節）。

**集団的デコヒーレンス.** チャネルのクラウス作用素が回転 $`e^{-i\theta\,\mathbf n\cdot\mathbf S}`$
（$`\theta\in\mathbb R`$、$`\mathbf n\in\mathbb R^3`$ は単位ベクトル）の定数倍であるとき、そのチャネルを **集団的** と呼ぶ。ここで
$`S_\alpha=\tfrac12\sum_j\sigma_{\alpha,j}`$、$`\sigma\in\{X,Y,Z\}`$ である。これは向きも揺らぐ場のモデルである（モデルであって、定理ではない）。
$`\mathbf n\cdot\mathbf S`$ は自己共役であるから、$`e^{-i\theta\,\mathbf n\cdot\mathbf S}`$ は $`\mathbf n\cdot\mathbf S`$ の多項式であり
（定理 1.5(a) の証明）、したがってこれらのクラウス作用素は $`1,S_x,S_y,S_z`$ が生成する $`*`$-環 $`\mathcal B`$ に属する。

<a id="r-10-4"></a>

**命題 10.4（一重項部分空間）.** $`K=\ker S_x\cap\ker S_y\cap\ker S_z`$ とする。クラウス作用素が
$`\mathcal B`$ に属する任意のチャネルは、$`K`$ に台をもつすべての状態を固定する。

*証明.* $`\mathcal B`$ は $`S_x,S_y,S_z`$ の語（および $`1`$）で張られる。長さ $`\ge1`$ の語は $`K`$ を
零に写すので、$`\mathcal B`$ の各元は $`K`$ 上でスカラーとして作用する。命題 10.1(b) を適用すればよい。$`\square`$

$`n=4`$ に対し、命題 4.13（[第4章](04_symmetries.md)、$`L=4`$ の場合）により
$`\mathcal B\cong M_5\oplus M_3\oplus M_1`$ であり、重複度は $`1,3,2`$（全スピン $`2,1,0`$）である。よって $`\dim\mathcal B=35`$ であり、
$`\mathcal B'\cong M_1\oplus M_3\oplus M_2`$ の次元は $`14`$ である。さらに $`\dim K=m_0=2`$ である：$`S_a`$ は直和成分ごとに作用するので、
$`K`$ は各直和成分との共通部分の和である。スピン $`S`$ の直和成分上では、それらは $`S^{(S)}_a\otimes1`$ として作用し（系 4.12）、
$`\sum_bu_b\otimes f_b`$（$`f_b`$ は $`\mathbb C^{m_S}`$ の基底）がそれらすべてによって零に写されるための必要十分条件は、各 $`u_b`$ が
すべての $`S^{(S)}_a`$ によって零に写されることである。そのような $`u`$ は既約なスピン $`S`$ の不変部分空間をなし（定理 4.11）、
$`S>0`$ ならば全体ではない（$`S^{(S)}_z\neq0`$ だから）ので $`0`$ であり、$`S=0`$ ならば $`\mathbb C^1`$ である。
命題 10.1(a) とその証明は、クラウス作用素を含む $`V(\bigoplus_cM_{k_c}\otimes1_{m_c})V^*`$ の形の任意の $`*`$-環に対して
そのまま成り立つ。ここではそれは $`\mathcal B`$ である：クラウス作用素が $`\mathcal B`$ に属するすべてのチャネルに対し、四つの
量子ビットは、デコヒーレンスフリーな量子ビット（$`k=1`$、$`m=2`$）とノイズレスなクォートリット（$`k=3`$、$`m=3`$）を担う。

*計算：* $`\dim\mathcal B'=14`$、$`\dim\mathcal B=35`$、$`\dim K=2`$、およびランダムな集団回転に対する命題 10.4 — [チュートリアル、10.1節](../../tutorial/ja/10_quantum_computing.md#s-10-1)。

<a id="s-10-2"></a>

## 10.2 誤り訂正：クニル–ラフラムの定理

**符号** とは、射影 $`P`$ をもつ零でない部分空間 $`\mathcal C\subseteq\mathbb C^D`$ のことである。$`\mathcal C`$ に台をもつ
すべての状態 $`\rho`$ に対し $`\mathcal R(\Phi(\rho))=\rho`$ となるチャネル $`\mathcal R`$（**回復**）が存在するとき、
$`\mathcal C`$ は $`\Phi`$ を **訂正する** という。

<a id="r-10-5"></a>

**定理 10.5（クニル–ラフラム）.** $`\mathcal C`$ が $`\Phi`$ を訂正するための必要十分条件は、次を満たす数 $`\lambda_{ij}`$
が存在することである。
```math
PE_i^*E_jP=\lambda_{ij}P\qquad\text{for all }i,j.
```

*証明.* *必要性.* $`\mathcal R`$ のクラウス作用素を $`R_m`$ とする。すると $`\mathcal R\circ\Phi`$ の
クラウス作用素は $`R_mE_i`$ である。単位ベクトル $`\psi\in\mathcal C`$ に対し
$`\sum_{m,i}R_mE_i|\psi\rangle\langle\psi|E_i^*R_m^*=|\psi\rangle\langle\psi|`$ であり、$`\phi\perp\psi`$ に対し
$`\phi`$ での期待値をとると $`\sum_{m,i}|\langle\phi,R_mE_i\psi\rangle|^2=0`$ となるから、$`R_mE_i\psi\in\mathbb C\psi`$ である。
したがって $`R_mE_iP=c_{mi}P`$ である：すべての $`x`$ に対し $`Tx\in\mathbb Cx`$ となる線形写像 $`T`$ は $`1`$ の定数倍である（
$`Tx=ax`$、$`Ty=by`$、$`T(x+y)=c(x+y)`$ で $`x,y`$ が一次独立ならば $`a=c=b`$）。$`\sum_mR_m^*R_m=1`$ であるから
$`PE_i^*E_jP=\sum_m(R_mE_iP)^*(R_mE_jP)=\sum_m\bar c_{mi}c_{mj}P`$ である。

*十分性.* $`\lambda_{ji}P=(PE_i^*E_jP)^*=\bar\lambda_{ij}P`$ であるから $`\Lambda=(\lambda_{ij})`$ はエルミートであり、
$`F=\sum_jc_jE_j`$ に対し $`\sum_{ij}\bar c_ic_j\lambda_{ij}P=PF^*FP\ge0`$ であるから正値である。$`U`$ をユニタリ、
$`D=\mathrm{diag}(d_k)`$、$`d_k\ge0`$ として $`\Lambda=UDU^*`$ と書き（定理 1.3）、$`F_k=\sum_iu_{ik}E_i`$ とおく。すると
$`\sum_kF_k\rho F_k^*=\sum_{ij}(UU^*)_{ij}E_i\rho E_j^*=\Phi(\rho)`$ かつ $`PF_k^*F_lP=(U^*\Lambda U)_{kl}P=d_k\delta_{kl}P`$ である。
$`d_k=0`$ ならば、すべての $`x`$ に対し $`\|F_kPx\|^2=0`$ であるから $`F_kP=0`$ である。$`d_k>0`$ ならば $`V_k=F_kP/\sqrt{d_k}`$ とおく。すると
$`V_k^*V_l=\delta_{kl}P`$ であるから $`V_kP=V_k`$ であり、$`Q_k=V_kV_k^*`$ は互いに直交する射影で、
$`Q=\sum_kQ_k`$ は射影である。回復
```math
\mathcal R(\sigma)=\sum_{k:d_k>0}V_k^*\sigma V_k+(1-Q)\sigma(1-Q)
```
は $`\sum_kV_kV_k^*+(1-Q)^2=1`$ であるからチャネルである。$`\rho=P\rho P`$ に対し $`\Phi(\rho)=\sum_{d_k>0}d_kV_k\rho V_k^*`$ であり、
$`V_l^*V_k=\delta_{kl}P`$ かつ $`(1-Q)V_k=(1-Q)Q_kV_k=0`$ であるから $`\mathcal R(\Phi(\rho))=\sum_kd_k\rho`$ となる。そして
$`\sum_i\lambda_{ii}P=P\sum_iE_i^*E_iP=P`$ であるから $`\sum_kd_k=\mathrm{Tr}\,\Lambda=1`$ である。$`\square`$

（Knill–Laflamme 1997；Nielsen–Chuang、定理 10.1。）$`\mathcal R`$ は **シンドローム**——互いに直交する
部分空間 $`Q_k\mathbb C^D`$ のどれが状態を含むか——を測定し、状態を回転させて符号に戻す。

**例：ビット反転符号.** $`\mathcal C=\mathrm{span}\{|000\rangle,|111\rangle\}`$ とし、ノイズを「高々一つのビット
反転」、$`E_0=\sqrt{p_0}\,1`$、$`E_j=\sqrt{p_j}\,X_j`$ とする。$`i\neq j`$ に対し、$`X_iX_j`$（$`X_0=1`$ とする）は一つまたは二つの
ビットを反転し、$`\mathcal C`$ を $`\mathcal C^\perp`$ に写す。したがって $`\Lambda=\mathrm{diag}(p_0,\dots,p_3)`$、$`V_k=X_kP`$、
$`Q=1`$、$`\mathcal R(\sigma)=\sum_{k=0}^3PX_k\sigma X_kP`$ である。各量子ビットが確率 $`p`$ で独立に
反転する場合には条件は成り立たない（クラウス作用素には $`X_1X_2`$ と $`X_3`$ の定数倍が含まれ、
$`PX_1X_2X_3P=|000\rangle\langle111|+|111\rangle\langle000|`$ は $`P`$ の定数倍ではない）。そして $`\mathcal R`$ は $`\ge2`$ 個の反転が起きたときには $`|0_L\rangle=|000\rangle`$ を $`|1_L\rangle`$ に写し、
そうでなければそのままにする：回復後の $`|0_L\rangle`$ の忠実度は
高々一つの反転が起きる確率であり、$`p=0.1`$ に対して $`(1-p)^3+3p(1-p)^2=1-3p^2+2p^3=0.972`$ である。

*計算：* この例での定理 10.5：行列 $`\Lambda`$、回復 $`\mathcal R`$、忠実度 $`0.972`$ — [チュートリアル、10.2節](../../tutorial/ja/10_quantum_computing.md#s-10-2)。*ライブラリでは：* クラウス作用素で与えられるチャネル `channels.kraus_channel` とその合成 `R @ Phi` — [チュートリアル、10.2節](../../tutorial/ja/10_quantum_computing.md#s-10-2)。

<a id="s-10-3"></a>

## 10.3 作用素環による誤り訂正

ハイゼンベルク描像を用いると、符号全体よりも少ないものだけを保護することができる。$`\mathfrak A`$ を
$`\mathcal C=P\mathbb C^D`$ 上のフォン・ノイマン環（$`P`$ を含む、$`\mathcal C`$ 上の作用素の $`*`$-環）とし、
これを保護すべきオブザーバブルとする。次を満たすチャネル
$`\mathcal R`$ が存在するとき、$`\mathfrak A`$ は $`\Phi`$ に対して **訂正可能** であるという。
```math
P\,(\mathcal R\circ\Phi)^*(A)\,P=A\qquad\text{for all }A\in\mathfrak A,
```
すなわち、$`\mathcal C`$ に台をもつすべての状態 $`\rho`$
とすべての $`A\in\mathfrak A`$ に対し $`\mathrm{Tr}(\mathcal R(\Phi(\rho))A)=\mathrm{Tr}(\rho A)`$ である（左辺は
$`\mathrm{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(A)P)`$ であり、これらすべての $`\rho`$ に対し $`\mathrm{Tr}(\rho Y)=0`$ を満たす作用素
$`Y\in B(\mathcal C)`$ は零である：すべての $`x`$ に対し $`\langle x,Yx\rangle=\langle Px,YPx\rangle=0`$ であり、系 1.4 の証明と同様に
偏極を用いればよい）。$`\mathfrak A=B(\mathcal C)`$ に対しては、これは10.2節の意味での訂正である：
$`\sigma=\mathcal R(\Phi(\rho))`$ が $`P\sigma P=\rho`$ を満たすならば（すべての $`A\in B(\mathcal C)`$ をとる）、$`\mathcal R\circ\Phi`$ はトレース保存であるから
$`\mathrm{Tr}(1-P)\sigma(1-P)=\mathrm{Tr}\,\sigma-\mathrm{Tr}\,P\sigma P=1-1=0`$ であり、したがって $`\sigma^{1/2}(1-P)=0`$ であり、
$`\sigma=\rho`$ となる。逆に、$`\mathcal R(\Phi(\rho))=\rho`$ ならばトレースの等式が成り立つ。反対の極端は可換環である：古典的なラベルだけが保護される。

<a id="r-10-6"></a>

**命題 10.6（セクターのラベルの訂正）.** 互いに直交する
射影 $`P_c\neq0`$ により $`P=\sum_{c=1}^CP_c`$ とし、$`\mathfrak A=\mathrm{span}\{P_1,\dots,P_C\}`$ とする。次は同値である：
1. $`N_1,\dots,N_C\ge0`$、$`\sum_cN_c=1`$ であって、すべての $`c`$ と
   $`P_c\mathbb C^D`$ に台をもつすべての状態 $`\rho`$ に対し $`\mathrm{Tr}(N_c\Phi(\rho))=1`$ となるものが存在する；
2. すべての $`i,j`$ と $`c\neq d`$ に対し $`P_cE_i^*E_jP_d=0`$、すなわち $`PE_i^*E_jP\in\mathfrak A'`$ である。ここで $`\mathfrak A'`$ は
   $`\mathfrak A`$ の $`B(\mathcal C)`$ における可換子環である；
3. $`\mathfrak A`$ は $`\Phi`$ に対して訂正可能である。

*証明.* $`T\in B(\mathcal C)`$ がすべての $`P_c`$ と可換であるための必要十分条件は、$`c\neq d`$ に対し $`P_cTP_d=0`$ となることである（系 10.3 と同様）。

(2)$`\Rightarrow`$(1)。部分空間 $`S_c=\mathrm{span}\{E_iP_cx\}`$ は互いに直交する。実際
$`\langle E_iP_cx,E_jP_dy\rangle=\langle x,P_cE_i^*E_jP_dy\rangle=0`$ である。$`c<C`$ に対し $`N_c`$ を $`S_c`$ の上への射影とし、
$`N_C=1-\sum_{c<C}N_c`$ とおくと、これも射影である。$`P_c\mathbb C^D`$ に台をもつ $`\rho`$ に対し、
$`\Phi(\rho)=\sum_iE_iP_c\rho P_cE_i^*`$ の値域は $`S_c`$ に含まれ、$`N_c`$ は $`S_c`$ 上で恒等写像である（$`c=C`$ については
$`S_C\perp S_{c'}`$、$`c'<C`$ であることによる）から、$`\mathrm{Tr}(N_c\Phi(\rho))=\mathrm{Tr}\,\Phi(\rho)=1`$ である。
(1)$`\Rightarrow`$(2)。単位ベクトル $`x\in P_c\mathbb C^D`$ と $`d\neq c`$ に対し、数
$`\mathrm{Tr}(N_{d'}\Phi(|x\rangle\langle x|))\ge0`$ の和は $`1`$ であり、$`d'=c`$ のものが $`1`$ であるから、
$`0=\mathrm{Tr}(N_d\Phi(|x\rangle\langle x|))=\sum_i\|N_d^{1/2}E_ix\|^2`$ である。よって $`N_dE_iP_c=0`$ であり、
$`d`$ について和をとると $`N_cE_iP_c=E_iP_c`$ である。ゆえに $`c\neq d`$ に対し $`P_cE_i^*E_jP_d=P_cE_i^*N_cE_jP_d=0`$ である。
(1),(2)$`\Rightarrow`$(3)。$`\omega_c=P_c/\mathrm{Tr}\,P_c`$ とし、$`\mathcal R(\sigma)=\sum_c\mathrm{Tr}(N_c\sigma)\,\omega_c`$
（「測定してから準備する」）とおく。これは、$`\omega_c`$ の固有対 $`(\mu,u)`$ と正規直交基底のベクトル $`f`$ にわたる
クラウス作用素 $`\sqrt\mu\,|u\rangle\langle f|N_c^{1/2}`$ をもつチャネルである。
$`\mathrm{Tr}(\mathcal R(\sigma)P_d)=\mathrm{Tr}(N_d\sigma)`$ であるから $`\mathcal R^*(P_d)=N_d`$ である。$`N_dE_iP_c=\delta_{cd}E_iP_c`$
（前段より）を用い、次に $`\sum_iE_i^*E_i=1`$ を用いると
```math
P(\mathcal R\circ\Phi)^*(P_d)P=P\,\Phi^*(N_d)\,P=\sum_i\sum_{c,c'}P_{c'}E_i^*N_dE_iP_c=\sum_i\sum_{c'}P_{c'}E_i^*E_iP_d=P_d\Big(\sum_iE_i^*E_i\Big)P_d=P_d .
```
(3)$`\Rightarrow`$(1)。$`c<C`$ に対し $`N_c=\mathcal R^*(P_c)`$、また $`N_C=\mathcal R^*(1-\sum_{c<C}P_c)`$ とおく。これらは
正値であり、和は $`\mathcal R^*(1)=1`$ である。$`P_c\mathbb C^D`$ に台をもつ $`\rho`$ と $`d<C`$ に対し
$`\mathrm{Tr}(N_d\Phi(\rho))=\mathrm{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(P_d)P)=\mathrm{Tr}(\rho P_d)=\delta_{cd}`$ であり、したがって
$`\mathrm{Tr}(N_C\Phi(\rho))=\delta_{cC}`$ でもある。$`\square`$

<a id="r-10-7"></a>

**定理 10.7（Bény–Kempf–Kribs）.** $`\mathcal C`$ 上のフォン・ノイマン環 $`\mathfrak A`$
が $`\Phi`$ に対して訂正可能であるための必要十分条件は、すべての $`i,j`$ に対し $`PE_i^*E_jP\in\mathfrak A'`$ となることである。ここで
$`\mathfrak A'`$ は $`\mathfrak A`$ の $`B(\mathcal C)`$ における可換子環である。

*証明.* これは[付録 D](D_gksl.md)で証明される定理 D.9 である。その証明は回復チャネルも構成する：それは
$`\mathfrak A\cong\bigoplus_cM_{k_c}`$ のセクターとシンドロームを測定し、定理 10.5 の証明と同様に各因子 $`\mathbb C^{k_c}`$ 上の
誤りを元に戻す。$`\square`$

$`\mathfrak A=B(\mathcal C)`$ に対しては $`\mathfrak A'=\mathbb CP`$ であり、これは定理 10.5 である。命題 10.6 は
可換な $`\mathfrak A=\mathrm{span}\{P_c\}`$ の場合である（系 D.10(a),(b)）。
$`\mathcal C=\mathbb C^k\otimes\mathbb C^m`$ 上の $`\mathfrak A=M_k\otimes1_m`$ に対しては作用素（部分系）
符号が得られる。10.1節のノイズレス部分系はそのような符号であり、$`\mathcal R=\mathrm{id}`$ によって訂正される
（系 D.10(c)、ただし二つのテンソル因子を入れ替える）。定理 2.6 により、すべての $`\mathfrak A`$ はある
$`\bigoplus_cM_{k_c}`$ と $`*`$-同型である。$`C\ge2`$ で、ある $`k_c\ge2`$ であるとき、訂正可能な $`\mathfrak A`$ はセクターのラベル $`c`$ を
各 $`M_{k_c}`$ の量子情報とともに保護する——**ハイブリッド** 符号である（系 D.10 の後の注意。クニル–ラフラム条件は
成り立たないが定理 10.7 が適用できる例は [チュートリアル、D.5節](../../tutorial/ja/D_gksl.md#s-D-5) で計算する）。
歴史的には：C. Bény, A. Kempf, D. W. Kribs, *Phys. Rev. Lett.* 98, 100502 (2007) および *Phys. Rev. A* 76,
042303 (2007)。作用素符号：Kribs–Laflamme–Poulin, *Phys. Rev. Lett.* 94, 180501 (2005)。

**例.** 独立な位相反転（量子ビットごとの確率 $`q`$）を受け、その後に
高々一つのビット反転を受ける繰り返し符号を考える。クラウス作用素は $`X_kZ^s`$（$`s\in\{0,1\}^3`$）の定数倍であり、各
$`E_i^*E_j`$ は高々二つのビットを反転するので、$`P_0=|000\rangle\langle000|`$、
$`P_1=|111\rangle\langle111|`$ について条件 2 が成り立つ：古典ビットは生き残り、多数決
$`N_0=\sum_{|b|\le1}|b\rangle\langle b|`$、$`N_1=1-N_0`$ で読み出される。量子ビットは生き残らない：$`PZ_1P`$ は $`P`$ の定数倍ではなく、
10.2節の回復 $`\mathcal R`$ の後、状態 $`|{+_L}\rangle=(|000\rangle+|111\rangle)/\sqrt2`$ は、奇数個の位相反転が起きたときに
ちょうど論理的な位相反転 $`Z_1Z_2Z_3`$ を被る。その忠実度は $`q=0.05`$ に対して
$`(1+(1-2q)^3)/2=0.8645`$ である。（チャネルは $`\Phi_{\mathrm{one}}\circ\mathcal Z`$ である。ここで $`\Phi_{\mathrm{one}}`$ は
ビット反転の例のチャネル「高々一つのビット反転」、$`\mathcal Z`$ は独立な位相反転である。
$`\mathcal C`$ 上ではすべての $`Z^s`$ が $`Z_1^{|s|}`$ として作用するので、$`\mathcal Z(\rho)=\sum_sq^{|s|}(1-q)^{3-|s|}Z_1^{|s|}\rho Z_1^{|s|}`$ は
再び $`\mathcal C`$ に台をもち、定理 10.5 により $`\mathcal R\circ\Phi_{\mathrm{one}}`$ はそのような状態の上で恒等写像である。
$`Z_1`$ は $`|{+_L}\rangle`$ を $`|{-_L}\rangle\perp|{+_L}\rangle`$ に写し、$`|s|`$ が偶数である確率は
$`\frac12\big((1-q+q)^3+(1-q-q)^3\big)`$ である。）

*計算：* この例での命題 10.6：条件 2、多数決 $`N_0,N_1`$、$`|{+_L}\rangle`$ の忠実度 $`0.8645`$ — [チュートリアル、10.3節](../../tutorial/ja/10_quantum_computing.md#s-10-3)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. 可換環を保護することは、古典的な情報（一つの
> ビット、超選択セクターのラベル）を保護することを意味し、これは量子ビットを保護するよりずっと安価である：セクターの射影と
> 可換なノイズは無害である（このとき $`c\neq d`$ に対し $`P_cE_i^*E_jP_d=E_i^*E_jP_cP_d=0`$ であり、これは命題 10.6 の条件 2 である）。

<a id="s-10-4"></a>

## 10.4 制約付きミキサー

量子交互作用素アンザッツ（Hadfield et al., *Algorithms* 12, 34, 2019）は
$`|\beta,\gamma\rangle=e^{-i\beta_pH_M}e^{-i\gamma_pH_C}\cdots e^{-i\beta_1H_M}e^{-i\gamma_1H_C}|b_0\rangle`$ を準備する。ここで $`H_C`$ は対角な
コスト・ハミルトニアン、$`H_M`$ は **ミキサー** である。制約 $`|b|=w`$（グラフの二分割、$`w`$ 銘柄の
資産からなるポートフォリオ）に対しては、$`|b_0|=w`$ ととり、リング上の **XY ミキサー** $`H_M=\tfrac14\sum_{(i,j)}(X_iX_j+Y_iY_j)`$ を用いる。

<a id="r-10-8"></a>

**命題 10.8.** $`H_M`$、$`H_C`$、$`e^{-i\beta H_M}`$、$`e^{-i\gamma H_C}`$ は
$`\mathcal M=\{P_0,\dots,P_n\}'=\bigoplus_wB(H_w)\cong\bigoplus_wM_{\binom nw}`$ に属し、したがって $`|\beta,\gamma\rangle\in H_w`$ である：
アンザッツが実行可能部分空間を離れることはない。

*証明.* 量子ビット $`i,j`$ 上で $`X\otimes X+Y\otimes Y=2(|01\rangle\langle10|+|10\rangle\langle01|)`$ である（$`Y|0\rangle=i|1\rangle`$、
$`Y|1\rangle=-i|0\rangle`$ を用いる）：$`X_iX_j+Y_iY_j`$ はビット $`i,j`$ が異なるときそれらを入れ替え、そうでないときは $`|b\rangle`$ を
零に写すので、各 $`H_w`$ をそれ自身に写す。対角な $`H_C`$ は対角な $`P_w`$ と可換であり、
指数関数も同様である（定理 1.5(d)）。$`\square`$

*ライブラリでは：* $`\mathcal M`$ は `SpinChain(n).algebra` であり、XY ミキサーは $`\Delta=0`$ とした `SpinChain.xxz` である — [チュートリアル、10.4節](../../tutorial/ja/10_quantum_computing.md#s-10-4)。*計算：* 6 量子ビット上のグラフ二分割に対する命題 10.8。ライブラリと稠密行列の両方で計算する — [チュートリアル、10.4節](../../tutorial/ja/10_quantum_computing.md#s-10-4)。

量子ビット $`i`$ の **振幅減衰**（確率 $`g`$ で $`|1\rangle\to|0\rangle`$）はクラウス作用素
$`K_0^{(i)}=|0\rangle\langle0|_i+\sqrt{1-g}\,|1\rangle\langle1|_i`$ と $`K_1^{(i)}=\sqrt g\,|0\rangle\langle1|_i`$ をもち、後者は
重みを下げる：これはセクター間の写像である。

*ライブラリでは：* `SpinChain.site_amplitude_damping`（`InterSectorChannel`） — [チュートリアル、10.4節](../../tutorial/ja/10_quantum_computing.md#s-10-4)。

<a id="r-10-9"></a>

**命題 10.9（漏れと事後選択）.** $`\rho_0`$ を $`H_w`$ に台をもつ状態、$`U_l`$ を
$`\mathcal M`$ のユニタリ、$`\mathcal D`$ をすべての量子ビットの減衰とし、$`\rho_l=\mathcal D(U_l\rho_{l-1}U_l^*)`$ とする。
このとき各 $`\rho_l`$ は重み $`\le w`$ のブロック対角であり、
```math
P_w\rho_pP_w=(1-g)^{wp}\;U\rho_0U^*,\qquad U=U_p\cdots U_1 :
```
である。すなわち、重み $`w`$ は確率 $`(1-g)^{wp}`$ で見いだされ、事後選択された状態はノイズのない状態である。

*証明.* $`\mathcal D`$ のクラウス作用素は $`K_s=\prod_iK^{(i)}_{s_i}`$、$`s\in\{0,1\}^n`$ である。$`K_s`$ は
$`H_v`$ を $`H_{v-|s|}`$ に写す、すなわち $`P_vK_s=K_sP_{v+|s|}`$ である（$`v>n`$ のとき $`P_v=0`$）。$`\sigma`$ が重み
$`\le w`$ のブロック対角ならば、$`P_vK_s\sigma K_s^*P_{v'}=K_sP_{v+|s|}\sigma P_{v'+|s|}K_s^*`$ は $`v\neq v'`$ のとき、および
$`v>w-|s|`$ のとき消える：$`\mathcal D(\sigma)`$ は重み $`\le w`$ のブロック対角であり、その重み $`w`$ のブロックは
$`s=0`$ からのみ生じる。$`K_0|b\rangle=(1-g)^{|b|/2}|b\rangle`$ であるから $`P_w\mathcal D(\sigma)P_w=(1-g)^wP_w\sigma P_w`$ である。$`U_l`$
は $`P_v`$ と可換である。$`l`$ に関する帰納法を用いればよい。$`\square`$

*計算：* 二分割回路での命題 10.9。事後選択には `DensityMatrix.condition_on` を用いる — [チュートリアル、10.4節](../../tutorial/ja/10_quantum_computing.md#s-10-4)。

**連続時間.** 同じことがリンドブラッド方程式（定理 8.6、[第8章](08_dynamics.md)）
```math
\dot\rho=-i[H,\rho]+\gamma_1\sum_i\Big(S^-_i\rho S^+_i-\tfrac12\{S^+_iS^-_i,\rho\}\Big)+\gamma_2\sum_k\Big(L_k\rho L_k^*-\tfrac12\{L_k^*L_k,\rho\}\Big)
```
についても成り立つ。ここで $`H\in\mathcal M`$、ジャンプ $`L_k\in\mathcal M`$ はセクター内のもので、例えば純粋位相緩和 $`L_k=S^z_k`$ である。
ここで $`\sum_iS^+_iS^-_i=\sum_i|1\rangle\langle1|_i=N`$ は **個数作用素** であり、$`NP_v=vP_v`$ である：これは中心に属する。

<a id="r-10-10"></a>

**命題 10.10.** $`\rho(0)`$ が $`H_w`$ に台をもつならば $`P_w\rho(t)P_w=e^{-\gamma_1wt}\sigma(t)`$ である。
ここで $`\sigma`$ は $`\gamma_1=0`$ とした同じ方程式の、$`\sigma(0)=\rho(0)`$ を満たす解である。したがって受理確率は
$`e^{-\gamma_1wt}`$ であり、事後選択された状態 $`\sigma(t)`$ は $`\gamma_1`$ に依存しない。

*証明.* $`B_{vv'}=P_v\rho P_{v'}`$ とおく。$`H,L_k,N`$ は $`P_v`$ と可換であり $`P_vS^-_i=S^-_iP_{v+1}`$ であるから、
$`\dot B_{vv'}=-i[H,B_{vv'}]+\gamma_1\sum_iS^-_iB_{v+1,v'+1}S^+_i-\tfrac{\gamma_1}2(v+v')B_{vv'}+\gamma_2\mathcal L_2(B_{vv'})`$
である。ここで $`\mathcal L_2`$ は位相緩和の部分である。$`\max(v,v')>w`$ のブロックは初期値が零の閉じた線形方程式系を
満たすので、それらは消える（補題 8.1 の一意性）。すると $`\dot B_{ww}=-i[H,B_{ww}]-\gamma_1wB_{ww}+\gamma_2\mathcal L_2(B_{ww})`$ であるから、
$`e^{\gamma_1wt}B_{ww}`$ は同じ初期値をもつ $`\sigma`$ の方程式の解であり、$`\sigma`$ に等しい
（補題 8.1）。その方程式はやはり GKSL 形であるから $`\mathrm{Tr}\,\sigma(t)=1`$ である（定理 8.6）。$`\square`$

命題 8.11 と比較せよ。そこでは中心に属さない $`\sum_jL_j^*L_j`$ が許されるが、セクター内のジャンプは許されない。ここでは
損失項 $`\gamma_1N`$ が中心に属するので、ジャンプなしの因子はスカラー $`e^{-\gamma_1wt}`$ となり、セクター
内部の位相緩和が許される。

$`\gamma_2=0`$ ならば $`\sigma(t)=e^{-iHt}\rho(0)e^{iHt}`$ であり、事後選択は完全である。

*計算：* 6 量子ビット上の XY ミキサーに対する命題 10.10。`dynamics.lindblad_evolve` を用いる — [チュートリアル、10.4節](../../tutorial/ja/10_quantum_computing.md#s-10-4)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. エネルギー緩和（$`T_1`$）は超伝導量子ビットの主要な誤りの一つであり、
> より低いセクターへのジャンプである。ハミング重みの測定——中心射影であり、セクター内部では
> 無害である（3.3節）——によってそれは検出され、事後選択によって受理率と引き換えに取り除かれる。
> 位相緩和はセクターの射影と可換であり、検出されない。

<a id="s-10-5"></a>

## 10.5 比較の対象

*ライブラリでは：* 他のソフトウェア（Qiskit、Cirq、PennyLane、QuTiP、stim）との比較 — [チュートリアル、10.5節](../../tutorial/ja/10_quantum_computing.md#s-10-5)。

<a id="x-10"></a>

## 演習問題

1. （定理 10.2 には単位性が必要である。6.4節と比較せよ。）一つの量子ビット上の振幅減衰
   $`E_0=|0\rangle\langle0|+\sqrt{1-g}\,|1\rangle\langle1|`$、$`E_1=\sqrt g\,|0\rangle\langle1|`$、$`0<g<1`$ について、
   $`\mathcal A'=\mathbb C1`$、$`\mathrm{Fix}(\Phi)=\mathbb C|0\rangle\langle0|`$、$`\mathrm{Fix}(\Phi^*)=\mathbb C1`$ であることを示せ。
2. $`\sum_w\binom nw^2=\binom{2n}n`$ を証明せよ。系 10.3 のチャネルで $`\theta=\pi/2`$ としたもの（これは系の
   仮定に反する）について、$`\mathcal A`$ と
   $`\mathrm{Fix}(\Phi)`$ は何か。$`n=4`$ について $`\dim\mathrm{Fix}(\Phi)`$ を確かめよ。
3. （誤り訂正の線形性.）有限個の作用素 $`E_a`$（チャネルのクラウス作用素とは
   限らない）について $`PE_a^*E_bP=\lambda_{ab}P`$ とし、定理 10.5 の証明と同様に $`\Lambda`$ から $`\mathcal R`$ を構成する。
   クラウス作用素が $`\mathrm{span}\{E_a\}`$ に属する任意のチャネルを $`\mathcal R`$ が訂正することを示せ。
   これより、ビット反転符号がコヒーレントな誤り $`e^{-i\phi X_1}`$ を訂正することを導け。

<details>
<summary>解答</summary>

1. $`E_1`$、$`E_1^*`$ とそれらの積からすべての行列単位が得られるので、$`\mathcal A=M_2`$ かつ $`\mathcal A'=\mathbb C1`$ である。
   $`T=(x_{kl})`$ に対し $`\Phi(T)=\begin{pmatrix}x_{00}+gx_{11}&\sqrt{1-g}\,x_{01}\\\sqrt{1-g}\,x_{10}&(1-g)x_{11}\end{pmatrix}=T`$
   となるための必要十分条件は $`x_{11}=x_{01}=x_{10}=0`$ であり、
   $`\Phi^*(T)=\begin{pmatrix}x_{00}&\sqrt{1-g}\,x_{01}\\\sqrt{1-g}\,x_{10}&(1-g)x_{11}+gx_{00}\end{pmatrix}`$
   が $`T`$ に等しいための必要十分条件は $`T=x_{00}1`$ である。（$`\Phi`$ は単位的でなく、$`\Phi^*`$ は単位的であるがトレース保存ではない。）

2. $`\{1,\dots,2n\}`$ の $`n`$ 元部分集合を、前半に含まれる元の個数 $`w`$ によって数えると
   $`\binom{2n}n=\sum_w\binom nw\binom n{n-w}=\sum_w\binom nw^2`$ である。$`\theta=\pi/2`$ のとき
   $`U=(-i)^n\sum_w(-1)^wP_w`$ は $`w`$ の偶奇しか見ないので、$`\mathcal A=\mathrm{span}\{P_{\mathrm{even}},P_{\mathrm{odd}}\}`$
   （$`P_{\mathrm{even}}=\sum_{w\text{ が偶数}}P_w`$）であり、定理 10.2 により $`\mathrm{Fix}(\Phi)`$ は偶奇に関して
   ブロック対角な作用素からなる：$`n=4`$ に対して次元は $`2\cdot(2^{n-1})^2=128`$ である。

3. $`F_k=\sum_au_{ak}E_a`$ とすると $`\mathrm{span}\{E_a\}=\mathrm{span}\{F_k\}`$ であるから、チャネルの各クラウス作用素は
   $`G_m=\sum_kh_{mk}F_k`$ と書け、$`G_mP=\sum_{d_k>0}h_{mk}\sqrt{d_k}V_k`$ である。$`V_l^*V_k=\delta_{lk}P`$ かつ $`(1-Q)V_k=0`$ であるから、
   $`\rho=P\rho P`$ に対し $`\mathcal R(G_m\rho G_m^*)=\sum_l|h_{ml}|^2d_l\,\rho`$ であり、したがって $`\mathcal R(\Phi(\rho))=c\rho`$ である。トレースをとれば
   $`c=1`$ である。$`E_a\in\{1,X_1,X_2,X_3\}`$ に対しては $`\Lambda=1`$、$`V_k=X_kP`$ であり、$`\mathcal R`$ は
   10.2節の回復である。$`e^{-i\phi X_1}=\cos\phi\,1-i\sin\phi\,X_1`$ はこの線形包に属する。

</details>

*計算：* 演習問題 1–3 の数値的検証 — [チュートリアル、演習問題](../../tutorial/ja/10_quantum_computing.md#x-10)。

**次へ：**[第11章](11_optimisation.md)：環のユニタリ群上の最適化と、トレース不等式。
