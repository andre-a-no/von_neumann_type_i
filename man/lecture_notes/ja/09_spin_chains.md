[English](../09_spin_chains.md) · [Español](../es/09_spin_chains.md) · [Français](../fr/09_spin_chains.md) · [Deutsch](../de/09_spin_chains.md) · [中文](../zh/09_spin_chains.md) · **日本語** · [Русский](../ru/09_spin_chains.md)

# 9. スピン鎖

$L$ 個のスピン $\frac12$ からなる鎖は $(\mathbb C^2)^{\otimes L}$ 上に住み、次元 $2^L$ の密な線形代数は
$L\approx14$ を超えると高くつく（密行列は $4^L$ 個の成分をもつ）。対称性はこの限界を押し上げる：磁化は空間を $L+1$ 個のセクターに分け
（第4章）、並進は各セクターをさらに分け、$SU(2)$ はブロックをもう一度簡約する。本章では XXZ 鎖に対して
これらのセクターを構成し、それらが基底状態について何を語るかを証明し、自由フェルミオンとベーテ仮説によって
ライブラリを検証し、さらに密行列（クリロフ法）と乱れのない鎖（乱れ）の先へ進む。

**本章で学ぶこと**：鎖の環がなぜ $\bigoplus_NM_{\binom LN}$ であるか、これらの環の塔がなぜ
パスカルの三角形をブラッテリ図とし、指数 $4\cos^2\frac\pi{2L+4}$ をもつか、運動量セクターの数え方、
マーシャルの符号則と、それが基底状態の運動量とスピンについて導く帰結（リーブ–マティス）、XX 鎖の
ジョルダン–ウィグナーによる解、二マグノンのベーテ方程式、ランチョス法がなぜ単調に収束するか、
そして準位統計と自由フェルミオンのエンタングルメントの測り方。本章の数値的検証は
[チュートリアル 9](../../tutorial/ja/09_spin_chains.md) にある。

本章を通じて、サイトは $0,\dots,L-1$ とし、$(\mathbb C^2)^{\otimes L}$ の基底は配置
$|b\rangle=|b_0\dots b_{L-1}\rangle$ からなり、上向きスピンに対して $b_i=1$ とする（ライブラリの規約：ラベルの
第 $i$ ビットがサイト $i$ である）。$s^a_i$（$a=x,y,z$）はサイト $i$ 上のスピン行列 $\frac12\sigma_a$、
$s^\pm_i=s^x_i\pm is^y_i$、$n_i=s^z_i+\frac12$、$S^a=\sum_is^a_i$ とし、$N=S^z+\frac L2$ は上向きスピンの
個数である。鎖は **開いた** 鎖（ボンド $(i,i+1)$、$i<L-1$）か **周期的** な鎖（ボンド $(L-1,0)$ も含む。このとき
$L\ge3$ で、添字は $L$ を法とする）である。$J>0$（反強磁性）、$\hbar=1$ とし、$\log$ は自然対数である。

<a id="s-9-1"></a>

## 9.1 XXZ 鎖とその $U(1)$ セクター

異方性 $\Delta\in\mathbb R$ と磁場 $h_i\in\mathbb R$ をもつ **XXZ 鎖** は
$$
H=J\sum_{(i,j)}\Big[\tfrac12\big(s^+_is^-_j+s^-_is^+_j\big)+\Delta\,s^z_is^z_j\Big]+\sum_ih_is^z_i ,
$$
である。和はボンドにわたってとる。$\frac12(s^+_is^-_j+s^-_is^+_j)=s^x_is^x_j+s^y_is^y_j$ は上向きスピンをボンドに沿って動かす。
$\Delta=1$ ではハイゼンベルク鎖 $J\sum\mathbf s_i\cdot\mathbf s_j$、$\Delta=0$ では XX 鎖である。

<a id="r-9-1"></a>

**命題 9.1.** $[H,S^z]=0$。したがって $H\in U(1)'\cong\bigoplus_{N=0}^LM_{\binom LN}(\mathbb C)$ であり、
ブロック $N$ は上向きスピンを $N$ 個もつ配置の張る空間に作用する。

*証明.* 対角項は対角な $S^z$ と可換である。$i\neq j$ のとき、$s^+_is^-_j$ は $l\ne i,j$ に対する
$s^z_l$ と可換であり、$[s^z,s^\pm]=\pm s^\pm$（4.5節）により
$[s^z_i+s^z_j,s^+_is^-_j]=[s^z_i,s^+_i]s^-_j+s^+_i[s^z_j,s^-_j]=s^+_is^-_j-s^+_is^-_j=0$ である。
$s^-_is^+_j$ については共役をとればよい。よって $H$ は $N=S^z+\frac L2$ と可換である。$N$ は整数の固有値
$0,\dots,L$ をもち、その固有空間は上向きスピンを $N$ 個もつ配置によって張られ、次元は $\binom LN$ である。残りは
$Q=N$ とした命題 4.10(2) である。$\square$

*計算：* 命題 9.1、および `SpinChain.xxz` が各セクターの配置（`chain.labels[N]`）から構成するブロック $H_N$ を、完全な $2^L\times2^L$ 行列と照合する — [チュートリアル、9.1節](../../tutorial/ja/09_spin_chains.md#s-9-1)。

**鎖の塔.** $\mathcal A_L=\bigoplus_NM_{\binom LN}$ を $L$ 個のサイトに作用させる。サイトを一つ加える写像
$x\mapsto x\otimes1_2$ は $\mathcal A_L$ を $\mathcal A_{L+1}$ の中に写す。$x\otimes1$ は
$S^z\otimes1+1\otimes s^z$ と可換だからである。これにより塔 $\mathcal A_0\subset\mathcal A_1\subset\mathcal A_2\subset\cdots$ が得られる。

<a id="r-9-2"></a>

**命題 9.2.** 包含 $\mathcal A_L\subset\mathcal A_{L+1}$ の包含行列は $(L+2)\times(L+1)$
行列 $\Lambda_L$ であり、$\lambda_{N,N}=\lambda_{N,N-1}=1$、他の成分はすべて $0$ である。ブロックの
大きさはパスカルの法則 $\binom{L+1}N=\binom LN+\binom L{N-1}$ に従い、塔のブラッテリ図は
パスカルの三角形であり、$\Lambda_L$ の二部グラフは $2L+3$ 個の頂点をもつ道である。

*証明.* $L+1$ 個のサイト上で上向きスピンを $N$ 個もつ配置は、最初の $L$ 個のサイトに上向きスピンが
$N$ 個あって新しいスピンが下向きのものと、$N-1$ 個あって新しいスピンが上向きのものである。この二つの部分空間の上で
$x\otimes1$ はそれぞれ $x_N$ と $x_{N-1}$ として作用するので、$\mathcal A_{L+1}$ のブロック $N$ は
$\mathcal A_L$ のブロック $N$ と $N-1$ を一度ずつ含む。次元を足し合わせるとパスカルの法則になる。二部グラフでは
$\mathcal A_{L+1}$ のセクター $N$ は $\mathcal A_L$ の $N-1$ と $N$ に結ばれ、これにより道
$0'-0-1'-1-\cdots-L'-L-(L+1)'$ が得られる（プライムは $\mathcal A_{L+1}$ のもの）。$\square$

ブラッテリ図の頂上からレベル $L$ のセクター $N$ に至る道は一つの配置 $b$ であり、その数は
$\binom LN$、すなわちブロックの大きさである。包含行列 $\Lambda$ をもつ連結な包含に対して、**指数** を
$\|\Lambda\|^2$（作用素ノルム）と **定義する**。$M_k\otimes1\subset M_k\otimes M_m$ に対しては包含行列は $1\times1$ 行列 $(m)$
であるから、指数は 7.2節の $m^2$ である。`Inclusion.index` がこれを計算する。（この名称は V. Jones の指数理論、
*Invent. Math.* 72 (1983) 1 から借りたものである。本講義では、有限次元の環の包含の指数とはこの定義のことであり、
それ以上の意味はない。）

*ライブラリでは：* 指数は `Inclusion.index` である — [チュートリアル、9.1節](../../tutorial/ja/09_spin_chains.md#s-9-1)。

<a id="r-9-3"></a>

**補題 9.3.** 頂点 $1,\dots,n$ をもつ道の隣接行列 $A_n$ は固有値
$2\cos\frac{\pi k}{n+1}$、$k=1,\dots,n$ と、固有ベクトル $v_k(j)=\sin\frac{\pi jk}{n+1}$ をもつ。したがって
$\|A_n\|=2\cos\frac\pi{n+1}$ である。

*証明.* $v_k(0)=v_k(n+1)=0$ とおく。$\sin(a-b)+\sin(a+b)=2\sin a\cos b$ により
$(A_nv_k)(j)=v_k(j-1)+v_k(j+1)=2\cos\frac{\pi k}{n+1}\,v_k(j)$ であり、$v_k(1)=\sin\frac{\pi k}{n+1}\neq0$ であるから $v_k\ne0$ である。
$n$ 個の固有値は相異なる（$\cos$ は $(0,\pi)$ 上で狭義単調減少）ので、これらで固有値はすべてであり、
対称行列 $A_n$ のノルムは固有値の絶対値の最大値で、$k=1$ と $k=n$ で達成される。$\square$

<a id="r-9-4"></a>

**系 9.4.** $[\mathcal A_{L+1}:\mathcal A_L]=\|\Lambda_L\|^2=4\cos^2\frac\pi{2L+4}$ であり、これは増加して $4$ に近づく。

*証明.* 頂点を命題 9.2 の道に沿って並べると、隣接行列
$B=\begin{pmatrix}0&\Lambda_L\\\Lambda_L^T&0\end{pmatrix}$ は $A_{2L+3}$ になる（置換による相似）。補題 1.1 により
$\|B\|^2=\|B^2\|$ かつ $\|\Lambda_L\Lambda_L^T\|=\|\Lambda_L^T\Lambda_L\|=\|\Lambda_L\|^2$ であり、$B^2=\Lambda_L\Lambda_L^T\oplus\Lambda_L^T\Lambda_L$ であるから
$\|B\|=\|\Lambda_L\|$ である。$n=2L+3$ とした補題 9.3 が値を与える。$\square$

*計算：* `inclusions.chain_tower` と `inclusions.bratteli_diagram` による、$L\le8$ に対する命題 9.2、補題 9.3、系 9.4 — [チュートリアル、9.1節](../../tutorial/ja/09_spin_chains.md#s-9-1)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. 相互作用が $z$ 軸まわりの回転で不変であれば、
> $S^z$ は保存される。$N$ は上向きスピンの数、あるいは 9.4節のジョルダン–ウィグナー写像の後では粒子の数を数える。
> 塔の中では、長い鎖の各電荷セクターに短い鎖の二つのセクターが流れ込む。

<a id="s-9-2"></a>

## 9.2 並進、運動量、マーシャルの符号則

周期的な鎖の上で、$T$ を並進 $T|b_0b_1\dots b_{L-1}\rangle=|b_{L-1}b_0\dots b_{L-2}\rangle$ とする。
これはサイト $i$ をサイト $i+1$ に移す。$T$ は置換行列であるからユニタリであり、$T^L=1$ かつ
$Ts^a_iT^*=s^a_{i+1}$ である。したがって $TS^zT^*=S^z$ であり、一様な磁場 $h_i=h$ に対しては $THT^*=H$ である：$T$ はボンドをボンドに写す。
したがってその固有値は $e^{ik}$ の形であり、**運動量** は $k=2\pi m/L$、$m=0,\dots,L-1$ である。

<a id="r-9-5"></a>

**命題 9.5（運動量セクター）.** $N$ を固定する。上向きスピンを $N$ 個もつ配置 $a$ に対し、$R_a$ を
その軌道 $\{T^ja\}$ の大きさとする。(a) $R_a$ は $L$ を割り切る。$mR_a\equiv0\pmod L$ を満たす $k=2\pi m/L$ に対し、ベクトル
$|a,k\rangle=\sum_{j=0}^{R_a-1}e^{-ikj}T^j|a\rangle$ は零でなく $T|a,k\rangle=e^{ik}|a,k\rangle$ を満たし、これら
$R_a$ 個のベクトルは軌道の張る空間の直交基底をなす。(b) セクター $N$ 上での $T$ の固有値 $e^{2\pi im/L}$ の
重複度は $d_{N,m}=\#\{\text{軌道 }:mR_a\equiv0\pmod L\}$ であり、
$S^z$ および $T$ と可換な任意の作用素——たとえば並進不変な XXZ 鎖——はブロック対角であり、$\bigoplus_N\bigoplus_mM_{d_{N,m}}$ に属する。

*証明.* (a) 巡回群 $\{T^j\}\cong\mathbb Z_L$ における $a$ の安定化群は部分群であり、$R_a$ はその
指数であるから $L$ の約数である。ベクトル $T^j|a\rangle$（$0\le j<R_a$）は正規直交であり、$T^{R_a}|a\rangle=|a\rangle$ である。
$e^{ikR_a}=1$ ならば、$j=R_a$ の項は $j=0$ の項に等しいので
$T|a,k\rangle=e^{ik}\sum_{j=1}^{R_a}e^{-ikj}T^j|a\rangle=e^{ik}|a,k\rangle$ である。$|a\rangle$ の係数は $1$ である。許される
$m\in\{0,\dots,L-1\}$ は $L/R_a$ の倍数で、$R_a$ 個あり、ユニタリな $T$ の相異なる固有値に対応するので、これらのベクトルは直交し、
軌道の張る $R_a$ 次元の空間を張る。(b) セクター $N$ はその軌道の張る空間の直交和であり、それぞれは $T$ で
不変である。(a) により、固有値 $e^{2\pi im/L}$ は、$mR_a\equiv0\pmod L$ ならば $a$ の軌道の張る空間にちょうど一回現れ、
そうでなければ現れない。$S^z$ および $T$ と可換な $X$ はそれらのスペクトル射影 $P_N$、$Q_m$ と可換である
（定理 1.5(d)）から、$X=\sum_{N,m}P_NQ_mXP_NQ_m$ である。$\square$

*計算：* $L=12$ に対する命題 9.5(b)、および `SpinChain.xxz_momentum` が構成する運動量ブロック — [チュートリアル、9.2節](../../tutorial/ja/09_spin_chains.md#s-9-2)。

<a id="r-9-6"></a>

**補題 9.6（ペロン–フロベニウス、対称な場合）.** $M$ を、$i\ne j$ に対して $M_{ij}\ge0$ であり、
そのグラフ（$i\sim j$ であるのは $M_{ij}>0$ のとき、かつそのときに限る）が連結な実対称行列とする。このとき最大固有値 $\lambda$ は単純であり、
すべての成分が $>0$ である固有ベクトルをもつ。

*証明.* $\lambda=\max_{\|x\|=1}\langle x,Mx\rangle$ であり、最大値は
$\lambda$ 固有空間の単位ベクトルでちょうど達成される（正規直交固有基底で展開せよ）。そのような $x$ に対し、絶対値のベクトル $|x|$ は
$\langle|x|,M|x|\rangle\ge\mathrm{Re}\langle x,Mx\rangle=\lambda$ を満たす。対角項は一致し、$i\ne j$ に対して
$M_{ij}|x_i||x_j|\ge\mathrm{Re}\,M_{ij}\bar x_ix_j$ だからである。よって $|x|$ は固有ベクトルである。$|x|_i=0$ ならば
$0=(\lambda-M_{ii})|x|_i=\sum_{j\ne i}M_{ij}|x|_j$ から $i$ のすべての隣接点 $j$ について $|x|_j=0$ が従い、
連結性により $x=0$ となる。したがって零でない固有ベクトルは零成分をもたない。しかし次元 $\ge2$ の固有空間は
$x_1=0$ を満たす零でないベクトルを含む。よって $\lambda$ は単純であり、$|x|>0$ は固有ベクトルである。$\square$

符号則を、任意の二部的な結合の配置に対して証明する。鎖はその特別な場合である。

<a id="r-9-7"></a>

**定理 9.7（マーシャルの符号則）.** サイト $\{0,\dots,L-1\}$ を二つの集合 $A\sqcup B$ に分け、
$$
H=\sum_{i<j}\Big[\tfrac{J_{ij}}2\big(s^+_is^-_j+s^-_is^+_j\big)+K_{ij}\,s^z_is^z_j\Big]+\sum_ih_is^z_i
$$
とする。ここで $K_{ij},h_i$ は実数、$J_{ij}$ は実数で、$i,j$ が異なる集合にあれば $J_{ij}\ge0$、同じ集合にあれば
$J_{ij}\le0$ であり、$J_{ij}\ne0$ のとき、かつそのときに限り辺 $i\sim j$ をもつサイト上のグラフ $G$ は連結であるとする。
$N_B(b)$ を $b$ の $B$ 上の上向きスピンの数とする。各セクター $N$ において $H_N$ の最小固有値は単純であり、
その固有ベクトルは、上向きスピンを $N$ 個もつすべての配置 $b$ に対して $(-1)^{N_B(b)}\psi(b)>0$ となるように選べる。

とくに、$J>0$ で $\Delta,h_i$ が任意の XXZ 鎖について、鎖が開いているか $L$ が偶数の周期的な鎖であれば、$A$ を
偶数番目のサイト、$B$ を奇数番目のサイトとしてこれが成り立つ：このときボンド上で $J_{ij}=J$、$K_{ij}=J\Delta$、
それ以外では $0$ であり、各ボンドは偶数サイトと奇数サイトを結び（環については $L$ が偶数だから）、$G$ は道
$0-1-\cdots-(L-1)$ を含む。

*証明.* $V|b\rangle=(-1)^{N_B(b)}|b\rangle$ とおく。$V$ は対角であり、$V=V^*=V^{-1}$ である。作用素 $s^+_is^-_j$ は
上向きスピンを $j$ から $i$ へ動かす。これは $i,j$ が異なる集合にあれば $N_B$ を $\pm1$ だけ変え、そうでなければ
変えないので、それに応じて $Vs^+_is^-_jV=\mp s^+_is^-_j$ であり、$s^-_is^+_j$ についても同様である。相異なる二つの
配置 $b,b'$ がこのような項で結ばれるのは、$b'$ が $b$ から二つのサイト $i,j$ の逆向きのスピンを交換して得られる
とき、かつそのときに限る。このとき $i,j$ は $b,b'$ によって決まる（両者が異なるサイト）。よって $M=-VH_NV$ の
非対角成分は $J_{ij}/2$（異なる集合）または $-J_{ij}/2$（同じ集合）、すなわちどちらの場合も $|J_{ij}|/2\ge0$ であり、
これが正であるのは $G$ において $i\sim j$ のとき、かつそのときに限る。それ以外の仕方で異なる配置の間では $0$ である。
$M$ のグラフは連結である。実際、サイト $i,j$ のスピンの交換を $(i\,j)$ と書く。相異なる $a,b,c$ に対して
$(a\,c)=(a\,b)(b\,c)(a\,b)$ であるから、$G$ における $i$ から $j$ への道に沿った帰納法により、任意の $(i\,j)$ は $G$ の辺に
沿った交換の合成であり、そのような各交換は配置を固定する（スピンが等しい）か、$M$ のグラフの辺に沿って動かす。
上向きスピンを $N$ 個もつ $b\ne b'$ に対しては、$b$ で上向き・$b'$ で下向きのサイト $i$ と、$b$ で下向き・$b'$ で
上向きのサイト $j$ が存在する。$(i\,j)b$ が $b'$ と異なるサイトは $b$ より二つ少ないので、有限回のこのような操作で
$b$ から $b'$ に至る。補題 9.6 により、$M$ の最大固有値——$H_N$ の最小固有値の符号を変えたもの——は単純で、
固有ベクトル $\varphi>0$ をもつ。$\psi=V\varphi$ とすればよい。$\square$

<a id="r-9-8"></a>

**系 9.8（基底状態の運動量）.** $L$ が偶数、$J>0$ で磁場が一様な周期的な鎖の上で、
セクター $N$ の基底状態の運動量は、$N$ が偶数なら $k=0$、$N$ が奇数なら $k=\pi$ である。

*証明.* $T$ は $H_N$ と可換であり、基底状態 $\psi=V\varphi$ は単純であるから、$T\psi=e^{ik}\psi$ である。ここで
$(T\psi)(b)=\psi(T^{-1}b)$ であり、$T^{-1}b$ はサイト $i$ にスピン $b_{i+1}$ をもつ。$L$ は偶数なので、奇数の $i$ に対して $i+1$ は偶数であり、
$N_B(T^{-1}b)=N-N_B(b)$、したがって $(VT\psi)(b)=(-1)^{N_B(b)}(-1)^{N-N_B(b)}\varphi(T^{-1}b)=(-1)^N\varphi(T^{-1}b)$ である。また
$VT\psi=e^{ik}\varphi$ である。$\varphi$ と $\varphi\circ T^{-1}$ はともに正であるから、$e^{ik}=(-1)^N$ である。$\square$

*計算：* ランダム磁場中の開いた鎖と周期的な鎖に対する定理 9.7、および $L=12$ に対する系 9.8 — [チュートリアル、9.2節](../../tutorial/ja/09_spin_chains.md#s-9-2)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. 反強磁性体では、格子が基底状態の符号を決める：
> スピンをボンドの向こうへ動かすと符号が変わる。物理学者はこれを用いて、二部格子の反強磁性体で負符号問題の
> ない量子モンテカルロ法を実行する。反強磁性的な結合が同じ集合の二つのサイトを結ぶフラストレートした鎖では、
> 定理 9.7 の仮定が成り立たない。

<a id="s-9-3"></a>

## 9.3 $\Delta=1$ における $SU(2)$

$\Delta=1$ かつ $h=0$ のとき、$H=J\sum\mathbf s_i\cdot\mathbf s_j$ は $S^x,S^y,S^z$ と可換である（命題 4.15）から、
$H\cong\bigoplus_S1_{2S+1}\otimes H_S$ であり、ブロックは $H_S\in M_{m_S}$、$m_S=\binom L{L/2-S}-\binom L{L/2-S-1}$
（命題 4.13）である。モジュール `su2` は最高ウェイトベクトルの上で $H_S$ を構成する（命題 4.15(b)）。
実数の $J_{ij}$ をもつ任意の $H=\sum_{i<j}J_{ij}\,\mathbf s_i\cdot\mathbf s_j$ についても同じことが成り立つ：命題 4.15(a) により、
各 $\mathbf s_i\cdot\mathbf s_j$ は $\{S_a\}'$ に属する。現れるスピン $S$ に対し、$E_0(S)$ を $H_S$ の最小固有値とする。

*ライブラリでは：* ブロック $H_S$ は `SU2Chain` が構成する — [チュートリアル、9.3節](../../tutorial/ja/09_spin_chains.md#s-9-3)。

<a id="r-9-9"></a>

**定理 9.9（リーブ–マティス）.** $H=\sum_{i<j}J_{ij}\,\mathbf s_i\cdot\mathbf s_j$ とし、$A,B,J_{ij}$ は定理 9.7 と同様
（したがって $K_{ij}=J_{ij}$、$h_i=0$）、$n_A=|A|\ge n_B=|B|$、$S_*=\frac{n_A-n_B}2$ とする。
(a) $S_*\le S\le\frac L2$（$S-\frac L2$ は整数）に対し、セクター $S^z=S$ の基底状態は全スピン $S$ をもつ。
$0\le M\le S_*$（$M-\frac L2$ は整数）に対し、セクター $S^z=M$ の基底状態は全スピン $S_*$ をもつ。
(b) したがって $S\ge S_*$ に対して $E_0(S)$ は $H_S$ の単純な固有値であり、$E_0(S_*)<E_0(S_*+1)<\dots<E_0(\frac L2)$、
かつ $S<S_*$ に対して $E_0(S)>E_0(S_*)$ である：$H$ の基底空間はスピン $S_*$ の単一の多重項である。
$J>0$ のハイゼンベルク鎖（開いているか、$L$ が偶数の周期的な鎖）では、$A$ と $B$ は偶数番目と奇数番目のサイトであり、
$S_*$ は $L$ が偶数なら $0$、$L$ が奇数なら $\frac12$（開いた鎖）であって、$E_0(S_*)<E_0(S_*+1)<\dots<E_0(\frac L2)$ が成り立つ。

*証明.* (a) $M\ge0$、$M-\frac L2\in\mathbb Z$、$s=\max(M,S_*)$ とし、$\psi$ をセクター $S^z=M$ の基底状態とする。
これは定理 9.7 により単純である：$\varphi>0$ として $\psi=V\varphi$。$\mathbf S^2$ は $H$ および $S^z$ と可換であるから、
セクターの一次元の基底空間を保つ：$\psi$ はある全スピン $S'$ をもつ。セクターのベクトル $\chi$ で、
$\mathbf S^2\chi=s(s+1)\chi$ かつ $\langle\psi,\chi\rangle\ne0$ となるものを構成する。$\mathbf S^2$ の相異なる固有値に属する
固有ベクトルは直交するので、これにより $S'=s$ が従う。$p=\frac L2-s$ とおく。$s\ge S_*$ は
$p\le\frac{n_A+n_B}2-\frac{n_A-n_B}2$ を意味するので、$0\le p\le n_B$ である。相異なる $a_1,\dots,a_p\in A$ と
$b_1,\dots,b_p\in B$ を選び、残りの $2s$ 個のサイトを $R$ とし、
$$
\chi_0=\bigotimes_{l=1}^p\tfrac1{\sqrt2}\big(|{\uparrow}_{a_l}{\downarrow}_{b_l}\rangle-|{\downarrow}_{a_l}{\uparrow}_{b_l}\rangle\big)\otimes|{\uparrow\cdots\uparrow}\rangle_R ,
\qquad \chi=(S^-)^{s-M}\chi_0 .
$$
とおく。各一重項は $s^\pm_{a_l}+s^\pm_{b_l}$ で消え、$R$ の上向きスピンは $s^+_i$ で消えるので、$S^+\chi_0=0$ である。
$S^z\chi_0=s\chi_0$ と合わせて $\mathbf S^2\chi_0=(S^-S^++(S^z)^2+S^z)\chi_0=s(s+1)\chi_0$（4.5節）を得、$\mathbf S^2$ は $S^-$ と
可換であるから $\mathbf S^2\chi=s(s+1)\chi$ である。$S^-=\sum_l(s^-_{a_l}+s^-_{b_l})+S^-_R$ と書くと、各項は可換で、
最初の項たちは一重項を消すので、$\chi$ は一重項と $(S^-_R)^{s-M}|{\uparrow\cdots\uparrow}\rangle_R=(s-M)!\sum_D|D\rangle$ の積である。
ここで和は、ちょうど $s-M$ 個のスピンが下向きである $R$ の配置にわたる（$s-M\le2s=|R|$ であるから、そのような配置は存在する）。
したがって $\chi\ne0$、$S^z\chi=M\chi$ であり、$\chi$ の係数は実数である。その符号：一重項は上向きスピンが $A$ に
あるとき $+$、$B$ にあるとき $-$ を寄与し、$R$ の部分は正である。$M=s$ ならば $R$ のサイトはすべて上向きであり、
$M<s$ ならば $s=S_*$、$p=n_B$ で $R\subseteq A$ である。どちらの場合も $N_B(b)$ は、上向きスピンが $B$ にある一重項の
数に定数 $|R\cap B|$（$M=s$ のとき）または $0$ を加えたものであるから、$(-1)^{N_B(b)}\chi(b)$ は $\chi(b)\ne0$ となる
すべての $b$ に対して同じ符号をもつ。よって $\langle\psi,\chi\rangle=\langle\varphi,V\chi\rangle=\pm\sum_b\varphi(b)|\chi(b)|\ne0$ である。

(b) 命題 4.15(a)、系 4.12、定理 4.11 により、スピン $S'$ の多重項はセクター $S^z=M$ と、$S'\ge|M|$ ならば
ちょうど一つのベクトルで交わり、そうでなければ交わらないので、セクター $S^z=M\ge0$ 上の $H$ のスペクトルは
$H_{S'}$（$S'\ge M$）のスペクトルの（重複度込みの）和集合である。(a) により、その最小固有値は単純で、スピン
$s=\max(M,S_*)$ に属する：よって $E_0(s)$ は $H_s$ において単純であり、$S'\ge M$、$S'\ne s$ なるすべての $S'$ に対して
$E_0(s)<E_0(S')$ である。$M=S\ge S_*$ に対してはこれは単純性と、$S'>S$ に対する $E_0(S)<E_0(S')$ を与え、
$M=S<S_*$ に対しては $E_0(S_*)<E_0(S)$ を与える。したがって $E_0(S_*)$ は $H$ の最小固有値であり、$H_{S_*}$ において
のみ達成され、そこで単純である。鎖については、偶数番目のサイトは奇数番目のサイトより $0$ 個（$L$ が偶数）または
$1$ 個（$L$ が奇数）多い。$\square$

定理 9.7 は W. Marshall（*Proc. Roy. Soc. A* 232 (1955) 48）に、定理 9.9 は E. Lieb と D. Mattis
（*J. Math. Phys.* 3 (1962) 749）にさかのぼる。偶数のハイゼンベルク環の **一重項–三重項ギャップ** $E_0(1)-E_0(0)$ は
定理 9.9 により正である。

*計算：* $L=12$ の $U(1)$ セクターによる定理 9.9 の検証、一重項–三重項ギャップ、および鎖ではない二部グラフ上での定理 9.7 と 9.9 — [チュートリアル、9.3節](../../tutorial/ja/09_spin_chains.md#s-9-3)。

**相関.** すべての偶数のハイゼンベルク環において、基底状態の相関 $\langle\mathbf s_0\cdot\mathbf s_r\rangle$ は
符号が交互に変わる。定理 9.7 と 9.9 がこれを強制する。$\psi$ を実で正規化された
基底状態とする。これは一重項であり（定理 9.9）、したがってセクター $S^z=0$ の基底状態である。$\mathbf S^2\psi=0$ から
$\|S^\pm\psi\|^2=\langle\psi,S^\mp S^\pm\psi\rangle=\langle\psi,(\mathbf S^2-(S^z)^2\mp S^z)\psi\rangle=0$ であるから $S^\pm\psi=0$ であり、
任意の $X$ に対して $\langle\psi,[S^+,[S^-,X]]\psi\rangle=0$ である（四つの項はどれも $\psi$ の隣に $S^\pm$ をもつ）。
$X=s^z_0s^z_r$、$r\ne0$ に対しては、関係式 $[s^\pm,s^z]=\mp s^\pm$、$[s^+,s^-]=2s^z$ から $[S^-,X]=s^-_0s^z_r+s^z_0s^-_r$ および
$[S^+,[S^-,X]]=4s^z_0s^z_r-s^+_0s^-_r-s^-_0s^+_r$ を得る。よって $\langle s^z_0s^z_r\rangle=\frac14\langle s^+_0s^-_r+s^-_0s^+_r\rangle$ であり、
$$
\langle\mathbf s_0\cdot\mathbf s_r\rangle=\langle s^z_0s^z_r\rangle+\tfrac12\langle s^+_0s^-_r+s^-_0s^+_r\rangle
=\tfrac34\sum_{b_0\ne b_r}\psi(b)\psi(b') ,
$$
となる。ここで $b'$ は $b$ からサイト $0$ と $r$ のスピンを交換して得られる。定理 9.7 により $\psi(b)\psi(b')$ は
$(-1)^{N_B(b)+N_B(b')}$ の符号をもち、$N_B(b')-N_B(b)$ は奇数の $r$ では $\pm1$、偶数の $r$ では $0$ である。和は空でないので、
$0<r<L$ に対して $\langle\mathbf s_0\cdot\mathbf s_r\rangle$ は零でなく、符号 $(-1)^r$ をもつ。

*計算：* `SU2Chain.spin_correlation` による $L=12$ の相関 $\langle\mathbf s_0\cdot\mathbf s_r\rangle$。これはウィグナー–エッカートの定理（定理 B.16）を用いる — [チュートリアル、9.3節](../../tutorial/ja/09_spin_chains.md#s-9-3)。

<a id="s-9-4"></a>

## 9.4 厳密な検証：自由フェルミオンとベーテ仮説

**ジョルダン–ウィグナー.** 開いた鎖の上で $c_j=\big(\prod_{l<j}(-1)^{n_l}\big)s^-_j$ とおく。

<a id="r-9-10"></a>

**定理 9.10（XX 鎖は自由である）.** (a) $\{c_i,c_j^*\}=\delta_{ij}$、$\{c_i,c_j\}=0$、$c_j^*c_j=n_j$。
(b) 開いた XX 鎖は $H=\frac J2\sum_{j<L-1}(s^+_js^-_{j+1}+s^-_js^+_{j+1})=\sum_{i,j}h_{ij}c_i^*c_j$ であり、
$h=\frac J2A_L$ である。(c) $u_k$ を $h$ の正規直交な実固有ベクトルで固有値が
$\varepsilon_k=J\cos\frac{\pi k}{L+1}$ のもの（補題 9.3）とし、$b_k=\sum_ju_k(j)c_j$ とおく。$K=\{k_1<\dots<k_N\}$ に対し、ベクトル
$b_K^*|0\rangle=b_{k_1}^*\cdots b_{k_N}^*|0\rangle$（$|0\rangle$ はすべてのスピンが下向きの状態）はセクター
$N$ の正規直交固有基底をなし、エネルギーは $\sum_{k\in K}\varepsilon_k$ である。

*証明.* (a) 異なるサイト上の作用素は可換であり、一つのサイト上では $(-1)^n$ は $s^\pm$ と反可換である。
$i<j$ のとき、$c_i$ と $c_j^*$（または $c_j$）はサイト $i$ で $s^-_i$ と $(-1)^{n_i}$ として出会い——一回の反交換——
それ以外では可換であるから、反可換である。$i=j$ のとき、ストリングの二乗は $1$ であり、
$c_jc_j^*+c_j^*c_j=s^-_js^+_j+s^+_js^-_j=1$、$c_j^*c_j=s^+_js^-_j=n_j$、$c_j^2=(s^-_j)^2=0$ である。$i>j$ の場合は
$i,j$ の役割を入れ替えれば（そして共役をとれば）従う。(b) $s^+_j\neq0$ となるのは $n_j=0$ の上だけであるから
$c_j^*c_{j+1}=s^+_j(-1)^{n_j}s^-_{j+1}=s^+_js^-_{j+1}$ である。共役をとればよい。(c) $\{b_k,b_l^*\}=\sum_ju_k(j)u_l(j)=\delta_{kl}$ かつ
$\{b_k,b_l\}=0$ であり、$h$ のスペクトル分解により $H=\sum_k\varepsilon_kb_k^*b_k$ である。正準反交換関係（CAR）から
$[b_l^*b_l,b_k^*]=\delta_{lk}b_k^*$、よって $[H,b_k^*]=\varepsilon_kb_k^*$ である。$H|0\rangle=0$ と合わせて
$Hb_K^*|0\rangle=(\sum_K\varepsilon_k)b_K^*|0\rangle$ を得る。$b_k|0\rangle=0$（各 $c_j$ はスピンを下げる）であるから、CAR を用いて消滅演算子を
右へ移すと、$N$ についての帰納法により $\langle b_{K'}^*0,b_K^*0\rangle=\delta_{KK'}$ を得る。各 $c_j^*$ は
$N$ を一つ上げるので、$b_K^*|0\rangle$ はセクター $|K|$ に属する：その次元の空間の中の $\binom LN$ 個の正規直交ベクトルである。$\square$

*計算：* $L=10$ の開いた XX 鎖の $1024$ 個の準位すべてを定理 9.10 と照合する — [チュートリアル、9.4節](../../tutorial/ja/09_spin_chains.md#s-9-4)。

**マグノン.** $h=0$ の環の上で、完全に偏極した状態 $|0\rangle$ のエネルギーは $E_F=J\Delta L/4$ である。$x$ にある一つの上向き
スピンは二つのボンドを壊す：$H|x\rangle=(E_F-J\Delta)|x\rangle+\frac J2(|x-1\rangle+|x+1\rangle)$（$L\ge3$）。

<a id="r-9-11"></a>

**命題 9.11（一マグノン）.** 環のセクター $N=1$ において、$|k\rangle=\sum_xe^{-ikx}|x\rangle$、
$k=2\pi m/L$ は $T|k\rangle=e^{ik}|k\rangle$ と $H|k\rangle=\big(E_F+J(\cos k-\Delta)\big)|k\rangle$ を満たす。

*証明.* $e^{ikL}=1$ を用いると $\sum_xe^{-ikx}(|x-1\rangle+|x+1\rangle)=\sum_y(e^{-ik(y+1)}+e^{-ik(y-1)})|y\rangle=2\cos k\,|k\rangle$ である。
同じ添字のつけかえで $T|k\rangle=e^{ik}|k\rangle$ が得られる。$\square$

二つのマグノンは散乱する。$x_1<x_2$ に対し、上向きスピンが $x_1,x_2$ にある状態を $|x_1,x_2\rangle$ と書く。隣り合う対は
二つのボンドを、離れた対は四つのボンドを壊す。

<a id="r-9-12"></a>

**命題 9.12（二マグノン；ベーテ方程式）.** $k_1,k_2\in\mathbb C$、$K=k_1+k_2$ とし、
$(A_{12},A_{21})\neq0$ が (i) $A_{12}(1+e^{iK}-2\Delta e^{ik_2})+A_{21}(1+e^{iK}-2\Delta e^{ik_1})=0$ と
(ii) $A_{21}=e^{ik_2L}A_{12}$、$A_{12}=e^{ik_1L}A_{21}$ を満たすとする。$\mathbb Z^2$ 上で
$\psi(x_1,x_2)=A_{12}e^{i(k_1x_1+k_2x_2)}+A_{21}e^{i(k_2x_1+k_1x_2)}$ とし、
$\Psi=\sum_{0\le x_1<x_2<L}\psi(x_1,x_2)|x_1,x_2\rangle$ とおく。$\Psi\neq0$ ならば
$H\Psi=\big(E_F+J(\cos k_1+\cos k_2-2\Delta)\big)\Psi$ である。分母が消えなければ、(i) と (ii) は次のことをいう：
$$
e^{ik_1L}=-\frac{1+e^{iK}-2\Delta e^{ik_1}}{1+e^{iK}-2\Delta e^{ik_2}},\qquad
e^{ik_2L}=-\frac{1+e^{iK}-2\Delta e^{ik_2}}{1+e^{iK}-2\Delta e^{ik_1}} .
$$

*証明.* $E=E_F+J(\cos k_1+\cos k_2-2\Delta)$ とおく。(1) すべての $(x_1,x_2)\in\mathbb Z^2$ に対し、各平面波について
$\psi(x_1\pm1,x_2)$ と $\psi(x_1,x_2\pm1)$ の和は $2(\cos k_1+\cos k_2)\psi(x_1,x_2)$ である。(2) 代入すると、
$\psi(x,x)+\psi(x+1,x+1)-2\Delta\psi(x,x+1)$ は (i) の左辺の $e^{iKx}$ 倍であるから、すべての $x$ に対して消える。
(3) (ii) により $\psi(x_2,x_1+L)=A_{12}e^{ik_2L}e^{i(k_2x_1+k_1x_2)}+A_{21}e^{ik_1L}e^{i(k_1x_1+k_2x_2)}=\psi(x_1,x_2)$ であり、
これを二回適用して $\psi(x_1+L,x_2+L)=\psi(x_1,x_2)$ を得る。配置 $\{y_1<y_2\}$ の
$D=\{x_1<x_2<x_1+L\}$ における代表元は $(y_1+nL,y_2+nL)$ と $(y_2+nL,y_1+(n+1)L)$ であるから、$\psi$ は配置の関数である。
$(x_1,x_2)\in D$ とする。$2\le x_2-x_1\le L-2$ ならば、四つのホッピングはすべて $D$ にとどまり、対角成分は $E_F-2J\Delta$ であるから、(1) により
$(H\Psi)(x_1,x_2)=(E_F-2J\Delta)\psi+\frac J2\cdot2(\cos k_1+\cos k_2)\psi=E\psi$ である。$x_2=x_1+1$ ならば、ホッピング
$(x_1+1,x_2)$、$(x_1,x_2-1)$ は阻まれ、対角成分は $E_F-J\Delta$ である。(1) と (2) により
$(H\Psi)(x_1,x_2)=(E_F-J\Delta)\psi+\frac J2\big(2(\cos k_1+\cos k_2)\psi-\psi(x_1,x_1)-\psi(x_1+1,x_1+1)\big)=E\psi$ である。
$x_2=x_1+L-1$ ならば、代表元 $(x_2,x_1+L)$ を用いる。これは前の場合に当たる（$L\ge3$）。$\square$

命題 9.12 は、解が固有ベクトルを与えることをいうだけであり、すべての固有ベクトルがこのようにして得られるとはいっていない。
$k_{1,2}=K/2\pm q$、$c=\cos\frac K2$ と書く。すると $1+e^{iK}-2\Delta e^{ik_{1,2}}=2e^{iK/2}(c-\Delta e^{\pm iq})$ であるから、
$A_{12}=c-\Delta e^{iq}$、$A_{21}=-(c-\Delta e^{-iq})$ は (i) を解く。(ii) と $(A_{12},A_{21})\ne0$ により両方の係数は
零でないので、(ii) の二つの式を掛け合わせると
$e^{iKL}=1$、$K=2\pi n/L$ となり、このとき (ii) は $F(q)=(-1)^ne^{iqL}(c-\Delta e^{-iq})+c-\Delta e^{iq}=0$ と同値である。
これは一変数の一つの方程式である。その根を、$c-\Delta e^{\pm iq}$ が互いに複素共役であるために
$e^{-iqL/2}F(q)$ が実数または純虚数となる $q\in(0,\pi)$ で探し（散乱状態）、また $F$ が実数となる直線
$q=i\kappa$、$\pi+i\kappa$（$\kappa>0$）の上で探す（束縛状態）。

偶数の $L$ と $\Delta\ne0$ では、ベーテ・ベクトルは完全ではない。$\Delta\ne0$ ならば、
$\Phi=\sum_x(-1)^x|x,x+1\rangle$ は $T\Phi=-\Phi$ と $H\Phi=(E_F-J\Delta)\Phi$ を満たす（$|x-1,x+1\rangle$ と
$|x,x+2\rangle$ へのホッピングは対ごとに打ち消し合う）。ベーテ・ベクトルは $\psi(x_1-1,x_2-1)=e^{-iK}\psi(x_1,x_2)$ で
あるから $T\Psi=e^{-iK}\Psi$ を満たす。$e^{-iK}\ne-1$ ならば、$\Psi$ と $\Phi$ はユニタリな $T$ の相異なる固有値に属する
固有ベクトルであるから $\Psi\perp\Phi$ である。$e^{iK}=-1$ ならば $\cos k_1+\cos k_2=2\cos\frac K2\cos\frac{k_1-k_2}2=0$ で
あり、$\Psi$ のエネルギーは $E_F-2J\Delta\ne E_F-J\Delta$ であるから、$H$ の相異なる固有値に属する固有ベクトルとして
$\Psi\perp\Phi$ である。よって $\Phi$ は命題 9.12 のすべてのベーテ・ベクトルと直交し、したがってそれらは偶数の環の
セクター $N=2$ を決して張らない。二マグノンの計算は
H. Bethe（*Z. Phys.* 71 (1931) 205）の仮説の $N=2$ の場合である。それ以外の場合は用いない。

*計算：* $L=9$、$\Delta=0.6$ に対する $F$ の根（セクター $N=2$ の $36$ 個の状態すべてを与える）と、$L=8$ に対するベクトル $\Phi$ — [チュートリアル、9.4節](../../tutorial/ja/09_spin_chains.md#s-9-4)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. マグノンはスピン波である：一つの反転したスピンが
> 運動量 $k$ をもって環全体に非局在化したものである。二つのマグノンは位相のずれを伴って散乱する（実数の $q$）か、
> 束縛される（複素数の $q$）。

<a id="s-9-5"></a>

## 9.5 ランチョス法とクリロフ法

セクター $N=L/2$ の次元は $L=16$ で $12870$、$L=20$ で $184756$ であるが、$H$ は疎であり（各行に $O(L)$ 個の非零
成分、`xxz_sparse`）、クリロフ法は積 $Hv$ だけを必要とする。$H$ をサイズ $n$ のエルミート行列とし、
$v\ne0$ とする。**クリロフ空間** は $\mathcal K_m=\mathrm{span}\{v,Hv,\dots,H^{m-1}v\}$ である。

<a id="r-9-13"></a>

**命題 9.13（ランチョス）.** $\dim\mathcal K_m=m$ と仮定する。$q_1=v/\|v\|$ とし、$j<m$ に対して $q_{j+1}$ を
$Hq_j$ の $\mathcal K_j$ に直交する成分を正規化したものとする。このとき $q_1,\dots,q_m$ は
$\mathcal K_m$ の正規直交基底である。$Q_m=(q_1\cdots q_m)$ とすると、$T_m=Q_m^*HQ_m$ は実対称三重対角行列で、対角成分は
$\alpha_j=\langle q_j,Hq_j\rangle$、非対角成分は $\beta_j=\langle q_{j+1},Hq_j\rangle>0$ である。また $HQ_m=Q_mT_m+r_me_m^T$ であり、ここで
$r_m\perp\mathcal K_m$ が消えるのは $\mathcal K_{m+1}=\mathcal K_m$ のとき、かつそのときに限る（そうでなければ $r_m=\beta_mq_{m+1}$、$q_{m+1}=r_m/\|r_m\|$）。

*証明.* $\mathrm{span}\{q_1,\dots,q_j\}=\mathcal K_j$ と仮定する（$j=1$ では正しい）。$Hq_i\in\mathcal K_{i+1}$ であるから
$H\mathcal K_j\subseteq\mathcal K_j+\mathbb CHq_j$ であり、したがって $\mathcal K_{j+1}=\mathcal K_j+\mathbb CH^jv\subseteq\mathcal K_j+\mathbb CHq_j\subseteq\mathcal K_{j+1}$、
すなわち $\mathcal K_{j+1}=\mathcal K_j+\mathbb CHq_j$ である。$j<m$ では次元が増えるので、$Hq_j$ の
$\mathcal K_j$ に直交する成分は零でなく、帰納法が進む。$j=m$ ではその成分 $r_m$ が消えるのは
$\mathcal K_{m+1}=\mathcal K_m$ のとき、かつそのときに限る。$i>j+1$ に対しては $q_i\perp\mathcal K_{j+1}\ni Hq_j$ であるから $(T_m)_{ij}=0$ であり、$T_m$ は
エルミートであるから $(T_m)_{ji}=0$ である。$\alpha_j$ は実数であり、$\beta_j=\langle q_{j+1},Hq_j\rangle$ はその成分のノルムで $>0$ である。よって
$j<m$ に対して $Hq_j=\beta_{j-1}q_{j-1}+\alpha_jq_j+\beta_jq_{j+1}$、また $Hq_m=\beta_{m-1}q_{m-1}+\alpha_mq_m+r_m$ である。$\square$

$T_m$ の固有値 $\theta^{(m)}_1\le\dots\le\theta^{(m)}_m$ を **リッツ値** という。

<a id="r-9-14"></a>

**定理 9.14（単調収束）.** $\lambda_1\le\dots\le\lambda_n$ を $H$ の固有値とし、
$\dim\mathcal K_m=m$ とする。このとき $\theta^{(m)}_j\ge\lambda_j$ であり、$\dim\mathcal K_{m+1}=m+1$ ならば $\theta^{(m)}_j\ge\theta^{(m+1)}_j$
（$j\le m$）である。$v$ が $\lambda_1$ の固有空間に零でない成分をもち、$\mathcal K_{m+1}=\mathcal K_m$ ならば
$\theta^{(m)}_1=\lambda_1$ である。

*証明.* $Q=Q_m$ に適用した系 A.9（圧縮に対するクーラン–フィッシャー、定理 A.8）により、
$\theta^{(m)}_j=\min_{U}\max_{x\in U,\|x\|=1}\langle x,Hx\rangle$ であり、最小は $j$ 次元の $U\subseteq\mathcal K_m=\mathrm{ran}\,Q_m$ に
わたってとる。同じ式を $U\subseteq\mathcal K_{m+1}$ にわたってとると $\theta^{(m+1)}_j$ が、$U\subseteq\mathbb C^n$ にわたってとると
$\lambda_j$ が得られる（定理 A.8）。より多くの部分空間にわたって最小化すれば値は小さくなる。
$\mathcal K_{m+1}=\mathcal K_m$ ならば $\mathcal K_m$ は $H$ 不変であり、任意の多項式 $p$ に対して $p(H)v$ を含み、
とくに $P_1v\ne0$ を含む（$P_1$ は $\lambda_1$ のスペクトル射影で、定理 1.5 と同様にラグランジュ補間により
$H$ の多項式である）。よって $x=P_1v/\|P_1v\|\in\mathcal K_m$ において $\theta^{(m)}_1\le\langle x,Hx\rangle=\lambda_1$ である。$\square$

小さな残差は近くに固有値があることを保証する：エルミートな $H$、数 $E$、$\|\psi\|=1$ に対し
$\|H\psi-E\psi\|^2=\sum_i|\lambda_i-E|^2|\langle u_i,\psi\rangle|^2\ge\min_i|\lambda_i-E|^2$（$u_i$ は正規直交固有基底）である。
しかしこの固有値が最小のものであることは保証しない。$L=4$ ではハイゼンベルク環の基底状態エネルギーは厳密に求まる：$\mathbf S_A=\mathbf s_0+\mathbf s_2$、
$\mathbf S_B=\mathbf s_1+\mathbf s_3$ とすると、環は
$H=J\,\mathbf S_A\cdot\mathbf S_B=\frac J2(\mathbf S^2-\mathbf S_A^2-\mathbf S_B^2)\ge\frac J2(0-2-2)=-2J$ である。$\mathbf S^2\ge0$ であり、
$\mathbf S_A^2,\mathbf S_B^2$ は固有値 $0,2$ をもつ（$L=2$ とした命題 4.13）からである。等号は二つの三重項の積の中のスピン $0$ の
ベクトルの上で成り立ち、そのようなベクトルは第4章の演習問題 2（$1\otimes1\cong2\oplus1\oplus0$）により存在する。

*計算：* `SpinChain.xxz_sparse` の疎なセクター上の `krylov.lanczos` による開いた XX 鎖での定理 9.14、および `krylov.ground_state` とその残差による $L=16$ までのハイゼンベルク環の基底状態 — [チュートリアル、9.5節](../../tutorial/ja/09_spin_chains.md#s-9-5)。

**時間発展.** $e^{-i\tau H}v$ に対しては $\mathcal K_m$ に射影する：$w_\tau=\|v\|\,Q_me^{-i\tau T_m}e_1$。

<a id="r-9-15"></a>

**命題 9.15.** $\dim\mathcal K_m=m$ とする。任意の $\tau$ に対して $\|w_\tau\|=\|v\|$ かつ $\langle w_\tau,Hw_\tau\rangle=\langle v,Hv\rangle$ である。
$\mathcal K_m$ が $H$ 不変ならば $w_\tau=e^{-i\tau H}v$ である。

*証明.* $Q_m$ は等長写像であり、$e^{-i\tau T_m}$ はユニタリである（$T_m$ は実対称）から、ノルムについての主張が従う。
$Q_m^*HQ_m=T_m$ は $e^{-i\tau T_m}$ と可換であるから、
$\langle w_\tau,Hw_\tau\rangle=\|v\|^2\langle e^{-i\tau T_m}e_1,T_me^{-i\tau T_m}e_1\rangle=\|v\|^2\langle e_1,T_me_1\rangle=\langle v,Hv\rangle$ である。
$\mathcal K_m$ が不変ならば $HQ_m=Q_mT_m$（命題 9.13、$r_m=0$）であるから、$H^jQ_m=Q_mT_m^j$ であり、指数級数を
足し合わせて（命題 8.2(1)）$e^{-i\tau H}Q_m=Q_me^{-i\tau T_m}$ を得る。
これを $e_1$ に作用させ、$Q_me_1=v/\|v\|$ を用いればよい。$\square$

*計算：* クリロフ空間でステップを進める `krylov.evolve` に対する命題 9.15 を、ネール状態からのクエンチで、密な伝播子と照合する — [チュートリアル、9.5節](../../tutorial/ja/09_spin_chains.md#s-9-5)。

<a id="s-9-6"></a>

## 9.6 乱れ：準位統計とエンタングルメント

**ランダム磁場ハイゼンベルク鎖** は、$[-W,W]$ 上の一様分布に従う独立な磁場 $h_i$ をもつ。これは $S^z$ を保存する
（命題 9.1）ので、準位統計はセクター $N$ ごとに計算する（定理 5.15 と 5.4節：セクターを混ぜると独立な準位を偽造して
しまう）。二つの磁場が異なれば $THT^*-H=\sum_i(h_{i-1}-h_i)s^z_i\ne0$ であるから、並進はもはや対称性ではない。

*計算：* ランダム磁場鎖の平均比 $\langle r\rangle$ を、GOE 行列および独立な準位と比較する（数値的観察）— [チュートリアル、9.6節](../../tutorial/ja/09_spin_chains.md#s-9-6)。

**自由フェルミオンのエンタングルメント.** 開いた XX 鎖の固有状態について、最初の $n$ 個のサイトの
エンタングルメント・エントロピーは $n\times n$ 行列によって決まる。

<a id="r-9-16"></a>

**定理 9.16（ペシェルの公式）.** $\Phi=b_{k_1}^*\cdots b_{k_N}^*|0\rangle$ を定理 9.10 のとおりとし、
$C_{ij}=\langle\Phi,c_i^*c_j\Phi\rangle$、$C_A=(C_{ij})_{i,j<n}$ とおく。このとき $0\le C_A\le1$ であり、サイト
$A=\{0,\dots,n-1\}$ のエンタングルメント・エントロピーは、$\nu$ が $C_A$ の固有値を重複度込みで動くとして、次で与えられる。
$$
S_A=\sum_\nu \eta(\nu),\qquad \eta(\nu)=-\nu\log\nu-(1-\nu)\log(1-\nu)\quad(0\log0=0).
$$

*証明.* 実ベクトル $x$ に対し $b(x)=\sum_jx_jc_j$ とおく。すると $b(x)^*b(y)^*=-b(y)^*b(x)^*$ かつ
$\{b(x),b(y)^*\}=\langle x,y\rangle$ である。$U\subseteq\mathbb R^L$ を $u_{k_a}$ の張る空間とし、$P_A$、$P_B=1-P_A$ を
$\{0,\dots,n-1\}$ とその補集合への座標射影とする。$U$ 上の形式 $\langle P_Ax,P_Ay\rangle$ を対角化すると、
$\langle P_A\phi_a,P_A\phi_b\rangle=\nu_a\delta_{ab}$ を満たす $U$ の正規直交基底 $\phi_1,\dots,\phi_N$ が得られる。ここで
$\nu_a=\|P_A\phi_a\|^2\in[0,1]$ である。したがって
$\langle P_B\phi_a,P_B\phi_b\rangle=(1-\nu_a)\delta_{ab}$ である。$(x_1,\dots,x_N)\mapsto b(x_1)^*\cdots b(x_N)^*|0\rangle$ は
多重線形かつ交代的であるから、$f_a=b(P_A\phi_a/\sqrt{\nu_a})$、$g_a=b(P_B\phi_a/\sqrt{1-\nu_a})$ として
$\Phi=\pm\prod_a\big(\sqrt{\nu_a}f_a^*+\sqrt{1-\nu_a}\,g_a^*\big)|0\rangle$ である（係数が零の項は省く）。
展開すると、部分集合 $I$ にわたる和として $\Phi=\sum_I\pm\prod_{a\in I}\sqrt{\nu_a}\prod_{a\notin I}\sqrt{1-\nu_a}\,F_I^*G_{I^c}^*|0\rangle$ となる。
ここで $F_I^*$、$G_{I^c}^*$ は順序づけられた積である。$j\ge n$ に対しては $c_j^*$ のストリングは $A$ 全体を含む：
$c_j^*=(-1)^{N_A}\otimes\tilde c_j^*$ であり、$\tilde c_j$ はサイト $\ge n$ のジョルダン–ウィグナー演算子である。一方
$i<n$ に対しては $c_i^*=\hat c_i^*\otimes1$ である。$(-1)^{N_A}$ は $|0_A\rangle$ の上で $1$ であるから、
$F_I^*G_{I^c}^*|0\rangle=\xi_I\otimes\chi_{I^c}$ であり、$\xi_I=\hat F_I^*|0_A\rangle$、$\chi_{I^c}=\tilde G_{I^c}^*|0_B\rangle$ は
定理 9.10 の CAR の計算により正規直交系である。これはシュミット分解（補題 7.7）であり、その
重み $p_I=\prod_{a\in I}\nu_a\prod_{a\notin I}(1-\nu_a)$（符号は $\xi_I$ に含める）は独立なベルヌーイ変数の分布であるから、
$\log p_I=\sum_{a\in I}\log\nu_a+\sum_{a\notin I}\log(1-\nu_a)$ と $\sum_{I\ni a}p_I=\nu_a$ により
$S_A=-\sum_Ip_I\log p_I=\sum_a\eta(\nu_a)$ である。最後に $C_{ij}=\sum_a\phi_a(i)\phi_a(j)$
（$c_i^*c_j$ を $b(\phi_a)$ で展開せよ）であるから、$C_A=\sum_a(P_A\phi_a)(P_A\phi_a)^T$ は
直交するベクトル $P_A\phi_a$ の上で固有値 $\nu_a$ をもち、それ以外では $0$ であり、$\eta(0)=0$ である。$\square$

定理 9.16 は I. Peschel（*J. Phys. A* 36 (2003) L205）による。

*計算：* $L=12$ の開いた XX 鎖の基底状態に対する定理 9.16 を `SpinChain.entanglement_entropy` と照合する — [チュートリアル、9.6節](../../tutorial/ja/09_spin_chains.md#s-9-6)。

<a id="s-9-7"></a>

## 9.7 比較の対象

*計算：* 他のソフトウェア（QuSpin、厳密対角化コード、テンソルネットワーク）との比較 — [チュートリアル、9.7節](../../tutorial/ja/09_spin_chains.md#s-9-7)。

<a id="x-9"></a>

## 演習問題

1. （運動量セクターを手で求める。）$L=4$、$N=2$ に対し、六つの配置の $T$ による軌道と
   次元 $d_{2,m}$（$m=0,\dots,3$）を求めよ。
2. （和則。）$\psi$ を $L$ が偶数のハイゼンベルク環の基底状態とする。
   $E_0=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle$ かつ $\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle=0$ を示せ。
3. （対称性が多すぎる。）セクター $N=L/2$ において、すべての係数が $1$ に等しいベクトルは
   $(S^-)^{L/2}|{\uparrow\cdots\uparrow}\rangle$ の定数倍であり、ハイゼンベルク環のエネルギー
   $JL/4$ の固有ベクトルであること、また任意のクリロフ空間 $\mathcal K_m(H,v)$ は $v$ を含む任意の $H$ 不変部分空間に含まれることを示せ。
   このベクトルから始めたとき、`krylov.ground_state` は何を返すか。

<details>
<summary>解答</summary>

1. $0011\to1001\to1100\to0110\to0011$ は $R=4$ の軌道であり、$0101\to1010\to0101$ は $R=2$ の軌道である。
   命題 9.5 により、前者はすべての $m$ に寄与し、後者は $2m\equiv0\pmod4$ を満たす $m$、すなわち $m=0,2$ に寄与する：
   $d_{2,m}=2,1,2,1$。

2. 基底状態は単純な一重項であり（定理 9.9）、$T$ は $H$ と可換であるから $T\psi=e^{ik}\psi$ である。
   $\mathbf s_i\cdot\mathbf s_{i+1}=T^i(\mathbf s_0\cdot\mathbf s_1)T^{*i}$ であるから、$\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle$ は $L$ 本のボンドすべてで等しく、
   $E_0=J\sum_i\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle$ である。$\mathbf S^2\psi=0$ であるから、
   再び並進不変性により $0=\langle\mathbf S^2\rangle=\sum_{i,j}\langle\mathbf s_i\cdot\mathbf s_j\rangle=L\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle$ である。

3. $S^-$ は上向きスピンの一つを下げる。完全に偏極した状態にこれを $L/2$ 回作用させると、上向きスピンを
   $L/2$ 個もつ各配置に $(L/2)!$ 通りの順序で到達し、どれも係数は $1$ である。$|{\uparrow\cdots\uparrow}\rangle$ では
   $L$ 本のボンドすべてが平行でエネルギーは $JL/4$ であり、$S^-$ は $H$ と可換であるから、このベクトルは同じ
   エネルギーの固有ベクトルである。$W$ が不変で $v\in W$ ならば、すべての $j$ に対して $H^jv\in W$ である。ここではすべての $m$ に対して $\mathcal K_m=\mathbb Cv$ であり、
   `ground_state` は $E=JL/4$、すなわちスペクトルの最上端を返す（命題 4.15(a) により各
   $\mathbf s_i\cdot\mathbf s_j=\frac12F_{ij}-\frac14\le\frac14$ である。$F_{ij}$ は固有値 $\pm1$ をもつ対称な置換行列だからである）。
   ランダムな初期ベクトルはこの罠を避ける：ガウス型の $v$ は密度をもつので、ほとんど確実に $\lambda_1$ の固有空間に
   零でない成分 $P_1v$ をもつ（系 A.7(c)。実の $v$ については $W=\mathbb R^n$ において真部分空間 $\{v:P_1v=0\}$ に適用する。
   実行列 $P_1\ne0$ は零でない列をもつので、これは真部分空間である）。すると定理 9.14 により、$H$ 不変なクリロフ空間
   $\mathcal K_m(H,v)$ は最小のリッツ値 $\lambda_1$ をもつが、ここでは不変な $\mathcal K_1=\mathbb Cv$ のリッツ値は $JL/4$ である。

</details>

*計算：* 演習問題 1–3 の数値的検証 — [チュートリアル、演習問題](../../tutorial/ja/09_spin_chains.md#x-9)。

**次へ：**[第10章](10_quantum_computing.md)ではセクターを用いて量子計算を扱う：デコヒーレンスフリー
部分空間とノイズレス部分系、作用素環による誤り訂正、そして保存されるハミング重みを尊重するミキサー——
本章の環 $\bigoplus_wM_{\binom nw}$ が新しい役割で現れる。
