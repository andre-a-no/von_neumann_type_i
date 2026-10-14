[English](../06_quantum_channels.md) · [Español](../es/06_quantum_channels.md) · [Français](../fr/06_quantum_channels.md) · [Deutsch](../de/06_quantum_channels.md) · [中文](../zh/06_quantum_channels.md) · **日本語** · [Русский](../ru/06_quantum_channels.md)

# 6. 量子チャネル

量子チャネルとは、完全正値かつトレース保存の線形写像である。ユニタリによる共役、中心の上への条件付き期待値、
等長作用素 $`V`$ により与えられる写像 $`\rho\mapsto\mathrm{Tr}_E(V\rho V^*)`$、そしてあるセクターから別のセクターへ
重みを移す写像などがその例である。それらの物理的な読み方（発展、結果が忘れられる測定、環境、粒子の損失）は
*物理的意味* の囲みの中でのみ論じる。本章では、
有限次元フォン・ノイマン環の間のそのような写像について標準的な記述——チョイ行列、クラウス作用素、
スタインスプリング拡張——を証明し、続いて直和に特有の事柄を調べる：セクター間の写像、第3章のトレースに関する双対、
そしてチャネルの古典的部分である。第2章と第3章ではすでにライブラリのチャネルを用いた。ここではその理論を述べる。

**本章で学ぶこと**：正値性だけでは不十分である理由、チョイ行列・クラウス作用素・超演算子・等長作用素の間を
行き来する方法、写像がセクターを変えるとき「トレース保存」がトレースに依存する理由、単位的チャネルが
どの作用素を固定するか、そしてランダムなチャネルが純粋状態に平均として何をするか。

本章の数値的検証は[チュートリアル 6](../../tutorial/ja/06_quantum_channels.md)にある。

本章を通じて、$`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$ および
$`\mathcal N=\bigoplus_{d=1}^{C'}M_{k'_d}(\mathbb C)`$ とし、それぞれの極小中心射影を $`P_c`$ と
$`Q_d`$ とする。$`A=\bigoplus_cA_c`$、$`M_{k'\times k}`$ は複素 $`k'\times k`$ 行列の空間、
$`E_{ij}=|i\rangle\langle j|`$ は行列単位、$`\mathrm{Tr}`$ は blunt トレース
$`\sum_c\mathrm{Tr}\,A_c`$、$`\mathrm{tr}_w=\sum_cw_c\mathrm{Tr}\,A_c`$ は忠実な
トレース（定理 3.1）である。ブロック行列は $`\mathbb C^n\otimes\mathbb C^k`$ に作用し、ブロックの添字を先に書く。

<a id="s-6-1"></a>

## 6.1 完全正値写像とチョイの定理

線形写像 $`\Phi:\mathcal M\to\mathcal N`$ は、$`A\ge0`$ ならば $`\Phi(A)\ge0`$ となるとき **正値** であるという。
成分が $`\mathcal M`$ に属する $`n\times n`$ 行列 $`[A_{ij}]`$ は再び行列環の直和をなす：
$`[A_{ij}]\mapsto([(A_{ij})_c]_{ij})_c`$ により $`M_n(\mathcal M)\cong\bigoplus_cM_n(M_{k_c})`$ であり、
したがって $`[A_{ij}]\ge0`$ とは、すべてのブロック行列 $`[(A_{ij})_c]\in M_{nk_c}`$ が正であることを意味する。
$`\Phi_n:[A_{ij}]\mapsto[\Phi(A_{ij})]`$ が正値のとき $`\Phi`$ は **$`n`$-正値** であるといい、すべての $`n`$ に対し
$`n`$-正値であるとき **完全正値**（CP）であるという。$`\mathrm{Tr}\,\Phi(A)=\mathrm{Tr}\,A`$ のとき
**トレース保存**（TP）、$`\Phi(1)=1`$ のとき **単位的** であるといい、CP かつ TP であるとき **チャネル** という。
任意の長方行列 $`T`$ に対し $`X\ge0`$ ならば $`TXT^*\ge0`$ であることを繰り返し用いる。実際
$`\langle y,TXT^*y\rangle=\langle T^*y,XT^*y\rangle`$ である。

<a id="r-6-1"></a>

**命題 6.1.** $`k\ge2`$ のとき、$`M_k`$ 上の転置 $`\theta(A)=A^T`$ は正値であるが
$`2`$-正値ではない。

*証明.* $`A\ge0`$ に対し $`\langle x,A^Tx\rangle=\sum_{ij}\bar x_iA_{ji}x_j=\langle\bar x,A\bar x\rangle\ge0`$ である。
$`\Omega=e_1\otimes e_1+e_2\otimes e_2\in\mathbb C^2\otimes\mathbb C^k`$ とすると
$`|\Omega\rangle\langle\Omega|=[E_{ij}]_{i,j=1}^2\ge0`$ であり、$`\theta_2`$ はこれを
$`F=[E_{ji}]_{i,j=1}^2`$ に写す。$`F`$ は $`e_1\otimes e_2-e_2\otimes e_1`$ をその $`-1`$ 倍に写す。$`\square`$

*計算：* 転置のチョイ行列は固有値 $`-1`$ をもつ — [チュートリアル、6.1節](../../tutorial/ja/06_quantum_channels.md#s-6-1)。

<a id="r-6-2"></a>

**定理 6.2（チョイ）.** 線形写像 $`\Phi:M_k\to M_{k'}`$ に対し、次は同値である：
(a) $`\Phi`$ は CP である。(b) $`\Phi`$ は $`k`$-正値である。(c) **チョイ行列**
```math
J(\Phi)=\sum_{i,j=1}^kE_{ij}\otimes\Phi(E_{ij})=[\Phi(E_{ij})]_{i,j=1}^k\in M_k(M_{k'})
```
は正である。(d) $`\Phi(A)=\sum_lK_lAK_l^*`$ を満たす $`K_1,\dots,K_r\in M_{k'\times k}`$（**クラウス作用素**）が存在し、
$`r=\mathrm{rank}\,J(\Phi)\le kk'`$ ととれる。

*証明.* (a)$`\Rightarrow`$(b) は自明である。(b)$`\Rightarrow`$(c)：
$`\Omega=\sum_ie_i\otimes e_i\in\mathbb C^k\otimes\mathbb C^k`$ に対し
$`|\Omega\rangle\langle\Omega|=[E_{ij}]\ge0`$ であり、$`J(\Phi)=\Phi_k([E_{ij}])`$ である。(c)$`\Rightarrow`$(d)：
定理 1.3 により、$`r=\mathrm{rank}\,J(\Phi)`$ として $`J(\Phi)=\sum_{l=1}^rv_lv_l^*`$ と書ける（$`v_l`$ は
非零固有値の固有ベクトルを固有値の平方根倍したもの）。任意の $`v\in\mathbb C^k\otimes\mathbb C^{k'}`$
は、ただ一つの $`K\in M_{k'\times k}`$ により $`v=\sum_ie_i\otimes Ke_i`$ と書ける（$`Ke_i`$ は $`v`$ の第 $`i`$ ブロック）。
また $`vv^*`$ の $`(i,j)`$ ブロックは $`Ke_ie_j^*K^*=KE_{ij}K^*`$ である。$`v_l`$ をこのように
$`K_l`$ を用いて書き、ブロックを比較すると $`\Phi(E_{ij})=\sum_lK_lE_{ij}K_l^*`$ を得るので、線形性により (d) が従う。
(d)$`\Rightarrow`$(a)：$`\Phi_n([A_{ij}])=\sum_l(1_n\otimes K_l)[A_{ij}](1_n\otimes K_l)^*\ge0`$。$`\square`$

<a id="r-6-3"></a>

**系 6.3.** $`\Phi(A)=\sum_{l=1}^sL_lAL_l^*`$ ならば $`s\ge\mathrm{rank}\,J(\Phi)`$ である。したがって
クラウス作用素の最小個数、すなわち **クラウス階数** は $`\mathrm{rank}\,J(\Phi)`$ である。

*証明.* 上の計算を逆にたどると、
$`w_l=\sum_ie_i\otimes L_le_i`$ として $`J(\Phi)=\sum_{l=1}^sw_lw_l^*`$ であり、これは階数 $`\le1`$ の行列 $`s`$ 個の和である。$`\square`$

*計算：* ランダムなチャネルのチョイ行列とクラウス階数 — [チュートリアル、6.1節](../../tutorial/ja/06_quantum_channels.md#s-6-1)。

直和の間の写像は、ブロック間の写像を成分とする行列である。$`\iota_c:M_{k_c}\to\mathcal M`$
を行列をブロック $`c`$ に置く写像、$`\pi_d:\mathcal N\to M_{k'_d}`$ をブロック $`d`$ を取り出す写像とし、
$`\Phi_{dc}=\pi_d\circ\Phi\circ\iota_c:M_{k_c}\to M_{k'_d}`$ とおく。$`A=\sum_c\iota_c(A_c)`$ であるから
$`\Phi(A)_d=\sum_c\Phi_{dc}(A_c)`$ である。

<a id="r-6-4"></a>

**命題 6.4.** $`\Phi:\mathcal M\to\mathcal N`$ が $`n`$-正値であるための必要十分条件は、すべての $`\Phi_{dc}`$ が $`n`$-正値であることである。したがって
$`\Phi`$ が CP であるための必要十分条件は、$`K^{dc}_l\in M_{k'_d\times k_c}`$（$`l=1,\dots,r`$）が存在して
```math
\Phi(A)_d=\sum_c\sum_{l=1}^rK^{dc}_l\,A_c\,K^{dc*}_l ,
```
となることであり、$`r=\max_{d,c}\mathrm{rank}\,J(\Phi_{dc})`$ で十分である。$`\Phi`$ が TP であるための必要十分条件は、
すべての $`c`$ に対し $`\sum_{d,l}K^{dc*}_lK^{dc}_l=1_{k_c}`$ となることである。

*証明.* すべての $`\Phi_{dc}`$ が $`n`$-正値で $`[A_{ij}]\ge0`$ ならば、各 $`X^c=[(A_{ij})_c]\ge0`$ であり、
$`\Phi_n([A_{ij}])`$ のブロック $`d`$ は $`\sum_c(\Phi_{dc})_n(X^c)\ge0`$ である。逆に、
$`M_n(M_{k_c})`$ の $`[X_{ij}]\ge0`$ に対し、$`M_n(\mathcal M)`$ において $`[\iota_c(X_{ij})]\ge0`$ である（そのブロック $`c`$ は $`[X_{ij}]`$ で、
他のブロックは零）。その像のブロック $`d`$ は $`[\Phi_{dc}(X_{ij})]`$ である。クラウス形は各組 $`(d,c)`$ に定理 6.2
を適用し、零作用素を補って共通の $`r`$ にそろえれば得られる。最後に、トレースの巡回性により
$`\mathrm{Tr}\,\Phi(A)=\sum_c\mathrm{Tr}\big(A_c\sum_{d,l}K^{dc*}_lK^{dc}_l\big)`$ であり、これがすべての $`A`$ に対し
$`\sum_c\mathrm{Tr}\,A_c`$ に等しいのは各和が $`1_{k_c}`$ であるときかつそのときに限る（行列単位をとればよい：
$`\mathrm{Tr}(E_{ij}X)=X_{ji}`$）。$`\square`$

線形写像 $`\Phi:M_k\to M_{k'}`$ がすべての $`X`$ に対し $`\Phi(X^*)=\Phi(X)^*`$ を満たすための必要十分条件は
$`J(\Phi)^*=J(\Phi)`$ である：$`J^*`$ の $`(i,j)`$ ブロックは $`\Phi(E_{ji})^*`$ であり、線形性により、これがすべての $`i,j`$ に対し
$`\Phi(E_{ij})`$ に等しいのは $`\Phi`$ が随伴を保つときかつそのときに限る。

*ライブラリでは：* セクター保存写像は `Channel`、命題 6.4 の一般の写像は `InterSectorChannel` であり、`Channel.choi()`、`Channel.superoperator()`、`Channel.from_superoperator` で表現の間を行き来する — [チュートリアル、6.1節](../../tutorial/ja/06_quantum_channels.md#s-6-1)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. 系は決して孤立していない：$`n`$ 準位の、手を触れられていない補助系と
> 相関しているかもしれない。物理的過程はこの組に $`\Phi_n=\mathrm{id}_n\otimes\Phi`$ として作用し、
> 状態を状態に写さなければならないので、$`\Phi`$ は完全正値でなければならない。転置は正値であるが、
> エンタングルした状態 $`|\Omega\rangle\langle\Omega|/2`$ を固有値 $`-1/2`$ をもつ作用素に写す。
> チョイ行列は最大エンタングル状態の像である：ただ一つの状態が完全正値性を決定する。

<a id="s-6-2"></a>

## 6.2 スタインスプリング拡張とクラウス作用素の自由度

$`M_k`$ から $`M_{k'}`$ への $`\Phi(A)=\sum_{l=1}^rK_lAK_l^*`$ に対し $`\Phi^*(B)=\sum_lK_l^*BK_l`$ とおく。
巡回性により $`\mathrm{Tr}(\Phi(A)B)=\mathrm{Tr}(A\Phi^*(B))`$ であるから、$`\Phi^*`$ はトレースに関する
$`\Phi`$ の双対であり（6.3節）、$`K_l`$ の選び方に依存しない（すべての $`A`$ に対し
$`\mathrm{Tr}(AX)=\mathrm{Tr}(AX')`$ ならば $`X=X'`$ である：$`A=E_{ji}`$ とおけばよい）。
$`\mathbb C^r\otimes\mathbb C^{k'}`$（**環境** $`\mathbb C^r`$ が先）上で
$`\mathrm{Tr}_E\,X=\sum_l(e_l^*\otimes1)X(e_l\otimes1)`$ とおく。

<a id="r-6-5"></a>

**定理 6.5（スタインスプリング、有限次元）.**
$`V=\sum_le_l\otimes K_l:\mathbb C^k\to\mathbb C^r\otimes\mathbb C^{k'}`$ とおく。このとき
```math
\Phi(A)=\mathrm{Tr}_E(VAV^*),\qquad\Phi^*(B)=V^*(1_r\otimes B)V,
```
であり、$`\Phi`$ が TP であるための必要十分条件は $`V`$ が等長作用素であることである。逆に、任意の $`V`$ は第一の式により CP 写像を定める。

*証明.* $`(e_l^*\otimes1)V=K_l`$ であるから $`\mathrm{Tr}_E(VAV^*)=\sum_lK_lAK_l^*`$ である。逆に
$`\mathrm{Tr}_E(VAV^*)`$ は $`K_l=(e_l^*\otimes1)V`$ をもつクラウス形である。次に
$`V^*(1\otimes B)V=\sum_{l,l'}\langle e_l,e_{l'}\rangle K_l^*BK_{l'}=\Phi^*(B)`$ であり、
$`V^*V=\sum_lK_l^*K_l`$ は $`\Phi`$ が TP であるときかつそのときに限り $`1`$ である（命題 6.4）。$`\square`$

したがって、ある $`V:\mathbb C^k\to\mathbb C^r\otimes\mathbb C^{k'}`$ により $`\Phi(A)=\mathrm{Tr}_E(VAV^*)`$ と
書けるような最小の $`r`$ はクラウス階数 $`\mathrm{rank}\,J(\Phi)`$ である：そのような $`V`$ はどれも $`r`$ 個の作用素をもつ
クラウス形を与えるから $`r\ge\mathrm{rank}\,J(\Phi)`$ であり（系 6.3）、$`r=\mathrm{rank}\,J(\Phi)`$ は
定理 6.2(d) により実現される。

*計算：* ランダムなチャネルに対する定理 6.5 — [チュートリアル、6.2節](../../tutorial/ja/06_quantum_channels.md#s-6-2)。

表現 $`\pi(B)=1_r\otimes B`$ を用いた式 $`\Phi^*(B)=V^*\pi(B)V`$ は、W. F. Stinespring（Proc. AMS 6 (1955) 211）が
この定理を述べた形であり、そこではここで証明した有限次元の場合よりはるかに一般的な設定が扱われている（歴史的な帰属；
さらに読むには：V. Paulsen, *Completely Bounded Maps and Operator Algebras*, Cambridge 2002, Ch. 4）。直和の間の
写像に対しては、$`\iota'_d`$ を $`\bigoplus_d\mathbb C^{k'_d}`$ の第 $`d`$ 直和成分の埋め込みとし、
$`V_c=\sum_{d,l}e_l\otimes\iota'_dK^{dc}_l:\mathbb C^{k_c}\to\mathbb C^r\otimes\bigoplus_d\mathbb C^{k'_d}`$ とおく。
すると $`(e_l^*\otimes1)V_c=\sum_d\iota'_dK^{dc}_l`$ であるから、
$`\sum_c\mathrm{Tr}_E(V_cA_cV_c^*)=\sum_{c,l}\big(\sum_d\iota'_dK^{dc}_l\big)A_c\big(\sum_{d'}\iota'_{d'}K^{d'c}_l\big)^*`$
の $`(d,d)`$ ブロックは $`\sum_{c,l}K^{dc}_lA_cK^{dc*}_l=\Phi(A)_d`$ である。また $`\iota'^*_d\iota'_{d'}=\delta_{dd'}1`$ であるから
$`V_c^*V_c=\sum_{d,l}K^{dc*}_lK^{dc}_l`$ となり、$`\Phi`$ が TP であるのはすべての $`V_c`$ が等長作用素であるときかつ
そのときに限る（命題 6.4）。`random_channel` と `random_inter_sector_channel` はこのようにしてチャネルを生成する
（6.5節）。

*ライブラリでは：* `channels.random_channel` と `channels.random_inter_sector_channel` — [チュートリアル、6.5節](../../tutorial/ja/06_quantum_channels.md#s-6-5)。

<a id="r-6-6"></a>

**補題 6.6.** $`X,Y\in M_{p\times n}`$ かつ $`XX^*=YY^*`$ ならば、あるユニタリ $`W\in M_n`$ により $`Y=XW`$ である。

*証明.* すべての $`z`$ に対し $`\|X^*z\|^2=\langle z,XX^*z\rangle=\|Y^*z\|^2`$ であるから、$`X^*z\mapsto Y^*z`$ は
$`X^*`$ の値域から $`Y^*`$ の値域の上へのwell-definedな線形等長写像である。定理 1.8 の証明と同様に、
これは $`UX^*=Y^*`$ を満たすユニタリ $`U`$ に拡張される。$`W=U^*`$ とすればよい。$`\square`$

<a id="r-6-7"></a>

**定理 6.7（ユニタリ自由度）.** $`K_1,\dots,K_n`$ と $`L_1,\dots,L_n`$ を $`M_{k'\times k}`$ の元とする
（短い方の族は零で補う）。このとき、すべての $`A`$ に対し $`\sum_lK_lAK_l^*=\sum_jL_jAL_j^*`$ となるための必要十分条件は、
あるユニタリ $`u\in M_n`$ により $`L_j=\sum_lu_{jl}K_l`$ となることである。

*証明.* $`L_j=\sum_lu_{jl}K_l`$ ならば
$`\sum_jL_jAL_j^*=\sum_{l,l'}(u^*u)_{l'l}K_lAK_{l'}^*=\sum_lK_lAK_l^*`$ である。逆に、等しい写像は
等しいチョイ行列をもつ：$`v_l=\sum_ie_i\otimes K_le_i`$、
$`w_j=\sum_ie_i\otimes L_je_i`$ として $`\sum_lv_lv_l^*=\sum_jw_jw_j^*`$ である。$`v_l`$ と $`w_j`$ を列とする行列に補題 6.6 を適用すると
$`w_j=\sum_lW_{lj}v_l`$ を得る。$`K\mapsto\sum_ie_i\otimes Ke_i`$ は線形かつ単射であるから
$`L_j=\sum_lW_{lj}K_l`$ であり、$`u=W^T`$ はユニタリである。$`\square`$

*計算：* $`n=3`$ での定理 6.7 — [チュートリアル、6.2節](../../tutorial/ja/06_quantum_channels.md#s-6-2)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. あらゆるチャネルは、系と固定された状態にある環境との等長的な相互作用の後に
> 環境を無視することとして得られる。環境を基底 $`e_l`$ で測定すると、確率 $`\mathrm{Tr}\,K_l\rho K_l^*`$ で
> 結果 $`l`$ が得られ、状態は $`K_l\rho K_l^*/\mathrm{Tr}(K_l\rho K_l^*)`$ となる：クラウス作用素は起こりうる「量子ジャンプ」である。
> 別の基底をとれば別のジャンプが得られるが、チャネルは同じである（定理 6.7）。また最小の環境の次元は、
> 定理 6.5 の後で示したとおり、クラウス階数に等しい。

<a id="s-6-3"></a>

## 6.3 トレースに関する双対

**定義.** $`\mathrm{tr}_w`$、$`\mathrm{tr}_v`$ をそれぞれ $`\mathcal M`$、
$`\mathcal N`$ 上の忠実なトレースとする。線形写像 $`\Phi:\mathcal M\to\mathcal N`$ の $`(w,v)`$ に関する **双対** とは、
```math
\mathrm{tr}_v\big(\Phi(x)\,y\big)=\mathrm{tr}_w\big(x\,\Phi^\sharp(y)\big)\qquad(x\in\mathcal M,\ y\in\mathcal N).
```
を満たす線形写像 $`\Phi^\sharp:\mathcal N\to\mathcal M`$ のことである。
これは存在して一意である：$`y`$ を固定すると $`x\mapsto\mathrm{tr}_v(\Phi(x)y)`$ は線形汎関数であるから、
定理 3.2 により、ただ一つの $`\rho=:\Phi^\sharp(y)`$ に対して $`\mathrm{tr}_w(\rho x)=\mathrm{tr}_w(x\rho)`$ に等しく、
一意性により $`\Phi^\sharp`$ は線形となる。blunt トレースの場合には $`\Phi^*`$ と書く。
$`z_w=\sum_cw_cP_c`$、$`z_v=\sum_dv_dQ_d`$ とおく（中心に属し、可逆である）。

<a id="r-6-8"></a>

**命題 6.8.** $`\Phi`$ がクラウス作用素 $`K^{dc}_l`$ をもつならば（命題 6.4）、
```math
\Phi^\sharp(y)_c=\sum_d\frac{v_d}{w_c}\sum_lK^{dc*}_l\,y_d\,K^{dc}_l,\qquad\text{i.e.}\qquad
\Phi^\sharp(y)=z_w^{-1}\,\Phi^*(z_v\,y).
```
特に $`\Phi^\sharp`$ は CP であり、クラウス作用素 $`(v_d/w_c)^{1/2}K^{dc*}_l`$ をもつ。

*証明.* 各ブロックにおけるトレースの巡回性により
```math
\mathrm{tr}_v(\Phi(x)y)=\sum_dv_d\sum_{c,l}\mathrm{Tr}\big(K^{dc}_lx_cK^{dc*}_ly_d\big)
=\sum_cw_c\mathrm{Tr}\Big(x_c\sum_d\frac{v_d}{w_c}\sum_lK^{dc*}_ly_dK^{dc}_l\Big),
```
であり、一意性からこの式が得られる。$`v=w=1`$ のときこれは $`\Phi^*(y)_c=\sum_{d,l}K^{dc*}_ly_dK^{dc}_l`$ となる。
因子 $`v_d`$ と $`1/w_c`$ は $`z_v`$ と $`z_w^{-1}`$ を掛けることにあたる。$`\square`$

すべての $`x`$ に対し $`\mathrm{tr}_v\,\Phi(x)=\mathrm{tr}_w(x)`$ のとき、$`\Phi`$ は **$`(w,v)`$ に関してトレース保存** であるという。

<a id="r-6-9"></a>

**系 6.9.** (a) $`\Phi`$ が $`(w,v)`$ に関して TP であるための必要十分条件は $`\Phi^\sharp(1)=1`$、すなわち
すべての $`c`$ に対し $`\sum_{d,l}(v_d/w_c)K^{dc*}_lK^{dc}_l=1_{k_c}`$ となることである。(b) $`\Phi`$ が単位的であるための必要十分条件は、$`\Phi^\sharp`$ が
$`(v,w)`$ に関して TP であることである。(c) $`\Phi^\sharp`$ の $`(v,w)`$ に関する双対は $`\Phi`$ である。(d)
$`\Phi:\mathcal M\to\mathcal M`$ がセクター保存で $`v=w`$ ならば、すべての $`w`$ に対し $`\Phi^\sharp=\Phi^*`$ であり、
$`\Phi`$ が一つのトレースに関して TP であることと、すべてのトレースに関して TP であることは同値である。

*証明.* (a) $`\mathrm{tr}_v\,\Phi(x)=\mathrm{tr}_w(x\Phi^\sharp(1))`$ であり、
すべての $`x`$ に対し $`\mathrm{tr}_w(x(\Phi^\sharp(1)-1))=0`$ となるのは $`\Phi^\sharp(1)=1`$ のときかつそのときに限る（
$`x=(\Phi^\sharp(1)-1)^*`$ ととる。$`\mathrm{tr}_w`$ は忠実である）。あとは命題 6.8 を用いる。(b)
$`\mathrm{tr}_w\,\Phi^\sharp(y)=\mathrm{tr}_v(\Phi(1)y)`$ がすべての
$`y`$ に対し $`\mathrm{tr}_v(y)`$ に等しいのは $`\Phi(1)=1`$ のときかつそのときに限る。(c) トレースの性質により
$`\mathrm{tr}_w(\Phi^\sharp(y)x)=\mathrm{tr}_v(y\Phi(x))`$ である。一意性を用いる。(d) 現れるのは $`d=c`$ の組のみであり、そこでは $`v_d/w_c=1`$ である。$`\square`$

**ハイゼンベルク描像における密度.** $`\Psi:\mathcal N\to\mathcal M`$ を線形写像、$`\omega`$ を $`\mathrm{tr}_w`$ に関する
密度 $`\rho`$ をもつ $`\mathcal M`$ 上の線形汎関数（定理 3.2）、$`\Psi^\sharp:\mathcal M\to\mathcal N`$ を $`(v,w)`$ に関する
$`\Psi`$ の双対とする。このとき
$`\omega(\Psi(y))=\mathrm{tr}_w(\Psi(y)\rho)=\mathrm{tr}_v(y\Psi^\sharp(\rho))=\mathrm{tr}_v(\Psi^\sharp(\rho)y)`$
であるから、$`\omega\circ\Psi`$ は $`\mathrm{tr}_v`$ に関して密度 $`\Psi^\sharp(\rho)`$ をもつ。$`\Psi`$ が単位的ならば、
系 6.9(b)（$`\mathcal M`$ と $`\mathcal N`$ の役割を入れ替えたもの）により、$`\Psi^\sharp`$ は $`(w,v)`$ に関して TP である。

セクターを変える写像に対しては系 6.9(d) は成り立たない。$`\mathcal M=\mathcal N=M_1\oplus M_2`$（セクター $`0`$ と $`1`$）とし、
$`\Phi`$ を **粒子損失** と呼ばれる写像（下の囲みを参照）
$`\Phi(x)=(x_0+\mathrm{Tr}\,x_1)\oplus0`$ とする。そのクラウス作用素は $`K^{00}=1`$、$`K^{01}_l=\langle l|`$ である。
これは blunt トレースに関して TP である。$`w=v=(1,\tfrac12)`$ の $`\mathrm{tr}_w`$ に対しては $`\Phi^\sharp(1)=1\oplus2\cdot1_2`$
であり、$`\mathrm{tr}_v\,\Phi(x)=x_0+\mathrm{Tr}\,x_1\ne x_0+\tfrac12\mathrm{Tr}\,x_1=\mathrm{tr}_w(x)`$ である。

*ライブラリでは：* `Channel.adjoint()` が $`\Phi^*`$、`InterSectorChannel.adjoint(trace=...)` が $`\Phi^\sharp`$ である。命題 6.8、系 6.9(a)、粒子損失の検証 — [チュートリアル、6.3節](../../tutorial/ja/06_quantum_channels.md#s-6-3)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. トレースに依存しない対象は、オブザーバブルに作用するハイゼンベルク写像
> $`\Psi:\mathcal N\to\mathcal M`$ である：過程の前の状態 $`\omega`$ は過程の後には
> $`\omega\circ\Psi`$ となり、$`\Psi`$ は CP かつ単位的である。密度はトレースに依存し
> （定理 3.2）、密度は選んだトレースに関する $`\Psi`$ の双対によって変換される（系 6.9 の後の段落）。環が
> 重複度をもって作用する場合、物理的なトレースの重みは $`w_c=m_c`$ である（2.3節）。上の例では、セクター $`0`$ を
> 「粒子なし」、セクター $`1`$ を「二つのサイト上の一つの粒子」と読む。すると写像は粒子を取り除く。これは状態を数える
> blunt トレースに関しては確率を保存するが、セクターを等しく重みづける `Tr_norm` に関しては保存しない。

<a id="s-6-4"></a>

## 6.4 デコヒーレンス

ライブラリの標準的なチャネルは各ブロックの内部で作用する（$`p,\gamma\in[0,1]`$）：

| チャネル | $`\Phi(\rho)_c`$ | ブロック $`c`$ のクラウス作用素 |
|---|---|---|
| `dephasing_channel` | $`(1-p)\rho_c+p\mathrm{diag}\,\rho_c`$ | $`\sqrt{1-p}\,1`$, $`\sqrt p\,E_{ii}`$ |
| `depolarizing_channel` | $`(1-p)\rho_c+p\mathrm{Tr}(\rho_c)\,1/k_c`$ | $`\sqrt{1-p}\,1`$, $`\sqrt{p/k_c}\,E_{ij}`$ |
| `center_expectation` | $`\mathrm{Tr}(\rho_c)\,1/k_c`$ | $`\sqrt{1/k_c}\,E_{ij}`$ |
| `amplitude_damping_channel` | 準位 $`i>0`$ は $`0`$ へ減衰する | $`E_{00}+\sqrt{1-\gamma}\sum_{i>0}E_{ii}`$, $`\sqrt\gamma\,E_{0i}`$ |
| `random_mixed_unitary_channel` | $`\sum_jp_jU_j\rho U_j^*`$ | $`\sqrt{p_j}\,U_j`$ |

クラウス形は直接計算で確かめられる。例えば
$`\sum_{ij}\frac1kE_{ij}\rho E_{ji}=\frac1k\sum_{ij}\rho_{jj}E_{ii}`$ である。五つとも TP であり、
$`\gamma>0`$ の振幅減衰（ある $`k_c\ge2`$ のとき）を除いてすべて単位的である（命題 6.4：$`\sum_lK_l^*K_l`$ と
$`\sum_lK_lK_l^*`$ を計算する。サイズ $`k`$ のブロックにおける振幅減衰では、それらは $`1`$ と
$`1+(k-1)\gamma E_{00}-\gamma\sum_{i>0}E_{ii}`$ であり、混合ユニタリ・チャネルでは $`\sum_jp_j=1`$ である）。
具体的には、振幅減衰は $`\rho_{00}\mapsto\rho_{00}+\gamma\sum_{i>0}\rho_{ii}`$、
$`\rho_{ii}\mapsto(1-\gamma)\rho_{ii}`$、$`i\ne j`$ がともに正のとき $`\rho_{ij}\mapsto(1-\gamma)\rho_{ij}`$、
$`j>0`$ に対し $`\rho_{0j}\mapsto\sqrt{1-\gamma}\,\rho_{0j}`$、$`\rho_{j0}\mapsto\sqrt{1-\gamma}\,\rho_{j0}`$ と写す。中心の上への条件付き期待値 $`E`$
（命題 2.8）は $`p=1`$ の脱分極である。セクターの測定
$`X\mapsto\sum_cP_cXP_c`$ は $`\mathcal M`$ 上で恒等写像である（3.3節）のと異なり、$`E`$ は $`p_c`$ だけを残す。

*ライブラリでは：* 表の五つのチャネルは `channels` の関数である — [チュートリアル、6.4節](../../tutorial/ja/06_quantum_channels.md#s-6-4)。

<a id="r-6-10"></a>

**定理 6.10（単位的チャネルの不動点）.** $`K_l\in\mathcal M`$ として $`\Phi(X)=\sum_lK_lXK_l^*`$
が TP かつ単位的であるとする。このとき $`\Phi(X)=X`$ となるための必要十分条件は、$`X`$ がすべての $`K_l`$ および $`K_l^*`$ と可換であることである。特に
不動点全体は $`\mathcal M`$ のフォン・ノイマン部分環をなす。

*証明.* $`X`$ がすべての $`K_l`$ と可換ならば $`\Phi(X)=X\sum_lK_lK_l^*=X`$ である。逆に $`\Phi(X^*)=\Phi(X)^*`$
であるから、不動点の自己共役部分 $`(X+X^*)/2`$ と $`(X-X^*)/2i`$ も不動点である。そこで
$`X=X^*`$ とする。$`Y_l=K_lX-XK_l`$ とおくと
```math
\sum_lY_lY_l^*=\sum_l\big(K_lX^2K_l^*-K_lXK_l^*X-XK_lXK_l^*+XK_lK_l^*X\big)
=\Phi(X^2)-\Phi(X)X-X\Phi(X)+X\Phi(1)X=\Phi(X^2)-X^2,
```
となる。ここで $`Y_l^*=XK_l^*-K_l^*X`$、$`\Phi(X)=X`$、$`\Phi(1)=1`$ を用いた。$`\Phi`$ が TP であるからそのトレースは $`0`$ である。
$`\mathrm{Tr}\,Y_lY_l^*=\|Y_l\|_2^2\ge0`$ であるから、すべての $`Y_l=0`$ である。
共役をとれば、$`X`$ は $`K_l^*`$ とも可換である。このような $`X`$ は自己共役な集合の可換子環と $`\mathcal M`$ との共通部分をなす：
これは $`1`$ を含む $`*`$-環であり、したがってフォン・ノイマン環である
（定理 2.2）。$`\square`$

したがって、位相緩和（$`p>0`$）は各ブロックの対角部分をちょうど保ち、脱分極（$`p>0`$）と $`E`$ は
ちょうど中心を保つ。単位性がなければ定理は成り立たない：振幅減衰は $`E_{00}`$ を固定するが、
$`E_{00}`$ はジャンプ $`E_{01}`$ と可換でない。$`0<\gamma\le1`$ かつ $`k\ge2`$ のとき、ブロックにおけるその不動点の空間と
クラウス作用素の可換子環はいずれも次元 $`1`$ であるが、それぞれ $`\mathbb CE_{00}`$ と
$`\mathbb C1`$ である：上の具体形により、不動点は $`i\ne j`$ に対し $`\rho_{ij}=0`$、$`i>0`$ に対し
$`\rho_{ii}=0`$ をみたし（これらの成分には $`\sqrt{1-\gamma}\ne1`$ または $`1-\gamma\ne1`$ が掛かる）、$`E_{00}`$ の
任意のスカラー倍は不動点である。また、すべての $`E_{0i}`$、$`E_{i0}`$（$`i>0`$）と可換なものは、それらの積
$`E_{ij}=E_{i0}E_{0j}`$ および $`E_{00}=E_{01}E_{10}`$ とも可換であり、したがって $`M_k`$ と可換であるから、$`1`$ のスカラー倍である
（第2章の演習問題 1）。

*計算：* 五つのチャネルの不動点と可換子環 — [チュートリアル、6.4節](../../tutorial/ja/06_quantum_channels.md#s-6-4)。

<a id="r-6-11"></a>

**定理 6.11（単位的チャネルはエントロピーを増加させる）.** $`\Phi:\mathcal M\to\mathcal N`$ を正値、
TP かつ単位的とし、$`\rho`$ と $`\Phi(\rho)`$ を固有値 $`\lambda_i`$ と $`\mu_j`$ をもつ状態とする。
このとき、ある二重確率行列 $`\Delta`$ により $`\mu=\Delta\lambda`$ となり、したがって $`[0,1]`$ 上の任意の凹関数 $`f`$ に対し
$`\sum_jf(\mu_j)\ge\sum_if(\lambda_i)`$ である。特に $`S(\Phi(\rho))\ge S(\rho)`$ かつ
$`\mathrm{Tr}\,\Phi(\rho)^2\le\mathrm{Tr}\,\rho^2`$ である。

*証明.* 各ブロックを対角化する：$`\rho=\sum_i\lambda_i|e_i\rangle\langle e_i|`$、
$`\Phi(\rho)=\sum_j\mu_j|f_j\rangle\langle f_j|`$。ここで $`\bigoplus_c\mathbb C^{k_c}`$
と $`\bigoplus_d\mathbb C^{k'_d}`$ の正規直交基底は、各ベクトルがそれぞれ一つの直和成分に属するようにとる。すると
$`|e_i\rangle\langle e_i|\in\mathcal M`$ である。このとき $`\mu_j=\langle f_j,\Phi(\rho)f_j\rangle=\sum_i\Delta_{ji}\lambda_i`$
であり、$`\Delta_{ji}=\langle f_j,\Phi(|e_i\rangle\langle e_i|)f_j\rangle\ge0`$、
$`\sum_j\Delta_{ji}=\mathrm{Tr}\,\Phi(|e_i\rangle\langle e_i|)=1`$、
$`\sum_i\Delta_{ji}=\langle f_j,\Phi(1)f_j\rangle=1`$ である（したがって $`\Delta`$ は正方行列である：全成分の和をとればわかる）。
有限形のイェンセンの不等式 $`f(\sum_i\Delta_{ji}\lambda_i)\ge\sum_i\Delta_{ji}f(\lambda_i)`$（第3章の演習問題 2 の解答と同様に、
凹性の定義から項数に関する帰納法で従う）により、
$`\sum_jf(\mu_j)\ge\sum_j\sum_i\Delta_{ji}f(\lambda_i)=\sum_if(\lambda_i)`$ である。$`f(t)=-t\log t`$
（$`f(0)=0`$）と $`f(t)=-t^2`$ をとればよい。これらは補題 3.4 により $`[0,1]`$ 上で凹である：$`t\log t`$ と $`t^2`$ は
$`[0,1]`$ 上で連続であり、$`(0,1)`$ 上での二階導関数はそれぞれ $`1/t>0`$ と $`2>0`$ である。$`\square`$

*計算：* 単位的チャネルのもとでのエントロピーと純度、および振幅減衰によるエントロピーの減少 — [チュートリアル、6.4節](../../tutorial/ja/06_quantum_channels.md#s-6-4)。

セクター間の写像の古典的部分はその **遷移行列**
$`T_{dc}(\rho)=\mathrm{Tr}\,\Phi_{dc}(\rho_c)/\mathrm{Tr}\,\rho_c`$（$`\rho_c\ne0`$）である。ライブラリの
`transition_matrix()` は、すべての入力セクターの最大混合状態に対する $`T=T(1)`$ である。

<a id="r-6-12"></a>

**命題 6.12.** $`\Phi:\mathcal M\to\mathcal N`$ を正値かつ TP とし、$`\rho`$ を
$`p_c=\mathrm{Tr}\,\rho_c>0`$ なる状態、$`q_d=\mathrm{Tr}\,\Phi(\rho)_d`$ とする。このとき $`T(\rho)`$ は列
確率行列であり、$`q=T(\rho)p`$ である。中心の上への条件付き期待値 $`E_{\mathcal M}`$、$`E_{\mathcal N}`$（命題 2.8）を用いると
```math
E_{\mathcal N}\Phi E_{\mathcal M}(\rho)=\bigoplus_d\frac{(Tp)_d}{k'_d}\,1_{k'_d},\qquad T=T(1).
```

*証明.* $`\Phi_{dc}`$ は正値であるから $`T_{dc}\ge0`$ であり、
$`\sum_dT_{dc}=\mathrm{Tr}\,\Phi(\iota_c\rho_c)/\mathrm{Tr}\,\rho_c=1`$ である。次に
$`q_d=\sum_c\mathrm{Tr}\,\Phi_{dc}(\rho_c)=\sum_cT_{dc}(\rho)p_c`$ である。最後に
$`E_{\mathcal M}\rho=\bigoplus_cp_c1_{k_c}/k_c`$ であるから、
$`\Phi(E_{\mathcal M}\rho)_d=\sum_cp_c\Phi_{dc}(1_{k_c}/k_c)`$ のトレースは $`\sum_cT_{dc}p_c`$ であり、$`E_{\mathcal N}`$
はブロック $`d`$ をそのトレースと $`1/k'_d`$ の積で置き換える。$`\square`$

*ライブラリでは：* 遷移行列は `InterSectorChannel.transition_matrix()` である。命題 6.12 の検証 — [チュートリアル、6.4節](../../tutorial/ja/06_quantum_channels.md#s-6-4)。

**量子ビットのコヒーレンス.** $`M_2`$ 上で $`D_p`$ を位相緩和、$`A_\gamma`$ を振幅減衰とする。表と上の具体形により、
$`D_p`$ は $`\rho_{01}`$ と $`\rho_{10}`$ に $`1-p`$ を掛けて対角部分を固定し、$`A_\gamma`$ は $`\rho_{01}`$ と $`\rho_{10}`$ に
$`\sqrt{1-\gamma}`$ を掛けて、対角部分をコヒーレンスとは無関係に $`(\rho_{00}+\gamma\rho_{11},(1-\gamma)\rho_{11})`$ に写す。
したがって $`D_pA_\gamma=A_\gamma D_p`$ である：どちらもコヒーレンスに $`(1-p)\sqrt{1-\gamma}`$ を掛け、対角部分には
$`A_\gamma`$ として作用する。特に、$`T_1,T_\varphi>0`$、$`t\ge0`$、$`1-p=e^{-t/T_\varphi}`$、$`1-\gamma=e^{-t/T_1}`$ のとき、
合成はコヒーレンスに $`e^{-t/T_2}`$ を掛ける。ここで
```math
\frac1{T_2}=\frac1{2T_1}+\frac1{T_\varphi},
```
である。またこれらの族は $`t`$ について加法的に合成される：$`D_{p(t)}D_{p(s)}=D_{p(t+s)}`$、
$`A_{\gamma(t)}A_{\gamma(s)}=A_{\gamma(t+s)}`$。因子 $`e^{-t/T}`$ は掛け合わされ、トレースは保存されるからである。

*計算：* コヒーレンス因子と $`1/T_2=1/(2T_1)+1/T_\varphi`$ — [チュートリアル、6.4節](../../tutorial/ja/06_quantum_channels.md#s-6-4)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. 量子ビットに対し、
> 位相緩和は時定数 $`T_\varphi`$ の純粋位相緩和を、振幅減衰は時定数 $`T_1`$ での励起準位の減衰をモデル化する。
> すると上の計算は、その両方を受ける量子ビットのコヒーレンス時間 $`T_2`$ を与える。定理 6.10 は
> デコヒーレンスのない情報の代数的な形である：すべてのジャンプ作用素と可換なものは生き残る（第10章）。
> 命題 6.12 によれば、中心を通して見ると、セクター間の写像は電荷の上の古典的なマルコフ
> 行列である。

<a id="s-6-5"></a>

## 6.5 ランダムなチャネル

**ハール・ランダムな等長作用素** $`V:\mathbb C^k\to\mathbb C^N`$ とは、ライブラリ（`random_isometry`）と同様に、
$`G=QR`$ の因子 $`Q`$ のことである。ここで $`G`$ は $`N\times k`$ のジニブル行列（成分は独立で、密度
$`\pi^{-1}e^{-|z|^2}`$ をもつ）、$`Q^*Q=1`$、$`R`$ は対角成分が正の上三角行列である（CPU 上では、ライブラリは LAPACK の分解 $`G=Q'R'`$ をとり、$`Q'`$ の第 $`i`$ 列に
$`R'_{ii}/|R'_{ii}|`$ を掛ける。定理 5.3 と同様に、これはすべての $`R'_{ii}\ne0`$ のとき、すなわち $`G`$ の階数が $`k`$ のとき
この分解を与え、補題 6.13 の証明によりこれはほとんど確実に成り立つ。GPU 上の大きなバッチに対してはグラム・シュミットの
直交化を用い、これは直接この分解を与える）。
`random_channel(alg, r)` はこのような $`V_c:\mathbb C^{k_c}\to\mathbb C^r\otimes\mathbb C^{k_c}`$ を
各ブロックで独立に生成し、定理 6.5 のクラウス作用素を用いる。

*ライブラリでは：* `random_isometry` と `channels.random_channel` — [チュートリアル、6.5節](../../tutorial/ja/06_quantum_channels.md#s-6-5)。

<a id="r-6-13"></a>

**補題 6.13.** $`V`$ は $`U(N)`$ のハール・ユニタリの最初の $`k`$ 列と同じ分布をもつ。
特に、任意のユニタリ $`U\in M_N`$ に対し、$`UV`$ は $`V`$ と同じ分布をもつ。

*証明.* $`Z`$ を最初の $`k`$ 列が $`G`$ である $`N\times N`$ のジニブル行列とし、$`E\in M_{N\times k}`$ を
$`1_N`$ の最初の $`k`$ 列とする。定理 5.3 により、$`Z`$ はほとんど確実に可逆であり、$`Z=Q_ZR_Z`$ で
$`Q_Z`$ はハール分布に従う。$`R_Z`$ は上三角であるから、$`G=ZE=(Q_ZE)(E^*R_ZE)`$ は
求める型の分解であり、しかもそれはただ一つである：$`G=QR=Q'R'`$ ならば $`R^*R=G^*G=R'^*R'`$ であるから、$`R'R^{-1}`$
はユニタリかつ対角成分が正の上三角行列であり、したがって $`1`$ である（補題 5.2 の証明）。ゆえに $`V=Q_ZE`$ であり、
$`UQ_Z`$ は $`Q_Z`$ と同じ分布をもつから、$`UV=(UQ_Z)E`$ は $`V`$ と同じ分布をもつ。$`\square`$

<a id="r-6-14"></a>

**補題 6.14.** $`u`$ を $`\mathbb C^N`$ のランダムな単位ベクトルとし、その分布はすべての
ユニタリのもとで不変であるとする。このとき $`\mathbb E[u_a\bar u_b]=\delta_{ab}/N`$ であり、
```math
\mathbb E\big[u_a\bar u_b\,u_c\bar u_d\big]=\frac{\delta_{ab}\delta_{cd}+\delta_{ad}\delta_{cb}}{N(N+1)} .
```

*証明.* $`\mathrm{diag}(e^{i\theta_1},\dots,e^{i\theta_N})`$ のもとでの不変性により、
これらのモーメントにはそれぞれ $`e^{i(\theta_a-\theta_b)}`$ と $`e^{i(\theta_a-\theta_b+\theta_c-\theta_d)}`$ が掛かるので、
$`a=b`$、あるいはそれぞれ多重集合として $`\{a,c\}=\{b,d\}`$ でない限り零となる。置換行列のもとでの不変性により
$`\mathbb E|u_a|^2=1/N`$ であり、残るのは $`\alpha=\mathbb E|u_1|^4`$ と $`\beta=\mathbb E|u_1|^2|u_2|^2`$ だけである。
第一に、$`\sum_c|u_c|^2=1`$ から $`\alpha+(N-1)\beta=\mathbb E|u_1|^2=1/N`$ を得る。第二に、$`(u_1+u_2)/\sqrt2`$ は
あるユニタリ $`U`$ に対する $`Uu`$ の第一座標であるから、$`\mathbb E|u_1+u_2|^4=4\alpha`$ である。
$`s=|u_1|^2+|u_2|^2`$、$`z=u_1\bar u_2`$ とおくと $`|u_1+u_2|^4=s^2+2s(z+\bar z)+(z+\bar z)^2`$ である。
不変性 $`u_1\mapsto-u_1`$ と $`u_1\mapsto iu_1`$ から $`\mathbb E[s(z+\bar z)]=0`$ と
$`\mathbb Ez^2=0`$ が得られるので、$`4\alpha=(2\alpha+2\beta)+2\beta`$、すなわち $`\alpha=2\beta`$ である。したがって
$`\beta=1/(N(N+1))`$、$`\alpha=2\beta`$ であり、これはすべての場合に上の式を与える。（ここでは $`N\ge2`$ を用いた。
$`N=1`$ のときは $`|u_1|=1`$ であり、両辺とも $`1`$ に等しい。）$`\square`$

<a id="r-6-15"></a>

**定理 6.15（ランダムなチャネルの出力純度）.** $`\psi`$ をサイズ
$`k`$ のセクターにおける純粋状態とし、$`\Phi`$ を `random_channel` でクラウス階数 $`r`$ として生成する。このとき
```math
\mathbb E\mathrm{Tr}\,\Phi(\psi)^2=\frac{k+r}{kr+1},\qquad\mathbb E\,\Phi(\psi)=\frac1k1_k .
```

*証明.* $`\psi`$ のセクターにおいて $`\Phi(\psi)=\mathrm{Tr}_E|u\rangle\langle u|`$ であり、$`u=V\psi`$ は
単位ベクトルで、補題 6.13 によりその分布はユニタリ不変である。$`\mathbb C^r\otimes\mathbb C^k`$ の添字を
組 $`(e,x)`$ で表す。すると $`\Phi(\psi)_{xy}=\sum_eu_{ex}\bar u_{ey}`$ であり、その平均は $`r\delta_{xy}/(kr)`$ である。また
$`\mathrm{Tr}\,\Phi(\psi)^2=\sum_{x,y,e,f}u_{ex}\bar u_{ey}u_{fy}\bar u_{fx}`$ である。補題 6.14 を
$`N=kr`$ として用いると、第一項は $`x=y`$（添字の組は $`kr^2`$ 個）を、第二項は $`e=f`$（$`k^2r`$ 個）を強いるので、
平均は $`(kr^2+k^2r)/(kr(kr+1))`$ である。$`\square`$

*計算：* $`r=1,2,4`$ に対する定理 6.15 のモンテカルロ法による検証 — [チュートリアル、6.5節](../../tutorial/ja/06_quantum_channels.md#s-6-5)。

$`r=1`$ のときチャネルはユニタリであり、純度は $`1`$ である。$`r\to\infty`$ のとき純度は
最大混合状態の純度 $`1/k`$ に近づく。証明は $`u`$ の分布のユニタリ不変性しか用いていないので、
ランダムな純粋状態 $`u=g/\|g\|`$（$`g`$ はガウス的であり、その分布は定理 A.20(b) によりユニタリ不変である）の
縮約状態の平均純度に関するルブキンの結果も与える。これは第7章の定理 7.11 でガウス積分を用いて証明される。補題 6.14 の数
$`\alpha,\beta`$ は、ハール・ユニタリの成分に関する系 5.11 の数と同じであり、系 A.21 は $`u=g/\|g\|`$ に対して
それらを計算し直す。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. ランダムなチャネルは、次元 $`r`$ の環境に由来する
> 一般的で構造のない雑音のモデルである。すでに $`r=4`$ で、純粋なクトリット状態の純度は平均して
> $`7/13\approx0.54`$ まで下がる（定理 6.15）。雑音がなければ $`1`$、最小値は $`1/3`$ である（コーシー・シュワルツの不等式により
> $`\mathrm{Tr}\,\rho^2\ge(\mathrm{Tr}\,\rho)^2/k`$）。超選択された系では、雑音は各セクターで
> 別々に作用し、問題となる次元は全空間の次元ではなくセクターの次元である。

<a id="s-6-6"></a>

## 6.6 比較の対象

*ライブラリでは：* 他のソフトウェアとの比較 — [チュートリアル、6.6節](../../tutorial/ja/06_quantum_channels.md#s-6-6)。

<a id="x-6"></a>

## 演習問題

1. （ワーナー–ホレボ）$`M_k`$（$`k\ge2`$）上の $`\Phi(\rho)=\big(\mathrm{Tr}(\rho)1-\rho^T\big)/(k-1)`$ は、
   転置から作られているにもかかわらず、クラウス階数 $`k(k-1)/2`$ の単位的チャネルであることを示せ。
2. $`0<p\le1`$ に対し、$`M_k`$ 上の位相緩和チャネルのクラウス階数は $`k`$ であり、ライブラリが格納する
   クラウス作用素の個数より一つ少ないことを示せ。
3. $`X`$ を $`\mathbb C^k`$ 上の巡回シフト、$`Z=\mathrm{diag}(1,\omega,\dots,\omega^{k-1})`$、
   $`\omega=e^{2\pi i/k}`$ とする。$`\frac1{k^2}\sum_{a,b=0}^{k-1}X^aZ^b\rho Z^{-b}X^{-a}=\mathrm{Tr}(\rho)1/k`$ を示せ：
   $`M_k`$ の中心の上への条件付き期待値は混合ユニタリ・チャネルである。
4. 中間の環で同じトレースを用いるとき $`(\Psi\circ\Phi)^\sharp=\Phi^\sharp\circ\Psi^\sharp`$ であることを示せ。

*計算：* 四つの演習問題の数値的検証 — [チュートリアル、演習問題](../../tutorial/ja/06_quantum_channels.md#x-6)。

<details>
<summary>解答</summary>

1. $`F`$ を入れ替え作用素とすると、$`J(\theta)=\sum_{ij}E_{ij}\otimes E_{ji}=F`$、$`J(\rho\mapsto\mathrm{Tr}(\rho)1)=1\otimes1`$
   であるから、$`J(\Phi)=(1-F)/(k-1)=2P_-/(k-1)`$ である。ここで $`P_-=(1-F)/2`$ は
   反対称ベクトルの上への射影であり（$`F=F^*`$、$`F^2=1`$）、その階数は $`\mathrm{Tr}\,P_-=(k^2-\mathrm{Tr}\,F)/2`$、
   ただし $`\mathrm{Tr}\,F=\sum_{ij}\langle e_i\otimes e_j,e_j\otimes e_i\rangle=k`$ である。よって $`J\ge0`$ で階数は $`k(k-1)/2`$ であり（定理 6.2、
   系 6.3）、$`\mathrm{Tr}\,\Phi(\rho)=\mathrm{Tr}\,\rho`$、$`\Phi(1)=1`$ は直接確かめられる。

2. $`J=(1-p)|\Omega\rangle\langle\Omega|+p\sum_iE_{ii}\otimes E_{ii}`$ は
   $`e_i\otimes e_i`$ の張る空間の外では零であり、この $`k`$ 次元空間上では行列 $`(1-p)\mathbf 1\mathbf 1^T+p1`$
   （$`\mathbf 1`$ はすべての成分が $`1`$ のベクトル）であって、その固有値は $`p`$ と $`(1-p)k+p`$ で、すべて正である。よって
   $`\mathrm{rank}\,J=k`$ であるが、ライブラリは $`1+k_{\max}`$ 個の作用素を格納する。

3. $`\mathrm{Tr}((X^aZ^b)^*X^{a'}Z^{b'})=k\,\delta_{aa'}\delta_{bb'}`$ である：$`a=a'`$ でなければ
   $`Z^{-b}X^{a'-a}Z^{b'}`$ の対角成分は零であり、$`\mathrm{Tr}\,Z^{b'-b}=\sum_j\omega^{j(b'-b)}=k\delta_{bb'}`$ である。よって
   $`B_{ab}=X^aZ^b/\sqrt k`$ は $`\mathrm{Tr}(A^*B)`$ に関する $`M_k`$ の正規直交基底であり、任意の
   二つの正規直交基底はユニタリで結ばれる。定理 6.7 により、$`\sum_mB_m\rho B_m^*`$ はすべての正規直交基底について
   同じである。基底 $`E_{ab}`$ に対しては
   $`\sum_{ab}E_{ab}\rho E_{ba}=\mathrm{Tr}(\rho)1`$ である。$`k`$ で割ればよい。

4. すべての $`x`$ と $`z`$ に対し $`\mathrm{tr}(\Psi(\Phi(x))z)=\mathrm{tr}(\Phi(x)\Psi^\sharp(z))=\mathrm{tr}(x\,\Phi^\sharp(\Psi^\sharp(z)))`$
   である。双対の一意性から主張が従う。

</details>

**次へ：**[第7章](07_tensor_products.md)ではテンソル積とエンタングルメントを扱う：部分トレース——
異なる環の間のチャネルであり、その双対は埋め込み $`A\mapsto A\otimes1`$ である——、ペイジの
公式、そして電荷の融合。
