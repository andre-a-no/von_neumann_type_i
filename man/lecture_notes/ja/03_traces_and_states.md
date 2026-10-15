[English](../03_traces_and_states.md) · [Español](../es/03_traces_and_states.md) · [Français](../fr/03_traces_and_states.md) · [Deutsch](../de/03_traces_and_states.md) · [中文](../zh/03_traces_and_states.md) · **日本語** · [Русский](../ru/03_traces_and_states.md)

# 3. トレース、状態、エントロピー

$`M_n(\mathbb C)`$ 上のトレースは本質的に一つしかない。$`\bigoplus_cM_{k_c}`$ 上には
$`C`$ 個のパラメータをもつ族があり、その選択は重要である：密度、エントロピー、写像の双対はすべてそれに依存する。本章ではトレースを分類し、状態を密度で表し、基本的なエントロピー不等式を証明し、最後にもっともらしいが偽であるトレース不等式を取り上げる。

**本章で学ぶこと**：$`\bigoplus_cM_{k_c}`$ 上のトレースの分類と、異なるトレースに関する状態の密度の間の関係、状態の古典的部分、フォン・ノイマン・エントロピーのセクターへの分解、クラインの不等式とその帰結、そしてマイケルソン・コントラスト。本章の数値的な検算は
[チュートリアル 3](../../tutorial/ja/03_traces_and_states.md) にある。

本章を通じて、$`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$（すべての $`k_c\ge1`$）、
$`D=\sum_ck_c`$ とし、$`P_c`$ は極小中心射影、$`A=\bigoplus_cA_c`$ とする。

<a id="s-3-1"></a>

## 3.1 トレース

$`\mathcal M`$ 上の **トレース** とは、すべての
$`A,B`$ に対し $`\tau(AB)=\tau(BA)`$ を満たす線形汎関数 $`\tau`$ のことである。すべての $`A`$ に対し $`\tau(A^*A)\ge0`$ のとき **正** であるといい、さらに
$`\tau(A^*A)=0`$ ならば $`A=0`$ となるとき **忠実** であるという。

<a id="r-3-1"></a>

**定理 3.1.** $`\mathcal M`$ 上の任意のトレースは、$`w\in\mathbb C^C`$ により
$`\mathrm{tr}_w(A)=\sum_cw_c\mathrm{Tr}\,A_c`$ の形をもつ。それが正であるための必要十分条件はすべての $`w_c\ge0`$ であり、忠実であるための必要十分条件はすべての $`w_c>0`$ である。

*証明.* $`E^c_{ij}`$ を第 $`c`$ ブロックの行列単位とする。トレース $`\tau`$ に対し
$`\tau(E^c_{ij})=\tau(E^c_{i1}E^c_{1j})=\tau(E^c_{1j}E^c_{i1})=\delta_{ij}\tau(E^c_{11})`$ である。
$`w_c=\tau(E^c_{11})`$ とおくと、線形性により $`\tau(A)=\sum_c\sum_i(A_c)_{ii}w_c`$ である。逆に
任意の $`\mathrm{tr}_w`$ はトレースである。正値性：$`E^c_{11}=(E^c_{11})^*E^c_{11}`$ であるから
$`w_c=\tau\big((E^c_{11})^*E^c_{11}\big)\ge0`$ は必要であり、$`E^c_{11}\neq0`$ であるから忠実性のためには $`w_c>0`$ が必要である。
逆に $`\mathrm{tr}_w(A^*A)=\sum_cw_c\|A_c\|_2^2`$ であり、これはすべての $`w_c\ge0`$ ならば $`\ge0`$ であり、
すべての $`w_c>0`$ ならば $`A=0`$ のときに限り零となる。$`\square`$

ライブラリは三つの忠実なトレースを実装している。

| 名称 | 重み $`w_c`$ | $`1`$ での値 | 意味 |
|---|---|---|---|
| `Tr_blunt` | $`1`$ | $`D`$ | $`\bigoplus_c\mathbb C^{k_c}`$ のトレース（重複度 $`m_c=1`$） |
| `Tr_norm` | $`1/k_c`$ | $`C`$ | 各ブロックを正規化したもの |
| `tau_vN` | $`1/(Ck_c)`$ | $`1`$ | セクターの重みが等しく $`1/C`$ であるトレース状態 |

2.3節により、$`\mathcal M`$ が重複度 $`m_c`$ で作用する空間のトレースは
$`w_c=m_c`$ とした $`\mathrm{tr}_w`$ である。したがって忠実なトレースとは「各セクターをどれだけ数えるか」の選択である。

*ライブラリでは：* トレース `Tr_blunt`, `Tr_norm`, `tau_vN`、その $`1`$ での値、および定理 3.1 のトレース性 — [チュートリアル、3.1節](../../tutorial/ja/03_traces_and_states.md#s-3-1)。

<a id="s-3-2"></a>

## 3.2 状態と密度

**状態** とは、$`\mathcal M`$ 上の線形汎関数 $`\omega`$ であって、正
（$`\omega(A^*A)\ge0`$）かつ正規化されている（$`\omega(1)=1`$）もののことである。

<a id="r-3-2"></a>

**定理 3.2.** $`\mathrm{tr}_w`$ を忠実なトレースとする。任意の線形汎関数
$`\omega`$ に対し、すべての $`A`$ について $`\omega(A)=\mathrm{tr}_w(\rho_wA)`$ を満たす $`\rho_w\in\mathcal M`$ がただ一つ存在する。汎関数 $`\omega`$ が正であるための必要十分条件は $`\rho_w\ge0`$ であり、状態であるための必要十分条件はさらに
$`\mathrm{tr}_w(\rho_w)=1`$ となることである。二つの忠実なトレースに対し、密度はブロックごとに次の関係にある。
```math
w_c\,(\rho_w)_c=w'_c\,(\rho_{w'})_c .
```

*証明.* $`\langle X,Y\rangle_w=\mathrm{tr}_w(X^*Y)`$ は $`\mathcal M`$ 上の内積である
（定理 3.1）。有限次元空間
$`\mathcal M`$ におけるリースの表現定理により、$`\omega(A)=\langle R,A\rangle_w=\mathrm{tr}_w(R^*A)`$ を満たす $`R`$ がただ一つ存在する。
$`\rho_w=R^*`$ とおく。$`\rho_w\ge0`$ ならば
$`\omega(A^*A)=\sum_cw_c\mathrm{Tr}(A_c\,(\rho_w)_c\,A_c^*)\ge0`$ である。実際
$`A_c(\rho_w)_cA_c^*\ge0`$ であり（系 1.7 を $`T=(\rho_w)_c^{1/2}A_c^*`$ として適用）、正の行列のトレースは
$`\ge0`$ であり（その固有値は $`\ge0`$、系 1.4(3)）、さらに $`w_c>0`$ だからである。逆に
$`\omega`$ が正であるとし、$`x\in\mathbb C^{k_c}`$（$`x\neq0`$）をとり、$`A`$ をブロック $`c`$ における階数1の作用素
$`|x\rangle\langle x|`$ とする。すると $`A=A^*A/\|x\|^2`$ であり
$`0\le\omega(A)=w_c\langle x,(\rho_w)_cx\rangle`$ であるから $`(\rho_w)_c\ge0`$ である。正規化は
$`\omega(1)=\mathrm{tr}_w(\rho_w)`$ である。最後に、$`\sum_cw_c\mathrm{Tr}((\rho_w)_cA_c)`$ と
$`\sum_cw'_c\mathrm{Tr}((\rho_{w'})_cA_c)`$ はともにすべての $`A`$ に対して $`\omega(A)`$ に等しいから、$`\bigoplus_cw_c(\rho_w)_c`$ と
$`\bigoplus_cw'_c(\rho_{w'})_c`$ はいずれも blunt トレース（$`w_c=1`$）に関する $`\omega`$ の密度であり、一意性により
$`w_c(\rho_w)_c=w'_c(\rho_{w'})_c`$ である。$`\square`$

*ライブラリでは：* 状態は `DensityMatrix` であり、その blunt 密度によって格納される。定理 3.2 の変換は `DensityMatrix.density` と `DensityMatrix.from_density` である — [チュートリアル、3.2節](../../tutorial/ja/03_traces_and_states.md#s-3-2)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. $`\omega(A)`$ はオブザーバブル $`A`$ の期待値である。
> blunt トレースに関する密度 $`\rho`$ は、$`\bigoplus_c\mathbb C^{k_c}`$ 上の量子力学における通常の密度行列である。他のトレースに関する密度は、
> 同じ環が重複度をもって作用するとき（第2章）、たとえばスピン系の
> 対称部分と非対称部分の上で現れる。

<a id="s-3-3"></a>

## 3.3 状態の古典的部分

$`\omega`$ の中心への制限はセクター上の確率分布
$`p_c=\omega(P_c)=\mathrm{Tr}\,\rho_c`$ である。$`p_c>0`$ ならば、セクター $`c`$ における **条件付き状態** は
$`\mathbb C^{k_c}`$ 上の密度 $`\sigma_c=\rho_c/p_c`$ をもち、
```math
\rho=\bigoplus_cp_c\,\sigma_c
```
である（$`p_c=0`$ ならば $`\rho_c=0`$ である。実際、$`\rho_c\ge0`$ の固有値は $`\ge0`$ でその和は
$`\mathrm{Tr}\,\rho_c=0`$ であるから（系 1.4(3)）すべて $`0`$ であり、$`\rho_c`$ はユニタリ対角化可能である
（定理 1.3）。その項は省く）。
したがって $`\mathcal M`$ 上の任意の状態はセクターの状態の凸結合である。写像
$`\rho\mapsto\sum_cP_c\rho P_c`$ は $`\mathcal M`$ 上で恒等写像である。なぜなら
$`\mathcal M`$ の任意の元はブロック対角だからである。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. $`\mathcal M`$ 上の
> 状態は、各セクターにおける量子状態の古典的な混合である。写像
> $`\rho\mapsto\sum_cP_c\rho P_c`$ は「どのセクターか」の測定であり、$`\mathcal M`$ のいかなる状態も
> 乱さない。これが超選択則の操作的な内容である。

*計算：* セクター確率 $`p_c`$（`DensityMatrix.sector_probabilities`）と、$`\rho`$ を乱さない測定 $`\rho\mapsto\sum_cP_c\rho P_c`$（`channels.measurement_channel`） — [チュートリアル、3.3節](../../tutorial/ja/03_traces_and_states.md#s-3-3)。

<a id="s-3-4"></a>

## 3.4 エントロピー

（blunt）密度 $`\rho`$ をもつ状態の **フォン・ノイマン・エントロピー** は
$`S(\rho)=-\mathrm{Tr}\,\rho\log\rho=-\sum_i\lambda_i\log\lambda_i`$ である。ここで和は
$`\rho`$ の固有値にわたり、$`0\log0=0`$ とする。すなわち、$`f(t)=t\log t`$, $`f(0)=0`$ に対して $`\rho\log\rho:=f(\rho)`$ であり、
$`f(\rho)=U\,\mathrm{diag}(f(\lambda_i))\,U^*`$ の定義（1.3節）により $`\mathrm{Tr}\,f(\rho)=\sum_if(\lambda_i)`$ である。

<a id="r-3-3"></a>

**定理 3.3（セクターの混合のエントロピー）.** $`p_c`$ と $`\sigma_c`$ を上のとおりとする
（$`p_c=0`$ のセクターは省く）と、
```math
S(\rho)=H(p)+\sum_cp_c\,S(\sigma_c),\qquad H(p)=-\sum_cp_c\log p_c .
```

*証明.* $`\rho`$ の非零固有値（重複度込み）は、$`\sigma_c`$ の非零
固有値 $`\mu`$（$`p_c>0`$ のセクター）に対する数 $`p_c\mu`$ である。零固有値は両辺に $`0`$ を寄与する。したがって
$`-\sum_c\sum_{\mu>0}p_c\mu\log(p_c\mu)=-\sum_cp_c\log p_c\sum_{\mu>0}\mu-\sum_cp_c\sum_{\mu>0}\mu\log\mu`$
であり、$`\sum_{\mu>0}\mu=\mathrm{Tr}\,\sigma_c=1`$ である。$`\square`$

*計算：* ランダムな状態に対する定理 3.3、`DensityMatrix.von_neumann_entropy` を用いて — [チュートリアル、3.4節](../../tutorial/ja/03_traces_and_states.md#s-3-4)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. この公式は系についての不確かさを二つの部分に分ける：$`H(p)`$ は電荷（系がどのセクターにあるか）に
> ついての古典的な不確かさであり、$`\sum_cp_cS(\sigma_c)`$ はセクター内部の平均的な量子的不確かさである。超選択則の
> もとでは、前者は古典的な確率変数、すなわち電荷のエントロピーである。

他のトレースに関するエントロピーは $`S_w(\rho)=-\mathrm{tr}_w(\rho_w\log\rho_w)`$ である。ここで
同じ規約 $`\rho_w\log\rho_w:=f(\rho_w)`$ を用いる。
定理 3.2 により $`(\rho_w)_c=\rho_c/w_c`$ であるから、
```math
S_w(\rho)=-\sum_cw_c\mathrm{Tr}\Big(\frac{\rho_c}{w_c}\log\frac{\rho_c}{w_c}\Big)
=S(\rho)+\sum_cp_c\log w_c .
```
エントロピーは状態 *と* トレースの性質である：同じ状態でも
`Tr_blunt` と `tau_vN` に関しては異なるエントロピーをもつ。

*計算：* $`\tau_{\mathrm{vN}}`$ に関するエントロピーと $`S_w(\rho)`$ の公式 — [チュートリアル、3.4節](../../tutorial/ja/03_traces_and_states.md#s-3-4)。

<a id="s-3-5"></a>

## 3.5 相対エントロピーとクラインの不等式

密度 $`\rho,\sigma`$（$`\sigma>0`$）をもつ状態に対し、**相対エントロピー** を
$`S(\rho\|\sigma)=\mathrm{Tr}\,\rho\log\rho-\mathrm{Tr}\,\rho\log\sigma`$ と定める。ここで
$`\mathrm{Tr}\,\rho\log\rho:=\sum_i\lambda_i\log\lambda_i`$ は $`\rho`$ の固有値にわたる和であり、
$`0\log0=0`$ とする（各ブロックについても同様）。これはトレースに依存しない：定理 3.2 により
$`(\rho_w)_c=\rho_c/w_c`$、$`(\sigma_w)_c=\sigma_c/w_c`$ であるから、$`(\rho_w)_c`$ の固有値は $`\rho_c`$ の固有値を $`w_c`$ で割ったものであり、
$`\log(\sigma_c/w_c)=\log\sigma_c-(\log w_c)1`$ であり（$`\sigma_c=U\,\mathrm{diag}(s_i)\,U^*`$ ならば
$`\sigma_c/w_c=U\,\mathrm{diag}(s_i/w_c)\,U^*`$ であるから、定理 1.5(a) による）、
$`\mathrm{tr}_w(\rho_w\log\rho_w)`$ と $`\mathrm{tr}_w(\rho_w\log\sigma_w)`$ はそれぞれ
$`\mathrm{Tr}\,\rho\log\rho`$ と $`\mathrm{Tr}\,\rho\log\sigma`$ から同じ量 $`-\sum_cp_c\log w_c`$ だけずれる
（規約 $`0\log0=0`$ により、$`\rho`$ の対数はその台の上でのみ必要である）。したがって差は変わらない。

<a id="r-3-4"></a>

**補題 3.4（接線不等式）.** $`f:I\to\mathbb R`$ を区間 $`I`$ 上の凸関数とし、
$`b\in I`$ で微分可能（$`b`$ が $`I`$ の端点ならば片側微分可能）とする。このときすべての $`a\in I`$ に対し $`f(a)\ge f(b)+(a-b)f'(b)`$ であり、
$`f`$ が狭義凸、すなわち $`I`$ の $`x\neq y`$ と $`t\in(0,1)`$ に対し
$`f((1-t)x+ty)<(1-t)f(x)+tf(y)`$ ならば、$`a\neq b`$ に対して不等号は狭義である。さらに、$`f`$ が $`I`$ 上で連続であり、
$`I`$ の内部で2回微分可能で $`f''>0`$ ならば、$`f`$ は狭義凸である。

*証明.* $`t\in(0,1]`$ に対し、凸性から $`f(b+t(a-b))\le(1-t)f(b)+tf(a)`$、すなわち
$`\big(f(b+t(a-b))-f(b)\big)/t\le f(a)-f(b)`$ を得る。$`t\to0`$ とすれば不等式を得る。$`f`$
が狭義凸で $`a\neq b`$ のとき、$`m=(a+b)/2`$ とおくと
$`f(a)-f(b)>2\big(f(m)-f(b)\big)\ge2(m-b)f'(b)=(a-b)f'(b)`$ である。第一段は $`m`$ における狭義
凸性により、第二段は今証明した不等式を $`m`$ に適用したものである。最後の主張については、
$`I`$ の $`x<y`$、$`t\in(0,1)`$ をとり、$`z=(1-t)x+ty`$ とおく。平均値の定理により、
$`f(z)-f(x)=t(y-x)f'(\xi_1)`$ かつ $`f(y)-f(z)=(1-t)(y-x)f'(\xi_2)`$ を満たす
$`\xi_1\in(x,z)`$ と $`\xi_2\in(z,y)`$ が存在し、$`[\xi_1,\xi_2]`$ 上で $`f''>0`$ であるから、$`f'`$ に対する
平均値の定理により $`f'(\xi_1)<f'(\xi_2)`$ である。$`t(1-t)(y-x)>0`$ を掛けると
$`(1-t)\big(f(z)-f(x)\big)<t\big(f(y)-f(z)\big)`$、すなわち $`f(z)<(1-t)f(x)+tf(y)`$ を得る。$`\square`$

<a id="r-3-5"></a>

**定理 3.5（クラインの不等式）.** $`f:I\to\mathbb R`$ を区間 $`I`$ 上の凸関数とし、
$`X,Y`$ をスペクトルが $`I`$ に含まれる自己共役作用素とし、$`f`$ は $`Y`$ のすべての固有値で微分可能
（$`I`$ の端点では片側微分可能）とする。
このとき $`\mathrm{Tr}[f(X)-f(Y)-(X-Y)f'(Y)]\ge0`$ である。$`f`$ が狭義凸ならば、等号は
$`X=Y`$ のときに限り成り立つ。

*証明.* 正規直交固有基底により $`X=\sum_ia_i|u_i\rangle\langle u_i|`$、$`Y=\sum_jb_j|v_j\rangle\langle v_j|`$ と書き、
$`c_{ij}=|\langle u_i,v_j\rangle|^2`$ とおく。すると $`\sum_ic_{ij}=\sum_jc_{ij}=1`$ である。
このとき $`\mathrm{Tr}\,f(X)=\sum_{ij}c_{ij}f(a_i)`$、$`\mathrm{Tr}\,f(Y)=\sum_{ij}c_{ij}f(b_j)`$、
$`\mathrm{Tr}\,Xf'(Y)=\sum_{ij}c_{ij}a_if'(b_j)`$、$`\mathrm{Tr}\,Yf'(Y)=\sum_{ij}c_{ij}b_jf'(b_j)`$ である。
したがって左辺は $`\sum_{ij}c_{ij}[f(a_i)-f(b_j)-(a_i-b_j)f'(b_j)]`$ に等しく、これは補題 3.4 により
非負の項の和である。$`f`$ が狭義凸で和が零ならば、
$`c_{ij}>0`$ のときは常に $`a_i=b_j`$ であるから、すべての $`j`$ に対し $`Xv_j=\sum_ia_i\langle u_i,v_j\rangle u_i=b_jv_j`$ となり、
$`X=Y`$ である。$`\square`$

<a id="r-3-6"></a>

**系 3.6.** $`S(\rho\|\sigma)\ge0`$ であり、等号は $`\rho=\sigma`$ のときかつそのときに限り成り立つ。特に
$`S(\rho)\le\log D`$ であり、等号は $`\rho=1/D`$ のときかつそのときに限り成り立つ。

*証明.* 定理 3.5 を $`f(t)=t\log t`$, $`f(0)=0`$ に適用する。これは
$`I=[0,\infty)`$ 上で補題 3.4 により狭義凸であり（$`t\to0^+`$ のとき $`t\log t\to0`$ であるから $`f`$ は $`0`$ で連続であり、
$`(0,\infty)`$ 上で $`f''(t)=1/t>0`$）、
$`(0,\infty)\supseteq\mathrm{spec}\,\sigma`$ 上で微分可能で $`f'(t)=\log t+1`$ である。
$`X=\rho`$, $`Y=\sigma`$ とすると、左辺は
$`\mathrm{Tr}\,\rho\log\rho-\mathrm{Tr}\,\sigma\log\sigma-\mathrm{Tr}(\rho-\sigma)(\log\sigma+1) =S(\rho\|\sigma)`$ である。$`\mathrm{Tr}\,\rho=\mathrm{Tr}\,\sigma=1`$ だからである。第二の主張については
$`\sigma=1/D`$ ととる：$`0\le S(\rho\|1/D)=-S(\rho)+\log D`$。$`\square`$

*計算：* $`S(\rho)\le\log D`$ と最大混合状態 — [チュートリアル、3.5節](../../tutorial/ja/03_traces_and_states.md#s-3-5)。

<a id="r-3-7"></a>

**定理 3.7（連鎖律；中心への条件付き期待値に関する単調性）.** $`\rho,\sigma`$ を
$`\sigma>0`$ なる状態とし、$`p_c=\mathrm{Tr}\,\rho_c`$、$`q_c=\mathrm{Tr}\,\sigma_c>0`$、$`\hat\sigma_c=\sigma_c/q_c`$ とおき、
$`p_c>0`$ に対し $`\hat\rho_c=\rho_c/p_c`$ とおく。このとき、和を $`p_c>0`$ なる $`c`$ にわたってとると（その他の
セクターは $`0\log0=0`$ により $`0`$ を寄与する）、
```math
S(\rho\|\sigma)=\sum_cp_c\log\frac{p_c}{q_c}+\sum_cp_c\,S(\hat\rho_c\|\hat\sigma_c)
\ \ge\ \sum_cp_c\log\frac{p_c}{q_c}=S(E\rho\|E\sigma),
```
が成り立つ。ここで $`E`$ は中心の上への条件付き期待値である（命題 2.8）。

*証明.* $`p_c>0`$ に対し $`\mathrm{Tr}\,\rho_c\log\rho_c=p_c\log p_c+p_c\mathrm{Tr}\,\hat\rho_c\log\hat\rho_c`$
（固有値は $`p_c\mu`$、定理 3.3 と同様）かつ $`\log\sigma_c=(\log q_c)1+\log\hat\sigma_c`$ であるから（$`\hat\sigma_c=U\,\mathrm{diag}(s_i)\,U^*`$ と対角化すれば
$`\sigma_c=U\,\mathrm{diag}(q_cs_i)\,U^*`$ であり、定理 1.5(a) を用いる）、
$`\mathrm{Tr}\,\rho_c\log\sigma_c=p_c\log q_c+p_c\mathrm{Tr}\,\hat\rho_c\log\hat\sigma_c`$ であり、
$`\mathrm{Tr}\,\rho_c\log\rho_c-\mathrm{Tr}\,\rho_c\log\sigma_c=p_c\log(p_c/q_c)+p_cS(\hat\rho_c\|\hat\sigma_c)`$ である。
$`p_c=0`$ ならば $`\rho_c=0`$ であり、そのセクターは $`0`$ を寄与する。
$`c`$ について和をとる。第二の和は、系 3.6 を $`M_{k_c}`$ において状態 $`\hat\rho_c`$ と $`\hat\sigma_c>0`$ に適用すれば $`\ge0`$ である。最後に
$`E\rho=\bigoplus_c(p_c/k_c)1_{k_c}`$ と $`E\sigma=\bigoplus_c(q_c/k_c)1_{k_c}>0`$ は状態である。$`E\rho`$ は状態 $`\omega\circ E`$
（ただし $`\omega(A)=\mathrm{Tr}\,\rho A`$）の密度である。実際
$`\mathrm{Tr}(\rho E(A))=\sum_cp_c\mathrm{Tr}(A_c)/k_c=\mathrm{Tr}(E(\rho)A)`$ である。そして
$`S(E\rho\|E\sigma)=\sum_ck_c\frac{p_c}{k_c}\big(\log\frac{p_c}{k_c}-\log\frac{q_c}{k_c}\big)=\sum_cp_c\log\frac{p_c}{q_c}`$ である。ここで
ブロック $`c`$ における $`E\rho`$ の固有値は $`p_c/k_c`$（$`k_c`$ 重）であり、$`p_c=0`$ のセクターは $`0`$ を寄与する。$`\square`$

*計算：* ランダムな状態に対する定理 3.7、`states.relative_entropy` と `channels.center_expectation` を用いて — [チュートリアル、3.5節](../../tutorial/ja/03_traces_and_states.md#s-3-5)。

したがって、二つの状態をともに $`E`$ による像で置き換えると、相対エントロピーは減少することしかできない。残るのはセクター分布のカルバック・ライブラー情報量 $`\sum_cp_c\log(p_c/q_c)`$ である。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. 相対エントロピーは
> 二つの状態をどれだけよく識別できるかを測る。セクター内部の量子情報を忘れることは、二つの状態を
> 識別しにくくすることしかできない。残るのは電荷分布の古典的な識別可能性である。

<a id="s-3-6"></a>

## 3.6 マイケルソン・コントラスト

$`A\ge0`$, $`A\neq0`$ に対し、**マイケルソン・コントラスト** を
$`\Delta(A)=(\lambda_{\max}-\lambda_{\min})/(\lambda_{\max}+\lambda_{\min})\in[0,1]`$ と定める。ここで
両端の固有値はすべてのセクターにわたってとる（また $`\Delta(0)=0`$ とする）。これは $`A`$ が
単位元のスカラー倍からどれだけ離れているかを測る：$`\Delta(A)=0`$ $`\iff`$ $`A=\lambda1`$。

<a id="r-3-8"></a>

**命題 3.8.** $`\Phi:\mathcal M\to\mathcal M`$ が線形、正かつ単位的
（$`\Phi(1)=1`$）ならば、すべての $`A\ge0`$ に対し $`\Delta(\Phi(A))\le\Delta(A)`$ である。

*証明.* $`\lambda_{\min}1\le A\le\lambda_{\max}1`$ であり（$`A-\lambda_{\min}1`$ と $`\lambda_{\max}1-A`$ に対する系 1.4(3)）、
正値性、線形性、単位性から
$`\lambda_{\min}1\le\Phi(A)\le\lambda_{\max}1`$ を得る。したがって $`\Phi(A)`$ は自己共役であり（正の行列は自己共役である、系 1.4(3)）、
$`\Phi(A)\ge0`$ であり（$`\lambda_{\min}\ge0`$ だから）、系 1.4(3) により $`\lambda_{\min}(\Phi(A))\ge\lambda_{\min}`$ かつ
$`\lambda_{\max}(\Phi(A))\le\lambda_{\max}`$ である。$`\Phi(A)=0`$ ならば主張は自明である。そうでなければ $`A\neq0`$ であり、
$`\lambda_{\max}`$ と $`\lambda_{\max}(\Phi(A))`$ はともに $`>0`$ であり、$`r=\lambda_{\min}/\lambda_{\max}\in[0,1]`$ とおくと
$`\Delta=(1-r)/(1+r)`$ は $`r`$ について減少であり、$`r(\Phi(A))\ge r(A)`$ である。$`\square`$

*計算：* ランダムな混合ユニタリ・チャネルと $`E`$ に対する命題 3.8、`Operator.michelson_contrast` を用いて — [チュートリアル、3.6節](../../tutorial/ja/03_traces_and_states.md#s-3-6)。

<a id="s-3-7"></a>

## 3.7 もっともらしいが偽である不等式

$`M_k`$ の正の $`X,Y`$ とユニタリ $`U`$ に対し、$`|\mathrm{Tr}(XUY)|\le\mathrm{Tr}(XY)`$ は成り立つだろうか。これは
$`U=1`$ のときは成り立ち、コーシー・シュワルツ型の主張のように見える。しかし $`X`$
と $`Y`$ が可換でないときは常に偽である。

<a id="r-3-9"></a>

**命題 3.9.** $`X,Y\ge0`$ を $`M_k`$ の元とする。このとき $`\max_U|\mathrm{Tr}(XUY)|=\|YX\|_1\ge\mathrm{Tr}(XY)`$ である（最大値は
ユニタリ $`U\in M_k`$ 全体にわたってとる）。等号は $`XY=YX`$ のときかつそのときに限り成り立つ。

*証明.* $`\mathrm{Tr}(XUY)=\mathrm{Tr}(YX\,U)`$ であるから、定理 1.9 により最大値は $`\|YX\|_1`$ であり、
$`\|YX\|_1\ge|\mathrm{Tr}(YX)|=\mathrm{Tr}(XY)`$ である。ここで
$`X^{1/2}YX^{1/2}=(Y^{1/2}X^{1/2})^*(Y^{1/2}X^{1/2})\ge0`$（系 1.7）より
$`\mathrm{Tr}(XY)=\mathrm{Tr}(X^{1/2}YX^{1/2})\ge0`$ であることを用いた。$`XY=YX`$ ならば
$`X^{1/2}`$ は $`Y`$ と可換であるから（定理 1.5(d)）、$`YX=X^{1/2}YX^{1/2}\ge0`$ であり、
$`\|YX\|_1=\mathrm{Tr}(YX)`$ である（$`P\ge0`$ に対し $`\|P\|_1=\mathrm{Tr}\,P`$、1.4節）。逆に
$`\|YX\|_1=\mathrm{Tr}(YX)`$ と仮定する。$`|YX|=0`$ ならば $`YX=0=(YX)^*=XY`$ である。そうでなければ
$`YX=W|YX|`$（定理 1.8）と書き、$`R=|YX|^{1/2}\neq0`$ とおく。巡回性により
$`\mathrm{Tr}(YX)=\mathrm{Tr}(RWR)=\langle R,WR\rangle_{HS}`$ であり、一方
$`\|R\|_2^2=\|WR\|_2^2=\mathrm{Tr}|YX|=\|YX\|_1`$ である。よって $`\langle R,WR\rangle_{HS}=\|R\|_2\|WR\|_2`$ であり、
コーシー・シュワルツの不等式で等号が成り立つから、あるスカラー $`\lambda`$ に対し $`WR=\lambda R`$ であり、
$`\lambda\|R\|_2^2=\langle R,WR\rangle_{HS}=\|R\|_2^2`$ から $`\lambda=1`$ を得る。したがって $`W|YX|^{1/2}=|YX|^{1/2}`$ かつ $`YX=W|YX|=|YX|\ge0`$ である。正
行列は自己共役であるから、$`YX=(YX)^*=XY`$ である。$`\square`$

定理 1.9 により、$`YX=W|YX|`$ に対しユニタリ $`U=W^*`$ が最大値 $`\|YX\|_1`$ を達成する。ランダムなユニタリは平均的にはより悪い：ハール分布に従う $`U`$
（定理 5.3）に対し、命題 11.11 を $`A=YX`$ として適用すると
$`\big(\mathbb E|\mathrm{Tr}(XUY)|^2\big)^{1/2}=\|YX\|_2/\sqrt k\le\|YX\|_1/\sqrt k=\max_V|\mathrm{Tr}(XVY)|/\sqrt k`$
を得る（$`\sum_is_i^2\le(\sum_is_i)^2`$ だから）。また系 11.12 は、スペクトルが $`[a,b]`$（$`a>0`$）に含まれる
$`X,Y`$ に対し、一つのサンプルが $`\mathrm{Tr}(XY)`$ に達する確率は高々 $`b^4/(a^4k^2)`$ であることを示す
（これが意味をもつのは $`k>b^2/a^2`$ のときだけである）。

*計算：* $`M_6`$ の可換でない組に対する命題 3.9、最大値を与える $`W^*`$、およびハール・ユニタリ上のランダム探索 — [チュートリアル、3.7節](../../tutorial/ja/03_traces_and_states.md#s-3-7)。

<a id="s-3-8"></a>

## 3.8 比較の対象

*ライブラリでは：* 他のソフトウェアとの比較 — [チュートリアル、3.8節](../../tutorial/ja/03_traces_and_states.md#s-3-8)。

<a id="x-3"></a>

## 演習問題

1. トレース状態 $`\tau_{\mathrm{vN}}`$（セクターの重みが等しい）の
   $`\tau_{\mathrm{vN}}`$ に関する密度が $`1`$ であることを示し、その blunt 密度とエントロピー $`S`$ を計算せよ。
2. $`\mathrm{tr}_w`$ を忠実とする。$`S_w(\rho)\le\log\mathrm{tr}_w(1)`$ であり、等号が
   $`\rho_w=1/\mathrm{tr}_w(1)`$ のときかつそのときに限り成り立つことを示せ。
3. $`M_D`$ 上の（$`\mathcal M`$ 上ではない！）測定チャネル $`\rho\mapsto\sum_cP_c\rho P_c`$ は状態を乱すことを示し、二つのセクターの重ね合わせについて測定前後の $`S(\rho)`$ を計算せよ。

<details>
<summary>解答</summary>

1. $`\tau_{\mathrm{vN}}(A)=\tau_{\mathrm{vN}}(1\cdot A)`$ であるから、一意性により $`\rho_{\tau}=1`$ である。
   定理 3.2 により blunt 密度は $`\bigoplus_c1_{k_c}/(Ck_c)`$ であり、$`p_c=1/C`$、
   $`\sigma_c=1/k_c`$ であるから、定理 3.3 により $`S=\log C+\frac1C\sum_c\log k_c`$ である。

2. 和は $`p_c>0`$ なる $`c`$ にわたる。このとき
   $`S_w(\rho)=S(\rho)+\sum_cp_c\log w_c=-\sum_cp_c\log p_c+\sum_cp_cS(\sigma_c)+\sum_cp_c\log w_c    \le\sum_cp_c\log\frac{k_cw_c}{p_c}\le\log\sum_{p_c>0}k_cw_c\le\log\mathrm{tr}_w(1)`$ である。ここで
   $`S(\sigma_c)\le\log k_c`$（$`M_{k_c}`$ における系 3.6）と、有限形のイェンセンの不等式
   $`\sum_{c=1}^mp_c\log x_c\le\log\sum_{c=1}^mp_cx_c`$（$`p_c>0`$, $`\sum_cp_c=1`$, $`x_c>0`$）を
   $`x_c=k_cw_c/p_c`$ に適用して用いた。この不等式で等号が成り立つのはすべての $`x_c`$ が等しいときに限る。$`m`$ に関する帰納法で証明する：
   $`m=1`$ では示すことはない。$`m\ge2`$ のとき $`q=1-p_m\in(0,1)`$、
   $`\bar x=\sum_{c<m}(p_c/q)x_c`$ とおく。帰納法の仮定により
   $`\sum_{c<m}(p_c/q)\log x_c\le\log\bar x`$ であり、等号は $`x_1=\dots=x_{m-1}`$ のときに限る。また $`\log`$ の
   狭義凹性（2階導関数が $`1/t^2>0`$ である $`-\log`$ に対する補題 3.4）により
   $`q\log\bar x+p_m\log x_m\le\log(q\bar x+p_mx_m)`$
   であり、等号は $`\bar x=x_m`$ のときに限る。第一の不等式に $`q`$ を掛けて加えると
   $`\sum_cp_c\log x_c\le\log\sum_cp_cx_c`$ を得て、等号は $`x_1=\dots=x_{m-1}=\bar x=x_m`$ を強いる。
   すべてで等号が成り立つためには、すべての $`p_c>0`$（すべての $`k_cw_c>0`$ だから）、$`\sigma_c=1/k_c`$ かつ $`k_cw_c/p_c`$ が一定、すなわち
   $`p_c=k_cw_c/\mathrm{tr}_w(1)`$、すなわち $`(\rho_w)_c=\rho_c/w_c=1_{k_c}/\mathrm{tr}_w(1)`$ でなければならない。
   逆に $`\rho_w=1/\mathrm{tr}_w(1)`$ に対しては
   $`S_w(\rho)=-\mathrm{tr}_w(1)\cdot\frac1{\mathrm{tr}_w(1)}\log\frac1{\mathrm{tr}_w(1)}=\log\mathrm{tr}_w(1)`$ である。

3. $`M_D`$ において、二つの
   セクターの単位ベクトルの純粋な重ね合わせ $`\psi=(\psi_1\oplus\psi_2)/\sqrt2`$ は $`S=0`$ をもつ。測定後の状態は
   $`\frac12|\psi_1\rangle\langle\psi_1|\oplus\frac12|\psi_2\rangle\langle\psi_2|`$ であり $`S=\log2`$ である。しかし
   $`\mathcal M`$ に制限すると、$`\psi`$ のベクトル状態は混合と一致する：
   すべての $`A\in\mathcal M`$ に対し
   $`\langle\psi,A\psi\rangle=\tfrac12\langle\psi_1,A_1\psi_1\rangle+\tfrac12\langle\psi_2,A_2\psi_2\rangle`$
   である。$`\mathcal M`$ のいかなる元もコヒーレンスを検出しない（物理的には、これが超選択則の内容である）。

</details>

*計算：* 演習問題 1 の数値的検算 — [チュートリアル、演習問題](../../tutorial/ja/03_traces_and_states.md#x-3)。

**次へ：**[第4章](04_symmetries.md)では対称性と表現を扱う：物理系への群作用がどのようにしてセクターを生み出すか、そしてライブラリがスピン系や量子ビットレジスタに対応する環をどのように構成するか。
