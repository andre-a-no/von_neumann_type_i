[English](../13_jordan_algebras.md) · [Español](../es/13_jordan_algebras.md) · [Français](../fr/13_jordan_algebras.md) · [Deutsch](../de/13_jordan_algebras.md) · [中文](../zh/13_jordan_algebras.md) · **日本語** · [Русский](../ru/13_jordan_algebras.md)

# 13. ジョルダン代数：量子力学はなぜ *-環を用いるのか

二つのオブザーバブルの積 $`xy`$ は、$`x`$ と $`y`$ が可換でない限りオブザーバブルではないが、対称化された積
$`x\circ y=\frac12(xy+yx)`$ は常にオブザーバブルである。1933年にジョルダン（P. Jordan）は、この積を備えたオブザーバブルを量子力学の基本構造とすることを提案し、ジョルダン、フォン・ノイマン、ウィグナーはそこから生じる有限次元の代数を分類した。本章では有限次元の理論を完全な証明とともに展開する。証明は本章で与えるが、四つの長いもの（分類定理、ケッヒャー–ヴィンベルクの定理、カディソンの定理、および $`\bigoplus_cM_{k_c}`$ のジョルダン自己同型）は付録 E と F で与える。そして、ジョルダン構造が何を見落としているかを問う。答えは交換子の符号である：一つの系のジョルダン代数は交換子を符号を除いてしか決定しない（13.5節）が、複合系はその符号を必要とする：その符号がどの結合作用素が正であるかを決め、それから作られる PPT 判定がエンタングルメントを検出する（命題 13.12）。

**本章で学ぶこと**：ジョルダン積と、$`\circ`$ とは異なり正値性を保つ二次表現 $`P(y)x=yxy`$、スピン因子とそのスペクトル理論およびローレンツ錐、$`*`$-代数の自己共役部分ではないエルミート行列のジョルダン部分代数、Jordan–von Neumann–Wigner の分類、転置がなぜ $`*`$-自己同型ではなくジョルダン自己同型であるのか、そしてそれが複合系の上で何をするか、さらに第11章のトレース不等式をジョルダンの言葉で書き、スピン因子の中で証明すること。

本章の数値的検証は[チュートリアル 13](../../tutorial/ja/13_jordan_algebras.md)にある。

本章を通じて、$`\mathrm{Herm}_n(\mathbb K)`$ は $`\mathbb K=\mathbb R`$ または $`\mathbb C`$ 上の $`n\times n`$ エルミート行列のなす実ベクトル空間（13.4節では E.6節で構成される四元数 $`\mathbb H`$ と八元数 $`\mathbb O`$、およびそれらの上の $`\mathrm{Herm}_n`$（E.7節）も用いる）、$`M_{sa}`$ は $`\mathcal M=\bigoplus_cM_{k_c}(\mathbb C)`$ の自己共役部分、
$`X,Y,Z`$ はパウリ行列で $`\sigma=(X,Y,Z)`$ とし、したがって $`u\in\mathbb R^3`$ に対し $`u\cdot\sigma=u_1X+u_2Y+u_3Z`$ である。
$`[a,b]=ab-ba`$ とし、$`|u|`$ は $`u\in\mathbb R^m`$ のユークリッドノルムである。

<a id="s-13-1"></a>

## 13.1 オブザーバブルとジョルダン積

**定義.** **ジョルダン代数** とは、双線形かつ可換な積 $`x\circ y`$ をもつ実ベクトル空間 $`J`$ であって、
**ジョルダン恒等式** $`(x^2\circ y)\circ x=x^2\circ(y\circ x)`$（ここで $`x^2=x\circ x`$）を満たすもののことである。ある結合代数の部分空間で $`\frac12(xy+yx)`$ について閉じたものと同型であるとき、
**特殊** であるという。$`J`$ 上で
$`L(y)x=y\circ x`$ とおき、**二次表現** を $`P(y)=2L(y)^2-L(y^2)`$ で定める。

<a id="r-13-1"></a>

**補題 13.1.** $`\mathcal A`$ を実または複素の結合代数とし、$`x\circ y=\frac12(xy+yx)`$ とする。このとき
$`\circ`$ はジョルダン恒等式を満たし、ジョルダン冪 $`x^{\circ1}=x`$、$`x^{\circ(k+1)}=x\circ x^{\circ k}`$ は結合的な冪 $`x^k`$ に等しい。特に $`M_{sa}`$、$`\mathrm{Herm}_n(\mathbb C)`$ および
$`\mathrm{Herm}_n(\mathbb R)`$ は特殊ジョルダン代数である。

*証明.* 展開すると $`4(x^2\circ y)\circ x=x^2yx+yx^3+x^3y+xyx^2`$ および
$`4x^2\circ(y\circ x)=x^2yx+x^3y+yx^3+xyx^2`$ であり、同じ四つの項からなる。帰納法により
$`x\circ x^k=\frac12(x^{k+1}+x^{k+1})=x^{k+1}`$ である。$`x,y`$ がエルミートならば $`(xy+yx)^*=yx+xy`$ であるから、三つの空間は $`\circ`$ について閉じている。実対称な $`x,y`$ に対しては積も実である。$`\square`$

*計算：* $`M_2\oplus M_3`$ におけるジョルダン恒等式と冪結合性 — [チュートリアル、13.1節](../../tutorial/ja/13_jordan_algebras.md#s-13-1)。

<a id="r-13-2"></a>

**命題 13.2.** $`x,y,z\in M_{sa}`$ とする。

1. $`P(y)x=yxy`$ であり、**基本公式** $`P(P(y)x)=P(y)P(x)P(y)`$ が成り立つ。
2. $`\mathrm{Tr}\,P(y)x=\mathrm{Tr}\,y^2x`$。
3. $`x\ge0`$ ならば $`P(y)x\ge0`$ である。$`y`$ が可逆ならば、$`P(y)`$ は正錐からそれ自身への全単射であり、その
   逆写像は $`P(y^{-1})`$ である。
4. 正の元のジョルダン積は正とは限らない：$`p=\begin{pmatrix}1&0\\0&0\end{pmatrix}`$ と
   $`q=\frac12\begin{pmatrix}1&1\\1&1\end{pmatrix}`$ に対し $`\det(p\circ q)=-\frac1{16}`$ である。

*証明.* (1) $`2y\circ(y\circ x)=\frac12(y^2x+2yxy+xy^2)`$ および $`y^2\circ x=\frac12(y^2x+xy^2)`$ であるから、差をとればよい。次に結合律により
$`P(P(y)x)z=yxy\,z\,yxy=P(y)P(x)P(y)z`$ である。(2) トレースの巡回性。(3)
$`yxy=(x^{1/2}y)^*(x^{1/2}y)\ge0`$（系 1.7）であり、$`P(y)P(y^{-1})=P(y^{-1})P(y)=\mathrm{id}`$ である。(4)
$`p\circ q=\frac14\begin{pmatrix}2&1\\1&0\end{pmatrix}`$。$`\square`$

*計算：* $`M_2\oplus M_3`$ における命題 13.2、および例 (4) — [チュートリアル、13.1節](../../tutorial/ja/13_jordan_algebras.md#s-13-1)。*ライブラリでは：* `jordan.jordan_product`、`jordan.jordan_power`、`jordan.quadratic_representation`。`TypeIAlgebra` のセクターにわたってバッチ化されている — [チュートリアル、13.1節](../../tutorial/ja/13_jordan_algebras.md#s-13-1)。

$`P(y)`$ を $`2L(y)^2-L(y^2)`$ と書けば (1) の恒等式はすべて $`\circ`$ だけで表される。(1) の証明は行列の積の結合律しか用いないので、すべての特殊ジョルダン代数 $`J\subseteq\mathcal A`$（$`\mathcal A`$ は結合的）において、$`J`$ 上で $`P(y)x=yxy`$ と基本公式 $`P(P(y)x)=P(y)P(x)P(y)`$ を与える。（これは本講義でそれを用いるすべての代数を含む。特殊でない代数に対しては基本公式は用いない。）

$`\mathrm{Herm}_n(\mathbb C)`$（$`n\ge2`$）上の積 $`\alpha xy+\beta yx`$（$`\alpha,\beta\in\mathbb R`$）のうち、$`x,y`$ について対称で $`y=x`$ のとき $`x^2`$ に等しいものは $`\frac12(xy+yx)`$ だけである：対称性はすべての $`x,y`$ に対し $`(\alpha-\beta)(xy-yx)=0`$ を意味し、$`XZ\ne ZX`$（左上の $`2\times2`$ ブロックにおいて）により $`\alpha=\beta`$ となる。さらに $`x=y=1`$ に対する $`(\alpha+\beta)1=1`$ から $`\alpha=\beta=\frac12`$ を得る。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. $`x\circ y`$ はオブザーバブルの対称化された積であり（上で示したように、そのような組合せは一意である）、
> 再びオブザーバブルとなる。二次表現は **リューダース則** である：イエス／ノーの問い $`p`$ を測定し、結果「イエス」を残すと、
> 密度 $`\rho`$ は $`P(p)\rho/\mathrm{Tr}\,P(p)\rho=p\rho p/\mathrm{Tr}(p\rho)`$ に変わる。命題 13.2(4)
> は、$`\circ`$ そのものには正の元に対してこのような操作的意味がないことを述べている。

<a id="s-13-2"></a>

## 13.2 スピン因子とローレンツ錐

**定義.** $`m\ge1`$ に対し、**スピン因子** $`V_m=\mathbb R\oplus\mathbb R^m`$ は元 $`x=(a,u)`$、
積 $`(a,u)\circ(b,v)=(ab+u\cdot v,\;av+bu)`$ および単位元 $`1=(1,0)`$ をもつ。その **トレース** は $`\mathrm{tr}\,x=2a`$ である。

積は双線形かつ可換である。ジョルダン恒等式は行列表現から従う。

<a id="r-13-3"></a>

**補題 13.3（反可換な対称元）.** $`p=\lfloor m/2\rfloor`$、$`D=2^p`$（$`m=1`$ のときは $`D=2`$）とする。
$`(\mathbb C^2)^{\otimes p}`$ 上で、$`j\le p`$ に対し $`s_{2j-1}=Z^{\otimes(j-1)}\otimes X\otimes1`$、$`s_{2j}=Z^{\otimes(j-1)}\otimes Y\otimes1`$
とおき、$`m`$ が奇数ならば $`s_m=Z^{\otimes p}`$ とおく（$`m=1`$ のときは $`s_1=Z`$）。このとき $`s_i=s_i^*`$ かつ
$`s_is_j+s_js_i=2\delta_{ij}1`$ である。写像 $`\phi(a,u)=a1+\sum_iu_is_i`$ は単射であり、
$`\phi(x\circ y)=\phi(x)\circ\phi(y)`$ を満たす。したがって $`V_m`$ は特殊ジョルダン代数である。$`m\ge2`$ のとき、$`s_i`$ たちは $`*`$-代数として
$`M_D(\mathbb C)`$ を生成する。

*証明.* パウリ行列は二乗すると $`1`$ になり、互いに反可換である。$`Z`$ と $`X`$、$`Z`$ と $`Y`$ も同様である。二つの異なる文字列
$`s_i,s_j`$ は、ちょうど一つのサイトで反可換な仕方で異なる：$`s_{2j-1},s_{2j}`$ についてはサイト $`j`$ で（$`X`$ 対
$`Y`$）、添字 $`j<j'`$ をもつ文字列についてはサイト $`j`$ で（$`X`$ または $`Y`$ 対 $`Z`$）異なり、その他のサイトでは等しい因子か因子 $`1`$ をもつ。$`s_m=Z^{\otimes p}`$ についても同じことが成り立つ。よって $`i\ne j`$ に対し $`s_is_j=-s_js_i`$ であり、
$`s_i^2=1`$ である。ゆえに $`(u\cdot s)\circ(v\cdot s)=\sum_{ij}u_iv_j\,s_i\circ s_j=(u\cdot v)1`$ であり、これから
$`\phi(x)\circ\phi(y)=\phi(x\circ y)`$ が得られる。各 $`s_i`$ はトレースが零であり（トレース零の因子をもつテンソル積である）、
$`\mathrm{Tr}(s_is_j)=\mathrm{Tr}(s_i\circ s_j)=D\delta_{ij}`$ であるから、$`1,s_1,\dots,s_m`$ は
$`\mathrm{Tr}(ab)`$ について直交し、$`\phi`$ は単射である。ジョルダン恒等式は $`\mathrm{Herm}_D`$ で成り立ち（補題 13.1）、
$`\phi`$ によって引き戻される。最後に $`s_{2j-1}s_{2j}=1\otimes\cdots\otimes XY\otimes\cdots\otimes1=iZ_j`$（$`Z`$ はサイト
$`j`$ にある）であるから、生成される $`*`$-代数は $`Z_j`$ を含み、したがって $`X_j=Z_1\cdots Z_{j-1}s_{2j-1}`$ と $`Y_j`$ を含み、それゆえパウリ行列のすべての
テンソル積、すなわち $`M_D(\mathbb C)`$ の基底を含む。$`\square`$

*ライブラリでは：* $`s_i`$ は `jordan.anticommuting_symmetries(m)` であり、写像 $`\phi`$ を伴う $`V_m`$ はクラス `SpinFactor`（`to_matrix`）である — [チュートリアル、13.2節](../../tutorial/ja/13_jordan_algebras.md#s-13-2)。

符号を除けば、$`s_{2j-1}`$ と $`s_{2j}`$ は 9.4節のジョルダン–ウィグナー・フェルミオン $`c_j`$ の
マヨラナ演算子 $`c_j+c_j^*`$ と $`i(c_j-c_j^*)`$ である。ここで第9章のサイト $`j-1`$ を第 $`j`$ テンソル因子とする：そこでは $`s^a=\frac12\sigma_a`$、
$`n=s^z+\frac12`$、$`s^\pm=s^x\pm is^y`$ であるから、$`(-1)^n=1-2n=-Z`$、$`s^-=\frac12(X-iY)`$ であり、したがって
$`c_j=(-1)^{j-1}Z^{\otimes(j-1)}\otimes s^-\otimes1`$ である。$`s^-+s^+=X`$ かつ $`i(s^--s^+)=i(-iY)=Y`$ であるから
```math
s_{2j-1}=(-1)^{j-1}(c_j+c_j^*),\qquad s_{2j}=(-1)^{j-1}\,i(c_j-c_j^*).
```
である。

*計算：* $`p=3`$ に対するこれらの恒等式 — [チュートリアル、13.2節](../../tutorial/ja/13_jordan_algebras.md#s-13-2)。

<a id="r-13-4"></a>

**定理 13.4（$`V_m`$ のスペクトル理論）.** $`x=(a,u)\in V_m`$、$`\lambda_\pm=a\pm|u|`$ とする。

1. $`V_m`$ の冪等元は $`0`$、$`1`$ および $`|w|=1`$ なる $`c_w=\frac12(1,w)`$ である。また $`c_w\circ c_{w'}=0`$ となるのは $`w'=-w`$ のときかつそのときに限る。
2. $`u\ne0`$ ならば $`e_\pm=c_{\pm u/|u|}`$ として $`x=\lambda_+e_++\lambda_-e_-`$ であり、$`u=0`$ ならば $`x=a1`$ である。$`x=\mu_1c_w+\mu_2c_{-w}`$ ならば
   $`\{\mu_1,\mu_2\}=\{\lambda_+,\lambda_-\}`$ である。$`\lambda_\pm`$ を $`x`$ の **固有値** とよぶ。
3. $`f:\{\lambda_+,\lambda_-\}\to\mathbb R`$ に対し $`f(x)=f(\lambda_+)e_++f(\lambda_-)e_-`$（$`u=0`$ のときは $`f(a)1`$）とおく。このとき
   $`x^{\circ k}=\lambda_+^ke_++\lambda_-^ke_-`$、$`\phi(f(x))=f(\phi(x))`$ であり、$`\phi(x)`$ の固有値は $`\lambda_\pm`$ で、
   それぞれの重複度は $`D/2`$、そして $`\mathrm{Tr}\,\phi(x)=\frac D2\mathrm{tr}\,x`$ である。また
   $`\mathrm{tr}\,x=\lambda_++\lambda_-`$、$`\det x:=\lambda_+\lambda_-=a^2-|u|^2`$ である。
4. 平方の錐 $`Q=\{x^2:x\in V_m\}`$ は **ローレンツ錐** $`\{(a,u):a\ge|u|\}=\{x:\lambda_-\ge0\}=\{x:\phi(x)\ge0\}`$ である。
   これは **自己双対** である：$`y\in Q`$ であるための必要十分条件は、すべての $`x\in Q`$ に対し $`\mathrm{tr}(x\circ y)\ge0`$ となることである。
5. $`\phi(P(y)x)=\phi(y)\phi(x)\phi(y)`$。したがってすべての $`y\in V_m`$ に対し $`P(y)Q\subseteq Q`$ である。

*証明.* (1) $`(a,u)^2=(a^2+|u|^2,2au)`$ が $`(a,u)`$ に等しいのは、$`u=0`$ かつ $`a\in\{0,1\}`$ のとき、または $`a=\frac12`$ かつ $`|u|=\frac12`$ のときに限る。
$`c_w\circ c_{w'}=\frac14(1+w\cdot w',\,w+w')`$ が零になるのは $`w'=-w`$ のときかつそのときに限る。(2)
$`\lambda_+e_++\lambda_-e_-=(a,0)+|u|(0,u/|u|)=x`$ であり、$`\mu_1c_w+\mu_2c_{-w}=(\frac{\mu_1+\mu_2}2,\frac{\mu_1-\mu_2}2w)`$ であるから、
$`\mu_1+\mu_2=2a`$ かつ $`|\mu_1-\mu_2|=2|u|`$ である。(3) (1) により $`(\alpha e_++\beta e_-)\circ(\gamma e_++\delta e_-)=\alpha\gamma e_++\beta\delta e_-`$ であり、帰納法で冪の公式を得る。$`w=u/|u|`$ とすると $`(w\cdot s)^2=1`$ であるから、
$`\phi(e_\pm)=\frac12(1\pm w\cdot s)`$ は和が $`1`$、積が $`0`$ の直交射影で、それぞれのトレースは $`D/2`$ である。したがって
$`\phi(x)=\lambda_+\phi(e_+)+\lambda_-\phi(e_-)`$ はスペクトル分解であり、
$`f(\phi(x))=f(\lambda_+)\phi(e_+)+f(\lambda_-)\phi(e_-)=\phi(f(x))`$ である（定理 1.5）。(4) $`x^2`$ の固有値は
$`\lambda_\pm^2\ge0`$ である。逆に $`\lambda_-\ge0`$ ならば $`x=(\sqrt x)^2`$ である。$`\lambda_-\ge0`$ は $`a\ge|u|`$ と同値であり、(3) により
$`\phi(x)\ge0`$ とも同値である。自己双対性については、$`\mathrm{tr}(x\circ y)=2(ab+u\cdot v)`$ である。$`a\ge|u|`$ かつ $`b\ge|v|`$ ならば
$`ab\ge|u||v|\ge-u\cdot v`$ である。$`b<|v|`$ ならば、$`x=(|v|,-v)\in Q`$ が $`2(|v|b-|v|^2)<0`$ を与える（$`v=0`$、$`b<0`$ のときは $`x=1`$ をとる）。
(5) $`P(y)x=2y\circ(y\circ x)-y^2\circ x`$ であり $`\phi`$ は $`\circ`$ を保つので、命題 13.2(1) により $`\phi(P(y)x)=\phi(y)\phi(x)\phi(y)`$ である。
あとは (4) と命題 13.2(3) を適用すればよい。$`\square`$

*計算：* $`V_5`$ における定理 13.4：スペクトル分解、$`\phi(x)`$ の固有値、$`\phi(e^x)=e^{\phi(x)}`$、平方の錐 — [チュートリアル、13.2節](../../tutorial/ja/13_jordan_algebras.md#s-13-2)。*ライブラリでは：* `SpinFactor.spectral_decomposition`、`exp`、`sqrt`、`in_cone` — [チュートリアル、13.2節](../../tutorial/ja/13_jordan_algebras.md#s-13-2)。

<a id="r-13-5"></a>

**命題 13.5.** (a) $`m=3`$ のとき、$`\phi(a,u)=a1+u\cdot\sigma`$ はジョルダン同型
$`V_3\cong\mathrm{Herm}_2(\mathbb C)`$ である。(b) $`V_m`$ のジョルダン自己同型は写像 $`(a,u)\mapsto(a,Ou)`$、
$`O\in O(m)`$ である。(c) $`m\ge2`$ のとき、$`\phi(V_m)`$ の次元は $`m+1`$ であり、それが生成する $`*`$-代数の自己共役部分の次元は
$`D^2=4^{\lfloor m/2\rfloor}`$ である。両者が等しいのは $`m=3`$ のときに限る。

*証明.* (a) $`m=3`$ のとき $`(s_1,s_2,s_3)=(X,Y,Z)`$ であり、$`\phi`$ は次元 $`4`$ の空間の間の単射である。(b) 写像
$`(a,u)\mapsto(a,Ou)`$ は $`u\cdot v`$ を保ち、したがって $`\circ`$ を保つ。逆に $`\alpha`$ をジョルダン自己同型とする。すると
$`\alpha(1)`$ は $`V_m`$ の単位元であるから $`\alpha(1)=1`$ である。$`x^2`$ の公式により、$`x^2=1`$ となるのは $`x=\pm1`$ のとき、または
$`|w|=1`$ なる $`x=(0,w)`$ のときに限る。$`|w|=1`$ に対し $`\alpha(0,w)^2=1`$ かつ $`\alpha(0,w)\ne\pm1=\alpha(\pm1)`$ であるから、$`|w'|=1`$ なる $`w'`$ により $`\alpha(0,w)=(0,w')`$ となる。
線形な $`\ell,O`$ を用いて $`\alpha(0,u)=(\ell(u),Ou)`$ と書けば、これは $`\ell=0`$ かつ $`|Ou|=|u|`$、すなわち $`O\in O(m)`$ を意味する。
(c) 補題 13.3。$`m=2`$ および $`m\ge4`$ に対し $`4^{\lfloor m/2\rfloor}>m+1`$ である。$`\square`$

$`m=1`$ のとき、$`s_1=Z`$ は対角行列を生成し、その次元は $`2=m+1`$ であり、$`V_1\cong\mathbb R\oplus\mathbb R`$ は単純ではない。

*計算：* ランダムな $`O\in O(5)`$ に対する (b)、および `jordan.compare_with_star_algebra` による $`m\le6`$ に対する (c) — [チュートリアル、13.2節](../../tutorial/ja/13_jordan_algebras.md#s-13-2)。

**$`V_m`$ の状態.** $`V_m`$ の **状態** とは、$`\omega(Q)\ge0`$ かつ $`\omega(1)=1`$ を満たす線形汎関数 $`\omega`$ のことである。状態はちょうど、$`\rho=\frac12(1,r)`$、$`|r|\le1`$ を用いた写像 $`x\mapsto\mathrm{tr}(\rho\circ x)`$ である。したがって状態空間は次元 $`m`$ の球体である。*証明.* $`\langle x,y\rangle=\mathrm{tr}(x\circ y)=2(ab+u\cdot v)`$ は $`V_m`$ 上の内積であるから、すべての線形な $`\omega`$ はちょうど一つの $`\rho=(b,v)`$ により $`\mathrm{tr}(\rho\circ\cdot)`$ と表される（有限次元のリース表現：正規直交基底において $`\rho=\sum_k\omega(e_k)e_k`$）。定理 13.4(4) の自己双対性により、$`\omega(Q)\ge0`$ は $`\rho\in Q`$、すなわち $`b\ge|v|`$ と同値であり、また $`\omega(1)=\mathrm{tr}\,\rho=2b=1`$ である。よって $`r=2v`$、$`|r|\le1`$ として $`\rho=\frac12(1,r)`$ である。逆に、そのような $`\rho`$ はどれも状態を与える。$`\square`$

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. $`m=3`$ のとき、$`V_m`$ の状態空間は量子ビットのブロッホ球である（命題 13.5(a)）。ローレンツ錐は $`1+m`$ 次元の特殊相対論における未来の
> 光円錐であり、命題 13.5(b) は、ジョルダン自己同型が時間軸を固定する空間の直交変換であることを述べている。

錐 $`Q`$ は **二次錐計画**（SOCP）の制約集合である：すなわちローレンツ錐の直積上の線形最適化である。SOCP の内点法は通常、定理 13.4 の量（固有値、$`\det`$、そして $`P(y)`$）を用いて定式化される。これらのアルゴリズムは本講義では扱わない（参考文献：
F. Alizadeh, D. Goldfarb, *Second-order cone programming*, Math. Program. 95 (2003) 3–51）。

<a id="s-13-3"></a>

## 13.3 *-部分代数ではないジョルダン部分代数

ジョルダン部分代数 $`J\subseteq\mathrm{Herm}_n(\mathbb C)`$ は複素 $`*`$-代数 $`\mathcal A`$ を生成し、
$`J\subseteq\mathcal A_{sa}`$ である。命題 13.5(c) は $`J\ne\mathcal A_{sa}`$ となる例を与えた。ここではさらに二つの例を挙げる。

<a id="r-13-6"></a>

**命題 13.6.** (a) $`\mathrm{Herm}_n(\mathbb R)`$ は $`\mathrm{Herm}_n(\mathbb C)`$ の次元
$`n(n+1)/2`$ のジョルダン部分代数である。それが生成する $`*`$-代数は $`M_n(\mathbb C)`$ であり、その自己共役部分の次元は $`n^2`$ である。
(b) $`x,y\in\mathrm{Herm}_2(\mathbb C)`$ に対し、$`UxU^*`$ と $`UyU^*`$ がともに実対称となるユニタリ $`U`$ が存在する。$`1,x,y`$ が生成するジョルダン代数は $`\mathrm{span}\{1,x,y\}`$ である。$`x,y,1`$ が線形独立ならば、それは次元 $`3`$ の $`U^*\mathrm{Herm}_2(\mathbb R)U\cong V_2`$ であり、一方
$`1,x,y`$ が生成する $`*`$-代数は次元 $`4`$ の $`M_2(\mathbb C)`$ である。

*証明.* (a) 補題 13.1。$`i\ne j`$ に対し $`E_{ij}=E_{ii}(E_{ij}+E_{ji})`$ は実対称行列の積であるから、
生成される代数はすべての行列単位を含む。(b) $`x=a1+u\cdot\sigma`$、$`y=b1+v\cdot\sigma`$ と書き、単位ベクトル
$`n\perp u,v`$ を選ぶ。行列 $`n\cdot\sigma`$ と $`Y`$ はエルミートでトレース零、二乗すると $`1`$ であるから、ともに
固有値 $`1,-1`$ をもち、あるユニタリ $`U`$ に対し $`U(n\cdot\sigma)U^*=Y`$ となる（定理 1.3）。$`U`$ による共役は
$`\mathrm{Tr}`$ と $`\circ`$ を保つので、それは $`u\cdot\sigma`$ を $`(Ru)\cdot\sigma`$ に写す。ここで $`R`$ は線形で、
$`Ru\cdot Rv=\frac12\mathrm{Tr}\big((Ru\cdot\sigma)(Rv\cdot\sigma)\big)=u\cdot v`$ であるから $`R\in O(3)`$ であり、$`Rn=e_2`$ である。
よって $`Ru,Rv\perp e_2`$ であり、$`X,Z`$ が実であるから $`UxU^*=a1+(Ru)_1X+(Ru)_3Z`$ は実である。$`UyU^*`$ も同様である。次に
$`x\circ y=(ab+u\cdot v)1+av\cdot\sigma+bu\cdot\sigma\in\mathrm{span}\{1,x,y\}`$ であり、$`x^2,y^2`$ も同様であるから、この線形包は
閉じている。その次元が $`3`$ ならば、$`U`$ による共役は $`\mathrm{Herm}_2(\mathbb R)=\mathrm{span}\{1,X,Z\}`$ 全体であり、
これは $`(a,u)\mapsto a1+u_1X+u_2Z`$ を介して $`\cong V_2`$ であり（補題 13.3 の証明は $`s_i`$ たちが互いに反可換なエルミート対称元であることしか用いず、$`X,Z`$ はそのようなものである）、(a) により $`M_2(\mathbb C)`$ を生成する。$`\square`$

*計算：* $`n\le5`$ に対する命題 13.6(a)、およびランダムな $`2\times2`$ エルミート行列に対する (b)、`jordan.compare_with_star_algebra` を用いて — [チュートリアル、13.3節](../../tutorial/ja/13_jordan_algebras.md#s-13-3)。

一般の三つの $`2\times2`$ エルミート行列は $`\mathrm{Herm}_2(\mathbb C)`$ 全体を生成する。したがって $`\mathcal M`$ の $`M_2`$ ブロックの中では、二つのオブザーバブルが複素構造を明らかにすることは決してないが、三つ目のオブザーバブルはそれを明らかにしうる。

<a id="s-13-4"></a>

## 13.4 分類

有限次元のジョルダン代数は、$`\sum_ix_i^2=0`$ ならばすべての $`x_i=0`$ となるとき **形式的実** であるといい、単位元 $`1`$ と $`\langle x\circ y,z\rangle=\langle y,x\circ z\rangle`$ を満たす内積をもつとき **ユークリッド的** であるという。二つの概念は一致する（定理 E.5(4)。「ユークリッド的 $`\Rightarrow`$ 形式的実」の向きは一行で済む：$`0=\sum_i\langle x_i^2,1\rangle=\sum_i\langle x_i,x_i\rangle`$）。$`\langle x,y\rangle=\mathrm{Tr}(xy)`$ をもつ $`M_{sa}`$ はユークリッド的である：$`\mathrm{Tr}(x^2)=\sum_{ij}|x_{ij}|^2`$ であり、巡回性により
$`\mathrm{Tr}((x\circ y)z)=\frac12\mathrm{Tr}(xyz+yxz)=\mathrm{Tr}(y(x\circ z))`$ である。$`\langle x,y\rangle=\mathrm{tr}(x\circ y)=2(ab+u\cdot v)`$ をもつ $`V_m`$ も同様である：$`z=(c,w)`$ に対し
$`\mathrm{tr}((x\circ y)\circ z)=2(abc+c\,u\cdot v+a\,v\cdot w+b\,u\cdot w)`$ は $`x,y,z`$ について対称である。

<a id="r-13-7"></a>

**定理 13.7（Jordan–von Neumann–Wigner）.** 有限次元の形式的実ジョルダン代数はすべて単純イデアルの直和であり、単純なものはどれも次のうちちょうど一つと同型である：$`\mathbb R`$、
$`\mathbb K=\mathbb R,\mathbb C,\mathbb H`$ かつ $`r\ge3`$ に対する $`\mathrm{Herm}_r(\mathbb K)`$、スピン因子
$`V_m`$（$`m\ge2`$）、八元数上の $`3\times3`$ エルミート行列のなす **アルバート代数** $`\mathrm{Herm}_3(\mathbb O)`$（次元 $`27`$）。さらに
$`\mathbb K=\mathbb R,\mathbb C,\mathbb H,\mathbb O`$ に対し $`\mathrm{Herm}_1(\mathbb K)\cong\mathbb R`$ かつ $`\mathrm{Herm}_2(\mathbb K)\cong V_{\dim\mathbb K+1}`$ である。

*証明.* 付録 E がこれを完全に証明する：分解は定理 E.9、リストと「ちょうど一つ」は定理 E.27(1),(2)、階数 1 と 2 における同型は定理 E.24(3) である。$`\square`$

（P. Jordan, J. von Neumann, E. Wigner, *On an algebraic generalization of the quantum mechanical formalism*, Ann. of
Math. 35 (1934) 29–64。付録 E の証明は [FK94, Chapters II–V] に従う。）

アルバート代数を除くこれらすべては特殊である：$`\mathbb R=\mathrm{Herm}_1(\mathbb R)`$ と $`\mathrm{Herm}_r(\mathbb K)`$（$`\mathbb K=\mathbb R,\mathbb C`$）は補題 13.1 により、$`\mathrm{Herm}_r(\mathbb H)`$ は四元数を $`2\times2`$ 複素行列として書くことで $`\mathrm{Herm}_{2r}(\mathbb C)`$ に埋め込まれるから（定理 E.24(1)）、そして $`V_m`$ は補題 13.3 による。
アルバート代数はジョルダン代数であり（定理 E.24(2)）、**例外的** である：それはいかなる実または複素の結合代数のジョルダン部分代数とも同型でない
（定理 E.26。これはアルバート（A. A. Albert, *On a certain algebra of quantum mechanics*, Ann. of Math. 35 (1934) 65–73）が用いた種類の議論である：結合代数へのジョルダン埋め込みがあれば八元数の積が結合的になってしまうが、実際にはそうではない、補題 E.15(2)）。$`r\ge4`$ に対し $`\mathrm{Herm}_r(\mathbb O)`$ はジョルダン恒等式を満たさず（系 E.25）、$`r\ge4`$ 個の元からなるジョルダン枠をもつ単純代数は結合的な合成代数 $`D\ne\mathbb O`$ による $`\mathrm{Herm}_r(D)`$ である（定理 E.22）。これが八元数が $`r\le3`$ に対してのみ現れる理由である。したがって $`\mathrm{Herm}_3(\mathbb O)`$ は、積 $`\frac12(xy+yx)`$ をもつ実または複素行列によって（あるいはいかなる結合代数の中でも）忠実に表現することはできず、その中での計算には行列の積の代わりに八元数の積（またはその構造定数）を用いなければならない。

$`\Omega`$ をユークリッド空間 $`V`$（$`\ne\{0\}`$ である有限次元の実内積空間）内の空でない開凸錐とする。その **双対錐** は $`\Omega^*=\{y:\langle x,y\rangle>0\text{ すべての }x\in\bar\Omega\setminus\{0\}\text{ に対し}\}`$ である。$`\Omega^*=\Omega`$ のとき $`\Omega`$ は **自己双対** であるといい、$`g\Omega=\Omega`$ を満たす線形写像 $`g`$ が $`\Omega`$ 上に推移的に作用するとき **等質** であるといい、その両方を満たすとき **対称** であるという。

<a id="r-13-8"></a>

**定理 13.8（ケッヒャー–ヴィンベルク）.** ユークリッド空間 $`V`$ 内の空でない開凸錐 $`\Omega`$ が対称であるための必要十分条件は、$`(V,\circ,\langle\cdot,\cdot\rangle)`$ がユークリッド的ジョルダン代数となり、$`\Omega`$ がその平方の錐の内部となるような積 $`\circ`$ が $`V`$ 上に存在することである。

*証明.* 定理 F.12（「ならば」の向き：定理 F.5(2) と F.6(2)）。$`\square`$

（M. Koecher, 1957; E. B. Vinberg, 1960。付録 F の証明は [FK94, Chapters I–III] に従う。）定理 13.4(4) はローレンツ錐について自己双対性の閉じた形（$`y\in Q`$ であるための必要十分条件は、すべての $`x\in Q`$ に対し $`\mathrm{tr}(x\circ y)\ge0`$ となること）を証明し、演習問題 1 はその等質性を証明する。ユークリッド的ジョルダン代数の直和（内積の直交和をもつ）の平方の錐は、$`(x_1,\dots,x_k)^2=(x_1^2,\dots,x_k^2)`$ であるから、それぞれの平方の錐の直積である。したがって「ならば」の向きにより、半正定値計画の錐（$`\mathrm{Herm}_n(\mathbb R)^+`$ または $`\mathrm{Herm}_n(\mathbb C)^+`$ の直積。定理 1.6 によりこれらは平方の錐である。実の $`x\ge0`$ に対しては平方根も実である。その複素共役も $`x`$ の正の平方根だからである）、SOCP の錐（ローレンツ錐の直積）、そして線形計画の錐（$`\mathbb R_+^n`$、すなわち $`\mathbb R^n=\mathbb R\oplus\dots\oplus\mathbb R`$ の平方の錐）の内部は対称錐である。

[FK94]: J. Faraut, A. Korányi, *Analysis on Symmetric Cones*, Oxford University Press, 1994.

<a id="s-13-5"></a>

## 13.5 ジョルダン写像と *-写像

フォン・ノイマン環の **ジョルダン $`*`$-準同型** とは、$`\Phi(x^*)=\Phi(x)^*`$ と
$`\Phi(x\circ y)=\Phi(x)\circ\Phi(y)`$ を満たす複素線形写像 $`\Phi`$ のことである。自己共役部分の上ではそれはジョルダン準同型であり、自己共役部分の間のすべてのジョルダン
準同型は $`\Phi(x+iy)=\Phi(x)+i\Phi(y)`$（$`x,y`$ は自己共役）によってジョルダン $`*`$-準同型に拡張される。（この $`\Phi`$ は複素線形で $`\Phi(x-iy)=\Phi(x+iy)^*`$ を満たす。また $`\circ`$ は複素双線形で、自己共役元の組を自己共役元に写すから、$`(x+iy)\circ(x'+iy')=(x\circ x'-y\circ y')+i(x\circ y'+y\circ x')`$ は $`\Phi(x+iy)\circ\Phi(x'+iy')`$ に写される。）

<a id="r-13-9"></a>

**命題 13.9（転置）.** $`\mathcal M`$ 上で $`T(A)=\bigoplus_cA_c^T`$ とする。(a) $`T`$ は単位的なジョルダン
$`*`$-自己同型であり、$`T(AB)=T(B)T(A)`$ を満たす。$`T`$ はスペクトルを保つので、$`T`$ と $`T^{-1}=T`$ は正値である。(b) $`M_k`$
（$`k\ge2`$）上で、チョイ行列 $`\sum_{ij}E_{ij}\otimes T(E_{ij})`$ はスワップ $`F(e_a\otimes e_b)=e_b\otimes e_a`$ であり、その
固有値は $`+1`$（重複度 $`k(k+1)/2`$）と $`-1`$（重複度 $`k(k-1)/2`$）である。したがって $`T`$ は完全正値ではない。

*証明.* (a) $`(A^T)^*=(A^*)^T`$、$`(AB)^T=B^TA^T`$ であるから $`(A\circ B)^T=B^T\circ A^T=A^T\circ B^T`$ である。また
$`\det(A^T-\lambda)=\det(A-\lambda)`$ である。(b) $`\sum_{ij}E_{ij}e_a\otimes E_{ji}e_b=e_b\otimes e_a`$。$`F^2=1`$ であり、$`F`$ は対称テンソル上で $`+1`$、
反対称テンソル上で $`-1`$ であって、これらは空間全体を張る。定理 6.2（および命題 6.1）により $`T`$ は
CP ではない。$`\square`$

*計算：* $`M_3`$ 上の転置のチョイ行列による命題 13.9 — [チュートリアル、13.5節](../../tutorial/ja/13_jordan_algebras.md#s-13-5)。*ライブラリでは：* 転置 $`T`$ は `jordan.transpose_map` である — [チュートリアル、13.5節](../../tutorial/ja/13_jordan_algebras.md#s-13-5)。

したがってジョルダン自己同型は $`*`$-自己同型であるとは限らない。この現象を除けば、正値性がジョルダン構造を決定する：

<a id="r-13-10"></a>

**定理 13.10（カディソン、有限次元の場合）.** $`\mathcal M=\bigoplus_cM_{k_c}(\mathbb C)`$、
$`\mathcal N=\bigoplus_dM_{k'_d}(\mathbb C)`$ とし、$`\Phi:\mathcal M\to\mathcal N`$ を
$`\Phi(1)=1`$ を満たす線形全単射で、$`\Phi`$ と $`\Phi^{-1}`$ がともに正値であるものとする。このとき $`\Phi`$ はジョルダン $`*`$-同型である。逆に、ジョルダン $`*`$-同型 $`\mathcal M\to\mathcal N`$ はどれも単位的であり、それ自身とその逆はともに正値である。

*証明.* 定理 F.15。（考え方：$`\Phi`$ は凸集合 $`\{0\le a\le1\}`$ を $`\mathcal N`$ の対応する集合の上に写し、したがって端点、すなわち射影を射影に写す。$`\Phi`$ は $`p+q\le1`$ を保ち、これは射影に対しては $`pq=0`$ を意味する。ゆえにスペクトル分解をスペクトル分解に写し、$`\Phi(x^2)=\Phi(x)^2`$ となる。）
$`\square`$

カディソン（R. V. Kadison）はこの定理をすべての単位的 C\*-環に対して証明した（*A generalized Schwarz inequality and algebraic
invariants for operator algebras*, Ann. of Math. 56 (1952) 494–503）。本講義では有限次元の場合だけを証明し、用いる。因子環 $`M_n`$ については、ジョルダン自己同型は正確に知られている：
$`M_n(\mathbb C)`$ のすべてのジョルダン $`*`$-自己同型は、ユニタリ $`U`$ により $`A\mapsto UAU^*`$ または $`A\mapsto UA^TU^*`$ と表される
（定理 F.17、すべての $`n`$ について。この結果は N. Jacobson, C. E. Rickart, Trans. Amer. Math. Soc. 69 (1950) 479–502、および I. N. Herstein, Trans. Amer. Math. Soc. 81 (1956) 331–341 による）。前者の写像は乗法的、後者は反乗法的であるから、$`M_n`$ のジョルダン $`*`$-自己同型はすべての交換子を保つか、すべての交換子を反転させるかのいずれかである：$`\Phi([x,y])=\pm[\Phi(x),\Phi(y)]`$。$`\bigoplus_cM_{k_c}`$ 上では、ジョルダン $`*`$-同型は同じ大きさのブロックを置換し、各ブロック上に二つの仕方のいずれかで作用する（定理 F.18）。それは単位的で、両方向に正値であり（定理 13.10）、系 F.19 により、それが完全正値であるのは大きさ $`\ge2`$ のブロックをどれも転置しないとき、すなわち $`*`$-同型であるとき、かつそのときに限る。以下では $`n=2`$ の場合をスピン因子の言葉で直接証明する。

<a id="r-13-11"></a>

**命題 13.11.** $`\mathrm{Herm}_2(\mathbb C)`$ のすべてのジョルダン自己同型は、$`U\in U(2)`$ により $`A\mapsto UAU^*`$ または
$`A\mapsto UA^TU^*`$ と表される。

*証明.* 命題 13.5(a,b) により、ジョルダン自己同型は $`O\in O(3)`$ を用いて $`a1+u\cdot\sigma\mapsto a1+(Ou)\cdot\sigma`$ と作用する。
この意味での $`A\mapsto UAU^*`$ の行列を $`R_U\in O(3)`$ と書く。$`R_{UV}=R_UR_V`$ である。(i) $`\det R_U=1`$：
$`U=W\mathrm{diag}(e^{i\theta_1},e^{i\theta_2})W^*`$ は道
$`U_t=W\mathrm{diag}(e^{it\theta_1},e^{it\theta_2})W^*`$ によって $`1`$ と結ばれ、$`\det R_{U_t}\in\{\pm1\}`$ は $`t`$ に連続的に依存する。
(ii) すべての $`R\in SO(3)`$ はある $`R_U`$ である：命題 13.6 と同様に、
$`U_1ZU_1^*=(Re_3)\cdot\sigma`$、すなわち $`R_{U_1}e_3=Re_3`$ となる $`U_1`$ が存在する。すると $`R_{U_1}^{-1}R\in SO(3)`$ は $`e_3`$ を固定するので、$`e_3`$ のまわりのある角度 $`\beta`$ の回転である。直接計算により、$`U_2=\mathrm{diag}(e^{-i\beta/2},e^{i\beta/2})`$ に対し
$`U_2XU_2^*=\cos\beta\,X+\sin\beta\,Y`$ および $`U_2YU_2^*=-\sin\beta\,X+\cos\beta\,Y`$ となるので、$`R=R_{U_1U_2}`$ である。(iii)
転置は $`X,Z`$ を固定し $`Y\mapsto-Y`$ と写す：$`O_T=\mathrm{diag}(1,-1,1)`$、$`\det O_T=-1`$。よって
$`\det O=1`$ ならば $`O=R_U`$ であり、$`\det O=-1`$ ならば $`OO_T=R_U`$、すなわち自己同型は $`A\mapsto UA^TU^*`$ である。$`\square`$

*計算：* $`O_T=\mathrm{diag}(1,-1,1)`$ — [チュートリアル、13.5節](../../tutorial/ja/13_jordan_algebras.md#s-13-5)。

**複合系.** 一つの因子環のジョルダン構造は、$`T`$ と恒等写像を区別できない。二つの系に対しては、この曖昧さが目に見えるものになる：

<a id="r-13-12"></a>

**命題 13.12.** $`x,x'\in\mathrm{Herm}_n(\mathbb C)`$、$`y,y'\in\mathrm{Herm}_m(\mathbb C)`$ とする。
$`\mathbb C^n\otimes\mathbb C^m`$ に作用する $`\mathrm{Herm}_{nm}(\mathbb C)`$ において
```math
(x\otimes y)\circ(x'\otimes y')=(x\circ x')\otimes(y\circ y')+\tfrac14[x,x']\otimes[y,y'].
```
である。$`\theta=\mathrm{id}\otimes T`$（**部分転置**）、$`A\circ_\theta B=\theta(\theta(A)\circ\theta(B))`$ とおく。このとき
$`\circ_\theta`$ は $`\mathrm{Herm}_{nm}(\mathbb C)`$ 上のジョルダン積であって、$`\theta`$ を介して $`\circ`$ と同型であり、
局所オブザーバブル $`x\otimes1`$、$`1\otimes y`$ のすべての組の上で $`\circ`$ と一致し、
$`(x\otimes y)\circ_\theta(x'\otimes y')=(x\circ x')\otimes(y\circ y')-\frac14[x,x']\otimes[y,y']`$ を満たす。その平方の錐は
$`\theta(\mathrm{Herm}_{nm}^+)`$ である。$`n=m=2`$、$`\Omega=e_1\otimes e_1+e_2\otimes e_2`$ に対し、状態
$`\omega=\frac12|\Omega\rangle\langle\Omega|\ge0`$ はこの錐に属さない：$`\theta(\omega)=\frac12F`$ は固有値
$`-\frac12`$ をもつ。

*証明.* $`j=x\circ x'`$、$`c=\frac12[x,x']`$、$`k=y\circ y'`$、$`d=\frac12[y,y']`$ とおくと、$`xx'=j+c`$、$`x'x=j-c`$ であり、
$`y`$ についても同様であるから、$`\frac12(xx'\otimes yy'+x'x\otimes y'y)=\frac12\big((j+c)\otimes(k+d)+(j-c)\otimes(k-d)\big)=j\otimes k+c\otimes d`$ である。
$`\theta`$ は $`\mathrm{Herm}_{nm}`$ の実線形な対合である：それは $`x\otimes y`$ を $`x\otimes y^T`$ に写し、
$`x\otimes y`$ たちは $`\mathbb R`$ 上 $`\mathrm{Herm}_{nm}`$ を張る（$`\mathrm{Herm}_n,\mathrm{Herm}_m`$ の実基底の積 $`h_\alpha\otimes g_\beta`$ は $`M_{nm}`$ の複素基底をなし、エルミート行列はこの基底に関して実の係数をもつ）。したがって $`\circ_\theta`$ はジョルダン積であり、$`\theta`$ は同型である。$`(y^T\circ y'^T)^T=y\circ y'`$ かつ
$`[y^T,y'^T]^T=-[y,y']`$ であるから、$`(x\otimes y^T)\circ(x'\otimes y'^T)`$ の公式に $`\theta`$ を施すと
$`\circ_\theta`$ の公式が得られる。局所オブザーバブルの上では交換子の項は消える（$`[1,1]=0`$、$`[x,1]=0`$）。平方は
$`A\circ_\theta A=\theta(\theta(A)^2)`$ であり、$`\theta^2=\mathrm{id}`$ であるから、$`\omega`$ が $`\theta(\mathrm{Herm}_{nm}^+)`$ に属するのは $`\theta(\omega)\ge0`$ のときかつそのときに限る。
最後に $`\theta(\omega)=\frac12\sum_{ij}E_{ij}\otimes E_{ji}=\frac12F`$ である
（命題 13.9(b)。これは命題 6.1 の証明における計算である）。$`\square`$

*計算：* 積の公式、$`\circ_\theta`$、および $`\theta(\omega)`$ の固有値 $`-\frac12`$ — [チュートリアル、13.5節](../../tutorial/ja/13_jordan_algebras.md#s-13-5)。

二つの積は同じように良いジョルダン代数であり、一度に一つの部分系で測定できるものすべてについて一致するが、どの結合作用素が正であるか、すなわち複合系の状態が何であるかについては一致しない。両者を区別するのは交換子の項の符号であり、この符号は
$`*`$-代数 $`M_n\otimes M_m`$ によって定まり、因子のジョルダン代数によっては定まらない。（$`\circ_\theta`$ は
結合的な積 $`A\star B=\theta(\theta(A)\theta(B))`$ のジョルダン積であり、この積については $`(x\otimes y)\star(x'\otimes y')=xx'\otimes y'y`$ である：
すなわち代数 $`M_n\otimes M_m^{\rm op}`$ である。）錐 $`\theta(\mathrm{Herm}_{nm}^+)`$ こそが
**PPT 判定条件**（A. Peres, Phys. Rev. Lett. 77 (1996) 1413）を与えるものである。$`M_n\otimes M_m`$ の状態 $`\rho`$ は、$`p_k\ge0`$、$`a_k,b_k\ge0`$ により
$`\rho=\sum_kp_ka_k\otimes b_k`$（有限和）と書けるとき **分離可能** であるといい、そうでないとき **エンタングルしている** という。分離可能状態は両方の錐に属する。なぜなら
$`\theta(\rho)=\sum_kp_ka_k\otimes b_k^T`$ で $`b_k^T\ge0`$（命題 13.9(a)）であり、
$`a\otimes b=(a^{1/2}\otimes b^{1/2})^2\ge0`$ だからである。したがって $`\theta(\mathrm{Herm}_4^+)`$ に属さない $`\omega`$ は
エンタングルしている。（ジョルダン代数の複合についての参考文献。ここでは用いない：H. Hanche-Olsen, *On the structure and tensor
products of JC-algebras*, Canad. J. Math. 35 (1983) 1059–1074; H. Barnum, M. Graydon, A. Wilce,
*Composites and categories of Euclidean Jordan algebras*, Quantum 4 (2020) 359。）

**ダイナミクス.** 時間発展が結合的な積を強制するのではないかと期待するかもしれない。そうはならない：

<a id="r-13-13"></a>

**命題 13.13.** $`H\in\mathrm{Herm}_n(\mathbb C)`$、$`\delta(z)=i[H,z]`$ とする。(a) $`\delta`$ は
$`\mathrm{Herm}_n`$ をそれ自身に写し、$`\circ`$ の導分である。$`e^{t\delta}(z)=e^{itH}ze^{-itH}`$ はジョルダン
自己同型である。(b) $`x,y,z\in\mathrm{Herm}_n`$ に対し $`[L(x),L(y)]z=\frac14[[x,y],z]`$。(c) ある $`x_k,y_k\in\mathrm{Herm}_n`$ に対し $`\delta=4\sum_k[L(x_k),L(y_k)]`$ となる：ハミルトン力学はジョルダン積だけによって生成される。

*証明.* (a) $`(i[H,z])^*=-i[z,H]=i[H,z]`$ である。$`[H,\cdot]`$ は結合的な積の導分であり、したがって $`\circ`$ の導分でもある。
$`z(t)=e^{itH}ze^{-itH}`$ は $`z'=i[H,z]`$ と $`z(0)=z`$ を満たすので（命題 8.2(3)、$`\mathcal M=M_n`$ として）、補題 8.1 の一意性により $`z(t)=e^{t\delta}(z)`$ である。ユニタリ $`e^{itH}`$ による共役は積と随伴を保ち、したがって $`\circ`$ と $`\mathrm{Herm}_n`$ を保つ。
(b) $`4x\circ(y\circ z)=xyz+xzy+yzx+zyx`$ から、$`x,y`$ を入れ替えた同じ式を引く。(c) (b) により、
$`iH-\frac{\mathrm{Tr}\,H}n i1`$（これは $`iH`$ と同じ交換子をもつ）を $`[x_k,y_k]`$ の和として書けば十分である。$`i\ne j`$ と
$`a=E_{ij}+E_{ji}`$、$`b=i(E_{ij}-E_{ji})`$ に対し、$`[a,E_{ii}]=E_{ji}-E_{ij}`$、$`[b,E_{ii}]=-i(E_{ij}+E_{ji})`$、
$`[a,b]=2i(E_{jj}-E_{ii})`$ である。これらは $`\mathbb R`$ 上でトレース零の反エルミート行列全体を張る。$`\square`$

*計算：* ランダムな $`3\times3`$ 行列に対する (a) と (b) — [チュートリアル、13.5節](../../tutorial/ja/13_jordan_algebras.md#s-13-5)。

ジョルダン代数が含んでいないのは、オブザーバブルとそれが生成するダイナミクスの間の **対応** $`H\mapsto i[H,\cdot]`$ である。たとえば $`\mathrm{Herm}_n(\mathbb R)`$ では、導分
$`[L(x),L(y)]=\frac14[[x,y],\cdot]`$ は実反対称な $`[x,y]`$ から来るが、これは実対称行列の $`i`$ 倍ではない（演習問題 4）。この対応は $`\circ`$ によっては決まらない：$`\mathrm{Herm}_n(\mathbb C)`$ 上では $`H\mapsto i[H,\cdot]`$ と $`H\mapsto-i[H,\cdot]`$ がともにすべてのオブザーバブルを $`\circ`$ の導分に写し（命題 13.13(a)）、ジョルダン自己同型である転置はこの二つを入れ替える。なぜなら $`T(i[H,T(z)])=-i[H^T,z]`$ だからである。後者は反対代数（$`a\cdot_{\rm op}b=ba`$）の交換子の $`i`$ 倍である：$`-i[H,z]=i(H\cdot_{\rm op}z-z\cdot_{\rm op}H)`$——命題 13.12 の二つの符号である。（このような「力学的対応」についての参考文献。ここでは用いない：E. M. Alfsen, F. W. Shultz, *Geometry of
State Spaces of Operator Algebras*, Birkhäuser, 2003。）

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. スピンをもたない粒子に対しては、転置（実基底における複素共役）は **時間反転** である。
> それは一つの系のオブザーバブルの対称性であるが、反ユニタリであり、一つの部分系だけに施すと物理的な操作ではない：
> 命題 13.12 は、それがベル状態を負の固有値をもつ作用素に変えることを示している。量子力学が複素 $`*`$-代数を用いるのは、
> 複合系が対称化された積だけでなく、符号つきの交換子を必要とするからである。

<a id="s-13-6"></a>

## 13.6 トレース不等式に戻って

命題 11.9（実験 3）は、$`X=X^*`$ と $`Y\ge0`$ に対して $`\mathrm{Tr}|YXY|\le\mathrm{Tr}(Y|X|Y)`$ を述べている（ここで
$`X,Y`$ は $`\mathcal M`$ の作用素であり、パウリ行列ではない）。その証明は $`Y\ge0`$ を一度も用いていない：
$`T=X_\pm^{1/2}Y`$ による $`YX_\pm Y=T^*T`$ は $`Y=Y^*`$ しか必要とせず、補題 11.8 は $`YX_+Y`$ と $`YX_-Y`$ に適用される。したがってこの不等式はすべての $`Y=Y^*`$ に対して成り立ち、
ジョルダンの言葉では
```math
\mathrm{tr}|P(y)x|\le\mathrm{tr}\,P(y)|x|,
```
となる。これはトレースを一つ固定すれば、すべてのユークリッド的ジョルダン代数で意味をもつ。これを $`V_m`$ において、三つの性質だけを用いる議論によって証明する。

<a id="r-13-14"></a>

**命題 13.14（スピン因子における実験 3）.** $`x,y\in V_m`$ に対し
$`\mathrm{tr}|P(y)x|\le\mathrm{tr}\,P(y)|x|`$ であり、$`x\in Q`$ または $`-x\in Q`$ ならば等号が成り立つ。

*証明.* $`f_+(t)=\max(t,0)`$、$`f_-(t)=\max(-t,0)`$ として $`x_\pm=f_\pm(x)`$ とおく（定理 13.4(3)）：$`x_\pm\in Q`$、$`x=x_+-x_-`$、
$`|x|=x_++x_-`$ である。$`A=P(y)x_+`$、$`B=P(y)x_-`$ とおく。定理 13.4(5) によりこれらは $`Q`$ に属し、$`P(y)x=A-B`$ である。
$`P(y)x=\mu_1c_1+\mu_2c_2`$ をスペクトル分解とし（定理 13.4(2)。$`c_1=c_w`$、$`c_2=c_{-w}`$ であるから $`c_1+c_2=1`$、
$`c_1\circ c_2=0`$）、$`c`$ を $`\mu_i>0`$ なる $`c_i`$ の和とし、
$`t>0`$ のとき $`\mathrm{sgn}'(t)=1`$、それ以外のとき $`-1`$ として $`s=2c-1=\sum_i\mathrm{sgn}'(\mu_i)c_i`$ とおく。すると
$`s\circ P(y)x=\sum_i|\mu_i|c_i=|P(y)x|`$ であり（定理 13.4(2,3) により、
$`\mu_1=\mu_2`$ のときもこれは $`w`$ の選び方によらない）、$`c,1-c\in Q`$ である（冪等元は平方である）。自己双対性の易しい半分（定理 13.4(4)：
$`a,b\in Q`$ に対し $`\mathrm{tr}(a\circ b)\ge0`$）を $`\mathrm{tr}(A\circ(1-c))`$ と $`\mathrm{tr}(B\circ c)`$ に適用すると
```math
\mathrm{tr}|P(y)x|=\mathrm{tr}(A\circ c)-\mathrm{tr}(A\circ(1-c))-\mathrm{tr}(B\circ c)+\mathrm{tr}(B\circ(1-c))
\le\mathrm{tr}\,A+\mathrm{tr}\,B=\mathrm{tr}\,P(y)|x|.
```
を得る。$`x\in Q`$ ならば $`x_-=0`$、$`B=0`$、$`|P(y)x|=A`$ であるから、両辺はともに $`\mathrm{tr}\,A`$ に等しい。$`-x\in Q`$ のときも同様である。
$`\square`$

*計算：* $`m=2,4,7`$ に対する命題 13.14、および $`\phi`$ を通した $`Y=Y^*`$ に対する命題 11.9 の不等式 — [チュートリアル、13.6節](../../tutorial/ja/13_jordan_algebras.md#s-13-6)。

証明が用いるのは次のことだけである：(i) すべての $`y`$ に対し $`P(y)Q\subseteq Q`$。(ii) スペクトル定理：すべての $`z`$ は、和が $`1`$ となる直交冪等元 $`c_i`$ を用いて
$`\sum_i\mu_ic_i`$ と書け、$`f(z)=\sum_if(\mu_i)c_i`$ である（したがって $`x_\pm`$ と冪等元は平方である）。(iii) $`\mathrm{tr}`$ は線形であり、$`\mathrm{tr}(a\circ1)=\mathrm{tr}\,a`$ であり、$`a,b\in Q`$ に対し $`\mathrm{tr}(a\circ b)\ge0`$ である。これら三つは、単位元 $`e`$ をもつすべてのユークリッド的ジョルダン代数 $`V`$ において線形形式 $`\tau(x)=\langle x,e\rangle`$ について成り立ち、この形式については $`\tau(a\circ b)=\langle b,a\circ e\rangle=\langle a,b\rangle`$ である：(ii) は定理 F.3、(i) は定理 F.6(4)、そして $`a,b\in Q`$ に対する $`\langle a,b\rangle\ge0`$ は定理 F.5(1) である。系 F.7 がこの議論を実行する：すべてのユークリッド的ジョルダン代数において $`\tau(|P(y)x|)\le\tau(P(y)|x|)`$ であり、$`x\in Q`$ または $`-x\in Q`$ ならば等号が成り立つ。$`\langle x,y\rangle=\mathrm{Tr}(xy)`$ をもつ $`M_{sa}`$ 上では $`\tau=\mathrm{Tr}`$、$`P(y)x=yxy`$ である。$`(M_{sa},\circ)`$ の冪等元は射影であり（$`c\circ c=c^2`$）、射影に対しては $`c\circ c'=0`$ と $`cc'=0`$ は同値である（補題 F.14(2)）から、定理 F.3 の一意性により、$`x`$ のその分解は固有射影へのスペクトル分解であり、$`|x|`$ は行列の絶対値である。したがってこれはすべての $`Y=Y^*`$ に対して命題 11.9 を再証明する。$`\langle x,y\rangle=\mathrm{tr}(x\circ y)`$ をもつ $`V_m`$ 上では $`\tau=\mathrm{tr}`$ であり、定理 F.3 の分解は定理 13.4(2) の分解（$`u=0`$ のときは $`x=a1`$）であって、これが命題 13.14 である。階数 $`r`$ の単純ユークリッド的ジョルダン代数 $`J`$ において、系 E.28 のトレース $`\mathrm{tr}\,x=\frac r{\dim J}\mathrm{Tr}\,L(x)`$ はすべてのジョルダン枠 $`(c_i)`$ に対し $`\mathrm{tr}(\sum_i\mu_ic_i)=\sum_i\mu_i`$ を満たし（定理 E.7）、$`\langle x,y\rangle=\mathrm{tr}(x\circ y)`$ は結合的な内積である。これについては $`\tau=\mathrm{tr}`$ であるから、不等式はこのトレースについて成り立つ。これは $`\mathrm{Herm}_r(\mathbb H)`$ とアルバート代数に適用される。これらは単純ユークリッド的ジョルダン代数である（定理 E.24）。$`r\ge3`$ のとき、これらはどの $`\mathcal M=\bigoplus_cM_{k_c}`$ の $`M_{sa}`$ ともジョルダン同型でない：$`M_{sa}=\bigoplus_c\mathrm{Herm}_{k_c}(\mathbb C)`$ はイデアルの直和であるから、単純になるのはブロックが一つだけのときに限り、そして $`\mathrm{Herm}_k(\mathbb C)`$（$`k=1`$ では $`\cong\mathbb R`$、$`k=2`$ では $`\cong V_3`$）は、定理 13.7 の一意性により、$`\mathrm{Herm}_r(\mathbb H)`$（$`r\ge3`$）とも $`\mathrm{Herm}_3(\mathbb O)`$ とも同型でない。$`V_m`$ では、$`\phi(P(y)x)=\phi(y)\phi(x)\phi(y)`$ かつ $`\phi(|x|)=|\phi(x)|`$ であるから、この不等式は $`\phi`$ を通して命題 11.9（上で $`Y=Y^*`$ に拡張したもの）からも従う。

<a id="s-13-7"></a>

## 13.7 比較の対象

*ライブラリでは：* 他のソフトウェア（SageMath、Magma、SOCP ソルバー）との比較と、ライブラリが付け加えるもの — [チュートリアル、13.7節](../../tutorial/ja/13_jordan_algebras.md#s-13-7)。

<a id="x-13"></a>

## 演習問題

1. （**ローレンツ錐の等質性.**）$`y\in V_m`$、$`\det y\ne0`$ とする。$`y^{-1}=(a,-u)/\det y`$ が
   $`y\circ y^{-1}=1`$ を満たすこと、$`P(y)P(y^{-1})=\mathrm{id}`$ であること、$`P(y)`$ が $`Q`$ の内部をそれ自身の上に写すこと、そして内部の $`x`$ に対し
   $`P(x^{1/2})1=x`$ であることを示せ。$`Q`$ の内部が等質であることを導け。
2. （**等方的状態.**）$`p\in[0,1]`$ に対し、$`\mathbb C^2\otimes\mathbb C^2`$ 上で $`\rho_p=p\,\omega+(1-p)\frac14`$ とする
   （命題 13.12）。$`\theta(\rho_p)`$ のスペクトルを計算し、$`p>\frac13`$ のとき $`\rho_p`$ がエンタングルしていることを示せ。
3. （**境界点の積.**）$`\mathbb R^m`$ において $`|w|=|w'|=1`$ とする。$`(1,w)\circ(1,w')\in Q`$ であるための必要十分条件が
   $`w'=w`$ または $`w'=-w`$ であることを示せ：ローレンツ錐の二本の境界半直線のジョルダン積はほとんど決して錐に属さない。
4. （**実量子力学.**）$`x,y\in\mathrm{Herm}_n(\mathbb R)`$ に対し、$`\mathrm{Herm}_n(\mathbb R)`$ の導分
   $`[L(x),L(y)]`$ は $`K`$ を実反対称行列として $`z\mapsto[K,z]`$ であること、そして
   $`n\ge2`$ のとき、そのような写像は $`0`$ でない限り決して $`H\in\mathrm{Herm}_n(\mathbb R)`$ による $`z\mapsto i[H,z]`$ にはならないことを示せ。

*計算：* 演習問題 1–3 の数値的検算 — [チュートリアル、演習問題](../../tutorial/ja/13_jordan_algebras.md#x-13)。

<details>
<summary>解答</summary>

1. $`(a,u)\circ(a,-u)=(a^2-|u|^2,0)`$ であり、これから $`y\circ y^{-1}=1`$ を得る。$`\phi(y)`$ と $`\phi(y^{-1})`$ は可換であり（ともに
   $`\mathrm{span}\{1,u\cdot s\}`$ に属する）、したがって $`\phi(y)\phi(y^{-1})=\phi(y)\circ\phi(y^{-1})=1`$ であり、定理 13.4(5) により
   $`\phi(P(y)P(y^{-1})x)=\phi(x)`$ である。$`\phi`$ は単射である。$`P(y)`$ は $`Q`$ を $`Q`$ の中に写し、その逆 $`P(y^{-1})`$ も同様であるから、
   $`P(y)Q=Q`$ である。線形同相写像であるから、それは内部を内部の上に写す。最後に
   $`P(y)1=2y\circ y-y^2=y^2`$ であるから $`P(x^{1/2})1=x`$ であり、$`x`$ が内部の点ならば $`P(x^{1/2})`$ は可逆である
   （$`\det x^{1/2}=\sqrt{\lambda_+\lambda_-}>0`$。$`a-|u|`$ は連続であり、$`a=|u|`$ のとき $`(a-\varepsilon,u)\notin Q`$ であるから、$`Q`$ の内部は $`\{a>|u|\}=\{\lambda_->0\}`$ である）。内部の任意の二点は
   $`P(x'^{1/2})P(x^{1/2})^{-1}`$ によって結ばれる。

2. 命題 13.12 により $`\theta(\rho_p)=\frac p2F+\frac{1-p}4`$ である。$`F`$ は固有値 $`1`$（三重）と $`-1`$（一重）をもつので、
   スペクトルは $`\frac{1+p}4`$（三重）と $`\frac{1-3p}4`$ であり、後者が負となるのは $`p>\frac13`$ のときかつそのときに限る。そのとき PPT 判定条件により $`\rho_p`$ は分離可能でない。

3. $`t=w\cdot w'\in[-1,1]`$ として $`(1,w)\circ(1,w')=(1+t,\,w+w')`$ であり、$`|w+w'|=\sqrt{2+2t}`$ である。これが $`Q`$ に属するのは
   $`1+t\ge\sqrt{2(1+t)}`$ のとき、すなわち $`t=-1`$（積は $`0`$）または $`1+t\ge2`$、すなわち $`t=1`$ のときかつそのときに限る。

4. 命題 13.13(b) により $`K=\frac14[x,y]`$ として $`[L(x),L(y)]z=[K,z]`$ であり、$`K^T=\frac14[y^T,x^T]=-K`$ は実である。実の $`H`$
   と $`z`$ に対し $`i[H,z]`$ は純虚数であり、一方 $`[K,z]`$ は実であるから、両者が一致するのはともに零となるときに限る。したがって実量子力学では、ダイナミクスの生成元 $`K`$ はオブザーバブルではない。

</details>

**この先へ.** 本章は本講義の最終章であり、付録 E と F がその証明を補完する。さらに学ぶために：有限次元の理論の標準的な文献は [FK94] である。K. McCrimmon, *A Taste of Jordan Algebras* (Springer 2004) は代数的側面と歴史を扱い、H. Hanche-Olsen and E. Størmer, *Jordan Operator Algebras* (Pitman 1984) と Alfsen–Shultz (2003) は無限次元の理論と C\*-環との関係を扱う。操作的公理からの量子論の再構成における複合系の役割については、上で引用した Barnum と Wilce の論文を参照せよ。第11章を越えるトレース不等式については E. Carlen, *Trace inequalities and quantum entropy* (Contemp. Math. 529, 2010) を参照せよ。二次表現とローレンツ錐は錐最適化に再び現れる：Alizadeh–Goldfarb (2003) と、そこで概観されているジョルダン代数的内点法を参照せよ。
