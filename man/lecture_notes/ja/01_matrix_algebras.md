[English](../01_matrix_algebras.md) · [Español](../es/01_matrix_algebras.md) · [Français](../fr/01_matrix_algebras.md) · [Deutsch](../de/01_matrix_algebras.md) · [中文](../zh/01_matrix_algebras.md) · **日本語** · [Русский](../ru/01_matrix_algebras.md)

# 1. 行列環と関数計算

本章では記号を定め、以降のすべての章で用いる $`M_n(\mathbb C)`$ に関する事実、すなわちスペクトル定理、関数計算とその微分、極分解、および特異値によって定義されるノルムを証明する。これらの結果のライブラリのコードと数値的な検証は[チュートリアル 1](../../tutorial/ja/01_matrix_algebras.md)にある。

**規約.** 内積 $`\langle x,y\rangle=\sum_i\bar x_iy_i`$ は第二引数について線形とする。$`A^*`$ は共役転置であり、したがって $`\langle x,Ay\rangle=\langle A^*x,y\rangle`$ である。$`e_1,\dots,e_n`$ は標準基底、$`1`$ は単位行列、$`\mathrm{Tr}`$ は通常のトレースを表す。

<a id="s-1-1"></a>

## 1.1 環 $`M_n(\mathbb C)`$

$`M_n(\mathbb C)`$ は複素 $`n\times n`$ 行列のなす環であり、対合 $`A\mapsto A^*`$ と作用素ノルム $`\|A\|=\sup_{\|x\|=1}\|Ax\|`$ を備える。元 $`A`$ は次のように呼ばれる。

| 名称 | 定義 |
|---|---|
| 正規 | $`A^*A=AA^*`$ |
| 自己共役 | $`A=A^*`$ |
| 正、$`A\ge0`$ | すべての $`x`$ に対し $`\langle x,Ax\rangle\ge0`$ |
| ユニタリ | $`A^*A=AA^*=1`$ |
| 射影 | $`A=A^*=A^2`$ |

$`B-A\ge0`$ のとき $`A\le B`$ と書き、$`A\ge0`$ かつ $`A`$ が可逆のとき $`A>0`$ と書く。

*ライブラリでは：* 指定したスペクトルをもつ自己共役または正の作用素、`TypeIAlgebra.operator_from_eigenvalues` — [チュートリアル、1.5節](../../tutorial/ja/01_matrix_algebras.md#s-1-5)；作用素ノルム、`operator_norm` — [チュートリアル、1.6節](../../tutorial/ja/01_matrix_algebras.md#s-1-6)。

<a id="r-1-1"></a>

**補題 1.1（C\*-恒等式）.** $`\|A^*\|=\|A\|`$ かつ $`\|A^*A\|=\|A\|^2`$ である。

*証明.* $`\|Ax\|=\sup_{\|y\|=1}|\langle y,Ax\rangle|=\sup_{\|y\|=1}|\langle A^*y,x\rangle|`$
であるから、$`\|A\|=\sup_{\|x\|=\|y\|=1}|\langle A^*y,x\rangle|=\|A^*\|`$ である。次に、$`\|x\|=1`$ に対し
$`\|Ax\|^2=\langle x,A^*Ax\rangle\le\|A^*A\|`$ であるから、$`\|A\|^2\le\|A^*A\|\le\|A^*\|\,\|A\|=\|A\|^2`$ を得る。$`\square`$

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. 量子力学では、$`n`$ 個の区別可能な配置をもつ系は $`\mathbb C^n`$ によって記述される。自己共役行列は
> *オブザーバブル*（エネルギー、スピン、格子上の位置）であり、その固有値は測定の取りうる結果である。射影は yes/no 型の問いである。ユニタリ
> 行列は対称性および時間発展であり、ハミルトニアン $`H`$ に対しては $`e^{-itH}`$ である。
> トレース 1 の正行列は *状態*（密度行列、第3章）である。

<a id="s-1-2"></a>

## 1.2 スペクトル定理

<a id="r-1-2"></a>

**補題 1.2（シューア）.** 任意の $`A\in M_n(\mathbb C)`$ に対し、$`U^*AU`$ が上三角となるユニタリ $`U`$ が存在する。

*証明.* $`n`$ に関する帰納法による。$`n=1`$ のときは自明である。$`A`$ の特性多項式は根 $`\lambda`$ をもつ（代数学の基本定理）ので、ある単位ベクトル $`v`$ に対し $`Av=\lambda v`$ である。$`v`$ を正規直交基底 $`v,w_2,\dots,w_n`$ に延長し、これらを列とするユニタリを
$`V`$ とする。すると $`A'\in M_{n-1}(\mathbb C)`$ により
$`V^*AV=\begin{pmatrix}\lambda&*\\0&A'\end{pmatrix}`$ と書ける。
帰納法の仮定より、あるユニタリ $`U'`$ に対して $`U'^*A'U'`$ は上三角であり、
$`U=V\,(1\oplus U')`$ が求めるものである。$`\square`$

<a id="r-1-3"></a>

**定理 1.3（スペクトル定理）.** $`A`$ が正規であるための必要十分条件は、あるユニタリ $`U`$ により
$`A=U\,\mathrm{diag}(\lambda_1,\dots,\lambda_n)\,U^*`$ と書けることである。
$`\lambda_i`$ は重複度込みの $`A`$ の固有値である。

*証明.* $`\Lambda`$ を対角行列として $`A=U\Lambda U^*`$ ならば、
$`A^*A=U\bar\Lambda\Lambda U^*=U\Lambda\bar\Lambda U^*=AA^*`$ である。逆に $`A`$
を正規とし、$`T=U^*AU`$ を上三角とする（補題 1.2）。$`T`$ は正規である。
$`(1,1)`$ 成分を比較すると $`(T^*T)_{11}=|t_{11}|^2`$ かつ $`(TT^*)_{11}=\sum_j|t_{1j}|^2`$ であるから、
$`j>1`$ に対し $`t_{1j}=0`$ である。よって $`T=t_{11}\oplus T'`$（$`T'`$ は正規かつ上三角）となり、帰納法により $`T`$ は対角である。$`U^*AU`$ の対角成分は
$`\det(A-\lambda)=\det(T-\lambda)`$ の根からなる。$`\square`$

*ライブラリでは：* 定理 1.3 を逆向きに読んだもの、すなわち指定した固有値とハール分布に従うランダムな $`U`$ による $`A=U\Lambda U^*`$ が `TypeIAlgebra.operator_from_eigenvalues` である — [チュートリアル、1.5節](../../tutorial/ja/01_matrix_algebras.md#s-1-5)。

<a id="r-1-4"></a>

**系 1.4.** $`A`$ を固有値 $`\lambda_i`$ をもつ正規行列とする。このとき

1. $`\|A\|=\max_i|\lambda_i|`$；
2. $`A`$ が自己共役 $`\iff`$ すべての $`\lambda_i\in\mathbb R`$；ユニタリ $`\iff`$ すべての $`|\lambda_i|=1`$；
   射影 $`\iff`$ すべての $`\lambda_i\in\{0,1\}`$；
3. $`A\ge0`$ $`\iff`$ $`A`$ が自己共役かつすべての $`\lambda_i\ge0`$。さらに、任意の正行列は自己共役であり、したがって正規であるから、
   (3) は *すべての* 行列の中で正値性を特徴づけ、以下の関数計算は任意の $`A\ge0`$ に適用できる。

*証明.* $`y=U^*x`$ とおくと $`\|Ax\|^2=\sum_i|\lambda_i|^2|y_i|^2`$ かつ $`\|y\|=\|x\|`$
であり、これより (1) が従う。等号は $`|\lambda_i|`$ が最大となる $`i`$ について $`x=Ue_i`$ のとき成り立つ。(2) では $`A^*=U\bar\Lambda U^*`$、
$`A^*A=U|\Lambda|^2U^*`$、$`A^2=U\Lambda^2U^*`$ であり、各条件は $`\Lambda`$ について成分ごとに成り立つときかつそのときに限り成り立つ。(3) について、まず任意の正の $`A`$ は自己共役である。実際、すべての $`x`$ に対し
$`\langle x,(A-A^*)x\rangle=\langle x,Ax\rangle-\overline{\langle x,Ax\rangle}=0`$ であり、
すべての $`x`$ に対し $`\langle x,Dx\rangle=0`$ を満たす行列 $`D`$ は零である（偏極恒等式：
$`4\langle x,Dy\rangle=\sum_{k=0}^3 i^{-k}\langle x+i^ky,D(x+i^ky)\rangle`$）。すると
$`\langle x,Ax\rangle=\sum_i\lambda_i|y_i|^2`$ であり、これがすべての $`x`$ に対して $`\ge0`$ となるのはすべての $`\lambda_i\ge0`$ のときかつそのときに限る。$`\square`$

<a id="s-1-3"></a>

## 1.3 関数計算

$`A=U\Lambda U^*`$ を正規とし、$`f`$ を $`\mathrm{spec}\,A`$ 上の複素数値関数とする。
```math
f(A)=U\,\mathrm{diag}\big(f(\lambda_1),\dots,f(\lambda_n)\big)\,U^*.
```
と定義する。

<a id="r-1-5"></a>

**定理 1.5.** (a) $`f(A)`$ は $`U`$ の選び方によらない。(b) 多項式 $`p`$ に対し、$`p(A)`$ は $`A`$ の通常の多項式である。(c) $`f\mapsto f(A)`$ は線形かつ乗法的であり、$`\bar f(A)=f(A)^*`$ かつ $`\mathrm{spec}\,f(A)=f(\mathrm{spec}\,A)`$ である。
(d) $`B`$ が $`A`$ と可換ならば、$`B`$ は $`f(A)`$ と可換である。(e) $`g`$ が
$`f(\mathrm{spec}\,A)`$ 上で定義されていれば、$`g(f(A))=(g\circ f)(A)`$ である。

*証明.* (b) $`A^k=U\Lambda^kU^*`$ であるから、線形性により $`p(A)=Up(\Lambda)U^*`$ である。(a) スペクトルは有限集合であるから、ラグランジュ補間により
$`\mathrm{spec}\,A`$ 上で $`p=f`$ となる多項式 $`p`$ が存在する：$`\mu_1,\dots,\mu_r`$ を相異なる固有値として
$`p(t)=\sum_if(\mu_i)\prod_{j\neq i}(t-\mu_j)/(\mu_i-\mu_j)`$ ととればよい。すると $`A`$ を対角化する *任意の* ユニタリ $`U`$ に対し $`U\,f(\Lambda)\,U^*=U\,p(\Lambda)\,U^*=p(A)`$ であり、$`p(A)`$ は $`U`$ を用いずに定義される。(c) は対角行列 $`f(\Lambda)`$ について成り立ち、$`U`$ による共役で保たれる。
$`Uf(\Lambda)U^*`$ のスペクトルは $`f(\Lambda)`$ の対角成分である。(d) $`B`$ は $`A`$ の任意の多項式と可換であり、特に $`p(A)=f(A)`$ と可換である。(e) $`f(A)=Uf(\Lambda)U^*`$ は正規であり同じ $`U`$ で対角化されるから、(a) により $`g(f(A))=U\,g(f(\Lambda))\,U^*=(g\circ f)(A)`$ である。$`\square`$

<a id="r-1-6"></a>

**定理 1.6（平方根）.** $`A\ge0`$ に対し、
$`B^2=A`$ を満たす $`B\ge0`$ はちょうど一つ存在し、それは関数計算で定義される $`B=\sqrt A`$ である。

*証明.* 定理 1.5(c) と系 1.4(3) により $`\sqrt A\ge0`$ かつ $`(\sqrt A)^2=A`$ である。
$`B\ge0`$ が $`B^2=A`$ を満たすとし、$`B=V\,\mathrm{diag}(\mu_j)\,V^*`$、$`\mu_j\ge0`$ とする。すると
$`A=V\,\mathrm{diag}(\mu_j^2)\,V^*`$ であるから $`\mathrm{spec}\,A=\{\mu_j^2\}`$ である。$`p`$ を
$`\mathrm{spec}\,A`$ 上で $`p(t)=\sqrt t`$ となる多項式とする。すると
$`p(A)=\sqrt A`$ かつ $`p(A)=p(B^2)=V\,\mathrm{diag}(p(\mu_j^2))\,V^*=V\,\mathrm{diag}(\mu_j)\,V^*=B`$ である。$`\square`$

*計算：* `apply_function`、`sqrt`、`log`、`expm`、`power` による定理 1.5 と定理 1.6 の検証 — [チュートリアル、1.5節](../../tutorial/ja/01_matrix_algebras.md#s-1-5)。

<a id="r-1-7"></a>

**系 1.7.** $`A\ge0`$ $`\iff`$ ある $`T`$ により $`A=T^*T`$。（$`A\ge0`$ ならば $`T=\sqrt A`$ ととればよい。
逆に $`\langle x,T^*Tx\rangle=\|Tx\|^2\ge0`$ である。）

**微分.** 勾配法では $`f(A)`$ を $`A`$ について微分する。自己共役な $`A`$ に対しては、微分は $`f`$ の差分商からなる行列との成分ごとの積である。

<a id="r-1-11"></a>

**定理 1.11（ダレツキー–クレイン）.** $`A=U\,\mathrm{diag}(\lambda_1,\dots,\lambda_n)\,U^*`$ を自己共役、$`I`$ を $`\mathrm{spec}\,A`$ を含む開区間、$`f:I\to\mathbb C`$ を連続微分可能とする。
```math
f^{[1]}(x,y)=\frac{f(x)-f(y)}{x-y}\quad(x\ne y),\qquad f^{[1]}(x,x)=f'(x)
```
を $`f`$ の**差分商**（divided difference）、$`L_{ij}=f^{[1]}(\lambda_i,\lambda_j)`$ を**レヴナー行列**とする。このとき任意の自己共役な $`E`$ に対して
```math
\frac{d}{dt}\Big|_{t=0}f(A+tE)=U\big(L\circ U^*EU\big)U^*
```
が成り立つ。ここで $`(X\circ Y)_{ij}=X_{ij}Y_{ij}`$ は成分ごとの（シューア）積である。

*証明.* 微積分の基本定理により、すべての $`x,y\in I`$ に対して $`f^{[1]}(x,y)=\int_0^1f'(y+s(x-y))\,ds`$ である（$`x=y`$ では両辺とも $`f'(x)`$）。よって $`f^{[1]}`$ は $`I\times I`$ 上で連続である。$`U`$ の列を並べ替えると $`L`$ と $`U^*EU`$ の行と列が同じように置換され、右辺は変わらないので、$`\lambda_1\le\dots\le\lambda_n`$ と仮定してよい。$`[\lambda_1-\delta,\lambda_n+\delta]\subset I`$ となる $`\delta>0`$ をとる。$`0<|t|\,\|E\|<\delta`$ に対して $`A+tE=V\,\mathrm{diag}(\mu_1,\dots,\mu_n)\,V^*`$（$`\mu_1\le\dots\le\mu_n`$）と書く。すると $`|\mu_j-\lambda_j|\le|t|\,\|E\|`$（定理 A.10）であり、すべての $`\mu_j\in I`$ である。$`W=U^*V`$（ユニタリ）、$`\tilde E=U^*EU`$ とおく。$`U^*f(A)=\mathrm{diag}(f(\lambda_i))\,U^*`$ かつ $`f(A+tE)V=V\,\mathrm{diag}(f(\mu_j))`$ であるから
```math
\big(U^*(f(A+tE)-f(A))V\big)_{ij}=\big(f(\mu_j)-f(\lambda_i)\big)W_{ij}=f^{[1]}(\lambda_i,\mu_j)\,(\mu_j-\lambda_i)\,W_{ij}
```
であり、$`f(x)=x`$ に対する同じ計算から $`(\mu_j-\lambda_i)W_{ij}=t\,(\tilde EW)_{ij}`$ を得る。$`(L_t)_{ij}=f^{[1]}(\lambda_i,\mu_j)`$ とおけば、これは厳密な恒等式
```math
\frac1t\big(f(A+tE)-f(A)\big)=U\big((L_t\circ\tilde EW)\,W^*\big)U^*
```
である。$`t\to0`$ のとき、$`f^{[1]}`$ が連続で $`\mu_j\to\lambda_j`$ であるから、$`L_t\to L`$（成分ごと）。$`\lambda_i=\lambda_j`$ ならば $`D_{ij}=W_{ij}`$、そうでなければ $`D_{ij}=0`$ とし、$`R=W-D`$ とおく。$`\lambda_i\ne\lambda_j`$ かつ $`|t|\,\|E\|\le|\lambda_i-\lambda_j|/2`$ ならば、$`|\mu_j-\lambda_i|\ge|\lambda_i-\lambda_j|/2`$ かつ $`|(\tilde EW)_{ij}|\le\|\tilde EW\|=\|E\|`$ なので $`|W_{ij}|\le2|t|\,\|E\|/|\lambda_i-\lambda_j|`$ である。したがって $`R\to0`$ であり、$`DD^*=(W-R)(W-R)^*\to WW^*=1`$。$`\lambda_k=\lambda_j`$ ならば $`L_{ij}=L_{ik}`$ であり、$`\lambda_k=\lambda_j`$ でなければ $`D_{kj}=0`$ であるから、$`(L\circ\tilde ED)_{ij}=\sum_kL_{ik}\tilde E_{ik}D_{kj}=((L\circ\tilde E)D)_{ij}`$。ゆえに
```math
(L_t\circ\tilde EW)W^*-L\circ\tilde E=\big((L_t-L)\circ\tilde EW\big)W^*+(L\circ\tilde ER)W^*+(L\circ\tilde ED)R^*+(L\circ\tilde E)(DD^*-1)
```
であり、$`W`$、$`D`$、$`\tilde E`$ の成分はそれぞれ $`1`$、$`1`$、$`\|E\|`$ で抑えられるので、各項は $`0`$ に収束する。$`\square`$

実数値の $`f`$ に対しては $`L`$ は実行列であり、任意の $`G`$ に対して $`\mathrm{Tr}\big(G^*U(L\circ U^*EU)U^*\big)=\mathrm{Tr}(\Gamma^*E)`$、$`\Gamma=U(L\circ U^*GU)U^*`$ が成り立つ（基底 $`U`$ で両辺を $`i,j`$ についての和として書けばよい）。すなわち $`A\mapsto\mathrm{Re}\,\mathrm{Tr}(G^*f(A))`$ の方向 $`E`$ への微分は $`\mathrm{Re}\,\mathrm{Tr}(\Gamma^*E)`$ である。

*ライブラリでは：* `apply_function` は定理 1.11 によって（autograd で）微分される。レヴナー行列は `functional.divided_differences` である。

**丸め.** $`\lambda_i,\lambda_j`$ が近いとき、商 $`(f(\lambda_i)-f(\lambda_j))/(\lambda_i-\lambda_j)`$ は桁落ちにより精度を失い、$`\lambda_i=\lambda_j`$ では定義されない。このときは導関数のほうがよいデータである。

<a id="r-1-12"></a>

**命題 1.12（中点値）.** $`a<b`$、$`h=b-a`$ とし、$`f`$ を $`[a,b]`$ 上で 3 回連続微分可能とする（複素数値でもよい）。このとき
```math
f^{[1]}(a,b)-\frac{f'(a)+f'(b)}2=-\frac{h^2}2\int_0^1s(1-s)\,f'''(a+sh)\,ds,
\qquad
\Big|f^{[1]}(a,b)-\frac{f'(a)+f'(b)}2\Big|\le\frac{h^2}{12}\max_{[a,b]}|f'''|
```
であり、$`3`$ 次多項式に対して等号が成り立つ。

*証明.* $`g(s)=f'(a+sh)`$ とおくと $`f^{[1]}(a,b)=\frac1h\int_a^bf'=\int_0^1g`$ であり、$`\varphi(s)=s(1-s)/2`$ とおくと $`\varphi(0)=\varphi(1)=0`$、$`\varphi'(0)=\frac12=-\varphi'(1)`$、$`\varphi''=-1`$ である。2 回部分積分して $`\int_0^1\varphi g''=-\int_0^1\varphi'g'=-[\varphi'g]_0^1+\int_0^1\varphi''g=\frac{g(0)+g(1)}2-\int_0^1g`$。$`g''(s)=h^2f'''(a+sh)`$ なのでこれが恒等式であり、$`\int_0^1\varphi=\frac1{12}`$ から評価が従う。3 次多項式では $`f'''`$ は定数である。$`\square`$

<a id="r-1-13"></a>

**命題 1.13（商か中点値か）.** 単位丸め誤差 $`u\le1/20`$ の浮動小数点演算を仮定する。すなわち浮動小数点数 $`x,y`$ と $`\circ\in\{+,-,\times,/\}`$ に対し、計算結果は $`(x\circ y)(1+\theta)`$、$`|\theta|\le u`$ である（オーバーフローもアンダーフローもないとする）。$`a\ne b`$ を実の浮動小数点数、$`h=|a-b|`$、$`x=f^{[1]}(a,b)`$ とし、$`\hat f_a,\hat f_b`$ を $`|\hat f_a-f(a)|\le\delta_a`$、$`|\hat f_b-f(b)|\le\delta_b`$ を満たす計算値とする。

(a) $`\hat f_a-\hat f_b`$ と $`a-b`$ の計算された商 $`\hat q`$ は $`|\hat q-x|\le R:=(\delta_a+\delta_b)/h+5u|\hat q|`$ を満たす。

(b) $`|\hat m-x|\le T`$ とする（たとえば $`\hat m`$ を計算された中点値、$`T`$ を命題 1.12 の評価にその丸め誤差を加えたもの）。$`|\hat q-\hat m|\le2R`$ ならば $`y=\hat m`$、そうでなければ $`y=\hat q`$ とおく。このとき $`|y-x|\le\min(T,3R)`$。

(c) $`\delta_a+\delta_b\le\Delta`$ かつ $`T=Mh^2+T_0`$ ならば、任意の $`h`$ に対して $`|y-x|\le(3\Delta)^{2/3}M^{1/3}+T_0+15u|\hat q|`$。二つの誤差項 $`3\Delta/h`$ と $`Mh^2`$ は $`h_*=(3\Delta/M)^{1/3}`$ で釣り合う。

*証明.* (a) $`\tilde q=(\hat f_a-\hat f_b)/(a-b)`$ とおく。三つの演算から $`\hat q=\tilde q\,(1+\theta_1)(1+\theta_3)/(1+\theta_2)=\tilde q\,(1+\theta)`$、$`|\theta_i|\le u`$ となり、$`-3u\le\frac{(1-u)^2}{1+u}-1\le\theta\le\frac{(1+u)^2}{1-u}-1=\frac{3u+u^2}{1-u}\le4u`$。よって $`|\hat q-\tilde q|=|\theta|\,|\tilde q|\le\frac{4u}{1-4u}|\hat q|\le5u|\hat q|`$ であり、$`|\tilde q-x|=|(\hat f_a-f(a))-(\hat f_b-f(b))|/h\le(\delta_a+\delta_b)/h`$。(b) $`y=\hat m`$ ならば $`|y-x|\le T`$ かつ $`|y-x|\le|\hat m-\hat q|+|\hat q-x|\le3R`$。$`y=\hat q`$ ならば $`|y-x|\le R`$ であり、$`T\ge|\hat m-x|\ge|\hat m-\hat q|-|\hat q-x|>2R-R=R`$。(c) $`\min(T,3R)\le\min(Mh^2,3\Delta/h)+T_0+15u|\hat q|`$。$`h\ge h_*`$ では $`3\Delta/h\le3\Delta/h_*`$、$`h<h_*`$ では $`Mh^2<Mh_*^2`$ であり、$`3\Delta/h_*=Mh_*^2=(3\Delta)^{2/3}M^{1/3}`$。$`\square`$

この証明は複素数値の $`f`$ に対してもそのまま成り立つ（実部と虚部は別々に丸められ、$`a-b`$ は実数である。このとき $`\theta`$ は同じ評価をもつ複素数である）。規則 (b) はデータから計算できる $`R`$ を必要とするが、$`T`$ は必要としない。$`f'''`$ について何も知らなくても、その誤差は二つの評価の小さいほうの高々 3 倍である。固定した閾値（$`h\le h_*`$ のとき、かつそのときに限り中点値をとる）には $`M=\max|f'''|/12`$ が必要になる。ライブラリの以前の版の閾値 $`h\le\varepsilon^{1/3}\max(|a|,|b|)`$（$`\varepsilon=2u`$）はこの種のもので、$`M`$ を $`\Delta/\max(|a|,|b|)^3`$ 程度とみなすことにあたる。これは導関数が $`|x|`$ とともにスケールする $`x^p`$ や $`\log x`$ には正しいが、$`|x|`$ と無関係なスケールで変化する関数には誤りである。$`a\approx c=10^4`$ における $`f(x)=e^{x-c}`$ では、$`h\approx0.06`$ まで中点値をとり、相対誤差は $`h^2/12\approx3\cdot10^{-4}`$ に達した（そこでは命題 1.12 は因子 $`e^{h}`$ を除いて等式である）。$`e^{-i\tau x}`$ と大きな定数シフトをもつハミルトニアンでも同じことが起こる。

*ライブラリでは：* `functional.divided_differences` は固有値のすべての対に規則 (b) を適用する。その際 $`\delta_a=\varepsilon\,(|\hat f_a|+|a\hat f'_a|)`$ とし、$`5u|\hat q|`$ の代わりに $`4\varepsilon|\hat q|`$ を用いる。$`\delta_a`$ の選び方は $`f`$ の評価に関する仮定であって定理ではない。すなわち、高々 $`\varepsilon|a|`$ だけ摂動された引数において相対誤差が高々 $`\varepsilon`$ であるという仮定である（一次の近似。第二項は $`\tau x`$ を丸める $`e^{-i\tau x}`$ のような関数をカバーする）。同じ規則は正規行列の複素固有値にも適用されるが、これは証明の範囲外である。等しい固有値の対、および固有値ソルバーが分離できない対 $`h\le n\varepsilon\max_l|\lambda_l|`$ に対しても中点値をとる。そこでは命題 1.12 によりその誤差は高々 $`(n\varepsilon\max_l|\lambda_l|)^2\max|f'''|/12`$ であり、一方、指示関数のような不連続な $`f`$ に対しては商は意味をもたない。

*計算：* 縮退したスペクトルにおける定理 1.11 と差分商との比較、および $`c=10^4`$ の $`e^{x-c}`$ に対する命題 1.12 と 1.13 の検証 — [チュートリアル、1.5節](../../tutorial/ja/01_matrix_algebras.md#s-1-5)。

<a id="s-1-4"></a>

## 1.4 極分解と特異値

任意の $`A`$ に対し、行列 $`A^*A`$ は正である（系 1.7）。その平方根
$`|A|=(A^*A)^{1/2}`$ を $`A`$ の **絶対値** という。

<a id="r-1-8"></a>

**定理 1.8（極分解）.** 任意の $`A\in M_n(\mathbb C)`$ は、ユニタリ $`W`$ を用いて
$`A=W|A|`$ と書ける。

*証明.* 任意の $`x`$ に対し $`\||A|x\|^2=\langle x,|A|^2x\rangle=\langle x,A^*Ax\rangle=\|Ax\|^2`$ である。
したがって $`W_0:|A|x\mapsto Ax`$ は $`|A|`$ の値域
$`R_1`$ から $`A`$ の値域 $`R_2`$ の上への well-defined な線形等長写像である：$`|A|x=|A|x'`$ ならば
$`\|A(x-x')\|=\||A|(x-x')\|=0`$ である。特に $`\dim R_1=\dim R_2`$ であるから、直交補空間の次元も等しい。$`W_1:R_1^\perp\to R_2^\perp`$ をそれらの間の任意の等長写像とする（正規直交基底を正規直交基底に写せばよい）。
$`W=W_0\oplus W_1`$、すなわち $`u_1\in R_1`$, $`u_2\in R_1^\perp`$ に対し $`W(u_1+u_2)=W_0u_1+W_1u_2`$ とおく。
$`W_0u_1\in R_2\perp R_2^\perp\ni W_1u_2`$ であるから、ピタゴラスの定理によりすべての $`u`$ に対し $`\|Wu\|^2=\|u_1\|^2+\|u_2\|^2=\|u\|^2`$ である。
よってすべての $`u`$ に対し $`\langle u,(W^*W-1)u\rangle=0`$ であり、偏極（系 1.4(3) の証明）により $`W^*W=1`$ である。
すると $`W`$ は可逆で $`W^{-1}=W^*`$ であるから $`WW^*=1`$ でもある：$`W`$ はユニタリであり、
$`W|A|x=W_0|A|x=Ax`$ である。$`\square`$

*計算：* ジニブル行列に対する定理 1.8 の検証（$`|A|`$ は `abs` で計算する） — [チュートリアル、1.6節](../../tutorial/ja/01_matrix_algebras.md#s-1-6)。

$`|A|`$ の固有値 $`s_1\ge\dots\ge s_n\ge0`$ を $`A`$ の **特異値** という。本講義を通じて三つのノルムを用いる。
```math
\|A\|=s_1,\qquad \|A\|_2=\Big(\sum_is_i^2\Big)^{1/2}=(\mathrm{Tr}\,A^*A)^{1/2},
\qquad\|A\|_1=\sum_is_i=\mathrm{Tr}|A|.
```
第一の等式は、任意の $`x`$ に対し $`\|Ax\|=\||A|x\|`$ であるから系 1.4(1) により
$`\|A\|=\||A|\|=s_1`$ となることから成り立つ。第二の等式は
$`\mathrm{Tr}\,A^*A=\mathrm{Tr}|A|^2=\sum_is_i^2`$ による。

<a id="r-1-9"></a>

**定理 1.9（トレースノルムの双対性）.**
$`\|A\|_1=\max_U|\mathrm{Tr}(AU)|`$（$`U`$ はユニタリ行列全体を動く）であり、最大値は
$`A=W|A|`$ として $`U=W^*`$ で達成される。特に $`|\mathrm{Tr}\,A|\le\|A\|_1`$ である。

*証明.* $`\mathrm{Tr}(AW^*)=\mathrm{Tr}(W|A|W^*)=\mathrm{Tr}|A|`$ である。任意のユニタリ $`U`$ に対し、
トレースの巡回性と $`\langle X,Y\rangle_{HS}=\mathrm{Tr}(X^*Y)`$ に関するコーシー・シュワルツの不等式により
```math
|\mathrm{Tr}(AU)|=|\mathrm{Tr}(|A|^{1/2}\,UW\,|A|^{1/2})|
=|\langle|A|^{1/2},\,UW|A|^{1/2}\rangle_{HS}|
\le\||A|^{1/2}\|_2\,\|UW|A|^{1/2}\|_2=\mathrm{Tr}|A|,
```
となる。ここでユニタリ $`UW`$ に対し $`\|UW X\|_2=\|X\|_2`$ であること、および $`\||A|^{1/2}\|_2^2=\mathrm{Tr}|A|`$ を用いた。$`\square`$

特に、上の三つの量はノルムである。$`\|\cdot\|`$ は作用素ノルムであり、
$`\|\cdot\|_2`$ は内積 $`\langle X,Y\rangle_{HS}`$ から定まるノルムである。
$`\|\cdot\|_1`$ については、定理 1.9 により
$`\|A+B\|_1=\max_U|\mathrm{Tr}(AU)+\mathrm{Tr}(BU)|\le\|A\|_1+\|B\|_1`$ および
$`\|cA\|_1=\max_U|c\mathrm{Tr}(AU)|=|c|\,\|A\|_1`$ が成り立ち、また
$`\|A\|_1=0`$ ならば $`|A|=0`$、したがって $`A=W|A|=0`$ である。$`P\ge0`$ に対しては定理 1.6 により $`|P|=(P^2)^{1/2}=P`$ であるから、
$`\|P\|_1=\mathrm{Tr}\,P`$ である。

*計算：* 三つのノルム（`operator_norm`、`frobenius_norm`、`trace_norm`）と定理 1.9 — [チュートリアル、1.6節](../../tutorial/ja/01_matrix_algebras.md#s-1-6)。

<a id="s-1-5"></a>

## 1.5 ライブラリにおける作用素

*ライブラリでは：* `TypeIAlgebra([n], [n])` としての環 $`M_n(\mathbb C)`$、その作用素の格納方法、および関数計算 `apply_function` — [チュートリアル、1.5節](../../tutorial/ja/01_matrix_algebras.md#s-1-5)。

<a id="s-1-6"></a>

## 1.6 ノルムとトレースノルムの双対性

*計算：* ジニブル行列に対する 1.4節のノルム、極分解、双対性の検証 — [チュートリアル、1.6節](../../tutorial/ja/01_matrix_algebras.md#s-1-6)。

<a id="s-1-7"></a>

## 1.7 バッチ：多数の作用素を一度に

多数のランダムな作用素を一度に用いて推定する量の例として、**ガウス型ユニタリアンサンブル**（GUE）の厳密なモーメントを計算する。

$`G`$ を独立な成分をもつ行列で、各成分の実部と虚部は独立に
$`\sim N(0,\tfrac12)`$ に従うものとする。したがって $`\mathbb E|G_{ij}|^2=1`$ である。$`H=(G+G^*)/\sqrt2`$ とおく。

<a id="r-1-10"></a>

**命題 1.10.** $`\mathbb E\,\mathrm{Tr}\,H^4=2n^3+n`$ である。したがって
$`X=H/\sqrt n`$ に対し $`\ \mathbb E\,\tfrac1n\mathrm{Tr}\,X^4=2+\tfrac1{n^2}`$ である。

*証明.* *イッサーリス（ウィック）の公式*（定理 A.24、$`m=4`$ の場合）。$`G`$ の成分の実部と虚部は独立な中心化ガウス変数であるから、
一つの中心化ガウスベクトルの座標であり、$`H`$ と $`\bar H`$ の成分はそれらの複素線形結合である。
したがって定理 A.24 により、そのような任意の四つの変数 $`Z_1,\dots,Z_4`$（重複を許す）に対し
```math
\mathbb E[Z_1Z_2Z_3Z_4]=\mathbb E[Z_1Z_2]\,\mathbb E[Z_3Z_4]+\mathbb E[Z_1Z_3]\,\mathbb E[Z_2Z_4]+\mathbb E[Z_1Z_4]\,\mathbb E[Z_2Z_3],
```
すなわち $`\{1,2,3,4\}`$ の三通りの対合わせ（ペアリング）にわたる、対ごとの期待値の積の和（複素共役はとらない）が成り立つ。
付録 A の証明はモーメント母関数のべき級数をガウス分布のモーメントと比較するものであり、期待値の中での微分は行わない。

*モーメントの計算.* 対ごとの期待値は
$`\mathbb E[H_{ij}H_{kl}]=\tfrac12\mathbb E[(G_{ij}+\bar G_{ji})(G_{kl}+\bar G_{lk})]=\delta_{il}\delta_{jk}`$
である。実際、$`\mathbb E[G_{ij}G_{kl}]=0`$ かつ $`\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}`$ である（命題 A.16(a)。
$`(i,j)\neq(k,l)`$ ならば二つの成分は独立かつ中心化されており、$`G_{ij}=x+iy`$ に対しては
$`\mathbb E[G_{ij}^2]=\mathbb E[x^2-y^2]+2i\,\mathbb E[xy]=0`$ および $`\mathbb E|G_{ij}|^2=\tfrac12+\tfrac12=1`$ である）。
$`\mathrm{Tr}\,H^4=\sum_{ijkl}H_{ij}H_{jk}H_{kl}H_{li}`$ において、各対合わせは次を与える：
$`\{12\}\{34\}`$：$`\delta_{ik}\delta_{ki}`$、和は $`n^3`$（自由な添字は $`i=k`$, $`j`$, $`l`$）；
$`\{14\}\{23\}`$：$`\delta_{jl}\delta_{jl}`$、和は $`n^3`$；
$`\{13\}\{24\}`$：$`\delta_{il}\delta_{jk}\delta_{ji}\delta_{kl}`$、これは $`i=j=k=l`$ を強制し、和は $`n`$。
$`n\cdot n^2`$ で割ると第二の主張を得る。$`\square`$

二つの *非交差* な対合わせ $`\{12\}\{34\}`$ と $`\{14\}\{23\}`$（対を点 $`1,2,3,4`$ の上に互いに交わらない弧として描ける）は、すべての $`n`$ に対して主要項 $`2`$ を与え、交差する対合わせ $`\{13\}\{24\}`$ は $`1/n^2`$ の補正を与える。

*計算：* $`200\,000`$ 個の標本のバッチを用いたモンテカルロ法による命題 1.10 の検証 — [チュートリアル、1.7節](../../tutorial/ja/01_matrix_algebras.md#s-1-7)。

<a id="s-1-8"></a>

## 1.8 比較の対象

*ライブラリでは：* 比較のための参照実装 — [チュートリアル、1.8節](../../tutorial/ja/01_matrix_algebras.md#s-1-8)。

<a id="x-1"></a>

## 演習問題

1. 正規な $`A`$ と $`\mathrm{spec}\,A`$ 上の関数 $`f,g`$ に対し、$`(fg)(A)=f(A)g(A)`$ を示せ。
2. 定理 1.9 を用いて、$`\|A\|_1=\|A^*\|_1`$ および $`\|AB\|_1\le\|A\|\,\|B\|_1`$ を示せ。
3. $`A,B>0`$ の **幾何平均** を $`A\#B=A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}`$ と定める。
   $`A\#B`$ が $`XA^{-1}X=B`$ のただ一つの正の解 $`X`$ であることを証明し、そこから
   $`A\#B=B\#A`$ を導け。

<details>
<summary>解答</summary>

1. 定理 1.5(c) による。

2. $`|\mathrm{Tr}(A^*U)|=|\overline{\mathrm{Tr}(U^*A)}|=|\mathrm{Tr}(AU^*)|`$ であり、
   $`U\mapsto U^*`$ はユニタリ群の全単射であるから、定理 1.9 における最大値は一致する。第二の主張について、$`B=V|B|`$ とする。ユニタリ $`U`$ に対し
   $`|\mathrm{Tr}(ABU)|=|\mathrm{Tr}(|B|^{1/2}UAV|B|^{1/2})|    \le\||B|^{1/2}\|_2\,\|UAV|B|^{1/2}\|_2\le\|A\|\mathrm{Tr}|B|`$ である。ここで
   $`\|YX\|_2\le\|Y\|\,\|X\|_2`$（$`X`$ の各列に $`\|Yx\|\le\|Y\|\,\|x\|`$ を適用する）および
   ユニタリ $`U,V`$ に対する $`\|UAV\|=\|A\|`$（ユニタリはベクトルの長さを保つ）を用いた。

3. $`T=A^{-1/2}XA^{-1/2}`$ とおく。すると $`XA^{-1}X=B`$ $`\iff`$ $`T^2=A^{-1/2}BA^{-1/2}`$ であり、
   $`\langle x,Xx\rangle=\langle A^{1/2}x,TA^{1/2}x\rangle`$ かつ
   $`A^{1/2}`$ は可逆である（ここで $`A^{\pm1/2}`$ は $`\mathrm{spec}\,A\subset(0,\infty)`$ 上の $`t^{\pm1/2}`$ による関数計算で与えられ、
   定理 1.5(c) により $`A^{-1/2}A^{1/2}=1`$ である）から、$`X\ge0`$ $`\iff`$ $`T\ge0`$ である。また
   $`\langle x,A^{-1/2}BA^{-1/2}x\rangle=\langle A^{-1/2}x,BA^{-1/2}x\rangle\ge0`$ であるから $`A^{-1/2}BA^{-1/2}\ge0`$ である。
   定理 1.6 により $`T=(A^{-1/2}BA^{-1/2})^{1/2}`$ はただ一つの正の解であり、すなわち $`X=A\#B`$ である。$`XA^{-1}X=B`$ の任意の解 $`X`$ は可逆である。
   実際 $`\det(X)^2\det(A)^{-1}=\det B\neq0`$ である。よって $`XA^{-1}X=B`$ は
   $`XB^{-1}X=A`$ と同値であり（$`A^{-1}=X^{-1}BX^{-1}`$ の両辺の逆をとる）、一意性から
   $`A\#B=B\#A`$ を得る。

</details>

*計算：* 演習問題 1 と 3 の数値的な検証 — [チュートリアル、演習問題](../../tutorial/ja/01_matrix_algebras.md#x-1)。

**次へ：**[第2章](02_von_neumann_algebras.md) では $`M_n`$ を行列環の直和に置き換え、すべての有限次元フォン・ノイマン環がこの形であることを証明する。
