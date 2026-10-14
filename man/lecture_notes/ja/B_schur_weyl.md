[English](../B_schur_weyl.md) · [Español](../es/B_schur_weyl.md) · [Français](../fr/B_schur_weyl.md) · [Deutsch](../de/B_schur_weyl.md) · [中文](../zh/B_schur_weyl.md) · **日本語** · [Русский](../ru/B_schur_weyl.md)

# 付録 B. シューア・ワイル双対性とワインガルテン計算

ユニタリ群 $`U(d)`$ は $`V\mapsto V^{\otimes n}`$ によって $`(\mathbb C^d)^{\otimes n}`$ に作用し、対称群 $`S_n`$ はテンソル因子の
置換によって作用する。これら二つの作用素の族は、それぞれが他方の可換子環である（**シューア・ワイル双対性**）。
この付録ではこの双対性を完全に証明し、そこから任意の次数のハール・ユニタリのモーメントに対するワインガルテン公式を
導く。その途中で、第4、5、9、11章がここから引用する事実を証明する：すべての $`L`$ に対する定理 4.14(b)、全スピンの
重複度空間上の $`S_L`$ の表現の既約性、4.5節のカタラン数、定理 5.12、9.3節で用いる $`SU(2)`$ に対するウィグナー・
エッカルトの定理、そして第11章で用いる直交群のハール測度である。

道筋は初等的である：偏極恒等式により置換の可換子環をすべての $`A^{\otimes n}`$（$`A\in M_d(\mathbb C)`$）の張る空間と
同一視し、ケーリー変換の議論によりユニタリ $`V^{\otimes n}`$ が同じ空間を張ることを示し、二重可換子環定理（定理 2.2）に
よってこれを双対性に変える。ワインガルテン公式は、ハール平均が可換子環への直交射影であることから従う。

この付録の数値的検証は[チュートリアル、付録 B](../../tutorial/ja/B_schur_weyl.md)にある。

以下を通じて $`d,n\ge1`$ とし、$`\mathcal H_n=(\mathbb C^d)^{\otimes n}`$ に正規直交基底
$`e_{\mathbf i}=e_{i_1}\otimes\dots\otimes e_{i_n}`$ を与える。ここで多重指数 $`\mathbf i=(i_1,\dots,i_n)\in[d]^n`$ は写像
$`[n]\to[d]`$、$`m\mapsto i_m`$ とみなす（$`[d]=\{1,\dots,d\}`$）。$`\mathrm{End}\,\mathcal H_n=M_{d^n}(\mathbb C)`$ には
ヒルベルト・シュミット内積 $`\langle X,Y\rangle=\mathrm{Tr}(X^*Y)`$ を入れる。$`X\in M_d`$ と $`1\le m\le n`$ に対し
$`X^{(m)}=1\otimes\dots\otimes X\otimes\dots\otimes1`$（$`X`$ は $`m`$ 番目の位置）、$`d\Pi(X)=\sum_{m=1}^nX^{(m)}`$ とおく。また
```math
\mathcal T_n=\mathrm{span}\{A^{\otimes n}:A\in M_d(\mathbb C)\},\qquad
\mathcal U_n=\mathrm{span}\{V^{\otimes n}:V\in U(d)\}.
```
と書く。クロネッカー積の積は因子ごとに計算される：$`(A_1\otimes\dots\otimes A_n)(B_1\otimes\dots\otimes B_n)=A_1B_1\otimes\dots\otimes A_nB_n`$、
また $`(A_1\otimes\dots\otimes A_n)^*=A_1^*\otimes\dots\otimes A_n^*`$ である。置換は写像として合成し、
$`(\sigma\tau)(m)=\sigma(\tau(m))`$ とする。$`\#\text{cycles}(\sigma)`$ は $`\sigma`$ の巡回の個数を、不動点も含めて数える。

<a id="s-B-1"></a>

## B.1 テンソル因子の置換

$`\sigma\in S_n`$ に対し $`\pi(\sigma)e_{\mathbf i}=e_{\mathbf i\circ\sigma^{-1}}`$ とおく。すなわち、スロット $`m`$ の基底ベクトルはスロット
$`\sigma(m)`$ に移る。これは `decompose.permutation_representation` および第2章と第4章の規約である。

*ライブラリでは：* 作用素 $`\pi(\sigma)`$ は `decompose.permutation_representation` が構成する — [チュートリアル、B.3節](../../tutorial/ja/B_schur_weyl.md#s-B-3)。

<a id="r-B-1"></a>

**補題 B.1.** (a) $`\pi`$ は $`S_n`$ のユニタリ表現であり、すべての $`x_m\in\mathbb C^d`$ に対し
$`\pi(\sigma)(x_1\otimes\dots\otimes x_n)=x_{\sigma^{-1}(1)}\otimes\dots\otimes x_{\sigma^{-1}(n)}`$ である。
(b) すべての $`A_m\in M_d`$ に対し $`\pi(\sigma)(A_1\otimes\dots\otimes A_n)\pi(\sigma)^*=A_{\sigma^{-1}(1)}\otimes\dots\otimes A_{\sigma^{-1}(n)}`$。
特に $`\pi(\sigma)`$ はすべての $`A^{\otimes n}`$ およびすべての $`d\Pi(X)`$ と可換であり、
$`\pi(\sigma)X^{(m)}\pi(\sigma)^*=X^{(\sigma(m))}`$ である。
(c) $`\mathrm{Tr}\,\pi(\sigma)=d^{\#\text{cycles}(\sigma)}`$ かつ $`\langle\pi(\sigma),\pi(\tau)\rangle=d^{\#\text{cycles}(\sigma^{-1}\tau)}`$。
巡回の個数は共役と逆元をとる操作で不変であり、二つの置換が共役であるのは、それらが同じ巡回型をもつときかつそのときに限る。
(d) $`n!`$ 個の作用素 $`\pi(\sigma)`$ が一次独立であるのは $`n\le d`$ のときかつそのときに限る。

*証明.* (a) $`\mathbf i\mapsto\mathbf i\circ\sigma^{-1}`$ は $`[d]^n`$ の全単射であるから、$`\pi(\sigma)`$ は正規直交基底を置換し、
ユニタリである。$`\pi(\sigma)\pi(\tau)e_{\mathbf i}=\pi(\sigma)e_{\mathbf i\circ\tau^{-1}}=e_{\mathbf i\circ\tau^{-1}\circ\sigma^{-1}}=\pi(\sigma\tau)e_{\mathbf i}`$。
$`x_m=\sum_ax_m(a)e_a`$ と展開すると $`x_1\otimes\dots\otimes x_n=\sum_{\mathbf i}\prod_mx_m(i_m)\,e_{\mathbf i}`$ である。$`\pi(\sigma)`$ を作用させ
$`\mathbf j=\mathbf i\circ\sigma^{-1}`$ と置き換えると $`\sum_{\mathbf j}\prod_mx_m(j_{\sigma(m)})\,e_{\mathbf j}=\sum_{\mathbf j}\prod_{p}x_{\sigma^{-1}(p)}(j_p)\,e_{\mathbf j}`$
となり、これは主張された積ベクトルである。
(b) (a) により、積ベクトル上で：$`\pi(\sigma)^*=\pi(\sigma^{-1})`$ は $`x_{\sigma(m)}`$ をスロット $`m`$ に置き、作用素はそこで $`A_m`$ を
作用させ、$`\pi(\sigma)`$ はスロット $`m`$ をスロット $`\sigma(m)`$ に移す。したがってスロット $`p`$ には最終的に $`A_{\sigma^{-1}(p)}x_p`$ が入る。
積ベクトルは $`\mathcal H_n`$ を張る。$`X^{(m)}`$ については因子 $`X`$ がスロット $`m`$ からスロット $`\sigma(m)`$ に移るので、$`d\Pi(X)=\sum_mX^{(m)}`$ は不変である。
(c) $`\pi(\sigma)`$ は置換行列であり、$`e_{\mathbf i}`$ が固定されるのはすべての $`p`$ に対し $`i_{\sigma^{-1}(p)}=i_p`$ となるとき、すなわち
$`\mathbf i`$ が $`\sigma`$ の各巡回上で一定であるときかつそのときに限る。そのような $`\mathbf i`$ は $`d^{\#\text{cycles}(\sigma)}`$ 個ある。さらに
$`\langle\pi(\sigma),\pi(\tau)\rangle=\mathrm{Tr}(\pi(\sigma^{-1})\pi(\tau))=\mathrm{Tr}\,\pi(\sigma^{-1}\tau)`$ である。
$`\kappa\rho\kappa^{-1}`$ の巡回は $`\rho`$ の巡回の $`\kappa`$ による像であり、$`\rho^{-1}`$ は同じ巡回を逆向きにたどる。逆に
$`\rho,\rho'`$ が同じ巡回型をもつならば、それらの巡回を長さの順に並べて上下二行に書き、$`\kappa`$ を第一行の各要素をその
真下の要素に写す写像とする。このとき $`\kappa\rho\kappa^{-1}=\rho'`$ である。
(d) $`n\le d`$ ならば $`\mathbf u=(1,2,\dots,n)`$ とおく。ベクトル $`\pi(\sigma)e_{\mathbf u}=e_{\mathbf u\circ\sigma^{-1}}`$ は相異なる $`\sigma`$ に対し
相異なる基底ベクトルであるから、$`\sum_\sigma c_\sigma\pi(\sigma)=0`$ を $`e_{\mathbf u}`$ に作用させるとすべての $`c_\sigma=0`$ を得る。$`n>d`$ ならば
$`\mathcal A=\sum_\sigma\mathrm{sgn}(\sigma)\pi(\sigma)`$ とおく。任意の $`\mathbf i\in[d]^n`$ はある $`a\neq b`$ に対し $`i_a=i_b`$ を満たす。
互換 $`t=(a\,b)`$ に対し $`\pi(t)e_{\mathbf i}=e_{\mathbf i}`$ かつ $`\mathcal A\pi(t)=\sum_\sigma\mathrm{sgn}(\sigma)\pi(\sigma t)=-\mathcal A`$ である。
よってすべての $`\mathbf i`$ に対し $`\mathcal Ae_{\mathbf i}=\mathcal A\pi(t)e_{\mathbf i}=-\mathcal Ae_{\mathbf i}=0`$ であり、$`\mathcal A=0`$ は自明でない
一次関係である。$`\square`$

*計算：* 補題 B.1(a) と (d) — [チュートリアル、B.3節](../../tutorial/ja/B_schur_weyl.md#s-B-3)。

<a id="s-B-2"></a>

## B.2 偏極：置換の可換子環

<a id="r-B-2"></a>

**補題 B.2（偏極）.** $`A_1,\dots,A_n\in M_d(\mathbb C)`$ に対し、$`|\varepsilon|=\sum_m\varepsilon_m`$ として
```math
\sum_{\sigma\in S_n}A_{\sigma(1)}\otimes\dots\otimes A_{\sigma(n)}=
\sum_{\varepsilon\in\{0,1\}^n}(-1)^{n-|\varepsilon|}\Big(\sum_{m=1}^n\varepsilon_mA_m\Big)^{\otimes n}.
```
が成り立つ。

*証明.* 多重線形性により $`\big(\sum_m\varepsilon_mA_m\big)^{\otimes n}=\sum_f\varepsilon_{f(1)}\cdots\varepsilon_{f(n)}\,A_{f(1)}\otimes\dots\otimes A_{f(n)}`$
であり、和はすべての写像 $`f:[n]\to[n]`$ にわたる。$`r\ge1`$ に対し $`\varepsilon_m^r=\varepsilon_m`$ であるから、右辺における項 $`f`$ の係数は
```math
\sum_{\varepsilon}\prod_{m=1}^n(-1)^{1-\varepsilon_m}\prod_{m\in f([n])}\varepsilon_m
=\prod_{m\notin f([n])}\Big(\sum_{\varepsilon_m=0,1}(-1)^{1-\varepsilon_m}\Big)\prod_{m\in f([n])}\Big(\sum_{\varepsilon_m=0,1}(-1)^{1-\varepsilon_m}\varepsilon_m\Big).
```
である。第一の種類の因子は $`-1+1=0`$、第二の種類は $`1`$ である。したがって係数は $`f`$ が全射、すなわち置換ならば $`1`$、
そうでなければ $`0`$ である。$`\square`$

*計算：* $`n=3`$ に対する補題 B.2 — [チュートリアル、B.3節](../../tutorial/ja/B_schur_weyl.md#s-B-3)。

$`\mathrm{End}\,\mathcal H_n`$ 上の**対称化作用素**を $`\mathcal S(X)=\frac1{n!}\sum_\sigma\pi(\sigma)X\pi(\sigma)^*`$ とする。

<a id="r-B-3"></a>

**命題 B.3.** $`\mathcal S`$ は $`\mathrm{End}\,\mathcal H_n`$ から $`\pi(S_n)'`$ への直交射影であり、
```math
\pi(S_n)'=\mathcal T_n,\qquad\dim\pi(S_n)'=\binom{d^2+n-1}n .
```
である。

*証明.* $`\mathcal S`$ は表現 $`\pi`$ における有限群 $`S_n`$ のツワールであるから、$`\pi(S_n)'`$ への直交射影である
（$`\pi_1=\pi_2=\pi`$ とした定理 12.5；4.2節も参照）。補題 B.1(b) により $`\mathcal T_n\subseteq\pi(S_n)'`$ である。逆に
$`\mathrm{End}\,\mathcal H_n`$ は行列単位の積 $`A_1\otimes\dots\otimes A_n`$ で張られ、補題 B.1(b) と置き換え $`\sigma\to\sigma^{-1}`$ により
$`\mathcal S(A_1\otimes\dots\otimes A_n)=\frac1{n!}\sum_\sigma A_{\sigma(1)}\otimes\dots\otimes A_{\sigma(n)}`$ であり、これは補題 B.2 により
$`\mathcal T_n`$ に属する。よって $`\pi(S_n)'=\mathcal S(\mathrm{End}\,\mathcal H_n)\subseteq\mathcal T_n`$ である。

*次元.* $`P=[d]\times[d]`$ とし、$`\alpha\in P^n`$ に対し $`E_\alpha=E_{\alpha_1}\otimes\dots\otimes E_{\alpha_n}`$ とおく。ここで行列単位は
$`E_{(a,b)}=|e_a\rangle\langle e_b|`$ であり、これらは $`\mathrm{End}\,\mathcal H_n`$ の基底をなす。補題 B.1(b) により
$`\mathcal S(E_\alpha)=\frac1{n!}\sum_\sigma E_{\alpha\circ\sigma^{-1}}`$ であり、これは軌道
$`O_\alpha=\{\alpha\circ\sigma\}`$ に対して $`|O_\alpha|^{-1}\sum_{\beta\in O_\alpha}E_\beta`$ に等しい（軌道の各要素は同じ個数 $`n!/|O_\alpha|`$ の
$`\sigma`$ に対して現れる）。異なる軌道にわたる和は互いに素な基底要素の集合を含むから、$`\mathcal S`$ の像は軌道で添字づけられた
基底をもつ。二つの組が同じ軌道に属するのは、$`P`$ の各要素が両者に同じ回数だけ現れるときかつそのときに限る：軌道は
$`P`$ の $`N=d^2`$ 個の要素からなる大きさ $`n`$ の多重集合である。$`P`$ を $`1,\dots,N`$ と順序づけると、そのような多重集合は列
$`1\le c_1\le\dots\le c_n\le N`$ であり、$`c_m\mapsto c_m+m-1`$ は $`[N+n-1]`$ の狭義単調増加列、すなわち $`n`$ 元部分集合への
全単射である。その個数は $`\binom{N+n-1}n`$ である。$`\square`$

*計算：* 小さな $`d`$ と $`n`$ に対する $`\dim\pi(S_n)'=\binom{d^2+n-1}n`$ — [チュートリアル、B.3節](../../tutorial/ja/B_schur_weyl.md#s-B-3)。

<a id="s-B-3"></a>

## B.3 シューア・ワイル双対性

<a id="r-B-4"></a>

**補題 B.4（ユニタリ群上で消える多項式）.** $`p:M_d(\mathbb C)\to\mathbb C`$ を $`d^2`$ 個の成分の多項式で、
次数 $`n\ge1`$ の斉次なものとする。すべての $`V\in U(d)`$ に対し $`p(V)=0`$ ならば $`p=0`$ である。

*証明.* *ステップ1（ケーリー変換）.* $`A\in M_d`$ に対し $`q(A)=p\big((1+iA)\mathrm{adj}(1-iA)\big)`$ とおく。ここで
$`\mathrm{adj}`$ は余因子行列である。$`q`$ は $`A`$ の成分の多項式である。$`H`$ をエルミートとする。その固有値は実数である
（系 1.4）から $`1-iH`$ は可逆であり、$`\mathrm{adj}(1-iH)=\det(1-iH)(1-iH)^{-1}`$（クラメルの公式）、そして斉次性により
$`C=(1+iH)(1-iH)^{-1}`$ として $`q(H)=\det(1-iH)^n\,p(C)`$ である。行列 $`1\pm iH`$ とその逆行列は $`H`$ の関数であって互いに
可換であり（定理 1.5）、$`(1+iH)^*=1-iH`$ であるから
$`C^*C=(1+iH)^{-1}(1-iH)(1+iH)(1-iH)^{-1}=1`$ である。したがって $`q(H)=0`$ である。

*ステップ2（エルミート行列からすべての行列へ）.* エルミート行列のなす実ベクトル空間の基底 $`B_1,\dots,B_{d^2}`$ を選ぶ。
これは $`\mathbb C`$ 上の $`M_d(\mathbb C)`$ の基底でもある。実際、任意の $`A`$ は $`\frac{A+A^*}2+i\frac{A-A^*}{2i}`$ であり、次元が一致する。
$`r(t)=q(\sum_kt_kB_k)`$ は $`t\in\mathbb C^{d^2}`$ の多項式で、$`\mathbb R^{d^2}`$ 上で消える。このような多項式の係数はすべて $`0`$ である：
一変数ならば無限個の根をもつ。$`N`$ 変数ならば $`r=\sum_jr_j(t')t_N^j`$ と書く——各実数の $`t'`$ に対し $`t_N`$ の多項式は
$`\mathbb R`$ 上で消えるので、各 $`r_j`$ は $`\mathbb R^{N-1}`$ 上で消え、帰納法が適用できる。したがってすべての
$`A=\sum_kt_kB_k\in M_d`$ に対し $`q(A)=r(t)=0`$ である。

*ステップ3（変換の逆）.* $`B\in M_d`$ を $`1+B`$ が可逆なものとし、$`A=-i(1+B)^{-1}(B-1)`$ とおく。すると
$`1-iA=1-(1+B)^{-1}(B-1)=2(1+B)^{-1}`$ かつ $`1+iA=(1+B)^{-1}2B`$ である。$`B`$ は $`1+B`$ と可換であるから
$`(1+iA)(1-iA)^{-1}=(1+B)^{-1}B(1+B)=B`$ である。よって $`\det(1-iA)\neq0`$ のもとで $`0=q(A)=\det(1-iA)^np(B)`$ となり、$`p(B)=0`$ である。

*ステップ4.* 任意の $`B\neq0`$ と $`0<t<1/\|B\|`$ に対し、$`-1`$ は $`tB`$ の固有値ではない（$`tBv=-v`$ ならば、$`v=0`$ でない限り
$`\|v\|=t\|Bv\|\le t\|B\|\|v\|<\|v\|`$ となる）から、$`t^np(B)=p(tB)=0`$ である。最後に斉次性により $`p(0)=0`$ である。$`\square`$

<a id="r-B-5"></a>

**補題 B.5（ニュートンの恒等式）.** $`x_1,\dots,x_n`$ をある環の互いに可換な元とし、$`x_S=\prod_{m\in S}x_m`$、
$`e_k=\sum_{|S|=k}x_S`$（$`e_0=1`$）、$`p_j=\sum_mx_m^j`$ とおく。このとき $`1\le k\le n`$ に対し $`k\,e_k=\sum_{j=1}^k(-1)^{j-1}e_{k-j}\,p_j`$ である。

*証明.* $`1\le j\le k`$ に対し $`a_j=\sum x_Sx_m^j`$ とおく。和は $`|S|=k-j`$ かつ $`m\notin S`$ を満たす組 $`(S,m)`$ にわたる。また
$`a_{k+1}=0`$ とする。$`e_{k-j}p_j=\sum_{|S|=k-j}\sum_mx_Sx_m^j`$ において、$`m\notin S`$ の項は $`a_j`$ を与え、$`m\in S`$ の項は
$`S'=S\setminus\{m\}`$ を用いて $`x_{S'}x_m^{j+1}`$ と書け、$`a_{j+1}`$ を与える。よって交代和は畳み込まれて $`a_1`$ となり、
$`a_1=\sum_{|T|=k}|T|\,x_T=k\,e_k`$ である（各 $`k`$ 元集合 $`T`$ は $`m\in T`$ の $`k`$ 通りの選び方から生じる）。$`\square`$

*計算：* $`x_m=A^{(m)}`$、$`n=4`$ に対する補題 B.5。これにより $`d\Pi(A^j)`$ から $`A^{\otimes4}`$ が復元される — [チュートリアル、B.3節](../../tutorial/ja/B_schur_weyl.md#s-B-3)。

<a id="r-B-6"></a>

**定理 B.6（シューア・ワイル双対性）.** $`d,n\ge1`$ とする。
(a) $`\pi(S_n)'=\mathcal T_n=\mathcal U_n=\mathrm{span}\{V^{\otimes n}:V\in SU(d)\}`$。この空間は $`*`$-環であるから、
$`V^{\otimes n}`$（$`V\in U(d)`$）が生成する環である。
(b) 張る空間が $`M_d`$ となる任意の集合 $`\mathcal X\subseteq M_d`$ に対し、
$`\{d\Pi(X):X\in\mathcal X\}`$ が生成する単位元をもつ環は $`\mathcal T_n`$ である。
(c) $`\{V^{\otimes n}:V\in U(d)\}'=\{V^{\otimes n}:V\in SU(d)\}'=\{d\Pi(X):X\in M_d\}'=\mathrm{span}\,\pi(S_n)`$。

*証明.* (a) $`\pi(S_n)'=\mathcal T_n`$ は命題 B.3 であり、$`\mathcal U_n\subseteq\mathcal T_n`$ である。$`\mathcal U_n\neq\mathcal T_n`$ と仮定する。すると
ある $`A_0^{\otimes n}\notin\mathcal U_n`$ が存在する。$`Y\neq0`$ を $`A_0^{\otimes n}`$ の $`\mathcal U_n`$ に直交する成分とする。関数
$`p(A)=\langle Y,A^{\otimes n}\rangle`$ は $`A`$ の成分の $`n`$ 次斉次多項式であり（$`A^{\otimes n}`$ の各成分は $`A`$ の $`n`$ 個の成分の積である）、
$`U(d)`$ 上で消え、しかも $`p(A_0)=\|Y\|^2\neq0`$ である。これは補題 B.4 に矛盾する。
$`V\in U(d)`$ に対し $`\det V=e^{i\theta}`$ となる $`\theta`$ を選ぶと、$`W=e^{-i\theta/d}V\in SU(d)`$ かつ $`V^{\otimes n}=e^{in\theta/d}W^{\otimes n}`$
であるから、二つの張る空間は一致する。$`\mathcal U_n`$ は $`1=1^{\otimes n}`$、$`V^{\otimes n}W^{\otimes n}=(VW)^{\otimes n}`$、
$`(V^{\otimes n})^*=(V^*)^{\otimes n}`$ を含むので $`*`$-環であり、$`V^{\otimes n}`$ を含む任意の環はそれらの張る空間を含む。

(b) $`\mathcal A`$ をその環とする。$`d\Pi`$ は線形であるから、$`\mathcal A`$ はすべての $`X\in M_d`$ に対する $`d\Pi(X)`$ を含む。補題 B.1(b) により
各 $`d\Pi(X)`$ は単位元をもつ環 $`\pi(S_n)'=\mathcal T_n`$ に属するので、$`\mathcal A\subseteq\mathcal T_n`$ である。逆に $`A\in M_d`$ を固定し、
互いに可換な $`x_m=A^{(m)}`$ に補題 B.5 を適用する：$`e_n=A^{\otimes n}`$ かつ $`p_j=d\Pi(A^j)\in\mathcal A`$ であり、$`k`$ に関する
帰納法によりすべての $`e_k=\frac1k\sum_j(-1)^{j-1}e_{k-j}p_j`$ は $`\mathcal A`$ に属する。よって $`A^{\otimes n}\in\mathcal A`$ であり、
$`\mathcal T_n\subseteq\mathcal A`$ である。

(c) 集合の可換子環は、それが生成する環の可換子環に等しい。(a) と (b) により三つの集合はいずれも
$`\mathcal T_n=\pi(S_n)'`$ を生成し、その可換子環は $`\pi(S_n)''`$ である。$`\pi(S_n)`$ の張る空間は $`*`$-環である
（$`\pi(\sigma)\pi(\tau)=\pi(\sigma\tau)`$、$`\pi(\sigma)^*=\pi(\sigma^{-1})`$、$`\pi(e)=1`$）から、定理 2.2 により
$`\pi(S_n)''=\mathrm{span}\,\pi(S_n)`$ である。$`\square`$

*計算：* 小さな $`d`$ と $`n`$ に対する定理 B.6(a) と (c) の次元 — [チュートリアル、B.3節](../../tutorial/ja/B_schur_weyl.md#s-B-3)。

<a id="r-B-7"></a>

**系 B.7（定理 4.14(b)）.** $`L`$ 個のスピン $`\frac12`$ の全スピン $`S_a=\sum_is_a^{(i)}`$ に対し、$`S_x,S_y,S_z`$ が生成する環は
$`\{S_a\}''=\pi(S_L)'`$ であり、$`\{S_a\}'=\mathrm{span}\,\pi(S_L)`$ である。

*証明.* $`M_2=\mathrm{span}\{1,s_x,s_y,s_z\}`$（パウリ行列と $`1`$ は基底をなす）、$`d\Pi(1)=L\cdot1`$、$`d\Pi(s_a)=S_a`$ である。
定理 B.6(b) により、$`S_x,S_y,S_z`$ が生成する単位元をもつ環 $`\mathcal A`$ は $`\mathcal T_L=\pi(S_L)'`$ である。$`S_a`$ はエルミートであるから
$`\mathcal A`$ は $`*`$-環であり（生成元の積の随伴は逆順の積である）、定理 2.2 により
$`\mathcal A=\mathcal A''\supseteq\{S_a\}''`$ である。また $`\{S_a\}''`$ は $`S_a`$ を含む単位元をもつ環であるから $`\{S_a\}''\supseteq\mathcal A`$ である。
最後に $`\{S_a\}'=\{S_a\}'''=\pi(S_L)''=\mathrm{span}\,\pi(S_L)`$ である（2.1節と定理 B.6(c)）。$`\square`$

$`d=2`$ の定理 B.6(a) により、同じ環は回転 $`V^{\otimes L}`$（$`V\in SU(2)`$）で張られる。これには 4.5節で触れた
$`SU(2)`$ の指数型パラメータ表示は必要ない。

<a id="s-B-4"></a>

## B.4 $`(\mathbb C^d)^{\otimes n}`$ の分解

<a id="r-B-8"></a>

**定理 B.8.** $`C\ge1`$、$`\sum_ck_cm_c=d^n`$ を満たす整数 $`k_c,m_c\ge1`$、およびユニタリ
$`W:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathcal H_n`$ が存在して
```math
W^*\pi(\sigma)W=\bigoplus_{c=1}^C\rho_c(\sigma)\otimes1_{m_c},\qquad W^*V^{\otimes n}W=\bigoplus_{c=1}^C1_{k_c}\otimes\tau_c(V)
\qquad(\sigma\in S_n,\ V\in U(d)),
```
となる。ここで $`\rho_c`$ は $`\mathbb C^{k_c}`$ 上の $`S_n`$ の互いに同値でない既約ユニタリ表現であり、
$`\tau_c`$ は $`\mathbb C^{m_c}`$ 上の $`U(d)`$ の互いに同値でない既約ユニタリ表現である。さらに
$`W^*\mathrm{span}\,\pi(S_n)W=\bigoplus_cM_{k_c}\otimes1_{m_c}`$、$`W^*\mathcal U_nW=\bigoplus_c1_{k_c}\otimes M_{m_c}`$、
```math
\sum_ck_c^2=\dim\mathrm{span}\,\pi(S_n),\qquad\sum_cm_c^2=\binom{d^2+n-1}n,
```
であり、$`n\le d`$ ならば $`\sum_ck_c^2=n!`$ である。

*証明.* 有限群 $`S_n`$ に対する命題 2.7 と $`\mathcal M=\mathrm{span}\,\pi(S_n)=\pi(S_n)''`$ に対する定理 2.6 により、
$`W`$、述べた性質をもつ $`\rho_c`$、および $`W^*\pi(S_n)'W=\bigoplus_c1_{k_c}\otimes M_{m_c}`$ が得られる。定理 B.6(a) により
$`\pi(S_n)'=\mathcal U_n`$ である。よってある $`\tau_c(V)\in M_{m_c}`$ に対し $`W^*V^{\otimes n}W=\bigoplus_c1\otimes\tau_c(V)`$ であり、
$`V\mapsto V^{\otimes n}`$ がユニタリ表現であるから $`V\mapsto\tau_c(V)`$ もユニタリ表現である。$`V^{\otimes n}`$ は $`\mathcal U_n`$ を張るので、線形写像
$`\sum_la_lV_l^{\otimes n}\mapsto(\sum_la_l\tau_c(V_l))_c`$ は $`\bigoplus_cM_{m_c}`$ への全射である。$`\tau_c`$ の既約性と非同値性は
命題 2.7 の証明とまったく同様に従う。次元は定理 2.6、命題 B.3、補題 B.1(d) による。$`\square`$

*計算：* `decompose.from_group` で分解した $`(\mathbb C^3)^{\otimes3}`$ 上の $`S_3`$ に対する定理 B.8 — [チュートリアル、B.4節](../../tutorial/ja/B_schur_weyl.md#s-B-4)。

$`d=2`$ の場合、セクターは全スピンによって識別される：

<a id="r-B-9"></a>

**系 B.9（$`L`$ 個のスピン $`\frac12`$）.** $`V`$ を、$`L`$ 個のスピン $`\frac12`$ の全スピンに対して
$`V^*S_aV=\bigoplus_SS^{(S)}_a\otimes1_{m_S}`$ を満たすユニタリとする（系 4.12、命題 4.13）。このとき
$`V^*\pi(\sigma)V=\bigoplus_S1_{2S+1}\otimes\rho_S(\sigma)`$ であり、ここで $`\rho_S`$ は $`S_L`$ の互いに同値でない既約ユニタリ表現で、
その次元は
```math
m_S=\binom Lj-\binom L{j-1},\qquad j=\tfrac L2-S\in\{0,1,\dots,\lfloor L/2\rfloor\}.
```
である。特に $`\dim\mathrm{span}\,\pi(S_L)=\sum_Sm_S^2`$ かつ $`\dim\pi(S_L)'=\sum_S(2S+1)^2=\binom{L+3}3`$ である。

*証明.* 定理 4.3 により $`V^*\{S_a\}'V=\bigoplus_S1_{2S+1}\otimes M_{m_S}`$ であり、系 B.7 により $`\{S_a\}'=\mathrm{span}\,\pi(S_L)`$ である。よって
ユニタリ表現 $`\rho_S`$ により $`V^*\pi(\sigma)V=\bigoplus_S1\otimes\rho_S(\sigma)`$ となり、写像
$`\sum_\sigma a_\sigma\pi(\sigma)\mapsto(\sum_\sigma a_\sigma\rho_S(\sigma))_S`$ は $`\bigoplus_SM_{m_S}`$ への全射である。既約性と非同値性は
命題 2.7 の証明と同様に従う。$`m_S`$ の値は命題 4.13 による。定理 4.3 により
$`\dim\{S_a\}'=\sum_Sm_S^2`$ かつ $`\dim\{S_a\}''=\sum_S(2S+1)^2`$ である。系 B.7 により $`\{S_a\}'=\mathrm{span}\,\pi(S_L)`$ かつ
$`\{S_a\}''=\pi(S_L)'`$ であり、後者の次元は $`d=2`$ の命題 B.3 により $`\binom{4+L-1}L=\binom{L+3}3`$ である。$`\square`$

*計算：* $`L=5,6`$ に対する組 $`(m_S,2S+1)`$ — [チュートリアル、B.4節](../../tutorial/ja/B_schur_weyl.md#s-B-4)。

<a id="r-B-10"></a>

**系 B.10（カタラン数）.** $`L\ge1`$ に対し $`\displaystyle\sum_Sm_S^2=\binom{2L}L\frac1{L+1}`$ である。これは
$`(\mathbb C^2)^{\otimes L}`$ 上の $`\dim\{S_a\}'=\dim\mathrm{span}\,\pi(S_L)`$ である。

*証明.* $`N=2^L`$ とし、$`\mathrm{vec}:M_N\to\mathbb C^N\otimes\mathbb C^N`$ を線形全単射
$`|x\rangle\langle y|\mapsto x\otimes\bar y`$ とする。すると $`\mathrm{vec}(AX-XA)=(A\otimes1-1\otimes A^T)\mathrm{vec}\,X`$ である。実際
$`A|x\rangle\langle y|B=|Ax\rangle\langle B^*y|\mapsto Ax\otimes B^T\bar y`$ である。よって $`\{S_a\}'`$ は
$`S_a\otimes1-1\otimes S_a^T`$（$`a=x,y,z`$）の共通核の上に写される。ここで $`S_a^T=\sum_i(s_a^T)^{(i)}`$ である。$`s_x^T=s_x`$、$`s_y^T=-s_y`$、$`s_z^T=s_z`$ と
$`\sigma_ys_x\sigma_y=-s_x`$、$`\sigma_ys_y\sigma_y=s_y`$、$`\sigma_ys_z\sigma_y=-s_z`$ から $`\sigma_y(-s_a^T)\sigma_y=s_a`$ を得るので、ユニタリ
$`Y=1\otimes\sigma_y^{\otimes L}`$ は $`Y(S_a\otimes1-1\otimes S_a^T)Y^*=S_a\otimes1+1\otimes S_a`$ を与える。これは
$`(\mathbb C^2)^{\otimes2L}=\mathbb C^N\otimes\mathbb C^N`$ 上の $`2L`$ 個のスピン $`\frac12`$ の全スピンである。$`2L`$ 個のスピンに対する系 4.12 の
分解において、三つの成分すべてで零化されるベクトルは $`\mathbf S^2`$ で零化されるが、これはスピン
$`S\ge\frac12`$ の直和成分上で $`S(S+1)\ne0`$ である。逆に、各成分はスピン $`0`$ の直和成分上で消える。よって共通核はその直和成分であり、
その次元は $`m_0(2L)=\binom{2L}L-\binom{2L}{L-1}=\binom{2L}L\big(1-\frac L{L+1}\big)`$ である（命題 4.13）。次元の第一の等式は
定理 4.3、第二の等式は系 B.7 による。$`\square`$

*計算：* $`L=5,6`$ に対する系 B.10 — [チュートリアル、B.4節](../../tutorial/ja/B_schur_weyl.md#s-B-4)。

<a id="s-B-5"></a>

## B.5 ワインガルテン計算

<a id="r-B-11"></a>

**補題 B.11（ハール・ツワールは可換子環への直交射影である）.** $`U`$ を $`U(d)`$ 上でハール分布に従うとし、
$`X\in\mathrm{End}\,\mathcal H_n`$ に対し $`\Phi(X)=\mathbb E[U^{\otimes n}X(U^{\otimes n})^*]`$ とおく。このとき $`\Phi`$ は
$`\mathrm{End}\,\mathcal H_n`$ から $`\{V^{\otimes n}:V\in U(d)\}'=\mathrm{span}\,\pi(S_n)`$ への直交射影である。

*証明.* $`U^{\otimes n}X(U^{\otimes n})^*`$ の成分は $`U`$ と $`\bar U`$ の成分の多項式であり、それらの成分は絶対値が $`1`$ 以下であるから、
$`\Phi`$ はきちんと定義された線形写像である。$`V\in U(d)`$ に対し $`VU`$ は $`U`$ と同じ分布をもつから
$`V^{\otimes n}\Phi(X)(V^{\otimes n})^*=\mathbb E[(VU)^{\otimes n}X((VU)^{\otimes n})^*]=\Phi(X)`$ であり、像は可換子環に含まれる。$`Y`$ が
可換子環に属するならば $`U^{\otimes n}Y(U^{\otimes n})^*=Y`$ であり $`\Phi(Y)=Y`$ である。よって $`\Phi`$ は冪等で、その像は可換子環である。
補題 5.1 により $`U^*=U^{-1}`$ は $`U`$ と同じ分布をもつので、
$`\langle Y,\Phi(X)\rangle=\mathbb E\mathrm{Tr}\big(((U^*)^{\otimes n}YU^{\otimes n})^*X\big)=\langle\Phi(Y),X\rangle`$ であり、$`\Phi`$ は自己共役である。
自己共役な冪等写像はその像への直交射影である
（$`\langle X-\Phi(X),\Phi(Z)\rangle=\langle\Phi(X-\Phi(X)),Z\rangle=0`$）。定理 B.6(c) により可換子環は $`\mathrm{span}\,\pi(S_n)`$ である。$`\square`$

*計算：* ハール・ユニタリにわたるモンテカルロ法による補題 B.11 — [チュートリアル、B.5節](../../tutorial/ja/B_schur_weyl.md#s-B-5)。

この射影を、$`n>d`$ では一次従属となる張る族 $`\pi(\sigma)`$ を用いて書くために、エルミート行列 $`G`$ の
**擬似逆行列**を用いる：$`\lambda\neq0`$ に対し $`f(\lambda)=1/\lambda`$、$`f(0)=0`$ として $`G^+=f(G)`$ とする
（関数計算、定理 1.5）。これはエルミートであり、$`G`$ が可逆ならば $`G^+=G^{-1}`$、$`GG^+=\mathbf 1_{\neq0}(G)`$ は
$`(\ker G)^\perp`$ への射影、$`\lambda\neq0`$ で $`Gv=\lambda v`$ ならば $`G^+v=v/\lambda`$、そしてユニタリ $`R`$ に対し
$`(R^*GR)^+=R^*G^+R`$ である（$`G=Q\Lambda Q^*`$ ならば $`R^*GR=(R^*Q)\Lambda(R^*Q)^*`$）。

<a id="r-B-12"></a>

**補題 B.12.** $`v_1,\dots,v_N`$ を有限次元内積空間 $`\mathcal K`$ のベクトルとし、$`G_{ab}=\langle v_a,v_b\rangle`$ を
そのグラム行列、$`W=G^+`$ とする。このとき $`\mathrm{span}\{v_a\}`$ への直交射影は
$`x\mapsto\sum_{a,b}v_a\,W_{ab}\,\langle v_b,x\rangle`$ である。

*証明.* $`T:\mathbb C^N\to\mathcal K`$、$`Tc=\sum_ac_av_a`$ とする。すると $`T^*x=(\langle v_b,x\rangle)_b`$、$`G=T^*T`$ であり、問題の写像は
$`Q=TG^+T^*`$ である。$`\langle c,Gc\rangle=\|Tc\|^2`$ であるから $`\ker G=\ker T`$ であり、したがって $`\ker T`$ への射影 $`\Pi`$ に対し
$`TG^+G=T(1-\Pi)=T`$ である。よって $`QT=T`$ であり、$`Q`$ は $`\mathrm{ran}\,T`$ 上で恒等写像である。$`Q`$ は $`(\mathrm{ran}\,T)^\perp=\ker T^*`$ 上で消え、
$`Q^*=Q`$ である。よって $`Q`$ は $`\mathrm{ran}\,T`$ への直交射影である。$`\square`$

<a id="r-B-13"></a>

**定理 B.13（任意の次数のワインガルテン公式）.** $`n,d\ge1`$ とし、$`U`$ を $`U(d)`$ 上でハール分布に従うとし、$`G`$ を
$`n!\times n!`$ 行列 $`G_{\sigma\tau}=d^{\#\text{cycles}(\sigma^{-1}\tau)}`$（$`\sigma,\tau\in S_n`$）、$`W=G^+`$ とする。
(a) $`S_n`$ 上の関数 $`\mathrm{Wg}(\rho)=\mathrm{Wg}(\rho,d)=W_{e\rho}`$ により $`W_{\sigma\tau}=\mathrm{Wg}(\sigma^{-1}\tau)`$ であり、この関数は
$`\rho`$ の巡回型のみに依存する。特に $`\mathrm{Wg}(\rho^{-1})=\mathrm{Wg}(\rho)`$ である。
(b) $`G`$ が可逆であるのは $`n\le d`$ のときかつそのときに限り、そのとき $`W=G^{-1}`$ である。
(c) $`X\in\mathrm{End}\,\mathcal H_n`$ に対し
$`\mathbb E[U^{\otimes n}X(U^{\otimes n})^*]=\sum_{\sigma,\tau}\mathrm{Wg}(\sigma^{-1}\tau)\mathrm{Tr}(\pi(\tau)^*X)\,\pi(\sigma)`$。
(d) すべての $`\mathbf i,\mathbf j,\mathbf i',\mathbf j'\in[d]^n`$ に対し
```math
\mathbb E\big[U_{i_1j_1}\cdots U_{i_nj_n}\overline{U_{i'_1j'_1}\cdots U_{i'_nj'_n}}\big]
=\sum_{\sigma,\tau\in S_n}\prod_{m=1}^n\delta_{i_mi'_{\sigma(m)}}\delta_{j_mj'_{\tau(m)}}\;\mathrm{Wg}(\sigma\tau^{-1}),
```
すなわち和は $`\mathbf i=\mathbf i'\circ\sigma`$ かつ $`\mathbf j=\mathbf j'\circ\tau`$ を満たす組にわたる。
(e) 因子 $`U`$ と $`\bar U`$ の個数 $`p\neq q`$ が異なるモーメントは零である。

これは（$`k=d`$ とした）定理 5.12 であり、ここに証明された。

*証明.* (a) $`g(\rho)=d^{\#\text{cycles}(\rho)}`$ と書くと $`G_{\sigma\tau}=g(\sigma^{-1}\tau)`$ である。$`\Lambda_\kappa,R_\kappa`$ を
$`\mathbb C^{S_n}`$ 上の置換行列で $`\Lambda_\kappa e_\sigma=e_{\kappa\sigma}`$、$`R_\kappa e_\sigma=e_{\sigma\kappa}`$ を満たすものとする。すると
$`(\Lambda_\kappa^*G\Lambda_\kappa)_{\sigma\tau}=g(\sigma^{-1}\kappa^{-1}\kappa\tau)=G_{\sigma\tau}`$ かつ
$`(R_\kappa^*GR_\kappa)_{\sigma\tau}=g(\kappa^{-1}\sigma^{-1}\tau\kappa)=G_{\sigma\tau}`$ である（補題 B.1(c)）。擬似逆行列の性質により、
$`W`$ も同じ二つの不変性をもつ：$`W_{\kappa\sigma,\kappa\tau}=W_{\sigma\tau}`$ かつ
$`W_{\sigma\kappa,\tau\kappa}=W_{\sigma\tau}`$。第一のものに $`\kappa=\sigma^{-1}`$ を代入すると $`W_{\sigma\tau}=\mathrm{Wg}(\sigma^{-1}\tau)`$ を得る。第二のものは
$`\mathrm{Wg}(\kappa^{-1}\rho\kappa)=\mathrm{Wg}(\rho)`$ を与えるので、$`\mathrm{Wg}`$ は共役類の上で一定であり、共役類は巡回型である
（補題 B.1(c)）。$`\rho^{-1}`$ は $`\rho`$ と同じ巡回型をもつ。
(b) $`G`$ は $`\pi(\sigma)`$ たちのグラム行列であり（補題 B.1(c)）、補題 B.12 の証明の $`T`$ を用いると、$`G=T^*T`$ が
可逆であるのは $`T`$ が単射であるとき、すなわち $`\pi(\sigma)`$ たちが一次独立であるときかつそのときに限る。補題 B.1(d) によりそれは $`n\le d`$ と同値である。
(c) $`\mathcal K=\mathrm{End}\,\mathcal H_n`$、$`v_\sigma=\pi(\sigma)`$ として補題 B.11 と補題 B.12 を用いる。
(d) クロネッカー積の成分は $`\langle e_{\mathbf i},U^{\otimes n}e_{\mathbf j}\rangle=\prod_mU_{i_mj_m}`$ であるから、モーメントは
$`X=|e_{\mathbf j}\rangle\langle e_{\mathbf j'}|`$ に対する $`\langle e_{\mathbf i},\Phi(X)e_{\mathbf i'}\rangle`$ である。ここで
$`\mathrm{Tr}(\pi(\tau)^*X)=\langle\pi(\tau)e_{\mathbf j'},e_{\mathbf j}\rangle=\delta_{\mathbf j,\mathbf j'\circ\tau^{-1}}`$ かつ
$`\langle e_{\mathbf i},\pi(\sigma)e_{\mathbf i'}\rangle=\delta_{\mathbf i,\mathbf i'\circ\sigma^{-1}}`$ である。(c) によりモーメントは
$`\sum_{\sigma,\tau}\delta_{\mathbf i,\mathbf i'\circ\sigma^{-1}}\delta_{\mathbf j,\mathbf j'\circ\tau^{-1}}\mathrm{Wg}(\sigma^{-1}\tau)`$ である。
$`\sigma\to\sigma^{-1}`$、$`\tau\to\tau^{-1}`$ と置き換えれば公式を得る。
(e) 実数 $`\theta`$ に対し $`e^{i\theta}1\in U(d)`$ であるから $`e^{i\theta}U`$ は $`U`$ と同じ分布をもち、モーメントはそれ自身の
$`e^{i(p-q)\theta}`$ 倍に等しい。$`e^{i(p-q)\theta}\neq1`$ となる $`\theta`$ をとればよい。$`\square`$

$`n=2\le d`$ のとき $`G=\begin{pmatrix}d^2&d\\d&d^2\end{pmatrix}`$ であり、$`G^{-1}`$ から $`\mathrm{Wg}(e)=\frac1{d^2-1}`$、
$`\mathrm{Wg}(\tau)=-\frac1{d(d^2-1)}`$ を得る：定理 5.10 である。$`d=1`$ のときすべての $`G_{\sigma\tau}=1`$、$`W=\frac1{(n!)^2}G`$ であり、
位相 $`U=U_{11}`$ のモーメントは $`\sum_{\sigma,\tau}\frac1{(n!)^2}=1=|U_{11}|^{2n}`$ となり、当然そうあるべき値である。

*計算：* 擬似逆行列から $`\mathrm{Wg}`$ を計算した定理 B.13。$`G`$ が特異となる $`d=2`$ での次数 $`3`$ のモーメントを含む — [チュートリアル、B.5節](../../tutorial/ja/B_schur_weyl.md#s-B-5)。

<a id="r-B-14"></a>

**系 B.14.** すべての $`n,d\ge1`$ に対し $`\mathbb E|U_{11}|^{2n}=\dfrac{n!}{d(d+1)\cdots(d+n-1)}`$ である。

*証明.* 定理 B.13(d) ですべての添字を $`1`$ とすると、モーメントは $`\sum_{\sigma,\tau}\mathrm{Wg}(\sigma\tau^{-1})=\sum_{\sigma,\tau}W_{\sigma^{-1},\tau^{-1}}=\mathbf 1^TW\mathbf 1`$ である（(a)）。$`G`$ の各行和は
$`c_n(d)=\sum_{\rho\in S_n}d^{\#\text{cycles}(\rho)}\ge d^n>0`$ であるから、$`G\mathbf 1=c_n(d)\mathbf 1`$、$`W\mathbf 1=\mathbf 1/c_n(d)`$ であり、
モーメントは $`n!/c_n(d)`$ である。最後に $`c_n(x)=x(x+1)\cdots(x+n-1)`$ である（第5章の演習問題 2 の解答）。$`\square`$

これは第5章の演習問題 2 を、$`G`$ が特異となる $`n>d`$ の場合に拡張するものである。

*計算：* $`n=3`$、$`d=2`$ に対する系 B.14 — [チュートリアル、B.5節](../../tutorial/ja/B_schur_weyl.md#s-B-5)。

<a id="s-B-6"></a>

## B.6 テンソル作用素と $`SU(2)`$ に対するウィグナー・エッカルトの定理

9.3節ではウィグナー・エッカルトの定理からスピン相関を構成する。ここでは 4.5節のスピン表現に対してこれを証明する。
スピン $`k`$ とは、基底 $`|k,q\rangle`$、$`q=k,\dots,-k`$ をもつ $`\mathbb C^{2k+1}`$ 上の定理 4.11 の表現 $`S^{(k)}_a`$ である。
$`S^{(k)}_\pm|k,q\rangle=\sqrt{k(k+1)-q(q\pm1)}\,|k,q\pm1\rangle`$ であったことを思い出そう（$`S_-=S_+^*`$ についてはこれは $`S_+`$ の
公式から従う）。

**定義.** $`(S_a)`$ を $`\mathbb C^N`$ 上のスピン表現とし、$`k\in\{0,\frac12,1,\dots\}`$ とする。**ランク $`k`$ のテンソル作用素**とは、
$`N\times N`$ 行列の族 $`(T_q)_{q=-k,\dots,k}`$ で
```math
[S_z,T_q]=q\,T_q,\qquad[S_\pm,T_q]=\sqrt{k(k+1)-q(q\pm1)}\;T_{q\pm1}
```
を満たすものである（$`q=\pm k`$ では係数が消え、$`T_{\pm(k+1)}`$ は定義されない）。**ベクトル作用素**とは、
$`[S_a,T_b]=i\sum_c\varepsilon_{abc}T_c`$、すなわち $`[S_x,T_y]=iT_z`$、$`[S_x,T_z]=-iT_y`$ とその巡回置換、および
$`[S_a,T_a]=0`$ を満たす三つ組 $`(T_x,T_y,T_z)`$ である。

<a id="r-B-15"></a>

**補題 B.15.** (a) $`(T_a)`$ がベクトル作用素ならば、その**球面成分** $`T_{\pm1}=\mp\frac1{\sqrt2}(T_x\pm iT_y)`$、
$`T_0=T_z`$ はランク $`1`$ のテンソル作用素をなし、二つのベクトル作用素に対し
$`\sum_aT_aT'_a=\sum_q(-1)^qT_qT'_{-q}`$ である。
(b) $`L`$ 個のスピン $`\frac12`$ と全スピン $`S_a`$ に対し、各サイトのスピン $`(s^{(i)}_a)_a`$ および $`(S_a)_a`$ 自身はベクトル作用素である。
(c) $`(T_q)`$ がランク $`k`$ のテンソル作用素であるのは、写像 $`\Psi:\mathbb C^N\otimes\mathbb C^{2k+1}\to\mathbb C^N`$、
$`\Psi(v\otimes|k,q\rangle)=T_qv`$ が $`a=x,y,z`$ に対し $`\Psi\,(S_a\otimes1+1\otimes S^{(k)}_a)=S_a\Psi`$ を満たすときかつそのときに限る。

*証明.* (a) 関係式から $`[S_z,T_x\pm iT_y]=\pm(T_x\pm iT_y)`$、$`[S_\pm,T_z]=\mp(T_x\pm iT_y)`$、
$`[S_\pm,T_x\mp iT_y]=\pm2T_z`$、$`[S_\pm,T_x\pm iT_y]=0`$ を得る（たとえば $`[S_x+iS_y,T_x-iT_y]=-i[S_x,T_y]+i[S_y,T_x]=T_z+T_z`$）。
球面成分で書けば $`[S_z,T_q]=qT_q`$、$`[S_\pm,T_0]=\sqrt2\,T_{\pm1}`$、$`[S_\pm,T_{\mp1}]=\sqrt2\,T_0`$、$`[S_\pm,T_{\pm1}]=0`$ であり、
これは $`k=1`$ の定義関係である（$`\sqrt{2-q(q\pm1)}`$ は $`q=0`$ と $`q=\mp1`$ で $`\sqrt2`$、$`q=\pm1`$ で $`0`$）。
積については $`T_{+1}T'_{-1}+T_{-1}T'_{+1}=-\frac12\big((T_x+iT_y)(T'_x-iT'_y)+(T_x-iT_y)(T'_x+iT'_y)\big)=-(T_xT'_x+T_yT'_y)`$ である。
(b) $`[S_a,s^{(i)}_b]=[s^{(i)}_a,s^{(i)}_b]`$（他のサイトは可換）であり、スピン $`\frac12`$ の行列はスピンの関係式を満たす。
$`S_a`$ については $`i`$ について和をとる。(c) $`\Psi(S_a\otimes1+1\otimes S^{(k)}_a)(v\otimes|k,q\rangle)=T_qS_av+\sum_{q'}(S^{(k)}_a)_{q'q}T_{q'}v`$ であるから、
$`S_a`$ に対する絡関係は、すべての $`q`$ に対し $`[S_a,T_q]=\sum_{q'}(S^{(k)}_a)_{q'q}T_{q'}`$ を意味する。$`S_z,S_\pm`$ に対してはこれらは
定義関係であり、$`S_x,S_y,S_z`$ に対する関係式と $`S_z,S_\pm`$ に対する関係式は互いの一次結合である。$`\square`$

*計算：* 六つのスピン $`\frac12`$ のサイトスピンに対する補題 B.15(a) と (b) — [チュートリアル、B.6節](../../tutorial/ja/B_schur_weyl.md#s-B-6)。

<a id="r-B-16"></a>

**定理 B.16（ウィグナー・エッカルト）.** $`(S_a)`$ を $`\mathbb C^N`$ 上のスピン表現、$`V`$ を
$`V^*S_aV=\bigoplus_SS^{(S)}_a\otimes1_{m_S}`$ を満たすユニタリ（系 4.12）、$`|S,M;b\rangle=V(|S,M\rangle\otimes f_b)`$ を対応する正規直交
基底（$`f_b`$ は $`\mathbb C^{m_S}`$ の標準基底）、$`(T_q)`$ をランク $`k`$ のテンソル作用素とする。$`S,S'`$ は現れるとする（$`m_S,m_{S'}\ge1`$）。
(a) $`S'\notin\{|S-k|,|S-k|+1,\dots,S+k\}`$ ならば、すべての $`M',M,q,a,b`$ に対し $`\langle S',M';a|T_q|S,M;b\rangle=0`$ である。
(b) そうでなければ、スピン $`S\otimes`$ スピン $`k`$ からスピン $`S'`$ への絡作用素は一次元の空間をなす。そのような零でない
$`C=C^{S'}_{S,k}:\mathbb C^{2S+1}\otimes\mathbb C^{2k+1}\to\mathbb C^{2S'+1}`$ のそれぞれに対し、$`m_{S'}\times m_S`$ 行列
$`R=\langle S'\|T\|S\rangle`$、すなわち**換算行列要素**が一意に存在して
```math
\langle S',M';a|\,T_q\,|S,M;b\rangle=\langle S',M'|\,C\,(|S,M\rangle\otimes|k,q\rangle)\rangle\;R_{ab}\qquad\text{for all }M',M,q,a,b .
```
が成り立つ。
(c) この要素は $`M'=M+q`$ でない限り消える。(b) の場合、$`\langle S',S'|C(|S,S'-q\rangle\otimes|k,q\rangle)\rangle\neq0`$ となる $`q`$ が存在する。
したがって $`R`$、そしてそれとともに多重項 $`S`$ と $`S'`$ の間の $`T`$ のブロック全体は、$`S'`$ の最高ウェイトベクトルと $`S`$ の各コピーの
一つのベクトルとの間の行列要素 $`\langle S',S';a|T_q|S,S'-q;b\rangle`$ によって決まる。

*証明.* $`a\le m_{S'}`$、$`b\le m_S`$ を固定し、$`J:\mathbb C^{2S+1}\to\mathbb C^N`$、$`Jx=V(x\otimes f_b)`$（直和成分 $`S`$ の中）とし、
$`K:\mathbb C^N\to\mathbb C^{2S'+1}`$ を $`x\mapsto V(x\otimes f_a)`$（直和成分 $`S'`$ の中）の随伴とする。成分の添字を
$`\alpha\in\{x,y,z\}`$ と書く。$`V^*S_\alpha V`$ の形から $`S_\alpha J=JS^{(S)}_\alpha`$ であり、$`S'`$ に対する同じ関係の随伴をとると
（スピン行列はエルミートである）$`KS_\alpha=S^{(S')}_\alpha K`$ である。補題 B.15(c) の $`\Psi`$ を用いて $`t=K\Psi(J\otimes1)`$ とおく。
その行列要素は $`\langle S',M'|t(|S,M\rangle\otimes|k,q\rangle)\rangle=\langle S',M';a|T_q|S,M;b\rangle`$ であり、
```math
t\,(S^{(S)}_\alpha\otimes1+1\otimes S^{(k)}_\alpha)=K\Psi(S_\alpha\otimes1+1\otimes S^{(k)}_\alpha)(J\otimes1)=KS_\alpha\Psi(J\otimes1)=S^{(S')}_\alpha\,t,
```
であるから、$`t`$ は $`\rho=`$ スピン $`S\otimes`$ スピン $`k`$（スピン表現
$`S^{(S)}_\alpha\otimes1+1\otimes S^{(k)}_\alpha`$、4.5節）から $`\sigma=`$ スピン $`S'`$ への絡作用素である。随伴をとる操作は
$`\mathrm{Hom}(\rho,\sigma)`$ と $`\mathrm{Hom}(\sigma,\rho)`$ の間の線形全単射である（$`t\rho(\alpha)=\sigma(\alpha)t`$ の随伴は
$`\rho(\alpha)t^*=t^*\sigma(\alpha)`$ であり、関係する行列はすべてエルミートである）。系 4.12 と第4章の演習問題 2（その解答では
系 4.12 のみから証明される）により、$`\rho`$ は $`\bigoplus_j(\text{spin }j)`$ にユニタリ同値であり、各
$`j\in\{|S-k|,|S-k|+1,\dots,S+k\}`$ がちょうど一つずつ現れる。これらのスピンは既約で互いに同値でない（定理 4.11）から、これは
定理 4.3 の意味での分解であり、その最後の主張により、$`S'`$ がこれらの $`j`$ の一つならば $`\dim\mathrm{Hom}(\sigma,\rho)=1`$、
そうでなければ $`\mathrm{Hom}(\sigma,\rho)=0`$ である（このときスピン $`S'`$ はそのいずれとも同値でない、定理 4.11）。
これで (a) が得られ、(b) の場合には一意な数 $`R_{ab}`$ に対し $`t=R_{ab}C`$ となる。
(c) $`S_z`$ に対する絡関係は $`M'\langle S',M'|t(|S,M\rangle\otimes|k,q\rangle)\rangle=(M+q)\langle S',M'|t(|S,M\rangle\otimes|k,q\rangle)\rangle`$ を与える。
$`C\neq0`$ の像は既約なスピン $`S'`$ の零でない不変部分空間であるから全体に等しく、したがってある基底ベクトルに対し
$`\langle S',S'|C(|S,M\rangle\otimes|k,q\rangle)\rangle\neq0`$ である。このとき絡作用素 $`C`$ に対する同じ関係により $`M=S'-q`$ である。
この数で割ればよい。$`\square`$

*計算：* `decompose.from_generators` が計算した分解における、基底によらない形の定理 B.16(a) と (b) — [チュートリアル、B.6節](../../tutorial/ja/B_schur_weyl.md#s-B-6)。

補題 B.15(a) により、
$`\mathbf s^{(i)}\cdot\mathbf s^{(j)}=\sum_q(-1)^q s^{(i)}_qs^{(j)}_{-q}`$ のようなスカラー積の多重項 $`S`$ と $`S'`$ の間のブロックは
$`\sum_{S''}\sum_q(-1)^q(s^{(i)}_q\text{ の }S''\to S'\text{ ブロック})(s^{(j)}_{-q}\text{ の }S\to S''\text{ ブロック})`$ であり、各因子は
定理 B.16 により換算行列要素とクレブシュ・ゴルダン写像を通じて与えられる。

*ライブラリでは：* これが `su2.spin_correlation` の背後にある原理であり、第9章でそれを `'direct'` 法と比較して確かめる — [チュートリアル、9.3節](../../tutorial/ja/09_spin_chains.md#s-9-3)。

<a id="s-B-7"></a>

## B.7 本コースで用いるハール測度

本コースでは次の群のハール測度を用いる：$`U(k)`$（定理 5.3）、$`U(1)`$（$`k=1`$ の場合；第5章の演習問題 3）、
$`U(\mathcal M)\cong\prod_cU(k_c)`$（系 5.4）、$`SU(k)`$（5.1節、系 5.4 の後）、直交群
$`O(k)`$（5.1節、第11章の演習問題 2、および実環）、そして補題 C.2 の対角群 $`T`$ である。これらすべてに対する一意性は
補題 5.1 による。いずれもある $`U(n)`$ の閉部分群だからである。存在はそれぞれの場合に明示的な構成によって示され
（$`T`$ については補題 C.2）、$`O(k)`$ に対する構成は以下のとおりである。コンパクト群のハール測度の一般的な存在定理は、
本コースのどこでも用いない。

<a id="r-B-17"></a>

**定理 B.17（直交群のハール測度）.** $`Z`$ を、成分が独立に $`N(0,1)`$ に従う $`k\times k`$ 実行列とする。
このとき $`Z`$ は確率 1 で可逆であり、$`O`$ を直交行列、$`R`$ を対角成分が正の実上三角行列として $`Z=OR`$ と書け、
$`O`$ は $`O(k)`$ 上でハール分布に従う。$`Z=O'R'`$ が任意の実 QR 分解（$`O'`$ は直交、$`R'`$ は実上三角）ならば
$`O'\mathrm{diag}(\mathrm{sgn}\,R'_{ii})=O`$ である。ハール分布に従う $`O`$ と
$`X\in M_k(\mathbb C)`$ に対し $`\mathbb E[O]=0`$ かつ $`\mathbb E[OXO^T]=\frac{\mathrm{Tr}\,X}k1`$ であり、特に
$`\mathbb E[O_{ij}O_{lm}]=\delta_{il}\delta_{jm}/k`$ である。

*証明.* $`\det Z`$ は $`k^2`$ 個の実成分の多項式で恒等的に零ではなく（$`\det1=1`$）、$`Z`$ は密度をもつから、
ほとんど確実に $`\det Z\neq0`$ である（系 A.7(b)、定理 5.3 の証明と同様）。実の列にグラム・シュミットを適用すると実の $`q_j`$ と実の係数、
すなわち $`R_{jj}>0`$ をもつ実分解 $`Z=OR`$ が得られる。補題 5.2 の一意性によりこれは補題 5.2 の分解であり、したがって $`O`$ は
$`Z`$ に連続的に依存する。$`Z`$ の密度 $`(2\pi)^{-k^2/2}e^{-\|Z\|_2^2/2}`$ は $`V\in O(k)`$ に対する $`Z\mapsto VZ`$ のもとで不変である。
これは $`\|Z\|_2^2=\mathrm{Tr}(Z^TZ)`$ を保つ $`\mathbb R^{k^2}`$ の線形写像であるから直交変換であり、ヤコビアンは $`1`$ である。
よって $`VZ`$ は $`Z`$ と同じ分布をもつ（命題 A.14(a)）。$`VZ=(VO)R`$ であるから一意性により $`O(VZ)=VO(Z)`$ であり、したがって $`O`$ の分布は
コンパクト群 $`O(k)`$ 上で左不変であり、補題 5.1 が適用できる。符号の補正は定理 5.3 の証明の最後の部分であり、
$`\mathrm{sgn}\,R'_{ii}=R'_{ii}/|R'_{ii}|`$ である（$`Z`$ は可逆なので $`R'_{ii}\neq0`$）。$`-1\in O(k)`$ であるから $`-O`$ は $`O`$ と同じ分布をもち、
$`\mathbb E[O]=-\mathbb E[O]=0`$ である（成分の絶対値は $`1`$ 以下）。最後に、$`VO`$ は $`O`$ と同じ分布をもつから、$`T=\mathbb E[OXO^T]`$ は
$`V\in O(k)`$ に対し $`VTV^T=T`$ を満たす。$`V`$ を位置 $`i`$ に $`-1`$、それ以外に $`1`$ をもつ対角行列とすると、$`j\neq i`$ に対し $`T_{ij}=-T_{ij}`$
を得る。置換行列を用いると、すべての $`T_{ii}`$ が等しい（補題 5.6 と同様）。よって $`T=\lambda1`$ であり、$`k\lambda=\mathrm{Tr}\,T=\mathrm{Tr}\,X`$ である。
$`X=|e_j\rangle\langle e_m|`$ に対しては $`(OXO^T)_{il}=O_{ij}O_{lm}`$ である。$`\square`$

*計算：* モンテカルロ法による定理 B.17 の二次モーメント、および符号の補正が必要な理由 — [チュートリアル、B.7節](../../tutorial/ja/B_schur_weyl.md#s-B-7)。

**戻る：**[第4章](04_symmetries.md)（定理 4.14）、[第5章](05_random_matrices.md)（定理 5.12）、
[第9章](09_spin_chains.md)（相関）。
