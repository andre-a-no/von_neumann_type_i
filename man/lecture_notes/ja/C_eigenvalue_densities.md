[English](../C_eigenvalue_densities.md) · [Español](../es/C_eigenvalue_densities.md) · [Français](../fr/C_eigenvalue_densities.md) · [Deutsch](../de/C_eigenvalue_densities.md) · [中文](../zh/C_eigenvalue_densities.md) · **日本語** · [Русский](../ru/C_eigenvalue_densities.md)

# 付録 C. 固有値密度とペイジの公式

第5章はガウス・アンサンブルの固有値の同時密度を用い、第7章はランダムな純粋状態の平均エンタングルメント・エントロピーに
対するペイジの公式を用いる。この付録では両方を証明する。
道筋は次のとおりである：ワイル積分公式（「対角化」のヤコビアン）、ガウス・アンサンブル、
直交多項式による行列式型密度の1点関数、ウィシャート・ラゲール・アンサンブルと縮約状態上の誘導測度、そして最後に
ラゲール多項式による厳密な計算である。

**本付録で学ぶこと**：不変なランダム行列の固有値がなぜ
$`\prod_{i<j}|\lambda_i-\lambda_j|^\beta`$ のように反発するのか、そのとき固有ベクトルがなぜハール分布に従い固有値と
独立になるのか、正規直交多項式の行列式がどのように積分されるのか、そして
$`\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}`$ がいくつかの積分からどのように従うのか。

本講義の予備知識（測度論、可逆な導関数をもつ単射 $`C^1`$ 写像に対する
$`\mathbb R^N`$ における変数変換公式、逆関数定理と陰関数
定理）のほかに、第5章から補題 5.1（ハール測度の一意性と右不変性）、
補題 5.2 と定理 5.3（ガウス QR による $`U(n)`$ 上のハール測度の存在）、付録 B から
定理 B.17（$`O(n)`$ に対する同じ主張）、付録 A から命題 A.16(a)（複素
ガウスベクトルのユニタリ不変性）と定理 A.20(b)（$`|z|^2`$ の分布）、そして第7章から定理 7.11 の証明における
表示 $`\rho_A=GG^*/\mathrm{Tr}\,GG^*`$ を用いる。

**記号.** $`\beta\in\{1,2\}`$ とし、それに応じて $`\mathbb F=\mathbb R`$ または $`\mathbb C`$ とする。$`\mathcal H`$ は
実対称（$`\beta=1`$）または複素エルミート（$`\beta=2`$）な $`n\times n`$ 行列の空間であり、座標 $`h_{ii}`$、$`\mathrm{Re}\,h_{ij}`$ および（$`\beta=2`$ のとき）
$`\mathrm{Im}\,h_{ij}`$（$`i<j`$）に関するルベーグ測度 $`dH`$ をもつ。その次元は $`n+\beta n(n-1)/2`$ である。$`K=O(n)`$ または $`U(n)`$ とし、$`T`$ は
$`K`$ の対角元のなす群（符号、または位相）、$`\mathfrak k=\{X\in M_n(\mathbb F):X^*=-X\}`$、
$`\mathfrak t`$ はその対角部分、$`\mathfrak m`$ はその非対角部分（対角成分が零）であり、
座標 $`\mathrm{Re}\,X_{ij}`$、$`\mathrm{Im}\,X_{ij}`$（$`i<j`$）をもつ。$`\lambda\in\mathbb R^n`$ に対し
$`\Lambda=\mathrm{diag}(\lambda)`$、$`\Delta(\lambda)=\prod_{i<j}(\lambda_j-\lambda_i)`$ とする。$`C=\{\lambda:\lambda_1<\dots<\lambda_n\}`$ とおくと、
その上で $`\Delta>0`$ である。$`H\in\mathcal H`$ に対し、$`\lambda(H)`$ はその固有値を昇順に並べたベクトルである。
$`\|X\|_2^2=\mathrm{Tr}\,X^*X`$ とする。$`V\in K`$ に対し写像 $`H\mapsto VHV^*`$ は $`\mathcal H`$ 上で線形であり、
$`\|\cdot\|_2`$ を保つ。線形等長写像は $`|\det|=1`$ をみたすので、これは **$`dH`$ を保つ**。

*計算：* この付録の数値的検証 — [チュートリアル、付録 C](../../tutorial/ja/C_eigenvalue_densities.md)。

<a id="s-C-1"></a>

## C.1 単純スペクトルと固有フレーム

<a id="r-C-1"></a>

**補題 C.1.** (a) 相異なる $`n`$ 個の固有値をもつ行列の集合 $`\Omega\subseteq\mathcal H`$ は
開集合であり、$`\mathcal H\setminus\Omega`$ はルベーグ零集合である。(b) $`\Omega`$ 上で写像 $`H\mapsto\lambda(H)`$ は
$`C^\infty`$ であり、スペクトル射影 $`P_j(H)=\prod_{k\neq j}\frac{H-\lambda_k(H)}{\lambda_j(H)-\lambda_k(H)}`$ は連続である。

*証明.* (a) $`\lambda=\lambda(H)`$ として $`A_{ij}=\lambda_j^{i-1}`$ とおく。すると $`\det A=\Delta(\lambda)`$（ヴァンデルモンド）であり、
$`(AA^T)_{ij}=\sum_k\lambda_k^{i+j-2}=\mathrm{Tr}\,H^{i+j-2}`$ であるから、
$`\mathrm{disc}(H)=\det\big(\mathrm{Tr}\,H^{i+j-2}\big)_{i,j=1}^n=\Delta(\lambda(H))^2`$ は $`H`$ の座標の多項式で、
$`\mathrm{diag}(1,\dots,n)`$ において零でなく、$`\Omega=\{\mathrm{disc}\,\neq0\}`$ は開集合である。
$`\mathbb R^N`$ 上の零でない多項式 $`p`$ の零点集合が零集合であることは $`N`$ に関する帰納法で示される：$`N=1`$ のときそれは
有限集合である。そうでなければ、ある $`p_{j_0}\neq0`$ をもつ $`p=\sum_jp_j(x')x_N^j`$ と書く。$`p_{j_0}`$ の零点集合（帰納法の仮定により零集合）
の外の $`x'`$ に対し $`x_N\mapsto p(x',x_N)`$ の零点は有限個であり、零点集合の指示関数にトネリの定理を
適用すると測度 $`0`$ を得る。
(b) $`\chi(H,t)=\det(t-H)=\prod_i(t-\lambda_i)`$ は $`(H,t)`$ の実多項式である。$`H_0\in\Omega`$ において
$`\partial_t\chi(H_0,\lambda_j(H_0))=\prod_{k\neq j}(\lambda_j-\lambda_k)\neq0`$ であるから、陰関数定理により $`H_0`$ の近くで
$`\chi(H,r_j(H))=0`$ かつ $`r_j(H_0)=\lambda_j(H_0)`$ をみたす $`C^\infty`$ 関数 $`r_j`$ が得られる。連続性により $`H_0`$ の近くで
$`r_1<\dots<r_n`$ である：これらは $`n`$ 次多項式の相異なる $`n`$ 個の根であるから、順序も含めてそのすべての根であり、
$`r_j=\lambda_j`$ である。$`H=\sum_k\lambda_kQ_k`$ をスペクトル分解とすると、ラグランジュ補間により
$`P_j(H)=Q_j`$ である。連続性は $`\lambda`$ の連続性から従う。$`\square`$

$`H\in\Omega`$ に対し、**固有フレーム** とは $`H=V\Lambda V^*`$（$`\Lambda=\mathrm{diag}\,\lambda(H)`$）をみたす $`V\in K`$ のことである。$`V,V'`$
が固有フレームならば、$`M=V'^*V`$ は $`\Lambda`$ と可換であるから $`(\lambda_i-\lambda_j)M_{ij}=0`$ となり $`M\in T`$ である。逆に、任意の $`D\in T`$ に対し $`VD`$
は固有フレームである。**ボレル固有フレーム** $`F(H)`$ は、第 $`j`$ 列として
$`P_j(H)e_m/\|P_j(H)e_m\|`$ をとることで得られる。ここで $`m`$ は $`P_j(H)e_m\neq0`$ となる最小の添字である。これは、これらの添字が固定される
各ボレル集合の上で連続であり、その列は正規直交な固有ベクトルであり、$`H`$ が実ならば実である。

<a id="r-C-2"></a>

**補題 C.2（ハール測度）.** $`U(n)`$、$`O(n)`$、$`T`$ はハール測度をもち、すべての $`W\in K`$ と有界ボレル関数 $`f`$ に対し
$`\int_Kf(VW)\,dV=\int_Kf\,dV`$ である。

*証明.* $`U(n)`$：定理 5.3。$`O(n)`$：定理 B.17（実ガウス QR による構成：実ガウス行列はほとんど確実に
可逆であり、その QR 分解の直交因子はハール分布に従う）。$`T`$：$`\{\pm1\}^n`$ 上の一様測度、または
各円周上の正規化された弧長の積は、固定された元による乗法のもとで不変である。右不変性は補題 5.1 である。$`\square`$

<a id="r-C-3"></a>

**補題 C.3（不変アンサンブルの固有フレーム）.** $`H`$ を $`\mathcal H`$ のランダムな元で
$`P(H\in\Omega)=1`$ をみたし、その分布がすべての $`V\in K`$ に対し $`H\mapsto VHV^*`$ のもとで不変であるものとする。$`D`$ を $`T`$ 上一様で
$`H`$ と独立とし、$`\tilde U=F(H)D`$ とおく。このとき $`\tilde U`$ は $`H`$ の固有フレームであり、$`K`$ 上でハール分布に従い、
$`\lambda(H)`$ と独立である。

*証明.* $`V\in K`$ を固定する。$`F(VHV^*)`$ と $`VF(H)`$ はともに $`VHV^*`$ の固有フレームであるから、
$`F(VHV^*)=VF(H)E(H)`$ であり、$`E(H)=(VF(H))^*F(VHV^*)\in T`$ は $`H`$ のボレル関数である。$`\mathbb R^n`$ 上の有界ボレル関数
$`\varphi`$ と $`K`$ 上の有界ボレル関数 $`\psi`$ に対し、$`VHV^*`$ は $`H`$ と同じ分布をもち $`\lambda(VHV^*)=\lambda(H)`$ であるから、
```math
\mathbb E\big[\varphi(\lambda(H))\psi(\tilde U)\big]=\mathbb E\big[\varphi(\lambda(H))\psi(VF(H)E(H)D)\big]
=\mathbb E\big[\varphi(\lambda(H))\psi(V\tilde U)\big],
```
ここで2番目の等号では先に $`D`$ について積分する（フビニ）：固定された $`H`$ に対し $`E(H)D`$ は $`T`$ 上一様である。次に $`V`$
をハール分布に従い他のすべてと独立なものとする。フビニの定理と補題 C.2 により
$`\mathbb E[\varphi(\lambda(H))\psi(\tilde U)]=\mathbb E\big[\varphi(\lambda(H))\int_K\psi(V\tilde U)\,dV\big] =\mathbb E[\varphi(\lambda(H))]\int_K\psi\,dV`$ である。$`\varphi=1`$ とすると $`\tilde U`$ はハール分布に従う。積の形が
独立性である。$`\square`$

<a id="s-C-2"></a>

## C.2 ワイル積分公式

$`\varepsilon>0`$ に対し $`\mathfrak m_\varepsilon=\{X\in\mathfrak m:\|X\|_2<\varepsilon\}`$ とおき、
```math
\Psi(X,\lambda)=e^X\Lambda e^{-X},\qquad (X,\lambda)\in\mathfrak m_\varepsilon\times C .
```
とする。$`X\in\mathfrak k`$ に対し $`(e^X)^*=e^{-X}=(e^X)^{-1}`$ であるから $`e^X\in K`$ である。

<a id="r-C-4"></a>

**補題 C.4（座標近傍）.** $`\varepsilon>0`$ と $`j(0)=1`$ をみたす連続関数 $`j:\overline{\mathfrak m_\varepsilon}\to(0,\infty)`$
が存在して、次が成り立つ：(a) $`|\det D\Psi(X,\lambda)|=\Delta(\lambda)^\beta j(X)`$；(b) $`\Psi`$ は単射である；(c)
$`H\in\Omega`$ に対し、$`H\in\Psi(\mathfrak m_\varepsilon\times C)`$ であるための必要十分条件は $`F(H)\in S=\{e^XD:X\in\mathfrak m_\varepsilon,D\in T\}`$ であり、
$`S`$ は $`ST=S`$ をみたすボレル集合である。

*証明.* (a) $`\Psi`$ は $`C^\infty`$ である（べき級数）。$`Y\in\mathfrak m`$ に対し $`A_X(Y)=e^{-X}\frac{d}{dt}\big|_0e^{X+tY}`$ とおく。
$`W(t)=e^{X+tY}\in K`$ に対し $`W(t)^*W(t)=1`$ を微分すると $`A_X(Y)\in\mathfrak k`$ を得る。また
$`\frac{d}{dt}W^{-1}=-W^{-1}W'W^{-1}`$ を用いると、$`g=e^X`$ に対し
```math
D\Psi(X,\lambda)(Y,\mu)=g\big([A_X(Y),\Lambda]+\mathrm{diag}(\mu)\big)g^{-1}.
```
対角行列は $`\Lambda`$ と可換であるから、$`B_X=\pi\circ A_X:\mathfrak m\to\mathfrak m`$ として $`[A_X(Y),\Lambda]=[B_XY,\Lambda]`$ である。ここで
$`\pi`$ は $`\mathfrak t`$ に沿った $`\mathfrak k`$ から $`\mathfrak m`$ への射影である。したがって $`L_\lambda(Z,\mu)=[Z,\Lambda]+\mathrm{diag}(\mu)`$ として
$`D\Psi=\mathrm{Ad}_g\,\circ L_\lambda\circ(B_X\times1)`$ である。$`([Z,\Lambda])_{ij}=Z_{ij}(\lambda_j-\lambda_i)`$ であるから、座標において
$`L_\lambda`$ は $`\mathrm{Re}\,Z_{ij}`$（および $`\mathrm{Im}\,Z_{ij}`$）に $`\lambda_j-\lambda_i`$ を掛け、$`\mu`$ を対角線に写す：
$`|\det L_\lambda|=\Delta(\lambda)^\beta`$。$`\mathrm{Ad}_g`$ は $`dH`$ を保つ。$`j(X)=|\det B_X|`$ とおく：これは連続であり、
$`\frac d{dt}e^{tY}|_0=Y`$ であるから $`B_0=1`$ である。$`\overline{\mathfrak m_\varepsilon}`$ 上で $`j>0`$ となるように $`\varepsilon`$ を選ぶ。
(b) $`\Psi(X,\lambda)=\Psi(X',\lambda')=H`$ ならば $`\lambda=\lambda'=\lambda(H)`$ であり、$`e^X,e^{X'}`$ は $`H\in\Omega`$ の固有フレームであるから、
$`D\in T`$ を用いて $`e^X=e^{X'}D`$ と書ける。$`\mathfrak m\times\mathfrak t`$ 上で $`G(X,t)=e^Xe^t`$ とおく（$`\beta=1`$ のとき $`\mathfrak t=\{0\}`$）。$`DG(0,0)(Y,s)=Y+s`$ は単射であり、
$`Y`$ は非対角、$`s`$ は対角であるから $`\|Y+s\|_2^2=\|Y\|_2^2+\|s\|_2^2`$ となり、積ノルムに関して $`\|DG(0,0)v\|=\|v\|`$ である。$`m=1`$ とおく。$`\|z\|<\delta`$ に対し $`\|DG(z)-DG(0)\|\le m/2`$ となる $`\delta>0`$ を選ぶ。この球の中の
$`u,v`$ に対し、$`G(u)-G(v)=\int_0^1DG(v+s(u-v))(u-v)\,ds`$ から $`\|G(u)-G(v)\|\ge\frac m2\|u-v\|`$ を得るので、$`G`$ はそこで単射である。
（$`\delta`$ は $`n,\beta`$ のみに依存し、$`\varepsilon`$ はその後で選ぶ。）
さて、べき級数と
$`\|AB\|_2\le\|A\|_2\|B\|_2`$ により $`\|D-1\|_2=\|e^{X'}(D-1)\|_2=\|e^X-e^{X'}\|_2\le\|e^X-1\|_2+\|e^{X'}-1\|_2\le2(e^\varepsilon-1)`$ である。$`\beta=1`$ のとき $`D=\mathrm{diag}(\pm1)`$ であり、成分 $`-1`$ があれば
$`\|D-1\|_2\ge2`$ となる。$`\varepsilon<\log2`$ ならば $`D=1=e^t`$（$`t=0`$）である。$`\beta=2`$ のとき $`D=e^t`$、$`t=i\mathrm{diag}(\theta)`$、
$`\theta_j\in(-\pi,\pi]`$ であり、$`|\theta_j|\le\frac\pi2|e^{i\theta_j}-1|`$ である（$`|e^{i\theta}-1|=2|\sin\frac\theta2|`$ であり、$`[0,\pi/2]`$ 上で $`\sin x\ge2x/\pi`$ だから）。したがって
$`\|t\|_2\le\frac\pi2\|D-1\|_2\le\pi(e^\varepsilon-1)`$ である。$`\varepsilon`$ が小さければ（$`\varepsilon^2+\pi^2(e^\varepsilon-1)^2<\delta^2`$）、$`(X,0)`$ と $`(X',t)`$ は $`\delta`$-球に含まれ、
$`G(X,0)=e^X=e^{X'}D=G(X',t)`$ であるから $`X=X'`$ である。
(c) $`H=\Psi(X,\lambda')`$ ならば $`e^X`$ は固有フレームであるから $`F(H)\in e^XT\subseteq S`$ である。$`F(H)=e^XD`$ ならば
$`H=F(H)\Lambda F(H)^*=\Psi(X,\lambda(H))`$ である。$`S=G(\mathfrak m_\varepsilon\times\mathfrak t)`$（または $`e^{\mathfrak m_\varepsilon}T`$、$`T`$ は有限）は
$`\sigma`$-コンパクト集合の連続像であるから、$`\sigma`$-コンパクトでありボレル集合である。$`ST=S`$ は明らかである。$`\square`$

*計算：* $`n=2,3`$、$`\beta=1,2`$ に対する $`\Psi`$ のヤコビアン — [チュートリアル、C.2節](../../tutorial/ja/C_eigenvalue_densities.md#s-C-2)。

<a id="r-C-5"></a>

**定理 C.5（ワイル積分公式）.** 定数 $`c_{n,\beta}>0`$ が存在して、任意のボレル関数
$`f:\mathcal H\to[0,\infty]`$ に対し
```math
\int_{\mathcal H}f(H)\,dH=c_{n,\beta}\int_C\Delta(\lambda)^\beta\Big(\int_Kf(V\Lambda V^*)\,dV\Big)d\lambda .
```
が成り立つ。同値なことだが、$`\Phi(V,\lambda)=V\Lambda V^*`$ は $`K\times C`$ 上の測度 $`c_{n,\beta}\,dV\otimes\Delta(\lambda)^\beta d\lambda`$ を
$`\mathcal H`$ 上のルベーグ測度に写す。因子 $`\Delta(\lambda)^\beta`$ は座標近傍 $`\Psi`$ の $`X=0`$ におけるヤコビアンである
（補題 C.4(a)）。定数 $`c_{n,\beta}`$ はハール測度の正規化を反映するものであり、以下で明示的に必要になることは
ない。この公式から得られる密度はすべて、その全質量で正規化されるからである。

*証明.* *参照アンサンブル.* $`H_0`$ が密度 $`p_0(H)=Z_0^{-1}e^{-\mathrm{Tr}\,H^2/2}`$ をもつとする。
$`\mathrm{Tr}\,H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2`$ であるから、$`Z_0`$ はガウス積分の有限積である。$`p_0`$ と $`dH`$
は $`H\mapsto VHV^*`$ のもとで不変であり、したがって $`H_0`$ の分布も不変である。補題 C.1 により $`P(H_0\in\Omega)=1`$ である。補題 C.3 により、
$`\tilde U`$ はハール分布に従い、$`\lambda=\lambda(H_0)`$ と独立であり、$`H_0=\tilde U\Lambda\tilde U^*`$ である。
*$`\lambda`$ の分布.* $`A\subseteq C`$ をボレル集合とする。補題 C.4(c) と $`ST=S`$ により、
零事象を除いて $`\{H_0\in\Psi(\mathfrak m_\varepsilon\times A)\}=\{\lambda\in A,\ \tilde U\in S\}`$ であり、その確率は
$`h=\int_K1_S\,dV`$ として $`P(\lambda\in A)\,h`$ である。一方、$`\Psi`$ は開集合 $`\mathfrak m_\varepsilon\times C`$ 上で単射かつ $`C^1`$ で
ヤコビアンは零でなく、$`p_0(\Psi(X,\lambda))=p_0(\Lambda)`$ である。変数変換公式と補題 C.4(a) により
```math
P\big(H_0\in\Psi(\mathfrak m_\varepsilon\times A)\big)=\int_{\mathfrak m_\varepsilon}\!\int_Ap_0(\Lambda)\Delta(\lambda)^\beta j(X)\,d\lambda\,dX
=\kappa\int_Ap_0(\Lambda)\Delta(\lambda)^\beta\,d\lambda,\qquad \kappa=\int_{\mathfrak m_\varepsilon}j>0 .
```
$`A=C`$ とすると $`h=\kappa\int_Cp_0\Delta^\beta>0`$ である。ゆえに $`\lambda`$ は $`C`$ 上で密度 $`c\,p_0(\Lambda)\Delta(\lambda)^\beta`$（$`c=\kappa/h`$）をもつ。
*一般の $`f`$.* $`p_0>0`$ であるから $`\int f\,dH=\mathbb E\big[(f/p_0)(\tilde U\Lambda\tilde U^*)\big]`$ である。$`(\tilde U,\lambda)`$ の同時分布は
$`dV\otimes c\,p_0(\Lambda)\Delta^\beta d\lambda`$ であり、$`p_0(V\Lambda V^*)=p_0(\Lambda)`$ である。トネリの定理により $`c_{n,\beta}=c`$ として公式が得られる。
定数は公式自身によって決まる（$`f=p_0`$ とせよ）。$`\square`$

<a id="r-C-6"></a>

**系 C.6（不変な密度）.** $`H`$ が $`\mathcal H`$ 上の密度 $`p`$ をもち、すべての
$`V\in K`$ に対し $`p(VHV^*)=p(H)`$ であるとする。このとき $`\lambda(H)`$ は $`C`$ 上で密度 $`c_{n,\beta}\,p(\Lambda)\Delta(\lambda)^\beta`$ をもち、補題 C.3 の固有フレーム $`\tilde U`$
はハール分布に従い $`\lambda(H)`$ と独立であり、$`\mathbb R^n`$ 上の対称なボレル関数 $`g\ge0`$ に対し
```math
\mathbb E\,g(\lambda(H))=\frac{c_{n,\beta}}{n!}\int_{\mathbb R^n}g(x)\,p(\mathrm{diag}\,x)\,|\Delta(x)|^\beta\,dx .
```
が成り立つ。

*証明.* 補題 C.1 により $`P(H\in\Omega)=1`$ である。ボレル集合 $`A\subseteq C`$ に対し、定理 C.5 を $`f=1_A(\lambda(\cdot))\,p`$ に適用する。これについては
$`f(V\Lambda V^*)=1_A(\lambda)p(\Lambda)`$ である。$`H`$ の分布は不変であるから、補題 C.3 が適用できる。二つの
座標が一致する零集合を除けば、$`\mathbb R^n`$ は座標の置換による $`C`$ の $`n!`$ 個の像の非交和であり、
$`g`$、$`|\Delta|`$ および $`x\mapsto p(\mathrm{diag}\,x)`$ は置換不変である（置換行列は $`K`$ に属する）。$`\square`$

<a id="s-C-3"></a>

## C.3 ガウス・アンサンブル

<a id="r-C-7"></a>

**定理 C.7（GOE と GUE）.** $`H=(G+G^T)/2`$（$`G`$ は成分が独立な $`N(0,1)`$ の実行列；$`\beta=1`$、GOE）、または
$`H=(G+G^*)/\sqrt2`$（$`G`$ はジニブル行列；$`\beta=2`$、GUE）とする（5.3節と同様）。このとき $`H`$ は $`\mathcal H`$ 上で密度
$`(2\pi)^{-n/2}\pi^{-\beta n(n-1)/4}e^{-\mathrm{Tr}\,H^2/2}`$ をもち、その固有値は密度
```math
Z_{n,\beta}^{-1}\,\Delta(\lambda)^\beta\,e^{-\sum_i\lambda_i^2/2}\quad\text{on }C,\qquad Z_{n,\beta}=\int_C\Delta^\beta e^{-|\lambda|^2/2}\,d\lambda\in(0,\infty),
```
をもち、固有フレーム $`\tilde U`$ は $`K`$ 上でハール分布に従い固有値と独立である。

*証明.* GOE：$`h_{ii}=G_{ii}\sim N(0,1)`$ であり、$`i<j`$ に対し $`h_{ij}=(G_{ij}+G_{ji})/2\sim N(0,\frac12)`$ である。GUE：$`h_{ii}=\sqrt2\mathrm{Re}\,G_{ii}\sim N(0,1)`$ であり、
$`\mathrm{Re}\,h_{ij}=(\mathrm{Re}\,G_{ij}+\mathrm{Re}\,G_{ji})/\sqrt2`$ と $`\mathrm{Im}\,h_{ij}=(\mathrm{Im}\,G_{ij}-\mathrm{Im}\,G_{ji})/\sqrt2`$ は $`N(0,\frac12)`$ に従う。
相異なる座標は独立な変数の互いに素な組の関数であるから独立であり、
密度は $`\prod_i(2\pi)^{-1/2}e^{-h_{ii}^2/2}\prod_{i<j}\pi^{-\beta/2}e^{-|h_{ij}|^2}`$ である。
$`\mathrm{Tr}\,H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2`$ であるから、これは主張の密度である。これは不変であり、系 C.6 が適用できる。定数は
正規化により定まる。$`|\lambda_j-\lambda_i|\le(1+|\lambda_i|)(1+|\lambda_j|)`$ であるから
$`|\Delta(\lambda)|\le\prod_i(1+|\lambda_i|)^{n-1}`$ となり、$`Z_{n,\beta}<\infty`$ である。$`\square`$

<a id="r-C-8"></a>

**系 C.8（準位密度）.** サイズ $`n`$ の GOE または GUE に対し、$`\mathbb R^n`$ 上で $`q(x)=(n!Z_{n,\beta})^{-1}|\Delta(x)|^\beta e^{-|x|^2/2}`$ とし、
$`\rho(x)=n\int_{\mathbb R^{n-1}}q(x,y)\,dy`$ とおく。このとき有界ボレル関数
$`\varphi`$ に対し $`\mathbb E\sum_i\varphi(\lambda_i(H))=\int\varphi\rho`$ であり、$`\int\rho=n`$ であり、$`\rho`$ は連続で $`0\le\rho(x)\le C_n(1+|x|)^{\beta(n-1)}e^{-x^2/2}`$ をみたす。特に、
独立な GOE または GUE ブロックに対して定理 5.15 の仮定が成り立つ。

*証明.* 系 C.6 を $`g(x)=\sum_i\varphi(x_i)`$ に適用し（$`\varphi\ge0`$ の場合。一般には $`\varphi^\pm=\max(\pm\varphi,0)`$ に適用する）、$`q`$ の対称性を用いると、
$`\mathbb E\sum_i\varphi(\lambda_i)=n\int\varphi(x_1)q(x)\,dx=\int\varphi\rho`$ である（フビニ）。$`|\Delta|`$ の評価により
$`q(x,y)\le(n!Z_{n,\beta})^{-1}(1+|x|)^{\beta(n-1)}e^{-x^2/2}\prod_{i}(1+|y_i|)^{\beta(n-1)}e^{-y_i^2/2}`$ であり、$`y`$ について積分すると
主張の評価を得る。また $`\sup_x(1+|x|)^{\beta(n-1)}e^{-x^2/2}<\infty`$ であるから、$`x`$ によらない可積分な優関数が得られる：連続性は
優収束定理から従う。$`\square`$

*計算：* $`n=3`$ の GOE に対する定理 C.7 と系 C.8 — [チュートリアル、C.3節](../../tutorial/ja/C_eigenvalue_densities.md#s-C-3)。

<a id="r-C-9"></a>

**命題 C.9（最小間隔）.** サイズ $`n\ge2`$ の GOE または GUE に対し $`\delta=\min_i(\lambda_{i+1}-\lambda_i)`$ とおく。
（$`n,\beta`$ に依存する）定数 $`0<c\le C`$ が存在して、$`0<\varepsilon\le1`$ に対し $`c\,\varepsilon^{\beta+1}\le P(\delta<\varepsilon)\le C\varepsilon^{\beta+1}`$ が成り立つ。

*証明.* *上からの評価.* $`P(\delta<\varepsilon)\le\sum_iP(\lambda_{i+1}-\lambda_i<\varepsilon)`$ である。$`\{0<u=x_{i+1}-x_i<\varepsilon\}`$ 上で、
対 $`(i,i+1)`$ の因子を $`u^\beta`$ で、その他の因子を上と同様に評価し、$`1+|x_{i+1}|\le2(1+|x_i|)`$ と $`e^{-x_{i+1}^2/2}\le1`$ を用いる。
$`x_{i+1}=x_i+u`$ と置換すると
```math
P(\lambda_{i+1}-\lambda_i<\varepsilon)\le\frac{2^{\beta(n-1)}}{Z_{n,\beta}}\int_0^\varepsilon u^\beta du\int(1+|x|)^{2\beta(n-1)}e^{-\frac{x^2}2}dx
\Big(\int(1+|x|)^{\beta(n-1)}e^{-\frac{x^2}2}dx\Big)^{n-2}.
```
*下からの評価.* $`R_\varepsilon=\{x:0<x_1<1,\ 0<x_2-x_1<\varepsilon,\ 3k<x_k<3k+1\ (3\le k\le n)\}\subseteq C`$ とおく。$`R_\varepsilon`$ 上では
$`\delta\le x_2-x_1`$ であり、$`(1,2)`$ 以外のすべての対の距離は $`\ge1`$ であり、$`\varepsilon\le1`$ によらない $`M`$ に対し $`|x|^2\le M`$ である。したがって
$`P(\delta<\varepsilon)\ge P(\lambda\in R_\varepsilon)\ge Z_{n,\beta}^{-1}e^{-M/2}\int_0^\varepsilon u^\beta du`$ である。$`\square`$

GOE に対しては、これは5.4節の終わりで引用したオーダー $`\varepsilon^2`$ である。ハール分布に従う $`U\in O(3)`$ に対し $`\mathbb E|U_{11}|^2=\frac13`$ である：置換行列を左から掛けると
第1列の成分が置換され、ハール測度は保たれるので、三つの $`\mathbb E|U_{i1}|^2`$ は等しく、その和は $`1`$ である。

<a id="s-C-4"></a>

## C.4 直交多項式と1点関数

<a id="r-C-10"></a>

**補題 C.10（1点関数）.** $`I\subseteq\mathbb R`$ を区間とし、$`w\ge0`$ を $`I`$ 上のボレル関数で、
正の測度の集合上で正であり、すべての $`m`$ に対し $`\int_I|x|^mw<\infty`$ をみたすものとする。$`p_0,\dots,p_{n-1}`$ を実多項式で、$`\deg p_k=k`$、
最高次係数 $`\kappa_k\neq0`$ をもち、$`L^2(w)`$ において正規直交なものとする（これらは存在する：$`1,x,\dots,x^{n-1}`$ にグラム・シュミットの方法を適用すればよい。
零でない多項式の零点は有限個であるから、これらは $`L^2(w)`$ において一次独立である。各 $`p_k`$ の符号は任意に選べる）。$`I^n`$ 上で
$`Z_n=\int_{I^n}\Delta^2\prod_iw(x_i)\,dx`$ として $`q=Z_n^{-1}\Delta(x)^2\prod_iw(x_i)`$ とおく。このとき $`Z_n=n!\prod_k\kappa_k^{-2}`$ であり、
```math
n\int_{I^{n-1}}q(x,y)\,dy=K_n(x)\,w(x),\qquad K_n(x)=\sum_{k=0}^{n-1}p_k(x)^2=v(x)^TM^{-1}v(x),
```
ここで $`v(x)=(1,x,\dots,x^{n-1})`$、$`M_{ij}=\int_Ix^{i+j}w`$、$`0\le i,j<n`$ である。

*証明.* $`\Delta(x)=\det(x_j^{i-1})_{i,j}`$ である。各行にそれより上の行の定数倍を加えると $`x_j^{i-1}`$ が
$`p_{i-1}(x_j)/\kappa_{i-1}`$ に置き換わるので、$`\Delta(x)=(\prod_k\kappa_k)^{-1}\det(p_{i-1}(x_j))`$ であり、
```math
\Delta(x)^2=\Big(\prod_k\kappa_k\Big)^{-2}\sum_{\sigma,\tau\in S_n}\mathrm{sgn}\,\sigma\mathrm{sgn}\,\tau\prod_{j=1}^np_{\sigma(j)-1}(x_j)\,p_{\tau(j)-1}(x_j).
```
$`x_2,\dots,x_n`$ について $`\prod_{j\ge2}w(x_j)`$ に関して積分すると、因子 $`j`$ から $`\delta_{\sigma(j)\tau(j)}`$ が出るので、$`\sigma=\tau`$ の項
だけが残り、結果は $`(\prod_k\kappa_k)^{-2}\sum_\sigma p_{\sigma(1)-1}(x_1)^2=(\prod_k\kappa_k)^{-2}(n-1)!\,K_n(x_1)`$ となる。
さらに $`w(x_1)`$ に関して積分すると $`Z_n=(\prod_k\kappa_k)^{-2}n!`$ を得て、割り算をすると公式を得る。最後に、
$`p_k(x)=\sum_iC_{ki}x^i`$ とすると、正規直交性は $`C`$ を可逆として $`CMC^T=1`$ と書けるから、$`M^{-1}=C^TC`$ であり、
$`v^TM^{-1}v=|Cv|^2=\sum_kp_k^2`$ である。$`\square`$

系 C.6 により、GUE の順序づけない固有値は、
$`\mathbb R`$ 上の $`w(x)=e^{-x^2/2}`$ に対する補題 C.10 の密度 $`q`$ をもつ。したがって GUE の準位密度は $`K_n(x)e^{-x^2/2}`$ である。

*計算：* $`n=3`$ に対する GUE の準位密度 — [チュートリアル、C.4節](../../tutorial/ja/C_eigenvalue_densities.md#s-C-4)。

<a id="s-C-5"></a>

## C.5 ウィシャート行列と誘導測度

$`a\le b`$ とし、$`G`$ を $`a\times b`$ のジニブル行列（成分は独立で、その実部と虚部は独立に $`N(0,\frac12)`$ に従う。
密度は $`\pi^{-ab}e^{-\mathrm{Tr}\,GG^*}`$）とし、$`W=GG^*`$ とおく。$`T_+`$ を対角成分が正の $`a\times a`$ 下三角行列の集合とする。
これは座標 $`L_{ii}`$、$`\mathrm{Re}\,L_{ij}`$、$`\mathrm{Im}\,L_{ij}`$（$`i>j`$）において $`\mathbb R^{a^2}`$ の開部分集合である。また
$`P_a`$ を正定値行列の集合とする。これは $`a\times a`$ エルミート行列の中で開集合である。

<a id="r-C-11"></a>

**補題 C.11（バートレット分解）.** ほとんど確実に $`G`$ の行 $`r_1,\dots,r_a`$ は一次独立であり、
グラム・シュミットの方法により $`L\in T_+`$ と正規直交な行をもつ $`V`$ を用いて $`G=LV`$ と書けるので、$`W=LL^*`$ である。成分 $`L_{ij}`$（$`i\ge j`$）は
独立であり、$`L_{ij}`$（$`i>j`$）は標準複素ガウス変数であり、$`L_{ii}^2`$ は $`(0,\infty)`$ 上で密度 $`x^{b-i}e^{-x}/(b-i)!`$ をもつ。
したがって $`L`$ は $`T_+`$ 上で密度 $`c_L\prod_iL_{ii}^{2(b-i)+1}e^{-\mathrm{Tr}\,LL^*}`$ をもつ。

*証明.* $`\mathcal F_{i-1}=\sigma(r_1,\dots,r_{i-1})`$ とおき、（帰納法により）$`r_1,\dots,r_{i-1}`$ はほとんど確実に一次独立で、
グラム・シュミット・ベクトル $`v_1,\dots,v_{i-1}`$ をもつと仮定する。$`u_1,\dots,u_b`$ を
$`v_1,\dots,v_{i-1},e_1,\dots,e_b`$ からグラム・シュミットの方法（従属なベクトルは除く）で得られる正規直交基底とする。$`k<i`$ に対し $`u_k=v_k`$ であり、$`u_k`$ は
$`r_1,\dots,r_{i-1}`$ のボレル関数である。$`c_k=\langle u_k,r_i\rangle`$ とおく。$`r_i`$ は $`\mathcal F_{i-1}`$ と独立で、その密度 $`\pi^{-b}e^{-\|r\|^2}`$ は
ユニタリ変換のもとで不変である（命題 A.16(a)）から、$`\mathcal F_{i-1}`$ のもとで条件づけると、ベクトル $`c`$ は標準複素ガウスベクトルである（フビニ）。さて
$`j<i`$ に対し $`L_{ij}=c_j`$ であり、$`r_i-\sum_{j<i}c_jv_j=\sum_{k\ge i}c_ku_k`$ であるから、ほとんど確実に $`L_{ii}^2=\sum_{k=i}^b|c_k|^2>0`$ である。したがって
$`\mathcal F_{i-1}`$ のもとでの $`L`$ の第 $`i`$ 行の条件付き分布は固定された分布であり、第 $`i`$ 行はそれ以前の行と独立で、その座標
$`c_1,\dots,c_{i-1}`$ と $`L_{ii}`$ は独立である。帰納法により $`L`$ のすべての行は独立である。$`N=m=b-i+1`$ として定理 A.20(b) を用いると、
$`L_{ii}^2=\|(c_i,\dots,c_b)\|^2`$ は密度 $`s^{m-1}e^{-s}/(m-1)!`$ をもつ。すると $`L_{ii}=\sqrt s`$ の密度は
$`2t^{2(b-i)+1}e^{-t^2}/(b-i)!`$ であり、$`\mathrm{Tr}\,LL^*=\sum_{i\ge j}|L_{ij}|^2`$ である。$`\square`$

*計算：* $`W`$ のコレスキー因子に関する補題 C.11 — [チュートリアル、C.5節](../../tutorial/ja/C_eigenvalue_densities.md#s-C-5)。

<a id="r-C-12"></a>

**補題 C.12（コレスキー写像）.** $`\chi(L)=LL^*`$ は $`T_+`$ から $`P_a`$ の上への $`C^\infty`$ 全単射であり、
$`|\det D\chi(L)|=2^a\prod_{i=1}^aL_{ii}^{2(a-i)+1}`$ である。

*証明.* *全射：* $`W\in P_a`$ に対し、補題 5.2 を $`W^{1/2}=QR`$ に適用する。すると $`W=R^*R`$ であり、$`L=R^*\in T_+`$ である。*単射：*
$`LL^*=MM^*`$ ならば $`M^{-1}L=M^*(L^*)^{-1}`$ は下三角かつ上三角であるから、対角成分が正の対角行列であり、さらに
ユニタリであるから $`1`$ である。*ヤコビアン：* $`i\ge j`$ に対し $`W_{ij}=L_{ij}L_{jj}+\sum_{k<j}L_{ik}\overline{L_{jk}}`$ である。$`L`$ の座標と
$`W`$ の下三角部分の座標を列ごとに並べ、第 $`j`$ 列の中では対角成分を最初に置く。すると $`W_{ij}`$ は
$`L_{ij}`$ とそれより前の座標にのみ依存するので、$`D\chi`$ は三角行列であり、その対角成分は $`2L_{jj}`$（$`\partial W_{jj}/\partial L_{jj}`$）と、
各 $`i>j`$ に対して2回現れる $`L_{jj}`$（$`W_{ij}`$ の実部と虚部）である。その積は $`\prod_j2L_{jj}\cdot L_{jj}^{2(a-j)}`$ である。$`\square`$

<a id="r-C-13"></a>

**定理 C.13（ウィシャート・アンサンブルとラゲール・アンサンブル）.** $`W=GG^*`$ は $`P_a`$ 上で密度 $`c_W\det(W)^{b-a}e^{-\mathrm{Tr}\,W}`$ をもち
（その外では零）、その固有値は密度
```math
Z_{a,b}^{-1}\,\Delta(\lambda)^2\prod_{i=1}^a\lambda_i^{\,b-a}e^{-\lambda_i}\quad\text{on }C_+=\{0<\lambda_1<\dots<\lambda_a\}.
```
をもつ。

*証明.* ボレル関数 $`g\ge0`$ に対し、補題 C.11 と変数変換 $`W=\chi(L)`$（補題 C.12）により、
$`L=\chi^{-1}(W)`$ として $`\mathbb E\,g(W)=\int_{P_a}g(W)\,c_L\prod_iL_{ii}^{2(b-i)+1}e^{-\mathrm{Tr}\,W}\big/\big(2^a\prod_iL_{ii}^{2(a-i)+1}\big)dW`$ であり、
$`\prod_iL_{ii}^{2(b-a)}=\det(W)^{b-a}`$ である。この密度は $`W\mapsto VWV^*`$（$`V\in U(a)`$）のもとで不変であり、
$`\Lambda\in P_a`$ であるための必要十分条件は $`\lambda_1>0`$ である。$`\beta=2`$ として系 C.6 を適用すればよい。$`\square`$

<a id="r-C-14"></a>

**定理 C.14（誘導測度）.** $`\psi`$ を $`\mathbb C^a\otimes\mathbb C^b`$ のハール・ランダムな純粋状態（7.5節）とし、
$`s=\mathrm{Tr}\,W`$ として $`\rho_A=\mathrm{Tr}_2|\psi\rangle\langle\psi|=W/s`$（定理 7.11 の証明）とし、
$`\mu_1<\dots<\mu_a`$ を $`\rho_A`$ の固有値とする。このとき $`s`$ と $`\mu`$ は独立であり、$`s`$ は密度
$`s^{ab-1}e^{-s}/(ab-1)!`$ をもち、$`(\mu_1,\dots,\mu_{a-1})`$ は密度
```math
Z'^{-1}_{a,b}\,\Delta(\mu)^2\prod_{i=1}^a\mu_i^{\,b-a}\quad\text{on }D_a=\{0<\mu_1<\dots<\mu_{a-1}<\mu_a\},\qquad\mu_a=1-\sum_{i<a}\mu_i ,
```
をもつ。これは物理の文献における密度「$`\propto\delta(1-\sum_i\mu_i)\,\Delta(\mu)^2\prod_i\mu_i^{b-a}`$」である。

*証明.* $`\Theta(s,\mu_1,\dots,\mu_{a-1})=s\,(\mu_1,\dots,\mu_{a-1},1-\sum_{i<a}\mu_i)`$ は $`(0,\infty)\times D_a`$ から
$`C_+`$ の上への $`C^\infty`$ 全単射であり、その逆は $`\lambda\mapsto(\sum\lambda_i,\lambda_1/\sum\lambda_i,\dots)`$ である。$`D\Theta`$ のすべての行を最後の行に加えると、最後の行は
$`(1,0,\dots,0)`$ となるので、$`|\det D\Theta|=s^{a-1}`$ である。$`\Delta`$ は $`a(a-1)/2`$ 次の斉次式であるから、$`(s,\mu)`$ の密度は
```math
Z_{a,b}^{-1}\,s^{a(a-1)+a(b-a)+a-1}e^{-s}\cdot\Delta(\mu)^2\prod_i\mu_i^{b-a}=Z_{a,b}^{-1}\,s^{ab-1}e^{-s}\cdot\Delta(\mu)^2\prod_i\mu_i^{b-a},
```
であり、積の形をしている。これにより独立性が示され、正規化の後で両方の因子が特定される。$`\square`$

*計算：* $`a=2`$ に対する定理 C.14 — [チュートリアル、C.5節](../../tutorial/ja/C_eigenvalue_densities.md#s-C-5)。*ライブラリでは：* 誘導測度、すなわち `rank=b` とした `states.random_density_matrix` — [チュートリアル、C.5節](../../tutorial/ja/C_eigenvalue_densities.md#s-C-5)。

<a id="s-C-6"></a>

## C.6 ラゲール多項式とペイジの公式

<a id="r-C-15"></a>

**補題 C.15（ラゲール多項式）.** 整数 $`\alpha,k\ge0`$ に対し
$`L_k^{(\alpha)}(x)=\sum_{j=0}^k(-1)^j\binom{k+\alpha}{k-j}\frac{x^j}{j!}`$ とおき、$`L_{-1}^{(\alpha)}=0`$ とする。このとき：

(a) $`L_k^{(\alpha)}(x)\,x^\alpha e^{-x}=\frac1{k!}\frac{d^k}{dx^k}\big(x^{k+\alpha}e^{-x}\big)`$ である。$`L_k^{(\alpha)}`$ は次数 $`k`$ で、最高次係数は $`(-1)^k/k!`$ である。

(b) $`\frac{d}{dx}L_k^{(\alpha)}=-L_{k-1}^{(\alpha+1)}`$ であり、したがって $`i\le k`$ に対し $`\frac{d^i}{dx^i}L_k^{(\alpha)}=(-1)^iL_{k-i}^{(\alpha+i)}`$ である。また $`L_k^{(\alpha)}=L_k^{(\alpha+1)}-L_{k-1}^{(\alpha+1)}`$ である。

(c) $`m<k`$ に対し $`\int_0^\infty x^mL_k^{(\alpha)}x^\alpha e^{-x}dx=0`$ であり、$`\int_0^\infty L_k^{(\alpha)}L_l^{(\alpha)}x^\alpha e^{-x}dx=\delta_{kl}\frac{(k+\alpha)!}{k!}`$ である。

その結果、補題 C.10 と定理 C.13 により、$`W`$ の固有値は準位密度
（ボレル関数 $`\varphi\ge0`$、および $`\int|\varphi|\rho_{a,b}<\infty`$ をみたすすべての $`\varphi`$ に対し $`\mathbb E\sum_i\varphi(\lambda_i(W))=\int_0^\infty\varphi\,\rho_{a,b}`$）
```math
\rho_{a,b}(x)=\sum_{k=0}^{a-1}\frac{k!}{(k+\alpha)!}L_k^{(\alpha)}(x)^2\,x^\alpha e^{-x},\qquad\alpha=b-a .
```
をもつ。

*証明.* (a) ライプニッツの公式により、$`\frac{d^k}{dx^k}(x^{k+\alpha}e^{-x})=\sum_j\binom kj(-1)^j\frac{(k+\alpha)!}{(\alpha+j)!}x^{\alpha+j}e^{-x}`$（$`e^{-x}`$ に $`j`$ 回の微分が
かかる項）であり、$`\frac1{k!}\binom kj\frac{(k+\alpha)!}{(\alpha+j)!}=\binom{k+\alpha}{k-j}\frac1{j!}`$ である。(b) 項別に微分すると
$`-\sum_{i=0}^{k-1}(-1)^i\binom{k+\alpha}{k-1-i}\frac{x^i}{i!}=-L_{k-1}^{(\alpha+1)}`$ である。2番目の恒等式はパスカルの法則
$`\binom{k+\alpha+1}{k-j}-\binom{k+\alpha}{k-1-j}=\binom{k+\alpha}{k-j}`$ である。(c) (a) を $`x^m`$ に対して $`k`$ 回部分積分する：
境界項は消える。$`\frac{d^r}{dx^r}(x^{k+\alpha}e^{-x})`$（$`r<k`$）は因子 $`x^{k+\alpha-r}`$ をもち、$`\infty`$ で減衰するからである。
結果は $`\frac{(-1)^k}{k!}\int\frac{d^kx^m}{dx^k}x^{k+\alpha}e^{-x}dx`$ であり、これは $`m<k`$ ならば $`0`$、$`m=k`$ ならば $`(-1)^k(k+\alpha)!`$ である。$`l<k`$ に対し
$`L_l^{(\alpha)}`$ は $`x^m`$（$`m<k`$）の一次結合である。$`l=k`$ に対しては最高次係数 $`(-1)^k/k!`$ のみが効く。したがって
$`p_k=\big(k!/(k+\alpha)!\big)^{1/2}L_k^{(\alpha)}`$ は $`(0,\infty)`$ 上の $`w(x)=x^\alpha e^{-x}`$ に関して正規直交であり、系 C.6 により $`W`$ の順序づけない
固有値は補題 C.10 の密度 $`q`$（$`n=a`$）をもつ。よって系 C.8 の証明と同様に、
ボレル関数 $`\varphi\ge0`$ に対し $`\mathbb E\sum_i\varphi(\lambda_i)=a\int\varphi(x_1)q(x)\,dx=\int\varphi\,K_aw=\int\varphi\,\rho_{a,b}`$ であり、線形性により
$`\int|\varphi|\rho_{a,b}<\infty`$ をみたす $`\varphi`$ に対しても成り立つ。$`\square`$

*計算：* 準位密度 $`\rho_{a,b}`$ とヒストグラムとの比較 — [チュートリアル、C.6節](../../tutorial/ja/C_eigenvalue_densities.md#s-C-6)。

<a id="r-C-16"></a>

**補題 C.16（対数積分）.** $`I_m=\int_0^\infty x^m\log x\,e^{-x}dx`$、$`H_m=\sum_{i=1}^m\frac1i`$ とおく。このとき
$`I_m=m!\,(I_0+H_m)`$ である。整数 $`\alpha\ge0`$ と $`0\le k\le m`$ に対し
```math
A^{(\alpha)}_{k,m}:=\int_0^\infty\log x\,L_k^{(\alpha)}(x)L_m^{(\alpha)}(x)\,x^\alpha e^{-x}dx=
\begin{cases}I_{k+\alpha}/k!, & m=k,\\ -\dfrac{(k+\alpha)!}{k!\,(m-k)}, & m>k.\end{cases}
```
（$`I_0`$ は定理 C.17 で打ち消し合うので、その値は必要ない。）

*証明.* $`m\ge1`$ に対し、$`(-e^{-x})'=e^{-x}`$ を用いて部分積分すると $`I_m=mI_{m-1}+(m-1)!`$、すなわち
$`I_m/m!=I_{m-1}/(m-1)!+1/m`$ を得る。$`A_{k,m}`$ については、$`g=x^{m+\alpha}e^{-x}`$ として $`L_m^{(\alpha)}x^\alpha e^{-x}=\frac1{m!}g^{(m)}`$ と書き
（補題 C.15(a)）、$`f=\log x\,L_k^{(\alpha)}`$ として $`m`$ 回部分積分する。境界項
$`f^{(i)}g^{(m-1-i)}`$ は $`\infty`$ と $`0`$ で消える。原点の近くでは $`f^{(i)}=O(|\log x|+x^{-i})`$ かつ $`g^{(m-1-i)}=O(x^{\alpha+1+i})`$ だからである。したがって
$`A_{k,m}=\frac{(-1)^m}{m!}\int_0^\infty f^{(m)}g\,dx`$ である（被積分関数は $`0`$ で $`O(x^\alpha(1+|\log x|))`$）。ここで
$`f^{(m)}=\sum_{i\le k}\binom mi(\log x)^{(m-i)}\frac{d^i}{dx^i}L_k^{(\alpha)}`$ である。$`m-i\ge1`$ に対し $`(\log x)^{(m-i)}=(-1)^{m-i-1}(m-i-1)!\,x^{i-m}`$ であり、
補題 C.15(b,c) により
```math
\int_0^\infty x^{i-m}\frac{d^iL_k^{(\alpha)}}{dx^i}\,x^{m+\alpha}e^{-x}dx=(-1)^i\int_0^\infty L_{k-i}^{(\alpha+i)}\,x^{\alpha+i}e^{-x}dx=(-1)^i(\alpha+i)!\,\delta_{ik}.
```
$`m>k`$ ならば $`i=k`$ のみが寄与する：$`A_{k,m}=\frac{(-1)^m}{m!}\binom mk(-1)^{m-k-1}(m-k-1)!(-1)^k(\alpha+k)!=-\frac{(k+\alpha)!}{k!(m-k)}`$。
$`m=k`$ ならば $`i<k`$ の項は消え、$`i=k`$ の項は $`\log x\cdot(-1)^k`$ であるから、$`A_{k,k}=\frac{(-1)^k}{k!}(-1)^kI_{k+\alpha}`$ を得る。$`\square`$

*計算：* 数値積分による積分 $`A^{(\alpha)}_{k,m}`$ — [チュートリアル、C.6節](../../tutorial/ja/C_eigenvalue_densities.md#s-C-6)。

<a id="r-C-17"></a>

**定理 C.17（ペイジの公式）.** $`a\le b`$ として $`\mathbb C^a\otimes\mathbb C^b`$ のハール・ランダムな純粋状態に対し
```math
\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}.
```

*証明.* *帰着.* 定理 C.14 の記号で $`\lambda_i=s\mu_i`$ かつ $`\sum_i\mu_i=1`$ であるから、
$`\sum_i\lambda_i\log\lambda_i=s\log s-s\,S(\rho_A)`$ である。すべての項は可積分であり（$`|x\log x|\le1+x^2`$ であり、$`\mathbb Es^2<\infty`$、
$`[0,1]`$ 上で $`0\le-x\log x\le1/e`$ であるから $`0\le S(\rho_A)\le a/e`$）、補題 C.15 により $`\mathbb E\sum_i\lambda_i\log\lambda_i=\int x\log x\,\rho_{a,b}(x)\,dx`$ であり、$`s`$ は
$`S(\rho_A)`$ と独立であるから $`\mathbb E\sum_i\lambda_i\log\lambda_i=\mathbb E[s\log s]-\mathbb E[s]\,\mathbb E\,S(\rho_A)`$ である。$`s`$ の密度から
$`\mathbb Es=ab`$ と $`\mathbb E[s\log s]=I_{ab}/(ab-1)!=ab\,(I_0+H_{ab})`$（補題 C.16）を得る。ゆえに
```math
\mathbb E\,S(\rho_A)=I_0+H_{ab}-\frac1{ab}\int_0^\infty x\log x\,\rho_{a,b}(x)\,dx .
```
*積分.* 補題 C.15(b) により $`L_k^{(\alpha)}=L_k^{(\alpha+1)}-L_{k-1}^{(\alpha+1)}`$ であるから、$`x\cdot x^\alpha=x^{\alpha+1}`$ と
$`\alpha+1`$ に対する補題 C.16 を用いて
```math
\frac{k!}{(k+\alpha)!}\int_0^\infty x\log x\,L_k^{(\alpha)}(x)^2x^\alpha e^{-x}dx=\frac{k!}{(k+\alpha)!}\Big(A^{(\alpha+1)}_{k,k}-2A^{(\alpha+1)}_{k-1,k}+A^{(\alpha+1)}_{k-1,k-1}\Big)
```
```math
=(k+\alpha+1)(I_0+H_{k+\alpha+1})+2k+k\,(I_0+H_{k+\alpha}),
```
であり、これは $`k=0`$ に対しても成り立つ。$`\sum_{k=0}^{a-1}(2k+\alpha+1)=a(a-1)+a(b-a+1)=ab`$ であるから
```math
\int_0^\infty x\log x\,\rho_{a,b}\,dx=ab\,I_0+\Sigma,\qquad \Sigma=\sum_{k=0}^{a-1}\big[(k+\alpha+1)H_{k+\alpha+1}+kH_{k+\alpha}+2k\big].
```
*和.* $`m=\alpha+1=b-a+1`$、$`j=k+m\in\{m,\dots,b\}`$ とおく。$`(j-m)H_{j-1}=(j-m)H_j-1+\frac mj`$ であるから
$`\Sigma=\sum_{j=m}^b(2j-m)H_j-a+m(H_b-H_{m-1})+a(a-1)`$ である。$`F(j)=j(j+1-m)`$ とすると $`F(j)-F(j-1)=2j-m`$、$`F(m-1)=0`$、
$`F(b)=ab`$ であり、部分和分により
```math
\sum_{j=m}^b(2j-m)H_j=abH_b-\sum_{j=m}^{b-1}\frac{F(j)}{j+1}=abH_b-\sum_{j=m}^{b-1}\Big(j-m+\frac m{j+1}\Big)=abH_b-\frac{(a-1)(a-2)}2-m(H_b-H_m).
```
$`m(H_m-H_{m-1})=1`$ であるから $`\Sigma=abH_b-\frac{(a-1)(a-2)}2+1-a+a(a-1)=abH_b+\frac{a(a-1)}2`$ である。したがって
$`\mathbb E\,S(\rho_A)=H_{ab}-H_b-\frac{a-1}{2b}`$ である。$`\square`$

*計算：* 厳密な算術による $`\Sigma`$ の恒等式、および数値積分（$`a=2`$）とモンテカルロ法によるペイジの公式 — [チュートリアル、C.6節](../../tutorial/ja/C_eigenvalue_densities.md#s-C-6)。*ライブラリでは：* エントロピー $`S(\rho_A)`$ は `DensityMatrix.von_neumann_entropy` である — [チュートリアル、3.4節](../../tutorial/ja/03_traces_and_states.md#s-3-4)。

<a id="r-C-18"></a>

**系 C.18.** $`a\le b`$ に対し
$`\log a-\frac{a-1}{ab+1}-\frac{a-1}{2b}\le\mathbb E\,S(\rho_A)\le\log a-\frac{a-1}{2b}`$ であり、特に
$`\big|\mathbb E\,S(\rho_A)-\log a+\frac a{2b}\big|\le\frac1{2b}`$ である。2量子ビットでは $`\mathbb E\,S(\rho_A)=\frac13`$ である。

*証明.* $`1/x`$ は減少するので $`\int_k^{k+1}\frac{dx}x\le\frac1k\le\int_{k-1}^k\frac{dx}x`$ である。$`b<k\le ab`$ について和をとると
$`\log\frac{ab+1}{b+1}\le H_{ab}-H_b\le\log a`$ であり、$`\log\frac{ab+1}{b+1}=\log a-\log\big(1+\frac{a-1}{ab+1}\big)\ge\log a-\frac{a-1}{ab+1}`$ である。
最後に $`\frac{a-1}{ab+1}\le\frac1b`$ であり、$`a=b=2`$ とすると $`\frac13+\frac14-\frac14`$ を得る。$`\square`$

*計算：* $`a<12`$、$`b<40`$ に対する二つの評価 — [チュートリアル、C.6節](../../tutorial/ja/C_eigenvalue_densities.md#s-C-6)。

この公式は D. N. Page が予想し（Phys. Rev. Lett. 71, 1291, 1993）、S. K. Foong と S. Kanno
（Phys. Rev. Lett. 72, 1148, 1994）、J. Sánchez-Ruiz（Phys. Rev. E 52, 5653, 1995）および S. Sen（Phys. Rev. Lett. 77, 1, 1996）が証明した。
上の計算は Sen によるラゲール・アンサンブルへの帰着に従うが、クリストッフェル・ダルブーの公式を用いる代わりに
1点密度を項ごとに足し合わせている。誘導測度は K. Życzkowski と H.-J. Sommers
（J. Phys. A 34, 7111, 2001）によるものであり、補題 C.11 の分解は M. S. Bartlett（1933）に、定理 C.5 は H. Weyl による。

**戻る：**[第5章](05_random_matrices.md)（定理 C.7、系 C.8、命題 C.9 を用いる）、
[第7章](07_tensor_products.md)（定理 C.14、C.17 と系 C.18 を用いる）。
