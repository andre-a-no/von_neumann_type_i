[English](../A_analysis.md) · [Español](../es/A_analysis.md) · [Français](../fr/A_analysis.md) · [Deutsch](../de/A_analysis.md) · [中文](../zh/A_analysis.md) · **日本語** · [Русский](../ru/A_analysis.md)

# 付録 A. 解析学と確率論からの事実

各章では、本講義の予備知識を超えるが、それ自体は独立した定理というほどでもない、解析学と確率論のいくつかの
事実を用いる：多項式の零点集合は零集合である、エルミート行列の順序づけられた固有値は行列に連続的に依存する、
ガウスベクトルの非退化な線形像は密度をもつ、ガウスベクトルの方向は球面上で一様分布し、その長さと独立である、
ガウス分布のモーメントはペアリングにわたる和である、そして一様な点の間隔は正規化された指数変数である。
本付録はこれらを証明し、各章が番号で引用できるようにする。

**用いる予備知識.** 線形代数；実解析（極限、連続性、コンパクト性と $\mathbb R^n$ におけるボルツァーノ・ワイエルシュトラスの
定理、一変数のべき級数、連続関数のリーマン積分、変数変換公式）；測度と積分（ルベーグ測度、
直積測度とフビニ・トネリの定理、単調収束定理と優収束定理、押し出し測度、
密度、独立性）；初等確率論（密度
$(2\pi s^2)^{-1/2}e^{-y^2/2s^2}$ とモーメント母関数 $\mathbb Ee^{aY}=e^{a^2s^2/2}$ をもつ分布 $N(0,s^2)$）。
$\mathbb R^m$ の可逆な線形写像 $L$ とボレル集合 $E$ に対し $\lambda_m(LE)=|\det L|\,\lambda_m(E)$ である。ここで
$\lambda_m$ はルベーグ測度である。これは変数変換公式の線形の場合である。

全体を通じて、「可測」はボレル可測を意味し、$|x|$ は $x\in\mathbb R^n$ または
$\mathbb C^n$ のユークリッドノルム、$\|A\|=\sup_{|x|=1}|Ax|$ は作用素ノルム、$\|A\|_2=(\operatorname{Tr}A^*A)^{1/2}$ は
フロベニウスノルムである。

## A.1 測度と積分

**補題 A.1（ディンキン）.** $\mu,\nu$ を可測空間 $(\Omega,\mathcal F)$ 上の有限測度で
$\mu(\Omega)=\nu(\Omega)$ をみたすものとし、$\mathcal P\subseteq\mathcal F$ は有限個の共通部分について閉じていて、
$\sigma$-加法族として $\mathcal F$ を生成するとする。すべての $P\in\mathcal P$ に対し $\mu(P)=\nu(P)$ ならば、$\mu=\nu$ である。

*証明.* $\mathcal D\subseteq\mathcal F$ が *ディンキン族* であるとは、$\Omega\in\mathcal D$ であり、
$A\subseteq B$ がともに $\mathcal D$ に属するときは常に $B\setminus A\in\mathcal D$ であり、$\mathcal D$ の任意の増大列
$A_k$ に対し $\bigcup_kA_k\in\mathcal D$ となることをいう。族 $\{A\in\mathcal F:\mu(A)=\nu(A)\}$ はディンキン族である。
有限測度に対しては $\mu(B\setminus A)=\mu(B)-\mu(A)$ であり、測度は増大する和集合に沿って連続だからである。
$d$ を $\mathcal P$ を含むすべてのディンキン族の共通部分とする。これ自身もディンキン族である。
族 $\mathcal D_1=\{A\in d:\text{すべての }P\in\mathcal P\text{ に対し }A\cap P\in d\}$ はディンキン族であり（
$A\subseteq A'$ に対し $(A'\setminus A)\cap P=(A'\cap P)\setminus(A\cap P)$）、$\mathcal P$ を含むので、
$\mathcal D_1=d$ である。すると $\mathcal D_2=\{A\in d:\text{すべての }B\in d\text{ に対し }A\cap B\in d\}$ は
$\mathcal P$ を含むディンキン族であり（$\mathcal D_1=d$ だから）、したがって $\mathcal D_2=d$ である：$d$ は有限個の共通部分について閉じている。
ゆえに $d$ は $\sigma$-加法族である：$\Omega\setminus A\in d$、$A\cup B=\Omega\setminus\big((\Omega\setminus A)\cap(\Omega\setminus B)\big)\in d$ であり、
可算個の和集合は有限個の和集合の増大する和集合である。よって
$\mathcal F=\sigma(\mathcal P)\subseteq d\subseteq\{\mu=\nu\}$ である。$\square$

**系 A.2.** 距離空間 $M$ 上の二つの有限ボレル測度 $\mu,\nu$ が、すべての有界連続関数 $f:M\to\mathbb R$ に対し
$\int f\,d\mu=\int f\,d\nu$ をみたすならば、両者は等しい。

*証明.* $f=1$ とすれば $\mu(M)=\nu(M)$ である。空でない閉集合 $K$ に対し（$K=\emptyset$ ならば示すことはない）、連続関数
$f_k(x)=\max(0,1-k\operatorname{dist}(x,K))$ は $[0,1]$ に値をとり、$K$ の指示関数に単調減少しながら収束する
（$x\notin K$ ならば $\operatorname{dist}(x,K)>0$）。優収束定理により $\mu(K)=\nu(K)$ である。閉集合の全体は
有限個の共通部分について閉じていて、ボレル集合族を生成する。補題 A.1 を適用すればよい。$\square$

補題 5.1 はまさにこの議論を用いている。

**補題 A.3（ベクトル空間上のルベーグ測度）.** $W$ を $m$ 次元の実ベクトル空間とする。
線形全単射 $\varphi:\mathbb R^m\to W$ に対し、ボレル集合 $E\subseteq W$ について $\lambda_\varphi(E)=\lambda_m(\varphi^{-1}E)$ とおく。

(a) 別の線形全単射 $\psi$ に対し $\lambda_\psi=|\det(\psi^{-1}\varphi)|\,\lambda_\varphi$ である。したがって $W$ の零集合も、
$W$ 内のランダムベクトルが密度をもつかどうかも、$\varphi$ の選び方によらない。

(b) すべての線形全単射 $L:W\to W$ に対し $\lambda_\varphi(LE)=|\det L|\,\lambda_\varphi(E)$ である。$L$ が $W$ 上の
ある内積について直交ならば、$|\det L|=1$ であり、$L$ は $\lambda_\varphi$ を保つ。

(c) $W$ 内のランダムベクトル $X$ が $\lambda_\varphi$ に関して密度 $f$ をもつならば、$LX$ は密度
$y\mapsto f(L^{-1}y)/|\det L|$ をもつ。

*証明.* (a) $\lambda_\psi(E)=\lambda_m\big((\psi^{-1}\varphi)\varphi^{-1}E\big)=|\det(\psi^{-1}\varphi)|\,\lambda_m(\varphi^{-1}E)$。
(b) $\lambda_\varphi(LE)=\lambda_m\big((\varphi^{-1}L\varphi)\varphi^{-1}E\big)=|\det(\varphi^{-1}L\varphi)|\,\lambda_\varphi(E)$ であり、
$\det(\varphi^{-1}L\varphi)=\det L$ は基底によらない。正規直交基底において直交な $L$ の行列 $O$ は
$O^TO=1$ をみたすので、$(\det L)^2=1$ である。(c) (b) により、$L$ による $\lambda_\varphi$ の押し出しは
$|\det L|^{-1}\lambda_\varphi$ であるから、ボレル集合 $E$ に対し
$\mathbb P(LX\in E)=\int1_E(Lx)f(x)\,d\lambda_\varphi(x)=|\det L|^{-1}\int1_E(y)f(L^{-1}y)\,d\lambda_\varphi(y)$ である。$\square$

*規約.* 実ベクトル空間上の「ルベーグ測度」とは任意の $\lambda_\varphi$ のことである。(a) により、零集合と密度に関しては
この選択は決して問題にならない。$\mathbb C^n$ と $M_{a\times b}(\mathbb C)$ の上では、成分の実部と虚部を
座標として用いる。これらは実内積 $\operatorname{Re}\operatorname{Tr}(X^*Y)$ について正規直交であり、この内積について
$V,W$ をユニタリとする写像 $X\mapsto VXW$ はすべて直交である（$\operatorname{Tr}(X^*Y)$ を保つ）。
エルミート行列の空間 $\operatorname{Herm}_n$ の上では、座標
$h_{ii}$ と $\operatorname{Re}h_{ij},\operatorname{Im}h_{ij}$（$i<j$）を用いる。これらは実内積 $\operatorname{Tr}(XY)$ について
正規直交ではないが、(b) により、この内積について直交な写像 $H\mapsto UHU^*$（$U$ はユニタリ）は
やはりルベーグ測度を保つ。同様に、実対称行列の空間 $\operatorname{Sym}_n(\mathbb R)$ の上では座標
$h_{ii}$ と $h_{ij}$（$i<j$）を用い、$\operatorname{Tr}(XY)$ について直交な写像 $H\mapsto OHO^T$（$O\in O(n)$）は
ルベーグ測度を保つ。

**補題 A.4（ベクトル値関数の積分）.** $V$ を有限次元の実または複素ベクトル空間とし、$F:[a,b]\to V$ を
連続とする。$\int_a^bF(t)\,dt$ をある基底について座標ごとに（実部と虚部は別々に）定義する。このとき：

(a) 分点 $\tau_i\in[t_{i-1},t_i]$ をもつ分割 $a=t_0<\dots<t_N=b$ に対し、リーマン和
$\sum_iF(\tau_i)(t_i-t_{i-1})$ は $\max_i(t_i-t_{i-1})\to0$ のとき $\int_a^bF$ に収束する。特に積分は
基底によらない；

(b) 有限次元空間 $V'$ への任意の（実）線形写像 $L:V\to V'$ に対し $L\int_a^bF=\int_a^bLF$；

(c) $V$ 上の任意のノルムに対し $\|\int_a^bF\|\le\int_a^b\|F(t)\|\,dt$；

(d) すべての $F(t)$ が閉凸集合 $K$ に属するならば $\frac1{b-a}\int_a^bF\in K$ である。さらに $K$ が錐ならば
$\int_a^bF\in K$；

(e) $t\mapsto\int_a^tF$ は微分可能で導関数は $F(t)$ である。また $\int_a^bF=\int_a^cF+\int_c^bF$ であり、
$\int_a^bF(s+c)\,ds=\int_{a+c}^{b+c}F(s)\,ds$ である。

*証明.* $F$ の各座標は連続な実関数であるから、(a) と (e) は座標ごとに成り立つ。(b) は (a) から従う。
$L$ はリーマン和をリーマン和に写し、連続だからである。ノルムは連続である
（$|\|x\|-\|y\||\le\|x-y\|\le\sum_k|x_k-y_k|\,\|e_k\|$）から、$t\mapsto\|F(t)\|$ は連続であり、(c) は
(a) とリーマン和に対する三角不等式から従う。(d) については $N$ 等分の分割をとる：
リーマン和を $b-a$ で割ったものは凸結合 $\frac1N\sum_iF(\tau_i)\in K$ であり、$K$ は閉である。$\square$

これにより、補題 8.4(c)（CP 写像のリーマン和）と命題 8.12 のチェザロ平均
$\sigma_\tau=\frac1\tau\int_0^\tau T_s\rho\,ds$ が正当化される：(d) により $\sigma_\tau$ は状態であり、
(b)、(c)、(e) により $\|T_t\sigma_\tau-\sigma_\tau\|_1\le2t/\tau$ である。

## A.2 解析関数の零点集合

多重指数 $\alpha\in\mathbb N^n$ に対し $|\alpha|=\alpha_1+\dots+\alpha_n$、$y^\alpha=y_1^{\alpha_1}\cdots y_n^{\alpha_n}$ と書く。
$a\in\mathbb R^n$ と $0<\rho\le\infty$ に対し、$Q_\rho(a)=\{x\in\mathbb R^n:\max_i|x_i-a_i|<\rho\}$ を開立方体とする
（$Q_\infty(a)=\mathbb R^n$）。

**定義.** 関数 $f:Q_\rho(a)\to\mathbb R$ が **$Q_\rho(a)$ 上のべき級数** であるとは、
$f(x)=\sum_\alpha c_\alpha(x-a)^\alpha$ であって、すべての $0<r<\rho$ に対し $\sum_\alpha|c_\alpha|r^{|\alpha|}<\infty$ となることをいう。

多項式は $\mathbb R^n$ 上のべき級数である（$c_\alpha\neq0$ となるものは有限個）。多変数の場合、本講義が
定理 A.6 を適用するのは多項式に対してだけである。（連結開集合上で局所的にのみべき級数であるような関数に対しては、
さらに多変数の一致の定理が必要になるが、本講義ではそれを用いない。）半径
$r<\rho$ の閉立方体上では各項は $|c_\alpha|r^{|\alpha|}$ で抑えられるので、級数はそこで絶対かつ一様に収束する
（ワイエルシュトラスの M 判定法）：その和は項の順序によらず、$f$ は連続である。

**補題 A.5（一変数）.** $g(t)=\sum_kc_k(t-a)^k$ を $I=(a-\rho,a+\rho)$ 上のべき級数とする。

(a) $b\in I$ と $|t-b|<\rho_b=\rho-|b-a|$ に対し、
$d_j=\sum_{k\ge j}\binom kjc_k(b-a)^{k-j}$ として $g(t)=\sum_jd_j(t-b)^j$ であり、これは $(b-\rho_b,b+\rho_b)$ 上のべき級数である。

(b) ある $d_j\neq0$ ならば、$g$ は $b$ のある穴あき近傍で零点をもたない。

(c) $g$ の零点が $I$ のある点に集積するならば、$I$ 上で $g=0$ である。$I$ 上で $g\not\equiv0$ ならば、その零点集合は
可算であり、したがってルベーグ零集合である。$I$ 上で $g\equiv0$ ならば、すべての $c_k=0$ である。

*証明.* (a) $0\le s<\rho_b$ に対し
$\sum_k\sum_{j\le k}|c_k|\binom kjs^j|b-a|^{k-j}=\sum_k|c_k|(s+|b-a|)^k<\infty$ であるから、
$(t-a)^k=\big((t-b)+(b-a)\big)^k$ を二項展開して得られる二重級数は $|t-b|\le s$ で絶対収束し、
どちらの順序で和をとってもよい。先に $k$ について和をとると $d_j$ が得られ、$\sum_j|d_j|s^j<\infty$ である。
(b) $d_j\ne0$ となる最小の $j$ を $j_0$ とする。すると $g(t)=(t-b)^{j_0}h(t)$ であり、
$h(t)=\sum_{j\ge j_0}d_j(t-b)^{j-j_0}$ は $b$ の近くで連続（M 判定法）で $h(b)=d_{j_0}\neq0$ である。したがって
$b$ の近くで $h\neq0$ である。(c) $g$ がその近くで恒等的に零となるような $b\in I$ 全体の集合を $A$ とする。$A$ は開である。
零点が $b\in I$ に集積するならば、(b) によりすべての $d_j$ は零であり、(a) により $b\in A$ である。$b\in I$ が
$A$ の点の極限ならば、零点は $b$ に集積するので $b\in A$ である：$A$ は $I$ において閉である。$I$ は連結だから、$A$ は
空か $I$ である。$g\not\equiv0$ ならば $A=\emptyset$ であり、$I$ のどの点も零点の集積点ではない。すると
任意のコンパクト区間 $J\subset I$ は有限個の零点しか含まず（$J$ の無限部分集合はボルツァーノ・ワイエルシュトラスの定理により
$J$ 内に集積点をもつ）、$I$ はそのような $J$ の可算個の和集合である。$g\equiv0$ ならば、
$b=a$（このとき $d_j=c_j$）での (b) によりすべての $c_k=0$ である。$\square$

**定理 A.6（零点集合は零集合である）.** $f$ を立方体 $Q\subseteq\mathbb R^n$ 上のべき級数とする。$f$ が $Q$ の
ルベーグ測度正の部分集合上で零となるならば、$Q$ 上で $f=0$ である。同値なことに：$Q$ 上で $f\not\equiv0$ ならば、
$\{x\in Q:f(x)=0\}$ はルベーグ零集合である。

*証明.* $n$ に関する帰納法による。$n=1$ の場合は補題 A.5(c) である。$n\ge2$ とし、
$x'\in\mathbb R^{n-1}$ として $x=(x',t)$、$Q=Q'\times I$、$a=(a',a_n)$ と書き、項を $t-a_n$ のべきでまとめる：
$$
f(x',t)=\sum_kh_k(x')(t-a_n)^k,\qquad h_k(x')=\sum_{\alpha'}c_{(\alpha',k)}(x'-a')^{\alpha'} .
$$
これは許される：$x'\in Q'$ かつ $\max_i|x'_i-a_i|\le r<\rho$ に対し
$\sum_k\sum_{\alpha'}|c_{(\alpha',k)}||x'-a'|^{\alpha'}r^k\le\sum_\alpha|c_\alpha|r^{|\alpha|}<\infty$ である（
$|y|^{\alpha'}=\prod_i|y_i|^{\alpha'_i}$ とする）。同じ評価により、各 $h_k$ は $Q'$ 上のべき級数であり、
固定した $x'\in Q'$ に対し $g_{x'}(t)=f(x',t)$ は係数 $h_k(x')$ をもつ $I$ 上のべき級数である。

$Z=\{x\in Q:f(x)=0\}$ とする。$f$ は連続だから $Z$ は $Q$ において閉である。$\lambda_n(Z)>0$ と仮定する。トネリの
定理により、$Z_{x'}=\{t\in I:f(x',t)=0\}$ として $\lambda_n(Z)=\int_{Q'}\lambda_1(Z_{x'})\,dx'$ であるから、集合
$F=\{x'\in Q':\lambda_1(Z_{x'})>0\}$ は正の $(n-1)$ 次元測度をもつ。$x'\in F$ に対し、補題 A.5(c) により
$I$ 上で $g_{x'}\equiv0$ であり、したがってすべての $k$ に対し $h_k(x')=0$ である。よって各 $h_k$ は $F$ 上で零であり、帰納法の
仮定により $Q'$ 上で $h_k=0$ である。ゆえに $Q$ 上で $f(x',t)=\sum_kh_k(x')(t-a_n)^k=0$ である。$\square$

**系 A.7.** (a) $\mathbb R^n$ 上の恒等的に零ではない実多項式 $p$ が零となるのは、ルベーグ零集合上においてだけである。
零でない係数をもつ多項式は恒等的に零ではない。

(b) $W$ を有限次元の実ベクトル空間とし、$p:W\to\mathbb R$ を、ある基底の座標についての多項式（このときどの基底についても多項式である）で、
恒等的に零ではないものとする。$W$ 内のランダムベクトル $X$ が密度をもつならば、
$\mathbb P(p(X)=0)=0$ である。

(c) 真の線形部分空間 $U\subsetneq W$ は零集合であり、密度をもつランダムベクトルはほとんど確実に $U$ の外にある。
特に、$E\neq0$ が $\mathbb C^n$ の部分空間で、$v$ が $\mathbb C^n$ 上で密度をもつならば
（例えばガウス初期ベクトル）、$v$ はほとんど確実に $E$ において零でない成分をもつ。

(d) $X$ と $Y$ を有限次元実空間 $W_1$ と $W_2$ 内の独立なランダムベクトルとし、$Y$ は密度をもつとする。また
$N\subseteq W_1\times W_2$ を、その切り口 $N_x=\{y:(x,y)\in N\}$ がすべての $x$ に対し零集合であるようなボレル集合とする。
このとき $\mathbb P((X,Y)\in N)=0$ である。

*証明.* (a) 第一の主張は $Q=\mathbb R^n$ とした定理 A.6 である。第二の主張は $n$ に関する帰納法による：
$p(x',t)=\sum_kp_k(x')t^k$ と書く。ある $p_k$ は零でない係数をもつので、帰納法によりある $x'_0$ で $p_k(x'_0)\neq0$
となる（$n=1$ では $p_k$ は定数であり、$x'$ は存在しない）。すると $t\mapsto p(x'_0,t)$ は一変数の零でない多項式であり、その根は有限個である。
(b) 基底の変換は線形であるから $p$ は多項式のままであり、(a) と補題 A.3(a) により $\{p=0\}$ は零集合である。
すると $\mathbb P(p(X)=0)=\int_{\{p=0\}}f\,d\lambda=0$ である。(c) $U$ 上で零となる零でない線形汎関数 $\ell$ を選び、
$p=\ell$ に (b) を適用する。例については、$v\mapsto P_Ev$（$P_E$ は直交射影）は零でない実線形
写像であり、$\{P_Ev=0\}$ はその核、すなわち真の部分空間である。(d) $(X,Y)$ の分布は分布の直積であるから、
トネリの定理により $\mathbb P((X,Y)\in N)=\int\mathbb P(Y\in N_x)\,d\mathbb P_X(x)=0$ である。$\square$

本講義は系 A.7 を、ジニブル行列 $Z$ の $\det Z$（$|\det Z|^2$ は成分の実部と虚部についての多項式で、
$Z=1$ で $1$ に等しい；定理 5.3、1.6節）、命題 4.7(a) のグラム行列式、
ガウス行列の階数（第12章の演習問題 4）、およびランチョス法の初期ベクトル（第9章）に適用する。
(d) は命題 4.7(b) における「$Y$ は $K$ と独立であるから、トネリの定理により確率は $0$ である」という段階である。

```python
import math
import numpy as np

rng = np.random.default_rng(0)


def mc(samples, exact, label, kind="exact"):
    """Monte Carlo mean with standard error; the reference value must lie within 5 standard errors."""
    samples = np.asarray(samples, dtype=float)
    mean, err = samples.mean(), samples.std() / math.sqrt(samples.size)
    print(f"{label}: {mean:.4f} +- {err:.4f}   ({kind} {exact:.4f})")
    assert abs(mean - exact) <= 5 * err + 1e-12, label


# Corollary A.7: det Z = 0 has probability 0 for a real Gaussian 3 x 3 matrix, so by continuity of the measure
# P(|det Z| < eps) -> 0 as eps -> 0. Numerical observation (not a theorem of this appendix): the fractions
# below shrink by about a factor 10 per factor 10 in eps.
d = np.abs(np.linalg.det(rng.standard_normal((200_000, 3, 3))))
fractions = [(d < eps).mean() for eps in (1e-1, 1e-2, 1e-3)]
print("P(|det Z| < eps) for eps = 0.1, 0.01, 0.001:", fractions)
assert fractions[1] < 0.2 * fractions[0] and fractions[2] < 0.2 * fractions[1]
```

## A.3 エルミート行列の固有値

$A\in\operatorname{Herm}_n$ に対し、$\lambda_1(A)\le\dots\le\lambda_n(A)$ を重複度込みの固有値、
$v_1,\dots,v_n$ を固有ベクトルからなる正規直交基底とする（定理 1.3、系 1.4）。単位ベクトル $x$ に対し
$\langle x,Ax\rangle=\sum_i\lambda_i|\langle v_i,x\rangle|^2$ かつ $\sum_i|\langle v_i,x\rangle|^2=1$ であるから、
$\lambda_1(A)\le\langle x,Ax\rangle\le\lambda_n(A)$ である。また $\|A\|=\max_i|\lambda_i(A)|$ である（系 1.4）。

**定理 A.8（クーラント・フィッシャー）.** $A\in\operatorname{Herm}_n$ と $1\le j\le n$ に対し
$$
\lambda_j(A)=\min_{\dim U=j}\ \max_{x\in U,\,|x|=1}\langle x,Ax\rangle=\max_{\dim U=n-j+1}\ \min_{x\in U,\,|x|=1}\langle x,Ax\rangle ,
$$
ここで $U$ は与えられた次元の部分空間 $U\subseteq\mathbb C^n$ を動く。（内側の最大値と最小値は $U$ の単位球面が
コンパクトであるから存在する。外側の最大値と最小値は達成される。）

*証明.* $\dim U=j$ とし、$W=\operatorname{span}(v_j,\dots,v_n)$ とする。$\dim U+\dim W=n+1>n$ であるから、単位ベクトル
$x\in U\cap W$ が存在し、$\langle x,Ax\rangle=\sum_{i\ge j}\lambda_i|\langle v_i,x\rangle|^2\ge\lambda_j$ である。よって $U$ 上の最大値は
$\ge\lambda_j$ である。$U=\operatorname{span}(v_1,\dots,v_j)$ に対しては $\langle x,Ax\rangle=\sum_{i\le j}\lambda_i|\langle v_i,x\rangle|^2\le\lambda_j$ であり、
$x=v_j$ で等号が成り立つ。よって $U$ にわたる最小値は達成され、$\lambda_j$ に等しい。第二の公式については、
固有値が $\lambda_j(-A)=-\lambda_{n+1-j}(A)$ である $-A$ に第一の公式を適用する：
$-\lambda_j(A)=\lambda_{n+1-j}(-A)=\min_{\dim U=n+1-j}\max_x(-\langle x,Ax\rangle)=-\max_{\dim U=n+1-j}\min_x\langle x,Ax\rangle$。$\square$

**系 A.9（圧縮）.** $Q\in M_{n\times m}(\mathbb C)$ は $Q^*Q=1_m$ をみたし、$A\in\operatorname{Herm}_n$ とする。このとき
$\lambda_j(Q^*AQ)=\min\max\langle x,Ax\rangle$ であり（$j$ 次元部分空間 $U\subseteq\operatorname{ran}Q$ と単位ベクトル $x\in U$ にわたる）、
$$
\lambda_j(A)\le\lambda_j(Q^*AQ)\le\lambda_{j+n-m}(A)\qquad(1\le j\le m).
$$

*証明.* $U'\mapsto QU'$ は $\mathbb C^m$ の $j$ 次元部分空間から $\operatorname{ran}Q$ の $j$ 次元部分空間への全単射であり、
$|Qy|=|y|$、$\langle y,Q^*AQy\rangle=\langle Qy,AQy\rangle$ である。したがって $Q^*AQ$ に対する定理 A.8 がこの
公式を与える。より大きな族であるすべての $j$ 次元部分空間 $U\subseteq\mathbb C^n$ にわたって最小化すると $\lambda_j(A)$ が得られ、
これはしたがってより小さい。同様に、$\lambda_j(Q^*AQ)$ の max–min 公式は $(m-j+1)$ 次元の
$U\subseteq\operatorname{ran}Q$ を動き、その次元のすべての $U\subseteq\mathbb C^n$ にわたると、より大きな数
$\lambda_{n-(m-j+1)+1}(A)=\lambda_{j+n-m}(A)$ を与える。$\square$

**定理 A.10（ワイル）.** $A,E\in\operatorname{Herm}_n$ とすべての $j$ に対し
$$
\lambda_j(A)+\lambda_1(E)\le\lambda_j(A+E)\le\lambda_j(A)+\lambda_n(E),
\qquad\text{hence}\qquad
|\lambda_j(A+E)-\lambda_j(A)|\le\|E\|\le\|E\|_2 .
$$

*証明.* 単位ベクトル $x$ に対し $\langle x,(A+E)x\rangle\le\langle x,Ax\rangle+\lambda_n(E)$ である。$j$ 次元部分空間 $U$ の
単位ベクトルにわたる最大値をとり、次に $U$ にわたる最小値をとる（定理 A.8）：$\lambda_j(A+E)\le\lambda_j(A)+\lambda_n(E)$。これを
$A+E$ と $-E$ に適用すると $\lambda_j(A)\le\lambda_j(A+E)+\lambda_n(-E)=\lambda_j(A+E)-\lambda_1(E)$ を得る。ゆえに
$|\lambda_j(A+E)-\lambda_j(A)|\le\max(|\lambda_1(E)|,|\lambda_n(E)|)=\|E\|$ であり、$\|E\|^2=\max_i\lambda_i(E)^2\le\sum_i\lambda_i(E)^2=\|E\|_2^2$ である。$\square$

**系 A.11.** 各 $\lambda_j:\operatorname{Herm}_n\to\mathbb R$ は $\|\cdot\|$ について $1$-リプシッツであり、したがって連続かつ
ボレル可測である。ランダムなエルミート行列 $H$ に対し、$(\lambda_1(H),\dots,\lambda_n(H))$ はランダムベクトルであり、その間隔と
間隔比も（それらの分母が零でないボレル集合上で）ランダムベクトルである。独立なランダムエルミート行列の固有値ベクトルは独立である。
$n$ 個の相異なる固有値をもつ行列の全体は開集合をなす。

*証明.* $E=B-A$ として定理 A.10 を用いる。独立な確率変数の関数は独立である。
すべての $j$ に対し $\lambda_j(A)<\lambda_{j+1}(A)$ ならば、$\|B-A\|<\frac12\min_j(\lambda_{j+1}(A)-\lambda_j(A))$ に対して狭義の不等式は保たれる。$\square$

**補題 A.12（ラグランジュの公式）.** $A\in\operatorname{Herm}_n$ が相異なる固有値 $\mu_1<\dots<\mu_r$ と
スペクトル射影 $P_1,\dots,P_r$（固有空間への直交射影。したがって $A=\sum_i\mu_iP_i$、
$P_iP_l=\delta_{il}P_i$、$\sum_iP_i=1$）をもつとする。このとき任意の多項式 $p$ に対し $p(A)=\sum_ip(\mu_i)P_i$ であり、
$$
P_i=\prod_{l\neq i}\frac{A-\mu_l}{\mu_i-\mu_l} .
$$

*証明.* $k$ に関する帰納法により（$P_iP_l=\delta_{il}P_i$ を用いて）$A^k=\sum_i\mu_i^kP_i$ であり、線形性により $p(A)=\sum_ip(\mu_i)P_i$ である。
多項式 $\ell_i(t)=\prod_{l\neq i}(t-\mu_l)/(\mu_i-\mu_l)$ は $\ell_i(\mu_l)=\delta_{il}$ をみたすので、$\ell_i(A)=P_i$ である。$\square$

**命題 A.13（固有値とスペクトル射影の連続性）.** $S$ を有限次元実ベクトル空間の部分集合とし、
$A:S\to\operatorname{Herm}_n$ を連続とする。すべての $x\in S$ に対し $A(x)$ がちょうど $r$ 個の相異なる固有値
$\mu_1(x)<\dots<\mu_r(x)$ をもつと仮定する。このときそれらの重複度は $S$ 上で局所定数であり、$\mu_i$ と
$A(x)$ のスペクトル射影 $P_i(x)$ は $S$ 上で連続である。特に、単純スペクトルをもつエルミート行列の
開集合上で、固有値と階数1のスペクトル射影は連続である。

*証明.* $r=1$ ならば、$\mu_1=\operatorname{Tr}A/n$、$P_1=1$ として $A(x)=\mu_1(x)1$ である。$r\ge2$ とし、$x_0\in S$ を固定して
$d=\min_i(\mu_{i+1}(x_0)-\mu_i(x_0))>0$ とおく。$\lambda_j(A(x_0))=\mu_i(x_0)$ となる添字 $j$ の集合を $J_i$ とする。これは
連続する整数のブロックである。$A$ の連続性により、$x\in N$ に対し $\|A(x)-A(x_0)\|<d/2$ となるような $x_0$ の $S$ における近傍 $N$ が存在する。
定理 A.10 により、$j\in J_i$ に対し $\lambda_j(A(x))$ は区間 $I_i=(\mu_i(x_0)-d/2,\mu_i(x_0)+d/2)$ にあり、
これら $r$ 個の区間は互いに素である。よって $A(x)$ は各 $I_i$ に少なくとも一つの固有値をもち、相異なる固有値をちょうど $r$ 個もつので、
ちょうど一つをもつ：$\lambda_j(A(x))$ はすべての $j\in J_i$ について同じである。ゆえに $N$ 上で重複度は定数 $|J_i|$ であり、
$\mu_i(x)=\lambda_{\max J_i}(A(x))$ は連続である（系 A.11）。補題 A.12 により
$P_i(x)=\prod_{l\neq i}(A(x)-\mu_l(x))/(\mu_i(x)-\mu_l(x))$ であり、分母は $0$ から離れているので $N$ 上で連続である。
最後の主張は $r=n$、$S$ を系 A.11 の開集合とした場合である。$\square$

命題 4.7 では、$S$ はグラム行列式 $\delta$ が零でない集合であり、その上で $K$ はちょうど $N$ 個の相異なる
固有値をもつ。したがって $K\mapsto\Pi_i(K)$ はそこで連続であり、命題 4.7(b) の事象はボレル集合である。定理 5.15
は系 A.11 を用い、定理 9.14 は系 A.9 を用いる。

```python
def herm(n, batch):
    G = rng.standard_normal((batch, n, n)) + 1j * rng.standard_normal((batch, n, n))
    return (G + G.conj().transpose(0, 2, 1)) / 2


n = 6
A, E = herm(n, 2000), 0.1 * herm(n, 2000)
lam_A, lam_AE, lam_E = np.linalg.eigvalsh(A), np.linalg.eigvalsh(A + E), np.linalg.eigvalsh(E)   # ascending
normE = np.linalg.norm(E, ord=2, axis=(1, 2))                                    # operator norm
assert (np.abs(lam_AE - lam_A).max(axis=1) <= normE + 1e-12).all()               # Theorem A.10
assert (lam_AE >= lam_A + lam_E[:, :1] - 1e-12).all() and (lam_AE <= lam_A + lam_E[:, -1:] + 1e-12).all()

a0, j = A[0], 3                          # Theorem A.8: every 3-dim subspace has max <x, Ax> >= lambda_3
for _ in range(200):
    Q, _ = np.linalg.qr(rng.standard_normal((n, j)) + 1j * rng.standard_normal((n, j)))
    assert np.linalg.eigvalsh(Q.conj().T @ a0 @ Q)[-1] >= lam_A[0, j - 1] - 1e-12
m = 4                                    # Corollary A.9: interlacing for a compression to 4 dimensions
Q, _ = np.linalg.qr(rng.standard_normal((n, m)) + 1j * rng.standard_normal((n, m)))
theta = np.linalg.eigvalsh(Q.conj().T @ a0 @ Q)
assert (theta >= lam_A[0, :m] - 1e-12).all() and (theta <= lam_A[0, n - m:] + 1e-12).all()

U, _ = np.linalg.qr(rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)))
mu, mult = np.array([-1.0, 0.5, 2.0]), [2, 3, 1]          # Lemma A.12 with multiplicities 2, 3, 1
spec = np.repeat(mu, mult)
A1 = U @ np.diag(spec) @ U.conj().T
for i, m_i in enumerate(mult):
    P = np.eye(n, dtype=complex)
    for l in range(len(mu)):
        if l != i:
            P = P @ (A1 - mu[l] * np.eye(n)) / (mu[i] - mu[l])
    Vi = U[:, spec == mu[i]]
    assert np.abs(P - Vi @ Vi.conj().T).max() < 1e-10 and abs(np.trace(P).real - m_i) < 1e-10
```

## A.4 ガウスベクトル

$\mathbb R^m$ の **標準ガウスベクトル** とは、座標が独立な $N(0,1)$ である $g=(g_1,\dots,g_m)$ のことである。その分布は
座標の分布の直積であるから、トネリの定理により密度
$\varphi_m(x)=(2\pi)^{-m/2}e^{-|x|^2/2}$ をもつ。

**命題 A.14（直交不変性）.** (a) すべての直交行列 $O\in O(m)$ に対し、$Og$ は $g$ と同じ分布をもつ。

(b) $W$ を $m$ 次元の実内積空間とする。正規直交基底 $(w_i)$ に対し、
$\sum_ig_iw_i$ の分布は基底によらない。このようなベクトルを **$W$ の標準ガウスベクトル** と呼ぶ。その座標は
任意の正規直交基底において独立な $N(0,1)$ であり、その分布は $W$ の直交写像のもとで不変であり、
$W=W_1\oplus W_2$ が直交直和ならば、その成分は $W_1$ と $W_2$ の独立な標準ガウスベクトルである。

(c) $v\in\mathbb R^m$ に対し $\langle v,g\rangle$ は分布 $N(0,|v|^2)$ をもつ。特に $\sum_ia_ig_i\sim N(0,\sum_ia_i^2)$ である。

*証明.* (a) 補題 A.3(c) により、$Og$ は密度 $\varphi_m(O^{-1}y)/|\det O|=\varphi_m(y)$ をもつ。(b) $w'_i=\sum_jO_{ji}w_j$
が別の正規直交基底ならば、$O$ は直交であり $\sum_ig_iw'_i=\sum_j(Og)_jw_j$ である。これは (a) により
$\sum_jg_jw_j$ と同じ分布をもつ。$X=\sum_ig_iw_i$ の別の正規直交基底 $(w'_i)$ における座標は $O^Tg$ であり、
(a) によりやはり標準的である。直交写像 $L$ は $(w_i)$ を正規直交基底に写し、$LX=\sum_ig_iLw_i$ である。最後の
主張については $W_1\oplus W_2$ に適合した $(w_i)$ をとればよい。(c) $v=0$ ならば両辺とも $0$ における点質量である。$v\neq0$ ならば
$v/|v|$ を正規直交基底に拡張する。(b) により座標 $\langle v/|v|,g\rangle$ は $N(0,1)$ であるから、
$\langle v,g\rangle=|v|\,\langle v/|v|,g\rangle$ は密度 $(2\pi|v|^2)^{-1/2}e^{-y^2/2|v|^2}$ をもつ（$m=1$ とした補題 A.3(c)）。$\square$

**定義.** 有限次元実ベクトル空間 $W$ 内の **中心化ガウスベクトル** とは、線形写像
$A:\mathbb R^M\to W$ と $\mathbb R^M$ の標準ガウスベクトル $g$ により $X=Ag$ と表されるもののことである。座標（$W=\mathbb R^m$）では、その
**共分散** は $C=\mathbb E[XX^T]=AA^T$ である。$A$ が全射であるとき、同値なことに $C$ が正定値であるとき
（$\langle x,Cx\rangle=|A^Tx|^2$）、それは **非退化** であるという。

本講義のランダムベクトルのうち、独立なガウス成分から線形に構成されるもの――複素ガウス
ベクトル、ジニブル行列、GOE 行列と GUE 行列、それらのツワールによる像――はすべて中心化ガウスベクトルである。

**定理 A.15（線形像と密度）.** $X=Ag$ を $W$ 内の中心化ガウスベクトルとする。

(a) すべての線形写像 $L:W\to W'$ に対し、$LX$ は中心化ガウスベクトルである。$X$ が非退化で $L$ が全射ならば、$LX$ は非退化である。

(b) $X$ が非退化ならば、$W$ 上のルベーグ測度に関して密度をもつ。座標
$W=\mathbb R^m$ では
$$
f_X(x)=(2\pi)^{-m/2}(\det C)^{-1/2}\,e^{-\langle x,C^{-1}x\rangle/2}.
$$
$X$ が退化しているならば、それは真部分空間 $\operatorname{ran}A$ に含まれ、これは零集合であり（系 A.7(c)）、$X$ は密度をもたない。

(c) $X$ の分布はその共分散だけで決まる。

*証明.* (a) $LX=(LA)g$ であり、$A$ と $L$ が全射ならば $LA$ も全射である。(b) $K=\ker A$ とすると、
$\mathbb R^M=K^\perp\oplus K$ である。命題 A.14(b) により $g=g_1+g_2$ で、$g_1$ は
$K^\perp$ の標準ガウスベクトルであり、$X=Ag_1=Bg_1$ である。ここで $B=A|_{K^\perp}$ は単射であり、$A$ が全射であるから全単射
$K^\perp\to W$ である。したがって $\dim K^\perp=m$ である。$K^\perp$ の正規直交基底に $K$ の正規直交基底を補った基底において、$A$ の行列は
$(B\ 0)$ であるから、$C=AA^T=BB^T$ である。$g_1$ の座標は密度 $\varphi_m$ をもち、補題 A.3(c) により $X$ は
密度 $\varphi_m(B^{-1}x)/|\det B|$ をもつ。ここで $|B^{-1}x|^2=\langle x,(BB^T)^{-1}x\rangle$、$|\det B|=(\det C)^{1/2}$ である。
$X$ が退化していて密度 $f$ をもつとすると、$1=\mathbb P(X\in\operatorname{ran}A)=\int_{\operatorname{ran}A}f\,d\lambda=0$ となる。
(c) $\ker C=\ker A^T$（$AA^Tx=0$ ならば $|A^Tx|^2=0$）であるから、$\operatorname{ran}C=(\ker C)^\perp=(\ker A^T)^\perp=\operatorname{ran}A$
は $X$ を含む。$d=\operatorname{rank}C$ とし（$d=0$ ならば $X=0$）、$\operatorname{ran}C$ の正規直交基底 $p_1,\dots,p_d$ を
選び（これは $C$ だけに依存する）、行が $p_k^T$ である $d\times m$ 行列を $P$ とする。すると $PP^T=1_d$ であり、
$P^TP$ は $\operatorname{ran}C$ への直交射影である。このとき $PX=(PA)g$ であり、$PA$ は $\mathbb R^d$ への全射であるから、これは共分散
$PCP^T$ をもつ非退化な中心化ガウスベクトルである。(b) によりその分布は $C$ だけで決まり、したがって
$X=P^T(PX)$ の分布も $C$ だけで決まる。$\square$

例えば第4章では、ツワール $P$ は $M_n(\mathbb C)$ から $\mathcal A'$ への線形射影であるから、
$Z\mapsto P(Z)$ と $Z\mapsto P(Z+Z^*)$ はそれぞれ $\mathcal A'$ と $\mathcal A'_h$ への実線形全射である（
$H\in\mathcal A'_h$ に対し $P(\frac H2+\frac H2)=H$）。ジニブル行列 $Z$（命題 A.16 により非退化）に対し、定理 A.15 により $Y=P(Z)$
と $K=P(Z+Z^*)$ はこれらの空間上で密度をもつ。$\mathcal A'$ の基底 $(B_k)$ のガウス結合 $\sum_k\xi_kB_k$ は、
標準ガウスベクトル $\xi$ の線形全単射による像である。

**命題 A.16（複素ガウスベクトル、ジニブル行列）.** (a) $z\in\mathbb C^n$ の座標は独立で、各座標の実部と虚部は
独立な $N(0,\frac12)$ であるとする（**標準複素ガウスベクトル**）。
このとき $\sqrt2z$ は内積 $\operatorname{Re}\langle u,w\rangle$ をもつ $\mathbb C^n\cong\mathbb R^{2n}$ の標準ガウスベクトルである。
$z$ は密度 $\pi^{-n}e^{-|z|^2}$ をもつ。すべての $U\in U(n)$ に対し $Uz$ は $z$ と同じ分布をもち、したがって
任意の正規直交基底における座標 $\langle w_i,z\rangle$ もまた独立な標準複素ガウス変数である。
そして $\mathbb E[z_i\bar z_j]=\delta_{ij}$、$\mathbb E[z_iz_j]=0$ である。

(b) $a\times b$ ジニブル行列 $G$（成分が独立な標準複素ガウス変数）は $M_{a\times b}(\mathbb C)$ 上で密度
$\pi^{-ab}e^{-\|G\|_2^2}$ をもち、すべての $V\in U(a)$、$W\in U(b)$ に対し $VGW$ は $G$ と同じ分布をもつ。

(c) 成分が独立な $N(0,1)$ である $a\times b$ 実行列 $G$ は密度 $(2\pi)^{-ab/2}e^{-\|G\|_2^2/2}$ をもち、
すべての $O\in O(a)$、$O'\in O(b)$ に対し $OGO'$ は $G$ と同じ分布をもつ。

*証明.* (a) $2n$ 個の数 $\sqrt2\operatorname{Re}z_i,\sqrt2\operatorname{Im}z_i$ は独立な $N(0,1)$ であり、それらは
正規直交基底 $(e_i,ie_i)$ における $\sqrt2z$ の座標である。$L=2^{-1/2}$ とした補題 A.3(c) により、$z$ は密度
$(2\pi)^{-n}e^{-|\sqrt2z|^2/2}2^{n}=\pi^{-n}e^{-|z|^2}$ をもつ。ユニタリ $U$ は実線形で
$\operatorname{Re}\langle u,w\rangle$ を保つので直交である。命題 A.14(b) を適用する。基底 $(w_i)$ における座標は、
$w_i$ を列とするユニタリ $W$ により $W^*z$ である。最後に、$z_j=x_j+iy_j$ とすると
$\mathbb E[z_i\bar z_j]=\mathbb E[x_ix_j+y_iy_j]+i\,\mathbb E[y_ix_j-x_iy_j]=\delta_{ij}$、$\mathbb E[z_iz_j]=\mathbb E[x_ix_j-y_iy_j]+i\,\mathbb E[x_iy_j+y_ix_j]=0$ である。
(b) $|G|=\|G\|_2$ として $M_{a\times b}(\mathbb C)\cong\mathbb C^{ab}$ であり、$X\mapsto VXW$ は
$\operatorname{Tr}(X^*Y)$ についてユニタリである（補題 A.3 の後の規約）。(a) を適用する。(c) $\mathbb R^{ab}$ と命題 A.14 を用いて同様である。$\square$

**命題 A.17（GOE と GUE）.** $H=(G+G^T)/2$（$G$ は成分が独立な $N(0,1)$ である $n\times n$ 実行列）
（GOE）、または $H=(G+G^*)/\sqrt2$（$G$ はジニブル行列）（GUE）とする。このとき成分 $h_{ij}$（$i\le j$）は独立であり、
$h_{ii}\sim N(0,1)$ であり、$i<j$ に対して $h_{ij}\sim N(0,\frac12)$（GOE）、あるいは $\operatorname{Re}h_{ij},\operatorname{Im}h_{ij}$ が
独立な $N(0,\frac12)$（GUE）である。補題 A.3 の後で定めた座標における $H$ の密度は
$$
(2\pi)^{-n/2}\pi^{-n(n-1)/4}e^{-\operatorname{Tr}H^2/2}\ \ \text{(GOE)},\qquad
(2\pi)^{-n/2}\pi^{-n(n-1)/2}e^{-\operatorname{Tr}H^2/2}\ \ \text{(GUE)},
$$
であり、$OHO^T$（$O\in O(n)$）、あるいは $UHU^*$（$U\in U(n)$）は $H$ と同じ分布をもつ。

*証明.* GOE：$G$ は内積 $\operatorname{Tr}(X^TY)$ をもつ $M_n(\mathbb R)$ の標準ガウスベクトルであり、行列
$e_{ii}$、$(e_{ij}\pm e_{ji})/\sqrt2$（$i<j$）は正規直交基底をなす。命題 A.14(b) により、座標
$G_{ii}$、$(G_{ij}+G_{ji})/\sqrt2$、$(G_{ij}-G_{ji})/\sqrt2$ は独立な $N(0,1)$ である。そして $h_{ii}=G_{ii}$、
$h_{ij}=2^{-1/2}\cdot(G_{ij}+G_{ji})/\sqrt2$ である。GUE：$\sqrt2G$ は内積
$\operatorname{Re}\operatorname{Tr}(X^*Y)$ をもつ $M_n(\mathbb C)$ の標準ガウスベクトルであり（命題 A.16）、正規直交基底は $e_{ii}$、$ie_{ii}$、
$(e_{ij}\pm e_{ji})/\sqrt2$、$i(e_{ij}\pm e_{ji})/\sqrt2$ である。$\sqrt2G$ の座標は $\sqrt2\operatorname{Re}G_{ii}$、
$\sqrt2\operatorname{Im}G_{ii}$、$\operatorname{Re}(G_{ij}\pm G_{ji})$、$\operatorname{Im}(G_{ij}\pm G_{ji})$ であり、独立な $N(0,1)$ である。そして
$h_{ii}=\sqrt2\operatorname{Re}G_{ii}$、$\operatorname{Re}h_{ij}=2^{-1/2}\operatorname{Re}(G_{ij}+G_{ji})$、
$\operatorname{Im}h_{ij}=2^{-1/2}\operatorname{Im}(G_{ij}-G_{ji})$ である。どちらの場合も密度は一次元の密度の積、すなわち
$(2\pi)^{-1/2}e^{-h_{ii}^2/2}$ と、分散 $\frac12$ の各実座標 $t$ に対する $\pi^{-1/2}e^{-t^2}$ の積であり、
$\operatorname{Tr}H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2$ である。写像 $H\mapsto UHU^*$ は線形かつ可逆で、
$\operatorname{Tr}(XY)$ について直交であるから、補題 A.3(b)、(c) により $UHU^*$ は密度 $f(U^*HU)=f(H)$ をもつ。$f$ は
$\operatorname{Tr}H^2$ だけに依存するからである。$O$ についても同様である。$\square$

```python
Amat = rng.standard_normal((2, 4))                       # Theorem A.15: X = A g in R^2, A of rank 2
C = Amat @ Amat.T
X = rng.standard_normal((400_000, 4)) @ Amat.T
h = 1 / 400                                              # midpoint rule for the density on [0, 1]^2
grid = (np.arange(400) + 0.5) * h
xx = np.stack(np.meshgrid(grid, grid, indexing='ij'), -1).reshape(-1, 2)
dens = np.exp(-0.5 * np.einsum('ni,ij,nj->n', xx, np.linalg.inv(C), xx)) / (2 * math.pi * math.sqrt(np.linalg.det(C)))
mc(((X >= 0) & (X <= 1)).all(axis=1), dens.sum() * h * h, "P(X in [0,1]^2) against the integral of f_X",
   kind="midpoint rule")

nn, B = 4, 100_000                                       # Proposition A.17
G = rng.standard_normal((B, nn, nn))
goe = (G + G.transpose(0, 2, 1)) / 2
Gc = (rng.standard_normal((B, nn, nn)) + 1j * rng.standard_normal((B, nn, nn))) / math.sqrt(2)
gue = (Gc + Gc.conj().transpose(0, 2, 1)) / math.sqrt(2)
mc(np.einsum('bij,bji->b', goe, goe), nn * (nn + 1) / 2, "E Tr H^2, GOE")      # n + 2 * n(n-1)/2 * 1/2
mc(np.einsum('bij,bji->b', gue, gue).real, nn ** 2, "E Tr H^2, GUE")           # n + 2 * n(n-1)/2 * 1
O, _ = np.linalg.qr(rng.standard_normal((nn, nn)))
rot = O @ goe @ O.T                                                             # has the law of the GOE
mc(rot[:, 0, 0] ** 4, 3.0, "E (O H O^T)_11^4")
mc(rot[:, 0, 1] ** 2, 0.5, "E (O H O^T)_12^2")
```

## A.5 ガウスベクトルの方向

$S^{n-1}=\{x\in\mathbb R^n:|x|=1\}$、$V_n=\lambda_n(\{|x|\le1\})$ とし、ボレル集合 $B\subseteq S^{n-1}$ に対し
$C(B)=\{x\in\mathbb R^n:0<|x|\le1,\ x/|x|\in B\}$ とおく。これはボレル集合である（$x\mapsto(|x|,x/|x|)$ は $\mathbb R^n\setminus0$ 上で連続）。

**定理 A.18（極座標）.** ボレル集合 $A\subseteq(0,\infty)$ と $B\subseteq S^{n-1}$ に対し
$$
\lambda_n\big(\{x\neq0:|x|\in A,\ x/|x|\in B\}\big)=n\,\lambda_n(C(B))\int_Ar^{n-1}\,dr .
$$
したがって、ボレル関数 $h:(0,\infty)\to[0,\infty]$ と $\phi:S^{n-1}\to[0,\infty]$ に対し、$\mu(B)=\lambda_n(C(B))$ として
$$
\int_{\mathbb R^n}h(|x|)\,\phi(x/|x|)\,dx=n\int_0^\infty h(r)r^{n-1}\,dr\int_{S^{n-1}}\phi\,d\mu .
$$

*証明.* $B$ を固定し、$\nu_1(A)$ と $\nu_2(A)$ を両辺とする。これらは $(0,\infty)$ 上のボレル測度である。拡大
$x\mapsto rx$ の行列式は $r^n$ であり、$rx/|rx|=x/|x|$ であるから $C(B)$ を $\{0<|x|\le r,\ x/|x|\in B\}$ の上へ写す。よって
$\nu_1((0,r])=r^n\lambda_n(C(B))=\nu_2((0,r])$ である。$R>0$ に対し、$\nu_1,\nu_2$ の $(0,R]$ への制限は有限であり
（高々 $\lambda_n(\{|x|\le R\})$）、区間 $(0,r]$（$r\le R$）上で一致する。これらの区間は共通部分について閉じていて、
$(0,R]$ を含み、$(0,R]$ のボレル集合族を生成する（それらから $(s,r]=(0,r]\setminus(0,s]$ が得られ、$(0,R]$ の
相対開集合はすべてそのような区間の可算個の和集合である）。補題 A.1 により $(0,R]$ 上で $\nu_1=\nu_2$ である。$R\to\infty$ とする。
積分公式は $h=1_A$、$\phi=1_B$ に対して成り立ち（いま証明した等式である；$\{0\}$ は零集合）、したがって線形性により非負単関数
$h,\phi$ に対して成り立ち、単関数 $h_k\uparrow h$、$\phi_k\uparrow\phi$ に沿った単調収束定理により一般の $h,\phi$ に対して成り立つ。$\square$

**系 A.19（動径積分）.** $(0,\infty)$ 上のボレル関数 $h\ge0$ に対し、
$\int_{\mathbb R^n}h(|x|)\,dx=s_n\int_0^\infty h(r)r^{n-1}\,dr$ である。ここで $s_n=nV_n=(2\pi)^{n/2}/I_n$、
$I_n=\int_0^\infty r^{n-1}e^{-r^2/2}dr$ である。$I_1=\sqrt{\pi/2}$、$I_2=1$、$I_{n+2}=nI_n$ であるから、$s_1=2$、$s_2=2\pi$、$s_3=4\pi$、
$s_4=2\pi^2$ である。一般に $s_n=2\pi^{n/2}/\Gamma(n/2)$ である。

*証明.* この公式は $\phi=1$、$\mu(S^{n-1})=V_n$ とした定理 A.18 である。$h(r)=e^{-r^2/2}$ に対し左辺は
$\big(\int_{\mathbb R}e^{-t^2/2}dt\big)^n=(2\pi)^{n/2}$（フビニ）であり、これから $s_n$ が得られる。$I_1=\frac12\sqrt{2\pi}$、$I_2=1$ であり、
部分積分により $I_{n+2}=\big[-r^ne^{-r^2/2}\big]_0^\infty+n\int_0^\infty r^{n-1}e^{-r^2/2}dr=nI_n$ を得る。置換
$r=\sqrt{2u}$ により $I_n=2^{n/2-1}\Gamma(n/2)$ を得る。$\square$

$n=2$ と $n=3$ の場合、これらは命題 5.13 の極座標と球座標である。

**定義.** $S^{n-1}$ 上の **一様分布** とは $\sigma_{n-1}(B)=\lambda_n(C(B))/V_n$ のことである。
$\mathbb C^N\cong\mathbb R^{2N}$ の単位球面上では、それは $\sigma_{2N-1}$ である。

**定理 A.20（動径的な分布の極分解）.** $X$ を、$f(|x|)$ の形の密度をもつ $\mathbb R^n$ 内のランダムベクトルとする。
このとき $X\neq0$ はほとんど確実に成り立ち、$|X|$ は $(0,\infty)$ 上で密度 $s_nr^{n-1}f(r)$ をもち、$X/|X|$ は分布
$\sigma_{n-1}$ をもち、$|X|$ と $X/|X|$ は独立である。さらに、すべての $O\in O(n)$ に対し $\sigma_{n-1}(OB)=\sigma_{n-1}(B)$ である。
特に：

(a) $\mathbb R^n$ の標準ガウスベクトル $g$ に対し、$g/|g|\sim\sigma_{n-1}$ であり、$|g|$ と独立である；

(b) 標準複素ガウスベクトル $z\in\mathbb C^N$ に対し、$\psi=z/|z|\sim\sigma_{2N-1}$ であり、これは $U(N)$ のもとで不変な分布で、
$|z|$ と独立である。また $|z|^2$ は $(0,\infty)$ 上で密度 $t^{N-1}e^{-t}/(N-1)!$ をもち、したがって $\mathbb E|z|^4=N(N+1)$ である。

*証明.* $h=f\cdot1_A$、$\phi=1_B$ とした定理 A.18 により
$\mathbb P(|X|\in A,\ X/|X|\in B)=n\,\lambda_n(C(B))\int_Af(r)r^{n-1}dr$ である。$A=(0,\infty)$ かつ $B=S^{n-1}$ とするとこれは
$1=nV_n\int_0^\infty f(r)r^{n-1}dr$ である（そして $\mathbb P(X=0)=0$）。$A=(0,\infty)$ とすると $\mathbb P(X/|X|\in B)=\sigma_{n-1}(B)$ を、
$B=S^{n-1}$ とすると $\mathbb P(|X|\in A)=s_n\int_Af(r)r^{n-1}dr$ を得て、一般の場合は両者の積である。
$O\in O(n)$ に対し $OC(B)=C(OB)$ かつ $\lambda_n(OE)=\lambda_n(E)$ である（補題 A.3(b)）。(a) $\varphi_n(x)=(2\pi)^{-n/2}e^{-|x|^2/2}$ は動径的である。
(b) $g$ を $\mathbb R^{2N}$ の標準ガウスベクトルとして $z=g/\sqrt2$ であるから（命題 A.16）、$z/|z|=g/|g|$、$|z|=|g|/\sqrt2$ である。ユニタリは
$\mathbb R^{2N}$ 上で直交である。$|g|$ の密度は $s_{2N}(2\pi)^{-N}r^{2N-1}e^{-r^2/2}$（$s_{2N}=2\pi^N/(N-1)!$）であり、
置換 $t=r^2/2$ により $|z|^2$ の密度が得られる。すると $\mathbb E|z|^4=\int_0^\infty t^{N+1}e^{-t}dt/(N-1)!=(N+1)N$ である。$\square$

補題 7.9 は (b) である。$N=1$ のとき、$\sigma_1$ は正規化された弧長、すなわち $[0,2\pi)$ 上の
$\frac1{2\pi}\lambda_1$ の $\theta\mapsto e^{i\theta}$ による押し出しである：ボレル集合 $B\subseteq S^1$ に対し
$\Theta(B)=\{\theta\in(0,2\pi):e^{i\theta}\in B\}$ とおく。写像 $(r,\theta)\mapsto re^{i\theta}$ は開集合
$(0,1)\times(0,2\pi)$ から、開単位円板から線分 $[0,1)$ を除いたものの上への $C^1$ 全単射で、ヤコビアンは $r$ である。また $C(B)$ は
$(0,1)\times\Theta(B)$ の像と、その線分と単位円周の部分集合、すなわち零集合だけ異なる。変数変換公式により
$\lambda_2(C(B))=\int_{\Theta(B)}\int_0^1r\,dr\,d\theta=\frac12\lambda_1(\Theta(B))$ であり、$V_2=\pi$ とあわせて
$\sigma_1(B)=\lambda_1(\Theta(B))/2\pi$ である（点 $\theta=0$ は零集合）。

**系 A.21（一様な単位ベクトルのモーメント）.** $u\sim\sigma_{2N-1}$ が $\mathbb C^N$ の単位球面上にあるならば、
$\mathbb E[u_a\bar u_b]=\delta_{ab}/N$、$\mathbb E|u_1|^4=\frac2{N(N+1)}$、$\mathbb E|u_1|^2|u_2|^2=\frac1{N(N+1)}$（$N\ge2$）である。
$u\sim\sigma_{n-1}$ が $S^{n-1}\subset\mathbb R^n$ 上にあるならば、$\mathbb Eu_1^2=\frac1n$、$\mathbb Eu_1^4=\frac3{n(n+2)}$、
$\mathbb Eu_1^2u_2^2=\frac1{n(n+2)}$ である。

*証明.* 定理 A.20(b) のように $u=z/|z|$ とする。$z_a\bar z_b=|z|^2u_a\bar u_b$ であり $|z|$ と $u$ は独立であるから、
$\delta_{ab}=\mathbb E[z_a\bar z_b]=\mathbb E|z|^2\,\mathbb E[u_a\bar u_b]=N\,\mathbb E[u_a\bar u_b]$ である（命題 A.16(a) により
$\mathbb E|z|^2=N$）。同様に $\mathbb E|z_1|^4=N(N+1)\mathbb E|u_1|^4$ かつ
$\mathbb E|z_1|^2|z_2|^2=N(N+1)\mathbb E|u_1|^2|u_2|^2$ であり、ここで独立性により $\mathbb E|z_1|^2|z_2|^2=1$、また独立な
$x,y\sim N(0,\frac12)$ に対し $\mathbb E|z_1|^4=\mathbb E(x^2+y^2)^2=\frac34+2\cdot\frac14+\frac34=2$ である。（
$Y\sim N(0,s^2)$ に対しては、$y^3\cdot ye^{-y^2/2s^2}=-s^2y^3\frac d{dy}e^{-y^2/2s^2}$ であるから、部分積分により
$\mathbb EY^4=3s^2\,\mathbb EY^2=3s^4$ を得る。）実の場合は $\mathbb E|g|^2=n$、
$\mathbb E|g|^4=s_n(2\pi)^{-n/2}I_{n+4}=I_{n+4}/I_n=n(n+2)$（系 A.19）、$\mathbb Eg_1^4=3$、$\mathbb Eg_1^2g_2^2=1$ であり、
$|g|^2$、$|g|^4$ は $u=g/|g|$ と独立である（定理 A.20(a)）。$\square$

これらは補題 6.14 と一致する。補題 6.14 はこれらをユニタリ不変性だけから導いている。逆に、不変性は分布を決定する：

**命題 A.22（不変な分布の一意性）.** $u$ を $\mathbb C^N$ 内のランダムな単位ベクトルで、その分布がすべてのユニタリのもとで
不変であるものとする。このとき $u\sim\sigma_{2N-1}$ である。$N=1$ の場合：円周上の回転不変なボレル確率測度はすべて
正規化された弧長である。

*証明.* $Q$ を $U(N)$ 上でハール分布に従い（定理 5.3）、$u$ と独立なものとする。固定した単位ベクトル $v$ に対し、
$We_1=v$ となる $W\in U(N)$ を選ぶ。補題 5.1 により $QW$ はハール分布に従うので、$Qv=(QW)e_1$ は $Qe_1$ と同じ分布をもつ。
定理 5.3 では $Q$ はジニブル行列のグラム・シュミット因子であるから（補題 5.2）、その第1列 $z$ に対し $Qe_1=z/|z|$ であり、
これは定理 A.20(b) により分布 $\sigma_{2N-1}$ をもつ。$u$ で条件づけると（フビニ、独立性）、すべてのボレル集合 $B$ に対し
$\mathbb P(Qu\in B)=\sigma_{2N-1}(B)$ である。代わりに $Q$ で条件づけ、固定した各ユニタリ $U$ に対し $Uu$ が $u$ と同じ分布をもつことを用いると、
$\mathbb P(Qu\in B)=\mathbb P(u\in B)$ である。円周上の確率測度は恒等写像の分布であり、
$U(1)$ は回転として作用する。$\square$

```python
N, B = 3, 400_000                                   # Theorem A.20 (b), Corollary A.21, in C^3
z = (rng.standard_normal((B, N)) + 1j * rng.standard_normal((B, N))) / math.sqrt(2)
R2 = (np.abs(z) ** 2).sum(axis=1)
u = z / np.sqrt(R2)[:, None]
t = 2.0
p_gamma = 1 - math.exp(-t) * sum(t ** k / math.factorial(k) for k in range(N))   # P(|z|^2 <= t) for Gamma(N)
mc(R2 <= t, p_gamma, "P(|z|^2 <= 2)")
mc(R2 ** 2, N * (N + 1), "E|z|^4")
mc(np.abs(u[:, 0]) ** 4, 2 / (N * (N + 1)), "E|u_1|^4")
mc(np.abs(u[:, 0] * u[:, 1]) ** 2, 1 / (N * (N + 1)), "E|u_1|^2 |u_2|^2")
mc((R2 <= t) * np.abs(u[:, 0]) ** 4, p_gamma * 2 / (N * (N + 1)), "E[1{|z|^2 <= 2} |u_1|^4], independence")
mc((u[:, 0] * u[:, 1].conj()).real, 0.0, "E Re u_1 conj(u_2)")

x3 = rng.uniform(-1, 1, (B, 3))                     # Corollary A.19: V_3 = s_3 / 3 = 4 pi / 3
mc((x3 ** 2).sum(axis=1) <= 1, math.pi / 6, "P(uniform point of [-1,1]^3 in the unit ball)")
g3 = rng.standard_normal((B, 3))
u3 = g3 / np.linalg.norm(g3, axis=1, keepdims=True)
mc(u3[:, 0] ** 4, 3 / 15, "E u_1^4 on S^2")
phase = np.angle(z[:, 0])                           # sigma_1 is normalised arc length
mc((phase >= 0) & (phase < 1), 1 / (2 * math.pi), "P(arg z_1 in [0, 1))")
```

## A.6 ガウス分布のモーメント：イッサーリス・ウィックの公式

**補題 A.23.** $Y$ を共分散 $C=(c_{ij})$ をもつ $\mathbb R^n$ 内の中心化ガウスベクトルとする。このとき
$c_{ij}=\mathbb E[Y_iY_j]$ であり、すべての $t\in\mathbb R^n$ に対し $\mathbb Ee^{\langle t,Y\rangle}=e^{\langle t,Ct\rangle/2}$ かつ
$\mathbb Ee^{|\langle t,Y\rangle|}<\infty$ であり、すべてのモーメント $\mathbb E|Y_i|^k$ は有限である。

*証明.* $Y=Ag$ と書く。すると $\mathbb E[Y_iY_j]=\sum_{k,l}A_{ik}A_{jl}\mathbb E[g_kg_l]=(AA^T)_{ij}$ である。$a=A^Tt$ とすると
$\langle t,Y\rangle=\sum_ka_kg_k$ であり、独立性により
$\mathbb Ee^{\langle t,Y\rangle}=\prod_k\mathbb Ee^{a_kg_k}=\prod_ke^{a_k^2/2}=e^{|A^Tt|^2/2}=e^{\langle t,Ct\rangle/2}$ である。最後に
$e^{|s|}\le e^s+e^{-s}$ かつ $|Y_i|^k\le k!\,e^{|Y_i|}$ である。$\square$

$\{1,\dots,2p\}$ の **ペアリング** とは、$p$ 個の2元ブロックへの分割のことである。その個数は $(2p-1)!!=(2p)!/(2^pp!)$ である
（$1$ の相手は $2p-1$ 通りに選べ、そのあと残りの $2p-2$ 個の元をペアにする；$p$ に関する帰納法）。

**定理 A.24（イッサーリス・ウィック）.** $Y_1,\dots,Y_m$ を中心化ガウスベクトルの座標とする（重複を
許す）。$m$ が奇数ならば $\mathbb E[Y_1\cdots Y_m]=0$ であり、$m=2p$ ならば
$$
\mathbb E[Y_1\cdots Y_{2p}]=\sum_\pi\prod_{\{a,b\}\in\pi}\mathbb E[Y_aY_b],
$$
ここで和は $\{1,\dots,2p\}$ のすべてのペアリング $\pi$ にわたる。同じ公式は、一つの中心化ガウスベクトルの座標の複素線形結合である
複素確率変数 $Z_1,\dots,Z_m$（例えばジニブル行列や GUE 行列の成分とその複素共役）に対しても、
対の因子を $\mathbb E[Z_aZ_b]$（共役をとらない）として成り立つ。

*証明.* ベクトル $(Y_1,\dots,Y_m)$ は与えられたガウスベクトルの線形写像による像であるから、中心化
ガウスベクトルである（定理 A.15(a)）。その共分散を $C$ とする。$t\in\mathbb R^m$ を固定し $S=\langle t,Y\rangle$ とおく。$s\in\mathbb R$ に対し
$\sum_k(sS)^k/k!$ の部分和は $e^{|sS|}$ で抑えられ、これは補題 A.23 により可積分であるから、優
収束定理と補題 A.23 により
$$
\sum_{k\ge0}\frac{s^k}{k!}\mathbb E[S^k]=\mathbb Ee^{sS}=e^{s^2q/2}=\sum_{p\ge0}\frac{s^{2p}q^p}{2^pp!},\qquad q=\langle t,Ct\rangle .
$$
両辺は A.2節の意味で $\mathbb R$ 上のべき級数である：右辺についてはこれは明らかであり、左辺については、すべての $\rho>0$ に対し
$\sum_k\frac{\rho^k}{k!}|\mathbb E[S^k]|\le\mathbb E\sum_k\frac{\rho^k|S|^k}{k!}=\mathbb Ee^{\rho|S|}<\infty$ である
（単調収束定理と補題 A.23）。よってそれらの係数は一致する（差に補題 A.5(c) を適用する）。ゆえにすべての
$t\in\mathbb R^m$ に対し、$\mathbb E\langle t,Y\rangle^m/m!$ は $m$ が奇数ならば $0$ に、$m=2p$ ならば $\langle t,Ct\rangle^p/(2^pp!)$ に等しい。
両辺は $t$ の多項式であるから、それらの係数は一致する（差に系 A.7(a) を適用する）。
$\mathbb E\langle t,Y\rangle^m=\sum_{i_1,\dots,i_m}t_{i_1}\cdots t_{i_m}\mathbb E[Y_{i_1}\cdots Y_{i_m}]$ において、単項式 $t_1t_2\cdots t_m$ は
$(1,\dots,m)$ の置換であるような $m!$ 個の列 $(i_1,\dots,i_m)$ から生じるので、
$\mathbb E\langle t,Y\rangle^m/m!$ におけるその係数は $\mathbb E[Y_1\cdots Y_m]$ である。したがって $m$ が奇数ならばそれは $0$ である。$m=2p$ ならば
$q^p=\sum c_{i_1j_1}\cdots c_{i_pj_p}\,t_{i_1}t_{j_1}\cdots t_{i_p}t_{j_p}$ であり、和は $p$ 個の順序対からなるすべての列にわたる。単項式
$t_1\cdots t_{2p}$ が現れるのは、順序なし対 $\{i_k,j_k\}$ がペアリング $\pi$ をなすときであり、かつそのときに限る。各ペアリングは
$p!\,2^p$ 個の列（対の順序と、各対の中の順序）から生じ、$C$ は対称であるから、それぞれが $\prod_{\{a,b\}\in\pi}c_{ab}$ を寄与する。
$c_{ab}=\mathbb E[Y_aY_b]$ であるから、$2^pp!$ で割ると公式が得られる。
複素の主張については $Z_a=\sum_i\alpha_{ai}Y_i$ と書く。両辺は $(Z_1,\dots,Z_m)$ について多重線形であり――左辺は
$\mathbb E$ の線形性により、右辺は各 $Z_a$ が各項のちょうど一つの因子に現れるから――、各
$Z_a$ が $Y_i$ のいずれかであるときに一致する。展開すればよい。$\square$

例えば標準複素ガウス変数 $z$ に対し $\mathbb E|z|^4=\mathbb E[z\bar zz\bar z]=1\cdot1+0+1\cdot1=2$ である。命題 1.10
と補題 7.10 は $m=4$ の場合を用いる。定理 A.24 は期待値の中で微分することなくそれを証明する。

```python
def pairings(idx):
    """All pairings of the list idx, as lists of pairs."""
    if not idx:
        yield []
        return
    for i in range(1, len(idx)):
        for rest in pairings(idx[1:i] + idx[i + 1:]):
            yield [(idx[0], idx[i])] + rest


def wick(cov, idx):                                       # sum over pairings of the positions of idx
    return sum(math.prod(cov[idx[a], idx[b]] for a, b in p) for p in pairings(list(range(len(idx)))))


assert [len(list(pairings(list(range(2 * p))))) for p in range(1, 5)] == [1, 3, 15, 105]    # (2p-1)!!
Amat = rng.standard_normal((3, 3)) / math.sqrt(3)
cov = Amat @ Amat.T
Y = rng.standard_normal((1_000_000, 3)) @ Amat.T
for idx in ((0, 1, 2, 0), (0, 0, 0, 0), (0, 1, 1, 2, 2, 0)):
    mc(np.prod(Y[:, list(idx)], axis=1), wick(cov, idx), f"E Y{idx}")
mc(np.prod(Y[:, [0, 1, 2]], axis=1), 0.0, "E Y(0, 1, 2), odd")
M = np.array([[1, 1j], [1, -1j]])                  # (z, conj z) = M (x, y) with x, y independent N(0, 1/2)
cz = M @ M.T / 2                                     # E[Z_a Z_b] = [[0, 1], [1, 0]]
assert abs(wick(cz, (0, 1, 0, 1)) - 2) < 1e-12                           # E|z|^4 = 2
```

## A.7 指数変数と一様な間隔

確率変数 $e$ が $(0,\infty)$ 上で密度 $e^{-y}$ をもつとき、$e$ は **指数変数** である（$e\sim\operatorname{Exp}(1)$）という。
$e/\beta$ はレート $\beta$ の指数変数である。$\Delta_n=\{x\in\mathbb R^n:0<x_1<\dots<x_n<1\}$ とおく。

**定理 A.25.** $e_0,\dots,e_n$（$n\ge1$）を独立な $\operatorname{Exp}(1)$ 変数とし、$T=e_0+\dots+e_n$、
$k=1,\dots,n$ に対し $X_k=(e_0+\dots+e_{k-1})/T$ とする。このとき $(X_1,\dots,X_n)$ は $\Delta_n$ 上で密度 $n!$ をもち、$T$ は
$(0,\infty)$ 上で密度 $t^ne^{-t}/n!$ をもち、$(X_1,\dots,X_n)$ と $T$ は独立である。特に
$(e_0,\dots,e_n)/T$ は $T$ と独立である。

*証明.* 写像 $\Phi(x,t)=t\,(x_1,\ x_2-x_1,\ \dots,\ x_n-x_{n-1},\ 1-x_n)$ は開集合
$\Delta_n\times(0,\infty)$ から $(0,\infty)^{n+1}$ の上への全単射であり、その逆写像は $y\mapsto\big(((y_0+\dots+y_{k-1})/\sum_iy_i)_{k=1}^n,\ \sum_iy_i\big)$ である。
どちらも $C^1$ である。$\Phi$ のヤコビ行列において、$\partial/\partial t$ の列の成分の和は $1$ であり、各列
$\partial/\partial x_j$（成分は $t$ と $-t$）の成分の和は $0$ である。最後の行をすべての行の和で置き換えると（これは
行列式を変えない）行 $(0,\dots,0,1)$ が得られ、残りの $n\times n$ ブロックは対角成分が $t$ の下三角行列である。
したがって $|\det D\Phi|=t^n$ である。ベクトル $e=(e_0,\dots,e_n)$ は密度 $e^{-\sum_iy_i}$ をもち、
$\sum_i\Phi(x,t)_i=t$ である。変数変換公式により、ボレル集合 $E\subseteq\Delta_n\times(0,\infty)$ に対し
$$
\mathbb P\big((X,T)\in E\big)=\mathbb P\big(e\in\Phi(E)\big)=\int_{\Phi(E)}e^{-\sum_iy_i}dy=\int_Ee^{-t}t^n\,dx\,dt .
$$
密度 $e^{-t}t^n=n!\cdot\frac{t^ne^{-t}}{n!}$ は積の形であり、$\int_0^\infty t^ne^{-t}dt=n!$ であるから
$\lambda_n(\Delta_n)=1/n!$ である。フビニの定理により $\mathbb P(X\in A,T\in A')=\mathbb P(X\in A)\mathbb P(T\in A')$ である。最後に $(e_0,\dots,e_n)/T$ は
$(X_1,X_2-X_1,\dots,1-X_n)$ であり、$X$ の関数である。$\square$

**系 A.26（一様な間隔）.** $U_1,\dots,U_n$ を独立で $(0,1)$ 上一様分布に従うものとし、その順序統計量を
$U_{(1)}<\dots<U_{(n)}$、間隔を $s_0=U_{(1)}$、$s_k=U_{(k+1)}-U_{(k)}$（$1\le k\le n-1$）、$s_n=1-U_{(n)}$ とする。このとき：

(a) 値が一致する確率は $0$ であり、$(U_{(1)},\dots,U_{(n)})$ は $\Delta_n$ 上で密度 $n!$ をもつ；

(b) $(s_0,\dots,s_n)$ は、$e_i$、$T$ を定理 A.25 のとおりとして $(e_0,\dots,e_n)/T$ と同じ分布をもつ；

(c) $\mathbb Es_i=\frac1{n+1}$、$\mathbb Es_i^2=\frac2{(n+1)(n+2)}$ であり、$i\neq j$ に対し $\mathbb E[s_is_j]=\frac1{(n+1)(n+2)}$ である；

(d) $i\neq j$ に対し $s_i/s_j$ は $e_i/e_j$ と同じ分布をもち、$\mathbb P(e_i/e_j\le t)=\frac t{1+t}$（$t>0$）である。したがって
$0\le i\le n-1$ に対し、間隔比 $\min(s_i,s_{i+1})/\max(s_i,s_{i+1})$ は命題 5.14(a) で求めた
$\min(e_i,e_{i+1})/\max(e_i,e_{i+1})$ と同じ分布、すなわち $[0,1]$ 上で $\mathbb P(r\le t)=\frac{2t}{1+t}$、平均 $2\ln2-1$ の分布をもつ。

区間 $[a,b]$ 上の $n$ 個の独立な一様な点に対しては、アフィン写像 $u\mapsto a+(b-a)u$ はすべての間隔を
$b-a$ 倍する：間隔を $b-a$ で割ったものは (b) と (c) をみたし、(d) の比は変わらない。

*証明.* (a) 二つの座標が一致する集合は有限個の超平面の和集合であり、零集合である（系 A.7(c)）。その外では
$U_{\pi(1)}<\dots<U_{\pi(n)}$ となる置換 $\pi$ がちょうど一つ存在するので、ボレル集合 $E\subseteq\Delta_n$ に対し事象
$\{(U_{(1)},\dots,U_{(n)})\in E\}$ は、零集合を除いて、$\pi$ にわたる
$\{(U_{\pi(1)},\dots,U_{\pi(n)})\in E\}$ の非交和である。それぞれの確率は $\lambda_n(E)$ である。$(U_1,\dots,U_n)$ は密度
$1_{(0,1)^n}$ をもち、座標の置換は $\lambda_n$ と $(0,1)^n$ をともに保つからである。よって確率は $n!\,\lambda_n(E)$ である。
(b) (a) と定理 A.25 により、ベクトル $(U_{(k)})_k$ と $(X_k)_k$ は同じ分布をもち、間隔はどちらについても同じ線形
関数である。(c) $w=(e_0,\dots,e_n)/T$ とする。これは $T$ と独立であり（定理 A.25）、$\mathbb ET=n+1$、
$\mathbb ET^2=\int_0^\infty t^{n+2}e^{-t}dt/n!=(n+1)(n+2)$ である。$e_i=Tw_i$ から $1=\mathbb Ee_i=\mathbb ET\,\mathbb Ew_i$、
$2=\mathbb Ee_i^2=\mathbb ET^2\,\mathbb Ew_i^2$、$1=\mathbb E[e_ie_j]=\mathbb ET^2\,\mathbb E[w_iw_j]$ である。(d) $s_i/s_j=w_i/w_j=e_i/e_j$ であり、
トネリの定理により $\mathbb P(e_i\le te_j)=\int_0^\infty e^{-y}(1-e^{-ty})\,dy=1-\frac1{1+t}$ である。間隔比は
$s_i/s_{i+1}$ の関数（すなわち $x=s_i/s_{i+1}$ における $\min(x,1/x)$ であり、すべての $e_k>0$ であるからほとんど確実に定義される）であるから、
$e_i/e_{i+1}$ の同じ関数と同じ分布をもち、それは独立な $\operatorname{Exp}(1)$ 変数 $e_i,e_{i+1}$ に対する
命題 5.14(a) の比 $r$ である。$\square$

特に、「ポアソン準位」$\lambda_k=e_0+\dots+e_{k-1}$、$k=1,\dots,n+1$（独立な指数
間隔 $e_0,\dots,e_n$）と、整列した $n$ 個の一様な点およびその $n+1$ 個の間隔 $s_0,\dots,s_n$ を考える。
連続する間隔比のベクトル $r_i=\min(d_i,d_{i+1})/\max(d_i,d_{i+1})$、$0\le i<n$ は、
$d=(e_0,\dots,e_n)$ に対しても $d=(s_0,\dots,s_n)$ に対しても同じ同時分布をもつ。すべての $d_i$ を $T$ で割っても比は変わらず、
$(e_0,\dots,e_n)/T$ は $(s_0,\dots,s_n)$ と同じ分布をもつからである（(b)）。そして定理 A.25 により、最初の $n$ 個のポアソン準位を
第 $n+1$ の準位で割ったものは、整列した $n$ 個の一様な点とちょうど同じ分布をもつ。

```python
n, B = 5, 200_000                                       # Corollary A.26
levels = np.sort(rng.uniform(0, 1, (B, n)), axis=1)
s = np.diff(np.concatenate([np.zeros((B, 1)), levels, np.ones((B, 1))], axis=1), axis=1)    # s_0, ..., s_n
e = rng.exponential(1.0, (B, n + 1))
w = e / e.sum(axis=1, keepdims=True)                    # normalised exponentials
for name, sp in (("uniform spacings", s), ("normalised exponentials", w)):
    mc(sp[:, 2], 1 / (n + 1), f"{name}: E s_2")
    mc(sp[:, 2] ** 2, 2 / ((n + 1) * (n + 2)), f"{name}: E s_2^2")
    mc(sp[:, 0] * sp[:, n], 1 / ((n + 1) * (n + 2)), f"{name}: E s_0 s_n")
    for tt in (0.5, 2.0):
        mc(sp[:, 1] / sp[:, 2] <= tt, tt / (1 + tt), f"{name}: P(s_1 / s_2 <= {tt})")
r = np.minimum(s[:, 1:-2], s[:, 2:-1]) / np.maximum(s[:, 1:-2], s[:, 2:-1])    # ratios of interior spacings
mc(r.mean(axis=1), 2 * math.log(2) - 1, "<r>, uniform levels")
T = e.sum(axis=1)
mc(w[:, 0] * T ** 2, n + 2, "E[w_0 T^2] = E w_0 E T^2 (independence)")
```

[目次](README.md)**に戻る**。
