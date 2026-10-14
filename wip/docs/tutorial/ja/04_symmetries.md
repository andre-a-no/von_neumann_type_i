[English](../04_symmetries.md) · [Español](../es/04_symmetries.md) · [Français](../fr/04_symmetries.md) · [Deutsch](../de/04_symmetries.md) · [中文](../zh/04_symmetries.md) · **日本語** · [Русский](../ru/04_symmetries.md)

# 4. 対称性と表現

第2章では、対称性の群がフォン・ノイマン環（対称性作用素の張る空間）とその可換子環（対称なオブザーバブル）を生み出すことを示し（命題 2.7）、一つの例を手で分解した。本章ではその背後にある表現論を与え、それをライブラリが任意の行列のリストに対して実行するアルゴリズムに変え、スピン系の二つの対称性に応用する：$U(1)$（粒子数あるいは磁化）と $SU(2)$（スピンの回転）である。

**本章で学ぶこと**：完全可約性とシューアの補題、既約表現への分解が可換子環と生成される環をどのように決定するか、指標がどのように重複度を数えるか、可換子環のランダムな元がなぜ分解を数値的に見つけるか、電荷セクターとは何か、既約なスピン表現と $L$ 個のスピン $\frac12$ の分解、そして $SU(2)$ 不変なハミルトニアンのブロックが $U(1)$ 不変なものよりはるかに小さい理由。

本章を通じて、$S'$ は $S\subseteq M_n(\mathbb C)$ の可換子環である（第2章）。$S^*=S$ ならば、$S''$ は $S$ **が生成するフォン・ノイマン環** である：これは $S$ を含む $*$-環であり、$S$ を含む任意の $*$-環 $\mathcal B\supseteq S$ は
$\mathcal B=\mathcal B''\supseteq S''$ を満たす（定理 2.2）。$\mathcal A'_h=\{K\in\mathcal A':K=K^*\}$ は実ベクトル空間である。$\mathbb C^k\otimes\mathbb C^m$ には第7章のクロネッカー順序を用いる（$e_a\otimes f_j$ は $0$ から数えて $am+j$ 番目の基底ベクトル）。したがって $1_k\otimes x$ は $x$ の $k$ 個のコピーからなるブロック対角行列である。$|G|$ は有限群の位数であり、$\hbar=1$ とする。

## 4.1 表現とシューアの補題

**定義.** $X$ を対合 $x\mapsto x^*$（$x^{**}=x$）をもつ集合とする。$\mathbb C^k$ 上の $X$ の
**表現** とは、$\rho(x^*)=\rho(x)^*$ を満たす写像 $\rho:X\to M_k(\mathbb C)$ のことである。部分空間 $K$ は、すべての $x$ に対し $\rho(x)K\subseteq K$ のとき **不変** であるという。$k\ge1$ かつ $0,\mathbb C^k$ のみが不変部分空間であるとき、$\rho$ は **既約** であるという。$\rho$ から $\mathbb C^{k'}$ 上の表現 $\rho'$ への **絡作用素** とは、すべての $x$ に対し $T\rho(x)=\rho'(x)T$ を満たす $k'\times k$ 行列 $T$ のことである。それらはベクトル空間 $\operatorname{Hom}(\rho,\rho')$ をなし、$\operatorname{Hom}(\rho,\rho)=\rho(X)'$ である。可逆な絡作用素が存在するとき、$\rho$ と $\rho'$ は **同値** であるといい、$\rho\cong\rho'$ と書く。

三つの場合が現れる：$g^*=g^{-1}$ とした群 $G$ で、$\rho$ が準同型 $G\to U(k)$ であることを要求する場合（**ユニタリ表現**；このとき $\rho(g^{-1})=\rho(g)^{-1}=\rho(g)^*$）；随伴をもつ $*$-環で、$\rho$ が単位的 $*$-準同型である場合；そして $x^*=x$ を満たす生成元のラベル、たとえば
$\{x,y,z\}$ で、スピン行列のようなエルミート行列の組で表現される場合（4.5節）である。$K$ が不変で、$Q\in M_{k\times d}$ の列が $K$ の正規直交基底をなすならば $\rho(x)Q=QQ^*\rho(x)Q$ であり、したがって
$\rho_K(x)=Q^*\rho(x)Q$ はふたたび表現である（$\rho$ が準同型ならば準同型）。これを $\rho$ の $K$ への **制限** という。

**補題 4.1（完全可約性）.** $\rho$ を $\mathbb C^n$ 上の $X$ の表現とする。
(a) $K$ が不変ならば $K^\perp$ も不変である。(b) $\mathbb C^n$ は、各々の上で $\rho$ が既約であるような不変部分空間の直交直和である。

*証明.* (a) $y\in K^\perp$ と $z\in K$ に対し、$\rho(x^*)z\in K$ であるから $\langle z,\rho(x)y\rangle=\langle\rho(x^*)z,y\rangle=0$ である。(b) $n$ に関する帰納法。次元が最小の零でない不変部分空間 $K_1$ は既約である。(a) により $K_1^\perp$ は不変であり、$K_1^\perp$ への制限の不変部分空間は $\rho$ の不変部分空間であるから、帰納法の仮定がそれに適用できる。$\square$

有限群ではユニタリ性は何の代償も要しない：$\rho(g)$ が単に可逆であるだけでも、$P=\sum_g\rho(g)^*\rho(g)$ は
$P\ge\rho(e)^*\rho(e)=1$ を満たすので $P>0$ であり、$\rho(h)^*P\rho(h)=\sum_g\rho(gh)^*\rho(gh)=P$ である（$g\to gh^{-1}$ と置換する）。
したがって $(P^{1/2}\rho(h)P^{-1/2})^*(P^{1/2}\rho(h)P^{-1/2})=P^{-1/2}\rho(h)^*P\rho(h)P^{-1/2}=1$ であり、
$h\mapsto P^{1/2}\rho(h)P^{-1/2}$ は $\rho$ と同値なユニタリ表現である（ワイルのトリック）。

**補題 4.2（シューアの補題）.** $\rho,\rho'$ を $\mathbb C^k,\mathbb C^{k'}$ 上の $X$ の既約表現とし、$T\in\operatorname{Hom}(\rho,\rho')$ とする。
(a) $T=0$ であるか $T$ は可逆である。(b) $\operatorname{Hom}(\rho,\rho)=\mathbb C1$。
(c) $T^*\in\operatorname{Hom}(\rho',\rho)$ であり、$T\neq0$ ならば $\lambda>0$ とユニタリ $U$ により $T=\lambda U$ と書ける。
したがって $\dim\operatorname{Hom}(\rho,\rho')$ は $\rho\cong\rho'$ ならば $1$、そうでなければ $0$ であり、同値な既約表現はユニタリ同値である。

*証明.* (a) と (b) は命題 2.7 の後で群に対して証明したが、その証明が用いるのは、$\ker T$ が $\rho$ で不変であること（$Tv=0$ ならば $T\rho(x)v=\rho'(x)Tv=0$）と、$T$ の像が $\rho'$ で不変であることだけである：よって $\ker T\in\{0,\mathbb C^k\}$ かつ $\operatorname{ran}T\in\{0,\mathbb C^{k'}\}$ である。(b) については、$\lambda$ を $T$ の固有値とすると $T-\lambda1$ は可逆でない絡作用素であり、したがって $0$ である。(c) $T\rho(x^*)=\rho'(x^*)T$ の随伴は $\rho(x)T^*=T^*\rho'(x)$ である。すると $T^*T\in\operatorname{Hom}(\rho,\rho)$ であるから、$T\neq0$ ならば $\lambda>0$ により $T^*T=\lambda^21$ となり（$T^*T\ge0$ は零でない）、$U=T/\lambda$ は (a) により正方行列で $U^*U=1$ を満たす。最後に、$T_1,T_2\in\operatorname{Hom}(\rho,\rho')$ が零でなければ $T_1^{-1}T_2\in\operatorname{Hom}(\rho,\rho)=\mathbb C1$ である。$\square$

すべての $\rho(x)$ が互いに可換ならば（アーベル群）、既約な $\rho$ は一次元である：各 $\rho(x)$ は
$\operatorname{Hom}(\rho,\rho)=\mathbb C1$ に属するので、すべての部分空間が不変となる。

## 4.2 等型分解と可換子環

**定理 4.3（等型分解）.** $\rho$ を $\mathbb C^n$ 上の $X$ の表現とする。
$\mathbb C^{k_c}$ 上の $X$ の互いに非同値な既約表現 $\rho_1,\dots,\rho_C$、$\sum_ck_cm_c=n$ を満たす整数 $m_c\ge1$、およびユニタリ
$V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n$ が存在して、すべての $x$ に対し
$V^*\rho(x)V=\bigoplus_c\rho_c(x)\otimes1_{m_c}$ となる。そのような任意の分解について
$$
V^*\rho(X)'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),\qquad
V^*\rho(X)''V=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},
$$
であり、$m_c=\dim\operatorname{Hom}(\rho_c,\rho)$ である。$X$ の既約表現 $\sigma$ が
$\operatorname{Hom}(\sigma,\rho)\neq0$ を満たすのは、ちょうど一つの $c$ に対して $\sigma\cong\rho_c$ となるときかつそのときに限る。

*証明.* *存在.* 補題 4.1 により、和が $\mathbb C^n$ で、制限 $\sigma_i(x)=Q_i^*\rho(x)Q_i$ が既約であるような、互いに直交する不変部分空間 $W_1,\dots,W_M$ が得られる。$W_i$ を制限の同値性によって類に分け、各類 $c$ において基準となる $\rho_c=\sigma_{i_0}$ を選ぶ。その類の他の各 $i$ について、補題 4.2(c) により $\sigma_i(x)=U_i\rho_c(x)U_i^*$ を満たすユニタリ $U_i$ が得られる。$Q_i$ を $Q_iU_i$ に置き換えれば $Q_i^*\rho(x)Q_i=\rho_c(x)$ となる。類 $c$ の元に $j=0,\dots,m_c-1$ と番号をつけ、それらの基底を $Q_{c,j}$ と書き、$c$ 番目の直和成分上で $V(e_a\otimes f_j)=Q_{c,j}e_a$ とおく。$W_i$ は不変で互いに直交するから $Q_{c,j}^*\rho(x)Q_{d,l}=\delta_{cd}\delta_{jl}\rho_c(x)$ であり、これが主張された形である。

*可換子環.* ここで $V$ を定理の主張にあるような任意の分解とし、$Q_{c,j}$ を（$c$ 番目の直和成分における）$n\times k_c$ 等長写像
$e_a\mapsto V(e_a\otimes f_j)$ とする。いま構成した $V$ に対しては、これらは上の $Q_{c,j}$ である。恒等式
$V^*\rho(x)V=\bigoplus_c\rho_c(x)\otimes1_{m_c}$ は $\rho(x)Q_{d,l}=Q_{d,l}\rho_d(x)$ を意味し、その随伴を $x^*$ に適用すると
$Q_{c,j}^*\rho(x)=\rho_c(x)Q_{c,j}^*$ を得る。$V$ はユニタリであるから、$Y\in M_n$ がすべての $\rho(x)$ と可換であることは、各ブロック
$Y_{cj,dl}=Q_{c,j}^*YQ_{d,l}$ が $Y_{cj,dl}\rho_d(x)=\rho_c(x)Y_{cj,dl}$ を満たすこと、すなわち
$\operatorname{Hom}(\rho_d,\rho_c)$ に属することと同値である。これは $c\neq d$ ならば $0$、$c=d$ ならば $\mathbb C1$ である（補題 4.2）。よって
$Y\in\rho(X)'$ であることと $Y_{cj,dl}=\delta_{cd}(y_c)_{jl}1_{k_c}$、すなわち $V^*YV=\bigoplus_c1_{k_c}\otimes y_c$ であることは同値である。

*二重可換子環.* $Z$ が $\rho(X)'$ と可換ならば、$Z$ は直和成分への射影と可換であり（$y_c=1$、$d\neq c$ に対し $y_d=0$ ととる）、したがって $V^*ZV=\bigoplus_cZ_c$ で、$Z_c$ は
$1_{k_c}\otimes M_{m_c}$ と可換である。補題 2.5 のステップ5（二つのテンソル因子を入れ替えたもの）により
$Z_c\in M_{k_c}\otimes1_{m_c}$ である。逆にこれらの作用素は $\rho(X)'$ と可換である。

*重複度.* 既約な $\sigma$ と $T\in\operatorname{Hom}(\sigma,\rho)$ に対し、同じ議論により
$Q_{c,j}^*T\in\operatorname{Hom}(\sigma,\rho_c)$ であり、これは $\sigma\cong\rho_c$ ならば可逆な $t_c$ による $\mathbb Ct_{c}$、そうでなければ
$0$ である。また $T=\sum_{c,j}Q_{c,j}Q_{c,j}^*T$ である。$\rho_c$ は互いに非同値であるから、$\sigma$ はそのうち高々一つ、たとえば $\rho_c$ と同値であり（そうでなければ $T=0$）、このとき $T=\sum_ja_jQ_{c,j}t_c$ である。逆に各 $Q_{c,j}t_c$ は $\operatorname{Hom}(\sigma,\rho)$ に属し、これら $m_c$ 個の写像は互いに直交する像をもつので、
$\dim\operatorname{Hom}(\sigma,\rho)=m_c$ である。$\square$

群の場合、$\rho(G)''$ は $\rho(G)$ の張る空間であり（2.4節）、定理 4.3 は命題 2.7 を逆向きに読んだものである：そこでは構造定理が既約表現を生み出し、ここでは既約表現が構造を生み出す。いずれにおいても $\dim\rho(X)'=\sum_cm_c^2$、$\dim\rho(X)''=\sum_ck_c^2$ である。$\rho(X)''$ の極小中心射影 $P_c=\sum_jQ_{c,j}Q_{c,j}^*$ の像が $\rho_c$ の **等型成分** である。
定理 4.3 により、任意のエルミートな $H\in\rho(X)'$ はエルミートな $h_c\in M_{m_c}$ により $V(\bigoplus_c1_{k_c}\otimes h_c)V^*$ と
書け、$1_{k_c}\otimes h_c$ は $h_c$ の $k_c$ 個のコピーからなるブロック対角行列である：$H$ のスペクトルは $m_c\times m_c$
ブロック $h_c$ のスペクトルの和集合であり、$h_c$ の各固有値は $k_c$ 回繰り返される。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. $\rho(X)$ が対称性で $H\in\rho(X)'$ がハミルトニアンならば、
> いま述べた $k_c$ 重に縮退した準位が対称性の **多重項** であり、対角化しなければならないブロックの大きさは、多重項の大きさではなく重複度 $m_c$ である。

有限群では、何も分解することなく重複度を計算できる。

**定理 4.4（シューアの直交関係）.** $\rho,\rho'$ を $\mathbb C^k,\mathbb C^{k'}$ 上の有限群 $G$ の既約ユニタリ表現とする。このとき
$$
\frac1{|G|}\sum_{g\in G}\rho(g)_{ij}\,\overline{\rho'(g)_{i'j'}}=
\begin{cases}\delta_{ii'}\delta_{jj'}/k&\text{if }\rho=\rho',\\0&\text{if }\rho\not\cong\rho'.\end{cases}
$$

*証明.* $k\times k'$ 行列 $Y$ に対し $T=\frac1{|G|}\sum_g\rho(g)Y\rho'(g)^*$ とおく。$h\in G$ について、$g\to hg$ と置換し
$\rho'(hg)^*=\rho'(g)^*\rho'(h)^{-1}$ を用いると $\rho(h)T=T\rho'(h)$ を得るので、
$T\in\operatorname{Hom}(\rho',\rho)$ である。$\rho\not\cong\rho'$ ならば $T=0$ であり、$\rho=\rho'$ ならば
$k\lambda=\operatorname{Tr}T=\operatorname{Tr}Y$ を満たす $\lambda$ により $T=\lambda1$ である（補題 4.2）。$Y=|j\rangle\langle j'|$ に対し、成分 $T_{ii'}$ は左辺であり、$\operatorname{Tr}Y=\delta_{jj'}$ である。$\square$

表現の **指標** とは $\chi_\rho(g)=\operatorname{Tr}\rho(g)$ のことである。同値な表現は等しい指標をもつ。実際
$\operatorname{Tr}(T\rho(g)T^{-1})=\operatorname{Tr}\rho(g)$ である。
$\langle\varphi,\psi\rangle_G=\frac1{|G|}\sum_g\overline{\varphi(g)}\psi(g)$ と書く。

**系 4.5（指標）.** $G$ を有限群とする。
(a) 既約ユニタリな $\rho,\rho'$ に対し、$\langle\chi_\rho,\chi_{\rho'}\rangle_G$ は $\rho\cong\rho'$ ならば $1$、そうでなければ $0$ である。
(b) 定理 4.3 のように分解されたユニタリ表現 $\pi$ と既約な $\rho$ に対し、$\pi$ における $\rho$ の重複度は
$$
m_\rho=\frac1{|G|}\sum_{g\in G}\overline{\chi_\rho(g)}\,\chi_\pi(g)
$$
である（$\rho\cong\rho_c$ ならば $m_\rho=m_c$、$\rho$ がどの $\rho_c$ とも同値でなければ $m_\rho=0$）。
(c) $\dim\pi(G)'=\sum_cm_c^2=\frac1{|G|}\sum_{g\in G}|\chi_\pi(g)|^2$。特に、$\pi$ が既約であるのはこの数が $1$ のときかつそのときに限る。

*証明.* (a) $\rho'$ を同値な表現（同じ指標をもつ）に置き換える。定理 4.4 を
$i=j$、$i'=j'$ について和をとり、複素共役をとる：結果は $\rho=\rho'$ ならば $\sum_{i,i'}\delta_{ii'}/k=1$、
$\rho\not\cong\rho'$ ならば $0$ である。(b) トレースはユニタリ不変であり $\operatorname{Tr}(a\otimes1_m)=m\operatorname{Tr}a$ であるから、
$\chi_\pi=\sum_cm_c\chi_{\rho_c}$ である。(a) を適用すればよい。(c) 同じ展開と (a) により
$\langle\chi_\pi,\chi_\pi\rangle_G=\sum_cm_c^2$ であり、これは定理 4.3 により $\dim\pi(G)'$ である。$\square$

(c) の第二の証明は、ライブラリが可換子環からサンプリングする方法そのものである。**ツワール**
$P(Y)=\frac1{|G|}\sum_g\pi(g)Y\pi(g)^*$ は像が $\pi(G)'$ の冪等写像であり（第2章、演習問題 3）、
$\langle Y,Z\rangle_{HS}=\operatorname{Tr}(Y^*Z)$ に関して自己共役である（$g\to g^{-1}$ と置換する）。したがって $\pi(G)'$ の上への直交射影であり、$\dim\pi(G)'$ はそのトレースである。写像 $Y\mapsto AYB$ のトレースは
$\operatorname{Tr}A\operatorname{Tr}B$ であり（その対角成分は $A_{ii}B_{jj}$）、これから $\frac1{|G|}\sum_g|\chi_\pi(g)|^2$ が得られる。

2.4節の例を指標を用いてやり直そう。置換行列のトレースはそれが固定する基底ベクトルの個数であり、$\pi(\sigma)$ が $|b_1b_2b_3\rangle$ を固定するのは $b$ が $\sigma$ の各巡回上で一定であるときかつそのときに限るので（補題 B.1(c)）、$\chi_\pi(\sigma)=2^{\#\text{cycles}(\sigma)}$ である。$S_3$ の三つのユニタリ表現が手元にある。$P_\sigma e_i=e_{\sigma(i)}$ を $\mathbb C^3$ 上の置換表現（置換行列への準同型）とする。**自明** 表現は $\sigma\mapsto1$ である。**符号** は
$\operatorname{sgn}\sigma=\det P_\sigma$ であり、行列式は乗法的で置換行列の行列式は $\pm1$ であるから、これは $\{\pm1\}\subseteq U(1)$ への準同型である。直線 $\mathbb C(1,1,1)$ は不変であり、したがって $E=(1,1,1)^\perp$ も不変である（補題 4.1(a)）。**標準** 表現は $P$ の $E$ への制限である。
$\mathbb C(1,1,1)\oplus E$ に適合した正規直交基底では $P_\sigma=1\oplus\operatorname{std}(\sigma)$ であるから、
$\chi_{\rm std}(\sigma)=\operatorname{Tr}P_\sigma-1$、すなわち固定点の個数から $1$ を引いたものである。$S_3$ の六つの元は、
単位元、三つの互換、二つの $3$-巡回であり、巡回の個数はそれぞれ $3,2,1$、固定点の個数はそれぞれ $3,1,0$ である。互換は
$1_3$ の二つの列を入れ替えるのでその符号は $-1$ であり、$3$-巡回は二つの互換の積（$(1\,2\,3)=(1\,3)(1\,2)$）であるからその符号は $+1$ である。こうして次の表を得る。

| | 単位元（$1$） | 互換（$3$） | $3$-巡回（$2$） |
|---|---|---|---|
| $\chi_{\rm triv}$ | $1$ | $1$ | $1$ |
| $\operatorname{sgn}$ | $1$ | $-1$ | $1$ |
| $\chi_{\rm std}$ | $2$ | $0$ | $-1$ |
| $\chi_\pi$ | $8$ | $4$ | $2$ |

（したがって $\operatorname{sgn}\sigma=(-1)^{3-\#\text{cycles}(\sigma)}$ であり、すべての指標は実数値である。）$\langle\varphi,\psi\rangle_G$ を共役類ごとに和をとって計算すると、
$\langle\chi_{\rm triv},\chi_{\rm triv}\rangle_G=\frac{1+3+2}6=1$、$\langle\operatorname{sgn},\operatorname{sgn}\rangle_G=\frac{1+3+2}6=1$、
$\langle\chi_{\rm std},\chi_{\rm std}\rangle_G=\frac{4+0+2}6=1$ であるから、系 4.5(c) によりこれら三つの表現は既約であり、指標が異なるので互いに非同値である。系 4.5(b) により $\pi$ における重複度は
$$
m_{\rm triv}=\tfrac{8+3\cdot4+2\cdot2}6=4,\qquad m_{\rm sgn}=\tfrac{8-3\cdot4+2\cdot2}6=0,\qquad
m_{\rm std}=\tfrac{2\cdot8+0-2\cdot2}6=2 .
$$
である。したがって自明表現と標準表現は、$m_c=4$ および $m_c=2$ をもつ $\pi$ の構成成分 $\rho_c$（定理 4.3）と同値である。これらで
$4\cdot1+2\cdot2=8=\sum_ck_cm_c$ 次元が尽くされるので、他の構成成分はない：$k=(1,2)$、$m=(4,2)$ であり、
$\dim\pi(S_3)'=4^2+2^2=\frac16(64+3\cdot16+2\cdot4)=20$ である（系 4.5(c)）。コードはこの表を計算し直す：$P_\sigma$ を
$|100\rangle,|010\rangle,|001\rangle$ の張る不変部分空間への $\pi(\sigma)$ の制限として求め（その上で $\pi(\sigma)$ はただ一つの $1$ の位置を置換する。補題 B.1(a)）、
$\#\text{cycles}(\sigma)=\dim\ker(P_\sigma-1)$（各巡回上で一定なベクトル）を用いる。

```python
import math
import torch
import torch_vn_algebra.decompose as dc

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    return (torch.as_tensor(X) - torch.as_tensor(Y)).abs().max().item() < tol

G = dc.permutation_representation(3, 2)                  # (6, 8, 8): pi(sigma) permutes three qubits
chi_pi = G.diagonal(dim1=-2, dim2=-1).sum(-1).real       # characters Tr pi(sigma)
P3 = G[:, [4, 2, 1]][:, :, [4, 2, 1]].real              # on |100>, |010>, |001>: S_3 permuting C^3
cycles = 3 - torch.linalg.matrix_rank(P3 - torch.eye(3, dtype=torch.float64)).double()   # dim ker(P - 1)
assert close(chi_pi, 2.0 ** cycles)
irreps = {'trivial': torch.ones(6, dtype=torch.float64), 'sign': (-1.0) ** (3 - cycles),
          'standard': P3.diagonal(dim1=-2, dim2=-1).sum(-1) - 1}                 # fixed points - 1
chis = torch.stack(list(irreps.values()))
assert close(chis @ chis.T / 6, torch.eye(3, dtype=torch.float64))    # Corollary 4.5(a): all irreducible
mult = chis @ chi_pi / 6                                               # Corollary 4.5(b)
print(dict(zip(irreps, mult.tolist())), (chi_pi ** 2).mean().item())
assert close(mult, torch.tensor([4.0, 0.0, 2.0], dtype=torch.float64))
assert dc.commutant_dimension(G) == 20                                # Corollary 4.5(c): M_4 (+) M_2
```

これは2.4節の分解であり、今度は表現ごとに六つの数からなる表から得られた。

## 4.3 生成元から環へ：ランダム化アルゴリズム

定理 4.3 は分解を記述するが、それを計算はしない：その証明には既約な不変部分空間が必要である。Murota、Kanno、Kojima、Kojima（*Japan J. Indust. Appl. Math.* 27 (2010) 125–160。Maehara and Murota, *ibid.* 263–293 も参照）は、*可換子環のランダムな元* の固有空間がそれらを見つけることに気づいた。$\mathcal A$ を与えられた行列 $S_1,\dots,S_J$ とその随伴が生成するフォン・ノイマン環とし、その分解を $V^*\mathcal AV=\bigoplus_cM_{k_c}\otimes1_{m_c}$、
$V^*\mathcal A'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$ とする（定理 2.6）。$V$ を用いるのは解析だけである。

**補題 4.6.** $K\in\mathcal A'_h$ とし、$V^*KV=\bigoplus_c1_{k_c}\otimes\kappa_c$ とする。$K$ の各固有空間は
$\mathcal A$ 不変である。すべての $\kappa_c$ を合わせた $N=\sum_cm_c$ 個の固有値が互いに相異なるならば、$K$ の固有空間はちょうど $N$ 個の部分空間 $V(\mathbb C^{k_c}\otimes e)$（$c$ 番目の直和成分において、$e$ は $\kappa_c$ の固有ベクトル）であり、$\mathcal A$ はそれぞれの上に $M_{k_c}$ として既約に作用する。

*証明.* $Kv=\lambda v$ かつ $a\in\mathcal A$ ならば $Kav=aKv=\lambda av$ である。正規直交な $e_{cj}$ により
$\kappa_c=\sum_j\mu_{cj}|e_{cj}\rangle\langle e_{cj}|$ と書く（定理 1.3）。すると
$V^*KV=\sum_{c,j}\mu_{cj}R_{cj}$ であり、ここで $R_{cj}$ は $\mathbb C^{k_c}\otimes e_{cj}$ の上への射影である。
$R_{cj}$ は互いに直交して和が $1$ であり、$\mu_{cj}$ は相異なるので、これらがスペクトル射影である。
$\mathbb C^{k_c}\otimes e_{cj}$ 上で元 $V^*aV=\bigoplus_da_d\otimes1$ は $a_c$ として作用し、$a_c$ は
$M_{k_c}$ 全体を動く。$M_{k_c}$ は $0$ と $\mathbb C^{k_c}$ 以外に不変部分空間をもたない。$\square$

**命題 4.7（一般性）.** $K$ と $Y$ を、それぞれ $\mathcal A'_h$ と
$\mathcal A'$ の独立なランダムな元で、その分布がこれらの実ベクトル空間上のルベーグ測度に関して密度をもつもの（補題 A.3）とする。このとき確率 $1$ で：(a) 補題 4.6 の仮定が成り立つ。(b) 像が同じ直和成分に含まれる $K$ の任意の二つの相異なるスペクトル射影 $\Pi_i,\Pi_j$ に対し $\Pi_iY\Pi_j\neq0$ である。

*証明.* (a) エルミート行列 $\kappa=\bigoplus_c\kappa_c\in M_N$ は $K$ に実線形に依存する。$\langle\cdot,\cdot\rangle_{HS}$ に関する $1,\kappa,\dots,\kappa^{N-1}$ のグラム行列式
$\delta(K)=\det[\operatorname{Tr}\kappa^{i+j}]_{i,j=0}^{N-1}$ は $K$ の座標の実多項式であり、これらの冪が線形独立であるときかつそのときに限り零でない。
$\kappa=U\operatorname{diag}(\lambda)U^*$ に対して $p(\kappa)=U\operatorname{diag}(p(\lambda_l))U^*$ であるから、次数 $<N$ の零でない多項式が $\kappa$ を零化するのは、それが $\operatorname{spec}\kappa$ 上で消えるときかつそのときに限り、それが可能なのは $\kappa$ の相異なる固有値が $N$ 個未満のときかつそのときに限る（それらの上で $\prod_\mu(t-\mu)$ をとる）。よって $\delta(K)\neq0$ は補題 4.6 の仮定であり、
$\delta\not\equiv0$ である：実際 $V(\bigoplus_c1\otimes\operatorname{diag}(\text{相異なる実数}))V^*\in\mathcal A'_h$ において零でない。したがって
系 A.7(b) により $\mathbb P(\delta(K)=0)=0$ である。

(b) 集合 $O=\{\delta\neq0\}\subseteq\mathcal A'_h$ は開集合である（$\delta$ は連続）。$K\in O$ に対し、$K$ はちょうど $N$ 個の相異なる固有値
$\mu_1(K)<\dots<\mu_N(K)$ をもつ（$\kappa$ の固有値で、重複度は $k_c$）ので、固有値の昇順に番号をつけたそのスペクトル射影
$\Pi_1(K),\dots,\Pi_N(K)$ は $O$ 上で連続である（命題 A.13）。$l\neq l'$ を固定し、$\mathcal N_{ll'}$ を、ある $c$ に対して
$P_c\Pi_l(K)=\Pi_l(K)$ かつ $P_c\Pi_{l'}(K)=\Pi_{l'}(K)$（$P_c$ は中心射影）であり、かつ $\Pi_l(K)Y\Pi_{l'}(K)=0$ であるような
$(K,Y)\in O\times\mathcal A'$ 全体の集合とする。これはボレル集合 $O\times\mathcal A'$ 上の連続な方程式で定義されるので、ボレル集合である。二つの像が直和成分 $c$ に含まれるような $K\in O$ を固定する。補題 4.6 により、
$V^*\Pi_lV$ と $V^*\Pi_{l'}V$ は、$\kappa_c$ の正規直交な固有ベクトル $e,e'$ により、それぞれ $\mathbb C^{k_c}\otimes e$ と $\mathbb C^{k_c}\otimes e'$ の上への射影であり、
$V^*YV=\bigoplus_d1\otimes y_d$ に対して $\Pi_lY\Pi_{l'}=0$ であることは $\langle e,y_ce'\rangle=0$ と同値である。写像
$Y\mapsto\langle e,y_ce'\rangle$ は実線形写像 $\mathcal A'\to\mathbb C$ であって零でない（
$V(0\oplus\dots\oplus1\otimes|e\rangle\langle e'|\oplus\dots\oplus0)V^*$ における値は $1$）ので、その核は真部分空間であり、零集合である（系 A.7(c)）。
よってすべての切り口 $\{Y:(K,Y)\in\mathcal N_{ll'}\}$ は零集合であり（二つの像が異なる直和成分に含まれるならば空集合）、$Y$ は $K$ と独立であるから、
系 A.7(d) により $\mathbb P((K,Y)\in\mathcal N_{ll'})=0$ である。(b) の事象が成り立たないのは
$\{\delta(K)=0\}\cup\bigcup_{l\neq l'}\{(K,Y)\in\mathcal N_{ll'}\}$ の上だけであり、これは確率 $0$ の事象の有限個の和集合である。$\square$

**補題 4.8（圧縮；コピーの整列）.** (a) $W$ を正規直交基底 $Q\in M_{n\times d}$ をもつ $\mathcal A$ 不変部分空間とし、
$\sigma(a)=Q^*aQ$ をその制限とする。このとき $\sigma(\mathcal A)$ は $*$-環であり、
$\{Q^*YQ:Y\in\mathcal A'\}=\sigma(\mathcal A)'$ かつ $\{Q^*KQ:K\in\mathcal A'_h\}=\sigma(\mathcal A)'_h$ である。$\mathcal A$ が $W$ 上に既約に作用するのは、
すべての $Y\in\mathcal A'$ に対して $Q^*YQ\in\mathbb C1$ となるときかつそのときに限る。そうでないとき、$X$ を $\mathcal A'$ または
$\mathcal A'_h$ の密度をもつランダムな元とすると、確率 $1$ で $Q^*XQ\notin\mathbb C1$ である。

(b) $W_i,W_j$ を、$\mathcal A$ がその上に既約に作用する不変部分空間とし、その正規直交基底を $Q_i,Q_j$、制限を $\sigma_i(a)=Q_i^*aQ_i$ とする。また $Y\in\mathcal A'$ とする。
このとき $T=Q_i^*YQ_j\in\operatorname{Hom}(\sigma_j,\sigma_i)$ である。$T\neq0$ ならば $\lambda>0$ とユニタリ $U$ により $T=\lambda U$ と書け、$Q_j'=Q_jU^*$ はすべての $a\in\mathcal A$ に対し $Q_j'^*aQ_j'=\sigma_i(a)$ を満たす。$W_i,W_j$ が異なる直和成分に含まれるならば（$W_i\subseteq P_c\mathbb C^n$、$W_j\subseteq P_d\mathbb C^n$、$c\neq d$）、$T=0$ である。

*証明.* (a) 定理 4.3 の証明と同様に（$W$ は $a$ と $a^*$ で不変である）、$a\in\mathcal A$ に対し $aQ=Q\sigma(a)$ かつ
$Q^*a=\sigma(a)Q^*$ である。したがって $\sigma(ab)=Q^*aQ\sigma(b)=\sigma(a)\sigma(b)$、$\sigma(a^*)=\sigma(a)^*$、
$\sigma(1)=Q^*Q=1_d$ であるから、$\sigma$ は単位的 $*$-準同型であり、$\sigma(\mathcal A)$ は $*$-環である。また $Y\in\mathcal A'$ に対し $Q^*YQ\sigma(a)=Q^*YaQ=Q^*aYQ=\sigma(a)Q^*YQ$ である。
逆に、$T\in\sigma(\mathcal A)'$ に対し $Y=QTQ^*$ とおくと、$aY=Q\sigma(a)TQ^*=QT\sigma(a)Q^*=Ya$ であるから $Y\in\mathcal A'$ であり、
$Q^*YQ=T$ で、$T$ がエルミートならば $Y$ もエルミートである。$\sigma$ が既約ならば $\sigma(\mathcal A)'=\operatorname{Hom}(\sigma,\sigma)=\mathbb C1$ である
（補題 4.2(b)）。$E\subseteq\mathbb C^d$ が $0,\mathbb C^d$ 以外の $\sigma$ 不変部分空間ならば、$E^\perp$ もそうであり（補題
4.1(a)）、$E$ の上への直交射影は $\sigma(\mathcal A)'_h$ のスカラーでない元である。この場合
$\mathbb C1\cap\sigma(\mathcal A)'$ と $\mathbb R1\cap\sigma(\mathcal A)'_h$ は真の実部分空間であり、全射な実線形写像 $X\mapsto Q^*XQ$ によるそれらの逆像は
$\mathcal A'$ と $\mathcal A'_h$ の真部分空間であって、これらは零集合である（系 A.7(c)）。

(b) (a) と同様に $aQ_j=Q_j\sigma_j(a)$ かつ $Q_i^*a=\sigma_i(a)Q_i^*$ であるから、
$T\sigma_j(a)=Q_i^*YaQ_j=Q_i^*aYQ_j=\sigma_i(a)T$ である。補題 4.2 により $T=\lambda U$ かつ
$U\sigma_j(a)U^*=\sigma_i(a)$ であり、後者は $Q_j'^*aQ_j'$ である。最後に $P_c\in\mathcal A\cap\mathcal A'$、$Q_i=P_cQ_i$、
$P_cQ_j=0$ であるから、$T=Q_i^*P_cYQ_j=Q_i^*YP_cQ_j=0$ である。$\square$

以上により、アルゴリズム（`decompose.wedderburn`、フロントエンドは `from_group` と `from_generators`）は次のとおりである：

1. *可換子環からのサンプリング.* 有限群の全元に対しては、ツワール $P$ と独立なジニブル行列 $Z,Z'$ により $Y=P(Z)$、$K=P(Z'+Z'^*)$ とする。
   $P$ は $M_n$ から $\mathcal A'$ の上への線形射影で $P(Z)^*=P(Z^*)$ を満たすので、$Z\mapsto P(Z)$ と $Z'\mapsto P(Z'+Z'^*)$ は $\mathcal A'$ の上および
   $\mathcal A'_h$ の上への実線形全射である（$H\in\mathcal A'_h$ に対し $P(\frac H2+\frac H2)=H$）。ジニブル行列は $M_n$ の非退化な
   中心化ガウスベクトルである（命題 A.16）から、$Y$ と $K$ は $\mathcal A'$ と $\mathcal A'_h$ の非退化な中心化ガウスベクトルであり、
   密度をもつ（定理 A.15）。生成元に対しては、
   $\mathcal A'=\{X:[S_j,X]=[S_j^*,X]=0\}$ は $M_n$ 上の線形写像の核である（ライブラリはエルミート部分
   $\frac12(S_j+S_j^*)$、$\frac1{2i}(S_j-S_j^*)$ を用いる。これらは $S_j,S_j^*$ と同じ複素空間を張るので、同じ可換子環をもつ）。その基底 $(B_r)_{r\le D}$ をとり、
   $Y=\sum_rc_rB_r$、$K=Y'+Y'^*$ とする。ここで $Y'$ は $Y$ の独立なコピー、$c_r$ は独立な標準複素ガウス変数である：これらは
   $\mathbb C^D$ の非退化なガウスベクトルの、$\mathcal A'$ の上への実線形全単射 $c\mapsto\sum_rc_rB_r$ による像、および $\mathcal A'_h$ の上への実線形全射
   $c\mapsto X+X^*$（$X=\sum_rc_rB_r$）による像であり、したがってやはり密度をもつ。どのランダムな元も新たに抽出されたもので、それ以前のものとは独立である。
   ライブラリはそれらを単位ノルムに規格化するが、これはどの固有空間も、どのブロックが消えるかも変えない。
2. *分割.* $K$ の固有空間は既約である（補題 4.6、命題 4.7(a)）。それらの次元が $k_c$ である。
3. *グループ化と整列.* 等しい次元の固有空間の間のブロック $Q_i^*YQ_j$ は、ちょうど同じ直和成分の内部で零でない（命題 4.7(b)、
   補題 4.8(b)）。「零でないブロックで結ばれる」という関係の類がセクターを与え、その大きさが $m_c$ を与え、極分解の因子が各コピーを類ごとに一つの基準コピーに整列させる
   （補題 4.8(b)。整列されたコピーに整列されたコピーは基準に整列されている）。
4. クロネッカー順序で $V$ を *組み立て*、$V^*aV=\bigoplus_ca_c\otimes1_{m_c}$ となるようにする。

したがって、厳密な算術で $K$ の固有空間が厳密に分離されるならば、アルゴリズムは確率 $1$ で定理 4.3 のような分解を返す：補題 4.6 により直和成分 $c$ は次元 $k_c$ の固有空間を $m_c$ 個含み、
整列の後の構成は定理 4.3 の証明の存在部分そのものである。浮動小数点演算では、ライブラリは許容誤差より近い固有値を等しいものとして扱うので、計算された「固有空間」$W$ は $K$ のいくつかの
固有空間の和であるかもしれない。それでも $W$ は $\mathcal A$ 不変であり（補題 4.6）、ライブラリは新しいランダムな
$X\in\mathcal A'_h$ をそれに圧縮する。補題 4.8(a) により、$\mathcal A$ が $W$ 上に既約に作用するならば $Q^*XQ$ は $1$ の倍数であり、
そうでなければ確率 $1$ で $1$ の倍数ではない。その場合、その固有空間は $\sigma(\mathcal A)$ で不変であり（$\sigma(\mathcal A)$ に適用した補題 4.6）、
$W$ をより小さな不変部分空間に分割する。ライブラリはすべての圧縮がスカラーになるまでこれを繰り返す。浮動小数点演算で何が起ころうと、結果は最後に証明書によって確認される：

**命題 4.9（証明書）.** $\mathcal A$ を $S_1,\dots,S_J$ とその随伴が生成する環とし、
ユニタリ $V$ により $\mathcal B=V(\bigoplus_cM_{k_c}\otimes1_{m_c})V^*$ とする。すべての $S_j\in\mathcal B$ かつ
$\dim\mathcal A'=\sum_cm_c^2$ ならば、$\mathcal A=\mathcal B$ である。

*証明.* $\mathcal B$ はすべての $S_j,S_j^*$ を含む $*$-環であるから $\mathcal A\subseteq\mathcal B$ かつ
$\mathcal B'\subseteq\mathcal A'$ である。$\mathbb C^n$ 上の $\mathcal B$ の表現に定理 4.3 を適用すると（ブロック写像 $b\mapsto b_c$ は既約であり、かつ非同値である。$c$ 番目の中心射影はそのうちの一つによって $1$ に、他によって $0$ に写されるからである）、$\dim\mathcal B'=\sum_cm_c^2=\dim\mathcal A'$ を得る。よって $\mathcal B'=\mathcal A'$ であり、
$\mathcal B=\mathcal B''=\mathcal A''=\mathcal A$ である。$\square$

`Decomposition.verify()` は、生成元（および可換子環の新しい元）の $\mathcal B$（あるいは $\mathcal B'$）からの相対距離を報告し、$\sum_cm_c^2$ を $\dim\mathcal A'$（群に対しては指標から、系 4.5(c)。生成元に対しては核から）と比較する。$S_3$ と三つの量子ビットに対してこれを実行し、それから中身を見てみる：一つのランダムな $K$ の固有空間と、一つのランダムな $Y$ の絡作用素である。

```python
D = dc.from_group(G)                                     # the twirl samples pi(S_3)'
print(D)
assert D.k_factors == (1, 2) and D.multiplicities == (4, 2)            # Section 2.4: k = (1, 2), m = (4, 2)
assert D.verify()['ok']
rho = D.to_algebra(G).matrix                             # (6, 2, 2, 2): the blocks rho_c(sigma), padded
assert close(rho[:, 0, 0, 0], torch.ones(6))                           # sector 1: trivial representation
assert close(rho[:, 1].diagonal(dim1=-2, dim2=-1).sum(-1), chis[2])    # sector 2: the standard one
schur = torch.einsum('gij,gkl->ijkl', rho[:, 1], rho[:, 1].conj()) / 6
assert close(schur.reshape(4, 4), torch.eye(4) / 2)                    # Theorem 4.4 with k = 2
assert close(D.from_algebra(D.to_algebra(G)), G)

Z = torch.randn(8, 8, dtype=dt)
K = dc.twirl(G, Z + Z.mH)                                # a random element of pi(S_3)'_h
w, Q = torch.linalg.eigh(K)
cuts = [0] + [i + 1 for i in range(7) if w[i + 1] - w[i] > 1e-8] + [8]
spaces = [Q[:, a:b] for a, b in zip(cuts[:-1], cuts[1:])]
assert sorted(W.shape[1] for W in spaces) == [1, 1, 1, 1, 2, 2]      # Lemma 4.6: m_c pieces of size k_c
assert all(close(G @ W, W @ (W.mH @ G @ W)) for W in spaces)           # invariant subspaces
Y = dc.twirl(G, torch.randn(8, 8, dtype=dt))
Wa, Wb = [W for W in spaces if W.shape[1] == 2]
s = torch.linalg.svdvals(Wa.mH @ Y @ Wb)                               # Lemma 4.8: a multiple of a unitary
assert s[1] > 1e-3 and abs(s[0] - s[1]) < 1e-10
assert all(close(W.mH @ Y @ Wa, 0) for W in spaces if W.shape[1] == 1)   # different summands: zero
```

計算量を支配するのは可換子環である：ツワールのサンプル一つに $n\times n$ 行列の積が $|G|$ 回かかる。生成元に対しては、
ライブラリはまず $\mathcal A$ のエルミートな元 $R=U\operatorname{diag}(\lambda)U^*$（生成元のエルミート部分のランダムな実線形結合に、
そのような二つの結合の対称化された積の定数倍を加えたもの）を対角化する。任意の $X\in\mathcal A'$ は $R$ と可換であり、したがってそのスペクトル射影とも可換である（定理 1.5(d)）から、
$U^*XU$ は $R$ の固有空間（次元 $d_\lambda$）に関してブロック対角である。これをあらかじめ課すことにより、核の問題の未知数は $n^2$ 個ではなく
$N_R=\sum_\lambda d_\lambda^2$ 個になる。いくつかの固有値をまとめてもブロックが大きくなるだけなので、縮小された問題は同じ解集合 $\mathcal A'$ をもつ。
分解において $V^*RV=\bigoplus_cr_c\otimes1_{m_c}$（エルミートな $r_c\in M_{k_c}$）である。$\lambda$ の $r_c$ における重複度を
$\mu_c(\lambda)$ とすると $d_\lambda=\sum_c\mu_c(\lambda)m_c$ であり、整数 $\mu\ge0$ に対して $\mu^2\ge\mu$ であるから
$$
N_R=\sum_\lambda\Big(\sum_c\mu_c(\lambda)m_c\Big)^2\ge\sum_{\lambda,c}\mu_c(\lambda)^2m_c^2\ge\sum_{\lambda,c}\mu_c(\lambda)m_c^2=\sum_ck_cm_c^2,
$$
であり、等号はすべての $r_c$ の固有値が単純かつ互いに相異なるときかつそのときに限り成り立つ。$L=3,4,5$ 個のスピンの全スピンについて、下のコードは等号を見いだす（*数値的観察*）。
最小値でさえ大きい：八つのスピン $\frac12$（$n=256$）では $\sum_ck_cm_c^2=4900$ であり、縮小された連立方程式は $(n^2+N_R)\times N_R$ 行列で、
$3.5\cdot10^8$ 個の複素成分をもつ。そのため4.6節ではスピン鎖に対して、答えをあらかじめ知っているモジュール `su2` を用いる。

```python
PLACEHOLDER
```

## 4.4 $U(1)$ と電荷セクター

**命題 4.10.** $Q$ を整数の固有値 $q$（**電荷**）をもつ $\mathbb C^n$ 上のエルミート行列とし、そのスペクトル射影を $P_q$、$d_q=\operatorname{rank}P_q$、$U_\theta=e^{i\theta Q}$ とする。
$U(1)'=\{U_\theta:\theta\in\mathbb R\}'$ と書き、$U(1)''$ でその可換子環を表す。

1. $\theta\mapsto U_\theta$ は $U(1)=\mathbb R/2\pi\mathbb Z$ のユニタリ表現であり、$\mathbb C^n$ 上の $U(1)$ の任意の連続な
   ユニタリ表現はこの形をもつ。
2. $X\in M_n$ に対し次は同値である：(i) すべての $\theta$ に対し $XU_\theta=U_\theta X$。(ii) $[X,Q]=0$。
   (iii) $X=\sum_qP_qXP_q$。したがって $U(1)'\cong\bigoplus_qM_{d_q}(\mathbb C)$（ブロック対角作用素）であり、
   $U(1)''=\operatorname{span}\{P_q\}$ である：すべての $k_c=1$ で、重複度は $d_q$ である。
3. （*保存電荷.*）エルミートな $H$ に対し、すべての $t$ について $e^{itH}Qe^{-itH}=Q$ であることは $[H,Q]=0$、すなわち
   $H\in U(1)'$ と同値である。このとき $\rho(t)=e^{-itH}\rho e^{itH}$ のすべてのセクターの重み $\operatorname{Tr}(P_q\rho(t))$ は一定である。

*証明.* (1) $U_\theta=\sum_qe^{iq\theta}P_q$ はユニタリであり、$U_\theta U_{\theta'}=U_{\theta+\theta'}$（定理 1.5(c)）、
$U_{2\pi}=1$ である。逆に、連続な表現に対して $(U_\theta)_{\theta\ge0}$ は行列の連続半群であるから、補題 8.5（第8章で証明する。その証明は有限次元であることしか用いない）により $\theta\ge0$ に対し $U_\theta=e^{\theta B}$ であり、$\theta<0$ に対しては $U_\theta=U_{-\theta}^{-1}=e^{\theta B}$ である。
$U_\theta^*U_\theta=1$ を $0$ で微分すると（補題 8.1）$B^*+B=0$ となるので、エルミートな $Q$ により $B=iQ$ と書ける。指数級数
$e^{i\theta Q}$ は関数計算に等しい（$H=-Q$ とした命題 8.2(1)）から、$U_{2\pi}=\sum_qe^{2\pi iq}P_q=1$ より $Q$ のすべての固有値 $q$ について
$e^{2\pi iq}=1$、$q\in\mathbb Z$ が従う。(2) (i)$\Rightarrow$(ii)：$\theta=0$ で微分する。
(ii)$\Rightarrow$(iii)：$X$ は $P_q=f(Q)$ と可換である（定理 1.5(d)）から、$q\neq q'$ に対し $P_qXP_{q'}=XP_qP_{q'}=0$ である。
(iii)$\Rightarrow$(i)：$XU_\theta=\sum_qe^{iq\theta}P_qXP_q=U_\theta X$。ブロック対角作用素は（$P_q$ の像に適合した正規直交基底において）
$\bigoplus_qM_{d_q}\otimes1_1$ をなし、その可換子環は $\bigoplus_q1_{d_q}\otimes M_1=\operatorname{span}\{P_q\}$ である：それらと可換な作用素はすべての $P_q$ と可換であるから
ブロック対角であり、各ブロックは $M_{d_q}$ と可換であるからスカラーである（$m=1$ とした補題 2.5 のステップ5）。
(3) $\frac d{dt}e^{itH}Qe^{-itH}=ie^{itH}[H,Q]e^{-itH}$（命題 8.2）がすべての $t$ で消えることは $[H,Q]=0$ と同値である。
このとき $e^{itH}$ は $P_q$ と可換であり（定理 1.5(d) を二回）、
$\operatorname{Tr}(P_q\rho(t))=\operatorname{Tr}(e^{itH}P_qe^{-itH}\rho)=\operatorname{Tr}(P_q\rho)$ である。$\square$

したがって $U(1)$ の既約表現は指標 $\theta\mapsto e^{iq\theta}$、$q\in\mathbb Z$ であり（任意のアーベル群と同様に一次元）、電荷 $q$ の重複度は $d_q$ である。電荷はテンソル積のもとで足し合わされる：$u_k$ と $w_l$ をそれぞれ固有値 $\alpha_k$ と $\beta_l$ をもつ $Q_1$ と $Q_2$ の正規直交固有基底とすると、
正規直交基底 $u_k\otimes w_l$ は $Q_1\otimes1+1\otimes Q_2$ の固有ベクトル（固有値 $\alpha_k+\beta_l$）からなり、同時に
$e^{i\theta Q_1}\otimes e^{i\theta Q_2}$ の固有ベクトル（固有値 $e^{i\theta\alpha_k}e^{i\theta\beta_l}$）からなる。したがって関数計算の定義により
$$
e^{i\theta(Q_1\otimes1+1\otimes Q_2)}=e^{i\theta Q_1}\otimes e^{i\theta Q_2},
$$
であり、これは $P_{q_1}\otimes P_{q_2}$ の像の上に電荷 $q_1+q_2$ で作用する。これが第7章の融合テンソル積の融合則 $(q_1,q_2)\mapsto q_1+q_2$ である。

$n$ 個の量子ビットに対しては、電荷は **個数作用素** $\hat N=\sum_jn_j$ である。ここで $n_j=|1\rangle\langle1|$ は量子ビット
$j$ 上に作用する：ハミング重み $|b|$ により $\hat N|b\rangle=|b|\,|b\rangle$ である。そのセクターの次元は
$\binom nw$ であるから、ハミング重みを保存する作用素全体は
$$
U(1)'\cong\bigoplus_{w=0}^nM_{\binom nw}(\mathbb C),
$$
をなす。これは `chains.SpinChain(n).algebra` であり、電荷 $0,1$ をもつ $M_1\oplus M_1$ の $n$ 個のコピーの融合積である
（ビット $1$ = 上向きスピン、$\hat N=S^z+\frac n2$。一方 `decompose.total_spin` では上向きスピンは $\mathbb C^2$ の *第一の* 基底ベクトルであり、これは基底の番号を付け替えるがセクターは変えない）。任意の XXZ ハミルトニアン
$\sum_{\langle ij\rangle}J(s^x_is^x_j+s^y_is^y_j+\Delta s^z_is^z_j)$ は $\hat N$ を保存する：$s^z_is^z_j$ は基底 $|b\rangle$ において対角であり、ホッピング項
$s^x_is^x_j+s^y_is^y_j=\frac12(s^+_is^-_j+s^-_is^+_j)$（$s^\pm=s^x\pm is^y$ を展開すると交差項は打ち消し合う）は $|b\rangle$ を $0$ に、あるいは
上向きスピンがサイト $i,j$ の一方から他方へ移った配置に写すので、$\hat N$ の各固有空間をそれ自身に写すからである。

```python
import functools
from torch_vn_algebra import SpinChain

L = 4
Sx, Sy, Sz = dc.total_spin(L)                              # (16, 16) each, site 0 most significant
N_op = Sz + L / 2 * torch.eye(2 ** L, dtype=dt)            # number of up spins
DU = dc.from_generators(N_op)                              # the algebra generated by U(1)
print(DU)
chain = SpinChain(L, boundary='periodic')                  # U(1)' = (+)_w M_binom(L, w)
assert set(DU.k_factors) == {1} and sorted(DU.multiplicities) == sorted(chain.algebra.k_factors)
assert sorted(chain.algebra.k_factors) == sorted(math.comb(L, w) for w in range(L + 1))

sx, sy, sz = dc.spin_matrices(0.5)
def two_site(a, b, i, j, n=L):                             # a on site i, b on site j of n sites
    ops = [a if k == i else b if k == j else torch.eye(2, dtype=dt) for k in range(n)]
    return functools.reduce(torch.kron, ops)

H = sum(two_site(sx, sx, i, j) + two_site(sy, sy, i, j) + 0.5 * two_site(sz, sz, i, j)
        for i, j in chain.bonds)                           # XXZ ring, Delta = 0.5
assert close(H @ N_op, N_op @ H)                           # [H, N] = 0 ...
assert close(DU.from_commutant(DU.to_commutant(H)), H)     # ... so H is block diagonal (Proposition 4.10)
Hc = chain.xxz(J=1.0, Delta=0.5)                           # the same H, built sector by sector
E_sec = torch.cat([torch.linalg.eigvalsh(Hc.matrix[0, N, :k, :k]) for N, k in enumerate(chain.algebra.k_factors)])
assert close(torch.sort(E_sec).values, torch.linalg.eigvalsh(H))
```

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. $U(1)$ 対称性とは、粒子数、電荷、あるいは磁化の保存である。命題 4.10(3) はこの設定におけるネーターの定理である：電荷が保存されるのはハミルトニアンが対称なときかつそのときに限り、そのときダイナミクスはセクターを決して離れない（定理 8.3）。さらにセクター間のコヒーレンスを測定できるものが何もなければ、セクターは超選択セクターである（第3章）。

## 4.5 $SU(2)$：スピン

**定義.** $\mathbb C^d$ 上の **スピン表現** とは、
$$
[S_x,S_y]=iS_z,\qquad[S_y,S_z]=iS_x,\qquad[S_z,S_x]=iS_y,
$$
を満たすエルミート行列の三つ組 $S_x,S_y,S_z$ のこと、すなわちこれらの関係をもつラベル $\{x,y,z\}$ の表現（4.1節）のことである。基本的な例は
$\mathbb C^2$ 上の $s_a=\frac12\sigma_a$（パウリ行列）である。$(S_a)$ と $(T_a)$ がスピン表現ならば
$S_a\otimes1+1\otimes T_a$ もそうである（交差項は可換）。したがって $L$ 個のスピン $\frac12$ の **全スピン** $S_a=\sum_is_a^{(i)}$
（`decompose.total_spin(L)`）はスピン表現である。回転 $e^{-i\theta\,\mathbf n\cdot\mathbf S}$（$\theta\in\mathbb R$、
$\mathbf n\in\mathbb R^3$ は単位ベクトル、関数計算）がすべての $\theta,\mathbf n$ について
$X$ と可換であることは、$X$ が三つの $S_a$ と可換であることと同値であり（$\theta=0$ で微分する、命題 8.2(1)。逆は定理 1.5(d)）、したがって「$SU(2)$ 作用」の可換子環は $\{S_x,S_y,S_z\}'$ であって、定理 2.2 により回転と $S_a$ は同じ環を生成する（どちらの集合も随伴をとる操作で閉じている）。$L$ 個のスピン $\frac12$ に対し、回転とは $V^{\otimes L}$、$V\in SU(2)$ のことである。実際、$h=\mathbf n\cdot\mathbf s$ として $\mathbf n\cdot\mathbf S=\sum_ih^{(i)}$ であり、4.4節の電荷の場合と同様に
（$L$ に関する帰納法）$e^{-i\theta\,\mathbf n\cdot\mathbf S}=(e^{-i\theta h})^{\otimes L}$ で、$e^{-i\theta h}\in SU(2)$ である（$h$ はトレースが零で $h^2=\frac14$ であるから、その固有値は
$e^{\mp i\theta/2}$ である）。逆に、$V\in SU(2)$ は絶対値 $1$ で積が $1$ の固有値をもつ正規行列であるから、ユニタリ $W$ により
$V=W\operatorname{diag}(e^{i\varphi},e^{-i\varphi})W^*$ と書け、トレースが零のエルミート行列 $A=W\operatorname{diag}(\varphi,-\varphi)W^*$ により $V=e^{iA}$ となる。トレースが零のエルミートな $2\times2$ 行列はパウリ行列の実線形結合
$A=2\,\mathbf c\cdot\mathbf s$ であるから、$\theta=2|\mathbf c|$、$\mathbf n=-\mathbf c/|\mathbf c|$ として $V=e^{-i\theta\,\mathbf n\cdot\mathbf s}$ である
（$\mathbf c=0$ ならば $\theta=0$）。定理 B.6(a),(c) と系 B.7（$d=2$）は、$V^{\otimes L}$、$V\in SU(2)$ が $S_a$ の生成する環を張ることを直接示す。

$S_\pm=S_x\pm iS_y$ とおくと $S_-=S_+^*$ である。また $\mathbf S^2=S_x^2+S_y^2+S_z^2$ とおく。関係式から
$$
[S_z,S_\pm]=\pm S_\pm,\qquad[S_+,S_-]=2S_z,\qquad\mathbf S^2=S_-S_++S_z^2+S_z,
$$
が得られ、$\mathbf S^2$ はすべての $S_a$ と可換である（たとえば
$[\mathbf S^2,S_z]=S_x[S_x,S_z]+[S_x,S_z]S_x+S_y[S_y,S_z]+[S_y,S_z]S_y=-i\{S_x,S_y\}+i\{S_y,S_x\}=0$）。

**定理 4.11（既約スピン表現）.** 任意の $s\in\{0,\frac12,1,\frac32,\dots\}$ に対し、
`decompose.spin_matrices(s)` の行列——基底 $|m\rangle$、$m=s,s-1,\dots,-s$ をもつ $\mathbb C^{2s+1}$ 上の
$$
S_z|m\rangle=m|m\rangle,\qquad S_+|m\rangle=\sqrt{s(s+1)-m(m+1)}\;|m+1\rangle,\qquad S_-=S_+^*,
$$
——は $\mathbf S^2=s(s+1)1$ を満たす既約なスピン表現（**スピン $s$**）をなす。任意の既約スピン表現は、これらのうちちょうど一つとユニタリ同値である。

*証明.* *スピン $s$ はスピン表現である.* $[S_z,S_+]=S_+$ は明らかである。$S_+S_-|m\rangle=(s(s+1)-(m-1)m)|m\rangle$
かつ $S_-S_+|m\rangle=(s(s+1)-m(m+1))|m\rangle$ であるから $[S_+,S_-]=2S_z$ である。これら二つの関係とその随伴は、三つの定義関係と同値である。$\mathbf S^2|m\rangle=(s(s+1)-m(m+1)+m^2+m)|m\rangle=s(s+1)|m\rangle$。
*既約性：* 不変な $W\neq0$ は $S_z$ で不変であり、$S_z$ の固有値は単純であるから、$W$ はある
$|m\rangle$ を含む。$S_\pm$ の係数は $m=\pm s$ でしか消えないので、$S_\pm$ を作用させればすべての $|m'\rangle$ が得られる。

*分類.* $(S_a)$ を $\mathbb C^d$ 上で既約とし、$\lambda$ を $S_z$ の最大固有値、
$v_0$ を単位固有ベクトルとする。$S_zS_+v_0=S_+(S_z+1)v_0=(\lambda+1)S_+v_0$ であるから、最大性により $S_+v_0=0$ である。
$v_j=S_-^jv_0$ とおく。すると $S_zv_j=(\lambda-j)v_j$ であり、$S_+S_-=S_-S_++2S_z$ を用いた帰納法により
$$
S_+v_j=j(2\lambda-j+1)\,v_{j-1}
$$
である（$S_+v_j=S_-S_+v_{j-1}+2S_zv_{j-1}=[(j-1)(2\lambda-j+2)+2(\lambda-j+1)]v_{j-1}$）。零でない $v_j$ は相異なる固有値に属する
$S_z$ の固有ベクトルであるから、最後の零でないもの $v_r$ が存在する。$v_0,\dots,v_r$ の張る空間は
$S_z,S_\pm$ で不変であるから $\mathbb C^d$ に等しく、$d=r+1$ である。$0=S_+v_{r+1}=(r+1)(2\lambda-r)v_r$ から
$\lambda=r/2=:s$ を得る。最後に $\|v_j\|^2=\langle v_{j-1},S_+v_j\rangle=j(2s-j+1)\|v_{j-1}\|^2$ であるから、単位ベクトル
$|s-j\rangle=v_j/\|v_j\|$ は $m=s-j$ として $S_+|m\rangle=\sqrt{(s-m)(s+m+1)}\,|m+1\rangle$ を満たし、これは上の公式である。異なる $s$ は異なる次元をもつ。$\square$

**系 4.12.** $\mathbb C^n$ 上の任意のスピン表現は
$\bigoplus_s(\text{spin }s)\otimes1_{m_s}$ にユニタリ同値であり、$d_M=\dim\ker(S_z-M)$ とおくと
$$
m_s=d_s-d_{s+1}.
$$
である。

*証明.* 定理 4.3 と定理 4.11 による。スピン $s'$ の各コピーは、$|M|\le s'$ かつ $s'-M\in\mathbb Z$ ならば固有値 $M$ に属する $S_z$ の固有ベクトルをちょうど一つ含み、そうでなければ含まない。よって
$d_M=\sum_{s'\ge|M|,\,s'-M\in\mathbb Z}m_{s'}$ であり、連続する二つの値の差が $m_s$ である。$\square$

**命題 4.13（$L$ 個のスピン $\frac12$）.** $(\mathbb C^2)^{\otimes L}$ 上の全スピンがスピン $S$ を含むのは
$S\in\{\frac L2,\frac L2-1,\dots\}$、$S\ge0$ のときかつそのときに限り、その重複度は
$$
m_S=\binom L{\frac L2-S}-\binom L{\frac L2-S-1}\qquad\Big(\tbinom L{-1}=0\Big),
$$
である。したがって $S_x,S_y,S_z$ が生成する環は $\bigoplus_SM_{2S+1}\otimes1_{m_S}$ であり、その可換子環は
$\bigoplus_S1_{2S+1}\otimes M_{m_S}$ である。

*証明.* $S_z=\sum_is_z^{(i)}$ は積基底において対角であり、$u$ 個の上向きスピンと $L-u$ 個の下向きスピンの積の上で固有値 $u-\frac L2$ をもつから、$\frac L2+M\in\{0,\dots,L\}$ ならば $d_M=\binom L{L/2+M}=\binom L{L/2-M}$、そうでなければ $0$ である。系 4.12 と定理 4.3 を適用すればよい。値が正であるのは、
$j<\frac{L-1}2$ に対して $\binom L{j+1}/\binom Lj=\frac{L-j}{j+1}>1$ だからである。$\square$

$L=2$ の場合：$\frac12\otimes\frac12=1\oplus0$（**クレブシュ・ゴルダン**）。三重項 $S_-^j|{\uparrow\uparrow}\rangle$、$j=0,1,2$ は
二つのスピンの入れ替え $F$ のもとで対称であり、一重項
$\frac1{\sqrt2}(|{\uparrow\downarrow}\rangle-|{\downarrow\uparrow}\rangle)$ は反対称で、$S_+$ と $S_z$ によって零化される。
$L=4$ の場合：$m_2=1$、$m_1=4-1=3$、$m_0=6-4=2$ であるから
$$
\{S_x,S_y,S_z\}''\cong M_5\oplus M_3\oplus M_1\ \text{ with multiplicities }1,3,2,\qquad
\{S_x,S_y,S_z\}'\cong M_1\oplus M_3\oplus M_2,
$$
であり、次元は $35$ と $14$ である——これは第10章で見いだされた数である。ライブラリは三つの生成元だけから同じものを見つける（セクターは $(k_c,m_c)$ の順に並べられる）：

```python
S4 = dc.total_spin(4)
assert close(S4[0] @ S4[1] - S4[1] @ S4[0], 1j * S4[2])                # [S_x, S_y] = i S_z
D4 = dc.from_generators(S4)
print(D4)
assert D4.k_factors == (1, 3, 5) and D4.multiplicities == (2, 3, 1)     # M_1, M_3, M_5 with m = 2, 3, 1
assert D4.verify()['ok'] and D4.dims == (35, 14)
m_S = {S: math.comb(4, 2 - S) - (math.comb(4, 1 - S) if S < 2 else 0) for S in (0, 1, 2)}
assert [m_S[(k - 1) // 2] for k in D4.k_factors] == list(D4.multiplicities)    # Proposition 4.13
casimir = D4.to_algebra((S4 @ S4).sum(0))                                  # S^2 = S(S+1) on spin S
for c, k in enumerate(D4.k_factors):
    assert close(casimir.matrix[0, c, :k, :k], (k - 1) / 2 * (k + 1) / 2 * torch.eye(k))
```

どの作用素がすべての回転と可換だろうか。サイトの置換はそうである：$\pi(\sigma)$ による共役は
$s^{(i)}_a$ を $s^{(\sigma(i))}_a$ に写す（補題 B.1(b)）ので、$\sum_is^{(i)}_a$ を固定する。逆に、それ以外には何もない：

**定理 4.14（シューア・ワイル双対性）.** $\pi$ を、$(\mathbb C^2)^{\otimes L}$ のテンソル因子を置換する $S_L$ の表現とする。(a) $\operatorname{span}\pi(S_L)\subseteq\{S_x,S_y,S_z\}'$。(b) 等号が成り立つ。同値なことに（定理 2.2）、全スピンが生成する環は $\pi(S_L)'$ である。

*証明.* (a) は今証明したところである。(b) は系 B.7、すなわち定理 B.6 の $d=2$ の場合である。定理 B.6 はすべての $d$ と $L$ について、
$(\mathbb C^d)^{\otimes L}$ と $U(d)$ の対角作用に対する双対性を証明している。$\square$

（さらに読むには：W. Fulton and J. Harris, *Representation Theory*, Springer 1991, Lecture 6; R. Goodman and
N. Wallach, *Symmetry, Representations, and Invariants*, Springer 2009, Section 9.1。）$L\le5$ について、コードは二つの空間を数値的に比較する：
(a) により、それらの次元が等しければ十分である。

```python
for n in range(2, 6):
    comm = dc.commutant_basis(dc.total_spin(n)).reshape(-1, 4 ** n)          # a basis of {S_a}'
    perm = dc.permutation_representation(n, 2).reshape(-1, 4 ** n)          # all n! permutations
    dims = [torch.linalg.matrix_rank(A).item() for A in (comm, perm, torch.cat([comm, perm]))]
    print(n, dims)
    assert dims[0] == dims[1] == dims[2]                                     # span pi(S_n) = {S_a}'
Dp = dc.from_group(dc.permutation_representation(4, 2))                   # pi(S_4)'' = {S_a}'
print(Dp)
assert sorted(zip(Dp.multiplicities, Dp.k_factors)) == sorted(zip(D4.k_factors, D4.multiplicities))
```

出力される次元は $2,5,14,42$ であり、カタラン数 $\binom{2L}L/(L+1)$ である。系 B.10 はすべての $L$ について
$\dim\{S_a\}'=\sum_Sm_S^2=\binom{2L}L/(L+1)$ を証明している。最後の行は双対性が働いている様子である：
$\pi(S_4)''$ はブロック $k=(1,2,3)$ と重複度 $(5,1,3)$ をもつ——全スピンの組 $(k_c,m_c)$ の役割を入れ替えたものである。一般に（系 B.9）スピン $S$ の重複度空間は次元 $m_S$ の $S_L$ の既約表現を担い、これらは互いに非同値である。ここでは $S=2,1,0$ に対してその次元は $1,3,2$ である。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. スピン $S$ の多重項は、任意の回転不変なハミルトニアンに対して等しいエネルギーをもつ $2S+1$ 個の状態からなり、それらは $S^z=M$ によって区別される。二つのスピン $\frac12$ は三重項と一重項をなし、
> 四つのスピンは一つの五重項、三つの三重項、二つの一重項をなす。二つの一重項の重複度空間と三つの三重項の重複度空間が、第10章のデコヒーレンスフリーな量子ビットとノイズレスなクォートリットである。

## 4.6 スピン系

$L$ 個のスピン $\frac12$ 上の **ハイゼンベルク模型** は、実数の結合定数 $J_{ij}$ により $H=\sum_{i<j}J_{ij}\,\mathbf s^{(i)}\cdot\mathbf s^{(j)}$ で与えられる。

**命題 4.15.** (a)（*ディラックの恒等式*）$\mathbf s^{(i)}\cdot\mathbf s^{(j)}=\frac12F_{ij}-\frac14$。ここで $F_{ij}$
はスピン $i$ と $j$ を入れ替える。したがって任意のハイゼンベルク・ハミルトニアンは $\operatorname{span}\pi(S_L)\subseteq\{S_a\}'$ に属し、
命題 4.13 の分解において $V^*HV=\bigoplus_S1_{2S+1}\otimes H_S$（ブロック $H_S\in M_{m_S}$）となる：$H_S$ の各固有値は $H$ の $(2S+1)$ 重に縮退した準位である。
(b) **最高ウェイトベクトル** の空間 $\ker S_+\cap\ker(S_z-S)$ は次元 $m_S$ をもち、$H$ で不変であり、
$H$ はその上に（基底のユニタリな変換を除いて）$H_S$ として作用する。

*証明.* (a) 二つのスピン上では $\mathbf s^{(1)}\cdot\mathbf s^{(2)}=\frac12(\mathbf S^2-\mathbf s^{(1)2}-\mathbf s^{(2)2})=\frac12\mathbf S^2-\frac34$
である（$s=\frac12$ の定理 4.11）。三重項のベクトル $S_-^j|{\uparrow\uparrow}\rangle$ は $|{\uparrow\uparrow}\rangle$、
$|{\downarrow\uparrow}\rangle+|{\uparrow\downarrow}\rangle$、$2|{\downarrow\downarrow}\rangle$ であり、一重項とあわせて $\mathbb C^2\otimes\mathbb C^2$ の直交基底をなす。
三重項上では $\mathbf S^2=2$ かつ $F=1$ である：$\mathbf S^2=S_-S_++S_z^2+S_z$ は $|{\uparrow\uparrow}\rangle$ 上で $0+1+1=2$ であり、
$\mathbf S^2$ と $F$ は $S_-$ と可換で、$F$ は $|{\uparrow\uparrow}\rangle$ を固定する。一重項上では
$\mathbf S^2=0$（$S_+$ と $S_z$ によって零化される）かつ $F=-1$ である。両辺はこれら四つのベクトル上で一致する。$L$ 個のうちのサイト $i,j$ については恒等作用素とのテンソル積をとる。残りは定理 4.14(a) と定理 4.3 である。(b) より多くを証明する：系 4.12 のように分解された任意のスピン表現と任意の $\mu\in\frac12\mathbb Z$ に対し、
$\ker S_+\cap\ker(S_z-\mu)$ は、$\mu\ge0$ ならば $|\mu,\mu\rangle\otimes\mathbb C^{m_\mu}$（の系 4.12 のユニタリによる像）であり、$\mu<0$ ならば $0$ である。実際、その基底において $S_z=\mu$ をもつベクトルは、
$S'-\mu\in\mathbb Z$ を満たす $S'\ge|\mu|$ にわたる和 $\sum_{S'}|S',\mu\rangle\otimes w_{S'}$ であり、定理 4.11 により $S_+$ はそれを
$\sum_{S'}\sqrt{(S'-\mu)(S'+\mu+1)}\,|S',\mu+1\rangle\otimes w_{S'}$ に写す。ここで $S'+\mu+1\ge1$ である。よってそれが消えるのは、
$S'\neq\mu$ なるすべての $S'$ に対して $w_{S'}=0$ となるときかつそのときに限る。ハイゼンベルクの $H$ に対し、$1\otimes H_S$ は $|S,S\rangle\otimes\mathbb C^{m_S}$ 上に $H_S$ として作用する。$\square$

(b) は、モジュール `su2`（スピン鎖、第9章）が数値的な分解なしに $H_S$ を構成する方法である：それは一度に一サイトずつ加えながら、
$j\otimes\frac12$ のクレブシュ・ゴルダン係数を用いて最高ウェイトベクトルの正規直交基底を構成する。これらの係数を最高ウェイトからの下降によって導こう。定理 4.11
と $S_-=S_+^*$ により、スピン $j$ は
$$
S_\pm|j,m\rangle=\sqrt{(j\mp m)(j\pm m+1)}\;|j,m\pm1\rangle ,
$$
を満たす。$|m|>j$ に対しては $|j,m\rangle=0$ とおく（これは公式と整合的である。その係数は $m=\pm j$ で消える）。

**定理 4.16（$j\otimes\frac12$ のクレブシュ・ゴルダン係数）.** $j\in\{\frac12,1,\frac32,\dots\}$ とし、
$\mathbb C^{2j+1}\otimes\mathbb C^2$ 上で $S_a=S^{(j)}_a\otimes1+1\otimes s_a$ とする（基底 $|j,m\rangle$ をもつ定理 4.11 のスピン $j$ を、
$|{\uparrow}\rangle=|\frac12,\frac12\rangle$、$|{\downarrow}\rangle=|\frac12,-\frac12\rangle$ をもつスピン $\frac12$ と結合したもの）。また
$\alpha_M=\sqrt{(j+M+\frac12)/(2j+1)}$、$\beta_M=\sqrt{(j-M+\frac12)/(2j+1)}$ とおく。$J=j\pm\frac12$ と $M=J,J-1,\dots,-J$ に対し
$$
\begin{aligned}
|j+\tfrac12,M\rangle\!\rangle&=\alpha_M\,|j,M-\tfrac12\rangle\otimes|{\uparrow}\rangle+\beta_M\,|j,M+\tfrac12\rangle\otimes|{\downarrow}\rangle,\\
|j-\tfrac12,M\rangle\!\rangle&=-\beta_M\,|j,M-\tfrac12\rangle\otimes|{\uparrow}\rangle+\alpha_M\,|j,M+\tfrac12\rangle\otimes|{\downarrow}\rangle .
\end{aligned}
$$
とおく。
(a) これら $2(2j+1)$ 個のベクトルは $\mathbb C^{2j+1}\otimes\mathbb C^2$ の正規直交基底である。
(b) $|j+\frac12,j+\frac12\rangle\!\rangle=|j,j\rangle\otimes|{\uparrow}\rangle$ であり、$|j-\frac12,j-\frac12\rangle\!\rangle$ は、絶対値 $1$ の因子を除いて、
$S_zx=(j-\frac12)x$ かつ $S_+x=0$ を満たす唯一の単位ベクトル $x$ である。この因子は
$|j,j\rangle\otimes|{\downarrow}\rangle$ の係数が正であることによって定まる（コンドン・ショートレーの規約）。
(c) $S_z|J,M\rangle\!\rangle=M|J,M\rangle\!\rangle$ かつ $S_-|J,M\rangle\!\rangle=\sqrt{(J+M)(J-M+1)}\,|J,M-1\rangle\!\rangle$ である：状態 $|J,M\rangle\!\rangle$
は $|J,J\rangle\!\rangle$ から下降と規格化によって得られる。
(d) したがって $a=x,y,z$ に対し $S_a|J,M\rangle\!\rangle=\sum_{M'}(S^{(J)}_a)_{M'M}|J,M'\rangle\!\rangle$ である：ユニタリ $|J,M\rangle\mapsto|J,M\rangle\!\rangle$ は
スピン $(j+\frac12)\oplus$ スピン $(j-\frac12)$ を $S_a$ に写すので、$j\otimes\frac12\cong(j+\frac12)\oplus(j-\frac12)$ であり、
クレブシュ・ゴルダン係数 $\langle j,m;\frac12,\sigma|J,M\rangle$（積基底における $|J,M\rangle\!\rangle$ の座標）は上に示した
$\pm\alpha_M,\beta_M$ である。

*証明.* (a) 上の式の基底ベクトルが $0$ ならば、その係数も $0$ である：$|j,M-\frac12\rangle=0$ となるのは
$M=-j-\frac12$ のときだけで、このとき $\alpha_M=0$ であり、$|j,M+\frac12\rangle=0$ となるのは $M=j+\frac12$ のときだけで、このとき $\beta_M=0$ である。よって各ベクトルのノルムの二乗は
$\alpha_M^2+\beta_M^2=1$ である。$M$ が異なるベクトルは異なる積基底ベクトルを含み、$M$ が等しいときの
内積は $-\alpha_M\beta_M+\beta_M\alpha_M=0$ である。それらは $(2j+2)+2j=2(2j+1)$ 個ある。
(b) 第一の主張は $\alpha_{j+1/2}=1$、$\beta_{j+1/2}=0$ である。固有空間 $S_z=j-\frac12$ は
$|j,j-1\rangle\otimes|{\uparrow}\rangle$ と $|j,j\rangle\otimes|{\downarrow}\rangle$ で張られ、$s_+|{\uparrow}\rangle=0$、$s_+|{\downarrow}\rangle=|{\uparrow}\rangle$ により
$$
S_+\big(p\,|j,j-1\rangle\otimes|{\uparrow}\rangle+q\,|j,j\rangle\otimes|{\downarrow}\rangle\big)=(p\sqrt{2j}+q)\,|j,j\rangle\otimes|{\uparrow}\rangle ,
$$
である。これが消えるのは $q=-\sqrt{2j}\,p$ のときかつそのときに限る。$|p|^2+|q|^2=1$ と $q>0$（したがって $p<0$）のもとでは、これは $p=-\beta_{j-1/2}=-1/\sqrt{2j+1}$、
$q=\alpha_{j-1/2}=\sqrt{2j/(2j+1)}$ である。(c) $S_z$：どちらの項も $S_z=(M\mp\frac12)\pm\frac12=M$ をもつ。$S_-$ については、
$x_M=p\,|j,M-\frac12\rangle\otimes|{\uparrow}\rangle+q\,|j,M+\frac12\rangle\otimes|{\downarrow}\rangle$ とおく。$s_-|{\uparrow}\rangle=|{\downarrow}\rangle$、$s_-|{\downarrow}\rangle=0$ であるから
$$
S_-x_M=p\sqrt{(j+M-\tfrac12)(j-M+\tfrac32)}\;|j,M-\tfrac32\rangle\otimes|{\uparrow}\rangle
+\Big(p+q\sqrt{(j+M+\tfrac12)(j-M+\tfrac12)}\Big)|j,M-\tfrac12\rangle\otimes|{\downarrow}\rangle .
$$
である。$\beta_M\sqrt{(j+M+\frac12)(j-M+\frac12)}=(j-M+\frac12)\alpha_M$ と $\alpha_M\sqrt{(j+M+\frac12)(j-M+\frac12)}=(j+M+\frac12)\beta_M$ を用いる。
$J=j+\frac12$ に対しては $(p,q)=(\alpha_M,\beta_M)$、$c=\sqrt{(J+M)(J-M+1)}=\sqrt{(j+M+\frac12)(j-M+\frac32)}$ であり、第一の係数は
$\sqrt{(j+M+\frac12)(j+M-\frac12)(j-M+\frac32)/(2j+1)}=c\,\alpha_{M-1}$、第二の係数は $\alpha_M(j-M+\frac32)=c\,\beta_{M-1}$ である。
$J=j-\frac12$ に対しては $(p,q)=(-\beta_M,\alpha_M)$、$c=\sqrt{(j+M-\frac12)(j-M+\frac12)}$ であり、第一の係数は
$-\sqrt{(j-M+\frac12)(j+M-\frac12)(j-M+\frac32)/(2j+1)}=-c\,\beta_{M-1}$、第二の係数は $\beta_M(j+M-\frac12)=c\,\alpha_{M-1}$ である。
いずれの場合も $M>-J$ に対して $S_-|J,M\rangle\!\rangle=c\,|J,M-1\rangle\!\rangle$ である。$M=-J$ に対しては両方の係数が消える（$J=j+\frac12$ ならば $p=\alpha_M=0$、
$J=j-\frac12$ ならば $j+M-\frac12=0$）ので、$S_-|J,-J\rangle\!\rangle=0$ である。(d) (a) と (c) により、正規直交基底
$|J,M\rangle\!\rangle$ における $S_z$ と $S_-$ の行列はスピン $(j+\frac12)\oplus$ スピン $(j-\frac12)$ のもの（定理 4.11
と上の公式）である。$S_+=S_-^*$ の行列も同様であり、したがって $S_x=\frac12(S_++S_-)$ と $S_y=\frac1{2i}(S_+-S_-)$ の行列も同様である。$\square$

演習問題 2 はこの分解を（係数は除いて）$j_1\otimes j_2$ に拡張する。ライブラリにとって重要なのは定理 4.16 の二つの
最高ウェイトベクトルだけであり、それらは既約であるか否かを問わず任意のスピン表現において機能する：

**系 4.17（一度に一サイトずつの最高ウェイトベクトル）.** $(T_a)$ を $E$ 上のスピン表現とし、
$\mathcal H_\mu(E)=\ker T_+\cap\ker(T_z-\mu)$ とおき、その正規直交基底を $B_\mu$ とする（$\mathcal H_\mu(E)=0$ ならば $B_\mu=\emptyset$）。また
$E\otimes\mathbb C^2$ 上で $S_a=T_a\otimes1+1\otimes s_a$ とする。任意の $S\in\{0,\frac12,1,\dots\}$ に対し、ベクトル
$$
v\otimes|{\uparrow}\rangle\quad(v\in B_{S-1/2}),\qquad
\sqrt{\tfrac{2S+1}{2S+2}}\;v\otimes|{\downarrow}\rangle-\tfrac1{\sqrt{(2S+1)(2S+2)}}\;T_-v\otimes|{\uparrow}\rangle\quad(v\in B_{S+1/2})
$$
は $\ker S_+\cap\ker(S_z-S)$ の正規直交基底をなす（$S=0$ では第一の族は空である）。一つのスピン
（$E=\mathbb C^2$、$B_{1/2}=\{|{\uparrow}\rangle\}$、他のすべての $B_\mu$ は空）から始めて $L-1$ 個のスピンを加えていくと、
$L$ 個のスピン $\frac12$ の最高ウェイトベクトルの空間の、次元 $m_S$ の正規直交基底が得られる（命題 4.15(b)）。

*証明.* $a=\sqrt{(2S+1)/(2S+2)}$、$b=-1/\sqrt{(2S+1)(2S+2)}$、$u(v)=a\,v\otimes|{\downarrow}\rangle+b\,T_-v\otimes|{\uparrow}\rangle$ と書く。
$v\in\mathcal H_{S+1/2}(E)$ に対し $T_+T_-v=(T_-T_++2T_z)v=(2S+1)v$ であり、$a+(2S+1)b=0$ である。*最高ウェイト：* $S_z$ はどちらの族の上でも $S$ である
（$v$ 上で $T_z=S\mp\frac12$、$T_-v$ 上で $S-\frac12$）。$S_+(v\otimes|{\uparrow}\rangle)=T_+v\otimes|{\uparrow}\rangle=0$ であり、
$S_+u(v)=a\,T_+v\otimes|{\downarrow}\rangle+a\,v\otimes|{\uparrow}\rangle+b\,T_+T_-v\otimes|{\uparrow}\rangle=(a+(2S+1)b)\,v\otimes|{\uparrow}\rangle=0$ である。
*正規直交性：* $v\in B_{S-1/2}$ に対し $\langle v\otimes|{\uparrow}\rangle,u(v')\rangle=b\langle v,T_-v'\rangle=b\langle T_+v,v'\rangle=0$ であり、
$\langle u(v),u(v')\rangle=a^2\langle v,v'\rangle+b^2\langle v,T_+T_-v'\rangle=(a^2+(2S+1)b^2)\langle v,v'\rangle=\langle v,v'\rangle$ である。
*張ること：* $S_zw=Sw$ かつ $S_+w=0$ とし、$w=x\otimes|{\uparrow}\rangle+y\otimes|{\downarrow}\rangle$ と書く。すると $T_zx=(S-\frac12)x$、
$T_zy=(S+\frac12)y$ であり、$0=S_+w=(T_+x+y)\otimes|{\uparrow}\rangle+T_+y\otimes|{\downarrow}\rangle$ であるから $T_+y=0$、すなわち $y\in\mathcal H_{S+1/2}(E)$
で、$T_+x=-y$ である。すると $x'=x+\frac1{2S+1}T_-y$ は $T_+x'=-y+y=0$ を満たすので $x'\in\mathcal H_{S-1/2}(E)$ であり、これは $S=0$ では $0$ である
（命題 4.15(b) の証明）。$b/a=-1/(2S+1)$ であるから、$w=x'\otimes|{\uparrow}\rangle+a^{-1}u(y)$（$u$ を線形に拡張する）は
二つの族の張る空間に属する。$L$ 個のスピンについては、$n=2,\dots,L$ に対して $E=(\mathbb C^2)^{\otimes(n-1)}$ とその全スピン $T_a$ にこれを適用する。
一つのスピン上では $\mathcal H_{1/2}=\mathbb C|{\uparrow}\rangle$ であり、他のすべての $\mathcal H_\mu=0$ である（スピンは $\frac12$ だけである）。$\square$

第二の族は、$j=S+\frac12$ とした定理 4.16 のベクトル $|S,S\rangle\!\rangle$ を、$T_-^kv$ の張るスピン $j$ のコピーの中で書いたものであり
（定理 4.11 の証明、$|j,j-1\rangle=T_-v/\sqrt{2j}$）、第一の族は $j=S-\frac12$ とした
$|j+\frac12,j+\frac12\rangle\!\rangle$ である。各基底ベクトルに中間スピンの列
$S_1=\frac12,S_2,\dots,S_L=S$、$|S_n-S_{n-1}|=\frac12$ でラベルをつけると、`su2.SU2Chain` が構成するヤマノウチ基底が得られる（そのソースは
まさにこれらの係数を用いる。下のコードは系 4.17 の結論だけを確かめる）。それは `SpinChain` の $U(1)$ セクター $\hat N=\frac L2+S$ の内部にある（`SpinChain` のサイトの順序は
テンソル因子の付け替えであり、補題 B.1(b) により全スピンと可換である）。コードは定理 4.16 を
密な行列とライブラリの `su2.clebsch_gordan` に対して確かめ、それから八つのスピンの基底について系 4.17 の結論を確かめる。また命題 4.15(a) を `D4` を用いて四つのスピンの環で確かめ、
八つのスピンについてブロックの大きさを比較する。

```python
def swap(i, j, n):
    p = list(range(n))
    p[i], p[j] = j, i
    return dc.permutation_representation(n, 2, [p])[0]

ring = [(i, (i + 1) % 4) for i in range(4)]
H4 = sum(swap(i, j, 4) / 2 - torch.eye(16, dtype=dt) / 4 for i, j in ring)        # Dirac's identity
assert close(H4, sum(two_site(a, a, i, j) for a in (sx, sy, sz) for i, j in ring))
h = D4.to_commutant(H4)                                   # blocks H_S of sizes m_S = 2, 3, 1
assert close(D4.from_commutant(h), H4)                    # H4 lies in {S_a}'
levels = torch.cat([torch.linalg.eigvalsh(h.matrix[0, c, :m, :m]).repeat_interleave(k)
                    for c, (k, m) in enumerate(zip(D4.k_factors, D4.multiplicities))])
assert close(torch.sort(levels).values, torch.linalg.eigvalsh(H4))   # level of H_S repeated 2S+1 times

from torch_vn_algebra import su2
u1, sym = SpinChain(8, boundary='periodic'), su2.SU2Chain(8, boundary='periodic')
H_u1, H_su2 = u1.xxz(J=1.0, Delta=1.0), sym.heisenberg(J=1.0)
E_u1 = torch.sort(torch.cat([torch.linalg.eigvalsh(H_u1.matrix[0, N, :k, :k])
                             for N, k in enumerate(u1.algebra.k_factors)])).values
E_su2, S_lev, _ = sym.full_spectrum(H_su2, return_quantum_numbers=True)   # level of block S, 2S+1 times
E_su2, S_lev = E_su2[0], S_lev[0]
assert close(E_u1, E_su2)
assert S_lev[0] == 0 and (S_lev[1:4] == 1).all() and E_su2[4] - E_su2[3] > 1e-3     # singlet, then triplet
assert close(sym.restrict(H_u1).matrix, H_su2.matrix)      # Proposition 4.15(b)
flops = lambda ks: sum(k ** 3 for k in ks)
print(u1.algebra.k_factors, sym.multiplicities, flops(u1.algebra.k_factors), flops(sym.multiplicities))
print(f"ground state energy {E_su2[0].item():.6f}, first excited {E_su2[1].item():.6f} (x3)")
```

八つのスピンの環では、$U(1)$ のブロックの大きさは $1,8,28,56,70,56,28,8,1$、$SU(2)$ のブロックの大きさは
$S=0,\dots,4$ に対して $m_S=14,28,20,7,1$ である。対角化のコスト $\sum k^3$ は $739162$ から $33040$ へ、
$22$ 分の1に減る。環の次元 $\sum k^2$ は $\binom{16}8=12870$ から $1430$ に減る。基底状態はエネルギー $-3.651093$ の一重項であり、第一励起準位は $-3.128419$ の三重項である。

> **物理的意味.** ハイゼンベルク模型は、磁性絶縁体における局在電子スピンの交換相互作用を記述する。$J_{ij}>0$ は反平行な隣接スピンを好む。そのスペクトルは
> 多重項からなり、磁場 $h\,S^z$ は各多重項を $2S+1$ 個の等間隔な準位に分裂させる（演習問題 4）。
> 二部格子の二つの副格子の間の結合が $J_{ij}>0$、副格子内の結合が $J_{ij}\le0$ で、
> 副格子の大きさが等しければ、基底状態は一重項である（リーブ・マティスの定理、E. Lieb and D. Mattis,
> *J. Math. Phys.* 3 (1962) 749。証明なしに述べる）。上の八つのスピンの環はその例である。

## 4.7 比較の対象

有限群については、GAP と SageMath が指標表を厳密に計算し、系 4.5 により重複度が得られる。RepLAB（MATLAB；D. Rosset, F. Montealegre-Mora, J.-D. Bancal）は有限群およびコンパクト群の表現を数値的に分解し、これも可換子環のランダムな元を用いる。
SymbolicWedderburn.jl（Julia）は半正定値計画の対称性による簡約のためにウェダーバーン分解を計算する。物理では、QuSpin や他の厳密対角化コードが
$U(1)$ と格子の対称性によってブロック対角化を行い、ITensor と TeNPy はテンソルネットワークにおいてアーベル的な電荷を扱い、
TensorKit.jl は融合木によって $SU(2)$ を扱う。ライブラリの `decompose` は小さな行列（$n$ は数百程度まで）のための密な手法である。任意の行列のリスト（群の元または生成元）を受け付け、結果を証明書によって確認し（命題 4.9）、環とその可換子環を `to_algebra`、`from_algebra`、`to_commutant`、`from_commutant` をもつ `TypeIAlgebra`
オブジェクトとして返すので、第3章、第5–8章の内容がただちに適用できる。スピン $\frac12$ の鎖については、`su2.SU2Chain` が $L=16$ までのすべてのスピンの値について $SU(2)$ のブロックを直接与える。

## 演習問題

1. $(\mathbb C^d)^{\otimes L}$ 上の $S_L$ の表現について、
   $\chi_\pi(\sigma)=d^{\#\text{cycles}(\sigma)}$ かつ $\dim\pi(S_L)'=\binom{d^2+L-1}L$ であることを示せ。$20$（2.4節）と
   $35$（四つのスピンに対する $=\dim\{S_a\}''$）を確かめよ。
2. （*クレブシュ・ゴルダン.*）スピン表現 $S^{(j_1)}_a\otimes1+1\otimes S^{(j_2)}_a$ に対して
   $j_1\otimes j_2\cong\bigoplus_jj$、$j=|j_1-j_2|,|j_1-j_2|+1,\dots,j_1+j_2$ を示せ——これは第7章で触れた $SU(2)$ の非アーベル的な融合則である——。また
   $1\otimes\frac12$ と $1\otimes1$ を `from_generators` で確かめよ。
3. $\theta/2\pi$ を無理数とする。単一の $U_\theta=e^{i\theta Q}$ に対して $\{U_\theta\}'=U(1)'$ であるが、
   $n$ 個の量子ビットに対して $\{e^{i\pi\hat N}\}'\cong M_{2^{n-1}}\oplus M_{2^{n-1}}$ であることを示せ。
4. （*ゼーマン分裂.*）$H$ を $L\ge1$ 個のスピン上のハイゼンベルク・ハミルトニアンとし、$H_h=H+hS_z$ とする。
   $H_h$ は $S_z$ および $\mathbf S^2$ と可換であるが、$h\neq0$ のとき $S_x$ とは可換でないこと、そしてその固有値が
   $H_S$ の固有値 $E$ に対して $E+hM$、$M=-S,\dots,S$ であることを示せ。八つのスピンの環でこれを確かめよ。

<details>
<summary>解答</summary>

1. $d=2$ の場合と同様に、$\pi(\sigma)$ は置換行列であり、$|b_1\dots b_L\rangle$ を固定するのは $b$ が各巡回上で一定であるときかつそのときに限る。そのような $b$ は $d^{\#\text{cycles}}$ 個ある。系 4.5(c) により
   $\dim\pi(S_L)'=\frac1{L!}\sum_\sigma(d^2)^{\#\text{cycles}(\sigma)}=\frac{x(x+1)\cdots(x+L-1)}{L!}$（$x=d^2$）
   であり（第5章、演習問題 2）、これは $\binom{d^2+L-1}L$ である：$\binom63=20$、$\binom74=35$。$35$ が
   $\dim\{S_a\}''$ でもあることは、$L=4$ の定理 4.14 である。

```python
for n, value in ((3, 20), (4, 35)):
    assert dc.commutant_dimension(dc.permutation_representation(n, 2)) == math.comb(n + 3, n) == value
```

2. $j_1\ge j_2$ とする。$S_z$ は $|m_1\rangle\otimes|m_2\rangle$ 上で固有値 $m_1+m_2$ をもつので、
   $M\in j_1+j_2+\mathbb Z$ に対し、$d_M$ は $|M-m_2|\le j_1$ を満たす $m_2\in\{-j_2,\dots,j_2\}$ の個数である：
   $|M|\le j_1-j_2$ ならば $2j_2+1$ 個すべて、$j_1-j_2\le|M|\le j_1+j_2$ ならば $j_1+j_2-|M|+1$ 個（$M\ge0$ に対しては条件は
   $m_2\ge M-j_1$）、それを超えれば $0$ 個である。系 4.12 により、$m_j=d_j-d_{j+1}$ は $j_1-j_2\le j\le j_1+j_2$ に対して $1$、
   $j<j_1-j_2$ に対して $(2j_2+1)-(2j_2+1)=0$ である。

```python
def coupled(a, b):                                        # S^(a) (x) 1 + 1 (x) S^(b)
    Sa, Sb = dc.spin_matrices(a), dc.spin_matrices(b)
    ea, eb = torch.eye(Sa.shape[-1], dtype=dt), torch.eye(Sb.shape[-1], dtype=dt)
    return torch.stack([torch.kron(x, eb) + torch.kron(ea, y) for x, y in zip(Sa, Sb)])
for a, b, ks in ((1.0, 0.5, (2, 4)), (1.0, 1.0, (1, 3, 5))):  # 1 (x) 1/2 = 1/2 (+) 3/2, 1 (x) 1 = 0 (+) 1 (+) 2
    D_cg = dc.from_generators(coupled(a, b))
    assert D_cg.k_factors == ks and set(D_cg.multiplicities) == {1}
```

3. $U_\theta=\sum_qe^{iq\theta}P_q$ であり、$e^{iq\theta}=e^{iq'\theta}$ であることは $(q-q')\theta\in2\pi\mathbb Z$ と同値で、
   $\theta/2\pi$ が無理数ならばこれは $q=q'$ を強いる。よって $U_\theta$ の固有空間は $P_q$ の像であり、作用素が
   正規行列 $U_\theta$ と可換であることは、それが $U_\theta$ の固有空間を保つことと同値である（定理 1.5(d) と、2.2節と同様のラグランジュ補間）。すなわち $X=\sum_qP_qXP_q$ と同値であり、命題 4.10 によりこれは $U(1)'$ である。
   $\theta=\pi$ に対しては、固有値 $(-1)^w$ はハミング重みの偶奇しか見ない：次元
   $2^{n-1}$ の二つの固有空間があり（$n\ge1$ のとき。偶数個および奇数個の部分集合）、可換子環は $M_{2^{n-1}}\oplus M_{2^{n-1}}$ である。

```python
for theta, m in ((1.0, (1, 1, 4, 4, 6)), (math.pi, (8, 8))):
    D_theta = dc.from_generators(torch.linalg.matrix_exp(1j * theta * N_op))
    assert D_theta.multiplicities == m and set(D_theta.k_factors) == {1}
```

4. $S_z$ と $\mathbf S^2$ は $H$（命題 4.15(a)）および $S_z$ と可換であるから、$H_h$ と可換である。また
   $[H_h,S_x]=h[S_z,S_x]=ihS_y\neq0$ である。命題 4.13 の分解において
   $V^*S_zV=\bigoplus_SS^{(S)}_z\otimes1_{m_S}$ かつ $V^*HV=\bigoplus_S1_{2S+1}\otimes H_S$（命題 4.15）であるから、
   スピン $S$ の直和成分上で $V^*H_hV=hS^{(S)}_z\otimes1+1\otimes H_S$ となる。これは可換な作用素の和であり、
   固有ベクトルの積によって対角化され、固有値は $hM+E$ である。

```python
assert not close(H4 @ Sx - Sx @ H4 + 0.3 * (Sz @ Sx - Sx @ Sz), 0)          # [H4 + 0.3 S_z, S_x] != 0
E_field = sym.full_spectrum(sym.heisenberg(J=1.0, h=0.3))[0]                # E + h M on every multiplet
H_field = u1.xxz(J=1.0, Delta=1.0, h=0.3)
E_u1_field = torch.cat([torch.linalg.eigvalsh(H_field.matrix[0, N, :k, :k])
                        for N, k in enumerate(u1.algebra.k_factors)])
assert close(torch.sort(E_u1_field).values, E_field)
```

</details>

**次へ：**[第5章](05_random_matrices.md) では環の中でランダム行列を抽出する：$\prod_cU(k_c)$ のハール・ユニタリ、ワインガルテン計算（その2次の場合は $\{V\otimes V\}$ の可換子環、すなわちシューア・ワイル双対性の最小の例に基づく）、そしてセクターごとの準位統計である。
