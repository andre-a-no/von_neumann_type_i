[English](../04_symmetries.md) · [Español](../es/04_symmetries.md) · [Français](../fr/04_symmetries.md) · [Deutsch](../de/04_symmetries.md) · [中文](../zh/04_symmetries.md) · **日本語** · [Русский](../ru/04_symmetries.md)

# 4. 対称性と表現

第2章では、対称性の群がフォン・ノイマン環（対称性作用素の張る空間）とその可換子環（対称なオブザーバブル）を生み出すことを示し（命題 2.7）、一つの例を手で分解した。本章ではその背後にある表現論を与え、それをライブラリが任意の行列のリストに対して実行するアルゴリズムに変え、スピン系の二つの対称性に応用する：$`U(1)`$（粒子数あるいは磁化）と $`SU(2)`$（スピンの回転）である。

**本章で学ぶこと**：完全可約性とシューアの補題、既約表現への分解が可換子環と生成される環をどのように決定するか、指標がどのように重複度を数えるか、可換子環のランダムな元がなぜ分解を数値的に見つけるか、電荷セクターとは何か、既約なスピン表現と $`L`$ 個のスピン $`\frac12`$ の分解、$`SU(2)`$ 不変なハミルトニアンのブロックが $`U(1)`$ 不変なものよりはるかに小さい理由、複素共役が実行列によって生成される環のセクターをどのように実型・複素型・四元数型に分類するか、そして雑音を含む生成元の数値分解をどの許容誤差が正しくし、その出力がどのように自らを保証するか。本章の数値的な検算は
[チュートリアル 4](../../tutorial/ja/04_symmetries.md) にある。

本章を通じて、$`S'`$ は $`S\subseteq M_n(\mathbb C)`$ の可換子環である（第2章）。$`S^*=S`$ ならば、$`S''`$ は $`S`$ **が生成するフォン・ノイマン環** である：これは $`S`$ を含む $`*`$-環であり、$`S`$ を含む任意の $`*`$-環 $`\mathcal B\supseteq S`$ は
$`\mathcal B=\mathcal B''\supseteq S''`$ を満たす（定理 2.2）。$`\mathcal A'_h=\{K\in\mathcal A':K=K^*\}`$ は実ベクトル空間である。$`\mathbb C^k\otimes\mathbb C^m`$ には第7章のクロネッカー順序を用いる（$`e_a\otimes f_j`$ は $`0`$ から数えて $`am+j`$ 番目の基底ベクトル）。したがって $`1_k\otimes x`$ は $`x`$ の $`k`$ 個のコピーからなるブロック対角行列である。$`|G|`$ は有限群の位数であり、$`\hbar=1`$ とする。

<a id="s-4-1"></a>

## 4.1 表現とシューアの補題

**定義.** $`X`$ を対合 $`x\mapsto x^*`$（$`x^{**}=x`$）をもつ集合とする。$`\mathbb C^k`$ 上の $`X`$ の
**表現** とは、$`\rho(x^*)=\rho(x)^*`$ を満たす写像 $`\rho:X\to M_k(\mathbb C)`$ のことである。部分空間 $`K`$ は、すべての $`x`$ に対し $`\rho(x)K\subseteq K`$ のとき **不変** であるという。$`k\ge1`$ かつ $`0,\mathbb C^k`$ のみが不変部分空間であるとき、$`\rho`$ は **既約** であるという。$`\rho`$ から $`\mathbb C^{k'}`$ 上の表現 $`\rho'`$ への **絡作用素** とは、すべての $`x`$ に対し $`T\rho(x)=\rho'(x)T`$ を満たす $`k'\times k`$ 行列 $`T`$ のことである。それらはベクトル空間 $`\mathrm{Hom}(\rho,\rho')`$ をなし、$`\mathrm{Hom}(\rho,\rho)=\rho(X)'`$ である。可逆な絡作用素が存在するとき、$`\rho`$ と $`\rho'`$ は **同値** であるといい、$`\rho\cong\rho'`$ と書く。

三つの場合が現れる：$`g^*=g^{-1}`$ とした群 $`G`$ で、$`\rho`$ が準同型 $`G\to U(k)`$ であることを要求する場合（**ユニタリ表現**；このとき $`\rho(g^{-1})=\rho(g)^{-1}=\rho(g)^*`$）；随伴をもつ $`*`$-環で、$`\rho`$ が単位的 $`*`$-準同型である場合；そして $`x^*=x`$ を満たす生成元のラベル、たとえば
$`\{x,y,z\}`$ で、スピン行列のようなエルミート行列の組で表現される場合（4.5節）である。$`K`$ が不変で、$`Q\in M_{k\times d}`$ の列が $`K`$ の正規直交基底をなすならば $`\rho(x)Q=QQ^*\rho(x)Q`$ であり、したがって
$`\rho_K(x)=Q^*\rho(x)Q`$ はふたたび表現である（$`\rho`$ が準同型ならば準同型）。これを $`\rho`$ の $`K`$ への **制限** という。

<a id="r-4-1"></a>

**補題 4.1（完全可約性）.** $`\rho`$ を $`\mathbb C^n`$ 上の $`X`$ の表現とする。
(a) $`K`$ が不変ならば $`K^\perp`$ も不変である。(b) $`\mathbb C^n`$ は、各々の上で $`\rho`$ が既約であるような不変部分空間の直交直和である。

*証明.* (a) $`y\in K^\perp`$ と $`z\in K`$ に対し、$`\rho(x^*)z\in K`$ であるから $`\langle z,\rho(x)y\rangle=\langle\rho(x^*)z,y\rangle=0`$ である。(b) $`n`$ に関する帰納法。次元が最小の零でない不変部分空間 $`K_1`$ は既約である。(a) により $`K_1^\perp`$ は不変であり、$`K_1^\perp`$ への制限の不変部分空間は $`\rho`$ の不変部分空間であるから、帰納法の仮定がそれに適用できる。$`\square`$

有限群ではユニタリ性は何の代償も要しない：$`\rho(g)`$ が単に可逆であるだけでも、$`P=\sum_g\rho(g)^*\rho(g)`$ は
$`P\ge\rho(e)^*\rho(e)=1`$ を満たすので $`P>0`$ であり、$`\rho(h)^*P\rho(h)=\sum_g\rho(gh)^*\rho(gh)=P`$ である（$`g\to gh^{-1}`$ と置換する）。
したがって $`(P^{1/2}\rho(h)P^{-1/2})^*(P^{1/2}\rho(h)P^{-1/2})=P^{-1/2}\rho(h)^*P\rho(h)P^{-1/2}=1`$ であり、
$`h\mapsto P^{1/2}\rho(h)P^{-1/2}`$ は $`\rho`$ と同値なユニタリ表現である（ワイルのトリック）。

<a id="r-4-2"></a>

**補題 4.2（シューアの補題）.** $`\rho,\rho'`$ を $`\mathbb C^k,\mathbb C^{k'}`$ 上の $`X`$ の既約表現とし、$`T\in\mathrm{Hom}(\rho,\rho')`$ とする。
(a) $`T=0`$ であるか $`T`$ は可逆である。(b) $`\mathrm{Hom}(\rho,\rho)=\mathbb C1`$。
(c) $`T^*\in\mathrm{Hom}(\rho',\rho)`$ であり、$`T\neq0`$ ならば $`\lambda>0`$ とユニタリ $`U`$ により $`T=\lambda U`$ と書ける。
したがって $`\dim\mathrm{Hom}(\rho,\rho')`$ は $`\rho\cong\rho'`$ ならば $`1`$、そうでなければ $`0`$ であり、同値な既約表現はユニタリ同値である。

*証明.* (a) と (b) は命題 2.7 の後で群に対して証明したが、その証明が用いるのは、$`\ker T`$ が $`\rho`$ で不変であること（$`Tv=0`$ ならば $`T\rho(x)v=\rho'(x)Tv=0`$）と、$`T`$ の像が $`\rho'`$ で不変であることだけである：よって $`\ker T\in\{0,\mathbb C^k\}`$ かつ $`\mathrm{ran}\,T\in\{0,\mathbb C^{k'}\}`$ である。(b) については、$`\lambda`$ を $`T`$ の固有値とすると $`T-\lambda1`$ は可逆でない絡作用素であり、したがって $`0`$ である。(c) $`T\rho(x^*)=\rho'(x^*)T`$ の随伴は $`\rho(x)T^*=T^*\rho'(x)`$ である。すると $`T^*T\in\mathrm{Hom}(\rho,\rho)`$ であるから、$`T\neq0`$ ならば $`\lambda>0`$ により $`T^*T=\lambda^21`$ となり（$`T^*T\ge0`$ は零でない）、$`U=T/\lambda`$ は (a) により正方行列で $`U^*U=1`$ を満たす。最後に、$`T_1,T_2\in\mathrm{Hom}(\rho,\rho')`$ が零でなければ $`T_1^{-1}T_2\in\mathrm{Hom}(\rho,\rho)=\mathbb C1`$ である。$`\square`$

すべての $`\rho(x)`$ が互いに可換ならば（アーベル群）、既約な $`\rho`$ は一次元である：各 $`\rho(x)`$ は
$`\mathrm{Hom}(\rho,\rho)=\mathbb C1`$ に属するので、すべての部分空間が不変となる。

<a id="s-4-2"></a>

## 4.2 等型分解と可換子環

<a id="r-4-3"></a>

**定理 4.3（等型分解）.** $`\rho`$ を $`\mathbb C^n`$ 上の $`X`$ の表現とする。
$`\mathbb C^{k_c}`$ 上の $`X`$ の互いに非同値な既約表現 $`\rho_1,\dots,\rho_C`$、$`\sum_ck_cm_c=n`$ を満たす整数 $`m_c\ge1`$、およびユニタリ
$`V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n`$ が存在して、すべての $`x`$ に対し
$`V^*\rho(x)V=\bigoplus_c\rho_c(x)\otimes1_{m_c}`$ となる。そのような任意の分解について
```math
V^*\rho(X)'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),\qquad
V^*\rho(X)''V=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},
```
であり、$`m_c=\dim\mathrm{Hom}(\rho_c,\rho)`$ である。$`X`$ の既約表現 $`\sigma`$ が
$`\mathrm{Hom}(\sigma,\rho)\neq0`$ を満たすのは、ちょうど一つの $`c`$ に対して $`\sigma\cong\rho_c`$ となるときかつそのときに限る。

*証明.* *存在.* 補題 4.1 により、和が $`\mathbb C^n`$ で、制限 $`\sigma_i(x)=Q_i^*\rho(x)Q_i`$ が既約であるような、互いに直交する不変部分空間 $`W_1,\dots,W_M`$ が得られる。$`W_i`$ を制限の同値性によって類に分け、各類 $`c`$ において基準となる $`\rho_c=\sigma_{i_0}`$ を選ぶ。その類の他の各 $`i`$ について、補題 4.2(c) により $`\sigma_i(x)=U_i\rho_c(x)U_i^*`$ を満たすユニタリ $`U_i`$ が得られる。$`Q_i`$ を $`Q_iU_i`$ に置き換えれば $`Q_i^*\rho(x)Q_i=\rho_c(x)`$ となる。類 $`c`$ の元に $`j=0,\dots,m_c-1`$ と番号をつけ、それらの基底を $`Q_{c,j}`$ と書き、$`c`$ 番目の直和成分上で $`V(e_a\otimes f_j)=Q_{c,j}e_a`$ とおく。$`W_i`$ は不変で互いに直交するから $`Q_{c,j}^*\rho(x)Q_{d,l}=\delta_{cd}\delta_{jl}\rho_c(x)`$ であり、これが主張された形である。

*可換子環.* ここで $`V`$ を定理の主張にあるような任意の分解とし、$`Q_{c,j}`$ を（$`c`$ 番目の直和成分における）$`n\times k_c`$ 等長写像
$`e_a\mapsto V(e_a\otimes f_j)`$ とする。いま構成した $`V`$ に対しては、これらは上の $`Q_{c,j}`$ である。恒等式
$`V^*\rho(x)V=\bigoplus_c\rho_c(x)\otimes1_{m_c}`$ は $`\rho(x)Q_{d,l}=Q_{d,l}\rho_d(x)`$ を意味し、その随伴を $`x^*`$ に適用すると
$`Q_{c,j}^*\rho(x)=\rho_c(x)Q_{c,j}^*`$ を得る。$`V`$ はユニタリであるから、$`Y\in M_n`$ がすべての $`\rho(x)`$ と可換であることは、各ブロック
$`Y_{cj,dl}=Q_{c,j}^*YQ_{d,l}`$ が $`Y_{cj,dl}\rho_d(x)=\rho_c(x)Y_{cj,dl}`$ を満たすこと、すなわち
$`\mathrm{Hom}(\rho_d,\rho_c)`$ に属することと同値である。これは $`c\neq d`$ ならば $`0`$、$`c=d`$ ならば $`\mathbb C1`$ である（補題 4.2）。よって
$`Y\in\rho(X)'`$ であることと $`Y_{cj,dl}=\delta_{cd}(y_c)_{jl}1_{k_c}`$、すなわち $`V^*YV=\bigoplus_c1_{k_c}\otimes y_c`$ であることは同値である。

*二重可換子環.* $`Z`$ が $`\rho(X)'`$ と可換ならば、$`Z`$ は直和成分への射影と可換であり（$`y_c=1`$、$`d\neq c`$ に対し $`y_d=0`$ ととる）、したがって $`V^*ZV=\bigoplus_cZ_c`$ で、$`Z_c`$ は
$`1_{k_c}\otimes M_{m_c}`$ と可換である。補題 2.5 のステップ5（二つのテンソル因子を入れ替えたもの）により
$`Z_c\in M_{k_c}\otimes1_{m_c}`$ である。逆にこれらの作用素は $`\rho(X)'`$ と可換である。

*重複度.* 既約な $`\sigma`$ と $`T\in\mathrm{Hom}(\sigma,\rho)`$ に対し、同じ議論により
$`Q_{c,j}^*T\in\mathrm{Hom}(\sigma,\rho_c)`$ であり、これは $`\sigma\cong\rho_c`$ ならば可逆な $`t_c`$ による $`\mathbb Ct_{c}`$、そうでなければ
$`0`$ である。また $`T=\sum_{c,j}Q_{c,j}Q_{c,j}^*T`$ である。$`\rho_c`$ は互いに非同値であるから、$`\sigma`$ はそのうち高々一つ、たとえば $`\rho_c`$ と同値であり（そうでなければ $`T=0`$）、このとき $`T=\sum_ja_jQ_{c,j}t_c`$ である。逆に各 $`Q_{c,j}t_c`$ は $`\mathrm{Hom}(\sigma,\rho)`$ に属し、これら $`m_c`$ 個の写像は互いに直交する像をもつので、
$`\dim\mathrm{Hom}(\sigma,\rho)=m_c`$ である。$`\square`$

*ライブラリでは：* 定理 4.3 のような分解は `decompose.from_group` と `decompose.from_generators` で計算する — [チュートリアル、4.3節](../../tutorial/ja/04_symmetries.md#s-4-3)。

群の場合、$`\rho(G)''`$ は $`\rho(G)`$ の張る空間であり（2.4節）、定理 4.3 は命題 2.7 を逆向きに読んだものである：そこでは構造定理が既約表現を生み出し、ここでは既約表現が構造を生み出す。いずれにおいても $`\dim\rho(X)'=\sum_cm_c^2`$、$`\dim\rho(X)''=\sum_ck_c^2`$ である。$`\rho(X)''`$ の極小中心射影 $`P_c=\sum_jQ_{c,j}Q_{c,j}^*`$ の像が $`\rho_c`$ の **等型成分** である。
定理 4.3 により、任意のエルミートな $`H\in\rho(X)'`$ はエルミートな $`h_c\in M_{m_c}`$ により $`V(\bigoplus_c1_{k_c}\otimes h_c)V^*`$ と
書け、$`1_{k_c}\otimes h_c`$ は $`h_c`$ の $`k_c`$ 個のコピーからなるブロック対角行列である：$`H`$ のスペクトルは $`m_c\times m_c`$
ブロック $`h_c`$ のスペクトルの和集合であり、$`h_c`$ の各固有値は $`k_c`$ 回繰り返される。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. $`\rho(X)`$ が対称性で $`H\in\rho(X)'`$ がハミルトニアンならば、
> いま述べた $`k_c`$ 重に縮退した準位が対称性の **多重項** であり、対角化しなければならないブロックの大きさは、多重項の大きさではなく重複度 $`m_c`$ である。

有限群では、何も分解することなく重複度を計算できる。

<a id="r-4-4"></a>

**定理 4.4（シューアの直交関係）.** $`\rho,\rho'`$ を $`\mathbb C^k,\mathbb C^{k'}`$ 上の有限群 $`G`$ の既約ユニタリ表現とする。このとき
```math
\frac1{|G|}\sum_{g\in G}\rho(g)_{ij}\,\overline{\rho'(g)_{i'j'}}=
\begin{cases}\delta_{ii'}\delta_{jj'}/k&\text{if }\rho=\rho',\\0&\text{if }\rho\not\cong\rho'.\end{cases}
```

*証明.* $`k\times k'`$ 行列 $`Y`$ に対し $`T=\frac1{|G|}\sum_g\rho(g)Y\rho'(g)^*`$ とおく。$`h\in G`$ について、$`g\to hg`$ と置換し
$`\rho'(hg)^*=\rho'(g)^*\rho'(h)^{-1}`$ を用いると $`\rho(h)T=T\rho'(h)`$ を得るので、
$`T\in\mathrm{Hom}(\rho',\rho)`$ である。$`\rho\not\cong\rho'`$ ならば $`T=0`$ であり、$`\rho=\rho'`$ ならば
$`k\lambda=\mathrm{Tr}\,T=\mathrm{Tr}\,Y`$ を満たす $`\lambda`$ により $`T=\lambda1`$ である（補題 4.2）。$`Y=|j\rangle\langle j'|`$ に対し、成分 $`T_{ii'}`$ は左辺であり、$`\mathrm{Tr}\,Y=\delta_{jj'}`$ である。$`\square`$

表現の **指標** とは $`\chi_\rho(g)=\mathrm{Tr}\,\rho(g)`$ のことである。同値な表現は等しい指標をもつ。実際
$`\mathrm{Tr}(T\rho(g)T^{-1})=\mathrm{Tr}\,\rho(g)`$ である。
$`\langle\varphi,\psi\rangle_G=\frac1{|G|}\sum_g\overline{\varphi(g)}\psi(g)`$ と書く。

<a id="r-4-5"></a>

**系 4.5（指標）.** $`G`$ を有限群とする。
(a) 既約ユニタリな $`\rho,\rho'`$ に対し、$`\langle\chi_\rho,\chi_{\rho'}\rangle_G`$ は $`\rho\cong\rho'`$ ならば $`1`$、そうでなければ $`0`$ である。
(b) 定理 4.3 のように分解されたユニタリ表現 $`\pi`$ と既約な $`\rho`$ に対し、$`\pi`$ における $`\rho`$ の重複度は
```math
m_\rho=\frac1{|G|}\sum_{g\in G}\overline{\chi_\rho(g)}\,\chi_\pi(g)
```
である（$`\rho\cong\rho_c`$ ならば $`m_\rho=m_c`$、$`\rho`$ がどの $`\rho_c`$ とも同値でなければ $`m_\rho=0`$）。
(c) $`\dim\pi(G)'=\sum_cm_c^2=\frac1{|G|}\sum_{g\in G}|\chi_\pi(g)|^2`$。特に、$`\pi`$ が既約であるのはこの数が $`1`$ のときかつそのときに限る。

*証明.* (a) $`\rho'`$ を同値な表現（同じ指標をもつ）に置き換える。定理 4.4 を
$`i=j`$、$`i'=j'`$ について和をとり、複素共役をとる：結果は $`\rho=\rho'`$ ならば $`\sum_{i,i'}\delta_{ii'}/k=1`$、
$`\rho\not\cong\rho'`$ ならば $`0`$ である。(b) トレースはユニタリ不変であり $`\mathrm{Tr}(a\otimes1_m)=m\mathrm{Tr}\,a`$ であるから、
$`\chi_\pi=\sum_cm_c\chi_{\rho_c}`$ である。(a) を適用すればよい。(c) 同じ展開と (a) により
$`\langle\chi_\pi,\chi_\pi\rangle_G=\sum_cm_c^2`$ であり、これは定理 4.3 により $`\dim\pi(G)'`$ である。$`\square`$

*計算：* 三つの量子ビットを置換する $`S_3`$ に対する系 4.5 — [チュートリアル、4.2節](../../tutorial/ja/04_symmetries.md#s-4-2)；$`S_3`$ の標準表現に対する定理 4.4 — [チュートリアル、4.3節](../../tutorial/ja/04_symmetries.md#s-4-3)。

(c) の第二の証明は、ライブラリが可換子環からサンプリングする方法そのものである。**ツワール**
$`P(Y)=\frac1{|G|}\sum_g\pi(g)Y\pi(g)^*`$ は像が $`\pi(G)'`$ の冪等写像であり（第2章、演習問題 3）、
$`\langle Y,Z\rangle_{HS}=\mathrm{Tr}(Y^*Z)`$ に関して自己共役である（$`g\to g^{-1}`$ と置換する）。したがって $`\pi(G)'`$ の上への直交射影であり、$`\dim\pi(G)'`$ はそのトレースである。写像 $`Y\mapsto AYB`$ のトレースは
$`\mathrm{Tr}\,A\mathrm{Tr}\,B`$ であり（その対角成分は $`A_{ii}B_{jj}`$）、これから $`\frac1{|G|}\sum_g|\chi_\pi(g)|^2`$ が得られる。

*ライブラリでは：* ツワールは `decompose.twirl` — [チュートリアル、4.3節](../../tutorial/ja/04_symmetries.md#s-4-3)；$`\dim\pi(G)'`$ は `decompose.commutant_dimension` — [チュートリアル、4.2節](../../tutorial/ja/04_symmetries.md#s-4-2)。

2.4節の例を指標を用いてやり直そう。置換行列のトレースはそれが固定する基底ベクトルの個数であり、$`\pi(\sigma)`$ が $`|b_1b_2b_3\rangle`$ を固定するのは $`b`$ が $`\sigma`$ の各巡回上で一定であるときかつそのときに限るので（補題 B.1(c)）、$`\chi_\pi(\sigma)=2^{\#\text{cycles}(\sigma)}`$ である。$`S_3`$ の三つのユニタリ表現が手元にある。$`P_\sigma e_i=e_{\sigma(i)}`$ を $`\mathbb C^3`$ 上の置換表現（置換行列への準同型）とする。**自明** 表現は $`\sigma\mapsto1`$ である。**符号** は
$`\mathrm{sgn}\,\sigma=\det P_\sigma`$ であり、行列式は乗法的で置換行列の行列式は $`\pm1`$ であるから、これは $`\{\pm1\}\subseteq U(1)`$ への準同型である。直線 $`\mathbb C(1,1,1)`$ は不変であり、したがって $`E=(1,1,1)^\perp`$ も不変である（補題 4.1(a)）。**標準** 表現は $`P`$ の $`E`$ への制限である。
$`\mathbb C(1,1,1)\oplus E`$ に適合した正規直交基底では $`P_\sigma=1\oplus\mathrm{std}(\sigma)`$ であるから、
$`\chi_{\rm std}(\sigma)=\mathrm{Tr}\,P_\sigma-1`$、すなわち固定点の個数から $`1`$ を引いたものである。$`S_3`$ の六つの元は、
単位元、三つの互換、二つの $`3`$-巡回であり、巡回の個数はそれぞれ $`3,2,1`$、固定点の個数はそれぞれ $`3,1,0`$ である。互換は
$`1_3`$ の二つの列を入れ替えるのでその符号は $`-1`$ であり、$`3`$-巡回は二つの互換の積（$`(1\,2\,3)=(1\,3)(1\,2)`$）であるからその符号は $`+1`$ である。こうして次の表を得る。

| | 単位元（$`1`$） | 互換（$`3`$） | $`3`$-巡回（$`2`$） |
|---|---|---|---|
| $`\chi_{\rm triv}`$ | $`1`$ | $`1`$ | $`1`$ |
| $`\mathrm{sgn}`$ | $`1`$ | $`-1`$ | $`1`$ |
| $`\chi_{\rm std}`$ | $`2`$ | $`0`$ | $`-1`$ |
| $`\chi_\pi`$ | $`8`$ | $`4`$ | $`2`$ |

（したがって $`\mathrm{sgn}\,\sigma=(-1)^{3-\#\text{cycles}(\sigma)}`$ であり、すべての指標は実数値である。）$`\langle\varphi,\psi\rangle_G`$ を共役類ごとに和をとって計算すると、
$`\langle\chi_{\rm triv},\chi_{\rm triv}\rangle_G=\frac{1+3+2}6=1`$、$`\langle\mathrm{sgn},\mathrm{sgn}\,\rangle_G=\frac{1+3+2}6=1`$、
$`\langle\chi_{\rm std},\chi_{\rm std}\rangle_G=\frac{4+0+2}6=1`$ であるから、系 4.5(c) によりこれら三つの表現は既約であり、指標が異なるので互いに非同値である。系 4.5(b) により $`\pi`$ における重複度は
```math
m_{\rm triv}=\tfrac{8+3\cdot4+2\cdot2}6=4,\qquad m_{\rm sgn}=\tfrac{8-3\cdot4+2\cdot2}6=0,\qquad
m_{\rm std}=\tfrac{2\cdot8+0-2\cdot2}6=2 .
```
である。したがって自明表現と標準表現は、$`m_c=4`$ および $`m_c=2`$ をもつ $`\pi`$ の構成成分 $`\rho_c`$（定理 4.3）と同値である。これらで
$`4\cdot1+2\cdot2=8=\sum_ck_cm_c`$ 次元が尽くされるので、他の構成成分はない：$`k=(1,2)`$、$`m=(4,2)`$ であり、
$`\dim\pi(S_3)'=4^2+2^2=\frac16(64+3\cdot16+2\cdot4)=20`$ である（系 4.5(c)）。

これは2.4節の分解であり、今度は表現ごとに六つの数からなる表から得られた。

*計算：* この表を $`\pi`$ から計算し直す — [チュートリアル、4.2節](../../tutorial/ja/04_symmetries.md#s-4-2)。

<a id="s-4-3"></a>

## 4.3 生成元から環へ：ランダム化アルゴリズム

定理 4.3 は分解を記述するが、それを計算はしない：その証明には既約な不変部分空間が必要である。Murota、Kanno、Kojima、Kojima（*Japan J. Indust. Appl. Math.* 27 (2010) 125–160。Maehara and Murota, *ibid.* 263–293 も参照）は、*可換子環のランダムな元* の固有空間がそれらを見つけることに気づいた。$`\mathcal A`$ を与えられた行列 $`S_1,\dots,S_J`$ とその随伴が生成するフォン・ノイマン環とし、その分解を $`V^*\mathcal AV=\bigoplus_cM_{k_c}\otimes1_{m_c}`$、
$`V^*\mathcal A'V=\bigoplus_c1_{k_c}\otimes M_{m_c}`$ とする（定理 2.6）。$`V`$ を用いるのは解析だけである。

<a id="r-4-6"></a>

**補題 4.6.** $`K\in\mathcal A'_h`$ とし、$`V^*KV=\bigoplus_c1_{k_c}\otimes\kappa_c`$ とする。$`K`$ の各固有空間は
$`\mathcal A`$ 不変である。すべての $`\kappa_c`$ を合わせた $`N=\sum_cm_c`$ 個の固有値が互いに相異なるならば、$`K`$ の固有空間はちょうど $`N`$ 個の部分空間 $`V(\mathbb C^{k_c}\otimes e)`$（$`c`$ 番目の直和成分において、$`e`$ は $`\kappa_c`$ の固有ベクトル）であり、$`\mathcal A`$ はそれぞれの上に $`M_{k_c}`$ として既約に作用する。

*証明.* $`Kv=\lambda v`$ かつ $`a\in\mathcal A`$ ならば $`Kav=aKv=\lambda av`$ である。正規直交な $`e_{cj}`$ により
$`\kappa_c=\sum_j\mu_{cj}|e_{cj}\rangle\langle e_{cj}|`$ と書く（定理 1.3）。すると
$`V^*KV=\sum_{c,j}\mu_{cj}R_{cj}`$ であり、ここで $`R_{cj}`$ は $`\mathbb C^{k_c}\otimes e_{cj}`$ の上への射影である。
$`R_{cj}`$ は互いに直交して和が $`1`$ であり、$`\mu_{cj}`$ は相異なるので、これらがスペクトル射影である。
$`\mathbb C^{k_c}\otimes e_{cj}`$ 上で元 $`V^*aV=\bigoplus_da_d\otimes1`$ は $`a_c`$ として作用し、$`a_c`$ は
$`M_{k_c}`$ 全体を動く。$`M_{k_c}`$ は $`0`$ と $`\mathbb C^{k_c}`$ 以外に不変部分空間をもたない。$`\square`$

<a id="r-4-7"></a>

**命題 4.7（一般性）.** $`K`$ と $`Y`$ を、それぞれ $`\mathcal A'_h`$ と
$`\mathcal A'`$ の独立なランダムな元で、その分布がこれらの実ベクトル空間上のルベーグ測度に関して密度をもつもの（補題 A.3）とする。このとき確率 $`1`$ で：(a) 補題 4.6 の仮定が成り立つ。(b) 像が同じ直和成分に含まれる $`K`$ の任意の二つの相異なるスペクトル射影 $`\Pi_i,\Pi_j`$ に対し $`\Pi_iY\Pi_j\neq0`$ である。

*証明.* (a) エルミート行列 $`\kappa=\bigoplus_c\kappa_c\in M_N`$ は $`K`$ に実線形に依存する。$`\langle\cdot,\cdot\rangle_{HS}`$ に関する $`1,\kappa,\dots,\kappa^{N-1}`$ のグラム行列式
$`\delta(K)=\det[\mathrm{Tr}\,\kappa^{i+j}]_{i,j=0}^{N-1}`$ は $`K`$ の座標の実多項式であり、これらの冪が線形独立であるときかつそのときに限り零でない。
$`\kappa=U\mathrm{diag}(\lambda)U^*`$ に対して $`p(\kappa)=U\mathrm{diag}(p(\lambda_l))U^*`$ であるから、次数 $`<N`$ の零でない多項式が $`\kappa`$ を零化するのは、それが $`\mathrm{spec}\,\kappa`$ 上で消えるときかつそのときに限り、それが可能なのは $`\kappa`$ の相異なる固有値が $`N`$ 個未満のときかつそのときに限る（それらの上で $`\prod_\mu(t-\mu)`$ をとる）。よって $`\delta(K)\neq0`$ は補題 4.6 の仮定であり、
$`\delta\not\equiv0`$ である：実際 $`V(\bigoplus_c1\otimes\mathrm{diag}(\text{相異なる実数}))V^*\in\mathcal A'_h`$ において零でない。したがって
系 A.7(b) により $`\mathbb P(\delta(K)=0)=0`$ である。

(b) 集合 $`O=\{\delta\neq0\}\subseteq\mathcal A'_h`$ は開集合である（$`\delta`$ は連続）。$`K\in O`$ に対し、$`K`$ はちょうど $`N`$ 個の相異なる固有値
$`\mu_1(K)<\dots<\mu_N(K)`$ をもつ（$`\kappa`$ の固有値で、重複度は $`k_c`$）ので、固有値の昇順に番号をつけたそのスペクトル射影
$`\Pi_1(K),\dots,\Pi_N(K)`$ は $`O`$ 上で連続である（命題 A.13）。$`l\neq l'`$ を固定し、$`\mathcal N_{ll'}`$ を、ある $`c`$ に対して
$`P_c\Pi_l(K)=\Pi_l(K)`$ かつ $`P_c\Pi_{l'}(K)=\Pi_{l'}(K)`$（$`P_c`$ は中心射影）であり、かつ $`\Pi_l(K)Y\Pi_{l'}(K)=0`$ であるような
$`(K,Y)\in O\times\mathcal A'`$ 全体の集合とする。これはボレル集合 $`O\times\mathcal A'`$ 上の連続な方程式で定義されるので、ボレル集合である。二つの像が直和成分 $`c`$ に含まれるような $`K\in O`$ を固定する。補題 4.6 により、
$`V^*\Pi_lV`$ と $`V^*\Pi_{l'}V`$ は、$`\kappa_c`$ の正規直交な固有ベクトル $`e,e'`$ により、それぞれ $`\mathbb C^{k_c}\otimes e`$ と $`\mathbb C^{k_c}\otimes e'`$ の上への射影であり、
$`V^*YV=\bigoplus_d1\otimes y_d`$ に対して $`\Pi_lY\Pi_{l'}=0`$ であることは $`\langle e,y_ce'\rangle=0`$ と同値である。写像
$`Y\mapsto\langle e,y_ce'\rangle`$ は実線形写像 $`\mathcal A'\to\mathbb C`$ であって零でない（
$`V(0\oplus\dots\oplus1\otimes|e\rangle\langle e'|\oplus\dots\oplus0)V^*`$ における値は $`1`$）ので、その核は真部分空間であり、零集合である（系 A.7(c)）。
よってすべての切り口 $`\{Y:(K,Y)\in\mathcal N_{ll'}\}`$ は零集合であり（二つの像が異なる直和成分に含まれるならば空集合）、$`Y`$ は $`K`$ と独立であるから、
系 A.7(d) により $`\mathbb P((K,Y)\in\mathcal N_{ll'})=0`$ である。(b) の事象が成り立たないのは
$`\{\delta(K)=0\}\cup\bigcup_{l\neq l'}\{(K,Y)\in\mathcal N_{ll'}\}`$ の上だけであり、これは確率 $`0`$ の事象の有限個の和集合である。$`\square`$

<a id="r-4-8"></a>

**補題 4.8（圧縮；コピーの整列）.** (a) $`W`$ を正規直交基底 $`Q\in M_{n\times d}`$ をもつ $`\mathcal A`$ 不変部分空間とし、
$`\sigma(a)=Q^*aQ`$ をその制限とする。このとき $`\sigma(\mathcal A)`$ は $`*`$-環であり、
$`\{Q^*YQ:Y\in\mathcal A'\}=\sigma(\mathcal A)'`$ かつ $`\{Q^*KQ:K\in\mathcal A'_h\}=\sigma(\mathcal A)'_h`$ である。$`\mathcal A`$ が $`W`$ 上に既約に作用するのは、
すべての $`Y\in\mathcal A'`$ に対して $`Q^*YQ\in\mathbb C1`$ となるときかつそのときに限る。そうでないとき、$`X`$ を $`\mathcal A'`$ または
$`\mathcal A'_h`$ の密度をもつランダムな元とすると、確率 $`1`$ で $`Q^*XQ\notin\mathbb C1`$ である。

(b) $`W_i,W_j`$ を、$`\mathcal A`$ がその上に既約に作用する不変部分空間とし、その正規直交基底を $`Q_i,Q_j`$、制限を $`\sigma_i(a)=Q_i^*aQ_i`$ とする。また $`Y\in\mathcal A'`$ とする。
このとき $`T=Q_i^*YQ_j\in\mathrm{Hom}(\sigma_j,\sigma_i)`$ である。$`T\neq0`$ ならば $`\lambda>0`$ とユニタリ $`U`$ により $`T=\lambda U`$ と書け、$`Q_j'=Q_jU^*`$ はすべての $`a\in\mathcal A`$ に対し $`Q_j'^*aQ_j'=\sigma_i(a)`$ を満たす。$`W_i,W_j`$ が異なる直和成分に含まれるならば（$`W_i\subseteq P_c\mathbb C^n`$、$`W_j\subseteq P_d\mathbb C^n`$、$`c\neq d`$）、$`T=0`$ である。

*証明.* (a) 定理 4.3 の証明と同様に（$`W`$ は $`a`$ と $`a^*`$ で不変である）、$`a\in\mathcal A`$ に対し $`aQ=Q\sigma(a)`$ かつ
$`Q^*a=\sigma(a)Q^*`$ である。したがって $`\sigma(ab)=Q^*aQ\sigma(b)=\sigma(a)\sigma(b)`$、$`\sigma(a^*)=\sigma(a)^*`$、
$`\sigma(1)=Q^*Q=1_d`$ であるから、$`\sigma`$ は単位的 $`*`$-準同型であり、$`\sigma(\mathcal A)`$ は $`*`$-環である。また $`Y\in\mathcal A'`$ に対し $`Q^*YQ\sigma(a)=Q^*YaQ=Q^*aYQ=\sigma(a)Q^*YQ`$ である。
逆に、$`T\in\sigma(\mathcal A)'`$ に対し $`Y=QTQ^*`$ とおくと、$`aY=Q\sigma(a)TQ^*=QT\sigma(a)Q^*=Ya`$ であるから $`Y\in\mathcal A'`$ であり、
$`Q^*YQ=T`$ で、$`T`$ がエルミートならば $`Y`$ もエルミートである。$`\sigma`$ が既約ならば $`\sigma(\mathcal A)'=\mathrm{Hom}(\sigma,\sigma)=\mathbb C1`$ である
（補題 4.2(b)）。$`E\subseteq\mathbb C^d`$ が $`0,\mathbb C^d`$ 以外の $`\sigma`$ 不変部分空間ならば、$`E^\perp`$ もそうであり（補題
4.1(a)）、$`E`$ の上への直交射影は $`\sigma(\mathcal A)'_h`$ のスカラーでない元である。この場合
$`\mathbb C1\cap\sigma(\mathcal A)'`$ と $`\mathbb R1\cap\sigma(\mathcal A)'_h`$ は真の実部分空間であり、全射な実線形写像 $`X\mapsto Q^*XQ`$ によるそれらの逆像は
$`\mathcal A'`$ と $`\mathcal A'_h`$ の真部分空間であって、これらは零集合である（系 A.7(c)）。

(b) (a) と同様に $`aQ_j=Q_j\sigma_j(a)`$ かつ $`Q_i^*a=\sigma_i(a)Q_i^*`$ であるから、
$`T\sigma_j(a)=Q_i^*YaQ_j=Q_i^*aYQ_j=\sigma_i(a)T`$ である。補題 4.2 により $`T=\lambda U`$ かつ
$`U\sigma_j(a)U^*=\sigma_i(a)`$ であり、後者は $`Q_j'^*aQ_j'`$ である。最後に $`P_c\in\mathcal A\cap\mathcal A'`$、$`Q_i=P_cQ_i`$、
$`P_cQ_j=0`$ であるから、$`T=Q_i^*P_cYQ_j=Q_i^*YP_cQ_j=0`$ である。$`\square`$

*計算：* $`\pi(S_3)'`$ のランダムな元に対する補題 4.6 と補題 4.8(b) — [チュートリアル、4.3節](../../tutorial/ja/04_symmetries.md#s-4-3)。

以上により、アルゴリズム（`decompose.wedderburn`、フロントエンドは `from_group` と `from_generators`）は次のとおりである：

1. *可換子環からのサンプリング.* 有限群の全元に対しては、ツワール $`P`$ と独立なジニブル行列 $`Z,Z'`$ により $`Y=P(Z)`$、$`K=P(Z'+Z'^*)`$ とする。
   $`P`$ は $`M_n`$ から $`\mathcal A'`$ の上への線形射影で $`P(Z)^*=P(Z^*)`$ を満たすので、$`Z\mapsto P(Z)`$ と $`Z'\mapsto P(Z'+Z'^*)`$ は $`\mathcal A'`$ の上および
   $`\mathcal A'_h`$ の上への実線形全射である（$`H\in\mathcal A'_h`$ に対し $`P(\frac H2+\frac H2)=H`$）。ジニブル行列は $`M_n`$ の非退化な
   中心化ガウスベクトルである（命題 A.16）から、$`Y`$ と $`K`$ は $`\mathcal A'`$ と $`\mathcal A'_h`$ の非退化な中心化ガウスベクトルであり、
   密度をもつ（定理 A.15）。生成元に対しては、
   $`\mathcal A'=\{X:[S_j,X]=[S_j^*,X]=0\}`$ は $`M_n`$ 上の線形写像の核である（ライブラリはエルミート部分
   $`\frac12(S_j+S_j^*)`$、$`\frac1{2i}(S_j-S_j^*)`$ を用いる。これらは $`S_j,S_j^*`$ と同じ複素空間を張るので、同じ可換子環をもつ）。その基底 $`(B_r)_{r\le D}`$ をとり、
   $`Y=\sum_rc_rB_r`$、$`K=Y'+Y'^*`$ とする。ここで $`Y'`$ は $`Y`$ の独立なコピー、$`c_r`$ は独立な標準複素ガウス変数である：これらは
   $`\mathbb C^D`$ の非退化なガウスベクトルの、$`\mathcal A'`$ の上への実線形全単射 $`c\mapsto\sum_rc_rB_r`$ による像、および $`\mathcal A'_h`$ の上への実線形全射
   $`c\mapsto X+X^*`$（$`X=\sum_rc_rB_r`$）による像であり、したがってやはり密度をもつ。どのランダムな元も新たに抽出されたもので、それ以前のものとは独立である。
   ライブラリはそれらを単位ノルムに規格化するが、これはどの固有空間も、どのブロックが消えるかも変えない。
2. *分割.* $`K`$ の固有空間は既約である（補題 4.6、命題 4.7(a)）。それらの次元が $`k_c`$ である。
3. *グループ化と整列.* 等しい次元の固有空間の間のブロック $`Q_i^*YQ_j`$ は、ちょうど同じ直和成分の内部で零でない（命題 4.7(b)、
   補題 4.8(b)）。「零でないブロックで結ばれる」という関係の類がセクターを与え、その大きさが $`m_c`$ を与え、極分解の因子が各コピーを類ごとに一つの基準コピーに整列させる
   （補題 4.8(b)。整列されたコピーに整列されたコピーは基準に整列されている）。
4. クロネッカー順序で $`V`$ を *組み立て*、$`V^*aV=\bigoplus_ca_c\otimes1_{m_c}`$ となるようにする。

したがって、厳密な算術で $`K`$ の固有空間が厳密に分離されるならば、アルゴリズムは確率 $`1`$ で定理 4.3 のような分解を返す：補題 4.6 により直和成分 $`c`$ は次元 $`k_c`$ の固有空間を $`m_c`$ 個含み、
整列の後の構成は定理 4.3 の証明の存在部分そのものである。浮動小数点演算では、ライブラリは許容誤差より近い固有値を等しいものとして扱うので、計算された「固有空間」$`W`$ は $`K`$ のいくつかの
固有空間の和であるかもしれない。それでも $`W`$ は $`\mathcal A`$ 不変であり（補題 4.6）、ライブラリは新しいランダムな
$`X\in\mathcal A'_h`$ をそれに圧縮する。補題 4.8(a) により、$`\mathcal A`$ が $`W`$ 上に既約に作用するならば $`Q^*XQ`$ は $`1`$ の倍数であり、
そうでなければ確率 $`1`$ で $`1`$ の倍数ではない。その場合、その固有空間は $`\sigma(\mathcal A)`$ で不変であり（$`\sigma(\mathcal A)`$ に適用した補題 4.6）、
$`W`$ をより小さな不変部分空間に分割する。ライブラリはすべての圧縮がスカラーになるまでこれを繰り返す。浮動小数点演算で何が起ころうと、結果は最後に証明書によって確認される：

<a id="r-4-9"></a>

**命題 4.9（証明書）.** $`\mathcal A`$ を $`S_1,\dots,S_J`$ とその随伴が生成する環とし、
ユニタリ $`V`$ により $`\mathcal B=V(\bigoplus_cM_{k_c}\otimes1_{m_c})V^*`$ とする。すべての $`S_j\in\mathcal B`$ かつ
$`\dim\mathcal A'=\sum_cm_c^2`$ ならば、$`\mathcal A=\mathcal B`$ である。

*証明.* $`\mathcal B`$ はすべての $`S_j,S_j^*`$ を含む $`*`$-環であるから $`\mathcal A\subseteq\mathcal B`$ かつ
$`\mathcal B'\subseteq\mathcal A'`$ である。$`\mathbb C^n`$ 上の $`\mathcal B`$ の表現に定理 4.3 を適用すると（ブロック写像 $`b\mapsto b_c`$ は既約であり、かつ非同値である。$`c`$ 番目の中心射影はそのうちの一つによって $`1`$ に、他によって $`0`$ に写されるからである）、$`\dim\mathcal B'=\sum_cm_c^2=\dim\mathcal A'`$ を得る。よって $`\mathcal B'=\mathcal A'`$ であり、
$`\mathcal B=\mathcal B''=\mathcal A''=\mathcal A`$ である。$`\square`$

*ライブラリでは：* アルゴリズム `decompose.wedderburn`（フロントエンド `from_group`、`from_generators`）、命題 4.9 の証明書 `Decomposition.verify()`、およびそれらの計算コスト — [チュートリアル、4.3節](../../tutorial/ja/04_symmetries.md#s-4-3)。

4.9節では浮動小数点の部分を定量化する：雑音を含む生成元に対してどの許容誤差が正しい可換子環を与えるか（定理 4.25）、ランダムな元の
固有値ギャップがどれほど小さいか（命題 4.27）、計算された分解が厳密な分解にどれほど近いか（定理 4.28）、そして出力がそれ自身で何を
保証するか（命題 4.26）。

<a id="s-4-4"></a>

## 4.4 $`U(1)`$ と電荷セクター

<a id="r-4-10"></a>

**命題 4.10.** $`Q`$ を整数の固有値 $`q`$（**電荷**）をもつ $`\mathbb C^n`$ 上のエルミート行列とし、そのスペクトル射影を $`P_q`$、$`d_q=\mathrm{rank}\,P_q`$、$`U_\theta=e^{i\theta Q}`$ とする。
$`U(1)'=\{U_\theta:\theta\in\mathbb R\}'`$ と書き、$`U(1)''`$ でその可換子環を表す。

1. $`\theta\mapsto U_\theta`$ は $`U(1)=\mathbb R/2\pi\mathbb Z`$ のユニタリ表現であり、$`\mathbb C^n`$ 上の $`U(1)`$ の任意の連続な
   ユニタリ表現はこの形をもつ。
2. $`X\in M_n`$ に対し次は同値である：(i) すべての $`\theta`$ に対し $`XU_\theta=U_\theta X`$。(ii) $`[X,Q]=0`$。
   (iii) $`X=\sum_qP_qXP_q`$。したがって $`U(1)'\cong\bigoplus_qM_{d_q}(\mathbb C)`$（ブロック対角作用素）であり、
   $`U(1)''=\mathrm{span}\{P_q\}`$ である：すべての $`k_c=1`$ で、重複度は $`d_q`$ である。
3. （*保存電荷.*）エルミートな $`H`$ に対し、すべての $`t`$ について $`e^{itH}Qe^{-itH}=Q`$ であることは $`[H,Q]=0`$、すなわち
   $`H\in U(1)'`$ と同値である。このとき $`\rho(t)=e^{-itH}\rho e^{itH}`$ のすべてのセクターの重み $`\mathrm{Tr}(P_q\rho(t))`$ は一定である。

*証明.* (1) $`U_\theta=\sum_qe^{iq\theta}P_q`$ はユニタリであり、$`U_\theta U_{\theta'}=U_{\theta+\theta'}`$（定理 1.5(c)）、
$`U_{2\pi}=1`$ である。逆に、連続な表現に対して $`(U_\theta)_{\theta\ge0}`$ は行列の連続半群であるから、補題 8.5（第8章で証明する。その証明は有限次元であることしか用いない）により $`\theta\ge0`$ に対し $`U_\theta=e^{\theta B}`$ であり、$`\theta<0`$ に対しては $`U_\theta=U_{-\theta}^{-1}=e^{\theta B}`$ である。
$`U_\theta^*U_\theta=1`$ を $`0`$ で微分すると（補題 8.1）$`B^*+B=0`$ となるので、エルミートな $`Q`$ により $`B=iQ`$ と書ける。指数級数
$`e^{i\theta Q}`$ は関数計算に等しい（$`H=-Q`$ とした命題 8.2(1)）から、$`U_{2\pi}=\sum_qe^{2\pi iq}P_q=1`$ より $`Q`$ のすべての固有値 $`q`$ について
$`e^{2\pi iq}=1`$、$`q\in\mathbb Z`$ が従う。(2) (i)$`\Rightarrow`$(ii)：$`\theta=0`$ で微分する。
(ii)$`\Rightarrow`$(iii)：$`X`$ は $`P_q=f(Q)`$ と可換である（定理 1.5(d)）から、$`q\neq q'`$ に対し $`P_qXP_{q'}=XP_qP_{q'}=0`$ である。
(iii)$`\Rightarrow`$(i)：$`XU_\theta=\sum_qe^{iq\theta}P_qXP_q=U_\theta X`$。ブロック対角作用素は（$`P_q`$ の像に適合した正規直交基底において）
$`\bigoplus_qM_{d_q}\otimes1_1`$ をなし、その可換子環は $`\bigoplus_q1_{d_q}\otimes M_1=\mathrm{span}\{P_q\}`$ である：それらと可換な作用素はすべての $`P_q`$ と可換であるから
ブロック対角であり、各ブロックは $`M_{d_q}`$ と可換であるからスカラーである（$`m=1`$ とした補題 2.5 のステップ5）。
(3) $`\frac d{dt}e^{itH}Qe^{-itH}=ie^{itH}[H,Q]e^{-itH}`$（命題 8.2）がすべての $`t`$ で消えることは $`[H,Q]=0`$ と同値である。
このとき $`e^{itH}`$ は $`P_q`$ と可換であり（定理 1.5(d) を二回）、
$`\mathrm{Tr}(P_q\rho(t))=\mathrm{Tr}(e^{itH}P_qe^{-itH}\rho)=\mathrm{Tr}(P_q\rho)`$ である。$`\square`$

したがって $`U(1)`$ の既約表現は指標 $`\theta\mapsto e^{iq\theta}`$、$`q\in\mathbb Z`$ であり（任意のアーベル群と同様に一次元）、電荷 $`q`$ の重複度は $`d_q`$ である。電荷はテンソル積のもとで足し合わされる：$`u_k`$ と $`w_l`$ をそれぞれ固有値 $`\alpha_k`$ と $`\beta_l`$ をもつ $`Q_1`$ と $`Q_2`$ の正規直交固有基底とすると、
正規直交基底 $`u_k\otimes w_l`$ は $`Q_1\otimes1+1\otimes Q_2`$ の固有ベクトル（固有値 $`\alpha_k+\beta_l`$）からなり、同時に
$`e^{i\theta Q_1}\otimes e^{i\theta Q_2}`$ の固有ベクトル（固有値 $`e^{i\theta\alpha_k}e^{i\theta\beta_l}`$）からなる。したがって関数計算の定義により
```math
e^{i\theta(Q_1\otimes1+1\otimes Q_2)}=e^{i\theta Q_1}\otimes e^{i\theta Q_2},
```
であり、これは $`P_{q_1}\otimes P_{q_2}`$ の像の上に電荷 $`q_1+q_2`$ で作用する。これが第7章の融合テンソル積の融合則 $`(q_1,q_2)\mapsto q_1+q_2`$ である。

$`n`$ 個の量子ビットに対しては、電荷は **個数作用素** $`\hat N=\sum_jn_j`$ である。ここで $`n_j=|1\rangle\langle1|`$ は量子ビット
$`j`$ 上に作用する：ハミング重み $`|b|`$ により $`\hat N|b\rangle=|b|\,|b\rangle`$ である。そのセクターの次元は
$`\binom nw`$ であるから、ハミング重みを保存する作用素全体は
```math
U(1)'\cong\bigoplus_{w=0}^nM_{\binom nw}(\mathbb C),
```
をなす。これは電荷 $`0,1`$ をもつ $`M_1\oplus M_1`$ の $`n`$ 個のコピーの融合積である
（ビット $`1`$ = 上向きスピン、$`\hat N=S^z+\frac n2`$）。任意の XXZ ハミルトニアン
$`\sum_{\langle ij\rangle}J(s^x_is^x_j+s^y_is^y_j+\Delta s^z_is^z_j)`$ は $`\hat N`$ を保存する：$`s^z_is^z_j`$ は基底 $`|b\rangle`$ において対角であり、ホッピング項
$`s^x_is^x_j+s^y_is^y_j=\frac12(s^+_is^-_j+s^-_is^+_j)`$（$`s^\pm=s^x\pm is^y`$ を展開すると交差項は打ち消し合う）は $`|b\rangle`$ を $`0`$ に、あるいは
上向きスピンがサイト $`i,j`$ の一方から他方へ移った配置に写すので、$`\hat N`$ の各固有空間をそれ自身に写すからである。

*計算：* 四つのスピンに対する命題 4.10：$`\hat N`$ だけから求めた電荷セクターと、セクターごとに対角化した XXZ 環 — [チュートリアル、4.4節](../../tutorial/ja/04_symmetries.md#s-4-4)。*ライブラリでは：* $`n`$ 個の量子ビットに対する $`U(1)'`$ は `SpinChain(n).algebra`、XXZ ハミルトニアンは `SpinChain.xxz` — [チュートリアル、4.4節](../../tutorial/ja/04_symmetries.md#s-4-4)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. $`U(1)`$ 対称性とは、粒子数、電荷、あるいは磁化の保存である。命題 4.10(3) はこの設定におけるネーターの定理である：電荷が保存されるのはハミルトニアンが対称なときかつそのときに限り、そのときダイナミクスはセクターを決して離れない（定理 8.3）。さらにセクター間のコヒーレンスを測定できるものが何もなければ、セクターは超選択セクターである（第3章）。

<a id="s-4-5"></a>

## 4.5 $`SU(2)`$：スピン

**定義.** $`\mathbb C^d`$ 上の **スピン表現** とは、
```math
[S_x,S_y]=iS_z,\qquad[S_y,S_z]=iS_x,\qquad[S_z,S_x]=iS_y,
```
を満たすエルミート行列の三つ組 $`S_x,S_y,S_z`$ のこと、すなわちこれらの関係をもつラベル $`\{x,y,z\}`$ の表現（4.1節）のことである。基本的な例は
$`\mathbb C^2`$ 上の $`s_a=\frac12\sigma_a`$（パウリ行列）である。$`(S_a)`$ と $`(T_a)`$ がスピン表現ならば
$`S_a\otimes1+1\otimes T_a`$ もそうである（交差項は可換）。したがって $`L`$ 個のスピン $`\frac12`$ の **全スピン** $`S_a=\sum_is_a^{(i)}`$
（`decompose.total_spin(L)`）はスピン表現である。回転 $`e^{-i\theta\,\mathbf n\cdot\mathbf S}`$（$`\theta\in\mathbb R`$、
$`\mathbf n\in\mathbb R^3`$ は単位ベクトル、関数計算）がすべての $`\theta,\mathbf n`$ について
$`X`$ と可換であることは、$`X`$ が三つの $`S_a`$ と可換であることと同値であり（$`\theta=0`$ で微分する、命題 8.2(1)。逆は定理 1.5(d)）、したがって「$`SU(2)`$ 作用」の可換子環は $`\{S_x,S_y,S_z\}'`$ であって、定理 2.2 により回転と $`S_a`$ は同じ環を生成する（どちらの集合も随伴をとる操作で閉じている）。$`L`$ 個のスピン $`\frac12`$ に対し、回転とは $`V^{\otimes L}`$、$`V\in SU(2)`$ のことである。実際、$`h=\mathbf n\cdot\mathbf s`$ として $`\mathbf n\cdot\mathbf S=\sum_ih^{(i)}`$ であり、4.4節の電荷の場合と同様に
（$`L`$ に関する帰納法）$`e^{-i\theta\,\mathbf n\cdot\mathbf S}=(e^{-i\theta h})^{\otimes L}`$ で、$`e^{-i\theta h}\in SU(2)`$ である（$`h`$ はトレースが零で $`h^2=\frac14`$ であるから、その固有値は
$`e^{\mp i\theta/2}`$ である）。逆に、$`V\in SU(2)`$ は絶対値 $`1`$ で積が $`1`$ の固有値をもつ正規行列であるから、ユニタリ $`W`$ により
$`V=W\mathrm{diag}(e^{i\varphi},e^{-i\varphi})W^*`$ と書け、トレースが零のエルミート行列 $`A=W\mathrm{diag}(\varphi,-\varphi)W^*`$ により $`V=e^{iA}`$ となる。トレースが零のエルミートな $`2\times2`$ 行列はパウリ行列の実線形結合
$`A=2\,\mathbf c\cdot\mathbf s`$ であるから、$`\theta=2|\mathbf c|`$、$`\mathbf n=-\mathbf c/|\mathbf c|`$ として $`V=e^{-i\theta\,\mathbf n\cdot\mathbf s}`$ である
（$`\mathbf c=0`$ ならば $`\theta=0`$）。定理 B.6(a),(c) と系 B.7（$`d=2`$）は、$`V^{\otimes L}`$、$`V\in SU(2)`$ が $`S_a`$ の生成する環を張ることを直接示す。

$`S_\pm=S_x\pm iS_y`$ とおくと $`S_-=S_+^*`$ である。また $`\mathbf S^2=S_x^2+S_y^2+S_z^2`$ とおく。関係式から
```math
[S_z,S_\pm]=\pm S_\pm,\qquad[S_+,S_-]=2S_z,\qquad\mathbf S^2=S_-S_++S_z^2+S_z,
```
が得られ、$`\mathbf S^2`$ はすべての $`S_a`$ と可換である（たとえば
$`[\mathbf S^2,S_z]=S_x[S_x,S_z]+[S_x,S_z]S_x+S_y[S_y,S_z]+[S_y,S_z]S_y=-i\{S_x,S_y\}+i\{S_y,S_x\}=0`$）。

<a id="r-4-11"></a>

**定理 4.11（既約スピン表現）.** 任意の $`s\in\{0,\frac12,1,\frac32,\dots\}`$ に対し、
`decompose.spin_matrices(s)` の行列——基底 $`|m\rangle`$、$`m=s,s-1,\dots,-s`$ をもつ $`\mathbb C^{2s+1}`$ 上の
```math
S_z|m\rangle=m|m\rangle,\qquad S_+|m\rangle=\sqrt{s(s+1)-m(m+1)}\;|m+1\rangle,\qquad S_-=S_+^*,
```
——は $`\mathbf S^2=s(s+1)1`$ を満たす既約なスピン表現（**スピン $`s`$**）をなす。任意の既約スピン表現は、これらのうちちょうど一つとユニタリ同値である。

*証明.* *スピン $`s`$ はスピン表現である.* $`[S_z,S_+]=S_+`$ は明らかである。$`S_+S_-|m\rangle=(s(s+1)-(m-1)m)|m\rangle`$
かつ $`S_-S_+|m\rangle=(s(s+1)-m(m+1))|m\rangle`$ であるから $`[S_+,S_-]=2S_z`$ である。これら二つの関係とその随伴は、三つの定義関係と同値である。$`\mathbf S^2|m\rangle=(s(s+1)-m(m+1)+m^2+m)|m\rangle=s(s+1)|m\rangle`$。
*既約性：* 不変な $`W\neq0`$ は $`S_z`$ で不変であり、$`S_z`$ の固有値は単純であるから、$`W`$ はある
$`|m\rangle`$ を含む。$`S_\pm`$ の係数は $`m=\pm s`$ でしか消えないので、$`S_\pm`$ を作用させればすべての $`|m'\rangle`$ が得られる。

*分類.* $`(S_a)`$ を $`\mathbb C^d`$ 上で既約とし、$`\lambda`$ を $`S_z`$ の最大固有値、
$`v_0`$ を単位固有ベクトルとする。$`S_zS_+v_0=S_+(S_z+1)v_0=(\lambda+1)S_+v_0`$ であるから、最大性により $`S_+v_0=0`$ である。
$`v_j=S_-^jv_0`$ とおく。すると $`S_zv_j=(\lambda-j)v_j`$ であり、$`S_+S_-=S_-S_++2S_z`$ を用いた帰納法により
```math
S_+v_j=j(2\lambda-j+1)\,v_{j-1}
```
である（$`S_+v_j=S_-S_+v_{j-1}+2S_zv_{j-1}=[(j-1)(2\lambda-j+2)+2(\lambda-j+1)]v_{j-1}`$）。零でない $`v_j`$ は相異なる固有値に属する
$`S_z`$ の固有ベクトルであるから、最後の零でないもの $`v_r`$ が存在する。$`v_0,\dots,v_r`$ の張る空間は
$`S_z,S_\pm`$ で不変であるから $`\mathbb C^d`$ に等しく、$`d=r+1`$ である。$`0=S_+v_{r+1}=(r+1)(2\lambda-r)v_r`$ から
$`\lambda=r/2=:s`$ を得る。最後に $`\|v_j\|^2=\langle v_{j-1},S_+v_j\rangle=j(2s-j+1)\|v_{j-1}\|^2`$ であるから、単位ベクトル
$`|s-j\rangle=v_j/\|v_j\|`$ は $`m=s-j`$ として $`S_+|m\rangle=\sqrt{(s-m)(s+m+1)}\,|m+1\rangle`$ を満たし、これは上の公式である。異なる $`s`$ は異なる次元をもつ。$`\square`$

*ライブラリでは：* スピン $`s`$ は `decompose.spin_matrices(s)` であり、[チュートリアル、4.6節](../../tutorial/ja/04_symmetries.md#s-4-6) で定理 4.11 と照合する；$`L`$ 個のスピン $`\frac12`$ の全スピンは `decompose.total_spin(L)` — [チュートリアル、4.5節](../../tutorial/ja/04_symmetries.md#s-4-5)。

<a id="r-4-12"></a>

**系 4.12.** $`\mathbb C^n`$ 上の任意のスピン表現は
$`\bigoplus_s(\text{spin }s)\otimes1_{m_s}`$ にユニタリ同値であり、$`d_M=\dim\ker(S_z-M)`$ とおくと
```math
m_s=d_s-d_{s+1}.
```
である。

*証明.* 定理 4.3 と定理 4.11 による。スピン $`s'`$ の各コピーは、$`|M|\le s'`$ かつ $`s'-M\in\mathbb Z`$ ならば固有値 $`M`$ に属する $`S_z`$ の固有ベクトルをちょうど一つ含み、そうでなければ含まない。よって
$`d_M=\sum_{s'\ge|M|,\,s'-M\in\mathbb Z}m_{s'}`$ であり、連続する二つの値の差が $`m_s`$ である。$`\square`$

<a id="r-4-13"></a>

**命題 4.13（$`L`$ 個のスピン $`\frac12`$）.** $`(\mathbb C^2)^{\otimes L}`$ 上の全スピンがスピン $`S`$ を含むのは
$`S\in\{\frac L2,\frac L2-1,\dots\}`$、$`S\ge0`$ のときかつそのときに限り、その重複度は
```math
m_S=\binom L{\frac L2-S}-\binom L{\frac L2-S-1}\qquad\Big(\tbinom L{-1}=0\Big),
```
である。したがって $`S_x,S_y,S_z`$ が生成する環は $`\bigoplus_SM_{2S+1}\otimes1_{m_S}`$ であり、その可換子環は
$`\bigoplus_S1_{2S+1}\otimes M_{m_S}`$ である。

*証明.* $`S_z=\sum_is_z^{(i)}`$ は積基底において対角であり、$`u`$ 個の上向きスピンと $`L-u`$ 個の下向きスピンの積の上で固有値 $`u-\frac L2`$ をもつから、$`\frac L2+M\in\{0,\dots,L\}`$ ならば $`d_M=\binom L{L/2+M}=\binom L{L/2-M}`$、そうでなければ $`0`$ である。系 4.12 と定理 4.3 を適用すればよい。値が正であるのは、
$`j<\frac{L-1}2`$ に対して $`\binom L{j+1}/\binom Lj=\frac{L-j}{j+1}>1`$ だからである。$`\square`$

$`L=2`$ の場合：$`\frac12\otimes\frac12=1\oplus0`$（**クレブシュ・ゴルダン**）。三重項 $`S_-^j|{\uparrow\uparrow}\rangle`$、$`j=0,1,2`$ は
二つのスピンの入れ替え $`F`$ のもとで対称であり、一重項
$`\frac1{\sqrt2}(|{\uparrow\downarrow}\rangle-|{\downarrow\uparrow}\rangle)`$ は反対称で、$`S_+`$ と $`S_z`$ によって零化される。
$`L=4`$ の場合：$`m_2=1`$、$`m_1=4-1=3`$、$`m_0=6-4=2`$ であるから
```math
\{S_x,S_y,S_z\}''\cong M_5\oplus M_3\oplus M_1\ \text{ with multiplicities }1,3,2,\qquad
\{S_x,S_y,S_z\}'\cong M_1\oplus M_3\oplus M_2,
```
であり、次元は $`35`$ と $`14`$ である——これは第10章で見いだされた数である。

*計算：* 三つの生成元だけから求めた $`L=4`$ の分解と、命題 4.13 およびスピン $`S`$ 上の $`\mathbf S^2=S(S+1)`$ — [チュートリアル、4.5節](../../tutorial/ja/04_symmetries.md#s-4-5)。

どの作用素がすべての回転と可換だろうか。サイトの置換はそうである：$`\pi(\sigma)`$ による共役は
$`s^{(i)}_a`$ を $`s^{(\sigma(i))}_a`$ に写す（補題 B.1(b)）ので、$`\sum_is^{(i)}_a`$ を固定する。逆に、それ以外には何もない：

<a id="r-4-14"></a>

**定理 4.14（シューア・ワイル双対性）.** $`\pi`$ を、$`(\mathbb C^2)^{\otimes L}`$ のテンソル因子を置換する $`S_L`$ の表現とする。(a) $`\mathrm{span}\,\pi(S_L)\subseteq\{S_x,S_y,S_z\}'`$。(b) 等号が成り立つ。同値なことに（定理 2.2）、全スピンが生成する環は $`\pi(S_L)'`$ である。

*証明.* (a) は今証明したところである。(b) は系 B.7、すなわち定理 B.6 の $`d=2`$ の場合である。定理 B.6 はすべての $`d`$ と $`L`$ について、
$`(\mathbb C^d)^{\otimes L}`$ と $`U(d)`$ の対角作用に対する双対性を証明している。$`\square`$

（さらに読むには：W. Fulton and J. Harris, *Representation Theory*, Springer 1991, Lecture 6; R. Goodman and
N. Wallach, *Symmetry, Representations, and Invariants*, Springer 2009, Section 9.1。）

*計算：* $`L\le5`$ に対する定理 4.14、$`\{S_a\}'`$ の次元 $`2,5,14,42`$（系 B.10）、および $`\pi(S_4)''`$（その組 $`(k_c,m_c)`$ は全スピンの組の役割を入れ替えたもの、系 B.9） — [チュートリアル、4.5節](../../tutorial/ja/04_symmetries.md#s-4-5)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. スピン $`S`$ の多重項は、任意の回転不変なハミルトニアンに対して等しいエネルギーをもつ $`2S+1`$ 個の状態からなり、それらは $`S^z=M`$ によって区別される。二つのスピン $`\frac12`$ は三重項と一重項をなし、
> 四つのスピンは一つの五重項、三つの三重項、二つの一重項をなす。二つの一重項の重複度空間と三つの三重項の重複度空間が、第10章のデコヒーレンスフリーな量子ビットとノイズレスなクォートリットである。

<a id="s-4-6"></a>

## 4.6 スピン系

$`L`$ 個のスピン $`\frac12`$ 上の **ハイゼンベルク模型** は、実数の結合定数 $`J_{ij}`$ により $`H=\sum_{i<j}J_{ij}\,\mathbf s^{(i)}\cdot\mathbf s^{(j)}`$ で与えられる。

<a id="r-4-15"></a>

**命題 4.15.** (a)（*ディラックの恒等式*）$`\mathbf s^{(i)}\cdot\mathbf s^{(j)}=\frac12F_{ij}-\frac14`$。ここで $`F_{ij}`$
はスピン $`i`$ と $`j`$ を入れ替える。したがって任意のハイゼンベルク・ハミルトニアンは $`\mathrm{span}\,\pi(S_L)\subseteq\{S_a\}'`$ に属し、
命題 4.13 の分解において $`V^*HV=\bigoplus_S1_{2S+1}\otimes H_S`$（ブロック $`H_S\in M_{m_S}`$）となる：$`H_S`$ の各固有値は $`H`$ の $`(2S+1)`$ 重に縮退した準位である。
(b) **最高ウェイトベクトル** の空間 $`\ker S_+\cap\ker(S_z-S)`$ は次元 $`m_S`$ をもち、$`H`$ で不変であり、
$`H`$ はその上に（基底のユニタリな変換を除いて）$`H_S`$ として作用する。

*証明.* (a) 二つのスピン上では $`\mathbf s^{(1)}\cdot\mathbf s^{(2)}=\frac12(\mathbf S^2-\mathbf s^{(1)2}-\mathbf s^{(2)2})=\frac12\mathbf S^2-\frac34`$
である（$`s=\frac12`$ の定理 4.11）。三重項のベクトル $`S_-^j|{\uparrow\uparrow}\rangle`$ は $`|{\uparrow\uparrow}\rangle`$、
$`|{\downarrow\uparrow}\rangle+|{\uparrow\downarrow}\rangle`$、$`2|{\downarrow\downarrow}\rangle`$ であり、一重項とあわせて $`\mathbb C^2\otimes\mathbb C^2`$ の直交基底をなす。
三重項上では $`\mathbf S^2=2`$ かつ $`F=1`$ である：$`\mathbf S^2=S_-S_++S_z^2+S_z`$ は $`|{\uparrow\uparrow}\rangle`$ 上で $`0+1+1=2`$ であり、
$`\mathbf S^2`$ と $`F`$ は $`S_-`$ と可換で、$`F`$ は $`|{\uparrow\uparrow}\rangle`$ を固定する。一重項上では
$`\mathbf S^2=0`$（$`S_+`$ と $`S_z`$ によって零化される）かつ $`F=-1`$ である。両辺はこれら四つのベクトル上で一致する。$`L`$ 個のうちのサイト $`i,j`$ については恒等作用素とのテンソル積をとる。残りは定理 4.14(a) と定理 4.3 である。(b) より多くを証明する：系 4.12 のように分解された任意のスピン表現と任意の $`\mu\in\frac12\mathbb Z`$ に対し、
$`\ker S_+\cap\ker(S_z-\mu)`$ は、$`\mu\ge0`$ ならば $`|\mu,\mu\rangle\otimes\mathbb C^{m_\mu}`$（の系 4.12 のユニタリによる像）であり、$`\mu<0`$ ならば $`0`$ である。実際、その基底において $`S_z=\mu`$ をもつベクトルは、
$`S'-\mu\in\mathbb Z`$ を満たす $`S'\ge|\mu|`$ にわたる和 $`\sum_{S'}|S',\mu\rangle\otimes w_{S'}`$ であり、定理 4.11 により $`S_+`$ はそれを
$`\sum_{S'}\sqrt{(S'-\mu)(S'+\mu+1)}\,|S',\mu+1\rangle\otimes w_{S'}`$ に写す。ここで $`S'+\mu+1\ge1`$ である。よってそれが消えるのは、
$`S'\neq\mu`$ なるすべての $`S'`$ に対して $`w_{S'}=0`$ となるときかつそのときに限る。ハイゼンベルクの $`H`$ に対し、$`1\otimes H_S`$ は $`|S,S\rangle\otimes\mathbb C^{m_S}`$ 上に $`H_S`$ として作用する。$`\square`$

*計算：* 四つのスピンの環と八つのスピンの環における命題 4.15 — [チュートリアル、4.6節](../../tutorial/ja/04_symmetries.md#s-4-6)。

(b) は、モジュール `su2`（スピン鎖、第9章）が数値的な分解なしに $`H_S`$ を構成する方法である：それは一度に一サイトずつ加えながら、
$`j\otimes\frac12`$ のクレブシュ・ゴルダン係数を用いて最高ウェイトベクトルの正規直交基底を構成する。これらの係数を最高ウェイトからの下降によって導こう。定理 4.11
と $`S_-=S_+^*`$ により、スピン $`j`$ は
```math
S_\pm|j,m\rangle=\sqrt{(j\mp m)(j\pm m+1)}\;|j,m\pm1\rangle ,
```
を満たす。$`|m|>j`$ に対しては $`|j,m\rangle=0`$ とおく（これは公式と整合的である。その係数は $`m=\pm j`$ で消える）。

<a id="r-4-16"></a>

**定理 4.16（$`j\otimes\frac12`$ のクレブシュ・ゴルダン係数）.** $`j\in\{\frac12,1,\frac32,\dots\}`$ とし、
$`\mathbb C^{2j+1}\otimes\mathbb C^2`$ 上で $`S_a=S^{(j)}_a\otimes1+1\otimes s_a`$ とする（基底 $`|j,m\rangle`$ をもつ定理 4.11 のスピン $`j`$ を、
$`|{\uparrow}\rangle=|\frac12,\frac12\rangle`$、$`|{\downarrow}\rangle=|\frac12,-\frac12\rangle`$ をもつスピン $`\frac12`$ と結合したもの）。また
$`\alpha_M=\sqrt{(j+M+\frac12)/(2j+1)}`$、$`\beta_M=\sqrt{(j-M+\frac12)/(2j+1)}`$ とおく。$`J=j\pm\frac12`$ と $`M=J,J-1,\dots,-J`$ に対し
```math
\begin{aligned}
|j+\tfrac12,M\rangle\!\rangle&=\alpha_M\,|j,M-\tfrac12\rangle\otimes|{\uparrow}\rangle+\beta_M\,|j,M+\tfrac12\rangle\otimes|{\downarrow}\rangle,\\
|j-\tfrac12,M\rangle\!\rangle&=-\beta_M\,|j,M-\tfrac12\rangle\otimes|{\uparrow}\rangle+\alpha_M\,|j,M+\tfrac12\rangle\otimes|{\downarrow}\rangle .
\end{aligned}
```
とおく。
(a) これら $`2(2j+1)`$ 個のベクトルは $`\mathbb C^{2j+1}\otimes\mathbb C^2`$ の正規直交基底である。
(b) $`|j+\frac12,j+\frac12\rangle\!\rangle=|j,j\rangle\otimes|{\uparrow}\rangle`$ であり、$`|j-\frac12,j-\frac12\rangle\!\rangle`$ は、絶対値 $`1`$ の因子を除いて、
$`S_zx=(j-\frac12)x`$ かつ $`S_+x=0`$ を満たす唯一の単位ベクトル $`x`$ である。この因子は
$`|j,j\rangle\otimes|{\downarrow}\rangle`$ の係数が正であることによって定まる（コンドン・ショートレーの規約）。
(c) $`S_z|J,M\rangle\!\rangle=M|J,M\rangle\!\rangle`$ かつ $`S_-|J,M\rangle\!\rangle=\sqrt{(J+M)(J-M+1)}\,|J,M-1\rangle\!\rangle`$ である：状態 $`|J,M\rangle\!\rangle`$
は $`|J,J\rangle\!\rangle`$ から下降と規格化によって得られる。
(d) したがって $`a=x,y,z`$ に対し $`S_a|J,M\rangle\!\rangle=\sum_{M'}(S^{(J)}_a)_{M'M}|J,M'\rangle\!\rangle`$ である：ユニタリ $`|J,M\rangle\mapsto|J,M\rangle\!\rangle`$ は
スピン $`(j+\frac12)\oplus`$ スピン $`(j-\frac12)`$ を $`S_a`$ に写すので、$`j\otimes\frac12\cong(j+\frac12)\oplus(j-\frac12)`$ であり、
クレブシュ・ゴルダン係数 $`\langle j,m;\frac12,\sigma|J,M\rangle`$（積基底における $`|J,M\rangle\!\rangle`$ の座標）は上に示した
$`\pm\alpha_M,\beta_M`$ である。

*証明.* (a) 上の式の基底ベクトルが $`0`$ ならば、その係数も $`0`$ である：$`|j,M-\frac12\rangle=0`$ となるのは
$`M=-j-\frac12`$ のときだけで、このとき $`\alpha_M=0`$ であり、$`|j,M+\frac12\rangle=0`$ となるのは $`M=j+\frac12`$ のときだけで、このとき $`\beta_M=0`$ である。よって各ベクトルのノルムの二乗は
$`\alpha_M^2+\beta_M^2=1`$ である。$`M`$ が異なるベクトルは異なる積基底ベクトルを含み、$`M`$ が等しいときの
内積は $`-\alpha_M\beta_M+\beta_M\alpha_M=0`$ である。それらは $`(2j+2)+2j=2(2j+1)`$ 個ある。
(b) 第一の主張は $`\alpha_{j+1/2}=1`$、$`\beta_{j+1/2}=0`$ である。固有空間 $`S_z=j-\frac12`$ は
$`|j,j-1\rangle\otimes|{\uparrow}\rangle`$ と $`|j,j\rangle\otimes|{\downarrow}\rangle`$ で張られ、$`s_+|{\uparrow}\rangle=0`$、$`s_+|{\downarrow}\rangle=|{\uparrow}\rangle`$ により
```math
S_+\big(p\,|j,j-1\rangle\otimes|{\uparrow}\rangle+q\,|j,j\rangle\otimes|{\downarrow}\rangle\big)=(p\sqrt{2j}+q)\,|j,j\rangle\otimes|{\uparrow}\rangle ,
```
である。これが消えるのは $`q=-\sqrt{2j}\,p`$ のときかつそのときに限る。$`|p|^2+|q|^2=1`$ と $`q>0`$（したがって $`p<0`$）のもとでは、これは $`p=-\beta_{j-1/2}=-1/\sqrt{2j+1}`$、
$`q=\alpha_{j-1/2}=\sqrt{2j/(2j+1)}`$ である。(c) $`S_z`$：どちらの項も $`S_z=(M\mp\frac12)\pm\frac12=M`$ をもつ。$`S_-`$ については、
$`x_M=p\,|j,M-\frac12\rangle\otimes|{\uparrow}\rangle+q\,|j,M+\frac12\rangle\otimes|{\downarrow}\rangle`$ とおく。$`s_-|{\uparrow}\rangle=|{\downarrow}\rangle`$、$`s_-|{\downarrow}\rangle=0`$ であるから
```math
S_-x_M=p\sqrt{(j+M-\tfrac12)(j-M+\tfrac32)}\;|j,M-\tfrac32\rangle\otimes|{\uparrow}\rangle
+\Big(p+q\sqrt{(j+M+\tfrac12)(j-M+\tfrac12)}\Big)|j,M-\tfrac12\rangle\otimes|{\downarrow}\rangle .
```
である。$`\beta_M\sqrt{(j+M+\frac12)(j-M+\frac12)}=(j-M+\frac12)\alpha_M`$ と $`\alpha_M\sqrt{(j+M+\frac12)(j-M+\frac12)}=(j+M+\frac12)\beta_M`$ を用いる。
$`J=j+\frac12`$ に対しては $`(p,q)=(\alpha_M,\beta_M)`$、$`c=\sqrt{(J+M)(J-M+1)}=\sqrt{(j+M+\frac12)(j-M+\frac32)}`$ であり、第一の係数は
$`\sqrt{(j+M+\frac12)(j+M-\frac12)(j-M+\frac32)/(2j+1)}=c\,\alpha_{M-1}`$、第二の係数は $`\alpha_M(j-M+\frac32)=c\,\beta_{M-1}`$ である。
$`J=j-\frac12`$ に対しては $`(p,q)=(-\beta_M,\alpha_M)`$、$`c=\sqrt{(j+M-\frac12)(j-M+\frac12)}`$ であり、第一の係数は
$`-\sqrt{(j-M+\frac12)(j+M-\frac12)(j-M+\frac32)/(2j+1)}=-c\,\beta_{M-1}`$、第二の係数は $`\beta_M(j+M-\frac12)=c\,\alpha_{M-1}`$ である。
いずれの場合も $`M>-J`$ に対して $`S_-|J,M\rangle\!\rangle=c\,|J,M-1\rangle\!\rangle`$ である。$`M=-J`$ に対しては両方の係数が消える（$`J=j+\frac12`$ ならば $`p=\alpha_M=0`$、
$`J=j-\frac12`$ ならば $`j+M-\frac12=0`$）ので、$`S_-|J,-J\rangle\!\rangle=0`$ である。(d) (a) と (c) により、正規直交基底
$`|J,M\rangle\!\rangle`$ における $`S_z`$ と $`S_-`$ の行列はスピン $`(j+\frac12)\oplus`$ スピン $`(j-\frac12)`$ のもの（定理 4.11
と上の公式）である。$`S_+=S_-^*`$ の行列も同様であり、したがって $`S_x=\frac12(S_++S_-)`$ と $`S_y=\frac1{2i}(S_+-S_-)`$ の行列も同様である。$`\square`$

*計算：* $`j\le2`$ に対する定理 4.16 を、密な行列および `su2.clebsch_gordan` と照合 — [チュートリアル、4.6節](../../tutorial/ja/04_symmetries.md#s-4-6)。

演習問題 2 はこの分解を（係数は除いて）$`j_1\otimes j_2`$ に拡張する。ライブラリにとって重要なのは定理 4.16 の二つの
最高ウェイトベクトルだけであり、それらは既約であるか否かを問わず任意のスピン表現において機能する：

<a id="r-4-17"></a>

**系 4.17（一度に一サイトずつの最高ウェイトベクトル）.** $`(T_a)`$ を $`E`$ 上のスピン表現とし、
$`\mathcal H_\mu(E)=\ker T_+\cap\ker(T_z-\mu)`$ とおき、その正規直交基底を $`B_\mu`$ とする（$`\mathcal H_\mu(E)=0`$ ならば $`B_\mu=\emptyset`$）。また
$`E\otimes\mathbb C^2`$ 上で $`S_a=T_a\otimes1+1\otimes s_a`$ とする。任意の $`S\in\{0,\frac12,1,\dots\}`$ に対し、ベクトル
```math
v\otimes|{\uparrow}\rangle\quad(v\in B_{S-1/2}),\qquad
\sqrt{\tfrac{2S+1}{2S+2}}\;v\otimes|{\downarrow}\rangle-\tfrac1{\sqrt{(2S+1)(2S+2)}}\;T_-v\otimes|{\uparrow}\rangle\quad(v\in B_{S+1/2})
```
は $`\ker S_+\cap\ker(S_z-S)`$ の正規直交基底をなす（$`S=0`$ では第一の族は空である）。一つのスピン
（$`E=\mathbb C^2`$、$`B_{1/2}=\{|{\uparrow}\rangle\}`$、他のすべての $`B_\mu`$ は空）から始めて $`L-1`$ 個のスピンを加えていくと、
$`L`$ 個のスピン $`\frac12`$ の最高ウェイトベクトルの空間の、次元 $`m_S`$ の正規直交基底が得られる（命題 4.15(b)）。

*証明.* $`a=\sqrt{(2S+1)/(2S+2)}`$、$`b=-1/\sqrt{(2S+1)(2S+2)}`$、$`u(v)=a\,v\otimes|{\downarrow}\rangle+b\,T_-v\otimes|{\uparrow}\rangle`$ と書く。
$`v\in\mathcal H_{S+1/2}(E)`$ に対し $`T_+T_-v=(T_-T_++2T_z)v=(2S+1)v`$ であり、$`a+(2S+1)b=0`$ である。*最高ウェイト：* $`S_z`$ はどちらの族の上でも $`S`$ である
（$`v`$ 上で $`T_z=S\mp\frac12`$、$`T_-v`$ 上で $`S-\frac12`$）。$`S_+(v\otimes|{\uparrow}\rangle)=T_+v\otimes|{\uparrow}\rangle=0`$ であり、
$`S_+u(v)=a\,T_+v\otimes|{\downarrow}\rangle+a\,v\otimes|{\uparrow}\rangle+b\,T_+T_-v\otimes|{\uparrow}\rangle=(a+(2S+1)b)\,v\otimes|{\uparrow}\rangle=0`$ である。
*正規直交性：* $`v\in B_{S-1/2}`$ に対し $`\langle v\otimes|{\uparrow}\rangle,u(v')\rangle=b\langle v,T_-v'\rangle=b\langle T_+v,v'\rangle=0`$ であり、
$`\langle u(v),u(v')\rangle=a^2\langle v,v'\rangle+b^2\langle v,T_+T_-v'\rangle=(a^2+(2S+1)b^2)\langle v,v'\rangle=\langle v,v'\rangle`$ である。
*張ること：* $`S_zw=Sw`$ かつ $`S_+w=0`$ とし、$`w=x\otimes|{\uparrow}\rangle+y\otimes|{\downarrow}\rangle`$ と書く。すると $`T_zx=(S-\frac12)x`$、
$`T_zy=(S+\frac12)y`$ であり、$`0=S_+w=(T_+x+y)\otimes|{\uparrow}\rangle+T_+y\otimes|{\downarrow}\rangle`$ であるから $`T_+y=0`$、すなわち $`y\in\mathcal H_{S+1/2}(E)`$
で、$`T_+x=-y`$ である。すると $`x'=x+\frac1{2S+1}T_-y`$ は $`T_+x'=-y+y=0`$ を満たすので $`x'\in\mathcal H_{S-1/2}(E)`$ であり、これは $`S=0`$ では $`0`$ である
（命題 4.15(b) の証明）。$`b/a=-1/(2S+1)`$ であるから、$`w=x'\otimes|{\uparrow}\rangle+a^{-1}u(y)`$（$`u`$ を線形に拡張する）は
二つの族の張る空間に属する。$`L`$ 個のスピンについては、$`n=2,\dots,L`$ に対して $`E=(\mathbb C^2)^{\otimes(n-1)}`$ とその全スピン $`T_a`$ にこれを適用する。
一つのスピン上では $`\mathcal H_{1/2}=\mathbb C|{\uparrow}\rangle`$ であり、他のすべての $`\mathcal H_\mu=0`$ である（スピンは $`\frac12`$ だけである）。$`\square`$

第二の族は、$`j=S+\frac12`$ とした定理 4.16 のベクトル $`|S,S\rangle\!\rangle`$ を、$`T_-^kv`$ の張るスピン $`j`$ のコピーの中で書いたものであり
（定理 4.11 の証明、$`|j,j-1\rangle=T_-v/\sqrt{2j}`$）、第一の族は $`j=S-\frac12`$ とした
$`|j+\frac12,j+\frac12\rangle\!\rangle`$ である。各基底ベクトルに中間スピンの列
$`S_1=\frac12,S_2,\dots,S_L=S`$、$`|S_n-S_{n-1}|=\frac12`$ でラベルをつけると、ヤマノウチ基底が得られる。

*ライブラリでは：* `su2.SU2Chain` はこれらの係数を用いて、`SpinChain` の $`U(1)`$ セクター $`\hat N=\frac L2+S`$ の内部にヤマノウチ基底を構成する；系 4.17 を八つのスピンについて確かめ、ハイゼンベルク環の $`U(1)`$ ブロックと $`SU(2)`$ ブロックおよびその計算コストを示す — [チュートリアル、4.6節](../../tutorial/ja/04_symmetries.md#s-4-6)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. ハイゼンベルク模型は、磁性絶縁体における局在電子スピンの交換相互作用を記述する。$`J_{ij}>0`$ は反平行な隣接スピンを好む。そのスペクトルは
> 多重項からなり、磁場 $`h\,S^z`$ は各多重項を $`2S+1`$ 個の等間隔な準位に分裂させる（演習問題 4）。
> 定理 9.9（二部的な結合。偶数長の反強磁性の鎖と環では一重項の基底状態）は、E. Lieb と D. Mattis
> （*J. Math. Phys.* 3 (1962) 749）にさかのぼる。

<a id="s-4-7"></a>

## 4.7 実代数：実型、複素型、四元数型

対称性はしばしば実行列で与えられる：$`\mathbb R^3`$ の回転、置換、ある基底で実であるハミルトニアンなど。それらが生成する複素代数は定理 2.6 によって分解されるが、その中の *実* 行列は実代数をなし、その単純な成分は $`M_k(\mathbb R)`$、$`M_k(\mathbb C)`$、$`M_k(\mathbb H)`$ の三種類である。
これらを複素分解ともう一つのデータ、すなわち複素共役がセクターに何をするかから導く。

$`a\in M_n(\mathbb C)`$ に対し、$`\bar a`$ を成分ごとの複素共役とする（ベクトルについても同様であり、集合については $`\bar S=\{\bar s:s\in S\}`$ とする）。
共役は反線形かつ全単射であり、$`\overline{ab}=\bar a\bar b`$、$`\overline{a^*}=\bar a^*`$、$`\bar{\bar a}=a`$ である；$`\bar a=a`$ であるのは $`a`$ が実であるときかつそのときに限り、そのとき $`a^*=a^T`$ である。行列（またはベクトル）の集合
$`\mathcal B`$ に対し $`\mathcal B_{\mathbb R}=\mathcal B\cap M_n(\mathbb R)`$（または
$`\mathcal B\cap\mathbb R^n`$）とおく。内積は $`\langle u,v\rangle=\sum_i\bar u_iv_i`$ であるから、
$`\langle\bar u,\bar v\rangle=\overline{\langle u,v\rangle}=\langle v,u\rangle`$ である。

<a id="r-4-18"></a>

**補題 4.18（実部）.** (a) $`\mathcal B`$ を $`M_n(\mathbb C)`$ または $`\mathbb C^n`$ の複素部分空間で
$`\bar{\mathcal B}=\mathcal B`$ を満たすものとする。このとき $`\mathcal B=\mathcal B_{\mathbb R}\oplus i\mathcal B_{\mathbb R}`$ かつ
$`\dim_{\mathbb R}\mathcal B_{\mathbb R}=\dim_{\mathbb C}\mathcal B`$ である。
(b) $`\mathcal A`$ を $`\bar{\mathcal A}=\mathcal A`$ を満たす $`\mathbb C^n`$ 上の $`*`$-環とする。このとき $`\mathcal A_{\mathbb R}`$ は
単位元をもち転置で閉じた実代数である；$`\overline{\mathcal A'}=\mathcal A'`$ である；そして $`(\mathcal A')_{\mathbb R}`$ は
すべての $`a\in\mathcal A_{\mathbb R}`$ に対し $`Xa=aX`$ を満たす実行列 $`X`$ 全体である。
(c) $`S\subseteq M_n(\mathbb C)`$ が $`S^*=S`$ と $`\bar S=S`$ を満たすならば——たとえば実行列 $`S_j`$ による $`S=\{S_1,\dots,S_J,S_1^T,\dots,S_J^T\}`$、
あるいは実直交行列からなる群——$`\overline{S''}=S''`$ である。

*証明.* (a) $`b\in\mathcal B`$ に対し、$`\mathrm{Re}\,b=\frac12(b+\bar b)`$ と $`\mathrm{Im}\,b=\frac1{2i}(b-\bar b)`$ は
$`\mathcal B`$ の実な元であって $`b=\mathrm{Re}\,b+i\mathrm{Im}\,b`$ であり、
$`\mathcal B_{\mathbb R}\cap i\mathcal B_{\mathbb R}`$ の元は実かつ純虚であるから $`0`$ である。したがって実空間 $`\mathcal B_{\mathbb R}`$ の基底 $`b_1,\dots,b_r`$ は
$`\mathbb C`$ 上で $`\mathcal B`$ を張り、しかも $`\mathbb C`$ 上線形独立である：実数 $`\alpha_j,\beta_j`$ により $`\sum_j(\alpha_j+i\beta_j)b_j=0`$ ならば、実部と虚部から
$`\sum_j\alpha_jb_j=\sum_j\beta_jb_j=0`$ が得られる。(b) $`\mathcal A_{\mathbb R}`$ は $`1`$ を含み、積について閉じ、
$`a\mapsto a^T=a^*`$ についても閉じている。$`Y`$ がすべての $`t\in T`$ と可換であることは $`\bar Y`$ がすべての $`\bar t`$ と可換であることと同値であるから、任意の集合 $`T`$ に対し $`\overline{T'}=\bar T'`$
であり、$`\overline{\mathcal A'}=\bar{\mathcal A}'=\mathcal A'`$ である。$`\mathcal A_{\mathbb R}`$ と可換な実行列 $`X`$ は
$`\mathcal A_{\mathbb R}+i\mathcal A_{\mathbb R}=\mathcal A`$ と可換である（(a) による）。(c) $`\overline{T'}=\bar T'`$ を二回適用する：$`\overline{S''}=\bar S''=S''`$。$`\square`$

**記号の約束.** 本節の残りでは、$`\mathcal A`$ は $`\bar{\mathcal A}=\mathcal A`$ を満たす $`\mathbb C^n`$ 上の $`*`$-環であり、$`V`$ は定理 2.6 のような分解である。$`V_c`$ は $`c`$ 番目の直和成分に属する $`V`$ の $`n\times k_cm_c`$ の列ブロックであり、したがって $`V_c^*V_d=\delta_{cd}1`$、
$`V_cV_c^*=P_c`$（極小中心射影）であって、
```math
aV_c=V_c\big(\pi_c(a)\otimes1_{m_c}\big)\quad(a\in\mathcal A),\qquad XV_c=V_c\big(1_{k_c}\otimes x_c\big)\quad(X\in\mathcal A'),
```
ここで $`\pi_c(a)\in M_{k_c}`$、$`x_c\in M_{m_c}`$ である。各 $`\pi_c`$ は $`\mathcal A`$ から
$`M_{k_c}`$ の上への単位的 $`*`$-準同型である：$`x\in M_{k_c}`$ に対し、$`\mathcal A`$ の元 $`V_c(x\otimes1)V_c^*`$ は $`\pi_c=x`$ をもち、$`e\neq c`$ に対し $`\pi_e=0`$ をもつ。
最後に $`\Omega_2=\begin{pmatrix}0&-1\\1&0\end{pmatrix}`$、$`\Omega_{2r}=1_r\otimes\Omega_2`$（$`\Omega_2`$ の $`r`$ 個のコピーからなるブロック対角行列）とおく。

部分空間 $`E\subseteq\mathbb C^k`$ 上の写像 $`\theta:E\to E`$ は、
$`\theta(\alpha u+\beta v)=\bar\alpha\theta u+\bar\beta\theta v`$ かつ $`\langle\theta u,\theta v\rangle=\langle v,u\rangle`$ であるとき **反ユニタリ** であるという。
ユニタリ $`J\in U(k)`$ に対し、$`\theta_Jv=J\bar v`$ は $`\mathbb C^k`$ 上で反ユニタリであり
（$`\langle J\bar u,J\bar v\rangle=\langle\bar u,\bar v\rangle=\langle v,u\rangle`$）、$`\theta_J^2v=J\bar Jv`$ である。

<a id="r-4-19"></a>

**補題 4.19（反ユニタリ写像）.** $`\theta`$ を $`E\neq0`$ 上の反ユニタリ写像で $`\theta^2=\varepsilon\,\mathrm{id}`$、
$`\varepsilon\in\{1,-1\}`$ を満たすものとする。
(a) $`\varepsilon=1`$ ならば、$`E`$ は $`\theta w=w`$ を満たすベクトル $`w`$ からなる正規直交基底をもつ。
(b) $`\varepsilon=-1`$ ならば、すべての $`u\in E`$ に対し $`\langle u,\theta u\rangle=0`$ である；任意の単位ベクトル $`u`$ に対し $`\mathrm{span}\{u,\theta u\}`$
は $`\theta`$ 不変な平面である；そして $`E`$ は正規直交基底 $`w_1,\theta w_1,\dots,w_r,\theta w_r`$ をもつ。特に
$`\dim E=2r`$ は偶数である。
(c) $`J\in U(k)`$ が $`J\bar J=\varepsilon1`$ を満たすとする。$`\varepsilon=1`$ ならば、$`J\bar W=W`$、すなわち
$`J=WW^T`$ を満たすユニタリ $`W`$ が存在する。$`\varepsilon=-1`$ ならば、$`k`$ は偶数であり、$`J\bar W=W\Omega_k`$、すなわち $`J=W\Omega_kW^T`$ を満たすユニタリ $`W`$ が存在する。

*証明.* (a) $`F=\{v\in E:\theta v=v\}`$ とおく。これは実部分空間である。任意の $`v\in E`$ は
$`v_1=\frac12(v+\theta v)`$、$`v_2=\frac1{2i}(v-\theta v)`$ により $`v_1+iv_2`$ と書け、どちらも $`F`$ に属する
（$`\theta v_2=\frac1{-2i}(\theta v-v)=v_2`$）。$`u,v\in F`$ に対し $`\langle u,v\rangle=\langle\theta u,\theta v\rangle=\langle v,u\rangle`$ は
実数であるから、$`\langle\cdot,\cdot\rangle`$ は $`F`$ 上の実内積である；$`w_1,\dots,w_s`$ を $`F`$ の正規直交基底とする
（グラム–シュミット）。$`w_j`$ は $`E`$ において正規直交であるから $`\mathbb C`$ 上線形独立であり、$`E=F+iF`$ であるから $`\mathbb C`$ 上で $`E`$ を張る。(b) $`\langle\theta x,\theta y\rangle=\langle y,x\rangle`$ を $`x=\theta u`$、$`y=u`$ に適用すると
$`\langle u,\theta u\rangle=\langle\theta^2u,\theta u\rangle=-\langle u,\theta u\rangle`$ である。よって単位ベクトル $`u`$ に対し $`\theta u`$ は $`u`$ に直交する
単位ベクトルであり、$`E_1=\mathrm{span}\{u,\theta u\}`$ は不変である：
$`\theta(\alpha u+\beta\theta u)=\bar\alpha\theta u-\bar\beta u`$。$`E`$ における $`E_1`$ の直交補空間も不変である：$`x\perp E_1`$ と $`y\in E_1`$ に対し $`y=-\theta(\theta y)`$ であり、
$`\theta y\in E_1`$ であるから $`\langle\theta x,y\rangle=-\langle\theta x,\theta(\theta y)\rangle=-\langle\theta y,x\rangle=0`$ である。
$`\dim E`$ に関する帰納法で基底が得られる。(c) $`\theta_J`$ は $`\mathbb C^k`$ 上の反ユニタリ写像で $`\theta_J^2=J\bar J=\varepsilon`$ である。$`W`$ を、
(a) の基底、あるいは (b) の基底 $`w_1,\theta_Jw_1,w_2,\theta_Jw_2,\dots`$ を列とするユニタリとする。$`J\bar W`$ の
列は $`W`$ の列の $`\theta_J`$ による像である。場合 (a) ではこれは $`J\bar W=W`$ を意味する。場合
(b) では像は $`\theta_Jw_j`$ と $`\theta_J^2w_j=-w_j`$ であり、これらは $`W\Omega_k`$ の第 $`2j-1`$ 列と第 $`2j`$ 列である
（$`\Omega_2e_1=e_2`$、$`\Omega_2e_2=-e_1`$ であるから）。$`\bar W^{-1}=\overline{W^*}=W^T`$ であるから、$`J=WW^T`$、それぞれ $`J=W\Omega_kW^T`$ である。$`\square`$

<a id="r-4-20"></a>

**定理 4.20（セクターの共役）.** (a) $`\sigma\circ\sigma=\mathrm{id}`$ を満たすセクターの置換 $`\sigma`$ と、ユニタリ
$`J_c\in U(k_c)`$、$`M_c\in U(m_c)`$ が存在して、$`\bar P_c=P_{\sigma(c)}`$、
$`k_{\sigma(c)}=k_c`$、$`m_{\sigma(c)}=m_c`$ であり、すべての $`c`$ とすべての $`a\in\mathcal A`$ に対し
```math
\overline{V_c}=V_{\sigma(c)}\,(J_c\otimes M_c),\qquad\pi_{\sigma(c)}(\bar a)=J_c\,\overline{\pi_c(a)}\,J_c^*.
```
となる。$`J_c\otimes M_c`$ は一意であり、$`J_c,M_c`$ は $`|\lambda|=1`$ による $`(J_c,M_c)\mapsto(\lambda J_c,\bar\lambda M_c)`$ を除いて一意である。
(b) $`\sigma(c)=c`$ ならば、$`J_c\bar J_c=\varepsilon_c1_{k_c}`$ かつ
$`M_c\bar M_c=\varepsilon_c1_{m_c}`$ となる $`\varepsilon_c\in\{1,-1\}`$ が存在する。$`\varepsilon_c=-1`$ ならば $`k_c`$ と $`m_c`$ は偶数である。
(c) $`\sigma`$ と $`\varepsilon_c`$ は $`V`$ の選び方によらない。
(d) $`y\in M_{m_c}`$ に対し $`\overline{V_c(1\otimes y)V_c^*}=V_{\sigma(c)}\big(1\otimes M_c\bar yM_c^*\big)V_{\sigma(c)}^*`$ である。したがって、
$`X\in\mathcal A'`$ のブロックが $`x_c`$ ならば、$`\bar X\in\mathcal A'`$ のブロックは $`(\bar X)_{\sigma(c)}=M_c\bar x_cM_c^*`$ である：
可換子環の上では、$`M_c`$ は $`J_c`$ が $`\mathcal A`$ の上で果たす役割を果たす。

*証明.* (a) 共役は $`\mathcal A`$ を $`\mathcal A`$ の上へ写し、和、積、随伴を保つので、
中心 $`Z(\mathcal A)`$ をそれ自身の上へ写し、$`Z(\mathcal A)`$ の零でない射影全体をそれ自身の上へ全単射に写し、
順序 $`p\le q`$（$`pq=p`$）を保つ。したがって極小なものを置換し、それらは $`P_c`$ である（補題 2.3）：
$`\bar P_c=P_{\sigma(c)}`$ であり、$`\bar{\bar P}_c=P_c`$ であるから $`\sigma\circ\sigma=\mathrm{id}`$ である。$`d=\sigma(c)`$ とおく。
$`\overline{V_c}`$ の値域は $`P_c`$ の値域の共役であり、それは $`\bar P_c=P_d`$ の値域である；よって
$`T=V_d^*\overline{V_c}`$ として $`\overline{V_c}=P_d\overline{V_c}=V_dT`$ であり、
$`T^*T=\overline{V_c}^*P_d\overline{V_c}=\overline{V_c^*V_c}=1`$ である。$`k_cm_c=\mathrm{rank}\,P_c=\mathrm{rank}\,\bar P_c=k_dm_d`$ であるから、
等長写像 $`T`$ は正方であり、したがってユニタリである。$`aV_c=V_c(\pi_c(a)\otimes1)`$ の共役をとり、
$`\bar aV_d=V_d(\pi_d(\bar a)\otimes1)`$ を用いると $`V_d(\pi_d(\bar a)\otimes1)T=V_dT(\overline{\pi_c(a)}\otimes1)`$ を得る；左から
$`V_d^*`$ を掛けると
```math
\big(\pi_d(\bar a)\otimes1_{m_d}\big)\,T=T\,\big(\overline{\pi_c(a)}\otimes1_{m_c}\big)\qquad(a\in\mathcal A).
```
となる。$`k_d\times k_c`$ ブロック $`T_{jl}`$ により $`T=\sum_{j,l}T_{jl}\otimes|j\rangle\langle l|`$ と書く。上の式は
$`\pi_d(\bar a)T_{jl}=T_{jl}\overline{\pi_c(a)}`$ を意味する：各 $`T_{jl}`$ は $`\rho(a)=\overline{\pi_c(a)}`$ から
$`\rho'(a)=\pi_d(\bar a)`$ への絡作用素である。これらは随伴をもつ $`\mathcal A`$ の 4.1節の意味での表現であり
（$`\rho(a^*)=\rho(a)^*`$、$`\rho'(a^*)=\rho'(a)^*`$）、その像はそれぞれ $`M_{k_c}`$ 全体と
$`M_{k_d}`$ 全体であるから既約である。$`T\neq0`$ であるから、ある $`T_{jl}\neq0`$ であり、したがって $`k_c=k_d`$ であって、あるユニタリ $`J_c`$ により $`\mathrm{Hom}(\rho,\rho')=\mathbb CJ_c`$ である
（補題 4.2）。よって $`T_{jl}=\mu_{jl}J_c`$、すなわち $`M_c=(\mu_{jl})`$ として $`T=J_c\otimes M_c`$ である；$`1=T^*T=1\otimes M_c^*M_c`$ かつ
$`m_c=k_cm_c/k_c=m_d`$ であるから $`M_c`$ はユニタリであり、上の式は $`\pi_d(\bar a)`$ の公式になる。$`T`$ は
$`V`$ によって決まる；$`J'`$ をユニタリとして $`T=J'\otimes M'`$ とも書けるならば、ブロック $`T_{jl}=M'_{jl}J'`$ は $`J'\in\mathrm{Hom}(\rho,\rho')`$ を示すので、
$`|\lambda|=1`$ により $`J'=\lambda J_c`$ であり、そのとき $`M'=\bar\lambda M_c`$ である。

(b) $`\sigma(c)=c`$ とし、$`J=J_c`$、$`M=M_c`$ とおく。$`\overline{V_c}=V_c(J\otimes M)`$ の共役をとると
$`V_c=\overline{V_c}(\bar J\otimes\bar M)=V_c(J\bar J\otimes M\bar M)`$ であるから、$`J\bar J\otimes M\bar M=1`$ である。(a) の公式を
二回適用し、$`\overline{J^*}=J^T`$ を用いると
```math
\pi_c(a)=\pi_c(\bar{\bar a})=J\,\overline{J\,\overline{\pi_c(a)}\,J^*}\,J^*=(J\bar J)\,\pi_c(a)\,(J\bar J)^*,
```
であるから、ユニタリ $`J\bar J`$ は $`\pi_c(\mathcal A)=M_{k_c}`$ と可換であり、スカラー $`\lambda1`$ である（補題 2.5 のステップ5、
$`m=1`$）、$`|\lambda|=1`$。すると $`\bar JJ=\overline{J\bar J}=\bar\lambda1`$ であり、$`\lambda J=(J\bar J)J=J(\bar JJ)=\bar\lambda J`$ であるから、
$`\lambda=\bar\lambda\in\{1,-1\}`$ である；これを $`\varepsilon_c`$ とよぶ。これは $`J\mapsto\lambda J`$ で変わらず、
$`\varepsilon_c1\otimes M\bar M=1`$ から $`M\bar M=\varepsilon_c1`$ が得られる。$`\varepsilon_c=-1`$ ならば
$`(-1)^{k_c}=\det(J\bar J)=\det J\cdot\overline{\det J}=|\det J|^2>0`$ であるから $`k_c`$ は偶数であり、$`m_c`$ も同様である。

(c) $`P_c`$ は $`V`$ によらない（定理 2.6）から、$`\sigma`$ もよらない。$`V'`$ を別の分解とし、
$`V'_cV_c'^*=P_c`$ となるように番号を付け（定理 2.6）、準同型を $`\pi'_c`$、データを $`J'_c,M'_c`$ とし、
$`\sigma(c)=c`$ とする。このとき $`Y=V_c^*V'_c`$ として $`V'_c=P_cV'_c=V_cY`$ であり、$`Y`$ は正方の等長写像であるからユニタリであって、
$`a\in\mathcal A`$ に対し $`(\pi_c(a)\otimes1)Y=Y(\pi'_c(a)\otimes1)`$ である。$`\rho,\rho'`$ の代わりに既約表現 $`\pi'_c,\pi_c`$ を用いた (a) の議論により、ユニタリ $`u,w`$ によって $`Y=u\otimes w`$ となる。したがって
$`\overline{V'_c}=\overline{V_c}(\bar u\otimes\bar w)=V_c(J_c\bar u\otimes M_c\bar w)=V'_c(u^*J_c\bar u\otimes w^*M_c\bar w)`$ であるから、
(a) の一意性により $`|\lambda|=1`$ として $`J'_c=\lambda u^*J_c\bar u`$ であり、$`\bar uu^T=\overline{uu^*}=1`$ であるから
$`J'_c\bar J'_c=u^*J_c\bar uu^T\bar J_cu=u^*J_c\bar J_cu=\varepsilon_c1`$ である。

(d) (a) により $`\overline{V_c(1\otimes y)V_c^*}=V_{\sigma(c)}(J_c\otimes M_c)(1\otimes\bar y)(J_c\otimes M_c)^*V_{\sigma(c)}^*`$ であり、
$`(J_c\otimes M_c)(1\otimes\bar y)(J_c\otimes M_c)^*=J_cJ_c^*\otimes M_c\bar yM_c^*=1\otimes M_c\bar yM_c^*`$ である。
$`X=\sum_cV_c(1\otimes x_c)V_c^*`$ について $`c`$ にわたって和をとればよい。$`\square`$

セクター $`c`$ は、$`\sigma(c)=c`$ かつ $`\varepsilon_c=1`$ のとき **実型**、$`\sigma(c)=c`$ かつ
$`\varepsilon_c=-1`$ のとき **四元数型**、$`\sigma(c)\neq c`$ のとき **複素型** であるという；定理 4.20(c) により、これは $`\mathcal A`$
と $`P_c`$ の性質である。定理 4.20(b) により、$`\sigma(c)=c`$ であって $`k_c`$ または $`m_c`$ が奇数であるセクターは実型である。

*ライブラリでは：* $`\sigma`$、$`J_c`$、$`M_c`$、$`\varepsilon_c`$ とセクターの型は `Decomposition.real_structure()` で計算する — [チュートリアル、4.7節](../../tutorial/ja/04_symmetries.md#s-4-7)。

<a id="r-4-21"></a>

**定理 4.21（実代数）.** $`\Gamma`$ を、$`\sigma(c)=c`$ を満たすすべてのセクター $`c`$ と、$`\sigma(c)\neq c`$ を満たす各対
$`\{c,\sigma(c)\}`$ から一つずつ選んだセクターを含む集合とする。$`\sigma(c)=c`$ に対し $`R_c=\{x\in M_{k_c}(\mathbb C):xJ_c=J_c\bar x\}`$ とおく。このとき
$`a\mapsto(\pi_c(a))_{c\in\Gamma}`$ は実代数の同型であって、転置を随伴に写し、
```math
\mathcal A_{\mathbb R}\;\cong\bigoplus_{c\in\Gamma,\ \sigma(c)=c}R_c\;\oplus\bigoplus_{c\in\Gamma,\ \sigma(c)\neq c}M_{k_c}(\mathbb C),
```
となる。ここで最後の直和成分は実代数とみなす（**複素型**）。さらに：
(a)（*実型*）$`\varepsilon_c=1`$ ならば、$`J_c\bar W_J=W_J`$、$`M_c\bar W_M=W_M`$ により $`J_c=W_JW_J^T`$、$`M_c=W_MW_M^T`$ と書く
（補題 4.19(c)）。このとき $`x\mapsto W_J^*xW_J`$ は $`R_c`$ を $`M_{k_c}(\mathbb R)`$ の上へ同型に写す。
$`O_c=V_c(W_J\otimes W_M)`$ の列は実かつ正規直交であり、$`\mathbb C`$ 上で $`P_c`$ の値域を張り、$`a\in\mathcal A_{\mathbb R}`$ に対し実行列
$`\alpha_c(a)=W_J^*\pi_c(a)W_J`$ により $`aO_c=O_c(\alpha_c(a)\otimes1_{m_c})`$ である。
(b)（*四元数型*）$`\varepsilon_c=-1`$ ならば、$`J_c\bar W_J=W_J\Omega_{k_c}`$ により $`J_c=W_J\Omega_{k_c}W_J^T`$ と書く
（補題 4.19(c)）。このとき $`x\mapsto\mathbb M^{-1}(W_J^*xW_J)`$ は実代数の同型 $`R_c\cong M_{k_c/2}(\mathbb H)`$ であり、
随伴を四元数的随伴 $`(q^*)_{ij}=\bar q_{ji}`$ に写す；ここで $`\mathbb M:M_r(\mathbb H)\to M_{2r}(\mathbb C)`$ は
各成分 $`q`$ を補題 E.15(2) の $`2\times2`$ ブロック $`M(q)`$ で置き換える写像であり、その像は
$`\{y\in M_{2r}(\mathbb C):y\Omega_{2r}=\Omega_{2r}\bar y\}`$ である。
特に $`\dim_{\mathbb R}\mathcal A_{\mathbb R}=\sum_ck_c^2=\dim_{\mathbb C}\mathcal A`$ である。

*証明.* 実な $`a`$ に対しては $`\bar a=a`$ であり、定理 4.20(a) は $`\pi_{\sigma(c)}(a)=J_c\overline{\pi_c(a)}J_c^*`$ となる。
$`\sigma(c)=c`$ に対しては、これは $`\pi_c(a)\in R_c`$ を意味する。よって $`\Phi(a)=(\pi_c(a))_{c\in\Gamma}`$ は $`\mathcal A_{\mathbb R}`$ を右辺の中へ写す；
これは実代数の準同型であり、$`\Phi(a^T)=\Phi(a^*)=\Phi(a)^*`$ である。$`\Phi(a)=0`$ ならば、同じ公式により $`c\in\Gamma`$ に対し
$`\pi_{\sigma(c)}(a)=0`$ となるので、すべての $`\pi_e(a)=0`$ であり $`a=0`$ である。*全射性:* $`c\in\Gamma`$ とし、$`\sigma(c)=c`$ ならば $`x\in R_c`$、
そうでなければ $`x\in M_{k_c}`$ とし、$`a=V_c(x\otimes1)V_c^*\in\mathcal A`$ とおく。セクター
$`\sigma(e)`$ に対する定理 4.20(a) により $`\pi_e(\bar a)=J_{\sigma(e)}\overline{\pi_{\sigma(e)}(a)}J_{\sigma(e)}^*`$ であり、これは $`\sigma(e)=c`$ でない限り $`0`$ である。
$`\sigma(c)=c`$ ならば $`\pi_c(\bar a)=J_c\bar xJ_c^*=x`$（$`R_c`$ の定義）であるから、実行列 $`\mathrm{Re}\,a=\frac12(a+\bar a)\in\mathcal A_{\mathbb R}`$
の $`\Phi(\mathrm{Re}\,a)`$ は位置 $`c`$ で $`x`$、それ以外で $`0`$ である。$`\sigma(c)\neq c`$ ならば $`\pi_c(\bar a)=0`$ であり、$`2\mathrm{Re}\,a`$ は
$`\pi_c=x`$ と、$`e\notin\{c,\sigma(c)\}`$ に対し $`\pi_e=0`$ をもつ；$`\sigma(c)\notin\Gamma`$ であるから、$`\Phi(2\mathrm{Re}\,a)`$ は位置
$`c`$ で $`x`$、それ以外で $`0`$ である。これらの元の和が右辺全体を与える。

(a) $`xW_JW_J^T=W_JW_J^T\bar x`$ に左から $`W_J^*`$、右から $`\bar W_J`$ を掛ける；
$`W_J^T\bar W_J=\overline{W_J^*W_J}=1`$ であるから、これは $`W_J^*xW_J=W_J^T\bar x\bar W_J=\overline{W_J^*xW_J}`$ となり、この手順は
逆にたどれる。よって $`R_c=W_JM_{k_c}(\mathbb R)W_J^*`$ であり、$`x\mapsto W_J^*xW_J`$ は随伴と両立する実代数の同型である。次に
$`\overline{O_c}=\overline{V_c}(\bar W_J\otimes\bar W_M)=V_c(J_c\bar W_J\otimes M_c\bar W_M)=V_c(W_J\otimes W_M)=O_c`$、$`O_c^*O_c=1`$ であり、
$`aO_c=V_c(\pi_c(a)W_J\otimes W_M)=O_c(W_J^*\pi_c(a)W_J\otimes1)`$ であって、$`a\in\mathcal A_{\mathbb R}`$ に対し $`W_J^*\pi_c(a)W_J`$ は実である。

(b) $`J_c=W_J\Omega W_J^T`$（$`\Omega=\Omega_{k_c}`$）について同じ計算をすると、$`xJ_c=J_c\bar x`$ であることは $`y=W_J^*xW_J`$ が
$`y\Omega=\Omega\bar y`$ を満たすことと同値である。$`y`$ を $`2\times2`$ ブロック $`y_{ij}`$ に分ける；$`\Omega`$ はブロック対角であるから、$`y\Omega=\Omega\bar y`$ は
すべての $`i,j`$ に対する $`y_{ij}\Omega_2=\Omega_2\bar y_{ij}`$ と同値であり、$`y_{ij}=\begin{pmatrix}\alpha&\beta\\\gamma&\delta\end{pmatrix}`$ に対してこれは
$`\begin{pmatrix}\beta&-\alpha\\\delta&-\gamma\end{pmatrix}=\begin{pmatrix}-\bar\gamma&-\bar\delta\\\bar\alpha&\bar\beta\end{pmatrix}`$、すなわち
$`\delta=\bar\alpha`$ かつ $`\beta=-\bar\gamma`$ である：$`M(a,b)=\begin{pmatrix}a&-b\\\bar b&\bar a\end{pmatrix}`$ として
$`y_{ij}=\begin{pmatrix}\alpha&-\bar\gamma\\\gamma&\bar\alpha\end{pmatrix}=M(\alpha,\bar\gamma)`$ である。よって $`\mathbb M`$ の像は主張どおりである。補題 E.15(2) により、$`M`$ は
$`M(\bar q)=M(q)^*`$ を満たす単射な実線形準同型 $`\mathbb H\to M_2(\mathbb C)`$ である；行列の積は成分の和と積しか用いないから、
$`\mathbb M`$ は $`\mathbb M(q^*)=\mathbb M(q)^*`$ を満たす単射な実代数の準同型である（定理 E.24(1) の証明と同様）。したがって $`x\mapsto\mathbb M^{-1}(W_J^*xW_J)`$ が
主張の同型である。*次元:* どちらの場合も $`\dim_{\mathbb R}R_c=k_c^2`$ であり（$`2r=k_c`$ として $`\dim_{\mathbb R}M_r(\mathbb H)=4r^2`$）、
対は $`\dim_{\mathbb R}M_{k_c}(\mathbb C)=2k_c^2=k_c^2+k_{\sigma(c)}^2`$ を寄与する。$`\square`$

可換子環についても、因子を入れ替えて同じことが成り立つ：$`(\mathcal A')_{\mathbb R}`$ は同じセクターにわたる、同じ型の
$`M_{m_c}(\mathbb R)`$、$`M_{m_c/2}(\mathbb H)`$、$`M_{m_c}(\mathbb C)`$ の直和である。実際、
$`\mathcal A'`$ は $`\overline{\mathcal A'}=\mathcal A'`$ を満たす $`*`$-環であり（補題 4.18(b)）、上の証明は三つの事実しか用いておらず、
それらは $`\pi_c`$ の代わりに $`X\mapsto x_c`$、$`J_c`$ の代わりに $`M_c`$ として $`\mathcal A'`$ について成り立つ：$`X\mapsto x_c`$ は
$`M_{m_c}`$ の上への単位的 $`*`$-準同型である（$`V_c(1\otimes y)V_c^*\in\mathcal A'`$ はセクター $`c`$ でブロック $`y`$、それ以外で $`0`$ をもつ）、
$`(\bar X)_{\sigma(c)}=M_c\bar x_cM_c^*`$（定理 4.20(d)）、そして $`\sigma(c)=c`$ に対し $`M_c\bar M_c=\varepsilon_c1`$
（定理 4.20(b)）。

*ライブラリでは：* 定理 4.21 に適合した実直交基底 `Decomposition.real_basis()` — [チュートリアル、4.7節](../../tutorial/ja/04_symmetries.md#s-4-7)。

三つの例を挙げる。(i) $`\mathbb R^3`$ の回転の生成元 $`(L_a)_{bc}=-\epsilon_{abc}`$ が
$`\mathbb R^3\otimes\mathbb R^2`$ 上に $`L_a\otimes1_2`$ として作用する場合：セクターは一つで $`k=3`$、$`m=2`$ である（$`iL_x,iL_y,iL_z`$ は $`\mathbb C^3`$ 上のスピン表現であり、そこで $`iL_z`$ は単純固有値
$`1,0,-1`$ をもつので、系 4.12 によりスピン $`1`$ である）；$`k`$ が奇数であるから
これは実型であり、$`\mathcal A_{\mathbb R}\cong M_3(\mathbb R)`$ である。(ii) $`\mathbb R^2`$ の回転の生成元
$`\Omega_2`$ が $`\mathbb R^2\otimes\mathbb R^2`$ 上に $`\Omega_2\otimes1_2`$ として作用する場合：$`\Omega_2^2=-1`$ かつ
$`\Omega_2^*=-\Omega_2`$ であるから、環は $`\mathrm{span}\{1,\Omega_2\otimes1\}`$ である；$`\pm i`$ に対する $`\Omega_2`$ の固有ベクトル $`(1,\mp i)/\sqrt2`$ は
互いに共役であり、$`k=1`$、$`m=2`$ の二つのセクターを与え、それらは共役で入れ替わる：複素型の対であり、
$`\mathcal A_{\mathbb R}=\{\alpha+\beta\,\Omega_2\otimes1\}\cong\mathbb C`$ である。
(iii) 四元数 $`q=(a,b)\in\mathbb H`$（E.6節）による $`\mathbb H=\mathbb R^4`$ 上の左乗法 $`L_q`$。座標は
$`(\mathrm{Re}\,a,\mathrm{Im}\,a,\mathrm{Re}\,b,\mathrm{Im}\,b)`$ とする。$`q\mapsto L_q`$ は単射な実代数の準同型であり
（$`\mathbb H`$ は結合的、補題 E.15(2)）、$`L_{\bar q}=L_q^T`$ である。なぜなら、
$`\mathbb R^4`$ のユークリッドノルムであるノルム $`N(a,b)=|a|^2+|b|^2`$ に関して $`\langle qx,y\rangle=\langle x,\bar qy\rangle`$（補題 E.14(2)）だからである。したがって $`L_I,L_J,L_K`$（$`I=(i,0)`$、$`J=(0,1)`$、$`K=(0,i)`$）とそれらの転置が生成する実代数は
$`\{L_q:q\in\mathbb H\}\cong\mathbb H`$ である；複素の環 $`\mathcal A`$ は生成元の積の複素線形包
$`\{L_q\}+i\{L_q\}`$ であるから、$`\mathcal A_{\mathbb R}=\{L_q\}\cong\mathbb H`$ である（補題 4.18(a)）。これだけで分解が決まる：
$`\sum_ck_c^2=\dim_{\mathbb C}\mathcal A=\dim_{\mathbb R}\mathcal A_{\mathbb R}=4`$（定理 4.21）であり、すべての $`k_c=1`$ ではない。もしそうなら
$`\mathcal A`$、したがって $`\mathcal A_{\mathbb R}\cong\mathbb H`$ が可換になってしまうからである；よってセクターは一つで、$`k=2`$ かつ
$`\sigma(c)=c`$ であり、これは実型ではない。もし実型なら定理 4.21(a) により $`\mathcal A_{\mathbb R}\cong M_2(\mathbb R)`$ となり、これは
零因子をもつが、$`\mathbb H`$ は零因子をもたない（補題 E.15(3)）からである：これは四元数型であり、$`\mathbb R^4`$ 上で $`m=4/2=2`$ である。（ここで単位四元数は $`SU(2)`$ として作用する：$`M(q)`$ は行列式
$`|a|^2+|b|^2=1`$ のユニタリであり、すべての $`U=\begin{pmatrix}\alpha&\beta\\\gamma&\delta\end{pmatrix}\in SU(2)`$ はこの形である。
実際 $`U^*=U^{-1}=\begin{pmatrix}\delta&-\beta\\-\gamma&\alpha\end{pmatrix}`$ から $`\delta=\bar\alpha`$、$`\beta=-\bar\gamma`$ が得られ、$`U=M(\alpha,\bar\gamma)`$ となる。）
$`q`$ が $`1_2\otimes L_q`$ として作用する $`\mathbb H^2=\mathbb R^8`$ 上では、同じ議論により $`k=2`$、$`m=4`$ の
四元数型のセクターが一つ得られ、$`\mathcal A_{\mathbb R}\cong M_1(\mathbb H)=\mathbb H`$ である。

*計算：* 三つの例を `Decomposition.real_structure()` と `Decomposition.real_basis()` で — [チュートリアル、4.7節](../../tutorial/ja/04_symmetries.md#s-4-7)。

実行列による群の表現論において、次はシューアの補題の実形である；これは三つの型が実表現そのものから読み取れることを示す。

<a id="r-4-22"></a>

**定理 4.22（実シューア補題）.** $`U\subseteq\mathbb R^n`$ を、すべての
$`a\in\mathcal A_{\mathbb R}`$ に対し $`aU\subseteq U`$ を満たす零でない部分空間で、**既約**、すなわち $`0`$ と $`U`$ がこの性質をもつ唯一の部分空間であるものとする。$`D_U`$ を、
すべての $`a\in\mathcal A_{\mathbb R}`$、$`u\in U`$ に対し $`Tau=aTu`$ を満たす実線形写像 $`T:U\to U`$ 全体とする。
(a) $`D_U`$ は単位元をもち零因子をもたない結合的実代数であり、したがって $`\mathbb R`$、$`\mathbb C`$
または $`\mathbb H`$ と同型である（定理 E.19）。
(b) $`\sigma`$ のちょうど一つの軌道 $`\mathcal O=\{c,\sigma(c)\}`$ について、$`U`$ は $`P_{\mathcal O}=\sum_{e\in\mathcal O}P_e`$ の値域に含まれ、
$`c`$ が実型、複素型、四元数型であるのに応じて $`\dim_{\mathbb R}D_U=1,2,4`$、$`\dim_{\mathbb R}U=k_c,2k_c,2k_c`$ である。よってそれぞれ
$`D_U\cong\mathbb R,\mathbb C,\mathbb H`$ である。
(c) 任意の軌道 $`\mathcal O`$ に対し、$`P_{\mathcal O}`$ の値域の中にそのような $`U`$ が存在する。$`\mathbb R^n`$ 自身が既約であるのは、
$`\sigma`$ がただ一つの軌道 $`\{c,\sigma(c)\}`$ をもち、$`m_c=1`$（実型と複素型）、それぞれ $`m_c=2`$（四元数型）であるときかつそのときに限る。

*証明.* (a) $`T\in D_U`$ に対し $`\ker T`$ と $`TU`$ は不変部分空間であるから、$`T=0`$ であるか $`T`$ は全単射であり、後者のとき
$`T^{-1}\in D_U`$ である。$`D_U`$ は合成について閉じ、合成は結合的であり、$`\mathrm{id}_U\neq0`$ を含み、$`ST=0`$ かつ
$`S\neq0`$ ならば $`T=S^{-1}ST=0`$ である。

(b) *不変部分空間と実射影.* $`U'\subseteq\mathbb R^n`$ が不変ならば、$`\mathbb R^n`$ におけるその直交補空間も不変である
（$`y\perp U'`$、$`u\in U'`$、$`a\in\mathcal A_{\mathbb R}`$ に対し、$`a^T\in\mathcal A_{\mathbb R}`$ であるから $`u\cdot ay=a^Tu\cdot y=0`$）。
したがって $`U'`$ への直交射影 $`e`$（実対称な冪等行列）は
$`\mathcal A_{\mathbb R}`$ と可換である；補題 4.18(b) により $`e\in(\mathcal A')_{\mathbb R}`$ である。逆に、任意の
射影 $`e\in(\mathcal A')_{\mathbb R}`$ に対し $`e\mathbb R^n`$ は不変である。よって $`U'\mapsto e`$ は不変部分空間全体から
$`(\mathcal A')_{\mathbb R}`$ の射影全体の集合 $`\mathcal P`$ の上への包含関係を保つ全単射であり、$`U`$ が既約であることはその射影 $`e`$
が $`\mathcal P`$ の極小な零でない元であることと同値である。値域 $`e\mathbb C^n`$（その実ベクトルは
$`e\mathbb R^n`$）に補題 4.18(a) を適用すると $`\dim_{\mathbb R}U=\mathrm{rank}\,e`$ である。

*$`D_U\cong e(\mathcal A')_{\mathbb R}e`$.* $`T\in D_U`$ に対し、$`\hat T`$ を実 $`n\times n`$ 行列 $`x\mapsto T(ex)`$ とする。
$`a\in\mathcal A_{\mathbb R}`$ に対し $`\hat Tax=T(eax)=T(aex)=aT(ex)=a\hat Tx`$ であるから、$`\hat T\in(\mathcal A')_{\mathbb R}`$（補題 4.18(b)）かつ
$`\hat T=e\hat Te`$ である。逆に $`X\in(\mathcal A')_{\mathbb R}`$ に対し、$`eXe`$ を $`U`$ に制限したものは $`D_U`$ に属する。これらの写像は
互いに逆であり、乗法的である。さらに $`e\mathcal A'e`$ は共役で閉じており、その実な元は
$`e(\mathcal A')_{\mathbb R}e`$ の元である（$`eXe`$ が実ならば、$`\mathrm{Re}\,X\in(\mathcal A')_{\mathbb R}`$ により $`eXe=e\mathrm{Re}(X)e`$）から、
補題 4.18(a) により $`\dim_{\mathbb R}D_U=\dim_{\mathbb C}e\mathcal A'e`$ である。

*分解における射影.* $`\mathcal A'`$ の射影 $`e`$ は射影
$`f_c\in M_{m_c}`$ により $`e=\sum_cV_c(1\otimes f_c)V_c^*`$ と書ける（定理 2.6）。このとき $`\mathrm{rank}\,e=\sum_ck_c\mathrm{rank}\,f_c`$ であり、
$`e\mathcal A'e=\sum_cV_c(1\otimes f_cM_{m_c}f_c)V_c^*`$ の複素次元は $`\sum_c(\mathrm{rank}\,f_c)^2`$ であり、
定理 4.20(d) により、$`\bar e=e`$ はすべての $`c`$ に対する $`f_{\sigma(c)}=M_c\bar f_cM_c^*`$ と同値である。

*極小元.* $`e`$ を $`\mathcal P`$ において極小とする。$`\sigma`$ の各軌道 $`\mathcal O`$ に対し、$`P_{\mathcal O}`$ は
$`Z(\mathcal A)\subseteq\mathcal A'`$ の実な中心射影であり（定理 4.20(a) により $`\bar P_{\mathcal O}=P_{\mathcal O}`$）、
$`e`$ と可換であるから、$`eP_{\mathcal O}\in\mathcal P`$ かつ $`eP_{\mathcal O}\le e`$ である；$`e=\sum_{\mathcal O}eP_{\mathcal O}`$ であるから、
ちょうど一つの軌道 $`\mathcal O`$ について $`eP_{\mathcal O}=e`$ である。すなわち $`U=e\mathbb R^n`$ は $`P_{\mathcal O}`$ の値域に含まれ、
$`c\notin\mathcal O`$ に対し $`f_c=0`$ である。$`\mathcal O=\{c,d\}`$、
$`d\neq c`$ ならば、$`f_d=M_c\bar f_cM_c^*`$ は $`f_c`$ と同じ階数をもつ；階数一の射影 $`g\le f_c`$ に対し、
$`Y=V_c(1\otimes g)V_c^*`$ として $`g'=Y+\bar Y`$ は $`\mathcal A'`$ の二つの直交する射影（セクター $`c`$ と $`d`$ における）の和であり、実であって、
$`g'\le e`$ である（$`Y\le e`$ かつ $`\bar Y\le\bar e=e`$ であるから）；よって $`g'=e`$、$`\mathrm{rank}\,f_c=\mathrm{rank}\,f_d=1`$、$`\dim_{\mathbb C}e\mathcal A'e=2`$、$`\mathrm{rank}\,e=2k_c`$ である。
$`\mathcal O=\{c\}`$ ならば、$`f=f_c`$、$`M=M_c`$、$`\theta v=M\bar v`$ とおく。これは $`\mathbb C^{m_c}`$ 上の反ユニタリ写像であり、
$`\theta^2=M\bar M=\varepsilon_c`$ である（定理 4.20(b)）。$`fM=M\bar f`$ から $`\theta(fv)=M\bar f\bar v=fM\bar v=f\theta v`$ であるので、
$`f`$ の値域 $`F`$ は $`\theta`$ 不変である。$`F_1\subseteq F`$ が射影 $`g`$ をもつ $`\theta`$ 不変部分空間ならば、$`M\bar gM^*`$ は
$`M\bar F_1=\theta F_1=F_1`$ への射影である（$`\theta`$ は単射）から、$`g=M\bar gM^*`$ であり、$`V_c(1\otimes g)V_c^*\in\mathcal P`$
は $`e`$ の下にある。極小性により、$`F`$ は $`0`$ と $`F`$ 以外に $`\theta`$ 不変部分空間をもたない。$`\varepsilon_c=1`$ ならば、$`F`$ における $`\theta`$ の
固定単位ベクトル（補題 4.19(a)）が一つを張るので、$`\mathrm{rank}\,f=1`$、$`\dim_{\mathbb C}e\mathcal A'e=1`$、
$`\mathrm{rank}\,e=k_c`$ である。$`\varepsilon_c=-1`$ ならば、単位ベクトル $`u\in F`$ の平面 $`\mathrm{span}\{u,\theta u\}`$ が一つであるので
（補題 4.19(b)）、$`\mathrm{rank}\,f=2`$、$`\dim_{\mathbb C}e\mathcal A'e=4`$、$`\mathrm{rank}\,e=2k_c`$ である。$`\mathbb R`$、
$`\mathbb C`$、$`\mathbb H`$ の次元は $`1,2,4`$ であるから、(a) によりそれぞれ $`D_U\cong\mathbb R,\mathbb C,\mathbb H`$ である。

(c) $`P_{\mathcal O}\in\mathcal P`$ は零でなく、その下にある $`\mathcal P\setminus\{0\}`$ の元のうち階数が最小のものは極小である。$`\mathbb R^n`$
が既約であることは $`e=1`$ が極小であることと同値である；そのブロックは $`f_c=1_{m_c}`$ であるから、(b) の証明により、これには軌道がただ一つであること、
および三つの場合にそれぞれ $`m_c=1,1,2`$ であることが必要である。逆に、これらの条件の下では、零でない $`e'\in\mathcal P`$ は
その軌道のセクターについて $`f'_c\neq0`$ をもち（$`f'_{\sigma(c)}=M_c\bar f'_cM_c^*`$）、したがって $`m_c=1`$ ならば $`f'_c=1`$ である；四元数型の
場合には $`f'_c`$ の値域は $`\theta`$ 不変であり、したがって偶数次元である（補題 4.19(b)）から、$`f'_c=1_2`$ である。よって $`e'=1`$ である。$`\square`$

実直交行列で作用する群 $`G`$ に対し、$`\mathcal A=\mathrm{span}\,\pi(G)`$ は共役で閉じており
（補題 4.18(c)）、$`\mathcal A_{\mathbb R}`$ は $`\pi(G)`$ の実線形包であり（補題 4.18(a)：$`\mathcal A_{\mathbb R}`$ は $`\mathcal A`$ の元の
実部からなる）、$`D_U`$ は $`G`$ の既約な実表現の可換子環である：それは
$`\mathbb R`$、$`\mathbb C`$ または $`\mathbb H`$ である。(b) の証明により、複素化 $`U+iU=e\mathbb C^n=V_c(\mathbb C^{k_c}\otimes F)`$
（対の場合はその共役を加えたもの）は、三つの場合にそれぞれ $`\pi_c`$、$`\pi_c\oplus\pi_{\sigma(c)}`$、$`\pi_c\oplus\pi_c`$ を担う：すなわち
既約であるか、二つの非同値な既約表現の和であるか、二つの同値な既約表現の和である。

*計算：* 三つの例の既約な実加群 $`\mathbb R^3`$、$`\mathbb R^2`$、$`\mathbb R^4`$ に対する定理 4.22 — [チュートリアル、4.7節](../../tutorial/ja/04_symmetries.md#s-4-7)。

<a id="r-4-23"></a>

**命題 4.23（可換子環の実な元）.** $`X\in\mathcal A'`$ を実とし、そのブロックを $`x_c`$ とする。このとき
すべての $`c`$ に対し $`x_{\sigma(c)}=M_c\bar x_cM_c^*`$ である。さらに $`X=X^T`$ ならば：
(a) 実型の $`c`$ に対し、$`y_c=W_M^*x_cW_M`$（$`W_M`$ は定理 4.21(a) のとおり）は実対称行列であり、
定理 4.21(a) の実で正規直交な列 $`O_c`$ について $`XO_c=O_c(1_{k_c}\otimes y_c)`$ である；
(b) 四元数型の $`c`$ に対し、$`x_c`$ のすべての固有値は偶数の重複度をもち、したがって $`P_c`$ の値域上の $`X`$ のすべての固有値は
$`2k_c`$ で割り切れる重複度をもつ（**クラマース縮退**）；
(c) 複素型の $`c`$ に対し、$`x_c`$ と $`x_{\sigma(c)}`$ は同じ固有値を同じ重複度でもつ。
$`X=X^T`$ がなくても、(a) の $`y_c`$ は実である。

*証明.* $`X=\bar X`$ であるから、定理 4.20(d) により $`x_{\sigma(c)}=M_c\bar x_cM_c^*`$ である。$`X=X^T`$ ならば $`X^*=X`$ であり、
各 $`x_c`$ はエルミートである。(a) $`M_c=W_MW_M^T`$、$`M_c^*=\bar W_MW_M^*`$、$`W_M^T\bar W_M=1`$ により
$`\bar y_c=W_M^T\bar x_c\bar W_M=W_M^TM_c^*x_cM_c\bar W_M=W_M^*x_cW_M=y_c`$ であり、$`y_c`$ はエルミートであるから実対称である。さらに
$`XO_c=XV_c(W_J\otimes W_M)=V_c(W_J\otimes x_cW_M)=O_c(1\otimes y_c)`$ である。(b) $`\theta v=M_c\bar v`$ とおく。これは
$`\theta^2=-1`$ を満たす反ユニタリ写像である（定理 4.20(b)）。$`x_cM_c=M_c\bar x_c`$ から $`x_c\theta v=M_c\bar x_c\bar v=\theta(x_cv)`$ であるので、$`x_c`$ の固有値
$`\lambda\in\mathbb R`$ に対し固有空間 $`E_\lambda`$ は $`\theta E_\lambda\subseteq E_\lambda`$ を満たし（$`x_c\theta v=\theta(\lambda v)=\lambda\theta v`$）、
補題 4.19(b) により $`\dim E_\lambda`$ は偶数である。$`P_c`$ の値域上で $`X`$ は $`V_c(1_{k_c}\otimes x_c)V_c^*`$ として作用し、その固有値は
$`x_c`$ の固有値をそれぞれ $`k_c`$ 回繰り返したものである。(c) $`\bar x_c`$ は $`x_c`$ の固有値 $`\bar\lambda=\lambda`$ を、
共役な固有ベクトルとともにもち、$`x_{\sigma(c)}`$ は $`\bar x_c`$ とユニタリ同値である。$`\square`$

四元数型の場合、補題 4.19(c) のように $`M_c=W_M\Omega_{m_c}W_M^T`$ とすると、$`y_c=W_M^*x_cW_M`$ は
$`y_c\Omega=\Omega\bar y_c`$ を満たす（$`J_c`$ の代わりに $`M_c`$ とした定理 4.21(b) の計算）ので、ある四元数行列 $`h`$ により $`y_c=\mathbb M(h)`$ である；
$`X=X^T`$ ならば、$`\mathbb M(h^*)=y_c^*=y_c`$ かつ $`\mathbb M`$ は単射であるから $`h^*=h`$ である。

*計算：* 三つの例のそれぞれにおける可換子環のランダムな実対称元に対する命題 4.23 — [チュートリアル、4.7節](../../tutorial/ja/04_symmetries.md#s-4-7)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. ハミルトニアンが対称性と可換であり、
> 対称性が実となる基底において実であるならば——これは時間反転対称性の効果であり、そこでは時間反転は
> 複素共役として作用する——、命題 4.23 はそのブロックがどのような形かを述べる：実型のセクターでは実対称、複素型の対では同じスペクトルをもつ相手を伴う複素
> エルミート、四元数型のセクターでは二重に縮退した準位をもつ四元数エルミートである。最後のものは、平方が
> $`-1`$ である奇数個のスピン $`\frac12`$ の時間反転の場合である（クラマースの定理）。ダイソンの *三重の道*（threefold way；F. J. Dyson, *J. Math. Phys.* 3 (1962) 1199）は、各種類の
> 一般的なハミルトニアンのブロックを対応するクラスのランダム行列——5.3節の GOE と GUE、
> および本講義では定義しないガウス型シンプレクティック集団（GSE）——でモデル化し、それに応じて各
> ブロックの準位統計を予言する。

<a id="s-4-8"></a>

## 4.8 比較の対象

*ライブラリでは：* 他のソフトウェアとの比較 — [チュートリアル、4.8節](../../tutorial/ja/04_symmetries.md#s-4-8)。

<a id="s-4-9"></a>

## 4.9 数値分解の安定性

4.3節では厳密な演算のもとでアルゴリズムを解析した。実際には生成元は近似的にしかわからない（測定値であったり丸められていたりする）。
そして対称性は許容誤差によって判定される：可換子環は、交換子写像の特異値のうち閾値 $`\theta`$ を下回るものに対応する特異ベクトルで
張られる。本節では、どの閾値が正しいか、計算された対象が厳密な対象からどれだけ離れているか、そして出力だけから何が保証できるかを述べる。

本節を通じて $`H_1,\dots,H_J`$ はエルミートであり、$`\mathcal A=\{H_j\}''`$ とする（ライブラリはまず生成元をそのエルミート部分に置き換える。
それらは同じ環を生成する、4.3節）。**交換子写像** は
```math
L:M_n(\mathbb C)\to M_n(\mathbb C)^J,\qquad LX=([H_1,X],\dots,[H_J,X])
```
であり、$`M_n`$ にはノルム $`\|X\|_2`$、$`M_n^J`$ には $`\|(X_j)_j\|_2=(\sum_j\|X_j\|_2^2)^{1/2}`$ を入れる。$`\ker L=\mathcal A'`$ の次元は
$`d=\sum_cm_c^2`$ である。有限次元複素内積空間の間の線形写像 $`T`$（定義域の次元を $`p`$ とする）に対し、$`\sigma_1(T)\le\dots\le\sigma_p(T)`$ で
$`T^*T\ge0`$ の固有値の平方根を *昇順* に並べたものを表す（正方行列ではこれは1.4節の特異値を逆順に並べたものである）。
$`|Tx|^2=\langle x,T^*Tx\rangle`$ だから $`\ker T=\ker T^*T`$ かつ $`\sigma_p(T)=\|T\|`$ である。交換子写像の **ギャップ** とは
$`\gamma=\sigma_{d+1}(L)`$、すなわちその最小の非零特異値である（$`d=n^2`$ なら $`\gamma=\infty`$）。$`\|AX\|_2\le\|A\|\|X\|_2`$（$`X`$ の各列に
$`|Ax|\le\|A\||x|`$ を適用する）と $`\|XA\|_2=\|A^*X^*\|_2\le\|X\|_2\|A\|`$ を用いる。したがって $`\|[E,X]\|_2\le2\|E\|\|X\|_2`$ である。

<a id="r-4-24"></a>

**補題 4.24（摂動）.** (a) $`T,\tilde T`$ を同じ有限次元複素内積空間の間の線形写像とし、$`\Delta=\|\tilde T-T\|`$ とする。このとき
すべての $`i`$ について $`|\sigma_i(\tilde T)-\sigma_i(T)|\le\Delta`$ である。さらに $`\sigma_d(T)=0<\gamma=\sigma_{d+1}(T)`$ かつ $`\Delta<\gamma/2`$ なら、
$`\Pi`$ を $`\ker T`$ への直交射影、$`\tilde\Pi`$ を $`\tilde T^*\tilde T`$ の固有値 $`\le\Delta^2`$ に対するスペクトル射影とすると、$`\tilde\Pi`$ の
階数は $`d`$ であり、$`\|\tilde\Pi-\Pi\|\le2\Delta/\gamma`$ である。

(b) $`A,\tilde A\in\mathrm{Herm}_n`$、$`\|\tilde A-A\|\le\Delta`$ とし、$`\mu`$ を $`A`$ の重複度 $`r`$ の固有値、$`P`$ をそのスペクトル射影とし、
$`\mu`$ は $`A`$ の他の固有値から距離 $`g>2\Delta`$ 離れているとする。このとき $`\tilde A`$ は $`[\mu-\Delta,\mu+\Delta]`$ にちょうど $`r`$ 個の
固有値をもち、他の固有値はすべて $`\mu`$ から距離 $`\ge g-\Delta`$ にあり、これら $`r`$ 個の固有値に対する $`\tilde A`$ のスペクトル射影 $`\tilde P`$ は
$`\|\tilde P-P\|\le2\Delta/g`$ を満たす。

(c) $`P,\tilde P`$ を階数 $`k`$ の直交射影で $`\|\tilde P-P\|\le\varepsilon\le\frac12`$ とし、$`\tilde Q\in M_{n\times k}`$ の列が
$`\mathrm{ran}\,\tilde P`$ の正規直交基底をなすとする。このとき $`S=\tilde Q^*P\tilde Q`$ は正定値であり、$`Q=P\tilde QS^{-1/2}`$ の列は
$`\mathrm{ran}\,P`$ の正規直交基底をなし、$`\|Q-\tilde Q\|\le2\varepsilon`$ である。

(d) 可逆な $`M\in M_k`$ に対し $`\mathrm{pol}(M)=M(M^*M)^{-1/2}`$ はユニタリである（定理 1.8 のユニタリ）。$`\lambda>0`$、$`U`$ ユニタリで
$`\|M-\lambda U\|\le\zeta<\lambda`$ ならば、$`M`$ は可逆で $`\|\mathrm{pol}(M)-U\|\le2\zeta/\lambda`$ である。

*証明.* (a) 定理 A.8 を $`T^*T`$ に適用し、$`t\mapsto t^{1/2}`$ が増加であることから
$`\sigma_i(T)=\min_{\dim U=i}\max_{x\in U,|x|=1}|Tx|`$ である。単位ベクトル $`x`$ に対し $`|\tilde Tx|\le|Tx|+\Delta`$；$`U`$ 上で最大をとり、次に最小を
とれば $`\sigma_i(\tilde T)\le\sigma_i(T)+\Delta`$ を得て、$`T,\tilde T`$ を入れ替えれば逆向きの不等式を得る。よって $`i\le d`$ では
$`\sigma_i(\tilde T)\le\Delta`$、$`i>d`$ では $`\sigma_i(\tilde T)\ge\gamma-\Delta>\Delta`$ であり、$`\tilde\Pi`$ の階数は $`d`$ である。単位ベクトル
$`x\in\mathrm{ran}\,\tilde\Pi`$ に対し、$`\tilde T^*\tilde T`$ の固有ベクトルで展開すると $`|\tilde Tx|^2=\langle x,\tilde T^*\tilde Tx\rangle\le\Delta^2`$ となるので
$`|Tx|\le2\Delta`$；一方 $`Tx=T(1-\Pi)x`$ で $`(1-\Pi)x\perp\ker T^*T`$ であり、そこでは $`T^*T`$ の固有値は $`\ge\gamma^2`$ だから
$`|Tx|\ge\gamma|(1-\Pi)x|`$。よって $`\|(1-\Pi)\tilde\Pi\|\le2\Delta/\gamma`$。単位ベクトル $`x\in\ker T`$ に対しては $`|\tilde Tx|\le\Delta`$ であり、
$`|\tilde Tx|^2=\langle\tilde\Pi x,\tilde T^*\tilde T\tilde\Pi x\rangle+\langle(1-\tilde\Pi)x,\tilde T^*\tilde T(1-\tilde\Pi)x\rangle\ge(\gamma-\Delta)^2|(1-\tilde\Pi)x|^2`$
（$`\tilde\Pi`$ は $`\tilde T^*\tilde T`$ と可換であり、$`\mathrm{ran}(1-\tilde\Pi)`$ 上でのその固有値は $`\ge(\gamma-\Delta)^2`$）だから
$`\|(1-\tilde\Pi)\Pi\|\le\Delta/(\gamma-\Delta)\le2\Delta/\gamma`$。最後に、
$`(\Pi-\tilde\Pi)x=\Pi(1-\tilde\Pi)(1-\tilde\Pi)x-(1-\Pi)\tilde\Pi\,\tilde\Pi x`$ は直交する二つのベクトルの和であり、
$`\|\Pi(1-\tilde\Pi)\|=\|(1-\tilde\Pi)\Pi\|`$（随伴、補題 1.1）だから
$`|(\Pi-\tilde\Pi)x|^2\le(2\Delta/\gamma)^2(|(1-\tilde\Pi)x|^2+|\tilde\Pi x|^2)=(2\Delta/\gamma)^2|x|^2`$。

(b) 固有値を昇順に番号付ける。定理 A.10 により、$`\lambda_i(A)=\mu`$ となる $`r`$ 個の固有値 $`\lambda_i(\tilde A)`$ は $`[\mu-\Delta,\mu+\Delta]`$
にあり、他は $`\mu`$ から距離 $`\ge g-\Delta>\Delta`$ にある。よって $`(\tilde A-\mu)^2`$ の小さい方から $`r`$ 個の固有値（いずれも $`\le\Delta^2`$）は
これら $`r`$ 個の固有値に対応し、$`T=A-\mu`$、$`\tilde T=\tilde A-\mu`$（このとき $`\ker T=\mathrm{ran}\,P`$、$`\sigma_{r+1}(T)=g`$）として (a) を
適用すれば評価を得る。

(c) $`(1-P)\tilde Q=(1-P)\tilde P\tilde Q=(1-P)(\tilde P-P)\tilde Q`$ のノルムは $`\le\varepsilon`$ だから、$`S=1-\tilde Q^*(1-P)\tilde Q`$ の固有値は
$`[1-\varepsilon^2,1]`$ にあり（定理 A.10）、$`S^{-1/2}`$ が存在する（定理 1.5）。$`Q^*Q=S^{-1/2}\tilde Q^*P\tilde QS^{-1/2}=1`$ であり、
$`\mathrm{ran}\,Q\subseteq\mathrm{ran}\,P`$ の次元は $`k`$ である。$`Q-\tilde Q=\tilde Q(S^{-1/2}-1)-(1-P)\tilde QS^{-1/2}`$ より
$`\|Q-\tilde Q\|\le(1-\varepsilon^2)^{-1/2}-1+\varepsilon(1-\varepsilon^2)^{-1/2}=f(\varepsilon)`$、ただし $`f(\varepsilon)=((1+\varepsilon)/(1-\varepsilon))^{1/2}-1`$。
その導関数 $`(1-\varepsilon)^{-1}(1-\varepsilon^2)^{-1/2}`$ は増加なので $`f`$ は凸で $`f(0)=0`$ であり、
$`f(\varepsilon)\le2\varepsilon f(\frac12)=2\varepsilon(\sqrt3-1)<2\varepsilon`$。

(d) $`M^*M>0`$ であり、$`\mathrm{pol}(M)^*\mathrm{pol}(M)=(M^*M)^{-1/2}M^*M(M^*M)^{-1/2}=1`$；$`\mathrm{pol}(M)|M|=M`$ だから、これは定理 1.8 の
ユニタリである。$`M=\lambda U(1+G)`$、$`\|G\|=\|M-\lambda U\|/\lambda\le\zeta/\lambda<1`$ と書く。(a)（$`T=1`$、$`\tilde T=1+G`$ とする）により
$`R=((1+G)^*(1+G))^{1/2}`$ の固有値は $`[1-\|G\|,1+\|G\|]`$ にあるので、$`R`$、$`1+G`$、$`M`$ は可逆である；$`M^*M=\lambda^2R^2`$ だから
$`\mathrm{pol}(M)=UW`$、ただし $`W=(1+G)R^{-1}`$ はユニタリであり、$`\|W-1\|\le\|W-WR\|+\|WR-1\|=\|1-R\|+\|G\|\le2\|G\|`$。$`\square`$

<a id="r-4-25"></a>

**定理 4.25（許容誤差のもとでの可換子環）.** $`\tilde H_1,\dots,\tilde H_J`$ をエルミートとし、$`\eta=(\sum_j\|\tilde H_j-H_j\|^2)^{1/2}`$、
$`\tilde L`$ をその交換子写像とする。$`\theta\ge0`$ に対し、$`\tilde L^*\tilde L`$ の固有値 $`\le\theta^2`$ の固有ベクトルが張る空間を $`\mathcal C_\theta`$ とし、
これを **許容誤差** $`\theta`$ **の可換子環** と呼ぶ（$`\mathcal C_0=\{\tilde H_j\}'`$）。

(a) すべての $`i`$ について $`|\sigma_i(\tilde L)-\sigma_i(L)|\le2\eta`$；特に $`i\le d`$ では $`\sigma_i(\tilde L)\le2\eta`$ であり、
$`\sigma_{d+1}(\tilde L)\ge\gamma-2\eta`$ である。

(b) $`2\eta<\theta<\gamma-2\eta`$ ならば、$`\dim\mathcal C_\theta=d`$、$`\mathcal C_\theta`$ は $`X\mapsto X^*`$ で閉じており、$`\mathcal C_\theta`$ への直交射影
$`\tilde\Pi`$ と $`\mathcal A'`$ への直交射影 $`\Pi`$ は $`\|\tilde\Pi-\Pi\|\le4\eta/\gamma`$ を満たす。

(c) $`\theta\ge\gamma+2\eta`$ ならば $`\dim\mathcal C_\theta>d`$ である。$`\theta<2\eta`$ では $`\dim\mathcal C_\theta<d`$ が起こりうる：
$`H_1=\mathrm{diag}(1,1,0)`$、$`\tilde H_1=H_1+\eta(E_{12}+E_{21})`$ に対し、$`d=5`$、$`\gamma=1`$ であり、$`\tilde L`$ の特異値は
$`0,0,0,2\eta,2\eta,1-\eta,1-\eta,1+\eta,1+\eta`$ で、$`\theta<2\eta`$ のとき $`\dim\mathcal C_\theta=3`$ である。

したがって許容誤差は雑音レベルの2倍を超え、かつギャップから雑音レベルの2倍を引いたものより小さくなければならない；このような
$`\theta`$ が存在するのは $`\eta<\gamma/4`$ のとき、かつそのときに限り、そのとき計算された可換子環は正しい次元をもち、真の可換子環から
$`4\eta/\gamma`$ 以内にある。評価 $`2\eta`$ は (c) の例で達成される。

*証明.* (a) $`\|(\tilde L-L)X\|_2^2=\sum_j\|[\tilde H_j-H_j,X]\|_2^2\le4\eta^2\|X\|_2^2`$ であり、$`i\le d`$ では $`\sigma_i(L)=0`$；補題 4.24(a) を
適用する。(b) (a) により $`i\le d`$ では $`\sigma_i(\tilde L)\le2\eta<\theta`$、$`i>d`$ では $`\sigma_i(\tilde L)\ge\gamma-2\eta>\theta`$ だから、
$`\mathcal C_\theta`$ は $`\Delta=2\eta<\gamma/2`$ とした補題 4.24(a) の射影 $`\tilde\Pi`$ の像である。エルミートな $`H`$ に対し
$`\langle Y,[H,X]\rangle_{HS}=\mathrm{Tr}(Y^*HX-Y^*XH)=\langle[H,Y],X\rangle_{HS}`$ だから $`\tilde L^*\tilde LX=\sum_j[\tilde H_j,[\tilde H_j,X]]`$ であり、
随伴をとると $`(\tilde L^*\tilde LX)^*=\sum_j[[X^*,\tilde H_j],\tilde H_j]=\tilde L^*\tilde L(X^*)`$。よって $`X\mapsto X^*`$ は $`\tilde L^*\tilde L`$ の各固有空間
（その固有値は実数）をそれ自身に写す。(c) (a) により $`\sigma_{d+1}(\tilde L)\le\gamma+2\eta\le\theta`$。例では
$`\mathbb C^3=\mathbb C^2\oplus\mathbb C`$ に対し $`X=\begin{pmatrix}X_{11}&x_{12}\\x_{21}&x_{22}\end{pmatrix}`$、$`\tilde H_1=(1+\eta\sigma)\oplus0`$、
$`\sigma=\begin{pmatrix}0&1\\1&0\end{pmatrix}`$ と書く。すると $`[\tilde H_1,X]`$ のブロックは $`\eta[\sigma,X_{11}]`$、$`(1+\eta\sigma)x_{12}`$、
$`-x_{21}(1+\eta\sigma)`$、$`0`$ であり、これは直交分解である。$`X_{11}\mapsto[\sigma,X_{11}]`$ は $`1,\sigma`$ 上で消え、正規直交な組
$`\mathrm{diag}(1,-1)/\sqrt2`$、$`\begin{pmatrix}0&-i\\i&0\end{pmatrix}/\sqrt2`$ を互いに相手の $`\mp2i`$ 倍に写す；$`1+\eta\sigma`$ の固有値は
$`1\pm\eta`$ である。これで特異値が得られる；$`\eta=0`$ とすれば $`d=5`$、$`\gamma=1`$ であり、$`\theta<2\eta`$ では三つの零だけが残る。$`\square`$

ライブラリでは $`\theta`$ は `tol` に $`\max(1,\sigma_{\max})`$ を掛けたものである（モジュールの docstring）；核は未知数を縮約したのちに計算され、
この縮約は厳密な生成元に対しては厳密である（4.3節）。報告される特異値は縮約された写像のものである。定理 4.25 は $`\eta`$ と $`\gamma`$ を
必要とするが、実際にはこれらは未知である。それでも出力は保証できる：

<a id="r-4-26"></a>

**命題 4.26（事後証明書）.** $`\tilde H_1,\dots,\tilde H_J`$ をエルミート、$`\tilde V`$ をユニタリ、$`\sum_ck_cm_c=n`$ とする；
$`\tilde{\mathcal B}=\tilde V(\bigoplus_cM_{k_c}\otimes1_{m_c})\tilde V^*`$、$`E`$ を $`\langle\cdot,\cdot\rangle_{HS}`$ に関する $`M_n`$ から $`\tilde{\mathcal B}`$ への
直交射影、$`H_j=E(\tilde H_j)`$、$`\tilde V^*H_j\tilde V=\bigoplus_ca_{j,c}\otimes1_{m_c}`$、$`r=(\sum_j\|\tilde H_j-H_j\|_2^2)^{1/2}`$（**後退誤差**）、
$`d=\sum_cm_c^2`$ とし、$`\gamma_{\mathcal B}`$ を $`H_1,\dots,H_J`$ の交換子写像の $`\sigma_{d+1}`$ とする。

(a) $`\gamma_{\mathcal B}`$ は、$`k_c\ge2`$ である $`c`$ についての $`\sigma_2(\Phi_{cc})`$ と、$`c\neq e`$ についての $`\sigma_1(\Phi_{ce})`$ の最小値である
（どちらもなければ $`\infty`$）。ここで $`\Phi_{ce}:M_{k_c\times k_e}\to M_{k_c\times k_e}^J`$、$`Y\mapsto(a_{j,c}Y-Ya_{j,e})_j`$。

(b) $`\gamma_{\mathcal B}>0`$ ならば $`H_j`$ は $`\tilde{\mathcal B}`$ を生成する：$`\tilde V`$ と $`(k_c,m_c)`$ は、入力から距離 $`r`$ にあるエルミート行列が
生成する環の厳密な分解である。

(c) $`2r<\gamma_{\mathcal B}`$ ならば、すべての $`\theta\in(2r,\gamma_{\mathcal B}-2r)`$ について入力の許容誤差 $`\theta`$ の可換子環は次元 $`d`$ をもち、
$`\tilde{\mathcal B}'`$ から $`4r/\gamma_{\mathcal B}`$ 以内にある；また $`(\sum_j\|H^0_j-\tilde H_j\|^2)^{1/2}<\gamma_{\mathcal B}/2-r`$ を満たすエルミート
行列 $`H^0_j`$ が生成する環の可換子環の次元は $`\le d`$ である。

*証明.* $`\tilde{\mathcal B}`$ は $`X\mapsto X^*`$ で閉じており、これは $`\mathrm{Re}\langle\cdot,\cdot\rangle_{HS}`$ を保つ反線形写像だから、
$`\tilde{\mathcal B}^\perp`$ をそれ自身に写し、$`E`$ と可換である；よって $`H_j`$ と $`a_{j,c}`$ はエルミートであり、$`\|\tilde H_j-H_j\|\le\|\tilde H_j-H_j\|_2`$。
(a) $`\tilde V^*X\tilde V`$ をブロック $`X_{ce}\in M_{k_cm_c\times k_em_e}=M_{k_c\times k_e}\otimes M_{m_c\times m_e}`$（クロネッカー順序）に分ける；
$`X\mapsto\tilde V^*X\tilde V`$ は $`\|\cdot\|_2`$ を保ち、これらのブロックは直交し、$`\tilde V^*[H_j,X]\tilde V`$ の $`(c,e)`$ ブロックは
$`(a_{j,c}\otimes1)X_{ce}-X_{ce}(a_{j,e}\otimes1)`$ であり、$`X_{ce}=Y\otimes Z`$ に対しては $`\Phi_{ce}(Y)_j\otimes Z`$ である。よって交換子写像は
$`(c,e)`$ ブロック上で $`\Phi_{ce}\otimes\mathrm{id}`$ として作用し、その $`T^*T`$ は $`\bigoplus_{c,e}\Phi_{ce}^*\Phi_{ce}\otimes1`$ であり、その特異値は
各 $`\Phi_{ce}`$ の特異値をそれぞれ $`m_cm_e`$ 回繰り返したものである。$`\Phi_{cc}(1)=0`$ だから少なくとも $`\sum_cm_c^2=d`$ 個が零であり、
ちょうど $`d`$ 個であるのは、すべての $`c`$ で $`\sigma_2(\Phi_{cc})>0`$ かつ $`c\neq e`$ で $`\sigma_1(\Phi_{ce})>0`$ のとき、かつそのときに限る；
そのとき $`\sigma_{d+1}`$ はこれらの数の最小値であり、そうでなければ両辺とも $`0`$ である。(b) $`H_j\in\tilde{\mathcal B}`$ は次元 $`d`$ の
$`\tilde{\mathcal B}'`$ と可換である（定理 4.3、命題 4.9 の証明と同様）から、$`\tilde{\mathcal B}'\subseteq\{H_j\}'`$；$`\gamma_{\mathcal B}>0`$ ならば零となる
特異値は高々 $`d`$ 個だから、$`\{H_j\}'=\tilde{\mathcal B}'`$ かつ $`\{H_j\}''=\tilde{\mathcal B}''=\tilde{\mathcal B}`$（定理 2.2）。(c) 厳密な生成元 $`H_j`$
（ギャップ $`\gamma_{\mathcal B}`$、(b) により可換子環は $`\tilde{\mathcal B}'`$）と、距離 $`\eta\le r`$ にある摂動された $`\tilde H_j`$ に定理 4.25 を
適用する。$`H^0_j`$ については $`(\sum_j\|H^0_j-H_j\|^2)^{1/2}<\gamma_{\mathcal B}/2`$ であり、定理 4.25(a) によりその交換子写像の
$`\sigma_{d+1}>0`$ である。$`\square`$

こうして出力は自らを保証する：後退誤差 $`r`$ が小さく、ギャップ $`\gamma_{\mathcal B}`$ が $`2r`$ よりはるかに大きければ、計算された構造は
わずかに変更した生成元に対して厳密であり、$`(2r,\gamma_{\mathcal B}-2r)`$ のどの閾値も同じ可換子環の次元を与え、より多くの対称性をもつ
モデルは入力から $`\gamma_{\mathcal B}/2-r`$ 以内には存在しない。

*ライブラリでは：* 証明書は `Decomposition.certificate()` である：(a) のようにブロックから求めた $`r`$ と $`\gamma_{\mathcal B}`$、窓
$`(2r,\gamma_{\mathcal B}-2r)`$、実際に用いた閾値の両側の特異値を返し、その閾値が窓に入っていなければ警告する — [チュートリアル、4.9節](../../tutorial/ja/04_symmetries.md#s-4-9)。

残りのステップ（ランダムな $`K`$ による分割とランダムな $`Y`$ によるグループ化）には定量的な一般性が必要である。ライブラリのサンプルは
スケールし直したガウスベクトルである（4.3節、ステップ 1）：可換子環の正規直交基底 $`(B_r)`$ と独立な標準複素ガウス変数 $`c_r`$ により、
$`Y=\sum_rc_rB_r`$ は $`\mathcal A'`$ の **標準複素ガウスベクトル** であり（その分布は基底によらない、命題 A.16(a)）、独立なコピー $`Y'`$ に対し
$`K=Y'+Y'^*=\sqrt2\,G`$ である。ここで $`G`$ は内積 $`\langle X,Y\rangle=\mathrm{Tr}(XY)`$ をもつ実空間 $`\mathcal A'_h`$ の **標準ガウスベクトル** である
（$`\sqrt2\,Y'`$ は $`\mathrm{Re}\langle\cdot,\cdot\rangle_{HS}`$ をもつ実空間 $`\mathcal A'`$ の標準ガウスベクトルであり、$`X\mapsto\frac12(X+X^*)`$ は
$`\mathcal A'_h`$ への直交射影である；命題 A.16(a)、A.14(b)）。群平均によるサンプリングでも同じことが成り立つ：$`Z+Z^*=\sqrt2\,H`$ で、
$`H`$ は GUE 行列、すなわち $`\mathrm{Herm}_n`$ の標準ガウスベクトルであり、その $`\mathcal A'_h`$ への直交射影が群平均である。また Ginibre
行列の群平均は $`\mathcal A'`$ の標準複素ガウスベクトルである。

<a id="r-4-27"></a>

**命題 4.27（可換子環のランダムな元）.** $`G`$ と $`Z`$ を独立とし、それぞれ $`\mathcal A'_h`$ の標準ガウスベクトル、$`\mathcal A'`$ の
標準複素ガウスベクトルとする。

(a) $`V^*GV=\bigoplus_c1_{k_c}\otimes k_c^{-1/2}H_c`$、$`V^*ZV=\bigoplus_c1_{k_c}\otimes k_c^{-1/2}Z_c`$ であり、$`H_c`$ はサイズ $`m_c`$ の GUE
行列、$`Z_c`$ はサイズ $`m_c`$ の Ginibre 行列で、すべて独立である；$`\mathbb E\|G\|_2^2=\mathbb E\|Z\|_2^2=d`$。

(b) $`\kappa_c=k_c^{-1/2}H_c`$ とし、$`\delta_\times=\min|\lambda_i(\kappa_c)-\lambda_j(\kappa_e)|`$（$`c\neq e`$ とすべての $`i,j`$ にわたる最小）とする。
$`x>0`$ に対し
```math
\mathbb P(\delta_\times<x)\le\sqrt{\tfrac2\pi}\;x\sum_{c<e}m_cm_e\Big(\frac1{k_cm_c}+\frac1{k_em_e}\Big)^{-1/2}.
```

(c) $`m_c\ge2`$ のとき $`\delta_c`$ を $`\kappa_c`$ の最小準位間隔とする。このとき $`k_c^{1/2}x\le1`$ に対し $`\mathbb P(\delta_c<x)\le C_{m_c}k_c^{3/2}x^3`$
であり、$`C_{m_c}`$ は命題 C.9 の定数（$`\beta=2`$）である；$`m_c=2`$ ではすべての $`x`$ について $`C_2=1/(6\sqrt\pi)`$ でよい。

(d) $`X`$ を次元 $`D`$ の実内積空間の標準ガウスベクトル、$`W`$ を $`X`$ と独立な確率変数、$`f\ge0`$ を $`t>0`$ に対し
$`f(w,tx)=tf(w,x)`$ を満たすボレル関数とする。このとき $`\mathbb P(f(W,X)<\varepsilon|X|)\le2\,\mathbb P(f(W,X)<\varepsilon(2D)^{1/2})`$ である。
$`X=G`$ では $`D=d`$；$`X=Z`$ では $`(2d)^{1/2}`$ で同じことが成り立つ。

(e) $`G`$ が補題 4.6 の仮定を満たし、$`\Pi,\Pi'`$ が同じ $`\kappa_c`$ の相異なる二つの固有値に対する $`G`$ のスペクトル射影ならば、$`G`$ を
与えたときの条件付きで $`\|\Pi Z\Pi'\|_2^2`$ は $`(0,\infty)`$ 上の密度 $`e^{-t}`$ をもつ：$`\mathbb P(\|\Pi Z\Pi'\|_2<x)=1-e^{-x^2}\le x^2`$。

したがって、$`\delta(G)=\min(\delta_\times,\min_c\delta_c)`$（$`G`$ の相異なる固有値が $`N=\sum_cm_c`$ 個のときその間の最小距離、そうでなければ
$`0`$）と、同じ直和成分内の相異なるスペクトル射影の組にわたる $`\beta(G,Z)=\min\|\Pi_l Z\Pi_{l'}\|_2`$ に対し
```math
\mathbb P\big(\delta(G)<\varepsilon\|G\|_2\big)\le\frac{4\sqrt d}{\sqrt\pi}\,\varepsilon\sum_{c<e}\frac{m_cm_e}{(\frac1{k_cm_c}+\frac1{k_em_e})^{1/2}}
+2\sum_{c:m_c\ge2}C_{m_c}(2dk_c)^{3/2}\varepsilon^3=:p_\delta(\varepsilon),\qquad
\mathbb P\big(\beta(G,Z)<\beta'\|Z\|_2\big)\le4d\beta'^2\sum_cm_c(m_c-1)=:p_\beta(\beta')
```
が成り立つ（前者は $`\varepsilon(2dk_c)^{1/2}\le1`$ のとき）。オーダー $`\varepsilon`$ は改善できない：異なるセクターの準位は反発しない（定理 5.15）が、
同じセクターの準位は反発する（オーダー $`\varepsilon^3`$）。

*証明.* (a) 行列 $`k_c^{-1/2}V(1\otimes b)V^*`$（$`b`$ は $`\mathrm{Tr}(xy)`$ に関する $`\mathrm{Herm}_{m_c}`$ の正規直交基底を動き、$`c=1,\dots,C`$）は
$`\mathcal A'_h`$ の正規直交基底をなす（$`\mathrm{Tr}((1\otimes b)(1\otimes b'))=k_c\mathrm{Tr}(bb')`$）から、命題 A.14(b) により
$`\kappa_c=k_c^{-1/2}\sum_bg_bb`$ であり、係数 $`g_b`$ は独立に $`N(0,1)`$ に従う。基底 $`e_{ii}`$、$`(e_{ij}+e_{ji})/\sqrt2`$、$`i(e_{ij}-e_{ji})/\sqrt2`$
に対し $`\sum_bg_bb`$ の成分は独立で $`h_{ii}\sim N(0,1)`$、$`\mathrm{Re}\,h_{ij},\mathrm{Im}\,h_{ij}\sim N(0,\frac12)`$：これは GUE 行列である（命題 A.17）。
同様に $`\mathcal A'`$ の正規直交基底 $`k_c^{-1/2}V(1\otimes e_{jj'})V^*`$ から、成分が独立な標準複素ガウス変数である $`Z_c`$ を得る。最後に
$`\mathbb E\|G\|_2^2=\sum\mathbb Eg_b^2=d`$ であり、$`Z`$ も同様である。

(b) $`c\neq e`$ を固定し、$`s=(\frac1{k_cm_c}+\frac1{k_em_e})^{1/2}`$ とおき、$`u\in\mathcal A'_h`$ を直和成分 $`c`$ でブロック $`\alpha1`$、直和成分 $`e`$
で $`-\beta1`$、他では $`0`$ とする（$`\alpha=1/(k_cm_cs)`$、$`\beta=1/(k_em_es)`$）。このとき $`\|u\|_2^2=\alpha^2k_cm_c+\beta^2k_em_e=1`$、
$`\alpha+\beta=s`$ である。命題 A.14(b) により $`G=G_0+tu`$ で、$`t=\langle u,G\rangle\sim N(0,1)`$ は $`G_0`$ と独立である。ブロックは
$`\kappa_c=\kappa^0_c+t\alpha`$、$`\kappa_e=\kappa^0_e-t\beta`$ だから $`\lambda_i(\kappa_c)-\lambda_j(\kappa_e)=a+ts`$ であり、
$`a=\lambda_i(\kappa^0_c)-\lambda_j(\kappa^0_e)`$ は $`G_0`$ の関数である（系 A.11）。$`a`$ を固定すると $`\{t:|a+ts|<x\}`$ は長さ $`2x/s`$ の区間で、
その確率は $`\le2x/(s\sqrt{2\pi})`$；独立性（Fubini）により、これが $`\mathbb P(|\lambda_i(\kappa_c)-\lambda_j(\kappa_e)|<x)`$ を抑える。
$`c<e`$、$`i\le m_c`$、$`j\le m_e`$ について和をとる。

(c) $`\delta_c<x`$ であるのは $`H_c`$ の最小準位間隔が $`<k_c^{1/2}x`$ のとき、かつそのときに限る；命題 C.9。$`m_c=2`$ では間隔の密度は
$`s^2e^{-s^2/4}/(2\sqrt\pi)\le s^2/(2\sqrt\pi)`$（命題 5.13 の証明）であり、その $`(0,y)`$ 上の積分は $`\le y^3/(6\sqrt\pi)`$ である。

(d) 正規直交座標で $`X`$ は $`\mathbb R^D`$ の標準ガウスベクトルである；定理 A.20(a) により $`|X|`$ は $`X/|X|`$ と独立であり、したがって
$`(W,X/|X|)`$ と独立である。マルコフの不等式により $`\mathbb P(|X|^2\ge2D)\le\mathbb E|X|^2/2D=\frac12`$。$`f(W,X)=|X|f(W,X/|X|)`$ だから
$`\mathbb P(f(W,X)<\varepsilon(2D)^{1/2})\ge\mathbb P(f(W,X/|X|)<\varepsilon)\,\mathbb P(|X|^2<2D)\ge\frac12\mathbb P(f(W,X)<\varepsilon|X|)`$。$`Z`$ については
$`\sqrt2Z`$ が次元 $`2d`$ の標準ガウスベクトルであり、$`\sqrt2Z`$ に対する $`\varepsilon(4d)^{1/2}`$ は $`Z`$ に対する $`\varepsilon(2d)^{1/2}`$ である。

(e) (a) と補題 4.6 により、$`V^*\Pi V`$ と $`V^*\Pi'V`$ は $`H_c`$ の正規直交固有ベクトル $`e,e'`$（$`G`$ のみに依存する）について
$`\mathbb C^{k_c}\otimes e`$ と $`\mathbb C^{k_c}\otimes e'`$ への射影である。すると $`\|\Pi Z\Pi'\|_2^2=k_c|\langle e,k_c^{-1/2}Z_ce'\rangle|^2=|\langle e,Z_ce'\rangle|^2`$
であり、$`e,e'`$ を固定すると $`\langle e,Z_ce'\rangle`$ は単位ベクトル $`|e\rangle\langle e'|`$ に沿った $`Z_c`$ の座標、すなわち標準複素ガウス変数で
あり（命題 A.16(a),(b)）、その絶対値の2乗は密度 $`e^{-t}`$ をもつ（定理 A.20(b)、$`N=1`$）。$`Z`$ は $`G`$ と独立だから、Fubini の定理により
条件付きの主張が得られる。*帰結.* 自明な $`W`$ に対して (d)、次に (b)、(c) と和集合の評価を $`x=\varepsilon(2d)^{1/2}`$ で用いる；また
$`W=G`$ に対して (d)、次に (e) と $`\sum_cm_c(m_c-1)`$ 個の順序対についての和集合の評価を $`x^2=2d\beta'^2`$ で用いる。系 A.11 と命題 A.13
により、これらの事象はボレルである。$`\square`$

アルゴリズムでは、可換子環のサンプルは $`\mathcal A'`$ からではなく計算された可換子環から取られる；それらを $`\mathcal A'`$ への射影と比較する。

<a id="r-4-28"></a>

**定理 4.28（分解の安定性）.** $`\tilde{\mathcal A}'\subseteq M_n`$ を $`X\mapsto X^*`$ で閉じた次元 $`d`$ の部分空間とし、その直交射影 $`\tilde\Pi`$
は $`\|\tilde\Pi-\Pi\|\le\varepsilon_0<1`$ を満たすとする（定理 4.25(b) により、窓の中の許容誤差での可換子環に対しては $`\varepsilon_0=4\eta/\gamma`$）。
$`\tilde K\in\tilde{\mathcal A}'_h`$、$`\tilde Y\in\tilde{\mathcal A}'`$、$`K=\Pi(\tilde K)\in\mathcal A'_h`$、$`Y=\Pi(\tilde Y)\in\mathcal A'`$ とする。$`K`$ は補題 4.6 の仮定を
満たし、固有値 $`\mu_1<\dots<\mu_N`$、直和成分 $`c(l)`$ 内の階数 $`k_{c(l)}`$ のスペクトル射影 $`P_l`$ をもつとし、$`\delta=\min_l(\mu_{l+1}-\mu_l)`$ とする；
$`\beta=\min\{\|P_lYP_{l'}\|_2:l\neq l',\,c(l)=c(l')\}`$（そのような組がなければ $`\infty`$）、$`h=\max_cm_c-1`$、そして
```math
\rho=\varepsilon_0\|\tilde K\|_2,\qquad\varepsilon_1=\frac{2\rho}\delta,\qquad\zeta=(\varepsilon_0+2\varepsilon_1)\|\tilde Y\|_2,\qquad
\xi=\frac{k_{\max}^{1/2}(\varepsilon_0+4\varepsilon_1)\|\tilde Y\|_2}\beta
```
とおく。

(a) $`4\rho<\delta`$ ならば、任意の $`\tau\in[2\rho,\delta-2\rho)`$ について、$`\tilde K`$ の昇順の固有値を隣り合う二つの差が $`\tau`$ を超える
ところで切ると $`N`$ 個のグループができる；$`l`$ 番目は $`k_{c(l)}`$ 個の元をもち、そのスペクトル射影 $`\tilde P_l`$ は
$`\|\tilde P_l-P_l\|\le\varepsilon_1`$ を満たす。

(b) さらに $`2\zeta<\beta`$ ならば、任意の $`\theta'\in[\zeta,\beta-\zeta)`$ について、同じ大きさの二つのグループ $`l\neq l'`$ が
$`\|\tilde P_l\tilde Y\tilde P_{l'}\|_2>\theta'`$ を満たすのは $`c(l)=c(l')`$ のとき、かつそのときに限る。よってステップ 3 の連結な類は直和成分で
あり、アルゴリズムは厳密な $`(k_c,m_c)`$ を返す。

(c) さらに $`3^h\xi\le1`$ ならば、$`\tilde V`$ をステップ 3–4 のように $`\mathrm{ran}\,\tilde P_l`$ の正規直交基底から組み立て、各直和成分の
コピーを深さ $`\le h`$ の任意の木に沿って極分解の因子で整列させたものとする。このとき、すべての $`a\in\mathcal A`$ について
$`\hat V^*a\hat V=\bigoplus_c\rho_c(a)\otimes1_{m_c}`$ を満たすユニタリ $`\hat V`$（定理 4.3 の意味の分解）が存在し、
$`\|\tilde V-\hat V\|\le N^{1/2}(2\varepsilon_1+3^h\xi)`$ である。

(d) ライブラリと同様に、$`\tilde K`$ を $`\tilde{\mathcal A}'_h`$ の標準ガウスベクトルの正の定数倍、$`\tilde Y`$ をそれと独立な $`\tilde{\mathcal A}'`$ の
標準複素ガウスベクトルの正の定数倍とする。（ライブラリはこれらを正規化するが、$`\tilde K`$ と $`\tilde Y`$ を正の定数倍し、$`\tau`$ と $`\theta'`$ をそれに合わせてスケールしても、(a)–(c) の仮定と結論は変わらない。）組からなる任意のボレル集合 $`\mathcal E`$ に対し、命題 4.27 と同じ定数倍 $`c,c'`$ について
$`\mathbb P((K,Y)\in\mathcal E)\le(1-\varepsilon_0^2)^{-3d/2}\,\mathbb P((cG,c'Z)\in\mathcal E)`$ が成り立つ。したがって $`\varepsilon,\beta'\in(0,1]`$ かつ
$`\varepsilon_0\le\varepsilon\beta'/(14\cdot3^hk_{\max}^{1/2})`$ ならば、確率 $`\le(1-\varepsilon_0^2)^{-3d/2}(p_\delta(\varepsilon)+p_\beta(\beta'))`$ の事象を除いて
(a)–(c) の仮定が成り立ち、$`\|\tilde V-\hat V\|\le17\cdot3^h(Nk_{\max})^{1/2}\varepsilon_0/(\varepsilon\beta')`$ である。$`\varepsilon_0^2\le1/(3d)`$ ならば
その因子は $`\le2`$ である。

*証明.* $`\Pi`$ は $`X\mapsto X^*`$ と可換である（命題 4.26 の $`E`$ と同様）から $`K`$ はエルミートである。$`X\in\mathrm{ran}\,\tilde\Pi`$ に対し
$`X-\Pi X=(1-\Pi)\tilde\Pi X=(\tilde\Pi-\Pi)\tilde\Pi X`$ だから $`\|X-\Pi X\|_2\le\varepsilon_0\|X\|_2`$ かつ $`\|\Pi X\|_2\le\|X\|_2`$。(a) よって
$`\|\tilde K-K\|\le\|\tilde K-K\|_2\le\rho`$。各 $`\mu_l`$ に補題 4.24(b)（$`g\ge\delta>4\rho`$、$`\Delta=\rho`$）を適用すると、$`\tilde K`$ は
$`[\mu_l-\rho,\mu_l+\rho]`$ にちょうど $`k_{c(l)}`$ 個の固有値をもつ；これらの区間は交わらず、一つの区間内で隣り合う固有値の差は
$`\le2\rho\le\tau`$、異なる区間の固有値の差は $`\ge\delta-2\rho>\tau`$ である；また $`\|\tilde P_l-P_l\|\le2\rho/\delta`$。(b)
$`\|AXB\|_2\le\|A\|\|X\|_2\|B\|`$ により
$`|\|\tilde P_l\tilde Y\tilde P_{l'}\|_2-\|P_lYP_{l'}\|_2|\le\|\tilde Y-Y\|_2+\|(\tilde P_l-P_l)Y\|_2+\|Y(\tilde P_{l'}-P_{l'})\|_2\le\zeta`$ であり、一方
$`c(l)\neq c(l')`$ なら $`P_lYP_{l'}=0`$（補題 4.8(b)）、そうでなければ $`\|P_lYP_{l'}\|_2\ge\beta`$ である。よって一つの直和成分のグループは
すべて互いに結ばれ、他のどのグループとも結ばれず、直和成分 $`c`$ は大きさ $`k_c`$ のグループを $`m_c`$ 個与える。

(c) $`\tilde Q_l`$ を $`\mathrm{ran}\,\tilde P_l`$ の基底、$`Q_l`$ を補題 4.24(c) による $`\mathrm{ran}\,P_l`$ の基底とする（$`\varepsilon_1<\frac12`$）。すると
$`\|\tilde Q_l-Q_l\|\le2\varepsilon_1`$。補題 4.8(b) により、一つの直和成分 $`c`$ 内の $`l,l'`$ に対し $`T_{ll'}=Q_l^*YQ_{l'}=\lambda_{ll'}U_{ll'}`$ であり、
$`U_{ll'}`$ はユニタリ、$`\lambda_{ll'}=k_c^{-1/2}\|P_lYP_{l'}\|_2\ge\beta k_{\max}^{-1/2}`$ である。また $`\tilde T_{ll'}=\tilde Q_l^*\tilde Y\tilde Q_{l'}`$ は
$`\|\tilde T_{ll'}-T_{ll'}\|\le\varepsilon_0\|\tilde Y\|_2+4\varepsilon_1\|Y\|_2\le\xi\lambda_{ll'}`$ を満たす。アルゴリズムは根で $`\tilde W_{l_0}=1`$ とし、
$`l`$ から $`l'`$ への辺に沿って $`\tilde W_{l'}=\mathrm{pol}(\tilde W_l^*\tilde T_{ll'})^*`$ とする；$`T_{ll'}`$ から同じように $`W_l`$ を定めると
$`W_{l'}=U_{ll'}^*W_l`$ である。補題 4.8(b) と木に沿った帰納法により、直和成分 $`c`$ のコピーの基底 $`Q_lW_l`$ は同じ行列
$`(Q_lW_l)^*a(Q_lW_l)`$ をもつ。よって $`\tilde Q_l\tilde W_l`$ が $`\tilde V`$ の中で並んでいるのと同じように $`Q_lW_l`$ を並べると、定理 4.3 の
証明の存在部分のような分解 $`\hat V`$ が得られる。$`e_l=\|\tilde W_l-W_l\|`$ とおく。すると
$`\|\tilde W_l^*\tilde T_{ll'}-W_l^*T_{ll'}\|\le e_l\|\tilde T_{ll'}\|+\|\tilde T_{ll'}-T_{ll'}\|\le\lambda_{ll'}(e_l(1+\xi)+\xi)`$ であり、これが
$`<\lambda_{ll'}`$ ならば、補題 4.24(d)（$`W_l^*T_{ll'}=\lambda_{ll'}W_l^*U_{ll'}`$ として）により $`e_{l'}\le2e_l(1+\xi)+2\xi\le3e_l+2\xi`$。深さ $`t`$ に
関する帰納法で $`e_l\le\xi(3^t-1)`$：帰納のステップは $`3\xi(3^t-1)+2\xi=\xi(3^{t+1}-1)`$ であり、$`t\le h-1`$ では
$`e_l(1+\xi)+\xi\le\xi3^t(1+\xi)\le\frac13\cdot\frac43<1`$（$`\xi\le3^{-h}\le\frac13`$ による）だから適用できる。よって
$`\|\tilde Q_l\tilde W_l-Q_lW_l\|\le2\varepsilon_1+3^h\xi`$；行列 $`\tilde V-\hat V`$ はこれら $`N`$ 個の列ブロックからなり、
$`|\sum_lA_lx_l|\le\sum_l\|A_l\||x_l|\le(\sum_l\|A_l\|^2)^{1/2}|x|`$ により評価を得る。

(d) $`\tilde{\mathcal A}'_h`$ と $`\mathcal A'_h`$ の正規直交座標で（どちらも実次元 $`d`$：随伴で閉じた部分空間は $`\mathcal S_h\oplus i\mathcal S_h`$）、
$`\Pi`$ は線形写像 $`B:\tilde{\mathcal A}'_h\to\mathcal A'_h`$ に制限され、$`(1-\varepsilon_0^2)|x|^2\le|Bx|^2\le|x|^2`$ を満たす。よって $`B`$ は全単射で、
$`|B^{-1}y|\ge|y|`$、$`|\det B|=\det(B^TB)^{1/2}\ge(1-\varepsilon_0^2)^{d/2}`$（定理 A.8）。補題 A.3(c) により $`K/c=B(\tilde K/c)`$ は密度
$`\varphi_d(B^{-1}y)/|\det B|\le(1-\varepsilon_0^2)^{-d/2}\varphi_d(y)`$ をもつ。ここで $`\varphi_d`$ は $`G`$ の密度である。次元 $`2d`$ の実空間
$`\tilde{\mathcal A}'`$ 上の同じ議論で $`Y`$ に対する因子 $`(1-\varepsilon_0^2)^{-d}`$ を得て、$`K,Y`$ は独立である。命題 4.27 の二つの事象（$`K`$ と
$`Y`$ について；これらは定数倍で不変）の補集合上で、$`(1-\varepsilon_0^2)^{1/2}\ge\frac{\sqrt3}2`$ を用いると：$`\delta\ge\varepsilon\|K\|_2\ge\frac{\sqrt3}2\varepsilon\|\tilde K\|_2`$
だから $`\varepsilon_1\le\frac4{\sqrt3}\frac{\varepsilon_0}\varepsilon<\frac12`$；$`\beta\ge\beta'\|Y\|_2\ge\frac{\sqrt3}2\beta'\|\tilde Y\|_2`$ であり、
$`2\zeta\le2(1+\frac8{\sqrt3\varepsilon})\varepsilon_0\|\tilde Y\|_2\le11.24\,\frac{\varepsilon_0}\varepsilon\|\tilde Y\|_2\le\frac{11.24}{14}\beta'\|\tilde Y\|_2<\beta`$；そして
$`3^h\xi\le3^hk_{\max}^{1/2}\frac{2}{\sqrt3\beta'}(1+\frac{16}{\sqrt3\varepsilon})\varepsilon_0\le11.83\cdot3^hk_{\max}^{1/2}\frac{\varepsilon_0}{\varepsilon\beta'}\le\frac{11.83}{14}<1`$。$`\|\tilde V-\hat V\|`$ の評価は
$`2\varepsilon_1\le4.62\,\varepsilon_0/\varepsilon`$、$`3^h\xi`$ の評価、および $`4.62+11.83<17`$ から従う。最後に、$`(1-\frac1m)^{-m}`$ は $`m\ge2`$ で減少するので
$`(1-\frac1{3d})^{-3d/2}\le(\frac23)^{-3/2}<2`$。$`\square`$

まとめると：雑音がギャップに比べて小さく、$`\eta\le\gamma\varepsilon\beta'/(56\cdot3^hk_{\max}^{1/2})`$ であり、閾値が定理 4.25 と 4.28 の窓に
入っていれば、アルゴリズムは確率 $`O(\varepsilon+\beta'^2)`$ の事象を除いて $`(k_c,m_c)`$ を復元し、厳密な分解から $`O(\eta/\gamma)`$ 以内の $`V`$ を
与える。定数は最悪の場合の評価であり、典型的な誤差ははるかに小さい。ライブラリは $`\eta`$ と $`\gamma`$ を知りえないので、可換子環の
サンプルの相対誤差を $`\nu=\max(\sigma_d(\tilde L),\epsilon_{\rm mach}\max(1,\sigma_{\max}))/\sigma_{d+1}(\tilde L)`$ で見積もり（補題 4.24(a) により
$`\varepsilon_0`$ と同じオーダー、$`\sigma_d(\tilde L)\le\Delta`$ だから）、クラスタリングの下限閾値と結合の閾値を $`\min(10^3\nu,\nu^{1/2})`$ に
置く（後者は少なくとも $`\epsilon_{\rm mach}^{1/2}`$）。これは雑音と構造の単位スケールの間である；その後、命題 4.26 が結果を検査する。

*計算：* 雑音を含む生成元での定理 4.25（特異値のギャップと許容誤差の比較、および (c) の例）、証明書についての命題 4.26、
モンテカルロ法による命題 4.27 — [チュートリアル、4.9節](../../tutorial/ja/04_symmetries.md#s-4-9)。

<a id="x-4"></a>

## 演習問題

1. $`(\mathbb C^d)^{\otimes L}`$ 上の $`S_L`$ の表現について、
   $`\chi_\pi(\sigma)=d^{\#\text{cycles}(\sigma)}`$ かつ $`\dim\pi(S_L)'=\binom{d^2+L-1}L`$ であることを示せ。$`20`$（2.4節）と
   $`35`$（四つのスピンに対する $`=\dim\{S_a\}''`$）を確かめよ。
2. （*クレブシュ・ゴルダン.*）スピン表現 $`S^{(j_1)}_a\otimes1+1\otimes S^{(j_2)}_a`$ に対して
   $`j_1\otimes j_2\cong\bigoplus_jj`$、$`j=|j_1-j_2|,|j_1-j_2|+1,\dots,j_1+j_2`$ を示せ——これは第7章で触れた $`SU(2)`$ の非アーベル的な融合則である——。
3. $`\theta/2\pi`$ を無理数とする。単一の $`U_\theta=e^{i\theta Q}`$ に対して $`\{U_\theta\}'=U(1)'`$ であるが、
   $`n`$ 個の量子ビットに対して $`\{e^{i\pi\hat N}\}'\cong M_{2^{n-1}}\oplus M_{2^{n-1}}`$ であることを示せ。
4. （*ゼーマン分裂.*）$`H`$ を $`L\ge1`$ 個のスピン上のハイゼンベルク・ハミルトニアンとし、$`H_h=H+hS_z`$ とする。
   $`H_h`$ は $`S_z`$ および $`\mathbf S^2`$ と可換であるが、$`h\neq0`$ のとき $`S_x`$ とは可換でないこと、そしてその固有値が
   $`H_S`$ の固有値 $`E`$ に対して $`E+hM`$、$`M=-S,\dots,S`$ であることを示せ。
5. （*四元数群.*）$`I=(i,0)`$、$`J=(0,1)`$、$`K=(0,i)`$（E.6節）として $`Q_8=\{\pm1,\pm I,\pm J,\pm K\}\subseteq\mathbb H`$ とし、
   $`\pi`$ を $`\mathbb R^8`$ 上のその正則表現 $`\pi(g)e_h=e_{gh}`$ とする。$`\mathcal A=\mathrm{span}\,\pi(Q_8)`$ は
   $`k=m=1`$ の実型のセクターを四つと、$`k=m=2`$ の四元数型のセクターを一つもち、したがって
   $`\mathcal A_{\mathbb R}\cong\mathbb R^4\oplus\mathbb H`$ であることを示せ。

<details>
<summary>解答</summary>

1. 補題 B.1(c) により、$`\pi(\sigma)`$ は置換行列であり、$`|b_1\dots b_L\rangle`$ を固定するのは $`b`$ が $`\sigma`$ の各巡回上で一定であるときかつそのときに限る。そのような $`b`$ は $`d^{\#\text{cycles}(\sigma)}`$ 個ある。系 4.5(c) により
   $`\dim\pi(S_L)'=\frac1{L!}\sum_\sigma(d^2)^{\#\text{cycles}(\sigma)}=\frac{x(x+1)\cdots(x+L-1)}{L!}`$（$`x=d^2`$）
   であり（恒等式 $`\sum_\sigma x^{\#\text{cycles}(\sigma)}=x(x+1)\cdots(x+L-1)`$ は第5章、演習問題 2 の解答で証明される）、これは $`\binom{d^2+L-1}L`$ である
   （命題 B.3 は多重集合を数えることで同じ次元を証明する）：$`\binom63=20`$、$`\binom74=35`$。$`35`$ が
   $`\dim\{S_a\}''`$ でもあることは、$`L=4`$ の定理 4.14(b) である：$`\{S_a\}''=\pi(S_4)'`$。

2. $`S_a=S^{(j_1)}_a\otimes1+1\otimes S^{(j_2)}_a`$ はスピン表現であり（4.5節）、ユニタリ
   $`x\otimes y\mapsto y\otimes x`$ はそれを $`(j_2,j_1)`$ に対するものと絡ませる。よって $`j_1\ge j_2`$ としてよい。$`S_z`$ は
   $`|m_1\rangle\otimes|m_2\rangle`$ 上で固有値 $`m_1+m_2`$ をもつので、その固有値は $`j_1+j_2+\mathbb Z`$ に属し、そのような $`M`$ に対して $`d_M`$ は
   $`|M-m_2|\le j_1`$ を満たす $`m_2\in\{-j_2,\dots,j_2\}`$ の個数である（このとき $`m_1=M-m_2\in j_1+\mathbb Z`$ は自動的に成り立つ）。
   $`M\notin j_1+j_2+\mathbb Z`$ に対しては $`d_M=0`$ である。$`M\ge0`$ に対しては条件は $`M-j_1\le m_2\le M+j_1`$ であり、$`M+j_1\ge j_2`$ であるから、$`m_2`$ は
   $`\max(-j_2,M-j_1)`$ から $`j_2`$ まで動く：$`0\le M\le j_1-j_2`$ に対して $`d_M=2j_2+1`$、
   $`j_1-j_2\le M\le j_1+j_2`$ に対して $`d_M=j_1+j_2-M+1`$、$`M>j_1+j_2`$ に対して $`d_M=0`$ である。系 4.12 により、スピン $`j\ge0`$ の重複度は
   $`m_j=d_j-d_{j+1}`$ である。これは $`j\notin j_1+j_2+\mathbb Z`$ または $`j>j_1+j_2`$ ならば $`0`$ であり、$`m_{j_1+j_2}=1-0=1`$ である。
   $`j_1-j_2\le j<j_1+j_2`$（$`j\in j_1+j_2+\mathbb Z`$ であるから $`j+1\le j_1+j_2`$）に対しては
   $`m_j=(j_1+j_2-j+1)-(j_1+j_2-j)=1`$ であり、$`0\le j<j_1-j_2`$（したがって $`j+1\le j_1-j_2`$）に対しては $`m_j=(2j_2+1)-(2j_2+1)=0`$ である。
   よって $`j_1\otimes j_2\cong\bigoplus_{j=j_1-j_2}^{j_1+j_2}j`$ であり、各スピンはちょうど一回現れる。$`j_2=\frac12`$ に対しては、これは
   定理 4.16(d) であり、そちらはさらに係数も与える。

3. $`U_\theta=\sum_qe^{iq\theta}P_q`$ であり、$`e^{iq\theta}=e^{iq'\theta}`$ であることは $`(q-q')\theta\in2\pi\mathbb Z`$ と同値で、
   $`\theta/2\pi`$ が無理数ならばこれは $`q=q'`$ を強いる。よって $`U_\theta`$ の固有空間は $`P_q`$ の像であり、作用素が
   正規行列 $`U_\theta`$ と可換であることは、それが $`U_\theta`$ の固有空間を保つことと同値である（定理 1.5(d) と、2.2節と同様のラグランジュ補間）。すなわち $`X=\sum_qP_qXP_q`$ と同値であり、命題 4.10 によりこれは $`U(1)'`$ である。この議論が用いるのは、現れる電荷 $`q`$ に対する
   $`e^{iq\theta}`$ が相異なることだけである。
   $`\theta=\pi`$ に対しては、固有値 $`(-1)^w`$ はハミング重みの偶奇しか見ない：次元
   $`2^{n-1}`$ の二つの固有空間があり（$`n\ge1`$ のとき。重みが偶数の配置と奇数の配置であり、最初のビットを反転することがそれらの間の全単射を与える）、可換子環は $`M_{2^{n-1}}\oplus M_{2^{n-1}}`$ である。

4. $`S_z`$ と $`\mathbf S^2`$ は $`H`$（命題 4.15(a)）および $`S_z`$ と可換であるから、$`H_h`$ と可換である。また
   $`[H_h,S_x]=h[S_z,S_x]=ihS_y\neq0`$ である。命題 4.13 の分解において
   $`V^*S_zV=\bigoplus_SS^{(S)}_z\otimes1_{m_S}`$ かつ $`V^*HV=\bigoplus_S1_{2S+1}\otimes H_S`$（命題 4.15）であるから、
   スピン $`S`$ の直和成分上で $`V^*H_hV=hS^{(S)}_z\otimes1+1\otimes H_S`$ となる。これは可換な作用素の和であり、
   固有ベクトルの積によって対角化され、固有値は $`hM+E`$ である。

5. 定理 E.19 の証明により $`I^2=J^2=K^2=-1`$ かつ $`IJ=-JI=K`$、$`JK=-KJ=I`$、$`KI=-IK=J`$ であるから、$`Q_8`$ は群であり、
   $`h\mapsto gh`$ は $`Q_8`$ の全単射である：$`\pi(g)`$ は置換行列であり、実行列によるユニタリ表現であって、
   $`\mathcal A=\pi(Q_8)''`$ は共役で閉じている（補題 4.18(c)）。$`\chi_\pi(g)`$ は $`gh=h`$ を満たす $`h`$ の個数、すなわち $`g=1`$ なら $`8`$、
   それ以外では $`0`$ であるから、系 4.5(b) により既約ユニタリ表現 $`\rho`$ の重複度は
   $`m_\rho=\frac18\overline{\chi_\rho(1)}\cdot8=\dim\rho`$ である。*一次元表現:* 符号を除けば $`I,J,K`$ は $`\mathbb Z_2\times\mathbb Z_2`$ の単位元以外の元と同じように
   掛け合わされるので、$`s,t\in\{\pm1\}`$ に対し $`\chi(\pm1)=1`$、$`\chi(\pm I)=s`$、
   $`\chi(\pm J)=t`$、$`\chi(\pm K)=st`$ は準同型である；これらは四つの相異なる指標であり、したがって四つの互いに非同値な既約
   表現であって、それぞれ $`m=1`$ である。*二次元表現:* $`\rho(g)=M(g)`$（補題 E.15(2)）はユニタリ表現であり
   （$`M(g)^*M(g)=M(\bar gg)=1`$）、$`\chi_\rho(\pm1)=\pm2`$、それ以外では $`\chi_\rho=0`$ である（$`\mathrm{Tr}\,M(a,b)=2\mathrm{Re}\,a`$）。
   よって $`\langle\chi_\rho,\chi_\rho\rangle_{Q_8}=\frac{4+4}8=1`$ であり、$`\rho`$ は既約で（系 4.5(c)）、$`m=2`$ である。定理 2.6 のような
   分解 $`V`$ に対し、$`\pi_c\circ\pi`$ は既約で互いに非同値である（命題 2.7）から、$`V`$ は定理 4.3 のような分解であり、上の五つの表現は
   いずれも $`m_\rho\ge1`$ をもつので、ちょうど一つの $`\pi_c\circ\pi`$ と同値であり、異なる表現は異なる $`c`$ に対応する（系 4.5(b)）。
   $`4\cdot1\cdot1+2\cdot2=8=\sum_ck_cm_c`$ であるから、他のセクターはない。定理 4.20(c) により、これらの表現に対して
   $`\pi_c(\pi(g))=\rho_c(g)`$ となる分解を用いてよい（補題 4.2(c) により $`\rho_c=u\,\pi_c(\pi(\cdot))\,u^*`$ を満たすユニタリ $`u`$ が得られる；
   $`V_c`$ を $`V_c(u^*\otimes1)`$ で置き換える）。$`\pi(g)`$ は実であるから、定理 4.20(a) により
   $`\pi_{\sigma(c)}(\pi(g))=J_c\overline{\rho_c(g)}J_c^*`$ である。指標 $`\chi`$ に対してはこれは $`\chi(g)`$ であるから、$`\pi_{\sigma(c)}\circ\pi\cong\pi_c\circ\pi`$、
   $`\sigma(c)=c`$（命題 2.7：異なるセクターは非同値な表現を担う）であり、$`J_c\bar J_c=|J_c|^2=1`$：実型である。
   $`\rho=M`$ に対しては $`\sigma(c)=c`$（$`k=2`$ のセクターはこれだけである）であり、$`\overline{M(g)}=\Omega_2^*M(g)\Omega_2`$（定理 4.21(b) の証明により $`M(q)\Omega_2=\Omega_2\overline{M(q)}`$）であるから、
   $`M(g)=(J_c\Omega_2^*)M(g)(J_c\Omega_2^*)^*`$ である。行列 $`M(1)=1`$、$`M(I)=\mathrm{diag}(i,-i)`$、$`M(J)=\Omega_2`$、
   $`M(K)=\begin{pmatrix}0&-i\\-i&0\end{pmatrix}`$ は線形独立であるから $`M_2(\mathbb C)`$ を張る；したがって
   $`J_c\Omega_2^*=\lambda1`$、$`J_c=\lambda\Omega_2`$ であり、$`J_c\bar J_c=|\lambda|^2\Omega_2^2=-1`$：四元数型である。
   定理 4.21 により $`\mathcal A_{\mathbb R}\cong\mathbb R^4\oplus M_1(\mathbb H)`$ であり、その次元は $`8=|Q_8|`$ である。

</details>

*計算：* 演習問題 1–5 の数値的な検証 — [チュートリアル、演習問題](../../tutorial/ja/04_symmetries.md#x-4)。

**次へ：**[第5章](05_random_matrices.md) では環の中でランダム行列を抽出する：$`\prod_cU(k_c)`$ のハール・ユニタリ、ワインガルテン計算（その2次の場合は $`\{V\otimes V\}`$ の可換子環、すなわちシューア・ワイル双対性の最小の例に基づく）、そしてセクターごとの準位統計である。
