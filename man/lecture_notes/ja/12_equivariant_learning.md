[English](../12_equivariant_learning.md) · [Español](../es/12_equivariant_learning.md) · [Français](../fr/12_equivariant_learning.md) · [Deutsch](../de/12_equivariant_learning.md) · [中文](../zh/12_equivariant_learning.md) · **日本語** · [Русский](../ru/12_equivariant_learning.md)

# 12. 同変機械学習

ニューラルネットワークの層が **同変** であるとは、入力を対称変換で変換してから層を適用しても、層を適用してから出力を変換しても同じ結果になることをいう：グラフの頂点の番号を付け替えると、その特徴量の番号も同じように付け替わる。線形層についてはこれは可換性の関係であるから、同変な層は可換子環をなし、第2章により、それは基底の取り替えを除いて環 $\bigoplus_cM_{m_c}(\mathbb C)$ である。本章では、それらのパラメータを数え、`torch_vn_algebra.decompose` を用いてそれらを構成・学習し、さらに環の構造を用いて、構成からして同変であるユニタリ層、正値層、スペクトル層を作る。

**本章で学ぶこと**：二つの表現の間の同変線形写像が、重複度で与えられる大きさのブロック行列であるのはなぜか。それらを軌道あるいは指標によってどう数えるか（DeepSets とグラフネットワークのベル数）。群上の平均化が同変写像への直交射影であるのはなぜか。関数計算とユニタリのパラメータ化が可換子環にブロックごとにどう作用するか。そして不変ベクトル、不変二次形式、およびあらゆる次数の不変多項式（モリーン（Molien）の公式）をどう数えるか。

層はライブラリを用いて[チュートリアル 12](../../tutorial/ja/12_equivariant_learning.md)で構成・学習され、同チュートリアルは本章の結果の数値的検証も行う。

本章を通じて、$G$ は有限群であり、$\mathbb C^n$ 上の $G$ の **表現** とはユニタリ値をとる準同型
$\pi:G\to M_n(\mathbb C)$ のことである（これは制約ではない：可逆な $\pi(g)$ は内積
$\sum_g\langle\pi(g)x,\pi(g)y\rangle$ に関してユニタリである。4.1節）。その **指標** は $\chi_\pi(g)=\mathrm{Tr}\,\pi(g)$ である。既約表現、同値な表現、および既約表現 $\rho$ の $\pi$ における **重複度** $m_\rho(\pi)$（$\rho$ が現れなければ $0$）は 2.4節のとおりであり、同節では $m_\rho(\pi)$ が分解によらないこととシューアの補題が証明されている。[第4章](04_symmetries.md) ではその両方を詳しく扱う（4.1–4.2節）。長方形行列の上で $\langle S,T\rangle=\mathrm{Tr}(S^*T)$ は **ヒルベルト・シュミット** 内積である。

<a id="s-12-1"></a>

## 12.1 同変線形写像

**定義.** $\pi_1,\pi_2$ を $\mathbb C^{n_1}$、$\mathbb C^{n_2}$ 上の $G$ の表現とする。行列
$T\in M_{n_2\times n_1}(\mathbb C)$ がすべての $g\in G$ に対して $T\pi_1(g)=\pi_2(g)T$ を満たすとき、$T$ は **同変** である（**絡作用素** である）という。それらはベクトル空間 $\mathrm{Hom}_G(\pi_1,\pi_2)$ をなし、$\pi_1=\pi_2=\pi$ のときこれは可換子環 $\mathrm{Hom}_G(\pi,\pi)=\pi(G)'$ である。

機械学習で最もよく現れる表現は座標を置換するものである。$G$ が有限集合 $\Omega$ に作用するとき、$\mathbb C^\Omega$ 上の **置換表現** は $\pi(g)e_x=e_{gx}$ である。例：集合の $n$ 個の点を置換する $S_n$（$\Omega=[n]=\{1,\dots,n\}$）、およびグラフの頂点の番号を付け替える $S_n$。後者は $\Omega=[n]^2$ に作用し、すなわち隣接行列に $X\mapsto P_\sigma XP_\sigma^T$ によって作用する。

<a id="r-12-1"></a>

**補題 12.1（軌道）.** $\pi_i$ を $\mathbb C^{\Omega_i}$ 上の $G$ の置換表現とする。このとき
$T\in\mathrm{Hom}_G(\pi_1,\pi_2)$ であるための必要十分条件は、すべての $g$、$y\in\Omega_2$、$x\in\Omega_1$ に対して $T_{gy,gx}=T_{yx}$ となることである。$\Omega_2\times\Omega_1$ 上の $G$ の軌道の指示行列は基底をなし、$\dim\mathrm{Hom}_G(\pi_1,\pi_2)$ はこれらの軌道の個数である。

*証明.* $\pi_2(g)T\pi_1(g)^*e_x=\pi_2(g)Te_{g^{-1}x}=\sum_yT_{y,g^{-1}x}e_{gy}$ であるから、$\pi_2(g)T\pi_1(g)^*$ の $(gy,x)$ 成分は $T_{y,g^{-1}x}$ である。よって $T$ が同変であること、すなわちすべての $g$ に対し $\pi_2(g)T\pi_1(g)^*=T$ であることは、すべての $g,y,x$ に対し $T_{gy,x}=T_{y,g^{-1}x}$ であることと同値であり、これは（$x$ を $gx$ で置き換えれば）述べた条件、すなわち $T$ が軌道上で一定であることにほかならない。軌道の指示行列は互いに交わらない台をもつので線形独立であり、軌道上で一定な関数はすべてそれらの線形結合である。$\square$

*計算：* 頂点特徴量と辺特徴量に作用する $S_4$ に対する補題 12.1、軌道の数え上げによる — [チュートリアル、12.1節](../../tutorial/ja/12_equivariant_learning.md#s-12-1)。

<a id="r-12-2"></a>

**系 12.2.** (a)（*DeepSets.*）$\mathbb C^n$ 上の $S_n$（$n\ge2$）に対し、可換子環は $\mathrm{span}\{1,J\}$（$J$ はすべての成分が $1$ の行列）であり、$\mathbb C^n=\mathbb C\mathbf 1\oplus\mathbf 1^\perp$ は互いに非同値な二つの既約表現、すなわち自明表現と次元 $n-1$ の **標準** 表現の和である。(b)（*グラフネットワーク.*）$\mathbb C^{[n]^p}$ から $\mathbb C^{[n]^q}$ への（$n$ 頂点上の $p$ 階および $q$ 階のテンソルの間の）同変写像の空間の次元は、$p+q$ 個の元からなる集合を高々 $n$ 個のブロックに分ける分割の個数に等しい。$n\ge p+q$ のとき、その個数は **ベル数** $B_{p+q}$、すなわち $p+q$ 個の元からなる集合のすべての分割の個数である（$B_2=2$、$B_3=5$、$B_4=15$）。

*証明.* (b) 二つの組 $x,y\in[n]^{p+q}$ が同じ軌道にあるための必要十分条件は、それらが同じ **等号パターン**、すなわち $i\sim j\Leftrightarrow x_i=x_j$ による位置の集合 $\{1,\dots,p+q\}$ の分割をもつことである。パターンが不変であることは明らかである。逆に $x,y$ が同じパターンをもてば、$x_i\mapsto y_i$ は $\{x_i\}$ から $[n]$ への矛盾なく定義された単射であり、これは $\sigma x=y$ を満たす置換 $\sigma$ に拡張される。高々 $n$ 個のブロックをもつ分割はすべてパターンである（各ブロックに相異なる値を与えればよい）。あとは補題 12.1 により結論が従う。$n\ge p+q$ のときはブロック数の制限は意味をもたない。ブロック数ごとに数えると：$B_2=1+1$；$B_3=1+3+1$（$3$ 元集合の二つのブロックへの分割は、その一元ブロックで定まる）；$B_4=1+7+6+1$（二つのブロック：大きさ $3+1$ のものが $4$ 個、$2+2$ のものが $3$ 個；三つのブロック：二元ブロックの選び方が $\binom42=6$ 通り）。(a) $p=q=1$ のとき、$[n]^2$ 上の軌道は対角線とその補集合の二つであり、その指示行列は $1$ と $J-1$ である。よって $\dim\pi(S_n)'=2=\sum_cm_c^2$（定理 2.6、命題 2.7）であり、これは $m_c=1$ の二つのセクター、すなわち互いに非同値な二つの既約な直和成分を強制し、$\pi(S_n)'$ はそれらの射影 $P_1,P_2$ で張られる。不変な直線 $\mathbb C\mathbf 1$ への射影は $\pi(S_n)'$ に属する（補題 2.1）から、それは $P_1$ または $P_2$ である：自明表現を担うこの直線が一つのセクターであり、$\mathbf 1^\perp$ がもう一つのセクターである。$\square$

DeepSets 層 $x\mapsto\lambda x+\gamma\,(\mathbf 1^Tx)\mathbf 1$（Zaheer et al. 2017）は (a) であり、グラフネットワークのための Maron et al. (2019) の基底は (b) である。

*計算：* $n=4$ の系 12.2（$B_2=2$、$B_3=5$、$B_4=15$）、および `dc.from_group` による $S_4$ の頂点特徴量と辺特徴量の分解（辺特徴量：次元 $k=(1,2,3,3)$、重複度 $m=(2,1,1,3)$ のセクター） — [チュートリアル、12.1節](../../tutorial/ja/12_equivariant_learning.md#s-12-1)。

二つの表現の場合も、すべては一つの可換子環に帰着する。

<a id="r-12-3"></a>

**補題 12.3（直和）.** $\pi=\pi_1\oplus\pi_2$ を $\mathbb C^{n_1}\oplus\mathbb C^{n_2}$ 上の表現とし、等長な埋め込みを $Q_i$、射影を $P_i=Q_iQ_i^*$ とする。このとき $P_1,P_2\in\pi(G)'$、
$\mathrm{Hom}_G(\pi_1,\pi_2)=Q_2^*\,\pi(G)'\,Q_1$ であり、$T\mapsto Q_2TQ_1^*$ は
$\mathrm{Hom}_G(\pi_1,\pi_2)$ から $P_2\,\pi(G)'P_1$ の上への線形全単射である。

*証明.* $\pi(g)$ はブロック対角であるから $\pi(g)Q_i=Q_i\pi_i(g)$ かつ $Q_i^*\pi(g)=\pi_i(g)Q_i^*$ であり、したがって
$P_i\pi(g)=Q_i\pi_i(g)Q_i^*=\pi(g)P_i$ である。$X\in\pi(G)'$ ならば
$Q_2^*XQ_1\pi_1(g)=Q_2^*X\pi(g)Q_1=Q_2^*\pi(g)XQ_1=\pi_2(g)Q_2^*XQ_1$ である。$T$ が同変ならば
$\pi(g)Q_2TQ_1^*=Q_2\pi_2(g)TQ_1^*=Q_2T\pi_1(g)Q_1^*=Q_2TQ_1^*\pi(g)$ であるから $Q_2TQ_1^*\in\pi(G)'$ であり、これは
$P_2(Q_2TQ_1^*)P_1$ に等しく、$Q_2^*(Q_2TQ_1^*)Q_1=T$ である。逆に、$X\in\pi(G)'$ による任意の $Y=P_2XP_1$ は
$\pi(G)'$ に属し、$Y=Q_2(Q_2^*YQ_1)Q_1^*$ である。$\square$

*計算：* 頂点特徴量と辺特徴量に対する角ブロック $Q_2^*\,\pi(G)'\,Q_1$ — [チュートリアル、12.2節](../../tutorial/ja/12_equivariant_learning.md#s-12-2)。

<a id="r-12-4"></a>

**定理 12.4（同変写像はブロック行列である）.** $\rho_1,\dots,\rho_C$ を $\pi_1\oplus\pi_2$ に現れる互いに非同値な既約表現とし、$k_c=\dim\rho_c$、$a_c=m_{\rho_c}(\pi_1)$、$b_c=m_{\rho_c}(\pi_2)$ とする。$V_i^*\pi_i(g)V_i=\bigoplus_c\rho_c(g)\otimes1$ を満たすユニタリ
$V_1:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{a_c}\to\mathbb C^{n_1}$ と
$V_2:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{b_c}\to\mathbb C^{n_2}$ が存在して
$$
\mathrm{Hom}_G(\pi_1,\pi_2)=\Big\{V_2\Big(\bigoplus_{c=1}^C1_{k_c}\otimes y_c\Big)V_1^*:\ y_c\in M_{b_c\times a_c}(\mathbb C)\Big\},
$$
となり、$(y_c)_c\mapsto T$ は単射である。したがって $\mathrm{Hom}_G(\pi_1,\pi_2)\cong\bigoplus_cM_{b_c\times a_c}(\mathbb C)$ であり、
$$
\dim\mathrm{Hom}_G(\pi_1,\pi_2)=\sum_\rho m_\rho(\pi_1)\,m_\rho(\pi_2),
$$
である。ここで和は同値を除いたすべての既約表現にわたる。$\pi_1=\pi_2$ のときこれは
$\pi(G)'\cong\bigoplus_cM_{m_c}$（命題 2.7）であり、$\mathrm{Hom}_G(\pi_1,\pi_2)=0$ であるための必要十分条件は、両方に現れる既約表現が存在しないことである。

*証明.* 定理 2.6 と命題 2.7 を $\pi=\pi_1\oplus\pi_2$ に適用する：
$V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$ かつ $V^*\pi(G)'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$ を満たすユニタリ $V$ が存在する。重複度の一意性（2.4節）により、$\pi$ に現れる既約表現はどれもちょうど一つの $\rho_c$ と同値であるから、$\rho_c$ たちは定理の主張における同値類の代表である。他の代表はユニタリ $w_c$ により $\rho_c'=w_c\rho_cw_c^*$ と表され（補題 4.2(c)）、$V$ を $V\bigoplus_c(w_c^*\otimes1)$ で置き換えると、一つ目の標準形は $\bigoplus_c\rho'_c(g)\otimes1_{m_c}$ に変わり、二つ目は変わらない。$w_c^*\otimes1$ は $1\otimes M_{m_c}$ と可換だからである。よって命題 2.7 の $\rho_c$ を用いてよい。
$X\mapsto V^*XV$ は $*$-同型であり、$P_1\in\pi(G)'$（補題 12.3）であるから、射影 $p_c$ により $V^*P_1V=\bigoplus_c1\otimes p_c$ と書ける。$u_c^*p_cu_c=1_{a_c}\oplus0_{b_c}$ となるユニタリ $u_c$ を選び（定理 1.3）、ここではさしあたり
$a_c=\mathrm{rank}\,p_c$、$b_c=m_c-a_c$ とおき、$V$ を $V\bigoplus_c(1\otimes u_c)$ で置き換える。$1\otimes u_c$ は $\rho_c(g)\otimes1$ と可換であるから、これはどちらの標準形も変えない。$\iota_c:\mathbb C^{a_c}\to\mathbb C^{m_c}$ と
$\kappa_c:\mathbb C^{b_c}\to\mathbb C^{m_c}$ をそれぞれ最初の $a_c$ 個の座標と最後の $b_c$ 個の座標の埋め込みとし、
$$
V_1=Q_1^*V\bigoplus_c(1\otimes\iota_c),\qquad V_2=Q_2^*V\bigoplus_c(1\otimes\kappa_c).
$$
とおく。$W_1=V\bigoplus_c(1\otimes\iota_c)$ は値域 $\mathrm{range}\,P_1=\mathrm{range}\,Q_1$ をもつ等長作用素であるから
$V_1^*V_1=W_1^*P_1W_1=1$ であり、$\sum_ck_ca_c=\mathrm{rank}\,P_1=n_1$ であるから $V_1$ はユニタリである。さらに
$\pi_1(g)V_1=Q_1^*\pi(g)W_1=Q_1^*V\big(\bigoplus_c\rho_c(g)\otimes1_{m_c}\big)\bigoplus_c(1\otimes\iota_c)
=V_1\bigoplus_c\rho_c(g)\otimes1_{a_c}$ である。よって互いに非同値な既約表現 $\rho_c$ により $\pi_1\cong\bigoplus_c\rho_c\otimes1_{a_c}$ であり、重複度の一意性（2.4節）により $a_c=m_{\rho_c}(\pi_1)$ である。$V_2$ と $b_c$ についても同様である。最後に、$X=V\big(\bigoplus_c1\otimes x_c\big)V^*\in\pi(G)'$ に対して
$P_2XP_1=V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*$ である。ここで $y_c=\kappa_c^*x_c\iota_c$ は $x_c$ の左下の
$b_c\times a_c$ 隅であり、任意の行列である。補題 12.3 により、
$\mathrm{Hom}_G(\pi_1,\pi_2)=Q_2^*P_2\,\pi(G)'P_1Q_1$ は行列
$Q_2^*V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*Q_1=V_2\big(\bigoplus_c1\otimes y_c\big)V_1^*$ からなり、$V_1,V_2$ はユニタリであるから $(y_c)\mapsto T$ は単射である。$\pi_1\oplus\pi_2$ に現れない既約表現 $\rho$ については $m_\rho(\pi_1)m_\rho(\pi_2)=0$ である。$\square$

*ライブラリでは：* アルゴリズムとしての定理 12.4 の証明、`dc.from_group`、`D.commutant_algebra()`、`D.from_commutant` の上に構成された層 `EquivariantLinear` と `EquivariantMap`。$S_4$ の頂点特徴量から辺特徴量へは $\dim\mathrm{Hom}_G=1\cdot2+1\cdot3=5=B_3$ — [チュートリアル、12.2節](../../tutorial/ja/12_equivariant_learning.md#s-12-2)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. 定理 12.4 は **選択則** の線形代数的な核心であり、ウィグナー・エッカートの定理のスカラーの場合である：対称性と可換な作用素は同じ対称性の型の状態しか結びつけることができず、一つの型の中では $1_{k_c}\otimes y_c$ として、すなわち多重項の中での位置を区別しない「換算行列要素」$y_c$ によって作用する。（ベクトルのように、それ自身が非自明な表現に従って変換する作用素には、さらにクレブシュ–ゴルダン係数が必要である。）対称なハミルトニアンは $\pi_1=\pi_2$ の場合であり、その縮退度 $k_c$ は既約表現の次元である。

<a id="s-12-2"></a>

## 12.2 ライブラリによる同変層

*ライブラリでは：* 定理 12.4 のブロックでパラメータ化された同変層 — [チュートリアル、12.2節](../../tutorial/ja/12_equivariant_learning.md#s-12-2)。

層を同変にする第二の方法は、群上で平均することである。**ツワール** は次で定義される。
$$
\mathcal P(W)=\frac1{|G|}\sum_{g\in G}\pi_2(g)\,W\,\pi_1(g)^*,\qquad W\in M_{n_2\times n_1}(\mathbb C).
$$

<a id="r-12-5"></a>

**定理 12.5（ツワールは直交射影である）.** $\mathcal P$ は $M_{n_2\times n_1}(\mathbb C)$ から $\mathrm{Hom}_G(\pi_1,\pi_2)$ へのヒルベルト・シュミット内積に関する直交射影であり、
$$
\dim\mathrm{Hom}_G(\pi_1,\pi_2)=\frac1{|G|}\sum_{g\in G}\overline{\chi_{\pi_1}(g)}\,\chi_{\pi_2}(g).
$$
である。

*証明.* $h\in G$ に対し $\pi_2(h)\mathcal P(W)\pi_1(h)^*=\frac1{|G|}\sum_g\pi_2(hg)W\pi_1(hg)^*=\mathcal P(W)$ である。
$g\mapsto hg$ は $G$ の全単射だからである。よって $\mathcal P(W)$ は同変である。$W$ が同変ならばすべての項が $W$ に等しいので $\mathcal P(W)=W$ である。したがって $\mathcal P^2=\mathcal P$ であり、その値域は $\mathrm{Hom}_G(\pi_1,\pi_2)$ である。$\mathcal P$ は自己共役である：$\langle S,\mathcal P(W)\rangle=\frac1{|G|}\sum_g\mathrm{Tr}\big((\pi_2(g)^*S\pi_1(g))^*W\big)$
$=\langle\mathcal P(S),W\rangle$。ここで $g\mapsto g^{-1}$ と置換し、$\pi_i(g^{-1})=\pi_i(g)^*$ を用いた。自己共役な冪等作用素はその値域への直交射影である。そのトレースは階数に等しい。$M_{n_2\times n_1}$ 上の線形写像
$W\mapsto AWB$ のトレースは $\mathrm{Tr}\,A\mathrm{Tr}\,B$ である（$AE_{ij}B$ における行列単位 $E_{ij}$ の係数は $A_{ii}B_{jj}$ である）。$A=\pi_2(g)$、$B=\pi_1(g)^*$ とすれば公式が得られる。$\square$

定理 12.4 と合わせると $\sum_\rho m_\rho(\pi_1)m_\rho(\pi_2)=|G|^{-1}\sum_g\overline{\chi_{\pi_1}(g)}\chi_{\pi_2}(g)$ である。
$\pi_1=\pi_2$ のとき、これは 4.2節の $\dim\pi(G)'$ に対する指標公式である。パラメータ化 $W\mapsto\mathcal P(W)$ はすべての同変写像に到達し、$\mathcal P(W)$ は $W$ に最も近い同変写像であり、$\mathcal P$ は自己共役であるから、$W$ に関する勾配は $T=\mathcal P(W)$ に関する勾配の $\mathcal P$ による像である。しかしそのパラメータのうち
$n_1n_2-\dim\mathrm{Hom}_G$ 個が核に浪費される。

*計算：* 辺特徴量に対する定理 12.5、および頂点特徴量から辺特徴量への次元公式 — [チュートリアル、12.2節](../../tutorial/ja/12_equivariant_learning.md#s-12-2)。*ライブラリでは：* ツワール `dc.twirl`、分解から計算した同じ射影 `D.conditional_expectation(W, 'commutant')`、および指標公式 `dc.commutant_dimension` — [チュートリアル、12.2節](../../tutorial/ja/12_equivariant_learning.md#s-12-2)。

以上の層は線形である。ネットワークはそれらを非線形性と交互に重ねるので、非線形性もまた同変でなければならない。置換表現に対しては、各座標に関数を適用する操作は同変である。$(\pi(g)x)_{gy}=x_y$ だからである。他の表現ではそうではなく、ノルムに作用する非線形性や、$G$ が置換する基底の係数に作用する非線形性が用いられる（演習問題 3）。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. ツワールは可換子環への条件付き期待値である（第2章の演習問題 3）：結合をすべての対称操作について平均すると、ちょうどその対称な部分が残る。ハミルトニアンや相互作用を対称群にわたる和によって「対称化」することは、同じ射影である。

<a id="s-12-3"></a>

## 12.3 制約付きの層とスペクトル層

同変層 $x\mapsto Tx$（$T\in\pi(G)'$）の重み $w=(w_c)_c$ は $*$-環 $\bigoplus_cM_{m_c}$ の中にある（定理 12.4）ので、第1章と第11章のすべての構成がそれらに適用できる。$w\in\bigoplus_cM_{m_c}$ に対して $\Phi(w)=V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*$ と書く。

<a id="r-12-6"></a>

**命題 12.6.** (a) $\Phi$ は $\bigoplus_cM_{m_c}$ から $\pi(G)'$ の上への $*$-同型であり、$\Phi(w)$ が自己共役、ユニタリ、あるいは正であるための必要十分条件は、すべての $w_c$ がそうであることである。(b) 正規な $w$ と
$\bigcup_c\mathrm{spec}\,w_c$ 上の関数 $f$ に対し $f(\Phi(w))=\Phi\big((f(w_c))_c\big)$ である。特に、すべての正規な $X\in\pi(G)'$ に対し $f(X)\in\pi(G)'$ である。(c) 任意に固定したユニタリ $u_0\in\bigoplus_cM_{m_c}$ に対し、同変なユニタリは、$s$ が
$\bigoplus_cM_{m_c}$ 全体を動くときの $\Phi(u_0e^{s-s^*})$ にほかならない。特に、それらは歪自己共役な $S\in\pi(G)'$ の指数関数 $e^S$ にほかならない。

*証明.* (a) $(1\otimes x)(1\otimes y)=1\otimes xy$ かつ $(1\otimes x)^*=1\otimes x^*$ であるから、$\Phi$ は線形、単位的、乗法的で $*$ を保つ。また単射であり、定理 2.6 により $\pi(G)'$ の上への写像である。よって $\Phi(w)^*=\Phi(w)$ と
$w^*=w$ は同値であり、$\Phi(w)^*\Phi(w)=1=\Phi(w)\Phi(w)^*$ と $w^*w=1=ww^*$ は同値である。すべての $w_c=t_c^*t_c$ ならば
$\Phi(w)=\Phi(t)^*\Phi(t)\ge0$ である（系 1.7）。逆に $\Phi(w)\ge0$ ならば、$v\in\mathbb C^{m_c}$ とセクター $c$ のベクトル $\xi=V(e_1\otimes v)$ に対して $0\le\langle\xi,\Phi(w)\xi\rangle=\langle v,w_cv\rangle$ である。
(b) $w_c=u_c\Lambda_cu_c^*$ とする（定理 1.3）。このときユニタリ $U=V\bigoplus_c(1\otimes u_c)$ により $\Phi(w)=U\big(\bigoplus_c1\otimes\Lambda_c\big)U^*$ であり、これは $w_c$ の固有値がそれぞれ $k_c$ 回繰り返される対角化である。定義と定理 1.5(a) により $f(\Phi(w))=U\big(\bigoplus_c1\otimes f(\Lambda_c)\big)U^*=\Phi\big((f(w_c))_c\big)$ である。
(c) (a) により、同変なユニタリは $\bigoplus_cM_{m_c}$ のユニタリ $u$ による $\Phi(u)$ であり、補題 11.1 によりそれらの $u$ は
$u_0e^{s-s^*}$ にほかならない。$u_0=1$ とすると (b) により $\Phi(e^{s-s^*})=e^{\Phi(s-s^*)}$ であり、$\Phi(s-s^*)$ は
$\pi(G)'$ の歪自己共役な元全体を動く。$\square$

*ライブラリでは：* `D.commutant_algebra()` 上の `optimize.UnitaryParam` と `PositiveParam` によるユニタリおよび正値の同変層、`apply_function` によるスペクトル層。同変なユニタリ層の学習、および $e^{-tL}$ とその $t$ に関する導関数に対する (b) — [チュートリアル、12.3節](../../tutorial/ja/12_equivariant_learning.md#s-12-3)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. 対称なハミルトニアン $H\in\pi(G)'$ は対称な時間発展
> $e^{-itH}\in\pi(G)'$ を生成する（第8章）。それはセクターを決して混ぜず、重複度空間にのみ作用する。同変なユニタリ層はそのような時間発展であり、与えられるのではなく学習されたものである。

<a id="s-12-4"></a>

## 12.4 不変量

ベクトル $x$ がすべての $g$ に対して $\pi(g)x=x$ を満たすとき **不変** であるといい、$\mathbb C^n$ 上の関数 $q$ が
$q(\pi(g)x)=q(x)$ を満たすとき不変であるという。不変な特徴量は、不変なネットワーク（たとえばグラフ分類器）が出力できるものである。

<a id="r-12-7"></a>

**命題 12.7.** (a) 不変ベクトルは直交射影
$P_0=|G|^{-1}\sum_g\pi(g)$ の値域であり、その次元は $m_{\mathrm{triv}}(\pi)=|G|^{-1}\sum_g\chi_\pi(g)$ である。また自明表現が $\rho_{c_0}$ として現れるならば、$P_0$ は $\pi(G)''$ の極小中心射影 $P_{c_0}$ である。
(b) エルミート形式 $q_A(x)=\langle x,Ax\rangle$（$A=A^*$）が不変であるための必要十分条件は $A\in\pi(G)'$ であり、不変なエルミート形式は次元 $\sum_cm_c^2$ の実ベクトル空間をなす。
(c) 座標 $x_1,\dots,x_n$ に関する $2$ 次の斉次不変多項式は、次元
$$
\frac1{|G|}\sum_{g\in G}\frac{\chi_\pi(g)^2+\chi_\pi(g^2)}2 .
$$
の空間をなす。置換表現に対しては、これは $\Omega$ の二つの元からなる多重集合 $\{x,y\}$ 上の $G$ の軌道の個数である。

*証明.* (a) 不変ベクトルとは $\mathrm{Hom}_G(\mathrm{triv},\pi)$ のことである。ここで $\mathrm{triv}$ は
$\mathbb C=M_{1\times1}$ 上の自明表現である。これに対してツワールは $\mathcal P(x)=P_0x$ であるから、定理 12.5 から射影と次元が、定理 12.4 から $m_{\mathrm{triv}}(\pi)$ が得られる。定理 2.6 の基底で、不変ベクトルを $\xi=\sum_c\xi_c$、
$\xi_c=\sum_jv_{c,j}\otimes e_j\in\mathbb C^{k_c}\otimes\mathbb C^{m_c}$ と書く。$\pi(g)$ はセクターを保つので、すべての $g$ に対し $\rho_c(g)v_{c,j}=v_{c,j}$ である。零でない $v_{c,j}$ は自明表現を担う不変な直線を張るから、そのとき $\rho_c$ は既約であることから自明であり、$c=c_0$ である。よって不変ベクトルはセクター $c_0$ のベクトル、すなわち $\pi(g)$ が $1$ として作用するところのベクトルであり、$P_{c_0}$ の値域をなす。したがって $P_{c_0}$ は $P_0$ である。
(b) $q_A(\pi(g)x)=\langle x,\pi(g)^*A\pi(g)x\rangle$ であり、すべての $x$ に対し $\langle x,Dx\rangle=0$ を満たす行列 $D$ は零である（偏極、系 1.4 の証明）。よって $q_A$ が不変であることは、すべての $g$ に対し $\pi(g)^*A\pi(g)=A$、すなわち $A\in\pi(G)'$ であることと同値である。命題 12.6(a) により、$\pi(G)'$ の自己共役元はエルミートな $w_c$ による $\Phi(w)$ であり、それらは次元 $\sum_cm_c^2$ の実空間をなす。
(c) 二次多項式 $p(x)=x^TBx$ はその対称行列 $B$ で定まる。$(g\cdot p)(x)=p(\pi(g)^{-1}x)$ とおくと
$g\cdot(h\cdot p)=(gh)\cdot p$ であり、$g\cdot p$ の行列は $\sigma(g)(B)=A^TBA$（$A=\pi(g)^{-1}$）である。平均
$M=|G|^{-1}\sum_g\sigma(g)$ は $\sigma(h)M=M$ を満たし、不変な $B$ に対し $MB=B$ を満たすから、（定理 12.5 の証明と同様に）値域が不変多項式全体である冪等写像であり、不変多項式の空間の次元は
$\mathrm{Tr}\,M=|G|^{-1}\sum_g\mathrm{Tr}\,\sigma(g)$ である（冪等写像は値域 $\oplus$ 核に適合した基底で $1\oplus0$ であるから、そのトレースは階数に等しい）。$A$ はユニタリであるから、
$\Lambda=\mathrm{diag}(\lambda)$ により $A=R\Lambda R^{-1}$ と書ける（定理 1.3）。対称行列
$B_{ij}=(R^{-1})^T(e_ie_j^T+e_je_i^T)R^{-1}$（$i\le j$）は基底をなし、
$A^TB_{ij}A=(R^{-1})^T\Lambda(e_ie_j^T+e_je_i^T)\Lambda R^{-1}=\lambda_i\lambda_jB_{ij}$ である。よって
$\mathrm{Tr}\,\sigma(g)=\sum_{i\le j}\lambda_i\lambda_j=\frac12\big((\mathrm{Tr}\,A)^2+\mathrm{Tr}\,A^2\big)
=\frac12\big(\chi_\pi(g^{-1})^2+\chi_\pi(g^{-2})\big)$ であり、$g$ の代わりに $g^{-1}$ について和をとれば公式が得られる。置換表現に対しては $G$ は単項式 $x_ax_b$ を置換するので、多項式が不変であるための必要十分条件は、その係数が単項式の軌道上で一定であることである。$\square$

*計算：* $S_4$ の辺特徴量に対する (a)、(b) および (c) の二通りの数え上げ — [チュートリアル、12.4節](../../tutorial/ja/12_equivariant_learning.md#s-12-4)。

グラフに対しては、不変な辺特徴量は対角線とその補集合の指示ベクトルで張られ（$B_2=2$）、不変な線形読み出しは対角成分の和と非対角成分の和である。頂点特徴量 $\mathbb C^4$ の上では、(c) は二つの頂点からなる多重集合上の $S_4$ の軌道を数える：それは二つある。

$\mathbb C^4$ 上では、それらは $\sum_ix_i^2$ と $\sum_{i<j}x_ix_j$ である。(c) の証明はあらゆる次数で通用する。

<a id="r-12-8"></a>

**命題 12.8（モリーン（Molien）の公式）.** $\mathbb C[x]_d$ を $x_1,\dots,x_n$ の $d$ 次斉次多項式の空間とし、$\mathbb C[x]^G_d$ をその中の不変なもの、すなわちすべての $g$ に対し $p(\pi(g)x)=p(x)$ を満たすもののなす部分空間とする。$|t|<1$ に対し
$$
\sum_{d\ge0}\dim\mathbb C[x]^G_d\,t^d=\frac1{|G|}\sum_{g\in G}\frac1{\det\big(1-t\pi(g)\big)} .
$$
が成り立つ。置換表現に対しては、$\dim\mathbb C[x]^G_d$ は $d$ 次の単項式（$\Omega$ の $d$ 個の元からなる多重集合）上の $G$ の軌道の個数である。$t^2$ の係数は命題 12.7(c) である。

*証明.* 可逆行列 $M$ に対し $\phi_M(p)(x)=p(Mx)$（代入）とおく。$\phi_M$ は線形で $\mathbb C[x]_d$ をそれ自身に写し、$\phi_N\circ\phi_M=\phi_{MN}$ を満たすので、$\phi_M$ は逆写像 $\phi_{M^{-1}}$ をもつ全単射である。よって
$\sigma(g)=\phi_{\pi(g)^{-1}}$ は $\sigma(g)\sigma(h)=\sigma(gh)$ を満たし、$p$ が不変であることは、すべての $g$ に対し $\sigma(g)p=p$ であることと同値である。命題 12.7(c) の証明と同様に、$\mathbb C[x]_d$ 上の $M_d=|G|^{-1}\sum_g\sigma(g)$ は値域が $\mathbb C[x]^G_d$ の冪等写像であるから、$\sigma_d(g)$ を $\sigma(g)$ の $\mathbb C[x]_d$ への制限として $\dim\mathbb C[x]^G_d=\mathrm{Tr}\,M_d=|G|^{-1}\sum_g\mathrm{Tr}\,\sigma_d(g)$ である。$g$ を固定し、ユニタリ $A=\pi(g)^{-1}$ を $A=R\Lambda R^{-1}$、$\Lambda=\mathrm{diag}(\lambda_1,\dots,\lambda_n)$、$|\lambda_i|=1$ と書く（定理 1.3）。一次形式
$\ell_i(x)=(R^{-1}x)_i$ を用いると、積 $\ell^\alpha=\prod_i\ell_i^{\alpha_i}=\phi_{R^{-1}}(x^\alpha)$（$|\alpha|=d$）は $\mathbb C[x]_d$ の基底であり（単項式基底の全単射 $\phi_{R^{-1}}$ による像）、$R^{-1}A=\Lambda R^{-1}$ であるから
$$
(\sigma(g)\ell^\alpha)(x)=\prod_i\big((R^{-1}Ax)_i\big)^{\alpha_i}=\prod_i\big(\lambda_i\ell_i(x)\big)^{\alpha_i}=\lambda^\alpha\ell^\alpha(x).
$$
となる。よって $\mathrm{Tr}\,\sigma_d(g)=\sum_{|\alpha|=d}\lambda^\alpha$ である。$|t|<1$ に対し、$n$ 個の絶対収束級数 $(1-t\lambda_i)^{-1}=\sum_{a\ge0}t^a\lambda_i^a$ の積は展開して次数ごとにまとめ直すことができる：
$$
\frac1{\det(1-tA)}=\prod_{i=1}^n\frac1{1-t\lambda_i}=\sum_\alpha t^{|\alpha|}\lambda^\alpha=\sum_{d\ge0}t^d\mathrm{Tr}\,\sigma_d(g).
$$
$g$ について平均し、$g$ を $g^{-1}$ で置き換える（$A=\pi(g)$ となるように）と公式が得られる。$\dim\mathbb C[x]^G_d\le\dim\mathbb C[x]_d=\binom{n+d-1}d$ であるから両辺は収束する。置換表現に対しては
$(\pi(g)^{-1}x)_y=x_{gy}$ であるから、$\sigma(g)$ は単項式 $\prod_jx_{y_j}$ を単項式 $\prod_jx_{gy_j}$ に写す：$G$ は $d$ 次の単項式を置換し、多項式が不変であるための必要十分条件はその係数が軌道上で一定であることであり、各軌道の単項式の和が $\mathbb C[x]^G_d$ の基底をなす。$\square$

（これは 1897 年のモリーンの定理である。さらに学ぶには B. Sturmfels, *Algorithms in Invariant Theory*, Springer 2008, 2.2節を参照。）$\mathbb C^4$ 上の $S_4$ では、$d\le4$ に対する個数 $1,1,2,3,5$ は $d$ を高々四つの部分に分ける分割の個数である（多重集合の軌道は、その重複度を並べ替えたものによって定まる）。

*計算：* 頂点特徴量と辺特徴量に対し、$t^4$ までのモリーン級数を単項式の軌道と比較 — [チュートリアル、12.4節](../../tutorial/ja/12_equivariant_learning.md#s-12-4)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. ランダウの相転移の理論では、自由エネルギーは実の秩序変数 $x$ の不変多項式であり、$G$ は $x$ に実直交行列として作用する。その二次の部分は $x^TBx$ であり、$B$ は（$\pi(g)^T=\pi(g)^{-1}$ であるから）可換子環に属するので、命題 12.6(a) によりセクターについてブロック対角である。そのうえでランダウの理論は、対称相は $B$ の固有値が最初に符号を変えるセクターで不安定になり、一般的にはそのセクターの既約表現が秩序相の対称性を決定する、と仮定する。

<a id="s-12-5"></a>

## 12.5 比較の対象

*ライブラリでは：* 他のソフトウェア（escnn、e3nn、EMLP）との比較と、ライブラリが付け加えるもの — [チュートリアル、12.5節](../../tutorial/ja/12_equivariant_learning.md#s-12-5)。

<a id="x-12"></a>

## 演習問題

1. （**チャネル付きの DeepSets.**）$S_n$（$n\ge2$）は $\mathbb C^n\otimes\mathbb C^d$ に $P_\sigma\otimes1_d$ によって作用する。$\mathbb C^n\otimes\mathbb C^{d_1}$ から $\mathbb C^n\otimes\mathbb C^{d_2}$ への同変写像は $A,B\in M_{d_2\times d_1}$ による $1\otimes A+J\otimes B$ にほかならないことを示せ。
2. （**実の層.**）$\pi_1,\pi_2$ が実行列をもつとする（たとえば置換表現）。実の同変行列は次元 $\dim_{\mathbb C}\mathrm{Hom}_G(\pi_1,\pi_2)$ の実ベクトル空間をなし、$T\mapsto\mathrm{Re}\,T$ は $\mathrm{Hom}_G(\pi_1,\pi_2)$ をその上へ写すことを示せ。
3. （**非線形性.**）(a) 関数 $f:\mathbb C\to\mathbb C$ を各座標に適用する操作は、すべての置換表現と可換であることを示せ。(b) $f(t)=\max(t,0)$ に対して、それが $\mathbb R^2$ の $\pi/2$ 回転と可換でないことを示せ。(c) $x\mapsto\varphi(\|x\|)x$ はすべてのユニタリ表現と可換であることを示せ。
4. （**層を決定するのに何個のサンプルが必要か？**）$T=\Phi(w)\in\pi(G)'$、$x_1,\dots,x_N\in\mathbb C^n$ とする。セクター $c$ における $V^*x_i$ を $X_{i,c}\in M_{k_c\times m_c}$ により $\sum_{a,j}(X_{i,c})_{aj}e_a\otimes e_j$ と書き、$S_c\in M_{Nk_c\times m_c}$ を行列 $X_{1,c},\dots,X_{N,c}$ を縦に積み重ねたものとする。$Tx_1=\dots=Tx_N=0$ を満たす同変な $T$ は次元 $\sum_cm_c(m_c-\mathrm{rank}\,S_c)$ の空間をなすことを示せ。独立な標準複素ガウス・ベクトル $x_i$ に対しては、これはほとんど確実に $\sum_cm_c\max(0,m_c-Nk_c)$ である。同変な $T$ が $Tx_1,\dots,Tx_N$ によって（同変写像の中で）決定されるための必要十分条件はこの次元が $0$ であることを導き、辺の表現では二つのガウス・サンプルがほとんど確実に同変層を決定し、一つのサンプルでは決して決定されないこと、一方、任意の $16\times16$ 行列は十六個未満のサンプルでは決して決定されないことを導け。

*計算：* 演習問題 1、3、4 の数値的検算 — [チュートリアル、演習問題](../../tutorial/ja/12_equivariant_learning.md#x-12)。

<details>
<summary>解答</summary>

1. $P_\sigma J=J=JP_\sigma$ であるから、これらの写像は同変である。それらは線形独立である：
   $1\otimes A+J\otimes B=0$ ならば、非対角ブロックから $B=0$、次いで $A=0$ が得られる。系 12.2(a) と定理 2.6 により $P_\sigma\otimes1_d\cong(\mathrm{triv}\otimes1_d)\oplus(\mathrm{std}\otimes1_d)$ であるから、二つの既約表現はともに重複度 $d$ をもち、定理 12.4 により $\dim\mathrm{Hom}_G=2d_1d_2$ となる。これは $1\otimes A+J\otimes B$ 全体のなす空間の次元である。

2. $T$ が同変ならば、$\overline{\pi_i(g)}=\pi_i(g)$ であるから $\bar T$ も同変である。よって
   $\mathrm{Re}\,T=(T+\bar T)/2$ と $\mathrm{Im}\,T=(T-\bar T)/2i$ は実かつ同変であり、したがって写像 $T\mapsto\mathrm{Re}\,T$ は実の同変行列全体の上への写像であり（その上では恒等写像である）、$H_{\mathbb R}$ を実の同変行列全体として $\mathrm{Hom}_G(\pi_1,\pi_2)=H_{\mathbb R}\oplus iH_{\mathbb R}$ である。すると $H_{\mathbb R}$ の実基底は $\mathrm{Hom}_G$ の複素基底である：それは全体を張り、実の $T_j$ に関する複素の線形関係 $\sum_j(\alpha_j+i\beta_j)T_j=0$ は実部と虚部に分かれる。補題 12.1 では軌道の指示行列がそのような実基底である。（定理 12.4 のブロック $y_c$ は実のデータに対しても複素である。ライブラリが計算する $V$ は実とは限らないからである。$\mathrm{Re}\,T$ を複素の $y_c$ でパラメータ化すると、必要な数の二倍の実パラメータを用いることになる。）

3. (a) $(\pi(g)f(x))_{gy}=f(x)_y=f(x_y)=f((\pi(g)x)_{gy})$。(b) $R(x_1,x_2)=(-x_2,x_1)$ と $x=(1,-1)$ に対し、
   $f(Rx)=f(1,1)=(1,1)$ であるが $Rf(x)=R(1,0)=(0,1)$ である。(c) $\|\pi(g)x\|=\|x\|$ であるから
   $\varphi(\|\pi(g)x\|)\pi(g)x=\pi(g)\varphi(\|x\|)x$ である。

4. $(1\otimes w_c)(e_a\otimes e_j)=\sum_le_a\otimes(w_c)_{lj}e_l$ であるから、$V^*Tx_i$ はセクター $c$ において係数行列
   $X_{i,c}w_c^T$ をもつ。よってすべての $i$ に対して $Tx_i=0$ であることは、すべての $c$ に対して $S_cw_c^T=0$ であること、すなわち $w_c^T$ のすべての列が次元 $m_c-\mathrm{rank}\,S_c$ の空間 $\ker S_c$ に属することと同値である：この性質をもつ $w_c$ は次元 $m_c(m_c-\mathrm{rank}\,S_c)$ の空間をなし、$\Phi$ は単射である。ガウス分布に従う $x_i$ に対しては、$V^*x_i$ の座標は $V$ の列からなる正規直交基底に関する $x_i$ の座標であり、やはり独立な標準複素ガウス変数であって（命題 A.16(a)）、$i$ についても独立である。よって $S_c$ は $Nk_c\times m_c$ のジニブル行列であり、密度をもつ（命題 A.16(b)）。その階数は高々 $\min(Nk_c,m_c)$ であり、その大きさの左上の正方小行列が特異でない限り等号が成り立つ。その小行列の $|\det|^2$ は成分の実部と虚部の実多項式であって恒等的に零ではない（小行列が単位行列のとき $1$ である）から、系 A.7(b) により階数はほとんど確実に $\min(Nk_c,m_c)$ である。辺の表現では $(k,m)=(1,2),(2,1),(3,1),(3,3)$（計算された分解）であり、寄与するのはセクター $(1,2)$ のみである：$N=1$ では $2\cdot(2-1)=2$、$N=2$ では $0$ である。任意の $16\times16$ 行列は自明群の場合である（セクターは一つで $k=1$、$m=16$）：$\mathrm{rank}\,S\le N$ であるから、$N<16$ 個のサンプルをどう選んでも次元は $16(16-\mathrm{rank}\,S)\ge16(16-N)>0$ であり、$N=16$ 個のガウス・サンプルではほとんど確実に $0$ である。$T\mapsto(Tx_1,\dots,Tx_N)$ は線形であるから、二つの同変写像がサンプル上で一致するための必要十分条件はその差がいま計算した空間に属することであり、これから判定条件が得られる。$N=1$ では、セクター $(1,2)$ において $\mathrm{rank}\,S_c\le Nk_c=1<m_c=2$ であるから、*どの* $x_1$ に対してもこの空間の次元は少なくとも $2$ である。

</details>

**この先へ.** 最終章である[第13章](13_jordan_algebras.md)では、量子力学がジョルダン代数ではなく *-環を用いる理由を論じる。[第4章](04_symmetries.md) は、ここで用いた表現論、すなわちシューアの補題、指標、そして `dc.from_group` の背後にある乱択アルゴリズムを展開する（4.1–4.3節）。[第9章](09_spin_chains.md)は同じセクターをスピン鎖に適用し、[第5章](05_random_matrices.md) の $\prod_cU(m_c)$ 上のハール測度は同変ユニタリ層の自然なランダム初期化を与える。数学については J.-P. Serre, *Linear Representations of Finite Groups* (Springer 1977) と W. Fulton and J. Harris, *Representation Theory* (Springer 1991) を、同変学習については M. Bronstein, J. Bruna, T. Cohen and P. Veličković, *Geometric Deep Learning* (2021) と M. Weiler et al., *Equivariant and Coordinate Independent Convolutional Networks* (2023) を、分解アルゴリズムについては `torch_vn_algebra/decompose.py` の参考文献を参照せよ。
