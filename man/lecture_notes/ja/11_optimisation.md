[English](../11_optimisation.md) · [Español](../es/11_optimisation.md) · [Français](../fr/11_optimisation.md) · [Deutsch](../de/11_optimisation.md) · [中文](../zh/11_optimisation.md) · **日本語** · [Русский](../ru/11_optimisation.md)

# 11. 最適化とトレース不等式

トレース不等式とは、作用素の汎関数の極値についての主張である。3.7節では、もっともらしいが偽であるトレース不等式（命題 3.9）を示した。サンプリングの実行（[チュートリアル、3.7節](../../tutorial/ja/03_traces_and_states.md#s-3-7)）は、数値的にはその破れを見つけられなかった。本章では、高次元でサンプリングが失敗する理由を説明し、それを厳密に
パラメータ化された制約集合上の最適化で置き換える。極値が閉じた形で分かっている場合にはそれを証明し、
最適化の結果を定理と照合できるようにする。

**本章で学ぶこと**：ユニタリと、指定されたコントラストをもつスペクトルをどのようにパラメータ化するか、
そしてなぜすべての許容される作用素に到達できるのか。フォン・ノイマンのトレース不等式とそれが与える厳密な極値。
等号条件つきのトレース不等式。ハール分布のサンプルが平均二乗の意味で極値から遠く離れている理由。そしてマルチスタート探索が何を証明するのか。本章の数値的な検算とライブラリの最適化器は
[チュートリアル 11](../../tutorial/ja/11_optimisation.md)にある。

本章を通じて、$\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$（$D=\sum_ck_c\ge2$）とし、$\mathrm{Tr}$ は
blunt トレース（3.1節）、$\|A\|_1=\sum_c\|A_c\|_1$ は $A$ を
$M_D(\mathbb C)$ のブロック対角行列とみなしたときのトレースノルムとする。$A\in M_k(\mathbb C)$ の特異値を $s_1(A)\ge\dots\ge s_k(A)$ とし、$x^\downarrow$ と
$x^\uparrow$ は実ベクトル $x$ を減少順および増加順に並べ替えたものとする。$\Delta$ はマイケルソン・
コントラスト（3.6節）であり、$\delta\in[0,1)$ に対し
$$
\ell(\delta)=\frac{1-\delta}{1+\delta},\qquad 1-\ell(\delta)=\frac{2\delta}{1+\delta};
$$
とおく。$\lambda_{\max}=1$ の正作用素のコントラストは $\ell(\lambda_{\min})$ であるから、そのコントラストが $\delta$ であるための必要十分条件は
$\lambda_{\min}=\ell(\delta)$ である。$\ell$ は $[0,1]$ の狭義減少な対合だからである。例として解説論文の三つの汎関数
（*From sampling to optimisation* の節）を扱う。

| | 汎関数 | 変数 |
|---|---|---|
| 実験 1 | $z=\lvert\mathrm{Tr}(XUY)\rvert-\mathrm{Tr}(XY)$ | $X,Y\ge0$、$U\in\mathcal M$ はユニタリ |
| 実験 2 | $z=\mathrm{Tr}\lvert XY\rvert-\mathrm{Tr}(\lvert X\rvert\,\lvert Y\rvert)$ | $X=UX_0$、$Y=VY_0$、$X_0,Y_0\ge0$、$U,V$ はユニタリ |
| 実験 3 | $z=\mathrm{Tr}(Y\lvert X\rvert Y)-\mathrm{Tr}\lvert YXY\rvert$ | $X=X^*$、$Y\ge0$ |

ここで $|X|$、$|Y|$ は $\lambda_{\max}=1$ に正規化し、それらのコントラストをパラメータとする。**包絡線** とは
$\overline z(\delta)=\sup\{z:\Delta(|X|)=\delta\}$ と $\underline z(\delta)=\inf\{z:\Delta(|X|)=\delta\}$ のことであり、その他の
変数はすべて自由に動かす。

<a id="s-11-1"></a>

## 11.1 制約集合とそのパラメータ化

**定義.** $\delta\in[0,1)$ と $\ell=\ell(\delta)$ に対し、**コントラスト殻** を次で定める。
$$
\mathcal P_\delta=\{X\in\mathcal M:X\ge0,\ \lambda_{\max}(X)=1,\ \Delta(X)=\delta\}
=\{X=X^*:\mathrm{spec}\,X\subseteq[\ell,1],\ 1,\ell\in\mathrm{spec}\,X\}.
$$
$\mathcal M$ の **スロット** とは $D$ 個の組 $(c,i)$（$i\le k_c$）のことである。$\lambda\in\mathbb R^D$ に対し
$\mathrm{diag}(\lambda)=\bigoplus_c\mathrm{diag}(\lambda_{c,1},\dots,\lambda_{c,k_c})$ とおく。

**パラメータ化.** コントラスト $\delta$ を指定すると、探索は各スタートごとに二つの異なるスロット
$s_\pm=(c_\pm,i_\pm)$（まず $s_+$ を $D$ 個のスロットから一様に、次に $s_-$ を残りの $D-1$ 個から一様に選ぶので、
相異なるスロットの順序対はどれも確率 $1/(D(D-1))$ をもつ）とハール分布に従うランダムなユニタリ $U_0\in\mathcal M$
（補題 11.10）を抽出し、$\theta\in\mathbb R^D$ と
$G=\bigoplus_cG_c$ を次のように写す。
$$
X(\theta,G)=U\mathrm{diag}(\lambda)U^*,\qquad U=U_0\,e^{G-G^*},\qquad
\lambda_{s_+}=1,\quad\lambda_{s_-}=\ell,\quad\lambda_s=\ell+(1-\ell)\tfrac{1+\sin\theta_s}2\ \text{otherwise}.
$$
ユニタリだけなら $G\mapsto U_0e^{G-G^*}$ でパラメータ化される。シグモイドとは異なり、$(1+\sin\theta)/2$ は有限の $\theta$ で $0$ と $1$
に到達する。本章の極値をとるスペクトル（系 11.7、命題 11.14）は $\ell$ と $1$ にあるので、これは都合がよい。
$[a,b]$ 内の自由なコントラストはパラメータが一つ増えるだけであり、$\delta=a+(b-a)(1+\sin\eta)/2$ とすればよい（例えば $[a,b]=[0,0.999]$）。この制限のもとで達成される値は制限なしでも
達成されるから、それは依然として制限なしの上限の下界（制限なしの下限の上界）である。

*ライブラリでは：* このパラメータ化は `optimize.PositiveParam` であり、$G\mapsto U_0e^{G-G^*}$ は `optimize.UnitaryParam` である — [チュートリアル、11.1節](../../tutorial/ja/11_optimisation.md#s-11-1)。

<a id="r-11-1"></a>

**補題 11.1.** 任意のユニタリ $U_0\in\mathcal M$ に対し、$G\mapsto U_0e^{G-G^*}$ は $\mathcal M$ のユニタリ群の上への
写像である。

*証明.* $U$ をユニタリとする。定理 1.3（ブロックごと）と系 1.4(2) により、
$W_c$ をユニタリ、$\varphi_{c,j}\in\mathbb R$ として $U_0^*U=\bigoplus_cW_c\mathrm{diag}(e^{i\varphi_{c,j}})_jW_c^*$ と書ける。このとき
$S=\bigoplus_cW_c\mathrm{diag}(i\varphi_{c,j})_jW_c^*$ は $S^*=-S$ と $e^S=U_0^*U$ を満たし（級数
$\sum_nS^n/n!$ は固有基底において項ごとに計算される）、$G=S/2$ とすれば $G-G^*=S$ である。$\square$

<a id="r-11-2"></a>

**命題 11.2（パラメータ化は全射である）.** $\delta$、$s_\pm$、$U_0$ を固定すると、
$(\theta,G)\mapsto X(\theta,G)$ の像は $\{X\in\mathcal P_\delta:1\in\mathrm{spec}\,X_{c_+},\ \ell\in\mathrm{spec}\,X_{c_-}\}$ であり、
$\mathcal P_\delta$ は相異なるスロットのすべての組にわたるこれらの像の和集合である。

*証明.* 像はこの集合に含まれる：$X$ の固有値は $\lambda_s\in[\ell,1]$ であり、$\lambda_{s_\pm}$ はブロック
$c_\pm$ にある。逆に、$X_c=W_c\mathrm{diag}(\mu_c)W_c^*$ と対角化し、各 $\mu_c$ を並べ替えて（$W_c$ を
$W_c\Pi_c$ で置き換える。$\Pi_c$ は置換行列）$\mu_{s_+}=1$ かつ $\mu_{s_-}=\ell$ となるようにする。これは可能である：$c_+=c_-$
かつ $\delta>0$ ならば、$1\ne\ell$ は同じブロックの二つの固有値であって異なる位置を占める。$\delta=0$ ならば、
すべての固有値は $1$ である。$\ell<1$ ならば、残りのスロットに対して $\theta_s=\arcsin\big(2(\mu_s-\ell)/(1-\ell)-1\big)$ とおけば
$\lambda=\mu$ となる（$\ell=1$ ならば $\theta$ は任意）。補題 11.1 により $U_0e^{G-G^*}=\bigoplus_cW_c$ を満たす $G$ が得られる。
*和集合:* $X\in\mathcal P_\delta$ に対し、固有値 $1$、$\ell$ を含むブロック $c_\pm$ と、その中の相異なるスロットを選ぶ。
$c_+=c_-$ かつ $\delta>0$ ならば、そのブロックは二つの異なる固有値をもつから $k_{c_+}\ge2$ である。$\delta=0$ ならば $X=1$ であり、
任意の二つのスロットでよい。$\square$

*計算：* 命題 11.2 の構成的な検証：証明のパラメータが目標の作用素を再現する — [チュートリアル、11.1節](../../tutorial/ja/11_optimisation.md#s-11-1)。

スロットのどの組も正の確率で抽出されるから、殻のすべての点は正の確率であるスタートの像に含まれ、
（厳密な算術では）制約はすべての反復で成り立つ。
自己共役な $X$（実験 3）に対しては $X=U\mathrm{diag}(\sigma\lambda)U^*$ を用いる。ここで符号 $\sigma\in\{\pm1\}^D$ はスタートごとに抽出され
（各ブロックで、負号の個数は $\{0,\dots,k_c\}$ 上で一様であり、その位置は一様にランダムであるから、どの符号パターンも
正の確率をもつ）、11.5節の符号反転によってのみ変更される。このとき
$|X|=U\mathrm{diag}(\lambda)U^*$ であり（定理 1.5(e)、$\lambda>0$ であるから）、符号とスロットにわたる和集合は
$\{X=X^*:|X|\in\mathcal P_\delta\}$ である：そのような $X$ に対し $X_c=W_c\mathrm{diag}(\nu_c)W_c^*$ と対角化すると、
$|\nu_{c,i}|\ge\ell>0$ として $|X|_c=W_c\mathrm{diag}(|\nu_c|)W_c^*$ となるので、$\mu=|\nu|$ と
$\sigma_s=\nu_s/|\nu_s|$ に命題 11.2 の証明を適用すればよい。意味をもつのは $S=G-G^*$ だけである。

*ライブラリでは：* `optimize.SelfAdjointParam` — [チュートリアル、11.1節](../../tutorial/ja/11_optimisation.md#s-11-1)。

チャート $S\mapsto U_0e^S$ が良い振る舞いをするのは $S=0$ の近くだけである。$\mathfrak u=\{S\in\mathcal M:S^*=-S\}$、
$\mathcal M_h=\{H\in\mathcal M:H^*=H\}$ とおく。これらは実部分空間であり、$\mathcal M=\mathfrak u\oplus\mathcal M_h$ である
（$X=\frac{X-X^*}2+\frac{X+X^*}2$）。*$0$ の近くで、$S\mapsto U_0e^S$ は $\mathfrak u$ の開球 $N\ni0$ からユニタリ群における
$U_0$ の近傍の上への同相写像である。* 実際、$\Psi(S,H)=U_0e^Se^H$ は $C^1$ 写像
$\mathfrak u\times\mathcal M_h\to\mathcal M$ であり（級数 $\sum_{m,n\ge0}U_0S^mH^n/(m!\,n!)$ と、その各項の実座標に関する
偏導関数の級数は有界集合上で一様収束する）、
$e^{tS}e^{tH}=1+t(S+H)+O(t^2)$ により、その $0$ における微分は $(S,H)\mapsto U_0(S+H)$ であって、これは
$\mathcal M$ の上への全単射である。逆関数定理により、$\mathfrak u$ の開球 $N\ni0$ と $\mathcal M_h$ の開球 $N_h\ni0$ が存在して、
$\Psi$ は $N\times N_h$ から開集合 $O\subseteq\mathcal M$ の上への同相写像となる。
$\Psi(S,H)$ がユニタリならば $e^H=e^{-S}U_0^*\Psi(S,H)$ もユニタリである（固有値が純虚数である正規行列 $S$ の固有基底では、
$e^{\pm S}$ は絶対値 $1$ の成分をもつ対角行列であり、互いに逆である）。しかし $e^H$ は正の固有値 $e^{h_i}$（$h_i$ は $H$ の
固有値）をもつから、$e^{h_i}=1$（系 1.4(2)）であり $H=0$ である。
したがって $O\cap U(\mathcal M)=U_0e^N$ である。ここで $U(\mathcal M)$ は $\mathcal M$ のユニタリ群である。これは
$U(\mathcal M)$ における $U_0$ の近傍であり、$S\mapsto U_0e^S$ は連続な全単射 $N\to U_0e^N$ で、その逆写像、すなわち
$\Psi^{-1}$ の第一成分は連続である。$0$ から遠いところではチャートは単射でない：任意の射影 $P\in\mathcal M$ に対し
$e^{2\pi iP}=1=e^0$ である。そのため探索は各ラウンドの後にチャートの **再中心化** を行う（$U_0\leftarrow U_0e^{G-G^*}$、
$G\leftarrow0$）。他の二つのチャートも $0$ において同じ微分 $S\mapsto U_0S$ をもつ。ケーリー・チャート
$S\mapsto U_0C(-S/2)$（演習 1）：$(1+tS/2)(1-tS/2)^{-1}=(1+tS/2)(1+tS/2+O(t^2))=1+tS+O(t^2)$。極分解レトラクション、
すなわち $S$ を $A=U_0(1+S)$ のユニタリ因子に写す写像（定理 1.8）：$A^*A=(1-S)(1+S)=1+S^*S\ge1$ は可逆であるから $|A|$ も
可逆であり、ユニタリ因子は $A|A|^{-1}$ である。$S$ の代わりに $tS$ をとると $|A|^{-1}=(1+t^2S^*S)^{-1/2}=1+O(t^2)$
（$S^*S$ の固有基底で計算する）であるから、この因子は $U_0(1+tS)+O(t^2)$ である。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. マイケルソン・コントラストは干渉縞の可視度 $(I_{\max}-I_{\min})/(I_{\max}+I_{\min})$ である。
> $\lambda_{\max}=1$ をもつ $0\le X\le1$（フィルター、減衰器、検出器の効率）に対しては、
> それはダイナミックレンジ $\lambda_{\min}/\lambda_{\max}$、すなわち条件数の逆数を定める。

<a id="s-11-2"></a>

## 11.2 フォン・ノイマンのトレース不等式

成分が $c_{ij}\ge0$ で、すべての行和と列和が $1$ である行列 $(c_{ij})$ を **二重確率行列** という。ユニタリ $G$ に対し、
$G$ の行と列は単位ベクトルであるから、$(|G_{ij}|^2)$ と $(|G_{ji}|^2)$ は二重確率行列である。

<a id="r-11-3"></a>

**補題 11.3（並べ替え）.** $p_1\ge\dots\ge p_k$、$q_1\ge\dots\ge q_k$ を実数とし、$(c_{ij})$ を二重確率行列とする。
このとき $\sum_ip_iq_{k+1-i}\le\sum_{i,j}p_iq_jc_{ij}\le\sum_ip_iq_i$ である。

*証明.* *上界.* $m\le k$ に対し $r_j=\sum_{i\le m}c_{ij}\in[0,1]$ とおくと、$\sum_jr_j=m$ である。
$j\le m$ に対し $q_j\ge q_m$、$j>m$ に対し $q_j\le q_m$ であるから、
$$
\sum_{j\le m}q_j-\sum_jq_jr_j=\sum_{j\le m}q_j(1-r_j)-\sum_{j>m}q_jr_j\ge q_m\Big(m-\sum_jr_j\Big)=0 .
$$
$x_i=\sum_jq_jc_{ij}$ とおくと、これは $\sum_{i\le m}x_i\le\sum_{j\le m}q_j$ を意味し、$m=k$ では等号が成り立つ。
アーベルの総和法と $p_m-p_{m+1}\ge0$ により、
$$
\sum_ip_ix_i=p_k\sum_{i\le k}x_i+\sum_{m<k}(p_m-p_{m+1})\sum_{i\le m}x_i
\le p_k\sum_{j\le k}q_j+\sum_{m<k}(p_m-p_{m+1})\sum_{j\le m}q_j=\sum_ip_iq_i .
$$
*下界.* 上界を $p$、減少列 $q'_j=-q_{k+1-j}$、および二重確率行列
$c'_{ij}=c_{i,k+1-j}$ に適用する：$-\sum_{ij}p_iq_jc_{ij}=\sum_{ij}p_iq'_jc'_{ij}\le\sum_ip_iq'_i=-\sum_ip_iq_{k+1-i}$。$\square$

任意の $A\in M_k(\mathbb C)$ は **特異値分解** $A=W\Sigma V^*$（$W,V$ はユニタリ、
$\Sigma=\mathrm{diag}(s_i(A))$）をもつ：$A=W'|A|$（定理 1.8）、$|A|=V\Sigma V^*$（定理 1.3）と書き、$W=W'V$ とすればよい。

<a id="r-11-4"></a>

**定理 11.4（フォン・ノイマンのトレース不等式）.** $A,B\in M_k(\mathbb C)$ に対し
$|\mathrm{Tr}(AB)|\le\sum_is_i(A)\,s_i(B)$ である。

*証明.* $A=W_A\Sigma_AV_A^*$、$B=W_B\Sigma_BV_B^*$、$a_i=s_i(A)$、$b_j=s_j(B)$、$G=V_A^*W_B$、$H=V_B^*W_A$ とおく。
巡回性により $\mathrm{Tr}(AB)=\mathrm{Tr}(\Sigma_AG\Sigma_BH)=\sum_{i,j}a_iG_{ij}b_jH_{ji}$ であるから、
$|g||h|\le\frac12(|g|^2+|h|^2)$ と $a_ib_j\ge0$ により
$$
|\mathrm{Tr}(AB)|\le\sum_{i,j}a_ib_j\,d_{ij},\qquad d_{ij}=\tfrac12\big(|G_{ij}|^2+|H_{ji}|^2\big).
$$
$(d_{ij})$ は二重確率行列であり、補題 11.3 により $\sum_{ij}a_ib_jd_{ij}\le\sum_ia_ib_i$ である。$\square$

<a id="r-11-5"></a>

**系 11.5.** (a) 固有値 $p,q$ をもつ自己共役な $P,Q\in M_k(\mathbb C)$ に対し、
$\sum_ip^\downarrow_iq^\uparrow_i\le\mathrm{Tr}(PQ)\le\sum_ip^\downarrow_iq^\downarrow_i$。
(b) $P,Q\ge0$ に対し、$\mathrm{Tr}(PQ)\le\|PQ\|_1\le\sum_ip^\downarrow_iq^\downarrow_i$。

*証明.* (a) 定理 3.5 の証明と同様に、減少順に並べた固有基底に対する二重確率行列
$c_{ij}=|\langle u_i,v_j\rangle|^2$ を用いて $\mathrm{Tr}(PQ)=\sum_{ij}p_iq_jc_{ij}$ と書ける。補題 11.3 を適用する（その下界は
$\sum_ip^\downarrow_iq^\downarrow_{k+1-i}=\sum_ip^\downarrow_iq^\uparrow_i$ である）。(b) 第一の不等式は
定理 1.9 である。定理 1.9 と定理 11.4 により $\|PQ\|_1=\max_U|\mathrm{Tr}(P\,QU)|\le\sum_is_i(P)s_i(QU)$ であり、
$|QU|=(U^*Q^2U)^{1/2}=U^*QU$（定理 1.6）であるから $s_i(QU)=q^\downarrow_i$ である。$\square$

*計算：* ランダム行列での定理 11.4 と系 11.5 — [チュートリアル、11.2節](../../tutorial/ja/11_optimisation.md#s-11-2)。

したがって、スペクトル $x,y$ をもつ $M_k(\mathbb C)$ の $X,Y\ge0$ に対し、命題 3.9 の破れ $\|YX\|_1-\mathrm{Tr}(XY)$ は
高々 $\sum_ix^\downarrow_i(y^\downarrow_i-y^\uparrow_i)$ である（$(Y,X)$ に対する系 11.5(b) と (a) の下界による）。
実験 2 については、この数がちょうど厳密な答えである。

<a id="r-11-6"></a>

**命題 11.6（スペクトルを固定した実験 2）.** 実験 2 では $|X|=X_0$、$|Y|=Y_0$ であり、
$z=\|X_0VY_0\|_1-\mathrm{Tr}(X_0Y_0)$ である。$(X_0)_c$ のスペクトル $x_c$ と $(Y_0)_c$ のスペクトル $y_c$ を固定し、
固有基底と $V$ を動かすとき、
$$
\sup z=-\inf z=\sum_c\sum_ix^\downarrow_{c,i}\big(y^\downarrow_{c,i}-y^\uparrow_{c,i}\big),
$$
であり、どちらも達成される。

*証明.* $X^*X=X_0U^*UX_0=X_0^2$ であるから $|X|=X_0$（定理 1.6）。同様に $|Y|=Y_0$ であり、$|UA|=|A|$ により
$\mathrm{Tr}|XY|=\|X_0VY_0\|_1$ である。すべての変数はブロックごとに作用し、$z$ はブロックにわたる和であるから、$C=1$ としてよい。
$Y'=VY_0V^*\ge0$（スペクトル $y$）とおく。$|AV|=V^*|A|V$ であるから $\|X_0VY_0\|_1=\|X_0Y'V\|_1=\|X_0Y'\|_1$ である。系 11.5
を組 $(X_0,Y')$ と $(X_0,Y_0)$ に適用し、$\|X_0Y'\|_1\ge\mathrm{Tr}(X_0Y')$ を用いると、
$$
\sum_ix^\downarrow_iy^\uparrow_i-\sum_ix^\downarrow_iy^\downarrow_i\le\|X_0Y'\|_1-\mathrm{Tr}(X_0Y_0)\le\sum_ix^\downarrow_iy^\downarrow_i-\sum_ix^\downarrow_iy^\uparrow_i .
$$
を得る。$J$ を逆順の置換 $Je_i=e_{k+1-i}$ とする。このとき $J\mathrm{diag}(w)=\mathrm{diag}(w_k,\dots,w_1)J$ であり、
$\mathrm{diag}(w)J$ のトレースノルムは $\sum_i|w_i|$ である。$X_0=\mathrm{diag}(x^\downarrow)$、$V=J$ とすると、
$Y_0=\mathrm{diag}(y^\uparrow)$ に対しては $X_0VY_0=\mathrm{diag}(x^\downarrow_iy^\downarrow_i)J$ かつ
$\mathrm{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\uparrow_i$（上界）であり、$Y_0=\mathrm{diag}(y^\downarrow)$ に対しては
$X_0VY_0=\mathrm{diag}(x^\downarrow_iy^\uparrow_i)J$ かつ $\mathrm{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\downarrow_i$
（下界）である。$\square$

<a id="r-11-7"></a>

**系 11.7（実験 2 の包絡線）.** $\ell_X,\ell_Y\in(0,1]$ とし、ある $k_c\ge2$ とする。
$\lambda_{\max}=1$、$\lambda_{\min}(X_0)=\ell_X$、$\lambda_{\min}(Y_0)=\ell_Y$ を満たすすべての $X_0,Y_0\ge0$ とすべてのユニタリにわたって、
$$
\sup z=-\inf z=\sum_c\Big\lfloor\frac{k_c}2\Big\rfloor(1-\ell_X)(1-\ell_Y).
$$

*証明.* 一つのブロックにおいて $x_i=x^\downarrow_{c,i}$、$y_i=y^\downarrow_{c,i}$、$k=k_c$ と書くと、$y^\uparrow_{c,i}=y_{k+1-i}$ である。
$i$ と $k+1-i$ を対にすると $x_i(y_i-y_{k+1-i})+x_{k+1-i}(y_{k+1-i}-y_i)=(x_i-x_{k+1-i})(y_i-y_{k+1-i})$ であり、
$k$ が奇数のとき中央の添字は $x_i(y_i-y_i)=0$ を寄与する。したがって
$$
\sum_ix_i(y_i-y_{k+1-i})=\sum_{i\le\lfloor k/2\rfloor}(x_i-x_{k+1-i})(y_i-y_{k+1-i})\le\Big\lfloor\frac k2\Big\rfloor(1-\ell_X)(1-\ell_Y),
$$
である。$x_i\in[\ell_X,1]$、$y_i\in[\ell_Y,1]$ は減少列であるから、二つの因子はそれぞれ $[0,1-\ell_X]$ と $[0,1-\ell_Y]$ に
含まれるからである。すべてのブロックにおいて、大きい方から $\lfloor k_c/2\rfloor$ 個の
固有値が $1$ で、小さい方から $\lfloor k_c/2\rfloor$ 個が $\ell_X$（それぞれ $\ell_Y$）であり、中央のもの（$k_c$ が奇数のとき）は
$[\ell_X,1]$（それぞれ $[\ell_Y,1]$）内で任意であるとき、等号が成り立つ。ある $k_c\ge2$ であるから、これらのスペクトルの全体での極値は $1$ と $\ell_X$（$\ell_Y$）である。命題 11.6
により結論を得る。$\square$

ここで $\Delta(|X|)=\delta$ に対し $1-\ell_X=2\delta/(1+\delta)$ である。$k_c$ が奇数のとき中央の固有値は関係しない。

*計算：* $k=3$、$C=2$ に対する命題 11.6 の極値配置と系 11.7 の包絡線 — [チュートリアル、11.2節](../../tutorial/ja/11_optimisation.md#s-11-2)；マルチスタート探索で求めた包絡線 — [チュートリアル、11.5節](../../tutorial/ja/11_optimisation.md#s-11-5)。

<a id="s-11-3"></a>

## 11.3 等号条件つきのトレース不等式

$X=X^*$ に対し、$f_+(t)=\max(t,0)$、$f_-(t)=\max(-t,0)$ として $X_\pm=f_\pm(X)$ とおく。このとき $X_\pm\ge0$、$X=X_+-X_-$、
$|X|=X_++X_-$、$X_+X_-=0$ である（定理 1.5(c)。$|X|=(X^2)^{1/2}=h(X)$、
$h(t)=\sqrt{t^2}=|t|=f_+(t)+f_-(t)$ については 1.5(e)。正値性は系 1.4(3)）。

<a id="r-11-8"></a>

**補題 11.8.** $M_n(\mathbb C)$ の $A,B\ge0$ に対し $\|A-B\|_1\le\mathrm{Tr}\,A+\mathrm{Tr}\,B$ であり、等号は
$AB=0$ のときかつそのときに限り成り立つ。

*証明.* $H=A-B$ とし、$(0,\infty)$ 上で $g=1$、それ以外で $g=0$ として $Q=g(H)$ とおく。$g=\bar g=g^2$ であるから $Q$ は射影である（定理 1.5(c)）。$|t|=t(2g(t)-1)$ であるから
$|H|=H(2Q-1)$ であり（定理 1.5(c),(e)）、射影 $R=Q,1-Q$ に対して $\mathrm{Tr}(AR)=\mathrm{Tr}(RAR)$（巡回性）を
用いると、
$$
\|H\|_1=\mathrm{Tr}\big(A(2Q-1)\big)-\mathrm{Tr}\big(B(2Q-1)\big)
=\mathrm{Tr}\,A+\mathrm{Tr}\,B-2\mathrm{Tr}\big((1-Q)A(1-Q)\big)-2\mathrm{Tr}(QBQ).
$$
である。引かれているトレースは $\ge0$ であり、$\mathrm{Tr}(RAR)=\|A^{1/2}R\|_2^2$ が零となるのは $AR=0$ のときかつそのときに限る。したがって等号が成り立つのは
$A(1-Q)=0$ かつ $BQ=0$ のときかつそのときに限り、そのとき $AB=AQB=A(BQ)^*=0$ である。逆に $AB=0$ ならば $BA=0$ であり、
$(A-B)^2=A^2+B^2=(A+B)^2$ であるから、定理 1.6 により $|A-B|=A+B$ である。$\square$

<a id="r-11-9"></a>

**命題 11.9（実験 3）.** $\mathcal M$ の $X=X^*$ と $Y\ge0$ に対し
$\mathrm{Tr}|YXY|\le\mathrm{Tr}(Y|X|Y)$ であり、等号は $YX_+Y^2X_-Y=0$ のときかつそのときに限り成り立つ。特に $XY=YX$ ならば等号が成り立つ。

*証明.* $M_D(\mathbb C)$ において $A=YX_+Y$ と $B=YX_-Y$ は正であり（系 1.7 を $T=X_\pm^{1/2}Y$ として適用）、$YXY=A-B$
かつ $\mathrm{Tr}(Y|X|Y)=\mathrm{Tr}\,A+\mathrm{Tr}\,B$ である。補題 11.8 を適用する。$XY=YX$ ならば $Y$ は
$X_\pm$ と可換であり（定理 1.5(d)）、$AB=Y^4X_+X_-=0$ である。$\square$

したがって実験 3 では $z\ge0$ であり、$\underline z=0$ である（$Y=1$ とする）。等号条件は可換性より弱い：
$X=\mathrm{diag}(1,-1,0)$ と、行が $(1,0,a)$、$(0,1,0)$、$(a,0,1)$（$|a|<1$）である $Y$ に対し、
$AB=(Ye_1)\langle Ye_1,Ye_2\rangle(Ye_2)^*=0$ であるが、$a\ne0$ ならば $XY\ne YX$ である。

*計算：* ランダムな組での命題 11.9 と、可換でないのに等号が成り立つ例 — [チュートリアル、11.3節](../../tutorial/ja/11_optimisation.md#s-11-3)。

<a id="s-11-4"></a>

## 11.4 サンプリング対最適化

$|\mathrm{Tr}(AU)|$ は、その最大値 $\|A\|_1$ と比べて典型的にどれほどの大きさか。我々は $U$ を、
**ジニブル行列** $Z$（成分は独立で密度
$\pi^{-1}e^{-|z|^2}$ をもつ。Mezzadri の処方）の QR 分解 $Z=QR$（$R_{ii}>0$）の因子 $Q$ として、各ブロックで独立にサンプリングする。

<a id="r-11-10"></a>

**補題 11.10.** $Q$ は $U(k)$ 上でハール分布に従う。特に、任意の固定されたユニタリ $V\in M_k(\mathbb C)$ に対し、$VQ$ は
$Q$ と同じ分布をもつ。独立なブロックをもつユニタリ $\bigoplus_cQ_c$ は、$\mathcal M$ のユニタリ群上でハール分布に従う。

*証明.* これは定理 5.3（左不変性はその証明の同変性の段階である。分布は一意性の補題 5.1 によってハール測度と
同定され、その測度論的な段階は系 A.2 である）と系 5.4 である。$\square$

*ライブラリでは：* このサンプラーは `TypeIAlgebra.random_unitary_operator` である — [チュートリアル、11.4節](../../tutorial/ja/11_optimisation.md#s-11-4)。

<a id="r-11-11"></a>

**命題 11.11.** ハール分布に従う $U\in U(k)$ と $A\in M_k(\mathbb C)$ に対し、$\mathbb E\mathrm{Tr}(AU)=0$ かつ
$\mathbb E|\mathrm{Tr}(AU)|^2=\|A\|_2^2/k$ である。

*証明.* 系 5.5 により $\mathbb EU_{ij}=0$ である。命題 5.7（$C=1$）を $X=|e_j\rangle\langle e_m|$ に適用すると、
1次のワインガルテン公式
$\mathbb E[U_{ij}\bar U_{lm}]=\mathbb E(UXU^*)_{il}=\delta_{il}\delta_{jm}/k$ を得る。$k\ge2$ に対しては、これは
2次の公式（定理 5.10）において $i_2=i'_2$、$j_2=j'_2$ とおき両者について和をとることからも従う。
$\sum_{i_2,j_2}|U_{i_2j_2}|^2=k$ だからである。$\mathrm{Tr}(AU)=\sum_{ij}A_{ji}U_{ij}$ を用いると、
$\mathbb E|\mathrm{Tr}(AU)|^2=\sum_{ijlm}A_{ji}\bar A_{ml}\mathbb E[U_{ij}\bar U_{lm}]=\frac1k\sum_{ij}|A_{ji}|^2$ である。$\square$

<a id="r-11-12"></a>

**系 11.12.** $M_k(\mathbb C)$ の $X,Y\ge0$ のスペクトルが $[a,b]$（$a>0$）に含まれるとし、$U$ をハール分布に従うとする。このとき
$$
\frac{\big(\mathbb E|\mathrm{Tr}(XUY)|^2\big)^{1/2}}{\max_V|\mathrm{Tr}(XVY)|}=\frac{\|YX\|_2}{\sqrt k\,\|YX\|_1}\le\min\Big(\frac1{\sqrt k},\frac{b^2}{a^2k}\Big),
\qquad
\mathbb P\big(|\mathrm{Tr}(XUY)|\ge\mathrm{Tr}(XY)\big)\le\frac{b^4}{a^4k^2}.
$$

*証明.* 等式は $A=YX$ に対する命題 11.11（$\mathrm{Tr}(XUY)=\mathrm{Tr}(AU)$ であるから）と命題 3.9 である。
$\|A\|_2^2=\sum_is_i^2\le(\sum_is_i)^2=\|A\|_1^2$（1.4節）であるから、比は $\le1/\sqrt k$ である。
$\|Av\|=\||A|v\|$（1.4節）であるから、$|A|$ を対角化すると $s_k(A)=\min_{\|v\|=1}\|Av\|$ および
$s_1(A)=\max_{\|v\|=1}\|Av\|$ を得る。また単位ベクトル $v$ に対し $\|Xv\|^2=\langle v,X^2v\rangle\in[a^2,b^2]$ であり
（$X^2-a^2$ と $b^2-X^2$ に対する系 1.4(3)）、$Y$ についても同様である。よって $a^2\le s_k(A)\le s_1(A)\le b^2$、$\|A\|_2^2\le kb^4$、
$\|A\|_1\ge ka^2$ である。最後に $\mathrm{Tr}(XY)=\mathrm{Tr}(X^{1/2}YX^{1/2})\ge a\mathrm{Tr}\,X\ge ka^2$ であり、
$|\mathrm{Tr}(XUY)|^2$ に対するマルコフの不等式により、確率は $(kb^4/k)/(ka^2)^2$ で抑えられる。$\square$

したがって、$N$ 個のハール・サンプルが $|\mathrm{Tr}(XUY)|\le\mathrm{Tr}(XY)$ の破れを見つける確率は高々
$Nb^4/(a^4k^2)$ である（$N$ 個の事象の和）。破れは可換でないすべての組に対して存在する（命題 3.9）にもかかわらず、である。
この上界が意味をもつのは $k>\sqrt N\,b^2/a^2$ のときだけである（$b/a=6$、$N=2000$ では $k>1609$）。

*計算：* モンテカルロ法による命題 11.11 の検証と、$k\le32$ に対する系 11.12 の比と $2000$ 個のサンプルの最良値との比較 — [チュートリアル、11.4節](../../tutorial/ja/11_optimisation.md#s-11-4)。

<a id="r-11-13"></a>

**補題 11.13.** $k\ge2$ ならば、任意の $A\in M_k(\mathbb C)$ に対し $\min_U|\mathrm{Tr}(AU)|=0$ である。

*証明.* $A=W\Sigma V^*$ とし、$P$ を巡回シフト $Pe_i=e_{(i\bmod k)+1}$ とする。$k\ge2$ ならばその対角成分は零である。
$U=VPW^*$ に対し $\mathrm{Tr}(AU)=\mathrm{Tr}(\Sigma P)=0$ である。$\square$

したがって $\inf_Uz=-\mathrm{Tr}(XY)$、$\sup_Uz=\|YX\|_1-\mathrm{Tr}(XY)$ である（命題 3.9）。命題 11.11 により $|\mathrm{Tr}(XUY)|$ の二乗平均平方根は $\|YX\|_2/\sqrt k$ である。
系 11.12 の証明における評価（スペクトルは $[a,b]$ 内）により $\mathbb E|\mathrm{Tr}(XUY)|^2\le b^4$ かつ
$\mathrm{Tr}(XY)\ge ka^2$ であるから、マルコフの不等式により $\varepsilon>0$ に対し
$\mathbb P\big(z\ge-(1-\varepsilon)\mathrm{Tr}(XY)\big)\le b^4/(\varepsilon^2a^4k^2)$ が成り立つ：$k$ が大きいとき、
$z$ は高い確率でその下限に近い。

*計算：* 勾配上昇法と勾配降下法による $U$ についての両極値と、$20\,000$ 個のハール・サンプルとの比較 — [チュートリアル、11.4節](../../tutorial/ja/11_optimisation.md#s-11-4)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. 系 11.12 は、統計力学で **典型性** と呼ばれるものの、
> 弱いが証明可能な一例である：高次元では、ほとんどすべての状態や時間発展が巨視的な量にほぼ同じ値を与えると
> 期待される。鋭い不等式は、微調整された配置——揃った固有基底、2準位のスペクトル——によって決まり、
> ランダムなユニタリ（ランダム回路、ランダムな散乱行列）はそうした配置をめったに生み出さない。

<a id="s-11-5"></a>

## 11.5 包絡線のためのマルチスタート探索

**マルチスタート探索** は、それぞれ固有のスロット、符号、$U_0$、コントラストをもつ多数の独立なスタートから、
勾配法をラウンドに分けて実行する。各ラウンドの後にユニタリ・チャートを再中心化し、スタートごとに、反復点における
目的関数の最良値を返す。自己共役なパラメータ化の符号は離散的である：ラウンドの間に、ランダムなスロットでの単一の
符号反転によって変更され、各反転は目的関数を改善する場合に保持される。

*ライブラリでは：* `optimize.extremize` — [チュートリアル、11.5節](../../tutorial/ja/11_optimisation.md#s-11-5)。

厳密な算術では、すべての反復は許容的である（命題 11.2 とその後の注意）から、返される値はすべて許容点における
$z$ の値である。浮動小数点ではそれは丸め誤差を除いてその値である。したがって最大化は上限の
**下界** を、最小化は下限の上界を与え、それ以上のものではない——局所探索は局所極値で止まりうる。
逆向きの評価には、系 11.7 や命題 11.9 のような証明か、凸緩和（11.6節）が必要である。計算された包絡線は
(i) 閉じた形と照合し、(ii) 小さな解の直和のような手で構成した配置と照合し（包絡線はそれを下回ってはならない）、
(iii) 構造化されたサンプル（2準位のスペクトル、置換基底）と照合する。後者の値もまた達成される値である。

<a id="r-11-14"></a>

**命題 11.14（実験 1 の下側包絡線）.** すべての $k_c\ge2$ とする。$X\in\mathcal P_\delta$、
$\lambda_{\max}(Y)=1$ を満たす $Y\ge0$、およびユニタリ $U\in\mathcal M$ にわたって $\ \min z=-\big(D-2\delta/(1+\delta)\big)$ である。

*証明.* $Y\le1$ から $X^{1/2}YX^{1/2}\le X$ が従うので、$\mathrm{Tr}(XY)\le\mathrm{Tr}\,X\le(D-1)+\ell(\delta)$ である。
$X$ の一つの固有値は $\ell(\delta)$ であり、他は $\le1$ だからである。よって $z\ge-\mathrm{Tr}(XY)\ge-(D-2\delta/(1+\delta))$ である。
等号は、$X$ の固有値が一度だけ $\ell(\delta)$ でそれ以外は $1$、$Y=1$、かつ $U=\bigoplus_cU_c$ が
$\mathrm{Tr}(X_cU_c)=0$ を満たす（補題 11.13）とき成り立つ。$\square$

型 (ii) の検証には直和を使うことができる：実験 3 の $z$ は直和について加法的であり、直和のコントラストは各成分のコントラストの
最大値である（演習 3）から、$\Delta(Y)$ を固定された区間 $[0,b]$ 内で自由とした実験 3 の包絡線について、
$$
\overline z_{M_2\oplus M_2}(\delta)\ \ge\ \overline z_{M_2}(\delta)+\sup_{\delta'\le\delta}\overline z_{M_2}(\delta');
$$
特に、$M_2$ 上でコントラスト $\delta$ と $\delta'\le\delta$ において達成される二つの値の和は、$M_2\oplus M_2$ 上でコントラスト
$\delta$ において達成される値である。

*計算：* マルチスタート探索による $M_2\oplus M_2$ 上の命題 11.14、$M_3\oplus M_3$ 上の系 11.7、命題 11.9 の検証と、型 (ii) の検証 — [チュートリアル、11.5節](../../tutorial/ja/11_optimisation.md#s-11-5)。

<a id="s-11-6"></a>

## 11.6 比較の対象

*ライブラリでは：* 他のソフトウェアとの比較。その中には、凸緩和によって上限を上から抑える SDP ソルバーも含まれる — [チュートリアル、11.6節](../../tutorial/ja/11_optimisation.md#s-11-6)。

<a id="x-11"></a>

## 演習問題

1. （**ケーリー変換.**）$S^*=-S$ に対し、$C(S)=(1-S)(1+S)^{-1}$ がユニタリであり、$S\mapsto C(S)$ が
   固有値 $-1$ をもたないユニタリ全体の上への全単射であることを示せ。したがって指数写像（補題 11.1）とは異なり、ケーリー・チャート
   $U_0C(S)$ は $-U_0$ を取りこぼす。
2. （**実ハール測度.**）$O(k)$ 上でハール分布に従う $O$（実ジニブル行列の QR）と実行列 $A$ に対し、
   $\mathbb E\mathrm{Tr}(AO)^2=\|A\|_2^2/k$ を示し、実環で確かめよ。
3. （**直和.**）実験 3 の $z$ が直和について加法的であること、および $\lambda_{\max}=1$ をもつ $A,A'\ge0$ に対し
   $\Delta(A\oplus A')=\max(\Delta(A),\Delta(A'))$ であることを示せ。11.5節の直和に関する
   不等式を導け（$\Delta(Y)$ は固定された区間 $[0,b]$ 内で自由とする）。

<details>
<summary>解答</summary>

1. $S$ は正規で $iS$ は自己共役であるから、$S$ の固有値は純虚数であり（系 1.4(2)）、$1\pm S$ は可逆である。
   これらは可換であり、$(1+S)^*=1-S$ である。したがって
   $C^*C=(1-S)^{-1}(1+S)(1-S)(1+S)^{-1}=1$ であり、$1+C=2(1+S)^{-1}$ は可逆であるから $-1\notin\mathrm{spec}\,C$ である。
   逆に、$1+U$ が可逆であるユニタリ $U$ に対し、$S=(1-U)(1+U)^{-1}$ は歪自己共役である：
   $S^*=(1+U^*)^{-1}(1-U^*)=(U+1)^{-1}U\,U^*(U-1)=(U+1)^{-1}(U-1)=-S$（$U\pm1$ は可換であるから）。可換な $C$ と $S$ に対し、
   $C(1+S)=1-S$ は $S(1+C)=1-C$ と同値である（どちらも $C+CS+S=1$ を意味する）。さて $C(S)$ は $S$ と可換であり、
   $S=(1-U)(1+U)^{-1}$ は $U$ と可換である。よって $U=C(S)$ と $S=(1-U)(1+U)^{-1}$ は同値である：二つの写像は互いに逆である。

2. 定理 B.17（補題 C.2 も参照）により、実ガウス行列の符号補正された QR 因子 $O$（ライブラリは実環に対してこれを
   サンプリングする）は $O(k)$ 上でハール分布に従うから、$V\in O(k)$ に対し $VO$ は $O$ と同じ分布をもつ。
   $\mathrm{diag}(\pm1)$ を左から掛けることにより $i\ne l$ に対して $\mathbb E[O_{ij}O_{lm}]=0$ を得る。置換により
   $\mathbb E[O_{ij}O_{im}]$ は $i$ に依存しないことが分かり、$\sum_iO_{ij}O_{im}=\delta_{jm}$ である。したがって
   $\mathbb E[O_{ij}O_{lm}]=\delta_{il}\delta_{jm}/k$ であり、$\mathbb E\mathrm{Tr}(AO)^2=\sum_{ij}A_{ji}^2/k$ である。

3. ブロック対角作用素のトレースとトレースノルムはブロックにわたる和であり、$|X\oplus X'|=|X|\oplus|X'|$ である。
   $\lambda_{\min}(A\oplus A')=\min(\lambda_{\min}(A),\lambda_{\min}(A'))$ であり、$\Delta$ は $\lambda_{\min}$ について減少する。
   したがって $(X,Y)$ が $\delta$ で許容的、$(X',Y')$ が $\delta'\le\delta$ で許容的ならば、$(X\oplus X',Y\oplus Y')$ は
   $\delta$ で許容的であり、$z=z(X,Y)+z(X',Y')$ である。両方の組と $\delta'\le\delta$ にわたって上限をとればよい。

</details>

*計算：* 演習問題 1 と 2 の数値的検証 — [チュートリアル、演習問題](../../tutorial/ja/11_optimisation.md#x-11)。

**次へ：**[第12章](12_equivariant_learning.md)では同変機械学習を扱う：対称性と可換な層は
可換子環（第2章）の元であり、同じバッチ化された微分可能な作用素から構成される。
