[English](../11_optimisation.md) · [Español](../es/11_optimisation.md) · [Français](../fr/11_optimisation.md) · [Deutsch](../de/11_optimisation.md) · [中文](../zh/11_optimisation.md) · **日本語** · [Русский](../ru/11_optimisation.md)

# 11. 最適化とトレース不等式

トレース不等式とは、作用素の汎関数の極値についての主張である。3.7節では、もっともらしいが偽であるトレース不等式と、
その破れを見逃すサンプリングによる検証を示した。本章では、高次元でサンプリングが失敗する理由を説明し、それを厳密に
パラメータ化された制約集合上の最適化で置き換える。極値が閉じた形で分かっている場合にはそれを証明し、
最適化の結果を定理と照合できるようにする。

**本章で学ぶこと**：ライブラリがユニタリと、指定されたコントラストをもつスペクトルをどのようにパラメータ化するか、
そしてなぜすべての許容される作用素に到達できるのか。フォン・ノイマンのトレース不等式とそれが与える厳密な極値。
等号条件つきのトレース不等式。ハール分布のサンプルが極値から遠く離れたままである理由。そしてマルチスタート探索が何を証明するのか。

本章を通じて、$\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$（$D=\sum_ck_c\ge2$）とし、$\operatorname{Tr}$ は
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
| 実験 1 | $z=\lvert\operatorname{Tr}(XUY)\rvert-\operatorname{Tr}(XY)$ | $X,Y\ge0$、$U\in\mathcal M$ はユニタリ |
| 実験 2 | $z=\operatorname{Tr}\lvert XY\rvert-\operatorname{Tr}(\lvert X\rvert\,\lvert Y\rvert)$ | $X=UX_0$、$Y=VY_0$、$X_0,Y_0\ge0$、$U,V$ はユニタリ |
| 実験 3 | $z=\operatorname{Tr}(Y\lvert X\rvert Y)-\operatorname{Tr}\lvert YXY\rvert$ | $X=X^*$、$Y\ge0$ |

ここで $|X|$、$|Y|$ は $\lambda_{\max}=1$ に正規化し、それらのコントラストをパラメータとする。**包絡線** とは
$\overline z(\delta)=\sup\{z:\Delta(|X|)=\delta\}$ と $\underline z(\delta)=\inf\{z:\Delta(|X|)=\delta\}$ のことであり、その他の
変数はすべて自由に動かす。

## 11.1 制約集合とそのパラメータ化

**定義.** $\delta\in[0,1)$ と $\ell=\ell(\delta)$ に対し、**コントラスト殻** を次で定める。
$$
\mathcal P_\delta=\{X\in\mathcal M:X\ge0,\ \lambda_{\max}(X)=1,\ \Delta(X)=\delta\}
=\{X=X^*:\operatorname{spec}X\subseteq[\ell,1],\ 1,\ell\in\operatorname{spec}X\}.
$$
$\mathcal M$ の **スロット** とは $D$ 個の組 $(c,i)$（$i\le k_c$）のことである。$\lambda\in\mathbb R^D$ に対し
$\operatorname{diag}(\lambda)=\bigoplus_c\operatorname{diag}(\lambda_{c,1},\dots,\lambda_{c,k_c})$ とおく。

コントラスト $\delta$ を指定した `optimize.PositiveParam` は、各スタートごとに二つの異なるスロット
$s_\pm=(c_\pm,i_\pm)$ とハール分布に従うランダムなユニタリ $U_0\in\mathcal M$ を抽出し、$\theta\in\mathbb R^D$ と
$G=\bigoplus_cG_c$ を次のように写す。
$$
X(\theta,G)=U\operatorname{diag}(\lambda)U^*,\qquad U=U_0\,e^{G-G^*},\qquad
\lambda_{s_+}=1,\quad\lambda_{s_-}=\ell,\quad\lambda_s=\ell+(1-\ell)\tfrac{1+\sin\theta_s}2\ \text{otherwise}.
$$
`UnitaryParam` は写像 $G\mapsto U_0e^{G-G^*}$ 単独である。シグモイドとは異なり、$(1+\sin\theta)/2$ は有限の $\theta$ で $0$ と $1$
に到達する。極値をとるスペクトルはしばしば $\ell$ と $1$ にあるので、これは重要である。$[a,b]$ 内の自由なコントラスト
はパラメータが一つ増えるだけであり、$\delta=a+(b-a)(1+\sin\eta)/2$ とすればよい。

**補題 11.1.** 任意のユニタリ $U_0\in\mathcal M$ に対し、$G\mapsto U_0e^{G-G^*}$ は $\mathcal M$ のユニタリ群の上への
写像である。

*証明.* $U$ をユニタリとする。定理 1.3（ブロックごと）と系 1.4(2) により、
$W_c$ をユニタリ、$\varphi_{c,j}\in\mathbb R$ として $U_0^*U=\bigoplus_cW_c\operatorname{diag}(e^{i\varphi_{c,j}})_jW_c^*$ と書ける。このとき
$S=\bigoplus_cW_c\operatorname{diag}(i\varphi_{c,j})_jW_c^*$ は $S^*=-S$ と $e^S=U_0^*U$ を満たし（級数
$\sum_nS^n/n!$ は固有基底において項ごとに計算される）、$G=S/2$ とすれば $G-G^*=S$ である。$\square$

**命題 11.2（パラメータ化は全射である）.** $\delta$、$s_\pm$、$U_0$ を固定すると、
$(\theta,G)\mapsto X(\theta,G)$ の像は $\{X\in\mathcal P_\delta:1\in\operatorname{spec}X_{c_+},\ \ell\in\operatorname{spec}X_{c_-}\}$ であり、
$\mathcal P_\delta$ は相異なるスロットのすべての組にわたるこれらの像の和集合である。

*証明.* 像はこの集合に含まれる：$X$ の固有値は $\lambda_s\in[\ell,1]$ であり、$\lambda_{s_\pm}$ はブロック
$c_\pm$ にある。逆に、$X_c=W_c\operatorname{diag}(\mu_c)W_c^*$ と対角化し、各 $\mu_c$ を並べ替えて（$W_c$ を
$W_c\Pi_c$ で置き換える。$\Pi_c$ は置換行列）$\mu_{s_+}=1$ かつ $\mu_{s_-}=\ell$ となるようにする。これは可能である：$c_+=c_-$
かつ $\delta>0$ ならば、$1\ne\ell$ は同じブロックの二つの固有値であって異なる位置を占める。$\delta=0$ ならば、
すべての固有値は $1$ である。$\ell<1$ ならば、残りのスロットに対して $\theta_s=\arcsin\big(2(\mu_s-\ell)/(1-\ell)-1\big)$ とおけば
$\lambda=\mu$ となる（$\ell=1$ ならば $\theta$ は任意）。補題 11.1 により $U_0e^{G-G^*}=\bigoplus_cW_c$ を満たす $G$ が得られる。
*和集合:* $X\in\mathcal P_\delta$ に対し、固有値 $1$、$\ell$ を含むブロック $c_\pm$ と、その中の相異なるスロットを選ぶ。
$c_+=c_-$ かつ $\delta>0$ ならば、そのブロックは二つの異なる固有値をもつから $k_{c_+}\ge2$ である。$\delta=0$ ならば $X=1$ であり、
任意の二つのスロットでよい。$\square$

したがって、ランダムなスロットによってマルチスタート探索は殻全体を覆うことができ、制約はすべての反復で厳密に成り立つ。
`SelfAdjointParam` は $X=U\operatorname{diag}(\sigma\lambda)U^*$ を用いる。ここで符号 $\sigma\in\{\pm1\}^D$ はスタートごとに抽出され
（11.5節の符号反転によってのみ変更される）、したがって $|X|=U\operatorname{diag}(\lambda)U^*$ であり、（同じ証明により）
符号とスロットにわたる和集合は $\{X=X^*:|X|\in\mathcal P_\delta\}$ である。意味をもつのは $S=G-G^*$ だけである。歪自己共役な $S$ 上の写像
$S\mapsto U_0e^S$ の $S=0$ における微分は $S\mapsto U_0S$ であり、これはユニタリ群の接空間
$U_0\cdot\{S^*=-S\}$ への全単射であるから、この写像は $S=0$ の近くで局所微分同相である。$0$ から遠いところでは単射ではなく、
$S$ の二つの固有値が $2\pi i$ の零でない整数倍だけ異なるところで微分は退化する。そのため `extremize`
は各ラウンドの後に **再中心化** を行う（$U_0\leftarrow U_0e^{G-G^*}$、$G\leftarrow0$）。ケーリー・チャート
$S\mapsto U_0C(-S/2)$（演習 1）と極分解レトラクション（$S$ を $U_0(1+S)$ のユニタリ因子に写す。定理 1.8）は、
$0$ において同じ微分をもつ。

コードは制約を検証し、さらに命題 11.2 を構成的に検証する：証明のパラメータが目標の作用素を再現する。

```python
import math
import scipy.linalg
import torch
from torch_vn_algebra import TypeIAlgebra, optimize as opt

torch.manual_seed(0)
alg = TypeIAlgebra([3, 3], [3, 3], precision='double')          # M_3 (+) M_3, D = 6
ell = lambda delta: (1 - delta) / (1 + delta)                     # lambda_min for contrast delta
d = torch.tensor([0.0, 0.25, 0.5, 0.9], dtype=torch.float64)
x = opt.PositiveParam(alg, 4, delta=d).operator()                 # four starts, one contrast each
assert torch.allclose(x.michelson_contrast, d, atol=1e-12)        # the constraint holds exactly
assert torch.allclose(x.lambda_max, torch.ones(4, dtype=torch.float64), atol=1e-12)

P = opt.PositiveParam(alg, 1, delta=torch.tensor([0.5], dtype=torch.float64))
lo = ell(0.5)
mu = lo + (1 - lo) * torch.rand(1, 2, 3, dtype=torch.float64)    # target spectrum in [lo, 1] ...
(cp, ip), (cm, im) = P.spectrum.pin_max[0].tolist(), P.spectrum.pin_min[0].tolist()
mu[0, cp, ip], mu[0, cm, im] = 1.0, lo                            # ... with 1 and lo in the drawn slots
W = alg.random_unitary_operator(1).matrix                         # target eigenvectors
target = W @ torch.diag_embed(mu.to(W.dtype)) @ W.mH
with torch.no_grad():                                             # the parameters of the proof
    P.spectrum.theta.copy_(torch.asin((2 * (mu - lo) / (1 - lo) - 1).clamp(-1, 1)))
    V = P.unitary.base.mH @ W                                     # U0^* W = exp(S), G = S / 2
    S = torch.stack([torch.from_numpy(scipy.linalg.logm(V[0, c].numpy())) for c in range(2)])
    P.unitary.G.copy_(S[None] / 2)
err = (P.operator().matrix - target).abs().max().item()
print(f"slots of 1 and lo: {(cp, ip)}, {(cm, im)};  reconstruction error {err:.1e}")
assert err < 1e-10
```

> **物理的意味.** マイケルソン・コントラストは干渉縞の可視度 $(I_{\max}-I_{\min})/(I_{\max}+I_{\min})$ である。
> $\lambda_{\max}=1$ をもつ $0\le X\le1$（フィルター、減衰器、検出器の効率）に対しては、
> それはダイナミックレンジ $\lambda_{\min}/\lambda_{\max}$、すなわち条件数の逆数を定める。

## 11.2 フォン・ノイマンのトレース不等式

成分が $c_{ij}\ge0$ で、すべての行和と列和が $1$ である行列 $(c_{ij})$ を **二重確率行列** という。ユニタリ $G$ に対し、
$G$ の行と列は単位ベクトルであるから、$(|G_{ij}|^2)$ と $(|G_{ji}|^2)$ は二重確率行列である。

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
$\Sigma=\operatorname{diag}(s_i(A))$）をもつ：$A=W'|A|$（定理 1.8）、$|A|=V\Sigma V^*$（定理 1.3）と書き、$W=W'V$ とすればよい。

**定理 11.4（フォン・ノイマンのトレース不等式）.** $A,B\in M_k(\mathbb C)$ に対し
$|\operatorname{Tr}(AB)|\le\sum_is_i(A)\,s_i(B)$ である。

*証明.* $A=W_A\Sigma_AV_A^*$、$B=W_B\Sigma_BV_B^*$、$a_i=s_i(A)$、$b_j=s_j(B)$、$G=V_A^*W_B$、$H=V_B^*W_A$ とおく。
巡回性により $\operatorname{Tr}(AB)=\operatorname{Tr}(\Sigma_AG\Sigma_BH)=\sum_{i,j}a_iG_{ij}b_jH_{ji}$ であるから、
$|g||h|\le\frac12(|g|^2+|h|^2)$ と $a_ib_j\ge0$ により
$$
|\operatorname{Tr}(AB)|\le\sum_{i,j}a_ib_j\,d_{ij},\qquad d_{ij}=\tfrac12\big(|G_{ij}|^2+|H_{ji}|^2\big).
$$
$(d_{ij})$ は二重確率行列であり、補題 11.3 により $\sum_{ij}a_ib_jd_{ij}\le\sum_ia_ib_i$ である。$\square$

**系 11.5.** (a) 固有値 $p,q$ をもつ自己共役な $P,Q\in M_k(\mathbb C)$ に対し、
$\sum_ip^\downarrow_iq^\uparrow_i\le\operatorname{Tr}(PQ)\le\sum_ip^\downarrow_iq^\downarrow_i$。
(b) $P,Q\ge0$ に対し、$\operatorname{Tr}(PQ)\le\|PQ\|_1\le\sum_ip^\downarrow_iq^\downarrow_i$。

*証明.* (a) 定理 3.5 の証明と同様に、減少順に並べた固有基底に対する二重確率行列
$c_{ij}=|\langle u_i,v_j\rangle|^2$ を用いて $\operatorname{Tr}(PQ)=\sum_{ij}p_iq_jc_{ij}$ と書ける。補題 11.3 を適用する（その下界は
$\sum_ip^\downarrow_iq^\downarrow_{k+1-i}=\sum_ip^\downarrow_iq^\uparrow_i$ である）。(b) 第一の不等式は
定理 1.9 である。定理 1.9 と定理 11.4 により $\|PQ\|_1=\max_U|\operatorname{Tr}(P\,QU)|\le\sum_is_i(P)s_i(QU)$ であり、
$|QU|=(U^*Q^2U)^{1/2}=U^*QU$（定理 1.6）であるから $s_i(QU)=q^\downarrow_i$ である。$\square$

したがって命題 3.9 における破れは高々 $\sum_ix^\downarrow_i(y^\downarrow_i-y^\uparrow_i)$ である。実験 2 については、この
数がちょうど厳密な答えである。

**命題 11.6（スペクトルを固定した実験 2）.** 実験 2 では $|X|=X_0$、$|Y|=Y_0$ であり、
$z=\|X_0VY_0\|_1-\operatorname{Tr}(X_0Y_0)$ である。$(X_0)_c$ のスペクトル $x_c$ と $(Y_0)_c$ のスペクトル $y_c$ を固定し、
固有基底と $V$ を動かすとき、
$$
\sup z=-\inf z=\sum_c\sum_ix^\downarrow_{c,i}\big(y^\downarrow_{c,i}-y^\uparrow_{c,i}\big),
$$
であり、どちらも達成される。

*証明.* $X^*X=X_0U^*UX_0=X_0^2$ であるから $|X|=X_0$（定理 1.6）。同様に $|Y|=Y_0$ であり、$|UA|=|A|$ により
$\operatorname{Tr}|XY|=\|X_0VY_0\|_1$ である。すべての変数はブロックごとに作用し、$z$ はブロックにわたる和であるから、$C=1$ としてよい。
$Y'=VY_0V^*\ge0$（スペクトル $y$）とおく。$|AV|=V^*|A|V$ であるから $\|X_0VY_0\|_1=\|X_0Y'V\|_1=\|X_0Y'\|_1$ である。系 11.5
を組 $(X_0,Y')$ と $(X_0,Y_0)$ に適用し、$\|X_0Y'\|_1\ge\operatorname{Tr}(X_0Y')$ を用いると、
$$
\sum_ix^\downarrow_iy^\uparrow_i-\sum_ix^\downarrow_iy^\downarrow_i\le\|X_0Y'\|_1-\operatorname{Tr}(X_0Y_0)\le\sum_ix^\downarrow_iy^\downarrow_i-\sum_ix^\downarrow_iy^\uparrow_i .
$$
を得る。$J$ を逆順の置換 $Je_i=e_{k+1-i}$ とする。このとき $J\operatorname{diag}(w)=\operatorname{diag}(w_k,\dots,w_1)J$ であり、
$\operatorname{diag}(w)J$ のトレースノルムは $\sum_i|w_i|$ である。$X_0=\operatorname{diag}(x^\downarrow)$、$V=J$ とすると、
$Y_0=\operatorname{diag}(y^\uparrow)$ に対しては $X_0VY_0=\operatorname{diag}(x^\downarrow_iy^\downarrow_i)J$ かつ
$\operatorname{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\uparrow_i$（上界）であり、$Y_0=\operatorname{diag}(y^\downarrow)$ に対しては
$X_0VY_0=\operatorname{diag}(x^\downarrow_iy^\uparrow_i)J$ かつ $\operatorname{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\downarrow_i$
（下界）である。$\square$

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
である。$x_i\in[\ell_X,1]$、$y_i\in[\ell_Y,1]$ だからである。すべてのブロックにおいて、大きい方から $\lfloor k_c/2\rfloor$ 個の
固有値が $1$ で、小さい方から $\lfloor k_c/2\rfloor$ 個が $\ell_X$（それぞれ $\ell_Y$）であり、中央のものは
任意であるとき、等号が成り立つ。ある $k_c\ge2$ であるから、これらのスペクトルの全体での極値は $1$ と $\ell_X$（$\ell_Y$）である。命題 11.6
により結論を得る。$\square$

ここで $\Delta(|X|)=\delta$ に対し $1-\ell_X=2\delta/(1+\delta)$ である。$k_c$ が奇数のとき中央の固有値は関係しない。コードは
$k=3$、$C=2$ について定理 11.4、系 11.5、および命題 11.6 の極値配置を検証する。

```python
M5 = TypeIAlgebra([5], [5], precision='double')
A = M5.operator(torch.randn(2000, 1, 5, 5, dtype=torch.complex128))         # Ginibre, no structure
B = M5.operator(torch.randn(2000, 1, 5, 5, dtype=torch.complex128))
sA, sB = torch.linalg.svdvals(A.matrix[:, 0]), torch.linalg.svdvals(B.matrix[:, 0])   # decreasing
assert ((A @ B).trace.abs() <= (sA * sB).sum(-1) + 1e-12).all()             # Theorem 11.4
Ph, Qh = (M5.operator_from_eigenvalues(lambda k: torch.randn(2000, k), batch_size=2000, force_self_adjoint=True)
          for _ in range(2))
p, q = Ph.eigenvalues()[0], Qh.eigenvalues()[0]                             # increasing order
t = (Ph @ Qh).trace.real
assert (t <= (p * q).sum(-1) + 1e-12).all() and (t >= (p * q.flip(-1)).sum(-1) - 1e-12).all()   # Cor. 11.5(a)
Pp, Qp = Ph @ Ph, Qh @ Qh                                                   # positive, spectra p^2, q^2
nuc = (Pp @ Qp).trace_norm()
assert ((Pp @ Qp).trace.real <= nuc + 1e-12).all()                          # Cor. 11.5(b)
assert (nuc <= ((p ** 2).sort(-1)[0] * (q ** 2).sort(-1)[0]).sum(-1) + 1e-10).all()

alg3 = TypeIAlgebra([3, 3], [3, 3], precision='double')
x_dn = torch.tensor([1.0, 0.6, ell(0.5)], dtype=torch.float64)              # Delta = 0.5, middle value 0.6
y_dn = torch.tensor([1.0, 0.7, ell(0.8)], dtype=torch.float64)              # Delta = 0.8, middle value 0.7
J = torch.eye(3, dtype=torch.float64).flip(0)                               # reversal permutation
blocks = lambda m: alg3.from_blocks([m, m])                                 # the same block in both sectors

def z_exp2(X0, V, Y0):
    return ((X0 @ V @ Y0).trace_norm() - (X0 @ Y0).trace.real).item()

X0 = blocks(torch.diag(x_dn))
z_max = z_exp2(X0, blocks(J), blocks(torch.diag(y_dn.flip(0))))
z_min = z_exp2(X0, blocks(J), blocks(torch.diag(y_dn)))
exact = 2 * (x_dn * (y_dn - y_dn.flip(0))).sum().item()                     # Proposition 11.6
envelope = 2 * (1 - ell(0.5)) * (1 - ell(0.8))                              # Corollary 11.7, floor(3/2) = 1
print(f"Exp. 2, k = 3, C = 2: z_max = {z_max:.6f}, z_min = {z_min:.6f}, formula {envelope:.6f}")
assert abs(z_max - exact) < 1e-12 and abs(z_min + exact) < 1e-12 and abs(exact - envelope) < 1e-12
```

## 11.3 等号条件つきのトレース不等式

$X=X^*$ に対し、$f_+(t)=\max(t,0)$、$f_-(t)=\max(-t,0)$ として $X_\pm=f_\pm(X)$ とおく。このとき $X_\pm\ge0$、$X=X_+-X_-$、
$|X|=X_++X_-$、$X_+X_-=0$ である（定理 1.5(c)）。

**補題 11.8.** $M_n(\mathbb C)$ の $A,B\ge0$ に対し $\|A-B\|_1\le\operatorname{Tr}A+\operatorname{Tr}B$ であり、等号は
$AB=0$ のときかつそのときに限り成り立つ。

*証明.* $H=A-B$ とし、$(0,\infty)$ 上で $g=1$、それ以外で $g=0$ として $Q=g(H)$ とおく。$|t|=t(2g(t)-1)$ であるから $|H|=H(2Q-1)$ であり、
$$
\|H\|_1=\operatorname{Tr}\big(A(2Q-1)\big)-\operatorname{Tr}\big(B(2Q-1)\big)
=\operatorname{Tr}A+\operatorname{Tr}B-2\operatorname{Tr}\big((1-Q)A(1-Q)\big)-2\operatorname{Tr}(QBQ).
$$
である。引かれているトレースは $\ge0$ であり、$\operatorname{Tr}(RAR)=\|A^{1/2}R\|_2^2$ が零となるのは $AR=0$ のときかつそのときに限る。したがって等号が成り立つのは
$A(1-Q)=0$ かつ $BQ=0$ のときかつそのときに限り、そのとき $AB=AQB=A(BQ)^*=0$ である。逆に $AB=0$ ならば $BA=0$ であり、
$(A-B)^2=A^2+B^2=(A+B)^2$ であるから、定理 1.6 により $|A-B|=A+B$ である。$\square$

**命題 11.9（実験 3）.** $\mathcal M$ の $X=X^*$ と $Y\ge0$ に対し
$\operatorname{Tr}|YXY|\le\operatorname{Tr}(Y|X|Y)$ であり、等号は $YX_+Y^2X_-Y=0$ のときかつそのときに限り成り立つ。特に $XY=YX$ ならば等号が成り立つ。

*証明.* $M_D(\mathbb C)$ において $A=YX_+Y$ と $B=YX_-Y$ は正であり（系 1.7 を $T=X_\pm^{1/2}Y$ として適用）、$YXY=A-B$
かつ $\operatorname{Tr}(Y|X|Y)=\operatorname{Tr}A+\operatorname{Tr}B$ である。補題 11.8 を適用する。$XY=YX$ ならば $Y$ は
$X_\pm$ と可換であり（定理 1.5(d)）、$AB=Y^4X_+X_-=0$ である。$\square$

したがって実験 3 では $z\ge0$ であり、$\underline z=0$ である（$Y=1$ とする）。等号条件は可換性より弱い：
$X=\operatorname{diag}(1,-1,0)$ と、行が $(1,0,a)$、$(0,1,0)$、$(a,0,1)$（$|a|<1$）である $Y$ に対し、
$AB=(Ye_1)\langle Ye_1,Ye_2\rangle(Ye_2)^*=0$ であるが、$a\ne0$ ならば $XY\ne YX$ である。

```python
def z_exp3(X, Y):
    return (Y @ X.abs() @ Y).trace.real - (Y @ X @ Y).trace_norm()

Xs = alg3.operator_from_eigenvalues(lambda k: torch.randn(500, k), batch_size=500, force_self_adjoint=True)
Yp = alg3.operator_from_eigenvalues(lambda k: torch.rand(500, k), batch_size=500, force_positive=True)
assert (z_exp3(Xs, Yp) >= -1e-12).all()                                     # Proposition 11.9

M3 = TypeIAlgebra([3], [3], precision='double')
Xe = M3.operator(torch.diag(torch.tensor([1.0, -1.0, 0.0]))[None, None])
Ye = M3.operator(torch.tensor([[1.0, 0, 0.5], [0, 1, 0], [0.5, 0, 1]])[None, None])
print(f"equality without commutation: z = {z_exp3(Xe, Ye).item():.1e}, "
      f"||XY - YX|| = {(Xe @ Ye - Ye @ Xe).operator_norm().item():.2f}")
assert abs(z_exp3(Xe, Ye).item()) < 1e-12 and (Xe @ Ye - Ye @ Xe).operator_norm().item() > 0.1
```

## 11.4 サンプリング対最適化

$|\operatorname{Tr}(AU)|$ は、その最大値 $\|A\|_1$ と比べて典型的にどれほどの大きさか。ライブラリは $U$ を、
**ジニブル行列** $Z$（成分は独立で密度
$\pi^{-1}e^{-|z|^2}$ をもつ。Mezzadri の処方）の QR 分解 $Z=QR$（$R_{ii}>0$）の因子 $Q$ として、各ブロックで独立にサンプリングする。

**補題 11.10.** $Q$ は $U(k)$ 上でハール分布に従う。特に、任意の固定されたユニタリ $V\in M_k(\mathbb C)$ に対し、$VQ$ は
$Q$ と同じ分布をもつ。独立なブロックをもつユニタリ $\bigoplus_cQ_c$ は、$\mathcal M$ のユニタリ群上でハール分布に従う。

*証明.* これは定理 5.3（左不変性はその証明の同変性の段階である）と系 5.4 である。$\square$

（測度の集中とこれらの群についてさらに詳しくは、E. Meckes, *The Random Matrix Theory of the Classical
Compact Groups*, Cambridge 2019 を参照。）

**命題 11.11.** ハール分布に従う $U\in U(k)$ と $A\in M_k(\mathbb C)$ に対し、$\mathbb E\operatorname{Tr}(AU)=0$ かつ
$\mathbb E|\operatorname{Tr}(AU)|^2=\|A\|_2^2/k$ である。

*証明.* 系 5.5 により $\mathbb EU_{ij}=0$ である。命題 5.7（$C=1$）を $X=|e_j\rangle\langle e_m|$ に適用すると、
1次のワインガルテン公式
$\mathbb E[U_{ij}\bar U_{lm}]=\mathbb E(UXU^*)_{il}=\delta_{il}\delta_{jm}/k$ を得る。$k\ge2$ に対しては、これは
2次の公式（定理 5.10）において $i_2=i'_2$、$j_2=j'_2$ とおき両者について和をとることからも従う。
$\sum_{i_2,j_2}|U_{i_2j_2}|^2=k$ だからである。$\operatorname{Tr}(AU)=\sum_{ij}A_{ji}U_{ij}$ を用いると、
$\mathbb E|\operatorname{Tr}(AU)|^2=\sum_{ijlm}A_{ji}\bar A_{ml}\mathbb E[U_{ij}\bar U_{lm}]=\frac1k\sum_{ij}|A_{ji}|^2$ である。$\square$

**系 11.12.** $M_k(\mathbb C)$ の $X,Y\ge0$ のスペクトルが $[a,b]$（$a>0$）に含まれるとし、$U$ をハール分布に従うとする。このとき
$$
\frac{\big(\mathbb E|\operatorname{Tr}(XUY)|^2\big)^{1/2}}{\max_V|\operatorname{Tr}(XVY)|}=\frac{\|YX\|_2}{\sqrt k\,\|YX\|_1}\le\min\Big(\frac1{\sqrt k},\frac{b^2}{a^2k}\Big),
\qquad
\mathbb P\big(|\operatorname{Tr}(XUY)|\ge\operatorname{Tr}(XY)\big)\le\frac{b^4}{a^4k^2}.
$$

*証明.* 等式は $A=YX$ に対する命題 11.11 と命題 3.9 である。$\sum_is_i^2\le(\sum_is_i)^2$ から
$1/\sqrt k$ を得る。$s_1(A)\le\|Y\|\|X\|\le b^2$ かつ $s_k(A)=\min_{\|v\|=1}\|YXv\|\ge a^2$ であるから、$\|A\|_2^2\le kb^4$ かつ
$\|A\|_1\ge ka^2$ である。最後に $\operatorname{Tr}(XY)=\operatorname{Tr}(X^{1/2}YX^{1/2})\ge a\operatorname{Tr}X\ge ka^2$ であり、
$|\operatorname{Tr}(XUY)|^2$ に対するマルコフの不等式により、確率は $(kb^4/k)/(ka^2)^2$ で抑えられる。$\square$

したがって、$N$ 個のハール・サンプルが $|\operatorname{Tr}(XUY)|\le\operatorname{Tr}(XY)$ の破れを見つける確率は高々
$Nb^4/(a^4k^2)$ である。破れは可換でないすべての組に対して存在するにもかかわらず、である。実際はもっと極端である：$f(U)=|\operatorname{Tr}(AU)|$
は $\|U-V\|_2$ に関してリプシッツ定数 $L=\|A\|_2$ のリプシッツ関数であるから、絶対定数 $c>0$ に対し
$\mathbb P(f\ge\mathbb Ef+t)\le e^{-ckt^2/L^2}$ が成り立つ（測度の集中。証明なしで述べる。Meckes, 前掲書, 第5章）。
ここでは指数は $k^2a^4/b^4$ のオーダーである。コードは命題 11.11 を検証し、$[0.2,1.2]$ 上で一様なスペクトル（$b/a=6$）に対して、
系 11.12 の比、$\sup z$、および $2000$ 個のサンプルのうち最良の $z$ を出力する。後二者はいずれも $\operatorname{Tr}(XY)$ に対する相対値である。

```python
def random_pair(alg_k):
    pos = lambda k: 0.2 + torch.rand(1, k)                                  # spectrum in [a, b] = [0.2, 1.2]
    return tuple(alg_k.operator_from_eigenvalues(pos, force_positive=True) for _ in range(2))

M6 = TypeIAlgebra([6], [6], precision='double')
Xf, Yf = random_pair(M6)
N = 20_000
t2 = (Yf @ Xf @ M6.random_unitary_operator(N)).trace.abs() ** 2
m2, err, exact2 = t2.mean().item(), t2.std().item() / math.sqrt(N), (Yf @ Xf).frobenius_norm().item() ** 2 / 6
print(f"E|Tr(AU)|^2 = {m2:.4f} +- {err:.4f}   (exact ||A||_2^2 / k = {exact2:.4f})")
assert abs(m2 - exact2) < 5 * err                                           # Proposition 11.11

print(" k   rms/max   sup z/TrXY   best of 2000 samples")
for k in (2, 4, 8, 16, 32):
    Mk = TypeIAlgebra([k], [k], precision='double')
    Xk, Yk = random_pair(Mk)
    trXY, nuc = (Xk @ Yk).trace.real.item(), (Yk @ Xk).trace_norm().item()
    ratio = (Yk @ Xk).frobenius_norm().item() / (math.sqrt(k) * nuc)
    best = (Xk @ Mk.random_unitary_operator(2000) @ Yk).trace.abs().max().item()
    print(f"{k:2d}   {ratio:.4f}    {(nuc - trXY) / trXY:.2e}     {(best - trXY) / trXY:+.3f}")
    assert ratio <= min(1 / math.sqrt(k), 36 / k) + 1e-12                   # Corollary 11.12
    assert nuc > trXY and best < nuc and (k < 4 or best < trXY)            # from k = 4 on: no violation found
```

比は $1/k$ のように減少し、$k=4$ 以降では最良のサンプルは $z=0$ にさえ届かない。極値はハール測度が
きわめて小さい集合であり、勾配法はそこに直接到達する。組 `Xf, Yf` については、$U$ に関する両方の極値が分かっている。

**補題 11.13.** $k\ge2$ ならば、任意の $A\in M_k(\mathbb C)$ に対し $\min_U|\operatorname{Tr}(AU)|=0$ である。

*証明.* $A=W\Sigma V^*$ とし、$P$ を巡回シフト $Pe_i=e_{(i\bmod k)+1}$ とする。$k\ge2$ ならばその対角成分は零である。
$U=VPW^*$ に対し $\operatorname{Tr}(AU)=\operatorname{Tr}(\Sigma P)=0$ である。$\square$

したがって $\inf_Uz=-\operatorname{Tr}(XY)$、$\sup_Uz=\|YX\|_1-\operatorname{Tr}(XY)$ である。`extremize` を一回呼び出せば両方が求まる。
$12$ 個のスタートが $z$ を最大化し、$12$ 個が $-z$ を最大化する。サンプリングは下限には近づく。下限は
典型的な値のすぐそばにあるからである（$|\operatorname{Tr}(XUY)|$ は典型的に $\|YX\|_2/\sqrt k$ のオーダー）。しかし上限には近づかない。

```python
trXY = (Xf @ Yf).trace.real
sup_exact, inf_exact = ((Yf @ Xf).trace_norm() - trXY).item(), -trXY.item()
U6 = opt.UnitaryParam(M6, 24)                                               # 24 starts
pm = torch.tensor([1.0, -1.0], dtype=torch.float64).repeat_interleave(12)   # 12 maximise z, 12 maximise -z
best = opt.extremize(lambda: pm * ((Xf @ U6.operator() @ Yf).trace.abs() - trXY), [U6],
                     maximize=True, rounds=6, steps=30).reshape(2, 12).max(1)[0]
best_up, best_down = best[0].item(), -best[1].item()
z_mc = (Xf @ M6.random_unitary_operator(N) @ Yf).trace.abs() - trXY
print(f"sup z: exact {sup_exact:.6e}, ascent {best_up:.6e}, best sample {z_mc.max().item():+.4f}")
print(f"inf z: exact {inf_exact:.6f}, descent {best_down:.6f}, worst sample {z_mc.min().item():+.6f}")
assert best_up <= sup_exact + 1e-12 and best_down >= inf_exact - 1e-12    # attained values: one-sided
assert sup_exact - best_up < 1e-3 * sup_exact and best_down - inf_exact < 1e-3 * abs(inf_exact)
assert z_mc.min().item() - inf_exact < 1e-2 * abs(inf_exact)               # sampling finds the typical side ...
assert z_mc.max().item() < 0                                               # ... but misses the violation
```

> **物理的意味.** 測度の集中は、統計力学における **典型性** の数学的な形である：高次元では、ほとんどすべての状態や
> 時間発展が巨視的な量に同じ値を与える。鋭い不等式は、微調整された配置——揃った固有基底、2準位のスペクトル——によって決まり、
> ランダムなユニタリ（ランダム回路、ランダムな散乱行列）はそうした配置を事実上決して生み出さない。

## 11.5 包絡線のためのマルチスタート探索

`optimize.extremize(objective, params, maximize)` は Adam を `steps` ステップずつ `rounds` ラウンド実行し（学習率は
`lr * decay**round`）、各ラウンドの後にユニタリ・チャートを再中心化し、**スタートごとに** 達成された最良の値を返す。
バッチの各要素は、それぞれ固有のスロット、符号、$U_0$、コントラストをもつ独立なスタートである。
`SelfAdjointParam` の符号は離散的である：ラウンドの間に、`discrete_step` はスタートごとに単一の符号反転を八回、
一つずつ順に試し、目的関数を改善する反転をそのつど保持する。

すべての反復は許容的である（命題 11.2）から、返される値はすべて達成された値である：最大化は上限の
**下界** を、最小化は下限の上界を与え、それ以上のものではない——局所探索は局所極値で止まりうる。
逆向きの評価には、系 11.7 や命題 11.9 のような証明か、凸緩和（11.6節）が必要である。計算された包絡線は
(i) 閉じた形と照合し、(ii) 小さな解の直和のような手で構成した配置と照合し（包絡線はそれを下回ってはならない）、
(iii) 構造化されたサンプル（2準位のスペクトル、置換基底）と照合する。後者はハール・サンプルとは異なり、
どの次元でも極値に近づく（`scripts/bounds_checks.py`）。

**命題 11.14（実験 1 の下側包絡線）.** すべての $k_c\ge2$ とする。$X\in\mathcal P_\delta$、
$\lambda_{\max}(Y)=1$ を満たす $Y\ge0$、およびユニタリ $U\in\mathcal M$ にわたって $\ \min z=-\big(D-2\delta/(1+\delta)\big)$ である。

*証明.* $Y\le1$ から $X^{1/2}YX^{1/2}\le X$ が従うので、$\operatorname{Tr}(XY)\le\operatorname{Tr}X\le(D-1)+\ell(\delta)$ である。
$X$ の一つの固有値は $\ell(\delta)$ であり、他は $\le1$ だからである。よって $z\ge-\operatorname{Tr}(XY)\ge-(D-2\delta/(1+\delta))$ である。
等号は、$X$ の固有値が一度だけ $\ell(\delta)$ でそれ以外は $1$、$Y=1$、かつ $U=\bigoplus_cU_c$ が
$\operatorname{Tr}(X_cU_c)=0$ を満たす（補題 11.13）とき成り立つ。$\square$

コードは $M_2\oplus M_2$ 上で命題 11.14 を、$M_3\oplus M_3$ 上で系 11.7 を（両方のコントラストを指定、
$k=3$ は奇数）、そして命題 11.9 を検証する。上限と下限は、上と同様に一回の呼び出しで求める。キーワード `is_self_adjoint` は
トレースノルムのアルゴリズム（特異値、あるいは $YXY$ に対しては固有値）を、各ステップで判定することなく選択する。最後の
検証は型 (ii) である：実験 3 の $z$ は直和について加法的であり、直和のコントラストは各成分のコントラストの
最大値である（演習 3）から、
$$
\overline z_{M_2\oplus M_2}(\delta)\ \ge\ \overline z_{M_2}(\delta)+\max_{\delta'\le\delta}\overline z_{M_2}(\delta').
$$

```python
torch.manual_seed(1)
S = 6                                                                       # starts per grid point
grid = torch.tensor([0.0, 0.5, 0.9], dtype=torch.float64)
alg22 = TypeIAlgebra([2, 2], [2, 2], precision='double')                    # D = 4

# Exp. 1, lower envelope along Delta(X); Delta(Y) free in [0, 0.999]
X1, Y1, U1 = (opt.PositiveParam(alg22, 3 * S, delta=grid.repeat_interleave(S)),
              opt.PositiveParam(alg22, 3 * S), opt.UnitaryParam(alg22, 3 * S))
def z1():
    x, y = X1.operator(), Y1.operator()
    return (x @ U1.operator() @ y).trace.abs() - (x @ y).trace.real
inf1 = opt.extremize(z1, [X1, Y1, U1], maximize=False, rounds=4, steps=30).reshape(3, S).min(1)[0]
exact1 = -(4 - 2 * grid / (1 + grid))                                       # Proposition 11.14
print("Exp. 1  inf:", [f"{v:.4f}" for v in inf1.tolist()], " exact:", [f"{v:.4f}" for v in exact1.tolist()])
assert (inf1 >= exact1 - 1e-12).all() and (inf1 - exact1).max() < 4e-3

# Exp. 2 on M_3 (+) M_3, Delta(X0) = 0.5, Delta(Y0) = 0.8: sup and inf in one call
sign = torch.tensor([1.0, -1.0], dtype=torch.float64).repeat_interleave(S)
X2, Y2, V2 = (opt.PositiveParam(alg3, 2 * S, delta=torch.full((2 * S,), 0.5, dtype=torch.float64)),
              opt.PositiveParam(alg3, 2 * S, delta=torch.full((2 * S,), 0.8, dtype=torch.float64)),
              opt.UnitaryParam(alg3, 2 * S))
def z2():
    x0, y0 = X2.operator(), Y2.operator()
    M = alg3.operator((x0 @ V2.operator() @ y0).matrix, is_self_adjoint=False)
    return sign * (M.trace_norm() - (x0 @ y0).trace.real)
best2 = opt.extremize(z2, [X2, Y2, V2], maximize=True, rounds=4, steps=30).reshape(2, S).max(1)[0]
sup2, inf2 = best2[0].item(), -best2[1].item()
print(f"Exp. 2  sup {sup2:.6f}, inf {inf2:.6f}, exact +-{envelope:.6f}")
assert sup2 <= envelope + 1e-12 and inf2 >= -envelope - 1e-12               # Corollary 11.7: one-sided ...
assert envelope - sup2 < 1e-3 * envelope and inf2 + envelope < 1e-3 * envelope   # ... and reached

# Exp. 3: maximise sign * z at the contrasts `deltas` of |X|; Delta(Y) free in [0, 0.999]
def exp3(alg_, deltas, sign):
    X3, Y3 = opt.SelfAdjointParam(alg_, len(deltas), delta=deltas), opt.PositiveParam(alg_, len(deltas))
    def z3():
        y = Y3.operator()
        M = alg_.operator((y @ X3.operator() @ y).matrix, is_self_adjoint=True)
        return sign * ((y @ X3.abs_operator() @ y).trace.real - M.trace_norm())
    return opt.extremize(z3, [X3, Y3], maximize=True, rounds=4, steps=30).reshape(2, S).max(1)[0]

best3 = exp3(alg22, torch.full((2 * S,), 0.5, dtype=torch.float64), sign)
sup_M22, inf3 = best3[0].item(), -best3[1].item()
s_M2 = exp3(TypeIAlgebra([2], [2], precision='double'),
            torch.tensor([0.0, 0.5], dtype=torch.float64).repeat_interleave(S), 1.0)
direct_sum = (s_M2[1] + s_M2.max()).item()                                 # delta = 0.5, best delta' <= 0.5
conj = [(1 - dl) * (2 * 0.999 / 1.999) ** 2 for dl in (0.0, 0.5)]          # unproven closed form, d = 0.999
print(f"Exp. 3  inf {inf3:.1e};  sup for M_2 at delta = 0, 0.5: {s_M2[0]:.5f}, {s_M2[1]:.5f}"
      f"  (closed form {conj[0]:.5f}, {conj[1]:.5f})")
print(f"        sup for M_2 (+) M_2 at 0.5: search {sup_M22:.5f}, direct sum {direct_sum:.5f}")
assert inf3 >= -1e-12 and inf3 < 1e-6                                       # Proposition 11.9
assert sup_M22 >= direct_sum - 5e-3 * direct_sum                           # (ii), up to 0.5 %
```

$M_2\oplus M_2$ 上の探索は直和を約 $0.2\,\%$ 下回る：探索は完全には収束しておらず、
$\overline z_{M_2\oplus M_2}(0.5)$ の最良の下界は大きい方の数である——これを検出することが検証 (ii)
の目的である。$M_2$ に対する値は $(1-\delta)\big(2d/(1+d)\big)^2=(1-\delta)(1-\ell(d))^2$ と $10^{-4}$ まで一致する。ここで
$d=0.999$ は $Y$ に許される最大のコントラストである。これは解説論文が（実環において）数値的に観察したが
証明できなかった閉じた形であり、コードはこれを assert しない。より大きなブロックについては、実験 1 と
実験 3 の上側包絡線は、我々の知る限り閉じた形では知られていない。その場合、計算された曲線は下界と予想であって、
それ以上のものではない。

## 11.6 比較の対象

**Manopt**、**pymanopt**、**geoopt** は行列多様体（シュティーフェル多様体、直交群、ユニタリ群）を、
QR、極分解、ケーリーなどのレトラクションや、リーマン的ソルバー（共役勾配法、信頼領域法、リーマン
Adam）とともに提供しており、これらは再中心化したチャートでの Adam より効率的であることが多い。SDP ソルバーと組み合わせた **cvxpy** は、
作用素変数の凸問題を大域的最適性まで解き、ソルバーの許容誤差の範囲で双対証明書を与える。
半正定値緩和は、局所探索では得られない方向の厳密な限界を与える。**NCTSSOS**、**ncpol2sdpa**、**NCSOStools** は、非可換
多項式の固有値とトレースをすべての次元で一度に評価する（0.3節）。ライブラリが加えるのは、$\bigoplus_cM_{k_c}$ のスペクトル制約と
コントラスト制約の厳密なパラメータ化（極値はどのセクターにもありうる）、ブロック構造を
保つユニタリ、そして多数のスタートと制約値に対する、CPU または GPU 上での一回のバッチ計算である。ライブラリは極値の
候補と上限の下界を生み出す。残りは証明か緩和が補わなければならない。

## 演習問題

1. （**ケーリー変換.**）$S^*=-S$ に対し、$C(S)=(1-S)(1+S)^{-1}$ がユニタリであり、$S\mapsto C(S)$ が
   固有値 $-1$ をもたないユニタリ全体の上への全単射であることを示せ。したがって指数写像（補題 11.1）とは異なり、ケーリー・チャート
   $U_0C(S)$ は $-U_0$ を取りこぼす。
2. （**実ハール測度.**）$O(k)$ 上でハール分布に従う $O$（実ジニブル行列の QR）と実行列 $A$ に対し、
   $\mathbb E\operatorname{Tr}(AO)^2=\|A\|_2^2/k$ を示し、実環で確かめよ。
3. （**直和.**）実験 3 の $z$ が直和について加法的であること、および $\lambda_{\max}=1$ をもつ $A,A'\ge0$ に対し
   $\Delta(A\oplus A')=\max(\Delta(A),\Delta(A'))$ であることを示せ。11.5節の直和に関する
   不等式を導け（$\Delta(Y)$ は固定された区間 $[0,b]$ 内で自由とする）。

<details>
<summary>解答</summary>

1. $S$ の固有値は純虚数であるから $1\pm S$ は可逆である。これらは可換であり、$(1+S)^*=1-S$ である。したがって
   $C^*C=(1-S)^{-1}(1+S)(1-S)(1+S)^{-1}=1$ であり、$1+C=2(1+S)^{-1}$ は可逆であるから $-1\notin\operatorname{spec}C$ である。
   逆に、$1+U$ が可逆であるユニタリ $U$ に対し、$S=(1-U)(1+U)^{-1}$ は歪自己共役である：
   $S^*=(1+U^*)^{-1}(1-U^*)=(U+1)^{-1}U\,U^*(U-1)=-S$。$C$ は $S$ と可換であるから、$C(1+S)=1-S$ は
   $S(1+C)=1-C$、すなわち $S=(1-C)(1+C)^{-1}$ と同値である：二つの写像は互いに逆である。

```python
Sk = torch.randn(4, 4, dtype=torch.complex128)
Sk = Sk - Sk.mH                                                   # skew-adjoint
I4 = torch.eye(4, dtype=torch.complex128)
Cay = (I4 - Sk) @ torch.linalg.inv(I4 + Sk)
assert torch.allclose(Cay.mH @ Cay, I4, atol=1e-12)
assert torch.allclose((I4 - Cay) @ torch.linalg.inv(I4 + Cay), Sk, atol=1e-10)
assert torch.linalg.eigvals(Cay).add(1).abs().min() > 1e-6        # -1 is never an eigenvalue
```

2. 定理 5.3 とその証明は、実ガウス行列（密度 $\propto e^{-\|Z\|_2^2/2}$）、直交
   群、直交行列 $V$ に対しても成り立つ（定理 5.3 の後の注意）から、$VO$ は $O$ と同じ分布をもつ。
   $\operatorname{diag}(\pm1)$ を左から掛けることにより $i\ne l$ に対して $\mathbb E[O_{ij}O_{lm}]=0$ を得る。置換により
   $\mathbb E[O_{ij}O_{im}]$ は $i$ に依存しないことが分かり、$\sum_iO_{ij}O_{im}=\delta_{jm}$ である。したがって
   $\mathbb E[O_{ij}O_{lm}]=\delta_{il}\delta_{jm}/k$ であり、$\mathbb E\operatorname{Tr}(AO)^2=\sum_{ij}A_{ji}^2/k$ である。

```python
Mr = TypeIAlgebra([5], [5], complex_valued=False, precision='double')
Ar = Mr.operator(torch.randn(1, 1, 5, 5, dtype=torch.float64))
tr2 = (Ar @ Mr.random_unitary_operator(20_000)).trace ** 2
assert abs(tr2.mean().item() - Ar.frobenius_norm().item() ** 2 / 5) < 5 * tr2.std().item() / math.sqrt(20_000)
```

3. ブロック対角作用素のトレースとトレースノルムはブロックにわたる和であり、$|X\oplus X'|=|X|\oplus|X'|$ である。
   $\lambda_{\min}(A\oplus A')=\min(\lambda_{\min}(A),\lambda_{\min}(A'))$ であり、$\Delta$ は $\lambda_{\min}$ について減少する。
   したがって $(X,Y)$ が $\delta$ で許容的、$(X',Y')$ が $\delta'\le\delta$ で許容的ならば、$(X\oplus X',Y\oplus Y')$ は
   $\delta$ で許容的であり、$z=z(X,Y)+z(X',Y')$ である。両方の組と $\delta'\le\delta$ にわたって上限をとればよい。

</details>

**次へ：**第12章（準備中）では同変機械学習を扱う：対称性と可換な層は
可換子環（第2章）の元であり、同じバッチ化された微分可能な作用素から構成される。
