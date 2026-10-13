[English](../01_matrix_algebras.md) · [Español](../es/01_matrix_algebras.md) · [Français](../fr/01_matrix_algebras.md) · [Deutsch](../de/01_matrix_algebras.md) · [中文](../zh/01_matrix_algebras.md) · **日本語** · [Русский](../ru/01_matrix_algebras.md)

# 1. 行列環と関数計算

本章では記号を定め、以降のすべての章で用いる $M_n(\mathbb C)$ に関する事実、すなわちスペクトル定理、関数計算、極分解、および特異値によって定義されるノルムを証明する。各結果の後には、対応するライブラリのコードと数値的な検証を置く。

**本章で学ぶこと**：指定したスペクトルをもつ作用素の構成法、作用素への関数の適用法、作用素のバッチの格納法、そして計算を厳密な結果と照合して検証する方法。

**規約.** 内積 $\langle x,y\rangle=\sum_i\bar x_iy_i$ は第二引数について線形とする。$A^*$ は共役転置であり、したがって $\langle x,Ay\rangle=\langle A^*x,y\rangle$ である。$e_1,\dots,e_n$ は標準基底、$1$ は単位行列、$\operatorname{Tr}$ は通常のトレースを表す。

## 1.1 環 $M_n(\mathbb C)$

$M_n(\mathbb C)$ は複素 $n\times n$ 行列のなす環であり、対合 $A\mapsto A^*$ と作用素ノルム $\|A\|=\sup_{\|x\|=1}\|Ax\|$ を備える。元 $A$ は次のように呼ばれる。

| 名称 | 定義 |
|---|---|
| 正規 | $A^*A=AA^*$ |
| 自己共役 | $A=A^*$ |
| 正、$A\ge0$ | すべての $x$ に対し $\langle x,Ax\rangle\ge0$ |
| ユニタリ | $A^*A=AA^*=1$ |
| 射影 | $A=A^*=A^2$ |

$B-A\ge0$ のとき $A\le B$ と書き、$A\ge0$ かつ $A$ が可逆のとき $A>0$ と書く。

**補題 1.1（C\*-恒等式）.** $\|A^*\|=\|A\|$ かつ $\|A^*A\|=\|A\|^2$ である。

*証明.* $\|Ax\|=\sup_{\|y\|=1}|\langle y,Ax\rangle|=\sup_{\|y\|=1}|\langle A^*y,x\rangle|$
であるから、$\|A\|=\sup_{\|x\|=\|y\|=1}|\langle A^*y,x\rangle|=\|A^*\|$ である。次に、$\|x\|=1$ に対し
$\|Ax\|^2=\langle x,A^*Ax\rangle\le\|A^*A\|$ であるから、$\|A\|^2\le\|A^*A\|\le\|A^*\|\,\|A\|=\|A\|^2$ を得る。$\square$

> **物理的意味.** 量子力学では、$n$ 個の区別可能な配置をもつ系は $\mathbb C^n$ によって記述される。自己共役行列は
> *オブザーバブル*（エネルギー、スピン、格子上の位置）であり、その固有値は測定の取りうる結果である。射影は yes/no 型の問いである。ユニタリ
> 行列は対称性および時間発展であり、ハミルトニアン $H$ に対しては $e^{-itH}$ である。
> トレース 1 の正行列は *状態*（密度行列、第3章）である。

## 1.2 スペクトル定理

**補題 1.2（シューア）.** 任意の $A\in M_n(\mathbb C)$ に対し、$U^*AU$ が上三角となるユニタリ $U$ が存在する。

*証明.* $n$ に関する帰納法による。$n=1$ のときは自明である。$A$ の特性多項式は根 $\lambda$ をもつ（代数学の基本定理）ので、ある単位ベクトル $v$ に対し $Av=\lambda v$ である。$v$ を正規直交基底 $v,w_2,\dots,w_n$ に延長し、これらを列とするユニタリを
$V$ とする。すると $A'\in M_{n-1}(\mathbb C)$ により
$V^*AV=\begin{pmatrix}\lambda&*\\0&A'\end{pmatrix}$ と書ける。
帰納法の仮定より、あるユニタリ $U'$ に対して $U'^*A'U'$ は上三角であり、
$U=V\,(1\oplus U')$ が求めるものである。$\square$

**定理 1.3（スペクトル定理）.** $A$ が正規であるための必要十分条件は、あるユニタリ $U$ により
$A=U\,\mathrm{diag}(\lambda_1,\dots,\lambda_n)\,U^*$ と書けることである。
$\lambda_i$ は重複度込みの $A$ の固有値である。

*証明.* $\Lambda$ を対角行列として $A=U\Lambda U^*$ ならば、
$A^*A=U\bar\Lambda\Lambda U^*=U\Lambda\bar\Lambda U^*=AA^*$ である。逆に $A$
を正規とし、$T=U^*AU$ を上三角とする（補題 1.2）。$T$ は正規である。
$(1,1)$ 成分を比較すると $(T^*T)_{11}=|t_{11}|^2$ かつ $(TT^*)_{11}=\sum_j|t_{1j}|^2$ であるから、
$j>1$ に対し $t_{1j}=0$ である。よって $T=t_{11}\oplus T'$（$T'$ は正規かつ上三角）となり、帰納法により $T$ は対角である。$U^*AU$ の対角成分は
$\det(A-\lambda)=\det(T-\lambda)$ の根からなる。$\square$

**系 1.4.** $A$ を固有値 $\lambda_i$ をもつ正規行列とする。このとき

1. $\|A\|=\max_i|\lambda_i|$；
2. $A$ が自己共役 $\iff$ すべての $\lambda_i\in\mathbb R$；ユニタリ $\iff$ すべての $|\lambda_i|=1$；
   射影 $\iff$ すべての $\lambda_i\in\{0,1\}$；
3. $A\ge0$ $\iff$ $A$ が自己共役かつすべての $\lambda_i\ge0$。

*証明.* $y=U^*x$ とおくと $\|Ax\|^2=\sum_i|\lambda_i|^2|y_i|^2$ かつ $\|y\|=\|x\|$
であり、これより (1) が従う。等号は $|\lambda_i|$ が最大となる $i$ について $x=Ue_i$ のとき成り立つ。(2) では $A^*=U\bar\Lambda U^*$、
$A^*A=U|\Lambda|^2U^*$、$A^2=U\Lambda^2U^*$ であり、各条件は $\Lambda$ について成分ごとに成り立つときかつそのときに限り成り立つ。(3) について、まず任意の正の $A$ は自己共役である。実際、すべての $x$ に対し
$\langle x,(A-A^*)x\rangle=\langle x,Ax\rangle-\overline{\langle x,Ax\rangle}=0$ であり、
すべての $x$ に対し $\langle x,Dx\rangle=0$ を満たす行列 $D$ は零である（偏極恒等式：
$4\langle x,Dy\rangle=\sum_{k=0}^3 i^{-k}\langle x+i^ky,D(x+i^ky)\rangle$）。すると
$\langle x,Ax\rangle=\sum_i\lambda_i|y_i|^2$ であり、これがすべての $x$ に対して $\ge0$ となるのはすべての $\lambda_i\ge0$ のときかつそのときに限る。$\square$

## 1.3 関数計算

$A=U\Lambda U^*$ を正規とし、$f$ を $\operatorname{spec}A$ 上の複素数値関数とする。
$$
f(A)=U\,\mathrm{diag}\big(f(\lambda_1),\dots,f(\lambda_n)\big)\,U^*.
$$
と定義する。

**定理 1.5.** (a) $f(A)$ は $U$ の選び方によらない。(b) 多項式 $p$ に対し、$p(A)$ は $A$ の通常の多項式である。(c) $f\mapsto f(A)$ は線形かつ乗法的であり、$\bar f(A)=f(A)^*$ かつ $\operatorname{spec}f(A)=f(\operatorname{spec}A)$ である。
(d) $B$ が $A$ と可換ならば、$B$ は $f(A)$ と可換である。(e) $g$ が
$f(\operatorname{spec}A)$ 上で定義されていれば、$g(f(A))=(g\circ f)(A)$ である。

*証明.* (b) $A^k=U\Lambda^kU^*$ であるから、線形性により $p(A)=Up(\Lambda)U^*$ である。(a) スペクトルは有限集合であるから、ラグランジュ補間により
$\operatorname{spec}A$ 上で $p=f$ となる多項式 $p$ が存在する。すると $A$ を対角化する *任意の* ユニタリ $U$ に対し $U\,f(\Lambda)\,U^*=U\,p(\Lambda)\,U^*=p(A)$ であり、$p(A)$ は $U$ を用いずに定義される。(c) は対角行列 $f(\Lambda)$ について成り立ち、$U$ による共役で保たれる。
$Uf(\Lambda)U^*$ のスペクトルは $f(\Lambda)$ の対角成分である。(d) $B$ は $A$ の任意の多項式と可換であり、特に $p(A)=f(A)$ と可換である。(e) $f(A)=Uf(\Lambda)U^*$ は正規であり同じ $U$ で対角化されるから、(a) により $g(f(A))=U\,g(f(\Lambda))\,U^*=(g\circ f)(A)$ である。$\square$

**定理 1.6（平方根）.** $A\ge0$ に対し、
$B^2=A$ を満たす $B\ge0$ はちょうど一つ存在し、それは関数計算で定義される $B=\sqrt A$ である。

*証明.* 定理 1.5(c) と系 1.4(3) により $\sqrt A\ge0$ かつ $(\sqrt A)^2=A$ である。
$B\ge0$ が $B^2=A$ を満たすとし、$B=V\,\mathrm{diag}(\mu_j)\,V^*$、$\mu_j\ge0$ とする。すると
$A=V\,\mathrm{diag}(\mu_j^2)\,V^*$ であるから $\operatorname{spec}A=\{\mu_j^2\}$ である。$p$ を
$\operatorname{spec}A$ 上で $p(t)=\sqrt t$ となる多項式とする。すると
$p(A)=\sqrt A$ かつ $p(A)=p(B^2)=V\,\mathrm{diag}(p(\mu_j^2))\,V^*=V\,\mathrm{diag}(\mu_j)\,V^*=B$ である。$\square$

**系 1.7.** $A\ge0$ $\iff$ ある $T$ により $A=T^*T$。（$A\ge0$ ならば $T=\sqrt A$ ととればよい。
逆に $\langle x,T^*Tx\rangle=\|Tx\|^2\ge0$ である。）

## 1.4 極分解と特異値

任意の $A$ に対し、行列 $A^*A$ は正である（系 1.7）。その平方根
$|A|=(A^*A)^{1/2}$ を $A$ の **絶対値** という。

**定理 1.8（極分解）.** 任意の $A\in M_n(\mathbb C)$ は、ユニタリ $W$ を用いて
$A=W|A|$ と書ける。

*証明.* 任意の $x$ に対し $\||A|x\|^2=\langle x,|A|^2x\rangle=\langle x,A^*Ax\rangle=\|Ax\|^2$ である。
したがって $W_0:|A|x\mapsto Ax$ は $|A|$ の値域
$R_1$ から $A$ の値域 $R_2$ の上への well-defined な線形等長写像である：$|A|x=|A|x'$ ならば
$\|A(x-x')\|=\||A|(x-x')\|=0$ である。特に $\dim R_1=\dim R_2$ であるから、直交補空間の次元も等しい。$W_1:R_1^\perp\to R_2^\perp$ をそれらの間の任意の等長写像とする（正規直交基底を正規直交基底に写せばよい）。すると
$W=W_0\oplus W_1$ はユニタリであり、$W|A|x=W_0|A|x=Ax$ である。$\square$

$|A|$ の固有値 $s_1\ge\dots\ge s_n\ge0$ を $A$ の **特異値** という。本講義を通じて三つのノルムを用いる。
$$
\|A\|=s_1,\qquad \|A\|_2=\Big(\sum_is_i^2\Big)^{1/2}=(\operatorname{Tr}A^*A)^{1/2},
\qquad\|A\|_1=\sum_is_i=\operatorname{Tr}|A|.
$$
第一の等式は、任意の $x$ に対し $\|Ax\|=\||A|x\|$ であるから系 1.4(1) により
$\|A\|=\||A|\|=s_1$ となることから成り立つ。第二の等式は
$\operatorname{Tr}A^*A=\operatorname{Tr}|A|^2=\sum_is_i^2$ による。

**定理 1.9（トレースノルムの双対性）.**
$\|A\|_1=\max_U|\operatorname{Tr}(AU)|$（$U$ はユニタリ行列全体を動く）であり、最大値は
$A=W|A|$ として $U=W^*$ で達成される。特に $|\operatorname{Tr}A|\le\|A\|_1$ である。

*証明.* $\operatorname{Tr}(AW^*)=\operatorname{Tr}(W|A|W^*)=\operatorname{Tr}|A|$ である。任意のユニタリ $U$ に対し、
トレースの巡回性と $\langle X,Y\rangle_{HS}=\operatorname{Tr}(X^*Y)$ に関するコーシー・シュワルツの不等式により
$$
|\operatorname{Tr}(AU)|=|\operatorname{Tr}(|A|^{1/2}\,UW\,|A|^{1/2})|
=|\langle|A|^{1/2},\,UW|A|^{1/2}\rangle_{HS}|
\le\||A|^{1/2}\|_2\,\|UW|A|^{1/2}\|_2=\operatorname{Tr}|A|,
$$
となる。ここでユニタリ $UW$ に対し $\|UW X\|_2=\|X\|_2$ であること、および $\||A|^{1/2}\|_2^2=\operatorname{Tr}|A|$ を用いた。$\square$

特に、上の三つの量はノルムである。$\|\cdot\|$ は作用素ノルムであり、
$\|\cdot\|_2$ は内積 $\langle X,Y\rangle_{HS}$ から定まるノルムである。
$\|\cdot\|_1$ については、定理 1.9 により
$\|A+B\|_1=\max_U|\operatorname{Tr}(AU)+\operatorname{Tr}(BU)|\le\|A\|_1+\|B\|_1$ であり、また
$\|A\|_1=0$ ならば $|A|=0$、したがって $A=W|A|=0$ である。$P\ge0$ に対しては定理 1.6 により $|P|=(P^2)^{1/2}=P$ であるから、
$\|P\|_1=\operatorname{Tr}P$ である。

## 1.5 ライブラリにおける作用素

ライブラリは環 $\bigoplus_cM_{k_c}(\mathbb C)$（第2章）を扱う。単一のブロックは $C=1$ の場合であり、`TypeIAlgebra([n], [n])` が $M_n(\mathbb C)$ である。（二つのリストはそれぞれ名目上のブロックサイズと実際のブロックサイズであり、当面は両者は等しい。）

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra

torch.manual_seed(0)
M3 = TypeIAlgebra([3], [3], precision='double')       # the algebra M_3(C), complex128

# a positive operator with eigenvalues 1, 2, 4 and Haar-random eigenvectors
A = M3.operator_from_eigenvalues(lambda d: torch.tensor([[1.0, 2.0, 4.0]]),
                                 force_positive=True)
print(A.matrix.shape)          # (batch, C, k_max, k_max) = (1, 1, 3, 3)
print(A.eigenvalues()[0])      # one tensor per block, eigenvalues in ascending order
```

各作用素は形状 `(batch, C, k_max, k_max)` の一つのテンソルとして格納される。すなわち作用素のバッチであって、各作用素は `C` 個のブロックのリストである。ここでは `batch = C = 1` である。この構成はまさに定理 1.3 を逆向きに読んだものである：$U$ をハール分布に従うユニタリとして $A=U\Lambda U^*$。フラグ `force_positive=True` は固有値に対して検査される（系 1.4(3)）ので、ライブラリは改めて検査することなく `A` が正であることを知っている。

自己共役作用素に対する関数計算は `apply_function` である（各ブロックを `eigh` で対角化する）。関数 $f$ は要素ごとに作用しなければならない：同じサイズ $k$ のすべてのセクターの固有値を、形 `(batch, n·k)` の一つのテンソルとして受け取り（$n$ はそのようなセクターの個数）、同じ形のテンソルを返す。よく使う関数には専用のメソッドがある。定理 1.5 と 1.6 を検証しよう。

```python
def close(X, Y, tol=1e-10):
    """Maximal entry of |X - Y| for two operators (or tensors) is below tol."""
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

R = A.sqrt()
assert close(R @ R, A)                                  # Theorem 1.6
assert close(A.log().expm(), A)                         # exp(log A) = A, Theorem 1.5(e)
assert close(A.apply_function(lambda w: w ** 3), A @ A @ A)   # Theorem 1.5(b)
assert close(A.power(-1), A.inverse())                  # f(t) = 1/t
spec = A.apply_function(torch.exp).eigenvalues()[0]     # spec f(A) = f(spec A)
assert torch.allclose(spec, torch.exp(torch.tensor([[1.0, 2.0, 4.0]], dtype=spec.dtype)))
```

`R @ R` はこの時点ではまだ何も計算しない。積、和、関数は計算の手順（レシピ）を組み立てるだけであり、`.matrix` やトレースなどの数値が要求されたときに評価される。

## 1.6 ノルムとトレースノルムの双対性

一般の（正規とは限らない）行列は `alg.operator(tensor)` によって導入する。三つのノルムを PyTorch で直接計算した特異値と照合し、さらに極分解（定理 1.8）と双対性（定理 1.9）を検証する。

```python
G = torch.randn(1, 1, 3, 3, dtype=torch.complex128)    # Ginibre matrix: no symmetry at all
X = M3.operator(G)
s = torch.linalg.svdvals(G[0, 0])                       # singular values, descending

assert abs(X.operator_norm().item() - s[0].item()) < 1e-12
assert abs(X.frobenius_norm().item() - s.norm().item()) < 1e-12
assert abs(X.trace_norm().item() - s.sum().item()) < 1e-12

absX = X.abs()
W = X @ absX.inverse()                                  # |X| is invertible (almost surely, for Ginibre)
assert close(W.adjoint() @ W, M3.identity())            # W is unitary
assert abs((X @ W.adjoint()).trace.item() - s.sum().item()) < 1e-10    # Tr(X W*) = ||X||_1
assert abs(X.trace.item()) <= X.trace_norm().item()     # |Tr X| <= ||X||_1
```

## 1.7 バッチ：多数の作用素を一度に

バッチ軸こそがこのライブラリの存在理由である。ランダム行列理論、モンテカルロ法、マルチスタート最適化には数千個の作用素が必要であり、バッチ化された一回の呼び出しは Python のループよりはるかに高速である。例として、**ガウス型ユニタリアンサンブル**（GUE）の厳密なモーメントを検証する。

$G$ を独立な成分をもつ行列で、各成分の実部と虚部は独立に
$\sim N(0,\tfrac12)$ に従うものとする。したがって $\mathbb E|G_{ij}|^2=1$ である。$H=(G+G^*)/\sqrt2$ とおく。

**命題 1.10.** $\mathbb E\,\operatorname{Tr}H^4=2n^3+n$ である。したがって
$X=H/\sqrt n$ に対し $\ \mathbb E\,\tfrac1n\operatorname{Tr}X^4=2+\tfrac1{n^2}$ である。

*証明.* *イッサーリス（ウィック）の公式.* $Y_1,\dots,Y_4$ を同時ガウス型で中心化された実確率変数（重複を許す）とし、$c_{ij}=\mathbb E[Y_iY_j]$ とおく。それらの同時積率母関数は
$q(t)=\sum_{ij}c_{ij}t_it_j$ として $\mathbb E\,e^{\sum_it_iY_i}=e^{q(t)/2}$ である（ガウス分布の積率母関数であり、確率論の標準的事実である）。期待値の下での微分は正当化される：積率母関数はすべての $t$ で有限であり、$|t_i|\le1$ では $e^{\sum_it_iY_i}$ の4階までのすべての偏導関数が $(1+\sum_i|Y_i|)^4e^{\sum_i|Y_i|}$ で抑えられ、$e^{\sum_i|Y_i|}\le\sum_{s\in\{\pm1\}^4}e^{\sum_is_iY_i}$ よりこれは可積分である。したがって優収束定理が適用できる。$t=0$ において $\partial^4/\partial t_1\cdots\partial t_4$ を作用させる。左辺からは
$\mathbb E[Y_1Y_2Y_3Y_4]$ が得られる。右辺では指数級数のうち項 $q^2/8$ のみが $t_1t_2t_3t_4$ を含み、その係数は $c_{12}c_{34}+c_{13}c_{24}+c_{14}c_{23}$ である。
したがって $\mathbb E[Y_1Y_2Y_3Y_4]$ は、三通りの対合わせ（ペアリング）にわたる、対ごとの期待値の積の和である。多重線形性により、同時ガウス型の実確率変数の複素線形結合である複素確率変数（$H$ の成分など）についても同じことが成り立つ。

*モーメントの計算.* 対ごとの期待値は
$\mathbb E[H_{ij}H_{kl}]=\tfrac12\mathbb E[(G_{ij}+\bar G_{ji})(G_{kl}+\bar G_{lk})]=\delta_{il}\delta_{jk}$
である。実際、$\mathbb E[G_{ij}G_{kl}]=0$ かつ $\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}$ である。
$\operatorname{Tr}H^4=\sum_{ijkl}H_{ij}H_{jk}H_{kl}H_{li}$ において、各対合わせは次を与える：
$\{12\}\{34\}$：$\delta_{ik}\delta_{ki}$、和は $n^3$（自由な添字は $i=k$, $j$, $l$）；
$\{14\}\{23\}$：$\delta_{jl}\delta_{jl}$、和は $n^3$；
$\{13\}\{24\}$：$\delta_{il}\delta_{jk}\delta_{ji}\delta_{kl}$、これは $i=j=k=l$ を強制し、和は $n$。
$n\cdot n^2$ で割ると第二の主張を得る。$\square$

二つの *平面的* な対合わせは半円則のカタラン数 $2$ を与え、交差する対合わせは $1/n^2$ の補正を与える。

```python
n, B = 4, 200_000
G = torch.complex(torch.randn(B, 1, n, n), torch.randn(B, 1, n, n)).to(torch.complex128) / math.sqrt(2)
M4 = TypeIAlgebra([n], [n], precision='double')
X = M4.operator((G + G.conj().transpose(-2, -1)) / math.sqrt(2 * n), is_self_adjoint=True)

X2 = X @ X
m4 = (X2 @ X2).trace.real / n                           # shape (B,): one value per sample
mean, err = m4.mean().item(), m4.std().item() / math.sqrt(B)
print(f"E tr X^4 / n = {mean:.4f} +- {err:.4f}   (exact {2 + 1 / n**2:.4f})")
assert abs(mean - (2 + 1 / n**2)) < 5 * err
```

モンテカルロ法による検証は統計的検定である：推定値は厳密値と標準誤差の数倍以内で一致しなければならない。これが本チュートリアルにおける既知の結果との照合すべてに共通する型である：厳密な主張、数値的推定、そして明示的な許容誤差。

## 1.8 比較の対象

単一の行列については、`numpy.linalg` と `scipy.linalg`（`sqrtm`, `expm`, `logm`,
`polar`）が参照実装であり、上の検証はそれらを用いても繰り返すことができる。本ライブラリは、後で重要になる三つのものを付け加える：バッチ化（1.7節）、ブロック構造を担う作用素（第2章）、そしてすべての演算を通じた自動微分である。

## 演習問題

1. ランダムな自己共役行列 $A$ について、$f=\sin$, $g=\exp$ に対し $(fg)(A)=f(A)g(A)$ を確かめよ。
2. 定理 1.9 を用いて、$\|A\|_1=\|A^*\|_1$ および $\|AB\|_1\le\|A\|\,\|B\|_1$ を示せ。
3. $A,B>0$ の **幾何平均** を $A\#B=A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}$ と定める。
   $A\#B$ が $XA^{-1}X=B$ のただ一つの正の解 $X$ であることを証明し、そこから
   $A\#B=B\#A$ を導け。また両者を数値的に確かめよ。

<details>
<summary>解答</summary>

1. 定理 1.5(c) による。

```python
Hs = M3.operator_from_eigenvalues(lambda d: torch.randn(1, d), force_self_adjoint=True)
assert close(Hs.apply_function(lambda w: torch.sin(w) * torch.exp(w)),
             Hs.apply_function(torch.sin) @ Hs.apply_function(torch.exp))
```

2. $|\operatorname{Tr}(A^*U)|=|\overline{\operatorname{Tr}(U^*A)}|=|\operatorname{Tr}(AU^*)|$ であり、
   $U\mapsto U^*$ はユニタリ群の全単射であるから、定理 1.9 における最大値は一致する。第二の主張について、$B=V|B|$ とする。ユニタリ $U$ に対し
   $|\operatorname{Tr}(ABU)|=|\operatorname{Tr}(|B|^{1/2}UAV|B|^{1/2})|
   \le\||B|^{1/2}\|_2\,\|UAV|B|^{1/2}\|_2\le\|A\|\operatorname{Tr}|B|$ である。ここで
   $\|YX\|_2\le\|Y\|\,\|X\|_2$（$X$ の各列に $\|Yx\|\le\|Y\|\,\|x\|$ を適用する）および
   ユニタリ $U,V$ に対する $\|UAV\|=\|A\|$（ユニタリはベクトルの長さを保つ）を用いた。

3. $T=A^{-1/2}XA^{-1/2}$ とおく。すると $XA^{-1}X=B$ $\iff$ $T^2=A^{-1/2}BA^{-1/2}$ であり、
   $\langle x,Xx\rangle=\langle A^{1/2}x,TA^{1/2}x\rangle$ かつ
   $A^{1/2}$ は可逆であるから、$X\ge0$ $\iff$ $T\ge0$ である。定理 1.6 により $T=(A^{-1/2}BA^{-1/2})^{1/2}$ はただ一つの正の解であり、すなわち $X=A\#B$ である。$XA^{-1}X=B$ の任意の解 $X$ は可逆である。
   実際 $\det(X)^2\det(A)^{-1}=\det B\neq0$ である。よって $XA^{-1}X=B$ は
   $XB^{-1}X=A$ と同値であり（$A^{-1}=X^{-1}BX^{-1}$ の両辺の逆をとる）、一意性から
   $A\#B=B\#A$ を得る。

```python
def geometric_mean(A, B):
    a, ai = A.sqrt(), A.power(-0.5)
    return a @ (ai @ B @ ai).sqrt() @ a

P = M3.operator_from_eigenvalues(lambda d: 0.5 + torch.rand(1, d), force_positive=True)
Q = M3.operator_from_eigenvalues(lambda d: 0.5 + torch.rand(1, d), force_positive=True)
GM = geometric_mean(P, Q)
assert close(GM, geometric_mean(Q, P), 1e-9)
assert close(GM @ P.inverse() @ GM, Q, 1e-9)
```

</details>

**次へ：**[第2章](02_von_neumann_algebras.md) では $M_n$ を行列環の直和に置き換え、すべての有限次元フォン・ノイマン環がこの形であることを証明する。
