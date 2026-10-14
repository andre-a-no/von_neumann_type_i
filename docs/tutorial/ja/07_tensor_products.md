[English](../07_tensor_products.md) · [Español](../es/07_tensor_products.md) · [Français](../fr/07_tensor_products.md) · [Deutsch](../de/07_tensor_products.md) · [中文](../zh/07_tensor_products.md) · **日本語** · [Русский](../ru/07_tensor_products.md)

# 7. テンソル積とエンタングルメント

本章では、セクターをもつ二つの環のテンソル積を調べる（物理的には、二つの系を合わせたものの環である。*物理的意味* の囲みを参照）。セクターがある場合には二つの積がある：各部分がそれぞれ自分の電荷を保つ通常の積と、全電荷だけが保存される *融合* 積である。両者の構造を証明し、部分トレースを埋め込みの双対として定義し、エンタングルメントを調べる。

**本章で学ぶこと**：セクターがどのように掛け合わされ、ライブラリがそれらをどう並べるか、部分トレースがなぜチャネルであり、正規化すれば有限指数の条件付き期待値であるか、電荷を保存する状態の縮約状態は何か、エンタングルメントが粒子数の部分と配置の部分にどう分かれるか、そしてランダムな状態の平均純度と平均エントロピー。

本章を通じて、$\mathcal M_1=\bigoplus_{c=1}^{C_1}M_{k_c}(\mathbb C)$ は
$H_1=\bigoplus_c\mathbb C^{k_c}$ に極小中心射影 $P_c$ をもって作用し、
$\mathcal M_2=\bigoplus_{d=1}^{C_2}M_{l_d}(\mathbb C)$ は $H_2=\bigoplus_d\mathbb C^{l_d}$ に
極小中心射影 $Q_d$ をもって作用する。$\operatorname{Tr}$ は作用素が作用する空間のトレースであり、$\mathcal M_i$ 上では第3章の blunt トレースである。
$\mathbb C^k\otimes\mathbb C^l$ 上の $X$ に対し $X_{(a,b),(a',b')}=\langle e_a\otimes e_b,X\,e_{a'}\otimes e_{b'}\rangle$ と書き、
$e_a\otimes e_b\mapsto e_{(a-1)l+b}$ によって $\mathbb C^k\otimes\mathbb C^l$ を $\mathbb C^{kl}$ と同一視する
（クロネッカー規約）。したがって $(A\otimes B)_{(a,b),(a',b')}=A_{aa'}B_{bb'}$ である。

## 7.1 二つの環のテンソル積

$H_1$ 上の $*$-環 $\mathcal M_1$ と $H_2$ 上の $*$-環 $\mathcal M_2$ に対し、**テンソル積**
$\mathcal M_1\otimes\mathcal M_2$ とは、$H_1\otimes H_2$ 上の作用素 $A\otimes B$（$A\in\mathcal M_1$, $B\in\mathcal M_2$）の線形包である。
$(A\otimes B)(A'\otimes B')=AA'\otimes BB'$、$(A\otimes B)^*=A^*\otimes B^*$ かつ $1\otimes1=1$ であるから、
これは $*$-環であり、したがってフォン・ノイマン環である（定理 2.2）。

**定理 7.1.** (a) $H_1\otimes H_2=\bigoplus_{c,d}\mathbb C^{k_c}\otimes\mathbb C^{l_d}$ は射影
$P_c\otimes Q_d$ の値域の直交和であり、
$$
\mathcal M_1\otimes\mathcal M_2=\bigoplus_{c,d}M_{k_c}\otimes M_{l_d}\cong\bigoplus_{c,d}M_{k_cl_d}(\mathbb C):
$$
これは $\mathbb C^{k_c}\otimes\mathbb C^{l_d}$ 上のブロックが任意である作用素 $X=\sum_{c,d}(P_c\otimes Q_d)X(P_c\otimes Q_d)$ 全体からなる。
(b) その極小中心射影は $C_1C_2$ 個の射影 $P_c\otimes Q_d$ であり、$\dim(\mathcal M_1\otimes\mathcal M_2)=\dim\mathcal M_1\dim\mathcal M_2$ である。
(c) $\mathcal M_1$ 上のトレース $\operatorname{tr}_w$ と $\mathcal M_2$ 上のトレース $\operatorname{tr}_v$（定理 3.1）に対し、
$(w\otimes v)_{cd}=w_cv_d$ として $\operatorname{tr}_w(A)\operatorname{tr}_v(B)=\operatorname{tr}_{w\otimes v}(A\otimes B)$ である。
特に、積の `Tr_blunt`、`Tr_norm`、`tau_vN` は、因子の同じトレースの積である。

*証明.* (a) $P_c\otimes Q_d$ は互いに直交する射影で和は $1\otimes1$ であり、
$P_c\otimes Q_d$ の値域は $\mathbb C^{k_c}\otimes\mathbb C^{l_d}$ である。$P_c$, $Q_d$ は中心元であるから
$(P_c\otimes Q_d)(A\otimes B)=AP_c\otimes BQ_d=(A\otimes B)(P_c\otimes Q_d)$ となり、
$\mathcal M_1\otimes\mathcal M_2$ の任意の元はブロック対角である。逆に、ブロック $c$ と $d$ の行列単位
$E^c_{aa'}\in\mathcal M_1$ と $F^d_{bb'}\in\mathcal M_2$ を用いると、
$\mathbb C^{k_c}\otimes\mathbb C^{l_d}$ の内部で $E^c_{aa'}\otimes F^d_{bb'}=|e_a\otimes e_b\rangle\langle e_{a'}\otimes e_{b'}|$ である。
これらは $M_{k_cl_d}(\mathbb C)$ の行列単位であるから、ブロック $(c,d)$ に台をもつ任意の作用素は
$\mathcal M_1\otimes\mathcal M_2$ に属する。(b) 全行列環の直和の中心は各直和成分への射影で張られる（第2章、演習問題 1）。また
$\sum_{c,d}k_c^2l_d^2=\sum_ck_c^2\sum_dl_d^2$ である。(c)
$\operatorname{Tr}(X\otimes Y)=\sum_{a,b}X_{aa}Y_{bb}=\operatorname{Tr}X\operatorname{Tr}Y$ であるから
$\operatorname{tr}_{w\otimes v}(A\otimes B)=\sum_{c,d}w_cv_d\operatorname{Tr}A_c\operatorname{Tr}B_d=\operatorname{tr}_w(A)\operatorname{tr}_v(B)$ である。
三つのトレースの重み $1$, $1/k_c$, $1/(C_1k_c)$（3.1節）を掛け合わせると $1$,
$1/(k_cl_d)$, $1/(C_1C_2k_cl_d)$ となり、これは積の上の同じトレースの重みである。$\square$

`tensor_product(alg1, alg2)` はこの環を構成する。セクター $(c,d)$ の番号は $s=cC_2+d$（コードでは番号は
$0$ から始まる）、電荷は $(q_c,r_d)$ であり、その内部で $e_a\otimes e_b$ の番号は $a\,l_d+b$ である。
`kron(A, B, alg12)` は $A\otimes B$ を計算する。セクターの順序は $H_1\otimes H_2$ のクロネッカー順序と基底の置換だけ異なる。この置換は下の `kron_rows` が明示的に構成し、それを用いて全行列のクロネッカー積と比較する。

```python
import math
import numpy as np
import torch
from torch_vn_algebra import (TypeIAlgebra, DensityMatrix, states, composite,
                              tensor_product, fused_tensor_product, kron, partial_trace)

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def full(X):
    """An operator of (+)_c M_{k_c} as a block-diagonal matrix on (+)_c C^{k_c}, shape (batch, D, D)."""
    return torch.stack([torch.block_diag(*[x[c, :k, :k] for c, k in enumerate(X.algebra.k_factors)])
                        for x in X.matrix])

def kron_rows(alg12):
    """For every sector of a product algebra, the positions of its basis vectors in H1 (x) H2 (Kronecker order)."""
    k1, k2 = alg12.tensor_factors[0].k_factors, alg12.tensor_factors[1].k_factors
    return [torch.tensor([(sum(k1[:c]) + a) * sum(k2) + sum(k2[:d]) + b
                          for c, d, off in pairs for a in range(k1[c]) for b in range(k2[d])])
            for pairs in alg12.product_blocks]              # pair blocks (c, d, offset), by offset

def full12(X):
    """An operator of a (plain or fused) product algebra as a matrix on H1 (x) H2 in Kronecker order."""
    f1, f2 = X.algebra.tensor_factors
    N = sum(f1.k_factors) * sum(f2.k_factors)
    out = torch.zeros(X.matrix.shape[0], N, N, dtype=X.matrix.dtype)
    for s, r in enumerate(kron_rows(X.algebra)):
        out[:, r[:, None], r[None, :]] = X.matrix[:, s, :len(r), :len(r)]
    return out

def random_operator(alg):
    return alg.from_blocks([torch.randn(k, k, dtype=dt) for k in alg.k_factors])

alg1 = TypeIAlgebra([1, 2], [1, 2], precision='double')            # C (+) M_2
alg2 = TypeIAlgebra([2, 3], [2, 3], precision='double')            # M_2 (+) M_3
alg12 = tensor_product(alg1, alg2)
print(alg12.k_factors, alg12.charges, sum(k * k for k in alg12.k_factors))     # 4 sectors, dimension 65
assert sum(k * k for k in alg12.k_factors) == 5 * 13                           # Theorem 7.1(b)

A, A2, B, B2 = random_operator(alg1), random_operator(alg1), random_operator(alg2), random_operator(alg2)
AB = kron(A, B, alg12)
assert close(full12(AB)[0], torch.kron(full(A)[0], full(B)[0]))               # Theorem 7.1(a)
assert close(AB @ kron(A2, B2, alg12), kron(A @ A2, B @ B2, alg12))
for tr in ('Tr_blunt', 'Tr_norm', 'tau_vN'):                                   # Theorem 7.1(c)
    assert (getattr(AB, tr)() - getattr(A, tr)() * getattr(B, tr)()).abs().max() < 1e-10
```

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. 二つの系がそれぞれ自分の超選択則に従うとき――粒子を交換できない二つの箱――、結合系は
> 電荷の *組* ごとに一つのセクターをもち、異なる組を結ぶオブザーバブルは存在しない。結合セクターは重みの積で数えられる
> （定理 7.1(c)）。

## 7.2 部分トレース

**定義.** $\rho\in\mathcal M_1\otimes\mathcal M_2$ とする。**部分トレース**
$\operatorname{Tr}_2\rho\in\mathcal M_1$ とは、次を満たす元である。
$$
\operatorname{Tr}\big((\operatorname{Tr}_2\rho)\,A\big)=\operatorname{Tr}\big(\rho\,(A\otimes1)\big)\qquad\text{for all }A\in\mathcal M_1 .
$$
線形汎関数 $A\mapsto\operatorname{Tr}(\rho(A\otimes1))$ と $\mathcal M_1$ の忠実なトレース $\operatorname{Tr}$ に定理 3.2 を適用すれば、
これはただ一つ存在する：部分トレースは、blunt トレースに関する埋め込み $A\mapsto A\otimes1$ の双対である。
$\operatorname{Tr}_1$ も同様に定義する。状態 $\rho$ に対し、$\operatorname{Tr}_2\rho$ と
$\operatorname{Tr}_1\rho$ を **縮約状態** という。全行列環 $\mathcal M_i=B(H_i)$ に対しては、これは
**ヒルベルト空間の部分トレース** $\operatorname{Tr}_{H_2}:B(H_1\otimes H_2)\to B(H_1)$ である。

**命題 7.2.** (a) $\rho=\bigoplus_{c,d}\rho_{cd}$ に対し
$$
(\operatorname{Tr}_2\rho)_c=\sum_d\operatorname{tr}_{l_d}\rho_{cd},\qquad
(\operatorname{tr}_lX)_{aa'}=\sum_{b=1}^lX_{(a,b),(a',b)}\quad(X\text{ on }\mathbb C^k\otimes\mathbb C^l),
$$
であり、$\operatorname{Tr}_2\rho=\operatorname{Tr}_{H_2}\rho$ である。ここで $(\operatorname{Tr}_{H_2}X)_{ii'}=\sum_jX_{(i,j),(i',j)}$ である。
(b) $\operatorname{Tr}_2(A\otimes B)=(\operatorname{Tr}B)\,A$。
(c) $A,A'\in\mathcal M_1$ に対し $\operatorname{Tr}_2\big((A\otimes1)X(A'\otimes1)\big)=A\,(\operatorname{Tr}_2X)\,A'$。
(d) $\operatorname{Tr}(\operatorname{Tr}_2X)=\operatorname{Tr}X$ であり、$\operatorname{Tr}_2$ は状態を状態に写す。

*証明.* (a) $A=\bigoplus_cA_c$ とすると $A\otimes1=\bigoplus_{c,d}A_c\otimes1_{l_d}$ であり、
$$
\operatorname{Tr}\big(\rho(A\otimes1)\big)=\sum_{c,d}\sum_{a,a',b}(\rho_{cd})_{(a,b),(a',b)}(A_c)_{a'a}
=\sum_c\operatorname{Tr}\Big(\Big(\sum_d\operatorname{tr}_{l_d}\rho_{cd}\Big)A_c\Big),
$$
であるから、右辺の元は定義を満たし、$\operatorname{Tr}_2\rho$ である。
$B(H_i)$ について同じ計算をすれば $\operatorname{Tr}_{H_2}$ の公式を得る。ブロック対角な $\rho$ では、$i,i'$ が異なるブロックにある成分
$\rho_{(i,j),(i',j)}$ は零であるから、$\operatorname{Tr}_{H_2}\rho$ は (a) のブロックをもつ。(b) $\operatorname{tr}_l(X\otimes Y)=(\operatorname{Tr}Y)X$。(c) $B\in\mathcal M_1$ に対し
$\operatorname{Tr}\big(\operatorname{Tr}_2((A\otimes1)X(A'\otimes1))B\big)=\operatorname{Tr}\big(X(A'BA\otimes1)\big)
=\operatorname{Tr}\big(A(\operatorname{Tr}_2X)A'B\big)$ である。あとは一意性を用いる。(d) $A=1$ ととる。$\rho$ が状態ならば、
$A^*A\otimes1=(A\otimes1)^*(A\otimes1)$ であるから $A\mapsto\operatorname{Tr}(\rho(A\otimes1))$ は正であり、
$A=1$ で $1$ をとるので、その密度 $\operatorname{Tr}_2\rho$ は状態である（定理 3.2）。$\square$

チャネルは第6章の主題である。ここでは次の定義だけが必要である：線形写像
$\Phi$ が **完全正** であるとは、すべての $n$ に対し $\mathrm{id}_n\otimes\Phi$ が
$M_n(\mathbb C)\otimes(\cdot)$ 上で正であることをいい、さらにトレースを保つとき **チャネル** という。

**定理 7.3.** $j$ を $H_2$ の標準基底にわたって動かし、$K_j=1\otimes\langle j|:H_1\otimes H_2\to H_1$,
$K_j(x\otimes y)=\langle j,y\rangle x$ とする。このとき
$$
\operatorname{Tr}_{H_2}X=\sum_jK_jXK_j^*,\qquad\sum_jK_j^*K_j=1,\qquad\sum_jK_j^*AK_j=A\otimes1 .
$$
である。したがって $\operatorname{Tr}_2$ は $\mathcal M_1\otimes\mathcal M_2$ から $\mathcal M_1$ へのチャネルであり、
その blunt トレースに関する双対は $A\mapsto A\otimes1$ である。

*証明.* $K_j^*x=x\otimes e_j$ であるから $\langle e_a,K_jXK_j^*e_{a'}\rangle=X_{(a,j),(a',j)}$ であり、$j$ について和をとると
命題 7.2(a) を得る。さらに $\sum_jK_j^*K_j(x\otimes y)=\sum_j\langle j,y\rangle x\otimes e_j=x\otimes y$
かつ $\sum_jK_j^*AK_j=\sum_jA\otimes|j\rangle\langle j|=A\otimes1$ である。$Y\ge0$ に対し、
$\langle x,TYT^*x\rangle=\langle T^*x,YT^*x\rangle\ge0$ であるから
$(\mathrm{id}_n\otimes\operatorname{Tr}_{H_2})(Y)=\sum_j(1_n\otimes K_j)Y(1_n\otimes K_j)^*\ge0$ である。これは部分環
$\mathcal M_1\otimes\mathcal M_2$ の上でもそのまま成り立つ。最後に $\operatorname{Tr}\sum_jK_jXK_j^*A=\operatorname{Tr}X\sum_jK_j^*AK_j$ であり、
これは $A=1$ のときトレースの保存を、一般には双対性を与える。$\square$

`partial_trace(X, keep=1)` は各ブロック $\rho_{cd}$ を $(k_c,l_d,k_c,l_d)$ の形に並べ替えてトレースをとる。
`composite.partial_trace_channel(alg12, keep=1)` は定理 7.3 のチャネルであり、そのクラウス作用素は
セクター $(c,d)$ からセクター $c$ への $1_{k_c}\otimes\langle j|$ である。両者を定義、独立な `numpy.einsum`、
および双対性と照合する。

```python
rho = states.random_density_matrix(alg12, batch_size=3)
rho1, rho2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
assert isinstance(rho1, DensityMatrix)                                      # states give states
R = full12(rho).numpy().reshape(3, 3, 5, 3, 5)                              # indices (batch, i, j, i', j')
assert np.abs(full(rho1).numpy() - np.einsum('xijkj->xik', R)).max() < 1e-12
assert np.abs(full(rho2).numpy() - np.einsum('xijil->xjl', R)).max() < 1e-12
one2 = alg2.identity()
assert (states.expectation(rho, kron(A, one2, alg12)) - states.expectation(rho1, A)).abs().max() < 1e-12
assert close(partial_trace(AB, keep=1), A * B.trace.item())                 # Proposition 7.2(b)

Ptr = composite.partial_trace_channel(alg12, keep=1)
print(Ptr)
assert Ptr.is_trace_preserving()
assert close(Ptr(rho), rho1)
assert close(Ptr.adjoint('blunt')(A), kron(A, one2, alg12))                 # Theorem 7.3: the dual
```

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. 縮約状態 $\operatorname{Tr}_2\rho$ は、第一の系しか測定できない観測者が知りうるすべてであり、
> 定義により、それはオブザーバブル $A\otimes1$ のあらゆる期待値を再現する。縮約状態が等しい結合状態どうしは、相関において異なる。

**条件付き期待値としての部分トレース.** $\mathcal M_2=M_m(\mathbb C)$ とすると、
$\mathcal M_1\otimes M_m$ のセクターは $\mathcal M_1$ のセクターと同じである（定理 7.1）。
$\iota(A)=A\otimes1_m$、$\mathcal E=\frac1m\operatorname{Tr}_2:\mathcal M_1\otimes M_m\to\mathcal M_1$ とおき、
$e_{ab}$ を $M_m$ の行列単位とする。

**命題 7.4.** (a) $\iota\circ\mathcal E$ は $\mathcal M_1\otimes M_m$ から $\mathcal M_1\otimes1$ の上への条件付き期待値であり
（2.6節）、$\mathcal M_1\otimes M_m$ の任意の忠実なトレース $\operatorname{tr}_w$ に対し、$\operatorname{tr}_w$ を保つものはこれただ一つである。
(b) $m^2$ 個の元 $u_{ab}=\sqrt m\,(1\otimes e_{ab})$ は、すべての $X$ に対し $X=\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*X))$ を満たし、
$\sum_{a,b}u_{ab}u_{ab}^*=m^2\,1$ である。(c) すべての $X$ に対し
$X=\sum_kv_k\,\iota(\mathcal E(v_k^*X))$ を満たす任意の族 $(v_k)$ は $\sum_kv_kv_k^*=m^2\,1$ を満たす。

*証明.* (a) $\mathcal E$ は正であり（定理 7.3）、$\iota$ も正である
（$T^*T\otimes1=(T\otimes1)^*(T\otimes1)$）。命題 7.2(b) により $\mathcal E(A\otimes1)=A$ であり、
双加群性は命題 7.2(c) である。$\mathcal M_1$ 上の同じ重み $w_c$ をもつトレースも $\operatorname{tr}_w$ と書くと、
$\operatorname{tr}_w(Y\otimes1)=m\operatorname{tr}_w(Y)$ である。$\operatorname{Tr}_2$ の定義を $AP_c$ に適用し、$w_c$ を掛けて $c$ について和をとると
$\operatorname{tr}_w(X(A\otimes1))=\operatorname{tr}_w((\operatorname{Tr}_2X)A)$ を得て、$A=1$ とすると
$\operatorname{tr}_w(\iota\mathcal E(X))=\operatorname{tr}_w(X)$ を得る。$F=\iota\circ F_0$ が $\operatorname{tr}_w$ を保つ $\mathcal M_1\otimes1$ の上への別の条件付き期待値ならば、$F(X(A\otimes1))=F(X)(A\otimes1)$ であるから
$\operatorname{tr}_w(X(A\otimes1))=\operatorname{tr}_w(F(X)(A\otimes1))=m\operatorname{tr}_w(F_0(X)A)$ であるから、
すべての $A\in\mathcal M_1$ に対し $\operatorname{tr}_w((\operatorname{Tr}_2X-mF_0(X))A)=0$ となり、忠実性により $F_0=\mathcal E$ である
（$A=(\operatorname{Tr}_2X-mF_0(X))^*$ ととる）。(b) $X_{ab}\in\mathcal M_1$ により $X=\sum_{a,b}X_{ab}\otimes e_{ab}$ と書く
（定理 7.1(a)）。命題 7.2(b) により
$\mathcal E(u_{ab}^*X)=\frac{\sqrt m}m\sum_{b'}\operatorname{Tr}_2(X_{ab'}\otimes e_{bb'})=X_{ab}/\sqrt m$ であるから
$u_{ab}\,\iota(\mathcal E(u_{ab}^*X))=X_{ab}\otimes e_{ab}$ であり、また
$\sum_{a,b}u_{ab}u_{ab}^*=m\sum_{a,b}1\otimes e_{aa}=m^21$ である。(c) 命題 7.2(a) により $\operatorname{Tr}_2(Y^*)=(\operatorname{Tr}_2Y)^*$ であるから、
$Y=u_{ab}$ の展開の随伴は
$u_{ab}^*=\sum_k\iota(\mathcal E(u_{ab}^*v_k))\,v_k^*$ であり、(b) により
$\sum_{a,b}u_{ab}u_{ab}^*=\sum_k\big(\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*v_k))\big)v_k^*=\sum_kv_kv_k^*$ である。$\square$

(c) のような族を $\mathcal E$ の *準基底* といい、(c) によりその選び方によらず (b) により存在する $\sum_kv_kv_k^*$ が
$\mathcal E$ の **指数** $\operatorname{Ind}\mathcal E=m^2$ である（これは定義であり、名称は Y. Watatani, Mem. Amer. Math.
Soc. 83, no. 424, 1990 にならう）。$\bigoplus_cM_{k_c}$ の $\bigoplus_sM_{n_s}$ への包含の **包含行列** は、成分
$\Lambda_{sc}$ として、部分環のブロック $c$ が大きい環のブロック $s$ の中に現れるコピーの個数をもつ
（第9章で環の塔に用いる）。因子の場合、$M_k\otimes1\subset M_k\otimes M_m$ である：$\mathbb C^{km}$ の基底
$e_a\otimes e_b$ をまず $b$ で並べると、$A\otimes1_m=\operatorname{diag}(A,\dots,A)$（$m$ 個のコピー）となるので
$\Lambda=(m)$ であり、$m^2=\|\Lambda\|^2=\dim(M_k\otimes M_m)/\dim M_k$ である。（発展的な文献で、ここでは用いない：
部分因子の指数については V. F. R. Jones, Invent. Math. 72, 1, 1983、包含行列については
F. M. Goodman, P. de la Harpe and V. F. R. Jones, *Coxeter graphs and towers of algebras*, Springer,
1989。）一般の $\mathcal M_1$ に対しては、包含は
$C_1$ 個の包含 $M_{k_c}\otimes1\subset M_{k_c}\otimes M_m$ の直和であり、$\Lambda=m1_{C_1}$ である。モジュール
`torch_vn_algebra.inclusions` は包含とその条件付き期待値および指数を実装している。
`from_tensor_product(alg1, m)` はこの包含を構成する。

```python
from torch_vn_algebra import inclusions
m = 3
inc = inclusions.from_tensor_product(alg1, m)                # alg1 (x) 1 c alg1 (x) M_3
X = random_operator(inc.M)
for trace in ('Tr_blunt', 'Tr_norm', 'tau_vN'):
    E = inc.conditional_expectation(trace.replace('Tr_', ''))
    assert close(E(X), partial_trace(X, keep=1) * (1 / m))   # Proposition 7.4(a): the same E for every trace
    assert abs(getattr(inc.embed(E(X)), trace)() - getattr(X, trace)()).max() < 1e-10
    assert close(E.index(), inc.M.identity() * m**2)         # Proposition 7.4(c)

def unit(a, b):
    e = torch.zeros(m, m, dtype=dt)
    e[a, b] = 1
    return math.sqrt(m) * kron(alg1.identity(), inc.M.tensor_factors[1].from_blocks([e]), inc.M)

u = [unit(a, b) for a in range(m) for b in range(m)]         # the quasi-basis of Proposition 7.4(b)
assert close(sum((v @ inc.embed(E(v.adjoint() @ X))).matrix for v in u), X)
assert close(sum((v @ v.adjoint()).matrix for v in u), inc.M.identity() * m**2)
print(inc.matrix.tolist(), inc.component_indices())           # Lambda = 3 * 1_2, index 9 per component
```

## 7.3 電荷の融合

いま $\mathcal M_1$ のセクター $c$ が電荷 $q_c$ を、$\mathcal M_2$ のセクター $d$ が電荷
$r_d$ をもち、**融合則** $f$ が全電荷 $f(q_c,r_d)$ を与えるとする――粒子数なら加法、フェルミオン・パリティなら $2$ を法とする加法、あるいは他の任意の写像である。
$$
H_q=\bigoplus_{f(q_c,r_d)=q}\mathbb C^{k_c}\otimes\mathbb C^{l_d},\qquad
\Pi_q=\sum_{f(q_c,r_d)=q}P_c\otimes Q_d,\qquad n_q=\sum_{f(q_c,r_d)=q}k_cl_d .
$$
とおく。**融合積** $\mathcal M_1\otimes_f\mathcal M_2=\{\Pi_q\}'$ は、
$H_1\otimes H_2$ 上の作用素のうち全電荷を保存するものからなる。

**命題 7.5.** (a) $\mathcal M_1\otimes_f\mathcal M_2=\bigoplus_qB(H_q)\cong\bigoplus_qM_{n_q}(\mathbb C)$ であり、
極小中心射影は $\Pi_q$ である。(b) $\mathcal M_1\otimes\mathcal M_2\subseteq\mathcal M_1\otimes_f\mathcal M_2$ であり、
等号はすべての $H_q$ がただ一つの組 $(c,d)$ しか含まないときかつそのときに限り成り立つ。(c) 電荷が実数で
$f$ が加法ならば、$N=N_1\otimes1+1\otimes N_2$, $N_1=\sum_cq_cP_c$, $N_2=\sum_dr_dQ_d$ として
$\mathcal M_1\otimes_f\mathcal M_2=\{N\}'$ である。

*証明.* (a) 各組はちょうど一つの $H_q$ に属するから、$\Pi_q$ は互いに直交する射影で和は $1$ であり、
$X$ がそれらすべてと可換であるための必要十分条件は、$B(H_q)$ のブロックを任意として $X=\sum_q\Pi_qX\Pi_q$ となることである。中心は定理 7.1(b) と同様である。(b)
$\mathcal M_1\otimes\mathcal M_2$ の任意の元は各 $P_c\otimes Q_d$ と可換であり、したがってその和 $\Pi_q$ とも可換である。
各 $H_q$ が一つの $\mathbb C^{k_c}\otimes\mathbb C^{l_d}$ ならば、両方の環は
$\bigoplus_{c,d}B(\mathbb C^{k_c}\otimes\mathbb C^{l_d})$ である。$H_q$ が二つの組を含むならば、$B(H_q)$ は
一方の組のブロックを他方へ写す作用素を含み、それらはブロック対角でない。(c)
$N=\sum_{c,d}(q_c+r_d)P_c\otimes Q_d=\sum_qq\,\Pi_q$ であるから、零でない $\Pi_q$ は $N$ のスペクトル射影である。
それらは $N$ の多項式であり（定理 1.5(a) の証明）、$N$ はそれらの線形結合であるから、
$X$ が $N$ と可換であるための必要十分条件は、すべての $\Pi_q$ と可換であることである。$\square$

`fused_tensor_product(alg1, alg2, fuse)`（既定の `fuse` は加法）は、セクターを融合後の電荷で並べ、
セクター内の組 $(c,d)$ を辞書式に並べ、それぞれに基底 $e_a\otimes e_b$ を用いる。
例：ハードコア・ボソン。$L$ サイト上の $(\mathbb C^2)^{\otimes L}$ に占有数
$N=\sum_in^{(i)}$, $n=\operatorname{diag}(0,1)$ を考えると、固有空間 $N=q$ は $q$ 個のサイトが占有された $\binom Lq$ 個の積
基底ベクトルで張られるので、粒子数を保存する作用素 $\{N\}'$ は
$\bigoplus_qM_{\binom Lq}$ をなす（命題 4.10(2)）：2サイトでは $M_1\oplus M_2\oplus M_1$（電荷 $0,1,2$）、
3サイトでは $M_1\oplus M_3\oplus M_3\oplus M_1$ である。命題 7.5(c) により、融合積は基底の順序を除いて5サイトの全占有数に対する
$\{N\}'$ であり、実際
$n_q=\sum_c\binom2c\binom3{q-c}=\binom5q$ である（ヴァンデルモンドの恒等式：5サイトの $q$ 元部分集合は、最初の2サイトのうちの
$c$ 個と最後の3サイトのうちの $q-c$ 個からなる）。

```python
bos1 = TypeIAlgebra([1, 2, 1], [1, 2, 1], precision='double', charges=[0, 1, 2])
bos2 = TypeIAlgebra([1, 3, 3, 1], [1, 3, 3, 1], precision='double', charges=[0, 1, 2, 3])
F = fused_tensor_product(bos1, bos2)                         # U(1): total particle number
P = tensor_product(bos1, bos2)                               # particle numbers conserved separately
print(F.charges, F.k_factors, sum(k * k for k in F.k_factors), sum(k * k for k in P.k_factors))
assert F.k_factors == [math.comb(5, q) for q in range(6)]
assert sum(k * k for k in F.k_factors) == 252 and sum(k * k for k in P.k_factors) == 120

N1, N2 = full(bos1.central([0.0, 1.0, 2.0]))[0], full(bos2.central([0.0, 1.0, 2.0, 3.0]))[0]
N = torch.kron(N1, torch.eye(8, dtype=dt)) + torch.kron(torch.eye(4, dtype=dt), N2)
X = full12(F.operator(torch.randn(1, F.C, F.k_max, F.k_max, dtype=dt)))[0]
assert close(X @ N, N @ X)                                   # Proposition 7.5(c)
A, B = random_operator(bos1), random_operator(bos2)
assert close(full12(kron(A, B, F)), full12(kron(A, B, P)))   # Proposition 7.5(b)
```

融合環の次元は $\sum_q\binom5q^2=\binom{10}5=252$（再びヴァンデルモンドの恒等式）、通常の積の次元は
$(1+4+1)(1+9+9+1)=120$ である（定理 7.1(b)）。定理 7.1(a) と命題 7.5(b) により、通常の積は融合積の元のうち組 $(c,d)$
についてブロック対角なものからなる。残りの $132$ 次元は、組のブロック $(c,d)$ を同じ $H_q$ の別の組のブロック
$(c',d')$, $c'\ne c$（したがって $d'\ne d$）へ写す作用素で張られる：これらは全電荷を保存しつつ各部分の電荷を変える、
すなわちここでは二つの部分の間で粒子を移動させる。

$\mathcal M_1\otimes_f\mathcal M_2$ 上の状態は組のブロックの間にコヒーレンスをもちうるので、命題
7.2 は適用できない。$B(H_1)$ 上で $E_1(X)=\sum_cP_cXP_c$ とおく（2.6節の条件付き期待値）。これは $\mathcal M_1$ の上への写像であり、
$A\in\mathcal M_1$ に対し $\operatorname{Tr}(E_1(X)A)=\operatorname{Tr}(XA)$ である。$A=\sum_cP_cAP_c$ だからである。

**定理 7.6.** $\rho$ を $\mathcal M_1\otimes_f\mathcal M_2$ 上の状態 $\omega$ の密度とする。
その $\mathcal M_1$ への制限 $A\mapsto\omega(A\otimes1)$ の密度は
$$
\rho_1=E_1\big(\operatorname{Tr}_{H_2}\rho\big).
$$
である。$\operatorname{Tr}_{H_2}\rho$ のブロック $(c,c')$ は $\sum_d\operatorname{tr}_{l_d}\rho_{(c,d),(c',d)}$ であり、和は
$f(q_c,r_d)=f(q_{c'},r_d)$ を満たす $d$ にわたる。ここで $\rho_{(c,d),(c',d)}$ は $\rho$ の組 $(c',d)$ から組 $(c,d)$ へのブロックである。
したがって $E_1$ がすべての $\operatorname{Tr}_{H_2}\rho$ に自明に作用するための必要十分条件は、
$f(q_c,r_d)=f(q_{c'},r_d)$ となる $c\neq c'$ と $d$ が存在しないことである。たとえば $f$ が加法で
$q_c$ が互いに異なる場合がそうである。

*証明.* $A\otimes1\in\mathcal M_1\otimes_f\mathcal M_2$（命題 7.5(b)）であり、
$\omega(A\otimes1)=\operatorname{Tr}(\rho(A\otimes1))=\operatorname{Tr}\big((\operatorname{Tr}_{H_2}\rho)A\big)
=\operatorname{Tr}\big(E_1(\operatorname{Tr}_{H_2}\rho)A\big)$ である。$E_1(\operatorname{Tr}_{H_2}\rho)\in\mathcal M_1$ であるから、
定理 3.2 の一意性によりこれが密度である。ブロックの公式は $\operatorname{Tr}_{H_2}$ に対する命題 7.2(a) で、
$i$ をブロック $c$ に、$i'$ をブロック $c'$ に、$j$ をブロック $d$ にとったものである。
$\rho$ の $(c',d)$ と $(c,d)$ の間のブロックは、両方の組が同じ $H_q$ に属さない限り零である。
$f(q_c,r_d)=f(q_{c'},r_d)=q$ となる $c\neq c'$ と $d$ が存在するならば、単位ベクトル $x\in\mathbb C^{k_c}$,
$x'\in\mathbb C^{k_{c'}}$, $y\in\mathbb C^{l_d}$ をとり、$\psi=(x+x')\otimes y/\sqrt2\in H_q$ とおく。すると
$|\psi\rangle\langle\psi|\in B(H_q)$ であり、$\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|=\frac12|x+x'\rangle\langle x+x'|$
は零でないブロック $\frac12|x\rangle\langle x'|$ をもち、$E_1$ はそれを取り除く。$\square$

定理 7.6 のブロック公式により $(\rho_1)_c=\sum_d\operatorname{tr}_{l_d}\rho_{(c,d),(c,d)}$ である。ライブラリの
二つの実装はいずれも、任意の融合則に対してこの写像を計算するように作られている：`partial_trace` は
対角の組のブロック $\rho_{(c,d),(c,d)}$ だけのトレースをとり、`partial_trace_channel` は一つのセクター内で
$c$ が同じ組に別々のクラウス添字を用い、その双対がやはり $A\mapsto A\otimes1$ となるようにしている。コードは
両者を、下の二つの融合則と残す因子の両方の選び方について、`einsum` で独立に計算した
$E_1(\operatorname{Tr}_{H_2}\rho)$ と照合し、双対を $A\otimes1$ と照合する。
我々の例では、加法の場合 $E_1$ は何もしない。パリティ則 $f(q,r)=q+r\bmod2$ の場合には、`bos1` の電荷 $0$
と $2$ はどの $r$ とも同じように融合し、$E_1$ は本物のコヒーレンスを取り除く。

```python
G = fused_tensor_product(bos1, bos2, fuse=lambda q, r: (q + r) % 2)       # only the parity is conserved
print(G.k_factors)

def block_mask(alg):
    """Ones in the diagonal blocks of (+)_c C^{k_c}: X * block_mask(alg) = sum_c P_c X P_c."""
    return full(alg.operator(torch.ones(1, alg.C, alg.k_max, alg.k_max, dtype=dt)))

for alg, coherent in ((F, False), (G, True)):
    rho = states.random_density_matrix(alg, batch_size=2)
    R = full12(rho).reshape(2, 4, 8, 4, 8)
    for keep, X_H, alg_k in ((1, torch.einsum('xijkj->xik', R), bos1), (2, torch.einsum('xijil->xjl', R), bos2)):
        E_X = X_H * block_mask(alg_k)                                     # E(Tr_H rho), Theorem 7.6
        assert close(full(partial_trace(rho, keep=keep)), E_X)
        assert close(full(composite.partial_trace_channel(alg, keep=keep)(rho)), E_X)
        lost = (X_H - E_X).abs().max().item()                             # coherences removed by E
        print(f"keep={keep}: largest coherence removed: {lost:.3e}")
        assert (lost > 1e-3) if coherent else (lost == 0)
    assert close(composite.partial_trace_channel(alg, keep=1).adjoint('blunt')(A), kron(A, bos2.identity(), alg))
```

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. 超伝導体の BCS 記述では電子数のパリティだけが保存されるので、全体の状態はある
> 領域の粒子数 $0$ と $2$ を重ね合わせうる。オブザーバブルがその領域の粒子数を保存する観測者には、
> このコヒーレンスは見えない：その観測者にとって縮約状態は定理 7.6 のブロック対角化（ピンチング）されたもの（この囲みの数学的内容）、すなわち粒子数の混合である。

## 7.4 純粋状態のエンタングルメント

**補題 7.7（シュミット分解）.** 任意の単位ベクトル $\psi\in\mathbb C^a\otimes\mathbb C^b$ は、
$s_i>0$, $\sum_is_i^2=1$ と正規直交系
$u_i\in\mathbb C^a$, $v_i\in\mathbb C^b$ により $\psi=\sum_{i=1}^rs_i\,u_i\otimes v_i$ と書ける。このとき $\operatorname{Tr}_2|\psi\rangle\langle\psi|=\sum_is_i^2|u_i\rangle\langle u_i|$
かつ $\operatorname{Tr}_1|\psi\rangle\langle\psi|=\sum_is_i^2|v_i\rangle\langle v_i|$ である：二つの縮約状態は
同じ非零固有値 $s_i^2$ をもつ。

*証明.* $\Psi$ を $a\times b$ 行列 $\Psi_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,\psi\rangle$ とする。
スペクトル定理により、正規直交基底 $u_i$ と $\lambda_i=\|\Psi^*u_i\|^2\ge0$ を用いて
$\Psi\Psi^*=\sum_{i=1}^a\lambda_i|u_i\rangle\langle u_i|$ と書ける。ここで $\lambda_i>0$ となるのがちょうど $i\le r$ であるように並べておく。
$i\le r$ に対し $s_i=\sqrt{\lambda_i}$, $w_i=\Psi^*u_i/s_i$ とおくと、
$\langle w_i,w_j\rangle=\langle u_i,\Psi\Psi^*u_j\rangle/(s_is_j)=\delta_{ij}$ であり、$i>r$ に対しては $\Psi^*u_i=0$ である。
したがって $\Psi=\sum_i|u_i\rangle\langle u_i|\Psi=\sum_{i\le r}|u_i\rangle\langle\Psi^*u_i|=\sum_{i\le r}s_i|u_i\rangle\langle w_i|$、
すなわち $\Psi_{\alpha\beta}=\sum_is_i(u_i)_\alpha\overline{(w_i)_\beta}$ である：これは $v_i=\bar w_i$ としたときの主張であり、
$\sum_is_i^2=\operatorname{Tr}\Psi\Psi^*=\|\psi\|^2=1$ である。最後に
$|\psi\rangle\langle\psi|=\sum_{i,j}s_is_j|u_i\rangle\langle u_j|\otimes|v_i\rangle\langle v_j|$ であり、
命題 7.2(b) により縮約状態が得られる。$\square$

$\psi$ の **エンタングルメント・エントロピー** は
$E(\psi)=S(\operatorname{Tr}_2|\psi\rangle\langle\psi|)=S(\operatorname{Tr}_1|\psi\rangle\langle\psi|)=-\sum_is_i^2\log s_i^2$ である。
これが零となるのは $r=1$ のとき（各項 $-s_i^2\log s_i^2$ は $\ge0$ であり、$0$ となるのは $s_i^2=1$ のときに限る）、すなわち
$\psi=u\otimes v$ のとき（このとき $\Psi=uv^T$ の階数は $1$ であり、$r=\operatorname{rank}\Psi\Psi^*=\operatorname{rank}\Psi$）かつそのときに限る。また
$E(\psi)\le\log\min(a,b)$ であり、等号は小さい方の縮約状態が最大混合状態であるときかつそのときに限り成り立つ
（系 3.6 を両方の縮約状態に適用する）。

電荷がある場合、全電荷が確定した単位ベクトル $\psi$、すなわち融合積の一つのセクター $H_q$ に属するものをとる（異なる $H_q$ の
重ね合わせは、$\mathcal M_1\otimes_f\mathcal M_2$ 上で混合と同じ状態を与える。第3章、演習問題 3）。
$\psi_{cd}\in\mathbb C^{k_c}\otimes\mathbb C^{l_d}$ として $\psi=\sum_{(c,d)}\psi_{cd}$ と書く。

**命題 7.8（対称性分解エンタングルメント）.** $f$ を加法とし、$q_c$ は互いに異なり、$r_d$ も互いに異なるとする。
したがって各 $c$ に対し $q_c+r_d=q$ を満たす $d=d(c)$ は高々一つである。$\psi_c=\psi_{c\,d(c)}$（または $0$）、$p_c=\|\psi_c\|^2$ とおき、
$p_c>0$ に対し $\hat\psi_c=\psi_c/\sqrt{p_c}$ とおく。このとき $\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|$ は $\mathcal M_1$ に属し
（$N_1$ と可換である）、$\bigoplus_cp_c\operatorname{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|$ に等しい。
$\mathcal M_2$ 上の縮約状態は $d(c)$ においてセクター確率 $p_c$ をもつ。そして
$$
E(\psi)=H(p)+\sum_cp_c\,E(\hat\psi_c),\qquad H(p)=-\sum_cp_c\log p_c .
$$
である。

*証明.* 定理 7.6 のブロックの公式により、$\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|$ のブロック $(c,c')$ は
$q_c+r_d=q=q_{c'}+r_d$ を満たす $d$ にわたる $\operatorname{tr}_{l_d}|\psi_{cd}\rangle\langle\psi_{c'd}|$ の和であり、この条件から
$c=c'$ かつ $d=d(c)$ が従う。よって縮約状態は
$\bigoplus_c\operatorname{Tr}_2|\psi_c\rangle\langle\psi_c|$ であり、ブロックのトレースは $p_c$ である。因子の役割を入れ替えれば
$\mathcal M_2$ についての主張を得る。補題 7.7 により $E(\psi)$ はこの状態のエントロピーであり、
$\sigma_c=\operatorname{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|$, $S(\sigma_c)=E(\hat\psi_c)$ として定理 3.3 を適用すれば公式を得る。$\square$

$H(p)$ を **粒子数エンタングルメント**、$\sum_cp_cE(\hat\psi_c)$ を **配置エンタングルメント** という。
5サイト上の $2$ 粒子のランダムな状態について、両側の部分トレースと $4\times8$ 行列 $\Psi$ の特異値から
$E(\psi)$ を計算する。

```python
s = F.charges.index(2)                                      # the sector of two particles, n_q = 10
n = F.k_factors[s]
psi = torch.randn(n, dtype=dt)
psi = psi / psi.norm()
mat = torch.zeros(1, F.C, F.k_max, F.k_max, dtype=dt)
mat[0, s, :n, :n] = torch.outer(psi, psi.conj())
rho = DensityMatrix(F, matrix=mat)
r1, r2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)

psi_H = torch.zeros(4 * 8, dtype=dt)
psi_H[kron_rows(F)[s]] = psi                                # psi as a vector of H1 (x) H2
lam = torch.linalg.svdvals(psi_H.reshape(4, 8)) ** 2        # squared Schmidt coefficients
lam = lam[lam > 1e-14]
E_psi = -(lam * lam.log()).sum().item()
assert abs(r1.von_neumann_entropy().item() - E_psi) < 1e-10
assert abs(r2.von_neumann_entropy().item() - E_psi) < 1e-10  # Lemma 7.7

def entropy(sigma):
    w = torch.linalg.eigvalsh(sigma).clamp(min=1e-300)
    return -(w * w.log()).sum().item()

p = r1.sector_probabilities()[0]                            # c = 0, 1, 2 particles in the first part
assert torch.allclose(r2.sector_probabilities()[0, [2, 1, 0]], p)   # d(c) = 2 - c
H_num = -(p * p.log()).sum().item()
S_conf = sum(p[c].item() * entropy(r1.matrix[0, c, :k, :k] / p[c]) for c, k in enumerate(bos1.k_factors))
print(f"E = {E_psi:.4f} = H(p) {H_num:.4f} + configurational {S_conf:.4f}")
assert abs(E_psi - H_num - S_conf) < 1e-10                  # Proposition 7.8
```

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. 粒子数の超選択則に縛られた当事者たちは、配置の部分だけを資源として使うことができる：これが
> H. M. Wiseman と J. A. Vaccaro の *粒子のエンタングルメント* である（Phys. Rev. Lett. 91, 097902, 2003）。$H(p)$ は冷却原子において
> *粒子数エントロピー* として測定されている（A. Lukin et al., Science 364, 256, 2019）。

## 7.5 ランダムな状態とペイジの公式

$g\in\mathbb C^N$ の座標は独立で、その実部と虚部は独立に $N(0,\tfrac12)$ に従うとする（標準複素ガウスベクトル、
命題 A.16(a)）。したがって $\mathbb E|g_i|^2=1$ であり、$g\neq0$ はほとんど確実に成り立つ（定理 A.20）。
**ハール・ランダムな純粋状態** とは $\psi=g/\|g\|$ のことである。その分布はユニタリ不変である：すべてのユニタリ $U$ に対し
$Ug$ は $g$ と同じ分布に従い（命題 A.16(a)）、$Ug/\|Ug\|=U\psi$ である。命題 A.22 により、これは単位球面上の唯一の
ユニタリ不変な分布であり、その証明で示したように、ハール・ユニタリ（定理 5.3）の第1列の分布である。

**補題 7.9.** $\|g\|$ と $\psi=g/\|g\|$ は独立である。

*証明.*（これは定理 A.20(b) の一部である。）$h\ge0$ と $\phi\ge0$ をそれぞれ $(0,\infty)$ と
$\mathbb C^N\cong\mathbb R^{2N}$ の単位球面 $S^{2N-1}$ 上のボレル関数とする。$g$ は密度
$\pi^{-N}e^{-|x|^2}$ をもつ（命題 A.16(a)）から、極座標（$n=2N$ とした定理 A.18 と、そこで定義された測度 $\mu$）により
$$
\mathbb E\big[h(\|g\|)\phi(\psi)\big]=\pi^{-N}\!\int_{\mathbb R^{2N}}\!h(|x|)e^{-|x|^2}\phi\Big(\frac x{|x|}\Big)dx
=2N\pi^{-N}\int_0^\infty h(t)e^{-t^2}t^{2N-1}dt\int_{S^{2N-1}}\phi\,d\mu=:\alpha(h)\beta(\phi).
$$
を得る。$\phi=1$、次に $h=1$、次に両方をとると、$\mathbb Eh(\|g\|)=\alpha(h)\beta(1)$,
$\mathbb E\phi(\psi)=\alpha(1)\beta(\phi)$ および $\alpha(1)\beta(1)=1$ を得るから、
$\mathbb E[h(\|g\|)\phi(\psi)]=\mathbb Eh(\|g\|)\,\mathbb E\phi(\psi)$ である。$\square$

**補題 7.10.** $\mathbb E\|g\|^4=N(N+1)$ である。$G$ が $g_i$ と同じ分布に従う独立な成分をもつ $a\times b$ 行列
（*ジニブル* 行列）ならば、$\mathbb E\operatorname{Tr}(GG^*)^2=ab(a+b)$ である。

*証明.* $g_i=x+iy$ で $x,y$ が独立に $N(0,\frac12)$ に従うとすると、$\mathbb Ex^4=3(\mathbb Ex^2)^2=\frac34$
（定理 A.24：三つのペアリング）であるから、
$\mathbb E|g_i|^4=\mathbb Ex^4+2\mathbb Ex^2\mathbb Ey^2+\mathbb Ey^4=\tfrac34+\tfrac12+\tfrac34=2$ であり、
$\mathbb E\|g\|^4=\sum_i\mathbb E|g_i|^4+\sum_{i\ne j}\mathbb E|g_i|^2\mathbb E|g_j|^2=2N+N(N-1)$ である。次に、
$\operatorname{Tr}(GG^*)^2=\sum_{i,i',j,j'}G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}$ である。成分とその複素共役は、一つの中心化
ガウスベクトル（すべての成分の実部と虚部）の座標の複素線形結合であり、
$\mathbb E[G_{ij}G_{kl}]=\mathbb E[\bar G_{ij}\bar G_{kl}]=0$, $\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}$ を満たす
（命題 A.16(a)）。イッサーリス・ウィックの公式（定理 A.24 の複素形、四つの因子の三つのペアリング）により
$\mathbb E[G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}]=\delta_{ii'}\delta_{i'i}+0+\delta_{jj'}\delta_{jj'}=\delta_{ii'}+\delta_{jj'}$ である。
$i,i'\le a$ と $j,j'\le b$ について和をとると $ab^2+a^2b$ を得る。$\square$

**定理 7.11（E. Lubkin, J. Math. Phys. 19, 1028, 1978）.** ハール・ランダムな
$\psi\in\mathbb C^a\otimes\mathbb C^b$ と $\rho_A=\operatorname{Tr}_2|\psi\rangle\langle\psi|$ に対し、
$$
\mathbb E\operatorname{Tr}\rho_A^2=\frac{a+b}{ab+1}.
$$

*証明.* $G_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,g\rangle$ とおくと、これはジニブル行列である。$\psi$ の行列は
$\Psi=G/\|g\|$ であり、命題 7.2(a) により $(\rho_A)_{\alpha\alpha'}=\sum_\beta\Psi_{\alpha\beta}\bar\Psi_{\alpha'\beta}$、すなわち
$\rho_A=\Psi\Psi^*=GG^*/\|g\|^2$ である。したがって
$\operatorname{Tr}(GG^*)^2=\|g\|^4\operatorname{Tr}\rho_A^2$ であり、$\operatorname{Tr}\rho_A^2\in[0,1]$ は
$\psi$ の関数である。補題 7.9 により $\mathbb E\operatorname{Tr}(GG^*)^2=\mathbb E\|g\|^4\,\mathbb E\operatorname{Tr}\rho_A^2$ であり、補題 7.10 により
$\mathbb E\operatorname{Tr}\rho_A^2=ab(a+b)/(ab(ab+1))$ である。$\square$

第6章では、同じ平均をユニタリ不変な単位ベクトルの4次モーメント（補題 6.14）から、ランダムなチャネルの出力純度として
得る（定理 6.15 とその後の注意）。

**系 7.12.** $\log\frac{ab+1}{a+b}\le\mathbb E\,S(\rho_A)\le\log\min(a,b)$。

*証明.* 系 3.6 により $S(\rho_A)\le\log a$ であり、補題 7.7 により
$S(\rho_A)=S(\operatorname{Tr}_1|\psi\rangle\langle\psi|)\le\log b$ である。凹関数 $\log$ に対し、
$\rho_A$ の非零固有値 $\lambda_i>0$ を重みとして有限形のイェンセンの不等式（第3章の演習問題 2 の解答）を用いると
$\sum_i\lambda_i\log\lambda_i\le\log\sum_i\lambda_i^2$、すなわち $X=\operatorname{Tr}\rho_A^2\in(0,1]$ として $S(\rho_A)\ge-\log X$ を得る。
$m=\mathbb EX=\frac{a+b}{ab+1}$ とおく（定理 7.11）。$(0,\infty)$ 上の凸関数 $-\log$ に点 $m$ で補題 3.4 を適用すると
$-\log X\ge-\log m-(X-m)/m$ を得る。両辺は有界であり（コーシー・シュワルツにより $X\ge1/a$）、期待値をとると
$\mathbb E[-\log X]\ge-\log m$ を得る。$\square$

$a$ を固定して $b\to\infty$ とすると上下の評価は一致する：ランダムな純粋状態の小さな部分の平均エントロピーは、その最大値
$\log a$ に近づく。厳密な平均はペイジの公式であり、付録 C で証明される。

**定理 7.13（ペイジの公式）.** $a\le b$ に対し、
$$
\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}.
$$
さらに $\log a-\frac{a-1}{ab+1}-\frac{a-1}{2b}\le\mathbb E\,S(\rho_A)\le\log a-\frac{a-1}{2b}$ であり、したがって
$\big|\mathbb E\,S(\rho_A)-\log a+\frac a{2b}\big|\le\frac1{2b}$ である。2量子ビットでは $\mathbb E\,S(\rho_A)=\frac13$ である。

*証明.* これは定理 C.17 と系 C.18 である。付録 C の証明は
$\rho_A=GG^*/\operatorname{Tr}GG^*$（定理 7.11 の証明）から出発し、ウィシャート行列 $GG^*$ の固有値の同時密度（定理 C.13）と
$\rho_A$ の固有値の同時密度（定理 C.14）を計算し、$x\log x$ をラゲール準位密度に対して積分する
（補題 C.15 と C.16）。$\square$

この公式は D. N. Page が予想し（Phys. Rev. Lett. 71, 1291, 1993）、S. K. Foong と
S. Kanno（Phys. Rev. Lett. 72, 1148, 1994）、J. Sánchez-Ruiz（Phys. Rev. E 52, 5653, 1995）および
S. Sen（Phys. Rev. Lett. 77, 1, 1996）が証明した。

$a>b$ のときは $a$ と $b$ の役割が入れ替わる。実際、$\mathbb C^a\otimes\mathbb C^b$ から $\mathbb C^b\otimes\mathbb C^a$ の上への
ユニタリ $e_\alpha\otimes e_\beta\mapsto e_\beta\otimes e_\alpha$ は $g$ の座標を置換するので、ハール・ランダムな純粋状態 $\psi$ を
$\mathbb C^b\otimes\mathbb C^a$ のハール・ランダムな純粋状態 $\psi'$ に写す。その行列は $\Psi^T$ であるから
$\operatorname{Tr}_2|\psi'\rangle\langle\psi'|=\Psi^T\bar\Psi=\operatorname{Tr}_1|\psi\rangle\langle\psi|$ であり（命題 7.2(a)）、
補題 7.7 によりそのエントロピーは $S(\rho_A)$ である。したがって $a'=\min(a,b)$, $b'=\max(a,b)$ とすると、定理 7.13 は、
平均エンタングルメント・エントロピーがその最大値 $\log a'$ に少なくとも $\frac{a'-1}{2b'}$、高々
$\frac{a'-1}{2b'}+\frac{a'-1}{a'b'+1}$ だけ足りないことを述べている。後者は1ナット未満である：$b'=1$ ならば $0$ であり、
そうでなければ $\frac{a'-1}{2b'}\le\frac12$ かつ $\frac{a'-1}{a'b'+1}<\frac1{b'}\le\frac12$ である。ゆえに $ab$ を固定すると、
$a$ の関数として、平均は $\log\min(a,b)=\min(\log a,\log(ab)-\log a)$ との差が1ナット以内である。

$\rho_A$ の分布は **誘導測度** である（この名称は K. Życzkowski and H.-J. Sommers,
J. Phys. A 34, 7111, 2001 による。その固有値の密度は定理 C.14 が与える）。
`states.random_density_matrix(alg, rank=r)` は独立な $k_c\times r$ のジニブル・ブロック $G_c$ から
$G_cG_c^*$ を作り、$\sum_c\operatorname{Tr}G_cG_c^*$ で割る。ブロックが一つならば、これは上の $\rho_A$ の分布で
$b=r$ としたものである（定理 7.11 の証明）。ブロックが複数ならば、これは $H_1\otimes\mathbb C^r$ のハール・ランダムな
純粋状態 $\psi=g/\|g\|$ のベクトル状態の $\mathcal M_1$ への制限の密度の分布である。実際、
$G_{i\beta}=\langle e_i\otimes e_\beta,g\rangle$ を $D\times r$ ジニブル行列（$D=\sum_ck_c$）とし、$G_c$ をブロック $c$ に属するその行とする：
これらは独立な $k_c\times r$ ジニブル行列である。$\mathcal M_2=M_r$（セクターは一つ）と定数の融合則をとると、
全電荷はただ一つで $\Pi=1$ であり、$\mathcal M_1\otimes_f\mathcal M_2=B(H_1\otimes\mathbb C^r)$ である。定理 7.6 により、
制限は密度 $E_1(\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|)$ をもち、$\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|=GG^*/\|g\|^2$
（命題 7.2(a)、定理 7.11 の証明と同様）である。その対角ブロックは $G_cG_c^*/\|g\|^2$ であり、
$\|g\|^2=\sum_c\operatorname{Tr}G_cG_c^*$ である。

**命題 7.14.** $\bigoplus_cM_{k_c}$ 上で `random_density_matrix(alg, rank=r)` により抽出された $\rho$ と
$D=\sum_ck_c$ に対し、
$$
\mathbb E\operatorname{Tr}\rho^2=\frac{\sum_ck_c(k_c+r)}{D(Dr+1)} .
$$

*証明.* $g\in\mathbb C^{Dr}$ をすべての成分を集めたものとして $\operatorname{Tr}\rho^2=\sum_c\operatorname{Tr}(G_cG_c^*)^2/\|g\|^4$ であり、
これは $g/\|g\|$ の関数で $[0,1]$ に値をとる。定理 7.11 と同様に、補題 7.9 により
$\sum_c\mathbb E\operatorname{Tr}(G_cG_c^*)^2=\mathbb E\|g\|^4\,\mathbb E\operatorname{Tr}\rho^2$ であり、補題 7.10 により
$\mathbb E\operatorname{Tr}(G_cG_c^*)^2=k_cr(k_c+r)$ および $\mathbb E\|g\|^4=Dr(Dr+1)$ であるから、
$\mathbb E\operatorname{Tr}\rho^2=\sum_ck_cr(k_c+r)/(Dr(Dr+1))$ である。$\square$

`partial_trace` で縮約した $\mathbb C^3\otimes\mathbb C^6$ の $10^4$ 個の純粋状態と誘導測度について定理 7.11 と 7.13 を、
$M_1\oplus M_2\oplus M_3$, $r=2$ について命題 7.14 を確かめる。

```python
a, b, n_samples = 3, 6, 10_000
g = torch.randn(n_samples, a * b, dtype=dt)                 # complex Gaussian: E|g_i|^2 = 1
psi = g / g.norm(dim=-1, keepdim=True)
AB_alg = tensor_product(TypeIAlgebra([a], [a], precision='double'), TypeIAlgebra([b], [b], precision='double'))
pure = DensityMatrix(AB_alg, matrix=(psi[:, :, None] * psi.conj()[:, None, :])[:, None],
                     validate=False)                        # pure by construction
rho_A = partial_trace(pure, keep=1)
purity_exact = (a + b) / (a * b + 1)
page = sum(1 / k for k in range(b + 1, a * b + 1)) - (a - 1) / (2 * b)

def mean_se(x):
    return x.mean().item(), x.std().item() / math.sqrt(x.numel())

for name, rho_s in (('reduced pure states', rho_A),
                    ('induced measure', states.random_density_matrix(rho_A.algebra, 40_000, rank=b))):
    (m_p, e_p), (m_S, e_S) = mean_se(rho_s.purity()), mean_se(rho_s.von_neumann_entropy())
    print(f"{name}: E Tr rho^2 = {m_p:.4f} +- {e_p:.4f} (exact {purity_exact:.4f}), "
          f"E S = {m_S:.4f} +- {e_S:.4f} (Page {page:.4f})")
    assert abs(m_p - purity_exact) < 5 * e_p and abs(m_S - page) < 5 * e_S     # Theorems 7.11, 7.13
    assert math.log((a * b + 1) / (a + b)) < m_S < math.log(a)                # Corollary 7.12
assert math.log(a) - (a - 1) / (a * b + 1) - (a - 1) / (2 * b) <= page <= math.log(a) - (a - 1) / (2 * b)
assert abs(sum(1 / k for k in (3, 4)) - 1 / 4 - 1 / 3) < 1e-15                # two qubits: E S = 1/3

alg3 = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')
m_p, e_p = mean_se(states.random_density_matrix(alg3, 40_000, rank=2).purity())
print(f"E Tr rho^2 = {m_p:.4f} +- {e_p:.4f} (exact {26 / 78:.4f})")
assert abs(m_p - 26 / 78) < 5 * e_p                                          # Proposition 7.14
```

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. ランダムな純粋状態の平均エンタングルメント・エントロピーを、
> 全次元を固定して部分の大きさに対してプロットしたもの（定理 7.13 の後で述べた曲線）を *ペイジ曲線* という。物理学では、
> 情報がホーキング放射によってブラックホールからどのように出ていくか、そして閉じた系の部分がなぜ熱的に見えるかを
> 記述するために用いられる。これらは解釈であり、本講義の結果ではない。

## 7.6 比較の対象

単一のヒルベルト空間上では、QuTiP（`Qobj.ptrace`, `entropy_vn`, `rand_ket`, ジニブル分布と階数を指定する `rand_dm`）と
`numpy`（上のように reshape と `einsum`）が通常の積について部分トレース、エントロピー、誘導状態を計算する
（上の `einsum` による検証は、そのような計算をライブラリと比較している）。保存電荷をもつテンソルネットワークのライブラリ――ITensor（量子数によるブロック疎テンソル）、
TeNPy（電荷付きの脚、電荷ごとのエンタングルメント・スペクトル）、そして非可換対称性に対しては TensorKit.jl（融合木）――は、
大きな系に対して同じ仕方でセクターを融合し、対称性分解エンタングルメントを与える。ライブラリが付け加えるのは代数的な層である：
積はセクターとトレースをもつフォン・ノイマン環であり、融合則はラベルの任意の写像であり、部分トレースは
各トレースに関する双対をもつチャネルで、融合積に対しては正しい制限を与え（定理 7.3, 7.6）、そしてそのすべてがバッチ化されている。
テンソルネットワークは備えておらず、`fused_tensor_product` が受け取るのは、可換な対称性の場合のように電荷の組ごとに結果が一つの規則である。
$SU(2)$ の場合のような非可換な規則――スピンの組 $j_1,j_2$ が複数の全スピン $|j_1-j_2|,\dots,j_1+j_2$ に合成される
（第4章、演習問題 2）――は、ここでの融合積ではない（ライブラリはスピン $\frac12$ 鎖に対する $SU(2)$ を別に扱う、第4章）。

## 演習問題

1. $\mathcal M_1\otimes\mathcal M_2$ と $\mathcal M_1$ の両方の上の $\tau_{\mathrm{vN}}$ に関する
   $\operatorname{Tr}_2$ の双対を計算し、
   `partial_trace_channel(alg12).adjoint('tau_vN')` で確かめよ。
2. （劣加法性.）$\rho$ を $\mathcal M_1\otimes\mathcal M_2$ 上の状態で、縮約状態 $\rho_1,\rho_2$ が可逆であるものとする。
   $S(\rho)\le S(\rho_1)+S(\rho_2)$ であり、等号は
   $\rho=\rho_1\otimes\rho_2$ のときかつそのときに限り成り立つことを示せ。
3. `site` を電荷 $0,1$ をもつ $\mathbb C\oplus\mathbb C$（空のサイトまたは占有されたサイト）とし、
   加法によってそれ自身と融合する。$\psi=(|01\rangle+|10\rangle)/\sqrt2$（二つのサイトで共有された一つの粒子）について
   $E(\psi)$ とその二つの部分を計算し、
   $(|00\rangle+|11\rangle)/\sqrt2$ のベクトル状態の融合環への制限を記述せよ。

<details>
<summary>解答</summary>

1. 重みは $\mathcal M_1$ 上で $1/(C_1k_c)$、積の上で $1/(C_1C_2k_cl_d)$ である。
   $A\in\mathcal M_1$ に対し、命題 7.2(a) により
   $\tau_{\mathrm{vN}}((\operatorname{Tr}_2X)A)=\sum_c\frac1{C_1k_c}\sum_d\operatorname{Tr}(X_{cd}(A_c\otimes1))
   =\sum_{c,d}\frac1{C_1C_2k_cl_d}\operatorname{Tr}\big(X_{cd}\,C_2l_d(A_c\otimes1)\big)$ であるから、双対は
   $\mathcal M_2$ の中心元 $Z=\bigoplus_dC_2l_d1_{l_d}$ を用いて
   $A\mapsto\bigoplus_{c,d}C_2l_d\,A_c\otimes1_{l_d}=A\otimes Z$ である。これが単位的であるのは $Z=1$、すなわち
   $\mathcal M_2=\mathbb C$ のときに限り、この例では単位的でない（$Z=4\cdot1_2\oplus6\cdot1_3$）。`Tr_norm` に対しては同じ計算により
   $A\otimes\bigoplus_dl_d1_{l_d}$ を得て、これが単位的であるのはすべての $l_d=1$ のときに限る。blunt トレースは常に
   $\operatorname{Tr}_2$ を $A\mapsto A\otimes1$ の双対にする。一般には、同じ計算により、$\mathcal M_1$ 上の $\operatorname{tr}_w$ と
   積の上の重み $w_{cd}=w_c$ がそうする。

```python
dual_tau = composite.partial_trace_channel(alg12, keep=1).adjoint('tau_vN')
A = random_operator(alg1)
Z = alg2.central([alg2.C * l for l in alg2.k_factors])
assert close(dual_tau(A), kron(A, Z, alg12))
```

2. $\rho_1\otimes\rho_2$ は状態であり（定理 7.1(c)）、可逆で、固有ベクトル $u_i\otimes v_j$ と
   固有値 $\lambda_i\mu_j$ をもつから、$\log(\rho_1\otimes\rho_2)=\log\rho_1\otimes1+1\otimes\log\rho_2$ である
   （定理 1.5(a)）。部分トレースの定義により
   $\operatorname{Tr}\rho(\log\rho_1\otimes1)=\operatorname{Tr}\rho_1\log\rho_1$ であり、第二の因子についても同様である。
   したがって $S(\rho\|\rho_1\otimes\rho_2)=S(\rho_1)+S(\rho_2)-S(\rho)$ であり、これは $\ge0$ で、等号は
   $\rho=\rho_1\otimes\rho_2$ のときかつそのときに限り成り立つ（系 3.6）。この差が相互情報量である。

```python
rho = states.random_density_matrix(alg12, batch_size=4)
rho1, rho2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
mutual = rho1.von_neumann_entropy() + rho2.von_neumann_entropy() - rho.von_neumann_entropy()
assert (mutual > 0).all()
assert torch.allclose(mutual, states.relative_entropy(rho, kron(rho1, rho2, alg12)))
```

3. セクターは $N=0,1,2$ で、大きさは $1,2,1$ である。$\psi$ は $N=1$ に属し、$\psi_0=|01\rangle/\sqrt2$,
   $\psi_1=|10\rangle/\sqrt2$ であるから $p=(\frac12,\frac12)$ であり、ブロックが1次元であるから
   $E(\hat\psi_c)=0$ である：命題 7.8 により $E(\psi)=H(p)=\log2$ はすべて粒子数エンタングルメントであり、
   配置の部分は $0$ である。ベクトル $(|00\rangle+|11\rangle)/\sqrt2$ は $N=0$ と
   $N=2$ を重ね合わせる。ブロック対角な融合環の上では、そのベクトル状態は積状態の混合
   $\frac12|00\rangle\langle00|+\frac12|11\rangle\langle11|$ である（第3章、演習問題 3）。

```python
site = TypeIAlgebra([1, 1], [1, 1], precision='double', charges=[0, 1])
S2 = fused_tensor_product(site, site)
assert S2.k_factors == [1, 2, 1]
mat = torch.zeros(1, 3, 2, 2, dtype=dt)
mat[0, 1] = 0.5                                       # |psi><psi| with psi = (|01> + |10>)/sqrt 2
r1 = partial_trace(DensityMatrix(S2, matrix=mat), keep=1)
assert torch.allclose(r1.sector_probabilities()[0], torch.tensor([0.5, 0.5], dtype=torch.float64))
assert abs(r1.von_neumann_entropy().item() - math.log(2)) < 1e-12
```

</details>

**次へ：**[第8章](08_dynamics.md)ではダイナミクスを扱う：セクターをもつ環におけるシュレーディンガー方程式とリンドブラッド方程式、セクター間の跳躍、そして事後選択。
