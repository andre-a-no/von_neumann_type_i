[English](../13_jordan_algebras.md) · [Español](../es/13_jordan_algebras.md) · [Français](../fr/13_jordan_algebras.md) · [Deutsch](../de/13_jordan_algebras.md) · [中文](../zh/13_jordan_algebras.md) · **日本語** · [Русский](../ru/13_jordan_algebras.md)

# 13. ジョルダン代数：量子力学はなぜ *-環を用いるのか

二つのオブザーバブルの積 $xy$ は、$x$ と $y$ が可換でない限りオブザーバブルではないが、対称化された積
$x\circ y=\frac12(xy+yx)$ は常にオブザーバブルである。1933年にジョルダン（P. Jordan）は、この積を備えたオブザーバブルを量子力学の基本構造とすることを提案し、ジョルダン、フォン・ノイマン、ウィグナーはそこから生じる有限次元の代数を分類した。本章では有限次元の理論を、証明が短い場合には完全な証明とともに展開し、ジョルダン構造が何を見落としているかを問う。答えは交換子の符号である：一つの系のジョルダン代数は交換子を符号を除いてしか決定しない（13.5節）が、複合系はその符号を必要とし、しかもその符号が状態がエンタングルしているかどうかを決める。

**本章で学ぶこと**：ジョルダン積と、$\circ$ とは異なり正値性を保つ二次表現 $P(y)x=yxy$、スピン因子とそのスペクトル理論およびローレンツ錐、$*$-代数の自己共役部分ではないエルミート行列のジョルダン部分代数、Jordan–von Neumann–Wigner の分類、転置がなぜ $*$-自己同型ではなくジョルダン自己同型であるのか、そしてそれが複合系の上で何をするか、さらに第11章のトレース不等式をジョルダンの言葉で書き、スピン因子の中で証明すること。

本章を通じて、$\operatorname{Herm}_n(\mathbb K)$ は $\mathbb K=\mathbb R$ または $\mathbb C$ 上の $n\times n$ エルミート行列のなす実ベクトル空間、$M_{sa}$ は $\mathcal M=\bigoplus_cM_{k_c}(\mathbb C)$ の自己共役部分、
$X,Y,Z$ はパウリ行列で $\sigma=(X,Y,Z)$ とし、したがって $u\in\mathbb R^3$ に対し $u\cdot\sigma=u_1X+u_2Y+u_3Z$ である。
$[a,b]=ab-ba$ とし、$|u|$ は $u\in\mathbb R^m$ のユークリッドノルムである。

## 13.1 オブザーバブルとジョルダン積

**定義.** **ジョルダン代数** とは、双線形かつ可換な積 $x\circ y$ をもつ実ベクトル空間 $J$ であって、
**ジョルダン恒等式** $(x^2\circ y)\circ x=x^2\circ(y\circ x)$（ここで $x^2=x\circ x$）を満たすもののことである。ある結合代数の部分空間で $\frac12(xy+yx)$ について閉じたものと同型であるとき、
**特殊** であるという。$J$ 上で
$L(y)x=y\circ x$ とおき、**二次表現** を $P(y)=2L(y)^2-L(y^2)$ で定める。

**補題 13.1.** $\mathcal A$ を実または複素の結合代数とし、$x\circ y=\frac12(xy+yx)$ とする。このとき
$\circ$ はジョルダン恒等式を満たし、ジョルダン冪 $x^{\circ1}=x$、$x^{\circ(k+1)}=x\circ x^{\circ k}$ は結合的な冪 $x^k$ に等しい。特に $M_{sa}$、$\operatorname{Herm}_n(\mathbb C)$ および
$\operatorname{Herm}_n(\mathbb R)$ は特殊ジョルダン代数である。

*証明.* 展開すると $4(x^2\circ y)\circ x=x^2yx+yx^3+x^3y+xyx^2$ および
$4x^2\circ(y\circ x)=x^2yx+x^3y+yx^3+xyx^2$ であり、同じ四つの項からなる。帰納法により
$x\circ x^k=\frac12(x^{k+1}+x^{k+1})=x^{k+1}$ である。$x,y$ がエルミートならば $(xy+yx)^*=yx+xy$ であるから、三つの空間は $\circ$ について閉じている。実対称な $x,y$ に対しては積も実である。$\square$

**命題 13.2.** $x,y,z\in M_{sa}$ とする。

1. $P(y)x=yxy$ であり、**基本公式** $P(P(y)x)=P(y)P(x)P(y)$ が成り立つ。
2. $\operatorname{Tr}P(y)x=\operatorname{Tr}y^2x$。
3. $x\ge0$ ならば $P(y)x\ge0$ である。$y$ が可逆ならば、$P(y)$ は正錐からそれ自身への全単射であり、その
   逆写像は $P(y^{-1})$ である。
4. 正の元のジョルダン積は正とは限らない：$p=\begin{pmatrix}1&0\\0&0\end{pmatrix}$ と
   $q=\frac12\begin{pmatrix}1&1\\1&1\end{pmatrix}$ に対し $\det(p\circ q)=-\frac1{16}$ である。

*証明.* (1) $2y\circ(y\circ x)=\frac12(y^2x+2yxy+xy^2)$ および $y^2\circ x=\frac12(y^2x+xy^2)$ であるから、差をとればよい。次に結合律により
$P(P(y)x)z=yxy\,z\,yxy=P(y)P(x)P(y)z$ である。(2) トレースの巡回性。(3)
$yxy=(x^{1/2}y)^*(x^{1/2}y)\ge0$（系 1.7）であり、$P(y)P(y^{-1})=P(y^{-1})P(y)=\mathrm{id}$ である。(4)
$p\circ q=\frac14\begin{pmatrix}2&1\\1&0\end{pmatrix}$。$\square$

$P(y)$ を $2L(y)^2-L(y^2)$ と書けば (1) の恒等式はすべて $\circ$ だけで表されるので、すべての特殊ジョルダン代数で成り立つ。基本公式はすべてのジョルダン代数で成り立つ（K. McCrimmon, *A Taste of Jordan Algebras*, Springer 2004。ここではその証明を再現しない）。モジュール `torch_vn_algebra.jordan` はこれらの写像を `TypeIAlgebra` の作用素の上に、遅延評価で、かつすべてのセクターにわたってバッチ化して実装している：

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, SpinFactor
from torch_vn_algebra import jordan as jd

torch.manual_seed(0)
dt = torch.complex128
jp, P = jd.jordan_product, jd.quadratic_representation

def close(A, B, tol=1e-10):
    A = A.matrix if hasattr(A, 'matrix') else A
    B = B.matrix if hasattr(B, 'matrix') else B
    return (A - B).abs().max().item() < tol

alg = TypeIAlgebra([2, 3], [2, 3], precision='double')           # M_2 (+) M_3
def sa(B=4):
    return alg.operator_from_eigenvalues(lambda d: torch.randn(B, d), batch_size=B, force_self_adjoint=True)
x, y, z = sa(), sa(), sa()
x2 = jp(x, x)
assert close(jp(jp(x2, y), x), jp(x2, jp(y, x)))                  # Jordan identity (Lemma 13.1)
assert close(jd.jordan_power(x, 5), x @ x @ x @ x @ x)            # power associativity
assert close(P(y, x), 2 * jp(y, jp(y, x)) - jp(jp(y, y), x))      # P(y) = 2L(y)^2 - L(y^2)
assert close(P(P(y, x), z), P(y, P(x, P(y, z))))                  # fundamental formula
assert close(P(y, x).trace, (jp(y, y) @ x).trace)                 # Tr P(y)x = Tr y^2 x
pos = alg.operator_from_eigenvalues(lambda d: torch.rand(4, d), batch_size=4, force_positive=True)
assert all((e > -1e-12).all() for e in P(y, pos).eigenvalues())  # Proposition 13.2(3)

M2 = TypeIAlgebra([2], [2], precision='double')
p = M2.operator(torch.tensor([[[[1.0, 0], [0, 0]]]]), is_self_adjoint=True, is_positive=True)
q = M2.operator(torch.tensor([[[[0.5, 0.5], [0.5, 0.5]]]]), is_self_adjoint=True, is_positive=True)
print(jp(p, q).eigenvalues()[0])                                  # (1 -+ sqrt 2)/4: not positive
assert abs(torch.linalg.det(jp(p, q).matrix[0, 0]).real.item() + 1 / 16) < 1e-14
```

> **物理的意味.** $x\circ y$ はオブザーバブルの対称化された積である：それは $x,y$ について対称で $y=x$ のとき $x^2$ に帰着する唯一の組合せ
> $\alpha xy+\beta yx$ であり、再びオブザーバブルとなる。二次表現は **リューダース則** である：イエス／ノーの問い $p$ を測定し、結果「イエス」を残すと、
> 密度 $\rho$ は $P(p)\rho/\operatorname{Tr}P(p)\rho=p\rho p/\operatorname{Tr}(p\rho)$ に変わる。命題 13.2(4)
> は、$\circ$ そのものには正の元に対してこのような操作的意味がないことを述べている。

## 13.2 スピン因子とローレンツ錐

**定義.** $m\ge1$ に対し、**スピン因子** $V_m=\mathbb R\oplus\mathbb R^m$ は元 $x=(a,u)$、
積 $(a,u)\circ(b,v)=(ab+u\cdot v,\;av+bu)$ および単位元 $1=(1,0)$ をもつ。その **トレース** は $\operatorname{tr}x=2a$ である。

積は双線形かつ可換である。ジョルダン恒等式は行列表現から従う。

**補題 13.3（反可換な対称元）.** $p=\lfloor m/2\rfloor$、$D=2^p$（$m=1$ のときは $D=2$）とする。
$(\mathbb C^2)^{\otimes p}$ 上で、$j\le p$ に対し $s_{2j-1}=Z^{\otimes(j-1)}\otimes X\otimes1$、$s_{2j}=Z^{\otimes(j-1)}\otimes Y\otimes1$
とおき、$m$ が奇数ならば $s_m=Z^{\otimes p}$ とおく（$m=1$ のときは $s_1=Z$）。このとき $s_i=s_i^*$ かつ
$s_is_j+s_js_i=2\delta_{ij}1$ である。写像 $\phi(a,u)=a1+\sum_iu_is_i$ は単射であり、
$\phi(x\circ y)=\phi(x)\circ\phi(y)$ を満たす。したがって $V_m$ は特殊ジョルダン代数である。$m\ge2$ のとき、$s_i$ たちは $*$-代数として
$M_D(\mathbb C)$ を生成する。

*証明.* パウリ行列は二乗すると $1$ になり、互いに反可換である。$Z$ と $X$、$Z$ と $Y$ も同様である。二つの異なる文字列
$s_i,s_j$ は、ちょうど一つのサイトで反可換な仕方で異なる：$s_{2j-1},s_{2j}$ についてはサイト $j$ で（$X$ 対
$Y$）、添字 $j<j'$ をもつ文字列についてはサイト $j$ で（$X$ または $Y$ 対 $Z$）異なり、その他のサイトでは等しい因子か因子 $1$ をもつ。$s_m=Z^{\otimes p}$ についても同じことが成り立つ。よって $i\ne j$ に対し $s_is_j=-s_js_i$ であり、
$s_i^2=1$ である。ゆえに $(u\cdot s)\circ(v\cdot s)=\sum_{ij}u_iv_j\,s_i\circ s_j=(u\cdot v)1$ であり、これから
$\phi(x)\circ\phi(y)=\phi(x\circ y)$ が得られる。各 $s_i$ はトレースが零であり（トレース零の因子をもつテンソル積である）、
$\operatorname{Tr}(s_is_j)=\operatorname{Tr}(s_i\circ s_j)=D\delta_{ij}$ であるから、$1,s_1,\dots,s_m$ は
$\operatorname{Tr}(ab)$ について直交し、$\phi$ は単射である。ジョルダン恒等式は $\operatorname{Herm}_D$ で成り立ち（補題 13.1）、
$\phi$ によって引き戻される。最後に $s_{2j-1}s_{2j}=1\otimes\cdots\otimes XY\otimes\cdots\otimes1=iZ_j$（$Z$ はサイト
$j$ にある）であるから、生成される $*$-代数は $Z_j$ を含み、したがって $X_j=Z_1\cdots Z_{j-1}s_{2j-1}$ と $Y_j$ を含み、それゆえパウリ行列のすべての
テンソル積、すなわち $M_D(\mathbb C)$ の基底を含む。$\square$

符号を除けば、$s_{2j-1}$ と $s_{2j}$ は 9.4節のジョルダン–ウィグナー・フェルミオン $c_j$ の
マヨラナ演算子 $c_j+c_j^*$ と $i(c_j-c_j^*)$ である。

**定理 13.4（$V_m$ のスペクトル理論）.** $x=(a,u)\in V_m$、$\lambda_\pm=a\pm|u|$ とする。

1. $V_m$ の冪等元は $0$、$1$ および $|w|=1$ なる $c_w=\frac12(1,w)$ である。また $c_w\circ c_{w'}=0$ となるのは $w'=-w$ のときかつそのときに限る。
2. $u\ne0$ ならば $e_\pm=c_{\pm u/|u|}$ として $x=\lambda_+e_++\lambda_-e_-$ であり、$u=0$ ならば $x=a1$ である。$x=\mu_1c_w+\mu_2c_{-w}$ ならば
   $\{\mu_1,\mu_2\}=\{\lambda_+,\lambda_-\}$ である。$\lambda_\pm$ を $x$ の **固有値** とよぶ。
3. $f:\{\lambda_+,\lambda_-\}\to\mathbb R$ に対し $f(x)=f(\lambda_+)e_++f(\lambda_-)e_-$（$u=0$ のときは $f(a)1$）とおく。このとき
   $x^{\circ k}=\lambda_+^ke_++\lambda_-^ke_-$、$\phi(f(x))=f(\phi(x))$ であり、$\phi(x)$ の固有値は $\lambda_\pm$ で、
   それぞれの重複度は $D/2$、そして $\operatorname{Tr}\phi(x)=\frac D2\operatorname{tr}x$ である。また
   $\operatorname{tr}x=\lambda_++\lambda_-$、$\det x:=\lambda_+\lambda_-=a^2-|u|^2$ である。
4. 平方の錐 $Q=\{x^2:x\in V_m\}$ は **ローレンツ錐** $\{(a,u):a\ge|u|\}=\{x:\lambda_-\ge0\}=\{x:\phi(x)\ge0\}$ である。
   これは **自己双対** である：$y\in Q$ であるための必要十分条件は、すべての $x\in Q$ に対し $\operatorname{tr}(x\circ y)\ge0$ となることである。
5. $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$。したがってすべての $y\in V_m$ に対し $P(y)Q\subseteq Q$ である。

*証明.* (1) $(a,u)^2=(a^2+|u|^2,2au)$ が $(a,u)$ に等しいのは、$u=0$ かつ $a\in\{0,1\}$ のとき、または $a=\frac12$ かつ $|u|=\frac12$ のときに限る。
$c_w\circ c_{w'}=\frac14(1+w\cdot w',\,w+w')$ が零になるのは $w'=-w$ のときかつそのときに限る。(2)
$\lambda_+e_++\lambda_-e_-=(a,0)+|u|(0,u/|u|)=x$ であり、$\mu_1c_w+\mu_2c_{-w}=(\frac{\mu_1+\mu_2}2,\frac{\mu_1-\mu_2}2w)$ であるから、
$\mu_1+\mu_2=2a$ かつ $|\mu_1-\mu_2|=2|u|$ である。(3) (1) により $(\alpha e_++\beta e_-)\circ(\gamma e_++\delta e_-)=\alpha\gamma
e_++\beta\delta e_-$ であり、帰納法で冪の公式を得る。$w=u/|u|$ とすると $(w\cdot s)^2=1$ であるから、
$\phi(e_\pm)=\frac12(1\pm w\cdot s)$ は和が $1$、積が $0$ の直交射影で、それぞれのトレースは $D/2$ である。したがって
$\phi(x)=\lambda_+\phi(e_+)+\lambda_-\phi(e_-)$ はスペクトル分解であり、
$f(\phi(x))=f(\lambda_+)\phi(e_+)+f(\lambda_-)\phi(e_-)=\phi(f(x))$ である（定理 1.5）。(4) $x^2$ の固有値は
$\lambda_\pm^2\ge0$ である。逆に $\lambda_-\ge0$ ならば $x=(\sqrt x)^2$ である。$\lambda_-\ge0$ は $a\ge|u|$ と同値であり、(3) により
$\phi(x)\ge0$ とも同値である。自己双対性については、$\operatorname{tr}(x\circ y)=2(ab+u\cdot v)$ である。$a\ge|u|$ かつ $b\ge|v|$ ならば
$ab\ge|u||v|\ge-u\cdot v$ である。$b<|v|$ ならば、$x=(|v|,-v)\in Q$ が $2(|v|b-|v|^2)<0$ を与える（$v=0$、$b<0$ のときは $x=1$ をとる）。
(5) $P(y)x=2y\circ(y\circ x)-y^2\circ x$ であり $\phi$ は $\circ$ を保つので、命題 13.2(1) により $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$ である。
あとは (4) と命題 13.2(3) を適用すればよい。$\square$

**命題 13.5.** (a) $m=3$ のとき、$\phi(a,u)=a1+u\cdot\sigma$ はジョルダン同型
$V_3\cong\operatorname{Herm}_2(\mathbb C)$ である。(b) $V_m$ のジョルダン自己同型は写像 $(a,u)\mapsto(a,Ou)$、
$O\in O(m)$ である。(c) $m\ge2$ のとき、$\phi(V_m)$ の次元は $m+1$ であり、それが生成する $*$-代数の自己共役部分の次元は
$D^2=4^{\lfloor m/2\rfloor}$ である。両者が等しいのは $m=3$ のときに限る。

*証明.* (a) $m=3$ のとき $(s_1,s_2,s_3)=(X,Y,Z)$ であり、$\phi$ は次元 $4$ の空間の間の単射である。(b) 写像
$(a,u)\mapsto(a,Ou)$ は $u\cdot v$ を保ち、したがって $\circ$ を保つ。逆に $\alpha$ をジョルダン自己同型とする。すると
$\alpha(1)$ は $V_m$ の単位元であるから $\alpha(1)=1$ である。$x^2$ の公式により、$x^2=1$ となるのは $x=\pm1$ のとき、または
$|w|=1$ なる $x=(0,w)$ のときに限る。$|w|=1$ に対し $\alpha(0,w)^2=1$ かつ $\alpha(0,w)\ne\pm1=\alpha(\pm1)$ であるから、$|w'|=1$ なる $w'$ により $\alpha(0,w)=(0,w')$ となる。
線形な $\ell,O$ を用いて $\alpha(0,u)=(\ell(u),Ou)$ と書けば、これは $\ell=0$ かつ $|Ou|=|u|$、すなわち $O\in O(m)$ を意味する。
(c) 補題 13.3。$m=2$ および $m\ge4$ に対し $4^{\lfloor m/2\rfloor}>m+1$ である。$\square$

```python
V = SpinFactor(5)
xs, ys = V.random(1000), V.random(1000)                          # (1000, 6): batches of (a, u)
x2 = V.square(xs)
assert (V.product(V.product(x2, ys), xs) - V.product(x2, V.product(ys, xs))).abs().max() < 1e-12
lam, e = V.spectral_decomposition(xs)                            # (lambda-, lambda+), (e-, e+)
assert (lam[:, 1, None] * e[:, 1] + lam[:, 0, None] * e[:, 0] - xs).abs().max() < 1e-12
Xm, Ym, Dm = V.to_matrix(xs), V.to_matrix(ys), V.representation_dim  # phi into Herm_4(C)
assert (V.to_matrix(V.product(xs, ys)) - 0.5 * (Xm @ Ym + Ym @ Xm)).abs().max() < 1e-12
assert (torch.linalg.eigvalsh(Xm) - lam.repeat_interleave(Dm // 2, dim=-1)).abs().max() < 1e-12
assert (V.to_matrix(V.exp(xs)) - torch.linalg.matrix_exp(Xm)).abs().max() < 1e-9       # phi(f(x)) = f(phi(x))
cone = V.random(1000, cone=True)                                 # squares
assert V.in_cone(cone).all() and (torch.linalg.eigvalsh(V.to_matrix(cone)) > -1e-12).all()
assert (V.square(V.sqrt(cone)) - cone).abs().max() < 1e-10
O, _ = torch.linalg.qr(torch.randn(5, 5, dtype=torch.float64))  # Proposition 13.5(b)
rot = lambda w: torch.cat([w[..., :1], w[..., 1:] @ O.T], dim=-1)
assert (rot(V.product(xs, ys)) - V.product(rot(xs), rot(ys))).abs().max() < 1e-12

for m in range(1, 7):                                            # Proposition 13.5(c)
    r = jd.compare_with_star_algebra(jd.anticommuting_symmetries(m))
    print(f"m = {m}: dim V_m = {r['jordan_dim']}, dim of the self-adjoint part of the *-algebra = {r['star_dim']}")
    assert r['jordan_dim'] == m + 1 and r['equal'] == (m in (1, 3))
```

$m=1$ のとき、$s_1=Z$ は対角行列を生成し、その次元は $2=m+1$ であり、$V_1\cong\mathbb R\oplus\mathbb R$ は単純ではない。表示される次元は、$m=2$ で $3,4$、$m=4$ で $5,16$、$m=6$ で $7,64$ である。

> **物理的意味.** 自己双対性により、$V_m$ の状態（正で正規化された線形汎関数）は
> $\rho=\frac12(1,r)\in Q$、$|r|\le1$ を用いて $x\mapsto\operatorname{tr}(\rho\circ x)$ と表される：状態空間は次元
> $m$ の球体である。$m=3$ のとき、それは量子ビットのブロッホ球である（命題 13.5(a)）。ローレンツ錐は $1+m$ 次元の特殊相対論における未来の
> 光円錐であり、命題 13.5(b) は、ジョルダン自己同型が時間軸を固定する空間の直交変換であることを述べている。

錐 $Q$ は **二次錐計画**（SOCP）の制約集合である：すなわちローレンツ錐の直積上の線形最適化である。SOCP の内点法はまさに定理 13.4 を用いる：固有値、$\det$、そしてスケーリングのための
$P(y)$ である（F. Alizadeh, D. Goldfarb, *Second-order cone programming*, Math. Program. 95 (2003) 3–51）。

## 13.3 *-部分代数ではないジョルダン部分代数

ジョルダン部分代数 $J\subseteq\operatorname{Herm}_n(\mathbb C)$ は複素 $*$-代数 $\mathcal A$ を生成し、
$J\subseteq\mathcal A_{sa}$ である。命題 13.5(c) は $J\ne\mathcal A_{sa}$ となる例を与えた。ここではさらに二つの例を挙げる。

**命題 13.6.** (a) $\operatorname{Herm}_n(\mathbb R)$ は $\operatorname{Herm}_n(\mathbb C)$ の次元
$n(n+1)/2$ のジョルダン部分代数である。それが生成する $*$-代数は $M_n(\mathbb C)$ であり、その自己共役部分の次元は $n^2$ である。
(b) $x,y\in\operatorname{Herm}_2(\mathbb C)$ に対し、$UxU^*$ と $UyU^*$ がともに実対称となるユニタリ $U$ が存在する。$1,x,y$ が生成するジョルダン代数は $\operatorname{span}\{1,x,y\}$ である。$x,y,1$ が線形独立ならば、それは次元 $3$ の $U^*\operatorname{Herm}_2(\mathbb R)U\cong V_2$ であり、一方
$x,y$ が生成する $*$-代数は次元 $4$ の $M_2(\mathbb C)$ である。

*証明.* (a) 補題 13.1。$i\ne j$ に対し $E_{ij}=E_{ii}(E_{ij}+E_{ji})$ は実対称行列の積であるから、
生成される代数はすべての行列単位を含む。(b) $x=a1+u\cdot\sigma$、$y=b1+v\cdot\sigma$ と書き、単位ベクトル
$n\perp u,v$ を選ぶ。行列 $n\cdot\sigma$ と $Y$ はエルミートでトレース零、二乗すると $1$ であるから、ともに
固有値 $1,-1$ をもち、あるユニタリ $U$ に対し $U(n\cdot\sigma)U^*=Y$ となる（定理 1.3）。$U$ による共役は
$\operatorname{Tr}$ と $\circ$ を保つので、それは $u\cdot\sigma$ を $(Ru)\cdot\sigma$ に写す。ここで $R$ は線形で、
$Ru\cdot Rv=\frac12\operatorname{Tr}\big((Ru\cdot\sigma)(Rv\cdot\sigma)\big)=u\cdot v$ であるから $R\in O(3)$ であり、$Rn=e_2$ である。
よって $Ru,Rv\perp e_2$ であり、$X,Z$ が実であるから $UxU^*=a1+(Ru)_1X+(Ru)_3Z$ は実である。$UyU^*$ も同様である。次に
$x\circ y=(ab+u\cdot v)1+av\cdot\sigma+bu\cdot\sigma\in\operatorname{span}\{1,x,y\}$ であり、$x^2,y^2$ も同様であるから、この線形包は
閉じている。その次元が $3$ ならば、$U$ による共役は $\operatorname{Herm}_2(\mathbb R)=\operatorname{span}\{1,X,Z\}$ 全体であり、
これは補題 13.3 により $\cong V_2$ であり（$X,Z$ は反可換）、(a) により $M_2(\mathbb C)$ を生成する。$\square$

```python
for n in range(2, 6):                                            # Proposition 13.6(a), a random real pair
    G = torch.randn(2, n, n, dtype=torch.float64)
    r = jd.compare_with_star_algebra(G + G.mT)
    assert (r['jordan_dim'], r['star_dim']) == (n * (n + 1) // 2, n * n)

H = torch.randn(3, 2, 2, dtype=dt)
H = H + H.mH                                                     # three random Hermitian 2 x 2 matrices
r2, r3 = jd.compare_with_star_algebra(H[:2]), jd.compare_with_star_algebra(H)
print("two generators:", r2['jordan_dim'], r2['star_dim'], "  three:", r3['jordan_dim'], r3['star_dim'])
assert (r2['jordan_dim'], r2['star_dim'], r3['jordan_dim']) == (3, 4, 4)

V3 = SpinFactor(3)
u, v = V3.from_matrix(H[0])[1:], V3.from_matrix(H[1])[1:]        # Bloch vectors
nrm = torch.linalg.cross(u, v)
nrm = nrm / nrm.norm()
_, W1 = torch.linalg.eigh(V3.to_matrix(torch.cat([torch.zeros(1, dtype=torch.float64), nrm])))
_, W2 = torch.linalg.eigh(jd.anticommuting_symmetries(3)[1])     # Y
Um = W2 @ W1.mH                                                  # U (n . sigma) U* = Y
assert all((Um @ A @ Um.mH).imag.abs().max() < 1e-12 for A in H[:2])
```

一般の三つの $2\times2$ エルミート行列は $\operatorname{Herm}_2(\mathbb C)$ 全体を生成する：両方の列で次元は $4$ である。したがって $\mathcal M$ の $M_2$ ブロックの中では、二つのオブザーバブルが複素構造を明らかにすることは決してないが、三つ目のオブザーバブルはそれを明らかにしうる。

## 13.4 分類

ジョルダン代数は、$\sum_ix_i^2=0$ ならばすべての $x_i=0$ となるとき **形式的実** であるという。有限次元では、これは
**ユークリッド的** であること、すなわち $\langle x\circ y,z\rangle=\langle y,x\circ z\rangle$ を満たす内積をもつことと同値である
[FK94]。$M_{sa}$（$\operatorname{Tr}(xy)$ による）と $V_m$（$\operatorname{tr}(x\circ y)$ による）はユークリッド的である。

**定理 13.7（Jordan–von Neumann–Wigner。証明なしに述べる）.** 有限次元の形式的実ジョルダン代数はすべて単純なものの直和であり、単純なものはどれも次のうちちょうど一つと同型である：$\mathbb R$、
$\mathbb K=\mathbb R,\mathbb C,\mathbb H$ かつ $r\ge3$ に対する $\operatorname{Herm}_r(\mathbb K)$、スピン因子
$V_m$（$m\ge2$）、八元数上の $3\times3$ エルミート行列のなす **アルバート代数** $\operatorname{Herm}_3(\mathbb O)$（次元 $27$）。階数 $2$ では、
$\mathbb K=\mathbb R,\mathbb C,\mathbb H,\mathbb O$ に対し $\operatorname{Herm}_2(\mathbb K)\cong V_{\dim\mathbb K+1}$ である。（P. Jordan, J. von Neumann, E. Wigner, *On an algebraic
generalization of the quantum mechanical formalism*, Ann. of Math. 35 (1934) 29–64。現代的な証明は [FK94,
Chapter V] にある。）

アルバート代数を除くこれらすべては特殊である：$\operatorname{Herm}_r(\mathbb H)$ は四元数を $2\times2$ 複素行列として書くことで
$\operatorname{Herm}_{2r}(\mathbb C)$ に埋め込まれ、$V_m$ は補題 13.3 により埋め込まれる。
アルバート代数は **例外的** である：それはいかなる結合代数のジョルダン部分代数とも同型でない
（A. A. Albert, *On a certain algebra of quantum mechanics*, Ann. of Math. 35 (1934) 65–73）。短く述べられる理由は次のとおりである：グレニー（C. M. Glennie）は、すべての特殊ジョルダン代数で成り立つが
$\operatorname{Herm}_3(\mathbb O)$ では成り立たない多項式恒等式を見出した（*Some identities valid in special Jordan algebras but not valid in all Jordan
algebras*, Pacific J. Math. 16 (1966) 47–59。McCrimmon, *A Taste of Jordan Algebras* も参照）。八元数は結合的でなく、$\operatorname{Herm}_r(\mathbb O)$ がジョルダン恒等式を満たすのは $r\le3$ のときに限る。これがリストがそこで止まる理由である。したがって $\operatorname{Herm}_3(\mathbb O)$ は、積 $\frac12(xy+yx)$ をもつ実または複素行列によって（あるいはいかなる結合代数の中でも）忠実に表現することはできず、ソフトウェアは八元数の成分を直接扱わなければならない。

**定理 13.8（ケッヒャー–ヴィンベルク。証明なしに述べる）.** ユークリッド的ジョルダン代数の平方の錐の内部は、ちょうど **対称錐**、すなわちユークリッド空間内の開凸錐であって自己双対かつ
等質である（その線形自己同型が推移的に作用する）もの、に一致する。[FK94, Chapter III.]

定理 13.4(4) と演習問題 1 は、ローレンツ錐についてこの二つの性質を証明する。したがって半正定値計画の錐
（$\operatorname{Herm}_n(\mathbb R)^+$ または $\operatorname{Herm}_n(\mathbb C)^+$ の直積）、SOCP の錐、そして線形計画の錐
（$\mathbb R_+^n$、すなわち $\mathbb R^n=\mathbb R\oplus\dots\oplus\mathbb R$ の錐）はすべてこの形をしており、
内点法はすべてのユークリッド的ジョルダン代数に対して一度に書くことができる。

[FK94]: J. Faraut, A. Korányi, *Analysis on Symmetric Cones*, Oxford University Press, 1994.

## 13.5 ジョルダン写像と *-写像

フォン・ノイマン環の **ジョルダン $*$-準同型** とは、$\Phi(x^*)=\Phi(x)^*$ と
$\Phi(x\circ y)=\Phi(x)\circ\Phi(y)$ を満たす複素線形写像 $\Phi$ のことである。自己共役部分の上ではそれはジョルダン準同型であり、自己共役部分の間のすべてのジョルダン
準同型は $\Phi(x+iy)=\Phi(x)+i\Phi(y)$ によってジョルダン $*$-準同型に拡張される。

**命題 13.9（転置）.** $\mathcal M$ 上で $T(A)=\bigoplus_cA_c^T$ とする。(a) $T$ は単位的なジョルダン
$*$-自己同型であり、$T(AB)=T(B)T(A)$ を満たす。$T$ はスペクトルを保つので、$T$ と $T^{-1}=T$ は正値である。(b) $M_k$
（$k\ge2$）上で、チョイ行列 $\sum_{ij}E_{ij}\otimes T(E_{ij})$ はスワップ $F(e_a\otimes e_b)=e_b\otimes e_a$ であり、その
固有値は $+1$（重複度 $k(k+1)/2$）と $-1$（重複度 $k(k-1)/2$）である。したがって $T$ は完全正値ではない。

*証明.* (a) $(A^T)^*=(A^*)^T$、$(AB)^T=B^TA^T$ であるから $(A\circ B)^T=B^T\circ A^T=A^T\circ B^T$ である。また
$\det(A^T-\lambda)=\det(A-\lambda)$ である。(b) $\sum_{ij}E_{ij}e_a\otimes E_{ji}e_b=e_b\otimes e_a$。$F^2=1$ であり、$F$ は対称テンソル上で $+1$、
反対称テンソル上で $-1$ であって、これらは空間全体を張る。定理 6.2（および命題 6.1）により $T$ は
CP ではない。$\square$

したがってジョルダン自己同型は $*$-自己同型であるとは限らない。この現象を除けば、正値性がジョルダン構造を決定する：

**定理 13.10（カディソン。証明なしに述べる）.** $\Phi$ を単位的 C\*-環の間の線形全単射で、
$\Phi(1)=1$ を満たし、$\Phi$ と $\Phi^{-1}$ がともに正値であるものとする。このとき $\Phi$ はジョルダン $*$-同型である。（R. V. Kadison,
*A generalized Schwarz inequality and algebraic invariants for operator algebras*, Ann. of Math. 56 (1952) 494–503,
Corollary 5。*Isometries of operator algebras*, Ann. of Math. 54 (1951) 325–338 において、彼はすでに、ある単位的 C\*-環から別の単位的 C\*-環の上への線形等距写像 $\Phi$ は、ユニタリ $\Phi(1)$ とジョルダン $*$-同型の積であることを示していた。）

逆は直ちに従う：ジョルダン $*$-同型は自己共役元の平方 $h^2=h\circ h$ を平方に写し、
したがって正錐を正錐の上に写す。因子環 $M_n$ については、ジョルダン自己同型は正確に知られている：
$M_n(\mathbb C)$ のすべてのジョルダン $*$-自己同型は、ユニタリ $U$ により $A\mapsto UAU^*$ または $A\mapsto UA^TU^*$ と表される
（N. Jacobson, C. E. Rickart, *Jordan homomorphisms of rings*, Trans. Amer. Math. Soc. 69 (1950) 479–502; I. N.
Herstein, *Jordan homomorphisms*, Trans. Amer. Math. Soc. 81 (1956) 331–341）。ここでは $n=2$ の場合を証明する。

**命題 13.11.** $\operatorname{Herm}_2(\mathbb C)$ のすべてのジョルダン自己同型は、$U\in U(2)$ により $A\mapsto UAU^*$ または
$A\mapsto UA^TU^*$ と表される。

*証明.* 命題 13.5(a,b) により、ジョルダン自己同型は $O\in O(3)$ を用いて $a1+u\cdot\sigma\mapsto a1+(Ou)\cdot\sigma$ と作用する。
この意味での $A\mapsto UAU^*$ の行列を $R_U\in O(3)$ と書く。$R_{UV}=R_UR_V$ である。(i) $\det R_U=1$：
$U=W\operatorname{diag}(e^{i\theta_1},e^{i\theta_2})W^*$ は道
$U_t=W\operatorname{diag}(e^{it\theta_1},e^{it\theta_2})W^*$ によって $1$ と結ばれ、$\det R_{U_t}\in\{\pm1\}$ は $t$ に連続的に依存する。
(ii) すべての $R\in SO(3)$ はある $R_U$ である：命題 13.6 と同様に、
$U_1ZU_1^*=(Re_3)\cdot\sigma$、すなわち $R_{U_1}e_3=Re_3$ となる $U_1$ が存在する。すると $R_{U_1}^{-1}R\in SO(3)$ は $e_3$ を固定するので、$e_3$ のまわりのある角度 $\beta$ の回転である。直接計算により、$U_2=\operatorname{diag}(e^{-i\beta/2},e^{i\beta/2})$ に対し
$U_2XU_2^*=\cos\beta\,X+\sin\beta\,Y$ および $U_2YU_2^*=-\sin\beta\,X+\cos\beta\,Y$ となるので、$R=R_{U_1U_2}$ である。(iii)
転置は $X,Z$ を固定し $Y\mapsto-Y$ と写す：$O_T=\operatorname{diag}(1,-1,1)$、$\det O_T=-1$。よって
$\det O=1$ ならば $O=R_U$ であり、$\det O=-1$ ならば $OO_T=R_U$、すなわち自己同型は $A\mapsto UA^TU^*$ である。$\square$

**複合系.** 一つの因子環のジョルダン構造は、$T$ と恒等写像を区別できない。二つの系に対しては、この曖昧さが目に見えるものになる：

**命題 13.12.** $x,x'\in\operatorname{Herm}_n(\mathbb C)$、$y,y'\in\operatorname{Herm}_m(\mathbb C)$ とする。
$\mathbb C^n\otimes\mathbb C^m$ に作用する $\operatorname{Herm}_{nm}(\mathbb C)$ において
$$
(x\otimes y)\circ(x'\otimes y')=(x\circ x')\otimes(y\circ y')+\tfrac14[x,x']\otimes[y,y'].
$$
である。$\theta=\mathrm{id}\otimes T$（**部分転置**）、$A\circ_\theta B=\theta(\theta(A)\circ\theta(B))$ とおく。このとき
$\circ_\theta$ は $\operatorname{Herm}_{nm}(\mathbb C)$ 上のジョルダン積であって、$\theta$ を介して $\circ$ と同型であり、
局所オブザーバブル $x\otimes1$、$1\otimes y$ のすべての組の上で $\circ$ と一致し、
$(x\otimes y)\circ_\theta(x'\otimes y')=(x\circ x')\otimes(y\circ y')-\frac14[x,x']\otimes[y,y']$ を満たす。その平方の錐は
$\theta(\operatorname{Herm}_{nm}^+)$ である。$n=m=2$、$\Omega=e_1\otimes e_1+e_2\otimes e_2$ に対し、状態
$\omega=\frac12|\Omega\rangle\langle\Omega|\ge0$ はこの錐に属さない：$\theta(\omega)=\frac12F$ は固有値
$-\frac12$ をもつ。

*証明.* $j=x\circ x'$、$c=\frac12[x,x']$、$k=y\circ y'$、$d=\frac12[y,y']$ とおくと、$xx'=j+c$、$x'x=j-c$ であり、
$y$ についても同様であるから、$\frac12(xx'\otimes yy'+x'x\otimes y'y)=\frac12\big((j+c)\otimes(k+d)+(j-c)\otimes(k-d)\big)=j\otimes k+c\otimes d$ である。
$\theta$ は $\operatorname{Herm}_{nm}$ の実線形な対合である：それは $x\otimes y$ を $x\otimes y^T$ に写し、
$x\otimes y$ たちは $\mathbb R$ 上 $\operatorname{Herm}_{nm}$ を張る（$\operatorname{Herm}_n,\operatorname{Herm}_m$ の実基底の積 $h_\alpha\otimes g_\beta$ は $M_{nm}$ の複素基底をなし、エルミート行列はこの基底に関して実の係数をもつ）。したがって $\circ_\theta$ はジョルダン積であり、$\theta$ は同型である。$(y^T\circ y'^T)^T=y\circ y'$ かつ
$[y^T,y'^T]^T=-[y,y']$ であるから、$(x\otimes y^T)\circ(x'\otimes y'^T)$ の公式に $\theta$ を施すと
$\circ_\theta$ の公式が得られる。局所オブザーバブルの上では交換子の項は消える（$[1,1]=0$、$[x,1]=0$）。平方は
$A\circ_\theta A=\theta(\theta(A)^2)$ であり、$\theta^2=\mathrm{id}$ であるから、$\omega$ が $\theta(\operatorname{Herm}_{nm}^+)$ に属するのは $\theta(\omega)\ge0$ のときかつそのときに限る。
最後に $\theta(\omega)=\frac12\sum_{ij}E_{ij}\otimes E_{ji}=\frac12F$ である
（命題 13.9(b)。これは命題 6.1 の証明における計算である）。$\square$

二つの積は同じように良いジョルダン代数であり、一度に一つの部分系で測定できるものすべてについて一致するが、どの結合作用素が正であるか、すなわち複合系の状態が何であるかについては一致しない。両者を区別するのは交換子の項の符号であり、この符号は
$*$-代数 $M_n\otimes M_m$ によって定まり、因子のジョルダン代数によっては定まらない。（$\circ_\theta$ は
結合的な積 $A\star B=\theta(\theta(A)\theta(B))$ のジョルダン積であり、この積については $(x\otimes y)\star(x'\otimes y')=xx'\otimes y'y$ である：
すなわち代数 $M_n\otimes M_m^{\rm op}$ である。）錐 $\theta(\operatorname{Herm}_{nm}^+)$ こそが
**PPT 判定条件**（A. Peres, Phys. Rev. Lett. 77 (1996) 1413）を与えるものである：分離可能状態
$\rho=\sum_kp_ka_k\otimes b_k$（$p_k\ge0$、$a_k,b_k\ge0$）は両方の錐に属する。なぜなら
$\theta(\rho)=\sum_kp_ka_k\otimes b_k^T$ で $b_k^T\ge0$（命題 13.9(a)）であり、
$a\otimes b=(a^{1/2}\otimes b^{1/2})^2\ge0$ だからである。したがって $\theta(\operatorname{Herm}_4^+)$ に属さない $\omega$ は
エンタングルしている。ジョルダン代数のどのような複合が存在するかについての結果は、H. Hanche-Olsen, *On the structure and tensor
products of JC-algebras*, Canad. J. Math. 35 (1983) 1059–1074、および H. Barnum, M. Graydon, A. Wilce,
*Composites and categories of Euclidean Jordan algebras*, Quantum 4 (2020) 359 を参照せよ。

**ダイナミクス.** 時間発展が結合的な積を強制するのではないかと期待するかもしれない。そうはならない：

**命題 13.13.** $H\in\operatorname{Herm}_n(\mathbb C)$、$\delta(z)=i[H,z]$ とする。(a) $\delta$ は
$\operatorname{Herm}_n$ をそれ自身に写し、$\circ$ の導分である。$e^{t\delta}(z)=e^{itH}ze^{-itH}$ はジョルダン
自己同型である。(b) $x,y,z\in\operatorname{Herm}_n$ に対し $[L(x),L(y)]z=\frac14[[x,y],z]$。(c) ある $x_k,y_k\in\operatorname{Herm}_n$ に対し $\delta=4\sum_k[L(x_k),L(y_k)]$ となる：ハミルトン力学はジョルダン積だけによって生成される。

*証明.* (a) $(i[H,z])^*=-i[z,H]=i[H,z]$ である。$[H,\cdot]$ は結合的な積の導分であり、したがって $\circ$ の導分でもある。
(b) $4x\circ(y\circ z)=xyz+xzy+yzx+zyx$ から、$x,y$ を入れ替えた同じ式を引く。(c) (b) により、
$iH-\frac{\operatorname{Tr}H}n i1$（これは $iH$ と同じ交換子をもつ）を $[x_k,y_k]$ の和として書けば十分である。$i\ne j$ と
$a=E_{ij}+E_{ji}$、$b=i(E_{ij}-E_{ji})$ に対し、$[a,E_{ii}]=E_{ji}-E_{ij}$、$[b,E_{ii}]=-i(E_{ij}+E_{ji})$、
$[a,b]=2i(E_{jj}-E_{ii})$ である。これらは $\mathbb R$ 上でトレース零の反エルミート行列全体を張る。$\square$

ジョルダン代数が含んでいないのは、オブザーバブルとそれが生成するダイナミクスの間の **対応** $H\mapsto i[H,\cdot]$ である。たとえば $\operatorname{Herm}_n(\mathbb R)$ では、導分
$[L(x),L(y)]=\frac14[[x,y],\cdot]$ は実反対称な $[x,y]$ から来るが、これは実対称行列の $i$ 倍ではない。アルフセン（E. M. Alfsen）とシュルツ（F. W. Shultz）は、JBW-代数（形式的実ジョルダン代数の無限次元版。有限次元では二つの概念は一致する）がフォン・ノイマン環の自己共役部分であるための必要十分条件は、それがそのような「力学的対応」を許すことであることを示し、また C\*-環の状態空間を幾何学的に特徴づけた
（*Geometry of State Spaces of Operator Algebras*, Birkhäuser, 2003。証明なしに述べる）。対応は一意ではない：$\operatorname{Herm}_n(\mathbb C)$ 上では $H\mapsto i[H,\cdot]$ と $H\mapsto-i[H,\cdot]$ がともに条件を満たす。転置はこの二つを入れ替える。なぜなら $T(i[H,T(z)])=-i[H^T,z]$ だからである。そしてこれらはそれぞれ $M_n$ とその反対代数に属する——命題 13.12 の二つの符号である。

```python
x, y = sa(), sa()
T = jd.transpose_map
assert close(T(jp(x, y)), jp(T(x), T(y))) and close(T(x @ y), T(y) @ T(x))   # Proposition 13.9(a)
assert not close(T(x @ y), T(x) @ T(y))

k = 3
Mk = TypeIAlgebra([k], [k], precision='double')
units = torch.eye(k * k, dtype=dt).reshape(k * k, 1, k, k)                     # E_ij, batch index i*k + j
TE = T(Mk.operator(units)).matrix[:, 0]
choi = TE.reshape(k, k, k, k).permute(0, 2, 1, 3).reshape(k * k, k * k)       # block (i, j) = T(E_ij)
swap = torch.eye(k * k, dtype=dt).reshape(k, k, k, k).permute(0, 1, 3, 2).reshape(k * k, k * k)
w = torch.linalg.eigvalsh(choi)
assert torch.equal(choi, swap) and (w < 0).sum() == k * (k - 1) // 2          # Proposition 13.9(b)

sig = jd.anticommuting_symmetries(3)                                           # X, Y, Z
O_T = torch.stack([V3.from_matrix(s.T)[1:] for s in sig], dim=1)               # transpose on Bloch vectors
assert torch.allclose(O_T, torch.diag(torch.tensor([1.0, -1.0, 1.0], dtype=torch.float64)))

def ptranspose(R, n, m):                                                       # theta = id (x) T
    return R.reshape(n, m, n, m).permute(0, 3, 2, 1).reshape(n * m, n * m)
def herm(n):
    A = torch.randn(n, n, dtype=dt)
    return A + A.mH
jm = lambda A, B: 0.5 * (A @ B + B @ A)
cm = lambda A, B: A @ B - B @ A
a, a2, b, b2 = herm(2), herm(2), herm(3), herm(3)
local = torch.kron(jm(a, a2), jm(b, b2))
assert close(jm(torch.kron(a, b), torch.kron(a2, b2)), local + 0.25 * torch.kron(cm(a, a2), cm(b, b2)))
theta = lambda R: ptranspose(R, 2, 3)
assert close(theta(jm(theta(torch.kron(a, b)), theta(torch.kron(a2, b2)))),
             local - 0.25 * torch.kron(cm(a, a2), cm(b, b2)))                  # Proposition 13.12
Om = torch.tensor([1.0, 0, 0, 1.0], dtype=dt) / math.sqrt(2)
omega = torch.outer(Om, Om.conj())                                             # Bell state
print("spectrum of the partial transpose of the Bell state:", torch.linalg.eigvalsh(ptranspose(omega, 2, 2)))
assert abs(torch.linalg.eigvalsh(ptranspose(omega, 2, 2))[0] + 0.5) < 1e-12

Hm, z1, z2 = herm(3), herm(3), herm(3)
delta = lambda A: 1j * cm(Hm, A)
assert close(delta(jm(z1, z2)), jm(delta(z1), z2) + jm(z1, delta(z2)))         # Proposition 13.13(a)
assert close(jm(z1, jm(z2, Hm)) - jm(z2, jm(z1, Hm)), 0.25 * cm(cm(z1, z2), Hm))   # (b)
```

> **物理的意味.** スピンをもたない粒子に対しては、転置（実基底における複素共役）は **時間反転** である。
> それは一つの系のオブザーバブルの対称性であるが、反ユニタリであり、一つの部分系だけに施すと物理的な操作ではない：
> 命題 13.12 は、それがベル状態を負の固有値をもつ作用素に変えることを示している。量子力学が複素 $*$-代数を用いるのは、
> 複合系が対称化された積だけでなく、符号つきの交換子を必要とするからである。

## 13.6 トレース不等式に戻って

命題 11.9（実験 3）は、$X=X^*$ と $Y\ge0$ に対して $\operatorname{Tr}|YXY|\le\operatorname{Tr}(Y|X|Y)$ を述べている（ここで
$X,Y$ は $\mathcal M$ の作用素であり、パウリ行列ではない）。その証明は $Y\ge0$ を一度も用いていない：
$T=X_\pm^{1/2}Y$ による $YX_\pm Y=T^*T$ は $Y=Y^*$ しか必要とせず、補題 11.8 は $YX_+Y$ と $YX_-Y$ に適用される。したがってこの不等式はすべての $Y=Y^*$ に対して成り立ち、
ジョルダンの言葉では
$$
\operatorname{tr}|P(y)x|\le\operatorname{tr}P(y)|x|,
$$
となる。これはすべてのユークリッド的ジョルダン代数で意味をもつ。これを $V_m$ において、三つの性質だけを用いる議論によって証明する。

**命題 13.14（スピン因子における実験 3）.** $x,y\in V_m$ に対し
$\operatorname{tr}|P(y)x|\le\operatorname{tr}P(y)|x|$ であり、$x\in Q$ または $-x\in Q$ ならば等号が成り立つ。

*証明.* $f_+(t)=\max(t,0)$、$f_-(t)=\max(-t,0)$ として $x_\pm=f_\pm(x)$ とおく（定理 13.4(3)）：$x_\pm\in Q$、$x=x_+-x_-$、
$|x|=x_++x_-$ である。$A=P(y)x_+$、$B=P(y)x_-$ とおく。定理 13.4(5) によりこれらは $Q$ に属し、$P(y)x=A-B$ である。
$P(y)x=\mu_1c_1+\mu_2c_2$ をスペクトル分解とし（定理 13.4(2)。$c_1=c_w$、$c_2=c_{-w}$ であるから $c_1+c_2=1$、
$c_1\circ c_2=0$）、$c$ を $\mu_i>0$ なる $c_i$ の和とし、
$t>0$ のとき $\operatorname{sgn}'(t)=1$、それ以外のとき $-1$ として $s=2c-1=\sum_i\operatorname{sgn}'(\mu_i)c_i$ とおく。すると
$s\circ P(y)x=\sum_i|\mu_i|c_i=|P(y)x|$ であり（定理 13.4(2,3) により、
$\mu_1=\mu_2$ のときもこれは $w$ の選び方によらない）、$c,1-c\in Q$ である（冪等元は平方である）。自己双対性の易しい半分（定理 13.4(4)：
$a,b\in Q$ に対し $\operatorname{tr}(a\circ b)\ge0$）を $\operatorname{tr}(A\circ(1-c))$ と $\operatorname{tr}(B\circ c)$ に適用すると
$$
\operatorname{tr}|P(y)x|=\operatorname{tr}(A\circ c)-\operatorname{tr}(A\circ(1-c))-\operatorname{tr}(B\circ c)+\operatorname{tr}(B\circ(1-c))
\le\operatorname{tr}A+\operatorname{tr}B=\operatorname{tr}P(y)|x|.
$$
を得る。$x\in Q$ ならば $x_-=0$、$B=0$、$|P(y)x|=A$ であるから、両辺はともに $\operatorname{tr}A$ に等しい。$-x\in Q$ のときも同様である。
$\square$

証明が用いるのは次のことだけである：(i) すべての $y$ に対し $P(y)Q\subseteq Q$。(ii) スペクトル定理：すべての $z$ は、和が $1$ となる直交冪等元 $c_i$ を用いて
$\sum_i\mu_ic_i$ と書け、$f(z)=\sum_if(\mu_i)c_i$ である（したがって $x_\pm$ と冪等元は平方である）。(iii) $\operatorname{tr}$ は線形であり、$\operatorname{tr}c_i\ge0$ として $\operatorname{tr}(\sum_i\mu_ic_i)=\sum_i\mu_i\operatorname{tr}c_i$ であり、$a,b\in Q$ に対し $\operatorname{tr}(a\circ b)\ge0$ である。階数 $r$ のユークリッド的ジョルダン代数では、(ii) は $r$ 個の原始冪等元 $c_i$（ジョルダン枠）について成り立ち、トレースは
$\operatorname{tr}z=\sum_i\mu_i$ であり、形式 $\operatorname{tr}(x\circ y)$ は結合的、すなわち
$\operatorname{tr}((x\circ y)\circ z)=\operatorname{tr}(x\circ(y\circ z))$ であり、そして (i) が成り立つ [FK94, Chapters II–III]。
これらの事実はここでは証明しない。すると (iii) の最後の部分は (i) から従う：$P(a^{1/2})=2L(a^{1/2})^2-L(a)$ を展開し、結合性を用いると $\operatorname{tr}(a\circ b)=\operatorname{tr}P(a^{1/2})b$ となり、これは $Q$ の元のトレースであるから
$\ge0$ である。したがってこの不等式はすべてのユークリッド的ジョルダン代数で成り立つ——複素 $*$-代数の自己共役部分とジョルダン同型でない
$\operatorname{Herm}_r(\mathbb H)$ やアルバート代数も含めて。$M_{sa}$ では、$\operatorname{Tr}$ と $P(y)x=yxy$ を用いた同じ議論が命題 11.9 を再証明する。$V_m$ では、$\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$ かつ $\phi(|x|)=|\phi(x)|$ であるから、この不等式は $\phi$ を通して命題 11.9 からも従う。コードは両方を確かめる：

```python
for m in (2, 4, 7):
    Vm = SpinFactor(m)
    xs, ys = Vm.random(20_000), Vm.random(20_000)
    lhs = Vm.trace(Vm.apply_function(Vm.quadratic_representation(ys, xs), torch.abs))
    rhs = Vm.trace(Vm.quadratic_representation(ys, Vm.apply_function(xs, torch.abs)))
    definite = Vm.determinant(xs) >= 0                                 # x or -x in the cone
    print(f"m = {m}: min gap {(rhs - lhs).min().item():.1e}, max gap {(rhs - lhs).max().item():.1f}, "
          f"definite x: {definite.double().mean().item():.3f}")
    assert (rhs - lhs >= -1e-12 * rhs).all()                            # Proposition 13.14
    assert ((rhs - lhs)[definite].abs() <= 1e-12 * rhs[definite]).all()
    Xs, Ys = Vm.to_matrix(xs[:200]), Vm.to_matrix(ys[:200])
    tn = torch.linalg.svdvals(Ys @ Xs @ Ys).sum(-1)                     # Tr |phi(y) phi(x) phi(y)|
    assert (tn - Vm.representation_dim / 2 * lhs[:200]).abs().max() < 1e-10 * tn.max()
```

差 $\operatorname{tr}P(y)|x|-\operatorname{tr}|P(y)x|$ は定符号の $x$ に対しては零となり、典型的な不定符号の $x$ に対しては正である（これらのガウス・サンプルについて、表示される最大値は $10^1$–$10^2$ の程度である）。

## 13.7 比較の対象

我々の知る限り、ジョルダン代数のための主流の数値ライブラリは存在しない。本章の演算は通常、行列を用いて手で書かれる。代数的な側面は計算機代数がカバーしている：SageMath の `JordanAlgebra` は特殊ジョルダン代数、対称双線形形式のジョルダン代数（スピン因子を含む）、そして八元数代数から 27 次元の例外ジョルダン代数を構成する。Magma にもジョルダン代数のコンストラクタがあり、例外的なものも含まれる。これらは厳密な代数的道具であって、バッチ化された数値的道具ではない。ローレンツ錐は最適化では標準的である：SOCP ソルバー（ECOS、SCS、Clarabel、MOSEK。cvxpy の `SOC` 制約を通じて利用できる）は、$V_m$ のスペクトル分解（錐への射影のため）またはその二次表現（内点法のネステロフ–トッド・スケーリングのため）を用い、半正定値計画のソルバーは $\operatorname{Herm}_n$ について同じことを行う。本ライブラリが付け加えるのは、$\bigoplus_cM_{k_c}$ 上のバッチ化された微分可能な作用素としてのジョルダン積、二次表現、転置、関数計算と行列表現を備えたスピン因子、そして与えられた行列が生成するジョルダン代数の乱択計算と、それらが生成する $*$-代数との比較である。

## 演習問題

1. （**ローレンツ錐の等質性.**）$y\in V_m$、$\det y\ne0$ とする。$y^{-1}=(a,-u)/\det y$ が
   $y\circ y^{-1}=1$ を満たすこと、$P(y)P(y^{-1})=\mathrm{id}$ であること、$P(y)$ が $Q$ の内部をそれ自身の上に写すこと、そして内部の $x$ に対し
   $P(x^{1/2})1=x$ であることを示せ。$Q$ の内部が等質であることを導け。
2. （**等方的状態.**）$p\in[0,1]$ に対し、$\mathbb C^2\otimes\mathbb C^2$ 上で $\rho_p=p\,\omega+(1-p)\frac14$ とする
   （命題 13.12）。$\theta(\rho_p)$ のスペクトルを計算し、$p>\frac13$ のとき $\rho_p$ がエンタングルしていることを示せ。
3. （**境界点の積.**）$\mathbb R^m$ において $|w|=|w'|=1$ とする。$(1,w)\circ(1,w')\in Q$ であるための必要十分条件が
   $w'=w$ または $w'=-w$ であることを示せ：ローレンツ錐の二本の境界半直線のジョルダン積はほとんど決して錐に属さない。
4. （**実量子力学.**）$x,y\in\operatorname{Herm}_n(\mathbb R)$ に対し、$\operatorname{Herm}_n(\mathbb R)$ の導分
   $[L(x),L(y)]$ は $K$ を実反対称行列として $z\mapsto[K,z]$ であること、そして
   $n\ge2$ のとき、そのような写像は $0$ でない限り決して $H\in\operatorname{Herm}_n(\mathbb R)$ による $z\mapsto i[H,z]$ にはならないことを示せ。

<details>
<summary>解答</summary>

1. $(a,u)\circ(a,-u)=(a^2-|u|^2,0)$ であり、これから $y\circ y^{-1}=1$ を得る。$\phi(y)$ と $\phi(y^{-1})$ は可換であり（ともに
   $\operatorname{span}\{1,u\cdot s\}$ に属する）、したがって $\phi(y)\phi(y^{-1})=\phi(y)\circ\phi(y^{-1})=1$ であり、定理 13.4(5) により
   $\phi(P(y)P(y^{-1})x)=\phi(x)$ である。$\phi$ は単射である。$P(y)$ は $Q$ を $Q$ の中に写し、その逆 $P(y^{-1})$ も同様であるから、
   $P(y)Q=Q$ である。線形同相写像であるから、それは内部を内部の上に写す。最後に
   $P(y)1=2y\circ y-y^2=y^2$ であるから $P(x^{1/2})1=x$ であり、$x$ が内部の点ならば $P(x^{1/2})$ は可逆である
   （$\det x^{1/2}=\sqrt{\lambda_+\lambda_-}>0$）。内部の任意の二点は
   $P(x'^{1/2})P(x^{1/2})^{-1}$ によって結ばれる。

```python
V4 = SpinFactor(4)
xi = V4.random(100, cone=True)                                   # interior points (almost surely)
r = V4.sqrt(xi)
assert (V4.quadratic_representation(r, V4.unit((100,))) - xi).abs().max() < 1e-10
yv = V4.random(100)
back = V4.quadratic_representation(V4.inverse(yv), V4.quadratic_representation(yv, xi))
assert (back - xi).abs().max() < 1e-8 * xi.abs().max()
assert V4.in_cone(V4.quadratic_representation(yv, xi)).all()
```

2. 命題 13.12 により $\theta(\rho_p)=\frac p2F+\frac{1-p}4$ である。$F$ は固有値 $1$（三重）と $-1$（一重）をもつので、
   スペクトルは $\frac{1+p}4$（三重）と $\frac{1-3p}4$ であり、後者が負となるのは $p>\frac13$ のときかつそのときに限る。そのとき PPT 判定条件により $\rho_p$ は分離可能でない。（$2\times2$ 系に対しては逆も成り立つ：M., P. and R. Horodecki, Phys.
   Lett. A 223 (1996) 1。ここでは証明しない。）

```python
for pp in (0.2, 1 / 3, 0.5):
    rho_p = pp * omega + (1 - pp) * torch.eye(4, dtype=dt) / 4
    lam_min = torch.linalg.eigvalsh(ptranspose(rho_p, 2, 2))[0].item()
    assert abs(lam_min - (1 - 3 * pp) / 4) < 1e-12
```

3. $t=w\cdot w'\in[-1,1]$ として $(1,w)\circ(1,w')=(1+t,\,w+w')$ であり、$|w+w'|=\sqrt{2+2t}$ である。これが $Q$ に属するのは
   $1+t\ge\sqrt{2(1+t)}$ のとき、すなわち $t=-1$（積は $0$）または $1+t\ge2$、すなわち $t=1$ のときかつそのときに限る。

```python
V3c = SpinFactor(3)
wv = torch.nn.functional.normalize(torch.randn(2, 1000, 3, dtype=torch.float64), dim=-1)
prod = V3c.product(V3c.element(1.0, wv[0]), V3c.element(1.0, wv[1]))
assert not V3c.in_cone(prod).any()
assert V3c.in_cone(V3c.product(V3c.element(1.0, wv[0]), V3c.element(1.0, wv[0]))).all()
```

4. 命題 13.13(b) により $K=\frac14[x,y]$ として $[L(x),L(y)]z=[K,z]$ であり、$K^T=\frac14[y^T,x^T]=-K$ は実である。実の $H$
   と $z$ に対し $i[H,z]$ は純虚数であり、一方 $[K,z]$ は実であるから、両者が一致するのはともに零となるときに限る。したがって実量子力学では、ダイナミクスの生成元 $K$ はオブザーバブルではない。

</details>

**この先へ.** 本章は本講義の最終章である。有限次元の理論の標準的な文献は [FK94] である。K. McCrimmon, *A Taste of Jordan Algebras* (Springer 2004) は代数的側面と歴史を扱い、H. Hanche-Olsen and E. Størmer, *Jordan Operator Algebras* (Pitman 1984) と Alfsen–Shultz (2003) は無限次元の理論と C\*-環との関係を扱う。操作的公理からの量子論の再構成における複合系の役割については、上で引用した Barnum と Wilce の論文を参照せよ。第11章を越えるトレース不等式については E. Carlen, *Trace inequalities and quantum entropy* (Contemp. Math. 529, 2010) を参照せよ。二次表現とローレンツ錐は錐最適化に再び現れる：Alizadeh–Goldfarb (2003) と、そこで概観されているジョルダン代数的内点法を参照せよ。
