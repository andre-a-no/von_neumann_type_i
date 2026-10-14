[English](../E_jordan_classification.md) · [Español](../es/E_jordan_classification.md) · [Français](../fr/E_jordan_classification.md) · [Deutsch](../de/E_jordan_classification.md) · [中文](../zh/E_jordan_classification.md) · **日本語** · [Русский](../ru/E_jordan_classification.md)

# 付録 E. ユークリッド的ジョルダン代数の分類

第13章では、有限次元の形式的実ジョルダン代数の Jordan–von Neumann–Wigner による分類を用いた。
本付録では、ジョルダン恒等式と有限次元の線形代数だけからこれを完全に証明する。道筋は
古典的なものである（J. Faraut, A. Korányi, *Analysis on Symmetric Cones*, Chapters II–V；K. McCrimmon, *A Taste of
Jordan Algebras*）：冪結合性とスペクトル補題により冪等元と結合的で正値なトレース形式が得られる；
ジョルダン枠に関するピアス分解は代数を部分 $J_{ij}$ に分け、それらの積は
厳密な規則に従う；階数 $r\ge3$ の単純代数では、これらの規則がピアス空間 $J_{13}$ を合成代数にし、それは
フルヴィッツの定理により $\mathbb R,\mathbb C,\mathbb H$ または $\mathbb O$ であり、代数全体はその上の
$\operatorname{Herm}_r$ となる。E.6節では実結合的可除代数に関するフロベニウスの定理も証明し、
E.8節ではアルバート代数 $\operatorname{Herm}_3(\mathbb O)$ が、いかなる結合代数にも埋め込めない
ジョルダン代数であることを証明する。

**本付録で学ぶこと**：冪結合性、トレース形式の正値性、単位元の存在が形式的実性だけから従う理由、
ピアス規則と、それらが単純代数の座標に数のような掛け算を強いる仕組み、
ケーリー–ディクソンの倍化によるフルヴィッツの定理の証明、そして八元数が新しい単純代数を与えるのは階数 $3$ のときだけである理由。

**規約.** 代数はすべて、双線形な積をもつ有限次元実ベクトル空間とする。**ジョルダン代数**
は 13.1節と同様である：$x\circ y=y\circ x$ かつ $(x^2\circ y)\circ x=x^2\circ(y\circ x)$、$x^2=x\circ x$ であり、単位元は
仮定しない。$L(x)y=x\circ y$、$[A,B]=AB-BA$ とし、$\operatorname{Tr}$ は線形写像のトレースである。**部分代数** とは
$\circ$ について閉じた部分空間のこと、**イデアル** とは $J\circ I\subseteq I$ を満たす部分空間 $I$ のことであり、
$J\ne0$ でそのイデアルが $0$ と $J$ だけであるとき $J$ は **単純** であるという。**冪等元** とは $c^2=c$ を満たす $c$ のことであり、
$c\circ c'=0$ のとき冪等元 $c,c'$ は **直交する** という。$\sum_ix_i^2=0$ ならばすべての $x_i=0$ となるとき $J$ は **形式的実** であるといい、
単位元 $1$ と $\langle x\circ y,z\rangle=\langle y,x\circ z\rangle$ を満たす内積をもつとき **ユークリッド的** であるという。
多項式恒等式は次のように **線形化** する：$p(t)=\sum_kp_kt^k$ がベクトル係数をもち、すべての
実数 $t$ に対し $p(t)=0$ ならば、すべての $p_k=0$ である（$\deg p+1$ 個の相異なる点で値をとり、ファンデルモンド行列の逆を掛ければよい）。

## E.1 線形化されたジョルダン恒等式と冪結合性

**補題 E.1.** すべてのジョルダン代数において、すべての $x,y,z,w$ に対し：

1. $[L(x),L(x^2)]=0$ かつ $[L(x^2),L(z)]+2[L(x\circ z),L(x)]=0$；
2. $L(x^2\circ y)=2L(x\circ y)L(x)+L(x^2)L(y)-2L(x)L(y)L(x)$；
3. $L((x\circ w)\circ y)=L(w\circ y)L(x)+L(x\circ y)L(w)+L(x\circ w)L(y)-L(w)L(y)L(x)-L(x)L(y)L(w)$；
4. $L((x\circ w)\circ y)-L(x\circ(w\circ y))=[L(w),[L(x),L(y)]]$；特に
   $\operatorname{Tr}L((x\circ w)\circ y)=\operatorname{Tr}L(x\circ(w\circ y))$。

*証明.* (1) ジョルダン恒等式は $L(x)L(x^2)y=L(x^2)L(x)y$ ということである。$x$ を $x+tz$ で置き換えると、$t$ の係数は
$[L(z),L(x^2)]+2[L(x),L(x\circ z)]=0$ である。(2) ジョルダン恒等式そのものにおいて $x$ を $x+tz$ で置き換えると、$t$ の係数は
$$
2((x\circ z)\circ y)\circ x+(x^2\circ y)\circ z=2(x\circ z)\circ(y\circ x)+x^2\circ(y\circ z).
$$
となる。$z$ の線形写像として読めば、四つの項はそれぞれ $2L(x)L(y)L(x)z$、$L(x^2\circ y)z$、$2L(x\circ y)L(x)z$、
$L(x^2)L(y)z$ である。(3) (2) において $x$ を $x+sw$ で置き換え、$s$ の係数をとって $2$ で割る。(4) (3) を
そのままの形で一度、$x$ と $y$ を入れ替えてもう一度書く。右辺の最初の三項は一致し、
$(y\circ w)\circ x=x\circ(w\circ y)$ であるから、差は
$-L(w)L(y)L(x)-L(x)L(y)L(w)+L(w)L(x)L(y)+L(y)L(x)L(w)=[L(w),[L(x),L(y)]]$ という交換子であり、そのトレースは $0$ である。
$\square$

**補題 E.2（冪結合性）.** $x^1=x$、$x^{k+1}=x\circ x^k$ とおく。このときすべての
$p,q\ge1$ に対し $x^p\circ x^q=x^{p+q}$ であり、作用素 $L(x^k)$ は互いに可換である。したがって、定数項をもたない実多項式 $f,g$ と
$f(x):=\sum_kf_kx^k$ に対し $(fg)(x)=f(x)\circ g(x)$ であり、$\mathbb R[x]=\{f(x)\}$ は結合的かつ可換な部分代数である。

*証明.* 補題 E.1(1) により $L(x)$ と $L(x^2)$ は可換である。$\mathcal A$ をそれらが生成する可換な作用素の
代数とする。$n\ge2$ に関する帰納法で次を示す：$p+q\le n$ ならば $x^p\circ x^q=x^{p+q}$ であり、
$k\le n$ に対し $L(x^k)\in\mathcal A$ である。$n=2$ のときこれは定義である。$n$ について成り立つと仮定し、$p+q=n+1$ とする。$q=1$ ならば
可換性により $x^p\circ x=x^{n+1}$ である。$q\ge2$ ならば
$x^p\circ x^q=L(x^p)L(x)x^{q-1}=L(x)L(x^p)x^{q-1}=L(x)x^{n}=x^{n+1}$ である。次に、補題 E.1(2) を $y=x^{n-1}$ に適用し、
$x^2\circ x^{n-1}=x^{n+1}$ を用いると、$L(x^{n+1})=2L(x^n)L(x)+L(x^2)L(x^{n-1})-2L(x)^2L(x^{n-1})\in\mathcal A$ を得る。残りは
双線形性から従う。$\square$

## E.2 冪等元とピアス分解

**補題 E.3（スペクトル補題）.** $J$ を形式的実とし、$x\in J$ とする。このとき $k\ge0$、相異なる 零でない実数
$\lambda_1,\dots,\lambda_k$、および 零でない互いに直交する冪等元 $c_1,\dots,c_k\in\mathbb R[x]$ が存在して
$$
x=\sum_{i=1}^k\lambda_ic_i,\qquad x^n=\sum_{i=1}^k\lambda_i^nc_i\quad(n\ge1),\qquad\mathbb R[x]=\operatorname{span}\{c_1,\dots,c_k\}.
$$
となる。特に、零でない形式的実ジョルダン代数はすべて 零でない冪等元を含む。

*証明.* 項が一つの場合の形式的実性は、$z^2=0$ ならば $z=0$ ということである。$K=\{q\in\mathbb R[t]:(tq)(x)=0\}$ とおく。
補題 E.2 により $K$ は $\mathbb R[t]$ のイデアルである：$q\in K$ かつ $s=s_0+s_1$、$s_1(0)=0$ ならば
$(tqs)(x)=s_0(tq)(x)+(tq)(x)\circ s_1(x)=0$ である。$x,x^2,\dots,x^{\dim J+1}$ は線形従属であるから $K\ne0$ であり、
モニック多項式 $\mu$ により $K=\mu\,\mathbb R[t]$ と書ける。(i) $\mu(0)\neq0$：$\mu=th$ ならば $w=(th)(x)$ は
$\mu h\in K$ であるから $w^2=(t\cdot\mu h)(x)=0$ を満たす；よって $w=0$、$h\in K$ となり $\mu\mid h$ であるが、これは次数から不可能である。(ii) $\mu$ は
重複因子をもたない：$\deg g\ge1$ として $\mu=g^2h$ ならば、$w=(tgh)(x)$ は $w^2=(t\cdot\mu\cdot th)(x)=0$ を満たすので $gh\in K$ となるが、
これは次数から不可能である。(iii) $\mu$ は実でない根をもたない：$\beta\ne0$ として $\mu=((t-\alpha)^2+\beta^2)h$ ならば、
$a=(t(t-\alpha)h)(x)$ と $b=(\beta th)(x)$ は $a^2+b^2=(t\cdot\mu\cdot th)(x)=0$ を満たすので $b=0$ かつ $h\in K$ となり、不可能である。
したがって、相異なる 零でない実数 $\lambda_i$ により $\mu=\prod_{i=1}^k(t-\lambda_i)$ である。
$$
p_i(t)=\frac t{\lambda_i}\prod_{j\ne i}\frac{t-\lambda_j}{\lambda_i-\lambda_j},\qquad c_i=p_i(x).
$$
とおく。多項式 $p_ip_l-\delta_{il}p_i$ は、$t$ とすべての $\lambda_j$ で消える多項式との積であり、したがって $t$ と
$K$ の元との積である；よって $c_i\circ c_l=\delta_{il}c_i$ である。同様に $t-\sum_i\lambda_ip_i=t\big(1-\sum_i\lambda_ip_i/t\big)$ であって
括弧の中はすべての $\lambda_j$ で消えるので、$x=\sum_i\lambda_ic_i$ である。$c_i=0$ ならば $p_i/t\in K$ であり、その次数は $k-1<\deg\mu$ である；
よって $c_i\ne0$ である。$x^n$ の公式は直交性から従い、それにより $\mathbb R[x]\subseteq\operatorname{span}\{c_i\}$ がわかる。
$\square$

**命題 E.4（一つの冪等元に関するピアス分解）.** $c$ をジョルダン代数 $J$ の冪等元とする。

1. $L(c)(2L(c)-1)(L(c)-1)=0$。したがって
   $J_\alpha(c)=\{x:c\circ x=\alpha x\}$ として $J=J_1(c)\oplus J_{1/2}(c)\oplus J_0(c)$ であり、$c\in J_1(c)$ である。
2. $z\in J_0(c)\cup J_1(c)$ に対し $[L(c),L(z)]=0$ であり、$y\in J_{1/2}(c)$ に対し $[L(c),L(y^2)]=0$ である。
3. $J_1\circ J_1\subseteq J_1$、$J_0\circ J_0\subseteq J_0$、$J_1\circ J_0=0$、$(J_1+J_0)\circ J_{1/2}\subseteq J_{1/2}$
   かつ $J_{1/2}\circ J_{1/2}\subseteq J_1+J_0$（いずれも $c$ に関して）。
4. $c'$ が $c$ に直交する冪等元ならば、$c+c'$ は冪等元であり、$[L(c),L(c')]=0$ かつ
   $J_1(c')\subseteq J_0(c)$ である。

*証明.* (1) 補題 E.1(2) を $x=y=c$ に適用すると $L(c)=3L(c)^2-2L(c)^3$ を得る。単根のみをもつ多項式で
零化される線形写像は対角化可能であり、その固有値はそれらの根のうちにある。(2) 補題 E.1(1) を $x=c$ に適用すると
$[L(c),L(z)]=2[L(c),L(c\circ z)]$、すなわち $z\in J_\alpha(c)$ に対し $(1-2\alpha)[L(c),L(z)]=0$ を得る。補題 E.1(1) を $x=y$、
$z=c$ に適用すると $[L(y^2),L(c)]+2[L(y\circ c),L(y)]=[L(y^2),L(c)]+[L(y),L(y)]=[L(y^2),L(c)]=0$ を得る。(3) (2) により、$z\in J_0\cup J_1$ ならば $L(z)$ は
$L(c)$ の各固有空間をそれ自身に写す。これにより一番目、二番目、四番目の規則と
$J_1\circ J_0\subseteq J_1\cap J_0=0$ が得られる。$y\in J_{1/2}$ に対し $L(y^2)$ は $L(c)$ と可換であるから、$c\circ y^2=y^2\circ c$ は
$c$ と同様に $L(c)$ の固有値 $1$ の固有ベクトルである；$y^2=a_1+a_{1/2}+a_0$ と書くと $c\circ y^2=a_1+\frac12a_{1/2}$ が $J_1$ に属するので、
$a_{1/2}=0$ である。偏極化（$y\to y+y'$）により最後の規則が得られる。(4) $(c+c')^2=c+2c\circ c'+c'=c+c'$ であり、
$c'\in J_0(c)$ であるから、(2) により $L(c)$ と $L(c')$ は可換で同時対角化可能である；特に $J_1(c')$ は
共通の固有ベクトルで張られる。$L(c)$ に対する固有値 $\alpha$ と $L(c')$ に対する固有値 $1$ をもつ
共通の固有ベクトルの上で、$L(c+c')$ は固有値 $\alpha+1\in\{0,\frac12,1\}$ をもつので
$\alpha=0$ である。$\square$

## E.3 トレース形式；ユークリッド的 = 形式的実

**定理 E.5.** $J$ を形式的実とし、$\tau(x,y)=\operatorname{Tr}L(x\circ y)$ とする。

1. $\tau$ は対称、双線形かつ **結合的** である：$\tau(x\circ y,z)=\tau(x,y\circ z)$。
2. $\tau$ は正定値である；零でないすべての冪等元 $c$ に対し $\tau(c,c)=\operatorname{Tr}L(c)\ge1$ である。
3. $J$ は単位元をもつ。
4. $J$ は内積 $\tau$ に関してユークリッド的である。逆に、すべてのユークリッド的ジョルダン代数は形式的実である。

*証明.* (1) 補題 E.1(4)。(2) 命題 E.4(1) により、
$0\ne c\in J_1(c)$ であるから $\operatorname{Tr}L(c)=\dim J_1(c)+\frac12\dim J_{1/2}(c)\ge1$ である。補題 E.3 のように $x=\sum_i\lambda_ic_i$ と書くと $x^2=\sum_i\lambda_i^2c_i$ であり、
$\tau(x,x)=\operatorname{Tr}L(x^2)=\sum_i\lambda_i^2\operatorname{Tr}L(c_i)$ は $k=0$、すなわち $x=0$ でない限り $>0$ である。
(3) $J\neq0$ とし、$e$ を $\dim J_1(e)$ が最大となる 零でない冪等元（補題 E.3）とする。$J_0(e)\ne0$ ならば、それは
形式的実な部分代数（命題 E.4(3)）であり、零でない冪等元 $f$ を含む；$f$ は $e$ に直交するので、
$e+f$ は冪等元であり、$J_1(e)\subseteq J_1(e+f)$（$f\circ J_1(e)\subseteq J_0\circ J_1=0$ だから）かつ
$f\in J_1(e+f)\setminus J_1(e)$ となって、最大性に矛盾する。よって $J_0(e)=0$ である。すると $y\in J_{1/2}(e)$ に対し
$y^2\in J_1(e)$ であり、結合性により
$\frac12\tau(y,y)=\tau(e\circ y,y)=\tau(e,y^2)=\operatorname{Tr}L(e\circ y^2)=\operatorname{Tr}L(y^2)=\tau(y,y)$ となるので、$y=0$ である。
したがって $J=J_1(e)$ であり、$e$ は単位元である。(4) (1)–(3) により $\tau(x\circ y,z)=\tau(y,x\circ z)$ である。逆に、$J$ がユークリッド的で
$\sum_ix_i^2=0$ ならば、$0=\sum_i\langle x_i\circ x_i,1\rangle=\sum_i\langle x_i,x_i\rangle$ である。$\square$

したがって「形式的実」と「ユークリッド的」は有限次元代数の同じクラスである。以下、$J$ は形式的実であり、
$1$ はその単位元、$\tau$ はそのトレース形式である。以下のコードは、構造定数で与えた $\operatorname{Herm}_3(\mathbb C)$ において、補題 E.1(4)、補題 E.3（冪等元 $p_i(x)$ は
ジョルダン冪だけで計算する）、命題 E.4、定理 E.5 を検証する；
それを構成するのに用いるケーリー–ディクソン積 `cd` は E.6節のものである。

```python
import math
import torch

torch.manual_seed(0)
dt = torch.float64

def conj(x):                                      # Cayley-Dickson conjugation (a, b)* = (a*, -b) on R^(2^k)
    if x.shape[-1] == 1:
        return x
    h = x.shape[-1] // 2
    return torch.cat([conj(x[..., :h]), -x[..., h:]], -1)

def cd(x, y):                                     # (a, b)(c, d) = (ac - d*b, da + bc*): R, C, H, O for 1, 2, 4, 8
    if x.shape[-1] == 1:
        return x * y
    h = x.shape[-1] // 2
    a, b, c, d = x[..., :h], x[..., h:], y[..., :h], y[..., h:]
    return torch.cat([cd(a, c) - cd(conj(d), b), cd(d, a) + cd(b, conj(c))], -1)

def herm(n, d):
    """Herm_n(D), D = R, C, H, O for d = 1, 2, 4, 8; coordinates: diagonal, then the entries a_ij, i < j."""
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    def to_mat(v):
        X = torch.zeros(v.shape[:-1] + (n, n, d), dtype=dt)
        X[..., range(n), range(n), 0] = v[..., :n]
        for p, (i, j) in enumerate(pairs):
            X[..., i, j, :] = v[..., n + d * p:n + d * (p + 1)]
            X[..., j, i, :] = conj(X[..., i, j, :])
        return X
    def from_mat(X):
        return torch.cat([X[..., range(n), range(n), 0]] + [X[..., i, j, :] for i, j in pairs], -1)
    def matmul(X, Y):
        return torch.stack([torch.stack([sum(cd(X[..., i, k, :], Y[..., k, j, :]) for k in range(n))
                                         for j in range(n)], -2) for i in range(n)], -3)
    def prod(x, y):
        X, Y = to_mat(x), to_mat(y)
        return from_mat(0.5 * (matmul(X, Y) + matmul(Y, X)))
    return n + d * len(pairs), prod

def spin(m):                                      # V_m: (a, u) o (b, v) = (ab + u.v, av + bu)
    def prod(x, y):
        return torch.cat([x[..., :1] * y[..., :1] + (x[..., 1:] * y[..., 1:]).sum(-1, keepdim=True),
                          x[..., :1] * y[..., 1:] + y[..., :1] * x[..., 1:]], -1)
    return m + 1, prod

def structure(dim, prod):                         # C[a, b] = e_a o e_b
    E = torch.eye(dim, dtype=dt)
    return prod(E[:, None, :].expand(dim, dim, dim), E[None, :, :].expand(dim, dim, dim))

mul = lambda C, x, y: torch.einsum('...a,...b,abc->...c', x, y, C)
Lop = lambda C, x: torch.einsum('...a,abc->...cb', x, C)          # matrix of L(x)
comm = lambda A, B: A @ B - B @ A

C3 = structure(*herm(3, 2))                       # Herm_3(C), dimension 9
x, y, w = torch.randn(3, 9, dtype=dt)
lhs = Lop(C3, mul(C3, mul(C3, x, w), y)) - Lop(C3, mul(C3, x, mul(C3, w, y)))
assert torch.allclose(lhs, comm(Lop(C3, w), comm(Lop(C3, x), Lop(C3, y))))     # Lemma E.1(4)

def jpow(C, x, k):                                # Jordan powers x^1 = x, x^(k+1) = x o x^k
    p = x
    for _ in range(k - 1):
        p = mul(C, x, p)
    return p

def coords(Z):                                    # complex Hermitian 3 x 3 -> coordinates of Herm_3(C)
    return torch.cat([Z.diagonal().real] + [torch.stack([Z[i, j].real, Z[i, j].imag]) for i, j in ((0, 1), (0, 2), (1, 2))])

G = torch.randn(3, 3, dtype=torch.complex128)
Z = G + G.mH
lam = torch.linalg.eigvalsh(Z)                    # distinct and non-zero (almost surely)
xz = coords(Z)
powers = torch.stack([jpow(C3, xz, k) for k in range(1, 4)])  # x, x^2, x^3
cs = []
for i in range(3):                                # Lagrange idempotents p_i(x) of Lemma E.3
    others = [lam[j] for j in range(3) if j != i]
    coef = torch.tensor([others[0] * others[1], -(others[0] + others[1]), 1.0], dtype=dt)
    den = lam[i] * (lam[i] - others[0]) * (lam[i] - others[1])
    cs.append((coef[:, None] * powers).sum(0) / den)
for i in range(3):
    for j in range(3):
        assert torch.allclose(mul(C3, cs[i], cs[j]), cs[i] if i == j else torch.zeros(9, dtype=dt), atol=1e-10)
assert torch.allclose(sum(l * c for l, c in zip(lam, cs)), xz)

c = cs[0] + cs[1]                                 # the idempotent c_1 + c_2: Proposition E.4
Lc = Lop(C3, c)
assert torch.allclose(Lc @ (2 * Lc - torch.eye(9, dtype=dt)) @ (Lc - torch.eye(9, dtype=dt)), torch.zeros(9, 9, dtype=dt), atol=1e-10)
ev = torch.linalg.eigvals(Lc).real.sort().values
print("eigenvalues of L(c):", ev.numpy().round(6))
tau = lambda C: torch.einsum('abc,c->ab', C, torch.einsum('cbb->c', Lop(C, torch.eye(C.shape[0], dtype=dt))))  # Tr L(e_a o e_b)
T3 = tau(C3)
assert torch.linalg.eigvalsh(T3).min() > 0                                     # Theorem E.5(2)
assert torch.allclose(mul(C3, x, y) @ T3 @ w, x @ T3 @ mul(C3, y, w))         # Theorem E.5(1)
```

## E.4 ジョルダン枠とスペクトル定理

零でない冪等元は、零でない二つの直交する冪等元の和でないとき **原始的** であるという。**ジョルダン枠**
とは、互いに直交する原始冪等元の族 $(c_1,\dots,c_r)$ であって $\sum_ic_i=1$ を満たすもののことである。

**補題 E.6.** (1) $c$ が原始的ならば $J_1(c)=\mathbb Rc$ である。(2) 互いに直交する 零でない冪等元は
線形独立である；零でない冪等元はすべて互いに直交する原始冪等元の和である；そして $J$ は
ジョルダン枠をもつ。

*証明.* (1) $J_1(c)$ は $c$ を単位元とする形式的実な部分代数である（命題 E.4(3)）。$x\in J_1(c)$ に対し、補題 E.3 により
$c_i\in\mathbb R[x]\subseteq J_1(c)$ として $x=\sum_i\lambda_ic_i$ と書く。すると $c-c_i$ は $c_i$ に直交する冪等元である
（$(c-c_i)^2=c-2c_i+c_i$ かつ $c_i\circ(c-c_i)=0$）から、原始性により $c-c_i=0$ である。相異なる $c_i$ は直交するので、
それは高々一つしかなく、$x\in\mathbb Rc$ である。(2) $\sum_ia_if_i=0$ ならば、$f_j$ との積をとって $a_jf_j=0$ を得る。よって
冪等元 $e$ は、互いに直交する 零でない冪等元の可能な限り多くの個数の和として書ける。そのうちの一つ $f$ が、
$g,g'$ を直交する 零でない冪等元として $g+g'$ であったとすると、他のすべての成員 $f'$ に対し $g=f\circ g\in J_1(f)\subseteq J_0(f')$
（命題 E.4(4)）であり、$g'$ についても同様であるから、$f$ を $g,g'$ で置き換えることができる。$e=1$ ととればよい。$\square$

**定理 E.7（スペクトル定理）.** すべての $x\in J$ に対し、ジョルダン枠 $(c_1,\dots,c_r)$ と実数
$\lambda_1,\dots,\lambda_r$（相異なるとは限らない）が存在して $x=\sum_i\lambda_ic_i$ となる。

*証明.* 補題 E.3 により $x=\sum_{i=1}^k\lambda_ie_i$ と書き、$e_0=1-\sum_ie_i$、$\lambda_0=0$ を付け加える；$e_0$ は
すべての $e_i$ に直交する冪等元である。零でない各 $e_i$ を、互いに直交する原始冪等元に分解する
（補題 E.6(2)）。$e_i$ の一片 $f$ は $l\ne i$ に対し $J_1(e_i)\subseteq J_0(e_l)$ に属し、$e_l$ の一片は
$J_1(e_l)$ に属する；$J_0(e_l)\circ J_1(e_l)=0$ であるから、すべての片は互いに直交する。それらの和は $\sum_ie_i=1$ である。$\square$

$J$ のすべてのジョルダン枠が同じ個数の元をもつことは、系 E.28 で証明する。

**定理 E.8（枠に関するピアス分解）.** $(c_1,\dots,c_r)$ をジョルダン枠とし、$t_i=\tau(c_i,c_i)$、
$i\ne j$ に対し $J_{ij}=J_{ji}=J_{1/2}(c_i)\cap J_{1/2}(c_j)$ とおく。

1. $L(c_k)$ たちは可換であり、$J=\bigoplus_i\mathbb Rc_i\oplus\bigoplus_{i<j}J_{ij}$ であって、
   $x\in J_{ij}$ に対し $c_k\circ x=\frac12(\delta_{ki}+\delta_{kj})x$ である。
2. 相異なる $i,j,k$ に対し $J_{ij}\circ J_{jk}\subseteq J_{ik}$ である；また $\{i,j\}\cap\{k,l\}=\emptyset$ ならば $J_{ij}\circ J_{kl}=0$ である。
3. $J_{ij}$ 上の正定値二次形式 $N$ で $x^2=N(x)(c_i+c_j)$ となるものが存在する；その極形式を
   $\langle\cdot,\cdot\rangle$ とすると、$x,y\in J_{ij}$ に対し $x\circ y=\langle x,y\rangle(c_i+c_j)$ かつ $\tau(x,y)=(t_i+t_j)\langle x,y\rangle$ である。$J_{ij}\ne0$ ならば $t_i=t_j$ である。

*証明.* (1) $L(c_k)$ たちは命題 E.4(4) により可換であり、対角化可能で、和は $L(1)=\mathrm{id}$ である。共通の
固有ベクトルは $\sum_k\alpha_k=1$ を満たす固有値 $\alpha_k\in\{0,\frac12,1\}$ をもつ：一つの $\alpha_i=1$ で他が
消えるか（このときそのベクトルは $J_1(c_i)=\mathbb Rc_i$ に属する（補題 E.6(1)））、ある一組について $\alpha_i=\alpha_j=\frac12$ で
他が消えるかのいずれかである。よって $J_{ij}$ はちょうど第二の種類の共通固有空間であり、$J$ は共通
固有空間の直和である。(2) $x\in J_{ij}$、$y\in J_{jk}$ とする。$c_i$ に関して $x\in J_{1/2}$、$y\in J_0$ であるから
$x\circ y\in J_{1/2}(c_i)$ である；$c_k$ に関しても同様に $x\circ y\in J_{1/2}(c_k)$ である。したがって $x\circ y\in J_{ik}$ である。
代わりに $y\in J_{kl}$ ならば、$x\circ y$ は $p=i,j,k,l$ に対し $J_{1/2}(c_p)$ に属する；四つの固有値がすべて
$\frac12$ である共通固有ベクトルは存在しないので、$x\circ y=0$ である。(3) $x,y\in J_{ij}$ に対し、$p=i,j$ について $x\circ y\in J_0(c_p)+J_1(c_p)$ であり、
$k\notin\{i,j\}$ について $x\circ y\in J_0(c_k)$ である（命題 E.4(3)）；(1) により $x\circ y\in\mathbb Rc_i+\mathbb Rc_j$ である。
$x^2=\alpha c_i+\beta c_j$ と書く。$y=c_i$ としたジョルダン恒等式から $(x^2\circ c_i)\circ x=\frac12\alpha x$ と
$x^2\circ(c_i\circ x)=\frac12x^2\circ x=\frac14(\alpha+\beta)x$ が得られるので、$x\ne0$ ならば $\alpha=\beta$ である；$N(x)=\alpha$ とおく。
$\tau(c_i,1)=\operatorname{Tr}L(c_i)=t_i$ であるから $\tau(x,x)=\tau(x^2,1)=N(x)(t_i+t_j)$ である；よって $N$ は二次形式であり、
定理 E.5(2) により正定値であって、偏極化により $x\circ y$ と $\tau$ の公式が得られる。最後に
$N(x)t_i=\tau(x^2,c_i)=\tau(x,x\circ c_i)=\frac12\tau(x,x)$ であり、$j$ についても同じであるから、ある $x\ne0$ があれば $t_i=t_j$ である。$\square$

## E.5 単純代数と階数 2

**定理 E.9.** $J$ は単純イデアルの直和 $J^{(1)}\oplus\dots\oplus J^{(s)}$ であり、
$p\ne q$ に対し $J^{(p)}\circ J^{(q)}=0$ である。各 $J^{(p)}$ は形式的実ジョルダン代数である。

*証明.* $I$ がイデアルならば、$\tau$ に関するその直交補空間 $I^\perp$ もイデアルである：
$\tau(J\circ I^\perp,I)=\tau(I^\perp,J\circ I)=0$。したがって $J=I\oplus I^\perp$ かつ $I\circ I^\perp\subseteq I\cap I^\perp=0$ であるから、
$J$ は代数 $I$ と $I^\perp$ の直和であり、一方の直和成分のイデアルは $J$ のイデアルである。$I$ を
次元が最小の 零でないイデアルとすると $I$ は単純であり、形式的実代数 $I^\perp$ に次元に関する帰納法を
適用して結論を得る。$\square$

**補題 E.10.** $(c_i)$ をジョルダン枠、$i,j,k,l$ を相異なる添字とし、$u\in J_{ij}$、$w\in J_{jk}$、$z\in J_{kl}$ とする。このとき
(1) $u\circ(u\circ w)=\frac14N(u)w$ であり、(2) $(u\circ w)\circ z=u\circ(w\circ z)$ である。

*証明.* (1) $x=u$ とした補題 E.1(1) を $c_j$ に適用する。$L(u^2)=N(u)L(c_i+c_j)$、$c_j\circ w=\frac12w$、
$c_j\circ u=\frac12u$、$c_j\circ(u\circ w)=0$（$u\circ w\in J_{ik}$ だから）を用いると、
$[L(u^2),L(w)]c_j=\frac14N(u)w-\frac12N(u)w$ かつ $2[L(u\circ w),L(u)]c_j=(u\circ w)\circ u$ であり、それらの和は消える。(2)
$(x,w,y)$ を $(u,z,w)$ で置き換えた補題 E.1(3) を $c_j$ に適用する。ここで $u\circ z=0$（定理 E.8(2)）、$c_j\circ u=\frac12u$、
$c_j\circ z=0$ であるから、$0=\frac12(z\circ w)\circ u-\frac12z\circ(w\circ u)$ である。$\square$

**命題 E.11.** $(c_1,\dots,c_r)$ を $J$ のジョルダン枠とし、$J_{ij}\ne0$ のとき $i\neq j$ を辺で結ぶ。

1. $J$ が単純であるための必要十分条件は、このグラフが連結であることである。
2. $J$ が単純で $r\ge2$ ならば、すべての $J_{ij}$（$i\ne j$）は同じ次元 $d\ge1$ をもち、すべての $t_i$ は
   $t=1+\frac{(r-1)d}2$ に等しく、$\dim J=r+d\,\frac{r(r-1)}2=rt$ である。

*証明.* (1) グラフが連結でなければ、$S$ を連結成分とし、$I$ を $i,j\in S$ に対する $c_i$ と $J_{ij}$ の
張る空間、$I'$ を補集合に対する同様の空間とする。すると $J=I\oplus I'$（$i\in S$、$j\notin S$ のすべての $J_{ij}$ は消える）であり、
定理 E.8 により $I\circ I\subseteq I$ かつ $I\circ I'=0$ である；よって $I$ は $\ne0,J$ なるイデアルである。逆に、グラフが
連結で $I\ne0$ がイデアルであるとする。$I$ は可換で対角化可能な $L(c_k)$ たちで不変であるから、共通固有空間との
共通部分の和である。$0\ne x\in I\cap J_{ij}$ ならば $c_i\circ x^2=N(x)c_i\in I$ である；よってある $c_i$ が
$I$ に属する。$c_i\in I$ かつ $0\ne x\in J_{ij}$ ならば、$x=2c_i\circ x\in I$ かつ $c_j=c_j\circ x^2/N(x)\in I$ である。
連結性によりすべての $c_j\in I$ となるので $1\in I$ であり、$I=J$ である。(2) $u\in J_{ij}$ を $N(u)=1$ を満たすものとし、$k\notin\{i,j\}$ とする。
定理 E.8(2) により $2L(u)$ は $J_{jk}$ を $J_{ik}$ に、$J_{ik}$ を $J_{jk}$ に写し、補題 E.10(1)（$i,j$ を
入れ替えたものも含めて）によりその平方は両方の上で恒等写像である。よって $J_{ij}\ne0$ ならば $\dim J_{ik}=\dim J_{jk}$ である。したがって辺 $ij$、$jk$ が
あれば辺 $ik$ がある；すると連結なグラフは完全グラフであり、添字を一つずつ変えていけば、すべての $\dim J_{ij}$
が等しいことがわかる。最後に、$L(c_i)$ は $\mathbb Rc_i$ 上で固有値 $1$、$r-1$ 個の空間 $J_{ij}$ 上で $\frac12$、それ以外で $0$
をもつので、$t_i=\operatorname{Tr}L(c_i)=1+\frac{(r-1)d}2$ である。$\square$

本節の残りでは、$J$ は単純、$(c_1,\dots,c_r)$ は $r\ge2$ のジョルダン枠であり、$x\in J_{ij}$、
$y\in J_{jk}$（$i,j,k$ は相異なる）に対し
$$
m(x,y)=2\,x\circ y\in J_{ik}.
$$
と書く。定理 E.8(3) と命題 E.11(2) により、すべての $J_{ij}$ 上で $\tau(x,y)=2t\langle x,y\rangle$ である。

**補題 E.12.** 相異なる $i,j,k$ と $x\in J_{ij}$、$y\in J_{jk}$、$z\in J_{ik}$ に対し：

1. $N(m(x,y))=N(x)N(y)$ かつ $m(x,m(x,y))=N(x)y$；
2. $\langle m(x,y),z\rangle=\langle y,m(x,z)\rangle=\langle x,m(y,z)\rangle$；
3. $l\notin\{i,j,k\}$ と $v\in J_{kl}$ に対し $m(m(x,y),v)=m(x,m(y,v))$。

*証明.* (2) 三つの数は $\tau(x\circ y,z)/t$、$\tau(y,x\circ z)/t$、$\tau(x,y\circ z)/t$ であり、
結合性と可換性により等しい。(1) 二番目の公式は補題 E.10(1) である。すると
$2tN(m(x,y))=4\tau(x\circ y,x\circ y)=4\tau(y,x\circ(x\circ y))=N(x)\tau(y,y)=2tN(x)N(y)$ である。(3) 補題 E.10(2)。$\square$

**定理 E.13（階数 1 と 2）.** $J$ を $r$ 個の元からなるジョルダン枠をもつ単純代数とする。$r=1$ ならば
$J=\mathbb Rc_1\cong\mathbb R$ である。$r=2$ ならば、$d=\dim J_{12}\ge1$ として $J\cong V_{d+1}$ である（スピン因子は 13.2節のとおり）。

*証明.* $r=1$：定理 E.8(1)。$r=2$：$W=\mathbb Rf\oplus J_{12}$、$f=c_1-c_2$ として $J=\mathbb R1\oplus W$ である。
定理 E.8 により、$x,y\in J_{12}$ に対し $f\circ f=1$、$f\circ x=\frac12x-\frac12x=0$、$x\circ y=\langle x,y\rangle1$ である。よって
$w,w'\in W$ に対し $w\circ w'=B(w,w')1$ であり、ここで $B(\alpha f+x,\beta f+y)=\alpha\beta+\langle x,y\rangle$ は内積である。
$U:W\to\mathbb R^{d+1}$ を等長写像とすると、$a1+w\mapsto(a,Uw)$ は
$(a1+w)\circ(b1+w')=(ab+B(w,w'))1+aw'+bw$ を $V_{d+1}$ の積に写す。$\square$

## E.6 合成代数：フルヴィッツの定理とフロベニウスの定理

**合成代数** とは、単位元 $1$ と正定値二次形式 $N$ をもつ代数 $A$ であって
$N(xy)=N(x)N(y)$ を満たすもののことである。極形式を $\langle\cdot,\cdot\rangle$ と書き（$N(x)=\langle x,x\rangle$）、
$\operatorname{Re}x=\langle x,1\rangle$、$\bar x=2\operatorname{Re}(x)1-x$ とおく。$N(1)=1$ であるから（補題 E.14(1)。その
証明は共役を用いない）、これは直線 $\mathbb R1$ に関する鏡映であり、したがって $x\mapsto\bar x$ は等長な
対合で、$\operatorname{Re}\bar x=\operatorname{Re}x$ である。合成代数の同型とは、等長な代数
同型のことである。**結合子** とは $(x,y,z)=(xy)z-x(yz)$ のことである。

**補題 E.14.** 合成代数において、すべての $x,y,z,w$ に対し：

1. $N(1)=1$、$\langle xy,xz\rangle=N(x)\langle y,z\rangle$、$\langle yx,zx\rangle=\langle y,z\rangle N(x)$ および
   $\langle xy,wz\rangle+\langle wy,xz\rangle=2\langle x,w\rangle\langle y,z\rangle$；
2. $\langle xy,z\rangle=\langle y,\bar xz\rangle=\langle x,z\bar y\rangle$；
3. $\bar x(xy)=N(x)y=(yx)\bar x$；したがって $\bar xx=x\bar x=N(x)1$、$x^2-2\operatorname{Re}(x)x+N(x)1=0$、
   $xy+yx=2\operatorname{Re}(x)y+2\operatorname{Re}(y)x-2\langle x,y\rangle1$、
   $\bar x(wy)+\bar w(xy)=2\langle x,w\rangle y$ および $(yx)\bar w+(yw)\bar x=2\langle x,w\rangle y$；
4. $\overline{xy}=\bar y\bar x$、$\operatorname{Re}(xy)=\operatorname{Re}(yx)$ および
   $\operatorname{Re}((xy)z)=\operatorname{Re}(x(yz))$；
5. $x(xy)=x^2y$ かつ $(yx)x=yx^2$；結合子は交代的である（二つの引数を入れ替えると符号が変わる）。

*証明.* (1) $N(1)=N(1)^2>0$。$N(xy)=N(x)N(y)$ を $y$ について、それぞれ $x$ について偏極化し、次に $\langle xy,xz\rangle=N(x)\langle y,z\rangle$ を $x$ について偏極化する。
(2) (1) の最後の恒等式で $w=1$ とすると $\langle xy,z\rangle+\langle y,xz\rangle=2\operatorname{Re}(x)\langle y,z\rangle$、
すなわち $\langle xy,z\rangle=\langle y,\bar xz\rangle$ となる；$\langle yx,zx\rangle=\langle y,z\rangle N(x)$ を $x$ について偏極化し、
二つの変数の一方を $1$ とすると $\langle yx,z\rangle=\langle y,z\bar x\rangle$ を得る。(3) (2) と (1) により、
すべての $z$ に対し $\langle\bar x(xy),z\rangle=\langle xy,xz\rangle=N(x)\langle y,z\rangle$ であり、右側についても同様である。$y=1$ とおき、
$\bar x=2\operatorname{Re}(x)1-x$ を代入する；そして偏極化する。(4) (2) と共役の等長性により
$\langle\bar y\bar x,z\rangle=\langle\bar x,yz\rangle=\langle\bar x\bar z,y\rangle=\langle\bar z,xy\rangle=\langle z,\overline{xy}\rangle$ である。
次に $\operatorname{Re}(xy)=\langle y,\bar x\rangle=\langle\bar y,x\rangle=\operatorname{Re}(yx)$ であり、
$\operatorname{Re}((xy)z)=\langle xy,\bar z\rangle=\langle y,\bar x\bar z\rangle=\langle yz,\bar x\rangle=\operatorname{Re}(x(yz))$ である。
(5) (3) により $x(xy)=(2\operatorname{Re}(x)-\bar x)(xy)=(2\operatorname{Re}(x)x-N(x))y=x^2y$ であり、右側についても同様である。よって
$(x,x,y)=(y,x,x)=0$ である；偏極化すると、結合子は最初の二つまたは最後の二つの引数を入れ替えると符号を変え、
これらの互換はすべての置換を生成する。$\square$

**ケーリー–ディクソン構成.** 単位元、対合 $b\mapsto\bar b$、二次形式 $N$ をもつ代数 $B$ に対し、
$\operatorname{CD}(B)=B\times B$ を
$$
(a,b)(c,d)=(ac-\bar db,\;da+b\bar c),\qquad N(a,b)=N(a)+N(b),\qquad\overline{(a,b)}=(\bar a,-b).
$$
で定める。$\mathbb R$（$\bar x=x$、$N(x)=x^2$）から出発して $\mathbb C=\operatorname{CD}(\mathbb R)$、
$\mathbb H=\operatorname{CD}(\mathbb C)$、$\mathbb O=\operatorname{CD}(\mathbb H)$ を定義する。次元はそれぞれ $1,2,4,8$ であり、各倍化は
直前に定義した対合と形式を用いる（下の補題 E.15(1) により、それらはその合成代数の
共役とノルムである）。
$\mathbb C$ では $(a,b)(c,d)=(ac-bd,ad+bc)$ である：これは複素数 $a+bi$ である。$B$ は部分代数 $B\times0$ である。

**補題 E.15.** (1) $B$ が結合的な合成代数ならば、$\operatorname{CD}(B)$ は単位元 $(1,0)$ と
共役 $\overline{(a,b)}=(\bar a,-b)$ をもつ合成代数である。(2) $\mathbb C$ は可換かつ結合的である。
$M(a,b)=\begin{pmatrix}a&-b\\\bar b&\bar a\end{pmatrix}$ は $M(\bar q)=M(q)^*$ を満たす単射な代数準同型 $\mathbb H\to M_2(\mathbb C)$ である；
したがって $\mathbb H$ は結合的である；それは可換ではない。$\mathbb O$ は結合的でない。(3)
$\mathbb R,\mathbb C,\mathbb H,\mathbb O$ は合成代数であり、$\mathbb R,\mathbb C,\mathbb H$ は可除代数である：
すべての $x\ne0$ は両側逆元 $\bar x/N(x)$ をもつので、$xy=0$ ならば $x=0$ または $y=0$ である。

*証明.* (1) $N((a,b)(c,d))=N(a)N(c)+N(d)N(b)+N(d)N(a)+N(b)N(c)-2\langle ac,\bar db\rangle+2\langle da,b\bar c\rangle$ であり、
補題 E.14(2) により $\langle ac,\bar db\rangle=\langle d(ac),b\rangle$、$\langle da,b\bar c\rangle=\langle(da)c,b\rangle$ であって、これらは
$B$ が結合的なら一致する。$(1,0)$ は単位元であり、$2\langle(a,b),(1,0)\rangle(1,0)-(a,b)=(\bar a,-b)$ である。(2) 積を
展開し $\mathbb C$ の可換性を用いると $M(a,b)M(c,d)=M(ac-\bar db,da+b\bar c)$ である；$M$ は線形かつ単射であり、
$M(\bar a,-b)=M(a,b)^*$ である。$\mathbb H$ において $(i,0)(0,1)=(0,i)$ であるが $(0,1)(i,0)=(0,-i)$ である。$\mathbb O=\mathbb H\times\mathbb H$ において
$p=((i,0),0)$、$q=((0,1),0)$、$s=(0,1)$ とおくと、$(pq)s=(0,(0,i))$ かつ $p(qs)=(0,(0,-i))$ である。(3) (1)、(2) と帰納法による。
$x\ne0$ ならば補題 E.14(3) により $x(\bar x/N(x))=(\bar x/N(x))x=1$ である；結合代数では、このとき $xy=0$ から
$y=(\bar x/N(x))(xy)=0$ が得られる。$\square$

**補題 E.16（倍化）.** $A$ を合成代数、$B\subseteq A$ を $1\in B$ なる部分代数とし、$i\in A$ は
$i\perp B$ かつ $N(i)=1$ を満たすとする。$a,b,c,d\in B$ に対し $Bi\perp B$、$N(a+bi)=N(a)+N(b)$ であり、
$$
(a+bi)(c+di)=(ac-\bar db)+(da+b\bar c)\,i .
$$
である。よって $B\oplus Bi$ は次元 $2\dim B$ の部分代数であり、$h:B\to D$ が合成代数の同型ならば $a+bi\mapsto(h(a),h(b))$ は同型
$B\oplus Bi\cong\operatorname{CD}(D)$ である。さらに $B$ は
結合的である；$B\oplus Bi$ が結合的ならば $B$ は可換である；$B\oplus Bi$ が可換ならば $B=\mathbb R1$ である。

*証明.* $B$ は $b$ とともに $\bar b=2\operatorname{Re}(b)1-b$ を含む。補題 E.14 により $\langle bi,a\rangle=\langle i,\bar ba\rangle=0$、
$N(bi)=N(b)$、$\bar i=-i$ であり、$c\in B$ に対し $ci+ic=2\operatorname{Re}(c)i$、すなわち $ci=i\bar c$ である。
補題 E.14(3) の恒等式から順に次が得られる：$(x,w,y)=(\bar a,i,b)$ とすると $a(ib)=i(\bar ab)$；$ib=\bar bi$ かつ
$i(\bar ab)=\overline{\bar ab}\,i=(\bar ba)i$ であるから、$b$ を $\bar b$ で置き換えて $a(bi)=(ba)i$。$(y,x,w)=(a,i,\bar c)$ とすると
$(ai)c=(a\bar c)i$。$x=i$ とすると $i(iy)=-y$。$(x,w,y)=(ai,i,\bar d)$ とすると、$\overline{ai}=-ai$（$ai\perp1$ だから）かつ
$\langle ai,i\rangle=\operatorname{Re}a$（補題 E.14(1)）であるから
$$
-(ai)(i\bar d)-i((ai)\bar d)=2\operatorname{Re}(a)\bar d .
$$
となる。ここで $i\bar d=di$ かつ $i((ai)\bar d)=i((ad)i)=i(i\,\overline{ad})=-\bar d\bar a$ であるから、
$(ai)(di)=\bar d\bar a-2\operatorname{Re}(a)\bar d=-\bar da$ である。これら四つの規則で $(a+bi)(c+di)$ を展開すると公式が得られ、
$N$ と積は $h$ のもとで $\operatorname{CD}(D)$ のものと一致する。$N$ は
$B\oplus Bi\subseteq A$ 上で乗法的であるから、補題 E.15(1) の証明の計算により、すべての
$b\in B$ に対し $\langle d(ac)-(da)c,b\rangle=0$ となり、$B$ は結合的である。$B\oplus Bi$ が結合的ならば $(ab)i=a(bi)=(ba)i$ であるから $ab=ba$ である。可換ならば、
$ai=ia=\bar ai$ であるから $a=\bar a\in\mathbb R1$ である。$\square$

**定理 E.17（フルヴィッツ）.** すべての合成代数は、$\mathbb R,\mathbb C,\mathbb H,\mathbb O$ のうちちょうど一つと同型である。

*証明.* $B_0=\mathbb R1\cong\mathbb R$ とおく。$B_k\ne A$ である限り、$N(i_{k+1})=1$ なる $i_{k+1}\perp B_k$ を選び、
$B_{k+1}=B_k\oplus B_ki_{k+1}$ とおく；補題 E.16 により $B_k\cong\operatorname{CD}^k(\mathbb R)$ の次元は $2^k$ である。$B_4$ が
存在したとすると、補題 E.16 により $B_3$ は結合的、したがって $B_2$ は可換、したがって $B_1=\mathbb R1$ となるが、これは次元が
$1\ne2$ である。よってある $k\le3$ に対し $A=B_k$ である。四つの代数は次元が異なる。$\square$

**補題 E.18.** $A$ を合成代数とする。(1) $A$ の任意の二元は、ある結合的な部分代数に含まれる。
(2)（**ムーファン恒等式**）すべての $u,x,y$ に対し $(ux)(yu)=u((xy)u)=(u(xy))u$ である。

*証明.* (1) $p,q\in A$ とする；$p\notin\mathbb R1$ と仮定してよい（そうでなければ両者を入れ替え、両方とも実なら $\mathbb R1$ を用いる）。
$B=\operatorname{span}\{1,p\}$ は部分代数（$p^2=2\operatorname{Re}(p)p-N(p)1$）であり、次元 $2$ の合成代数であるから、
$B\cong\mathbb C$ である（定理 E.17）。$q\in B$ ならば証明は終わる；そうでなければ $i$ を $q$ の $B$ に直交する成分を
正規化したものとする。補題 E.16 により $B\oplus Bi\cong\operatorname{CD}(\mathbb C)=\mathbb H$ は結合的で、$p,q$ を含む。
(2) 結合子は交代的である（補題 E.14(5)）から $(u,xy,u)=0$ であり、これが二番目の等式である。$X=(u,x,y)$ とおく。
補題 E.14(4) により $\operatorname{Re}X=0$ であり、$\langle u,X\rangle=\langle u\bar y,ux\rangle-N(u)\langle1,xy\rangle
=N(u)(\langle\bar y,x\rangle-\langle\bar x,y\rangle)=0$ である；よって補題 E.14(3) の偏極化された二次恒等式から
$uX+Xu=2\operatorname{Re}(u)X$、すなわち $\bar uX=Xu$ が得られる。定義を展開すると、すべての代数において
$$
(wx,y,z)-(w,xy,z)+(w,x,yz)=w(x,y,z)+(w,x,y)z .
$$
であることがわかる。$(w,x,y,z)=(u,u,x,y)$ とし、$(u^2,x,y)=2\operatorname{Re}(u)X$ を用いると：$(u,ux,y)=2\operatorname{Re}(u)X-uX=\bar uX=Xu$。
$(w,x,y,z)=(u,x,y,u)$ とし、$(ux,y,u)=(u,ux,y)=Xu$ と $(x,y,u)=X$（巡回置換）を用いると：$(u,x,yu)=uX+Xu-Xu=uX$。したがって
$(ux)(yu)-u((xy)u)=(u,x,yu)-u(x,y,u)=uX-uX=0$ である。$\square$

**定理 E.19（フロベニウス）.** 単位元 $1\ne0$ をもち零因子をもたない有限次元の結合的実代数 $A$ はすべて、
$\mathbb R$、$\mathbb C$ または $\mathbb H$ と同型である。

*証明.* (i) $x\in A$ に対し冪 $1,x,x^2,\dots$ は線形従属である；$\mu$ を $\mu(x)=0$ を満たす次数最小のモニックな実多項式とする。
$\deg f,\deg g\ge1$ として $\mu=fg$ ならば、$f(x)g(x)=0$ から $f(x)=0$ または $g(x)=0$ となる；よって $\mu$ は
既約であり、したがって（代数学の基本定理により）$\mu=t-\alpha$ または $\mu=(t-\alpha)^2+\beta^2$、$\beta>0$ である。
$V=\{v\in A:v^2\in(-\infty,0]\cdot1\}$ とおく。すると $A=\mathbb R1+V$ であり（第二の場合には $x-\alpha\in V$）、
$\mathbb R1\cap V=0$ である。(ii) $V$ は部分空間である：スカラー倍で閉じており、$u,v\in V$ に対し $u+v\in V$ を示す。
$u,v$ が線形従属ならば $u+v$ はそのどちらかの倍数である；そこで独立とする。すると $1,u,v$ は独立である：
たとえば $\alpha\ne0$ として $1=\alpha u+\gamma v$ ならば、$u^2=(1-2\gamma v+\gamma^2v^2)/\alpha^2\in\mathbb R1$ から $\gamma=0$
となり、$u\in\mathbb R1\cap V=0$ となる。(i) により、実数 $p_\pm,q_\pm$ を用いて $(u\pm v)^2=p_\pm(u\pm v)+q_\pm1$ と書ける；足し合わせると
$2u^2+2v^2=(p_++p_-)u+(p_+-p_-)v+(q_++q_-)1$ となるので $p_\pm=0$ である。$q_+>0$ ならば
$(u+v-\sqrt{q_+})(u+v+\sqrt{q_+})=0$ となって $u+v\in\mathbb R1$ となるが、これは不可能である；よって $u+v\in V$ である。(iii) $u,v\in V$ に対し
$uv+vu=(u+v)^2-u^2-v^2=-2\langle u,v\rangle1$ は $\langle u,u\rangle=-u^2\ge0$ を満たす対称双線形形式を定め、
$\langle u,u\rangle=0$ となるのは $u^2=0$、すなわち $u=0$ のときに限る：これは内積である。直交する $u,v\in V$ は反可換である。(iv) $V=0$ ならば
$A=\mathbb R$ である。$\dim V=1$ ならば、$i^2=-1$ として $V=\mathbb Ri$ であり、$A\cong\mathbb C$ である。$\dim V\ge2$ ならば、正規直交な
$i,j\in V$ をとり $k=ij$ とおく。すると $k^2=ijij=-i^2j^2=-1$ であるから $k\in V$ であり、$ik+ki=-j+j=0$、$jk+kj=i-i=0$ である：$i,j,k$ は
正規直交で、$ij=-ji=k$、$jk=-kj=i$、$ki=-ik=j$ である。$i,j,k$ に直交する単位ベクトル $l\in V$ があったとすると、
$lk=(li)j=-i(lj)=(ij)l=kl$ かつ $lk=-kl$ であるから $kl=0$ となり、不可能である。したがって $A=\operatorname{span}\{1,i,j,k\}$ である。$\mathbb H$ において、
元 $I=(i,0)$、$J=(0,1)$、$K=(0,i)$ は同じ乗積表をもつ（補題 E.15(2) により $IJ=K=-JI$ であり、
$I^2=J^2=-1$、$K=IJ$ により残りは結合性から決まる）から、$A\cong\mathbb H$ である。$\square$

以下のコードは、ケーリー–ディクソンの倍化により $\mathbb O$ を構成し、ランダムな八元数について補題 E.14 の恒等式、補題 E.18 のムーファン恒等式、
補題 E.15 の非結合性を検証し、さらに準同型 $M$ を検証する：

```python
o = torch.randn(4, 1000, 8, dtype=dt)             # random octonions
a, b, cc, u = o
N = lambda z: (z * z).sum(-1)
ip = lambda p, q: (p * q).sum(-1)
assert torch.allclose(N(cd(a, b)), N(a) * N(b))                               # composition
assert torch.allclose(ip(cd(a, b), cc), ip(b, cd(conj(a), cc)))               # Lemma E.14(2)
assert torch.allclose(ip(cd(a, b), cc), ip(a, cd(cc, conj(b))))
assert torch.allclose(conj(cd(a, b)), cd(conj(b), conj(a)))                   # Lemma E.14(4)
assert torch.allclose(cd(a, cd(a, b)), cd(cd(a, a), b))                       # alternative
assert torch.allclose(cd(cd(u, a), cd(b, u)), cd(u, cd(cd(a, b), u)))         # Moufang
assoc = cd(cd(a, b), cc) - cd(a, cd(b, cc))
assert assoc[..., 0].abs().max() < 1e-12 and assoc.abs().max() > 1          # Re-associative, not associative
e = torch.eye(8, dtype=dt)
print("(e1 e2) e4 =", cd(cd(e[1], e[2]), e[4]).tolist(), "\ne1 (e2 e4) =", cd(e[1], cd(e[2], e[4])).tolist())

def M(q):                                         # H -> M_2(C): (a, b) -> [[a, -b], [b*, a*]]
    za, zb = torch.complex(q[..., 0], q[..., 1]), torch.complex(q[..., 2], q[..., 3])
    return torch.stack([torch.stack([za, -zb], -1), torch.stack([zb.conj(), za.conj()], -1)], -2)
p4, q4 = torch.randn(2, 1000, 4, dtype=dt)
assert torch.allclose(M(cd(p4, q4)), M(p4) @ M(q4)) and torch.allclose(M(conj(p4)), M(p4).mH)
```

## E.7 座標化

$D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}$、$d=\dim D$、$n\ge1$ に対し、
$\operatorname{Herm}_n(D)=\{X\in M_n(D):X_{ji}=\overline{X_{ij}}\}$ に積 $X\circ Y=\frac12(XY+YX)$ を入れる。ここで
$(XY)_{ik}=\sum_jX_{ij}Y_{jk}$ である（これを定義するのに $D$ の結合性は必要ない）。対角成分は $\bar x=x$ を満たすので
実数である。$E_{ij}$ を行列単位とし、$i\ne j$ と $a\in D$ に対し $a[ij]=aE_{ij}+\bar aE_{ji}$ とおく；したがって
$a[ji]=\bar a[ij]$ である。

**補題 E.20.** $\operatorname{Herm}_n(D)$ は $\circ$ について閉じており、単位元 $\sum_iE_{ii}$ と次元
$n+d\frac{n(n-1)}2$ をもち、$a,b\in D$ と、最初の二式では任意の $i,j,k$、残りの式では相異なる $i,j,k,l$ に対し（$a[ij]$ では常に $i\ne j$）
$$
E_{ii}\circ E_{jj}=\delta_{ij}E_{ii},\quad E_{kk}\circ a[ij]=\tfrac12(\delta_{ki}+\delta_{kj})a[ij],\quad
a[ij]\circ b[ij]=\langle a,b\rangle(E_{ii}+E_{jj}),\quad a[ij]\circ b[jk]=\tfrac12(ab)[ik],\quad a[ij]\circ b[kl]=0 .
$$
が成り立つ。さらに $(X|Y)=\operatorname{Re}\operatorname{tr}(XY)$ は $(X\circ Y|Z)=(Y|X\circ Z)$ を満たす内積である。

*証明.* 行列の積により $a[ij]\,b[ij]=a\bar bE_{ii}+\bar abE_{jj}$ であり、$a\bar b+b\bar a=\bar ab+\bar ba=2\langle a,b\rangle$ である
（$x\bar x=\bar xx=N(x)$ を偏極化する）；$a[ij]\,b[jk]=abE_{ik}$、$b[jk]\,a[ij]=\bar b\bar aE_{ki}=\overline{ab}E_{ki}$ であり、
その他の積は直ちにわかる。基底の積はすべて $\operatorname{Herm}_n(D)$ に属するので、すべての積がそうである。次に、
$X\ne0$ に対し $(X|X)=\sum_{i,j}\operatorname{Re}(X_{ij}\overline{X_{ij}})=\sum_{i,j}N(X_{ij})>0$ である。
$T(X,Y,Z)=\operatorname{Re}\operatorname{tr}((XY)Z)=\sum_{i,j,k}\operatorname{Re}((X_{ij}Y_{jk})Z_{ki})$ とおく。補題 E.14(4) により
$\operatorname{Re}\operatorname{tr}(AB)=\operatorname{Re}\operatorname{tr}(BA)$、$T(X,Y,Z)=T(Z,X,Y)$ であり、共役をとって
（$\operatorname{Re}\bar p=\operatorname{Re}p$、$\overline{(pq)r}=\bar r(\bar q\bar p)$）$X,Y,Z$ がエルミートであることを用いると
$T(X,Y,Z)=\sum\operatorname{Re}(Z_{ik}(Y_{kj}X_{ji}))=T(Z,Y,X)$ である。よって $T$ は対称であり、
$(X\circ Y|Z)=\frac12(T(X,Y,Z)+T(Y,X,Z))=T(X,Y,Z)$ は $X,Y,Z$ について対称である。$\square$

**補題 E.21（向き）.** $P_{12},P_{13},P_{23}$ を内積空間、$P_{ji}=P_{ij}$ とし、$\{1,2,3\}$ の相異なる元からなる
各順序三つ組 $(i,j,k)$ に対し $m:P_{ij}\times P_{jk}\to P_{ik}$ を双線形写像とし、
$x\in P_{ij}$、$y\in P_{jk}$、$z\in P_{ik}$ に対し $m(x,y)=m(y,x)$（右辺は $(k,j,i)$ に対する写像）および $\langle m(x,y),z\rangle=\langle y,m(x,z)\rangle$ が成り立つとする。$D$ を合成代数とし、$\iota_{ij}:D\to P_{ij}$（$i<j$）を
$m(\iota_{12}a,\iota_{23}b)=\iota_{13}(ab)$ を満たす線形な等長全単射とする；$\iota_{ji}(a)=\iota_{ij}(\bar a)$ とおく。このとき
すべての順序三つ組とすべての $a,b\in D$ に対し $m(\iota_{ij}a,\iota_{jk}b)=\iota_{ik}(ab)$ である。

仮定は次の場合に成り立つ：(i) 単純な $J$ のジョルダン枠の三つの添字について、ピアス空間と $m=2\circ$ をとる場合
（可換性と補題 E.12(2)）；(ii) $\operatorname{Herm}_3(D)$ において、$P_{ij}=\{a[ij]\}$、
$\langle a[ij],b[ij]\rangle=\langle a,b\rangle$（共役は等長であるから矛盾なく定義される）、$m=2\circ$ をとる場合。実際
$m(a[ij],b[jk])=(ab)[ik]$ であり、$c[ik]$ に対し $\langle ab,c\rangle=\langle b,\bar ac\rangle$ である（補題 E.14(2)）。

*証明.* 主張が $(i,j,k)$ について成り立てば、$(k,j,i)$ についても成り立つ：
$m(\iota_{kj}a,\iota_{ji}b)=m(\iota_{ij}\bar b,\iota_{jk}\bar a)=\iota_{ik}(\bar b\bar a)=\iota_{ik}(\overline{ab})=\iota_{ki}(ab)$。
よって $(1,2,3)$ のもとで $(2,1,3)$ と $(1,3,2)$ を扱えば十分である。すべての $c\in D$ に対し、随伴性、$(1,2,3)$、
補題 E.14(2) を用いて
$$
\langle m(\iota_{21}a,\iota_{13}b),\iota_{23}c\rangle=\langle\iota_{13}b,m(\iota_{12}\bar a,\iota_{23}c)\rangle=\langle b,\bar ac\rangle=\langle ab,c\rangle=\langle\iota_{23}(ab),\iota_{23}c\rangle,
$$
$$
\langle m(\iota_{13}a,\iota_{32}b),\iota_{12}c\rangle=\langle\iota_{13}a,m(\iota_{12}c,\iota_{23}\bar b)\rangle=\langle a,c\bar b\rangle=\langle ab,c\rangle=\langle\iota_{12}(ab),\iota_{12}c\rangle,
$$
であり、$\iota_{23},\iota_{12}$ は全射である。（二行目では、$m(x,y)=m(y,x)$ の後で、$x$ の役割にある
$y=\iota_{32}b$ に対して随伴性を用いている。）$\square$

**定理 E.22（座標化）.** $J$ を、$r\ge3$ のジョルダン枠 $(c_1,\dots,c_r)$ をもつ単純代数とする。このとき
次元 $d$ の合成代数 $D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}$ と、$\Phi(E_{ii})=c_i$ を満たす代数同型
$\Phi:\operatorname{Herm}_r(D)\to J$ が存在する。$r\ge4$ ならば $D$ は結合的、すなわち $D\ne\mathbb O$ である。

*証明.* *ステップ1（一つの三角形）.* $N=1$ を満たす $e_{12}\in J_{12}$、$e_{23}\in J_{23}$ を選ぶ（命題 E.11）。
$D'=J_{13}$ 上に $p\cdot q=m(m(p,e_{23}),m(e_{12},q))$ を定義する。ここで積 $J_{13}\times J_{32}\to J_{12}$、
$J_{21}\times J_{13}\to J_{23}$、$J_{12}\times J_{23}\to J_{13}$ を用いる。補題 E.12(1) により $1'=m(e_{12},e_{23})$ は単位元である：
$p\cdot1'=m(m(p,e_{23}),e_{23})=p$ かつ $1'\cdot q=m(e_{12},m(e_{12},q))=q$；また $N(p\cdot q)=N(p)N(q)$ である。よって $(D',N)$ は
合成代数であり、定理 E.17 により $D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}$ と同型 $h:D\to D'$ が存在する。
$\iota_{13}=h$、$\iota_{12}=m(h(\cdot),e_{23})$、$\iota_{23}=m(e_{12},h(\cdot))$ とおく；これらは等長全単射であり
（補題 E.12(1)）、$m(\iota_{12}a,\iota_{23}b)=h(a)\cdot h(b)=\iota_{13}(ab)$ である。補題 E.21 により、$k=3$ に対し次が成り立つ：

> $(\mathrm T_k)$ 等長全単射 $\iota_{ij}:D\to J_{ij}$、$i<j\le k$ が存在して、$\iota_{ji}(a)=\iota_{ij}(\bar a)$ とおくとき、
> 相異なるすべての $i,j,l\le k$ とすべての $a,b\in D$ に対し $m(\iota_{ij}a,\iota_{jl}b)=\iota_{il}(ab)$ となる。

*ステップ2（添字の追加）.* $4\le k\le r$ として $(\mathrm T_{k-1})$ を仮定する。$N(v)=1$ を満たす $v\in J_{1k}$ を選び、
$2\le i<k$ に対し $\iota_{ik}(a)=m(\iota_{i1}a,v)$、また $\iota_{1k}(a)=m(\iota_{12}a,\iota_{2k}1)$ とおく。これらは
補題 E.12(1) により等長全単射である。補題 E.12(3)（毎回四つの相異なる添字について）と $(\mathrm T_{k-1})$ を用いると、相異なる
$i,j\in\{2,\dots,k-1\}$ に対し
$$
m(\iota_{ij}a,\iota_{jk}b)=m(\iota_{ij}a,m(\iota_{j1}b,v))=m(m(\iota_{ij}a,\iota_{j1}b),v)=m(\iota_{i1}(ab),v)=\iota_{ik}(ab).
$$
である。特に $\iota_{2k}(b)=m(\iota_{23}b,\iota_{3k}1)$ であり、したがって
$$
m(\iota_{12}a,\iota_{2k}b)=m(m(\iota_{12}a,\iota_{23}b),\iota_{3k}1)=m(\iota_{13}(ab),\iota_{3k}1)
=m(m(\iota_{12}(ab),\iota_{23}1),\iota_{3k}1)=m(\iota_{12}(ab),\iota_{2k}1)=\iota_{1k}(ab).
$$
である。$3\le j<k$ に対し $\iota_{1j}a=m(\iota_{12}a,\iota_{2j}1)$ であるから、
$m(\iota_{1j}a,\iota_{jk}b)=m(\iota_{12}a,m(\iota_{2j}1,\iota_{jk}b))=m(\iota_{12}a,\iota_{2k}b)=\iota_{1k}(ab)$ である。よってこの関係は
すべての三つ組 $i<j<k$ について成り立ち、各三角形 $\{i,j,k\}$ に補題 E.21 を適用すれば $(\mathrm T_k)$ が得られる。

*ステップ3（結合性）.* $r\ge4$ ならば、$(\mathrm T_4)$ と補題 E.12(3) により
$\iota_{14}((ab)c)=m(m(\iota_{12}a,\iota_{23}b),\iota_{34}c)=m(\iota_{12}a,m(\iota_{23}b,\iota_{34}c))=\iota_{14}(a(bc))$ である。

*ステップ4（同型）.* $(\mathrm T_r)$ を用いて線形写像
$\Phi(\sum_i\alpha_iE_{ii}+\sum_{i<j}a_{ij}[ij])=\sum_i\alpha_ic_i+\sum_{i<j}\iota_{ij}(a_{ij})$ を定義する；これは
定理 E.8(1) により全単射であり、すべての $i\ne j$ に対し $\Phi(a[ij])=\iota_{ij}(a)$ である。定理 E.8 と $(\mathrm T_r)$ により
$c_i\circ c_j=\delta_{ij}c_i$、$c_k\circ\iota_{ij}(a)=\frac12(\delta_{ki}+\delta_{kj})\iota_{ij}(a)$、
$\iota_{ij}(a)\circ\iota_{ij}(b)=\langle a,b\rangle(c_i+c_j)$、$\iota_{ij}(a)\circ\iota_{jk}(b)=\frac12\iota_{ik}(ab)$ および
$\iota_{ij}(a)\circ\iota_{kl}(b)=0$ である：これは補題 E.20 の乗積表である。双線形性により $\Phi(X\circ Y)=\Phi(X)\circ\Phi(Y)$ である。$\square$

## E.8 アルバート代数

**補題 E.23（$\operatorname{Herm}_3(D)$ の自己同型）.** $D$ を合成代数とし、$f_{12},f_{23},f_{13}$ を
すべての $a,b$ に対し $f_{13}(ab)=f_{12}(a)f_{23}(b)$ を満たす $D$ の等長な線形全単射とする。このとき
$$
\Phi\Big(\sum_i\alpha_iE_{ii}+\sum_{i<j}a_{ij}[ij]\Big)=\sum_i\alpha_iE_{ii}+\sum_{i<j}f_{ij}(a_{ij})[ij]
$$
は $(\operatorname{Herm}_3(D),\circ)$ の自己同型である。$N(u)=1$ を満たすすべての $u\in D$ に対し、
$L_ua=ua$、$R_ua=au$、$B_ua=uau$ として $(f_{12},f_{23},f_{13})=(L_u,R_u,B_u)$ はそのような三つ組である。

*証明.* 例 (ii) において $\iota_{ij}(a)=f_{ij}(a)[ij]$ に補題 E.21 を適用する：すると、すべての $i\ne j$ に対し $\Phi(a[ij])=\iota_{ij}(a)$ であり、
相異なるすべての $i,j,k$ に対し $2\iota_{ij}(a)\circ\iota_{jk}(b)=\iota_{ik}(ab)$ である。$\Phi$ は $E_{ii}$ を固定し、$f_{ij}$ は
等長であるから、$\Phi$ は補題 E.20 の乗積表のすべての積を保ち、したがって $\circ$ を保つ。三つ組については、$L_u,R_u,B_u$
は等長（$N(ua)=N(a)$）かつ全単射であり、$u(ab)u=(ua)(bu)$ はムーファン恒等式（補題 E.18(2)）である。$\square$

**定理 E.24.** (1) $D\in\{\mathbb R,\mathbb C,\mathbb H\}$ と $n\ge1$ に対し、$\operatorname{Herm}_n(D)$ はユークリッド的ジョルダン
代数であり、$\operatorname{Herm}_n(\mathbb H)$ は $\operatorname{Herm}_{2n}(\mathbb C)$ のジョルダン部分代数と同型である。
(2) **アルバート代数** $\operatorname{Herm}_3(\mathbb O)$（次元 $27$）と $\operatorname{Herm}_2(\mathbb O)$ は
ユークリッド的ジョルダン代数である。(3) これらのそれぞれにおいて、$(E_{11},\dots,E_{nn})$ はジョルダン枠であり、そのピアス空間
$J_{ij}=\{a[ij]:a\in D\}$ の次元は $d=\dim D$ である；その代数は単純である；$D=\mathbb R,\mathbb C,\mathbb H,\mathbb O$ に対し $\operatorname{Herm}_1(D)\cong\mathbb R$ かつ
$\operatorname{Herm}_2(D)\cong V_{d+1}$ である。

*証明.* (1) $M_n(D)$ は結合的実代数である（補題 E.15(2)）から、補題 E.20 と補題 13.1 により
$\operatorname{Herm}_n(D)$ はジョルダン代数であり、補題 E.20 によりユークリッド的である。各成分 $X_{ij}\in\mathbb H$ を
$2\times2$ ブロック $M(X_{ij})$ で置き換えることは単射な代数準同型 $M_n(\mathbb H)\to M_{2n}(\mathbb C)$ であり（補題 E.15(2)）、
$M(\bar q)=M(q)^*$ であるからエルミート行列をエルミート行列に写す；よって $\circ$ を保つ。

(2) $\operatorname{Herm}_3(\mathbb O)$ 上で $J(x,y)=(x^2\circ y)\circ x-x^2\circ(y\circ x)$ とおく；自己同型 $\Phi$ は
$J(\Phi x,\Phi y)=\Phi J(x,y)$ を満たす。$x=\sum_i\alpha_iE_{ii}+p[12]+q[23]+s[13]$ とする。補題 E.23 の形の自己同型を
三つ合成して、$(1,3)$ 成分と $(2,3)$ 成分を実にする。(a) $s\notin\mathbb R$ ならば、
$v\perp1$、$N(v)=1$ として $s=|s|(\cos\theta+v\sin\theta)$ と書く。すると $v^2=-1$ である（補題 E.14(3)）；
$u=\cos\frac\theta2-v\sin\frac\theta2$ に対し、可換な部分代数 $\operatorname{span}\{1,v\}\cong\mathbb C$ の中で計算すると
$B_us=u^2s=|s|$ である。(b) これで $s$ は実である。$N(u)=1$ かつ $u\perp q$ を満たす $u\perp1$ を選ぶ（$\dim\mathbb O=8$ であるから可能）；
すると $B_us=su^2=-s$ は実であり、$\bar q\in\operatorname{span}\{1,q\}$ であるから $\operatorname{Re}(qu)=\langle u,\bar q\rangle=0$ である。
(c) これで $q\perp1$ である；$q\ne0$ ならば $u=q/|q|$ ととる：$qu=q^2/|q|=-|q|$ と $B_us=-s$ は実である。よって、$\Phi x$ が
$(1,2)$ 成分の $p'$ を除いて実の成分をもつような自己同型 $\Phi$ が存在する。$J(\Phi x,\cdot)$ は線形であるから、
$y=E_{ii}$ と $y=c[ij]$、$c\in\mathbb O$ に対し $J(\Phi x,y)=0$ を示せば十分である。補題 E.18(1) により $p'$ と $c$ はある結合的
部分代数 $Q$ に含まれる；すると $\Phi x,y\in\operatorname{Herm}_3(Q)$ であり、これは
結合代数 $M_3(Q)$ の中にある $\operatorname{Herm}_3(\mathbb O)$ の部分代数であって、補題 13.1 により $J(\Phi x,y)=0$ である。したがってすべての $y$ に対し $J(x,\Phi^{-1}y)=\Phi^{-1}J(\Phi x,y)=0$ である：
ジョルダン恒等式が成り立ち、補題 E.20 によりこの代数はユークリッド的である。第三行と第三列が零である行列は、
単位元 $E_{11}+E_{22}$ をもつ部分代数（補題 E.20）をなし、これが $\operatorname{Herm}_2(\mathbb O)$ である。

(3) 補題 E.20 により $E_{ii}$ たちは和が $1$ の互いに直交する冪等元である；$L(E_{kk})$ は $E_{jj}$ に
$\delta_{kj}$ 倍として、$a[ij]$ に $\frac12(\delta_{ki}+\delta_{kj})$ 倍として作用するので、$J_1(E_{ii})=\mathbb RE_{ii}$ かつ $J_{1/2}(E_{ii})\cap
J_{1/2}(E_{jj})=\{a[ij]\}$ である。$E_{ii}$ は原始的である：$E_{ii}=g+g'$ が直交する 零でない冪等元の和ならば、
$g=E_{ii}\circ g\in\mathbb RE_{ii}$ であるから $g=E_{ii}$、$g'=0$ である。命題 E.11 のグラフは完全グラフであるから、この代数は
単純である。$n=1,2$ に対しては定理 E.13 を適用する。$\square$

**系 E.25.** $n\ge4$ に対し、$\operatorname{Herm}_n(\mathbb O)$ はジョルダン代数ではない。

*証明.* もしそうならば補題 E.10(2) が適用できる：その証明はジョルダン恒等式と補題 E.20 の積
（$c_k=E_{kk}$ として）しか用いない。すると、すべての $a,b,c\in\mathbb O$ に対し $\frac14((ab)c)[14]=(a[12]\circ b[23])\circ c[34]=a[12]\circ(b[23]\circ c[34])=\frac14(a(bc))[14]$
となり、補題 E.15(2) に矛盾する。$\square$

**定理 E.26（アルバート）.** $D$ を合成代数とする。結合的実代数 $A$ と、
$\varphi(x\circ y)=\frac12(\varphi(x)\varphi(y)+\varphi(y)\varphi(x))$ を満たす単射な線形写像 $\varphi:\operatorname{Herm}_3(D)\to A$ が存在するならば、
$D$ は結合的である。特に、アルバート代数は **例外的** である：それは $M_N(\mathbb C)$ のような、いかなる結合代数の
ジョルダン部分代数とも同型でない。

*証明.* $e_i=\varphi(E_{ii})$、$e=e_1+e_2+e_3=\varphi(1)$ とおく。すると $e_i^2=e_i$ であり、$i\ne j$ に対し $e_ie_j+e_je_i=0$ である；
左から、それぞれ右から $e_i$ を掛けると $e_ie_j+e_ie_je_i=0=e_ie_je_i+e_je_i$ となるので、
$e_ie_j=e_je_i=0$ である。$x=\varphi(X)$ に対し、$ex+xe=2x$ から $exe=ex=xe=x$ が得られるので、
$x_{pq}=e_pxe_q$ として $x=\sum_{p,q}x_{pq}$ である。$X=a[ij]$ に対し、関係 $e_i\circ x=e_j\circ x=\frac12x$ と $e_k\circ x=0$（$k\ne i,j$）は
$\frac12([p=i]+[q=i])x_{pq}=\frac12x_{pq}$ などを意味し、これにより $x_{ij}$ と $x_{ji}$ だけが残る。順序対に対し
$\alpha_{ij}(a)=e_i\varphi(a[ij])e_j$ と定義する。相異なる $i,j,k$ に対し、$\varphi(a[ij])\varphi(b[jk])=\alpha_{ij}(a)\alpha_{jk}(b)$ は
$\varphi(a[ij])\varphi(b[jk])+\varphi(b[jk])\varphi(a[ij])=\varphi((ab)[ik])$ の $e_iAe_k$ に属する唯一の項であるから、
$$
\alpha_{ij}(a)\,\alpha_{jk}(b)=\alpha_{ik}(ab).
$$
である。$\varphi(1[13])^2=\varphi(E_{11}+E_{33})$ から $\alpha_{13}(1)\alpha_{31}(1)=e_1$、$\alpha_{31}(1)\alpha_{13}(1)=e_3$ を得る；
$\varphi(a[13])^2=N(a)(e_1+e_3)$ から $\alpha_{13}(a)\alpha_{31}(\bar a)=N(a)e_1$ を得るので、$\alpha_{13}$ は単射である。ここで
$$
\alpha_{13}(ay)=\alpha_{12}(a)\alpha_{23}(y)=\alpha_{13}(a)\alpha_{32}(1)\alpha_{21}(1)\alpha_{13}(y)=\alpha_{13}(a)\alpha_{31}(1)\alpha_{13}(y),
$$
であるから、$\pi=\alpha_{31}(1)$ とおき $A$ における結合性を用いると、
$\alpha_{13}((ab)c)=\alpha_{13}(a)\pi\alpha_{13}(b)\pi\alpha_{13}(c)=\alpha_{13}(a(bc))$ である。よって $(ab)c=a(bc)$ である。補題 E.15(2) により、
$\mathbb O$ は結合的でない。$\square$

アルバートのもともとの議論（A. A. Albert, Ann. of Math. 35 (1934) 65–73）はこの種のものである。以下のコードは、ランダムな元について
$\operatorname{Herm}_3(\mathbb O)$ のジョルダン恒等式と、そのトレース形式の結合性および正値性を検証し、
$\operatorname{Herm}_4(\mathbb O)$ ではジョルダン恒等式が成り立たないことを数値的に示し（系 E.25）、
補題 E.23 の自己同型を検証し、さらにランダムな元に対して定理 E.24(2) の証明の三つのステップ (a)–(c) を実行する：

```python
n3, prod3 = herm(3, 8)
X, Y, Z = torch.randn(3, 200, 27, dtype=dt)
X2 = prod3(X, X)
assert (prod3(prod3(X2, Y), X) - prod3(X2, prod3(Y, X))).abs().max() < 1e-10    # Jordan identity, Theorem E.24(2)
n4, prod4 = herm(4, 8)
X4, Y4 = torch.randn(2, 200, n4, dtype=dt)
X42 = prod4(X4, X4)
defect4 = (prod4(prod4(X42, Y4), X4) - prod4(X42, prod4(Y4, X4))).abs().max().item()
print("Herm_4(O): max |Jordan defect| =", defect4)
assert defect4 > 1                                                               # Corollary E.25
retr = lambda v: v[..., :3].sum(-1)               # Re tr
assert torch.allclose(retr(prod3(prod3(X, Y), Z)), retr(prod3(X, prod3(Y, Z))))
sq = (X[:, :3] ** 2).sum(-1) + 2 * (X[:, 3:] ** 2).sum(-1)
assert torch.allclose(retr(prod3(X, X)), sq)                                    # (X|X) = sum of N(X_ij)

def moufang_auto(v, uu):                          # (L_u, R_u, B_u) on the entries a_12, a_23, a_13
    a12, a13, a23 = v[..., 3:11], v[..., 11:19], v[..., 19:27]   # pairs order (0,1),(0,2),(1,2)
    uu = uu.expand_as(a12)
    return torch.cat([v[..., :3], cd(uu, a12), cd(cd(uu, a13), uu), cd(a23, uu)], -1)
uu = torch.randn(8, dtype=dt); uu = uu / uu.norm()
assert torch.allclose(moufang_auto(prod3(X, Y), uu), prod3(moufang_auto(X, uu), moufang_auto(Y, uu)))

x = torch.randn(27, dtype=dt)                     # steps (a)-(c) of the proof of Theorem E.24(2)
s = x[11:19]
im = torch.cat([torch.zeros(1, dtype=dt), s[1:]])
v, th = im / im.norm(), torch.atan2(im.norm(), s[0]).item()
x1 = moufang_auto(x, math.cos(th / 2) * e[0] - math.sin(th / 2) * v)            # (a): (1,3) entry -> |s|
imq = torch.cat([torch.zeros(1, dtype=dt), x1[20:27]])
g = torch.cat([torch.zeros(1, dtype=dt), torch.randn(7, dtype=dt)])
g = g - (g @ imq) / (imq @ imq) * imq                                             # u perpendicular to 1 and q
x2 = moufang_auto(x1, g / g.norm())                                               # (b): (2,3) entry imaginary
x3 = moufang_auto(x2, x2[19:27] / x2[19:27].norm())                              # (c): (2,3) entry real
print("entries (1,3), (2,3) after (a)-(c):", x3[11:19].numpy().round(3), x3[19:27].numpy().round(3))
assert x3[12:19].abs().max() < 1e-12 and x3[20:27].abs().max() < 1e-12
assert abs(abs(x3[11].item()) - s.norm().item()) < 1e-12
```

## E.9 分類

**定理 E.27（Jordan–von Neumann–Wigner）.**

1. 有限次元の形式的実（同値に、ユークリッド的）ジョルダン代数はすべて単純イデアルの直和であり、
   単純なものはどれも、ある $n\ge1$ に対する $\operatorname{Herm}_n(\mathbb R)$、$\operatorname{Herm}_n(\mathbb C)$、
   $\operatorname{Herm}_n(\mathbb H)$、$m\ge2$ のスピン因子 $V_m$、またはアルバート代数
   $\operatorname{Herm}_3(\mathbb O)$ のいずれかと同型である。
2. 代数 $\mathbb R$；$V_m$（$m\ge2$）；$\operatorname{Herm}_n(D)$（$D=\mathbb R,\mathbb C,\mathbb H$、$n\ge3$）；
   $\operatorname{Herm}_3(\mathbb O)$ は単純かつ形式的実で、互いに同型でない。よって単純な形式的
   実ジョルダン代数はすべて、それらのうちちょうど一つと同型である。さらに $D=\mathbb R,\mathbb C,\mathbb H,\mathbb O$ に対し $\operatorname{Herm}_1(D)\cong\mathbb R$ かつ
   $\operatorname{Herm}_2(D)\cong V_{\dim D+1}$ である。
3. 単純な形式的実ジョルダン代数では、すべてのジョルダン枠は同じ個数 $r$ の元をもち、そのすべてのピアス
   空間 $J_{ij}$ は同じ次元 $d$ をもつ：$\mathbb R$ では $r=1$、$V_m$ では $(r,d)=(2,m-1)$、
   $\operatorname{Herm}_n(D)$、$n\ge2$ では $(n,\dim D)$ である。

*証明.* (1) 定理 E.9 により、$J$ は単純な形式的実イデアルの直和である。そのようなイデアルはある大きさ
$r$ のジョルダン枠をもつ（補題 E.6(2)）；それは $r=1$ ならば $\cong\mathbb R=\operatorname{Herm}_1(\mathbb R)$、$r=2$ ならば $\cong V_{d+1}$
（定理 E.13）であり、それ以外では $\cong\operatorname{Herm}_r(D)$ であって、$r\ge4$ ならば $D\ne\mathbb O$、$r=3$ ならば $D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}$
である（定理 E.22）。

(2) $\operatorname{Herm}_n(D)$ については定理 E.24 を参照。$V_m$ は補題 13.3 によりジョルダン代数であり、
$\langle(a,u),(b,v)\rangle=ab+u\cdot v$ に関してユークリッド的である。実際 $\langle x\circ y,z\rangle=abc+c\,u\cdot v+a\,v\cdot w+b\,u\cdot w$ は
$x=(a,u)$、$y=(b,v)$、$z=(c,w)$ について対称である。単位ベクトル $w$ に対し $c_{\pm w}=\frac12(1,\pm w)$ はジョルダン枠であり
（定理 13.4(1)；直交する 零でない二つの冪等元の和 $c_{w'}+c_{-w'}$ は $1$ であるから、$c_{\pm w}$ は原始的である）、
$J_{12}=\{(0,v):v\perp w\}$ の次元は $m-1\ge1$ である（次元の和は $m+1$）；よって命題 E.11 により $V_m$
は単純である。代数を区別するために、$s(J)$ を互いに直交する 零でない冪等元の最大個数とする
（補題 E.6(2) により有限）；$s$ と $\dim$ は同型で不変である。明らかに $s(\mathbb R)=1$ である。
$V_m$ では 零でない冪等元は $1$ と $c_w$ たちであり、$1\circ c_w=c_w$ であって、$c_w\circ c_{w'}=0$ となるのは $w'=-w$ のときに限る
（定理 13.4(1)）；よって $s(V_m)=2$ である。$D\in\{\mathbb R,\mathbb C,\mathbb H\}$ に対し $s(\operatorname{Herm}_n(D))=n$ である：
$E_{ii}$ たちにより $s\ge n$ である。逆に、直交する冪等元 $p_k\ne0$ は $k\ne l$ に対し $p_kp_l=0$ を満たす（定理 E.26 の
証明と同様）ので、$D^n$ の右 $D$-部分加群 $V_k=p_kD^n$ は零でなく独立である：$v_k\in V_k$ として $\sum_kv_k=0$ ならば、
$v_l=p_l\sum_kv_k=0$ である。零でない $V_k$ はある $v\ne0$ を含み、したがって実次元 $d$ の $vD$ を含む
（$D$ は零因子をもたない）；よって $sd\le nd$ である。$\operatorname{Herm}_3(\mathbb O)$ については、$E_{ii}$ たちにより $s\ge3$ である；そして $s$ 個の
互いに直交する 零でない冪等元は、$1$ からそれらの和を引いたものとあわせて、定理 E.7 の証明と同様に細分すると
$r'\ge s$ 個の元からなるジョルダン枠になる。命題 E.11(2) と定理 E.22 により、$d'\ge1$ として $27=r'+d'\frac{r'(r'-1)}2$ であり、
$r'\ge4$ ならば $d'\in\{1,2,4\}$ である。$r'=4,5,6$ に対し右辺はそれぞれ $\{10,16,28\}$、$\{15,25,45\}$、$\{21,36,66\}$ に属し、
$r'\ge7$ に対しては $27$ を超える；よって $s\le r'=3$ である。組 $(s,\dim)$ は $(1,1)$、$(2,m+1)$、$n\ge3$ に対する $(n,\frac{n(n+1)}2)$、
$(n,n^2)$、$(n,n(2n-1))$、そして $(3,27)$ であり、$\frac{n(n+1)}2<n^2<n(2n-1)$ かつ $6,9,15\ne27$ であるから、すべて異なる。
階数 $1$ と $2$ における同型は定理 E.24(3) である。

(3) $F$ を単純な $J$ の $r_F$ 個の元からなるジョルダン枠とすると、(1) の証明を $F$ で行うことにより、$J$ と
(2) の階数 $r_F$ のモデルとの同型が得られ、そのモデルの $s$ は $r_F$ である。よってすべての枠について $r_F=s(J)$ であり、
$r\ge2$ に対し $d=(\dim J-r)/\binom r2$ である（命題 E.11(2)）。$\square$

**系 E.28（階数とトレース）.** すべての形式的実ジョルダン代数において、すべてのジョルダン枠は同じ個数の
元をもつ。この数を **階数** という。$J$ が階数 $r$ の単純代数ならば、$\operatorname{tr}x=\frac r{\dim J}\operatorname{Tr}L(x)$ はすべての原始冪等元 $c$ に対し
$\operatorname{tr}c=1$ を満たし、したがって定理 E.7 のすべての分解に対し $\operatorname{tr}(\sum_i\lambda_ic_i)=\sum_i\lambda_i$ であって、
$\operatorname{tr}(x\circ y)$ は結合的な内積である。

*証明.* 定理 E.9 のように $J=\bigoplus_pJ^{(p)}$ とし、$c$ を原始的とする。その成分 $c^{(p)}\in J^{(p)}$ は
互いに直交する冪等元である（$c^2=\sum_p(c^{(p)})^2$）から、ちょうど一つが 零でない。したがって $J$ のジョルダン枠の
成員のうち $J^{(p)}$ に属するものは $J^{(p)}$ のジョルダン枠をなし（それらの和は $1$ の成分 $1^{(p)}$ である）、$J$ の階数は
$J^{(p)}$ たちの階数の和である（定理 E.27(3)）。単純な $J$ では、原始的な $c$ はあるジョルダン枠に属する
（定理 E.7 の証明と同様に $1-c$ を細分する）から、命題 E.11(2) により $\operatorname{Tr}L(c)=t=\dim J/r$ である（$r=1$ のときは
$J=\mathbb Rc$ で $\operatorname{Tr}L(c)=1$）。最後に $\operatorname{tr}(x\circ y)=\frac r{\dim J}\tau(x,y)$ である（定理 E.5）。$\square$

以下のコードは、各型について、ジョルダン枠 $(E_{ii})$ または $(c_w,c_{-w})$、$L(c_i)$ たちの共通固有空間としての
ピアス空間の次元（定理 E.8）、命題 E.11(2) の次元公式、そしてトレース
形式 $\tau$ の正値性（定理 E.5）を計算する：

```python
def peirce(C, frame):
    """dimensions of the joint Peirce spaces J_ij = {x : c_i o x = c_j o x = x/2}"""
    I = torch.eye(C.shape[0], dtype=dt)
    Ls = [Lop(C, f) for f in frame]
    return [C.shape[0] - torch.linalg.matrix_rank(torch.cat([Ls[i] - I / 2, Ls[j] - I / 2])).item()
            for i in range(len(frame)) for j in range(i + 1, len(frame))]

cases = [(f"Herm_{n}({'RCHO'[k]})", herm(n, 2 ** k), n, 2 ** k) for n, k in ((2, 0), (3, 0), (4, 0), (3, 1), (4, 1), (3, 2), (2, 3), (3, 3))]
cases += [(f"V_{m}", spin(m), 2, m - 1) for m in (2, 3, 5)]
for name, (dim, prod), r, d in cases:
    C = structure(dim, prod)
    if name.startswith('V'):
        frame = [0.5 * torch.eye(dim, dtype=dt)[0] + s * 0.5 * torch.eye(dim, dtype=dt)[1] for s in (1, -1)]
    else:
        frame = list(torch.eye(dim, dtype=dt)[:r])
    for i, f in enumerate(frame):
        assert torch.allclose(mul(C, f, f), f) and all(torch.allclose(mul(C, f, g), torch.zeros(dim, dtype=dt)) for g in frame[i + 1:])
    dims = peirce(C, frame)
    t = torch.einsum('cbb->c', Lop(C, torch.eye(dim, dtype=dt)))
    gram = torch.einsum('abc,c->ab', C, t)
    print(f"{name:10s} dim {dim:2d}  rank {r}  Peirce dimensions {dims}  min eig of trace form {torch.linalg.eigvalsh(gram).min().item():.2f}")
    assert dims == [d] * (r * (r - 1) // 2) and dim == r + d * r * (r - 1) // 2
    assert torch.linalg.eigvalsh(gram).min() > 0
```

## E.10 本付録が第13章のために証明すること

| 第13章の主張 | ここでの証明 |
|---|---|
| 形式的実 $\Leftrightarrow$ ユークリッド的（13.4節） | 定理 E.5 |
| 定理 13.7：単純イデアルへの分解、リスト、「ちょうど一つ」 | 定理 E.9、E.27 |
| $\mathbb K=\mathbb R,\mathbb C,\mathbb H,\mathbb O$ に対する $\operatorname{Herm}_2(\mathbb K)\cong V_{\dim\mathbb K+1}$ | 定理 E.24(3) |
| $\operatorname{Herm}_r(\mathbb H)$ は $\operatorname{Herm}_{2r}(\mathbb C)$ に埋め込まれる | 定理 E.24(1) |
| アルバート代数はジョルダン代数であり、例外的である | 定理 E.24(2)、E.26 |
| $\operatorname{Herm}_r(\mathbb O)$ がジョルダン代数であるのは $r\le3$ のときに限る | 定理 E.24(2)、系 E.25 |
| ジョルダン枠によるスペクトル定理、階数、$\operatorname{tr}z=\sum_i\mu_i$、$\operatorname{tr}(x\circ y)$ の結合性（13.6節） | 定理 E.7、系 E.28 |

**注記.** この分類は P. Jordan, J. von Neumann, E. Wigner, *On an algebraic generalization of the
quantum mechanical formalism*, Ann. of Math. 35 (1934) 29–64 によるものであり、$\operatorname{Herm}_3(\mathbb O)$ の例外性は
A. A. Albert（同じ巻、65–73）による。上の証明は Faraut–Korányi, *Analysis on Symmetric
Cones*（Oxford 1994）の Chapters II–V の構成に従う；フルヴィッツの定理（Nachr. Ges. Wiss. Göttingen 1898, 309–316）の倍化による証明は
J. H. Conway, D. A. Smith, *On Quaternions and Octonions*（2003）および T. A. Springer, F. D. Veldkamp, *Octonions, Jordan
Algebras and Exceptional Groups*（2000）にあるとおりであり、フロベニウスの定理（1878）の証明は R. S. Palais によるもの（Amer. Math. Monthly
75 (1968) 366–368）である。さらに読むべき文献（本講義のどこでも用いない）：多項式恒等式による例外性への別のアプローチとして C. M. Glennie, Pacific J. Math. 16 (1966)
47–59、そして歴史と一般論については McCrimmon の *A Taste of Jordan
Algebras*（2004）。
