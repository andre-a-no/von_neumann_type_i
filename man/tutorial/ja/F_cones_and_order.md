[English](../F_cones_and_order.md) · [Español](../es/F_cones_and_order.md) · [Français](../fr/F_cones_and_order.md) · [Deutsch](../de/F_cones_and_order.md) · [中文](../zh/F_cones_and_order.md) · **日本語** · [Русский](../ru/F_cones_and_order.md)

# チュートリアル、付録 F. 対称錐、順序同型とジョルダン写像

このチュートリアルは講義ノートの[付録 F](../../lecture_notes/ja/F_cones_and_order.md)に対応する。付録 F は、ユークリッド的ジョルダン代数の平方の錐が対称錐であること、
ケッヒャー–ヴィンベルクの定理（[定理 F.12](../../lecture_notes/ja/F_cones_and_order.md#r-F-12)）、そしてカディソンの定理と、環 $\bigoplus_cM_{k_c}(\mathbb C)$ のジョルダン $*$-同型の形
（[定理 F.15](../../lecture_notes/ja/F_cones_and_order.md#r-F-15) と [F.18](../../lecture_notes/ja/F_cones_and_order.md#r-F-18)）を証明する。ここではこれらの結果を、行列、ライブラリのスピン因子 $V_4$、および
正定値行列の錐の上で数値的に確かめる。

<a id="s-F-2"></a>

## F.2 平方の錐は対称錐である

*講義ノート：* [F.2節](../../lecture_notes/ja/F_cones_and_order.md#s-F-2)；[定理 F.5](../../lecture_notes/ja/F_cones_and_order.md#r-F-5)、[定理 F.6](../../lecture_notes/ja/F_cones_and_order.md#r-F-6)、[系 F.7](../../lecture_notes/ja/F_cones_and_order.md#r-F-7)。

$V_4$ では、コードはライブラリの `SpinFactor`（[チュートリアル、13.2節](../../tutorial/ja/13_jordan_algebras.md#s-13-2)）を用いる：`random`（$Q$ の元には `cone=True` を指定）、
`inner`, `spectral_decomposition`, `in_cone`, `product`, `quadratic_representation`, `exp`, `log`, `apply_function`。

コードは $\mathrm{Herm}_4(\mathbb C)$ とローレンツ錐 $Q\subseteq V_4$ に対して定理 F.5 を（負のペアリングは
(1), $\langle y,c_i\rangle=\lambda_i\|c_i\|^2<0$ の証拠 $c_i$ である）、$V_4$ において定理 F.6 を（座標
$(a,u)$ において作用素 $L(b)$ と $P(y)$ は $5\times5$ 対称行列である）、そして
$\mathrm{Herm}_5(\mathbb R)$ において系 F.7 を検証する：

```python
import math
import torch
from torch_vn_algebra import SpinFactor

torch.manual_seed(0)
dt = torch.complex128
f64 = torch.float64
close = lambda A, B, tol=1e-10: (A - B).abs().max().item() < tol

G = torch.randn(2, 500, 4, 4, dtype=dt)
X, Y = G[0] @ G[0].mH, G[1] @ G[1].mH                             # squares in Herm_4(C)
assert ((X @ Y).diagonal(dim1=-2, dim2=-1).sum(-1).real > -1e-10).all()   # <Q, Q> >= 0
Z = torch.randn(500, 4, 4, dtype=dt)
Z = Z + Z.mH
w, U = torch.linalg.eigh(Z)
c1 = U[:, :, :1] @ U[:, :, :1].mH                                  # eigenprojection of the smallest eigenvalue
pair = (Z @ c1).diagonal(dim1=-2, dim2=-1).sum(-1).real            # <Z, c1> = lambda_min ||c1||^2
assert close(pair, w[:, 0]) and (pair < 0).float().mean() > 0.9   # witnesses of Z not in Q

V = SpinFactor(4)                                                  # <x, y> = tr(x o y) = V.inner
x, y = V.random(2000, cone=True), V.random(2000, cone=True)
assert (V.inner(x, y) > -1e-12).all()
z = V.random(2000)
lam, idem = V.spectral_decomposition(z)                            # (lambda-, lambda+), (c-, c+); ||c-||^2 = 1
assert close(V.inner(z, idem[:, 0]), lam[:, 0])
assert ((lam[:, 0] < 0) == ~V.in_cone(z, tol=0.0)).all()

I5 = torch.eye(5, dtype=f64)
Lmat = lambda b: V.product(b.expand(5, 5), I5).T                   # columns b o e_k
Pmat = lambda b: V.quadratic_representation(b.expand(5, 5), I5).T
for _ in range(20):
    b, yi, yg = V.random(), V.random(cone=True), V.random()
    E = torch.linalg.matrix_exp(Lmat(b))
    assert close(E @ V.unit(), V.exp(b))                           # Theorem F.6(1)
    assert V.in_cone((E @ x.T).T, tol=1e-9).all()
    assert close(Pmat(yi), torch.linalg.matrix_exp(2 * Lmat(V.log(yi))))      # Theorem F.6(3)
    lam_g, idem_g = V.spectral_decomposition(yg)
    s = (torch.where(lam_g >= 0, 1.0, -1.0)[:, None] * idem_g).sum(0)
    assert close(Pmat(yg), Pmat(s) @ Pmat(V.apply_function(yg, torch.abs)))  # P(y) = P(s) P(|y|)
    assert V.in_cone((Pmat(yg) @ x.T).T, tol=1e-9).all()                     # Theorem F.6(4)

xs = torch.randn(300, 5, 5, dtype=f64)
xs, ys = xs + xs.mT, torch.randn(300, 5, 5, dtype=f64)
ys = ys + ys.mT                                                    # Corollary F.7 in Herm_5(R), tau = Tr
absm = lambda A: (lambda w, U: U @ torch.diag_embed(w.abs()) @ U.mT)(*torch.linalg.eigh(A))
lhs = absm(ys @ xs @ ys).diagonal(dim1=-2, dim2=-1).sum(-1)
rhs = (ys @ absm(xs) @ ys).diagonal(dim1=-2, dim2=-1).sum(-1)
assert (rhs - lhs > -1e-9 * rhs).all()
```

<a id="s-F-3"></a>

## F.3 ケッヒャー–ヴィンベルクの定理

*講義ノート：* [F.3節](../../lecture_notes/ja/F_cones_and_order.md#s-F-3)；[命題 F.11](../../lecture_notes/ja/F_cones_and_order.md#r-F-11)、[定理 F.12](../../lecture_notes/ja/F_cones_and_order.md#r-F-12)、[例 F.13](../../lecture_notes/ja/F_cones_and_order.md#r-F-13)。

コードは $\mathrm{Herm}_3(\mathbb R)$ の正規直交座標において $f(x)=-\frac{n+1}2\log\det x$ を自動微分することで $e$ における
$D^3f$ を計算し、命題 F.11 と定理 F.12 の主張を検証する：

```python
from torch.autograd.functional import jacobian, hessian

n = 3
pairs = [(i, j) for i in range(n) for j in range(i, n)]
Nd = len(pairs)                                                    # N = 6
Bs = torch.zeros(Nd, n, n, dtype=f64)                              # orthonormal basis for Tr(xy)
for k, (i, j) in enumerate(pairs):
    Bs[k, i, j] = Bs[k, j, i] = 1.0 if i == j else 1 / math.sqrt(2)
mat = lambda v: torch.einsum('k,kij->ij', v, Bs)
vec = lambda A: torch.einsum('kij,ij->k', Bs, A)
f = lambda v: -(n + 1) / 2 * torch.logdet(mat(v))                  # log phi up to a constant (Example F.13)
kappa = math.sqrt((n + 1) / 2)
e = vec(kappa * torch.eye(n, dtype=f64))
assert close(-jacobian(f, e), e)                                   # e* = e
assert close(hessian(f, e), torch.eye(Nd, dtype=f64))             # H(e) = I (Proposition F.11(2))
T = jacobian(lambda v: hessian(f, v, create_graph=True), e)        # T[a, b, c] = D^3 f(e)[B_a, B_b, B_c]
prod = lambda u, v: torch.einsum('a,b,abc->c', u, v, -0.5 * T)     # Theorem F.12
u, v = torch.randn(2, Nd, dtype=f64)
assert close(prod(e, v), v)                                        # unit
assert close(prod(u, v), vec(mat(u) @ mat(v) + mat(v) @ mat(u)) / (2 * kappa))
u2 = prod(u, u)
assert close(prod(u, prod(u2, v)), prod(u2, prod(u, v)))           # Jordan identity
```

<a id="s-F-4"></a>

## F.4 順序同型：有限次元におけるカディソンの定理

*講義ノート：* [F.4節](../../lecture_notes/ja/F_cones_and_order.md#s-F-4)；[定理 F.15](../../lecture_notes/ja/F_cones_and_order.md#r-F-15)、[定理 F.18](../../lecture_notes/ja/F_cones_and_order.md#r-F-18)。

[講義ノート、定理 F.15](../../lecture_notes/ja/F_cones_and_order.md#r-F-15)により、$\Phi$ と $\Phi^{-1}$ がともに正値であるような単位的線形全単射 $\Phi$ はジョルダン $*$-同型である。
コードは、[定理 F.18](../../lecture_notes/ja/F_cones_and_order.md#r-F-18) の形をもつそのような $M_2\oplus M_2\oplus M_3$ の写像を構成する：それは二つの $M_2$ ブロックを入れ替え、
二つのブロックを転置し、ユニタリで共役をとる。そして、$\Phi$ が単位的で両方向に正値であること、ジョルダン積を保つこと、
乗法的でないことを確かめる：

```python
def blocks_rand(ks, psd=False):
    out = []
    for k in ks:
        A = torch.randn(k, k, dtype=dt)
        out.append(A @ A.mH if psd else A + A.mH)
    return out

ks = (2, 2, 3)                                                     # M = M_2 (+) M_2 (+) M_3
Us = [torch.linalg.qr(torch.randn(k, k, dtype=dt))[0] for k in ks]
def Phi(x):                                                        # swap the M_2 blocks, transpose two blocks
    x1, x2, x3 = x
    return [Us[0] @ x2.mT @ Us[0].mH, Us[1] @ x1 @ Us[1].mH, Us[2] @ x3.mT @ Us[2].mH]
def Phi_inv(y):
    y1, y2, y3 = y
    return [Us[1].mH @ y2 @ Us[1], (Us[0].mH @ y1 @ Us[0]).mT, (Us[2].mH @ y3 @ Us[2]).mT]
is_psd = lambda x: all(torch.linalg.eigvalsh(b)[0] > -1e-10 for b in x)
one = [torch.eye(k, dtype=dt) for k in ks]
assert all(close(a, b) for a, b in zip(Phi(one), one))            # unital
for _ in range(50):
    p = blocks_rand(ks, psd=True)
    assert is_psd(Phi(p)) and is_psd(Phi_inv(p))                   # positive in both directions
x, y = blocks_rand(ks), blocks_rand(ks)
assert all(close(a, b) for a, b in zip(Phi_inv(Phi(x)), x))
jordan = lambda a, b: [0.5 * (s @ t + t @ s) for s, t in zip(a, b)]
assert all(close(a, b) for a, b in zip(Phi(jordan(x, y)), jordan(Phi(x), Phi(y))))   # Theorem F.15
assert not close(Phi([s @ t for s, t in zip(x, y)])[0], (Phi(x)[0] @ Phi(y)[0]))      # not multiplicative
```

<a id="s-F-5"></a>

## F.5 $\bigoplus_cM_{k_c}$ のジョルダン $*$-同型

*講義ノート：* [F.5節](../../lecture_notes/ja/F_cones_and_order.md#s-F-5)；[定理 F.17](../../lecture_notes/ja/F_cones_and_order.md#r-F-17)。

コードは $M_3$ のジョルダン自己同型 $x\mapsto Ux^TU^*$ を構成し、証明の手順（[講義ノート、定理 F.17](../../lecture_notes/ja/F_cones_and_order.md#r-F-17)）によって
$U$ を（位相を除いて）復元する：

```python
k = 3
E = lambda i, j: torch.zeros(k, k, dtype=dt).index_put_((torch.tensor(i), torch.tensor(j)), torch.tensor(1.0, dtype=dt))
U0 = torch.linalg.qr(torch.randn(k, k, dtype=dt))[0]
Psi = lambda A: U0 @ A.mT @ U0.mH                                  # the unknown Jordan *-automorphism

W = torch.stack([torch.linalg.eigh(Psi(E(i, i)))[1][:, -1] for i in range(k)], dim=1)   # Psi(E_ii) = xi_i xi_i^*
Psi1 = lambda A: W.mH @ Psi(A) @ W
assert all(close(Psi1(E(i, i)), E(i, i)) for i in range(k))
transposed = Psi1(E(0, 1))[1, 0].abs().item() > 0.5               # epsilon = - : Psi1(E_01) is a multiple of E_10
Psi2 = (lambda A: Psi1(A.mT)) if transposed else Psi1
alpha = torch.stack([Psi2(E(i, 0))[i, 0] for i in range(k)])       # u_i = alpha_{i1}
for i in range(k):
    for j in range(k):
        assert close(Psi2(E(i, j)), alpha[i] * alpha[j].conj() * E(i, j))
Ufound = W @ torch.diag(alpha)
X = torch.randn(k, k, dtype=dt)
assert transposed and close(Psi(X), Ufound @ X.mT @ Ufound.mH)     # Theorem F.17, second form
phase = (U0.mH @ Ufound)[0, 0]
assert close(U0.mH @ Ufound, phase * torch.eye(k, dtype=dt)) and abs(abs(phase) - 1) < 1e-12
print("U recovered up to the phase", complex(phase))
```

<a id="x-F"></a>

## 演習問題

講義ノートの[演習問題 2](../../lecture_notes/ja/F_cones_and_order.md#x-F)：$0<t<1$ に対し、$M_n$ 上の $\Psi_t(x)=(1-t)x+\frac tn\mathrm{Tr}(x)1$ は、逆写像が正値でない
単位的な正値全単射であり、射影を射影に写さない。$n=3$、$t=\frac14$、階数1の射影 $p$ についてこれを数値的に確かめよ。

<details>
<summary>解答</summary>

```python
nn, t = 3, 0.25
p = E(0, 0)
Psi_t = lambda A: (1 - t) * A + t / nn * torch.trace(A) * torch.eye(nn, dtype=dt)
Psi_t_inv = lambda A: (A - t / nn * torch.trace(A) * torch.eye(nn, dtype=dt)) / (1 - t)
assert close(Psi_t_inv(Psi_t(X)), X)
assert abs(torch.linalg.eigvalsh(Psi_t_inv(p))[0].item() + t / (nn * (1 - t))) < 1e-12
assert not close(Psi_t(p) @ Psi_t(p), Psi_t(p))                   # not a projection
```

</details>
