[English](../E_jordan_classification.md) · [Español](../es/E_jordan_classification.md) · [Français](../fr/E_jordan_classification.md) · [Deutsch](../de/E_jordan_classification.md) · [中文](../zh/E_jordan_classification.md) · **日本語** · [Русский](../ru/E_jordan_classification.md)

# チュートリアル、付録 E. ユークリッド的ジョルダン代数の分類

このチュートリアルは講義ノートの[付録 E](../../lecture_notes/ja/E_jordan_classification.md)に対応する。付録 E は、形式的実ジョルダン代数の Jordan–von Neumann–Wigner による分類
（[定理 E.27](../../lecture_notes/ja/E_jordan_classification.md#r-E-27)）を証明する。コードはこれらの代数を素の PyTorch で、構造定数とケーリー–ディクソン積から構成し、
証明の主要なステップをランダムな元について検証する。ライブラリにおける作用素のジョルダン積は[チュートリアル 13](13_jordan_algebras.md)で扱う。

<a id="s-E-3"></a>

## E.3 トレース形式；ユークリッド的 = 形式的実

*講義ノート：* [E.3節](../../lecture_notes/ja/E_jordan_classification.md#s-E-3)；[補題 E.1](../../lecture_notes/ja/E_jordan_classification.md#r-E-1)、[補題 E.3](../../lecture_notes/ja/E_jordan_classification.md#r-E-3)、[命題 E.4](../../lecture_notes/ja/E_jordan_classification.md#r-E-4)、[定理 E.5](../../lecture_notes/ja/E_jordan_classification.md#r-E-5)。

形式的実ジョルダン代数では、トレース形式 $\tau(x,y)=\mathrm{Tr}\,L(x\circ y)$ は結合的かつ正定値である
（[講義ノート、定理 E.5](../../lecture_notes/ja/E_jordan_classification.md#r-E-5)）。以下のコードは、構造定数で与えた $\mathrm{Herm}_3(\mathbb C)$（[E.7節](../../lecture_notes/ja/E_jordan_classification.md#s-E-7)の代数 $\mathrm{Herm}_n(D)$ の一つ）において、
[補題 E.1(4)](../../lecture_notes/ja/E_jordan_classification.md#r-E-1)、[補題 E.3](../../lecture_notes/ja/E_jordan_classification.md#r-E-3)（冪等元 $p_i(x)$ はジョルダン冪だけで計算する）、[命題 E.4](../../lecture_notes/ja/E_jordan_classification.md#r-E-4)、定理 E.5 を検証する；
それを構成するのに用いるケーリー–ディクソン積 `cd` は[E.6節](../../lecture_notes/ja/E_jordan_classification.md#s-E-6)のものである。

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

<a id="s-E-6"></a>

## E.6 合成代数：フルヴィッツの定理とフロベニウスの定理

*講義ノート：* [E.6節](../../lecture_notes/ja/E_jordan_classification.md#s-E-6)；[補題 E.14](../../lecture_notes/ja/E_jordan_classification.md#r-E-14)、[補題 E.15](../../lecture_notes/ja/E_jordan_classification.md#r-E-15)、[補題 E.18](../../lecture_notes/ja/E_jordan_classification.md#r-E-18)。

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

<a id="s-E-8"></a>

## E.8 アルバート代数

*講義ノート：* [E.8節](../../lecture_notes/ja/E_jordan_classification.md#s-E-8)；[補題 E.20](../../lecture_notes/ja/E_jordan_classification.md#r-E-20)、[補題 E.23](../../lecture_notes/ja/E_jordan_classification.md#r-E-23)、[定理 E.24](../../lecture_notes/ja/E_jordan_classification.md#r-E-24)、[系 E.25](../../lecture_notes/ja/E_jordan_classification.md#r-E-25)。

以下のコードは、ランダムな元について、アルバート代数 $\mathrm{Herm}_3(\mathbb O)$ のジョルダン恒等式
（[講義ノート、定理 E.24(2)](../../lecture_notes/ja/E_jordan_classification.md#r-E-24)）と、そのトレース形式の結合性および正値性（[補題 E.20](../../lecture_notes/ja/E_jordan_classification.md#r-E-20)）を検証し、
$\mathrm{Herm}_4(\mathbb O)$ ではジョルダン恒等式が成り立たないことを数値的に示し（[系 E.25](../../lecture_notes/ja/E_jordan_classification.md#r-E-25)）、
[補題 E.23](../../lecture_notes/ja/E_jordan_classification.md#r-E-23) の自己同型を検証し、さらにランダムな元に対して定理 E.24(2) の証明の三つのステップ (a)–(c) を実行する：

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

<a id="s-E-9"></a>

## E.9 分類

*講義ノート：* [E.9節](../../lecture_notes/ja/E_jordan_classification.md#s-E-9)；[定理 E.5](../../lecture_notes/ja/E_jordan_classification.md#r-E-5)、[定理 E.8](../../lecture_notes/ja/E_jordan_classification.md#r-E-8)、[命題 E.11](../../lecture_notes/ja/E_jordan_classification.md#r-E-11)、[定理 E.27](../../lecture_notes/ja/E_jordan_classification.md#r-E-27)。

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

**次へ：** [チュートリアル、付録 F. 対称錐、順序同型とジョルダン写像](F_cones_and_order.md)。
