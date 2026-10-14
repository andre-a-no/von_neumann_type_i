[English](../C_eigenvalue_densities.md) · [Español](../es/C_eigenvalue_densities.md) · [Français](../fr/C_eigenvalue_densities.md) · [Deutsch](../de/C_eigenvalue_densities.md) · [中文](../zh/C_eigenvalue_densities.md) · **日本語** · [Русский](../ru/C_eigenvalue_densities.md)

# チュートリアル、付録 C. 固有値密度とペイジの公式

このチュートリアルは講義ノートの[付録 C](../../lecture_notes/ja/C_eigenvalue_densities.md)に対応する。付録 C は、ワイル積分公式からガウス・アンサンブルと
ウィシャート・アンサンブルの固有値密度を導き、ペイジの公式（[定理 C.17](../../lecture_notes/ja/C_eigenvalue_densities.md#r-C-17)）を証明する。ここでは、ワイル積分の
座標近傍のヤコビアン、GOE・GUE・ラゲール・アンサンブルの準位密度、バートレット分解、誘導測度、そしてペイジの公式を、
モンテカルロ法、数値積分、厳密な算術によって確かめる。

<a id="s-C-2"></a>

## C.2 ワイル積分公式

*講義ノート：* [C.2節](../../lecture_notes/ja/C_eigenvalue_densities.md#s-C-2)；[補題 C.4](../../lecture_notes/ja/C_eigenvalue_densities.md#r-C-4)、[定理 C.5](../../lecture_notes/ja/C_eigenvalue_densities.md#r-C-5)。

以下の検証では、$`n=2,3`$ と両方の $`\beta`$ について $`\Psi`$ を数値的に微分する：$`X=0`$ においてヤコビアンは $`\Delta(\lambda)^\beta`$ であり、
固定した $`X\neq0`$ においてはその $`\Delta(\lambda)^\beta`$ に対する比は $`\lambda`$ に依存しない（補題 C.4(a)）。

```python
import math
import torch
from fractions import Fraction
from torch_vn_algebra import TypeIAlgebra, states

torch.manual_seed(0)
dt, f64 = torch.complex128, torch.float64

def mc(samples, exact, name):
    """Monte Carlo mean of real samples with its standard error; must agree within 5 standard errors."""
    mean, err = samples.mean().item(), samples.std().item() / math.sqrt(samples.numel())
    print(f"{name}: {mean:.4f} +- {err:.4f}   (theory {exact:.4f})")
    assert abs(mean - exact) <= 5 * err + 1e-12
    return mean

def vandermonde(lam):
    """Delta(lam) = prod_{i<j} (lam_j - lam_i) along the last axis."""
    n = lam.shape[-1]
    return math.prod(lam[..., j] - lam[..., i] for i in range(n) for j in range(i + 1, n))

def coords(H, beta):
    """Coordinates of dH: h_ii, Re h_ij (i < j) and, for beta = 2, Im h_ij."""
    iu = torch.triu_indices(H.shape[-1], H.shape[-1], 1)
    off = H[..., iu[0], iu[1]]
    return torch.cat([H.diagonal(dim1=-2, dim2=-1).real, off.real] + ([off.imag] if beta == 2 else []), dim=-1)

def from_m(c, n, beta):
    """The element X of m (X* = -X, zero diagonal) with coordinates Re X_ij (, Im X_ij), i < j."""
    iu = torch.triu_indices(n, n, 1)
    p = iu.shape[1]
    X = torch.zeros(n, n, dtype=dt)
    X[iu[0], iu[1]] = c[:p] + 1j * c[p:] if beta == 2 else c[:p].to(dt)
    return X - X.mH

def Psi(p, n, beta):                                       # (X, lambda) -> e^X Lambda e^{-X}
    d = beta * n * (n - 1) // 2
    E = torch.linalg.matrix_exp(from_m(p[:d], n, beta))
    return coords(E @ torch.diag(p[d:]).to(dt) @ E.mH, beta)

def jacobian(p, n, beta, h=1e-6):
    J = torch.stack([(Psi(p + h * e, n, beta) - Psi(p - h * e, n, beta)) / (2 * h)
                     for e in torch.eye(len(p), dtype=f64)], dim=1)
    return torch.linalg.det(J).abs().item()

for n in (2, 3):
    for beta in (1, 2):
        d = beta * n * (n - 1) // 2
        X = 0.4 * torch.randn(d, dtype=f64)
        ratios = []
        for trial in range(3):
            lam = torch.randn(n, dtype=f64).sort().values
            D = vandermonde(lam).item() ** beta
            assert abs(jacobian(torch.cat([torch.zeros(d, dtype=f64), lam]), n, beta) / D - 1) < 1e-6
            ratios.append(jacobian(torch.cat([X, lam]), n, beta) / D)
        print(f"n = {n}, beta = {beta}: j(X) = {ratios[0]:.6f} for a random X")
        assert max(ratios) - min(ratios) < 1e-6              # j(X) does not depend on lambda
```

<a id="s-C-3"></a>

## C.3 ガウス・アンサンブル

*講義ノート：* [C.3節](../../lecture_notes/ja/C_eigenvalue_densities.md#s-C-3)；[定理 C.7](../../lecture_notes/ja/C_eigenvalue_densities.md#r-C-7)、[系 C.8](../../lecture_notes/ja/C_eigenvalue_densities.md#r-C-8)。

$`n=3`$ について定理 C.7 を検証する：GOE の準位密度を、同時密度を格子上で他の二つの固有値について積分して計算し、
各ビンの中の固有値の平均個数と比較する。また、固有ベクトルが固有値と無相関であることも検証する
（$`|\tilde U_{11}|^2`$ は位相 $`D`$ によらないので、`eigh` が返すどの固有フレームからも読み取れる）。ハール分布に従う
$`U\in O(3)`$ に対し $`\mathbb E|U_{11}|^2=\frac13`$ である（[講義ノート、C.3節](../../lecture_notes/ja/C_eigenvalue_densities.md#s-C-3)）。

```python
B, n = 200_000, 3
G = torch.randn(B, n, n, dtype=f64)
goe_lam, goe_vec = torch.linalg.eigh((G + G.mT) / 2)                    # Section 5.3 normalisation
edges = torch.arange(-4.0, 4.01, 1.0, dtype=f64)

h = 0.04                                                                # grid for the other two eigenvalues
xs = torch.arange(-6, 6 + h / 2, h, dtype=f64)
Y1, Y2 = torch.meshgrid(xs, xs, indexing='ij')
joint = lambda x1, x2, x3: (vandermonde(torch.stack([x1, x2, x3], -1)).abs()
                            * torch.exp(-(x1 ** 2 + x2 ** 2 + x3 ** 2) / 2))
marg = torch.stack([joint(x + 0 * Y1, Y1, Y2).sum() * h * h for x in xs])   # int |Delta| e^{-|x|^2/2} dy dz
rho_goe = 3 * marg / (marg.sum() * h)                                   # Corollary C.8, normalised by n! Z
for lo, hi in zip(edges[:-1], edges[1:]):
    sel = (xs >= lo - 1e-9) & (xs <= hi + 1e-9)
    count = ((goe_lam >= lo) & (goe_lam < hi)).sum(-1).double()
    mc(count, torch.trapezoid(rho_goe[sel], xs[sel]).item(), f"GOE(3): E #eigenvalues in [{lo:.0f}, {hi:.0f})")

u11 = goe_vec[:, 0, 0] ** 2                                             # |U_11|^2, Haar: mean 1/3
mc(u11, 1 / 3, "GOE(3): E |U_11|^2")
mc((u11 - 1 / 3) * goe_lam[:, 0], 0.0, "GOE(3): E (|U_11|^2 - 1/3) lambda_1")   # independence (Theorem C.7)
```

<a id="s-C-4"></a>

## C.4 直交多項式と1点関数

*講義ノート：* [C.4節](../../lecture_notes/ja/C_eigenvalue_densities.md#s-C-4)；[補題 C.10](../../lecture_notes/ja/C_eigenvalue_densities.md#r-C-10)。

[補題 C.10](../../lecture_notes/ja/C_eigenvalue_densities.md#r-C-10) と講義ノートの系 C.6 により、GUE の準位密度は $`K_n(x)e^{-x^2/2}`$ であり、
$`K_n=v^TM^{-1}v`$ はモーメントの行列 $`M`$ から計算される。$`n=3`$ について、モーメント
$`\int x^{2m}e^{-x^2/2}dx=\sqrt{2\pi}(2m-1)(2m-3)\cdots1`$ を用いてこれを検証する：

```python
def one_point(x, n, moment, w):
    """K_n(x) w(x) = v(x)^T M^{-1} v(x) w(x), M the matrix of moments (Lemma C.10)."""
    M = torch.tensor([[moment(i + j) for j in range(n)] for i in range(n)], dtype=f64)
    v = torch.stack([x ** i for i in range(n)], dim=-1)
    return (v @ torch.linalg.inv(M) * v).sum(-1) * w(x)

gauss_moment = lambda m: math.sqrt(2 * math.pi) * math.prod(range(m - 1, 0, -2)) if m % 2 == 0 else 0.0
Gc = torch.randn(B, n, n, dtype=dt)
gue_lam = torch.linalg.eigvalsh((Gc + Gc.mH) / math.sqrt(2))
for lo, hi in zip(edges[:-1], edges[1:]):
    x = torch.linspace(lo, hi, 401, dtype=f64)
    exact = torch.trapezoid(one_point(x, n, gauss_moment, lambda t: torch.exp(-t ** 2 / 2)), x).item()
    mc(((gue_lam >= lo) & (gue_lam < hi)).sum(-1).double(), exact, f"GUE(3): E #eigenvalues in [{lo:.0f}, {hi:.0f})")
```

<a id="s-C-5"></a>

## C.5 ウィシャート行列と誘導測度

*講義ノート：* [C.5節](../../lecture_notes/ja/C_eigenvalue_densities.md#s-C-5)；[補題 C.11](../../lecture_notes/ja/C_eigenvalue_densities.md#r-C-11)（バートレット）、[補題 C.12](../../lecture_notes/ja/C_eigenvalue_densities.md#r-C-12)（コレスキー）、[定理 C.14](../../lecture_notes/ja/C_eigenvalue_densities.md#r-C-14)（誘導測度）。

ライブラリでは、`states.random_density_matrix(TypeIAlgebra([a], [a]), N, rank=b)` が、$`G`$ を $`a\times b`$ の
ジニブル行列として $`N`$ 個の状態 $`\rho_A=GG^*/\mathrm{Tr}\,GG^*`$ をサンプリングする。これは[定理 C.14](../../lecture_notes/ja/C_eigenvalue_densities.md#r-C-14) の
誘導測度である。`rho.eigenvalues()` は各ブロックの固有値を昇順で返す。

補題 C.11 を `torch.linalg.cholesky` で計算したコレスキー因子（補題 C.12 により、これは一意的な $`L\in T_+`$ である）について検証し、
定理 C.14 を $`a=2`$ について検証する。このとき $`\mu_1\in(0,\frac12)`$ は密度 $`\propto(1-2t)^2(t(1-t))^{b-2}`$ をもつ：

```python
a, b, B = 3, 5, 100_000
G = torch.randn(B, a, b, dtype=dt)                        # E|G_ij|^2 = 1
W = G @ G.mH
L = torch.linalg.cholesky(W)                              # W = L L*, L in T_+
for i in range(a):                                        # 0-based i: L_ii^2 ~ Gamma(b - i) (Lemma C.11, index i + 1)
    mc(L[:, i, i].real ** 2, b - i, f"E L_{i + 1}{i + 1}^2")
    mc(L[:, i, i].real ** 4, (b - i) * (b - i + 1), f"E L_{i + 1}{i + 1}^4")
mc(L[:, 2, 0].abs() ** 2 * L[:, 1, 1].real ** 2, 1.0 * (b - 1), "E |L_31|^2 L_22^2")   # independence

for b2 in (2, 4):
    rho = states.random_density_matrix(TypeIAlgebra([2], [2], precision='double'), 100_000, rank=b2)
    mu1 = rho.eigenvalues()[0][:, 0]                      # smallest eigenvalue of rho_A
    t = torch.linspace(0, 0.5, 20_001, dtype=f64)
    dens = (1 - 2 * t) ** 2 * (t * (1 - t)) ** (b2 - 2)
    dens = dens / torch.trapezoid(dens, t)
    for lo, hi in ((0.0, 0.1), (0.1, 0.25), (0.25, 0.4)):
        sel = (t >= lo) & (t <= hi)
        mc(((mu1 >= lo) & (mu1 < hi)).double(), torch.trapezoid(dens[sel], t[sel]).item(),
           f"a = 2, b = {b2}: P({lo} <= mu_1 < {hi})")
```

<a id="s-C-6"></a>

## C.6 ラゲール多項式とペイジの公式

*講義ノート：* [C.6節](../../lecture_notes/ja/C_eigenvalue_densities.md#s-C-6)；[補題 C.15](../../lecture_notes/ja/C_eigenvalue_densities.md#r-C-15)、[補題 C.16](../../lecture_notes/ja/C_eigenvalue_densities.md#r-C-16)、[定理 C.17](../../lecture_notes/ja/C_eigenvalue_densities.md#r-C-17)（ペイジの公式）、[系 C.18](../../lecture_notes/ja/C_eigenvalue_densities.md#r-C-18)。

ラゲール準位密度をヒストグラムと比較し、補題 C.16 の値 $`A^{(\alpha)}_{k,m}`$ と中間
結果 $`\int x\log x\,\rho_{a,b}=ab\,I_0+\Sigma`$ を数値積分で（$`I_0`$ 自身も数値積分で計算して）検証し、
$`\Sigma`$ に関する恒等式を厳密な算術で検証し、ペイジの公式を $`a=2`$ に対する誘導密度の直接積分（厳密な有理数値）
およびモンテカルロと比較し、最後に系 C.18 を検証する：

```python
def laguerre(k, al, x):
    return sum((-1) ** j * math.comb(k + al, k - j) * x ** j / math.factorial(j) for j in range(k + 1))

def rho_laguerre(x, a, b):                                # Lemma C.15
    al = b - a
    return sum(math.factorial(k) / math.factorial(k + al) * laguerre(k, al, x) ** 2 for k in range(a)) * x ** al * torch.exp(-x)

lam = torch.linalg.eigvalsh(W)                            # a = 3, b = 5 from above
x = torch.linspace(0, 80, 160_001, dtype=f64)
r = rho_laguerre(x, a, b)
assert abs(torch.trapezoid(r, x).item() - a) < 1e-8 and abs(torch.trapezoid(x * r, x).item() - a * b) < 1e-8
for lo, hi in ((0, 1), (1, 2), (2, 4), (4, 8), (8, 12), (12, 20)):
    sel = (x >= lo) & (x <= hi)
    mc(((lam >= lo) & (lam < hi)).sum(-1).double(), torch.trapezoid(r[sel], x[sel]).item(),
       f"Laguerre(3, 5): E #eigenvalues in [{lo}, {hi})")

def Hn(m):
    return sum(Fraction(1, i) for i in range(1, m + 1))

def page(a, b):
    return Hn(a * b) - Hn(b) - Fraction(a - 1, 2 * b)

vv = torch.linspace(0, 12, 1_200_001, dtype=f64)[1:]      # x = v^2 smooths log x at 0
def log_integral(g):                                      # int_0^oo log(x) g(x) dx
    xv = vv ** 2
    return torch.trapezoid(xv.log() * g(xv) * 2 * vv, vv).item()
I0 = log_integral(lambda t: torch.exp(-t))
print(f"I_0 = {I0:.10f} (quadrature)")
for al in range(3):                                       # Lemma C.16
    for k in range(4):
        for m in range(k, 5):
            A = log_integral(lambda t: laguerre(k, al, t) * laguerre(m, al, t) * t ** al * torch.exp(-t))
            Im = math.factorial(k + al) * (I0 + float(Hn(k + al)))
            assert abs(A - (Im / math.factorial(k) if m == k else -math.factorial(k + al) / (math.factorial(k) * (m - k)))) < 1e-6
for a_, b_ in ((2, 2), (2, 3), (3, 3), (3, 5), (4, 6)):
    al = b_ - a_
    Sigma = sum((k + al + 1) * Hn(k + al + 1) + k * Hn(k + al) + 2 * k for k in range(a_))
    assert Sigma == a_ * b_ * Hn(b_) + Fraction(a_ * (a_ - 1), 2)            # proof of Theorem C.17, exactly
    integral = log_integral(lambda t: t * rho_laguerre(t, a_, b_))
    assert abs(integral - (a_ * b_ * I0 + float(Sigma))) < 1e-6

print("a = 2:", ", ".join(f"b = {b_}: {page(2, b_)}" for b_ in range(2, 7)))
assert page(2, 2) == Fraction(1, 3)
for b_ in range(2, 7):                                    # a = 2: integrate the density of Theorem C.14
    v = torch.linspace(0, math.sqrt(0.5), 200_001, dtype=f64)[1:]
    t = v ** 2                                            # mu_1 = v^2 smooths the logarithm at 0
    dens = (1 - 2 * t) ** 2 * (t * (1 - t)) ** (b_ - 2) * 2 * v
    S = -(t * t.log() + (1 - t) * (1 - t).log())
    assert abs((torch.trapezoid(dens * S, v) / torch.trapezoid(dens, v)).item() - float(page(2, b_))) < 1e-9

for a_, b_ in ((2, 2), (2, 3), (3, 3), (3, 5)):
    rho = states.random_density_matrix(TypeIAlgebra([a_], [a_], precision='double'), 40_000, rank=b_)
    mc(rho.von_neumann_entropy(), float(page(a_, b_)), f"E S(rho_A), a = {a_}, b = {b_}")   # Theorem C.17

for a_ in range(1, 12):
    for b_ in range(a_, 40):
        p = float(page(a_, b_))
        assert math.log(a_) - (a_ - 1) / (a_ * b_ + 1) - (a_ - 1) / (2 * b_) - 1e-12 <= p <= math.log(a_) - (a_ - 1) / (2 * b_) + 1e-12
        assert abs(p - math.log(a_) + a_ / (2 * b_)) <= 1 / (2 * b_) + 1e-12  # Corollary C.18
```

**次へ：** [チュートリアル、付録 D. 量子動力学的半群の生成子](D_gksl.md)。
