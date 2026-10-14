[English](../C_eigenvalue_densities.md) · [Español](../es/C_eigenvalue_densities.md) · [Français](../fr/C_eigenvalue_densities.md) · [Deutsch](../de/C_eigenvalue_densities.md) · **中文** · [日本語](../ja/C_eigenvalue_densities.md) · [Русский](../ru/C_eigenvalue_densities.md)

# 教程，附录 C. 特征值密度与 Page 公式

本教程配合讲义的[附录 C](../../lecture_notes/zh/C_eigenvalue_densities.md)，该附录由 Weyl 积分公式推导出高斯系综与 Wishart 系综的特征值密度，并证明了
Page 公式（[定理 C.17](../../lecture_notes/zh/C_eigenvalue_densities.md#r-C-17)）。我们用蒙特卡罗方法、数值积分和精确算术检验 Weyl 坐标卡的 Jacobi 行列式、GOE、
GUE 以及 Laguerre 系综的能级密度、Bartlett 分解、诱导测度和 Page 公式。

<a id="s-C-2"></a>

## C.2 Weyl 积分公式

*讲义：* [C.2 节](../../lecture_notes/zh/C_eigenvalue_densities.md#s-C-2)；[引理 C.4](../../lecture_notes/zh/C_eigenvalue_densities.md#r-C-4)、[定理 C.5](../../lecture_notes/zh/C_eigenvalue_densities.md#r-C-5)。

下面的检验对 $n=2,3$ 和两种 $\beta$ 数值地对 $\Psi$ 求导：在 $X=0$ 处 Jacobi 行列式是 $\Delta(\lambda)^\beta$，而
在固定的 $X\neq0$ 处它与 $\Delta(\lambda)^\beta$ 之比不依赖于 $\lambda$（引理 C.4(a)）。

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

## C.3 高斯系综

*讲义：* [C.3 节](../../lecture_notes/zh/C_eigenvalue_densities.md#s-C-3)；[定理 C.7](../../lecture_notes/zh/C_eigenvalue_densities.md#r-C-7)、[推论 C.8](../../lecture_notes/zh/C_eigenvalue_densities.md#r-C-8)。

我们对 $n=3$ 检验定理 C.7：GOE 的能级密度由联合密度在网格上对另外两个特征值积分算出，并与每个区间内
特征值的平均个数比较；我们还检验特征向量与特征值不相关（$|\tilde U_{11}|^2$ 不依赖于相位 $D$，所以可以从
`eigh` 返回的任何特征标架读出）。对服从 Haar 分布的 $U\in O(3)$，$\mathbb E|U_{11}|^2=\frac13$
（[讲义，C.3 节](../../lecture_notes/zh/C_eigenvalue_densities.md#s-C-3)）。

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

## C.4 正交多项式与单点函数

*讲义：* [C.4 节](../../lecture_notes/zh/C_eigenvalue_densities.md#s-C-4)；[引理 C.10](../../lecture_notes/zh/C_eigenvalue_densities.md#r-C-10)。

由[引理 C.10](../../lecture_notes/zh/C_eigenvalue_densities.md#r-C-10) 与讲义的推论 C.6，GUE 的能级密度是 $K_n(x)e^{-x^2/2}$，其中
$K_n=v^TM^{-1}v$ 由矩构成的矩阵 $M$ 算出。我们对 $n=3$ 用矩
$\int x^{2m}e^{-x^2/2}dx=\sqrt{2\pi}(2m-1)(2m-3)\cdots1$ 来检验这一点：

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

## C.5 Wishart 矩阵与诱导测度

*讲义：* [C.5 节](../../lecture_notes/zh/C_eigenvalue_densities.md#s-C-5)；[引理 C.11](../../lecture_notes/zh/C_eigenvalue_densities.md#r-C-11)（Bartlett），[引理 C.12](../../lecture_notes/zh/C_eigenvalue_densities.md#r-C-12)（Cholesky），[定理 C.14](../../lecture_notes/zh/C_eigenvalue_densities.md#r-C-14)（诱导测度）。

在本库中，`states.random_density_matrix(TypeIAlgebra([a], [a]), N, rank=b)` 抽取 $N$ 个态
$\rho_A=GG^*/\mathrm{Tr}\,GG^*$，其中 $G$ 是 $a\times b$ Ginibre 矩阵，也就是[定理 C.14](../../lecture_notes/zh/C_eigenvalue_densities.md#r-C-14) 的诱导测度；
`rho.eigenvalues()` 按升序返回每个块的特征值。

我们在 `torch.linalg.cholesky` 计算出的 Cholesky 因子上检验引理 C.11（由引理 C.12，它就是唯一的 $L\in T_+$），
并对 $a=2$ 检验定理 C.14，此时 $\mu_1\in(0,\frac12)$ 的密度 $\propto(1-2t)^2(t(1-t))^{b-2}$：

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

## C.6 Laguerre 多项式与 Page 公式

*讲义：* [C.6 节](../../lecture_notes/zh/C_eigenvalue_densities.md#s-C-6)；[引理 C.15](../../lecture_notes/zh/C_eigenvalue_densities.md#r-C-15)、[引理 C.16](../../lecture_notes/zh/C_eigenvalue_densities.md#r-C-16)、[定理 C.17](../../lecture_notes/zh/C_eigenvalue_densities.md#r-C-17)（Page 公式），[推论 C.18](../../lecture_notes/zh/C_eigenvalue_densities.md#r-C-18)。

我们检验：Laguerre 能级密度与直方图的比较；引理 C.16 的值 $A^{(\alpha)}_{k,m}$ 以及中间
结果 $\int x\log x\,\rho_{a,b}=ab\,I_0+\Sigma$（用数值积分，其中 $I_0$ 本身也用数值积分计算）；$\Sigma$ 的恒等式
（精确算术）；Page 公式，与 $a=2$ 时诱导密度的直接积分（精确有理数值）以及蒙特卡罗结果比较；
最后是推论 C.18：

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

**下一章：** [教程，附录 D. 量子动力学半群的生成元](D_gksl.md)。
