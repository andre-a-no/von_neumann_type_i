[English](../A_analysis.md) · [Español](../es/A_analysis.md) · [Français](../fr/A_analysis.md) · [Deutsch](../de/A_analysis.md) · **中文** · [日本語](../ja/A_analysis.md) · [Русский](../ru/A_analysis.md)

# 教程，附录 A. 分析与概率中的若干事实

本教程配合讲义的[附录 A](../../lecture_notes/zh/A_analysis.md)，该附录证明了各章按编号引用的分析与概率中的事实：多项式的零点集、
特征值不等式、高斯向量及其方向、高斯矩以及均匀间距。这里的检验只使用 NumPy，不使用本库：每个代码块对随机向量
或随机矩阵抽样，并与精确的结论比较。第一个代码块中定义的辅助函数 `mc` 断言蒙特卡罗均值与精确值之差在 $`5`$ 个
标准误差以内。

<a id="s-A-2"></a>

## A.2 解析函数的零点集

*讲义：* [A.2 节](../../lecture_notes/zh/A_analysis.md#s-A-2)；[推论 A.7](../../lecture_notes/zh/A_analysis.md#r-A-7)。

由[推论 A.7](../../lecture_notes/zh/A_analysis.md#r-A-7)，不恒为零的多项式在具有密度的随机向量处取零的概率为 $`0`$。对 $`\det Z`$，其中 $`Z`$
是元素相互独立且服从 $`N(0,1)`$ 的 $`3\times3`$ 实矩阵，这给出当 $`\varepsilon\to0`$ 时
$`\mathbb P(|\det Z|<\varepsilon)\to0`$。

```python
import math
import numpy as np

rng = np.random.default_rng(0)


def mc(samples, exact, label, kind="exact"):
    """Monte Carlo mean with standard error; the reference value must lie within 5 standard errors."""
    samples = np.asarray(samples, dtype=float)
    mean, err = samples.mean(), samples.std() / math.sqrt(samples.size)
    print(f"{label}: {mean:.4f} +- {err:.4f}   ({kind} {exact:.4f})")
    assert abs(mean - exact) <= 5 * err + 1e-12, label


# Corollary A.7: det Z = 0 has probability 0 for a real Gaussian 3 x 3 matrix, so by continuity of the measure
# P(|det Z| < eps) -> 0 as eps -> 0. Numerical observation (not a theorem of this appendix): the fractions
# below shrink by about a factor 10 per factor 10 in eps.
d = np.abs(np.linalg.det(rng.standard_normal((200_000, 3, 3))))
fractions = [(d < eps).mean() for eps in (1e-1, 1e-2, 1e-3)]
print("P(|det Z| < eps) for eps = 0.1, 0.01, 0.001:", fractions)
assert fractions[1] < 0.2 * fractions[0] and fractions[2] < 0.2 * fractions[1]
```

<a id="s-A-3"></a>

## A.3 Hermite 矩阵的特征值

*讲义：* [A.3 节](../../lecture_notes/zh/A_analysis.md#s-A-3)；[定理 A.8](../../lecture_notes/zh/A_analysis.md#r-A-8)、[推论 A.9](../../lecture_notes/zh/A_analysis.md#r-A-9)、[定理 A.10](../../lecture_notes/zh/A_analysis.md#r-A-10)、[引理 A.12](../../lecture_notes/zh/A_analysis.md#r-A-12)。

该代码块在随机 Hermite 矩阵上检验 Weyl 不等式
$`\lambda_j(A)+\lambda_1(E)\le\lambda_j(A+E)\le\lambda_j(A)+\lambda_n(E)`$（[定理 A.10](../../lecture_notes/zh/A_analysis.md#r-A-10)）、对 $`j`$ 维 $`U`$ 的
Courant–Fischer 界 $`\max_{x\in U,|x|=1}\langle x,Ax\rangle\ge\lambda_j(A)`$（[定理 A.8](../../lecture_notes/zh/A_analysis.md#r-A-8)）、交错不等式
$`\lambda_j(A)\le\lambda_j(Q^*AQ)\le\lambda_{j+n-m}(A)`$（[推论 A.9](../../lecture_notes/zh/A_analysis.md#r-A-9)），以及谱投影的 Lagrange 公式
（[引理 A.12](../../lecture_notes/zh/A_analysis.md#r-A-12)）。

```python
def herm(n, batch):
    G = rng.standard_normal((batch, n, n)) + 1j * rng.standard_normal((batch, n, n))
    return (G + G.conj().transpose(0, 2, 1)) / 2


n = 6
A, E = herm(n, 2000), 0.1 * herm(n, 2000)
lam_A, lam_AE, lam_E = np.linalg.eigvalsh(A), np.linalg.eigvalsh(A + E), np.linalg.eigvalsh(E)   # ascending
normE = np.linalg.norm(E, ord=2, axis=(1, 2))                                    # operator norm
assert (np.abs(lam_AE - lam_A).max(axis=1) <= normE + 1e-12).all()               # Theorem A.10
assert (lam_AE >= lam_A + lam_E[:, :1] - 1e-12).all() and (lam_AE <= lam_A + lam_E[:, -1:] + 1e-12).all()

a0, j = A[0], 3                          # Theorem A.8: every 3-dim subspace has max <x, Ax> >= lambda_3
for _ in range(200):
    Q, _ = np.linalg.qr(rng.standard_normal((n, j)) + 1j * rng.standard_normal((n, j)))
    assert np.linalg.eigvalsh(Q.conj().T @ a0 @ Q)[-1] >= lam_A[0, j - 1] - 1e-12
m = 4                                    # Corollary A.9: interlacing for a compression to 4 dimensions
Q, _ = np.linalg.qr(rng.standard_normal((n, m)) + 1j * rng.standard_normal((n, m)))
theta = np.linalg.eigvalsh(Q.conj().T @ a0 @ Q)
assert (theta >= lam_A[0, :m] - 1e-12).all() and (theta <= lam_A[0, n - m:] + 1e-12).all()

U, _ = np.linalg.qr(rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)))
mu, mult = np.array([-1.0, 0.5, 2.0]), [2, 3, 1]          # Lemma A.12 with multiplicities 2, 3, 1
spec = np.repeat(mu, mult)
A1 = U @ np.diag(spec) @ U.conj().T
for i, m_i in enumerate(mult):
    P = np.eye(n, dtype=complex)
    for l in range(len(mu)):
        if l != i:
            P = P @ (A1 - mu[l] * np.eye(n)) / (mu[i] - mu[l])
    Vi = U[:, spec == mu[i]]
    assert np.abs(P - Vi @ Vi.conj().T).max() < 1e-10 and abs(np.trace(P).real - m_i) < 1e-10
```

<a id="s-A-4"></a>

## A.4 高斯向量

*讲义：* [A.4 节](../../lecture_notes/zh/A_analysis.md#s-A-4)；[定理 A.15](../../lecture_notes/zh/A_analysis.md#r-A-15)、[命题 A.17](../../lecture_notes/zh/A_analysis.md#r-A-17)。

[定理 A.15](../../lecture_notes/zh/A_analysis.md#r-A-15) 给出非退化中心化高斯向量 $`X`$ 的密度 $`f_X`$；该代码块把 $`\mathbb P(X\in[0,1]^2)`$
与用中点法则算出的 $`f_X`$ 的积分比较。由[命题 A.17](../../lecture_notes/zh/A_analysis.md#r-A-17)，$`\mathbb E\mathrm{Tr}\,H^2`$ 对 GOE 为
$`n(n+1)/2`$，对 GUE 为 $`n^2`$，并且 $`OHO^T`$ 与 GOE 矩阵 $`H`$ 同分布；该代码块检验其元素的两个矩。

```python
Amat = rng.standard_normal((2, 4))                       # Theorem A.15: X = A g in R^2, A of rank 2
C = Amat @ Amat.T
X = rng.standard_normal((400_000, 4)) @ Amat.T
h = 1 / 400                                              # midpoint rule for the density on [0, 1]^2
grid = (np.arange(400) + 0.5) * h
xx = np.stack(np.meshgrid(grid, grid, indexing='ij'), -1).reshape(-1, 2)
dens = np.exp(-0.5 * np.einsum('ni,ij,nj->n', xx, np.linalg.inv(C), xx)) / (2 * math.pi * math.sqrt(np.linalg.det(C)))
mc(((X >= 0) & (X <= 1)).all(axis=1), dens.sum() * h * h, "P(X in [0,1]^2) against the integral of f_X",
   kind="midpoint rule")

nn, B = 4, 100_000                                       # Proposition A.17
G = rng.standard_normal((B, nn, nn))
goe = (G + G.transpose(0, 2, 1)) / 2
Gc = (rng.standard_normal((B, nn, nn)) + 1j * rng.standard_normal((B, nn, nn))) / math.sqrt(2)
gue = (Gc + Gc.conj().transpose(0, 2, 1)) / math.sqrt(2)
mc(np.einsum('bij,bji->b', goe, goe), nn * (nn + 1) / 2, "E Tr H^2, GOE")      # n + 2 * n(n-1)/2 * 1/2
mc(np.einsum('bij,bji->b', gue, gue).real, nn ** 2, "E Tr H^2, GUE")           # n + 2 * n(n-1)/2 * 1
O, _ = np.linalg.qr(rng.standard_normal((nn, nn)))
rot = O @ goe @ O.T                                                             # has the law of the GOE
mc(rot[:, 0, 0] ** 4, 3.0, "E (O H O^T)_11^4")
mc(rot[:, 0, 1] ** 2, 0.5, "E (O H O^T)_12^2")
```

<a id="s-A-5"></a>

## A.5 高斯向量的方向

*讲义：* [A.5 节](../../lecture_notes/zh/A_analysis.md#s-A-5)；[推论 A.19](../../lecture_notes/zh/A_analysis.md#r-A-19)、[定理 A.20](../../lecture_notes/zh/A_analysis.md#r-A-20)、[推论 A.21](../../lecture_notes/zh/A_analysis.md#r-A-21)。

对标准复高斯向量 $`z\in\mathbb C^3`$，[定理 A.20](../../lecture_notes/zh/A_analysis.md#r-A-20)(b) 指出 $`|z|^2`$ 有密度
$`t^{N-1}e^{-t}/(N-1)!`$，并且 $`z/|z|`$ 在球面上均匀分布且与 $`|z|`$ 独立；[推论 A.21](../../lecture_notes/zh/A_analysis.md#r-A-21) 给出均匀单位向量的矩。
该代码块检验这些结论、单位球的体积 $`V_3=4\pi/3`$（[推论 A.19](../../lecture_notes/zh/A_analysis.md#r-A-19)），以及 $`z_1`$ 的相位是均匀的
（[A.5 节](../../lecture_notes/zh/A_analysis.md#s-A-5)，定理 A.20 之后）。

```python
N, B = 3, 400_000                                   # Theorem A.20 (b), Corollary A.21, in C^3
z = (rng.standard_normal((B, N)) + 1j * rng.standard_normal((B, N))) / math.sqrt(2)
R2 = (np.abs(z) ** 2).sum(axis=1)
u = z / np.sqrt(R2)[:, None]
t = 2.0
p_gamma = 1 - math.exp(-t) * sum(t ** k / math.factorial(k) for k in range(N))   # P(|z|^2 <= t) for Gamma(N)
mc(R2 <= t, p_gamma, "P(|z|^2 <= 2)")
mc(R2 ** 2, N * (N + 1), "E|z|^4")
mc(np.abs(u[:, 0]) ** 4, 2 / (N * (N + 1)), "E|u_1|^4")
mc(np.abs(u[:, 0] * u[:, 1]) ** 2, 1 / (N * (N + 1)), "E|u_1|^2 |u_2|^2")
mc((R2 <= t) * np.abs(u[:, 0]) ** 4, p_gamma * 2 / (N * (N + 1)), "E[1{|z|^2 <= 2} |u_1|^4], independence")
mc((u[:, 0] * u[:, 1].conj()).real, 0.0, "E Re u_1 conj(u_2)")

x3 = rng.uniform(-1, 1, (B, 3))                     # Corollary A.19: V_3 = s_3 / 3 = 4 pi / 3
mc((x3 ** 2).sum(axis=1) <= 1, math.pi / 6, "P(uniform point of [-1,1]^3 in the unit ball)")
g3 = rng.standard_normal((B, 3))
u3 = g3 / np.linalg.norm(g3, axis=1, keepdims=True)
mc(u3[:, 0] ** 4, 3 / 15, "E u_1^4 on S^2")
phase = np.angle(z[:, 0])                           # sigma_1 is normalised arc length
mc((phase >= 0) & (phase < 1), 1 / (2 * math.pi), "P(arg z_1 in [0, 1))")
```

<a id="s-A-6"></a>

## A.6 高斯矩：Isserlis–Wick 公式

*讲义：* [A.6 节](../../lecture_notes/zh/A_analysis.md#s-A-6)；[定理 A.24](../../lecture_notes/zh/A_analysis.md#r-A-24)。

Isserlis–Wick 公式（[定理 A.24](../../lecture_notes/zh/A_analysis.md#r-A-24)）把中心化高斯向量的矩写成对 $`(2p-1)!!`$ 个配对求和的协方差乘积之和。
该代码块枚举这些配对，把公式与 $`\mathbb R^3`$ 中高斯向量的蒙特卡罗矩比较，并由定理的复形式对标准复高斯变量 $`z`$
重新得到 $`\mathbb E|z|^4=2`$。

```python
def pairings(idx):
    """All pairings of the list idx, as lists of pairs."""
    if not idx:
        yield []
        return
    for i in range(1, len(idx)):
        for rest in pairings(idx[1:i] + idx[i + 1:]):
            yield [(idx[0], idx[i])] + rest


def wick(cov, idx):                                       # sum over pairings of the positions of idx
    return sum(math.prod(cov[idx[a], idx[b]] for a, b in p) for p in pairings(list(range(len(idx)))))


assert [len(list(pairings(list(range(2 * p))))) for p in range(1, 5)] == [1, 3, 15, 105]    # (2p-1)!!
Amat = rng.standard_normal((3, 3)) / math.sqrt(3)
cov = Amat @ Amat.T
Y = rng.standard_normal((1_000_000, 3)) @ Amat.T
for idx in ((0, 1, 2, 0), (0, 0, 0, 0), (0, 1, 1, 2, 2, 0)):
    mc(np.prod(Y[:, list(idx)], axis=1), wick(cov, idx), f"E Y{idx}")
mc(np.prod(Y[:, [0, 1, 2]], axis=1), 0.0, "E Y(0, 1, 2), odd")
M = np.array([[1, 1j], [1, -1j]])                  # (z, conj z) = M (x, y) with x, y independent N(0, 1/2)
cz = M @ M.T / 2                                     # E[Z_a Z_b] = [[0, 1], [1, 0]]
assert abs(wick(cz, (0, 1, 0, 1)) - 2) < 1e-12                           # E|z|^4 = 2
```

<a id="s-A-7"></a>

## A.7 指数变量与均匀间距

*讲义：* [A.7 节](../../lecture_notes/zh/A_analysis.md#s-A-7)；[定理 A.25](../../lecture_notes/zh/A_analysis.md#r-A-25)、[推论 A.26](../../lecture_notes/zh/A_analysis.md#r-A-26)。

由[推论 A.26](../../lecture_notes/zh/A_analysis.md#r-A-26)，$`n`$ 个排好序的均匀点的 $`n+1`$ 个间距与归一化指数变量 $`(e_0,\dots,e_n)/T`$ 同分布，
具有 (c) 部分的矩、$`\mathbb P(s_i/s_j\le t)=t/(1+t)`$，以及平均间距比 $`2\ln2-1`$。该代码块对 $`n=5`$ 的两族检验这些
结论，并通过 $`\mathbb E[w_0T^2]=\mathbb Ew_0\,\mathbb ET^2=n+2`$ 检验 $`(e_0,\dots,e_n)/T`$ 与 $`T`$ 的独立性
（[定理 A.25](../../lecture_notes/zh/A_analysis.md#r-A-25)）。

```python
n, B = 5, 200_000                                       # Corollary A.26
levels = np.sort(rng.uniform(0, 1, (B, n)), axis=1)
s = np.diff(np.concatenate([np.zeros((B, 1)), levels, np.ones((B, 1))], axis=1), axis=1)    # s_0, ..., s_n
e = rng.exponential(1.0, (B, n + 1))
w = e / e.sum(axis=1, keepdims=True)                    # normalised exponentials
for name, sp in (("uniform spacings", s), ("normalised exponentials", w)):
    mc(sp[:, 2], 1 / (n + 1), f"{name}: E s_2")
    mc(sp[:, 2] ** 2, 2 / ((n + 1) * (n + 2)), f"{name}: E s_2^2")
    mc(sp[:, 0] * sp[:, n], 1 / ((n + 1) * (n + 2)), f"{name}: E s_0 s_n")
    for tt in (0.5, 2.0):
        mc(sp[:, 1] / sp[:, 2] <= tt, tt / (1 + tt), f"{name}: P(s_1 / s_2 <= {tt})")
r = np.minimum(s[:, 1:-2], s[:, 2:-1]) / np.maximum(s[:, 1:-2], s[:, 2:-1])    # ratios of interior spacings
mc(r.mean(axis=1), 2 * math.log(2) - 1, "<r>, uniform levels")
T = e.sum(axis=1)
mc(w[:, 0] * T ** 2, n + 2, "E[w_0 T^2] = E w_0 E T^2 (independence)")
```

**下一章：** [教程，附录 B. Schur–Weyl 对偶与 Weingarten 演算](B_schur_weyl.md)。
