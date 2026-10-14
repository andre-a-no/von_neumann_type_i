[English](../05_random_matrices.md) · [Español](../es/05_random_matrices.md) · [Français](../fr/05_random_matrices.md) · [Deutsch](../de/05_random_matrices.md) · **中文** · [日本語](../ja/05_random_matrices.md) · [Русский](../ru/05_random_matrices.md)

# 教程 5. 代数中的随机矩阵

本教程配合讲义的[第 5 章](../../lecture_notes/zh/05_random_matrices.md)。该章构造了 $\mathcal M=\bigoplus_cM_{k_c}$ 的酉群上的 Haar 测度
（[推论 5.4](../../lecture_notes/zh/05_random_matrices.md#r-5-4)），计算了它的一阶矩和二阶矩（[命题 5.7](../../lecture_notes/zh/05_random_matrices.md#r-5-7)、[定理 5.10](../../lecture_notes/zh/05_random_matrices.md#r-5-10)），并说明了为什么能级统计
必须逐个扇区计算（[定理 5.15](../../lecture_notes/zh/05_random_matrices.md#r-5-15)）。我们用本库抽取 Haar 酉矩阵和高斯系综，并用蒙特卡罗方法检验这些
结果：均值与精确值之差必须在 $5$ 个标准误差之内。没有定理支撑的结果都标为数值观察。

在本章中，$\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$，$D=\sum_ck_c$，$P_c$ 为
极小中心投影，$E$ 是命题 2.8 中到中心上的条件期望。$U(k)$ 是 $k\times k$ 酉矩阵构成的群，
$U(\mathcal M)=\{U\in\mathcal M:U^*U=1\}$ 是 $\mathcal M$ 的酉群。随机变量定义在一个
概率空间上，期望记为 $\mathbb E$；随机矩阵的期望按元素逐个取。*标准复高斯*变量的实部与虚部
相互独立且 $\sim N(0,\tfrac12)$，从而 $\mathbb E|z|^2=1$（按照 `torch.randn` 的文档，这正是它对复数
dtype 抽取的分布；第一个代码块检验了二阶矩）；*Ginibre 矩阵*的元素是相互独立的标准复高斯变量。

<a id="s-5-1"></a>

## 5.1 Haar 测度

*讲义：* [5.1 节](../../lecture_notes/zh/05_random_matrices.md#s-5-1)；[定理 5.3](../../lecture_notes/zh/05_random_matrices.md#r-5-3)（Ginibre 矩阵 QR 分解的酉因子在用 $\Lambda=\mathrm{diag}(R'_{ii}/|R'_{ii}|)$ 固定相位后服从 Haar 分布）、[推论 5.4](../../lecture_notes/zh/05_random_matrices.md#r-5-4)（$U(\mathcal M)$ 上 Haar 酉元的各块相互独立且服从 Haar 分布）和[推论 5.5](../../lecture_notes/zh/05_random_matrices.md#r-5-5)（$\mathbb EU_{ij}=0$，$\mathbb E|U_{ij}|^2=1/k$）。$SU(k)$ 的 Haar 测度是 $e^{-i\theta/k}U$ 的分布，其中 $\theta=\arg\det U$。

若不作[定理 5.3](../../lecture_notes/zh/05_random_matrices.md#r-5-3) 中的修正 $\Lambda$，$Q'$ 的分布依赖于 QR 例程的符号约定。
对 `torch.linalg.qr` 所调用的例程（LAPACK），下面的代码发现
$\mathbb E\mathrm{Re}\,\mathrm{Tr}\,Q'$ 离 $0$ 有许多个标准误差，而 Haar 酉矩阵满足
$\mathbb EU=0$（推论 5.5）：所以这个 $Q'$ 不服从 Haar 分布。这是关于该例程的陈述，
由代码检验，而不是定理。

本库用 `random_unitary_operator` 抽取 $U(\mathcal M)$ 上的 Haar 酉元：对每种块大小，
一次批量 QR 分解同时处理所有扇区和所有样本（在 GPU 上则用 Gram–Schmidt，由引理 5.2 中的
唯一性，它给出相同的 $Q$）。测度 `'haar'` 就是在每个块中应用定理 5.3；对复代数，
`'haar_su'` 如上乘以 $e^{-i\theta/k}$；`'coe'` 是习题 3 中的系综 $S=W^TW$（本库的其他测度
在本课程中不使用）。我们检验推论 5.4 和 5.5，
并检验本库的采样器恰好就是定理 5.3 中经过修正的 QR：

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, SpinChain, channels

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def mc(samples, exact, name):
    """Monte Carlo mean of real samples with its standard error; must agree within 5 standard errors."""
    mean, err = samples.mean().item(), samples.std().item() / math.sqrt(samples.numel())
    print(f"{name}: {mean:.4f} +- {err:.4f}   (theory {exact:.4f})")
    assert abs(mean - exact) <= 5 * err + 1e-12
    return mean

alg = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')        # M_1 (+) M_2 (+) M_3, D = 6
B = 100_000
U = alg.random_unitary_operator(B)                                   # Haar on U(1) x U(2) x U(3)
assert close(U.adjoint() @ U, alg.identity())                        # unitary elements of the algebra
u = [U.matrix[:, c, :k, :k] for c, k in enumerate(alg.k_factors)]   # the blocks U_c
for c, k in enumerate(alg.k_factors):                                # Corollary 5.5
    mc(u[c][:, 0, 0].real, 0.0, f"E Re (U_{c})_11")
    mc(u[c][:, 0, 0].abs() ** 2, 1 / k, f"E |(U_{c})_11|^2")
mc(u[1][:, 0, 0].abs() ** 2 * u[2][:, 0, 0].abs() ** 2, 1 / 6, "E |(U_1)_11|^2 |(U_2)_11|^2")   # independent blocks

z = torch.randn(B, dtype=dt)                                         # torch.randn draws standard complex Gaussians:
mc(z.real ** 2, 0.5, "E (Re z)^2")                                   # Re z, Im z ~ N(0, 1/2), uncorrelated
mc(z.imag ** 2, 0.5, "E (Im z)^2")
mc(z.real * z.imag, 0.0, "E Re z Im z")

SU = alg.random_unitary_operator(4, measure='haar_su').matrix
assert torch.allclose(torch.linalg.det(SU[:, 2, :3, :3]), torch.ones(4, dtype=dt))

torch.manual_seed(1)
Z = torch.randn(B, 2, 2, dtype=dt)                                   # Ginibre matrices
torch.manual_seed(1)
V = alg.random_unitary(2, batch_size=B)                              # the library, same random numbers
Q, R = torch.linalg.qr(Z)                                            # LAPACK: diag(R) real, either sign
d = R.diagonal(dim1=-2, dim2=-1)
assert close(V, Q * (d / d.abs()).unsqueeze(-2))                     # Theorem 5.3: V = Q Lambda
raw = Q.diagonal(dim1=-2, dim2=-1).sum(-1).real                      # Re Tr Q without the correction
print(f"E Re Tr Q without the phase correction: {raw.mean():.4f} +- {raw.std() / math.sqrt(B):.4f}")
assert abs(raw.mean()) > 20 * raw.std() / math.sqrt(B)               # far from 0: not Haar
```

*数值观察。* 若不作修正，代码对 $k=2$ 得到
$\mathbb E\mathrm{Re}\,\mathrm{Tr}\,Q'\approx-0.84$，而 Haar 酉矩阵满足 $\mathbb EU=0$。

<a id="s-5-2"></a>

## 5.2 矩与 Weingarten 演算

*讲义：* [5.2 节](../../lecture_notes/zh/05_random_matrices.md#s-5-2)；[命题 5.7](../../lecture_notes/zh/05_random_matrices.md#r-5-7)（$\mathbb E[UXU^*]=E(X)$）、[引理 5.8](../../lecture_notes/zh/05_random_matrices.md#r-5-8)（所有 $V\otimes V$ 的交换子是 $\mathrm{span}\{1,F\}$，$F$ 为交换算子）、[定理 5.10](../../lecture_notes/zh/05_random_matrices.md#r-5-10)（二阶 Weingarten 公式）和[推论 5.11](../../lecture_notes/zh/05_random_matrices.md#r-5-11)。

我们将命题 5.7 与 `channels.center_expectation` 对照检验，用线性代数检验引理 5.8（三个抽样得到的
$V\otimes V$ 的交换子作为零空间计算，它包含 $F$ 且是 $2$ 维的，所以由
引理 5.8 它已经就是 $\mathrm{span}\{1,F\}$），并用上面样本中的
$U(3)$ 块通过蒙特卡罗方法检验定理 5.10 和推论 5.11：

```python
X = alg.from_blocks([torch.randn(k, k, dtype=dt) for k in alg.k_factors])      # a fixed element
T = (U @ X @ U.adjoint()).matrix                                    # B samples of U X U*
mean, err = T.mean(0), T.std(0) / math.sqrt(B)
EX = channels.center_expectation(alg)(X).matrix[0]
random = err > 1e-12                                                 # entries that are not constant
z = ((mean - EX).abs()[random] / err[random]).max().item()
print(f"twirl over U(M) vs E(X): largest deviation {z:.2f} standard errors")
assert ((mean - EX).abs() <= 5 * err + 1e-12).all()                 # Proposition 5.7

def commutant_dim(ops, tol=1e-9):
    """Dimension of {X : XA = AX for all A in ops}, as the null space of X -> (AX - XA)_A."""
    d = ops[0].shape[0]
    eye = torch.eye(d, dtype=dt)
    ops = [A.contiguous() for A in ops]                              # torch.kron needs contiguous input
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    return int((torch.linalg.svdvals(L) < tol).sum())

k = 3
F = torch.zeros(k * k, k * k, dtype=dt)
for i in range(k):
    for j in range(k):
        F[j * k + i, i * k + j] = 1                                  # F (e_i (x) e_j) = e_j (x) e_i
VV = [torch.kron(W, W) for W in alg.random_unitary(k, batch_size=3).contiguous()]
assert all(close(F @ W, W @ F) for W in VV)
assert commutant_dim(VV) == 2                                        # Lemma 5.8: span{1, F}

u3 = u[2]                                                            # Haar on U(3)
mc(u3[:, 0, 0].abs() ** 4, 2 / (k * (k + 1)), "E|U_11|^4")
mc(u3[:, 0, 0].abs() ** 2 * u3[:, 0, 1].abs() ** 2, 1 / (k * (k + 1)), "E|U_11|^2 |U_12|^2")
mc((u3[:, 0, 0] * u3[:, 1, 1] * (u3[:, 0, 1] * u3[:, 1, 0]).conj()).real,
   -1 / ((k - 1) * k * (k + 1)), "E U_11 U_22 conj(U_12 U_21)")

n_s = 20_000                                                         # Theorem 5.10 for a random X
X2 = torch.randn(k * k, k * k, dtype=dt)
trX, trFX = X2.diagonal().sum(), (F @ X2).diagonal().sum()
alpha, beta = (trX - trFX / k) / (k ** 2 - 1), (trFX - trX / k) / (k ** 2 - 1)
UU = torch.einsum('bij,bkl->bikjl', u3[:n_s], u3[:n_s]).reshape(n_s, k * k, k * k)   # U (x) U
T2 = UU @ X2 @ UU.mH
mean, err = T2.mean(0), T2.std(0) / math.sqrt(n_s)
assert ((mean - alpha * torch.eye(k * k, dtype=dt) - beta * F).abs() <= 5 * err).all()
```

<a id="s-5-3"></a>

## 5.3 能级统计

*讲义：* [5.3 节](../../lecture_notes/zh/05_random_matrices.md#s-5-3)：GOE $H=(G+G^T)/2$，GUE $H=(G+G^*)/\sqrt2$，以及能级间距 $s_i=\lambda_{i+1}-\lambda_i$ 的间距比 $r_i=\min(s_i,s_{i+1})/\max(s_i,s_{i+1})$；[命题 5.13](../../lecture_notes/zh/05_random_matrices.md#r-5-13)（$2\times2$ 矩阵的 Wigner 猜测）和[命题 5.14](../../lecture_notes/zh/05_random_matrices.md#r-5-14)（平均间距比对独立能级为 $2\ln2-1$，对 $3\times3$ 的 GOE 与 GUE 分别为 $4-2\sqrt3$ 和 $\frac{2\sqrt3}\pi-\frac12$，其中 $\mathbb Er=I_1/I_0$，$I_j=\int_0^1x^jf_\beta(x)\,dx$）。

函数 `SpinChain.level_spacing_ratio`（位于 `chains` 中）返回每一行能级的平均间距比。
我们检验命题 5.13 和 5.14，并记录 $n=64$ 时 GOE 的数值：

```python
ratio = SpinChain.level_spacing_ratio
n2 = 200_000
alg_r = TypeIAlgebra([2], [2], complex_valued=False, precision='double')     # real symmetric 2 x 2
alg_c = TypeIAlgebra([2], [2], precision='double')
G = torch.randn(n2, 1, 2, 2, dtype=torch.float64)
Gc = torch.randn(n2, 1, 2, 2, dtype=dt)
goe2 = alg_r.operator((G + G.mT) / 2, is_self_adjoint=True)
gue2 = alg_c.operator((Gc + Gc.mH) / math.sqrt(2), is_self_adjoint=True)
cdf = {1: lambda t: 1 - math.exp(-math.pi * t ** 2 / 4),
       2: lambda t: math.erf(2 * t / math.sqrt(math.pi)) - 4 * t / math.pi * math.exp(-4 * t ** 2 / math.pi)}
for beta, H, mean_s in ((1, goe2, math.sqrt(math.pi)), (2, gue2, 4 / math.sqrt(math.pi))):
    lam = H.eigenvalues()[0]                                         # (n2, 2), ascending
    S = (lam[:, 1] - lam[:, 0]) / mean_s
    mc(S, 1.0, f"beta={beta}: E S")
    for t in (0.25, 1.0, 2.0):
        mc((S <= t).double(), cdf[beta](t), f"beta={beta}: P(S <= {t})")

levels = torch.rand(4000, 50, dtype=torch.float64)                   # independent levels
mc(ratio(levels), 2 * math.log(2) - 1, "<r>, independent uniform levels")      # Proposition 5.14(a)

n3 = 100_000                                                         # Proposition 5.14(b): 3 x 3
G = torch.randn(n3, 1, 3, 3, dtype=torch.float64)
Gc = torch.randn(n3, 1, 3, 3, dtype=dt)
goe3 = TypeIAlgebra([3], [3], complex_valued=False, precision='double').operator((G + G.mT) / 2, is_self_adjoint=True)
gue3 = TypeIAlgebra([3], [3], precision='double').operator((Gc + Gc.mH) / math.sqrt(2), is_self_adjoint=True)
h = 1e-4
xs = (torch.arange(int(1 / h), dtype=torch.float64) + 0.5) * h      # midpoint rule on [0, 1]
for beta, H, exact in ((1, goe3, 4 - 2 * math.sqrt(3)), (2, gue3, 2 * math.sqrt(3) / math.pi - 0.5)):
    f = (xs + xs ** 2) ** beta / (1 + xs + xs ** 2) ** (1 + 1.5 * beta)
    assert abs((xs * f).sum() / f.sum() - exact) < 1e-7              # E r = I_1 / I_0
    r3 = ratio(H.eigenvalues()[0])                                   # one ratio per matrix
    mc(r3, exact, f"beta={beta}: <r>, 3 x 3")
    mc((r3 <= 0.5).double(), (f * (xs <= 0.5)).sum().item() / f.sum().item(), f"beta={beta}: P(r <= 1/2), 3 x 3")

n = 64                                                               # numerical observation, no theorem
alg64 = TypeIAlgebra([n], [n], complex_valued=False, precision='double')
G = torch.randn(800, 1, n, n, dtype=torch.float64)
lam = alg64.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()[0]
r_goe = ratio(lam[:, n // 4: 3 * n // 4])                            # 800 samples of the mean ratio
print(f"<r>, GOE(64), middle half of the spectrum: {r_goe.mean():.4f} +- {r_goe.std() / math.sqrt(800):.4f}")

def same_mean(a, b, name):
    """Two independent Monte Carlo samples of the same law: means agree within 5 combined standard errors."""
    ea, eb = a.std().item() / math.sqrt(a.numel()), b.std().item() / math.sqrt(b.numel())
    print(f"{name}: {a.mean():.4f} vs {b.mean():.4f} +- {math.hypot(ea, eb):.4f}")
    assert abs(a.mean().item() - b.mean().item()) <= 5 * math.hypot(ea, eb)
```

*数值观察。* 代码对 $n=64$ 的 GOE 谱的中间一半得到 $\langle r\rangle\approx0.53$，低于
$3\times3$ 的精确值 $4-2\sqrt3\approx0.536$，而远高于独立能级的值
$2\ln2-1\approx0.386$。这个数是计算的输出，而不是
定理；下面它作为同样大小 $64$ 的 GOE 块的参考值，这些块的间距比具有
完全相同的分布。

<a id="s-5-4"></a>

## 5.4 为何不能混合扇区

*讲义：* [5.4 节](../../lecture_notes/zh/05_random_matrices.md#s-5-4)；由[定理 5.15](../../lecture_notes/zh/05_random_matrices.md#r-5-15)，不同的独立块的能级之间并不排斥：对两个 GOE 块，合并谱的最小间距 $\delta$ 满足 $P(\delta<\varepsilon)$ 至少是 $\varepsilon$ 的量级，而在一个 GOE 块内部它至多是 $\varepsilon^2$ 的量级。

扇区越多，这一效应越强。我们不对多个扇区证明极限定理；下面的
实验是一个数值观察。在本库中，各块从一开始就是分开的：
`eigenvalues()` 对每个扇区返回一个谱，混合它们需要显式地调用 `torch.cat`。我们把 $C$
个大小为 $64$ 的独立 GOE 块放进一个代数，并比较每个扇区的平均间距比与
合并谱的平均间距比。每个扇区都是大小为 $64$ 的 GOE 矩阵，所以它的间距比与上面的 `r_goe`
具有完全相同的分布，代码检验了这一点；合并谱的行为则只是被观察到：

```python
k, merged_r = 64, {}
for C in (1, 2, 4, 8, 16):
    alg_C = TypeIAlgebra([k] * C, [k] * C, complex_valued=False, precision='double')
    G = torch.randn(800 // C, C, k, k, dtype=torch.float64)
    lam = alg_C.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()    # one spectrum per sector
    per_sector = torch.stack([ratio(l[:, k // 4: 3 * k // 4]) for l in lam])
    merged = torch.cat(lam, dim=-1).sort(dim=-1).values                       # the spectrum on C^D
    r_m = ratio(merged[:, k * C // 4: 3 * k * C // 4])
    merged_r[C] = r_m.mean().item()
    print(f"C = {C:2d}: <r> per sector {per_sector.mean():.4f}, merged {merged_r[C]:.4f}"
          f" +- {r_m.std() / math.sqrt(len(r_m)):.4f}")
    same_mean(per_sector.flatten(), r_goe, f"C = {C:2d}: sectors vs GOE(64)")   # every sector is GOE(64)
assert merged_r[1] > merged_r[2] > merged_r[4] > merged_r[16]       # numerical observation
assert abs(merged_r[16] - (2 * math.log(2) - 1)) < 0.01              # numerical observation
```

*数值观察。* 每个扇区都有 $\langle r\rangle\approx0.53$，即 5.3 节中 GOE(64) 的值
（这一吻合是对一个事实的检验：各扇区都是 GOE(64) 矩阵）。合并谱在两个扇区时降到
$0.42$，在 $8$ 和 $16$ 个扇区时降到约 $0.39$，与独立能级的值 $2\ln2-1$ 相差不到 $0.01$
（命题 5.14(a)）。随着 $C$ 增大合并谱趋近于独立能级，这是代码观察到的现象；本课程并不证明它。

在物理模型中也会发生同样的事情。带随机场的 Heisenberg 链（第 9 章）守恒
磁化强度，所以它的哈密顿量属于 $\bigoplus_{N=0}^LM_{\binom LN}$（命题 9.1）。它的各块*不是*
相互独立的——它们都依赖于同一组随机场——所以定理 5.15 不能照原样应用；
下面的实验显示出同样的效应，作为一个数值观察。

```python
chain = SpinChain(10, complex_valued=False)                         # sectors N = 0, ..., 10
lam = chain.random_field_heisenberg(W=0.5, batch_size=20).eigenvalues()
k5 = chain.sector_dim(5)                                            # 252
r_sector = ratio(lam[5][:, k5 // 4: 3 * k5 // 4]).mean().item()
spectrum = torch.cat(lam, dim=-1).sort(dim=-1).values               # all 1024 levels mixed
r_all = ratio(spectrum[:, 256:768]).mean().item()
print(f"Heisenberg chain, L = 10: <r> in sector N = 5: {r_sector:.4f}, all sectors mixed: {r_all:.4f}")
assert r_sector > 0.49 and r_all < 0.42                             # numerical observation
```

*数值观察。* 在扇区 $N=5$ 内部，代码得到 $\langle r\rangle\approx0.53$，接近上面
GOE(64) 的值，而所有扇区混合后的谱给出 $0.39$，接近 $2\ln2-1$。

<a id="s-5-5"></a>

## 5.5 可与之比较的工具

Haar 酉矩阵：`scipy.stats.unitary_group` 和 `ortho_group`，QuTiP 的 Haar 随机酉矩阵，
以及任何 `numpy`/`torch` 的 QR 再加上相位修正。
Weingarten 演算：符号计算软件包 IntU（Mathematica；Z. Puchała, J. A. Miszczak）、RTNI
（Mathematica 与 Python；M. Fukuda, R. König, I. Nechita, J. Phys. A 52 (2019), 425303）以及
Python 软件包 Haarpy 能符号地计算任意阶的矩；定理 5.10 是它们的二阶情形，
而定理 5.12（在附录 B 中证明）是它们所实现的一般公式。能级统计通常用几行 `numpy` 代码计算，
而 QuSpin 等精确对角化程序逐个扇区地提供谱。本库
增加的是把代数的酉群作为一个对象：在批量调用中抽取所有大小的独立 Haar 块，
把 twirl 等同于条件期望
`channels.center_expectation`，并按扇区返回谱，使得混合扇区永远不是
默认行为。

<a id="x-5"></a>

## 习题

讲义中[习题](../../lecture_notes/zh/05_random_matrices.md#x-5)的数值部分：
1. 对 $U(\mathcal M)$ 上服从 Haar 分布的 $U$，$\mathcal M=M_1\oplus M_2\oplus M_3$，检验 $\mathbb E|\mathrm{Tr}\,U|^2=C=3$ 和
   $\mathbb E|\mathrm{Tr}\,U|^4=17$。
2. 对 $k=3$ 数值计算 $3$ 阶 Weingarten 矩阵（[定理 5.12](../../lecture_notes/zh/05_random_matrices.md#r-5-12)），并检验
   $\mathbb E|U_{11}|^6=\frac1{10}$。
3. **圆正交系综**（`measure='coe'`）是 $S=W^TW$，其中 $W$ 在 $U(k)$ 上服从 Haar 分布。检验本库的 COE
   就是 $W^TW$，$S$ 是酉且对称的，对 $m=1,2,3$ 有 $\mathbb E\mathrm{Tr}\,S^m=0$，并且本征相位的平均密度是均匀的。
   对 $k=32$ 的 COE 计算其本征相位的 $\langle r\rangle$，并与 5.3 节中 GOE 的值比较（一个数值观察）。

<details>
<summary>解答</summary>

*习题 1.*

```python
tr = U.trace                                                        # Tr U = sum_c Tr U_c, shape (B,)
mc(tr.abs() ** 2, 3.0, "E|Tr U|^2")
mc(tr.abs() ** 4, 17.0, "E|Tr U|^4")
```

*习题 2.*

```python
import itertools

def n_cycles(p):                                                    # number of cycles of a permutation
    seen, count = set(), 0
    for i in range(len(p)):
        count += i not in seen
        while i not in seen:
            seen.add(i)
            i = p[i]
    return count

perms = list(itertools.permutations(range(3)))                      # S_3; s.index(j) = s^{-1}(j)
Gram = torch.tensor([[3.0 ** n_cycles([s.index(p[i]) for i in range(3)]) for p in perms] for s in perms],
                    dtype=torch.float64)                            # k^{#cycles(s^{-1} p)}, k = 3
Wg3 = torch.linalg.inv(Gram)                                        # the Weingarten matrix, n = k = 3
assert abs(Wg3.sum().item() - 6 / (3 * 4 * 5)) < 1e-12
mc(u3[:, 0, 0].abs() ** 6, 0.1, "E|U_11|^6")
```

*习题 3.*

   *数值观察。* 代码对 $k=32$ 的 COE 得到 $\langle r\rangle\approx0.53$，与 5.3 节中 GOE(64) 的值
   相差不到 $0.015$。这一吻合是计算的输出；本课程中没有任何
   定理把 COE 的间距统计与 GOE 的间距统计联系起来（关于
   圆系综，见 F. J. Dyson, J. Math. Phys. 3 (1962), 140–156）。

```python
alg_s = TypeIAlgebra([32], [32], precision='double')
torch.manual_seed(2)
S = alg_s.random_unitary_operator(400, measure='coe').matrix[:, 0]
torch.manual_seed(2)
W = alg_s.random_unitary(32, batch_size=400)                         # same random numbers
assert close(S, W.mT @ W)                                            # the library's COE is W^T W
assert close(S, S.mT) and close(S @ S.mH, torch.eye(32, dtype=dt))
for m in (1, 2, 3):                                                  # E Tr S^m = 0
    trm = torch.linalg.matrix_power(S, m).diagonal(dim1=-2, dim2=-1).sum(-1)
    mc(trm.real, 0.0, f"E Re Tr S^{m}")
    mc(trm.imag, 0.0, f"E Im Tr S^{m}")
phases = torch.linalg.eigvals(S).angle().sort(dim=-1).values         # in (-pi, pi]
mc(((phases > 0) & (phases <= 1)).sum(-1).double(), 32 / (2 * math.pi), "E #eigenphases in (0, 1]")
s = torch.diff(torch.cat([phases, phases[:, :1] + 2 * math.pi], dim=-1), dim=-1)   # around the circle
r = torch.minimum(s, s.roll(-1, -1)) / torch.maximum(s, s.roll(-1, -1))
r_coe = r.mean(-1)
print(f"<r>, COE(32) eigenphases: {r_coe.mean():.4f}; GOE(64), Section 5.3: {r_goe.mean():.4f}")
assert abs(r_coe.mean() - r_goe.mean()) < 0.015                      # numerical observation
```

</details>

**下一章：** [教程 6](06_quantum_channels.md)：库中的量子信道。
