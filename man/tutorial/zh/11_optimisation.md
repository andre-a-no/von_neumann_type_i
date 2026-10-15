[English](../11_optimisation.md) · [Español](../es/11_optimisation.md) · [Français](../fr/11_optimisation.md) · [Deutsch](../de/11_optimisation.md) · **中文** · [日本語](../ja/11_optimisation.md) · [Русский](../ru/11_optimisation.md)

# 教程 11. 优化与迹不等式

本教程配合讲义的[第 11 章](../../lecture_notes/zh/11_optimisation.md)，该章证明了三个迹泛函的极值，并解释了为什么对酉算子采样找不到这些极值。
这里用这些定理来检验本库中约束集的参数化（`optimize.PositiveParam`、`UnitaryParam`、`SelfAdjointParam`）与
多起点优化器（`optimize.extremize`）。记 $`\Delta`$ 为 Michelson 对比度，$`\ell(\delta)=(1-\delta)/(1+\delta)`$，三个泛函为：

| | 泛函 | 变量 |
|---|---|---|
| 实验 1 | $`z=\lvert\mathrm{Tr}(XUY)\rvert-\mathrm{Tr}(XY)`$ | $`X,Y\ge0`$，$`U\in\mathcal M`$ 为酉算子 |
| 实验 2 | $`z=\mathrm{Tr}\lvert XY\rvert-\mathrm{Tr}(\lvert X\rvert\,\lvert Y\rvert)`$ | $`X=UX_0`$，$`Y=VY_0`$，$`X_0,Y_0\ge0`$，$`U,V`$ 为酉算子 |
| 实验 3 | $`z=\mathrm{Tr}(Y\lvert X\rvert Y)-\mathrm{Tr}\lvert YXY\rvert`$ | $`X=X^*`$，$`Y\ge0`$ |

其中 $`|X|`$、$`|Y|`$ 归一化为 $`\lambda_{\max}=1`$，它们的对比度作为参数。**包络**为
$`\overline z(\delta)=\sup\{z:\Delta(|X|)=\delta\}`$ 与 $`\underline z(\delta)=\inf\{z:\Delta(|X|)=\delta\}`$，其余
变量均自由变化。

<a id="s-11-1"></a>

## 11.1 约束集及其参数化

*讲义：* [11.1 节](../../lecture_notes/zh/11_optimisation.md#s-11-1)；[引理 11.1](../../lecture_notes/zh/11_optimisation.md#r-11-1)、[命题 11.2](../../lecture_notes/zh/11_optimisation.md#r-11-2)。

具有给定对比度 $`\delta`$ 的 `optimize.PositiveParam` 实现了讲义[11.1 节](../../lecture_notes/zh/11_optimisation.md#s-11-1)中的参数化：它对每个起点
抽取两个不同的槽位 $`s_\pm`$（`spectrum.pin_max`、`spectrum.pin_min`；每个有序对的概率都是 $`1/(D(D-1))`$）和一个
Haar 随机酉算子 $`U_0`$（`unitary.base`，由 `random_unitary_operator` 抽取，[引理 11.10](../../lecture_notes/zh/11_optimisation.md#r-11-10)），并把参数 $`\theta`$
（`spectrum.theta`）与 $`G`$（`unitary.G`）映为
```math
X(\theta,G)=U\mathrm{diag}(\lambda)U^*,\qquad U=U_0\,e^{G-G^*},\qquad
\lambda_{s_+}=1,\quad\lambda_{s_-}=\ell(\delta),\quad\lambda_s=\ell+(1-\ell)\tfrac{1+\sin\theta_s}2\ \text{otherwise}.
```
`UnitaryParam` 就是单独的映射 $`G\mapsto U_0e^{G-G^*}`$，而 `SelfAdjointParam` 再加上符号，$`X=U\mathrm{diag}(\sigma\lambda)U^*`$。
由[命题 11.2](../../lecture_notes/zh/11_optimisation.md#r-11-2)，其像是对比度壳层 $`\mathcal P_\delta`$ 中在所抽槽位上取值 $`1`$ 与 $`\ell(\delta)`$ 的算子全体。
$`[a,b]`$ 中的自由对比度只是多一个参数；本库的默认值是 $`[a,b]=[0,0.999]`$，下面的代码保留这一默认值。

代码检验约束，并以构造的方式检验命题 11.2：证明中的参数重现了一个目标算子。

```python
import math
import scipy.linalg
import torch
from torch_vn_algebra import TypeIAlgebra, optimize as opt

torch.manual_seed(0)
alg = TypeIAlgebra([3, 3], [3, 3], precision='double')          # M_3 (+) M_3, D = 6
ell = lambda delta: (1 - delta) / (1 + delta)                     # lambda_min for contrast delta
d = torch.tensor([0.0, 0.25, 0.5, 0.9], dtype=torch.float64)
x = opt.PositiveParam(alg, 4, delta=d).operator()                 # four starts, one contrast each
assert torch.allclose(x.michelson_contrast, d, atol=1e-12)        # the constraint holds exactly
assert torch.allclose(x.lambda_max, torch.ones(4, dtype=torch.float64), atol=1e-12)

P = opt.PositiveParam(alg, 1, delta=torch.tensor([0.5], dtype=torch.float64))
lo = ell(0.5)
mu = lo + (1 - lo) * torch.rand(1, 2, 3, dtype=torch.float64)    # target spectrum in [lo, 1] ...
(cp, ip), (cm, im) = P.spectrum.pin_max[0].tolist(), P.spectrum.pin_min[0].tolist()
mu[0, cp, ip], mu[0, cm, im] = 1.0, lo                            # ... with 1 and lo in the drawn slots
W = alg.random_unitary_operator(1).matrix                         # target eigenvectors
target = W @ torch.diag_embed(mu.to(W.dtype)) @ W.mH
with torch.no_grad():                                             # the parameters of the proof
    P.spectrum.theta.copy_(torch.asin((2 * (mu - lo) / (1 - lo) - 1).clamp(-1, 1)))
    V = P.unitary.base.mH @ W                                     # U0^* W = exp(S), G = S / 2
    S = torch.stack([torch.from_numpy(scipy.linalg.logm(V[0, c].numpy())) for c in range(2)])
    P.unitary.G.copy_(S[None] / 2)
err = (P.operator().matrix - target).abs().max().item()
print(f"slots of 1 and lo: {(cp, ip)}, {(cm, im)};  reconstruction error {err:.1e}")
assert err < 1e-10
```

<a id="s-11-2"></a>

## 11.2 冯·诺依曼迹不等式

*讲义：* [11.2 节](../../lecture_notes/zh/11_optimisation.md#s-11-2)；[定理 11.4](../../lecture_notes/zh/11_optimisation.md#r-11-4)、[推论 11.5](../../lecture_notes/zh/11_optimisation.md#r-11-5)、[命题 11.6](../../lecture_notes/zh/11_optimisation.md#r-11-6)、[推论 11.7](../../lecture_notes/zh/11_optimisation.md#r-11-7)。

由[讲义，定理 11.4](../../lecture_notes/zh/11_optimisation.md#r-11-4)（冯·诺依曼迹不等式），$`|\mathrm{Tr}(AB)|\le\sum_is_i(A)s_i(B)`$；由[推论 11.5](../../lecture_notes/zh/11_optimisation.md#r-11-5)，
对特征值为 $`p,q`$ 的自伴矩阵 $`P,Q`$ 有 $`\sum_ip^\downarrow_iq^\uparrow_i\le\mathrm{Tr}(PQ)\le\sum_ip^\downarrow_iq^\downarrow_i`$，对 $`P,Q\ge0`$ 有 $`\mathrm{Tr}(PQ)\le\|PQ\|_1\le\sum_ip^\downarrow_iq^\downarrow_i`$。
对固定谱时的实验 2，[命题 11.6](../../lecture_notes/zh/11_optimisation.md#r-11-6) 给出 $`\sup z=-\inf z=\sum_c\sum_ix^\downarrow_{c,i}(y^\downarrow_{c,i}-y^\uparrow_{c,i})`$，它在 $`X_0=\mathrm{diag}(x^\downarrow)`$、$`V=J`$（逆序置换）以及 $`Y_0=\mathrm{diag}(y^\uparrow)`$
（相应地 $`\mathrm{diag}(y^\downarrow)`$）处取到；[推论 11.7](../../lecture_notes/zh/11_optimisation.md#r-11-7) 给出包络 $`\sup z=-\inf z=\sum_c\lfloor k_c/2\rfloor(1-\ell_X)(1-\ell_Y)`$。

这里对 $`\Delta(|X|)=\delta`$ 有 $`1-\ell_X=2\delta/(1+\delta)`$；当 $`k_c`$ 为奇数时，中间的特征值无关紧要。代码
对 $`k=3`$、$`C=2`$ 检验定理 11.4、推论 11.5 以及命题 11.6 的极值构型。

```python
M5 = TypeIAlgebra([5], [5], precision='double')
A = M5.operator(torch.randn(2000, 1, 5, 5, dtype=torch.complex128))         # Ginibre, no structure
B = M5.operator(torch.randn(2000, 1, 5, 5, dtype=torch.complex128))
sA, sB = torch.linalg.svdvals(A.matrix[:, 0]), torch.linalg.svdvals(B.matrix[:, 0])   # decreasing
assert ((A @ B).trace.abs() <= (sA * sB).sum(-1) + 1e-12).all()             # Theorem 11.4
Ph, Qh = (M5.operator_from_eigenvalues(lambda k: torch.randn(2000, k), batch_size=2000, force_self_adjoint=True)
          for _ in range(2))
p, q = Ph.eigenvalues()[0], Qh.eigenvalues()[0]                             # increasing order
t = (Ph @ Qh).trace.real
assert (t <= (p * q).sum(-1) + 1e-12).all() and (t >= (p * q.flip(-1)).sum(-1) - 1e-12).all()   # Cor. 11.5(a)
Pp, Qp = Ph @ Ph, Qh @ Qh                                                   # positive, spectra p^2, q^2
nuc = (Pp @ Qp).trace_norm()
assert ((Pp @ Qp).trace.real <= nuc + 1e-12).all()                          # Cor. 11.5(b)
assert (nuc <= ((p ** 2).sort(-1)[0] * (q ** 2).sort(-1)[0]).sum(-1) + 1e-10).all()

alg3 = TypeIAlgebra([3, 3], [3, 3], precision='double')
x_dn = torch.tensor([1.0, 0.6, ell(0.5)], dtype=torch.float64)              # Delta = 0.5, middle value 0.6
y_dn = torch.tensor([1.0, 0.7, ell(0.8)], dtype=torch.float64)              # Delta = 0.8, middle value 0.7
J = torch.eye(3, dtype=torch.float64).flip(0)                               # reversal permutation
blocks = lambda m: alg3.from_blocks([m, m])                                 # the same block in both sectors

def z_exp2(X0, V, Y0):
    return ((X0 @ V @ Y0).trace_norm() - (X0 @ Y0).trace.real).item()

X0 = blocks(torch.diag(x_dn))
z_max = z_exp2(X0, blocks(J), blocks(torch.diag(y_dn.flip(0))))
z_min = z_exp2(X0, blocks(J), blocks(torch.diag(y_dn)))
exact = 2 * (x_dn * (y_dn - y_dn.flip(0))).sum().item()                     # Proposition 11.6
envelope = 2 * (1 - ell(0.5)) * (1 - ell(0.8))                              # Corollary 11.7, floor(3/2) = 1
print(f"Exp. 2, k = 3, C = 2: z_max = {z_max:.6f}, z_min = {z_min:.6f}, formula {envelope:.6f}")
assert abs(z_max - exact) < 1e-12 and abs(z_min + exact) < 1e-12 and abs(exact - envelope) < 1e-12
```

<a id="s-11-3"></a>

## 11.3 一个带有取等条件的迹不等式

*讲义：* [11.3 节](../../lecture_notes/zh/11_optimisation.md#s-11-3)；[引理 11.8](../../lecture_notes/zh/11_optimisation.md#r-11-8)、[命题 11.9](../../lecture_notes/zh/11_optimisation.md#r-11-9)。

由[讲义，命题 11.9](../../lecture_notes/zh/11_optimisation.md#r-11-9)，在实验 3 中 $`z=\mathrm{Tr}(Y|X|Y)-\mathrm{Tr}|YXY|\ge0`$，等号成立当且仅当 $`YX_+Y^2X_-Y=0`$，特别地，当 $`XY=YX`$ 时成立。
对 $`X=\mathrm{diag}(1,-1,0)`$ 以及行为 $`(1,0,a)`$, $`(0,1,0)`$, $`(a,0,1)`$, $`-1<a<1`$ 的 $`Y`$，虽然当 $`a\ne0`$ 时 $`XY\ne YX`$，等号仍然成立。代码对两者都作了检验
（取 $`a=0.5`$）；`X.abs()` 即 $`|X|`$，`trace_norm()` 即 $`\|\cdot\|_1`$。

```python
def z_exp3(X, Y):
    return (Y @ X.abs() @ Y).trace.real - (Y @ X @ Y).trace_norm()

Xs = alg3.operator_from_eigenvalues(lambda k: torch.randn(500, k), batch_size=500, force_self_adjoint=True)
Yp = alg3.operator_from_eigenvalues(lambda k: torch.rand(500, k), batch_size=500, force_positive=True)
assert (z_exp3(Xs, Yp) >= -1e-12).all()                                     # Proposition 11.9

M3 = TypeIAlgebra([3], [3], precision='double')
Xe = M3.operator(torch.diag(torch.tensor([1.0, -1.0, 0.0]))[None, None])
Ye = M3.operator(torch.tensor([[1.0, 0, 0.5], [0, 1, 0], [0.5, 0, 1]])[None, None])
print(f"equality without commutation: z = {z_exp3(Xe, Ye).item():.1e}, "
      f"||XY - YX|| = {(Xe @ Ye - Ye @ Xe).operator_norm().item():.2f}")
assert abs(z_exp3(Xe, Ye).item()) < 1e-12 and (Xe @ Ye - Ye @ Xe).operator_norm().item() > 0.1
```

<a id="s-11-4"></a>

## 11.4 采样与优化

*讲义：* [11.4 节](../../lecture_notes/zh/11_optimisation.md#s-11-4)；[引理 11.10](../../lecture_notes/zh/11_optimisation.md#r-11-10)、[命题 11.11](../../lecture_notes/zh/11_optimisation.md#r-11-11)、[推论 11.12](../../lecture_notes/zh/11_optimisation.md#r-11-12)、[引理 11.13](../../lecture_notes/zh/11_optimisation.md#r-11-13)。

与其最大值 $`\|A\|_1`$ 相比，$`|\mathrm{Tr}(AU)|`$ 的典型大小是多少？本库把 $`U`$ 取为
**Ginibre 矩阵** $`Z`$（元素相互独立，密度为
$`\pi^{-1}e^{-|z|^2}`$；Mezzadri 的方法）的 QR 分解 $`Z=QR`$（$`R_{ii}>0`$）中的因子 $`Q`$，在每个块中独立抽取。

由[讲义，命题 11.11](../../lecture_notes/zh/11_optimisation.md#r-11-11)，对 Haar 分布的 $`U\in U(k)`$ 有 $`\mathbb E|\mathrm{Tr}(AU)|^2=\|A\|_2^2/k`$；由[推论 11.12](../../lecture_notes/zh/11_optimisation.md#r-11-12)，
对谱含于 $`[a,b]`$、$`a>0`$ 的 $`X,Y\ge0`$，$`|\mathrm{Tr}(XUY)|`$ 的均方根至多为它关于 $`U`$ 的最大值的 $`\min(1/\sqrt k,b^2/(a^2k))`$ 倍，并且 $`\mathbb P\big(|\mathrm{Tr}(XUY)|\ge\mathrm{Tr}(XY)\big)\le b^4/(a^4k^2)`$。

因此，$`N`$ 个 Haar 样本找到 $`|\mathrm{Tr}(XUY)|\le\mathrm{Tr}(XY)`$ 的反例的概率至多为
$`Nb^4/(a^4k^2)`$（$`N`$ 个事件的并），尽管对每一对不可交换的矩阵这样的反例都存在（命题 3.9）。这个界
只有在 $`k>\sqrt N\,b^2/a^2`$ 时才有信息量（对 $`b/a=6`$ 与 $`N=2000`$：$`k>1609`$）；在数值上，样本在小得多的 $`k`$
下就已经漏掉了反例（见下文）。代码用 Monte Carlo 检验命题 11.11（对这次固定种子的运行，样本均值与精确值
之差在五个标准误之内），并对在 $`[0.2,1.2]`$ 中均匀分布的谱（$`b/a=6`$）打印
推论 11.12 中的比值、$`\sup z`$ 以及 $`2000`$ 个样本中最好的 $`z`$，后两者都相对于 $`\mathrm{Tr}(XY)`$。

```python
def random_pair(alg_k):
    pos = lambda k: 0.2 + torch.rand(1, k)                                  # spectrum in [a, b] = [0.2, 1.2]
    return tuple(alg_k.operator_from_eigenvalues(pos, force_positive=True) for _ in range(2))

M6 = TypeIAlgebra([6], [6], precision='double')
Xf, Yf = random_pair(M6)
N = 20_000
t2 = (Yf @ Xf @ M6.random_unitary_operator(N)).trace.abs() ** 2
m2, err, exact2 = t2.mean().item(), t2.std().item() / math.sqrt(N), (Yf @ Xf).frobenius_norm().item() ** 2 / 6
print(f"E|Tr(AU)|^2 = {m2:.4f} +- {err:.4f}   (exact ||A||_2^2 / k = {exact2:.4f})")
assert abs(m2 - exact2) < 5 * err                                           # Proposition 11.11

print(" k   rms/max   sup z/TrXY   best of 2000 samples")
for k in (2, 4, 8, 16, 32):
    Mk = TypeIAlgebra([k], [k], precision='double')
    Xk, Yk = random_pair(Mk)
    trXY, nuc = (Xk @ Yk).trace.real.item(), (Yk @ Xk).trace_norm().item()
    ratio = (Yk @ Xk).frobenius_norm().item() / (math.sqrt(k) * nuc)
    best = (Xk @ Mk.random_unitary_operator(2000) @ Yk).trace.abs().max().item()
    print(f"{k:2d}   {ratio:.4f}    {(nuc - trXY) / trXY:.2e}     {(best - trXY) / trXY:+.3f}")
    assert ratio <= min(1 / math.sqrt(k), 36 / k) + 1e-12                   # Corollary 11.12
    assert nuc > trXY and best < nuc and (k < 4 or best < trXY)            # from k = 4 on: no violation found
```

*数值观察*（打印出的表格；只对这次固定种子的运行作断言）：该比值始终低于已证明的界
$`\min(1/\sqrt k,b^2/(a^2k))`$，并且 $`k`$ 每翻一倍它大约减半，即按 $`1/k`$ 下降；而从 $`k=4`$ 起，$`2000`$ 个样本中
最好的那个甚至达不到 $`z=0`$。与此相反，在下面的运行中梯度搜索接近了两个极值（这同样是数值观察）。
对于矩阵对 `Xf, Yf`，关于 $`U`$ 的两个极值都是已知的：

由[讲义，引理 11.13](../../lecture_notes/zh/11_optimisation.md#r-11-13)，对每个 $`A\in M_k(\mathbb C)`$、$`k\ge2`$ 有 $`\min_U|\mathrm{Tr}(AU)|=0`$。

所以 $`\inf_Uz=-\mathrm{Tr}(XY)`$，$`\sup_Uz=\|YX\|_1-\mathrm{Tr}(XY)`$（命题 3.9）。调用一次 `extremize`
就能同时逼近二者：$`12`$ 个起点最大化 $`z`$，$`12`$ 个起点最大化 $`-z`$（在这次运行中相对精度达到 $`10^{-3}`$，
下面有断言）。由命题 11.11，$`|\mathrm{Tr}(XUY)|`$ 的均方根为 $`\|YX\|_2/\sqrt k`$；利用推论 11.12 的证明中的
界（谱含于 $`[a,b]`$），$`\mathbb E|\mathrm{Tr}(XUY)|^2\le b^4`$ 且 $`\mathrm{Tr}(XY)\ge ka^2`$，所以 Markov 不等式给出
$`\mathbb P\big(z\ge-(1-\varepsilon)\mathrm{Tr}(XY)\big)\le b^4/(\varepsilon^2a^4k^2)`$（$`\varepsilon>0`$）：对大的 $`k`$，
$`z`$ 以高概率接近其下确界。
*数值观察*（对这次固定种子的运行作断言）：采样到达了下确界的 $`1\,\%`$ 以内，却连 $`z=0`$ 都达不到，
而上确界是正的：

```python
trXY = (Xf @ Yf).trace.real
sup_exact, inf_exact = ((Yf @ Xf).trace_norm() - trXY).item(), -trXY.item()
U6 = opt.UnitaryParam(M6, 24)                                               # 24 starts
pm = torch.tensor([1.0, -1.0], dtype=torch.float64).repeat_interleave(12)   # 12 maximise z, 12 maximise -z
best = opt.extremize(lambda: pm * ((Xf @ U6.operator() @ Yf).trace.abs() - trXY), [U6],
                     maximize=True, rounds=6, steps=30).reshape(2, 12).max(1)[0]
best_up, best_down = best[0].item(), -best[1].item()
z_mc = (Xf @ M6.random_unitary_operator(N) @ Yf).trace.abs() - trXY
print(f"sup z: exact {sup_exact:.6e}, ascent {best_up:.6e}, best sample {z_mc.max().item():+.4f}")
print(f"inf z: exact {inf_exact:.6f}, descent {best_down:.6f}, worst sample {z_mc.min().item():+.6f}")
assert best_up <= sup_exact + 1e-12 and best_down >= inf_exact - 1e-12    # attained values: one-sided
assert sup_exact - best_up < 1e-3 * sup_exact and best_down - inf_exact < 1e-3 * abs(inf_exact)
assert z_mc.min().item() - inf_exact < 1e-2 * abs(inf_exact)               # sampling finds the typical side ...
assert z_mc.max().item() < 0                                               # ... but misses the violation
```

<a id="s-11-5"></a>

## 11.5 包络的多起点搜索

*讲义：* [11.5 节](../../lecture_notes/zh/11_optimisation.md#s-11-5)；[命题 11.14](../../lecture_notes/zh/11_optimisation.md#r-11-14)、[推论 11.7](../../lecture_notes/zh/11_optimisation.md#r-11-7)、[命题 11.9](../../lecture_notes/zh/11_optimisation.md#r-11-9)、[命题 11.15](../../lecture_notes/zh/11_optimisation.md#r-11-15)、[推论 11.16](../../lecture_notes/zh/11_optimisation.md#r-11-16)。

`optimize.extremize(objective, params, maximize)` 运行 Adam，共 `rounds` 轮，每轮 `steps` 步（学习率为
`lr * decay**round`），损失函数为 $`\mp\sum_b\text{objective}_b`$，每轮之后对酉坐标卡重新居中，并**对每个起点**
返回目标函数在梯度迭代点处以及在最终点处（计算一次）所取各值中的最大值（`maximize=False` 时为最小值）。
批中的每个元素都是一个独立的起点，有自己的槽位、符号、$`U_0`$ 与对比度：Adam 按坐标逐个作用，所以如果
目标函数的第 $`b`$ 个分量只依赖于第 $`b`$ 组参数（下面所有的目标函数都是如此），各起点之间互不影响。
`SelfAdjointParam` 的符号是离散的：在两轮之间，`discrete_step` 对每个起点依次尝试八次单个符号的
翻转（每次都在一个均匀随机的槽位上），并保留每一次使目标函数改进的翻转。

在精确算术中，每一步迭代都是容许的（命题 11.2 及其后的说明），所以返回的每个值都是 $`z`$ 在某个容许点处的值；
在浮点运算中，它在舍入误差范围内等于该值（下面的断言允许 $`10^{-12}`$）。因此最大化给出上确界的
**下界**，最小化给出下确界的上界，仅此而已——局部搜索可能停在局部极值处。另一个方向需要证明，
如推论 11.7 与命题 11.9，或者凸松弛（11.6 节）。计算得到的包络要接受以下检验：(i) 与闭式表达比较；(ii) 与
手工构造的构型比较，例如较小解的直和，包络不得低于它们；(iii) 与
结构化样本（两能级谱、置换基）比较，它们的值是另外一些可取到的值。

由[讲义，命题 11.14](../../lecture_notes/zh/11_optimisation.md#r-11-14)，若所有 $`k_c\ge2`$，则实验 1 的下包络为 $`\underline z(\delta)=-\big(D-2\delta/(1+\delta)\big)`$。

代码在 $`M_2\oplus M_2`$ 上检验命题 11.14，在 $`M_3\oplus M_3`$ 上检验推论 11.7（两个对比度都给定，
$`k=3`$ 为奇数），并检验命题 11.9；与上面一样，上确界和下确界在同一次调用中求出。单侧断言（没有任何返回值
越过已证明的极值）是这些定理的推论；而断言搜索到达极值的给定容差以内，则是这次固定种子运行的数值观察。
关键字 `is_self_adjoint` 选择
迹范数的算法（奇异值，或对 $`YXY`$ 用特征值），而无需在每一步都做检验。最后一项
检验属于类型 (ii)：实验 3 的 $`z`$ 对直和是可加的，而直和的对比度是各直和项对比度中
最大者（习题 3），所以对于实验 3 的包络（如代码中那样，$`\Delta(Y)`$ 在 $`[0,0.999]`$ 中自由变化），
```math
\overline z_{M_2\oplus M_2}(\delta)\ \ge\ \overline z_{M_2}(\delta)+\sup_{\delta'\le\delta}\overline z_{M_2}(\delta');
```
特别地，在 $`M_2`$ 上于对比度 $`\delta`$ 与 $`\delta'\le\delta`$ 处取到的两个值之和，是在 $`M_2\oplus M_2`$ 上于对比度
$`\delta`$ 处取到的值。

```python
torch.manual_seed(1)
S = 6                                                                       # starts per grid point
grid = torch.tensor([0.0, 0.5, 0.9], dtype=torch.float64)
alg22 = TypeIAlgebra([2, 2], [2, 2], precision='double')                    # D = 4

# Exp. 1, lower envelope along Delta(X); Delta(Y) free in [0, 0.999]
X1, Y1, U1 = (opt.PositiveParam(alg22, 3 * S, delta=grid.repeat_interleave(S)),
              opt.PositiveParam(alg22, 3 * S), opt.UnitaryParam(alg22, 3 * S))
def z1():
    x, y = X1.operator(), Y1.operator()
    return (x @ U1.operator() @ y).trace.abs() - (x @ y).trace.real
inf1 = opt.extremize(z1, [X1, Y1, U1], maximize=False, rounds=4, steps=30).reshape(3, S).min(1)[0]
exact1 = -(4 - 2 * grid / (1 + grid))                                       # Proposition 11.14
print("Exp. 1  inf:", [f"{v:.4f}" for v in inf1.tolist()], " exact:", [f"{v:.4f}" for v in exact1.tolist()])
assert (inf1 >= exact1 - 1e-12).all() and (inf1 - exact1).max() < 4e-3

# Exp. 2 on M_3 (+) M_3, Delta(X0) = 0.5, Delta(Y0) = 0.8: sup and inf in one call
sign = torch.tensor([1.0, -1.0], dtype=torch.float64).repeat_interleave(S)
X2, Y2, V2 = (opt.PositiveParam(alg3, 2 * S, delta=torch.full((2 * S,), 0.5, dtype=torch.float64)),
              opt.PositiveParam(alg3, 2 * S, delta=torch.full((2 * S,), 0.8, dtype=torch.float64)),
              opt.UnitaryParam(alg3, 2 * S))
def z2():
    x0, y0 = X2.operator(), Y2.operator()
    M = alg3.operator((x0 @ V2.operator() @ y0).matrix, is_self_adjoint=False)
    return sign * (M.trace_norm() - (x0 @ y0).trace.real)
best2 = opt.extremize(z2, [X2, Y2, V2], maximize=True, rounds=4, steps=30).reshape(2, S).max(1)[0]
sup2, inf2 = best2[0].item(), -best2[1].item()
print(f"Exp. 2  sup {sup2:.6f}, inf {inf2:.6f}, exact +-{envelope:.6f}")
assert sup2 <= envelope + 1e-12 and inf2 >= -envelope - 1e-12               # Corollary 11.7: one-sided ...
assert envelope - sup2 < 1e-3 * envelope and inf2 + envelope < 1e-3 * envelope   # ... and reached

# Exp. 3: maximise sign * z at the contrasts `deltas` of |X|; Delta(Y) free in [0, 0.999]
def exp3(alg_, deltas, sign):
    X3, Y3 = opt.SelfAdjointParam(alg_, len(deltas), delta=deltas), opt.PositiveParam(alg_, len(deltas))
    def z3():
        y = Y3.operator()
        M = alg_.operator((y @ X3.operator() @ y).matrix, is_self_adjoint=True)
        return sign * ((y @ X3.abs_operator() @ y).trace.real - M.trace_norm())
    return opt.extremize(z3, [X3, Y3], maximize=True, rounds=4, steps=30).reshape(2, S).max(1)[0]

best3 = exp3(alg22, torch.full((2 * S,), 0.5, dtype=torch.float64), sign)
sup_M22, inf3 = best3[0].item(), -best3[1].item()
s_M2 = exp3(TypeIAlgebra([2], [2], precision='double'),
            torch.tensor([0.0, 0.5], dtype=torch.float64).repeat_interleave(S), 1.0)
direct_sum = (s_M2[1] + s_M2.max()).item()                                 # delta = 0.5, best delta' <= 0.5
closed = [(1 - dl) * (2 * 0.999 / 1.999) ** 2 for dl in (0.0, 0.5)]        # Corollary 11.16, d = 0.999
print(f"Exp. 3  inf {inf3:.1e};  sup for M_2 at delta = 0, 0.5: {s_M2[0]:.5f}, {s_M2[1]:.5f}"
      f"  (closed form {closed[0]:.5f}, {closed[1]:.5f})")
print(f"        sup for M_2 (+) M_2 at 0.5: search {sup_M22:.5f}, direct sum {direct_sum:.5f}")
assert inf3 >= -1e-12 and inf3 < 1e-6                                       # Proposition 11.9
assert sup_M22 >= direct_sum - 5e-3 * direct_sum                           # (ii), up to 0.5 %
assert (s_M2 <= torch.tensor(closed, dtype=torch.float64) + 1e-12).all()   # Corollary 11.16: one-sided
assert sup_M22 <= (2 - 0.5) * (2 * 0.999 / 1.999) ** 2 + 1e-12
```

*数值观察*（这次运行打印出的值）。在 $`M_2\oplus M_2`$ 上的搜索比直和低约 $`0.2\,\%`$：它尚未完全收敛，而
$`\overline z_{M_2\oplus M_2}(0.5)`$ 的最佳下界是较大的那个数——发现这类情况正是检验 (ii)
的目的。$`M_2`$ 的值与 $`(1-\delta)\big(2d/(1+d)\big)^2=(1-\delta)(1-\ell(d))^2`$ 吻合到 $`10^{-4}`$，其中
$`d=0.999`$ 是 $`Y`$ 所允许的最大对比度。这个闭式最早由配套论文在数值上观察到，它就是
[推论 11.16](../../lecture_notes/zh/11_optimisation.md#r-11-16)（见下文）的上包络；代码只断言已证明的方向（没有值超过它）。
由同一推论，直和是最优的：$`\overline z_{M_2\oplus M_2}(0.5)=1.5\,(2d/(1+d))^2`$。对更大的块，我们不知道实验 1 与
实验 3 的上包络的闭式；此时计算得到的曲线只是下界（可取到的值），仅此而已。

由[讲义，命题 11.15](../../lecture_notes/zh/11_optimisation.md#r-11-15)，在 $`M_2`$ 上，对 $`X`$（实验 3 中为 $`|X|`$）的谱
$`x_1\ge x_2`$ 与 $`Y`$ 的谱 $`y_1\ge y_2`$，$`\sup z`$ 在实验 1 中为
$`(x_1-x_2)^2(y_1-y_2)^2/(4(x_1+x_2)(y_1+y_2))`$，在实验 3 中为 $`\frac{2x_1x_2}{x_1+x_2}(y_1-y_2)^2`$。设
$`e_1`$ 为 $`X`$ 对应于模为 $`x_1`$ 的特征值的本征向量，$`f_1`$ 为 $`Y`$ 对应于 $`y_1`$ 的本征向量，则当
$`\gamma=|\langle e_1,f_1\rangle|^2`$ 等于 $`\gamma_1=\frac14+\frac S{2q}`$（实验 1，且 $`|\mathrm{Tr}(XUY)|=\|YX\|_1`$）或
$`\gamma_3=\frac Sq`$（实验 3，$`X=\mathrm{diag}(x_1,-x_2)`$）时上确界被取到，其中 $`S=x_1y_2+x_2y_1`$，
$`q=(x_1+x_2)(y_1+y_2)`$。由
[推论 11.16](../../lecture_notes/zh/11_optimisation.md#r-11-16)，对 $`C`$ 个大小为 $`2`$ 的块，若 $`\Delta(|X|)=\delta`$、
$`\Delta(Y)=\delta'`$，则实验 1 有 $`\sup z=C\delta^2\delta'^2/((1+\delta)(1+\delta'))`$ 与
$`\inf z=-(2C-1+\ell(\delta)\ell(\delta'))`$，实验 3 有 $`\sup z=(C-\delta)\big(2\delta'/(1+\delta')\big)^2`$。两个
上确界都随 $`\delta'`$ 递增：当 $`\Delta(Y)`$ 在 $`[0,d]`$ 中自由变化时，上包络就是 $`\delta'=d`$ 处的值，
当 $`d\to1`$ 时它们趋于 $`C\delta^2/(2(1+\delta))`$ 与 $`C-\delta`$，即奇异 $`Y`$
（$`\Delta(Y)=1`$）时的值。

代码对 $`C=1,2`$ 构造证明中的实极大点（`from_blocks`），检验对比度与数值，在随机复矩阵对上检验命题 11.15
的界，在每个块中用秩一投影取到 $`d\to1`$ 的极限，并对 $`M_2`$ 上的实验 1 沿 $`\Delta(X)`$ 运行多起点搜索
（$`\Delta(Y)`$ 在 $`[0,0.999]`$ 中自由），其结果必须不超过闭式（一个定理），且在这次运行中与之相差不到
$`10^{-3}`$。

```python
def rot(gamma):                                                             # real basis f with |<e_1, f_1>|^2 = gamma
    t = math.acos(math.sqrt(gamma))
    return torch.tensor([[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]], dtype=torch.float64)

def gammas(x1, x2, y1, y2):                                                 # gamma_1 (Exp. 1), gamma_3 (Exp. 3)
    S, q = x1 * y2 + x2 * y1, (x1 + x2) * (y1 + y2)
    return 0.25 + S / (2 * q), S / q

env1 = lambda C, dl, dp: C * dl ** 2 * dp ** 2 / ((1 + dl) * (1 + dp))     # Corollary 11.16(a)
env3 = lambda C, dl, dp: (C - dl) * (2 * dp / (1 + dp)) ** 2               # Corollary 11.16(b)

def maximisers(C, dl, dp):
    """Exp. 1: (X, U, Y) and Exp. 3: (X, Y) attaining Corollary 11.16 on C blocks of size 2."""
    a2 = TypeIAlgebra([2] * C, [2] * C, complex_valued=False, precision='double')
    lx, ly = ell(dl), ell(dp)
    R = rot(gammas(1.0, lx, 1.0, ly)[0])
    X0 = torch.diag(torch.tensor([1.0, lx], dtype=torch.float64))
    Y0 = R @ torch.diag(torch.tensor([1.0, ly], dtype=torch.float64)) @ R.T
    W, _, Vh = torch.linalg.svd(Y0 @ X0)                                    # Y X = W Sigma V^T, U = V W^T
    exp1 = [a2.from_blocks([m] * C) for m in (X0, Vh.T @ W.T, Y0)]
    xs = [lx] + [1.0] * (C - 1)                                             # |X|: {1, ell} in block 1, {1, 1} else
    Rs = [rot(gammas(1.0, x2, 1.0, ly)[1]) for x2 in xs]
    X3 = a2.from_blocks([torch.diag(torch.tensor([1.0, -x2], dtype=torch.float64)) for x2 in xs])
    Y3 = a2.from_blocks([Rc @ torch.diag(torch.tensor([1.0, ly], dtype=torch.float64)) @ Rc.T for Rc in Rs])
    return exp1, (X3, Y3)

print(" C  delta  delta'   Exp. 1 z / formula      Exp. 3 z / formula")
for C in (1, 2):
    for dl, dp in ((0.5, 0.8), (0.9, 0.3), (0.4, 1.0)):
        (X, U, Y), (X3, Y3) = maximisers(C, dl, dp)
        assert abs(X.michelson_contrast.item() - dl) < 1e-12 and abs(Y.michelson_contrast.item() - dp) < 1e-12
        assert abs(X3.abs().michelson_contrast.item() - dl) < 1e-12 and abs(Y3.michelson_contrast.item() - dp) < 1e-12
        z1v = ((X @ U @ Y).trace.abs() - (X @ Y).trace.real).item()
        z3v = ((Y3 @ X3.abs() @ Y3).trace.real - (Y3 @ X3 @ Y3).trace_norm()).item()
        print(f" {C}   {dl:.1f}    {dp:.1f}    {z1v:.6f} / {env1(C, dl, dp):.6f}    {z3v:.6f} / {env3(C, dl, dp):.6f}")
        assert abs(z1v - env1(C, dl, dp)) < 1e-12 and abs(z3v - env3(C, dl, dp)) < 1e-12
print("limits d -> 1 at delta = 0.4:", env1(1, 0.4, 1.0), 0.4 ** 2 / (2 * 1.4), env3(1, 0.4, 1.0), 1 - 0.4)

# Proposition 11.15 on random complex pairs in M_2: no pair exceeds the bounds
M2c = TypeIAlgebra([2], [2], precision='double')
x = torch.rand(4000, 2, dtype=torch.float64).sort(-1, descending=True)[0]
y = torch.rand(4000, 2, dtype=torch.float64).sort(-1, descending=True)[0]
V, Wc = M2c.random_unitary_operator(4000).matrix, M2c.random_unitary_operator(4000).matrix
Xr = M2c.operator(V @ torch.diag_embed(x.to(V.dtype)).unsqueeze(1) @ V.mH, is_positive=True)
Xi = M2c.operator(V @ torch.diag_embed((x * torch.tensor([1.0, -1.0])).to(V.dtype)).unsqueeze(1) @ V.mH,
                  is_self_adjoint=True)
Yr = M2c.operator(Wc @ torch.diag_embed(y.to(V.dtype)).unsqueeze(1) @ Wc.mH, is_positive=True)
(x1, x2), (y1, y2) = x.unbind(-1), y.unbind(-1)
sup1 = (Yr @ Xr).trace_norm() - (Xr @ Yr).trace.real                        # sup over U (Theorem 1.9)
z3r = (Yr @ Xi.abs() @ Yr).trace.real - M2c.operator((Yr @ Xi @ Yr).matrix, is_self_adjoint=True).trace_norm()
assert (sup1 <= (x1 - x2) ** 2 * (y1 - y2) ** 2 / (4 * (x1 + x2) * (y1 + y2)) + 1e-12).all()
assert (z3r <= 2 * x1 * x2 / (x1 + x2) * (y1 - y2) ** 2 + 1e-12).all()

# multistart search for Exp. 1 on M_2 along Delta(X), Delta(Y) free in [0, 0.999]
torch.manual_seed(2)
M2r = TypeIAlgebra([2], [2], complex_valued=False, precision='double')
gx = torch.tensor([0.3, 0.6, 0.9], dtype=torch.float64)
Xp, Yp1, Up = opt.PositiveParam(M2r, 3 * S, delta=gx.repeat_interleave(S)), opt.PositiveParam(M2r, 3 * S), \
    opt.UnitaryParam(M2r, 3 * S)
def z1m():
    xm, ym = Xp.operator(), Yp1.operator()
    return (xm @ Up.operator() @ ym).trace.abs() - (xm @ ym).trace.real
sup1m = opt.extremize(z1m, [Xp, Yp1, Up], maximize=True, rounds=4, steps=30).reshape(3, S).max(1)[0]
exact1m = env1(1, gx, 0.999)                                                # Corollary 11.16(c)
print("Exp. 1 sup on M_2:", [f"{v:.5f}" for v in sup1m.tolist()], " closed form:", [f"{v:.5f}" for v in exact1m.tolist()])
assert (sup1m <= exact1m + 1e-12).all() and ((exact1m - sup1m) / exact1m).max() < 1e-3
```

<a id="s-11-6"></a>

## 11.6 可与之比较的工具

**Manopt**、**pymanopt** 与 **geoopt** 提供矩阵流形（Stiefel 流形、正交群与酉群），
带有 QR、极分解或 Cayley 等回缩，以及黎曼求解器（共轭梯度、信赖域、黎曼
Adam），可替代在重新居中的坐标卡中运行的 Adam（这里不作性能比较）。**cvxpy** 把算子变量的凸问题
交给 SDP 求解器，后者返回一个最优值以及一个对偶证书，二者都只在求解器容差范围内成立；
松弛——在更大的可行集上的凸问题——从上方界定上确界，这是局部搜索无法给出的方向。**NCTSSOS**、**ncpol2sdpa** 与 **NCSOStools** 在所有维数上同时界定非交换
多项式的特征值和迹（0.3 节）。本库补充的是：$`\bigoplus_cM_{k_c}`$ 的谱约束与
对比度约束的精确参数化（极值可以位于任何扇区）、遵从
块结构的酉算子，以及在 CPU 或 GPU 上对许多起点和约束值进行的一次批量计算。它给出候选的
极值点以及上确界的下界；其余部分必须由证明或松弛来补足。

<a id="x-11"></a>

## 习题

讲义[习题 1 与 2](../../lecture_notes/zh/11_optimisation.md#x-11) 的数值验证：反自伴 $`S`$ 的 Cayley 变换 $`C(S)=(1-S)(1+S)^{-1}`$ 是酉的，没有特征值 $`-1`$，
并且 $`S=(1-C)(1+C)^{-1}`$（习题 1）；对在 $`O(k)`$ 上服从 Haar 分布的 $`O`$（用实代数的 `random_unitary_operator` 抽取，
`complex_valued=False`），$`\mathbb E\mathrm{Tr}(AO)^2=\|A\|_2^2/k`$（习题 2）。

<details>
<summary>解答</summary>

```python
Sk = torch.randn(4, 4, dtype=torch.complex128)
Sk = Sk - Sk.mH                                                   # skew-adjoint
I4 = torch.eye(4, dtype=torch.complex128)
Cay = (I4 - Sk) @ torch.linalg.inv(I4 + Sk)
assert torch.allclose(Cay.mH @ Cay, I4, atol=1e-12)
assert torch.allclose((I4 - Cay) @ torch.linalg.inv(I4 + Cay), Sk, atol=1e-10)
assert torch.linalg.eigvals(Cay).add(1).abs().min() > 1e-6        # -1 is never an eigenvalue
```

```python
Mr = TypeIAlgebra([5], [5], complex_valued=False, precision='double')
Ar = Mr.operator(torch.randn(1, 1, 5, 5, dtype=torch.float64))
tr2 = (Ar @ Mr.random_unitary_operator(20_000)).trace ** 2
assert abs(tr2.mean().item() - Ar.frobenius_norm().item() ** 2 / 5) < 5 * tr2.std().item() / math.sqrt(20_000)
```

</details>

**下一章：** [教程 12](12_equivariant_learning.md)：等变层作为交换子的元素在本库中的实现。
