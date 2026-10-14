[English](../03_traces_and_states.md) · [Español](../es/03_traces_and_states.md) · [Français](../fr/03_traces_and_states.md) · [Deutsch](../de/03_traces_and_states.md) · **中文** · [日本語](../ja/03_traces_and_states.md) · [Русский](../ru/03_traces_and_states.md)

# 教程 3. 迹、态与熵

本教程配合讲义的[第 3 章](../../lecture_notes/zh/03_traces_and_states.md)。该章对
$\mathcal M=\bigoplus_cM_{k_c}$ 上的迹进行分类（[定理 3.1](../../lecture_notes/zh/03_traces_and_states.md#r-3-1)），用密度表示态
（[定理 3.2](../../lecture_notes/zh/03_traces_and_states.md#r-3-2)），并证明基本的熵不等式。我们使用本库的三种迹，
在它们之间换算密度，并在 $M_2\oplus M_3\oplus M_5$ 上数值地检验这些定理。

<a id="s-3-1"></a>

## 3.1 迹

*讲义：* [3.1 节](../../lecture_notes/zh/03_traces_and_states.md#s-3-1)；由[定理 3.1](../../lecture_notes/zh/03_traces_and_states.md#r-3-1)，每个迹都是 $\mathrm{tr}_w(A)=\sum_cw_c\mathrm{Tr}\,A_c$。

本库实现了三种忠实迹：

| 名称 | 权重 $w_c$ | 在 $1$ 上的值 | 含义 |
|---|---|---|---|
| `Tr_blunt` | $1$ | $D$ | $\bigoplus_c\mathbb C^{k_c}$ 的迹（重数 $m_c=1$） |
| `Tr_norm` | $1/k_c$ | $C$ | 每个块都归一化 |
| `tau_vN` | $1/(Ck_c)$ | $1$ | 各扇区权重均为 $1/C$ 的迹态 |

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, DensityMatrix, states, channels

torch.manual_seed(0)
alg = TypeIAlgebra([2, 3, 5], [2, 3, 5], precision='double')    # M_2 (+) M_3 (+) M_5, D = 10, C = 3
U = alg.random_unitary_operator(4)
A = alg.operator_from_eigenvalues(lambda d: torch.randn(4, d), batch_size=4, force_self_adjoint=True)
B = U @ A                                                       # a batch of four non-normal operators
one = alg.identity()

for f, value in ((lambda X: X.Tr_blunt(), 10), (lambda X: X.Tr_norm(), 3), (lambda X: X.tau_vN(), 1)):
    assert abs(f(one).item() - value) < 1e-12                   # normalisation
    assert (f(A @ B) - f(B @ A)).abs().max() < 1e-10            # trace property
assert (A.Tr_norm() - 3 * A.tau_vN()).abs().max() < 1e-12       # Tr_norm = C tau_vN
```

<a id="s-3-2"></a>

## 3.2 态与密度

*讲义：* [3.2 节](../../lecture_notes/zh/03_traces_and_states.md#s-3-2)；[定理 3.2](../../lecture_notes/zh/03_traces_and_states.md#r-3-2) 对每个忠实迹给出态的密度 $\rho_w$，且 $w_c(\rho_w)_c=w'_c(\rho_{w'})_c$。

本库通过态关于 `Tr_blunt` 的密度 $\rho=\rho_{\mathrm{blunt}}$ 来存储态
（类 `DensityMatrix`）；`density(trace)` 和 `DensityMatrix.from_density`
负责换算。态的取值与约定无关：

```python
rho = states.random_density_matrix(alg, batch_size=4)          # Hilbert-Schmidt random states
for trace, tr in (('blunt', lambda X: X.Tr_blunt()), ('norm', lambda X: X.Tr_norm()),
                  ('tau_vN', lambda X: X.tau_vN())):
    d = alg.operator(rho.density(trace))
    assert (tr(d @ A) - (rho @ A).Tr_blunt()).abs().max() < 1e-10       # omega(A) for all traces
    assert (tr(d) - 1).abs().max() < 1e-10                               # tr_w(rho_w) = 1
    back = DensityMatrix.from_density(alg, rho.density(trace), trace=trace)
    assert (back.matrix - rho.matrix).abs().max() < 1e-12
```

<a id="s-3-3"></a>

## 3.3 态的经典部分

*讲义：* [3.3 节](../../lecture_notes/zh/03_traces_and_states.md#s-3-3)：扇区概率 $p_c=\mathrm{Tr}\,\rho_c$，以及在 $\mathcal M$ 上为恒等映射的映射 $\rho\mapsto\sum_cP_c\rho P_c$。

```python
p = rho.sector_probabilities()                                   # shape (batch, C)
assert torch.allclose(p.sum(-1), torch.ones(4, dtype=p.dtype))
P = [alg.central([1.0 if c == d else 0.0 for d in range(alg.C)]) for c in range(alg.C)]
measure_sector = channels.measurement_channel(P)
assert (measure_sector(rho).matrix - rho.matrix).abs().max() < 1e-12   # no disturbance
```

<a id="s-3-4"></a>

## 3.4 熵

*讲义：* [3.4 节](../../lecture_notes/zh/03_traces_and_states.md#s-3-4)；[定理 3.3](../../lecture_notes/zh/03_traces_and_states.md#r-3-3) 给出 $S(\rho)=H(p)+\sum_cp_cS(\sigma_c)$，而相对于 $\mathrm{tr}_w$ 的熵为 $S_w(\rho)=S(\rho)+\sum_cp_c\log w_c$。

```python
S = rho.von_neumann_entropy()
sigma = [rho.matrix[:, c, :k, :k] / p[:, c, None, None] for c, k in enumerate(alg.k_factors)]
S_sec = [-(w * torch.log(w)).sum(-1) for w in (torch.linalg.eigvalsh(s) for s in sigma)]
H = -(p * torch.log(p)).sum(-1)
assert torch.allclose(S, H + sum(p[:, c] * S_sec[c] for c in range(alg.C)))    # Theorem 3.3

w_tau = torch.tensor([1 / (alg.C * k) for k in alg.k_factors], dtype=torch.float64)
d_tau = alg.operator(rho.density('tau_vN'), is_positive=True)
S_tau = -(d_tau @ d_tau.log()).tau_vN().real
assert torch.allclose(S_tau, S + (p * torch.log(w_tau)).sum(-1))
```

<a id="s-3-5"></a>

## 3.5 相对熵与 Klein 不等式

*讲义：* [3.5 节](../../lecture_notes/zh/03_traces_and_states.md#s-3-5)；[推论 3.6](../../lecture_notes/zh/03_traces_and_states.md#r-3-6)（$S(\rho)\le\log D$，当 $\rho=1/D$ 时取等号）与[定理 3.7](../../lecture_notes/zh/03_traces_and_states.md#r-3-7)（对到中心上的条件期望 $E$，$S(\rho\|\sigma)\ge S(E\rho\|E\sigma)=\sum_cp_c\log(p_c/q_c)$）。

```python
sigma_state = states.random_density_matrix(alg, batch_size=4)
q = sigma_state.sector_probabilities()
D_rel = states.relative_entropy(rho, sigma_state)
kl = (p * torch.log(p / q)).sum(-1)
E = channels.center_expectation(alg)
assert (D_rel >= kl - 1e-12).all()                                       # Theorem 3.7
assert torch.allclose(states.relative_entropy(E(rho), E(sigma_state)), kl)
assert (rho.von_neumann_entropy() <= math.log(10) + 1e-12).all()         # Corollary 3.6
assert abs(states.maximally_mixed_state(alg).von_neumann_entropy().item() - math.log(10)) < 1e-12
```

<a id="s-3-6"></a>

## 3.6 Michelson 对比度

*讲义：* [3.6 节](../../lecture_notes/zh/03_traces_and_states.md#s-3-6)；由[命题 3.8](../../lecture_notes/zh/03_traces_and_states.md#r-3-8)，正的保单位映射不会增大 Michelson 对比度；本库以 `Operator.michelson_contrast` 计算它。

```python
A_pos = alg.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(4, d), batch_size=4, force_positive=True)
Phi = channels.random_mixed_unitary_channel(alg, 3, batch_size=4)       # unital
assert (Phi(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
assert (E(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
```

<a id="s-3-7"></a>

## 3.7 一个看似合理却不成立的不等式

*讲义：* [3.7 节](../../lecture_notes/zh/03_traces_and_states.md#s-3-7)；由[命题 3.9](../../lecture_notes/zh/03_traces_and_states.md#r-3-9)，对 $M_k$ 中的 $X,Y\ge0$，$\max_U|\mathrm{Tr}(XUY)|=\|YX\|_1\ge\mathrm{Tr}(XY)$，等号成立当且仅当 $XY=YX$。

本库借助极分解求出最大值点：对 $YX=W|YX|$，酉矩阵 $U=W^*$ 达到 $\|YX\|_1$
（[讲义，定理 1.9](../../lecture_notes/zh/01_matrix_algebras.md#r-1-9)）。在按 Haar 分布的酉矩阵上随机搜索
（本库正是这样对 $U$ 采样的）平均而言更差
（[讲义，3.7 节](../../lecture_notes/zh/03_traces_and_states.md#s-3-7)）；那里回顾的概率上界 $b^4/(a^4k^2)$ 在下面的运行中
什么也没说，因为 $b/a=6$ 且 $k=6$。代码对一对 $k=6$ 的不可交换矩阵计算精确的上确界以及 $20000$ 个 Haar
样本中的最佳值，二者都减去 $\mathrm{Tr}(XY)$。由命题 3.9，没有样本能超过上确界；最佳样本低于上确界
多少，是本次运行的*数值观察*，而不是定理。

```python
Mk = TypeIAlgebra([6], [6], precision='double')
X = Mk.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
Y = Mk.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
assert (X @ Y - Y @ X).matrix.abs().max().item() > 1e-6             # X and Y do not commute
YX = Y @ X
W = YX @ YX.abs().inverse()                                          # polar part of YX
gap = (X @ W.adjoint() @ Y).trace.abs() - (X @ Y).trace.real
assert abs(gap.item() - (YX.trace_norm() - (X @ Y).trace.real).item()) < 1e-10   # Theorem 1.9
assert gap.item() > 0                                  # Proposition 3.9: the inequality fails
samples = Mk.random_unitary_operator(20_000)
best = (X @ samples @ Y).trace.abs().max() - (X @ Y).trace.real
print(f"sup - Tr(XY) = {gap.item():.4f}; numerical observation: best of 20000 Haar samples "
      f"- Tr(XY) = {best.item():.4f}")
assert best.item() <= gap.item() + 1e-10               # Proposition 3.9: no sample beats the sup
```

<a id="s-3-8"></a>

## 3.8 可与之比较的工具

QuTiP（`entropy_vn`、`entropy_relative`）以及其他量子信息软件包（例如
`toqito`）用单个 Hilbert 空间的迹计算该空间上密度矩阵的熵与相对熵；在
$\bigoplus_c\mathbb C^{k_c}$ 上，这就是这里的 blunt 迹（至多相差对数的底，不同软件包的取法不同）。它们不提供的是定理 3.1
中的迹族以及定理 3.2 中的换算，而一旦代数以一定重数作用，或者各扇区必须赋予不同权重，
这些就变得必不可少。

<a id="x-3"></a>

## 习题

讲义的[习题 1](../../lecture_notes/zh/03_traces_and_states.md#x-3)：迹态 $\tau_{\mathrm{vN}}$
（`states.tracial_state`）关于 $\tau_{\mathrm{vN}}$ 的密度为 $1$，其熵为
$S=\log C+\frac1C\sum_c\log k_c$。对 $M_2\oplus M_3\oplus M_5$ 检验这两点。

<details>
<summary>解答</summary>

```python
tau = states.tracial_state(alg)
assert (tau.density('tau_vN') - one.matrix).abs().max() < 1e-12
assert abs(tau.von_neumann_entropy().item() - (math.log(3) + sum(math.log(k) for k in (2, 3, 5)) / 3)) < 1e-12
```

</details>

**下一章：** [教程 4. 对称性与表示](04_symmetries.md)
