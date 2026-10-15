[English](../08_dynamics.md) · [Español](../es/08_dynamics.md) · [Français](../fr/08_dynamics.md) · [Deutsch](../de/08_dynamics.md) · **中文** · [日本語](../ja/08_dynamics.md) · [Русский](../ru/08_dynamics.md)

# 教程 8. 动力学

本教程配合讲义的[第 8 章](../../lecture_notes/zh/08_dynamics.md)。该章证明 $`\mathcal M`$ 中的哈密顿量保持扇区权重不变
（[定理 8.3](../../lecture_notes/zh/08_dynamics.md#r-8-3)），GKSL 形式的生成元产生量子动力学半群（[定理 8.6](../../lecture_notes/zh/08_dynamics.md#r-8-6)），
并讨论扇区间的跳跃、后选择与稳态。我们用 `dynamics` 的精确求解器和 Runge–Kutta 求解器求解 Schrödinger、
von Neumann 和 Lindblad 方程，并将结果与 SciPy 以及 $`M_D`$ 上的 GKSL 生成元对照。

<a id="s-8-1"></a>

## 8.1 封闭系统

*讲义：* [8.1 节](../../lecture_notes/zh/08_dynamics.md#s-8-1)；[命题 8.2](../../lecture_notes/zh/08_dynamics.md#r-8-2)（传播子 $`U(t)=e^{-itH}`$）、[定理 8.3](../../lecture_notes/zh/08_dynamics.md#r-8-3)（荷守恒）。

由命题 8.2(1)，演化分解为 $`C`$ 个大小为 $`k_c`$ 的独立问题，每个块一个。
本库对每种块大小做一次特征分解来计算 $`U(t)`$（`dynamics.propagator`），而
`dynamics.von_neumann` 在本征基中一次性把密度演化到许多时刻，在该基中
$`\rho(t)_{ij}=e^{-it(\lambda_i-\lambda_j)}\rho_{ij}`$。我们用 `scipy.linalg.expm` 加以检验：

```python
import math
import numpy as np
import scipy.integrate
import scipy.linalg as sla
import torch
from torch_vn_algebra import TypeIAlgebra, DensityMatrix, InterSectorChannel, channels, dynamics, states

torch.manual_seed(0)
f64, c128 = torch.float64, torch.complex128
alg = TypeIAlgebra([2, 3, 4], [2, 3, 4], precision='double')     # M_2 (+) M_3 (+) M_4, D = 9
H = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
t = 0.7
U = dynamics.propagator(H, t)                                      # exp(-itH), block by block
for c, k in enumerate(alg.k_factors):                              # independent reference
    ref = sla.expm(-1j * t * H.matrix[0, c, :k, :k].numpy())
    assert np.abs(U.matrix[0, c, :k, :k].numpy() - ref).max() < 1e-12
assert (U.adjoint() @ U - alg.identity()).matrix.abs().max() < 1e-12

rho0 = states.random_density_matrix(alg)
times = [0.0, 0.5, 1.0, 2.0]
traj = dynamics.von_neumann(rho0, H, times)                        # exact, all times at once
for s, rho_s in zip(times, traj):
    U_s = dynamics.propagator(H, s)
    assert (rho_s.matrix - (U_s @ rho0 @ U_s.adjoint()).matrix).abs().max() < 1e-12
    assert (rho_s.sector_probabilities() - rho0.sector_probabilities()).abs().max() < 1e-12   # Theorem 8.3
A = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
A_t = U.adjoint() @ A @ U                                          # Heisenberg picture
assert abs(rho0.expectation(A_t) - dynamics.von_neumann(rho0, H, [t])[0].expectation(A)) < 1e-12
```

单位向量 $`\psi=\bigoplus_c\psi_c`$ 在 $`\mathcal M`$ 上定义态 $`\bigoplus_c|\psi_c\rangle\langle\psi_c|`$，它按 $`U(t)`$ 演化
（[讲义，8.1 节](../../lecture_notes/zh/08_dynamics.md#s-8-1)），所以 `dynamics.schrodinger`（向量采用块布局 `(batch, C, k_max)`）与 `von_neumann`
在它上面必须一致。对含时的 $`H(t)`$，
本库使用固定步长的经典 Runge–Kutta 方法（`dynamics.rk4`、`schrodinger_rk4`）：每一步把
$`\psi`$ 乘以 $`H(t)\in\mathcal M`$ 若干取值的一个多项式，所以它从不在扇区之间转移振幅，
但它不是酉的：对常数 $`H`$，步长为 $`h`$ 的一步把 $`H`$ 的属于特征值 $`\lambda`$ 的特征向量乘以
$`T(-ih\lambda)`$，其中 $`T(z)=\sum_{k=0}^4z^k/k!`$，而对 $`0<|x|<2\sqrt2`$ 有
$`|T(ix)|^2=1-\frac{x^6}{72}+\frac{x^8}{576}<1`$（展开即得）。所以 RK4 并不严格守恒扇区权重：对常数 $`H`$，
一个在某个满足 $`0<|h\lambda|<2\sqrt2`$ 的特征向量上系数非零的扇区分量在每一步都会损失权重。
它们在一次运行中变化多少是一个数值观察（见下文）。
（对于一个大的扇区，`krylov.evolve` 只需要矩阵–向量乘积；见[教程 9，9.5 节](../../tutorial/zh/09_spin_chains.md#s-9-5)。）

```python
rows = torch.tensor([[1, 1, 0, 0], [1, 1, 1, 0], [1, 1, 1, 1]], dtype=c128)   # entries of C^2 (+) C^3 (+) C^4
psi0 = torch.randn(1, alg.C, alg.k_max, dtype=c128) * rows
psi0 = psi0 / psi0.norm()
psi = dynamics.schrodinger(psi0, H, times)                         # shape (T, batch, C, k_max)
p_vec = (psi.abs() ** 2).sum(-1)                                   # ||P_c psi(t)||^2
assert (p_vec - p_vec[0]).abs().max() < 1e-12
mix = DensityMatrix(alg, matrix=psi0[..., :, None] * psi0[..., None, :].conj())   # the state of psi0 on M
outer = psi[..., :, None] * psi[..., None, :].conj()
assert all((o - r.matrix).abs().max() < 1e-12 for o, r in zip(outer, dynamics.von_neumann(mix, H, times)))

V = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
H_of_t = lambda s: H.matrix + math.cos(3 * s) * V.matrix          # driven, in M at every time
psi_rk = dynamics.schrodinger_rk4(psi0, H_of_t, times, substeps=50, algebra=alg)
drift = ((psi_rk.abs() ** 2).sum(-1) - p_vec[0]).abs().max().item()
print(f"RK4, time-dependent H: largest change of a sector weight {drift:.1e}")
assert drift < 1e-8 and (dynamics.schrodinger_rk4(psi0, H, times, substeps=50) - psi).abs().max() < 1e-7

```

*数值观察*：在这次运行中，扇区权重在 RK4 下的变化为 $`1.1\cdot10^{-9}`$（打印出的值），而在精确
求解器下的变化小于断言所用的容差 $`10^{-12}`$；对精确求解器，定理 8.3 给出严格的守恒。

<a id="s-8-2"></a>

## 8.2 量子动力学半群

*讲义：* [8.2 节](../../lecture_notes/zh/08_dynamics.md#s-8-2)；[定理 8.6](../../lecture_notes/zh/08_dynamics.md#r-8-6)（GKSL 生成元与 Dyson 展开）、[定理 8.7](../../lecture_notes/zh/08_dynamics.md#r-8-7)。

本库有四个 Lindblad 例程：`lindblad_rhs` 以函数形式返回右端项，供
`rk4` 使用；`lindblad_evolve` 用 RK4 积分；`lindblad_superoperator` 对每个扇区返回 $`\mathcal L`$ 作用在
按行展开的 $`\mathrm{vec}\,\rho_c`$ 上的 $`k_c^2\times k_c^2`$ 矩阵；
`lindblad_channel` 逐扇区精确计算 $`e^{t\mathcal L}`$，结果为一个 `Channel`。我们对
$`M_3`$ 上任意的跳跃检验定理 8.6：$`e^{t\mathcal L}`$ 的 Choi 矩阵是正的（定理
6.2），该映射是 TP 的，半群律成立，RK4 与指数一致，并且跳跃部分
$`e^{t\mathcal L}-S_0(t)=\sum_{n\ge1}S_n(t)`$ 本身是 CP 的：

```python
M3 = TypeIAlgebra([3], [3], precision='double')
H3 = M3.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
jumps = [M3.operator(torch.randn(1, 1, 3, 3, dtype=c128) / 2) for _ in range(2)]     # arbitrary L_j
T = lambda s: dynamics.lindblad_channel(M3, H3, jumps, t=s)        # exp(sL), exact
S = lambda s: T(s).superoperator()[0][0]                           # its 9 x 9 matrix on row-major vec
choi = lambda X: X.reshape(3, 3, 3, 3).permute(2, 0, 3, 1).reshape(9, 9)   # J_{(i,a),(j,b)} = S_{(a,b),(i,j)}
assert (choi(S(0.8)) - T(0.8).choi()[0][0]).abs().max() < 1e-12
assert torch.linalg.eigvalsh(choi(S(0.8))).min() > -1e-12         # completely positive
assert T(0.8).trace_preservation_error().max() < 1e-12             # trace preserving
assert (S(0.3) @ S(0.5) - S(0.8)).abs().max() < 1e-12              # semigroup

rho3 = states.random_density_matrix(M3)
rk = dynamics.lindblad_evolve(rho3, H3, jumps, [0.0, 0.8], substeps=40)[-1]
assert (rk.matrix - T(0.8)(rho3).matrix).abs().max() < 1e-7        # RK4 against exp(tL)
assert abs(dynamics.lindblad_rhs(H3, jumps)(0.0, rho3.matrix)[0, 0].trace()) < 1e-12   # Tr L(rho) = 0

H_eff = H3.matrix - 0.5j * sum(L.adjoint().matrix @ L.matrix for L in jumps)
S0 = channels.Channel(M3, torch.linalg.matrix_exp(-0.8j * H_eff)[:, None]).superoperator()[0][0]
jump_part = torch.linalg.eigvalsh(choi(S(0.8) - S0))
print(f"smallest eigenvalue of the Choi matrix of the jump part: {jump_part.min().item():.1e}")
assert jump_part.min() > -1e-12                                    # sum_{n >= 1} S_n(t) is CP
```

<a id="s-8-3"></a>

## 8.3 有扇区的代数中的 Lindblad 方程

*讲义：* [8.3 节](../../lecture_notes/zh/08_dynamics.md#s-8-3)：扇区内部与扇区之间都有跳跃的生成元 $`\mathcal L`$，由 `dynamics.lindblad_rhs` 实现；[命题 8.8](../../lecture_notes/zh/08_dynamics.md#r-8-8)、[命题 8.9](../../lecture_notes/zh/08_dynamics.md#r-8-9)（扇区权重）。

**例：粒子丢失**（[讲义，8.3 节](../../lecture_notes/zh/08_dynamics.md#s-8-3)）。两个玻色模式，至多两个粒子，
$`\mathcal M=M_1\oplus M_2\oplus M_3`$（扇区 $`N=0,1,2`$），Bose–Hubbard 二聚体
$`H=J(a_1^*a_2+a_2^*a_1)+\frac U2\sum_jn_j(n_j-1)`$，并且每个模式以速率 $`\gamma`$ 丢失粒子：这是一个
`InterSectorChannel`，其 Kraus 算子为 $`a_1,a_2`$，限制在 $`N\to N-1`$ 上。由
[命题 8.9(2)](../../lecture_notes/zh/08_dynamics.md#r-8-9)，扇区权重服从 $`p'=Qp`$，其中对速率 $`1`$ 有 $`Q=R-\mathrm{diag}(\sum_dR_{dc})_c`$，
$`R`$ = `transition_matrix()`；从 $`N=2`$ 出发，它们服从二项分布 $`p_2=q^2`$，$`p_1=2q(1-q)`$，$`p_0=(1-q)^2`$，
$`q=e^{-\gamma t}`$，无论隧穿与相互作用如何。代码检验
$`\Phi^*(1)=0\oplus1\oplus2`$，即在所用的块上 $`a_1^*a_1+a_2^*a_2=N`$，并用两个公式检验 $`p(t)`$。

```python
fock = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double', charges=[0, 1, 2])   # N = 0, 1, 2
r2 = math.sqrt(2)
a1 = {(0, 1): [[1.0, 0.0]], (1, 2): [[r2, 0.0, 0.0], [0.0, 1.0, 0.0]]}           # a_1 : N -> N - 1
a2 = {(0, 1): [[0.0, 1.0]], (1, 2): [[0.0, 1.0, 0.0], [0.0, 0.0, r2]]}           # a_2 : N -> N - 1
blocks = lambda a: {dc: torch.tensor(K, dtype=f64) for dc, K in a.items()}
loss = InterSectorChannel.from_blocks(fock, fock, {dc: [blocks(a1)[dc], blocks(a2)[dc]] for dc in a1})
J, U_int = 1.0, 0.8                                                # hopping and interaction
H_bh = fock.from_blocks([torch.zeros(1, 1, dtype=f64), torch.tensor([[0.0, J], [J, 0.0]], dtype=f64),
                         torch.tensor([[U_int, r2 * J, 0.0], [r2 * J, 0.0, r2 * J], [0.0, r2 * J, U_int]], dtype=f64)])

assert (loss.adjoint()(fock.identity()).matrix - fock.central([0.0, 1.0, 2.0]).matrix).abs().max() < 1e-12
R = loss.transition_matrix()[0]                                    # r_dc (columns: from sector c)
Q = R - torch.diag(R.sum(0))                                       # Proposition 8.9(2), for rate 1
print(Q)

gamma = 0.5
psi2 = torch.randn(3, dtype=c128)
psi2 = psi2 / psi2.norm()
m0 = torch.zeros(1, 3, 3, 3, dtype=c128)
m0[0, 2] = torch.outer(psi2, psi2.conj())                          # a pure state with N = 2
ts = [0.0, 0.5, 1.0, 2.0, 4.0]
traj = dynamics.lindblad_evolve(DensityMatrix(fock, matrix=m0), H_bh, [loss], ts, rates=[gamma], substeps=100)
for s, rho_s in zip(ts, traj):
    q = math.exp(-gamma * s)
    p = rho_s.sector_probabilities()[0]
    assert (p - torch.tensor([(1 - q) ** 2, 2 * q * (1 - q), q ** 2], dtype=f64)).abs().max() < 1e-8
    assert (p - torch.from_numpy(sla.expm(gamma * s * Q.numpy())[:, 2])).abs().max() < 1e-8
neg = min(w.min().item() for rho_s in traj for w in rho_s.eigenvalues())
print(f"smallest eigenvalue of rho(t) under RK4: {neg:.1e}")
assert neg > -1e-7                                                 # positive up to the RK4 error

sm = torch.tensor([[0.0, 1.0], [0.0, 0.0]], dtype=c128)           # |g><e|; H = 0, rate 3, one step h = 1
decay = lambda s, X: 3.0 * (sm @ X @ sm.mH - 0.5 * (sm.mH @ sm @ X + X @ sm.mH @ sm))
step = dynamics.rk4(decay, torch.diag(torch.tensor([0.0, 1.0], dtype=c128)), [0.0, 1.0], substeps=1)[-1]
assert (step - torch.diag(torch.tensor([-3 / 8, 11 / 8], dtype=c128))).abs().max() < 1e-14   # not positive
```

RK4 不能保证保持正性：对 $`\rho'=\mathcal L\rho`$，步长为 $`h`$ 的一步是
$`\rho\mapsto T(h\mathcal L)\rho`$，其中 $`T`$ 是 8.1 节的多项式（代入四个阶段即得）；它是
保迹的（因为对所有 $`X`$ 有 $`\mathrm{Tr}\,\mathcal L(X)=0`$），但不一定是正的。例如，在以
$`|g\rangle,|e\rangle`$ 为基的 $`M_2`$ 上，取 $`H=0`$ 和一个速率为 $`\gamma`$ 的跳跃 $`|g\rangle\langle e|`$，则 $`\mathcal L`$ 把
$`\mathrm{diag}(a,b)`$ 映为 $`\gamma b\mathrm{diag}(1,-1)`$，所以 $`h\gamma=3`$ 的一步把 $`\mathrm{diag}(0,1)`$ 映为
$`\mathrm{diag}(1-T(-3),T(-3))=\mathrm{diag}(-\frac38,\frac{11}8)`$。*数值观察*：在这次运行中，
$`\rho(t)`$ 的最小特征值为 $`\approx-2\cdot10^{-8}`$（打印出的值，步长至多 $`0.02`$）。`lindblad_channel` 是精确的，
但只处理扇区内部的跳跃。最后，我们对扇区之间的一个随机映射检验命题 8.8(2)，
并对粒子丢失检验命题 8.8(3)，与用 numpy 构造的 $`M_6`$ 上的 GKSL 生成元进行比较：

```python
Phi = channels.random_inter_sector_channel(fock, fock, kraus_rank=2)
rho = states.random_density_matrix(fock)
on_M = dynamics.lindblad_rhs(None, [Phi], algebra=fock)(0.0, rho.matrix)[0].numpy()
ks, K = fock.k_factors, Phi.kraus[0].numpy()                       # K[l, d, c] = K_l^{dc}, zero padded
embed = lambda Y: sla.block_diag(*[Y[c, :k, :k] for c, k in enumerate(ks)])    # (+)_c Y_c in M_6
Ls = [np.block([[K[l, d, c, :kd, :kc] for c, kc in enumerate(ks)] for d, kd in enumerate(ks)])
      for l in range(2)]                                           # L_l = sum_{d,c} J_l^{dc}
X = embed(rho.matrix[0].numpy())
full = sum(L @ X @ L.conj().T - 0.5 * (L.conj().T @ L @ X + X @ L.conj().T @ L) for L in Ls)   # GKSL on M_6
E_full = full * embed(np.ones((3, 3, 3)))                          # E(X) = sum_d P_d X P_d
assert np.abs(E_full - embed(on_M)).max() < 1e-12                  # Proposition 8.8(2)
print(f"coherences between sectors created on M_6: {np.abs(full - E_full).max():.2f}")

zero = lambda d, c: np.zeros((ks[d], ks[c]))
L_loss = [np.block([[blocks(a)[(d, c)].numpy() if (d, c) in a else zero(d, c) for c in range(3)] for d in range(3)])
          for a in (a1, a2)]                                       # a_1, a_2 on M_6: one block per block row and column
full_loss = sum(L @ X @ L.T - 0.5 * (L.T @ L @ X + X @ L.T @ L) for L in L_loss)
loss_on_M = dynamics.lindblad_rhs(None, [loss], algebra=fock)(0.0, rho.matrix)[0].numpy()
assert np.abs(full_loss - embed(loss_on_M)).max() < 1e-12         # Proposition 8.8(3): no compression needed
```

<a id="s-8-4"></a>

## 8.4 后选择与无跳跃演化

*讲义：* [8.4 节](../../lecture_notes/zh/08_dynamics.md#s-8-4)；[命题 8.10](../../lecture_notes/zh/08_dynamics.md#r-8-10)（无跳跃概率 $`w_0`$）、[命题 8.11](../../lecture_notes/zh/08_dynamics.md#r-8-11)（对"没有粒子丢失"的后选择）。

本库的 `states.lueders_update(rho, P)` 返回 $`P\rho P/\mathrm{Tr}(P\rho P)`$ 和
$`\mathrm{Tr}(\rho P)`$。现在让两个模式以不同的速率 $`\gamma_1>\gamma_2`$ 泄漏；则
$`A=\gamma_1n_1+\gamma_2n_2`$ 不是中心元（在扇区 $`N=2`$ 中，它在基 $`|20\rangle,|11\rangle,|02\rangle`$ 下是
$`\mathrm{diag}(2\gamma_1,\gamma_1+\gamma_2,2\gamma_2)`$），命题 8.9(2) 不适用，而 $`N=2`$ 上后选择的态是
由 $`H_{\mathrm{eff}}`$ 给出的归一化演化，其反厄米部分 $`-\frac i2A`$ 在该扇区中不是 $`1`$ 的倍数。我们对 $`n=0`$ 检验命题 8.11，并通过用 Simpson 公式计算
$`S_1`$ 的积分，对 $`n=1`$ 检验它：

```python
g1, g2 = 0.6, 0.2                                                  # mode 1 leaks faster
loss1 = InterSectorChannel.from_blocks(fock, fock, {dc: [K] for dc, K in blocks(a1).items()})
loss2 = InterSectorChannel.from_blocks(fock, fock, {dc: [K] for dc, K in blocks(a2).items()})
m1 = torch.zeros(1, 3, 3, 3, dtype=c128)
m1[0, 2, 0, 0] = 1.0                                               # |20>: both particles in mode 1
ts = [0.0, 0.5, 1.0, 2.0]
traj = dynamics.lindblad_evolve(DensityMatrix(fock, matrix=m1), H_bh, [loss1, loss2], ts, rates=[g1, g2],
                                substeps=100)

n_1 = [np.diag(x) for x in ([0.0], [1.0, 0.0], [2.0, 1.0, 0.0])]   # n_1 in each sector; n_2 = N - n_1
H_eff = [H_bh.matrix[0, N, :N + 1, :N + 1].numpy() - 0.5j * (g1 * n + g2 * (N * np.eye(N + 1) - n))
         for N, n in enumerate(n_1)]

def S0(s, N, X):                                  # e^{-isH_eff} X e^{isH_eff^*} in sector N, s a time or an array
    E = torch.linalg.matrix_exp(torch.from_numpy(-1j * np.asarray(s)[..., None, None] * H_eff[N])).numpy()
    return E @ X @ E.conj().swapaxes(-1, -2)

rho_20 = m1[0, 2].numpy()
for s, rho_s in zip(ts, traj):
    w0 = np.trace(S0(s, 2, rho_20)).real                           # probability of no loss
    post, prob = states.lueders_update(rho_s, fock.central([0.0, 0.0, 1.0]))   # post-select N = 2
    assert abs(prob.item() - w0) < 1e-8                            # Proposition 8.11, n = 0
    assert np.abs(post.matrix[0, 2].numpy() - S0(s, 2, rho_20) / w0).max() < 1e-7
closed = dynamics.von_neumann(DensityMatrix(fock, matrix=m1), H_bh, [2.0])[0]   # no losses at all
dist = 0.5 * np.abs(np.linalg.eigvalsh(S0(2.0, 2, rho_20) / w0 - closed.matrix[0, 2].numpy())).sum()
print(f"t = 2: P(no loss) = {w0:.3f} (without hopping: {math.exp(-2 * g1 * 2):.3f}); "
      f"trace distance post-selected/closed system: {dist:.3f}")

K1, K2 = blocks(a1)[(1, 2)].numpy(), blocks(a2)[(1, 2)].numpy()   # S_1(2) rho = int_0^2 S_0(2-s) Phi S_0(s) rho ds
jump = lambda X: g1 * K1 @ X @ K1.T + g2 * K2 @ X @ K2.T
s_grid = np.linspace(0.0, 2.0, 201)
S1 = scipy.integrate.simpson(S0(2.0 - s_grid, 1, jump(S0(s_grid, 2, rho_20))), x=s_grid, axis=0)
assert np.abs(S1 - traj[-1].matrix[0, 1, :2, :2].numpy()).max() < 1e-8     # Proposition 8.11, n = 1
```

*数值观察*：在 $`t=2`$ 时，扇区 $`N=2`$ 的权重 $`p_2(2)=w_0(2)`$（"没有损失"）为 $`0.230`$（见上面的
打印输出）。若没有隧穿（$`J=0`$），它将恰好是 $`e^{-2\gamma_1t}=e^{-2.4}=0.091`$，因为这时 $`|20\rangle`$ 是
$`H_{\mathrm{eff}}`$ 的属于特征值 $`U-i\gamma_1`$ 的特征向量，且 $`w_0(t)=|e^{-it(U-i\gamma_1)}|^2`$（命题 8.10）；所以
这里扇区权重确实依赖于哈密顿量。后选择的态与封闭系统的态之间的迹距离为
$`0.135`$（见上面的打印输出）：以 $`N=2`$ 为条件会改变态，这种改变超出了酉演化。

<a id="s-8-5"></a>

## 8.5 稳态与弛豫

*讲义：* [8.5 节](../../lecture_notes/zh/08_dynamics.md#s-8-5)；[推论 8.13](../../lecture_notes/zh/08_dynamics.md#r-8-13)、[命题 8.14](../../lecture_notes/zh/08_dynamics.md#r-8-14)（自发辐射）。

由[命题 8.14](../../lecture_notes/zh/08_dynamics.md#r-8-14)，$`M_2`$ 上的自发辐射（$`H=\omega|e\rangle\langle e|`$，一个速率为 $`\gamma`$ 的跳跃
$`|g\rangle\langle e|`$）满足 $`e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}}`$，其谱为
$`\{0,-\gamma,-\gamma/2\pm i\omega\}`$，谱隙为 $`\gamma/2`$。我们用
`lindblad_superoperator` 检验谱，用 `lindblad_channel` 检验信道与迹距离，并对
$`M_2\oplus M_3`$ 扇区内部的随机跳跃检验[推论 8.13](../../lecture_notes/zh/08_dynamics.md#r-8-13)：

```python
q2 = TypeIAlgebra([2], [2], precision='double')                    # basis |g>, |e>
omega, gamma = 1.3, 0.4
H_at = q2.from_blocks([torch.tensor([[0.0, 0.0], [0.0, omega]], dtype=f64)])
sigma_minus = q2.from_blocks([torch.tensor([[0.0, 1.0], [0.0, 0.0]], dtype=f64)])   # |g><e|
gen = dynamics.lindblad_superoperator(q2, H_at, [sigma_minus], rates=[gamma])[0][0]
spec = sorted(torch.linalg.eigvals(gen).tolist(), key=lambda z: (z.real, z.imag))
print([f"{z.real:+.2f}{z.imag:+.2f}i" for z in spec])
assert np.allclose(spec, [-gamma, -gamma / 2 - 1j * omega, -gamma / 2 + 1j * omega, 0.0], atol=1e-12)

psi_at = torch.tensor([0.6, 0.8], dtype=c128)                      # rho_ee(0) = 0.64, rho_ge(0) = 0.48
rho_at = DensityMatrix(q2, matrix=torch.outer(psi_at, psi_at.conj())[None, None])
ground = DensityMatrix(q2, matrix=torch.tensor([[[[1.0, 0.0], [0.0, 0.0]]]], dtype=c128))
for s in (0.5, 2.0, 10.0):
    ch = dynamics.lindblad_channel(q2, H_at, [sigma_minus], t=s, rates=[gamma])
    UA = channels.unitary_channel(dynamics.propagator(H_at, s)).compose(
        channels.amplitude_damping_channel(q2, 1 - math.exp(-gamma * s)))
    assert (ch.superoperator()[0] - UA.superoperator()[0]).abs().max() < 1e-12
    r = ch(rho_at).matrix[0, 0]
    assert abs(r[1, 1].real - 0.64 * math.exp(-gamma * s)) < 1e-12
    assert abs(r[0, 1] - 0.48 * np.exp((1j * omega - gamma / 2) * s)) < 1e-12
    exact = math.sqrt(0.64 ** 2 * math.exp(-2 * gamma * s) + 0.48 ** 2 * math.exp(-gamma * s))
    assert abs(states.trace_distance(ch(rho_at), ground).item() - exact) < 1e-12

two = TypeIAlgebra([2, 3], [2, 3], precision='double')
H_two = two.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
L_two = [two.operator(torch.randn(1, 2, 3, 3, dtype=c128)) for _ in range(2)]
kernel = [int((torch.linalg.eigvals(G_c[0]).abs() < 1e-10).sum())
          for G_c in dynamics.lindblad_superoperator(two, H_two, L_two)]
print("stationary states per sector:", kernel)
assert all(dim >= 1 for dim in kernel)                             # Corollary 8.13
```

*数值观察*：对这些随机跳跃，代码在每个扇区中恰好找到一个稳态。
由逐扇区应用的定理 D.6（[讲义，8.5 节](../../lecture_notes/zh/08_dynamics.md#s-8-5)），只要每个扇区都有一个忠实的稳态，并且与 $`H_c`$ 以及所有
$`L_{j,c}`$、$`L_{j,c}^*`$ 可交换的矩阵只有 $`1`$ 的倍数，这一点就成立；代码检验这两个假设（在浮点运算中，对上面所用的速率 $`g_j=1`$），
然后断言其结论。

```python
for c, k in enumerate(two.k_factors):                              # Theorem D.6 in every sector
    G_c = dynamics.lindblad_superoperator(two, H_two, L_two)[c][0].numpy()   # on row-major vec(rho_c)
    w_c, V_c = np.linalg.eig(G_c)
    sigma_c = V_c[:, np.argmin(abs(w_c))].reshape(k, k)
    sigma_c = sigma_c / np.trace(sigma_c)
    assert np.abs(G_c @ sigma_c.reshape(-1)).max() < 1e-12            # stationary
    assert np.abs(sigma_c - sigma_c.conj().T).max() < 1e-10 and np.linalg.eigvalsh(sigma_c).min() > 1e-3   # faithful
    ops = [H_two.matrix[0, c, :k, :k].numpy()] + [L.matrix[0, c, :k, :k].numpy() for L in L_two]
    ops += [A.conj().T for A in ops[1:]]
    comm = np.concatenate([np.kron(A, np.eye(k)) - np.kron(np.eye(k), A.T) for A in ops])   # X -> [A, X]
    assert int((np.linalg.svd(comm, compute_uv=False) < 1e-9).sum()) == 1   # commutant = C 1
assert kernel == [1, 1]                                            # the conclusion of Theorem D.6
```

<a id="s-8-6"></a>

## 8.6 可与之比较的工具

QuTiP（`sesolve`；`mesolve` 以自适应 ODE 求解器求解 Lindblad 方程；`mcsolve` 用于
量子轨迹；`steadystate`；`liouvillian`）、dynamiqs（JAX；在 GPU 上成批且可微的
Schrödinger 与 Lindblad 求解器）、QuantumOptics.jl 以及 SciPy（`scipy.linalg.expm`，以及对稀疏的
向量化生成元计算 $`e^{t\mathcal L}\rho`$ 的 `scipy.sparse.linalg.expm_multiply`）在单个 Hilbert 空间上
求解同样的方程，对于上面的问题它们必须给出相同的数值。
在 $`\mathcal M`$ 上，需要把 $`D\times D`$ 的块对角矩阵传给它们；这时向量化的生成元
位于一个 $`D^2`$ 维空间上（稀疏，但各扇区在一次求解中耦合在一起），除非用户
手动把问题拆分到各扇区。本库内建了块结构：逐扇区的精确传播子
与 Lindblad 信道（$`k_c^2\times k_c^2`$ 超算子），作为
`InterSectorChannel` 的扇区间跳跃以及命题 8.8 的压缩，成批的哈密顿量与跳跃
算子（用于无序平均），以及自动微分（后两者 dynamiqs 也有）。它
没有自适应步长控制（RK4 是固定步长的，并且只是近似保正），没有轨迹
求解器，并且 `lindblad_channel` 不接受扇区之间的跳跃；对于这些情形，以 $`M_D`$ 上的 QuTiP 或 dynamiqs
为参照。除 RK4 外，`lindblad_evolve` 还有两种时间步进方法，其每一步对每个步长都是信道（8.7 节）：二阶分裂
格式 `method='kraus'`，以及通过 $`\mathcal L`$ 在 $`\mathcal M`$ 上的矩阵（适用于小扇区）也处理扇区之间跳跃的 `method='exact'`。

<a id="s-8-7"></a>

## 8.7 本身是信道的时间步

*讲义：* [8.7 节](../../lecture_notes/zh/08_dynamics.md#s-8-7)；[命题 8.15](../../lecture_notes/zh/08_dynamics.md#r-8-15)（RK4 不是正的）、[定理 8.16](../../lecture_notes/zh/08_dynamics.md#r-8-16)（每一步都是信道的分裂格式）。

`lindblad_evolve` 有三种时间步进方法，由 `method` 选择：`'rk4'`（默认值，保持不变）是步进 $`R_h=T_4(h\mathcal L)`$；
`'kraus'` 是定理 8.16 的分裂格式 $`\tilde\Psi_h`$，由 $`E_h=e^{-ihH_{\mathrm{eff}}/2}`$、跳跃部分 $`\Phi_{\mathrm{tot}}`$（扇区内部与
扇区之间的跳跃）以及修正 $`M_h^{-1/2}`$ 构成；`'exact'` 是 $`e^{h\mathcal L}`$，由 $`\mathcal L`$ 在 $`\bigoplus_cM_{k_c}`$ 上的
矩阵计算（补零到 $`C k_{\max}^2`$ 个坐标；对常数生成元，每种区间长度计算一次指数，且不使用 `substeps`）。我们先在
$`\omega=\gamma=1`$、$`h=1`$ 的二能级系统上检验命题 8.15：RK4 的一步等于用 `lindblad_superoperator` 的矩阵算出的
$`T_4(h\mathcal L)`$，$`\tau-|a|^2<0`$，并且纯态 $`\cos\frac\theta2|g\rangle+\sin\frac\theta2|e\rangle`$（$`\theta=0.2746`$）被映为
具有负特征值的矩阵；另外两种方法给出的都是态。

```python
T4 = lambda z: 1 + z + z ** 2 / 2 + z ** 3 / 6 + z ** 4 / 24
om, gm, h = 1.0, 1.0, 1.0                                          # omega = gamma: Proposition 8.15(3)
qb = TypeIAlgebra([2], [2], precision='double')                    # basis |g>, |e>
H_qb = qb.from_blocks([torch.tensor([[0.0, 0.0], [0.0, om]], dtype=f64)])
L_qb = [qb.from_blocks([torch.tensor([[0.0, 1.0], [0.0, 0.0]], dtype=f64)])]   # |g><e|
step = lambda rho, method: dynamics.lindblad_evolve(rho, H_qb, L_qb, [0.0, h], rates=[gm], substeps=1,
                                                    method=method)[-1]
gen = dynamics.lindblad_superoperator(qb, H_qb, L_qb, rates=[gm])[0][0].numpy()
R_h = sum(np.linalg.matrix_power(h * gen, n) / math.factorial(n) for n in range(5))   # T_4(hL)
X = states.random_density_matrix(qb)
assert np.abs(step(X, 'rk4').matrix[0, 0].numpy().reshape(-1) - R_h @ X.matrix[0, 0].numpy().reshape(-1)).max() < 1e-14

tau, a = T4(-h * gm), T4(h * (1j * om - gm / 2))
u = gm ** 2 + 4 * om ** 2
c = [-gm * u * (4 * om ** 2 - 3 * gm ** 2) / 384, u ** 2 * (4 * om ** 2 - 5 * gm ** 2) / 4608, gm * u ** 3 / 9216,
     -u ** 4 / 147456]
assert abs(tau - abs(a) ** 2 - h ** 5 * sum(ci * h ** i for i, ci in enumerate(c))) < 1e-15   # Proposition 8.15(3)
print(f"tau - |a|^2 = {tau - abs(a) ** 2:.4f}")

th = 0.2746
psi_th = torch.tensor([math.cos(th / 2), math.sin(th / 2)], dtype=c128)
pure = DensityMatrix(qb, matrix=torch.outer(psi_th, psi_th.conj())[None, None])
low = {m: torch.linalg.eigvalsh(step(pure, m).matrix[0, 0]).min().item() for m in dynamics.LINDBLAD_METHODS}
print({m: f"{v:.2e}" for m, v in low.items()})
assert low['rk4'] < -8e-5 and low['kraus'] > 0 and low['exact'] > 0
x_star = max(r.real for r in np.roots([1, -4, 12, -24]) if abs(r.imag) < 1e-12)
print(f"x_* = {x_star:.4f}")
```

*数值观察*：这里 $`\tau-|a|^2=-0.0091`$，RK4 一步之后的最小特征值为 $`-8.5\cdot10^{-5}`$（打印出的值）；在四位数字的精度下，
$`\theta=0.2746`$ 是使它最小的角度。命题 8.15 讨论的是一个量子比特；在带扇区的代数中，同样的缺陷出现在扇区分量 $`R_{dc}`$
的 Choi 矩阵中，而这一步是 CP 的当且仅当这些矩阵都是正的（命题 6.4）。下面的函数 `step_map` 由 `lindblad_evolve` 在厄米
矩阵上的作用重建一步（输出被取为自伴的，所以在 $`E_{ij}+E_{ji}`$ 与 $`i(E_{ij}-E_{ji})`$ 上求值，再按线性延拓）。对 8.3 节的
粒子丢失，我们检验定理 8.16(1) 与命题 8.8(1)：对从 $`0.1`$ 到 $`1000`$ 的步长，`'kraus'` 与 `'exact'` 的每一步都是 CP 且
TP 的，而 RK4 的一步在 $`h=0.1`$ 时就已经不是 CP 的。

```python
def step_map(alg, H, jumps, rates, h, method):
    """{(d, c): matrix of the sector component R_dc of one step, on row-major vec}"""
    ks_, out = alg.k_factors, {}
    for c, kc in enumerate(ks_):
        E = torch.zeros(kc * kc, alg.C, alg.k_max, alg.k_max, dtype=c128)
        for m in range(kc * kc):
            E[m, c, m // kc, m % kc] = 1.0
        Y = [dynamics.lindblad_evolve(alg.operator(Z), H, jumps, [0.0, h], rates=rates, substeps=1,
                                      method=method)[-1].matrix for Z in (E + E.mH, 1j * (E - E.mH))]
        RE = 0.5 * (Y[0] - 1j * Y[1])                                  # R(E_ij) by linearity
        for d, kd in enumerate(ks_):
            out[d, c] = RE[:, d, :kd, :kd].reshape(kc * kc, kd * kd).T.numpy()
    return out

def choi_min(alg, S):                                              # smallest eigenvalue of all J(R_dc)
    ks_ = alg.k_factors
    return min(np.linalg.eigvalsh(S[d, c].reshape(ks_[d], ks_[d], ks_[c], ks_[c]).transpose(2, 0, 3, 1)
                                  .reshape(ks_[c] * ks_[d], -1)).min() for d, c in S)

def tp_error(alg, S):                                              # Tr sum_d R_dc(X) = Tr X
    ks_ = alg.k_factors
    return max(np.abs(sum(np.trace(S[d, c].reshape(kd, kd, -1), axis1=0, axis2=1) for d, kd in enumerate(ks_))
                      - np.eye(kc).reshape(-1)).max() for c, kc in enumerate(ks_))

for h in (0.1, 3.0, 1000.0):
    for m in dynamics.LINDBLAD_METHODS:
        S = step_map(fock, H_bh, [loss], [0.5], h, m)
        print(f"h = {h:6}, {m:5}: smallest Choi eigenvalue {choi_min(fock, S):+.1e}, trace error {tp_error(fock, S):.0e}")
        if m != 'rk4':
            assert choi_min(fock, S) > -1e-12 and tp_error(fock, S) < 1e-12        # a channel
assert choi_min(fock, step_map(fock, H_bh, [loss], [0.5], 0.1, 'rk4')) < -5e-5
```

*数值观察*：对 $`h=0.1`$，RK4 一步的某个 Choi 矩阵的最小特征值为 $`-9.7\cdot10^{-5}`$；对 $`h=3`$（此时 $`h\mathcal L`$ 的某些
特征值 $`z`$ 满足 $`|T_4(z)|>1`$，RK4 不稳定），它为 $`-9.6\cdot10^{2}`$（见上面的打印）；另外两种方法在舍入误差范围内仍是
信道。对 $`h=0.1`$，所有特征值都满足 $`|T_4(z)|\le1`$：正性的丧失不是一种不稳定性。用 $`M_h^{-1/2}`$ 的重新归一化是由
$`M_h`$ 的特征值算出的；当一个扇区内的衰减率不同时，$`M_h`$ 的条件数可能随 $`h`$ 指数增长；当它超过
$`\varepsilon^{-1/4}`$（$`\varepsilon`$ 为机器精度）时，本库递归地把这一步分成两个步长为 $`h/2`$ 的步，使算出的一步在舍入
误差范围内保迹。

最后看精度。`'exact'` 与 $`M_6`$ 上带跳跃 $`\sqrt\gamma a_1,\sqrt\gamma a_2`$ 的生成元 $`\hat{\mathcal L}`$ 的
$`e^{t\hat{\mathcal L}}`$ 一致（命题 8.8(1) 与 (3)，用 SciPy 由 8.3 节的矩阵 `L_loss` 算出），而另外两种方法到它的迹距离
分别像 $`h^2`$（定理 8.16(3)）和 $`h^4`$ 那样下降：

```python
rho_f = states.random_density_matrix(fock)
ref = dynamics.lindblad_evolve(rho_f, H_bh, [loss], [0.0, 2.0], rates=[0.5], method='exact')[-1]
H6, I6 = embed(H_bh.matrix[0].numpy()), np.eye(6)
G6 = -1j * (np.kron(H6, I6) - np.kron(I6, H6.T))
for L6 in L_loss:
    L6 = math.sqrt(0.5) * L6
    G6 = G6 + np.kron(L6, L6.conj()) - 0.5 * (np.kron(L6.T @ L6, I6) + np.kron(I6, (L6.T @ L6).T))
exact6 = sla.expm(2.0 * G6) @ embed(rho_f.matrix[0].numpy()).reshape(-1)
assert np.abs(exact6 - embed(ref.matrix[0].numpy()).reshape(-1)).max() < 1e-12

for m, subs in (('kraus', (4, 8, 16, 32)), ('rk4', (4, 8, 16, 32))):
    err = [states.trace_distance(dynamics.lindblad_evolve(rho_f, H_bh, [loss], [0.0, 2.0], rates=[0.5], substeps=n,
                                                          method=m)[-1], ref).item() for n in subs]
    ratios = [e0 / e1 for e0, e1 in zip(err[:-1], err[1:])]
    print(m, [f"{e:.1e}" for e in err], "ratios", [f"{r:.1f}" for r in ratios])
    assert all((3.5 < r < 4.5) if m == 'kraus' else (12 < r < 20) for r in ratios)
```

*数值观察*：把 $`h`$ 减半，`'kraus'` 的误差约缩小为 $`1/4`$，RK4 的误差约缩小为 $`1/16`$（打印出的比值）：分裂格式以较低的
阶换取了对每个步长的正性。对含时的 $`H`$ 或含时的跳跃，`'kraus'` 与 `'exact'` 把生成元冻结在每一步的中点，所以由定理
8.16(1)，每一步仍是信道；这种冻结中点格式的阶在讲义中没有证明。

<a id="x-8"></a>

## 习题

讲义的[习题 2 和 4](../../lecture_notes/zh/08_dynamics.md#x-8)：8.4 节扇区 $`N=1`$ 中 $`H_{\mathrm{eff}}`$ 的特征值
$`-\frac i4(\gamma_1+\gamma_2)\pm\sqrt{J^2-\delta^2}`$，$`\delta=(\gamma_1-\gamma_2)/4`$，在 $`J=|\delta|`$ 处只有一个大小为 $`2`$ 的
Jordan 块；以及 $`M_2`$ 上纯退相位（$`L=\sigma_z`$）的谱 $`\{0,0,-2\gamma\pm i\omega\}`$。用数值方法检验这两者。

<details>
<summary>解答</summary>

```python
lam = np.linalg.eigvals(H_eff[1])
delta = (g1 - g2) / 4
root = np.sqrt(J ** 2 - delta ** 2 + 0j)
assert np.allclose(sorted(lam, key=np.real), sorted(-0.25j * (g1 + g2) + np.array([root, -root]), key=np.real))
H_ep = np.array([[-0.5j * g1, delta], [delta, -0.5j * g2]])         # J = delta: the exceptional point
Y = H_ep + 0.25j * (g1 + g2) * np.eye(2)                           # H_ep minus its only eigenvalue
assert np.abs(Y @ Y).max() < 1e-15 and np.linalg.matrix_rank(Y) == 1   # one Jordan block of size 2
```

```python
deph = dynamics.lindblad_superoperator(q2, H_at, [q2.from_blocks([torch.diag(torch.tensor([1.0, -1.0], dtype=f64))])],
                                       rates=[gamma])[0][0]
ev = torch.linalg.eigvals(deph)
assert int((ev.abs() < 1e-12).sum()) == 2
assert np.allclose(sorted(ev.tolist(), key=lambda z: (z.real, z.imag))[:2],
                   [-2 * gamma - 1j * omega, -2 * gamma + 1j * omega], atol=1e-12)
```

</details>

**下一章：** [教程 9. 自旋链](09_spin_chains.md)。
