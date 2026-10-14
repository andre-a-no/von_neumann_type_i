[English](../06_quantum_channels.md) · [Español](../es/06_quantum_channels.md) · [Français](../fr/06_quantum_channels.md) · [Deutsch](../de/06_quantum_channels.md) · **中文** · [日本語](../ja/06_quantum_channels.md) · [Русский](../ru/06_quantum_channels.md)

# 教程 6. 量子信道

本教程配合讲义的[第 6 章](../../lecture_notes/zh/06_quantum_channels.md)。该章证明了矩阵代数直和之间完全正映射的各种描述——Choi 矩阵、
Kraus 算子、Stinespring 扩张——并研究对偶、不动点、转移矩阵和随机信道。这里我们看本库如何存储信道，
并数值地检验这些结果。

在本章中，$\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$，
$\mathcal N=\bigoplus_{d=1}^{C'}M_{k'_d}(\mathbb C)$，其极小中心投影分别为 $P_c$ 和
$Q_d$；$A=\bigoplus_cA_c$；$M_{k'\times k}$ 是复 $k'\times k$ 矩阵构成的空间；
$E_{ij}=|i\rangle\langle j|$ 为矩阵单位；$\mathrm{Tr}$ 是 blunt 迹
$\sum_c\mathrm{Tr}\,A_c$，而 $\mathrm{tr}_w=\sum_cw_c\mathrm{Tr}\,A_c$ 是一个忠实
迹（定理 3.1）。分块矩阵作用在 $\mathbb C^n\otimes\mathbb C^k$ 上，块指标在前。

<a id="s-6-1"></a>

## 6.1 完全正映射与 Choi 定理

*讲义：* [6.1 节](../../lecture_notes/zh/06_quantum_channels.md#s-6-1)：[命题 6.1](../../lecture_notes/zh/06_quantum_channels.md#r-6-1)（转置）、[定理 6.2](../../lecture_notes/zh/06_quantum_channels.md#r-6-2)（Choi）、[推论 6.3](../../lecture_notes/zh/06_quantum_channels.md#r-6-3)（Kraus 秩）和[命题 6.4](../../lecture_notes/zh/06_quantum_channels.md#r-6-4)（直和之间的映射）。

本库有两个类。`Channel` 存储**保扇区**的映射（$C'=C$，$k'=k$，
对 $d\ne c$ 有 $\Phi_{dc}=0$），其 Kraus 算子 $K_l=\bigoplus_cK^{cc}_l$ 属于 $\mathcal M$，以形状为
`(batch, r, C, k_max, k_max)` 的张量存储；`InterSectorChannel` 以形状
`(batch, r, C', C, k'_max, k_max)` 存储[命题 6.4](../../lecture_notes/zh/06_quantum_channels.md#r-6-4) 中的一般映射。`Channel.choi()` 对每个扇区返回
$J(\Phi_{cc})$，而 `Channel.superoperator()` 返回 $\Phi_{cc}$ 作用在按行展开的向量上的矩阵 $S_c$，
$\mathrm{vec}(X)_{(a,b)}=X_{ab}$。二者以不同的顺序包含同样的数
$\Phi(E_{ij})_{ab}=J_{(i,a),(j,b)}=S_{(a,b),(i,j)}$。`Channel.from_superoperator`
在某个 Choi 矩阵不是厄米矩阵，即 $\|J-J^*\|>2\cdot10^{-10}s$（双精度下）时抛出 `ValueError`，其中
$s=\max(1,\|J\|)$，范数取谱范数；否则它把 $J$ 替换为其厄米部分 $(J+J^*)/2$，
像[定理 6.2](../../lecture_notes/zh/06_quantum_channels.md#r-6-2) 的证明那样将其对角化，若某个特征值小于 $-10^{-10}s$ 则抛出 `ValueError`，
并把其余的负特征值置为 $0$。由于 $J^*=J$ 当且仅当对所有 $X$ 有 $\Phi(X^*)=\Phi(X)^*$
（[讲义，6.1 节](../../lecture_notes/zh/06_quantum_channels.md#s-6-1)），在这些舍入容差的范围内，它拒绝的恰好是不保持伴随的映射，
以及其余映射中恰好那些不是 CP 的映射；例如 $X\mapsto iX$ 和转置都会被拒绝。

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, states, channels
from torch_vn_algebra.channels import Channel, InterSectorChannel

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def unit(i, j, k):
    X = torch.zeros(k, k, dtype=dt)
    X[i, j] = 1
    return X

def choi(f, k):
    """Choi matrix sum_ij E_ij (x) f(E_ij) of a linear map f on k x k matrices."""
    return sum(torch.kron(unit(i, j, k), f(unit(i, j, k)).contiguous()) for i in range(k) for j in range(k))

def to_superoperator(J, k):
    """S[(a,b),(i,j)] = J[(i,a),(j,b)], with a leading batch axis."""
    return J.reshape(k, k, k, k).permute(1, 3, 0, 2).reshape(1, k * k, k * k)

ev = torch.linalg.eigvalsh(choi(lambda X: X.T, 2))
print(ev)                                                 # the Choi matrix of the transpose is the swap
assert abs(ev[0] + 1) < 1e-12                             # Proposition 6.1

alg = TypeIAlgebra([2, 3], [2, 3], precision='double')    # M_2 (+) M_3
Phi = channels.random_channel(alg, kraus_rank=2)          # Kraus tensor of shape (1, 2, 2, 3, 3)

def block(Phi, c):
    """The block Phi_cc of a sector-preserving map, as a function on k_c x k_c matrices."""
    ks = Phi.algebra.k_factors
    def f(X):
        blocks = [torch.zeros(k, k, dtype=dt) for k in ks]
        blocks[c] = X
        return Phi(Phi.algebra.from_blocks(blocks)).matrix[0, c, :ks[c], :ks[c]]
    return f

J, S = Phi.choi(), Phi.superoperator()                    # one matrix per sector, (1, k_c^2, k_c^2)
rho = states.random_density_matrix(alg)
out = Phi(rho)
for c, k in enumerate(alg.k_factors):
    assert close(J[c][0], choi(block(Phi, c), k))         # the convention of Theorem 6.2
    w = torch.linalg.eigvalsh(J[c][0])
    assert w.min() > -1e-12 and int((w > 1e-10).sum()) == 2      # J >= 0, Kraus rank 2 (Corollary 6.3)
    assert close(S[c][0], to_superoperator(J[c][0], k)[0])
    assert close(S[c][0] @ rho.matrix[0, c, :k, :k].reshape(-1), out.matrix[0, c, :k, :k].reshape(-1))
assert Phi.is_trace_preserving()
assert close(Channel.from_superoperator(alg, S)(rho), out)       # back from the superoperator
try:
    Channel.from_superoperator(alg, [to_superoperator(choi(lambda X: X.T, k), k) for k in (2, 3)])
    raise AssertionError("the transpose was accepted")
except ValueError as e:
    print(e)                                              # superoperator is not completely positive
try:
    Channel.from_superoperator(alg, [1j * torch.eye(k * k, dtype=dt)[None] for k in (2, 3)])
    raise AssertionError("X -> iX was accepted")
except ValueError as e:
    print(e)                                              # superoperator does not preserve Hermiticity
```

<a id="s-6-2"></a>

## 6.2 Stinespring 扩张与 Kraus 算子的自由度

*讲义：* [6.2 节](../../lecture_notes/zh/06_quantum_channels.md#s-6-2)：[定理 6.5](../../lecture_notes/zh/06_quantum_channels.md#r-6-5)（Stinespring）和[定理 6.7](../../lecture_notes/zh/06_quantum_channels.md#r-6-7)（酉自由度）。对 6.1 节中信道的块 $c=1$，代码检验 $V$ 是等距，$\Phi(\rho)=\mathrm{Tr}_E(V\rho V^*)$ 和 $\Phi^*(B)=V^*(1\otimes B)V$ 成立，并且 Kraus 算子 $\sum_lu_{jl}K_l$ 给出同一个信道。

在代码中，$V$ 的环境指标是 Kraus 张量的第一个指标：

```python
K = Phi.kraus[0, :, 1, :3, :3]                            # the two Kraus operators of block c = 1
V = K.reshape(2 * 3, 3)                                   # V = sum_l e_l (x) K_l : C^3 -> C^2 (x) C^3
assert close(V.conj().T @ V, torch.eye(3, dtype=dt))      # TP <=> V is an isometry
r1 = rho.matrix[0, 1, :3, :3]
TrE = torch.einsum('lalb->ab', (V @ r1 @ V.conj().T).reshape(2, 3, 2, 3))
assert close(TrE, out.matrix[0, 1, :3, :3])               # Phi(rho) = Tr_E(V rho V^*)
B1 = torch.randn(3, 3, dtype=dt)
Bop = alg.from_blocks([torch.zeros(2, 2, dtype=dt), B1])
assert close(V.conj().T @ torch.kron(torch.eye(2, dtype=dt), B1) @ V,
             Phi.adjoint()(Bop).matrix[0, 1, :3, :3])    # Phi^*(B) = V^*(1 (x) B)V

u = alg.random_unitary(3)                                 # Theorem 6.7 with n = 3: pad with a zero operator
Kpad = torch.cat([Phi.kraus, torch.zeros_like(Phi.kraus[:, :1])], dim=1)
PhiL = Channel(alg, torch.einsum('jl,blcxy->bjcxy', u, Kpad))
assert close(PhiL(rho), out) and all(close(a, b) for a, b in zip(PhiL.choi(), J))
```

<a id="s-6-3"></a>

## 6.3 关于迹的对偶

*讲义：* [6.3 节](../../lecture_notes/zh/06_quantum_channels.md#s-6-3)：关于 $(w,v)$ 的对偶 $\Phi^\sharp$、[命题 6.8](../../lecture_notes/zh/06_quantum_channels.md#r-6-8) 和[推论 6.9](../../lecture_notes/zh/06_quantum_channels.md#r-6-9)。

在本库中，`Channel.adjoint()` 是 $\Phi^*$，而 `InterSectorChannel.adjoint(trace=...)` 是对两个代数
使用同一个具名迹（`'blunt'`、`'norm'`、`'tau_vN'`）时的 $\Phi^\sharp$；`tau_vN` 的权重
$1/(Ck_c)$ 和 $1/(C'k'_d)$ 还会受到扇区个数的影响。下面对一个随机映射
$M_2\oplus M_3\to M_1\oplus M_2\oplus M_2$ 和三种迹检验[命题 6.8](../../lecture_notes/zh/06_quantum_channels.md#r-6-8)，并对 $M_1\oplus M_2$ 上的粒子丢失
$\Phi(x)=(x_0+\mathrm{Tr}\,x_1)\oplus0$ 检验[推论 6.9(a)](../../lecture_notes/zh/06_quantum_channels.md#r-6-9)；它对 blunt 迹是保迹的，对 `Tr_norm` 则不是
（[讲义，6.3 节](../../lecture_notes/zh/06_quantum_channels.md#s-6-3)）：

```python
alg_out = TypeIAlgebra([1, 2, 2], [1, 2, 2], precision='double')
Psi = channels.random_inter_sector_channel(alg, alg_out, kraus_rank=2)     # Kraus tensor (1, 2, 3, 2, 2, 3)
x = alg.operator(torch.randn(1, 2, 3, 3, dtype=dt))
y = alg_out.operator(torch.randn(1, 3, 2, 2, dtype=dt))
weights = {'blunt': lambda k, C: 1.0, 'norm': lambda k, C: 1 / k, 'tau_vN': lambda k, C: 1 / (C * k)}
trace = {'blunt': lambda X: X.Tr_blunt(), 'norm': lambda X: X.Tr_norm(), 'tau_vN': lambda X: X.tau_vN()}

def z(a, t):                                              # the central element sum_c w_c P_c of trace t
    return a.central(torch.tensor([weights[t](k, a.C) for k in a.k_factors], dtype=torch.float64))

for t in trace:
    dual = Psi.adjoint(trace=t)
    assert (trace[t](Psi(x) @ y) - trace[t](x @ dual(y))).abs().max() < 1e-10   # the definition
    assert close(dual(y), z(alg, t).inverse() @ Psi.adjoint()(z(alg_out, t) @ y))   # Proposition 6.8
dev = [(Psi.adjoint(trace=t)(alg_out.identity()) - alg.identity()).matrix.abs().max().item() for t in trace]
print(dev)                         # TP only for the blunt traces: the other duals are not unital
assert Psi.is_trace_preserving() and dev[0] < 1e-12 and min(dev[1:]) > 0.1

alg01 = TypeIAlgebra([1, 2], [1, 2], precision='double')
loss = InterSectorChannel.from_blocks(alg01, alg01, {(0, 0): [torch.ones(1, 1)],
                                                     (0, 1): [torch.tensor([[1.0, 0.0]]), torch.tensor([[0.0, 1.0]])]})
assert loss.is_trace_preserving()
assert close(loss.adjoint(trace='norm')(alg01.identity()), alg01.central([1.0, 2.0]))
s = states.random_density_matrix(alg01)
assert abs(loss(s).Tr_norm() - (s @ alg01.central([1.0, 2.0])).Tr_norm()) < 1e-12   # Corollary 6.9(a)
print(s.Tr_norm().real.item(), loss(s).Tr_norm().real.item())      # Tr_norm is not conserved
```

<a id="s-6-4"></a>

## 6.4 退相干

*讲义：* [6.4 节](../../lecture_notes/zh/06_quantum_channels.md#s-6-4)：[定理 6.10](../../lecture_notes/zh/06_quantum_channels.md#r-6-10)（保单位信道的不动点）、[定理 6.11](../../lecture_notes/zh/06_quantum_channels.md#r-6-11)（保单位信道使熵增加）和[命题 6.12](../../lecture_notes/zh/06_quantum_channels.md#r-6-12)（转移矩阵）。

本库的标准信道在每个块内部作用（$p,\gamma\in[0,1]$）：

| 信道 | $\Phi(\rho)_c$ | 块 $c$ 中的 Kraus 算子 |
|---|---|---|
| `dephasing_channel` | $(1-p)\rho_c+p\mathrm{diag}\,\rho_c$ | $\sqrt{1-p}\,1$, $\sqrt p\,E_{ii}$ |
| `depolarizing_channel` | $(1-p)\rho_c+p\mathrm{Tr}(\rho_c)\,1/k_c$ | $\sqrt{1-p}\,1$, $\sqrt{p/k_c}\,E_{ij}$ |
| `center_expectation` | $\mathrm{Tr}(\rho_c)\,1/k_c$ | $\sqrt{1/k_c}\,E_{ij}$ |
| `amplitude_damping_channel` | 能级 $i>0$ 衰减到 $0$ | $E_{00}+\sqrt{1-\gamma}\sum_{i>0}E_{ii}$, $\sqrt\gamma\,E_{0i}$ |
| `random_mixed_unitary_channel` | $\sum_jp_jU_j\rho U_j^*$ | $\sqrt{p_j}\,U_j$ |

对表中 $M_2\oplus M_3$ 上的每个信道，代码计算其不动点空间（$S_c-1$ 的零空间）的维数和其 Kraus 算子的
交换子的维数：对保单位信道二者一致（[定理 6.10](../../lecture_notes/zh/06_quantum_channels.md#r-6-10)），对振幅阻尼二者维数相同但是不同的空间。
然后它在 $100$ 个随机态上检验[定理 6.11](../../lecture_notes/zh/06_quantum_channels.md#r-6-11)，并对 6.3 节的映射 `Psi` 和粒子丢失检验[命题 6.12](../../lecture_notes/zh/06_quantum_channels.md#r-6-12)。

```python
def fixed_dims(Phi):
    """dim {X in M_{k_c} : Phi_cc(X) = X} for every sector: the null space of S_c - 1."""
    return [int((torch.linalg.svdvals(S_c[0] - torch.eye(k * k, dtype=dt)) < 1e-9).sum())
            for S_c, k in zip(Phi.superoperator(), Phi.algebra.k_factors)]

def commutant_dims(Phi):
    """dim of the commutant of {K_l, K_l^*} in every block (row-major vec, as in Chapter 2)."""
    out = []
    for c, k in enumerate(Phi.algebra.k_factors):
        I, Ks = torch.eye(k, dtype=dt), Phi.kraus[0, :, c, :k, :k]
        ops = [A.contiguous() for A in torch.cat([Ks, Ks.mH])]
        L = torch.cat([torch.kron(A, I) - torch.kron(I, A.T.contiguous()) for A in ops])
        out.append(k * k - int((torch.linalg.svdvals(L) > 1e-9).sum()))
    return out

zoo = {'dephasing': channels.dephasing_channel(alg, 0.3),
       'depolarizing': channels.depolarizing_channel(alg, 0.3),
       'centre': channels.center_expectation(alg),
       'mixed unitary': channels.random_mixed_unitary_channel(alg, 2),
       'amplitude damping': channels.amplitude_damping_channel(alg, 0.4)}
for name, ch in zoo.items():
    print(f"{name:18s} unital {ch.is_unital()!s:5s} fixed {fixed_dims(ch)} commutant {commutant_dims(ch)}")
    assert ch.is_trace_preserving()
    if ch.is_unital():
        assert fixed_dims(ch) == commutant_dims(ch)                         # Theorem 6.10
assert fixed_dims(zoo['dephasing']) == [2, 3] and fixed_dims(zoo['depolarizing']) == [1, 1]
E00 = alg.from_blocks([unit(0, 0, 2), unit(0, 0, 3)])
E01 = alg.from_blocks([unit(0, 1, 2), unit(0, 1, 3)])
assert close(zoo['amplitude damping'](E00), E00) and not close(E00 @ E01, E01 @ E00)   # equal dims, different spaces
assert fixed_dims(zoo['amplitude damping']) == commutant_dims(zoo['amplitude damping']) == [1, 1]

rhos = states.random_density_matrix(alg, batch_size=100)
for name in ('dephasing', 'depolarizing', 'centre', 'mixed unitary'):                 # Theorem 6.11
    assert (zoo[name](rhos).von_neumann_entropy() >= rhos.von_neumann_entropy() - 1e-12).all()
    assert (zoo[name](rhos).purity() <= rhos.purity() + 1e-12).all()
mm = states.maximally_mixed_state(alg)
S_after = channels.amplitude_damping_channel(alg, 1.0)(mm).von_neumann_entropy().item()
print(mm.von_neumann_entropy().item(), S_after)          # not unital: the entropy drops
assert abs(S_after + 0.4 * math.log(0.4) + 0.6 * math.log(0.6)) < 1e-12

T = Psi.transition_matrix()[0]                            # (C', C) = (3, 2)
p, q = rho.sector_probabilities()[0], Psi(rho).sector_probabilities()[0]
assert torch.allclose(T.sum(0), torch.ones(2, dtype=T.dtype))                       # Proposition 6.12
assert torch.allclose(q, Psi.transition_matrix(rho)[0] @ p)
classical = channels.center_expectation(alg_out) @ Psi @ channels.center_expectation(alg)
Tp = T @ p
assert close(classical(rho), alg_out.central([Tp[d] / k for d, k in enumerate(alg_out.k_factors)]))
assert torch.allclose(loss.transition_matrix()[0], torch.tensor([[1.0, 1.0], [0.0, 0.0]], dtype=T.dtype))
```

在完全振幅阻尼下，$M_2\oplus M_3$ 的最大混合态的熵从
$\log5\approx1.609$ 降到 $H(\tfrac25,\tfrac35)\approx0.673$：非保单位的信道可以起到纯化作用。

在 $M_2$ 上，退相位 $D_p$ 与振幅阻尼 $A_\gamma$ 可交换，它们的复合把相干项乘以 $(1-p)\sqrt{1-\gamma}$；
取 $1-p=e^{-t/T_\varphi}$ 和 $1-\gamma=e^{-t/T_1}$，就得到 $1/T_2=1/(2T_1)+1/T_\varphi$（[讲义，6.4 节](../../lecture_notes/zh/06_quantum_channels.md#s-6-4)，
*量子比特的相干项*）。代码检验这些恒等式以及复合律 $D_{p(t)}D_{p(s)}=D_{p(t+s)}$、$A_{\gamma(t)}A_{\gamma(s)}=A_{\gamma(t+s)}$：

```python
q2 = TypeIAlgebra([2], [2], precision='double')
p_, g_ = 0.3, 0.4
Dp, Ag = channels.dephasing_channel(q2, p_), channels.amplitude_damping_channel(q2, g_)
r2 = states.random_density_matrix(q2)
assert close(Dp(Ag(r2)), Ag(Dp(r2)))                                         # D_p A_gamma = A_gamma D_p
coh = Dp(Ag(r2)).matrix[0, 0, 0, 1] / r2.matrix[0, 0, 0, 1]
assert abs(coh - (1 - p_) * math.sqrt(1 - g_)) < 1e-12                       # coherence factor
T1, Tphi, t1, t2 = 2.0, 3.0, 0.7, 1.1
gam = lambda t: 1 - math.exp(-t / T1)
assert close(channels.amplitude_damping_channel(q2, gam(t1))(channels.amplitude_damping_channel(q2, gam(t2))(r2)),
             channels.amplitude_damping_channel(q2, gam(t1 + t2))(r2))       # A_gamma(t) A_gamma(s) = A_gamma(t+s)
pdeph = lambda t: 1 - math.exp(-t / Tphi)
assert close(channels.dephasing_channel(q2, pdeph(t1))(channels.dephasing_channel(q2, pdeph(t2))(r2)),
             channels.dephasing_channel(q2, pdeph(t1 + t2))(r2))           # D_p(t) D_p(s) = D_p(t+s)
both = channels.dephasing_channel(q2, pdeph(t1))(channels.amplitude_damping_channel(q2, gam(t1))(r2))
coh_t = (both.matrix[0, 0, 0, 1] / r2.matrix[0, 0, 0, 1]).item()
assert abs(coh_t - math.exp(-t1 * (1 / (2 * T1) + 1 / Tphi))) < 1e-12       # 1/T_2 = 1/(2T_1) + 1/T_phi
```

<a id="s-6-5"></a>

## 6.5 随机信道

*讲义：* [6.5 节](../../lecture_notes/zh/06_quantum_channels.md#s-6-5)：[引理 6.13](../../lecture_notes/zh/06_quantum_channels.md#r-6-13) 和[定理 6.15](../../lecture_notes/zh/06_quantum_channels.md#r-6-15)（随机信道的输出纯度）。

与本库（`random_isometry`）一样，**Haar 随机等距** $V:\mathbb C^k\to\mathbb C^N$ 定义为
$G=QR$ 中的因子 $Q$，其中 $G$ 是 $N\times k$ 的 Ginibre 矩阵（各元素独立，密度为
$\pi^{-1}e^{-|z|^2}$），$Q^*Q=1$，$R$ 是对角元为正的上三角矩阵（在 CPU 上，本库先取 LAPACK
分解 $G=Q'R'$，再把 $Q'$ 的第 $i$ 列乘以 $R'_{ii}/|R'_{ii}|$；与定理 5.3 中一样，只要所有
$R'_{ii}\ne0$，即只要 $G$ 的秩为 $k$，这就给出上述分解，而由引理 6.13 的证明，这几乎必然成立；
对 GPU 上的大批量，本库使用 Gram–Schmidt 正交化，它直接给出该分解）。
`random_channel(alg, r)` 在每个块中独立地抽取这样的 $V_c:\mathbb C^{k_c}\to\mathbb C^r\otimes\mathbb C^{k_c}$，
并使用定理 6.5 中的 Kraus 算子。

由[讲义，定理 6.15](../../lecture_notes/zh/06_quantum_channels.md#r-6-15)，对大小为 $k$ 的扇区中的纯态 $\psi$ 和由 `random_channel` 以 Kraus 秩 $r$ 抽取的 $\Phi$，
$\mathbb E\mathrm{Tr}\,\Phi(\psi)^2=(k+r)/(kr+1)$ 且 $\mathbb E\,\Phi(\psi)=1_k/k$。我们用 $20\,000$ 个样本估计二者，
并给出均值的标准误差：

```python
B = 20_000
for r in (1, 2, 4):
    Phi_r = channels.random_channel(alg, kraus_rank=r, batch_size=B)
    for c, k in enumerate(alg.k_factors):
        psi = torch.randn(k, 1, dtype=dt)
        psi = psi / psi.norm()
        out_r = Phi_r(alg.from_blocks([psi @ psi.conj().T if d == c else torch.zeros(kk, kk, dtype=dt)
                                       for d, kk in enumerate(alg.k_factors)]))     # a batch of B outputs
        purity = (out_r @ out_r).trace.real
        mean, err = purity.mean().item(), purity.std().item() / math.sqrt(B)
        exact = (k + r) / (k * r + 1)
        print(f"k = {k}, r = {r}: E Tr Phi(psi)^2 = {mean:.4f} +- {err:.4f}   (exact {exact:.4f})")
        assert abs(mean - exact) < 5 * err + 1e-12
        avg = out_r.matrix.mean(0)[c, :k, :k]
        assert (avg - torch.eye(k, dtype=dt) / k).abs().max() < 0.02          # E Phi(psi) = 1/k
```

<a id="s-6-6"></a>

## 6.6 可与之比较的工具

QuTiP（`to_choi`、`to_kraus`、`to_super`、`to_stinespring`）、Qiskit（`qiskit.quantum_info`：
类 `Kraus`、`Choi`、`SuperOp`、`Stinespring` 以及 `random_quantum_channel`；Qiskit Aer 中的
`amplitude_damping_error` 等噪声模型）、toqito（`toqito.channel_ops`、`toqito.channel_props`）以及
QETLAB（`ChoiMatrix`、`KrausOperators`、`ApplyMap`、`RandomSuperoperator`）可以在 6.1–6.2 节的
各种表示之间转换、检验完全正性，并且（大多数）可以抽取随机信道，
适用于完整矩阵代数之间的映射 $M_n\to M_m$。
它们的约定各不相同（按列或按行堆叠、Choi 矩阵中因子的顺序），所以
在比较数值之前要先比较约定。$\mathcal M$ 上的保扇区信道在那里可以
当作 $M_D$（$D=\sum_ck_c$）上的信道处理，其超算子的大小为 $D^2\times D^2$，
而不是每个扇区 $k_c^2\times k_c^2$。本库增加的是块结构本身：按扇区并成批存储的
Kraus 算子、不同代数之间的映射（`InterSectorChannel`）、
关于第 3 章各种迹的对偶，以及扇区之间的转移矩阵。

<a id="x-6"></a>

## 习题

讲义中[习题](../../lecture_notes/zh/06_quantum_channels.md#x-6)的数值检验：

1. （Werner–Holevo）用 `Channel.from_superoperator` 检验 $M_3$ 上的 $\Phi(\rho)=\big(\mathrm{Tr}(\rho)1-\rho^T\big)/(k-1)$
   是 Kraus 秩为 $3$ 的保单位信道。
2. 检验 $M_2\oplus M_3$ 上的退相位信道在块 $c$ 中的 Kraus 秩为 $k_c$，而本库存储了 $1+k_{\max}$ 个 Kraus 算子。
3. 检验以 $X^aZ^b/k$ 为 Kraus 算子的混合酉信道就是 $M_3$ 上的 `center_expectation`。
4. 对本库的三种迹检验 $(\Psi\circ\Phi)^\sharp=\Phi^\sharp\circ\Psi^\sharp$。

<details>
<summary>解答</summary>

*习题 1.*

```python
k = 3
M3 = TypeIAlgebra([k], [k], precision='double')
J_WH = choi(lambda X: (X.trace() * torch.eye(k, dtype=dt) - X.T) / (k - 1), k)
w = torch.linalg.eigvalsh(J_WH)
assert w.min() > -1e-12 and int((w > 1e-10).sum()) == k * (k - 1) // 2
WH = Channel.from_superoperator(M3, [to_superoperator(J_WH, k)])
assert WH.is_trace_preserving() and WH.is_unital()
r3 = states.random_density_matrix(M3).matrix[0, 0]
assert close(WH(M3.operator(r3[None, None])).matrix[0, 0], (torch.eye(k, dtype=dt) - r3.T) / (k - 1))
```

*习题 2.*

```python
ranks = [int((torch.linalg.eigvalsh(J_c[0]) > 1e-10).sum()) for J_c in zoo['dephasing'].choi()]
print(ranks, zoo['dephasing'].kraus_rank)                 # [2, 3] versus 4 stored operators
assert ranks == alg.k_factors and zoo['dephasing'].kraus_rank == 1 + max(alg.k_factors)
```

*习题 3.*

```python
Xs = torch.roll(torch.eye(k, dtype=dt), 1, dims=0)
Zs = torch.diag(torch.exp(2j * math.pi * torch.arange(k, dtype=torch.float64) / k)).to(dt)
W = [torch.linalg.matrix_power(Xs, a) @ torch.linalg.matrix_power(Zs, b) / k for a in range(k) for b in range(k)]
twirl = channels.kraus_channel([M3.operator(Wab[None, None]) for Wab in W])
r3op = states.random_density_matrix(M3)
assert close(twirl(r3op), channels.center_expectation(M3)(r3op))
```

*习题 4.*

```python
Psi2 = channels.random_inter_sector_channel(alg_out, alg, kraus_rank=2)
for t in trace:
    assert close((Psi2 @ Psi).adjoint(trace=t)(x), Psi.adjoint(trace=t)(Psi2.adjoint(trace=t)(x)))
```

</details>

**下一章：** [教程 7](07_tensor_products.md)：库中的张量积、偏迹与纠缠。
