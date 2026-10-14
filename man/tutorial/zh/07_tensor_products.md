[English](../07_tensor_products.md) · [Español](../es/07_tensor_products.md) · [Français](../fr/07_tensor_products.md) · [Deutsch](../de/07_tensor_products.md) · **中文** · [日本語](../ja/07_tensor_products.md) · [Русский](../ru/07_tensor_products.md)

# 教程 7. 张量积与纠缠

本教程配合讲义的[第 7 章](../../lecture_notes/zh/07_tensor_products.md)。该章证明了两个带扇区的代数的普通张量积与融合张量积的结构，
把偏迹定义为一个嵌入的对偶，并研究纠缠与随机态。这里我们用本库构造这些积，计算偏迹和纠缠熵，
并数值地检验这些定理。

在本章中，$\mathcal M_1=\bigoplus_{c=1}^{C_1}M_{k_c}(\mathbb C)$ 作用在
$H_1=\bigoplus_c\mathbb C^{k_c}$ 上，极小中心投影为 $P_c$；
$\mathcal M_2=\bigoplus_{d=1}^{C_2}M_{l_d}(\mathbb C)$ 作用在 $H_2=\bigoplus_d\mathbb C^{l_d}$ 上，
极小中心投影为 $Q_d$。$\mathrm{Tr}$ 是算子所作用的空间上的迹；在 $\mathcal M_i$ 上它就是第 3 章的
blunt 迹。对 $\mathbb C^k\otimes\mathbb C^l$ 上的 $X$，记
$X_{(a,b),(a',b')}=\langle e_a\otimes e_b,X\,e_{a'}\otimes e_{b'}\rangle$，
并通过 $e_a\otimes e_b\mapsto e_{(a-1)l+b}$ 把 $\mathbb C^k\otimes\mathbb C^l$ 与 $\mathbb C^{kl}$ 等同
（Kronecker 约定），从而 $(A\otimes B)_{(a,b),(a',b')}=A_{aa'}B_{bb'}$。

<a id="s-7-1"></a>

## 7.1 两个代数的张量积

*讲义：* [7.1 节](../../lecture_notes/zh/07_tensor_products.md#s-7-1)：[定理 7.1](../../lecture_notes/zh/07_tensor_products.md#r-7-1)（$\mathcal M_1\otimes\mathcal M_2$ 的结构）。

`tensor_product(alg1, alg2)` 构造这个代数。扇区 $(c,d)$ 的编号为 $s=cC_2+d$（代码中编号
从 $0$ 开始），荷为 $(q_c,r_d)$，在扇区内部 $e_a\otimes e_b$ 的编号为 $a\,l_d+b$。
`kron(A, B, alg12)` 计算 $A\otimes B$。扇区顺序与 $H_1\otimes H_2$ 的 Kronecker 顺序相差一个基的置换，
下面的 `kron_rows` 显式地构造这个置换；借助它，我们与全矩阵的 Kronecker 积进行比较。

```python
import math
import numpy as np
import torch
from torch_vn_algebra import (TypeIAlgebra, DensityMatrix, states, composite,
                              tensor_product, fused_tensor_product, kron, partial_trace)

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def full(X):
    """An operator of (+)_c M_{k_c} as a block-diagonal matrix on (+)_c C^{k_c}, shape (batch, D, D)."""
    return torch.stack([torch.block_diag(*[x[c, :k, :k] for c, k in enumerate(X.algebra.k_factors)])
                        for x in X.matrix])

def kron_rows(alg12):
    """For every sector of a product algebra, the positions of its basis vectors in H1 (x) H2 (Kronecker order)."""
    k1, k2 = alg12.tensor_factors[0].k_factors, alg12.tensor_factors[1].k_factors
    return [torch.tensor([(sum(k1[:c]) + a) * sum(k2) + sum(k2[:d]) + b
                          for c, d, off in pairs for a in range(k1[c]) for b in range(k2[d])])
            for pairs in alg12.product_blocks]              # pair blocks (c, d, offset), by offset

def full12(X):
    """An operator of a (plain or fused) product algebra as a matrix on H1 (x) H2 in Kronecker order."""
    f1, f2 = X.algebra.tensor_factors
    N = sum(f1.k_factors) * sum(f2.k_factors)
    out = torch.zeros(X.matrix.shape[0], N, N, dtype=X.matrix.dtype)
    for s, r in enumerate(kron_rows(X.algebra)):
        out[:, r[:, None], r[None, :]] = X.matrix[:, s, :len(r), :len(r)]
    return out

def random_operator(alg):
    return alg.from_blocks([torch.randn(k, k, dtype=dt) for k in alg.k_factors])

alg1 = TypeIAlgebra([1, 2], [1, 2], precision='double')            # C (+) M_2
alg2 = TypeIAlgebra([2, 3], [2, 3], precision='double')            # M_2 (+) M_3
alg12 = tensor_product(alg1, alg2)
print(alg12.k_factors, alg12.charges, sum(k * k for k in alg12.k_factors))     # 4 sectors, dimension 65
assert sum(k * k for k in alg12.k_factors) == 5 * 13                           # Theorem 7.1(b)

A, A2, B, B2 = random_operator(alg1), random_operator(alg1), random_operator(alg2), random_operator(alg2)
AB = kron(A, B, alg12)
assert close(full12(AB)[0], torch.kron(full(A)[0], full(B)[0]))               # Theorem 7.1(a)
assert close(AB @ kron(A2, B2, alg12), kron(A @ A2, B @ B2, alg12))
for tr in ('Tr_blunt', 'Tr_norm', 'tau_vN'):                                   # Theorem 7.1(c)
    assert (getattr(AB, tr)() - getattr(A, tr)() * getattr(B, tr)()).abs().max() < 1e-10
```

<a id="s-7-2"></a>

## 7.2 偏迹

*讲义：* [7.2 节](../../lecture_notes/zh/07_tensor_products.md#s-7-2)：$\mathrm{Tr}_2$ 的定义、[命题 7.2](../../lecture_notes/zh/07_tensor_products.md#r-7-2)、[定理 7.3](../../lecture_notes/zh/07_tensor_products.md#r-7-3)（偏迹是一个信道）和[命题 7.4](../../lecture_notes/zh/07_tensor_products.md#r-7-4)（指标为 $m^2$ 的条件期望）。

`partial_trace(X, keep=1)` 把每个块 $\rho_{cd}$ 重排为形状 $(k_c,l_d,k_c,l_d)$ 再求迹；
`composite.partial_trace_channel(alg12, keep=1)` 是定理 7.3 中的信道，其 Kraus 算子为从扇区 $(c,d)$
到扇区 $c$ 的 $1_{k_c}\otimes\langle j|$。我们把两者与定义、一个独立的 `numpy.einsum` 以及对偶关系进行核对：

```python
rho = states.random_density_matrix(alg12, batch_size=3)
rho1, rho2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
assert isinstance(rho1, DensityMatrix)                                      # states give states
R = full12(rho).numpy().reshape(3, 3, 5, 3, 5)                              # indices (batch, i, j, i', j')
assert np.abs(full(rho1).numpy() - np.einsum('xijkj->xik', R)).max() < 1e-12
assert np.abs(full(rho2).numpy() - np.einsum('xijil->xjl', R)).max() < 1e-12
one2 = alg2.identity()
assert (states.expectation(rho, kron(A, one2, alg12)) - states.expectation(rho1, A)).abs().max() < 1e-12
assert close(partial_trace(AB, keep=1), A * B.trace.item())                 # Proposition 7.2(b)

Ptr = composite.partial_trace_channel(alg12, keep=1)
print(Ptr)
assert Ptr.is_trace_preserving()
assert close(Ptr(rho), rho1)
assert close(Ptr.adjoint('blunt')(A), kron(A, one2, alg12))                 # Theorem 7.3: the dual
```

模块 `torch_vn_algebra.inclusions` 实现了包含、它们的条件期望和指标；
`from_tensor_product(alg1, m)` 构造[命题 7.4](../../lecture_notes/zh/07_tensor_products.md#r-7-4) 中的包含 $\mathcal M_1\otimes1\subset\mathcal M_1\otimes M_m$，
其包含矩阵为 $\Lambda=m1_{C_1}$，指标为 $m^2$（[讲义，7.2 节](../../lecture_notes/zh/07_tensor_products.md#s-7-2)）：

```python
from torch_vn_algebra import inclusions
m = 3
inc = inclusions.from_tensor_product(alg1, m)                # alg1 (x) 1 c alg1 (x) M_3
X = random_operator(inc.M)
for trace in ('Tr_blunt', 'Tr_norm', 'tau_vN'):
    E = inc.conditional_expectation(trace.replace('Tr_', ''))
    assert close(E(X), partial_trace(X, keep=1) * (1 / m))   # Proposition 7.4(a): the same E for every trace
    assert abs(getattr(inc.embed(E(X)), trace)() - getattr(X, trace)()).max() < 1e-10
    assert close(E.index(), inc.M.identity() * m**2)         # Proposition 7.4(c)

def unit(a, b):
    e = torch.zeros(m, m, dtype=dt)
    e[a, b] = 1
    return math.sqrt(m) * kron(alg1.identity(), inc.M.tensor_factors[1].from_blocks([e]), inc.M)

u = [unit(a, b) for a in range(m) for b in range(m)]         # the quasi-basis of Proposition 7.4(b)
assert close(sum((v @ inc.embed(E(v.adjoint() @ X))).matrix for v in u), X)
assert close(sum((v @ v.adjoint()).matrix for v in u), inc.M.identity() * m**2)
print(inc.matrix.tolist(), inc.component_indices())           # Lambda = 3 * 1_2, index 9 per component
```

<a id="s-7-3"></a>

## 7.3 荷的融合

*讲义：* [7.3 节](../../lecture_notes/zh/07_tensor_products.md#s-7-3)：[命题 7.5](../../lecture_notes/zh/07_tensor_products.md#r-7-5)（融合积）和[定理 7.6](../../lecture_notes/zh/07_tensor_products.md#r-7-6)（融合积的约化态）。

`fused_tensor_product(alg1, alg2, fuse)`（默认的 `fuse` 为加法）按融合后的荷对扇区排序，扇区内部的
对 $(c,d)$ 按字典序排列，每个对都采用基 $e_a\otimes e_b$。
例子：硬核玻色子（[讲义，7.3 节](../../lecture_notes/zh/07_tensor_products.md#s-7-3)）。两个格点上保持粒子数的算子构成 $M_1\oplus M_2\oplus M_1$（荷 $0,1,2$），
三个格点上构成 $M_1\oplus M_3\oplus M_3\oplus M_1$；由[命题 7.5(c)](../../lecture_notes/zh/07_tensor_products.md#r-7-5)，它们的融合积就是五个格点的总占据数所对应的
$\{N\}'$，且 $n_q=\binom5q$。

```python
bos1 = TypeIAlgebra([1, 2, 1], [1, 2, 1], precision='double', charges=[0, 1, 2])
bos2 = TypeIAlgebra([1, 3, 3, 1], [1, 3, 3, 1], precision='double', charges=[0, 1, 2, 3])
F = fused_tensor_product(bos1, bos2)                         # U(1): total particle number
P = tensor_product(bos1, bos2)                               # particle numbers conserved separately
print(F.charges, F.k_factors, sum(k * k for k in F.k_factors), sum(k * k for k in P.k_factors))
assert F.k_factors == [math.comb(5, q) for q in range(6)]
assert sum(k * k for k in F.k_factors) == 252 and sum(k * k for k in P.k_factors) == 120

N1, N2 = full(bos1.central([0.0, 1.0, 2.0]))[0], full(bos2.central([0.0, 1.0, 2.0, 3.0]))[0]
N = torch.kron(N1, torch.eye(8, dtype=dt)) + torch.kron(torch.eye(4, dtype=dt), N2)
X = full12(F.operator(torch.randn(1, F.C, F.k_max, F.k_max, dtype=dt)))[0]
assert close(X @ N, N @ X)                                   # Proposition 7.5(c)
A, B = random_operator(bos1), random_operator(bos2)
assert close(full12(kron(A, B, F)), full12(kron(A, B, P)))   # Proposition 7.5(b)
```

融合代数的维数是 $\sum_q\binom5q^2=\binom{10}5=252$（又是 Vandermonde 恒等式），普通积的维数是
$(1+4+1)(1+9+9+1)=120$（定理 7.1(b)）。由定理 7.1(a) 和命题 7.5(b)，普通积由融合积中关于对 $(c,d)$
块对角的元素组成；其余 $132$ 维由这样的算子张成：它们把一个对块 $(c,d)$ 映到同一个 $H_q$ 中另一个对块
$(c',d')$，$c'\ne c$（因而 $d'\ne d$）。它们在保持总荷守恒的同时改变每个部分的荷，即在这里把粒子在两部分之间移动。

由定理 7.6 的块公式，$(\rho_1)_c=\sum_d\mathrm{tr}_{l_d}\,\rho_{(c,d),(c,d)}$。本库的两种实现都被设计为
对任意融合规则计算这个映射：`partial_trace` 只对对角的对块 $\rho_{(c,d),(c,d)}$ 求迹，
而 `partial_trace_channel` 对同一扇区中 $c$ 相同的对使用各自独立的 Kraus 指标，使得它的对偶仍然是
$A\mapsto A\otimes1$。代码对下面两种融合规则以及保留因子的两种选择，把两者与用 `einsum` 独立计算的
$E_1(\mathrm{Tr}_{H_2}\,\rho)$ 进行核对，并把对偶与 $A\otimes1$ 进行核对。在我们的例子中，取加法时 $E_1$ 什么也不做；取宇称规则 $f(q,r)=q+r\bmod2$ 时，
`bos1` 的荷 $0$ 和 $2$ 与每个 $r$ 的融合结果都相同，于是 $E_1$ 去掉了真正的相干性：

```python
G = fused_tensor_product(bos1, bos2, fuse=lambda q, r: (q + r) % 2)       # only the parity is conserved
print(G.k_factors)

def block_mask(alg):
    """Ones in the diagonal blocks of (+)_c C^{k_c}: X * block_mask(alg) = sum_c P_c X P_c."""
    return full(alg.operator(torch.ones(1, alg.C, alg.k_max, alg.k_max, dtype=dt)))

for alg, coherent in ((F, False), (G, True)):
    rho = states.random_density_matrix(alg, batch_size=2)
    R = full12(rho).reshape(2, 4, 8, 4, 8)
    for keep, X_H, alg_k in ((1, torch.einsum('xijkj->xik', R), bos1), (2, torch.einsum('xijil->xjl', R), bos2)):
        E_X = X_H * block_mask(alg_k)                                     # E(Tr_H rho), Theorem 7.6
        assert close(full(partial_trace(rho, keep=keep)), E_X)
        assert close(full(composite.partial_trace_channel(alg, keep=keep)(rho)), E_X)
        lost = (X_H - E_X).abs().max().item()                             # coherences removed by E
        print(f"keep={keep}: largest coherence removed: {lost:.3e}")
        assert (lost > 1e-3) if coherent else (lost == 0)
    assert close(composite.partial_trace_channel(alg, keep=1).adjoint('blunt')(A), kron(A, bos2.identity(), alg))
```

<a id="s-7-4"></a>

## 7.4 纯态的纠缠

*讲义：* [7.4 节](../../lecture_notes/zh/07_tensor_products.md#s-7-4)：[引理 7.7](../../lecture_notes/zh/07_tensor_products.md#r-7-7)（Schmidt 分解）和[命题 7.8](../../lecture_notes/zh/07_tensor_products.md#r-7-8)（对称性分辨纠缠）：$E(\psi)=H(p)+\sum_cp_cE(\hat\psi_c)$，其中 $p_c$ 是第一部分的扇区 $c$ 的概率。

$H(p)$ 称为**数纠缠**，$\sum_cp_cE(\hat\psi_c)$ 称为**构型纠缠**。对五个格点上 $2$ 个粒子的随机态，
我们分别由两侧的偏迹以及 $4\times8$ 矩阵 $\Psi$ 的奇异值计算 $E(\psi)$：

```python
s = F.charges.index(2)                                      # the sector of two particles, n_q = 10
n = F.k_factors[s]
psi = torch.randn(n, dtype=dt)
psi = psi / psi.norm()
mat = torch.zeros(1, F.C, F.k_max, F.k_max, dtype=dt)
mat[0, s, :n, :n] = torch.outer(psi, psi.conj())
rho = DensityMatrix(F, matrix=mat)
r1, r2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)

psi_H = torch.zeros(4 * 8, dtype=dt)
psi_H[kron_rows(F)[s]] = psi                                # psi as a vector of H1 (x) H2
lam = torch.linalg.svdvals(psi_H.reshape(4, 8)) ** 2        # squared Schmidt coefficients
lam = lam[lam > 1e-14]
E_psi = -(lam * lam.log()).sum().item()
assert abs(r1.von_neumann_entropy().item() - E_psi) < 1e-10
assert abs(r2.von_neumann_entropy().item() - E_psi) < 1e-10  # Lemma 7.7

def entropy(sigma):
    w = torch.linalg.eigvalsh(sigma).clamp(min=1e-300)
    return -(w * w.log()).sum().item()

p = r1.sector_probabilities()[0]                            # c = 0, 1, 2 particles in the first part
assert torch.allclose(r2.sector_probabilities()[0, [2, 1, 0]], p)   # d(c) = 2 - c
H_num = -(p * p.log()).sum().item()
S_conf = sum(p[c].item() * entropy(r1.matrix[0, c, :k, :k] / p[c]) for c, k in enumerate(bos1.k_factors))
print(f"E = {E_psi:.4f} = H(p) {H_num:.4f} + configurational {S_conf:.4f}")
assert abs(E_psi - H_num - S_conf) < 1e-10                  # Proposition 7.8
```

<a id="s-7-5"></a>

## 7.5 随机态与 Page 公式

*讲义：* [7.5 节](../../lecture_notes/zh/07_tensor_products.md#s-7-5)：[定理 7.11](../../lecture_notes/zh/07_tensor_products.md#r-7-11)（Lubkin）、[推论 7.12](../../lecture_notes/zh/07_tensor_products.md#r-7-12)、[定理 7.13](../../lecture_notes/zh/07_tensor_products.md#r-7-13)（Page 公式）和[命题 7.14](../../lecture_notes/zh/07_tensor_products.md#r-7-14)（诱导测度）。

我们在 $\mathbb C^3\otimes\mathbb C^6$ 的 $10^4$ 个经 `partial_trace` 约化的纯态上以及诱导测度上检验定理 7.11
和 7.13，并对 $M_1\oplus M_2\oplus M_3$、$r=2$ 检验命题 7.14：

```python
a, b, n_samples = 3, 6, 10_000
g = torch.randn(n_samples, a * b, dtype=dt)                 # complex Gaussian: E|g_i|^2 = 1
psi = g / g.norm(dim=-1, keepdim=True)
AB_alg = tensor_product(TypeIAlgebra([a], [a], precision='double'), TypeIAlgebra([b], [b], precision='double'))
pure = DensityMatrix(AB_alg, matrix=(psi[:, :, None] * psi.conj()[:, None, :])[:, None],
                     validate=False)                        # pure by construction
rho_A = partial_trace(pure, keep=1)
purity_exact = (a + b) / (a * b + 1)
page = sum(1 / k for k in range(b + 1, a * b + 1)) - (a - 1) / (2 * b)

def mean_se(x):
    return x.mean().item(), x.std().item() / math.sqrt(x.numel())

for name, rho_s in (('reduced pure states', rho_A),
                    ('induced measure', states.random_density_matrix(rho_A.algebra, 40_000, rank=b))):
    (m_p, e_p), (m_S, e_S) = mean_se(rho_s.purity()), mean_se(rho_s.von_neumann_entropy())
    print(f"{name}: E Tr rho^2 = {m_p:.4f} +- {e_p:.4f} (exact {purity_exact:.4f}), "
          f"E S = {m_S:.4f} +- {e_S:.4f} (Page {page:.4f})")
    assert abs(m_p - purity_exact) < 5 * e_p and abs(m_S - page) < 5 * e_S     # Theorems 7.11, 7.13
    assert math.log((a * b + 1) / (a + b)) < m_S < math.log(a)                # Corollary 7.12
assert math.log(a) - (a - 1) / (a * b + 1) - (a - 1) / (2 * b) <= page <= math.log(a) - (a - 1) / (2 * b)
assert abs(sum(1 / k for k in (3, 4)) - 1 / 4 - 1 / 3) < 1e-15                # two qubits: E S = 1/3

alg3 = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')
m_p, e_p = mean_se(states.random_density_matrix(alg3, 40_000, rank=2).purity())
print(f"E Tr rho^2 = {m_p:.4f} +- {e_p:.4f} (exact {26 / 78:.4f})")
assert abs(m_p - 26 / 78) < 5 * e_p                                          # Proposition 7.14
```

<a id="s-7-6"></a>

## 7.6 可与之比较的工具

在单个 Hilbert 空间上，QuTiP（`Qobj.ptrace`、`entropy_vn`、`rand_ket`、带 Ginibre 分布和秩参数的
`rand_dm`）和 `numpy`（如上所示，一次 reshape 加一次 `einsum`）可以对普通积计算偏迹、熵和诱导态
（上面的 `einsum` 检验就是把这样的计算与本库进行比较）。带守恒荷的张量网络库——ITensor（量子数块稀疏张量）、TeNPy（带荷的
腿、按荷分解的纠缠谱），以及针对非阿贝尔对称性的 TensorKit.jl（融合树）——对大系统以同样的方式
融合扇区，并给出对称性分辨纠缠。本库增加的是代数层面：积是带有扇区和迹的冯·诺依曼代数，
融合规则可以是标签的任意映射，偏迹是一个信道，对每种迹都有其对偶，并且对融合积给出正确的
限制（定理 7.3、7.6），而且这一切都是批量化的。本库没有张量网络，并且 `fused_tensor_product`
接受的规则对每一对荷只给出一个结果，正如阿贝尔对称性那样；非阿贝尔规则，例如 $SU(2)$ 的规则（一对自旋 $j_1,j_2$
组合成若干个总自旋 $|j_1-j_2|,\dots,j_1+j_2$，第 4 章习题 2），在这里不是融合积（本库对自旋 $\frac12$ 链上的 $SU(2)$
另行处理，见第 4 章）。

<a id="x-7"></a>

## 习题

讲义中[习题](../../lecture_notes/zh/07_tensor_products.md#x-7)的数值检验：

1. 检验 $\mathrm{Tr}_2$ 关于 $\tau_{\mathrm{vN}}$ 的对偶
   `partial_trace_channel(alg12).adjoint('tau_vN')` 是 $A\mapsto A\otimes Z$，其中 $Z=\bigoplus_dC_2l_d1_{l_d}$。
2. （次可加性。）在 `alg12` 的随机态上检验互信息
   $S(\rho_1)+S(\rho_2)-S(\rho)$ 是正的，并且等于 $S(\rho\|\rho_1\otimes\rho_2)$。
3. 对两个格点共享的一个粒子 $\psi=(|01\rangle+|10\rangle)/\sqrt2$，检验
   $E(\psi)=H(p)=\log2$。

<details>
<summary>解答</summary>

*习题 1.*

```python
dual_tau = composite.partial_trace_channel(alg12, keep=1).adjoint('tau_vN')
A = random_operator(alg1)
Z = alg2.central([alg2.C * l for l in alg2.k_factors])
assert close(dual_tau(A), kron(A, Z, alg12))
```

*习题 2.*

```python
rho = states.random_density_matrix(alg12, batch_size=4)
rho1, rho2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
mutual = rho1.von_neumann_entropy() + rho2.von_neumann_entropy() - rho.von_neumann_entropy()
assert (mutual > 0).all()
assert torch.allclose(mutual, states.relative_entropy(rho, kron(rho1, rho2, alg12)))
```

*习题 3.*

```python
site = TypeIAlgebra([1, 1], [1, 1], precision='double', charges=[0, 1])
S2 = fused_tensor_product(site, site)
assert S2.k_factors == [1, 2, 1]
mat = torch.zeros(1, 3, 2, 2, dtype=dt)
mat[0, 1] = 0.5                                       # |psi><psi| with psi = (|01> + |10>)/sqrt 2
r1 = partial_trace(DensityMatrix(S2, matrix=mat), keep=1)
assert torch.allclose(r1.sector_probabilities()[0], torch.tensor([0.5, 0.5], dtype=torch.float64))
assert abs(r1.von_neumann_entropy().item() - math.log(2)) < 1e-12
```

</details>

**下一章：** [教程 8](08_dynamics.md)：库中的 Schrödinger 方程与 Lindblad 方程。
