[English](../02_von_neumann_algebras.md) · [Español](../es/02_von_neumann_algebras.md) · [Français](../fr/02_von_neumann_algebras.md) · [Deutsch](../de/02_von_neumann_algebras.md) · **中文** · [日本語](../ja/02_von_neumann_algebras.md) · [Русский](../ru/02_von_neumann_algebras.md)

# 教程 2. 有限维冯·诺依曼代数

本教程配合讲义的[第 2 章](../../lecture_notes/zh/02_von_neumann_algebras.md)。该章证明：$\mathbb C^n$ 上的每个冯·诺依曼代数，
在相差一个酉的基变换的意义下，都是
$\bigoplus_cM_{k_c}(\mathbb C)\otimes1_{m_c}$（[定理 2.6](../../lecture_notes/zh/02_von_neumann_algebras.md#r-2-6)）。
我们以一个对称群为例数值地计算交换子，并看本库如何存储这样的代数、它的中心以及到中心上的条件期望。

<a id="s-2-4"></a>

## 2.4 例子：对称群及其交换子

*讲义：* [2.4 节](../../lecture_notes/zh/02_von_neumann_algebras.md#s-2-4)；[定理 2.2](../../lecture_notes/zh/02_von_neumann_algebras.md#r-2-2)、[定理 2.6](../../lecture_notes/zh/02_von_neumann_algebras.md#r-2-6)、[命题 2.7](../../lecture_notes/zh/02_von_neumann_algebras.md#r-2-7)。

讲义的 [2.4 节](../../lecture_notes/zh/02_von_neumann_algebras.md#s-2-4)分解了置换 $(\mathbb C^2)^{\otimes3}$ 三个张量因子的
$S_3$ 的表示 $\pi$：$\pi(S_3)'\cong M_4\oplus M_2$ 的维数为
$20$，$\pi(S_3)''\cong M_1\oplus M_2$ 的维数为 $5$。下面的代码对这两个维数做数值上的确认，
并确认 $\pi(S_3)'$ 的交换子又是 $\pi(S_3)''$
（[定理 2.2](../../lecture_notes/zh/02_von_neumann_algebras.md#r-2-2)）：集合 $S$ 的交换子是线性映射
$X\mapsto(AX-XA)_{A\in S}$ 的零空间，这里在浮点运算中通过计数大于某个容差的奇异值来计算其维数。

```python
import itertools
import torch

torch.manual_seed(0)
dt = torch.complex128

def permutation_operator(perm, n_qubits=3):
    """Unitary on (C^2)^{(x)3} that permutes the tensor factors."""
    dim = 2 ** n_qubits
    P = torch.zeros(dim, dim, dtype=dt)
    for bits in itertools.product([0, 1], repeat=n_qubits):
        src = int(''.join(map(str, bits)), 2)
        dst = int(''.join(str(bits[perm[i]]) for i in range(n_qubits)), 2)
        P[dst, src] = 1
    return P

group = [permutation_operator(p) for p in itertools.permutations(range(3))]

def commutant_basis(ops, tol=1e-10):
    """Orthonormal basis (Hilbert-Schmidt) of {X : X A = A X for all A in ops}."""
    d = ops[0].shape[0]
    eye = torch.eye(d, dtype=dt)
    # vec(AX - XA) = (A (x) 1 - 1 (x) A^T) vec(X) for row-major vec
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    _, s, Vh = torch.linalg.svd(L)
    null = Vh[(s > tol).sum():].conj().resolve_conj()   # rows span the null space
    return [v.reshape(d, d) for v in null]

def span_dim(ops, tol=1e-10):
    return int((torch.linalg.svdvals(torch.stack([A.reshape(-1) for A in ops])) > tol).sum())

comm = commutant_basis(group)
print(len(comm), span_dim(group), len(commutant_basis(comm)))
assert len(comm) == 20                                   # Step 1: pi(S3)'  = M_4 (+) M_2
assert span_dim(group) == 5                              # Step 2: pi(S3)'' = M_1 (+) M_2
assert len(commutant_basis(comm)) == 5                   # double commutant theorem
```

重数在谱中是可见的（[讲义，2.4 节](../../lecture_notes/zh/02_von_neumann_algebras.md#s-2-4)，由
[定理 2.6](../../lecture_notes/zh/02_von_neumann_algebras.md#r-2-6) 得出）：
每个自伴的 $K\in\pi(S_3)'$ 在对称子空间 $F$ 的补空间上的特征值都具有偶数重数；而对 $K=E(Z+Z^*)$，
其中 $Z$ 是 Ginibre 矩阵、$E$ 是群平均，以概率一，
$K$ 在 $F$ 上有 $4$ 个单重特征值，在其补空间上有 $2$ 个二重简并特征值。代码检验了这些
简并对（它们对每个 $K$ 都成立），并对所抽取的 $K$ 检验两对特征值与 $F$ 上的四个特征值互不相同；
这里利用了到 $F$ 上的正交投影 $P_F=\frac1{|G|}\sum_g\pi(g)$。

```python
def group_average(X):
    return sum(g @ X @ g.conj().T for g in group) / len(group)

H0 = torch.randn(8, 8, dtype=dt)                        # Ginibre: Re, Im independent N(0, 1/2)
Hsym = group_average(H0 + H0.conj().T)                  # K = E(Z + Z^*), a random symmetric Hamiltonian
assert all(torch.allclose(g @ Hsym, Hsym @ g) for g in group)

P_sym = sum(group) / len(group)                         # P_F, the projection onto the symmetric subspace
w, Q = torch.linalg.eigh(P_sym)
Q_rest = Q[:, w.real < 0.5]                             # orthonormal basis of the complement (dim 4)
e_rest = torch.linalg.eigvalsh(Q_rest.conj().T @ Hsym @ Q_rest)
print(e_rest)
assert Q_rest.shape[1] == 4
assert torch.allclose(e_rest[0], e_rest[1]) and torch.allclose(e_rest[2], e_rest[3])   # even multiplicities
e_sym = torch.linalg.eigvalsh(Q[:, w.real > 0.5].conj().T @ Hsym @ Q[:, w.real > 0.5])   # eigenvalues on F
assert e_sym.shape[0] == 4 and e_rest[2] - e_rest[1] > 1e-6                           # this sample: generic
assert (e_sym[1:] - e_sym[:-1]).min() > 1e-6                                           # 4 simple eigenvalues on F
assert (e_sym[:, None] - e_rest[None, :]).abs().min() > 1e-6                           # all six distinct
```

<a id="s-2-5"></a>

## 2.5 库中的代数

*讲义：* [2.3 节](../../lecture_notes/zh/02_von_neumann_algebras.md#s-2-3)、[定理 2.6](../../lecture_notes/zh/02_von_neumann_algebras.md#r-2-6)；中心：[2.2 节](../../lecture_notes/zh/02_von_neumann_algebras.md#s-2-2)、[引理 2.3](../../lecture_notes/zh/02_von_neumann_algebras.md#r-2-3)。

本库在 $\bigoplus_c\mathbb C^{k_c}$ 上存储定理 2.6 中的抽象代数 $\bigoplus_cM_{k_c}(\mathbb C)$，
不带重数：一个算子就是它的块 $A_c$ 组成的列表，各块填充到共同的大小。上例中的交换子就是
`TypeIAlgebra([4, 2], [4, 2])`。

```python
from torch_vn_algebra import TypeIAlgebra, channels

alg = TypeIAlgebra([4, 2], [4, 2], precision='double')    # M_4 (+) M_2
print(alg.k_factors, alg.C)

A = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
print(A.matrix.shape)              # (1, 2, 4, 4): block 1 occupies the top-left 2 x 2 corner
assert A.matrix[0, 1, 2:, :].abs().max() == 0 and A.matrix[0, 1, :, 2:].abs().max() == 0
```

中心由元素 $\sum_cz_cP_c$ 组成；`alg.central` 构造这些元素，而它们与一切元素可交换：

```python
def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

P0, P1 = alg.central([1.0, 0.0]), alg.central([0.0, 1.0])   # minimal central projections
Z = alg.central([2.0, -1.0])
assert P0.is_projection and P1.is_projection
assert close(P0 + P1, alg.identity())
assert close(Z @ A, A @ Z)
assert close(A @ P0 + A @ P1, A)                         # A = sum_c A P_c
```

<a id="s-2-6"></a>

## 2.6 条件期望

*讲义：* [2.6 节](../../lecture_notes/zh/02_von_neumann_algebras.md#s-2-6)。

由[讲义，命题 2.8](../../lecture_notes/zh/02_von_neumann_algebras.md#r-2-8)，$E(A)=\sum_c\frac{\mathrm{Tr}\,A_c}{k_c}\,P_c$
是从 $\bigoplus_cM_{k_c}$ 到其中心上的条件期望：它是正的，在中心上是恒等映射，
是中心上的双模映射，并且保持每个迹 $\sum_cw_c\mathrm{Tr}\,A_c$。

在本库中，这个映射是信道 `channels.center_expectation(alg)`。把它作用于算子，得到其
经典（按扇区平均的）部分：

```python
E = channels.center_expectation(alg)
EA = E(A)
assert close(E(EA), EA)                                  # idempotent
assert close(E(Z @ A @ P1), Z @ EA @ P1)                 # bimodule property
A2 = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
assert close(EA @ A2, A2 @ EA)                           # E(A) commutes with a fresh element
for c, k in enumerate(alg.k_factors):                    # E(A)_c = (Tr A_c / k_c) 1, so Tr E(A)_c = Tr A_c
    trA = A.matrix[0, c].diagonal().sum()
    assert close(EA.matrix[0, c, :k, :k], trA / k * torch.eye(k, dtype=dt))
    assert abs(EA.matrix[0, c].diagonal().sum() - trA) < 1e-10
```

<a id="s-2-7"></a>

## 2.7 可与之比较的工具

2.4 节的例子可以用任何线性代数软件包检验；
表示论软件（GAP、SageMath）能以符号方式给出分解
$\pi\cong\bigoplus\rho_\lambda\otimes1_{m_\lambda}$。物理中的严格对角化代码（QuSpin、QuTiP
的对称性模块）使用同样的块结构，逐个扇区地对角化哈密顿量。本库以这一结构为
出发点：一旦块已知，后续所有计算——函数、迹、态、信道、动力学——都停留在代数内部。

<a id="x-2"></a>

## 习题

讲义的[习题 2](../../lecture_notes/zh/02_von_neumann_algebras.md#x-2)：对 $\mathbb C^3$ 上的循环移位 $S$，
$Se_i=e_{i+1\bmod3}$，$\{S\}'$ 与 $\{S\}''$ 都等于 $\mathrm{span}\{1,S,S^2\}$，即循环矩阵全体。
用数值方法检验其维数。

<details>
<summary>解答</summary>

```python
S = torch.roll(torch.eye(3, dtype=dt), 1, dims=0)
assert len(commutant_basis([S])) == 3
assert len(commutant_basis(commutant_basis([S]))) == 3
assert span_dim([torch.eye(3, dtype=dt), S, S @ S]) == 3
```

</details>

**下一章：** [教程 3. 迹、态与熵](03_traces_and_states.md)
