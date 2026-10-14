[English](../04_symmetries.md) · [Español](../es/04_symmetries.md) · [Français](../fr/04_symmetries.md) · [Deutsch](../de/04_symmetries.md) · **中文** · [日本語](../ja/04_symmetries.md) · [Русский](../ru/04_symmetries.md)

# 教程 4. 对称性与表示

本教程配合讲义的[第 4 章](../../lecture_notes/zh/04_symmetries.md)。该章发展了对称性代数背后的表示论，把它变成一个随机化算法，
用以分解由任意一组矩阵生成的冯·诺依曼代数（[命题 4.9](../../lecture_notes/zh/04_symmetries.md#r-4-9)），将其应用于 $`U(1)`$ 和 $`SU(2)`$，
并把由实矩阵生成的代数的扇区分为实型、复型与四元数型（[定理 4.20](../../lecture_notes/zh/04_symmetries.md#r-4-20)）。我们计算特征标，
对群和对生成元运行分解，用 `SpinChain` 和 `su2.SU2Chain` 逐个扇区地构造自旋链，并计算一个分解的实结构。

<a id="s-4-2"></a>

## 4.2 同型分解与交换子

*讲义：* [4.2 节](../../lecture_notes/zh/04_symmetries.md#s-4-2)；[定理 4.3](../../lecture_notes/zh/04_symmetries.md#r-4-3)、[定理 4.4](../../lecture_notes/zh/04_symmetries.md#r-4-4)、[推论 4.5](../../lecture_notes/zh/04_symmetries.md#r-4-5)。

讲义的 [4.2 节](../../lecture_notes/zh/04_symmetries.md#s-4-2)手工计算了 $`S_3`$ 的三个不可约表示——平凡表示、符号表示以及
$`(1,1,1)^\perp\subseteq\mathbb C^3`$ 上的标准表示——的特征标，以及置换 $`(\mathbb C^2)^{\otimes3}`$ 的三个因子的表示 $`\pi`$ 的特征标
$`\chi_\pi(\sigma)=2^{\#\text{cycles}(\sigma)}`$（见下表），并由[推论 4.5](../../lecture_notes/zh/04_symmetries.md#r-4-5)得到重数 $`m_{\rm triv}=4`$，$`m_{\rm sgn}=0`$，
$`m_{\rm std}=2`$ 以及 $`\dim\pi(S_3)'=20`$。代码重新计算这张表：它把 $`P_\sigma`$ 求为 $`\pi(\sigma)`$ 在
$`|100\rangle,|010\rangle,|001\rangle`$ 张成的不变子空间上的限制（$`\pi(\sigma)`$ 在其上置换唯一的那个 $`1`$ 的位置，引理 B.1(a)），
并且 $`\#\text{cycles}(\sigma)=\dim\ker(P_\sigma-1)`$（在各轮换上为常数的向量）。

| | 恒等元（$`1`$） | 对换（$`3`$） | $`3`$-轮换（$`2`$） |
|---|---|---|---|
| $`\chi_{\rm triv}`$ | $`1`$ | $`1`$ | $`1`$ |
| $`\mathrm{sgn}`$ | $`1`$ | $`-1`$ | $`1`$ |
| $`\chi_{\rm std}`$ | $`2`$ | $`0`$ | $`-1`$ |
| $`\chi_\pi`$ | $`8`$ | $`4`$ | $`2`$ |

```python
import math
import torch
import torch_vn_algebra.decompose as dc

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    return (torch.as_tensor(X) - torch.as_tensor(Y)).abs().max().item() < tol

G = dc.permutation_representation(3, 2)                  # (6, 8, 8): pi(sigma) permutes three qubits
chi_pi = G.diagonal(dim1=-2, dim2=-1).sum(-1).real       # characters Tr pi(sigma)
P3 = G[:, [4, 2, 1]][:, :, [4, 2, 1]].real              # on |100>, |010>, |001>: S_3 permuting C^3
cycles = 3 - torch.linalg.matrix_rank(P3 - torch.eye(3, dtype=torch.float64)).double()   # dim ker(P - 1)
assert close(chi_pi, 2.0 ** cycles)
irreps = {'trivial': torch.ones(6, dtype=torch.float64), 'sign': (-1.0) ** (3 - cycles),
          'standard': P3.diagonal(dim1=-2, dim2=-1).sum(-1) - 1}                 # fixed points - 1
chis = torch.stack(list(irreps.values()))
assert close(chis @ chis.T / 6, torch.eye(3, dtype=torch.float64))    # Corollary 4.5(a): all irreducible
mult = chis @ chi_pi / 6                                               # Corollary 4.5(b)
print(dict(zip(irreps, mult.tolist())), (chi_pi ** 2).mean().item())
assert close(mult, torch.tensor([4.0, 0.0, 2.0], dtype=torch.float64))
assert dc.commutant_dimension(G) == 20                                # Corollary 4.5(c): M_4 (+) M_2
```

这就是 2.4 节的分解，现在由每个表示的一张六个数的表得到。

<a id="s-4-3"></a>

## 4.3 从生成元到代数：一个随机化算法

*讲义：* [4.3 节](../../lecture_notes/zh/04_symmetries.md#s-4-3)；[引理 4.6](../../lecture_notes/zh/04_symmetries.md#r-4-6)、[命题 4.7](../../lecture_notes/zh/04_symmetries.md#r-4-7)、[引理 4.8](../../lecture_notes/zh/04_symmetries.md#r-4-8)、[命题 4.9](../../lecture_notes/zh/04_symmetries.md#r-4-9)。

讲义中的随机化算法就是 `decompose.wedderburn`，其前端为 `from_group`（有限群的全部元素；用 twirl `decompose.twirl`
对交换子采样）和 `from_generators`（任意一组矩阵；交换子是一个零空间，`decompose.commutant_basis`）。它把 $`\mathbb C^n`$
分解为一个随机 $`K\in\mathcal A'_h`$ 的特征子空间，借助一个随机 $`Y\in\mathcal A'`$ 的块对它们分组并对齐，然后返回一个
`Decomposition`，其中含有酉矩阵 `V`、诸 $`k_c`$（`k_factors`）和诸 $`m_c`$（`multiplicities`）。在浮点运算中，它把相差小于
某个容差的特征值视为相等，并把交换子的新随机元素压缩到所得子空间上，以进一步分裂这些子空间
（[引理 4.8(a)](../../lecture_notes/zh/04_symmetries.md#r-4-8)）。

`Decomposition.verify()` 报告生成元（以及交换子的新元素）到 $`\mathcal B`$（相应地 $`\mathcal B'`$）的相对距离，
并将 $`\sum_cm_c^2`$ 与 $`\dim\mathcal A'`$ 比较（对群由特征标得到，推论 4.5(c)；对生成元由零空间得到）。我们
对 $`S_3`$ 与三个量子比特运行它，然后看看内部：一个随机 $`K`$ 的特征子空间和一个随机 $`Y`$ 的缠结算子。

```python
D = dc.from_group(G)                                     # the twirl samples pi(S_3)'
print(D)
assert D.k_factors == (1, 2) and D.multiplicities == (4, 2)            # Section 2.4: k = (1, 2), m = (4, 2)
assert D.verify()['ok']
rho = D.to_algebra(G).matrix                             # (6, 2, 2, 2): the blocks rho_c(sigma), padded
assert close(rho[:, 0, 0, 0], torch.ones(6))                           # sector 1: trivial representation
assert close(rho[:, 1].diagonal(dim1=-2, dim2=-1).sum(-1), chis[2])    # sector 2: the standard one
schur = torch.einsum('gij,gkl->ijkl', rho[:, 1], rho[:, 1].conj()) / 6
assert close(schur.reshape(4, 4), torch.eye(4) / 2)                    # Theorem 4.4 with k = 2
assert close(D.from_algebra(D.to_algebra(G)), G)

Z = torch.randn(8, 8, dtype=dt)
K = dc.twirl(G, Z + Z.mH)                                # a random element of pi(S_3)'_h
w, Q = torch.linalg.eigh(K)
cuts = [0] + [i + 1 for i in range(7) if w[i + 1] - w[i] > 1e-8] + [8]
spaces = [Q[:, a:b] for a, b in zip(cuts[:-1], cuts[1:])]
assert sorted(W.shape[1] for W in spaces) == [1, 1, 1, 1, 2, 2]      # Lemma 4.6: m_c pieces of size k_c
assert all(close(G @ W, W @ (W.mH @ G @ W)) for W in spaces)           # invariant subspaces
Y = dc.twirl(G, torch.randn(8, 8, dtype=dt))
Wa, Wb = [W for W in spaces if W.shape[1] == 2]
s = torch.linalg.svdvals(Wa.mH @ Y @ Wb)                               # Lemma 4.8(b): a multiple of a unitary
assert s[1] > 1e-3 and abs(s[0] - s[1]) < 1e-10
assert all(close(W.mH @ Y @ Wa, 0) for W in spaces if W.shape[1] == 1)   # different summands: zero
```

主要开销在交换子：一个 twirl 样本需要 $`|G|`$ 次 $`n\times n`$ 矩阵乘积。对生成元，本库首先对角化 $`\mathcal A`$ 的
一个厄米元素 $`R=U\mathrm{diag}(\lambda)U^*`$（生成元厄米部分的一个随机实线性组合，加上两个这样的组合的对称化
乘积的一个倍数）。每个 $`X\in\mathcal A'`$ 都与 $`R`$ 可交换，因而与它的谱投影可交换（定理 1.5(d)），所以 $`U^*XU`$
关于 $`R`$ 的特征子空间（维数为 $`d_\lambda`$）是分块对角的；这一点被预先施加，于是零空间问题只有
$`N_R=\sum_\lambda d_\lambda^2`$ 个未知数，而不是 $`n^2`$ 个。把若干特征值归并在一起只会使块变大，所以约化后的问题有
相同的解集 $`\mathcal A'`$。在分解中，$`V^*RV=\bigoplus_cr_c\otimes1_{m_c}`$，其中 $`r_c\in M_{k_c}`$ 厄米；若 $`\lambda`$ 在
$`r_c`$ 中的重数为 $`\mu_c(\lambda)`$，则 $`d_\lambda=\sum_c\mu_c(\lambda)m_c`$，并且由于对整数 $`\mu\ge0`$ 有 $`\mu^2\ge\mu`$，
```math
N_R=\sum_\lambda\Big(\sum_c\mu_c(\lambda)m_c\Big)^2\ge\sum_{\lambda,c}\mu_c(\lambda)^2m_c^2\ge\sum_{\lambda,c}\mu_c(\lambda)m_c^2=\sum_ck_cm_c^2,
```
等号成立当且仅当所有 $`r_c`$ 的特征值都是单重的且两两不同。对 $`L=3,4,5`$ 个自旋的总自旋，下面的代码发现等号
成立（*数值观察*）。即使是最小值也很大：对八个自旋 $`\frac12`$（$`n=256`$），$`\sum_ck_cm_c^2=4900`$，而约化后的
方程组是一个 $`(n^2+N_R)\times N_R`$ 矩阵，有 $`3.5\cdot10^8`$ 个复元素。因此对自旋链，4.6 节使用模块 `su2`，它事先
就知道答案。

```python
for n in (3, 4, 5):
    Dn = dc.from_generators(dc.total_spin(n))
    _, info = dc.commutant_basis(dc.total_spin(n), return_info=True)
    N_min = sum(k * m * m for k, m in zip(Dn.k_factors, Dn.multiplicities))
    assert info['unknowns'] >= N_min                                  # N_R >= sum_c k_c m_c^2
    print(n, info['unknowns'], N_min)                                # observed: equality
```

<a id="s-4-4"></a>

## 4.4 $`U(1)`$ 与荷扇区

*讲义：* [4.4 节](../../lecture_notes/zh/04_symmetries.md#s-4-4)；[命题 4.10](../../lecture_notes/zh/04_symmetries.md#r-4-10)。

对 $`n`$ 个量子比特，荷是数算子 $`\hat N=\sum_jn_j`$，其中 $`n_j=|1\rangle\langle1|`$ 作用在量子比特 $`j`$ 上，保持它守恒的算子构成
$`U(1)'\cong\bigoplus_{w=0}^nM_{\binom nw}(\mathbb C)`$（[讲义，4.4 节](../../lecture_notes/zh/04_symmetries.md#s-4-4)），这就是
`chains.SpinChain(n).algebra`（比特 $`1`$ = 自旋向上，$`\hat N=S^z+\frac n2`$；而在 `decompose.total_spin` 中，自旋向上是
$`\mathbb C^2`$ 的*第一个*基向量，这只是给基重新编号，并不改变扇区）。每个 XXZ 哈密顿量都属于 $`U(1)'`$。代码仅由
$`\hat N`$ 找出各扇区，并逐个扇区对角化一个四自旋的 XXZ 环。

```python
import functools
from torch_vn_algebra import SpinChain

L = 4
Sx, Sy, Sz = dc.total_spin(L)                              # (16, 16) each, site 0 most significant
N_op = Sz + L / 2 * torch.eye(2 ** L, dtype=dt)            # number of up spins
DU = dc.from_generators(N_op)                              # the algebra generated by U(1)
print(DU)
chain = SpinChain(L, boundary='periodic')                  # U(1)' = (+)_w M_binom(L, w)
assert set(DU.k_factors) == {1} and sorted(DU.multiplicities) == sorted(chain.algebra.k_factors)
assert sorted(chain.algebra.k_factors) == sorted(math.comb(L, w) for w in range(L + 1))

sx, sy, sz = dc.spin_matrices(0.5)
def two_site(a, b, i, j, n=L):                             # a on site i, b on site j of n sites
    ops = [a if k == i else b if k == j else torch.eye(2, dtype=dt) for k in range(n)]
    return functools.reduce(torch.kron, ops)

H = sum(two_site(sx, sx, i, j) + two_site(sy, sy, i, j) + 0.5 * two_site(sz, sz, i, j)
        for i, j in chain.bonds)                           # XXZ ring, Delta = 0.5
assert close(H @ N_op, N_op @ H)                           # [H, N] = 0 ...
assert close(DU.from_commutant(DU.to_commutant(H)), H)     # ... so H is block diagonal (Proposition 4.10)
Hc = chain.xxz(J=1.0, Delta=0.5)                           # the same H, built sector by sector
E_sec = torch.cat([torch.linalg.eigvalsh(Hc.matrix[0, N, :k, :k]) for N, k in enumerate(chain.algebra.k_factors)])
assert close(torch.sort(E_sec).values, torch.linalg.eigvalsh(H))
```

<a id="s-4-5"></a>

## 4.5 $`SU(2)`$：自旋

*讲义：* [4.5 节](../../lecture_notes/zh/04_symmetries.md#s-4-5)；[定理 4.11](../../lecture_notes/zh/04_symmetries.md#r-4-11)、[命题 4.13](../../lecture_notes/zh/04_symmetries.md#r-4-13)、[定理 4.14](../../lecture_notes/zh/04_symmetries.md#r-4-14)。自旋 $`s`$ 即 `decompose.spin_matrices(s)`，$`L`$ 个自旋 $`\frac12`$ 的总自旋即 `decompose.total_spin(L)`。

对 $`L=4`$，[讲义，命题 4.13](../../lecture_notes/zh/04_symmetries.md#r-4-13) 给出 $`\{S_x,S_y,S_z\}''\cong M_5\oplus M_3\oplus M_1`$，重数为 $`1,3,2`$，
维数为 $`35`$，而交换子 $`\{S_x,S_y,S_z\}'\cong M_1\oplus M_3\oplus M_2`$，维数为 $`14`$。本库仅由三个生成元就得到同样的结果
（扇区按 $`(k_c,m_c)`$ 排序），并且在自旋 $`S`$ 上 $`\mathbf S^2=S(S+1)`$（[定理 4.11](../../lecture_notes/zh/04_symmetries.md#r-4-11)）：

```python
S4 = dc.total_spin(4)
assert close(S4[0] @ S4[1] - S4[1] @ S4[0], 1j * S4[2])                # [S_x, S_y] = i S_z
D4 = dc.from_generators(S4)
print(D4)
assert D4.k_factors == (1, 3, 5) and D4.multiplicities == (2, 3, 1)     # M_1, M_3, M_5 with m = 2, 3, 1
assert D4.verify()['ok'] and D4.dims == (35, 14)
m_S = {S: math.comb(4, 2 - S) - (math.comb(4, 1 - S) if S < 2 else 0) for S in (0, 1, 2)}
assert [m_S[(k - 1) // 2] for k in D4.k_factors] == list(D4.multiplicities)    # Proposition 4.13
casimir = D4.to_algebra((S4 @ S4).sum(0))                                  # S^2 = S(S+1) on spin S
for c, k in enumerate(D4.k_factors):
    assert close(casimir.matrix[0, c, :k, :k], (k - 1) / 2 * (k + 1) / 2 * torch.eye(k))
```

由[讲义，定理 4.14](../../lecture_notes/zh/04_symmetries.md#r-4-14)（Schur–Weyl 对偶），格点的置换张成总自旋的交换子：
$`\mathrm{span}\,\pi(S_L)=\{S_x,S_y,S_z\}'`$。对 $`L\le5`$，代码数值地比较这两个空间：由定理的 (a) 部分，
$`\mathrm{span}\,\pi(S_L)\subseteq\{S_x,S_y,S_z\}'`$，只需它们维数相同。

```python
for n in range(2, 6):
    comm = dc.commutant_basis(dc.total_spin(n)).reshape(-1, 4 ** n)          # a basis of {S_a}'
    perm = dc.permutation_representation(n, 2).reshape(-1, 4 ** n)          # all n! permutations
    dims = [torch.linalg.matrix_rank(A).item() for A in (comm, perm, torch.cat([comm, perm]))]
    print(n, dims)
    assert dims[0] == dims[1] == dims[2]                                     # span pi(S_n) = {S_a}'
Dp = dc.from_group(dc.permutation_representation(4, 2))                   # pi(S_4)'' = {S_a}'
print(Dp)
assert sorted(zip(Dp.multiplicities, Dp.k_factors)) == sorted(zip(D4.k_factors, D4.multiplicities))
```

打印出的维数是 $`2,5,14,42`$，即 Catalan 数 $`\binom{2L}L/(L+1)`$；推论 B.10 对每个 $`L`$ 证明了
$`\dim\{S_a\}'=\sum_Sm_S^2=\binom{2L}L/(L+1)`$。最后一行展示了对偶的作用：$`\pi(S_4)''`$
的块为 $`k=(1,2,3)`$，重数为 $`(5,1,3)`$——正是总自旋的数对 $`(k_c,m_c)`$ 互换角色。一般地（推论 B.9），自旋 $`S`$ 的
重数空间承载 $`S_L`$ 的一个维数为 $`m_S`$ 的不可约表示，并且这些表示两两不等价；这里对 $`S=2,1,0`$ 其维数分别为 $`1,3,2`$。

<a id="s-4-6"></a>

## 4.6 自旋系统

*讲义：* [4.6 节](../../lecture_notes/zh/04_symmetries.md#s-4-6)；[命题 4.15](../../lecture_notes/zh/04_symmetries.md#r-4-15)、[定理 4.16](../../lecture_notes/zh/04_symmetries.md#r-4-16)、[推论 4.17](../../lecture_notes/zh/04_symmetries.md#r-4-17)。

由[讲义，推论 4.17](../../lecture_notes/zh/04_symmetries.md#r-4-17)，$`L`$ 个自旋 $`\frac12`$ 的最高权向量的规范正交基可以用
[定理 4.16](../../lecture_notes/zh/04_symmetries.md#r-4-16)的 Clebsch–Gordan 系数逐个格点地构造。用中间自旋序列 $`S_1=\frac12,S_2,\dots,S_L=S`$
标记每个基向量，就得到 `su2.SU2Chain` 所构造的 Yamanouchi 基（其源代码恰好使用这些系数；下面的代码只检验推论 4.17
的结论），它位于 `SpinChain` 的 $`U(1)`$ 扇区 $`\hat N=\frac L2+S`$ 之中（`SpinChain` 的格点顺序是张量因子的重新标记，由
引理 B.1(b)，它与总自旋可交换）。代码用稠密矩阵以及本库的 `su2.clebsch_gordan` 检验定理 4.16，然后对八个自旋的基
检验推论 4.17 的结论。我们还用 `D4` 在四个自旋的环上检验[命题 4.15(a)](../../lecture_notes/zh/04_symmetries.md#r-4-15)，并比较八个自旋
时的块大小。

```python
from torch_vn_algebra import su2

def ket(j, m):                                            # |j, m> of Theorem 4.11 (zero if |m| > j)
    v = torch.zeros(int(2 * j + 1), dtype=dt)
    if abs(m) <= j:
        v[int(j - m)] = 1
    return v

up, dn = ket(0.5, 0.5), ket(0.5, -0.5)
for j in (0.5, 1.0, 1.5, 2.0):                            # Theorem 4.16
    Sj = dc.spin_matrices(j)
    sp_j = [math.sqrt(j * (j + 1) - (j - k) * (j - k + 1)) for k in range(1, int(2 * j) + 1)]
    assert close(Sj[0] + 1j * Sj[1], torch.diag(torch.tensor(sp_j, dtype=dt), 1))     # S_+ of Theorem 4.11
    assert close(Sj[2], torch.diag(torch.tensor([j - k for k in range(int(2 * j) + 1)], dtype=dt)))
    Sc = torch.stack([torch.kron(a, torch.eye(2, dtype=dt)) + torch.kron(torch.eye(len(a), dtype=dt), b)
                      for a, b in zip(Sj, dc.spin_matrices(0.5))])
    Sp, Sm = Sc[0] + 1j * Sc[1], Sc[0] - 1j * Sc[1]
    cols = []
    for J in (j + 0.5, j - 0.5):
        for M in [J - k for k in range(int(2 * J + 1))]:
            al, be = math.sqrt((j + M + 0.5) / (2 * j + 1)), math.sqrt((j - M + 0.5) / (2 * j + 1))
            p, q = (al, be) if J > j else (-be, al)
            x = p * torch.kron(ket(j, M - 0.5), up) + q * torch.kron(ket(j, M + 0.5), dn)
            assert close(Sc[2] @ x, M * x)
            if M == J:
                assert close(Sp @ x, 0)                               # highest weight
            else:
                assert close(Sm @ cols[-1], math.sqrt((J + M + 1) * (J - M)) * x)   # lowering, (c)
            for s2, e in ((0.5, up), (-0.5, dn)):                     # the library's Racah formula agrees
                cg = su2.clebsch_gordan(j, M - s2, 0.5, s2, J, M) if abs(M - s2) <= j else 0.0
                assert abs(cg - torch.vdot(torch.kron(ket(j, M - s2), e), x).real.item()) < 1e-12
            cols.append(x)
    U = torch.stack(cols, 1)
    assert U.shape == (2 * len(Sj[0]), 2 * len(Sj[0])) and close(U.mH @ U, torch.eye(U.shape[1]))   # (a)

def swap(i, j, n):
    p = list(range(n))
    p[i], p[j] = j, i
    return dc.permutation_representation(n, 2, [p])[0]

ring = [(i, (i + 1) % 4) for i in range(4)]
H4 = sum(swap(i, j, 4) / 2 - torch.eye(16, dtype=dt) / 4 for i, j in ring)        # Dirac's identity
assert close(H4, sum(two_site(a, a, i, j) for a in (sx, sy, sz) for i, j in ring))
h = D4.to_commutant(H4)                                   # blocks H_S of sizes m_S = 2, 3, 1
assert close(D4.from_commutant(h), H4)                    # H4 lies in {S_a}'
levels = torch.cat([torch.linalg.eigvalsh(h.matrix[0, c, :m, :m]).repeat_interleave(k)
                    for c, (k, m) in enumerate(zip(D4.k_factors, D4.multiplicities))])
assert close(torch.sort(levels).values, torch.linalg.eigvalsh(H4))   # level of H_S repeated 2S+1 times

u1, sym = SpinChain(8, boundary='periodic'), su2.SU2Chain(8, boundary='periodic')
for S, m in zip(sym.spins, sym.multiplicities):           # Corollary 4.17: orthonormal highest-weight bases
    B, N = sym.hw_basis(S), int(4 + S)                     # columns in the U(1) sector N = L/2 + S
    assert B.shape == (math.comb(8, N), m) and close(B.T @ B, torch.eye(m, dtype=B.dtype))
    assert N == 8 or close(sym.raise_(B, N), 0)            # S_+ B = 0 (S_z = S holds in sector N)
H_u1, H_su2 = u1.xxz(J=1.0, Delta=1.0), sym.heisenberg(J=1.0)
E_u1 = torch.sort(torch.cat([torch.linalg.eigvalsh(H_u1.matrix[0, N, :k, :k])
                             for N, k in enumerate(u1.algebra.k_factors)])).values
E_su2, S_lev, _ = sym.full_spectrum(H_su2, return_quantum_numbers=True)   # level of block S, 2S+1 times
E_su2, S_lev = E_su2[0], S_lev[0]
assert close(E_u1, E_su2)
assert S_lev[0] == 0 and (S_lev[1:4] == 1).all() and E_su2[4] - E_su2[3] > 1e-3     # singlet, then triplet
assert E_su2[1] - E_su2[0] > 1e-3 and close(E_su2[1:4], E_su2[1].expand(3))       # simple singlet (Theorem 9.9)
assert close(sym.restrict(H_u1).matrix, H_su2.matrix)      # Proposition 4.15(b)
flops = lambda ks: sum(k ** 3 for k in ks)
print(u1.algebra.k_factors, sym.multiplicities, flops(u1.algebra.k_factors), flops(sym.multiplicities))
print(f"ground state energy {E_su2[0].item():.6f}, first excited {E_su2[1].item():.6f} (x3)")
```

对八个自旋的环，$`U(1)`$ 块的大小为 $`\binom8w=1,8,28,56,70,56,28,8,1`$（4.4 节），$`SU(2)`$ 块的大小为
$`m_S=14,28,20,7,1`$（$`S=0,\dots,4`$，命题 4.13）。对角化的开销 $`\sum k^3`$ 从 $`739162`$ 降到 $`33040`$，减小为约
$`1/22`$；代数的维数 $`\sum k^2`$ 从 $`\sum_w\binom8w\binom8{8-w}=\binom{16}8=12870`$（从 $`16`$ 个元素中选 $`8`$ 个，其中 $`w`$ 个
取自前 $`8`$ 个）降到 $`\sum_Sm_S^2=1430=\binom{16}8/9`$（推论 B.10）。由定理 9.9（环的长度为偶数且 $`J>0`$），基态是
唯一的单态，且 $`E_0(0)<E_0(1)`$。*数值观察*（由上面的代码计算）：基态能量为 $`-3.651093`$，第一激发能级是位于
$`-3.128419`$ 的三重态。

<a id="s-4-7"></a>

## 4.7 实代数：实型、复型与四元数型

*讲义：* [4.7 节](../../lecture_notes/zh/04_symmetries.md#s-4-7)。对在复共轭下封闭的 $`*`$-代数 $`\mathcal A`$，以及其列块 $`V_c`$ 属于各扇区的
分解 $`V`$，[定理 4.20](../../lecture_notes/zh/04_symmetries.md#r-4-20) 给出满足 $`\bar P_c=P_{\sigma(c)}`$ 的扇区置换 $`\sigma`$、满足
$`\overline{V_c}=V_{\sigma(c)}(J_c\otimes M_c)`$ 的酉矩阵 $`J_c`$、$`M_c`$，以及当 $`\sigma(c)=c`$ 时满足 $`J_c\bar J_c=\varepsilon_c1`$ 的符号
$`\varepsilon_c`$；类型与实代数 $`\mathcal A_{\mathbb R}`$ 见[定理 4.21](../../lecture_notes/zh/04_symmetries.md#r-4-21)，实 Schur 引理是
[定理 4.22](../../lecture_notes/zh/04_symmetries.md#r-4-22)，交换子的实元素由[命题 4.23](../../lecture_notes/zh/04_symmetries.md#r-4-23)刻画。

`Decomposition.real_structure()` 对 `wedderburn`（`from_generators`、`from_group`）返回的分解计算 $`\sigma`$、$`J_c`$、
$`M_c`$ 以及类型。它检查所存储的生成元的共轭落在 $`V(\bigoplus_cM_{k_c}\otimes1)V^*`$ 中（相对残差
`'conjugation'`），并构造酉矩阵 $`P=V^*\bar V`$，其第 $`c`$ 个列块为 $`V^*\overline{V_c}`$：由定理 4.20(a)，它只在扇区
$`\sigma(c)`$ 的行上非零，并在那里等于 $`J_c\otimes M_c`$。本库取承载该块最大部分的扇区作为 $`\sigma(c)`$，要求
$`(k_{\sigma(c)},m_{\sigma(c)})=(k_c,m_c)`$，并界定其余部分（`'leak'`）；它把该块重排为
$`k_c^2\times m_c^2`$ 矩阵，其元素为 $`T_{(a,b),(j,l)}=T_{(a,j),(b,l)}`$，对 $`T=J\otimes M`$ 这就是秩一矩阵
$`\mathrm{vec}\,J\,(\mathrm{vec}\,M)^T`$，由其主奇异向量读出 $`J_c`$ 与 $`M_c`$（取
$`\|J_c\|_F^2=k_c`$），并记录该块与 $`J_c\otimes M_c`$ 的距离（`'kron'`）。若这三个残差之一超过容差、形状不同，或
配对不构成对合，它就抛出 `ValueError`（该代数在共轭下不封闭）。对 $`\sigma(c)=c`$，它把 $`\varepsilon_c`$ 设为
$`\mathrm{Re}\,\mathrm{Tr}(J_c\bar J_c)/k_c`$ 的符号——在精确算术中它恰为 $`\pm1`$，所以两个候选值相差 $`2`$——并把
这些扇区上 $`\|J_c\bar J_c-\varepsilon_c\|_F/\sqrt{k_c}`$ 与 $`\|M_c\bar M_c-\varepsilon_c\|_F/\sqrt{m_c}`$ 的最大值记为
`'sign'`；若某个 $`\varepsilon_c=-1`$ 而 $`k_c`$ 或 $`m_c`$ 为奇数（由定理 4.20(b) 这不可能），则把它提高到至少 $`1`$。结果
具有 `types`、`partner`（$`\sigma`$）、`J`、`M`、`residuals` 和 `ok`，后者表示所有残差都不超过容差（默认是分解的
接受阈值）；此时 $`|\mathrm{Tr}(J_c\bar J_c)/k_c-\varepsilon_c|\le\|J_c\bar J_c-\varepsilon_c\|_F/\sqrt{k_c}`$
（Cauchy–Schwarz）也不超过容差。

`Decomposition.real_basis()` 返回一个实正交 $`n\times n`$ 矩阵 $`O`$，以及描述其各列块的列表（类型、扇区、
列范围）：对实型扇区，是定理 4.21(a) 的列 $`O_c=V_c(W_J\otimes W_M)`$，使得 $`O_c^TaO_c=\alpha_c(a)\otimes1_{m_c}`$，
其中 $`\alpha_c(a)`$ 为实矩阵；对四元数型扇区，是列 $`V_cW`$，其中 $`J_c\otimes M_c=WW^T`$ 且
$`(J_c\otimes M_c)\bar W=W`$（引理 4.19(c)，取 $`\varepsilon=1`$，因为 $`(J_c\otimes M_c)\overline{(J_c\otimes M_c)}=\varepsilon_c^2=1`$），
由定理 4.21(a) 的计算，它们是实的；对一对 $`c,\sigma(c)`$（$`c<\sigma(c)`$），是列 $`\sqrt2\mathrm{Re}\,V_c`$，后接
$`\sqrt2\mathrm{Im}\,V_c`$。这些列是规范正交的：$`V_c^*\overline{V_c}=V_c^*V_{\sigma(c)}(J_c\otimes M_c)=0`$ 与
$`V_c^T\overline{V_c}=\overline{V_c^*V_c}=1`$ 及其共轭给出 $`(\mathrm{Re}\,V_c)^T\mathrm{Re}\,V_c=(\mathrm{Im}\,V_c)^T\mathrm{Im}\,V_c=\frac12`$ 与
$`(\mathrm{Re}\,V_c)^T\mathrm{Im}\,V_c=0`$；不同块的列落在 $`\sigma`$ 的不同轨道 $`\mathcal O`$ 的正交投影
$`P_{\mathcal O}=\sum_{e\in\mathcal O}P_e`$ 的像中（因为 $`\overline{V_c}=V_{\sigma(c)}(J_c\otimes M_c)`$），所以 $`O`$ 是正交的。对实的 $`a`$，
取 $`aV_c=V_cx`$（$`x=\pi_c(a)\otimes1_{m_c}`$）的实部与虚部，得
```math
O^TaO=\begin{pmatrix}\mathrm{Re}\,x&\mathrm{Im}\,x\\-\mathrm{Im}\,x&\mathrm{Re}\,x\end{pmatrix}\quad\text{on these columns.}
```
对称酉矩阵 $`J=X+iY`$ 的满足 $`WW^T=J`$ 的 Takagi 因子 $`W`$，由 $`X+tY`$（$`t`$ 随机）的正交对角化计算得到
（$`X`$ 与 $`Y`$ 是可交换的实对称矩阵，因为 $`J\bar J=JJ^*=1`$ 的虚部是 $`YX-XY`$）；本库依据残差
$`\|WW^T-J\|_F/\|J\|_F`$ 与 $`\|\mathrm{Im}(V_cW)\|_F/\sqrt{k_cm_c}`$ 接受结果，若其中之一超过
$`\max(\text{tolerance},10^3\,\epsilon_{\rm mach}\,n)`$，则发出警告（`DecompositionWarning`）。

讲义的 [4.7 节](../../lecture_notes/zh/04_symmetries.md#s-4-7)详细讨论了三个由实矩阵生成的代数：(i) $`\mathbb R^3\otimes\mathbb R^2`$ 上的
$`L_a\otimes1_2`$，其中 $`(L_a)_{bc}=-\epsilon_{abc}`$ 是 $`\mathbb R^3`$ 的转动的生成元：一个实型扇区，$`k=3`$，$`m=2`$；(ii)
$`\mathbb R^2\otimes\mathbb R^2`$ 上的 $`\Omega_2\otimes1_2`$：一对复型扇区，$`k=1`$，$`m=2`$；(iii) 四元数在
$`\mathbb H^2=\mathbb R^8`$ 上的左乘 $`1_2\otimes L_q`$：一个四元数型扇区，$`k=2`$，$`m=4`$。代码计算每个代数的分解与实结构，
并检验定理 4.20(a),(b)，以及实型（定理 4.21(a)）和复型扇区对的 $`O^TaO`$ 的分块形式。

```python
def qmul(p, q):                                           # H = CD(C), Section E.6: (a, b)(c, d)
    (a1, b1), (a2, b2) = p, q
    return (a1 * a2 - b2.conjugate() * b1, b2 * a1 + b1 * a2.conjugate())

def q_left(q):                                            # left multiplication by q on H = R^4
    cols = []                                             # coordinates (Re a, Im a, Re b, Im b)
    for e in torch.eye(4, dtype=torch.float64).tolist():
        u, v = qmul(q, (complex(e[0], e[1]), complex(e[2], e[3])))
        cols.append([u.real, u.imag, v.real, v.imag])
    return torch.tensor(cols, dtype=torch.float64).T

L_so3 = torch.zeros(3, 3, 3, dtype=torch.float64)         # (L_a)_bc = -epsilon_abc
for i, (b1, b2) in enumerate([(1, 2), (2, 0), (0, 1)]):
    L_so3[i, b1, b2], L_so3[i, b2, b1] = -1.0, 1.0
R_rot = torch.tensor([[[0.0, -1.0], [1.0, 0.0]]])         # Omega_2
L_H = torch.stack([q_left(u) for u in ((1j, 0), (0, 1), (0, 1j))])            # L_I, L_J, L_K
assert close(L_H[0] @ L_H[1], L_H[2]) and close(q_left((2 - 1j, 3j)).T, q_left((2 + 1j, -3j)))   # IJ = K, L_qbar = L_q^T
eye2 = torch.eye(2, dtype=torch.float64)
real_examples = {'R': (torch.kron(L_so3, eye2[None]), (3,), (2,)),            # L_a (x) 1_2 on R^3 (x) R^2
                 'C': (torch.kron(R_rot, eye2[None]), (1, 1), (2, 2)),        # Omega_2 (x) 1_2
                 'H': (torch.kron(eye2[None], L_H), (2,), (4,))}              # L_q on H^2 = R^8
dec_real = {}
for t, (gens, ks, ms) in real_examples.items():
    D_t = dc.from_generators(gens)
    rs_t, (O_t, blocks_t) = D_t.real_structure(), D_t.real_basis()
    print(t, D_t, rs_t)
    assert D_t.verify()['ok'] and rs_t.ok and (D_t.k_factors, D_t.multiplicities) == (ks, ms)
    assert set(rs_t.types) == {t} and close(O_t.T @ O_t, torch.eye(D_t.n))     # real orthogonal basis
    for c, d in enumerate(rs_t.partner):                                      # Theorem 4.20
        o_c, o_d, w = D_t.offsets[c], D_t.offsets[d], ks[c] * ms[c]
        assert close(D_t.V[:, o_c:o_c + w].conj(), D_t.V[:, o_d:o_d + w] @ torch.kron(rs_t.J[c], rs_t.M[c]))
        if d == c:                                                            # J conj(J) = M conj(M) = +-1
            eps = 1.0 if t == 'R' else -1.0
            assert close(rs_t.J[c] @ rs_t.J[c].conj(), eps * torch.eye(ks[c]))
            assert close(rs_t.M[c] @ rs_t.M[c].conj(), eps * torch.eye(ms[c]))
    dec_real[t] = D_t, rs_t, O_t

D_t, rs_t, O_t = dec_real['R']
for g in real_examples['R'][0]:                           # Theorem 4.21(a): O^T a O = alpha (x) 1_2, alpha real
    G_t = O_t.T @ g @ O_t
    assert close(G_t, torch.kron(G_t[0::2, 0::2], eye2))
D_t, rs_t, O_t = dec_real['C']                            # the pair: [[Re x, Im x], [-Im x, Re x]]
x_t = torch.kron(D_t.to_algebra(real_examples['C'][0]).matrix[0, 0, :1, :1], torch.eye(2, dtype=dt))
assert close(O_t.T @ real_examples['C'][0][0] @ O_t,
             torch.cat([torch.cat([x_t.real, x_t.imag], 1), torch.cat([-x_t.imag, x_t.real], 1)]))
```

代码在三个例子的不可约实模 $`\mathbb R^3`$、$`\mathbb R^2`$、$`\mathbb R^4`$ 上验证定理 4.22（本库找到的类型与重数——
$`m=1`$ 的实型、$`m=1`$ 的复型对、$`m=2`$ 的四元数型——由定理 4.22(c) 使它们不可约），其中实交换子作为实线性方程
$`Xs=sX`$，$`Xs^T=s^TX`$ 的零空间来计算；然后对每个例子中交换子的一个随机实对称元素验证命题 4.23。

```python
def real_commutant_dim(gens):                             # dim_R {X real : X s = s X, X s^T = s^T X}
    gens = torch.cat([gens, gens.mT])
    eye_n = torch.eye(gens.shape[-1], dtype=torch.float64)
    eqs = torch.cat([torch.kron(eye_n, s) - torch.kron(s.T.contiguous(), eye_n) for s in gens])
    return gens.shape[-1] ** 2 - torch.linalg.matrix_rank(eqs).item()

for gens, types, ms, dim_D in ((L_so3, ('R',), (1,), 1), (R_rot, ('C', 'C'), (1, 1), 2), (L_H, ('H',), (2,), 4)):
    D_irr = dc.from_generators(gens)                      # Theorem 4.22(c): R^3, R^2, R^4 are irreducible
    assert D_irr.real_structure().types == types and D_irr.multiplicities == ms
    assert real_commutant_dim(gens) == dim_D              # Theorem 4.22(b): D_U = R, C, H

gen_X = torch.Generator().manual_seed(1)
for t, (gens, ks, ms) in real_examples.items():           # Proposition 4.23
    D_t, rs_t, O_t = dec_real[t]
    B_t = dc.commutant_basis(gens)
    X_t = torch.einsum('r,rij->ij', torch.randn(len(B_t), dtype=dt, generator=gen_X), B_t).real
    X_t = X_t + X_t.T                                     # real symmetric, in A' (Lemma 4.18(b))
    assert close(X_t @ gens, gens @ X_t)
    x_t = D_t.to_commutant(X_t).matrix[0]
    ev = [torch.linalg.eigvalsh(x_t[c, :m, :m]) for c, m in enumerate(ms)]
    if t == 'R':                                          # (a): O^T X O = 1_3 (x) y, y real symmetric
        Y_t = O_t.T @ X_t @ O_t
        assert close(Y_t, torch.kron(torch.eye(3, dtype=torch.float64), Y_t[:2, :2]))
    if t == 'C':                                          # (c): the two sectors have the same spectrum
        assert close(ev[0], ev[1])
    if t == 'H':                                          # (b): Kramers pairs, each level 2 k_c = 4 times in X
        assert close(ev[0][0::2], ev[0][1::2])
        assert close(torch.linalg.eigvalsh(X_t), ev[0].repeat_interleave(2))
        print('quaternionic block x_c:', [round(v, 4) for v in ev[0].tolist()])
```

<a id="s-4-8"></a>

## 4.8 可与之比较的工具

对有限群，GAP 和 SageMath 精确计算特征标表，然后由推论 4.5 得到重数；RepLAB（MATLAB；D. Rosset,
F. Montealegre-Mora, J.-D. Bancal）数值地分解有限群和紧群的表示，同样使用交换子的随机元素；
SymbolicWedderburn.jl（Julia）为半定规划的对称性约化计算 Wedderburn 分解。在物理中，QuSpin 和其他精确
对角化程序对 $`U(1)`$ 和晶格对称性进行块对角化，ITensor 和 TeNPy 在张量网络中携带阿贝尔荷，TensorKit.jl
用融合树处理 $`SU(2)`$。本库的 `decompose` 是面向小矩阵（$`n`$ 至多几百）的稠密方法；它接受任意矩阵列表
（群元素或生成元），验证其结果（命题 4.9），并把代数及其交换子作为带有 `to_algebra`、`from_algebra`、
`to_commutant`、`from_commutant` 的 `TypeIAlgebra` 对象返回，使第 3 章、第 5–8 章的内容可以直接应用。对
自旋 $`\frac12`$ 链，`su2.SU2Chain` 直接给出 $`SU(2)`$ 块（推论 4.17）；它最多接受 $`L=18`$ 个自旋，这一上限由基所需的内存决定。

<a id="x-4"></a>

## 习题

讲义的[习题](../../lecture_notes/zh/04_symmetries.md#x-4)的数值检验（下面的解答按此顺序排列）：
1. 对 $`d=2`$ 和 $`L=3,4`$ 验证 $`\dim\pi(S_L)'=\binom{d^2+L-1}L`$（即 $`20`$ 和 $`35`$）；
2. 用 `from_generators` 验证 $`1\otimes\frac12\cong\frac12\oplus\frac32`$ 和 $`1\otimes1\cong0\oplus1\oplus2`$；
3. 四个量子比特上 $`e^{i\theta\hat N}`$ 的交换子，取 $`\theta=1`$（荷 $`0,\dots,4`$ 给出互不相同的 $`e^{iq}`$，因为
   $`0<|q-q'|\le4<2\pi`$）和 $`\theta=\pi`$；
4. Zeeman 分裂：八个自旋的 Heisenberg 环的 $`H+hS_z`$ 的谱，分别用 `su2.SU2Chain` 和在各 $`U(1)`$ 扇区中计算；
5. $`Q_8`$ 的正则表示的扇区，用 `real_structure` 计算。

<details>
<summary>解答</summary>

```python
for n, value in ((3, 20), (4, 35)):
    assert dc.commutant_dimension(dc.permutation_representation(n, 2)) == math.comb(n + 3, n) == value
```

```python
def coupled(a, b):                                        # S^(a) (x) 1 + 1 (x) S^(b)
    Sa, Sb = dc.spin_matrices(a), dc.spin_matrices(b)
    ea, eb = torch.eye(Sa.shape[-1], dtype=dt), torch.eye(Sb.shape[-1], dtype=dt)
    return torch.stack([torch.kron(x, eb) + torch.kron(ea, y) for x, y in zip(Sa, Sb)])
for a, b, ks in ((1.0, 0.5, (2, 4)), (1.0, 1.0, (1, 3, 5))):  # 1 (x) 1/2 = 1/2 (+) 3/2, 1 (x) 1 = 0 (+) 1 (+) 2
    D_cg = dc.from_generators(coupled(a, b))
    assert D_cg.k_factors == ks and set(D_cg.multiplicities) == {1}
```

```python
for theta, m in ((1.0, (1, 1, 4, 4, 6)), (math.pi, (8, 8))):   # theta = 1: e^{iq}, q = 0..4, distinct as 4 < 2 pi
    D_theta = dc.from_generators(torch.linalg.matrix_exp(1j * theta * N_op))
    assert D_theta.multiplicities == m and set(D_theta.k_factors) == {1}
```

```python
assert not close(H4 @ Sx - Sx @ H4 + 0.3 * (Sz @ Sx - Sx @ Sz), 0)          # [H4 + 0.3 S_z, S_x] != 0
E_field = sym.full_spectrum(sym.heisenberg(J=1.0, h=0.3))[0]                # E + h M on every multiplet
H_field = u1.xxz(J=1.0, Delta=1.0, h=0.3)
E_u1_field = torch.cat([torch.linalg.eigvalsh(H_field.matrix[0, N, :k, :k])
                        for N, k in enumerate(u1.algebra.k_factors)])
assert close(torch.sort(E_u1_field).values, E_field)
```

```python
Q8 = [(s * a, s * b) for a, b in ((1, 0), (1j, 0), (0, 1), (0, 1j)) for s in (1, -1)]
key = lambda q: tuple(round(v) for z in q for v in (complex(z).real, complex(z).imag))
pos = {key(q): i for i, q in enumerate(Q8)}
reg = torch.zeros(8, 8, 8, dtype=torch.float64)           # pi(g) e_h = e_{gh}
for i, g in enumerate(Q8):
    for j, h in enumerate(Q8):
        reg[i, pos[key(qmul(g, h))], j] = 1.0
D_Q8 = dc.from_group(reg)
rs_Q8 = D_Q8.real_structure()
assert rs_Q8.ok and sorted(zip(rs_Q8.types, D_Q8.k_factors, D_Q8.multiplicities)) == [('H', 2, 2)] + [('R', 1, 1)] * 4
assert D_Q8.dims[0] == 8                                  # dim_R A_R = dim_C A = |Q_8| (Theorem 4.21)
```

</details>

**下一章：** [教程 5](05_random_matrices.md)：代数中的随机矩阵——库中的 Haar 酉矩阵、Weingarten 演算与能级统计。
