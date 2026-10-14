[English](../09_spin_chains.md) · [Español](../es/09_spin_chains.md) · [Français](../fr/09_spin_chains.md) · [Deutsch](../de/09_spin_chains.md) · **中文** · [日本語](../ja/09_spin_chains.md) · [Русский](../ru/09_spin_chains.md)

# 教程 9. 自旋链

本教程配合讲义的[第 9 章](../../lecture_notes/zh/09_spin_chains.md)。该章构造 XXZ 链的扇区——磁化（[命题 9.1](../../lecture_notes/zh/09_spin_chains.md#r-9-1)）、
动量（[命题 9.5](../../lecture_notes/zh/09_spin_chains.md#r-9-5)）和总自旋——并证明它们关于基态说明了什么（[定理 9.7](../../lecture_notes/zh/09_spin_chains.md#r-9-7)、[定理 9.9](../../lecture_notes/zh/09_spin_chains.md#r-9-9)）。
我们用 `SpinChain` 和 `SU2Chain` 构造链，用自由费米子和 Bethe 拟设检验本库，对稠密矩阵无法处理的大扇区使用
`krylov` 的 Lanczos 与 Krylov 方法，并计算能级统计和纠缠。

在本章中，格点为 $`0,\dots,L-1`$，$`(\mathbb C^2)^{\otimes L}`$ 的基由构型
$`|b\rangle=|b_0\dots b_{L-1}\rangle`$ 组成，自旋向上时 $`b_i=1`$（本库的约定：标号的第 $`i`$ 位
对应格点 $`i`$），$`s^a_i`$（$`a=x,y,z`$）是格点 $`i`$ 上的自旋矩阵 $`\frac12\sigma_a`$，
$`s^\pm_i=s^x_i\pm is^y_i`$，$`n_i=s^z_i+\frac12`$，$`S^a=\sum_is^a_i`$，而 $`N=S^z+\frac L2`$ 是向上自旋的
个数。链是**开链**（键 $`(i,i+1)`$，$`i<L-1`$）或**周期链**（还有键 $`(L-1,0)`$；此时
$`L\ge3`$，指标取模 $`L`$）。$`J>0`$（反铁磁），$`\hbar=1`$，$`\log`$ 为自然对数。

<a id="s-9-1"></a>

## 9.1 XXZ 链及其 $`U(1)`$ 扇区

*讲义：* [9.1 节](../../lecture_notes/zh/09_spin_chains.md#s-9-1)；[命题 9.1](../../lecture_notes/zh/09_spin_chains.md#r-9-1)（$`H`$ 保持 $`S^z`$ 守恒）、[命题 9.2](../../lecture_notes/zh/09_spin_chains.md#r-9-2)、[引理 9.3](../../lecture_notes/zh/09_spin_chains.md#r-9-3) 和[推论 9.4](../../lecture_notes/zh/09_spin_chains.md#r-9-4)（链的塔及其指标）。

`SpinChain(L).xxz` 直接由每个扇区的构型（`chain.labels[N]`）构造这些块。
我们把它们与用 Kronecker 积构造的完整 $`2^L\times2^L`$ 矩阵作比较：

```python
import functools, itertools, math, time
import numpy as np
import torch
from scipy.optimize import brentq
from torch_vn_algebra import SpinChain, SU2Chain, DensityMatrix, inclusions, krylov

torch.manual_seed(0)
f64, c128 = torch.float64, torch.complex128
close = lambda a, b, tol=1e-10: bool((torch.as_tensor(a) - torch.as_tensor(b)).abs().max() < tol)
def op(L, ops):                             # product of one-site operators {site: 2 x 2}; site 0 = lowest bit
    return functools.reduce(torch.kron, [ops.get(j, torch.eye(2, dtype=f64)) for j in reversed(range(L))])

sp = torch.tensor([[0., 0.], [1., 0.]], dtype=f64)                    # s^+ : |0> = down -> |1> = up
sm, sz = sp.T.contiguous(), torch.diag(torch.tensor([-0.5, 0.5], dtype=f64))
L = 6
chain = SpinChain(L, boundary='periodic')
h = 0.3 * torch.randn(L, dtype=f64)                                  # random fields
H_full = sum(0.5 * (op(L, {i: sp, j: sm}) + op(L, {i: sm, j: sp})) + 0.7 * op(L, {i: sz, j: sz})
             for i, j in chain.bonds) + sum(h[i] * op(L, {i: sz}) for i in range(L))   # J = 1, Delta = 0.7
Sz = sum(op(L, {i: sz}) for i in range(L))
assert close(H_full @ Sz, Sz @ H_full)                               # Proposition 9.1
H = chain.xxz(J=1.0, Delta=0.7, h=h)                                 # an element of (+)_N M_binom(6, N)
for N, k in enumerate(chain.algebra.k_factors):
    lab = chain.labels[N]                                            # configurations of sector N
    assert close(H.matrix[0, N, :k, :k], H_full[lab][:, lab])
```

链的塔 $`\mathcal A_0\subset\mathcal A_1\subset\cdots`$：`inclusions.chain_tower` 构造这些包含，
`inclusions.bratteli_diagram` 构造它们的 Bratteli 图，`Inclusion.index` 给出指标 $`\|\Lambda\|^2`$。代码检验
Pascal 法则、[引理 9.3](../../lecture_notes/zh/09_spin_chains.md#r-9-3) 以及 $`[\mathcal A_{L+1}:\mathcal A_L]=4\cos^2\frac\pi{2L+4}`$（[推论 9.4](../../lecture_notes/zh/09_spin_chains.md#r-9-4)）。

```python
tower = inclusions.chain_tower(8)                                    # A_0 c A_1 c ... c A_8
print(inclusions.bratteli_diagram(tower).levels[:6])                 # Pascal's triangle
for L_, inc in enumerate(tower):
    sizes = [math.comb(L_ + 1, N) for N in range(L_ + 2)]
    assert inc.M.k_factors == sizes and inc.is_connected
    assert list(inc.matrix @ np.array(inc.N.k_factors)) == sizes      # unitality = Pascal's rule
    n = 2 * L_ + 3
    A_path = np.eye(n, k=1) + np.eye(n, k=-1)
    assert abs(np.linalg.eigvalsh(A_path).max() - 2 * math.cos(math.pi / (n + 1))) < 1e-12   # Lemma 9.3
    assert abs(inc.index() - 4 * math.cos(math.pi / (2 * L_ + 4)) ** 2) < 1e-12           # Corollary 9.4
```

<a id="s-9-2"></a>

## 9.2 平移、动量与 Marshall 符号规则

*讲义：* [9.2 节](../../lecture_notes/zh/09_spin_chains.md#s-9-2)；[命题 9.5](../../lecture_notes/zh/09_spin_chains.md#r-9-5)（动量扇区）、[定理 9.7](../../lecture_notes/zh/09_spin_chains.md#r-9-7)（Marshall 符号规则）、[推论 9.8](../../lecture_notes/zh/09_spin_chains.md#r-9-8)（基态的动量）。

`xxz_momentum` 在归一化的基 $`|a,k\rangle`$ 中构造这些块，每条轨道取一个代表元 $`a`$。我们
独立地计数轨道，并把谱与稠密扇区作比较：

```python
L = 12
ring = SpinChain(L, boundary='periodic')
rot = lambda x: ((x << 1) | (x >> (L - 1))) & ((1 << L) - 1)        # T on labels: site i -> i + 1
for N in range(L + 1):
    orbits = {frozenset(itertools.accumulate(range(L - 1), lambda x, _: rot(x), initial=a))
              for a in ring.labels[N].tolist()}                       # {a, Ta, ..., T^{L-1} a}
    d = [sum((m * len(o)) % L == 0 for o in orbits) for m in range(L)]
    assert d == ring.momentum_dims(N)                                   # Proposition 9.5(b)
print(ring.momentum_dims(6))
Hk = ring.xxz_momentum(J=1.0, Delta=1.0, h=0.2, sector=6)               # 12 blocks instead of one of size 924
E_k = torch.sort(torch.cat([e[0] for e in Hk.eigenvalues()])).values
assert close(E_k, torch.linalg.eigvalsh(ring.xxz(1.0, 1.0, h=0.2, sector=6).matrix[0, 0]))
```

含 $`924`$ 个态的扇区 $`N=6`$ 分裂成大小为 $`80,75,78,76,78,75,80,75,78,76,78,75`$ 的块。

由 Marshall 符号规则（[定理 9.7](../../lecture_notes/zh/09_spin_chains.md#r-9-7)），对带随机场的开链和周期链，
$`(-1)^{N_B(b)}\psi(b)`$ 在每个扇区中只有一种符号；由[推论 9.8](../../lecture_notes/zh/09_spin_chains.md#r-9-8)，环的扇区 $`N`$ 的最低能级
对偶数 $`N`$ 位于动量块 $`k=0`$ 中，对奇数 $`N`$ 位于 $`k=\pi`$ 中：

```python
for bc in ('open', 'periodic'):
    c = SpinChain(10, bc, complex_valued=False)
    for N in range(11):
        w, V = torch.linalg.eigh(c.xxz(1.0, 0.7, h=0.3 * torch.randn(10, dtype=f64), sector=N).matrix[0, 0])
        x = (-1.0) ** c.bits(N)[:, 1::2].sum(-1) * V[:, 0]                  # (-1)^{N_B(b)} psi(b)
        assert (x * x.sum().sign()).min() > 0 and (len(w) == 1 or w[1] - w[0] > 1e-8)   # Theorem 9.7
for N in range(1, L):
    lowest = [e[0, 0].item() for e in ring.xxz_momentum(1.0, 0.5, sector=N).eigenvalues()]
    m0 = int(np.argmin(lowest))
    assert m0 == (0 if N % 2 == 0 else L // 2)                              # Corollary 9.8: k = 0 or pi
    assert sorted(lowest)[1] - lowest[m0] > 1e-8
```

<a id="s-9-3"></a>

## 9.3 $`\Delta=1`$ 时的 $`SU(2)`$

*讲义：* [9.3 节](../../lecture_notes/zh/09_spin_chains.md#s-9-3)；[定理 9.9](../../lecture_notes/zh/09_spin_chains.md#r-9-9)（Lieb–Mattis）。对 $`\Delta=1`$ 和 $`h=0`$，$`H\cong\bigoplus_S1_{2S+1}\otimes H_S`$，而 `SU2Chain` 在最高权向量上构造块 $`H_S`$。

我们比较十二个自旋时的块大小，用 $`U(1)`$ 扇区检验定理 9.9，计算能隙，并在一个不是链的二部图上
（$`n_A=5`$，$`n_B=3`$，$`S_*=1`$，随机耦合）检验定理 9.7 和 9.9：

```python
sym = SU2Chain(12, boundary='periodic', complex_valued=False)
H_su2 = sym.heisenberg(J=1.0)                                       # (+)_S M_{m_S}, S = 0, ..., 6
u1_sizes = [math.comb(12, N) for N in range(13)]
print(sym.multiplicities, max(u1_sizes), sum(k ** 3 for k in u1_sizes), sum(m ** 3 for m in sym.multiplicities))
E0_S = sym.lowest_energies(H_su2)[0]                                # E_0(S), S = 0, ..., 6
E0_N = [torch.linalg.eigvalsh(ring.xxz(1.0, 1.0, sector=6 + S).matrix[0, 0])[0] for S in range(7)]
assert close(E0_S, torch.stack(E0_N))                               # Theorem 9.9: sector S^z = S has spin S
assert (E0_S[1:] > E0_S[:-1] + 1e-8).all()                          # E_0(S) strictly increasing
gap = {}
for n in (6, 8, 10, 12):
    s = SU2Chain(n, 'periodic', complex_valued=False)
    e = [torch.linalg.eigvalsh(s.heisenberg(spin=S).matrix[0, 0])[0].item() for S in (0, 1)]
    gap[n] = n * (e[1] - e[0])
print({n: round(g, 4) for n, g in gap.items()})
assert 0 < gap[6] < gap[8] < gap[10] < gap[12]                      # > 0: Theorem 9.9; increase: numerical
n8, A_set = 8, (0, 1, 2, 3, 4)                                      # A = {0, ..., 4}, B = {5, 6, 7}
edges = [(0, 5), (1, 5), (1, 6), (2, 6), (3, 6), (4, 7), (0, 7), (0, 1), (2, 3), (5, 6)]   # connected
g8 = torch.Generator().manual_seed(2)
J8 = [(1 if (i in A_set) != (j in A_set) else -1) * (0.5 + torch.rand(1, generator=g8).item()) for i, j in edges]
dot = lambda i, j: 0.5 * (op(n8, {i: sp, j: sm}) + op(n8, {i: sm, j: sp})) + op(n8, {i: sz, j: sz})
H8 = sum(Jij * dot(i, j) for Jij, (i, j) in zip(J8, edges))          # J_ij > 0 between, < 0 within A, B
Sp8, Sz8 = sum(op(n8, {i: sp}) for i in range(n8)), sum(op(n8, {i: sz}) for i in range(n8))
S2 = Sp8.T @ Sp8 + Sz8 @ Sz8 + Sz8                                  # S^2 = S^- S^+ + (S^z)^2 + S^z
n_up = torch.tensor([bin(b).count('1') for b in range(2 ** n8)])
n_B = torch.tensor([bin(b >> 5).count('1') for b in range(2 ** n8)])  # up spins on B
for N in range(4, 9):                                               # sectors S^z = M = N - 4 >= 0
    idx = torch.nonzero(n_up == N).flatten()
    w, U = torch.linalg.eigh(H8[idx][:, idx])
    x, spin = (-1.0) ** n_B[idx] * U[:, 0], max(N - 4, 1)            # expected spin max(M, S_*)
    assert (x * x.sum().sign()).min() > 0 and (len(w) == 1 or w[1] - w[0] > 1e-8)          # Theorem 9.7
    assert abs(U[:, 0] @ S2[idx][:, idx] @ U[:, 0] - spin * (spin + 1)) < 1e-10           # Theorem 9.9(a)
w8 = torch.linalg.eigvalsh(H8)
assert w8[2] - w8[0] < 1e-10 and w8[3] - w8[0] > 1e-8                # ground space: one triplet, S_* = 1
```

$`SU(2)`$ 块的大小为 $`m_S=132,297,275,154,54,11,1`$，而最大的 $`U(1)`$ 块为 $`924`$，
对角化的开销 $`\sum k^3`$ 从 $`2.0\cdot10^9`$ 降到 $`5.3\cdot10^7`$。*数值观察*：$`L`$ 乘以能隙增长得很慢，
对 $`L=6,\dots,12`$ 依次为 $`4.108,4.181,4.232,4.270`$；关于它在 $`L\to\infty`$ 时的行为，我们不作任何证明。

**关联。** 各格点的自旋构成矢量算子，其球分量构成秩一张量算子（引理 B.15(a)、(b)）。由
**Wigner–Eckart 定理**（定理 B.16），它们在两个多重态之间的矩阵元等于一个 Clebsch–Gordan 系数乘以一个
与磁量子数无关的**约化矩阵元**，并且 $`\mathbf s_i\cdot\mathbf s_j=\sum_q(-1)^qs^{(i)}_qs^{(j)}_{-q}`$（引理 B.15(a)）。
`spin_correlation` 利用这一点仅从最高权向量出发在每个 $`H_S`$ 上构造 $`\mathbf s_i\cdot\mathbf s_j`$；`'direct'` 方法则改为
把键算子压缩到最高权向量上，由于 $`\mathbf s_i\cdot\mathbf s_j\in\{S_a\}'`$，由命题 4.15(b) 这是正确的。第一种方法的
实现是否正确，这里不作证明；代码检验两者彼此一致，并与按定义在 $`U(1)`$ 扇区中显式写出的
基态上计算的值一致：

```python
gs = sym.ground_state(H_su2)
S0, coef = gs.spin.item(), gs.coefficients[0, :sym.multiplicity(gs.spin.item())]
pairs = [(0, r) for r in range(12)]
corr_we = sym.expectation(sym.spin_correlation(pairs), S0, coef)
corr_dir = sym.expectation(sym.spin_correlation(pairs, method='direct'), S0, coef)
psi = sym.embed(S0, coef)                                            # the ground state in the sector N = 6
lab, z = sym.u1.labels[6], sym.u1.bits(6).to(f64) - 0.5
def s_dot_s(i, j):                                                   # <psi| s_i . s_j |psi>, from the definition
    if i == j:
        return 0.75
    flip = ((lab >> i) & 1) != ((lab >> j) & 1)
    partner = sym.u1.index_of[lab[flip] ^ ((1 << i) | (1 << j))]
    return ((psi ** 2 * z[:, i] * z[:, j]).sum() + 0.5 * (psi[flip] * psi[partner]).sum()).item()
assert S0 == 0 and close(corr_we, corr_dir)
assert close(corr_we, torch.tensor([s_dot_s(0, r) for r in range(12)], dtype=f64))
assert ((-1.0) ** torch.arange(1, 12) * corr_we[1:] > 0).all()       # signs (-1)^r
print([round(x, 4) for x in corr_we.tolist()])
```

关联 $`\langle\mathbf s_0\cdot\mathbf s_r\rangle=0.75,-0.4489,0.1878,-0.1658,0.1209,-0.1224,0.1070,\dots`$
符号交替，这正是定理 9.7 和 9.9 对每个偶数环所强制的（[讲义，9.3 节](../../lecture_notes/zh/09_spin_chains.md#s-9-3)）。

<a id="s-9-4"></a>

## 9.4 精确检验：自由费米子与 Bethe 拟设

*讲义：* [9.4 节](../../lecture_notes/zh/09_spin_chains.md#s-9-4)；由[定理 9.10](../../lecture_notes/zh/09_spin_chains.md#r-9-10)，开 XX 链扇区 $`N`$ 的能量是 $`N`$ 个互不相同的 $`\varepsilon_k=J\cos\frac{\pi k}{L+1}`$ 之和；[命题 9.12](../../lecture_notes/zh/09_spin_chains.md#r-9-12)（双磁振子，Bethe 方程）。

```python
L = 10
xx = SpinChain(L)
eps = [math.cos(math.pi * k / (L + 1)) for k in range(1, L + 1)]
for N in range(L + 1):
    E = torch.linalg.eigvalsh(xx.xxz(J=1.0, Delta=0.0, sector=N).matrix[0, 0])
    E_ff = sorted(sum(K) for K in itertools.combinations(eps, N))     # sums over N-element subsets
    assert close(E, torch.tensor(E_ff, dtype=f64))                    # all 1024 levels, Theorem 9.10
    assert abs(E[0].item() - sum(sorted(eps)[:N])) < 1e-12           # ground state: the N lowest modes
```

由[命题 9.12](../../lecture_notes/zh/09_spin_chains.md#r-9-12)，$`F(q)=(-1)^ne^{iqL}(c-\Delta e^{-iq})+c-\Delta e^{iq}`$（其中
$`K=2\pi n/L`$，$`c=\cos\frac K2`$）的每个根 $`q`$ 都给出环的扇区 $`N=2`$ 的一个特征向量 $`\Psi`$
（[讲义，9.4 节](../../lecture_notes/zh/09_spin_chains.md#s-9-4)）。代码在 $`q\in(0,\pi)`$ 中（散射态）以及在直线
$`q=i\kappa`$、$`\pi+i\kappa`$（$`\kappa>0`$）上（束缚态）寻找根，对 $`L=9`$、$`\Delta=0.6`$ 为每个根构造 $`\Psi`$，
并与含 $`36`$ 个态的扇区 $`N=2`$ 作比较：

```python
L, Delta = 9, 0.6
ring9 = SpinChain(L, boundary='periodic')
H2 = ring9.xxz(J=1.0, Delta=Delta, sector=2).matrix[0, 0]
x1, x2 = np.triu_indices(L, 1)                                       # all pairs x1 < x2
rows = ring9.index_of[torch.tensor((1 << x1) | (1 << x2))]
def roots(g, grid):                                                  # sign changes of g on the grid, refined
    v = np.array([g(t) for t in grid])
    return [brentq(g, grid[i], grid[i + 1], xtol=1e-15) for i in np.nonzero(v[:-1] * v[1:] < 0)[0]]
E_bethe, vecs, n_bound = [], [], 0
for n in range(L):
    K, c = 2 * np.pi * n / L, math.cos(np.pi * n / L)
    F = lambda q: (-1) ** n * np.exp(1j * q * L) * (c - Delta * np.exp(-1j * q)) + c - Delta * np.exp(1j * q)
    G = lambda q: (lambda z: z.real + z.imag)(np.exp(-0.5j * q * L) * F(q))    # real, or imaginary
    bound = [s + 1j * t for s in (0, np.pi)                          # q = i kappa and pi + i kappa
             for t in roots(lambda t: F(s + 1j * t).real, np.linspace(1e-3, 20, 4000))]
    n_bound += len(bound)
    for q in roots(G, (np.arange(1, 2000) + 0.5 ** 0.5) * np.pi / 2000) + bound:
        psi = np.exp(0.5j * K * (x1 + x2)) * ((c - Delta * np.exp(1j * q)) * np.exp(1j * q * (x1 - x2))
                                              - (c - Delta * np.exp(-1j * q)) * np.exp(-1j * q * (x1 - x2)))
        v = torch.zeros(36, dtype=c128)
        v[rows] = torch.tensor(psi / np.linalg.norm(psi))
        E = 9 * Delta / 4 + 2 * (c * np.cos(q) - Delta)                 # E_F + J (cos k1 + cos k2 - 2 Delta)
        assert abs(E.imag) < 1e-12 and close(H2 @ v, E.real * v, 1e-9)   # Proposition 9.12
        E_bethe.append(E.real)
        vecs.append(v)
assert torch.linalg.matrix_rank(torch.stack(vecs)) == 36             # a complete eigenbasis
assert close(torch.sort(torch.tensor(E_bethe)).values, torch.linalg.eigvalsh(H2))
print(len(E_bethe), n_bound)
ring8, Phi = SpinChain(8, 'periodic'), torch.zeros(28, dtype=c128)     # even L: Phi = sum_x (-1)^x |x, x+1>
for x in range(8):
    Phi[ring8.index_of[(1 << x) | (1 << (x + 1) % 8)]] = (-1) ** x
assert close(ring8.xxz(1.0, Delta, sector=2).matrix[0, 0] @ Phi, (8 * Delta / 4 - Delta) * Phi)   # E_F - J Delta
```

*数值观察*：对 $`L=9`$ 和 $`\Delta=0.6`$，搜索找到 $`36`$ 个根，即 $`34`$ 个散射态和 $`2`$ 个束缚态，它们的向量
是线性无关的特征向量（在舍入误差范围内）：全部 $`36`$ 个本征态都来自 Bethe 方程。对其他 $`L`$，完备性没有得到证明，
而对偶数 $`L`$ 和 $`\Delta\ne0`$ 它不成立（[讲义，9.4 节](../../lecture_notes/zh/09_spin_chains.md#s-9-4)）：向量 $`\Phi=\sum_x(-1)^x|x,x+1\rangle`$ 与所有 Bethe 向量正交，
代码的最后几行对 $`L=8`$ 检验了 $`H\Phi=(E_F-J\Delta)\Phi`$。

<a id="s-9-5"></a>

## 9.5 Lanczos 与 Krylov 方法

*讲义：* [9.5 节](../../lecture_notes/zh/09_spin_chains.md#s-9-5)；[命题 9.13](../../lecture_notes/zh/09_spin_chains.md#r-9-13)（Lanczos：$`T_m=Q_m^*HQ_m`$ 是三对角的）、[定理 9.14](../../lecture_notes/zh/09_spin_chains.md#r-9-14)（Ritz 值，即 $`T_m`$ 的特征值，满足 $`\theta^{(m)}_j\ge\lambda_j`$ 并随 $`m`$ 递减）、[命题 9.15](../../lecture_notes/zh/09_spin_chains.md#r-9-15)（Krylov 时间演化）。

扇区 $`N=L/2`$ 的维数在 $`L=16`$ 时为 $`12870`$，在 $`L=20`$ 时为 $`184756`$，但 $`H`$ 是稀疏的（每行 $`O(L)`$ 个非零
元素，`xxz_sparse`），而 Krylov 方法只需要乘积 $`Hv`$。设 $`H`$ 是 $`n`$ 阶厄米矩阵，
$`v\ne0`$。**Krylov 空间**为 $`\mathcal K_m=\mathrm{span}\{v,Hv,\dots,H^{m-1}v\}`$。

定理 9.14 并没有说明 $`\theta^{(m)}_1`$ 趋近 $`\lambda_1`$ 的快慢；我们只在下面观察它。在浮点运算中
$`q_j`$ 会失去正交性；`krylov.lanczos` 进行完全再正交化，代码检验计算出的 $`q_j`$ 在舍入误差范围内
规范正交且 $`Q_m^*HQ_m=T_m`$。对开放 XX 链，扇区 $`N`$ 的精确基态能量是 $`N`$ 个最低的 $`\varepsilon_k`$ 之和
（定理 9.10(c)）：

```python
L, N = 12, 6
xx = SpinChain(L, complex_valued=False)
Hs = xx.xxz_sparse(J=1.0, Delta=0.0, sector=N)                      # sparse: diagonal + hopping
alpha, beta, Q = krylov.lanczos(Hs.matvec, torch.randn(1, xx.sector_dim(N), dtype=f64), 40)
T = torch.diag(alpha[0]) + torch.diag(beta[0], 1) + torch.diag(beta[0], -1)
assert close(Q[0] @ Q[0].T, torch.eye(40, dtype=f64)) and close(Q[0] @ Hs.to_dense()[0] @ Q[0].T, T)
ritz = torch.tensor([torch.linalg.eigvalsh(T[:m, :m])[0].item() for m in range(1, 41)])   # T_m = leading block
E_xx = sum(sorted(math.cos(math.pi * k / (L + 1)) for k in range(1, L + 1))[:N])
assert (ritz[1:] <= ritz[:-1] + 1e-12).all() and (ritz >= E_xx - 1e-12).all()             # Theorem 9.14
print([f"{(ritz[m - 1] - E_xx).item():.1e}" for m in (5, 10, 20, 40)])
```

*数值观察*：最低 Ritz 值的误差单调下降（由定理 9.14，这是必然的），在 $`5`$、$`10`$ 和 $`20`$ 步之后分别为
$`0.42`$、$`3.1\cdot10^{-3}`$ 和 $`3.1\cdot10^{-6}`$，在 $`40`$ 步之后达到舍入误差水平，而维数为 $`924`$。
`krylov.ground_state` 反复重启 Lanczos（除非给定 `v0`，否则从一个高斯随机向量出发），直到其归一化
Ritz 对 $`(E,\psi)`$ 的残差 $`\|H\psi-E\psi\|`$ 低于 $`10^{-10}\max(1,|E|)`$（至多重启 $`20`$ 次，否则发出警告），
并返回该残差。小的残差保证附近有一个特征值，但不保证它是最低的（[讲义，9.5 节](../../lecture_notes/zh/09_spin_chains.md#s-9-5)），
所以我们与稠密对角化作比较，这是我们唯一的参照（在这里计算，
所以这些值是数值结果）：直到 $`L=12`$ 用扇区 $`N=L/2`$，对 $`L=16`$ 用维数为 $`810`$ 的动量块 $`k=0`$，
由推论 9.8 它包含扇区 $`N=8`$ 的基态。对 $`L=4`$，值 $`-2J`$ 是精确的（[讲义，9.5 节](../../lecture_notes/zh/09_spin_chains.md#s-9-5)）。

```python
E_dense = {}                                                         # numerical references, dense
for n in (4, 6, 8, 10, 12):
    E_dense[n] = torch.linalg.eigvalsh(SpinChain(n, 'periodic').heisenberg(sector=n // 2).matrix[0, 0])[0].item()
    print(n, round(E_dense[n] / n, 6))
assert abs(E_dense[4] + 2) < 1e-12                                   # exact value for L = 4
ring16 = SpinChain(16, 'periodic')
k0 = ring16.momentum_dims(8)[0]                                      # dimension of the block k = 0
H16 = ring16.xxz_momentum(J=1.0, Delta=1.0, sector=8).matrix[0, 0, :k0, :k0]
E_dense[16] = torch.linalg.eigvalsh(H16)[0].item()                   # ground state of N = 8 (Corollary 9.8)
t0 = time.time()
E16, psi16, res = krylov.ground_state(SpinChain(16, 'periodic', complex_valued=False).xxz_sparse(sector=8),
                                      generator=torch.Generator().manual_seed(0))
print(f"L = 16: E0 = {E16.item():.10f}, residual {res.item():.1e}, {time.time() - t0:.2f} s; k = 0 block: {k0}")
assert abs(E16.item() - E_dense[16]) < 1e-9 and res.item() < 1e-9
```

对 $`L=16`$（维数 $`12870`$），Lanczos 在单个 CPU 线程上用时远少于一秒。*数值观察*：每个格点的能量
$`-0.5,-0.467129,-0.456387,-0.451545,-0.448949`$（$`L=4,\dots,12`$）随 $`L`$ 增大；关于其极限我们不作任何证明。

**时间演化。** 对 $`e^{-i\tau H}v`$，我们投影到 $`\mathcal K_m`$ 上：$`w_\tau=\|v\|\,Q_me^{-i\tau T_m}e_1`$。

由[命题 9.15](../../lecture_notes/zh/09_spin_chains.md#r-9-15)，对每个 $`\tau`$ 有 $`\|w_\tau\|=\|v\|`$ 和 $`\langle w_\tau,Hw_\tau\rangle=\langle v,Hv\rangle`$，并且若
$`\mathcal K_m`$ 是 $`H`$-不变的，则 $`w_\tau=e^{-i\tau H}v`$。

所以每个 Krylov 步都保持范数和能量，不论精确与否；除非 $`\mathcal K_m`$ 是不变的，命题 9.15 对精度什么也
没有说。`krylov.evolve` 把时间区间细分为短步，每一步都在 $`m=30`$ 的 Krylov 空间中完成，并在一个启发式
误差指标较大时发出警告；我们不为它证明任何误差界，而是把它与由稠密特征分解（定理 1.5）计算出的
$`e^{-itH}`$ 作比较。我们在 Heisenberg 环上对 Néel 态 $`|{\uparrow\downarrow\uparrow\downarrow}\cdots\rangle`$ 作淬火，跟踪交错磁化
$`m_s=\frac2L\sum_i(-1)^i\langle s^z_i\rangle`$，并对 $`L=12`$ 与稠密传播子作比较：

```python
times, m_s = [0.0, 0.5, 1.0, 2.0, 4.0], {}
stag = lambda c, N, out: [((p[0].abs() ** 2) @ (c.bits(N).to(f64) - 0.5) @ torch.tensor(
    [(-1.0) ** i for i in range(c.L)], dtype=f64) * 2 / c.L).item() for p in out]
for n in (12, 16):
    c = SpinChain(n, 'periodic')
    Hn = c.xxz_sparse(J=1.0, Delta=1.0, sector=n // 2)
    psi0 = c.vector_in_sector('10' * (n // 2))[None]
    out = krylov.evolve(Hn, psi0, times)                               # (times, batch, dim)
    energy = torch.stack([(p.conj() * Hn.matvec(p)).sum(-1).real for p in out])
    assert close(out.norm(dim=-1), 1.0) and close(energy, energy[0])   # Proposition 9.15
    m_s[n] = stag(c, n // 2, out)
    print(n, [round(x, 4) for x in m_s[n]])
    if n == 12:
        w, U = torch.linalg.eigh(c.xxz(1.0, 1.0, sector=6).matrix[0, 0])
        for t, p in zip(times, out):
            assert close(p[0], U @ (torch.exp(-1j * t * w) * (U.mH @ psi0[0])))
print(f"{max(abs(a - b) for a, b in zip(m_s[12][:4], m_s[16][:4])):.1e}")   # L = 12 vs 16 up to t = 2
```

*数值观察*：$`m_s`$ 从 $`1`$ 经过 $`0.770`$ 和 $`0.279`$ 衰减到 $`t=2`$ 时的 $`-0.183`$，两种尺寸的结果相差在
$`1.7\cdot10^{-5}`$ 以内；在 $`t=4`$ 时它们不同（$`0.026`$ 对 $`0.035`$）。

<a id="s-9-6"></a>

## 9.6 无序：能级统计与纠缠

*讲义：* [9.6 节](../../lecture_notes/zh/09_spin_chains.md#s-9-6)；[定理 9.16](../../lecture_notes/zh/09_spin_chains.md#r-9-16)（Peschel 公式）。

**随机场 Heisenberg 链**的外场 $`h_i`$ 相互独立，在 $`[-W,W]`$ 上均匀分布。它保持 $`S^z`$ 守恒
（命题 9.1），所以能级统计按扇区 $`N`$ 分别计算（定理 5.15 和 5.4 节：混合扇区会伪造出独立的能级）。
若两个外场不同，则 $`THT^*-H=\sum_i(h_{i-1}-h_i)s^z_i\ne0`$，所以平移不再是对称性。我们在扇区 $`N=L/2`$
的中间一半能级中计算 $`\langle r\rangle`$，给出关于独立实现的一个标准误差（每个实现给出一个平均比值），
并与两个参照值比较：$`2\log2-1\approx0.3863`$，即独立指数间距的精确平均比值（命题 5.14(a)），以及在这里
抽样得到的、与 $`L=10`$ 时该扇区大小相同（$`252`$）的 GOE 矩阵（按 5.3 节的定义；代码抽样 $`G+G^T`$，它是
GOE 矩阵的两倍，这不改变任何比值）的平均比值。这些比较是数值观察，而不是定理：关于随机链的谱，
我们不作任何证明。

```python
gen = torch.Generator().manual_seed(1)
r_mean = {}
for n, R in ((8, 200), (10, 100)):
    c = SpinChain(n, 'periodic', complex_valued=False)
    k = c.sector_dim(n // 2)
    for W in (0.5, 1, 2, 4, 8, 16):
        E = torch.linalg.eigvalsh(c.random_field_heisenberg(W, R, sector=n // 2, generator=gen).matrix[:, 0])
        r = SpinChain.level_spacing_ratio(E[:, k // 4: 3 * k // 4])     # one value per realisation
        r_mean[n, W] = (r.mean().item(), (r.std() / math.sqrt(R)).item())
    print(n, '  '.join(f"W={W}: {m:.3f}+-{s:.3f}" for (n_, W), (m, s) in r_mean.items() if n_ == n))
X = torch.randn(100, k, k, dtype=f64, generator=gen)                  # 100 GOE matrices of size k = 252
r_goe = SpinChain.level_spacing_ratio(torch.linalg.eigvalsh(X + X.mT)[:, k // 4: 3 * k // 4])
print(f"GOE({k}): {r_goe.mean().item():.3f}+-{(r_goe.std() / 10).item():.3f}, Poisson: {2 * math.log(2) - 1:.4f}")
assert abs(r_mean[10, 0.5][0] - r_goe.mean().item()) < 0.015               # numerical observation
assert abs(r_mean[10, 16][0] - (2 * math.log(2) - 1)) < 0.015               # numerical observation
assert r_mean[10, 1][0] > r_mean[10, 2][0] > r_mean[10, 4][0] > r_mean[10, 8][0]
```

*数值观察*：对 $`L=10`$，平均比值从 $`0.530`$（$`W=0.5`$，接近抽样得到的 GOE 值 $`0.533`$）经过 $`0.484`$（$`W=2`$）
和 $`0.410`$（$`W=4`$）变到 $`0.384`$（$`W=16`$，接近 $`2\log2-1`$），标准误差为 $`0.003`$–$`0.004`$。对 $`L=8`$，$`W=0.5`$
处的值（$`0.507`$）低于 $`W=1`$ 处的值（$`0.525`$）。在 $`W=0`$ 时，环与 $`T`$（9.2 节）以及总自旋（9.3 节）可交换，
所以它的扇区 $`N=L/2`$ 分裂成若干块；外场只通过 $`\sum_ih_is^z_i`$ 耦合这些块，其算子范数至多为 $`\frac{WL}2`$。

> **物理意义**（物理解释，而非本课程的定理）。物理学家把 $`\langle r\rangle`$ 的变化解读为从像 GOE 的能级那样
> 相互排斥的能级（弱无序）到表现得像相互独立的能级的交叉过渡，并将后者归因于本征态的
> **多体局域化**（强无序）。这一交叉过渡随 $`L`$ 漂移；它在 $`L\to\infty`$ 时位于何处，以及局域相究竟是否
> 存在，都是物理文献中的开放问题（D. Luitz, N. Laflorencie, F. Alet, *Phys. Rev. B* 91 (2015) 081103；
> J. Šuntajs, J. Bonča, T. Prosen, L. Vidmar, *Phys. Rev. E* 102 (2020) 062144）。$`L=8`$、$`W=0.5`$ 处的下凹被解读为
> $`W=0`$ 时额外块结构的残留（参见 5.4 节）。

**自由费米子的纠缠。** 对开放 XX 链的本征态，前 $`n`$ 个格点的纠缠熵
由一个 $`n\times n`$ 矩阵决定。

`entanglement_entropy` 逐个扇区计算约化态（它与块的磁化可交换，
命题 7.8）及其熵。扇区 $`N`$ 的基态是由 $`N`$ 个最低模构成的 $`\Phi`$，
由于各 $`\varepsilon_k`$ 互不相同，它是单的：

```python
L, N = 12, 6
xx = SpinChain(L)
w, V = torch.linalg.eigh(xx.xxz(J=1.0, Delta=0.0, sector=N).matrix[0, 0])
rho = DensityMatrix(xx.sector_algebra(N), matrix=torch.outer(V[:, 0], V[:, 0].conj())[None, None])
eps_k, u = np.linalg.eigh(np.eye(L, k=1) / 2 + np.eye(L, k=-1) / 2)  # h = A_L / 2, ascending
C = u[:, :N] @ u[:, :N].T                                            # correlation matrix of the N lowest modes
bern = lambda nu: -(nu * np.log(nu) + (1 - nu) * np.log(1 - nu))
S_peschel = [bern(np.clip(np.linalg.eigvalsh(C[:n, :n]), 1e-15, 1 - 1e-15)).sum() for n in range(1, L)]
S_lib = [xx.entanglement_entropy(rho, n, sector=N).item() for n in range(1, L)]
assert np.allclose(S_lib, S_peschel, atol=1e-9)                       # Theorem 9.16
print([round(s, 4) for s in S_lib[:6]])
```

*数值观察*：熵 $`0.6931,0.5199,0.7456,0.6283,0.7719,0.6560`$（$`n=1,\dots,6`$）随 $`n`$ 的奇偶性振荡。
定理 9.16 归功于 I. Peschel（*J. Phys. A* 36 (2003) L205）。

<a id="s-9-7"></a>

## 9.7 可与之比较的工具

QuSpin（Python）在磁化、动量、宇称和自旋反演的扇区中构造自旋链哈密顿量，
并使用 SciPy 的稀疏本征求解器和矩阵指数；专门的精确对角化程序
（XDiag、SPINPACK、HΦ）针对大型稀疏问题作了优化，借助格点对称性和大型计算机，
这类程序已达到约 40–50 个自旋；ITensor、TeNPy 和 TensorKit.jl（带有 $`SU(2)`$ 融合树）用
矩阵乘积态处理长链。本库所增加的东西不多而具体：把链作为
`TypeIAlgebra` $`\bigoplus_NM_{\binom LN}`$，使第 3–8 章的态、信道、熵和动力学
可以逐扇区应用；带有指标的包含塔；由最高权向量得到的 $`SU(2)`$ 块以及
Wigner–Eckart 关联；一次调用处理一批无序实现；以及从不离开扇区的约化态。
它没有反射或自旋翻转扇区，动量扇区仅适用于均匀 XXZ 链，稠密
方法到 $`L\approx16`$，Krylov 方法到 $`L\approx24`$（在大型 GPU 上约 $`26`$），并且没有张量网络。

<a id="x-9"></a>

## 习题

讲义的[习题 1–3](../../lecture_notes/zh/09_spin_chains.md#x-9)：$`L=4`$、$`N=2`$ 时动量扇区的维数 $`d_{2,m}=2,1,2,1`$；
Heisenberg 环基态的求和规则 $`E_0=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle`$ 和 $`\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle=0`$；
以及从扇区 $`N=L/2`$ 中所有系数都为 $`1`$ 的向量出发的 `krylov.ground_state`，它返回谱的顶端 $`E=JL/4`$。
用数值方法检验它们。

<details>
<summary>解答</summary>

```python
assert SpinChain(4, 'periodic').momentum_dims(2) == [2, 1, 2, 1]
```

```python
assert abs(12 * corr_we[1].item() - gs.energy.item()) < 1e-10 and abs(corr_we.sum().item()) < 1e-10
```

```python
c10 = SpinChain(10, 'periodic', complex_valued=False)
E_sym, _, _ = krylov.ground_state(c10.xxz_sparse(1.0, 1.0, sector=5), v0=torch.ones(1, 252, dtype=f64))
assert abs(E_sym.item() - 2.5) < 1e-10 and E_sym.item() > E_dense[10] + 7  # far from E_0
```

</details>

**下一章：** [教程 10. 量子计算](10_quantum_computing.md)。
