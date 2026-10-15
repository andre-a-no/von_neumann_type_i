[English](../10_quantum_computing.md) · [Español](../es/10_quantum_computing.md) · [Français](../fr/10_quantum_computing.md) · [Deutsch](../de/10_quantum_computing.md) · **中文** · [日本語](../ja/10_quantum_computing.md) · [Русский](../ru/10_quantum_computing.md)

# 教程 10. 量子计算

本教程配合讲义的[第 10 章](../../lecture_notes/zh/10_quantum_computing.md)，该章用量子计算机上噪声的 Kraus 算子所生成的 $`*`$-代数来描述噪声。
我们计算不动点代数与无噪声子系统（[定理 10.2](../../lecture_notes/zh/10_quantum_computing.md#r-10-2)），在重复码上检验 Knill–Laflamme 定理（[定理 10.5](../../lecture_notes/zh/10_quantum_computing.md#r-10-5)）
以及扇区标记的纠错（[命题 10.6](../../lecture_notes/zh/10_quantum_computing.md#r-10-6)），并在汉明重量扇区的代数 $`\bigoplus_wM_{\binom nw}`$ 中模拟带振幅阻尼与
后选择的 XY 混合算子拟设（[命题 10.8](../../lecture_notes/zh/10_quantum_computing.md#r-10-8)–[10.10](../../lecture_notes/zh/10_quantum_computing.md#r-10-10)）。

<a id="s-10-1"></a>

## 10.1 噪声代数、无噪声子系统与不动点

*讲义：* [10.1 节](../../lecture_notes/zh/10_quantum_computing.md#s-10-1)；[定理 10.2](../../lecture_notes/zh/10_quantum_computing.md#r-10-2)、[推论 10.3](../../lecture_notes/zh/10_quantum_computing.md#r-10-3)、[命题 10.4](../../lecture_notes/zh/10_quantum_computing.md#r-10-4)。

由[讲义，推论 10.3](../../lecture_notes/zh/10_quantum_computing.md#r-10-3)，$`n`$ 个量子比特上的集体退相位 $`\Phi(\rho)=\tfrac12(U\rho U^*+U^*\rho U)`$，
$`U=e^{-i\theta Z_{\mathrm{tot}}}`$，其相互作用代数为 $`\mathcal A=\mathrm{span}\{P_0,\dots,P_n\}`$
（$`P_w`$ 是到由汉明重量为 $`w`$ 的基向量张成的子空间 $`H_w`$ 上的投影），并且
$`\mathrm{Fix}(\Phi)=\mathcal A'=\bigoplus_wB(H_w)`$（[定理 10.2](../../lecture_notes/zh/10_quantum_computing.md#r-10-2)）。

交换子是 $`T\mapsto(AT-TA)_A`$ 的零空间（第 2 章），不动点则是
$`\sum_iE_i\otimes\bar E_i-1`$（即 $`\Phi-\mathrm{id}`$ 在按行优先向量化的矩阵上的矩阵）的零空间。对 $`n=4`$，
推论中的三个空间的维数分别为 $`5`$、$`70=1+16+36+16+1`$ 和 $`70`$：

```python
import itertools
import math
import torch

torch.manual_seed(0)
dt = torch.complex128
I2, X = torch.eye(2, dtype=dt), torch.tensor([[0, 1], [1, 0]], dtype=dt)
Y, Z = torch.tensor([[0, -1j], [1j, 0]], dtype=dt), torch.tensor([[1, 0], [0, -1]], dtype=dt)

def on(op, j, n):                             # op on qubit j of n (j = 0 is the leftmost factor)
    out = torch.ones(1, 1, dtype=dt)
    for i in range(n):
        out = torch.kron(out, op if i == j else I2)
    return out

def null_space(L, tol=1e-9):                  # orthonormal columns spanning {x : L x = 0}
    w, V = torch.linalg.eigh(L.mH @ L)
    return V[:, w < tol]

def commutant_basis(ops):
    d, ops = ops[0].shape[0], [A.resolve_conj().contiguous() for A in ops]
    eye = torch.eye(d, dtype=dt)
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    return [v.reshape(d, d) for v in null_space(L).T]

def fixed_points(kraus):                      # vec(E X E^*) = (E (x) conj E) vec(X)
    d = kraus[0].shape[0]
    S = sum(torch.kron(E.resolve_conj(), E.conj().resolve_conj()) for E in kraus)
    return [v.reshape(d, d) for v in null_space(S - torch.eye(d * d, dtype=dt)).T]

n, theta = 4, 0.3
U = torch.linalg.matrix_exp(-1j * theta * sum(on(Z, j, n) for j in range(n)))
kraus = [U / math.sqrt(2), U.mH / math.sqrt(2)]                  # unital: U U^* = U^* U = 1
comm = commutant_basis(kraus + [E.mH for E in kraus])             # A'
A_int, fix = commutant_basis(comm), fixed_points(kraus)           # A = A'' (Theorem 2.2), Fix(Phi)
print(len(A_int), len(comm), len(fix))
assert (len(A_int), len(comm), len(fix)) == (5, 70, 70)
span = torch.linalg.svdvals(torch.stack([F.reshape(-1) for F in comm + fix]))
assert int((span > 1e-9).sum()) == 70                             # Fix(Phi) = A' (Theorem 10.2)
```

$`\bigoplus_wM_{\binom nw}`$ 是 `SpinChain`（第 9 章）的代数，其扇区由向上自旋的个数，即汉明重量
标记。在那里，Kraus 算子是中心元，
$`U=\sum_we^{-i\theta(n-2w)}P_w`$，而信道是恒等映射：本库存储的正是无噪声的部分。

```python
from torch_vn_algebra import SpinChain, TypeIAlgebra, DensityMatrix, channels, states, dynamics

M = SpinChain(n).algebra                                          # (+)_w M_binom(4, w)
print(M.k_factors)
assert sum(k * k for k in M.k_factors) == len(fix)
phases = [torch.exp(torch.tensor(-1j * theta * (n - 2 * w), dtype=dt)) for w in range(n + 1)]
Phi_M = channels.kraus_channel([M.central([z / math.sqrt(2) for z in phases]),
                                M.central([z.conj() / math.sqrt(2) for z in phases])])
rho = states.random_density_matrix(M, batch_size=3)
assert (Phi_M(rho).matrix - rho.matrix).abs().max() < 1e-12      # no state of M is disturbed
```

集体转动 $`e^{-i\theta\,\mathbf n\cdot\mathbf S}`$（$`S_\alpha=\tfrac12\sum_j\sigma_{\alpha,j}`$）属于由 $`1,S_x,S_y,S_z`$ 生成的 $`*`$-代数
$`\mathcal B`$。对 $`n=4`$，讲义的[10.1 节](../../lecture_notes/zh/10_quantum_computing.md#s-10-1)表明 $`\dim\mathcal B=35`$、$`\dim\mathcal B'=14`$，且对
$`K=\ker S_x\cap\ker S_y\cap\ker S_z`$ 有 $`\dim K=2`$；代码求得相同的数值，并检验[命题 10.4](../../lecture_notes/zh/10_quantum_computing.md#r-10-4)：以随机集体转动为
Kraus 算子的信道保持支撑在 $`K`$ 中的态不变，而会改变一般的态。

```python
S = [sum(on(P, j, n) for j in range(n)) / 2 for P in (X, Y, Z)]
comm_S = commutant_basis(S)
print(len(comm_S), len(commutant_basis(comm_S)))
assert len(comm_S) == 14 and len(commutant_basis(comm_S)) == 35
K = null_space(torch.cat(S))                                      # columns: a basis of K
assert K.shape[1] == 2
axes = torch.nn.functional.normalize(torch.randn(3, 3, dtype=torch.float64), dim=1)
rotations = [torch.linalg.matrix_exp(-1j * t * sum(c * Sa for c, Sa in zip(ax, S))) / math.sqrt(3)
             for t, ax in zip(torch.rand(3, dtype=torch.float64) * math.pi, axes)]
apply = lambda kr, r: sum(E @ r @ E.mH for E in kr)
for G in (K @ torch.randn(2, 2, dtype=dt), torch.randn(16, 16, dtype=dt)):  # supported in K / generic
    r = G @ G.mH / (G @ G.mH).trace()
    assert ((apply(rotations, r) - r).abs().max() < 1e-12) == (G.shape[1] == 2)
```

<a id="s-10-2"></a>

## 10.2 纠错：Knill–Laflamme 定理

*讲义：* [10.2 节](../../lecture_notes/zh/10_quantum_computing.md#s-10-2)；[定理 10.5](../../lecture_notes/zh/10_quantum_computing.md#r-10-5)。

讲义[10.2 节](../../lecture_notes/zh/10_quantum_computing.md#s-10-2)中的比特翻转例子：对 $`\mathcal C=\mathrm{span}\{|000\rangle,|111\rangle\}`$ 和"至多一个比特翻转"，
[定理 10.5](../../lecture_notes/zh/10_quantum_computing.md#r-10-5) 的 Knill–Laflamme 条件成立，且 $`\Lambda=\mathrm{diag}(p_0,\dots,p_3)`$，而恢复
$`\mathcal R(\sigma)=\sum_{k=0}^3PX_k\sigma X_kP`$ 能还原码中的每个态。若每个量子比特以概率 $`p`$ 独立地翻转，则该条件不成立，
恢复后 $`|0_L\rangle`$ 的保真度为 $`1-3p^2+2p^3=0.972`$（$`p=0.1`$）。辅助函数 `kl_matrix` 返回 $`\Lambda`$，若条件不成立则返回
`None`。在本库中，$`M_8`$ 是只有一个块的代数，而 `R @ Phi` 就是复合 $`\mathcal R\circ\Phi`$：

```python
M8 = TypeIAlgebra([8], [8], precision='double')
op8 = lambda T: M8.operator(T.reshape(1, 1, 8, 8))
state8 = lambda v: DensityMatrix(M8, matrix=torch.outer(v, v.conj()).reshape(1, 1, 8, 8))
ket = lambda bits: torch.eye(8, dtype=dt)[int(bits, 2)]
def kl_matrix(kraus, P, tol=1e-12):
    """(lambda_ij) with P E_i^* E_j P = lambda_ij P, or None if the Knill-Laflamme condition fails."""
    lam = torch.tensor([[(P @ a.mH @ b @ P).trace() / P.trace() for b in kraus] for a in kraus])
    ok = all((P @ a.mH @ b @ P - lam[i, j] * P).abs().max() < tol
             for i, a in enumerate(kraus) for j, b in enumerate(kraus))
    return lam if ok else None
def independent(op, prob):                    # op on each qubit independently with probability prob
    return [math.sqrt(prob ** sum(s) * (1 - prob) ** (3 - sum(s)))
            * torch.linalg.multi_dot([on(op, j, 3) if s[j] else torch.eye(8, dtype=dt) for j in range(3)])
            for s in itertools.product([0, 1], repeat=3)]

P = torch.outer(ket('000'), ket('000')) + torch.outer(ket('111'), ket('111'))
flips = [torch.eye(8, dtype=dt)] + [on(X, j, 3) for j in range(3)]
E_one = [math.sqrt(pk) * Xk for pk, Xk in zip((0.7, 0.1, 0.1, 0.1), flips)]
assert torch.allclose(kl_matrix(E_one, P), torch.diag(torch.tensor([0.7, 0.1, 0.1, 0.1], dtype=dt)))
Phi_one = channels.kraus_channel([op8(E) for E in E_one])
R = channels.kraus_channel([op8(P @ Xk) for Xk in flips])        # syndrome measurement + correction
assert R.is_trace_preserving()
a, b = torch.nn.functional.normalize(torch.randn(2, dtype=dt), dim=0)
rho_L = state8(a * ket('000') + b * ket('111'))
assert abs((R @ Phi_one)(rho_L).fidelity(rho_L).item() - 1) < 1e-12       # Theorem 10.5
p = 0.1                                                           # independent flips: KL fails
E_ind = independent(X, p)
assert kl_matrix(E_ind, P) is None
zero_L = state8(ket('000'))
F = (R @ channels.kraus_channel([op8(E) for E in E_ind]))(zero_L).fidelity(zero_L).item()
print(f"fidelity after recovery {F:.4f}, exact {1 - 3 * p**2 + 2 * p**3:.4f}")
assert abs(F - (1 - 3 * p**2 + 2 * p**3)) < 1e-12
```

<a id="s-10-3"></a>

## 10.3 算子代数纠错

*讲义：* [10.3 节](../../lecture_notes/zh/10_quantum_computing.md#s-10-3)；[命题 10.6](../../lecture_notes/zh/10_quantum_computing.md#r-10-6)。

讲义[10.3 节](../../lecture_notes/zh/10_quantum_computing.md#s-10-3)中的例子：重复码先受到独立的相位翻转（每个量子比特的概率为 $`q`$），再受到至多一个比特翻转。
量子比特不受保护（Knill–Laflamme 条件不成立），但经典比特受到保护：[命题 10.6](../../lecture_notes/zh/10_quantum_computing.md#r-10-6) 的条件 2 对
$`P_0=|000\rangle\langle000|`$、$`P_1=|111\rangle\langle111|`$ 成立，多数表决 $`N_0=\sum_{|b|\le1}|b\rangle\langle b|`$，
$`N_1=1-N_0`$ 读出该比特。经过 10.2 节的恢复后，$`|{+_L}\rangle`$ 的保真度为
$`(1+(1-2q)^3)/2=0.8645`$（$`q=0.05`$ 时）。

```python
q = 0.05
Phi2 = Phi_one @ channels.kraus_channel([op8(E) for E in independent(Z, q)])   # dephasing, then a flip
E2 = [Phi2.kraus[0, i, 0] for i in range(Phi2.kraus_rank)]        # its 32 Kraus operators
P0, P1 = torch.outer(ket('000'), ket('000')), torch.outer(ket('111'), ket('111'))
assert kl_matrix(E2, P) is None                                   # the qubit is not protected
assert max((P0 @ e.mH @ f @ P1).abs().max().item() for e in E2 for f in E2) == 0      # condition 2
N0 = torch.diag(torch.tensor([float(bin(i).count('1') <= 1) for i in range(8)], dtype=dt))
for N, Pd in ((N0, P0), (torch.eye(8, dtype=dt) - N0, P1)):       # P Phi^*(N_d) P = P_d
    assert (P @ sum(e.mH @ N @ e for e in E2) @ P - Pd).abs().max() < 1e-12
votes = states.born_probabilities(Phi2(rho_L), [op8(N0), op8(torch.eye(8, dtype=dt) - N0)])[0]
assert torch.allclose(votes, torch.stack([abs(a) ** 2, abs(b) ** 2]))      # the classical bit survives
plus_L = state8((ket('000') + ket('111')) / math.sqrt(2))
F_plus = (R @ Phi2)(plus_L).fidelity(plus_L).item()
print(f"fidelity of |+_L> {F_plus:.4f}, exact {(1 + (1 - 2 * q)**3) / 2:.4f}")
assert abs(F_plus - (1 + (1 - 2 * q) ** 3) / 2) < 1e-12
```

<a id="s-10-4"></a>

## 10.4 带约束的混合算子

*讲义：* [10.4 节](../../lecture_notes/zh/10_quantum_computing.md#s-10-4)；[命题 10.8](../../lecture_notes/zh/10_quantum_computing.md#r-10-8)、[命题 10.9](../../lecture_notes/zh/10_quantum_computing.md#r-10-9)、[命题 10.10](../../lecture_notes/zh/10_quantum_computing.md#r-10-10)。

量子交替算子拟设交替作用 $`e^{-i\gamma H_C}`$（$`H_C`$ 是对角的代价哈密顿量）与 $`e^{-i\beta H_M}`$
（$`H_M=\tfrac14\sum_{(i,j)}(X_iX_j+Y_iY_j)`$ 是环上的 XY 混合算子）。由[讲义，命题 10.8](../../lecture_notes/zh/10_quantum_computing.md#r-10-8)，它们都属于
$`\mathcal M=\{P_0,\dots,P_n\}'\cong\bigoplus_wM_{\binom nw}`$，因此拟设永远不会离开其初态所在的汉明重量扇区。

在本库中，$`\mathcal M`$ 就是 `SpinChain(n).algebra`（比特 $`1`$ = 自旋向上，所以 $`S^z_i=-\tfrac12Z_i`$），而 $`H_M`$
是它的 XXZ 哈密顿量 $`J\sum_{(i,j)}\big(\tfrac12(S^+_iS^-_j+S^-_iS^+_j)+\Delta S^z_iS^z_j\big)`$，取 $`J=1`$，
$`\Delta=0`$，因为 $`S^+_iS^-_j+S^-_iS^+_j=\tfrac12(X_iX_j+Y_iY_j)`$；我们将它与稠密的
$`64\times64`$ 矩阵进行对照。对于 $`H_C`$，我们取由边 $`(2,3)`$、$`(5,0)`$ 连接的两个三角形 $`\{0,1,2\}`$、$`\{3,4,5\}`$
上的图二分问题：$`H_C=\sum_{(i,j)\in E}\tfrac12(1-Z_iZ_j)`$（$`Z_iZ_j=4S^z_iS^z_j`$）计算被切割的
边数。最优的二分切割 $`2`$ 条边（即分开两个三角形；任何其他由三个顶点构成的集合都会把两个
三角形都分割开，而一个被分割的三角形有两条被切割的边）；我们从切割 $`6`$ 条边的 $`|101010\rangle`$ 出发。最后，整个
拟设用稠密的 $`64\times64`$ 矩阵重新计算一遍，与本库无关。

```python
n, w = 6, 3
chain = SpinChain(n, 'periodic')
M = chain.algebra                                                 # (+)_w M_binom(6, w), padded to 20
H_M = chain.xxz(J=1.0, Delta=0.0)
H_dense = sum(on(X, i, n) @ on(X, j, n) + on(Y, i, n) @ on(Y, j, n) for i, j in chain.bonds) / 4
spec = torch.sort(torch.cat([e[0] for e in H_M.eigenvalues()]))[0]
assert torch.allclose(spec, torch.linalg.eigvalsh(H_dense))       # same spectrum, sector by sector
edges = chain.bonds + [(0, 2), (3, 5)]                            # two triangles joined by (2,3), (5,0)
H_C = M.zero()
for i, j in edges:
    H_C = H_C + (M.identity() - chain.sz(i) @ chain.sz(j) * 4.0) * 0.5
layers = [dynamics.propagator(H_M, beta) @ dynamics.propagator(H_C, gamma)
          for gamma, beta in ((0.6, 1.0), (1.4, 1.4))]
rho0 = ideal = chain.basis_state('101010')
for U_l in layers:
    ideal = channels.unitary_channel(U_l)(ideal)
assert abs(ideal.sector_probabilities()[0, w].item() - 1) < 1e-12           # Proposition 10.8

HC_dense = sum(torch.eye(64, dtype=dt) - on(Z, i, n) @ on(Z, j, n) for i, j in edges) / 2
psi = torch.eye(64, dtype=dt)[int('101010', 2)]                   # the same circuit on 6 full qubits
for gamma, beta in ((0.6, 1.0), (1.4, 1.4)):
    psi = torch.linalg.matrix_exp(-1j * beta * H_dense) @ torch.linalg.matrix_exp(-1j * gamma * HC_dense) @ psi
weight = torch.tensor([bin(b).count('1') for b in range(64)])
assert abs(psi[weight == w].norm().item() - 1) < 1e-12           # Proposition 10.8, dense
assert abs((psi.conj() @ HC_dense @ psi).real.item() - states.expectation(ideal, H_C).real.item()) < 1e-10
```

量子比特 $`i`$ 的**振幅阻尼**（以概率 $`g`$ 发生 $`|1\rangle\to|0\rangle`$）具有 Kraus 算子
$`K_0^{(i)}=|0\rangle\langle0|_i+\sqrt{1-g}\,|1\rangle\langle1|_i`$ 和 $`K_1^{(i)}=\sqrt g\,|0\rangle\langle1|_i`$，后者降低
重量：这是扇区之间的映射，即 `chain.site_amplitude_damping(i, g)`（一个 `InterSectorChannel`）。

由[讲义，命题 10.9](../../lecture_notes/zh/10_quantum_computing.md#r-10-9)，若线路的每一层之后都对每个量子比特施加阻尼，则经过 $`p`$ 层后仍以概率
$`(1-g)^{wp}`$ 测得初态的重量 $`w`$，而后选择得到的态就是无噪声的态。

`InterSectorChannel` 只保留其输出的对角块；这里不存在其他块（命题 10.9 的证明逐字适用于单个量子比特的阻尼，
此时 $`s\in\{0,1\}`$），所以没有任何损失。取 $`g=0.05`$、$`p=2`$，接受概率为 $`0.95^6=0.7351`$：

```python
g = 0.05
damping = [chain.site_amplitude_damping(i, g) for i in range(n)]
noisy = rho0
for U_l in layers:
    noisy = channels.unitary_channel(U_l)(noisy)
    for D in damping:                                             # one qubit after the other
        noisy = D(noisy)
print(noisy.sector_probabilities()[0])                            # weight leaks to w = 2, 1, 0
P_w = M.central([1.0 if v == w else 0.0 for v in range(n + 1)])
post, accept = noisy.condition_on(P_w)                            # measure the weight, keep w = 3
cut = lambda r: states.expectation(r, H_C).real.item()
print(f"P(weight 3) {accept.item():.4f}, fidelity {noisy.fidelity(ideal).item():.4f}, "
      f"post-selected {post.fidelity(ideal).item():.4f}; "
      f"cut: ideal {cut(ideal):.4f}, noisy {cut(noisy):.4f}, post-selected {cut(post):.4f}")
assert abs(accept.item() - (1 - g) ** (w * len(layers))) < 1e-12              # Proposition 10.9
assert abs(post.fidelity(ideal).item() - 1) < 1e-10 and abs(cut(post) - cut(ideal)) < 1e-10
```

*数值观察*（打印出的值）：不可行部分使期望切割数从 $`4.3789`$ 降为 $`4.3123`$，这是一种虚假的改进，
会误导角度的优化器；后选择之后，该值是精确的（命题 10.9）。

由[讲义，命题 10.10](../../lecture_notes/zh/10_quantum_computing.md#r-10-10)，在带振幅阻尼（速率 $`\gamma_1`$）与退相位 $`L_k=S^z_k`$（速率 $`\gamma_2`$）的
Lindblad 方程下，从 $`H_w`$ 出发的态以概率 $`e^{-\gamma_1wt}`$ 保持重量 $`w`$，且后选择的态不依赖于 $`\gamma_1`$；
当 $`\gamma_2=0`$ 时后选择是完美的。配套论文
（`examples/paper/ex9_qaoa_xy_mixer.py`）让 $`|111000\rangle`$ 在 $`H_M`$ 下演化到 $`T=2`$；我们重现其中
$`\gamma_1=\gamma_2=0.1`$ 那一行：$`P(w=3)=e^{-0.6}=0.5488`$（命题 10.10），而作为下面代码的*数值观察*，保真度为 $`0.4805`$，
后选择后的保真度为 $`0.8755`$；由命题 10.10，后者与 $`\gamma_1=0`$ 时相同。

```python
T = 2.0
jumps = [chain.lowering(i) for i in range(n)] + [chain.sz(i) for i in range(n)]
rho0 = chain.basis_state('111000')                                # the initial state of the paper
ideal_T = dynamics.von_neumann(rho0, H_M, [T])[0]
post_T = {}
for g1, g2 in ((0.1, 0.0), (0.0, 0.1), (0.1, 0.1)):          # the last one is printed
    rho_T = dynamics.lindblad_evolve(rho0, H_M, jumps, [0.0, T], rates=[g1] * n + [g2] * n, substeps=80)[-1]
    post_T[g1, g2], accept_T = rho_T.condition_on(P_w)
    assert abs(accept_T.item() - math.exp(-w * g1 * T)) < 1e-8                 # Proposition 10.10
print(f"P(weight 3) {accept_T.item():.4f}, fidelity {rho_T.fidelity(ideal_T).item():.4f}, "
      f"post-selected {post_T[0.1, 0.1].fidelity(ideal_T).item():.4f}")
assert abs(post_T[0.1, 0.0].fidelity(ideal_T).item() - 1) < 1e-6               # gamma_2 = 0: perfect
assert post_T[0.1, 0.1].trace_distance(post_T[0.0, 0.1]).item() < 1e-6       # independent of gamma_1
```

<a id="s-10-5"></a>

## 10.5 可与之比较的工具

Qiskit（`qiskit.quantum_info`、Aer 噪声模型、带自定义混合算子的 `QAOAAnsatz`）、Cirq、PennyLane
（`qml.qaoa.xy_mixer`、`default.mixed`）以及 QuTiP/`qutip-qip` 在 $`n`$ 个完整的量子比特上模拟上述全部内容，
是处理线路的合适工具；stim 对带 Pauli 型噪声的 Clifford 线路进行采样，可处理数千个
量子比特（振幅阻尼不是 Pauli 信道：Pauli 信道是保单位的，而振幅阻尼不是）。据我们所知，这些软件包中没有任何函数能
计算信道的相互作用代数、它的交换子或它的整个不动点代数（QuTiP 的
`steadystate` 返回 Lindblad 生成元的一个稳态）；关于矩阵 $`*`$-代数的数值块分解，见
Murota et al., *Japan J. Indust. Appl. Math.* 27 (2010)。本库增加的是：在
$`\bigoplus_wM_{\binom nw}`$ 中的存储、扇区之间的映射以及作为中心投影的后选择，并且都是批量化的。
一个扇区需要 $`\binom nw^2`$ 个元素，而不是 $`4^n`$ 个（对 $`n=16`$，
$`w=8`$，分别为 $`1.7\cdot10^8`$ 与 $`4.3\cdot10^9`$），但整个代数会被填充到 $`(n+1)\binom n{\lfloor n/2\rfloor}^2`$ 个元素（$`2.8\cdot10^9`$），
节省有限。它不是线路模拟器，也没有稳定子形式体系。

<a id="x-10"></a>

## 习题

讲义[习题 1–3](../../lecture_notes/zh/10_quantum_computing.md#x-10) 的数值验证：对单个量子比特上的振幅阻尼，$`\mathcal A'`$、$`\mathrm{Fix}(\Phi)`$ 与
$`\mathrm{Fix}(\Phi^*)`$ 都是一维的（习题 1）；对[推论 10.3](../../lecture_notes/zh/10_quantum_computing.md#r-10-3) 中取 $`\theta=\pi/2`$、$`n=4`$ 的信道，
$`\dim\mathrm{Fix}(\Phi)=128`$（习题 2）；比特翻转码的恢复 $`\mathcal R`$ 能纠正相干误差 $`e^{-i\phi X_1}`$（习题 3）。

<details>
<summary>解答</summary>

```python
E_ad = [torch.tensor([[1, 0], [0, math.sqrt(0.7)]], dtype=dt), torch.tensor([[0, math.sqrt(0.3)], [0, 0]], dtype=dt)]
fix_ad = fixed_points(E_ad)
assert len(fix_ad) == 1 and len(commutant_basis(E_ad + [E.mH for E in E_ad])) == 1
assert torch.allclose(fix_ad[0] / fix_ad[0][0, 0], torch.tensor([[1, 0], [0, 0]], dtype=dt))
assert len(fixed_points([E.mH for E in E_ad])) == 1                     # Fix(Phi^*) = C 1
```

```python
U2 = torch.linalg.matrix_exp(-1j * (math.pi / 2) * sum(on(Z, j, 4) for j in range(4)))
assert len(fixed_points([U2 / math.sqrt(2), U2.mH / math.sqrt(2)])) == 128
assert sum(math.comb(4, v) ** 2 for v in range(5)) == math.comb(8, 4) == 70
```

```python
for phi in (0.2, 1.0, 2.5):
    G1 = torch.linalg.matrix_exp(-1j * phi * on(X, 0, 3))
    out = (R @ channels.kraus_channel([op8(G1)]))(rho_L)
    assert (out.matrix - rho_L.matrix).abs().max() < 1e-12       # corrected exactly
```

</details>

**下一章：** [教程 11](11_optimisation.md)：代数的酉群上的优化，以及迹不等式。
