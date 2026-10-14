[English](../12_equivariant_learning.md) · [Español](../es/12_equivariant_learning.md) · [Français](../fr/12_equivariant_learning.md) · [Deutsch](../de/12_equivariant_learning.md) · **中文** · [日本語](../ja/12_equivariant_learning.md) · [Русский](../ru/12_equivariant_learning.md)

# 教程 12. 等变机器学习

本教程配合讲义的[第 12 章](../../lecture_notes/zh/12_equivariant_learning.md)。该章证明：有限群的两个表示之间的等变线性映射是分块矩阵
（[定理 12.4](../../lecture_notes/zh/12_equivariant_learning.md#r-12-4)），在群上取平均就是到这些映射上的正交投影（[定理 12.5](../../lecture_notes/zh/12_equivariant_learning.md#r-12-5)），并对不变量进行计数。
这里我们用 `torch_vn_algebra.decompose` 构造并训练等变层，对象是作用在图的顶点特征与边特征上的 $`S_4`$，
并数值地检验这些结果。

在代码中，一批向量是形状为 `(batch, n)` 的张量，矩阵 $`T`$ 以 `x @ T.T` 的方式作用于其上。

<a id="s-12-1"></a>

## 12.1 等变线性映射

*讲义：* [12.1 节](../../lecture_notes/zh/12_equivariant_learning.md#s-12-1)：[引理 12.1](../../lecture_notes/zh/12_equivariant_learning.md#r-12-1)（轨道）与[推论 12.2](../../lecture_notes/zh/12_equivariant_learning.md#r-12-2)（DeepSets、Bell 数）。

DeepSets 层（Zaheer et al. 2017）就是[讲义，推论 12.2](../../lecture_notes/zh/12_equivariant_learning.md#r-12-2)(a)，Maron et al.（2019）关于图网络的基就是 (b)。
下面的代码用穷举法计数轨道，与数值计算出的交换子作比较（[引理 12.1](../../lecture_notes/zh/12_equivariant_learning.md#r-12-1)），
并分解 $`S_4`$ 的这两个表示：

```python
import itertools
import torch
import torch_vn_algebra.decompose as dc

torch.manual_seed(0)
dt = torch.complex128
n = 4
perms = list(itertools.permutations(range(n)))                   # S_4, |G| = 24

def perm_matrix(p):
    """pi(sigma) e_x = e_{sigma(x)} on C^n."""
    P = torch.zeros(n, n, dtype=dt)
    P[list(p), list(range(n))] = 1
    return P

G_node = torch.stack([perm_matrix(p) for p in perms])            # S_4 on C^4: vertex features
G_edge = torch.stack([torch.kron(P, P) for P in G_node])         # S_4 on C^16: X -> P X P^T, row-major

def n_orbits(tuples):
    """Number of orbits of S_4 on tuples of vertices: the canonical form is the smallest image."""
    return len({min(tuple(p[i] for i in t) for p in perms) for t in tuples})

bell = {r: n_orbits(itertools.product(range(n), repeat=r)) for r in (2, 3, 4)}
print(bell)
assert bell == {2: 2, 3: 5, 4: 15}                               # Corollary 12.2(b)
for G, r in ((G_node, 2), (G_edge, 4)):                          # Lemma 12.1: dim pi(G)' = number of orbits
    assert dc.commutant_basis(G, group=True).shape[0] == dc.commutant_dimension(G) == bell[r]
for X in dc.commutant_basis(G_node, group=True):                 # Corollary 12.2(a): X = a 1 + b (J - 1)
    off = X[~torch.eye(n, dtype=torch.bool)]
    assert (X.diagonal() - X[0, 0]).abs().max() < 1e-12 and (off - off[0]).abs().max() < 1e-12

D_node, D_edge = dc.from_group(G_node), dc.from_group(G_edge)
print(D_node)
print(D_edge)
assert (D_node.k_factors, D_node.multiplicities) == ((1, 3), (1, 1))
assert (D_edge.k_factors, D_edge.multiplicities) == ((1, 2, 3, 3), (2, 1, 1, 3))
```

所以，在 `dc.from_group` 计算出的分解中（断言检验了它的输出），$`\mathbb C^{16}`$ 包含两次维数为 $`1`$ 的
表示——即平凡表示，由命题 12.7(a) 与引理 12.1，它的重数等于 $`[4]^2`$ 上的轨道个数 $`2`$（12.4 节的代码
识别出这个扇区）——一次维数为 $`2`$ 的表示，以及两个维数为 $`3`$ 的表示，分别出现一次和三次：
$`16=1\cdot2+2\cdot1+3\cdot1+3\cdot3`$，而 $`15=2^2+1^2+1^2+3^2`$ 与引理 12.1 的轨道计数以及
$`\dim\pi(G)'=\sum_cm_c^2`$（定理 2.6）相符。对于两个表示，一切都归结为一个交换子：

由[讲义，引理 12.3](../../lecture_notes/zh/12_equivariant_learning.md#r-12-3)，对 $`\pi=\pi_1\oplus\pi_2`$ 有 $`\mathrm{Hom}_G(\pi_1,\pi_2)=Q_2^*\,\pi(G)'\,Q_1`$；而由
[定理 12.4](../../lecture_notes/zh/12_equivariant_learning.md#r-12-4)，在适当的基下它就是 $`\bigoplus_cM_{b_c\times a_c}(\mathbb C)`$，维数为 $`\sum_\rho m_\rho(\pi_1)\,m_\rho(\pi_2)`$。

<a id="s-12-2"></a>

## 12.2 用库构造等变层

*讲义：* [12.2 节](../../lecture_notes/zh/12_equivariant_learning.md#s-12-2)：[定理 12.5](../../lecture_notes/zh/12_equivariant_learning.md#r-12-5)（twirl）；这些层遵循[引理 12.3](../../lecture_notes/zh/12_equivariant_learning.md#r-12-3)与[定理 12.4](../../lecture_notes/zh/12_equivariant_learning.md#r-12-4)的证明。

一旦知道了定理 2.6 的分解，定理 12.4 的证明就是一个算法；`dc.from_group` 用 4.3 节的随机化方法
计算这个分解（在精确算术中以概率 1 正确；这里由断言检验计算出的输出）。当 $`\pi_1=\pi_2`$ 时，它就是定理 2.6 本身：权重是块 $`w_c\in M_{m_c}`$，即抽象交换子
`D.commutant_algebra()` $`=\bigoplus_cM_{m_c}`$ 的一个算子，而层是
$`x\mapsto V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*x`$，由 `D.from_commutant` 计算。下面的模块是本章的
代码，不属于库的一部分。

```python
class EquivariantLinear(torch.nn.Module):
    """x -> T x, T = V ((+)_c 1_{k_c} (x) w_c) V^* in pi(G)'; one weight block w_c of size m_c x m_c per sector."""
    def __init__(self, D):
        super().__init__()
        self.D, self.CA = D, D.commutant_algebra()
        self.w = torch.nn.ParameterList(torch.nn.Parameter(torch.randn(m, m, dtype=dt) / m ** 0.5)
                                        for m in D.multiplicities)

    def matrix(self):
        return self.D.from_commutant(self.CA.from_blocks(list(self.w)))[0]

    def forward(self, x):                                        # x: (batch, n)
        return x @ self.matrix().T

def equivariance_error(layer, G_in, G_out):
    """max_g |layer(pi_in(g) x) - pi_out(g) layer(x)| on a random batch."""
    x = torch.randn(8, G_in.shape[-1], dtype=dt)
    with torch.no_grad():
        return max((layer(x @ g.T) - layer(x) @ h.T).abs().max().item() for g, h in zip(G_in, G_out))

edge_layer = EquivariantLinear(D_edge)
n_params = sum(p.numel() for p in edge_layer.parameters())
err = equivariance_error(edge_layer, G_edge, G_edge)
print(f"edge -> edge: {n_params} parameters, equivariance error {err:.1e}")
assert n_params == D_edge.commutant_dim == bell[4] and err < 1e-12
```

当 $`\pi_1\ne\pi_2`$ 时，逐行照搬证明：分解 $`\pi_1\oplus\pi_2`$，用 `to_commutant` 计算 $`P_1`$ 的块 $`p_c`$，
将它们对角化，并保留左下角。$`p_c`$ 的迹是 $`a_c=m_c(\pi_1)`$。

对 $`n=4`$ 时的推论 12.2(b)，在计算出的分解中（下面打印并断言）：在 $`\pi_1\oplus\pi_2`$ 的四个扇区（维数
为 $`1,2,3,3`$；第一个是平凡表示，因为由推论 12.2(a)，$`\pi_1`$ 唯一的一维成分是平凡表示）中，顶点特征
$`\pi_1`$ 有 $`a=(1,0,0,1)`$，边特征 $`\pi_2`$ 有 $`b=(2,1,1,3)`$，所以
$`\dim\mathrm{Hom}_G=1\cdot2+1\cdot3=5=B_3`$，正如推论 12.2(b) 所证明的。

```python
class EquivariantMap(torch.nn.Module):
    """x -> T x with T in Hom_G(pi_1, pi_2) (Theorem 12.4): one weight block y_c of size b_c x a_c per sector."""
    def __init__(self, G1, G2):
        super().__init__()
        self.n1 = n1 = G1.shape[-1]
        self.D = dc.from_group(torch.stack([torch.block_diag(g1, g2) for g1, g2 in zip(G1, G2)]))
        P1 = torch.block_diag(torch.eye(n1, dtype=dt), torch.zeros(G2.shape[-1], G2.shape[-1], dtype=dt))
        p = self.D.to_commutant(P1).matrix[0]                    # blocks p_c of P_1, padded
        self.u, self.a = [], []
        for c, m in enumerate(self.D.multiplicities):
            lam, u = torch.linalg.eigh(p[c, :m, :m])             # eigenvalues 0, ..., 0, 1, ..., 1
            self.a.append(round(lam.sum().item()))               # a_c = rank p_c = m_c(pi_1)
            self.u.append(u.flip(-1))                            # u_c^* p_c u_c = 1_{a_c} (+) 0_{b_c}
        self.y = torch.nn.ParameterList(torch.nn.Parameter(torch.randn(m - a, a, dtype=dt))
                                        for m, a in zip(self.D.multiplicities, self.a))

    def matrix(self):
        blocks = []
        for u, a, y in zip(self.u, self.a, self.y):
            x = torch.zeros(u.shape, dtype=dt)
            x[a:, :a] = y                                        # kappa_c y_c iota_c^*
            blocks.append(u @ x @ u.mH)
        T = self.D.from_commutant(self.D.commutant_algebra().from_blocks(blocks))[0]
        return T[self.n1:, :self.n1]                             # Q_2^* (...) Q_1

    def forward(self, x):
        return x @ self.matrix().T

node_to_edge = EquivariantMap(G_node, G_edge)
n_params = sum(p.numel() for p in node_to_edge.parameters())
err = equivariance_error(node_to_edge, G_node, G_edge)
print(f"node -> edge: sectors k = {node_to_edge.D.k_factors}, a = {node_to_edge.a}, "
      f"{n_params} parameters, equivariance error {err:.1e}")
assert n_params == bell[3] and err < 1e-12
b = tuple(m - a for m, a in zip(node_to_edge.D.multiplicities, node_to_edge.a))
assert node_to_edge.D.k_factors == (1, 2, 3, 3) and node_to_edge.a == [1, 0, 0, 1] and b == (2, 1, 1, 3)
corners = node_to_edge.D.basis('commutant')[:, n:, :n].reshape(-1, 16 * n)   # Lemma 12.3: Q_2^* pi(G)' Q_1
assert torch.linalg.matrix_rank(corners).item() == bell[3]
```

使层等变的第二种方法是在群上对它取平均。**twirl** 定义为
```math
\mathcal P(W)=\frac1{|G|}\sum_{g\in G}\pi_2(g)\,W\,\pi_1(g)^*,\qquad W\in M_{n_2\times n_1}(\mathbb C).
```

由[讲义，定理 12.5](../../lecture_notes/zh/12_equivariant_learning.md#r-12-5)，$`\mathcal P`$ 是到 $`\mathrm{Hom}_G(\pi_1,\pi_2)`$ 上关于 Hilbert–Schmidt 内积的正交投影，并且
$`\dim\mathrm{Hom}_G(\pi_1,\pi_2)=|G|^{-1}\sum_g\overline{\chi_{\pi_1}(g)}\chi_{\pi_2}(g)`$；当 $`\pi_1=\pi_2`$ 时，`dc.commutant_dimension` 计算的正是这个
特征标公式。参数化 $`W\mapsto\mathcal P(W)`$ 有 $`n_1n_2-\dim\mathrm{Hom}_G`$ 个参数浪费在核上（对边特征为 $`256-15`$ 个），并且对稠密的 $`\pi(g)`$，
一次求值的代价为 $`O(|G|n_1n_2(n_1+n_2))`$，而对 $`S_n`$ 有 $`|G|=n!`$。库从分解出发计算同一个投影，即
`D.conditional_expectation(W, 'commutant')`，只需四次 $`n\times n`$ 矩阵乘积。
代码检验 `dc.twirl` 与 `D.conditional_expectation` 一致，$`\mathcal P`$ 在 $`\pi(G)'`$ 上是恒等映射
且是自伴的，并检验从顶点特征到边特征的维数公式。

```python
W = torch.randn(16, 16, dtype=dt)
PW = dc.twirl(G_edge, W)
S = torch.randn(16, 16, dtype=dt)
hs = lambda A, B: (A.conj() * B).sum()                           # Hilbert-Schmidt inner product
T_layer = edge_layer.matrix().detach()
assert (PW - D_edge.conditional_expectation(W, 'commutant')[0]).abs().max() < 1e-12
assert (dc.twirl(G_edge, T_layer) - T_layer).abs().max() < 1e-12           # P = id on pi(G)'
assert abs(hs(S, PW) - hs(dc.twirl(G_edge, S), W)) < 1e-10                  # P is self-adjoint
assert abs(hs(T_layer, W - PW)) < 1e-10                                     # W - P(W) is orthogonal to pi(G)'
chi_node, chi_edge = (G.diagonal(dim1=-2, dim2=-1).sum(-1) for G in (G_node, G_edge))
assert abs((chi_node.conj() * chi_edge).mean() - bell[3]) < 1e-12           # Theorem 12.5: 5
```

<a id="s-12-3"></a>

## 12.3 带约束的层与谱层

*讲义：* [12.3 节](../../lecture_notes/zh/12_equivariant_learning.md#s-12-3)：[命题 12.6](../../lecture_notes/zh/12_equivariant_learning.md#r-12-6)。

`EquivariantLinear` 的权重 $`w=(w_c)_c`$ 位于一个 $`*`$-代数中，所以第 1 章与第 11 章的每个构造都适用于
它们。对 $`w\in\bigoplus_cM_{m_c}`$ 记 $`\Phi(w)=V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*`$。

由[讲义，命题 12.6](../../lecture_notes/zh/12_equivariant_learning.md#r-12-6)，(a) $`\Phi`$ 是从 $`\bigoplus_cM_{m_c}`$ 到 $`\pi(G)'`$ 上的 $`*`$-同构，并且 $`\Phi(w)`$ 是
自伴的、酉的或正的，当且仅当每个 $`w_c`$ 都是；(b) 对正规的 $`w`$ 有 $`f(\Phi(w))=\Phi\big((f(w_c))_c\big)`$；(c)
等变酉算子恰为 $`\Phi(u_0e^{s-s^*})`$。

所以正交（酉）等变层就是在 `D.commutant_algebra()` 上的 `optimize.UnitaryParam` 再接上 $`\Phi`$；由 (c)，
这个参数化能取到每个等变酉算子，并且在精确算术中，梯度方法的每一步迭代都是酉的且等变的（代码在训练
之后以 $`10^{-12}`$ 的精度检验这两点）。类似地，
`PositiveParam` 给出正的等变层，其归一化如第 11 章（谱在 $`[\ell,1]`$ 中且最大值为 $`1`$；由 (a)，一般的
正层是 $`\Phi(t^*t)`$）；而自伴等变生成元 $`L`$ 的谱层 $`f(L)`$，例如扩散 $`e^{-tL}`$，需要对各块 $`w_c`$ 作特征
分解，代价的量级为 $`\sum_cm_c^3`$（这里为 $`37`$），而不是 $`n^3=4096`$；库对 `apply_function` 求导（对其反向
传播我们不作任何证明；下面的代码检验了一个导数）。下面的代码用 Adam 把一个等变酉算子拟合到 $`32`$ 对
输入–输出数据上——这是本次运行的*数值观察*，而不是定理：损失下降超过 $`10^3`$ 倍，这由断言检验，我们
不声称任何收敛性结论——并将 (b) 以及关于 $`t`$ 的梯度与稠密矩阵和公式 $`\frac d{dt}e^{-tL}=-Le^{-tL}`$ 对照
检验（对 $`L=U\mathrm{diag}(\lambda)U^*`$，由定理 1.5 有 $`e^{-tL}=U\mathrm{diag}(e^{-t\lambda_i})U^*`$，
且 $`U`$ 不依赖于 $`t`$）。

```python
from torch_vn_algebra import optimize as opt

CA = D_edge.commutant_algebra()                                  # M_2 (+) M_1 (+) M_1 (+) M_3
U_target = D_edge.from_commutant(CA.random_unitary_operator(1))[0]
X = torch.randn(32, 16, dtype=dt)
Y = X @ U_target.T

U = opt.UnitaryParam(CA, 1)                                      # u_0 Haar-random, s = 0
adam = torch.optim.Adam(U.parameters(), lr=0.05)
for step in range(150):
    U_full = D_edge.from_commutant(U.operator())[0]
    loss = (X @ U_full.T - Y).abs().pow(2).mean()
    if step == 0:
        loss0 = loss.item()
    adam.zero_grad()
    loss.backward()
    adam.step()
U_full = D_edge.from_commutant(U.operator())[0].detach()
unit_err = (U_full.mH @ U_full - torch.eye(16, dtype=dt)).abs().max().item()
equi_err = max((g @ U_full - U_full @ g).abs().max().item() for g in G_edge)
print(f"loss {loss0:.3f} -> {loss.item():.1e};  |U*U - 1| = {unit_err:.1e}, |[U, pi(g)]| = {equi_err:.1e}")
assert loss.item() < 1e-3 * loss0 and unit_err < 1e-12 and equi_err < 1e-12

h = torch.randn(1, CA.C, CA.k_max, CA.k_max, dtype=dt)
h = CA.operator(h + h.mH, is_self_adjoint=True)                  # a self-adjoint generator h = (h_c)
L = D_edge.from_commutant(h)[0]
t = torch.tensor(0.7, dtype=torch.float64, requires_grad=True)
E = D_edge.from_commutant(h.apply_function(lambda w: torch.exp(-t * w)))[0]   # e^{-tL}, in the small algebra
Z = torch.randn(16, 16, dtype=dt)
(Z.conj() * E).sum().real.backward()                             # d/dt Re <Z, e^{-tL}>
dE = -L @ E.detach()
assert (E - torch.linalg.matrix_exp(-0.7 * L)).abs().max() < 1e-12          # Proposition 12.6(b)
assert abs(t.grad.item() - (Z.conj() * dE).sum().real.item()) < 1e-10
```

<a id="s-12-4"></a>

## 12.4 不变量

*讲义：* [12.4 节](../../lecture_notes/zh/12_equivariant_learning.md#s-12-4)：[命题 12.7](../../lecture_notes/zh/12_equivariant_learning.md#r-12-7)（不变向量与二次型）与[命题 12.8](../../lecture_notes/zh/12_equivariant_learning.md#r-12-8)（Molien 公式）。

对于图，不变边特征由对角线及其补集的示性向量张成（$`B_2=2`$），不变的线性读出是对角元之和与非对角元
之和。下面的代码从块 $`\rho_c(g)`$（`to_algebra`）中找出平凡扇区，并检验 (a)、(b) 以及 (c) 的两种计数。

```python
rho = D_edge.to_algebra(G_edge).matrix                           # blocks rho_c(g), shape (24, C, 3, 3)
c0 = [c for c, k in enumerate(D_edge.k_factors) if k == 1 and (rho[:, c, 0, 0] - 1).abs().max() < 1e-12]
P0 = G_edge.mean(0)
assert c0 == [0] and (P0 - D_edge.central_projections()[0]).abs().max() < 1e-12      # (a)
assert abs(chi_edge.mean() - D_edge.multiplicities[0]) < 1e-12
diag = torch.eye(n, dtype=dt).reshape(16)
assert (P0 @ diag - diag).abs().max() < 1e-12 and (P0 @ (1 - diag) - (1 - diag)).abs().max() < 1e-12

w = [torch.randn(m, m, dtype=dt) for m in D_edge.multiplicities]
A = D_edge.from_commutant(CA.from_blocks([b + b.mH for b in w]))[0]                   # (b)
x = torch.randn(16, dtype=dt)
q = lambda v: (v.conj() @ A @ v).real
assert (A - A.mH).abs().max() < 1e-12
assert max(abs(q(g @ x) - q(x)).item() for g in G_edge) < 1e-12

def degree2(G):                                                  # (c): characters
    chi, chi_sq = (H.diagonal(dim1=-2, dim2=-1).sum(-1) for H in (G, G @ G))
    return ((chi ** 2 + chi_sq) / 2).mean().real.item()

pairs = list(itertools.product(range(n), repeat=2))
multisets = [(a, b) for i, a in enumerate(pairs) for b in pairs[i:]]
orbits2 = len({min(tuple(sorted(tuple(p[v] for v in e) for e in (a, b))) for p in perms) for a, b in multisets})
print(f"invariant quadratic polynomials: vertices {degree2(G_node):.4f}, edges {degree2(G_edge):.4f} "
      f"(orbits of S_4 on monomials: {orbits2})")
assert abs(degree2(G_node) - 2) < 1e-12 and abs(degree2(G_edge) - orbits2) < 1e-12 and orbits2 == 11
```

在 $`\mathbb C^4`$ 上，它们是 $`\sum_ix_i^2`$ 与 $`\sum_{i<j}x_ix_j`$。(c) 的证明对每个次数都适用：

由[讲义，命题 12.8](../../lecture_notes/zh/12_equivariant_learning.md#r-12-8)（Molien 公式），
$`\sum_{d\ge0}\dim\mathbb C[x]^G_d\,t^d=|G|^{-1}\sum_g\det\big(1-t\pi(g)\big)^{-1}`$，并且对置换表示，
$`\dim\mathbb C[x]^G_d`$ 等于 $`G`$ 在 $`d`$ 次单项式上的轨道个数。

（这是 1897 年的 Molien 定理；延伸阅读参见 B. Sturmfels, *Algorithms in Invariant Theory*, Springer 2008，
2.2 节。）下面的代码把右边展开到 $`t^4`$，对每个特征值乘上一个几何级数，并与单项式的轨道作比较；对作用在
$`\mathbb C^4`$ 上的 $`S_4`$，计数 $`1,1,2,3,5`$ 就是把 $`d`$ 分拆为至多四个部分的分拆个数（多重集的轨道由其
排序后的重数决定）。

```python
def molien(G, d_max):
    """Coefficients of |G|^{-1} sum_g det(1 - t pi(g))^{-1} up to t^d_max: products of geometric series."""
    lam = torch.linalg.eigvals(G)                                # (|G|, n), |lambda| = 1
    c = torch.zeros(G.shape[0], d_max + 1, dtype=dt)
    c[:, 0] = 1
    for i in range(lam.shape[1]):                                # multiply by 1 / (1 - t lambda_i)
        for d in range(1, d_max + 1):
            c[:, d] = c[:, d] + lam[:, i] * c[:, d - 1]
    return c.mean(0)

def monomial_orbits(G, d):
    """Orbits of G on the monomials of degree d, i.e. on the multisets of d points of Omega."""
    act = [g.real.argmax(0).tolist() for g in G]                 # pi(g) e_x = e_{act[x]}
    mons = itertools.combinations_with_replacement(range(G.shape[-1]), d)
    return len({min(tuple(sorted(a[x] for x in m)) for a in act) for m in mons})

for G in (G_node, G_edge):
    series, counts = molien(G, 4), [monomial_orbits(G, d) for d in range(5)]
    print("Molien series:", [f"{v.real:.4f}" for v in series], " orbits of monomials:", counts)
    assert (series - torch.tensor(counts, dtype=dt)).abs().max() < 1e-9          # Proposition 12.8
assert [monomial_orbits(G_node, d) for d in range(5)] == [1, 1, 2, 3, 5] and monomial_orbits(G_edge, 2) == orbits2
```

<a id="s-12-5"></a>

## 12.5 可与之比较的工具

针对特定群的等变层已有成熟的库提供。**escnn**（Cesa、Lang 与 Weiler 2022，e2cnn 的后继者）由不可约
表示和核的调和基为 $`O(2)`$ 与 $`O(3)`$ 的子群构造可操控卷积；**e3nn**（Geiger 与 Smidt 2022）由球谐函数
与 Clebsch–Gordan 张量积构造 $`E(3)`$-等变网络。**EMLP**（Finzi、Welling 与 Wilson 2021）处理任意矩阵群：
它对离散生成元施加 $`\pi_2(g)T=T\pi_1(g)`$，对 Lie 代数生成元施加 $`d\pi_2(X)T=Td\pi_1(X)`$，并把解空间的
一组基作为数值零空间算出——这正是 `dc.commutant_basis` 所解的线性问题。群卷积（Cohen 与 Welling 2016）
是引理 12.1 中 $`\Omega_1=\Omega_2=G`$ 的情形，其中 $`G`$ 通过左乘作用在自身上：$`(y,x)`$ 与 $`(y',x')`$ 属于同一
轨道当且仅当 $`x^{-1}y=x'^{-1}y'`$（若如此，$`g=x'x^{-1}`$ 把一对映为另一对），所以等变映射就是
$`(Tv)(y)=\sum_xf(x^{-1}y)v(x)`$，即与 $`G`$ 上的函数 $`f`$ 的卷积。Kondor 与 Trivedi（2018）研究了紧群的齐性空间
上的类似问题，本课程不讨论这一问题。对置换而言，DeepSets（Zaheer et al. 2017）与 Maron et al.（2019）的基就是推论 12.2。

本库不提供网格上的卷积、连续空间上的可操控核、$`\mathbb R`$ 上的分解（实型、复型与四元数型；习题 2
处理实矩阵），也不适用于大的 $`n`$：$`V`$ 是稠密的 $`n\times n`$ 矩阵，并且对由全部元素给出的群，分解的
代价为每个随机样本 $`O(|G|n^3)`$（Lie 代数生成元通过 `dc.from_generators` 处理；4.3 节）。它所增加的
是解空间的**结构**：交换子的一组基只是一列矩阵，而分解把交换子与 $`*`$-代数 $`\bigoplus_cM_{m_c}`$
等同起来。于是乘积、伴随、酉性、正性和函数演算都在小矩阵上逐块作用（命题 12.6），库可以对它们求导，
而 `to_commutant`/`from_commutant` 在两种图景之间转换。

<a id="x-12"></a>

## 习题

讲义[习题](../../lecture_notes/zh/12_equivariant_learning.md#x-12) 1、3 和 4 的数值检验：

- *习题 1*（带通道的 DeepSets）：用 `EquivariantMap` 验证从
  $`\mathbb C^n\otimes\mathbb C^{d_1}`$ 到 $`\mathbb C^n\otimes\mathbb C^{d_2}`$ 的等变映射恰为 $`1\otimes A+J\otimes B`$。
- *习题 3*（非线性）：逐坐标的 ReLU 与置换表示可交换，但与 $`\mathbb R^2`$ 上转角为 $`\pi/2`$ 的旋转
  不可交换；$`x\mapsto\varphi(\|x\|)x`$ 与每个酉表示都可交换。
- *习题 4*：对边表示，一个高斯样本会留下一个 $`2`$ 维的等变层空间无法确定，两个样本则能确定该层。

<details>
<summary>解答</summary>

*习题 1.*

```python
G_in = torch.stack([torch.kron(P, torch.eye(2, dtype=dt)) for P in G_node])     # d_1 = 2
G_out = torch.stack([torch.kron(P, torch.eye(3, dtype=dt)) for P in G_node])    # d_2 = 3
deepsets = EquivariantMap(G_in, G_out)
assert sum(p.numel() for p in deepsets.parameters()) == 2 * 2 * 3
assert equivariance_error(deepsets, G_in, G_out) < 1e-12
T = deepsets.matrix().detach().reshape(n, 3, n, 2)                # blocks T[i, :, j, :] = delta_ij A + B
A_, B_ = T[0, :, 0] - T[0, :, 1], T[0, :, 1]
assert all((T[i, :, j] - (A_ if i == j else 0) - B_).abs().max() < 1e-12 for i in range(n) for j in range(n))
```

*习题 3.*

```python
relu = lambda v: v.clamp(min=0)
R = torch.tensor([[0.0, -1.0], [1.0, 0.0]], dtype=torch.float64)
v = torch.tensor([1.0, -1.0], dtype=torch.float64)
assert (relu(R @ v) - R @ relu(v)).abs().max() == 1                             # (b)
xr = torch.randn(16, dtype=torch.float64)
assert all(torch.equal(relu(g.real @ xr), g.real @ relu(xr)) for g in G_edge)   # (a): exact
gate = lambda v: torch.tanh(v.norm()) * v
assert max((gate(g @ x) - g @ gate(x)).abs().max().item() for g in G_edge) < 1e-12   # (c)
```

*习题 4.*

```python
basis = D_edge.basis('commutant')                                # (15, 16, 16), Hilbert-Schmidt orthonormal
for N, lost in ((1, 2), (2, 0)):
    Xs = torch.randn(16, N, dtype=dt)
    M = torch.stack([(B @ Xs).reshape(-1) for B in basis], dim=1)               # w -> (T x_1, ..., T x_N)
    assert torch.linalg.matrix_rank(M).item() == 15 - lost
```

</details>

**下一章：** [教程 13. 若尔当代数：量子力学为何使用 *-代数](13_jordan_algebras.md)。
