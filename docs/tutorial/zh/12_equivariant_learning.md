[English](../12_equivariant_learning.md) · [Español](../es/12_equivariant_learning.md) · [Français](../fr/12_equivariant_learning.md) · [Deutsch](../de/12_equivariant_learning.md) · **中文** · [日本語](../ja/12_equivariant_learning.md) · [Русский](../ru/12_equivariant_learning.md)

# 12. 等变机器学习

神经网络的一层称为**等变的**，如果先用一个对称变换作用于输入再应用该层，与先应用该层再变换输出
得到相同的结果：对图的顶点重新编号，会以同样的方式对其特征重新编号。对线性层而言这是一个交换
关系，因此等变层构成一个交换子，由第 2 章，在相差一个基变换的意义下它就是代数 $\bigoplus_cM_{m_c}(\mathbb C)$。
本章计算它们的参数个数，用 `torch_vn_algebra.decompose` 构造并训练它们，并利用代数结构得到按构造
即为等变的酉层、正层与谱层。

**你将学到**为什么两个表示之间的等变线性映射是分块矩阵，其块的大小由重数给出；如何通过轨道或
特征标来计数它们（DeepSets 以及图网络的 Bell 数）；为什么在群上取平均就是到等变映射上的正交投影；
函数演算与酉参数化如何在交换子上逐块作用；以及如何计数不变向量与不变二次型。

在本章中，$G$ 是有限群，$G$ 在 $\mathbb C^n$ 上的**表示**是取酉值的同态
$\pi:G\to M_n(\mathbb C)$（这不是限制：可逆的 $\pi(g)$ 对内积
$\sum_g\langle\pi(g)x,\pi(g)y\rangle$ 是酉的，见 4.1 节）。它的**特征标**是 $\chi_\pi(g)=\operatorname{Tr}\pi(g)$。不可约
表示、等价表示，以及不可约表示 $\rho$ 在 $\pi$ 中的**重数** $m_\rho(\pi)$（若 $\rho$ 不出现则为 $0$）
均如 2.4 节所述，该节证明了 $m_\rho(\pi)$ 不依赖于分解，并证明了 Schur 引理；[第 4 章](04_symmetries.md)
对二者有详细讨论（4.1–4.2 节）。在长方矩阵上，$\langle S,T\rangle=\operatorname{Tr}(S^*T)$ 是
**Hilbert–Schmidt** 内积。在代码中，一批向量是形状为 `(batch, n)` 的张量，矩阵 $T$ 以 `x @ T.T`
的方式作用于其上。

## 12.1 等变线性映射

**定义。** 设 $\pi_1,\pi_2$ 是 $G$ 在 $\mathbb C^{n_1}$、$\mathbb C^{n_2}$ 上的表示。若对所有
$g\in G$ 有 $T\pi_1(g)=\pi_2(g)T$，则称矩阵 $T\in M_{n_2\times n_1}(\mathbb C)$ 是**等变的**（或称
**缠结算子**）。它们构成向量空间 $\operatorname{Hom}_G(\pi_1,\pi_2)$；当 $\pi_1=\pi_2=\pi$ 时，它就是交换子
$\operatorname{Hom}_G(\pi,\pi)=\pi(G)'$。

机器学习中最常见的表示是置换坐标的表示。设 $G$ 作用在有限集 $\Omega$ 上；$\mathbb C^\Omega$ 上的
**置换表示**为 $\pi(g)e_x=e_{gx}$。例子：$S_n$ 置换一个集合的 $n$ 个点（$\Omega=[n]=\{1,\dots,n\}$），
以及 $S_n$ 对图的顶点重新编号，它作用在 $\Omega=[n]^2$ 上，即通过 $X\mapsto P_\sigma XP_\sigma^T$ 作用在
邻接矩阵上。

**引理 12.1（轨道）。** 设 $\pi_i$ 是 $G$ 在 $\mathbb C^{\Omega_i}$ 上的置换表示。则
$T\in\operatorname{Hom}_G(\pi_1,\pi_2)$ 当且仅当对所有 $g$、$y\in\Omega_2$、$x\in\Omega_1$ 有 $T_{gy,gx}=T_{yx}$。
$G$ 在 $\Omega_2\times\Omega_1$ 上各轨道的示性矩阵构成一组基，并且
$\dim\operatorname{Hom}_G(\pi_1,\pi_2)$ 等于这些轨道的个数。

*证明。* $\pi_2(g)T\pi_1(g)^*e_x=\pi_2(g)Te_{g^{-1}x}=\sum_yT_{y,g^{-1}x}e_{gy}$，所以
$\pi_2(g)T\pi_1(g)^*$ 的 $(gy,x)$ 元为 $T_{y,g^{-1}x}$。因此 $T$ 等变，即对所有 $g$ 有
$\pi_2(g)T\pi_1(g)^*=T$，当且仅当对所有 $g,y,x$ 有 $T_{gy,x}=T_{y,g^{-1}x}$，而这（把 $x$ 换成 $gx$）
正是所述条件：$T$ 在轨道上为常数。各轨道的示性矩阵支撑互不相交，因而线性无关，并且每个在轨道上
为常数的函数都是它们的线性组合。$\square$

**推论 12.2.** (a)（*DeepSets。*）对作用在 $\mathbb C^n$ 上的 $S_n$，$n\ge2$，交换子为
$\operatorname{span}\{1,J\}$，其中 $J$ 是全 1 矩阵，并且 $\mathbb C^n=\mathbb C\mathbf 1\oplus\mathbf 1^\perp$
是两个不等价不可约表示之和：平凡表示，以及维数为 $n-1$ 的**标准**表示。(b)（*图网络。*）从
$\mathbb C^{[n]^p}$ 到 $\mathbb C^{[n]^q}$（$n$ 个顶点上的 $p$ 阶与 $q$ 阶张量）的等变映射空间的维数，
等于把一个 $p+q$ 元集合划分为至多 $n$ 块的划分个数；当 $n\ge p+q$ 时，这个数是 **Bell 数** $B_{p+q}$
（$B_2=2$，$B_3=5$，$B_4=15$）。

*证明。* (b) 两个元组 $x,y\in[n]^{p+q}$ 属于同一轨道，当且仅当它们有相同的**相等模式**，即位置集
$\{1,\dots,p+q\}$ 按 $i\sim j\Leftrightarrow x_i=x_j$ 所得的划分。模式显然是不变的。反之，若 $x,y$
有相同的模式，则 $x_i\mapsto y_i$ 是从 $\{x_i\}$ 到 $[n]$ 的一个定义良好的单射，它可以延拓为满足
$\sigma x=y$ 的置换 $\sigma$。每个至多有 $n$ 块的划分都是某个模式（给各块赋以不同的值）。由引理 12.1
即得结论。(a) 对 $p=q=1$，$[n]^2$ 上有两个轨道：对角线及其补集，示性矩阵为 $1$ 与 $J-1$。所以
$\dim\pi(S_n)'=2=\sum_cm_c^2$（定理 2.6，命题 2.7），这迫使有两个 $m_c=1$ 的扇区，即两个不等价的
不可约直和项，并且 $\pi(S_n)'$ 由它们的投影 $P_1,P_2$ 张成。到不变直线 $\mathbb C\mathbf 1$ 上的投影
属于 $\pi(S_n)'$（引理 2.1），所以它是 $P_1$ 或 $P_2$：承载平凡表示的这条直线是一个扇区，
$\mathbf 1^\perp$ 是另一个。$\square$

DeepSets 层 $x\mapsto\lambda x+\gamma\,(\mathbf 1^Tx)\mathbf 1$（Zaheer et al. 2017）就是 (a)；
Maron et al.（2019）关于图网络的基就是 (b)。下面的代码用穷举法计数轨道，与数值计算出的交换子作
比较，并分解 $S_4$ 的这两个表示：

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

所以 $\mathbb C^{16}$ 包含两次平凡表示、一次维数为 $2$ 的表示，以及两个维数为 $3$ 的表示，分别出现
一次和三次：$16=1\cdot2+2\cdot1+3\cdot1+3\cdot3$，$15=2^2+1^2+1^2+3^2$。对于两个表示，一切都
归结为一个交换子：

**引理 12.3（直和）。** 设 $\pi=\pi_1\oplus\pi_2$ 作用在 $\mathbb C^{n_1}\oplus\mathbb C^{n_2}$ 上，等距
嵌入为 $Q_i$，投影为 $P_i=Q_iQ_i^*$。则 $P_1,P_2\in\pi(G)'$，
$\operatorname{Hom}_G(\pi_1,\pi_2)=Q_2^*\,\pi(G)'\,Q_1$，并且 $T\mapsto Q_2TQ_1^*$ 是从
$\operatorname{Hom}_G(\pi_1,\pi_2)$ 到 $P_2\,\pi(G)'P_1$ 上的线性双射。

*证明。* 由于 $\pi(g)$ 是分块对角的，$\pi(g)Q_i=Q_i\pi_i(g)$ 且 $Q_i^*\pi(g)=\pi_i(g)Q_i^*$；因此
$P_i\pi(g)=Q_i\pi_i(g)Q_i^*=\pi(g)P_i$。若 $X\in\pi(G)'$，则
$Q_2^*XQ_1\pi_1(g)=Q_2^*X\pi(g)Q_1=Q_2^*\pi(g)XQ_1=\pi_2(g)Q_2^*XQ_1$。若 $T$ 等变，则
$\pi(g)Q_2TQ_1^*=Q_2\pi_2(g)TQ_1^*=Q_2T\pi_1(g)Q_1^*=Q_2TQ_1^*\pi(g)$，所以 $Q_2TQ_1^*\in\pi(G)'$，它等于
$P_2(Q_2TQ_1^*)P_1$，且 $Q_2^*(Q_2TQ_1^*)Q_1=T$。反之，每个 $Y=P_2XP_1$（$X\in\pi(G)'$）都属于
$\pi(G)'$，且 $Y=Q_2(Q_2^*YQ_1)Q_1^*$。$\square$

**定理 12.4（等变映射是分块矩阵）。** 设 $\rho_1,\dots,\rho_C$ 是在 $\pi_1\oplus\pi_2$ 中出现的两两
不等价的不可约表示，$k_c=\dim\rho_c$，$a_c=m_{\rho_c}(\pi_1)$，$b_c=m_{\rho_c}(\pi_2)$。存在酉算子
$V_1:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{a_c}\to\mathbb C^{n_1}$ 与
$V_2:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{b_c}\to\mathbb C^{n_2}$，满足
$V_i^*\pi_i(g)V_i=\bigoplus_c\rho_c(g)\otimes1$，使得
$$
\operatorname{Hom}_G(\pi_1,\pi_2)=\Big\{V_2\Big(\bigoplus_{c=1}^C1_{k_c}\otimes y_c\Big)V_1^*:\ y_c\in M_{b_c\times a_c}(\mathbb C)\Big\},
$$
并且 $(y_c)_c\mapsto T$ 是单射。因此 $\operatorname{Hom}_G(\pi_1,\pi_2)\cong\bigoplus_cM_{b_c\times a_c}(\mathbb C)$，
$$
\dim\operatorname{Hom}_G(\pi_1,\pi_2)=\sum_\rho m_\rho(\pi_1)\,m_\rho(\pi_2),
$$
其中求和遍及等价意义下的所有不可约表示；当 $\pi_1=\pi_2$ 时，这就是
$\pi(G)'\cong\bigoplus_cM_{m_c}$（命题 2.7），并且 $\operatorname{Hom}_G(\pi_1,\pi_2)=0$ 当且仅当没有
不可约表示同时出现在两者之中。

*证明。* 对 $\pi=\pi_1\oplus\pi_2$ 应用定理 2.6 与命题 2.7：存在酉算子 $V$，使得
$V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$ 且 $V^*\pi(G)'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$。由于
$X\mapsto V^*XV$ 是 $*$-同构且 $P_1\in\pi(G)'$（引理 12.3），有 $V^*P_1V=\bigoplus_c1\otimes p_c$，其中
$p_c$ 为投影。选取酉矩阵 $u_c$ 使 $u_c^*p_cu_c=1_{a_c}\oplus0_{b_c}$（定理 1.3），这里暂时令
$a_c=\operatorname{rank}p_c$，$b_c=m_c-a_c$，并把 $V$ 换成 $V\bigoplus_c(1\otimes u_c)$；这不改变两个
标准形，因为 $1\otimes u_c$ 与 $\rho_c(g)\otimes1$ 可交换。设 $\iota_c:\mathbb C^{a_c}\to\mathbb C^{m_c}$ 与
$\kappa_c:\mathbb C^{b_c}\to\mathbb C^{m_c}$ 分别是前 $a_c$ 个坐标与后 $b_c$ 个坐标的嵌入，并令
$$
V_1=Q_1^*V\bigoplus_c(1\otimes\iota_c),\qquad V_2=Q_2^*V\bigoplus_c(1\otimes\kappa_c).
$$
$W_1=V\bigoplus_c(1\otimes\iota_c)$ 是值域为 $\operatorname{range}P_1=\operatorname{range}Q_1$ 的等距算子，所以
$V_1^*V_1=W_1^*P_1W_1=1$，又因为 $\sum_ck_ca_c=\operatorname{rank}P_1=n_1$，$V_1$ 是酉的。此外
$\pi_1(g)V_1=Q_1^*\pi(g)W_1=Q_1^*V\big(\bigoplus_c\rho_c(g)\otimes1_{m_c}\big)\bigoplus_c(1\otimes\iota_c)
=V_1\bigoplus_c\rho_c(g)\otimes1_{a_c}$。所以 $\pi_1\cong\bigoplus_c\rho_c\otimes1_{a_c}$，其中各 $\rho_c$
为不等价的不可约表示，由重数的唯一性（2.4 节）得 $a_c=m_{\rho_c}(\pi_1)$；对 $V_2$ 与 $b_c$ 同理。
最后，对 $X=V\big(\bigoplus_c1\otimes x_c\big)V^*\in\pi(G)'$，有
$P_2XP_1=V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*$，其中 $y_c=\kappa_c^*x_c\iota_c$ 是 $x_c$ 左下角的
$b_c\times a_c$ 子块，它是任意矩阵。由引理 12.3，
$\operatorname{Hom}_G(\pi_1,\pi_2)=Q_2^*P_2\,\pi(G)'P_1Q_1$ 由矩阵
$Q_2^*V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*Q_1=V_2\big(\bigoplus_c1\otimes y_c\big)V_1^*$ 构成，并且由于
$V_1,V_2$ 是酉的，$(y_c)\mapsto T$ 是单射。不在 $\pi_1\oplus\pi_2$ 中出现的不可约表示 $\rho$ 满足
$m_\rho(\pi_1)m_\rho(\pi_2)=0$。$\square$

对 $n=4$ 时的推论 12.2(b)：在 $\pi_1\oplus\pi_2$ 的四个扇区（平凡表示、维数为 $2$ 的表示、两个维数为
$3$ 的表示）中，顶点特征 $\pi_1$ 有 $a=(1,0,0,1)$，边特征 $\pi_2$ 有 $b=(2,1,1,3)$，所以
$\dim\operatorname{Hom}_G=1\cdot2+1\cdot3=5=B_3$。

> **物理意义。** 定理 12.4 是**选择定则**的线性代数核心，也是 Wigner–Eckart 定理的标量情形：与
> 对称性可交换的算子只能连接对称类型相同的态，并且在同一类型内部它以 $1_{k_c}\otimes y_c$ 的方式作用，
> 即通过一个对多重态内部的位置不敏感的“约化矩阵元” $y_c$ 作用。（本身按非平凡表示变换的算子，例如
> 矢量算子，还需要 Clebsch–Gordan 系数。）对称的哈密顿量是 $\pi_1=\pi_2$ 的情形；它的简并度 $k_c$
> 就是各不可约表示的维数。

## 12.2 用库构造等变层

一旦知道了定理 2.6 的分解，定理 12.4 的证明就是一个算法；`dc.from_group` 用 4.3 节的随机化方法
计算这个分解。当 $\pi_1=\pi_2$ 时，它就是定理 2.6 本身：权重是块 $w_c\in M_{m_c}$，即抽象交换子
`D.commutant_algebra()` $=\bigoplus_cM_{m_c}$ 的一个算子，而层是
$x\mapsto V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*x$，由 `D.from_commutant` 计算。下面的模块是本章的
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

当 $\pi_1\ne\pi_2$ 时，逐行照搬证明：分解 $\pi_1\oplus\pi_2$，用 `to_commutant` 计算 $P_1$ 的块 $p_c$，
将它们对角化，并保留左下角。$p_c$ 的迹是 $a_c=m_c(\pi_1)$。

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
corners = node_to_edge.D.basis('commutant')[:, n:, :n].reshape(-1, 16 * n)   # Lemma 12.3: Q_2^* pi(G)' Q_1
assert torch.linalg.matrix_rank(corners).item() == bell[3]
```

使层等变的第二种方法是在群上对它取平均。**twirl** 定义为
$$
\mathcal P(W)=\frac1{|G|}\sum_{g\in G}\pi_2(g)\,W\,\pi_1(g)^*,\qquad W\in M_{n_2\times n_1}(\mathbb C).
$$

**定理 12.5（twirl 是正交投影）。** $\mathcal P$ 是 $M_{n_2\times n_1}(\mathbb C)$ 到
$\operatorname{Hom}_G(\pi_1,\pi_2)$ 上关于 Hilbert–Schmidt 内积的正交投影，并且
$$
\dim\operatorname{Hom}_G(\pi_1,\pi_2)=\frac1{|G|}\sum_{g\in G}\overline{\chi_{\pi_1}(g)}\,\chi_{\pi_2}(g).
$$

*证明。* 对 $h\in G$，$\pi_2(h)\mathcal P(W)\pi_1(h)^*=\frac1{|G|}\sum_g\pi_2(hg)W\pi_1(hg)^*=\mathcal P(W)$，因为
$g\mapsto hg$ 是 $G$ 的双射；所以 $\mathcal P(W)$ 是等变的。若 $W$ 等变，则每一项都等于 $W$，所以
$\mathcal P(W)=W$。因此 $\mathcal P^2=\mathcal P$，其值域为 $\operatorname{Hom}_G(\pi_1,\pi_2)$。它是自伴的：
$\langle S,\mathcal P(W)\rangle=\frac1{|G|}\sum_g\operatorname{Tr}\big((\pi_2(g)^*S\pi_1(g))^*W\big)$
$=\langle\mathcal P(S),W\rangle$，这里作了代换 $g\mapsto g^{-1}$ 并用到 $\pi_i(g^{-1})=\pi_i(g)^*$。自伴
幂等算子就是到其值域上的正交投影。它的迹等于它的秩。$M_{n_2\times n_1}$ 上的线性映射
$W\mapsto AWB$ 的迹为 $\operatorname{Tr}A\operatorname{Tr}B$（矩阵单位 $E_{ij}$ 在 $AE_{ij}B$ 中的系数为
$A_{ii}B_{jj}$）；取 $A=\pi_2(g)$，$B=\pi_1(g)^*$ 即得该公式。$\square$

结合定理 12.4，$\sum_\rho m_\rho(\pi_1)m_\rho(\pi_2)=|G|^{-1}\sum_g\overline{\chi_{\pi_1}(g)}\chi_{\pi_2}(g)$；
当 $\pi_1=\pi_2$ 时，这就是 4.2 节中 $\dim\pi(G)'$ 的特征标公式，`dc.commutant_dimension` 计算的正是它。
参数化 $W\mapsto\mathcal P(W)$ 能取到每一个等变映射，$\mathcal P(W)$ 是离 $W$ 最近的等变映射，并且由于
$\mathcal P$ 是自伴的，关于 $W$ 的梯度就是关于 $T=\mathcal P(W)$ 的梯度在 $\mathcal P$ 下的像。但有
$n_1n_2-\dim\operatorname{Hom}_G$ 个参数浪费在核上（对边特征为 $256-15$ 个），并且对稠密的 $\pi(g)$，
一次求值的代价为 $O(|G|n_1n_2(n_1+n_2))$，而对 $S_n$ 有 $|G|=n!$。库从分解出发计算同一个投影，即
`D.conditional_expectation(W, 'commutant')`，只需四次 $n\times n$ 矩阵乘积。

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

上面的层都是线性的。网络让它们与非线性交替出现，而这些非线性也必须是等变的。对置换表示，逐坐标
作用的函数是等变的，因为 $(\pi(g)x)_{gy}=x_y$；对其他表示则不然，此时人们使用作用于范数的非线性，
或者作用于一组被 $G$ 置换的基的系数上的非线性（习题 3）。

> **物理意义。** twirl 是到交换子上的条件期望（第 2 章习题 3）：把一个耦合在所有对称操作上取平均，
> 恰好保留它的对称部分。通过对对称群求和来“对称化”哈密顿量或相互作用，就是同一个投影。

## 12.3 带约束的层与谱层

`EquivariantLinear` 的权重 $w=(w_c)_c$ 位于一个 $*$-代数中，所以第 1 章与第 11 章的每个构造都适用于
它们。对 $w\in\bigoplus_cM_{m_c}$ 记 $\Phi(w)=V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*$。

**命题 12.6.** (a) $\Phi$ 是从 $\bigoplus_cM_{m_c}$ 到 $\pi(G)'$ 上的 $*$-同构，并且 $\Phi(w)$ 是自伴的、
酉的或正的，当且仅当每个 $w_c$ 都是。(b) 对正规的 $w$ 以及
$\bigcup_c\operatorname{spec}w_c$ 上的函数 $f$，有 $f(\Phi(w))=\Phi\big((f(w_c))_c\big)$；特别地，对每个正规的
$X\in\pi(G)'$ 有 $f(X)\in\pi(G)'$。(c) 对任意固定的酉元 $u_0\in\bigoplus_cM_{m_c}$，等变酉算子恰为
$\Phi(u_0e^{s-s^*})$，其中 $s$ 取遍 $\bigoplus_cM_{m_c}$；特别地，它们恰为反自伴的 $S\in\pi(G)'$ 的指数 $e^S$。

*证明。* (a) 由于 $(1\otimes x)(1\otimes y)=1\otimes xy$ 且 $(1\otimes x)^*=1\otimes x^*$，$\Phi$ 是线性、
保单位、可乘且保持 $*$ 的；它是单射，并由定理 2.6 映满 $\pi(G)'$。因此 $\Phi(w)^*=\Phi(w)$ 当且仅当
$w^*=w$，且 $\Phi(w)^*\Phi(w)=1=\Phi(w)\Phi(w)^*$ 当且仅当 $w^*w=1=ww^*$。若所有 $w_c=t_c^*t_c$，则
$\Phi(w)=\Phi(t)^*\Phi(t)\ge0$（推论 1.7）。反之，若 $\Phi(w)\ge0$，则对 $v\in\mathbb C^{m_c}$ 以及扇区 $c$
中的向量 $\xi=V(e_1\otimes v)$，有 $0\le\langle\xi,\Phi(w)\xi\rangle=\langle v,w_cv\rangle$。
(b) 设 $w_c=u_c\Lambda_cu_c^*$（定理 1.3）。则 $\Phi(w)=U\big(\bigoplus_c1\otimes\Lambda_c\big)U^*$，其中酉算子
$U=V\bigoplus_c(1\otimes u_c)$；这是一个对角化，$w_c$ 的特征值各重复 $k_c$ 次；由定义与定理 1.5(a)，
$f(\Phi(w))=U\big(\bigoplus_c1\otimes f(\Lambda_c)\big)U^*=\Phi\big((f(w_c))_c\big)$。
(c) 由 (a)，等变酉算子是 $\Phi(u)$，其中 $u$ 是 $\bigoplus_cM_{m_c}$ 中的酉元，而由引理 11.1，这些
$u$ 恰为 $u_0e^{s-s^*}$。取 $u_0=1$，(b) 给出 $\Phi(e^{s-s^*})=e^{\Phi(s-s^*)}$，而 $\Phi(s-s^*)$ 取遍
$\pi(G)'$ 的反自伴元素。$\square$

所以正交（酉）等变层就是在 `D.commutant_algebra()` 上的 `optimize.UnitaryParam` 再接上 $\Phi$；由 (c)，
这个参数化能取到每个等变酉算子，并且梯度方法的每一步迭代都精确地是酉的、精确地是等变的。类似地，
`PositiveParam` 给出正的等变层，其归一化如第 11 章（谱在 $[\ell,1]$ 中且最大值为 $1$；由 (a)，一般的
正层是 $\Phi(t^*t)$）；而自伴等变生成元 $L$ 的谱层 $f(L)$，例如扩散 $e^{-tL}$，所需特征分解的代价为
$\sum_cm_c^3$（这里为 $37$），而不是 $n^3=4096$，其梯度由 Daleckii–Krein 公式给出（`apply_function`）。
下面的代码用 $32$ 对输入–输出数据拟合一个等变酉算子，并将 (b) 以及关于 $t$ 的导数
$\frac d{dt}e^{-tL}=-Le^{-tL}$ 与稠密矩阵的结果对照检验。

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

> **物理意义。** 对称的哈密顿量 $H\in\pi(G)'$ 生成对称的时间演化
> $e^{-itH}\in\pi(G)'$（第 8 章），它从不混合扇区，只作用在重数空间上。等变酉层就是这样一种演化，
> 只不过它是学习得到的，而不是给定的。

## 12.4 不变量

若对所有 $g$ 有 $\pi(g)x=x$，则称向量 $x$ 是**不变的**；若 $q(\pi(g)x)=q(x)$，则称 $\mathbb C^n$ 上的
函数 $q$ 是不变的。不变特征就是不变网络（例如图分类器）所能输出的东西。

**命题 12.7.** (a) 不变向量构成正交投影 $P_0=|G|^{-1}\sum_g\pi(g)$ 的值域；其维数为
$m_{\mathrm{triv}}(\pi)=|G|^{-1}\sum_g\chi_\pi(g)$，并且若平凡表示作为 $\rho_{c_0}$ 出现，则 $P_0$ 是
$\pi(G)''$ 的极小中心投影 $P_{c_0}$。
(b) Hermite 型 $q_A(x)=\langle x,Ax\rangle$（$A=A^*$）是不变的，当且仅当 $A\in\pi(G)'$；不变 Hermite 型
构成维数为 $\sum_cm_c^2$ 的实向量空间。
(c) 关于坐标 $x_1,\dots,x_n$ 的 $2$ 次齐次不变多项式构成维数为
$$
\frac1{|G|}\sum_{g\in G}\frac{\chi_\pi(g)^2+\chi_\pi(g^2)}2 .
$$
的空间。对置换表示，它等于 $G$ 在由 $\Omega$ 的两个元素组成的多重集 $\{x,y\}$ 上的轨道个数。

*证明。* (a) 不变向量就是 $\operatorname{Hom}_G(\mathrm{triv},\pi)$，其中 $\mathrm{triv}$ 是
$\mathbb C=M_{1\times1}$ 上的平凡表示；对它而言 twirl 为 $\mathcal P(x)=P_0x$，所以定理 12.5 给出投影与
维数，定理 12.4 给出 $m_{\mathrm{triv}}(\pi)$。在定理 2.6 的基下，把不变向量写成 $\xi=\sum_c\xi_c$，其中
$\xi_c=\sum_jv_{c,j}\otimes e_j\in\mathbb C^{k_c}\otimes\mathbb C^{m_c}$。由于 $\pi(g)$ 保持各扇区，对所有 $g$ 有
$\rho_c(g)v_{c,j}=v_{c,j}$。非零的 $v_{c,j}$ 张成一条承载平凡表示的不变直线，于是 $\rho_c$ 作为不可约
表示是平凡的，且 $c=c_0$。因此不变向量就是扇区 $c_0$ 中的向量，$\pi(g)$ 在那里作用为 $1$：它们构成
$P_{c_0}$ 的值域，因此 $P_{c_0}$ 就是 $P_0$。
(b) $q_A(\pi(g)x)=\langle x,\pi(g)^*A\pi(g)x\rangle$，而对所有 $x$ 满足 $\langle x,Dx\rangle=0$ 的矩阵 $D$ 必为零
（极化，见推论 1.4 的证明）。所以 $q_A$ 不变当且仅当对所有 $g$ 有 $\pi(g)^*A\pi(g)=A$，即 $A\in\pi(G)'$。
由命题 12.6(a)，$\pi(G)'$ 的自伴元素是各 $w_c$ 为 Hermite 矩阵的 $\Phi(w)$，它们构成维数为
$\sum_cm_c^2$ 的实空间。
(c) 二次多项式 $p(x)=x^TBx$ 由其对称矩阵 $B$ 决定。令 $(g\cdot p)(x)=p(\pi(g)^{-1}x)$；则
$g\cdot(h\cdot p)=(gh)\cdot p$，且 $g\cdot p$ 的矩阵为 $\sigma(g)(B)=A^TBA$，$A=\pi(g)^{-1}$。平均
$M=|G|^{-1}\sum_g\sigma(g)$ 满足 $\sigma(h)M=M$，并且对不变的 $B$ 有 $MB=B$，所以（与定理 12.5 的证明
一样）它是幂等的，其值域为不变多项式，而不变多项式空间的维数为
$\operatorname{Tr}M=|G|^{-1}\sum_g\operatorname{Tr}\sigma(g)$。$A$ 是酉的，因此 $A=R\Lambda R^{-1}$，其中
$\Lambda=\operatorname{diag}(\lambda)$（定理 1.3）。对称矩阵
$B_{ij}=(R^{-1})^T(e_ie_j^T+e_je_i^T)R^{-1}$（$i\le j$）构成一组基，并且
$A^TB_{ij}A=(R^{-1})^T\Lambda(e_ie_j^T+e_je_i^T)\Lambda R^{-1}=\lambda_i\lambda_jB_{ij}$。所以
$\operatorname{Tr}\sigma(g)=\sum_{i\le j}\lambda_i\lambda_j=\frac12\big((\operatorname{Tr}A)^2+\operatorname{Tr}A^2\big)
=\frac12\big(\chi_\pi(g^{-1})^2+\chi_\pi(g^{-2})\big)$，再把对 $g$ 的求和换成对 $g^{-1}$ 的求和即得该公式。
对置换表示，$G$ 置换单项式 $x_ax_b$，所以多项式不变当且仅当其系数在单项式的轨道上为常数。$\square$

对于图，不变边特征由对角线及其补集的示性向量张成（$B_2=2$），不变的线性读出是对角元之和与非对角元
之和。下面的代码从块 $\rho_c(g)$（`to_algebra`）中找出平凡扇区，并检验 (a)、(b) 以及 (c) 的两种计数。

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

在 $\mathbb C^4$ 上，它们是 $\sum_ix_i^2$ 与 $\sum_{i<j}x_ix_j$。所有次数的不变多项式由 Molien 公式
$\sum_d\dim\mathbb C[x]^G_dt^d=|G|^{-1}\sum_g\det(1-t\pi(g))^{-1}$ 计数，其 $t^2$ 的系数就是 (c)（此处不加
证明地引用；B. Sturmfels, *Algorithms in Invariant Theory*, Springer 2008，定理 2.2.1）。

> **物理意义。** 在 Landau 相变理论中，自由能是实序参量 $x$ 的不变多项式，$G$ 通过实正交矩阵
> 作用在 $x$ 上。它的二次部分是 $x^TBx$，其中 $B$ 属于交换子（因为 $\pi(g)^T=\pi(g)^{-1}$），所以由
> 命题 12.6(a)，它在各扇区上是分块对角的：对称相在 $B$ 的某个特征值首先改变符号的那个扇区中变得
> 不稳定，并且（在一般情形下）该扇区的不可约表示决定了有序相的对称性。

## 12.5 可与之比较的工具

针对特定群的等变层已有成熟的库提供。**escnn**（Cesa、Lang 与 Weiler 2022，e2cnn 的后继者）由不可约
表示和核的调和基为 $O(2)$ 与 $O(3)$ 的子群构造可操控卷积；**e3nn**（Geiger 与 Smidt 2022）由球谐函数
与 Clebsch–Gordan 张量积构造 $E(3)$-等变网络。**EMLP**（Finzi、Welling 与 Wilson 2021）处理任意矩阵群：
它对离散生成元施加 $\pi_2(g)T=T\pi_1(g)$，对 Lie 代数生成元施加 $d\pi_2(X)T=Td\pi_1(X)$，并把解空间的
一组基作为数值零空间算出——这正是 `dc.commutant_basis` 所解的线性问题。群卷积（Cohen 与 Welling 2016）
是正则表示的情形，而 Kondor 与 Trivedi（2018）证明了，对于紧群齐性空间上的特征映射，等变线性层就是
广义卷积。对置换而言，DeepSets（Zaheer et al. 2017）与 Maron et al.（2019）的基就是推论 12.2。

本库不提供网格上的卷积、连续空间上的可操控核、$\mathbb R$ 上的分解（实型、复型与四元数型；习题 2
处理实矩阵），也不适用于大的 $n$：$V$ 是稠密的 $n\times n$ 矩阵，并且对由全部元素给出的群，分解的
代价为每个随机样本 $O(|G|n^3)$（Lie 代数生成元通过 `dc.from_generators` 处理；4.3 节）。它所增加的
是解空间的**结构**：交换子的一组基只是一列矩阵，而分解把交换子与 $*$-代数 $\bigoplus_cM_{m_c}$
等同起来。于是乘积、伴随、酉性、正性和函数演算都在小矩阵上逐块作用，并带有梯度（命题 12.6），而
`to_commutant`/`from_commutant` 在两种图景之间转换。

## 习题

1. （**带通道的 DeepSets。**）$S_n$（$n\ge2$）通过 $P_\sigma\otimes1_d$ 作用在 $\mathbb C^n\otimes\mathbb C^d$ 上。
   证明从 $\mathbb C^n\otimes\mathbb C^{d_1}$ 到 $\mathbb C^n\otimes\mathbb C^{d_2}$ 的等变映射恰为
   $1\otimes A+J\otimes B$，其中 $A,B\in M_{d_2\times d_1}$，并用 `EquivariantMap` 加以验证。
2. （**实的层。**）设 $\pi_1,\pi_2$ 的矩阵是实的（例如置换表示）。证明实的等变矩阵构成维数为
   $\dim_{\mathbb C}\operatorname{Hom}_G(\pi_1,\pi_2)$ 的实向量空间，并且 $T\mapsto\operatorname{Re}T$ 把
   $\operatorname{Hom}_G(\pi_1,\pi_2)$ 映满这个空间。
3. （**非线性。**）(a) 证明对每个坐标应用函数 $f:\mathbb C\to\mathbb C$ 与每个置换表示都可交换。
   (b) 证明对 $f(t)=\max(t,0)$，它与 $\mathbb R^2$ 上转角为 $\pi/2$ 的旋转不可交换。(c) 证明
   $x\mapsto\varphi(\|x\|)x$ 与每个酉表示都可交换。
4. （**多少个样本能确定一个层？**）设 $T=\Phi(w)\in\pi(G)'$，$x_1,\dots,x_N\in\mathbb C^n$。把 $V^*x_i$
   在扇区 $c$ 中的部分写成 $\sum_{a,j}(X_{i,c})_{aj}e_a\otimes e_j$，其中 $X_{i,c}\in M_{k_c\times m_c}$，并设
   $S_c\in M_{Nk_c\times m_c}$ 是把矩阵 $X_{1,c},\dots,X_{N,c}$ 上下堆叠所得的矩阵。证明满足
   $Tx_1=\dots=Tx_N=0$ 的等变 $T$ 构成维数为 $\sum_cm_c(m_c-\operatorname{rank}S_c)$ 的空间；对高斯的 $x_i$，
   它几乎必然等于 $\sum_cm_c\max(0,m_c-Nk_c)$。对边表示，两个样本就能确定一个等变层，而任意的
   $16\times16$ 矩阵需要十六个样本。

<details>
<summary>解答</summary>

1. $P_\sigma J=J=JP_\sigma$，所以这些映射是等变的。它们线性无关：若
   $1\otimes A+J\otimes B=0$，则非对角块给出 $B=0$，进而 $A=0$。由推论 12.2(a) 与
   定理 2.6，$P_\sigma\otimes1_d\cong(\mathrm{triv}\otimes1_d)\oplus(\mathrm{std}\otimes1_d)$，所以两个不可约
   表示的重数都是 $d$，而定理 12.4 给出 $\dim\operatorname{Hom}_G=2d_1d_2$，这正是由形如
   $1\otimes A+J\otimes B$ 的映射构成的空间的维数。

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

2. 若 $T$ 等变，则 $\bar T$ 也等变，因为 $\overline{\pi_i(g)}=\pi_i(g)$。因此
   $\operatorname{Re}T=(T+\bar T)/2$ 与 $\operatorname{Im}T=(T-\bar T)/2i$ 是实的且等变，所以映射
   $T\mapsto\operatorname{Re}T$ 映满实的等变矩阵（在它们上面它是恒等映射），并且
   $\operatorname{Hom}_G(\pi_1,\pi_2)=H_{\mathbb R}\oplus iH_{\mathbb R}$，其中 $H_{\mathbb R}$ 是实的等变
   矩阵全体。于是 $H_{\mathbb R}$ 的一组实基就是 $\operatorname{Hom}_G$ 的一组复基：它张成整个空间，而对
   实的 $T_j$，复线性关系 $\sum_j(\alpha_j+i\beta_j)T_j=0$ 可以拆成实部与虚部。在引理 12.1 中，轨道的
   示性矩阵就是这样一组实基。（即使数据是实的，定理 12.4 的块 $y_c$ 也是复的：一般而言 $V$ 是复的，
   而用 $\operatorname{Re}T$ 参数化所用的实参数是所需的两倍。）

3. (a) $(\pi(g)f(x))_{gy}=f(x)_y=f(x_y)=f((\pi(g)x)_{gy})$。(b) 对 $R(x_1,x_2)=(-x_2,x_1)$ 与 $x=(1,-1)$，
   $f(Rx)=f(1,1)=(1,1)$，但 $Rf(x)=R(1,0)=(0,1)$。(c) $\|\pi(g)x\|=\|x\|$，所以
   $\varphi(\|\pi(g)x\|)\pi(g)x=\pi(g)\varphi(\|x\|)x$。

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

4. $(1\otimes w_c)(e_a\otimes e_j)=\sum_le_a\otimes(w_c)_{lj}e_l$，所以 $V^*Tx_i$ 在扇区 $c$ 中的系数矩阵为
   $X_{i,c}w_c^T$。因此对所有 $i$ 有 $Tx_i=0$，当且仅当对所有 $c$ 有 $S_cw_c^T=0$，即当且仅当 $w_c^T$ 的
   每一列都属于 $\ker S_c$，这是一个维数为 $m_c-\operatorname{rank}S_c$ 的空间：具有这一性质的 $w_c$ 构成
   维数为 $m_c(m_c-\operatorname{rank}S_c)$ 的空间，而 $\Phi$ 是单射。对高斯的 $x_i$，$V^*x_i$ 是高斯的
   （$V$ 是酉的），所以 $S_c$ 的元素是独立的高斯变量，并且几乎必然有秩 $\min(Nk_c,m_c)$：对于该尺寸的
   左上角方子矩阵，$|\det|^2$ 是各元素实部与虚部的实多项式，且不恒为零（当该子矩阵为单位矩阵时它
   等于 $1$），而这样的多项式只在 Lebesgue 零测集上为零（在 4.3 节中证明）。对边表示，
   $(k,m)=(1,2),(2,1),(3,1),(3,3)$，只有扇区 $(1,2)$ 有贡献：$N=1$ 时为 $2\cdot(2-1)=2$，$N=2$ 时为 $0$。

```python
basis = D_edge.basis('commutant')                                # (15, 16, 16), Hilbert-Schmidt orthonormal
for N, lost in ((1, 2), (2, 0)):
    Xs = torch.randn(16, N, dtype=dt)
    M = torch.stack([(B @ Xs).reshape(-1) for B in basis], dim=1)               # w -> (T x_1, ..., T x_N)
    assert torch.linalg.matrix_rank(M).item() == 15 - lost
```

</details>

**下一步。** 这是本课程的最后一章。[第 4 章](04_symmetries.md)发展了这里用到的表示论：Schur 引理、
特征标，以及 `dc.from_group` 背后的随机化算法（4.1–4.3 节）。[第 9 章](09_spin_chains.md)把同样的扇区应用于
自旋链，而[第 5 章](05_random_matrices.md)中 $\prod_cU(m_c)$ 上的 Haar 测度给出了等变酉层的自然随机
初始化。数学方面参见 J.-P. Serre, *Linear Representations of Finite Groups* (Springer 1977) 与
W. Fulton and J. Harris, *Representation Theory* (Springer 1991)；等变学习方面参见 M. Bronstein, J. Bruna,
T. Cohen and P. Veličković, *Geometric Deep Learning* (2021) 与 M. Weiler et al., *Equivariant and Coordinate
Independent Convolutional Networks* (2023)；分解算法方面参见 `torch_vn_algebra/decompose.py` 中的参考文献。
