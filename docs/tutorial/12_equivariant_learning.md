**English** · [Español](es/12_equivariant_learning.md) · [Français](fr/12_equivariant_learning.md) · [Deutsch](de/12_equivariant_learning.md) · [中文](zh/12_equivariant_learning.md) · [日本語](ja/12_equivariant_learning.md) · [Русский](ru/12_equivariant_learning.md)

# 12. Equivariant machine learning

A layer of a neural network is **equivariant** if transforming the input by a symmetry and then applying the
layer gives the same result as applying the layer and then transforming the output: relabelling the vertices of a
graph relabels its features in the same way. For linear layers this is a commutation relation, so the equivariant
layers form a commutant, which by Chapter 2 is, up to a change of basis, the algebra $\bigoplus_cM_{m_c}(\mathbb C)$.
This chapter counts their parameters, builds and trains them with `torch_vn_algebra.decompose`, and uses the
algebra structure for unitary, positive and spectral layers that are equivariant by construction.

**You will learn** why the equivariant linear maps between two representations are block matrices of
sizes given by multiplicities; how to count them by orbits or by characters (DeepSets and the Bell numbers of
graph networks); why averaging over the group is the orthogonal projection onto the equivariant maps; how
functional calculus and unitary parametrisations act blockwise on the commutant; and how invariant vectors
and invariant quadratic forms are counted.

Throughout, $G$ is a finite group, and a **representation** of $G$ on $\mathbb C^n$ is a homomorphism
$\pi:G\to M_n(\mathbb C)$ with unitary values (no restriction: invertible $\pi(g)$ are unitary for the inner product
$\sum_g\langle\pi(g)x,\pi(g)y\rangle$, Section 4.1). Its **character** is $\chi_\pi(g)=\operatorname{Tr}\pi(g)$. Irreducible and
equivalent representations and the **multiplicity** $m_\rho(\pi)$ of an irreducible $\rho$ in $\pi$ ($0$ if $\rho$
does not occur) are as in Section 2.4, which proves that $m_\rho(\pi)$ does not depend on the decomposition and
proves Schur's lemma; [Chapter 4](04_symmetries.md) treats both at length (Sections 4.1–4.2). On rectangular matrices
$\langle S,T\rangle=\operatorname{Tr}(S^*T)$ is the **Hilbert–Schmidt** inner product. In the code a batch of vectors
is a tensor of shape `(batch, n)`, on which a matrix $T$ acts as `x @ T.T`.

## 12.1 Equivariant linear maps

**Definition.** Let $\pi_1,\pi_2$ be representations of $G$ on $\mathbb C^{n_1}$, $\mathbb C^{n_2}$. A matrix
$T\in M_{n_2\times n_1}(\mathbb C)$ is **equivariant** (an **intertwiner**) if $T\pi_1(g)=\pi_2(g)T$ for all
$g\in G$. They form a vector space $\operatorname{Hom}_G(\pi_1,\pi_2)$; for $\pi_1=\pi_2=\pi$ it is the
commutant, $\operatorname{Hom}_G(\pi,\pi)=\pi(G)'$.

The commonest representations in machine learning permute coordinates. Let $G$ act on a finite set $\Omega$;
the **permutation representation** on $\mathbb C^\Omega$ is $\pi(g)e_x=e_{gx}$. Examples: $S_n$ permuting the
$n$ points of a set ($\Omega=[n]=\{1,\dots,n\}$), and $S_n$ relabelling the vertices of a graph, acting on
$\Omega=[n]^2$, i.e. on the adjacency matrix by $X\mapsto P_\sigma XP_\sigma^T$.

**Lemma 12.1 (orbits).** Let $\pi_i$ be the permutation representations of $G$ on $\mathbb C^{\Omega_i}$. Then
$T\in\operatorname{Hom}_G(\pi_1,\pi_2)$ iff $T_{gy,gx}=T_{yx}$ for all $g$, $y\in\Omega_2$, $x\in\Omega_1$. The
indicator matrices of the orbits of $G$ on $\Omega_2\times\Omega_1$ are a basis, and
$\dim\operatorname{Hom}_G(\pi_1,\pi_2)$ is the number of these orbits.

*Proof.* $\pi_2(g)T\pi_1(g)^*e_x=\pi_2(g)Te_{g^{-1}x}=\sum_yT_{y,g^{-1}x}e_{gy}$, so the $(gy,x)$ entry of
$\pi_2(g)T\pi_1(g)^*$ is $T_{y,g^{-1}x}$. Hence $T$ is equivariant, i.e. $\pi_2(g)T\pi_1(g)^*=T$ for all $g$, iff
$T_{gy,x}=T_{y,g^{-1}x}$ for all $g,y,x$, which (replace $x$ by $gx$) is the stated condition: $T$ is constant
on orbits. The indicator matrices of the orbits have disjoint supports, so they are linearly independent, and
every function constant on orbits is a combination of them. $\square$

**Corollary 12.2.** (a) (*DeepSets.*) For $S_n$ on $\mathbb C^n$, $n\ge2$, the commutant is
$\operatorname{span}\{1,J\}$, $J$ the all-ones matrix, and $\mathbb C^n=\mathbb C\mathbf 1\oplus\mathbf 1^\perp$ is
a sum of two inequivalent irreducible representations, the trivial one and the **standard** one of dimension
$n-1$. (b) (*Graph networks.*) The equivariant maps from $\mathbb C^{[n]^p}$ to $\mathbb C^{[n]^q}$ (tensors
of order $p$ and $q$ on $n$ vertices) are counted by the partitions of a set of $p+q$ elements into at most $n$
blocks; for $n\ge p+q$ their number is the **Bell number** $B_{p+q}$ ($B_2=2$, $B_3=5$, $B_4=15$).

*Proof.* (b) Two tuples $x,y\in[n]^{p+q}$ lie in the same orbit iff they have the same **equality pattern**,
the partition of the positions $\{1,\dots,p+q\}$ by $i\sim j\Leftrightarrow x_i=x_j$. The pattern is clearly
invariant. Conversely, if $x,y$ have the same pattern, $x_i\mapsto y_i$ is a well-defined injection of
$\{x_i\}$ into $[n]$, which extends to a permutation $\sigma$ with $\sigma x=y$. Every partition with at most
$n$ blocks is a pattern (give the blocks distinct values). Lemma 12.1 concludes. (a) For $p=q=1$ there are two
orbits on $[n]^2$, the diagonal and its complement, with indicators $1$ and $J-1$. So
$\dim\pi(S_n)'=2=\sum_cm_c^2$ (Theorem 2.6, Proposition 2.7), which forces two sectors with $m_c=1$, i.e. two
inequivalent irreducible summands, and $\pi(S_n)'$ is spanned by their projections $P_1,P_2$. The projection onto the
invariant line $\mathbb C\mathbf 1$ lies in $\pi(S_n)'$ (Lemma 2.1), so it is $P_1$ or $P_2$: the line, which
carries the trivial representation, is one sector, and $\mathbf 1^\perp$ is the other. $\square$

The DeepSets layer $x\mapsto\lambda x+\gamma\,(\mathbf 1^Tx)\mathbf 1$ (Zaheer et al. 2017) is (a); the bases of
Maron et al. (2019) for graph networks are (b). The code counts the orbits by brute force, compares with the
commutant computed numerically, and decomposes both representations of $S_4$:

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

So $\mathbb C^{16}$ contains the trivial representation twice, one of dimension $2$ once, and two of dimension $3$
once and three times: $16=1\cdot2+2\cdot1+3\cdot1+3\cdot3$ and $15=2^2+1^2+1^2+3^2$. For two representations,
everything reduces to one commutant:

**Lemma 12.3 (direct sum).** Let $\pi=\pi_1\oplus\pi_2$ on $\mathbb C^{n_1}\oplus\mathbb C^{n_2}$, with the
isometric inclusions $Q_i$ and projections $P_i=Q_iQ_i^*$. Then $P_1,P_2\in\pi(G)'$,
$\operatorname{Hom}_G(\pi_1,\pi_2)=Q_2^*\,\pi(G)'\,Q_1$, and $T\mapsto Q_2TQ_1^*$ is a linear bijection of
$\operatorname{Hom}_G(\pi_1,\pi_2)$ onto $P_2\,\pi(G)'P_1$.

*Proof.* $\pi(g)Q_i=Q_i\pi_i(g)$ and $Q_i^*\pi(g)=\pi_i(g)Q_i^*$, because $\pi(g)$ is block diagonal; hence
$P_i\pi(g)=Q_i\pi_i(g)Q_i^*=\pi(g)P_i$. If $X\in\pi(G)'$, then
$Q_2^*XQ_1\pi_1(g)=Q_2^*X\pi(g)Q_1=Q_2^*\pi(g)XQ_1=\pi_2(g)Q_2^*XQ_1$. If $T$ is equivariant, then
$\pi(g)Q_2TQ_1^*=Q_2\pi_2(g)TQ_1^*=Q_2T\pi_1(g)Q_1^*=Q_2TQ_1^*\pi(g)$, so $Q_2TQ_1^*\in\pi(G)'$, it equals
$P_2(Q_2TQ_1^*)P_1$, and $Q_2^*(Q_2TQ_1^*)Q_1=T$. Conversely every $Y=P_2XP_1$ with $X\in\pi(G)'$ lies in
$\pi(G)'$ and $Y=Q_2(Q_2^*YQ_1)Q_1^*$. $\square$

**Theorem 12.4 (equivariant maps are block matrices).** Let $\rho_1,\dots,\rho_C$ be the pairwise inequivalent
irreducible representations occurring in $\pi_1\oplus\pi_2$, $k_c=\dim\rho_c$, $a_c=m_{\rho_c}(\pi_1)$,
$b_c=m_{\rho_c}(\pi_2)$. There are unitaries $V_1:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{a_c}\to\mathbb C^{n_1}$
and $V_2:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{b_c}\to\mathbb C^{n_2}$ with
$V_i^*\pi_i(g)V_i=\bigoplus_c\rho_c(g)\otimes1$, such that
$$
\operatorname{Hom}_G(\pi_1,\pi_2)=\Big\{V_2\Big(\bigoplus_{c=1}^C1_{k_c}\otimes y_c\Big)V_1^*:\ y_c\in M_{b_c\times a_c}(\mathbb C)\Big\},
$$
and $(y_c)_c\mapsto T$ is injective. Hence $\operatorname{Hom}_G(\pi_1,\pi_2)\cong\bigoplus_cM_{b_c\times a_c}(\mathbb C)$,
$$
\dim\operatorname{Hom}_G(\pi_1,\pi_2)=\sum_\rho m_\rho(\pi_1)\,m_\rho(\pi_2),
$$
the sum running over all irreducible representations up to equivalence; for $\pi_1=\pi_2$ this is
$\pi(G)'\cong\bigoplus_cM_{m_c}$ (Proposition 2.7), and $\operatorname{Hom}_G(\pi_1,\pi_2)=0$ iff no
irreducible representation occurs in both.

*Proof.* Apply Theorem 2.6 and Proposition 2.7 to $\pi=\pi_1\oplus\pi_2$: there is a unitary $V$ with
$V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$ and $V^*\pi(G)'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$. Since
$X\mapsto V^*XV$ is a $*$-isomorphism and $P_1\in\pi(G)'$ (Lemma 12.3), $V^*P_1V=\bigoplus_c1\otimes p_c$ with
projections $p_c$. Choose unitaries $u_c$ with $u_c^*p_cu_c=1_{a_c}\oplus0_{b_c}$ (Theorem 1.3), where now
$a_c=\operatorname{rank}p_c$, $b_c=m_c-a_c$, and replace $V$ by $V\bigoplus_c(1\otimes u_c)$; this changes neither
normal form, as $1\otimes u_c$ commutes with $\rho_c(g)\otimes1$. Let $\iota_c:\mathbb C^{a_c}\to\mathbb C^{m_c}$ and
$\kappa_c:\mathbb C^{b_c}\to\mathbb C^{m_c}$ be the inclusions of the first $a_c$ and the last $b_c$ coordinates, and
$$
V_1=Q_1^*V\bigoplus_c(1\otimes\iota_c),\qquad V_2=Q_2^*V\bigoplus_c(1\otimes\kappa_c).
$$
$W_1=V\bigoplus_c(1\otimes\iota_c)$ is an isometry with range $\operatorname{range}P_1=\operatorname{range}Q_1$, so
$V_1^*V_1=W_1^*P_1W_1=1$, and $V_1$ is unitary because $\sum_ck_ca_c=\operatorname{rank}P_1=n_1$. Moreover
$\pi_1(g)V_1=Q_1^*\pi(g)W_1=Q_1^*V\big(\bigoplus_c\rho_c(g)\otimes1_{m_c}\big)\bigoplus_c(1\otimes\iota_c)
=V_1\bigoplus_c\rho_c(g)\otimes1_{a_c}$. So $\pi_1\cong\bigoplus_c\rho_c\otimes1_{a_c}$ with inequivalent
irreducible $\rho_c$, and $a_c=m_{\rho_c}(\pi_1)$ by the uniqueness of multiplicities (Section 2.4); likewise for
$V_2$ and $b_c$. Finally, for $X=V\big(\bigoplus_c1\otimes x_c\big)V^*\in\pi(G)'$,
$P_2XP_1=V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*$, where $y_c=\kappa_c^*x_c\iota_c$ is the lower left
$b_c\times a_c$ corner of $x_c$, an arbitrary matrix. By Lemma 12.3,
$\operatorname{Hom}_G(\pi_1,\pi_2)=Q_2^*P_2\,\pi(G)'P_1Q_1$ consists of the matrices
$Q_2^*V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*Q_1=V_2\big(\bigoplus_c1\otimes y_c\big)V_1^*$, and
$(y_c)\mapsto T$ is injective since $V_1,V_2$ are unitary. An irreducible $\rho$ that does not occur in
$\pi_1\oplus\pi_2$ has $m_\rho(\pi_1)m_\rho(\pi_2)=0$. $\square$

For Corollary 12.2(b) with $n=4$: vertex features $\pi_1$ have $a=(1,0,0,1)$ in the four sectors of
$\pi_1\oplus\pi_2$ (trivial, dimension $2$, the two of dimension $3$), edge features $\pi_2$ have
$b=(2,1,1,3)$, so $\dim\operatorname{Hom}_G=1\cdot2+1\cdot3=5=B_3$.

> **Physical meaning.** Theorem 12.4 is the linear-algebra core of **selection rules** and the scalar case of
> the Wigner–Eckart theorem: an operator that commutes with the symmetry can only connect states of the same
> symmetry type, and within a type it acts as $1_{k_c}\otimes y_c$, i.e. by a "reduced matrix element" $y_c$ that
> is blind to the position inside the multiplet. (Operators that themselves transform by a non-trivial
> representation, such as vectors, need in addition the Clebsch–Gordan coefficients.) A symmetric Hamiltonian
> is the case $\pi_1=\pi_2$; its degeneracies $k_c$ are the dimensions of the irreducible representations.

## 12.2 An equivariant layer with the library

The proof of Theorem 12.4 is an algorithm, once the decomposition of Theorem 2.6 is known; `dc.from_group`
computes it by the randomised method of Section 4.3. For $\pi_1=\pi_2$ it is just Theorem 2.6: the weights are blocks
$w_c\in M_{m_c}$, an operator of the abstract commutant `D.commutant_algebra()` $=\bigoplus_cM_{m_c}$, and the
layer is $x\mapsto V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*x$, computed by `D.from_commutant`. The module
below is chapter code, not part of the library.

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

For $\pi_1\ne\pi_2$ the proof is followed line by line: decompose $\pi_1\oplus\pi_2$, compute the blocks $p_c$ of
$P_1$ with `to_commutant`, diagonalise them, and keep the lower left corners. The trace of $p_c$ is
$a_c=m_c(\pi_1)$.

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

The second way to make a layer equivariant is to average it over the group. The **twirl** is
$$
\mathcal P(W)=\frac1{|G|}\sum_{g\in G}\pi_2(g)\,W\,\pi_1(g)^*,\qquad W\in M_{n_2\times n_1}(\mathbb C).
$$

**Theorem 12.5 (the twirl is the orthogonal projection).** $\mathcal P$ is the Hilbert–Schmidt orthogonal
projection of $M_{n_2\times n_1}(\mathbb C)$ onto $\operatorname{Hom}_G(\pi_1,\pi_2)$, and
$$
\dim\operatorname{Hom}_G(\pi_1,\pi_2)=\frac1{|G|}\sum_{g\in G}\overline{\chi_{\pi_1}(g)}\,\chi_{\pi_2}(g).
$$

*Proof.* For $h\in G$, $\pi_2(h)\mathcal P(W)\pi_1(h)^*=\frac1{|G|}\sum_g\pi_2(hg)W\pi_1(hg)^*=\mathcal P(W)$, as
$g\mapsto hg$ is a bijection of $G$; so $\mathcal P(W)$ is equivariant. If $W$ is equivariant, every term equals
$W$, so $\mathcal P(W)=W$. Hence $\mathcal P^2=\mathcal P$ with range $\operatorname{Hom}_G(\pi_1,\pi_2)$. It is
self-adjoint: $\langle S,\mathcal P(W)\rangle=\frac1{|G|}\sum_g\operatorname{Tr}\big((\pi_2(g)^*S\pi_1(g))^*W\big)$
$=\langle\mathcal P(S),W\rangle$, where we substituted $g\mapsto g^{-1}$ and used $\pi_i(g^{-1})=\pi_i(g)^*$. A
self-adjoint idempotent is the orthogonal projection onto its range. Its trace is its rank. The linear map
$W\mapsto AWB$ on $M_{n_2\times n_1}$ has trace $\operatorname{Tr}A\operatorname{Tr}B$ (the coefficient of the
matrix unit $E_{ij}$ in $AE_{ij}B$ is $A_{ii}B_{jj}$); with $A=\pi_2(g)$, $B=\pi_1(g)^*$ this gives the formula. $\square$

With Theorem 12.4, $\sum_\rho m_\rho(\pi_1)m_\rho(\pi_2)=|G|^{-1}\sum_g\overline{\chi_{\pi_1}(g)}\chi_{\pi_2}(g)$;
for $\pi_1=\pi_2$ this is the character formula for $\dim\pi(G)'$ of Section 4.2, which `dc.commutant_dimension`
evaluates. The parametrisation $W\mapsto\mathcal P(W)$ reaches every
equivariant map, $\mathcal P(W)$ is the equivariant map closest to $W$, and as $\mathcal P$ is self-adjoint, the
gradient with respect to $W$ is $\mathcal P$ of the gradient with respect to $T=\mathcal P(W)$. But
$n_1n_2-\dim\operatorname{Hom}_G$ parameters are wasted on the kernel ($256-15$ for edge features), and an
evaluation costs $O(|G|n_1n_2(n_1+n_2))$ for dense $\pi(g)$, with $|G|=n!$ for $S_n$. The library computes the same
projection from the decomposition, `D.conditional_expectation(W, 'commutant')`, with four $n\times n$ products.

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

The layers above are linear. A network alternates them with nonlinearities, and these must be equivariant too.
For a permutation representation, a function applied to each coordinate is equivariant, since
$(\pi(g)x)_{gy}=x_y$; for other representations it is not, and one uses nonlinearities that act on norms or on
the coefficients of a basis that $G$ permutes (Exercise 3).

> **Physical meaning.** The twirl is the conditional expectation onto the commutant (Exercise 3 of Chapter 2):
> averaging a coupling over all symmetry operations keeps exactly its symmetric part. "Symmetrising" a
> Hamiltonian or an interaction by summing over the symmetry group is the same projection.

## 12.3 Constrained and spectral layers

The weights $w=(w_c)_c$ of `EquivariantLinear` live in a $*$-algebra, so every construction of Chapters 1 and 11
applies to them. Write $\Phi(w)=V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*$ for $w\in\bigoplus_cM_{m_c}$.

**Proposition 12.6.** (a) $\Phi$ is a $*$-isomorphism of $\bigoplus_cM_{m_c}$ onto $\pi(G)'$, and $\Phi(w)$ is
self-adjoint, unitary or positive iff every $w_c$ is. (b) For normal $w$ and a function $f$ on
$\bigcup_c\operatorname{spec}w_c$, $f(\Phi(w))=\Phi\big((f(w_c))_c\big)$; in particular $f(X)\in\pi(G)'$ for every
normal $X\in\pi(G)'$. (c) The equivariant unitaries are the $\Phi(u_0e^{s-s^*})$, $s$ ranging over
$\bigoplus_cM_{m_c}$, for any fixed unitary $u_0\in\bigoplus_cM_{m_c}$; in particular they are the exponentials $e^S$
of the skew-adjoint $S\in\pi(G)'$.

*Proof.* (a) $\Phi$ is linear, unital, multiplicative and $*$-preserving because $(1\otimes x)(1\otimes y)=1\otimes xy$
and $(1\otimes x)^*=1\otimes x^*$; it is injective, and onto $\pi(G)'$ by Theorem 2.6. Hence $\Phi(w)^*=\Phi(w)$ iff
$w^*=w$, and $\Phi(w)^*\Phi(w)=1=\Phi(w)\Phi(w)^*$ iff $w^*w=1=ww^*$. If all $w_c=t_c^*t_c$, then
$\Phi(w)=\Phi(t)^*\Phi(t)\ge0$ (Corollary 1.7). Conversely, if $\Phi(w)\ge0$, then for $v\in\mathbb C^{m_c}$ and
the vector $\xi=V(e_1\otimes v)$ in sector $c$, $0\le\langle\xi,\Phi(w)\xi\rangle=\langle v,w_cv\rangle$.
(b) Let $w_c=u_c\Lambda_cu_c^*$ (Theorem 1.3). Then $\Phi(w)=U\big(\bigoplus_c1\otimes\Lambda_c\big)U^*$ with the
unitary $U=V\bigoplus_c(1\otimes u_c)$, a diagonalisation with the eigenvalues of $w_c$ repeated $k_c$ times; by
definition and Theorem 1.5(a), $f(\Phi(w))=U\big(\bigoplus_c1\otimes f(\Lambda_c)\big)U^*=\Phi\big((f(w_c))_c\big)$.
(c) By (a), the equivariant unitaries are $\Phi(u)$ with $u$ unitary in $\bigoplus_cM_{m_c}$, and these are the
$u_0e^{s-s^*}$ by Lemma 11.1. With $u_0=1$, (b) gives $\Phi(e^{s-s^*})=e^{\Phi(s-s^*)}$, and $\Phi(s-s^*)$ runs over
the skew-adjoint elements of $\pi(G)'$. $\square$

So an orthogonal (unitary) equivariant layer is `optimize.UnitaryParam` on `D.commutant_algebra()` followed by
$\Phi$; by (c) this parametrisation reaches every equivariant unitary, and every iterate of a gradient method is
exactly unitary and exactly equivariant. Likewise `PositiveParam` gives positive equivariant layers, normalised as
in Chapter 11 (spectrum in $[\ell,1]$ with maximum $1$; by (a), a general positive layer is $\Phi(t^*t)$), and a
spectral layer $f(L)$ of a self-adjoint equivariant generator $L$, such as the diffusion $e^{-tL}$, needs
eigendecompositions of cost $\sum_cm_c^3$ ($37$ here) instead of $n^3=4096$, with gradients by the
Daleckii–Krein formula (`apply_function`). The code fits an equivariant unitary to $32$ input–output pairs, and checks (b)
and its gradient with respect to $t$, $\frac d{dt}e^{-tL}=-Le^{-tL}$, against dense matrices.

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

> **Physical meaning.** A symmetric Hamiltonian $H\in\pi(G)'$ generates a symmetric time evolution
> $e^{-itH}\in\pi(G)'$ (Chapter 8), which never mixes sectors and acts on the multiplicity spaces only. An
> equivariant unitary layer is such an evolution, learned instead of given.

## 12.4 Invariants

A vector $x$ is **invariant** if $\pi(g)x=x$ for all $g$; a function $q$ on $\mathbb C^n$ is invariant if
$q(\pi(g)x)=q(x)$. Invariant features are what an invariant network (a graph classifier, say) can output.

**Proposition 12.7.** (a) The invariant vectors are the range of the orthogonal projection
$P_0=|G|^{-1}\sum_g\pi(g)$; their dimension is $m_{\mathrm{triv}}(\pi)=|G|^{-1}\sum_g\chi_\pi(g)$, and if the trivial
representation occurs, as $\rho_{c_0}$, then $P_0$ is the minimal central projection $P_{c_0}$ of $\pi(G)''$.
(b) A Hermitian form $q_A(x)=\langle x,Ax\rangle$, $A=A^*$, is invariant iff $A\in\pi(G)'$; the invariant Hermitian
forms are a real vector space of dimension $\sum_cm_c^2$.
(c) The homogeneous invariant polynomials of degree $2$ in the coordinates $x_1,\dots,x_n$ form a space of dimension
$$
\frac1{|G|}\sum_{g\in G}\frac{\chi_\pi(g)^2+\chi_\pi(g^2)}2 .
$$
For a permutation representation it is the number of orbits of $G$ on the multisets $\{x,y\}$ of two elements of $\Omega$.

*Proof.* (a) Invariant vectors are $\operatorname{Hom}_G(\mathrm{triv},\pi)$, with $\mathrm{triv}$ the trivial
representation on $\mathbb C=M_{1\times1}$; for it the twirl is $\mathcal P(x)=P_0x$, so Theorem 12.5 gives the
projection and the dimension, and Theorem 12.4 gives $m_{\mathrm{triv}}(\pi)$. In the basis of Theorem 2.6, write an
invariant vector as $\xi=\sum_c\xi_c$ with $\xi_c=\sum_jv_{c,j}\otimes e_j\in\mathbb C^{k_c}\otimes\mathbb C^{m_c}$.
Since $\pi(g)$ preserves the sectors, $\rho_c(g)v_{c,j}=v_{c,j}$ for all $g$. A non-zero $v_{c,j}$ spans an invariant
line carrying the trivial representation, so then $\rho_c$, being irreducible, is trivial and $c=c_0$. Hence the
invariant vectors are the vectors of sector $c_0$, where $\pi(g)$ acts as $1$: the range of $P_{c_0}$, which is
therefore $P_0$.
(b) $q_A(\pi(g)x)=\langle x,\pi(g)^*A\pi(g)x\rangle$, and a matrix $D$ with $\langle x,Dx\rangle=0$ for all $x$ vanishes
(polarisation, proof of Corollary 1.4). So $q_A$ is invariant iff $\pi(g)^*A\pi(g)=A$ for all $g$, i.e. $A\in\pi(G)'$.
By Proposition 12.6(a) the self-adjoint elements of $\pi(G)'$ are the $\Phi(w)$ with Hermitian $w_c$, a real space of
dimension $\sum_cm_c^2$.
(c) A quadratic polynomial $p(x)=x^TBx$ is determined by its symmetric matrix $B$. Put $(g\cdot p)(x)=p(\pi(g)^{-1}x)$;
then $g\cdot(h\cdot p)=(gh)\cdot p$, and $g\cdot p$ has the matrix $\sigma(g)(B)=A^TBA$, $A=\pi(g)^{-1}$. The
average $M=|G|^{-1}\sum_g\sigma(g)$ satisfies $\sigma(h)M=M$ and $MB=B$ for invariant $B$, so (as in the proof of
Theorem 12.5) it is idempotent with range the invariant polynomials, and their number is
$\operatorname{Tr}M=|G|^{-1}\sum_g\operatorname{Tr}\sigma(g)$. $A$ is unitary, hence $A=R\Lambda R^{-1}$ with
$\Lambda=\operatorname{diag}(\lambda)$ (Theorem 1.3). The symmetric matrices
$B_{ij}=(R^{-1})^T(e_ie_j^T+e_je_i^T)R^{-1}$, $i\le j$, form a basis, and
$A^TB_{ij}A=(R^{-1})^T\Lambda(e_ie_j^T+e_je_i^T)\Lambda R^{-1}=\lambda_i\lambda_jB_{ij}$. So
$\operatorname{Tr}\sigma(g)=\sum_{i\le j}\lambda_i\lambda_j=\frac12\big((\operatorname{Tr}A)^2+\operatorname{Tr}A^2\big)
=\frac12\big(\chi_\pi(g^{-1})^2+\chi_\pi(g^{-2})\big)$, and summing over $g^{-1}$ instead of $g$ gives the formula.
For a permutation representation, $G$ permutes the monomials $x_ax_b$, so a polynomial is invariant iff its
coefficients are constant on the orbits of monomials. $\square$

For graphs, the invariant edge features are spanned by the indicators of the diagonal and of its complement
($B_2=2$), and the invariant linear readouts are the sums of the diagonal and of the off-diagonal entries. The code
finds the trivial sector from the blocks $\rho_c(g)$ (`to_algebra`) and checks (a), (b) and the two counts of (c).

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

On $\mathbb C^4$ they are $\sum_ix_i^2$ and $\sum_{i<j}x_ix_j$. All degrees are counted by Molien's formula
$\sum_d\dim\mathbb C[x]^G_dt^d=|G|^{-1}\sum_g\det(1-t\pi(g))^{-1}$, whose coefficient of $t^2$ is (c) (stated without
proof; B. Sturmfels, *Algorithms in Invariant Theory*, Springer 2008, Theorem 2.2.1).

> **Physical meaning.** In Landau's theory of phase transitions the free energy is an invariant polynomial in a
> real order parameter $x$, on which $G$ acts by real orthogonal matrices. Its quadratic part is $x^TBx$ with $B$ in
> the commutant (as $\pi(g)^T=\pi(g)^{-1}$), so by Proposition 12.6(a) it is block diagonal over the sectors: the
> symmetric phase becomes unstable in the sector where an eigenvalue of $B$ first changes sign, and (generically)
> the irreducible representation of that sector determines the symmetry of the ordered phase.

## 12.5 What to compare with

Equivariant layers for specific groups are provided by mature libraries. **escnn** (Cesa, Lang and Weiler 2022,
successor of e2cnn) builds steerable convolutions for subgroups of $O(2)$ and $O(3)$ from irreducible
representations and harmonic bases of kernels; **e3nn** (Geiger and Smidt 2022) builds $E(3)$-equivariant networks
from spherical harmonics and Clebsch–Gordan tensor products. **EMLP** (Finzi, Welling and Wilson 2021) handles any
matrix group: it imposes $\pi_2(g)T=T\pi_1(g)$ for discrete generators and $d\pi_2(X)T=Td\pi_1(X)$ for Lie algebra
generators and computes a basis of the solutions as a numerical null space — the linear problem of
`dc.commutant_basis`. Group convolutions (Cohen and Welling 2016) are the case of the regular representation,
and Kondor and Trivedi (2018) show that, for feature maps on homogeneous spaces of a compact group, equivariant
linear layers are generalised convolutions. For permutations, DeepSets (Zaheer et al. 2017) and the bases of
Maron et al. (2019) are Corollary 12.2.

The library does not provide convolutions on grids, steerable kernels on continuous spaces, decompositions over
$\mathbb R$ (real, complex and quaternionic types; Exercise 2 treats real matrices), or large $n$: $V$ is a dense
$n\times n$ matrix and the decomposition costs $O(|G|n^3)$ per random sample for a group given by all its elements
(Lie algebra generators go through `dc.from_generators`; Section 4.3). What it adds is the
**structure** of the solution space: a basis of the commutant is a list of matrices, whereas the decomposition
identifies the commutant with the $*$-algebra $\bigoplus_cM_{m_c}$. Products, adjoints, unitarity, positivity and
the functional calculus then act blockwise on small matrices, with gradients (Proposition 12.6), and
`to_commutant`/`from_commutant` translate between the two pictures.

## Exercises

1. (**DeepSets with channels.**) $S_n$ ($n\ge2$) acts on $\mathbb C^n\otimes\mathbb C^d$ by $P_\sigma\otimes1_d$. Show
   that the equivariant maps from $\mathbb C^n\otimes\mathbb C^{d_1}$ to $\mathbb C^n\otimes\mathbb C^{d_2}$ are
   $1\otimes A+J\otimes B$ with $A,B\in M_{d_2\times d_1}$, and check it with `EquivariantMap`.
2. (**Real layers.**) Let $\pi_1,\pi_2$ have real matrices (e.g. permutation representations). Show that the real
   equivariant matrices form a real vector space of dimension $\dim_{\mathbb C}\operatorname{Hom}_G(\pi_1,\pi_2)$, and
   that $T\mapsto\operatorname{Re}T$ maps $\operatorname{Hom}_G(\pi_1,\pi_2)$ onto it.
3. (**Nonlinearities.**) (a) Show that applying a function $f:\mathbb C\to\mathbb C$ to every coordinate commutes
   with every permutation representation. (b) Show that it does not commute with the rotation by $\pi/2$ of
   $\mathbb R^2$ for $f(t)=\max(t,0)$. (c) Show that $x\mapsto\varphi(\|x\|)x$ commutes with every unitary
   representation.
4. (**How many samples determine a layer?**) Let $T=\Phi(w)\in\pi(G)'$ and $x_1,\dots,x_N\in\mathbb C^n$. Write
   $V^*x_i$ in sector $c$ as $\sum_{a,j}(X_{i,c})_{aj}e_a\otimes e_j$ with $X_{i,c}\in M_{k_c\times m_c}$, and let
   $S_c\in M_{Nk_c\times m_c}$ be the matrices $X_{1,c},\dots,X_{N,c}$ stacked. Show that the equivariant $T$ with
   $Tx_1=\dots=Tx_N=0$ form a space of dimension $\sum_cm_c(m_c-\operatorname{rank}S_c)$; for Gaussian $x_i$ this is
   $\sum_cm_c\max(0,m_c-Nk_c)$ almost surely. For the edge representation, two samples determine an equivariant
   layer, while an arbitrary $16\times16$ matrix needs sixteen.

<details>
<summary>Solutions</summary>

1. $P_\sigma J=J=JP_\sigma$, so these maps are equivariant. They are linearly independent: if
   $1\otimes A+J\otimes B=0$, the off-diagonal blocks give $B=0$ and then $A=0$. By Corollary 12.2(a) and
   Theorem 2.6, $P_\sigma\otimes1_d\cong(\mathrm{triv}\otimes1_d)\oplus(\mathrm{std}\otimes1_d)$, so both irreducible
   representations have multiplicity $d$, and Theorem 12.4 gives $\dim\operatorname{Hom}_G=2d_1d_2$, the dimension of
   the space of the $1\otimes A+J\otimes B$.

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

2. If $T$ is equivariant, so is $\bar T$, since $\overline{\pi_i(g)}=\pi_i(g)$. Hence
   $\operatorname{Re}T=(T+\bar T)/2$ and $\operatorname{Im}T=(T-\bar T)/2i$ are real and equivariant, so the map
   $T\mapsto\operatorname{Re}T$ is onto the real equivariant matrices (it is the identity on them), and
   $\operatorname{Hom}_G(\pi_1,\pi_2)=H_{\mathbb R}\oplus iH_{\mathbb R}$ with $H_{\mathbb R}$ the real equivariant
   matrices. A real basis of $H_{\mathbb R}$ is then a complex basis of $\operatorname{Hom}_G$: it spans, and a
   complex relation $\sum_j(\alpha_j+i\beta_j)T_j=0$ with real $T_j$ splits into real and imaginary parts. In
   Lemma 12.1 the orbit indicators are such a real basis. (The blocks $y_c$ of Theorem 12.4 are complex even for real
   data: $V$ is complex in general, and the parametrisation by $\operatorname{Re}T$ uses twice as many real parameters
   as needed.)

3. (a) $(\pi(g)f(x))_{gy}=f(x)_y=f(x_y)=f((\pi(g)x)_{gy})$. (b) For $R(x_1,x_2)=(-x_2,x_1)$ and $x=(1,-1)$,
   $f(Rx)=f(1,1)=(1,1)$ but $Rf(x)=R(1,0)=(0,1)$. (c) $\|\pi(g)x\|=\|x\|$, so
   $\varphi(\|\pi(g)x\|)\pi(g)x=\pi(g)\varphi(\|x\|)x$.

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

4. $(1\otimes w_c)(e_a\otimes e_j)=\sum_le_a\otimes(w_c)_{lj}e_l$, so $V^*Tx_i$ has, in sector $c$, the coefficient
   matrix $X_{i,c}w_c^T$. Hence $Tx_i=0$ for all $i$ iff $S_cw_c^T=0$ for all $c$, i.e. iff every column of $w_c^T$ lies in
   $\ker S_c$, a space of dimension $m_c-\operatorname{rank}S_c$: the $w_c$ with this property form a space of
   dimension $m_c(m_c-\operatorname{rank}S_c)$, and $\Phi$ is injective. For Gaussian $x_i$, $V^*x_i$ is Gaussian
   ($V$ is unitary), so $S_c$ has independent Gaussian entries and rank $\min(Nk_c,m_c)$ almost surely: for the
   leading square submatrix of that size, $|\det|^2$ is a real polynomial in the real and imaginary parts of the
   entries, not identically zero (it is $1$ when the submatrix is the identity), and such a polynomial vanishes
   only on a Lebesgue null set (proved in Section 4.3). For the edge representation, $(k,m)=(1,2),(2,1),(3,1),(3,3)$,
   only the sector $(1,2)$ contributes: $2\cdot(2-1)=2$ for $N=1$, and $0$ for $N=2$.

```python
basis = D_edge.basis('commutant')                                # (15, 16, 16), Hilbert-Schmidt orthonormal
for N, lost in ((1, 2), (2, 0)):
    Xs = torch.randn(16, N, dtype=dt)
    M = torch.stack([(B @ Xs).reshape(-1) for B in basis], dim=1)               # w -> (T x_1, ..., T x_N)
    assert torch.linalg.matrix_rank(M).item() == 15 - lost
```

</details>

**Where to go from here.** The last chapter, [Chapter 13](13_jordan_algebras.md), asks why quantum mechanics uses *-algebras rather than Jordan algebras. [Chapter 4](04_symmetries.md) develops the
representation theory used here: Schur's lemma, characters, and the randomised algorithm behind `dc.from_group`
(Sections 4.1–4.3). [Chapter 9](09_spin_chains.md) applies the same sectors to spin chains, and the Haar measure on
$\prod_cU(m_c)$ of [Chapter 5](05_random_matrices.md) gives the natural random initialisation of an equivariant
unitary layer. For the mathematics see J.-P. Serre, *Linear Representations of Finite Groups* (Springer 1977) and
W. Fulton and J. Harris, *Representation Theory* (Springer 1991); for equivariant learning, M. Bronstein, J. Bruna,
T. Cohen and P. Veličković, *Geometric Deep Learning* (2021), and M. Weiler et al., *Equivariant and Coordinate
Independent Convolutional Networks* (2023); for the decomposition algorithm, the references in
`torch_vn_algebra/decompose.py`.
