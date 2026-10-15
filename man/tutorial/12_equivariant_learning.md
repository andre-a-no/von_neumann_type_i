**English** · [Español](es/12_equivariant_learning.md) · [Français](fr/12_equivariant_learning.md) · [Deutsch](de/12_equivariant_learning.md) · [中文](zh/12_equivariant_learning.md) · [日本語](ja/12_equivariant_learning.md) · [Русский](ru/12_equivariant_learning.md)

# Tutorial 12. Equivariant machine learning

This tutorial accompanies [Chapter 12](../lecture_notes/12_equivariant_learning.md) of the lecture notes, which proves that the equivariant linear
maps between two representations of a finite group are block matrices ([Theorem 12.4](../lecture_notes/12_equivariant_learning.md#r-12-4)), that
averaging over the group is the orthogonal projection onto them ([Theorem 12.5](../lecture_notes/12_equivariant_learning.md#r-12-5)), and counts invariants.
Here we build and train equivariant layers with `torch_vn_algebra.decompose`, for $`S_4`$ acting on the vertex and
edge features of graphs, and check these results numerically.

In the code a batch of vectors is a tensor of shape `(batch, n)`, on which a matrix $`T`$ acts as `x @ T.T`.

<a id="s-12-1"></a>

## 12.1 Equivariant linear maps

*Lecture notes:* [Section 12.1](../lecture_notes/12_equivariant_learning.md#s-12-1): [Lemma 12.1](../lecture_notes/12_equivariant_learning.md#r-12-1) (orbits) and [Corollary 12.2](../lecture_notes/12_equivariant_learning.md#r-12-2) (DeepSets, Bell numbers).

The DeepSets layer (Zaheer et al. 2017) is [lecture notes, Corollary 12.2](../lecture_notes/12_equivariant_learning.md#r-12-2)(a), the bases of Maron et al. (2019) for graph networks are (b).
The code counts the orbits by brute force, compares with the commutant computed numerically ([Lemma 12.1](../lecture_notes/12_equivariant_learning.md#r-12-1)),
and decomposes both representations of $`S_4`$:

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

So, in the decomposition computed by `dc.from_group` (the asserts check its output), $`\mathbb C^{16}`$ contains a
representation of dimension $`1`$ twice — the trivial one, whose multiplicity is the number $`2`$ of orbits on $`[4]^2`$ by
Proposition 12.7(a) and Lemma 12.1 (the code of Section 12.4 identifies this sector) — one of dimension $`2`$ once, and
two of dimension $`3`$ once and three times:
$`16=1\cdot2+2\cdot1+3\cdot1+3\cdot3`$, and $`15=2^2+1^2+1^2+3^2`$ agrees with the orbit count of Lemma 12.1 and
$`\dim\pi(G)'=\sum_cm_c^2`$ (Theorem 2.6). For two representations, everything reduces to one commutant:

By [lecture notes, Lemma 12.3](../lecture_notes/12_equivariant_learning.md#r-12-3), $`\mathrm{Hom}_G(\pi_1,\pi_2)=Q_2^*\,\pi(G)'\,Q_1`$ for $`\pi=\pi_1\oplus\pi_2`$, and by
[Theorem 12.4](../lecture_notes/12_equivariant_learning.md#r-12-4) it is, in suitable bases, $`\bigoplus_cM_{b_c\times a_c}(\mathbb C)`$, of dimension $`\sum_\rho m_\rho(\pi_1)\,m_\rho(\pi_2)`$.

<a id="s-12-2"></a>

## 12.2 An equivariant layer with the library

*Lecture notes:* [Section 12.2](../lecture_notes/12_equivariant_learning.md#s-12-2): [Theorem 12.5](../lecture_notes/12_equivariant_learning.md#r-12-5) (the twirl); the layers follow the proofs of [Lemma 12.3](../lecture_notes/12_equivariant_learning.md#r-12-3) and [Theorem 12.4](../lecture_notes/12_equivariant_learning.md#r-12-4).

The proof of Theorem 12.4 is an algorithm, once the decomposition of Theorem 2.6 is known; `dc.from_group`
computes it by the randomised method of Section 4.3 (correct with probability one in exact arithmetic; here the
asserts check the computed output). For $`\pi_1=\pi_2`$ it is just Theorem 2.6: the weights are blocks
$`w_c\in M_{m_c}`$, an operator of the abstract commutant `D.commutant_algebra()` $`=\bigoplus_cM_{m_c}`$, and the
layer is $`x\mapsto V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*x`$, computed by `D.from_commutant`. The module
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

For $`\pi_1\ne\pi_2`$ the proof is followed line by line: decompose $`\pi_1\oplus\pi_2`$, compute the blocks $`p_c`$ of
$`P_1`$ with `to_commutant`, diagonalise them, and keep the lower left corners. The trace of $`p_c`$ is
$`a_c=m_c(\pi_1)`$.

For Corollary 12.2(b) with $`n=4`$, in the computed decomposition (printed and asserted below): vertex features
$`\pi_1`$ have $`a=(1,0,0,1)`$ in the four sectors of $`\pi_1\oplus\pi_2`$ (dimensions $`1,2,3,3`$; the first is the trivial
one, since by Corollary 12.2(a) the only one-dimensional constituent of $`\pi_1`$ is trivial), edge features $`\pi_2`$
have $`b=(2,1,1,3)`$, so $`\dim\mathrm{Hom}_G=1\cdot2+1\cdot3=5=B_3`$, as Corollary 12.2(b) proves.

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

The second way to make a layer equivariant is to average it over the group. The **twirl** is
```math
\mathcal P(W)=\frac1{|G|}\sum_{g\in G}\pi_2(g)\,W\,\pi_1(g)^*,\qquad W\in M_{n_2\times n_1}(\mathbb C).
```

By [lecture notes, Theorem 12.5](../lecture_notes/12_equivariant_learning.md#r-12-5), $`\mathcal P`$ is the Hilbert–Schmidt orthogonal projection onto $`\mathrm{Hom}_G(\pi_1,\pi_2)`$, and
$`\dim\mathrm{Hom}_G(\pi_1,\pi_2)=|G|^{-1}\sum_g\overline{\chi_{\pi_1}(g)}\chi_{\pi_2}(g)`$; for $`\pi_1=\pi_2`$, `dc.commutant_dimension` evaluates this
character formula. The parametrisation $`W\mapsto\mathcal P(W)`$ wastes $`n_1n_2-\dim\mathrm{Hom}_G`$ parameters on the kernel ($`256-15`$ for edge features), and an
evaluation costs $`O(|G|n_1n_2(n_1+n_2))`$ for dense $`\pi(g)`$, with $`|G|=n!`$ for $`S_n`$. The library computes the same
projection from the decomposition, `D.conditional_expectation(W, 'commutant')`, with four $`n\times n`$ products.
The code checks that `dc.twirl` and `D.conditional_expectation` agree, that $`\mathcal P`$ is the identity on $`\pi(G)'`$
and self-adjoint, and the dimension formula from vertex to edge features.

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

## 12.3 Constrained and spectral layers

*Lecture notes:* [Section 12.3](../lecture_notes/12_equivariant_learning.md#s-12-3): [Proposition 12.6](../lecture_notes/12_equivariant_learning.md#r-12-6).

The weights $`w=(w_c)_c`$ of `EquivariantLinear` live in a $`*`$-algebra, so every construction of Chapters 1 and 11
applies to them. Write $`\Phi(w)=V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*`$ for $`w\in\bigoplus_cM_{m_c}`$.

By [lecture notes, Proposition 12.6](../lecture_notes/12_equivariant_learning.md#r-12-6), (a) $`\Phi`$ is a $`*`$-isomorphism of $`\bigoplus_cM_{m_c}`$ onto $`\pi(G)'`$, and $`\Phi(w)`$ is
self-adjoint, unitary or positive iff every $`w_c`$ is; (b) $`f(\Phi(w))=\Phi\big((f(w_c))_c\big)`$ for normal $`w`$; (c) the
equivariant unitaries are the $`\Phi(u_0e^{s-s^*})`$.

So an orthogonal (unitary) equivariant layer is `optimize.UnitaryParam` on `D.commutant_algebra()` followed by
$`\Phi`$; by (c) this parametrisation reaches every equivariant unitary, and in exact arithmetic every iterate of a
gradient method is unitary and equivariant (the code checks both to $`10^{-12}`$ after training). Likewise
`PositiveParam` gives positive equivariant layers, normalised as
in Chapter 11 (spectrum in $`[\ell,1]`$ with maximum $`1`$; by (a), a general positive layer is $`\Phi(t^*t)`$), and a
spectral layer $`f(L)`$ of a self-adjoint equivariant generator $`L`$, such as the diffusion $`e^{-tL}`$, needs
eigendecompositions of the blocks $`w_c`$, of cost of order $`\sum_cm_c^3`$ ($`37`$ here) instead of $`n^3=4096`$; the library
differentiates through `apply_function` (we prove nothing about its backward pass; the code checks one derivative
below). The code fits an equivariant unitary to $`32`$ input–output pairs with Adam — a numerical observation of this
run, not a theorem: the loss decreases by more than a factor $`10^3`$, which the assert checks, and no convergence
statement is claimed — and checks (b) and the gradient with respect to $`t`$ against dense matrices and the formula
$`\frac d{dt}e^{-tL}=-Le^{-tL}`$ (for $`L=U\mathrm{diag}(\lambda)U^*`$, $`e^{-tL}=U\mathrm{diag}(e^{-t\lambda_i})U^*`$ by
Theorem 1.5, and $`U`$ does not depend on $`t`$).

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

## 12.4 Invariants

*Lecture notes:* [Section 12.4](../lecture_notes/12_equivariant_learning.md#s-12-4): [Proposition 12.7](../lecture_notes/12_equivariant_learning.md#r-12-7) (invariant vectors and quadratic forms) and [Proposition 12.8](../lecture_notes/12_equivariant_learning.md#r-12-8) (Molien's formula).

For graphs, the invariant edge features are spanned by the indicators of the diagonal and of its complement
($`B_2=2`$), and the invariant linear readouts are the sums of the diagonal and of the off-diagonal entries. The code
finds the trivial sector from the blocks $`\rho_c(g)`$ (`to_algebra`) and checks (a), (b) and the two counts of (c).

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

On $`\mathbb C^4`$ they are $`\sum_ix_i^2`$ and $`\sum_{i<j}x_ix_j`$. The proof of (c) works in every degree:

By [lecture notes, Proposition 12.8](../lecture_notes/12_equivariant_learning.md#r-12-8) (Molien's formula),
$`\sum_{d\ge0}\dim\mathbb C[x]^G_d\,t^d=|G|^{-1}\sum_g\det\big(1-t\pi(g)\big)^{-1}`$, and for a permutation representation
$`\dim\mathbb C[x]^G_d`$ is the number of orbits of $`G`$ on the monomials of degree $`d`$.

(This is Molien's theorem of 1897; see B. Sturmfels, *Algorithms in Invariant Theory*, Springer 2008, Section 2.2,
for further reading.) The code expands the right side up to $`t^4`$, multiplying by one geometric series per
eigenvalue, and compares with the orbits of monomials; for $`S_4`$ on $`\mathbb C^4`$ the counts $`1,1,2,3,5`$ are the
partitions of $`d`$ into at most four parts (the orbits of multisets are determined by their multiplicities, sorted).

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

## 12.5 What to compare with

Equivariant layers for specific groups are provided by mature libraries. **escnn** (Cesa, Lang and Weiler 2022,
successor of e2cnn) builds steerable convolutions for subgroups of $`O(2)`$ and $`O(3)`$ from irreducible
representations and harmonic bases of kernels; **e3nn** (Geiger and Smidt 2022) builds $`E(3)`$-equivariant networks
from spherical harmonics and Clebsch–Gordan tensor products. **EMLP** (Finzi, Welling and Wilson 2021) handles any
matrix group: it imposes $`\pi_2(g)T=T\pi_1(g)`$ for discrete generators and $`d\pi_2(X)T=Td\pi_1(X)`$ for Lie algebra
generators and computes a basis of the solutions as a numerical null space — the linear problem of
`dc.commutant_basis`. Group convolutions (Cohen and Welling 2016) are the case $`\Omega_1=\Omega_2=G`$ of Lemma 12.1,
with $`G`$ acting on itself by left multiplication: $`(y,x)`$ and $`(y',x')`$ lie in the same orbit iff
$`x^{-1}y=x'^{-1}y'`$ (if so, $`g=x'x^{-1}`$ maps one pair to the other), so the equivariant maps are
$`(Tv)(y)=\sum_xf(x^{-1}y)v(x)`$, the convolutions with functions $`f`$ on $`G`$. Kondor and Trivedi (2018) study the
analogous question for homogeneous spaces of compact groups, which this course does not treat. For permutations,
DeepSets (Zaheer et al. 2017) and the bases of Maron et al. (2019) are Corollary 12.2.

The library does not provide convolutions on grids, steerable kernels on continuous spaces, decompositions over
$`\mathbb R`$ (real, complex and quaternionic types; Exercise 2 treats real matrices), or large $`n`$: $`V`$ is a dense
$`n\times n`$ matrix and the decomposition costs $`O(|G|n^3)`$ per random sample for a group given by all its elements
(Lie algebra generators go through `dc.from_generators`; Section 4.3). What it adds is the
**structure** of the solution space: a basis of the commutant is a list of matrices, whereas the decomposition
identifies the commutant with the $`*`$-algebra $`\bigoplus_cM_{m_c}`$. Products, adjoints, unitarity, positivity and
the functional calculus then act blockwise on small matrices (Proposition 12.6), the library differentiates
through them, and `to_commutant`/`from_commutant` translate between the two pictures.

<a id="x-12"></a>

## Exercises

[Exercises](../lecture_notes/12_equivariant_learning.md#x-12) 1, 3 and 4 of the lecture notes, numerically:

- *Exercise 1* (DeepSets with channels): check with `EquivariantMap` that the equivariant maps from
  $`\mathbb C^n\otimes\mathbb C^{d_1}`$ to $`\mathbb C^n\otimes\mathbb C^{d_2}`$ are the $`1\otimes A+J\otimes B`$.
- *Exercise 3* (nonlinearities): ReLU on the coordinates commutes with the permutation representation but not
  with the rotation by $`\pi/2`$ of $`\mathbb R^2`$; $`x\mapsto\varphi(\|x\|)x`$ commutes with every unitary representation.
- *Exercise 4*: for the edge representation, one Gaussian sample leaves a $`2`$-dimensional space of equivariant
  layers undetermined, two samples determine the layer.

<details>
<summary>Solutions</summary>

*Exercise 1.*

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

*Exercise 3.*

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

*Exercise 4.*

```python
basis = D_edge.basis('commutant')                                # (15, 16, 16), Hilbert-Schmidt orthonormal
for N, lost in ((1, 2), (2, 0)):
    Xs = torch.randn(16, N, dtype=dt)
    M = torch.stack([(B @ Xs).reshape(-1) for B in basis], dim=1)               # w -> (T x_1, ..., T x_N)
    assert torch.linalg.matrix_rank(M).item() == 15 - lost
```

</details>

**Next:** [Tutorial 13. Jordan algebras: why quantum mechanics uses *-algebras](13_jordan_algebras.md).
