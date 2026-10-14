[English](../12_equivariant_learning.md) · [Español](../es/12_equivariant_learning.md) · **Français** · [Deutsch](../de/12_equivariant_learning.md) · [中文](../zh/12_equivariant_learning.md) · [日本語](../ja/12_equivariant_learning.md) · [Русский](../ru/12_equivariant_learning.md)

# Tutoriel 12. Apprentissage automatique équivariant

Ce tutoriel accompagne le [chapitre 12](../../lecture_notes/fr/12_equivariant_learning.md) des notes de cours, qui démontre que les applications linéaires
équivariantes entre deux représentations d'un groupe fini sont des matrices par blocs ([théorème 12.4](../../lecture_notes/fr/12_equivariant_learning.md#r-12-4)), que
la moyenne sur le groupe est la projection orthogonale sur ces applications ([théorème 12.5](../../lecture_notes/fr/12_equivariant_learning.md#r-12-5)), et qui compte les invariants.
Nous construisons et entraînons ici des couches équivariantes avec `torch_vn_algebra.decompose`, pour $S_4$ agissant sur les
caractéristiques de sommets et d'arêtes de graphes, et vérifions numériquement ces résultats.

Dans le code, un lot de vecteurs est un tenseur de forme `(batch, n)`, sur lequel une matrice $T$ agit par `x @ T.T`.

<a id="s-12-1"></a>

## 12.1 Applications linéaires équivariantes

*Notes de cours :* [section 12.1](../../lecture_notes/fr/12_equivariant_learning.md#s-12-1) : [lemme 12.1](../../lecture_notes/fr/12_equivariant_learning.md#r-12-1) (orbites) et [corollaire 12.2](../../lecture_notes/fr/12_equivariant_learning.md#r-12-2) (DeepSets, nombres de Bell).

La couche DeepSets (Zaheer et al. 2017) est [notes de cours, corollaire 12.2](../../lecture_notes/fr/12_equivariant_learning.md#r-12-2)(a), les bases de Maron et al. (2019) pour les réseaux de graphes sont (b).
Le code compte les orbites par force brute, les compare au commutant calculé numériquement ([lemme 12.1](../../lecture_notes/fr/12_equivariant_learning.md#r-12-1))
et décompose les deux représentations de $S_4$ :

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

Ainsi, dans la décomposition calculée par `dc.from_group` (les asserts vérifient son résultat), $\mathbb C^{16}$
contient une représentation de dimension $1$ deux fois — la triviale, dont la multiplicité est le nombre $2$
d'orbites sur $[4]^2$ d'après la proposition 12.7(a) et le lemme 12.1 (le code de la section 12.4 identifie ce
secteur) —, une représentation de dimension $2$ une fois, et deux représentations de dimension $3$, l'une une fois
et l'autre trois fois : $16=1\cdot2+2\cdot1+3\cdot1+3\cdot3$, et $15=2^2+1^2+1^2+3^2$ concorde avec le décompte des
orbites du lemme 12.1 et avec $\dim\pi(G)'=\sum_cm_c^2$ (théorème 2.6). Pour deux représentations, tout se ramène à
un seul commutant :

D'après [notes de cours, lemme 12.3](../../lecture_notes/fr/12_equivariant_learning.md#r-12-3), $\mathrm{Hom}_G(\pi_1,\pi_2)=Q_2^*\,\pi(G)'\,Q_1$ pour $\pi=\pi_1\oplus\pi_2$, et d'après le
[théorème 12.4](../../lecture_notes/fr/12_equivariant_learning.md#r-12-4), c'est, dans des bases convenables, $\bigoplus_cM_{b_c\times a_c}(\mathbb C)$, de dimension $\sum_\rho m_\rho(\pi_1)\,m_\rho(\pi_2)$.

<a id="s-12-2"></a>

## 12.2 Une couche équivariante avec la bibliothèque

*Notes de cours :* [section 12.2](../../lecture_notes/fr/12_equivariant_learning.md#s-12-2) : [théorème 12.5](../../lecture_notes/fr/12_equivariant_learning.md#r-12-5) (le twirl) ; les couches suivent les démonstrations du [lemme 12.3](../../lecture_notes/fr/12_equivariant_learning.md#r-12-3) et du [théorème 12.4](../../lecture_notes/fr/12_equivariant_learning.md#r-12-4).

La démonstration du théorème 12.4 est un algorithme, dès que l'on connaît la décomposition du théorème 2.6 ;
`dc.from_group` la calcule par la méthode randomisée de la section 4.3 (correcte avec probabilité un en
arithmétique exacte ; ici, les asserts vérifient le résultat calculé). Pour $\pi_1=\pi_2$, c'est simplement le
théorème 2.6 : les poids sont des blocs $w_c\in M_{m_c}$, un opérateur du commutant abstrait
`D.commutant_algebra()` $=\bigoplus_cM_{m_c}$, et la couche est $x\mapsto V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*x$,
calculée par `D.from_commutant`. Le module ci-dessous est du code du chapitre, il ne fait pas partie de la
bibliothèque.

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

Pour $\pi_1\ne\pi_2$, on suit la démonstration ligne à ligne : on décompose $\pi_1\oplus\pi_2$, on calcule les blocs
$p_c$ de $P_1$ avec `to_commutant`, on les diagonalise et l'on garde les coins inférieurs gauches. La trace de $p_c$
est $a_c=m_c(\pi_1)$.

Pour le corollaire 12.2(b) avec $n=4$, dans la décomposition calculée (affichée et vérifiée ci-dessous) : les
caractéristiques de sommets $\pi_1$ ont $a=(1,0,0,1)$ dans les quatre secteurs de $\pi_1\oplus\pi_2$ (de dimensions
$1,2,3,3$ ; le premier est le secteur trivial, car d'après le corollaire 12.2(a) le seul constituant de dimension un
de $\pi_1$ est trivial), les caractéristiques d'arêtes $\pi_2$ ont $b=(2,1,1,3)$, donc
$\dim\mathrm{Hom}_G=1\cdot2+1\cdot3=5=B_3$, comme le démontre le corollaire 12.2(b).

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

La seconde façon de rendre une couche équivariante est de la moyenner sur le groupe. Le **twirl** est
$$
\mathcal P(W)=\frac1{|G|}\sum_{g\in G}\pi_2(g)\,W\,\pi_1(g)^*,\qquad W\in M_{n_2\times n_1}(\mathbb C).
$$

D'après [notes de cours, théorème 12.5](../../lecture_notes/fr/12_equivariant_learning.md#r-12-5), $\mathcal P$ est la projection orthogonale, pour le produit scalaire de Hilbert–Schmidt, sur $\mathrm{Hom}_G(\pi_1,\pi_2)$, et
$\dim\mathrm{Hom}_G(\pi_1,\pi_2)=|G|^{-1}\sum_g\overline{\chi_{\pi_1}(g)}\chi_{\pi_2}(g)$ ; pour $\pi_1=\pi_2$, `dc.commutant_dimension` évalue cette
formule des caractères. Le paramétrage $W\mapsto\mathcal P(W)$ gaspille $n_1n_2-\dim\mathrm{Hom}_G$ paramètres dans le noyau ($256-15$ pour les caractéristiques d'arêtes), et une
évaluation coûte $O(|G|n_1n_2(n_1+n_2))$ pour des $\pi(g)$ denses, avec $|G|=n!$ pour $S_n$. La bibliothèque calcule la même
projection à partir de la décomposition, `D.conditional_expectation(W, 'commutant')`, avec quatre produits $n\times n$.
Le code vérifie que `dc.twirl` et `D.conditional_expectation` coïncident, que $\mathcal P$ est l'identité sur $\pi(G)'$
et autoadjointe, ainsi que la formule de dimension des caractéristiques de sommets vers les caractéristiques d'arêtes.

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

## 12.3 Couches contraintes et spectrales

*Notes de cours :* [section 12.3](../../lecture_notes/fr/12_equivariant_learning.md#s-12-3) : [proposition 12.6](../../lecture_notes/fr/12_equivariant_learning.md#r-12-6).

Les poids $w=(w_c)_c$ de `EquivariantLinear` vivent dans une $*$-algèbre, de sorte que toutes les constructions des
chapitres 1 et 11 s'y appliquent. Notons $\Phi(w)=V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*$ pour
$w\in\bigoplus_cM_{m_c}$.

D'après [notes de cours, proposition 12.6](../../lecture_notes/fr/12_equivariant_learning.md#r-12-6), (a) $\Phi$ est un $*$-isomorphisme de $\bigoplus_cM_{m_c}$ sur $\pi(G)'$, et $\Phi(w)$ est
autoadjoint, unitaire ou positif si et seulement si chaque $w_c$ l'est ; (b) $f(\Phi(w))=\Phi\big((f(w_c))_c\big)$ pour $w$ normal ; (c) les
unitaires équivariants sont les $\Phi(u_0e^{s-s^*})$.

Ainsi une couche équivariante orthogonale (unitaire) est `optimize.UnitaryParam` sur `D.commutant_algebra()` suivi de
$\Phi$ ; d'après (c), ce paramétrage atteint tout unitaire équivariant, et en arithmétique exacte chaque itéré d'une
méthode de gradient est unitaire et équivariant (le code vérifie les deux à $10^{-12}$ près après l'entraînement).
De même, `PositiveParam` fournit des couches équivariantes positives,
normalisées comme au chapitre 11 (spectre dans $[\ell,1]$ de maximum $1$ ; d'après (a), une couche positive générale
est $\Phi(t^*t)$), et une couche spectrale $f(L)$ d'un générateur équivariant autoadjoint $L$, comme la diffusion
$e^{-tL}$, demande des diagonalisations des blocs $w_c$, d'un coût de l'ordre de $\sum_cm_c^3$ (ici $37$) au lieu de
$n^3=4096$ ; la bibliothèque différencie à travers `apply_function` (nous ne démontrons rien sur sa passe arrière ;
le code vérifie une dérivée ci-dessous). Le code ajuste un unitaire équivariant à $32$ couples entrée–sortie avec
Adam — une observation numérique de cette exécution, pas un théorème : la perte diminue d'un facteur supérieur à
$10^3$, ce que vérifie l'assert, et aucun énoncé de convergence n'est revendiqué — et confronte (b) et le gradient par
rapport à $t$ aux matrices denses et à la formule $\frac d{dt}e^{-tL}=-Le^{-tL}$ (pour
$L=U\mathrm{diag}(\lambda)U^*$, $e^{-tL}=U\mathrm{diag}(e^{-t\lambda_i})U^*$ d'après le théorème 1.5, et
$U$ ne dépend pas de $t$).

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

*Notes de cours :* [section 12.4](../../lecture_notes/fr/12_equivariant_learning.md#s-12-4) : [proposition 12.7](../../lecture_notes/fr/12_equivariant_learning.md#r-12-7) (vecteurs invariants et formes quadratiques) et [proposition 12.8](../../lecture_notes/fr/12_equivariant_learning.md#r-12-8) (formule de Molien).

Pour les graphes, les caractéristiques d'arêtes invariantes sont engendrées par les indicatrices de la diagonale et de
son complémentaire ($B_2=2$), et les lectures (*readouts*) linéaires invariantes sont la somme des coefficients
diagonaux et celle des coefficients hors diagonale. Le code trouve le secteur trivial à partir des blocs
$\rho_c(g)$ (`to_algebra`) et vérifie (a), (b) et les deux décomptes de (c).

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

Sur $\mathbb C^4$, ce sont $\sum_ix_i^2$ et $\sum_{i<j}x_ix_j$. La démonstration de (c) fonctionne en tout degré :

D'après [notes de cours, proposition 12.8](../../lecture_notes/fr/12_equivariant_learning.md#r-12-8) (formule de Molien),
$\sum_{d\ge0}\dim\mathbb C[x]^G_d\,t^d=|G|^{-1}\sum_g\det\big(1-t\pi(g)\big)^{-1}$, et pour une représentation de permutation,
$\dim\mathbb C[x]^G_d$ est le nombre d'orbites de $G$ sur les monômes de degré $d$.

(C'est le théorème de Molien de 1897 ; pour aller plus loin, voir B. Sturmfels, *Algorithms in Invariant Theory*,
Springer 2008, section 2.2.) Le code développe le membre de droite jusqu'à $t^4$, en multipliant par une série
géométrique par valeur propre, et le compare aux orbites de monômes ; pour $S_4$ sur $\mathbb C^4$, les nombres
$1,1,2,3,5$ sont les partitions de $d$ en au plus quatre parts (les orbites de multiensembles sont déterminées par
leurs multiplicités, triées).

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

## 12.5 Points de comparaison

Des bibliothèques matures fournissent des couches équivariantes pour des groupes particuliers. **escnn** (Cesa, Lang
et Weiler 2022, successeur de e2cnn) construit des convolutions orientables (*steerable*) pour des sous-groupes de
$O(2)$ et $O(3)$ à partir de représentations irréductibles et de bases harmoniques de noyaux ; **e3nn** (Geiger et
Smidt 2022) construit des réseaux $E(3)$-équivariants à partir d'harmoniques sphériques et de produits tensoriels de
Clebsch–Gordan. **EMLP** (Finzi, Welling et Wilson 2021) traite n'importe quel groupe de matrices : il impose
$\pi_2(g)T=T\pi_1(g)$ pour des générateurs discrets et $d\pi_2(X)T=Td\pi_1(X)$ pour des générateurs de l'algèbre de
Lie, et calcule une base des solutions comme un noyau numérique — le problème linéaire de `dc.commutant_basis`. Les
convolutions de groupe (Cohen et Welling 2016) sont le cas $\Omega_1=\Omega_2=G$ du lemme 12.1, $G$ agissant sur
lui-même par multiplication à gauche : $(y,x)$ et $(y',x')$ sont dans la même orbite si et seulement si
$x^{-1}y=x'^{-1}y'$ (dans ce cas, $g=x'x^{-1}$ envoie une paire sur l'autre), donc les applications équivariantes sont
les $(Tv)(y)=\sum_xf(x^{-1}y)v(x)$, les convolutions par des fonctions $f$ sur $G$. Kondor et Trivedi (2018) étudient
la question analogue pour les espaces homogènes de groupes compacts, que ce cours ne traite pas. Pour les permutations, DeepSets (Zaheer et al. 2017) et
les bases de Maron et al. (2019) sont le corollaire 12.2.

La bibliothèque ne fournit ni convolutions sur des grilles, ni noyaux orientables sur des espaces continus, ni
décompositions sur $\mathbb R$ (types réel, complexe et quaternionique ; l'exercice 2 traite les matrices réelles), et
ne vise pas les grands $n$ : $V$ est une matrice dense $n\times n$ et la décomposition coûte $O(|G|n^3)$ par
échantillon aléatoire pour un groupe donné par tous ses éléments (les générateurs d'algèbres de Lie passent par
`dc.from_generators` ; section 4.3). Ce qu'elle apporte, c'est la **structure** de l'espace des solutions : une base du
commutant est une liste de matrices, tandis que la décomposition identifie le commutant à la $*$-algèbre
$\bigoplus_cM_{m_c}$. Les produits, les adjoints, l'unitarité, la positivité et le calcul fonctionnel agissent alors
bloc par bloc sur de petites matrices (proposition 12.6), la bibliothèque différencie à travers eux, et
`to_commutant`/`from_commutant` font la traduction entre les deux points de vue.

<a id="x-12"></a>

## Exercices

[Exercices](../../lecture_notes/fr/12_equivariant_learning.md#x-12) 1, 3 et 4 des notes de cours, numériquement :

- *Exercice 1* (DeepSets avec canaux) : vérifier avec `EquivariantMap` que les applications équivariantes de
  $\mathbb C^n\otimes\mathbb C^{d_1}$ dans $\mathbb C^n\otimes\mathbb C^{d_2}$ sont les $1\otimes A+J\otimes B$.
- *Exercice 3* (non-linéarités) : ReLU appliquée aux coordonnées commute avec la représentation de permutation, mais pas
  avec la rotation d'angle $\pi/2$ de $\mathbb R^2$ ; $x\mapsto\varphi(\|x\|)x$ commute avec toute représentation unitaire.
- *Exercice 4* : pour la représentation des arêtes, un échantillon gaussien laisse indéterminé un espace de couches
  équivariantes de dimension $2$ ; deux échantillons déterminent la couche.

<details>
<summary>Solutions</summary>

*Exercice 1.*

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

*Exercice 3.*

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

*Exercice 4.*

```python
basis = D_edge.basis('commutant')                                # (15, 16, 16), Hilbert-Schmidt orthonormal
for N, lost in ((1, 2), (2, 0)):
    Xs = torch.randn(16, N, dtype=dt)
    M = torch.stack([(B @ Xs).reshape(-1) for B in basis], dim=1)               # w -> (T x_1, ..., T x_N)
    assert torch.linalg.matrix_rank(M).item() == 15 - lost
```

</details>

**Suite :** [Tutoriel 13. Algèbres de Jordan : pourquoi la mécanique quantique utilise des *-algèbres](13_jordan_algebras.md).
