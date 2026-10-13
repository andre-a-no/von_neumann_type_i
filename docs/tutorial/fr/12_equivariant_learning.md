[English](../12_equivariant_learning.md) · [Español](../es/12_equivariant_learning.md) · **Français** · [Deutsch](../de/12_equivariant_learning.md) · [中文](../zh/12_equivariant_learning.md) · [日本語](../ja/12_equivariant_learning.md) · [Русский](../ru/12_equivariant_learning.md)

# 12. Apprentissage automatique équivariant

Une couche d'un réseau de neurones est **équivariante** si transformer l'entrée par une symétrie puis appliquer la
couche donne le même résultat qu'appliquer la couche puis transformer la sortie : renuméroter les sommets d'un
graphe renumérote ses caractéristiques de la même façon. Pour des couches linéaires, c'est une relation de
commutation, de sorte que les couches équivariantes forment un commutant, qui d'après le chapitre 2 est, à un
changement de base près, l'algèbre $\bigoplus_cM_{m_c}(\mathbb C)$. Ce chapitre compte leurs paramètres, les
construit et les entraîne avec `torch_vn_algebra.decompose`, et se sert de la structure d'algèbre pour obtenir des
couches unitaires, positives et spectrales qui sont équivariantes par construction.

**Vous apprendrez** pourquoi les applications linéaires équivariantes entre deux représentations sont des matrices
par blocs dont les tailles sont données par les multiplicités ; comment les compter par les orbites ou par les
caractères (DeepSets et les nombres de Bell des réseaux de graphes) ; pourquoi la moyenne sur le groupe est la
projection orthogonale sur les applications équivariantes ; comment le calcul fonctionnel et les paramétrages
unitaires agissent bloc par bloc sur le commutant ; et comment on compte les vecteurs invariants et les formes
quadratiques invariantes.

Dans tout ce chapitre, $G$ est un groupe fini, et une **représentation** de $G$ sur $\mathbb C^n$ est un
homomorphisme $\pi:G\to M_n(\mathbb C)$ à valeurs unitaires (ce n'est pas une restriction : des $\pi(g)$ inversibles
sont unitaires pour le produit scalaire $\sum_g\langle\pi(g)x,\pi(g)y\rangle$, section 4.1). Son **caractère** est
$\chi_\pi(g)=\operatorname{Tr}\pi(g)$. Les représentations irréductibles et équivalentes et la **multiplicité**
$m_\rho(\pi)$ d'une $\rho$ irréductible dans $\pi$ ($0$ si $\rho$ n'apparaît pas) sont comme à la section 2.4, qui
démontre que $m_\rho(\pi)$ ne dépend pas de la décomposition ainsi que le lemme de Schur ; le
[chapitre 4](04_symmetries.md) traite les deux en détail (sections 4.1–4.2). Sur les matrices rectangulaires,
$\langle S,T\rangle=\operatorname{Tr}(S^*T)$ est le produit scalaire de **Hilbert–Schmidt**. Dans le code, un lot de
vecteurs est un tenseur de forme `(batch, n)`, sur lequel une matrice $T$ agit par `x @ T.T`.

## 12.1 Applications linéaires équivariantes

**Définition.** Soient $\pi_1,\pi_2$ des représentations de $G$ sur $\mathbb C^{n_1}$, $\mathbb C^{n_2}$. Une matrice
$T\in M_{n_2\times n_1}(\mathbb C)$ est **équivariante** (un **opérateur d'entrelacement**) si $T\pi_1(g)=\pi_2(g)T$
pour tout $g\in G$. Ces matrices forment un espace vectoriel $\operatorname{Hom}_G(\pi_1,\pi_2)$ ; pour
$\pi_1=\pi_2=\pi$, c'est le commutant, $\operatorname{Hom}_G(\pi,\pi)=\pi(G)'$.

Les représentations les plus courantes en apprentissage automatique permutent les coordonnées. Supposons que $G$ agisse
sur un ensemble fini $\Omega$ ; la **représentation de permutation** sur $\mathbb C^\Omega$ est $\pi(g)e_x=e_{gx}$.
Exemples : $S_n$ permutant les $n$ points d'un ensemble ($\Omega=[n]=\{1,\dots,n\}$), et $S_n$ renumérotant les
sommets d'un graphe, qui agit sur $\Omega=[n]^2$, c'est-à-dire sur la matrice d'adjacence par
$X\mapsto P_\sigma XP_\sigma^T$.

**Lemme 12.1 (orbites).** Soient $\pi_i$ les représentations de permutation de $G$ sur $\mathbb C^{\Omega_i}$. Alors
$T\in\operatorname{Hom}_G(\pi_1,\pi_2)$ si et seulement si $T_{gy,gx}=T_{yx}$ pour tous $g$, $y\in\Omega_2$,
$x\in\Omega_1$. Les matrices indicatrices des orbites de $G$ sur $\Omega_2\times\Omega_1$ forment une base, et
$\dim\operatorname{Hom}_G(\pi_1,\pi_2)$ est le nombre de ces orbites.

*Démonstration.* $\pi_2(g)T\pi_1(g)^*e_x=\pi_2(g)Te_{g^{-1}x}=\sum_yT_{y,g^{-1}x}e_{gy}$, donc le coefficient $(gy,x)$
de $\pi_2(g)T\pi_1(g)^*$ est $T_{y,g^{-1}x}$. Par conséquent, $T$ est équivariante, c'est-à-dire
$\pi_2(g)T\pi_1(g)^*=T$ pour tout $g$, si et seulement si $T_{gy,x}=T_{y,g^{-1}x}$ pour tous $g,y,x$, ce qui (en
remplaçant $x$ par $gx$) est la condition de l'énoncé : $T$ est constante sur les orbites. Les matrices indicatrices
des orbites ont des supports disjoints, donc elles sont linéairement indépendantes, et toute fonction constante sur
les orbites en est une combinaison. $\square$

**Corollaire 12.2.** (a) (*DeepSets.*) Pour $S_n$ sur $\mathbb C^n$, $n\ge2$, le commutant est
$\operatorname{span}\{1,J\}$, où $J$ est la matrice dont tous les coefficients valent un, et
$\mathbb C^n=\mathbb C\mathbf 1\oplus\mathbf 1^\perp$ est la somme de deux représentations irréductibles non
équivalentes, la triviale et la représentation **standard**, de dimension $n-1$. (b) (*Réseaux de graphes.*) Les
applications équivariantes de $\mathbb C^{[n]^p}$ dans $\mathbb C^{[n]^q}$ (tenseurs d'ordre $p$ et $q$ sur $n$
sommets) sont comptées par les partitions d'un ensemble à $p+q$ éléments en au plus $n$ blocs ; pour $n\ge p+q$, leur
nombre est le **nombre de Bell** $B_{p+q}$ ($B_2=2$, $B_3=5$, $B_4=15$).

*Démonstration.* (b) Deux uplets $x,y\in[n]^{p+q}$ sont dans la même orbite si et seulement s'ils ont le même
**motif d'égalités**, la partition des positions $\{1,\dots,p+q\}$ définie par $i\sim j\Leftrightarrow x_i=x_j$. Le
motif est clairement invariant. Réciproquement, si $x,y$ ont le même motif, $x_i\mapsto y_i$ est une injection bien
définie de $\{x_i\}$ dans $[n]$, qui se prolonge en une permutation $\sigma$ telle que $\sigma x=y$. Toute partition
en au plus $n$ blocs est un motif (donner aux blocs des valeurs distinctes). Le lemme 12.1 permet de conclure.
(a) Pour $p=q=1$, il y a deux orbites sur $[n]^2$, la diagonale et son complémentaire, d'indicatrices $1$ et $J-1$.
Donc $\dim\pi(S_n)'=2=\sum_cm_c^2$ (théorème 2.6, proposition 2.7), ce qui impose deux secteurs avec $m_c=1$,
c'est-à-dire deux facteurs irréductibles non équivalents, et $\pi(S_n)'$ est engendré par leurs projections
$P_1,P_2$. La projection sur la droite invariante $\mathbb C\mathbf 1$ appartient à $\pi(S_n)'$ (lemme 2.1), donc c'est
$P_1$ ou $P_2$ : la droite, qui porte la représentation triviale, est un secteur, et $\mathbf 1^\perp$ est
l'autre. $\square$

La couche DeepSets $x\mapsto\lambda x+\gamma\,(\mathbf 1^Tx)\mathbf 1$ (Zaheer et al. 2017) est (a) ; les bases de
Maron et al. (2019) pour les réseaux de graphes sont (b). Le code compte les orbites par force brute, les compare au
commutant calculé numériquement et décompose les deux représentations de $S_4$ :

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

Ainsi $\mathbb C^{16}$ contient la représentation triviale deux fois, une représentation de dimension $2$ une fois, et
deux représentations de dimension $3$, l'une une fois et l'autre trois fois : $16=1\cdot2+2\cdot1+3\cdot1+3\cdot3$ et
$15=2^2+1^2+1^2+3^2$. Pour deux représentations, tout se ramène à un seul commutant :

**Lemme 12.3 (somme directe).** Soit $\pi=\pi_1\oplus\pi_2$ sur $\mathbb C^{n_1}\oplus\mathbb C^{n_2}$, avec les
inclusions isométriques $Q_i$ et les projections $P_i=Q_iQ_i^*$. Alors $P_1,P_2\in\pi(G)'$,
$\operatorname{Hom}_G(\pi_1,\pi_2)=Q_2^*\,\pi(G)'\,Q_1$, et $T\mapsto Q_2TQ_1^*$ est une bijection linéaire de
$\operatorname{Hom}_G(\pi_1,\pi_2)$ sur $P_2\,\pi(G)'P_1$.

*Démonstration.* $\pi(g)Q_i=Q_i\pi_i(g)$ et $Q_i^*\pi(g)=\pi_i(g)Q_i^*$, car $\pi(g)$ est diagonale par blocs ; donc
$P_i\pi(g)=Q_i\pi_i(g)Q_i^*=\pi(g)P_i$. Si $X\in\pi(G)'$, alors
$Q_2^*XQ_1\pi_1(g)=Q_2^*X\pi(g)Q_1=Q_2^*\pi(g)XQ_1=\pi_2(g)Q_2^*XQ_1$. Si $T$ est équivariante, alors
$\pi(g)Q_2TQ_1^*=Q_2\pi_2(g)TQ_1^*=Q_2T\pi_1(g)Q_1^*=Q_2TQ_1^*\pi(g)$, donc $Q_2TQ_1^*\in\pi(G)'$, elle est égale à
$P_2(Q_2TQ_1^*)P_1$, et $Q_2^*(Q_2TQ_1^*)Q_1=T$. Réciproquement, tout $Y=P_2XP_1$ avec $X\in\pi(G)'$ appartient à
$\pi(G)'$ et $Y=Q_2(Q_2^*YQ_1)Q_1^*$. $\square$

**Théorème 12.4 (les applications équivariantes sont des matrices par blocs).** Soient $\rho_1,\dots,\rho_C$ les
représentations irréductibles deux à deux non équivalentes qui apparaissent dans $\pi_1\oplus\pi_2$,
$k_c=\dim\rho_c$, $a_c=m_{\rho_c}(\pi_1)$, $b_c=m_{\rho_c}(\pi_2)$. Il existe des unitaires
$V_1:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{a_c}\to\mathbb C^{n_1}$ et
$V_2:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{b_c}\to\mathbb C^{n_2}$ avec
$V_i^*\pi_i(g)V_i=\bigoplus_c\rho_c(g)\otimes1$, tels que
$$
\operatorname{Hom}_G(\pi_1,\pi_2)=\Big\{V_2\Big(\bigoplus_{c=1}^C1_{k_c}\otimes y_c\Big)V_1^*:\ y_c\in M_{b_c\times a_c}(\mathbb C)\Big\},
$$
et $(y_c)_c\mapsto T$ est injective. Par conséquent $\operatorname{Hom}_G(\pi_1,\pi_2)\cong\bigoplus_cM_{b_c\times a_c}(\mathbb C)$,
$$
\dim\operatorname{Hom}_G(\pi_1,\pi_2)=\sum_\rho m_\rho(\pi_1)\,m_\rho(\pi_2),
$$
la somme portant sur toutes les représentations irréductibles à équivalence près ; pour $\pi_1=\pi_2$, on retrouve
$\pi(G)'\cong\bigoplus_cM_{m_c}$ (proposition 2.7), et $\operatorname{Hom}_G(\pi_1,\pi_2)=0$ si et seulement si aucune
représentation irréductible n'apparaît dans les deux.

*Démonstration.* Appliquons le théorème 2.6 et la proposition 2.7 à $\pi=\pi_1\oplus\pi_2$ : il existe un unitaire
$V$ tel que $V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$ et $V^*\pi(G)'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$. Comme
$X\mapsto V^*XV$ est un $*$-isomorphisme et que $P_1\in\pi(G)'$ (lemme 12.3), $V^*P_1V=\bigoplus_c1\otimes p_c$ avec
des projections $p_c$. Choisissons des unitaires $u_c$ tels que $u_c^*p_cu_c=1_{a_c}\oplus0_{b_c}$ (théorème 1.3), où
désormais $a_c=\operatorname{rank}p_c$, $b_c=m_c-a_c$, et remplaçons $V$ par $V\bigoplus_c(1\otimes u_c)$ ; cela ne
modifie aucune des deux formes normales, car $1\otimes u_c$ commute avec $\rho_c(g)\otimes1$. Soient
$\iota_c:\mathbb C^{a_c}\to\mathbb C^{m_c}$ et $\kappa_c:\mathbb C^{b_c}\to\mathbb C^{m_c}$ les inclusions des $a_c$
premières et des $b_c$ dernières coordonnées, et
$$
V_1=Q_1^*V\bigoplus_c(1\otimes\iota_c),\qquad V_2=Q_2^*V\bigoplus_c(1\otimes\kappa_c).
$$
$W_1=V\bigoplus_c(1\otimes\iota_c)$ est une isométrie d'image $\operatorname{range}P_1=\operatorname{range}Q_1$, donc
$V_1^*V_1=W_1^*P_1W_1=1$, et $V_1$ est unitaire puisque $\sum_ck_ca_c=\operatorname{rank}P_1=n_1$. De plus
$\pi_1(g)V_1=Q_1^*\pi(g)W_1=Q_1^*V\big(\bigoplus_c\rho_c(g)\otimes1_{m_c}\big)\bigoplus_c(1\otimes\iota_c)
=V_1\bigoplus_c\rho_c(g)\otimes1_{a_c}$. Ainsi $\pi_1\cong\bigoplus_c\rho_c\otimes1_{a_c}$ avec des $\rho_c$
irréductibles non équivalentes, et $a_c=m_{\rho_c}(\pi_1)$ par l'unicité des multiplicités (section 2.4) ; de même
pour $V_2$ et $b_c$. Enfin, pour $X=V\big(\bigoplus_c1\otimes x_c\big)V^*\in\pi(G)'$,
$P_2XP_1=V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*$, où $y_c=\kappa_c^*x_c\iota_c$ est le coin inférieur
gauche de taille $b_c\times a_c$ de $x_c$, une matrice arbitraire. D'après le lemme 12.3,
$\operatorname{Hom}_G(\pi_1,\pi_2)=Q_2^*P_2\,\pi(G)'P_1Q_1$ est formé des matrices
$Q_2^*V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*Q_1=V_2\big(\bigoplus_c1\otimes y_c\big)V_1^*$, et
$(y_c)\mapsto T$ est injective puisque $V_1,V_2$ sont unitaires. Une $\rho$ irréductible qui n'apparaît pas dans
$\pi_1\oplus\pi_2$ vérifie $m_\rho(\pi_1)m_\rho(\pi_2)=0$. $\square$

Pour le corollaire 12.2(b) avec $n=4$ : les caractéristiques de sommets $\pi_1$ ont $a=(1,0,0,1)$ dans les quatre
secteurs de $\pi_1\oplus\pi_2$ (la triviale, celle de dimension $2$, les deux de dimension $3$), les
caractéristiques d'arêtes $\pi_2$ ont $b=(2,1,1,3)$, donc $\dim\operatorname{Hom}_G=1\cdot2+1\cdot3=5=B_3$.

> **Interprétation physique.** Le théorème 12.4 est le noyau d'algèbre linéaire des **règles de sélection** et le
> cas scalaire du théorème de Wigner–Eckart : un opérateur qui commute avec la symétrie ne peut relier que des états
> de même type de symétrie, et à l'intérieur d'un type il agit comme $1_{k_c}\otimes y_c$, c'est-à-dire par un
> « élément de matrice réduit » $y_c$ aveugle à la position à l'intérieur du multiplet. (Les opérateurs qui se
> transforment eux-mêmes selon une représentation non triviale, comme les opérateurs vectoriels, requièrent en outre
> les coefficients de Clebsch–Gordan.) Un hamiltonien symétrique est le cas $\pi_1=\pi_2$ ; ses dégénérescences $k_c$
> sont les dimensions des représentations irréductibles.

## 12.2 Une couche équivariante avec la bibliothèque

La démonstration du théorème 12.4 est un algorithme, dès que l'on connaît la décomposition du théorème 2.6 ;
`dc.from_group` la calcule par la méthode randomisée de la section 4.3. Pour $\pi_1=\pi_2$, c'est simplement le
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

La seconde façon de rendre une couche équivariante est de la moyenner sur le groupe. Le **twirl** est
$$
\mathcal P(W)=\frac1{|G|}\sum_{g\in G}\pi_2(g)\,W\,\pi_1(g)^*,\qquad W\in M_{n_2\times n_1}(\mathbb C).
$$

**Théorème 12.5 (le twirl est la projection orthogonale).** $\mathcal P$ est la projection orthogonale, pour le
produit scalaire de Hilbert–Schmidt, de $M_{n_2\times n_1}(\mathbb C)$ sur $\operatorname{Hom}_G(\pi_1,\pi_2)$, et
$$
\dim\operatorname{Hom}_G(\pi_1,\pi_2)=\frac1{|G|}\sum_{g\in G}\overline{\chi_{\pi_1}(g)}\,\chi_{\pi_2}(g).
$$

*Démonstration.* Pour $h\in G$, $\pi_2(h)\mathcal P(W)\pi_1(h)^*=\frac1{|G|}\sum_g\pi_2(hg)W\pi_1(hg)^*=\mathcal P(W)$,
car $g\mapsto hg$ est une bijection de $G$ ; donc $\mathcal P(W)$ est équivariante. Si $W$ est équivariante, chaque
terme est égal à $W$, donc $\mathcal P(W)=W$. Par conséquent $\mathcal P^2=\mathcal P$, d'image
$\operatorname{Hom}_G(\pi_1,\pi_2)$. Elle est autoadjointe :
$\langle S,\mathcal P(W)\rangle=\frac1{|G|}\sum_g\operatorname{Tr}\big((\pi_2(g)^*S\pi_1(g))^*W\big)$
$=\langle\mathcal P(S),W\rangle$, où l'on a substitué $g\mapsto g^{-1}$ et utilisé $\pi_i(g^{-1})=\pi_i(g)^*$. Un
idempotent autoadjoint est la projection orthogonale sur son image. Sa trace est son rang. L'application linéaire
$W\mapsto AWB$ sur $M_{n_2\times n_1}$ a pour trace $\operatorname{Tr}A\operatorname{Tr}B$ (le coefficient de l'unité
matricielle $E_{ij}$ dans $AE_{ij}B$ est $A_{ii}B_{jj}$) ; avec $A=\pi_2(g)$, $B=\pi_1(g)^*$, cela donne la
formule. $\square$

Avec le théorème 12.4,
$\sum_\rho m_\rho(\pi_1)m_\rho(\pi_2)=|G|^{-1}\sum_g\overline{\chi_{\pi_1}(g)}\chi_{\pi_2}(g)$ ; pour $\pi_1=\pi_2$, c'est
la formule des caractères pour $\dim\pi(G)'$ de la section 4.2, que calcule `dc.commutant_dimension`. Le paramétrage
$W\mapsto\mathcal P(W)$ atteint toute application équivariante, $\mathcal P(W)$ est l'application équivariante la plus
proche de $W$, et, comme $\mathcal P$ est autoadjointe, le gradient par rapport à $W$ est l'image par $\mathcal P$ du
gradient par rapport à $T=\mathcal P(W)$. Mais $n_1n_2-\dim\operatorname{Hom}_G$ paramètres sont gaspillés dans le
noyau ($256-15$ pour les caractéristiques d'arêtes), et une évaluation coûte $O(|G|n_1n_2(n_1+n_2))$ pour des $\pi(g)$
denses, avec $|G|=n!$ pour $S_n$. La bibliothèque calcule la même projection à partir de la décomposition,
`D.conditional_expectation(W, 'commutant')`, avec quatre produits $n\times n$.

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

Les couches ci-dessus sont linéaires. Un réseau les alterne avec des non-linéarités, qui doivent elles aussi être
équivariantes. Pour une représentation de permutation, une fonction appliquée à chaque coordonnée est équivariante,
puisque $(\pi(g)x)_{gy}=x_y$ ; pour d'autres représentations, elle ne l'est pas, et l'on utilise des non-linéarités
qui agissent sur des normes ou sur les coefficients dans une base que $G$ permute (exercice 3).

> **Interprétation physique.** Le twirl est l'espérance conditionnelle sur le commutant (exercice 3 du chapitre 2) :
> moyenner un couplage sur toutes les opérations de symétrie en conserve exactement la partie symétrique.
> « Symétriser » un hamiltonien ou une interaction en sommant sur le groupe de symétrie est la même projection.

## 12.3 Couches contraintes et spectrales

Les poids $w=(w_c)_c$ de `EquivariantLinear` vivent dans une $*$-algèbre, de sorte que toutes les constructions des
chapitres 1 et 11 s'y appliquent. Notons $\Phi(w)=V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*$ pour
$w\in\bigoplus_cM_{m_c}$.

**Proposition 12.6.** (a) $\Phi$ est un $*$-isomorphisme de $\bigoplus_cM_{m_c}$ sur $\pi(G)'$, et $\Phi(w)$ est
autoadjoint, unitaire ou positif si et seulement si chaque $w_c$ l'est. (b) Pour $w$ normal et une fonction $f$ sur
$\bigcup_c\operatorname{spec}w_c$, $f(\Phi(w))=\Phi\big((f(w_c))_c\big)$ ; en particulier $f(X)\in\pi(G)'$ pour tout
$X\in\pi(G)'$ normal. (c) Les unitaires équivariants sont les $\Phi(u_0e^{s-s^*})$, $s$ parcourant
$\bigoplus_cM_{m_c}$, pour n'importe quel unitaire fixé $u_0\in\bigoplus_cM_{m_c}$ ; en particulier, ce sont les
exponentielles $e^S$ des $S\in\pi(G)'$ antiautoadjoints.

*Démonstration.* (a) $\Phi$ est linéaire, unifère, multiplicative et compatible avec l'involution, car
$(1\otimes x)(1\otimes y)=1\otimes xy$ et $(1\otimes x)^*=1\otimes x^*$ ; elle est injective, et surjective sur
$\pi(G)'$ d'après le théorème 2.6. Donc $\Phi(w)^*=\Phi(w)$ si et seulement si $w^*=w$, et
$\Phi(w)^*\Phi(w)=1=\Phi(w)\Phi(w)^*$ si et seulement si $w^*w=1=ww^*$. Si tous les $w_c=t_c^*t_c$, alors
$\Phi(w)=\Phi(t)^*\Phi(t)\ge0$ (corollaire 1.7). Réciproquement, si $\Phi(w)\ge0$, alors pour $v\in\mathbb C^{m_c}$ et
le vecteur $\xi=V(e_1\otimes v)$ du secteur $c$, $0\le\langle\xi,\Phi(w)\xi\rangle=\langle v,w_cv\rangle$.
(b) Soit $w_c=u_c\Lambda_cu_c^*$ (théorème 1.3). Alors $\Phi(w)=U\big(\bigoplus_c1\otimes\Lambda_c\big)U^*$ avec
l'unitaire $U=V\bigoplus_c(1\otimes u_c)$, une diagonalisation où les valeurs propres de $w_c$ sont répétées $k_c$
fois ; par définition et d'après le théorème 1.5(a),
$f(\Phi(w))=U\big(\bigoplus_c1\otimes f(\Lambda_c)\big)U^*=\Phi\big((f(w_c))_c\big)$.
(c) D'après (a), les unitaires équivariants sont les $\Phi(u)$ avec $u$ unitaire dans $\bigoplus_cM_{m_c}$, et ce sont
les $u_0e^{s-s^*}$ d'après le lemme 11.1. Avec $u_0=1$, (b) donne $\Phi(e^{s-s^*})=e^{\Phi(s-s^*)}$, et $\Phi(s-s^*)$
parcourt les éléments antiautoadjoints de $\pi(G)'$. $\square$

Ainsi une couche équivariante orthogonale (unitaire) est `optimize.UnitaryParam` sur `D.commutant_algebra()` suivi de
$\Phi$ ; d'après (c), ce paramétrage atteint tout unitaire équivariant, et chaque itéré d'une méthode de gradient est
exactement unitaire et exactement équivariant. De même, `PositiveParam` fournit des couches équivariantes positives,
normalisées comme au chapitre 11 (spectre dans $[\ell,1]$ de maximum $1$ ; d'après (a), une couche positive générale
est $\Phi(t^*t)$), et une couche spectrale $f(L)$ d'un générateur équivariant autoadjoint $L$, comme la diffusion
$e^{-tL}$, demande des diagonalisations de coût $\sum_cm_c^3$ (ici $37$) au lieu de $n^3=4096$, avec des gradients
donnés par la formule de Daleckii–Krein (`apply_function`). Le code ajuste un unitaire équivariant à $32$ couples
entrée–sortie, et confronte (b) et son gradient par rapport à $t$, $\frac d{dt}e^{-tL}=-Le^{-tL}$, aux matrices
denses.

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

> **Interprétation physique.** Un hamiltonien symétrique $H\in\pi(G)'$ engendre une évolution temporelle symétrique
> $e^{-itH}\in\pi(G)'$ (chapitre 8), qui ne mélange jamais les secteurs et n'agit que sur les espaces de
> multiplicité. Une couche unitaire équivariante est une telle évolution, apprise au lieu d'être donnée.

## 12.4 Invariants

Un vecteur $x$ est **invariant** si $\pi(g)x=x$ pour tout $g$ ; une fonction $q$ sur $\mathbb C^n$ est invariante si
$q(\pi(g)x)=q(x)$. Les caractéristiques invariantes sont ce qu'un réseau invariant (un classifieur de graphes, par
exemple) peut renvoyer.

**Proposition 12.7.** (a) Les vecteurs invariants forment l'image de la projection orthogonale
$P_0=|G|^{-1}\sum_g\pi(g)$ ; leur dimension est $m_{\mathrm{triv}}(\pi)=|G|^{-1}\sum_g\chi_\pi(g)$, et si la
représentation triviale apparaît, en tant que $\rho_{c_0}$, alors $P_0$ est la projection centrale minimale $P_{c_0}$
de $\pi(G)''$.
(b) Une forme hermitienne $q_A(x)=\langle x,Ax\rangle$, $A=A^*$, est invariante si et seulement si $A\in\pi(G)'$ ; les
formes hermitiennes invariantes forment un espace vectoriel réel de dimension $\sum_cm_c^2$.
(c) Les polynômes invariants homogènes de degré $2$ en les coordonnées $x_1,\dots,x_n$ forment un espace de dimension
$$
\frac1{|G|}\sum_{g\in G}\frac{\chi_\pi(g)^2+\chi_\pi(g^2)}2 .
$$
Pour une représentation de permutation, c'est le nombre d'orbites de $G$ sur les multiensembles $\{x,y\}$ de deux
éléments de $\Omega$.

*Démonstration.* (a) Les vecteurs invariants sont $\operatorname{Hom}_G(\mathrm{triv},\pi)$, où $\mathrm{triv}$ est la
représentation triviale sur $\mathbb C=M_{1\times1}$ ; pour celle-ci, le twirl est $\mathcal P(x)=P_0x$, de sorte que
le théorème 12.5 donne la projection et la dimension, et le théorème 12.4 donne $m_{\mathrm{triv}}(\pi)$. Dans la base
du théorème 2.6, écrivons un vecteur invariant $\xi=\sum_c\xi_c$ avec
$\xi_c=\sum_jv_{c,j}\otimes e_j\in\mathbb C^{k_c}\otimes\mathbb C^{m_c}$. Comme $\pi(g)$ préserve les secteurs,
$\rho_c(g)v_{c,j}=v_{c,j}$ pour tout $g$. Un $v_{c,j}$ non nul engendre une droite invariante qui porte la
représentation triviale ; alors $\rho_c$, étant irréductible, est triviale et $c=c_0$. Les vecteurs invariants sont
donc les vecteurs du secteur $c_0$, où $\pi(g)$ agit comme $1$ : l'image de $P_{c_0}$, qui est par conséquent $P_0$.
(b) $q_A(\pi(g)x)=\langle x,\pi(g)^*A\pi(g)x\rangle$, et une matrice $D$ telle que $\langle x,Dx\rangle=0$ pour tout
$x$ est nulle (polarisation, démonstration du corollaire 1.4). Donc $q_A$ est invariante si et seulement si
$\pi(g)^*A\pi(g)=A$ pour tout $g$, c'est-à-dire $A\in\pi(G)'$. D'après la proposition 12.6(a), les éléments
autoadjoints de $\pi(G)'$ sont les $\Phi(w)$ avec des $w_c$ hermitiennes, un espace réel de dimension $\sum_cm_c^2$.
(c) Un polynôme quadratique $p(x)=x^TBx$ est déterminé par sa matrice symétrique $B$. Posons
$(g\cdot p)(x)=p(\pi(g)^{-1}x)$ ; alors $g\cdot(h\cdot p)=(gh)\cdot p$, et $g\cdot p$ a pour matrice
$\sigma(g)(B)=A^TBA$, $A=\pi(g)^{-1}$. La moyenne $M=|G|^{-1}\sum_g\sigma(g)$ vérifie $\sigma(h)M=M$ et $MB=B$ pour $B$
invariante, donc (comme dans la démonstration du théorème 12.5) elle est idempotente, d'image les polynômes
invariants, et la dimension de ceux-ci est $\operatorname{Tr}M=|G|^{-1}\sum_g\operatorname{Tr}\sigma(g)$. $A$ est
unitaire, donc $A=R\Lambda R^{-1}$ avec $\Lambda=\operatorname{diag}(\lambda)$ (théorème 1.3). Les matrices symétriques
$B_{ij}=(R^{-1})^T(e_ie_j^T+e_je_i^T)R^{-1}$, $i\le j$, forment une base, et
$A^TB_{ij}A=(R^{-1})^T\Lambda(e_ie_j^T+e_je_i^T)\Lambda R^{-1}=\lambda_i\lambda_jB_{ij}$. Donc
$\operatorname{Tr}\sigma(g)=\sum_{i\le j}\lambda_i\lambda_j=\frac12\big((\operatorname{Tr}A)^2+\operatorname{Tr}A^2\big)
=\frac12\big(\chi_\pi(g^{-1})^2+\chi_\pi(g^{-2})\big)$, et sommer sur $g^{-1}$ au lieu de $g$ donne la formule. Pour une
représentation de permutation, $G$ permute les monômes $x_ax_b$, donc un polynôme est invariant si et seulement si ses
coefficients sont constants sur les orbites de monômes. $\square$

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

Sur $\mathbb C^4$, ce sont $\sum_ix_i^2$ et $\sum_{i<j}x_ix_j$. Tous les degrés sont comptés par la formule de Molien
$\sum_d\dim\mathbb C[x]^G_dt^d=|G|^{-1}\sum_g\det(1-t\pi(g))^{-1}$, dont le coefficient de $t^2$ est (c) (énoncée sans
démonstration ; B. Sturmfels, *Algorithms in Invariant Theory*, Springer 2008, théorème 2.2.1).

> **Interprétation physique.** Dans la théorie de Landau des transitions de phase, l'énergie libre est un polynôme
> invariant en un paramètre d'ordre réel $x$, sur lequel $G$ agit par des matrices orthogonales réelles. Sa partie
> quadratique est $x^TBx$ avec $B$ dans le commutant (car $\pi(g)^T=\pi(g)^{-1}$), donc, d'après la
> proposition 12.6(a), elle est diagonale par blocs selon les secteurs : la phase symétrique devient instable dans le
> secteur où une valeur propre de $B$ change de signe en premier, et (génériquement) la représentation irréductible de
> ce secteur détermine la symétrie de la phase ordonnée.

## 12.5 Points de comparaison

Des bibliothèques matures fournissent des couches équivariantes pour des groupes particuliers. **escnn** (Cesa, Lang
et Weiler 2022, successeur de e2cnn) construit des convolutions orientables (*steerable*) pour des sous-groupes de
$O(2)$ et $O(3)$ à partir de représentations irréductibles et de bases harmoniques de noyaux ; **e3nn** (Geiger et
Smidt 2022) construit des réseaux $E(3)$-équivariants à partir d'harmoniques sphériques et de produits tensoriels de
Clebsch–Gordan. **EMLP** (Finzi, Welling et Wilson 2021) traite n'importe quel groupe de matrices : il impose
$\pi_2(g)T=T\pi_1(g)$ pour des générateurs discrets et $d\pi_2(X)T=Td\pi_1(X)$ pour des générateurs de l'algèbre de
Lie, et calcule une base des solutions comme un noyau numérique — le problème linéaire de `dc.commutant_basis`. Les
convolutions de groupe (Cohen et Welling 2016) sont le cas de la représentation régulière, et Kondor et Trivedi (2018)
montrent que, pour des cartes de caractéristiques sur des espaces homogènes d'un groupe compact, les couches
linéaires équivariantes sont des convolutions généralisées. Pour les permutations, DeepSets (Zaheer et al. 2017) et
les bases de Maron et al. (2019) sont le corollaire 12.2.

La bibliothèque ne fournit ni convolutions sur des grilles, ni noyaux orientables sur des espaces continus, ni
décompositions sur $\mathbb R$ (types réel, complexe et quaternionique ; l'exercice 2 traite les matrices réelles), et
ne vise pas les grands $n$ : $V$ est une matrice dense $n\times n$ et la décomposition coûte $O(|G|n^3)$ par
échantillon aléatoire pour un groupe donné par tous ses éléments (les générateurs d'algèbres de Lie passent par
`dc.from_generators` ; section 4.3). Ce qu'elle apporte, c'est la **structure** de l'espace des solutions : une base du
commutant est une liste de matrices, tandis que la décomposition identifie le commutant à la $*$-algèbre
$\bigoplus_cM_{m_c}$. Les produits, les adjoints, l'unitarité, la positivité et le calcul fonctionnel agissent alors
bloc par bloc sur de petites matrices, avec des gradients (proposition 12.6), et `to_commutant`/`from_commutant`
font la traduction entre les deux points de vue.

## Exercices

1. (**DeepSets avec canaux.**) $S_n$ ($n\ge2$) agit sur $\mathbb C^n\otimes\mathbb C^d$ par $P_\sigma\otimes1_d$.
   Montrer que les applications équivariantes de $\mathbb C^n\otimes\mathbb C^{d_1}$ dans
   $\mathbb C^n\otimes\mathbb C^{d_2}$ sont les $1\otimes A+J\otimes B$ avec $A,B\in M_{d_2\times d_1}$, et le
   vérifier avec `EquivariantMap`.
2. (**Couches réelles.**) Supposons que $\pi_1,\pi_2$ soient à matrices réelles (par exemple des représentations de
   permutation). Montrer que les matrices équivariantes réelles forment un espace vectoriel réel de dimension
   $\dim_{\mathbb C}\operatorname{Hom}_G(\pi_1,\pi_2)$, et que $T\mapsto\operatorname{Re}T$ envoie
   $\operatorname{Hom}_G(\pi_1,\pi_2)$ sur cet espace.
3. (**Non-linéarités.**) (a) Montrer qu'appliquer une fonction $f:\mathbb C\to\mathbb C$ à chaque coordonnée commute
   avec toute représentation de permutation. (b) Montrer que cela ne commute pas avec la rotation d'angle $\pi/2$ de
   $\mathbb R^2$ pour $f(t)=\max(t,0)$. (c) Montrer que $x\mapsto\varphi(\|x\|)x$ commute avec toute représentation
   unitaire.
4. (**Combien d'échantillons déterminent une couche ?**) Soient $T=\Phi(w)\in\pi(G)'$ et
   $x_1,\dots,x_N\in\mathbb C^n$. Écrivons $V^*x_i$ dans le secteur $c$ sous la forme
   $\sum_{a,j}(X_{i,c})_{aj}e_a\otimes e_j$ avec $X_{i,c}\in M_{k_c\times m_c}$, et soit $S_c\in M_{Nk_c\times m_c}$ la
   matrice obtenue en empilant $X_{1,c},\dots,X_{N,c}$. Montrer que les $T$ équivariantes telles que
   $Tx_1=\dots=Tx_N=0$ forment un espace de dimension $\sum_cm_c(m_c-\operatorname{rank}S_c)$ ; pour des $x_i$
   gaussiens, cette dimension est presque sûrement $\sum_cm_c\max(0,m_c-Nk_c)$. Pour la représentation des arêtes,
   deux échantillons déterminent une couche équivariante, alors qu'une matrice $16\times16$ quelconque en demande
   seize.

<details>
<summary>Solutions</summary>

1. $P_\sigma J=J=JP_\sigma$, donc ces applications sont équivariantes. Elles sont linéairement indépendantes : si
   $1\otimes A+J\otimes B=0$, les blocs hors diagonale donnent $B=0$, puis $A=0$. D'après le corollaire 12.2(a) et le
   théorème 2.6, $P_\sigma\otimes1_d\cong(\mathrm{triv}\otimes1_d)\oplus(\mathrm{std}\otimes1_d)$, donc les deux
   représentations irréductibles ont multiplicité $d$, et le théorème 12.4 donne $\dim\operatorname{Hom}_G=2d_1d_2$,
   la dimension de l'espace des $1\otimes A+J\otimes B$.

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

2. Si $T$ est équivariante, $\bar T$ l'est aussi, puisque $\overline{\pi_i(g)}=\pi_i(g)$. Par conséquent
   $\operatorname{Re}T=(T+\bar T)/2$ et $\operatorname{Im}T=(T-\bar T)/2i$ sont réelles et équivariantes, de sorte que
   l'application $T\mapsto\operatorname{Re}T$ est surjective sur les matrices équivariantes réelles (elle est
   l'identité sur celles-ci), et $\operatorname{Hom}_G(\pi_1,\pi_2)=H_{\mathbb R}\oplus iH_{\mathbb R}$, où
   $H_{\mathbb R}$ désigne les matrices équivariantes réelles. Une base réelle de $H_{\mathbb R}$ est alors une base
   complexe de $\operatorname{Hom}_G$ : elle l'engendre, et une relation complexe $\sum_j(\alpha_j+i\beta_j)T_j=0$ avec
   des $T_j$ réelles se scinde en parties réelle et imaginaire. Dans le lemme 12.1, les indicatrices des orbites
   forment une telle base réelle. (Les blocs $y_c$ du théorème 12.4 sont complexes même pour des données réelles :
   $V$ est complexe en général, et le paramétrage par $\operatorname{Re}T$ utilise deux fois plus de paramètres réels
   que nécessaire.)

3. (a) $(\pi(g)f(x))_{gy}=f(x)_y=f(x_y)=f((\pi(g)x)_{gy})$. (b) Pour $R(x_1,x_2)=(-x_2,x_1)$ et $x=(1,-1)$,
   $f(Rx)=f(1,1)=(1,1)$ mais $Rf(x)=R(1,0)=(0,1)$. (c) $\|\pi(g)x\|=\|x\|$, donc
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

4. $(1\otimes w_c)(e_a\otimes e_j)=\sum_le_a\otimes(w_c)_{lj}e_l$, donc $V^*Tx_i$ a, dans le secteur $c$, la matrice
   de coefficients $X_{i,c}w_c^T$. Ainsi $Tx_i=0$ pour tout $i$ si et seulement si $S_cw_c^T=0$ pour tout $c$,
   c'est-à-dire si et seulement si chaque colonne de $w_c^T$ appartient à $\ker S_c$, un espace de dimension
   $m_c-\operatorname{rank}S_c$ : les $w_c$ ayant cette propriété forment un espace de dimension
   $m_c(m_c-\operatorname{rank}S_c)$, et $\Phi$ est injective. Pour des $x_i$ gaussiens, $V^*x_i$ est gaussien ($V$
   est unitaire), donc $S_c$ a des coefficients gaussiens indépendants et est presque sûrement de rang
   $\min(Nk_c,m_c)$ : pour la sous-matrice carrée supérieure gauche de cette taille, $|\det|^2$ est un polynôme réel
   en les parties réelles et imaginaires des coefficients, non identiquement nul (il vaut $1$ lorsque la sous-matrice
   est l'identité), et un tel polynôme ne s'annule que sur un ensemble de mesure de Lebesgue nulle (démontré à la
   section 4.3). Pour la représentation des arêtes, $(k,m)=(1,2),(2,1),(3,1),(3,3)$, seul le secteur $(1,2)$
   contribue : $2\cdot(2-1)=2$ pour $N=1$, et $0$ pour $N=2$.

```python
basis = D_edge.basis('commutant')                                # (15, 16, 16), Hilbert-Schmidt orthonormal
for N, lost in ((1, 2), (2, 0)):
    Xs = torch.randn(16, N, dtype=dt)
    M = torch.stack([(B @ Xs).reshape(-1) for B in basis], dim=1)               # w -> (T x_1, ..., T x_N)
    assert torch.linalg.matrix_rank(M).item() == 15 - lost
```

</details>

**Pour aller plus loin.** Ce chapitre est le dernier du cours. Le [chapitre 4](04_symmetries.md) développe la théorie
des représentations utilisée ici : le lemme de Schur, les caractères et l'algorithme randomisé sur lequel repose
`dc.from_group` (sections 4.1–4.3). Le [chapitre 9](09_spin_chains.md) applique les mêmes secteurs aux chaînes de spins,
et la mesure de Haar sur $\prod_cU(m_c)$ du [chapitre 5](05_random_matrices.md) fournit l'initialisation aléatoire
naturelle d'une couche unitaire équivariante. Pour les mathématiques, voir J.-P. Serre, *Linear Representations of
Finite Groups* (Springer 1977) et W. Fulton et J. Harris, *Representation Theory* (Springer 1991) ; pour
l'apprentissage équivariant, M. Bronstein, J. Bruna, T. Cohen et P. Veličković, *Geometric Deep Learning* (2021), et
M. Weiler et al., *Equivariant and Coordinate Independent Convolutional Networks* (2023) ; pour l'algorithme de
décomposition, les références de `torch_vn_algebra/decompose.py`.
