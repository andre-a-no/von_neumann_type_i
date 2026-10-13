[English](../07_tensor_products.md) · [Español](../es/07_tensor_products.md) · **Français** · [Deutsch](../de/07_tensor_products.md) · [中文](../zh/07_tensor_products.md) · [日本語](../ja/07_tensor_products.md) · [Русский](../ru/07_tensor_products.md)

# 7. Produits tensoriels et intrication

Deux systèmes quantiques pris ensemble sont décrits par le produit tensoriel de leurs algèbres.
En présence de secteurs, il y a deux produits : le produit simple, où chaque partie conserve sa
propre charge, et le produit *fusionné*, où seule la charge totale est conservée. Ce chapitre
établit la structure des deux, définit la trace partielle comme la duale d'un plongement et
étudie l'intrication.

**Vous apprendrez** comment les secteurs se multiplient et comment la bibliothèque les ordonne,
pourquoi la trace partielle est un canal et, une fois normalisée, une espérance conditionnelle
d'indice fini, quel est l'état réduit d'un état qui conserve la charge, comment l'intrication se
décompose en une partie de nombre et une partie configurationnelle, ainsi que la pureté et
l'entropie moyennes des états aléatoires.

Dans tout ce chapitre, $\mathcal M_1=\bigoplus_{c=1}^{C_1}M_{k_c}(\mathbb C)$ agit sur
$H_1=\bigoplus_c\mathbb C^{k_c}$, avec les projecteurs centraux minimaux $P_c$, et
$\mathcal M_2=\bigoplus_{d=1}^{C_2}M_{l_d}(\mathbb C)$ agit sur $H_2=\bigoplus_d\mathbb C^{l_d}$, avec
les projecteurs centraux minimaux $Q_d$. $\operatorname{Tr}$ désigne la trace de l'espace sur lequel
agit un opérateur ; sur $\mathcal M_i$, c'est la trace brute du chapitre 3. Pour $X$ sur
$\mathbb C^k\otimes\mathbb C^l$, on note $X_{(a,b),(a',b')}=\langle e_a\otimes e_b,X\,e_{a'}\otimes e_{b'}\rangle$,
et l'on identifie $\mathbb C^k\otimes\mathbb C^l$ à $\mathbb C^{kl}$ par $e_a\otimes e_b\mapsto e_{(a-1)l+b}$
(la convention de Kronecker), de sorte que $(A\otimes B)_{(a,b),(a',b')}=A_{aa'}B_{bb'}$.

## 7.1 Le produit tensoriel de deux algèbres

Pour des $*$-algèbres $\mathcal M_1$ sur $H_1$ et $\mathcal M_2$ sur $H_2$, le **produit tensoriel**
$\mathcal M_1\otimes\mathcal M_2$ est le sous-espace engendré par les opérateurs $A\otimes B$ sur
$H_1\otimes H_2$, avec $A\in\mathcal M_1$, $B\in\mathcal M_2$. Comme
$(A\otimes B)(A'\otimes B')=AA'\otimes BB'$, $(A\otimes B)^*=A^*\otimes B^*$ et $1\otimes1=1$, c'est
une $*$-algèbre, donc une algèbre de von Neumann (théorème 2.2).

**Théorème 7.1.** (a) $H_1\otimes H_2=\bigoplus_{c,d}\mathbb C^{k_c}\otimes\mathbb C^{l_d}$, somme
orthogonale des images des projecteurs $P_c\otimes Q_d$, et
$$
\mathcal M_1\otimes\mathcal M_2=\bigoplus_{c,d}M_{k_c}\otimes M_{l_d}\cong\bigoplus_{c,d}M_{k_cl_d}(\mathbb C):
$$
elle se compose de tous les opérateurs $X=\sum_{c,d}(P_c\otimes Q_d)X(P_c\otimes Q_d)$, avec des
blocs arbitraires sur les $\mathbb C^{k_c}\otimes\mathbb C^{l_d}$. (b) Ses projecteurs centraux minimaux
sont les $C_1C_2$ projecteurs $P_c\otimes Q_d$, et
$\dim(\mathcal M_1\otimes\mathcal M_2)=\dim\mathcal M_1\dim\mathcal M_2$.
(c) Pour des traces $\operatorname{tr}_w$ sur $\mathcal M_1$ et $\operatorname{tr}_v$ sur $\mathcal M_2$ (théorème 3.1),
$\operatorname{tr}_w(A)\operatorname{tr}_v(B)=\operatorname{tr}_{w\otimes v}(A\otimes B)$ avec
$(w\otimes v)_{cd}=w_cv_d$. En particulier, `Tr_blunt`, `Tr_norm` et `tau_vN` du produit sont les
produits des mêmes traces des facteurs.

*Démonstration.* (a) Les $P_c\otimes Q_d$ sont des projecteurs deux à deux orthogonaux de somme
$1\otimes1$, et l'image de $P_c\otimes Q_d$ est $\mathbb C^{k_c}\otimes\mathbb C^{l_d}$. Comme $P_c$ et
$Q_d$ sont centraux, $(P_c\otimes Q_d)(A\otimes B)=AP_c\otimes BQ_d=(A\otimes B)(P_c\otimes Q_d)$, donc
tout élément de $\mathcal M_1\otimes\mathcal M_2$ est diagonal par blocs. Réciproquement, avec les
unités matricielles $E^c_{aa'}\in\mathcal M_1$ et $F^d_{bb'}\in\mathcal M_2$ des blocs $c$ et $d$, on a
$E^c_{aa'}\otimes F^d_{bb'}=|e_a\otimes e_b\rangle\langle e_{a'}\otimes e_{b'}|$ dans
$\mathbb C^{k_c}\otimes\mathbb C^{l_d}$ ; ce sont les unités matricielles de $M_{k_cl_d}(\mathbb C)$, donc
tout opérateur porté par le bloc $(c,d)$ appartient à $\mathcal M_1\otimes\mathcal M_2$. (b) Le centre
d'une somme directe d'algèbres de matrices pleines est engendré par les projecteurs sur les termes
de la somme (chapitre 2, exercice 1), et $\sum_{c,d}k_c^2l_d^2=\sum_ck_c^2\sum_dl_d^2$. (c)
$\operatorname{Tr}(X\otimes Y)=\sum_{a,b}X_{aa}Y_{bb}=\operatorname{Tr}X\operatorname{Tr}Y$, donc
$\operatorname{tr}_{w\otimes v}(A\otimes B)=\sum_{c,d}w_cv_d\operatorname{Tr}A_c\operatorname{Tr}B_d=\operatorname{tr}_w(A)\operatorname{tr}_v(B)$.
Les poids $1$, $1/k_c$, $1/(C_1k_c)$ des trois traces (section 3.1) se multiplient pour donner $1$,
$1/(k_cl_d)$, $1/(C_1C_2k_cl_d)$, les poids des mêmes traces sur le produit. $\square$

`tensor_product(alg1, alg2)` construit cette algèbre. Le secteur $(c,d)$ a pour indice $s=cC_2+d$
(indices à partir de $0$ dans le code) et pour charge $(q_c,r_d)$, et, à l'intérieur de ce
secteur, $e_a\otimes e_b$ a pour indice $a\,l_d+b$. `kron(A, B, alg12)` calcule $A\otimes B$. L'ordre
des secteurs diffère de l'ordre de Kronecker de $H_1\otimes H_2$ par une permutation de la base, que
`kron_rows` construit explicitement ci-dessous ; elle permet de comparer avec des produits de
Kronecker de matrices pleines.

```python
import math
import numpy as np
import torch
from torch_vn_algebra import (TypeIAlgebra, DensityMatrix, states, composite,
                              tensor_product, fused_tensor_product, kron, partial_trace)

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def full(X):
    """An operator of (+)_c M_{k_c} as a block-diagonal matrix on (+)_c C^{k_c}, shape (batch, D, D)."""
    return torch.stack([torch.block_diag(*[x[c, :k, :k] for c, k in enumerate(X.algebra.k_factors)])
                        for x in X.matrix])

def kron_rows(alg12):
    """For every sector of a product algebra, the positions of its basis vectors in H1 (x) H2 (Kronecker order)."""
    k1, k2 = alg12.tensor_factors[0].k_factors, alg12.tensor_factors[1].k_factors
    return [torch.tensor([(sum(k1[:c]) + a) * sum(k2) + sum(k2[:d]) + b
                          for c, d, off in pairs for a in range(k1[c]) for b in range(k2[d])])
            for pairs in alg12.product_blocks]              # pair blocks (c, d, offset), by offset

def full12(X):
    """An operator of a (plain or fused) product algebra as a matrix on H1 (x) H2 in Kronecker order."""
    f1, f2 = X.algebra.tensor_factors
    N = sum(f1.k_factors) * sum(f2.k_factors)
    out = torch.zeros(X.matrix.shape[0], N, N, dtype=X.matrix.dtype)
    for s, r in enumerate(kron_rows(X.algebra)):
        out[:, r[:, None], r[None, :]] = X.matrix[:, s, :len(r), :len(r)]
    return out

def random_operator(alg):
    return alg.from_blocks([torch.randn(k, k, dtype=dt) for k in alg.k_factors])

alg1 = TypeIAlgebra([1, 2], [1, 2], precision='double')            # C (+) M_2
alg2 = TypeIAlgebra([2, 3], [2, 3], precision='double')            # M_2 (+) M_3
alg12 = tensor_product(alg1, alg2)
print(alg12.k_factors, alg12.charges, sum(k * k for k in alg12.k_factors))     # 4 sectors, dimension 65
assert sum(k * k for k in alg12.k_factors) == 5 * 13                           # Theorem 7.1(b)

A, A2, B, B2 = random_operator(alg1), random_operator(alg1), random_operator(alg2), random_operator(alg2)
AB = kron(A, B, alg12)
assert close(full12(AB)[0], torch.kron(full(A)[0], full(B)[0]))               # Theorem 7.1(a)
assert close(AB @ kron(A2, B2, alg12), kron(A @ A2, B @ B2, alg12))
for tr in ('Tr_blunt', 'Tr_norm', 'tau_vN'):                                   # Theorem 7.1(c)
    assert (getattr(AB, tr)() - getattr(A, tr)() * getattr(B, tr)()).abs().max() < 1e-10
```

> **Interprétation physique.** Si chacun de deux systèmes obéit à sa propre règle de
> supersélection — deux boîtes qui ne peuvent pas échanger de particules —, le système composé a un
> secteur pour chaque *couple* de charges, et aucune observable ne relie des couples différents. Un
> secteur composé compte avec le produit des poids.

## 7.2 La trace partielle

**Définition.** Soit $\rho\in\mathcal M_1\otimes\mathcal M_2$. La **trace partielle**
$\operatorname{Tr}_2\rho\in\mathcal M_1$ est l'élément tel que
$$
\operatorname{Tr}\big((\operatorname{Tr}_2\rho)\,A\big)=\operatorname{Tr}\big(\rho\,(A\otimes1)\big)\qquad\text{for all }A\in\mathcal M_1 .
$$
Elle existe et est unique d'après le théorème 3.2, appliqué à la forme linéaire
$A\mapsto\operatorname{Tr}(\rho(A\otimes1))$ et à la trace fidèle $\operatorname{Tr}$ de $\mathcal M_1$ : la
trace partielle est la duale, relativement aux traces brutes, du plongement $A\mapsto A\otimes1$.
$\operatorname{Tr}_1$ se définit de la même façon. Pour les algèbres pleines $\mathcal M_i=B(H_i)$, c'est la
**trace partielle de l'espace de Hilbert** $\operatorname{Tr}_{H_2}:B(H_1\otimes H_2)\to B(H_1)$.

**Proposition 7.2.** (a) Pour $\rho=\bigoplus_{c,d}\rho_{cd}$,
$$
(\operatorname{Tr}_2\rho)_c=\sum_d\operatorname{tr}_{l_d}\rho_{cd},\qquad
(\operatorname{tr}_lX)_{aa'}=\sum_{b=1}^lX_{(a,b),(a',b)}\quad(X\text{ on }\mathbb C^k\otimes\mathbb C^l),
$$
et $\operatorname{Tr}_2\rho=\operatorname{Tr}_{H_2}\rho$, où $(\operatorname{Tr}_{H_2}X)_{ii'}=\sum_jX_{(i,j),(i',j)}$.
(b) $\operatorname{Tr}_2(A\otimes B)=(\operatorname{Tr}B)\,A$.
(c) $\operatorname{Tr}_2\big((A\otimes1)X(A'\otimes1)\big)=A\,(\operatorname{Tr}_2X)\,A'$ pour $A,A'\in\mathcal M_1$.
(d) $\operatorname{Tr}(\operatorname{Tr}_2X)=\operatorname{Tr}X$, et $\operatorname{Tr}_2$ envoie les états sur des états.

*Démonstration.* (a) Avec $A=\bigoplus_cA_c$, on a $A\otimes1=\bigoplus_{c,d}A_c\otimes1_{l_d}$ et
$$
\operatorname{Tr}\big(\rho(A\otimes1)\big)=\sum_{c,d}\sum_{a,a',b}(\rho_{cd})_{(a,b),(a',b)}(A_c)_{a'a}
=\sum_c\operatorname{Tr}\Big(\Big(\sum_d\operatorname{tr}_{l_d}\rho_{cd}\Big)A_c\Big),
$$
donc l'élément de droite vérifie la définition et est égal à $\operatorname{Tr}_2\rho$. Le même calcul
pour $B(H_i)$ donne la formule de $\operatorname{Tr}_{H_2}$ ; pour $\rho$ diagonale par blocs, les
coefficients $\rho_{(i,j),(i',j)}$ avec $i,i'$ dans des blocs différents sont nuls, donc
$\operatorname{Tr}_{H_2}\rho$ a les blocs de (a). (b) $\operatorname{tr}_l(X\otimes Y)=(\operatorname{Tr}Y)X$. (c) Pour
$B\in\mathcal M_1$,
$\operatorname{Tr}\big(\operatorname{Tr}_2((A\otimes1)X(A'\otimes1))B\big)=\operatorname{Tr}\big(X(A'BA\otimes1)\big)
=\operatorname{Tr}\big(A(\operatorname{Tr}_2X)A'B\big)$ ; on conclut par unicité. (d) On prend $A=1$. Si $\rho$ est
un état, $A\mapsto\operatorname{Tr}(\rho(A\otimes1))$ est positive, puisque $A^*A\otimes1=(A\otimes1)^*(A\otimes1)$,
et vaut $1$ en $A=1$, donc sa densité $\operatorname{Tr}_2\rho$ est un état (théorème 3.2). $\square$

Les canaux font l'objet du chapitre 6 ; ici, seule la définition suivante est nécessaire : une
application linéaire $\Phi$ est **complètement positive** si $\mathrm{id}_n\otimes\Phi$ est positive
sur $M_n(\mathbb C)\otimes(\cdot)$ pour tout $n$, et c'est un **canal** si elle préserve de plus la
trace.

**Théorème 7.3.** Faisons parcourir à $j$ la base canonique de $H_2$ et posons
$K_j=1\otimes\langle j|:H_1\otimes H_2\to H_1$, $K_j(x\otimes y)=\langle j,y\rangle x$. Alors
$$
\operatorname{Tr}_{H_2}X=\sum_jK_jXK_j^*,\qquad\sum_jK_j^*K_j=1,\qquad\sum_jK_j^*AK_j=A\otimes1 .
$$
Par conséquent, $\operatorname{Tr}_2$ est un canal de $\mathcal M_1\otimes\mathcal M_2$ dans $\mathcal M_1$, et sa
duale relativement aux traces brutes est $A\mapsto A\otimes1$.

*Démonstration.* $K_j^*x=x\otimes e_j$, donc $\langle e_a,K_jXK_j^*e_{a'}\rangle=X_{(a,j),(a',j)}$ ; la
somme sur $j$ donne la proposition 7.2(a). De plus,
$\sum_jK_j^*K_j(x\otimes y)=\sum_j\langle j,y\rangle x\otimes e_j=x\otimes y$
et $\sum_jK_j^*AK_j=\sum_jA\otimes|j\rangle\langle j|=A\otimes1$. Pour $Y\ge0$,
$(\mathrm{id}_n\otimes\operatorname{Tr}_{H_2})(Y)=\sum_j(1_n\otimes K_j)Y(1_n\otimes K_j)^*\ge0$, car
$\langle x,TYT^*x\rangle=\langle T^*x,YT^*x\rangle\ge0$ ; cela reste vrai sur la sous-algèbre
$\mathcal M_1\otimes\mathcal M_2$. Enfin,
$\operatorname{Tr}\sum_jK_jXK_j^*A=\operatorname{Tr}X\sum_jK_j^*AK_j$, ce qui, pour $A=1$, est la préservation de
la trace et, en général, la dualité. $\square$

`partial_trace(X, keep=1)` remet chaque bloc $\rho_{cd}$ sous la forme $(k_c,l_d,k_c,l_d)$ et prend
la trace ; `composite.partial_trace_channel(alg12, keep=1)` est le canal du théorème 7.3, avec les
opérateurs de Kraus $1_{k_c}\otimes\langle j|$ du secteur $(c,d)$ vers le secteur $c$. On vérifie les
deux par rapport à la définition, à un `numpy.einsum` indépendant et à la dualité :

```python
rho = states.random_density_matrix(alg12, batch_size=3)
rho1, rho2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
assert isinstance(rho1, DensityMatrix)                                      # states give states
R = full12(rho).numpy().reshape(3, 3, 5, 3, 5)                              # indices (batch, i, j, i', j')
assert np.abs(full(rho1).numpy() - np.einsum('xijkj->xik', R)).max() < 1e-12
assert np.abs(full(rho2).numpy() - np.einsum('xijil->xjl', R)).max() < 1e-12
one2 = alg2.identity()
assert (states.expectation(rho, kron(A, one2, alg12)) - states.expectation(rho1, A)).abs().max() < 1e-12
assert close(partial_trace(AB, keep=1), A * B.trace.item())                 # Proposition 7.2(b)

Ptr = composite.partial_trace_channel(alg12, keep=1)
print(Ptr)
assert Ptr.is_trace_preserving()
assert close(Ptr(rho), rho1)
assert close(Ptr.adjoint('blunt')(A), kron(A, one2, alg12))                 # Theorem 7.3: the dual
```

> **Interprétation physique.** $\operatorname{Tr}_2\rho$ est l'**état réduit** : tout ce qu'apprend un
> observateur qui ne peut mesurer que le premier système — la définition dit qu'il reproduit toutes
> les valeurs moyennes des observables $A\otimes1$. Des états composés ayant les mêmes états réduits
> diffèrent par leurs corrélations.

**La trace partielle comme espérance conditionnelle.** Soit $\mathcal M_2=M_m(\mathbb C)$, de sorte
que les secteurs de $\mathcal M_1\otimes M_m$ sont ceux de $\mathcal M_1$ (théorème 7.1) ; soient
$\iota(A)=A\otimes1_m$ et $\mathcal E=\frac1m\operatorname{Tr}_2:\mathcal M_1\otimes M_m\to\mathcal M_1$, et
soient $e_{ab}$ les unités matricielles de $M_m$.

**Proposition 7.4.** (a) $\iota\circ\mathcal E$ est une espérance conditionnelle de
$\mathcal M_1\otimes M_m$ sur $\mathcal M_1\otimes1$ (section 2.6) et, pour toute trace fidèle
$\operatorname{tr}_w$ de $\mathcal M_1\otimes M_m$, c'est la seule qui préserve $\operatorname{tr}_w$. (b) Les
$m^2$ éléments $u_{ab}=\sqrt m\,(1\otimes e_{ab})$ vérifient
$X=\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*X))$ pour tout $X$, et
$\sum_{a,b}u_{ab}u_{ab}^*=m^2\,1$. (c) Toute famille $(v_k)$ telle que
$X=\sum_kv_k\,\iota(\mathcal E(v_k^*X))$ pour tout $X$ vérifie $\sum_kv_kv_k^*=m^2\,1$.

*Démonstration.* (a) $\mathcal E$ est positive (théorème 7.3), et $\iota$ aussi
($T^*T\otimes1=(T\otimes1)^*(T\otimes1)$) ; $\mathcal E(A\otimes1)=A$ d'après la proposition 7.2(b), et
la propriété de bimodule est la proposition 7.2(c). Notons encore $\operatorname{tr}_w$ la trace de mêmes
poids $w_c$ sur $\mathcal M_1$ ; alors $\operatorname{tr}_w(Y\otimes1)=m\operatorname{tr}_w(Y)$. La définition de
$\operatorname{Tr}_2$ appliquée à $AP_c$, multipliée par $w_c$ et sommée sur $c$, donne
$\operatorname{tr}_w(X(A\otimes1))=\operatorname{tr}_w((\operatorname{Tr}_2X)A)$, et $A=1$ donne
$\operatorname{tr}_w(\iota\mathcal E(X))=\operatorname{tr}_w(X)$. Si $F=\iota\circ F_0$ est une autre espérance
conditionnelle sur $\mathcal M_1\otimes1$ préservant $\operatorname{tr}_w$, alors
$\operatorname{tr}_w(X(A\otimes1))=\operatorname{tr}_w(F(X)(A\otimes1))=m\operatorname{tr}_w(F_0(X)A)$, donc
$\operatorname{tr}_w((\operatorname{Tr}_2X-mF_0(X))A)=0$ pour tout $A\in\mathcal M_1$, et $F_0=\mathcal E$ par fidélité
(prendre $A=(\operatorname{Tr}_2X-mF_0(X))^*$). (b) Écrivons $X=\sum_{a,b}X_{ab}\otimes e_{ab}$ avec
$X_{ab}\in\mathcal M_1$ (théorème 7.1(a)). D'après la proposition 7.2(b),
$\mathcal E(u_{ab}^*X)=\frac{\sqrt m}m\sum_{b'}\operatorname{Tr}_2(X_{ab'}\otimes e_{bb'})=X_{ab}/\sqrt m$, donc
$u_{ab}\,\iota(\mathcal E(u_{ab}^*X))=X_{ab}\otimes e_{ab}$ ; et
$\sum_{a,b}u_{ab}u_{ab}^*=m\sum_{a,b}1\otimes e_{aa}=m^21$. (c) $\operatorname{Tr}_2(Y^*)=(\operatorname{Tr}_2Y)^*$
d'après la proposition 7.2(a), donc l'adjoint du développement de $Y=u_{ab}$ est
$u_{ab}^*=\sum_k\iota(\mathcal E(u_{ab}^*v_k))\,v_k^*$, et, d'après (b),
$\sum_{a,b}u_{ab}u_{ab}^*=\sum_k\big(\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*v_k))\big)v_k^*=\sum_kv_kv_k^*$. $\square$

Une famille comme en (b) est une *quasi-base* de $\mathcal E$, et $\sum_kv_kv_k^*$, qui d'après (c)
n'en dépend pas, est l'**indice de Watatani** $\operatorname{Ind}\mathcal E=m^2$ (Y. Watatani, Mem. Amer.
Math. Soc. 83, no. 424, 1990). Pour un facteur, $M_k\otimes1\subset M_k\otimes M_m$ contient
$\mathbb C^k$ avec la multiplicité $m$ dans $\mathbb C^{km}$ : sa matrice d'inclusion est
$\Lambda=(m)$, et $m^2=\|\Lambda\|^2=\dim(M_k\otimes M_m)/\dim M_k$ est son indice de Jones
$[\mathcal M:\mathcal N]$ (V. F. R. Jones, Invent. Math. 72, 1, 1983 ; F. M. Goodman, P. de la Harpe
et V. F. R. Jones, *Coxeter graphs and towers of algebras*, Springer, 1989). Pour $\mathcal M_1$
quelconque, l'inclusion est la somme directe des $C_1$ inclusions
$M_{k_c}\otimes1\subset M_{k_c}\otimes M_m$, $\Lambda=m1_{C_1}$. Le module
`torch_vn_algebra.inclusions` implémente les inclusions, leurs espérances conditionnelles et leurs
indices ; `from_tensor_product(alg1, m)` construit celle-ci :

```python
from torch_vn_algebra import inclusions
m = 3
inc = inclusions.from_tensor_product(alg1, m)                # alg1 (x) 1 c alg1 (x) M_3
X = random_operator(inc.M)
for trace in ('Tr_blunt', 'Tr_norm', 'tau_vN'):
    E = inc.conditional_expectation(trace.replace('Tr_', ''))
    assert close(E(X), partial_trace(X, keep=1) * (1 / m))   # Proposition 7.4(a): the same E for every trace
    assert abs(getattr(inc.embed(E(X)), trace)() - getattr(X, trace)()).max() < 1e-10
    assert close(E.index(), inc.M.identity() * m**2)         # Proposition 7.4(c)

def unit(a, b):
    e = torch.zeros(m, m, dtype=dt)
    e[a, b] = 1
    return math.sqrt(m) * kron(alg1.identity(), inc.M.tensor_factors[1].from_blocks([e]), inc.M)

u = [unit(a, b) for a in range(m) for b in range(m)]         # the quasi-basis of Proposition 7.4(b)
assert close(sum((v @ inc.embed(E(v.adjoint() @ X))).matrix for v in u), X)
assert close(sum((v @ v.adjoint()).matrix for v in u), inc.M.identity() * m**2)
print(inc.matrix.tolist(), inc.component_indices())           # Lambda = 3 * 1_2, index 9 per component
```

## 7.3 Fusion des charges

Supposons maintenant que le secteur $c$ de $\mathcal M_1$ porte une charge $q_c$, le secteur $d$ de
$\mathcal M_2$ une charge $r_d$, et qu'une **règle de fusion** $f$ donne la charge totale
$f(q_c,r_d)$ — l'addition pour des nombres de particules, l'addition modulo $2$ pour la parité
fermionique, ou toute autre application. Posons
$$
H_q=\bigoplus_{f(q_c,r_d)=q}\mathbb C^{k_c}\otimes\mathbb C^{l_d},\qquad
\Pi_q=\sum_{f(q_c,r_d)=q}P_c\otimes Q_d,\qquad n_q=\sum_{f(q_c,r_d)=q}k_cl_d .
$$
Le **produit fusionné** $\mathcal M_1\otimes_f\mathcal M_2=\{\Pi_q\}'$ est formé des opérateurs sur
$H_1\otimes H_2$ qui conservent la charge totale.

**Proposition 7.5.** (a) $\mathcal M_1\otimes_f\mathcal M_2=\bigoplus_qB(H_q)\cong\bigoplus_qM_{n_q}(\mathbb C)$,
avec les projecteurs centraux minimaux $\Pi_q$. (b) $\mathcal M_1\otimes\mathcal M_2\subseteq\mathcal M_1\otimes_f\mathcal M_2$,
avec égalité si et seulement si chaque $H_q$ ne contient qu'un seul couple $(c,d)$. (c) Si les
charges sont des nombres réels et si $f$ est l'addition, alors $\mathcal M_1\otimes_f\mathcal M_2=\{N\}'$
avec $N=N_1\otimes1+1\otimes N_2$, $N_1=\sum_cq_cP_c$, $N_2=\sum_dr_dQ_d$.

*Démonstration.* (a) Chaque couple appartient à exactement un $H_q$, donc les $\Pi_q$ sont des
projecteurs deux à deux orthogonaux de somme $1$, et $X$ commute avec tous si et seulement si
$X=\sum_q\Pi_qX\Pi_q$, avec des blocs arbitraires dans $B(H_q)$ ; le centre est comme dans le
théorème 7.1(b). (b) Tout élément de $\mathcal M_1\otimes\mathcal M_2$ commute avec chaque
$P_c\otimes Q_d$, donc avec leurs sommes $\Pi_q$. Si chaque $H_q$ est un seul
$\mathbb C^{k_c}\otimes\mathbb C^{l_d}$, les deux algèbres sont égales à
$\bigoplus_{c,d}B(\mathbb C^{k_c}\otimes\mathbb C^{l_d})$ ; si $H_q$ contient deux couples, $B(H_q)$
contient des opérateurs qui envoient le bloc d'un couple dans celui de l'autre, et qui ne sont pas
diagonaux par blocs. (c) $N=\sum_{c,d}(q_c+r_d)P_c\otimes Q_d=\sum_qq\,\Pi_q$, donc les $\Pi_q$ non
nuls sont les projecteurs spectraux de $N$. Ce sont des polynômes en $N$ (démonstration du
théorème 1.5(a)) et $N$ en est une combinaison linéaire, donc $X$ commute avec $N$ si et seulement
s'il commute avec chaque $\Pi_q$. $\square$

`fused_tensor_product(alg1, alg2, fuse)` (par défaut, `fuse` est l'addition) ordonne ses secteurs
selon la charge fusionnée, et les couples $(c,d)$ à l'intérieur d'un secteur dans l'ordre
lexicographique, chacun avec la base $e_a\otimes e_b$. Exemple : les bosons de cœur dur. Les
opérateurs qui conservent le nombre de particules de deux sites forment
$M_1\oplus M_2\oplus M_1$ (charges $0,1,2$), ceux de trois sites $M_1\oplus M_3\oplus M_3\oplus M_1$ ;
le produit fusionné est celui de cinq sites, à l'ordre de la base près :
$n_q=\sum_c\binom2c\binom3{q-c}=\binom5q$.

```python
bos1 = TypeIAlgebra([1, 2, 1], [1, 2, 1], precision='double', charges=[0, 1, 2])
bos2 = TypeIAlgebra([1, 3, 3, 1], [1, 3, 3, 1], precision='double', charges=[0, 1, 2, 3])
F = fused_tensor_product(bos1, bos2)                         # U(1): total particle number
P = tensor_product(bos1, bos2)                               # particle numbers conserved separately
print(F.charges, F.k_factors, sum(k * k for k in F.k_factors), sum(k * k for k in P.k_factors))
assert F.k_factors == [math.comb(5, q) for q in range(6)]

N1, N2 = full(bos1.central([0.0, 1.0, 2.0]))[0], full(bos2.central([0.0, 1.0, 2.0, 3.0]))[0]
N = torch.kron(N1, torch.eye(8, dtype=dt)) + torch.kron(torch.eye(4, dtype=dt), N2)
X = full12(F.operator(torch.randn(1, F.C, F.k_max, F.k_max, dtype=dt)))[0]
assert close(X @ N, N @ X)                                   # Proposition 7.5(c)
A, B = random_operator(bos1), random_operator(bos2)
assert close(full12(kron(A, B, F)), full12(kron(A, B, P)))   # Proposition 7.5(b)
```

L'algèbre fusionnée est de dimension $252$, le produit simple de dimension $120$ ; la différence
est formée des opérateurs qui déplacent des particules d'une partie à l'autre, comme les termes
de saut.

Un état de $\mathcal M_1\otimes_f\mathcal M_2$ peut présenter des cohérences entre les blocs de
couples différents, si bien que la proposition 7.2 ne s'applique pas. Soit $E_1(X)=\sum_cP_cXP_c$
sur $B(H_1)$ (l'espérance conditionnelle de la section 2.6, Interprétation physique) ; son image
est $\mathcal M_1$, et $\operatorname{Tr}(E_1(X)A)=\operatorname{Tr}(XA)$ pour $A\in\mathcal M_1$, puisque
$A=\sum_cP_cAP_c$.

**Théorème 7.6.** Soit $\rho$ la densité d'un état $\omega$ sur $\mathcal M_1\otimes_f\mathcal M_2$. La
densité de sa restriction $A\mapsto\omega(A\otimes1)$ à $\mathcal M_1$ est
$$
\rho_1=E_1\big(\operatorname{Tr}_{H_2}\rho\big).
$$
Le bloc $(c,c')$ de $\operatorname{Tr}_{H_2}\rho$ est $\sum_d\operatorname{tr}_{l_d}\rho_{(c,d),(c',d)}$, la somme
portant sur les $d$ tels que $f(q_c,r_d)=f(q_{c'},r_d)$, où $\rho_{(c,d),(c',d)}$ est le bloc de
$\rho$ allant du couple $(c',d)$ au couple $(c,d)$. Par conséquent, $E_1$ agit trivialement sur tous
les $\operatorname{Tr}_{H_2}\rho$ si et seulement s'il n'existe pas de $c\neq c'$ et de $d$ tels que
$f(q_c,r_d)=f(q_{c'},r_d)$ ; c'est le cas par exemple si $f$ est l'addition et si les $q_c$ sont
deux à deux distincts.

*Démonstration.* $A\otimes1\in\mathcal M_1\otimes_f\mathcal M_2$ (proposition 7.5(b)) et
$\omega(A\otimes1)=\operatorname{Tr}(\rho(A\otimes1))=\operatorname{Tr}\big((\operatorname{Tr}_{H_2}\rho)A\big)
=\operatorname{Tr}\big(E_1(\operatorname{Tr}_{H_2}\rho)A\big)$. Comme $E_1(\operatorname{Tr}_{H_2}\rho)\in\mathcal M_1$, c'est
la densité, par l'unicité dans le théorème 3.2. La formule par blocs est la proposition 7.2(a)
pour $\operatorname{Tr}_{H_2}$, avec $i$ dans le bloc $c$, $i'$ dans le bloc $c'$ et $j$ dans le bloc $d$ ; le
bloc de $\rho$ entre $(c',d)$ et $(c,d)$ est nul sauf si les deux couples appartiennent au même
$H_q$. S'il existe $c\neq c'$ et $d$ tels que $f(q_c,r_d)=f(q_{c'},r_d)=q$, prenons des vecteurs
unitaires $x\in\mathbb C^{k_c}$, $x'\in\mathbb C^{k_{c'}}$, $y\in\mathbb C^{l_d}$ et
$\psi=(x+x')\otimes y/\sqrt2\in H_q$. Alors $|\psi\rangle\langle\psi|\in B(H_q)$, et
$\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|=\frac12|x+x'\rangle\langle x+x'|$ possède le bloc non nul
$\frac12|x\rangle\langle x'|$, que $E_1$ supprime. $\square$

Les deux implémentations de la bibliothèque calculent cette application, pour toute règle de
fusion : `partial_trace` ne prend la trace que des blocs diagonaux de couples
$\rho_{(c,d),(c,d)}$, et `partial_trace_channel` utilise des indices de Kraus distincts pour les
couples d'un même secteur ayant le même $c$ ; sa duale est de nouveau $A\mapsto A\otimes1$. Avec
l'addition, $E_1$ n'a aucun effet dans notre exemple ; avec la règle de parité
$f(q,r)=q+r\bmod2$, les charges $0$ et $2$ de `bos1` fusionnent de la même manière avec tout $r$, et
$E_1$ supprime de véritables cohérences :

```python
G = fused_tensor_product(bos1, bos2, fuse=lambda q, r: (q + r) % 2)       # only the parity is conserved
print(G.k_factors)

def block_mask(alg):
    """Ones in the diagonal blocks of (+)_c C^{k_c}: X * block_mask(alg) = sum_c P_c X P_c."""
    return full(alg.operator(torch.ones(1, alg.C, alg.k_max, alg.k_max, dtype=dt)))

for alg, coherent in ((F, False), (G, True)):
    rho = states.random_density_matrix(alg, batch_size=2)
    R = full12(rho).reshape(2, 4, 8, 4, 8)
    for keep, X_H, alg_k in ((1, torch.einsum('xijkj->xik', R), bos1), (2, torch.einsum('xijil->xjl', R), bos2)):
        E_X = X_H * block_mask(alg_k)                                     # E(Tr_H rho), Theorem 7.6
        assert close(full(partial_trace(rho, keep=keep)), E_X)
        assert close(full(composite.partial_trace_channel(alg, keep=keep)(rho)), E_X)
        lost = (X_H - E_X).abs().max().item()                             # coherences removed by E
        print(f"keep={keep}: largest coherence removed: {lost:.3e}")
        assert (lost > 1e-3) if coherent else (lost == 0)
    assert close(composite.partial_trace_channel(alg, keep=1).adjoint('blunt')(A), kron(A, bos2.identity(), alg))
```

> **Interprétation physique.** Dans la description BCS d'un supraconducteur, seule la parité du
> nombre d'électrons est conservée, de sorte que l'état global peut superposer les nombres de
> particules $0$ et $2$ d'une région. Un observateur dont les observables conservent le nombre de
> particules de la région ne peut pas voir cette cohérence : pour lui, l'état réduit est l'état
> rendu diagonal par blocs du théorème 7.6, un mélange de nombres de particules.

## 7.4 Intrication des états purs

**Lemme 7.7 (décomposition de Schmidt).** Tout vecteur unitaire $\psi\in\mathbb C^a\otimes\mathbb C^b$
peut s'écrire $\psi=\sum_{i=1}^rs_i\,u_i\otimes v_i$ avec $s_i>0$, $\sum_is_i^2=1$ et des familles
orthonormées $u_i\in\mathbb C^a$, $v_i\in\mathbb C^b$. Alors
$\operatorname{Tr}_2|\psi\rangle\langle\psi|=\sum_is_i^2|u_i\rangle\langle u_i|$
et $\operatorname{Tr}_1|\psi\rangle\langle\psi|=\sum_is_i^2|v_i\rangle\langle v_i|$ : les deux états réduits ont
les mêmes valeurs propres non nulles $s_i^2$.

*Démonstration.* Soit $\Psi$ la matrice $a\times b$ définie par $\Psi_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,\psi\rangle$.
D'après le théorème spectral, $\Psi\Psi^*=\sum_{i=1}^a\lambda_i|u_i\rangle\langle u_i|$ avec une base
orthonormée $u_i$ et $\lambda_i=\|\Psi^*u_i\|^2\ge0$, ordonnés de sorte que $\lambda_i>0$ exactement
pour $i\le r$. Posons $s_i=\sqrt{\lambda_i}$ et $w_i=\Psi^*u_i/s_i$ pour $i\le r$ ; alors
$\langle w_i,w_j\rangle=\langle u_i,\Psi\Psi^*u_j\rangle/(s_is_j)=\delta_{ij}$, et $\Psi^*u_i=0$ pour $i>r$.
Donc $\Psi=\sum_i|u_i\rangle\langle u_i|\Psi=\sum_{i\le r}|u_i\rangle\langle\Psi^*u_i|=\sum_{i\le r}s_i|u_i\rangle\langle w_i|$,
c'est-à-dire $\Psi_{\alpha\beta}=\sum_is_i(u_i)_\alpha\overline{(w_i)_\beta}$ : c'est l'assertion avec
$v_i=\bar w_i$, et $\sum_is_i^2=\operatorname{Tr}\Psi\Psi^*=\|\psi\|^2=1$. Enfin,
$|\psi\rangle\langle\psi|=\sum_{i,j}s_is_j|u_i\rangle\langle u_j|\otimes|v_i\rangle\langle v_j|$, et la
proposition 7.2(b) donne les états réduits. $\square$

L'**entropie d'intrication** de $\psi$ est
$E(\psi)=S(\operatorname{Tr}_2|\psi\rangle\langle\psi|)=S(\operatorname{Tr}_1|\psi\rangle\langle\psi|)=-\sum_is_i^2\log s_i^2$.
Elle est nulle si et seulement si $r=1$, c'est-à-dire si et seulement si $\psi=u\otimes v$, et
$E(\psi)\le\log\min(a,b)$, avec égalité si et seulement si l'état réduit sur le plus petit facteur
est maximalement mélangé (corollaire 3.6).

En présence de charges, un état pur appartient à un secteur $H_q$ du produit fusionné ; écrivons
$\psi=\sum_{(c,d)}\psi_{cd}$ avec $\psi_{cd}\in\mathbb C^{k_c}\otimes\mathbb C^{l_d}$.

**Proposition 7.8 (intrication résolue en symétrie).** Supposons que $f$ soit l'addition, que les
$q_c$ soient deux à deux distincts et que les $r_d$ soient deux à deux distincts, de sorte que,
pour chaque $c$, au plus un $d=d(c)$ vérifie $q_c+r_d=q$. Posons $\psi_c=\psi_{c\,d(c)}$ (ou $0$),
$p_c=\|\psi_c\|^2$ et, pour $p_c>0$, $\hat\psi_c=\psi_c/\sqrt{p_c}$. Alors
$\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|$ appartient à $\mathcal M_1$ (il commute avec $N_1$) et vaut
$\bigoplus_cp_c\operatorname{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|$ ; l'état réduit sur $\mathcal M_2$ a la
probabilité de secteur $p_c$ en $d(c)$ ; et
$$
E(\psi)=H(p)+\sum_cp_c\,E(\hat\psi_c),\qquad H(p)=-\sum_cp_c\log p_c .
$$

*Démonstration.* D'après la formule par blocs du théorème 7.6, le bloc $(c,c')$ de
$\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|$ est la somme des $\operatorname{tr}_{l_d}|\psi_{cd}\rangle\langle\psi_{c'd}|$
sur les $d$ tels que $q_c+r_d=q=q_{c'}+r_d$, ce qui impose $c=c'$ et $d=d(c)$. L'état réduit est
donc $\bigoplus_c\operatorname{Tr}_2|\psi_c\rangle\langle\psi_c|$, de traces de blocs $p_c$ ; en échangeant les
rôles des facteurs, on obtient l'énoncé pour $\mathcal M_2$. D'après le lemme 7.7, $E(\psi)$ est
l'entropie de cet état, et le théorème 3.3 avec $\sigma_c=\operatorname{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|$,
$S(\sigma_c)=E(\hat\psi_c)$, donne la formule. $\square$

$H(p)$ est l'**intrication de nombre** et $\sum_cp_cE(\hat\psi_c)$ l'**intrication
configurationnelle**. Pour un état aléatoire de $2$ particules sur les cinq sites, on calcule
$E(\psi)$ à partir des traces partielles des deux côtés et des valeurs singulières de la matrice
$\Psi$ de taille $4\times8$ :

```python
s = F.charges.index(2)                                      # the sector of two particles, n_q = 10
n = F.k_factors[s]
psi = torch.randn(n, dtype=dt)
psi = psi / psi.norm()
mat = torch.zeros(1, F.C, F.k_max, F.k_max, dtype=dt)
mat[0, s, :n, :n] = torch.outer(psi, psi.conj())
rho = DensityMatrix(F, matrix=mat)
r1, r2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)

psi_H = torch.zeros(4 * 8, dtype=dt)
psi_H[kron_rows(F)[s]] = psi                                # psi as a vector of H1 (x) H2
lam = torch.linalg.svdvals(psi_H.reshape(4, 8)) ** 2        # squared Schmidt coefficients
lam = lam[lam > 1e-14]
E_psi = -(lam * lam.log()).sum().item()
assert abs(r1.von_neumann_entropy().item() - E_psi) < 1e-10
assert abs(r2.von_neumann_entropy().item() - E_psi) < 1e-10  # Lemma 7.7

def entropy(sigma):
    w = torch.linalg.eigvalsh(sigma).clamp(min=1e-300)
    return -(w * w.log()).sum().item()

p = r1.sector_probabilities()[0]                            # c = 0, 1, 2 particles in the first part
assert torch.allclose(r2.sector_probabilities()[0, [2, 1, 0]], p)   # d(c) = 2 - c
H_num = -(p * p.log()).sum().item()
S_conf = sum(p[c].item() * entropy(r1.matrix[0, c, :k, :k] / p[c]) for c, k in enumerate(bos1.k_factors))
print(f"E = {E_psi:.4f} = H(p) {H_num:.4f} + configurational {S_conf:.4f}")
assert abs(E_psi - H_num - S_conf) < 1e-10                  # Proposition 7.8
```

> **Interprétation physique.** Des parties soumises à la règle de supersélection du nombre de
> particules ne peuvent utiliser comme ressource que la partie configurationnelle : c'est
> l'*intrication de particules* de H. M. Wiseman et J. A. Vaccaro (Phys. Rev. Lett. 91, 097902,
> 2003). $H(p)$ est mesurée dans les atomes froids sous le nom d'*entropie de nombre* (A. Lukin et
> al., Science 364, 256, 2019).

## 7.5 États aléatoires et formule de Page

Soit $g\in\mathbb C^N$ à coordonnées indépendantes, dont les parties réelle et imaginaire sont
indépendantes de loi $N(0,\tfrac12)$, de sorte que $\mathbb E|g_i|^2=1$. Un **état pur aléatoire de
Haar** est $\psi=g/\|g\|$ ; sa loi est invariante par les unitaires, comme le sont la densité
$\pi^{-N}e^{-\|g\|^2}$ et la mesure de Lebesgue (c'est l'unique loi invariante sur la sphère, fait
que nous n'utiliserons pas).

**Lemme 7.9.** $\|g\|$ et $\psi=g/\|g\|$ sont indépendants.

*Démonstration.* Soient $h\ge0$ et $\phi\ge0$ mesurables, respectivement sur $(0,\infty)$ et sur la
sphère unité $S^{2N-1}$. L'intégration en coordonnées polaires dans $\mathbb R^{2N}$, $\sigma$
désignant la mesure de surface, donne
$$
\mathbb E\big[h(\|g\|)\phi(\psi)\big]=\pi^{-N}\int_0^\infty h(t)e^{-t^2}t^{2N-1}dt\int_{S^{2N-1}}\phi\,d\sigma=:\alpha(h)\beta(\phi).
$$
En prenant $\phi=1$, puis $h=1$, puis les deux, on obtient $\mathbb Eh(\|g\|)=\alpha(h)\beta(1)$,
$\mathbb E\phi(\psi)=\alpha(1)\beta(\phi)$ et $\alpha(1)\beta(1)=1$, d'où
$\mathbb E[h(\|g\|)\phi(\psi)]=\mathbb Eh(\|g\|)\,\mathbb E\phi(\psi)$. $\square$

**Lemme 7.10.** $\mathbb E\|g\|^4=N(N+1)$. Si $G$ est une matrice $a\times b$ à coefficients
indépendants distribués comme les $g_i$ (une matrice *de Ginibre*), alors
$\mathbb E\operatorname{Tr}(GG^*)^2=ab(a+b)$.

*Démonstration.* Avec $g_i=x+iy$, $\mathbb E|g_i|^4=\mathbb Ex^4+2\mathbb Ex^2\mathbb Ey^2+\mathbb Ey^4=\tfrac34+\tfrac12+\tfrac34=2$,
donc $\mathbb E\|g\|^4=\sum_i\mathbb E|g_i|^4+\sum_{i\ne j}\mathbb E|g_i|^2\mathbb E|g_j|^2=2N+N(N-1)$. Ensuite,
$\operatorname{Tr}(GG^*)^2=\sum_{i,i',j,j'}G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}$. Les coefficients et leurs
conjugués sont des combinaisons linéaires complexes de gaussiennes réelles indépendantes, avec
$\mathbb E[G_{ij}G_{kl}]=0$ et $\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}$, donc la formule
d'Isserlis (démonstration de la proposition 1.10) donne
$\mathbb E[G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}]=\delta_{ii'}+\delta_{jj'}$. En sommant, on obtient
$ab^2+a^2b$. $\square$

**Théorème 7.11 (E. Lubkin, J. Math. Phys. 19, 1028, 1978).** Pour $\psi\in\mathbb C^a\otimes\mathbb C^b$
aléatoire de Haar et $\rho_A=\operatorname{Tr}_2|\psi\rangle\langle\psi|$,
$$
\mathbb E\operatorname{Tr}\rho_A^2=\frac{a+b}{ab+1}.
$$

*Démonstration.* Soit $G_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,g\rangle$, une matrice de Ginibre ;
la matrice de $\psi$ est $\Psi=G/\|g\|$ et, d'après la démonstration du lemme 7.7,
$\rho_A=\Psi\Psi^*=GG^*/\|g\|^2$. Donc $\operatorname{Tr}(GG^*)^2=\|g\|^4\operatorname{Tr}\rho_A^2$, où
$\operatorname{Tr}\rho_A^2\in[0,1]$ est une fonction de $\psi$. D'après le lemme 7.9,
$\mathbb E\operatorname{Tr}(GG^*)^2=\mathbb E\|g\|^4\,\mathbb E\operatorname{Tr}\rho_A^2$, et d'après le lemme 7.10,
$\mathbb E\operatorname{Tr}\rho_A^2=ab(a+b)/(ab(ab+1))$. $\square$

Le chapitre 6 obtient la même moyenne à partir des moments d'ordre quatre d'un vecteur de Haar,
comme pureté de sortie d'un canal aléatoire.

**Corollaire 7.12.** $\log\frac{ab+1}{a+b}\le\mathbb E\,S(\rho_A)\le\log\min(a,b)$.

*Démonstration.* D'après le corollaire 3.6, $S(\rho_A)\le\log a$, et d'après le lemme 7.7,
$S(\rho_A)=S(\operatorname{Tr}_1|\psi\rangle\langle\psi|)\le\log b$. L'inégalité de Jensen pour la fonction
concave $\log$, avec les poids $\lambda_i>0$, valeurs propres non nulles de $\rho_A$, donne
$\sum_i\lambda_i\log\lambda_i\le\log\sum_i\lambda_i^2$, c'est-à-dire $S(\rho_A)\ge-\log\operatorname{Tr}\rho_A^2$ ;
l'inégalité de Jensen pour la fonction convexe $-\log$ et le théorème 7.11 donnent
$\mathbb E[-\log\operatorname{Tr}\rho_A^2]\ge-\log\mathbb E\operatorname{Tr}\rho_A^2$. $\square$

Pour $a$ fixé et $b\to\infty$, les bornes se rejoignent : une petite partie d'un état pur aléatoire
est presque maximalement mélangée. La moyenne exacte est donnée par la formule de Page, que nous
énonçons sans démonstration.

**Théorème 7.13 (formule de Page).** Pour $a\le b$,
$$
\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}.
$$
Elle a été conjecturée par D. N. Page (Phys. Rev. Lett. 71, 1291, 1993) et démontrée par
S. K. Foong et S. Kanno (Phys. Rev. Lett. 72, 1148, 1994), J. Sánchez-Ruiz (Phys. Rev. E 52, 5653,
1995) et S. Sen (Phys. Rev. Lett. 77, 1, 1996). Pour $1\ll a\le b$, elle vaut
$\approx\log a-\frac a{2b}$ ; pour deux qubits, elle vaut $\frac13$.

La loi de $\rho_A$ est la **mesure induite** (K. Życzkowski et H.-J. Sommers, J. Phys. A 34,
7111, 2001). `states.random_density_matrix(alg, rank=r)` forme $G_cG_c^*$ avec des blocs de Ginibre
$G_c$ de taille $k_c\times r$ et divise par $\sum_c\operatorname{Tr}G_cG_c^*$. Pour un seul bloc, c'est la loi
de $\rho_A$ ci-dessus avec $b=r$. Pour plusieurs blocs, c'est la restriction à $\mathcal M_1$ d'un
état pur aléatoire de Haar de $H_1\otimes\mathbb C^r$ : il suffit d'appliquer le théorème 7.6 avec
$\mathcal M_2=M_r$ et une règle de fusion constante, de sorte que
$\mathcal M_1\otimes_f\mathcal M_2=B(H_1\otimes\mathbb C^r)$, avec pour $G_c$ les lignes de $G$ du bloc $c$.

**Proposition 7.14.** Pour $\rho$ tiré par `random_density_matrix(alg, rank=r)` sur
$\bigoplus_cM_{k_c}$, avec $D=\sum_ck_c$,
$$
\mathbb E\operatorname{Tr}\rho^2=\frac{\sum_ck_c(k_c+r)}{D(Dr+1)} .
$$

*Démonstration.* $\operatorname{Tr}\rho^2=\sum_c\operatorname{Tr}(G_cG_c^*)^2/\|g\|^4$, où $g\in\mathbb C^{Dr}$ est la
collection de tous les coefficients ; comme dans le théorème 7.11, les lemmes 7.9 et 7.10 donnent
$\sum_ck_cr(k_c+r)/(Dr(Dr+1))$. $\square$

On vérifie les théorèmes 7.11 et 7.13 sur $10^4$ états purs de $\mathbb C^3\otimes\mathbb C^6$ réduits
avec `partial_trace` et sur la mesure induite, ainsi que la proposition 7.14 pour
$M_1\oplus M_2\oplus M_3$, $r=2$ :

```python
a, b, n_samples = 3, 6, 10_000
g = torch.randn(n_samples, a * b, dtype=dt)                 # complex Gaussian: E|g_i|^2 = 1
psi = g / g.norm(dim=-1, keepdim=True)
AB_alg = tensor_product(TypeIAlgebra([a], [a], precision='double'), TypeIAlgebra([b], [b], precision='double'))
pure = DensityMatrix(AB_alg, matrix=(psi[:, :, None] * psi.conj()[:, None, :])[:, None],
                     validate=False)                        # pure by construction
rho_A = partial_trace(pure, keep=1)
purity_exact = (a + b) / (a * b + 1)
page = sum(1 / k for k in range(b + 1, a * b + 1)) - (a - 1) / (2 * b)

def mean_se(x):
    return x.mean().item(), x.std().item() / math.sqrt(x.numel())

for name, rho_s in (('reduced pure states', rho_A),
                    ('induced measure', states.random_density_matrix(rho_A.algebra, 40_000, rank=b))):
    (m_p, e_p), (m_S, e_S) = mean_se(rho_s.purity()), mean_se(rho_s.von_neumann_entropy())
    print(f"{name}: E Tr rho^2 = {m_p:.4f} +- {e_p:.4f} (exact {purity_exact:.4f}), "
          f"E S = {m_S:.4f} +- {e_S:.4f} (Page {page:.4f})")
    assert abs(m_p - purity_exact) < 5 * e_p and abs(m_S - page) < 5 * e_S     # Theorems 7.11, 7.13
    assert math.log((a * b + 1) / (a + b)) < m_S < math.log(a)                # Corollary 7.12

alg3 = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')
m_p, e_p = mean_se(states.random_density_matrix(alg3, 40_000, rank=2).purity())
print(f"E Tr rho^2 = {m_p:.4f} +- {e_p:.4f} (exact {26 / 78:.4f})")
assert abs(m_p - 26 / 78) < 5 * e_p                                          # Proposition 7.14
```

> **Interprétation physique.** Un état pur typique est presque maximalement intriqué entre toute
> partie d'au plus la moitié du système et le reste. Quand la partie grandit, l'intrication croît
> comme $\log a$ puis décroît au-delà de la moitié : c'est la *courbe de Page*, qui décrit comment
> l'information quitte un trou noir dans le rayonnement de Hawking, et pourquoi les parties d'un
> système fermé paraissent thermiques.

## 7.6 Points de comparaison

Sur un seul espace de Hilbert, QuTiP (`Qobj.ptrace`, `entropy_vn`, `rand_ket`, `rand_dm` avec la
distribution de Ginibre et un rang) et `numpy` (un `reshape` et un `einsum`, comme ci-dessus)
calculent traces partielles, entropies et états induits, et coïncident avec la bibliothèque sur
les produits simples. Les bibliothèques de réseaux de tenseurs à charges conservées — ITensor
(tenseurs creux par blocs indexés par des nombres quantiques), TeNPy (pattes chargées, spectres
d'intrication par charge) et, pour les symétries non abéliennes, TensorKit.jl (arbres de fusion) —
fusionnent les secteurs de la même manière pour de grands systèmes et donnent l'intrication
résolue en symétrie. La bibliothèque ajoute la couche algébrique : les produits sont des algèbres
de von Neumann avec secteurs et traces, la règle de fusion est une application quelconque des
étiquettes, la trace partielle est un canal avec sa duale pour chaque trace et la restriction
correcte pour les produits fusionnés (théorèmes 7.3 et 7.6), et tout cela fonctionne par lots. Elle
n'a pas de réseaux de tenseurs, et `fused_tensor_product` prend une règle avec un seul résultat par
couple de charges, comme pour les symétries abéliennes ; une règle non abélienne comme
$j_1\otimes j_2=\bigoplus_jj$ pour $SU(2)$, avec ses coefficients de Clebsch–Gordan, n'est pas ici un
produit fusionné (la bibliothèque traite $SU(2)$ séparément pour les chaînes de spins $\frac12$,
chapitre 4).

## Exercices

1. Calculer la duale de $\operatorname{Tr}_2$ relativement à $\tau_{\mathrm{vN}}$ sur
   $\mathcal M_1\otimes\mathcal M_2$ et sur $\mathcal M_1$, et la vérifier avec
   `partial_trace_channel(alg12).adjoint('tau_vN')`.
2. (Sous-additivité.) Soit $\rho$ un état sur $\mathcal M_1\otimes\mathcal M_2$ dont les états réduits
   $\rho_1,\rho_2$ sont inversibles. Montrer que $S(\rho)\le S(\rho_1)+S(\rho_2)$, avec égalité si
   et seulement si $\rho=\rho_1\otimes\rho_2$.
3. Soit `site` l'algèbre $\mathbb C\oplus\mathbb C$ avec les charges $0,1$ (un site vide ou occupé),
   fusionnée avec elle-même par l'addition. Calculer $E(\psi)$ et ses deux parties pour
   $\psi=(|01\rangle+|10\rangle)/\sqrt2$, une particule partagée entre deux sites, et décrire la
   restriction de l'état vectoriel de $(|00\rangle+|11\rangle)/\sqrt2$ à l'algèbre fusionnée.

<details>
<summary>Solutions</summary>

1. Les poids sont $1/(C_1k_c)$ sur $\mathcal M_1$ et $1/(C_1C_2k_cl_d)$ sur le produit. Pour
   $A\in\mathcal M_1$, d'après la proposition 7.2(a),
   $\tau_{\mathrm{vN}}((\operatorname{Tr}_2X)A)=\sum_c\frac1{C_1k_c}\sum_d\operatorname{Tr}(X_{cd}(A_c\otimes1))
   =\sum_{c,d}\frac1{C_1C_2k_cl_d}\operatorname{Tr}\big(X_{cd}\,C_2l_d(A_c\otimes1)\big)$, donc la duale est
   $A\mapsto\bigoplus_{c,d}C_2l_d\,A_c\otimes1_{l_d}=A\otimes Z$, avec l'élément central
   $Z=\bigoplus_dC_2l_d1_{l_d}$ de $\mathcal M_2$. Elle n'est pas unifère : des trois traces de la
   bibliothèque, seules les traces brutes font de $\operatorname{Tr}_2$ la duale de $A\mapsto A\otimes1$ (en
   général, c'est le cas de $\operatorname{tr}_w$ sur $\mathcal M_1$ avec les poids $w_{cd}=w_c$ sur le
   produit, par le même calcul).

```python
dual_tau = composite.partial_trace_channel(alg12, keep=1).adjoint('tau_vN')
A = random_operator(alg1)
Z = alg2.central([alg2.C * l for l in alg2.k_factors])
assert close(dual_tau(A), kron(A, Z, alg12))
```

2. $\rho_1\otimes\rho_2$ est un état (théorème 7.1(c)), inversible, de vecteurs propres
   $u_i\otimes v_j$ et de valeurs propres $\lambda_i\mu_j$, donc
   $\log(\rho_1\otimes\rho_2)=\log\rho_1\otimes1+1\otimes\log\rho_2$ (théorème 1.5(a)). Par définition
   de la trace partielle, $\operatorname{Tr}\rho(\log\rho_1\otimes1)=\operatorname{Tr}\rho_1\log\rho_1$, et de même
   pour le second facteur. Ainsi $S(\rho\|\rho_1\otimes\rho_2)=S(\rho_1)+S(\rho_2)-S(\rho)$, qui est
   $\ge0$, avec égalité si et seulement si $\rho=\rho_1\otimes\rho_2$ (corollaire 3.6). Cette
   différence est l'information mutuelle.

```python
rho = states.random_density_matrix(alg12, batch_size=4)
rho1, rho2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
mutual = rho1.von_neumann_entropy() + rho2.von_neumann_entropy() - rho.von_neumann_entropy()
assert (mutual > 0).all()
assert torch.allclose(mutual, states.relative_entropy(rho, kron(rho1, rho2, alg12)))
```

3. Les secteurs sont $N=0,1,2$, de tailles $1,2,1$. $\psi$ appartient à $N=1$, avec
   $\psi_0=|01\rangle/\sqrt2$, $\psi_1=|10\rangle/\sqrt2$, donc $p=(\frac12,\frac12)$ et, les blocs
   étant de dimension un, $E(\hat\psi_c)=0$ : d'après la proposition 7.8, $E(\psi)=H(p)=\log2$ est
   entièrement de l'intrication de nombre, et la partie configurationnelle utilisable est nulle.
   Le vecteur $(|00\rangle+|11\rangle)/\sqrt2$ superpose $N=0$ et $N=2$ ; sur l'algèbre fusionnée,
   diagonale par blocs, son état vectoriel est le mélange
   $\frac12|00\rangle\langle00|+\frac12|11\rangle\langle11|$ d'états produits (chapitre 3,
   exercice 3).

```python
site = TypeIAlgebra([1, 1], [1, 1], precision='double', charges=[0, 1])
S2 = fused_tensor_product(site, site)
assert S2.k_factors == [1, 2, 1]
mat = torch.zeros(1, 3, 2, 2, dtype=dt)
mat[0, 1] = 0.5                                       # |psi><psi| with psi = (|01> + |10>)/sqrt 2
r1 = partial_trace(DensityMatrix(S2, matrix=mat), keep=1)
assert torch.allclose(r1.sector_probabilities()[0], torch.tensor([0.5, 0.5], dtype=torch.float64))
assert abs(r1.von_neumann_entropy().item() - math.log(2)) < 1e-12
```

</details>

**Suite :** le [chapitre 8](08_dynamics.md) traite de la dynamique : les équations de Schrödinger et
de Lindblad dans une algèbre à secteurs, les sauts entre secteurs et la post-sélection.
