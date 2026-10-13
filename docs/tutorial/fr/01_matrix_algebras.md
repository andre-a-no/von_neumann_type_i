[English](../01_matrix_algebras.md) · [Español](../es/01_matrix_algebras.md) · **Français** · [Deutsch](../de/01_matrix_algebras.md) · [中文](../zh/01_matrix_algebras.md) · [日本語](../ja/01_matrix_algebras.md) · [Русский](../ru/01_matrix_algebras.md)

# 1. Algèbres de matrices et calcul fonctionnel

Ce chapitre fixe les notations et démontre les propriétés de $M_n(\mathbb C)$ dont se
servent tous les chapitres suivants : le théorème spectral, le calcul fonctionnel, la
décomposition polaire et les normes définies par les valeurs singulières. Chaque résultat
est suivi du code correspondant de la bibliothèque et d'une vérification numérique.

**Vous apprendrez** à construire des opérateurs de spectre prescrit, à appliquer une
fonction à un opérateur, à manipuler des lots d'opérateurs, et à confronter un calcul à un
résultat exact.

**Conventions.** Les produits scalaires $\langle x,y\rangle=\sum_i\bar x_iy_i$ sont linéaires
par rapport au second argument. $A^*$ désigne la transposée conjuguée, de sorte que
$\langle x,Ay\rangle=\langle A^*x,y\rangle$. On note $e_1,\dots,e_n$ la base canonique, $1$
la matrice identité et $\operatorname{Tr}$ la trace usuelle.

## 1.1 L'algèbre $M_n(\mathbb C)$

$M_n(\mathbb C)$ est l'algèbre des matrices complexes $n\times n$, munie de l'involution
$A\mapsto A^*$ et de la norme d'opérateur $\|A\|=\sup_{\|x\|=1}\|Ax\|$. Un élément $A$ est
dit

| nom | définition |
|---|---|
| normal | $A^*A=AA^*$ |
| autoadjoint | $A=A^*$ |
| positif, $A\ge0$ | $\langle x,Ax\rangle\ge0$ pour tout $x$ |
| unitaire | $A^*A=AA^*=1$ |
| projecteur | $A=A^*=A^2$ |

On écrit $A\le B$ si $B-A\ge0$, et $A>0$ si $A\ge0$ et $A$ est inversible.

**Lemme 1.1 (identité C\*).** $\|A^*\|=\|A\|$ et $\|A^*A\|=\|A\|^2$.

*Démonstration.* $\|Ax\|=\sup_{\|y\|=1}|\langle y,Ax\rangle|=\sup_{\|y\|=1}|\langle A^*y,x\rangle|$,
donc $\|A\|=\sup_{\|x\|=\|y\|=1}|\langle A^*y,x\rangle|=\|A^*\|$. Ensuite,
$\|Ax\|^2=\langle x,A^*Ax\rangle\le\|A^*A\|$ pour $\|x\|=1$, d'où $\|A\|^2\le\|A^*A\|\le\|A^*\|\,\|A\|=\|A\|^2$. $\square$

> **Interprétation physique.** En mécanique quantique, un système possédant $n$
> configurations distinguables est décrit par $\mathbb C^n$. Les matrices autoadjointes
> sont les *observables* (énergie, spin, position sur un réseau) ; leurs valeurs propres
> sont les résultats possibles d'une mesure. Les projecteurs sont les questions à réponse
> oui/non. Les matrices unitaires sont les symétries et les évolutions temporelles,
> $e^{-itH}$ pour un hamiltonien $H$. Les matrices positives de trace un sont les *états*
> (matrices densité, chapitre 3).

## 1.2 Le théorème spectral

**Lemme 1.2 (Schur).** Pour tout $A\in M_n(\mathbb C)$, il existe une matrice unitaire $U$
telle que $U^*AU$ soit triangulaire supérieure.

*Démonstration.* Par récurrence sur $n$ ; le cas $n=1$ est trivial. Le polynôme
caractéristique de $A$ admet une racine $\lambda$ (théorème de d'Alembert-Gauss), donc
$Av=\lambda v$ pour un certain vecteur unitaire $v$. On complète $v$ en une base
orthonormée $v,w_2,\dots,w_n$ et l'on note $V$ la matrice unitaire ayant ces vecteurs pour
colonnes. Alors $V^*AV=\begin{pmatrix}\lambda&*\\0&A'\end{pmatrix}$ avec
$A'\in M_{n-1}(\mathbb C)$. Par hypothèse de récurrence, $U'^*A'U'$ est triangulaire
supérieure pour une certaine matrice unitaire $U'$, et $U=V\,(1\oplus U')$ convient. $\square$

**Théorème 1.3 (théorème spectral).** $A$ est normale si et seulement si
$A=U\,\mathrm{diag}(\lambda_1,\dots,\lambda_n)\,U^*$ pour une matrice unitaire $U$. Les
$\lambda_i$ sont alors les valeurs propres de $A$, comptées avec multiplicité.

*Démonstration.* Si $A=U\Lambda U^*$ avec $\Lambda$ diagonale, alors
$A^*A=U\bar\Lambda\Lambda U^*=U\Lambda\bar\Lambda U^*=AA^*$. Réciproquement, soit $A$
normale et $T=U^*AU$ triangulaire supérieure (lemme 1.2) ; $T$ est normale. Comparons les
coefficients d'indice $(1,1)$ : $(T^*T)_{11}=|t_{11}|^2$ et $(TT^*)_{11}=\sum_j|t_{1j}|^2$,
donc $t_{1j}=0$ pour $j>1$. Ainsi $T=t_{11}\oplus T'$ avec $T'$ normale et triangulaire
supérieure, et une récurrence montre que $T$ est diagonale. La diagonale de $U^*AU$ est
formée des racines de $\det(A-\lambda)=\det(T-\lambda)$. $\square$

**Corollaire 1.4.** Soit $A$ normale, de valeurs propres $\lambda_i$. Alors

1. $\|A\|=\max_i|\lambda_i|$ ;
2. $A$ est autoadjointe si et seulement si tous les $\lambda_i\in\mathbb R$ ; unitaire si et
   seulement si tous les $|\lambda_i|=1$ ; un projecteur si et seulement si tous les
   $\lambda_i\in\{0,1\}$ ;
3. $A\ge0$ si et seulement si $A$ est autoadjointe et tous les $\lambda_i\ge0$.

*Démonstration.* Avec $y=U^*x$, on a $\|Ax\|^2=\sum_i|\lambda_i|^2|y_i|^2$ et $\|y\|=\|x\|$,
ce qui donne (1), avec égalité pour $x=Ue_i$, $|\lambda_i|$ maximal. Pour (2),
$A^*=U\bar\Lambda U^*$, $A^*A=U|\Lambda|^2U^*$ et $A^2=U\Lambda^2U^*$, et les conditions
sont vérifiées si et seulement si elles le sont coefficient par coefficient pour $\Lambda$.
Pour (3), remarquons d'abord que tout $A$ positif est autoadjoint, puisque
$\langle x,(A-A^*)x\rangle=\langle x,Ax\rangle-\overline{\langle x,Ax\rangle}=0$ pour tout $x$,
et qu'une matrice $D$ telle que $\langle x,Dx\rangle=0$ pour tout $x$ est nulle (identité de
polarisation : $4\langle x,Dy\rangle=\sum_{k=0}^3 i^{-k}\langle x+i^ky,D(x+i^ky)\rangle$).
Ensuite $\langle x,Ax\rangle=\sum_i\lambda_i|y_i|^2$, quantité positive pour tout $x$ si et
seulement si tous les $\lambda_i\ge0$. $\square$

## 1.3 Le calcul fonctionnel

Soit $A=U\Lambda U^*$ normale et $f$ une fonction complexe définie sur $\operatorname{spec}A$.
On pose
$$
f(A)=U\,\mathrm{diag}\big(f(\lambda_1),\dots,f(\lambda_n)\big)\,U^*.
$$

**Théorème 1.5.** (a) $f(A)$ ne dépend pas du choix de $U$. (b) Pour un polynôme $p$,
$p(A)$ est le polynôme usuel en $A$. (c) L'application $f\mapsto f(A)$ est linéaire et
multiplicative, $\bar f(A)=f(A)^*$, et $\operatorname{spec}f(A)=f(\operatorname{spec}A)$.
(d) Si $B$ commute avec $A$, alors $B$ commute avec $f(A)$. (e) Si $g$ est définie sur
$f(\operatorname{spec}A)$, alors $g(f(A))=(g\circ f)(A)$.

*Démonstration.* (b) $A^k=U\Lambda^kU^*$, donc $p(A)=Up(\Lambda)U^*$ par linéarité. (a) Le
spectre est fini, donc l'interpolation de Lagrange fournit un polynôme $p$ égal à $f$ sur
$\operatorname{spec}A$. Alors $U\,f(\Lambda)\,U^*=U\,p(\Lambda)\,U^*=p(A)$ pour *toute*
matrice unitaire $U$ diagonalisant $A$, et $p(A)$ se définit sans recourir à $U$. (c) est
vrai pour les matrices diagonales $f(\Lambda)$ et se conserve par conjugaison par $U$ ; le
spectre de $Uf(\Lambda)U^*$ est la diagonale de $f(\Lambda)$. (d) $B$ commute avec tout
polynôme en $A$, en particulier avec $p(A)=f(A)$. (e) $f(A)=Uf(\Lambda)U^*$ est normale et
diagonalisée par la même matrice $U$, donc d'après (a) $g(f(A))=U\,g(f(\Lambda))\,U^*=(g\circ f)(A)$. $\square$

**Théorème 1.6 (racine carrée).** Pour $A\ge0$, il existe un unique $B\ge0$ tel que
$B^2=A$, à savoir $B=\sqrt A$ défini par le calcul fonctionnel.

*Démonstration.* $\sqrt A\ge0$ et $(\sqrt A)^2=A$ d'après le théorème 1.5(c) et le
corollaire 1.4(3). Soit $B\ge0$ tel que $B^2=A$, $B=V\,\mathrm{diag}(\mu_j)\,V^*$,
$\mu_j\ge0$. Alors $A=V\,\mathrm{diag}(\mu_j^2)\,V^*$, donc $\operatorname{spec}A=\{\mu_j^2\}$.
Soit $p$ un polynôme tel que $p(t)=\sqrt t$ sur $\operatorname{spec}A$. Alors
$p(A)=\sqrt A$ et $p(A)=p(B^2)=V\,\mathrm{diag}(p(\mu_j^2))\,V^*=V\,\mathrm{diag}(\mu_j)\,V^*=B$. $\square$

**Corollaire 1.7.** $A\ge0$ si et seulement si $A=T^*T$ pour une certaine matrice $T$. (Si
$A\ge0$, on prend $T=\sqrt A$ ; réciproquement, $\langle x,T^*Tx\rangle=\|Tx\|^2\ge0$.)

## 1.4 Décomposition polaire et valeurs singulières

Pour toute matrice $A$, la matrice $A^*A$ est positive (corollaire 1.7) ; sa racine carrée
$|A|=(A^*A)^{1/2}$ est la **valeur absolue** (ou module) de $A$.

**Théorème 1.8 (décomposition polaire).** Toute matrice $A\in M_n(\mathbb C)$ s'écrit
$A=W|A|$ avec $W$ unitaire.

*Démonstration.* Pour tout $x$, $\||A|x\|^2=\langle x,|A|^2x\rangle=\langle x,A^*Ax\rangle=\|Ax\|^2$.
Par conséquent $W_0:|A|x\mapsto Ax$ est une isométrie linéaire bien définie de l'image
$R_1$ de $|A|$ sur l'image $R_2$ de $A$ : si $|A|x=|A|x'$, alors
$\|A(x-x')\|=\||A|(x-x')\|=0$. En particulier $\dim R_1=\dim R_2$, de sorte que les
supplémentaires orthogonaux ont eux aussi même dimension ; soit $W_1:R_1^\perp\to R_2^\perp$
une isométrie quelconque entre eux (qui envoie une base orthonormée sur une base
orthonormée). Alors $W=W_0\oplus W_1$ est unitaire et $W|A|x=W_0|A|x=Ax$. $\square$

Les valeurs propres $s_1\ge\dots\ge s_n\ge0$ de $|A|$ sont les **valeurs singulières** de
$A$. Trois normes seront utilisées tout au long du cours :
$$
\|A\|=s_1,\qquad \|A\|_2=\Big(\sum_is_i^2\Big)^{1/2}=(\operatorname{Tr}A^*A)^{1/2},
\qquad\|A\|_1=\sum_is_i=\operatorname{Tr}|A|.
$$
La première identité est vraie car $\|Ax\|=\||A|x\|$ pour tout $x$, donc
$\|A\|=\||A|\|=s_1$ d'après le corollaire 1.4(1) ; la deuxième car
$\operatorname{Tr}A^*A=\operatorname{Tr}|A|^2=\sum_is_i^2$.

**Théorème 1.9 (dualité de la norme trace).**
$\|A\|_1=\max_{U}|\operatorname{Tr}(AU)|$ (maximum sur les matrices
unitaires $U$), et ce maximum est atteint en $U=W^*$, où $A=W|A|$. En particulier
$|\operatorname{Tr}A|\le\|A\|_1$.

*Démonstration.* $\operatorname{Tr}(AW^*)=\operatorname{Tr}(W|A|W^*)=\operatorname{Tr}|A|$. Pour
toute matrice unitaire $U$, par cyclicité de la trace et par l'inégalité de Cauchy–Schwarz
pour $\langle X,Y\rangle_{HS}=\operatorname{Tr}(X^*Y)$,
$$
|\operatorname{Tr}(AU)|=|\operatorname{Tr}(|A|^{1/2}\,UW\,|A|^{1/2})|
=|\langle|A|^{1/2},\,UW|A|^{1/2}\rangle_{HS}|
\le\||A|^{1/2}\|_2\,\|UW|A|^{1/2}\|_2=\operatorname{Tr}|A|,
$$
puisque $\|UW X\|_2=\|X\|_2$ pour $UW$ unitaire et $\||A|^{1/2}\|_2^2=\operatorname{Tr}|A|$. $\square$

En particulier, les trois quantités ci-dessus sont bien des normes. $\|\cdot\|$ est une
norme d'opérateur ; $\|\cdot\|_2$ est la norme associée au produit scalaire
$\langle X,Y\rangle_{HS}$ ; quant à $\|\cdot\|_1$, le théorème 1.9 donne
$\|A+B\|_1=\max_U|\operatorname{Tr}(AU)+\operatorname{Tr}(BU)|\le\|A\|_1+\|B\|_1$, tandis que
$\|A\|_1=0$ impose $|A|=0$ et $A=W|A|=0$. Pour $P\ge0$, on a $|P|=(P^2)^{1/2}=P$ d'après le
théorème 1.6, donc $\|P\|_1=\operatorname{Tr}P$.

## 1.5 Les opérateurs dans la bibliothèque

La bibliothèque travaille avec des algèbres $\bigoplus_cM_{k_c}(\mathbb C)$ (chapitre 2). Un
bloc unique correspond au cas $C=1$ : `TypeIAlgebra([n], [n])` est $M_n(\mathbb C)$. (Les
deux listes sont les tailles nominales et les tailles effectives des blocs ; pour l'instant
elles coïncident.)

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra

torch.manual_seed(0)
M3 = TypeIAlgebra([3], [3], precision='double')       # the algebra M_3(C), complex128

# a positive operator with eigenvalues 1, 2, 4 and Haar-random eigenvectors
A = M3.operator_from_eigenvalues(lambda d: torch.tensor([[1.0, 2.0, 4.0]]),
                                 force_positive=True)
print(A.matrix.shape)          # (batch, C, k_max, k_max) = (1, 1, 3, 3)
print(A.eigenvalues()[0])      # one tensor per block, eigenvalues in ascending order
```

Chaque opérateur est stocké sous la forme d'un unique tenseur de forme
`(batch, C, k_max, k_max)` : un lot d'opérateurs, chacun étant une liste de `C` blocs. Ici
`batch = C = 1`. La construction n'est autre que le théorème 1.3 lu à l'envers :
$A=U\Lambda U^*$ avec $U$ distribuée selon la mesure de Haar. L'option
`force_positive=True` est contrôlée sur les valeurs propres (corollaire 1.4(3)), si bien que
la bibliothèque sait que `A` est positif sans avoir à le tester de nouveau.

Le calcul fonctionnel pour les opérateurs autoadjoints est assuré par `apply_function` (qui
diagonalise chaque bloc au moyen de `eigh`). La fonction $f$ doit agir coefficient par
coefficient : elle reçoit les valeurs propres de tous les secteurs d'une même taille $k$ sous la
forme d'un seul tenseur de forme `(batch, n·k)`, où $n$ est le nombre de ces secteurs, et renvoie
un tenseur de même forme. Les fonctions usuelles disposent de leurs propres méthodes. Vérifions
les théorèmes 1.5 et 1.6 :

```python
def close(X, Y, tol=1e-10):
    """Maximal entry of |X - Y| for two operators (or tensors) is below tol."""
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

R = A.sqrt()
assert close(R @ R, A)                                  # Theorem 1.6
assert close(A.log().expm(), A)                         # exp(log A) = A, Theorem 1.5(e)
assert close(A.apply_function(lambda w: w ** 3), A @ A @ A)   # Theorem 1.5(b)
assert close(A.power(-1), A.inverse())                  # f(t) = 1/t
spec = A.apply_function(torch.exp).eigenvalues()[0]     # spec f(A) = f(spec A)
assert torch.allclose(spec, torch.exp(torch.tensor([[1.0, 2.0, 4.0]], dtype=spec.dtype)))
```

`R @ R` ne calcule encore rien : produits, sommes et fonctions construisent une recette,
qui n'est évaluée que lorsqu'on demande `.matrix` ou un nombre tel qu'une trace.

## 1.6 Normes et dualité de la norme trace

Une matrice quelconque (non normale) s'introduit par `alg.operator(tensor)`. Nous
comparons les trois normes aux valeurs singulières calculées directement par PyTorch, puis
nous vérifions la décomposition polaire (théorème 1.8) et la dualité (théorème 1.9) :

```python
G = torch.randn(1, 1, 3, 3, dtype=torch.complex128)    # Ginibre matrix: no symmetry at all
X = M3.operator(G)
s = torch.linalg.svdvals(G[0, 0])                       # singular values, descending

assert abs(X.operator_norm().item() - s[0].item()) < 1e-12
assert abs(X.frobenius_norm().item() - s.norm().item()) < 1e-12
assert abs(X.trace_norm().item() - s.sum().item()) < 1e-12

absX = X.abs()
W = X @ absX.inverse()                                  # |X| is invertible (almost surely, for Ginibre)
assert close(W.adjoint() @ W, M3.identity())            # W is unitary
assert abs((X @ W.adjoint()).trace.item() - s.sum().item()) < 1e-10    # Tr(X W*) = ||X||_1
assert abs(X.trace.item()) <= X.trace_norm().item()     # |Tr X| <= ||X||_1
```

## 1.7 Lots : de nombreux opérateurs à la fois

L'axe des lots est la raison d'être de la bibliothèque : la théorie des matrices
aléatoires, les méthodes de Monte-Carlo et l'optimisation à départs multiples requièrent
des milliers d'opérateurs, et un seul appel vectorisé est bien plus rapide qu'une boucle
Python. À titre d'exemple, vérifions un moment exact de l'**ensemble gaussien unitaire**
(GUE).

Soit $G$ une matrice à coefficients indépendants, dont les parties réelles et imaginaires
sont indépendantes et de loi $N(0,\tfrac12)$, de sorte que $\mathbb E|G_{ij}|^2=1$, et soit
$H=(G+G^*)/\sqrt2$.

**Proposition 1.10.** $\mathbb E\,\operatorname{Tr}H^4=2n^3+n$. Par conséquent, pour
$X=H/\sqrt n$, $\ \mathbb E\,\tfrac1n\operatorname{Tr}X^4=2+\tfrac1{n^2}$.

*Démonstration.* *Formule d'Isserlis (ou de Wick).* Soient $Y_1,\dots,Y_4$ des variables
aléatoires réelles centrées conjointement gaussiennes (les répétitions sont permises), et
$c_{ij}=\mathbb E[Y_iY_j]$. Leur fonction génératrice des moments conjointe est
$\mathbb E\,e^{\sum_it_iY_i}=e^{q(t)/2}$ avec $q(t)=\sum_{ij}c_{ij}t_it_j$ (fonction
génératrice des moments d'une gaussienne, résultat classique de probabilités). La dérivation sous l'espérance est légitime : la fonction
génératrice des moments est finie pour tout $t$, et pour $|t_i|\le1$ toutes les dérivées partielles
de $e^{\sum_it_iY_i}$ jusqu'à l'ordre $4$ sont majorées par $(1+\sum_i|Y_i|)^4e^{\sum_i|Y_i|}$, qui est
intégrable car $e^{\sum_i|Y_i|}\le\sum_{s\in\{\pm1\}^4}e^{\sum_is_iY_i}$ ; le théorème de convergence
dominée s'applique donc. Appliquons
$\partial^4/\partial t_1\cdots\partial t_4$ en $t=0$ : à gauche, on obtient
$\mathbb E[Y_1Y_2Y_3Y_4]$ ; à droite, seul le terme $q^2/8$ de la série exponentielle
contient $t_1t_2t_3t_4$, avec le coefficient $c_{12}c_{34}+c_{13}c_{24}+c_{14}c_{23}$.
Ainsi $\mathbb E[Y_1Y_2Y_3Y_4]$ est la somme, sur les trois appariements, des produits des
espérances des paires. Par multilinéarité, il en va de même pour des variables complexes
qui sont des combinaisons linéaires complexes de variables réelles conjointement
gaussiennes, comme les coefficients de $H$.

*Le moment.* Les espérances des paires valent
$\mathbb E[H_{ij}H_{kl}]=\tfrac12\mathbb E[(G_{ij}+\bar G_{ji})(G_{kl}+\bar G_{lk})]=\delta_{il}\delta_{jk}$,
car $\mathbb E[G_{ij}G_{kl}]=0$ et $\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}$.
Dans $\operatorname{Tr}H^4=\sum_{ijkl}H_{ij}H_{jk}H_{kl}H_{li}$, les appariements donnent :
$\{12\}\{34\}$ : $\delta_{ik}\delta_{ki}$, de somme $n^3$ (indices libres $i=k$, $j$, $l$) ;
$\{14\}\{23\}$ : $\delta_{jl}\delta_{jl}$, de somme $n^3$ ;
$\{13\}\{24\}$ : $\delta_{il}\delta_{jk}\delta_{ji}\delta_{kl}$, qui impose $i=j=k=l$, de somme $n$.
En divisant par $n\cdot n^2$, on obtient le second énoncé. $\square$

Les deux appariements *planaires* donnent le nombre de Catalan $2$ de la loi du demi-cercle ;
l'appariement croisé donne la correction en $1/n^2$.

```python
n, B = 4, 200_000
G = torch.complex(torch.randn(B, 1, n, n), torch.randn(B, 1, n, n)).to(torch.complex128) / math.sqrt(2)
M4 = TypeIAlgebra([n], [n], precision='double')
X = M4.operator((G + G.conj().transpose(-2, -1)) / math.sqrt(2 * n), is_self_adjoint=True)

X2 = X @ X
m4 = (X2 @ X2).trace.real / n                           # shape (B,): one value per sample
mean, err = m4.mean().item(), m4.std().item() / math.sqrt(B)
print(f"E tr X^4 / n = {mean:.4f} +- {err:.4f}   (exact {2 + 1 / n**2:.4f})")
assert abs(mean - (2 + 1 / n**2)) < 5 * err
```

Une vérification par Monte-Carlo est un test statistique : l'estimation doit coïncider avec
la valeur exacte à quelques écarts-types près. C'est le schéma de toutes les vérifications
par rapport à des résultats connus dans ce cours : un énoncé exact, une estimation
numérique et une tolérance explicite.

## 1.8 Points de comparaison

Pour des matrices isolées, `numpy.linalg` et `scipy.linalg` (`sqrtm`, `expm`, `logm`,
`polar`) sont les implémentations de référence, et les vérifications ci-dessus peuvent être
refaites avec elles. La bibliothèque apporte trois éléments qui compteront par la suite :
le traitement par lots (section 1.7), des opérateurs porteurs d'une structure par blocs
(chapitre 2), et la différentiation automatique à travers toutes les opérations.

## Exercices

1. Vérifier que $(fg)(A)=f(A)g(A)$ pour $f=\sin$, $g=\exp$ et $A$ autoadjointe aléatoire.
2. Montrer que $\|A\|_1=\|A^*\|_1$ et $\|AB\|_1\le\|A\|\,\|B\|_1$ à l'aide du théorème 1.9.
3. La **moyenne géométrique** de $A,B>0$ est $A\#B=A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}$.
   Démontrer que $A\#B$ est l'unique solution positive $X$ de $XA^{-1}X=B$, en déduire
   $A\#B=B\#A$, et vérifier numériquement ces deux résultats.

<details>
<summary>Solutions</summary>

1. D'après le théorème 1.5(c).

```python
Hs = M3.operator_from_eigenvalues(lambda d: torch.randn(1, d), force_self_adjoint=True)
assert close(Hs.apply_function(lambda w: torch.sin(w) * torch.exp(w)),
             Hs.apply_function(torch.sin) @ Hs.apply_function(torch.exp))
```

2. $|\operatorname{Tr}(A^*U)|=|\overline{\operatorname{Tr}(U^*A)}|=|\operatorname{Tr}(AU^*)|$, et
   $U\mapsto U^*$ est une bijection du groupe unitaire, donc les maxima du théorème 1.9
   coïncident. Pour la seconde assertion, écrivons $B=V|B|$. Pour $U$ unitaire,
   $|\operatorname{Tr}(ABU)|=|\operatorname{Tr}(|B|^{1/2}UAV|B|^{1/2})|
   \le\||B|^{1/2}\|_2\,\|UAV|B|^{1/2}\|_2\le\|A\|\operatorname{Tr}|B|$, en utilisant
   $\|YX\|_2\le\|Y\|\,\|X\|_2$ (appliquer $\|Yx\|\le\|Y\|\,\|x\|$ à chaque colonne de $X$) et
   $\|UAV\|=\|A\|$ pour $U,V$ unitaires (elles conservent la norme des vecteurs).

3. Posons $T=A^{-1/2}XA^{-1/2}$. Alors $XA^{-1}X=B$ si et seulement si
   $T^2=A^{-1/2}BA^{-1/2}$, et $X\ge0$ si et seulement si $T\ge0$, puisque
   $\langle x,Xx\rangle=\langle A^{1/2}x,TA^{1/2}x\rangle$ et que $A^{1/2}$ est inversible.
   D'après le théorème 1.6, $T=(A^{-1/2}BA^{-1/2})^{1/2}$ est l'unique solution positive,
   c'est-à-dire $X=A\#B$. Toute solution $X$ de $XA^{-1}X=B$ est inversible, car
   $\det(X)^2\det(A)^{-1}=\det B\neq0$. Ainsi $XA^{-1}X=B$ équivaut à $XB^{-1}X=A$ (à partir
   de $A^{-1}=X^{-1}BX^{-1}$, on inverse les deux membres), et l'unicité donne $A\#B=B\#A$.

```python
def geometric_mean(A, B):
    a, ai = A.sqrt(), A.power(-0.5)
    return a @ (ai @ B @ ai).sqrt() @ a

P = M3.operator_from_eigenvalues(lambda d: 0.5 + torch.rand(1, d), force_positive=True)
Q = M3.operator_from_eigenvalues(lambda d: 0.5 + torch.rand(1, d), force_positive=True)
GM = geometric_mean(P, Q)
assert close(GM, geometric_mean(Q, P), 1e-9)
assert close(GM @ P.inverse() @ GM, Q, 1e-9)
```

</details>

**Suite :** le [chapitre 2](02_von_neumann_algebras.md) remplace $M_n$ par une somme directe
d'algèbres de matrices et démontre que toute algèbre de von Neumann de dimension finie est
de cette forme.
