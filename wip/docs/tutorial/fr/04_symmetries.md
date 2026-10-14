[English](../04_symmetries.md) · [Español](../es/04_symmetries.md) · **Français** · [Deutsch](../de/04_symmetries.md) · [中文](../zh/04_symmetries.md) · [日本語](../ja/04_symmetries.md) · [Русский](../ru/04_symmetries.md)

# 4. Symétries et représentations

Le chapitre 2 a montré qu'un groupe de symétrie produit une algèbre de von Neumann, le sous-espace engendré
par les opérateurs de symétrie, ainsi que son commutant, les observables symétriques (proposition 2.7), et il
a décomposé un exemple à la main. Ce chapitre fournit la théorie des représentations sous-jacente, la
transforme en un algorithme que la bibliothèque exécute sur n'importe quelle liste de matrices, et l'applique
aux deux symétries des systèmes de spins : $U(1)$ (un nombre de particules ou une aimantation) et $SU(2)$ (les
rotations des spins).

**Vous apprendrez** la réductibilité complète et le lemme de Schur, comment la décomposition en irréductibles
détermine le commutant et l'algèbre engendrée, comment les caractères comptent les multiplicités, pourquoi un
élément aléatoire du commutant trouve numériquement la décomposition, ce que sont les secteurs de charge, les
représentations de spin irréductibles et la décomposition de $L$ spins $\frac12$, et pourquoi les hamiltoniens
invariants par $SU(2)$ ont des blocs bien plus petits que les hamiltoniens invariants par $U(1)$.

Dans tout ce chapitre, $S'$ est le commutant de $S\subseteq M_n(\mathbb C)$ (chapitre 2). Si $S^*=S$, $S''$ est
l'**algèbre de von Neumann engendrée par** $S$ : c'est une $*$-algèbre contenant $S$, et toute $*$-algèbre
$\mathcal B\supseteq S$ vérifie $\mathcal B=\mathcal B''\supseteq S''$ (théorème 2.2). $\mathcal A'_h=\{K\in\mathcal A':K=K^*\}$
est un espace vectoriel réel. $\mathbb C^k\otimes\mathbb C^m$ est muni de l'ordre de Kronecker du chapitre 7
($e_a\otimes f_j$ est le vecteur de base $am+j$, en numérotant à partir de $0$), de sorte que $1_k\otimes x$ est
diagonale par blocs avec $k$ copies de $x$. $|G|$ désigne l'ordre d'un groupe fini, et $\hbar=1$.

## 4.1 Représentations et lemme de Schur

**Définition.** Soit $X$ un ensemble muni d'une involution $x\mapsto x^*$ ($x^{**}=x$). Une
**représentation** de $X$ sur $\mathbb C^k$ est une application $\rho:X\to M_k(\mathbb C)$ telle que
$\rho(x^*)=\rho(x)^*$. Un sous-espace $K$ est **invariant** si $\rho(x)K\subseteq K$ pour tout $x$ ; $\rho$ est
**irréductible** si $k\ge1$ et si $0,\mathbb C^k$ sont les seuls sous-espaces invariants. Un
**opérateur d'entrelacement** de $\rho$ vers une représentation $\rho'$ sur $\mathbb C^{k'}$ est une matrice $T$
de taille $k'\times k$ telle que $T\rho(x)=\rho'(x)T$ pour tout $x$ ; ils forment un espace vectoriel
$\operatorname{Hom}(\rho,\rho')$, et $\operatorname{Hom}(\rho,\rho)=\rho(X)'$. $\rho$ et $\rho'$ sont **équivalentes**,
$\rho\cong\rho'$, s'il existe un opérateur d'entrelacement inversible.

Trois cas se présentent : un groupe $G$ avec $g^*=g^{-1}$, où l'on exige que $\rho$ soit un morphisme
$G\to U(k)$ (une **représentation unitaire** ; alors $\rho(g^{-1})=\rho(g)^{-1}=\rho(g)^*$) ; une $*$-algèbre munie
de son adjoint, où $\rho$ est un $*$-morphisme unifère ; et des étiquettes de générateurs avec $x^*=x$, par exemple
$\{x,y,z\}$, représentées par des familles de matrices hermitiennes telles que les matrices de spin
(section 4.5). Si $K$ est invariant et si les colonnes de $Q\in M_{k\times d}$ forment une base orthonormée de $K$,
alors $\rho(x)Q=QQ^*\rho(x)Q$, de sorte que $\rho_K(x)=Q^*\rho(x)Q$ est encore une représentation (un morphisme si
$\rho$ en est un), la **restriction** de $\rho$ à $K$.

**Lemme 4.1 (réductibilité complète).** Soit $\rho$ une représentation de $X$ sur $\mathbb C^n$.
(a) Si $K$ est invariant, $K^\perp$ l'est aussi. (b) $\mathbb C^n$ est somme directe orthogonale de sous-espaces
invariants sur chacun desquels $\rho$ est irréductible.

*Démonstration.* (a) Pour $y\in K^\perp$ et $z\in K$, $\langle z,\rho(x)y\rangle=\langle\rho(x^*)z,y\rangle=0$,
puisque $\rho(x^*)z\in K$. (b) Récurrence sur $n$. Un sous-espace invariant non nul $K_1$ de dimension minimale
est irréductible ; $K_1^\perp$ est invariant d'après (a), et les sous-espaces invariants de la restriction à
$K_1^\perp$ sont des sous-espaces invariants de $\rho$, si bien que l'hypothèse de récurrence s'y applique. $\square$

Pour un groupe fini, l'unitarité ne coûte rien : si les $\rho(g)$ sont seulement inversibles,
$P=\sum_g\rho(g)^*\rho(g)$ vérifie $P\ge\rho(e)^*\rho(e)=1$, donc $P>0$, et $\rho(h)^*P\rho(h)=\sum_g\rho(gh)^*\rho(gh)=P$
(substituer $g\to gh^{-1}$) ; par conséquent
$(P^{1/2}\rho(h)P^{-1/2})^*(P^{1/2}\rho(h)P^{-1/2})=P^{-1/2}\rho(h)^*P\rho(h)P^{-1/2}=1$, et
$h\mapsto P^{1/2}\rho(h)P^{-1/2}$ est une représentation unitaire équivalente à $\rho$ (astuce de Weyl).

**Lemme 4.2 (lemme de Schur).** Soient $\rho,\rho'$ des représentations irréductibles de $X$ sur
$\mathbb C^k,\mathbb C^{k'}$ et $T\in\operatorname{Hom}(\rho,\rho')$.
(a) $T=0$ ou $T$ est inversible. (b) $\operatorname{Hom}(\rho,\rho)=\mathbb C1$.
(c) $T^*\in\operatorname{Hom}(\rho',\rho)$, et si $T\neq0$, alors $T=\lambda U$ avec $\lambda>0$ et $U$ unitaire.
Par conséquent $\dim\operatorname{Hom}(\rho,\rho')$ vaut $1$ si $\rho\cong\rho'$ et $0$ sinon, et des représentations
irréductibles équivalentes sont unitairement équivalentes.

*Démonstration.* (a) et (b) ont été démontrés pour les groupes après la proposition 2.7, et la démonstration
n'utilise que l'invariance de $\ker T$ par $\rho$ ($T\rho(x)v=\rho'(x)Tv=0$ si $Tv=0$) et celle de l'image de $T$
par $\rho'$ : ainsi $\ker T\in\{0,\mathbb C^k\}$ et $\operatorname{ran}T\in\{0,\mathbb C^{k'}\}$ ; et pour (b),
$T-\lambda1$, où $\lambda$ est une valeur propre de $T$, est un opérateur d'entrelacement non inversible, donc
nul. (c) L'adjoint de $T\rho(x^*)=\rho'(x^*)T$ est $\rho(x)T^*=T^*\rho'(x)$. Alors $T^*T\in\operatorname{Hom}(\rho,\rho)$,
donc $T^*T=\lambda^21$ avec $\lambda>0$ si $T\neq0$ ($T^*T\ge0$ est non nul), et $U=T/\lambda$ est carrée d'après
(a), avec $U^*U=1$. Enfin, si $T_1,T_2\in\operatorname{Hom}(\rho,\rho')$ sont non nuls, alors
$T_1^{-1}T_2\in\operatorname{Hom}(\rho,\rho)=\mathbb C1$. $\square$

Si tous les $\rho(x)$ commutent entre eux (un groupe abélien), une $\rho$ irréductible est de dimension un :
chaque $\rho(x)$ appartient à $\operatorname{Hom}(\rho,\rho)=\mathbb C1$, donc tout sous-espace est invariant.

## 4.2 La décomposition isotypique et le commutant

**Théorème 4.3 (décomposition isotypique).** Soit $\rho$ une représentation de $X$ sur $\mathbb C^n$.
Il existe des représentations irréductibles deux à deux non équivalentes $\rho_1,\dots,\rho_C$ de $X$ sur
$\mathbb C^{k_c}$, des entiers $m_c\ge1$ avec $\sum_ck_cm_c=n$, et un opérateur unitaire
$V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n$ tels que
$V^*\rho(x)V=\bigoplus_c\rho_c(x)\otimes1_{m_c}$ pour tout $x$. Pour toute décomposition de ce type,
$$
V^*\rho(X)'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),\qquad
V^*\rho(X)''V=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},
$$
et $m_c=\dim\operatorname{Hom}(\rho_c,\rho)$ ; une représentation irréductible $\sigma$ de $X$ vérifie
$\operatorname{Hom}(\sigma,\rho)\neq0$ si et seulement si $\sigma\cong\rho_c$ pour exactement un $c$.

*Démonstration.* *Existence.* Le lemme 4.1 fournit des sous-espaces invariants deux à deux orthogonaux
$W_1,\dots,W_M$ de somme $\mathbb C^n$, à restrictions irréductibles $\sigma_i(x)=Q_i^*\rho(x)Q_i$. Regroupons les
$W_i$ en classes de restrictions équivalentes et choisissons dans chaque classe $c$ une référence
$\rho_c=\sigma_{i_0}$. Pour tout autre $i$ de la classe, le lemme 4.2(c) fournit un unitaire $U_i$ tel que
$\sigma_i(x)=U_i\rho_c(x)U_i^*$ ; en remplaçant $Q_i$ par $Q_iU_i$, on obtient $Q_i^*\rho(x)Q_i=\rho_c(x)$. Numérotons
les membres de la classe $c$ par $j=0,\dots,m_c-1$, notons $Q_{c,j}$ leurs bases, et posons
$V(e_a\otimes f_j)=Q_{c,j}e_a$ sur le $c$-ième terme de la somme. Comme les $W_i$ sont invariants et
orthogonaux, $Q_{c,j}^*\rho(x)Q_{d,l}=\delta_{cd}\delta_{jl}\rho_c(x)$, ce qui est la forme annoncée.

*Commutant.* Soit maintenant $V$ une décomposition quelconque comme dans l'énoncé, et soit $Q_{c,j}$ l'isométrie
$n\times k_c$ $e_a\mapsto V(e_a\otimes f_j)$ (dans le $c$-ième terme de la somme) ; pour le $V$ que l'on vient de
construire, ce sont les $Q_{c,j}$ ci-dessus. L'identité $V^*\rho(x)V=\bigoplus_c\rho_c(x)\otimes1_{m_c}$ signifie
$\rho(x)Q_{d,l}=Q_{d,l}\rho_d(x)$, et son adjointe appliquée à $x^*$ donne $Q_{c,j}^*\rho(x)=\rho_c(x)Q_{c,j}^*$. Comme
$V$ est unitaire, $Y\in M_n$ commute avec tous les $\rho(x)$ si et seulement si chaque
bloc $Y_{cj,dl}=Q_{c,j}^*YQ_{d,l}$ vérifie $Y_{cj,dl}\rho_d(x)=\rho_c(x)Y_{cj,dl}$, c'est-à-dire appartient à
$\operatorname{Hom}(\rho_d,\rho_c)$, qui vaut $0$ pour $c\neq d$ et $\mathbb C1$ pour $c=d$ (lemme 4.2). Ainsi
$Y\in\rho(X)'$ si et seulement si $Y_{cj,dl}=\delta_{cd}(y_c)_{jl}1_{k_c}$, c'est-à-dire
$V^*YV=\bigoplus_c1_{k_c}\otimes y_c$.

*Bicommutant.* Si $Z$ commute avec $\rho(X)'$, il commute avec les projecteurs sur les termes de la somme
(prendre $y_c=1$ et $y_d=0$ pour $d\neq c$), donc $V^*ZV=\bigoplus_cZ_c$, et $Z_c$ commute avec
$1_{k_c}\otimes M_{m_c}$ ; d'après le lemme 2.5, étape 5 (en échangeant les deux facteurs tensoriels),
$Z_c\in M_{k_c}\otimes1_{m_c}$. Réciproquement, ces opérateurs commutent avec $\rho(X)'$.

*Multiplicités.* Pour $T\in\operatorname{Hom}(\sigma,\rho)$ avec $\sigma$ irréductible, le même argument montre que
$Q_{c,j}^*T\in\operatorname{Hom}(\sigma,\rho_c)$, qui vaut $\mathbb Ct_{c}$ avec $t_c$ inversible si $\sigma\cong\rho_c$ et
$0$ sinon ; et $T=\sum_{c,j}Q_{c,j}Q_{c,j}^*T$. Comme les $\rho_c$ sont deux à deux non équivalentes, $\sigma$ est
équivalente à au plus l'une d'elles, disons $\rho_c$ (sinon $T=0$), et alors $T=\sum_ja_jQ_{c,j}t_c$ ;
réciproquement, chaque $Q_{c,j}t_c$ appartient à $\operatorname{Hom}(\sigma,\rho)$, et ces $m_c$ applications ont des
images orthogonales, donc $\dim\operatorname{Hom}(\sigma,\rho)=m_c$. $\square$

Pour un groupe, $\rho(G)''$ est le sous-espace engendré par $\rho(G)$ (section 2.4), et le théorème 4.3 est la
proposition 2.7 lue dans l'autre sens : là, le théorème de structure produisait les représentations
irréductibles, ici ce sont elles qui produisent la structure ; dans les deux cas, $\dim\rho(X)'=\sum_cm_c^2$ et
$\dim\rho(X)''=\sum_ck_c^2$. L'image de $P_c=\sum_jQ_{c,j}Q_{c,j}^*$, un projecteur central minimal de
$\rho(X)''$, est la **composante isotypique** de $\rho_c$.
D'après le théorème 4.3, tout $H\in\rho(X)'$ hermitien s'écrit $V(\bigoplus_c1_{k_c}\otimes h_c)V^*$ avec des
$h_c\in M_{m_c}$ hermitiennes, et $1_{k_c}\otimes h_c$ est diagonale par blocs avec $k_c$ copies de $h_c$ : le spectre de
$H$ est la réunion des spectres des blocs $h_c$ de taille $m_c\times m_c$, chaque valeur propre de $h_c$ étant répétée
$k_c$ fois.

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Si $\rho(X)$ est une symétrie
> et $H\in\rho(X)'$ l'hamiltonien, les niveaux $k_c$ fois dégénérés que l'on vient de décrire sont les **multiplets**
> de la symétrie, et les blocs qu'il faut diagonaliser ont pour tailles les multiplicités $m_c$, et non les tailles des
> multiplets.

Pour un groupe fini, les multiplicités peuvent se calculer sans rien décomposer.

**Théorème 4.4 (orthogonalité de Schur).** Soient $\rho,\rho'$ des représentations unitaires irréductibles d'un
groupe fini $G$ sur $\mathbb C^k,\mathbb C^{k'}$. Alors
$$
\frac1{|G|}\sum_{g\in G}\rho(g)_{ij}\,\overline{\rho'(g)_{i'j'}}=
\begin{cases}\delta_{ii'}\delta_{jj'}/k&\text{if }\rho=\rho',\\0&\text{if }\rho\not\cong\rho'.\end{cases}
$$

*Démonstration.* Pour une matrice $Y$ de taille $k\times k'$, posons $T=\frac1{|G|}\sum_g\rho(g)Y\rho'(g)^*$. Pour
$h\in G$, la substitution $g\to hg$ et l'égalité $\rho'(hg)^*=\rho'(g)^*\rho'(h)^{-1}$ donnent $\rho(h)T=T\rho'(h)$,
donc $T\in\operatorname{Hom}(\rho',\rho)$. Si $\rho\not\cong\rho'$, $T=0$ ; si $\rho=\rho'$, $T=\lambda1$ avec
$k\lambda=\operatorname{Tr}T=\operatorname{Tr}Y$ (lemme 4.2). Pour $Y=|j\rangle\langle j'|$, le coefficient $T_{ii'}$ est
le membre de gauche, et $\operatorname{Tr}Y=\delta_{jj'}$. $\square$

Le **caractère** d'une représentation est $\chi_\rho(g)=\operatorname{Tr}\rho(g)$ ; des représentations équivalentes
ont le même caractère, puisque $\operatorname{Tr}(T\rho(g)T^{-1})=\operatorname{Tr}\rho(g)$. On note
$\langle\varphi,\psi\rangle_G=\frac1{|G|}\sum_g\overline{\varphi(g)}\psi(g)$.

**Corollaire 4.5 (caractères).** Soit $G$ fini.
(a) Pour $\rho,\rho'$ unitaires irréductibles : $\langle\chi_\rho,\chi_{\rho'}\rangle_G$ vaut $1$ si $\rho\cong\rho'$ et $0$ sinon.
(b) Pour une représentation unitaire $\pi$ décomposée comme dans le théorème 4.3 et une $\rho$ irréductible, la
multiplicité de $\rho$ dans $\pi$ est
$$
m_\rho=\frac1{|G|}\sum_{g\in G}\overline{\chi_\rho(g)}\,\chi_\pi(g)
$$
($m_\rho=m_c$ si $\rho\cong\rho_c$, et $m_\rho=0$ si $\rho$ n'est équivalente à aucune $\rho_c$).
(c) $\dim\pi(G)'=\sum_cm_c^2=\frac1{|G|}\sum_{g\in G}|\chi_\pi(g)|^2$. En particulier, $\pi$ est irréductible si et
seulement si ce nombre vaut $1$.

*Démonstration.* (a) Remplaçons $\rho'$ par une représentation équivalente (de même caractère). Sommons le
théorème 4.4 sur $i=j$, $i'=j'$ et conjuguons : on obtient $\sum_{i,i'}\delta_{ii'}/k=1$ si $\rho=\rho'$, et $0$ si
$\rho\not\cong\rho'$. (b) La trace est unitairement invariante et $\operatorname{Tr}(a\otimes1_m)=m\operatorname{Tr}a$, donc
$\chi_\pi=\sum_cm_c\chi_{\rho_c}$ ; on applique (a). (c) Par le même développement et (a),
$\langle\chi_\pi,\chi_\pi\rangle_G=\sum_cm_c^2$, qui est $\dim\pi(G)'$ d'après le théorème 4.3. $\square$

Une seconde démonstration de (c) correspond à la manière dont la bibliothèque échantillonne les commutants. Le
**twirl** $P(Y)=\frac1{|G|}\sum_g\pi(g)Y\pi(g)^*$ est idempotent, d'image $\pi(G)'$ (chapitre 2, exercice 3), et
autoadjoint pour $\langle Y,Z\rangle_{HS}=\operatorname{Tr}(Y^*Z)$ (substituer $g\to g^{-1}$) ; c'est donc le projecteur
orthogonal sur $\pi(G)'$, et $\dim\pi(G)'$ est sa trace. L'application $Y\mapsto AYB$ a pour trace
$\operatorname{Tr}A\operatorname{Tr}B$ (ses coefficients diagonaux sont les $A_{ii}B_{jj}$), ce qui donne
$\frac1{|G|}\sum_g|\chi_\pi(g)|^2$.

Reprenons l'exemple de la section 2.4 à l'aide des caractères. Une matrice de permutation a pour trace le nombre
de vecteurs de base qu'elle fixe, et $\pi(\sigma)$ fixe $|b_1b_2b_3\rangle$ si et seulement si $b$ est constant sur
les cycles de $\sigma$ (lemme B.1(c)), donc $\chi_\pi(\sigma)=2^{\#\text{cycles}(\sigma)}$. Trois représentations unitaires
de $S_3$ sont à portée de main. Soit $P_\sigma e_i=e_{\sigma(i)}$ la représentation par permutations sur $\mathbb C^3$ (un
homomorphisme à valeurs dans les matrices de permutation). La représentation **triviale** est $\sigma\mapsto1$. La
**signature** est $\operatorname{sgn}\sigma=\det P_\sigma$, un homomorphisme à valeurs dans $\{\pm1\}\subseteq U(1)$ car le
déterminant est multiplicatif et le déterminant d'une matrice de permutation vaut $\pm1$. La droite $\mathbb C(1,1,1)$
est invariante, donc $E=(1,1,1)^\perp$ l'est aussi (lemme 4.1(a)) ; la représentation **standard** est la restriction
de $P$ à $E$. Dans une base orthonormée adaptée à $\mathbb C(1,1,1)\oplus E$, $P_\sigma=1\oplus\operatorname{std}(\sigma)$,
donc $\chi_{\rm std}(\sigma)=\operatorname{Tr}P_\sigma-1$, le nombre de points fixes moins un. Les six éléments de $S_3$
sont l'identité, trois transpositions et deux $3$-cycles, avec $3,2,1$ cycles et $3,1,0$ points fixes ; une
transposition échange deux colonnes de $1_3$, donc sa signature vaut $-1$, et un $3$-cycle est un produit de deux
transpositions ($(1\,2\,3)=(1\,3)(1\,2)$), donc sa signature vaut $+1$. On obtient la table

| | identité ($1$) | transpositions ($3$) | $3$-cycles ($2$) |
|---|---|---|---|
| $\chi_{\rm triv}$ | $1$ | $1$ | $1$ |
| $\operatorname{sgn}$ | $1$ | $-1$ | $1$ |
| $\chi_{\rm std}$ | $2$ | $0$ | $-1$ |
| $\chi_\pi$ | $8$ | $4$ | $2$ |

(ainsi $\operatorname{sgn}\sigma=(-1)^{3-\#\text{cycles}(\sigma)}$ ; tous les caractères sont réels). Avec
$\langle\varphi,\psi\rangle_G$ sommé classe par classe :
$\langle\chi_{\rm triv},\chi_{\rm triv}\rangle_G=\frac{1+3+2}6=1$, $\langle\operatorname{sgn},\operatorname{sgn}\rangle_G=\frac{1+3+2}6=1$ et
$\langle\chi_{\rm std},\chi_{\rm std}\rangle_G=\frac{4+0+2}6=1$, donc les trois représentations sont irréductibles
d'après le corollaire 4.5(c), et elles sont deux à deux non équivalentes car leurs caractères diffèrent. Le
corollaire 4.5(b) donne les multiplicités dans $\pi$ :
$$
m_{\rm triv}=\tfrac{8+3\cdot4+2\cdot2}6=4,\qquad m_{\rm sgn}=\tfrac{8-3\cdot4+2\cdot2}6=0,\qquad
m_{\rm std}=\tfrac{2\cdot8+0-2\cdot2}6=2 .
$$
Ainsi la représentation triviale et la représentation standard sont équivalentes à des constituants $\rho_c$ de $\pi$
(théorème 4.3) avec $m_c=4$ et $m_c=2$ ; elles rendent compte de $4\cdot1+2\cdot2=8=\sum_ck_cm_c$ dimensions, il n'y a
donc pas d'autre constituant : $k=(1,2)$, $m=(4,2)$, et $\dim\pi(S_3)'=4^2+2^2=\frac16(64+3\cdot16+2\cdot4)=20$
(corollaire 4.5(c)). Le code recalcule la table : il obtient $P_\sigma$ comme restriction de $\pi(\sigma)$ au
sous-espace invariant engendré par $|100\rangle,|010\rangle,|001\rangle$ (sur lequel $\pi(\sigma)$ permute les positions
de l'unique $1$, lemme B.1(a)), et $\#\text{cycles}(\sigma)=\dim\ker(P_\sigma-1)$ (les vecteurs constants sur les
cycles).

```python
import math
import torch
import torch_vn_algebra.decompose as dc

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    return (torch.as_tensor(X) - torch.as_tensor(Y)).abs().max().item() < tol

G = dc.permutation_representation(3, 2)                  # (6, 8, 8): pi(sigma) permutes three qubits
chi_pi = G.diagonal(dim1=-2, dim2=-1).sum(-1).real       # characters Tr pi(sigma)
P3 = G[:, [4, 2, 1]][:, :, [4, 2, 1]].real              # on |100>, |010>, |001>: S_3 permuting C^3
cycles = 3 - torch.linalg.matrix_rank(P3 - torch.eye(3, dtype=torch.float64)).double()   # dim ker(P - 1)
assert close(chi_pi, 2.0 ** cycles)
irreps = {'trivial': torch.ones(6, dtype=torch.float64), 'sign': (-1.0) ** (3 - cycles),
          'standard': P3.diagonal(dim1=-2, dim2=-1).sum(-1) - 1}                 # fixed points - 1
chis = torch.stack(list(irreps.values()))
assert close(chis @ chis.T / 6, torch.eye(3, dtype=torch.float64))    # Corollary 4.5(a): all irreducible
mult = chis @ chi_pi / 6                                               # Corollary 4.5(b)
print(dict(zip(irreps, mult.tolist())), (chi_pi ** 2).mean().item())
assert close(mult, torch.tensor([4.0, 0.0, 2.0], dtype=torch.float64))
assert dc.commutant_dimension(G) == 20                                # Corollary 4.5(c): M_4 (+) M_2
```

C'est la décomposition de la section 2.4, obtenue cette fois à partir d'une table de six nombres par représentation.

## 4.3 Des générateurs à l'algèbre : un algorithme randomisé

Le théorème 4.3 décrit la décomposition mais ne la calcule pas : sa démonstration a besoin des sous-espaces
invariants irréductibles. Murota, Kanno, Kojima et Kojima (*Japan J. Indust. Appl. Math.* 27 (2010) 125–160 ;
voir aussi Maehara et Murota, *ibid.* 263–293) ont observé que les sous-espaces propres d'un *élément aléatoire du
commutant* les trouvent. Soit $\mathcal A$ l'algèbre de von Neumann engendrée par des matrices données
$S_1,\dots,S_J$ et leurs adjointes, avec une décomposition $V^*\mathcal AV=\bigoplus_cM_{k_c}\otimes1_{m_c}$,
$V^*\mathcal A'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$ (théorème 2.6) ; seule l'analyse utilise $V$.

**Lemme 4.6.** Soit $K\in\mathcal A'_h$ avec $V^*KV=\bigoplus_c1_{k_c}\otimes\kappa_c$. Tout sous-espace propre de
$K$ est $\mathcal A$-invariant. Si les $N=\sum_cm_c$ valeurs propres de l'ensemble des $\kappa_c$ sont deux à deux
distinctes, les sous-espaces propres de $K$ sont les $N$ sous-espaces $V(\mathbb C^{k_c}\otimes e)$ (dans le
$c$-ième terme de la somme directe, $e$ étant un vecteur propre de $\kappa_c$), et $\mathcal A$ agit irréductiblement sur chacun
d'eux, comme $M_{k_c}$.

*Démonstration.* Si $Kv=\lambda v$ et $a\in\mathcal A$, alors $Kav=aKv=\lambda av$. Écrivons
$\kappa_c=\sum_j\mu_{cj}|e_{cj}\rangle\langle e_{cj}|$ avec des $e_{cj}$ orthonormés (théorème 1.3). Alors
$V^*KV=\sum_{c,j}\mu_{cj}R_{cj}$, où $R_{cj}$ est le projecteur sur $\mathbb C^{k_c}\otimes e_{cj}$ ; les $R_{cj}$
sont orthogonaux de somme $1$ et les $\mu_{cj}$ distincts, donc ce sont les projecteurs spectraux. Sur
$\mathbb C^{k_c}\otimes e_{cj}$, l'élément $V^*aV=\bigoplus_da_d\otimes1$ agit comme $a_c$, et les $a_c$ parcourent
$M_{k_c}$, qui n'a pas d'autres sous-espaces invariants que $0$ et $\mathbb C^{k_c}$. $\square$

**Proposition 4.7 (généricité).** Soient $K$ et $Y$ des éléments aléatoires indépendants de $\mathcal A'_h$ et de
$\mathcal A'$ dont les lois ont des densités par rapport à la mesure de Lebesgue sur ces espaces vectoriels réels
(lemme A.3). Presque sûrement : (a) l'hypothèse du lemme 4.6 est satisfaite ; (b) $\Pi_iY\Pi_j\neq0$ pour deux
projecteurs spectraux distincts quelconques $\Pi_i,\Pi_j$ de $K$ dont les images sont dans le même terme de la somme
directe.

*Démonstration.* (a) La matrice hermitienne $\kappa=\bigoplus_c\kappa_c\in M_N$ dépend de $K$ de façon $\mathbb R$-linéaire.
Le déterminant de Gram $\delta(K)=\det[\operatorname{Tr}\kappa^{i+j}]_{i,j=0}^{N-1}$ de $1,\kappa,\dots,\kappa^{N-1}$ pour
$\langle\cdot,\cdot\rangle_{HS}$ est un polynôme réel en les coordonnées de $K$, non nul si et seulement si ces
puissances sont linéairement indépendantes. Comme $p(\kappa)=U\operatorname{diag}(p(\lambda_l))U^*$ pour
$\kappa=U\operatorname{diag}(\lambda)U^*$, un polynôme non nul de degré $<N$ annule $\kappa$ si et seulement s'il
s'annule sur $\operatorname{spec}\kappa$, ce qui est possible si et seulement si $\kappa$ a moins de $N$ valeurs
propres distinctes (prendre $\prod_\mu(t-\mu)$ sur celles-ci). Ainsi $\delta(K)\neq0$ est l'hypothèse du lemme 4.6,
et $\delta\not\equiv0$ : il est non nul en $V(\bigoplus_c1\otimes\operatorname{diag}(\text{réels distincts}))V^*\in\mathcal A'_h$.
Donc $\mathbb P(\delta(K)=0)=0$ d'après le corollaire A.7(b).

(b) L'ensemble $O=\{\delta\neq0\}\subseteq\mathcal A'_h$ est ouvert ($\delta$ est continue). Pour $K\in O$, $K$ a exactement
$N$ valeurs propres distinctes $\mu_1(K)<\dots<\mu_N(K)$ (celles de $\kappa$, de multiplicités $k_c$), donc ses projecteurs
spectraux $\Pi_1(K),\dots,\Pi_N(K)$, numérotés par valeur propre croissante, sont continus sur $O$ (proposition A.13).
Fixons $l\neq l'$ et soit $\mathcal N_{ll'}$ l'ensemble des $(K,Y)\in O\times\mathcal A'$ tels que $P_c\Pi_l(K)=\Pi_l(K)$ et
$P_c\Pi_{l'}(K)=\Pi_{l'}(K)$ pour un certain $c$ ($P_c$ les projecteurs centraux) et $\Pi_l(K)Y\Pi_{l'}(K)=0$ ; il est
défini par des équations continues sur le borélien $O\times\mathcal A'$, donc borélien. Fixons $K\in O$ pour lequel les
deux images sont dans le terme $c$. D'après le lemme 4.6, $V^*\Pi_lV$ et $V^*\Pi_{l'}V$ projettent sur
$\mathbb C^{k_c}\otimes e$ et $\mathbb C^{k_c}\otimes e'$, avec des vecteurs propres orthonormés $e,e'$ de $\kappa_c$, et pour
$V^*YV=\bigoplus_d1\otimes y_d$, $\Pi_lY\Pi_{l'}=0$ si et seulement si $\langle e,y_ce'\rangle=0$. L'application
$Y\mapsto\langle e,y_ce'\rangle$ est $\mathbb R$-linéaire de $\mathcal A'$ dans $\mathbb C$ et non nulle (sa valeur en
$V(0\oplus\dots\oplus1\otimes|e\rangle\langle e'|\oplus\dots\oplus0)V^*$ vaut $1$), donc son noyau est un sous-espace
propre, un ensemble négligeable (corollaire A.7(c)). Ainsi toute section $\{Y:(K,Y)\in\mathcal N_{ll'}\}$ est négligeable
(elle est vide si les images sont dans des termes différents) et, comme $Y$ est indépendant de $K$,
$\mathbb P((K,Y)\in\mathcal N_{ll'})=0$ d'après le corollaire A.7(d). L'événement de (b) n'est en défaut que sur
$\{\delta(K)=0\}\cup\bigcup_{l\neq l'}\{(K,Y)\in\mathcal N_{ll'}\}$, une réunion finie d'événements de probabilité nulle.
$\square$

**Lemme 4.8 (compressions ; alignement des copies).** (a) Soit $W$ un sous-espace $\mathcal A$-invariant de base
orthonormée $Q\in M_{n\times d}$ et $\sigma(a)=Q^*aQ$ la restriction. Alors $\sigma(\mathcal A)$ est une $*$-algèbre,
$\{Q^*YQ:Y\in\mathcal A'\}=\sigma(\mathcal A)'$ et $\{Q^*KQ:K\in\mathcal A'_h\}=\sigma(\mathcal A)'_h$. $\mathcal A$ agit
irréductiblement sur $W$ si et seulement si $Q^*YQ\in\mathbb C1$ pour tout $Y\in\mathcal A'$. Sinon, si $X$ est un élément
aléatoire de $\mathcal A'$ ou de $\mathcal A'_h$ à densité, alors $Q^*XQ\notin\mathbb C1$ presque sûrement.

(b) Soient $W_i,W_j$ des sous-espaces invariants sur lesquels $\mathcal A$ agit
irréductiblement, de bases orthonormées $Q_i,Q_j$ et de restrictions $\sigma_i(a)=Q_i^*aQ_i$, et soit
$Y\in\mathcal A'$. Alors $T=Q_i^*YQ_j\in\operatorname{Hom}(\sigma_j,\sigma_i)$. Si $T\neq0$, alors $T=\lambda U$ avec
$\lambda>0$, $U$ unitaire, et $Q_j'=Q_jU^*$ vérifie $Q_j'^*aQ_j'=\sigma_i(a)$ pour tout $a\in\mathcal A$. Si
$W_i,W_j$ sont dans des termes différents ($W_i\subseteq P_c\mathbb C^n$, $W_j\subseteq P_d\mathbb C^n$, $c\neq d$),
alors $T=0$.

*Démonstration.* (a) Comme dans la démonstration du théorème 4.3 ($W$ est invariant par $a$ et par $a^*$),
$aQ=Q\sigma(a)$ et $Q^*a=\sigma(a)Q^*$ pour $a\in\mathcal A$. Donc $\sigma(ab)=Q^*aQ\sigma(b)=\sigma(a)\sigma(b)$,
$\sigma(a^*)=\sigma(a)^*$ et $\sigma(1)=Q^*Q=1_d$, si bien que $\sigma$ est un $*$-homomorphisme unifère et $\sigma(\mathcal A)$
une $*$-algèbre ; et pour $Y\in\mathcal A'$, $Q^*YQ\sigma(a)=Q^*YaQ=Q^*aYQ=\sigma(a)Q^*YQ$. Réciproquement, pour
$T\in\sigma(\mathcal A)'$ posons $Y=QTQ^*$ ; alors $aY=Q\sigma(a)TQ^*=QT\sigma(a)Q^*=Ya$, donc $Y\in\mathcal A'$, avec
$Q^*YQ=T$, et $Y$ est hermitien si $T$ l'est. Si $\sigma$ est irréductible,
$\sigma(\mathcal A)'=\operatorname{Hom}(\sigma,\sigma)=\mathbb C1$ (lemme 4.2(b)). Si $E\subseteq\mathbb C^d$ est un
sous-espace $\sigma$-invariant autre que $0,\mathbb C^d$, $E^\perp$ l'est aussi (lemme 4.1(a)), et le projecteur
orthogonal sur $E$ est un élément non scalaire de $\sigma(\mathcal A)'_h$. Dans ce cas, $\mathbb C1\cap\sigma(\mathcal A)'$
et $\mathbb R1\cap\sigma(\mathcal A)'_h$ sont des sous-espaces réels propres, leurs images réciproques par les applications
$\mathbb R$-linéaires surjectives $X\mapsto Q^*XQ$ sont des sous-espaces propres de $\mathcal A'$ et de $\mathcal A'_h$, et
ceux-ci sont négligeables (corollaire A.7(c)).

(b) Comme en (a), $aQ_j=Q_j\sigma_j(a)$ et $Q_i^*a=\sigma_i(a)Q_i^*$,
donc $T\sigma_j(a)=Q_i^*YaQ_j=Q_i^*aYQ_j=\sigma_i(a)T$. Le lemme 4.2 donne $T=\lambda U$ et
$U\sigma_j(a)U^*=\sigma_i(a)$, qui vaut $Q_j'^*aQ_j'$. Enfin $P_c\in\mathcal A\cap\mathcal A'$, $Q_i=P_cQ_i$ et
$P_cQ_j=0$, donc $T=Q_i^*P_cYQ_j=Q_i^*YP_cQ_j=0$. $\square$

L'algorithme (`decompose.wedderburn`, avec les interfaces `from_group` et `from_generators`) est alors le suivant :

1. *Échantillonner le commutant.* Pour tous les éléments d'un groupe fini, $Y=P(Z)$ et $K=P(Z'+Z'^*)$ avec le twirl
   $P$ et des matrices de Ginibre indépendantes $Z,Z'$. $P$ est une projection linéaire de $M_n$ sur $\mathcal A'$ avec
   $P(Z)^*=P(Z^*)$, donc $Z\mapsto P(Z)$ et $Z'\mapsto P(Z'+Z'^*)$ sont des surjections $\mathbb R$-linéaires sur
   $\mathcal A'$ et sur $\mathcal A'_h$ (pour $H\in\mathcal A'_h$, $P(\frac H2+\frac H2)=H$) ; une matrice de Ginibre est un
   vecteur gaussien centré non dégénéré de $M_n$ (proposition A.16), donc $Y$ et $K$ sont des vecteurs gaussiens centrés
   non dégénérés de $\mathcal A'$ et de $\mathcal A'_h$ et ont des densités (théorème A.15). Pour des générateurs,
   $\mathcal A'=\{X:[S_j,X]=[S_j^*,X]=0\}$ est le noyau d'une application linéaire de $M_n$ (la bibliothèque utilise les
   parties hermitiennes $\frac12(S_j+S_j^*)$, $\frac1{2i}(S_j-S_j^*)$, qui engendrent le même espace complexe que
   $S_j,S_j^*$ et ont donc le même commutant) ; avec une base $(B_r)_{r\le D}$ de celui-ci, $Y=\sum_rc_rB_r$ et
   $K=Y'+Y'^*$, où $Y'$ est une copie indépendante de $Y$ et les $c_r$ des gaussiennes complexes standard
   indépendantes : ce sont des images d'un vecteur gaussien non dégénéré de $\mathbb C^D$ par la bijection
   $\mathbb R$-linéaire $c\mapsto\sum_rc_rB_r$ sur $\mathcal A'$ et par la surjection $\mathbb R$-linéaire $c\mapsto X+X^*$
   ($X=\sum_rc_rB_r$) sur $\mathcal A'_h$, qui ont donc là encore des densités. Chaque élément aléatoire est un nouveau
   tirage, indépendant des précédents ; la bibliothèque les renormalise à norme un, ce qui ne change aucun sous-espace
   propre ni l'annulation d'aucun bloc.
2. *Scinder.* Les sous-espaces propres de $K$ sont irréductibles (lemme 4.6, proposition 4.7(a)) ; leurs dimensions
   sont les $k_c$.
3. *Regrouper et aligner.* Les blocs $Q_i^*YQ_j$ entre sous-espaces propres de même dimension sont non nuls
   exactement à l'intérieur d'un même terme (proposition 4.7(b), lemme 4.8(b)) ; les classes de la relation « reliés
   par un bloc non nul » donnent les secteurs, leurs tailles les $m_c$, et les facteurs polaires alignent les copies sur
   une copie de référence par classe (lemme 4.8(b) ; une copie alignée sur une copie alignée est alignée sur la
   référence).
4. *Assembler* $V$ dans l'ordre de Kronecker, de sorte que $V^*aV=\bigoplus_ca_c\otimes1_{m_c}$.

Ainsi, en arithmétique exacte, les sous-espaces propres de $K$ étant séparés exactement, l'algorithme renvoie presque
sûrement une décomposition comme dans le théorème 4.3 : d'après le lemme 4.6, le terme $c$ contient $m_c$ sous-espaces
propres de dimension $k_c$, et après alignement la construction est la partie existence de la démonstration du
théorème 4.3. En virgule flottante, la bibliothèque traite comme égales des valeurs propres plus proches qu'une
tolérance, si bien qu'un « sous-espace propre » calculé $W$ peut être une somme de plusieurs sous-espaces propres de
$K$ ; il reste $\mathcal A$-invariant (lemme 4.6), et la bibliothèque y comprime un nouvel élément aléatoire
$X\in\mathcal A'_h$. D'après le lemme 4.8(a), $Q^*XQ$ est un multiple de $1$ dès que $\mathcal A$ agit irréductiblement sur
$W$, et sinon presque sûrement pas ; ses sous-espaces propres, invariants par $\sigma(\mathcal A)$ (lemme 4.6 appliqué à
$\sigma(\mathcal A)$), scindent alors $W$ en sous-espaces invariants plus petits, et la bibliothèque recommence jusqu'à ce
que toute compression soit scalaire. Quoi qu'il arrive en virgule flottante, le résultat est certifié à la fin :

**Proposition 4.9 (certificat).** Soit $\mathcal A$ engendrée par $S_1,\dots,S_J$ et leurs adjointes, et
$\mathcal B=V(\bigoplus_cM_{k_c}\otimes1_{m_c})V^*$ avec $V$ unitaire. Si tous les $S_j\in\mathcal B$ et
$\dim\mathcal A'=\sum_cm_c^2$, alors $\mathcal A=\mathcal B$.

*Démonstration.* $\mathcal B$ est une $*$-algèbre contenant tous les $S_j,S_j^*$, donc $\mathcal A\subseteq\mathcal B$
et $\mathcal B'\subseteq\mathcal A'$. Le théorème 4.3, appliqué à la représentation de $\mathcal B$ sur $\mathbb C^n$
(les applications de bloc $b\mapsto b_c$ sont irréductibles, et non équivalentes, car le $c$-ième projecteur central
est envoyé sur $1$ par l'une et sur $0$ par les autres), donne $\dim\mathcal B'=\sum_cm_c^2=\dim\mathcal A'$. Donc
$\mathcal B'=\mathcal A'$ et $\mathcal B=\mathcal B''=\mathcal A''=\mathcal A$. $\square$

`Decomposition.verify()` renvoie les distances relatives des générateurs (et de nouveaux éléments du commutant) à
$\mathcal B$ (resp. $\mathcal B'$), et compare $\sum_cm_c^2$ à $\dim\mathcal A'$ (obtenue à partir des caractères pour
un groupe, corollaire 4.5(c) ; à partir du noyau pour des générateurs). Nous l'exécutons sur $S_3$ et trois qubits,
puis nous regardons à l'intérieur : les sous-espaces propres d'un $K$ aléatoire et les opérateurs d'entrelacement
d'un $Y$ aléatoire.

```python
D = dc.from_group(G)                                     # the twirl samples pi(S_3)'
print(D)
assert D.k_factors == (1, 2) and D.multiplicities == (4, 2)            # Section 2.4: k = (1, 2), m = (4, 2)
assert D.verify()['ok']
rho = D.to_algebra(G).matrix                             # (6, 2, 2, 2): the blocks rho_c(sigma), padded
assert close(rho[:, 0, 0, 0], torch.ones(6))                           # sector 1: trivial representation
assert close(rho[:, 1].diagonal(dim1=-2, dim2=-1).sum(-1), chis[2])    # sector 2: the standard one
schur = torch.einsum('gij,gkl->ijkl', rho[:, 1], rho[:, 1].conj()) / 6
assert close(schur.reshape(4, 4), torch.eye(4) / 2)                    # Theorem 4.4 with k = 2
assert close(D.from_algebra(D.to_algebra(G)), G)

Z = torch.randn(8, 8, dtype=dt)
K = dc.twirl(G, Z + Z.mH)                                # a random element of pi(S_3)'_h
w, Q = torch.linalg.eigh(K)
cuts = [0] + [i + 1 for i in range(7) if w[i + 1] - w[i] > 1e-8] + [8]
spaces = [Q[:, a:b] for a, b in zip(cuts[:-1], cuts[1:])]
assert sorted(W.shape[1] for W in spaces) == [1, 1, 1, 1, 2, 2]      # Lemma 4.6: m_c pieces of size k_c
assert all(close(G @ W, W @ (W.mH @ G @ W)) for W in spaces)           # invariant subspaces
Y = dc.twirl(G, torch.randn(8, 8, dtype=dt))
Wa, Wb = [W for W in spaces if W.shape[1] == 2]
s = torch.linalg.svdvals(Wa.mH @ Y @ Wb)                               # Lemma 4.8(b): a multiple of a unitary
assert s[1] > 1e-3 and abs(s[0] - s[1]) < 1e-10
assert all(close(W.mH @ Y @ Wa, 0) for W in spaces if W.shape[1] == 1)   # different summands: zero
```

Le commutant domine le coût : un échantillon du twirl coûte $|G|$ produits de matrices $n\times n$. Pour des
générateurs, la bibliothèque diagonalise d'abord un élément hermitien $R=U\operatorname{diag}(\lambda)U^*$ de $\mathcal A$
(une combinaison réelle aléatoire des parties hermitiennes des générateurs, plus un multiple d'un produit symétrisé de
deux telles combinaisons). Tout $X\in\mathcal A'$ commute avec $R$, donc avec ses projecteurs spectraux (théorème
1.5(d)), si bien que $U^*XU$ est diagonale par blocs relativement aux sous-espaces propres de $R$, de dimensions
$d_\lambda$ ; on l'impose d'avance, et le problème de noyau a $N_R=\sum_\lambda d_\lambda^2$ inconnues au lieu de $n^2$.
Regrouper plusieurs valeurs propres ne fait qu'agrandir les blocs, donc le problème réduit a le même ensemble de
solutions $\mathcal A'$. Dans la décomposition, $V^*RV=\bigoplus_cr_c\otimes1_{m_c}$ avec des $r_c\in M_{k_c}$
hermitiennes ; si $\lambda$ a la multiplicité $\mu_c(\lambda)$ dans $r_c$, alors $d_\lambda=\sum_c\mu_c(\lambda)m_c$, et comme
$\mu^2\ge\mu$ pour les entiers $\mu\ge0$,
$$
N_R=\sum_\lambda\Big(\sum_c\mu_c(\lambda)m_c\Big)^2\ge\sum_{\lambda,c}\mu_c(\lambda)^2m_c^2\ge\sum_{\lambda,c}\mu_c(\lambda)m_c^2=\sum_ck_cm_c^2,
$$
avec égalité si et seulement si les valeurs propres de toutes les $r_c$ sont simples et deux à deux distinctes. Pour le
spin total de $L=3,4,5$ spins, le code ci-dessous trouve l'égalité (observation numérique). Même le minimum est
grand : pour huit spins $\frac12$ ($n=256$), $\sum_ck_cm_c^2=4900$, et le système réduit, une matrice
$(n^2+N_R)\times N_R$, a $3.5\cdot10^8$ coefficients complexes. C'est pourquoi, pour les chaînes de spins, la
section 4.6 utilise le module `su2`, qui connaît la réponse à l'avance.

```python
for n in (3, 4, 5):
    Dn = dc.from_generators(dc.total_spin(n))
    _, info = dc.commutant_basis(dc.total_spin(n), return_info=True)
    N_min = sum(k * m * m for k, m in zip(Dn.k_factors, Dn.multiplicities))
    assert info['unknowns'] >= N_min                                  # N_R >= sum_c k_c m_c^2
    print(n, info['unknowns'], N_min)                                # observed: equality
```

## 4.4 $U(1)$ et secteurs de charge

**Proposition 4.10.** Soit $Q$ une matrice hermitienne sur $\mathbb C^n$ à valeurs propres entières $q$ (les
**charges**), de projecteurs spectraux $P_q$, avec $d_q=\operatorname{rank}P_q$, et soit $U_\theta=e^{i\theta Q}$ ; on
note $U(1)'=\{U_\theta:\theta\in\mathbb R\}'$ et $U(1)''$ son commutant.

1. $\theta\mapsto U_\theta$ est une représentation unitaire de $U(1)=\mathbb R/2\pi\mathbb Z$, et toute représentation
   unitaire continue de $U(1)$ sur $\mathbb C^n$ est de cette forme.
2. Pour $X\in M_n$, les assertions suivantes sont équivalentes : (i) $XU_\theta=U_\theta X$ pour tout $\theta$ ;
   (ii) $[X,Q]=0$ ; (iii) $X=\sum_qP_qXP_q$. Par conséquent $U(1)'\cong\bigoplus_qM_{d_q}(\mathbb C)$ (les opérateurs
   diagonaux par blocs) et $U(1)''=\operatorname{span}\{P_q\}$ : tous les $k_c$ valent $1$, et les multiplicités sont
   les $d_q$.
3. (*Charge conservée.*) Pour $H$ hermitien, $e^{itH}Qe^{-itH}=Q$ pour tout $t$ si et seulement si $[H,Q]=0$,
   c'est-à-dire si et seulement si $H\in U(1)'$. Alors chaque poids de secteur $\operatorname{Tr}(P_q\rho(t))$ de
   $\rho(t)=e^{-itH}\rho e^{itH}$ est constant.

*Démonstration.* (1) $U_\theta=\sum_qe^{iq\theta}P_q$ est unitaire, $U_\theta U_{\theta'}=U_{\theta+\theta'}$
(théorème 1.5(c)) et $U_{2\pi}=1$. Réciproquement, pour une représentation continue, $(U_\theta)_{\theta\ge0}$ est un
semi-groupe continu de matrices, donc $U_\theta=e^{\theta B}$ pour $\theta\ge0$ d'après le lemme 8.5 (démontré au
chapitre 8 ; la démonstration n'utilise que la dimension finie), et pour $\theta<0$,
$U_\theta=U_{-\theta}^{-1}=e^{\theta B}$. En dérivant $U_\theta^*U_\theta=1$ en $0$ (lemme 8.1), on obtient $B^*+B=0$,
donc $B=iQ$ avec $Q$ hermitienne. La série exponentielle $e^{i\theta Q}$ coïncide avec le calcul fonctionnel
(proposition 8.2(1) avec $H=-Q$), donc $U_{2\pi}=\sum_qe^{2\pi iq}P_q=1$ impose $e^{2\pi iq}=1$, $q\in\mathbb Z$, pour toute
valeur propre $q$ de $Q$.
(2) (i)$\Rightarrow$(ii) : dériver en $\theta=0$. (ii)$\Rightarrow$(iii) : $X$ commute avec $P_q=f(Q)$
(théorème 1.5(d)), donc $P_qXP_{q'}=XP_qP_{q'}=0$ pour $q\neq q'$. (iii)$\Rightarrow$(i) :
$XU_\theta=\sum_qe^{iq\theta}P_qXP_q=U_\theta X$. Les opérateurs diagonaux par blocs forment
$\bigoplus_qM_{d_q}\otimes1_1$ (dans une base orthonormée adaptée aux images des $P_q$), dont le commutant est
$\bigoplus_q1_{d_q}\otimes M_1=\operatorname{span}\{P_q\}$ : un opérateur qui commute avec eux commute avec chaque $P_q$,
donc il est diagonal par blocs, et chaque bloc commute avec $M_{d_q}$, donc est un scalaire (lemme 2.5, étape 5, avec
$m=1$). (3) $\frac d{dt}e^{itH}Qe^{-itH}=ie^{itH}[H,Q]e^{-itH}$ (proposition 8.2) s'annule pour tout $t$ si
et seulement si $[H,Q]=0$. Alors $e^{itH}$ commute avec $P_q$ (théorème 1.5(d), deux fois), et
$\operatorname{Tr}(P_q\rho(t))=\operatorname{Tr}(e^{itH}P_qe^{-itH}\rho)=\operatorname{Tr}(P_q\rho)$. $\square$

Les représentations irréductibles de $U(1)$ sont donc les caractères $\theta\mapsto e^{iq\theta}$, $q\in\mathbb Z$ (de
dimension un, comme pour tout groupe abélien), et la multiplicité de la charge $q$ est $d_q$. Les charges
s'additionnent dans les produits tensoriels : si $u_k$ et $w_l$ sont des bases propres orthonormées de $Q_1$ et $Q_2$, de
valeurs propres $\alpha_k$ et $\beta_l$, la base orthonormée $u_k\otimes w_l$ est formée de vecteurs propres de
$Q_1\otimes1+1\otimes Q_2$ (valeur propre $\alpha_k+\beta_l$) et de $e^{i\theta Q_1}\otimes e^{i\theta Q_2}$ (valeur propre
$e^{i\theta\alpha_k}e^{i\theta\beta_l}$), donc, par définition du calcul fonctionnel,
$$
e^{i\theta(Q_1\otimes1+1\otimes Q_2)}=e^{i\theta Q_1}\otimes e^{i\theta Q_2},
$$
et cet opérateur agit sur l'image de $P_{q_1}\otimes P_{q_2}$ avec la charge $q_1+q_2$. C'est la règle de fusion
$(q_1,q_2)\mapsto q_1+q_2$ des produits tensoriels fusionnés du chapitre 7.

Pour $n$ qubits, la charge est l'**opérateur nombre** $\hat N=\sum_jn_j$, où $n_j=|1\rangle\langle1|$ agit sur le
qubit $j$ : $\hat N|b\rangle=|b|\,|b\rangle$, avec le poids de Hamming $|b|$. Ses secteurs sont de dimension
$\binom nw$, de sorte que les opérateurs qui conservent le poids de Hamming forment
$$
U(1)'\cong\bigoplus_{w=0}^nM_{\binom nw}(\mathbb C),
$$
qui est `chains.SpinChain(n).algebra`, le produit fusionné de $n$ copies de $M_1\oplus M_1$ avec les charges $0,1$
(bit $1$ = spin haut, $\hat N=S^z+\frac n2$ ; dans `decompose.total_spin`, le spin haut est au contraire le *premier*
vecteur de base de $\mathbb C^2$, ce qui renomme la base mais pas les secteurs). Tout hamiltonien XXZ
$\sum_{\langle ij\rangle}J(s^x_is^x_j+s^y_is^y_j+\Delta s^z_is^z_j)$ conserve $\hat N$ : $s^z_is^z_j$ est diagonal dans la base
$|b\rangle$, et le terme de saut $s^x_is^x_j+s^y_is^y_j=\frac12(s^+_is^-_j+s^-_is^+_j)$ (développer $s^\pm=s^x\pm is^y$ ;
les termes croisés se compensent) envoie $|b\rangle$ sur $0$ ou sur la configuration dans laquelle un spin haut est passé
de l'un des sites $i,j$ à l'autre, donc il envoie chaque sous-espace propre de $\hat N$ dans lui-même.

```python
import functools
from torch_vn_algebra import SpinChain

L = 4
Sx, Sy, Sz = dc.total_spin(L)                              # (16, 16) each, site 0 most significant
N_op = Sz + L / 2 * torch.eye(2 ** L, dtype=dt)            # number of up spins
DU = dc.from_generators(N_op)                              # the algebra generated by U(1)
print(DU)
chain = SpinChain(L, boundary='periodic')                  # U(1)' = (+)_w M_binom(L, w)
assert set(DU.k_factors) == {1} and sorted(DU.multiplicities) == sorted(chain.algebra.k_factors)
assert sorted(chain.algebra.k_factors) == sorted(math.comb(L, w) for w in range(L + 1))

sx, sy, sz = dc.spin_matrices(0.5)
def two_site(a, b, i, j, n=L):                             # a on site i, b on site j of n sites
    ops = [a if k == i else b if k == j else torch.eye(2, dtype=dt) for k in range(n)]
    return functools.reduce(torch.kron, ops)

H = sum(two_site(sx, sx, i, j) + two_site(sy, sy, i, j) + 0.5 * two_site(sz, sz, i, j)
        for i, j in chain.bonds)                           # XXZ ring, Delta = 0.5
assert close(H @ N_op, N_op @ H)                           # [H, N] = 0 ...
assert close(DU.from_commutant(DU.to_commutant(H)), H)     # ... so H is block diagonal (Proposition 4.10)
Hc = chain.xxz(J=1.0, Delta=0.5)                           # the same H, built sector by sector
E_sec = torch.cat([torch.linalg.eigvalsh(Hc.matrix[0, N, :k, :k]) for N, k in enumerate(chain.algebra.k_factors)])
assert close(torch.sort(E_sec).values, torch.linalg.eigvalsh(H))
```

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). La symétrie $U(1)$ est la conservation d'un nombre de particules, d'une charge ou
> d'une aimantation. La proposition 4.10(3) est le théorème de Noether dans ce cadre : la charge est conservée si et
> seulement si le hamiltonien est symétrique, et alors la dynamique ne quitte jamais un secteur (théorème 8.3). Si
> de plus rien ne peut mesurer les cohérences entre secteurs, les secteurs sont des secteurs de supersélection
> (chapitre 3).

## 4.5 $SU(2)$ : le spin

**Définition.** Une **représentation de spin** sur $\mathbb C^d$ est un triplet de matrices hermitiennes
$S_x,S_y,S_z$ telles que
$$
[S_x,S_y]=iS_z,\qquad[S_y,S_z]=iS_x,\qquad[S_z,S_x]=iS_y,
$$
c'est-à-dire une représentation des étiquettes $\{x,y,z\}$ (section 4.1) vérifiant ces relations. L'exemple de base
est $s_a=\frac12\sigma_a$ sur $\mathbb C^2$ (matrices de Pauli). Si $(S_a)$ et $(T_a)$ sont des représentations de
spin, $S_a\otimes1+1\otimes T_a$ en est une aussi (les termes croisés commutent) ; par conséquent le **spin total**
$S_a=\sum_is_a^{(i)}$ de $L$ spins $\frac12$ (`decompose.total_spin(L)`) en est une. Les rotations
$e^{-i\theta\,\mathbf n\cdot\mathbf S}$ ($\theta\in\mathbb R$, $\mathbf n\in\mathbb R^3$ un vecteur unitaire, calcul fonctionnel)
commutent avec $X$ pour tous $\theta,\mathbf n$ si et seulement si $X$ commute avec les trois $S_a$ (dériver en
$\theta=0$, proposition 8.2(1) ; réciproquement, théorème 1.5(d)), de sorte que le commutant de « l'action de $SU(2)$ »
est $\{S_x,S_y,S_z\}'$, et d'après le théorème 2.2 les rotations et les $S_a$ engendrent la même algèbre (les deux
ensembles sont stables par adjonction). Pour $L$ spins $\frac12$, les rotations sont les $V^{\otimes L}$, $V\in SU(2)$.
En effet, $\mathbf n\cdot\mathbf S=\sum_ih^{(i)}$ avec $h=\mathbf n\cdot\mathbf s$ et, comme pour les charges à la
section 4.4 (récurrence sur $L$), $e^{-i\theta\,\mathbf n\cdot\mathbf S}=(e^{-i\theta h})^{\otimes L}$, avec
$e^{-i\theta h}\in SU(2)$ (ses valeurs propres sont $e^{\mp i\theta/2}$, car $h$ est de trace nulle avec $h^2=\frac14$).
Réciproquement, un $V\in SU(2)$ est normal, de valeurs propres de module $1$ et de produit $1$, donc
$V=W\operatorname{diag}(e^{i\varphi},e^{-i\varphi})W^*$ avec $W$ unitaire, et $V=e^{iA}$ pour la matrice hermitienne de trace
nulle $A=W\operatorname{diag}(\varphi,-\varphi)W^*$. Une matrice hermitienne $2\times2$ de trace nulle est une combinaison
réelle $A=2\,\mathbf c\cdot\mathbf s$ des matrices de Pauli, donc $V=e^{-i\theta\,\mathbf n\cdot\mathbf s}$ avec
$\theta=2|\mathbf c|$, $\mathbf n=-\mathbf c/|\mathbf c|$ (et $\theta=0$ si $\mathbf c=0$). Le théorème B.6(a),(c) et le
corollaire B.7 (avec $d=2$) montrent directement que les $V^{\otimes L}$, $V\in SU(2)$, engendrent linéairement
l'algèbre engendrée par les $S_a$.

Posons $S_\pm=S_x\pm iS_y$, de sorte que $S_-=S_+^*$, et $\mathbf S^2=S_x^2+S_y^2+S_z^2$. Les relations donnent
$$
[S_z,S_\pm]=\pm S_\pm,\qquad[S_+,S_-]=2S_z,\qquad\mathbf S^2=S_-S_++S_z^2+S_z,
$$
et $\mathbf S^2$ commute avec chaque $S_a$ (par exemple
$[\mathbf S^2,S_z]=S_x[S_x,S_z]+[S_x,S_z]S_x+S_y[S_y,S_z]+[S_y,S_z]S_y=-i\{S_x,S_y\}+i\{S_y,S_x\}=0$).

**Théorème 4.11 (représentations de spin irréductibles).** Pour tout $s\in\{0,\frac12,1,\frac32,\dots\}$, les
matrices de `decompose.spin_matrices(s)` — sur $\mathbb C^{2s+1}$ de base $|m\rangle$, $m=s,s-1,\dots,-s$,
$$
S_z|m\rangle=m|m\rangle,\qquad S_+|m\rangle=\sqrt{s(s+1)-m(m+1)}\;|m+1\rangle,\qquad S_-=S_+^*,
$$
forment une représentation de spin irréductible (le **spin $s$**) avec $\mathbf S^2=s(s+1)1$. Toute représentation
de spin irréductible est unitairement équivalente à exactement l'une d'elles.

*Démonstration.* *Le spin $s$ est une représentation de spin.* $[S_z,S_+]=S_+$ est clair.
$S_+S_-|m\rangle=(s(s+1)-(m-1)m)|m\rangle$ et $S_-S_+|m\rangle=(s(s+1)-m(m+1))|m\rangle$, donc $[S_+,S_-]=2S_z$ ; ces
deux relations et leurs adjointes équivalent aux trois relations de définition.
$\mathbf S^2|m\rangle=(s(s+1)-m(m+1)+m^2+m)|m\rangle=s(s+1)|m\rangle$. *Irréductibilité :* un $W\neq0$ invariant est
invariant par $S_z$, dont les valeurs propres sont simples, donc il contient un certain $|m\rangle$ ; les coefficients
de $S_\pm$ ne s'annulent qu'en $m=\pm s$, donc en appliquant $S_\pm$ on obtient tous les $|m'\rangle$.

*Classification.* Soit $(S_a)$ irréductible sur $\mathbb C^d$, $\lambda$ la plus grande valeur propre de $S_z$ et
$v_0$ un vecteur propre unitaire. $S_zS_+v_0=S_+(S_z+1)v_0=(\lambda+1)S_+v_0$, donc $S_+v_0=0$ par maximalité.
Posons $v_j=S_-^jv_0$ ; alors $S_zv_j=(\lambda-j)v_j$ et, par récurrence, en utilisant $S_+S_-=S_-S_++2S_z$,
$$
S_+v_j=j(2\lambda-j+1)\,v_{j-1}
$$
($S_+v_j=S_-S_+v_{j-1}+2S_zv_{j-1}=[(j-1)(2\lambda-j+2)+2(\lambda-j+1)]v_{j-1}$). Les $v_j$ non nuls sont des vecteurs
propres de $S_z$ pour des valeurs propres distinctes, il y a donc un dernier vecteur non nul, $v_r$. Le sous-espace
engendré par $v_0,\dots,v_r$ est invariant par $S_z,S_\pm$, donc égal à $\mathbb C^d$, et $d=r+1$. De
$0=S_+v_{r+1}=(r+1)(2\lambda-r)v_r$ on tire $\lambda=r/2=:s$. Enfin
$\|v_j\|^2=\langle v_{j-1},S_+v_j\rangle=j(2s-j+1)\|v_{j-1}\|^2$, de sorte que les vecteurs unitaires
$|s-j\rangle=v_j/\|v_j\|$ vérifient $S_+|m\rangle=\sqrt{(s-m)(s+m+1)}\,|m+1\rangle$ avec $m=s-j$, ce qui est la formule
ci-dessus. Des valeurs différentes de $s$ donnent des dimensions différentes. $\square$

**Corollaire 4.12.** Toute représentation de spin sur $\mathbb C^n$ est unitairement équivalente à
$\bigoplus_s(\text{spin }s)\otimes1_{m_s}$, et, avec $d_M=\dim\ker(S_z-M)$,
$$
m_s=d_s-d_{s+1}.
$$

*Démonstration.* Théorème 4.3 et théorème 4.11. Chaque copie du spin $s'$ contient exactement un vecteur propre de
$S_z$ pour la valeur propre $M$ si $|M|\le s'$ et $s'-M\in\mathbb Z$, et aucun sinon ; donc
$d_M=\sum_{s'\ge|M|,\,s'-M\in\mathbb Z}m_{s'}$, et la différence de deux valeurs consécutives est $m_s$. $\square$

**Proposition 4.13 ($L$ spins $\frac12$).** Le spin total sur $(\mathbb C^2)^{\otimes L}$ contient le spin $S$ si et
seulement si $S\in\{\frac L2,\frac L2-1,\dots\}$, $S\ge0$, avec la multiplicité
$$
m_S=\binom L{\frac L2-S}-\binom L{\frac L2-S-1}\qquad\Big(\tbinom L{-1}=0\Big),
$$
de sorte que l'algèbre engendrée par $S_x,S_y,S_z$ est $\bigoplus_SM_{2S+1}\otimes1_{m_S}$ et son commutant
$\bigoplus_S1_{2S+1}\otimes M_{m_S}$.

*Démonstration.* $S_z=\sum_is_z^{(i)}$ est diagonal dans la base produit, de valeur propre $u-\frac L2$ sur un produit
de $u$ spins haut et $L-u$ spins bas, donc $d_M=\binom L{L/2+M}=\binom L{L/2-M}$ si
$\frac L2+M\in\{0,\dots,L\}$ et $0$ sinon. On applique le corollaire 4.12 et le théorème 4.3. Les valeurs sont
strictement positives car $\binom L{j+1}/\binom Lj=\frac{L-j}{j+1}>1$ pour $j<\frac{L-1}2$. $\square$

Pour $L=2$ : $\frac12\otimes\frac12=1\oplus0$ (**Clebsch–Gordan**) : le triplet $S_-^j|{\uparrow\uparrow}\rangle$,
$j=0,1,2$, symétrique sous l'échange $F$ des deux spins, et le singulet
$\frac1{\sqrt2}(|{\uparrow\downarrow}\rangle-|{\downarrow\uparrow}\rangle)$, antisymétrique et annulé par $S_+$ et $S_z$.
Pour $L=4$ : $m_2=1$, $m_1=4-1=3$, $m_0=6-4=2$, donc
$$
\{S_x,S_y,S_z\}''\cong M_5\oplus M_3\oplus M_1\ \text{ with multiplicities }1,3,2,\qquad
\{S_x,S_y,S_z\}'\cong M_1\oplus M_3\oplus M_2,
$$
de dimensions $35$ et $14$ — les nombres trouvés au chapitre 10. La bibliothèque retrouve la même chose à partir
des trois seuls générateurs (secteurs triés selon $(k_c,m_c)$) :

```python
S4 = dc.total_spin(4)
assert close(S4[0] @ S4[1] - S4[1] @ S4[0], 1j * S4[2])                # [S_x, S_y] = i S_z
D4 = dc.from_generators(S4)
print(D4)
assert D4.k_factors == (1, 3, 5) and D4.multiplicities == (2, 3, 1)     # M_1, M_3, M_5 with m = 2, 3, 1
assert D4.verify()['ok'] and D4.dims == (35, 14)
m_S = {S: math.comb(4, 2 - S) - (math.comb(4, 1 - S) if S < 2 else 0) for S in (0, 1, 2)}
assert [m_S[(k - 1) // 2] for k in D4.k_factors] == list(D4.multiplicities)    # Proposition 4.13
casimir = D4.to_algebra((S4 @ S4).sum(0))                                  # S^2 = S(S+1) on spin S
for c, k in enumerate(D4.k_factors):
    assert close(casimir.matrix[0, c, :k, :k], (k - 1) / 2 * (k + 1) / 2 * torch.eye(k))
```

Quels opérateurs commutent avec toutes les rotations ? Les permutations des sites, en tout cas : la conjugaison par
$\pi(\sigma)$ envoie $s^{(i)}_a$ sur $s^{(\sigma(i))}_a$ (lemme B.1(b)), donc elle fixe $\sum_is^{(i)}_a$. Réciproquement,
rien d'autre ne le fait :

**Théorème 4.14 (dualité de Schur–Weyl).** Soit $\pi$ la représentation de $S_L$ qui permute les facteurs de
$(\mathbb C^2)^{\otimes L}$. (a) $\operatorname{span}\pi(S_L)\subseteq\{S_x,S_y,S_z\}'$. (b) L'égalité a lieu ;
de manière équivalente (théorème 2.2), l'algèbre engendrée par le spin total est $\pi(S_L)'$.

*Démonstration.* La partie (a) vient d'être démontrée. La partie (b) est le corollaire B.7, le cas $d=2$ du
théorème B.6, qui démontre la dualité pour $(\mathbb C^d)^{\otimes L}$ et l'action diagonale de $U(d)$ pour tous $d$ et
$L$. $\square$

(Pour aller plus loin : W. Fulton et J. Harris, *Representation Theory*, Springer 1991, leçon 6 ; R. Goodman et
N. Wallach, *Symmetry, Representations, and Invariants*, Springer 2009, section 9.1.) Pour $L\le5$, le code compare
numériquement les deux espaces : d'après (a), il suffit qu'ils aient la même dimension.

```python
for n in range(2, 6):
    comm = dc.commutant_basis(dc.total_spin(n)).reshape(-1, 4 ** n)          # a basis of {S_a}'
    perm = dc.permutation_representation(n, 2).reshape(-1, 4 ** n)          # all n! permutations
    dims = [torch.linalg.matrix_rank(A).item() for A in (comm, perm, torch.cat([comm, perm]))]
    print(n, dims)
    assert dims[0] == dims[1] == dims[2]                                     # span pi(S_n) = {S_a}'
Dp = dc.from_group(dc.permutation_representation(4, 2))                   # pi(S_4)'' = {S_a}'
print(Dp)
assert sorted(zip(Dp.multiplicities, Dp.k_factors)) == sorted(zip(D4.k_factors, D4.multiplicities))
```

Les dimensions affichées sont $2,5,14,42$, les nombres de Catalan $\binom{2L}L/(L+1)$ ; le corollaire B.10 démontre
$\dim\{S_a\}'=\sum_Sm_S^2=\binom{2L}L/(L+1)$ pour tout $L$. La dernière ligne montre la
dualité à l'œuvre : $\pi(S_4)''$ a des blocs $k=(1,2,3)$ de multiplicités $(5,1,3)$ — les couples $(k_c,m_c)$ du
spin total, rôles échangés. En général (corollaire B.9), l'espace de multiplicité du spin $S$ porte une
représentation irréductible de $S_L$, de dimension $m_S$, et ces représentations sont deux à deux non équivalentes ;
ici les dimensions sont $1,3,2$ pour $S=2,1,0$.

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Un multiplet de spin $S$ est formé de $2S+1$ états de même énergie pour tout
> hamiltonien invariant par rotation, distingués par $S^z=M$. Deux spins $\frac12$ forment un triplet et un
> singulet ; quatre spins forment un quintuplet, trois triplets et deux singulets. Les espaces de multiplicité des
> deux singulets et des trois triplets sont le qubit sans décohérence et le qutrit sans bruit du chapitre 10.

## 4.6 Systèmes de spins

Le **modèle de Heisenberg** sur $L$ spins $\frac12$ est $H=\sum_{i<j}J_{ij}\,\mathbf s^{(i)}\cdot\mathbf s^{(j)}$ avec
des couplages réels $J_{ij}$.

**Proposition 4.15.** (a) (*Identité de Dirac*) $\mathbf s^{(i)}\cdot\mathbf s^{(j)}=\frac12F_{ij}-\frac14$, où $F_{ij}$
échange les spins $i$ et $j$. Par conséquent tout hamiltonien de Heisenberg appartient à
$\operatorname{span}\pi(S_L)\subseteq\{S_a\}'$, et $V^*HV=\bigoplus_S1_{2S+1}\otimes H_S$ avec des blocs $H_S\in M_{m_S}$
dans la décomposition de la proposition 4.13 : chaque valeur propre de $H_S$ est un niveau $(2S+1)$ fois dégénéré
de $H$.
(b) L'espace des **vecteurs de plus haut poids** $\ker S_+\cap\ker(S_z-S)$ est de dimension $m_S$, invariant par
$H$, et $H$ y agit comme $H_S$ (à un changement de base unitaire près).

*Démonstration.* (a) Sur deux spins,
$\mathbf s^{(1)}\cdot\mathbf s^{(2)}=\frac12(\mathbf S^2-\mathbf s^{(1)2}-\mathbf s^{(2)2})=\frac12\mathbf S^2-\frac34$
(théorème 4.11 avec $s=\frac12$). Les vecteurs du triplet $S_-^j|{\uparrow\uparrow}\rangle$ sont $|{\uparrow\uparrow}\rangle$,
$|{\downarrow\uparrow}\rangle+|{\uparrow\downarrow}\rangle$, $2|{\downarrow\downarrow}\rangle$, et avec le singulet ils forment une base
orthogonale de $\mathbb C^2\otimes\mathbb C^2$. Sur le triplet, $\mathbf S^2=2$ et $F=1$ : $\mathbf S^2=S_-S_++S_z^2+S_z$ vaut
$0+1+1=2$ sur $|{\uparrow\uparrow}\rangle$, et $\mathbf S^2$ et $F$ commutent avec $S_-$, tandis que $F$ fixe
$|{\uparrow\uparrow}\rangle$. Sur le singulet, $\mathbf S^2=0$ (il est annihilé par $S_+$ et $S_z$) et $F=-1$. Les deux
membres coïncident sur ces quatre vecteurs. Pour des sites $i,j$ parmi $L$, on tensorise avec l'identité. Le reste
découle du théorème 4.14(a) et du théorème 4.3. (b) Nous démontrons davantage : pour toute représentation de spin
décomposée comme dans le corollaire 4.12 et tout $\mu\in\frac12\mathbb Z$, $\ker S_+\cap\ker(S_z-\mu)$ est (l'image par
l'unitaire du corollaire 4.12 de) $|\mu,\mu\rangle\otimes\mathbb C^{m_\mu}$ si $\mu\ge0$, et $0$ si $\mu<0$. En effet, dans
cette base, un vecteur avec $S_z=\mu$ s'écrit $\sum_{S'}|S',\mu\rangle\otimes w_{S'}$, la somme portant sur les $S'\ge|\mu|$
avec $S'-\mu\in\mathbb Z$, et d'après le théorème 4.11, $S_+$ l'envoie sur
$\sum_{S'}\sqrt{(S'-\mu)(S'+\mu+1)}\,|S',\mu+1\rangle\otimes w_{S'}$, où $S'+\mu+1\ge1$ ; il s'annule donc si et seulement si
$w_{S'}=0$ pour tout $S'\neq\mu$. Pour le $H$ de Heisenberg, $1\otimes H_S$ agit sur $|S,S\rangle\otimes\mathbb C^{m_S}$
comme $H_S$. $\square$

La partie (b) est la manière dont le module `su2` (chaînes de spins, chapitre 9) construit $H_S$ sans décomposition
numérique : il construit une base orthonormée des vecteurs de plus haut poids, en ajoutant un site après l'autre,
avec les coefficients de Clebsch–Gordan de $j\otimes\frac12$. Nous les obtenons en descendant à partir du plus haut
poids. D'après le théorème 4.11 et $S_-=S_+^*$, le spin $j$ vérifie
$$
S_\pm|j,m\rangle=\sqrt{(j\mp m)(j\pm m+1)}\;|j,m\pm1\rangle ,
$$
et nous posons $|j,m\rangle=0$ pour $|m|>j$ (ce qui est cohérent avec la formule, dont le coefficient s'annule en
$m=\pm j$).

**Théorème 4.16 (coefficients de Clebsch–Gordan de $j\otimes\frac12$).** Soient $j\in\{\frac12,1,\frac32,\dots\}$,
$S_a=S^{(j)}_a\otimes1+1\otimes s_a$ sur $\mathbb C^{2j+1}\otimes\mathbb C^2$ (le spin $j$ du théorème 4.11, de base
$|j,m\rangle$, couplé à un spin $\frac12$ avec $|{\uparrow}\rangle=|\frac12,\frac12\rangle$, $|{\downarrow}\rangle=|\frac12,-\frac12\rangle$),
et posons $\alpha_M=\sqrt{(j+M+\frac12)/(2j+1)}$, $\beta_M=\sqrt{(j-M+\frac12)/(2j+1)}$. Pour $J=j\pm\frac12$ et
$M=J,J-1,\dots,-J$, soit
$$
\begin{aligned}
|j+\tfrac12,M\rangle\!\rangle&=\alpha_M\,|j,M-\tfrac12\rangle\otimes|{\uparrow}\rangle+\beta_M\,|j,M+\tfrac12\rangle\otimes|{\downarrow}\rangle,\\
|j-\tfrac12,M\rangle\!\rangle&=-\beta_M\,|j,M-\tfrac12\rangle\otimes|{\uparrow}\rangle+\alpha_M\,|j,M+\tfrac12\rangle\otimes|{\downarrow}\rangle .
\end{aligned}
$$
(a) Ces $2(2j+1)$ vecteurs forment une base orthonormée de $\mathbb C^{2j+1}\otimes\mathbb C^2$.
(b) $|j+\frac12,j+\frac12\rangle\!\rangle=|j,j\rangle\otimes|{\uparrow}\rangle$, et $|j-\frac12,j-\frac12\rangle\!\rangle$ est, à un facteur de
module un près, l'unique vecteur unitaire $x$ tel que $S_zx=(j-\frac12)x$ et $S_+x=0$ ; le facteur est fixé en imposant
un coefficient positif devant $|j,j\rangle\otimes|{\downarrow}\rangle$ (la convention de Condon–Shortley).
(c) $S_z|J,M\rangle\!\rangle=M|J,M\rangle\!\rangle$ et $S_-|J,M\rangle\!\rangle=\sqrt{(J+M)(J-M+1)}\,|J,M-1\rangle\!\rangle$ : les états
$|J,M\rangle\!\rangle$ s'obtiennent à partir de $|J,J\rangle\!\rangle$ par descente et normalisation.
(d) Par conséquent $S_a|J,M\rangle\!\rangle=\sum_{M'}(S^{(J)}_a)_{M'M}|J,M'\rangle\!\rangle$ pour $a=x,y,z$ : l'unitaire
$|J,M\rangle\mapsto|J,M\rangle\!\rangle$ envoie spin $(j+\frac12)\oplus$ spin $(j-\frac12)$ sur $S_a$, donc
$j\otimes\frac12\cong(j+\frac12)\oplus(j-\frac12)$, et les coefficients de Clebsch–Gordan
$\langle j,m;\frac12,\sigma|J,M\rangle$ (les coordonnées de $|J,M\rangle\!\rangle$ dans la base produit) sont $\pm\alpha_M,\beta_M$
comme ci-dessus.

*Démonstration.* (a) Si un vecteur de base figurant dans la formule est nul, son coefficient l'est aussi :
$|j,M-\frac12\rangle=0$ seulement pour $M=-j-\frac12$, où $\alpha_M=0$, et $|j,M+\frac12\rangle=0$ seulement pour
$M=j+\frac12$, où $\beta_M=0$. Chaque vecteur a donc pour norme $\alpha_M^2+\beta_M^2=1$. Des vecteurs de $M$ différents
font intervenir des vecteurs différents de la base produit, et pour un même $M$ le produit scalaire vaut
$-\alpha_M\beta_M+\beta_M\alpha_M=0$. Il y en a $(2j+2)+2j=2(2j+1)$.
(b) La première assertion est $\alpha_{j+1/2}=1$, $\beta_{j+1/2}=0$. Le sous-espace propre $S_z=j-\frac12$ est engendré
par $|j,j-1\rangle\otimes|{\uparrow}\rangle$ et $|j,j\rangle\otimes|{\downarrow}\rangle$ et, avec $s_+|{\uparrow}\rangle=0$,
$s_+|{\downarrow}\rangle=|{\uparrow}\rangle$,
$$
S_+\big(p\,|j,j-1\rangle\otimes|{\uparrow}\rangle+q\,|j,j\rangle\otimes|{\downarrow}\rangle\big)=(p\sqrt{2j}+q)\,|j,j\rangle\otimes|{\uparrow}\rangle ,
$$
qui s'annule si et seulement si $q=-\sqrt{2j}\,p$ ; avec $|p|^2+|q|^2=1$ et $q>0$ (donc $p<0$), cela donne
$p=-\beta_{j-1/2}=-1/\sqrt{2j+1}$, $q=\alpha_{j-1/2}=\sqrt{2j/(2j+1)}$. (c) $S_z$ : les deux termes ont
$S_z=(M\mp\frac12)\pm\frac12=M$. Pour $S_-$, soit
$x_M=p\,|j,M-\frac12\rangle\otimes|{\uparrow}\rangle+q\,|j,M+\frac12\rangle\otimes|{\downarrow}\rangle$ ; comme
$s_-|{\uparrow}\rangle=|{\downarrow}\rangle$ et $s_-|{\downarrow}\rangle=0$,
$$
S_-x_M=p\sqrt{(j+M-\tfrac12)(j-M+\tfrac32)}\;|j,M-\tfrac32\rangle\otimes|{\uparrow}\rangle
+\Big(p+q\sqrt{(j+M+\tfrac12)(j-M+\tfrac12)}\Big)|j,M-\tfrac12\rangle\otimes|{\downarrow}\rangle .
$$
Nous utilisons $\beta_M\sqrt{(j+M+\frac12)(j-M+\frac12)}=(j-M+\frac12)\alpha_M$ et
$\alpha_M\sqrt{(j+M+\frac12)(j-M+\frac12)}=(j+M+\frac12)\beta_M$.
Pour $J=j+\frac12$, $(p,q)=(\alpha_M,\beta_M)$ et $c=\sqrt{(J+M)(J-M+1)}=\sqrt{(j+M+\frac12)(j-M+\frac32)}$ : le premier
coefficient vaut $\sqrt{(j+M+\frac12)(j+M-\frac12)(j-M+\frac32)/(2j+1)}=c\,\alpha_{M-1}$, le second
$\alpha_M(j-M+\frac32)=c\,\beta_{M-1}$.
Pour $J=j-\frac12$, $(p,q)=(-\beta_M,\alpha_M)$ et $c=\sqrt{(j+M-\frac12)(j-M+\frac12)}$ : le premier coefficient vaut
$-\sqrt{(j-M+\frac12)(j+M-\frac12)(j-M+\frac32)/(2j+1)}=-c\,\beta_{M-1}$, le second $\beta_M(j+M-\frac12)=c\,\alpha_{M-1}$.
Dans les deux cas $S_-|J,M\rangle\!\rangle=c\,|J,M-1\rangle\!\rangle$ pour $M>-J$. Pour $M=-J$, les deux coefficients
s'annulent ($p=\alpha_M=0$ si $J=j+\frac12$ ; $j+M-\frac12=0$ si $J=j-\frac12$), donc $S_-|J,-J\rangle\!\rangle=0$. (d) D'après
(a) et (c), dans la base orthonormée $|J,M\rangle\!\rangle$, les matrices de $S_z$ et $S_-$ sont celles de spin
$(j+\frac12)\oplus$ spin $(j-\frac12)$ (théorème 4.11 et la formule ci-dessus) ; il en va de même de la matrice de
$S_+=S_-^*$, donc de celles de $S_x=\frac12(S_++S_-)$ et $S_y=\frac1{2i}(S_+-S_-)$. $\square$

L'exercice 2 étend la décomposition (mais pas les coefficients) à $j_1\otimes j_2$. Pour la bibliothèque, seuls
comptent les deux vecteurs de plus haut poids du théorème 4.16, et ils fonctionnent dans toute représentation de spin,
irréductible ou non :

**Corollaire 4.17 (vecteurs de plus haut poids site par site).** Soient $(T_a)$ une représentation de spin sur $E$,
$\mathcal H_\mu(E)=\ker T_+\cap\ker(T_z-\mu)$, muni d'une base orthonormée $B_\mu$ ($B_\mu=\emptyset$ si $\mathcal H_\mu(E)=0$),
et $S_a=T_a\otimes1+1\otimes s_a$ sur $E\otimes\mathbb C^2$. Pour tout $S\in\{0,\frac12,1,\dots\}$, les vecteurs
$$
v\otimes|{\uparrow}\rangle\quad(v\in B_{S-1/2}),\qquad
\sqrt{\tfrac{2S+1}{2S+2}}\;v\otimes|{\downarrow}\rangle-\tfrac1{\sqrt{(2S+1)(2S+2)}}\;T_-v\otimes|{\uparrow}\rangle\quad(v\in B_{S+1/2})
$$
forment une base orthonormée de $\ker S_+\cap\ker(S_z-S)$ (la première famille est vide pour $S=0$). En partant d'un
spin ($E=\mathbb C^2$, $B_{1/2}=\{|{\uparrow}\rangle\}$, tous les autres $B_\mu$ vides) et en ajoutant $L-1$ spins, on obtient
des bases orthonormées des espaces de vecteurs de plus haut poids de $L$ spins $\frac12$, de dimensions $m_S$
(proposition 4.15(b)).

*Démonstration.* Notons $a=\sqrt{(2S+1)/(2S+2)}$, $b=-1/\sqrt{(2S+1)(2S+2)}$,
$u(v)=a\,v\otimes|{\downarrow}\rangle+b\,T_-v\otimes|{\uparrow}\rangle$. Pour $v\in\mathcal H_{S+1/2}(E)$,
$T_+T_-v=(T_-T_++2T_z)v=(2S+1)v$, et $a+(2S+1)b=0$. *Plus haut poids :* $S_z$ vaut $S$ sur les deux familles
($T_z=S\mp\frac12$ sur $v$ et $S-\frac12$ sur $T_-v$) ; $S_+(v\otimes|{\uparrow}\rangle)=T_+v\otimes|{\uparrow}\rangle=0$, et
$S_+u(v)=a\,T_+v\otimes|{\downarrow}\rangle+a\,v\otimes|{\uparrow}\rangle+b\,T_+T_-v\otimes|{\uparrow}\rangle=(a+(2S+1)b)\,v\otimes|{\uparrow}\rangle=0$.
*Orthonormalité :* $\langle v\otimes|{\uparrow}\rangle,u(v')\rangle=b\langle v,T_-v'\rangle=b\langle T_+v,v'\rangle=0$ pour
$v\in B_{S-1/2}$, et
$\langle u(v),u(v')\rangle=a^2\langle v,v'\rangle+b^2\langle v,T_+T_-v'\rangle=(a^2+(2S+1)b^2)\langle v,v'\rangle=\langle v,v'\rangle$.
*Génération :* soit $S_zw=Sw$ et $S_+w=0$, et écrivons $w=x\otimes|{\uparrow}\rangle+y\otimes|{\downarrow}\rangle$. Alors
$T_zx=(S-\frac12)x$, $T_zy=(S+\frac12)y$ et $0=S_+w=(T_+x+y)\otimes|{\uparrow}\rangle+T_+y\otimes|{\downarrow}\rangle$, donc $T_+y=0$,
c'est-à-dire $y\in\mathcal H_{S+1/2}(E)$, et $T_+x=-y$. Alors $x'=x+\frac1{2S+1}T_-y$ vérifie $T_+x'=-y+y=0$, donc
$x'\in\mathcal H_{S-1/2}(E)$, qui est nul pour $S=0$ (démonstration de la proposition 4.15(b)). Comme $b/a=-1/(2S+1)$,
$w=x'\otimes|{\uparrow}\rangle+a^{-1}u(y)$, avec $u$ prolongée linéairement, appartient au sous-espace engendré par les deux
familles. Pour $L$ spins, on applique ceci avec $E=(\mathbb C^2)^{\otimes(n-1)}$ et son spin total $T_a$, $n=2,\dots,L$ ;
sur un spin, $\mathcal H_{1/2}=\mathbb C|{\uparrow}\rangle$ et tous les autres $\mathcal H_\mu=0$ (le seul spin est $\frac12$).
$\square$

La seconde famille est formée des vecteurs $|S,S\rangle\!\rangle$ du théorème 4.16 avec $j=S+\frac12$, écrits dans la copie du
spin $j$ engendrée par les $T_-^kv$ (démonstration du théorème 4.11, $|j,j-1\rangle=T_-v/\sqrt{2j}$), et la première
famille est $|j+\frac12,j+\frac12\rangle\!\rangle$ avec $j=S-\frac12$. Étiqueter chaque vecteur de base par sa suite de spins
intermédiaires $S_1=\frac12,S_2,\dots,S_L=S$, $|S_n-S_{n-1}|=\frac12$, donne la base de Yamanouchi que construit
`su2.SU2Chain` (son code source utilise exactement ces coefficients ; le code ci-dessous ne vérifie que les conclusions
du corollaire 4.17), à l'intérieur du secteur $U(1)$ $\hat N=\frac L2+S$ de `SpinChain` (dont l'ordre des sites est un
réétiquetage des facteurs tensoriels, qui commute avec le spin total d'après le lemme B.1(b)). Le code vérifie le
théorème 4.16 sur des matrices denses et par comparaison avec `su2.clebsch_gordan` de la bibliothèque, puis les
conclusions du corollaire 4.17 pour les bases de huit spins. Nous vérifions aussi la proposition 4.15(a) sur un
anneau de quatre spins avec `D4`, et comparons les tailles des blocs pour huit spins.

```python
from torch_vn_algebra import su2

def ket(j, m):                                            # |j, m> of Theorem 4.11 (zero if |m| > j)
    v = torch.zeros(int(2 * j + 1), dtype=dt)
    if abs(m) <= j:
        v[int(j - m)] = 1
    return v

up, dn = ket(0.5, 0.5), ket(0.5, -0.5)
for j in (0.5, 1.0, 1.5, 2.0):                            # Theorem 4.16
    Sj = dc.spin_matrices(j)
    sp_j = [math.sqrt(j * (j + 1) - (j - k) * (j - k + 1)) for k in range(1, int(2 * j) + 1)]
    assert close(Sj[0] + 1j * Sj[1], torch.diag(torch.tensor(sp_j, dtype=dt), 1))     # S_+ of Theorem 4.11
    assert close(Sj[2], torch.diag(torch.tensor([j - k for k in range(int(2 * j) + 1)], dtype=dt)))
    Sc = torch.stack([torch.kron(a, torch.eye(2, dtype=dt)) + torch.kron(torch.eye(len(a), dtype=dt), b)
                      for a, b in zip(Sj, dc.spin_matrices(0.5))])
    Sp, Sm = Sc[0] + 1j * Sc[1], Sc[0] - 1j * Sc[1]
    cols = []
    for J in (j + 0.5, j - 0.5):
        for M in [J - k for k in range(int(2 * J + 1))]:
            al, be = math.sqrt((j + M + 0.5) / (2 * j + 1)), math.sqrt((j - M + 0.5) / (2 * j + 1))
            p, q = (al, be) if J > j else (-be, al)
            x = p * torch.kron(ket(j, M - 0.5), up) + q * torch.kron(ket(j, M + 0.5), dn)
            assert close(Sc[2] @ x, M * x)
            if M == J:
                assert close(Sp @ x, 0)                               # highest weight
            else:
                assert close(Sm @ cols[-1], math.sqrt((J + M + 1) * (J - M)) * x)   # lowering, (c)
            for s2, e in ((0.5, up), (-0.5, dn)):                     # the library's Racah formula agrees
                cg = su2.clebsch_gordan(j, M - s2, 0.5, s2, J, M) if abs(M - s2) <= j else 0.0
                assert abs(cg - torch.vdot(torch.kron(ket(j, M - s2), e), x).real.item()) < 1e-12
            cols.append(x)
    U = torch.stack(cols, 1)
    assert U.shape == (2 * len(Sj[0]), 2 * len(Sj[0])) and close(U.mH @ U, torch.eye(U.shape[1]))   # (a)

def swap(i, j, n):
    p = list(range(n))
    p[i], p[j] = j, i
    return dc.permutation_representation(n, 2, [p])[0]

ring = [(i, (i + 1) % 4) for i in range(4)]
H4 = sum(swap(i, j, 4) / 2 - torch.eye(16, dtype=dt) / 4 for i, j in ring)        # Dirac's identity
assert close(H4, sum(two_site(a, a, i, j) for a in (sx, sy, sz) for i, j in ring))
h = D4.to_commutant(H4)                                   # blocks H_S of sizes m_S = 2, 3, 1
assert close(D4.from_commutant(h), H4)                    # H4 lies in {S_a}'
levels = torch.cat([torch.linalg.eigvalsh(h.matrix[0, c, :m, :m]).repeat_interleave(k)
                    for c, (k, m) in enumerate(zip(D4.k_factors, D4.multiplicities))])
assert close(torch.sort(levels).values, torch.linalg.eigvalsh(H4))   # level of H_S repeated 2S+1 times

u1, sym = SpinChain(8, boundary='periodic'), su2.SU2Chain(8, boundary='periodic')
for S, m in zip(sym.spins, sym.multiplicities):           # Corollary 4.17: orthonormal highest-weight bases
    B, N = sym.hw_basis(S), int(4 + S)                     # columns in the U(1) sector N = L/2 + S
    assert B.shape == (math.comb(8, N), m) and close(B.T @ B, torch.eye(m, dtype=B.dtype))
    assert N == 8 or close(sym.raise_(B, N), 0)            # S_+ B = 0 (S_z = S holds in sector N)
H_u1, H_su2 = u1.xxz(J=1.0, Delta=1.0), sym.heisenberg(J=1.0)
E_u1 = torch.sort(torch.cat([torch.linalg.eigvalsh(H_u1.matrix[0, N, :k, :k])
                             for N, k in enumerate(u1.algebra.k_factors)])).values
E_su2, S_lev, _ = sym.full_spectrum(H_su2, return_quantum_numbers=True)   # level of block S, 2S+1 times
E_su2, S_lev = E_su2[0], S_lev[0]
assert close(E_u1, E_su2)
assert S_lev[0] == 0 and (S_lev[1:4] == 1).all() and E_su2[4] - E_su2[3] > 1e-3     # singlet, then triplet
assert E_su2[1] - E_su2[0] > 1e-3 and close(E_su2[1:4], E_su2[1].expand(3))       # simple singlet (Theorem 9.9)
assert close(sym.restrict(H_u1).matrix, H_su2.matrix)      # Proposition 4.15(b)
flops = lambda ks: sum(k ** 3 for k in ks)
print(u1.algebra.k_factors, sym.multiplicities, flops(u1.algebra.k_factors), flops(sym.multiplicities))
print(f"ground state energy {E_su2[0].item():.6f}, first excited {E_su2[1].item():.6f} (x3)")
```

Pour l'anneau de huit spins, les blocs de $U(1)$ ont pour tailles $\binom8w=1,8,28,56,70,56,28,8,1$ (section 4.4), les
blocs de $SU(2)$ les tailles $m_S=14,28,20,7,1$ pour $S=0,\dots,4$ (proposition 4.13). Le coût de la diagonalisation,
$\sum k^3$, passe de $739162$ à $33040$, soit un facteur $22$ ; la dimension de l'algèbre, $\sum k^2$, passe de
$\sum_w\binom8w\binom8{8-w}=\binom{16}8=12870$ (choisir $8$ éléments parmi $16$, dont $w$ parmi les $8$ premiers) à
$\sum_Sm_S^2=1430=\binom{16}8/9$ (corollaire B.10). D'après le théorème 9.9 (l'anneau est de longueur paire et $J>0$),
l'état fondamental est un singulet unique et $E_0(0)<E_0(1)$. *Observation numérique* (calculée par le code
ci-dessus) : l'énergie fondamentale vaut $-3.651093$, et le premier niveau excité est un triplet à $-3.128419$.

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Le modèle de Heisenberg décrit
> l'interaction d'échange entre spins électroniques localisés dans les isolants magnétiques ; $J_{ij}>0$ favorise les
> voisins antiparallèles. Son spectre est formé de multiplets, et un champ magnétique $h\,S^z$ scinde chaque multiplet
> en $2S+1$ niveaux équidistants (exercice 4). Le théorème 9.9 (couplages bipartis ; pour les chaînes et anneaux
> antiferromagnétiques de longueur paire, un état fondamental singulet) remonte à E. Lieb et D. Mattis
> (*J. Math. Phys.* 3 (1962) 749).

## 4.7 Points de comparaison

Pour les groupes finis, GAP et SageMath calculent exactement les tables de caractères, et le corollaire 4.5 donne
alors les multiplicités ; RepLAB (MATLAB ; D. Rosset, F. Montealegre-Mora, J.-D. Bancal) décompose numériquement les
représentations des groupes finis et compacts, lui aussi à l'aide d'éléments aléatoires du commutant, et
SymbolicWedderburn.jl (Julia) calcule des décompositions de Wedderburn pour la réduction par symétrie de programmes
semi-définis. En physique, QuSpin et d'autres codes de diagonalisation exacte diagonalisent par blocs les symétries
$U(1)$ et de réseau, ITensor et TeNPy gèrent des charges abéliennes dans les réseaux de tenseurs, et TensorKit.jl
traite $SU(2)$ avec des arbres de fusion. Le module `decompose` de la bibliothèque est une méthode dense pour de
petites matrices ($n$ jusqu'à quelques centaines) ; il accepte n'importe quelle liste de matrices (éléments de groupe
ou générateurs), certifie son résultat (proposition 4.9) et renvoie l'algèbre et son commutant sous forme d'objets
`TypeIAlgebra` munis de `to_algebra`, `from_algebra`, `to_commutant`, `from_commutant`, de sorte que les chapitres 3
et 5–8 s'appliquent immédiatement. Pour les chaînes de spins $\frac12$, `su2.SU2Chain` fournit directement les blocs
de $SU(2)$ (corollaire 4.17) ; il accepte jusqu'à $L=18$ spins, une limite fixée par la mémoire des bases.

## Exercices

1. Montrer que, pour la représentation de $S_L$ sur $(\mathbb C^d)^{\otimes L}$,
   $\chi_\pi(\sigma)=d^{\#\text{cycles}(\sigma)}$ et $\dim\pi(S_L)'=\binom{d^2+L-1}L$. Vérifier $20$ (section 2.4) et
   $35$ ($=\dim\{S_a\}''$ pour quatre spins).
2. (*Clebsch–Gordan.*) Montrer que $j_1\otimes j_2\cong\bigoplus_jj$, $j=|j_1-j_2|,|j_1-j_2|+1,\dots,j_1+j_2$, pour
   la représentation de spin $S^{(j_1)}_a\otimes1+1\otimes S^{(j_2)}_a$ — la règle de fusion non abélienne de
   $SU(2)$ que mentionne le chapitre 7 — et vérifier $1\otimes\frac12$ et $1\otimes1$ avec `from_generators`.
3. Soit $\theta/2\pi$ irrationnel. Montrer que $\{U_\theta\}'=U(1)'$ pour un seul $U_\theta=e^{i\theta Q}$, mais que
   $\{e^{i\pi\hat N}\}'\cong M_{2^{n-1}}\oplus M_{2^{n-1}}$ pour $n$ qubits.
4. (*Effet Zeeman.*) Soit $H$ un hamiltonien de Heisenberg sur $L\ge1$ spins et $H_h=H+hS_z$. Montrer que $H_h$
   commute avec $S_z$ et $\mathbf S^2$ mais, pour $h\neq0$, pas avec $S_x$, et que ses valeurs propres sont
   $E+hM$, $M=-S,\dots,S$, pour les valeurs propres $E$ de $H_S$. Le vérifier pour l'anneau de huit spins.

<details>
<summary>Solutions</summary>

1. D'après le lemme B.1(c), $\pi(\sigma)$ est une matrice de permutation qui fixe $|b_1\dots b_L\rangle$ si et seulement
   si $b$ est constant sur les cycles de $\sigma$, et il y a $d^{\#\text{cycles}(\sigma)}$ tels $b$. D'après le
   corollaire 4.5(c),
   $\dim\pi(S_L)'=\frac1{L!}\sum_\sigma(d^2)^{\#\text{cycles}(\sigma)}=\frac{x(x+1)\cdots(x+L-1)}{L!}$ avec $x=d^2$
   (l'identité $\sum_\sigma x^{\#\text{cycles}(\sigma)}=x(x+1)\cdots(x+L-1)$ est démontrée dans la solution du chapitre 5,
   exercice 2), ce qui vaut $\binom{d^2+L-1}L$ (la proposition B.3 démontre la même dimension en dénombrant des
   multiensembles) : $\binom63=20$ et $\binom74=35$. Que $35$ soit aussi $\dim\{S_a\}''$, c'est le théorème 4.14(b)
   pour $L=4$ : $\{S_a\}''=\pi(S_4)'$.

```python
for n, value in ((3, 20), (4, 35)):
    assert dc.commutant_dimension(dc.permutation_representation(n, 2)) == math.comb(n + 3, n) == value
```

2. $S_a=S^{(j_1)}_a\otimes1+1\otimes S^{(j_2)}_a$ est une représentation de spin (section 4.5), et l'unitaire
   $x\otimes y\mapsto y\otimes x$ l'entrelace avec celle de $(j_2,j_1)$ ; on peut donc supposer $j_1\ge j_2$. $S_z$ a la
   valeur propre $m_1+m_2$ sur $|m_1\rangle\otimes|m_2\rangle$, donc ses valeurs propres sont dans $j_1+j_2+\mathbb Z$, et pour
   un tel $M$, $d_M$ est le nombre de $m_2\in\{-j_2,\dots,j_2\}$ tels que $|M-m_2|\le j_1$ (alors $m_1=M-m_2\in j_1+\mathbb Z$
   automatiquement) ; $d_M=0$ pour $M\notin j_1+j_2+\mathbb Z$. Pour $M\ge0$, la condition s'écrit
   $M-j_1\le m_2\le M+j_1$, et $M+j_1\ge j_2$, donc $m_2$ va de $\max(-j_2,M-j_1)$ à $j_2$ : $d_M=2j_2+1$ pour
   $0\le M\le j_1-j_2$, $d_M=j_1+j_2-M+1$ pour $j_1-j_2\le M\le j_1+j_2$, et $d_M=0$ pour $M>j_1+j_2$. D'après le
   corollaire 4.12, la multiplicité du spin $j\ge0$ est $m_j=d_j-d_{j+1}$. Elle vaut $0$ si $j\notin j_1+j_2+\mathbb Z$ ou
   $j>j_1+j_2$ ; $m_{j_1+j_2}=1-0=1$ ; pour $j_1-j_2\le j<j_1+j_2$ (avec $j\in j_1+j_2+\mathbb Z$, donc $j+1\le j_1+j_2$),
   $m_j=(j_1+j_2-j+1)-(j_1+j_2-j)=1$ ; et pour $0\le j<j_1-j_2$ (donc $j+1\le j_1-j_2$), $m_j=(2j_2+1)-(2j_2+1)=0$.
   Par conséquent $j_1\otimes j_2\cong\bigoplus_{j=j_1-j_2}^{j_1+j_2}j$, chaque spin une fois. Pour $j_2=\frac12$, c'est le
   théorème 4.16(d), qui donne en outre les coefficients.

```python
def coupled(a, b):                                        # S^(a) (x) 1 + 1 (x) S^(b)
    Sa, Sb = dc.spin_matrices(a), dc.spin_matrices(b)
    ea, eb = torch.eye(Sa.shape[-1], dtype=dt), torch.eye(Sb.shape[-1], dtype=dt)
    return torch.stack([torch.kron(x, eb) + torch.kron(ea, y) for x, y in zip(Sa, Sb)])
for a, b, ks in ((1.0, 0.5, (2, 4)), (1.0, 1.0, (1, 3, 5))):  # 1 (x) 1/2 = 1/2 (+) 3/2, 1 (x) 1 = 0 (+) 1 (+) 2
    D_cg = dc.from_generators(coupled(a, b))
    assert D_cg.k_factors == ks and set(D_cg.multiplicities) == {1}
```

3. $U_\theta=\sum_qe^{iq\theta}P_q$, et $e^{iq\theta}=e^{iq'\theta}$ si et seulement si $(q-q')\theta\in2\pi\mathbb Z$, ce
   qui, pour $\theta/2\pi$ irrationnel, impose $q=q'$. Les sous-espaces propres de $U_\theta$ sont donc les images des
   $P_q$, et un opérateur commute avec la matrice normale $U_\theta$ si et seulement s'il préserve ses sous-espaces
   propres (théorème 1.5(d) et interpolation de Lagrange comme à la section 2.2), c'est-à-dire si et seulement si
   $X=\sum_qP_qXP_q$ ; c'est $U(1)'$ d'après la proposition 4.10. L'argument utilise seulement que les $e^{iq\theta}$
   pour les charges $q$ qui apparaissent sont distincts ; le code prend $\theta=1$ et les charges $0,\dots,4$ de quatre
   qubits, pour lesquelles $0<|q-q'|\le4<2\pi$ suffit.
   Pour $\theta=\pi$, la valeur propre $(-1)^w$ ne voit que la parité du poids de Hamming : deux sous-espaces propres de
   dimension $2^{n-1}$ (pour $n\ge1$, les configurations de poids pair et de poids impair ; inverser le premier bit est
   une bijection entre elles), et le commutant est $M_{2^{n-1}}\oplus M_{2^{n-1}}$.

```python
for theta, m in ((1.0, (1, 1, 4, 4, 6)), (math.pi, (8, 8))):   # theta = 1: e^{iq}, q = 0..4, distinct as 4 < 2 pi
    D_theta = dc.from_generators(torch.linalg.matrix_exp(1j * theta * N_op))
    assert D_theta.multiplicities == m and set(D_theta.k_factors) == {1}
```

4. $S_z$ et $\mathbf S^2$ commutent avec $H$ (proposition 4.15(a)) et avec $S_z$, donc avec $H_h$ ; et
   $[H_h,S_x]=h[S_z,S_x]=ihS_y\neq0$. Dans la décomposition de la proposition 4.13,
   $V^*S_zV=\bigoplus_SS^{(S)}_z\otimes1_{m_S}$ et $V^*HV=\bigoplus_S1_{2S+1}\otimes H_S$ (proposition 4.15), donc sur le
   terme de spin $S$, $V^*H_hV=hS^{(S)}_z\otimes1+1\otimes H_S$, une somme d'opérateurs qui commutent, diagonalisée
   par les produits de vecteurs propres, de valeurs propres $hM+E$.

```python
assert not close(H4 @ Sx - Sx @ H4 + 0.3 * (Sz @ Sx - Sx @ Sz), 0)          # [H4 + 0.3 S_z, S_x] != 0
E_field = sym.full_spectrum(sym.heisenberg(J=1.0, h=0.3))[0]                # E + h M on every multiplet
H_field = u1.xxz(J=1.0, Delta=1.0, h=0.3)
E_u1_field = torch.cat([torch.linalg.eigvalsh(H_field.matrix[0, N, :k, :k])
                        for N, k in enumerate(u1.algebra.k_factors)])
assert close(torch.sort(E_u1_field).values, E_field)
```

</details>

**Suite :** le [chapitre 5](05_random_matrices.md) tire des matrices aléatoires dans une algèbre : des unitaires de
Haar de $\prod_cU(k_c)$, le calcul de Weingarten — dont le cas d'ordre deux repose sur le commutant de
$\{V\otimes V\}$, le plus petit exemple de dualité de Schur–Weyl — et les statistiques de niveaux secteur par
secteur.
