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
représentations de spin irréductibles et la décomposition de $L$ spins $\frac12$, pourquoi les hamiltoniens
invariants par $SU(2)$ ont des blocs bien plus petits que les hamiltoniens invariants par $U(1)$, et comment la
conjugaison complexe répartit les secteurs d'une algèbre engendrée par des matrices réelles en type réel, complexe et
quaternionique. Les vérifications numériques de ce chapitre se trouvent dans le
[Tutoriel 4](../../tutorial/fr/04_symmetries.md).

Dans tout ce chapitre, $S'$ est le commutant de $S\subseteq M_n(\mathbb C)$ (chapitre 2). Si $S^*=S$, $S''$ est
l'**algèbre de von Neumann engendrée par** $S$ : c'est une $*$-algèbre contenant $S$, et toute $*$-algèbre
$\mathcal B\supseteq S$ vérifie $\mathcal B=\mathcal B''\supseteq S''$ (théorème 2.2). $\mathcal A'_h=\{K\in\mathcal A':K=K^*\}$
est un espace vectoriel réel. $\mathbb C^k\otimes\mathbb C^m$ est muni de l'ordre de Kronecker du chapitre 7
($e_a\otimes f_j$ est le vecteur de base $am+j$, en numérotant à partir de $0$), de sorte que $1_k\otimes x$ est
diagonale par blocs avec $k$ copies de $x$. $|G|$ désigne l'ordre d'un groupe fini, et $\hbar=1$.

<a id="s-4-1"></a>

## 4.1 Représentations et lemme de Schur

**Définition.** Soit $X$ un ensemble muni d'une involution $x\mapsto x^*$ ($x^{**}=x$). Une
**représentation** de $X$ sur $\mathbb C^k$ est une application $\rho:X\to M_k(\mathbb C)$ telle que
$\rho(x^*)=\rho(x)^*$. Un sous-espace $K$ est **invariant** si $\rho(x)K\subseteq K$ pour tout $x$ ; $\rho$ est
**irréductible** si $k\ge1$ et si $0,\mathbb C^k$ sont les seuls sous-espaces invariants. Un
**opérateur d'entrelacement** de $\rho$ vers une représentation $\rho'$ sur $\mathbb C^{k'}$ est une matrice $T$
de taille $k'\times k$ telle que $T\rho(x)=\rho'(x)T$ pour tout $x$ ; ils forment un espace vectoriel
$\mathrm{Hom}(\rho,\rho')$, et $\mathrm{Hom}(\rho,\rho)=\rho(X)'$. $\rho$ et $\rho'$ sont **équivalentes**,
$\rho\cong\rho'$, s'il existe un opérateur d'entrelacement inversible.

Trois cas se présentent : un groupe $G$ avec $g^*=g^{-1}$, où l'on exige que $\rho$ soit un morphisme
$G\to U(k)$ (une **représentation unitaire** ; alors $\rho(g^{-1})=\rho(g)^{-1}=\rho(g)^*$) ; une $*$-algèbre munie
de son adjoint, où $\rho$ est un $*$-morphisme unifère ; et des étiquettes de générateurs avec $x^*=x$, par exemple
$\{x,y,z\}$, représentées par des familles de matrices hermitiennes telles que les matrices de spin
(section 4.5). Si $K$ est invariant et si les colonnes de $Q\in M_{k\times d}$ forment une base orthonormée de $K$,
alors $\rho(x)Q=QQ^*\rho(x)Q$, de sorte que $\rho_K(x)=Q^*\rho(x)Q$ est encore une représentation (un morphisme si
$\rho$ en est un), la **restriction** de $\rho$ à $K$.

<a id="r-4-1"></a>

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

<a id="r-4-2"></a>

**Lemme 4.2 (lemme de Schur).** Soient $\rho,\rho'$ des représentations irréductibles de $X$ sur
$\mathbb C^k,\mathbb C^{k'}$ et $T\in\mathrm{Hom}(\rho,\rho')$.
(a) $T=0$ ou $T$ est inversible. (b) $\mathrm{Hom}(\rho,\rho)=\mathbb C1$.
(c) $T^*\in\mathrm{Hom}(\rho',\rho)$, et si $T\neq0$, alors $T=\lambda U$ avec $\lambda>0$ et $U$ unitaire.
Par conséquent $\dim\mathrm{Hom}(\rho,\rho')$ vaut $1$ si $\rho\cong\rho'$ et $0$ sinon, et des représentations
irréductibles équivalentes sont unitairement équivalentes.

*Démonstration.* (a) et (b) ont été démontrés pour les groupes après la proposition 2.7, et la démonstration
n'utilise que l'invariance de $\ker T$ par $\rho$ ($T\rho(x)v=\rho'(x)Tv=0$ si $Tv=0$) et celle de l'image de $T$
par $\rho'$ : ainsi $\ker T\in\{0,\mathbb C^k\}$ et $\mathrm{ran}\,T\in\{0,\mathbb C^{k'}\}$ ; et pour (b),
$T-\lambda1$, où $\lambda$ est une valeur propre de $T$, est un opérateur d'entrelacement non inversible, donc
nul. (c) L'adjoint de $T\rho(x^*)=\rho'(x^*)T$ est $\rho(x)T^*=T^*\rho'(x)$. Alors $T^*T\in\mathrm{Hom}(\rho,\rho)$,
donc $T^*T=\lambda^21$ avec $\lambda>0$ si $T\neq0$ ($T^*T\ge0$ est non nul), et $U=T/\lambda$ est carrée d'après
(a), avec $U^*U=1$. Enfin, si $T_1,T_2\in\mathrm{Hom}(\rho,\rho')$ sont non nuls, alors
$T_1^{-1}T_2\in\mathrm{Hom}(\rho,\rho)=\mathbb C1$. $\square$

Si tous les $\rho(x)$ commutent entre eux (un groupe abélien), une $\rho$ irréductible est de dimension un :
chaque $\rho(x)$ appartient à $\mathrm{Hom}(\rho,\rho)=\mathbb C1$, donc tout sous-espace est invariant.

<a id="s-4-2"></a>

## 4.2 La décomposition isotypique et le commutant

<a id="r-4-3"></a>

**Théorème 4.3 (décomposition isotypique).** Soit $\rho$ une représentation de $X$ sur $\mathbb C^n$.
Il existe des représentations irréductibles deux à deux non équivalentes $\rho_1,\dots,\rho_C$ de $X$ sur
$\mathbb C^{k_c}$, des entiers $m_c\ge1$ avec $\sum_ck_cm_c=n$, et un opérateur unitaire
$V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n$ tels que
$V^*\rho(x)V=\bigoplus_c\rho_c(x)\otimes1_{m_c}$ pour tout $x$. Pour toute décomposition de ce type,
$$
V^*\rho(X)'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),\qquad
V^*\rho(X)''V=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},
$$
et $m_c=\dim\mathrm{Hom}(\rho_c,\rho)$ ; une représentation irréductible $\sigma$ de $X$ vérifie
$\mathrm{Hom}(\sigma,\rho)\neq0$ si et seulement si $\sigma\cong\rho_c$ pour exactement un $c$.

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
$\mathrm{Hom}(\rho_d,\rho_c)$, qui vaut $0$ pour $c\neq d$ et $\mathbb C1$ pour $c=d$ (lemme 4.2). Ainsi
$Y\in\rho(X)'$ si et seulement si $Y_{cj,dl}=\delta_{cd}(y_c)_{jl}1_{k_c}$, c'est-à-dire
$V^*YV=\bigoplus_c1_{k_c}\otimes y_c$.

*Bicommutant.* Si $Z$ commute avec $\rho(X)'$, il commute avec les projecteurs sur les termes de la somme
(prendre $y_c=1$ et $y_d=0$ pour $d\neq c$), donc $V^*ZV=\bigoplus_cZ_c$, et $Z_c$ commute avec
$1_{k_c}\otimes M_{m_c}$ ; d'après le lemme 2.5, étape 5 (en échangeant les deux facteurs tensoriels),
$Z_c\in M_{k_c}\otimes1_{m_c}$. Réciproquement, ces opérateurs commutent avec $\rho(X)'$.

*Multiplicités.* Pour $T\in\mathrm{Hom}(\sigma,\rho)$ avec $\sigma$ irréductible, le même argument montre que
$Q_{c,j}^*T\in\mathrm{Hom}(\sigma,\rho_c)$, qui vaut $\mathbb Ct_{c}$ avec $t_c$ inversible si $\sigma\cong\rho_c$ et
$0$ sinon ; et $T=\sum_{c,j}Q_{c,j}Q_{c,j}^*T$. Comme les $\rho_c$ sont deux à deux non équivalentes, $\sigma$ est
équivalente à au plus l'une d'elles, disons $\rho_c$ (sinon $T=0$), et alors $T=\sum_ja_jQ_{c,j}t_c$ ;
réciproquement, chaque $Q_{c,j}t_c$ appartient à $\mathrm{Hom}(\sigma,\rho)$, et ces $m_c$ applications ont des
images orthogonales, donc $\dim\mathrm{Hom}(\sigma,\rho)=m_c$. $\square$

*Dans la bibliothèque :* une décomposition comme dans le théorème 4.3 est calculée par `decompose.from_group` et `decompose.from_generators` — [tutoriel, section 4.3](../../tutorial/fr/04_symmetries.md#s-4-3).

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

<a id="r-4-4"></a>

**Théorème 4.4 (orthogonalité de Schur).** Soient $\rho,\rho'$ des représentations unitaires irréductibles d'un
groupe fini $G$ sur $\mathbb C^k,\mathbb C^{k'}$. Alors
$$
\frac1{|G|}\sum_{g\in G}\rho(g)_{ij}\,\overline{\rho'(g)_{i'j'}}=
\begin{cases}\delta_{ii'}\delta_{jj'}/k&\text{if }\rho=\rho',\\0&\text{if }\rho\not\cong\rho'.\end{cases}
$$

*Démonstration.* Pour une matrice $Y$ de taille $k\times k'$, posons $T=\frac1{|G|}\sum_g\rho(g)Y\rho'(g)^*$. Pour
$h\in G$, la substitution $g\to hg$ et l'égalité $\rho'(hg)^*=\rho'(g)^*\rho'(h)^{-1}$ donnent $\rho(h)T=T\rho'(h)$,
donc $T\in\mathrm{Hom}(\rho',\rho)$. Si $\rho\not\cong\rho'$, $T=0$ ; si $\rho=\rho'$, $T=\lambda1$ avec
$k\lambda=\mathrm{Tr}\,T=\mathrm{Tr}\,Y$ (lemme 4.2). Pour $Y=|j\rangle\langle j'|$, le coefficient $T_{ii'}$ est
le membre de gauche, et $\mathrm{Tr}\,Y=\delta_{jj'}$. $\square$

Le **caractère** d'une représentation est $\chi_\rho(g)=\mathrm{Tr}\,\rho(g)$ ; des représentations équivalentes
ont le même caractère, puisque $\mathrm{Tr}(T\rho(g)T^{-1})=\mathrm{Tr}\,\rho(g)$. On note
$\langle\varphi,\psi\rangle_G=\frac1{|G|}\sum_g\overline{\varphi(g)}\psi(g)$.

<a id="r-4-5"></a>

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
$\rho\not\cong\rho'$. (b) La trace est unitairement invariante et $\mathrm{Tr}(a\otimes1_m)=m\mathrm{Tr}\,a$, donc
$\chi_\pi=\sum_cm_c\chi_{\rho_c}$ ; on applique (a). (c) Par le même développement et (a),
$\langle\chi_\pi,\chi_\pi\rangle_G=\sum_cm_c^2$, qui est $\dim\pi(G)'$ d'après le théorème 4.3. $\square$

*Calcul :* le corollaire 4.5 pour $S_3$ permutant trois qubits — [tutoriel, section 4.2](../../tutorial/fr/04_symmetries.md#s-4-2) ; le théorème 4.4 pour la représentation standard de $S_3$ — [tutoriel, section 4.3](../../tutorial/fr/04_symmetries.md#s-4-3).

Une seconde démonstration de (c) correspond à la manière dont la bibliothèque échantillonne les commutants. Le
**twirl** $P(Y)=\frac1{|G|}\sum_g\pi(g)Y\pi(g)^*$ est idempotent, d'image $\pi(G)'$ (chapitre 2, exercice 3), et
autoadjoint pour $\langle Y,Z\rangle_{HS}=\mathrm{Tr}(Y^*Z)$ (substituer $g\to g^{-1}$) ; c'est donc le projecteur
orthogonal sur $\pi(G)'$, et $\dim\pi(G)'$ est sa trace. L'application $Y\mapsto AYB$ a pour trace
$\mathrm{Tr}\,A\mathrm{Tr}\,B$ (ses coefficients diagonaux sont les $A_{ii}B_{jj}$), ce qui donne
$\frac1{|G|}\sum_g|\chi_\pi(g)|^2$.

*Dans la bibliothèque :* le twirl est `decompose.twirl` — [tutoriel, section 4.3](../../tutorial/fr/04_symmetries.md#s-4-3) ; $\dim\pi(G)'$ est `decompose.commutant_dimension` — [tutoriel, section 4.2](../../tutorial/fr/04_symmetries.md#s-4-2).

Reprenons l'exemple de la section 2.4 à l'aide des caractères. Une matrice de permutation a pour trace le nombre
de vecteurs de base qu'elle fixe, et $\pi(\sigma)$ fixe $|b_1b_2b_3\rangle$ si et seulement si $b$ est constant sur
les cycles de $\sigma$ (lemme B.1(c)), donc $\chi_\pi(\sigma)=2^{\#\text{cycles}(\sigma)}$. Trois représentations unitaires
de $S_3$ sont à portée de main. Soit $P_\sigma e_i=e_{\sigma(i)}$ la représentation par permutations sur $\mathbb C^3$ (un
homomorphisme à valeurs dans les matrices de permutation). La représentation **triviale** est $\sigma\mapsto1$. La
**signature** est $\mathrm{sgn}\,\sigma=\det P_\sigma$, un homomorphisme à valeurs dans $\{\pm1\}\subseteq U(1)$ car le
déterminant est multiplicatif et le déterminant d'une matrice de permutation vaut $\pm1$. La droite $\mathbb C(1,1,1)$
est invariante, donc $E=(1,1,1)^\perp$ l'est aussi (lemme 4.1(a)) ; la représentation **standard** est la restriction
de $P$ à $E$. Dans une base orthonormée adaptée à $\mathbb C(1,1,1)\oplus E$, $P_\sigma=1\oplus\mathrm{std}(\sigma)$,
donc $\chi_{\rm std}(\sigma)=\mathrm{Tr}\,P_\sigma-1$, le nombre de points fixes moins un. Les six éléments de $S_3$
sont l'identité, trois transpositions et deux $3$-cycles, avec $3,2,1$ cycles et $3,1,0$ points fixes ; une
transposition échange deux colonnes de $1_3$, donc sa signature vaut $-1$, et un $3$-cycle est un produit de deux
transpositions ($(1\,2\,3)=(1\,3)(1\,2)$), donc sa signature vaut $+1$. On obtient la table

| | identité ($1$) | transpositions ($3$) | $3$-cycles ($2$) |
|---|---|---|---|
| $\chi_{\rm triv}$ | $1$ | $1$ | $1$ |
| $\mathrm{sgn}$ | $1$ | $-1$ | $1$ |
| $\chi_{\rm std}$ | $2$ | $0$ | $-1$ |
| $\chi_\pi$ | $8$ | $4$ | $2$ |

(ainsi $\mathrm{sgn}\,\sigma=(-1)^{3-\#\text{cycles}(\sigma)}$ ; tous les caractères sont réels). Avec
$\langle\varphi,\psi\rangle_G$ sommé classe par classe :
$\langle\chi_{\rm triv},\chi_{\rm triv}\rangle_G=\frac{1+3+2}6=1$, $\langle\mathrm{sgn},\mathrm{sgn}\,\rangle_G=\frac{1+3+2}6=1$ et
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
(corollaire 4.5(c)).

C'est la décomposition de la section 2.4, obtenue cette fois à partir d'une table de six nombres par représentation.

*Calcul :* cette table, recalculée à partir de $\pi$ — [tutoriel, section 4.2](../../tutorial/fr/04_symmetries.md#s-4-2).

<a id="s-4-3"></a>

## 4.3 Des générateurs à l'algèbre : un algorithme randomisé

Le théorème 4.3 décrit la décomposition mais ne la calcule pas : sa démonstration a besoin des sous-espaces
invariants irréductibles. Murota, Kanno, Kojima et Kojima (*Japan J. Indust. Appl. Math.* 27 (2010) 125–160 ;
voir aussi Maehara et Murota, *ibid.* 263–293) ont observé que les sous-espaces propres d'un *élément aléatoire du
commutant* les trouvent. Soit $\mathcal A$ l'algèbre de von Neumann engendrée par des matrices données
$S_1,\dots,S_J$ et leurs adjointes, avec une décomposition $V^*\mathcal AV=\bigoplus_cM_{k_c}\otimes1_{m_c}$,
$V^*\mathcal A'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$ (théorème 2.6) ; seule l'analyse utilise $V$.

<a id="r-4-6"></a>

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

<a id="r-4-7"></a>

**Proposition 4.7 (généricité).** Soient $K$ et $Y$ des éléments aléatoires indépendants de $\mathcal A'_h$ et de
$\mathcal A'$ dont les lois ont des densités par rapport à la mesure de Lebesgue sur ces espaces vectoriels réels
(lemme A.3). Presque sûrement : (a) l'hypothèse du lemme 4.6 est satisfaite ; (b) $\Pi_iY\Pi_j\neq0$ pour deux
projecteurs spectraux distincts quelconques $\Pi_i,\Pi_j$ de $K$ dont les images sont dans le même terme de la somme
directe.

*Démonstration.* (a) La matrice hermitienne $\kappa=\bigoplus_c\kappa_c\in M_N$ dépend de $K$ de façon $\mathbb R$-linéaire.
Le déterminant de Gram $\delta(K)=\det[\mathrm{Tr}\,\kappa^{i+j}]_{i,j=0}^{N-1}$ de $1,\kappa,\dots,\kappa^{N-1}$ pour
$\langle\cdot,\cdot\rangle_{HS}$ est un polynôme réel en les coordonnées de $K$, non nul si et seulement si ces
puissances sont linéairement indépendantes. Comme $p(\kappa)=U\mathrm{diag}(p(\lambda_l))U^*$ pour
$\kappa=U\mathrm{diag}(\lambda)U^*$, un polynôme non nul de degré $<N$ annule $\kappa$ si et seulement s'il
s'annule sur $\mathrm{spec}\,\kappa$, ce qui est possible si et seulement si $\kappa$ a moins de $N$ valeurs
propres distinctes (prendre $\prod_\mu(t-\mu)$ sur celles-ci). Ainsi $\delta(K)\neq0$ est l'hypothèse du lemme 4.6,
et $\delta\not\equiv0$ : il est non nul en $V(\bigoplus_c1\otimes\mathrm{diag}(\text{réels distincts}))V^*\in\mathcal A'_h$.
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

<a id="r-4-8"></a>

**Lemme 4.8 (compressions ; alignement des copies).** (a) Soit $W$ un sous-espace $\mathcal A$-invariant de base
orthonormée $Q\in M_{n\times d}$ et $\sigma(a)=Q^*aQ$ la restriction. Alors $\sigma(\mathcal A)$ est une $*$-algèbre,
$\{Q^*YQ:Y\in\mathcal A'\}=\sigma(\mathcal A)'$ et $\{Q^*KQ:K\in\mathcal A'_h\}=\sigma(\mathcal A)'_h$. $\mathcal A$ agit
irréductiblement sur $W$ si et seulement si $Q^*YQ\in\mathbb C1$ pour tout $Y\in\mathcal A'$. Sinon, si $X$ est un élément
aléatoire de $\mathcal A'$ ou de $\mathcal A'_h$ à densité, alors $Q^*XQ\notin\mathbb C1$ presque sûrement.

(b) Soient $W_i,W_j$ des sous-espaces invariants sur lesquels $\mathcal A$ agit
irréductiblement, de bases orthonormées $Q_i,Q_j$ et de restrictions $\sigma_i(a)=Q_i^*aQ_i$, et soit
$Y\in\mathcal A'$. Alors $T=Q_i^*YQ_j\in\mathrm{Hom}(\sigma_j,\sigma_i)$. Si $T\neq0$, alors $T=\lambda U$ avec
$\lambda>0$, $U$ unitaire, et $Q_j'=Q_jU^*$ vérifie $Q_j'^*aQ_j'=\sigma_i(a)$ pour tout $a\in\mathcal A$. Si
$W_i,W_j$ sont dans des termes différents ($W_i\subseteq P_c\mathbb C^n$, $W_j\subseteq P_d\mathbb C^n$, $c\neq d$),
alors $T=0$.

*Démonstration.* (a) Comme dans la démonstration du théorème 4.3 ($W$ est invariant par $a$ et par $a^*$),
$aQ=Q\sigma(a)$ et $Q^*a=\sigma(a)Q^*$ pour $a\in\mathcal A$. Donc $\sigma(ab)=Q^*aQ\sigma(b)=\sigma(a)\sigma(b)$,
$\sigma(a^*)=\sigma(a)^*$ et $\sigma(1)=Q^*Q=1_d$, si bien que $\sigma$ est un $*$-homomorphisme unifère et $\sigma(\mathcal A)$
une $*$-algèbre ; et pour $Y\in\mathcal A'$, $Q^*YQ\sigma(a)=Q^*YaQ=Q^*aYQ=\sigma(a)Q^*YQ$. Réciproquement, pour
$T\in\sigma(\mathcal A)'$ posons $Y=QTQ^*$ ; alors $aY=Q\sigma(a)TQ^*=QT\sigma(a)Q^*=Ya$, donc $Y\in\mathcal A'$, avec
$Q^*YQ=T$, et $Y$ est hermitien si $T$ l'est. Si $\sigma$ est irréductible,
$\sigma(\mathcal A)'=\mathrm{Hom}(\sigma,\sigma)=\mathbb C1$ (lemme 4.2(b)). Si $E\subseteq\mathbb C^d$ est un
sous-espace $\sigma$-invariant autre que $0,\mathbb C^d$, $E^\perp$ l'est aussi (lemme 4.1(a)), et le projecteur
orthogonal sur $E$ est un élément non scalaire de $\sigma(\mathcal A)'_h$. Dans ce cas, $\mathbb C1\cap\sigma(\mathcal A)'$
et $\mathbb R1\cap\sigma(\mathcal A)'_h$ sont des sous-espaces réels propres, leurs images réciproques par les applications
$\mathbb R$-linéaires surjectives $X\mapsto Q^*XQ$ sont des sous-espaces propres de $\mathcal A'$ et de $\mathcal A'_h$, et
ceux-ci sont négligeables (corollaire A.7(c)).

(b) Comme en (a), $aQ_j=Q_j\sigma_j(a)$ et $Q_i^*a=\sigma_i(a)Q_i^*$,
donc $T\sigma_j(a)=Q_i^*YaQ_j=Q_i^*aYQ_j=\sigma_i(a)T$. Le lemme 4.2 donne $T=\lambda U$ et
$U\sigma_j(a)U^*=\sigma_i(a)$, qui vaut $Q_j'^*aQ_j'$. Enfin $P_c\in\mathcal A\cap\mathcal A'$, $Q_i=P_cQ_i$ et
$P_cQ_j=0$, donc $T=Q_i^*P_cYQ_j=Q_i^*YP_cQ_j=0$. $\square$

*Calcul :* le lemme 4.6 et le lemme 4.8(b) pour un élément aléatoire de $\pi(S_3)'$ — [tutoriel, section 4.3](../../tutorial/fr/04_symmetries.md#s-4-3).

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

<a id="r-4-9"></a>

**Proposition 4.9 (certificat).** Soit $\mathcal A$ engendrée par $S_1,\dots,S_J$ et leurs adjointes, et
$\mathcal B=V(\bigoplus_cM_{k_c}\otimes1_{m_c})V^*$ avec $V$ unitaire. Si tous les $S_j\in\mathcal B$ et
$\dim\mathcal A'=\sum_cm_c^2$, alors $\mathcal A=\mathcal B$.

*Démonstration.* $\mathcal B$ est une $*$-algèbre contenant tous les $S_j,S_j^*$, donc $\mathcal A\subseteq\mathcal B$
et $\mathcal B'\subseteq\mathcal A'$. Le théorème 4.3, appliqué à la représentation de $\mathcal B$ sur $\mathbb C^n$
(les applications de bloc $b\mapsto b_c$ sont irréductibles, et non équivalentes, car le $c$-ième projecteur central
est envoyé sur $1$ par l'une et sur $0$ par les autres), donne $\dim\mathcal B'=\sum_cm_c^2=\dim\mathcal A'$. Donc
$\mathcal B'=\mathcal A'$ et $\mathcal B=\mathcal B''=\mathcal A''=\mathcal A$. $\square$

*Dans la bibliothèque :* l'algorithme, `decompose.wedderburn` (interfaces `from_group`, `from_generators`), le certificat de la proposition 4.9, `Decomposition.verify()`, et leur coût — [tutoriel, section 4.3](../../tutorial/fr/04_symmetries.md#s-4-3).

<a id="s-4-4"></a>

## 4.4 $U(1)$ et secteurs de charge

<a id="r-4-10"></a>

**Proposition 4.10.** Soit $Q$ une matrice hermitienne sur $\mathbb C^n$ à valeurs propres entières $q$ (les
**charges**), de projecteurs spectraux $P_q$, avec $d_q=\mathrm{rank}\,P_q$, et soit $U_\theta=e^{i\theta Q}$ ; on
note $U(1)'=\{U_\theta:\theta\in\mathbb R\}'$ et $U(1)''$ son commutant.

1. $\theta\mapsto U_\theta$ est une représentation unitaire de $U(1)=\mathbb R/2\pi\mathbb Z$, et toute représentation
   unitaire continue de $U(1)$ sur $\mathbb C^n$ est de cette forme.
2. Pour $X\in M_n$, les assertions suivantes sont équivalentes : (i) $XU_\theta=U_\theta X$ pour tout $\theta$ ;
   (ii) $[X,Q]=0$ ; (iii) $X=\sum_qP_qXP_q$. Par conséquent $U(1)'\cong\bigoplus_qM_{d_q}(\mathbb C)$ (les opérateurs
   diagonaux par blocs) et $U(1)''=\mathrm{span}\{P_q\}$ : tous les $k_c$ valent $1$, et les multiplicités sont
   les $d_q$.
3. (*Charge conservée.*) Pour $H$ hermitien, $e^{itH}Qe^{-itH}=Q$ pour tout $t$ si et seulement si $[H,Q]=0$,
   c'est-à-dire si et seulement si $H\in U(1)'$. Alors chaque poids de secteur $\mathrm{Tr}(P_q\rho(t))$ de
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
$\bigoplus_q1_{d_q}\otimes M_1=\mathrm{span}\{P_q\}$ : un opérateur qui commute avec eux commute avec chaque $P_q$,
donc il est diagonal par blocs, et chaque bloc commute avec $M_{d_q}$, donc est un scalaire (lemme 2.5, étape 5, avec
$m=1$). (3) $\frac d{dt}e^{itH}Qe^{-itH}=ie^{itH}[H,Q]e^{-itH}$ (proposition 8.2) s'annule pour tout $t$ si
et seulement si $[H,Q]=0$. Alors $e^{itH}$ commute avec $P_q$ (théorème 1.5(d), deux fois), et
$\mathrm{Tr}(P_q\rho(t))=\mathrm{Tr}(e^{itH}P_qe^{-itH}\rho)=\mathrm{Tr}(P_q\rho)$. $\square$

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
qui est le produit fusionné de $n$ copies de $M_1\oplus M_1$ avec les charges $0,1$
(bit $1$ = spin haut, $\hat N=S^z+\frac n2$). Tout hamiltonien XXZ
$\sum_{\langle ij\rangle}J(s^x_is^x_j+s^y_is^y_j+\Delta s^z_is^z_j)$ conserve $\hat N$ : $s^z_is^z_j$ est diagonal dans la base
$|b\rangle$, et le terme de saut $s^x_is^x_j+s^y_is^y_j=\frac12(s^+_is^-_j+s^-_is^+_j)$ (développer $s^\pm=s^x\pm is^y$ ;
les termes croisés se compensent) envoie $|b\rangle$ sur $0$ ou sur la configuration dans laquelle un spin haut est passé
de l'un des sites $i,j$ à l'autre, donc il envoie chaque sous-espace propre de $\hat N$ dans lui-même.

*Calcul :* la proposition 4.10 pour quatre spins : les secteurs de charge trouvés à partir de $\hat N$ seul, et un anneau XXZ diagonalisé secteur par secteur — [tutoriel, section 4.4](../../tutorial/fr/04_symmetries.md#s-4-4). *Dans la bibliothèque :* $U(1)'$ pour $n$ qubits est `SpinChain(n).algebra`, et l'hamiltonien XXZ est `SpinChain.xxz` — [tutoriel, section 4.4](../../tutorial/fr/04_symmetries.md#s-4-4).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). La symétrie $U(1)$ est la conservation d'un nombre de particules, d'une charge ou
> d'une aimantation. La proposition 4.10(3) est le théorème de Noether dans ce cadre : la charge est conservée si et
> seulement si le hamiltonien est symétrique, et alors la dynamique ne quitte jamais un secteur (théorème 8.3). Si
> de plus rien ne peut mesurer les cohérences entre secteurs, les secteurs sont des secteurs de supersélection
> (chapitre 3).

<a id="s-4-5"></a>

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
$V=W\mathrm{diag}(e^{i\varphi},e^{-i\varphi})W^*$ avec $W$ unitaire, et $V=e^{iA}$ pour la matrice hermitienne de trace
nulle $A=W\mathrm{diag}(\varphi,-\varphi)W^*$. Une matrice hermitienne $2\times2$ de trace nulle est une combinaison
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

<a id="r-4-11"></a>

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

*Dans la bibliothèque :* le spin $s$ est `decompose.spin_matrices(s)`, vérifié par rapport au théorème 4.11 dans [tutoriel, section 4.6](../../tutorial/fr/04_symmetries.md#s-4-6) ; le spin total de $L$ spins $\frac12$ est `decompose.total_spin(L)` — [tutoriel, section 4.5](../../tutorial/fr/04_symmetries.md#s-4-5).

<a id="r-4-12"></a>

**Corollaire 4.12.** Toute représentation de spin sur $\mathbb C^n$ est unitairement équivalente à
$\bigoplus_s(\text{spin }s)\otimes1_{m_s}$, et, avec $d_M=\dim\ker(S_z-M)$,
$$
m_s=d_s-d_{s+1}.
$$

*Démonstration.* Théorème 4.3 et théorème 4.11. Chaque copie du spin $s'$ contient exactement un vecteur propre de
$S_z$ pour la valeur propre $M$ si $|M|\le s'$ et $s'-M\in\mathbb Z$, et aucun sinon ; donc
$d_M=\sum_{s'\ge|M|,\,s'-M\in\mathbb Z}m_{s'}$, et la différence de deux valeurs consécutives est $m_s$. $\square$

<a id="r-4-13"></a>

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
de dimensions $35$ et $14$ — les nombres trouvés au chapitre 10.

*Calcul :* la décomposition pour $L=4$ à partir des trois seuls générateurs, avec la proposition 4.13 et $\mathbf S^2=S(S+1)$ sur le spin $S$ — [tutoriel, section 4.5](../../tutorial/fr/04_symmetries.md#s-4-5).

Quels opérateurs commutent avec toutes les rotations ? Les permutations des sites, en tout cas : la conjugaison par
$\pi(\sigma)$ envoie $s^{(i)}_a$ sur $s^{(\sigma(i))}_a$ (lemme B.1(b)), donc elle fixe $\sum_is^{(i)}_a$. Réciproquement,
rien d'autre ne le fait :

<a id="r-4-14"></a>

**Théorème 4.14 (dualité de Schur–Weyl).** Soit $\pi$ la représentation de $S_L$ qui permute les facteurs de
$(\mathbb C^2)^{\otimes L}$. (a) $\mathrm{span}\,\pi(S_L)\subseteq\{S_x,S_y,S_z\}'$. (b) L'égalité a lieu ;
de manière équivalente (théorème 2.2), l'algèbre engendrée par le spin total est $\pi(S_L)'$.

*Démonstration.* La partie (a) vient d'être démontrée. La partie (b) est le corollaire B.7, le cas $d=2$ du
théorème B.6, qui démontre la dualité pour $(\mathbb C^d)^{\otimes L}$ et l'action diagonale de $U(d)$ pour tous $d$ et
$L$. $\square$

(Pour aller plus loin : W. Fulton et J. Harris, *Representation Theory*, Springer 1991, leçon 6 ; R. Goodman et
N. Wallach, *Symmetry, Representations, and Invariants*, Springer 2009, section 9.1.)

*Calcul :* le théorème 4.14 pour $L\le5$, les dimensions $2,5,14,42$ de $\{S_a\}'$ (corollaire B.10), et $\pi(S_4)''$, dont les couples $(k_c,m_c)$ sont ceux du spin total, rôles échangés (corollaire B.9) — [tutoriel, section 4.5](../../tutorial/fr/04_symmetries.md#s-4-5).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Un multiplet de spin $S$ est formé de $2S+1$ états de même énergie pour tout
> hamiltonien invariant par rotation, distingués par $S^z=M$. Deux spins $\frac12$ forment un triplet et un
> singulet ; quatre spins forment un quintuplet, trois triplets et deux singulets. Les espaces de multiplicité des
> deux singulets et des trois triplets sont le qubit sans décohérence et le qutrit sans bruit du chapitre 10.

<a id="s-4-6"></a>

## 4.6 Systèmes de spins

Le **modèle de Heisenberg** sur $L$ spins $\frac12$ est $H=\sum_{i<j}J_{ij}\,\mathbf s^{(i)}\cdot\mathbf s^{(j)}$ avec
des couplages réels $J_{ij}$.

<a id="r-4-15"></a>

**Proposition 4.15.** (a) (*Identité de Dirac*) $\mathbf s^{(i)}\cdot\mathbf s^{(j)}=\frac12F_{ij}-\frac14$, où $F_{ij}$
échange les spins $i$ et $j$. Par conséquent tout hamiltonien de Heisenberg appartient à
$\mathrm{span}\,\pi(S_L)\subseteq\{S_a\}'$, et $V^*HV=\bigoplus_S1_{2S+1}\otimes H_S$ avec des blocs $H_S\in M_{m_S}$
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

*Calcul :* la proposition 4.15 sur un anneau de quatre spins et sur un anneau de huit spins — [tutoriel, section 4.6](../../tutorial/fr/04_symmetries.md#s-4-6).

La partie (b) est la manière dont le module `su2` (chaînes de spins, chapitre 9) construit $H_S$ sans décomposition
numérique : il construit une base orthonormée des vecteurs de plus haut poids, en ajoutant un site après l'autre,
avec les coefficients de Clebsch–Gordan de $j\otimes\frac12$. Nous les obtenons en descendant à partir du plus haut
poids. D'après le théorème 4.11 et $S_-=S_+^*$, le spin $j$ vérifie
$$
S_\pm|j,m\rangle=\sqrt{(j\mp m)(j\pm m+1)}\;|j,m\pm1\rangle ,
$$
et nous posons $|j,m\rangle=0$ pour $|m|>j$ (ce qui est cohérent avec la formule, dont le coefficient s'annule en
$m=\pm j$).

<a id="r-4-16"></a>

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

*Calcul :* le théorème 4.16 pour $j\le2$, sur des matrices denses et par comparaison avec `su2.clebsch_gordan` — [tutoriel, section 4.6](../../tutorial/fr/04_symmetries.md#s-4-6).

L'exercice 2 étend la décomposition (mais pas les coefficients) à $j_1\otimes j_2$. Pour la bibliothèque, seuls
comptent les deux vecteurs de plus haut poids du théorème 4.16, et ils fonctionnent dans toute représentation de spin,
irréductible ou non :

<a id="r-4-17"></a>

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
intermédiaires $S_1=\frac12,S_2,\dots,S_L=S$, $|S_n-S_{n-1}|=\frac12$, donne la base de Yamanouchi.

*Dans la bibliothèque :* `su2.SU2Chain` construit la base de Yamanouchi avec ces coefficients, à l'intérieur du secteur $U(1)$ $\hat N=\frac L2+S$ de `SpinChain` ; le corollaire 4.17 est vérifié pour huit spins, avec les blocs $U(1)$ et $SU(2)$ de l'anneau de Heisenberg et leur coût — [tutoriel, section 4.6](../../tutorial/fr/04_symmetries.md#s-4-6).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Le modèle de Heisenberg décrit
> l'interaction d'échange entre spins électroniques localisés dans les isolants magnétiques ; $J_{ij}>0$ favorise les
> voisins antiparallèles. Son spectre est formé de multiplets, et un champ magnétique $h\,S^z$ scinde chaque multiplet
> en $2S+1$ niveaux équidistants (exercice 4). Le théorème 9.9 (couplages bipartis ; pour les chaînes et anneaux
> antiferromagnétiques de longueur paire, un état fondamental singulet) remonte à E. Lieb et D. Mattis
> (*J. Math. Phys.* 3 (1962) 749).

<a id="s-4-7"></a>

## 4.7 Algèbres réelles : type réel, complexe et quaternionique

Les symétries sont souvent données par des matrices réelles : les rotations de $\mathbb R^3$, les permutations, des
hamiltoniens réels dans une certaine base. L'algèbre complexe qu'elles engendrent se décompose par le théorème 2.6,
mais les matrices *réelles* qu'elle contient forment une algèbre réelle, dont les composantes simples sont de trois
sortes, $M_k(\mathbb R)$, $M_k(\mathbb C)$ et $M_k(\mathbb H)$. Nous les déduisons de la décomposition complexe et
d'une donnée supplémentaire : l'effet de la conjugaison complexe sur les secteurs.

Pour $a\in M_n(\mathbb C)$, soit $\bar a$ la conjuguée complexe coefficient par coefficient (de même pour les
vecteurs, et $\bar S=\{\bar s:s\in S\}$ pour les ensembles). La conjugaison est antilinéaire et bijective,
$\overline{ab}=\bar a\bar b$, $\overline{a^*}=\bar a^*$ et $\bar{\bar a}=a$ ; $\bar a=a$ si et seulement si $a$ est
réelle, et alors $a^*=a^T$. Pour un ensemble $\mathcal B$ de matrices (ou de vecteurs), soit
$\mathcal B_{\mathbb R}=\mathcal B\cap M_n(\mathbb R)$ (ou $\mathcal B\cap\mathbb R^n$). Le produit scalaire est
$\langle u,v\rangle=\sum_i\bar u_iv_i$, de sorte que
$\langle\bar u,\bar v\rangle=\overline{\langle u,v\rangle}=\langle v,u\rangle$.

<a id="r-4-18"></a>

**Lemme 4.18 (parties réelles).** (a) Soit $\mathcal B$ un sous-espace complexe de $M_n(\mathbb C)$ ou de
$\mathbb C^n$ tel que $\bar{\mathcal B}=\mathcal B$. Alors
$\mathcal B=\mathcal B_{\mathbb R}\oplus i\mathcal B_{\mathbb R}$ et
$\dim_{\mathbb R}\mathcal B_{\mathbb R}=\dim_{\mathbb C}\mathcal B$.
(b) Soit $\mathcal A$ une $*$-algèbre sur $\mathbb C^n$ telle que $\bar{\mathcal A}=\mathcal A$. Alors
$\mathcal A_{\mathbb R}$ est une algèbre réelle unifère, stable par transposition ;
$\overline{\mathcal A'}=\mathcal A'$ ; et $(\mathcal A')_{\mathbb R}$ est l'ensemble des $X$ réelles telles que
$Xa=aX$ pour tout $a\in\mathcal A_{\mathbb R}$.
(c) Si $S\subseteq M_n(\mathbb C)$ vérifie $S^*=S$ et $\bar S=S$ — par exemple $S=\{S_1,\dots,S_J,S_1^T,\dots,S_J^T\}$
avec des $S_j$ réelles, ou un groupe de matrices orthogonales réelles —, alors $\overline{S''}=S''$.

*Démonstration.* (a) Pour $b\in\mathcal B$, $\mathrm{Re}\,b=\frac12(b+\bar b)$ et
$\mathrm{Im}\,b=\frac1{2i}(b-\bar b)$ sont des éléments réels de $\mathcal B$ avec
$b=\mathrm{Re}\,b+i\mathrm{Im}\,b$, et un élément de $\mathcal B_{\mathbb R}\cap i\mathcal B_{\mathbb R}$ est
réel et imaginaire pur, donc nul. Une base $b_1,\dots,b_r$ de l'espace réel $\mathcal B_{\mathbb R}$ engendre donc
$\mathcal B$ sur $\mathbb C$, et elle est libre sur $\mathbb C$ : si $\sum_j(\alpha_j+i\beta_j)b_j=0$ avec
$\alpha_j,\beta_j$ réels, les parties réelle et imaginaire donnent $\sum_j\alpha_jb_j=\sum_j\beta_jb_j=0$. (b)
$\mathcal A_{\mathbb R}$ contient $1$ et est stable par produit, ainsi que par $a\mapsto a^T=a^*$. $Y$ commute avec
tous les $t\in T$ si et seulement si $\bar Y$ commute avec tous les $\bar t$, donc $\overline{T'}=\bar T'$ pour tout
ensemble $T$, et $\overline{\mathcal A'}=\bar{\mathcal A}'=\mathcal A'$. Une $X$ réelle qui commute avec
$\mathcal A_{\mathbb R}$ commute avec $\mathcal A_{\mathbb R}+i\mathcal A_{\mathbb R}=\mathcal A$ (d'après (a)). (c)
Appliquer deux fois $\overline{T'}=\bar T'$ : $\overline{S''}=\bar S''=S''$. $\square$

**Notations fixées.** Dans toute la suite de cette section, $\mathcal A$ est une $*$-algèbre sur $\mathbb C^n$ telle
que $\bar{\mathcal A}=\mathcal A$, et $V$ est une décomposition comme dans le théorème 2.6. $V_c$ est le bloc de
colonnes $n\times k_cm_c$ de $V$ qui correspond au $c$-ième terme, de sorte que $V_c^*V_d=\delta_{cd}1$ et
$V_cV_c^*=P_c$, les projecteurs centraux minimaux, et
$$
aV_c=V_c\big(\pi_c(a)\otimes1_{m_c}\big)\quad(a\in\mathcal A),\qquad XV_c=V_c\big(1_{k_c}\otimes x_c\big)\quad(X\in\mathcal A'),
$$
avec $\pi_c(a)\in M_{k_c}$ et $x_c\in M_{m_c}$. Chaque $\pi_c$ est un $*$-morphisme unifère de $\mathcal A$ sur
$M_{k_c}$ : pour $x\in M_{k_c}$, l'élément $V_c(x\otimes1)V_c^*$ de $\mathcal A$ a $\pi_c=x$ et $\pi_e=0$ pour
$e\neq c$. Enfin $\Omega_2=\begin{pmatrix}0&-1\\1&0\end{pmatrix}$ et $\Omega_{2r}=1_r\otimes\Omega_2$ (diagonale par
blocs avec $r$ copies de $\Omega_2$).

Une application $\theta:E\to E$ sur un sous-espace $E\subseteq\mathbb C^k$ est **antiunitaire** si
$\theta(\alpha u+\beta v)=\bar\alpha\theta u+\bar\beta\theta v$ et
$\langle\theta u,\theta v\rangle=\langle v,u\rangle$. Pour une matrice unitaire $J\in U(k)$, $\theta_Jv=J\bar v$ est
antiunitaire sur $\mathbb C^k$ ($\langle J\bar u,J\bar v\rangle=\langle\bar u,\bar v\rangle=\langle v,u\rangle$), avec
$\theta_J^2v=J\bar Jv$.

<a id="r-4-19"></a>

**Lemme 4.19 (applications antiunitaires).** Soit $\theta$ antiunitaire sur $E\neq0$ avec
$\theta^2=\varepsilon\,\mathrm{id}$, $\varepsilon\in\{1,-1\}$.
(a) Si $\varepsilon=1$, $E$ admet une base orthonormée de vecteurs $w$ tels que $\theta w=w$.
(b) Si $\varepsilon=-1$, alors $\langle u,\theta u\rangle=0$ pour tout $u\in E$ ; pour tout vecteur unitaire $u$,
$\mathrm{span}\{u,\theta u\}$ est un plan $\theta$-invariant ; et $E$ admet une base orthonormée
$w_1,\theta w_1,\dots,w_r,\theta w_r$. En particulier, $\dim E=2r$ est pair.
(c) Soit $J\in U(k)$ avec $J\bar J=\varepsilon1$. Si $\varepsilon=1$, il existe une matrice unitaire $W$ telle que
$J\bar W=W$, c'est-à-dire $J=WW^T$. Si $\varepsilon=-1$, $k$ est pair et il existe une matrice unitaire $W$ telle que
$J\bar W=W\Omega_k$, c'est-à-dire $J=W\Omega_kW^T$.

*Démonstration.* (a) Soit $F=\{v\in E:\theta v=v\}$, un sous-espace réel. Tout $v\in E$ s'écrit $v_1+iv_2$ avec
$v_1=\frac12(v+\theta v)$ et $v_2=\frac1{2i}(v-\theta v)$, et tous deux sont dans $F$
($\theta v_2=\frac1{-2i}(\theta v-v)=v_2$). Pour $u,v\in F$,
$\langle u,v\rangle=\langle\theta u,\theta v\rangle=\langle v,u\rangle$ est réel, donc $\langle\cdot,\cdot\rangle$ est
un produit scalaire réel sur $F$ ; soit $w_1,\dots,w_s$ une base orthonormée de $F$ (Gram–Schmidt). Les $w_j$ sont
orthonormés dans $E$, donc libres sur $\mathbb C$, et ils engendrent $E$ sur $\mathbb C$ puisque $E=F+iF$. (b) En
appliquant $\langle\theta x,\theta y\rangle=\langle y,x\rangle$ à $x=\theta u$, $y=u$ :
$\langle u,\theta u\rangle=\langle\theta^2u,\theta u\rangle=-\langle u,\theta u\rangle$. Ainsi, pour un vecteur
unitaire $u$, $\theta u$ est un vecteur unitaire orthogonal à $u$, et $E_1=\mathrm{span}\{u,\theta u\}$ est
invariant : $\theta(\alpha u+\beta\theta u)=\bar\alpha\theta u-\bar\beta u$. L'orthogonal de $E_1$ dans $E$ est
invariant lui aussi : pour $x\perp E_1$ et $y\in E_1$, $y=-\theta(\theta y)$ et
$\langle\theta x,y\rangle=-\langle\theta x,\theta(\theta y)\rangle=-\langle\theta y,x\rangle=0$, puisque
$\theta y\in E_1$. Une récurrence sur $\dim E$ donne la base. (c) $\theta_J$ est antiunitaire sur $\mathbb C^k$ avec
$\theta_J^2=J\bar J=\varepsilon$. Soit $W$ la matrice unitaire dont les colonnes sont la base de (a), resp. la base
$w_1,\theta_Jw_1,w_2,\theta_Jw_2,\dots$ de (b). Les colonnes de $J\bar W$ sont les images par $\theta_J$ des colonnes
de $W$. Dans le cas (a), cela dit $J\bar W=W$. Dans le cas (b), les images sont $\theta_Jw_j$ et $\theta_J^2w_j=-w_j$,
qui sont les colonnes $2j-1$ et $2j$ de $W\Omega_k$ (car $\Omega_2e_1=e_2$, $\Omega_2e_2=-e_1$). Comme
$\bar W^{-1}=\overline{W^*}=W^T$, $J=WW^T$, resp. $J=W\Omega_kW^T$. $\square$

<a id="r-4-20"></a>

**Théorème 4.20 (conjugaison des secteurs).** (a) Il existe une permutation $\sigma$ des secteurs avec
$\sigma\circ\sigma=\mathrm{id}$ et des matrices unitaires $J_c\in U(k_c)$, $M_c\in U(m_c)$ telles que
$\bar P_c=P_{\sigma(c)}$, $k_{\sigma(c)}=k_c$, $m_{\sigma(c)}=m_c$ et, pour tout $c$ et tout $a\in\mathcal A$,
$$
\overline{V_c}=V_{\sigma(c)}\,(J_c\otimes M_c),\qquad\pi_{\sigma(c)}(\bar a)=J_c\,\overline{\pi_c(a)}\,J_c^*.
$$
$J_c\otimes M_c$ est unique, et $J_c,M_c$ sont uniques à $(J_c,M_c)\mapsto(\lambda J_c,\bar\lambda M_c)$ près, avec
$|\lambda|=1$.
(b) Si $\sigma(c)=c$, il existe $\varepsilon_c\in\{1,-1\}$ tel que $J_c\bar J_c=\varepsilon_c1_{k_c}$ et
$M_c\bar M_c=\varepsilon_c1_{m_c}$. Si $\varepsilon_c=-1$, alors $k_c$ et $m_c$ sont pairs.
(c) $\sigma$ et les $\varepsilon_c$ ne dépendent pas du choix de $V$.
(d) Pour $y\in M_{m_c}$,
$\overline{V_c(1\otimes y)V_c^*}=V_{\sigma(c)}\big(1\otimes M_c\bar yM_c^*\big)V_{\sigma(c)}^*$. Par conséquent, si
$X\in\mathcal A'$ a les blocs $x_c$, alors $\bar X\in\mathcal A'$ a les blocs
$(\bar X)_{\sigma(c)}=M_c\bar x_cM_c^*$ : sur le commutant, $M_c$ joue le rôle que $J_c$ joue sur $\mathcal A$.

*Démonstration.* (a) La conjugaison envoie $\mathcal A$ sur $\mathcal A$ et préserve sommes, produits et adjoints ;
elle envoie donc le centre $Z(\mathcal A)$ sur lui-même, et l'ensemble des projecteurs non nuls de $Z(\mathcal A)$
bijectivement sur lui-même, en préservant l'ordre $p\le q$ ($pq=p$). Elle permute donc les projecteurs minimaux, qui
sont les $P_c$ (lemme 2.3) : $\bar P_c=P_{\sigma(c)}$, et $\sigma\circ\sigma=\mathrm{id}$ car $\bar{\bar P}_c=P_c$.
Soit $d=\sigma(c)$. L'image de $\overline{V_c}$ est la conjuguée de l'image de $P_c$, c'est-à-dire l'image de
$\bar P_c=P_d$ ; donc $\overline{V_c}=P_d\overline{V_c}=V_dT$ avec $T=V_d^*\overline{V_c}$, et
$T^*T=\overline{V_c}^*P_d\overline{V_c}=\overline{V_c^*V_c}=1$. Comme
$k_cm_c=\mathrm{rank}\,P_c=\mathrm{rank}\,\bar P_c=k_dm_d$, l'isométrie $T$ est carrée, donc unitaire. En
conjuguant $aV_c=V_c(\pi_c(a)\otimes1)$ et en utilisant $\bar aV_d=V_d(\pi_d(\bar a)\otimes1)$, on obtient
$V_d(\pi_d(\bar a)\otimes1)T=V_dT(\overline{\pi_c(a)}\otimes1)$ ; en multipliant par $V_d^*$,
$$
\big(\pi_d(\bar a)\otimes1_{m_d}\big)\,T=T\,\big(\overline{\pi_c(a)}\otimes1_{m_c}\big)\qquad(a\in\mathcal A).
$$
Écrivons $T=\sum_{j,l}T_{jl}\otimes|j\rangle\langle l|$ avec des blocs $T_{jl}$ de taille $k_d\times k_c$. L'égalité
ci-dessus dit $\pi_d(\bar a)T_{jl}=T_{jl}\overline{\pi_c(a)}$ : chaque $T_{jl}$ est un opérateur d'entrelacement de
$\rho(a)=\overline{\pi_c(a)}$ vers $\rho'(a)=\pi_d(\bar a)$. Ce sont des représentations de $\mathcal A$ muni de son
adjoint au sens de la section 4.1 ($\rho(a^*)=\rho(a)^*$, $\rho'(a^*)=\rho'(a)^*$), et elles sont irréductibles,
puisque leurs images sont $M_{k_c}$ et $M_{k_d}$ tout entiers. Comme $T\neq0$, un $T_{jl}$ est non nul, donc $k_c=k_d$
et $\mathrm{Hom}(\rho,\rho')=\mathbb CJ_c$ pour une matrice unitaire $J_c$ (lemme 4.2). Ainsi
$T_{jl}=\mu_{jl}J_c$, c'est-à-dire $T=J_c\otimes M_c$ avec $M_c=(\mu_{jl})$ ; $1=T^*T=1\otimes M_c^*M_c$ et
$m_c=k_cm_c/k_c=m_d$, donc $M_c$ est unitaire, et l'égalité devient la formule pour $\pi_d(\bar a)$. $T$ est
déterminée par $V$ ; si l'on a aussi $T=J'\otimes M'$ avec $J'$ unitaire, les blocs $T_{jl}=M'_{jl}J'$ montrent que
$J'\in\mathrm{Hom}(\rho,\rho')$, donc $J'=\lambda J_c$ avec $|\lambda|=1$, et alors $M'=\bar\lambda M_c$.

(b) Soit $\sigma(c)=c$ et $J=J_c$, $M=M_c$. En conjuguant $\overline{V_c}=V_c(J\otimes M)$, on obtient
$V_c=\overline{V_c}(\bar J\otimes\bar M)=V_c(J\bar J\otimes M\bar M)$, donc $J\bar J\otimes M\bar M=1$. En appliquant
deux fois la formule de (a), avec $\overline{J^*}=J^T$,
$$
\pi_c(a)=\pi_c(\bar{\bar a})=J\,\overline{J\,\overline{\pi_c(a)}\,J^*}\,J^*=(J\bar J)\,\pi_c(a)\,(J\bar J)^*,
$$
donc la matrice unitaire $J\bar J$ commute avec $\pi_c(\mathcal A)=M_{k_c}$ et est un scalaire $\lambda1$ (lemme 2.5,
étape 5, avec $m=1$), $|\lambda|=1$. Alors $\bar JJ=\overline{J\bar J}=\bar\lambda1$ et
$\lambda J=(J\bar J)J=J(\bar JJ)=\bar\lambda J$, donc $\lambda=\bar\lambda\in\{1,-1\}$ ; notons-le $\varepsilon_c$. Il
ne change pas par $J\mapsto\lambda J$, et $\varepsilon_c1\otimes M\bar M=1$ donne $M\bar M=\varepsilon_c1$. Si
$\varepsilon_c=-1$, alors $(-1)^{k_c}=\det(J\bar J)=\det J\cdot\overline{\det J}=|\det J|^2>0$, donc $k_c$ est pair,
et de même $m_c$.

(c) Les $P_c$ ne dépendent pas de $V$ (théorème 2.6), donc $\sigma$ non plus. Soit $V'$ une autre décomposition,
numérotée de sorte que $V'_cV_c'^*=P_c$ (théorème 2.6), avec des morphismes $\pi'_c$ et des données $J'_c,M'_c$, et
soit $\sigma(c)=c$. Alors $V'_c=P_cV'_c=V_cY$ avec $Y=V_c^*V'_c$, une isométrie carrée, donc unitaire, et
$(\pi_c(a)\otimes1)Y=Y(\pi'_c(a)\otimes1)$ pour $a\in\mathcal A$. L'argument de (a), avec les représentations
irréductibles $\pi'_c,\pi_c$ à la place de $\rho,\rho'$, donne $Y=u\otimes w$ avec des matrices unitaires $u,w$. Par
conséquent
$\overline{V'_c}=\overline{V_c}(\bar u\otimes\bar w)=V_c(J_c\bar u\otimes M_c\bar w)=V'_c(u^*J_c\bar u\otimes w^*M_c\bar w)$,
donc $J'_c=\lambda u^*J_c\bar u$ avec $|\lambda|=1$ par l'unicité dans (a), et, comme $\bar uu^T=\overline{uu^*}=1$,
$J'_c\bar J'_c=u^*J_c\bar uu^T\bar J_cu=u^*J_c\bar J_cu=\varepsilon_c1$.

(d) D'après (a),
$\overline{V_c(1\otimes y)V_c^*}=V_{\sigma(c)}(J_c\otimes M_c)(1\otimes\bar y)(J_c\otimes M_c)^*V_{\sigma(c)}^*$, et
$(J_c\otimes M_c)(1\otimes\bar y)(J_c\otimes M_c)^*=J_cJ_c^*\otimes M_c\bar yM_c^*=1\otimes M_c\bar yM_c^*$. Sommer
sur $c$ pour $X=\sum_cV_c(1\otimes x_c)V_c^*$. $\square$

Un secteur $c$ est de **type réel** si $\sigma(c)=c$ et $\varepsilon_c=1$, de **type quaternionique** si $\sigma(c)=c$
et $\varepsilon_c=-1$, et de **type complexe** si $\sigma(c)\neq c$ ; d'après le théorème 4.20(c), c'est une propriété
de $\mathcal A$ et de $P_c$. Un secteur avec $\sigma(c)=c$ et $k_c$ impair ou $m_c$ impair est de type réel, d'après
le théorème 4.20(b).

*Dans la bibliothèque :* $\sigma$, $J_c$, $M_c$, $\varepsilon_c$ et les types des secteurs sont calculés par `Decomposition.real_structure()` — [tutoriel, section 4.7](../../tutorial/fr/04_symmetries.md#s-4-7).

<a id="r-4-21"></a>

**Théorème 4.21 (l'algèbre réelle).** Soit $\Gamma$ un ensemble de secteurs contenant tout secteur $c$ tel que
$\sigma(c)=c$ et un secteur de chaque paire $\{c,\sigma(c)\}$ avec $\sigma(c)\neq c$. Pour $\sigma(c)=c$, soit
$R_c=\{x\in M_{k_c}(\mathbb C):xJ_c=J_c\bar x\}$. Alors $a\mapsto(\pi_c(a))_{c\in\Gamma}$ est un isomorphisme
d'algèbres réelles, qui envoie la transposée sur l'adjoint,
$$
\mathcal A_{\mathbb R}\;\cong\bigoplus_{c\in\Gamma,\ \sigma(c)=c}R_c\;\oplus\bigoplus_{c\in\Gamma,\ \sigma(c)\neq c}M_{k_c}(\mathbb C),
$$
les derniers termes étant vus comme des algèbres réelles (**type complexe**). De plus :
(a) (*type réel*) Si $\varepsilon_c=1$, écrivons $J_c=W_JW_J^T$ et $M_c=W_MW_M^T$ avec $J_c\bar W_J=W_J$,
$M_c\bar W_M=W_M$ (lemme 4.19(c)). Alors $x\mapsto W_J^*xW_J$ envoie $R_c$ isomorphiquement sur $M_{k_c}(\mathbb R)$.
Les colonnes de $O_c=V_c(W_J\otimes W_M)$ sont réelles et orthonormées, elles engendrent l'image de $P_c$ sur
$\mathbb C$, et $aO_c=O_c(\alpha_c(a)\otimes1_{m_c})$ avec la matrice réelle $\alpha_c(a)=W_J^*\pi_c(a)W_J$ pour
$a\in\mathcal A_{\mathbb R}$.
(b) (*type quaternionique*) Si $\varepsilon_c=-1$, écrivons $J_c=W_J\Omega_{k_c}W_J^T$ avec
$J_c\bar W_J=W_J\Omega_{k_c}$ (lemme 4.19(c)). Alors $x\mapsto\mathbb M^{-1}(W_J^*xW_J)$ est un isomorphisme
d'algèbres réelles $R_c\cong M_{k_c/2}(\mathbb H)$ qui envoie l'adjoint sur l'adjoint quaternionique
$(q^*)_{ij}=\bar q_{ji}$ ; ici $\mathbb M:M_r(\mathbb H)\to M_{2r}(\mathbb C)$ remplace chaque coefficient $q$ par le
bloc $2\times2$ $M(q)$ du lemme E.15(2), et son image est $\{y\in M_{2r}(\mathbb C):y\Omega_{2r}=\Omega_{2r}\bar y\}$.
En particulier, $\dim_{\mathbb R}\mathcal A_{\mathbb R}=\sum_ck_c^2=\dim_{\mathbb C}\mathcal A$.

*Démonstration.* Pour $a$ réelle, $\bar a=a$, et le théorème 4.20(a) s'écrit
$\pi_{\sigma(c)}(a)=J_c\overline{\pi_c(a)}J_c^*$. Pour $\sigma(c)=c$, cela dit que $\pi_c(a)\in R_c$. Ainsi
$\Phi(a)=(\pi_c(a))_{c\in\Gamma}$ envoie $\mathcal A_{\mathbb R}$ dans le membre de droite ; c'est un homomorphisme
d'algèbres réelles, et $\Phi(a^T)=\Phi(a^*)=\Phi(a)^*$. Si $\Phi(a)=0$, la même formule donne $\pi_{\sigma(c)}(a)=0$
pour $c\in\Gamma$, donc tous les $\pi_e(a)$ sont nuls et $a=0$. *Surjectivité :* soient $c\in\Gamma$, $x\in R_c$ si
$\sigma(c)=c$ et $x\in M_{k_c}$ sinon, et $a=V_c(x\otimes1)V_c^*\in\mathcal A$. D'après le théorème 4.20(a) pour le
secteur $\sigma(e)$, $\pi_e(\bar a)=J_{\sigma(e)}\overline{\pi_{\sigma(e)}(a)}J_{\sigma(e)}^*$, qui est nul sauf si
$\sigma(e)=c$. Si $\sigma(c)=c$, $\pi_c(\bar a)=J_c\bar xJ_c^*=x$ (définition de $R_c$), donc la matrice réelle
$\mathrm{Re}\,a=\frac12(a+\bar a)\in\mathcal A_{\mathbb R}$ a $\Phi(\mathrm{Re}\,a)=x$ dans la composante $c$
et $0$ ailleurs. Si $\sigma(c)\neq c$, $\pi_c(\bar a)=0$, et $2\mathrm{Re}\,a$ a $\pi_c=x$ et $\pi_e=0$ pour
$e\notin\{c,\sigma(c)\}$ ; comme $\sigma(c)\notin\Gamma$, $\Phi(2\mathrm{Re}\,a)$ vaut $x$ dans la composante $c$
et $0$ ailleurs. Les sommes de ces éléments donnent tout le membre de droite.

(a) Multiplions $xW_JW_J^T=W_JW_J^T\bar x$ par $W_J^*$ à gauche et par $\bar W_J$ à droite ; comme
$W_J^T\bar W_J=\overline{W_J^*W_J}=1$, cela devient $W_J^*xW_J=W_J^T\bar x\bar W_J=\overline{W_J^*xW_J}$, et les
étapes sont réversibles. Donc $R_c=W_JM_{k_c}(\mathbb R)W_J^*$, et $x\mapsto W_J^*xW_J$ est un isomorphisme d'algèbres
réelles compatible avec l'adjoint. Ensuite
$\overline{O_c}=\overline{V_c}(\bar W_J\otimes\bar W_M)=V_c(J_c\bar W_J\otimes M_c\bar W_M)=V_c(W_J\otimes W_M)=O_c$,
$O_c^*O_c=1$, et $aO_c=V_c(\pi_c(a)W_J\otimes W_M)=O_c(W_J^*\pi_c(a)W_J\otimes1)$, avec $W_J^*\pi_c(a)W_J$ réelle pour
$a\in\mathcal A_{\mathbb R}$.

(b) Le même calcul avec $J_c=W_J\Omega W_J^T$ ($\Omega=\Omega_{k_c}$) montre que $xJ_c=J_c\bar x$ si et seulement si
$y=W_J^*xW_J$ vérifie $y\Omega=\Omega\bar y$. Écrivons $y$ en blocs $2\times2$ $y_{ij}$ ; comme $\Omega$ est diagonale
par blocs, $y\Omega=\Omega\bar y$ si et seulement si $y_{ij}\Omega_2=\Omega_2\bar y_{ij}$ pour tous $i,j$, et pour
$y_{ij}=\begin{pmatrix}\alpha&\beta\\\gamma&\delta\end{pmatrix}$ cela s'écrit
$\begin{pmatrix}\beta&-\alpha\\\delta&-\gamma\end{pmatrix}=\begin{pmatrix}-\bar\gamma&-\bar\delta\\\bar\alpha&\bar\beta\end{pmatrix}$,
c'est-à-dire $\delta=\bar\alpha$ et $\beta=-\bar\gamma$ :
$y_{ij}=\begin{pmatrix}\alpha&-\bar\gamma\\\gamma&\bar\alpha\end{pmatrix}=M(\alpha,\bar\gamma)$ avec
$M(a,b)=\begin{pmatrix}a&-b\\\bar b&\bar a\end{pmatrix}$. L'image de $\mathbb M$ est donc celle annoncée. D'après le
lemme E.15(2), $M$ est un homomorphisme $\mathbb R$-linéaire injectif $\mathbb H\to M_2(\mathbb C)$ avec
$M(\bar q)=M(q)^*$ ; comme un produit matriciel n'utilise que des sommes et des produits de coefficients, $\mathbb M$
est un homomorphisme injectif d'algèbres réelles avec $\mathbb M(q^*)=\mathbb M(q)^*$ (comme dans la démonstration du
théorème E.24(1)). Par conséquent $x\mapsto\mathbb M^{-1}(W_J^*xW_J)$ est l'isomorphisme annoncé. *Dimensions :*
$\dim_{\mathbb R}R_c=k_c^2$ dans les deux cas ($\dim_{\mathbb R}M_r(\mathbb H)=4r^2$ avec $2r=k_c$), et une paire
contribue $\dim_{\mathbb R}M_{k_c}(\mathbb C)=2k_c^2=k_c^2+k_{\sigma(c)}^2$. $\square$

Il en va de même pour le commutant, les facteurs étant échangés : $(\mathcal A')_{\mathbb R}$ est la somme directe de
$M_{m_c}(\mathbb R)$, $M_{m_c/2}(\mathbb H)$ et $M_{m_c}(\mathbb C)$ sur les mêmes secteurs, avec le même type. En
effet, $\mathcal A'$ est une $*$-algèbre telle que $\overline{\mathcal A'}=\mathcal A'$ (lemme 4.18(b)), et la
démonstration ci-dessus n'a utilisé que trois faits, qui valent pour $\mathcal A'$ avec $X\mapsto x_c$ à la place de
$\pi_c$ et $M_c$ à la place de $J_c$ : $X\mapsto x_c$ est un $*$-morphisme unifère sur $M_{m_c}$
($V_c(1\otimes y)V_c^*\in\mathcal A'$ a le bloc $y$ dans le secteur $c$ et $0$ ailleurs),
$(\bar X)_{\sigma(c)}=M_c\bar x_cM_c^*$ (théorème 4.20(d)), et $M_c\bar M_c=\varepsilon_c1$ pour $\sigma(c)=c$
(théorème 4.20(b)).

*Dans la bibliothèque :* une base orthogonale réelle adaptée au théorème 4.21, `Decomposition.real_basis()` — [tutoriel, section 4.7](../../tutorial/fr/04_symmetries.md#s-4-7).

Trois exemples. (i) Les générateurs $(L_a)_{bc}=-\epsilon_{abc}$ des rotations
de $\mathbb R^3$, agissant sur $\mathbb R^3\otimes\mathbb R^2$ par $L_a\otimes1_2$ : il y a un secteur, $k=3$, $m=2$ ($iL_x,iL_y,iL_z$ est une représentation de spin sur $\mathbb C^3$ dans laquelle
$iL_z$ a les valeurs propres simples $1,0,-1$, c'est donc le spin $1$ d'après le corollaire 4.12) ; il est de type réel car $k$ est impair, et $\mathcal A_{\mathbb R}\cong M_3(\mathbb R)$. (ii)
Le générateur des rotations $\Omega_2$ de $\mathbb R^2$, agissant sur $\mathbb R^2\otimes\mathbb R^2$ par
$\Omega_2\otimes1_2$ : comme $\Omega_2^2=-1$ et $\Omega_2^*=-\Omega_2$, l'algèbre est
$\mathrm{span}\{1,\Omega_2\otimes1\}$ ; les vecteurs propres $(1,\mp i)/\sqrt2$ de $\Omega_2$ pour $\pm i$,
conjugués l'un de l'autre, donnent deux secteurs avec $k=1$, $m=2$, échangés par la conjugaison : une paire de type
complexe, et $\mathcal A_{\mathbb R}=\{\alpha+\beta\,\Omega_2\otimes1\}\cong\mathbb C$. (iii) La multiplication à
gauche $L_q$ par un quaternion $q=(a,b)\in\mathbb H$ (section E.6) sur $\mathbb H=\mathbb R^4$, avec les coordonnées
$(\mathrm{Re}\,a,\mathrm{Im}\,a,\mathrm{Re}\,b,\mathrm{Im}\,b)$. $q\mapsto L_q$ est un homomorphisme
injectif d'algèbres réelles ($\mathbb H$ est associative, lemme E.15(2)), et $L_{\bar q}=L_q^T$, car
$\langle qx,y\rangle=\langle x,\bar qy\rangle$ (lemme E.14(2)) pour la norme $N(a,b)=|a|^2+|b|^2$, qui est la norme
euclidienne de $\mathbb R^4$. L'algèbre réelle engendrée par $L_I,L_J,L_K$ ($I=(i,0)$, $J=(0,1)$, $K=(0,i)$) et leurs
transposées est donc $\{L_q:q\in\mathbb H\}\cong\mathbb H$ ; l'algèbre complexe $\mathcal A$ est l'espace engendré sur
$\mathbb C$ par les produits des générateurs, $\{L_q\}+i\{L_q\}$, donc $\mathcal A_{\mathbb R}=\{L_q\}\cong\mathbb H$
(lemme 4.18(a)). Cela fixe déjà la décomposition :
$\sum_ck_c^2=\dim_{\mathbb C}\mathcal A=\dim_{\mathbb R}\mathcal A_{\mathbb R}=4$ (théorème 4.21), et les $k_c$ ne
valent pas tous $1$, car sinon $\mathcal A$, donc $\mathcal A_{\mathbb R}\cong\mathbb H$, serait commutative ; il y a
donc un seul secteur, avec $k=2$ et $\sigma(c)=c$, et il n'est pas de type réel, car le théorème 4.21(a) donnerait
alors $\mathcal A_{\mathbb R}\cong M_2(\mathbb R)$, qui a des diviseurs de zéro, alors que $\mathbb H$ n'en a pas
(lemme E.15(3)) : il est de type quaternionique, avec $m=4/2=2$ sur $\mathbb R^4$. (Les quaternions de norme $1$
agissent ici comme $SU(2)$ : $M(q)$ est unitaire de déterminant $|a|^2+|b|^2=1$, et toute
$U=\begin{pmatrix}\alpha&\beta\\\gamma&\delta\end{pmatrix}\in SU(2)$ est de cette forme, car
$U^*=U^{-1}=\begin{pmatrix}\delta&-\beta\\-\gamma&\alpha\end{pmatrix}$ donne $\delta=\bar\alpha$, $\beta=-\bar\gamma$,
donc $U=M(\alpha,\bar\gamma)$.) Sur $\mathbb H^2=\mathbb R^8$, où $q$ agit par $1_2\otimes L_q$, le même argument
donne un secteur avec $k=2$, $m=4$, de type quaternionique, $\mathcal A_{\mathbb R}\cong M_1(\mathbb H)=\mathbb H$.

*Calcul :* les trois exemples, avec `Decomposition.real_structure()` et `Decomposition.real_basis()` — [tutoriel, section 4.7](../../tutorial/fr/04_symmetries.md#s-4-7).

Dans la théorie des représentations d'un groupe par des matrices réelles, l'énoncé suivant est la forme réelle du
lemme de Schur ; il montre que les trois types se lisent sur la représentation réelle elle-même.

<a id="r-4-22"></a>

**Théorème 4.22 (lemme de Schur réel).** Soit $U\subseteq\mathbb R^n$ un sous-espace non nul tel que $aU\subseteq U$
pour tout $a\in\mathcal A_{\mathbb R}$ et qui est **irréductible** : $0$ et $U$ sont ses seuls sous-espaces ayant
cette propriété. Soit $D_U$ l'ensemble des applications $\mathbb R$-linéaires $T:U\to U$ telles que $Tau=aTu$ pour
tous $a\in\mathcal A_{\mathbb R}$, $u\in U$.
(a) $D_U$ est une algèbre réelle associative unifère sans diviseurs de zéro, donc isomorphe à $\mathbb R$, $\mathbb C$
ou $\mathbb H$ (théorème E.19).
(b) Pour exactement une orbite $\mathcal O=\{c,\sigma(c)\}$ de $\sigma$, $U$ est contenu dans l'image de
$P_{\mathcal O}=\sum_{e\in\mathcal O}P_e$, et $\dim_{\mathbb R}D_U=1,2,4$ et $\dim_{\mathbb R}U=k_c,2k_c,2k_c$ selon
que $c$ est de type réel, complexe ou quaternionique. Ainsi $D_U\cong\mathbb R,\mathbb C,\mathbb H$ respectivement.
(c) Pour toute orbite $\mathcal O$, il existe un tel $U$ dans l'image de $P_{\mathcal O}$. $\mathbb R^n$ lui-même est
irréductible si et seulement si $\sigma$ a une seule orbite $\{c,\sigma(c)\}$ et $m_c=1$ (types réel et complexe),
resp. $m_c=2$ (type quaternionique).

*Démonstration.* (a) Pour $T\in D_U$, $\ker T$ et $TU$ sont des sous-espaces invariants, donc $T=0$ ou $T$ est
bijective, et alors $T^{-1}\in D_U$. $D_U$ est stable par composition, qui est associative, contient
$\mathrm{id}_U\neq0$, et $ST=0$ avec $S\neq0$ impose $T=S^{-1}ST=0$.

(b) *Sous-espaces invariants et projecteurs réels.* Si $U'\subseteq\mathbb R^n$ est invariant, son orthogonal dans
$\mathbb R^n$ l'est aussi (pour $y\perp U'$, $u\in U'$ et $a\in\mathcal A_{\mathbb R}$, $u\cdot ay=a^Tu\cdot y=0$, car
$a^T\in\mathcal A_{\mathbb R}$), donc le projecteur orthogonal $e$ sur $U'$, un idempotent réel symétrique, commute
avec $\mathcal A_{\mathbb R}$ ; d'après le lemme 4.18(b), $e\in(\mathcal A')_{\mathbb R}$. Réciproquement,
$e\mathbb R^n$ est invariant pour tout projecteur $e\in(\mathcal A')_{\mathbb R}$. Ainsi $U'\mapsto e$ est une
bijection croissante (pour l'inclusion) des sous-espaces invariants sur l'ensemble $\mathcal P$ des projecteurs de
$(\mathcal A')_{\mathbb R}$, et $U$ est irréductible si et seulement si son projecteur $e$ est un élément non nul
minimal de $\mathcal P$. D'après le lemme 4.18(a) appliqué à l'image $e\mathbb C^n$ (dont les vecteurs réels forment
$e\mathbb R^n$), $\dim_{\mathbb R}U=\mathrm{rank}\,e$.

*$D_U\cong e(\mathcal A')_{\mathbb R}e$.* Pour $T\in D_U$, soit $\hat T$ la matrice réelle $n\times n$
$x\mapsto T(ex)$. Pour $a\in\mathcal A_{\mathbb R}$, $\hat Tax=T(eax)=T(aex)=aT(ex)=a\hat Tx$, donc
$\hat T\in(\mathcal A')_{\mathbb R}$ (lemme 4.18(b)) et $\hat T=e\hat Te$. Réciproquement, pour
$X\in(\mathcal A')_{\mathbb R}$, la restriction de $eXe$ à $U$ appartient à $D_U$. Ces applications sont inverses
l'une de l'autre et multiplicatives. De plus, $e\mathcal A'e$ est stable par conjugaison, et ses éléments réels sont
les $e(\mathcal A')_{\mathbb R}e$ (si $eXe$ est réelle, $eXe=e\mathrm{Re}(X)e$ avec
$\mathrm{Re}\,X\in(\mathcal A')_{\mathbb R}$), donc $\dim_{\mathbb R}D_U=\dim_{\mathbb C}e\mathcal A'e$ d'après le
lemme 4.18(a).

*Projecteurs dans la décomposition.* Un projecteur $e\in\mathcal A'$ s'écrit $e=\sum_cV_c(1\otimes f_c)V_c^*$ avec des
projecteurs $f_c\in M_{m_c}$ (théorème 2.6). Alors $\mathrm{rank}\,e=\sum_ck_c\mathrm{rank}\,f_c$,
$e\mathcal A'e=\sum_cV_c(1\otimes f_cM_{m_c}f_c)V_c^*$ est de dimension complexe $\sum_c(\mathrm{rank}\,f_c)^2$
et, d'après le théorème 4.20(d), $\bar e=e$ si et seulement si $f_{\sigma(c)}=M_c\bar f_cM_c^*$ pour tout $c$.

*Éléments minimaux.* Soit $e$ minimal dans $\mathcal P$. Pour chaque orbite $\mathcal O$ de $\sigma$, $P_{\mathcal O}$
est un projecteur central réel ($\bar P_{\mathcal O}=P_{\mathcal O}$ d'après le théorème 4.20(a)) dans
$Z(\mathcal A)\subseteq\mathcal A'$, qui commute avec $e$, donc $eP_{\mathcal O}\in\mathcal P$ et
$eP_{\mathcal O}\le e$ ; comme $e=\sum_{\mathcal O}eP_{\mathcal O}$, une orbite $\mathcal O$ et une seule vérifie
$eP_{\mathcal O}=e$, c'est-à-dire que $U=e\mathbb R^n$ est dans l'image de $P_{\mathcal O}$, et $f_c=0$ pour
$c\notin\mathcal O$. Si $\mathcal O=\{c,d\}$ avec $d\neq c$, alors $f_d=M_c\bar f_cM_c^*$ a le rang de $f_c$ ; pour un
projecteur de rang un $g\le f_c$, $g'=Y+\bar Y$ avec $Y=V_c(1\otimes g)V_c^*$ est une somme de deux projecteurs
orthogonaux de $\mathcal A'$ (dans les secteurs $c$ et $d$), réelle, et $g'\le e$ (car $Y\le e$ et
$\bar Y\le\bar e=e$) ; donc $g'=e$, $\mathrm{rank}\,f_c=\mathrm{rank}\,f_d=1$,
$\dim_{\mathbb C}e\mathcal A'e=2$ et $\mathrm{rank}\,e=2k_c$. Si $\mathcal O=\{c\}$, soient $f=f_c$, $M=M_c$ et
$\theta v=M\bar v$, antiunitaire sur $\mathbb C^{m_c}$ avec $\theta^2=M\bar M=\varepsilon_c$ (théorème 4.20(b)). De
$fM=M\bar f$ on tire $\theta(fv)=M\bar f\bar v=fM\bar v=f\theta v$, donc l'image $F$ de $f$ est $\theta$-invariante.
Si $F_1\subseteq F$ est $\theta$-invariant, de projecteur $g$, alors $M\bar gM^*$ est le projecteur sur
$M\bar F_1=\theta F_1=F_1$ ($\theta$ est injective), donc $g=M\bar gM^*$ et $V_c(1\otimes g)V_c^*\in\mathcal P$ est
majoré par $e$. Par minimalité, $F$ n'a pas d'autres sous-espaces $\theta$-invariants que $0$ et $F$. Si
$\varepsilon_c=1$, un vecteur unitaire de $F$ fixé par $\theta$ (lemme 4.19(a)) en engendre un, donc
$\mathrm{rank}\,f=1$, $\dim_{\mathbb C}e\mathcal A'e=1$ et $\mathrm{rank}\,e=k_c$. Si $\varepsilon_c=-1$, le
plan $\mathrm{span}\{u,\theta u\}$ d'un vecteur unitaire $u\in F$ en est un (lemme 4.19(b)), donc
$\mathrm{rank}\,f=2$, $\dim_{\mathbb C}e\mathcal A'e=4$ et $\mathrm{rank}\,e=2k_c$. Comme $\mathbb R$,
$\mathbb C$, $\mathbb H$ sont de dimensions $1,2,4$, (a) donne $D_U\cong\mathbb R,\mathbb C,\mathbb H$ respectivement.

(c) $P_{\mathcal O}\in\mathcal P$ est non nul, et un élément de $\mathcal P\setminus\{0\}$ majoré par lui et de rang
minimal est minimal. $\mathbb R^n$ est irréductible si et seulement si $e=1$ est minimal ; ses blocs sont
$f_c=1_{m_c}$, donc, d'après la démonstration de (b), cela exige une seule orbite et $m_c=1,1,2$ dans les trois cas.
Réciproquement, sous ces conditions, un $e'\in\mathcal P$ non nul a $f'_c\neq0$ pour les secteurs de l'orbite
($f'_{\sigma(c)}=M_c\bar f'_cM_c^*$), donc $f'_c=1$ si $m_c=1$ ; dans le cas quaternionique, l'image de $f'_c$ est
$\theta$-invariante, donc de dimension paire (lemme 4.19(b)), d'où $f'_c=1_2$. Ainsi $e'=1$. $\square$

Pour un groupe $G$ agissant par des matrices orthogonales réelles, $\mathcal A=\mathrm{span}\,\pi(G)$ est stable
par conjugaison (lemme 4.18(c)), $\mathcal A_{\mathbb R}$ est l'espace engendré sur $\mathbb R$ par $\pi(G)$ (lemme
4.18(a) : $\mathcal A_{\mathbb R}$ est formée des parties réelles des éléments de $\mathcal A$), et $D_U$ est le
commutant d'une représentation réelle irréductible de $G$ : c'est $\mathbb R$, $\mathbb C$ ou $\mathbb H$. D'après la
démonstration de (b), le complexifié $U+iU=e\mathbb C^n=V_c(\mathbb C^{k_c}\otimes F)$ (plus son conjugué dans le cas
d'une paire) porte $\pi_c$, $\pi_c\oplus\pi_{\sigma(c)}$ ou $\pi_c\oplus\pi_c$ dans les trois cas : il est
irréductible, somme de deux représentations irréductibles non équivalentes, ou somme de deux représentations
équivalentes.

*Calcul :* le théorème 4.22 sur les modules réels irréductibles $\mathbb R^3$, $\mathbb R^2$, $\mathbb R^4$ des trois exemples — [tutoriel, section 4.7](../../tutorial/fr/04_symmetries.md#s-4-7).

<a id="r-4-23"></a>

**Proposition 4.23 (éléments réels du commutant).** Soit $X\in\mathcal A'$ réelle, de blocs $x_c$. Alors
$x_{\sigma(c)}=M_c\bar x_cM_c^*$ pour tout $c$. Si de plus $X=X^T$ :
(a) pour $c$ de type réel, $y_c=W_M^*x_cW_M$ ($W_M$ comme dans le théorème 4.21(a)) est une matrice réelle symétrique,
et $XO_c=O_c(1_{k_c}\otimes y_c)$ avec les colonnes réelles orthonormées $O_c$ du théorème 4.21(a) ;
(b) pour $c$ de type quaternionique, toute valeur propre de $x_c$ est de multiplicité paire, de sorte que toute valeur
propre de $X$ sur l'image de $P_c$ a une multiplicité divisible par $2k_c$ (**dégénérescence de Kramers**) ;
(c) pour $c$ de type complexe, $x_c$ et $x_{\sigma(c)}$ ont les mêmes valeurs propres avec les mêmes multiplicités.
Sans l'hypothèse $X=X^T$, $y_c$ dans (a) est encore réelle.

*Démonstration.* $X=\bar X$, donc le théorème 4.20(d) donne $x_{\sigma(c)}=M_c\bar x_cM_c^*$. Si $X=X^T$, alors
$X^*=X$ et chaque $x_c$ est hermitienne. (a) Avec $M_c=W_MW_M^T$, $M_c^*=\bar W_MW_M^*$ et $W_M^T\bar W_M=1$ :
$\bar y_c=W_M^T\bar x_c\bar W_M=W_M^TM_c^*x_cM_c\bar W_M=W_M^*x_cW_M=y_c$, et $y_c$ est hermitienne, donc réelle
symétrique. De plus $XO_c=XV_c(W_J\otimes W_M)=V_c(W_J\otimes x_cW_M)=O_c(1\otimes y_c)$. (b) Soit
$\theta v=M_c\bar v$, antiunitaire avec $\theta^2=-1$ (théorème 4.20(b)). De $x_cM_c=M_c\bar x_c$ on tire
$x_c\theta v=M_c\bar x_c\bar v=\theta(x_cv)$, donc, pour une valeur propre $\lambda\in\mathbb R$ de $x_c$, le
sous-espace propre $E_\lambda$ vérifie $\theta E_\lambda\subseteq E_\lambda$
($x_c\theta v=\theta(\lambda v)=\lambda\theta v$), et $\dim E_\lambda$ est paire d'après le lemme 4.19(b). Sur l'image
de $P_c$, $X$ agit comme $V_c(1_{k_c}\otimes x_c)V_c^*$, dont les valeurs propres sont celles de $x_c$, chacune
répétée $k_c$ fois. (c) $\bar x_c$ a les valeurs propres $\bar\lambda=\lambda$ de $x_c$, avec les vecteurs propres
conjugués, et $x_{\sigma(c)}$ est unitairement équivalente à $\bar x_c$. $\square$

Dans le cas quaternionique, $y_c=W_M^*x_cW_M$, avec $M_c=W_M\Omega_{m_c}W_M^T$ comme dans le lemme 4.19(c), vérifie
$y_c\Omega=\Omega\bar y_c$ (le calcul du théorème 4.21(b) avec $M_c$ à la place de $J_c$), donc $y_c=\mathbb M(h)$
pour une matrice quaternionique $h$ ; si $X=X^T$, alors $h^*=h$, puisque $\mathbb M(h^*)=y_c^*=y_c$ et que $\mathbb M$
est injective.

*Calcul :* la proposition 4.23 pour un élément symétrique réel aléatoire du commutant dans chacun des trois exemples — [tutoriel, section 4.7](../../tutorial/fr/04_symmetries.md#s-4-7).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Si un hamiltonien commute avec
> une symétrie et est réel dans une base où la symétrie est réelle — l'effet d'une symétrie par renversement du temps,
> qui y agit comme la conjugaison complexe —, la proposition 4.23 dit à quoi ressemblent ses blocs : réels symétriques
> dans les secteurs de type réel, hermitiens complexes avec un partenaire de même spectre dans les paires de type
> complexe, et hermitiens quaternioniques avec des niveaux doublement dégénérés dans les secteurs de type
> quaternionique, comme pour le renversement du temps d'un nombre impair de spins $\frac12$, dont le carré vaut $-1$
> (théorème de Kramers). La *voie triple* de Dyson (*threefold way* ; F. J. Dyson, *J. Math. Phys.* 3 (1962) 1199)
> modélise les blocs d'un hamiltonien générique de chaque sorte par des matrices aléatoires de la classe
> correspondante — le GOE et le GUE de la section 5.3, et l'ensemble symplectique gaussien (GSE), que ce cours ne
> définit pas — et prédit en conséquence la statistique des niveaux de chaque bloc.

<a id="s-4-8"></a>

## 4.8 Points de comparaison

*Dans la bibliothèque :* la comparaison avec d'autres logiciels — [tutoriel, section 4.8](../../tutorial/fr/04_symmetries.md#s-4-8).

<a id="x-4"></a>

## Exercices

1. Montrer que, pour la représentation de $S_L$ sur $(\mathbb C^d)^{\otimes L}$,
   $\chi_\pi(\sigma)=d^{\#\text{cycles}(\sigma)}$ et $\dim\pi(S_L)'=\binom{d^2+L-1}L$. Vérifier $20$ (section 2.4) et
   $35$ ($=\dim\{S_a\}''$ pour quatre spins).
2. (*Clebsch–Gordan.*) Montrer que $j_1\otimes j_2\cong\bigoplus_jj$, $j=|j_1-j_2|,|j_1-j_2|+1,\dots,j_1+j_2$, pour
   la représentation de spin $S^{(j_1)}_a\otimes1+1\otimes S^{(j_2)}_a$ — la règle de fusion non abélienne de
   $SU(2)$ que mentionne le chapitre 7.
3. Soit $\theta/2\pi$ irrationnel. Montrer que $\{U_\theta\}'=U(1)'$ pour un seul $U_\theta=e^{i\theta Q}$, mais que
   $\{e^{i\pi\hat N}\}'\cong M_{2^{n-1}}\oplus M_{2^{n-1}}$ pour $n$ qubits.
4. (*Effet Zeeman.*) Soit $H$ un hamiltonien de Heisenberg sur $L\ge1$ spins et $H_h=H+hS_z$. Montrer que $H_h$
   commute avec $S_z$ et $\mathbf S^2$ mais, pour $h\neq0$, pas avec $S_x$, et que ses valeurs propres sont
   $E+hM$, $M=-S,\dots,S$, pour les valeurs propres $E$ de $H_S$.
5. (*Le groupe des quaternions.*) Soit $Q_8=\{\pm1,\pm I,\pm J,\pm K\}\subseteq\mathbb H$ avec $I=(i,0)$, $J=(0,1)$,
   $K=(0,i)$ (section E.6), et soit $\pi$ sa représentation régulière sur $\mathbb R^8$, $\pi(g)e_h=e_{gh}$. Montrer
   que $\mathcal A=\mathrm{span}\,\pi(Q_8)$ a quatre secteurs de type réel avec $k=m=1$ et un de type
   quaternionique avec $k=m=2$, de sorte que $\mathcal A_{\mathbb R}\cong\mathbb R^4\oplus\mathbb H$.

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

3. $U_\theta=\sum_qe^{iq\theta}P_q$, et $e^{iq\theta}=e^{iq'\theta}$ si et seulement si $(q-q')\theta\in2\pi\mathbb Z$, ce
   qui, pour $\theta/2\pi$ irrationnel, impose $q=q'$. Les sous-espaces propres de $U_\theta$ sont donc les images des
   $P_q$, et un opérateur commute avec la matrice normale $U_\theta$ si et seulement s'il préserve ses sous-espaces
   propres (théorème 1.5(d) et interpolation de Lagrange comme à la section 2.2), c'est-à-dire si et seulement si
   $X=\sum_qP_qXP_q$ ; c'est $U(1)'$ d'après la proposition 4.10. L'argument utilise seulement que les $e^{iq\theta}$
   pour les charges $q$ qui apparaissent sont distincts.
   Pour $\theta=\pi$, la valeur propre $(-1)^w$ ne voit que la parité du poids de Hamming : deux sous-espaces propres de
   dimension $2^{n-1}$ (pour $n\ge1$, les configurations de poids pair et de poids impair ; inverser le premier bit est
   une bijection entre elles), et le commutant est $M_{2^{n-1}}\oplus M_{2^{n-1}}$.

4. $S_z$ et $\mathbf S^2$ commutent avec $H$ (proposition 4.15(a)) et avec $S_z$, donc avec $H_h$ ; et
   $[H_h,S_x]=h[S_z,S_x]=ihS_y\neq0$. Dans la décomposition de la proposition 4.13,
   $V^*S_zV=\bigoplus_SS^{(S)}_z\otimes1_{m_S}$ et $V^*HV=\bigoplus_S1_{2S+1}\otimes H_S$ (proposition 4.15), donc sur le
   terme de spin $S$, $V^*H_hV=hS^{(S)}_z\otimes1+1\otimes H_S$, une somme d'opérateurs qui commutent, diagonalisée
   par les produits de vecteurs propres, de valeurs propres $hM+E$.

5. D'après la démonstration du théorème E.19, $I^2=J^2=K^2=-1$ et $IJ=-JI=K$, $JK=-KJ=I$, $KI=-IK=J$, donc $Q_8$ est
   un groupe, et $h\mapsto gh$ est une bijection de $Q_8$ : les $\pi(g)$ sont des matrices de permutation, une
   représentation unitaire par des matrices réelles, et $\mathcal A=\pi(Q_8)''$ est stable par conjugaison (lemme
   4.18(c)). $\chi_\pi(g)$ est le nombre de $h$ tels que $gh=h$, c'est-à-dire $8$ pour $g=1$ et $0$ sinon, donc,
   d'après le corollaire 4.5(b), une représentation unitaire irréductible $\rho$ a la multiplicité
   $m_\rho=\frac18\overline{\chi_\rho(1)}\cdot8=\dim\rho$. *Dimension un :* aux signes près, $I,J,K$ se multiplient
   comme les éléments non neutres de $\mathbb Z_2\times\mathbb Z_2$, donc, pour $s,t\in\{\pm1\}$, $\chi(\pm1)=1$,
   $\chi(\pm I)=s$, $\chi(\pm J)=t$, $\chi(\pm K)=st$ est un morphisme ; ce sont quatre caractères distincts, donc
   quatre représentations irréductibles non équivalentes, chacune avec $m=1$. *Dimension deux :* $\rho(g)=M(g)$ (lemme
   E.15(2)) est une représentation unitaire ($M(g)^*M(g)=M(\bar gg)=1$) avec $\chi_\rho(\pm1)=\pm2$ et $\chi_\rho=0$
   ailleurs ($\mathrm{Tr}\,M(a,b)=2\mathrm{Re}\,a$), donc
   $\langle\chi_\rho,\chi_\rho\rangle_{Q_8}=\frac{4+4}8=1$ et $\rho$ est irréductible (corollaire 4.5(c)), avec $m=2$.
   Pour une décomposition $V$ comme dans le théorème 2.6, les $\pi_c\circ\pi$ sont irréductibles et deux à deux non
   équivalentes (proposition 2.7), donc $V$ est une décomposition comme dans le théorème 4.3, et chacune des cinq
   représentations ci-dessus, ayant $m_\rho\ge1$, est équivalente à exactement une $\pi_c\circ\pi$, des
   représentations différentes correspondant à des $c$ différents (corollaire 4.5(b)). Comme
   $4\cdot1\cdot1+2\cdot2=8=\sum_ck_cm_c$, il n'y a pas d'autres secteurs. D'après le théorème 4.20(c), on peut
   utiliser une décomposition avec $\pi_c(\pi(g))=\rho_c(g)$ pour ces représentations (le lemme 4.2(c) donne des
   matrices unitaires $u$ avec $\rho_c=u\,\pi_c(\pi(\cdot))\,u^*$ ; remplacer $V_c$ par $V_c(u^*\otimes1)$). Comme
   $\pi(g)$ est réelle, le théorème 4.20(a) donne $\pi_{\sigma(c)}(\pi(g))=J_c\overline{\rho_c(g)}J_c^*$. Pour un
   caractère $\chi$, cela vaut $\chi(g)$, donc $\pi_{\sigma(c)}\circ\pi\cong\pi_c\circ\pi$, $\sigma(c)=c$ (proposition
   2.7 : des secteurs différents portent des représentations non équivalentes), et $J_c\bar J_c=|J_c|^2=1$ : type
   réel. Pour $\rho=M$, $\sigma(c)=c$ (le seul secteur avec $k=2$), et $\overline{M(g)}=\Omega_2^*M(g)\Omega_2$
   (d'après la démonstration du théorème 4.21(b), $M(q)\Omega_2=\Omega_2\overline{M(q)}$), donc
   $M(g)=(J_c\Omega_2^*)M(g)(J_c\Omega_2^*)^*$. Les matrices $M(1)=1$, $M(I)=\mathrm{diag}(i,-i)$,
   $M(J)=\Omega_2$, $M(K)=\begin{pmatrix}0&-i\\-i&0\end{pmatrix}$ sont linéairement indépendantes, donc elles
   engendrent $M_2(\mathbb C)$ ; par conséquent $J_c\Omega_2^*=\lambda1$, $J_c=\lambda\Omega_2$ et
   $J_c\bar J_c=|\lambda|^2\Omega_2^2=-1$ : type quaternionique. D'après le théorème 4.21,
   $\mathcal A_{\mathbb R}\cong\mathbb R^4\oplus M_1(\mathbb H)$, de dimension $8=|Q_8|$.

</details>

*Calcul :* les exercices 1 à 5, numériquement — [tutoriel, Exercices](../../tutorial/fr/04_symmetries.md#x-4).

**Suite :** le [chapitre 5](05_random_matrices.md) tire des matrices aléatoires dans une algèbre : des unitaires de
Haar de $\prod_cU(k_c)$, le calcul de Weingarten — dont le cas d'ordre deux repose sur le commutant de
$\{V\otimes V\}$, le plus petit exemple de dualité de Schur–Weyl — et les statistiques de niveaux secteur par
secteur.
