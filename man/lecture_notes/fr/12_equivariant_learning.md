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
unitaires agissent bloc par bloc sur le commutant ; et comment on compte les vecteurs invariants, les formes
quadratiques invariantes et les polynômes invariants de tout degré (formule de Molien).

Les couches sont construites et entraînées avec la bibliothèque dans le [Tutoriel 12](../../tutorial/fr/12_equivariant_learning.md), qui vérifie aussi numériquement les résultats de ce chapitre.

Dans tout ce chapitre, $G$ est un groupe fini, et une **représentation** de $G$ sur $\mathbb C^n$ est un
homomorphisme $\pi:G\to M_n(\mathbb C)$ à valeurs unitaires (ce n'est pas une restriction : des $\pi(g)$ inversibles
sont unitaires pour le produit scalaire $\sum_g\langle\pi(g)x,\pi(g)y\rangle$, section 4.1). Son **caractère** est
$\chi_\pi(g)=\mathrm{Tr}\,\pi(g)$. Les représentations irréductibles et équivalentes et la **multiplicité**
$m_\rho(\pi)$ d'une $\rho$ irréductible dans $\pi$ ($0$ si $\rho$ n'apparaît pas) sont comme à la section 2.4, qui
démontre que $m_\rho(\pi)$ ne dépend pas de la décomposition ainsi que le lemme de Schur ; le
[chapitre 4](04_symmetries.md) traite les deux en détail (sections 4.1–4.2). Sur les matrices rectangulaires,
$\langle S,T\rangle=\mathrm{Tr}(S^*T)$ est le produit scalaire de **Hilbert–Schmidt**.

<a id="s-12-1"></a>

## 12.1 Applications linéaires équivariantes

**Définition.** Soient $\pi_1,\pi_2$ des représentations de $G$ sur $\mathbb C^{n_1}$, $\mathbb C^{n_2}$. Une matrice
$T\in M_{n_2\times n_1}(\mathbb C)$ est **équivariante** (un **opérateur d'entrelacement**) si $T\pi_1(g)=\pi_2(g)T$
pour tout $g\in G$. Ces matrices forment un espace vectoriel $\mathrm{Hom}_G(\pi_1,\pi_2)$ ; pour
$\pi_1=\pi_2=\pi$, c'est le commutant, $\mathrm{Hom}_G(\pi,\pi)=\pi(G)'$.

Les représentations les plus courantes en apprentissage automatique permutent les coordonnées. Supposons que $G$ agisse
sur un ensemble fini $\Omega$ ; la **représentation de permutation** sur $\mathbb C^\Omega$ est $\pi(g)e_x=e_{gx}$.
Exemples : $S_n$ permutant les $n$ points d'un ensemble ($\Omega=[n]=\{1,\dots,n\}$), et $S_n$ renumérotant les
sommets d'un graphe, qui agit sur $\Omega=[n]^2$, c'est-à-dire sur la matrice d'adjacence par
$X\mapsto P_\sigma XP_\sigma^T$.

<a id="r-12-1"></a>

**Lemme 12.1 (orbites).** Soient $\pi_i$ les représentations de permutation de $G$ sur $\mathbb C^{\Omega_i}$. Alors
$T\in\mathrm{Hom}_G(\pi_1,\pi_2)$ si et seulement si $T_{gy,gx}=T_{yx}$ pour tous $g$, $y\in\Omega_2$,
$x\in\Omega_1$. Les matrices indicatrices des orbites de $G$ sur $\Omega_2\times\Omega_1$ forment une base, et
$\dim\mathrm{Hom}_G(\pi_1,\pi_2)$ est le nombre de ces orbites.

*Démonstration.* $\pi_2(g)T\pi_1(g)^*e_x=\pi_2(g)Te_{g^{-1}x}=\sum_yT_{y,g^{-1}x}e_{gy}$, donc le coefficient $(gy,x)$
de $\pi_2(g)T\pi_1(g)^*$ est $T_{y,g^{-1}x}$. Par conséquent, $T$ est équivariante, c'est-à-dire
$\pi_2(g)T\pi_1(g)^*=T$ pour tout $g$, si et seulement si $T_{gy,x}=T_{y,g^{-1}x}$ pour tous $g,y,x$, ce qui (en
remplaçant $x$ par $gx$) est la condition de l'énoncé : $T$ est constante sur les orbites. Les matrices indicatrices
des orbites ont des supports disjoints, donc elles sont linéairement indépendantes, et toute fonction constante sur
les orbites en est une combinaison. $\square$

*Calcul :* le lemme 12.1 pour $S_4$ sur les caractéristiques de sommets et d'arêtes, en comptant les orbites — [tutoriel, section 12.1](../../tutorial/fr/12_equivariant_learning.md#s-12-1).

<a id="r-12-2"></a>

**Corollaire 12.2.** (a) (*DeepSets.*) Pour $S_n$ sur $\mathbb C^n$, $n\ge2$, le commutant est
$\mathrm{span}\{1,J\}$, où $J$ est la matrice dont tous les coefficients valent un, et
$\mathbb C^n=\mathbb C\mathbf 1\oplus\mathbf 1^\perp$ est la somme de deux représentations irréductibles non
équivalentes, la triviale et la représentation **standard**, de dimension $n-1$. (b) (*Réseaux de graphes.*) Les
applications équivariantes de $\mathbb C^{[n]^p}$ dans $\mathbb C^{[n]^q}$ (tenseurs d'ordre $p$ et $q$ sur $n$
sommets) sont comptées par les partitions d'un ensemble à $p+q$ éléments en au plus $n$ blocs ; pour $n\ge p+q$, leur
nombre est le **nombre de Bell** $B_{p+q}$, le nombre de toutes les partitions d'un ensemble à $p+q$ éléments
($B_2=2$, $B_3=5$, $B_4=15$).

*Démonstration.* (b) Deux uplets $x,y\in[n]^{p+q}$ sont dans la même orbite si et seulement s'ils ont le même
**motif d'égalités**, la partition des positions $\{1,\dots,p+q\}$ définie par $i\sim j\Leftrightarrow x_i=x_j$. Le
motif est clairement invariant. Réciproquement, si $x,y$ ont le même motif, $x_i\mapsto y_i$ est une injection bien
définie de $\{x_i\}$ dans $[n]$, qui se prolonge en une permutation $\sigma$ telle que $\sigma x=y$. Toute partition
en au plus $n$ blocs est un motif (donner aux blocs des valeurs distinctes). Le lemme 12.1 permet de conclure.
Pour $n\ge p+q$, la borne sur le nombre de blocs est vide. En comptant selon le nombre de blocs : $B_2=1+1$ ;
$B_3=1+3+1$ (une partition d'un ensemble à $3$ éléments en deux blocs est déterminée par son singleton) ;
$B_4=1+7+6+1$ (deux blocs : $4$ de tailles $3+1$ et $3$ de tailles $2+2$ ; trois blocs : les $\binom42=6$ choix de
la paire). (a) Pour $p=q=1$, il y a deux orbites sur $[n]^2$, la diagonale et son complémentaire, d'indicatrices $1$ et $J-1$.
Donc $\dim\pi(S_n)'=2=\sum_cm_c^2$ (théorème 2.6, proposition 2.7), ce qui impose deux secteurs avec $m_c=1$,
c'est-à-dire deux facteurs irréductibles non équivalents, et $\pi(S_n)'$ est engendré par leurs projections
$P_1,P_2$. La projection sur la droite invariante $\mathbb C\mathbf 1$ appartient à $\pi(S_n)'$ (lemme 2.1), donc c'est
$P_1$ ou $P_2$ : la droite, qui porte la représentation triviale, est un secteur, et $\mathbf 1^\perp$ est
l'autre. $\square$

La couche DeepSets $x\mapsto\lambda x+\gamma\,(\mathbf 1^Tx)\mathbf 1$ (Zaheer et al. 2017) est (a) ; les bases de
Maron et al. (2019) pour les réseaux de graphes sont (b).

*Calcul :* le corollaire 12.2 pour $n=4$ ($B_2=2$, $B_3=5$, $B_4=15$), et la décomposition des caractéristiques de sommets et d'arêtes de $S_4$ par `dc.from_group` (caractéristiques d'arêtes : secteurs de dimensions $k=(1,2,3,3)$ avec les multiplicités $m=(2,1,1,3)$) — [tutoriel, section 12.1](../../tutorial/fr/12_equivariant_learning.md#s-12-1).

Pour deux représentations, tout se ramène à un seul commutant :

<a id="r-12-3"></a>

**Lemme 12.3 (somme directe).** Soit $\pi=\pi_1\oplus\pi_2$ sur $\mathbb C^{n_1}\oplus\mathbb C^{n_2}$, avec les
inclusions isométriques $Q_i$ et les projections $P_i=Q_iQ_i^*$. Alors $P_1,P_2\in\pi(G)'$,
$\mathrm{Hom}_G(\pi_1,\pi_2)=Q_2^*\,\pi(G)'\,Q_1$, et $T\mapsto Q_2TQ_1^*$ est une bijection linéaire de
$\mathrm{Hom}_G(\pi_1,\pi_2)$ sur $P_2\,\pi(G)'P_1$.

*Démonstration.* $\pi(g)Q_i=Q_i\pi_i(g)$ et $Q_i^*\pi(g)=\pi_i(g)Q_i^*$, car $\pi(g)$ est diagonale par blocs ; donc
$P_i\pi(g)=Q_i\pi_i(g)Q_i^*=\pi(g)P_i$. Si $X\in\pi(G)'$, alors
$Q_2^*XQ_1\pi_1(g)=Q_2^*X\pi(g)Q_1=Q_2^*\pi(g)XQ_1=\pi_2(g)Q_2^*XQ_1$. Si $T$ est équivariante, alors
$\pi(g)Q_2TQ_1^*=Q_2\pi_2(g)TQ_1^*=Q_2T\pi_1(g)Q_1^*=Q_2TQ_1^*\pi(g)$, donc $Q_2TQ_1^*\in\pi(G)'$, elle est égale à
$P_2(Q_2TQ_1^*)P_1$, et $Q_2^*(Q_2TQ_1^*)Q_1=T$. Réciproquement, tout $Y=P_2XP_1$ avec $X\in\pi(G)'$ appartient à
$\pi(G)'$ et $Y=Q_2(Q_2^*YQ_1)Q_1^*$. $\square$

*Calcul :* les coins $Q_2^*\,\pi(G)'\,Q_1$ pour les caractéristiques de sommets et d'arêtes — [tutoriel, section 12.2](../../tutorial/fr/12_equivariant_learning.md#s-12-2).

<a id="r-12-4"></a>

**Théorème 12.4 (les applications équivariantes sont des matrices par blocs).** Soient $\rho_1,\dots,\rho_C$ les
représentations irréductibles deux à deux non équivalentes qui apparaissent dans $\pi_1\oplus\pi_2$,
$k_c=\dim\rho_c$, $a_c=m_{\rho_c}(\pi_1)$, $b_c=m_{\rho_c}(\pi_2)$. Il existe des unitaires
$V_1:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{a_c}\to\mathbb C^{n_1}$ et
$V_2:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{b_c}\to\mathbb C^{n_2}$ avec
$V_i^*\pi_i(g)V_i=\bigoplus_c\rho_c(g)\otimes1$, tels que
$$
\mathrm{Hom}_G(\pi_1,\pi_2)=\Big\{V_2\Big(\bigoplus_{c=1}^C1_{k_c}\otimes y_c\Big)V_1^*:\ y_c\in M_{b_c\times a_c}(\mathbb C)\Big\},
$$
et $(y_c)_c\mapsto T$ est injective. Par conséquent $\mathrm{Hom}_G(\pi_1,\pi_2)\cong\bigoplus_cM_{b_c\times a_c}(\mathbb C)$,
$$
\dim\mathrm{Hom}_G(\pi_1,\pi_2)=\sum_\rho m_\rho(\pi_1)\,m_\rho(\pi_2),
$$
la somme portant sur toutes les représentations irréductibles à équivalence près ; pour $\pi_1=\pi_2$, on retrouve
$\pi(G)'\cong\bigoplus_cM_{m_c}$ (proposition 2.7), et $\mathrm{Hom}_G(\pi_1,\pi_2)=0$ si et seulement si aucune
représentation irréductible n'apparaît dans les deux.

*Démonstration.* Appliquons le théorème 2.6 et la proposition 2.7 à $\pi=\pi_1\oplus\pi_2$ : il existe un unitaire
$V$ tel que $V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$ et $V^*\pi(G)'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$.
D'après l'unicité des multiplicités (section 2.4), toute représentation irréductible apparaissant dans $\pi$ est
équivalente à exactement une $\rho_c$, de sorte que les $\rho_c$ sont des représentants des classes de l'énoncé.
D'autres représentants quelconques sont $\rho_c'=w_c\rho_cw_c^*$ avec des unitaires $w_c$ (lemme 4.2(c)), et
remplacer $V$ par $V\bigoplus_c(w_c^*\otimes1)$ transforme la première forme normale en
$\bigoplus_c\rho'_c(g)\otimes1_{m_c}$ et laisse la seconde inchangée, car $w_c^*\otimes1$ commute avec
$1\otimes M_{m_c}$ ; on peut donc utiliser les $\rho_c$ de la proposition 2.7. Comme
$X\mapsto V^*XV$ est un $*$-isomorphisme et que $P_1\in\pi(G)'$ (lemme 12.3), $V^*P_1V=\bigoplus_c1\otimes p_c$ avec
des projections $p_c$. Choisissons des unitaires $u_c$ tels que $u_c^*p_cu_c=1_{a_c}\oplus0_{b_c}$ (théorème 1.3), où
désormais $a_c=\mathrm{rank}\,p_c$, $b_c=m_c-a_c$, et remplaçons $V$ par $V\bigoplus_c(1\otimes u_c)$ ; cela ne
modifie aucune des deux formes normales, car $1\otimes u_c$ commute avec $\rho_c(g)\otimes1$. Soient
$\iota_c:\mathbb C^{a_c}\to\mathbb C^{m_c}$ et $\kappa_c:\mathbb C^{b_c}\to\mathbb C^{m_c}$ les inclusions des $a_c$
premières et des $b_c$ dernières coordonnées, et
$$
V_1=Q_1^*V\bigoplus_c(1\otimes\iota_c),\qquad V_2=Q_2^*V\bigoplus_c(1\otimes\kappa_c).
$$
$W_1=V\bigoplus_c(1\otimes\iota_c)$ est une isométrie d'image $\mathrm{range}\,P_1=\mathrm{range}\,Q_1$, donc
$V_1^*V_1=W_1^*P_1W_1=1$, et $V_1$ est unitaire puisque $\sum_ck_ca_c=\mathrm{rank}\,P_1=n_1$. De plus
$\pi_1(g)V_1=Q_1^*\pi(g)W_1=Q_1^*V\big(\bigoplus_c\rho_c(g)\otimes1_{m_c}\big)\bigoplus_c(1\otimes\iota_c)
=V_1\bigoplus_c\rho_c(g)\otimes1_{a_c}$. Ainsi $\pi_1\cong\bigoplus_c\rho_c\otimes1_{a_c}$ avec des $\rho_c$
irréductibles non équivalentes, et $a_c=m_{\rho_c}(\pi_1)$ par l'unicité des multiplicités (section 2.4) ; de même
pour $V_2$ et $b_c$. Enfin, pour $X=V\big(\bigoplus_c1\otimes x_c\big)V^*\in\pi(G)'$,
$P_2XP_1=V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*$, où $y_c=\kappa_c^*x_c\iota_c$ est le coin inférieur
gauche de taille $b_c\times a_c$ de $x_c$, une matrice arbitraire. D'après le lemme 12.3,
$\mathrm{Hom}_G(\pi_1,\pi_2)=Q_2^*P_2\,\pi(G)'P_1Q_1$ est formé des matrices
$Q_2^*V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*Q_1=V_2\big(\bigoplus_c1\otimes y_c\big)V_1^*$, et
$(y_c)\mapsto T$ est injective puisque $V_1,V_2$ sont unitaires. Une $\rho$ irréductible qui n'apparaît pas dans
$\pi_1\oplus\pi_2$ vérifie $m_\rho(\pi_1)m_\rho(\pi_2)=0$. $\square$

*Dans la bibliothèque :* la démonstration du théorème 12.4 en tant qu'algorithme, les couches `EquivariantLinear` et `EquivariantMap` construites sur `dc.from_group`, `D.commutant_algebra()` et `D.from_commutant` ; des caractéristiques de sommets vers les caractéristiques d'arêtes de $S_4$, $\dim\mathrm{Hom}_G=1\cdot2+1\cdot3=5=B_3$ — [tutoriel, section 12.2](../../tutorial/fr/12_equivariant_learning.md#s-12-2).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Le théorème 12.4 est le noyau
> d'algèbre linéaire des **règles de sélection** et le cas scalaire du théorème de Wigner–Eckart : un opérateur qui commute avec la symétrie ne peut relier que des états
> de même type de symétrie, et à l'intérieur d'un type il agit comme $1_{k_c}\otimes y_c$, c'est-à-dire par un
> « élément de matrice réduit » $y_c$ aveugle à la position à l'intérieur du multiplet. (Les opérateurs qui se
> transforment eux-mêmes selon une représentation non triviale, comme les opérateurs vectoriels, requièrent en outre
> les coefficients de Clebsch–Gordan.) Un hamiltonien symétrique est le cas $\pi_1=\pi_2$ ; ses dégénérescences $k_c$
> sont les dimensions des représentations irréductibles.

<a id="s-12-2"></a>

## 12.2 Une couche équivariante avec la bibliothèque

*Dans la bibliothèque :* des couches équivariantes paramétrées par les blocs du théorème 12.4 — [tutoriel, section 12.2](../../tutorial/fr/12_equivariant_learning.md#s-12-2).

La seconde façon de rendre une couche équivariante est de la moyenner sur le groupe. Le **twirl** est
$$
\mathcal P(W)=\frac1{|G|}\sum_{g\in G}\pi_2(g)\,W\,\pi_1(g)^*,\qquad W\in M_{n_2\times n_1}(\mathbb C).
$$

<a id="r-12-5"></a>

**Théorème 12.5 (le twirl est la projection orthogonale).** $\mathcal P$ est la projection orthogonale, pour le
produit scalaire de Hilbert–Schmidt, de $M_{n_2\times n_1}(\mathbb C)$ sur $\mathrm{Hom}_G(\pi_1,\pi_2)$, et
$$
\dim\mathrm{Hom}_G(\pi_1,\pi_2)=\frac1{|G|}\sum_{g\in G}\overline{\chi_{\pi_1}(g)}\,\chi_{\pi_2}(g).
$$

*Démonstration.* Pour $h\in G$, $\pi_2(h)\mathcal P(W)\pi_1(h)^*=\frac1{|G|}\sum_g\pi_2(hg)W\pi_1(hg)^*=\mathcal P(W)$,
car $g\mapsto hg$ est une bijection de $G$ ; donc $\mathcal P(W)$ est équivariante. Si $W$ est équivariante, chaque
terme est égal à $W$, donc $\mathcal P(W)=W$. Par conséquent $\mathcal P^2=\mathcal P$, d'image
$\mathrm{Hom}_G(\pi_1,\pi_2)$. Elle est autoadjointe :
$\langle S,\mathcal P(W)\rangle=\frac1{|G|}\sum_g\mathrm{Tr}\big((\pi_2(g)^*S\pi_1(g))^*W\big)$
$=\langle\mathcal P(S),W\rangle$, où l'on a substitué $g\mapsto g^{-1}$ et utilisé $\pi_i(g^{-1})=\pi_i(g)^*$. Un
idempotent autoadjoint est la projection orthogonale sur son image. Sa trace est son rang. L'application linéaire
$W\mapsto AWB$ sur $M_{n_2\times n_1}$ a pour trace $\mathrm{Tr}\,A\mathrm{Tr}\,B$ (le coefficient de l'unité
matricielle $E_{ij}$ dans $AE_{ij}B$ est $A_{ii}B_{jj}$) ; avec $A=\pi_2(g)$, $B=\pi_1(g)^*$, cela donne la
formule. $\square$

Avec le théorème 12.4,
$\sum_\rho m_\rho(\pi_1)m_\rho(\pi_2)=|G|^{-1}\sum_g\overline{\chi_{\pi_1}(g)}\chi_{\pi_2}(g)$ ; pour $\pi_1=\pi_2$, c'est
la formule des caractères pour $\dim\pi(G)'$ de la section 4.2. Le paramétrage
$W\mapsto\mathcal P(W)$ atteint toute application équivariante, $\mathcal P(W)$ est l'application équivariante la plus
proche de $W$, et, comme $\mathcal P$ est autoadjointe, le gradient par rapport à $W$ est l'image par $\mathcal P$ du
gradient par rapport à $T=\mathcal P(W)$. Mais $n_1n_2-\dim\mathrm{Hom}_G$ de ses paramètres sont gaspillés dans le
noyau.

*Calcul :* le théorème 12.5 sur les caractéristiques d'arêtes et la formule de dimension des caractéristiques de sommets vers les caractéristiques d'arêtes — [tutoriel, section 12.2](../../tutorial/fr/12_equivariant_learning.md#s-12-2). *Dans la bibliothèque :* le twirl `dc.twirl`, la même projection calculée à partir de la décomposition, `D.conditional_expectation(W, 'commutant')`, et la formule des caractères `dc.commutant_dimension` — [tutoriel, section 12.2](../../tutorial/fr/12_equivariant_learning.md#s-12-2).

Les couches ci-dessus sont linéaires. Un réseau les alterne avec des non-linéarités, qui doivent elles aussi être
équivariantes. Pour une représentation de permutation, une fonction appliquée à chaque coordonnée est équivariante,
puisque $(\pi(g)x)_{gy}=x_y$ ; pour d'autres représentations, elle ne l'est pas, et l'on utilise des non-linéarités
qui agissent sur des normes ou sur les coefficients dans une base que $G$ permute (exercice 3).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Le twirl est l'espérance conditionnelle sur le commutant (exercice 3 du chapitre 2) :
> moyenner un couplage sur toutes les opérations de symétrie en conserve exactement la partie symétrique.
> « Symétriser » un hamiltonien ou une interaction en sommant sur le groupe de symétrie est la même projection.

<a id="s-12-3"></a>

## 12.3 Couches contraintes et spectrales

Les poids $w=(w_c)_c$ d'une couche équivariante $x\mapsto Tx$, $T\in\pi(G)'$, vivent dans la $*$-algèbre
$\bigoplus_cM_{m_c}$ (théorème 12.4), de sorte que toutes les constructions des chapitres 1 et 11 s'y appliquent.
Notons $\Phi(w)=V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*$ pour $w\in\bigoplus_cM_{m_c}$.

<a id="r-12-6"></a>

**Proposition 12.6.** (a) $\Phi$ est un $*$-isomorphisme de $\bigoplus_cM_{m_c}$ sur $\pi(G)'$, et $\Phi(w)$ est
autoadjoint, unitaire ou positif si et seulement si chaque $w_c$ l'est. (b) Pour $w$ normal et une fonction $f$ sur
$\bigcup_c\mathrm{spec}\,w_c$, $f(\Phi(w))=\Phi\big((f(w_c))_c\big)$ ; en particulier $f(X)\in\pi(G)'$ pour tout
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

*Dans la bibliothèque :* des couches équivariantes unitaires et positives avec `optimize.UnitaryParam` et `PositiveParam` sur `D.commutant_algebra()`, des couches spectrales avec `apply_function` ; l'entraînement d'une couche unitaire équivariante, et (b) pour $e^{-tL}$ avec sa dérivée en $t$ — [tutoriel, section 12.3](../../tutorial/fr/12_equivariant_learning.md#s-12-3).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Un hamiltonien symétrique
> $H\in\pi(G)'$ engendre une évolution temporelle symétrique $e^{-itH}\in\pi(G)'$ (chapitre 8), qui ne mélange
> jamais les secteurs et n'agit que sur les espaces de multiplicité. Une couche unitaire équivariante est une telle
> évolution, apprise au lieu d'être donnée.

<a id="s-12-4"></a>

## 12.4 Invariants

Un vecteur $x$ est **invariant** si $\pi(g)x=x$ pour tout $g$ ; une fonction $q$ sur $\mathbb C^n$ est invariante si
$q(\pi(g)x)=q(x)$. Les caractéristiques invariantes sont ce qu'un réseau invariant (un classifieur de graphes, par
exemple) peut renvoyer.

<a id="r-12-7"></a>

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

*Démonstration.* (a) Les vecteurs invariants sont $\mathrm{Hom}_G(\mathrm{triv},\pi)$, où $\mathrm{triv}$ est la
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
invariants, et la dimension de ceux-ci est $\mathrm{Tr}\,M=|G|^{-1}\sum_g\mathrm{Tr}\,\sigma(g)$ (un
idempotent est $1\oplus0$ dans une base adaptée à image $\oplus$ noyau, donc sa trace est son rang). $A$ est
unitaire, donc $A=R\Lambda R^{-1}$ avec $\Lambda=\mathrm{diag}(\lambda)$ (théorème 1.3). Les matrices symétriques
$B_{ij}=(R^{-1})^T(e_ie_j^T+e_je_i^T)R^{-1}$, $i\le j$, forment une base, et
$A^TB_{ij}A=(R^{-1})^T\Lambda(e_ie_j^T+e_je_i^T)\Lambda R^{-1}=\lambda_i\lambda_jB_{ij}$. Donc
$\mathrm{Tr}\,\sigma(g)=\sum_{i\le j}\lambda_i\lambda_j=\frac12\big((\mathrm{Tr}\,A)^2+\mathrm{Tr}\,A^2\big)
=\frac12\big(\chi_\pi(g^{-1})^2+\chi_\pi(g^{-2})\big)$, et sommer sur $g^{-1}$ au lieu de $g$ donne la formule. Pour une
représentation de permutation, $G$ permute les monômes $x_ax_b$, donc un polynôme est invariant si et seulement si ses
coefficients sont constants sur les orbites de monômes. $\square$

*Calcul :* (a), (b) et les deux décomptes de (c) pour les caractéristiques d'arêtes de $S_4$ — [tutoriel, section 12.4](../../tutorial/fr/12_equivariant_learning.md#s-12-4).

Pour les graphes, les caractéristiques d'arêtes invariantes sont engendrées par les indicatrices de la diagonale et de
son complémentaire ($B_2=2$), et les lectures (*readouts*) linéaires invariantes sont la somme des coefficients
diagonaux et celle des coefficients hors diagonale. Sur les caractéristiques de sommets $\mathbb C^4$, (c) compte les
orbites de $S_4$ sur les multiensembles de deux sommets : il y en a deux.

Sur $\mathbb C^4$, ce sont $\sum_ix_i^2$ et $\sum_{i<j}x_ix_j$. La démonstration de (c) fonctionne en tout degré :

<a id="r-12-8"></a>

**Proposition 12.8 (formule de Molien).** Soient $\mathbb C[x]_d$ l'espace des polynômes homogènes de degré $d$ en
$x_1,\dots,x_n$ et $\mathbb C[x]^G_d$ le sous-espace des polynômes invariants, $p(\pi(g)x)=p(x)$ pour tout $g$. Pour
$|t|<1$,
$$
\sum_{d\ge0}\dim\mathbb C[x]^G_d\,t^d=\frac1{|G|}\sum_{g\in G}\frac1{\det\big(1-t\pi(g)\big)} .
$$
Pour une représentation de permutation, $\dim\mathbb C[x]^G_d$ est le nombre d'orbites de $G$ sur les monômes de
degré $d$ (les multiensembles de $d$ éléments de $\Omega$). Le coefficient de $t^2$ est la proposition 12.7(c).

*Démonstration.* Pour une matrice inversible $M$, posons $\phi_M(p)(x)=p(Mx)$ (substitution). $\phi_M$ est linéaire,
envoie $\mathbb C[x]_d$ dans lui-même, et $\phi_N\circ\phi_M=\phi_{MN}$, donc $\phi_M$ est bijective d'inverse
$\phi_{M^{-1}}$. Par conséquent $\sigma(g)=\phi_{\pi(g)^{-1}}$ vérifie $\sigma(g)\sigma(h)=\sigma(gh)$, et $p$ est
invariant si et seulement si $\sigma(g)p=p$ pour tout $g$. Comme dans la démonstration de la proposition 12.7(c),
$M_d=|G|^{-1}\sum_g\sigma(g)$ sur $\mathbb C[x]_d$ est idempotente d'image $\mathbb C[x]^G_d$, donc
$\dim\mathbb C[x]^G_d=\mathrm{Tr}\,M_d=|G|^{-1}\sum_g\mathrm{Tr}\,\sigma_d(g)$, où $\sigma_d(g)$ est la restriction
de $\sigma(g)$ à $\mathbb C[x]_d$. Fixons $g$ et écrivons l'unitaire $A=\pi(g)^{-1}$ sous la forme $A=R\Lambda R^{-1}$,
$\Lambda=\mathrm{diag}(\lambda_1,\dots,\lambda_n)$, $|\lambda_i|=1$ (théorème 1.3). Avec les formes linéaires
$\ell_i(x)=(R^{-1}x)_i$, les produits $\ell^\alpha=\prod_i\ell_i^{\alpha_i}=\phi_{R^{-1}}(x^\alpha)$, $|\alpha|=d$, forment
une base de $\mathbb C[x]_d$ (l'image de la base des monômes par la bijection $\phi_{R^{-1}}$), et comme
$R^{-1}A=\Lambda R^{-1}$,
$$
(\sigma(g)\ell^\alpha)(x)=\prod_i\big((R^{-1}Ax)_i\big)^{\alpha_i}=\prod_i\big(\lambda_i\ell_i(x)\big)^{\alpha_i}=\lambda^\alpha\ell^\alpha(x).
$$
Donc $\mathrm{Tr}\,\sigma_d(g)=\sum_{|\alpha|=d}\lambda^\alpha$. Pour $|t|<1$, le produit des $n$ séries absolument
convergentes $(1-t\lambda_i)^{-1}=\sum_{a\ge0}t^a\lambda_i^a$ peut être développé et regroupé par degré :
$$
\frac1{\det(1-tA)}=\prod_{i=1}^n\frac1{1-t\lambda_i}=\sum_\alpha t^{|\alpha|}\lambda^\alpha=\sum_{d\ge0}t^d\mathrm{Tr}\,\sigma_d(g).
$$
La moyenne sur $g$, en remplaçant $g$ par $g^{-1}$ (de sorte que $A=\pi(g)$), donne la formule ; les deux membres
convergent, car $\dim\mathbb C[x]^G_d\le\dim\mathbb C[x]_d=\binom{n+d-1}d$. Pour une représentation de permutation,
$(\pi(g)^{-1}x)_y=x_{gy}$, donc $\sigma(g)$ envoie le monôme $\prod_jx_{y_j}$ sur le monôme $\prod_jx_{gy_j}$ : $G$
permute les monômes de degré $d$, un polynôme est invariant si et seulement si ses coefficients sont constants sur
les orbites, et les sommes sur les orbites forment une base de $\mathbb C[x]^G_d$. $\square$

(C'est le théorème de Molien de 1897 ; pour aller plus loin, voir B. Sturmfels, *Algorithms in Invariant Theory*,
Springer 2008, section 2.2.) Pour $S_4$ sur $\mathbb C^4$, les nombres $1,1,2,3,5$ pour $d\le4$ sont les partitions de
$d$ en au plus quatre parts (les orbites de multiensembles sont déterminées par leurs multiplicités, triées).

*Calcul :* la série de Molien jusqu'à $t^4$ confrontée aux orbites de monômes, pour les caractéristiques de sommets et d'arêtes — [tutoriel, section 12.4](../../tutorial/fr/12_equivariant_learning.md#s-12-4).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Dans la théorie de Landau des
> transitions de phase, l'énergie libre est un polynôme invariant en un paramètre d'ordre réel $x$, sur lequel $G$
> agit par des matrices orthogonales réelles. Sa partie quadratique est $x^TBx$ avec $B$ dans le commutant (car
> $\pi(g)^T=\pi(g)^{-1}$), donc, d'après la proposition 12.6(a), elle est diagonale par blocs selon les secteurs. La
> théorie de Landau postule alors que la phase symétrique devient instable dans le secteur où une valeur propre de
> $B$ change de signe en premier, et que, génériquement, la représentation irréductible de ce secteur détermine la
> symétrie de la phase ordonnée.

<a id="s-12-5"></a>

## 12.5 Points de comparaison

*Dans la bibliothèque :* la comparaison avec d'autres logiciels (escnn, e3nn, EMLP) et ce qu'apporte la bibliothèque — [tutoriel, section 12.5](../../tutorial/fr/12_equivariant_learning.md#s-12-5).

<a id="x-12"></a>

## Exercices

1. (**DeepSets avec canaux.**) $S_n$ ($n\ge2$) agit sur $\mathbb C^n\otimes\mathbb C^d$ par $P_\sigma\otimes1_d$.
   Montrer que les applications équivariantes de $\mathbb C^n\otimes\mathbb C^{d_1}$ dans
   $\mathbb C^n\otimes\mathbb C^{d_2}$ sont les $1\otimes A+J\otimes B$ avec $A,B\in M_{d_2\times d_1}$.
2. (**Couches réelles.**) Supposons que $\pi_1,\pi_2$ soient à matrices réelles (par exemple des représentations de
   permutation). Montrer que les matrices équivariantes réelles forment un espace vectoriel réel de dimension
   $\dim_{\mathbb C}\mathrm{Hom}_G(\pi_1,\pi_2)$, et que $T\mapsto\mathrm{Re}\,T$ envoie
   $\mathrm{Hom}_G(\pi_1,\pi_2)$ sur cet espace.
3. (**Non-linéarités.**) (a) Montrer qu'appliquer une fonction $f:\mathbb C\to\mathbb C$ à chaque coordonnée commute
   avec toute représentation de permutation. (b) Montrer que cela ne commute pas avec la rotation d'angle $\pi/2$ de
   $\mathbb R^2$ pour $f(t)=\max(t,0)$. (c) Montrer que $x\mapsto\varphi(\|x\|)x$ commute avec toute représentation
   unitaire.
4. (**Combien d'échantillons déterminent une couche ?**) Soient $T=\Phi(w)\in\pi(G)'$ et
   $x_1,\dots,x_N\in\mathbb C^n$. Écrivons $V^*x_i$ dans le secteur $c$ sous la forme
   $\sum_{a,j}(X_{i,c})_{aj}e_a\otimes e_j$ avec $X_{i,c}\in M_{k_c\times m_c}$, et soit $S_c\in M_{Nk_c\times m_c}$ la
   matrice obtenue en empilant $X_{1,c},\dots,X_{N,c}$. Montrer que les $T$ équivariantes telles que
   $Tx_1=\dots=Tx_N=0$ forment un espace de dimension $\sum_cm_c(m_c-\mathrm{rank}\,S_c)$ ; pour des vecteurs
   $x_i$ gaussiens complexes standard indépendants, cette dimension est presque sûrement
   $\sum_cm_c\max(0,m_c-Nk_c)$. En déduire qu'une $T$ équivariante est déterminée par $Tx_1,\dots,Tx_N$ (parmi les
   applications équivariantes) si et seulement si cette dimension est $0$, et que, pour la représentation des arêtes,
   deux échantillons gaussiens déterminent presque sûrement une couche équivariante et qu'un seul échantillon ne le
   fait jamais, alors qu'une matrice $16\times16$ quelconque n'est jamais déterminée par moins de seize
   échantillons.

*Calcul :* les exercices 1, 3 et 4, numériquement — [tutoriel, Exercices](../../tutorial/fr/12_equivariant_learning.md#x-12).

<details>
<summary>Solutions</summary>

1. $P_\sigma J=J=JP_\sigma$, donc ces applications sont équivariantes. Elles sont linéairement indépendantes : si
   $1\otimes A+J\otimes B=0$, les blocs hors diagonale donnent $B=0$, puis $A=0$. D'après le corollaire 12.2(a) et le
   théorème 2.6, $P_\sigma\otimes1_d\cong(\mathrm{triv}\otimes1_d)\oplus(\mathrm{std}\otimes1_d)$, donc les deux
   représentations irréductibles ont multiplicité $d$, et le théorème 12.4 donne $\dim\mathrm{Hom}_G=2d_1d_2$,
   la dimension de l'espace des $1\otimes A+J\otimes B$.

2. Si $T$ est équivariante, $\bar T$ l'est aussi, puisque $\overline{\pi_i(g)}=\pi_i(g)$. Par conséquent
   $\mathrm{Re}\,T=(T+\bar T)/2$ et $\mathrm{Im}\,T=(T-\bar T)/2i$ sont réelles et équivariantes, de sorte que
   l'application $T\mapsto\mathrm{Re}\,T$ est surjective sur les matrices équivariantes réelles (elle est
   l'identité sur celles-ci), et $\mathrm{Hom}_G(\pi_1,\pi_2)=H_{\mathbb R}\oplus iH_{\mathbb R}$, où
   $H_{\mathbb R}$ désigne les matrices équivariantes réelles. Une base réelle de $H_{\mathbb R}$ est alors une base
   complexe de $\mathrm{Hom}_G$ : elle l'engendre, et une relation complexe $\sum_j(\alpha_j+i\beta_j)T_j=0$ avec
   des $T_j$ réelles se scinde en parties réelle et imaginaire. Dans le lemme 12.1, les indicatrices des orbites
   forment une telle base réelle. (Les blocs $y_c$ du théorème 12.4 sont complexes même pour des données réelles,
   car la $V$ calculée par la bibliothèque n'est pas nécessairement réelle ; paramétrer $\mathrm{Re}\,T$ par des
   $y_c$ complexes utilise deux fois plus de paramètres réels que nécessaire.)

3. (a) $(\pi(g)f(x))_{gy}=f(x)_y=f(x_y)=f((\pi(g)x)_{gy})$. (b) Pour $R(x_1,x_2)=(-x_2,x_1)$ et $x=(1,-1)$,
   $f(Rx)=f(1,1)=(1,1)$ mais $Rf(x)=R(1,0)=(0,1)$. (c) $\|\pi(g)x\|=\|x\|$, donc
   $\varphi(\|\pi(g)x\|)\pi(g)x=\pi(g)\varphi(\|x\|)x$.

4. $(1\otimes w_c)(e_a\otimes e_j)=\sum_le_a\otimes(w_c)_{lj}e_l$, donc $V^*Tx_i$ a, dans le secteur $c$, la matrice
   de coefficients $X_{i,c}w_c^T$. Ainsi $Tx_i=0$ pour tout $i$ si et seulement si $S_cw_c^T=0$ pour tout $c$,
   c'est-à-dire si et seulement si chaque colonne de $w_c^T$ appartient à $\ker S_c$, un espace de dimension
   $m_c-\mathrm{rank}\,S_c$ : les $w_c$ ayant cette propriété forment un espace de dimension
   $m_c(m_c-\mathrm{rank}\,S_c)$, et $\Phi$ est injective. Pour des $x_i$ gaussiens, les coordonnées de $V^*x_i$
   sont celles de $x_i$ dans la base orthonormée des colonnes de $V$, de nouveau des gaussiennes complexes standard
   indépendantes (proposition A.16(a)), et indépendantes en $i$ ; donc $S_c$ est une matrice de Ginibre
   $Nk_c\times m_c$, qui a une densité (proposition A.16(b)). Son rang est au plus $\min(Nk_c,m_c)$, avec égalité sauf
   si la sous-matrice carrée supérieure gauche de cette taille est singulière ; le $|\det|^2$ de cette sous-matrice est
   un polynôme réel en les parties réelles et imaginaires des coefficients, non identiquement nul (il vaut $1$ lorsque
   la sous-matrice est l'identité), donc, d'après le corollaire A.7(b), le rang est presque sûrement
   $\min(Nk_c,m_c)$. Pour la représentation des arêtes, $(k,m)=(1,2),(2,1),(3,1),(3,3)$ (la décomposition calculée),
   seul le secteur $(1,2)$ contribue : $2\cdot(2-1)=2$ pour $N=1$, et $0$ pour $N=2$. Une matrice $16\times16$
   quelconque est le cas du groupe trivial (un secteur, $k=1$, $m=16$) : la dimension est
   $16(16-\mathrm{rank}\,S)\ge16(16-N)>0$ pour tout choix de $N<16$ échantillons, car
   $\mathrm{rank}\,S\le N$, et presque sûrement $0$ pour $N=16$ échantillons gaussiens. Comme
   $T\mapsto(Tx_1,\dots,Tx_N)$ est linéaire, deux applications équivariantes coïncident sur les échantillons si et
   seulement si leur différence appartient à l'espace que l'on vient de calculer, ce qui donne le critère. Pour $N=1$,
   l'espace est de dimension au moins $2$ pour *tout* $x_1$, car $\mathrm{rank}\,S_c\le Nk_c=1<m_c=2$ dans le
   secteur $(1,2)$.

</details>

**Pour aller plus loin.** Le dernier chapitre, le [chapitre 13](13_jordan_algebras.md), explique pourquoi la mécanique quantique utilise des *-algèbres plutôt que des algèbres de Jordan. Le [chapitre 4](04_symmetries.md) développe la théorie
des représentations utilisée ici : le lemme de Schur, les caractères et l'algorithme randomisé sur lequel repose
`dc.from_group` (sections 4.1–4.3). Le [chapitre 9](09_spin_chains.md) applique les mêmes secteurs aux chaînes de spins,
et la mesure de Haar sur $\prod_cU(m_c)$ du [chapitre 5](05_random_matrices.md) fournit l'initialisation aléatoire
naturelle d'une couche unitaire équivariante. Pour les mathématiques, voir J.-P. Serre, *Linear Representations of
Finite Groups* (Springer 1977) et W. Fulton et J. Harris, *Representation Theory* (Springer 1991) ; pour
l'apprentissage équivariant, M. Bronstein, J. Bruna, T. Cohen et P. Veličković, *Geometric Deep Learning* (2021), et
M. Weiler et al., *Equivariant and Coordinate Independent Convolutional Networks* (2023) ; pour l'algorithme de
décomposition, les références de `torch_vn_algebra/decompose.py`.
