[English](../08_dynamics.md) · [Español](../es/08_dynamics.md) · **Français** · [Deutsch](../de/08_dynamics.md) · [中文](../zh/08_dynamics.md) · [日本語](../ja/08_dynamics.md) · [Русский](../ru/08_dynamics.md)

# 8. Dynamique

En mécanique quantique, un système fermé est modélisé par un groupe à un paramètre d'unitaires ; un
système ouvert — couplé à un environnement que l'on n'observe pas — par un semi-groupe de canaux. Dans une algèbre
à secteurs, les deux présentent une particularité invisible sur une algèbre de matrices unique : un
hamiltonien de $`\mathcal M`$ ne peut jamais changer le secteur, de sorte que les charges sont
conservées, alors qu'un système ouvert peut sauter d'un secteur à l'autre, par exemple en perdant
une particule. Ce chapitre démontre les faits fondamentaux sur ces deux types d'évolution, écrit
l'équation de Lindblad avec des sauts entre secteurs et montre que « aucune particule n'a été
perdue » est une mesure de la charge dont l'état post-sélectionné évolue selon un hamiltonien non
hermitien.

**Vous apprendrez** pourquoi les poids des secteurs sont des constantes du mouvement, pourquoi les
générateurs de forme de Gorini–Kossakowski–Sudarshan–Lindblad (GKSL) produisent des canaux, comment
les sauts entre secteurs entrent dans l'équation de Lindblad et quand les poids des secteurs
obéissent à une équation de taux classique, ce que signifient l'évolution sans saut et la
post-sélection, et comment trouver les états stationnaires et les taux de relaxation. Les
vérifications numériques de ce chapitre se trouvent dans le [Tutoriel 8](../../tutorial/fr/08_dynamics.md).

Dans tout ce chapitre, $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$ agit sur
$`\mathcal H=\bigoplus_c\mathbb C^{k_c}`$, $`D=\sum_ck_c`$, les $`P_c`$ sont les projecteurs centraux
minimaux, les états sont donnés par leurs densités brutes $`\rho`$, de poids des secteurs
$`p_c=\mathrm{Tr}\,\rho_c`$ (chapitre 3), $`[A,B]=AB-BA`$, $`\{A,B\}=AB+BA`$ et $`\hbar=1`$. Un
**superopérateur** est une application linéaire sur $`\mathcal M`$ ou sur $`M_D(\mathbb C)`$ ; toutes
les normes de superopérateurs sont des normes d'opérateur, donc sous-multiplicatives. Pour une
application linéaire $`B`$ sur un espace de dimension finie, $`e^{tB}=\sum_{n\ge0}t^nB^n/n!`$.

<a id="s-8-1"></a>

## 8.1 Systèmes fermés

Un **hamiltonien** est un élément autoadjoint $`H\in\mathcal M`$. Un vecteur évolue selon l'**équation
de Schrödinger** $`i\psi'(t)=H\psi(t)`$, une densité selon l'**équation de von Neumann**
$`\rho'(t)=-i[H,\rho(t)]`$. Ce sont deux équations linéaires à coefficients constants sur un espace de
dimension finie, et le lemme suivant est tout ce dont nous avons besoin de la théorie des équations
différentielles.

<a id="r-8-1"></a>

**Lemme 8.1.** Soit $`B`$ une application linéaire sur un espace vectoriel normé de dimension finie
$`X`$. Alors $`\frac d{dt}e^{tB}=Be^{tB}=e^{tB}B`$, $`e^{(s+t)B}=e^{sB}e^{tB}`$,
$`\|e^{tB}\|\le e^{|t|\,\|B\|}`$, et pour tout $`x_0\in X`$, l'application $`x(t)=e^{tB}x_0`$ est l'unique
$`x`$ dérivable sur un intervalle $`I\ni0`$ telle que $`x'=Bx`$ et $`x(0)=x_0`$.

*Démonstration.* La série converge absolument, avec la majoration de la norme, et, en tant que série
entière en $`t`$, elle se dérive terme à terme, ce qui donne la première formule. Ensuite
$`\frac d{dt}(e^{tB}e^{-tB})=Be^{tB}e^{-tB}-e^{tB}Be^{-tB}=0`$, donc $`e^{tB}e^{-tB}=1`$. Si $`x'=Bx`$,
alors $`\frac d{dt}(e^{-tB}x(t))=e^{-tB}(x'-Bx)=0`$, donc $`e^{-tB}x(t)=x_0`$ et $`x(t)=e^{tB}x_0`$.
Enfin, $`t\mapsto e^{(s+t)B}x_0`$ et $`t\mapsto e^{tB}e^{sB}x_0`$ résolvent $`y'=By`$ avec la même valeur
initiale. $`\square`$

<a id="r-8-2"></a>

**Proposition 8.2.** Soit $`H=H^*\in\mathcal M`$ et $`U(t)=e^{-itH}`$, défini par le calcul fonctionnel
de la section 1.3 avec $`f_t(\lambda)=e^{-it\lambda}`$.

1. $`U(t)=\bigoplus_ce^{-itH_c}`$ est un unitaire de $`\mathcal M`$, $`U(s+t)=U(s)U(t)`$ et
   $`U'(t)=-iHU(t)`$ ; $`U(t)`$ est égal à la série exponentielle de $`-itH`$.
2. $`\psi(t)=U(t)\psi_0`$ est l'unique solution de $`i\psi'=H\psi`$ telle que $`\psi(0)=\psi_0`$, et
   $`\rho(t)=U(t)\rho_0U(t)^*`$ est l'unique solution dans $`\mathcal M`$ de $`\rho'=-i[H,\rho]`$ telle
   que $`\rho(0)=\rho_0`$. Si $`\rho_0`$ est un état, chaque $`\rho(t)`$ l'est aussi.
3. (*Représentation de Heisenberg*) $`\mathrm{Tr}(\rho(t)A)=\mathrm{Tr}(\rho_0A(t))`$ avec
   $`A(t)=U(t)^*AU(t)`$, et $`A'=i[H,A]`$.

*Démonstration.* (1) Diagonalisons chaque bloc, $`H_c=V_c\Lambda_cV_c^*`$ ; alors $`V=\bigoplus_cV_c`$
diagonalise $`H`$, et $`U(t)=Ve^{-it\Lambda}V^*`$ est diagonal par blocs (d'après le théorème 1.5(a), le
résultat ne dépend pas de l'unitaire diagonalisant). Comme $`|f_t|^2=1`$ et $`f_sf_t=f_{s+t}`$, le
théorème 1.5(c) donne l'unitarité et la loi de groupe. En dérivant $`e^{-it\Lambda}`$ coefficient par
coefficient, on obtient $`U'=V(-i\Lambda)e^{-it\Lambda}V^*=-iHU`$, donc $`U(t)`$ est la série $`e^{-itH}`$
d'après le lemme 8.1 ($`X=M_D`$, $`B`$ la multiplication à gauche par $`-iH`$). (2)
$`(U\psi_0)'=-iHU\psi_0`$, et $`(U\rho_0U^*)'=-iHU\rho_0U^*+U\rho_0U^*(iH)`$ car $`(U')^*=iU^*H`$.
L'unicité est le lemme 8.1 avec $`B=-iH`$ sur $`\mathcal H`$, resp. $`B=-i[H,\cdot\,]`$ sur $`\mathcal M`$.
Positivité : $`\langle x,U\rho_0U^*x\rangle=\langle U^*x,\rho_0U^*x\rangle\ge0`$ ; la trace, par
cyclicité. (3) Cyclicité de la trace, et $`(U^*AU)'=iU^*HAU-iU^*AHU=i[H,A(t)]`$. $`\square`$

*Dans la bibliothèque :* $`U(t)=e^{-itH}`$, bloc par bloc, est `dynamics.propagator` ; les solutions de la proposition 8.2(2) sont `dynamics.schrodinger` et `dynamics.von_neumann` — [tutoriel, section 8.1](../../tutorial/fr/08_dynamics.md#s-8-1).

<a id="r-8-3"></a>

**Théorème 8.3 (conservation de la charge).** Soit $`t\mapsto H(t)`$, définie sur un intervalle $`I`$, à
valeurs dans les éléments autoadjoints de $`\mathcal M`$.

1. Si $`\psi:I\to\mathcal H`$ est dérivable et vérifie $`i\psi'=H(t)\psi`$, alors chaque
   $`\|P_c\psi(t)\|`$ est constant.
2. Si $`\rho:I\to\mathcal M`$ est dérivable et vérifie $`\rho'=-i[H(t),\rho]`$, alors
   $`\mathrm{Tr}(Z\rho(t))`$ est constant pour tout $`Z`$ central ; en particulier, chaque poids
   de secteur $`p_c(t)=\mathrm{Tr}(P_c\rho(t))`$ l'est.
3. Si $`H`$ est constant et si $`A\in\mathcal M`$ commute avec $`H`$, alors $`\mathrm{Tr}(A\rho(t))`$
   est constant.

*Démonstration.* (1) Le produit scalaire est antilinéaire en le premier argument, donc
$`\frac d{dt}\langle\psi,P_c\psi\rangle=\langle-iH\psi,P_c\psi\rangle+\langle\psi,-iP_cH\psi\rangle =i\langle\psi,HP_c\psi\rangle-i\langle\psi,P_cH\psi\rangle=0`$,
puisque $`P_c`$ est central. (2)
$`\frac d{dt}\mathrm{Tr}(Z\rho)=-i\mathrm{Tr}(Z[H,\rho])=-i\mathrm{Tr}([Z,H]\rho)=0`$
par cyclicité. (3) Le même calcul avec $`A`$. $`\square`$

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Un hamiltonien
> de $`\mathcal M`$ commute avec les charges (nombre de particules, aimantation, nombre quantique de
> symétrie), si bien que la probabilité de chaque charge est une constante du mouvement, même pour des
> hamiltoniens forcés, dépendant du temps.

Sur $`\mathcal M`$, un vecteur unitaire $`\psi=\bigoplus_c\psi_c`$ définit l'état du mélange
$`\bigoplus_c|\psi_c\rangle\langle\psi_c|`$ : tout $`A\in\mathcal M`$ est diagonal par blocs, donc
$`\langle\psi,A\psi\rangle=\sum_c\langle\psi_c,A_c\psi_c\rangle=\mathrm{Tr}\big(A\bigoplus_c|\psi_c\rangle\langle\psi_c|\big)`$
(comme au chapitre 3, exercice 3). Si $`\psi(t)`$ résout $`i\psi'=H\psi`$, alors
$`\bigoplus_c|\psi_c(t)\rangle\langle\psi_c(t)|=U(t)\big(\bigoplus_c|\psi_c(0)\rangle\langle\psi_c(0)|\big)U(t)^*`$, puisque $`U(t)`$ est
diagonal par blocs.

*Calcul :* le théorème 8.3 et la concordance des évolutions de $`\psi`$ et de $`\bigoplus_c|\psi_c\rangle\langle\psi_c|`$, avec les solveurs exacts et avec le solveur de Runge–Kutta `dynamics.schrodinger_rk4` pour un hamiltonien forcé — [tutoriel, section 8.1](../../tutorial/fr/08_dynamics.md#s-8-1).

<a id="s-8-2"></a>

## 8.2 Semi-groupes dynamiques quantiques

Les applications CP et TP et les canaux sont définis à la section 6.1 ([chapitre
6](06_quantum_channels.md)). Outre le théorème 6.2 et la proposition 6.4 (une application
$`A\mapsto\sum_lK_lAK_l^*`$, y compris entre secteurs, est CP), nous avons besoin de trois propriétés
de stabilité.

<a id="r-8-4"></a>

**Lemme 8.4.** (a) Les sommes, les multiples positifs et les composées d'applications CP sont CP.
Soit $`\Phi`$ CP sur $`\mathcal A=M_D`$ ou sur $`\mathcal A=\mathcal M`$, et soit $`\mathcal N\subseteq\mathcal A`$
égale à $`\mathcal M`$ (si $`\mathcal A=M_D`$) ou à l'un des coins $`P_c\mathcal MP_c\cong M_{k_c}`$ (d'unité
$`P_c`$). Si $`\Phi`$ envoie $`\mathcal N`$ dans elle-même, alors $`\Phi|_{\mathcal N}`$ est CP. (b) Si des applications CP $`\Phi_m`$ convergent vers $`\Phi`$
(c'est-à-dire $`\Phi_m(A)\to\Phi(A)`$ pour tout $`A`$), alors $`\Phi`$ est CP. (c) Si $`s\mapsto\Psi(s)`$
est une application continue de $`[a,b]`$ dans les applications CP, alors $`\int_a^b\Psi(s)\,ds`$ est
CP.

*Démonstration.* (a) $`\Phi\mapsto\Phi_n`$ est linéaire et $`(\Phi\Psi)_n=\Phi_n\Psi_n`$ ; une somme de
matrices positives est positive. Pour la restriction : $`M_n(\mathcal N)\subseteq M_n(\mathcal A)`$,
et la positivité est la même condition dans les deux. En effet, d'après la section 6.1, un élément
$`[A_{ij}]`$ de $`M_n(\mathcal M)`$ est positif si et seulement si chaque matrice par blocs
$`[(A_{ij})_c]`$ l'est ; pour $`\mathcal A=M_D`$, c'est la positivité de $`[A_{ij}]\in M_{nD}`$, car après
une permutation de la base $`[A_{ij}]`$ est la somme directe de ces blocs (tous les $`A_{ij}`$ sont
diagonaux par blocs) ; et un élément de $`M_n(P_c\mathcal MP_c)`$ a pour unique bloc non nul
$`[(A_{ij})_c]`$, dont la positivité est la positivité dans $`M_n(M_{k_c})`$. (b) Les matrices positives
forment un fermé ($`\langle x,Xx\rangle\ge0`$ pour chaque $`x`$ est une condition fermée), et
$`(\Phi_m)_n(X)\to\Phi_n(X)`$. (c) Les sommes de Riemann sont CP d'après (a) et convergent vers
l'intégrale (lemme A.4(a)) ; on applique (b). $`\square`$

**Définition.** Un **semi-groupe dynamique quantique** sur $`\mathcal M`$ est une famille
$`(T_t)_{t\ge0}`$ de canaux $`\mathcal M\to\mathcal M`$ telle que $`T_0=\mathrm{id}`$, $`T_{s+t}=T_sT_t`$ et
$`t\mapsto T_t`$ est continue.

<a id="r-8-5"></a>

**Lemme 8.5.** Tout semi-groupe continu $`(T_t)_{t\ge0}`$ de superopérateurs ($`T_0=\mathrm{id}`$,
$`T_{s+t}=T_sT_t`$) s'écrit $`T_t=e^{t\mathcal L}`$ pour un unique superopérateur $`\mathcal L`$, son
**générateur**, et $`\mathcal L=\lim_{t\to0}(T_t-\mathrm{id})/t`$. Il en va de même, avec la même
démonstration, des semi-groupes continus d'applications linéaires sur un espace de dimension finie
quelconque (par exemple de matrices).

*Démonstration.* (Les intégrales de fonctions continues à valeurs superopérateurs sont celles du
lemme A.4.) $`\|V_h-\mathrm{id}\|\le\max_{0\le s\le h}\|T_s-\mathrm{id}\|\to0`$ pour
$`V_h=\frac1h\int_0^hT_s\,ds`$ (lemme A.4(c)), donc $`\det V_h\neq0`$ pour un certain $`h`$ petit, fixé, par
continuité du déterminant. Or $`T_tV_h=\frac1h\int_0^hT_{t+s}\,ds=\frac1h\int_t^{t+h}T_s\,ds`$ (lemme
A.4(b),(e)) est dérivable en $`t`$, de
dérivée $`\frac1h(T_{t+h}-T_t)=T_t\frac{T_h-\mathrm{id}}h`$. Par conséquent $`T_t`$ est dérivable et
$`T_t'=T_t\mathcal L`$ avec $`\mathcal L=\frac{T_h-\mathrm{id}}hV_h^{-1}`$, et $`T_t=e^{t\mathcal L}`$
d'après le lemme 8.1 (avec pour $`B`$ la composition à droite par $`\mathcal L`$). Unicité :
$`\mathcal L=T'(0)`$. $`\square`$

<a id="r-8-6"></a>

**Théorème 8.6 (générateurs GKSL).** Soient $`H=H^*`$ et $`L_1,\dots,L_J`$ dans $`M_D`$, et
```math
\mathcal L(\rho)=-i[H,\rho]+\sum_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
```
Alors $`(e^{t\mathcal L})_{t\ge0}`$ est un semi-groupe dynamique quantique sur $`M_D`$. Plus
précisément, posons $`H_{\mathrm{eff}}=H-\frac i2\sum_jL_j^*L_j`$,
$`\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)`$ et
$`\Phi(\rho)=\sum_jL_j\rho L_j^*`$, de sorte que $`\mathcal L=\mathcal L_0+\Phi`$, et définissons
```math
S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},\qquad
S_n(t)=\int_0^tS_0(t-s)\,\Phi\,S_{n-1}(s)\,ds\quad(n\ge1).
```
Alors chaque $`S_n(t)`$ est CP, et $`e^{t\mathcal L}=\sum_{n\ge0}S_n(t)`$ (le **développement de
Dyson**), avec $`\|S_n(t)\|\le e^{t\|\mathcal L_0\|}(t\|\Phi\|)^n/n!`$.

*Démonstration.* *Étape 1.* $`X(t)=e^{-itH_{\mathrm{eff}}}\rho e^{itH_{\mathrm{eff}}^*}`$ vérifie
$`X'=\mathcal L_0(X)`$ et $`X(0)=\rho`$, donc $`S_0(t)=e^{t\mathcal L_0}`$ d'après le lemme 8.1 ; il est
CP, avec l'unique opérateur de Kraus $`e^{-itH_{\mathrm{eff}}}`$. *Étape 2.* Notons
$`a=\|\mathcal L_0\|`$, $`b=\|\Phi\|`$. Par récurrence, les $`S_n`$ sont continus et
$`\|S_n(t)\|\le\int_0^te^{(t-s)a}\,b\,e^{sa}\frac{(bs)^{n-1}}{(n-1)!}ds=e^{ta}\frac{(bt)^n}{n!}`$,
donc la série converge uniformément sur tout $`[0,T]`$. *Étape 3.* Comme
$`S_n(t)=S_0(t)\int_0^te^{-s\mathcal L_0}\Phi S_{n-1}(s)\,ds`$, la règle de dérivation d'un produit et
le théorème fondamental de l'analyse (lemme A.4(e)) donnent $`S_n'=\mathcal L_0S_n+\Phi S_{n-1}`$ ($`n\ge1`$) et
$`S_0'=\mathcal L_0S_0`$. Grâce aux majorations, la série des dérivées converge uniformément sur
$`[0,T]`$, donc $`Y(t)=\sum_nS_n(t)`$ est dérivable avec $`Y'=\mathcal L_0Y+\Phi Y=\mathcal LY`$,
$`Y(0)=\mathrm{id}`$, et $`Y(t)=e^{t\mathcal L}`$ d'après le lemme 8.1. *Étape 4.* $`\Phi`$ est CP
(théorème 6.2), donc, par récurrence, l'intégrande $`s\mapsto S_0(t-s)\Phi S_{n-1}(s)`$ est une
famille continue d'applications CP et $`S_n(t)`$ est CP (lemme 8.4) ; il en va de même des sommes
partielles et de leur limite $`e^{t\mathcal L}`$. *Étape 5.* Par cyclicité,
$`\mathrm{Tr}\,\mathcal L(\rho)=\sum_j\big(\mathrm{Tr}(L_j\rho L_j^*)-\mathrm{Tr}(L_j^*L_j\rho)\big)=0`$,
donc
$`\frac d{dt}\mathrm{Tr}\,e^{t\mathcal L}\rho=\mathrm{Tr}\,\mathcal L(e^{t\mathcal L}\rho)=0`$ :
chaque $`e^{t\mathcal L}`$ est TP. La loi de semi-groupe et la continuité découlent du lemme 8.1.
$`\square`$

*Calcul :* le théorème 8.6 pour des sauts arbitraires sur $`M_3`$, y compris la complète positivité de la partie due aux sauts $`\sum_{n\ge1}S_n(t)`$ — [tutoriel, section 8.2](../../tutorial/fr/08_dynamics.md#s-8-2). *Dans la bibliothèque :* $`e^{t\mathcal L}`$ est `dynamics.lindblad_channel`, et sa matrice $`\mathcal L`$ est `dynamics.lindblad_superoperator` — [tutoriel, section 8.2](../../tutorial/fr/08_dynamics.md#s-8-2).

<a id="r-8-7"></a>

**Théorème 8.7 (Gorini–Kossakowski–Sudarshan ; Lindblad).** (1) Réciproquement, le générateur de
tout semi-groupe dynamique quantique sur $`M_D`$ a la forme du théorème 8.6, avec $`J\le D^2-1`$ sauts. On
peut choisir $`\mathrm{Tr}\,H=0`$, $`\mathrm{Tr}\,L_j=0`$ et $`\mathrm{Tr}(L_j^*L_k)=0`$ pour
$`j\neq k`$ ; avec ces normalisations, $`H`$ et $`\Phi=\sum_jL_j\cdot L_j^*`$ sont déterminés par
$`\mathcal L`$, et les $`L_j`$ le sont à $`L'_j=\sum_ku_{jk}L_k`$ près, avec une unitaire $`u`$ (après avoir
complété par des zéros la famille la plus courte).

(2) Le générateur de tout semi-groupe dynamique quantique sur $`\mathcal M=\bigoplus_cM_{k_c}`$ est
```math
\mathcal L(\rho)=\Phi(\rho)-\tfrac12\{\Phi^*(1),\rho\}-i[H,\rho]
```
avec une application CP $`\Phi:\mathcal M\to\mathcal M`$ (opérateurs de Kraus entre les secteurs,
proposition 6.4) et $`H=H^*\in\mathcal M`$ : c'est la forme de la section 8.3, avec un seul saut
$`\Phi_1=\Phi`$ entre secteurs et le taux $`1`$.

*Démonstration.* (1) est le théorème D.4 avec les assertions d'unicité du théorème D.5(b),(c), et (2)
est le théorème D.7, tous deux démontrés dans l'[annexe D](D_gksl.md). $`\square`$

L'étape clé (théorème D.3) est que la matrice de Choi du générateur d'un semi-groupe d'applications CP
est positive sur l'orthogonal du vecteur maximalement intriqué ; la seule propriété de semi-groupe,
avec des applications $`T_t`$ positives, ne suffirait pas (l'exemple qui suit le théorème D.3). Pour (2),
les applications $`T_t\circ E`$, où $`E(X)=\sum_cP_cXP_c`$ est l'espérance conditionnelle, sont CP sur
$`M_D`$ ; la matrice de Choi de la dérivée $`\mathcal L\circ E`$ en $`t=0`$ est positive sur l'orthogonal
des vecteurs $`(1\otimes P_c)\Omega`$, ce qui donne une décomposition de $`\mathcal L\circ E`$ (lemme
D.2), que $`E`$ comprime ensuite sur $`\mathcal M`$. En représentation de Heisenberg (corollaire D.8), le
générateur d'un semi-groupe continu d'applications CP unifères sur $`\mathcal M`$ est
$`A\mapsto\Psi(A)+K^*A+AK`$ avec $`\Psi`$ CP et $`K\in\mathcal M`$. Historiquement, (1) est dû à V. Gorini,
A. Kossakowski, E. C. G. Sudarshan [*J. Math. Phys.* 17, 821 (1976)] et G. Lindblad [*Commun. Math.
Phys.* 48, 119 (1976)], et la forme en représentation de Heisenberg, pour des semi-groupes continus en
norme sur des algèbres de von Neumann quelconques, à E. Christensen et D. E. Evans [*J. London Math.
Soc.* 20, 358 (1979)].

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Un système qui
> interagit faiblement avec un grand environnement dont les corrélations décroissent rapidement
> (l'approximation de Born–Markov) devrait perdre la mémoire de son passé : son état à l'instant
> $`s+t`$ s'obtient alors à partir de l'état à l'instant $`s`$ par un canal qui ne dépend que de $`t`$.
> *Si* c'est le cas, le théorème 8.7 affirme que la dynamique est la somme d'une partie hamiltonienne
> et de « sauts » $`L_j`$ (émission de photons, perte de particules, déphasage), à des taux fixés par la
> normalisation des $`L_j`$.

<a id="s-8-3"></a>

## 8.3 Équations de Lindblad dans une algèbre à secteurs

Sur $`\mathcal M`$, il y a deux sortes de sauts. Un saut **à l'intérieur des secteurs** est un
opérateur $`L_j\in\mathcal M`$. Un saut **entre secteurs** est une application CP
$`\Phi:\mathcal M\to\mathcal M`$ d'opérateurs de Kraus $`K^{dc}_l:\mathbb C^{k_c}\to\mathbb C^{k_d}`$
(proposition 6.4, un `InterSectorChannel`) : $`\Phi(\rho)_d=\sum_{c,l}K^{dc}_l\rho_cK^{dc*}_l`$.
D'après la proposition 6.8, son application duale brute vérifie
$`\Phi^*(1)_c=\sum_{d,l}K^{dc*}_lK^{dc}_l`$, élément diagonal par blocs de $`\mathcal M`$. Avec des taux
$`g_j,g_m>0`$, le générateur de `lindblad_rhs` est
```math
\mathcal L(\rho)=-i[H,\rho]+\sum_jg_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big)
+\sum_mg_m\Big(\Phi_m(\rho)-\tfrac12\{\Phi_m^*(1),\rho\}\Big).
```
D'après la proposition 8.8(1) ci-dessous, il engendre un semi-groupe dynamique quantique sur
$`\mathcal M`$ ; réciproquement, d'après le théorème 8.7(2) (théorème D.7), tout semi-groupe dynamique
quantique sur $`\mathcal M`$ a un générateur de cette forme.

*Dans la bibliothèque :* un saut entre secteurs est un `InterSectorChannel`, et $`\mathcal L(\rho)`$ est `dynamics.lindblad_rhs` — [tutoriel, section 8.3](../../tutorial/fr/08_dynamics.md#s-8-3).

<a id="r-8-8"></a>

**Proposition 8.8.** Soit $`J^{dc}_l\in M_D`$ la matrice ayant le bloc $`K^{dc}_l`$ en position $`(d,c)`$
et des zéros ailleurs.

1. $`\mathcal L`$ envoie $`\mathcal M`$ dans $`\mathcal M`$ et est la restriction à $`\mathcal M`$ du
   générateur GKSL $`\hat{\mathcal L}`$ sur $`M_D`$ de sauts $`\sqrt{g_j}L_j`$ et $`\sqrt{g_m}J^{dc}_{m,l}`$
   (pour tous $`m,l,d,c`$). Par conséquent, $`(e^{t\mathcal L})_{t\ge0}`$ est un semi-groupe dynamique
   quantique sur $`\mathcal M`$.
2. Soit au contraire $`\tilde{\mathcal L}`$ le générateur GKSL sur $`M_D`$ dans lequel les blocs d'un
   même indice de Kraus sont sommés en un seul saut, $`L_{m,l}=\sum_{d,c}J^{dc}_{m,l}`$ (avec le même
   hamiltonien, les mêmes sauts $`\sqrt{g_j}L_j`$ et les sauts $`\sqrt{g_m}L_{m,l}`$). Alors
   $`\mathcal L(\rho)=E\big(\tilde{\mathcal L}(\rho)\big)`$ pour $`\rho\in\mathcal M`$, où
   $`E(X)=\sum_dP_dXP_d`$ est l'espérance conditionnelle de $`M_D`$ sur $`\mathcal M`$ ; en général,
   $`\tilde{\mathcal L}(\rho)\notin\mathcal M`$.
3. Si chaque $`L_{m,l}`$ a au plus un bloc non nul dans chaque ligne de blocs et dans chaque colonne de
   blocs (par exemple s'il envoie seulement le secteur $`c`$ dans le secteur $`c-1`$, comme les
   opérateurs d'annihilation de l'exemple ci-dessous), alors $`\tilde{\mathcal L}(\rho)=\mathcal L(\rho)`$
   pour $`\rho\in\mathcal M`$.

*Démonstration.* Comme $`\rho=\sum_cP_c\rho P_c`$ et $`J^{dc}P_{c'}=\delta_{cc'}J^{dc}`$,
$`J^{dc}\rho J^{dc*}`$ est le bloc $`K^{dc}\rho_cK^{dc*}`$ en position $`(d,d)`$, et $`J^{dc*}J^{dc}`$ est
$`K^{dc*}K^{dc}`$ en position $`(c,c)`$. En sommant sur $`d,l`$, on obtient $`\Phi(\rho)`$ et $`\Phi^*(1)`$,
ce qui démontre l'identité de (1), et tous les termes sont dans $`\mathcal M`$. Ainsi
$`\hat{\mathcal L}`$ et toutes ses puissances envoient $`\mathcal M`$ dans $`\mathcal M`$,
$`e^{t\hat{\mathcal L}}|_{\mathcal M}=e^{t\mathcal L}`$, et c'est un canal de $`\mathcal M`$ d'après le
théorème 8.6 et le lemme 8.4(a). (2) $`L\rho L^*=\sum_{d,d',c}J^{dc}\rho J^{d'c*}`$ pour
$`L=\sum_{d,c}J^{dc}`$, et $`E`$ garde les termes $`d=d'`$. Pour $`B\in\mathcal M`$, $`P_dXBP_d=P_dXP_dB`$,
donc $`E(XB)=E(X)B`$ et $`E(BX)=BE(X)`$ ; ainsi $`E\{L^*L,\rho\}=\{E(L^*L),\rho\}`$, et
$`E(L^*L)=\sum_cP_cL^*LP_c=\sum_{c,d}J^{dc*}J^{dc}`$ car $`J^{dc*}J^{d'c}=0`$ pour $`d\ne d'`$. Enfin
$`E(-i[H,\rho])=-i[H,\rho]`$, et les termes des sauts $`\sqrt{g_j}L_j\in\mathcal M`$ sont dans
$`\mathcal M`$. (3) Soit $`L=\sum_{(d,c)\in S}J^{dc}`$, où l'ensemble $`S`$ de paires contient au plus une
paire de $`d`$ donné et au plus une de $`c`$ donné. Pour $`\rho\in\mathcal M`$,
$`J^{dc}\rho J^{d'c'*}=J^{dc}P_c\rho P_{c'}J^{d'c'*}`$ est nul sauf si $`c=c'`$, et alors $`d=d'`$ ; et
$`J^{dc*}J^{d'c'}`$ est nul sauf si $`d=d'`$, et alors $`c=c'`$. Donc $`L\rho L^*`$ et $`L^*L`$ sont diagonaux par
blocs, tous les termes de $`\tilde{\mathcal L}(\rho)`$ sont dans $`\mathcal M`$, et
$`\tilde{\mathcal L}(\rho)=E(\tilde{\mathcal L}(\rho))=\mathcal L(\rho)`$ d'après (2). $`\square`$

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). La partie (2)
> dit que l'équation de Lindblad sur $`\mathcal M`$ est ce qui reste d'une équation sur $`M_D`$ lorsque
> les cohérences entre secteurs créées par les sauts sont effacées immédiatement — c'est le contenu
> d'une règle de supersélection. Pour la perte de particules, rien n'est effacé (partie (3)) :
> l'équation de Lindblad usuelle avec les opérateurs d'annihilation comme sauts est déjà une équation
> sur $`\mathcal M`$.

*Calcul :* la proposition 8.8(2) pour une application aléatoire entre les secteurs et la proposition 8.8(3) pour la perte de particules, comparées à des générateurs GKSL sur $`M_6`$ — [tutoriel, section 8.3](../../tutorial/fr/08_dynamics.md#s-8-3).

<a id="r-8-9"></a>

**Proposition 8.9 (poids des secteurs).** Soit $`\rho(t)=e^{t\mathcal L}\rho_0`$ et $`p(t)`$ ses poids
des secteurs.

1. S'il n'y a pas de sauts entre secteurs ($`K^{dc}_{m,l}=0`$ pour tous $`m,l`$ et $`d\ne c`$), chaque
   $`p_c(t)`$ est constant.
2. Si, pour tout $`m`$ et toute paire $`d\ne c`$, $`\sum_lK^{dc*}_{m,l}K^{dc}_{m,l}=r_{m,dc}1_{k_c}`$ avec
   des nombres $`r_{m,dc}`$, alors $`p'=Qp`$ avec $`Q_{dc}=\sum_mg_mr_{m,dc}`$ pour $`d\ne c`$ et
   $`Q_{cc}=-\sum_{d\ne c}Q_{dc}`$. Par conséquent, $`p(t)=e^{tQ}p(0)`$ ne dépend ni de $`H`$, ni des
   $`L_j`$, ni des états à l'intérieur des secteurs.

*Démonstration.* $`p_d'=\mathrm{Tr}(P_d\mathcal L(\rho))`$ d'après le lemme 8.1. Le terme hamiltonien s'annule comme
dans le théorème 8.3. Pour $`L_j\in\mathcal M`$,
$`\mathrm{Tr}(P_dL_j\rho L_j^*)-\mathrm{Tr}(P_dL_j^*L_j\rho) =\mathrm{Tr}(L_j^*P_dL_j\rho)-\mathrm{Tr}(L_j^*L_jP_d\rho)=0`$,
puisque $`P_d`$ est central. Pour un saut $`\Phi`$ entre secteurs,
$`\mathrm{Tr}(P_d\Phi(\rho))=\sum_c\mathrm{Tr}\big(\sum_lK^{dc*}_lK^{dc}_l\rho_c\big)`$ et
$`\mathrm{Tr}(P_d\tfrac12\{\Phi^*(1),\rho\})=\sum_{d'}\mathrm{Tr}\big(\sum_lK^{d'd*}_lK^{d'd}_l\rho_d\big)`$ ;
les termes $`c=d`$ et $`d'=d`$ se compensent. Cela démontre (1), et sous l'hypothèse de (2) les termes
restants sont $`\sum_{c\ne d}r_{dc}p_c-\sum_{d'\ne d}r_{d'd}p_d`$. Le lemme 8.1 donne
$`p(t)=e^{tQ}p(0)`$. $`\square`$

Sous l'hypothèse de (2), pour $`d\ne c`$, $`r_{m,dc}`$ est le coefficient $`T_{dc}(\rho)`$ de la matrice
de transition de $`\Phi_m`$ (définie avant la proposition 6.12), pour tout état $`\rho`$ tel que $`p_c>0`$,
puisque
$`\mathrm{Tr}\,\Phi_{m,dc}(\rho_c)=\mathrm{Tr}\big(\sum_lK^{dc*}_{m,l}K^{dc}_{m,l}\rho_c\big)=r_{m,dc}p_c`$.
Ainsi, pour un seul saut de taux $`1`$, $`Q=R-\mathrm{diag}(\sum_dR_{dc})_c`$ avec $`R`$ = `transition_matrix()`
(la diagonale de $`R`$ se compense).

*Dans la bibliothèque :* $`R`$ est `InterSectorChannel.transition_matrix` — [tutoriel, section 8.3](../../tutorial/fr/08_dynamics.md#s-8-3).

**Exemple : perte de particules.** Deux modes bosoniques avec au plus deux particules ; les secteurs
sont $`N=0,1,2`$, de bases $`|00\rangle`$ ; $`|10\rangle,|01\rangle`$ ;
$`|20\rangle,|11\rangle,|02\rangle`$, de sorte que $`\mathcal M=M_1\oplus M_2\oplus M_3`$. Le dimère de
Bose–Hubbard $`H=J(a_1^*a_2+a_2^*a_1)+\frac U2\sum_jn_j(n_j-1)`$ commute avec $`N`$, et chaque mode perd
des particules au taux $`\gamma`$ : un saut entre secteurs d'opérateurs de Kraus $`a_1,a_2`$, restreints
à $`N\to N-1`$. Comme $`a_1^*a_1+a_2^*a_2=N`$, la proposition 8.9(2) s'applique avec
$`Q_{N-1,N}=\gamma N`$, et en partant de $`N=2`$,
```math
p_2(t)=q^2,\qquad p_1(t)=2q(1-q),\qquad p_0(t)=(1-q)^2,\qquad q=e^{-\gamma t},
```
la loi binomiale de paramètres $`2`$ et $`q`$ (la loi du nombre de survivants si chaque particule survit
indépendamment avec probabilité $`q`$), quels que soient les effets de l'effet tunnel entre modes
(hopping) et de l'interaction : les deux membres résolvent $`p'=Qp`$ (dériver) avec $`p(0)=(0,0,1)`$, donc
ils coïncident d'après le lemme 8.1.

*Calcul :* l'exemple, avec $`\Phi^*(1)=0\oplus1\oplus2`$ et $`p(t)`$ comparé aux deux formules — [tutoriel, section 8.3](../../tutorial/fr/08_dynamics.md#s-8-3).

<a id="s-8-4"></a>

## 8.4 Post-sélection et évolution sans saut

Regroupons les sauts du générateur de la section 8.3 en
$`\Phi_{\mathrm{tot}}(\rho)=\sum_jg_jL_j\rho L_j^*+\sum_mg_m\Phi_m(\rho)`$, avec
$`A=\Phi_{\mathrm{tot}}^*(1)=\sum_jg_jL_j^*L_j+\sum_mg_m\Phi_m^*(1)\ge0`$ dans $`\mathcal M`$, et
$`H_{\mathrm{eff}}=H-\frac i2A`$. Appliquons le théorème 8.6 au générateur $`\hat{\mathcal L}`$ sur $`M_D`$
de la proposition 8.8(1) : d'après la démonstration de la proposition 8.8, sa partie de sauts envoie
$`\mathcal M`$ dans $`\mathcal M`$, où elle est égale à $`\Phi_{\mathrm{tot}}`$, et sa somme
$`\sum\hat L^*\hat L`$ sur les sauts est $`A`$, donc son $`H_{\mathrm{eff}}`$ est
$`H-\frac i2A\in\mathcal M`$. Par conséquent ses $`S_n(t)`$ envoient $`\mathcal M`$ dans $`\mathcal M`$ (par
récurrence sur $`n`$), et pour $`\rho\in\mathcal M`$
```math
\rho(t)=\sum_{n\ge0}S_n(t)\rho,\qquad S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},
\qquad S_n(t)=\int_0^tS_0(t-s)\,\Phi_{\mathrm{tot}}\,S_{n-1}(s)\,ds .
```

<a id="r-8-10"></a>

**Proposition 8.10.** Pour un état $`\rho`$, posons $`w_n(t)=\mathrm{Tr}\,S_n(t)\rho`$. Alors
$`w_n(t)\ge0`$, $`\sum_nw_n(t)=1`$, et $`w_0`$ est décroissante avec
$`w_0'(t)=-\mathrm{Tr}\big(A\,S_0(t)\rho\big)`$. Pour un état pur
$`\rho=|\psi\rangle\langle\psi|`$, $`w_0(t)=\|e^{-itH_{\mathrm{eff}}}\psi\|^2`$.

*Démonstration.* $`S_n(t)`$ est CP, donc positive, et $`\sum_nw_n=\mathrm{Tr}\,\rho(t)=1`$. Avec
$`X=S_0(t)\rho\ge0`$,
$`w_0'=\mathrm{Tr}(-iH_{\mathrm{eff}}X+iXH_{\mathrm{eff}}^*)=-\mathrm{Tr}(AX)`$, et
$`\mathrm{Tr}(AX)=\mathrm{Tr}(X^{1/2}AX^{1/2})\ge0`$. La dernière assertion est
$`\mathrm{Tr}|\phi\rangle\langle\phi|=\|\phi\|^2`$. $`\square`$

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Supposons qu'un détecteur enregistre chaque saut — chaque photon
> émis, chaque atome perdu. Alors $`w_n(t)`$ est la probabilité d'exactement $`n`$ clics jusqu'à
> l'instant $`t`$, $`S_n(t)\rho/w_n(t)`$ l'état sachant qu'il y a eu $`n`$ clics, et $`S_0(t)\rho/w_0(t)`$
> l'état sachant qu'il n'y a eu **aucun clic** : un résultat nul modifie lui aussi l'état, par
> l'intermédiaire du $`H_{\mathrm{eff}}`$ non hermitien. Cette interprétation en termes de processus
> de comptage [E. B. Davies, *Quantum Theory of Open Systems*, 1976 ; M. D. Srinivas et E. B.
> Davies, *Optica Acta* 28, 981 (1981)] est à la base de la méthode des trajectoires quantiques
> (fonctions d'onde Monte-Carlo) [J. Dalibard, Y. Castin, K. Mølmer, *Phys. Rev. Lett.* 68, 580
> (1992)] ; nous n'en utilisons ici que la partie mathématique.

Pour des pertes, les termes de Dyson sont les composantes de secteur de $`\rho(t)`$ (dans le langage de
l'encadré : compter les clics revient à mesurer la charge). La **mise à jour de Lüders** d'un état
$`\rho`$ par un projecteur $`P`$ tel que $`\mathrm{Tr}(P\rho)>0`$ est l'état
$`P\rho P/\mathrm{Tr}(P\rho P)`$ (noter que $`\mathrm{Tr}(P\rho P)=\mathrm{Tr}(P\rho)`$).

*Dans la bibliothèque :* la mise à jour de Lüders est `states.lueders_update` — [tutoriel, section 8.4](../../tutorial/fr/08_dynamics.md#s-8-4).

<a id="r-8-11"></a>

**Proposition 8.11.** Supposons les secteurs indexés par $`N=0,\dots,C-1`$, sans sauts à l'intérieur
des secteurs, chaque saut abaissant la charge d'une unité ($`K^{dc}_{m,l}=0`$ sauf si $`d=c-1`$), et
soit $`\rho=P_{N_0}\rho P_{N_0}`$. Alors $`S_n(t)\rho=P_{N_0-n}\,\rho(t)\,P_{N_0-n}`$ pour
$`0\le n\le N_0`$ et $`S_n(t)\rho=0`$ pour $`n>N_0`$. En particulier, $`p_{N_0}(t)=w_0(t)>0`$, et la mise à
jour de Lüders de
$`\rho(t)`$ par $`P_{N_0}`$ (« aucune particule perdue ») est l'état sans saut :
```math
\frac{P_{N_0}\rho(t)P_{N_0}}{p_{N_0}(t)}=\frac{e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*}}{w_0(t)},
\qquad p_{N_0}(t)=w_0(t).
```

*Démonstration.* $`S_0(s)`$ envoie les opérateurs à support dans le secteur $`N`$ sur des opérateurs à
support dans le secteur $`N`$, car $`H_{\mathrm{eff}}\in\mathcal M`$ ; $`\Phi_{\mathrm{tot}}`$ envoie le
secteur $`N\ge1`$ sur le secteur $`N-1`$ et annule le secteur $`0`$. Par récurrence, $`S_n(s)\rho`$ est à
support dans le secteur $`N_0-n`$ pour tout $`s`$, et égal à $`0`$ pour $`n>N_0`$. Comprimer
$`\rho(t)=\sum_nS_n(t)\rho`$ par $`P_{N_0-n}`$ conserve exactement le $`n`$-ième terme. Enfin, $`w_0(t)>0`$, car $`e^{-itH_{\mathrm{eff}}}`$ est
inversible (lemme 8.1), donc $`S_0(t)\rho`$ est une matrice positive non nulle. $`\square`$

*Calcul :* la proposition 8.11 pour $`n=0`$ et $`n=1`$, pour deux modes qui fuient à des taux différents — [tutoriel, section 8.4](../../tutorial/fr/08_dynamics.md#s-8-4).

<a id="s-8-5"></a>

## 8.5 États stationnaires et relaxation

Un état $`\sigma`$ est **stationnaire** si $`e^{t\mathcal L}\sigma=\sigma`$ pour tout $`t\ge0`$, ce qui
équivaut à $`\mathcal L(\sigma)=0`$ (dériver en $`t=0`$ ; réciproquement, utiliser la série).

<a id="r-8-12"></a>

**Proposition 8.12.** Tout semi-groupe dynamique quantique sur $`\mathcal M`$ admet un état
stationnaire.

*Démonstration.* Fixons un état $`\rho`$ et posons $`\sigma_\tau=\frac1\tau\int_0^\tau T_s\rho\,ds`$.
C'est un état d'après le lemme A.4(d) (les états forment un convexe fermé). D'après le lemme
A.4(b),(e),
$`T_t\sigma_\tau-\sigma_\tau=\frac1\tau\big(\int_\tau^{\tau+t}-\int_0^t\big)T_s\rho\,ds`$, dont la
norme trace est au plus $`2t/\tau`$ d'après le lemme A.4(c), puisque $`\|T_s\rho\|_1=1`$. Les états forment un fermé borné d'un
espace de dimension finie, donc une suite $`\sigma_{\tau_m}`$, $`\tau_m\to\infty`$, converge vers un
état $`\sigma`$ (Bolzano–Weierstrass), et
$`T_t\sigma-\sigma=\lim_m(T_t\sigma_{\tau_m}-\sigma_{\tau_m})=0`$. $`\square`$

<a id="r-8-13"></a>

**Corollaire 8.13.** S'il n'y a pas de sauts entre secteurs (comme dans la proposition 8.9(1)), les
états stationnaires engendrent un espace de dimension au moins $`C`$ : chaque secteur porte son propre
état stationnaire.

*Démonstration.* $`\mathcal L`$ envoie $`P_c\mathcal MP_c\cong M_{k_c}`$ dans lui-même ($`H`$ et les $`L_j`$
sont dans $`\mathcal M`$, et $`\Phi_m(X)=\sum_lK^{cc}_{m,l}X_cK^{cc*}_{m,l}`$ et $`\Phi_m^*(1)`$ sont
diagonaux par blocs), donc chaque $`e^{t\mathcal L}`$ aussi. Sa restriction est CP (lemme 8.4(a)) et TP,
c'est donc un semi-groupe dynamique quantique sur $`M_{k_c}`$, qui admet un état stationnaire
$`\sigma_c`$ d'après la proposition 8.12 ; en tant qu'élément de $`\mathcal M`$, $`\sigma_c`$ est un état
stationnaire de $`e^{t\mathcal L}`$. Les $`\sigma_c`$ sont dans des blocs différents, donc linéairement
indépendants. $`\square`$

Le théorème D.6 (le cas de dimension finie d'un critère de A. Frigerio) donne une condition
suffisante d'unicité : si un générateur GKSL sur $`M_k`$ de données $`(H;L_1,\dots,L_r)`$ a un état stationnaire
fidèle $`\sigma>0`$ et si les seules matrices qui commutent avec $`H`$ et avec tous les $`L_j`$, $`L_j^*`$
sont les multiples de $`1`$, alors $`\sigma`$ est l'unique état stationnaire. Sans sauts entre secteurs,
la restriction de $`\mathcal L`$ à $`P_c\mathcal MP_c\cong M_{k_c}`$ est le générateur GKSL de données
$`(H_c;\sqrt{g_1}L_{1,c},\sqrt{g_2}L_{2,c},\dots)`$, de sorte que le théorème s'applique secteur par
secteur. L'exercice 4 montre que l'unicité peut être en défaut même dans
un facteur : la seconde hypothèse y est en défaut.

*Calcul :* le corollaire 8.13 et les hypothèses du théorème D.6 pour des sauts aléatoires à l'intérieur des secteurs de $`M_2\oplus M_3`$ — [tutoriel, section 8.5](../../tutorial/fr/08_dynamics.md#s-8-5).

<a id="r-8-14"></a>

**Proposition 8.14 (émission spontanée).** Sur $`M_2`$, de base $`|g\rangle,|e\rangle`$, soient
$`H=\omega|e\rangle\langle e|`$, un saut $`L=|g\rangle\langle e|`$ et un taux $`\gamma>0`$. Alors
```math
\rho_{ee}(t)=e^{-\gamma t}\rho_{ee}(0),\qquad\rho_{ge}(t)=e^{(i\omega-\gamma/2)t}\rho_{ge}(0),
\qquad e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}},
```
où $`\mathcal U_t(\rho)=e^{-itH}\rho e^{itH}`$ et $`\mathcal A_p`$ est le canal d'amortissement
d'amplitude de paramètre $`p`$ (section 6.4). Le spectre de $`\mathcal L`$ est
$`\{0,-\gamma,-\gamma/2\pm i\omega\}`$, $`|g\rangle\langle g|`$ est l'unique état stationnaire, et
$`\tfrac12\|\rho(t)-|g\rangle\langle g|\|_1=\big(\rho_{ee}(0)^2e^{-2\gamma t}+|\rho_{ge}(0)|^2e^{-\gamma t}\big)^{1/2}`$.

*Démonstration.* Ici $`L^*L=|e\rangle\langle e|`$ et $`L\rho L^*=\rho_{ee}|g\rangle\langle g|`$, donc,
coefficient par coefficient, $`\rho_{ee}'=-\gamma\rho_{ee}`$, $`\rho_{gg}'=\gamma\rho_{ee}`$ et
$`\rho_{ge}'=(i\omega-\gamma/2)\rho_{ge}`$, ce qui donne la solution. Dans la base
$`(\rho_{gg},\rho_{ge},\rho_{eg},\rho_{ee})`$, la matrice de $`\mathcal L`$ est triangulaire supérieure,
de diagonale $`(0,-\gamma/2+i\omega,-\gamma/2-i\omega,-\gamma)`$. $`\mathcal L(\sigma)=0`$ impose
$`\sigma_{ee}=\sigma_{ge}=0`$. Le canal d'amortissement d'amplitude multiplie $`\rho_{ee}`$ par $`1-p`$ et
$`\rho_{ge}`$ par $`\sqrt{1-p}`$ (ses opérateurs de Kraus sont
$`|g\rangle\langle g|+\sqrt{1-p}\,|e\rangle\langle e|`$ et $`\sqrt p\,|g\rangle\langle e|`$), et
$`\mathcal U_t`$ multiplie $`\rho_{ge}`$ par $`e^{i\omega t}`$ ; tous deux sont TP, ce qui fixe
$`\rho_{gg}`$. Enfin, $`\rho(t)-|g\rangle\langle g|`$ est autoadjoint, de trace nulle et de diagonale
$`(-\rho_{ee}(t),\rho_{ee}(t))`$, donc ses valeurs propres sont
$`\pm(\rho_{ee}(t)^2+|\rho_{ge}(t)|^2)^{1/2}`$. $`\square`$

La cohérence décroît à un taux moitié de celui de la population : avec $`T_1=1/\gamma`$ et
$`T_2=2/\gamma`$, c'est la relation $`\frac1{T_2}=\frac1{2T_1}+\frac1{T_\varphi}`$ de la section 6.4 sans
le terme de déphasage. Le plus petit taux de décroissance non nul, le **trou spectral**
$`\min\{-\mathrm{Re}\,\mu:\mu\ne0\text{ valeur propre de }\mathcal L\}=\gamma/2`$, est ici le taux
de relaxation : d'après la dernière formule, la distance à $`|g\rangle\langle g|`$ est au plus
$`e^{-\gamma t/2}`$ (puisque $`|\rho_{ge}|^2\le\rho_{gg}\rho_{ee}`$ par positivité,
$`\rho_{ee}^2+|\rho_{ge}|^2\le\rho_{ee}\le1`$), et $`e^{\gamma t/2}`$ fois cette distance tend vers
$`|\rho_{ge}(0)|`$ quand $`t\to\infty`$.

*Calcul :* le spectre de $`\mathcal L`$, la factorisation $`e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}}`$ et la distance de trace de la proposition 8.14, avec `dynamics.lindblad_superoperator` et `dynamics.lindblad_channel` — [tutoriel, section 8.5](../../tutorial/fr/08_dynamics.md#s-8-5).

<a id="s-8-6"></a>

## 8.6 Points de comparaison

*Calcul :* la comparaison avec d'autres logiciels (QuTiP, dynamiqs, QuantumOptics.jl, SciPy) — [tutoriel, section 8.6](../../tutorial/fr/08_dynamics.md#s-8-6).

<a id="s-8-7"></a>

## 8.7 Pas de temps qui sont des canaux

Un solveur numérique remplace $`e^{t\mathcal L}`$ par $`n`$ pas d'une application $`R_h`$, $`h=t/n`$. Chaque
$`e^{h\mathcal L}`$ est un canal de $`\mathcal M`$ (proposition 8.8(1)) ; si chaque pas est lui aussi un canal, les densités
calculées sont des états, quel que soit $`h`$. Pour $`\rho'=\mathcal L(\rho)`$, la méthode de Runge–Kutta classique (RK4)
est le pas $`\rho\mapsto\rho+\frac h6(k_1+2k_2+2k_3+k_4)`$ avec $`k_1=\mathcal L(\rho)`$, $`k_2=\mathcal L(\rho+\frac h2k_1)`$,
$`k_3=\mathcal L(\rho+\frac h2k_2)`$ et $`k_4=\mathcal L(\rho+hk_3)`$. En insérant chaque étage dans le suivant,
$`k_2=(\mathcal L+\frac h2\mathcal L^2)\rho`$, $`k_3=(\mathcal L+\frac h2\mathcal L^2+\frac{h^2}4\mathcal L^3)\rho`$ et
$`k_4=(\mathcal L+h\mathcal L^2+\frac{h^2}2\mathcal L^3+\frac{h^3}4\mathcal L^4)\rho`$, donc
```math
R_h=T_4(h\mathcal L),\qquad T_4(z)=1+z+\frac{z^2}2+\frac{z^3}6+\frac{z^4}{24}.
```
$`R_h`$ est TP, car $`\mathrm{Tr}\,\mathcal L(X)=0`$ pour tout $`X`$ (étape 5 de la démonstration du théorème 8.6, avec la
proposition 8.8(1)), et envoie les éléments autoadjoints sur des éléments autoadjoints, car $`\mathcal L(X^*)=\mathcal L(X)^*`$.
Elle diffère de $`e^{h\mathcal L}`$ par les termes d'ordre $`\ge5`$ de la série exponentielle, mais elle n'est pas
forcément positive.

<a id="r-8-15"></a>

**Proposition 8.15 (RK4 n'est pas positif).** Sur $`M_2`$, de base $`|g\rangle,|e\rangle`$, soient $`\mathcal L`$ le
générateur de la proposition 8.14 ($`H=\omega|e\rangle\langle e|`$ avec $`\omega\in\mathbb R`$, un saut $`|g\rangle\langle e|`$ de
taux $`\gamma>0`$), $`\lambda=i\omega-\gamma/2`$, $`h>0`$, $`\tau=T_4(-h\gamma)`$ et $`a=T_4(h\lambda)`$.

1. $`R_h(\rho)=\begin{pmatrix}\rho_{gg}+(1-\tau)\rho_{ee}&a\rho_{ge}\\\bar a\rho_{eg}&\tau\rho_{ee}\end{pmatrix}`$, et les
   assertions suivantes sont équivalentes : $`R_h`$ est positive ; $`R_h`$ est CP ; $`\tau\le1`$ et $`|a|^2\le\tau`$.
2. $`\tau\le1`$ si et seulement si $`h\gamma\le x_*`$, où $`x_*=2.7853\ldots`$ est l'unique racine réelle de
   $`x^3-4x^2+12x-24`$. Pour $`h\gamma>x_*`$, l'état $`|e\rangle\langle e|`$ est envoyé sur $`\mathrm{diag}(1-\tau,\tau)`$, dont
   le premier coefficient est négatif.
3. Avec $`u=\gamma^2+4\omega^2`$,
   ```math
   \tau-|a|^2=h^5\Big(-\frac{\gamma u(4\omega^2-3\gamma^2)}{384}+\frac{u^2(4\omega^2-5\gamma^2)}{4608}\,h
   +\frac{\gamma u^3}{9216}\,h^2-\frac{u^4}{147456}\,h^3\Big).
   ```
   Par conséquent, si $`2\omega>\sqrt3\,\gamma`$, il existe $`h_0>0`$ tel que $`R_h`$ n'est pas positive pour $`0<h<h_0`$ ; si
   $`\omega=\gamma`$, $`R_h`$ n'est positive pour aucun $`h>0`$ ; si $`\omega=0`$, $`R_h`$ est CP si et seulement si
   $`h\gamma\le x_*`$.

*Démonstration.* (1) D'après la démonstration de la proposition 8.14, $`\mathcal L(|g\rangle\langle g|)=0`$,
$`\mathcal L(|e\rangle\langle e|)=\gamma(|g\rangle\langle g|-|e\rangle\langle e|)`$, $`\mathcal L(|g\rangle\langle e|)=\lambda|g\rangle\langle e|`$
et $`\mathcal L(|e\rangle\langle g|)=\bar\lambda|e\rangle\langle g|`$. Par récurrence,
$`\mathcal L^n(|e\rangle\langle e|)=(-\gamma)^n(|e\rangle\langle e|-|g\rangle\langle g|)`$ pour $`n\ge1`$, donc
$`R_h(|e\rangle\langle e|)=|e\rangle\langle e|+(\tau-1)(|e\rangle\langle e|-|g\rangle\langle g|)`$, tandis que $`R_h`$ fixe
$`|g\rangle\langle g|`$ et multiplie $`|g\rangle\langle e|`$ par $`T_4(h\lambda)=a`$ et $`|e\rangle\langle g|`$ par
$`T_4(h\bar\lambda)=\bar a`$ ($`T_4`$ est à coefficients réels). C'est la formule. Si $`\tau\le1`$ et $`|a|^2\le\tau`$, les
opérateurs de Kraus $`|g\rangle\langle g|+\bar a|e\rangle\langle e|`$, $`\sqrt{1-\tau}\,|g\rangle\langle e|`$ et
$`\sqrt{\tau-|a|^2}\,|e\rangle\langle e|`$ donnent $`R_h`$ (comparer les quatre coefficients), donc $`R_h`$ est CP
(théorème 6.2), en particulier positive. Réciproquement, soit $`R_h`$ positive. Alors $`R_h(|e\rangle\langle e|)\ge0`$ donne
$`1-\tau\ge0`$. Pour $`\psi=c|g\rangle+s|e\rangle`$ avec $`c,s>0`$, $`c^2+s^2=1`$, la matrice $`R_h(|\psi\rangle\langle\psi|)\ge0`$ a
pour déterminant $`s^2\big(c^2(\tau-|a|^2)+\tau(1-\tau)s^2\big)\ge0`$ ; en faisant tendre $`s`$ vers $`0`$, on obtient
$`\tau-|a|^2\ge0`$.

(2) Avec $`x=h\gamma`$, $`\tau-1=\frac x{24}q(x)`$, $`q(x)=x^3-4x^2+12x-24`$. Comme $`q'(x)=3x^2-8x+12>0`$ (son
discriminant vaut $`64-144<0`$), $`q`$ est strictement croissante, et $`q(2)=-8<0<3=q(3)`$ ; donc $`q`$ a exactement une
racine réelle $`x_*\in(2,3)`$, et $`\tau\le1`$ si et seulement si $`x\le x_*`$.

(3) $`\tau-|a|^2=T_4(-h\gamma)-T_4(h\lambda)T_4(h\bar\lambda)`$ est un polynôme en $`h`$ de degré $`8`$. Le coefficient de
$`h^n`$ dans $`T_4(h\lambda)T_4(h\bar\lambda)`$ est $`\sum\lambda^j\bar\lambda^l/(j!\,l!)`$ sur $`j+l=n`$, $`0\le j,l\le4`$. Pour
$`n\le4`$, tous ces couples apparaissent, et la somme vaut $`(\lambda+\bar\lambda)^n/n!=(-\gamma)^n/n!`$ par la formule du
binôme, le coefficient de $`h^n`$ dans $`T_4(-h\gamma)`$ ; donc les coefficients de $`h^0,\dots,h^4`$ dans $`\tau-|a|^2`$
s'annulent. Pour $`n=5,\dots,8`$, en fonction de $`\sigma=\lambda+\bar\lambda=-\gamma`$ et $`\pi=\lambda\bar\lambda=u/4`$, les
sommes valent $`\frac{\pi(\lambda^3+\bar\lambda^3)}{24}+\frac{\pi^2\sigma}{12}=\frac{\pi\sigma^3-\pi^2\sigma}{24}`$,
$`\frac{\pi^2(\lambda^2+\bar\lambda^2)}{48}+\frac{\pi^3}{36}=\frac{\pi^2\sigma^2}{48}-\frac{\pi^3}{72}`$, $`\frac{\pi^3\sigma}{144}`$ et
$`\frac{\pi^4}{576}`$ (avec $`\lambda^2+\bar\lambda^2=\sigma^2-2\pi`$, $`\lambda^3+\bar\lambda^3=\sigma^3-3\pi\sigma`$), et
$`T_4(-h\gamma)`$ n'a pas de tels termes ; en substituant $`\sigma`$ et $`\pi`$, on obtient la formule. Si
$`2\omega>\sqrt3\gamma`$, le coefficient de $`h^5`$ est négatif, donc $`\tau-|a|^2<0`$ pour tout $`h>0`$ petit, et $`R_h`$ n'est
pas positive d'après (1). Si $`\omega=\gamma`$, alors $`u=5\gamma^2`$, et avec $`y=5h\gamma`$ la formule devient
$`\tau-|a|^2=-\frac{5(h\gamma)^5}{147456}\,p(y)`$, $`p(y)=y^3-16y^2+32y+384`$. Or $`p(y)>0`$ pour $`y>0`$ :
$`p(y)=y(y-8)^2+(384-32y)`$ est une somme de deux termes positifs ou nuls, non tous deux nuls, pour $`y\le12`$ ; pour
$`12\le y\le16`$, $`y(y-8)^2\ge192`$ (croissant pour $`y\ge8`$) et $`384-32y\ge-128`$ ; et $`p(y)=y^2(y-16)+32y+384>0`$ pour
$`y\ge16`$. Donc $`\tau<|a|^2`$ pour tout $`h>0`$. Si $`\omega=0`$, alors $`u=\gamma^2`$ et, avec $`x=h\gamma`$,
$`\tau-|a|^2=x^5\big(\frac1{128}-\frac{5x}{4608}+\frac{x^2}{9216}-\frac{x^3}{147456}\big)\ge
x^5\big(\frac1{128}-\frac{15}{4608}-\frac{27}{147456}\big)>0`$ pour $`0<x\le3`$ ; avec (1) et (2), $`R_h`$ est CP si et
seulement si $`x\le x_*`$. $`\square`$

Ainsi, pour un système à deux niveaux faiblement amorti ($`2\omega>\sqrt3\,\gamma`$), le pas de RK4 n'est pas positif, si
petit que soit le pas, alors que ni l'amortissement seul ($`\omega=0`$, d'après (3)) ni la rotation seule n'ont ce défaut pour
de petits pas : pour $`\gamma=0`$, la démonstration de (1) s'applique avec $`\tau=1`$, et
$`|T_4(ix)|^2=1-\frac{x^6}{72}+\frac{x^8}{576}\le1`$ pour $`x^2\le8`$ (développer). La violation est petite,
$`\tau-|a|^2=O(h^5)`$, mais ce n'est pas un effet d'arrondi, et la conservation de la trace par $`R_h`$ n'y remédie pas. La
théorie de la stabilité de RK4 (la région où $`|T_4(z)|\le1`$) ne la voit pas non plus : ici $`h\lambda`$ et $`-h\gamma`$ sont
dans cette région pour $`h`$ petit.

*Calcul :* la proposition 8.15 avec `lindblad_evolve(..., method='rk4')` : le pas comme $`T_4(h\mathcal L)`$, un état pur
envoyé sur une matrice ayant une valeur propre négative, et une matrice de Choi ayant une valeur propre négative pour la
perte de particules de la section 8.3 — [tutoriel, section 8.7](../../tutorial/fr/08_dynamics.md#s-8-7).

Le schéma suivant est un canal pour toute taille de pas. Il utilise la décomposition
$`\mathcal L=\mathcal L_0+\Phi_{\mathrm{tot}}`$ de la section 8.4, $`\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)`$
avec $`H_{\mathrm{eff}}=H-\frac i2A`$ et $`A=\Phi_{\mathrm{tot}}^*(1)`$, et le dual $`\Theta^*`$ d'un superopérateur $`\Theta`$ sur
$`\mathcal M`$ pour la trace, $`\mathrm{Tr}(\Theta(X)Z)=\mathrm{Tr}(X\Theta^*(Z))`$. Pour les superopérateurs $`\Theta`$ sur
$`\mathcal M`$, on note $`\|\Theta\|_1=\max\{\|\Theta(X)\|_1:X\in\mathcal M,\ \|X\|_1\le1\}`$, avec la norme de trace de la
section 1.4 ; cette norme est sous-multiplicative.

<a id="r-8-16"></a>

**Théorème 8.16 (un schéma de décomposition dont les pas sont des canaux).** Pour $`h\ge0`$, soient
$`E_h=e^{-ihH_{\mathrm{eff}}/2}`$, $`S_h(\rho)=E_h\rho E_h^*`$ et
```math
\Psi_h=S_h\circ\Big(\mathrm{id}+h\Phi_{\mathrm{tot}}+\frac{h^2}2\Phi_{\mathrm{tot}}^2\Big)\circ S_h,\qquad M_h=\Psi_h^*(1),\qquad
\tilde\Psi_h(\rho)=\Psi_h\big(M_h^{-1/2}\rho M_h^{-1/2}\big).
```

1. $`M_h\in\mathcal M`$ et $`M_h\ge e^{ihH_{\mathrm{eff}}^*}e^{-ihH_{\mathrm{eff}}}`$, qui est définie positive ; $`\tilde\Psi_h`$
   est un canal de $`\mathcal M`$ pour tout $`h\ge0`$.
2. (*Erreur locale.*) Il existe une constante $`C`$, dépendant de $`\mathcal L`$, telle que
   $`\|\tilde\Psi_h-e^{h\mathcal L}\|_1\le Ch^3`$ pour $`0\le h\le1`$.
3. (*Erreur globale.*) Pour tout état $`\rho`$, tout $`t>0`$ et tout entier $`n\ge t`$, avec $`h=t/n`$,
   $`\|\tilde\Psi_h^n(\rho)-e^{t\mathcal L}\rho\|_1\le Cth^2`$.

*Démonstration.* *Étape 0.* Si $`\Theta`$ est une application positive et TP sur $`\mathcal M`$ et $`Y=Y^*\in\mathcal M`$,
alors $`\|\Theta(Y)\|_1\le\|Y\|_1`$. En effet, écrivons $`Y=Y_+-Y_-`$, où $`Y_\pm\ge0`$ sont les parties de la
décomposition spectrale aux valeurs propres positives et aux valeurs propres négatives. Les valeurs singulières de $`Y`$
sont les valeurs absolues de ses valeurs propres ($`|Y|=(Y^2)^{1/2}`$), donc $`\|Y\|_1=\mathrm{Tr}\,Y_++\mathrm{Tr}\,Y_-`$,
et, par l'inégalité triangulaire et $`\|P\|_1=\mathrm{Tr}\,P`$ pour $`P\ge0`$ (section 1.4),
$`\|\Theta(Y)\|_1\le\mathrm{Tr}\,\Theta(Y_+)+\mathrm{Tr}\,\Theta(Y_-)=\|Y\|_1`$. Un $`X\in\mathcal M`$ quelconque s'écrit
$`Y_1+iY_2`$ avec les autoadjoints $`Y_1=\frac12(X+X^*)`$, $`Y_2=\frac1{2i}(X-X^*)`$ de norme de trace $`\le\|X\|_1`$
(car $`\|X^*\|_1=\|X\|_1`$, chapitre 1, exercice 2) ; donc $`\|\Theta\|_1\le2`$.

*Étape 1 : (1).* $`H_{\mathrm{eff}}\in\mathcal M`$ (section 8.4), donc $`E_h\in\mathcal M`$ (tous les termes de sa série
exponentielle sont dans le sous-espace $`\mathcal M`$). $`\Phi_{\mathrm{tot}}`$ envoie $`\mathcal M`$ dans $`\mathcal M`$ et est CP,
d'opérateurs de Kraus $`\sqrt{g_j}L_j`$ et $`\sqrt{g_m}K^{dc}_{m,l}`$ (proposition 6.4). Donc
$`\mathrm{id}+h\Phi_{\mathrm{tot}}+\frac{h^2}2\Phi_{\mathrm{tot}}^2`$ est CP (lemme 8.4(a)), $`S_h`$ est CP d'opérateur de Kraus
$`E_h`$, et $`\Psi_h`$ est CP. Par cyclicité de la trace, $`S_h^*(Z)=E_h^*ZE_h`$, et
$`\Psi_h^*=S_h^*\circ(\mathrm{id}+h\Phi_{\mathrm{tot}}^*+\frac{h^2}2(\Phi_{\mathrm{tot}}^*)^2)\circ S_h^*`$, où
$`\Phi_{\mathrm{tot}}^*(Z)=\sum_jg_jL_j^*ZL_j+\sum_mg_m\Phi_m^*(Z)`$ (exercice 1) est positive (proposition 6.8 : les
opérateurs de Kraus de $`\Phi_m^*`$ sont les $`K^{dc*}_{m,l}`$). Avec $`W=S_h^*(1)=E_h^*E_h\ge0`$, le facteur central envoie
$`W`$ sur $`W`$ plus un élément positif, et $`S_h^*`$ préserve l'ordre, donc
$`M_h\ge S_h^*(W)=(E_h^2)^*E_h^2=e^{ihH_{\mathrm{eff}}^*}e^{-ihH_{\mathrm{eff}}}`$ ($`E_h^2=e^{-ihH_{\mathrm{eff}}}`$ d'après le
lemme 8.1). Cette matrice est définie positive : $`\langle x,(E_h^2)^*E_h^2x\rangle=\|E_h^2x\|^2>0`$ pour $`x\ne0`$, car
$`E_h^2`$ a pour inverse $`e^{ihH_{\mathrm{eff}}}`$ (lemme 8.1). Toutes les valeurs propres de l'autoadjoint $`M_h\in\mathcal M`$
sont donc positives, et $`N_h=M_h^{-1/2}\in\mathcal M`$ (calcul fonctionnel, section 1.3) vérifie $`N_hM_hN_h=1`$.
L'application $`\rho\mapsto N_h\rho N_h`$ est CP, donc $`\tilde\Psi_h`$ aussi (lemme 8.4(a)), et
$`\mathrm{Tr}\,\tilde\Psi_h(\rho)=\mathrm{Tr}(N_h\rho N_hM_h)=\mathrm{Tr}(\rho N_hM_hN_h)=\mathrm{Tr}\,\rho`$.

*Étape 2 : coefficients de Taylor.* Posons $`\alpha=\frac12\mathcal L_0`$ et $`\beta=\Phi_{\mathrm{tot}}`$, superopérateurs sur
$`\mathcal M`$. Comme à l'étape 1 de la démonstration du théorème 8.6, $`X(s)=E_s\rho E_s^*`$ résout $`X'=\alpha(X)`$ (dériver
$`E_s=e^{-isH_{\mathrm{eff}}/2}`$ grâce au lemme 8.1), donc $`S_h=e^{h\alpha}`$ d'après le lemme 8.1, et
$`\Psi_h=e^{h\alpha}(1+h\beta+\frac{h^2}2\beta^2)e^{h\alpha}`$, tandis que $`e^{h\mathcal L}=e^{h(2\alpha+\beta)}`$. En
développant les séries absolument convergentes, $`\Psi_h=\sum_{n\ge0}h^nP_n`$ et $`e^{h\mathcal L}=\sum_{n\ge0}h^nQ_n`$ avec
$`\sum_n\|P_n\|_1\le e^{2\|\alpha\|_1}(1+\|\beta\|_1+\frac12\|\beta\|_1^2)`$ et $`\sum_n\|Q_n\|_1\le e^{\|2\alpha+\beta\|_1}`$. Ici
$`P_0=Q_0=\mathrm{id}`$, $`P_1=2\alpha+\beta=Q_1`$, et
$`P_2=\frac{\alpha^2}2+\alpha\beta+\alpha^2+\beta\alpha+\frac{\beta^2}2+\frac{\alpha^2}2=\frac12(2\alpha+\beta)^2=Q_2`$. Par
conséquent, pour $`0\le h\le1`$, $`\|\Psi_h-e^{h\mathcal L}\|_1\le C_1h^3`$ avec $`C_1=\sum_{n\ge3}(\|P_n\|_1+\|Q_n\|_1)`$, et
$`\|\Psi_h\|_1\le K=\sum_n\|P_n\|_1`$.

*Étape 3 : (2).* Soit $`x`$ un vecteur propre unitaire de $`M_h`$, de valeur propre $`\mu`$. Comme $`M_h`$ est diagonale par
blocs, on peut choisir $`x`$ dans un seul secteur, de sorte que $`X=|x\rangle\langle x|\in\mathcal M`$, $`\|X\|_1=1`$. Comme
$`e^{h\mathcal L}`$ est TP, $`\mu-1=\mathrm{Tr}(XM_h)-\mathrm{Tr}\,X=\mathrm{Tr}\,\Psi_h(X)-\mathrm{Tr}\,e^{h\mathcal L}(X)`$, et
$`|\mathrm{Tr}\,Y|\le\|Y\|_1`$ (section 1.4) donne $`|\mu-1|\le C_1h^3`$. Soit $`h_1=\min\{1,(2C_1)^{-1/3}\}`$. Pour $`h\le h_1`$,
tous les $`\mu\ge\frac12`$, et $`N_h`$ a les valeurs propres $`\mu^{-1/2}\le\sqrt2`$ avec
$`|\mu^{-1/2}-1|=\frac{|1-\mu|}{\sqrt\mu(1+\sqrt\mu)}\le|\mu-1|`$, car $`\sqrt\mu(1+\sqrt\mu)\ge2^{-1/2}(1+2^{-1/2})>1`$. La
norme d'opérateur d'une matrice autoadjointe est sa plus grande valeur propre en valeur absolue, donc
$`\|N_h-1\|\le C_1h^3`$ et $`\|N_h\|\le\sqrt2`$. Or
```math
\tilde\Psi_h(X)-e^{h\mathcal L}(X)=\Psi_h\big((N_h-1)XN_h+X(N_h-1)\big)+\big(\Psi_h-e^{h\mathcal L}\big)(X),
```
et $`\|AB\|_1\le\|A\|\,\|B\|_1`$, $`\|BA\|_1\le\|A\|\,\|B\|_1`$ (chapitre 1, exercice 2 ; la seconde en passant aux adjoints)
donnent $`\|\tilde\Psi_h-e^{h\mathcal L}\|_1\le\big((1+\sqrt2)K+1\big)C_1h^3`$ pour $`h\le h_1`$. Pour $`h_1<h\le1`$, les deux
applications sont positives et TP, donc $`\|\tilde\Psi_h-e^{h\mathcal L}\|_1\le4\le4h^3/h_1^3`$ d'après l'étape 0. On prend
$`C=\max\{((1+\sqrt2)K+1)C_1,\,4h_1^{-3}\}`$.

*Étape 4 : (3).* La somme
$`\tilde\Psi_h^n-e^{nh\mathcal L}=\sum_{k=0}^{n-1}\tilde\Psi_h^{\,n-1-k}\big(\tilde\Psi_h-e^{h\mathcal L}\big)e^{kh\mathcal L}`$
est télescopique. Appliqué à $`\rho`$, le $`k`$-ième terme est $`\tilde\Psi_h^{\,n-1-k}(Y_k)`$ avec
$`Y_k=(\tilde\Psi_h-e^{h\mathcal L})(e^{kh\mathcal L}\rho)`$, qui est autoadjoint (les applications CP préservent
l'autoadjonction, par la forme de Kraus) avec $`\|Y_k\|_1\le Ch^3`$ d'après (2), car l'état $`e^{kh\mathcal L}\rho`$ est de norme
de trace $`1`$. D'après (1) et l'étape 0, $`\|\tilde\Psi_h^{\,n-1-k}(Y_k)\|_1\le Ch^3`$, et les $`n`$ termes ont une somme
au plus égale à $`nCh^3=Cth^2`$. $`\square`$

C'est à l'étape 4 que le fait que les pas soient des canaux paie : les erreurs commises à chaque pas ne sont jamais
amplifiées, quel que soit $`n`$, et aucune condition de stabilité sur $`h`$ n'est nécessaire. En développant la définition,
les opérateurs de Kraus de $`\tilde\Psi_h`$ sont $`E_h^2N_h`$, $`\sqrt h\,E_hJE_hN_h`$ et $`\frac h{\sqrt2}E_hJ'JE_hN_h`$, où
$`J,J'`$ parcourent les opérateurs de Kraus de $`\Phi_{\mathrm{tot}}`$ : un demi-pas de l'évolution sans saut de la
section 8.4, au plus deux sauts, l'autre demi-pas, et la correction $`N_h`$ des probabilités, avec $`\|N_h-1\|\le C_1h^3`$
d'après l'étape 3.

*Dans la bibliothèque :* $`\tilde\Psi_h`$ est `lindblad_evolve(..., method='kraus')` ; `method='exact'` calcule
$`e^{h\mathcal L}`$ lui-même à partir de la matrice de $`\mathcal L`$ sur $`\mathcal M`$ (de taille $`\sum_ck_c^2`$, sauts entre
secteurs compris), dont les pas sont des canaux d'après la proposition 8.8(1) ; le `method='rk4'` par défaut est
$`R_h`$ — [tutoriel, section 8.7](../../tutorial/fr/08_dynamics.md#s-8-7).

<a id="x-8"></a>

## Exercices

1. (Représentation de Heisenberg) Montrer que
   $`\mathrm{Tr}(\mathcal L(\rho)A)=\mathrm{Tr}(\rho\,\mathcal L^*(A))`$ pour le
   générateur de la section 8.3, avec
   $`\mathcal L^*(A)=i[H,A]+\sum_jg_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})+\sum_mg_m(\Phi_m^*(A)-\frac12\{\Phi_m^*(1),A\})`$,
   que $`\mathcal L^*(1)=0`$, et que $`\mathcal L^*(P_c)=0`$ s'il n'y a pas de sauts entre secteurs.
2. (Point exceptionnel) Dans le secteur $`N=1`$ de la section 8.4,
   $`H_{\mathrm{eff}}=\begin{pmatrix}-i\gamma_1/2&J\\J&-i\gamma_2/2\end{pmatrix}`$. Trouver ses
   valeurs propres et montrer que, pour $`J=|\gamma_1-\gamma_2|/4>0`$, il n'est pas diagonalisable.
3. On suppose que $`N_0`$ bosons répartis sur $`m`$ modes sont perdus chacun au taux $`\gamma`$ (sauts
   $`a_1,\dots,a_m`$ restreints à $`N\to N-1`$). Montrer que $`p_N(t)=\binom{N_0}Nq^N(1-q)^{N_0-N}`$ avec
   $`q=e^{-\gamma t}`$, pour tout hamiltonien de $`\mathcal M`$.
4. (Déphasage pur) Sur $`M_2`$, soient $`H=\omega|e\rangle\langle e|`$ et
   $`L=\sigma_z=\mathrm{diag}(1,-1)`$ de taux $`\gamma`$. Montrer que les populations sont
   constantes, que $`\rho_{ge}(t)=e^{(i\omega-2\gamma)t}\rho_{ge}(0)`$, et que les états stationnaires
   engendrent l'espace de dimension deux des matrices diagonales, bien que $`C=1`$.

<details>
<summary>Solutions</summary>

1. Par cyclicité, $`\mathrm{Tr}(L\rho L^*A)=\mathrm{Tr}(\rho L^*AL)`$,
   $`\mathrm{Tr}(\{L^*L,\rho\}A)=\mathrm{Tr}(\rho\{L^*L,A\})`$ et
   $`\mathrm{Tr}(-i[H,\rho]A)=\mathrm{Tr}(\rho\,i[H,A])`$ ; pour les sauts entre secteurs,
   $`\mathrm{Tr}(\Phi(\rho)A)=\mathrm{Tr}(\rho\Phi^*(A))`$ est la définition de $`\Phi^*`$
   (section 6.3). Alors
   $`\mathcal L^*(1)=\sum_jg_j(L_j^*L_j-L_j^*L_j)+\sum_mg_m(\Phi_m^*(1)-\Phi_m^*(1))=0`$, ce qui est
   encore une fois la préservation de la trace. Sans sauts entre secteurs, $`P_c`$ commute avec $`H`$ et
   avec tous les $`L_j`$, donc $`L_j^*P_cL_j=L_j^*L_jP_c=\frac12\{L_j^*L_j,P_c\}`$, chaque $`\Phi_m`$ n'a que les blocs $`K^{cc}_{m,l}`$, donc
   $`\Phi_m^*(P_c)=\Phi_m^*(1)P_c=\frac12\{\Phi_m^*(1),P_c\}`$, et $`\mathcal L^*(P_c)=0`$ ; alors $`p_c'=\mathrm{Tr}(\rho\mathcal L^*(P_c))=0`$, ce qui fournit
   une seconde démonstration de la proposition 8.9(1).

2. $`H_{\mathrm{eff}}=-\frac i4(\gamma_1+\gamma_2)1+\begin{pmatrix}-i\delta&J\\J&i\delta\end{pmatrix}`$
   avec $`\delta=(\gamma_1-\gamma_2)/4`$ ; la seconde matrice est de trace $`0`$ et de déterminant
   $`\delta^2-J^2`$, donc les valeurs propres sont
   $`-\frac i4(\gamma_1+\gamma_2)\pm\sqrt{J^2-\delta^2}`$ : deux taux de décroissance pour
   $`J<|\delta|`$, un taux et deux fréquences pour $`J>|\delta|`$. Pour $`J=|\delta|>0`$, la seconde
   matrice $`Y`$ vérifie $`Y^2=(J^2-\delta^2)1=0`$ mais $`Y\ne0`$, donc $`H_{\mathrm{eff}}`$ a pour unique
   valeur propre $`-\frac i4(\gamma_1+\gamma_2)`$ sans être un multiple de $`1`$ : il n'est pas
   diagonalisable.

3. Sur le secteur $`N`$, $`\sum_ja_j^*a_j=N1`$, et la seule paire de secteurs est $`(N-1,N)`$, donc la
   proposition 8.9(2) donne $`p_N'=-\gamma Np_N+\gamma(N+1)p_{N+1}`$. Avec $`q'=-\gamma q`$, la dérivée
   du terme binomial est $`-\gamma Np_N+\gamma(N_0-N)\binom{N_0}Nq^{N+1}(1-q)^{N_0-N-1}`$, et
   $`(N+1)\binom{N_0}{N+1}=(N_0-N)\binom{N_0}N`$ montre que le second terme vaut $`\gamma(N+1)p_{N+1}`$.
   En $`t=0`$, c'est la masse de Dirac en $`N_0`$ ; l'unicité est le lemme 8.1.

4. $`L\rho L^*-\frac12\{L^*L,\rho\}=\sigma_z\rho\sigma_z-\rho`$ multiplie les coefficients hors
   diagonale par $`-2`$ et annule la diagonale, et $`-i[H,\rho]_{ge}=i\omega\rho_{ge}`$. Ainsi
   $`\mathcal L(\sigma)=0`$ si et seulement si $`\sigma_{ge}=0`$ (pour $`\gamma>0`$), et le spectre est
   $`\{0,0,-2\gamma\pm i\omega\}`$. Avec $`U_t=e^{-itH}`$,
   $`e^{t\mathcal L}(\rho)=\frac{1+e^{-2\gamma t}}2U_t\rho U_t^*+\frac{1-e^{-2\gamma t}}2U_t\sigma_z\rho\sigma_zU_t^*`$ :
   comme $`\sigma_z`$ commute avec $`U_t`$ et $`\sigma_z^2=1`$, la dérivation montre que le membre de droite
   $`R(t)`$ vérifie $`R'=-i[H,R]+\gamma(\sigma_zR\sigma_z-R)=\mathcal L(R)`$ et $`R(0)=\rho`$, donc le lemme
   8.1 s'applique. Pour $`t>0`$, c'est un canal unifère dont les points fixes sont les matrices qui
   commutent avec ses opérateurs de Kraus (théorème 6.10), c'est-à-dire les matrices diagonales.

</details>

*Calcul :* les exercices 2 et 4, numériquement — [tutoriel, Exercices](../../tutorial/fr/08_dynamics.md#x-8).

**Suite :** le [chapitre 9](09_spin_chains.md) traite des chaînes de spins : secteurs d'aimantation et
d'impulsion, évolution temporelle par les méthodes de Lanczos et de Krylov pour les grands secteurs,
et désordre.
