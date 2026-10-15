[English](../C_eigenvalue_densities.md) · [Español](../es/C_eigenvalue_densities.md) · **Français** · [Deutsch](../de/C_eigenvalue_densities.md) · [中文](../zh/C_eigenvalue_densities.md) · [日本語](../ja/C_eigenvalue_densities.md) · [Русский](../ru/C_eigenvalue_densities.md)

# Annexe C. Densités de valeurs propres et formule de Page

Le chapitre 5 utilise la densité jointe des valeurs propres des ensembles gaussiens, et le chapitre 7 utilise
la formule de Page pour l'entropie d'intrication moyenne d'un état pur aléatoire. Cette annexe démontre les deux.
Le chemin est le suivant : la formule d'intégration de Weyl (le jacobien de la « diagonalisation »), les ensembles
gaussiens, la fonction à un point d'une densité déterminantale par les polynômes orthogonaux, l'ensemble de
Wishart–Laguerre et la mesure induite sur les états réduits, et enfin un calcul exact avec les polynômes de
Laguerre.

**Vous apprendrez** pourquoi les valeurs propres d'une matrice aléatoire invariante se repoussent comme
$`\prod_{i<j}|\lambda_i-\lambda_j|^\beta`$, pourquoi les vecteurs propres sont alors distribués selon la loi de Haar et
indépendants des valeurs propres, comment s'intègre un déterminant de polynômes orthonormés, et comment
$`\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}`$ découle d'une poignée d'intégrales.

Outre les prérequis du cours (théorie de la mesure, la formule du changement de variables dans
$`\mathbb R^N`$ pour les applications $`C^1`$ injectives de différentielle inversible, les théorèmes d'inversion locale et des
fonctions implicites), nous utilisons du chapitre 5 le lemme 5.1 (unicité et invariance à droite de la mesure de Haar),
le lemme 5.2 et le théorème 5.3 (existence de la mesure de Haar sur $`U(n)`$ par la QR gaussienne), de l'annexe B
le théorème B.17 (la même chose pour $`O(n)`$), de l'annexe A la proposition A.16(a) (invariance unitaire des vecteurs
gaussiens complexes) et le théorème A.20(b) (la loi de $`|z|^2`$), et du chapitre 7 la
représentation $`\rho_A=GG^*/\mathrm{Tr}\,GG^*`$ de la démonstration du théorème 7.11.

**Notations.** $`\beta\in\{1,2\}`$ et, en conséquence, $`\mathbb F=\mathbb R`$ ou $`\mathbb C`$. $`\mathcal H`$ est
l'espace des matrices $`n\times n`$ réelles symétriques ($`\beta=1`$) ou complexes hermitiennes ($`\beta=2`$), muni de la
mesure de Lebesgue $`dH`$ dans les coordonnées $`h_{ii}`$, $`\mathrm{Re}\,h_{ij}`$ et (pour $`\beta=2`$)
$`\mathrm{Im}\,h_{ij}`$, $`i<j`$ ; sa dimension est $`n+\beta n(n-1)/2`$. $`K=O(n)`$ ou $`U(n)`$, $`T`$ est le
groupe des éléments diagonaux de $`K`$ (signes, ou phases), $`\mathfrak k=\{X\in M_n(\mathbb F):X^*=-X\}`$,
$`\mathfrak t`$ sa partie diagonale et $`\mathfrak m`$ sa partie hors diagonale (diagonale nulle), avec les
coordonnées $`\mathrm{Re}\,X_{ij}`$, $`\mathrm{Im}\,X_{ij}`$, $`i<j`$. Pour $`\lambda\in\mathbb R^n`$,
$`\Lambda=\mathrm{diag}(\lambda)`$ et $`\Delta(\lambda)=\prod_{i<j}(\lambda_j-\lambda_i)`$ ; $`C=\{\lambda:\lambda_1<\dots<\lambda_n\}`$,
sur lequel $`\Delta>0`$. Pour $`H\in\mathcal H`$, $`\lambda(H)`$ est le vecteur de ses valeurs propres rangées par ordre
croissant. $`\|X\|_2^2=\mathrm{Tr}\,X^*X`$. Pour $`V\in K`$ l'application $`H\mapsto VHV^*`$ est linéaire sur $`\mathcal H`$
et conserve $`\|\cdot\|_2`$ ; une isométrie linéaire vérifie $`|\det|=1`$, donc elle **conserve $`dH`$**.

*Calcul :* les vérifications numériques de cette annexe — [tutoriel, annexe C](../../tutorial/fr/C_eigenvalue_densities.md).

<a id="s-C-1"></a>

## C.1 Spectres simples et repères propres

<a id="r-C-1"></a>

**Lemme C.1.** (a) L'ensemble $`\Omega\subseteq\mathcal H`$ des matrices à $`n`$ valeurs propres distinctes est
ouvert, et $`\mathcal H\setminus\Omega`$ est Lebesgue-négligeable. (b) Sur $`\Omega`$ l'application $`H\mapsto\lambda(H)`$ est
$`C^\infty`$, et les projecteurs spectraux $`P_j(H)=\prod_{k\neq j}\frac{H-\lambda_k(H)}{\lambda_j(H)-\lambda_k(H)}`$ sont continus.

*Démonstration.* (a) Soit $`A_{ij}=\lambda_j^{i-1}`$ avec $`\lambda=\lambda(H)`$. Alors $`\det A=\Delta(\lambda)`$ (Vandermonde) et
$`(AA^T)_{ij}=\sum_k\lambda_k^{i+j-2}=\mathrm{Tr}\,H^{i+j-2}`$, donc
$`\mathrm{disc}(H)=\det\big(\mathrm{Tr}\,H^{i+j-2}\big)_{i,j=1}^n=\Delta(\lambda(H))^2`$ est un polynôme
en les coordonnées de $`H`$, non nul en $`\mathrm{diag}(1,\dots,n)`$, et $`\Omega=\{\mathrm{disc}\,\neq0\}`$ est ouvert.
L'ensemble des zéros d'un polynôme non nul $`p`$ sur $`\mathbb R^N`$ est négligeable, par récurrence sur $`N`$ : pour $`N=1`$ il est
fini ; sinon on écrit $`p=\sum_jp_j(x')x_N^j`$ avec un $`p_{j_0}\neq0`$ ; pour $`x'`$ hors de l'ensemble des zéros de
$`p_{j_0}`$ (négligeable par récurrence) $`x_N\mapsto p(x',x_N)`$ a un nombre fini de zéros, et le théorème de Tonelli appliqué
à l'indicatrice de l'ensemble des zéros donne la mesure $`0`$.
(b) $`\chi(H,t)=\det(t-H)=\prod_i(t-\lambda_i)`$ est un polynôme réel en $`(H,t)`$. En $`H_0\in\Omega`$,
$`\partial_t\chi(H_0,\lambda_j(H_0))=\prod_{k\neq j}(\lambda_j-\lambda_k)\neq0`$, donc le théorème des fonctions implicites fournit des
fonctions $`C^\infty`$ $`r_j`$ près de $`H_0`$ avec $`\chi(H,r_j(H))=0`$ et $`r_j(H_0)=\lambda_j(H_0)`$. Par continuité
$`r_1<\dots<r_n`$ près de $`H_0`$ : ce sont $`n`$ racines distinctes d'un polynôme de degré $`n`$, donc toutes ses racines, dans
l'ordre, et $`r_j=\lambda_j`$. Si $`H=\sum_k\lambda_kQ_k`$ est la décomposition spectrale, l'interpolation de Lagrange donne
$`P_j(H)=Q_j`$ ; la continuité découle de celle de $`\lambda`$. $`\square`$

Pour $`H\in\Omega`$, un **repère propre** est un $`V\in K`$ tel que $`H=V\Lambda V^*`$, $`\Lambda=\mathrm{diag}\,\lambda(H)`$. Si $`V,V'`$
sont des repères propres, $`M=V'^*V`$ commute avec $`\Lambda`$, donc $`(\lambda_i-\lambda_j)M_{ij}=0`$ et $`M\in T`$ ; réciproquement $`VD`$
est un repère propre pour tout $`D\in T`$. On obtient un **repère propre borélien** $`F(H)`$ en prenant pour $`j`$-ième colonne
$`P_j(H)e_m/\|P_j(H)e_m\|`$, où $`m`$ est le plus petit indice tel que $`P_j(H)e_m\neq0`$ : il est continu sur chacun des
boréliens où ces indices sont fixés, ses colonnes sont des vecteurs propres orthonormés, et il est réel lorsque $`H`$ l'est.

<a id="r-C-2"></a>

**Lemme C.2 (mesures de Haar).** $`U(n)`$, $`O(n)`$ et $`T`$ possèdent des mesures de Haar, et $`\int_Kf(VW)\,dV=\int_Kf\,dV`$
pour tout $`W\in K`$ et toute $`f`$ borélienne bornée.

*Démonstration.* $`U(n)`$ : théorème 5.3. $`O(n)`$ : théorème B.17 (la construction QR gaussienne réelle : une matrice gaussienne réelle est
presque sûrement inversible, et le facteur orthogonal de sa décomposition QR est distribué selon la loi de Haar). $`T`$ : la mesure uniforme sur $`\{\pm1\}^n`$, ou le produit des longueurs d'arc normalisées sur les
cercles, est invariante par multiplication par un élément fixé. L'invariance à droite est le lemme 5.1. $`\square`$

<a id="r-C-3"></a>

**Lemme C.3 (repères propres des ensembles invariants).** Soit $`H`$ un élément aléatoire de $`\mathcal H`$ avec
$`P(H\in\Omega)=1`$ dont la loi est invariante par $`H\mapsto VHV^*`$ pour tout $`V\in K`$. Soient $`D`$ uniforme sur $`T`$ et
indépendant de $`H`$, et $`\tilde U=F(H)D`$. Alors $`\tilde U`$ est un repère propre de $`H`$, il est distribué selon la loi de Haar sur
$`K`$, et il est indépendant de $`\lambda(H)`$.

*Démonstration.* Fixons $`V\in K`$. $`F(VHV^*)`$ et $`VF(H)`$ sont tous deux des repères propres de $`VHV^*`$, donc
$`F(VHV^*)=VF(H)E(H)`$ avec $`E(H)=(VF(H))^*F(VHV^*)\in T`$, fonction borélienne de $`H`$. Pour des fonctions boréliennes bornées
$`\varphi`$ sur $`\mathbb R^n`$ et $`\psi`$ sur $`K`$, comme $`VHV^*`$ a la loi de $`H`$ et $`\lambda(VHV^*)=\lambda(H)`$,
```math
\mathbb E\big[\varphi(\lambda(H))\psi(\tilde U)\big]=\mathbb E\big[\varphi(\lambda(H))\psi(VF(H)E(H)D)\big]
=\mathbb E\big[\varphi(\lambda(H))\psi(V\tilde U)\big],
```
où la deuxième étape intègre d'abord en $`D`$ (Fubini) : pour $`H`$ fixé, $`E(H)D`$ est uniforme sur $`T`$. Soit maintenant $`V`$
distribué selon la loi de Haar et indépendant de tout le reste ; par Fubini et le lemme C.2,
$`\mathbb E[\varphi(\lambda(H))\psi(\tilde U)]=\mathbb E\big[\varphi(\lambda(H))\int_K\psi(V\tilde U)\,dV\big] =\mathbb E[\varphi(\lambda(H))]\int_K\psi\,dV`$. Avec $`\varphi=1`$, $`\tilde U`$ suit la loi de Haar ; la forme produit est
l'indépendance. $`\square`$

<a id="s-C-2"></a>

## C.2 La formule d'intégration de Weyl

Pour $`\varepsilon>0`$ soient $`\mathfrak m_\varepsilon=\{X\in\mathfrak m:\|X\|_2<\varepsilon\}`$ et
```math
\Psi(X,\lambda)=e^X\Lambda e^{-X},\qquad (X,\lambda)\in\mathfrak m_\varepsilon\times C .
```
Pour $`X\in\mathfrak k`$, $`(e^X)^*=e^{-X}=(e^X)^{-1}`$, donc $`e^X\in K`$.

<a id="r-C-4"></a>

**Lemme C.4 (une carte).** Il existe $`\varepsilon>0`$ et une fonction continue $`j:\overline{\mathfrak m_\varepsilon}\to(0,\infty)`$
avec $`j(0)=1`$ tels que : (a) $`|\det D\Psi(X,\lambda)|=\Delta(\lambda)^\beta j(X)`$ ; (b) $`\Psi`$ est injective ; (c) pour
$`H\in\Omega`$, $`H\in\Psi(\mathfrak m_\varepsilon\times C)`$ si et seulement si $`F(H)\in S=\{e^XD:X\in\mathfrak m_\varepsilon,D\in T\}`$, et
$`S`$ est un borélien vérifiant $`ST=S`$.

*Démonstration.* (a) $`\Psi`$ est $`C^\infty`$ (séries entières). Pour $`Y\in\mathfrak m`$ posons $`A_X(Y)=e^{-X}\frac{d}{dt}\big|_0e^{X+tY}`$.
En dérivant $`W(t)^*W(t)=1`$ pour $`W(t)=e^{X+tY}\in K`$ on obtient $`A_X(Y)\in\mathfrak k`$, et avec
$`\frac{d}{dt}W^{-1}=-W^{-1}W'W^{-1}`$, pour $`g=e^X`$,
```math
D\Psi(X,\lambda)(Y,\mu)=g\big([A_X(Y),\Lambda]+\mathrm{diag}(\mu)\big)g^{-1}.
```
Les matrices diagonales commutent avec $`\Lambda`$, donc $`[A_X(Y),\Lambda]=[B_XY,\Lambda]`$ avec $`B_X=\pi\circ A_X:\mathfrak m\to\mathfrak m`$ et
$`\pi`$ la projection de $`\mathfrak k`$ sur $`\mathfrak m`$ parallèlement à $`\mathfrak t`$. Ainsi $`D\Psi=\mathrm{Ad}_g\,\circ L_\lambda\circ(B_X\times1)`$
avec $`L_\lambda(Z,\mu)=[Z,\Lambda]+\mathrm{diag}(\mu)`$. Comme $`([Z,\Lambda])_{ij}=Z_{ij}(\lambda_j-\lambda_i)`$, dans les coordonnées
$`L_\lambda`$ multiplie $`\mathrm{Re}\,Z_{ij}`$ (et $`\mathrm{Im}\,Z_{ij}`$) par $`\lambda_j-\lambda_i`$ et envoie $`\mu`$ sur la diagonale :
$`|\det L_\lambda|=\Delta(\lambda)^\beta`$. $`\mathrm{Ad}_g`$ conserve $`dH`$. Posons $`j(X)=|\det B_X|`$ : elle est continue, et
$`B_0=1`$ car $`\frac d{dt}e^{tY}|_0=Y`$ ; on choisit $`\varepsilon`$ tel que $`j>0`$ sur $`\overline{\mathfrak m_\varepsilon}`$.
(b) Si $`\Psi(X,\lambda)=\Psi(X',\lambda')=H`$, alors $`\lambda=\lambda'=\lambda(H)`$, et $`e^X,e^{X'}`$ sont des repères propres de $`H\in\Omega`$,
donc $`e^X=e^{X'}D`$ avec $`D\in T`$. Soit $`G(X,t)=e^Xe^t`$ sur $`\mathfrak m\times\mathfrak t`$ (pour $`\beta=1`$, $`\mathfrak t=\{0\}`$) ; $`DG(0,0)(Y,s)=Y+s`$ est injective,
et comme $`Y`$ est hors diagonale et $`s`$ diagonale, $`\|Y+s\|_2^2=\|Y\|_2^2+\|s\|_2^2`$, donc $`\|DG(0,0)v\|=\|v\|`$ pour la norme produit ; posons $`m=1`$. On choisit $`\delta>0`$ tel que $`\|DG(z)-DG(0)\|\le m/2`$ pour $`\|z\|<\delta`$ ; pour $`u,v`$ dans cette
boule, $`G(u)-G(v)=\int_0^1DG(v+s(u-v))(u-v)\,ds`$ donne $`\|G(u)-G(v)\|\ge\frac m2\|u-v\|`$, donc $`G`$ y est injective.
($`\delta`$ ne dépend que de $`n,\beta`$ ; $`\varepsilon`$ est choisi ensuite.)
Or $`\|D-1\|_2=\|e^{X'}(D-1)\|_2=\|e^X-e^{X'}\|_2\le\|e^X-1\|_2+\|e^{X'}-1\|_2\le2(e^\varepsilon-1)`$, par la série entière et
$`\|AB\|_2\le\|A\|_2\|B\|_2`$. Pour $`\beta=1`$, $`D=\mathrm{diag}(\pm1)`$ et un coefficient $`-1`$ donnerait
$`\|D-1\|_2\ge2`$ ; pour $`\varepsilon<\log2`$, $`D=1=e^t`$ avec $`t=0`$. Pour $`\beta=2`$, $`D=e^t`$ avec $`t=i\mathrm{diag}(\theta)`$,
$`\theta_j\in(-\pi,\pi]`$, et $`|\theta_j|\le\frac\pi2|e^{i\theta_j}-1|`$ (car $`|e^{i\theta}-1|=2|\sin\frac\theta2|`$ et $`\sin x\ge2x/\pi`$ sur $`[0,\pi/2]`$), donc
$`\|t\|_2\le\frac\pi2\|D-1\|_2\le\pi(e^\varepsilon-1)`$. Pour $`\varepsilon`$ petit ($`\varepsilon^2+\pi^2(e^\varepsilon-1)^2<\delta^2`$), $`(X,0)`$ et $`(X',t)`$ sont dans la boule de rayon $`\delta`$ et
$`G(X,0)=e^X=e^{X'}D=G(X',t)`$, donc $`X=X'`$.
(c) Si $`H=\Psi(X,\lambda')`$, alors $`e^X`$ est un repère propre, donc $`F(H)\in e^XT\subseteq S`$. Si $`F(H)=e^XD`$, alors
$`H=F(H)\Lambda F(H)^*=\Psi(X,\lambda(H))`$. $`S=G(\mathfrak m_\varepsilon\times\mathfrak t)`$ (ou $`e^{\mathfrak m_\varepsilon}T`$, $`T`$ fini) est une
image continue d'un ensemble $`\sigma`$-compact, donc $`\sigma`$-compact et borélien ; $`ST=S`$ est clair. $`\square`$

*Calcul :* le jacobien de $`\Psi`$ pour $`n=2,3`$ et $`\beta=1,2`$ — [tutoriel, section C.2](../../tutorial/fr/C_eigenvalue_densities.md#s-C-2).

<a id="r-C-5"></a>

**Théorème C.5 (formule d'intégration de Weyl).** Il existe une constante $`c_{n,\beta}>0`$ telle que pour toute fonction borélienne
$`f:\mathcal H\to[0,\infty]`$,
```math
\int_{\mathcal H}f(H)\,dH=c_{n,\beta}\int_C\Delta(\lambda)^\beta\Big(\int_Kf(V\Lambda V^*)\,dV\Big)d\lambda .
```
De manière équivalente, $`\Phi(V,\lambda)=V\Lambda V^*`$ envoie la mesure $`c_{n,\beta}\,dV\otimes\Delta(\lambda)^\beta d\lambda`$ sur $`K\times C`$
sur la mesure de Lebesgue de $`\mathcal H`$. Le facteur $`\Delta(\lambda)^\beta`$ est le jacobien de la carte $`\Psi`$ en $`X=0`$
(lemme C.4(a)) ; la constante $`c_{n,\beta}`$ reflète la normalisation de la mesure de Haar et n'est jamais nécessaire
explicitement dans la suite, puisque toute densité obtenue par la formule est normalisée par sa masse totale.

*Démonstration.* *Ensemble de référence.* Soit $`H_0`$ de densité $`p_0(H)=Z_0^{-1}e^{-\mathrm{Tr}\,H^2/2}`$ ; comme
$`\mathrm{Tr}\,H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2`$, $`Z_0`$ est un produit fini d'intégrales gaussiennes. $`p_0`$ et $`dH`$
sont invariants par $`H\mapsto VHV^*`$, donc la loi de $`H_0`$ l'est aussi ; $`P(H_0\in\Omega)=1`$ par le lemme C.1. Par le lemme C.3,
$`\tilde U`$ suit la loi de Haar, est indépendant de $`\lambda=\lambda(H_0)`$, et $`H_0=\tilde U\Lambda\tilde U^*`$.
*La loi de $`\lambda`$.* Soit $`A\subseteq C`$ borélien. Par le lemme C.4(c) et $`ST=S`$,
$`\{H_0\in\Psi(\mathfrak m_\varepsilon\times A)\}=\{\lambda\in A,\ \tilde U\in S\}`$ à un événement négligeable près, dont la probabilité est
$`P(\lambda\in A)\,h`$ avec $`h=\int_K1_S\,dV`$. D'autre part $`\Psi`$ est injective et $`C^1`$ de jacobien non nul sur
l'ouvert $`\mathfrak m_\varepsilon\times C`$, et $`p_0(\Psi(X,\lambda))=p_0(\Lambda)`$ ; la formule du changement de variables et le lemme C.4(a) donnent
```math
P\big(H_0\in\Psi(\mathfrak m_\varepsilon\times A)\big)=\int_{\mathfrak m_\varepsilon}\!\int_Ap_0(\Lambda)\Delta(\lambda)^\beta j(X)\,d\lambda\,dX
=\kappa\int_Ap_0(\Lambda)\Delta(\lambda)^\beta\,d\lambda,\qquad \kappa=\int_{\mathfrak m_\varepsilon}j>0 .
```
Avec $`A=C`$, $`h=\kappa\int_Cp_0\Delta^\beta>0`$. Donc $`\lambda`$ a pour densité $`c\,p_0(\Lambda)\Delta(\lambda)^\beta`$ sur $`C`$, $`c=\kappa/h`$.
*$`f`$ quelconque.* Comme $`p_0>0`$, $`\int f\,dH=\mathbb E\big[(f/p_0)(\tilde U\Lambda\tilde U^*)\big]`$. La loi jointe de $`(\tilde U,\lambda)`$ est
$`dV\otimes c\,p_0(\Lambda)\Delta^\beta d\lambda`$, et $`p_0(V\Lambda V^*)=p_0(\Lambda)`$ ; le théorème de Tonelli donne la formule avec $`c_{n,\beta}=c`$.
La constante est déterminée par la formule elle-même (prendre $`f=p_0`$). $`\square`$

<a id="r-C-6"></a>

**Corollaire C.6 (densités invariantes).** Soit $`H`$ de densité $`p`$ sur $`\mathcal H`$ avec $`p(VHV^*)=p(H)`$ pour tout
$`V\in K`$. Alors $`\lambda(H)`$ a pour densité $`c_{n,\beta}\,p(\Lambda)\Delta(\lambda)^\beta`$ sur $`C`$ ; le repère propre $`\tilde U`$ du lemme C.3
est distribué selon la loi de Haar et indépendant de $`\lambda(H)`$ ; et pour toute fonction borélienne symétrique $`g\ge0`$ sur $`\mathbb R^n`$,
```math
\mathbb E\,g(\lambda(H))=\frac{c_{n,\beta}}{n!}\int_{\mathbb R^n}g(x)\,p(\mathrm{diag}\,x)\,|\Delta(x)|^\beta\,dx .
```

*Démonstration.* $`P(H\in\Omega)=1`$ par le lemme C.1. Pour $`A\subseteq C`$ borélien, on applique le théorème C.5 à $`f=1_A(\lambda(\cdot))\,p`$, pour laquelle
$`f(V\Lambda V^*)=1_A(\lambda)p(\Lambda)`$. La loi de $`H`$ est invariante, donc le lemme C.3 s'applique. À l'ensemble négligeable près où deux
coordonnées coïncident, $`\mathbb R^n`$ est la réunion disjointe des $`n!`$ images de $`C`$ par les permutations des coordonnées, et
$`g`$, $`|\Delta|`$ et $`x\mapsto p(\mathrm{diag}\,x)`$ sont invariantes par permutation (une matrice de permutation appartient à $`K`$). $`\square`$

<a id="s-C-3"></a>

## C.3 Les ensembles gaussiens

<a id="r-C-7"></a>

**Théorème C.7 (GOE et GUE).** Soit $`H=(G+G^T)/2`$ avec $`G`$ réelle à coefficients indépendants de loi $`N(0,1)`$ ($`\beta=1`$, GOE), ou
$`H=(G+G^*)/\sqrt2`$ avec $`G`$ de Ginibre ($`\beta=2`$, GUE), comme dans la section 5.3. Alors $`H`$ a pour densité
$`(2\pi)^{-n/2}\pi^{-\beta n(n-1)/4}e^{-\mathrm{Tr}\,H^2/2}`$ sur $`\mathcal H`$, ses valeurs propres ont pour densité
```math
Z_{n,\beta}^{-1}\,\Delta(\lambda)^\beta\,e^{-\sum_i\lambda_i^2/2}\quad\text{on }C,\qquad Z_{n,\beta}=\int_C\Delta^\beta e^{-|\lambda|^2/2}\,d\lambda\in(0,\infty),
```
et le repère propre $`\tilde U`$ est distribué selon la loi de Haar sur $`K`$ et indépendant des valeurs propres.

*Démonstration.* GOE : $`h_{ii}=G_{ii}\sim N(0,1)`$ et $`h_{ij}=(G_{ij}+G_{ji})/2\sim N(0,\frac12)`$ pour $`i<j`$. GUE : $`h_{ii}=\sqrt2\mathrm{Re}\,G_{ii}\sim N(0,1)`$,
$`\mathrm{Re}\,h_{ij}=(\mathrm{Re}\,G_{ij}+\mathrm{Re}\,G_{ji})/\sqrt2`$ et $`\mathrm{Im}\,h_{ij}=(\mathrm{Im}\,G_{ij}-\mathrm{Im}\,G_{ji})/\sqrt2`$ suivent la loi $`N(0,\frac12)`$.
Des coordonnées distinctes sont fonctions de familles disjointes de variables indépendantes, donc elles sont indépendantes, et la
densité est $`\prod_i(2\pi)^{-1/2}e^{-h_{ii}^2/2}\prod_{i<j}\pi^{-\beta/2}e^{-|h_{ij}|^2}`$, qui est celle de l'énoncé car
$`\mathrm{Tr}\,H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2`$. Elle est invariante, et le corollaire C.6 s'applique ; la constante est fixée par
la normalisation. $`Z_{n,\beta}<\infty`$ car $`|\lambda_j-\lambda_i|\le(1+|\lambda_i|)(1+|\lambda_j|)`$, donc
$`|\Delta(\lambda)|\le\prod_i(1+|\lambda_i|)^{n-1}`$. $`\square`$

<a id="r-C-8"></a>

**Corollaire C.8 (densité de niveaux).** Pour le GOE ou le GUE de taille $`n`$, soient $`q(x)=(n!Z_{n,\beta})^{-1}|\Delta(x)|^\beta e^{-|x|^2/2}`$ sur
$`\mathbb R^n`$ et $`\rho(x)=n\int_{\mathbb R^{n-1}}q(x,y)\,dy`$. Alors $`\mathbb E\sum_i\varphi(\lambda_i(H))=\int\varphi\rho`$ pour toute fonction borélienne bornée
$`\varphi`$, $`\int\rho=n`$, et $`\rho`$ est continue avec $`0\le\rho(x)\le C_n(1+|x|)^{\beta(n-1)}e^{-x^2/2}`$. En particulier
l'hypothèse du théorème 5.15 est satisfaite pour des blocs GOE ou GUE indépendants.

*Démonstration.* Par le corollaire C.6 avec $`g(x)=\sum_i\varphi(x_i)`$ (pour $`\varphi\ge0`$ ; en général on l'applique à $`\varphi^\pm=\max(\pm\varphi,0)`$) et la symétrie de $`q`$,
$`\mathbb E\sum_i\varphi(\lambda_i)=n\int\varphi(x_1)q(x)\,dx=\int\varphi\rho`$ (Fubini). Par la majoration de $`|\Delta|`$,
$`q(x,y)\le(n!Z_{n,\beta})^{-1}(1+|x|)^{\beta(n-1)}e^{-x^2/2}\prod_{i}(1+|y_i|)^{\beta(n-1)}e^{-y_i^2/2}`$, ce qui donne la majoration après
intégration en $`y`$ et, puisque $`\sup_x(1+|x|)^{\beta(n-1)}e^{-x^2/2}<\infty`$, un majorant intégrable indépendant de $`x`$ : la continuité
découle de la convergence dominée. $`\square`$

*Calcul :* le théorème C.7 et le corollaire C.8 pour le GOE avec $`n=3`$ — [tutoriel, section C.3](../../tutorial/fr/C_eigenvalue_densities.md#s-C-3).

<a id="r-C-9"></a>

**Proposition C.9 (plus petit espacement).** Pour le GOE ou le GUE de taille $`n\ge2`$ soit $`\delta=\min_i(\lambda_{i+1}-\lambda_i)`$. Il existe
des constantes $`0<c\le C`$ (dépendant de $`n,\beta`$) telles que $`c\,\varepsilon^{\beta+1}\le P(\delta<\varepsilon)\le C\varepsilon^{\beta+1}`$ pour $`0<\varepsilon\le1`$.

*Démonstration.* *Majoration.* $`P(\delta<\varepsilon)\le\sum_iP(\lambda_{i+1}-\lambda_i<\varepsilon)`$. Sur $`\{0<u=x_{i+1}-x_i<\varepsilon\}`$ on majore le facteur de
la paire $`(i,i+1)`$ par $`u^\beta`$, les autres comme ci-dessus, $`1+|x_{i+1}|\le2(1+|x_i|)`$ et $`e^{-x_{i+1}^2/2}\le1`$ ; en substituant
$`x_{i+1}=x_i+u`$,
```math
P(\lambda_{i+1}-\lambda_i<\varepsilon)\le\frac{2^{\beta(n-1)}}{Z_{n,\beta}}\int_0^\varepsilon u^\beta du\int(1+|x|)^{2\beta(n-1)}e^{-\frac{x^2}2}dx
\Big(\int(1+|x|)^{\beta(n-1)}e^{-\frac{x^2}2}dx\Big)^{n-2}.
```
*Minoration.* Soit $`R_\varepsilon=\{x:0<x_1<1,\ 0<x_2-x_1<\varepsilon,\ 3k<x_k<3k+1\ (3\le k\le n)\}\subseteq C`$. Sur $`R_\varepsilon`$,
$`\delta\le x_2-x_1`$, toute paire autre que $`(1,2)`$ est à distance $`\ge1`$, et $`|x|^2\le M`$ avec $`M`$ indépendant de $`\varepsilon\le1`$ ; donc
$`P(\delta<\varepsilon)\ge P(\lambda\in R_\varepsilon)\ge Z_{n,\beta}^{-1}e^{-M/2}\int_0^\varepsilon u^\beta du`$. $`\square`$

Pour le GOE, c'est l'ordre $`\varepsilon^2`$ cité à la fin de la section 5.4. Pour $`U\in O(3)`$ de Haar, $`\mathbb E|U_{11}|^2=\frac13`$ : la multiplication à gauche par une matrice de permutation permute les
coefficients de la première colonne et conserve la mesure de Haar, donc les trois $`\mathbb E|U_{i1}|^2`$ sont égaux, et leur somme vaut $`1`$.

<a id="s-C-4"></a>

## C.4 Polynômes orthogonaux et fonction à un point

<a id="r-C-10"></a>

**Lemme C.10 (fonction à un point).** Soient $`I\subseteq\mathbb R`$ un intervalle et $`w\ge0`$ une fonction borélienne sur $`I`$, strictement positive sur un
ensemble de mesure non nulle, avec $`\int_I|x|^mw<\infty`$ pour tout $`m`$. Soient $`p_0,\dots,p_{n-1}`$ des polynômes réels, $`\deg p_k=k`$
de coefficient dominant $`\kappa_k\neq0`$, orthonormés dans $`L^2(w)`$ (ils existent : on applique Gram–Schmidt à $`1,x,\dots,x^{n-1}`$, qui sont
indépendants dans $`L^2(w)`$ car un polynôme non nul n'a qu'un nombre fini de zéros ; le signe de chaque $`p_k`$ peut être choisi arbitrairement). Soient $`q=Z_n^{-1}\Delta(x)^2\prod_iw(x_i)`$ sur
$`I^n`$ avec $`Z_n=\int_{I^n}\Delta^2\prod_iw(x_i)\,dx`$. Alors $`Z_n=n!\prod_k\kappa_k^{-2}`$ et
```math
n\int_{I^{n-1}}q(x,y)\,dy=K_n(x)\,w(x),\qquad K_n(x)=\sum_{k=0}^{n-1}p_k(x)^2=v(x)^TM^{-1}v(x),
```
où $`v(x)=(1,x,\dots,x^{n-1})`$ et $`M_{ij}=\int_Ix^{i+j}w`$, $`0\le i,j<n`$.

*Démonstration.* $`\Delta(x)=\det(x_j^{i-1})_{i,j}`$. Ajouter à chaque ligne des multiples des lignes situées au-dessus remplace $`x_j^{i-1}`$ par
$`p_{i-1}(x_j)/\kappa_{i-1}`$, donc $`\Delta(x)=(\prod_k\kappa_k)^{-1}\det(p_{i-1}(x_j))`$ et
```math
\Delta(x)^2=\Big(\prod_k\kappa_k\Big)^{-2}\sum_{\sigma,\tau\in S_n}\mathrm{sgn}\,\sigma\mathrm{sgn}\,\tau\prod_{j=1}^np_{\sigma(j)-1}(x_j)\,p_{\tau(j)-1}(x_j).
```
En intégrant en $`x_2,\dots,x_n`$ contre $`\prod_{j\ge2}w(x_j)`$, le facteur $`j`$ donne $`\delta_{\sigma(j)\tau(j)}`$, donc seul $`\sigma=\tau`$
subsiste, et le résultat est $`(\prod_k\kappa_k)^{-2}\sum_\sigma p_{\sigma(1)-1}(x_1)^2=(\prod_k\kappa_k)^{-2}(n-1)!\,K_n(x_1)`$.
En intégrant contre $`w(x_1)`$ on obtient $`Z_n=(\prod_k\kappa_k)^{-2}n!`$, et en divisant on obtient la formule. Enfin, si
$`p_k(x)=\sum_iC_{ki}x^i`$, l'orthonormalité s'écrit $`CMC^T=1`$ avec $`C`$ inversible, donc $`M^{-1}=C^TC`$ et
$`v^TM^{-1}v=|Cv|^2=\sum_kp_k^2`$. $`\square`$

Par le corollaire C.6, les valeurs propres non ordonnées du GUE ont la densité $`q`$ du lemme C.10 avec $`w(x)=e^{-x^2/2}`$ sur
$`\mathbb R`$, donc la densité de niveaux du GUE est $`K_n(x)e^{-x^2/2}`$.

*Calcul :* la densité de niveaux du GUE pour $`n=3`$ — [tutoriel, section C.4](../../tutorial/fr/C_eigenvalue_densities.md#s-C-4).

<a id="s-C-5"></a>

## C.5 Matrices de Wishart et mesure induite

Soient $`a\le b`$, $`G`$ une matrice de Ginibre $`a\times b`$ (coefficients indépendants, de parties réelle et imaginaire indépendantes de loi $`N(0,\frac12)`$,
de densité $`\pi^{-ab}e^{-\mathrm{Tr}\,GG^*}`$) et $`W=GG^*`$. Soit $`T_+`$ l'ensemble des matrices $`a\times a`$ triangulaires inférieures à
diagonale strictement positive, ouvert de $`\mathbb R^{a^2}`$ dans les coordonnées $`L_{ii}`$, $`\mathrm{Re}\,L_{ij}`$, $`\mathrm{Im}\,L_{ij}`$ ($`i>j`$), et
$`P_a`$ l'ensemble des matrices définies positives, ouvert dans l'espace des matrices hermitiennes $`a\times a`$.

<a id="r-C-11"></a>

**Lemme C.11 (décomposition de Bartlett).** Presque sûrement les lignes $`r_1,\dots,r_a`$ de $`G`$ sont linéairement indépendantes, et
Gram–Schmidt donne $`G=LV`$ avec $`L\in T_+`$ et $`V`$ à lignes orthonormées, donc $`W=LL^*`$. Les coefficients $`L_{ij}`$, $`i\ge j`$, sont
indépendants, $`L_{ij}`$ ($`i>j`$) est gaussien complexe standard et $`L_{ii}^2`$ a pour densité $`x^{b-i}e^{-x}/(b-i)!`$ sur $`(0,\infty)`$.
Par conséquent $`L`$ a pour densité $`c_L\prod_iL_{ii}^{2(b-i)+1}e^{-\mathrm{Tr}\,LL^*}`$ sur $`T_+`$.

*Démonstration.* Soit $`\mathcal F_{i-1}=\sigma(r_1,\dots,r_{i-1})`$ et supposons (récurrence) que $`r_1,\dots,r_{i-1}`$ sont p.s. indépendants, de
vecteurs de Gram–Schmidt $`v_1,\dots,v_{i-1}`$. Soit $`u_1,\dots,u_b`$ la base orthonormée obtenue par Gram–Schmidt à partir de
$`v_1,\dots,v_{i-1},e_1,\dots,e_b`$ (en écartant les vecteurs dépendants) ; $`u_k=v_k`$ pour $`k<i`$, et les $`u_k`$ sont des fonctions boréliennes de
$`r_1,\dots,r_{i-1}`$. Posons $`c_k=\langle u_k,r_i\rangle`$. Comme $`r_i`$ est indépendant de $`\mathcal F_{i-1}`$ et que sa densité $`\pi^{-b}e^{-\|r\|^2}`$ est
invariante par les unitaires (proposition A.16(a)), conditionnellement à $`\mathcal F_{i-1}`$ le vecteur $`c`$ est gaussien complexe standard (Fubini). Or
$`L_{ij}=c_j`$ pour $`j<i`$ et $`r_i-\sum_{j<i}c_jv_j=\sum_{k\ge i}c_ku_k`$, donc $`L_{ii}^2=\sum_{k=i}^b|c_k|^2>0`$ p.s. La loi conditionnelle de
la ligne $`i`$ de $`L`$ sachant $`\mathcal F_{i-1}`$ est donc une loi fixe, si bien que la ligne $`i`$ est indépendante des lignes précédentes, avec des coordonnées
$`c_1,\dots,c_{i-1}`$ et $`L_{ii}`$ indépendantes ; par récurrence toutes les lignes de $`L`$ sont indépendantes. Par le théorème A.20(b) avec $`N=m=b-i+1`$,
$`L_{ii}^2=\|(c_i,\dots,c_b)\|^2`$ a pour densité $`s^{m-1}e^{-s}/(m-1)!`$. La densité de $`L_{ii}=\sqrt s`$ est alors
$`2t^{2(b-i)+1}e^{-t^2}/(b-i)!`$, et $`\mathrm{Tr}\,LL^*=\sum_{i\ge j}|L_{ij}|^2`$. $`\square`$

*Calcul :* le lemme C.11 sur le facteur de Cholesky de $`W`$ — [tutoriel, section C.5](../../tutorial/fr/C_eigenvalue_densities.md#s-C-5).

<a id="r-C-12"></a>

**Lemme C.12 (application de Cholesky).** $`\chi(L)=LL^*`$ est une bijection $`C^\infty`$ de $`T_+`$ sur $`P_a`$ avec
$`|\det D\chi(L)|=2^a\prod_{i=1}^aL_{ii}^{2(a-i)+1}`$.

*Démonstration.* *Surjectivité :* pour $`W\in P_a`$ on applique le lemme 5.2 à $`W^{1/2}=QR`$ ; alors $`W=R^*R`$ et $`L=R^*\in T_+`$. *Injectivité :* si
$`LL^*=MM^*`$, alors $`M^{-1}L=M^*(L^*)^{-1}`$ est triangulaire inférieure et supérieure, donc diagonale à diagonale strictement positive, et
unitaire, donc égale à $`1`$. *Jacobien :* pour $`i\ge j`$, $`W_{ij}=L_{ij}L_{jj}+\sum_{k<j}L_{ik}\overline{L_{jk}}`$. Ordonnons les coordonnées de $`L`$ et du
triangle inférieur de $`W`$ colonne par colonne et, dans la colonne $`j`$, la diagonale en premier. Alors $`W_{ij}`$ ne dépend que de
$`L_{ij}`$ et de coordonnées antérieures, donc $`D\chi`$ est triangulaire, de coefficients diagonaux $`2L_{jj}`$ ($`\partial W_{jj}/\partial L_{jj}`$) et,
pour chaque $`i>j`$, deux fois $`L_{jj}`$ (parties réelle et imaginaire de $`W_{ij}`$). Le produit vaut $`\prod_j2L_{jj}\cdot L_{jj}^{2(a-j)}`$. $`\square`$

<a id="r-C-13"></a>

**Théorème C.13 (ensembles de Wishart et de Laguerre).** $`W=GG^*`$ a pour densité $`c_W\det(W)^{b-a}e^{-\mathrm{Tr}\,W}`$ sur $`P_a`$
(nulle ailleurs), et ses valeurs propres ont pour densité
```math
Z_{a,b}^{-1}\,\Delta(\lambda)^2\prod_{i=1}^a\lambda_i^{\,b-a}e^{-\lambda_i}\quad\text{on }C_+=\{0<\lambda_1<\dots<\lambda_a\}.
```

*Démonstration.* Pour $`g\ge0`$ borélienne, le lemme C.11 et le changement de variables $`W=\chi(L)`$ (lemme C.12) donnent
$`\mathbb E\,g(W)=\int_{P_a}g(W)\,c_L\prod_iL_{ii}^{2(b-i)+1}e^{-\mathrm{Tr}\,W}\big/\big(2^a\prod_iL_{ii}^{2(a-i)+1}\big)dW`$ avec
$`L=\chi^{-1}(W)`$, et $`\prod_iL_{ii}^{2(b-a)}=\det(W)^{b-a}`$. Cette densité est invariante par $`W\mapsto VWV^*`$, $`V\in U(a)`$, et
$`\Lambda\in P_a`$ si et seulement si $`\lambda_1>0`$ ; on applique le corollaire C.6 avec $`\beta=2`$. $`\square`$

<a id="r-C-14"></a>

**Théorème C.14 (mesure induite).** Soient $`\psi`$ un état pur aléatoire de Haar de $`\mathbb C^a\otimes\mathbb C^b`$ (section 7.5),
$`\rho_A=\mathrm{Tr}_2|\psi\rangle\langle\psi|=W/s`$ avec $`s=\mathrm{Tr}\,W`$ (démonstration du théorème 7.11), et
$`\mu_1<\dots<\mu_a`$ les valeurs propres de $`\rho_A`$. Alors $`s`$ et $`\mu`$ sont indépendants, $`s`$ a pour densité
$`s^{ab-1}e^{-s}/(ab-1)!`$, et $`(\mu_1,\dots,\mu_{a-1})`$ a pour densité
```math
Z'^{-1}_{a,b}\,\Delta(\mu)^2\prod_{i=1}^a\mu_i^{\,b-a}\quad\text{on }D_a=\{0<\mu_1<\dots<\mu_{a-1}<\mu_a\},\qquad\mu_a=1-\sum_{i<a}\mu_i ,
```
qui est la densité « $`\propto\delta(1-\sum_i\mu_i)\,\Delta(\mu)^2\prod_i\mu_i^{b-a}`$ » de la littérature physique.

*Démonstration.* $`\Theta(s,\mu_1,\dots,\mu_{a-1})=s\,(\mu_1,\dots,\mu_{a-1},1-\sum_{i<a}\mu_i)`$ est une bijection $`C^\infty`$ de $`(0,\infty)\times D_a`$ sur
$`C_+`$, de réciproque $`\lambda\mapsto(\sum\lambda_i,\lambda_1/\sum\lambda_i,\dots)`$. Ajouter toutes les lignes de $`D\Theta`$ à la dernière transforme celle-ci en
$`(1,0,\dots,0)`$, donc $`|\det D\Theta|=s^{a-1}`$. Comme $`\Delta`$ est homogène de degré $`a(a-1)/2`$, la densité de $`(s,\mu)`$ est
```math
Z_{a,b}^{-1}\,s^{a(a-1)+a(b-a)+a-1}e^{-s}\cdot\Delta(\mu)^2\prod_i\mu_i^{b-a}=Z_{a,b}^{-1}\,s^{ab-1}e^{-s}\cdot\Delta(\mu)^2\prod_i\mu_i^{b-a},
```
un produit, ce qui démontre l'indépendance et identifie les deux facteurs après normalisation. $`\square`$

*Calcul :* le théorème C.14 pour $`a=2`$ — [tutoriel, section C.5](../../tutorial/fr/C_eigenvalue_densities.md#s-C-5). *Dans la bibliothèque :* la mesure induite, `states.random_density_matrix` avec `rank=b` — [tutoriel, section C.5](../../tutorial/fr/C_eigenvalue_densities.md#s-C-5).

<a id="s-C-6"></a>

## C.6 Polynômes de Laguerre et formule de Page

<a id="r-C-15"></a>

**Lemme C.15 (polynômes de Laguerre).** Pour des entiers $`\alpha,k\ge0`$ soit
$`L_k^{(\alpha)}(x)=\sum_{j=0}^k(-1)^j\binom{k+\alpha}{k-j}\frac{x^j}{j!}`$, et $`L_{-1}^{(\alpha)}=0`$. Alors :

(a) $`L_k^{(\alpha)}(x)\,x^\alpha e^{-x}=\frac1{k!}\frac{d^k}{dx^k}\big(x^{k+\alpha}e^{-x}\big)`$ ; $`L_k^{(\alpha)}`$ est de degré $`k`$ et de coefficient dominant $`(-1)^k/k!`$.

(b) $`\frac{d}{dx}L_k^{(\alpha)}=-L_{k-1}^{(\alpha+1)}`$, donc $`\frac{d^i}{dx^i}L_k^{(\alpha)}=(-1)^iL_{k-i}^{(\alpha+i)}`$ pour $`i\le k`$ ; et $`L_k^{(\alpha)}=L_k^{(\alpha+1)}-L_{k-1}^{(\alpha+1)}`$.

(c) $`\int_0^\infty x^mL_k^{(\alpha)}x^\alpha e^{-x}dx=0`$ pour $`m<k`$, et $`\int_0^\infty L_k^{(\alpha)}L_l^{(\alpha)}x^\alpha e^{-x}dx=\delta_{kl}\frac{(k+\alpha)!}{k!}`$.

Par conséquent, d'après le lemme C.10 et le théorème C.13, les valeurs propres de $`W`$ ont pour densité de niveaux
($`\mathbb E\sum_i\varphi(\lambda_i(W))=\int_0^\infty\varphi\,\rho_{a,b}`$ pour toute $`\varphi\ge0`$ borélienne, et pour toute $`\varphi`$ telle que $`\int|\varphi|\rho_{a,b}<\infty`$)
```math
\rho_{a,b}(x)=\sum_{k=0}^{a-1}\frac{k!}{(k+\alpha)!}L_k^{(\alpha)}(x)^2\,x^\alpha e^{-x},\qquad\alpha=b-a .
```

*Démonstration.* (a) Par la formule de Leibniz, $`\frac{d^k}{dx^k}(x^{k+\alpha}e^{-x})=\sum_j\binom kj(-1)^j\frac{(k+\alpha)!}{(\alpha+j)!}x^{\alpha+j}e^{-x}`$ (avec $`j`$ dérivations
portant sur $`e^{-x}`$), et $`\frac1{k!}\binom kj\frac{(k+\alpha)!}{(\alpha+j)!}=\binom{k+\alpha}{k-j}\frac1{j!}`$. (b) En dérivant terme à terme,
$`-\sum_{i=0}^{k-1}(-1)^i\binom{k+\alpha}{k-1-i}\frac{x^i}{i!}=-L_{k-1}^{(\alpha+1)}`$ ; la seconde identité est la règle de Pascal
$`\binom{k+\alpha+1}{k-j}-\binom{k+\alpha}{k-1-j}=\binom{k+\alpha}{k-j}`$. (c) On intègre (a) par parties $`k`$ fois contre $`x^m`$ : les
termes de bord s'annulent, car $`\frac{d^r}{dx^r}(x^{k+\alpha}e^{-x})`$ contient le facteur $`x^{k+\alpha-r}`$, $`r<k`$, et décroît en $`\infty`$. Le
résultat est $`\frac{(-1)^k}{k!}\int\frac{d^kx^m}{dx^k}x^{k+\alpha}e^{-x}dx`$, qui vaut $`0`$ pour $`m<k`$ et $`(-1)^k(k+\alpha)!`$ pour $`m=k`$. Pour $`l<k`$,
$`L_l^{(\alpha)}`$ est une combinaison des $`x^m`$, $`m<k`$ ; pour $`l=k`$ seul le coefficient dominant $`(-1)^k/k!`$ compte. Ainsi les
$`p_k=\big(k!/(k+\alpha)!\big)^{1/2}L_k^{(\alpha)}`$ sont orthonormés pour $`w(x)=x^\alpha e^{-x}`$ sur $`(0,\infty)`$, et les valeurs propres non ordonnées
de $`W`$ ont la densité $`q`$ du lemme C.10 (avec $`n=a`$) d'après le corollaire C.6, de sorte que, comme dans la démonstration du corollaire C.8,
$`\mathbb E\sum_i\varphi(\lambda_i)=a\int\varphi(x_1)q(x)\,dx=\int\varphi\,K_aw=\int\varphi\,\rho_{a,b}`$ pour toute $`\varphi\ge0`$ borélienne, et par linéarité pour
$`\varphi`$ telle que $`\int|\varphi|\rho_{a,b}<\infty`$. $`\square`$

*Calcul :* la densité de niveaux $`\rho_{a,b}`$ comparée à un histogramme — [tutoriel, section C.6](../../tutorial/fr/C_eigenvalue_densities.md#s-C-6).

<a id="r-C-16"></a>

**Lemme C.16 (intégrales logarithmiques).** Soient $`I_m=\int_0^\infty x^m\log x\,e^{-x}dx`$ et $`H_m=\sum_{i=1}^m\frac1i`$. Alors
$`I_m=m!\,(I_0+H_m)`$. Pour des entiers $`\alpha\ge0`$ et $`0\le k\le m`$,
```math
A^{(\alpha)}_{k,m}:=\int_0^\infty\log x\,L_k^{(\alpha)}(x)L_m^{(\alpha)}(x)\,x^\alpha e^{-x}dx=
\begin{cases}I_{k+\alpha}/k!, & m=k,\\ -\dfrac{(k+\alpha)!}{k!\,(m-k)}, & m>k.\end{cases}
```
($`I_0`$ se simplifie dans le théorème C.17, sa valeur n'est donc pas nécessaire.)

*Démonstration.* Pour $`m\ge1`$, l'intégration par parties avec $`(-e^{-x})'=e^{-x}`$ donne $`I_m=mI_{m-1}+(m-1)!`$, c'est-à-dire
$`I_m/m!=I_{m-1}/(m-1)!+1/m`$. Pour $`A_{k,m}`$, on écrit $`L_m^{(\alpha)}x^\alpha e^{-x}=\frac1{m!}g^{(m)}`$ avec $`g=x^{m+\alpha}e^{-x}`$
(lemme C.15(a)) et $`f=\log x\,L_k^{(\alpha)}`$, et on intègre par parties $`m`$ fois. Les termes de bord
$`f^{(i)}g^{(m-1-i)}`$ s'annulent en $`\infty`$ et en $`0`$, où $`f^{(i)}=O(|\log x|+x^{-i})`$ et $`g^{(m-1-i)}=O(x^{\alpha+1+i})`$. Ainsi
$`A_{k,m}=\frac{(-1)^m}{m!}\int_0^\infty f^{(m)}g\,dx`$ (l'intégrande est $`O(x^\alpha(1+|\log x|))`$ en $`0`$) avec
$`f^{(m)}=\sum_{i\le k}\binom mi(\log x)^{(m-i)}\frac{d^i}{dx^i}L_k^{(\alpha)}`$. Pour $`m-i\ge1`$, $`(\log x)^{(m-i)}=(-1)^{m-i-1}(m-i-1)!\,x^{i-m}`$,
et d'après le lemme C.15(b,c)
```math
\int_0^\infty x^{i-m}\frac{d^iL_k^{(\alpha)}}{dx^i}\,x^{m+\alpha}e^{-x}dx=(-1)^i\int_0^\infty L_{k-i}^{(\alpha+i)}\,x^{\alpha+i}e^{-x}dx=(-1)^i(\alpha+i)!\,\delta_{ik}.
```
Si $`m>k`$, seul $`i=k`$ contribue : $`A_{k,m}=\frac{(-1)^m}{m!}\binom mk(-1)^{m-k-1}(m-k-1)!(-1)^k(\alpha+k)!=-\frac{(k+\alpha)!}{k!(m-k)}`$.
Si $`m=k`$, les termes $`i<k`$ s'annulent et le terme $`i=k`$ vaut $`\log x\cdot(-1)^k`$, ce qui donne $`A_{k,k}=\frac{(-1)^k}{k!}(-1)^kI_{k+\alpha}`$. $`\square`$

*Calcul :* les intégrales $`A^{(\alpha)}_{k,m}`$ par quadrature — [tutoriel, section C.6](../../tutorial/fr/C_eigenvalue_densities.md#s-C-6).

<a id="r-C-17"></a>

**Théorème C.17 (formule de Page).** Pour un état pur aléatoire de Haar de $`\mathbb C^a\otimes\mathbb C^b`$ avec $`a\le b`$,
```math
\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}.
```

*Démonstration.* *Réduction.* Avec les notations du théorème C.14, $`\lambda_i=s\mu_i`$ et $`\sum_i\mu_i=1`$, donc
$`\sum_i\lambda_i\log\lambda_i=s\log s-s\,S(\rho_A)`$. Tous les termes sont intégrables ($`|x\log x|\le1+x^2`$, et $`\mathbb Es^2<\infty`$,
$`0\le S(\rho_A)\le a/e`$ car $`0\le-x\log x\le1/e`$ sur $`[0,1]`$), $`\mathbb E\sum_i\lambda_i\log\lambda_i=\int x\log x\,\rho_{a,b}(x)\,dx`$ d'après le lemme C.15, et $`s`$ est indépendant de
$`S(\rho_A)`$, donc $`\mathbb E\sum_i\lambda_i\log\lambda_i=\mathbb E[s\log s]-\mathbb E[s]\,\mathbb E\,S(\rho_A)`$. D'après la densité de $`s`$,
$`\mathbb Es=ab`$ et $`\mathbb E[s\log s]=I_{ab}/(ab-1)!=ab\,(I_0+H_{ab})`$ (lemme C.16). Donc
```math
\mathbb E\,S(\rho_A)=I_0+H_{ab}-\frac1{ab}\int_0^\infty x\log x\,\rho_{a,b}(x)\,dx .
```
*L'intégrale.* D'après le lemme C.15(b), $`L_k^{(\alpha)}=L_k^{(\alpha+1)}-L_{k-1}^{(\alpha+1)}`$, donc avec $`x\cdot x^\alpha=x^{\alpha+1}`$ et le lemme C.16
pour $`\alpha+1`$,
```math
\frac{k!}{(k+\alpha)!}\int_0^\infty x\log x\,L_k^{(\alpha)}(x)^2x^\alpha e^{-x}dx=\frac{k!}{(k+\alpha)!}\Big(A^{(\alpha+1)}_{k,k}-2A^{(\alpha+1)}_{k-1,k}+A^{(\alpha+1)}_{k-1,k-1}\Big)
```
```math
=(k+\alpha+1)(I_0+H_{k+\alpha+1})+2k+k\,(I_0+H_{k+\alpha}),
```
y compris pour $`k=0`$. Comme $`\sum_{k=0}^{a-1}(2k+\alpha+1)=a(a-1)+a(b-a+1)=ab`$,
```math
\int_0^\infty x\log x\,\rho_{a,b}\,dx=ab\,I_0+\Sigma,\qquad \Sigma=\sum_{k=0}^{a-1}\big[(k+\alpha+1)H_{k+\alpha+1}+kH_{k+\alpha}+2k\big].
```
*La somme.* Posons $`m=\alpha+1=b-a+1`$ et $`j=k+m\in\{m,\dots,b\}`$. Comme $`(j-m)H_{j-1}=(j-m)H_j-1+\frac mj`$,
$`\Sigma=\sum_{j=m}^b(2j-m)H_j-a+m(H_b-H_{m-1})+a(a-1)`$. Avec $`F(j)=j(j+1-m)`$ on a $`F(j)-F(j-1)=2j-m`$, $`F(m-1)=0`$,
$`F(b)=ab`$, et la sommation par parties donne
```math
\sum_{j=m}^b(2j-m)H_j=abH_b-\sum_{j=m}^{b-1}\frac{F(j)}{j+1}=abH_b-\sum_{j=m}^{b-1}\Big(j-m+\frac m{j+1}\Big)=abH_b-\frac{(a-1)(a-2)}2-m(H_b-H_m).
```
Comme $`m(H_m-H_{m-1})=1`$, $`\Sigma=abH_b-\frac{(a-1)(a-2)}2+1-a+a(a-1)=abH_b+\frac{a(a-1)}2`$. Donc
$`\mathbb E\,S(\rho_A)=H_{ab}-H_b-\frac{a-1}{2b}`$. $`\square`$

*Calcul :* l'identité pour $`\Sigma`$ en arithmétique exacte, et la formule de Page par quadrature ($`a=2`$) et par Monte Carlo — [tutoriel, section C.6](../../tutorial/fr/C_eigenvalue_densities.md#s-C-6). *Dans la bibliothèque :* l'entropie $`S(\rho_A)`$ est `DensityMatrix.von_neumann_entropy` — [tutoriel, section 3.4](../../tutorial/fr/03_traces_and_states.md#s-3-4).

<a id="r-C-18"></a>

**Corollaire C.18.** Pour $`a\le b`$,
$`\log a-\frac{a-1}{ab+1}-\frac{a-1}{2b}\le\mathbb E\,S(\rho_A)\le\log a-\frac{a-1}{2b}`$ ; en particulier
$`\big|\mathbb E\,S(\rho_A)-\log a+\frac a{2b}\big|\le\frac1{2b}`$. Pour deux qubits, $`\mathbb E\,S(\rho_A)=\frac13`$.

*Démonstration.* Comme $`1/x`$ est décroissante, $`\int_k^{k+1}\frac{dx}x\le\frac1k\le\int_{k-1}^k\frac{dx}x`$ ; en sommant sur $`b<k\le ab`$,
$`\log\frac{ab+1}{b+1}\le H_{ab}-H_b\le\log a`$, et $`\log\frac{ab+1}{b+1}=\log a-\log\big(1+\frac{a-1}{ab+1}\big)\ge\log a-\frac{a-1}{ab+1}`$.
Enfin $`\frac{a-1}{ab+1}\le\frac1b`$, et $`a=b=2`$ donne $`\frac13+\frac14-\frac14`$. $`\square`$

*Calcul :* les deux bornes pour $`a<12`$, $`b<40`$ — [tutoriel, section C.6](../../tutorial/fr/C_eigenvalue_densities.md#s-C-6).

La formule a été conjecturée par D. N. Page (Phys. Rev. Lett. 71, 1291, 1993) et démontrée par S. K. Foong et S. Kanno
(Phys. Rev. Lett. 72, 1148, 1994), J. Sánchez-Ruiz (Phys. Rev. E 52, 5653, 1995) et S. Sen (Phys. Rev. Lett. 77, 1, 1996) ;
le calcul ci-dessus suit la réduction de Sen à l'ensemble de Laguerre, mais somme la densité à un point terme à terme
au lieu d'utiliser la formule de Christoffel–Darboux. La mesure induite est celle de K. Życzkowski et H.-J. Sommers
(J. Phys. A 34, 7111, 2001) ; la décomposition du lemme C.11 est due à M. S. Bartlett (1933), et le théorème C.5 à H. Weyl.

**Retour :** le [chapitre 5](05_random_matrices.md) utilise le théorème C.7, le corollaire C.8 et la proposition C.9 ;
le [chapitre 7](07_tensor_products.md) utilise les théorèmes C.14 et C.17 et le corollaire C.18.
