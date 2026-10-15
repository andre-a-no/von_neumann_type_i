[English](../D_gksl.md) · [Español](../es/D_gksl.md) · **Français** · [Deutsch](../de/D_gksl.md) · [中文](../zh/D_gksl.md) · [日本語](../ja/D_gksl.md) · [Русский](../ru/D_gksl.md)

# Annexe D. Générateurs des semi-groupes dynamiques quantiques

Le chapitre 8 a démontré qu'un générateur de la forme de Gorini–Kossakowski–Sudarshan–Lindblad (GKSL) engendre un
semi-groupe dynamique quantique (théorème 8.6) et a utilisé la réciproque (théorème 8.7) et sa version pour
$`\mathcal M=\bigoplus_cM_{k_c}`$. Cette annexe démontre les deux en dimension finie, ainsi que la
liberté dans la forme GKSL et un critère d'unicité de l'état stationnaire. L'outil est la
matrice de Choi du chapitre 6 : un générateur n'est en général pas complètement positif, mais sa matrice de Choi est
positive sur l'orthogonal du vecteur maximalement intriqué. Le même type d'argument — positivité
d'une famille d'opérateurs, puis une décomposition — démontre dans la dernière section le théorème de correction
d'erreurs par algèbres d'opérateurs énoncé au chapitre 10 (théorème 10.7), avec un canal de récupération explicite.

**Vous apprendrez** pourquoi le générateur d'un semi-groupe complètement positif est conditionnellement complètement
positif, comment lire $`H`$ et les opérateurs de saut sur la matrice de Choi d'un générateur, en quel sens
ils sont uniques, pourquoi tout générateur sur une algèbre à secteurs a la forme utilisée à la section 8.3, quand
l'état stationnaire est unique, et quelles algèbres d'observables un canal de récupération peut protéger.

Dans toute cette annexe, $`M_n=M_n(\mathbb C)`$, les $`E_{ij}=|i\rangle\langle j|`$ sont ses unités matricielles, un
**superopérateur** est une application linéaire sur $`M_n`$ (ou sur $`\mathcal M`$), et
```math
\Omega=\sum_{i=1}^ne_i\otimes e_i\in\mathbb C^n\otimes\mathbb C^n,\qquad
J(\Psi)=\sum_{i,j=1}^nE_{ij}\otimes\Psi(E_{ij})=[\Psi(E_{ij})]_{i,j=1}^n
```
sont le vecteur maximalement intriqué et la **matrice de Choi** d'un superopérateur $`\Psi`$ (théorème 6.2 ;
ici $`\Psi`$ n'est pas nécessairement positif). Pour des vecteurs $`v,w`$, $`|v\rangle\langle w|`$ est l'opérateur
$`x\mapsto\langle w,x\rangle v`$. Un superopérateur $`\Psi`$ **préserve l'hermiticité** si
$`\Psi(A^*)=\Psi(A)^*`$ pour tout $`A`$. Les semi-groupes dynamiques quantiques et leurs générateurs sont ceux de la
section 8.2 (lemme 8.5) : $`T_t=e^{t\mathcal L}`$ avec des canaux $`T_t`$.

*Calcul :* les vérifications numériques de cette annexe — [tutoriel, annexe D](../../tutorial/fr/D_gksl.md).

<a id="s-D-1"></a>

## D.1 Matrices de Choi des briques élémentaires

<a id="r-D-1"></a>

**Lemme D.1.** Soit $`\Psi`$ un superopérateur sur $`M_n`$.

(a) $`\Psi\mapsto J(\Psi)`$ est une bijection linéaire de l'espace des superopérateurs sur $`M_n`$ sur $`M_n(M_n)=M_{n^2}`$.
(b) Tout $`v\in\mathbb C^n\otimes\mathbb C^n`$ s'écrit $`v=(1\otimes K)\Omega`$ pour exactement un $`K\in M_n`$, et
$`\langle(1\otimes K)\Omega,(1\otimes L)\Omega\rangle=\mathrm{Tr}(K^*L)`$. En particulier,
$`\langle\Omega,(1\otimes K)\Omega\rangle=\mathrm{Tr}\,K`$.
(c) Pour $`K,L\in M_n`$, le superopérateur $`A\mapsto KAL^*`$ a pour matrice de Choi
$`|(1\otimes K)\Omega\rangle\langle(1\otimes L)\Omega|`$. En particulier, $`J(\mathrm{id})=|\Omega\rangle\langle\Omega|`$.
(d) $`\Psi`$ préserve l'hermiticité si et seulement si $`J(\Psi)`$ est autoadjointe. Toute application CP préserve l'hermiticité.

*Démonstration.* (a) Le bloc $`(i,j)`$ de $`J(\Psi)`$ est $`\Psi(E_{ij})`$ et les $`E_{ij}`$ forment une base de $`M_n`$, donc
$`J`$ est linéaire et injective, et toute matrice par blocs $`[X_{ij}]`$ est la matrice de Choi du superopérateur
défini par $`E_{ij}\mapsto X_{ij}`$. (b) $`(1\otimes K)\Omega=\sum_ie_i\otimes Ke_i`$ a pour $`i`$-ième bloc $`Ke_i`$, la
$`i`$-ième colonne de $`K`$, qui détermine $`K`$, et tout $`v`$ s'obtient ainsi. Ensuite
$`\langle\sum_ie_i\otimes Ke_i,\sum_je_j\otimes Le_j\rangle=\sum_i\langle Ke_i,Le_i\rangle=\mathrm{Tr}(K^*L)`$ ;
pour la dernière assertion, on remplace $`(K,L)`$ par $`(1,K)`$. (c) Le bloc $`(i,j)`$ du membre de droite est
$`Ke_i(Le_j)^*=KE_{ij}L^*`$. (d) Le bloc $`(i,j)`$ de $`J(\Psi)^*`$ est $`\Psi(E_{ji})^*`$. Si $`\Psi`$
préserve l'hermiticité, $`\Psi(E_{ji})^*=\Psi(E_{ji}^*)=\Psi(E_{ij})`$. Réciproquement, si
$`\Psi(E_{ji})^*=\Psi(E_{ij})`$ pour tous $`i,j`$ et $`A=\sum a_{ij}E_{ij}`$, alors
$`\Psi(A)^*=\sum\bar a_{ij}\Psi(E_{ij})^*=\sum\bar a_{ij}\Psi(E_{ji})=\Psi(A^*)`$. Une application CP a une forme de Kraus
(théorème 6.2), et $`(\sum_lK_lAK_l^*)^*=\sum_lK_lA^*K_l^*`$. $`\square`$

<a id="r-D-2"></a>

**Lemme D.2 (décomposition).** Soient $`q_1,\dots,q_C`$ des vecteurs orthonormés de
$`\mathbb C^n\otimes\mathbb C^n`$, écrits $`q_c=(1\otimes Q_c)\Omega`$, soit $`W=\mathrm{span}\{q_c\}`$, et soit
$`\Psi`$ un superopérateur qui préserve l'hermiticité, avec $`\langle v,J(\Psi)v\rangle\ge0`$ pour tout $`v\perp W`$.
Posons $`J=J(\Psi)`$ et $`\xi_c=Jq_c-\frac12\sum_d\langle q_d,Jq_c\rangle q_d`$, et soit $`K_c\in M_n`$ défini par
$`(1\otimes K_c)\Omega=\xi_c`$. Alors il existe $`L_1,\dots,L_r\in M_n`$, où $`r=\mathrm{rank}(P_W^\perp JP_W^\perp)`$
et $`P_W^\perp`$ est le projecteur sur $`W^\perp`$, tels que
```math
\Psi(A)=\sum_{j=1}^rL_jAL_j^*+\sum_{c=1}^C\big(K_cAQ_c^*+Q_cAK_c^*\big)\qquad(A\in M_n),
```
$`\mathrm{Tr}(Q_c^*L_j)=0`$ pour tous $`c,j`$, et $`\mathrm{Tr}(L_j^*L_k)=0`$ pour $`j\neq k`$.

*Démonstration.* $`J`$ est autoadjointe d'après le lemme D.1(d). Soient $`P=\sum_c|q_c\rangle\langle q_c|`$ le projecteur sur $`W`$
et $`P^\perp=P_W^\perp=1-P`$. En développant $`J=(P+P^\perp)J(P+P^\perp)`$ et en utilisant $`PJP^\perp=PJ-PJP`$ et
$`P^\perp JP=JP-PJP`$,
```math
J=P^\perp JP^\perp+\big(PJ+JP-PJP\big).
```
Le premier terme est autoadjoint, et $`\langle v,P^\perp JP^\perp v\rangle=\langle P^\perp v,JP^\perp v\rangle\ge0`$
car $`P^\perp v\perp W`$ ; donc $`P^\perp JP^\perp\ge0`$. D'après le théorème 1.3,
$`P^\perp JP^\perp=\sum_{j=1}^r|l_j\rangle\langle l_j|`$, où les $`l_j`$ sont les vecteurs propres des valeurs propres non nulles
$`\mu_j`$, multipliés par $`\sqrt{\mu_j}`$ : ils sont deux à deux orthogonaux et appartiennent à l'image de
$`P^\perp JP^\perp`$, qui est orthogonale à $`W`$. Écrivons $`l_j=(1\otimes L_j)\Omega`$ ; d'après le lemme D.1(b),
$`\mathrm{Tr}(Q_c^*L_j)=\langle q_c,l_j\rangle=0`$ et $`\mathrm{Tr}(L_j^*L_k)=\langle l_j,l_k\rangle=0`$ pour $`j\ne k`$,
et d'après le lemme D.1(c), $`P^\perp JP^\perp`$ est la matrice de Choi de $`A\mapsto\sum_jL_jAL_j^*`$. Pour le second terme,
$`JP=\sum_c|Jq_c\rangle\langle q_c|`$ et $`PJP=\sum_{c,d}\langle q_d,Jq_c\rangle|q_d\rangle\langle q_c|`$, donc
```math
\sum_c|\xi_c\rangle\langle q_c|=JP-\tfrac12PJP,\qquad\sum_c|q_c\rangle\langle\xi_c|=PJ-\tfrac12PJP
```
(la seconde est l'adjointe de la première ; $`PJP`$ est autoadjointe), et leur somme vaut $`PJ+JP-PJP`$. D'après le
lemme D.1(c), $`|\xi_c\rangle\langle q_c|`$ et $`|q_c\rangle\langle\xi_c|`$ sont les matrices de Choi de
$`A\mapsto K_cAQ_c^*`$ et $`A\mapsto Q_cAK_c^*`$. Ainsi les deux membres de l'identité annoncée ont pour matrice de Choi $`J`$,
et le lemme D.1(a) donne l'identité. $`\square`$

<a id="s-D-2"></a>

## D.2 Positivité complète conditionnelle

Un superopérateur $`\Psi`$ sur $`M_n`$ est **conditionnellement complètement positif** (CCP) si
$`\langle v,J(\Psi)v\rangle\ge0`$ pour tout $`v\perp\Omega`$. D'après le lemme D.1(b), ces $`v`$ sont les vecteurs
$`(1\otimes K)\Omega`$ avec $`\mathrm{Tr}\,K=0`$.

<a id="r-D-3"></a>

**Théorème D.3.** Pour un superopérateur $`\mathcal L`$ sur $`M_n`$, les assertions suivantes sont équivalentes :

1. $`e^{t\mathcal L}`$ est CP pour tout $`t\ge0`$ ;
2. $`\mathcal L`$ préserve l'hermiticité et est conditionnellement complètement positif ;
3. il existe une application CP $`\Phi`$ sur $`M_n`$ et $`G\in M_n`$ tels que $`\mathcal L(A)=\Phi(A)+GA+AG^*`$ pour tout $`A`$.

Dans 3, on peut prendre $`\Phi(A)=\sum_{j=1}^rL_jAL_j^*`$ avec $`r\le n^2-1`$, $`\mathrm{Tr}\,L_j=0`$ et
$`\mathrm{Tr}(L_j^*L_k)=0`$ pour $`j\neq k`$, et $`G`$ avec $`\mathrm{Tr}\,G\in\mathbb R`$.

*Démonstration.* (1)$`\Rightarrow`$(2). D'après le lemme 8.1, $`(e^{t\mathcal L}-\mathrm{id})/t\to\mathcal L`$ quand $`t\to0`$, et $`J`$ est
linéaire sur un espace de dimension finie, donc continue ; ainsi
$`J(e^{t\mathcal L})=J(\mathrm{id})+t\,J\big((e^{t\mathcal L}-\mathrm{id})/t\big)`$, où la dernière matrice de Choi tend vers
$`J(\mathcal L)`$. Pour $`v\perp\Omega`$, $`\langle v,J(\mathrm{id})v\rangle=|\langle\Omega,v\rangle|^2=0`$ (lemme D.1(c)), donc
pour $`t>0`$, comme $`J(e^{t\mathcal L})\ge0`$ (théorème 6.2),
```math
0\le\frac1t\langle v,J(e^{t\mathcal L})v\rangle=\big\langle v,J\big((e^{t\mathcal L}-\mathrm{id})/t\big)v\big\rangle
\;\xrightarrow[t\to0]{}\;\langle v,J(\mathcal L)v\rangle .
```
Chaque $`e^{t\mathcal L}`$ est CP, donc préserve l'hermiticité (lemme D.1(d)), d'où
$`\mathcal L(A^*)=\lim_t(e^{t\mathcal L}(A)^*-A^*)/t=\mathcal L(A)^*`$.

(2)$`\Rightarrow`$(3). On applique le lemme D.2 avec $`C=1`$ et $`q_1=\Omega/\sqrt n`$, c'est-à-dire $`Q_1=1/\sqrt n`$ :
$`\mathcal L(A)=\sum_jL_jAL_j^*+GA+AG^*`$ avec $`G=K_1/\sqrt n`$, et $`\Phi=\sum_jL_j\cdot L_j^*`$ est CP
(théorème 6.2). De plus, $`\mathrm{Tr}\,L_j=\sqrt n\mathrm{Tr}(Q_1^*L_j)=0`$, $`r\le\mathrm{rank}\,P_W^\perp=n^2-1`$,
et d'après le lemme D.1(b),
$`\mathrm{Tr}\,G=\langle\Omega,\xi_1\rangle/\sqrt n=\langle q_1,\xi_1\rangle=\frac12\langle q_1,J(\mathcal L)q_1\rangle\in\mathbb R`$,
puisque $`J(\mathcal L)`$ est autoadjointe.

(3)$`\Rightarrow`$(1). Les étapes 1 à 4 de la démonstration du théorème 8.6 n'utilisent que le fait que le générateur est
$`\mathcal L_0+\Phi`$ avec $`\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)`$ pour une certaine matrice
$`H_{\mathrm{eff}}`$ et $`\Phi`$ CP. Avec $`H_{\mathrm{eff}}=iG`$, on a $`\mathcal L_0(\rho)=G\rho+\rho G^*`$, et ces
étapes montrent que $`e^{t\mathcal L}=\sum_{m\ge0}S_m(t)`$ est CP. $`\square`$

*Calcul :* le théorème D.3 sur un générateur GKSL sur $`M_3`$, et l'exemple $`\theta-\mathrm{id}`$ ci-dessous — [tutoriel, section D.2](../../tutorial/fr/D_gksl.md#s-D-2). *Dans la bibliothèque :* la matrice d'un générateur GKSL, `dynamics.lindblad_superoperator` — [tutoriel, section 8.2](../../tutorial/fr/08_dynamics.md#s-8-2).

La positivité seule ne suffirait pas. Soient $`n\ge2`$, $`\theta`$ la transposition sur $`M_n`$ et $`\mathcal L=\theta-\mathrm{id}`$. Comme $`\theta^2=\mathrm{id}`$,
```math
e^{t\mathcal L}=e^{-t}\big(\cosh t\,\mathrm{id}+\sinh t\,\theta\big)=\tfrac{1+e^{-2t}}2\,\mathrm{id}+\tfrac{1-e^{-2t}}2\,\theta ,
```
combinaison convexe d'applications positives qui préservent la trace (proposition 6.1). Mais $`J(\theta)=[E_{ji}]=\sum_{i,j}E_{ij}\otimes E_{ji}`$ est l'opérateur d'échange $`e_k\otimes e_l\mapsto e_l\otimes e_k`$, qui envoie donc
$`v=e_1\otimes e_2-e_2\otimes e_1`$, orthogonal à $`\Omega`$, sur $`-v`$ (comme dans la démonstration de la proposition 6.1) ; avec
$`\|v\|^2=2`$ et $`\langle v,J(\mathrm{id})v\rangle=|\langle\Omega,v\rangle|^2=0`$, on obtient
$`\langle v,J(\mathcal L)v\rangle=-2`$ : $`\mathcal L`$ n'est pas CCP, et de fait
$`\langle v,J(e^{t\mathcal L})v\rangle=-(1-e^{-2t})<0`$ pour tout $`t>0`$. Un semi-groupe d'applications positives qui préservent
la trace n'est pas nécessairement formé de canaux.

<a id="s-D-3"></a>

## D.3 Le théorème GKSL

<a id="r-D-4"></a>

**Théorème D.4 (Gorini–Kossakowski–Sudarshan ; Lindblad ; théorème 8.7).** Un superopérateur $`\mathcal L`$ sur
$`M_n`$ est le générateur d'un semi-groupe dynamique quantique sur $`M_n`$ si et seulement s'il existe $`H=H^*`$ et
$`L_1,\dots,L_r`$ dans $`M_n`$ tels que
```math
\mathcal L(\rho)=-i[H,\rho]+\sum_{j=1}^r\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
```
On peut choisir $`r\le n^2-1`$, $`\mathrm{Tr}\,H=0`$, $`\mathrm{Tr}\,L_j=0`$ et $`\mathrm{Tr}(L_j^*L_k)=0`$ pour $`j\neq k`$.

*Démonstration.* Le sens « si » est le théorème 8.6. Réciproquement, soit $`(T_t)`$ un semi-groupe dynamique quantique de générateur
$`\mathcal L`$, de sorte que $`T_t=e^{t\mathcal L}`$ (lemme 8.5). Comme tout $`T_t`$ est CP, le théorème D.3 donne
$`\mathcal L(\rho)=\sum_jL_j\rho L_j^*+G\rho+\rho G^*`$ avec $`r\le n^2-1`$ opérateurs $`L_j`$ de trace nulle, deux à deux orthogonaux, et
$`\mathrm{Tr}\,G\in\mathbb R`$. Posons $`H=\frac i2(G-G^*)`$ et $`R=\frac12(G+G^*)`$, tous deux autoadjoints ; alors
$`G=-iH+R`$, $`G^*=iH+R`$ et $`G\rho+\rho G^*=-i[H,\rho]+\{R,\rho\}`$. De plus,
$`\mathrm{Tr}\,H=\frac i2(\mathrm{Tr}\,G-\overline{\mathrm{Tr}\,G})=0`$. Comme tout $`T_t`$ est TP,
$`\mathrm{Tr}\,\mathcal L(\rho)=\frac d{dt}\mathrm{Tr}\,T_t(\rho)\big|_{t=0}=0`$ (lemme 8.1), et par cyclicité
```math
0=\mathrm{Tr}\,\mathcal L(\rho)=\mathrm{Tr}\Big(\rho\,\Big(\sum_jL_j^*L_j+2R\Big)\Big)\qquad\text{for all }\rho .
```
Avec $`\rho=X=\sum_jL_j^*L_j+2R`$, qui est autoadjoint, $`\mathrm{Tr}(X^*X)=0`$, donc $`X=0`$ et
$`R=-\frac12\sum_jL_j^*L_j`$, ce qui est la forme GKSL. $`\square`$

*Calcul :* la démonstration du théorème D.4 exécutée comme un algorithme, qui retrouve $`H`$ et les sauts à partir de $`\mathcal L`$ — [tutoriel, section D.3](../../tutorial/fr/D_gksl.md#s-D-3).

<a id="r-D-5"></a>

**Théorème D.5 (liberté dans la forme GKSL).** Soient $`(H;L_1,\dots,L_r)`$ et $`(H';L'_1,\dots,L'_s)`$ des données GKSL
($`H,H'`$ autoadjoints) de générateurs $`\mathcal L`$, $`\mathcal L'`$ comme dans le théorème D.4.

(a) $`\mathcal L'=\mathcal L`$ dans chacun des cas suivants :
(i) $`L'_j=L_j`$ et $`H'=H+c1`$, $`c\in\mathbb R`$ ;
(ii) $`H'=H`$ et, après avoir complété par des zéros la famille la plus courte, $`L'_j=\sum_ku_{jk}L_k`$ avec une matrice unitaire $`u`$ ;
(iii) $`L'_j=L_j+a_j1`$ avec $`a_j\in\mathbb C`$, et $`H'=H-\frac i2\sum_j(\bar a_jL_j-a_jL_j^*)`$, qui est autoadjoint.

(b) Si $`\mathrm{Tr}\,H=\mathrm{Tr}\,H'=0`$, si tous les $`L_j,L'_j`$ sont de trace nulle et si $`\mathcal L=\mathcal L'`$, alors $`H=H'`$ et,
après complétion par des zéros, $`L'_j=\sum_ku_{jk}L_k`$ avec une matrice unitaire $`u`$.

(c) Par conséquent, deux données GKSL ont le même générateur si et seulement si l'une s'obtient de l'autre par un nombre fini
d'opérations (i)–(iii). En particulier, l'opérateur $`H`$ de la forme normalisée du théorème D.4 et sa
partie dissipative $`\Phi(\rho)=\sum_jL_j\rho L_j^*`$ sont déterminés par $`\mathcal L`$.

*Démonstration.* (a) (i) $`[c1,\rho]=0`$. (ii) D'après le théorème 6.7, $`\Phi'=\sum_jL'_j\cdot L_j'^*`$ est égale à
$`\Phi=\sum_jL_j\cdot L_j^*`$ ; la duale ne dépend pas des opérateurs de Kraus (section 6.2), donc
$`\sum_jL_j'^*L'_j=\Phi^*(1)=\sum_jL_j^*L_j`$. (iii) Pour un seul saut $`L`$ et $`a\in\mathbb C`$, le développement donne
```math
(L+a)\rho(L+a)^*-\tfrac12\{(L+a)^*(L+a),\rho\}-\Big(L\rho L^*-\tfrac12\{L^*L,\rho\}\Big)=\tfrac12[B,\rho],
\qquad B=\bar aL-aL^* ,
```
car les termes $`\bar aL\rho+a\rho L^*`$ du premier produit se combinent avec
$`-\frac12(aL^*\rho+a\rho L^*+\bar aL\rho+\bar a\rho L)`$ et les termes en $`|a|^2`$ se compensent. Comme $`B^*=-B`$, la
matrice $`\frac i2B`$ est autoadjointe, donc $`H'=H-\frac i2\sum_jB_j`$ est autoadjoint, et
$`-i[H'-H,\rho]=-\frac12\sum_j[B_j,\rho]`$ compense la variation du dissipateur.

(b) Soient $`\Phi=\sum_jL_j\cdot L_j^*`$ et $`G=-iH-\frac12\Phi^*(1)`$, de sorte que $`\mathcal L(\rho)=\Phi(\rho)+G\rho+\rho G^*`$.
D'après le lemme D.1(c),
$`J(\mathcal L)=\sum_j|l_j\rangle\langle l_j|+|g\rangle\langle\Omega|+|\Omega\rangle\langle g|`$ avec $`l_j=(1\otimes L_j)\Omega`$ et
$`g=(1\otimes G)\Omega`$, et $`l_j\perp\Omega`$ car $`\mathrm{Tr}\,L_j=0`$ (lemme D.1(b)). En notant $`P^\perp`$ le projecteur
sur $`\Omega^\perp`$, on a donc $`P^\perp J(\mathcal L)P^\perp=\sum_j|l_j\rangle\langle l_j|=J(\Phi)`$, et de même
$`P^\perp J(\mathcal L')P^\perp=J(\Phi')`$. Ainsi $`J(\Phi)=J(\Phi')`$, $`\Phi=\Phi'`$ (lemme D.1(a)), et le théorème 6.7 fournit la
matrice unitaire $`u`$. Alors $`\Phi^*(1)=\Phi'^*(1)`$, donc $`[H-H',\rho]=0`$ pour tout $`\rho`$. Une matrice $`X`$ qui commute avec toutes
les $`E_{ij}`$ est un multiple de $`1`$ (comparer les coefficients de $`XE_{ij}`$ et de $`E_{ij}X`$), donc $`H-H'=\lambda1`$, et
$`\lambda=0`$ en prenant les traces.

(c) D'après (a), les opérations (i)–(iii) ne changent pas le générateur, et l'inverse de chaque opération est une opération du
même type ((iii) avec $`-a_j`$, puisque $`\bar a_jL'_j-a_jL_j'^*=\bar a_jL_j-a_jL_j^*`$). Toute donnée se ramène à la
forme normalisée par (iii) avec $`a_j=-\mathrm{Tr}\,L_j/n`$, suivie de (i) avec $`c=-\mathrm{Tr}\,H''/n`$, où $`H''`$ est
l'hamiltonien après (iii) ; deux données normalisées de même générateur sont reliées par (ii), d'après (b). $`\square`$

*Calcul :* le théorème D.5 : la forme normalisée par (iii) et (i), et l'unicité dans (b) — [tutoriel, section D.3](../../tutorial/fr/D_gksl.md#s-D-3).

<a id="r-D-6"></a>

**Théorème D.6 (unicité de l'état stationnaire).** Soit $`\mathcal L`$ le générateur GKSL du théorème D.4
de données $`(H;L_1,\dots,L_r)`$. Si $`\mathcal L`$ a un état stationnaire fidèle $`\sigma`$ ($`\sigma>0`$) et si les seules
matrices qui commutent avec $`H`$ et avec tous les $`L_j`$, $`L_j^*`$ sont les multiples de $`1`$, alors $`\ker\mathcal L=\mathbb C\sigma`$ :
$`\sigma`$ est le seul état stationnaire. De plus, $`\ker\mathcal L^*=\mathbb C1`$.

*Démonstration.* La duale pour l'accouplement par la trace, $`\mathrm{Tr}(\mathcal L(\rho)A)=\mathrm{Tr}(\rho\,\mathcal L^*(A))`$, est
$`\mathcal L^*(A)=i[H,A]+\sum_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})`$ : par cyclicité,
$`\mathrm{Tr}(-i[H,\rho]A)=\mathrm{Tr}(\rho\,i[H,A])`$, $`\mathrm{Tr}(L\rho L^*A)=\mathrm{Tr}(\rho L^*AL)`$ et
$`\mathrm{Tr}(\{L^*L,\rho\}A)=\mathrm{Tr}(\rho\{L^*L,A\})`$ (c'est l'exercice 1 du chapitre 8).
Un développement direct, dans lequel les termes hamiltoniens se compensent parce que $`A\mapsto i[H,A]`$ est une dérivation, donne
pour tout $`X`$
```math
\mathcal L^*(X^*X)-X^*\mathcal L^*(X)-\mathcal L^*(X^*)X=\sum_j[L_j,X]^*[L_j,X]:
```
pour un seul saut, parmi les six termes d'anticommutateur, quatre se compensent deux à deux et les deux autres,
$`\frac12X^*L^*LX+\frac12X^*L^*LX`$, s'additionnent, de sorte que le membre de gauche vaut $`L^*X^*XL-X^*L^*XL-L^*X^*LX+X^*L^*LX`$, c'est-à-dire
$`(XL-LX)^*(XL-LX)`$. Soit $`X\in\ker\mathcal L^*`$. Alors $`\mathcal L^*(X^*)=\mathcal L^*(X)^*=0`$, car $`\mathcal L^*`$
préserve l'hermiticité, donc $`\mathcal L^*(X^*X)=\sum_j[L_j,X]^*[L_j,X]=:Y\ge0`$, et
$`\mathrm{Tr}(\sigma Y)=\mathrm{Tr}(\mathcal L(\sigma)X^*X)=0`$. Comme $`\sigma>0`$,
$`\mathrm{Tr}(\sigma Y)=\mathrm{Tr}(\sigma^{1/2}Y\sigma^{1/2})`$, avec $`\sigma^{1/2}Y\sigma^{1/2}\ge0`$, ne s'annule que si
$`\sigma^{1/2}Y\sigma^{1/2}=0`$, c'est-à-dire $`Y=0`$ ; en prenant les traces, $`\sum_j\mathrm{Tr}([L_j,X]^*[L_j,X])=0`$ avec des termes
positifs, donc $`[L_j,X]=0`$ pour tout $`j`$. Le même argument appliqué à $`X^*\in\ker\mathcal L^*`$
donne $`[L_j,X^*]=0`$, c'est-à-dire $`[L_j^*,X]=0`$. Alors $`X`$ commute avec $`L_j^*L_j`$, la partie dissipative de $`\mathcal L^*(X)`$
s'annule, et $`0=\mathcal L^*(X)=i[H,X]`$. Par hypothèse, $`X\in\mathbb C1`$ ; réciproquement, $`\mathcal L^*(1)=0`$. Donc
$`\ker\mathcal L^*=\mathbb C1`$. Dans les bases $`(E_{ij})`$ et $`(E_{ji})`$ de $`M_n`$, qui sont duales pour l'accouplement
$`(\rho,A)\mapsto\mathrm{Tr}(\rho A)`$, la matrice de $`\mathcal L^*`$ est la transposée de celle de $`\mathcal L`$ ; elles ont donc le
même rang, et $`\dim\ker\mathcal L=\dim\ker\mathcal L^*=1`$. Un état de $`\mathbb C\sigma`$ de trace $`1`$ est $`\sigma`$. $`\square`$

C'est le cas de dimension finie d'un critère de Frigerio [*Commun. Math. Phys.* 63, 269 (1978)],
utilisé à la section 8.5. Les deux hypothèses sont nécessaires. Dans l'exercice 4 du chapitre 8, l'état $`\tfrac12\,1`$ est fidèle
et stationnaire, mais $`H`$ et $`\sigma_z`$ commutent avec toutes les matrices diagonales, et les états stationnaires engendrent
l'espace de dimension deux des matrices diagonales. Inversement, sur $`M_3`$, soient $`H=0`$, $`L_1=|e_1\rangle\langle e_3|`$,
$`L_2=|e_2\rangle\langle e_3|`$. La $`*`$-algèbre engendrée par $`L_1,L_2`$ contient $`L_1L_2^*=E_{12}`$, $`L_1^*L_1=E_{33}`$ et donc
toutes les unités matricielles, de sorte que seuls les multiples de $`1`$ commutent avec les $`L_j,L_j^*`$ ; mais tout $`\rho=P\rho P`$ avec $`P=E_{11}+E_{22}`$
est stationnaire, puisque $`L_j\rho=0=\rho L_j^*`$ et $`\{E_{33},\rho\}=0`$, donc $`\dim\ker\mathcal L\ge4`$ ; d'après le théorème D.6, aucun état
stationnaire n'est fidèle.

*Calcul :* le théorème D.6 pour le générateur de la section D.2, et le second exemple — [tutoriel, section D.3](../../tutorial/fr/D_gksl.md#s-D-3).

<a id="s-D-4"></a>

## D.4 Algèbres à secteurs

Soit $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}\subseteq M_D`$ l'algèbre des matrices diagonales par blocs, $`D=\sum_ck_c`$, de
projecteurs centraux minimaux $`P_c`$, et soit $`E(X)=\sum_cP_cXP_c`$ l'espérance conditionnelle de $`M_D`$ sur
$`\mathcal M`$. Comme dans la démonstration de la proposition 8.8, $`E(B)=B`$, $`E(XB)=E(X)B`$ et $`E(BX)=BE(X)`$ pour $`B\in\mathcal M`$,
et $`E(X^*)=E(X)^*`$. Pour une application CP $`\Phi:\mathcal M\to\mathcal M`$, $`\Phi^*`$ désigne sa duale brute (proposition 6.8).

<a id="r-D-7"></a>

**Théorème D.7.** Soit $`\mathcal L`$ le générateur d'un semi-groupe dynamique quantique $`(T_t)`$ sur $`\mathcal M`$. Il existe
une application CP $`\Phi:\mathcal M\to\mathcal M`$ et $`H=H^*\in\mathcal M`$ tels que
```math
\mathcal L(\rho)=\Phi(\rho)-\tfrac12\{\Phi^*(1),\rho\}-i[H,\rho]\qquad(\rho\in\mathcal M).
```
D'après la proposition 6.4, $`\Phi(\rho)_d=\sum_{c,l}K^{dc}_l\rho_cK^{dc*}_l`$ avec des opérateurs de Kraus entre les secteurs ;
réciproquement, tout $`\mathcal L`$ de cette forme engendre un semi-groupe dynamique quantique sur $`\mathcal M`$ (proposition 8.8(1)).
Ainsi les générateurs de la section 8.3 (où $`\Phi=\sum_jg_jL_j\cdot L_j^*+\sum_mg_m\Phi_m`$) sont exactement les
générateurs des semi-groupes dynamiques quantiques sur $`\mathcal M`$.

*Démonstration.* *Étape 1 : une application CP sur $`M_D`$.* D'après la proposition 6.4,
$`T_t(\rho)_d=\sum_{c,l}K^{dc}_l(t)\rho_cK^{dc}_l(t)^*`$. Soit $`N^{dc}_l(t)\in M_D`$ la matrice ayant le bloc $`K^{dc}_l(t)`$ en
position $`(d,c)`$ et des zéros ailleurs. Pour $`X\in M_D`$, $`N^{dc}_lXN^{dc*}_l`$ est le bloc
$`K^{dc}_lX_{cc}K^{dc*}_l`$ en position $`(d,d)`$, où $`X_{cc}`$ est le $`c`$-ième bloc diagonal, c'est-à-dire le $`c`$-ième bloc
de $`E(X)`$. Par conséquent, $`\hat T_t=T_t\circ E=\sum_{d,c,l}N^{dc}_l\cdot N^{dc*}_l`$ est CP sur $`M_D`$ (théorème 6.2), avec
$`\hat T_0=E`$. D'après le lemme 8.5, $`(\hat T_t-E)/t=\big((T_t-\mathrm{id})/t\big)\circ E\to\mathcal L\circ E`$.

*Étape 2 : positivité hors des vecteurs des secteurs.* Soient $`p_c=(1\otimes P_c)\Omega\in\mathbb C^D\otimes\mathbb C^D`$ et
$`W=\mathrm{span}\{p_c\}`$. D'après le lemme D.1(c), $`J(E)=\sum_c|p_c\rangle\langle p_c|`$, donc pour $`v\perp W`$ et $`t>0`$,
$`0\le\frac1t\langle v,J(\hat T_t)v\rangle=\langle v,J((\hat T_t-E)/t)v\rangle\to\langle v,J(\mathcal L\circ E)v\rangle`$.
Tout $`T_t`$ préserve l'hermiticité (par sa forme de Kraus), donc $`\mathcal L`$ aussi (comme dans le théorème D.3), ainsi que $`E`$
et $`\mathcal L\circ E`$.

*Étape 3 : décomposition.* D'après le lemme D.1(b), $`\langle p_c,p_d\rangle=\mathrm{Tr}(P_cP_d)=\delta_{cd}k_c`$, donc les
$`q_c=p_c/\sqrt{k_c}=(1\otimes Q_c)\Omega`$, avec $`Q_c=P_c/\sqrt{k_c}`$, sont orthonormés. Le lemme D.2 appliqué à
$`\Psi=\mathcal L\circ E`$ fournit une application CP $`\hat\Phi=\sum_j\hat L_j\cdot\hat L_j^*`$ sur $`M_D`$ et des $`K_c\in M_D`$ tels que
```math
\mathcal L(E(A))=\hat\Phi(A)+\sum_c\big(K_cAP_c+P_cAK_c^*\big)/\sqrt{k_c}\qquad(A\in M_D).
```
Pour $`A\in\mathcal M`$, $`E(A)=A`$ et $`AP_c=P_cA`$, donc la somme vaut $`\tilde GA+A\tilde G^*`$ avec
$`\tilde G=\sum_cK_cP_c/\sqrt{k_c}`$.

*Étape 4 : compression.* Pour $`A\in\mathcal M`$, $`\mathcal L(A)\in\mathcal M`$, donc
$`\mathcal L(A)=E(\mathcal L(A))=E(\hat\Phi(A))+GA+AG^*`$ avec $`G=E(\tilde G)\in\mathcal M`$. L'application
$`\Phi=E\circ\hat\Phi`$ sur $`\mathcal M`$ s'écrit
$`\Phi(A)_d=\sum_{c,j}(P_d\hat L_jP_c)A_c(P_d\hat L_jP_c)^*`$ : elle est CP d'après la proposition 6.4, et ses opérateurs de Kraus sont les
blocs des $`\hat L_j`$.

*Étape 5 : préservation de la trace.* Comme dans le théorème D.4, $`\mathrm{Tr}\,\mathcal L(A)=0`$ pour $`A\in\mathcal M`$, et
$`\mathrm{Tr}\,\Phi(A)=\mathrm{Tr}(A\Phi^*(1))`$, donc $`\mathrm{Tr}(AX)=0`$ pour tout $`A\in\mathcal M`$, où
$`X=\Phi^*(1)+G+G^*\in\mathcal M`$ est autoadjoint. Avec $`A=X`$, on obtient $`X=0`$ (la trace brute est fidèle, théorème 3.1).
Avec $`H=\frac i2(G-G^*)\in\mathcal M`$, on a $`G=-iH-\frac12\Phi^*(1)`$ et $`G\rho+\rho G^*=-i[H,\rho]-\frac12\{\Phi^*(1),\rho\}`$. $`\square`$

Les vecteurs $`p_c`$ remplacent $`\Omega=\sum_cp_c`$ parce que $`\hat T_0=E`$ et non $`\mathrm{id}`$ : $`J(\mathcal L\circ E)`$ n'est pas nécessairement
positive sur $`\Omega^\perp`$, et dans l'exemple ci-dessous elle ne l'est pas.

*Calcul :* le théorème D.7 sur $`M_1\oplus M_2\subseteq M_3`$, en suivant les cinq étapes de la démonstration ; $`J(\mathcal L\circ E)`$ n'y est pas positive sur $`\Omega^\perp`$ — [tutoriel, section D.4](../../tutorial/fr/D_gksl.md#s-D-4). *Dans la bibliothèque :* les générateurs sur une algèbre à secteurs, `dynamics.lindblad_rhs` avec `channels.random_inter_sector_channel` — [tutoriel, section 8.3](../../tutorial/fr/08_dynamics.md#s-8-3).

<a id="r-D-8"></a>

**Corollaire D.8 (représentation de Heisenberg).** Soit $`(S_t)_{t\ge0}`$ une famille continue d'applications CP unifères sur
$`\mathcal M`$ avec $`S_0=\mathrm{id}`$ et $`S_{s+t}=S_sS_t`$. Son générateur est $`A\mapsto\Psi(A)+K^*A+AK`$ avec
$`\Psi`$ CP sur $`\mathcal M`$, $`K\in\mathcal M`$ et $`\Psi(1)+K+K^*=0`$.

*Démonstration.* Soient $`T_t=S_t^*`$ (duales brutes). Elles sont CP (proposition 6.8) et TP (corollaire 6.9(b)),
$`T_0=\mathrm{id}`$, et $`T_{s+t}=(S_tS_s)^*=S_s^*S_t^*=T_sT_t`$, car
$`\mathrm{Tr}(S_tS_s(x)y)=\mathrm{Tr}(xS_s^*S_t^*(y))`$ et $`S_tS_s=S_{s+t}`$. Comme $`\Psi\mapsto\Psi^*`$ est linéaire sur un
espace de dimension finie, $`t\mapsto T_t`$ est continue et $`(T_t-\mathrm{id})/t\to\mathcal L_S^*`$, où $`\mathcal L_S`$ est le
générateur de $`(S_t)`$. Le théorème D.7 donne $`\mathcal L_S^*(\rho)=\Phi(\rho)+G\rho+\rho G^*`$ avec $`G+G^*=-\Phi^*(1)`$, et par
cyclicité, $`\mathrm{Tr}((G\rho+\rho G^*)A)=\mathrm{Tr}(\rho(G^*A+AG))`$, donc d'après le corollaire 6.9(c),
$`\mathcal L_S(A)=\Phi^*(A)+G^*A+AG`$. On prend $`\Psi=\Phi^*`$ (CP d'après la proposition 6.8) et $`K=G`$. $`\square`$

C'est le cas de dimension finie du théorème de Christensen et Evans cité après le théorème 8.7.

*Calcul :* le générateur du corollaire D.8 comme dual de $`\mathcal L`$ — [tutoriel, section D.4](../../tutorial/fr/D_gksl.md#s-D-4).

<a id="s-D-5"></a>

## D.5 Correction d'erreurs par algèbres d'opérateurs

Nous reprenons le cadre de la section 10.3. $`\Phi(\rho)=\sum_{i=1}^rE_i\rho E_i^*`$ est un canal sur $`M_D`$, $`\mathcal C=P\mathbb C^D`$
un code, $`B(\mathcal C)=\{X\in M_D:X=PXP\}`$ l'ensemble des opérateurs sur $`\mathcal C`$, $`\mathfrak A\subseteq B(\mathcal C)`$ une $`*`$-algèbre
d'unité $`P`$ (une algèbre de von Neumann sur $`\mathcal C`$ d'après le théorème 2.2), et $`\mathfrak A'`$ son commutant dans
$`B(\mathcal C)`$. $`\mathfrak A`$ est **corrigeable** pour $`\Phi`$ s'il existe un canal $`\mathcal R`$ sur $`M_D`$ tel que
$`P(\mathcal R\circ\Phi)^*(A)P=A`$ pour tout $`A\in\mathfrak A`$. Si $`\mathcal R`$ a pour opérateurs de Kraus les $`R_m`$, alors
$`\mathcal R\circ\Phi`$ a pour opérateurs de Kraus les $`R_mE_i`$, et $`\sum_mR_m^*R_m=1`$.

<a id="r-D-9"></a>

**Théorème D.9 (Bény–Kempf–Kribs ; théorème 10.7).** $`\mathfrak A`$ est corrigeable pour $`\Phi`$ si et seulement si
$`PE_i^*E_jP\in\mathfrak A'`$ pour tous $`i,j`$. Dans ce cas, le canal $`\mathcal R`$ construit dans la démonstration corrige $`\mathfrak A`$.

*Démonstration.* *Condition nécessaire.* Supposons que $`\mathcal R`$ corrige $`\mathfrak A`$ et posons $`W_{mi}=R_mE_iP`$, de sorte que
$`\sum_{m,i}W_{mi}^*XW_{mi}=P(\mathcal R\circ\Phi)^*(X)P`$ pour tout $`X\in M_D`$. Ceci vaut $`X`$ pour $`X\in\mathfrak A`$, et $`P`$ pour
$`X=1`$, car $`(\mathcal R\circ\Phi)^*`$ est unifère (corollaire 6.9(a)). Soit $`p\in\mathfrak A`$ un projecteur ; alors $`p\le P`$
et $`\sum W_{mi}^*(1-p)W_{mi}=P-p`$. Pour $`x=px`$,
$`\sum_{m,i}\|(1-p)W_{mi}x\|^2=\langle x,(P-p)x\rangle=0`$, et pour $`x=(P-p)x`$,
$`\sum_{m,i}\|pW_{mi}x\|^2=\langle x,px\rangle=0`$. Par conséquent $`(1-p)W_{mi}p=0`$ et $`pW_{mi}(P-p)=0`$, et comme
$`W_{mi}=W_{mi}P`$,
```math
W_{mi}p=pW_{mi}p=pW_{mi}P=pW_{mi}.
```
Tout $`A\in\mathfrak A`$ autoadjoint est une combinaison linéaire de projecteurs de $`\mathfrak A`$ : d'après le théorème 1.3 sur $`\mathcal C`$,
$`A=\sum_k\lambda_kp_k`$ avec des $`\lambda_k`$ distincts et des projecteurs orthogonaux $`p_k`$, $`\sum_kp_k=P`$, et
$`p_k=\prod_{l\neq k}(A-\lambda_lP)/(\lambda_k-\lambda_l)\in\mathfrak A`$ (appliquer les deux membres aux vecteurs propres dans $`\mathcal C`$ et à $`\mathcal C^\perp`$ ; le produit vide vaut $`P`$). Tout élément
de $`\mathfrak A`$ s'écrit $`A_1+iA_2`$ avec $`A_1,A_2\in\mathfrak A`$ autoadjoints. Ainsi $`W_{mi}A=AW_{mi}`$ pour tout $`A\in\mathfrak A`$ et,
en passant aux adjoints, $`W_{mi}^*A=AW_{mi}^*`$. Comme $`\sum_mW_{mi}^*W_{mj}=PE_i^*(\sum_mR_m^*R_m)E_jP=PE_i^*E_jP`$,
```math
A\,PE_i^*E_jP=\sum_mW_{mi}^*AW_{mj}=PE_i^*E_jP\,A\qquad(A\in\mathfrak A),
```
c'est-à-dire $`PE_i^*E_jP\in\mathfrak A'`$.

*Condition suffisante.* D'après le théorème 2.6 appliqué à $`\mathfrak A`$ sur $`\mathcal C`$, il existe une base orthonormée $`|c;a,b\rangle`$
($`c=1,\dots,C`$, $`a=1,\dots,k_c`$, $`b=1,\dots,m_c`$) de $`\mathcal C`$ — les images des $`e_a\otimes f_b`$ du $`c`$-ième terme par
l'unitaire $`V`$ — dans laquelle $`\mathfrak A`$ est formée des opérateurs $`\sum_{c,a,a',b}(A_c)_{aa'}|c;a,b\rangle\langle c;a',b|`$,
$`A_c\in M_{k_c}`$, et $`\mathfrak A'`$ des opérateurs $`\sum_{c,a,b,b'}(Y_c)_{bb'}|c;a,b\rangle\langle c;a,b'|`$, $`Y_c\in M_{m_c}`$.
L'hypothèse fournit donc des matrices $`\Lambda^c_{ij}\in M_{m_c}`$ telles que
```math
\langle c;a,b|E_i^*E_j|d;a',b'\rangle=\delta_{cd}\,\delta_{aa'}\,(\Lambda^c_{ij})_{bb'} .
```
*Étape 1.* Pour chaque $`c`$ et chaque couple $`\alpha=(i,b)`$, soient $`F^c_\alpha:\mathbb C^{k_c}\to\mathbb C^D`$,
$`F^c_\alpha e_a=E_i|c;a,b\rangle`$, et $`\Lambda^c_{\alpha\beta}=(\Lambda^c_{ij})_{bb'}`$ pour $`\beta=(j,b')`$. La formule ci-dessus s'écrit
$`F^{c*}_\alpha F^d_\beta=\delta_{cd}\Lambda^c_{\alpha\beta}1_{k_c}`$ : une condition de Knill–Laflamme pour les « erreurs » $`F^c_\alpha`$ sur le
code $`\mathbb C^{k_c}`$, avec en outre l'orthogonalité entre les secteurs.

*Étape 2.* Comme dans la démonstration du théorème 10.5, $`\Lambda^c`$ est autoadjointe, car
$`\Lambda^c_{\beta\alpha}1=F^{c*}_\beta F^c_\alpha=(F^{c*}_\alpha F^c_\beta)^*=\bar\Lambda^c_{\alpha\beta}1`$, et positive, car
$`\sum\bar z_\alpha z_\beta\Lambda^c_{\alpha\beta}1=(\sum z_\alpha F^c_\alpha)^*(\sum z_\beta F^c_\beta)\ge0`$. Écrivons
$`\Lambda^c=U^c\mathrm{diag}(d^c_\kappa)U^{c*}`$ avec $`U^c`$ unitaire et $`d^c_\kappa\ge0`$, et
$`G^c_\kappa=\sum_\alpha u^c_{\alpha\kappa}F^c_\alpha`$. Alors $`G^{c*}_\kappa G^d_\lambda=\delta_{cd}\delta_{\kappa\lambda}d^c_\kappa1`$, donc
$`G^c_\kappa=0`$ si $`d^c_\kappa=0`$, et sinon $`V^c_\kappa=G^c_\kappa/\sqrt{d^c_\kappa}`$ est une isométrie avec
$`V^{c*}_\kappa V^d_\lambda=\delta_{cd}\delta_{\kappa\lambda}1`$. En inversant, $`F^c_\alpha=\sum_\kappa\bar u^c_{\alpha\kappa}G^c_\kappa =\sum_{\kappa:d^c_\kappa>0}\bar u^c_{\alpha\kappa}\sqrt{d^c_\kappa}\,V^c_\kappa`$, d'où
```math
V^{c'*}_\kappa F^d_\beta=\delta_{c'd}\,\bar u^d_{\beta\kappa}\sqrt{d^d_\kappa}\,1,\qquad
F^{c*}_\alpha V^{c'}_\kappa=\delta_{cc'}\,u^c_{\alpha\kappa}\sqrt{d^c_\kappa}\,1 .
```
*Étape 3 : la récupération.* Les $`V^c_\kappa V^{c*}_\kappa`$ ($`d^c_\kappa>0`$) sont des projecteurs deux à deux orthogonaux ; soit $`Q`$ leur
somme. Soient $`\iota_c:\mathbb C^{k_c}\to\mathbb C^D`$, $`\iota_ce_a=|c;a,1\rangle`$, une isométrie, et
```math
\mathcal R(\sigma)=\sum_{c}\sum_{\kappa:d^c_\kappa>0}\iota_cV^{c*}_\kappa\,\sigma\,V^c_\kappa\iota_c^*+(1-Q)\sigma(1-Q).
```
C'est un canal : $`\sum_{c,\kappa}V^c_\kappa\iota_c^*\iota_cV^{c*}_\kappa+(1-Q)^2=Q+1-Q=1`$. Il mesure le secteur $`c`$ et
le syndrome $`\kappa`$, annule l'erreur sur $`\mathbb C^{k_c}`$ et écrit le résultat dans la copie $`b=1`$.

*Étape 4 : correction.* $`\mathcal R^*(X)=\sum_{c,\kappa}V^c_\kappa\iota_c^*X\iota_cV^{c*}_\kappa+(1-Q)X(1-Q)`$, et pour
$`A\in\mathfrak A`$ de blocs $`A_c`$, $`\iota_c^*A\iota_c=A_c`$. Comme chaque $`F^d_\beta`$ est une combinaison des $`V^d_\kappa`$,
$`(1-Q)F^d_\beta=0`$. Avec l'étape 2, pour $`\alpha=(i,b)`$, $`\beta=(j,b')`$,
```math
\langle c;a,b|E_i^*\mathcal R^*(A)E_j|d;a',b'\rangle=\sum_{c',\kappa}\langle e_a,F^{c*}_\alpha V^{c'}_\kappa A_{c'}V^{c'*}_\kappa F^d_\beta e_{a'}\rangle
=\delta_{cd}\sum_\kappa u^c_{\alpha\kappa}d^c_\kappa\bar u^c_{\beta\kappa}(A_c)_{aa'}=\delta_{cd}(\Lambda^c_{ij})_{bb'}(A_c)_{aa'} .
```
Prenons $`j=i`$ et sommons sur $`i`$ : d'après la formule de l'hypothèse avec $`\sum_iE_i^*E_i=1`$, $`\sum_i\Lambda^c_{ii}=1_{m_c}`$, donc
```math
\langle c;a,b|(\mathcal R\circ\Phi)^*(A)|d;a',b'\rangle=\delta_{cd}\,\delta_{bb'}(A_c)_{aa'}=\langle c;a,b|A|d;a',b'\rangle .
```
$`P(\mathcal R\circ\Phi)^*(A)P`$ et $`A`$ appartiennent tous deux à $`B(\mathcal C)`$ et ont la même matrice dans la base $`|c;a,b\rangle`$ ; ils sont donc
égaux. $`\square`$

*Calcul :* le théorème D.9 et la récupération de sa démonstration sur un code hybride dans $`\mathbb C^8`$ — [tutoriel, section D.5](../../tutorial/fr/D_gksl.md#s-D-5).

<a id="r-D-10"></a>

**Corollaire D.10.** (a) Pour $`\mathfrak A=B(\mathcal C)`$, $`\mathfrak A'=\mathbb CP`$, et le théorème D.9 est le théorème de
Knill–Laflamme 10.5 ($`C=1`$, $`m_1=1`$ ; l'étape 3 est la récupération du théorème 10.5). (b) Pour $`\mathfrak A=\mathrm{span}\{P_1,\dots,P_C\}`$ avec
des projecteurs orthogonaux $`P_c\ne0`$, $`\sum_cP_c=P`$, c'est l'équivalence de 2 et 3 dans la proposition 10.6. (c) (*Sous-systèmes
sans bruit.*) Avec les notations de la proposition 10.1(a), soient $`\mathcal C=V(\mathbb C^{k_c}\otimes\mathbb C^{m_c})`$ le $`c`$-ième
terme pour l'algèbre d'interaction $`\mathcal A`$ de $`\Phi`$, et $`\mathfrak A=V(1_{k_c}\otimes M_{m_c})V^*`$ sur $`\mathcal C`$, une
algèbre de la forme $`M_{m}\otimes1_{k}`$ après échange des deux facteurs. Alors $`\mathfrak A`$ est corrigée par
$`\mathcal R=\mathrm{id}`$, et $`PE_i^*E_jP=V(E_{i,c}^*E_{j,c}\otimes1_{m_c})V^*\in\mathfrak A'`$ ($`P=VV^*`$).

*Démonstration.* (a) $`\mathfrak A'=\mathbb CP`$ est le cas $`C=1`$, $`m_1=1`$ de la description de $`\mathfrak A'`$ ci-dessus, et
$`PE_i^*E_jP\in\mathbb CP`$ est la condition de Knill–Laflamme. La corrigeabilité de $`B(\mathcal C)`$ est la correction au sens de la
section 10.2 : si $`\mathcal R(\Phi(\rho))=\rho`$ pour tous les états $`\rho`$ à support dans $`\mathcal C`$, alors
$`\mathrm{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(A)P)=\mathrm{Tr}(\rho A)`$ pour chacun d'eux, et ces $`\rho`$ engendrent $`B(\mathcal C)`$ ; la
réciproque est montrée à la section 10.3. Avec $`C=1`$, $`m_1=1`$, les opérateurs de Kraus $`\iota_1V^{1*}_\kappa`$ de l'étape 3 sont ceux
du théorème 10.5, car $`\iota_1\iota_1^*=P`$. (b) $`T\in B(\mathcal C)`$ commute avec tous les $`P_c`$ si et seulement si
$`P_cTP_d=0`$ pour $`c\ne d`$, ce qui est la condition 2 de la proposition 10.6. (c) Tout $`A\in\mathfrak A`$, prolongé par $`0`$ aux
autres termes, appartient à $`\mathcal A'=V(\bigoplus_c1_{k_c}\otimes M_{m_c})V^*`$ (théorème 2.6), donc $`\Phi^*(A)=A`$ d'après la proposition
10.1(a), et $`PAP=A`$. La formule pour $`PE_i^*E_jP`$ découle de $`V^*E_iV=\bigoplus_cE_{i,c}\otimes1_{m_c}`$, et
$`E_{i,c}^*E_{j,c}\otimes1_{m_c}`$ commute avec $`1_{k_c}\otimes M_{m_c}`$. $`\square`$

D'après le théorème 2.6, toute $`\mathfrak A`$ est $`*`$-isomorphe à une algèbre $`\bigoplus_cM_{k_c}`$. Lorsque $`C\ge2`$ et que $`k_c\ge2`$ pour au moins un $`c`$, une
$`\mathfrak A`$ corrigeable protège l'étiquette de secteur $`c`$ (information classique) en même temps que l'information quantique de
chaque $`M_{k_c}`$ ; de tels codes sont appelés codes **hybrides**.

Les résultats de cette annexe sont dus à V. Gorini, A. Kossakowski et E. C. G. Sudarshan
[*J. Math. Phys.* 17, 821 (1976)] et G. Lindblad [*Commun. Math. Phys.* 48, 119 (1976)] (théorèmes D.3–D.5),
E. Christensen et D. E. Evans [*J. London Math. Soc.* 20, 358 (1979)] (corollaire D.8 dans des algèbres de von Neumann
générales), A. Frigerio (théorème D.6), et C. Bény, A. Kempf et D. W. Kribs [*Phys. Rev. Lett.* 98, 100502
(2007) ; *Phys. Rev. A* 76, 042303 (2007)] (théorème D.9) ; les démonstrations ci-dessus sont complètes pour les énoncés
de dimension finie utilisés dans le cours.
