[English](../06_quantum_channels.md) · [Español](../es/06_quantum_channels.md) · **Français** · [Deutsch](../de/06_quantum_channels.md) · [中文](../zh/06_quantum_channels.md) · [日本語](../ja/06_quantum_channels.md) · [Русский](../ru/06_quantum_channels.md)

# 6. Canaux quantiques

Un canal quantique est une application linéaire complètement positive qui préserve la trace. La
conjugaison par une unitaire, l'espérance conditionnelle sur le centre, l'application
$`\rho\mapsto\mathrm{Tr}_E(V\rho V^*)`$ donnée par une isométrie $`V`$ et les applications qui
transfèrent du poids d'un secteur à un autre en sont des exemples ; leur lecture physique (une
évolution, une mesure dont on oublie le résultat, un environnement, la perte d'une particule) n'est
discutée que dans les encadrés *Interprétation physique*. Ce chapitre démontre les descriptions standard de ces applications — matrices de
Choi, opérateurs de Kraus, dilatation de Stinespring — pour les applications entre algèbres de von
Neumann de dimension finie, puis étudie ce que les sommes directes ont de particulier : les
applications entre secteurs, les duales relativement aux traces du chapitre 3, et la partie
classique d'un canal. Les chapitres 2 et 3 ont déjà utilisé les canaux de la bibliothèque ; en voici
la théorie.

**Vous apprendrez** pourquoi la positivité ne suffit pas, comment passer des matrices de Choi aux
opérateurs de Kraus, aux superopérateurs et aux isométries et inversement, pourquoi « préserver la
trace » dépend de la trace dès qu'une application change de secteur, quels opérateurs un canal unifère
laisse fixes, et ce qu'un canal aléatoire fait en moyenne à un état pur.

Les vérifications numériques de ce chapitre se trouvent dans le [Tutoriel 6](../../tutorial/fr/06_quantum_channels.md).

Dans tout ce chapitre, $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$ et
$`\mathcal N=\bigoplus_{d=1}^{C'}M_{k'_d}(\mathbb C)`$, de projecteurs centraux minimaux $`P_c`$ et
$`Q_d`$ ; $`A=\bigoplus_cA_c`$ ; $`M_{k'\times k}`$ est l'espace des matrices complexes $`k'\times k`$ ;
les $`E_{ij}=|i\rangle\langle j|`$ sont les unités matricielles ; $`\mathrm{Tr}`$ est la trace brute
$`\sum_c\mathrm{Tr}\,A_c`$, et $`\mathrm{tr}_w=\sum_cw_c\mathrm{Tr}\,A_c`$ une trace fidèle
(théorème 3.1). Les matrices par blocs agissent sur $`\mathbb C^n\otimes\mathbb C^k`$, l'indice de
bloc venant en premier.

<a id="s-6-1"></a>

## 6.1 Applications complètement positives et théorème de Choi

Une application linéaire $`\Phi:\mathcal M\to\mathcal N`$ est **positive** si $`A\ge0`$ entraîne
$`\Phi(A)\ge0`$. Les matrices $`n\times n`$ $`[A_{ij}]`$ à coefficients dans $`\mathcal M`$ forment à nouveau
une somme directe d'algèbres de matrices, $`M_n(\mathcal M)\cong\bigoplus_cM_n(M_{k_c})`$ via
$`[A_{ij}]\mapsto([(A_{ij})_c]_{ij})_c`$, de sorte que $`[A_{ij}]\ge0`$ signifie que chaque matrice par
blocs $`[(A_{ij})_c]\in M_{nk_c}`$ est positive. $`\Phi`$ est **$`n`$-positive** si
$`\Phi_n:[A_{ij}]\mapsto[\Phi(A_{ij})]`$ est positive, et **complètement positive** (CP) si elle est
$`n`$-positive pour tout $`n`$. Elle **préserve la trace** (TP) si
$`\mathrm{Tr}\,\Phi(A)=\mathrm{Tr}\,A`$, elle est **unifère** si $`\Phi(1)=1`$, et c'est un
**canal** si elle est CP et TP. Nous utiliserons constamment le fait que $`X\ge0`$ entraîne
$`TXT^*\ge0`$ pour toute matrice rectangulaire $`T`$, puisque
$`\langle y,TXT^*y\rangle=\langle T^*y,XT^*y\rangle`$.

<a id="r-6-1"></a>

**Proposition 6.1.** Pour $`k\ge2`$, la transposition $`\theta(A)=A^T`$ sur $`M_k`$ est positive mais non
$`2`$-positive.

*Démonstration.* Pour $`A\ge0`$, $`\langle x,A^Tx\rangle=\sum_{ij}\bar x_iA_{ji}x_j=\langle\bar x,A\bar x\rangle\ge0`$.
Avec $`\Omega=e_1\otimes e_1+e_2\otimes e_2\in\mathbb C^2\otimes\mathbb C^k`$,
$`|\Omega\rangle\langle\Omega|=[E_{ij}]_{i,j=1}^2\ge0`$, et $`\theta_2`$ l'envoie sur
$`F=[E_{ji}]_{i,j=1}^2`$, qui envoie $`e_1\otimes e_2-e_2\otimes e_1`$ sur son opposé. $`\square`$

*Calcul :* la matrice de Choi de la transposition a la valeur propre $`-1`$ — [tutoriel, section 6.1](../../tutorial/fr/06_quantum_channels.md#s-6-1).

<a id="r-6-2"></a>

**Théorème 6.2 (Choi).** Pour une application linéaire $`\Phi:M_k\to M_{k'}`$, les assertions suivantes
sont équivalentes : (a) $`\Phi`$ est CP ; (b) $`\Phi`$ est $`k`$-positive ; (c) la **matrice de Choi**
```math
J(\Phi)=\sum_{i,j=1}^kE_{ij}\otimes\Phi(E_{ij})=[\Phi(E_{ij})]_{i,j=1}^k\in M_k(M_{k'})
```
est positive ; (d) il existe $`K_1,\dots,K_r\in M_{k'\times k}`$ (**opérateurs de Kraus**) tels que
$`\Phi(A)=\sum_lK_lAK_l^*`$, et l'on peut prendre $`r=\mathrm{rank}\,J(\Phi)\le kk'`$.

*Démonstration.* (a)$`\Rightarrow`$(b) est trivial. (b)$`\Rightarrow`$(c) : pour
$`\Omega=\sum_ie_i\otimes e_i\in\mathbb C^k\otimes\mathbb C^k`$,
$`|\Omega\rangle\langle\Omega|=[E_{ij}]\ge0`$ et $`J(\Phi)=\Phi_k([E_{ij}])`$. (c)$`\Rightarrow`$(d) :
d'après le théorème 1.3, $`J(\Phi)=\sum_{l=1}^rv_lv_l^*`$ avec $`r=\mathrm{rank}\,J(\Phi)`$ (les
vecteurs propres des valeurs propres non nulles, multipliés par leurs racines carrées). Tout
$`v\in\mathbb C^k\otimes\mathbb C^{k'}`$ s'écrit $`v=\sum_ie_i\otimes Ke_i`$ pour exactement une
$`K\in M_{k'\times k}`$ ($`Ke_i`$ est le $`i`$-ième bloc de $`v`$), et le bloc $`(i,j)`$ de $`vv^*`$ est
$`Ke_ie_j^*K^*=KE_{ij}K^*`$. En écrivant ainsi chaque $`v_l`$ avec $`K_l`$ et en comparant les blocs, on
obtient $`\Phi(E_{ij})=\sum_lK_lE_{ij}K_l^*`$, d'où (d) par linéarité.
(d)$`\Rightarrow`$(a) : $`\Phi_n([A_{ij}])=\sum_l(1_n\otimes K_l)[A_{ij}](1_n\otimes K_l)^*\ge0`$. $`\square`$

<a id="r-6-3"></a>

**Corollaire 6.3.** Si $`\Phi(A)=\sum_{l=1}^sL_lAL_l^*`$, alors $`s\ge\mathrm{rank}\,J(\Phi)`$. Ainsi, le
nombre minimal d'opérateurs de Kraus, le **rang de Kraus**, est $`\mathrm{rank}\,J(\Phi)`$.

*Démonstration.* En lisant le dernier calcul à rebours, $`J(\Phi)=\sum_{l=1}^sw_lw_l^*`$ avec
$`w_l=\sum_ie_i\otimes L_le_i`$, somme de $`s`$ matrices de rang $`\le1`$. $`\square`$

*Calcul :* la matrice de Choi et le rang de Kraus d'un canal aléatoire — [tutoriel, section 6.1](../../tutorial/fr/06_quantum_channels.md#s-6-1).

Une application entre sommes directes est une matrice d'applications entre blocs. Soient
$`\iota_c:M_{k_c}\to\mathcal M`$ l'application qui place une matrice dans le bloc $`c`$,
$`\pi_d:\mathcal N\to M_{k'_d}`$ celle qui extrait le bloc $`d`$, et
$`\Phi_{dc}=\pi_d\circ\Phi\circ\iota_c:M_{k_c}\to M_{k'_d}`$. Comme $`A=\sum_c\iota_c(A_c)`$,
$`\Phi(A)_d=\sum_c\Phi_{dc}(A_c)`$.

<a id="r-6-4"></a>

**Proposition 6.4.** $`\Phi:\mathcal M\to\mathcal N`$ est $`n`$-positive si et seulement si chaque
$`\Phi_{dc}`$ l'est. Par conséquent, $`\Phi`$ est CP si et seulement s'il existe
$`K^{dc}_l\in M_{k'_d\times k_c}`$, $`l=1,\dots,r`$, tels que
```math
\Phi(A)_d=\sum_c\sum_{l=1}^rK^{dc}_l\,A_c\,K^{dc*}_l ,
```
et $`r=\max_{d,c}\mathrm{rank}\,J(\Phi_{dc})`$ suffit. $`\Phi`$ est TP si et seulement si
$`\sum_{d,l}K^{dc*}_lK^{dc}_l=1_{k_c}`$ pour tout $`c`$.

*Démonstration.* Si toutes les $`\Phi_{dc}`$ sont $`n`$-positives et $`[A_{ij}]\ge0`$, alors chaque
$`X^c=[(A_{ij})_c]\ge0`$, et le bloc $`d`$ de $`\Phi_n([A_{ij}])`$ est $`\sum_c(\Phi_{dc})_n(X^c)\ge0`$.
Réciproquement, pour $`[X_{ij}]\ge0`$ dans $`M_n(M_{k_c})`$, on a $`[\iota_c(X_{ij})]\ge0`$ dans
$`M_n(\mathcal M)`$ (son bloc $`c`$ est $`[X_{ij}]`$, les autres sont nuls), et le bloc $`d`$ de son image est
$`[\Phi_{dc}(X_{ij})]`$. La forme de Kraus découle du théorème 6.2 appliqué à chaque couple $`(d,c)`$,
complétée par des opérateurs nuls jusqu'à un $`r`$ commun. Enfin, par cyclicité,
$`\mathrm{Tr}\,\Phi(A)=\sum_c\mathrm{Tr}\big(A_c\sum_{d,l}K^{dc*}_lK^{dc}_l\big)`$, qui vaut
$`\sum_c\mathrm{Tr}\,A_c`$ pour tout $`A`$ si et seulement si chaque somme vaut $`1_{k_c}`$ (prendre des
unités matricielles : $`\mathrm{Tr}(E_{ij}X)=X_{ji}`$). $`\square`$

Une application linéaire $`\Phi:M_k\to M_{k'}`$ vérifie $`\Phi(X^*)=\Phi(X)^*`$ pour tout $`X`$ si et seulement
si $`J(\Phi)^*=J(\Phi)`$ : le bloc $`(i,j)`$ de $`J^*`$ est $`\Phi(E_{ji})^*`$, qui est égal à $`\Phi(E_{ij})`$ pour
tous $`i,j`$ si et seulement si $`\Phi`$ préserve les adjoints, par linéarité.

*Dans la bibliothèque :* les applications qui préservent les secteurs sont des `Channel`, les applications générales de la proposition 6.4 des `InterSectorChannel` ; `Channel.choi()`, `Channel.superoperator()` et `Channel.from_superoperator` font passer d'une représentation à l'autre — [tutoriel, section 6.1](../../tutorial/fr/06_quantum_channels.md#s-6-1).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Un système n'est jamais seul : il peut être corrélé à une ancilla
> non touchée à $`n`$ niveaux. Un processus physique agit comme $`\Phi_n=\mathrm{id}_n\otimes\Phi`$ sur la
> paire et doit envoyer des états sur des états ; $`\Phi`$ doit donc être complètement positive. La
> transposition est positive, mais elle envoie l'état intriqué $`|\Omega\rangle\langle\Omega|/2`$ sur un
> opérateur de valeur propre $`-1/2`$. La matrice de Choi est l'image de l'état maximalement intriqué :
> un seul état décide de la complète positivité.

<a id="s-6-2"></a>

## 6.2 Dilatation de Stinespring et liberté dans le choix des opérateurs de Kraus

Pour $`\Phi(A)=\sum_{l=1}^rK_lAK_l^*`$ de $`M_k`$ dans $`M_{k'}`$, posons $`\Phi^*(B)=\sum_lK_l^*BK_l`$ ; par
cyclicité, $`\mathrm{Tr}(\Phi(A)B)=\mathrm{Tr}(A\Phi^*(B))`$, donc $`\Phi^*`$ est la duale de
$`\Phi`$ relativement à la trace (section 6.3) et ne dépend pas du choix des $`K_l`$ (si
$`\mathrm{Tr}(AX)=\mathrm{Tr}(AX')`$ pour tout $`A`$, alors $`X=X'`$ : prendre $`A=E_{ji}`$). Sur
$`\mathbb C^r\otimes\mathbb C^{k'}`$ (l'**environnement** $`\mathbb C^r`$ en premier), posons
$`\mathrm{Tr}_E\,X=\sum_l(e_l^*\otimes1)X(e_l\otimes1)`$.

<a id="r-6-5"></a>

**Théorème 6.5 (Stinespring, dimension finie).** Soit
$`V=\sum_le_l\otimes K_l:\mathbb C^k\to\mathbb C^r\otimes\mathbb C^{k'}`$. Alors
```math
\Phi(A)=\mathrm{Tr}_E(VAV^*),\qquad\Phi^*(B)=V^*(1_r\otimes B)V,
```
et $`\Phi`$ est TP si et seulement si $`V`$ est une isométrie. Réciproquement, tout $`V`$ définit par la
première formule une application CP.

*Démonstration.* $`(e_l^*\otimes1)V=K_l`$, donc $`\mathrm{Tr}_E(VAV^*)=\sum_lK_lAK_l^*`$ ;
réciproquement, $`\mathrm{Tr}_E(VAV^*)`$ est une forme de Kraus avec $`K_l=(e_l^*\otimes1)V`$.
Ensuite, $`V^*(1\otimes B)V=\sum_{l,l'}\langle e_l,e_{l'}\rangle K_l^*BK_{l'}=\Phi^*(B)`$, et
$`V^*V=\sum_lK_l^*K_l`$, qui vaut $`1`$ si et seulement si $`\Phi`$ est TP (proposition 6.4). $`\square`$

Ainsi, le plus petit $`r`$ pour lequel $`\Phi(A)=\mathrm{Tr}_E(VAV^*)`$ avec un certain
$`V:\mathbb C^k\to\mathbb C^r\otimes\mathbb C^{k'}`$ est le rang de Kraus $`\mathrm{rank}\,J(\Phi)`$ :
tout tel $`V`$ donne une forme de Kraus à $`r`$ opérateurs, donc $`r\ge\mathrm{rank}\,J(\Phi)`$
(corollaire 6.3), et $`r=\mathrm{rank}\,J(\Phi)`$ est atteint d'après le théorème 6.2(d).

*Calcul :* le théorème 6.5 pour un canal aléatoire — [tutoriel, section 6.2](../../tutorial/fr/06_quantum_channels.md#s-6-2).

La formule $`\Phi^*(B)=V^*\pi(B)V`$ avec la représentation $`\pi(B)=1_r\otimes B`$ est la forme sous
laquelle W. F. Stinespring (Proc. AMS 6 (1955) 211) a énoncé le théorème, dans un cadre bien plus
général que le cadre de dimension finie démontré ici (attribution historique ; pour aller plus loin :
V. Paulsen, *Completely Bounded Maps and Operator Algebras*, Cambridge 2002, chap. 4). Pour une
application entre sommes directes, soient $`\iota'_d`$ l'injection du $`d`$-ième facteur de
$`\bigoplus_d\mathbb C^{k'_d}`$ et
$`V_c=\sum_{d,l}e_l\otimes\iota'_dK^{dc}_l:\mathbb C^{k_c}\to\mathbb C^r\otimes\bigoplus_d\mathbb C^{k'_d}`$.
Alors $`(e_l^*\otimes1)V_c=\sum_d\iota'_dK^{dc}_l`$, de sorte que le bloc $`(d,d)`$ de
$`\sum_c\mathrm{Tr}_E(V_cA_cV_c^*)=\sum_{c,l}\big(\sum_d\iota'_dK^{dc}_l\big)A_c\big(\sum_{d'}\iota'_{d'}K^{d'c}_l\big)^*`$
est $`\sum_{c,l}K^{dc}_lA_cK^{dc*}_l=\Phi(A)_d`$ ; et, comme $`\iota'^*_d\iota'_{d'}=\delta_{dd'}1`$,
$`V_c^*V_c=\sum_{d,l}K^{dc*}_lK^{dc}_l`$, donc $`\Phi`$ est TP si et seulement si chaque $`V_c`$ est une
isométrie (proposition 6.4). C'est ainsi que `random_channel` et `random_inter_sector_channel` tirent
des canaux (section 6.5).

*Dans la bibliothèque :* `channels.random_channel` et `channels.random_inter_sector_channel` — [tutoriel, section 6.5](../../tutorial/fr/06_quantum_channels.md#s-6-5).

<a id="r-6-6"></a>

**Lemme 6.6.** Si $`X,Y\in M_{p\times n}`$ et $`XX^*=YY^*`$, alors $`Y=XW`$ pour une unitaire $`W\in M_n`$.

*Démonstration.* $`\|X^*z\|^2=\langle z,XX^*z\rangle=\|Y^*z\|^2`$ pour tout $`z`$, donc $`X^*z\mapsto Y^*z`$
est une isométrie linéaire bien définie de l'image de $`X^*`$ sur l'image de $`Y^*`$. Comme dans la
démonstration du théorème 1.8, elle se prolonge en une unitaire $`U`$ telle que $`UX^*=Y^*`$ ; on prend
$`W=U^*`$. $`\square`$

<a id="r-6-7"></a>

**Théorème 6.7 (liberté unitaire).** Soient $`K_1,\dots,K_n`$ et $`L_1,\dots,L_n`$ dans $`M_{k'\times k}`$
(on complète la famille la plus courte par des zéros). Alors $`\sum_lK_lAK_l^*=\sum_jL_jAL_j^*`$ pour
tout $`A`$ si et seulement si $`L_j=\sum_lu_{jl}K_l`$ pour une unitaire $`u\in M_n`$.

*Démonstration.* Si $`L_j=\sum_lu_{jl}K_l`$, alors
$`\sum_jL_jAL_j^*=\sum_{l,l'}(u^*u)_{l'l}K_lAK_{l'}^*=\sum_lK_lAK_l^*`$. Réciproquement, deux
applications égales ont des matrices de Choi égales, $`\sum_lv_lv_l^*=\sum_jw_jw_j^*`$ avec
$`v_l=\sum_ie_i\otimes K_le_i`$ et $`w_j=\sum_ie_i\otimes L_je_i`$. Le lemme 6.6, appliqué aux matrices de
colonnes $`v_l`$ et $`w_j`$, donne $`w_j=\sum_lW_{lj}v_l`$ ; comme $`K\mapsto\sum_ie_i\otimes Ke_i`$ est
linéaire et injective, $`L_j=\sum_lW_{lj}K_l`$, et $`u=W^T`$ est unitaire. $`\square`$

*Calcul :* le théorème 6.7 avec $`n=3`$ — [tutoriel, section 6.2](../../tutorial/fr/06_quantum_channels.md#s-6-2).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Tout canal est une interaction isométrique du système avec un
> environnement dans un état fixé, suivie de l'oubli de l'environnement. Mesurer l'environnement dans
> la base $`e_l`$ donne le résultat $`l`$ avec probabilité $`\mathrm{Tr}\,K_l\rho K_l^*`$ et laisse
> l'état $`K_l\rho K_l^*/\mathrm{Tr}(K_l\rho K_l^*)`$ : les opérateurs de Kraus sont les « sauts
> quantiques » possibles. Une autre base donne d'autres sauts, mais le même canal (théorème 6.7), et le
> plus petit environnement est de dimension égale au rang de Kraus, comme on l'a vu après le
> théorème 6.5.

<a id="s-6-3"></a>

## 6.3 Applications duales relativement à une trace

**Définition.** Soient $`\mathrm{tr}_w`$, $`\mathrm{tr}_v`$ des traces fidèles sur $`\mathcal M`$,
$`\mathcal N`$. La **duale** d'une application linéaire $`\Phi:\mathcal M\to\mathcal N`$ relativement à
$`(w,v)`$ est l'application linéaire $`\Phi^\sharp:\mathcal N\to\mathcal M`$ telle que
```math
\mathrm{tr}_v\big(\Phi(x)\,y\big)=\mathrm{tr}_w\big(x\,\Phi^\sharp(y)\big)\qquad(x\in\mathcal M,\ y\in\mathcal N).
```
Elle existe et est unique : pour $`y`$ fixé, $`x\mapsto\mathrm{tr}_v(\Phi(x)y)`$ est une forme
linéaire, donc, d'après le théorème 3.2, elle est égale à
$`\mathrm{tr}_w(\rho x)=\mathrm{tr}_w(x\rho)`$ pour exactement un $`\rho=:\Phi^\sharp(y)`$, et
l'unicité rend $`\Phi^\sharp`$ linéaire. Pour les traces brutes, on écrit $`\Phi^*`$. Posons
$`z_w=\sum_cw_cP_c`$ et $`z_v=\sum_dv_dQ_d`$ (centraux et inversibles).

<a id="r-6-8"></a>

**Proposition 6.8.** Si $`\Phi`$ a pour opérateurs de Kraus les $`K^{dc}_l`$ (proposition 6.4), alors
```math
\Phi^\sharp(y)_c=\sum_d\frac{v_d}{w_c}\sum_lK^{dc*}_l\,y_d\,K^{dc}_l,\qquad\text{i.e.}\qquad
\Phi^\sharp(y)=z_w^{-1}\,\Phi^*(z_v\,y).
```
En particulier, $`\Phi^\sharp`$ est CP, d'opérateurs de Kraus $`(v_d/w_c)^{1/2}K^{dc*}_l`$.

*Démonstration.* Par cyclicité de la trace dans chaque bloc,
```math
\mathrm{tr}_v(\Phi(x)y)=\sum_dv_d\sum_{c,l}\mathrm{Tr}\big(K^{dc}_lx_cK^{dc*}_ly_d\big)
=\sum_cw_c\mathrm{Tr}\Big(x_c\sum_d\frac{v_d}{w_c}\sum_lK^{dc*}_ly_dK^{dc}_l\Big),
```
et l'unicité donne la formule. Pour $`v=w=1`$, elle s'écrit $`\Phi^*(y)_c=\sum_{d,l}K^{dc*}_ly_dK^{dc}_l`$ ;
les facteurs $`v_d`$ et $`1/w_c`$ sont la multiplication par $`z_v`$ et par $`z_w^{-1}`$. $`\square`$

$`\Phi`$ **préserve la trace relativement à $`(w,v)`$** si
$`\mathrm{tr}_v\,\Phi(x)=\mathrm{tr}_w(x)`$ pour tout $`x`$.

<a id="r-6-9"></a>

**Corollaire 6.9.** (a) $`\Phi`$ est TP relativement à $`(w,v)`$ si et seulement si $`\Phi^\sharp(1)=1`$,
c'est-à-dire $`\sum_{d,l}(v_d/w_c)K^{dc*}_lK^{dc}_l=1_{k_c}`$ pour tout $`c`$. (b) $`\Phi`$ est unifère si
et seulement si $`\Phi^\sharp`$ est TP relativement à $`(v,w)`$. (c) La duale de $`\Phi^\sharp`$
relativement à $`(v,w)`$ est $`\Phi`$. (d) Si $`\Phi:\mathcal M\to\mathcal M`$ préserve les secteurs et
$`v=w`$, alors $`\Phi^\sharp=\Phi^*`$ pour tout $`w`$, et $`\Phi`$ est TP pour une trace si et seulement si
elle l'est pour toutes.

*Démonstration.* (a) $`\mathrm{tr}_v\,\Phi(x)=\mathrm{tr}_w(x\Phi^\sharp(1))`$, et
$`\mathrm{tr}_w(x(\Phi^\sharp(1)-1))=0`$ pour tout $`x`$ si et seulement si $`\Phi^\sharp(1)=1`$ (prendre
$`x=(\Phi^\sharp(1)-1)^*`$ ; $`\mathrm{tr}_w`$ est fidèle) ; on utilise ensuite la proposition 6.8.
(b) $`\mathrm{tr}_w\,\Phi^\sharp(y)=\mathrm{tr}_v(\Phi(1)y)`$ est égal à
$`\mathrm{tr}_v(y)`$ pour tout $`y`$ si et seulement si $`\Phi(1)=1`$. (c)
$`\mathrm{tr}_w(\Phi^\sharp(y)x)=\mathrm{tr}_v(y\Phi(x))`$ par la propriété de trace ; on
utilise l'unicité. (d) Seuls les couples $`d=c`$ interviennent, pour lesquels $`v_d/w_c=1`$. $`\square`$

**Densités en représentation de Heisenberg.** Soient $`\Psi:\mathcal N\to\mathcal M`$ linéaire,
$`\omega`$ une forme linéaire sur $`\mathcal M`$ de densité $`\rho`$ pour $`\mathrm{tr}_w`$
(théorème 3.2), et $`\Psi^\sharp:\mathcal M\to\mathcal N`$ la duale de $`\Psi`$ relativement à $`(v,w)`$.
Alors
$`\omega(\Psi(y))=\mathrm{tr}_w(\Psi(y)\rho)=\mathrm{tr}_v(y\Psi^\sharp(\rho))=\mathrm{tr}_v(\Psi^\sharp(\rho)y)`$,
donc $`\omega\circ\Psi`$ a pour densité $`\Psi^\sharp(\rho)`$ pour $`\mathrm{tr}_v`$ ; si $`\Psi`$ est
unifère, $`\Psi^\sharp`$ est TP relativement à $`(w,v)`$ d'après le corollaire 6.9(b) (en échangeant les
rôles de $`\mathcal M`$ et $`\mathcal N`$).

Pour les applications qui changent de secteur, le corollaire 6.9(d) est en défaut. Soient
$`\mathcal M=\mathcal N=M_1\oplus M_2`$ (secteurs $`0`$ et $`1`$) et $`\Phi`$ l'application appelée **perte de
particule** (voir l'encadré ci-dessous), $`\Phi(x)=(x_0+\mathrm{Tr}\,x_1)\oplus0`$,
d'opérateurs de Kraus $`K^{00}=1`$, $`K^{01}_l=\langle l|`$. Elle est TP pour les traces brutes. Pour
$`\mathrm{tr}_w`$ avec $`w=v=(1,\tfrac12)`$, $`\Phi^\sharp(1)=1\oplus2\cdot1_2`$ et
$`\mathrm{tr}_v\,\Phi(x)=x_0+\mathrm{Tr}\,x_1\ne x_0+\tfrac12\mathrm{Tr}\,x_1=\mathrm{tr}_w(x)`$.

*Dans la bibliothèque :* `Channel.adjoint()` est $`\Phi^*`$ et `InterSectorChannel.adjoint(trace=...)` est $`\Phi^\sharp`$ ; vérifications de la proposition 6.8, du corollaire 6.9(a) et de la perte de particule — [tutoriel, section 6.3](../../tutorial/fr/06_quantum_channels.md#s-6-3).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). L'objet indépendant de la trace est l'application de Heisenberg
> $`\Psi:\mathcal N\to\mathcal M`$ sur les observables : un état $`\omega`$ avant le processus devient
> $`\omega\circ\Psi`$ après, et $`\Psi`$ est CP et unifère. Les densités dépendent de la trace
> (théorème 3.2), et les densités se transforment par la duale de $`\Psi`$ pour les traces choisies
> (le paragraphe qui suit le corollaire 6.9). Si les algèbres agissent avec des multiplicités, la
> trace physique a pour poids $`w_c=m_c`$ (section 2.3). Dans l'exemple, le secteur $`0`$ se lit
> « aucune particule » et le secteur $`1`$ « une particule sur deux sites » ; l'application retire
> alors la particule. Elle conserve la probabilité pour la trace brute, qui compte les états, mais pas
> pour `Tr_norm`, qui pondère également les secteurs.

<a id="s-6-4"></a>

## 6.4 Décohérence

Les canaux standard de la bibliothèque agissent à l'intérieur de chaque bloc ($`p,\gamma\in[0,1]`$) :

| canal | $`\Phi(\rho)_c`$ | opérateurs de Kraus dans le bloc $`c`$ |
|---|---|---|
| `dephasing_channel` | $`(1-p)\rho_c+p\mathrm{diag}\,\rho_c`$ | $`\sqrt{1-p}\,1`$, $`\sqrt p\,E_{ii}`$ |
| `depolarizing_channel` | $`(1-p)\rho_c+p\mathrm{Tr}(\rho_c)\,1/k_c`$ | $`\sqrt{1-p}\,1`$, $`\sqrt{p/k_c}\,E_{ij}`$ |
| `center_expectation` | $`\mathrm{Tr}(\rho_c)\,1/k_c`$ | $`\sqrt{1/k_c}\,E_{ij}`$ |
| `amplitude_damping_channel` | les niveaux $`i>0`$ se désexcitent vers $`0`$ | $`E_{00}+\sqrt{1-\gamma}\sum_{i>0}E_{ii}`$, $`\sqrt\gamma\,E_{0i}`$ |
| `random_mixed_unitary_channel` | $`\sum_jp_jU_j\rho U_j^*`$ | $`\sqrt{p_j}\,U_j`$ |

Les formes de Kraus se vérifient par un calcul direct, par exemple
$`\sum_{ij}\frac1kE_{ij}\rho E_{ji}=\frac1k\sum_{ij}\rho_{jj}E_{ii}`$. Les cinq canaux sont TP, et tous
sauf l'amortissement d'amplitude avec $`\gamma>0`$ (si un $`k_c\ge2`$) sont unifères (proposition 6.4 : calculer
$`\sum_lK_l^*K_l`$ et $`\sum_lK_lK_l^*`$ ; pour l'amortissement d'amplitude dans un bloc de taille $`k`$, ils
valent $`1`$ et $`1+(k-1)\gamma E_{00}-\gamma\sum_{i>0}E_{ii}`$, et $`\sum_jp_j=1`$ pour le canal unitaire
mixte). Explicitement, l'amortissement d'amplitude envoie
$`\rho_{00}\mapsto\rho_{00}+\gamma\sum_{i>0}\rho_{ii}`$, $`\rho_{ii}\mapsto(1-\gamma)\rho_{ii}`$ et
$`\rho_{ij}\mapsto(1-\gamma)\rho_{ij}`$ pour $`i\ne j`$ tous deux positifs, et
$`\rho_{0j}\mapsto\sqrt{1-\gamma}\,\rho_{0j}`$, $`\rho_{j0}\mapsto\sqrt{1-\gamma}\,\rho_{j0}`$ pour $`j>0`$.
L'espérance conditionnelle $`E`$ sur le centre (proposition 2.8) est la dépolarisation avec $`p=1`$. Contrairement à la mesure du secteur,
$`X\mapsto\sum_cP_cXP_c`$, qui est l'identité sur $`\mathcal M`$ (section 3.3), elle ne conserve que les
$`p_c`$.

*Dans la bibliothèque :* les cinq canaux du tableau sont des fonctions de `channels` — [tutoriel, section 6.4](../../tutorial/fr/06_quantum_channels.md#s-6-4).

<a id="r-6-10"></a>

**Théorème 6.10 (points fixes des canaux unifères).** Soit $`\Phi(X)=\sum_lK_lXK_l^*`$, avec
$`K_l\in\mathcal M`$, TP et unifère. Alors $`\Phi(X)=X`$ si et seulement si $`X`$ commute avec tous les
$`K_l`$ et $`K_l^*`$. En particulier, les points fixes forment une sous-algèbre de von Neumann de
$`\mathcal M`$.

*Démonstration.* Si $`X`$ commute avec tous les $`K_l`$, alors $`\Phi(X)=X\sum_lK_lK_l^*=X`$.
Réciproquement, $`\Phi(X^*)=\Phi(X)^*`$, donc les parties autoadjointes $`(X+X^*)/2`$ et $`(X-X^*)/2i`$
d'un point fixe sont des points fixes ; supposons $`X=X^*`$. Avec $`Y_l=K_lX-XK_l`$,
```math
\sum_lY_lY_l^*=\sum_l\big(K_lX^2K_l^*-K_lXK_l^*X-XK_lXK_l^*+XK_lK_l^*X\big)
=\Phi(X^2)-\Phi(X)X-X\Phi(X)+X\Phi(1)X=\Phi(X^2)-X^2,
```
en utilisant $`Y_l^*=XK_l^*-K_l^*X`$, $`\Phi(X)=X`$ et $`\Phi(1)=1`$. Sa trace est nulle car $`\Phi`$ est
TP. Comme $`\mathrm{Tr}\,Y_lY_l^*=\|Y_l\|_2^2\ge0`$, chaque $`Y_l=0`$ ; en passant aux adjoints, $`X`$ commute aussi avec $`K_l^*`$. Ces $`X`$ forment le commutant d'un
ensemble autoadjoint, intersecté avec $`\mathcal M`$ : une $`*`$-algèbre contenant $`1`$, donc une algèbre
de von Neumann (théorème 2.2). $`\square`$

Ainsi, le déphasage ($`p>0`$) conserve exactement la diagonale de chaque bloc, et la dépolarisation
($`p>0`$) ainsi que $`E`$ conservent exactement le centre. Sans l'hypothèse unifère, le théorème est en
défaut : l'amortissement d'amplitude fixe $`E_{00}`$, qui ne commute pas avec le saut $`E_{01}`$. Pour
$`0<\gamma\le1`$ et $`k\ge2`$, ses points fixes dans un bloc et le commutant de ses opérateurs de Kraus
sont tous deux de dimension $`1`$, mais ce sont $`\mathbb CE_{00}`$ et $`\mathbb C1`$ : d'après la forme
explicite ci-dessus, un point fixe a $`\rho_{ij}=0`$ pour $`i\ne j`$ et $`\rho_{ii}=0`$ pour $`i>0`$ (ces
coefficients sont multipliés par $`\sqrt{1-\gamma}\ne1`$ ou $`1-\gamma\ne1`$), et tout multiple de $`E_{00}`$
est fixe ; et ce qui commute avec tous les $`E_{0i}`$, $`E_{i0}`$ ($`i>0`$) commute avec leurs produits
$`E_{ij}=E_{i0}E_{0j}`$ et $`E_{00}=E_{01}E_{10}`$, donc avec $`M_k`$, et c'est donc un multiple de $`1`$
(chapitre 2, exercice 1).

*Calcul :* points fixes et commutants des cinq canaux — [tutoriel, section 6.4](../../tutorial/fr/06_quantum_channels.md#s-6-4).

<a id="r-6-11"></a>

**Théorème 6.11 (les canaux unifères augmentent l'entropie).** Soit $`\Phi:\mathcal M\to\mathcal N`$
positive, TP et unifère, et soient $`\rho`$ et $`\Phi(\rho)`$ des états de valeurs propres $`\lambda_i`$ et
$`\mu_j`$. Alors $`\mu=\Delta\lambda`$ pour une matrice bistochastique $`\Delta`$, de sorte que
$`\sum_jf(\mu_j)\ge\sum_if(\lambda_i)`$ pour toute fonction $`f`$ concave sur $`[0,1]`$. En particulier,
$`S(\Phi(\rho))\ge S(\rho)`$ et $`\mathrm{Tr}\,\Phi(\rho)^2\le\mathrm{Tr}\,\rho^2`$.

*Démonstration.* Diagonalisons chaque bloc : $`\rho=\sum_i\lambda_i|e_i\rangle\langle e_i|`$ et
$`\Phi(\rho)=\sum_j\mu_j|f_j\rangle\langle f_j|`$, avec des bases orthonormées de
$`\bigoplus_c\mathbb C^{k_c}`$ et de $`\bigoplus_d\mathbb C^{k'_d}`$ dont chaque vecteur appartient à un
seul facteur, de sorte que $`|e_i\rangle\langle e_i|\in\mathcal M`$. Alors
$`\mu_j=\langle f_j,\Phi(\rho)f_j\rangle=\sum_i\Delta_{ji}\lambda_i`$ avec
$`\Delta_{ji}=\langle f_j,\Phi(|e_i\rangle\langle e_i|)f_j\rangle\ge0`$,
$`\sum_j\Delta_{ji}=\mathrm{Tr}\,\Phi(|e_i\rangle\langle e_i|)=1`$ et
$`\sum_i\Delta_{ji}=\langle f_j,\Phi(1)f_j\rangle=1`$ (donc $`\Delta`$ est carrée : sommer tous les
coefficients). Par l'inégalité de Jensen finie $`f(\sum_i\Delta_{ji}\lambda_i)\ge\sum_i\Delta_{ji}f(\lambda_i)`$
(par récurrence sur le nombre de termes à partir de la définition de la concavité, comme dans la
solution de l'exercice 2 du chapitre 3),
$`\sum_jf(\mu_j)\ge\sum_j\sum_i\Delta_{ji}f(\lambda_i)=\sum_if(\lambda_i)`$. On prend $`f(t)=-t\log t`$
($`f(0)=0`$) et $`f(t)=-t^2`$, qui sont concaves sur $`[0,1]`$ d'après le lemme 3.4 : $`t\log t`$ et $`t^2`$
sont continues sur $`[0,1]`$, de dérivées secondes $`1/t>0`$ et $`2>0`$ sur $`(0,1)`$. $`\square`$

*Calcul :* l'entropie et la pureté sous les canaux unifères, et la baisse d'entropie sous l'amortissement d'amplitude — [tutoriel, section 6.4](../../tutorial/fr/06_quantum_channels.md#s-6-4).

La partie classique d'une application entre secteurs est sa **matrice de transition**
$`T_{dc}(\rho)=\mathrm{Tr}\,\Phi_{dc}(\rho_c)/\mathrm{Tr}\,\rho_c`$ ($`\rho_c\ne0`$) ; la méthode
`transition_matrix()` de la bibliothèque renvoie $`T=T(1)`$, pour l'état maximalement mélangé de chaque
secteur d'entrée.

<a id="r-6-12"></a>

**Proposition 6.12.** Soient $`\Phi:\mathcal M\to\mathcal N`$ positive et TP, $`\rho`$ un état tel que
$`p_c=\mathrm{Tr}\,\rho_c>0`$, et $`q_d=\mathrm{Tr}\,\Phi(\rho)_d`$. Alors $`T(\rho)`$ est
stochastique en colonnes et $`q=T(\rho)p`$. Avec les espérances conditionnelles $`E_{\mathcal M}`$,
$`E_{\mathcal N}`$ sur les centres (proposition 2.8),
```math
E_{\mathcal N}\Phi E_{\mathcal M}(\rho)=\bigoplus_d\frac{(Tp)_d}{k'_d}\,1_{k'_d},\qquad T=T(1).
```

*Démonstration.* $`T_{dc}\ge0`$ car $`\Phi_{dc}`$ est positive, et
$`\sum_dT_{dc}=\mathrm{Tr}\,\Phi(\iota_c\rho_c)/\mathrm{Tr}\,\rho_c=1`$. Ensuite,
$`q_d=\sum_c\mathrm{Tr}\,\Phi_{dc}(\rho_c)=\sum_cT_{dc}(\rho)p_c`$. Enfin,
$`E_{\mathcal M}\rho=\bigoplus_cp_c1_{k_c}/k_c`$, de sorte que la trace de
$`\Phi(E_{\mathcal M}\rho)_d=\sum_cp_c\Phi_{dc}(1_{k_c}/k_c)`$ vaut $`\sum_cT_{dc}p_c`$, et $`E_{\mathcal N}`$
remplace le bloc $`d`$ par sa trace multipliée par $`1/k'_d`$. $`\square`$

*Dans la bibliothèque :* la matrice de transition est `InterSectorChannel.transition_matrix()` ; une vérification de la proposition 6.12 — [tutoriel, section 6.4](../../tutorial/fr/06_quantum_channels.md#s-6-4).

**Cohérences d'un qubit.** Sur $`M_2`$, soient $`D_p`$ le déphasage et $`A_\gamma`$ l'amortissement
d'amplitude. D'après le tableau et la forme explicite ci-dessus, $`D_p`$ multiplie $`\rho_{01}`$ et
$`\rho_{10}`$ par $`1-p`$ et fixe la diagonale, tandis que $`A_\gamma`$ multiplie $`\rho_{01}`$ et $`\rho_{10}`$
par $`\sqrt{1-\gamma}`$ et envoie la diagonale sur $`(\rho_{00}+\gamma\rho_{11},(1-\gamma)\rho_{11})`$,
indépendamment des cohérences. Par conséquent $`D_pA_\gamma=A_\gamma D_p`$ : les deux composées
multiplient les cohérences par $`(1-p)\sqrt{1-\gamma}`$ et agissent sur la diagonale comme $`A_\gamma`$. En
particulier, pour $`T_1,T_\varphi>0`$, $`t\ge0`$, $`1-p=e^{-t/T_\varphi}`$ et $`1-\gamma=e^{-t/T_1}`$, la composée
multiplie les cohérences par $`e^{-t/T_2}`$ avec
```math
\frac1{T_2}=\frac1{2T_1}+\frac1{T_\varphi},
```
et ces familles se composent additivement en $`t`$ : $`D_{p(t)}D_{p(s)}=D_{p(t+s)}`$ et
$`A_{\gamma(t)}A_{\gamma(s)}=A_{\gamma(t+s)}`$, puisque les facteurs $`e^{-t/T}`$ se multiplient et que la
trace est préservée.

*Calcul :* le facteur de cohérence et $`1/T_2=1/(2T_1)+1/T_\varphi`$ — [tutoriel, section 6.4](../../tutorial/fr/06_quantum_channels.md#s-6-4).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Pour un qubit,
> le déphasage modélise le déphasage pur de temps $`T_\varphi`$ et l'amortissement d'amplitude la
> désexcitation du niveau excité de temps $`T_1`$ ; le calcul ci-dessus donne alors le temps de
> cohérence $`T_2`$ d'un qubit soumis aux deux. Le théorème 6.10 est
> la forme algébrique de l'information sans décohérence : ce qui commute avec tous les opérateurs de
> saut survit (chapitre 10). La proposition 6.12 dit que, vue à travers le centre, une application
> entre secteurs est une matrice de Markov classique sur les charges.

<a id="s-6-5"></a>

## 6.5 Canaux aléatoires

Une **isométrie aléatoire de Haar** $`V:\mathbb C^k\to\mathbb C^N`$ est, comme dans la bibliothèque
(`random_isometry`), le facteur $`Q`$ de $`G=QR`$, où $`G`$ est une matrice de Ginibre $`N\times k`$
(coefficients indépendants de densité $`\pi^{-1}e^{-|z|^2}`$), $`Q^*Q=1`$, et $`R`$ est triangulaire
supérieure à diagonale positive (sur le CPU, la bibliothèque prend une factorisation LAPACK
$`G=Q'R'`$ et multiplie la $`i`$-ième colonne de $`Q'`$ par $`R'_{ii}/|R'_{ii}|`$, ce qui donne cette
factorisation, comme dans le théorème 5.3, dès que tous les $`R'_{ii}\ne0`$, c'est-à-dire dès que $`G`$
est de rang $`k`$, ce qui a lieu presque sûrement d'après la démonstration du lemme 6.13 ; pour de
grands lots sur un GPU, elle utilise Gram–Schmidt, qui la donne directement).
`random_channel(alg, r)` tire de telles
$`V_c:\mathbb C^{k_c}\to\mathbb C^r\otimes\mathbb C^{k_c}`$ indépendamment dans chaque bloc et utilise
les opérateurs de Kraus du théorème 6.5.

*Dans la bibliothèque :* `random_isometry` et `channels.random_channel` — [tutoriel, section 6.5](../../tutorial/fr/06_quantum_channels.md#s-6-5).

<a id="r-6-13"></a>

**Lemme 6.13.** $`V`$ a la loi des $`k`$ premières colonnes d'une unitaire de Haar de $`U(N)`$. En
particulier, pour toute unitaire $`U\in M_N`$, $`UV`$ a la même loi que $`V`$.

*Démonstration.* Soient $`Z`$ une matrice de Ginibre $`N\times N`$ dont les $`k`$ premières colonnes sont
$`G`$, et $`E\in M_{N\times k}`$ les $`k`$ premières colonnes de $`1_N`$. D'après le théorème 5.3, $`Z`$ est
presque sûrement inversible et $`Z=Q_ZR_Z`$ avec $`Q_Z`$ distribuée selon la loi de Haar. Comme $`R_Z`$ est
triangulaire supérieure, $`G=ZE=(Q_ZE)(E^*R_ZE)`$ est une factorisation du type voulu, et c'est la
seule : si $`G=QR=Q'R'`$, alors $`R^*R=G^*G=R'^*R'`$, donc $`R'R^{-1}`$ est unitaire, triangulaire
supérieure à diagonale positive, donc égale à $`1`$ (démonstration du lemme 5.2). Ainsi $`V=Q_ZE`$, et
$`UV=(UQ_Z)E`$ a la loi de $`V`$ car $`UQ_Z`$ a celle de $`Q_Z`$. $`\square`$

<a id="r-6-14"></a>

**Lemme 6.14.** Soit $`u`$ un vecteur unitaire aléatoire de $`\mathbb C^N`$ dont la loi est invariante par
toutes les unitaires. Alors $`\mathbb E[u_a\bar u_b]=\delta_{ab}/N`$ et
```math
\mathbb E\big[u_a\bar u_b\,u_c\bar u_d\big]=\frac{\delta_{ab}\delta_{cd}+\delta_{ad}\delta_{cb}}{N(N+1)} .
```

*Démonstration.* L'invariance par $`\mathrm{diag}(e^{i\theta_1},\dots,e^{i\theta_N})`$ multiplie
les moments par $`e^{i(\theta_a-\theta_b)}`$ et $`e^{i(\theta_a-\theta_b+\theta_c-\theta_d)}`$, si bien
qu'ils sont nuls sauf si $`a=b`$, resp. $`\{a,c\}=\{b,d\}`$ en tant que multiensembles. Par invariance par
les matrices de permutation, $`\mathbb E|u_a|^2=1/N`$, et il ne reste que $`\alpha=\mathbb E|u_1|^4`$ et
$`\beta=\mathbb E|u_1|^2|u_2|^2`$. D'abord, $`\sum_c|u_c|^2=1`$ donne
$`\alpha+(N-1)\beta=\mathbb E|u_1|^2=1/N`$. Ensuite, $`(u_1+u_2)/\sqrt2`$ est la première coordonnée de
$`Uu`$ pour une unitaire $`U`$, donc $`\mathbb E|u_1+u_2|^4=4\alpha`$. Avec $`s=|u_1|^2+|u_2|^2`$ et
$`z=u_1\bar u_2`$, $`|u_1+u_2|^4=s^2+2s(z+\bar z)+(z+\bar z)^2`$ ; les invariances $`u_1\mapsto-u_1`$ et
$`u_1\mapsto iu_1`$ donnent $`\mathbb E[s(z+\bar z)]=0`$ et $`\mathbb Ez^2=0`$, donc
$`4\alpha=(2\alpha+2\beta)+2\beta`$, c'est-à-dire $`\alpha=2\beta`$. Par conséquent $`\beta=1/(N(N+1))`$ et
$`\alpha=2\beta`$, ce qui est la formule dans tous les cas. (On a utilisé $`N\ge2`$ ; pour $`N=1`$,
$`|u_1|=1`$ et les deux membres valent $`1`$.) $`\square`$

<a id="r-6-15"></a>

**Théorème 6.15 (pureté de sortie d'un canal aléatoire).** Soient $`\psi`$ un état pur dans un secteur
de taille $`k`$ et $`\Phi`$ tiré par `random_channel` avec un rang de Kraus $`r`$. Alors
```math
\mathbb E\mathrm{Tr}\,\Phi(\psi)^2=\frac{k+r}{kr+1},\qquad\mathbb E\,\Phi(\psi)=\frac1k1_k .
```

*Démonstration.* Dans le secteur de $`\psi`$, $`\Phi(\psi)=\mathrm{Tr}_E|u\rangle\langle u|`$ avec
$`u=V\psi`$, un vecteur unitaire dont la loi est invariante par les unitaires d'après le lemme 6.13.
Indexons $`\mathbb C^r\otimes\mathbb C^k`$ par des couples $`(e,x)`$. Alors
$`\Phi(\psi)_{xy}=\sum_eu_{ex}\bar u_{ey}`$, de moyenne $`r\delta_{xy}/(kr)`$, et
$`\mathrm{Tr}\,\Phi(\psi)^2=\sum_{x,y,e,f}u_{ex}\bar u_{ey}u_{fy}\bar u_{fx}`$. D'après le lemme 6.14
avec $`N=kr`$, le premier terme impose $`x=y`$ ($`kr^2`$ uplets d'indices) et le second $`e=f`$ ($`k^2r`$
uplets), de sorte que la moyenne vaut $`(kr^2+k^2r)/(kr(kr+1))`$. $`\square`$

*Calcul :* le théorème 6.15 par Monte-Carlo pour $`r=1,2,4`$ — [tutoriel, section 6.5](../../tutorial/fr/06_quantum_channels.md#s-6-5).

Pour $`r=1`$, le canal est unitaire et la pureté vaut $`1`$ ; quand $`r\to\infty`$, elle tend vers $`1/k`$, la
pureté de l'état maximalement mélangé. La démonstration n'a utilisé que l'invariance unitaire de la
loi de $`u`$ ; elle donne donc aussi la pureté moyenne de Lubkin de l'état réduit d'un état pur
aléatoire $`u=g/\|g\|`$, $`g`$ gaussien (dont la loi est unitairement invariante d'après le
théorème A.20(b)), que le chapitre 7, théorème 7.11, démontre à l'aide d'intégrales gaussiennes ; les
nombres $`\alpha,\beta`$ du lemme 6.14 sont ceux du corollaire 5.11 pour les coefficients d'une unitaire
de Haar, et le corollaire A.21 les recalcule pour $`u=g/\|g\|`$.

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Un canal aléatoire modélise un bruit générique, sans structure, dû à un
> environnement de dimension $`r`$. Dès $`r=4`$, la pureté d'un état pur de qutrit tombe en moyenne à
> $`7/13\approx0.54`$ (théorème 6.15), contre $`1`$ sans bruit et le minimum $`1/3`$ (par Cauchy–Schwarz,
> $`\mathrm{Tr}\,\rho^2\ge(\mathrm{Tr}\,\rho)^2/k`$). Dans un système à supersélection, le
> bruit agit séparément dans chaque secteur, et la dimension pertinente est celle du secteur, non celle
> de l'espace entier.

<a id="s-6-6"></a>

## 6.6 Points de comparaison

*Dans la bibliothèque :* la comparaison avec d'autres logiciels — [tutoriel, section 6.6](../../tutorial/fr/06_quantum_channels.md#s-6-6).

<a id="x-6"></a>

## Exercices

1. (Werner–Holevo) Montrer que $`\Phi(\rho)=\big(\mathrm{Tr}(\rho)1-\rho^T\big)/(k-1)`$ sur $`M_k`$,
   $`k\ge2`$, est un canal unifère de rang de Kraus $`k(k-1)/2`$, bien qu'il soit construit à partir de
   la transposition.
2. Montrer que pour $`0<p\le1`$ le canal de déphasage sur $`M_k`$ a pour rang de Kraus $`k`$, soit un de
   moins que le nombre d'opérateurs de Kraus stockés par la bibliothèque.
3. Soient $`X`$ le décalage cyclique sur $`\mathbb C^k`$ et $`Z=\mathrm{diag}(1,\omega,\dots,\omega^{k-1})`$,
   $`\omega=e^{2\pi i/k}`$. Montrer que $`\frac1{k^2}\sum_{a,b=0}^{k-1}X^aZ^b\rho Z^{-b}X^{-a}=\mathrm{Tr}(\rho)1/k`$ :
   l'espérance conditionnelle sur le centre de $`M_k`$ est un canal unitaire mixte.
4. Montrer que $`(\Psi\circ\Phi)^\sharp=\Phi^\sharp\circ\Psi^\sharp`$ lorsqu'on utilise la même trace sur
   l'algèbre intermédiaire.

*Calcul :* les quatre exercices, numériquement — [tutoriel, Exercices](../../tutorial/fr/06_quantum_channels.md#x-6).

<details>
<summary>Solutions</summary>

1. En notant $`F`$ l'opérateur d'échange, $`J(\theta)=\sum_{ij}E_{ij}\otimes E_{ji}=F`$ et $`J(\rho\mapsto\mathrm{Tr}(\rho)1)=1\otimes1`$,
   donc $`J(\Phi)=(1-F)/(k-1)=2P_-/(k-1)`$, où $`P_-=(1-F)/2`$ est le projecteur sur les vecteurs
   antisymétriques ($`F=F^*`$, $`F^2=1`$), de rang $`\mathrm{Tr}\,P_-=(k^2-\mathrm{Tr}\,F)/2`$ avec
   $`\mathrm{Tr}\,F=\sum_{ij}\langle e_i\otimes e_j,e_j\otimes e_i\rangle=k`$. Ainsi $`J\ge0`$, de rang $`k(k-1)/2`$ (théorème 6.2,
   corollaire 6.3), et $`\mathrm{Tr}\,\Phi(\rho)=\mathrm{Tr}\,\rho`$, $`\Phi(1)=1`$ directement.

2. $`J=(1-p)|\Omega\rangle\langle\Omega|+p\sum_iE_{ii}\otimes E_{ii}`$ est nulle en dehors du sous-espace
   engendré par les $`e_i\otimes e_i`$, et sur cet espace de dimension $`k`$ c'est la matrice
   $`(1-p)\mathbf 1\mathbf 1^T+p1`$ ($`\mathbf 1`$ le vecteur dont toutes les composantes valent $`1`$), de
   valeurs propres $`p`$ et $`(1-p)k+p`$, toutes strictement positives. Ainsi
   $`\mathrm{rank}\,J=k`$, tandis que la bibliothèque stocke $`1+k_{\max}`$ opérateurs.

3. $`\mathrm{Tr}((X^aZ^b)^*X^{a'}Z^{b'})=k\,\delta_{aa'}\delta_{bb'}`$ : $`Z^{-b}X^{a'-a}Z^{b'}`$ a une
   diagonale nulle sauf si $`a=a'`$, et $`\mathrm{Tr}\,Z^{b'-b}=\sum_j\omega^{j(b'-b)}=k\delta_{bb'}`$.
   Les $`B_{ab}=X^aZ^b/\sqrt k`$ forment donc une base orthonormée de $`M_k`$ pour
   $`\mathrm{Tr}(A^*B)`$, et deux bases orthonormées quelconques se déduisent l'une de l'autre par
   une unitaire ; d'après le théorème 6.7, $`\sum_mB_m\rho B_m^*`$ est le même pour toutes les bases
   orthonormées. Pour la base des $`E_{ab}`$, il vaut $`\sum_{ab}E_{ab}\rho E_{ba}=\mathrm{Tr}(\rho)1`$.
   Il reste à diviser par $`k`$.

4. $`\mathrm{tr}(\Psi(\Phi(x))z)=\mathrm{tr}(\Phi(x)\Psi^\sharp(z))=\mathrm{tr}(x\,\Phi^\sharp(\Psi^\sharp(z)))`$
   pour tous $`x`$ et $`z`$ ; l'unicité de la duale donne le résultat.

</details>

**Suite :** le [chapitre 7](07_tensor_products.md) traite des produits tensoriels et de l'intrication : la
trace partielle — un canal entre algèbres différentes dont la duale est le plongement
$`A\mapsto A\otimes1`$ —, la formule de Page et la fusion des charges.
