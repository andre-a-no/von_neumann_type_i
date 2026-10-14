[English](../03_traces_and_states.md) · [Español](../es/03_traces_and_states.md) · **Français** · [Deutsch](../de/03_traces_and_states.md) · [中文](../zh/03_traces_and_states.md) · [日本語](../ja/03_traces_and_states.md) · [Русский](../ru/03_traces_and_states.md)

# 3. Traces, états et entropie

Sur $`M_n(\mathbb C)`$, il n'existe pour l'essentiel qu'une seule trace. Sur
$`\bigoplus_cM_{k_c}`$, on dispose d'une famille à $`C`$ paramètres, et le choix importe :
densités, entropies et applications duales en dépendent toutes. Ce chapitre classifie les
traces, représente les états par des densités, démontre les inégalités d'entropie
fondamentales, et se termine par une inégalité de traces plausible qui est fausse.

**Vous apprendrez** la classification des traces sur $`\bigoplus_cM_{k_c}`$ et la relation entre les densités
d'un état relativement à différentes traces, la partie classique d'un état, la décomposition de l'entropie de von Neumann
selon les secteurs, l'inégalité de Klein et ses conséquences, ainsi que le contraste de
Michelson. Les vérifications numériques de ce chapitre se trouvent dans le
[Tutoriel 3](../../tutorial/fr/03_traces_and_states.md).

Dans tout ce chapitre, $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$ avec tous les
$`k_c\ge1`$, $`D=\sum_ck_c`$, les $`P_c`$ sont les projecteurs centraux minimaux, et
$`A=\bigoplus_cA_c`$.

<a id="s-3-1"></a>

## 3.1 Traces

Une **trace** sur $`\mathcal M`$ est une forme linéaire $`\tau`$ telle que $`\tau(AB)=\tau(BA)`$
pour tous $`A,B`$. Elle est **positive** si $`\tau(A^*A)\ge0`$ pour tout $`A`$, et **fidèle** si,
de plus, $`\tau(A^*A)=0`$ entraîne $`A=0`$.

<a id="r-3-1"></a>

**Théorème 3.1.** Toute trace sur $`\mathcal M`$ est de la forme
$`\mathrm{tr}_w(A)=\sum_cw_c\mathrm{Tr}\,A_c`$ avec $`w\in\mathbb C^C`$. Elle est positive
si et seulement si tous les $`w_c\ge0`$, et fidèle si et seulement si tous les $`w_c>0`$.

*Démonstration.* Soient $`E^c_{ij}`$ les unités matricielles du $`c`$-ième bloc. Pour une trace
$`\tau`$,
$`\tau(E^c_{ij})=\tau(E^c_{i1}E^c_{1j})=\tau(E^c_{1j}E^c_{i1})=\delta_{ij}\tau(E^c_{11})`$.
Posons $`w_c=\tau(E^c_{11})`$ ; par linéarité, $`\tau(A)=\sum_c\sum_i(A_c)_{ii}w_c`$.
Réciproquement, toute $`\mathrm{tr}_w`$ est une trace. Positivité :
$`E^c_{11}=(E^c_{11})^*E^c_{11}`$, donc la condition $`w_c=\tau\big((E^c_{11})^*E^c_{11}\big)\ge0`$ est
nécessaire, et $`w_c>0`$ est nécessaire pour la fidélité, puisque $`E^c_{11}\neq0`$. Réciproquement,
$`\mathrm{tr}_w(A^*A)=\sum_cw_c\|A_c\|_2^2`$, qui est $`\ge0`$ si tous les $`w_c\ge0`$, et ne
s'annule que pour $`A=0`$ si tous les $`w_c>0`$. $`\square`$

La bibliothèque implémente trois traces fidèles :

| nom | poids $`w_c`$ | valeur en $`1`$ | signification |
|---|---|---|---|
| `Tr_blunt` | $`1`$ | $`D`$ | trace de $`\bigoplus_c\mathbb C^{k_c}`$ (multiplicités $`m_c=1`$) |
| `Tr_norm` | $`1/k_c`$ | $`C`$ | chaque bloc normalisé |
| `tau_vN` | $`1/(Ck_c)`$ | $`1`$ | l'état tracial de poids égaux $`1/C`$ sur les secteurs |

D'après la section 2.3, la trace d'un espace sur lequel $`\mathcal M`$ agit avec les
multiplicités $`m_c`$ est $`\mathrm{tr}_w`$ avec $`w_c=m_c`$. Choisir une trace fidèle, c'est
donc choisir « combien compte chaque secteur ».

*Dans la bibliothèque :* les traces `Tr_blunt`, `Tr_norm`, `tau_vN`, leurs valeurs en $`1`$ et la propriété de trace du théorème 3.1 — [tutoriel, section 3.1](../../tutorial/fr/03_traces_and_states.md#s-3-1).

<a id="s-3-2"></a>

## 3.2 États et densités

Un **état** est une forme linéaire $`\omega`$ sur $`\mathcal M`$ positive
($`\omega(A^*A)\ge0`$) et normalisée ($`\omega(1)=1`$).

<a id="r-3-2"></a>

**Théorème 3.2.** Soit $`\mathrm{tr}_w`$ une trace fidèle. Pour toute forme linéaire
$`\omega`$, il existe un unique $`\rho_w\in\mathcal M`$ tel que
$`\omega(A)=\mathrm{tr}_w(\rho_wA)`$ pour tout $`A`$. La forme $`\omega`$ est positive si et
seulement si $`\rho_w\ge0`$, et c'est un état si et seulement si, de plus,
$`\mathrm{tr}_w(\rho_w)=1`$. Pour deux traces fidèles, les densités sont reliées bloc par
bloc par
```math
w_c\,(\rho_w)_c=w'_c\,(\rho_{w'})_c .
```

*Démonstration.* $`\langle X,Y\rangle_w=\mathrm{tr}_w(X^*Y)`$ est un produit scalaire sur
$`\mathcal M`$ (théorème 3.1). D'après le théorème de représentation de Riesz dans l'espace
de dimension finie $`\mathcal M`$, il existe un unique $`R`$ tel que
$`\omega(A)=\langle R,A\rangle_w=\mathrm{tr}_w(R^*A)`$ ; on pose $`\rho_w=R^*`$. Si
$`\rho_w\ge0`$, alors $`\omega(A^*A)=\sum_cw_c\mathrm{Tr}(A_c\,(\rho_w)_c\,A_c^*)\ge0`$,
puisque $`A_c(\rho_w)_cA_c^*\ge0`$ (corollaire 1.7 avec $`T=(\rho_w)_c^{1/2}A_c^*`$), qu'une
matrice positive a une trace $`\ge0`$ (ses valeurs propres sont $`\ge0`$, corollaire 1.4(3)) et que
$`w_c>0`$. Réciproquement, si $`\omega`$ est positive, prenons $`x\in\mathbb C^{k_c}`$ et pour $`A`$
l'opérateur de rang un $`|x\rangle\langle x|`$ dans le bloc $`c`$ ($`x\neq0`$) ; alors $`A=A^*A/\|x\|^2`$ et
$`0\le\omega(A)=w_c\langle x,(\rho_w)_cx\rangle`$, donc $`(\rho_w)_c\ge0`$. La normalisation
s'écrit $`\omega(1)=\mathrm{tr}_w(\rho_w)`$. Enfin,
$`\sum_cw_c\mathrm{Tr}((\rho_w)_cA_c)`$ et $`\sum_cw'_c\mathrm{Tr}((\rho_{w'})_cA_c)`$
sont tous deux égaux à $`\omega(A)`$ pour tout $`A`$, donc $`\bigoplus_cw_c(\rho_w)_c`$ et
$`\bigoplus_cw'_c(\rho_{w'})_c`$ sont toutes deux des densités de $`\omega`$ relativement à la trace brute
($`w_c=1`$), et par unicité $`w_c(\rho_w)_c=w'_c(\rho_{w'})_c`$. $`\square`$

*Dans la bibliothèque :* un état est une `DensityMatrix`, stockée par sa densité brute ; les conversions du théorème 3.2 sont `DensityMatrix.density` et `DensityMatrix.from_density` — [tutoriel, section 3.2](../../tutorial/fr/03_traces_and_states.md#s-3-2).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). $`\omega(A)`$ est la valeur moyenne de l'observable $`A`$. La
> densité $`\rho`$ relative à la trace brute est la matrice densité usuelle de la mécanique
> quantique sur $`\bigoplus_c\mathbb C^{k_c}`$. Les densités relatives à d'autres traces
> apparaissent lorsque la même algèbre agit avec des multiplicités (chapitre 2), par
> exemple sur la partie symétrique et la partie non symétrique d'un système de spins.

<a id="s-3-3"></a>

## 3.3 La partie classique d'un état

La restriction de $`\omega`$ au centre est une loi de probabilité sur les secteurs,
$`p_c=\omega(P_c)=\mathrm{Tr}\,\rho_c`$. Si $`p_c>0`$, l'**état conditionnel** dans le
secteur $`c`$ a pour densité $`\sigma_c=\rho_c/p_c`$ sur $`\mathbb C^{k_c}`$, et
```math
\rho=\bigoplus_cp_c\,\sigma_c
```
(si $`p_c=0`$, alors $`\rho_c=0`$, puisque les valeurs propres de $`\rho_c\ge0`$ sont $`\ge0`$ et de
somme $`\mathrm{Tr}\,\rho_c=0`$, corollaire 1.4(3), donc toutes nulles, et que $`\rho_c`$ est
unitairement diagonalisable, théorème 1.3 ; le terme est omis).
Ainsi, tout état sur $`\mathcal M`$ est une combinaison convexe d'états des secteurs.
L'application $`\rho\mapsto\sum_cP_c\rho P_c`$ est l'identité sur $`\mathcal M`$, car tout élément de
$`\mathcal M`$ est diagonal par blocs.

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Un état
> sur $`\mathcal M`$ est un mélange classique d'états quantiques dans les secteurs. L'application
> $`\rho\mapsto\sum_cP_c\rho P_c`$ est la mesure de « quel secteur », et elle ne perturbe aucun
> état de $`\mathcal M`$. C'est là le contenu opérationnel d'une règle de supersélection.

*Calcul :* les probabilités des secteurs $`p_c`$ (`DensityMatrix.sector_probabilities`) et la mesure $`\rho\mapsto\sum_cP_c\rho P_c`$ (`channels.measurement_channel`), qui ne perturbe pas $`\rho`$ — [tutoriel, section 3.3](../../tutorial/fr/03_traces_and_states.md#s-3-3).

<a id="s-3-4"></a>

## 3.4 Entropie

L'**entropie de von Neumann** d'un état de densité (brute) $`\rho`$ est
$`S(\rho)=-\mathrm{Tr}\,\rho\log\rho=-\sum_i\lambda_i\log\lambda_i`$, la somme portant sur
les valeurs propres de $`\rho`$, avec la convention $`0\log0=0`$ ; autrement dit,
$`\rho\log\rho:=f(\rho)`$ pour $`f(t)=t\log t`$, $`f(0)=0`$, et
$`\mathrm{Tr}\,f(\rho)=\sum_if(\lambda_i)`$ par définition de $`f(\rho)=U\,\mathrm{diag}(f(\lambda_i))\,U^*`$
(section 1.3).

<a id="r-3-3"></a>

**Théorème 3.3 (entropie d'un mélange de secteurs).** Avec $`p_c`$ et $`\sigma_c`$ comme
ci-dessus (les secteurs tels que $`p_c=0`$ étant omis),
```math
S(\rho)=H(p)+\sum_cp_c\,S(\sigma_c),\qquad H(p)=-\sum_cp_c\log p_c .
```

*Démonstration.* Les valeurs propres non nulles de $`\rho`$ (avec multiplicité) sont les nombres $`p_c\mu`$, où
$`\mu`$ parcourt les valeurs propres non nulles des $`\sigma_c`$ (secteurs tels que $`p_c>0`$) ;
les valeurs propres nulles contribuent pour $`0`$ aux deux membres. Par conséquent
$`-\sum_c\sum_{\mu>0}p_c\mu\log(p_c\mu)=-\sum_cp_c\log p_c\sum_{\mu>0}\mu-\sum_cp_c\sum_{\mu>0}\mu\log\mu`$,
et $`\sum_{\mu>0}\mu=\mathrm{Tr}\,\sigma_c=1`$. $`\square`$

*Calcul :* le théorème 3.3 pour des états aléatoires, avec `DensityMatrix.von_neumann_entropy` — [tutoriel, section 3.4](../../tutorial/fr/03_traces_and_states.md#s-3-4).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). La formule sépare l'incertitude sur le système en deux parts :
> $`H(p)`$, l'incertitude classique sur la charge (le secteur où se trouve le système), et
> $`\sum_cp_cS(\sigma_c)`$, l'incertitude quantique moyenne à l'intérieur des secteurs. Sous une règle
> de supersélection, la première part est l'entropie d'une variable aléatoire classique : la
> charge.

L'entropie relativement à une autre trace est $`S_w(\rho)=-\mathrm{tr}_w(\rho_w\log\rho_w)`$,
avec la même convention $`\rho_w\log\rho_w:=f(\rho_w)`$. D'après le théorème 3.2, $`(\rho_w)_c=\rho_c/w_c`$, donc
```math
S_w(\rho)=-\sum_cw_c\mathrm{Tr}\Big(\frac{\rho_c}{w_c}\log\frac{\rho_c}{w_c}\Big)
=S(\rho)+\sum_cp_c\log w_c .
```
L'entropie est une propriété de l'état *et* de la trace : un même état a des entropies
différentes relativement à `Tr_blunt` et à `tau_vN`.

*Calcul :* l'entropie relativement à $`\tau_{\mathrm{vN}}`$ et la formule de $`S_w(\rho)`$ — [tutoriel, section 3.4](../../tutorial/fr/03_traces_and_states.md#s-3-4).

<a id="s-3-5"></a>

## 3.5 Entropie relative et inégalité de Klein

Pour des états de densités $`\rho,\sigma`$, avec $`\sigma>0`$, l'**entropie relative** est
$`S(\rho\|\sigma)=\mathrm{Tr}\,\rho\log\rho-\mathrm{Tr}\,\rho\log\sigma`$, où
$`\mathrm{Tr}\,\rho\log\rho:=\sum_i\lambda_i\log\lambda_i`$, la somme portant sur les
valeurs propres de $`\rho`$, avec $`0\log0=0`$ (et de même dans chaque bloc). Elle ne dépend
pas de la trace : d'après le théorème 3.2, $`(\rho_w)_c=\rho_c/w_c`$ et $`(\sigma_w)_c=\sigma_c/w_c`$, donc
les valeurs propres de $`(\rho_w)_c`$ sont celles de $`\rho_c`$ divisées par $`w_c`$,
$`\log(\sigma_c/w_c)=\log\sigma_c-(\log w_c)1`$ (si $`\sigma_c=U\,\mathrm{diag}(s_i)\,U^*`$, alors
$`\sigma_c/w_c=U\,\mathrm{diag}(s_i/w_c)\,U^*`$, et théorème 1.5(a)), et
$`\mathrm{tr}_w(\rho_w\log\rho_w)`$ comme $`\mathrm{tr}_w(\rho_w\log\sigma_w)`$ diffèrent de
$`\mathrm{Tr}\,\rho\log\rho`$ et de $`\mathrm{Tr}\,\rho\log\sigma`$ de la même quantité
$`-\sum_cp_c\log w_c`$ (le logarithme de $`\rho`$ n'intervient que sur son support, par la convention
$`0\log0=0`$). La différence est donc inchangée.

<a id="r-3-4"></a>

**Lemme 3.4 (inégalité de la tangente).** Soit $`f:I\to\mathbb R`$ convexe sur un intervalle $`I`$
et dérivable en $`b\in I`$ (à droite ou à gauche si $`b`$ est une extrémité de $`I`$). Alors
$`f(a)\ge f(b)+(a-b)f'(b)`$ pour tout $`a\in I`$, avec inégalité stricte pour $`a\neq b`$ si $`f`$ est
strictement convexe, c'est-à-dire si $`f((1-t)x+ty)<(1-t)f(x)+tf(y)`$ pour $`x\neq y`$ dans $`I`$ et
$`t\in(0,1)`$. De plus, si $`f`$ est continue sur $`I`$ et deux fois dérivable avec $`f''>0`$ sur
l'intérieur de $`I`$, alors $`f`$ est strictement convexe.

*Démonstration.* Pour $`t\in(0,1]`$, la convexité donne $`f(b+t(a-b))\le(1-t)f(b)+tf(a)`$,
c'est-à-dire $`\big(f(b+t(a-b))-f(b)\big)/t\le f(a)-f(b)`$ ; en faisant tendre $`t`$ vers $`0`$,
on obtient l'inégalité. Si $`f`$ est strictement convexe et $`a\neq b`$, posons $`m=(a+b)/2`$ ;
alors $`f(a)-f(b)>2\big(f(m)-f(b)\big)\ge2(m-b)f'(b)=(a-b)f'(b)`$, la première étape par
stricte convexité en $`m`$, la seconde par l'inégalité qui vient d'être démontrée, appliquée
à $`m`$. Pour la dernière assertion, soient $`x<y`$ dans $`I`$, $`t\in(0,1)`$ et $`z=(1-t)x+ty`$. D'après
le théorème des accroissements finis, il existe $`\xi_1\in(x,z)`$ et $`\xi_2\in(z,y)`$ tels que
$`f(z)-f(x)=t(y-x)f'(\xi_1)`$ et $`f(y)-f(z)=(1-t)(y-x)f'(\xi_2)`$, et $`f'(\xi_1)<f'(\xi_2)`$ d'après
le théorème des accroissements finis appliqué à $`f'`$ sur $`[\xi_1,\xi_2]`$, puisque $`f''>0`$ sur cet
intervalle. En multipliant par $`t(1-t)(y-x)>0`$, on obtient
$`(1-t)\big(f(z)-f(x)\big)<t\big(f(y)-f(z)\big)`$, c'est-à-dire $`f(z)<(1-t)f(x)+tf(y)`$. $`\square`$

<a id="r-3-5"></a>

**Théorème 3.5 (inégalité de Klein).** Soit $`f:I\to\mathbb R`$ convexe sur un intervalle
$`I`$, et $`X,Y`$ autoadjoints de spectres contenus dans $`I`$, $`f`$ étant dérivable en chaque
valeur propre de $`Y`$ (à droite ou à gauche en une extrémité de $`I`$). Alors $`\mathrm{Tr}[f(X)-f(Y)-(X-Y)f'(Y)]\ge0`$. Si $`f`$ est
strictement convexe, l'égalité n'a lieu que pour $`X=Y`$.

*Démonstration.* Écrivons $`X=\sum_ia_i|u_i\rangle\langle u_i|`$ et
$`Y=\sum_jb_j|v_j\rangle\langle v_j|`$ dans des bases propres orthonormées, et posons
$`c_{ij}=|\langle u_i,v_j\rangle|^2`$, de sorte que $`\sum_ic_{ij}=\sum_jc_{ij}=1`$. Alors
$`\mathrm{Tr}\,f(X)=\sum_{ij}c_{ij}f(a_i)`$, $`\mathrm{Tr}\,f(Y)=\sum_{ij}c_{ij}f(b_j)`$ et
$`\mathrm{Tr}\,Xf'(Y)=\sum_{ij}c_{ij}a_if'(b_j)`$, $`\mathrm{Tr}\,Yf'(Y)=\sum_{ij}c_{ij}b_jf'(b_j)`$.
Le membre de gauche vaut donc $`\sum_{ij}c_{ij}[f(a_i)-f(b_j)-(a_i-b_j)f'(b_j)]`$, somme de
termes positifs d'après le lemme 3.4. Si $`f`$ est strictement convexe et que
la somme est nulle, alors $`a_i=b_j`$ dès que $`c_{ij}>0`$, donc
$`Xv_j=\sum_ia_i\langle u_i,v_j\rangle u_i=b_jv_j`$ pour tout $`j`$, et $`X=Y`$. $`\square`$

<a id="r-3-6"></a>

**Corollaire 3.6.** $`S(\rho\|\sigma)\ge0`$, avec égalité si et seulement si $`\rho=\sigma`$.
En particulier $`S(\rho)\le\log D`$, avec égalité si et seulement si $`\rho=1/D`$.

*Démonstration.* On applique le théorème 3.5 à $`f(t)=t\log t`$, $`f(0)=0`$, qui est
strictement convexe sur $`I=[0,\infty)`$ d'après le lemme 3.4 ($`f`$ est continue en $`0`$, car
$`t\log t\to0`$ quand $`t\to0^+`$, et $`f''(t)=1/t>0`$ sur $`(0,\infty)`$) et dérivable sur $`(0,\infty)\supseteq\mathrm{spec}\,\sigma`$ avec
$`f'(t)=\log t+1`$, et à $`X=\rho`$, $`Y=\sigma`$ : le membre de gauche vaut
$`\mathrm{Tr}\,\rho\log\rho-\mathrm{Tr}\,\sigma\log\sigma-\mathrm{Tr}(\rho-\sigma)(\log\sigma+1) =S(\rho\|\sigma)`$, puisque $`\mathrm{Tr}\,\rho=\mathrm{Tr}\,\sigma=1`$. Pour la seconde
assertion, on prend $`\sigma=1/D`$ : $`0\le S(\rho\|1/D)=-S(\rho)+\log D`$. $`\square`$

*Calcul :* $`S(\rho)\le\log D`$ et l'état maximalement mélangé — [tutoriel, section 3.5](../../tutorial/fr/03_traces_and_states.md#s-3-5).

<a id="r-3-7"></a>

**Théorème 3.7 (règle de la chaîne ; monotonie sous l'espérance conditionnelle sur le centre).** Soient
$`\rho,\sigma`$ des états avec $`\sigma>0`$, $`p_c=\mathrm{Tr}\,\rho_c`$,
$`q_c=\mathrm{Tr}\,\sigma_c>0`$, $`\hat\sigma_c=\sigma_c/q_c`$ et, pour $`p_c>0`$,
$`\hat\rho_c=\rho_c/p_c`$. Alors, les sommes portant sur les $`c`$ tels que $`p_c>0`$ (les autres
secteurs contribuent pour $`0`$ en vertu de $`0\log0=0`$),
```math
S(\rho\|\sigma)=\sum_cp_c\log\frac{p_c}{q_c}+\sum_cp_c\,S(\hat\rho_c\|\hat\sigma_c)
\ \ge\ \sum_cp_c\log\frac{p_c}{q_c}=S(E\rho\|E\sigma),
```
où $`E`$ est l'espérance conditionnelle sur le centre (proposition 2.8).

*Démonstration.* Pour $`p_c>0`$,
$`\mathrm{Tr}\,\rho_c\log\rho_c=p_c\log p_c+p_c\mathrm{Tr}\,\hat\rho_c\log\hat\rho_c`$
(valeurs propres $`p_c\mu`$, comme dans le théorème 3.3) et
$`\log\sigma_c=(\log q_c)1+\log\hat\sigma_c`$ (on diagonalise $`\hat\sigma_c=U\,\mathrm{diag}(s_i)\,U^*`$, de
sorte que $`\sigma_c=U\,\mathrm{diag}(q_cs_i)\,U^*`$, et on utilise le théorème 1.5(a)), donc
$`\mathrm{Tr}\,\rho_c\log\sigma_c=p_c\log q_c+p_c\mathrm{Tr}\,\hat\rho_c\log\hat\sigma_c`$ et
$`\mathrm{Tr}\,\rho_c\log\rho_c-\mathrm{Tr}\,\rho_c\log\sigma_c=p_c\log(p_c/q_c)+p_cS(\hat\rho_c\|\hat\sigma_c)`$.
Si $`p_c=0`$, alors $`\rho_c=0`$ et le secteur contribue pour $`0`$.
On somme sur $`c`$ ; la seconde somme est $`\ge0`$ d'après le corollaire 3.6, appliqué dans
$`M_{k_c}`$ aux états $`\hat\rho_c`$ et $`\hat\sigma_c>0`$. Enfin
$`E\rho=\bigoplus_c(p_c/k_c)1_{k_c}`$ et $`E\sigma=\bigoplus_c(q_c/k_c)1_{k_c}>0`$ sont des états ; $`E\rho`$ est
la densité de l'état $`\omega\circ E`$ (où $`\omega(A)=\mathrm{Tr}\,\rho A`$), puisque
$`\mathrm{Tr}(\rho E(A))=\sum_cp_c\mathrm{Tr}(A_c)/k_c=\mathrm{Tr}(E(\rho)A)`$. Et
$`S(E\rho\|E\sigma)=\sum_ck_c\frac{p_c}{k_c}\big(\log\frac{p_c}{k_c}-\log\frac{q_c}{k_c}\big)=\sum_cp_c\log\frac{p_c}{q_c}`$,
car les valeurs propres de $`E\rho`$ dans le bloc $`c`$ sont $`p_c/k_c`$ ($`k_c`$ fois) et les secteurs
tels que $`p_c=0`$ contribuent pour $`0`$. $`\square`$

*Calcul :* le théorème 3.7 pour des états aléatoires, avec `states.relative_entropy` et `channels.center_expectation` — [tutoriel, section 3.5](../../tutorial/fr/03_traces_and_states.md#s-3-5).

Ainsi, remplacer les deux états par leurs images par $`E`$ ne peut que diminuer l'entropie
relative ; ce qui subsiste est la divergence de Kullback–Leibler $`\sum_cp_c\log(p_c/q_c)`$ entre
les lois des secteurs.

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours).
> L'entropie relative mesure à quel point deux états peuvent être distingués. Oublier
> l'information quantique contenue dans les secteurs ne peut que rendre deux états plus
> difficiles à distinguer ; ce qui subsiste est la discernabilité classique des lois de la charge.

<a id="s-3-6"></a>

## 3.6 Le contraste de Michelson

Pour $`A\ge0`$, $`A\neq0`$, le **contraste de Michelson** est
$`\Delta(A)=(\lambda_{\max}-\lambda_{\min})/(\lambda_{\max}+\lambda_{\min})\in[0,1]`$, les valeurs
propres extrêmes étant prises sur l'ensemble des secteurs (et $`\Delta(0)=0`$). Il mesure
l'écart entre $`A`$ et un multiple de l'identité : $`\Delta(A)=0`$ si et seulement si
$`A=\lambda1`$.

<a id="r-3-8"></a>

**Proposition 3.8.** Si $`\Phi:\mathcal M\to\mathcal M`$ est linéaire, positive et unifère
($`\Phi(1)=1`$), alors $`\Delta(\Phi(A))\le\Delta(A)`$ pour tout $`A\ge0`$.

*Démonstration.* $`\lambda_{\min}1\le A\le\lambda_{\max}1`$ (corollaire 1.4(3) pour
$`A-\lambda_{\min}1`$ et $`\lambda_{\max}1-A`$), et la positivité, la linéarité ainsi que
$`\Phi(1)=1`$ donnent $`\lambda_{\min}1\le\Phi(A)\le\lambda_{\max}1`$. Par conséquent $`\Phi(A)`$ est
autoadjoint (une matrice positive l'est, corollaire 1.4(3)), $`\Phi(A)\ge0`$ (car
$`\lambda_{\min}\ge0`$), et d'après le corollaire 1.4(3) $`\lambda_{\min}(\Phi(A))\ge\lambda_{\min}`$ et
$`\lambda_{\max}(\Phi(A))\le\lambda_{\max}`$. Si $`\Phi(A)=0`$, l'assertion est triviale. Sinon $`A\neq0`$,
$`\lambda_{\max}`$ et $`\lambda_{\max}(\Phi(A))`$ sont tous deux $`>0`$, et en posant
$`r=\lambda_{\min}/\lambda_{\max}\in[0,1]`$,
$`\Delta=(1-r)/(1+r)`$ est une fonction décroissante de $`r`$, et $`r(\Phi(A))\ge r(A)`$. $`\square`$

*Calcul :* la proposition 3.8 pour un canal aléatoire à mélange d'unitaires et pour $`E`$, avec `Operator.michelson_contrast` — [tutoriel, section 3.6](../../tutorial/fr/03_traces_and_states.md#s-3-6).

<a id="s-3-7"></a>

## 3.7 Une inégalité plausible mais fausse

Pour $`X,Y`$ positifs dans $`M_k`$ et $`U`$ unitaire, a-t-on
$`|\mathrm{Tr}(XUY)|\le\mathrm{Tr}(XY)`$ ? C'est vrai pour $`U=1`$, et cela ressemble à
un énoncé de type Cauchy–Schwarz. C'est pourtant faux dès que $`X`$ et $`Y`$ ne commutent pas.

<a id="r-3-9"></a>

**Proposition 3.9.** Soient $`X,Y\ge0`$ dans $`M_k`$. Alors $`\max_U|\mathrm{Tr}(XUY)|=\|YX\|_1\ge\mathrm{Tr}(XY)`$, le maximum
portant sur les unitaires $`U\in M_k`$, avec égalité si et seulement si $`XY=YX`$.

*Démonstration.* $`\mathrm{Tr}(XUY)=\mathrm{Tr}(YX\,U)`$, donc le maximum vaut
$`\|YX\|_1`$ d'après le théorème 1.9, et $`\|YX\|_1\ge|\mathrm{Tr}(YX)|=\mathrm{Tr}(XY)`$,
puisque $`\mathrm{Tr}(XY)=\mathrm{Tr}(X^{1/2}YX^{1/2})\ge0`$, car
$`X^{1/2}YX^{1/2}=(Y^{1/2}X^{1/2})^*(Y^{1/2}X^{1/2})\ge0`$ (corollaire 1.7). Si $`XY=YX`$, alors
$`X^{1/2}`$ commute avec $`Y`$ (théorème 1.5(d)), donc $`YX=X^{1/2}YX^{1/2}\ge0`$ et
$`\|YX\|_1=\mathrm{Tr}(YX)`$ (pour $`P\ge0`$, $`\|P\|_1=\mathrm{Tr}\,P`$, section 1.4).
Réciproquement, supposons $`\|YX\|_1=\mathrm{Tr}(YX)`$. Si $`|YX|=0`$, alors
$`YX=0=(YX)^*=XY`$. Sinon, écrivons $`YX=W|YX|`$ (théorème 1.8) et $`R=|YX|^{1/2}\neq0`$. Par
cyclicité, $`\mathrm{Tr}(YX)=\mathrm{Tr}(RWR)=\langle R,WR\rangle_{HS}`$, tandis que
$`\|R\|_2^2=\|WR\|_2^2=\mathrm{Tr}|YX|=\|YX\|_1`$. Donc $`\langle R,WR\rangle_{HS}=\|R\|_2\|WR\|_2`$ :
il y a égalité dans l'inégalité de Cauchy–Schwarz, donc $`WR=\lambda R`$ pour un scalaire
$`\lambda`$, et $`\lambda\|R\|_2^2=\langle R,WR\rangle_{HS}=\|R\|_2^2`$ donne $`\lambda=1`$. Ainsi
$`W|YX|^{1/2}=|YX|^{1/2}`$ et $`YX=W|YX|=|YX|\ge0`$. Une matrice positive est autoadjointe, donc
$`YX=(YX)^*=XY`$. $`\square`$

D'après le théorème 1.9, pour $`YX=W|YX|`$, l'unitaire $`U=W^*`$ atteint le maximum $`\|YX\|_1`$. Un unitaire
aléatoire fait moins bien en moyenne : pour un $`U`$ distribué selon la loi de Haar
(théorème 5.3), la proposition 11.11 avec $`A=YX`$ donne
$`\big(\mathbb E|\mathrm{Tr}(XUY)|^2\big)^{1/2}=\|YX\|_2/\sqrt k\le\|YX\|_1/\sqrt k=\max_V|\mathrm{Tr}(XVY)|/\sqrt k`$
(car $`\sum_is_i^2\le(\sum_is_i)^2`$), et le corollaire 11.12 montre que, pour $`X,Y`$ de spectres
contenus dans $`[a,b]`$, $`a>0`$, un seul tirage atteint $`\mathrm{Tr}(XY)`$ avec une probabilité au
plus égale à $`b^4/(a^4k^2)`$ (ce qui n'est informatif que pour $`k>b^2/a^2`$).

*Calcul :* la proposition 3.9 pour un couple qui ne commute pas dans $`M_6`$, l'unitaire optimal $`W^*`$ et une recherche aléatoire parmi des unitaires de Haar — [tutoriel, section 3.7](../../tutorial/fr/03_traces_and_states.md#s-3-7).

<a id="s-3-8"></a>

## 3.8 Points de comparaison

*Dans la bibliothèque :* la comparaison avec d'autres logiciels — [tutoriel, section 3.8](../../tutorial/fr/03_traces_and_states.md#s-3-8).

<a id="x-3"></a>

## Exercices

1. Montrer que la densité de l'état tracial $`\tau_{\mathrm{vN}}`$ (poids égaux sur les
   secteurs) relativement à $`\tau_{\mathrm{vN}}`$ est $`1`$, puis calculer sa densité brute et
   son entropie $`S`$.
2. Soit $`\mathrm{tr}_w`$ fidèle. Montrer que $`S_w(\rho)\le\log\mathrm{tr}_w(1)`$, avec
   égalité si et seulement si $`\rho_w=1/\mathrm{tr}_w(1)`$.
3. Montrer que le canal de mesure $`\rho\mapsto\sum_cP_c\rho P_c`$ sur $`M_D`$ (et non sur
   $`\mathcal M`$ !) perturbe bel et bien les états, et calculer $`S(\rho)`$ avant et après la
   mesure pour une superposition de deux secteurs.

<details>
<summary>Solutions</summary>

1. $`\tau_{\mathrm{vN}}(A)=\tau_{\mathrm{vN}}(1\cdot A)`$, donc $`\rho_{\tau}=1`$ par unicité.
   D'après le théorème 3.2, la densité brute est $`\bigoplus_c1_{k_c}/(Ck_c)`$, avec
   $`p_c=1/C`$ et $`\sigma_c=1/k_c`$, donc $`S=\log C+\frac1C\sum_c\log k_c`$ d'après le
   théorème 3.3.

2. Les sommes portent sur les $`c`$ tels que $`p_c>0`$. Alors
   $`S_w(\rho)=S(\rho)+\sum_cp_c\log w_c=-\sum_cp_c\log p_c+\sum_cp_cS(\sigma_c)+\sum_cp_c\log w_c    \le\sum_cp_c\log\frac{k_cw_c}{p_c}\le\log\sum_{p_c>0}k_cw_c\le\log\mathrm{tr}_w(1)`$, en
   utilisant $`S(\sigma_c)\le\log k_c`$ (corollaire 3.6 dans $`M_{k_c}`$) et l'inégalité de Jensen finie
   $`\sum_{c=1}^mp_c\log x_c\le\log\sum_{c=1}^mp_cx_c`$ ($`p_c>0`$, $`\sum_cp_c=1`$, $`x_c>0`$), avec
   égalité seulement si tous les $`x_c`$ sont égaux, appliquée à $`x_c=k_cw_c/p_c`$. Preuve par
   récurrence sur $`m`$ : pour $`m=1`$, il n'y a rien à démontrer. Pour $`m\ge2`$, posons
   $`q=1-p_m\in(0,1)`$ et $`\bar x=\sum_{c<m}(p_c/q)x_c`$ ; par l'hypothèse de récurrence,
   $`\sum_{c<m}(p_c/q)\log x_c\le\log\bar x`$, avec égalité seulement si $`x_1=\dots=x_{m-1}`$, et par
   stricte concavité de $`\log`$ (lemme 3.4 pour $`-\log`$, dont la dérivée seconde est $`1/t^2>0`$)
   $`q\log\bar x+p_m\log x_m\le\log(q\bar x+p_mx_m)`$,
   avec égalité seulement si $`\bar x=x_m`$. En multipliant la première inégalité par $`q`$ et en
   additionnant, $`\sum_cp_c\log x_c\le\log\sum_cp_cx_c`$, et l'égalité impose
   $`x_1=\dots=x_{m-1}=\bar x=x_m`$. L'égalité partout impose que tous les $`p_c`$ soient $`>0`$ (car
   tous les $`k_cw_c>0`$), que
   $`\sigma_c=1/k_c`$ et que $`k_cw_c/p_c`$ soit constant, c'est-à-dire
   $`p_c=k_cw_c/\mathrm{tr}_w(1)`$, soit encore $`(\rho_w)_c=\rho_c/w_c=1_{k_c}/\mathrm{tr}_w(1)`$.
   Réciproquement, pour $`\rho_w=1/\mathrm{tr}_w(1)`$,
   $`S_w(\rho)=-\mathrm{tr}_w(1)\cdot\frac1{\mathrm{tr}_w(1)}\log\frac1{\mathrm{tr}_w(1)}=\log\mathrm{tr}_w(1)`$.

3. Dans $`M_D`$, une superposition pure $`\psi=(\psi_1\oplus\psi_2)/\sqrt2`$ de vecteurs
   unitaires de deux secteurs a pour entropie $`S=0`$ ; après la mesure, l'état est
   $`\frac12|\psi_1\rangle\langle\psi_1|\oplus\frac12|\psi_2\rangle\langle\psi_2|`$, d'entropie
   $`S=\log2`$. Restreint à $`\mathcal M`$, en revanche, l'état vectoriel associé à $`\psi`$
   coïncide avec le mélange :
   $`\langle\psi,A\psi\rangle=\tfrac12\langle\psi_1,A_1\psi_1\rangle+\tfrac12\langle\psi_2,A_2\psi_2\rangle`$
   pour tout $`A\in\mathcal M`$. Aucun élément de $`\mathcal M`$ ne détecte la cohérence
   (physiquement, c'est le contenu d'une règle de supersélection).

</details>

*Calcul :* l'exercice 1, numériquement — [tutoriel, Exercices](../../tutorial/fr/03_traces_and_states.md#x-3).

**Suite :** le [chapitre 4](04_symmetries.md) traite des symétries et des représentations :
comment l'action d'un groupe sur un système physique fait apparaître les secteurs, et
comment la bibliothèque construit les algèbres correspondantes pour les systèmes de spins
et les registres de qubits.
