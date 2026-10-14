[English](../11_optimisation.md) · [Español](../es/11_optimisation.md) · **Français** · [Deutsch](../de/11_optimisation.md) · [中文](../zh/11_optimisation.md) · [日本語](../ja/11_optimisation.md) · [Русский](../ru/11_optimisation.md)

# 11. Optimisation et inégalités de traces

Une inégalité de traces est un énoncé sur les valeurs extrêmes d'une fonctionnelle d'opérateurs. La section 3.7 en
a présenté une qui est plausible mais fausse (proposition 3.9) ; un calcul par échantillonnage ([tutoriel, section 3.7](../../tutorial/fr/03_traces_and_states.md#s-3-7)) n'a pas trouvé numériquement la violation. Ce chapitre explique pourquoi l'échantillonnage échoue en grande dimension et le remplace par une optimisation sur
des ensembles de contraintes paramétrés de façon exacte ; lorsque les extrêmes sont connus sous forme close, nous
les démontrons, de sorte que l'optimiseur est confronté à un théorème.

**Vous apprendrez** comment on paramètre les unitaires et les spectres de contraste prescrit, et
pourquoi tout opérateur admissible est atteint ; l'inégalité de traces de von Neumann et les extrêmes exacts
qu'elle fournit ; une inégalité de traces avec son cas d'égalité ; pourquoi les échantillons de Haar sont, en moyenne
quadratique, loin des extrêmes ; et ce que démontre une recherche multi-départs.
Les vérifications numériques de ce chapitre, ainsi que l'optimiseur de la bibliothèque, se trouvent dans le
[Tutoriel 11](../../tutorial/fr/11_optimisation.md).

Dans tout ce chapitre, $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$ avec $`D=\sum_ck_c\ge2`$, $`\mathrm{Tr}`$ est
la trace brute (section 3.1) et $`\|A\|_1=\sum_c\|A_c\|_1`$, la norme trace de $`A`$ vue comme matrice diagonale par blocs
dans $`M_D(\mathbb C)`$. Les valeurs singulières de $`A\in M_k(\mathbb C)`$ sont $`s_1(A)\ge\dots\ge s_k(A)`$ ;
$`x^\downarrow`$ et $`x^\uparrow`$ désignent les réarrangements décroissant et croissant d'un vecteur réel $`x`$. $`\Delta`$
est le contraste de Michelson (section 3.6), et pour $`\delta\in[0,1)`$
```math
\ell(\delta)=\frac{1-\delta}{1+\delta},\qquad 1-\ell(\delta)=\frac{2\delta}{1+\delta};
```
un opérateur positif avec $`\lambda_{\max}=1`$ a pour contraste $`\ell(\lambda_{\min})`$ ; il a donc le contraste $`\delta`$
si et seulement si $`\lambda_{\min}=\ell(\delta)`$, puisque $`\ell`$ est une involution strictement décroissante de
$`[0,1]`$. Les exemples sont les trois fonctionnelles de l'article qui accompagne la bibliothèque (section *From
sampling to optimisation*) :

| | fonctionnelle | variables |
|---|---|---|
| Exp. 1 | $`z=\lvert\mathrm{Tr}(XUY)\rvert-\mathrm{Tr}(XY)`$ | $`X,Y\ge0`$, $`U\in\mathcal M`$ unitaire |
| Exp. 2 | $`z=\mathrm{Tr}\lvert XY\rvert-\mathrm{Tr}(\lvert X\rvert\,\lvert Y\rvert)`$ | $`X=UX_0`$, $`Y=VY_0`$, $`X_0,Y_0\ge0`$, $`U,V`$ unitaires |
| Exp. 3 | $`z=\mathrm{Tr}(Y\lvert X\rvert Y)-\mathrm{Tr}\lvert YXY\rvert`$ | $`X=X^*`$, $`Y\ge0`$ |

avec $`|X|`$, $`|Y|`$ normalisés par $`\lambda_{\max}=1`$ et leurs contrastes comme paramètres. Les **enveloppes** sont
$`\overline z(\delta)=\sup\{z:\Delta(|X|)=\delta\}`$ et $`\underline z(\delta)=\inf\{z:\Delta(|X|)=\delta\}`$, toutes les
autres variables étant libres.

<a id="s-11-1"></a>

## 11.1 Ensembles de contraintes et leurs paramétrages

**Définition.** Pour $`\delta\in[0,1)`$ et $`\ell=\ell(\delta)`$, la **couche de contraste** est
```math
\mathcal P_\delta=\{X\in\mathcal M:X\ge0,\ \lambda_{\max}(X)=1,\ \Delta(X)=\delta\}
=\{X=X^*:\mathrm{spec}\,X\subseteq[\ell,1],\ 1,\ell\in\mathrm{spec}\,X\}.
```
Les **emplacements** de $`\mathcal M`$ sont les $`D`$ couples $`(c,i)`$, $`i\le k_c`$ ; pour $`\lambda\in\mathbb R^D`$,
$`\mathrm{diag}(\lambda)=\bigoplus_c\mathrm{diag}(\lambda_{c,1},\dots,\lambda_{c,k_c})`$.

**Le paramétrage.** Pour un contraste prescrit $`\delta`$, une recherche tire pour chaque départ deux emplacements distincts
$`s_\pm=(c_\pm,i_\pm)`$ ($`s_+`$ uniformément parmi les $`D`$ emplacements, puis $`s_-`$ uniformément parmi les $`D-1`$ autres,
de sorte que chaque couple ordonné d'emplacements distincts a la probabilité $`1/(D(D-1))`$) et un unitaire
$`U_0\in\mathcal M`$ selon la loi de Haar (lemme 11.10), et envoie $`\theta\in\mathbb R^D`$ et
$`G=\bigoplus_cG_c`$ sur
```math
X(\theta,G)=U\mathrm{diag}(\lambda)U^*,\qquad U=U_0\,e^{G-G^*},\qquad
\lambda_{s_+}=1,\quad\lambda_{s_-}=\ell,\quad\lambda_s=\ell+(1-\ell)\tfrac{1+\sin\theta_s}2\ \text{otherwise}.
```
Les unitaires seuls sont paramétrés par $`G\mapsto U_0e^{G-G^*}`$. Contrairement à une sigmoïde, $`(1+\sin\theta)/2`$
atteint $`0`$ et $`1`$ pour des $`\theta`$ finis, ce qui est commode car les spectres extrémaux de ce chapitre
(corollaire 11.7, proposition 11.14) se trouvent en $`\ell`$ et $`1`$. Un contraste libre dans $`[a,b]`$ est un paramètre de
plus, $`\delta=a+(b-a)(1+\sin\eta)/2`$, par exemple avec $`[a,b]=[0,0.999]`$. Une valeur atteinte avec cette restriction l'est aussi sans elle, et reste donc une borne
inférieure pour un supremum non restreint (une borne supérieure pour un infimum non restreint).

*Dans la bibliothèque :* ce paramétrage est `optimize.PositiveParam`, et $`G\mapsto U_0e^{G-G^*}`$ est `optimize.UnitaryParam` — [tutoriel, section 11.1](../../tutorial/fr/11_optimisation.md#s-11-1).

<a id="r-11-1"></a>

**Lemme 11.1.** Pour tout unitaire $`U_0\in\mathcal M`$, l'application $`G\mapsto U_0e^{G-G^*}`$ est surjective sur le
groupe unitaire de $`\mathcal M`$.

*Démonstration.* Soit $`U`$ unitaire. D'après le théorème 1.3 (bloc par bloc) et le corollaire 1.4(2),
$`U_0^*U=\bigoplus_cW_c\mathrm{diag}(e^{i\varphi_{c,j}})_jW_c^*`$ avec $`W_c`$ unitaire et $`\varphi_{c,j}\in\mathbb R`$.
Alors $`S=\bigoplus_cW_c\mathrm{diag}(i\varphi_{c,j})_jW_c^*`$ vérifie $`S^*=-S`$ et $`e^S=U_0^*U`$ (la série
$`\sum_nS^n/n!`$ se calcule terme à terme dans la base propre), et $`G=S/2`$ donne $`G-G^*=S`$. $`\square`$

<a id="r-11-2"></a>

**Proposition 11.2 (le paramétrage est surjectif).** Pour $`\delta`$, $`s_\pm`$ et $`U_0`$ fixés, l'image de
$`(\theta,G)\mapsto X(\theta,G)`$ est $`\{X\in\mathcal P_\delta:1\in\mathrm{spec}\,X_{c_+},\ \ell\in\mathrm{spec}\,X_{c_-}\}`$,
et $`\mathcal P_\delta`$ est la réunion de ces images sur tous les couples d'emplacements distincts.

*Démonstration.* L'image est contenue dans cet ensemble : $`X`$ a des valeurs propres $`\lambda_s\in[\ell,1]`$, avec
$`\lambda_{s_\pm}`$ dans les blocs $`c_\pm`$. Réciproquement, diagonalisons $`X_c=W_c\mathrm{diag}(\mu_c)W_c^*`$ et
permutons chaque $`\mu_c`$ (en remplaçant $`W_c`$ par $`W_c\Pi_c`$, où $`\Pi_c`$ est une matrice de permutation) de sorte
que $`\mu_{s_+}=1`$ et $`\mu_{s_-}=\ell`$. C'est possible : si $`c_+=c_-`$ et $`\delta>0`$, alors $`1\ne\ell`$ sont deux valeurs
propres du même bloc et occupent des positions différentes ; si $`\delta=0`$, toutes les valeurs propres valent $`1`$.
Si $`\ell<1`$, posons $`\theta_s=\arcsin\big(2(\mu_s-\ell)/(1-\ell)-1\big)`$ pour les autres emplacements, de sorte que
$`\lambda=\mu`$ (si $`\ell=1`$, tout $`\theta`$ convient). Le lemme 11.1 fournit $`G`$ avec $`U_0e^{G-G^*}=\bigoplus_cW_c`$.
*Réunion :* pour $`X\in\mathcal P_\delta`$, choisissons des blocs $`c_\pm`$ contenant les valeurs propres $`1`$, $`\ell`$ et
des emplacements distincts dans ces blocs ; si $`c_+=c_-`$ et $`\delta>0`$, le bloc a deux valeurs propres différentes,
donc $`k_{c_+}\ge2`$ ; si $`\delta=0`$, $`X=1`$ et deux emplacements quelconques conviennent. $`\square`$

*Calcul :* la proposition 11.2 de façon constructive : les paramètres de la démonstration reproduisent un opérateur cible — [tutoriel, section 11.1](../../tutorial/fr/11_optimisation.md#s-11-1).

Comme chaque couple d'emplacements est tiré avec une probabilité positive, chaque point de la couche se trouve, avec
une probabilité positive, dans l'image d'un départ, et (en arithmétique exacte) la contrainte est satisfaite à chaque
itéré. Pour $`X`$ autoadjoint (Exp. 3), on utilise $`X=U\mathrm{diag}(\sigma\lambda)U^*`$ avec des signes $`\sigma\in\{\pm1\}^D`$
tirés pour chaque départ (dans chaque bloc, le nombre de signes moins est uniforme dans $`\{0,\dots,k_c\}`$ et leurs
positions sont uniformément aléatoires, de sorte que chaque configuration de signes a une probabilité positive) et
modifiés uniquement par les changements de signe de la section 11.5. Alors $`|X|=U\mathrm{diag}(\lambda)U^*`$
(théorème 1.5(e), car $`\lambda>0`$), et la réunion sur les signes et les emplacements est
$`\{X=X^*:|X|\in\mathcal P_\delta\}`$ : pour un tel $`X`$, diagonalisons $`X_c=W_c\mathrm{diag}(\nu_c)W_c^*`$, de sorte
que $`|X|_c=W_c\mathrm{diag}(|\nu_c|)W_c^*`$ avec $`|\nu_{c,i}|\ge\ell>0`$, et appliquons la démonstration de la
proposition 11.2 à $`\mu=|\nu|`$ avec $`\sigma_s=\nu_s/|\nu_s|`$. Seul $`S=G-G^*`$ compte.

*Dans la bibliothèque :* `optimize.SelfAdjointParam` — [tutoriel, section 11.1](../../tutorial/fr/11_optimisation.md#s-11-1).

La carte $`S\mapsto U_0e^S`$ ne se comporte bien qu'au voisinage de $`S=0`$. Soient $`\mathfrak u=\{S\in\mathcal M:S^*=-S\}`$
et $`\mathcal M_h=\{H\in\mathcal M:H^*=H\}`$, sous-espaces réels avec $`\mathcal M=\mathfrak u\oplus\mathcal M_h`$
($`X=\frac{X-X^*}2+\frac{X+X^*}2`$). *Au voisinage de $`0`$, $`S\mapsto U_0e^S`$ est un homéomorphisme d'une boule ouverte
$`N\ni0`$ de $`\mathfrak u`$ sur un voisinage de $`U_0`$ dans le groupe unitaire.* En effet, $`\Psi(S,H)=U_0e^Se^H`$ est une
application $`C^1`$ $`\mathfrak u\times\mathcal M_h\to\mathcal M`$ (la série $`\sum_{m,n\ge0}U_0S^mH^n/(m!\,n!)`$ et la série
des dérivées partielles de ses termes par rapport aux coordonnées réelles convergent uniformément sur les parties
bornées), et $`e^{tS}e^{tH}=1+t(S+H)+O(t^2)`$ montre que sa différentielle en $`0`$ est $`(S,H)\mapsto U_0(S+H)`$, une
bijection sur $`\mathcal M`$. D'après le théorème d'inversion locale, il existe des boules ouvertes $`N\ni0`$ dans
$`\mathfrak u`$ et $`N_h\ni0`$ dans $`\mathcal M_h`$ telles que $`\Psi`$ soit un homéomorphisme de $`N\times N_h`$ sur un ouvert
$`O\subseteq\mathcal M`$. Si $`\Psi(S,H)`$ est unitaire, alors $`e^H=e^{-S}U_0^*\Psi(S,H)`$ l'est aussi (dans une base propre
de la matrice normale $`S`$, dont les valeurs propres sont imaginaires, $`e^{\pm S}`$ sont diagonales, à coefficients de
module $`1`$, et inverses l'une de l'autre) ; mais $`e^H`$ a les valeurs propres positives $`e^{h_i}`$, où les $`h_i`$ sont
les valeurs propres de $`H`$, donc $`e^{h_i}=1`$ (corollaire 1.4(2)) et $`H=0`$. Par conséquent $`O\cap U(\mathcal M)=U_0e^N`$,
où $`U(\mathcal M)`$ est le groupe unitaire de $`\mathcal M`$ : c'est un voisinage de $`U_0`$ dans $`U(\mathcal M)`$, et
$`S\mapsto U_0e^S`$ est une bijection continue $`N\to U_0e^N`$ dont la réciproque, la première composante de
$`\Psi^{-1}`$, est continue. Loin de $`0`$, la carte n'est pas injective : $`e^{2\pi iP}=1=e^0`$ pour tout projecteur
$`P\in\mathcal M`$. C'est pourquoi une recherche **recentre** la carte après chaque tour ($`U_0\leftarrow U_0e^{G-G^*}`$,
$`G\leftarrow0`$). Deux autres cartes ont la même différentielle $`S\mapsto U_0S`$ en $`0`$. La carte de Cayley
$`S\mapsto U_0C(-S/2)`$ (exercice 1) : $`(1+tS/2)(1-tS/2)^{-1}=(1+tS/2)(1+tS/2+O(t^2))=1+tS+O(t^2)`$. La rétraction
polaire, $`S\mapsto`$ le facteur unitaire de $`A=U_0(1+S)`$ (théorème 1.8) : $`A^*A=(1-S)(1+S)=1+S^*S\ge1`$ est inversible,
donc $`|A|`$ est inversible et le facteur unitaire est $`A|A|^{-1}`$ ; avec $`tS`$ au lieu de $`S`$,
$`|A|^{-1}=(1+t^2S^*S)^{-1/2}=1+O(t^2)`$ (dans une base propre de $`S^*S`$), donc le facteur vaut $`U_0(1+tS)+O(t^2)`$.

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Le contraste de Michelson est
> la visibilité $`(I_{\max}-I_{\min})/(I_{\max}+I_{\min})`$ de franges d'interférence. Pour $`0\le X\le1`$ avec
> $`\lambda_{\max}=1`$ — un filtre, un atténuateur, l'efficacité d'un détecteur —, il fixe la dynamique
> $`\lambda_{\min}/\lambda_{\max}`$, l'inverse d'un conditionnement.

<a id="s-11-2"></a>

## 11.2 L'inégalité de traces de von Neumann

Une matrice $`(c_{ij})`$ avec $`c_{ij}\ge0`$ dont toutes les sommes de lignes et de colonnes valent $`1`$ est dite
**bistochastique** ; pour un unitaire $`G`$, $`(|G_{ij}|^2)`$ et $`(|G_{ji}|^2)`$ sont bistochastiques, car les lignes et
les colonnes de $`G`$ sont des vecteurs unitaires.

<a id="r-11-3"></a>

**Lemme 11.3 (réarrangement).** Soient $`p_1\ge\dots\ge p_k`$, $`q_1\ge\dots\ge q_k`$ réels et $`(c_{ij})`$
bistochastique. Alors $`\sum_ip_iq_{k+1-i}\le\sum_{i,j}p_iq_jc_{ij}\le\sum_ip_iq_i`$.

*Démonstration.* *Majoration.* Pour $`m\le k`$, soit $`r_j=\sum_{i\le m}c_{ij}\in[0,1]`$, de sorte que $`\sum_jr_j=m`$.
Comme $`q_j\ge q_m`$ pour $`j\le m`$ et $`q_j\le q_m`$ pour $`j>m`$,
```math
\sum_{j\le m}q_j-\sum_jq_jr_j=\sum_{j\le m}q_j(1-r_j)-\sum_{j>m}q_jr_j\ge q_m\Big(m-\sum_jr_j\Big)=0 .
```
Avec $`x_i=\sum_jq_jc_{ij}`$, cela signifie $`\sum_{i\le m}x_i\le\sum_{j\le m}q_j`$, avec égalité pour $`m=k`$. Par
sommation d'Abel et $`p_m-p_{m+1}\ge0`$,
```math
\sum_ip_ix_i=p_k\sum_{i\le k}x_i+\sum_{m<k}(p_m-p_{m+1})\sum_{i\le m}x_i
\le p_k\sum_{j\le k}q_j+\sum_{m<k}(p_m-p_{m+1})\sum_{j\le m}q_j=\sum_ip_iq_i .
```
*Minoration.* Appliquons la majoration à $`p`$, à la suite décroissante $`q'_j=-q_{k+1-j}`$ et à la matrice
bistochastique $`c'_{ij}=c_{i,k+1-j}`$ :
$`-\sum_{ij}p_iq_jc_{ij}=\sum_{ij}p_iq'_jc'_{ij}\le\sum_ip_iq'_i=-\sum_ip_iq_{k+1-i}`$. $`\square`$

Toute $`A\in M_k(\mathbb C)`$ admet une **décomposition en valeurs singulières** $`A=W\Sigma V^*`$, avec $`W,V`$ unitaires
et $`\Sigma=\mathrm{diag}(s_i(A))`$ : écrivons $`A=W'|A|`$ (théorème 1.8), $`|A|=V\Sigma V^*`$ (théorème 1.3) et
$`W=W'V`$.

<a id="r-11-4"></a>

**Théorème 11.4 (inégalité de traces de von Neumann).** Pour $`A,B\in M_k(\mathbb C)`$,
$`|\mathrm{Tr}(AB)|\le\sum_is_i(A)\,s_i(B)`$.

*Démonstration.* Soient $`A=W_A\Sigma_AV_A^*`$, $`B=W_B\Sigma_BV_B^*`$, $`a_i=s_i(A)`$, $`b_j=s_j(B)`$, et $`G=V_A^*W_B`$,
$`H=V_B^*W_A`$. Par cyclicité, $`\mathrm{Tr}(AB)=\mathrm{Tr}(\Sigma_AG\Sigma_BH)=\sum_{i,j}a_iG_{ij}b_jH_{ji}`$,
donc, d'après $`|g||h|\le\frac12(|g|^2+|h|^2)`$ et $`a_ib_j\ge0`$,
```math
|\mathrm{Tr}(AB)|\le\sum_{i,j}a_ib_j\,d_{ij},\qquad d_{ij}=\tfrac12\big(|G_{ij}|^2+|H_{ji}|^2\big).
```
$`(d_{ij})`$ est bistochastique, et le lemme 11.3 donne $`\sum_{ij}a_ib_jd_{ij}\le\sum_ia_ib_i`$. $`\square`$

<a id="r-11-5"></a>

**Corollaire 11.5.** (a) Pour $`P,Q\in M_k(\mathbb C)`$ autoadjointes de valeurs propres $`p,q`$ :
$`\sum_ip^\downarrow_iq^\uparrow_i\le\mathrm{Tr}(PQ)\le\sum_ip^\downarrow_iq^\downarrow_i`$.
(b) Pour $`P,Q\ge0`$ : $`\mathrm{Tr}(PQ)\le\|PQ\|_1\le\sum_ip^\downarrow_iq^\downarrow_i`$.

*Démonstration.* (a) Comme dans la démonstration du théorème 3.5, $`\mathrm{Tr}(PQ)=\sum_{ij}p_iq_jc_{ij}`$ avec la
matrice bistochastique $`c_{ij}=|\langle u_i,v_j\rangle|^2`$ pour des bases propres ordonnées de façon décroissante ;
appliquons le lemme 11.3 (sa minoration est $`\sum_ip^\downarrow_iq^\downarrow_{k+1-i}=\sum_ip^\downarrow_iq^\uparrow_i`$).
(b) La première inégalité est le théorème 1.9. D'après le théorème 1.9 et le théorème 11.4,
$`\|PQ\|_1=\max_U|\mathrm{Tr}(P\,QU)|\le\sum_is_i(P)s_i(QU)`$, et $`s_i(QU)=q^\downarrow_i`$ car
$`|QU|=(U^*Q^2U)^{1/2}=U^*QU`$ (théorème 1.6). $`\square`$

*Calcul :* le théorème 11.4 et le corollaire 11.5 sur des matrices aléatoires — [tutoriel, section 11.2](../../tutorial/fr/11_optimisation.md#s-11-2).

Ainsi, pour $`X,Y\ge0`$ dans $`M_k(\mathbb C)`$ de spectres $`x,y`$, la violation $`\|YX\|_1-\mathrm{Tr}(XY)`$ de la
proposition 3.9 est au plus $`\sum_ix^\downarrow_i(y^\downarrow_i-y^\uparrow_i)`$ (corollaire 11.5(b) pour $`(Y,X)`$ et la
minoration de (a)). Pour Exp. 2, ce nombre est la réponse exacte :

<a id="r-11-6"></a>

**Proposition 11.6 (Exp. 2 à spectres fixés).** Dans Exp. 2, $`|X|=X_0`$, $`|Y|=Y_0`$ et
$`z=\|X_0VY_0\|_1-\mathrm{Tr}(X_0Y_0)`$. Si les spectres $`x_c`$ de $`(X_0)_c`$ et $`y_c`$ de $`(Y_0)_c`$ sont fixés et
que les bases propres et $`V`$ varient, alors
```math
\sup z=-\inf z=\sum_c\sum_ix^\downarrow_{c,i}\big(y^\downarrow_{c,i}-y^\uparrow_{c,i}\big),
```
et les deux sont atteints.

*Démonstration.* $`X^*X=X_0U^*UX_0=X_0^2`$, donc $`|X|=X_0`$ (théorème 1.6) ; de même $`|Y|=Y_0`$, et $`|UA|=|A|`$ donne
$`\mathrm{Tr}|XY|=\|X_0VY_0\|_1`$. Toutes les variables agissent bloc par bloc et $`z`$ est une somme sur les blocs ;
on peut donc supposer $`C=1`$. Posons $`Y'=VY_0V^*\ge0`$ (de spectre $`y`$) ; comme $`|AV|=V^*|A|V`$,
$`\|X_0VY_0\|_1=\|X_0Y'V\|_1=\|X_0Y'\|_1`$. Le corollaire 11.5 appliqué aux couples $`(X_0,Y')`$ et $`(X_0,Y_0)`$, avec
$`\|X_0Y'\|_1\ge\mathrm{Tr}(X_0Y')`$, donne
```math
\sum_ix^\downarrow_iy^\uparrow_i-\sum_ix^\downarrow_iy^\downarrow_i\le\|X_0Y'\|_1-\mathrm{Tr}(X_0Y_0)\le\sum_ix^\downarrow_iy^\downarrow_i-\sum_ix^\downarrow_iy^\uparrow_i .
```
Soit $`J`$ la permutation de retournement, $`Je_i=e_{k+1-i}`$ ; alors $`J\mathrm{diag}(w)=\mathrm{diag}(w_k,\dots,w_1)J`$,
et la norme trace de $`\mathrm{diag}(w)J`$ vaut $`\sum_i|w_i|`$. Avec $`X_0=\mathrm{diag}(x^\downarrow)`$ et $`V=J`$ :
pour $`Y_0=\mathrm{diag}(y^\uparrow)`$, $`X_0VY_0=\mathrm{diag}(x^\downarrow_iy^\downarrow_i)J`$ et
$`\mathrm{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\uparrow_i`$ (la borne supérieure) ; pour $`Y_0=\mathrm{diag}(y^\downarrow)`$,
$`X_0VY_0=\mathrm{diag}(x^\downarrow_iy^\uparrow_i)J`$ et $`\mathrm{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\downarrow_i`$
(la borne inférieure). $`\square`$

<a id="r-11-7"></a>

**Corollaire 11.7 (enveloppe d'Exp. 2).** Soient $`\ell_X,\ell_Y\in(0,1]`$, et supposons $`k_c\ge2`$ pour au moins un $`c`$. Sur tous les
$`X_0,Y_0\ge0`$ avec $`\lambda_{\max}=1`$, $`\lambda_{\min}(X_0)=\ell_X`$, $`\lambda_{\min}(Y_0)=\ell_Y`$, et tous les unitaires,
```math
\sup z=-\inf z=\sum_c\Big\lfloor\frac{k_c}2\Big\rfloor(1-\ell_X)(1-\ell_Y).
```

*Démonstration.* Dans un bloc, écrivons $`x_i=x^\downarrow_{c,i}`$, $`y_i=y^\downarrow_{c,i}`$, $`k=k_c`$, de sorte que
$`y^\uparrow_{c,i}=y_{k+1-i}`$. En associant $`i`$ à $`k+1-i`$,
$`x_i(y_i-y_{k+1-i})+x_{k+1-i}(y_{k+1-i}-y_i)=(x_i-x_{k+1-i})(y_i-y_{k+1-i})`$, et pour $`k`$ impair l'indice du milieu
contribue $`x_i(y_i-y_i)=0`$. Par conséquent
```math
\sum_ix_i(y_i-y_{k+1-i})=\sum_{i\le\lfloor k/2\rfloor}(x_i-x_{k+1-i})(y_i-y_{k+1-i})\le\Big\lfloor\frac k2\Big\rfloor(1-\ell_X)(1-\ell_Y),
```
puisque $`x_i\in[\ell_X,1]`$, $`y_i\in[\ell_Y,1]`$ sont décroissants, de sorte que les deux facteurs sont dans
$`[0,1-\ell_X]`$ et $`[0,1-\ell_Y]`$. L'égalité a lieu si, dans chaque bloc, les $`\lfloor k_c/2\rfloor`$ plus grandes
valeurs propres valent $`1`$ et les $`\lfloor k_c/2\rfloor`$ plus petites valent $`\ell_X`$ (resp. $`\ell_Y`$), celle du milieu
($`k_c`$ impair) étant arbitraire dans $`[\ell_X,1]`$ (resp. $`[\ell_Y,1]`$) ; ces spectres ont pour extrêmes globaux $`1`$ et $`\ell_X`$ ($`\ell_Y`$) parce que $`k_c\ge2`$ pour au moins un $`c`$.
La proposition 11.6 permet de conclure. $`\square`$

Ici $`1-\ell_X=2\delta/(1+\delta)`$ pour $`\Delta(|X|)=\delta`$ ; pour $`k_c`$ impair, la valeur propre du milieu est sans
importance.

*Calcul :* les configurations extrémales de la proposition 11.6 et l'enveloppe du corollaire 11.7 pour $`k=3`$, $`C=2`$ — [tutoriel, section 11.2](../../tutorial/fr/11_optimisation.md#s-11-2) ; l'enveloppe trouvée par une recherche multi-départs — [tutoriel, section 11.5](../../tutorial/fr/11_optimisation.md#s-11-5).

<a id="s-11-3"></a>

## 11.3 Une inégalité de traces avec cas d'égalité

Pour $`X=X^*`$, soit $`X_\pm=f_\pm(X)`$ avec $`f_+(t)=\max(t,0)`$, $`f_-(t)=\max(-t,0)`$. Alors $`X_\pm\ge0`$, $`X=X_+-X_-`$,
$`|X|=X_++X_-`$ et $`X_+X_-=0`$ (théorème 1.5(c), et 1.5(e) pour $`|X|=(X^2)^{1/2}=h(X)`$ avec
$`h(t)=\sqrt{t^2}=|t|=f_+(t)+f_-(t)`$ ; corollaire 1.4(3) pour la positivité).

<a id="r-11-8"></a>

**Lemme 11.8.** Pour $`A,B\ge0`$ dans $`M_n(\mathbb C)`$, $`\|A-B\|_1\le\mathrm{Tr}\,A+\mathrm{Tr}\,B`$, avec égalité si et
seulement si $`AB=0`$.

*Démonstration.* Soient $`H=A-B`$ et $`Q=g(H)`$ avec $`g=1`$ sur $`(0,\infty)`$ et $`g=0`$ ailleurs ; $`Q`$ est un projecteur,
car $`g=\bar g=g^2`$ (théorème 1.5(c)). Comme $`|t|=t(2g(t)-1)`$, on a $`|H|=H(2Q-1)`$ (théorème 1.5(c),(e)) et, en utilisant
$`\mathrm{Tr}(AR)=\mathrm{Tr}(RAR)`$ pour les projecteurs $`R=Q,1-Q`$ (cyclicité),
```math
\|H\|_1=\mathrm{Tr}\big(A(2Q-1)\big)-\mathrm{Tr}\big(B(2Q-1)\big)
=\mathrm{Tr}\,A+\mathrm{Tr}\,B-2\mathrm{Tr}\big((1-Q)A(1-Q)\big)-2\mathrm{Tr}(QBQ).
```
Les traces soustraites sont $`\ge0`$, et $`\mathrm{Tr}(RAR)=\|A^{1/2}R\|_2^2`$ s'annule si et seulement si $`AR=0`$ ;
l'égalité a donc lieu si et seulement si $`A(1-Q)=0`$ et $`BQ=0`$, et alors $`AB=AQB=A(BQ)^*=0`$. Réciproquement, si
$`AB=0`$, alors $`BA=0`$, $`(A-B)^2=A^2+B^2=(A+B)^2`$, et le théorème 1.6 donne $`|A-B|=A+B`$. $`\square`$

<a id="r-11-9"></a>

**Proposition 11.9 (Exp. 3).** Pour $`X=X^*`$ et $`Y\ge0`$ dans $`\mathcal M`$,
$`\mathrm{Tr}|YXY|\le\mathrm{Tr}(Y|X|Y)`$, avec égalité si et seulement si $`YX_+Y^2X_-Y=0`$, en particulier si
$`XY=YX`$.

*Démonstration.* Dans $`M_D(\mathbb C)`$, $`A=YX_+Y`$ et $`B=YX_-Y`$ sont positifs (corollaire 1.7 avec $`T=X_\pm^{1/2}Y`$),
$`YXY=A-B`$ et $`\mathrm{Tr}(Y|X|Y)=\mathrm{Tr}\,A+\mathrm{Tr}\,B`$ ; appliquons le lemme 11.8. Si $`XY=YX`$, alors
$`Y`$ commute avec $`X_\pm`$ (théorème 1.5(d)), et $`AB=Y^4X_+X_-=0`$. $`\square`$

Ainsi $`z\ge0`$ dans Exp. 3, et $`\underline z=0`$ (prendre $`Y=1`$). La condition d'égalité est plus faible que la
commutation : pour $`X=\mathrm{diag}(1,-1,0)`$ et $`Y`$ de lignes $`(1,0,a)`$, $`(0,1,0)`$, $`(a,0,1)`$, $`|a|<1`$, on a
$`AB=(Ye_1)\langle Ye_1,Ye_2\rangle(Ye_2)^*=0`$, alors que $`XY\ne YX`$ pour $`a\ne0`$.

*Calcul :* la proposition 11.9 sur des couples aléatoires, et l'exemple d'égalité sans commutation — [tutoriel, section 11.3](../../tutorial/fr/11_optimisation.md#s-11-3).

<a id="s-11-4"></a>

## 11.4 Échantillonnage contre optimisation

Quelle est la taille typique de $`|\mathrm{Tr}(AU)|`$, comparée à son maximum $`\|A\|_1`$ ? Nous tirons $`U`$
comme le facteur $`Q`$ de la décomposition QR $`Z=QR`$, $`R_{ii}>0`$, d'une **matrice de Ginibre** $`Z`$ (coefficients
indépendants de densité $`\pi^{-1}e^{-|z|^2}`$ ; recette de Mezzadri), indépendamment dans chaque bloc.

<a id="r-11-10"></a>

**Lemme 11.10.** $`Q`$ est distribué selon la loi de Haar sur $`U(k)`$ ; en particulier, pour tout unitaire fixé
$`V\in M_k(\mathbb C)`$, $`VQ`$ a la même loi que $`Q`$. L'unitaire $`\bigoplus_cQ_c`$ à blocs indépendants est distribué
selon la loi de Haar sur le groupe unitaire de $`\mathcal M`$.

*Démonstration.* C'est le théorème 5.3 (l'invariance à gauche est l'étape d'équivariance de sa démonstration ; la
loi est identifiée à la mesure de Haar par le lemme d'unicité 5.1, dont l'étape de théorie de la mesure est le
corollaire A.2) et le corollaire 5.4. $`\square`$

*Dans la bibliothèque :* cet échantillonneur est `TypeIAlgebra.random_unitary_operator` — [tutoriel, section 11.4](../../tutorial/fr/11_optimisation.md#s-11-4).

<a id="r-11-11"></a>

**Proposition 11.11.** Pour $`U\in U(k)`$ de loi de Haar et $`A\in M_k(\mathbb C)`$ : $`\mathbb E\mathrm{Tr}(AU)=0`$ et
$`\mathbb E|\mathrm{Tr}(AU)|^2=\|A\|_2^2/k`$.

*Démonstration.* $`\mathbb EU_{ij}=0`$ d'après le corollaire 5.5. La proposition 5.7 (avec $`C=1`$) pour
$`X=|e_j\rangle\langle e_m|`$ donne la formule de Weingarten d'ordre un,
$`\mathbb E[U_{ij}\bar U_{lm}]=\mathbb E(UXU^*)_{il}=\delta_{il}\delta_{jm}/k`$ ; pour $`k\ge2`$, elle découle aussi de la
formule d'ordre deux, le théorème 5.10, en posant $`i_2=i'_2`$, $`j_2=j'_2`$ et en sommant sur les deux indices, puisque
$`\sum_{i_2,j_2}|U_{i_2j_2}|^2=k`$. Avec $`\mathrm{Tr}(AU)=\sum_{ij}A_{ji}U_{ij}`$,
$`\mathbb E|\mathrm{Tr}(AU)|^2=\sum_{ijlm}A_{ji}\bar A_{ml}\mathbb E[U_{ij}\bar U_{lm}]=\frac1k\sum_{ij}|A_{ji}|^2`$. $`\square`$

<a id="r-11-12"></a>

**Corollaire 11.12.** Soient $`X,Y\ge0`$ dans $`M_k(\mathbb C)`$ à spectres dans $`[a,b]`$, $`a>0`$, et $`U`$ de loi de Haar.
Alors
```math
\frac{\big(\mathbb E|\mathrm{Tr}(XUY)|^2\big)^{1/2}}{\max_V|\mathrm{Tr}(XVY)|}=\frac{\|YX\|_2}{\sqrt k\,\|YX\|_1}\le\min\Big(\frac1{\sqrt k},\frac{b^2}{a^2k}\Big),
\qquad
\mathbb P\big(|\mathrm{Tr}(XUY)|\ge\mathrm{Tr}(XY)\big)\le\frac{b^4}{a^4k^2}.
```

*Démonstration.* L'égalité résulte de la proposition 11.11 pour $`A=YX`$ (car
$`\mathrm{Tr}(XUY)=\mathrm{Tr}(AU)`$) et de la proposition 3.9. Comme
$`\|A\|_2^2=\sum_is_i^2\le(\sum_is_i)^2=\|A\|_1^2`$ (section 1.4), le rapport est $`\le1/\sqrt k`$. Comme
$`\|Av\|=\||A|v\|`$ (section 1.4), en diagonalisant $`|A|`$ on obtient $`s_k(A)=\min_{\|v\|=1}\|Av\|`$ et
$`s_1(A)=\max_{\|v\|=1}\|Av\|`$ ; et $`\|Xv\|^2=\langle v,X^2v\rangle\in[a^2,b^2]`$ pour $`v`$ unitaire (corollaire 1.4(3)
pour $`X^2-a^2`$ et $`b^2-X^2`$), de même pour $`Y`$. Par conséquent $`a^2\le s_k(A)\le s_1(A)\le b^2`$, $`\|A\|_2^2\le kb^4`$ et
$`\|A\|_1\ge ka^2`$. Enfin
$`\mathrm{Tr}(XY)=\mathrm{Tr}(X^{1/2}YX^{1/2})\ge a\mathrm{Tr}\,X\ge ka^2`$, et l'inégalité de Markov pour
$`|\mathrm{Tr}(XUY)|^2`$ majore la probabilité par $`(kb^4/k)/(ka^2)^2`$. $`\square`$

Ainsi, $`N`$ échantillons de Haar trouvent une violation de $`|\mathrm{Tr}(XUY)|\le\mathrm{Tr}(XY)`$ avec une
probabilité d'au plus $`Nb^4/(a^4k^2)`$ (une réunion de $`N`$ événements), bien qu'il en existe une pour tout couple qui
ne commute pas (proposition 3.9). La borne n'est informative que pour $`k>\sqrt N\,b^2/a^2`$ (pour $`b/a=6`$ et $`N=2000`$ :
$`k>1609`$).

*Calcul :* la proposition 11.11 par Monte-Carlo, et le rapport du corollaire 11.12 comparé au meilleur de $`2000`$ échantillons pour $`k\le32`$ — [tutoriel, section 11.4](../../tutorial/fr/11_optimisation.md#s-11-4).

<a id="r-11-13"></a>

**Lemme 11.13.** Si $`k\ge2`$, alors $`\min_U|\mathrm{Tr}(AU)|=0`$ pour toute $`A\in M_k(\mathbb C)`$.

*Démonstration.* Soient $`A=W\Sigma V^*`$ et $`P`$ le décalage cyclique, $`Pe_i=e_{(i\bmod k)+1}`$, dont la diagonale est
nulle pour $`k\ge2`$. Pour $`U=VPW^*`$, $`\mathrm{Tr}(AU)=\mathrm{Tr}(\Sigma P)=0`$. $`\square`$

Ainsi $`\inf_Uz=-\mathrm{Tr}(XY)`$ et $`\sup_Uz=\|YX\|_1-\mathrm{Tr}(XY)`$ (proposition 3.9). D'après la proposition 11.11, la moyenne quadratique de
$`|\mathrm{Tr}(XUY)|`$ vaut $`\|YX\|_2/\sqrt k`$ ; avec les bornes de la démonstration du corollaire 11.12 (spectres
dans $`[a,b]`$), $`\mathbb E|\mathrm{Tr}(XUY)|^2\le b^4`$ et $`\mathrm{Tr}(XY)\ge ka^2`$, donc l'inégalité de Markov
donne $`\mathbb P\big(z\ge-(1-\varepsilon)\mathrm{Tr}(XY)\big)\le b^4/(\varepsilon^2a^4k^2)`$ pour $`\varepsilon>0`$ :
pour $`k`$ grand, $`z`$ est proche de son infimum avec grande probabilité.

*Calcul :* les deux extrêmes sur $`U`$ par montée et descente de gradient, comparés à $`20\,000`$ échantillons de Haar — [tutoriel, section 11.4](../../tutorial/fr/11_optimisation.md#s-11-4).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Le corollaire 11.12 est un
> exemple faible, mais démontrable, de ce que la mécanique statistique appelle **typicité** : en grande dimension,
> presque tout état ou toute évolution est censé donner à peu près la même valeur d'une grandeur macroscopique. Une
> inégalité optimale se joue sur des configurations finement ajustées — bases propres alignées, spectres à deux
> niveaux — qu'un unitaire aléatoire (un circuit aléatoire, une matrice de diffusion aléatoire) produit rarement.

<a id="s-11-5"></a>

## 11.5 Recherche multi-départs des enveloppes

Une **recherche multi-départs** exécute une méthode de gradient, par tours, depuis de nombreux départs
indépendants, chacun avec ses propres emplacements, signes, $`U_0`$ et contraste ; elle recentre les cartes unitaires après
chaque tour et renvoie, pour chaque départ, la meilleure valeur de l'objectif aux itérés. Les signes du paramétrage
autoadjoint sont discrets ; entre les tours, ils sont modifiés par des changements de signe individuels en des
emplacements aléatoires, chacun étant conservé s'il améliore l'objectif.

*Dans la bibliothèque :* `optimize.extremize` — [tutoriel, section 11.5](../../tutorial/fr/11_optimisation.md#s-11-5).

En arithmétique exacte, chaque itéré est admissible (la proposition 11.2 et les remarques qui la suivent), donc chaque
valeur renvoyée est la valeur de $`z`$ en un point admissible ; en virgule flottante, c'est cette valeur aux erreurs d'arrondi près. Par conséquent, une maximisation donne une **borne
inférieure** du supremum, une minimisation une borne supérieure de l'infimum, et rien de plus — une recherche locale
peut s'arrêter en un extremum local. L'autre sens exige une démonstration, comme le corollaire 11.7 et la
proposition 11.9, ou une relaxation convexe (section 11.6). Une enveloppe calculée est confrontée (i) à des formes
closes ; (ii) à des configurations construites à la main, telles que des sommes directes de solutions plus petites,
en dessous desquelles elle ne doit pas tomber ; (iii) à des échantillons structurés (spectres à deux niveaux, bases
de permutation), dont les valeurs sont d'autres valeurs atteintes.

<a id="r-11-14"></a>

**Proposition 11.14 (enveloppe inférieure d'Exp. 1).** Supposons tous les $`k_c\ge2`$. Sur $`X\in\mathcal P_\delta`$,
$`Y\ge0`$ avec $`\lambda_{\max}(Y)=1`$ et les unitaires $`U\in\mathcal M`$, $`\ \min z=-\big(D-2\delta/(1+\delta)\big)`$.

*Démonstration.* $`Y\le1`$ donne $`X^{1/2}YX^{1/2}\le X`$, donc $`\mathrm{Tr}(XY)\le\mathrm{Tr}\,X\le(D-1)+\ell(\delta)`$,
puisqu'une valeur propre de $`X`$ vaut $`\ell(\delta)`$ et que les autres sont $`\le1`$ ; d'où
$`z\ge-\mathrm{Tr}(XY)\ge-(D-2\delta/(1+\delta))`$. L'égalité a lieu pour $`X`$ de valeurs propres $`\ell(\delta)`$ une
fois et $`1`$ sinon, $`Y=1`$, et $`U=\bigoplus_cU_c`$ avec $`\mathrm{Tr}(X_cU_c)=0`$ (lemme 11.13). $`\square`$

Les vérifications de type (ii) peuvent utiliser des sommes directes : le $`z`$ d'Exp. 3 est additif sur les sommes directes, et le contraste d'une somme directe est le plus grand
contraste de ses termes (exercice 3), donc, pour les enveloppes d'Exp. 3 avec $`\Delta(Y)`$ libre dans un intervalle fixé $`[0,b]`$,
```math
\overline z_{M_2\oplus M_2}(\delta)\ \ge\ \overline z_{M_2}(\delta)+\sup_{\delta'\le\delta}\overline z_{M_2}(\delta');
```
en particulier, la somme de deux valeurs atteintes sur $`M_2`$ aux contrastes $`\delta`$ et $`\delta'\le\delta`$ est une
valeur atteinte sur $`M_2\oplus M_2`$ au contraste $`\delta`$.

*Calcul :* la proposition 11.14 sur $`M_2\oplus M_2`$, le corollaire 11.7 sur $`M_3\oplus M_3`$ et la proposition 11.9 par recherche multi-départs, avec une vérification de type (ii) — [tutoriel, section 11.5](../../tutorial/fr/11_optimisation.md#s-11-5).

<a id="s-11-6"></a>

## 11.6 Points de comparaison

*Dans la bibliothèque :* la comparaison avec d'autres logiciels, parmi lesquels des solveurs SDP, dont les relaxations convexes majorent un supremum — [tutoriel, section 11.6](../../tutorial/fr/11_optimisation.md#s-11-6).

<a id="x-11"></a>

## Exercices

1. (**Transformation de Cayley.**) Pour $`S^*=-S`$, montrer que $`C(S)=(1-S)(1+S)^{-1}`$ est unitaire et que
   $`S\mapsto C(S)`$ est une bijection sur les unitaires sans valeur propre $`-1`$. Ainsi, contrairement à l'application
   exponentielle (lemme 11.1), une carte de Cayley $`U_0C(S)`$ manque $`-U_0`$.
2. (**Mesure de Haar réelle.**) Pour $`O`$ distribué selon la loi de Haar sur $`O(k)`$ (QR d'une matrice de Ginibre
   réelle) et $`A`$ réelle, montrer que $`\mathbb E\mathrm{Tr}(AO)^2=\|A\|_2^2/k`$, et le vérifier avec une algèbre
   réelle.
3. (**Sommes directes.**) Montrer que le $`z`$ d'Exp. 3 est additif sur les sommes directes et que
   $`\Delta(A\oplus A')=\max(\Delta(A),\Delta(A'))`$ pour $`A,A'\ge0`$ avec $`\lambda_{\max}=1`$. En déduire l'inégalité de
   somme directe de la section 11.5 (avec $`\Delta(Y)`$ libre dans un intervalle fixé $`[0,b]`$).

<details>
<summary>Solutions</summary>

1. $`S`$ est normal et $`iS`$ est autoadjoint, donc $`S`$ a des valeurs propres imaginaires pures (corollaire 1.4(2)) et
   $`1\pm S`$ sont inversibles ; ils commutent, et $`(1+S)^*=1-S`$. Par conséquent
   $`C^*C=(1-S)^{-1}(1+S)(1-S)(1+S)^{-1}=1`$, et $`1+C=2(1+S)^{-1}`$ est inversible, donc $`-1\notin\mathrm{spec}\,C`$.
   Réciproquement, pour un unitaire $`U`$ avec $`1+U`$ inversible, $`S=(1-U)(1+U)^{-1}`$ est antiautoadjoint :
   $`S^*=(1+U^*)^{-1}(1-U^*)=(U+1)^{-1}U\,U^*(U-1)=(U+1)^{-1}(U-1)=-S`$, car $`U\pm1`$ commutent. Pour $`C`$ et $`S`$ qui
   commutent, $`C(1+S)=1-S`$ équivaut à $`S(1+C)=1-C`$ (les deux disent $`C+CS+S=1`$). Or $`C(S)`$ commute avec $`S`$, et
   $`S=(1-U)(1+U)^{-1}`$ commute avec $`U`$ ; donc $`U=C(S)`$ si et seulement si $`S=(1-U)(1+U)^{-1}`$ : les deux applications
   sont réciproques l'une de l'autre.

2. D'après le théorème B.17 (voir aussi le lemme C.2), le facteur QR à signes corrigés $`O`$ d'une matrice gaussienne
   réelle, que la bibliothèque tire pour une algèbre réelle, est distribué selon la loi de Haar sur $`O(k)`$, donc $`VO`$
   a la même loi que $`O`$ pour $`V\in O(k)`$. La multiplication à gauche par $`\mathrm{diag}(\pm1)`$ donne
   $`\mathbb E[O_{ij}O_{lm}]=0`$ pour $`i\ne l`$, les permutations montrent que $`\mathbb E[O_{ij}O_{im}]`$ ne dépend pas de
   $`i`$, et $`\sum_iO_{ij}O_{im}=\delta_{jm}`$ ; donc $`\mathbb E[O_{ij}O_{lm}]=\delta_{il}\delta_{jm}/k`$ et
   $`\mathbb E\mathrm{Tr}(AO)^2=\sum_{ij}A_{ji}^2/k`$.

3. Les traces et les normes trace d'opérateurs diagonaux par blocs sont des sommes sur les blocs, et
   $`|X\oplus X'|=|X|\oplus|X'|`$ ; $`\lambda_{\min}(A\oplus A')=\min(\lambda_{\min}(A),\lambda_{\min}(A'))`$, et $`\Delta`$
   décroît avec $`\lambda_{\min}`$. Si donc $`(X,Y)`$ est admissible en $`\delta`$ et $`(X',Y')`$ en $`\delta'\le\delta`$, alors
   $`(X\oplus X',Y\oplus Y')`$ est admissible en $`\delta`$, avec $`z=z(X,Y)+z(X',Y')`$ ; on passe au supremum sur les deux
   couples et sur $`\delta'\le\delta`$.

</details>

*Calcul :* les exercices 1 et 2, numériquement — [tutoriel, Exercices](../../tutorial/fr/11_optimisation.md#x-11).

**Suite :** le [chapitre 12](12_equivariant_learning.md) traite de l'apprentissage automatique équivariant : les couches qui
commutent avec une symétrie sont des éléments d'un commutant (chapitre 2), construits à partir des mêmes opérateurs
différentiables et traités par lots.
