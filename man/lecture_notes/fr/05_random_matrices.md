[English](../05_random_matrices.md) · [Español](../es/05_random_matrices.md) · **Français** · [Deutsch](../de/05_random_matrices.md) · [中文](../zh/05_random_matrices.md) · [日本語](../ja/05_random_matrices.md) · [Русский](../ru/05_random_matrices.md)

# 5. Matrices aléatoires dans une algèbre

La théorie des matrices aléatoires cherche à quoi ressemble un opérateur « typique ». Dans une
algèbre $`\bigoplus_cM_{k_c}`$, la réponse doit respecter les secteurs : les unitaires aléatoires
naturels sont ceux de l'algèbre, et non ceux de la grande algèbre de matrices $`M_D`$, et les
hamiltoniens aléatoires naturels ont des blocs indépendants. Ce chapitre construit la mesure de
Haar sur le groupe unitaire de l'algèbre, calcule exactement ses moments d'ordre un et deux, puis
montre ce qui se passe mal dans les statistiques spectrales lorsqu'on ignore les secteurs.

**Vous apprendrez** pourquoi la décomposition QR d'une matrice gaussienne fournit un unitaire de
Haar (et pourquoi il faut fixer les phases), que la moyenne sur le groupe unitaire de l'algèbre est
l'espérance conditionnelle sur le centre, la formule de Weingarten d'ordre deux, la surmise de
Wigner, le rapport d'espacements et sa valeur exacte pour les matrices gaussiennes $`3\times3`$, et
pourquoi les statistiques de niveaux calculées sur un spectre où les secteurs sont mélangés
ressemblent à celles de niveaux indépendants, d'autant plus que les secteurs mélangés sont nombreux.

*Calcul :* les vérifications numériques de ce chapitre — [Tutoriel 5](../../tutorial/fr/05_random_matrices.md).

Dans tout ce chapitre, $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$ avec $`D=\sum_ck_c`$, les
$`P_c`$ sont les projecteurs centraux minimaux, et $`E`$ est l'espérance conditionnelle sur le centre
de la proposition 2.8. $`U(k)`$ désigne le groupe des matrices unitaires $`k\times k`$ et
$`U(\mathcal M)=\{U\in\mathcal M:U^*U=1\}`$ le groupe unitaire de $`\mathcal M`$. Les variables
aléatoires sont définies sur un espace probabilisé d'espérance $`\mathbb E`$ ; l'espérance d'une
matrice aléatoire se prend coefficient par coefficient. Une variable *gaussienne complexe
standard* a des parties réelle et imaginaire indépendantes de loi $`\sim N(0,\tfrac12)`$, de sorte
que $`\mathbb E|z|^2=1`$ ; une *matrice de Ginibre* a des coefficients indépendants gaussiens
complexes standard.

<a id="s-5-1"></a>

## 5.1 Mesure de Haar

Soit $`G`$ un groupe compact de matrices (un sous-groupe fermé d'un certain $`U(n)`$). Une mesure de
probabilité borélienne $`\mu`$ sur $`G`$ est **invariante à gauche** si $`\mu(gB)=\mu(B)`$ pour tout
$`g\in G`$ et tout borélien $`B`$, ou, de manière équivalente, si $`\int f(gh)\,d\mu(h)=\int f\,d\mu`$
pour toute fonction borélienne bornée $`f`$ ; l'**invariance à droite** se définit avec $`Bg`$ et
$`f(hg)`$.

<a id="r-5-1"></a>

**Lemme 5.1 (unicité).** Si $`\mu`$ est une mesure de probabilité borélienne invariante à gauche et
$`\nu`$ une mesure de probabilité borélienne invariante à droite sur $`G`$, alors $`\mu=\nu`$. Par
conséquent, $`G`$ porte au plus une mesure de probabilité borélienne invariante à gauche, et si elle
existe, elle est aussi invariante à droite et invariante par $`g\mapsto g^{-1}`$.

*Démonstration.* Soit $`f:G\to\mathbb R`$ continue (donc bornée). D'après le théorème de Fubini,
$`\int\!\big(\int f(gh)\,d\mu(h)\big)d\nu(g)=\int\!\big(\int f(gh)\,d\nu(g)\big)d\mu(h)`$.
L'intégrale intérieure de gauche vaut $`\int f\,d\mu`$ par invariance à gauche, celle de droite vaut
$`\int f\,d\nu`$ par invariance à droite ; donc $`\int f\,d\mu=\int f\,d\nu`$. Deux mesures de
probabilité boréliennes sur un espace métrique qui ont les mêmes intégrales sur les fonctions
continues sont égales (corollaire A.2 : les fonctions continues $`\max(0,1-n\mathrm{dist}(\cdot,K))`$
décroissent vers l'indicatrice d'un fermé $`K`$, et les fermés déterminent la mesure d'après le
lemme de Dynkin, lemme A.1). Si $`\mu_1,\mu_2`$ sont invariantes à gauche, alors
$`\nu(B)=\mu_2(B^{-1})`$ est invariante à droite, puisque $`(Bg)^{-1}=g^{-1}B^{-1}`$ ; donc
$`\mu_1=\nu=\mu_2`$, et $`\mu_1`$ est invariante à droite et par inversion. $`\square`$

**Définition.** La **mesure de Haar** de $`G`$ est sa mesure de probabilité borélienne invariante à
gauche. Un élément aléatoire de cette loi est dit **distribué selon la loi de Haar**.

Le lemme 5.1 ne dit rien de l'existence. Le cours n'a besoin que des mesures de Haar de $`U(k)`$,
$`U(1)`$, $`U(\mathcal M)`$, $`SU(k)`$, $`O(k)`$ et des groupes diagonaux $`T`$ du lemme C.2, et chacune
d'elles est construite explicitement : $`U(k)`$ au théorème 5.3, $`U(\mathcal M)`$ au corollaire 5.4,
$`SU(k)`$ juste après, $`O(k)`$ au théorème B.17 et $`T`$ au lemme C.2 (voir la section B.7).

<a id="r-5-2"></a>

**Lemme 5.2 (décomposition QR).** Toute $`Z\in M_k(\mathbb C)`$ inversible s'écrit $`Z=QR`$ avec $`Q`$
unitaire et $`R`$ triangulaire supérieure à diagonale positive. $`Q`$ et $`R`$ sont uniques et dépendent
continûment de $`Z`$.

*Démonstration.* Appliquons le procédé de Gram–Schmidt aux colonnes $`z_1,\dots,z_k`$ de $`Z`$, qui
sont linéairement indépendantes : $`q_j=w_j/\|w_j\|`$ avec
$`w_j=z_j-\sum_{i<j}\langle q_i,z_j\rangle q_i`$ ; ici $`w_j\neq0`$ car
$`z_j\notin\mathrm{span}\{z_1,\dots,z_{j-1}\}=\mathrm{span}\{q_1,\dots,q_{j-1}\}`$. Alors
$`z_j=\sum_{i<j}\langle q_i,z_j\rangle q_i+\|w_j\|q_j`$, c'est-à-dire $`Z=QR`$ avec $`R_{jj}=\|w_j\|>0`$,
et chaque étape est continue sur les matrices inversibles. Si $`Q_1R_1=Q_2R_2`$, alors
$`T=Q_2^*Q_1=R_2R_1^{-1}`$ est unitaire et triangulaire supérieure à diagonale positive. Son inverse
$`T^*`$ est triangulaire inférieure, et aussi triangulaire supérieure (l'inverse d'une matrice
triangulaire supérieure est triangulaire supérieure), donc $`T`$ est diagonale ; une matrice
unitaire diagonale à diagonale positive vaut $`1`$. $`\square`$

<a id="r-5-3"></a>

**Théorème 5.3 (unitaires de Haar à partir de matrices gaussiennes).** Soit $`Z`$ une matrice de
Ginibre $`k\times k`$. Alors $`Z`$ est inversible avec probabilité un, et le facteur $`Q`$ du lemme 5.2
est distribué selon la loi de Haar sur $`U(k)`$. De plus, si $`Z=Q'R'`$ est une décomposition QR
quelconque avec $`Q'`$ unitaire et $`R'`$ triangulaire supérieure (par exemple celle que renvoie
`torch.linalg.qr`), alors $`Q=Q'\Lambda`$ avec $`\Lambda=\mathrm{diag}(R'_{ii}/|R'_{ii}|)`$ est
distribué selon la loi de Haar.

*Démonstration.* *Inversibilité.* $`|\det Z|^2`$ est un polynôme réel en les $`2k^2`$ coordonnées
réelles de $`Z`$, qui n'est pas identiquement nul (il vaut $`1`$ en $`Z=1`$), et $`Z`$ a la densité
$`\pi^{-k^2}e^{-\|Z\|_2^2}`$ (proposition A.16(b)) ; donc $`P(\det Z=0)=0`$ d'après le
corollaire A.7(b).

*Invariance de la loi gaussienne.* Pour $`V\in U(k)`$, l'application $`Z\mapsto VZ`$ est linéaire sur
les réels et conserve $`\|Z\|_2`$ (section 1.4) ; c'est donc une application orthogonale de
$`M_k(\mathbb C)\cong\mathbb R^{2k^2}`$, qui conserve la mesure de Lebesgue (lemme A.3(b)) et la
densité : $`VZ`$ a même loi que $`Z`$ (proposition A.16(b)).

*Équivariance.* Si $`Z=QR`$, alors $`VZ=(VQ)R`$ est la décomposition du lemme 5.2 de $`VZ`$, donc
$`Q(VZ)=VQ(Z)`$. Comme $`Q`$ est continue sur les matrices inversibles (lemme 5.2), $`Q(Z)`$ est un
élément aléatoire de $`U(k)`$ et, pour tout borélien $`B\subseteq U(k)`$,
$`P(VQ(Z)\in B)=P(Q(VZ)\in B)=P(Q(Z)\in B)`$ : la loi de $`Q(Z)`$ est invariante à gauche et, d'après le
lemme 5.1, c'est la mesure de Haar (ce qui prouve aussi que la mesure de Haar de $`U(k)`$ existe).
Enfin $`Q'R'=(Q'\Lambda)(\Lambda^*R')`$, et $`\Lambda^*R'`$ est triangulaire supérieure de diagonale
$`|R'_{ii}|>0`$, donc $`Q'\Lambda=Q(Z)`$ par unicité. $`\square`$

La même démonstration avec des matrices gaussiennes réelles donne la mesure de Haar du groupe
orthogonal $`O(k)`$, avec la correction $`\mathrm{diag}(\mathrm{sgn}\,R'_{ii})`$ ; c'est ce
que détaille le théorème B.17. Sans la correction $`\Lambda`$, la loi de $`Q'`$ dépend de la convention
de signe de la routine QR et n'est pas nécessairement la mesure de Haar (la recette avec $`\Lambda`$
est due à F. Mezzadri, *How to generate random matrices from the classical compact groups*,
Notices AMS 54 (2007), 592–604). Par exemple, si la routine rend $`R'_{11}`$ réel de signe opposé à celui de
$`\mathrm{Re}\,Z_{11}`$ (comme le fait la QR de Householder de LAPACK), les premières colonnes $`z_1=R'_{11}q'_1`$ de $`Z=Q'R'`$ donnent
$`\mathrm{Re}\,Q'_{11}=-|\mathrm{Re}\,Z_{11}|/\|z_1\|<0`$ presque sûrement ($`\mathrm{Re}\,Z_{11}\sim N(0,\frac12)`$ s'annule avec probabilité $`0`$) ;
donc $`\mathbb E\,\mathrm{Re}\,Q'_{11}<0`$, alors que $`\mathbb EU_{11}=0`$ pour un unitaire de Haar (corollaire 5.5).

*Calcul :* l'échantillonneur de la bibliothèque `TypeIAlgebra.random_unitary` est la QR corrigée du théorème 5.3, et le $`Q'`$ non corrigé de `torch.linalg.qr` n'est pas distribué selon la loi de Haar — [tutoriel, section 5.1](../../tutorial/fr/05_random_matrices.md#s-5-1).

<a id="r-5-4"></a>

**Corollaire 5.4 (mesure de Haar de l'algèbre).** $`U(\mathcal M)=\{\bigoplus_cU_c:U_c\in U(k_c)\}\cong\prod_cU(k_c)`$,
et $`U=\bigoplus_cU_c`$ est distribué selon la loi de Haar sur $`U(\mathcal M)`$ si et seulement si les
blocs $`U_c`$ sont indépendants et chaque $`U_c`$ est distribué selon la loi de Haar sur $`U(k_c)`$.

*Démonstration.* $`U^*U=\bigoplus_cU_c^*U_c`$, donc $`U`$ est unitaire si et seulement si chaque bloc
l'est. Supposons les $`U_c`$ indépendants et distribués selon la loi de Haar. Pour
$`V=\bigoplus_cV_c\in U(\mathcal M)`$, les blocs $`V_cU_c`$ de $`VU`$ sont indépendants (fonctions de
variables indépendantes) et $`V_cU_c`$ a la loi de $`U_c`$, donc $`VU`$ a la loi de $`U`$. La loi de $`U`$
est donc invariante à gauche et, d'après le lemme 5.1, c'est la mesure de Haar du groupe compact
$`U(\mathcal M)`$. Réciproquement, par unicité, tout $`U`$ distribué selon la loi de Haar a cette
loi. $`\square`$

*Dans la bibliothèque :* les unitaires de Haar de $`U(\mathcal M)`$, à blocs indépendants, sont tirés par `TypeIAlgebra.random_unitary_operator` — [tutoriel, section 5.1](../../tutorial/fr/05_random_matrices.md#s-5-1).

Le groupe $`U(\mathcal M)`$ est bien plus petit que $`U(D)`$ : pour $`C\ge2`$, un unitaire de Haar $`U`$
de $`U(D)`$ appartient à $`\mathcal M`$ avec probabilité nulle. En effet, d'après le théorème 5.3 et le
lemme 5.1, $`U`$ a la loi du facteur $`Q`$ d'une matrice de Ginibre $`Z`$ de taille $`D\times D`$, et
d'après la construction de Gram–Schmidt du lemme 5.2, la première colonne de $`Q`$ est $`z/\|z\|`$ pour
la première colonne $`z`$ de $`Z`$, un vecteur gaussien complexe standard de $`\mathbb C^D`$ qui admet une
densité (proposition A.16(a)). Pour chaque $`i`$, $`\{z_i=0\}`$ est un sous-espace propre de
$`\mathbb C^D`$, donc $`P(z_i=0)=0`$ (corollaire A.7(c)), et tous les coefficients de la première
colonne sont presque sûrement non nuls. Un élément de $`\mathcal M`$ a des coefficients nuls dans la
première colonne en dehors des lignes du premier bloc, et ces lignes existent lorsque $`C\ge2`$.

La mesure de Haar de $`SU(k)`$ s'obtient dans le même esprit : si $`U`$ suit la loi de Haar sur $`U(k)`$
et $`\theta=\arg\det U\in(-\pi,\pi]`$ (une fonction borélienne de $`U`$), alors $`V=e^{-i\theta/k}U\in SU(k)`$,
et pour $`W\in SU(k)`$, $`\det(WU)=\det U`$ donne $`V(WU)=WV(U)`$ ; comme $`WU`$ a la loi de $`U`$, la loi de
$`V`$ est invariante à gauche sur $`SU(k)`$, et c'est donc la mesure de Haar de $`SU(k)`$ d'après le
lemme 5.1.

<a id="r-5-5"></a>

**Corollaire 5.5 (premiers moments des coefficients).** Si $`U`$ est distribué selon la loi de Haar
sur $`U(k)`$, alors $`\mathbb EU_{ij}=0`$ et $`\mathbb E|U_{ij}|^2=1/k`$.

*Démonstration.* $`-1\in U(k)`$, donc $`-U`$ a la loi de $`U`$ et $`\mathbb EU=-\mathbb EU`$. Pour une
matrice de permutation $`P`$, $`PU`$ a la loi de $`U`$, donc $`\mathbb E|U_{ij}|^2`$ ne dépend pas de $`i`$ ;
comme $`\sum_i|U_{ij}|^2=1`$, cette quantité vaut $`1/k`$. $`\square`$

*Calcul :* les corollaires 5.4 et 5.5 par Monte-Carlo, et la mesure de Haar de $`SU(k)`$ (`measure='haar_su'`) — [tutoriel, section 5.1](../../tutorial/fr/05_random_matrices.md#s-5-1).

<a id="s-5-2"></a>

## 5.2 Moments et calcul de Weingarten

<a id="r-5-6"></a>

**Lemme 5.6.** Si $`T\in M_k(\mathbb C)`$ commute avec tous les unitaires diagonaux et toutes les
matrices de permutation, alors $`T\in\mathbb C1`$. En particulier, le commutant de $`U(k)`$ est
$`\mathbb C1`$.

*Démonstration.* Pour $`\Delta=\mathrm{diag}(e^{i\theta_1},\dots,e^{i\theta_k})`$ avec des
$`e^{i\theta_j}`$ deux à deux distincts, $`\Delta T=T\Delta`$ donne
$`(e^{i\theta_i}-e^{i\theta_j})T_{ij}=0`$, donc $`T`$ est diagonale. Pour la matrice de permutation
$`P_\sigma e_j=e_{\sigma(j)}`$, $`P_\sigma TP_\sigma^{-1}=T`$ donne $`T_{\sigma(j)\sigma(j)}=T_{jj}`$, donc
tous les coefficients diagonaux sont égaux. $`\square`$

<a id="r-5-7"></a>

**Proposition 5.7 (le twirl de Haar est l'espérance conditionnelle sur le centre).** Soit $`U`$
distribué selon la loi de Haar sur $`U(\mathcal M)`$. Pour tout $`X\in\mathcal M`$,
```math
\mathbb E\,[UXU^*]=E(X)=\bigoplus_c\frac{\mathrm{Tr}\,X_c}{k_c}\,1_{k_c}.
```

*Démonstration.* $`(UXU^*)_c=U_cX_cU_c^*`$, et $`U_c`$ suit la loi de Haar sur $`U(k_c)`$
(corollaire 5.4). Posons $`T_c=\mathbb E[U_cX_cU_c^*]`$ (les coefficients sont bornés). Pour
$`V\in U(k_c)`$, $`VU_c`$ a la loi de $`U_c`$, donc $`VT_cV^*=\mathbb E[(VU_c)X_c(VU_c)^*]=T_c`$. D'après le
lemme 5.6, $`T_c=\lambda1`$, et
$`k_c\lambda=\mathrm{Tr}\,T_c=\mathbb E\mathrm{Tr}(U_cX_cU_c^*)=\mathrm{Tr}\,X_c`$. $`\square`$

La normalisation est celle de la proposition 2.8 et de `channels.center_expectation`. Si
$`X\in\mathcal M`$ commute avec $`U(\mathcal M)`$, alors $`UXU^*=X`$, donc $`X=E(X)`$ est central ; la
moyenne de Haar est l'analogue continu de la moyenne sur le groupe du chapitre 2, exercice 3. Si au
contraire $`U`$ suit la loi de Haar sur $`U(D)`$, la même démonstration (avec $`C=1`$) donne
$`\mathbb E[UXU^*]=\frac{\mathrm{Tr}\,X}{D}1`$ : les poids des secteurs sont perdus, et $`UXU^*`$
n'appartient même pas nécessairement à $`\mathcal M`$ (pour $`C\ge2`$, $`X=P_1`$ et la rotation $`U`$
d'angle $`\pi/4`$ dans le plan de $`e_1`$ et $`e_{k_1+1}`$, $`UP_1U^*`$ a le coefficient $`\frac12`$ en
position $`(1,k_1+1)`$).

*Calcul :* la proposition 5.7 par Monte-Carlo, comparée à `channels.center_expectation` — [tutoriel, section 5.2](../../tutorial/fr/05_random_matrices.md#s-5-2).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Un unitaire aléatoire de Haar de l'algèbre est une évolution
> aléatoire qui respecte la règle de supersélection. En moyennant sur lui, on brouille tout à
> l'intérieur de chaque secteur et on conserve exactement la partie classique : le poids de chaque
> charge. C'est le « twirl » de l'information quantique, et c'est l'espérance conditionnelle sur le
> centre du chapitre 2.

Pour les moments d'ordre deux, il nous faut le commutant de $`\{V\otimes V\}`$ sur
$`\mathbb C^k\otimes\mathbb C^k`$. Soit $`F`$ l'**échange** (*swap*), $`F(x\otimes y)=y\otimes x`$ ; il
est autoadjoint et unitaire, et
$`\mathrm{Tr}\,1=k^2`$, $`\mathrm{Tr}\,F=\sum_{i,j}\langle e_i\otimes e_j,e_j\otimes e_i\rangle=k`$, $`F^2=1`$.

<a id="r-5-8"></a>

**Lemme 5.8.** Soit $`k\ge2`$. Une matrice $`X`$ sur $`\mathbb C^k\otimes\mathbb C^k`$ commute avec
$`V\otimes V`$ pour tout $`V\in U(k)`$ si et seulement si $`X=a1+bF`$ avec $`a,b\in\mathbb C`$ ; de plus,
$`1`$ et $`F`$ sont linéairement indépendants.

*Démonstration.* $`F(V\otimes V)(x\otimes y)=Vy\otimes Vx=(V\otimes V)F(x\otimes y)`$, donc $`1`$ et $`F`$
commutent avec tous les $`V\otimes V`$. Réciproquement, supposons que $`X`$ commute avec eux, et
écrivons $`X_{ij,pq}=\langle e_i\otimes e_j,X(e_p\otimes e_q)\rangle`$.
*Unitaires diagonaux :* $`\Delta\otimes\Delta`$ multiplie $`e_i\otimes e_j`$ par
$`e^{i(\theta_i+\theta_j)}`$, donc $`(e^{i(\theta_i+\theta_j)}-e^{i(\theta_p+\theta_q)})X_{ij,pq}=0`$. Si
les multiensembles $`\{i,j\}`$ et $`\{p,q\}`$ diffèrent, un indice $`m`$ apparaît $`n\in\{\pm1,\pm2\}`$
fois de plus dans le premier ; avec $`\theta_m=\pi/2`$ et tous les autres $`\theta=0`$, le crochet
vaut $`e^{i(\theta_p+\theta_q)}(e^{in\pi/2}-1)\neq0`$, donc $`X_{ij,pq}=0`$. Par conséquent
$`X(e_i\otimes e_i)=a_ie_i\otimes e_i`$ et, pour $`i\neq j`$,
$`X(e_i\otimes e_j)=b_{ij}\,e_i\otimes e_j+c_{ij}\,e_j\otimes e_i`$.
*Permutations :* en appliquant $`X(P_\sigma\otimes P_\sigma)=(P_\sigma\otimes P_\sigma)X`$ à
$`e_i\otimes e_j`$, on obtient $`a_{\sigma(i)}=a_i`$, $`b_{\sigma(i)\sigma(j)}=b_{ij}`$,
$`c_{\sigma(i)\sigma(j)}=c_{ij}`$. Les permutations agissent transitivement sur les couples d'indices
distincts, donc $`a_i=a'`$, $`b_{ij}=b'`$, $`c_{ij}=c'`$ ; ainsi $`X=b'1+c'F`$ sur le sous-espace engendré
par les $`e_i\otimes e_j`$ avec $`i\ne j`$.
*Une rotation :* soit $`R`$ l'unitaire tel que $`Re_1=(e_1+e_2)/\sqrt2`$, $`Re_2=(-e_1+e_2)/\sqrt2`$, $`Re_m=e_m`$ pour $`m\ge3`$.
Alors $`(R\otimes R)(e_1\otimes e_1)=\frac12(e_1\otimes e_1+e_1\otimes e_2+e_2\otimes e_1+e_2\otimes e_2)`$, et
$`X(R\otimes R)(e_1\otimes e_1)=\frac12\big(a'e_1\otimes e_1+(b'+c')(e_1\otimes e_2+e_2\otimes e_1)+a'e_2\otimes e_2\big)`$
doit être égal à $`(R\otimes R)X(e_1\otimes e_1)=a'(R\otimes R)(e_1\otimes e_1)`$. En comparant les
coefficients de $`e_1\otimes e_2`$, on obtient $`a'=b'+c'`$, qui est la valeur de $`b'1+c'F`$ sur
$`e_i\otimes e_i`$. Donc $`X=b'1+c'F`$. Enfin, $`1`$ et $`F`$ envoient $`e_1\otimes e_2`$ sur les vecteurs
indépendants $`e_1\otimes e_2`$ et $`e_2\otimes e_1`$. $`\square`$

*Calcul :* le lemme 5.8 par algèbre linéaire, avec le commutant de trois $`V\otimes V`$ tirés au hasard — [tutoriel, section 5.2](../../tutorial/fr/05_random_matrices.md#s-5-2).

<a id="r-5-9"></a>

**Corollaire 5.9 (dualité de Schur–Weyl d'ordre deux).** Soient $`k\ge2`$ et $`P_\pm=\frac12(1\pm F)`$.
Ce sont les projecteurs sur les sous-espaces symétrique et antisymétrique
$`\mathrm{Sym}=\{x:Fx=x\}`$ et $`\mathrm{Alt}=\{x:Fx=-x\}`$, de dimensions $`k(k+1)/2`$ et $`k(k-1)/2`$.
La représentation $`V\mapsto V\otimes V`$ de $`U(k)`$ est irréductible sur $`\mathrm{Sym}`$ et sur
$`\mathrm{Alt}`$, et les deux représentations ne sont pas équivalentes.

*Démonstration.* $`F=F^*`$ et $`F^2=1`$ font de $`P_\pm`$ des projecteurs avec $`P_++P_-=1`$ et d'images
$`\mathrm{Sym}`$, $`\mathrm{Alt}`$ ; leurs dimensions sont $`\mathrm{Tr}\,P_\pm=(k^2\pm k)/2`$. Le
sous-espace $`\mathcal A`$ engendré par $`\{V\otimes V\}`$ est une $`*`$-algèbre ($`1\otimes1`$, produits
$`VW\otimes VW`$, adjoints $`V^*\otimes V^*`$) de commutant
$`\mathcal A'=\mathrm{span}\{1,F\}=\mathrm{span}\{P_+,P_-\}`$ d'après le lemme 5.8. Soit
$`K\subseteq\mathrm{Sym}`$ invariant par tous les $`V\otimes V`$. D'après le lemme 2.1, son projecteur
est dans $`\mathcal A'`$, il s'écrit donc $`xP_++yP_-`$ avec $`x,y\in\{0,1\}`$ (c'est un projecteur) ; il
s'annule sur $`\mathrm{Alt}`$, donc $`y=0`$, et $`K=0`$ ou $`K=\mathrm{Sym}`$. Le même argument s'applique
à $`\mathrm{Alt}`$. Si $`T:\mathrm{Sym}\to\mathrm{Alt}`$ entrelace les deux représentations, alors
$`TP_+`$ (vu comme opérateur sur $`\mathbb C^k\otimes\mathbb C^k`$) commute avec $`\mathcal A`$, donc
$`TP_+=xP_++yP_-`$ ; il s'annule sur $`\mathrm{Alt}`$, donc $`y=0`$, et son image est contenue dans
$`\mathrm{Alt}\cap\mathrm{Sym}=0`$, donc $`T=0`$. $`\square`$

<a id="r-5-10"></a>

**Théorème 5.10 (formule de Weingarten d'ordre deux).** Soient $`k\ge2`$, $`U`$ distribué selon la loi
de Haar sur $`U(k)`$, et $`X`$ une matrice sur $`\mathbb C^k\otimes\mathbb C^k`$. Alors
```math
\mathbb E\,[(U\otimes U)X(U\otimes U)^*]=\alpha\,1+\beta F,\qquad
\alpha=\frac{\mathrm{Tr}\,X-\mathrm{Tr}(FX)/k}{k^2-1},\quad
\beta=\frac{\mathrm{Tr}(FX)-\mathrm{Tr}\,X/k}{k^2-1}.
```
De manière équivalente, avec $`\mathrm{Wg}(e)=\frac1{k^2-1}`$, $`\mathrm{Wg}(\tau)=-\frac1{k(k^2-1)}`$
sur le groupe $`S_2=\{e,\tau\}`$,
```math
\mathbb E\big[U_{i_1j_1}U_{i_2j_2}\overline{U_{i'_1j'_1}U_{i'_2j'_2}}\big]
=\sum_{\sigma,\pi\in S_2}\delta_{i_1i'_{\sigma(1)}}\delta_{i_2i'_{\sigma(2)}}\,
\delta_{j_1j'_{\pi(1)}}\delta_{j_2j'_{\pi(2)}}\,\mathrm{Wg}(\sigma\pi^{-1}).
```

*Démonstration.* Posons $`T=\mathbb E[(U\otimes U)X(U\otimes U)^*]`$. Comme dans la proposition 5.7,
$`T`$ commute avec tous les $`V\otimes V`$, donc $`T=\alpha1+\beta F`$ d'après le lemme 5.8. En prenant
les traces, $`\mathrm{Tr}\,T=\mathrm{Tr}\,X`$ et, puisque $`F`$ commute avec $`U\otimes U`$,
$`\mathrm{Tr}(FT)=\mathbb E\mathrm{Tr}((U\otimes U)^*F(U\otimes U)X)=\mathrm{Tr}(FX)`$.
On obtient $`k^2\alpha+k\beta=\mathrm{Tr}\,X`$ et $`k\alpha+k^2\beta=\mathrm{Tr}(FX)`$, un
système linéaire de déterminant $`k^2(k^2-1)\neq0`$, dont la solution est celle de l'énoncé. Pour
les coefficients, prenons $`X=|e_{j_1}\otimes e_{j_2}\rangle\langle e_{j'_1}\otimes e_{j'_2}|`$ :
alors $`\langle e_{i_1}\otimes e_{i_2},T\,e_{i'_1}\otimes e_{i'_2}\rangle`$ est le membre de gauche,
$`\mathrm{Tr}\,X=\delta_{j_1j'_1}\delta_{j_2j'_2}`$ (le terme $`\pi=e`$),
$`\mathrm{Tr}(FX)=\delta_{j_1j'_2}\delta_{j_2j'_1}`$ (le terme $`\pi=\tau`$), et $`1`$, $`F`$
contribuent $`\delta_{i_1i'_1}\delta_{i_2i'_2}`$ (le terme $`\sigma=e`$) et
$`\delta_{i_1i'_2}\delta_{i_2i'_1}`$ ($`\sigma=\tau`$). En regroupant les termes,
$`\alpha=\sum_\pi\mathrm{Wg}(e\pi^{-1})[\pi]`$ et $`\beta=\sum_\pi\mathrm{Wg}(\tau\pi^{-1})[\pi]`$,
où $`[\pi]`$ désigne les deltas en $`j`$ associés à $`\pi`$. $`\square`$

La matrice $`\big(\mathrm{Wg}(\sigma^{-1}\pi)\big)_{\sigma,\pi}`$ est l'inverse de la matrice de Gram
$`\begin{pmatrix}k^2&k\\k&k^2\end{pmatrix}`$ de $`1`$ et $`F`$ pour
$`\langle A,B\rangle=\mathrm{Tr}(A^*B)`$ ; c'est sous cette forme que le résultat se
généralise.

*Calcul :* le théorème 5.10 par Monte-Carlo pour $`k=3`$ — [tutoriel, section 5.2](../../tutorial/fr/05_random_matrices.md#s-5-2).

<a id="r-5-11"></a>

**Corollaire 5.11.** Pour $`k\ge2`$ :
```math
\mathbb E|U_{11}|^4=\frac2{k(k+1)},\qquad
\mathbb E|U_{11}|^2|U_{12}|^2=\frac1{k(k+1)},\qquad
\mathbb E\big[U_{11}U_{22}\overline{U_{12}U_{21}}\big]=-\frac1{(k-1)k(k+1)}.
```

*Démonstration.* Dans le théorème 5.10, (i) tous les indices valent $`1`$ : chaque produit de deltas
vaut $`1`$, donc le moment vaut $`2(\mathrm{Wg}(e)+\mathrm{Wg}(\tau))=2\frac{k-1}{k(k^2-1)}`$.
(ii) $`i=(1,1)`$, $`j=(1,2)`$, $`i'=(1,1)`$, $`j'=(1,2)`$ : seul $`\pi=e`$ subsiste, ce qui donne
$`\mathrm{Wg}(e)+\mathrm{Wg}(\tau)`$. (iii) $`i=(1,2)`$, $`j=(1,2)`$, $`i'=(1,2)`$, $`j'=(2,1)`$ : seuls
$`\sigma=e`$, $`\pi=\tau`$ subsistent, ce qui donne $`\mathrm{Wg}(\tau)`$. $`\square`$

*Calcul :* les trois moments du corollaire 5.11 par Monte-Carlo pour $`k=3`$ — [tutoriel, section 5.2](../../tutorial/fr/05_random_matrices.md#s-5-2).

Le dernier moment est négatif : les coefficients d'un unitaire de Haar ne sont pas indépendants,
et les corrélations sont d'ordre relatif $`1/k`$. Pour l'algèbre, les coefficients de blocs
différents sont indépendants (corollaire 5.4), de sorte que les moments mixtes entre secteurs se
factorisent en moments de blocs individuels.

Pour les moments d'ordre quelconque, le même raisonnement s'applique dès que l'on connaît le
commutant de $`\{V^{\otimes n}\}`$ ; il est engendré par les permutations des $`n`$ facteurs tensoriels
(dualité de Schur–Weyl, théorème B.6(c)). L'annexe B mène ce programme à bien.

<a id="r-5-12"></a>

**Théorème 5.12 (formule de Weingarten d'ordre quelconque).** Soient $`n\ge1`$ et $`U`$ distribué selon
la loi de Haar sur $`U(k)`$. Alors
```math
\mathbb E\big[U_{i_1j_1}\cdots U_{i_nj_n}\overline{U_{i'_1j'_1}\cdots U_{i'_nj'_n}}\big]
=\sum_{\sigma,\pi\in S_n}\prod_{m=1}^n\delta_{i_mi'_{\sigma(m)}}\delta_{j_mj'_{\pi(m)}}\;\mathrm{Wg}(\sigma\pi^{-1},k),
```
où $`\mathrm{Wg}(\cdot,k)`$ est une fonction sur $`S_n`$ qui ne dépend que du type cyclique. Pour
$`n\le k`$, la matrice $`\big(\mathrm{Wg}(\sigma^{-1}\pi,k)\big)_{\sigma,\pi}`$ est l'inverse de la
matrice de Gram $`\big(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)}\big)_{\sigma,\pi}`$ ; pour $`n>k`$, cette
matrice est singulière et l'on prend son pseudo-inverse. Les moments où les nombres de facteurs
$`U`$ et $`\overline U`$ diffèrent sont nuls.

*Démonstration.* C'est le théorème B.13 avec $`d=k`$ : la partie (a) donne
$`W_{\sigma\pi}=\mathrm{Wg}(\sigma^{-1}\pi)`$ avec $`\mathrm{Wg}`$ fonction du type cyclique, la
partie (b) l'inversibilité de la matrice de Gram si et seulement si $`n\le k`$, la partie (d) la
formule des moments et la partie (e) la nullité des moments non équilibrés. $`\square`$

(La formule est due à D. Weingarten, J. Math. Phys. 19 (1978), 999–1001, et à B. Collins, Int.
Math. Res. Not. 2003, 953–982 ; voir aussi B. Collins, P. Śniady, Comm. Math. Phys. 264 (2006),
773–795.)

*Calcul :* la matrice de Weingarten d'ordre $`3`$ pour $`k=3`$ et $`\mathbb E|U_{11}|^6=\frac1{10}`$ (exercice 2) — [tutoriel, Exercices](../../tutorial/fr/05_random_matrices.md#x-5).

<a id="s-5-3"></a>

## 5.3 Statistiques de niveaux

Soit $`H`$ une matrice autoadjointe aléatoire $`n\times n`$ de valeurs propres
$`\lambda_1\le\dots\le\lambda_n`$. Ses **espacements** sont les $`s_i=\lambda_{i+1}-\lambda_i`$. Deux
ensembles sont classiques : le **GOE**, $`H=(G+G^T)/2`$ avec $`G`$ réelle à coefficients indépendants
de loi $`N(0,1)`$, et le **GUE**, $`H=(G+G^*)/\sqrt2`$ avec $`G`$ de Ginibre (section 1.7). Dans les deux
cas, les coefficients diagonaux suivent la loi $`N(0,1)`$, les coefficients au-dessus de la
diagonale suivent la loi $`N(0,\tfrac12)`$ dans le GOE et ont des parties réelle et imaginaire
indépendantes de loi $`N(0,\tfrac12)`$ dans le GUE, les coefficients $`h_{ij}`$, $`i\le j`$, sont
indépendants, et la densité de $`H`$ est proportionnelle à $`e^{-\mathrm{Tr}\,H^2/2}`$
(proposition A.17 ; la proportionnalité provient de $`\mathrm{Tr}\,H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2`$).

<a id="r-5-13"></a>

**Proposition 5.13 (surmise de Wigner).** Soient $`n=2`$, $`s=\lambda_2-\lambda_1`$ et
$`S=s/\mathbb Es`$. Pour le GOE ($`\beta=1`$) et le GUE ($`\beta=2`$), $`S`$ a pour densités
```math
p_1(S)=\frac\pi2\,S\,e^{-\pi S^2/4},\qquad p_2(S)=\frac{32}{\pi^2}\,S^2\,e^{-4S^2/\pi},
```
de fonctions de répartition $`P(S\le t)=1-e^{-\pi t^2/4}`$ et
$`P(S\le t)=\mathrm{erf}(2t/\sqrt\pi)-\frac{4t}\pi e^{-4t^2/\pi}`$. En particulier,
$`P(S\le\varepsilon)`$ est d'ordre $`\varepsilon^{\beta+1}`$ quand $`\varepsilon\to0`$.

*Démonstration.* Pour $`H=\begin{pmatrix}a&b\\\bar b&c\end{pmatrix}`$, les valeurs propres sont
$`\frac{a+c}2\pm\frac12\sqrt{(a-c)^2+4|b|^2}`$, donc $`s=\sqrt{(a-c)^2+4|b|^2}`$. D'après la
proposition A.17, $`a`$, $`c`$ et les coordonnées réelles de $`b`$ sont indépendantes, $`a,c\sim N(0,1)`$ et
chaque coordonnée réelle de $`b`$ suit la loi $`N(0,\frac12)`$ ; donc $`a-c\sim N(0,2)`$
(proposition A.14(c)) et $`2\mathrm{Re}\,b`$, $`2\mathrm{Im}\,b`$ suivent la loi $`N(0,2)`$, toutes
étant indépendantes. *GOE :* $`x=a-c`$ et $`y=2b`$ sont des variables $`N(0,2)`$ indépendantes, et en
coordonnées polaires (corollaire A.19 avec $`n=2`$, $`s_2=2\pi`$)
$`P(s\le t)=\int_{x^2+y^2\le t^2}\frac{e^{-(x^2+y^2)/4}}{4\pi}\,dx\,dy=\int_0^t\frac r2e^{-r^2/4}dr=1-e^{-t^2/4}`$.
Alors $`\mathbb Es=\int_0^\infty P(s>t)\,dt=\int_0^\infty e^{-t^2/4}dt=\sqrt\pi`$, et
$`P(S\le t)=P(s\le\sqrt\pi t)=1-e^{-\pi t^2/4}`$ ; il reste à dériver. *GUE :* $`x=a-c`$,
$`y=2\mathrm{Re}\,b`$, $`y'=2\mathrm{Im}\,b`$ sont des $`N(0,2)`$ indépendantes, et d'après le
théorème A.20 ($`n=3`$, $`s_3=4\pi`$) la densité de $`s=\sqrt{x^2+y^2+y'^2}`$ est
$`4\pi r^2(4\pi)^{-3/2}e^{-r^2/4}=\frac{r^2e^{-r^2/4}}{2\sqrt\pi}`$. Sa moyenne vaut
$`\frac1{2\sqrt\pi}\int_0^\infty r^3e^{-r^2/4}dr=\frac{8}{2\sqrt\pi}=\frac4{\sqrt\pi}`$. Le changement
de variable $`r=4S/\sqrt\pi`$ donne $`p_2`$, et l'intégration par parties de $`p_2`$
($`\int x^2e^{-ax^2}dx=-\frac{x}{2a}e^{-ax^2}+\frac1{2a}\int e^{-ax^2}dx`$ avec $`a=4/\pi`$) donne la
fonction de répartition. Le comportement pour $`\varepsilon`$ petit découle de
$`p_\beta(S)\sim cS^\beta`$. $`\square`$

*Calcul :* la proposition 5.13 par Monte-Carlo pour le GOE et le GUE $`2\times2`$ — [tutoriel, section 5.3](../../tutorial/fr/05_random_matrices.md#s-5-3).

L'exposant compte des dimensions : une matrice symétrique réelle $`2\times2`$ est dégénérée si et
seulement si deux conditions réelles sont satisfaites ($`a=c`$, $`b=0`$), une matrice hermitienne
complexe si et seulement si trois le sont. Pour tout $`n`$, les valeurs propres ordonnées
$`\lambda_1<\dots<\lambda_n`$ du GOE ($`\beta=1`$) et du GUE ($`\beta=2`$) ont la densité jointe
```math
Z_{n,\beta}^{-1}\prod_{i<j}(\lambda_j-\lambda_i)^\beta\,e^{-\sum_i\lambda_i^2/2}
```
(théorème C.7, démontré dans l'annexe C au moyen de la formule d'intégration de Weyl), et le plus
petit espacement $`\delta`$ vérifie $`c\,\varepsilon^{\beta+1}\le P(\delta<\varepsilon)\le C\varepsilon^{\beta+1}`$
pour $`0<\varepsilon\le1`$, avec des constantes $`0<c\le C`$ dépendant de $`n`$ et de $`\beta`$
(proposition C.9). La limite $`n\to\infty`$ n'est pas traitée dans ce cours (voir M. L. Mehta,
*Random Matrices*, 3e éd., Elsevier, 2004).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Wigner a
> proposé que les niveaux d'énergie des systèmes quantiques complexes (noyaux lourds) se comportent
> statistiquement comme les valeurs propres d'un hamiltonien aléatoire ayant la bonne symétrie :
> symétrique réel s'il existe une symétrie par renversement du temps ($`\beta=1`$), hermitien complexe
> sinon ($`\beta=2`$). Les niveaux se *repoussent* alors : deux niveaux de même symétrie ne se
> rapprochent presque jamais (von Neumann et Wigner, 1929). La répulsion est la signature typique du
> chaos quantique, des niveaux indépendants (sans répulsion) celle des systèmes intégrables ou
> localisés ; ce sont les conjectures de Bohigas–Giannoni–Schmit et de Berry–Tabor de la physique,
> qui ne sont ni démontrées ni utilisées dans ce cours.

Les espacements dépendent de la densité locale de niveaux, qui varie le long du spectre ; pour les
comparer à $`p_\beta`$, il faut procéder à un **dépliement** (*unfolding*), le changement de variable
$`\lambda\mapsto\bar N(\lambda)`$ par une approximation lisse $`\bar N`$ de la fonction de comptage, qui
rend l'espacement moyen égal à $`1`$. Le résultat dépend du choix de $`\bar N`$. Le rapport
d'espacements consécutifs évite ce choix.

**Définition.** Le **rapport d'espacements** est
$`r_i=\min(s_i,s_{i+1})/\max(s_i,s_{i+1})\in[0,1]`$ (V. Oganesyan, D. A. Huse, Phys. Rev. B 75
(2007), 155111). Il est invariant par $`\lambda\mapsto a\lambda+b`$, $`a>0`$, et ne dépend donc pas de
l'espacement moyen local ; c'est pourquoi on l'utilise sans dépliement.

*Dans la bibliothèque :* le rapport d'espacements moyen de chaque ligne de niveaux est `SpinChain.level_spacing_ratio` — [tutoriel, section 5.3](../../tutorial/fr/05_random_matrices.md#s-5-3).

Un rapport fait intervenir trois niveaux, de sorte que les plus petites matrices gaussiennes
possédant un rapport d'espacements sont de taille $`3\times3`$ ; pour elles, le rapport moyen se
calcule exactement, comme la surmise pour $`2\times2`$.

<a id="r-5-14"></a>

**Proposition 5.14 (rapports d'espacements moyens).** (a) *Niveaux indépendants.* Soient $`s,s'`$
indépendantes, de même loi exponentielle (comme deux espacements consécutifs des niveaux
$`e_0+\dots+e_{k-1}`$ construits à partir d'exponentielles indépendantes $`e_j`$, section A.7). Alors
$`r=\min(s,s')/\max(s,s')`$ vérifie $`P(r\le t)=\frac{2t}{1+t}`$ sur $`[0,1]`$, a pour densité
$`\frac2{(1+t)^2}`$, et
```math
\mathbb Er=2\ln2-1\approx0.3863 .
```
(b) *GOE et GUE $`3\times3`$.* Soit $`H`$ une matrice $`3\times3`$ du GOE ($`\beta=1`$) ou du GUE
($`\beta=2`$) de valeurs propres $`\lambda_1<\lambda_2<\lambda_3`$, d'espacements
$`s_1=\lambda_2-\lambda_1`$, $`s_2=\lambda_3-\lambda_2`$, et soit $`r=\min(s_1,s_2)/\max(s_1,s_2)`$. Alors
$`x=s_2/s_1`$ a une densité proportionnelle à
```math
f_\beta(x)=\frac{(x+x^2)^\beta}{(1+x+x^2)^{1+3\beta/2}}\qquad\text{on }(0,\infty),
```
$`r`$ a la densité $`f_\beta/I_0`$ sur $`(0,1)`$, où $`I_0=\int_0^1f_\beta`$, et
```math
\mathbb Er=4-2\sqrt3\approx0.5359\quad(\text{GOE}),\qquad
\mathbb Er=\frac{2\sqrt3}\pi-\frac12\approx0.6027\quad(\text{GUE}).
```

*Démonstration.* (a) Par changement d'échelle, on peut supposer le paramètre égal à $`1`$. Pour
$`x=s/s'`$, par indépendance et d'après le théorème de Tonelli,
$`P(x\le t)=\int_0^\infty e^{-y}(1-e^{-ty})\,dy=1-\frac1{1+t}=\frac t{1+t}`$.
Pour $`t<1`$, les événements $`\{x\le t\}`$ et $`\{1/x\le t\}`$ sont disjoints, $`x`$ et $`1/x`$ ont même loi
($`s,s'`$ sont échangeables), et $`\{r\le t\}`$ est leur réunion, donc $`P(r\le t)=\frac{2t}{1+t}`$.
Alors
$`\mathbb Er=\int_0^1\frac{2t}{(1+t)^2}dt=2\Big[\ln(1+t)+\frac1{1+t}\Big]_0^1=2\ln2-1`$.

(b) *La loi de $`x`$.* D'après le théorème C.7, le vecteur $`\lambda=(\lambda_1,\lambda_2,\lambda_3)`$ a
la densité $`Z^{-1}\Delta(\lambda)^\beta e^{-|\lambda|^2/2}`$ sur $`C=\{\lambda_1<\lambda_2<\lambda_3\}`$,
où $`\Delta(\lambda)=\prod_{i<j}(\lambda_j-\lambda_i)`$. L'application linéaire
```math
(m,s_1,s_2)\mapsto\Big(m-\tfrac{2s_1+s_2}3,\ m+\tfrac{s_1-s_2}3,\ m+\tfrac{s_1+2s_2}3\Big)
```
a pour déterminant $`1`$ et envoie $`\mathbb R\times(0,\infty)^2`$ sur $`C`$ ; son inverse est
$`m=\frac13\sum_i\lambda_i`$, $`s_1=\lambda_2-\lambda_1`$, $`s_2=\lambda_3-\lambda_2`$. Dans ces variables,
$`\Delta=s_1s_2(s_1+s_2)`$ et $`|\lambda|^2=3m^2+\frac23(s_1^2+s_1s_2+s_2^2)`$ (en développant). D'après
le lemme A.3(c) et le théorème de Fubini (en intégrant en $`m`$), $`(s_1,s_2)`$ a la densité
$`c\,\big(s_1s_2(s_1+s_2)\big)^\beta e^{-(s_1^2+s_1s_2+s_2^2)/3}`$ sur $`(0,\infty)^2`$ pour une
constante $`c>0`$. Posons $`q(x)=1+x+x^2`$. Pour tout borélien $`A\subseteq(0,\infty)`$, en substituant
$`s_2=xs_1`$ dans l'intégrale intérieure et en utilisant le théorème de Tonelli,
```math
P(x\in A)=c\int_0^\infty\!\!\int_A s_1^{3\beta+1}(x+x^2)^\beta e^{-s_1^2q(x)/3}\,dx\,ds_1
=c\,\gamma_\beta\int_A(x+x^2)^\beta\Big(\frac3{q(x)}\Big)^{1+3\beta/2}dx ,
```
où le changement de variable $`s_1=(3/q(x))^{1/2}u`$ donne le facteur, avec
$`\gamma_\beta=\int_0^\infty u^{3\beta+1}e^{-u^2}du\in(0,\infty)`$. Ainsi $`x`$ a la densité
$`f_\beta/\int_0^\infty f_\beta`$.

*Réduction à $`[0,1]`$.* $`f_\beta(1/x)=x^2f_\beta(x)`$ (multiplier numérateur et dénominateur par des
puissances de $`x`$). Avec $`x=1/y`$,
$`\int_{1/A}f_\beta(x)\,dx=\int_Af_\beta(1/y)y^{-2}dy=\int_Af_\beta`$ pour tout borélien
$`A\subseteq(0,1)`$, où $`1/A=\{1/y:y\in A\}`$ ; en particulier $`\int_0^\infty f_\beta=2I_0`$ avec
$`I_j=\int_0^1x^jf_\beta(x)\,dx`$. Comme $`r=\min(x,1/x)`$ et $`P(x=1)=0`$,
$`\{r\in A\}=\{x\in A\}\cup\{x\in1/A\}`$ (réunion disjointe), donc
```math
P(r\in A)=\frac{2\int_Af_\beta}{2I_0}=\int_A\frac{f_\beta}{I_0},\qquad \mathbb Er=\frac{I_1}{I_0} .
```
Comme $`x+x^2=q-1`$, $`f_\beta=g_\beta(q)`$ avec $`g_1(u)=u^{-3/2}-u^{-5/2}`$ et $`g_2(u)=(u-1)^2u^{-4}`$.
Écrivons $`x=\frac12(q'(x)-1)`$ avec $`q'=2x+1`$, et posons $`u=q(x)`$ (de $`q(0)=1`$ à $`q(1)=3`$) :
```math
I_1=\frac12\int_1^3g_\beta(u)\,du-\frac12I_0,\qquad \mathbb Er=\frac{\int_1^3g_\beta}{2I_0}-\frac12 .
```
Ici $`\int_1^3g_1=\big[-2u^{-1/2}+\frac23u^{-3/2}\big]_1^3=\frac43-\frac{16}{9\sqrt3}`$ et
$`\int_1^3g_2=\big[-u^{-1}+u^{-2}-\frac13u^{-3}\big]_1^3=\frac8{81}`$.

*Les intégrales $`I_0`$.* Posons $`L_a=\int_0^1q^{-a}dx`$. De $`(2x+1)^2=4q-3`$ on tire
$`\frac d{dx}\big[(2x+1)q^{-a}\big]=2q^{-a}-a(2x+1)^2q^{-a-1}=(2-4a)q^{-a}+3aq^{-a-1}`$, et en intégrant
sur $`[0,1]`$,
```math
3a\,L_{a+1}=3^{1-a}-1+(4a-2)L_a .
```
GOE : $`f_1=q^{-3/2}-q^{-5/2}`$. Avec $`a=\frac12`$, $`L_{3/2}=\frac23(\sqrt3-1)`$ ; avec $`a=\frac32`$,
$`L_{5/2}=\frac29\big(\frac1{\sqrt3}-1+4L_{3/2}\big)=\frac{18\sqrt3-22}{27}`$. Donc
$`I_0=L_{3/2}-L_{5/2}=\frac4{27}`$ et
$`\mathbb Er=\frac{27}8\big(\frac43-\frac{16}{9\sqrt3}\big)-\frac12=\frac92-2\sqrt3-\frac12=4-2\sqrt3`$.
GUE : $`f_2=q^{-2}-2q^{-3}+q^{-4}`$. Comme $`q=(x+\frac12)^2+\frac34`$,
$`L_1=\frac2{\sqrt3}\big[\arctan\frac{2x+1}{\sqrt3}\big]_0^1=\frac2{\sqrt3}\big(\frac\pi3-\frac\pi6\big)=\frac\pi{3\sqrt3}`$, et la
récurrence avec $`a=1,2,3`$ donne $`L_2=\frac23L_1=\frac{2\pi}{9\sqrt3}`$, $`L_3=L_2-\frac19`$, $`L_4=\frac19\big(10L_3-\frac89\big)`$. Donc
$`I_0=L_2-2L_3+L_4=\frac{2\pi}{81\sqrt3}`$ et $`\mathbb Er=\frac{81\sqrt3}{4\pi}\cdot\frac8{81}-\frac12=\frac{2\sqrt3}\pi-\frac12`$. $`\square`$

Les valeurs de (b) ont été obtenues par Y. Y. Atas, E. Bogomolny, O. Giraud et G. Roux (Phys. Rev.
Lett. 110 (2013), 084101). Pour $`n`$ plus grand, le rapport moyen n'est pas calculé dans ce cours.

La partie (a) vaut exactement pour $`n`$ niveaux uniformes indépendants sur un intervalle : leurs
$`n+1`$ espacements (y compris les deux lacunes aux extrémités) ont la loi de
$`(e_0,\dots,e_n)/\sum_je_j`$ pour des exponentielles indépendantes $`e_j`$, de sorte que chaque
rapport $`s_i/s_{i+1}`$ d'espacements consécutifs a la loi de $`e_i/e_{i+1}`$, et chaque rapport
d'espacements $`r_i`$ a la loi de $`r`$ dans (a) (théorème A.25 et corollaire A.26(b),(d)).

*Calcul :* la proposition 5.14 par Monte-Carlo et par quadrature de $`f_\beta`$, et le rapport moyen du GOE avec $`n=64`$, observé numériquement — [tutoriel, section 5.3](../../tutorial/fr/05_random_matrices.md#s-5-3).

<a id="s-5-4"></a>

## 5.4 Pourquoi il ne faut pas mélanger les secteurs

Soit $`H=\bigoplus_cH_c\in\mathcal M`$ un élément autoadjoint aléatoire à **blocs indépendants**, par
exemple des blocs GOE. En tant qu'opérateur sur $`\mathbb C^D`$, son spectre est la réunion des
spectres des blocs. À l'intérieur d'un bloc, les niveaux se repoussent ; nous montrons que les
niveaux de blocs différents ne se repoussent pas, et qu'avec de nombreux blocs indépendants le spectre
fusionné a les rapports d'espacements de niveaux indépendants.

Supposons que chaque bloc ait une **densité de niveaux** $`\rho_c`$ continue et bornée :
$`\mathbb E\sum_i\varphi(\lambda_i(H_c))=\int\varphi\,\rho_c`$ pour toute $`\varphi`$ borélienne bornée
(de sorte que $`\int\rho_c=k_c`$). Pour des blocs GOE et GUE, c'est le corollaire C.8, qui déduit
$`\rho_c`$ de la densité jointe du théorème C.7.

<a id="r-5-15"></a>

**Théorème 5.15 (pas de répulsion entre secteurs).** Pour $`c\neq d`$ et toute
$`\psi:\mathbb R^2\to\mathbb R`$ borélienne bornée,
```math
\mathbb E\sum_{i,j}\psi\big(\lambda_i(H_c),\lambda_j(H_d)\big)=\iint\psi(x,y)\,\rho_c(x)\rho_d(y)\,dx\,dy .
```
Par conséquent, si $`N_\varepsilon`$ est le nombre de paires de niveaux de blocs différents à
distance $`<\varepsilon`$, alors $`\mathbb EN_\varepsilon/(2\varepsilon)\to\sum_{c<d}\int\rho_c\rho_d`$
quand $`\varepsilon\to0`$, et le plus petit espacement $`\delta`$ du spectre fusionné vérifie
```math
P(\delta<\varepsilon)\ \ge\ \frac{\mathbb EN_\varepsilon}{\sum_{c<d}k_ck_d}
=\frac{2\varepsilon\sum_{c<d}\int\rho_c\rho_d}{\sum_{c<d}k_ck_d}+o(\varepsilon).
```

*Démonstration.* Les valeurs propres ordonnées sont des fonctions continues de la matrice
(corollaire A.11), donc $`\Lambda_c=(\lambda_i(H_c))_i`$ et $`\Lambda_d`$ sont des vecteurs aléatoires
indépendants. D'après le théorème de Fubini (leur loi jointe est le produit de leurs lois),
$`\mathbb E\sum_{i,j}\psi(\Lambda_{c,i},\Lambda_{d,j})=\mathbb E\sum_j g(\Lambda_{d,j})`$ avec
$`g(y)=\mathbb E\sum_i\psi(\Lambda_{c,i},y)=\int\psi(x,y)\rho_c(x)\,dx`$, qui est borélienne bornée ;
on applique la définition de $`\rho_d`$ à $`g`$. Avec $`\psi(x,y)=1\{|x-y|<\varepsilon\}`$,
```math
\frac{\mathbb EN_\varepsilon}{2\varepsilon}=\sum_{c<d}\int\rho_c(x)\Big(\frac1{2\varepsilon}\int_{x-\varepsilon}^{x+\varepsilon}\rho_d(y)\,dy\Big)dx ;
```
la parenthèse tend vers $`\rho_d(x)`$ par continuité et est majorée par $`\sup\rho_d`$, et $`\rho_c`$ est
intégrable, donc la convergence dominée donne la limite. Enfin $`N_\varepsilon\le\sum_{c<d}k_ck_d`$,
donc $`\mathbb EN_\varepsilon\le\big(\sum_{c<d}k_ck_d\big)P(N_\varepsilon\ge1)`$, et
$`N_\varepsilon\ge1`$ entraîne $`\delta<\varepsilon`$. $`\square`$

Ainsi, pour deux blocs GOE ou GUE, $`P(\delta<\varepsilon)`$ est au moins d'ordre $`\varepsilon`$ :
leurs densités de niveaux sont strictement positives partout (dans le corollaire C.8,
$`q(x,\cdot)>0`$ hors de l'ensemble négligeable où deux coordonnées de $`(x,y)`$ coïncident), donc
$`\int\rho_c\rho_d>0`$. À l'intérieur d'un bloc GOE $`2\times2`$, elle vaut
$`1-e^{-\varepsilon^2/4}\le\varepsilon^2/4`$ (proposition 5.13, espacement non normalisé), et à
l'intérieur d'un bloc GOE $`n\times n`$ elle est au plus une constante fois $`\varepsilon^2`$
(proposition C.9). Le test de répulsion des niveaux échoue pour le spectre fusionné, bien que
chaque bloc soit un GOE parfait.

Avec de nombreux secteurs, le spectre fusionné perd toute trace des blocs : ses rapports d'espacements deviennent ceux
de niveaux indépendants. Nous le démontrons, avec une vitesse, pour $`C`$ copies indépendantes $`H_1,\dots,H_C`$ d'une même
matrice autoadjointe aléatoire $`H`$ de taille $`k\times k`$, avec $`k`$ fixé et $`C\to\infty`$, pour les niveaux d'une fenêtre
où la densité de niveaux est strictement positive. Deux hypothèses sur $`H`$ sont utilisées :

- **(H1)** sa densité de niveaux $`\rho`$ est bornée et lipschitzienne : $`\rho\le\rho_{\max}`$ et $`|\rho(x)-\rho(y)|\le L|x-y|`$ ;
- **(H2)** son plus petit espacement $`\delta=\min_i(\lambda_{i+1}-\lambda_i)`$ ($`\delta=+\infty`$ si $`k=1`$) vérifie
  $`P(\delta\le\varepsilon)\le K\varepsilon^2`$ pour tout $`\varepsilon\ge0`$.

(H2) est une forme faible de répulsion des niveaux à l'intérieur du bloc ; avec $`\varepsilon=0`$, elle dit que le spectre de
$`H`$ est simple presque sûrement.

<a id="r-5-16"></a>

**Lemme 5.16 (blocs GOE et GUE).** Les matrices GOE et GUE de toute taille $`n\ge1`$ vérifient (H1) et (H2), et leur
densité de niveaux est strictement positive partout.

*Démonstration.* La positivité a été montrée après le théorème 5.15. *(H2).* Pour $`n=1`$, il n'y a rien à démontrer. Pour
$`n\ge2`$ et $`0\le\varepsilon<1`$, la proposition C.9 fournit une constante $`K_0`$ telle que
$`P(\delta\le\varepsilon)\le P(\delta<\varepsilon')\le K_0\varepsilon'^{\beta+1}\le K_0\varepsilon'^2`$ pour tout $`\varepsilon'\in(\varepsilon,1]`$ ; faire
$`\varepsilon'\downarrow\varepsilon`$. Avec $`K=\max(K_0,1)`$, la borne vaut aussi pour $`\varepsilon\ge1`$.
*(H1).* $`\rho`$ est bornée par le corollaire C.8. Dans ses notations, avec $`x\in\mathbb R`$ et $`y\in\mathbb R^{n-1}`$,
$`|\Delta(x,y)|=|\Delta(y)|\prod_i|x-y_i|`$, donc
```math
\rho(x)=\frac{n}{n!\,Z_{n,\beta}}\int_{\mathbb R^{n-1}}|\Delta(y)|^\beta e^{-|y|^2/2}\,g_y(x)\,dy,\qquad
g_y(x)=f_0(x)\prod_{i=1}^{n-1}f_i(x),
```
avec $`f_0(x)=e^{-x^2/(2n)}`$ et $`f_i(x)=|x-y_i|^\beta e^{-x^2/(2n)}`$ (pour $`n=1`$ il n'y a pas de $`y`$, et $`\rho=f_0/Z_{1,\beta}`$).
Soit $`a=(\beta+1)\sup_x(1+|x|)^{\beta+1}e^{-x^2/(2n)}<\infty`$. De $`|x-y_i|\le(1+|x|)(1+|y_i|)`$ on tire $`|f_i|\le a(1+|y_i|)^\beta`$, et
$`f_i`$ est continue, de classe $`C^1`$ hors de $`y_i`$, avec
$`|f_i'(x)|\le\big(\beta|x-y_i|^{\beta-1}+|x||x-y_i|^\beta\big)e^{-x^2/(2n)}\le a(1+|y_i|)^\beta`$. Par l'inégalité des accroissements finis de
chaque côté de $`y_i`$, $`f_i`$ est lipschitzienne de constante $`a(1+|y_i|)^\beta`$ ; de même $`f_0`$ de constante $`a`$, et
$`|f_0|\le1\le a`$. En écrivant la différence de deux produits comme une somme télescopique,
$`|g_y(x)-g_y(x')|\le\sum_i|f_i(x)-f_i(x')|\prod_{j\neq i}\sup|f_j|\le na^n\prod_i(1+|y_i|)^\beta\,|x-x'|`$. Donc $`\rho`$ est
lipschitzienne avec $`L=\frac{n^2a^n}{n!\,Z_{n,\beta}}\int|\Delta(y)|^\beta\prod_i(1+|y_i|)^\beta e^{-|y|^2/2}dy`$, qui est fini grâce à la
borne $`|\Delta(y)|\le\prod_i(1+|y_i|)^{n-2}`$ de la démonstration du théorème C.7. $`\square`$

Soit $`\Lambda_c`$ l'ensemble des niveaux de $`H_c`$ ; le **spectre fusionné de $`m`$ blocs** est $`\Lambda_1\cup\dots\cup\Lambda_m`$. Pour
un point $`x`$, soient $`L_m(x)`$ et $`R_m(x)`$ les distances de $`x`$ au niveau le plus proche du spectre fusionné de $`m`$
blocs au-dessous et au-dessus de $`x`$ ($`+\infty`$ s'il n'y en a pas) et, lorsque les deux sont finies,
$`r_m(x)=\min(L_m(x),R_m(x))/\max(L_m(x),R_m(x))`$. Pour un intervalle $`J`$ de longueur $`|J|`$, soient $`N(J)`$ le nombre de
niveaux de $`H`$ dans $`J`$, $`\rho(J)=\int_J\rho=\mathbb EN(J)`$ et $`\sigma(J)=P(N(J)=0)`$.

<a id="r-5-17"></a>

**Lemme 5.17 (les deux lacunes autour d'un point fixé).** Supposons (H1) et (H2), et soient $`\ell,\rho_*>0`$ et
$`x\in\mathbb R`$ avec $`\rho\ge\rho_*`$ sur $`[x-\ell,x+\ell]`$. Il existe une constante $`A`$, ne dépendant que de
$`k,K,L,\rho_{\max},\rho_*,\ell`$, telle que pour tout $`m\ge1`$ et $`0\le t\le1`$
```math
\Big|P\big(L_m(x),R_m(x)<\infty,\ r_m(x)\le t\big)-\frac{2t}{1+t}\Big|\le\frac Am .
```

*Démonstration.* Pour $`m\le2`$, tout $`A\ge2`$ convient ; soit $`m\ge3`$. *Un bloc.* Comme
$`1\{N\ge1\}\le N\le k\,1\{N\ge1\}`$, en prenant les espérances on obtient (1) : $`1-\rho(J)\le\sigma(J)\le1-\rho(J)/k\le e^{-\rho(J)/k}`$.
Si $`N(J)\ge2`$, deux niveaux sont dans $`J`$ et $`\delta\le|J|`$ ; comme $`N-1\{N\ge1\}\le k\,1\{N\ge2\}`$, (H2) donne (2) :
$`\sigma(J)\le1-\rho(J)+kK|J|^2`$. Pour $`c>0`$, d'après le théorème de Tonelli et (H2),
```math
\mathbb Ee^{-c\delta}=\int_0^1P\big(e^{-c\delta}>u\big)\,du\le\int_0^1\frac{K\ln^2u}{c^2}\,du=\frac{2K}{c^2}.\qquad(3)
```
*Simplicité.* D'après (H2) avec $`\varepsilon=0`$, chaque bloc a un spectre simple presque sûrement ; d'après le théorème 5.15
avec $`\psi(u,v)=1\{u=v\}`$ (la diagonale est négligeable), deux blocs n'ont presque sûrement aucun niveau commun ; et $`x`$
n'est presque sûrement pas un niveau, car $`\mathbb EN(\{x\})=0`$.

*Une formule.* Notons $`\Lambda`$ le spectre fusionné de $`m`$ blocs, $`L,R`$ pour $`L_m(x),R_m(x)`$, fixons $`t`$ et posons
$`J_b=[x-tb,x+b)`$ et $`\sigma_b=\sigma(J_b)`$ pour $`b>0`$ ; $`J_b`$ croît avec $`b`$, donc $`\sigma_b`$ est décroissante.
Presque sûrement
```math
1\{R<\infty,\ L>tR\}=\sum_{y\in\Lambda,\ y>x}1\{\Lambda\cap J_{y-x}=\{y\}\},
```
car pour $`y=x+R`$ les autres niveaux dans $`J_R`$ sont ceux de $`[x-tR,x)`$, qui existent si et seulement si $`L\le tR`$, et pour
$`y>x+R`$ le niveau $`x+R`$ est dans $`J_{y-x}`$. Séparons la somme selon le bloc de $`y`$. Les blocs sont indépendants et de
même loi, et, les niveaux du bloc de $`y`$ étant donnés, les $`m-1`$ autres blocs évitent $`J_{y-x}`$ avec probabilité
$`\sigma_{y-x}^{m-1}`$ ; d'après le théorème de Fubini, comme dans la démonstration du théorème 5.15,
```math
P(R<\infty,\ L>tR)=m\,\mathbb E\sum_{y\in\Lambda_1,\ y>x}1\{\Lambda_1\cap J_{y-x}=\{y\}\}\,\sigma_{y-x}^{m-1}=M-E_1,
```
où $`M=m\,\mathbb E\sum_{y\in\Lambda_1,y>x}\sigma_{y-x}^{m-1}=m\int_0^\infty\rho(x+b)\sigma_b^{m-1}db`$ (définition de $`\rho`$) et $`E_1\ge0`$
est la même espérance avec $`1\{\Lambda_1\cap J_{y-x}\neq\{y\}\}`$.

*$`E_1`$ est petit.* Sur cet événement, $`H_1`$ a un second niveau à distance $`\le|J_{y-x}|\le2(y-x)`$ de $`y`$, donc
$`y-x\ge\delta_1/2`$, où $`\delta_1`$ est le plus petit espacement de $`H_1`$. D'après (1), $`\sigma_b\le e^{-\rho_*\min(b,\ell)/k}`$, car
$`J_b\supseteq[x,x+\min(b,\ell))`$. Avec au plus $`k`$ termes, $`\mu=(m-1)\rho_*/k`$ et (3),
$`E_1\le mk\,\mathbb Ee^{-\mu\min(\delta_1/2,\ell)}\le mk\big(8Kk^2/((m-1)\rho_*)^2+e^{-(m-1)\rho_*\ell/k}\big)`$.

*$`M`$ est proche de $`\frac1{1+t}`$.* Soient $`\rho_0=\rho(x)`$, $`a=(1+t)\rho_0`$ et $`b_1=\min(\ell,1/a)\ge b_0=\min(\ell,1/(2\rho_{\max}))`$.
Alors $`m\int_0^{1/a}\rho_0(1-ab)^{m-1}db=\rho_0/a=\frac1{1+t}`$. Pour $`0<b\le b_1`$ :
$`|\rho(J_b)-ab|\le\int_{J_b}|\rho-\rho_0|\le L\frac{(1+t^2)b^2}2\le Lb^2`$, donc d'après (1) et (2) $`|\sigma_b-(1-ab)|\le K_1b^2`$ avec
$`K_1=L+4kK`$ ; les deux nombres sont dans $`[0,e^{-\rho_*b/k}]`$ ; et $`|\rho(x+b)-\rho_0|\le Lb`$. D'après
$`|u^n-v^n|\le n|u-v|\max(u,v)^{n-1}`$ pour $`u,v\ge0`$,
```math
\big|\rho(x+b)\sigma_b^{m-1}-\rho_0(1-ab)^{m-1}\big|\le\rho_{\max}(m-1)K_1b^2e^{-(m-2)\rho_*b/k}+Lb\,e^{-(m-1)\rho_*b/k}.
```
Avec $`\int_0^\infty b^je^{-\lambda b}db=j!/\lambda^{j+1}`$, $`m`$ fois l'intégrale de ceci sur $`(0,b_1]`$ est au plus
$`2m(m-1)\rho_{\max}K_1k^3/((m-2)\rho_*)^3+mLk^2/((m-1)\rho_*)^2`$. Les parties restantes sont
$`m\int_{b_1}^\infty\rho(x+b)\sigma_b^{m-1}db\le mk\,e^{-(m-1)\rho_*b_0/k}`$ (car $`\int\rho=k`$ et $`\sigma_b\le\sigma_{b_1}`$) et
$`m\int_{b_1}^{1/a}\rho_0(1-ab)^{m-1}db\le(1-ab_1)^m\le e^{-m\rho_*b_0}`$.

*Des deux côtés.* Le même argument pour les matrices $`-H_c`$ (densité de niveaux $`\rho(-\,\cdot)`$, mêmes hypothèses) au point
$`-x`$ donne $`P(L<\infty,\ R>tL)=M'-E_1'`$ avec les mêmes bornes. Pour $`t<1`$, les événements $`\{R<\infty,\ L\le tR\}`$ et
$`\{L<\infty,\ R\le tL\}`$ sont disjoints (ensemble ils donneraient $`L\le t^2L`$, alors que $`L>0`$), et leur réunion est
$`\{L,R<\infty,\ r_m(x)\le t\}`$, donc
```math
P\big(L,R<\infty,\ r_m(x)\le t\big)-\frac{2t}{1+t}=-P(R=\infty)-P(L=\infty)-\Big(M-\frac1{1+t}\Big)-\Big(M'-\frac1{1+t}\Big)+E_1+E_1'.
```
D'après (1), $`P(R=\infty)\le\sigma((x,x+\ell))^m\le e^{-m\rho_*\ell/k}`$, et de même pour $`L`$ ; pour $`t=1`$, la probabilité vaut
$`1-P(L=\infty\text{ ou }R=\infty)`$. Chaque terme est au plus une constante divisée par $`m`$, car
$`m(m-1)/(m-2)^3\le18/m`$, $`m/(m-1)^2\le4/m`$ et $`me^{-(m-1)c}\le4/(mc^2)`$ pour $`c>0`$ (d'après $`u^2e^{-u}\le1`$). $`\square`$

<a id="r-5-18"></a>

**Théorème 5.18 (de nombreux secteurs indépendants ressemblent à des niveaux indépendants).** Supposons (H1) et (H2).
Soient $`I=[\alpha,\gamma]`$ et $`\ell,\rho_*>0`$ avec $`\rho\ge\rho_*`$ sur $`[\alpha-\ell,\gamma+\ell]`$, soit $`C\ge2`$, et soit
$`\Lambda`$ le spectre fusionné de $`C`$ blocs (simple presque sûrement). Pour un niveau $`y`$ de voisins $`y_-<y<y_+`$ dans
$`\Lambda`$, soit $`r(y)=\min(y-y_-,y_+-y)/\max(y-y_-,y_+-y)`$. Soient $`N(I)`$ le nombre de niveaux dans $`I`$, $`N_t(I)`$ le
nombre de ceux qui ont deux voisins et $`r(y)\le t`$, et $`S(I)`$ la somme des $`r(y)`$ sur les niveaux de $`I`$ ayant deux
voisins. Alors $`\mathbb EN(I)=C\rho(I)`$, et il existe une constante $`A`$, ne dépendant que de $`k,K,L,\rho_{\max},\rho_*,\ell`$ et
$`\rho(I)`$, telle que pour $`0\le t\le1`$
```math
\Big|\frac{\mathbb EN_t(I)}{\mathbb EN(I)}-\frac{2t}{1+t}\Big|\le\frac AC,\qquad
\Big|\frac{\mathbb ES(I)}{\mathbb EN(I)}-(2\ln2-1)\Big|\le\frac AC .
```

Les rapports d'espérances sont ce qu'estime un calcul qui regroupe de nombreux échantillons : la proportion des niveaux de
$`I`$ ayant $`r\le t`$, et leur rapport d'espacements moyen. D'après le lemme 5.16, le théorème s'applique à des blocs GOE
ou GUE indépendants de toute taille fixée et à tout intervalle compact $`I`$.

*Démonstration.* $`\mathbb EN(I)=\sum_c\rho(I)`$. Les blocs ont la même loi, donc
$`\mathbb EN_t(I)=C\,\mathbb E\sum_{y\in\Lambda_1\cap I}1\{y\text{ a deux voisins},\ r(y)\le t\}`$. Soit $`m=C-1`$ et, pour
$`y\in\Lambda_1`$, soient $`L'(y),R'(y)`$ ses lacunes dans le spectre fusionné $`\Lambda'`$ de $`H_2,\dots,H_C`$ (comme dans le
lemme 5.17) et $`r'(y)`$ leur rapport. Les autres niveaux de $`H_1`$ sont à distance $`\ge\delta_1`$ de $`y`$ ; donc si
$`\max(L'(y),R'(y))<\delta_1`$, les voisins de $`y`$ dans $`\Lambda`$ sont ceux dans $`\Lambda'`$, et $`r(y)=r'(y)`$. Or $`\Lambda'`$ est
indépendant de $`H_1`$, et $`\max(L'(y),R'(y))\ge s`$ exige que $`\Lambda'`$ évite $`(y-s,y)`$ ou $`(y,y+s)`$, ce qui d'après (1) a
une probabilité au plus $`2e^{-m\rho_*\min(s,\ell)/k}`$ pour $`y\in I`$. D'après le théorème de Fubini (comme dans la
démonstration du théorème 5.15) et (3),
```math
\Big|\mathbb EN_t(I)-C\,\mathbb E\sum_{y\in\Lambda_1\cap I}p_m(y,t)\Big|\le2Ck\big(\mathbb Ee^{-m\rho_*\delta_1/k}+e^{-m\rho_*\ell/k}\big)
\le2Ck\Big(\frac{2Kk^2}{(m\rho_*)^2}+e^{-m\rho_*\ell/k}\Big),
```
où $`p_m(x,t)=P\big(L_m(x),R_m(x)<\infty,\ r_m(x)\le t\big)`$, fonction borélienne de $`x`$ d'après le théorème de Fubini.
D'après la définition de $`\rho`$, $`\mathbb E\sum_{y\in\Lambda_1\cap I}p_m(y,t)=\int_I\rho(x)p_m(x,t)\,dx`$, et le lemme 5.17
s'applique en tout $`x\in I`$ avec une même constante $`A_0`$. En divisant par $`C\rho(I)`$,
$`\big|\mathbb EN_t(I)/\mathbb EN(I)-\frac{2t}{1+t}\big|\le A_0/m+2k\big(2Kk^2/(m\rho_*)^2+e^{-m\rho_*\ell/k}\big)/\rho(I)`$, qui est au
plus une constante divisée par $`C`$ puisque $`m\ge C/2`$. Pour la moyenne, $`r(y)=\int_0^11\{r(y)>t\}\,dt`$, donc d'après le
théorème de Tonelli $`\mathbb ES(I)=\int_0^1\big(\mathbb EN_1(I)-\mathbb EN_t(I)\big)dt`$, tandis que
$`\int_0^1\big(1-\frac{2t}{1+t}\big)dt=\int_0^1\frac{1-t}{1+t}dt=2\ln2-1`$ (proposition 5.14(a)) ; la seconde déviation est donc au plus le
double de la borne de la première. $`\square`$

La démonstration montre le mécanisme. Un niveau $`y`$ d'un bloc voit les $`C-1`$ autres blocs comme une configuration
indépendante, dont les niveaux voisins de $`y`$ sont à des distances d'ordre $`1/(C\rho)`$ ; à de telles distances, les autres
niveaux du bloc de $`y`$ sont absents avec grande probabilité (c'est là qu'interviennent $`C\to\infty`$ et (H2)), de sorte que
$`r(y)`$ est le rapport dans lequel $`y`$ partage la lacune des autres blocs qui le contient. Pour une densité constante et au
plus un niveau par bloc près de $`y`$, la formule de la démonstration donne exactement la valeur des niveaux indépendants,
$`m\int_0^{1/a}\rho_0(1-ab)^{m-1}db=\frac1{1+t}`$, pour tout $`m`$ ; les erreurs proviennent de la variation de $`\rho`$ sur une
lacune et de deux niveaux d'un même bloc dans une lacune. La vitesse $`1/C`$ est celle que donne la démonstration ; pour
des blocs GOE, l'écart trouvé numériquement décroît comme $`1/C^2`$ (une observation numérique, dans le tutoriel). Que les
superpositions de nombreux processus ponctuels indépendants et clairsemés soient proches de processus de Poisson est un
thème classique des probabilités (B. Grigelionis, Theory Probab. Appl. 8 (1963), 177–182) ; le théorème 5.18 en est le
cas utile ici, avec une démonstration interne à ce cours. L'indépendance des blocs est essentielle à la démonstration : le
théorème ne s'applique pas aux secteurs de la chaîne de Heisenberg en champs aléatoires (chapitre 9), qui dépendent tous
des mêmes champs.

*Calcul :* le théorème 5.18 pour $`C`$ blocs GOE indépendants de taille $`64`$, avec une vérification de la vitesse $`1/C`$ ; et la chaîne de Heisenberg en champs aléatoires, dont les blocs ne sont pas indépendants, où le même effet est une observation numérique — [tutoriel, section 5.4](../../tutorial/fr/05_random_matrices.md#s-5-4).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Des niveaux
> portant des nombres quantiques conservés différents n'interagissent pas, ils se croisent donc
> librement ; seuls les niveaux de même symétrie se repoussent. Un spectre qui mélange des secteurs
> de symétrie présente par conséquent des statistiques de type Poisson même lorsque la dynamique
> dans chaque secteur est pleinement chaotique, et une symétrie non résolue peut être prise à tort
> pour de l'intégrabilité ou de la localisation. La règle de l'art — résoudre toutes les symétries
> avant de calculer des statistiques de niveaux — revient à dire qu'il faut calculer dans l'algèbre
> $`\bigoplus_cM_{k_c}`$, et non dans $`M_D`$.

<a id="s-5-5"></a>

## 5.5 Points de comparaison

*Dans la bibliothèque :* la comparaison avec d'autres logiciels — [tutoriel, section 5.5](../../tutorial/fr/05_random_matrices.md#s-5-5).

<a id="x-5"></a>

## Exercices

1. Soit $`U`$ distribué selon la loi de Haar sur $`U(\mathcal M)`$ et
   $`\mathrm{Tr}\,U=\sum_c\mathrm{Tr}\,U_c`$. Montrer que $`\mathbb E|\mathrm{Tr}\,U|^2=C`$
   et que $`\mathbb E|\mathrm{Tr}\,U|^4=\sum_cm(k_c)+2C(C-1)`$ avec $`m(1)=1`$, $`m(k)=2`$ pour
   $`k\ge2`$ ; comparer avec $`1`$ et $`2`$ pour un unitaire de Haar de $`U(D)`$, $`D\ge2`$, et évaluer pour
   $`\mathcal M=M_1\oplus M_2\oplus M_3`$.
2. À l'aide du théorème 5.12, montrer que $`\mathbb E|U_{11}|^{2n}=\frac{n!}{k(k+1)\cdots(k+n-1)}`$
   pour $`n\le k`$ (le corollaire B.14 étend ce résultat à $`n>k`$).
3. L'**ensemble orthogonal circulaire** est $`S=W^TW`$ avec $`W`$ distribué selon la
   loi de Haar sur $`U(k)`$. Montrer que $`S`$ est unitaire et symétrique, et que sa loi est invariante
   par $`S\mapsto V^TSV`$ pour tout $`V\in U(k)`$. Soient $`\lambda_1,\dots,\lambda_k`$ les valeurs propres
   de $`S`$ (avec multiplicité) et $`N_B(S)=\#\{j:\lambda_j\in B\}`$ pour les boréliens $`B`$ du cercle
   unité $`\mathbb T`$. Montrer que $`N_B(S)`$ est une variable aléatoire, que
   $`\nu(B)=\frac1k\mathbb EN_B(S)`$ est la longueur d'arc normalisée (la densité moyenne des phases
   propres est uniforme), et que $`\mathbb E\mathrm{Tr}\,S^m=0`$ pour tout entier $`m\neq0`$.

<details>
<summary>Solutions</summary>

1. Posons $`t_c=\mathrm{Tr}\,U_c`$ ; les $`t_c`$ sont indépendants (corollaire 5.4). Multiplier
   $`U_c`$ par $`e^{i\theta}`$ ne change pas sa loi et multiplie $`t_c`$ par $`e^{i\theta}`$ ; un moment
   $`\mathbb E[t_at_b\bar t_c\bar t_d]`$ est donc nul sauf si chaque secteur apparaît aussi souvent
   dans $`(a,b)`$ que dans $`(c,d)`$, c'est-à-dire si $`\{a,b\}=\{c,d\}`$ en tant que multiensembles, et
   de même $`\mathbb E[t_a\bar t_b]=0`$ pour $`a\neq b`$. D'après la proposition 5.7,
   $`\mathbb E[(U_c)_{ii}\overline{(U_c)_{jj}}]=\langle e_i,\mathbb E[U_ce_ie_j^*U_c^*]e_j\rangle=\delta_{ij}/k_c`$,
   donc $`\mathbb E|t_c|^2=1`$ et $`\mathbb E|\mathrm{Tr}\,U|^2=C`$. Dans
   $`\mathbb E|\sum_at_a|^4`$, les termes avec $`a=b`$ donnent $`\sum_a\mathbb E|t_a|^4`$, ceux avec
   $`a\neq b`$ donnent $`2\sum_{a\neq b}\mathbb E|t_a|^2\,\mathbb E|t_b|^2=2C(C-1)`$ (deux choix de
   $`(c,d)`$). Pour $`k=1`$, $`|t|=1`$. Pour $`k\ge2`$,
   $`\mathbb E|t|^4=\sum_{i,j,i',j'}\mathbb E[U_{ii}U_{jj}\overline{U_{i'i'}U_{j'j'}}]`$ ; dans le
   théorème 5.10, les deux termes avec $`\sigma=\pi`$ imposent $`\{i',j'\}=\{i,j\}`$ d'une seule façon
   et donnent chacun $`k^2\mathrm{Wg}(e)`$ après sommation, les deux termes avec $`\sigma\neq\pi`$
   imposent $`i=j=i'=j'`$ et donnent chacun $`k\,\mathrm{Wg}(\tau)`$ : au total
   $`2\big(\frac{k^2}{k^2-1}-\frac1{k^2-1}\big)=2`$. Le cas $`C=1`$, $`k=D`$ donne $`1`$ et $`2`$ pour
   $`U(D)`$. Pour $`M_1\oplus M_2\oplus M_3`$ : $`\mathbb E|\mathrm{Tr}\,U|^4=1+2+2+12=17`$, contre
   $`2`$ pour un unitaire de Haar de $`U(6)`$.

2. Lorsque tous les indices valent $`1`$, chaque delta vaut $`1`$, donc
   $`\mathbb E|U_{11}|^{2n}=\sum_{\sigma,\pi}\mathrm{Wg}(\sigma\pi^{-1})=\mathbf 1^TW\mathbf 1`$, où
   $`W=\Gamma^{-1}=(\mathrm{Wg}(\sigma^{-1}\pi))`$ et $`\Gamma=(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)})`$.
   On a utilisé ici que $`\mathrm{Wg}`$ ne dépend que du type cyclique, et que
   $`\sigma\pi^{-1}=\sigma(\pi^{-1}\sigma)\sigma^{-1}`$ est conjuguée à
   $`\pi^{-1}\sigma=(\sigma^{-1}\pi)^{-1}`$, qui a le type cyclique de $`\sigma^{-1}\pi`$. Chaque somme
   de ligne de $`\Gamma`$ vaut $`c_n(k)=\sum_{\rho\in S_n}k^{\#\mathrm{cycles}(\rho)}`$, donc
   $`\Gamma\mathbf 1=c_n(k)\mathbf 1`$, $`W\mathbf 1=\mathbf 1/c_n(k)`$ et
   $`\mathbf 1^TW\mathbf 1=n!/c_n(k)`$. Enfin $`c_n(x)=x(x+1)\cdots(x+n-1)`$ par récurrence : une
   permutation de $`\{1,\dots,n\}`$ s'obtient à partir d'une permutation de $`\{1,\dots,n-1\}`$ soit en
   ajoutant $`n`$ comme point fixe (un cycle de plus, facteur $`x`$), soit en insérant $`n`$ après l'un
   des $`n-1`$ éléments dans son cycle ($`n-1`$ façons, même nombre de cycles) ; c'est une bijection,
   donc $`c_n(x)=(x+n-1)c_{n-1}(x)`$.

3. $`S^*S=W^*\bar WW^TW=1`$ puisque $`W^T`$ est unitaire, et $`S^T=W^TW=S`$. Pour $`V\in U(k)`$, $`WV`$ a la
   loi de $`W`$ (invariance à droite, lemme 5.1), et $`(WV)^T(WV)=V^TSV`$. Avec $`V=e^{i\alpha/2}1`$ :
   $`e^{i\alpha}S`$ a la loi de $`S`$ pour tout réel $`\alpha`$.

   *Mesurabilité.* D'après le théorème 1.3, $`S=Y\mathrm{diag}(\lambda_1,\dots,\lambda_k)Y^*`$
   avec $`Y`$ unitaire. Si $`v`$ est un vecteur unitaire tel que $`Sv=\mu v`$, l'écriture
   $`v=\sum_jc_jYe_j`$ donne $`c_j(\lambda_j-\mu)=0`$, donc $`v`$ appartient à
   $`E_K=\mathrm{span}\{Ye_j:\lambda_j\in K\}`$ dès que $`\mu\in K`$ ; et $`\dim E_K=N_K(S)`$. Par
   conséquent, pour $`K\subseteq\mathbb T`$ et $`m\ge1`$ : $`N_K(S)\ge m`$ si et seulement si $`S`$ possède
   $`m`$ vecteurs propres orthonormés de valeurs propres dans $`K`$. Si $`K`$ est fermé, l'ensemble
   $`F_{K,m}`$ des unitaires ayant cette propriété est fermé dans $`U(k)`$ : si $`S_l\to S`$ et si $`S_l`$ a
   des vecteurs propres orthonormés $`v_i^{(l)}`$ de valeurs propres $`\mu_i^{(l)}\in K`$ ($`i\le m`$),
   une sous-suite vérifie $`v_i^{(l)}\to v_i`$ et $`\mu_i^{(l)}\to\mu_i\in K`$ pour tout $`i\le m`$
   (Bolzano–Weierstrass dans le compact $`(\text{sphère unité}\times K)^m`$), et les limites sont
   orthonormées avec $`Sv_i=\mu_iv_i`$.
   Soient $`0\le\theta_1(S)\le\dots\le\theta_k(S)<2\pi`$ les phases propres, de sorte que les valeurs
   propres sont les $`e^{i\theta_j(S)}`$ (corollaire 1.4 : $`|\lambda_j|=1`$). Pour $`0\le t<2\pi`$,
   l'arc $`K_t=\{e^{i\varphi}:0\le\varphi\le t\}`$ est fermé et
   $`\{\theta_j\le t\}=\{N_{K_t}\ge j\}=F_{K_t,j}`$, donc chaque $`\theta_j`$ est une fonction
   borélienne sur $`U(k)`$. Ainsi $`N_B(S)=\sum_j1_B(e^{i\theta_j(S)})`$ est une variable aléatoire,
   et $`\nu(B)=\frac1k\sum_jP(e^{i\theta_j(S)}\in B)`$ est une mesure de probabilité borélienne sur
   $`\mathbb T`$, la moyenne des lois des $`k`$ valeurs propres.

   *Uniformité.* Les valeurs propres de $`e^{i\alpha}S`$ sont les $`e^{i\alpha}\lambda_j`$, donc
   $`N_B(e^{i\alpha}S)=N_{e^{-i\alpha}B}(S)`$ et, comme $`e^{i\alpha}S`$ a la loi de $`S`$,
   $`\nu(e^{-i\alpha}B)=\frac1k\mathbb EN_B(e^{i\alpha}S)=\nu(B)`$. Ainsi $`\nu`$ est invariante par
   rotation et, d'après la proposition A.22 (avec $`N=1`$), c'est la longueur d'arc normalisée :
   $`\mathbb EN_B(S)=k\,|B|/2\pi`$ pour tout arc de longueur $`|B|`$. De même,
   $`\mathrm{Tr}\,S^m=\sum_j\lambda_j^m`$ (aussi pour $`m<0`$, puisque $`S^{-1}=S^*`$) est une
   fonction continue bornée de $`S`$, et
   $`\mathbb E\mathrm{Tr}\,S^m=\mathbb E\mathrm{Tr}(e^{i\alpha}S)^m=e^{im\alpha}\,\mathbb E\mathrm{Tr}\,S^m`$
   pour tout $`\alpha`$, donc cette quantité est nulle pour $`m\neq0`$ : ce sont $`k`$ fois les
   coefficients de Fourier de $`\nu`$. En particulier, la densité moyenne des phases propres est
   constante, et les espacements le long du cercle ne nécessitent aucun dépliement.

</details>

*Calcul :* les parties numériques des exercices 1–3 — [tutoriel, Exercices](../../tutorial/fr/05_random_matrices.md#x-5).

**Suite :** le [chapitre 6](06_quantum_channels.md) traite des canaux quantiques : représentations
de Kraus et de Stinespring, duaux relativement à une trace, et applications entre secteurs.
L'espérance conditionnelle $`E`$, rencontrée ici comme moyenne de Haar, en est le premier exemple.
