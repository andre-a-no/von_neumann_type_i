[English](../09_spin_chains.md) · [Español](../es/09_spin_chains.md) · **Français** · [Deutsch](../de/09_spin_chains.md) · [中文](../zh/09_spin_chains.md) · [日本語](../ja/09_spin_chains.md) · [Русский](../ru/09_spin_chains.md)

# 9. Chaînes de spins

Une chaîne de $`L`$ spins $`\frac12`$ vit sur $`(\mathbb C^2)^{\otimes L}`$, et l'algèbre linéaire dense en dimension
$`2^L`$ devient coûteuse au-delà de $`L\approx14`$ (une matrice dense a $`4^L`$ coefficients). Les symétries repoussent cette limite : l'aimantation découpe l'espace en $`L+1`$
secteurs (chapitre 4), les translations découpent à nouveau chaque secteur, et $`SU(2)`$ réduit encore les blocs.
Ce chapitre construit ces secteurs pour la chaîne XXZ, démontre ce qu'ils disent des états fondamentaux, confronte
la bibliothèque aux fermions libres et à l'ansatz de Bethe, et va au-delà des matrices denses (méthodes de Krylov)
et des chaînes propres (désordre).

**Vous apprendrez** pourquoi l'algèbre d'une chaîne est $`\bigoplus_NM_{\binom LN}`$ et pourquoi la tour de ces
algèbres a pour diagramme de Bratteli le triangle de Pascal et pour indice $`4\cos^2\frac\pi{2L+4}`$, comment on
dénombre les secteurs d'impulsion, la règle des signes de Marshall et ses conséquences pour l'impulsion et le spin
des états fondamentaux (Lieb–Mattis), la solution de Jordan–Wigner de la chaîne XX, les équations de Bethe à deux
magnons, pourquoi Lanczos converge de façon monotone, quelle est la précision d'un pas de temps de Krylov, et comment mesurer les statistiques de niveaux et
l'intrication des fermions libres. Les vérifications numériques de ce chapitre se trouvent dans le
[Tutoriel 9](../../tutorial/fr/09_spin_chains.md).

Dans tout ce chapitre, les sites sont $`0,\dots,L-1`$, la base de $`(\mathbb C^2)^{\otimes L}`$ est formée des
configurations $`|b\rangle=|b_0\dots b_{L-1}\rangle`$ avec $`b_i=1`$ pour un spin vers le haut (convention de la
bibliothèque : le bit $`i`$ d'une étiquette est le site $`i`$), $`s^a_i`$ ($`a=x,y,z`$) sont les matrices de spin
$`\frac12\sigma_a`$ sur le site $`i`$, $`s^\pm_i=s^x_i\pm is^y_i`$, $`n_i=s^z_i+\frac12`$, $`S^a=\sum_is^a_i`$, et
$`N=S^z+\frac L2`$ est le nombre de spins vers le haut. Une chaîne est **ouverte** (liaisons $`(i,i+1)`$, $`i<L-1`$) ou
**périodique** (avec aussi la liaison $`(L-1,0)`$ ; alors $`L\ge3`$ et les indices sont pris modulo $`L`$). $`J>0`$
(antiferromagnétique), $`\hbar=1`$, et $`\log`$ est le logarithme népérien.

<a id="s-9-1"></a>

## 9.1 La chaîne XXZ et ses secteurs $`U(1)`$

La **chaîne XXZ** d'anisotropie $`\Delta\in\mathbb R`$ et de champs $`h_i\in\mathbb R`$ est
```math
H=J\sum_{(i,j)}\Big[\tfrac12\big(s^+_is^-_j+s^-_is^+_j\big)+\Delta\,s^z_is^z_j\Big]+\sum_ih_is^z_i ,
```
la somme portant sur les liaisons ; $`\frac12(s^+_is^-_j+s^-_is^+_j)=s^x_is^x_j+s^y_is^y_j`$ déplace un spin vers le
haut le long de la liaison. Pour $`\Delta=1`$, c'est la chaîne de Heisenberg $`J\sum\mathbf s_i\cdot\mathbf s_j`$ ; pour
$`\Delta=0`$, la chaîne XX.

<a id="r-9-1"></a>

**Proposition 9.1.** $`[H,S^z]=0`$. Par conséquent $`H\in U(1)'\cong\bigoplus_{N=0}^LM_{\binom LN}(\mathbb C)`$, le bloc
$`N`$ agissant sur l'espace engendré par les configurations à $`N`$ spins vers le haut.

*Démonstration.* Les termes diagonaux commutent avec la matrice diagonale $`S^z`$. Pour $`i\neq j`$, $`s^+_is^-_j`$ commute
avec $`s^z_l`$, $`l\ne i,j`$, et $`[s^z_i+s^z_j,s^+_is^-_j]=[s^z_i,s^+_i]s^-_j+s^+_i[s^z_j,s^-_j]=s^+_is^-_j-s^+_is^-_j=0`$
d'après $`[s^z,s^\pm]=\pm s^\pm`$ (section 4.5) ; pour $`s^-_is^+_j`$, on passe aux adjoints. Par conséquent $`H`$
commute avec $`N=S^z+\frac L2`$, qui a les valeurs propres entières $`0,\dots,L`$, de sous-espaces propres engendrés par
les configurations à $`N`$ spins vers le haut, de dimension $`\binom LN`$ ; le reste est la proposition 4.10(2) avec
$`Q=N`$. $`\square`$

*Calcul :* la proposition 9.1, et les blocs $`H_N`$ que `SpinChain.xxz` construit à partir des configurations de chaque secteur (`chain.labels[N]`), comparés à la matrice complète $`2^L\times2^L`$ — [tutoriel, section 9.1](../../tutorial/fr/09_spin_chains.md#s-9-1).

**La tour des chaînes.** Soit $`\mathcal A_L=\bigoplus_NM_{\binom LN}`$, agissant sur $`L`$ sites. Ajouter un site,
$`x\mapsto x\otimes1_2`$, envoie $`\mathcal A_L`$ dans $`\mathcal A_{L+1}`$, car $`x\otimes1`$ commute avec
$`S^z\otimes1+1\otimes s^z`$. On obtient une tour $`\mathcal A_0\subset\mathcal A_1\subset\mathcal A_2\subset\cdots`$.

<a id="r-9-2"></a>

**Proposition 9.2.** L'inclusion $`\mathcal A_L\subset\mathcal A_{L+1}`$ a pour matrice d'inclusion la matrice
$`\Lambda_L`$ de taille $`(L+2)\times(L+1)`$ telle que $`\lambda_{N,N}=\lambda_{N,N-1}=1`$, tous les autres coefficients étant
nuls ; les tailles des blocs obéissent à la règle de Pascal $`\binom{L+1}N=\binom LN+\binom L{N-1}`$, le diagramme de
Bratteli de la tour est le triangle de Pascal, et le graphe biparti de $`\Lambda_L`$ est un chemin à $`2L+3`$ sommets.

*Démonstration.* Les configurations à $`N`$ spins vers le haut sur $`L+1`$ sites sont celles qui ont $`N`$ spins vers le
haut sur les $`L`$ premiers sites et le nouveau spin vers le bas, et celles qui en ont $`N-1`$ avec le nouveau spin vers
le haut. Sur ces deux sous-espaces, $`x\otimes1`$ agit comme $`x_N`$ et $`x_{N-1}`$, donc le bloc $`N`$ de $`\mathcal A_{L+1}`$
contient une fois chacun les blocs $`N`$ et $`N-1`$ de $`\mathcal A_L`$ ; les dimensions s'additionnent selon la règle de
Pascal. Dans le graphe biparti, le secteur $`N`$ de $`\mathcal A_{L+1}`$ est relié aux secteurs $`N-1`$ et $`N`$ de
$`\mathcal A_L`$, ce qui donne le chemin $`0'-0-1'-1-\cdots-L'-L-(L+1)'`$ (les primes désignent $`\mathcal A_{L+1}`$). $`\square`$

Un chemin allant du sommet du diagramme de Bratteli au secteur $`N`$ du niveau $`L`$ est une configuration $`b`$ ; il y
en a $`\binom LN`$, la taille du bloc. Pour une inclusion connexe de matrice d'inclusion $`\Lambda`$, nous **définissons**
l'**indice** comme $`\|\Lambda\|^2`$ (norme d'opérateur) ; pour $`M_k\otimes1\subset M_k\otimes M_m`$, la matrice
d'inclusion est la matrice $`1\times1`$ $`(m)`$, donc l'indice est le $`m^2`$ de la section 7.2. `Inclusion.index` le
calcule. (Le nom est emprunté à la théorie de l'indice de V. Jones, *Invent. Math.* 72 (1983) 1 ; dans ce cours,
l'indice d'une inclusion d'algèbres de dimension finie est cette définition et rien de plus.)

*Dans la bibliothèque :* l'indice est `Inclusion.index` — [tutoriel, section 9.1](../../tutorial/fr/09_spin_chains.md#s-9-1).

<a id="r-9-3"></a>

**Lemme 9.3.** La matrice d'adjacence $`A_n`$ du chemin de sommets $`1,\dots,n`$ a pour valeurs propres
$`2\cos\frac{\pi k}{n+1}`$, $`k=1,\dots,n`$, de vecteurs propres $`v_k(j)=\sin\frac{\pi jk}{n+1}`$ ; par conséquent
$`\|A_n\|=2\cos\frac\pi{n+1}`$.

*Démonstration.* Posons $`v_k(0)=v_k(n+1)=0`$. Alors $`(A_nv_k)(j)=v_k(j-1)+v_k(j+1)=2\cos\frac{\pi k}{n+1}\,v_k(j)`$
d'après $`\sin(a-b)+\sin(a+b)=2\sin a\cos b`$, et $`v_k\ne0`$ car $`v_k(1)=\sin\frac{\pi k}{n+1}\neq0`$. Les $`n`$ valeurs
propres sont distinctes ($`\cos`$ est strictement décroissant sur $`(0,\pi)`$), ce sont donc toutes les valeurs
propres, et la norme de la matrice symétrique $`A_n`$ est le plus grand module, atteint en $`k=1`$ et $`k=n`$. $`\square`$

<a id="r-9-4"></a>

**Corollaire 9.4.** $`[\mathcal A_{L+1}:\mathcal A_L]=\|\Lambda_L\|^2=4\cos^2\frac\pi{2L+4}`$, qui croît vers $`4`$.

*Démonstration.* En ordonnant les sommets le long du chemin de la proposition 9.2, la matrice d'adjacence
$`B=\begin{pmatrix}0&\Lambda_L\\\Lambda_L^T&0\end{pmatrix}`$ devient $`A_{2L+3}`$ (une similitude par permutation). D'après le
lemme 1.1, $`\|B\|^2=\|B^2\|`$ et $`\|\Lambda_L\Lambda_L^T\|=\|\Lambda_L^T\Lambda_L\|=\|\Lambda_L\|^2`$, et $`B^2=\Lambda_L\Lambda_L^T\oplus\Lambda_L^T\Lambda_L`$, donc
$`\|B\|=\|\Lambda_L\|`$ ; le lemme 9.3 avec $`n=2L+3`$ donne la valeur. $`\square`$

*Calcul :* la proposition 9.2, le lemme 9.3 et le corollaire 9.4 pour $`L\le8`$, avec `inclusions.chain_tower` et `inclusions.bratteli_diagram` — [tutoriel, section 9.1](../../tutorial/fr/09_spin_chains.md#s-9-1).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). $`S^z`$ est conservé dès que l'interaction est invariante par les rotations autour
> de l'axe $`z`$ ; $`N`$ compte les spins vers le haut, ou les particules après la transformation de Jordan–Wigner de
> la section 9.4. Dans la tour, chaque secteur de charge de la chaîne la plus longue est alimenté par deux
> secteurs de la plus courte.

<a id="s-9-2"></a>

## 9.2 Translations, impulsion et règle des signes de Marshall

Sur une chaîne périodique, soit $`T`$ la translation $`T|b_0b_1\dots b_{L-1}\rangle=|b_{L-1}b_0\dots b_{L-2}\rangle`$, qui
envoie le site $`i`$ sur le site $`i+1`$. $`T`$ est une matrice de permutation, donc unitaire, $`T^L=1`$ et
$`Ts^a_iT^*=s^a_{i+1}`$. Par conséquent $`TS^zT^*=S^z`$ et $`THT^*=H`$ pour un champ uniforme $`h_i=h`$ : $`T`$ envoie les
liaisons sur des liaisons. Ses valeurs propres sont donc de la forme $`e^{ik}`$, d'**impulsions** $`k=2\pi m/L`$, $`m=0,\dots,L-1`$.

<a id="r-9-5"></a>

**Proposition 9.5 (secteurs d'impulsion).** Fixons $`N`$. Pour une configuration $`a`$ à $`N`$ spins vers le haut, soit
$`R_a`$ la taille de son orbite $`\{T^ja\}`$. (a) $`R_a`$ divise $`L`$ ; pour $`k=2\pi m/L`$ avec $`mR_a\equiv0\pmod L`$, le
vecteur $`|a,k\rangle=\sum_{j=0}^{R_a-1}e^{-ikj}T^j|a\rangle`$ est non nul et $`T|a,k\rangle=e^{ik}|a,k\rangle`$, et ces $`R_a`$
vecteurs forment une base orthogonale de l'espace engendré par l'orbite. (b) La valeur propre $`e^{2\pi im/L}`$ de
$`T`$ sur le secteur $`N`$ a la multiplicité $`d_{N,m}=\#\{\text{orbites }:mR_a\equiv0\pmod L\}`$, et tout opérateur qui
commute avec $`S^z`$ et $`T`$ — par exemple la chaîne XXZ invariante par translation — est diagonal par blocs, dans
$`\bigoplus_N\bigoplus_mM_{d_{N,m}}`$.

*Démonstration.* (a) Le stabilisateur de $`a`$ dans le groupe cyclique $`\{T^j\}\cong\mathbb Z_L`$ est un sous-groupe, et
$`R_a`$ en est l'indice, un diviseur de $`L`$ ; les vecteurs $`T^j|a\rangle`$, $`0\le j<R_a`$, sont orthonormés et
$`T^{R_a}|a\rangle=|a\rangle`$. Si $`e^{ikR_a}=1`$, alors $`T|a,k\rangle=e^{ik}\sum_{j=1}^{R_a}e^{-ikj}T^j|a\rangle=e^{ik}|a,k\rangle`$,
puisque le terme $`j=R_a`$ est égal au terme $`j=0`$ ; le coefficient de $`|a\rangle`$ vaut $`1`$. Les
$`m\in\{0,\dots,L-1\}`$ admissibles sont les multiples de $`L/R_a`$, au nombre de $`R_a`$, avec des valeurs propres
distinctes de l'unitaire $`T`$ ; les vecteurs sont donc orthogonaux et engendrent l'espace de dimension $`R_a`$
engendré par l'orbite. (b) Le secteur $`N`$ est la somme orthogonale des espaces engendrés par ses orbites,
chacun invariant par $`T`$, et d'après (a) la valeur propre $`e^{2\pi im/L}`$ apparaît une fois dans l'espace engendré
par l'orbite de $`a`$ si $`mR_a\equiv0\pmod L`$, et pas du tout sinon. Un $`X`$
qui commute avec $`S^z`$ et $`T`$ commute avec leurs projecteurs spectraux $`P_N`$, $`Q_m`$ (théorème 1.5(d)), donc
$`X=\sum_{N,m}P_NQ_mXP_NQ_m`$. $`\square`$

*Calcul :* la proposition 9.5(b) pour $`L=12`$, et les blocs d'impulsion construits par `SpinChain.xxz_momentum` — [tutoriel, section 9.2](../../tutorial/fr/09_spin_chains.md#s-9-2).

<a id="r-9-6"></a>

**Lemme 9.6 (Perron–Frobenius, cas symétrique).** Soit $`M`$ une matrice réelle symétrique avec $`M_{ij}\ge0`$ pour
$`i\ne j`$, dont le graphe ($`i\sim j`$ si et seulement si $`M_{ij}>0`$) est connexe. Alors sa plus grande valeur propre
$`\lambda`$ est simple, avec un vecteur propre dont tous les coefficients sont $`>0`$.

*Démonstration.* $`\lambda=\max_{\|x\|=1}\langle x,Mx\rangle`$, et le maximum est atteint exactement sur les vecteurs
unitaires du sous-espace propre de $`\lambda`$ (développer dans une base orthonormée de vecteurs propres). Pour un tel
$`x`$, le vecteur $`|x|`$ des modules vérifie $`\langle|x|,M|x|\rangle\ge\mathrm{Re}\langle x,Mx\rangle=\lambda`$, car les
termes diagonaux coïncident et $`M_{ij}|x_i||x_j|\ge\mathrm{Re}\,M_{ij}\bar x_ix_j`$ pour $`i\ne j`$ ; donc $`|x|`$ est un
vecteur propre. Si $`|x|_i=0`$, alors $`0=(\lambda-M_{ii})|x|_i=\sum_{j\ne i}M_{ij}|x|_j`$ impose $`|x|_j=0`$ pour tous les
voisins $`j`$ de $`i`$, et par connexité $`x=0`$. Aucun vecteur propre non nul n'a donc de coefficient nul ; or un
sous-espace propre de dimension $`\ge2`$ contient un vecteur non nul avec $`x_1=0`$. Par conséquent $`\lambda`$ est simple,
et $`|x|>0`$ est un vecteur propre. $`\square`$

Nous démontrons la règle des signes pour toute disposition bipartie des couplages ; la chaîne en est un cas
particulier.

<a id="r-9-7"></a>

**Théorème 9.7 (règle des signes de Marshall).** Soit $`A\sqcup B`$ une partition des sites $`\{0,\dots,L-1\}`$ en deux
ensembles, et soit
```math
H=\sum_{i<j}\Big[\tfrac{J_{ij}}2\big(s^+_is^-_j+s^-_is^+_j\big)+K_{ij}\,s^z_is^z_j\Big]+\sum_ih_is^z_i
```
avec $`K_{ij},h_i`$ réels et des $`J_{ij}`$ réels tels que $`J_{ij}\ge0`$ si $`i,j`$ sont dans des ensembles différents,
$`J_{ij}\le0`$ s'ils sont dans le même ensemble, et que le graphe $`G`$ sur les sites, avec une arête $`i\sim j`$ si et
seulement si $`J_{ij}\ne0`$, soit connexe. Soit $`N_B(b)`$ le nombre de spins vers le haut de $`b`$ sur $`B`$. Dans chaque
secteur $`N`$, la plus petite valeur propre de $`H_N`$ est simple, et son vecteur propre peut être choisi tel que
$`(-1)^{N_B(b)}\psi(b)>0`$ pour toute configuration $`b`$ à $`N`$ spins vers le haut.

En particulier, cela vaut pour la chaîne XXZ avec $`J>0`$ et $`\Delta,h_i`$ quelconques si la chaîne est ouverte, ou
périodique avec $`L`$ pair, avec pour $`A`$ les sites pairs et pour $`B`$ les sites impairs : alors $`J_{ij}=J`$,
$`K_{ij}=J\Delta`$ sur les liaisons et $`0`$ sinon, toute liaison relie un site pair et un site impair (pour l'anneau,
parce que $`L`$ est pair), et $`G`$ contient le chemin $`0-1-\cdots-(L-1)`$.

*Démonstration.* Soit $`V|b\rangle=(-1)^{N_B(b)}|b\rangle`$ ; $`V`$ est diagonal et $`V=V^*=V^{-1}`$. L'opérateur $`s^+_is^-_j`$
déplace un spin vers le haut de $`j`$ à $`i`$ ; il change $`N_B`$ de $`\pm1`$ si $`i,j`$ sont dans des ensembles différents
et le laisse inchangé sinon, donc $`Vs^+_is^-_jV=\mp s^+_is^-_j`$ selon le cas, et de même pour $`s^-_is^+_j`$. Deux
configurations distinctes $`b,b'`$ sont reliées par un tel terme si et seulement si $`b'`$ s'obtient à partir de $`b`$ en
échangeant les spins opposés de deux sites $`i,j`$, et alors $`i,j`$ sont déterminés par $`b,b'`$ (les sites où elles
diffèrent). Ainsi $`M=-VH_NV`$ a pour coefficients hors diagonale $`J_{ij}/2`$ (ensembles différents) ou $`-J_{ij}/2`$
(même ensemble), c'est-à-dire $`|J_{ij}|/2\ge0`$ dans les deux cas, strictement positifs si et seulement si $`i\sim j`$
dans $`G`$, et $`0`$ entre configurations qui diffèrent autrement. Le graphe de $`M`$ est connexe. En effet, notons
$`(i\,j)`$ l'échange des spins des sites $`i,j`$ ; pour $`a,b,c`$ distincts, $`(a\,c)=(a\,b)(b\,c)(a\,b)`$, donc, par
récurrence le long d'un chemin de $`G`$ de $`i`$ à $`j`$, tout $`(i\,j)`$ est une composée d'échanges le long d'arêtes de
$`G`$, et chacun de ces échanges ou bien fixe une configuration (spins égaux), ou bien suit une arête du graphe de
$`M`$. Si $`b\ne b'`$ ont $`N`$ spins vers le haut, il existe des sites $`i`$, vers le haut dans $`b`$ et vers le bas dans
$`b'`$, et $`j`$, vers le bas dans $`b`$ et vers le haut dans $`b'`$ ; $`(i\,j)b`$ diffère de $`b'`$ en deux sites de moins que
$`b`$, donc un nombre fini de telles étapes mène de $`b`$ à $`b'`$. D'après le lemme 9.6, la plus grande valeur propre de
$`M`$ — l'opposée de la plus petite de $`H_N`$ — est simple, avec un vecteur propre $`\varphi>0`$ ; on prend
$`\psi=V\varphi`$. $`\square`$

<a id="r-9-8"></a>

**Corollaire 9.8 (impulsion de l'état fondamental).** Sur une chaîne périodique avec $`L`$ pair, $`J>0`$ et un champ
uniforme, l'état fondamental du secteur $`N`$ a l'impulsion $`k=0`$ si $`N`$ est pair et $`k=\pi`$ si $`N`$ est impair.

*Démonstration.* $`T`$ commute avec $`H_N`$ et l'état fondamental $`\psi=V\varphi`$ est simple, donc $`T\psi=e^{ik}\psi`$. Or
$`(T\psi)(b)=\psi(T^{-1}b)`$, et $`T^{-1}b`$ porte sur le site $`i`$ le spin $`b_{i+1}`$ ; comme $`L`$ est pair, $`i+1`$ est pair
pour $`i`$ impair, donc $`N_B(T^{-1}b)=N-N_B(b)`$ et
$`(VT\psi)(b)=(-1)^{N_B(b)}(-1)^{N-N_B(b)}\varphi(T^{-1}b)=(-1)^N\varphi(T^{-1}b)`$. Par ailleurs $`VT\psi=e^{ik}\varphi`$.
$`\varphi`$ et $`\varphi\circ T^{-1}`$ sont tous deux positifs, donc $`e^{ik}=(-1)^N`$. $`\square`$

*Calcul :* le théorème 9.7 pour des chaînes ouvertes et périodiques en champs aléatoires, et le corollaire 9.8 pour $`L=12`$ — [tutoriel, section 9.2](../../tutorial/fr/09_spin_chains.md#s-9-2).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Dans un
> antiferromagnétique, le réseau fixe les signes de l'état fondamental : déplacer un spin à travers une liaison
> change le signe. Les physiciens s'en servent pour faire tourner le Monte-Carlo quantique sans problème du signe
> sur les antiferromagnétiques bipartis ; pour les chaînes frustrées, où un couplage antiferromagnétique relie
> deux sites du même ensemble, l'hypothèse du théorème 9.7 n'est pas satisfaite.

<a id="s-9-3"></a>

## 9.3 $`SU(2)`$ en $`\Delta=1`$

Pour $`\Delta=1`$ et $`h=0`$, $`H=J\sum\mathbf s_i\cdot\mathbf s_j`$ commute avec $`S^x,S^y,S^z`$ (proposition 4.15), donc
$`H\cong\bigoplus_S1_{2S+1}\otimes H_S`$ avec des blocs $`H_S\in M_{m_S}`$, $`m_S=\binom L{L/2-S}-\binom L{L/2-S-1}`$
(proposition 4.13), et le module `su2` construit $`H_S`$ sur les vecteurs de plus haut poids (proposition 4.15(b)).
Il en va de même pour tout $`H=\sum_{i<j}J_{ij}\,\mathbf s_i\cdot\mathbf s_j`$ avec des $`J_{ij}`$ réels : d'après la
proposition 4.15(a), chaque $`\mathbf s_i\cdot\mathbf s_j`$ appartient à $`\{S_a\}'`$. Notons $`E_0(S)`$ la plus petite valeur propre
de $`H_S`$, pour les spins $`S`$ qui apparaissent.

*Dans la bibliothèque :* les blocs $`H_S`$ sont construits par `SU2Chain` — [tutoriel, section 9.3](../../tutorial/fr/09_spin_chains.md#s-9-3).

<a id="r-9-9"></a>

**Théorème 9.9 (Lieb–Mattis).** Soit $`H=\sum_{i<j}J_{ij}\,\mathbf s_i\cdot\mathbf s_j`$ avec $`A,B,J_{ij}`$ comme dans le
théorème 9.7 (donc $`K_{ij}=J_{ij}`$, $`h_i=0`$), $`n_A=|A|\ge n_B=|B|`$ et $`S_*=\frac{n_A-n_B}2`$.
(a) Pour $`S_*\le S\le\frac L2`$ (avec $`S-\frac L2`$ entier), l'état fondamental du secteur $`S^z=S`$ a le spin total $`S`$ ;
pour $`0\le M\le S_*`$ (avec $`M-\frac L2`$ entier), l'état fondamental du secteur $`S^z=M`$ a le spin total $`S_*`$.
(b) Par conséquent $`E_0(S)`$ est une valeur propre simple de $`H_S`$ pour $`S\ge S_*`$,
$`E_0(S_*)<E_0(S_*+1)<\dots<E_0(\frac L2)`$, et $`E_0(S)>E_0(S_*)`$ pour $`S<S_*`$ : l'espace fondamental de $`H`$ est un unique
multiplet, de spin $`S_*`$.
Pour la chaîne de Heisenberg avec $`J>0`$, ouverte ou périodique avec $`L`$ pair, $`A`$ et $`B`$ sont les sites pairs et
impairs, $`S_*=0`$ pour $`L`$ pair et $`S_*=\frac12`$ pour $`L`$ impair (chaîne ouverte), et
$`E_0(S_*)<E_0(S_*+1)<\dots<E_0(\frac L2)`$.

*Démonstration.* (a) Soient $`M\ge0`$ avec $`M-\frac L2\in\mathbb Z`$, $`s=\max(M,S_*)`$, et $`\psi`$ l'état fondamental du
secteur $`S^z=M`$, simple d'après le théorème 9.7 : $`\psi=V\varphi`$ avec $`\varphi>0`$. $`\mathbf S^2`$ commute avec $`H`$ et avec
$`S^z`$, donc il préserve l'espace fondamental de dimension un du secteur : $`\psi`$ a un spin total $`S'`$. Nous
construisons un vecteur $`\chi`$ du secteur avec $`\mathbf S^2\chi=s(s+1)\chi`$ et $`\langle\psi,\chi\rangle\ne0`$ ; comme des
vecteurs propres de $`\mathbf S^2`$ pour des valeurs propres distinctes sont orthogonaux, cela donne $`S'=s`$. Posons
$`p=\frac L2-s`$ ; alors $`0\le p\le n_B`$, car $`s\ge S_*`$ signifie $`p\le\frac{n_A+n_B}2-\frac{n_A-n_B}2`$. Choisissons des
$`a_1,\dots,a_p\in A`$ distincts et des $`b_1,\dots,b_p\in B`$ distincts, soit $`R`$ l'ensemble des $`2s`$ sites restants, et
soit
```math
\chi_0=\bigotimes_{l=1}^p\tfrac1{\sqrt2}\big(|{\uparrow}_{a_l}{\downarrow}_{b_l}\rangle-|{\downarrow}_{a_l}{\uparrow}_{b_l}\rangle\big)\otimes|{\uparrow\cdots\uparrow}\rangle_R ,
\qquad \chi=(S^-)^{s-M}\chi_0 .
```
Chaque singulet est annulé par $`s^\pm_{a_l}+s^\pm_{b_l}`$ et les spins vers le haut de $`R`$ par $`s^+_i`$, donc
$`S^+\chi_0=0`$ ; avec $`S^z\chi_0=s\chi_0`$, cela donne $`\mathbf S^2\chi_0=(S^-S^++(S^z)^2+S^z)\chi_0=s(s+1)\chi_0`$
(section 4.5), et $`\mathbf S^2\chi=s(s+1)\chi`$ puisque $`\mathbf S^2`$ commute avec $`S^-`$. En écrivant
$`S^-=\sum_l(s^-_{a_l}+s^-_{b_l})+S^-_R`$, les termes commutent et les premiers annulent les singulets, donc $`\chi`$ est le
produit des singulets par $`(S^-_R)^{s-M}|{\uparrow\cdots\uparrow}\rangle_R=(s-M)!\sum_D|D\rangle`$, la somme portant sur les
configurations de $`R`$ ayant exactement $`s-M`$ spins vers le bas (il en existe, car $`s-M\le2s=|R|`$) ; ainsi $`\chi\ne0`$,
$`S^z\chi=M\chi`$, et les coefficients de $`\chi`$ sont réels. Son signe : un singulet contribue $`+`$ quand son spin vers
le haut est sur $`A`$ et $`-`$ quand il est sur $`B`$, et la partie sur $`R`$ est positive. Si $`M=s`$, les sites de $`R`$ sont
tous vers le haut ; si $`M<s`$, alors $`s=S_*`$, $`p=n_B`$ et $`R\subseteq A`$. Dans les deux cas, $`N_B(b)`$ est le nombre de
singulets dont le spin vers le haut est sur $`B`$, plus une constante, $`|R\cap B|`$ (si $`M=s`$) ou $`0`$, donc
$`(-1)^{N_B(b)}\chi(b)`$ a le même signe pour tous les $`b`$ tels que $`\chi(b)\ne0`$. Par conséquent
$`\langle\psi,\chi\rangle=\langle\varphi,V\chi\rangle=\pm\sum_b\varphi(b)|\chi(b)|\ne0`$.

(b) D'après la proposition 4.15(a), le corollaire 4.12 et le théorème 4.11, un multiplet de spin $`S'`$ rencontre le
secteur $`S^z=M`$ en un vecteur si $`S'\ge|M|`$ et pas du tout sinon, donc le spectre de $`H`$ sur le secteur $`S^z=M\ge0`$
est la réunion des spectres des $`H_{S'}`$, $`S'\ge M`$, avec multiplicités. D'après (a), sa plus petite valeur propre
est simple et appartient au spin $`s=\max(M,S_*)`$ : ainsi $`E_0(s)`$ est simple dans $`H_s`$ et $`E_0(s)<E_0(S')`$ pour tout
$`S'\ge M`$, $`S'\ne s`$. Pour $`M=S\ge S_*`$, cela donne la simplicité et $`E_0(S)<E_0(S')`$ pour $`S'>S`$ ; pour $`M=S<S_*`$,
cela donne $`E_0(S_*)<E_0(S)`$. Ainsi $`E_0(S_*)`$ est la plus petite valeur propre de $`H`$, atteinte seulement dans
$`H_{S_*}`$, où elle est simple. Pour la chaîne, les sites pairs sont plus nombreux que les impairs de $`0`$ ($`L`$ pair)
ou de $`1`$ ($`L`$ impair). $`\square`$

Le théorème 9.7 remonte à W. Marshall (*Proc. Roy. Soc. A* 232 (1955) 48) et le théorème 9.9 à E. Lieb et
D. Mattis (*J. Math. Phys.* 3 (1962) 749). Le **gap singulet–triplet** $`E_0(1)-E_0(0)`$ de l'anneau de Heisenberg
pair est positif d'après le théorème 9.9.

*Calcul :* le théorème 9.9 à l'aide des secteurs $`U(1)`$ pour $`L=12`$, le gap singulet–triplet, et les théorèmes 9.7 et 9.9 sur un graphe biparti qui n'est pas une chaîne — [tutoriel, section 9.3](../../tutorial/fr/09_spin_chains.md#s-9-3).

**Corrélations.** Sur tout anneau de Heisenberg pair, les corrélations de l'état fondamental
$`\langle\mathbf s_0\cdot\mathbf s_r\rangle`$ alternent en signe ; les théorèmes 9.7 et 9.9 l'imposent. Soit $`\psi`$ l'état fondamental
réel normalisé, un singulet (théorème 9.9), donc l'état fondamental du secteur $`S^z=0`$. De $`\mathbf S^2\psi=0`$ on tire
$`\|S^\pm\psi\|^2=\langle\psi,S^\mp S^\pm\psi\rangle=\langle\psi,(\mathbf S^2-(S^z)^2\mp S^z)\psi\rangle=0`$, donc $`S^\pm\psi=0`$ et
$`\langle\psi,[S^+,[S^-,X]]\psi\rangle=0`$ pour tout $`X`$ (chacun des quatre termes a $`S^\pm`$ à côté de $`\psi`$). Pour
$`X=s^z_0s^z_r`$, $`r\ne0`$, les relations $`[s^\pm,s^z]=\mp s^\pm`$, $`[s^+,s^-]=2s^z`$ donnent $`[S^-,X]=s^-_0s^z_r+s^z_0s^-_r`$ et
$`[S^+,[S^-,X]]=4s^z_0s^z_r-s^+_0s^-_r-s^-_0s^+_r`$. Par conséquent
$`\langle s^z_0s^z_r\rangle=\frac14\langle s^+_0s^-_r+s^-_0s^+_r\rangle`$ et
```math
\langle\mathbf s_0\cdot\mathbf s_r\rangle=\langle s^z_0s^z_r\rangle+\tfrac12\langle s^+_0s^-_r+s^-_0s^+_r\rangle
=\tfrac34\sum_{b_0\ne b_r}\psi(b)\psi(b') ,
```
où $`b'`$ s'obtient à partir de $`b`$ en échangeant les spins des sites $`0`$ et $`r`$. D'après le théorème 9.7,
$`\psi(b)\psi(b')`$ a le signe de $`(-1)^{N_B(b)+N_B(b')}`$, et $`N_B(b')-N_B(b)`$ vaut $`\pm1`$ pour $`r`$ impair et $`0`$ pour $`r`$
pair ; la somme n'est pas vide, donc $`\langle\mathbf s_0\cdot\mathbf s_r\rangle`$ est non nul, de signe $`(-1)^r`$, pour
$`0<r<L`$.

*Calcul :* les corrélations $`\langle\mathbf s_0\cdot\mathbf s_r\rangle`$ pour $`L=12`$, avec `SU2Chain.spin_correlation`, qui utilise le théorème de Wigner–Eckart (théorème B.16) — [tutoriel, section 9.3](../../tutorial/fr/09_spin_chains.md#s-9-3).

<a id="s-9-4"></a>

## 9.4 Vérifications exactes : fermions libres et ansatz de Bethe

**Jordan–Wigner.** Sur une chaîne ouverte, posons $`c_j=\big(\prod_{l<j}(-1)^{n_l}\big)s^-_j`$.

<a id="r-9-10"></a>

**Théorème 9.10 (la chaîne XX est libre).** (a) $`\{c_i,c_j^*\}=\delta_{ij}`$, $`\{c_i,c_j\}=0`$ et $`c_j^*c_j=n_j`$.
(b) La chaîne XX ouverte s'écrit $`H=\frac J2\sum_{j<L-1}(s^+_js^-_{j+1}+s^-_js^+_{j+1})=\sum_{i,j}h_{ij}c_i^*c_j`$ avec
$`h=\frac J2A_L`$. (c) Soient $`u_k`$ des vecteurs propres réels orthonormés de $`h`$, de valeurs propres
$`\varepsilon_k=J\cos\frac{\pi k}{L+1}`$ (lemme 9.3), et $`b_k=\sum_ju_k(j)c_j`$. Pour $`K=\{k_1<\dots<k_N\}`$, les vecteurs
$`b_K^*|0\rangle=b_{k_1}^*\cdots b_{k_N}^*|0\rangle`$ ($`|0\rangle`$ : tous les spins vers le bas) forment une base orthonormée
de vecteurs propres du secteur $`N`$, d'énergies $`\sum_{k\in K}\varepsilon_k`$.

*Démonstration.* (a) Des opérateurs sur des sites différents commutent, et sur un même site $`(-1)^n`$ anticommute
avec $`s^\pm`$. Pour $`i<j`$, $`c_i`$ et $`c_j^*`$ (ou $`c_j`$) se rencontrent sur le site $`i`$ sous la forme $`s^-_i`$ et
$`(-1)^{n_i}`$ — une anticommutation — et commutent par ailleurs, donc ils anticommutent. Pour $`i=j`$, les chaînes de
signes sont de carré $`1`$ et $`c_jc_j^*+c_j^*c_j=s^-_js^+_j+s^+_js^-_j=1`$, $`c_j^*c_j=s^+_js^-_j=n_j`$, $`c_j^2=(s^-_j)^2=0`$ ; le cas $`i>j`$ s'en déduit
en échangeant les rôles de $`i,j`$ (et en passant aux adjoints). (b)
$`c_j^*c_{j+1}=s^+_j(-1)^{n_j}s^-_{j+1}=s^+_js^-_{j+1}`$, car $`s^+_j\neq0`$ seulement sur $`n_j=0`$ ; on passe aux adjoints.
(c) $`\{b_k,b_l^*\}=\sum_ju_k(j)u_l(j)=\delta_{kl}`$ et $`\{b_k,b_l\}=0`$, et $`H=\sum_k\varepsilon_kb_k^*b_k`$ par la
décomposition spectrale de $`h`$. D'après les CAR, $`[b_l^*b_l,b_k^*]=\delta_{lk}b_k^*`$, donc $`[H,b_k^*]=\varepsilon_kb_k^*`$ ;
avec $`H|0\rangle=0`$, cela donne $`Hb_K^*|0\rangle=(\sum_K\varepsilon_k)b_K^*|0\rangle`$. Comme $`b_k|0\rangle=0`$ (chaque $`c_j`$
abaisse un spin), faire passer les annihilateurs à droite à l'aide des CAR donne
$`\langle b_{K'}^*0,b_K^*0\rangle=\delta_{KK'}`$ par récurrence sur $`N`$. Chaque $`c_j^*`$ augmente $`N`$ de un, donc
$`b_K^*|0\rangle`$ appartient au secteur $`|K|`$ : $`\binom LN`$ vecteurs orthonormés dans un espace de cette dimension. $`\square`$

*Calcul :* les $`1024`$ niveaux de la chaîne XX ouverte avec $`L=10`$, comparés au théorème 9.10 — [tutoriel, section 9.4](../../tutorial/fr/09_spin_chains.md#s-9-4).

**Magnons.** Sur l'anneau avec $`h=0`$, l'état complètement polarisé $`|0\rangle`$ a l'énergie $`E_F=J\Delta L/4`$. Un spin
vers le haut en $`x`$ brise deux liaisons : $`H|x\rangle=(E_F-J\Delta)|x\rangle+\frac J2(|x-1\rangle+|x+1\rangle)`$ ($`L\ge3`$).

<a id="r-9-11"></a>

**Proposition 9.11 (un magnon).** Dans le secteur $`N=1`$ de l'anneau, $`|k\rangle=\sum_xe^{-ikx}|x\rangle`$,
$`k=2\pi m/L`$, vérifie $`T|k\rangle=e^{ik}|k\rangle`$ et $`H|k\rangle=\big(E_F+J(\cos k-\Delta)\big)|k\rangle`$.

*Démonstration.* $`\sum_xe^{-ikx}(|x-1\rangle+|x+1\rangle)=\sum_y(e^{-ik(y+1)}+e^{-ik(y-1)})|y\rangle=2\cos k\,|k\rangle`$, en
utilisant $`e^{ikL}=1`$ ; le même changement d'indice donne $`T|k\rangle=e^{ik}|k\rangle`$. $`\square`$

Deux magnons diffusent l'un sur l'autre. Notons $`|x_1,x_2\rangle`$, $`x_1<x_2`$, l'état aux spins vers le haut en
$`x_1,x_2`$ ; une paire de voisins brise deux liaisons, une paire séparée quatre.

<a id="r-9-12"></a>

**Proposition 9.12 (deux magnons ; équations de Bethe).** Soient $`k_1,k_2\in\mathbb C`$, $`K=k_1+k_2`$, et
$`(A_{12},A_{21})\neq0`$ tels que (i) $`A_{12}(1+e^{iK}-2\Delta e^{ik_2})+A_{21}(1+e^{iK}-2\Delta e^{ik_1})=0`$ et
(ii) $`A_{21}=e^{ik_2L}A_{12}`$, $`A_{12}=e^{ik_1L}A_{21}`$. Soient
$`\psi(x_1,x_2)=A_{12}e^{i(k_1x_1+k_2x_2)}+A_{21}e^{i(k_2x_1+k_1x_2)}`$ sur $`\mathbb Z^2`$ et
$`\Psi=\sum_{0\le x_1<x_2<L}\psi(x_1,x_2)|x_1,x_2\rangle`$. Si $`\Psi\neq0`$, alors
$`H\Psi=\big(E_F+J(\cos k_1+\cos k_2-2\Delta)\big)\Psi`$. Si les dénominateurs ne s'annulent pas, (i) et (ii)
s'écrivent
```math
e^{ik_1L}=-\frac{1+e^{iK}-2\Delta e^{ik_1}}{1+e^{iK}-2\Delta e^{ik_2}},\qquad
e^{ik_2L}=-\frac{1+e^{iK}-2\Delta e^{ik_2}}{1+e^{iK}-2\Delta e^{ik_1}} .
```

*Démonstration.* Soit $`E=E_F+J(\cos k_1+\cos k_2-2\Delta)`$. (1) Pour tout $`(x_1,x_2)\in\mathbb Z^2`$, chaque onde plane
donne des valeurs $`\psi(x_1\pm1,x_2)`$ et $`\psi(x_1,x_2\pm1)`$ de somme $`2(\cos k_1+\cos k_2)\psi(x_1,x_2)`$. (2) Par
substitution, $`\psi(x,x)+\psi(x+1,x+1)-2\Delta\psi(x,x+1)`$ vaut $`e^{iKx}`$ fois le membre de gauche de (i), donc
s'annule pour tout $`x`$. (3) D'après (ii),
$`\psi(x_2,x_1+L)=A_{12}e^{ik_2L}e^{i(k_2x_1+k_1x_2)}+A_{21}e^{ik_1L}e^{i(k_1x_1+k_2x_2)}=\psi(x_1,x_2)`$, et en appliquant
ceci deux fois, $`\psi(x_1+L,x_2+L)=\psi(x_1,x_2)`$. Les représentants d'une configuration $`\{y_1<y_2\}`$ dans
$`D=\{x_1<x_2<x_1+L\}`$ sont $`(y_1+nL,y_2+nL)`$ et $`(y_2+nL,y_1+(n+1)L)`$, donc $`\psi`$ est une fonction des
configurations. Soit $`(x_1,x_2)\in D`$. Si $`2\le x_2-x_1\le L-2`$, les quatre sauts restent dans $`D`$ et le terme
diagonal vaut $`E_F-2J\Delta`$, donc $`(H\Psi)(x_1,x_2)=(E_F-2J\Delta)\psi+\frac J2\cdot2(\cos k_1+\cos k_2)\psi=E\psi`$
d'après (1). Si $`x_2=x_1+1`$, les sauts $`(x_1+1,x_2)`$, $`(x_1,x_2-1)`$ sont bloqués et le terme diagonal vaut
$`E_F-J\Delta`$ ; d'après (1) et (2),
$`(H\Psi)(x_1,x_2)=(E_F-J\Delta)\psi+\frac J2\big(2(\cos k_1+\cos k_2)\psi-\psi(x_1,x_1)-\psi(x_1+1,x_1+1)\big)=E\psi`$.
Si $`x_2=x_1+L-1`$, on utilise le représentant $`(x_2,x_1+L)`$, qui est du type précédent ($`L\ge3`$). $`\square`$

La proposition 9.12 dit seulement que les solutions donnent des vecteurs propres ; elle ne dit pas que tous les
vecteurs propres s'obtiennent ainsi. Écrivons $`k_{1,2}=K/2\pm q`$ et $`c=\cos\frac K2`$. Alors
$`1+e^{iK}-2\Delta e^{ik_{1,2}}=2e^{iK/2}(c-\Delta e^{\pm iq})`$, donc $`A_{12}=c-\Delta e^{iq}`$,
$`A_{21}=-(c-\Delta e^{-iq})`$ résolvent (i). D'après (ii) et $`(A_{12},A_{21})\ne0`$, les deux coefficients sont non
nuls, donc le produit des deux équations (ii) donne $`e^{iKL}=1`$, $`K=2\pi n/L`$, et
(ii) équivaut alors à $`F(q)=(-1)^ne^{iqL}(c-\Delta e^{-iq})+c-\Delta e^{iq}=0`$, une équation à une variable. Nous
cherchons ses racines $`q\in(0,\pi)`$, où $`e^{-iqL/2}F(q)`$ est réel ou imaginaire pur parce que $`c-\Delta e^{\pm iq}`$
sont complexes conjugués (états de diffusion), et sur les droites $`q=i\kappa`$, $`\pi+i\kappa`$, $`\kappa>0`$, où $`F`$ est
réel (états liés).

Pour $`L`$ pair et $`\Delta\ne0`$, les vecteurs de Bethe ne sont pas complets. Si $`\Delta\ne0`$,
$`\Phi=\sum_x(-1)^x|x,x+1\rangle`$ vérifie $`T\Phi=-\Phi`$ et $`H\Phi=(E_F-J\Delta)\Phi`$ (les sauts vers $`|x-1,x+1\rangle`$ et
$`|x,x+2\rangle`$ se compensent deux à deux). Un vecteur de Bethe vérifie $`T\Psi=e^{-iK}\Psi`$, puisque
$`\psi(x_1-1,x_2-1)=e^{-iK}\psi(x_1,x_2)`$. Si $`e^{-iK}\ne-1`$, $`\Psi\perp\Phi`$ car ce sont des vecteurs propres de
l'unitaire $`T`$ pour des valeurs propres distinctes ; si $`e^{iK}=-1`$, alors
$`\cos k_1+\cos k_2=2\cos\frac K2\cos\frac{k_1-k_2}2=0`$, l'énergie de $`\Psi`$ vaut $`E_F-2J\Delta\ne E_F-J\Delta`$, et
$`\Psi\perp\Phi`$ en tant que vecteurs propres de $`H`$ pour des valeurs propres distinctes. Ainsi $`\Phi`$ est orthogonal à
tous les vecteurs de Bethe de la proposition 9.12, qui n'engendrent donc jamais le secteur $`N=2`$ d'un anneau pair.
Le calcul à deux magnons est le cas $`N=2`$ de l'ansatz
de H. Bethe (*Z. Phys.* 71 (1931) 205) ; nous n'utilisons aucun autre cas.

*Calcul :* les racines de $`F`$ pour $`L=9`$, $`\Delta=0.6`$, qui donnent les $`36`$ états du secteur $`N=2`$, et le vecteur $`\Phi`$ pour $`L=8`$ — [tutoriel, section 9.4](../../tutorial/fr/09_spin_chains.md#s-9-4).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Un magnon est une onde de
> spin : un spin retourné délocalisé sur l'anneau avec l'impulsion $`k`$. Deux magnons diffusent avec un déphasage
> ($`q`$ réel) ou forment un état lié ($`q`$ complexe).

<a id="s-9-5"></a>

## 9.5 Lanczos et méthodes de Krylov

Le secteur $`N=L/2`$ est de dimension $`12870`$ pour $`L=16`$ et $`184756`$ pour $`L=20`$, mais $`H`$ est creuse ($`O(L)`$
coefficients non nuls par ligne, `xxz_sparse`), et les méthodes de Krylov n'ont besoin que de produits $`Hv`$. Soient
$`H`$ hermitienne de taille $`n`$ et $`v\ne0`$. L'**espace de Krylov** est
$`\mathcal K_m=\mathrm{span}\{v,Hv,\dots,H^{m-1}v\}`$.

<a id="r-9-13"></a>

**Proposition 9.13 (Lanczos).** Supposons $`\dim\mathcal K_m=m`$. Soit $`q_1=v/\|v\|`$ et, pour $`j<m`$, soit $`q_{j+1}`$
la composante normalisée de $`Hq_j`$ orthogonale à $`\mathcal K_j`$. Alors $`q_1,\dots,q_m`$ est une base orthonormée de
$`\mathcal K_m`$ ; avec $`Q_m=(q_1\cdots q_m)`$, $`T_m=Q_m^*HQ_m`$ est réelle, symétrique et tridiagonale, de diagonale
$`\alpha_j=\langle q_j,Hq_j\rangle`$ et de sous-diagonale $`\beta_j=\langle q_{j+1},Hq_j\rangle>0`$ ; et $`HQ_m=Q_mT_m+r_me_m^T`$,
où $`r_m\perp\mathcal K_m`$ s'annule si et seulement si $`\mathcal K_{m+1}=\mathcal K_m`$ (sinon $`r_m=\beta_mq_{m+1}`$,
$`q_{m+1}=r_m/\|r_m\|`$).

*Démonstration.* Supposons $`\mathrm{span}\{q_1,\dots,q_j\}=\mathcal K_j`$ (vrai pour $`j=1`$). Alors
$`H\mathcal K_j\subseteq\mathcal K_j+\mathbb CHq_j`$, puisque $`Hq_i\in\mathcal K_{i+1}`$ ; donc
$`\mathcal K_{j+1}=\mathcal K_j+\mathbb CH^jv\subseteq\mathcal K_j+\mathbb CHq_j\subseteq\mathcal K_{j+1}`$, et
$`\mathcal K_{j+1}=\mathcal K_j+\mathbb CHq_j`$. Pour $`j<m`$ la dimension croît, donc la composante de $`Hq_j`$ orthogonale à
$`\mathcal K_j`$ est non nulle et la récurrence se poursuit ; pour $`j=m`$ cette composante, $`r_m`$, s'annule si et
seulement si $`\mathcal K_{m+1}=\mathcal K_m`$. Pour $`i>j+1`$, $`q_i\perp\mathcal K_{j+1}\ni Hq_j`$, donc $`(T_m)_{ij}=0`$, et
$`(T_m)_{ji}=0`$ puisque $`T_m`$ est hermitienne. $`\alpha_j`$ est réel, et $`\beta_j=\langle q_{j+1},Hq_j\rangle`$ est la norme de
la composante, $`>0`$. Ainsi $`Hq_j=\beta_{j-1}q_{j-1}+\alpha_jq_j+\beta_jq_{j+1}`$ pour $`j<m`$ et
$`Hq_m=\beta_{m-1}q_{m-1}+\alpha_mq_m+r_m`$. $`\square`$

Les valeurs propres $`\theta^{(m)}_1\le\dots\le\theta^{(m)}_m`$ de $`T_m`$ sont les **valeurs de Ritz**.

<a id="r-9-14"></a>

**Théorème 9.14 (convergence monotone).** Soient $`\lambda_1\le\dots\le\lambda_n`$ les valeurs propres de $`H`$ et
$`\dim\mathcal K_m=m`$. Alors $`\theta^{(m)}_j\ge\lambda_j`$, et $`\theta^{(m)}_j\ge\theta^{(m+1)}_j`$ si
$`\dim\mathcal K_{m+1}=m+1`$ ($`j\le m`$). Si $`v`$ a une composante non nulle dans le sous-espace propre de $`\lambda_1`$ et
si $`\mathcal K_{m+1}=\mathcal K_m`$, alors $`\theta^{(m)}_1=\lambda_1`$.

*Démonstration.* D'après le corollaire A.9 (Courant–Fischer pour les compressions, théorème A.8) appliqué à
$`Q=Q_m`$, $`\theta^{(m)}_j=\min_{U}\max_{x\in U,\|x\|=1}\langle x,Hx\rangle`$, le minimum portant sur les
$`U\subseteq\mathcal K_m=\mathrm{ran}\,Q_m`$ de dimension $`j`$. La même formule sur $`U\subseteq\mathcal K_{m+1}`$ donne
$`\theta^{(m+1)}_j`$ et sur $`U\subseteq\mathbb C^n`$ donne $`\lambda_j`$ (théorème A.8) ; minimiser sur davantage de sous-espaces donne une valeur plus petite. Si $`\mathcal K_{m+1}=\mathcal K_m`$, alors
$`\mathcal K_m`$ est $`H`$-invariant et contient $`p(H)v`$ pour tout polynôme $`p`$, en particulier $`P_1v\ne0`$ ($`P_1`$ le
projecteur spectral de $`\lambda_1`$, un polynôme en $`H`$ par interpolation de Lagrange, comme dans le théorème 1.5) ;
donc $`\theta^{(m)}_1\le\langle x,Hx\rangle=\lambda_1`$ en $`x=P_1v/\|P_1v\|\in\mathcal K_m`$. $`\square`$

Un petit résidu certifie une valeur propre proche : pour $`H`$ hermitienne, un nombre $`E`$ et
$`\|\psi\|=1`$, $`\|H\psi-E\psi\|^2=\sum_i|\lambda_i-E|^2|\langle u_i,\psi\rangle|^2\ge\min_i|\lambda_i-E|^2`$ ($`u_i`$ une base
orthonormée de vecteurs propres). Il ne certifie pas que cette valeur propre soit la plus petite. Pour $`L=4`$,
l'énergie de l'état fondamental de l'anneau de Heisenberg est exacte : avec
$`\mathbf S_A=\mathbf s_0+\mathbf s_2`$, $`\mathbf S_B=\mathbf s_1+\mathbf s_3`$, l'anneau s'écrit
$`H=J\,\mathbf S_A\cdot\mathbf S_B=\frac J2(\mathbf S^2-\mathbf S_A^2-\mathbf S_B^2)\ge\frac J2(0-2-2)=-2J`$, puisque $`\mathbf S^2\ge0`$ et que
$`\mathbf S_A^2,\mathbf S_B^2`$ ont les valeurs propres $`0,2`$ (proposition 4.13 avec $`L=2`$) ; l'égalité a lieu sur un vecteur
de spin $`0`$ du produit des deux triplets, qui existe d'après le chapitre 4, exercice 2
($`1\otimes1\cong2\oplus1\oplus0`$).

*Calcul :* le théorème 9.14 pour la chaîne XX ouverte, avec `krylov.lanczos` sur le secteur creux de `SpinChain.xxz_sparse`, et des états fondamentaux de l'anneau de Heisenberg jusqu'à $`L=16`$ avec `krylov.ground_state` et son résidu — [tutoriel, section 9.5](../../tutorial/fr/09_spin_chains.md#s-9-5).

**Évolution temporelle.** Pour $`e^{-i\tau H}v`$, on projette sur $`\mathcal K_m`$ : $`w_\tau=\|v\|\,Q_me^{-i\tau T_m}e_1`$.

<a id="r-9-15"></a>

**Proposition 9.15.** Soit $`\dim\mathcal K_m=m`$. Pour tout $`\tau`$, $`\|w_\tau\|=\|v\|`$ et
$`\langle w_\tau,Hw_\tau\rangle=\langle v,Hv\rangle`$. Si $`\mathcal K_m`$ est $`H`$-invariant, $`w_\tau=e^{-i\tau H}v`$.

*Démonstration.* $`Q_m`$ est une isométrie et $`e^{-i\tau T_m}`$ est unitaire ($`T_m`$ est réelle symétrique), ce qui donne
la norme. Comme $`Q_m^*HQ_m=T_m`$ commute avec $`e^{-i\tau T_m}`$,
$`\langle w_\tau,Hw_\tau\rangle=\|v\|^2\langle e^{-i\tau T_m}e_1,T_me^{-i\tau T_m}e_1\rangle=\|v\|^2\langle e_1,T_me_1\rangle=\langle v,Hv\rangle`$. Si
$`\mathcal K_m`$ est invariant, $`HQ_m=Q_mT_m`$ (proposition 9.13, $`r_m=0`$), donc $`H^jQ_m=Q_mT_m^j`$ et, en sommant
la série exponentielle (proposition 8.2(1)), $`e^{-i\tau H}Q_m=Q_me^{-i\tau T_m}`$ ; on applique à $`e_1`$, avec $`Q_me_1=v/\|v\|`$. $`\square`$

*Calcul :* la proposition 9.15 pour `krylov.evolve`, qui avance par pas dans des espaces de Krylov, sur une trempe (quench) à partir de l'état de Néel, comparée au propagateur dense — [tutoriel, section 9.5](../../tutorial/fr/09_spin_chains.md#s-9-5).

**Bornes d'erreur.** Sauf si $`\mathcal K_m`$ est invariant, la proposition 9.15 ne dit rien de la précision de $`w_\tau`$.
L'erreur dépend de la qualité avec laquelle les polynômes de degré $`<m`$ approchent $`e^{-i\tau x}`$ sur le spectre, et
les polynômes de Tchebychev fournissent l'approximation. Les **polynômes de Tchebychev** $`\mathcal T_0=1`$, $`\mathcal T_1(y)=y`$,
$`\mathcal T_{k+1}(y)=2y\mathcal T_k(y)-\mathcal T_{k-1}(y)`$ sont de degré $`k`$ et vérifient $`\mathcal T_k(\cos\vartheta)=\cos k\vartheta`$ (par
récurrence, à partir de $`\cos(k+1)\vartheta+\cos(k-1)\vartheta=2\cos\vartheta\cos k\vartheta`$), donc $`|\mathcal T_k|\le1`$ sur $`[-1,1]`$. On note
```math
S_m(x)=\sum_{j\ge m}\frac{x^j}{j!}\le\frac{x^m}{m!}\Big(1-\frac x{m+1}\Big)^{-1}\qquad(0\le x<m+1),
```
où la borne découle de $`m!/(m+i)!\le(m+1)^{-i}`$ et de la série géométrique.

<a id="r-9-17"></a>

**Lemme 9.17 (développement de Tchebychev de $`e^{-i\rho y}`$).** Pour $`\rho`$ réel et $`y\in[-1,1]`$,
```math
e^{-i\rho y}=a_0+2\sum_{k\ge1}a_k\mathcal T_k(y),\qquad a_k=(-i)^k\sum_{l\ge0}\frac{(-1)^l(\rho/2)^{k+2l}}{l!\,(k+l)!},\qquad|a_k|\le\frac{(|\rho|/2)^k}{k!}.
```
Par conséquent, pour tout $`m\ge1`$, le polynôme $`p=a_0+2\sum_{k=1}^{m-1}a_k\mathcal T_k`$ de degré $`\le m-1`$ vérifie
$`|e^{-i\rho y}-p(y)|\le2S_m(|\rho|/2)`$ sur $`[-1,1]`$.

(Les sommes qui définissent $`a_k`$ sont les fonctions de Bessel, $`a_k=(-i)^kJ_k(\rho)`$ ; la démonstration n'en a pas
besoin.)

*Démonstration.* *Développement.* Pour $`y=\cos\vartheta`$, $`\cos^j\vartheta=2^{-j}(e^{i\vartheta}+e^{-i\vartheta})^j=2^{-j}\sum_{l=0}^j\binom jle^{i(j-2l)\vartheta}`$, donc
```math
e^{-i\rho\cos\vartheta}=\sum_{j\ge0}\frac{(-i\rho)^j}{j!}\,2^{-j}\sum_{l=0}^j\binom jl\,e^{i(j-2l)\vartheta}.
```
La série double converge absolument (les modules de ses termes ont pour somme $`e^{|\rho|}`$), on peut donc la sommer
dans n'importe quel ordre, en particulier selon $`k=j-2l`$. Pour $`k\ge0`$, les termes avec $`j-2l=k`$ ont $`j=k+2l`$, et
comme $`\binom{k+2l}l/(k+2l)!=1/(l!\,(k+l)!)`$, leur somme est $`a_ke^{ik\vartheta}`$ ; pour $`k<0`$, $`\binom jl=\binom j{j-l}`$ montre que le
coefficient de $`e^{ik\vartheta}`$ est $`a_{|k|}`$. D'où
$`e^{-i\rho\cos\vartheta}=a_0+\sum_{k\ge1}a_k(e^{ik\vartheta}+e^{-ik\vartheta})=a_0+2\sum_{k\ge1}a_k\mathcal T_k(\cos\vartheta)`$.

*Borne.* Fixons $`k\ge0`$. Soit $`\vartheta`$ uniforme sur $`[0,2\pi]`$ ; soit $`\sigma=1`$ si $`k=0`$, et pour $`k\ge1`$ soit $`\sigma`$
indépendante de $`\vartheta`$, de densité $`k(1-t)^{k-1}`$ sur $`[0,1]`$. Posons $`Z=\sqrt\sigma\cos\vartheta\in[-1,1]`$. Alors
$`\mathbb E\cos^{2l}\vartheta=4^{-l}\binom{2l}l`$ (dans le développement binomial ci-dessus, seul le terme constant est de
moyenne non nulle), et $`\mathbb E\sigma^l=k\int_0^1t^l(1-t)^{k-1}dt=l!\,k!/(k+l)!`$ (intégrer $`l`$ fois par parties,
$`\int_0^1t^l(1-t)^{k-1}dt=\frac lk\int_0^1t^{l-1}(1-t)^kdt`$ ; pour $`k=0`$ les deux membres valent $`1`$). Par indépendance,
$`\mathbb EZ^{2l}=(2l)!\,k!/(4^l\,l!\,(k+l)!)`$. La série $`\cos(\rho Z)=\sum_l(-1)^l(\rho Z)^{2l}/(2l)!`$ converge uniformément,
puisque $`|Z|\le1`$, donc on peut en prendre l'espérance terme à terme :
```math
\frac{(\rho/2)^k}{k!}\,\mathbb E\cos(\rho Z)=\sum_{l\ge0}\frac{(-1)^l\rho^{2l}}{(2l)!}\cdot\frac{(\rho/2)^k\,(2l)!}{4^l\,l!\,(k+l)!}=\sum_{l\ge0}\frac{(-1)^l(\rho/2)^{k+2l}}{l!\,(k+l)!}=i^ka_k.
```
Comme $`|\cos|\le1`$, $`|a_k|\le(|\rho|/2)^k/k!`$.

*Troncature.* $`|e^{-i\rho y}-p(y)|=2\,|\sum_{k\ge m}a_k\mathcal T_k(y)|\le2\sum_{k\ge m}(|\rho|/2)^k/k!=2S_m(|\rho|/2)`$. $`\square`$

<a id="r-9-18"></a>

**Théorème 9.18 (erreur du propagateur de Krylov).** Soient $`\dim\mathcal K_m=m\ge2`$ et $`\|v\|=1`$, $`Q_m,T_m,r_m`$ comme
dans la proposition 9.13, $`\beta_m=\|r_m\|`$, $`\theta_1\le\dots\le\theta_m`$ les valeurs de Ritz, $`c_T=(\theta_1+\theta_m)/2`$ et
$`w_T=(\theta_m-\theta_1)/2`$. Soient $`w_\tau=Q_me^{-i\tau T_m}e_1`$ et $`\varphi(t)=\langle e_m,e^{-itT_m}e_1\rangle`$, et soient
$`c\in\mathbb R`$, $`w\ge0`$ avec $`\mathrm{spec}\,H\subseteq[c-w,c+w]`$. Pour $`\tau`$ réel :

(a) (*invariance par translation*) pour $`s`$ réel, les vecteurs de Lanczos de $`H+s`$ et $`v`$ sont ceux de $`H`$ et $`v`$,
$`T_m`$ devient $`T_m+s`$, $`r_m`$ ne change pas, et $`e^{-i\tau H}v`$ et $`w_\tau`$ sont tous deux multipliés par $`e^{-i\tau s}`$ ;
l'erreur $`\|e^{-i\tau H}v-w_\tau\|`$ et les bornes (b)–(d) ne changent pas ;

(b) (*a posteriori*) $`e^{-i\tau H}v-w_\tau=-i\int_0^\tau\varphi(t)\,e^{-i(\tau-t)H}r_m\,dt`$, d'où
$`\|e^{-i\tau H}v-w_\tau\|\le\beta_m\,\big|\int_0^\tau|\varphi(t)|\,dt\big|`$ ;

(c) (*borne calculable*) $`w_T>0`$, $`|\varphi(t)|\le2S_{m-1}(|t|w_T/2)`$, et
```math
\|e^{-i\tau H}v-w_\tau\|\le\frac{4\beta_m}{w_T}\,S_m\Big(\frac{|\tau|\,w_T}2\Big);
```

(d) (*a priori*) $`\beta_m\le w`$ et $`w_T\le w`$, et avec $`\rho=|\tau|\,w`$
```math
\|e^{-i\tau H}v-w_\tau\|\le\frac{4\beta_m}w\,S_m\Big(\frac\rho2\Big)\le4S_m\Big(\frac\rho2\Big)\le\frac4{m!}\Big(\frac\rho2\Big)^m\Big(1-\frac\rho{2m+2}\Big)^{-1}\qquad(\rho<2m+2).
```

*Démonstration.* (a) $`(H+s)^jv`$ est $`H^jv`$ plus une combinaison des $`H^iv`$, $`i<j`$, et réciproquement, donc
$`\mathcal K_j(H+s,v)=\mathcal K_j(H,v)`$ pour tout $`j`$. La composante de $`(H+s)q_j=Hq_j+sq_j`$ orthogonale à $`\mathcal K_j\ni q_j`$
est celle de $`Hq_j`$, donc la construction de la proposition 9.13 donne les mêmes $`q_j`$ ; alors $`Q_m^*(H+s)Q_m=T_m+s`$ et
$`(H+s)Q_m-Q_m(T_m+s)=HQ_m-Q_mT_m`$. D'après le théorème 1.5(c), $`e^{-i\tau(X+s)}=e^{-i\tau s}e^{-i\tau X}`$ pour $`X=H`$ et
$`X=T_m`$. Enfin $`\beta_m`$, $`|\varphi|`$ et $`w_T`$ ne changent pas, et $`\mathrm{spec}(H+s)\subseteq[c+s-w,c+s+w]`$.

(b) D'après la proposition 8.2(1) pour $`T_m`$, $`\frac d{dt}w_t=-iQ_mT_me^{-itT_m}e_1`$, et $`Q_mT_m=HQ_m-r_me_m^T`$
(proposition 9.13), donc $`\frac d{dt}w_t=-iHw_t+i\varphi(t)r_m`$. Comme $`\frac d{dt}e^{itH}=iHe^{itH}`$ (proposition 8.2(1)), la
règle du produit donne $`\frac d{dt}(e^{itH}w_t)=i\varphi(t)e^{itH}r_m`$. En intégrant de $`0`$ à $`\tau`$ (lemme A.4(e)), avec
$`w_0=Q_me_1=v`$, $`e^{i\tau H}w_\tau-v=i\int_0^\tau\varphi(t)e^{itH}r_m\,dt`$. On applique l'unitaire $`e^{-i\tau H}`$ (lemme A.4(b))
et on prend les normes (lemme A.4(c)).

(c) $`T_m`$ n'est pas un multiple de $`1`$ ($`\beta_1>0`$), donc ses valeurs propres ne sont pas toutes égales (théorème 1.3)
et $`w_T>0`$. Comme $`T_m`$ est tridiagonale, $`T_m^je_1`$ s'annule sous son coefficient $`j+1`$, donc
$`\langle e_m,p(T_m)e_1\rangle=0`$ pour tout polynôme $`p`$ de degré $`\le m-2`$. D'après le lemme 9.17, avec $`m-1`$ au lieu de
$`m`$ et $`\rho=tw_T`$, il existe un polynôme $`\tilde p`$ de degré $`\le m-2`$ avec $`|e^{-itw_Ty}-\tilde p(y)|\le2S_{m-1}(|t|w_T/2)`$
sur $`[-1,1]`$ ; alors $`p(x)=e^{-itc_T}\tilde p((x-c_T)/w_T)`$ est de degré $`\le m-2`$ et
$`|e^{-itx}-p(x)|\le2S_{m-1}(|t|w_T/2)`$ en chaque valeur de Ritz. D'après le théorème 1.5 et le corollaire 1.4(1),
$`|\varphi(t)|=|\langle e_m,(e^{-itT_m}-p(T_m))e_1\rangle|\le\max_i|e^{-it\theta_i}-p(\theta_i)|`$. La série $`S_{m-1}`$ est à termes
positifs et converge uniformément sur les intervalles bornés, on peut donc l'intégrer terme à terme :
$`\int_0^{|\tau|}2S_{m-1}(tw_T/2)\,dt=\frac4{w_T}S_m(|\tau|w_T/2)`$ ; avec (b), c'est la borne.

(d) $`Q_m(T_m-c)e_m=Q_mQ_m^*(H-c)q_m`$ est la projection orthogonale de $`(H-c)q_m`$ sur $`\mathcal K_m`$, et
$`r_m=(H-c)q_m-Q_m(T_m-c)e_m`$ est la composante orthogonale, donc $`\beta_m\le\|(H-c)q_m\|\le\|H-c\|\le w`$ (corollaire 1.4(1)).
Pour un vecteur propre unitaire $`y`$ de $`T_m`$, la valeur de Ritz $`\langle y,T_my\rangle=\langle Q_my,HQ_my\rangle`$ est dans
$`[c-w,c+w]`$, donc $`w_T\le w`$. Or $`\frac{4\beta_m}{w_T}S_m(x)=2\beta_m|\tau|\,S_m(x)/x`$ avec $`x=|\tau|w_T/2\le\rho/2`$, et
$`S_m(x)/x=\sum_{j\ge m}x^{j-1}/j!`$ est croissante ; cela donne la première inégalité, $`\beta_m\le w`$ la deuxième et la
borne sur $`S_m`$ la troisième. $`\square`$

La borne (c) n'utilise que les données du calcul de Lanczos, $`\beta_m`$ et les valeurs de Ritz. Les valeurs de Ritz sont
à l'intérieur du spectre (théorème 9.14) et peuvent beaucoup sous-estimer sa demi-largeur $`w`$, pour un vecteur de
départ ayant de petites composantes sur les vecteurs propres extrêmes ; mais aucune estimation de $`w`$ n'entre dans
(c), car la représentation (b) ne fait intervenir que $`T_m`$. D'après (d) et $`m!\ge(m/e)^m`$, l'erreur est au plus de
l'ordre de $`(e\rho/2m)^m`$ : elle est petite dès que $`\rho=|\tau|w`$ est nettement inférieur à $`2m/e\approx0.74\,m`$, et
les pas peuvent donc être longs. La quantité $`\beta_m|\varphi(\tau)|`$, que les versions antérieures de la bibliothèque
additionnaient comme estimation d'erreur, est, d'après la démonstration de (b), la norme du défaut
$`\frac d{d\tau}w_\tau+iHw_\tau`$ : l'intégrande de (b) à l'extrémité, un taux et non une erreur. Pour $`\tau`$ petit,
$`\varphi(t)\approx(-it)^{m-1}\beta_1\cdots\beta_{m-1}/(m-1)!`$, si bien qu'elle dépasse la borne (b) d'un facteur d'environ
$`m/|\tau|`$ ; pour $`\tau`$ grand, $`\varphi`$ oscille et peut s'annuler à l'extrémité sans que l'erreur s'annule.

<a id="r-9-19"></a>

**Corollaire 9.19 (règle de pas).** Soient $`\tau=\tau_1+\dots+\tau_N`$ avec $`\tau_k>0`$, $`\psi_0=v\ne0`$,
$`\psi_k=\|\psi_{k-1}\|\,Q_me^{-i\tau_kT_m}e_1`$ l'approximation de Krylov de $`e^{-i\tau_kH}\psi_{k-1}`$ dans
$`\mathcal K_m(H,\psi_{k-1})`$, et $`\varepsilon_k`$ la borne (c) de ce pas multipliée par $`\|\psi_{k-1}\|`$ ($`\varepsilon_k=0`$ si cet
espace de Krylov est de dimension $`<m`$). Alors $`\|e^{-i(\tau_1+\dots+\tau_k)H}v-\psi_k\|\le\varepsilon_1+\dots+\varepsilon_k`$ pour tout
$`k`$. En particulier, si $`\varepsilon_k\le\mathrm{tol}\cdot\|v\|\,\tau_k/\tau`$ à chaque pas, l'erreur après $`k`$ pas est au plus
$`\mathrm{tol}\cdot\|v\|\,(\tau_1+\dots+\tau_k)/\tau`$. Pour $`m\ge2`$, la borne d'un pas de longueur $`h`$, divisée par $`h`$, est
croissante en $`h`$ et tend vers $`0`$ quand $`h\to0`$ ; le plus long pas admissible est donc positif.

*Démonstration.* Soit $`d_k=e^{-i(\tau_1+\dots+\tau_k)H}v-\psi_k`$. Alors $`d_k=e^{-i\tau_kH}d_{k-1}+(e^{-i\tau_kH}\psi_{k-1}-\psi_k)`$, et
$`e^{-i\tau_kH}`$ est unitaire, donc $`\|d_k\|\le\|d_{k-1}\|+\varepsilon_k`$ d'après le théorème 9.18(c) appliqué à
$`\psi_{k-1}/\|\psi_{k-1}\|`$ ; si l'espace de Krylov est de dimension $`j<m`$, il est égal à $`\mathcal K_j`$ et invariant, et le
pas est exact d'après la proposition 9.15. Toutes les $`\|\psi_k\|=\|v\|`$ d'après la proposition 9.15. La borne d'un pas de
longueur $`h`$ divisée par $`h`$ vaut $`2\beta_m\|\psi_{k-1}\|S_m(x)/x`$ avec $`x=hw_T/2`$, croissante en $`h`$, et
$`S_m(x)/x\le\frac{x^{m-1}}{m!}(1-\frac x{m+1})^{-1}\to0`$ pour $`m\ge2`$. $`\square`$

*Dans la bibliothèque :* `krylov.expm_multiply_step` renvoie la borne (c) (`krylov.krylov_error_bound`) avec chaque pas.
`krylov.evolve` (par défaut `dt=None`) construit l'espace de Krylov en l'état courant et fait le plus long pas $`h`$,
jusqu'au prochain instant de sortie, dont la borne (c) est au plus $`\mathrm{tol}\cdot\|\psi_0\|\,h/(t_{\mathrm{end}}-t_0)`$ ; d'après le
corollaire 9.19, chaque état renvoyé est alors à distance au plus $`\mathrm{tol}\cdot\|\psi_0\|`$ de $`e^{-i(t-t_0)H}\psi_0`$. Cela vaut
en arithmétique exacte : les erreurs d'arrondi ne sont pas comprises (elles croissent avec le nombre de pas, et avec
$`\varepsilon\|H\|`$ par unité de temps pour une grande translation constante de $`H`$, que porte chaque produit
matrice–vecteur). D'après le théorème 9.18(a), les pas ne dépendent pas d'une translation constante. Avec un `dt` fixe,
la somme des bornes est comparée à $`\mathrm{tol}\cdot\|\psi_0\|`$, avec un avertissement si elle est plus grande.

*Calcul :* le lemme 9.17, le théorème 9.18 et le corollaire 9.19 comparés au propagateur dense, pour des spectres avec
des amas, une grande translation et une grande norme, et les pas de `krylov.evolve` — [tutoriel, section 9.5](../../tutorial/fr/09_spin_chains.md#s-9-5).

<a id="s-9-6"></a>

## 9.6 Désordre : statistiques de niveaux et intrication

La **chaîne de Heisenberg en champ aléatoire** a des champs indépendants $`h_i`$ uniformes sur $`[-W,W]`$. Elle conserve
$`S^z`$ (proposition 9.1), donc les statistiques de niveaux se calculent par secteur $`N`$ (théorème 5.15 et
section 5.4 : mélanger les secteurs simule des niveaux indépendants). Si deux champs diffèrent,
$`THT^*-H=\sum_i(h_{i-1}-h_i)s^z_i\ne0`$, donc les translations ne sont plus une symétrie.

*Calcul :* le rapport moyen $`\langle r\rangle`$ de la chaîne en champ aléatoire, comparé à des matrices du GOE et à des niveaux indépendants (observations numériques) — [tutoriel, section 9.6](../../tutorial/fr/09_spin_chains.md#s-9-6).

**Intrication des fermions libres.** Pour les états propres de la chaîne XX ouverte, l'entropie d'intrication des
$`n`$ premiers sites est déterminée par une matrice $`n\times n`$.

<a id="r-9-16"></a>

**Théorème 9.16 (formule de Peschel).** Soient $`\Phi=b_{k_1}^*\cdots b_{k_N}^*|0\rangle`$ comme dans le théorème 9.10,
$`C_{ij}=\langle\Phi,c_i^*c_j\Phi\rangle`$ et $`C_A=(C_{ij})_{i,j<n}`$. Alors $`0\le C_A\le1`$, et l'entropie d'intrication des
sites $`A=\{0,\dots,n-1\}`$ est, $`\nu`$ parcourant les valeurs propres de $`C_A`$ comptées avec multiplicité,
```math
S_A=\sum_\nu \eta(\nu),\qquad \eta(\nu)=-\nu\log\nu-(1-\nu)\log(1-\nu)\quad(0\log0=0).
```

*Démonstration.* Pour un vecteur réel $`x`$, posons $`b(x)=\sum_jx_jc_j`$ ; alors $`b(x)^*b(y)^*=-b(y)^*b(x)^*`$ et
$`\{b(x),b(y)^*\}=\langle x,y\rangle`$. Soient $`U\subseteq\mathbb R^L`$ l'espace engendré par les $`u_{k_a}`$, et $`P_A`$,
$`P_B=1-P_A`$ les projecteurs de coordonnées sur $`\{0,\dots,n-1\}`$ et sur son complémentaire. Diagonaliser la forme
$`\langle P_Ax,P_Ay\rangle`$ sur $`U`$ donne une base orthonormée $`\phi_1,\dots,\phi_N`$ de $`U`$ avec
$`\langle P_A\phi_a,P_A\phi_b\rangle=\nu_a\delta_{ab}`$, où $`\nu_a=\|P_A\phi_a\|^2\in[0,1]`$, d'où $`\langle P_B\phi_a,P_B\phi_b\rangle=(1-\nu_a)\delta_{ab}`$. Comme
$`(x_1,\dots,x_N)\mapsto b(x_1)^*\cdots b(x_N)^*|0\rangle`$ est multilinéaire et alternée,
$`\Phi=\pm\prod_a\big(\sqrt{\nu_a}f_a^*+\sqrt{1-\nu_a}\,g_a^*\big)|0\rangle`$ avec $`f_a=b(P_A\phi_a/\sqrt{\nu_a})`$,
$`g_a=b(P_B\phi_a/\sqrt{1-\nu_a})`$ (en omettant les termes de coefficient nul). En développant,
$`\Phi=\sum_I\pm\prod_{a\in I}\sqrt{\nu_a}\prod_{a\notin I}\sqrt{1-\nu_a}\,F_I^*G_{I^c}^*|0\rangle`$, la somme portant sur les
parties $`I`$, avec $`F_I^*`$, $`G_{I^c}^*`$ les produits ordonnés. Pour $`j\ge n`$, la chaîne de signes de $`c_j^*`$ contient
tout $`A`$ : $`c_j^*=(-1)^{N_A}\otimes\tilde c_j^*`$ avec des opérateurs de Jordan–Wigner $`\tilde c_j`$ des sites $`\ge n`$,
tandis que pour $`i<n`$, $`c_i^*=\hat c_i^*\otimes1`$. Comme $`(-1)^{N_A}`$ vaut $`1`$ sur $`|0_A\rangle`$,
$`F_I^*G_{I^c}^*|0\rangle=\xi_I\otimes\chi_{I^c}`$ avec $`\xi_I=\hat F_I^*|0_A\rangle`$, $`\chi_{I^c}=\tilde G_{I^c}^*|0_B\rangle`$, qui
sont des familles orthonormées d'après le calcul avec les CAR du théorème 9.10. C'est une décomposition de Schmidt
(lemme 7.7) de poids $`p_I=\prod_{a\in I}\nu_a\prod_{a\notin I}(1-\nu_a)`$ (les signes passent dans $`\xi_I`$), la loi de
variables de Bernoulli indépendantes, donc, comme $`\log p_I=\sum_{a\in I}\log\nu_a+\sum_{a\notin I}\log(1-\nu_a)`$ et
$`\sum_{I\ni a}p_I=\nu_a`$, $`S_A=-\sum_Ip_I\log p_I=\sum_a\eta(\nu_a)`$. Enfin $`C_{ij}=\sum_a\phi_a(i)\phi_a(j)`$ (développer
$`c_i^*c_j`$ sur les $`b(\phi_a)`$), donc $`C_A=\sum_a(P_A\phi_a)(P_A\phi_a)^T`$ a les valeurs propres $`\nu_a`$ sur les vecteurs
orthogonaux $`P_A\phi_a`$ et $`0`$ ailleurs, et $`\eta(0)=0`$. $`\square`$

Le théorème 9.16 est dû à I. Peschel (*J. Phys. A* 36 (2003) L205).

*Calcul :* le théorème 9.16 pour les états fondamentaux de la chaîne XX ouverte avec $`L=12`$, comparé à `SpinChain.entanglement_entropy` — [tutoriel, section 9.6](../../tutorial/fr/09_spin_chains.md#s-9-6).

<a id="s-9-7"></a>

## 9.7 Points de comparaison

*Calcul :* la comparaison avec d'autres logiciels (QuSpin, codes de diagonalisation exacte, réseaux de tenseurs) — [tutoriel, section 9.7](../../tutorial/fr/09_spin_chains.md#s-9-7).

<a id="x-9"></a>

## Exercices

1. (Secteurs d'impulsion à la main) Pour $`L=4`$ et $`N=2`$, trouver les orbites des six configurations sous $`T`$ et les
   dimensions $`d_{2,m}`$, $`m=0,\dots,3`$.
2. (Règles de somme) Soit $`\psi`$ l'état fondamental de l'anneau de Heisenberg avec $`L`$ pair. Montrer que
   $`E_0=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle`$ et $`\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle=0`$.
3. (Trop de symétrie) Montrer que, dans le secteur $`N=L/2`$, le vecteur dont tous les coefficients valent $`1`$ est un
   multiple de $`(S^-)^{L/2}|{\uparrow\cdots\uparrow}\rangle`$ et un vecteur propre de l'anneau de Heisenberg d'énergie
   $`JL/4`$, et que tout espace de Krylov $`\mathcal K_m(H,v)`$ est contenu dans chaque sous-espace $`H`$-invariant
   contenant $`v`$. Que renvoie `krylov.ground_state` lorsqu'on part de ce vecteur ?

<details>
<summary>Solutions</summary>

1. $`0011\to1001\to1100\to0110\to0011`$ est une orbite avec $`R=4`$, et $`0101\to1010\to0101`$ une orbite avec $`R=2`$.
   D'après la proposition 9.5, la première contribue à tout $`m`$, la seconde aux $`m`$ tels que $`2m\equiv0\pmod4`$,
   c'est-à-dire $`m=0,2`$ : $`d_{2,m}=2,1,2,1`$.

2. L'état fondamental est un singulet simple (théorème 9.9) et $`T`$ commute avec $`H`$, donc $`T\psi=e^{ik}\psi`$. Comme
   $`\mathbf s_i\cdot\mathbf s_{i+1}=T^i(\mathbf s_0\cdot\mathbf s_1)T^{*i}`$, $`\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle`$ est le même
   pour les $`L`$ liaisons, et $`E_0=J\sum_i\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle`$. Comme
   $`\mathbf S^2\psi=0`$,
   $`0=\langle\mathbf S^2\rangle=\sum_{i,j}\langle\mathbf s_i\cdot\mathbf s_j\rangle=L\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle`$, là encore
   par invariance par translation.

3. $`S^-`$ abaisse l'un des spins vers le haut ; en l'appliquant $`L/2`$ fois à l'état complètement polarisé, on
   atteint chaque configuration à $`L/2`$ spins vers le haut selon $`(L/2)!`$ ordres, tous avec le coefficient $`1`$.
   $`|{\uparrow\cdots\uparrow}\rangle`$ a ses $`L`$ liaisons parallèles, d'énergie $`JL/4`$, et $`S^-`$ commute avec $`H`$, donc
   le vecteur est un vecteur propre de même énergie. Si $`W`$ est invariant et $`v\in W`$, alors $`H^jv\in W`$ pour tout
   $`j`$. Ici $`\mathcal K_m=\mathbb Cv`$ pour tout $`m`$, et `ground_state` renvoie $`E=JL/4`$, le sommet du spectre (chaque
   $`\mathbf s_i\cdot\mathbf s_j=\frac12F_{ij}-\frac14\le\frac14`$ d'après la proposition 4.15(a), car $`F_{ij}`$ est une matrice de
   permutation symétrique de valeurs propres $`\pm1`$). Un vecteur de départ aléatoire évite ce piège : un $`v`$ gaussien
   a une densité, donc il a presque sûrement une composante non nulle $`P_1v`$ dans le sous-espace propre de
   $`\lambda_1`$ (corollaire A.7(c) ; pour un $`v`$ réel, on l'applique dans $`W=\mathbb R^n`$ au sous-espace strict
   $`\{v:P_1v=0\}`$, strict parce que la matrice réelle $`P_1\ne0`$ a une colonne non nulle), et alors, d'après le
   théorème 9.14, un espace de Krylov $`H`$-invariant $`\mathcal K_m(H,v)`$ a pour plus petite valeur de Ritz $`\lambda_1`$,
   alors qu'ici l'espace invariant $`\mathcal K_1=\mathbb Cv`$ a la valeur de Ritz $`JL/4`$.

</details>

*Calcul :* les exercices 1–3, numériquement — [tutoriel, Exercices](../../tutorial/fr/09_spin_chains.md#x-9).

**Suite :** le [chapitre 10](10_quantum_computing.md) traite du calcul quantique avec des secteurs : sous-espaces
sans décohérence et sous-systèmes sans bruit, correction d'erreurs dans les algèbres d'opérateurs, et mélangeurs
qui respectent un poids de Hamming conservé — l'algèbre $`\bigoplus_wM_{\binom nw}`$ de ce chapitre dans un nouveau
rôle.
