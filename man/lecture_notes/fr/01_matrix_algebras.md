[English](../01_matrix_algebras.md) · [Español](../es/01_matrix_algebras.md) · **Français** · [Deutsch](../de/01_matrix_algebras.md) · [中文](../zh/01_matrix_algebras.md) · [日本語](../ja/01_matrix_algebras.md) · [Русский](../ru/01_matrix_algebras.md)

# 1. Algèbres de matrices et calcul fonctionnel

Ce chapitre fixe les notations et démontre les propriétés de $`M_n(\mathbb C)`$ dont se
servent tous les chapitres suivants : le théorème spectral, le calcul fonctionnel et sa dérivée, la
décomposition polaire et les normes définies par les valeurs singulières. Le code de la
bibliothèque et les vérifications numériques de ces résultats se trouvent dans le
[Tutoriel 1](../../tutorial/fr/01_matrix_algebras.md).

**Conventions.** Les produits scalaires $`\langle x,y\rangle=\sum_i\bar x_iy_i`$ sont linéaires
par rapport au second argument. $`A^*`$ désigne la transposée conjuguée, de sorte que
$`\langle x,Ay\rangle=\langle A^*x,y\rangle`$. On note $`e_1,\dots,e_n`$ la base canonique, $`1`$
la matrice identité et $`\mathrm{Tr}`$ la trace usuelle.

<a id="s-1-1"></a>

## 1.1 L'algèbre $`M_n(\mathbb C)`$

$`M_n(\mathbb C)`$ est l'algèbre des matrices complexes $`n\times n`$, munie de l'involution
$`A\mapsto A^*`$ et de la norme d'opérateur $`\|A\|=\sup_{\|x\|=1}\|Ax\|`$. Un élément $`A`$ est
dit

| nom | définition |
|---|---|
| normal | $`A^*A=AA^*`$ |
| autoadjoint | $`A=A^*`$ |
| positif, $`A\ge0`$ | $`\langle x,Ax\rangle\ge0`$ pour tout $`x`$ |
| unitaire | $`A^*A=AA^*=1`$ |
| projecteur | $`A=A^*=A^2`$ |

On écrit $`A\le B`$ si $`B-A\ge0`$, et $`A>0`$ si $`A\ge0`$ et $`A`$ est inversible.

*Dans la bibliothèque :* les opérateurs de spectre prescrit, autoadjoints ou positifs, `TypeIAlgebra.operator_from_eigenvalues` — [tutoriel, section 1.5](../../tutorial/fr/01_matrix_algebras.md#s-1-5) ; la norme d'opérateur, `operator_norm` — [tutoriel, section 1.6](../../tutorial/fr/01_matrix_algebras.md#s-1-6).

<a id="r-1-1"></a>

**Lemme 1.1 (identité C\*).** $`\|A^*\|=\|A\|`$ et $`\|A^*A\|=\|A\|^2`$.

*Démonstration.* $`\|Ax\|=\sup_{\|y\|=1}|\langle y,Ax\rangle|=\sup_{\|y\|=1}|\langle A^*y,x\rangle|`$,
donc $`\|A\|=\sup_{\|x\|=\|y\|=1}|\langle A^*y,x\rangle|=\|A^*\|`$. Ensuite,
$`\|Ax\|^2=\langle x,A^*Ax\rangle\le\|A^*A\|`$ pour $`\|x\|=1`$, d'où $`\|A\|^2\le\|A^*A\|\le\|A^*\|\,\|A\|=\|A\|^2`$. $`\square`$

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). En
> mécanique quantique, un système possédant $`n`$ configurations distinguables est décrit par
> $`\mathbb C^n`$. Les matrices autoadjointes sont les *observables* (énergie, spin, position sur
> un réseau) ; leurs valeurs propres sont les résultats possibles d'une mesure. Les projecteurs sont les questions à réponse
> oui/non. Les matrices unitaires sont les symétries et les évolutions temporelles,
> $`e^{-itH}`$ pour un hamiltonien $`H`$. Les matrices positives de trace un sont les *états*
> (matrices densité, chapitre 3).

<a id="s-1-2"></a>

## 1.2 Le théorème spectral

<a id="r-1-2"></a>

**Lemme 1.2 (Schur).** Pour tout $`A\in M_n(\mathbb C)`$, il existe une matrice unitaire $`U`$
telle que $`U^*AU`$ soit triangulaire supérieure.

*Démonstration.* Par récurrence sur $`n`$ ; le cas $`n=1`$ est trivial. Le polynôme
caractéristique de $`A`$ admet une racine $`\lambda`$ (théorème de d'Alembert-Gauss), donc
$`Av=\lambda v`$ pour un certain vecteur unitaire $`v`$. On complète $`v`$ en une base
orthonormée $`v,w_2,\dots,w_n`$ et l'on note $`V`$ la matrice unitaire ayant ces vecteurs pour
colonnes. Alors $`V^*AV=\begin{pmatrix}\lambda&*\\0&A'\end{pmatrix}`$ avec
$`A'\in M_{n-1}(\mathbb C)`$. Par hypothèse de récurrence, $`U'^*A'U'`$ est triangulaire
supérieure pour une certaine matrice unitaire $`U'`$, et $`U=V\,(1\oplus U')`$ convient. $`\square`$

<a id="r-1-3"></a>

**Théorème 1.3 (théorème spectral).** $`A`$ est normale si et seulement si
$`A=U\,\mathrm{diag}(\lambda_1,\dots,\lambda_n)\,U^*`$ pour une matrice unitaire $`U`$. Les
$`\lambda_i`$ sont alors les valeurs propres de $`A`$, comptées avec multiplicité.

*Démonstration.* Si $`A=U\Lambda U^*`$ avec $`\Lambda`$ diagonale, alors
$`A^*A=U\bar\Lambda\Lambda U^*=U\Lambda\bar\Lambda U^*=AA^*`$. Réciproquement, soit $`A`$
normale et $`T=U^*AU`$ triangulaire supérieure (lemme 1.2) ; $`T`$ est normale. Comparons les
coefficients d'indice $`(1,1)`$ : $`(T^*T)_{11}=|t_{11}|^2`$ et $`(TT^*)_{11}=\sum_j|t_{1j}|^2`$,
donc $`t_{1j}=0`$ pour $`j>1`$. Ainsi $`T=t_{11}\oplus T'`$ avec $`T'`$ normale et triangulaire
supérieure, et une récurrence montre que $`T`$ est diagonale. La diagonale de $`U^*AU`$ est
formée des racines de $`\det(A-\lambda)=\det(T-\lambda)`$. $`\square`$

*Dans la bibliothèque :* le théorème 1.3 lu à l'envers, $`A=U\Lambda U^*`$ avec des valeurs propres prescrites et $`U`$ aléatoire de loi de Haar, est `TypeIAlgebra.operator_from_eigenvalues` — [tutoriel, section 1.5](../../tutorial/fr/01_matrix_algebras.md#s-1-5).

<a id="r-1-4"></a>

**Corollaire 1.4.** Soit $`A`$ normale, de valeurs propres $`\lambda_i`$. Alors

1. $`\|A\|=\max_i|\lambda_i|`$ ;
2. $`A`$ est autoadjointe si et seulement si tous les $`\lambda_i\in\mathbb R`$ ; unitaire si et
   seulement si tous les $`|\lambda_i|=1`$ ; un projecteur si et seulement si tous les
   $`\lambda_i\in\{0,1\}`$ ;
3. $`A\ge0`$ si et seulement si $`A`$ est autoadjointe et tous les $`\lambda_i\ge0`$. De plus, toute matrice
   positive est autoadjointe, donc normale, si bien que (3) caractérise la positivité parmi *toutes*
   les matrices et que le calcul fonctionnel ci-dessous s'applique à toute $`A\ge0`$.

*Démonstration.* Avec $`y=U^*x`$, on a $`\|Ax\|^2=\sum_i|\lambda_i|^2|y_i|^2`$ et $`\|y\|=\|x\|`$,
ce qui donne (1), avec égalité pour $`x=Ue_i`$, $`|\lambda_i|`$ maximal. Pour (2),
$`A^*=U\bar\Lambda U^*`$, $`A^*A=U|\Lambda|^2U^*`$ et $`A^2=U\Lambda^2U^*`$, et les conditions
sont vérifiées si et seulement si elles le sont coefficient par coefficient pour $`\Lambda`$.
Pour (3), remarquons d'abord que tout $`A`$ positif est autoadjoint, puisque
$`\langle x,(A-A^*)x\rangle=\langle x,Ax\rangle-\overline{\langle x,Ax\rangle}=0`$ pour tout $`x`$,
et qu'une matrice $`D`$ telle que $`\langle x,Dx\rangle=0`$ pour tout $`x`$ est nulle (identité de
polarisation : $`4\langle x,Dy\rangle=\sum_{k=0}^3 i^{-k}\langle x+i^ky,D(x+i^ky)\rangle`$).
Ensuite $`\langle x,Ax\rangle=\sum_i\lambda_i|y_i|^2`$, quantité positive pour tout $`x`$ si et
seulement si tous les $`\lambda_i\ge0`$. $`\square`$

<a id="s-1-3"></a>

## 1.3 Le calcul fonctionnel

Soit $`A=U\Lambda U^*`$ normale et $`f`$ une fonction complexe définie sur $`\mathrm{spec}\,A`$.
On pose
```math
f(A)=U\,\mathrm{diag}\big(f(\lambda_1),\dots,f(\lambda_n)\big)\,U^*.
```

<a id="r-1-5"></a>

**Théorème 1.5.** (a) $`f(A)`$ ne dépend pas du choix de $`U`$. (b) Pour un polynôme $`p`$,
$`p(A)`$ est le polynôme usuel en $`A`$. (c) L'application $`f\mapsto f(A)`$ est linéaire et
multiplicative, $`\bar f(A)=f(A)^*`$, et $`\mathrm{spec}\,f(A)=f(\mathrm{spec}\,A)`$.
(d) Si $`B`$ commute avec $`A`$, alors $`B`$ commute avec $`f(A)`$. (e) Si $`g`$ est définie sur
$`f(\mathrm{spec}\,A)`$, alors $`g(f(A))=(g\circ f)(A)`$.

*Démonstration.* (b) $`A^k=U\Lambda^kU^*`$, donc $`p(A)=Up(\Lambda)U^*`$ par linéarité. (a) Le
spectre est fini, donc l'interpolation de Lagrange fournit un polynôme $`p`$ égal à $`f`$ sur
$`\mathrm{spec}\,A`$ : si $`\mu_1,\dots,\mu_r`$ sont les valeurs propres distinctes, on prend
$`p(t)=\sum_if(\mu_i)\prod_{j\neq i}(t-\mu_j)/(\mu_i-\mu_j)`$. Alors $`U\,f(\Lambda)\,U^*=U\,p(\Lambda)\,U^*=p(A)`$ pour *toute*
matrice unitaire $`U`$ diagonalisant $`A`$, et $`p(A)`$ se définit sans recourir à $`U`$. (c) est
vrai pour les matrices diagonales $`f(\Lambda)`$ et se conserve par conjugaison par $`U`$ ; le
spectre de $`Uf(\Lambda)U^*`$ est la diagonale de $`f(\Lambda)`$. (d) $`B`$ commute avec tout
polynôme en $`A`$, en particulier avec $`p(A)=f(A)`$. (e) $`f(A)=Uf(\Lambda)U^*`$ est normale et
diagonalisée par la même matrice $`U`$, donc d'après (a) $`g(f(A))=U\,g(f(\Lambda))\,U^*=(g\circ f)(A)`$. $`\square`$

<a id="r-1-6"></a>

**Théorème 1.6 (racine carrée).** Pour $`A\ge0`$, il existe un unique $`B\ge0`$ tel que
$`B^2=A`$, à savoir $`B=\sqrt A`$ défini par le calcul fonctionnel.

*Démonstration.* $`\sqrt A\ge0`$ et $`(\sqrt A)^2=A`$ d'après le théorème 1.5(c) et le
corollaire 1.4(3). Soit $`B\ge0`$ tel que $`B^2=A`$, $`B=V\,\mathrm{diag}(\mu_j)\,V^*`$,
$`\mu_j\ge0`$. Alors $`A=V\,\mathrm{diag}(\mu_j^2)\,V^*`$, donc $`\mathrm{spec}\,A=\{\mu_j^2\}`$.
Soit $`p`$ un polynôme tel que $`p(t)=\sqrt t`$ sur $`\mathrm{spec}\,A`$. Alors
$`p(A)=\sqrt A`$ et $`p(A)=p(B^2)=V\,\mathrm{diag}(p(\mu_j^2))\,V^*=V\,\mathrm{diag}(\mu_j)\,V^*=B`$. $`\square`$

*Calcul :* les théorèmes 1.5 et 1.6, vérifiés avec `apply_function`, `sqrt`, `log`, `expm` et `power` — [tutoriel, section 1.5](../../tutorial/fr/01_matrix_algebras.md#s-1-5).

<a id="r-1-7"></a>

**Corollaire 1.7.** $`A\ge0`$ si et seulement si $`A=T^*T`$ pour une certaine matrice $`T`$. (Si
$`A\ge0`$, on prend $`T=\sqrt A`$ ; réciproquement, $`\langle x,T^*Tx\rangle=\|Tx\|^2\ge0`$.)

**Dérivées.** Les méthodes de gradient dérivent $`f(A)`$ par rapport à $`A`$. Pour $`A`$ autoadjointe, la dérivée est un
produit coefficient par coefficient avec une matrice de quotients différentiels de $`f`$.

<a id="r-1-11"></a>

**Théorème 1.11 (Daleckii–Krein).** Soient $`A=U\,\mathrm{diag}(\lambda_1,\dots,\lambda_n)\,U^*`$ autoadjointe, $`I`$ un
intervalle ouvert contenant $`\mathrm{spec}\,A`$ et $`f:I\to\mathbb C`$ continûment dérivable. Soient
```math
f^{[1]}(x,y)=\frac{f(x)-f(y)}{x-y}\quad(x\ne y),\qquad f^{[1]}(x,x)=f'(x)
```
la **différence divisée** de $`f`$ et $`L_{ij}=f^{[1]}(\lambda_i,\lambda_j)`$ la **matrice de Loewner**. Alors, pour toute
$`E`$ autoadjointe,
```math
\frac{d}{dt}\Big|_{t=0}f(A+tE)=U\big(L\circ U^*EU\big)U^*,
```
où $`(X\circ Y)_{ij}=X_{ij}Y_{ij}`$ est le produit coefficient par coefficient (de Schur).

*Démonstration.* Par le théorème fondamental de l'analyse, $`f^{[1]}(x,y)=\int_0^1f'(y+s(x-y))\,ds`$ pour tous
$`x,y\in I`$ (pour $`x=y`$ les deux membres valent $`f'(x)`$), donc $`f^{[1]}`$ est continue sur $`I\times I`$. Réordonner
les colonnes de $`U`$ permute de la même façon les lignes et les colonnes de $`L`$ et de $`U^*EU`$ et ne change pas le
membre de droite ; on peut donc supposer $`\lambda_1\le\dots\le\lambda_n`$. Soit $`\delta>0`$ avec
$`[\lambda_1-\delta,\lambda_n+\delta]\subset I`$. Pour $`0<|t|\,\|E\|<\delta`$, écrivons $`A+tE=V\,\mathrm{diag}(\mu_1,\dots,\mu_n)\,V^*`$
avec $`\mu_1\le\dots\le\mu_n`$ ; alors $`|\mu_j-\lambda_j|\le|t|\,\|E\|`$ (théorème A.10), donc tous les $`\mu_j\in I`$. Posons
$`W=U^*V`$ (unitaire) et $`\tilde E=U^*EU`$. Comme $`U^*f(A)=\mathrm{diag}(f(\lambda_i))\,U^*`$ et $`f(A+tE)V=V\,\mathrm{diag}(f(\mu_j))`$,
```math
\big(U^*(f(A+tE)-f(A))V\big)_{ij}=\big(f(\mu_j)-f(\lambda_i)\big)W_{ij}=f^{[1]}(\lambda_i,\mu_j)\,(\mu_j-\lambda_i)\,W_{ij},
```
et le même calcul pour $`f(x)=x`$ donne $`(\mu_j-\lambda_i)W_{ij}=t\,(\tilde EW)_{ij}`$. Avec
$`(L_t)_{ij}=f^{[1]}(\lambda_i,\mu_j)`$, c'est l'identité exacte
```math
\frac1t\big(f(A+tE)-f(A)\big)=U\big((L_t\circ\tilde EW)\,W^*\big)U^*.
```
Quand $`t\to0`$, $`L_t\to L`$ coefficient par coefficient, car $`f^{[1]}`$ est continue et $`\mu_j\to\lambda_j`$. Soient
$`D_{ij}=W_{ij}`$ si $`\lambda_i=\lambda_j`$, $`D_{ij}=0`$ sinon, et $`R=W-D`$. Si $`\lambda_i\ne\lambda_j`$ et
$`|t|\,\|E\|\le|\lambda_i-\lambda_j|/2`$, alors $`|\mu_j-\lambda_i|\ge|\lambda_i-\lambda_j|/2`$ et
$`|(\tilde EW)_{ij}|\le\|\tilde EW\|=\|E\|`$, donc $`|W_{ij}|\le2|t|\,\|E\|/|\lambda_i-\lambda_j|`$. Ainsi $`R\to0`$ et
$`DD^*=(W-R)(W-R)^*\to WW^*=1`$. Comme $`L_{ij}=L_{ik}`$ dès que $`\lambda_k=\lambda_j`$, et $`D_{kj}=0`$ sauf si
$`\lambda_k=\lambda_j`$, on a $`(L\circ\tilde ED)_{ij}=\sum_kL_{ik}\tilde E_{ik}D_{kj}=((L\circ\tilde E)D)_{ij}`$. Par conséquent
```math
(L_t\circ\tilde EW)W^*-L\circ\tilde E=\big((L_t-L)\circ\tilde EW\big)W^*+(L\circ\tilde ER)W^*+(L\circ\tilde ED)R^*+(L\circ\tilde E)(DD^*-1),
```
et chaque terme tend vers $`0`$, car les coefficients de $`W`$, $`D`$ et $`\tilde E`$ sont bornés par $`1`$, $`1`$ et $`\|E\|`$. $`\square`$

Pour $`f`$ réelle, $`L`$ est réelle, et $`\mathrm{Tr}\big(G^*U(L\circ U^*EU)U^*\big)=\mathrm{Tr}(\Gamma^*E)`$ avec
$`\Gamma=U(L\circ U^*GU)U^*`$, pour toute $`G`$ (écrire les deux membres comme sommes sur $`i,j`$ dans la base $`U`$) : la
dérivée de $`A\mapsto\mathrm{Re}\,\mathrm{Tr}(G^*f(A))`$ dans la direction $`E`$ est $`\mathrm{Re}\,\mathrm{Tr}(\Gamma^*E)`$.

*Dans la bibliothèque :* `apply_function` est dérivée (par autograd) au moyen du théorème 1.11 ; la matrice de
Loewner est `functional.divided_differences`.

**Arrondi.** Pour $`\lambda_i,\lambda_j`$ proches, le quotient $`(f(\lambda_i)-f(\lambda_j))/(\lambda_i-\lambda_j)`$ perd en
précision par compensation, et pour $`\lambda_i=\lambda_j`$ il n'est pas défini. Les dérivées sont alors les meilleures
données.

<a id="r-1-12"></a>

**Proposition 1.12 (valeur au point milieu).** Soient $`a<b`$, $`h=b-a`$ et $`f`$ trois fois continûment dérivable sur
$`[a,b]`$ (valeurs complexes permises). Alors
```math
f^{[1]}(a,b)-\frac{f'(a)+f'(b)}2=-\frac{h^2}2\int_0^1s(1-s)\,f'''(a+sh)\,ds,
\qquad
\Big|f^{[1]}(a,b)-\frac{f'(a)+f'(b)}2\Big|\le\frac{h^2}{12}\max_{[a,b]}|f'''|,
```
avec égalité pour les polynômes de degré $`3`$.

*Démonstration.* Soient $`g(s)=f'(a+sh)`$, de sorte que $`f^{[1]}(a,b)=\frac1h\int_a^bf'=\int_0^1g`$, et $`\varphi(s)=s(1-s)/2`$,
de sorte que $`\varphi(0)=\varphi(1)=0`$, $`\varphi'(0)=\frac12=-\varphi'(1)`$ et $`\varphi''=-1`$. Deux intégrations par parties
donnent $`\int_0^1\varphi g''=-\int_0^1\varphi'g'=-[\varphi'g]_0^1+\int_0^1\varphi''g=\frac{g(0)+g(1)}2-\int_0^1g`$. Comme
$`g''(s)=h^2f'''(a+sh)`$, c'est l'identité ; $`\int_0^1\varphi=\frac1{12}`$ donne la borne, et pour un polynôme de degré
$`3`$, $`f'''`$ est constante. $`\square`$

<a id="r-1-13"></a>

**Proposition 1.13 (quotient ou point milieu).** Supposons une arithmétique flottante d'unité d'arrondi $`u\le1/20`$ :
pour des nombres flottants $`x,y`$ et $`\circ\in\{+,-,\times,/\}`$, le résultat calculé est $`(x\circ y)(1+\theta)`$ avec
$`|\theta|\le u`$ (ni dépassement de capacité ni sous-dépassement). Soient $`a\ne b`$ des nombres flottants réels,
$`h=|a-b|`$, $`x=f^{[1]}(a,b)`$, et $`\hat f_a,\hat f_b`$ des valeurs calculées avec $`|\hat f_a-f(a)|\le\delta_a`$ et
$`|\hat f_b-f(b)|\le\delta_b`$.

(a) Le quotient calculé $`\hat q`$ de $`\hat f_a-\hat f_b`$ par $`a-b`$ vérifie
$`|\hat q-x|\le R:=(\delta_a+\delta_b)/h+5u|\hat q|`$.

(b) Soit $`|\hat m-x|\le T`$ (par exemple $`\hat m`$ la valeur au point milieu calculée et $`T`$ la borne de la
proposition 1.12 plus son erreur d'arrondi). Posons $`y=\hat m`$ si $`|\hat q-\hat m|\le2R`$, et $`y=\hat q`$ sinon. Alors
$`|y-x|\le\min(T,3R)`$.

(c) Si $`\delta_a+\delta_b\le\Delta`$ et $`T=Mh^2+T_0`$, alors $`|y-x|\le(3\Delta)^{2/3}M^{1/3}+T_0+15u|\hat q|`$ pour tout $`h`$.
Les deux termes d'erreur $`3\Delta/h`$ et $`Mh^2`$ s'équilibrent en $`h_*=(3\Delta/M)^{1/3}`$.

*Démonstration.* (a) Soit $`\tilde q=(\hat f_a-\hat f_b)/(a-b)`$. Les trois opérations donnent
$`\hat q=\tilde q\,(1+\theta_1)(1+\theta_3)/(1+\theta_2)=\tilde q\,(1+\theta)`$ avec $`|\theta_i|\le u`$, d'où
$`-3u\le\frac{(1-u)^2}{1+u}-1\le\theta\le\frac{(1+u)^2}{1-u}-1=\frac{3u+u^2}{1-u}\le4u`$. Donc
$`|\hat q-\tilde q|=|\theta|\,|\tilde q|\le\frac{4u}{1-4u}|\hat q|\le5u|\hat q|`$, et
$`|\tilde q-x|=|(\hat f_a-f(a))-(\hat f_b-f(b))|/h\le(\delta_a+\delta_b)/h`$. (b) Si $`y=\hat m`$, alors $`|y-x|\le T`$ et
$`|y-x|\le|\hat m-\hat q|+|\hat q-x|\le3R`$. Si $`y=\hat q`$, alors $`|y-x|\le R`$, et
$`T\ge|\hat m-x|\ge|\hat m-\hat q|-|\hat q-x|>2R-R=R`$. (c) $`\min(T,3R)\le\min(Mh^2,3\Delta/h)+T_0+15u|\hat q|`$ ; pour
$`h\ge h_*`$, $`3\Delta/h\le3\Delta/h_*`$, pour $`h<h_*`$, $`Mh^2<Mh_*^2`$, et $`3\Delta/h_*=Mh_*^2=(3\Delta)^{2/3}M^{1/3}`$. $`\square`$

La démonstration vaut mot pour mot pour $`f`$ à valeurs complexes (parties réelle et imaginaire sont arrondies
séparément, et $`a-b`$ est réel ; $`\theta`$ est alors complexe avec la même borne). La règle (b) a besoin de $`R`$, qui
se calcule à partir des données, mais pas de $`T`$ : sans rien savoir de $`f'''`$, son erreur est au plus trois fois la
plus petite des deux bornes. Un seuil fixe, valeur au point milieu si et seulement si $`h\le h_*`$, demanderait
$`M=\max|f'''|/12`$. Le seuil $`h\le\varepsilon^{1/3}\max(|a|,|b|)`$ ($`\varepsilon=2u`$) des versions antérieures de la
bibliothèque est de ce type, avec $`M`$ de l'ordre de $`\Delta/\max(|a|,|b|)^3`$ : juste pour $`x^p`$ ou $`\log x`$, dont les
dérivées suivent l'échelle de $`|x|`$, faux pour les fonctions qui varient à une échelle sans rapport avec $`|x|`$. Pour
$`f(x)=e^{x-c}`$ en $`a\approx c=10^4`$, il prenait la valeur au point milieu jusqu'à $`h\approx0.06`$, avec une erreur
relative atteignant $`h^2/12\approx3\cdot10^{-4}`$ (la proposition 1.12 y est une égalité à un facteur $`e^{h}`$ près) ; il
en va de même pour $`e^{-i\tau x}`$ et un hamiltonien comportant un grand décalage constant.

*Dans la bibliothèque :* `functional.divided_differences` applique la règle (b) à chaque paire de valeurs propres,
avec $`\delta_a=\varepsilon\,(|\hat f_a|+|a\hat f'_a|)`$ et $`4\varepsilon|\hat q|`$ au lieu de $`5u|\hat q|`$. Le choix de $`\delta_a`$
est une hypothèse sur l'évaluation de $`f`$, non un théorème : erreur relative au plus $`\varepsilon`$, en un argument
perturbé d'au plus $`\varepsilon|a|`$ (au premier ordre ; le second terme couvre des fonctions comme $`e^{-i\tau x}`$, qui
arrondissent $`\tau x`$). La même règle s'applique aux valeurs propres complexes des matrices normales, que la
démonstration ne couvre pas. La valeur au point milieu est aussi prise pour les égalités et pour les paires que le
solveur de valeurs propres ne sépare pas, $`h\le n\varepsilon\max_l|\lambda_l|`$ : son erreur y est au plus
$`(n\varepsilon\max_l|\lambda_l|)^2\max|f'''|/12`$ par la proposition 1.12, tandis que le quotient n'aurait pas de sens pour
une $`f`$ discontinue, comme une fonction indicatrice.

*Calcul :* le théorème 1.11 comparé à un quotient différentiel en un spectre dégénéré, et les propositions 1.12 et
1.13 pour $`e^{x-c}`$ avec $`c=10^4`$ — [tutoriel, section 1.5](../../tutorial/fr/01_matrix_algebras.md#s-1-5).

<a id="s-1-4"></a>

## 1.4 Décomposition polaire et valeurs singulières

Pour toute matrice $`A`$, la matrice $`A^*A`$ est positive (corollaire 1.7) ; sa racine carrée
$`|A|=(A^*A)^{1/2}`$ est la **valeur absolue** (ou module) de $`A`$.

<a id="r-1-8"></a>

**Théorème 1.8 (décomposition polaire).** Toute matrice $`A\in M_n(\mathbb C)`$ s'écrit
$`A=W|A|`$ avec $`W`$ unitaire.

*Démonstration.* Pour tout $`x`$, $`\||A|x\|^2=\langle x,|A|^2x\rangle=\langle x,A^*Ax\rangle=\|Ax\|^2`$.
Par conséquent $`W_0:|A|x\mapsto Ax`$ est une isométrie linéaire bien définie de l'image
$`R_1`$ de $`|A|`$ sur l'image $`R_2`$ de $`A`$ : si $`|A|x=|A|x'`$, alors
$`\|A(x-x')\|=\||A|(x-x')\|=0`$. En particulier $`\dim R_1=\dim R_2`$, de sorte que les
supplémentaires orthogonaux ont eux aussi même dimension ; soit $`W_1:R_1^\perp\to R_2^\perp`$
une isométrie quelconque entre eux (qui envoie une base orthonormée sur une base
orthonormée). Posons $`W=W_0\oplus W_1`$, c'est-à-dire $`W(u_1+u_2)=W_0u_1+W_1u_2`$ pour $`u_1\in R_1`$,
$`u_2\in R_1^\perp`$. Comme $`W_0u_1\in R_2\perp R_2^\perp\ni W_1u_2`$, le théorème de Pythagore donne
$`\|Wu\|^2=\|u_1\|^2+\|u_2\|^2=\|u\|^2`$ pour tout $`u`$. Par conséquent $`\langle u,(W^*W-1)u\rangle=0`$
pour tout $`u`$, donc $`W^*W=1`$ par polarisation (démonstration du corollaire 1.4(3)) ; alors $`W`$ est
inversible et $`W^{-1}=W^*`$, donc aussi $`WW^*=1`$ : $`W`$ est unitaire, et
$`W|A|x=W_0|A|x=Ax`$. $`\square`$

*Calcul :* le théorème 1.8 sur une matrice de Ginibre, avec $`|A|`$ calculé par `abs` — [tutoriel, section 1.6](../../tutorial/fr/01_matrix_algebras.md#s-1-6).

Les valeurs propres $`s_1\ge\dots\ge s_n\ge0`$ de $`|A|`$ sont les **valeurs singulières** de
$`A`$. Trois normes seront utilisées tout au long du cours :
```math
\|A\|=s_1,\qquad \|A\|_2=\Big(\sum_is_i^2\Big)^{1/2}=(\mathrm{Tr}\,A^*A)^{1/2},
\qquad\|A\|_1=\sum_is_i=\mathrm{Tr}|A|.
```
La première identité est vraie car $`\|Ax\|=\||A|x\|`$ pour tout $`x`$, donc
$`\|A\|=\||A|\|=s_1`$ d'après le corollaire 1.4(1) ; la deuxième car
$`\mathrm{Tr}\,A^*A=\mathrm{Tr}|A|^2=\sum_is_i^2`$.

<a id="r-1-9"></a>

**Théorème 1.9 (dualité de la norme trace).**
$`\|A\|_1=\max_{U}|\mathrm{Tr}(AU)|`$ (maximum sur les matrices
unitaires $`U`$), et ce maximum est atteint en $`U=W^*`$, où $`A=W|A|`$. En particulier
$`|\mathrm{Tr}\,A|\le\|A\|_1`$.

*Démonstration.* $`\mathrm{Tr}(AW^*)=\mathrm{Tr}(W|A|W^*)=\mathrm{Tr}|A|`$. Pour
toute matrice unitaire $`U`$, par cyclicité de la trace et par l'inégalité de Cauchy–Schwarz
pour $`\langle X,Y\rangle_{HS}=\mathrm{Tr}(X^*Y)`$,
```math
|\mathrm{Tr}(AU)|=|\mathrm{Tr}(|A|^{1/2}\,UW\,|A|^{1/2})|
=|\langle|A|^{1/2},\,UW|A|^{1/2}\rangle_{HS}|
\le\||A|^{1/2}\|_2\,\|UW|A|^{1/2}\|_2=\mathrm{Tr}|A|,
```
puisque $`\|UW X\|_2=\|X\|_2`$ pour $`UW`$ unitaire et $`\||A|^{1/2}\|_2^2=\mathrm{Tr}|A|`$. $`\square`$

En particulier, les trois quantités ci-dessus sont bien des normes. $`\|\cdot\|`$ est une
norme d'opérateur ; $`\|\cdot\|_2`$ est la norme associée au produit scalaire
$`\langle X,Y\rangle_{HS}`$ ; quant à $`\|\cdot\|_1`$, le théorème 1.9 donne
$`\|A+B\|_1=\max_U|\mathrm{Tr}(AU)+\mathrm{Tr}(BU)|\le\|A\|_1+\|B\|_1`$ et
$`\|cA\|_1=\max_U|c\mathrm{Tr}(AU)|=|c|\,\|A\|_1`$, tandis que
$`\|A\|_1=0`$ impose $`|A|=0`$ et $`A=W|A|=0`$. Pour $`P\ge0`$, on a $`|P|=(P^2)^{1/2}=P`$ d'après le
théorème 1.6, donc $`\|P\|_1=\mathrm{Tr}\,P`$.

*Calcul :* les trois normes (`operator_norm`, `frobenius_norm`, `trace_norm`) et le théorème 1.9 — [tutoriel, section 1.6](../../tutorial/fr/01_matrix_algebras.md#s-1-6).

<a id="s-1-5"></a>

## 1.5 Les opérateurs dans la bibliothèque

*Dans la bibliothèque :* l'algèbre $`M_n(\mathbb C)`$ sous la forme `TypeIAlgebra([n], [n])`, la façon dont ses opérateurs sont stockés, et le calcul fonctionnel `apply_function` — [tutoriel, section 1.5](../../tutorial/fr/01_matrix_algebras.md#s-1-5).

<a id="s-1-6"></a>

## 1.6 Normes et dualité de la norme trace

*Calcul :* les normes, la décomposition polaire et la dualité de la section 1.4, sur une matrice de Ginibre — [tutoriel, section 1.6](../../tutorial/fr/01_matrix_algebras.md#s-1-6).

<a id="s-1-7"></a>

## 1.7 Lots : de nombreux opérateurs à la fois

À titre d'exemple d'une quantité que l'on estime à l'aide de nombreux opérateurs aléatoires à la fois,
calculons un moment exact de l'**ensemble gaussien unitaire** (GUE).

Soit $`G`$ une matrice à coefficients indépendants, dont les parties réelles et imaginaires
sont indépendantes et de loi $`N(0,\tfrac12)`$, de sorte que $`\mathbb E|G_{ij}|^2=1`$, et soit
$`H=(G+G^*)/\sqrt2`$.

<a id="r-1-10"></a>

**Proposition 1.10.** $`\mathbb E\,\mathrm{Tr}\,H^4=2n^3+n`$. Par conséquent, pour
$`X=H/\sqrt n`$, $`\ \mathbb E\,\tfrac1n\mathrm{Tr}\,X^4=2+\tfrac1{n^2}`$.

*Démonstration.* *Formule d'Isserlis (ou de Wick)* (théorème A.24, cas $`m=4`$). Les parties réelles
et imaginaires des coefficients de $`G`$ sont des variables gaussiennes centrées indépendantes ; ce
sont donc les coordonnées d'un même vecteur gaussien centré, et les coefficients de $`H`$ et de
$`\bar H`$ en sont des combinaisons linéaires complexes. Le théorème A.24 donne alors, pour quatre
quelconques de ces variables $`Z_1,\dots,Z_4`$ (les répétitions sont permises),
```math
\mathbb E[Z_1Z_2Z_3Z_4]=\mathbb E[Z_1Z_2]\,\mathbb E[Z_3Z_4]+\mathbb E[Z_1Z_3]\,\mathbb E[Z_2Z_4]+\mathbb E[Z_1Z_4]\,\mathbb E[Z_2Z_3],
```
la somme, sur les trois appariements de $`\{1,2,3,4\}`$, des produits des espérances des paires
(sans conjugaison complexe). La démonstration de l'annexe A compare la série entière de la
fonction génératrice des moments aux moments gaussiens et ne dérive pas sous l'espérance.

*Le moment.* Les espérances des paires valent
$`\mathbb E[H_{ij}H_{kl}]=\tfrac12\mathbb E[(G_{ij}+\bar G_{ji})(G_{kl}+\bar G_{lk})]=\delta_{il}\delta_{jk}`$,
car $`\mathbb E[G_{ij}G_{kl}]=0`$ et $`\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}`$ (proposition A.16(a) ;
pour $`(i,j)\neq(k,l)`$ les coefficients sont indépendants et centrés, et pour $`G_{ij}=x+iy`$ on a
$`\mathbb E[G_{ij}^2]=\mathbb E[x^2-y^2]+2i\,\mathbb E[xy]=0`$ et $`\mathbb E|G_{ij}|^2=\tfrac12+\tfrac12=1`$).
Dans $`\mathrm{Tr}\,H^4=\sum_{ijkl}H_{ij}H_{jk}H_{kl}H_{li}`$, les appariements donnent :
$`\{12\}\{34\}`$ : $`\delta_{ik}\delta_{ki}`$, de somme $`n^3`$ (indices libres $`i=k`$, $`j`$, $`l`$) ;
$`\{14\}\{23\}`$ : $`\delta_{jl}\delta_{jl}`$, de somme $`n^3`$ ;
$`\{13\}\{24\}`$ : $`\delta_{il}\delta_{jk}\delta_{ji}\delta_{kl}`$, qui impose $`i=j=k=l`$, de somme $`n`$.
En divisant par $`n\cdot n^2`$, on obtient le second énoncé. $`\square`$

Les deux appariements *non croisés* $`\{12\}\{34\}`$ et $`\{14\}\{23\}`$ (les paires peuvent être
dessinées comme des arcs sans intersection au-dessus des points $`1,2,3,4`$) donnent le terme
principal $`2`$, pour tout $`n`$ ; l'appariement croisé $`\{13\}\{24\}`$ donne la correction en $`1/n^2`$.

*Calcul :* la proposition 1.10 par Monte-Carlo, sur un lot de $`200\,000`$ échantillons — [tutoriel, section 1.7](../../tutorial/fr/01_matrix_algebras.md#s-1-7).

<a id="s-1-8"></a>

## 1.8 Points de comparaison

*Dans la bibliothèque :* les implémentations de référence auxquelles se comparer — [tutoriel, section 1.8](../../tutorial/fr/01_matrix_algebras.md#s-1-8).

<a id="x-1"></a>

## Exercices

1. Montrer que $`(fg)(A)=f(A)g(A)`$ pour $`A`$ normale et des fonctions $`f,g`$ sur $`\mathrm{spec}\,A`$.
2. Montrer que $`\|A\|_1=\|A^*\|_1`$ et $`\|AB\|_1\le\|A\|\,\|B\|_1`$ à l'aide du théorème 1.9.
3. La **moyenne géométrique** de $`A,B>0`$ est $`A\#B=A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}`$.
   Démontrer que $`A\#B`$ est l'unique solution positive $`X`$ de $`XA^{-1}X=B`$, et en déduire
   $`A\#B=B\#A`$.

<details>
<summary>Solutions</summary>

1. D'après le théorème 1.5(c).

2. $`|\mathrm{Tr}(A^*U)|=|\overline{\mathrm{Tr}(U^*A)}|=|\mathrm{Tr}(AU^*)|`$, et
   $`U\mapsto U^*`$ est une bijection du groupe unitaire, donc les maxima du théorème 1.9
   coïncident. Pour la seconde assertion, écrivons $`B=V|B|`$. Pour $`U`$ unitaire,
   $`|\mathrm{Tr}(ABU)|=|\mathrm{Tr}(|B|^{1/2}UAV|B|^{1/2})|    \le\||B|^{1/2}\|_2\,\|UAV|B|^{1/2}\|_2\le\|A\|\mathrm{Tr}|B|`$, en utilisant
   $`\|YX\|_2\le\|Y\|\,\|X\|_2`$ (appliquer $`\|Yx\|\le\|Y\|\,\|x\|`$ à chaque colonne de $`X`$) et
   $`\|UAV\|=\|A\|`$ pour $`U,V`$ unitaires (elles conservent la norme des vecteurs).

3. Posons $`T=A^{-1/2}XA^{-1/2}`$. Alors $`XA^{-1}X=B`$ si et seulement si
   $`T^2=A^{-1/2}BA^{-1/2}`$, et $`X\ge0`$ si et seulement si $`T\ge0`$, puisque
   $`\langle x,Xx\rangle=\langle A^{1/2}x,TA^{1/2}x\rangle`$ et que $`A^{1/2}`$ est inversible (ici
   $`A^{\pm1/2}`$ sont données par le calcul fonctionnel avec $`t^{\pm1/2}`$ sur
   $`\mathrm{spec}\,A\subset(0,\infty)`$, donc $`A^{-1/2}A^{1/2}=1`$ d'après le théorème 1.5(c)). De
   plus $`A^{-1/2}BA^{-1/2}\ge0`$, puisque
   $`\langle x,A^{-1/2}BA^{-1/2}x\rangle=\langle A^{-1/2}x,BA^{-1/2}x\rangle\ge0`$.
   D'après le théorème 1.6, $`T=(A^{-1/2}BA^{-1/2})^{1/2}`$ est l'unique solution positive,
   c'est-à-dire $`X=A\#B`$. Toute solution $`X`$ de $`XA^{-1}X=B`$ est inversible, car
   $`\det(X)^2\det(A)^{-1}=\det B\neq0`$. Ainsi $`XA^{-1}X=B`$ équivaut à $`XB^{-1}X=A`$ (à partir
   de $`A^{-1}=X^{-1}BX^{-1}`$, on inverse les deux membres), et l'unicité donne $`A\#B=B\#A`$.

</details>

*Calcul :* les exercices 1 et 3, numériquement — [tutoriel, Exercices](../../tutorial/fr/01_matrix_algebras.md#x-1).

**Suite :** le [chapitre 2](02_von_neumann_algebras.md) remplace $`M_n`$ par une somme directe
d'algèbres de matrices et démontre que toute algèbre de von Neumann de dimension finie est
de cette forme.
