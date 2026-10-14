[English](../13_jordan_algebras.md) · [Español](../es/13_jordan_algebras.md) · **Français** · [Deutsch](../de/13_jordan_algebras.md) · [中文](../zh/13_jordan_algebras.md) · [日本語](../ja/13_jordan_algebras.md) · [Русский](../ru/13_jordan_algebras.md)

# 13. Algèbres de Jordan : pourquoi la mécanique quantique utilise des *-algèbres

Le produit $`xy`$ de deux observables n'est pas une observable, sauf si $`x`$ et $`y`$ commutent, mais le produit symétrisé
$`x\circ y=\frac12(xy+yx)`$ en est toujours une. En 1933, P. Jordan proposa de prendre les observables munies de ce
produit comme structure de base de la mécanique quantique, et Jordan, von Neumann et Wigner classifièrent les algèbres
de dimension finie qui apparaissent. Ce chapitre développe la théorie en dimension finie, avec des démonstrations
complètes ici ou, pour les quatre longues (la classification, le théorème de Koecher–Vinberg, le théorème de Kadison et
les automorphismes de Jordan de $`\bigoplus_cM_{k_c}`$), dans les annexes E et F, et se demande ce qui échappe à la
structure de Jordan. La réponse est le signe du commutateur : l'algèbre de Jordan d'un système ne détermine le
commutateur qu'au signe près (section 13.5), mais les systèmes composés ont besoin du signe : c'est lui qui décide
quels opérateurs conjoints sont positifs, et le test PPT construit à partir de lui détecte l'intrication
(proposition 13.12).

**Vous apprendrez** le produit de Jordan et la représentation quadratique $`P(y)x=yxy`$ qui, contrairement à $`\circ`$,
préserve la positivité ; les facteurs de spin, leur théorie spectrale et le cône de Lorentz ; des sous-algèbres de
Jordan de matrices hermitiennes qui ne sont pas des parties autoadjointes de $`*`$-algèbres ; la classification de
Jordan–von Neumann–Wigner ; pourquoi la transposée est un automorphisme de Jordan et non un $`*`$-automorphisme, et ce
qu'elle fait sur un système composé ; et une inégalité de traces du chapitre 11 écrite dans le langage de Jordan et
démontrée dans un facteur de spin.

Les vérifications numériques de ce chapitre se trouvent dans le [Tutoriel 13](../../tutorial/fr/13_jordan_algebras.md).

Dans tout ce chapitre, $`\mathrm{Herm}_n(\mathbb K)`$ est l'espace vectoriel réel des matrices hermitiennes
$`n\times n`$ sur $`\mathbb K=\mathbb R`$ ou $`\mathbb C`$ (la section 13.4 utilise aussi les quaternions $`\mathbb H`$ et les
octonions $`\mathbb O`$, construits à la section E.6, et $`\mathrm{Herm}_n`$ sur ceux-ci, section E.7), $`M_{sa}`$ la
partie autoadjointe de
$`\mathcal M=\bigoplus_cM_{k_c}(\mathbb C)`$, $`X,Y,Z`$ les matrices de Pauli et $`\sigma=(X,Y,Z)`$, de sorte que
$`u\cdot\sigma=u_1X+u_2Y+u_3Z`$ pour $`u\in\mathbb R^3`$. $`[a,b]=ab-ba`$, et $`|u|`$ est la norme euclidienne de
$`u\in\mathbb R^m`$.

<a id="s-13-1"></a>

## 13.1 Observables et produit de Jordan

**Définition.** Une **algèbre de Jordan** est un espace vectoriel réel $`J`$ muni d'un produit bilinéaire et commutatif
$`x\circ y`$ qui vérifie l'**identité de Jordan** $`(x^2\circ y)\circ x=x^2\circ(y\circ x)`$, où $`x^2=x\circ x`$. Elle est
**spéciale** si elle est isomorphe à un sous-espace d'une algèbre associative stable par $`\frac12(xy+yx)`$. Sur $`J`$,
$`L(y)x=y\circ x`$, et la **représentation quadratique** est $`P(y)=2L(y)^2-L(y^2)`$.

<a id="r-13-1"></a>

**Lemme 13.1.** Soit $`\mathcal A`$ une algèbre associative réelle ou complexe et $`x\circ y=\frac12(xy+yx)`$. Alors
$`\circ`$ vérifie l'identité de Jordan, et les puissances de Jordan $`x^{\circ1}=x`$, $`x^{\circ(k+1)}=x\circ x^{\circ k}`$
sont égales aux puissances associatives $`x^k`$. En particulier, $`M_{sa}`$, $`\mathrm{Herm}_n(\mathbb C)`$ et
$`\mathrm{Herm}_n(\mathbb R)`$ sont des algèbres de Jordan spéciales.

*Démonstration.* En développant, $`4(x^2\circ y)\circ x=x^2yx+yx^3+x^3y+xyx^2`$ et
$`4x^2\circ(y\circ x)=x^2yx+x^3y+yx^3+xyx^2`$, les mêmes quatre termes. Par récurrence,
$`x\circ x^k=\frac12(x^{k+1}+x^{k+1})=x^{k+1}`$. Si $`x,y`$ sont hermitiennes, $`(xy+yx)^*=yx+xy`$, donc les trois espaces
sont stables par $`\circ`$ ; pour $`x,y`$ symétriques réelles, le produit est réel. $`\square`$

*Calcul :* l'identité de Jordan et l'associativité des puissances dans $`M_2\oplus M_3`$ — [tutoriel, section 13.1](../../tutorial/fr/13_jordan_algebras.md#s-13-1).

<a id="r-13-2"></a>

**Proposition 13.2.** Soient $`x,y,z\in M_{sa}`$.

1. $`P(y)x=yxy`$, et la **formule fondamentale** $`P(P(y)x)=P(y)P(x)P(y)`$ est vérifiée.
2. $`\mathrm{Tr}\,P(y)x=\mathrm{Tr}\,y^2x`$.
3. Si $`x\ge0`$, alors $`P(y)x\ge0`$. Si $`y`$ est inversible, $`P(y)`$ est une bijection du cône positif sur lui-même,
   d'inverse $`P(y^{-1})`$.
4. Le produit de Jordan d'éléments positifs n'est pas nécessairement positif : pour $`p=\begin{pmatrix}1&0\\0&0\end{pmatrix}`$ et
   $`q=\frac12\begin{pmatrix}1&1\\1&1\end{pmatrix}`$, $`\det(p\circ q)=-\frac1{16}`$.

*Démonstration.* (1) $`2y\circ(y\circ x)=\frac12(y^2x+2yxy+xy^2)`$ et $`y^2\circ x=\frac12(y^2x+xy^2)`$ ; on soustrait.
Ensuite, $`P(P(y)x)z=yxy\,z\,yxy=P(y)P(x)P(y)z`$ par associativité. (2) Cyclicité de la trace. (3)
$`yxy=(x^{1/2}y)^*(x^{1/2}y)\ge0`$ (corollaire 1.7), et $`P(y)P(y^{-1})=P(y^{-1})P(y)=\mathrm{id}`$. (4)
$`p\circ q=\frac14\begin{pmatrix}2&1\\1&0\end{pmatrix}`$. $`\square`$

*Calcul :* la proposition 13.2 dans $`M_2\oplus M_3`$, et l'exemple (4) — [tutoriel, section 13.1](../../tutorial/fr/13_jordan_algebras.md#s-13-1). *Dans la bibliothèque :* `jordan.jordan_product`, `jordan.jordan_power` et `jordan.quadratic_representation`, par lots sur les secteurs d'une `TypeIAlgebra` — [tutoriel, section 13.1](../../tutorial/fr/13_jordan_algebras.md#s-13-1).

En écrivant $`P(y)`$ comme $`2L(y)^2-L(y^2)`$, chaque identité de (1) s'exprime avec $`\circ`$ seul. La démonstration de (1)
n'utilise que l'associativité du produit matriciel, donc dans toute algèbre de Jordan spéciale $`J\subseteq\mathcal A`$
($`\mathcal A`$ associative) elle donne $`P(y)x=yxy`$ et la formule fondamentale $`P(P(y)x)=P(y)P(x)P(y)`$ sur $`J`$. (Cela
couvre toutes les algèbres dans lesquelles le cours l'utilise ; la formule fondamentale n'est pas utilisée pour des
algèbres non spéciales.)

Parmi les produits $`\alpha xy+\beta yx`$ ($`\alpha,\beta\in\mathbb R`$) sur $`\mathrm{Herm}_n(\mathbb C)`$, $`n\ge2`$,
$`\frac12(xy+yx)`$ est le seul qui soit symétrique en $`x,y`$ et égal à $`x^2`$ pour $`y=x`$ : la symétrie signifie
$`(\alpha-\beta)(xy-yx)=0`$ pour tous $`x,y`$, et $`XZ\ne ZX`$ (dans le bloc $`2\times2`$ en haut à gauche) impose
$`\alpha=\beta`$ ; ensuite $`(\alpha+\beta)1=1`$ pour $`x=y=1`$ donne $`\alpha=\beta=\frac12`$.

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). $`x\circ y`$ est le produit
> symétrisé d'observables (la seule combinaison de ce type, comme on l'a vu ci-dessus), et c'est de nouveau une
> observable. La représentation quadratique est la **règle de Lüders** : mesurer une question oui/non $`p`$ et garder le
> résultat « oui » transforme une densité $`\rho`$ en $`P(p)\rho/\mathrm{Tr}\,P(p)\rho=p\rho p/\mathrm{Tr}(p\rho)`$.
> La proposition 13.2(4) dit que $`\circ`$ lui-même n'a pas de telle signification opérationnelle pour les éléments
> positifs.

<a id="s-13-2"></a>

## 13.2 Facteurs de spin et cône de Lorentz

**Définition.** Pour $`m\ge1`$, le **facteur de spin** $`V_m=\mathbb R\oplus\mathbb R^m`$ a pour éléments les $`x=(a,u)`$,
pour produit $`(a,u)\circ(b,v)=(ab+u\cdot v,\;av+bu)`$ et pour unité $`1=(1,0)`$. Sa **trace** est
$`\mathrm{tr}\,x=2a`$.

Le produit est bilinéaire et commutatif ; l'identité de Jordan découle d'une représentation matricielle.

<a id="r-13-3"></a>

**Lemme 13.3 (symétries anticommutantes).** Soient $`p=\lfloor m/2\rfloor`$ et $`D=2^p`$ ($`D=2`$ si $`m=1`$). Sur
$`(\mathbb C^2)^{\otimes p}`$, posons $`s_{2j-1}=Z^{\otimes(j-1)}\otimes X\otimes1`$, $`s_{2j}=Z^{\otimes(j-1)}\otimes Y\otimes1`$
pour $`j\le p`$, et $`s_m=Z^{\otimes p}`$ si $`m`$ est impair ($`s_1=Z`$ si $`m=1`$). Alors $`s_i=s_i^*`$ et
$`s_is_j+s_js_i=2\delta_{ij}1`$. L'application $`\phi(a,u)=a1+\sum_iu_is_i`$ est injective et vérifie
$`\phi(x\circ y)=\phi(x)\circ\phi(y)`$ ; par conséquent, $`V_m`$ est une algèbre de Jordan spéciale. Pour $`m\ge2`$, les
$`s_i`$ engendrent $`M_D(\mathbb C)`$ en tant que $`*`$-algèbre.

*Démonstration.* Les matrices de Pauli sont de carré $`1`$ et anticommutent deux à deux ; il en va de même de $`Z`$ avec
$`X`$ et $`Y`$. Deux chaînes distinctes $`s_i,s_j`$ diffèrent en exactement un site, de façon anticommutante : pour
$`s_{2j-1},s_{2j}`$ au site $`j`$ ($`X`$ contre $`Y`$) ; pour des chaînes d'indices $`j<j'`$, au site $`j`$ ($`X`$ ou $`Y`$ contre
$`Z`$), tandis qu'aux autres sites elles portent des facteurs égaux ou un facteur $`1`$ ; il en va de même pour
$`s_m=Z^{\otimes p}`$. Donc $`s_is_j=-s_js_i`$ pour $`i\ne j`$, et $`s_i^2=1`$. Ainsi
$`(u\cdot s)\circ(v\cdot s)=\sum_{ij}u_iv_j\,s_i\circ s_j=(u\cdot v)1`$, ce qui donne $`\phi(x)\circ\phi(y)=\phi(x\circ y)`$.
Chaque $`s_i`$ est de trace nulle (produit tensoriel comportant un facteur de trace nulle) et
$`\mathrm{Tr}(s_is_j)=\mathrm{Tr}(s_i\circ s_j)=D\delta_{ij}`$, de sorte que $`1,s_1,\dots,s_m`$ sont orthogonaux
pour $`\mathrm{Tr}(ab)`$ et $`\phi`$ est injective. L'identité de Jordan est vraie dans $`\mathrm{Herm}_D`$
(lemme 13.1) et se transporte par $`\phi`$. Enfin, $`s_{2j-1}s_{2j}=1\otimes\cdots\otimes XY\otimes\cdots\otimes1=iZ_j`$
($`Z`$ au site $`j`$), donc la $`*`$-algèbre engendrée contient $`Z_j`$, puis $`X_j=Z_1\cdots Z_{j-1}s_{2j-1}`$ et $`Y_j`$, donc
tous les produits tensoriels de matrices de Pauli, qui forment une base de $`M_D(\mathbb C)`$. $`\square`$

*Dans la bibliothèque :* les $`s_i`$ sont `jordan.anticommuting_symmetries(m)`, et $`V_m`$ avec l'application $`\phi`$ est la classe `SpinFactor` (`to_matrix`) — [tutoriel, section 13.2](../../tutorial/fr/13_jordan_algebras.md#s-13-2).

Au signe près, $`s_{2j-1}`$ et $`s_{2j}`$ sont les opérateurs de Majorana $`c_j+c_j^*`$ et $`i(c_j-c_j^*)`$ des fermions de
Jordan–Wigner $`c_j`$ de la section 9.4, le site $`j-1`$ du chapitre 9 étant le $`j`$-ième facteur tensoriel : là,
$`s^a=\frac12\sigma_a`$, $`n=s^z+\frac12`$ et $`s^\pm=s^x\pm is^y`$, donc $`(-1)^n=1-2n=-Z`$ et $`s^-=\frac12(X-iY)`$, d'où
$`c_j=(-1)^{j-1}Z^{\otimes(j-1)}\otimes s^-\otimes1`$. Comme $`s^-+s^+=X`$ et $`i(s^--s^+)=i(-iY)=Y`$,
```math
s_{2j-1}=(-1)^{j-1}(c_j+c_j^*),\qquad s_{2j}=(-1)^{j-1}\,i(c_j-c_j^*).
```

*Calcul :* ces identités pour $`p=3`$ — [tutoriel, section 13.2](../../tutorial/fr/13_jordan_algebras.md#s-13-2).

<a id="r-13-4"></a>

**Théorème 13.4 (théorie spectrale de $`V_m`$).** Soient $`x=(a,u)\in V_m`$ et $`\lambda_\pm=a\pm|u|`$.

1. Les idempotents de $`V_m`$ sont $`0`$, $`1`$ et $`c_w=\frac12(1,w)`$ avec $`|w|=1`$ ; et $`c_w\circ c_{w'}=0`$ si et seulement
   si $`w'=-w`$.
2. Si $`u\ne0`$, $`x=\lambda_+e_++\lambda_-e_-`$ avec $`e_\pm=c_{\pm u/|u|}`$ ; si $`u=0`$, $`x=a1`$. Si $`x=\mu_1c_w+\mu_2c_{-w}`$,
   alors $`\{\mu_1,\mu_2\}=\{\lambda_+,\lambda_-\}`$. On appelle $`\lambda_\pm`$ les **valeurs propres** de $`x`$.
3. Pour $`f:\{\lambda_+,\lambda_-\}\to\mathbb R`$, posons $`f(x)=f(\lambda_+)e_++f(\lambda_-)e_-`$ ($`f(a)1`$ si $`u=0`$).
   Alors $`x^{\circ k}=\lambda_+^ke_++\lambda_-^ke_-`$, $`\phi(f(x))=f(\phi(x))`$, les valeurs propres de $`\phi(x)`$ sont
   $`\lambda_\pm`$, chacune de multiplicité $`D/2`$, et $`\mathrm{Tr}\,\phi(x)=\frac D2\mathrm{tr}\,x`$. De plus,
   $`\mathrm{tr}\,x=\lambda_++\lambda_-`$ et $`\det x:=\lambda_+\lambda_-=a^2-|u|^2`$.
4. Le cône des carrés $`Q=\{x^2:x\in V_m\}`$ est le **cône de Lorentz** $`\{(a,u):a\ge|u|\}=\{x:\lambda_-\ge0\}=\{x:\phi(x)\ge0\}`$.
   Il est **autodual** : $`y\in Q`$ si et seulement si $`\mathrm{tr}(x\circ y)\ge0`$ pour tout $`x\in Q`$.
5. $`\phi(P(y)x)=\phi(y)\phi(x)\phi(y)`$ ; par conséquent, $`P(y)Q\subseteq Q`$ pour tout $`y\in V_m`$.

*Démonstration.* (1) $`(a,u)^2=(a^2+|u|^2,2au)`$ est égal à $`(a,u)`$ si et seulement si $`u=0`$, $`a\in\{0,1\}`$, ou
$`a=\frac12`$, $`|u|=\frac12`$. $`c_w\circ c_{w'}=\frac14(1+w\cdot w',\,w+w')`$ s'annule si et seulement si $`w'=-w`$. (2)
$`\lambda_+e_++\lambda_-e_-=(a,0)+|u|(0,u/|u|)=x`$, et $`\mu_1c_w+\mu_2c_{-w}=(\frac{\mu_1+\mu_2}2,\frac{\mu_1-\mu_2}2w)`$, donc
$`\mu_1+\mu_2=2a`$ et $`|\mu_1-\mu_2|=2|u|`$. (3) D'après (1), $`(\alpha e_++\beta e_-)\circ(\gamma e_++\delta e_-)=\alpha\gamma e_++\beta\delta e_-`$, ce qui donne les puissances par récurrence. Avec $`w=u/|u|`$, $`(w\cdot s)^2=1`$, donc
$`\phi(e_\pm)=\frac12(1\pm w\cdot s)`$ sont des projecteurs orthogonaux de somme $`1`$ et de produit $`0`$, chacun de trace
$`D/2`$ ; ainsi $`\phi(x)=\lambda_+\phi(e_+)+\lambda_-\phi(e_-)`$ est une décomposition spectrale, et
$`f(\phi(x))=f(\lambda_+)\phi(e_+)+f(\lambda_-)\phi(e_-)=\phi(f(x))`$ (théorème 1.5). (4) $`x^2`$ a pour valeurs propres
$`\lambda_\pm^2\ge0`$ ; réciproquement, si $`\lambda_-\ge0`$, alors $`x=(\sqrt x)^2`$. $`\lambda_-\ge0`$ si et seulement si
$`a\ge|u|`$, si et seulement si $`\phi(x)\ge0`$ d'après (3). Pour l'autodualité, $`\mathrm{tr}(x\circ y)=2(ab+u\cdot v)`$.
Si $`a\ge|u|`$ et $`b\ge|v|`$, alors $`ab\ge|u||v|\ge-u\cdot v`$. Si $`b<|v|`$, alors $`x=(|v|,-v)\in Q`$ donne
$`2(|v|b-|v|^2)<0`$ (pour $`v=0`$, $`b<0`$, prendre $`x=1`$). (5) $`P(y)x=2y\circ(y\circ x)-y^2\circ x`$ et $`\phi`$ préserve
$`\circ`$, donc $`\phi(P(y)x)=\phi(y)\phi(x)\phi(y)`$ d'après la proposition 13.2(1) ; on applique (4) et la
proposition 13.2(3). $`\square`$

*Calcul :* le théorème 13.4 dans $`V_5`$ : la décomposition spectrale, les valeurs propres de $`\phi(x)`$, $`\phi(e^x)=e^{\phi(x)}`$ et le cône des carrés — [tutoriel, section 13.2](../../tutorial/fr/13_jordan_algebras.md#s-13-2). *Dans la bibliothèque :* `SpinFactor.spectral_decomposition`, `exp`, `sqrt`, `in_cone` — [tutoriel, section 13.2](../../tutorial/fr/13_jordan_algebras.md#s-13-2).

<a id="r-13-5"></a>

**Proposition 13.5.** (a) Pour $`m=3`$, $`\phi(a,u)=a1+u\cdot\sigma`$ est un isomorphisme de Jordan
$`V_3\cong\mathrm{Herm}_2(\mathbb C)`$. (b) Les automorphismes de Jordan de $`V_m`$ sont les applications
$`(a,u)\mapsto(a,Ou)`$, $`O\in O(m)`$. (c) Pour $`m\ge2`$, $`\phi(V_m)`$ est de dimension $`m+1`$, et la partie autoadjointe de
la $`*`$-algèbre qu'il engendre est de dimension $`D^2=4^{\lfloor m/2\rfloor}`$ ; elles ne sont égales que pour $`m=3`$.

*Démonstration.* (a) Pour $`m=3`$, $`(s_1,s_2,s_3)=(X,Y,Z)`$, et $`\phi`$ est injective entre espaces de dimension $`4`$. (b)
Les applications $`(a,u)\mapsto(a,Ou)`$ préservent $`u\cdot v`$, donc $`\circ`$. Réciproquement, soit $`\alpha`$ un
automorphisme de Jordan. Alors $`\alpha(1)`$ est une unité de $`V_m`$, donc $`\alpha(1)=1`$. D'après la formule pour $`x^2`$,
$`x^2=1`$ si et seulement si $`x=\pm1`$ ou $`x=(0,w)`$ avec $`|w|=1`$. Pour $`|w|=1`$, $`\alpha(0,w)^2=1`$ et
$`\alpha(0,w)\ne\pm1=\alpha(\pm1)`$, donc $`\alpha(0,w)=(0,w')`$ avec $`|w'|=1`$. En écrivant $`\alpha(0,u)=(\ell(u),Ou)`$ avec
$`\ell,O`$ linéaires, cela signifie $`\ell=0`$ et $`|Ou|=|u|`$, c'est-à-dire $`O\in O(m)`$. (c) Lemme 13.3 ;
$`4^{\lfloor m/2\rfloor}>m+1`$ pour $`m=2`$ et $`m\ge4`$. $`\square`$

Pour $`m=1`$, $`s_1=Z`$ engendre les matrices diagonales, de dimension $`2=m+1`$, et $`V_1\cong\mathbb R\oplus\mathbb R`$
n'est pas simple.

*Calcul :* (b) pour un $`O\in O(5)`$ aléatoire, et (c) pour $`m\le6`$ avec `jordan.compare_with_star_algebra` — [tutoriel, section 13.2](../../tutorial/fr/13_jordan_algebras.md#s-13-2).

**États de $`V_m`$.** Un **état** de $`V_m`$ est une forme linéaire $`\omega`$ telle que $`\omega(Q)\ge0`$ et $`\omega(1)=1`$.
Les états sont exactement les applications $`x\mapsto\mathrm{tr}(\rho\circ x)`$ avec $`\rho=\frac12(1,r)`$,
$`|r|\le1`$ ; l'espace des états est donc une boule de dimension $`m`$. *Démonstration.*
$`\langle x,y\rangle=\mathrm{tr}(x\circ y)=2(ab+u\cdot v)`$ est un produit scalaire sur $`V_m`$, donc toute $`\omega`$
linéaire s'écrit $`\mathrm{tr}(\rho\circ\cdot)`$ pour exactement un $`\rho=(b,v)`$ (représentation de Riesz en
dimension finie : dans une base orthonormée, $`\rho=\sum_k\omega(e_k)e_k`$). Par l'autodualité du théorème 13.4(4),
$`\omega(Q)\ge0`$ si et seulement si $`\rho\in Q`$, c'est-à-dire $`b\ge|v|`$ ; et $`\omega(1)=\mathrm{tr}\,\rho=2b=1`$.
Donc $`\rho=\frac12(1,r)`$ avec $`r=2v`$, $`|r|\le1`$ ; réciproquement, tout tel $`\rho`$ donne un état. $`\square`$

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Pour $`m=3`$, l'espace des états
> de $`V_m`$ est la boule de Bloch d'un qubit (proposition 13.5(a)). Le cône de Lorentz est
> le cône de lumière futur de la relativité restreinte en dimension $`1+m`$, et la proposition 13.5(b) dit que les
> automorphismes de Jordan sont les applications orthogonales de l'espace qui fixent l'axe du temps.

Les cônes $`Q`$ sont les ensembles de contraintes de la **programmation conique du second ordre** (SOCP) : optimisation
linéaire sur des produits de cônes de Lorentz. Les grandeurs du théorème 13.4 (valeurs propres, $`\det`$ et $`P(y)`$) sont
celles dans lesquelles on formule habituellement les méthodes de points intérieurs pour la SOCP ; ces algorithmes ne
sont pas traités dans ce cours (pour aller plus loin : F. Alizadeh, D. Goldfarb, *Second-order cone programming*,
Math. Program. 95 (2003) 3–51).

<a id="s-13-3"></a>

## 13.3 Sous-algèbres de Jordan qui ne sont pas des *-sous-algèbres

Une sous-algèbre de Jordan $`J\subseteq\mathrm{Herm}_n(\mathbb C)`$ engendre une $`*`$-algèbre complexe $`\mathcal A`$
avec $`J\subseteq\mathcal A_{sa}`$. La proposition 13.5(c) a donné des exemples avec $`J\ne\mathcal A_{sa}`$ ; en voici
deux autres.

<a id="r-13-6"></a>

**Proposition 13.6.** (a) $`\mathrm{Herm}_n(\mathbb R)`$ est une sous-algèbre de Jordan de
$`\mathrm{Herm}_n(\mathbb C)`$ de dimension $`n(n+1)/2`$ ; la $`*`$-algèbre qu'elle engendre est $`M_n(\mathbb C)`$, dont
la partie autoadjointe est de dimension $`n^2`$. (b) Pour $`x,y\in\mathrm{Herm}_2(\mathbb C)`$, il existe une
unitaire $`U`$ telle que $`UxU^*`$ et $`UyU^*`$ soient symétriques réelles. L'algèbre de Jordan engendrée par $`1,x,y`$ est
$`\mathrm{span}\{1,x,y\}`$ ; si $`x,y,1`$ sont linéairement indépendants, c'est $`U^*\mathrm{Herm}_2(\mathbb R)U\cong V_2`$,
de dimension $`3`$, tandis que la $`*`$-algèbre engendrée par $`1,x,y`$ est $`M_2(\mathbb C)`$, de dimension $`4`$.

*Démonstration.* (a) Lemme 13.1. Pour $`i\ne j`$, $`E_{ij}=E_{ii}(E_{ij}+E_{ji})`$ est un produit de matrices symétriques
réelles, donc l'algèbre engendrée contient toutes les unités matricielles. (b) Écrivons $`x=a1+u\cdot\sigma`$,
$`y=b1+v\cdot\sigma`$ et choisissons un vecteur unitaire $`n\perp u,v`$. Les matrices $`n\cdot\sigma`$ et $`Y`$ sont
hermitiennes, de trace nulle et de carré $`1`$, donc toutes deux ont les valeurs propres $`1,-1`$, et
$`U(n\cdot\sigma)U^*=Y`$ pour une certaine unitaire $`U`$ (théorème 1.3). La conjugaison par $`U`$ préserve
$`\mathrm{Tr}`$ et $`\circ`$, donc elle envoie $`u\cdot\sigma`$ sur $`(Ru)\cdot\sigma`$ avec $`R`$ linéaire et
$`Ru\cdot Rv=\frac12\mathrm{Tr}\big((Ru\cdot\sigma)(Rv\cdot\sigma)\big)=u\cdot v`$ : $`R\in O(3)`$, et $`Rn=e_2`$.
Donc $`Ru,Rv\perp e_2`$, et $`UxU^*=a1+(Ru)_1X+(Ru)_3Z`$ est réelle, puisque $`X,Z`$ le sont ; de même pour $`UyU^*`$.
Ensuite, $`x\circ y=(ab+u\cdot v)1+av\cdot\sigma+bu\cdot\sigma\in\mathrm{span}\{1,x,y\}`$, et de même pour
$`x^2,y^2`$, donc l'espace engendré est stable. S'il est de dimension $`3`$, son conjugué par $`U`$ est
$`\mathrm{Herm}_2(\mathbb R)=\mathrm{span}\{1,X,Z\}`$ tout entier, qui est $`\cong V_2`$ via
$`(a,u)\mapsto a1+u_1X+u_2Z`$ (la démonstration du lemme 13.3 n'utilise que le fait que les $`s_i`$ sont des symétries
hermitiennes anticommutantes, et $`X,Z`$ en sont), et engendre $`M_2(\mathbb C)`$ d'après (a). $`\square`$

*Calcul :* la proposition 13.6(a) pour $`n\le5`$, et (b) pour des matrices hermitiennes $`2\times2`$ aléatoires, avec `jordan.compare_with_star_algebra` — [tutoriel, section 13.3](../../tutorial/fr/13_jordan_algebras.md#s-13-3).

Trois matrices hermitiennes $`2\times2`$ génériques engendrent $`\mathrm{Herm}_2(\mathbb C)`$ tout entier. Dans un
bloc $`M_2`$ de $`\mathcal M`$, deux observables ne révèlent donc jamais la structure complexe ; une troisième le peut.

<a id="s-13-4"></a>

## 13.4 La classification

Une algèbre de Jordan de dimension finie est **formellement réelle** si $`\sum_ix_i^2=0`$ entraîne que tous les $`x_i=0`$,
et **euclidienne** si elle possède une unité $`1`$ et un produit scalaire tel que
$`\langle x\circ y,z\rangle=\langle y,x\circ z\rangle`$. Les deux notions coïncident (théorème E.5(4) ; le sens
« euclidienne $`\Rightarrow`$ formellement réelle » tient en une ligne :
$`0=\sum_i\langle x_i^2,1\rangle=\sum_i\langle x_i,x_i\rangle`$). $`M_{sa}`$ muni de
$`\langle x,y\rangle=\mathrm{Tr}(xy)`$ est euclidienne : $`\mathrm{Tr}(x^2)=\sum_{ij}|x_{ij}|^2`$, et par
cyclicité $`\mathrm{Tr}((x\circ y)z)=\frac12\mathrm{Tr}(xyz+yxz)=\mathrm{Tr}(y(x\circ z))`$. Il en va de
même de $`V_m`$ muni de $`\langle x,y\rangle=\mathrm{tr}(x\circ y)=2(ab+u\cdot v)`$ : pour $`z=(c,w)`$,
$`\mathrm{tr}((x\circ y)\circ z)=2(abc+c\,u\cdot v+a\,v\cdot w+b\,u\cdot w)`$ est symétrique en $`x,y,z`$.

<a id="r-13-7"></a>

**Théorème 13.7 (Jordan–von Neumann–Wigner).** Toute algèbre de Jordan formellement réelle de dimension finie est
somme directe d'idéaux simples, et toute algèbre simple est isomorphe à exactement une des algèbres suivantes :
$`\mathbb R`$ ; $`\mathrm{Herm}_r(\mathbb K)`$ pour $`\mathbb K=\mathbb R,\mathbb C,\mathbb H`$ et $`r\ge3`$ ; les
facteurs de spin $`V_m`$, $`m\ge2`$ ; l'**algèbre d'Albert** $`\mathrm{Herm}_3(\mathbb O)`$ des matrices hermitiennes
$`3\times3`$ à coefficients octonions (dimension $`27`$). De plus, $`\mathrm{Herm}_1(\mathbb K)\cong\mathbb R`$ et
$`\mathrm{Herm}_2(\mathbb K)\cong V_{\dim\mathbb K+1}`$ pour $`\mathbb K=\mathbb R,\mathbb C,\mathbb H,\mathbb O`$.

*Démonstration.* L'annexe E le démontre entièrement : la décomposition est le théorème E.9, la liste et « exactement
une » sont le théorème E.27(1),(2), et les isomorphismes en rang un et deux sont le théorème E.24(3). $`\square`$

(P. Jordan, J. von Neumann, E. Wigner, *On an algebraic generalization of the quantum mechanical formalism*, Ann. of
Math. 35 (1934) 29–64 ; la démonstration de l'annexe E suit [FK94, chapitres II–V].)

Toutes sauf l'algèbre d'Albert sont spéciales : $`\mathbb R=\mathrm{Herm}_1(\mathbb R)`$ et
$`\mathrm{Herm}_r(\mathbb K)`$, $`\mathbb K=\mathbb R,\mathbb C`$, d'après le lemme 13.1,
$`\mathrm{Herm}_r(\mathbb H)`$ parce qu'écrire les quaternions comme des matrices complexes $`2\times2`$ la plonge
dans $`\mathrm{Herm}_{2r}(\mathbb C)`$ (théorème E.24(1)), et $`V_m`$ d'après le lemme 13.3. L'algèbre d'Albert est
une algèbre de Jordan (théorème E.24(2)) et elle est **exceptionnelle** : elle n'est isomorphe à aucune sous-algèbre de
Jordan d'une algèbre associative réelle ou complexe (théorème E.26, un argument du type de celui qu'utilisa A. A.
Albert, *On a certain algebra of quantum mechanics*, Ann. of Math. 35 (1934) 65–73 : un plongement de Jordan dans une
algèbre associative rendrait associatif le produit des octonions, ce qu'il n'est pas, lemme E.15(2)). Pour $`r\ge4`$,
$`\mathrm{Herm}_r(\mathbb O)`$ ne vérifie pas l'identité de Jordan (corollaire E.25), et une algèbre simple
possédant un repère de Jordan à $`r\ge4`$ éléments est $`\mathrm{Herm}_r(D)`$ pour une algèbre de composition
associative $`D\ne\mathbb O`$ (théorème E.22), ce qui explique que les octonions n'apparaissent que pour $`r\le3`$. Ainsi
$`\mathrm{Herm}_3(\mathbb O)`$ ne peut être représentée fidèlement par des matrices réelles ou complexes (ni dans
aucune algèbre associative) munies du produit $`\frac12(xy+yx)`$ ; un calcul dans cette algèbre doit utiliser le produit
des octonions (ou ses constantes de structure) au lieu d'un produit matriciel.

Soit $`\Omega`$ un cône convexe ouvert non vide d'un espace euclidien $`V`$ (un espace préhilbertien réel de dimension
finie $`\ne\{0\}`$). Son **cône dual** est
$`\Omega^*=\{y:\langle x,y\rangle>0\text{ pour tout }x\in\bar\Omega\setminus\{0\}\}`$ ; $`\Omega`$ est **autodual** si
$`\Omega^*=\Omega`$, **homogène** si les applications linéaires $`g`$ telles que $`g\Omega=\Omega`$ agissent transitivement
sur $`\Omega`$, et **symétrique** s'il est les deux à la fois.

<a id="r-13-8"></a>

**Théorème 13.8 (Koecher–Vinberg).** Un cône convexe ouvert non vide $`\Omega`$ d'un espace euclidien $`V`$ est symétrique
si et seulement s'il existe un produit $`\circ`$ sur $`V`$ tel que $`(V,\circ,\langle\cdot,\cdot\rangle)`$ soit une algèbre de
Jordan euclidienne et que $`\Omega`$ soit l'intérieur de son cône des carrés.

*Démonstration.* Théorème F.12 (« si » : théorèmes F.5(2) et F.6(2)). $`\square`$

(M. Koecher, 1957 ; E. B. Vinberg, 1960 ; la démonstration de l'annexe F suit [FK94, chapitres I–III].) Le théorème
13.4(4) démontre la forme fermée de l'autodualité pour le cône de Lorentz ($`y\in Q`$ si et seulement si
$`\mathrm{tr}(x\circ y)\ge0`$ pour tout $`x\in Q`$), et l'exercice 1 son homogénéité. Le cône des carrés d'une somme
directe d'algèbres de Jordan euclidiennes (munie de la somme orthogonale des produits scalaires) est le produit de leurs
cônes des carrés, puisque $`(x_1,\dots,x_k)^2=(x_1^2,\dots,x_k^2)`$. Ainsi, d'après le sens « si », les intérieurs des
cônes de la programmation semi-définie (produits de $`\mathrm{Herm}_n(\mathbb R)^+`$ ou de
$`\mathrm{Herm}_n(\mathbb C)^+`$, qui sont les cônes des carrés d'après le théorème 1.6 ; pour $`x\ge0`$ réel, la
racine carrée est réelle, puisque sa conjuguée complexe est aussi une racine carrée positive de $`x`$), de la SOCP
(produits de cônes de Lorentz) et de la programmation linéaire ($`\mathbb R_+^n`$, le cône des carrés de
$`\mathbb R^n=\mathbb R\oplus\dots\oplus\mathbb R`$) sont des cônes symétriques.

[FK94]: J. Faraut, A. Korányi, *Analysis on Symmetric Cones*, Oxford University Press, 1994.

<a id="s-13-5"></a>

## 13.5 Applications de Jordan et *-applications

Un **$`\ast`$-homomorphisme de Jordan** d'algèbres de von Neumann est une application $`\mathbb C`$-linéaire $`\Phi`$ telle que
$`\Phi(x^*)=\Phi(x)^*`$ et $`\Phi(x\circ y)=\Phi(x)\circ\Phi(y)`$ ; sur les parties autoadjointes, c'est un homomorphisme
de Jordan, et tout homomorphisme de Jordan entre parties autoadjointes s'étend en un tel $`*`$-homomorphisme par
$`\Phi(x+iy)=\Phi(x)+i\Phi(y)`$ ($`x,y`$ autoadjoints). (Ce $`\Phi`$ est $`\mathbb C`$-linéaire et vérifie
$`\Phi(x-iy)=\Phi(x+iy)^*`$ ; et comme $`\circ`$ est $`\mathbb C`$-bilinéaire et envoie les couples d'éléments autoadjoints
sur des éléments autoadjoints, $`(x+iy)\circ(x'+iy')=(x\circ x'-y\circ y')+i(x\circ y'+y\circ x')`$ est envoyé sur
$`\Phi(x+iy)\circ\Phi(x'+iy')`$.)

<a id="r-13-9"></a>

**Proposition 13.9 (la transposée).** Soit $`T(A)=\bigoplus_cA_c^T`$ sur $`\mathcal M`$. (a) $`T`$ est un
$`*`$-automorphisme de Jordan unifère avec $`T(AB)=T(B)T(A)`$ ; il préserve les spectres, donc $`T`$ et $`T^{-1}=T`$ sont
positifs. (b) Sur $`M_k`$, $`k\ge2`$, la matrice de Choi $`\sum_{ij}E_{ij}\otimes T(E_{ij})`$ est l'opérateur d'échange
$`F(e_a\otimes e_b)=e_b\otimes e_a`$, de valeurs propres $`+1`$ (multiplicité $`k(k+1)/2`$) et $`-1`$ (multiplicité
$`k(k-1)/2`$). Donc $`T`$ n'est pas complètement positif.

*Démonstration.* (a) $`(A^T)^*=(A^*)^T`$, $`(AB)^T=B^TA^T`$, d'où $`(A\circ B)^T=B^T\circ A^T=A^T\circ B^T`$ ;
$`\det(A^T-\lambda)=\det(A-\lambda)`$. (b) $`\sum_{ij}E_{ij}e_a\otimes E_{ji}e_b=e_b\otimes e_a`$. $`F^2=1`$, et $`F`$ vaut
$`+1`$ sur les tenseurs symétriques et $`-1`$ sur les tenseurs antisymétriques, qui engendrent l'espace. D'après le
théorème 6.2 (et la proposition 6.1), $`T`$ n'est pas CP. $`\square`$

*Calcul :* la proposition 13.9 avec la matrice de Choi de la transposée sur $`M_3`$ — [tutoriel, section 13.5](../../tutorial/fr/13_jordan_algebras.md#s-13-5). *Dans la bibliothèque :* la transposée $`T`$ est `jordan.transpose_map` — [tutoriel, section 13.5](../../tutorial/fr/13_jordan_algebras.md#s-13-5).

Un automorphisme de Jordan n'est donc pas nécessairement un $`*`$-automorphisme. À ce phénomène près, la positivité
détermine la structure de Jordan :

<a id="r-13-10"></a>

**Théorème 13.10 (Kadison, dimension finie).** Soient $`\mathcal M=\bigoplus_cM_{k_c}(\mathbb C)`$ et
$`\mathcal N=\bigoplus_dM_{k'_d}(\mathbb C)`$, et soit $`\Phi:\mathcal M\to\mathcal N`$ une bijection linéaire, avec
$`\Phi(1)=1`$, telle que $`\Phi`$ et $`\Phi^{-1}`$ soient positives. Alors $`\Phi`$ est un $`*`$-isomorphisme de Jordan.
Réciproquement, tout $`*`$-isomorphisme de Jordan $`\mathcal M\to\mathcal N`$ est unifère, et lui et son inverse sont
positifs.

*Démonstration.* Théorème F.15. (L'idée : $`\Phi`$ envoie l'ensemble convexe $`\{0\le a\le1\}`$ sur l'ensemble
correspondant de $`\mathcal N`$, donc les points extrémaux, qui sont les projections, sur des projections ; elle préserve
$`p+q\le1`$, ce qui pour des projections signifie $`pq=0`$ ; elle envoie donc les décompositions spectrales sur des
décompositions spectrales, et $`\Phi(x^2)=\Phi(x)^2`$.) $`\square`$

R. V. Kadison a démontré le théorème pour toutes les C\*-algèbres unifères (*A generalized Schwarz inequality and
algebraic invariants for operator algebras*, Ann. of Math. 56 (1952) 494–503) ; seul le cas de dimension finie est
démontré et utilisé dans ce cours. Pour les facteurs $`M_n`$, les automorphismes de Jordan sont connus exactement : tout
$`*`$-automorphisme de Jordan de $`M_n(\mathbb C)`$ est $`A\mapsto UAU^*`$ ou $`A\mapsto UA^TU^*`$ avec $`U`$ unitaire
(théorème F.17, pour tout $`n`$ ; le résultat est dû à N. Jacobson, C. E. Rickart, Trans. Amer. Math. Soc. 69 (1950)
479–502, et à I. N. Herstein, Trans. Amer. Math. Soc. 81 (1956) 331–341). Les premières applications sont
multiplicatives et les secondes antimultiplicatives, donc un $`*`$-automorphisme de Jordan de $`M_n`$ soit préserve tous
les commutateurs, soit les renverse tous : $`\Phi([x,y])=\pm[\Phi(x),\Phi(y)]`$. Sur $`\bigoplus_cM_{k_c}`$, un
$`*`$-isomorphisme de Jordan permute des blocs de même taille et agit sur chaque bloc de l'une des deux manières
(théorème F.18) ; il est unifère et positif dans les deux sens (théorème 13.10), et d'après le corollaire F.19 il est
complètement positif si et seulement s'il ne transpose aucun bloc de taille $`\ge2`$, si et seulement si c'est un
$`*`$-isomorphisme. Voici une démonstration directe du cas $`n=2`$ dans le langage des facteurs de spin.

<a id="r-13-11"></a>

**Proposition 13.11.** Tout automorphisme de Jordan de $`\mathrm{Herm}_2(\mathbb C)`$ est $`A\mapsto UAU^*`$ ou
$`A\mapsto UA^TU^*`$ avec $`U\in U(2)`$.

*Démonstration.* D'après la proposition 13.5(a,b), un automorphisme de Jordan agit par
$`a1+u\cdot\sigma\mapsto a1+(Ou)\cdot\sigma`$ avec $`O\in O(3)`$. Notons $`R_U\in O(3)`$ la matrice de $`A\mapsto UAU^*`$ en
ce sens ; $`R_{UV}=R_UR_V`$. (i) $`\det R_U=1`$ : $`U=W\mathrm{diag}(e^{i\theta_1},e^{i\theta_2})W^*`$ est relié à $`1`$
par le chemin $`U_t=W\mathrm{diag}(e^{it\theta_1},e^{it\theta_2})W^*`$, et $`\det R_{U_t}\in\{\pm1\}`$ dépend
continûment de $`t`$. (ii) Tout $`R\in SO(3)`$ est un $`R_U`$ : comme dans la proposition 13.6, il existe $`U_1`$ tel que
$`U_1ZU_1^*=(Re_3)\cdot\sigma`$, c'est-à-dire $`R_{U_1}e_3=Re_3`$. Alors $`R_{U_1}^{-1}R\in SO(3)`$ fixe $`e_3`$, c'est donc
une rotation d'un certain angle $`\beta`$ autour de $`e_3`$. Un calcul direct donne, pour
$`U_2=\mathrm{diag}(e^{-i\beta/2},e^{i\beta/2})`$, $`U_2XU_2^*=\cos\beta\,X+\sin\beta\,Y`$ et
$`U_2YU_2^*=-\sin\beta\,X+\cos\beta\,Y`$, donc $`R=R_{U_1U_2}`$. (iii) La transposée fixe $`X,Z`$ et envoie $`Y\mapsto-Y`$ :
$`O_T=\mathrm{diag}(1,-1,1)`$, $`\det O_T=-1`$. Donc $`O=R_U`$ si $`\det O=1`$, et $`OO_T=R_U`$, c'est-à-dire que
l'automorphisme est $`A\mapsto UA^TU^*`$, si $`\det O=-1`$. $`\square`$

*Calcul :* $`O_T=\mathrm{diag}(1,-1,1)`$ — [tutoriel, section 13.5](../../tutorial/fr/13_jordan_algebras.md#s-13-5).

**Systèmes composés.** La structure de Jordan d'un facteur ne permet pas de distinguer $`T`$ de l'identité. Pour deux
systèmes, cette ambiguïté devient visible :

<a id="r-13-12"></a>

**Proposition 13.12.** Soient $`x,x'\in\mathrm{Herm}_n(\mathbb C)`$, $`y,y'\in\mathrm{Herm}_m(\mathbb C)`$.
Dans $`\mathrm{Herm}_{nm}(\mathbb C)`$, agissant sur $`\mathbb C^n\otimes\mathbb C^m`$,
```math
(x\otimes y)\circ(x'\otimes y')=(x\circ x')\otimes(y\circ y')+\tfrac14[x,x']\otimes[y,y'].
```
Soient $`\theta=\mathrm{id}\otimes T`$ (la **transposée partielle**) et $`A\circ_\theta B=\theta(\theta(A)\circ\theta(B))`$.
Alors $`\circ_\theta`$ est un produit de Jordan sur $`\mathrm{Herm}_{nm}(\mathbb C)`$, isomorphe à $`\circ`$ via
$`\theta`$, qui coïncide avec $`\circ`$ sur toutes les paires d'observables locales $`x\otimes1`$, $`1\otimes y`$, et
$`(x\otimes y)\circ_\theta(x'\otimes y')=(x\circ x')\otimes(y\circ y')-\frac14[x,x']\otimes[y,y']`$. Son cône des carrés
est $`\theta(\mathrm{Herm}_{nm}^+)`$. Pour $`n=m=2`$ et $`\Omega=e_1\otimes e_1+e_2\otimes e_2`$, l'état
$`\omega=\frac12|\Omega\rangle\langle\Omega|\ge0`$ n'appartient pas à ce cône : $`\theta(\omega)=\frac12F`$ a la valeur
propre $`-\frac12`$.

*Démonstration.* Avec $`j=x\circ x'`$, $`c=\frac12[x,x']`$, $`k=y\circ y'`$, $`d=\frac12[y,y']`$ : $`xx'=j+c`$, $`x'x=j-c`$, et de
même pour $`y`$, donc
$`\frac12(xx'\otimes yy'+x'x\otimes y'y)=\frac12\big((j+c)\otimes(k+d)+(j-c)\otimes(k-d)\big)=j\otimes k+c\otimes d`$.
$`\theta`$ est une involution $`\mathbb R`$-linéaire de $`\mathrm{Herm}_{nm}`$ : elle envoie $`x\otimes y`$ sur
$`x\otimes y^T`$, et les $`x\otimes y`$ engendrent $`\mathrm{Herm}_{nm}`$ sur $`\mathbb R`$ (les produits
$`h_\alpha\otimes g_\beta`$ de bases réelles de $`\mathrm{Herm}_n,\mathrm{Herm}_m`$ forment une base complexe
de $`M_{nm}`$, et une matrice hermitienne y a des coefficients réels). Donc $`\circ_\theta`$ est un produit de Jordan et
$`\theta`$ un isomorphisme. Comme $`(y^T\circ y'^T)^T=y\circ y'`$ et $`[y^T,y'^T]^T=-[y,y']`$, appliquer $`\theta`$ à la
formule pour $`(x\otimes y^T)\circ(x'\otimes y'^T)`$ donne la formule pour $`\circ_\theta`$. Sur les observables locales,
le terme de commutateur s'annule ($`[1,1]=0`$, $`[x,1]=0`$). Les carrés sont $`A\circ_\theta A=\theta(\theta(A)^2)`$, et
$`\omega`$ appartient à $`\theta(\mathrm{Herm}_{nm}^+)`$ si et seulement si $`\theta(\omega)\ge0`$, puisque
$`\theta^2=\mathrm{id}`$. Enfin, $`\theta(\omega)=\frac12\sum_{ij}E_{ij}\otimes E_{ji}=\frac12F`$ (proposition 13.9(b) ;
c'est le calcul de la démonstration de la proposition 6.1). $`\square`$

*Calcul :* la formule du produit, $`\circ_\theta`$, et la valeur propre $`-\frac12`$ de $`\theta(\omega)`$ — [tutoriel, section 13.5](../../tutorial/fr/13_jordan_algebras.md#s-13-5).

Les deux produits sont des algèbres de Jordan aussi bonnes l'une que l'autre, ils coïncident sur tout ce qui peut être
mesuré sur un sous-système à la fois, et ils diffèrent sur la question de savoir quels opérateurs joints sont
positifs, c'est-à-dire sur ce que sont les états du système composé. Ce qui les distingue est le signe du terme de
commutateur, et ce signe est fixé par la $`*`$-algèbre $`M_n\otimes M_m`$ et non par les algèbres de Jordan des facteurs.
($`\circ_\theta`$ est le produit de Jordan du produit associatif $`A\star B=\theta(\theta(A)\theta(B))`$, pour lequel
$`(x\otimes y)\star(x'\otimes y')=xx'\otimes y'y`$ : l'algèbre $`M_n\otimes M_m^{\rm op}`$.) Le cône
$`\theta(\mathrm{Herm}_{nm}^+)`$ est ce qui donne le **critère PPT** (A. Peres, Phys. Rev. Lett. 77 (1996) 1413).
Un état $`\rho`$ de $`M_n\otimes M_m`$ est dit **séparable** si $`\rho=\sum_kp_ka_k\otimes b_k`$ (une somme finie) avec
$`p_k\ge0`$, $`a_k,b_k\ge0`$, et **intriqué** sinon. Un état séparable appartient aux deux cônes, puisque
$`\theta(\rho)=\sum_kp_ka_k\otimes b_k^T`$ avec $`b_k^T\ge0`$ (proposition 13.9(a)) et
$`a\otimes b=(a^{1/2}\otimes b^{1/2})^2\ge0`$. Donc $`\omega`$, qui n'appartient pas à $`\theta(\mathrm{Herm}_4^+)`$,
est intriqué. (Pour aller plus loin sur les composés d'algèbres de Jordan, non utilisé ici : H. Hanche-Olsen, *On the
structure and tensor products of JC-algebras*, Canad. J. Math. 35 (1983) 1059–1074 ; H. Barnum, M. Graydon,
A. Wilce, *Composites and categories of Euclidean Jordan algebras*, Quantum 4 (2020) 359.)

**Dynamique.** On pourrait espérer que l'évolution temporelle impose le produit associatif. Il n'en est rien :

<a id="r-13-13"></a>

**Proposition 13.13.** Soient $`H\in\mathrm{Herm}_n(\mathbb C)`$ et $`\delta(z)=i[H,z]`$. (a) $`\delta`$ envoie
$`\mathrm{Herm}_n`$ dans lui-même et est une dérivation de $`\circ`$ ; $`e^{t\delta}(z)=e^{itH}ze^{-itH}`$ est un
automorphisme de Jordan. (b) $`[L(x),L(y)]z=\frac14[[x,y],z]`$ pour $`x,y,z\in\mathrm{Herm}_n`$. (c)
$`\delta=4\sum_k[L(x_k),L(y_k)]`$ pour certains $`x_k,y_k\in\mathrm{Herm}_n`$ : la dynamique hamiltonienne est
engendrée par le seul produit de Jordan.

*Démonstration.* (a) $`(i[H,z])^*=-i[z,H]=i[H,z]`$ ; $`[H,\cdot]`$ est une dérivation du produit associatif, donc de
$`\circ`$. $`z(t)=e^{itH}ze^{-itH}`$ vérifie $`z'=i[H,z]`$ et $`z(0)=z`$ (proposition 8.2(3) avec $`\mathcal M=M_n`$), donc
$`z(t)=e^{t\delta}(z)`$ par l'unicité du lemme 8.1 ; la conjugaison par l'unitaire $`e^{itH}`$ préserve les produits et les
adjoints, donc $`\circ`$ et $`\mathrm{Herm}_n`$.
(b) $`4x\circ(y\circ z)=xyz+xzy+yzx+zyx`$ ; on soustrait la même expression avec $`x,y`$ échangés. (c) D'après
(b), il suffit d'écrire $`iH-\frac{\mathrm{Tr}\,H}n i1`$ (qui a les mêmes commutateurs que $`iH`$) comme somme de
$`[x_k,y_k]`$. Pour $`i\ne j`$ et $`a=E_{ij}+E_{ji}`$, $`b=i(E_{ij}-E_{ji})`$ : $`[a,E_{ii}]=E_{ji}-E_{ij}`$,
$`[b,E_{ii}]=-i(E_{ij}+E_{ji})`$ et $`[a,b]=2i(E_{jj}-E_{ii})`$ ; ces matrices engendrent sur $`\mathbb R`$ les matrices
antihermitiennes de trace nulle. $`\square`$

*Calcul :* (a) et (b) pour des matrices $`3\times3`$ aléatoires — [tutoriel, section 13.5](../../tutorial/fr/13_jordan_algebras.md#s-13-5).

Ce que l'algèbre de Jordan ne contient pas, c'est la **correspondance** $`H\mapsto i[H,\cdot]`$ entre une observable et
la dynamique qu'elle engendre. Dans $`\mathrm{Herm}_n(\mathbb R)`$, par exemple, les dérivations
$`[L(x),L(y)]=\frac14[[x,y],\cdot]`$ proviennent de matrices $`[x,y]`$ réelles antisymétriques, qui ne sont pas $`i`$ fois
une matrice symétrique réelle (exercice 4). La correspondance n'est pas déterminée par $`\circ`$ : sur
$`\mathrm{Herm}_n(\mathbb C)`$, $`H\mapsto i[H,\cdot]`$ et $`H\mapsto-i[H,\cdot]`$ envoient toutes deux chaque
observable sur une dérivation de $`\circ`$ (proposition 13.13(a)), et la transposée, un automorphisme de Jordan, les
échange, puisque $`T(i[H,T(z)])=-i[H^T,z]`$. La seconde est $`i`$ fois le commutateur de l'algèbre opposée
($`a\cdot_{\rm op}b=ba`$) : $`-i[H,z]=i(H\cdot_{\rm op}z-z\cdot_{\rm op}H)`$ — les deux signes de la proposition 13.12.
(Pour aller plus loin sur ces « correspondances dynamiques », non utilisé ici : E. M. Alfsen, F. W. Shultz, *Geometry
of State Spaces of Operator Algebras*, Birkhäuser, 2003.)

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Pour une particule sans spin,
> la transposée (la conjugaison complexe dans une base réelle) est le **renversement du temps**. C'est une symétrie des observables d'un système unique, mais elle est
> antiunitaire et, appliquée à un seul sous-système, ce n'est pas une opération physique : la proposition 13.12 montre
> qu'elle transforme l'état de Bell en un opérateur ayant une valeur propre négative. La mécanique quantique utilise des
> $`*`$-algèbres complexes parce que les systèmes composés ont besoin du commutateur avec son signe, et pas seulement du
> produit symétrisé.

<a id="s-13-6"></a>

## 13.6 Retour aux inégalités de traces

La proposition 11.9 (Exp. 3) énonce $`\mathrm{Tr}|YXY|\le\mathrm{Tr}(Y|X|Y)`$ pour $`X=X^*`$ et $`Y\ge0`$ (ici
$`X,Y`$ sont des opérateurs de $`\mathcal M`$, et non des matrices de Pauli). Sa démonstration n'utilise jamais $`Y\ge0`$ :
$`YX_\pm Y=T^*T`$ avec $`T=X_\pm^{1/2}Y`$ ne demande que $`Y=Y^*`$, et le lemme 11.8 est appliqué à $`YX_+Y`$ et $`YX_-Y`$. Elle
est donc vraie pour tout $`Y=Y^*`$ et s'écrit, dans le langage de Jordan,
```math
\mathrm{tr}|P(y)x|\le\mathrm{tr}\,P(y)|x|,
```
ce qui a un sens dans toute algèbre de Jordan euclidienne une fois une trace fixée. Nous la démontrons dans $`V_m`$ par un argument qui
n'utilise que trois propriétés.

<a id="r-13-14"></a>

**Proposition 13.14 (Exp. 3 dans un facteur de spin).** Pour $`x,y\in V_m`$,
$`\mathrm{tr}|P(y)x|\le\mathrm{tr}\,P(y)|x|`$, avec égalité si $`x\in Q`$ ou $`-x\in Q`$.

*Démonstration.* Soient $`x_\pm=f_\pm(x)`$ avec $`f_+(t)=\max(t,0)`$, $`f_-(t)=\max(-t,0)`$ (théorème 13.4(3)) :
$`x_\pm\in Q`$, $`x=x_+-x_-`$ et $`|x|=x_++x_-`$. Posons $`A=P(y)x_+`$ et $`B=P(y)x_-`$ ; ils appartiennent à $`Q`$ d'après le
théorème 13.4(5), et $`P(y)x=A-B`$. Soient $`P(y)x=\mu_1c_1+\mu_2c_2`$ une décomposition spectrale (théorème 13.4(2) ;
$`c_1=c_w`$, $`c_2=c_{-w}`$, donc $`c_1+c_2=1`$ et $`c_1\circ c_2=0`$), $`c`$ la somme des $`c_i`$ tels que $`\mu_i>0`$, et
$`s=2c-1=\sum_i\mathrm{sgn}'(\mu_i)c_i`$ avec $`\mathrm{sgn}'(t)=1`$ pour $`t>0`$ et $`-1`$ sinon. Alors
$`s\circ P(y)x=\sum_i|\mu_i|c_i=|P(y)x|`$ (d'après le théorème 13.4(2,3), cela ne dépend pas du choix de $`w`$ lorsque
$`\mu_1=\mu_2`$), et $`c,1-c\in Q`$ (les idempotents sont des carrés). Par la moitié facile de l'autodualité (théorème
13.4(4) : $`\mathrm{tr}(a\circ b)\ge0`$ pour $`a,b\in Q`$), appliquée à $`\mathrm{tr}(A\circ(1-c))`$ et à
$`\mathrm{tr}(B\circ c)`$,
```math
\mathrm{tr}|P(y)x|=\mathrm{tr}(A\circ c)-\mathrm{tr}(A\circ(1-c))-\mathrm{tr}(B\circ c)+\mathrm{tr}(B\circ(1-c))
\le\mathrm{tr}\,A+\mathrm{tr}\,B=\mathrm{tr}\,P(y)|x|.
```
Si $`x\in Q`$, alors $`x_-=0`$, $`B=0`$ et $`|P(y)x|=A`$, donc les deux membres valent $`\mathrm{tr}\,A`$ ; de même pour
$`-x\in Q`$. $`\square`$

*Calcul :* la proposition 13.14 pour $`m=2,4,7`$, et, via $`\phi`$, l'inégalité de la proposition 11.9 pour $`Y=Y^*`$ — [tutoriel, section 13.6](../../tutorial/fr/13_jordan_algebras.md#s-13-6).

La démonstration n'utilise que : (i) $`P(y)Q\subseteq Q`$ pour tout $`y`$ ; (ii) le théorème spectral : tout $`z`$ s'écrit
$`\sum_i\mu_ic_i`$ avec des idempotents orthogonaux $`c_i`$ de somme $`1`$, et $`f(z)=\sum_if(\mu_i)c_i`$ (de sorte que
$`x_\pm`$ et les idempotents sont des carrés) ; (iii) $`\mathrm{tr}`$ est linéaire,
$`\mathrm{tr}(a\circ1)=\mathrm{tr}\,a`$, et $`\mathrm{tr}(a\circ b)\ge0`$ pour $`a,b\in Q`$. Toutes trois sont
vraies dans toute algèbre de Jordan euclidienne $`V`$ d'unité $`e`$ pour la forme linéaire $`\tau(x)=\langle x,e\rangle`$,
pour laquelle $`\tau(a\circ b)=\langle b,a\circ e\rangle=\langle a,b\rangle`$ : (ii) est le théorème F.3, (i) est le
théorème F.6(4), et $`\langle a,b\rangle\ge0`$ pour $`a,b\in Q`$ est le théorème F.5(1). Le corollaire F.7 mène
l'argument à son terme : $`\tau(|P(y)x|)\le\tau(P(y)|x|)`$ dans toute algèbre de Jordan euclidienne, avec égalité si
$`x\in Q`$ ou $`-x\in Q`$. Sur $`M_{sa}`$ muni de $`\langle x,y\rangle=\mathrm{Tr}(xy)`$, $`\tau=\mathrm{Tr}`$ et
$`P(y)x=yxy`$ ; les idempotents de $`(M_{sa},\circ)`$ sont les projections ($`c\circ c=c^2`$), et pour des projections
$`c\circ c'=0`$ si et seulement si $`cc'=0`$ (lemme F.14(2)), donc, par l'unicité dans le théorème F.3, sa décomposition de
$`x`$ est la décomposition spectrale selon les projecteurs propres et $`|x|`$ est la valeur absolue matricielle. Cela
redémontre donc la proposition 11.9 pour tout $`Y=Y^*`$. Sur $`V_m`$ muni de $`\langle x,y\rangle=\mathrm{tr}(x\circ y)`$,
$`\tau=\mathrm{tr}`$, la décomposition du théorème F.3 est celle du théorème 13.4(2) (avec $`x=a1`$ si $`u=0`$), et
c'est la proposition 13.14. Dans une algèbre de Jordan euclidienne simple $`J`$ de rang $`r`$, la trace
$`\mathrm{tr}\,x=\frac r{\dim J}\mathrm{Tr}\,L(x)`$ du corollaire E.28 vérifie
$`\mathrm{tr}(\sum_i\mu_ic_i)=\sum_i\mu_i`$ pour tout repère de Jordan $`(c_i)`$ (théorème E.7), et
$`\langle x,y\rangle=\mathrm{tr}(x\circ y)`$ est un produit scalaire associatif ; pour celui-ci,
$`\tau=\mathrm{tr}`$, donc l'inégalité est vraie pour cette trace. Cela s'applique à
$`\mathrm{Herm}_r(\mathbb H)`$ et à l'algèbre d'Albert, qui sont des algèbres de Jordan euclidiennes simples
(théorème E.24). Pour $`r\ge3`$, elles ne sont Jordan-isomorphes à $`M_{sa}`$ pour aucun $`\mathcal M=\bigoplus_cM_{k_c}`$ :
$`M_{sa}=\bigoplus_c\mathrm{Herm}_{k_c}(\mathbb C)`$ est une somme directe d'idéaux, donc n'est simple qu'avec un
seul bloc, et $`\mathrm{Herm}_k(\mathbb C)`$ ($`\cong\mathbb R`$ pour $`k=1`$, $`\cong V_3`$ pour $`k=2`$) n'est isomorphe
ni à $`\mathrm{Herm}_r(\mathbb H)`$, $`r\ge3`$, ni à $`\mathrm{Herm}_3(\mathbb O)`$, par l'unicité dans le
théorème 13.7. Dans $`V_m`$, l'inégalité découle aussi de la proposition 11.9 (étendue ci-dessus à $`Y=Y^*`$) via $`\phi`$,
puisque $`\phi(P(y)x)=\phi(y)\phi(x)\phi(y)`$ et $`\phi(|x|)=|\phi(x)|`$.

<a id="s-13-7"></a>

## 13.7 Points de comparaison

*Dans la bibliothèque :* la comparaison avec d'autres logiciels (SageMath, Magma, solveurs SOCP) et ce qu'apporte la bibliothèque — [tutoriel, section 13.7](../../tutorial/fr/13_jordan_algebras.md#s-13-7).

<a id="x-13"></a>

## Exercices

1. (**Homogénéité du cône de Lorentz.**) Soit $`y\in V_m`$ avec $`\det y\ne0`$. Montrer que $`y^{-1}=(a,-u)/\det y`$
   vérifie $`y\circ y^{-1}=1`$, que $`P(y)P(y^{-1})=\mathrm{id}`$, que $`P(y)`$ envoie l'intérieur de $`Q`$ sur lui-même, et
   que $`P(x^{1/2})1=x`$ pour $`x`$ intérieur. En déduire que l'intérieur de $`Q`$ est homogène.
2. (**États isotropes.**) Pour $`p\in[0,1]`$, soit $`\rho_p=p\,\omega+(1-p)\frac14`$ sur $`\mathbb C^2\otimes\mathbb C^2`$
   (proposition 13.12). Calculer le spectre de $`\theta(\rho_p)`$ et montrer que $`\rho_p`$ est intriqué pour $`p>\frac13`$.
3. (**Produits de points du bord.**) Soient $`|w|=|w'|=1`$ dans $`\mathbb R^m`$. Montrer que $`(1,w)\circ(1,w')\in Q`$ si et
   seulement si $`w'=w`$ ou $`w'=-w`$ : le produit de Jordan de deux rayons du bord du cône de Lorentz n'est presque jamais
   dans le cône.
4. (**Mécanique quantique réelle.**) Montrer que, pour $`x,y\in\mathrm{Herm}_n(\mathbb R)`$, la dérivation
   $`[L(x),L(y)]`$ de $`\mathrm{Herm}_n(\mathbb R)`$ est $`z\mapsto[K,z]`$ avec $`K`$ réelle antisymétrique, et que pour
   $`n\ge2`$ une telle application n'est jamais de la forme $`z\mapsto i[H,z]`$ avec $`H\in\mathrm{Herm}_n(\mathbb R)`$,
   sauf si elle est nulle.

*Calcul :* les exercices 1–3, numériquement — [tutoriel, Exercices](../../tutorial/fr/13_jordan_algebras.md#x-13).

<details>
<summary>Solutions</summary>

1. $`(a,u)\circ(a,-u)=(a^2-|u|^2,0)`$, ce qui donne $`y\circ y^{-1}=1`$. $`\phi(y)`$ et $`\phi(y^{-1})`$ commutent (tous deux
   sont dans $`\mathrm{span}\{1,u\cdot s\}`$), donc $`\phi(y)\phi(y^{-1})=\phi(y)\circ\phi(y^{-1})=1`$ et, d'après le
   théorème 13.4(5), $`\phi(P(y)P(y^{-1})x)=\phi(x)`$ ; $`\phi`$ est injective. $`P(y)`$ envoie $`Q`$ dans $`Q`$, de même que son
   inverse $`P(y^{-1})`$, donc $`P(y)Q=Q`$ ; étant un homéomorphisme linéaire, il envoie l'intérieur sur l'intérieur.
   Enfin, $`P(y)1=2y\circ y-y^2=y^2`$, donc $`P(x^{1/2})1=x`$, et $`P(x^{1/2})`$ est inversible lorsque $`x`$ est intérieur
   ($`\det x^{1/2}=\sqrt{\lambda_+\lambda_-}>0`$ ; l'intérieur de $`Q`$ est $`\{a>|u|\}=\{\lambda_->0\}`$, puisque $`a-|u|`$
   est continue et que $`(a-\varepsilon,u)\notin Q`$ lorsque $`a=|u|`$). Deux points intérieurs quelconques sont reliés par
   $`P(x'^{1/2})P(x^{1/2})^{-1}`$.

2. $`\theta(\rho_p)=\frac p2F+\frac{1-p}4`$ d'après la proposition 13.12. $`F`$ a pour valeurs propres $`1`$ (trois fois) et
   $`-1`$ (une fois), donc le spectre est $`\frac{1+p}4`$ (trois fois) et $`\frac{1-3p}4`$, qui est négatif si et seulement
   si $`p>\frac13`$ ; $`\rho_p`$ n'est alors pas séparable d'après le critère PPT.

3. $`(1,w)\circ(1,w')=(1+t,\,w+w')`$ avec $`t=w\cdot w'\in[-1,1]`$ et $`|w+w'|=\sqrt{2+2t}`$. Il appartient à $`Q`$ si et
   seulement si $`1+t\ge\sqrt{2(1+t)}`$, c'est-à-dire si et seulement si $`t=-1`$ (le produit est $`0`$) ou $`1+t\ge2`$,
   c'est-à-dire $`t=1`$.

4. D'après la proposition 13.13(b), $`[L(x),L(y)]z=[K,z]`$ avec $`K=\frac14[x,y]`$, et $`K^T=\frac14[y^T,x^T]=-K`$ est
   réelle. Pour $`H`$ et $`z`$ réels, $`i[H,z]`$ est imaginaire pure, tandis que $`[K,z]`$ est réelle, donc elles ne
   coïncident que si toutes deux sont nulles. Ainsi, en mécanique quantique réelle, les générateurs de la dynamique,
   $`K`$, ne sont pas des observables.

</details>

**Pour aller plus loin.** Ce chapitre est le dernier du cours ; les annexes E et F complètent ses démonstrations.
Comme lectures complémentaires, la référence standard pour la théorie en dimension finie est [FK94] ; K. McCrimmon, *A Taste of Jordan Algebras* (Springer 2004), couvre l'algèbre et l'histoire, et
H. Hanche-Olsen et E. Størmer, *Jordan Operator Algebras* (Pitman 1984), ainsi qu'Alfsen–Shultz (2003), la théorie en
dimension infinie et sa relation avec les C\*-algèbres. Pour le rôle des systèmes composés dans les reconstructions de
la théorie quantique à partir d'axiomes opérationnels, voir les articles de Barnum et Wilce cités plus haut. Pour les
inégalités de traces au-delà du chapitre 11, E. Carlen, *Trace inequalities and quantum entropy* (Contemp. Math. 529,
2010). La représentation quadratique et le cône de Lorentz réapparaissent en optimisation conique : Alizadeh–Goldfarb
(2003) et les méthodes de points intérieurs jordaniennes qui y sont passées en revue.
