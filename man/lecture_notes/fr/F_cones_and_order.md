[English](../F_cones_and_order.md) · [Español](../es/F_cones_and_order.md) · **Français** · [Deutsch](../de/F_cones_and_order.md) · [中文](../zh/F_cones_and_order.md) · [日本語](../ja/F_cones_and_order.md) · [Русский](../ru/F_cones_and_order.md)

# Annexe F. Cônes symétriques, isomorphismes d'ordre et applications de Jordan

Le chapitre 13 a relié l'ordre des observables à leur produit de Jordan par trois théorèmes : le théorème de
Koecher–Vinberg (théorème 13.8), le théorème de Kadison (théorème 13.10) et la description des automorphismes de
Jordan de $`M_n(\mathbb C)`$ (la proposition 13.11 la démontre pour $`n=2`$). Cette annexe démontre le premier en toute
généralité et les deux autres pour les algèbres $`\bigoplus_cM_{k_c}(\mathbb C)`$ du cours, pour toutes les tailles de
blocs. Chemin faisant, elle démontre les faits sur les algèbres de Jordan euclidiennes dont la section 13.6 a besoin :
un théorème spectral, l'autodualité du cône des carrés et $`P(y)Q\subseteq Q`$.

**Vous apprendrez** pourquoi toute algèbre de Jordan euclidienne est à puissances associatives et possède un théorème
spectral ; pourquoi son cône des carrés est autodual et homogène ; comment on retrouve un produit de Jordan à partir du
seul cône autodual homogène, par la dérivée troisième du logarithme de sa fonction caractéristique ; pourquoi une
bijection linéaire unifère qui préserve la positivité dans les deux sens préserve les carrés ; et pourquoi toute
application de ce type de $`\bigoplus_cM_{k_c}`$ permute des blocs égaux et agit sur chaque bloc par $`UxU^*`$ ou par
$`Ux^TU^*`$.

*Calcul :* les vérifications numériques de cette annexe — [Tutoriel, annexe F](../../tutorial/fr/F_cones_and_order.md).

Outre l'algèbre linéaire, le calcul différentiel et des résultats des chapitres 1, 6 et 13 cités par leur numéro, la
section F.3 utilise l'intégration de Lebesgue (le lemme de Fatou, la convergence dominée et la dérivation sous le signe
intégrale, Fubini–Tonelli, la formule de changement de variables linéaire) ; rien d'autre n'est supposé. Le théorème
spectral de l'annexe E, avec les repères de Jordan, raffine le théorème F.3 ci-dessous ; rien ici n'en dépend. Un
**espace euclidien** est un espace vectoriel réel de dimension finie $`V\ne\{0\}`$ muni d'un produit scalaire
$`\langle\cdot,\cdot\rangle`$ ; $`\|x\|`$ est la norme (les barres simples $`|x|`$ d'un élément de $`V`$ désigneront la valeur
absolue du théorème F.3, jamais la norme), $`\|A\|`$ la norme d'opérateur d'une application linéaire $`A`$ de $`V`$, $`A^T`$ son
adjointe ($`\langle Ax,y\rangle=\langle x,A^Ty\rangle`$), et $`A\ge0`$ ($`A>0`$) signifie que $`A=A^T`$ est semi-définie positive
(définie positive). $`e^A=\exp A=\sum_{k\ge0}A^k/k!`$.

<a id="s-F-1"></a>

## F.1 Algèbres de Jordan euclidiennes

Comme à la section 13.1, une **algèbre de Jordan** est un espace vectoriel réel muni d'un produit bilinéaire et
commutatif $`\circ`$ tel que $`x\circ(x^2\circ y)=x^2\circ(x\circ y)`$, c'est-à-dire $`[L(x),L(x^2)]=0`$, où $`L(x)y=x\circ y`$ et
$`x^2=x\circ x`$ ; et $`P(x)=2L(x)^2-L(x^2)`$.

**Définition.** Une **algèbre de Jordan euclidienne** est une algèbre de Jordan $`V`$ de dimension finie munie d'une
unité $`e`$ ($`e\circ x=x`$) et d'un produit scalaire pour lequel tout $`L(x)`$ est symétrique : $`\langle x\circ y,z\rangle=\langle y,x\circ z\rangle`$.
Son **cône des carrés** est $`Q=\{x^2:x\in V\}`$.

Exemples : $`M_{sa}`$ (la partie autoadjointe de $`\bigoplus_cM_{k_c}(\mathbb C)`$) et $`\mathrm{Herm}_n(\mathbb R)`$
avec $`\langle x,y\rangle=\mathrm{Tr}(xy)`$, qui est réel et vérifie $`\mathrm{Tr}(x^2)=\sum_{ij}|x_{ij}|^2`$ (lemme 13.1 ;
$`\mathrm{Tr}((x\circ y)z)=\frac12\mathrm{Tr}(xyz+yxz)=\mathrm{Re}\,\mathrm{Tr}(xyz)`$, puisque
$`\mathrm{Tr}(yxz)=\overline{\mathrm{Tr}((yxz)^*)}=\overline{\mathrm{Tr}(zxy)}=\overline{\mathrm{Tr}(xyz)}`$,
est invariant par permutation circulaire et, pour la même raison, par $`x\leftrightarrow z`$, donc symétrique en
$`x,y,z`$) ; les facteurs de spin $`V_m`$ (algèbres de Jordan d'après le lemme 13.3) avec le produit scalaire
$`\langle x,y\rangle=\mathrm{tr}(x\circ y)=2(ab+u\cdot v)`$ pour $`x=(a,u)`$, $`y=(b,v)`$, puisque
$`\mathrm{tr}((x\circ y)\circ z)=2(abc+c\,u\cdot v+a\,v\cdot w+b\,u\cdot w)`$ pour $`z=(c,w)`$ est symétrique en
$`x,y,z`$ ; et $`\mathbb R^N`$ avec le produit coordonnée par coordonnée (exercice 1).

<a id="r-F-1"></a>

**Lemme F.1.** Dans toute algèbre de Jordan, pour tous $`a,b,d,x,y`$ :

1. $`[L(a),L(b\circ d)]+[L(b),L(d\circ a)]+[L(d),L(a\circ b)]=0`$ ;
2. $`L(x^2\circ y)=L(x^2)L(y)-2L(x)L(y)L(x)+2L(x\circ y)L(x)`$.

*Démonstration.* (1) Pour $`s,t,r\in\mathbb R`$ et $`x=sa+tb+rd`$, $`[L(x),L(x^2)]=0`$ est une identité polynomiale en
$`(s,t,r)`$, donc chaque coefficient est nul ; le coefficient de $`str`$ est le double du membre de gauche de (1). (2) Posons
$`a=b=x`$ dans (1) : $`2[L(x),L(x\circ d)]+[L(d),L(x^2)]=0`$. Appliqué à $`y`$, cela s'écrit
$`2x\circ((x\circ d)\circ y)-2(x\circ d)\circ(x\circ y)+d\circ(x^2\circ y)-x^2\circ(d\circ y)=0`$, c'est-à-dire
$`\big(2L(x)L(y)L(x)-2L(x\circ y)L(x)+L(x^2\circ y)-L(x^2)L(y)\big)d=0`$ pour tout $`d`$. $`\square`$

<a id="r-F-2"></a>

**Proposition F.2 (associativité des puissances).** Soient $`V`$ une algèbre de Jordan d'unité $`e`$, $`x\in V`$, $`x^0=e`$ et
$`x^{k+1}=x\circ x^k`$. Alors $`x^p\circ x^q=x^{p+q}`$ pour tous $`p,q\ge0`$, et tout $`L(x^p)`$ est un polynôme en les
opérateurs $`L(x)`$ et $`L(x^2)`$, qui commutent. Par conséquent $`p(x)\circ q(x)=(pq)(x)`$ pour des polynômes réels $`p,q`$, et
les opérateurs $`L(p(x))`$ commutent entre eux.

*Démonstration.* $`L(x)`$ et $`L(x^2)`$ commutent d'après l'identité de Jordan. Montrons par récurrence sur $`N\ge2`$ : (A$`_N`$)
$`x^i\circ x^j=x^{i+j}`$ pour $`i+j\le N`$, et (B$`_N`$) $`L(x^i)`$ est un polynôme en $`L(x),L(x^2)`$ pour $`i\le N`$. Les deux
sont clairs pour $`N=2`$. Supposons-les vraies pour $`N`$, et soit $`i+j=N+1`$ avec $`i,j\ge1`$ (pour $`i=0`$ ou $`j=0`$ il n'y a
rien à démontrer). Comme $`i\le N`$, $`L(x^i)`$ commute avec $`L(x)`$, donc d'après (A$`_N`$)
```math
x^i\circ x^j=L(x^i)L(x)x^{j-1}=L(x)L(x^i)x^{j-1}=x\circ x^N=x^{N+1}.
```
C'est (A$`_{N+1}`$). Le lemme F.1(2) avec $`y=x^{N-1}`$ et $`x^2\circ x^{N-1}=x^{N+1}`$ donne
$`L(x^{N+1})=L(x^2)L(x^{N-1})-2L(x)L(x^{N-1})L(x)+2L(x^N)L(x)`$, ce qui est (B$`_{N+1}`$). La dernière phrase s'ensuit par
bilinéarité. $`\square`$

<a id="r-F-3"></a>

**Théorème F.3 (théorème spectral).** Soient $`V`$ une algèbre de Jordan euclidienne et $`x\in V`$. Il existe un unique
$`k\ge1`$, d'uniques réels $`\lambda_1>\dots>\lambda_k`$ et d'uniques $`c_1,\dots,c_k\in V`$ non nuls tels que
```math
c_i\circ c_j=\delta_{ij}c_i,\qquad \textstyle\sum_ic_i=e,\qquad x=\sum_i\lambda_ic_i .
```
Les $`c_i`$ sont des polynômes en $`x`$, et $`p(x)=\sum_ip(\lambda_i)c_i`$ pour tout polynôme réel $`p`$.

On appelle $`(c_i)`$ un **système complet d'idempotents orthogonaux** et les $`\lambda_i`$ les **valeurs propres** de $`x`$,
et l'on pose $`f(x)=\sum_if(\lambda_i)c_i`$ pour $`f:\{\lambda_1,\dots,\lambda_k\}\to\mathbb R`$ ; par exemple $`|x|`$,
$`x_\pm=\max(\pm x,0)`$, $`\sqrt x`$ (si tous les $`\lambda_i\ge0`$), $`\exp x`$, et $`\log x`$ (si tous les $`\lambda_i>0`$).

*Démonstration.* Existence. $`\mathbb R[x]=\{p(x)\}`$ est un sous-espace qui contient $`e`$ et qui, d'après la proposition
F.2, vérifie $`p(x)\circ q(x)=(pq)(x)`$. Ainsi $`M=L(x)|_{\mathbb R[x]}`$ envoie $`\mathbb R[x]`$ dans lui-même, est symétrique
et vérifie $`p(M)q(x)=p(x)\circ q(x)`$. Un opérateur symétrique est diagonalisable à valeurs propres réelles ; soient
$`\lambda_1>\dots>\lambda_k`$ celles de $`M`$ et $`m(t)=\prod_i(t-\lambda_i)`$ son polynôme minimal. Alors $`m(x)=m(M)e=0`$. Avec
les polynômes de Lagrange $`\ell_i(t)=\prod_{j\ne i}\frac{t-\lambda_j}{\lambda_i-\lambda_j}`$, les polynômes
$`\ell_i\ell_j-\delta_{ij}\ell_i`$, $`\sum_i\ell_i-1`$ et $`\sum_i\lambda_i\ell_i-t`$ s'annulent en chaque $`\lambda_l`$, donc sont
des multiples de $`m`$ ; en les évaluant en $`x`$, on voit que les $`c_i=\ell_i(x)`$ ont les trois propriétés. Si $`c_i=0`$, alors
$`\ell_i(M)q(x)=\ell_i(x)\circ q(x)=0`$ pour tout $`q`$, donc $`\ell_i(M)=0`$ alors que $`\deg\ell_i<\deg m`$ ; ainsi $`c_i\ne0`$.
Enfin $`p(x)=\sum_ip(\lambda_i)c_i`$ : c'est vrai pour $`p=1`$ et $`p=t`$, et $`x\circ\sum_i\lambda_i^rc_i=\sum_i\lambda_i^{r+1}c_i`$.

Unicité. Soit $`x=\sum_j\mu_jd_j`$ avec des $`\mu_j`$ distincts et un système complet d'idempotents orthogonaux non nuls $`d_j`$.
Comme ci-dessus, $`p(x)=\sum_jp(\mu_j)d_j`$. Les $`d_j`$ sont linéairement indépendants (multiplier $`\sum_ja_jd_j=0`$ par
$`d_i`$), donc $`p(x)=0`$ si et seulement si $`p(\mu_j)=0`$ pour tout $`j`$. Ainsi $`m(x)=0`$ donne $`\{\mu_j\}\subseteq\{\lambda_i\}`$,
et $`\tilde m(t)=\prod_j(t-\mu_j)`$ vérifie $`\tilde m(M)q(x)=\tilde m(x)\circ q(x)=0`$, donc $`m`$ divise $`\tilde m`$ et $`\{\lambda_i\}\subseteq\{\mu_j\}`$.
Alors $`d_j=\ell_j(x)=c_j`$. $`\square`$

*Dans la bibliothèque :* la décomposition spectrale et les fonctions $`f(x)`$ dans un facteur de spin, `SpinFactor.spectral_decomposition`, `SpinFactor.apply_function`, `SpinFactor.exp`, `SpinFactor.log` — [tutoriel, section F.2](../../tutorial/fr/F_cones_and_order.md#s-F-2).

<a id="r-F-4"></a>

**Lemme F.4 (décomposition de Peirce).** Soient $`V`$ une algèbre de Jordan euclidienne et $`c`$ un idempotent ($`c^2=c`$).

1. $`L(c)(2L(c)-I)(L(c)-I)=0`$. Par conséquent $`V`$ est la somme orthogonale des sous-espaces propres
   $`V(c,\nu)=\{z:c\circ z=\nu z\}`$, $`\nu\in\{0,\frac12,1\}`$, et $`L(c)\ge0`$.
2. $`V(c,1)\circ V(c,0)=\{0\}`$ ; $`V(c,\nu)\circ V(c,\nu)\subseteq V(c,\nu)`$ et $`V(c,\nu)\circ V(c,\frac12)\subseteq V(c,\frac12)`$
   pour $`\nu\in\{0,1\}`$ ; $`V(c,\frac12)\circ V(c,\frac12)\subseteq V(c,0)+V(c,1)`$.
3. $`s=2c-e`$ vérifie $`s^2=e`$, et $`P(s)`$ est un automorphisme de Jordan de $`V`$ : il vaut $`+1`$ sur $`V(c,0)+V(c,1)`$ et $`-1`$
   sur $`V(c,\frac12)`$, et $`P(s)(z\circ w)=P(s)z\circ P(s)w`$.
4. Soit $`x=\sum_{i=1}^k\lambda_ic_i`$ comme dans le théorème F.3. Les $`L(c_i)`$ commutent, et $`V`$ est la somme orthogonale
   des sous-espaces $`V_{ii}=\{z:c_i\circ z=z\}`$ et $`V_{ij}=\{z:c_i\circ z=c_j\circ z=\frac12z\}`$ ($`i<j`$). Pour
   $`a=\sum_la_lc_l`$ ($`a_l\in\mathbb R`$), $`L(a)`$ agit sur $`V_{ij}`$ ($`i\le j`$) comme le scalaire $`\frac12(a_i+a_j)`$, et $`P(a)`$
   comme le scalaire $`a_ia_j`$.

*Démonstration.* (1) Le lemme F.1(2) avec $`x=y=c`$ donne $`L(c)=L(c)^2-2L(c)^3+2L(c)^2`$, c'est-à-dire
$`2L(c)^3-3L(c)^2+L(c)=0`$, ce qui est la factorisation annoncée. Les valeurs propres de l'opérateur symétrique $`L(c)`$
sont donc dans $`\{0,\frac12,1\}`$.

(2) Le lemme F.1(1) avec $`a=b=c`$, $`d=z\in V(c,\nu)`$ donne $`2\nu[L(c),L(z)]+[L(z),L(c)]=0`$, donc $`L(z)`$ commute avec $`L(c)`$
si $`\nu\in\{0,1\}`$ et préserve alors chaque $`V(c,\mu)`$. Cela donne la deuxième et la troisième inclusion, et
$`V(c,1)\circ V(c,0)\subseteq V(c,0)\cap V(c,1)=\{0\}`$. Pour $`z,w\in V(c,\frac12)`$, appliquons le lemme F.1(1) avec $`(a,b,d)=(z,w,c)`$
à $`c`$ : comme $`L(w\circ c)=\frac12L(w)`$ et $`L(c\circ z)=\frac12L(z)`$, les deux premiers termes se compensent, et
$`[L(c),L(u)]c=L(c)^2u-L(c)u=0`$ pour $`u=z\circ w`$. En écrivant $`u=u_1+u_{1/2}+u_0`$ selon (1), cela s'écrit $`-\frac14u_{1/2}=0`$.

(3) $`s^2=4c-4c+e=e`$, $`L(s)=2L(c)-I`$ vaut $`1,0,-1`$ sur $`V(c,1),V(c,\frac12),V(c,0)`$, et $`P(s)=2L(s)^2-I`$. Écrivons
$`z=z'+z''`$, $`w=w'+w''`$ avec $`z',w'\in V(c,0)+V(c,1)`$ et $`z'',w''\in V(c,\frac12)`$. D'après (2), $`z'\circ w'+z''\circ w''`$
appartient à $`V(c,0)+V(c,1)`$ et $`z'\circ w''+z''\circ w'`$ à $`V(c,\frac12)`$, donc
$`P(s)(z\circ w)=(z'-z'')\circ(w'-w'')=P(s)z\circ P(s)w`$.

(4) Les $`c_l`$ sont des polynômes en $`x`$, donc les $`L(c_l)`$ commutent (proposition F.2) ; étant symétriques, ils ont une
base orthonormée commune de vecteurs propres. Sur un vecteur propre commun, $`L(c_l)`$ a une valeur propre
$`\varepsilon_l\in\{0,\frac12,1\}`$ d'après (1), et $`\sum_l\varepsilon_l=1`$ car $`\sum_lL(c_l)=L(e)=I`$. Donc, ou bien un
$`\varepsilon_i=1`$ et les autres sont nuls, ou bien $`\varepsilon_i=\varepsilon_j=\frac12`$ pour une paire $`i<j`$ et les autres
sont nuls. Ainsi $`V`$ est la somme orthogonale des sous-espaces propres communs $`W_{ii}`$ (où $`L(c_i)=1`$, $`L(c_l)=0`$ pour
$`l\ne i`$) et $`W_{ij}`$ (où $`L(c_i)=L(c_j)=\frac12`$, $`L(c_l)=0`$ sinon), et $`W_{ij}=V_{ij}`$ : un vecteur tel que $`c_i\circ z=z`$
(resp. $`c_i\circ z=c_j\circ z=\frac12z`$) n'a pas de composante dans les autres sous-espaces propres communs, sur lesquels
$`L(c_i)`$ (resp. la paire $`(L(c_i),L(c_j))`$) prend une autre valeur. Comme
$`a^2=\sum_la_l^2c_l`$, sur $`V_{ij}`$ ($`i<j`$) on a $`L(a)=\frac12(a_i+a_j)`$, $`L(a^2)=\frac12(a_i^2+a_j^2)`$ et
$`P(a)=\frac12(a_i+a_j)^2-\frac12(a_i^2+a_j^2)=a_ia_j`$ ; sur $`V_{ii}`$, $`L(a)=a_i`$ et $`P(a)=2a_i^2-a_i^2`$. $`\square`$

<a id="s-F-2"></a>

## F.2 Le cône des carrés est symétrique

Soit $`\Omega`$ un **cône convexe** ouvert non vide d'un espace euclidien $`V`$ ($`tx\in\Omega`$ pour $`x\in\Omega`$, $`t>0`$). Son
**cône dual** (ouvert) est $`\Omega^*=\{y:\langle x,y\rangle>0\text{ pour tout }x\in\bar\Omega\setminus\{0\}\}`$ ; $`\Omega`$ est
**autodual** si $`\Omega^*=\Omega`$. $`G(\Omega)=\{g\in GL(V):g\Omega=\Omega\}`$, et $`\Omega`$ est **homogène** si
$`G(\Omega)`$ agit transitivement sur $`\Omega`$. Un **cône symétrique** est un cône autodual homogène. Si
$`\Omega_1\subseteq\Omega_2`$, alors $`\Omega_2^*\subseteq\Omega_1^*`$.

<a id="r-F-5"></a>

**Théorème F.5 (autodualité).** Soient $`V`$ une algèbre de Jordan euclidienne et $`y=\sum_i\lambda_ic_i`$ (théorème F.3).

1. Les conditions suivantes sont équivalentes : (i) $`y\in Q`$ ; (ii) tous les $`\lambda_i\ge0`$ ; (iii) $`L(y)\ge0`$ ; (iv)
   $`\langle y,x\rangle\ge0`$ pour tout $`x\in Q`$. En particulier, $`Q`$ est un cône convexe fermé.
2. Les conditions suivantes sont équivalentes : (i) $`y\in\Omega:=\mathrm{int}\,Q`$ ; (ii) tous les $`\lambda_i>0`$ ; (iii)
   $`L(y)>0`$ ; (iv) $`\langle y,x\rangle>0`$ pour tout $`x\in Q\setminus\{0\}`$. De plus $`e\in\Omega`$, $`\bar\Omega=Q`$ et
   $`\Omega^*=\Omega`$ : $`\Omega`$ est un cône convexe ouvert autodual.
3. Si $`x,y\in Q`$ et $`\langle x,y\rangle=0`$, alors $`x\circ y=0`$.
4. $`x^2=0`$ entraîne $`x=0`$, et $`Q\cap(-Q)=\{0\}`$.

*Démonstration.* On utilise $`\langle y,x^2\rangle=\langle y\circ x,x\rangle=\langle L(y)x,x\rangle`$ et
$`\langle x^2,e\rangle=\langle x,x\circ e\rangle=\|x\|^2`$ ; la seconde donne (4) avec (1)(ii), puisque les valeurs propres
de $`-y`$ sont les $`-\lambda_i`$ par unicité.

(1) (i)$`\Rightarrow`$(iii) : si $`y=x^2`$ avec $`x=\sum_i\mu_id_i`$, alors $`y=\sum_i\mu_i^2d_i`$ et
$`L(y)=\sum_i\mu_i^2L(d_i)\ge0`$ d'après le lemme F.4(1). (iii)$`\Leftrightarrow`$(iv) par la première identité.
(iii)$`\Rightarrow`$(ii) : $`0\le\langle L(y)c_i,c_i\rangle=\lambda_i\|c_i\|^2`$. (ii)$`\Rightarrow`$(i) :
$`y=(\sum_i\sqrt{\lambda_i}c_i)^2`$. D'après (iv), $`Q`$ est une intersection de demi-espaces fermés passant par $`0`$.

(2) (ii)$`\Rightarrow`$(iii) : $`L(y)=\sum_i\lambda_iL(c_i)\ge(\min_i\lambda_i)\sum_iL(c_i)=(\min_i\lambda_i)I`$.
(iii)$`\Rightarrow`$(iv) : un $`x\in Q`$ non nul est $`z^2`$ avec $`z\ne0`$, et $`\langle y,x\rangle=\langle L(y)z,z\rangle>0`$.
(iv)$`\Rightarrow`$(ii) : $`\langle y,c_i\rangle=\lambda_i\|c_i\|^2`$ et $`c_i=c_i^2\in Q\setminus\{0\}`$. (iii)$`\Rightarrow`$(i) : si
$`L(y)\ge mI`$ avec $`m>0`$, alors $`L(y')\ge(m-\|L\|\,\|y-y'\|)I`$, où $`\|L\|`$ est la norme de l'application linéaire
$`b\mapsto L(b)`$, donc un voisinage de $`y`$ est formé d'éléments tels que
$`L>0`$, qui sont dans $`Q`$ d'après (1). (i)$`\Rightarrow`$(iii) : $`L(y)\ge0`$ d'après (1) ; si $`\langle L(y)z,z\rangle=0`$ avec
$`z\ne0`$, alors $`x=z^2\ne0`$ et $`\langle y,x\rangle=0`$ ; pour $`\varepsilon>0`$ petit, $`y-\varepsilon x\in Q`$ mais
$`\langle y-\varepsilon x,x\rangle=-\varepsilon\|x\|^2<0`$, ce qui contredit (1)(iv). Ensuite, $`L(e)=I>0`$, et pour $`y\in Q`$,
$`y+\varepsilon e\in\Omega`$ tend vers $`y`$, donc $`\bar\Omega=Q`$ ; alors $`\Omega^*=\{y:\langle y,x\rangle>0\ \forall x\in Q\setminus\{0\}\}=\Omega`$
d'après (iv). $`\Omega=\{L>0\}`$ est convexe.

(3) Soit $`x=\sum_i\lambda_ic_i`$ avec $`\lambda_i\ge0`$. Alors $`0=\sum_i\lambda_i\langle c_i,y\rangle`$, tous les termes étant
$`\ge0`$ d'après (1). Si $`\lambda_i>0`$, alors $`0=\langle y,c_i\rangle=\langle L(y)c_i,c_i\rangle`$ avec $`L(y)\ge0`$, donc
$`L(y)c_i=0`$. Par conséquent $`x\circ y=\sum_i\lambda_i\,c_i\circ y=0`$. $`\square`$

*Calcul :* (1) dans $`\mathrm{Herm}_4(\mathbb C)`$ et pour le cône de Lorentz de $`V_4`$ (le test d'appartenance au cône `SpinFactor.in_cone`) — [tutoriel, section F.2](../../tutorial/fr/F_cones_and_order.md#s-F-2).

<a id="r-F-6"></a>

**Théorème F.6 (homogénéité et représentation quadratique).** Soient $`V`$ une algèbre de Jordan euclidienne, $`Q`$ et
$`\Omega`$ comme dans le théorème F.5.

1. Pour tout $`b\in V`$ : $`e^{L(b)}Q=Q`$, $`e^{L(b)}\Omega=\Omega`$ et $`e^{L(b)}e=\exp b`$.
2. Tout $`y\in\Omega`$ est égal à $`\exp(\log y)`$. Par conséquent $`\Omega`$ est homogène, et c'est un cône symétrique.
3. Pour $`y\in\Omega`$, $`P(y)=e^{2L(\log y)}`$ ; donc $`P(y)Q=Q`$, $`P(y)\Omega=\Omega`$ et $`P(y)e=y^2`$.
4. $`P(y)Q\subseteq Q`$ pour tout $`y\in V`$.

*Démonstration.* (1) Soit $`0\le s<1/\|L(b)\|`$ (n'importe quel $`s\ge0`$ si $`L(b)=0`$). L'opérateur symétrique $`I-sL(b)`$ a
des valeurs propres $`\ge1-s\|L(b)\|>0`$. Soient $`x\in Q`$ et $`z=(I-sL(b))^{-1}x=\sum_j\nu_jd_j`$ (théorème F.3). En utilisant
$`z\circ d_j=\nu_jd_j`$,
$`\langle d_l,d_j\rangle=\langle d_l,d_j\circ d_j\rangle=\langle d_l\circ d_j,d_j\rangle=\delta_{lj}\|d_j\|^2`$ et
$`\langle b,d_j\rangle=\langle b,d_j^2\rangle=\langle L(b)d_j,d_j\rangle\le\|L(b)\|\,\|d_j\|^2`$,
```math
0\le\langle x,d_j\rangle=\langle z-s\,b\circ z,d_j\rangle=\nu_j\|d_j\|^2-s\langle b,z\circ d_j\rangle
=\nu_j\big(\|d_j\|^2-s\langle b,d_j\rangle\big),
```
où l'inégalité de gauche est le théorème F.5(1)(iv). La parenthèse est $`\ge(1-s\|L(b)\|)\|d_j\|^2>0`$, donc tous les
$`\nu_j\ge0`$ et $`z\in Q`$. Ainsi, pour $`t\ge0`$ et $`n>t\|L(b)\|`$, $`(I-\frac tnL(b))^{-n}Q\subseteq Q`$. Quand $`n\to\infty`$, cet
opérateur tend vers $`e^{tL(b)}`$ (dans une base orthonormée de vecteurs propres de $`L(b)`$ :
$`(1-t\beta/n)^{-n}\to e^{t\beta}`$), et $`Q`$ est fermé, donc $`e^{tL(b)}Q\subseteq Q`$. Avec $`-b`$ à la place de $`b`$,
$`e^{-L(b)}Q\subseteq Q`$, donc $`e^{L(b)}Q=Q`$, et l'homéomorphisme $`e^{L(b)}`$ envoie $`\mathrm{int}\,Q`$ sur lui-même.
Enfin $`L(b)^re=b^r`$ (puissances de Jordan), donc $`e^{L(b)}e=\sum_rb^r/r!=\sum_ie^{\beta_i}d_i=\exp b`$ pour
$`b=\sum_i\beta_id_i`$, d'après le théorème F.3.

(2) Si $`y=\sum_i\lambda_ic_i`$ avec $`\lambda_i>0`$, alors $`\log y=\sum_i(\log\lambda_i)c_i`$ et $`\exp\log y=y`$. Ainsi
$`y=e^{L(\log y)}e`$ avec $`e^{L(\log y)}\in G(\Omega)`$ d'après (1) ; l'autodualité est le théorème F.5(2).

(3) D'après le lemme F.4(4) appliqué aux $`c_i`$ de $`y`$, $`P(y)`$ et $`e^{2L(\log y)}`$ agissent tous deux sur $`V_{ij}`$ comme
$`\lambda_i\lambda_j=e^{\log\lambda_i+\log\lambda_j}`$. On utilise (1), et $`P(y)e=2y\circ y-y^2=y^2`$.

(4) Soient $`y=\sum_i\lambda_ic_i`$, $`c=\sum_{\lambda_i\ge0}c_i`$ (un idempotent), $`s=2c-e=\sum_i\varepsilon_ic_i`$ avec
$`\varepsilon_i=\pm1`$, et $`|y|=\sum_i|\lambda_i|c_i`$. D'après le lemme F.4(4), $`P(s)`$, $`P(|y|)`$ et $`P(y)`$ agissent sur $`V_{ij}`$
comme $`\varepsilon_i\varepsilon_j`$, $`|\lambda_i\lambda_j|`$ et $`\lambda_i\lambda_j`$, donc $`P(y)=P(s)P(|y|)`$. Pour $`\eta>0`$,
$`|y|+\eta e\in\Omega`$, donc $`P(|y|+\eta e)Q\subseteq Q`$ d'après (3) ; $`P(\cdot)`$ est polynomial et $`Q`$ est fermé, donc
$`P(|y|)Q\subseteq Q`$ ($`\eta\to0`$). $`P(s)`$ est un automorphisme de Jordan (lemme F.4(3)), donc $`P(s)x^2=(P(s)x)^2\in Q`$. $`\square`$

*Calcul :* (1), (3) et (4) dans $`V_4`$, avec `SpinFactor.product` et `SpinFactor.quadratic_representation` — [tutoriel, section F.2](../../tutorial/fr/F_cones_and_order.md#s-F-2).

La section 13.6 a démontré $`\mathrm{tr}|P(y)x|\le\mathrm{tr}\,P(y)|x|`$ dans $`V_m`$ à partir de trois propriétés.
Elles sont vraies dans toute algèbre de Jordan euclidienne pour la forme linéaire $`\tau(x)=\langle x,e\rangle`$, qui est
$`\mathrm{Tr}`$ sur $`M_{sa}`$ et $`\mathrm{tr}`$ sur $`V_m`$ (avec les produits scalaires ci-dessus) :

<a id="r-F-7"></a>

**Corollaire F.7 (proposition 13.14 dans toute algèbre de Jordan euclidienne).** Pour $`x,y`$ dans une algèbre de Jordan
euclidienne, $`\tau(|P(y)x|)\le\tau(P(y)|x|)`$, avec égalité si $`x\in Q`$ ou $`-x\in Q`$.

*Démonstration.* $`x_\pm\in Q`$, $`x=x_+-x_-`$, $`|x|=x_++x_-`$ (théorème F.3). $`A=P(y)x_+`$ et $`B=P(y)x_-`$ appartiennent à $`Q`$
(théorème F.6(4)), et $`z=P(y)x=A-B`$. Soient $`z=\sum_i\mu_ic_i`$ et $`c=\sum_{\mu_i>0}c_i`$ ; alors $`|z|=z\circ(2c-e)`$, et
$`\tau(a\circ b)=\langle a,b\rangle`$. Comme $`\langle Q,Q\rangle\ge0`$ (théorème F.5(1)) et $`c,e-c\in Q`$,
```math
\tau(|z|)=\langle A,c\rangle-\langle A,e-c\rangle-\langle B,c\rangle+\langle B,e-c\rangle\le\tau(A)+\tau(B)=\tau(P(y)|x|).
```
Si $`x\in Q`$, alors $`x_-=0`$ et $`z=A\in Q`$, donc $`|z|=z=P(y)|x|`$ ; de même si $`-x\in Q`$. $`\square`$

*Calcul :* l'inégalité dans $`\mathrm{Herm}_5(\mathbb R)`$ — [tutoriel, section F.2](../../tutorial/fr/F_cones_and_order.md#s-F-2).

<a id="s-F-3"></a>

## F.3 Le théorème de Koecher–Vinberg

La section F.2 a démontré que l'intérieur du cône des carrés d'une algèbre de Jordan euclidienne est symétrique. La
réciproque reconstruit le produit de Jordan à partir du cône. Jusqu'au théorème F.12, $`\Omega`$ est un cône convexe
ouvert autodual d'un espace euclidien $`V`$ de dimension $`N`$, et $`dy`$ est la mesure de Lebesgue en coordonnées
orthonormées (les changements de coordonnées orthogonaux vérifient $`|\det|=1`$).

<a id="r-F-8"></a>

**Lemme F.8.** Soit $`\Omega`$ autodual.

1. $`0\notin\Omega`$ ; $`\langle x,y\rangle\ge0`$ pour $`x,y\in\bar\Omega`$ ; pour tout $`x\in\Omega`$, il existe $`\delta>0`$ tel que
   $`\langle x,y\rangle\ge\delta\|y\|`$ pour tout $`y\in\bar\Omega`$.
2. Si $`x\in\bar\Omega\setminus\Omega`$, il existe $`u\in\bar\Omega\setminus\{0\}`$ tel que $`\langle x,u\rangle=0`$. Si $`y\in\Omega`$,
   $`u\in\bar\Omega`$ et $`t\ge0`$, alors $`y+tu\in\Omega`$.
3. $`\Omega`$ est l'intérieur de $`\bar\Omega`$.
4. $`G(\Omega)`$ est un sous-groupe fermé de $`GL(V)`$, et $`g\in G(\Omega)`$ entraîne $`g^T\in G(\Omega)`$.

*Démonstration.* (1) Si $`0\in\Omega`$, alors $`\Omega`$ contient une boule centrée en $`0`$, donc $`\Omega=V`$ ; mais
$`V^*=\emptyset`$. Pour $`x,y\in\Omega`$, $`\langle x,y\rangle>0`$ car $`y\in\Omega^*`$ ; par continuité, $`\langle x,y\rangle\ge0`$
sur $`\bar\Omega`$. Si la boule fermée $`\bar B(x,\delta)`$ est contenue dans $`\Omega`$ et $`y\in\bar\Omega\setminus\{0\}`$, alors
$`\langle x-\delta y/\|y\|,y\rangle\ge0`$. (2) $`x\notin\Omega=\Omega^*`$ fournit $`u\in\bar\Omega\setminus\{0\}`$ tel que
$`\langle u,x\rangle\le0`$, et $`\langle u,x\rangle\ge0`$ d'après (1). Pour $`x'\in\bar\Omega\setminus\{0\}`$,
$`\langle x',y+tu\rangle\ge\langle x',y\rangle>0`$, donc $`y+tu\in\Omega^*=\Omega`$. (3) Si $`\bar B(y,\varepsilon)\subseteq\bar\Omega`$,
alors pour $`x\in\bar\Omega\setminus\{0\}`$, $`\langle x,y-\varepsilon x/\|x\|\rangle\ge0`$, donc $`\langle x,y\rangle\ge\varepsilon\|x\|>0`$
et $`y\in\Omega^*=\Omega`$. (4) Pour $`g\in G(\Omega)`$, $`g`$ envoie $`\bar\Omega\setminus\{0\}`$ sur lui-même, donc pour $`y\in\Omega`$
et $`x\in\bar\Omega\setminus\{0\}`$, $`\langle x,g^Ty\rangle=\langle gx,y\rangle>0`$ : $`g^T\Omega\subseteq\Omega`$. Le même argument
pour $`g^{-1}`$ donne $`g^T\Omega=\Omega`$. Si des $`g_k\in G(\Omega)`$ convergent vers $`g\in GL(V)`$, alors $`g\Omega\subseteq\bar\Omega`$
est ouvert, donc $`g\Omega\subseteq\Omega`$ d'après (3) ; de même $`g^{-1}=\lim g_k^{-1}`$ envoie $`\Omega`$ dans $`\Omega`$, donc
$`g\Omega=\Omega`$. $`\square`$

La **fonction caractéristique** de $`\Omega`$ est $`\varphi(x)=\int_{\Omega^*}e^{-\langle x,y\rangle}dy=\int_\Omega e^{-\langle x,y\rangle}dy`$.
Posons $`f=\log\varphi`$, soit $`H(x)`$ l'opérateur symétrique tel que $`\langle H(x)u,v\rangle=D^2f(x)[u,v]`$, et
$`x^*=-\nabla f(x)`$, c'est-à-dire $`Df(x)[u]=-\langle x^*,u\rangle`$.

<a id="r-F-9"></a>

**Lemme F.9 (fonction caractéristique).** Soit $`\Omega`$ autodual.

1. $`\varphi`$ est finie, positive et $`C^\infty`$ sur $`\Omega`$, et peut être dérivée sous le signe intégrale.
2. $`\varphi(tx)=t^{-N}\varphi(x)`$ pour $`t>0`$, et $`\varphi(gx)=|\det g|^{-1}\varphi(x)`$ pour $`g\in G(\Omega)`$.
3. $`H(x)>0`$ et $`x^*\in\Omega`$ pour $`x\in\Omega`$.
4. $`H(x)x=x^*`$ et $`D^3f(x)[x,u,v]=-2\langle H(x)u,v\rangle`$ ; pour $`g\in G(\Omega)`$, $`(gx)^*=g^{-T}x^*`$ et
   $`H(gx)=g^{-T}H(x)g^{-1}`$.
5. $`\varphi(x)\ge c\|x\|^{-N}`$ pour une constante $`c>0`$, et $`\varphi(x_k)\to\infty`$ si des $`x_k\in\Omega`$ convergent vers un
   point de $`\bar\Omega\setminus\Omega`$.

*Démonstration.* (1) Soient $`x\in\Omega`$ et $`\delta`$ comme dans le lemme F.8(1). Pour $`\|x'-x\|<\delta/2`$ et $`y\in\Omega`$,
$`\langle x',y\rangle\ge\frac\delta2\|y\|`$, donc les dérivées $`k`$-ièmes de l'intégrande en $`x'`$,
$`(-1)^k\prod_{i}\langle u_i,y\rangle e^{-\langle x',y\rangle}`$, sont majorées par $`\prod_i\|u_i\|\,\|y\|^ke^{-\delta\|y\|/2}`$, qui est
intégrable ($`\|y\|\ge N^{-1/2}\sum_i|y_i|`$ et Fubini). La convergence dominée permet de dériver sous le signe intégrale.
(2) On fait le changement de variables $`y=y'/t`$, qui envoie $`\Omega`$ sur $`\Omega`$ ; et $`y'=g^Ty`$, qui envoie $`\Omega`$ sur
$`\Omega`$ d'après le lemme F.8(4), avec $`\langle gx,y\rangle=\langle x,y'\rangle`$ et $`dy=|\det g|^{-1}dy'`$.
(3) Soit $`\mu_x`$ la mesure de probabilité $`e^{-\langle x,y\rangle}dy/\varphi(x)`$ sur $`\Omega`$. D'après (1),
$`Df(x)[u]=-\mathbb E\langle u,y\rangle`$ et $`D^2f(x)[u,u]=\mathbb E\langle u,y\rangle^2-(\mathbb E\langle u,y\rangle)^2`$, la
variance de $`\langle u,y\rangle`$. Elle ne s'annule que si $`\langle u,y\rangle`$ est égal à une constante $`\gamma`$
$`\mu_x`$-presque partout, c'est-à-dire Lebesgue-presque partout sur $`\Omega`$ (la densité de $`\mu_x`$ est strictement
positive) ; l'ensemble $`\{y\in\Omega:\langle u,y\rangle\ne\gamma\}`$ est ouvert et négligeable, donc vide, si bien que la
fonction linéaire $`\langle u,\cdot\rangle`$ est constante sur une boule, ce qui impose $`u=0`$. Et $`x^*=\mathbb E\,y`$, donc
$`\langle z,x^*\rangle=\mathbb E\langle z,y\rangle>0`$ pour $`z\in\bar\Omega\setminus\{0\}`$ : $`x^*\in\Omega^*=\Omega`$.
(4) D'après (2), $`f(tx)=f(x)-N\log t`$, donc $`Df(x)[x]=-N`$ sur $`\Omega`$. En dérivant en $`x`$ dans la direction $`u`$ :
$`D^2f(x)[x,u]+Df(x)[u]=0`$, c'est-à-dire $`\langle H(x)x,u\rangle=\langle x^*,u\rangle`$ ; puis dans la direction $`v`$ :
$`D^3f(x)[x,u,v]+2D^2f(x)[u,v]=0`$. Dériver une et deux fois $`f(gx)=f(x)-\log|\det g|`$ donne
$`\langle(gx)^*,gu\rangle=\langle x^*,u\rangle`$ et $`\langle H(gx)gu,gv\rangle=\langle H(x)u,v\rangle`$.
(5) $`\langle x,y\rangle\le\|x\|\|y\|\le1`$ pour $`\|y\|\le1/\|x\|`$, donc $`\varphi(x)\ge e^{-1}\mathrm{vol}(\Omega\cap B(0,1/\|x\|))=c\|x\|^{-N}`$
avec $`c=e^{-1}\mathrm{vol}(\Omega\cap B(0,1))>0`$. Soit $`x_k\to x_0\in\bar\Omega\setminus\Omega`$. D'après le lemme de Fatou,
$`\liminf\varphi(x_k)\ge\int_\Omega e^{-\langle x_0,y\rangle}dy`$. Prenons $`u`$ comme dans le lemme F.8(2), une boule
$`B(y_1,r)\subseteq\Omega`$ et les boules disjointes $`B_l=B(y_1+2lru/\|u\|,r)`$, $`l\ge0`$ ; elles sont contenues dans $`\Omega`$
d'après le lemme F.8(2), et sur $`B_l`$, $`\langle x_0,y\rangle=\langle x_0,y-2lru/\|u\|\rangle\le\|x_0\|(\|y_1\|+r)`$. L'intégrale est
donc infinie. $`\square`$

<a id="r-F-10"></a>

**Lemme F.10 (groupes à un paramètre).** Soient $`G`$ un sous-groupe fermé de $`GL(V)`$ et
$`\mathfrak g=\{Z:e^{tZ}\in G\text{ pour tout }t\in\mathbb R\}`$.

1. Si $`c:[0,s_0)\to G`$ vérifie $`c(s)=I+sZ+o(s)`$ quand $`s\to0^+`$, alors $`Z\in\mathfrak g`$.
2. $`\mathfrak g`$ est une algèbre de Lie : un sous-espace vectoriel stable par $`[Z,W]=ZW-WZ`$.

*Démonstration.* (1) Fixons $`t>0`$ et posons $`A_n=c(t/n)`$, $`B_n=e^{tZ/n}`$. Alors $`\|A_n-B_n\|=o(1/n)`$ et
$`\|A_n\|,\|B_n\|\le1+K/n`$ pour un certain $`K`$ et $`n`$ grand, donc
```math
\|A_n^n-e^{tZ}\|=\Big\|\sum_{r=0}^{n-1}A_n^r(A_n-B_n)B_n^{n-1-r}\Big\|\le n\,e^K\|A_n-B_n\|\to0 .
```
$`A_n^n\in G`$ et $`G`$ est fermé dans $`GL(V)`$, donc $`e^{tZ}\in G`$ ; de plus $`e^{-tZ}=(e^{tZ})^{-1}\in G`$. (2) Les multiples
sont clairs. On applique (1) à $`c(s)=e^{sZ}e^{sW}=I+s(Z+W)+O(s^2)`$, et à $`c(s)=e^{aZ}e^{aW}e^{-aZ}e^{-aW}`$ avec $`a=\sqrt s`$ :
en développant chaque facteur à l'ordre deux, les termes d'ordre $`a`$ se compensent et ceux d'ordre $`a^2`$ ont pour somme
$`ZW-WZ`$, donc $`c(s)=I+s[Z,W]+O(s^{3/2})`$. $`\square`$

<a id="r-F-11"></a>

**Proposition F.11 (le point base).** Soient $`\Omega`$ un cône symétrique, $`G=G(\Omega)`$ et $`\mathfrak g`$ comme dans le
lemme F.10.

1. $`x^{**}=x`$ et $`H(x^*)=H(x)^{-1}`$ pour $`x\in\Omega`$.
2. Il existe $`e\in\Omega`$ tel que $`e^*=e`$, et tout tel $`e`$ vérifie $`H(e)=I`$.
3. Pour $`x\in\Omega`$ et tout $`g\in G`$ tel que $`ge=x`$ : $`H(x)^{-1}=gg^T\in G`$.
4. Si $`Z\in\mathfrak g`$ et $`Ze=0`$, alors $`Z^T=-Z`$.

*Démonstration.* (1) L'application $`x\mapsto x^*`$ est lisse de $`\Omega`$ dans $`\Omega`$ (lemme F.9(3)), de dérivée $`-H(x)`$.
Comme $`g^{-T}\in G`$ (lemme F.8(4)), le lemme F.9(2,4) donne
$`f(gx)+f((gx)^*)=f(x)-\log|\det g|+f(g^{-T}x^*)=f(x)+f(x^*)`$, donc $`F(x)=f(x)+f(x^*)`$ est constante sur l'orbite
$`G x=\Omega`$. Par conséquent $`0=DF(x)[u]=-\langle x^*,u\rangle+\langle x^{**},H(x)u\rangle`$, c'est-à-dire
$`H(x)x^{**}=x^*=H(x)x`$ (lemme F.9(4)), et $`x^{**}=x`$ puisque $`H(x)`$ est inversible. En dérivant $`x^{**}=x`$, on obtient
$`H(x^*)H(x)=I`$.
(2) $`F_2(x)=f(x)+\frac12\|x\|^2`$ vérifie $`F_2(x)\ge\log c-N\log\|x\|+\frac12\|x\|^2`$ (lemme F.9(5)), qui tend vers $`\infty`$
quand $`\|x\|\to\infty`$. Donc $`S=\{x\in\Omega:F_2(x)\le F_2(x_1)\}`$ ($`x_1\in\Omega`$) est borné ; il est fermé dans $`V`$, car
un point limite dans $`\bar\Omega\setminus\Omega`$ rendrait $`f`$ non bornée sur $`S`$ (lemme F.9(5)). La fonction continue $`F_2`$
atteint son minimum sur le compact $`S`$ en un point $`e`$ de l'ouvert $`\Omega`$, où $`0=\nabla F_2(e)=-e^*+e`$. Alors
$`H(e)=H(e^*)=H(e)^{-1}`$ d'après (1), donc $`H(e)^2=I`$, et $`H(e)=I`$ puisque $`H(e)>0`$.
(3) $`H(x)=H(ge)=g^{-T}g^{-1}`$ d'après le lemme F.9(4) et (2), et $`gg^T\in G`$ d'après le lemme F.8(4).
(4) $`e^{tZ}\in G`$ fixe $`e`$, donc $`I=H(e)=e^{-tZ^T}H(e)e^{-tZ}=e^{-tZ^T}e^{-tZ}`$ ; on dérive en $`t=0`$. $`\square`$

*Calcul :* $`e^*=e`$ et $`H(e)=I`$ pour les matrices définies positives de $`\mathrm{Herm}_3(\mathbb R)`$ — [tutoriel, section F.3](../../tutorial/fr/F_cones_and_order.md#s-F-3).

<a id="r-F-12"></a>

**Théorème F.12 (Koecher–Vinberg).** Soit $`\Omega`$ un cône convexe ouvert non vide d'un espace euclidien $`V`$. Alors
$`\Omega`$ est un cône symétrique si et seulement s'il existe un produit $`\circ`$ sur $`V`$ tel que $`(V,\circ,\langle\cdot,\cdot\rangle)`$
soit une algèbre de Jordan euclidienne et que $`\Omega`$ soit l'intérieur de son cône des carrés. Pour un cône symétrique,
on peut prendre pour unité n'importe quel $`e\in\Omega`$ tel que $`e^*=e`$, et le produit
```math
\langle u\circ v,w\rangle=-\tfrac12\,D^3\log\varphi(e)[u,v,w].
```

*Démonstration.* Le « si » est donné par les théorèmes F.5(2) et F.6(2). Réciproquement, soient $`\Omega`$ symétrique et
$`e`$ comme dans la proposition F.11(2) ; définissons $`L(u)`$ par $`\langle L(u)v,w\rangle=-\frac12D^3f(e)[u,v,w]`$ et
$`u\circ v=L(u)v`$.

*Étape 1 : commutativité, symétrie, unité.* $`D^3f(e)`$ est symétrique en ses trois arguments, donc $`u\circ v=v\circ u`$ et
$`L(u)^T=L(u)`$. D'après le lemme F.9(4) et $`H(e)=I`$, $`\langle e\circ v,w\rangle=\langle H(e)v,w\rangle=\langle v,w\rangle`$ :
$`e`$ est une unité.

*Étape 2 : $`L(u)\in\mathfrak g`$.* $`R(x)=H(x)^{-1}`$ appartient à $`G`$ (proposition F.11(3)), est lisse, et $`R(e)=I`$. Sa
dérivée est $`DR(e)[u]=-H(e)^{-1}DH(e)[u]H(e)^{-1}=-DH(e)[u]=2L(u)`$, puisque $`\langle DH(e)[u]v,w\rangle=D^3f(e)[u,v,w]`$.
Le lemme F.10(1) pour $`c(s)=R(e+su)`$ ($`s\ge0`$ petit, de sorte que $`e+su\in\Omega`$), avec le sous-groupe fermé $`G`$
(lemme F.8(4)), donne $`2L(u)\in\mathfrak g`$, donc $`L(u)\in\mathfrak g`$.

*Étape 3 : unicité.* Si $`Y\in\mathfrak g`$ est symétrique, alors $`Y=L(Ye)`$ : l'élément symétrique $`Y-L(Ye)\in\mathfrak g`$
annule $`e`$ ($`L(w)e=w\circ e=w`$), donc il est aussi antisymétrique d'après la proposition F.11(4), donc égal à $`0`$.

*Étape 4 : dérivations intérieures.* Soit $`D=[L(x),L(y)]`$. Il appartient à $`\mathfrak g`$ (lemme F.10(2)), $`D^T=-D`$, et
$`De=x\circ y-y\circ x=0`$. Pour $`z\in V`$, $`[D,L(z)]\in\mathfrak g`$ est symétrique et $`[D,L(z)]e=Dz-L(z)De=Dz`$, donc
$`[D,L(z)]=L(Dz)`$ d'après l'étape 3 ; appliqué à $`w`$ : $`D(z\circ w)=Dz\circ w+z\circ Dw`$.

*Étape 5 : l'identité de Jordan.* Pour $`D=[L(y),L(x)]`$, $`Dx=y\circ x^2-x\circ(x\circ y)`$ et
$`D(x^2)=y\circ(x\circ x^2)-x\circ(y\circ x^2)`$. L'étape 4 donne $`D(x^2)=2x\circ Dx`$, c'est-à-dire
$`y\circ(x\circ x^2)=3x\circ(x^2\circ y)-2x\circ(x\circ(x\circ y))`$ pour tout $`y`$ :
```math
L(x\circ x^2)=3L(x)L(x^2)-2L(x)^3 .
```
Tous les $`L`$ sont symétriques, donc en transposant on obtient $`L(x\circ x^2)=3L(x^2)L(x)-2L(x)^3`$, et
$`[L(x),L(x^2)]=0`$. Ainsi $`V`$ est une algèbre de Jordan euclidienne ; soit $`\Omega_J`$ l'intérieur de son cône des carrés.

*Étape 6 : $`\Omega=\Omega_J`$.* Pour $`b\in V`$, $`e^{L(b)}\in G`$ d'après l'étape 2, donc $`\exp b=e^{L(b)}e\in\Omega`$
(théorème F.6(1)). D'après le théorème F.6(2), $`\Omega_J\subseteq\Omega`$. Les deux cônes sont autoduaux pour
$`\langle\cdot,\cdot\rangle`$ (théorème F.5(2)), donc $`\Omega=\Omega^*\subseteq\Omega_J^*=\Omega_J`$. $`\square`$

*Calcul :* le produit $`-\frac12D^3f(e)`$, son unité et l'identité de Jordan, par différentiation automatique, pour les matrices $`3\times3`$ définies positives — [tutoriel, section F.3](../../tutorial/fr/F_cones_and_order.md#s-F-3).

La démonstration suit la structure de J. Faraut, A. Korányi, *Analysis on Symmetric Cones* (Oxford 1994), chapitres
I–III ; là où ils utilisent la structure de groupe de Lie de $`G(\Omega)`$, le lemme F.10 et la proposition F.11(3)
fournissent directement les groupes à un paramètre nécessaires. M. Koecher (1957) et E. B. Vinberg (1960) ont démontré
le théorème.

<a id="r-F-13"></a>

**Exemple F.13 (le cône de la programmation semi-définie).** Soient $`V=\mathrm{Herm}_n(\mathbb R)`$, $`\langle x,y\rangle=\mathrm{Tr}(xy)`$,
$`N=n(n+1)/2`$ et $`\Omega`$ l'ensemble des matrices définies positives, symétrique d'après les théorèmes F.5 et F.6. Alors
$`\varphi(x)=\varphi(1)\det(x)^{-(n+1)/2}`$, $`e=\kappa1`$ avec $`\kappa=\sqrt{(n+1)/2}`$, et le produit du théorème F.12 est
$`u\circ_\kappa v=\frac1{2\kappa}(uv+vu)`$, donc $`x\mapsto x/\kappa`$ est un isomorphisme sur l'algèbre de Jordan usuelle.

*Démonstration.* La décomposition spectrale $`x=\sum_i\lambda_ip_i`$ d'une matrice symétrique selon ses projecteurs propres
est la décomposition du théorème F.3 (unicité ; $`p_i\circ p_j=p_ip_j`$ puisque les $`p_i`$ commutent), donc, d'après le
théorème F.5(2), l'intérieur du cône des carrés est bien l'ensemble des matrices définies positives. Pour $`a\in\Omega`$,
$`P(a)x=axa`$ (proposition 13.2(1)) appartient à $`G(\Omega)`$ (théorème F.6(3)). Dans une base orthonormée de vecteurs
propres de $`a`$, les matrices $`E_{ii}`$ et $`(E_{ij}+E_{ji})/\sqrt2`$ ($`i<j`$) forment une base orthonormée de vecteurs propres
de $`P(a)`$, de valeurs propres $`a_ia_j`$, donc $`\det P(a)=\prod_{i\le j}a_ia_j=\det(a)^{n+1}`$ (chaque $`a_i`$ apparaît $`n+1`$
fois). Le lemme F.9(2) donne $`\varphi(a^2)=\varphi(P(a)1)=\det(a)^{-(n+1)}\varphi(1)`$, et tout $`x\in\Omega`$ est $`a^2`$ avec
$`a=\sqrt x`$. Avec
$`D\log\det(x)[u]=\mathrm{Tr}(x^{-1}u)`$ (d'après $`\det(x+su)=\det x\,\det(1+sx^{-1}u)=\det x\,(1+s\mathrm{Tr}(x^{-1}u)+O(s^2))`$)
et $`D(x^{-1})[u]=-x^{-1}ux^{-1}`$ :
$`x^*=\kappa^2x^{-1}`$ ; donc $`e^*=e`$ signifie $`e^2=\kappa^2 1`$ avec $`e`$ définie positive, c'est-à-dire $`e=\kappa1`$ ; $`D^2f(x)[v,w]=\kappa^2\mathrm{Tr}(x^{-1}vx^{-1}w)`$, et
$`D^3f(e)[u,v,w]=-\kappa^2\kappa^{-3}\mathrm{Tr}(uvw+vuw)=-\frac2\kappa\mathrm{Tr}((u\circ v)w)`$. $`\square`$

*Calcul :* $`e=\kappa1`$ et $`u\circ_\kappa v=\frac1{2\kappa}(uv+vu)`$ pour $`n=3`$ — [tutoriel, section F.3](../../tutorial/fr/F_cones_and_order.md#s-F-3).

<a id="s-F-4"></a>

## F.4 Isomorphismes d'ordre : le théorème de Kadison en dimension finie

Soient $`\mathcal M=\bigoplus_cM_{k_c}(\mathbb C)`$ et $`\mathcal N=\bigoplus_dM_{k'_d}(\mathbb C)`$, vues comme des matrices
diagonales par blocs ; $`x\ge0`$ signifie $`x=x^*`$ à valeurs propres positives ou nulles (de façon équivalente, chaque bloc
est semi-défini positif), et $`[0,1]=\{a\in\mathcal M_{sa}:0\le a\le1\}`$. Un **projecteur** est un $`p=p^*=p^2`$. Une
application linéaire $`\Phi`$ est **positive** si $`x\ge0`$ entraîne $`\Phi(x)\ge0`$ ; un **$`\ast`$-isomorphisme de Jordan** est
une bijection linéaire telle que $`\Phi(x^*)=\Phi(x)^*`$ et $`\Phi(x\circ y)=\Phi(x)\circ\Phi(y)`$ (section 13.5). Un point $`a`$
d'un ensemble convexe $`K`$ est **extrémal** s'il n'est pas le milieu $`\frac12(b+b')`$ de deux points distincts $`b,b'\in K`$
(de façon équivalente, s'il n'est pas de la forme $`tb+(1-t)b'`$ avec $`b\ne b'`$ et $`0<t<1`$ : un tel $`a`$ est le milieu de
$`a\pm\varepsilon(b-b')\in K`$ pour $`\varepsilon>0`$ petit).

<a id="r-F-14"></a>

**Lemme F.14.** Dans $`\mathcal M`$ :

1. Les points extrémaux de l'ensemble convexe $`[0,1]`$ sont exactement les projecteurs.
2. Pour des projecteurs $`p,q`$ : $`pq=0`$ si et seulement si $`p+q\le1`$, si et seulement si $`p\circ q=0`$.
3. Pour des projecteurs $`p,q`$ : $`p\le q`$ si et seulement si $`pq=p`$, si et seulement si $`p\circ q=p`$. Les projecteurs non
   nuls minimaux sont les projecteurs de rang un (dans un bloc).
4. Un projecteur $`p`$ est central si et seulement si $`pxp+(1-p)x(1-p)=x`$ pour tout $`x\in\mathcal M`$. Les projecteurs
   centraux sont les sommes des unités de bloc $`P_c`$, et les projecteurs centraux non nuls minimaux sont les $`P_c`$.

*Démonstration.* (1) Soient $`p`$ un projecteur et $`p=\frac12(a+b)`$ avec $`a,b\in[0,1]`$. Pour $`\xi\in\ker p`$,
$`0=\langle\xi,a\xi\rangle+\langle\xi,b\xi\rangle`$ avec deux termes $`\ge0`$, donc $`\|a^{1/2}\xi\|^2=0`$ et $`a\xi=0`$ ; pour $`\xi`$
dans l'image de $`p`$, le même argument avec $`1-a,1-b`$ donne $`a\xi=\xi`$. Donc $`a=p`$. Réciproquement, si $`a\in[0,1]`$ n'est
pas un projecteur, il a une valeur propre $`\lambda\in(0,1)`$ ; avec le projecteur spectral $`q_\lambda\in\mathcal M`$ (un
polynôme en $`a`$) et $`\varepsilon=\min(\lambda,1-\lambda)`$, $`a\pm\varepsilon q_\lambda\in[0,1]`$ et $`a`$ est leur milieu.
(2) Si $`pq=0`$, alors $`qp=(pq)^*=0`$, $`p+q`$ est un projecteur, et $`1-(p+q)=(1-p-q)^2\ge0`$ ; de plus $`p\circ q=0`$. Si
$`p+q\le1`$, alors pour tout $`\xi`$, $`\|qp\xi\|^2=\langle p\xi,qp\xi\rangle\le\langle p\xi,(1-p)p\xi\rangle=0`$, donc $`qp=0=pq`$. Si
$`p\circ q=0`$, alors $`pq=-qp`$, donc $`pq=p(pq)=-pqp=(qp)p=qp`$, et $`pq=-pq=0`$.
(3) $`pq=p`$ entraîne $`qp=p`$, $`p\circ q=p`$ et $`q-p=(q-p)^2\ge0`$. Si $`p\circ q=p`$, en multipliant $`pq+qp=2p`$ par $`p`$ à gauche
et à droite, on obtient $`pq+pqp=2p=pqp+qp`$, donc $`pq=qp=p`$. Si $`p\le q`$, alors $`p+(1-q)\le1`$ et (2) donne $`p(1-q)=0`$. Un
projecteur de rang $`\ge2`$ domine strictement le projecteur $`\xi\xi^*`$ sur une droite de son image contenue dans un bloc
($`\xi`$ un vecteur unitaire de l'image d'un bloc de $`p`$). Un projecteur de rang un $`p`$ est contenu dans un seul bloc, et si
$`0\ne q\le p`$, alors $`q=pq`$ a son image contenue dans celle de $`p`$, donc $`q=p`$.
(4) $`x=pxp+px(1-p)+(1-p)xp+(1-p)x(1-p)`$, et $`px(1-p)+(1-p)xp=0`$ si et seulement si les deux termes sont nuls (multiplier
par $`p`$ et $`1-p`$), si et seulement si $`px=pxp=xp`$. Un élément central commute avec toutes les unités matricielles de
chaque bloc, donc c'est un scalaire $`z_c`$ sur le bloc $`c`$, et $`z_c\in\{0,1\}`$ pour un projecteur. $`\square`$

<a id="r-F-15"></a>

**Théorème F.15 (Kadison, dimension finie).** Soit $`\Phi:\mathcal M\to\mathcal N`$ une bijection linéaire telle que
$`\Phi(1)=1`$ et que $`\Phi`$ et $`\Phi^{-1}`$ soient positives. Alors $`\Phi`$ est un $`*`$-isomorphisme de Jordan.
Réciproquement, tout $`*`$-isomorphisme de Jordan est unifère, et lui-même ainsi que son inverse sont positifs.

*Démonstration.* Tout $`x\in\mathcal M_{sa}`$ s'écrit $`x_+-x_-`$ avec $`x_\pm\ge0`$, donc $`\Phi(x)`$ est autoadjoint ; en
écrivant $`x=h+ik`$ avec $`h,k`$ autoadjoints, $`\Phi(x^*)=\Phi(h)-i\Phi(k)=\Phi(x)^*`$. Il en va de même pour $`\Phi^{-1}`$, donc
$`\Phi`$ est une bijection $`\mathbb R`$-linéaire $`\mathcal M_{sa}\to\mathcal N_{sa}`$ telle que $`x\ge0\iff\Phi(x)\ge0`$. Comme
$`\Phi(1)=1`$, $`\Phi`$ envoie $`[0,1]`$ sur $`[0,1]`$. Une bijection linéaire entre ensembles convexes envoie les points
extrémaux sur les points extrémaux (si $`\Phi(a)=\frac12(b+b')`$, alors $`a=\frac12(\Phi^{-1}b+\Phi^{-1}b')`$), donc, d'après le
lemme F.14(1), $`\Phi`$ envoie les projecteurs sur les projecteurs et, d'après le lemme F.14(2),
$`pq=0\iff p+q\le1\iff\Phi(p)+\Phi(q)\le1\iff\Phi(p)\Phi(q)=0`$. Pour $`x=\sum_i\lambda_ip_i\in\mathcal M_{sa}`$ (décomposition
spectrale, $`p_ip_j=0`$ pour $`i\ne j`$), $`\Phi(x)^2=\sum_i\lambda_i^2\Phi(p_i)=\Phi(x^2)`$. La polarisation,
$`x\circ y=\frac12((x+y)^2-x^2-y^2)`$, donne $`\Phi(x\circ y)=\Phi(x)\circ\Phi(y)`$ sur $`\mathcal M_{sa}`$, et la bilinéarité
complexe sur $`\mathcal M=\mathcal M_{sa}+i\mathcal M_{sa}`$.

Réciproquement, $`\Phi(1)\circ\Phi(y)=\Phi(y)`$ pour tout $`y`$ et $`\Phi`$ est surjectif, donc $`\Phi(1)`$ est une unité pour
$`\circ`$, et $`\Phi(1)=\Phi(1)\circ1=1`$. Tout $`x\ge0`$ est $`h^2`$ avec $`h=\sqrt x`$ autoadjoint (théorème 1.6), et
$`\Phi(h^2)=\Phi(h)^2\ge0`$ puisque $`\Phi(h)`$ est autoadjoint ; $`\Phi^{-1}`$ est aussi un $`*`$-isomorphisme de Jordan. $`\square`$

*Calcul :* un isomorphisme d'ordre unifère de $`M_2\oplus M_2\oplus M_3`$ qui est un $`*`$-isomorphisme de Jordan mais n'est pas multiplicatif — [tutoriel, section F.4](../../tutorial/fr/F_cones_and_order.md#s-F-4).

Kadison a démontré le théorème pour toutes les C\*-algèbres unifères (R. V. Kadison, *A generalized Schwarz inequality
and algebraic invariants for operator algebras*, Ann. of Math. 56 (1952) 494–503) ; cette version n'est pas utilisée
dans le cours.

<a id="s-F-5"></a>

## F.5 $`*`$-isomorphismes de Jordan de $`\bigoplus_cM_{k_c}`$

<a id="r-F-16"></a>

**Lemme F.16.** Soit $`\Phi:\mathcal M\to\mathcal N`$ un $`*`$-isomorphisme de Jordan. Alors $`\Phi(1)=1`$,
$`\Phi(xyx)=\Phi(x)\Phi(y)\Phi(x)`$, et $`\Phi`$ envoie bijectivement les projecteurs de $`\mathcal M`$ sur ceux de $`\mathcal N`$,
en préservant et en reflétant $`pq=0`$, $`p\le q`$, la minimalité, la centralité et la minimalité parmi les projecteurs
centraux.

*Démonstration.* $`\Phi(1)=1`$ comme dans le théorème F.15. Dans toute algèbre associative,
$`xyx=2x\circ(x\circ y)-(x\circ x)\circ y`$, ce qui donne la deuxième assertion. $`p=p^*=p\circ p`$ est préservé, ainsi que les
caractérisations du lemme F.14(2,3,4), qui n'utilisent que $`\circ`$, le produit triple $`xyx`$, $`1`$ et des quantificateurs
portant sur toute l'algèbre ; $`\Phi^{-1}`$ est aussi un $`*`$-isomorphisme de Jordan. $`\square`$

<a id="r-F-17"></a>

**Théorème F.17 ($`*`$-automorphismes de Jordan de $`M_n`$).** Tout $`*`$-isomorphisme de Jordan
$`\Phi:M_n(\mathbb C)\to M_n(\mathbb C)`$ est $`x\mapsto UxU^*`$ ou $`x\mapsto Ux^TU^*`$ pour un unitaire $`U`$. Pour $`n\ge2`$, les
premières applications sont multiplicatives et les secondes antimultiplicatives, donc aucune application n'est des deux
formes.

*Démonstration.* Pour $`n=1`$, $`\Phi=\mathrm{id}`$. Soient $`n\ge2`$ et $`E_{ij}`$ les unités matricielles. D'après le lemme
F.16, les $`\Phi(E_{ii})`$ sont des projecteurs de rang un deux à deux orthogonaux, $`\Phi(E_{ii})=\xi_i\xi_i^*`$ avec des
$`\xi_i`$ orthonormés ; avec l'unitaire $`W`$, $`We_i=\xi_i`$, le $`*`$-isomorphisme de Jordan $`\Phi_1(x)=W^*\Phi(x)W`$ vérifie
$`\Phi_1(E_{ii})=E_{ii}`$.

Soit $`i\ne j`$. L'espace $`S_{ij}=\mathrm{span}\{E_{ij},E_{ji}\}`$ est $`\{x:pxp=x,\ E_{ii}xE_{ii}=E_{jj}xE_{jj}=0\}`$ avec
$`p=E_{ii}+E_{jj}`$, donc $`\Phi_1(S_{ij})\subseteq S_{ij}`$ d'après le lemme F.16. Écrivons $`\Phi_1(E_{ij})=\alpha E_{ij}+\beta E_{ji}`$.
De $`E_{ij}\circ E_{ij}=0`$ : $`0=\Phi_1(E_{ij})^2=\alpha\beta(E_{ii}+E_{jj})`$, donc exactement l'un des $`\alpha,\beta`$ est non nul
($`\Phi_1`$ est injectif). Soit $`\varepsilon(i,j)=+`$ si $`\Phi_1(E_{ij})\in\mathbb CE_{ij}`$ et $`-`$ si $`\Phi_1(E_{ij})\in\mathbb CE_{ji}`$ ;
$`\varepsilon(j,i)=\varepsilon(i,j)`$ puisque $`\Phi_1(E_{ji})=\Phi_1(E_{ij})^*`$. Pour $`i,j,k`$ distincts,
$`E_{ij}\circ E_{jk}=\frac12E_{ik}\ne0`$ ; si $`\varepsilon(i,j)=+`$ et $`\varepsilon(j,k)=-`$, alors
$`\Phi_1(E_{ij})\circ\Phi_1(E_{jk})\in\mathbb C\,E_{ij}\circ E_{kj}=\{0\}`$, ce qui est contradictoire ; de même pour $`(-,+)`$.
Ainsi $`\varepsilon(i,j)=\varepsilon(j,k)`$ dès que $`i,j,k`$ sont distincts, et $`\varepsilon`$ est constant : deux paires ou bien
ont un indice commun, ou bien sont reliées par une troisième paire. Si $`\varepsilon=-`$, on remplace $`\Phi_1`$ par
$`\Phi_1\circ T`$ ($`T`$ la transposée, un $`*`$-automorphisme de Jordan d'après la proposition 13.9), qui a $`\varepsilon=+`$ et
fixe toujours les $`E_{ii}`$.

Maintenant $`\Phi_1(E_{ij})=\alpha_{ij}E_{ij}`$ pour tous $`i,j`$, avec $`\alpha_{ii}=1`$ et $`\alpha_{ji}=\bar\alpha_{ij}`$. De
$`E_{ij}\circ E_{ji}=\frac12(E_{ii}+E_{jj})`$ : $`|\alpha_{ij}|^2=1`$ ; de $`E_{ij}\circ E_{jk}=\frac12E_{ik}`$ ($`i,j,k`$ distincts) :
$`\alpha_{ij}\alpha_{jk}=\alpha_{ik}`$. Avec $`u_i=\alpha_{i1}`$, cela donne $`\alpha_{ij}=\alpha_{i1}\alpha_{1j}=u_i\bar u_j`$ (aussi
pour $`i=j`$, $`i=1`$ ou $`j=1`$). Donc $`\Phi_1(x)=DxD^*`$ avec $`D=\mathrm{diag}(u_i)`$, et $`\Phi(x)=WDxD^*W^*`$ ou
$`\Phi(x)=WDx^TD^*W^*`$. Si les deux formes coïncidaient, $`T`$ serait multiplicative, mais $`(E_{12}E_{21})^T=E_{11}\ne E_{22}=E_{12}^TE_{21}^T`$.
$`\square`$

*Calcul :* retrouver $`U`$ à partir de $`x\mapsto Ux^TU^*`$ en suivant les étapes de la démonstration — [tutoriel, section F.5](../../tutorial/fr/F_cones_and_order.md#s-F-5).

Cela démontre le résultat de N. Jacobson, C. E. Rickart (Trans. Amer. Math. Soc. 69 (1950) 479–502) et I. N.
Herstein (Trans. Amer. Math. Soc. 81 (1956) 331–341) cité à la section 13.5, dans le cas de la dimension finie, pour tout
$`n`$ ; la proposition 13.11 est le cas $`n=2`$. Un automorphisme de Jordan de $`\mathrm{Herm}_n(\mathbb C)`$ s'étend en
un $`*`$-automorphisme de Jordan de $`M_n`$ par $`\Phi(x+iy)=\Phi(x)+i\Phi(y)`$, donc il est lui aussi de l'une des deux formes.

<a id="r-F-18"></a>

**Théorème F.18 (sommes directes).** Soit $`\Phi:\bigoplus_{c=1}^CM_{k_c}\to\bigoplus_{d=1}^{C'}M_{k'_d}`$ un
$`*`$-isomorphisme de Jordan. Alors $`C'=C`$, et il existe une bijection $`\sigma`$ telle que $`k'_{\sigma(c)}=k_c`$, des unitaires
$`U_c\in M_{k_c}`$ et un ensemble $`S`$ de blocs tels que
```math
\Phi(x)_{\sigma(c)}=U_cx_cU_c^*\ (c\notin S),\qquad \Phi(x)_{\sigma(c)}=U_cx_c^TU_c^*\ (c\in S).
```
Réciproquement, toute application de ce type est un $`*`$-isomorphisme de Jordan, et $`S`$ est unique aux blocs tels que
$`k_c=1`$ près.

*Démonstration.* D'après le lemme F.16 et le lemme F.14(4), $`\Phi(P_c)=P'_{\sigma(c)}`$ pour une bijection $`\sigma`$ entre les
projecteurs centraux minimaux, donc $`C=C'`$. Si $`x=P_cxP_c`$, alors $`\Phi(x)=P'_{\sigma(c)}\Phi(x)P'_{\sigma(c)}`$, donc $`\Phi`$
envoie le bloc $`c`$ dans le bloc $`\sigma(c)`$, et $`\Phi^{-1}`$ envoie le bloc $`\sigma(c)`$ dans le bloc $`c`$ ; par conséquent
$`k_c^2=k'^2_{\sigma(c)}`$, et la restriction est un $`*`$-isomorphisme de Jordan $`M_{k_c}\to M_{k_c}`$, de l'une des deux formes
d'après le théorème F.17. Comme $`x=\sum_cP_cxP_c`$, cela détermine $`\Phi`$. La réciproque se vérifie bloc par bloc
(proposition 13.9), et l'unicité de $`S`$ découle du théorème F.17. $`\square`$

*Calcul :* une application de cette forme qui échange deux blocs et en transpose deux — [tutoriel, section F.4](../../tutorial/fr/F_cones_and_order.md#s-F-4).

<a id="r-F-19"></a>

**Corollaire F.19 (isomorphismes d'ordre unifères).** Les bijections linéaires unifères $`\Phi:\mathcal M\to\mathcal N`$
telles que $`\Phi`$ et $`\Phi^{-1}`$ soient positives sont exactement les applications du théorème F.18. Une telle $`\Phi`$ est
$`2`$-positive si et seulement si elle est complètement positive, si et seulement si $`S`$ ne contient aucun bloc tel que
$`k_c\ge2`$, si et seulement si $`\Phi`$ est un $`*`$-isomorphisme.

*Démonstration.* La première phrase résulte des théorèmes F.15 et F.18. Si $`S`$ ne contient aucun bloc tel que
$`k_c\ge2`$, chaque application de bloc est $`x\mapsto U_cxU_c^*`$ (pour $`k_c=1`$, $`x^T=x`$), donc $`\Phi`$ est complètement
positive (proposition 6.4) et multiplicative. Si $`c\in S`$ avec $`k_c\ge2`$ et si $`\Phi`$ est $`2`$-positive, alors il en va
de même de l'application de bloc $`\Phi_{\sigma(c)c}(a)=U_ca^TU_c^*`$ (proposition 6.4), ainsi que de
$`a\mapsto U_c^*\Phi_{\sigma(c)c}(a)U_c=a^T`$, puisque $`[A_{ij}]\ge0`$ entraîne $`[U^*A_{ij}U]=(1\otimes U)^*[A_{ij}](1\otimes U)\ge0`$ ;
cela contredit la proposition 6.1. Une telle $`\Phi`$ n'est pas non plus multiplicative, d'après le théorème F.17. $`\square`$

Ainsi un isomorphisme d'ordre unifère est déterminé, aux unitaires près, par les blocs égaux qu'il échange et par les
blocs qu'il transpose, et, d'après le corollaire F.19, ce sont exactement les blocs transposés qui empêchent la complète
positivité.

<a id="x-F"></a>

## Exercices

1. Montrer que $`\mathbb R^N`$, muni du produit coordonnée par coordonnée et du produit scalaire usuel, est une algèbre de
   Jordan euclidienne dont le cône des carrés est $`\mathbb R^N_{\ge0}`$, et déduire du corollaire F.19 que les
   isomorphismes d'ordre unifères de $`\mathbb C^N=\bigoplus_{c=1}^NM_1`$ sont les permutations des coordonnées.
2. (**La positivité ne suffit pas.**) Pour $`0<t<1`$, soit $`\Psi_t(x)=(1-t)x+\frac tn\mathrm{Tr}(x)1`$ sur $`M_n`$,
   $`n\ge2`$. Montrer que $`\Psi_t`$ est une bijection positive unifère, que $`\Psi_t^{-1}`$ n'est pas positive et que $`\Psi_t`$
   n'envoie pas les projecteurs sur des projecteurs.

<details>
<summary>Solutions</summary>

1. Le produit est commutatif et associatif, donc de Jordan, d'unité $`(1,\dots,1)`$, et $`L(x)`$ est la matrice diagonale
   $`\mathrm{diag}(x)`$, qui est symétrique. Les carrés sont les vecteurs à coordonnées positives ou nulles. D'après
   le corollaire F.19 avec tous les $`k_c=1`$, $`\Phi(x)_{\sigma(c)}=U_cx_cU_c^*=x_c`$ (les unitaires de $`M_1`$ sont des
   phases), donc $`\Phi`$ permute les coordonnées ; réciproquement, les permutations sont des isomorphismes d'ordre
   unifères.

2. $`\Psi_t(1)=(1-t)1+t1=1`$ ; pour $`x\ge0`$, les deux termes sont $`\ge0`$. Comme $`\mathrm{Tr}\,\Psi_t(x)=\mathrm{Tr}\,x`$,
   $`\Psi_t^{-1}(y)=\frac1{1-t}\big(y-\frac tn\mathrm{Tr}(y)1\big)`$. Pour un projecteur de rang un $`p`$, $`\Psi_t^{-1}(p)`$
   a la valeur propre $`-\frac{t}{n(1-t)}<0`$ sur $`\ker p`$, et $`\Psi_t(p)`$ a la valeur propre $`\frac tn\notin\{0,1\}`$. C'est
   cohérent avec le théorème F.15 : c'est la positivité de l'inverse qui force les projecteurs à être envoyés sur des
   projecteurs.

</details>

*Calcul :* l'exercice 2, numériquement — [tutoriel, Exercices](../../tutorial/fr/F_cones_and_order.md#x-F).

**Retour :** [chapitre 13](13_jordan_algebras.md), sections 13.4–13.6 : le théorème F.12 démontre le théorème 13.8, le
théorème F.15 démontre le théorème 13.10 pour les algèbres $`\bigoplus_cM_{k_c}(\mathbb C)`$, le théorème F.17 étend la
proposition 13.11 à tout $`n`$, et le corollaire F.7 étend la proposition 13.14 à toute algèbre de Jordan euclidienne.
