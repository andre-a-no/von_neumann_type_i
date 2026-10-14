[English](../B_schur_weyl.md) · [Español](../es/B_schur_weyl.md) · **Français** · [Deutsch](../de/B_schur_weyl.md) · [中文](../zh/B_schur_weyl.md) · [日本語](../ja/B_schur_weyl.md) · [Русский](../ru/B_schur_weyl.md)

# Annexe B. Dualité de Schur–Weyl et calcul de Weingarten

Le groupe unitaire $U(d)$ agit sur $(\mathbb C^d)^{\otimes n}$ par $V\mapsto V^{\otimes n}$, et le groupe symétrique $S_n$
en permutant les facteurs tensoriels. Chacune des deux familles d'opérateurs est le commutant de l'autre
(**dualité de Schur–Weyl**). Cette annexe démontre complètement la dualité et en déduit la formule de Weingarten pour
les moments de tout ordre d'un unitaire de Haar. Chemin faisant, elle démontre les faits que les chapitres 4, 5, 9 et
11 empruntent ici : le théorème 4.14(b) pour tout $L$, l'irréductibilité des représentations de $S_L$ sur les espaces
de multiplicité du spin total, les nombres de Catalan de la section 4.5, le théorème 5.12, le théorème de
Wigner–Eckart pour $SU(2)$ utilisé à la section 9.3, et la mesure de Haar du groupe orthogonal utilisée au chapitre 11.

La démarche est élémentaire : une identité de polarisation identifie le commutant des permutations avec l'espace
engendré par tous les $A^{\otimes n}$, $A\in M_d(\mathbb C)$ ; un argument de transformée de Cayley montre que les
unitaires $V^{\otimes n}$ engendrent le même espace ; le théorème du bicommutant (théorème 2.2) transforme cela en la
dualité. La formule de Weingarten s'ensuit alors, car la moyenne de Haar est la projection orthogonale sur le
commutant.

Les vérifications numériques de cette annexe se trouvent dans [Tutoriel, annexe B](../../tutorial/fr/B_schur_weyl.md).

Dans toute l'annexe, $d,n\ge1$, $\mathcal H_n=(\mathbb C^d)^{\otimes n}$ avec la base orthonormée
$e_{\mathbf i}=e_{i_1}\otimes\dots\otimes e_{i_n}$, où un multi-indice $\mathbf i=(i_1,\dots,i_n)\in[d]^n$ est lu comme une
application $[n]\to[d]$, $m\mapsto i_m$ ($[d]=\{1,\dots,d\}$). $\mathrm{End}\,\mathcal H_n=M_{d^n}(\mathbb C)$ est muni du
produit scalaire de Hilbert–Schmidt $\langle X,Y\rangle=\mathrm{Tr}(X^*Y)$. Pour $X\in M_d$ et $1\le m\le n$,
$X^{(m)}=1\otimes\dots\otimes X\otimes\dots\otimes1$ ($X$ en position $m$) et $d\Pi(X)=\sum_{m=1}^nX^{(m)}$. Nous notons
$$
\mathcal T_n=\mathrm{span}\{A^{\otimes n}:A\in M_d(\mathbb C)\},\qquad
\mathcal U_n=\mathrm{span}\{V^{\otimes n}:V\in U(d)\}.
$$
Les produits de produits de Kronecker se calculent facteur par facteur, $(A_1\otimes\dots\otimes A_n)(B_1\otimes\dots\otimes B_n)=A_1B_1\otimes\dots\otimes A_nB_n$,
et $(A_1\otimes\dots\otimes A_n)^*=A_1^*\otimes\dots\otimes A_n^*$. Les permutations se composent comme des applications,
$(\sigma\tau)(m)=\sigma(\tau(m))$, et $\#\text{cycles}(\sigma)$ compte les cycles de $\sigma$, points fixes compris.

<a id="s-B-1"></a>

## B.1 Permutations des facteurs tensoriels

Pour $\sigma\in S_n$, soit $\pi(\sigma)e_{\mathbf i}=e_{\mathbf i\circ\sigma^{-1}}$ : le vecteur de base de la case $m$ passe dans
la case $\sigma(m)$. C'est la convention de `decompose.permutation_representation` et des chapitres 2 et 4.

*Dans la bibliothèque :* les opérateurs $\pi(\sigma)$ sont construits par `decompose.permutation_representation` — [tutoriel, section B.3](../../tutorial/fr/B_schur_weyl.md#s-B-3).

<a id="r-B-1"></a>

**Lemme B.1.** (a) $\pi$ est une représentation unitaire de $S_n$, et
$\pi(\sigma)(x_1\otimes\dots\otimes x_n)=x_{\sigma^{-1}(1)}\otimes\dots\otimes x_{\sigma^{-1}(n)}$ pour tous $x_m\in\mathbb C^d$.
(b) $\pi(\sigma)(A_1\otimes\dots\otimes A_n)\pi(\sigma)^*=A_{\sigma^{-1}(1)}\otimes\dots\otimes A_{\sigma^{-1}(n)}$ pour tous
$A_m\in M_d$. En particulier $\pi(\sigma)$ commute avec tout $A^{\otimes n}$ et tout $d\Pi(X)$, et
$\pi(\sigma)X^{(m)}\pi(\sigma)^*=X^{(\sigma(m))}$.
(c) $\mathrm{Tr}\,\pi(\sigma)=d^{\#\text{cycles}(\sigma)}$ et $\langle\pi(\sigma),\pi(\tau)\rangle=d^{\#\text{cycles}(\sigma^{-1}\tau)}$.
Le nombre de cycles est invariant par conjugaison et par inversion, et deux permutations sont conjuguées si et
seulement si elles ont le même type cyclique.
(d) Les $n!$ opérateurs $\pi(\sigma)$ sont linéairement indépendants si et seulement si $n\le d$.

*Démonstration.* (a) $\mathbf i\mapsto\mathbf i\circ\sigma^{-1}$ est une bijection de $[d]^n$, donc $\pi(\sigma)$ permute une base
orthonormée et est unitaire. $\pi(\sigma)\pi(\tau)e_{\mathbf i}=\pi(\sigma)e_{\mathbf i\circ\tau^{-1}}=e_{\mathbf i\circ\tau^{-1}\circ\sigma^{-1}}=\pi(\sigma\tau)e_{\mathbf i}$.
En développant $x_m=\sum_ax_m(a)e_a$, $x_1\otimes\dots\otimes x_n=\sum_{\mathbf i}\prod_mx_m(i_m)\,e_{\mathbf i}$ ; en appliquant
$\pi(\sigma)$ et en substituant $\mathbf j=\mathbf i\circ\sigma^{-1}$, on obtient $\sum_{\mathbf j}\prod_mx_m(j_{\sigma(m)})\,e_{\mathbf j}=\sum_{\mathbf j}\prod_{p}x_{\sigma^{-1}(p)}(j_p)\,e_{\mathbf j}$,
qui est le vecteur produit annoncé.
(b) D'après (a), sur un vecteur produit : $\pi(\sigma)^*=\pi(\sigma^{-1})$ place $x_{\sigma(m)}$ dans la case $m$, l'opérateur
y applique $A_m$, et $\pi(\sigma)$ envoie la case $m$ dans la case $\sigma(m)$ ; la case $p$ contient donc finalement
$A_{\sigma^{-1}(p)}x_p$. Les vecteurs produits engendrent $\mathcal H_n$. Pour $X^{(m)}$, le facteur $X$ passe de la case
$m$ à la case $\sigma(m)$, donc $d\Pi(X)$ est fixé.
(c) $\pi(\sigma)$ est une matrice de permutation, et $e_{\mathbf i}$ est fixé si et seulement si $i_{\sigma^{-1}(p)}=i_p$
pour tout $p$, c'est-à-dire si et seulement si $\mathbf i$ est constant sur les cycles de $\sigma$ : il y a
$d^{\#\text{cycles}(\sigma)}$ tels $\mathbf i$. Ensuite
$\langle\pi(\sigma),\pi(\tau)\rangle=\mathrm{Tr}(\pi(\sigma^{-1})\pi(\tau))=\mathrm{Tr}\,\pi(\sigma^{-1}\tau)$. Les cycles de
$\kappa\rho\kappa^{-1}$ sont les images par $\kappa$ des cycles de $\rho$, et $\rho^{-1}$ parcourt les mêmes cycles à
rebours. Réciproquement, si $\rho,\rho'$ ont le même type cyclique, on écrit leurs cycles l'un au-dessous de l'autre,
ordonnés par longueur, et l'on note $\kappa$ l'application qui envoie chaque élément de la première ligne sur
l'élément situé au-dessous ; alors $\kappa\rho\kappa^{-1}=\rho'$.
(d) Si $n\le d$, soit $\mathbf u=(1,2,\dots,n)$. Les vecteurs $\pi(\sigma)e_{\mathbf u}=e_{\mathbf u\circ\sigma^{-1}}$ sont des
vecteurs de base distincts pour des $\sigma$ distinctes, donc $\sum_\sigma c_\sigma\pi(\sigma)=0$ appliqué à $e_{\mathbf u}$
donne $c_\sigma=0$ pour tout $\sigma$. Si $n>d$, soit $\mathcal A=\sum_\sigma\mathrm{sgn}(\sigma)\pi(\sigma)$. Tout
$\mathbf i\in[d]^n$ vérifie $i_a=i_b$ pour certains $a\neq b$ ; pour la transposition $t=(a\,b)$, $\pi(t)e_{\mathbf i}=e_{\mathbf i}$
et $\mathcal A\pi(t)=\sum_\sigma\mathrm{sgn}(\sigma)\pi(\sigma t)=-\mathcal A$. Donc
$\mathcal Ae_{\mathbf i}=\mathcal A\pi(t)e_{\mathbf i}=-\mathcal Ae_{\mathbf i}=0$ pour tout $\mathbf i$ : $\mathcal A=0$ est une
relation linéaire non triviale. $\square$

*Calcul :* le lemme B.1(a) et (d) — [tutoriel, section B.3](../../tutorial/fr/B_schur_weyl.md#s-B-3).

<a id="s-B-2"></a>

## B.2 Polarisation : le commutant des permutations

<a id="r-B-2"></a>

**Lemme B.2 (polarisation).** Pour $A_1,\dots,A_n\in M_d(\mathbb C)$, avec $|\varepsilon|=\sum_m\varepsilon_m$,
$$
\sum_{\sigma\in S_n}A_{\sigma(1)}\otimes\dots\otimes A_{\sigma(n)}=
\sum_{\varepsilon\in\{0,1\}^n}(-1)^{n-|\varepsilon|}\Big(\sum_{m=1}^n\varepsilon_mA_m\Big)^{\otimes n}.
$$

*Démonstration.* Par multilinéarité, $\big(\sum_m\varepsilon_mA_m\big)^{\otimes n}=\sum_f\varepsilon_{f(1)}\cdots\varepsilon_{f(n)}\,A_{f(1)}\otimes\dots\otimes A_{f(n)}$,
la somme portant sur toutes les applications $f:[n]\to[n]$. Comme $\varepsilon_m^r=\varepsilon_m$ pour $r\ge1$, le
coefficient du terme $f$ dans le membre de droite est
$$
\sum_{\varepsilon}\prod_{m=1}^n(-1)^{1-\varepsilon_m}\prod_{m\in f([n])}\varepsilon_m
=\prod_{m\notin f([n])}\Big(\sum_{\varepsilon_m=0,1}(-1)^{1-\varepsilon_m}\Big)\prod_{m\in f([n])}\Big(\sum_{\varepsilon_m=0,1}(-1)^{1-\varepsilon_m}\varepsilon_m\Big).
$$
Le premier type de facteur vaut $-1+1=0$, le second vaut $1$. Le coefficient vaut donc $1$ si $f$ est surjective,
c'est-à-dire une permutation, et $0$ sinon. $\square$

*Calcul :* le lemme B.2 pour $n=3$ — [tutoriel, section B.3](../../tutorial/fr/B_schur_weyl.md#s-B-3).

Le **symétriseur** est $\mathcal S(X)=\frac1{n!}\sum_\sigma\pi(\sigma)X\pi(\sigma)^*$ sur $\mathrm{End}\,\mathcal H_n$.

<a id="r-B-3"></a>

**Proposition B.3.** $\mathcal S$ est la projection orthogonale de $\mathrm{End}\,\mathcal H_n$ sur $\pi(S_n)'$, et
$$
\pi(S_n)'=\mathcal T_n,\qquad\dim\pi(S_n)'=\binom{d^2+n-1}n .
$$

*Démonstration.* $\mathcal S$ est le twirl du groupe fini $S_n$ dans la représentation $\pi$, donc la projection
orthogonale sur $\pi(S_n)'$ (théorème 12.5 avec $\pi_1=\pi_2=\pi$ ; voir aussi la section 4.2). D'après le lemme B.1(b),
$\mathcal T_n\subseteq\pi(S_n)'$. Réciproquement, $\mathrm{End}\,\mathcal H_n$ est engendré par les produits
$A_1\otimes\dots\otimes A_n$ d'unités matricielles, et d'après le lemme B.1(b) et la substitution $\sigma\to\sigma^{-1}$,
$\mathcal S(A_1\otimes\dots\otimes A_n)=\frac1{n!}\sum_\sigma A_{\sigma(1)}\otimes\dots\otimes A_{\sigma(n)}$, qui appartient à
$\mathcal T_n$ d'après le lemme B.2. Donc $\pi(S_n)'=\mathcal S(\mathrm{End}\,\mathcal H_n)\subseteq\mathcal T_n$.

*Dimension.* Soient $P=[d]\times[d]$ et $E_\alpha=E_{\alpha_1}\otimes\dots\otimes E_{\alpha_n}$ pour $\alpha\in P^n$, avec les
unités matricielles $E_{(a,b)}=|e_a\rangle\langle e_b|$ : une base de $\mathrm{End}\,\mathcal H_n$. D'après le lemme B.1(b),
$\mathcal S(E_\alpha)=\frac1{n!}\sum_\sigma E_{\alpha\circ\sigma^{-1}}$, qui vaut $|O_\alpha|^{-1}\sum_{\beta\in O_\alpha}E_\beta$ pour l'orbite
$O_\alpha=\{\alpha\circ\sigma\}$ (chaque élément de l'orbite apparaît pour le même nombre $n!/|O_\alpha|$ de $\sigma$). Des
sommes sur des orbites différentes font intervenir des ensembles disjoints d'éléments de base, donc l'image de
$\mathcal S$ a une base indexée par les orbites. Deux uplets sont dans la même orbite si et seulement si chaque élément
de $P$ y apparaît le même nombre de fois : les orbites sont les multiensembles de taille $n$ formés à partir des
$N=d^2$ éléments de $P$. En numérotant $P$ par $1,\dots,N$, un tel multiensemble est une suite
$1\le c_1\le\dots\le c_n\le N$, et $c_m\mapsto c_m+m-1$ est une bijection sur les suites strictement croissantes dans
$[N+n-1]$, c'est-à-dire les parties à $n$ éléments ; il y en a $\binom{N+n-1}n$. $\square$

*Calcul :* $\dim\pi(S_n)'=\binom{d^2+n-1}n$ pour de petits $d$ et $n$ — [tutoriel, section B.3](../../tutorial/fr/B_schur_weyl.md#s-B-3).

<a id="s-B-3"></a>

## B.3 Dualité de Schur–Weyl

<a id="r-B-4"></a>

**Lemme B.4 (polynômes qui s'annulent sur le groupe unitaire).** Soit $p:M_d(\mathbb C)\to\mathbb C$ un polynôme en les
$d^2$ coefficients, homogène de degré $n\ge1$. Si $p(V)=0$ pour tout $V\in U(d)$, alors $p=0$.

*Démonstration.* *Étape 1 (transformée de Cayley).* Pour $A\in M_d$, soit $q(A)=p\big((1+iA)\mathrm{adj}(1-iA)\big)$,
où $\mathrm{adj}$ désigne la comatrice transposée ; $q$ est un polynôme en les coefficients de $A$. Soit $H$
hermitienne. Ses valeurs propres sont réelles (corollaire 1.4), donc $1-iH$ est inversible,
$\mathrm{adj}(1-iH)=\det(1-iH)(1-iH)^{-1}$ (règle de Cramer), et par homogénéité $q(H)=\det(1-iH)^n\,p(C)$ avec
$C=(1+iH)(1-iH)^{-1}$. Les matrices $1\pm iH$ et leurs inverses sont des fonctions de $H$ et commutent (théorème 1.5),
et $(1+iH)^*=1-iH$, donc $C^*C=(1+iH)^{-1}(1-iH)(1+iH)(1-iH)^{-1}=1$. Par conséquent $q(H)=0$.

*Étape 2 (des matrices hermitiennes à toutes les matrices).* Choisissons une base $B_1,\dots,B_{d^2}$ de l'espace
vectoriel réel des matrices hermitiennes. C'est aussi une base de $M_d(\mathbb C)$ sur $\mathbb C$, puisque toute $A$
s'écrit $\frac{A+A^*}2+i\frac{A-A^*}{2i}$ et que les dimensions coïncident. $r(t)=q(\sum_kt_kB_k)$ est un polynôme en
$t\in\mathbb C^{d^2}$ qui s'annule sur $\mathbb R^{d^2}$. Un tel polynôme a tous ses coefficients nuls : en une variable,
il a une infinité de racines ; pour $N$ variables, on écrit $r=\sum_jr_j(t')t_N^j$ — pour tout $t'$ réel, le polynôme
en $t_N$ s'annule sur $\mathbb R$, donc chaque $r_j$ s'annule sur $\mathbb R^{N-1}$, et l'on conclut par récurrence.
Donc $q(A)=r(t)=0$ pour toute $A=\sum_kt_kB_k\in M_d$.

*Étape 3 (inversion de la transformée).* Soit $B\in M_d$ avec $1+B$ inversible, et $A=-i(1+B)^{-1}(B-1)$. Alors
$1-iA=1-(1+B)^{-1}(B-1)=2(1+B)^{-1}$ et $1+iA=(1+B)^{-1}2B$ ; comme $B$ commute avec $1+B$,
$(1+iA)(1-iA)^{-1}=(1+B)^{-1}B(1+B)=B$. Donc $0=q(A)=\det(1-iA)^np(B)$ avec $\det(1-iA)\neq0$ : $p(B)=0$.

*Étape 4.* Pour tout $B\neq0$ et $0<t<1/\|B\|$, $-1$ n'est pas valeur propre de $tB$ (si $tBv=-v$, alors
$\|v\|=t\|Bv\|\le t\|B\|\|v\|<\|v\|$ sauf si $v=0$), donc $t^np(B)=p(tB)=0$. Enfin $p(0)=0$ par homogénéité. $\square$

<a id="r-B-5"></a>

**Lemme B.5 (identités de Newton).** Soient $x_1,\dots,x_n$ des éléments d'une algèbre qui commutent deux à deux,
$x_S=\prod_{m\in S}x_m$, $e_k=\sum_{|S|=k}x_S$ ($e_0=1$) et $p_j=\sum_mx_m^j$. Alors
$k\,e_k=\sum_{j=1}^k(-1)^{j-1}e_{k-j}\,p_j$ pour $1\le k\le n$.

*Démonstration.* Pour $1\le j\le k$, soit $a_j=\sum x_Sx_m^j$, la somme portant sur les couples $(S,m)$ avec
$|S|=k-j$ et $m\notin S$, et soit $a_{k+1}=0$. Dans $e_{k-j}p_j=\sum_{|S|=k-j}\sum_mx_Sx_m^j$, les termes avec
$m\notin S$ donnent $a_j$, et les termes avec $m\in S$, écrits avec $S'=S\setminus\{m\}$, valent $x_{S'}x_m^{j+1}$ et
donnent $a_{j+1}$. La somme alternée est donc télescopique et vaut $a_1$, et $a_1=\sum_{|T|=k}|T|\,x_T=k\,e_k$ (chaque
partie $T$ à $k$ éléments provient des $k$ choix de $m\in T$). $\square$

*Calcul :* le lemme B.5 pour $x_m=A^{(m)}$ et $n=4$, qui retrouve $A^{\otimes4}$ à partir des $d\Pi(A^j)$ — [tutoriel, section B.3](../../tutorial/fr/B_schur_weyl.md#s-B-3).

<a id="r-B-6"></a>

**Théorème B.6 (dualité de Schur–Weyl).** Soient $d,n\ge1$.
(a) $\pi(S_n)'=\mathcal T_n=\mathcal U_n=\mathrm{span}\{V^{\otimes n}:V\in SU(d)\}$. Cet espace est une $*$-algèbre, c'est
donc l'algèbre engendrée par les $V^{\otimes n}$, $V\in U(d)$.
(b) Pour tout ensemble $\mathcal X\subseteq M_d$ qui engendre $M_d$ comme espace vectoriel, l'algèbre unifère engendrée par
$\{d\Pi(X):X\in\mathcal X\}$ est $\mathcal T_n$.
(c) $\{V^{\otimes n}:V\in U(d)\}'=\{V^{\otimes n}:V\in SU(d)\}'=\{d\Pi(X):X\in M_d\}'=\mathrm{span}\,\pi(S_n)$.

*Démonstration.* (a) $\pi(S_n)'=\mathcal T_n$ est la proposition B.3, et $\mathcal U_n\subseteq\mathcal T_n$. Supposons
$\mathcal U_n\neq\mathcal T_n$. Alors un certain $A_0^{\otimes n}\notin\mathcal U_n$ ; soit $Y\neq0$ la composante de
$A_0^{\otimes n}$ orthogonale à $\mathcal U_n$. La fonction $p(A)=\langle Y,A^{\otimes n}\rangle$ est un polynôme homogène
de degré $n$ en les coefficients de $A$ (chaque coefficient de $A^{\otimes n}$ est un produit de $n$ coefficients de
$A$), elle s'annule sur $U(d)$, et $p(A_0)=\|Y\|^2\neq0$, ce qui contredit le lemme B.4.
Pour $V\in U(d)$, choisissons $\theta$ avec $\det V=e^{i\theta}$ ; alors $W=e^{-i\theta/d}V\in SU(d)$ et
$V^{\otimes n}=e^{in\theta/d}W^{\otimes n}$, donc les deux espaces engendrés coïncident. $\mathcal U_n$ contient
$1=1^{\otimes n}$, $V^{\otimes n}W^{\otimes n}=(VW)^{\otimes n}$ et $(V^{\otimes n})^*=(V^*)^{\otimes n}$, c'est donc une
$*$-algèbre, et toute algèbre contenant les $V^{\otimes n}$ contient l'espace qu'ils engendrent.

(b) Soit $\mathcal A$ cette algèbre. $d\Pi$ est linéaire, donc $\mathcal A$ contient $d\Pi(X)$ pour tout $X\in M_d$.
D'après le lemme B.1(b), tout $d\Pi(X)$ appartient à l'algèbre unifère $\pi(S_n)'=\mathcal T_n$, donc
$\mathcal A\subseteq\mathcal T_n$. Réciproquement, fixons $A\in M_d$ et appliquons le lemme B.5 aux $x_m=A^{(m)}$, qui
commutent deux à deux : alors $e_n=A^{\otimes n}$ et $p_j=d\Pi(A^j)\in\mathcal A$, et par récurrence sur $k$ tout
$e_k=\frac1k\sum_j(-1)^{j-1}e_{k-j}p_j$ appartient à $\mathcal A$. Donc $A^{\otimes n}\in\mathcal A$ et
$\mathcal T_n\subseteq\mathcal A$.

(c) Le commutant d'un ensemble est égal au commutant de l'algèbre qu'il engendre. D'après (a) et (b), les trois
ensembles engendrent $\mathcal T_n=\pi(S_n)'$, dont le commutant est $\pi(S_n)''$. L'espace engendré par $\pi(S_n)$ est
une $*$-algèbre ($\pi(\sigma)\pi(\tau)=\pi(\sigma\tau)$, $\pi(\sigma)^*=\pi(\sigma^{-1})$, $\pi(e)=1$), donc
$\pi(S_n)''=\mathrm{span}\,\pi(S_n)$ d'après le théorème 2.2. $\square$

*Calcul :* les dimensions du théorème B.6(a) et (c) pour de petits $d$ et $n$ — [tutoriel, section B.3](../../tutorial/fr/B_schur_weyl.md#s-B-3).

<a id="r-B-7"></a>

**Corollaire B.7 (théorème 4.14(b)).** Pour le spin total $S_a=\sum_is_a^{(i)}$ de $L$ spins $\frac12$, l'algèbre
engendrée par $S_x,S_y,S_z$ est $\{S_a\}''=\pi(S_L)'$, et $\{S_a\}'=\mathrm{span}\,\pi(S_L)$.

*Démonstration.* $M_2=\mathrm{span}\{1,s_x,s_y,s_z\}$ (les matrices de Pauli et $1$ forment une base),
$d\Pi(1)=L\cdot1$ et $d\Pi(s_a)=S_a$. D'après le théorème B.6(b), l'algèbre unifère $\mathcal A$ engendrée par
$S_x,S_y,S_z$ est $\mathcal T_L=\pi(S_L)'$. Comme les $S_a$ sont hermitiens, $\mathcal A$ est une $*$-algèbre (l'adjoint
d'un produit de générateurs est le produit dans l'ordre inverse), donc $\mathcal A=\mathcal A''\supseteq\{S_a\}''$
d'après le théorème 2.2, et $\{S_a\}''\supseteq\mathcal A$ car $\{S_a\}''$ est une algèbre unifère contenant les $S_a$.
Enfin $\{S_a\}'=\{S_a\}'''=\pi(S_L)''=\mathrm{span}\,\pi(S_L)$ (section 2.1 et théorème B.6(c)). $\square$

D'après le théorème B.6(a) pour $d=2$, la même algèbre est engendrée linéairement par les rotations $V^{\otimes L}$,
$V\in SU(2)$ ; cela ne nécessite pas la paramétrisation exponentielle de $SU(2)$ évoquée à la section 4.5.

<a id="s-B-4"></a>

## B.4 La décomposition de $(\mathbb C^d)^{\otimes n}$

<a id="r-B-8"></a>

**Théorème B.8.** Il existe $C\ge1$, des entiers $k_c,m_c\ge1$ avec $\sum_ck_cm_c=d^n$, et un unitaire
$W:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathcal H_n$ tels que
$$
W^*\pi(\sigma)W=\bigoplus_{c=1}^C\rho_c(\sigma)\otimes1_{m_c},\qquad W^*V^{\otimes n}W=\bigoplus_{c=1}^C1_{k_c}\otimes\tau_c(V)
\qquad(\sigma\in S_n,\ V\in U(d)),
$$
où les $\rho_c$ sont des représentations unitaires irréductibles de $S_n$ sur $\mathbb C^{k_c}$, deux à deux non
équivalentes, et les $\tau_c$ des représentations unitaires irréductibles de $U(d)$ sur $\mathbb C^{m_c}$, deux à deux
non équivalentes. De plus $W^*\mathrm{span}\,\pi(S_n)W=\bigoplus_cM_{k_c}\otimes1_{m_c}$,
$W^*\mathcal U_nW=\bigoplus_c1_{k_c}\otimes M_{m_c}$,
$$
\sum_ck_c^2=\dim\mathrm{span}\,\pi(S_n),\qquad\sum_cm_c^2=\binom{d^2+n-1}n,
$$
et $\sum_ck_c^2=n!$ si $n\le d$.

*Démonstration.* La proposition 2.7 pour le groupe fini $S_n$ et le théorème 2.6 pour
$\mathcal M=\mathrm{span}\,\pi(S_n)=\pi(S_n)''$ donnent $W$, les $\rho_c$ avec les propriétés énoncées, et
$W^*\pi(S_n)'W=\bigoplus_c1_{k_c}\otimes M_{m_c}$. D'après le théorème B.6(a), $\pi(S_n)'=\mathcal U_n$. Donc
$W^*V^{\otimes n}W=\bigoplus_c1\otimes\tau_c(V)$ pour certains $\tau_c(V)\in M_{m_c}$, et $V\mapsto\tau_c(V)$ est une
représentation unitaire puisque $V\mapsto V^{\otimes n}$ en est une. Comme les $V^{\otimes n}$ engendrent $\mathcal U_n$,
l'application linéaire $\sum_la_lV_l^{\otimes n}\mapsto(\sum_la_l\tau_c(V_l))_c$ est surjective sur $\bigoplus_cM_{m_c}$ ;
l'irréductibilité et la non-équivalence des $\tau_c$ s'en déduisent mot pour mot comme dans la démonstration de la
proposition 2.7. Les dimensions sont le théorème 2.6, la proposition B.3 et le lemme B.1(d). $\square$

*Calcul :* le théorème B.8 pour $S_3$ sur $(\mathbb C^3)^{\otimes3}$, décomposé par `decompose.from_group` — [tutoriel, section B.4](../../tutorial/fr/B_schur_weyl.md#s-B-4).

Pour $d=2$, les secteurs sont identifiés par le spin total :

<a id="r-B-9"></a>

**Corollaire B.9 ($L$ spins $\frac12$).** Soit $V$ un unitaire avec $V^*S_aV=\bigoplus_SS^{(S)}_a\otimes1_{m_S}$ pour le
spin total de $L$ spins $\frac12$ (corollaire 4.12, proposition 4.13). Alors
$V^*\pi(\sigma)V=\bigoplus_S1_{2S+1}\otimes\rho_S(\sigma)$, où les $\rho_S$ sont des représentations unitaires
irréductibles de $S_L$, deux à deux non équivalentes, de dimension
$$
m_S=\binom Lj-\binom L{j-1},\qquad j=\tfrac L2-S\in\{0,1,\dots,\lfloor L/2\rfloor\}.
$$
En particulier $\dim\mathrm{span}\,\pi(S_L)=\sum_Sm_S^2$ et $\dim\pi(S_L)'=\sum_S(2S+1)^2=\binom{L+3}3$.

*Démonstration.* D'après le théorème 4.3, $V^*\{S_a\}'V=\bigoplus_S1_{2S+1}\otimes M_{m_S}$, et
$\{S_a\}'=\mathrm{span}\,\pi(S_L)$ d'après le corollaire B.7. Donc $V^*\pi(\sigma)V=\bigoplus_S1\otimes\rho_S(\sigma)$ avec
une représentation unitaire $\rho_S$, et l'application $\sum_\sigma a_\sigma\pi(\sigma)\mapsto(\sum_\sigma a_\sigma\rho_S(\sigma))_S$
est surjective sur $\bigoplus_SM_{m_S}$ ; l'irréductibilité et la non-équivalence s'en déduisent comme dans la
démonstration de la proposition 2.7. Les valeurs de $m_S$ sont la proposition 4.13. D'après le théorème 4.3,
$\dim\{S_a\}'=\sum_Sm_S^2$ et $\dim\{S_a\}''=\sum_S(2S+1)^2$ ; d'après le corollaire B.7,
$\{S_a\}'=\mathrm{span}\,\pi(S_L)$ et $\{S_a\}''=\pi(S_L)'$, dont la dimension est $\binom{4+L-1}L=\binom{L+3}3$
d'après la proposition B.3 avec $d=2$. $\square$

*Calcul :* les couples $(m_S,2S+1)$ pour $L=5,6$ — [tutoriel, section B.4](../../tutorial/fr/B_schur_weyl.md#s-B-4).

<a id="r-B-10"></a>

**Corollaire B.10 (nombres de Catalan).** Pour $L\ge1$, $\displaystyle\sum_Sm_S^2=\binom{2L}L\frac1{L+1}$ ; c'est
$\dim\{S_a\}'=\dim\mathrm{span}\,\pi(S_L)$ sur $(\mathbb C^2)^{\otimes L}$.

*Démonstration.* Posons $N=2^L$ et soit $\mathrm{vec}:M_N\to\mathbb C^N\otimes\mathbb C^N$ la bijection linéaire
$|x\rangle\langle y|\mapsto x\otimes\bar y$. Alors $\mathrm{vec}(AX-XA)=(A\otimes1-1\otimes A^T)\mathrm{vec}\,X$, puisque
$A|x\rangle\langle y|B=|Ax\rangle\langle B^*y|\mapsto Ax\otimes B^T\bar y$. Ainsi $\{S_a\}'$ est envoyé sur le noyau commun
des $S_a\otimes1-1\otimes S_a^T$, $a=x,y,z$, où $S_a^T=\sum_i(s_a^T)^{(i)}$. De $s_x^T=s_x$, $s_y^T=-s_y$, $s_z^T=s_z$ et
$\sigma_ys_x\sigma_y=-s_x$, $\sigma_ys_y\sigma_y=s_y$, $\sigma_ys_z\sigma_y=-s_z$, on tire $\sigma_y(-s_a^T)\sigma_y=s_a$, donc
l'unitaire $Y=1\otimes\sigma_y^{\otimes L}$ donne $Y(S_a\otimes1-1\otimes S_a^T)Y^*=S_a\otimes1+1\otimes S_a$, le spin total de
$2L$ spins $\frac12$ sur $(\mathbb C^2)^{\otimes2L}=\mathbb C^N\otimes\mathbb C^N$. Dans la décomposition du corollaire 4.12
pour $2L$ spins, un vecteur annulé par les trois composantes est annulé par $\mathbf S^2$, qui vaut $S(S+1)\ne0$ sur
les termes de spin $S\ge\frac12$ ; réciproquement, les composantes s'annulent sur le terme de spin $0$. Le noyau
commun est donc ce terme, de dimension $m_0(2L)=\binom{2L}L-\binom{2L}{L-1}=\binom{2L}L\big(1-\frac L{L+1}\big)$
(proposition 4.13). La première égalité de dimensions est le théorème 4.3, la seconde le corollaire B.7. $\square$

*Calcul :* le corollaire B.10 pour $L=5,6$ — [tutoriel, section B.4](../../tutorial/fr/B_schur_weyl.md#s-B-4).

<a id="s-B-5"></a>

## B.5 Le calcul de Weingarten

<a id="r-B-11"></a>

**Lemme B.11 (le twirl de Haar est la projection orthogonale sur le commutant).** Soient $U$ distribué selon la loi
de Haar sur $U(d)$ et $\Phi(X)=\mathbb E[U^{\otimes n}X(U^{\otimes n})^*]$ pour $X\in\mathrm{End}\,\mathcal H_n$. Alors
$\Phi$ est la projection orthogonale de $\mathrm{End}\,\mathcal H_n$ sur
$\{V^{\otimes n}:V\in U(d)\}'=\mathrm{span}\,\pi(S_n)$.

*Démonstration.* Les coefficients de $U^{\otimes n}X(U^{\otimes n})^*$ sont des polynômes en les coefficients de $U$ et
de $\bar U$, qui sont bornés par $1$, donc $\Phi$ est une application linéaire bien définie. Pour $V\in U(d)$, $VU$ a
la loi de $U$, donc $V^{\otimes n}\Phi(X)(V^{\otimes n})^*=\mathbb E[(VU)^{\otimes n}X((VU)^{\otimes n})^*]=\Phi(X)$ : l'image
est contenue dans le commutant. Si $Y$ est dans le commutant, $U^{\otimes n}Y(U^{\otimes n})^*=Y$ et $\Phi(Y)=Y$. Donc
$\Phi$ est idempotente, d'image le commutant. D'après le lemme 5.1, $U^*=U^{-1}$ a la loi de $U$, donc
$\langle Y,\Phi(X)\rangle=\mathbb E\mathrm{Tr}\big(((U^*)^{\otimes n}YU^{\otimes n})^*X\big)=\langle\Phi(Y),X\rangle$ : $\Phi$ est
autoadjointe. Un idempotent autoadjoint est la projection orthogonale sur son image
($\langle X-\Phi(X),\Phi(Z)\rangle=\langle\Phi(X-\Phi(X)),Z\rangle=0$). Le commutant est $\mathrm{span}\,\pi(S_n)$ d'après
le théorème B.6(c). $\square$

*Calcul :* le lemme B.11 par Monte Carlo sur des unitaires de Haar — [tutoriel, section B.5](../../tutorial/fr/B_schur_weyl.md#s-B-5).

Pour écrire cette projection à l'aide de la famille génératrice $\pi(\sigma)$, qui est liée pour $n>d$, nous utilisons
le **pseudo-inverse** d'une matrice hermitienne $G$ : $G^+=f(G)$ avec $f(\lambda)=1/\lambda$ pour $\lambda\neq0$ et
$f(0)=0$ (calcul fonctionnel, théorème 1.5). Il est hermitien, $G^+=G^{-1}$ si $G$ est inversible,
$GG^+=\mathbf 1_{\neq0}(G)$ est la projection sur $(\ker G)^\perp$, $G^+v=v/\lambda$ si $Gv=\lambda v$ avec $\lambda\neq0$, et
$(R^*GR)^+=R^*G^+R$ pour $R$ unitaire (si $G=Q\Lambda Q^*$, alors $R^*GR=(R^*Q)\Lambda(R^*Q)^*$).

<a id="r-B-12"></a>

**Lemme B.12.** Soient $v_1,\dots,v_N$ des vecteurs d'un espace préhilbertien de dimension finie $\mathcal K$,
$G_{ab}=\langle v_a,v_b\rangle$ leur matrice de Gram et $W=G^+$. Alors la projection orthogonale sur
$\mathrm{span}\{v_a\}$ est $x\mapsto\sum_{a,b}v_a\,W_{ab}\,\langle v_b,x\rangle$.

*Démonstration.* Soit $T:\mathbb C^N\to\mathcal K$, $Tc=\sum_ac_av_a$. Alors $T^*x=(\langle v_b,x\rangle)_b$, $G=T^*T$, et
l'application considérée est $Q=TG^+T^*$. Comme $\langle c,Gc\rangle=\|Tc\|^2$, $\ker G=\ker T$, donc
$TG^+G=T(1-\Pi)=T$ pour la projection $\Pi$ sur $\ker T$. Donc $QT=T$ : $Q$ est l'identité sur $\mathrm{ran}\,T$. $Q$
s'annule sur $(\mathrm{ran}\,T)^\perp=\ker T^*$, et $Q^*=Q$. Ainsi $Q$ est la projection orthogonale sur
$\mathrm{ran}\,T$. $\square$

<a id="r-B-13"></a>

**Théorème B.13 (formule de Weingarten de tout ordre).** Soient $n,d\ge1$, $U$ distribué selon la loi de Haar sur
$U(d)$, $G$ la matrice $n!\times n!$ $G_{\sigma\tau}=d^{\#\text{cycles}(\sigma^{-1}\tau)}$ ($\sigma,\tau\in S_n$) et $W=G^+$.
(a) $W_{\sigma\tau}=\mathrm{Wg}(\sigma^{-1}\tau)$ pour la fonction $\mathrm{Wg}(\rho)=\mathrm{Wg}(\rho,d)=W_{e\rho}$ sur $S_n$, qui ne
dépend que du type cyclique de $\rho$ ; en particulier $\mathrm{Wg}(\rho^{-1})=\mathrm{Wg}(\rho)$.
(b) $G$ est inversible si et seulement si $n\le d$, et alors $W=G^{-1}$.
(c) Pour $X\in\mathrm{End}\,\mathcal H_n$,
$\mathbb E[U^{\otimes n}X(U^{\otimes n})^*]=\sum_{\sigma,\tau}\mathrm{Wg}(\sigma^{-1}\tau)\mathrm{Tr}(\pi(\tau)^*X)\,\pi(\sigma)$.
(d) Pour tous $\mathbf i,\mathbf j,\mathbf i',\mathbf j'\in[d]^n$,
$$
\mathbb E\big[U_{i_1j_1}\cdots U_{i_nj_n}\overline{U_{i'_1j'_1}\cdots U_{i'_nj'_n}}\big]
=\sum_{\sigma,\tau\in S_n}\prod_{m=1}^n\delta_{i_mi'_{\sigma(m)}}\delta_{j_mj'_{\tau(m)}}\;\mathrm{Wg}(\sigma\tau^{-1}),
$$
c'est-à-dire que la somme porte sur les couples tels que $\mathbf i=\mathbf i'\circ\sigma$ et $\mathbf j=\mathbf j'\circ\tau$.
(e) Les moments où les nombres $p\neq q$ de facteurs $U$ et $\bar U$ diffèrent sont nuls.

C'est le théorème 5.12 (avec $k=d$), désormais démontré.

*Démonstration.* (a) Notons $g(\rho)=d^{\#\text{cycles}(\rho)}$, de sorte que $G_{\sigma\tau}=g(\sigma^{-1}\tau)$, et soient
$\Lambda_\kappa,R_\kappa$ les matrices de permutation sur $\mathbb C^{S_n}$ telles que $\Lambda_\kappa e_\sigma=e_{\kappa\sigma}$,
$R_\kappa e_\sigma=e_{\sigma\kappa}$. Alors $(\Lambda_\kappa^*G\Lambda_\kappa)_{\sigma\tau}=g(\sigma^{-1}\kappa^{-1}\kappa\tau)=G_{\sigma\tau}$
et $(R_\kappa^*GR_\kappa)_{\sigma\tau}=g(\kappa^{-1}\sigma^{-1}\tau\kappa)=G_{\sigma\tau}$ (lemme B.1(c)). D'après les propriétés du
pseudo-inverse, $W$ possède les deux mêmes invariances : $W_{\kappa\sigma,\kappa\tau}=W_{\sigma\tau}$ et
$W_{\sigma\kappa,\tau\kappa}=W_{\sigma\tau}$. La première avec $\kappa=\sigma^{-1}$ donne $W_{\sigma\tau}=\mathrm{Wg}(\sigma^{-1}\tau)$ ;
la seconde donne $\mathrm{Wg}(\kappa^{-1}\rho\kappa)=\mathrm{Wg}(\rho)$, donc $\mathrm{Wg}$ est constante sur les classes de
conjugaison, qui sont les types cycliques (lemme B.1(c)) ; $\rho^{-1}$ a le type cyclique de $\rho$.
(b) $G$ est la matrice de Gram des $\pi(\sigma)$ (lemme B.1(c)), et avec $T$ comme dans la démonstration du lemme B.12,
$G=T^*T$ est inversible si et seulement si $T$ est injective, c'est-à-dire si et seulement si les $\pi(\sigma)$ sont
linéairement indépendants : si et seulement si $n\le d$ d'après le lemme B.1(d).
(c) Lemme B.11 et lemme B.12 avec $\mathcal K=\mathrm{End}\,\mathcal H_n$ et $v_\sigma=\pi(\sigma)$.
(d) Les coefficients d'un produit de Kronecker sont $\langle e_{\mathbf i},U^{\otimes n}e_{\mathbf j}\rangle=\prod_mU_{i_mj_m}$,
donc le moment vaut $\langle e_{\mathbf i},\Phi(X)e_{\mathbf i'}\rangle$ pour $X=|e_{\mathbf j}\rangle\langle e_{\mathbf j'}|$. Ici
$\mathrm{Tr}(\pi(\tau)^*X)=\langle\pi(\tau)e_{\mathbf j'},e_{\mathbf j}\rangle=\delta_{\mathbf j,\mathbf j'\circ\tau^{-1}}$ et
$\langle e_{\mathbf i},\pi(\sigma)e_{\mathbf i'}\rangle=\delta_{\mathbf i,\mathbf i'\circ\sigma^{-1}}$. D'après (c), le moment vaut
$\sum_{\sigma,\tau}\delta_{\mathbf i,\mathbf i'\circ\sigma^{-1}}\delta_{\mathbf j,\mathbf j'\circ\tau^{-1}}\mathrm{Wg}(\sigma^{-1}\tau)$ ;
la substitution $\sigma\to\sigma^{-1}$, $\tau\to\tau^{-1}$ donne la formule.
(e) Pour $\theta$ réel, $e^{i\theta}1\in U(d)$, donc $e^{i\theta}U$ a la loi de $U$, et le moment est égal à
$e^{i(p-q)\theta}$ fois lui-même. On prend $\theta$ tel que $e^{i(p-q)\theta}\neq1$. $\square$

Pour $n=2\le d$, $G=\begin{pmatrix}d^2&d\\d&d^2\end{pmatrix}$ et $G^{-1}$ donne $\mathrm{Wg}(e)=\frac1{d^2-1}$,
$\mathrm{Wg}(\tau)=-\frac1{d(d^2-1)}$ : c'est le théorème 5.10. Pour $d=1$, tous les $G_{\sigma\tau}$ valent $1$,
$W=\frac1{(n!)^2}G$, et le moment de la phase $U=U_{11}$ vaut $\sum_{\sigma,\tau}\frac1{(n!)^2}=1=|U_{11}|^{2n}$, comme il se
doit.

*Calcul :* le théorème B.13, avec $\mathrm{Wg}$ calculée à partir du pseudo-inverse, y compris les moments d'ordre $3$ pour $d=2$, où $G$ est singulière — [tutoriel, section B.5](../../tutorial/fr/B_schur_weyl.md#s-B-5).

<a id="r-B-14"></a>

**Corollaire B.14.** Pour tous $n,d\ge1$, $\mathbb E|U_{11}|^{2n}=\dfrac{n!}{d(d+1)\cdots(d+n-1)}$.

*Démonstration.* Avec tous les indices égaux à $1$ dans le théorème B.13(d), le moment vaut $\sum_{\sigma,\tau}\mathrm{Wg}(\sigma\tau^{-1})=\sum_{\sigma,\tau}W_{\sigma^{-1},\tau^{-1}}=\mathbf 1^TW\mathbf 1$
(partie (a)). Chaque somme de ligne de $G$ vaut $c_n(d)=\sum_{\rho\in S_n}d^{\#\text{cycles}(\rho)}\ge d^n>0$, donc
$G\mathbf 1=c_n(d)\mathbf 1$, $W\mathbf 1=\mathbf 1/c_n(d)$ et le moment vaut $n!/c_n(d)$. Enfin
$c_n(x)=x(x+1)\cdots(x+n-1)$ (chapitre 5, exercice 2, solution). $\square$

Cela étend l'exercice 2 du chapitre 5 au cas $n>d$, où $G$ est singulière.

*Calcul :* le corollaire B.14 pour $n=3$, $d=2$ — [tutoriel, section B.5](../../tutorial/fr/B_schur_weyl.md#s-B-5).

<a id="s-B-6"></a>

## B.6 Opérateurs tensoriels et théorème de Wigner–Eckart pour $SU(2)$

La section 9.3 construit des corrélations de spin à partir du théorème de Wigner–Eckart. Nous le démontrons pour les
représentations de spin de la section 4.5. Le spin $k$ est la représentation $S^{(k)}_a$ du théorème 4.11 sur
$\mathbb C^{2k+1}$, de base $|k,q\rangle$, $q=k,\dots,-k$ ; rappelons que
$S^{(k)}_\pm|k,q\rangle=\sqrt{k(k+1)-q(q\pm1)}\,|k,q\pm1\rangle$ (pour $S_-=S_+^*$, cela découle de la formule pour $S_+$).

**Définition.** Soient $(S_a)$ une représentation de spin sur $\mathbb C^N$ et $k\in\{0,\frac12,1,\dots\}$. Un
**opérateur tensoriel de rang $k$** est une famille $(T_q)_{q=-k,\dots,k}$ de matrices $N\times N$ telle que
$$
[S_z,T_q]=q\,T_q,\qquad[S_\pm,T_q]=\sqrt{k(k+1)-q(q\pm1)}\;T_{q\pm1}
$$
(le coefficient s'annule pour $q=\pm k$, où $T_{\pm(k+1)}$ n'est pas défini). Un **opérateur vectoriel** est un
triplet $(T_x,T_y,T_z)$ tel que $[S_a,T_b]=i\sum_c\varepsilon_{abc}T_c$, c'est-à-dire $[S_x,T_y]=iT_z$, $[S_x,T_z]=-iT_y$
et de même par permutation circulaire, et $[S_a,T_a]=0$.

<a id="r-B-15"></a>

**Lemme B.15.** (a) Si $(T_a)$ est un opérateur vectoriel, ses **composantes sphériques**
$T_{\pm1}=\mp\frac1{\sqrt2}(T_x\pm iT_y)$, $T_0=T_z$ forment un opérateur tensoriel de rang $1$, et pour deux opérateurs
vectoriels $\sum_aT_aT'_a=\sum_q(-1)^qT_qT'_{-q}$.
(b) Pour $L$ spins $\frac12$ et le spin total $S_a$, chaque spin de site $(s^{(i)}_a)_a$ ainsi que $(S_a)_a$ lui-même sont
des opérateurs vectoriels.
(c) $(T_q)$ est un opérateur tensoriel de rang $k$ si et seulement si l'application
$\Psi:\mathbb C^N\otimes\mathbb C^{2k+1}\to\mathbb C^N$, $\Psi(v\otimes|k,q\rangle)=T_qv$, vérifie
$\Psi\,(S_a\otimes1+1\otimes S^{(k)}_a)=S_a\Psi$ pour $a=x,y,z$.

*Démonstration.* (a) Des relations, on tire $[S_z,T_x\pm iT_y]=\pm(T_x\pm iT_y)$, $[S_\pm,T_z]=\mp(T_x\pm iT_y)$,
$[S_\pm,T_x\mp iT_y]=\pm2T_z$ et $[S_\pm,T_x\pm iT_y]=0$ (par exemple $[S_x+iS_y,T_x-iT_y]=-i[S_x,T_y]+i[S_y,T_x]=T_z+T_z$).
En composantes sphériques : $[S_z,T_q]=qT_q$, $[S_\pm,T_0]=\sqrt2\,T_{\pm1}$, $[S_\pm,T_{\mp1}]=\sqrt2\,T_0$,
$[S_\pm,T_{\pm1}]=0$, qui sont les relations de définition pour $k=1$ ($\sqrt{2-q(q\pm1)}$ vaut $\sqrt2$ pour $q=0$ et
$q=\mp1$, et $0$ pour $q=\pm1$). Pour le produit,
$T_{+1}T'_{-1}+T_{-1}T'_{+1}=-\frac12\big((T_x+iT_y)(T'_x-iT'_y)+(T_x-iT_y)(T'_x+iT'_y)\big)=-(T_xT'_x+T_yT'_y)$.
(b) $[S_a,s^{(i)}_b]=[s^{(i)}_a,s^{(i)}_b]$ (les autres sites commutent), et les matrices de spin $\frac12$ vérifient les
relations de spin ; pour $S_a$, on somme sur $i$. (c)
$\Psi(S_a\otimes1+1\otimes S^{(k)}_a)(v\otimes|k,q\rangle)=T_qS_av+\sum_{q'}(S^{(k)}_a)_{q'q}T_{q'}v$, donc la relation
d'entrelacement pour $S_a$ s'écrit $[S_a,T_q]=\sum_{q'}(S^{(k)}_a)_{q'q}T_{q'}$ pour tout $q$. Pour $S_z,S_\pm$, ce sont les
relations de définition, et les relations pour $S_x,S_y,S_z$ et pour $S_z,S_\pm$ sont combinaisons linéaires les unes
des autres. $\square$

*Calcul :* le lemme B.15(a) et (b) pour les spins de site de six spins $\frac12$ — [tutoriel, section B.6](../../tutorial/fr/B_schur_weyl.md#s-B-6).

<a id="r-B-16"></a>

**Théorème B.16 (Wigner–Eckart).** Soient $(S_a)$ une représentation de spin sur $\mathbb C^N$, $V$ un unitaire avec
$V^*S_aV=\bigoplus_SS^{(S)}_a\otimes1_{m_S}$ (corollaire 4.12), $|S,M;b\rangle=V(|S,M\rangle\otimes f_b)$ la base orthonormée
correspondante ($f_b$ la base canonique de $\mathbb C^{m_S}$), et $(T_q)$ un opérateur tensoriel de rang $k$. On suppose
que $S,S'$ apparaissent ($m_S,m_{S'}\ge1$).
(a) Si $S'\notin\{|S-k|,|S-k|+1,\dots,S+k\}$, alors $\langle S',M';a|T_q|S,M;b\rangle=0$ pour tous $M',M,q,a,b$.
(b) Sinon, les opérateurs d'entrelacement de spin $S\otimes$ spin $k$ vers spin $S'$ forment un espace de dimension un ;
pour tout tel $C=C^{S'}_{S,k}:\mathbb C^{2S+1}\otimes\mathbb C^{2k+1}\to\mathbb C^{2S'+1}$ non nul, il existe une unique
matrice $m_{S'}\times m_S$ $R=\langle S'\|T\|S\rangle$, l'**élément de matrice réduit**, telle que
$$
\langle S',M';a|\,T_q\,|S,M;b\rangle=\langle S',M'|\,C\,(|S,M\rangle\otimes|k,q\rangle)\rangle\;R_{ab}\qquad\text{for all }M',M,q,a,b .
$$
(c) L'élément est nul sauf si $M'=M+q$. Dans le cas (b), il existe un $q$ tel que
$\langle S',S'|C(|S,S'-q\rangle\otimes|k,q\rangle)\rangle\neq0$, de sorte que $R$, et avec lui tout le bloc de $T$ entre les
multiplets $S$ et $S'$, est déterminé par les éléments de matrice $\langle S',S';a|T_q|S,S'-q;b\rangle$ entre les
vecteurs de plus haut poids de $S'$ et un vecteur de chaque copie de $S$.

*Démonstration.* Fixons $a\le m_{S'}$, $b\le m_S$, et soient $J:\mathbb C^{2S+1}\to\mathbb C^N$, $Jx=V(x\otimes f_b)$ (dans
le terme $S$), et $K:\mathbb C^N\to\mathbb C^{2S'+1}$ l'adjoint de $x\mapsto V(x\otimes f_a)$ (dans le terme $S'$). Notons
$\alpha\in\{x,y,z\}$ l'indice de composante. D'après la forme de $V^*S_\alpha V$, $S_\alpha J=JS^{(S)}_\alpha$ et, en prenant
les adjoints dans la même relation pour $S'$ (les matrices de spin sont hermitiennes), $KS_\alpha=S^{(S')}_\alpha K$. Soit
$t=K\Psi(J\otimes1)$ avec $\Psi$ du lemme B.15(c) ; ses éléments de matrice sont
$\langle S',M'|t(|S,M\rangle\otimes|k,q\rangle)\rangle=\langle S',M';a|T_q|S,M;b\rangle$, et
$$
t\,(S^{(S)}_\alpha\otimes1+1\otimes S^{(k)}_\alpha)=K\Psi(S_\alpha\otimes1+1\otimes S^{(k)}_\alpha)(J\otimes1)=KS_\alpha\Psi(J\otimes1)=S^{(S')}_\alpha\,t,
$$
donc $t$ est un opérateur d'entrelacement de $\rho=$ spin $S\otimes$ spin $k$ (la représentation de spin
$S^{(S)}_\alpha\otimes1+1\otimes S^{(k)}_\alpha$, section 4.5) vers $\sigma=$ spin $S'$. Le passage à l'adjoint est une
bijection linéaire entre $\mathrm{Hom}(\rho,\sigma)$ et $\mathrm{Hom}(\sigma,\rho)$ (l'adjoint de
$t\rho(\alpha)=\sigma(\alpha)t$ est $\rho(\alpha)t^*=t^*\sigma(\alpha)$, puisque toutes les matrices en jeu sont
hermitiennes). D'après le corollaire 4.12 et l'exercice 2 du chapitre 4 (démontré dans sa solution à partir du seul
corollaire 4.12), $\rho$ est unitairement équivalente à $\bigoplus_j(\text{spin }j)$, une copie de chaque
$j\in\{|S-k|,|S-k|+1,\dots,S+k\}$. Ces spins sont irréductibles et deux à deux non équivalents (théorème 4.11), il
s'agit donc d'une décomposition comme dans le théorème 4.3, et son dernier énoncé donne
$\dim\mathrm{Hom}(\sigma,\rho)=1$ si $S'$ est l'un de ces $j$, et $\mathrm{Hom}(\sigma,\rho)=0$ sinon (le spin
$S'$ n'est alors équivalent à aucun d'eux, d'après le théorème 4.11). Cela donne (a), et dans le cas (b) $t=R_{ab}C$
pour un unique nombre $R_{ab}$.
(c) La relation d'entrelacement pour $S_z$ donne $M'\langle S',M'|t(|S,M\rangle\otimes|k,q\rangle)\rangle=(M+q)\langle S',M'|t(|S,M\rangle\otimes|k,q\rangle)\rangle$.
L'image de $C\neq0$ est un sous-espace invariant non nul du spin irréductible $S'$, donc l'espace tout entier, de sorte
que $\langle S',S'|C(|S,M\rangle\otimes|k,q\rangle)\rangle\neq0$ pour un certain vecteur de base, et alors $M=S'-q$ par la
même relation pour l'opérateur d'entrelacement $C$. On divise par ce nombre. $\square$

*Calcul :* le théorème B.16(a) et (b) sous une forme indépendante de la base, dans la décomposition calculée par `decompose.from_generators` — [tutoriel, section B.6](../../tutorial/fr/B_schur_weyl.md#s-B-6).

D'après le lemme B.15(a), le bloc entre les multiplets $S$ et $S'$ d'un produit scalaire tel que
$\mathbf s^{(i)}\cdot\mathbf s^{(j)}=\sum_q(-1)^q s^{(i)}_qs^{(j)}_{-q}$ est
$\sum_{S''}\sum_q(-1)^q(\text{bloc }S''\to S'\text{ de }s^{(i)}_q)(\text{bloc }S\to S''\text{ de }s^{(j)}_{-q})$, et chaque facteur est
donné par le théorème B.16 au moyen d'un élément de matrice réduit et d'une application de Clebsch–Gordan.

*Dans la bibliothèque :* c'est le principe sur lequel repose `su2.spin_correlation`, que le chapitre 9 vérifie par comparaison avec la méthode `'direct'` — [tutoriel, section 9.3](../../tutorial/fr/09_spin_chains.md#s-9-3).

<a id="s-B-7"></a>

## B.7 Mesures de Haar utilisées dans le cours

Le cours utilise les mesures de Haar des groupes suivants : $U(k)$ (théorème 5.3), $U(1)$ (le cas $k=1$ ; chapitre 5,
exercice 3), $U(\mathcal M)\cong\prod_cU(k_c)$ (corollaire 5.4), $SU(k)$ (section 5.1, après le corollaire 5.4), le groupe
orthogonal $O(k)$ (section 5.1, chapitre 11, exercice 2, et algèbres réelles), et les groupes diagonaux $T$ du
lemme C.2. L'unicité, pour tous ces groupes, est le lemme 5.1, puisque chacun est un sous-groupe fermé d'un certain
$U(n)$ ; l'existence est une construction explicite dans chaque cas (pour $T$, dans le lemme C.2), et la construction
pour $O(k)$ est la suivante. Aucun théorème général d'existence de la mesure de Haar d'un groupe compact n'est utilisé
nulle part dans le cours.

<a id="r-B-17"></a>

**Théorème B.17 (mesure de Haar du groupe orthogonal).** Soit $Z$ une matrice réelle $k\times k$ à coefficients
$N(0,1)$ indépendants. Alors $Z$ est inversible avec probabilité un, $Z=OR$ avec $O$ orthogonale et $R$ réelle
triangulaire supérieure à diagonale positive, et $O$ est distribuée selon la loi de Haar sur $O(k)$. Si $Z=O'R'$ est
une décomposition QR réelle quelconque ($O'$ orthogonale, $R'$ réelle triangulaire supérieure), alors
$O'\mathrm{diag}(\mathrm{sgn}\,R'_{ii})=O$. Pour $O$ distribuée selon la loi de Haar et $X\in M_k(\mathbb C)$,
$\mathbb E[O]=0$ et $\mathbb E[OXO^T]=\frac{\mathrm{Tr}\,X}k1$ ; en particulier
$\mathbb E[O_{ij}O_{lm}]=\delta_{il}\delta_{jm}/k$.

*Démonstration.* $\det Z$ est un polynôme en les $k^2$ coefficients réels, non identiquement nul ($\det1=1$), et $Z$
a une densité, donc $\det Z\neq0$ presque sûrement (corollaire A.7(b), comme dans la démonstration du théorème 5.3).
Gram–Schmidt appliqué à des colonnes réelles produit des $q_j$ réels et des coefficients réels, c'est-à-dire une
décomposition réelle $Z=OR$ avec $R_{jj}>0$ ; par l'unicité du lemme 5.2, c'est la décomposition du lemme 5.2, donc
$O$ dépend continûment de $Z$. La densité $(2\pi)^{-k^2/2}e^{-\|Z\|_2^2/2}$ de $Z$ est invariante par $Z\mapsto VZ$ pour
$V\in O(k)$, application linéaire de $\mathbb R^{k^2}$ qui préserve $\|Z\|_2^2=\mathrm{Tr}(Z^TZ)$, donc orthogonale
de jacobien $1$ ; ainsi $VZ$ a la loi de $Z$ (proposition A.14(a)). Comme $VZ=(VO)R$, l'unicité donne $O(VZ)=VO(Z)$,
donc la loi de $O$ est invariante à gauche sur le groupe compact $O(k)$, et le lemme 5.1 s'applique. La correction de
signe est la fin de la démonstration du théorème 5.3, avec $\mathrm{sgn}\,R'_{ii}=R'_{ii}/|R'_{ii}|$
($R'_{ii}\neq0$ car $Z$ est inversible). Comme $-1\in O(k)$, $-O$ a la loi de $O$, donc
$\mathbb E[O]=-\mathbb E[O]=0$ (les coefficients sont bornés par $1$). Enfin $T=\mathbb E[OXO^T]$ vérifie $VTV^T=T$ pour
$V\in O(k)$, puisque $VO$ a la loi de $O$. Avec pour $V$ la matrice diagonale ayant $-1$ en position $i$ et $1$
ailleurs, $T_{ij}=-T_{ij}$ pour $j\neq i$ ; avec des matrices de permutation, tous les $T_{ii}$ sont égaux (comme dans
le lemme 5.6). Donc $T=\lambda1$, et $k\lambda=\mathrm{Tr}\,T=\mathrm{Tr}\,X$. Pour $X=|e_j\rangle\langle e_m|$,
$(OXO^T)_{il}=O_{ij}O_{lm}$. $\square$

*Calcul :* les moments d'ordre deux du théorème B.17 par Monte Carlo, et pourquoi la correction de signe est nécessaire — [tutoriel, section B.7](../../tutorial/fr/B_schur_weyl.md#s-B-7).

**Retour :** [chapitre 4](04_symmetries.md) (théorème 4.14), [chapitre 5](05_random_matrices.md) (théorème 5.12),
[chapitre 9](09_spin_chains.md) (corrélations).
