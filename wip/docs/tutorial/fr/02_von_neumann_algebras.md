[English](../02_von_neumann_algebras.md) · [Español](../es/02_von_neumann_algebras.md) · **Français** · [Deutsch](../de/02_von_neumann_algebras.md) · [中文](../zh/02_von_neumann_algebras.md) · [日本語](../ja/02_von_neumann_algebras.md) · [Русский](../ru/02_von_neumann_algebras.md)

# 2. Algèbres de von Neumann de dimension finie

Une algèbre de von Neumann est une algèbre d'opérateurs stable par passage au « double
commutant ». En dimension finie, c'est une condition purement algébrique, et le théorème
principal de ce chapitre affirme qu'une telle algèbre est toujours, à un changement de
base unitaire près, une somme directe d'algèbres de matrices pleines, chacune répétée avec
une certaine multiplicité. C'est sur cette structure que repose la bibliothèque.

**Vous apprendrez** le théorème du bicommutant, le théorème de structure avec une
démonstration complète, la manière dont un groupe de symétrie engendre une telle algèbre,
ce qu'est le centre, et comment la bibliothèque représente tout cela.

Dans tout ce chapitre, une *$*$-algèbre sur $\mathbb C^n$* désigne un sous-espace
vectoriel complexe $\mathcal M\subseteq M_n(\mathbb C)$ stable par produit et par passage à
l'adjoint, et contenant l'identité $1$.

## 2.1 Commutants et théorème du bicommutant

Pour une partie $S\subseteq M_n(\mathbb C)$, le **commutant** est
$S'=\{X\in M_n(\mathbb C):XA=AX\text{ pour tout }A\in S\}$ (l'ensemble des $X$ tels que
$XA=AX$ pour tout $A\in S$). C'est toujours une algèbre unifère ; si $S$ est stable par
passage à l'adjoint ($S^*=S$), alors $S'$ est une $*$-algèbre, puisque $XA=AX$ pour tout
$A\in S$ entraîne $X^*A^*=A^*X^*$ pour tout $A\in S$. Tout élément de $S$ commute avec tout
élément de $S'$, donc $S\subseteq S''$ ; appliqué à $S'$, cela donne $S'\subseteq S'''$, et
en passant aux commutants dans $S\subseteq S''$ (ce qui renverse les inclusions), on obtient
$S'''\subseteq S'$. Par conséquent $S'=S'''$.

Lorsqu'une algèbre agit sur un sous-espace $E$ de $\mathbb C^n$, on identifie $E$ à
$\mathbb C^{\dim E}$ au moyen d'une base orthonormée, de sorte que tous les résultats
énoncés sur $\mathbb C^n$ s'y appliquent.

**Définition.** Une **algèbre de von Neumann** sur $\mathbb C^n$ est une $*$-algèbre
$\mathcal M$ telle que $\mathcal M=\mathcal M''$. Son **centre** est
$Z(\mathcal M)=\mathcal M\cap\mathcal M'$. Un **facteur** est une algèbre de von Neumann
telle que $Z(\mathcal M)=\mathbb C1$.

La définition est faite pour des opérateurs sur des espaces de Hilbert de dimension infinie, que ce
cours ne traite pas. En dimension finie, la condition $\mathcal M=\mathcal M''$ est automatique :

**Lemme 2.1.** Soit $\mathcal M$ une $*$-algèbre sur $\mathbb C^n$ et
$K\subseteq\mathbb C^n$ un sous-espace tel que $\mathcal MK\subseteq K$. Alors le projecteur
orthogonal $P$ sur $K$ appartient à $\mathcal M'$.

*Démonstration.* Pour $A\in\mathcal M$, on a $AK\subseteq K$ et $A^*K\subseteq K$ (car
$A^*\in\mathcal M$). La seconde inclusion donne $AK^\perp\subseteq K^\perp$ : pour
$y\in K^\perp$ et $x\in K$, $\langle x,Ay\rangle=\langle A^*x,y\rangle=0$. Ainsi $A$ envoie
$K$ dans $K$ et $K^\perp$ dans $K^\perp$, ce qui signifie $AP=PAP=PA$. $\square$

**Théorème 2.2 (théorème du bicommutant).** Toute $*$-algèbre $\mathcal M$ sur
$\mathbb C^n$ vérifie $\mathcal M=\mathcal M''$.

*Démonstration.* On a $\mathcal M\subseteq\mathcal M''$ ; soit $T\in\mathcal M''$.

*Étape 1 : un seul vecteur.* Pour $x\in\mathbb C^n$, le sous-espace $K=\mathcal Mx$ est
invariant par $\mathcal M$ et contient $x=1x$. D'après le lemme 2.1, son projecteur $P$
appartient à $\mathcal M'$, donc $TP=PT$ et $Tx=TPx=PTx\in K$ : il existe $A\in\mathcal M$
tel que $Tx=Ax$.

*Étape 2 : tous les vecteurs à la fois.* Sur $\mathbb C^n\otimes\mathbb C^n=(\mathbb C^n)^{\oplus n}$,
considérons la $*$-algèbre $\tilde{\mathcal M}=\{A\oplus\dots\oplus A:A\in\mathcal M\}$. Une
matrice par blocs $[X_{ij}]_{i,j=1}^n$ commute avec tous les $A\oplus\dots\oplus A$ si et
seulement si $X_{ij}A=AX_{ij}$ pour tous $i,j$, donc
$\tilde{\mathcal M}'=\{[X_{ij}]:X_{ij}\in\mathcal M'\}$. L'opérateur
$\tilde T=T\oplus\dots\oplus T$ commute avec toutes ces matrices $[X_{ij}]$, car
$TX_{ij}=X_{ij}T$. Par conséquent $\tilde T\in\tilde{\mathcal M}''$. Appliquons l'étape 1 à
$\tilde{\mathcal M}$, à $\tilde T$ et au vecteur $\xi=(e_1,\dots,e_n)$ : il existe
$A\in\mathcal M$ tel que $(Te_1,\dots,Te_n)=(Ae_1,\dots,Ae_n)$, c'est-à-dire
$T=A\in\mathcal M$. $\square$

Ainsi, en dimension finie, « algèbre de von Neumann » et « $*$-algèbre unifère de
matrices » désignent la même chose. Le théorème du bicommutant n'en reste pas moins le bon
point de vue : en physique, une algèbre d'observables est très souvent *donnée* comme un
commutant — celui des opérateurs qui commutent avec une symétrie (section 2.4).

## 2.2 Le centre

**Lemme 2.3.** Soit $\mathcal Z$ une $*$-algèbre commutative sur $\mathbb C^n$. Il existe
des projecteurs $P_1,\dots,P_C\in\mathcal Z$, non nuls, deux à deux orthogonaux, avec
$\sum_cP_c=1$, tels que $\mathcal Z=\operatorname{span}\{P_1,\dots,P_C\}$. Ce sont les
projecteurs non nuls minimaux de $\mathcal Z$, et ils sont uniques à l'ordre près.

*Démonstration.* Tout $z\in\mathcal Z$ est normal ($z^*\in\mathcal Z$ commute avec $z$).
Soit $z_1,\dots,z_r$ une famille génératrice de $\mathcal Z$. Décomposons $\mathbb C^n$ en
sous-espaces propres communs. Les sous-espaces propres de $z_1$ sont deux à deux orthogonaux
et engendrent $\mathbb C^n$ (théorème 1.3) ; ils sont invariants par $z_2,\dots,z_r$ et par
leurs adjoints (qui commutent avec $z_1$), donc chaque $z_j$ restreint à un sous-espace
propre de $z_1$ est encore normal, et l'on peut décomposer chaque sous-espace propre de
$z_1$ en sous-espaces propres de $z_2$, et ainsi de suite. On obtient une décomposition
orthogonale $\mathbb C^n=E_1\oplus\dots\oplus E_C$ en sous-espaces non nuls sur lesquels
chaque $z\in\mathcal Z$ agit comme un scalaire $\chi_c(z)$, où l'on a regroupé les
sous-espaces correspondant à la même fonction $\chi_c:\mathcal Z\to\mathbb C$, de sorte que
$\chi_c\neq\chi_d$ pour $c\neq d$. Soit $P_c$ le projecteur sur $E_c$. Alors
$z=\sum_c\chi_c(z)P_c$ pour tout $z\in\mathcal Z$.

Il reste à montrer que $P_c\in\mathcal Z$. Pour $d\neq c$, choisissons $z_{cd}\in\mathcal Z$
tel que $\chi_c(z_{cd})\neq\chi_d(z_{cd})$. Le produit
$\prod_{d\neq c}\big(z_{cd}-\chi_d(z_{cd})1\big)/\big(\chi_c(z_{cd})-\chi_d(z_{cd})\big)$
appartient à $\mathcal Z$, agit comme $1$ sur $E_c$ et comme $0$ sur chaque $E_d$,
$d\neq c$ : c'est $P_c$ (pour $C=1$, $P_1=1$ : le produit vide vaut $1$). Minimalité et unicité : un projecteur $\sum_c a_cP_c$ vérifie
$a_c\in\{0,1\}$, donc les projecteurs non nuls de $\mathcal Z$ sont les sommes des $P_c$
indexées par des parties non vides. $\square$

Appliqués à $\mathcal Z=Z(\mathcal M)$, les projecteurs $P_c$ sont les **projecteurs
centraux minimaux** de $\mathcal M$. Tout $A\in\mathcal M$ se décompose en
$A=\sum_cAP_c$ avec $AP_c=P_cAP_c$, donc
$\mathcal M=\bigoplus_c\mathcal MP_c$ (somme directe orthogonale d'algèbres agissant sur
$E_c=P_c\mathbb C^n$).

**Lemme 2.4.** Chaque $\mathcal M_c=\mathcal MP_c$, vue comme algèbre sur $E_c$ d'unité
$P_c$, est un facteur.

*Démonstration.* Soit $y\in\mathcal M_c$ commutant avec $\mathcal M_c$. Pour $A\in\mathcal M$,
$Ay=AP_cy=(AP_c)y=y(AP_c)=yP_cA=yA$, puisque $P_c$ est central et $y=yP_c=P_cy$. Donc
$y\in Z(\mathcal M)$, $y=\sum_d\chi_d(y)P_d$, et $y=yP_c=\chi_c(y)P_c$. $\square$

## 2.3 Le théorème de structure

**Lemme 2.5 (facteurs).** Soit $\mathcal N$ un facteur sur un espace $E$. Il existe des
entiers $k,m\ge1$ et un opérateur unitaire $V:\mathbb C^k\otimes\mathbb C^m\to E$ tels que
$V^*\mathcal NV=M_k(\mathbb C)\otimes1_m$. De plus $V^*\mathcal N'V=1_k\otimes M_m(\mathbb C)$.

*Démonstration.* *Étape 1 : un projecteur minimal.* $\mathcal N$ contient des projecteurs
non nuls (par exemple $1$). Soit $e\in\mathcal N$ un projecteur non nul de rang minimal.
Montrons que $e\mathcal Ne=\mathbb Ce$. Soit $x=x^*\in e\mathcal Ne$. Chaque projecteur
spectral $q$ de $x$ associé à une valeur propre non nulle $\lambda$ s'écrit $q=p(x)$, où
$p$ est le polynôme de Lagrange tel que $p(\lambda)=1$ et $p(\mu)=0$ pour
$\mu\in(\operatorname{spec}x\cup\{0\})\setminus\{\lambda\}$ (théorème 1.5(a),(b) : $p$ coïncide
sur $\operatorname{spec}x$ avec la fonction indicatrice de $\{\lambda\}$). Comme $p(0)=0$, on a
$p(t)=t\,r(t)$ pour un polynôme $r$, donc $q=x\,r(x)=r(x)\,x\in\mathcal N$, et $x=ex=xe$ donne
$q=eq=qe$ ; l'image de $q$ est donc contenue dans celle de $e$. Si $q\neq0$, la minimalité donne
$\operatorname{rank}q\ge\operatorname{rank}e$, de sorte que les deux images coïncident et $q=e$. Comme
les projecteurs spectraux associés à des valeurs propres distinctes sont orthogonaux, $x$ a au plus une
valeur propre non nulle $\lambda$, de projecteur spectral $e$, donc $x=\lambda e$ (et $x=0=0\cdot e$ s'il
n'y en a aucune). Un élément quelconque de $e\mathcal Ne$ s'écrit $x_1+ix_2$ avec
$x_{1,2}=x_{1,2}^*\in e\mathcal Ne$.

*Étape 2 : $e$ engendre $E$.* Soit $F=\operatorname{span}\mathcal Ne E$. Ce sous-espace est
invariant par $\mathcal N$, donc son projecteur $p_F$ appartient à $\mathcal N'$
(lemme 2.1). Il est aussi invariant par la $*$-algèbre $\mathcal N'$ : pour
$b\in\mathcal N'$, $bAe\xi=Aeb\xi\in F$ ; donc, d'après le lemme 2.1 appliqué à
$\mathcal N'$, $p_F\in\mathcal N''=\mathcal N$ (théorème 2.2). Par conséquent
$p_F\in Z(\mathcal N)=\mathbb C1$, et $p_F\neq0$ (car $eE\subseteq F$), donc $F=E$.

*Étape 3 : unités matricielles.* Appelons *$e$-isométrie* un élément $v\in\mathcal N$ tel
que $v^*v=e$ ; alors $vv^*$ est un projecteur de $\mathcal N$. (D'abord $v=ve$ : pour tout
$x$, $\|v(1-e)x\|^2=\langle(1-e)x,v^*v(1-e)x\rangle=\langle(1-e)x,e(1-e)x\rangle=0$. Donc
$vv^*v=ve=v$ et $(vv^*)^2=vev^*=vv^*$. Pour les $v_i$ ci-dessous, cela donne $e_iv_i=v_i$
et $v_i^*e_i=v_i^*$.) Choisissons des $e$-isométries $v_1=e,v_2,\dots,v_k$ d'images
$e_i=v_iv_i^*$ deux à deux orthogonales, avec $k$ maximal (les images sont des
sous-espaces orthogonaux de même dimension $\operatorname{rank}e$, donc
$k\le\dim E/\operatorname{rank}e$). Montrons que $\sum_ie_i=1$. Sinon,
$g=1-\sum_ie_i\neq0$ est un projecteur de $\mathcal N$. D'après l'étape 2, $gAe\neq0$ pour
un certain $A\in\mathcal N$ (si $gAe=0$ pour tout $A$, alors $gAe\xi=0$ pour tous $A,\xi$,
donc $g$ s'annule sur $F=E$). Posons $y=gAe$. Alors $y^*y\in e\mathcal Ne$ est positif et
non nul, donc $y^*y=\lambda e$ avec $\lambda>0$, et $v=y/\sqrt\lambda$ est une
$e$-isométrie d'image contenue dans $gE$, orthogonale à tous les $e_i$ — ce qui contredit
la maximalité.

Pour $i\neq j$, $v_i^*v_j=v_i^*e_ie_jv_j=0$. Posons $E_{ij}=v_iv_j^*$. Alors
$E_{ij}E_{kl}=v_iv_j^*v_kv_l^*=\delta_{jk}v_iev_l^*=\delta_{jk}E_{il}$,
$E_{ij}^*=E_{ji}$ et $\sum_iE_{ii}=1$. Pour $A\in\mathcal N$,
$A=\sum_{ij}e_iAe_j=\sum_{ij}v_i(v_i^*Av_j)v_j^*$ et $v_i^*Av_j\in e\mathcal Ne=\mathbb Ce$,
donc $A=\sum_{ij}a_{ij}E_{ij}$ avec $a_{ij}\in\mathbb C$. Ainsi
$\mathcal N=\operatorname{span}\{E_{ij}\}$.

*Étape 4 : l'opérateur unitaire.* Soit $m=\operatorname{rank}e$ et $f_1,\dots,f_m$ une base
orthonormée de $eE$. Les vecteurs $v_if_r$ sont orthonormés :
$\langle v_if_r,v_jf_s\rangle=\langle f_r,v_i^*v_jf_s\rangle=\delta_{ij}\langle f_r,ef_s\rangle=\delta_{ij}\delta_{rs}$.
Ils engendrent $E$, puisque $E=\bigoplus_ie_iE$ et $e_iE=v_ieE$. On définit
$V(|i\rangle\otimes f_r)=v_if_r$, où $|1\rangle,\dots,|k\rangle$ est la base canonique de
$\mathbb C^k$ et où l'on identifie $\mathbb C^m$ à $eE$ au moyen des $f_r$. Alors
$E_{ij}v_lf_r=\delta_{jl}v_if_r$, c'est-à-dire $V^*E_{ij}V=|i\rangle\langle j|\otimes1_m$, et
$V^*\mathcal NV=M_k(\mathbb C)\otimes1_m$.

*Étape 5 : le commutant.* Écrivons $X\in M_k\otimes M_m$ comme une matrice par blocs
$X=\sum_{ij}|i\rangle\langle j|\otimes X_{ij}$. Pour $a,b$ fixés,
$(|a\rangle\langle b|\otimes1)X=\sum_j|a\rangle\langle j|\otimes X_{bj}$ et
$X(|a\rangle\langle b|\otimes1)=\sum_i|i\rangle\langle b|\otimes X_{ia}$. S'ils sont égaux, la
comparaison des blocs en position $(a,j)$ avec $j\neq b$ donne $X_{bj}=0$, et en position
$(a,b)$ elle donne $X_{bb}=X_{aa}$. Comme $a,b$ sont arbitraires, $X=1_k\otimes Y$ avec
$Y=X_{11}$ ; réciproquement, tout $1_k\otimes Y$ commute avec $M_k\otimes1_m$. $\square$

**Théorème 2.6 (structure des algèbres de von Neumann de dimension finie).** Soit
$\mathcal M$ une algèbre de von Neumann sur $\mathbb C^n$. Il existe $C\ge1$, des entiers
$k_c,m_c\ge1$ avec $\sum_ck_cm_c=n$, et un opérateur unitaire
$V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n$ tels que
$$
V^*\mathcal MV=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},\qquad
V^*\mathcal M'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),
$$
et $V^*Z(\mathcal M)V=\operatorname{span}\{V^*P_1V,\dots,V^*P_CV\}$, où $P_1,\dots,P_C$ sont
les projecteurs centraux minimaux de $\mathcal M$ (lemme 2.3) et $V^*P_cV$ est le projecteur sur
le $c$-ième facteur de la somme directe. Les couples $(k_c,m_c)$ sont uniques à l'ordre
près : pour *tout* opérateur unitaire $V$ ayant la première propriété, $V(\cdot)V^*$ envoie les
projecteurs sur les facteurs de la somme directe sur les $P_c$ (voir la démonstration), de sorte que
$k_c^2=\dim\mathcal MP_c$ et $m_c=\operatorname{rank}P_c/k_c$ sont déterminés par $\mathcal M$. En particulier $\mathcal M\cong\bigoplus_cM_{k_c}(\mathbb C)$
en tant que $*$-algèbres, et $\dim\mathcal M=\sum_ck_c^2$, $\dim\mathcal M'=\sum_cm_c^2$.

*Démonstration.* D'après le lemme 2.3 et la remarque qui le suit,
$\mathcal M=\bigoplus_c\mathcal M_c$, où $\mathcal M_c$ agit sur $E_c$, et chaque
$\mathcal M_c$ est un facteur (lemme 2.4). Le lemme 2.5 fournit des opérateurs unitaires
$V_c$ tels que $V_c^*\mathcal M_cV_c=M_{k_c}\otimes1_{m_c}$ ; on pose
$V=\bigoplus_cV_c$. Pour le commutant : $X\in\mathcal M'$ commute avec chaque $P_c$, donc
$X=\sum_cP_cXP_c$, et $P_cXP_c$ commute avec $\mathcal M_c$ ; réciproquement, toute somme de
cette forme commute avec $\mathcal M$. Le lemme 2.5 identifie $\mathcal M_c'$ sur
$E_c$, et la conjugaison par $V$ donne la formule pour $V^*\mathcal M'V$.

*Unicité.* Soit $V$ un opérateur unitaire quelconque tel que $V^*\mathcal MV=\mathcal B:=\bigoplus_{c=1}^{C'}M_{k_c}\otimes1_{m_c}$,
et soit $Q_c$ le projecteur sur le $c$-ième facteur de la somme directe. Un élément $\bigoplus_cX_c$ de $\mathcal B$ ($X_c\in M_{k_c}\otimes1_{m_c}$)
est central si et seulement si chaque $X_c$ commute avec $M_{k_c}\otimes1_{m_c}$, c'est-à-dire (lemme 2.5, étape 5) $X_c\in(M_{k_c}\otimes1_{m_c})\cap(1_{k_c}\otimes M_{m_c})=\mathbb C1$
(un élément $a\otimes1=1\otimes b$ a tous ses blocs diagonaux égaux à $b$ et tous ses blocs non diagonaux nuls, donc $a$ est un scalaire).
Par conséquent $Z(\mathcal B)=\operatorname{span}\{Q_1,\dots,Q_{C'}\}$, et $Z(\mathcal M)=VZ(\mathcal B)V^*$ est engendré par les projecteurs
$VQ_cV^*$, non nuls, deux à deux orthogonaux et de somme $1$. Chacun d'eux est un projecteur non nul de $Z(\mathcal M)$, donc (lemme 2.3)
la somme des $P_d$ indexés par une partie non vide $S_c$ ; les $S_c$ sont disjointes (orthogonalité) et recouvrent tous les $d$ (somme $1$), et
comme $C'=\dim Z(\mathcal M)=C$, chaque $S_c$ est un singleton. Ainsi, après renumérotation, $VQ_cV^*=P_c$, d'où
$\dim\mathcal MP_c=\dim(M_{k_c}\otimes1_{m_c})=k_c^2$ et $\operatorname{rank}P_c=k_cm_c$. $\square$

Les nombres $k_c$ décrivent l'algèbre elle-même. Les **multiplicités** $m_c$ décrivent la
manière dont elle est plongée dans $M_n(\mathbb C)$ : elles sont invisibles pour la
structure algébrique, mais elles déterminent la restriction à $\mathcal M$ de la trace de
$\mathbb C^n$,
$$
\operatorname{Tr}_{\mathbb C^n}(A)=\sum_cm_c\operatorname{Tr}(A_c),\qquad A=V\Big(\bigoplus_cA_c\otimes1_{m_c}\Big)V^*.
$$
C'est pourquoi la bibliothèque travaille avec une famille de traces (chapitre 3). De même, pour $A$
autoadjoint, la restriction de $A$ à $V(\mathbb C^{k_c}\otimes\mathbb C^{m_c})$ a les valeurs propres du bloc
$A_c$, chacune avec $m_c$ fois sa multiplicité dans $A_c$, puisque $A_c\otimes1_{m_c}$ est formé de $m_c$ copies
de $A_c$.

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Les $P_c$ sont les **secteurs de supersélection** ou
> **secteurs de symétrie** : aucune observable de $\mathcal M$ ne relie des secteurs
> différents, de sorte qu'une mesure de « quel secteur » est compatible avec toute autre
> mesure. Le centre est la partie *classique* du système. La multiplicité $m_c$ est une
> dégénérescence de toute observable de $\mathcal M$.

## 2.4 Exemple : un groupe de symétrie et son commutant

Soit un groupe fini $G$ agissant sur $\mathbb C^n$ par une représentation unitaire $\pi$.
Les opérateurs qui commutent avec la symétrie forment $\pi(G)'$, qui est une algèbre de von
Neumann ; son commutant $\pi(G)''$ est le sous-espace engendré par les $\pi(g)$
(théorème 2.2 appliqué à ce sous-espace, qui est une $*$-algèbre : il contient $\pi(e)=1$,
les produits $\pi(g)\pi(h)=\pi(gh)$ et les adjoints $\pi(g)^*=\pi(g^{-1})$).

**Proposition 2.7.** Appliquons le théorème 2.6 à $\mathcal M=\pi(G)''$. Alors
$V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$, où chaque $\rho_c$ est une représentation
unitaire irréductible de $G$ sur $\mathbb C^{k_c}$, et où $\rho_c,\rho_d$ ne sont pas
équivalentes pour $c\neq d$. Par conséquent
$$
\pi\cong\bigoplus_c\rho_c\otimes1_{m_c},\qquad
\pi(G)''\cong\bigoplus_cM_{\dim\rho_c},\qquad
\pi(G)'\cong\bigoplus_cM_{m_c}:
$$
les blocs des opérateurs symétriques $\pi(G)'$ ont pour tailles les multiplicités des
représentations irréductibles, et les rôles de $k_c$ et $m_c$ s'échangent entre une algèbre
et son commutant.

*Démonstration.* $\pi(g)\in\mathcal M$, donc $V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$
pour certains $\rho_c(g)\in M_{k_c}$, et $g\mapsto\rho_c(g)$ est une représentation
unitaire puisque $\pi$ en est une. Comme les $\pi(g)$ engendrent $\mathcal M$, les
$\rho_c(g)$ engendrent $M_{k_c}$, et plus précisément l'application
$\sum_ga_g\pi(g)\mapsto\big(\sum_ga_g\rho_c(g)\big)_c$ est surjective sur
$\bigoplus_cM_{k_c}$. *Irréductibilité :* un sous-espace de $\mathbb C^{k_c}$ invariant par
tous les $\rho_c(g)$ est invariant par le sous-espace $M_{k_c}$ qu'ils engendrent, donc
égal à $0$ ou à $\mathbb C^{k_c}$. *Non-équivalence :* si $\rho_d(g)=T\rho_c(g)T^{-1}$ pour
tout $g$ et un certain $T$ inversible, avec $c\neq d$, alors tout élément $(x_{c'})_{c'}$
de l'image vérifie $x_d=Tx_cT^{-1}$, ce qui exclut l'élément tel que $x_c=1$, $x_d=0$ — en
contradiction avec la surjectivité. Les deux derniers isomorphismes découlent du
théorème 2.6. $\square$

La décomposition d'une représentation en irréductibles est unique à équivalence près, de sorte
que les $m_c$ sont *les* multiplicités dans toute décomposition. *Démonstration.* Soient
$\rho,\rho'$ irréductibles et $T$ un opérateur d'entrelacement, $T\rho(g)=\rho'(g)T$ pour tout $g$.
Le noyau et l'image de $T$ sont invariants, donc $T=0$ ou $T$ est inversible (lemme de Schur) ; si
$\rho'=\rho$, il en va de même de $T-\lambda1$ pour une valeur propre $\lambda$ de $T$, donc
$T=\lambda1$. Ainsi les opérateurs d'entrelacement de $\rho$ vers $\rho'$ forment un espace de
dimension $1$ si $\rho\cong\rho'$ (pour un opérateur d'entrelacement inversible $S$, $T\mapsto S^{-1}T$
les identifie aux opérateurs d'entrelacement de $\rho$ vers elle-même, qui forment $\mathbb C1$) et $0$ sinon
(un opérateur non nul serait inversible, donc une équivalence). Si $\pi\cong\bigoplus_j\sigma_j$ avec des $\sigma_j$
irréductibles (répétitions permises), un opérateur d'entrelacement d'une $\rho$ irréductible vers
$\bigoplus_j\sigma_j$ est une famille $(T_j)_j$ d'opérateurs d'entrelacement de $\rho$ vers $\sigma_j$, et la
composition avec l'équivalence ne change pas les dimensions ; les opérateurs d'entrelacement de $\rho$ vers
$\pi$ forment donc un espace de dimension $\#\{j:\sigma_j\cong\rho\}$, nombre qui ne dépend que de
$\pi$ et de $\rho$. Pour la décomposition de la proposition 2.7, dont les $\rho_c$ sont deux à deux
non équivalentes, ce nombre vaut $m_c$ si $\rho\cong\rho_c$ et $0$ si $\rho$ n'est équivalente à
aucune $\rho_c$. $\square$

Prenons $G=S_3$ agissant par permutation des trois facteurs de $(\mathbb C^2)^{\otimes3}$,
c'est-à-dire trois qubits ou trois spins $\tfrac12$ : $\pi(\sigma)|x\rangle=|\sigma\cdot x\rangle$ pour $x\in\{0,1\}^3$, où
$(\sigma\cdot x)_{\sigma(i)}=x_i$. Soit $F=\{v:\pi(g)v=v\text{ pour tout }g\}$ le **sous-espace symétrique**. La
décomposition est la suivante : $F$ (dimension $4$) porte la représentation triviale avec multiplicité $4$ ; son
supplémentaire orthogonal (dimension $4$) porte la **représentation standard** de dimension deux (la
représentation par permutations de $S_3$ sur les vecteurs de $\mathbb C^3$ de somme des coordonnées nulle) avec
multiplicité $2$ ; la représentation signature n'apparaît pas. Ainsi
$\pi(S_3)''\cong M_1\oplus M_2$ (dimension $5$) et $\pi(S_3)'\cong M_4\oplus M_2$
(dimension $20$).

*Démonstration.* Écrivons $V^*\pi(S_3)''V=\bigoplus_cM_{k_c}\otimes1_{m_c}$ (théorème 2.6), de sorte que
$\sum_ck_cm_c=8$ et $\sum_cm_c^2=\dim\pi(S_3)'$.

*Étape 1 : $\dim\pi(S_3)'=20$.* Comme $\pi(g)^*|g\cdot x\rangle=|x\rangle$, la matrice
$\pi(g)X\pi(g)^*$ a le coefficient $X_{xy}$ en position $(g\cdot x,g\cdot y)$. Par conséquent $X\in\pi(S_3)'$, c'est-à-dire
$\pi(g)X\pi(g)^*=X$ pour tout $g$, si et seulement si $(x,y)\mapsto X_{xy}$ est constante sur les orbites de $S_3$ agissant
sur les couples $(x,y)$ par $g\cdot(x,y)=(g\cdot x,g\cdot y)$ ; les matrices indicatrices des orbites forment une base de
$\pi(S_3)'$, dont la dimension est donc le nombre d'orbites. Un couple $(x,y)$ est un triplet
$\big((x_1,y_1),(x_2,y_2),(x_3,y_3)\big)$ d'éléments de l'ensemble à quatre éléments $\{0,1\}^2$, et $S_3$
permute les trois composantes ; une orbite est donc un multiensemble de taille $3$ dans un ensemble à quatre éléments : $4$ de la forme
$\{a,a,a\}$, $4\cdot3=12$ de la forme $\{a,a,b\}$ avec $b\neq a$, et $4$ de la forme $\{a,b,c\}$ à
éléments distincts, $20$ en tout.

*Étape 2 : $k=(1,2)$, $m=(4,2)$.* Comme $k_c\ge1$, $\sum_cm_c\le\sum_ck_cm_c=8$. Soit $M=\max_cm_c$.
Alors $20=\sum_cm_c^2\le M\sum_cm_c\le8M$ et $M^2\le20$, donc $M\in\{3,4\}$. Si $M=4$, les autres
$m_c$ ont une somme des carrés égale à $4$ : ce sont $(2)$ ou $(1,1,1,1)$. Si $M=3$, les autres $m_c$ ont une somme
des carrés égale à $11$ et une somme au plus égale à $5$ ; s'ils étaient tous $\le2$, leur somme des carrés serait au plus
$2\cdot5=10$, donc l'un d'eux vaut $3$ et les autres ont une somme des carrés égale à $2$ : ce sont $(1,1)$.
Ainsi $m$ est, à l'ordre près, $(4,2)$, $(4,1,1,1,1)$ ou $(3,3,1,1)$. Dans les deux derniers cas,
$\sum_cm_c=8=\sum_ck_cm_c$ impose $k_c=1$ pour tout $c$, et alors $\pi(S_3)''\cong\bigoplus_cM_1$ serait
commutative, alors que les transpositions $(12)$ et $(23)$ ne commutent pas :
$\pi((12))\pi((23))|001\rangle=|100\rangle$ mais $\pi((23))\pi((12))|001\rangle=|010\rangle$. Donc
$m=(4,2)$ (en numérotant les facteurs de sorte que $m_1=4$), et $4k_1+2k_2=8$ avec $k_1,k_2\ge1$ donne $k=(1,2)$. En particulier
$\dim\pi(S_3)''=\sum_ck_c^2=5$.

*Étape 3 : la représentation triviale et la signature.* $\rho_2$ ne fixe aucun vecteur non nul $u$ : la droite
$\mathbb Cu$ serait un sous-espace invariant de $\mathbb C^2$, contrairement à l'irréductibilité
(proposition 2.7). Écrivons $V^*v=v_1\oplus v_2$ avec $v_1\in\mathbb C^1\otimes\mathbb C^4$ et
$v_2=\sum_{r=1}^2u_r\otimes f_r\in\mathbb C^2\otimes\mathbb C^2$ ($f_1,f_2$ la base canonique). Alors
$v\in F$ si et seulement si $\rho_1(g)v_1=v_1$ et $\rho_2(g)u_r=u_r$ pour tous $g$ et $r$, donc tout $v\in F$ vérifie
$v_2=0$. Le vecteur non nul $|000\rangle$ appartient à $F$, donc le scalaire $\rho_1(g)$ fixe un vecteur non nul :
$\rho_1$ est triviale, $F=V(\mathbb C^1\otimes\mathbb C^4)$ est de dimension $4$, et son supplémentaire
orthogonal est $V(\mathbb C^2\otimes\mathbb C^2)$. La représentation signature n'est équivalente ni à
$\rho_1$ (car $\operatorname{sgn}(12)=-1$) ni à $\rho_2$ (les dimensions diffèrent) ; d'après l'unicité des
multiplicités démontrée plus haut, elle apparaît donc $0$ fois.

*Étape 4 : $\rho_2$ est la représentation standard.* Les vecteurs $e_1=|100\rangle$, $e_2=|010\rangle$,
$e_3=|001\rangle$ sont permutés par les $\pi(\sigma)$ comme $\sigma$ permute les positions, donc
$W=\{\sum_ix_ie_i:\sum_ix_i=0\}$ est invariant et $\pi$ y agit par la représentation standard.
$W\perp F$ : pour $v\in F$ et $w=\sum_ix_ie_i\in W$, $\langle v,w\rangle=\langle\pi(g)v,w\rangle=\langle v,\pi(g)^{-1}w\rangle$
pour tout $g$, donc $\langle v,w\rangle=\langle v,\frac16\sum_g\pi(g)w\rangle$, et
$\frac16\sum_g\pi(g)w=\frac13\big(\sum_ix_i\big)(e_1+e_2+e_3)=0$ (chaque $e_i$ est envoyé sur chaque $e_j$ par
exactement deux des six permutations). Par conséquent $V^*W\subseteq\mathbb C^2\otimes\mathbb C^2$. Les applications
$T_r:\mathbb C^2\otimes\mathbb C^2\to\mathbb C^2$, $u\otimes f\mapsto\langle f_r,f\rangle u$, vérifient
$T_r(\rho_2(g)\otimes1)=\rho_2(g)T_r$. Choisissons $r$ tel que $T_rV^*|_W\neq0$ (si les deux s'annulaient, $V^*W$ serait
nul). Alors $T_rV^*|_W$ est un opérateur d'entrelacement non nul de $\pi|_W$ vers $\rho_2$ ; son image est un sous-espace
invariant non nul de $\mathbb C^2$, donc égale à $\mathbb C^2$, et comme $\dim W=2$, il est bijectif. Ainsi $\rho_2$ est
équivalente à la représentation standard (qui, au passage, est irréductible). $\square$

Le code ci-dessous confirme numériquement les deux dimensions : le commutant d'une partie $S$ est le
noyau de l'application linéaire $X\mapsto(AX-XA)_{A\in S}$, dont la dimension est calculée ici en virgule
flottante, en comptant les valeurs singulières supérieures à une tolérance. C'est une vérification des valeurs
qui viennent d'être démontrées, non une partie de la démonstration.

```python
import itertools
import torch

torch.manual_seed(0)
dt = torch.complex128

def permutation_operator(perm, n_qubits=3):
    """Unitary on (C^2)^{(x)3} that permutes the tensor factors."""
    dim = 2 ** n_qubits
    P = torch.zeros(dim, dim, dtype=dt)
    for bits in itertools.product([0, 1], repeat=n_qubits):
        src = int(''.join(map(str, bits)), 2)
        dst = int(''.join(str(bits[perm[i]]) for i in range(n_qubits)), 2)
        P[dst, src] = 1
    return P

group = [permutation_operator(p) for p in itertools.permutations(range(3))]

def commutant_basis(ops, tol=1e-10):
    """Orthonormal basis (Hilbert-Schmidt) of {X : X A = A X for all A in ops}."""
    d = ops[0].shape[0]
    eye = torch.eye(d, dtype=dt)
    # vec(AX - XA) = (A (x) 1 - 1 (x) A^T) vec(X) for row-major vec
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    _, s, Vh = torch.linalg.svd(L)
    null = Vh[(s > tol).sum():].conj().resolve_conj()   # rows span the null space
    return [v.reshape(d, d) for v in null]

def span_dim(ops, tol=1e-10):
    return int((torch.linalg.svdvals(torch.stack([A.reshape(-1) for A in ops])) > tol).sum())

comm = commutant_basis(group)
print(len(comm), span_dim(group), len(commutant_basis(comm)))
assert len(comm) == 20                                   # Step 1: pi(S3)'  = M_4 (+) M_2
assert span_dim(group) == 5                              # Step 2: pi(S3)'' = M_1 (+) M_2
assert len(commutant_basis(comm)) == 5                   # double commutant theorem
```

Les multiplicités se lisent sur le spectre. Tout $K\in\pi(S_3)'$ autoadjoint est de la forme
$V^*KV=(1_1\otimes\kappa_1)\oplus(1_2\otimes\kappa_2)$ avec $\kappa_1\in M_4$, $\kappa_2\in M_2$ autoadjoints
(théorème 2.6). Sur $F$, ses valeurs propres sont celles de $\kappa_1$ ; sur le supplémentaire de $F$, ce sont celles
de $\kappa_2$, chacune répétée deux fois, puisque $1_2\otimes\kappa_2$ est formé de deux copies de $\kappa_2$. Ainsi *toute* valeur propre de
$K$ sur le supplémentaire de $F$ est de multiplicité paire. Pour un $K$ aléatoire, on a davantage. La moyenne sur le
groupe, $E(X)=\frac1{|G|}\sum_g\pi(g)X\pi(g)^*$, envoie $M_8(\mathbb C)$ sur le commutant et y est
l'identité (exercice 3) ; pour une matrice de Ginibre $Z$ (coefficients gaussiens complexes standard
indépendants), $K=E(Z+Z^*)$ a donc une densité sur l'espace vectoriel réel $\pi(S_3)'_h$ des éléments autoadjoints
(théorème A.15 et la remarque qui le suit). D'après la proposition 4.7(a) (une référence anticipée ; sa démonstration utilise
le théorème 2.6 et l'annexe A, pas ce paragraphe), presque sûrement les six
valeurs propres de $\kappa_1$ et $\kappa_2$ sont deux à deux distinctes : $K$ a alors $4$ valeurs propres simples sur $F$ et
$2$ valeurs propres doubles sur son supplémentaire. Le code vérifie les paires dégénérées (valables pour tout $K$)
et, pour le $K$ tiré, que les deux paires et les quatre valeurs propres sur $F$ sont distinctes, en utilisant que
$P_F=\frac1{|G|}\sum_g\pi(g)$ est le projecteur orthogonal sur $F$ ($\pi(h)P_F=P_F$ pour tout $h$, donc son
image est contenue dans $F$ ; $P_Fv=v$ pour $v\in F$ ; et $P_F^*=P_F$, par le changement de variable $g\to g^{-1}$).

```python
def group_average(X):
    return sum(g @ X @ g.conj().T for g in group) / len(group)

H0 = torch.randn(8, 8, dtype=dt)                        # Ginibre: Re, Im independent N(0, 1/2)
Hsym = group_average(H0 + H0.conj().T)                  # K = E(Z + Z^*), a random symmetric Hamiltonian
assert all(torch.allclose(g @ Hsym, Hsym @ g) for g in group)

P_sym = sum(group) / len(group)                         # P_F, the projection onto the symmetric subspace
w, Q = torch.linalg.eigh(P_sym)
Q_rest = Q[:, w.real < 0.5]                             # orthonormal basis of the complement (dim 4)
e_rest = torch.linalg.eigvalsh(Q_rest.conj().T @ Hsym @ Q_rest)
print(e_rest)
assert Q_rest.shape[1] == 4
assert torch.allclose(e_rest[0], e_rest[1]) and torch.allclose(e_rest[2], e_rest[3])   # even multiplicities
e_sym = torch.linalg.eigvalsh(Q[:, w.real > 0.5].conj().T @ Hsym @ Q[:, w.real > 0.5])   # eigenvalues on F
assert e_sym.shape[0] == 4 and e_rest[2] - e_rest[1] > 1e-6                           # this sample: generic
assert (e_sym[1:] - e_sym[:-1]).min() > 1e-6                                           # 4 simple eigenvalues on F
assert (e_sym[:, None] - e_rest[None, :]).abs().min() > 1e-6                           # all six distinct
```

## 2.5 L'algèbre dans la bibliothèque

La bibliothèque stocke l'algèbre abstraite $\bigoplus_cM_{k_c}(\mathbb C)$ du théorème 2.6
sur $\bigoplus_c\mathbb C^{k_c}$, sans les multiplicités : un opérateur est la liste de ses
blocs $A_c$, complétés par des zéros jusqu'à une taille commune. Le commutant de l'exemple
est `TypeIAlgebra([4, 2], [4, 2])`.

```python
from torch_vn_algebra import TypeIAlgebra, channels

alg = TypeIAlgebra([4, 2], [4, 2], precision='double')    # M_4 (+) M_2
print(alg.k_factors, alg.C)

A = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
print(A.matrix.shape)              # (1, 2, 4, 4): block 1 occupies the top-left 2 x 2 corner
assert A.matrix[0, 1, 2:, :].abs().max() == 0 and A.matrix[0, 1, :, 2:].abs().max() == 0
```

Le centre est formé des éléments $\sum_cz_cP_c$ ; `alg.central` les construit, et ils
commutent avec tout :

```python
def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

P0, P1 = alg.central([1.0, 0.0]), alg.central([0.0, 1.0])   # minimal central projections
Z = alg.central([2.0, -1.0])
assert P0.is_projection and P1.is_projection
assert close(P0 + P1, alg.identity())
assert close(Z @ A, A @ Z)
assert close(A @ P0 + A @ P1, A)                         # A = sum_c A P_c
```

## 2.6 Espérances conditionnelles

Une **espérance conditionnelle** de $\mathcal M$ sur une sous-algèbre de von Neumann
$\mathcal N\subseteq\mathcal M$ (une $*$-algèbre contenue dans $\mathcal M$, de même unité) est une
application linéaire $E:\mathcal M\to\mathcal N$ positive ($A\ge0\Rightarrow E(A)\ge0$), telle que $E(B)=B$ pour $B\in\mathcal N$, et qui est
un morphisme de $\mathcal N$-bimodules : $E(B_1AB_2)=B_1E(A)B_2$ pour
$B_1,B_2\in\mathcal N$.

**Proposition 2.8.** Pour $\mathcal M=\bigoplus_cM_{k_c}$ sur $\bigoplus_c\mathbb C^{k_c}$, avec $P_c$ le projecteur sur
le $c$-ième bloc, et $\mathcal N=Z(\mathcal M)=\operatorname{span}\{P_1,\dots,P_C\}$ (exercice 1),
$$
E(A)=\sum_c\frac{\operatorname{Tr}A_c}{k_c}\,P_c
$$
est une espérance conditionnelle. Elle préserve toute trace de la forme
$\sum_cw_c\operatorname{Tr}A_c$.

*Démonstration.* La linéarité est claire ; $A\ge0$ entraîne $\operatorname{Tr}A_c=\sum_i\langle e_i,Ae_i\rangle\ge0$
(somme sur les vecteurs $e_i$ de la base canonique du $c$-ième bloc), donc $E(A)\ge0$. Pour $B=\sum_cb_cP_c$, $E(B)=\sum_cb_ck_cP_c/k_c=B$. Pour
$B_i=\sum_cb^{(i)}_cP_c$, $\operatorname{Tr}(B_1AB_2)_c=b^{(1)}_cb^{(2)}_c\operatorname{Tr}A_c$,
ce qui donne la propriété de bimodule. Enfin
$\operatorname{Tr}(E(A)_c)=\operatorname{Tr}A_c$ pour tout $c$. $\square$

Dans la bibliothèque, cette application est le canal `channels.center_expectation(alg)`.
Appliquée aux opérateurs, elle en donne la partie classique (moyennée sur chaque secteur) :

```python
E = channels.center_expectation(alg)
EA = E(A)
assert close(E(EA), EA)                                  # idempotent
assert close(E(Z @ A @ P1), Z @ EA @ P1)                 # bimodule property
A2 = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
assert close(EA @ A2, A2 @ EA)                           # E(A) commutes with a fresh element
for c, k in enumerate(alg.k_factors):                    # E(A)_c = (Tr A_c / k_c) 1, so Tr E(A)_c = Tr A_c
    trA = A.matrix[0, c].diagonal().sum()
    assert close(EA.matrix[0, c, :k, :k], trA / k * torch.eye(k, dtype=dt))
    assert abs(EA.matrix[0, c].diagonal().sum() - trA) < 1e-10
```

La moyenne sur le groupe de la section 2.4 est une espérance conditionnelle sur le
commutant $\pi(G)'$ dans $M_n(\mathbb C)$, et elle préserve $\operatorname{Tr}$
(exercice 3). De même, dans l'algèbre de matrices complète $M_D(\mathbb C)$, $D=\sum_ck_c$, qui contient
$\mathcal M=\bigoplus_cM_{k_c}$ comme ses matrices diagonales par blocs, l'application $E_{\mathcal M}(X)=\sum_cP_cXP_c$
est une espérance conditionnelle sur $\mathcal M$ qui préserve $\operatorname{Tr}$. *Démonstration.* $P_cXP_c$
s'annule hors du $c$-ième bloc diagonal, donc $E_{\mathcal M}(X)\in\mathcal M$. Pour $X\ge0$,
$\langle v,P_cXP_cv\rangle=\langle P_cv,XP_cv\rangle\ge0$. Un élément $B\in\mathcal M$ est diagonal par blocs, donc
$BP_c=P_cB$ et $B=\sum_cP_cBP_c$ ; par conséquent $E_{\mathcal M}(B)=B$ et
$E_{\mathcal M}(B_1XB_2)=\sum_cB_1P_cXP_cB_2=B_1E_{\mathcal M}(X)B_2$ pour $B_1,B_2\in\mathcal M$. Enfin
$\sum_c\operatorname{Tr}(P_cXP_c)=\sum_c\operatorname{Tr}(P_cX)=\operatorname{Tr}X$. $\square$

Pour un **état** $\omega$ sur $\mathcal M$ (une forme linéaire telle que $\omega(A^*A)\ge0$ pour tout $A$ et $\omega(1)=1$ ;
les états sont étudiés au chapitre 3), la proposition 2.8 et la linéarité donnent
$\omega(E(A))=\sum_c\omega(P_c)\operatorname{Tr}(A_c)/k_c$, et $\omega\circ E=\omega$ sur le centre, puisque $E(z)=z$
sur celui-ci. (Les espérances conditionnelles ne sont pas uniques en général : sur
$\mathbb C1$, chaque état $\varphi$ en fournit une, $A\mapsto\varphi(A)1$ ; sur $M_2(\mathbb C)$, les états
$A\mapsto A_{11}$ et $A\mapsto A_{22}$ en fournissent deux différentes.) Ce sont toutes des
instances d'une même idée : projeter sur les opérateurs qui « ne voient pas » une certaine
structure.

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). L'application $E_{\mathcal M}$
> efface les cohérences entre secteurs de supersélection, ce qui se produit lorsque l'on mesure la
> charge (le secteur) sans enregistrer le résultat. L'espérance conditionnelle $E$ sur le centre va
> plus loin : elle ne garde que l'information classique, le poids $\omega(P_c)$ de chaque secteur, et
> remplace l'état à l'intérieur de chaque secteur par l'état maximalement mélangé (la trace
> normalisée $\operatorname{Tr}(A_c)/k_c$).

## 2.7 Points de comparaison

Les exemples de la section 2.4 peuvent être vérifiés avec n'importe quelle bibliothèque
d'algèbre linéaire ; les logiciels de théorie des représentations (GAP, SageMath)
fournissent la décomposition $\pi\cong\bigoplus\rho_\lambda\otimes1_{m_\lambda}$ de façon
symbolique. Les codes de diagonalisation exacte utilisés en physique (QuSpin, les modules
de symétrie de QuTiP) exploitent la même structure par blocs pour diagonaliser les
hamiltoniens secteur par secteur. La bibliothèque prend cette structure comme point de
départ : une fois les blocs connus, tous les calculs ultérieurs — fonctions, traces,
états, canaux, dynamique — restent à l'intérieur de l'algèbre.

## Exercices

1. Montrer que le centre de $M_k(\mathbb C)$ est $\mathbb C1$ (utiliser les unités
   matricielles $|i\rangle\langle j|$), et que $Z(\bigoplus_cM_{k_c})$ est de dimension $C$.
2. Soit $S$ la permutation circulaire de $\mathbb C^3$, $Se_i=e_{i+1\bmod3}$. Déterminer
   $\{S\}'$ et $\{S\}''$, et vérifier numériquement les dimensions.
3. Pour une représentation unitaire $\pi$ d'un groupe fini, montrer que
   $E(X)=\frac1{|G|}\sum_g\pi(g)X\pi(g)^*$ est une espérance conditionnelle de
   $M_n(\mathbb C)$ sur $\pi(G)'$ qui préserve $\operatorname{Tr}$.

<details>
<summary>Solutions</summary>

1. Si $X$ commute avec tous les $|i\rangle\langle j|$, alors
   $X|i\rangle\langle j|=|i\rangle\langle j|X$ ; en appliquant les deux membres à $e_j$ et
   en comparant, on obtient $Xe_i=X_{jj}e_i$ pour tous $i,j$, donc $X=X_{11}1$. Dans une
   somme directe, un élément $\bigoplus X_c$ est central si et seulement si chaque $X_c$
   l'est, donc le centre est $\{\bigoplus_cz_c1_{k_c}\}$, de dimension $C$.

2. $S$ est unitaire, de valeurs propres $1,\omega,\omega^2$, $\omega=e^{2\pi i/3}$, toutes
   simples : les vecteurs $u_j=\sum_i\omega^{-ij}e_i$ ($j=0,1,2$) vérifient $Su_j=\omega^ju_j$ et sont orthogonaux.
   Un opérateur qui commute avec $S$ préserve ses sous-espaces propres ; dans la
   base propre, il est donc diagonal, et toute matrice diagonale commute avec $S$ : $\{S\}'\cong\mathbb C^3$.
   Il contient $1,S,S^2$, qui dans la base propre s'écrivent $\operatorname{diag}(\lambda^0)$, $\operatorname{diag}(\lambda)$,
   $\operatorname{diag}(\lambda^2)$ pour $\lambda=(1,\omega,\omega^2)$, linéairement indépendants (Vandermonde : les
   valeurs propres sont distinctes) ; donc $\{S\}'=\operatorname{span}\{1,S,S^2\}$, les matrices circulantes. Alors
   $\{S\}''$ est lui aussi égal à $\{S\}'$, car dans la base propre de
   $S$ une matrice qui commute avec toutes les matrices diagonales est diagonale (comparer
   les coefficients de $XD$ et de $DX$ pour une matrice $D$ à coefficients diagonaux
   distincts).

```python
S = torch.roll(torch.eye(3, dtype=dt), 1, dims=0)
assert len(commutant_basis([S])) == 3
assert len(commutant_basis(commutant_basis([S]))) == 3
assert span_dim([torch.eye(3, dtype=dt), S, S @ S]) == 3
```

3. $E$ est linéaire et positive (somme d'applications $X\mapsto UXU^*$). Pour $h\in G$,
   $\pi(h)E(X)\pi(h)^*=\frac1{|G|}\sum_g\pi(hg)X\pi(hg)^*=E(X)$, donc $E(X)\in\pi(G)'$. Si
   $X\in\pi(G)'$, chaque terme vaut $X$, donc $E(X)=X$. Pour $B_1,B_2\in\pi(G)'$,
   $\pi(g)B_1XB_2\pi(g)^*=B_1\pi(g)X\pi(g)^*B_2$, ce qui donne la propriété de bimodule.
   Enfin $\operatorname{Tr}(\pi(g)X\pi(g)^*)=\operatorname{Tr}X$ pour tout $g$, donc
   $\operatorname{Tr}E(X)=\operatorname{Tr}X$.

</details>

**Suite :** le [chapitre 3](03_traces_and_states.md) classifie les traces sur
$\bigoplus_cM_{k_c}$, introduit les états et les densités, et démontre les inégalités
d'entropie et de traces utilisées par la suite.
