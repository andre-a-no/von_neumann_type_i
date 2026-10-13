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

En dimension infinie, $\mathcal M=\mathcal M''$ est une condition topologique (fermeture
pour la topologie faible des opérateurs). En dimension finie, elle est automatique :

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
$\mu\in(\operatorname{spec}x\cup\{0\})\setminus\{\lambda\}$ (démonstration du
théorème 1.5(a)). Comme $p(0)=0$ et $x=ex=xe$, on a $q\in\mathcal N$ et $q=eqe$, donc
l'image de $q$ est contenue dans celle de $e$. Par minimalité, $q=e$ ou $q=0$ ; ainsi $x$ a
au plus une valeur propre non nulle et le sous-espace propre associé est l'image de $e$,
donc $x=\lambda e$. Un élément quelconque de $e\mathcal Ne$ s'écrit $x_1+ix_2$ avec
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
près : les $P_c$ sont uniques (lemme 2.3), $k_c^2=\dim\mathcal MP_c$ et
$m_c=\operatorname{rank}P_c/k_c$. En particulier $\mathcal M\cong\bigoplus_cM_{k_c}(\mathbb C)$
en tant que $*$-algèbres, et $\dim\mathcal M=\sum_ck_c^2$, $\dim\mathcal M'=\sum_cm_c^2$.

*Démonstration.* D'après le lemme 2.3 et la remarque qui le suit,
$\mathcal M=\bigoplus_c\mathcal M_c$, où $\mathcal M_c$ agit sur $E_c$, et chaque
$\mathcal M_c$ est un facteur (lemme 2.4). Le lemme 2.5 fournit des opérateurs unitaires
$V_c$ tels que $V_c^*\mathcal M_cV_c=M_{k_c}\otimes1_{m_c}$ ; on pose
$V=\bigoplus_cV_c$. Pour le commutant : $X\in\mathcal M'$ commute avec chaque $P_c$, donc
$X=\sum_cP_cXP_c$, et $P_cXP_c$ commute avec $\mathcal M_c$ ; réciproquement, toute somme de
cette forme commute avec $\mathcal M$. Le lemme 2.5 identifie $\mathcal M_c'$ sur
$E_c$. $\square$

Les nombres $k_c$ décrivent l'algèbre elle-même. Les **multiplicités** $m_c$ décrivent la
manière dont elle est plongée dans $M_n(\mathbb C)$ : elles sont invisibles pour la
structure algébrique, mais elles déterminent la restriction à $\mathcal M$ de la trace de
$\mathbb C^n$,
$$
\operatorname{Tr}_{\mathbb C^n}(A)=\sum_cm_c\operatorname{Tr}(A_c),\qquad A=V\Big(\bigoplus_cA_c\otimes1_{m_c}\Big)V^*.
$$
C'est pourquoi la bibliothèque travaille avec une famille de traces (chapitre 3).

> **Interprétation physique.** Les $P_c$ sont les **secteurs de supersélection** ou
> **secteurs de symétrie** : aucune observable de $\mathcal M$ ne relie des secteurs
> différents, de sorte qu'une mesure de « quel secteur » est compatible avec toute autre
> mesure. Le centre est la partie *classique* du système. La multiplicité $m_c$ est une
> dégénérescence : chaque valeur propre du bloc $A_c$ apparaît dans $\mathbb C^n$ avec
> $m_c$ fois sa multiplicité dans $A_c$.

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
dimension $1$ si $\rho\cong\rho'$ et $0$ sinon. Si $\pi\cong\bigoplus_j\sigma_j$ avec des $\sigma_j$
irréductibles (répétitions permises), les opérateurs d'entrelacement d'une $\rho$ irréductible vers
$\pi$ forment donc un espace de dimension $\#\{j:\sigma_j\cong\rho\}$, nombre qui ne dépend que de
$\pi$ et de $\rho$. Pour la décomposition de la proposition 2.7, dont les $\rho_c$ sont deux à deux
non équivalentes, ce nombre vaut $m_c$ si $\rho\cong\rho_c$ et $0$ si $\rho$ n'est équivalente à
aucune $\rho_c$. $\square$

Prenons $G=S_3$ agissant par permutation des trois facteurs de $(\mathbb C^2)^{\otimes3}$,
c'est-à-dire trois qubits ou trois spins $\tfrac12$. La décomposition est la suivante : le
sous-espace symétrique (spin total $\tfrac32$, dimension $4$) porte la représentation triviale
avec multiplicité $4$ ; son supplémentaire (dimension $4$) porte la représentation standard de
dimension deux avec multiplicité $2$ ; la représentation signature n'apparaît pas. Ainsi
$\pi(S_3)''\cong M_1\oplus M_2$ (dimension $5$) et $\pi(S_3)'\cong M_4\oplus M_2$
(dimension $20$).

*Démonstration.* Écrivons $\pi(S_3)''\cong\bigoplus_cM_{k_c}\otimes1_{m_c}$ (théorème 2.6). Alors
$\sum_ck_cm_c=8$, $\sum_ck_c^2=\dim\pi(S_3)''=5$ et $\sum_cm_c^2=\dim\pi(S_3)'=20$ ; le code
ci-dessous calcule ces deux dimensions par algèbre linéaire exacte (le commutant d'une partie $S$
est le noyau de l'application linéaire $X\mapsto(AX-XA)_{A\in S}$). Le nombre $20$ se compte
aussi à la main : $X$ commute avec tous les $\pi(g)$ si et seulement si ses coefficients
$X_{xy}$, $x,y\in\{0,1\}^3$, sont constants sur les orbites de $S_3$ agissant sur les couples
$(x,y)$, c'est-à-dire sur les triplets d'éléments de $\{0,1\}^2$, et il y a $\binom63=20$ telles
orbites (multiensembles de taille $3$ dans un ensemble à $4$ éléments). De $\sum_ck_c^2=5$ on tire
$k=(1,2)$ ou $k=(1,1,1,1,1)$. Dans le premier cas, $m_1+2m_2=8$ et $m_1^2+m_2^2=20$ donnent
$m=(4,2)$. Dans le second, cinq entiers strictement positifs de somme $8$ et de somme des carrés
$20$ sont nécessairement $m=(4,1,1,1,1)$ ; mais alors $\pi(S_3)''$ serait commutative, alors que
les transpositions $(12)$ et $(23)$ ne commutent pas : $\pi((12))\pi((23))|001\rangle=|100\rangle$
mais $\pi((23))\pi((12))|001\rangle=|010\rangle$. Donc $k=(1,2)$ et $m=(4,2)$. La représentation
de dimension un $\rho_1$ est la triviale ou la signature. Ce n'est pas la signature : chaque
vecteur de base $|b_1b_2b_3\rangle$ a deux bits égaux, donc il est fixé par une transposition $t$,
et $A=\sum_g\operatorname{sgn}(g)\pi(g)$, qui vérifie $A\pi(t)=-A$, donne
$A|b_1b_2b_3\rangle=A\pi(t)|b_1b_2b_3\rangle=-A|b_1b_2b_3\rangle$ ; ainsi $A=0$, et la
représentation signature n'apparaît pas. Donc $\rho_1$ est triviale, ses $4$ copies engendrent le
sous-espace symétrique, et $\rho_2$ est une représentation irréductible de dimension deux
(proposition 2.7) : la représentation standard. $\square$

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
assert len(comm) == 20                                   # pi(S3)'  = M_4 (+) M_2
assert span_dim(group) == 5                              # pi(S3)'' = M_1 (+) M_2
assert len(commutant_basis(comm)) == 5                   # double commutant theorem
```

Les multiplicités se lisent sur le spectre. La moyenne sur le groupe,
$X\mapsto\frac1{|G|}\sum_g\pi(g)X\pi(g)^*$, envoie tout opérateur dans le commutant
(exercice 3). Un élément autoadjoint aléatoire de $\pi(S_3)'$ possède donc (presque
sûrement) $4$ valeurs propres simples sur le sous-espace symétrique et $2$ valeurs propres
doubles sur son supplémentaire, car le bloc de $\pi(S_3)'$ y est $1_2\otimes M_2$ ; le code
vérifie les paires dégénérées.

```python
def group_average(X):
    return sum(g @ X @ g.conj().T for g in group) / len(group)

H0 = torch.randn(8, 8, dtype=dt)
Hsym = group_average(H0 + H0.conj().T)                  # a random symmetric Hamiltonian
assert all(torch.allclose(g @ Hsym, Hsym @ g) for g in group)

P_sym = sum(group) / len(group)                         # projection onto the symmetric subspace
w, Q = torch.linalg.eigh(P_sym)
Q_rest = Q[:, w.real < 0.5]                             # orthonormal basis of the complement (dim 4)
e_rest = torch.linalg.eigvalsh(Q_rest.conj().T @ Hsym @ Q_rest)
print(e_rest)
assert Q_rest.shape[1] == 4
assert torch.allclose(e_rest[0], e_rest[1]) and torch.allclose(e_rest[2], e_rest[3])
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
$\mathcal N\subseteq\mathcal M$ est une application linéaire $E:\mathcal M\to\mathcal N$
positive ($A\ge0\Rightarrow E(A)\ge0$), telle que $E(B)=B$ pour $B\in\mathcal N$, et qui est
un morphisme de $\mathcal N$-bimodules : $E(B_1AB_2)=B_1E(A)B_2$ pour
$B_1,B_2\in\mathcal N$.

**Proposition 2.8.** Pour $\mathcal M=\bigoplus_cM_{k_c}$ et $\mathcal N=Z(\mathcal M)$,
$$
E(A)=\sum_c\frac{\operatorname{Tr}A_c}{k_c}\,P_c
$$
est une espérance conditionnelle. Elle préserve toute trace de la forme
$\sum_cw_c\operatorname{Tr}A_c$.

*Démonstration.* La linéarité est claire ; $A\ge0$ entraîne $\operatorname{Tr}A_c\ge0$,
donc $E(A)\ge0$. Pour $B=\sum_cb_cP_c$, $E(B)=\sum_cb_ck_cP_c/k_c=B$. Pour
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
for c in range(alg.C):                                   # Tr E(A)_c = Tr A_c
    assert abs(EA.matrix[0, c].diagonal().sum() - A.matrix[0, c].diagonal().sum()) < 1e-10
```

La moyenne sur le groupe de la section 2.4 est une espérance conditionnelle sur le
commutant $\pi(G)'$ dans $M_n(\mathbb C)$, et elle préserve $\operatorname{Tr}$
(exercice 3). (Les espérances conditionnelles ne sont pas uniques en général : sur
$\mathbb C1$, chaque état $\varphi$ en fournit une, $A\mapsto\varphi(A)1$.) Ce sont deux
instances d'une même idée : projeter sur les opérateurs qui « ne voient pas » une certaine
structure.

> **Interprétation physique.** Dans l'algèbre de matrices complète $M_D(\mathbb C)$, $D=\sum_ck_c$,
> l'application $X\mapsto\sum_cP_cXP_c$ est l'espérance conditionnelle sur $\mathcal M$ : elle
> efface les cohérences entre secteurs de supersélection, ce qui se produit lorsque l'on mesure la
> charge (le secteur) sans enregistrer le résultat. L'espérance conditionnelle $E$ sur le centre va
> plus loin : elle ne garde que l'information classique, le poids $p_c$ de chaque secteur, et
> remplace l'état à l'intérieur de chaque secteur par l'état maximalement mélangé. Un état
> $\omega$ et l'état $\omega\circ E$ coïncident sur toute observable qui ne demande que « quel
> secteur ? ».

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
   simples. Un opérateur qui commute avec $S$ préserve ses sous-espaces propres ; dans la
   base propre, il est donc diagonal : $\{S\}'\cong\mathbb C^3$, ce sont les matrices
   circulantes. Alors $\{S\}''$ est lui aussi $\cong\mathbb C^3$, car dans la base propre de
   $S$ une matrice qui commute avec toutes les matrices diagonales est diagonale (comparer
   les coefficients de $XD$ et de $DX$ pour une matrice $D$ à coefficients diagonaux
   distincts) ; il est engendré par $1,S,S^2$ (Vandermonde : les valeurs propres sont
   distinctes).

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
