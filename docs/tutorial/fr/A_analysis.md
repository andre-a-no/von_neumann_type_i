[English](../A_analysis.md) · [Español](../es/A_analysis.md) · **Français** · [Deutsch](../de/A_analysis.md) · [中文](../zh/A_analysis.md) · [日本語](../ja/A_analysis.md) · [Русский](../ru/A_analysis.md)

# Annexe A. Résultats d'analyse et de probabilités

Les chapitres utilisent quelques résultats d'analyse et de probabilités qui dépassent les prérequis du cours sans
être des théorèmes à part entière : les ensembles de zéros des polynômes sont négligeables, les valeurs propres ordonnées
d'une matrice hermitienne en dépendent continûment, les images linéaires non dégénérées de vecteurs gaussiens ont une
densité, la direction d'un vecteur gaussien est uniforme sur la sphère et indépendante de sa longueur, les moments
gaussiens sont des sommes sur les appariements, et les espacements de points uniformes sont des variables exponentielles
normalisées. Cette annexe les démontre, afin que les chapitres puissent les citer par leur numéro.

**Prérequis utilisés.** Algèbre linéaire ; analyse réelle (limites, continuité, compacité et le théorème de
Bolzano–Weierstrass dans $\mathbb R^n$, séries entières d'une variable, intégrale de Riemann d'une fonction continue,
formule du changement de variables) ; mesure et intégration (mesure de Lebesgue, mesures produits et Fubini–Tonelli,
convergence monotone et dominée, mesures images, densités, indépendance) ; probabilités élémentaires (la loi $N(0,s^2)$,
de densité $(2\pi s^2)^{-1/2}e^{-y^2/2s^2}$ et de fonction génératrice des moments $\mathbb Ee^{aY}=e^{a^2s^2/2}$). Pour
une application linéaire inversible $L$ de $\mathbb R^m$ et un borélien $E$, $\lambda_m(LE)=|\det L|\,\lambda_m(E)$, où
$\lambda_m$ est la mesure de Lebesgue ; c'est le cas linéaire de la formule du changement de variables.

Dans toute la suite, « mesurable » signifie borélien, $|x|$ est la norme euclidienne de $x\in\mathbb R^n$ ou
$\mathbb C^n$, $\|A\|=\sup_{|x|=1}|Ax|$ la norme d'opérateur et $\|A\|_2=(\operatorname{Tr}A^*A)^{1/2}$ la norme de
Frobenius.

## A.1 Mesures et intégrales

**Lemme A.1 (Dynkin).** Soient $\mu,\nu$ des mesures finies sur un espace mesurable $(\Omega,\mathcal F)$ avec
$\mu(\Omega)=\nu(\Omega)$, et soit $\mathcal P\subseteq\mathcal F$ stable par intersections finies et engendrant
$\mathcal F$ comme $\sigma$-algèbre. Si $\mu(P)=\nu(P)$ pour tout $P\in\mathcal P$, alors $\mu=\nu$.

*Démonstration.* Appelons *système de Dynkin* une famille $\mathcal D\subseteq\mathcal F$ telle que $\Omega\in\mathcal D$,
$B\setminus A\in\mathcal D$ dès que $A\subseteq B$ sont dans $\mathcal D$, et $\bigcup_kA_k\in\mathcal D$ pour toute suite croissante
$A_k\in\mathcal D$. La famille $\{A\in\mathcal F:\mu(A)=\nu(A)\}$ est un système de Dynkin, puisque
$\mu(B\setminus A)=\mu(B)-\mu(A)$ pour des mesures finies et que les mesures sont continues le long des réunions
croissantes. Soit $d$ l'intersection de tous les systèmes de Dynkin contenant $\mathcal P$, qui est elle-même un système de Dynkin.
La famille $\mathcal D_1=\{A\in d:A\cap P\in d\text{ pour tout }P\in\mathcal P\}$ est un système de Dynkin (pour
$A\subseteq A'$, $(A'\setminus A)\cap P=(A'\cap P)\setminus(A\cap P)$) et contient $\mathcal P$, donc
$\mathcal D_1=d$. Alors $\mathcal D_2=\{A\in d:A\cap B\in d\text{ pour tout }B\in d\}$ est un système de Dynkin qui
contient $\mathcal P$ (car $\mathcal D_1=d$), donc $\mathcal D_2=d$ : $d$ est stable par intersections finies.
Par conséquent $d$ est une $\sigma$-algèbre : $\Omega\setminus A\in d$, $A\cup B=\Omega\setminus\big((\Omega\setminus A)\cap(\Omega\setminus B)\big)\in d$,
et une réunion dénombrable est la réunion croissante de réunions finies. Ainsi
$\mathcal F=\sigma(\mathcal P)\subseteq d\subseteq\{\mu=\nu\}$. $\square$

**Corollaire A.2.** Deux mesures boréliennes finies $\mu,\nu$ sur un espace métrique $M$ telles que $\int f\,d\mu=\int f\,d\nu$
pour toute fonction continue bornée $f:M\to\mathbb R$ sont égales.

*Démonstration.* Avec $f=1$, $\mu(M)=\nu(M)$. Pour un fermé non vide $K$ (pour $K=\emptyset$ il n'y a rien à démontrer), les fonctions continues
$f_k(x)=\max(0,1-k\operatorname{dist}(x,K))$ sont à valeurs dans $[0,1]$ et décroissent vers l'indicatrice de $K$
(pour $x\notin K$, $\operatorname{dist}(x,K)>0$). Par convergence dominée, $\mu(K)=\nu(K)$. Les fermés sont
stables par intersections finies et engendrent les boréliens ; on applique le lemme A.1. $\square$

Le lemme 5.1 utilise exactement cet argument.

**Lemme A.3 (mesure de Lebesgue sur un espace vectoriel).** Soit $W$ un espace vectoriel réel de dimension $m$. Pour une
bijection linéaire $\varphi:\mathbb R^m\to W$, on pose $\lambda_\varphi(E)=\lambda_m(\varphi^{-1}E)$ pour tout borélien $E\subseteq W$.

(a) Pour une autre bijection linéaire $\psi$, $\lambda_\psi=|\det(\psi^{-1}\varphi)|\,\lambda_\varphi$. Ainsi, les ensembles négligeables de
$W$, et le fait qu'un vecteur aléatoire de $W$ ait une densité, ne dépendent pas du choix de $\varphi$.

(b) Pour toute bijection linéaire $L:W\to W$, $\lambda_\varphi(LE)=|\det L|\,\lambda_\varphi(E)$. Si $L$ est orthogonale
pour un produit scalaire sur $W$, alors $|\det L|=1$ et $L$ préserve $\lambda_\varphi$.

(c) Si un vecteur aléatoire $X$ de $W$ a la densité $f$ par rapport à $\lambda_\varphi$, alors $LX$ a la densité
$y\mapsto f(L^{-1}y)/|\det L|$.

*Démonstration.* (a) $\lambda_\psi(E)=\lambda_m\big((\psi^{-1}\varphi)\varphi^{-1}E\big)=|\det(\psi^{-1}\varphi)|\,\lambda_m(\varphi^{-1}E)$.
(b) $\lambda_\varphi(LE)=\lambda_m\big((\varphi^{-1}L\varphi)\varphi^{-1}E\big)=|\det(\varphi^{-1}L\varphi)|\,\lambda_\varphi(E)$, et
$\det(\varphi^{-1}L\varphi)=\det L$ ne dépend pas de la base. Dans une base orthonormée, une application orthogonale $L$ a une
matrice $O$ avec $O^TO=1$, donc $(\det L)^2=1$. (c) D'après (b), la mesure image de $\lambda_\varphi$ par $L$ est
$|\det L|^{-1}\lambda_\varphi$, donc, pour tout borélien $E$,
$\mathbb P(LX\in E)=\int1_E(Lx)f(x)\,d\lambda_\varphi(x)=|\det L|^{-1}\int1_E(y)f(L^{-1}y)\,d\lambda_\varphi(y)$. $\square$

*Conventions.* « Mesure de Lebesgue » sur un espace vectoriel réel désigne n'importe quelle $\lambda_\varphi$ ; d'après (a), le choix
n'importe jamais pour les ensembles négligeables ni pour les densités. Sur $\mathbb C^n$ et sur $M_{a\times b}(\mathbb C)$, on prend comme
coordonnées les parties réelles et imaginaires des coefficients ; elles sont orthonormées pour le produit scalaire réel
$\operatorname{Re}\operatorname{Tr}(X^*Y)$, pour lequel toute application $X\mapsto VXW$ avec $V,W$ unitaires est orthogonale
(elle préserve $\operatorname{Tr}(X^*Y)$). Sur les matrices hermitiennes $\operatorname{Herm}_n$, on prend les coordonnées
$h_{ii}$ et $\operatorname{Re}h_{ij},\operatorname{Im}h_{ij}$ ($i<j$) ; elles ne sont pas orthonormées pour le produit scalaire réel
$\operatorname{Tr}(XY)$, mais d'après (b) les applications $H\mapsto UHU^*$ ($U$ unitaire), qui sont orthogonales pour celui-ci,
préservent encore la mesure de Lebesgue. De même, sur les matrices symétriques réelles $\operatorname{Sym}_n(\mathbb R)$, on prend les
coordonnées $h_{ii}$ et $h_{ij}$ ($i<j$), et les applications $H\mapsto OHO^T$ ($O\in O(n)$), orthogonales pour $\operatorname{Tr}(XY)$,
préservent la mesure de Lebesgue.

**Lemme A.4 (intégrales de fonctions à valeurs vectorielles).** Soient $V$ un espace vectoriel réel ou complexe de dimension
finie et $F:[a,b]\to V$ continue. On définit $\int_a^bF(t)\,dt$ coordonnée par coordonnée dans une base (parties réelle et
imaginaire séparément). Alors :

(a) pour des subdivisions $a=t_0<\dots<t_N=b$ pointées par $\tau_i\in[t_{i-1},t_i]$, les sommes de Riemann
$\sum_iF(\tau_i)(t_i-t_{i-1})$ convergent vers $\int_a^bF$ lorsque $\max_i(t_i-t_{i-1})\to0$ ; en particulier, l'intégrale
ne dépend pas de la base ;

(b) $L\int_a^bF=\int_a^bLF$ pour toute application ($\mathbb R$-)linéaire $L:V\to V'$ vers un espace $V'$ de dimension finie ;

(c) $\|\int_a^bF\|\le\int_a^b\|F(t)\|\,dt$ pour toute norme sur $V$ ;

(d) si tous les $F(t)$ appartiennent à un convexe fermé $K$, alors $\frac1{b-a}\int_a^bF\in K$ ; si de plus $K$ est un cône,
$\int_a^bF\in K$ ;

(e) $t\mapsto\int_a^tF$ est dérivable de dérivée $F(t)$ ; $\int_a^bF=\int_a^cF+\int_c^bF$, et
$\int_a^bF(s+c)\,ds=\int_{a+c}^{b+c}F(s)\,ds$.

*Démonstration.* Chaque coordonnée de $F$ est une fonction réelle continue, donc (a) et (e) valent coordonnée par coordonnée ;
(b) découle de (a), car $L$ envoie les sommes de Riemann sur des sommes de Riemann et est continue. Une norme est continue
($|\|x\|-\|y\||\le\|x-y\|\le\sum_k|x_k-y_k|\,\|e_k\|$), donc $t\mapsto\|F(t)\|$ est continue, et (c) découle de
(a) et de l'inégalité triangulaire pour les sommes de Riemann. Pour (d), on prend les subdivisions en $N$ parties égales : les
sommes de Riemann divisées par $b-a$ sont des combinaisons convexes $\frac1N\sum_iF(\tau_i)\in K$, et $K$ est fermé. $\square$

Ceci justifie le lemme 8.4(c) (sommes de Riemann d'applications CP) et les moyennes de Cesàro
$\sigma_\tau=\frac1\tau\int_0^\tau T_s\rho\,ds$ de la proposition 8.12 : $\sigma_\tau$ est un état d'après (d), et
$\|T_t\sigma_\tau-\sigma_\tau\|_1\le2t/\tau$ d'après (b), (c) et (e).

## A.2 Ensembles de zéros des fonctions analytiques

Pour un multi-indice $\alpha\in\mathbb N^n$, on note $|\alpha|=\alpha_1+\dots+\alpha_n$ et $y^\alpha=y_1^{\alpha_1}\cdots y_n^{\alpha_n}$.
Pour $a\in\mathbb R^n$ et $0<\rho\le\infty$, soit $Q_\rho(a)=\{x\in\mathbb R^n:\max_i|x_i-a_i|<\rho\}$ le cube
ouvert ($Q_\infty(a)=\mathbb R^n$).

**Définition.** Une fonction $f:Q_\rho(a)\to\mathbb R$ est une **série entière sur $Q_\rho(a)$** si
$f(x)=\sum_\alpha c_\alpha(x-a)^\alpha$ avec $\sum_\alpha|c_\alpha|r^{|\alpha|}<\infty$ pour tout $0<r<\rho$.

Les polynômes sont des séries entières sur $\mathbb R^n$ (un nombre fini de $c_\alpha\neq0$) ; en plusieurs variables, le cours n'applique
le théorème A.6 qu'aux polynômes. (Une fonction qui n'est que localement une série entière sur un ouvert connexe exigerait en
outre le principe du prolongement analytique en plusieurs variables, que le cours n'utilise pas.) Sur le cube fermé de rayon
$r<\rho$, les termes sont majorés par $|c_\alpha|r^{|\alpha|}$, donc la série y converge absolument et uniformément
(critère de Weierstrass) : sa somme ne dépend pas de l'ordre des termes, et $f$ est continue.

**Lemme A.5 (une variable).** Soit $g(t)=\sum_kc_k(t-a)^k$ une série entière sur $I=(a-\rho,a+\rho)$.

(a) Pour $b\in I$ et $|t-b|<\rho_b=\rho-|b-a|$, on a $g(t)=\sum_jd_j(t-b)^j$ avec
$d_j=\sum_{k\ge j}\binom kjc_k(b-a)^{k-j}$, une série entière sur $(b-\rho_b,b+\rho_b)$.

(b) Si l'un des $d_j$ est non nul, alors $g$ n'a pas de zéro dans un voisinage épointé de $b$.

(c) Si les zéros de $g$ s'accumulent en un point de $I$, alors $g=0$ sur $I$. Si $g\not\equiv0$ sur $I$, son ensemble de zéros est
dénombrable, donc Lebesgue-négligeable. Si $g\equiv0$ sur $I$, alors tous les $c_k$ sont nuls.

*Démonstration.* (a) Pour $0\le s<\rho_b$,
$\sum_k\sum_{j\le k}|c_k|\binom kjs^j|b-a|^{k-j}=\sum_k|c_k|(s+|b-a|)^k<\infty$, donc la série double obtenue en
développant $(t-a)^k=\big((t-b)+(b-a)\big)^k$ par la formule du binôme converge absolument pour $|t-b|\le s$ et peut être
sommée dans un ordre quelconque ; en sommant d'abord sur $k$, on obtient les $d_j$, avec $\sum_j|d_j|s^j<\infty$.
(b) Soit $j_0$ le plus petit $j$ tel que $d_j\ne0$. Alors $g(t)=(t-b)^{j_0}h(t)$ avec
$h(t)=\sum_{j\ge j_0}d_j(t-b)^{j-j_0}$, qui est continue près de $b$ (critère de Weierstrass) avec $h(b)=d_{j_0}\neq0$ ; donc
$h\neq0$ près de $b$. (c) Soit $A$ l'ensemble des $b\in I$ au voisinage desquels $g$ est identiquement nulle ; $A$ est ouvert. Si les
zéros s'accumulent en $b\in I$, alors, d'après (b), tous les $d_j$ sont nuls et $b\in A$ d'après (a). Si $b\in I$ est limite de
points de $A$, alors les zéros s'accumulent en $b$, donc $b\in A$ : $A$ est fermé dans $I$. Comme $I$ est connexe, $A$ est
vide ou égal à $I$. Si $g\not\equiv0$, alors $A=\emptyset$, donc aucun point de $I$ n'est point d'accumulation de zéros ; tout
intervalle compact $J\subset I$ ne contient alors qu'un nombre fini de zéros (une partie infinie de $J$ a un
point d'accumulation dans $J$ d'après Bolzano–Weierstrass), et $I$ est réunion dénombrable de tels $J$. Si $g\equiv0$,
(b) en $b=a$ (où $d_j=c_j$) montre que tous les $c_k$ sont nuls. $\square$

**Théorème A.6 (les ensembles de zéros sont négligeables).** Soit $f$ une série entière sur un cube $Q\subseteq\mathbb R^n$. Si $f$ s'annule
sur une partie de $Q$ de mesure de Lebesgue strictement positive, alors $f=0$ sur $Q$. De manière équivalente : si $f\not\equiv0$ sur $Q$,
alors $\{x\in Q:f(x)=0\}$ est Lebesgue-négligeable.

*Démonstration.* Récurrence sur $n$ ; le cas $n=1$ est le lemme A.5(c). Soit $n\ge2$ ; écrivons $x=(x',t)$ avec
$x'\in\mathbb R^{n-1}$, $Q=Q'\times I$ et $a=(a',a_n)$, et regroupons les termes selon la puissance de $t-a_n$ :
$$
f(x',t)=\sum_kh_k(x')(t-a_n)^k,\qquad h_k(x')=\sum_{\alpha'}c_{(\alpha',k)}(x'-a')^{\alpha'} .
$$
C'est légitime : pour $x'\in Q'$ et $\max_i|x'_i-a'_i|\le r<\rho$,
$\sum_k\sum_{\alpha'}|c_{(\alpha',k)}||x'-a'|^{\alpha'}r^k\le\sum_\alpha|c_\alpha|r^{|\alpha|}<\infty$ (avec
$|y|^{\alpha'}=\prod_i|y_i|^{\alpha'_i}$). La même majoration montre que chaque $h_k$ est une série entière sur $Q'$ et que, pour
$x'\in Q'$ fixé, $g_{x'}(t)=f(x',t)$ est une série entière sur $I$ de coefficients $h_k(x')$.

Soit $Z=\{x\in Q:f(x)=0\}$, fermé dans $Q$ puisque $f$ est continue, et supposons $\lambda_n(Z)>0$. D'après le théorème de
Tonelli, $\lambda_n(Z)=\int_{Q'}\lambda_1(Z_{x'})\,dx'$ avec $Z_{x'}=\{t\in I:f(x',t)=0\}$, donc l'ensemble
$F=\{x'\in Q':\lambda_1(Z_{x'})>0\}$ est de mesure $(n-1)$-dimensionnelle strictement positive. Pour $x'\in F$, le lemme A.5(c) donne
$g_{x'}\equiv0$ sur $I$, puis $h_k(x')=0$ pour tout $k$. Ainsi chaque $h_k$ s'annule sur $F$ et, par l'hypothèse de
récurrence, $h_k=0$ sur $Q'$. Par conséquent $f(x',t)=\sum_kh_k(x')(t-a_n)^k=0$ sur $Q$. $\square$

**Corollaire A.7.** (a) Un polynôme réel $p$ sur $\mathbb R^n$ qui n'est pas identiquement nul ne s'annule que sur un ensemble
Lebesgue-négligeable. Un polynôme ayant un coefficient non nul n'est pas identiquement nul.

(b) Soient $W$ un espace vectoriel réel de dimension finie et $p:W\to\mathbb R$ un polynôme en les coordonnées d'une
base (c'en est alors un dans toute base), non identiquement nul. Si un vecteur aléatoire $X$ de $W$ a une densité,
alors $\mathbb P(p(X)=0)=0$.

(c) Un sous-espace linéaire propre $U\subsetneq W$ est négligeable, et un vecteur aléatoire à densité est en dehors de
$U$ presque sûrement. En particulier, si $E\neq0$ est un sous-espace de $\mathbb C^n$ et si $v$ a une densité sur $\mathbb C^n$
(par exemple un vecteur initial gaussien), alors $v$ a une composante non nulle dans $E$ presque sûrement.

(d) Soient $X$ et $Y$ des vecteurs aléatoires indépendants dans des espaces réels de dimension finie $W_1$ et $W_2$, $Y$ à
densité, et $N\subseteq W_1\times W_2$ un borélien dont les sections $N_x=\{y:(x,y)\in N\}$ sont négligeables pour tout
$x$. Alors $\mathbb P((X,Y)\in N)=0$.

*Démonstration.* (a) La première assertion est le théorème A.6 avec $Q=\mathbb R^n$. Pour la seconde, récurrence sur $n$ : on écrit
$p(x',t)=\sum_kp_k(x')t^k$ ; l'un des $p_k$ a un coefficient non nul, donc $p_k(x'_0)\neq0$ pour un certain $x'_0$ par
récurrence (pour $n=1$, les $p_k$ sont des constantes et il n'y a pas de $x'$) ; alors $t\mapsto p(x'_0,t)$ est un polynôme non nul d'une variable, qui n'a qu'un nombre fini de racines.
(b) Un changement de base est linéaire, donc $p$ reste un polynôme, et $\{p=0\}$ est négligeable d'après (a) et le lemme A.3(a) ;
alors $\mathbb P(p(X)=0)=\int_{\{p=0\}}f\,d\lambda=0$. (c) On choisit une forme linéaire non nulle $\ell$ qui s'annule sur $U$
et on applique (b) à $p=\ell$. Pour l'exemple, $v\mapsto P_Ev$ ($P_E$ la projection orthogonale) est une application
$\mathbb R$-linéaire non nulle, et $\{P_Ev=0\}$ est son noyau, un sous-espace propre. (d) La loi de $(X,Y)$ est le produit des lois, donc, d'après
le théorème de Tonelli, $\mathbb P((X,Y)\in N)=\int\mathbb P(Y\in N_x)\,d\mathbb P_X(x)=0$. $\square$

Le cours applique le corollaire A.7 à $\det Z$ pour une matrice de Ginibre $Z$ ($|\det Z|^2$ est un polynôme en les parties réelles
et imaginaires des coefficients, égal à $1$ en $Z=1$ ; théorème 5.3, section 1.6), au déterminant de Gram de la
proposition 4.7(a), au rang des matrices gaussiennes (chapitre 12, exercice 4) et aux vecteurs initiaux de la
méthode de Lanczos (chapitre 9). La partie (d) est l'étape « comme $Y$ est indépendant de $K$, le théorème de Tonelli donne
la probabilité $0$ » de la proposition 4.7(b).

```python
import math
import numpy as np

rng = np.random.default_rng(0)


def mc(samples, exact, label, kind="exact"):
    """Monte Carlo mean with standard error; the reference value must lie within 5 standard errors."""
    samples = np.asarray(samples, dtype=float)
    mean, err = samples.mean(), samples.std() / math.sqrt(samples.size)
    print(f"{label}: {mean:.4f} +- {err:.4f}   ({kind} {exact:.4f})")
    assert abs(mean - exact) <= 5 * err + 1e-12, label


# Corollary A.7: det Z = 0 has probability 0 for a real Gaussian 3 x 3 matrix, so by continuity of the measure
# P(|det Z| < eps) -> 0 as eps -> 0. Numerical observation (not a theorem of this appendix): the fractions
# below shrink by about a factor 10 per factor 10 in eps.
d = np.abs(np.linalg.det(rng.standard_normal((200_000, 3, 3))))
fractions = [(d < eps).mean() for eps in (1e-1, 1e-2, 1e-3)]
print("P(|det Z| < eps) for eps = 0.1, 0.01, 0.001:", fractions)
assert fractions[1] < 0.2 * fractions[0] and fractions[2] < 0.2 * fractions[1]
```

## A.3 Valeurs propres des matrices hermitiennes

Pour $A\in\operatorname{Herm}_n$, soient $\lambda_1(A)\le\dots\le\lambda_n(A)$ ses valeurs propres comptées avec multiplicité, et
$v_1,\dots,v_n$ une base orthonormée de vecteurs propres (théorème 1.3, corollaire 1.4). Pour un vecteur unitaire $x$,
$\langle x,Ax\rangle=\sum_i\lambda_i|\langle v_i,x\rangle|^2$ avec $\sum_i|\langle v_i,x\rangle|^2=1$, donc
$\lambda_1(A)\le\langle x,Ax\rangle\le\lambda_n(A)$ ; et $\|A\|=\max_i|\lambda_i(A)|$ (corollaire 1.4).

**Théorème A.8 (Courant–Fischer).** Pour $A\in\operatorname{Herm}_n$ et $1\le j\le n$,
$$
\lambda_j(A)=\min_{\dim U=j}\ \max_{x\in U,\,|x|=1}\langle x,Ax\rangle=\max_{\dim U=n-j+1}\ \min_{x\in U,\,|x|=1}\langle x,Ax\rangle ,
$$
sur les sous-espaces $U\subseteq\mathbb C^n$ de la dimension indiquée. (Les maxima et minima intérieurs existent car la
sphère unité de $U$ est compacte ; les extérieurs sont atteints.)

*Démonstration.* Soient $\dim U=j$ et $W=\operatorname{span}(v_j,\dots,v_n)$. Comme $\dim U+\dim W=n+1>n$, il existe un vecteur unitaire
$x\in U\cap W$, et $\langle x,Ax\rangle=\sum_{i\ge j}\lambda_i|\langle v_i,x\rangle|^2\ge\lambda_j$. Donc le maximum sur $U$
est $\ge\lambda_j$. Pour $U=\operatorname{span}(v_1,\dots,v_j)$, $\langle x,Ax\rangle=\sum_{i\le j}\lambda_i|\langle v_i,x\rangle|^2\le\lambda_j$,
avec égalité en $x=v_j$ ; donc le minimum sur $U$ est atteint et vaut $\lambda_j$. Pour la seconde formule,
on applique la première à $-A$, dont les valeurs propres sont $\lambda_j(-A)=-\lambda_{n+1-j}(A)$ :
$-\lambda_j(A)=\lambda_{n+1-j}(-A)=\min_{\dim U=n+1-j}\max_x(-\langle x,Ax\rangle)=-\max_{\dim U=n+1-j}\min_x\langle x,Ax\rangle$. $\square$

**Corollaire A.9 (compressions).** Soient $Q\in M_{n\times m}(\mathbb C)$ avec $Q^*Q=1_m$ et $A\in\operatorname{Herm}_n$. Alors
$\lambda_j(Q^*AQ)=\min\max\langle x,Ax\rangle$ sur les sous-espaces $U\subseteq\operatorname{ran}Q$ de dimension $j$ et les $x\in U$ unitaires, et
$$
\lambda_j(A)\le\lambda_j(Q^*AQ)\le\lambda_{j+n-m}(A)\qquad(1\le j\le m).
$$

*Démonstration.* $U'\mapsto QU'$ est une bijection des sous-espaces de dimension $j$ de $\mathbb C^m$ sur ceux de
$\operatorname{ran}Q$, $|Qy|=|y|$ et $\langle y,Q^*AQy\rangle=\langle Qy,AQy\rangle$ ; le théorème A.8 pour $Q^*AQ$ donne donc la
formule. En minimisant sur la famille plus grande de tous les $U\subseteq\mathbb C^n$ de dimension $j$, on obtient $\lambda_j(A)$, qui est
donc plus petit. De même, la formule max–min de $\lambda_j(Q^*AQ)$ porte sur les $U\subseteq\operatorname{ran}Q$ de dimension
$m-j+1$, et sur tous les $U\subseteq\mathbb C^n$ de cette dimension elle donne le nombre plus grand
$\lambda_{n-(m-j+1)+1}(A)=\lambda_{j+n-m}(A)$. $\square$

**Théorème A.10 (Weyl).** Pour $A,E\in\operatorname{Herm}_n$ et tout $j$,
$$
\lambda_j(A)+\lambda_1(E)\le\lambda_j(A+E)\le\lambda_j(A)+\lambda_n(E),
\qquad\text{hence}\qquad
|\lambda_j(A+E)-\lambda_j(A)|\le\|E\|\le\|E\|_2 .
$$

*Démonstration.* Pour $x$ unitaire, $\langle x,(A+E)x\rangle\le\langle x,Ax\rangle+\lambda_n(E)$. On prend le maximum sur les vecteurs unitaires d'un
sous-espace $U$ de dimension $j$, puis le minimum sur $U$ (théorème A.8) : $\lambda_j(A+E)\le\lambda_j(A)+\lambda_n(E)$. Appliqué à
$A+E$ et $-E$, cela donne $\lambda_j(A)\le\lambda_j(A+E)+\lambda_n(-E)=\lambda_j(A+E)-\lambda_1(E)$. Donc
$|\lambda_j(A+E)-\lambda_j(A)|\le\max(|\lambda_1(E)|,|\lambda_n(E)|)=\|E\|$, et $\|E\|^2=\max_i\lambda_i(E)^2\le\sum_i\lambda_i(E)^2=\|E\|_2^2$. $\square$

**Corollaire A.11.** Chaque $\lambda_j:\operatorname{Herm}_n\to\mathbb R$ est $1$-lipschitzienne pour $\|\cdot\|$, donc continue et
borélienne. Pour une matrice hermitienne aléatoire $H$, $(\lambda_1(H),\dots,\lambda_n(H))$ est un vecteur aléatoire, de même que ses
espacements et ses rapports d'espacements (sur le borélien où leurs dénominateurs sont non nuls) ; les vecteurs de valeurs propres de matrices hermitiennes aléatoires indépendantes sont indépendants. Les matrices ayant
$n$ valeurs propres distinctes forment un ouvert.

*Démonstration.* Théorème A.10 avec $E=B-A$. Des fonctions de variables aléatoires indépendantes sont indépendantes. Si
$\lambda_j(A)<\lambda_{j+1}(A)$ pour tout $j$, les inégalités strictes persistent pour $\|B-A\|<\frac12\min_j(\lambda_{j+1}(A)-\lambda_j(A))$. $\square$

**Lemme A.12 (formule de Lagrange).** Soit $A\in\operatorname{Herm}_n$ de valeurs propres distinctes $\mu_1<\dots<\mu_r$ et de
projections spectrales $P_1,\dots,P_r$ (projections orthogonales sur les sous-espaces propres, de sorte que $A=\sum_i\mu_iP_i$,
$P_iP_l=\delta_{il}P_i$, $\sum_iP_i=1$). Alors $p(A)=\sum_ip(\mu_i)P_i$ pour tout polynôme $p$, et
$$
P_i=\prod_{l\neq i}\frac{A-\mu_l}{\mu_i-\mu_l} .
$$

*Démonstration.* $A^k=\sum_i\mu_i^kP_i$ par récurrence sur $k$ (en utilisant $P_iP_l=\delta_{il}P_i$), et $p(A)=\sum_ip(\mu_i)P_i$ par linéarité.
Le polynôme $\ell_i(t)=\prod_{l\neq i}(t-\mu_l)/(\mu_i-\mu_l)$ vérifie $\ell_i(\mu_l)=\delta_{il}$, donc $\ell_i(A)=P_i$. $\square$

**Proposition A.13 (continuité des valeurs propres et des projections spectrales).** Soient $S$ une partie d'un espace vectoriel réel
de dimension finie et $A:S\to\operatorname{Herm}_n$ continue, et supposons que $A(x)$ ait exactement $r$ valeurs propres distinctes
$\mu_1(x)<\dots<\mu_r(x)$ pour tout $x\in S$. Alors leurs multiplicités sont localement constantes sur $S$, et les $\mu_i$ ainsi que
les projections spectrales $P_i(x)$ de $A(x)$ sont continues sur $S$. En particulier, sur l'ouvert des matrices hermitiennes
à spectre simple, les valeurs propres et les projections spectrales de rang un sont continues.

*Démonstration.* Si $r=1$, $A(x)=\mu_1(x)1$ avec $\mu_1=\operatorname{Tr}A/n$ et $P_1=1$. Soit $r\ge2$ ; fixons $x_0\in S$ et posons
$d=\min_i(\mu_{i+1}(x_0)-\mu_i(x_0))>0$. Soit $J_i$ l'ensemble des indices $j$ tels que $\lambda_j(A(x_0))=\mu_i(x_0)$, un bloc
d'entiers consécutifs. Par continuité de $A$, il existe un voisinage $N$ de $x_0$ dans $S$ tel que $\|A(x)-A(x_0)\|<d/2$ pour
$x\in N$ ; d'après le théorème A.10, $\lambda_j(A(x))$ appartient à l'intervalle $I_i=(\mu_i(x_0)-d/2,\mu_i(x_0)+d/2)$ pour $j\in J_i$, et
ces $r$ intervalles sont disjoints. Donc $A(x)$ a au moins une valeur propre distincte dans chaque $I_i$ et, en ayant exactement $r$,
exactement une : $\lambda_j(A(x))$ est la même pour tous les $j\in J_i$. Ainsi, sur $N$, les multiplicités sont les constantes $|J_i|$ et
$\mu_i(x)=\lambda_{\max J_i}(A(x))$ est continue (corollaire A.11). D'après le lemme A.12,
$P_i(x)=\prod_{l\neq i}(A(x)-\mu_l(x))/(\mu_i(x)-\mu_l(x))$, continue sur $N$ car les dénominateurs restent éloignés de $0$.
La dernière assertion est le cas $r=n$, $S$ étant l'ouvert du corollaire A.11. $\square$

Dans la proposition 4.7, $S$ est l'ensemble où le déterminant de Gram $\delta$ est non nul, sur lequel $K$ a exactement $N$ valeurs
propres distinctes ; ainsi $K\mapsto\Pi_i(K)$ y est continue et l'événement de la proposition 4.7(b) est borélien. Le théorème 5.15
utilise le corollaire A.11, et le théorème 9.14 utilise le corollaire A.9.

```python
def herm(n, batch):
    G = rng.standard_normal((batch, n, n)) + 1j * rng.standard_normal((batch, n, n))
    return (G + G.conj().transpose(0, 2, 1)) / 2


n = 6
A, E = herm(n, 2000), 0.1 * herm(n, 2000)
lam_A, lam_AE, lam_E = np.linalg.eigvalsh(A), np.linalg.eigvalsh(A + E), np.linalg.eigvalsh(E)   # ascending
normE = np.linalg.norm(E, ord=2, axis=(1, 2))                                    # operator norm
assert (np.abs(lam_AE - lam_A).max(axis=1) <= normE + 1e-12).all()               # Theorem A.10
assert (lam_AE >= lam_A + lam_E[:, :1] - 1e-12).all() and (lam_AE <= lam_A + lam_E[:, -1:] + 1e-12).all()

a0, j = A[0], 3                          # Theorem A.8: every 3-dim subspace has max <x, Ax> >= lambda_3
for _ in range(200):
    Q, _ = np.linalg.qr(rng.standard_normal((n, j)) + 1j * rng.standard_normal((n, j)))
    assert np.linalg.eigvalsh(Q.conj().T @ a0 @ Q)[-1] >= lam_A[0, j - 1] - 1e-12
m = 4                                    # Corollary A.9: interlacing for a compression to 4 dimensions
Q, _ = np.linalg.qr(rng.standard_normal((n, m)) + 1j * rng.standard_normal((n, m)))
theta = np.linalg.eigvalsh(Q.conj().T @ a0 @ Q)
assert (theta >= lam_A[0, :m] - 1e-12).all() and (theta <= lam_A[0, n - m:] + 1e-12).all()

U, _ = np.linalg.qr(rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)))
mu, mult = np.array([-1.0, 0.5, 2.0]), [2, 3, 1]          # Lemma A.12 with multiplicities 2, 3, 1
spec = np.repeat(mu, mult)
A1 = U @ np.diag(spec) @ U.conj().T
for i, m_i in enumerate(mult):
    P = np.eye(n, dtype=complex)
    for l in range(len(mu)):
        if l != i:
            P = P @ (A1 - mu[l] * np.eye(n)) / (mu[i] - mu[l])
    Vi = U[:, spec == mu[i]]
    assert np.abs(P - Vi @ Vi.conj().T).max() < 1e-10 and abs(np.trace(P).real - m_i) < 1e-10
```

## A.4 Vecteurs gaussiens

Un **vecteur gaussien standard** de $\mathbb R^m$ est $g=(g_1,\dots,g_m)$ à coordonnées indépendantes de loi $N(0,1)$. Sa loi
est le produit des lois des coordonnées, donc, d'après le théorème de Tonelli, il a la densité
$\varphi_m(x)=(2\pi)^{-m/2}e^{-|x|^2/2}$.

**Proposition A.14 (invariance orthogonale).** (a) $Og$ a la loi de $g$ pour toute matrice orthogonale $O\in O(m)$.

(b) Soit $W$ un espace euclidien réel de dimension $m$. Pour une base orthonormée $(w_i)$, la loi de
$\sum_ig_iw_i$ ne dépend pas de la base ; un tel vecteur est appelé **vecteur gaussien standard de $W$**. Ses coordonnées
dans toute base orthonormée sont indépendantes de loi $N(0,1)$, sa loi est invariante par les applications orthogonales de $W$, et si
$W=W_1\oplus W_2$ orthogonalement, ses composantes sont des vecteurs gaussiens standard indépendants de $W_1$ et de $W_2$.

(c) Pour $v\in\mathbb R^m$, $\langle v,g\rangle$ suit la loi $N(0,|v|^2)$ ; en particulier $\sum_ia_ig_i\sim N(0,\sum_ia_i^2)$.

*Démonstration.* (a) D'après le lemme A.3(c), $Og$ a la densité $\varphi_m(O^{-1}y)/|\det O|=\varphi_m(y)$. (b) Si $w'_i=\sum_jO_{ji}w_j$
est une autre base orthonormée, $O$ est orthogonale et $\sum_ig_iw'_i=\sum_j(Og)_jw_j$, qui a la loi de
$\sum_jg_jw_j$ d'après (a). Les coordonnées de $X=\sum_ig_iw_i$ dans une autre base orthonormée $(w'_i)$ sont $O^Tg$, de nouveau
standard d'après (a). Une application orthogonale $L$ envoie $(w_i)$ sur une base orthonormée, et $LX=\sum_ig_iLw_i$. Pour la dernière
assertion, on prend $(w_i)$ adaptée à $W_1\oplus W_2$. (c) Pour $v=0$, les deux membres sont la masse de Dirac en $0$. Pour $v\neq0$,
on complète $v/|v|$ en une base orthonormée ; d'après (b), la coordonnée $\langle v/|v|,g\rangle$ suit la loi $N(0,1)$, donc
$\langle v,g\rangle=|v|\,\langle v/|v|,g\rangle$ a la densité $(2\pi|v|^2)^{-1/2}e^{-y^2/2|v|^2}$ (lemme A.3(c) avec $m=1$). $\square$

**Définition.** Un **vecteur gaussien centré** d'un espace vectoriel réel de dimension finie $W$ est $X=Ag$ pour une application
linéaire $A:\mathbb R^M\to W$ et un vecteur gaussien standard $g$ de $\mathbb R^M$. En coordonnées ($W=\mathbb R^m$), sa
**covariance** est $C=\mathbb E[XX^T]=AA^T$. Il est **non dégénéré** si $A$ est surjective, ce qui équivaut à $C$ définie positive
($\langle x,Cx\rangle=|A^Tx|^2$).

Tout vecteur aléatoire du cours construit linéairement à partir de coefficients gaussiens indépendants — vecteurs gaussiens
complexes, matrices de Ginibre, matrices du GOE et du GUE, leurs images par des twirls — est un vecteur gaussien centré.

**Théorème A.15 (images linéaires et densités).** Soit $X=Ag$ un vecteur gaussien centré de $W$.

(a) Pour toute application linéaire $L:W\to W'$, $LX$ est un vecteur gaussien centré ; il est non dégénéré si $X$ l'est et si $L$ est surjective.

(b) Si $X$ est non dégénéré, il a une densité par rapport à la mesure de Lebesgue sur $W$ ; en coordonnées
$W=\mathbb R^m$,
$$
f_X(x)=(2\pi)^{-m/2}(\det C)^{-1/2}\,e^{-\langle x,C^{-1}x\rangle/2}.
$$
Si $X$ est dégénéré, il appartient au sous-espace propre $\operatorname{ran}A$, un ensemble négligeable (corollaire A.7(c)), et n'a pas de densité.

(c) La loi de $X$ ne dépend que de sa covariance.

*Démonstration.* (a) $LX=(LA)g$, et $LA$ est surjective si $A$ et $L$ le sont. (b) Soit $K=\ker A$, de sorte que
$\mathbb R^M=K^\perp\oplus K$. D'après la proposition A.14(b), $g=g_1+g_2$ avec $g_1$ un vecteur gaussien standard de
$K^\perp$, et $X=Ag_1=Bg_1$ avec $B=A|_{K^\perp}$, qui est injective et, $A$ étant surjective, bijective
$K^\perp\to W$ ; donc $\dim K^\perp=m$. Dans une base orthonormée de $K^\perp$ complétée par une base orthonormée de $K$, la matrice de $A$ est
$(B\ 0)$, donc $C=AA^T=BB^T$. Les coordonnées de $g_1$ ont la densité $\varphi_m$ et, d'après le lemme A.3(c), $X$ a
la densité $\varphi_m(B^{-1}x)/|\det B|$, avec $|B^{-1}x|^2=\langle x,(BB^T)^{-1}x\rangle$ et $|\det B|=(\det C)^{1/2}$. Si
$X$ était dégénéré et avait une densité $f$, on aurait $1=\mathbb P(X\in\operatorname{ran}A)=\int_{\operatorname{ran}A}f\,d\lambda=0$.
(c) $\ker C=\ker A^T$ (si $AA^Tx=0$, alors $|A^Tx|^2=0$), donc $\operatorname{ran}C=(\ker C)^\perp=(\ker A^T)^\perp=\operatorname{ran}A$
contient $X$. Soit $d=\operatorname{rank}C$ (si $d=0$, alors $X=0$) ; on choisit une base orthonormée $p_1,\dots,p_d$ de
$\operatorname{ran}C$ (qui ne dépend que de $C$) et on note $P$ la matrice $d\times m$ de lignes $p_k^T$, de sorte que $PP^T=1_d$ et que
$P^TP$ est la projection orthogonale sur $\operatorname{ran}C$. Alors $PX=(PA)g$ avec $PA$ surjective sur $\mathbb R^d$ : c'est un vecteur
gaussien centré non dégénéré de covariance $PCP^T$ ; d'après (b), sa loi ne dépend que de $C$, et il en va de même de la loi de
$X=P^T(PX)$. $\square$

Par exemple, au chapitre 4, le twirl $P$ est une projection linéaire de $M_n(\mathbb C)$ sur $\mathcal A'$, donc
$Z\mapsto P(Z)$ et $Z\mapsto P(Z+Z^*)$ sont des surjections $\mathbb R$-linéaires sur $\mathcal A'$ et sur $\mathcal A'_h$ (pour
$H\in\mathcal A'_h$, $P(\frac H2+\frac H2)=H$) ; pour une matrice de Ginibre $Z$ (non dégénérée d'après la proposition A.16), $Y=P(Z)$
et $K=P(Z+Z^*)$ ont des densités sur ces espaces d'après le théorème A.15. Une combinaison gaussienne $\sum_k\xi_kB_k$ d'une base
$(B_k)$ de $\mathcal A'$ est l'image du vecteur gaussien standard $\xi$ par une bijection linéaire.

**Proposition A.16 (vecteurs gaussiens complexes, matrices de Ginibre).** (a) Soit $z\in\mathbb C^n$ à coordonnées indépendantes
dont les parties réelle et imaginaire sont indépendantes de loi $N(0,\frac12)$ (un **vecteur gaussien complexe standard**).
Alors $\sqrt2z$ est le vecteur gaussien standard de $\mathbb C^n\cong\mathbb R^{2n}$ muni du produit scalaire
$\operatorname{Re}\langle u,w\rangle$ ; $z$ a la densité $\pi^{-n}e^{-|z|^2}$ ; $Uz$ a la loi de $z$ pour tout
$U\in U(n)$, donc les coordonnées $\langle w_i,z\rangle$ dans n'importe quelle base orthonormée sont de nouveau des gaussiennes
complexes standard indépendantes ; et $\mathbb E[z_i\bar z_j]=\delta_{ij}$, $\mathbb E[z_iz_j]=0$.

(b) Une matrice de Ginibre $G$ de taille $a\times b$ (coefficients gaussiens complexes standard indépendants) a la densité
$\pi^{-ab}e^{-\|G\|_2^2}$ sur $M_{a\times b}(\mathbb C)$, et $VGW$ a la loi de $G$ pour tous $V\in U(a)$, $W\in U(b)$.

(c) Une matrice réelle $G$ de taille $a\times b$ à coefficients indépendants de loi $N(0,1)$ a la densité $(2\pi)^{-ab/2}e^{-\|G\|_2^2/2}$, et
$OGO'$ a la loi de $G$ pour tous $O\in O(a)$, $O'\in O(b)$.

*Démonstration.* (a) Les $2n$ nombres $\sqrt2\operatorname{Re}z_i,\sqrt2\operatorname{Im}z_i$ sont indépendants de loi $N(0,1)$, et ce sont les
coordonnées de $\sqrt2z$ dans la base orthonormée $(e_i,ie_i)$. D'après le lemme A.3(c) avec $L=2^{-1/2}$, $z$ a la densité
$(2\pi)^{-n}e^{-|\sqrt2z|^2/2}2^{n}=\pi^{-n}e^{-|z|^2}$. Une matrice unitaire $U$ est $\mathbb R$-linéaire et préserve
$\operatorname{Re}\langle u,w\rangle$, donc elle est orthogonale ; on applique la proposition A.14(b). Les coordonnées dans la base $(w_i)$ sont
$W^*z$, où $W$ est la matrice unitaire de colonnes $w_i$. Enfin, avec $z_j=x_j+iy_j$,
$\mathbb E[z_i\bar z_j]=\mathbb E[x_ix_j+y_iy_j]+i\,\mathbb E[y_ix_j-x_iy_j]=\delta_{ij}$ et $\mathbb E[z_iz_j]=\mathbb E[x_ix_j-y_iy_j]+i\,\mathbb E[x_iy_j+y_ix_j]=0$.
(b) $M_{a\times b}(\mathbb C)\cong\mathbb C^{ab}$ avec $|G|=\|G\|_2$, et $X\mapsto VXW$ est unitaire pour
$\operatorname{Tr}(X^*Y)$ (conventions après le lemme A.3) ; on applique (a). (c) De même avec $\mathbb R^{ab}$ et la proposition A.14. $\square$

**Proposition A.17 (GOE et GUE).** Soit $H=(G+G^T)/2$ avec $G$ une matrice réelle $n\times n$ à coefficients indépendants de loi $N(0,1)$
(GOE), ou $H=(G+G^*)/\sqrt2$ avec $G$ une matrice de Ginibre (GUE). Alors les coefficients $h_{ij}$, $i\le j$, sont indépendants,
$h_{ii}\sim N(0,1)$ et, pour $i<j$, $h_{ij}\sim N(0,\frac12)$ (GOE), respectivement $\operatorname{Re}h_{ij},\operatorname{Im}h_{ij}$ sont
indépendants de loi $N(0,\frac12)$ (GUE). La densité de $H$ dans les coordonnées fixées après le lemme A.3 est
$$
(2\pi)^{-n/2}\pi^{-n(n-1)/4}e^{-\operatorname{Tr}H^2/2}\ \ \text{(GOE)},\qquad
(2\pi)^{-n/2}\pi^{-n(n-1)/2}e^{-\operatorname{Tr}H^2/2}\ \ \text{(GUE)},
$$
et $OHO^T$ ($O\in O(n)$), respectivement $UHU^*$ ($U\in U(n)$), a la loi de $H$.

*Démonstration.* GOE : $G$ est le vecteur gaussien standard de $M_n(\mathbb R)$ muni de $\operatorname{Tr}(X^TY)$, et les matrices
$e_{ii}$, $(e_{ij}\pm e_{ji})/\sqrt2$ ($i<j$) forment une base orthonormée. D'après la proposition A.14(b), les coordonnées
$G_{ii}$, $(G_{ij}+G_{ji})/\sqrt2$, $(G_{ij}-G_{ji})/\sqrt2$ sont indépendantes de loi $N(0,1)$. Or $h_{ii}=G_{ii}$ et
$h_{ij}=2^{-1/2}\cdot(G_{ij}+G_{ji})/\sqrt2$. GUE : $\sqrt2G$ est le vecteur gaussien standard de $M_n(\mathbb C)$ muni de
$\operatorname{Re}\operatorname{Tr}(X^*Y)$ (proposition A.16), avec la base orthonormée $e_{ii}$, $ie_{ii}$,
$(e_{ij}\pm e_{ji})/\sqrt2$, $i(e_{ij}\pm e_{ji})/\sqrt2$ ; les coordonnées de $\sqrt2G$ sont $\sqrt2\operatorname{Re}G_{ii}$,
$\sqrt2\operatorname{Im}G_{ii}$, $\operatorname{Re}(G_{ij}\pm G_{ji})$, $\operatorname{Im}(G_{ij}\pm G_{ji})$, indépendantes de loi $N(0,1)$. Or
$h_{ii}=\sqrt2\operatorname{Re}G_{ii}$, $\operatorname{Re}h_{ij}=2^{-1/2}\operatorname{Re}(G_{ij}+G_{ji})$ et
$\operatorname{Im}h_{ij}=2^{-1/2}\operatorname{Im}(G_{ij}-G_{ji})$. Dans les deux cas, la densité est le produit des
densités unidimensionnelles, $(2\pi)^{-1/2}e^{-h_{ii}^2/2}$ et $\pi^{-1/2}e^{-t^2}$ pour chaque coordonnée réelle $t$ de variance
$\frac12$, et $\operatorname{Tr}H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2$. L'application $H\mapsto UHU^*$ est linéaire, inversible et
orthogonale pour $\operatorname{Tr}(XY)$, donc, d'après le lemme A.3(b),(c), $UHU^*$ a la densité $f(U^*HU)=f(H)$, puisque $f$ ne dépend que de
$\operatorname{Tr}H^2$ ; de même pour $O$. $\square$

```python
Amat = rng.standard_normal((2, 4))                       # Theorem A.15: X = A g in R^2, A of rank 2
C = Amat @ Amat.T
X = rng.standard_normal((400_000, 4)) @ Amat.T
h = 1 / 400                                              # midpoint rule for the density on [0, 1]^2
grid = (np.arange(400) + 0.5) * h
xx = np.stack(np.meshgrid(grid, grid, indexing='ij'), -1).reshape(-1, 2)
dens = np.exp(-0.5 * np.einsum('ni,ij,nj->n', xx, np.linalg.inv(C), xx)) / (2 * math.pi * math.sqrt(np.linalg.det(C)))
mc(((X >= 0) & (X <= 1)).all(axis=1), dens.sum() * h * h, "P(X in [0,1]^2) against the integral of f_X",
   kind="midpoint rule")

nn, B = 4, 100_000                                       # Proposition A.17
G = rng.standard_normal((B, nn, nn))
goe = (G + G.transpose(0, 2, 1)) / 2
Gc = (rng.standard_normal((B, nn, nn)) + 1j * rng.standard_normal((B, nn, nn))) / math.sqrt(2)
gue = (Gc + Gc.conj().transpose(0, 2, 1)) / math.sqrt(2)
mc(np.einsum('bij,bji->b', goe, goe), nn * (nn + 1) / 2, "E Tr H^2, GOE")      # n + 2 * n(n-1)/2 * 1/2
mc(np.einsum('bij,bji->b', gue, gue).real, nn ** 2, "E Tr H^2, GUE")           # n + 2 * n(n-1)/2 * 1
O, _ = np.linalg.qr(rng.standard_normal((nn, nn)))
rot = O @ goe @ O.T                                                             # has the law of the GOE
mc(rot[:, 0, 0] ** 4, 3.0, "E (O H O^T)_11^4")
mc(rot[:, 0, 1] ** 2, 0.5, "E (O H O^T)_12^2")
```

## A.5 La direction d'un vecteur gaussien

Soient $S^{n-1}=\{x\in\mathbb R^n:|x|=1\}$, $V_n=\lambda_n(\{|x|\le1\})$ et, pour un borélien $B\subseteq S^{n-1}$,
$C(B)=\{x\in\mathbb R^n:0<|x|\le1,\ x/|x|\in B\}$, qui est borélien ($x\mapsto(|x|,x/|x|)$ est continue sur $\mathbb R^n\setminus0$).

**Théorème A.18 (coordonnées polaires).** Pour des boréliens $A\subseteq(0,\infty)$ et $B\subseteq S^{n-1}$,
$$
\lambda_n\big(\{x\neq0:|x|\in A,\ x/|x|\in B\}\big)=n\,\lambda_n(C(B))\int_Ar^{n-1}\,dr .
$$
Par conséquent, pour des fonctions boréliennes $h:(0,\infty)\to[0,\infty]$ et $\phi:S^{n-1}\to[0,\infty]$, avec $\mu(B)=\lambda_n(C(B))$,
$$
\int_{\mathbb R^n}h(|x|)\,\phi(x/|x|)\,dx=n\int_0^\infty h(r)r^{n-1}\,dr\int_{S^{n-1}}\phi\,d\mu .
$$

*Démonstration.* Fixons $B$ et notons $\nu_1(A)$ et $\nu_2(A)$ les deux membres, qui sont des mesures boréliennes sur $(0,\infty)$. L'homothétie
$x\mapsto rx$ a pour déterminant $r^n$ et envoie $C(B)$ sur $\{0<|x|\le r,\ x/|x|\in B\}$, puisque $rx/|rx|=x/|x|$ ; donc
$\nu_1((0,r])=r^n\lambda_n(C(B))=\nu_2((0,r])$. Pour $R>0$, les restrictions de $\nu_1,\nu_2$ à $(0,R]$ sont finies
(au plus $\lambda_n(\{|x|\le R\})$) et coïncident sur les intervalles $(0,r]$, $r\le R$, qui sont stables par intersection,
contiennent $(0,R]$ et engendrent les boréliens de $(0,R]$ (ils donnent $(s,r]=(0,r]\setminus(0,s]$, et toute partie relativement ouverte
de $(0,R]$ est réunion dénombrable de tels intervalles). D'après le lemme A.1, $\nu_1=\nu_2$ sur $(0,R]$ ; on fait tendre $R\to\infty$. La
formule intégrale vaut pour $h=1_A$, $\phi=1_B$ (l'identité qui vient d'être démontrée ; $\{0\}$ est négligeable), donc pour $h,\phi$ étagées
positives par linéarité, et pour $h,\phi$ quelconques par convergence monotone le long de fonctions étagées $h_k\uparrow h$, $\phi_k\uparrow\phi$. $\square$

**Corollaire A.19 (intégrales radiales).** Pour $h\ge0$ borélienne sur $(0,\infty)$,
$\int_{\mathbb R^n}h(|x|)\,dx=s_n\int_0^\infty h(r)r^{n-1}\,dr$ avec $s_n=nV_n=(2\pi)^{n/2}/I_n$,
$I_n=\int_0^\infty r^{n-1}e^{-r^2/2}dr$. Ici $I_1=\sqrt{\pi/2}$, $I_2=1$, $I_{n+2}=nI_n$, donc $s_1=2$, $s_2=2\pi$, $s_3=4\pi$,
$s_4=2\pi^2$ ; en général $s_n=2\pi^{n/2}/\Gamma(n/2)$.

*Démonstration.* La formule est le théorème A.18 avec $\phi=1$, $\mu(S^{n-1})=V_n$. Pour $h(r)=e^{-r^2/2}$, le membre de gauche vaut
$\big(\int_{\mathbb R}e^{-t^2/2}dt\big)^n=(2\pi)^{n/2}$ (Fubini), ce qui donne $s_n$. $I_1=\frac12\sqrt{2\pi}$, $I_2=1$, et
une intégration par parties donne $I_{n+2}=\big[-r^ne^{-r^2/2}\big]_0^\infty+n\int_0^\infty r^{n-1}e^{-r^2/2}dr=nI_n$ ; le changement de variable
$r=\sqrt{2u}$ donne $I_n=2^{n/2-1}\Gamma(n/2)$. $\square$

Pour $n=2$ et $n=3$, ce sont les coordonnées polaires et sphériques de la proposition 5.13.

**Définition.** La **loi uniforme** sur $S^{n-1}$ est $\sigma_{n-1}(B)=\lambda_n(C(B))/V_n$. Sur la sphère unité de
$\mathbb C^N\cong\mathbb R^{2N}$, c'est $\sigma_{2N-1}$.

**Théorème A.20 (décomposition polaire des lois radiales).** Soit $X$ un vecteur aléatoire de $\mathbb R^n$ ayant une densité de la
forme $f(|x|)$. Alors $X\neq0$ presque sûrement, $|X|$ a la densité $s_nr^{n-1}f(r)$ sur $(0,\infty)$, $X/|X|$ suit la loi
$\sigma_{n-1}$, et $|X|$ et $X/|X|$ sont indépendants. De plus, $\sigma_{n-1}(OB)=\sigma_{n-1}(B)$ pour tout $O\in O(n)$.
En particulier :

(a) pour un vecteur gaussien standard $g$ de $\mathbb R^n$, $g/|g|\sim\sigma_{n-1}$, indépendamment de $|g|$ ;

(b) pour un vecteur gaussien complexe standard $z\in\mathbb C^N$, $\psi=z/|z|\sim\sigma_{2N-1}$, une loi invariante par $U(N)$,
indépendante de $|z|$ ; et $|z|^2$ a la densité $t^{N-1}e^{-t}/(N-1)!$ sur $(0,\infty)$, donc $\mathbb E|z|^4=N(N+1)$.

*Démonstration.* D'après le théorème A.18 avec $h=f\cdot1_A$, $\phi=1_B$ :
$\mathbb P(|X|\in A,\ X/|X|\in B)=n\,\lambda_n(C(B))\int_Af(r)r^{n-1}dr$. Avec $A=(0,\infty)$ et $B=S^{n-1}$, cela donne
$1=nV_n\int_0^\infty f(r)r^{n-1}dr$ (et $\mathbb P(X=0)=0$). Le choix $A=(0,\infty)$ donne $\mathbb P(X/|X|\in B)=\sigma_{n-1}(B)$,
le choix $B=S^{n-1}$ donne $\mathbb P(|X|\in A)=s_n\int_Af(r)r^{n-1}dr$, et le cas général est le produit des deux.
Pour $O\in O(n)$, $OC(B)=C(OB)$ et $\lambda_n(OE)=\lambda_n(E)$ (lemme A.3(b)). (a) $\varphi_n(x)=(2\pi)^{-n/2}e^{-|x|^2/2}$ est radiale.
(b) $z=g/\sqrt2$ avec $g$ standard dans $\mathbb R^{2N}$ (proposition A.16), donc $z/|z|=g/|g|$ et $|z|=|g|/\sqrt2$ ; les matrices unitaires sont
orthogonales sur $\mathbb R^{2N}$. La densité de $|g|$ est $s_{2N}(2\pi)^{-N}r^{2N-1}e^{-r^2/2}$ avec $s_{2N}=2\pi^N/(N-1)!$, et le
changement de variable $t=r^2/2$ donne la densité de $|z|^2$ ; alors $\mathbb E|z|^4=\int_0^\infty t^{N+1}e^{-t}dt/(N-1)!=(N+1)N$. $\square$

Le lemme 7.9 est la partie (b). Pour $N=1$, $\sigma_1$ est la longueur d'arc normalisée, c'est-à-dire la mesure image de
$\frac1{2\pi}\lambda_1$ sur $[0,2\pi)$ par $\theta\mapsto e^{i\theta}$ : pour un borélien $B\subseteq S^1$, soit
$\Theta(B)=\{\theta\in(0,2\pi):e^{i\theta}\in B\}$. L'application $(r,\theta)\mapsto re^{i\theta}$ est une bijection de classe $C^1$ de l'ouvert
$(0,1)\times(0,2\pi)$ sur le disque unité ouvert privé du segment $[0,1)$, de jacobien $r$, et $C(B)$ diffère de
l'image de $(0,1)\times\Theta(B)$ par une partie de ce segment et du cercle unité, qui est négligeable. D'après la formule du changement
de variables, $\lambda_2(C(B))=\int_{\Theta(B)}\int_0^1r\,dr\,d\theta=\frac12\lambda_1(\Theta(B))$ et, avec $V_2=\pi$,
$\sigma_1(B)=\lambda_1(\Theta(B))/2\pi$ (le point $\theta=0$ est négligeable).

**Corollaire A.21 (moments d'un vecteur unitaire uniforme).** Si $u\sim\sigma_{2N-1}$ sur la sphère unité de $\mathbb C^N$, alors
$\mathbb E[u_a\bar u_b]=\delta_{ab}/N$, $\mathbb E|u_1|^4=\frac2{N(N+1)}$ et $\mathbb E|u_1|^2|u_2|^2=\frac1{N(N+1)}$ ($N\ge2$). Si
$u\sim\sigma_{n-1}$ sur $S^{n-1}\subset\mathbb R^n$, alors $\mathbb Eu_1^2=\frac1n$, $\mathbb Eu_1^4=\frac3{n(n+2)}$ et
$\mathbb Eu_1^2u_2^2=\frac1{n(n+2)}$ ($n\ge2$).

*Démonstration.* Prenons $u=z/|z|$ comme dans le théorème A.20(b). Comme $z_a\bar z_b=|z|^2u_a\bar u_b$ et que $|z|$ et $u$ sont indépendants,
$\delta_{ab}=\mathbb E[z_a\bar z_b]=\mathbb E|z|^2\,\mathbb E[u_a\bar u_b]=N\,\mathbb E[u_a\bar u_b]$ ($\mathbb E|z|^2=N$ d'après la
proposition A.16(a)). De même, $\mathbb E|z_1|^4=N(N+1)\mathbb E|u_1|^4$ et
$\mathbb E|z_1|^2|z_2|^2=N(N+1)\mathbb E|u_1|^2|u_2|^2$, où $\mathbb E|z_1|^2|z_2|^2=1$ par indépendance et
$\mathbb E|z_1|^4=\mathbb E(x^2+y^2)^2=\frac34+2\cdot\frac14+\frac34=2$ pour $x,y\sim N(0,\frac12)$ indépendants. (Pour
$Y\sim N(0,s^2)$, une intégration par parties donne $\mathbb EY^4=3s^2\,\mathbb EY^2=3s^4$, puisque
$y^3\cdot ye^{-y^2/2s^2}=-s^2y^3\frac d{dy}e^{-y^2/2s^2}$.) Dans le cas réel, $\mathbb E|g|^2=n$,
$\mathbb E|g|^4=s_n(2\pi)^{-n/2}I_{n+4}=I_{n+4}/I_n=n(n+2)$ (corollaire A.19), $\mathbb Eg_1^4=3$, $\mathbb Eg_1^2g_2^2=1$, et
$|g|^2$, $|g|^4$ sont indépendants de $u=g/|g|$ (théorème A.20(a)). $\square$

Ces valeurs concordent avec le lemme 6.14, qui les déduit de la seule invariance unitaire. Réciproquement, l'invariance détermine la loi :

**Proposition A.22 (unicité de la loi invariante).** Soit $u$ un vecteur unitaire aléatoire de $\mathbb C^N$ dont la loi est
invariante par toutes les matrices unitaires. Alors $u\sim\sigma_{2N-1}$. Pour $N=1$ : toute mesure de probabilité borélienne sur le
cercle invariante par rotation est la longueur d'arc normalisée.

*Démonstration.* Soit $Q$ de loi de Haar sur $U(N)$ (théorème 5.3) et indépendante de $u$. Pour un vecteur unitaire fixé $v$, on choisit
$W\in U(N)$ avec $We_1=v$. D'après le lemme 5.1, $QW$ suit la loi de Haar, donc $Qv=(QW)e_1$ a la loi de $Qe_1$. Dans le
théorème 5.3, $Q$ est le facteur de Gram–Schmidt d'une matrice de Ginibre (lemme 5.2), donc $Qe_1=z/|z|$ pour sa première colonne $z$, qui suit
la loi $\sigma_{2N-1}$ d'après le théorème A.20(b). En conditionnant par $u$ (Fubini, indépendance), $\mathbb P(Qu\in B)=\sigma_{2N-1}(B)$
pour tout borélien $B$. En conditionnant plutôt par $Q$, et en utilisant que $Uu$ a la loi de $u$ pour chaque matrice unitaire fixée $U$,
$\mathbb P(Qu\in B)=\mathbb P(u\in B)$. Une mesure de probabilité sur le cercle est la loi de l'application identité, et
$U(1)$ agit par rotations. $\square$

```python
N, B = 3, 400_000                                   # Theorem A.20 (b), Corollary A.21, in C^3
z = (rng.standard_normal((B, N)) + 1j * rng.standard_normal((B, N))) / math.sqrt(2)
R2 = (np.abs(z) ** 2).sum(axis=1)
u = z / np.sqrt(R2)[:, None]
t = 2.0
p_gamma = 1 - math.exp(-t) * sum(t ** k / math.factorial(k) for k in range(N))   # P(|z|^2 <= t) for Gamma(N)
mc(R2 <= t, p_gamma, "P(|z|^2 <= 2)")
mc(R2 ** 2, N * (N + 1), "E|z|^4")
mc(np.abs(u[:, 0]) ** 4, 2 / (N * (N + 1)), "E|u_1|^4")
mc(np.abs(u[:, 0] * u[:, 1]) ** 2, 1 / (N * (N + 1)), "E|u_1|^2 |u_2|^2")
mc((R2 <= t) * np.abs(u[:, 0]) ** 4, p_gamma * 2 / (N * (N + 1)), "E[1{|z|^2 <= 2} |u_1|^4], independence")
mc((u[:, 0] * u[:, 1].conj()).real, 0.0, "E Re u_1 conj(u_2)")

x3 = rng.uniform(-1, 1, (B, 3))                     # Corollary A.19: V_3 = s_3 / 3 = 4 pi / 3
mc((x3 ** 2).sum(axis=1) <= 1, math.pi / 6, "P(uniform point of [-1,1]^3 in the unit ball)")
g3 = rng.standard_normal((B, 3))
u3 = g3 / np.linalg.norm(g3, axis=1, keepdims=True)
mc(u3[:, 0] ** 4, 3 / 15, "E u_1^4 on S^2")
phase = np.angle(z[:, 0])                           # sigma_1 is normalised arc length
mc((phase >= 0) & (phase < 1), 1 / (2 * math.pi), "P(arg z_1 in [0, 1))")
```

## A.6 Moments gaussiens : la formule d'Isserlis–Wick

**Lemme A.23.** Soit $Y$ un vecteur gaussien centré de $\mathbb R^n$ de covariance $C=(c_{ij})$. Alors
$c_{ij}=\mathbb E[Y_iY_j]$, $\mathbb Ee^{\langle t,Y\rangle}=e^{\langle t,Ct\rangle/2}$ et $\mathbb Ee^{|\langle t,Y\rangle|}<\infty$ pour tout
$t\in\mathbb R^n$, et tous les moments $\mathbb E|Y_i|^k$ sont finis.

*Démonstration.* Écrivons $Y=Ag$. Alors $\mathbb E[Y_iY_j]=\sum_{k,l}A_{ik}A_{jl}\mathbb E[g_kg_l]=(AA^T)_{ij}$. Avec $a=A^Tt$,
$\langle t,Y\rangle=\sum_ka_kg_k$ et, par indépendance,
$\mathbb Ee^{\langle t,Y\rangle}=\prod_k\mathbb Ee^{a_kg_k}=\prod_ke^{a_k^2/2}=e^{|A^Tt|^2/2}=e^{\langle t,Ct\rangle/2}$. Enfin,
$e^{|s|}\le e^s+e^{-s}$ et $|Y_i|^k\le k!\,e^{|Y_i|}$. $\square$

Un **appariement** de $\{1,\dots,2p\}$ est une partition en $p$ blocs à deux éléments ; il y en a $(2p-1)!!=(2p)!/(2^pp!)$
(le partenaire de $1$ peut être choisi de $2p-1$ façons, puis on apparie les $2p-2$ éléments restants ; récurrence sur $p$).

**Théorème A.24 (Isserlis–Wick).** Soient $Y_1,\dots,Y_m$ des coordonnées d'un vecteur gaussien centré (les répétitions
sont permises). Si $m$ est impair, $\mathbb E[Y_1\cdots Y_m]=0$ ; si $m=2p$,
$$
\mathbb E[Y_1\cdots Y_{2p}]=\sum_\pi\prod_{\{a,b\}\in\pi}\mathbb E[Y_aY_b],
$$
la somme portant sur tous les appariements $\pi$ de $\{1,\dots,2p\}$. La même formule vaut pour des variables aléatoires complexes
$Z_1,\dots,Z_m$ qui sont des combinaisons linéaires complexes des coordonnées d'un même vecteur gaussien centré (par exemple les
coefficients d'une matrice de Ginibre ou du GUE et leurs conjugués complexes), avec $\mathbb E[Z_aZ_b]$ (sans conjugaison) dans les facteurs associés aux paires.

*Démonstration.* Le vecteur $(Y_1,\dots,Y_m)$ est l'image du vecteur gaussien donné par une application linéaire, donc il est gaussien
centré (théorème A.15(a)) ; soit $C$ sa covariance. Fixons $t\in\mathbb R^m$ et posons $S=\langle t,Y\rangle$. Pour $s\in\mathbb R$, les
sommes partielles de $\sum_k(sS)^k/k!$ sont majorées par $e^{|sS|}$, qui est intégrable d'après le lemme A.23, donc, par convergence
dominée et le lemme A.23,
$$
\sum_{k\ge0}\frac{s^k}{k!}\mathbb E[S^k]=\mathbb Ee^{sS}=e^{s^2q/2}=\sum_{p\ge0}\frac{s^{2p}q^p}{2^pp!},\qquad q=\langle t,Ct\rangle .
$$
Les deux membres sont des séries entières sur $\mathbb R$ au sens de la section A.2 : à droite c'est clair, et à gauche
$\sum_k\frac{\rho^k}{k!}|\mathbb E[S^k]|\le\mathbb E\sum_k\frac{\rho^k|S|^k}{k!}=\mathbb Ee^{\rho|S|}<\infty$ pour tout $\rho>0$
(convergence monotone et lemme A.23). Leurs coefficients coïncident donc (lemme A.5(c) appliqué à la différence). Par conséquent $\mathbb E\langle t,Y\rangle^m/m!$ vaut $0$ pour $m$ impair et $\langle t,Ct\rangle^p/(2^pp!)$ pour $m=2p$, pour tout
$t\in\mathbb R^m$. Les deux membres sont des polynômes en $t$, donc leurs coefficients coïncident (corollaire A.7(a) appliqué à la différence).
Dans $\mathbb E\langle t,Y\rangle^m=\sum_{i_1,\dots,i_m}t_{i_1}\cdots t_{i_m}\mathbb E[Y_{i_1}\cdots Y_{i_m}]$ ($|Y_{i_1}\cdots Y_{i_m}|\le\sum_k|Y_{i_k}|^m$, A.23), le monôme $t_1t_2\cdots t_m$
provient des $m!$ suites $(i_1,\dots,i_m)$ qui sont des permutations de $(1,\dots,m)$, donc son coefficient dans
$\mathbb E\langle t,Y\rangle^m/m!$ est $\mathbb E[Y_1\cdots Y_m]$. Pour $m$ impair, il est donc nul. Pour $m=2p$,
$q^p=\sum c_{i_1j_1}\cdots c_{i_pj_p}\,t_{i_1}t_{j_1}\cdots t_{i_p}t_{j_p}$ sur toutes les suites de $p$ couples ordonnés ; le monôme
$t_1\cdots t_{2p}$ apparaît exactement lorsque les paires non ordonnées $\{i_k,j_k\}$ forment un appariement $\pi$, et chaque appariement provient de
$p!\,2^p$ suites (ordre des paires, ordre à l'intérieur de chaque paire), chacune contribuant $\prod_{\{a,b\}\in\pi}c_{ab}$ puisque $C$ est
symétrique. En divisant par $2^pp!$, on obtient la formule, puisque $c_{ab}=\mathbb E[Y_aY_b]$.
Pour l'énoncé complexe, écrivons $Z_a=\sum_i\alpha_{ai}Y_i$. Les deux membres sont multilinéaires en $(Z_1,\dots,Z_m)$ — celui de gauche par
linéarité de $\mathbb E$, celui de droite parce que chaque $Z_a$ apparaît dans exactement un facteur de chaque terme — et ils coïncident lorsque chaque
$Z_a$ est l'un des $Y_i$ ; il suffit de développer. $\square$

Par exemple, $\mathbb E|z|^4=\mathbb E[z\bar zz\bar z]=1\cdot1+0+1\cdot1=2$ pour une gaussienne complexe standard $z$. La proposition 1.10
et le lemme 7.10 utilisent le cas $m=4$ ; le théorème A.24 le démontre sans dériver sous l'espérance.

```python
def pairings(idx):
    """All pairings of the list idx, as lists of pairs."""
    if not idx:
        yield []
        return
    for i in range(1, len(idx)):
        for rest in pairings(idx[1:i] + idx[i + 1:]):
            yield [(idx[0], idx[i])] + rest


def wick(cov, idx):                                       # sum over pairings of the positions of idx
    return sum(math.prod(cov[idx[a], idx[b]] for a, b in p) for p in pairings(list(range(len(idx)))))


assert [len(list(pairings(list(range(2 * p))))) for p in range(1, 5)] == [1, 3, 15, 105]    # (2p-1)!!
Amat = rng.standard_normal((3, 3)) / math.sqrt(3)
cov = Amat @ Amat.T
Y = rng.standard_normal((1_000_000, 3)) @ Amat.T
for idx in ((0, 1, 2, 0), (0, 0, 0, 0), (0, 1, 1, 2, 2, 0)):
    mc(np.prod(Y[:, list(idx)], axis=1), wick(cov, idx), f"E Y{idx}")
mc(np.prod(Y[:, [0, 1, 2]], axis=1), 0.0, "E Y(0, 1, 2), odd")
M = np.array([[1, 1j], [1, -1j]])                  # (z, conj z) = M (x, y) with x, y independent N(0, 1/2)
cz = M @ M.T / 2                                     # E[Z_a Z_b] = [[0, 1], [1, 0]]
assert abs(wick(cz, (0, 1, 0, 1)) - 2) < 1e-12                           # E|z|^4 = 2
```

## A.7 Variables exponentielles et espacements uniformes

Une variable aléatoire $e$ est **exponentielle** ($e\sim\operatorname{Exp}(1)$) si elle a la densité $e^{-y}$ sur $(0,\infty)$ ;
$e/\beta$ est exponentielle de paramètre $\beta$. Soit $\Delta_n=\{x\in\mathbb R^n:0<x_1<\dots<x_n<1\}$.

**Théorème A.25.** Soient $e_0,\dots,e_n$ ($n\ge1$) des variables indépendantes de loi $\operatorname{Exp}(1)$, $T=e_0+\dots+e_n$ et
$X_k=(e_0+\dots+e_{k-1})/T$ pour $k=1,\dots,n$. Alors $(X_1,\dots,X_n)$ a la densité $n!$ sur $\Delta_n$, $T$ a la densité
$t^ne^{-t}/n!$ sur $(0,\infty)$, et $(X_1,\dots,X_n)$ et $T$ sont indépendants. En particulier,
$(e_0,\dots,e_n)/T$ est indépendant de $T$.

*Démonstration.* L'application $\Phi(x,t)=t\,(x_1,\ x_2-x_1,\ \dots,\ x_n-x_{n-1},\ 1-x_n)$ est une bijection de l'ouvert
$\Delta_n\times(0,\infty)$ sur $(0,\infty)^{n+1}$, de réciproque $y\mapsto\big(((y_0+\dots+y_{k-1})/\sum_iy_i)_{k=1}^n,\ \sum_iy_i\big)$ ;
toutes deux sont de classe $C^1$. Dans la matrice jacobienne de $\Phi$, la colonne de $\partial/\partial t$ a pour somme des coefficients $1$ et chaque colonne
$\partial/\partial x_j$ (coefficients $t$ et $-t$) a pour somme des coefficients $0$. En remplaçant la dernière ligne par la somme de toutes les lignes (ce qui ne
change pas le déterminant), on obtient la ligne $(0,\dots,0,1)$, et le bloc $n\times n$ restant est triangulaire inférieur de
diagonale $t$ ; donc $|\det D\Phi|=t^n$. Le vecteur $e=(e_0,\dots,e_n)$ a la densité $e^{-\sum_iy_i}$, et
$\sum_i\Phi(x,t)_i=t$. D'après la formule du changement de variables, pour tout borélien $E\subseteq\Delta_n\times(0,\infty)$,
$$
\mathbb P\big((X,T)\in E\big)=\mathbb P\big(e\in\Phi(E)\big)=\int_{\Phi(E)}e^{-\sum_iy_i}dy=\int_Ee^{-t}t^n\,dx\,dt .
$$
La densité $e^{-t}t^n=n!\cdot\frac{t^ne^{-t}}{n!}$ est un produit, et $\int_0^\infty t^ne^{-t}dt=n!$, donc
$\lambda_n(\Delta_n)=1/n!$ ; d'après Fubini, $\mathbb P(X\in A,T\in A')=\mathbb P(X\in A)\mathbb P(T\in A')$. Enfin, $(e_0,\dots,e_n)/T$ est
$(X_1,X_2-X_1,\dots,1-X_n)$, une fonction de $X$. $\square$

**Corollaire A.26 (espacements uniformes).** Soient $U_1,\dots,U_n$ indépendantes et uniformes sur $(0,1)$, de statistiques d'ordre
$U_{(1)}<\dots<U_{(n)}$ et d'espacements $s_0=U_{(1)}$, $s_k=U_{(k+1)}-U_{(k)}$ ($1\le k\le n-1$), $s_n=1-U_{(n)}$. Alors :

(a) les ex æquo ont une probabilité nulle, et $(U_{(1)},\dots,U_{(n)})$ a la densité $n!$ sur $\Delta_n$ ;

(b) $(s_0,\dots,s_n)$ a la loi de $(e_0,\dots,e_n)/T$, avec $e_i$, $T$ comme dans le théorème A.25 ;

(c) $\mathbb Es_i=\frac1{n+1}$, $\mathbb Es_i^2=\frac2{(n+1)(n+2)}$ et $\mathbb E[s_is_j]=\frac1{(n+1)(n+2)}$ pour $i\neq j$ ;

(d) pour $i\neq j$, $s_i/s_j$ a la loi de $e_i/e_j$, avec $\mathbb P(e_i/e_j\le t)=\frac t{1+t}$ ($t>0$) ; donc, pour
$0\le i\le n-1$, le rapport d'espacements $\min(s_i,s_{i+1})/\max(s_i,s_{i+1})$ a la loi de $\min(e_i,e_{i+1})/\max(e_i,e_{i+1})$
obtenue à la proposition 5.14(a), $\mathbb P(r\le t)=\frac{2t}{1+t}$ sur $[0,1]$, de moyenne $2\ln2-1$.

Pour $n$ points uniformes indépendants sur un intervalle $[a,b]$, l'application affine $u\mapsto a+(b-a)u$ multiplie tous les espacements
par $b-a$ : les espacements divisés par $b-a$ vérifient (b) et (c), et les rapports de (d) sont inchangés.

*Démonstration.* (a) L'ensemble où deux coordonnées coïncident est une réunion finie d'hyperplans, donc négligeable (corollaire A.7(c)). En dehors
de cet ensemble, il existe exactement une permutation $\pi$ telle que $U_{\pi(1)}<\dots<U_{\pi(n)}$, donc, pour tout borélien $E\subseteq\Delta_n$, l'événement
$\{(U_{(1)},\dots,U_{(n)})\in E\}$ est, à un ensemble négligeable près, la réunion disjointe sur $\pi$ des événements
$\{(U_{\pi(1)},\dots,U_{\pi(n)})\in E\}$. Chacun a la probabilité $\lambda_n(E)$, car $(U_1,\dots,U_n)$ a la densité
$1_{(0,1)^n}$ et permuter les coordonnées préserve à la fois $\lambda_n$ et $(0,1)^n$. La probabilité vaut donc $n!\,\lambda_n(E)$.
(b) D'après (a) et le théorème A.25, les vecteurs $(U_{(k)})_k$ et $(X_k)_k$ ont la même loi, et les espacements sont la même fonction
linéaire de l'un comme de l'autre. (c) Soit $w=(e_0,\dots,e_n)/T$, indépendant de $T$ (théorème A.25), avec $\mathbb ET=n+1$ et
$\mathbb ET^2=\int_0^\infty t^{n+2}e^{-t}dt/n!=(n+1)(n+2)$. De $e_i=Tw_i$ on tire $1=\mathbb Ee_i=\mathbb ET\,\mathbb Ew_i$,
$2=\mathbb Ee_i^2=\mathbb ET^2\,\mathbb Ew_i^2$ et $1=\mathbb E[e_ie_j]=\mathbb ET^2\,\mathbb E[w_iw_j]$. (d) $s_i/s_j=w_i/w_j=e_i/e_j$ et, d'après
Tonelli, $\mathbb P(e_i\le te_j)=\int_0^\infty e^{-y}(1-e^{-ty})\,dy=1-\frac1{1+t}$. Le rapport d'espacements est une fonction de
$s_i/s_{i+1}$ (à savoir $\min(x,1/x)$ en $x=s_i/s_{i+1}$, définie presque sûrement puisque tous les $e_k>0$), donc il a la loi de la
même fonction de $e_i/e_{i+1}$, qui est le rapport $r$ de la proposition 5.14(a) pour les variables indépendantes de loi $\operatorname{Exp}(1)$
$e_i,e_{i+1}$. $\square$

En particulier, considérons les « niveaux de Poisson » $\lambda_k=e_0+\dots+e_{k-1}$, $k=1,\dots,n+1$ (espacements exponentiels
indépendants $e_0,\dots,e_n$), et $n$ points uniformes triés avec leurs $n+1$ espacements $s_0,\dots,s_n$. Le vecteur des
rapports d'espacements consécutifs $r_i=\min(d_i,d_{i+1})/\max(d_i,d_{i+1})$, $0\le i<n$, a la même loi jointe pour
$d=(e_0,\dots,e_n)$ et pour $d=(s_0,\dots,s_n)$, car les rapports ne changent pas lorsque tous les $d_i$ sont divisés par $T$, et
$(e_0,\dots,e_n)/T$ a la loi de $(s_0,\dots,s_n)$ (partie (b)) ; et, d'après le théorème A.25, les $n$ premiers niveaux de Poisson divisés par le $(n+1)$-ième ont exactement la loi de $n$ points
uniformes triés.

```python
n, B = 5, 200_000                                       # Corollary A.26
levels = np.sort(rng.uniform(0, 1, (B, n)), axis=1)
s = np.diff(np.concatenate([np.zeros((B, 1)), levels, np.ones((B, 1))], axis=1), axis=1)    # s_0, ..., s_n
e = rng.exponential(1.0, (B, n + 1))
w = e / e.sum(axis=1, keepdims=True)                    # normalised exponentials
for name, sp in (("uniform spacings", s), ("normalised exponentials", w)):
    mc(sp[:, 2], 1 / (n + 1), f"{name}: E s_2")
    mc(sp[:, 2] ** 2, 2 / ((n + 1) * (n + 2)), f"{name}: E s_2^2")
    mc(sp[:, 0] * sp[:, n], 1 / ((n + 1) * (n + 2)), f"{name}: E s_0 s_n")
    for tt in (0.5, 2.0):
        mc(sp[:, 1] / sp[:, 2] <= tt, tt / (1 + tt), f"{name}: P(s_1 / s_2 <= {tt})")
r = np.minimum(s[:, 1:-2], s[:, 2:-1]) / np.maximum(s[:, 1:-2], s[:, 2:-1])    # ratios of interior spacings
mc(r.mean(axis=1), 2 * math.log(2) - 1, "<r>, uniform levels")
T = e.sum(axis=1)
mc(w[:, 0] * T ** 2, n + 2, "E[w_0 T^2] = E w_0 E T^2 (independence)")
```

**Retour à la** [table des matières](README.md).
