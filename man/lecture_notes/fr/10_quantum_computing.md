[English](../10_quantum_computing.md) · [Español](../es/10_quantum_computing.md) · **Français** · [Deutsch](../de/10_quantum_computing.md) · [中文](../zh/10_quantum_computing.md) · [日本語](../ja/10_quantum_computing.md) · [Русский](../ru/10_quantum_computing.md)

# 10. Calcul quantique

Le bruit sur un ordinateur quantique est un canal, et un canal engendre une algèbre : la
$*$-algèbre de ses opérateurs de Kraus. Sa structure (chapitre 2) indique quelle information survit
intacte au bruit ; une condition sur les produits d'opérateurs de Kraus, faisant intervenir une
autre algèbre, indique quelle information une récupération peut restaurer. La dernière partie
traite des algorithmes qui conservent une charge : leurs états vivent dans l'algèbre des secteurs
de charge, et l'amortissement d'amplitude les fait passer d'un secteur à l'autre.

**Vous apprendrez** pourquoi un canal unifère et son dual ont les mêmes points fixes (le commutant
du bruit), le théorème de Knill–Laflamme avec une démonstration complète, comment ne protéger que
l'étiquette d'un secteur, et à quoi ressemblent le mélangeur XY de l'ansatz à opérateurs alternés
quantique, l'amortissement d'amplitude et la post-sélection dans l'algèbre $\bigoplus_wM_{\binom nw}$
des secteurs de poids de Hamming.
Les vérifications numériques de ce chapitre se trouvent dans le
[Tutoriel 10](../../tutorial/fr/10_quantum_computing.md).

Dans tout ce chapitre, un **canal** sur $M_D(\mathbb C)$ est une application $\Phi(\rho)=\sum_{i=1}^rE_i\rho E_i^*$
dont les **opérateurs de Kraus** vérifient $\sum_iE_i^*E_i=1$ (préservation de la trace) ; il est
**unifère** si de plus $\sum_iE_iE_i^*=1$. Son **dual** $\Phi^*(T)=\sum_iE_i^*TE_i$ vérifie
$\mathrm{Tr}(\Phi(\rho)T)=\mathrm{Tr}(\rho\,\Phi^*(T))$ (section 6.2). (D'après le théorème
6.2, toute application complètement positive est de cette forme ; ici les canaux sont *donnés* par
des opérateurs de Kraus.) L'**algèbre d'interaction** $\mathcal A$ de $\Phi$ est la $*$-algèbre
engendrée par $1$ et les $E_i$, l'espace engendré par tous les produits des $E_i$ et des $E_i^*$ ;
une matrice commute avec $\mathcal A$ si et seulement si elle commute avec tous les $E_i$ et $E_i^*$.
Un état $\rho$ est **à support** dans un sous-espace $K$ de projecteur $P_K$ si $\rho=P_K\rho P_K$. Sur
$n$ qubits, $D=2^n$ ; $X,Y,Z$ sont les matrices de Pauli, $Z=\mathrm{diag}(1,-1)$, $X_j$ est $X$ sur le qubit $j$ (numérotés $1,\dots,n$ ;
la section 10.4 suit la bibliothèque, sites $0,\dots,n-1$), $|b\rangle$, $b\in\{0,1\}^n$, la base
de calcul, $|b|$ le **poids de Hamming**, et $P_w$ le projecteur sur $H_w=\mathrm{span}\{|b\rangle:|b|=w\}$.

<a id="s-10-1"></a>

## 10.1 Algèbres de bruit, sous-systèmes sans bruit et points fixes

D'après le théorème 2.6, il existe un unitaire $V$ tel que $V^*\mathcal AV=\bigoplus_cM_{k_c}\otimes1_{m_c}$
et $V^*\mathcal A'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$ : chaque opérateur de Kraus n'agit que sur le
premier facteur tensoriel de chaque terme de la somme.

<a id="r-10-1"></a>

**Proposition 10.1 (sous-systèmes sans bruit ; sous-espaces sans décohérence).** (a) Écrivons
$V^*E_iV=\bigoplus_cE_{i,c}\otimes1_{m_c}$, et soit $\rho=V(\sigma\otimes\tau)V^*$, placé dans le terme
$c$, pour des états $\sigma$ sur $\mathbb C^{k_c}$ et $\tau$ sur $\mathbb C^{m_c}$. Alors
$\Phi(\rho)=V(\Phi_c(\sigma)\otimes\tau)V^*$ avec le canal $\Phi_c(\sigma)=\sum_iE_{i,c}\sigma E_{i,c}^*$ sur
$M_{k_c}$ : l'état $\tau$ n'est pas touché. De plus, $\Phi^*(T)=T$ pour tout $T\in\mathcal A'$.
(b) Si $K$ est un sous-espace tel que $E_iP_K=e_iP_K$, $e_i\in\mathbb C$, pour tout $i$, alors
$\Phi(\rho)=\rho$ pour tout état $\rho$ à support dans $K$.

*Démonstration.* (a) $\Phi(\rho)=V\sum_i(E_{i,c}\sigma E_{i,c}^*\otimes\tau)V^*$, et $\Phi_c$ préserve la
trace car le terme $c$ de $V^*(\sum_iE_i^*E_i)V=1$ est $(\sum_iE_{i,c}^*E_{i,c})\otimes1_{m_c}$. Pour
$T\in\mathcal A'$, $\Phi^*(T)=T\sum_iE_i^*E_i=T$. (b) $P_KE_i^*=(E_iP_K)^*=\bar e_iP_K$, donc
$E_i\rho E_i^*=E_iP_K\rho P_KE_i^*=|e_i|^2\rho$, et $\sum_i|e_i|^2P_K=P_K\sum_iE_i^*E_iP_K=P_K$. $\square$

Un terme avec $k_c=1$ relève du cas (b) : c'est un **sous-espace sans décohérence** (Zanardi–Rasetti
1997, Lidar–Chuang–Whaley 1998). Pour $k_c>1$, le facteur $\mathbb C^{m_c}$ est un **sous-système
sans bruit** (Knill–Laflamme–Viola 2000) : le bruit agit, mais seulement sur l'autre facteur. La
proposition 10.1 vaut pour tout canal ; pour les canaux unifères, les observables préservées sont
*exactement* celles de $\mathcal A'$ :

<a id="r-10-2"></a>

**Théorème 10.2 (points fixes des canaux unifères).** Soit $\Phi$ unifère. Alors
$$
\mathrm{Fix}(\Phi)=\{T\in M_D(\mathbb C):\Phi(T)=T\}=\mathcal A'=\mathrm{Fix}(\Phi^*).
$$
En particulier, $\mathrm{Fix}(\Phi)$ est une algèbre de von Neumann, $*$-isomorphe à $\bigoplus_cM_{m_c}$.

*Démonstration.* La première égalité est le théorème 6.10 pour l'algèbre $M_D(\mathbb C)$ (un seul
bloc) : $T$ est un point fixe si et seulement s'il commute avec tous les $E_i$ et $E_i^*$,
c'est-à-dire si et seulement si $T\in\mathcal A'$. Le dual $\Phi^*$ est encore un canal unifère : ses
opérateurs de Kraus $E_i^*$ vérifient $\sum_iE_iE_i^*=1$ (préservation de la trace) et
$\sum_iE_i^*E_i=1$ (unifère), et ils engendrent la même $*$-algèbre $\mathcal A$. Le théorème 6.10
appliqué à $\Phi^*$ donne donc $\mathrm{Fix}(\Phi^*)=\mathcal A'$. La dernière assertion est le
théorème 2.6. $\square$

*Calcul :* $\mathrm{Fix}(\Phi)=\mathcal A'$ pour le déphasage collectif, les commutants et les points fixes étant calculés comme des noyaux — [tutoriel, section 10.1](../../tutorial/fr/10_quantum_computing.md#s-10-1).

Pour tout canal, $\mathcal A'\subseteq\mathrm{Fix}(\Phi^*)$ (proposition 10.1(a)), mais sans
l'hypothèse unifère le théorème est en défaut : pour l'amortissement d'amplitude,
$\mathrm{Fix}(\Phi)=\mathbb C|0\rangle\langle0|$, alors que $\mathcal A'=\mathrm{Fix}(\Phi^*)=\mathbb C1$
(section 6.4 et exercice 1). (Historiquement : Kribs 2003, Arias–Gheondea–Gudder 2002.)

<a id="r-10-3"></a>

**Corollaire 10.3 (déphasage collectif).** Soient $U=e^{-i\theta Z_{\mathrm{tot}}}$,
$Z_{\mathrm{tot}}=\sum_jZ_j$, avec $e^{2i\theta d}\neq1$ pour $d=1,\dots,n$, et
$\Phi(\rho)=\tfrac12(U\rho U^*+U^*\rho U)$. Alors $\mathcal A=\mathrm{span}\{P_0,\dots,P_n\}$, chaque
$H_w$ est un sous-espace sans décohérence, et $\mathrm{Fix}(\Phi)=\bigoplus_wB(H_w)\cong\bigoplus_wM_{\binom nw}$,
de dimension $\sum_w\binom nw^2$.

*Démonstration.* $Z_{\mathrm{tot}}|b\rangle=(n-2|b|)|b\rangle$, donc $U=\sum_we^{-i\theta(n-2w)}P_w$, de
valeurs propres distinctes (leurs quotients sont les $e^{2i\theta(w-w')}$, $0<|w-w'|\le n$). Les
polynômes de Lagrange de la démonstration du théorème 1.5(a) donnent $P_w=p_w(U)\in\mathcal A$ ;
réciproquement, $U,U^*$ appartiennent à l'algèbre $\mathrm{span}\{P_w\}$. Donc
$\mathcal A=\mathrm{span}\{P_w\}$ ($k_c=1$, $m_c=\binom nw$), et la proposition 10.1(b) s'applique à
chaque $H_w$. Si $T$ commute avec tous les $P_w$, alors $P_vTP_{v'}=TP_vP_{v'}=0$ pour $v\neq v'$,
c'est-à-dire que $T$ est diagonal par blocs, et réciproquement ; le théorème 10.2 donne le reste.
$\square$

*Calcul :* les dimensions $5$, $70$, $70$ pour $n=4$ — [tutoriel, section 10.1](../../tutorial/fr/10_quantum_computing.md#s-10-1). *Dans la bibliothèque :* $\bigoplus_wM_{\binom nw}$ est `SpinChain(n).algebra`, sur laquelle ce canal est l'identité — [tutoriel, section 10.1](../../tutorial/fr/10_quantum_computing.md#s-10-1).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Un champ fluctuant, le même sur tous les qubits, décale la phase de
> $|b\rangle$ proportionnellement à $n-2|b|$. Les états d'un même poids de Hamming acquièrent une
> phase commune, de sorte que $|01\rangle,|10\rangle$ codent un qubit insensible à ce bruit (ions
> piégés : Kielpinski et al., *Science* 291, 2001), tandis que les cohérences entre poids
> s'amortissent : le bruit impose une règle de supersélection (section 3.3).

**Décohérence collective.** Nous disons qu'un canal est **collectif** si ses opérateurs de Kraus sont des
multiples de rotations $e^{-i\theta\,\mathbf n\cdot\mathbf S}$ ($\theta\in\mathbb R$, $\mathbf n\in\mathbb R^3$ un vecteur unitaire),
$S_\alpha=\tfrac12\sum_j\sigma_{\alpha,j}$, $\sigma\in\{X,Y,Z\}$ ; cela modélise un champ qui fluctue aussi en direction
(le modèle, pas un théorème). Comme $\mathbf n\cdot\mathbf S$ est autoadjoint, $e^{-i\theta\,\mathbf n\cdot\mathbf S}$ en est un
polynôme (démonstration du théorème 1.5(a)), donc ces opérateurs de Kraus appartiennent à la $*$-algèbre $\mathcal B$
engendrée par $1,S_x,S_y,S_z$.

<a id="r-10-4"></a>

**Proposition 10.4 (sous-espace singulet).** Soit $K=\ker S_x\cap\ker S_y\cap\ker S_z$. Tout canal
dont les opérateurs de Kraus appartiennent à $\mathcal B$ fixe tout état à support dans $K$.

*Démonstration.* $\mathcal B$ est engendrée linéairement par les mots en $S_x,S_y,S_z$ (et $1$) ; les
mots de longueur $\ge1$ annulent $K$, donc chaque élément de $\mathcal B$ agit sur $K$ comme un
scalaire. On applique la proposition 10.1(b). $\square$

Pour $n=4$, la proposition 4.13 ([chapitre 4](04_symmetries.md), cas $L=4$) donne
$\mathcal B\cong M_5\oplus M_3\oplus M_1$ avec les multiplicités $1,3,2$ (spin total $2,1,0$), donc $\dim\mathcal B=35$ et
$\mathcal B'\cong M_1\oplus M_3\oplus M_2$ est de dimension $14$. De plus $\dim K=m_0=2$ : les $S_a$ agissent terme par
terme, donc $K$ est la somme de ses intersections avec les termes ; sur le terme de spin $S$, ils agissent comme
$S^{(S)}_a\otimes1$ (corollaire 4.12), et
$\sum_bu_b\otimes f_b$ ($f_b$ une base de $\mathbb C^{m_S}$) est annulé par tous ces opérateurs si et seulement si chaque $u_b$
est annulé par tous les $S^{(S)}_a$. Ces $u$ forment un sous-espace invariant du spin $S$ irréductible (théorème 4.11), qui
n'est pas l'espace entier pour $S>0$ (car $S^{(S)}_z\neq0$), donc est $0$, tandis que pour $S=0$ c'est $\mathbb C^1$. La proposition 10.1(a) et sa démonstration valent mot pour mot pour toute $*$-algèbre de la forme
$V(\bigoplus_cM_{k_c}\otimes1_{m_c})V^*$ qui contient les opérateurs de Kraus, ici $\mathcal B$ : quatre qubits portent un
qubit sans décohérence ($k=1$, $m=2$) et un qutrit sans bruit ($k=3$, $m=3$) pour tout canal dont les opérateurs de Kraus
appartiennent à $\mathcal B$.

*Calcul :* $\dim\mathcal B'=14$, $\dim\mathcal B=35$, $\dim K=2$, et la proposition 10.4 pour des rotations collectives aléatoires — [tutoriel, section 10.1](../../tutorial/fr/10_quantum_computing.md#s-10-1).

<a id="s-10-2"></a>

## 10.2 Correction d'erreurs : le théorème de Knill–Laflamme

Un **code** est un sous-espace non nul $\mathcal C\subseteq\mathbb C^D$ de projecteur $P$. Il **corrige**
$\Phi$ s'il existe un canal $\mathcal R$ (la **récupération**) tel que $\mathcal R(\Phi(\rho))=\rho$ pour
tout état $\rho$ à support dans $\mathcal C$.

<a id="r-10-5"></a>

**Théorème 10.5 (Knill–Laflamme).** $\mathcal C$ corrige $\Phi$ si et seulement s'il existe des
nombres $\lambda_{ij}$ tels que
$$
PE_i^*E_jP=\lambda_{ij}P\qquad\text{for all }i,j.
$$

*Démonstration.* *Condition nécessaire.* Soient $R_m$ les opérateurs de Kraus de $\mathcal R$ ;
alors $\mathcal R\circ\Phi$ a pour opérateurs de Kraus les $R_mE_i$. Pour un vecteur unitaire
$\psi\in\mathcal C$, $\sum_{m,i}R_mE_i|\psi\rangle\langle\psi|E_i^*R_m^*=|\psi\rangle\langle\psi|$, et pour
$\phi\perp\psi$ la valeur moyenne en $\phi$ donne $\sum_{m,i}|\langle\phi,R_mE_i\psi\rangle|^2=0$, donc
$R_mE_i\psi\in\mathbb C\psi$. Ainsi $R_mE_iP=c_{mi}P$ : une application linéaire $T$ telle que
$Tx\in\mathbb Cx$ pour tout $x$ est un multiple de $1$ (si $Tx=ax$, $Ty=by$, $T(x+y)=c(x+y)$, avec $x,y$
indépendants, alors $a=c=b$). Comme $\sum_mR_m^*R_m=1$,
$PE_i^*E_jP=\sum_m(R_mE_iP)^*(R_mE_jP)=\sum_m\bar c_{mi}c_{mj}P$.

*Condition suffisante.* $\Lambda=(\lambda_{ij})$ est hermitienne, car
$\lambda_{ji}P=(PE_i^*E_jP)^*=\bar\lambda_{ij}P$, et positive, car
$\sum_{ij}\bar c_ic_j\lambda_{ij}P=PF^*FP\ge0$ pour $F=\sum_jc_jE_j$. Écrivons $\Lambda=UDU^*$, avec
$U$ unitaire, $D=\mathrm{diag}(d_k)$, $d_k\ge0$ (théorème 1.3), et $F_k=\sum_iu_{ik}E_i$. Alors
$\sum_kF_k\rho F_k^*=\sum_{ij}(UU^*)_{ij}E_i\rho E_j^*=\Phi(\rho)$ et $PF_k^*F_lP=(U^*\Lambda U)_{kl}P=d_k\delta_{kl}P$.
Si $d_k=0$, alors $\|F_kPx\|^2=0$ pour tout $x$, donc $F_kP=0$. Si $d_k>0$, posons
$V_k=F_kP/\sqrt{d_k}$ ; alors $V_k^*V_l=\delta_{kl}P$, donc $V_kP=V_k$, les $Q_k=V_kV_k^*$ sont des
projecteurs deux à deux orthogonaux et $Q=\sum_kQ_k$ est un projecteur. La récupération
$$
\mathcal R(\sigma)=\sum_{k:d_k>0}V_k^*\sigma V_k+(1-Q)\sigma(1-Q)
$$
est un canal, puisque $\sum_kV_kV_k^*+(1-Q)^2=1$. Pour $\rho=P\rho P$, $\Phi(\rho)=\sum_{d_k>0}d_kV_k\rho V_k^*$ ;
comme $V_l^*V_k=\delta_{kl}P$ et $(1-Q)V_k=(1-Q)Q_kV_k=0$, $\mathcal R(\Phi(\rho))=\sum_kd_k\rho$, et
$\sum_kd_k=\mathrm{Tr}\,\Lambda=1$ car $\sum_i\lambda_{ii}P=P\sum_iE_i^*E_iP=P$. $\square$

(Knill–Laflamme 1997 ; Nielsen–Chuang, théorème 10.1.) $\mathcal R$ mesure le **syndrome** — lequel
des sous-espaces orthogonaux $Q_k\mathbb C^D$ contient l'état — et ramène l'état dans le code par une
rotation.

**Exemple : le code de bit-flip.** $\mathcal C=\mathrm{span}\{|000\rangle,|111\rangle\}$ et « au plus
un bit inversé », $E_0=\sqrt{p_0}\,1$, $E_j=\sqrt{p_j}\,X_j$. Pour $i\neq j$, $X_iX_j$ (avec $X_0=1$)
inverse un ou deux bits et envoie $\mathcal C$ dans $\mathcal C^\perp$ ; donc
$\Lambda=\mathrm{diag}(p_0,\dots,p_3)$, $V_k=X_kP$, $Q=1$ et $\mathcal R(\sigma)=\sum_{k=0}^3PX_k\sigma X_kP$.
Si chaque qubit est inversé indépendamment avec probabilité $p$, la condition est en défaut (les
opérateurs de Kraus comprennent des multiples de $X_1X_2$ et de $X_3$, et
$PX_1X_2X_3P=|000\rangle\langle111|+|111\rangle\langle000|$ n'est pas un multiple de $P$), et $\mathcal R$ envoie $|0_L\rangle=|000\rangle$ sur $|1_L\rangle$ dès que $\ge2$ inversions se produisent et le
laisse en place sinon : la fidélité de $|0_L\rangle$ après récupération est la probabilité d'au plus
une inversion, $(1-p)^3+3p(1-p)^2=1-3p^2+2p^3=0.972$ pour $p=0.1$.

*Calcul :* le théorème 10.5 sur cet exemple : la matrice $\Lambda$, la récupération $\mathcal R$ et la fidélité $0.972$ — [tutoriel, section 10.2](../../tutorial/fr/10_quantum_computing.md#s-10-2). *Dans la bibliothèque :* les canaux donnés par des opérateurs de Kraus, `channels.kraus_channel`, et leur composée `R @ Phi` — [tutoriel, section 10.2](../../tutorial/fr/10_quantum_computing.md#s-10-2).

<a id="s-10-3"></a>

## 10.3 Correction d'erreurs par algèbres d'opérateurs

Le point de vue de Heisenberg permet de protéger moins que le code entier. Soit $\mathfrak A$ une
algèbre de von Neumann sur $\mathcal C=P\mathbb C^D$ (une $*$-algèbre d'opérateurs sur $\mathcal C$
contenant $P$), les observables à protéger. $\mathfrak A$ est **corrigeable** pour $\Phi$ s'il existe
un canal $\mathcal R$ tel que
$$
P\,(\mathcal R\circ\Phi)^*(A)\,P=A\qquad\text{for all }A\in\mathfrak A,
$$
c'est-à-dire $\mathrm{Tr}(\mathcal R(\Phi(\rho))A)=\mathrm{Tr}(\rho A)$ pour tous les états
$\rho$ à support dans $\mathcal C$ et tout $A\in\mathfrak A$ (le membre de gauche vaut
$\mathrm{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(A)P)$, et un opérateur $Y\in B(\mathcal C)$ tel que $\mathrm{Tr}(\rho Y)=0$
pour tous ces $\rho$ est nul : $\langle x,Yx\rangle=\langle Px,YPx\rangle=0$ pour tout $x$, et polarisation comme dans la
démonstration du corollaire 1.4). Pour $\mathfrak A=B(\mathcal C)$, c'est la correction de la section 10.2 : si
$\sigma=\mathcal R(\Phi(\rho))$ vérifie $P\sigma P=\rho$ (prendre tous les $A\in B(\mathcal C)$), alors, comme
$\mathcal R\circ\Phi$ préserve la trace,
$\mathrm{Tr}(1-P)\sigma(1-P)=\mathrm{Tr}\,\sigma-\mathrm{Tr}\,P\sigma P=1-1=0$, donc $\sigma^{1/2}(1-P)=0$ et
$\sigma=\rho$ ; réciproquement, $\mathcal R(\Phi(\rho))=\rho$ donne l'égalité des traces. L'extrême opposé est une
algèbre commutative : seule une étiquette classique est protégée.

<a id="r-10-6"></a>

**Proposition 10.6 (correction de l'étiquette d'un secteur).** Soient $P=\sum_{c=1}^CP_c$, avec des
projecteurs $P_c\neq0$ deux à deux orthogonaux, et $\mathfrak A=\mathrm{span}\{P_1,\dots,P_C\}$. Les
assertions suivantes sont équivalentes :
1. il existe $N_1,\dots,N_C\ge0$, $\sum_cN_c=1$, tels que $\mathrm{Tr}(N_c\Phi(\rho))=1$ pour tout $c$
   et tout état $\rho$ à support dans $P_c\mathbb C^D$ ;
2. $P_cE_i^*E_jP_d=0$ pour tous $i,j$ et $c\neq d$, c'est-à-dire $PE_i^*E_jP\in\mathfrak A'$, le
   commutant de $\mathfrak A$ dans $B(\mathcal C)$ ;
3. $\mathfrak A$ est corrigeable pour $\Phi$.

*Démonstration.* $T\in B(\mathcal C)$ commute avec tous les $P_c$ si et seulement si $P_cTP_d=0$ pour
$c\neq d$ (comme dans le corollaire 10.3).

(2)$\Rightarrow$(1). Les sous-espaces $S_c=\mathrm{span}\{E_iP_cx\}$ sont deux à deux orthogonaux,
car $\langle E_iP_cx,E_jP_dy\rangle=\langle x,P_cE_i^*E_jP_dy\rangle=0$. Soient $N_c$ le projecteur sur
$S_c$ pour $c<C$ et $N_C=1-\sum_{c<C}N_c$, qui est aussi un projecteur. Pour $\rho$ à support dans
$P_c\mathbb C^D$, $\Phi(\rho)=\sum_iE_iP_c\rho P_cE_i^*$ a son image dans $S_c$, et $N_c$ est l'identité
sur $S_c$ (pour $c=C$ parce que $S_C\perp S_{c'}$, $c'<C$), donc
$\mathrm{Tr}(N_c\Phi(\rho))=\mathrm{Tr}\,\Phi(\rho)=1$.
(1)$\Rightarrow$(2). Pour un vecteur unitaire $x\in P_c\mathbb C^D$ et $d\neq c$, les nombres
$\mathrm{Tr}(N_{d'}\Phi(|x\rangle\langle x|))\ge0$ ont pour somme $1$ et celui pour $d'=c$ vaut $1$,
donc $0=\mathrm{Tr}(N_d\Phi(|x\rangle\langle x|))=\sum_i\|N_d^{1/2}E_ix\|^2$ ; ainsi $N_dE_iP_c=0$
et, en sommant sur $d$, $N_cE_iP_c=E_iP_c$. Par conséquent $P_cE_i^*E_jP_d=P_cE_i^*N_cE_jP_d=0$ pour
$c\neq d$.
(1),(2)$\Rightarrow$(3). Soient $\omega_c=P_c/\mathrm{Tr}\,P_c$ et $\mathcal R(\sigma)=\sum_c\mathrm{Tr}(N_c\sigma)\,\omega_c$
(« mesurer, puis préparer »), un canal ayant pour opérateurs de Kraus les
$\sqrt\mu\,|u\rangle\langle f|N_c^{1/2}$ pour les couples propres $(\mu,u)$ de $\omega_c$ et les vecteurs
$f$ d'une base orthonormée. Comme $\mathrm{Tr}(\mathcal R(\sigma)P_d)=\mathrm{Tr}(N_d\sigma)$,
$\mathcal R^*(P_d)=N_d$. En utilisant $N_dE_iP_c=\delta_{cd}E_iP_c$ (étape précédente), puis $\sum_iE_i^*E_i=1$,
$$
P(\mathcal R\circ\Phi)^*(P_d)P=P\,\Phi^*(N_d)\,P=\sum_i\sum_{c,c'}P_{c'}E_i^*N_dE_iP_c=\sum_i\sum_{c'}P_{c'}E_i^*E_iP_d=P_d\Big(\sum_iE_i^*E_i\Big)P_d=P_d .
$$
(3)$\Rightarrow$(1). Soient $N_c=\mathcal R^*(P_c)$ pour $c<C$ et $N_C=\mathcal R^*(1-\sum_{c<C}P_c)$ ; ils
sont positifs et de somme $\mathcal R^*(1)=1$. Pour $\rho$ à support dans $P_c\mathbb C^D$ et $d<C$,
$\mathrm{Tr}(N_d\Phi(\rho))=\mathrm{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(P_d)P)=\mathrm{Tr}(\rho P_d)=\delta_{cd}$,
donc aussi $\mathrm{Tr}(N_C\Phi(\rho))=\delta_{cC}$. $\square$

<a id="r-10-7"></a>

**Théorème 10.7 (Bény–Kempf–Kribs).** Une algèbre de von Neumann
$\mathfrak A$ sur $\mathcal C$ est corrigeable pour $\Phi$ si et seulement si $PE_i^*E_jP\in\mathfrak A'$
pour tous $i,j$, où $\mathfrak A'$ est le commutant de $\mathfrak A$ dans $B(\mathcal C)$.

*Démonstration.* C'est le théorème D.9, démontré dans l'[annexe D](D_gksl.md), dont la démonstration construit aussi
un canal de récupération : il mesure le secteur de $\mathfrak A\cong\bigoplus_cM_{k_c}$ et un syndrome, et défait l'erreur
sur chaque facteur $\mathbb C^{k_c}$ comme dans la démonstration du théorème 10.5. $\square$

Pour $\mathfrak A=B(\mathcal C)$, $\mathfrak A'=\mathbb CP$ et l'on retrouve le théorème 10.5 ; la
proposition 10.6 est le cas de l'algèbre commutative $\mathfrak A=\mathrm{span}\{P_c\}$ (corollaire D.10(a),(b)). Pour
$\mathfrak A=M_k\otimes1_m$ sur $\mathcal C=\mathbb C^k\otimes\mathbb C^m$, on obtient les codes
d'opérateurs (codes de sous-systèmes) ; les sous-systèmes sans bruit de la section 10.1 sont de tels codes, corrigés par
$\mathcal R=\mathrm{id}$ (corollaire D.10(c), les deux facteurs tensoriels étant échangés). D'après le théorème 2.6, toute
$\mathfrak A$ est $*$-isomorphe à une algèbre $\bigoplus_cM_{k_c}$ ; lorsque $C\ge2$ et que $k_c\ge2$ pour au moins un $c$, une
$\mathfrak A$ corrigeable protège l'étiquette de secteur $c$ en même temps que l'information quantique de chaque $M_{k_c}$ :
c'est un code **hybride** (remarque suivant le corollaire D.10, avec un exemple où la condition de Knill–Laflamme est en
défaut mais où le théorème 10.7 s'applique). Historiquement : C. Bény, A. Kempf, D. W. Kribs, *Phys. Rev. Lett.* 98,
100502 (2007) et *Phys. Rev. A* 76, 042303 (2007) ; codes d'opérateurs : Kribs–Laflamme–Poulin,
*Phys. Rev. Lett.* 94, 180501 (2005).

**Exemple.** Le code de répétition soumis à des inversions de phase indépendantes (probabilité $q$
par qubit), suivies d'au plus une inversion de bit. Les opérateurs de Kraus sont des multiples des
$X_kZ^s$, $s\in\{0,1\}^3$, et chaque $E_i^*E_j$ inverse au plus deux bits, donc la condition 2 est
vérifiée pour $P_0=|000\rangle\langle000|$, $P_1=|111\rangle\langle111|$ : le bit classique survit, lu par
vote majoritaire, $N_0=\sum_{|b|\le1}|b\rangle\langle b|$, $N_1=1-N_0$. Le qubit, lui, ne survit pas :
$PZ_1P$ n'est pas un multiple de $P$, et après la récupération $\mathcal R$ de la section 10.2 l'état
$|{+_L}\rangle=(|000\rangle+|111\rangle)/\sqrt2$ subit l'inversion de phase logique $Z_1Z_2Z_3$ exactement
lorsqu'un nombre impair d'inversions de phase s'est produit ; sa fidélité vaut
$(1+(1-2q)^3)/2=0.8645$ pour $q=0.05$. (Le canal est $\Phi_{\mathrm{one}}\circ\mathcal Z$, où $\Phi_{\mathrm{one}}$ est le
canal « au plus une inversion de bit » de l'exemple des inversions de bit et $\mathcal Z$ les inversions de phase
indépendantes ; sur $\mathcal C$, chaque $Z^s$ agit comme $Z_1^{|s|}$, donc
$\mathcal Z(\rho)=\sum_sq^{|s|}(1-q)^{3-|s|}Z_1^{|s|}\rho Z_1^{|s|}$ est de nouveau à support dans $\mathcal C$, et
$\mathcal R\circ\Phi_{\mathrm{one}}$ est l'identité sur de tels états d'après le théorème 10.5.
$Z_1$ envoie $|{+_L}\rangle$ sur $|{-_L}\rangle\perp|{+_L}\rangle$, et $|s|$ est pair avec probabilité
$\frac12\big((1-q+q)^3+(1-q-q)^3\big)$.)

*Calcul :* la proposition 10.6 sur cet exemple : la condition 2, le vote majoritaire $N_0,N_1$ et la fidélité $0.8645$ de $|{+_L}\rangle$ — [tutoriel, section 10.3](../../tutorial/fr/10_quantum_computing.md#s-10-3).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Protéger une
> algèbre commutative, c'est protéger de l'information classique (un bit, l'étiquette d'un secteur de
> supersélection), ce qui est bien moins coûteux que de protéger un qubit : un bruit qui commute avec
> les projecteurs des secteurs est inoffensif (alors $P_cE_i^*E_jP_d=E_i^*E_jP_cP_d=0$ pour $c\neq d$,
> condition 2 de la proposition 10.6).

<a id="s-10-4"></a>

## 10.4 Mélangeurs sous contraintes

L'ansatz à opérateurs alternés quantique (Hadfield et al., *Algorithms* 12, 34, 2019) prépare
$|\beta,\gamma\rangle=e^{-i\beta_pH_M}e^{-i\gamma_pH_C}\cdots e^{-i\beta_1H_M}e^{-i\gamma_1H_C}|b_0\rangle$, où $H_C$ est
un hamiltonien de coût diagonal et $H_M$ un **mélangeur**. Pour une contrainte $|b|=w$ (bissection de
graphe, portefeuilles de $w$ actifs), on prend $|b_0|=w$ et le **mélangeur XY** sur un anneau,
$H_M=\tfrac14\sum_{(i,j)}(X_iX_j+Y_iY_j)$.

<a id="r-10-8"></a>

**Proposition 10.8.** $H_M$, $H_C$, $e^{-i\beta H_M}$ et $e^{-i\gamma H_C}$ appartiennent à
$\mathcal M=\{P_0,\dots,P_n\}'=\bigoplus_wB(H_w)\cong\bigoplus_wM_{\binom nw}$, donc $|\beta,\gamma\rangle\in H_w$ :
l'ansatz ne quitte jamais le sous-espace admissible.

*Démonstration.* Sur les qubits $i,j$, $X\otimes X+Y\otimes Y=2(|01\rangle\langle10|+|10\rangle\langle01|)$
(utiliser $Y|0\rangle=i|1\rangle$, $Y|1\rangle=-i|0\rangle$) : $X_iX_j+Y_iY_j$ échange les bits $i,j$ lorsqu'ils
diffèrent et annule $|b\rangle$ sinon, donc il envoie chaque $H_w$ dans lui-même. Un $H_C$ diagonal
commute avec les $P_w$ diagonaux, et les exponentielles aussi (théorème 1.5(d)). $\square$

*Dans la bibliothèque :* $\mathcal M$ est `SpinChain(n).algebra` et le mélangeur XY est `SpinChain.xxz` avec $\Delta=0$ — [tutoriel, section 10.4](../../tutorial/fr/10_quantum_computing.md#s-10-4). *Calcul :* la proposition 10.8 pour la bissection de graphe sur six qubits, avec la bibliothèque et avec des matrices denses — [tutoriel, section 10.4](../../tutorial/fr/10_quantum_computing.md#s-10-4).

L'**amortissement d'amplitude** du qubit $i$ ($|1\rangle\to|0\rangle$ avec probabilité $g$) a pour
opérateurs de Kraus $K_0^{(i)}=|0\rangle\langle0|_i+\sqrt{1-g}\,|1\rangle\langle1|_i$ et
$K_1^{(i)}=\sqrt g\,|0\rangle\langle1|_i$, qui abaisse le poids : c'est une application entre secteurs.

*Dans la bibliothèque :* `SpinChain.site_amplitude_damping`, un `InterSectorChannel` — [tutoriel, section 10.4](../../tutorial/fr/10_quantum_computing.md#s-10-4).

<a id="r-10-9"></a>

**Proposition 10.9 (fuite et post-sélection).** Soient $\rho_0$ à support dans $H_w$, $U_l$ des
unitaires de $\mathcal M$, $\mathcal D$ l'amortissement de tous les qubits, et
$\rho_l=\mathcal D(U_l\rho_{l-1}U_l^*)$. Alors chaque $\rho_l$ est diagonal par blocs, de poids $\le w$,
et
$$
P_w\rho_pP_w=(1-g)^{wp}\;U\rho_0U^*,\qquad U=U_p\cdots U_1 :
$$
le poids $w$ est obtenu avec probabilité $(1-g)^{wp}$, et l'état post-sélectionné est l'état sans
bruit.

*Démonstration.* Les opérateurs de Kraus de $\mathcal D$ sont les $K_s=\prod_iK^{(i)}_{s_i}$,
$s\in\{0,1\}^n$ ; $K_s$ envoie $H_v$ dans $H_{v-|s|}$, c'est-à-dire $P_vK_s=K_sP_{v+|s|}$ ($P_v=0$ pour
$v>n$). Si $\sigma$ est diagonal par blocs de poids $\le w$, alors
$P_vK_s\sigma K_s^*P_{v'}=K_sP_{v+|s|}\sigma P_{v'+|s|}K_s^*$ s'annule pour $v\neq v'$ et pour $v>w-|s|$ :
$\mathcal D(\sigma)$ est diagonal par blocs de poids $\le w$, et son bloc de poids $w$ provient du seul
$s=0$. Comme $K_0|b\rangle=(1-g)^{|b|/2}|b\rangle$, $P_w\mathcal D(\sigma)P_w=(1-g)^wP_w\sigma P_w$. Les $U_l$
commutent avec les $P_v$ ; on conclut par récurrence sur $l$. $\square$

*Calcul :* la proposition 10.9 sur le circuit de bissection, avec post-sélection par `DensityMatrix.condition_on` — [tutoriel, section 10.4](../../tutorial/fr/10_quantum_computing.md#s-10-4).

**Temps continu.** Il en va de même pour l'équation de Lindblad (théorème 8.6, [chapitre 8](08_dynamics.md))
$$
\dot\rho=-i[H,\rho]+\gamma_1\sum_i\Big(S^-_i\rho S^+_i-\tfrac12\{S^+_iS^-_i,\rho\}\Big)+\gamma_2\sum_k\Big(L_k\rho L_k^*-\tfrac12\{L_k^*L_k,\rho\}\Big)
$$
avec $H\in\mathcal M$ et des sauts $L_k\in\mathcal M$ à l'intérieur des secteurs, par exemple le
déphasage pur $L_k=S^z_k$. Ici $\sum_iS^+_iS^-_i=\sum_i|1\rangle\langle1|_i=N$, l'**opérateur nombre**,
$NP_v=vP_v$ : il est central.

<a id="r-10-10"></a>

**Proposition 10.10.** Si $\rho(0)$ est à support dans $H_w$, alors $P_w\rho(t)P_w=e^{-\gamma_1wt}\sigma(t)$,
où $\sigma$ est solution de la même équation avec $\gamma_1=0$ et $\sigma(0)=\rho(0)$. Le taux
d'acceptation vaut donc $e^{-\gamma_1wt}$, et l'état post-sélectionné $\sigma(t)$ ne dépend pas de
$\gamma_1$.

*Démonstration.* Soit $B_{vv'}=P_v\rho P_{v'}$. Comme $H,L_k,N$ commutent avec les $P_v$ et que
$P_vS^-_i=S^-_iP_{v+1}$,
$\dot B_{vv'}=-i[H,B_{vv'}]+\gamma_1\sum_iS^-_iB_{v+1,v'+1}S^+_i-\tfrac{\gamma_1}2(v+v')B_{vv'}+\gamma_2\mathcal L_2(B_{vv'})$,
où $\mathcal L_2$ est la partie de déphasage. Les blocs avec $\max(v,v')>w$ vérifient un système
linéaire fermé de valeurs initiales nulles, donc ils sont nuls (unicité dans le lemme 8.1). Alors
$\dot B_{ww}=-i[H,B_{ww}]-\gamma_1wB_{ww}+\gamma_2\mathcal L_2(B_{ww})$, de sorte que
$e^{\gamma_1wt}B_{ww}$ est solution de l'équation de $\sigma$ avec la même valeur initiale, et est égal
à $\sigma$ (lemme 8.1). Cette équation est encore de forme GKSL, donc $\mathrm{Tr}\,\sigma(t)=1$
(théorème 8.6). $\square$

Comparer avec la proposition 8.11, qui admet un $\sum_jL_j^*L_j$ non central mais aucun saut à
l'intérieur des secteurs : ici le terme de perte $\gamma_1N$ est central, si bien que le facteur sans
saut est le scalaire $e^{-\gamma_1wt}$ et que le déphasage à l'intérieur du secteur est permis.

Pour $\gamma_2=0$, $\sigma(t)=e^{-iHt}\rho(0)e^{iHt}$ : la post-sélection est parfaite.

*Calcul :* la proposition 10.10 pour le mélangeur XY sur six qubits, avec `dynamics.lindblad_evolve` — [tutoriel, section 10.4](../../tutorial/fr/10_quantum_computing.md#s-10-4).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). La relaxation d'énergie ($T_1$), l'une des principales erreurs des
> qubits supraconducteurs, est un saut vers un secteur inférieur. La mesure du poids de Hamming —
> un projecteur central, inoffensif à l'intérieur d'un secteur (section 3.3) — la détecte, et la
> post-sélection l'élimine au prix du taux d'acceptation. Le déphasage commute avec les projecteurs
> des secteurs et n'est pas détecté.

<a id="s-10-5"></a>

## 10.5 Points de comparaison

*Dans la bibliothèque :* la comparaison avec d'autres logiciels (Qiskit, Cirq, PennyLane, QuTiP, stim) — [tutoriel, section 10.5](../../tutorial/fr/10_quantum_computing.md#s-10-5).

<a id="x-10"></a>

## Exercices

1. (Le théorème 10.2 exige l'hypothèse unifère ; comparer avec la section 6.4.) Pour
   l'amortissement d'amplitude d'un qubit, $E_0=|0\rangle\langle0|+\sqrt{1-g}\,|1\rangle\langle1|$,
   $E_1=\sqrt g\,|0\rangle\langle1|$, $0<g<1$, montrer que
   $\mathcal A'=\mathbb C1$, $\mathrm{Fix}(\Phi)=\mathbb C|0\rangle\langle0|$ et $\mathrm{Fix}(\Phi^*)=\mathbb C1$.
2. Montrer que $\sum_w\binom nw^2=\binom{2n}n$. Pour le canal du corollaire 10.3 avec $\theta=\pi/2$
   (qui viole son hypothèse), que valent $\mathcal A$ et
   $\mathrm{Fix}(\Phi)$ ? Vérifier $\dim\mathrm{Fix}(\Phi)$ pour $n=4$.
3. (Linéarité de la correction d'erreurs.) Soit $PE_a^*E_bP=\lambda_{ab}P$ pour un nombre fini
   d'opérateurs $E_a$ (pas nécessairement les opérateurs de Kraus d'un canal), et construisons
   $\mathcal R$ à partir de $\Lambda$ comme dans la démonstration du théorème 10.5. Montrer que
   $\mathcal R$ corrige tout canal dont les opérateurs de Kraus appartiennent à
   $\mathrm{span}\{E_a\}$. En déduire que le code de bit-flip corrige l'erreur cohérente
   $e^{-i\phi X_1}$.

<details>
<summary>Solutions</summary>

1. $E_1$, $E_1^*$ et leurs produits donnent toutes les unités matricielles, donc $\mathcal A=M_2$ et
   $\mathcal A'=\mathbb C1$. Pour $T=(x_{kl})$, $\Phi(T)=\begin{pmatrix}x_{00}+gx_{11}&\sqrt{1-g}\,x_{01}\\\sqrt{1-g}\,x_{10}&(1-g)x_{11}\end{pmatrix}=T$
   si et seulement si $x_{11}=x_{01}=x_{10}=0$, et
   $\Phi^*(T)=\begin{pmatrix}x_{00}&\sqrt{1-g}\,x_{01}\\\sqrt{1-g}\,x_{10}&(1-g)x_{11}+gx_{00}\end{pmatrix}$,
   qui est égal à $T$ si et seulement si $T=x_{00}1$. ($\Phi$ n'est pas unifère, et $\Phi^*$ est
   unifère mais ne préserve pas la trace.)

2. Dénombrer les parties à $n$ éléments de $\{1,\dots,2n\}$ selon le nombre $w$ de leurs éléments
   situés dans la première moitié : $\binom{2n}n=\sum_w\binom nw\binom n{n-w}=\sum_w\binom nw^2$. Pour
   $\theta=\pi/2$, $U=(-i)^n\sum_w(-1)^wP_w$ ne voit que la parité de $w$, donc
   $\mathcal A=\mathrm{span}\{P_{\mathrm{even}},P_{\mathrm{odd}}\}$ avec $P_{\mathrm{even}}=\sum_{w\text{ even}}P_w$,
   et d'après le théorème 10.2, $\mathrm{Fix}(\Phi)$ est formé des opérateurs diagonaux par
   blocs pour la parité : dimension $2\cdot(2^{n-1})^2=128$ pour $n=4$.

3. Avec $F_k=\sum_au_{ak}E_a$, $\mathrm{span}\{E_a\}=\mathrm{span}\{F_k\}$, donc chaque opérateur de
   Kraus du canal s'écrit $G_m=\sum_kh_{mk}F_k$, et $G_mP=\sum_{d_k>0}h_{mk}\sqrt{d_k}V_k$. Comme
   $V_l^*V_k=\delta_{lk}P$ et $(1-Q)V_k=0$, $\mathcal R(G_m\rho G_m^*)=\sum_l|h_{ml}|^2d_l\,\rho$ pour
   $\rho=P\rho P$, donc $\mathcal R(\Phi(\rho))=c\rho$ ; en prenant les traces, $c=1$. Pour
   $E_a\in\{1,X_1,X_2,X_3\}$, $\Lambda=1$, $V_k=X_kP$, et $\mathcal R$ est la récupération de la
   section 10.2 ; $e^{-i\phi X_1}=\cos\phi\,1-i\sin\phi\,X_1$ appartient à l'espace engendré.

</details>

*Calcul :* les exercices 1–3, numériquement — [tutoriel, Exercices](../../tutorial/fr/10_quantum_computing.md#x-10).

**Suite :** [chapitre 11](11_optimisation.md) : optimisation sur les groupes unitaires d'une algèbre, et
inégalités de traces.
