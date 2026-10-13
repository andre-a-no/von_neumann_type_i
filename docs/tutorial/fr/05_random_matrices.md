[English](../05_random_matrices.md) · [Español](../es/05_random_matrices.md) · **Français** · [Deutsch](../de/05_random_matrices.md) · [中文](../zh/05_random_matrices.md) · [日本語](../ja/05_random_matrices.md) · [Русский](../ru/05_random_matrices.md)

# 5. Matrices aléatoires dans une algèbre

La théorie des matrices aléatoires cherche à quoi ressemble un opérateur « typique ». Dans une
algèbre $\bigoplus_cM_{k_c}$, la réponse doit respecter les secteurs : les unitaires aléatoires
naturels sont ceux de l'algèbre, et non ceux de la grande algèbre de matrices $M_D$, et les
hamiltoniens aléatoires naturels ont des blocs indépendants. Ce chapitre construit la mesure de
Haar sur le groupe unitaire de l'algèbre, calcule exactement ses moments d'ordre un et deux, puis
montre ce qui se passe mal dans les statistiques spectrales lorsqu'on ignore les secteurs.

**Vous apprendrez** pourquoi la décomposition QR d'une matrice gaussienne fournit un unitaire de
Haar (et pourquoi il faut fixer les phases), que la moyenne sur le groupe unitaire de l'algèbre est
l'espérance conditionnelle sur le centre, la formule de Weingarten d'ordre deux, la surmise de
Wigner et le rapport d'espacements, et pourquoi les statistiques de niveaux calculées sur un
spectre où les secteurs sont mélangés ressemblent à celles de niveaux indépendants.

Dans tout ce chapitre, $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ avec $D=\sum_ck_c$, les
$P_c$ sont les projecteurs centraux minimaux, et $E$ est l'espérance conditionnelle sur le centre
de la proposition 2.8. $U(k)$ désigne le groupe des matrices unitaires $k\times k$ et
$U(\mathcal M)=\{U\in\mathcal M:U^*U=1\}$ le groupe unitaire de $\mathcal M$. Les variables
aléatoires sont définies sur un espace probabilisé d'espérance $\mathbb E$ ; l'espérance d'une
matrice aléatoire se prend coefficient par coefficient. Une variable *gaussienne complexe
standard* a des parties réelle et imaginaire indépendantes de loi $\sim N(0,\tfrac12)$, de sorte
que $\mathbb E|z|^2=1$ (c'est ce que tire `torch.randn` pour un dtype complexe) ; une *matrice de
Ginibre* a des coefficients indépendants gaussiens complexes standard.

## 5.1 Mesure de Haar

Soit $G$ un groupe compact de matrices (un sous-groupe fermé d'un certain $U(n)$). Une mesure de
probabilité borélienne $\mu$ sur $G$ est **invariante à gauche** si $\mu(gB)=\mu(B)$ pour tout
$g\in G$ et tout borélien $B$, ou, de manière équivalente, si $\int f(gh)\,d\mu(h)=\int f\,d\mu$
pour toute fonction borélienne bornée $f$ ; l'**invariance à droite** se définit avec $Bg$ et
$f(hg)$.

**Lemme 5.1 (unicité).** Si $\mu$ est une mesure de probabilité borélienne invariante à gauche et
$\nu$ une mesure de probabilité borélienne invariante à droite sur $G$, alors $\mu=\nu$. Par
conséquent, $G$ porte au plus une mesure de probabilité borélienne invariante à gauche, et si elle
existe, elle est aussi invariante à droite et invariante par $g\mapsto g^{-1}$.

*Démonstration.* Soit $f:G\to\mathbb R$ continue (donc bornée). D'après le théorème de Fubini,
$\int\!\big(\int f(gh)\,d\mu(h)\big)d\nu(g)=\int\!\big(\int f(gh)\,d\nu(g)\big)d\mu(h)$.
L'intégrale intérieure de gauche vaut $\int f\,d\mu$ par invariance à gauche, celle de droite vaut
$\int f\,d\nu$ par invariance à droite ; donc $\int f\,d\mu=\int f\,d\nu$. Deux mesures de
probabilité boréliennes sur un espace métrique qui ont les mêmes intégrales sur les fonctions
continues sont égales : les fonctions continues $\max(0,1-n\operatorname{dist}(\cdot,K))$ décroissent
vers l'indicatrice d'un fermé $K$, donc les mesures coïncident sur les fermés (convergence dominée),
qui engendrent la tribu borélienne et sont stables par intersection (théorème $\pi$–$\lambda$ de
Dynkin). Si $\mu_1,\mu_2$ sont invariantes à gauche, alors $\nu(B)=\mu_2(B^{-1})$ est invariante à
droite, puisque $(Bg)^{-1}=g^{-1}B^{-1}$ ; donc $\mu_1=\nu=\mu_2$, et $\mu_1$ est invariante à
droite et par inversion. $\square$

**Définition.** La **mesure de Haar** de $G$ est sa mesure de probabilité borélienne invariante à
gauche. Un élément aléatoire de cette loi est dit **distribué selon la loi de Haar**.

Tout groupe compact possède une mesure de Haar (voir G. B. Folland, *A Course in Abstract Harmonic
Analysis*, 2e éd., CRC Press, 2016, section 2.2 ; énoncé ici sans démonstration). Pour $U(k)$,
nous n'avons pas besoin de ce théorème général : le théorème 5.3 construit la mesure.

**Lemme 5.2 (décomposition QR).** Toute $Z\in M_k(\mathbb C)$ inversible s'écrit $Z=QR$ avec $Q$
unitaire et $R$ triangulaire supérieure à diagonale positive. $Q$ et $R$ sont uniques et dépendent
continûment de $Z$.

*Démonstration.* Appliquons le procédé de Gram–Schmidt aux colonnes $z_1,\dots,z_k$ de $Z$, qui
sont linéairement indépendantes : $q_j=w_j/\|w_j\|$ avec
$w_j=z_j-\sum_{i<j}\langle q_i,z_j\rangle q_i$ ; ici $w_j\neq0$ car
$z_j\notin\operatorname{span}\{z_1,\dots,z_{j-1}\}=\operatorname{span}\{q_1,\dots,q_{j-1}\}$. Alors
$z_j=\sum_{i<j}\langle q_i,z_j\rangle q_i+\|w_j\|q_j$, c'est-à-dire $Z=QR$ avec $R_{jj}=\|w_j\|>0$,
et chaque étape est continue sur les matrices inversibles. Si $Q_1R_1=Q_2R_2$, alors
$T=Q_2^*Q_1=R_2R_1^{-1}$ est unitaire et triangulaire supérieure à diagonale positive. Son inverse
$T^*$ est triangulaire inférieure, et aussi triangulaire supérieure (l'inverse d'une matrice
triangulaire supérieure est triangulaire supérieure), donc $T$ est diagonale ; une matrice
unitaire diagonale à diagonale positive vaut $1$. $\square$

**Théorème 5.3 (unitaires de Haar à partir de matrices gaussiennes).** Soit $Z$ une matrice de
Ginibre $k\times k$. Alors $Z$ est inversible avec probabilité un, et le facteur $Q$ du lemme 5.2
est distribué selon la loi de Haar sur $U(k)$. De plus, si $Z=Q'R'$ est une décomposition QR
quelconque avec $Q'$ unitaire et $R'$ triangulaire supérieure (par exemple celle que renvoie
LAPACK, dont la diagonale est réelle de signe quelconque), alors $Q=Q'\Lambda$ avec
$\Lambda=\operatorname{diag}(R'_{ii}/|R'_{ii}|)$ est distribué selon la loi de Haar.

*Démonstration.* *Inversibilité.* $\det Z$ est un polynôme non nul en les $2k^2$ coordonnées
réelles de $Z$, et l'ensemble des zéros d'un polynôme réel non nul est de mesure de Lebesgue nulle
(récurrence sur le nombre de variables : écrivons $p=\sum_jp_j(x')x_n^j$ avec un $p_j\neq0$ ; pour
$x'$ hors de l'ensemble négligeable des zéros de $p_j$, il n'y a qu'un nombre fini de racines
$x_n$ ; on conclut par Fubini). $Z$ admet une densité, donc $P(\det Z=0)=0$.

*Invariance de la loi gaussienne.* La densité de $Z$ par rapport à la mesure de Lebesgue sur
$M_k(\mathbb C)\cong\mathbb R^{2k^2}$ est $\pi^{-k^2}e^{-\|Z\|_2^2}$. Pour $V\in U(k)$,
l'application $Z\mapsto VZ$ est linéaire sur les réels et conserve $\|Z\|_2$ (section 1.4) ; c'est donc
une application orthogonale de $\mathbb R^{2k^2}$, de jacobien $1$, qui conserve la densité : $VZ$
a même loi que $Z$.

*Équivariance.* Si $Z=QR$, alors $VZ=(VQ)R$ est la décomposition du lemme 5.2 de $VZ$, donc
$Q(VZ)=VQ(Z)$. Par conséquent, pour tout borélien $B\subseteq U(k)$,
$P(VQ(Z)\in B)=P(Q(VZ)\in B)=P(Q(Z)\in B)$ : la loi de $Q(Z)$ est invariante à gauche et, d'après le
lemme 5.1, c'est la mesure de Haar (ce qui prouve aussi que la mesure de Haar de $U(k)$ existe).
Enfin $Q'R'=(Q'\Lambda)(\Lambda^*R')$, et $\Lambda^*R'$ est triangulaire supérieure de diagonale
$|R'_{ii}|>0$, donc $Q'\Lambda=Q(Z)$ par unicité. $\square$

La même démonstration avec des matrices gaussiennes réelles donne la mesure de Haar du groupe
orthogonal $O(k)$. Sans la correction $\Lambda$, le résultat n'est *pas* distribué selon la loi de
Haar : la convention de LAPACK pour les signes des $R'_{ii}$ introduit un biais, que le code
ci-dessous met en évidence (F. Mezzadri, *How to generate random matrices from the classical
compact groups*, Notices AMS 54 (2007), 592–604).

**Corollaire 5.4 (mesure de Haar de l'algèbre).** $U(\mathcal M)=\{\bigoplus_cU_c:U_c\in U(k_c)\}\cong\prod_cU(k_c)$,
et $U=\bigoplus_cU_c$ est distribué selon la loi de Haar sur $U(\mathcal M)$ si et seulement si les
blocs $U_c$ sont indépendants et chaque $U_c$ est distribué selon la loi de Haar sur $U(k_c)$.

*Démonstration.* $U^*U=\bigoplus_cU_c^*U_c$, donc $U$ est unitaire si et seulement si chaque bloc
l'est. Supposons les $U_c$ indépendants et distribués selon la loi de Haar. Pour
$V=\bigoplus_cV_c\in U(\mathcal M)$, les blocs $V_cU_c$ de $VU$ sont indépendants (fonctions de
variables indépendantes) et $V_cU_c$ a la loi de $U_c$, donc $VU$ a la loi de $U$. La loi de $U$
est donc invariante à gauche et, d'après le lemme 5.1, c'est la mesure de Haar du groupe compact
$U(\mathcal M)$. Réciproquement, par unicité, tout $U$ distribué selon la loi de Haar a cette
loi. $\square$

Le groupe $U(\mathcal M)$ est bien plus petit que $U(D)$ : pour $C\ge2$, un unitaire de Haar de
$U(D)$ appartient à $\mathcal M$ avec probabilité nulle (d'après le théorème 5.3 et la construction
de Gram–Schmidt du lemme 5.2, sa première colonne est $z/\|z\|$ pour un vecteur gaussien $z$, de
sorte que tous ses coefficients sont presque sûrement non nuls). La mesure de Haar de $SU(k)$
s'obtient dans le même esprit : si $U$ suit la loi de Haar sur $U(k)$ et $\theta=\arg\det U$, alors
$V=e^{-i\theta/k}U\in SU(k)$, et pour $W\in SU(k)$, $\det(WU)=\det U$ donne $V(WU)=WV(U)$ ; comme
$WU$ a la loi de $U$, la loi de $V$ est invariante à gauche sur $SU(k)$, et c'est donc la mesure de
Haar de $SU(k)$ d'après le lemme 5.1.

**Corollaire 5.5 (premiers moments des coefficients).** Si $U$ est distribué selon la loi de Haar
sur $U(k)$, alors $\mathbb EU_{ij}=0$ et $\mathbb E|U_{ij}|^2=1/k$.

*Démonstration.* $-1\in U(k)$, donc $-U$ a la loi de $U$ et $\mathbb EU=-\mathbb EU$. Pour une
matrice de permutation $P$, $PU$ a la loi de $U$, donc $\mathbb E|U_{ij}|^2$ ne dépend pas de $i$ ;
comme $\sum_i|U_{ij}|^2=1$, cette quantité vaut $1/k$. $\square$

La bibliothèque tire des unitaires de Haar de $U(\mathcal M)$ avec `random_unitary_operator`, au
moyen d'une seule décomposition QR par lots par taille de bloc, pour tous les secteurs et tous les
échantillons à la fois (sur GPU, par Gram–Schmidt, qui donne le même $Q$ grâce à l'unicité du
lemme 5.2). La mesure `'haar'` correspond au théorème 5.3 dans chaque bloc ; pour une algèbre
complexe, `'haar_su'` multiplie par $e^{-i\theta/k}$ comme ci-dessus ; `'coe'` et `'cse'` sont les
ensembles circulaires de Dyson (le premier fait l'objet de l'exercice 3). Nous vérifions les
corollaires 5.4 et 5.5, et que l'échantillonneur de la bibliothèque est exactement la QR corrigée
du théorème 5.3 :

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, SpinChain, channels

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def mc(samples, exact, name):
    """Monte Carlo mean of real samples with its standard error; must agree within 5 standard errors."""
    mean, err = samples.mean().item(), samples.std().item() / math.sqrt(samples.numel())
    print(f"{name}: {mean:.4f} +- {err:.4f}   (theory {exact:.4f})")
    assert abs(mean - exact) <= 5 * err + 1e-12
    return mean

alg = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')        # M_1 (+) M_2 (+) M_3, D = 6
B = 100_000
U = alg.random_unitary_operator(B)                                   # Haar on U(1) x U(2) x U(3)
assert close(U.adjoint() @ U, alg.identity())                        # unitary elements of the algebra
u = [U.matrix[:, c, :k, :k] for c, k in enumerate(alg.k_factors)]   # the blocks U_c
for c, k in enumerate(alg.k_factors):                                # Corollary 5.5
    mc(u[c][:, 0, 0].real, 0.0, f"E Re (U_{c})_11")
    mc(u[c][:, 0, 0].abs() ** 2, 1 / k, f"E |(U_{c})_11|^2")
mc(u[1][:, 0, 0].abs() ** 2 * u[2][:, 0, 0].abs() ** 2, 1 / 6, "E |(U_1)_11|^2 |(U_2)_11|^2")   # independent blocks

SU = alg.random_unitary_operator(4, measure='haar_su').matrix
assert torch.allclose(torch.linalg.det(SU[:, 2, :3, :3]), torch.ones(4, dtype=dt))

torch.manual_seed(1)
Z = torch.randn(B, 2, 2, dtype=dt)                                   # Ginibre matrices
torch.manual_seed(1)
V = alg.random_unitary(2, batch_size=B)                              # the library, same random numbers
Q, R = torch.linalg.qr(Z)                                            # LAPACK: diag(R) real, either sign
d = R.diagonal(dim1=-2, dim2=-1)
assert close(V, Q * (d / d.abs()).unsqueeze(-2))                     # Theorem 5.3: V = Q Lambda
raw = Q.diagonal(dim1=-2, dim2=-1).sum(-1).real                      # Re Tr Q without the correction
print(f"E Re Tr Q without the phase correction: {raw.mean():.4f} +- {raw.std() / math.sqrt(B):.4f}")
assert abs(raw.mean()) > 20 * raw.std() / math.sqrt(B)               # far from 0: not Haar
```

Sans la correction, $\mathbb E\operatorname{Re}\operatorname{Tr}Q\approx-0.84$ pour $k=2$, alors qu'un
unitaire de Haar vérifie $\mathbb EU=0$.

## 5.2 Moments et calcul de Weingarten

**Lemme 5.6.** Si $T\in M_k(\mathbb C)$ commute avec tous les unitaires diagonaux et toutes les
matrices de permutation, alors $T\in\mathbb C1$. En particulier, le commutant de $U(k)$ est
$\mathbb C1$.

*Démonstration.* Pour $\Delta=\operatorname{diag}(e^{i\theta_1},\dots,e^{i\theta_k})$ avec des
$e^{i\theta_j}$ deux à deux distincts, $\Delta T=T\Delta$ donne
$(e^{i\theta_i}-e^{i\theta_j})T_{ij}=0$, donc $T$ est diagonale. Pour la matrice de permutation
$P_\sigma e_j=e_{\sigma(j)}$, $P_\sigma TP_\sigma^{-1}=T$ donne $T_{\sigma(j)\sigma(j)}=T_{jj}$, donc
tous les coefficients diagonaux sont égaux. $\square$

**Proposition 5.7 (le twirl de Haar est l'espérance conditionnelle sur le centre).** Soit $U$
distribué selon la loi de Haar sur $U(\mathcal M)$. Pour tout $X\in\mathcal M$,
$$
\mathbb E\,[UXU^*]=E(X)=\bigoplus_c\frac{\operatorname{Tr}X_c}{k_c}\,1_{k_c}.
$$

*Démonstration.* $(UXU^*)_c=U_cX_cU_c^*$, et $U_c$ suit la loi de Haar sur $U(k_c)$
(corollaire 5.4). Posons $T_c=\mathbb E[U_cX_cU_c^*]$ (les coefficients sont bornés). Pour
$V\in U(k_c)$, $VU_c$ a la loi de $U_c$, donc $VT_cV^*=\mathbb E[(VU_c)X_c(VU_c)^*]=T_c$. D'après le
lemme 5.6, $T_c=\lambda1$, et
$k_c\lambda=\operatorname{Tr}T_c=\mathbb E\operatorname{Tr}(U_cX_cU_c^*)=\operatorname{Tr}X_c$. $\square$

La normalisation est celle de la proposition 2.8 et de `channels.center_expectation`. Si
$X\in\mathcal M$ commute avec $U(\mathcal M)$, alors $UXU^*=X$, donc $X=E(X)$ est central ; la
moyenne de Haar est l'analogue continu de la moyenne sur le groupe du chapitre 2, exercice 3. Si au
contraire $U$ suit la loi de Haar sur $U(D)$, la même démonstration (avec $C=1$) donne
$\mathbb E[UXU^*]=\frac{\operatorname{Tr}X}{D}1$ : les poids des secteurs sont perdus, et $UXU^*$
n'est même pas dans $\mathcal M$.

> **Interprétation physique.** Un unitaire aléatoire de Haar de l'algèbre est une évolution
> aléatoire qui respecte la règle de supersélection. En moyennant sur lui, on brouille tout à
> l'intérieur de chaque secteur et on conserve exactement la partie classique : le poids de chaque
> charge. C'est le « twirl » de l'information quantique, et c'est l'espérance conditionnelle sur le
> centre du chapitre 2.

Pour les moments d'ordre deux, il nous faut le commutant de $\{V\otimes V\}$ sur
$\mathbb C^k\otimes\mathbb C^k$. Soit $F$ l'**échange** (*swap*), $F(x\otimes y)=y\otimes x$ ; il
est autoadjoint et unitaire, et
$\operatorname{Tr}1=k^2$, $\operatorname{Tr}F=\sum_{i,j}\langle e_i\otimes e_j,e_j\otimes e_i\rangle=k$, $F^2=1$.

**Lemme 5.8.** Soit $k\ge2$. Une matrice $X$ sur $\mathbb C^k\otimes\mathbb C^k$ commute avec
$V\otimes V$ pour tout $V\in U(k)$ si et seulement si $X=a1+bF$ avec $a,b\in\mathbb C$ ; de plus,
$1$ et $F$ sont linéairement indépendants.

*Démonstration.* $F(V\otimes V)(x\otimes y)=Vy\otimes Vx=(V\otimes V)F(x\otimes y)$, donc $1$ et $F$
commutent avec tous les $V\otimes V$. Réciproquement, supposons que $X$ commute avec eux, et
écrivons $X_{ij,pq}=\langle e_i\otimes e_j,X(e_p\otimes e_q)\rangle$.
*Unitaires diagonaux :* $\Delta\otimes\Delta$ multiplie $e_i\otimes e_j$ par
$e^{i(\theta_i+\theta_j)}$, donc $(e^{i(\theta_i+\theta_j)}-e^{i(\theta_p+\theta_q)})X_{ij,pq}=0$. Si
les multiensembles $\{i,j\}$ et $\{p,q\}$ diffèrent, un indice $m$ apparaît $n\in\{\pm1,\pm2\}$
fois de plus dans le premier ; avec $\theta_m=\pi/2$ et tous les autres $\theta=0$, le crochet
vaut $e^{i(\theta_p+\theta_q)}(e^{in\pi/2}-1)\neq0$, donc $X_{ij,pq}=0$. Par conséquent
$X(e_i\otimes e_i)=a_ie_i\otimes e_i$ et, pour $i\neq j$,
$X(e_i\otimes e_j)=b_{ij}\,e_i\otimes e_j+c_{ij}\,e_j\otimes e_i$.
*Permutations :* en appliquant $X(P_\sigma\otimes P_\sigma)=(P_\sigma\otimes P_\sigma)X$ à
$e_i\otimes e_j$, on obtient $a_{\sigma(i)}=a_i$, $b_{\sigma(i)\sigma(j)}=b_{ij}$,
$c_{\sigma(i)\sigma(j)}=c_{ij}$. Les permutations agissent transitivement sur les couples d'indices
distincts, donc $a_i=a'$, $b_{ij}=b'$, $c_{ij}=c'$ ; ainsi $X=b'1+c'F$ sur le sous-espace engendré
par les $e_i\otimes e_j$ avec $i\ne j$.
*Une rotation :* soit $R$ l'unitaire tel que $Re_1=(e_1+e_2)/\sqrt2$, $Re_2=(-e_1+e_2)/\sqrt2$, $Re_m=e_m$ pour $m\ge3$.
Alors $(R\otimes R)(e_1\otimes e_1)=\frac12(e_1\otimes e_1+e_1\otimes e_2+e_2\otimes e_1+e_2\otimes e_2)$, et
$X(R\otimes R)(e_1\otimes e_1)=\frac12\big(a'e_1\otimes e_1+(b'+c')(e_1\otimes e_2+e_2\otimes e_1)+a'e_2\otimes e_2\big)$
doit être égal à $(R\otimes R)X(e_1\otimes e_1)=a'(R\otimes R)(e_1\otimes e_1)$. En comparant les
coefficients de $e_1\otimes e_2$, on obtient $a'=b'+c'$, qui est la valeur de $b'1+c'F$ sur
$e_i\otimes e_i$. Donc $X=b'1+c'F$. Enfin, $1$ et $F$ envoient $e_1\otimes e_2$ sur les vecteurs
indépendants $e_1\otimes e_2$ et $e_2\otimes e_1$. $\square$

**Corollaire 5.9 (dualité de Schur–Weyl d'ordre deux).** Soient $k\ge2$ et $P_\pm=\frac12(1\pm F)$.
Ce sont les projecteurs sur les sous-espaces symétrique et antisymétrique
$\mathrm{Sym}=\{x:Fx=x\}$ et $\mathrm{Alt}=\{x:Fx=-x\}$, de dimensions $k(k+1)/2$ et $k(k-1)/2$.
La représentation $V\mapsto V\otimes V$ de $U(k)$ est irréductible sur $\mathrm{Sym}$ et sur
$\mathrm{Alt}$, et les deux représentations ne sont pas équivalentes.

*Démonstration.* $F=F^*$ et $F^2=1$ font de $P_\pm$ des projecteurs avec $P_++P_-=1$ et d'images
$\mathrm{Sym}$, $\mathrm{Alt}$ ; leurs dimensions sont $\operatorname{Tr}P_\pm=(k^2\pm k)/2$. Le
sous-espace $\mathcal A$ engendré par $\{V\otimes V\}$ est une $*$-algèbre ($1\otimes1$, produits
$VW\otimes VW$, adjoints $V^*\otimes V^*$) de commutant
$\mathcal A'=\operatorname{span}\{1,F\}=\operatorname{span}\{P_+,P_-\}$ d'après le lemme 5.8. Soit
$K\subseteq\mathrm{Sym}$ invariant par tous les $V\otimes V$. D'après le lemme 2.1, son projecteur
est dans $\mathcal A'$, il s'écrit donc $xP_++yP_-$ avec $x,y\in\{0,1\}$ (c'est un projecteur) ; il
s'annule sur $\mathrm{Alt}$, donc $y=0$, et $K=0$ ou $K=\mathrm{Sym}$. Le même argument s'applique
à $\mathrm{Alt}$. Si $T:\mathrm{Sym}\to\mathrm{Alt}$ entrelace les deux représentations, alors
$TP_+$ (vu comme opérateur sur $\mathbb C^k\otimes\mathbb C^k$) commute avec $\mathcal A$, donc
$TP_+=xP_++yP_-$ ; il s'annule sur $\mathrm{Alt}$, donc $y=0$, et son image est contenue dans
$\mathrm{Alt}\cap\mathrm{Sym}=0$, donc $T=0$. $\square$

**Théorème 5.10 (formule de Weingarten d'ordre deux).** Soient $k\ge2$, $U$ distribué selon la loi
de Haar sur $U(k)$, et $X$ une matrice sur $\mathbb C^k\otimes\mathbb C^k$. Alors
$$
\mathbb E\,[(U\otimes U)X(U\otimes U)^*]=\alpha\,1+\beta F,\qquad
\alpha=\frac{\operatorname{Tr}X-\operatorname{Tr}(FX)/k}{k^2-1},\quad
\beta=\frac{\operatorname{Tr}(FX)-\operatorname{Tr}X/k}{k^2-1}.
$$
De manière équivalente, avec $\mathrm{Wg}(e)=\frac1{k^2-1}$, $\mathrm{Wg}(\tau)=-\frac1{k(k^2-1)}$
sur le groupe $S_2=\{e,\tau\}$,
$$
\mathbb E\big[U_{i_1j_1}U_{i_2j_2}\overline{U_{i'_1j'_1}U_{i'_2j'_2}}\big]
=\sum_{\sigma,\pi\in S_2}\delta_{i_1i'_{\sigma(1)}}\delta_{i_2i'_{\sigma(2)}}\,
\delta_{j_1j'_{\pi(1)}}\delta_{j_2j'_{\pi(2)}}\,\mathrm{Wg}(\sigma\pi^{-1}).
$$

*Démonstration.* Posons $T=\mathbb E[(U\otimes U)X(U\otimes U)^*]$. Comme dans la proposition 5.7,
$T$ commute avec tous les $V\otimes V$, donc $T=\alpha1+\beta F$ d'après le lemme 5.8. En prenant
les traces, $\operatorname{Tr}T=\operatorname{Tr}X$ et, puisque $F$ commute avec $U\otimes U$,
$\operatorname{Tr}(FT)=\mathbb E\operatorname{Tr}((U\otimes U)^*F(U\otimes U)X)=\operatorname{Tr}(FX)$.
On obtient $k^2\alpha+k\beta=\operatorname{Tr}X$ et $k\alpha+k^2\beta=\operatorname{Tr}(FX)$, un
système linéaire de déterminant $k^2(k^2-1)\neq0$, dont la solution est celle de l'énoncé. Pour
les coefficients, prenons $X=|e_{j_1}\otimes e_{j_2}\rangle\langle e_{j'_1}\otimes e_{j'_2}|$ :
alors $\langle e_{i_1}\otimes e_{i_2},T\,e_{i'_1}\otimes e_{i'_2}\rangle$ est le membre de gauche,
$\operatorname{Tr}X=\delta_{j_1j'_1}\delta_{j_2j'_2}$ (le terme $\pi=e$),
$\operatorname{Tr}(FX)=\delta_{j_1j'_2}\delta_{j_2j'_1}$ (le terme $\pi=\tau$), et $1$, $F$
contribuent $\delta_{i_1i'_1}\delta_{i_2i'_2}$ (le terme $\sigma=e$) et
$\delta_{i_1i'_2}\delta_{i_2i'_1}$ ($\sigma=\tau$). En regroupant les termes,
$\alpha=\sum_\pi\mathrm{Wg}(e\pi^{-1})[\pi]$ et $\beta=\sum_\pi\mathrm{Wg}(\tau\pi^{-1})[\pi]$,
où $[\pi]$ désigne les deltas en $j$ associés à $\pi$. $\square$

La matrice $\big(\mathrm{Wg}(\sigma^{-1}\pi)\big)_{\sigma,\pi}$ est l'inverse de la matrice de Gram
$\begin{pmatrix}k^2&k\\k&k^2\end{pmatrix}$ de $1$ et $F$ pour
$\langle A,B\rangle=\operatorname{Tr}(A^*B)$ ; c'est sous cette forme que le résultat se
généralise.

**Corollaire 5.11.** Pour $k\ge2$ :
$$
\mathbb E|U_{11}|^4=\frac2{k(k+1)},\qquad
\mathbb E|U_{11}|^2|U_{12}|^2=\frac1{k(k+1)},\qquad
\mathbb E\big[U_{11}U_{22}\overline{U_{12}U_{21}}\big]=-\frac1{(k-1)k(k+1)}.
$$

*Démonstration.* Dans le théorème 5.10, (i) tous les indices valent $1$ : chaque produit de deltas
vaut $1$, donc le moment vaut $2(\mathrm{Wg}(e)+\mathrm{Wg}(\tau))=2\frac{k-1}{k(k^2-1)}$.
(ii) $i=(1,1)$, $j=(1,2)$, $i'=(1,1)$, $j'=(1,2)$ : seul $\pi=e$ subsiste, ce qui donne
$\mathrm{Wg}(e)+\mathrm{Wg}(\tau)$. (iii) $i=(1,2)$, $j=(1,2)$, $i'=(1,2)$, $j'=(2,1)$ : seuls
$\sigma=e$, $\pi=\tau$ subsistent, ce qui donne $\mathrm{Wg}(\tau)$. $\square$

Le dernier moment est négatif : les coefficients d'un unitaire de Haar ne sont pas indépendants,
et les corrélations sont d'ordre relatif $1/k$. Pour l'algèbre, les coefficients de blocs
différents sont indépendants (corollaire 5.4), de sorte que les moments mixtes entre secteurs se
factorisent en moments de blocs individuels.

Pour les moments d'ordre quelconque, le même raisonnement s'applique dès que l'on connaît le
commutant de $\{V^{\otimes n}\}$ ; il est engendré par les permutations des $n$ facteurs tensoriels
(dualité de Schur–Weyl, voir R. Goodman, N. R. Wallach, *Symmetry, Representations, and
Invariants*, Springer, 2009, chapitre 9). Nous énonçons le résultat sans démonstration.

**Théorème 5.12 (formule de Weingarten ; énoncé sans démonstration).** Soient $n\ge1$ et $U$
distribué selon la loi de Haar sur $U(k)$. Alors
$$
\mathbb E\big[U_{i_1j_1}\cdots U_{i_nj_n}\overline{U_{i'_1j'_1}\cdots U_{i'_nj'_n}}\big]
=\sum_{\sigma,\pi\in S_n}\prod_{m=1}^n\delta_{i_mi'_{\sigma(m)}}\delta_{j_mj'_{\pi(m)}}\;\mathrm{Wg}(\sigma\pi^{-1},k),
$$
où $\mathrm{Wg}(\cdot,k)$ est une fonction sur $S_n$ qui ne dépend que du type cyclique. Pour
$n\le k$, la matrice $\big(\mathrm{Wg}(\sigma^{-1}\pi,k)\big)_{\sigma,\pi}$ est l'inverse de la
matrice de Gram $\big(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)}\big)_{\sigma,\pi}$ ; pour $n>k$, cette
matrice est singulière et l'on prend son pseudo-inverse. Les moments où les nombres de facteurs
$U$ et $\overline U$ diffèrent sont nuls. (D. Weingarten, J. Math. Phys. 19 (1978), 999–1001 ;
B. Collins, Int. Math. Res. Not. 2003, 953–982 ; B. Collins, P. Śniady, Comm. Math. Phys. 264
(2006), 773–795.)

Nous vérifions la proposition 5.7 par rapport à `channels.center_expectation`, le commutant du
lemme 5.8 par algèbre linéaire exacte, et le théorème 5.10 et le corollaire 5.11 par Monte-Carlo
avec le bloc $U(3)$ des échantillons précédents :

```python
X = alg.from_blocks([torch.randn(k, k, dtype=dt) for k in alg.k_factors])      # a fixed element
T = (U @ X @ U.adjoint()).matrix                                    # B samples of U X U*
mean, err = T.mean(0), T.std(0) / math.sqrt(B)
EX = channels.center_expectation(alg)(X).matrix[0]
random = err > 1e-12                                                 # entries that are not constant
z = ((mean - EX).abs()[random] / err[random]).max().item()
print(f"twirl over U(M) vs E(X): largest deviation {z:.2f} standard errors")
assert ((mean - EX).abs() <= 5 * err + 1e-12).all()                 # Proposition 5.7

def commutant_dim(ops, tol=1e-9):
    """Dimension of {X : XA = AX for all A in ops}, as the null space of X -> (AX - XA)_A."""
    d = ops[0].shape[0]
    eye = torch.eye(d, dtype=dt)
    ops = [A.contiguous() for A in ops]                              # torch.kron needs contiguous input
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    return int((torch.linalg.svdvals(L) < tol).sum())

k = 3
F = torch.zeros(k * k, k * k, dtype=dt)
for i in range(k):
    for j in range(k):
        F[j * k + i, i * k + j] = 1                                  # F (e_i (x) e_j) = e_j (x) e_i
VV = [torch.kron(W, W) for W in alg.random_unitary(k, batch_size=3).contiguous()]
assert all(close(F @ W, W @ F) for W in VV)
assert commutant_dim(VV) == 2                                        # Lemma 5.8: span{1, F}

u3 = u[2]                                                            # Haar on U(3)
mc(u3[:, 0, 0].abs() ** 4, 2 / (k * (k + 1)), "E|U_11|^4")
mc(u3[:, 0, 0].abs() ** 2 * u3[:, 0, 1].abs() ** 2, 1 / (k * (k + 1)), "E|U_11|^2 |U_12|^2")
mc((u3[:, 0, 0] * u3[:, 1, 1] * (u3[:, 0, 1] * u3[:, 1, 0]).conj()).real,
   -1 / ((k - 1) * k * (k + 1)), "E U_11 U_22 conj(U_12 U_21)")

n_s = 20_000                                                         # Theorem 5.10 for a random X
X2 = torch.randn(k * k, k * k, dtype=dt)
trX, trFX = X2.diagonal().sum(), (F @ X2).diagonal().sum()
alpha, beta = (trX - trFX / k) / (k ** 2 - 1), (trFX - trX / k) / (k ** 2 - 1)
UU = torch.einsum('bij,bkl->bikjl', u3[:n_s], u3[:n_s]).reshape(n_s, k * k, k * k)   # U (x) U
T2 = UU @ X2 @ UU.mH
mean, err = T2.mean(0), T2.std(0) / math.sqrt(n_s)
assert ((mean - alpha * torch.eye(k * k, dtype=dt) - beta * F).abs() <= 5 * err).all()
```

## 5.3 Statistiques de niveaux

Soit $H$ une matrice autoadjointe aléatoire $n\times n$ de valeurs propres
$\lambda_1\le\dots\le\lambda_n$. Ses **espacements** sont les $s_i=\lambda_{i+1}-\lambda_i$. Deux
ensembles sont classiques : le **GOE**, $H=(G+G^T)/2$ avec $G$ réelle à coefficients indépendants
de loi $N(0,1)$, et le **GUE**, $H=(G+G^*)/\sqrt2$ avec $G$ de Ginibre (section 1.7). Dans les deux
cas, les coefficients diagonaux suivent la loi $N(0,1)$, les coefficients au-dessus de la
diagonale suivent la loi $N(0,\tfrac12)$ dans le GOE et ont des parties réelle et imaginaire
indépendantes de loi $N(0,\tfrac12)$ dans le GUE, et la densité de $H$ est
proportionnelle à $e^{-\operatorname{Tr}H^2/2}$.

**Proposition 5.13 (surmise de Wigner).** Soient $n=2$, $s=\lambda_2-\lambda_1$ et
$S=s/\mathbb Es$. Pour le GOE ($\beta=1$) et le GUE ($\beta=2$), $S$ a pour densités
$$
p_1(S)=\frac\pi2\,S\,e^{-\pi S^2/4},\qquad p_2(S)=\frac{32}{\pi^2}\,S^2\,e^{-4S^2/\pi},
$$
de fonctions de répartition $P(S\le t)=1-e^{-\pi t^2/4}$ et
$P(S\le t)=\operatorname{erf}(2t/\sqrt\pi)-\frac{4t}\pi e^{-4t^2/\pi}$. En particulier,
$P(S\le\varepsilon)$ est d'ordre $\varepsilon^{\beta+1}$ quand $\varepsilon\to0$.

*Démonstration.* Pour $H=\begin{pmatrix}a&b\\\bar b&c\end{pmatrix}$, les valeurs propres sont
$\frac{a+c}2\pm\frac12\sqrt{(a-c)^2+4|b|^2}$, donc $s=\sqrt{(a-c)^2+4|b|^2}$. *GOE :* $x=a-c$ et
$y=2b$ sont des variables $N(0,2)$ indépendantes, et en coordonnées polaires
$P(s\le t)=\int_{x^2+y^2\le t^2}\frac{e^{-(x^2+y^2)/4}}{4\pi}\,dx\,dy=\int_0^t\frac r2e^{-r^2/4}dr=1-e^{-t^2/4}$.
Alors $\mathbb Es=\int_0^\infty P(s>t)\,dt=\int_0^\infty e^{-t^2/4}dt=\sqrt\pi$, et
$P(S\le t)=P(s\le\sqrt\pi t)=1-e^{-\pi t^2/4}$ ; il reste à dériver. *GUE :* $x=a-c$,
$y=2\operatorname{Re}b$, $y'=2\operatorname{Im}b$ sont des $N(0,2)$ indépendantes, et en
coordonnées sphériques la densité de $s=\sqrt{x^2+y^2+y'^2}$ est
$4\pi r^2(4\pi)^{-3/2}e^{-r^2/4}=\frac{r^2e^{-r^2/4}}{2\sqrt\pi}$. Sa moyenne vaut
$\frac1{2\sqrt\pi}\int_0^\infty r^3e^{-r^2/4}dr=\frac{8}{2\sqrt\pi}=\frac4{\sqrt\pi}$. Le changement
de variable $r=4S/\sqrt\pi$ donne $p_2$, et l'intégration par parties de $p_2$
($\int x^2e^{-ax^2}dx=-\frac{x}{2a}e^{-ax^2}+\frac1{2a}\int e^{-ax^2}dx$ avec $a=4/\pi$) donne la
fonction de répartition. Le comportement pour $\varepsilon$ petit découle de
$p_\beta(S)\sim cS^\beta$. $\square$

L'exposant compte des dimensions : une matrice symétrique réelle $2\times2$ est dégénérée si et
seulement si deux conditions réelles sont satisfaites ($a=c$, $b=0$), une matrice hermitienne
complexe si et seulement si trois le sont. Pour tout $n$, la densité jointe des valeurs propres est
proportionnelle à $\prod_{i<j}|\lambda_i-\lambda_j|^\beta e^{-\sum_i\lambda_i^2/2}$, et quand
$n\to\infty$ la loi des espacements dans le cœur du spectre converge vers une limite
(Gaudin–Mehta) proche de la surmise, mais différente (M. L. Mehta, *Random Matrices*, 3e éd.,
Elsevier, 2004 ; énoncé sans démonstration).

> **Interprétation physique.** Wigner a proposé que les niveaux d'énergie des systèmes quantiques
> complexes (noyaux lourds) se comportent statistiquement comme les valeurs propres d'un
> hamiltonien aléatoire ayant la bonne symétrie : symétrique réel s'il existe une symétrie par
> renversement du temps ($\beta=1$), hermitien complexe sinon ($\beta=2$). Les niveaux se
> *repoussent* alors : deux niveaux de même symétrie ne se rapprochent presque jamais (von Neumann
> et Wigner, 1929). La répulsion est la signature typique du chaos quantique, des niveaux
> indépendants (sans répulsion) celle des systèmes intégrables ou localisés (conjectures de
> Bohigas–Giannoni–Schmit et de Berry–Tabor).

Les espacements dépendent de la densité locale de niveaux, qui varie le long du spectre ; pour les
comparer à $p_\beta$, il faut procéder à un **dépliement** (*unfolding*), le changement de variable
$\lambda\mapsto\bar N(\lambda)$ par une approximation lisse $\bar N$ de la fonction de comptage, qui
rend l'espacement moyen égal à $1$. Le choix de $\bar N$ est une source de biais. Le rapport
d'espacements consécutifs permet de l'éviter.

**Définition.** Le **rapport d'espacements** est
$r_i=\min(s_i,s_{i+1})/\max(s_i,s_{i+1})\in[0,1]$ (V. Oganesyan, D. A. Huse, Phys. Rev. B 75
(2007), 155111). Il est invariant par $\lambda\mapsto a\lambda+b$, $a>0$, et ne nécessite donc
aucun dépliement là où la densité est à peu près constante sur deux espacements.

**Proposition 5.14 (niveaux indépendants).** Soient $s,s'$ indépendantes, de même loi
exponentielle (comme deux espacements consécutifs d'un processus de Poisson). Alors
$r=\min(s,s')/\max(s,s')$ vérifie $P(r\le t)=\frac{2t}{1+t}$ sur $[0,1]$, a pour densité
$\frac2{(1+t)^2}$, et
$$
\mathbb Er=2\ln2-1\approx0.3863 .
$$

*Démonstration.* Par changement d'échelle, on peut supposer le paramètre égal à $1$. Pour
$x=s/s'$, $P(x\le t)=\mathbb E\,P(s\le ts'\mid s')=\mathbb E(1-e^{-ts'})=1-\frac1{1+t}=\frac t{1+t}$.
Pour $t<1$, les événements $\{x\le t\}$ et $\{1/x\le t\}$ sont disjoints, $x$ et $1/x$ ont même loi
($s,s'$ sont échangeables), et $\{r\le t\}$ est leur réunion, donc $P(r\le t)=\frac{2t}{1+t}$.
Alors
$\mathbb Er=\int_0^1\frac{2t}{(1+t)^2}dt=2\Big[\ln(1+t)+\frac1{1+t}\Big]_0^1=2\ln2-1$. $\square$

Il en va exactement de même pour $n$ niveaux uniformes indépendants sur un intervalle, car leurs
espacements ont la loi de $(e_1,\dots,e_{n+1})/\sum_je_j$ pour des exponentielles indépendantes
$e_j$, de sorte que chaque rapport $s_i/s_{i+1}$ a la loi de $e_i/e_{i+1}$ (W. Feller, *An Introduction to
Probability Theory and Its Applications*, vol. II, Wiley, 1971, chapitre I ; énoncé sans
démonstration). Pour le GOE et le GUE dans la limite $n\to\infty$, $\langle r\rangle=0.5307$ et
$0.5996$ (valeurs numériques ; les analogues $3\times3$ de la surmise donnent
$4-2\sqrt3\approx0.5359$ et $2\sqrt3/\pi-\frac12\approx0.6027$ ; Y. Y. Atas, E. Bogomolny,
O. Giraud, G. Roux, Phys. Rev. Lett. 110 (2013), 084101 ; énoncé sans démonstration).

La fonction `SpinChain.level_spacing_ratio` (dans `chains`) renvoie le rapport moyen de chaque
ligne de niveaux. Nous vérifions les propositions 5.13 et 5.14 ainsi que la valeur du GOE :

```python
ratio = SpinChain.level_spacing_ratio
n2 = 200_000
alg_r = TypeIAlgebra([2], [2], complex_valued=False, precision='double')     # real symmetric 2 x 2
alg_c = TypeIAlgebra([2], [2], precision='double')
G = torch.randn(n2, 1, 2, 2, dtype=torch.float64)
Gc = torch.randn(n2, 1, 2, 2, dtype=dt)
goe2 = alg_r.operator((G + G.mT) / 2, is_self_adjoint=True)
gue2 = alg_c.operator((Gc + Gc.mH) / math.sqrt(2), is_self_adjoint=True)
cdf = {1: lambda t: 1 - math.exp(-math.pi * t ** 2 / 4),
       2: lambda t: math.erf(2 * t / math.sqrt(math.pi)) - 4 * t / math.pi * math.exp(-4 * t ** 2 / math.pi)}
for beta, H, mean_s in ((1, goe2, math.sqrt(math.pi)), (2, gue2, 4 / math.sqrt(math.pi))):
    lam = H.eigenvalues()[0]                                         # (n2, 2), ascending
    S = (lam[:, 1] - lam[:, 0]) / mean_s
    mc(S, 1.0, f"beta={beta}: E S")
    for t in (0.25, 1.0, 2.0):
        mc((S <= t).double(), cdf[beta](t), f"beta={beta}: P(S <= {t})")

levels = torch.rand(4000, 50, dtype=torch.float64)                   # independent levels
mc(ratio(levels), 2 * math.log(2) - 1, "<r>, independent uniform levels")

n = 64
alg64 = TypeIAlgebra([n], [n], complex_valued=False, precision='double')
G = torch.randn(800, 1, n, n, dtype=torch.float64)
lam = alg64.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()[0]
mc(ratio(lam[:, n // 4: 3 * n // 4]), 0.5307, "<r>, GOE(64), middle half of the spectrum")
```

## 5.4 Pourquoi il ne faut pas mélanger les secteurs

Soit $H=\bigoplus_cH_c\in\mathcal M$ un élément autoadjoint aléatoire à **blocs indépendants**, par
exemple des blocs GOE. En tant qu'opérateur sur $\mathbb C^D$, son spectre est la réunion des
spectres des blocs. À l'intérieur d'un bloc, les niveaux se repoussent ; nous montrons que les
niveaux de blocs différents ne se repoussent pas.

Supposons que chaque bloc ait une **densité de niveaux** $\rho_c$ continue et bornée :
$\mathbb E\sum_i\varphi(\lambda_i(H_c))=\int\varphi\,\rho_c$ pour toute $\varphi$ borélienne bornée
(de sorte que $\int\rho_c=k_c$). Pour le GOE et le GUE, cela découle de la densité jointe citée
plus haut.

**Théorème 5.15 (pas de répulsion entre secteurs).** Pour $c\neq d$ et toute
$\psi:\mathbb R^2\to\mathbb R$ borélienne bornée,
$$
\mathbb E\sum_{i,j}\psi\big(\lambda_i(H_c),\lambda_j(H_d)\big)=\iint\psi(x,y)\,\rho_c(x)\rho_d(y)\,dx\,dy .
$$
Par conséquent, si $N_\varepsilon$ est le nombre de paires de niveaux de blocs différents à
distance $<\varepsilon$, alors $\mathbb EN_\varepsilon/(2\varepsilon)\to\sum_{c<d}\int\rho_c\rho_d$
quand $\varepsilon\to0$, et le plus petit espacement $\delta$ du spectre fusionné vérifie
$$
P(\delta<\varepsilon)\ \ge\ \frac{\mathbb EN_\varepsilon}{\sum_{c<d}k_ck_d}
=\frac{2\varepsilon\sum_{c<d}\int\rho_c\rho_d}{\sum_{c<d}k_ck_d}+o(\varepsilon).
$$

*Démonstration.* Les valeurs propres ordonnées sont des fonctions continues de la matrice, donc
$\Lambda_c=(\lambda_i(H_c))_i$ et $\Lambda_d$ sont des vecteurs aléatoires indépendants. D'après
Fubini (en conditionnant par $\Lambda_d$),
$\mathbb E\sum_{i,j}\psi(\Lambda_{c,i},\Lambda_{d,j})=\mathbb E\sum_j g(\Lambda_{d,j})$ avec
$g(y)=\mathbb E\sum_i\psi(\Lambda_{c,i},y)=\int\psi(x,y)\rho_c(x)\,dx$, qui est borélienne bornée ;
on applique la définition de $\rho_d$ à $g$. Avec $\psi(x,y)=1\{|x-y|<\varepsilon\}$,
$$
\frac{\mathbb EN_\varepsilon}{2\varepsilon}=\sum_{c<d}\int\rho_c(x)\Big(\frac1{2\varepsilon}\int_{x-\varepsilon}^{x+\varepsilon}\rho_d(y)\,dy\Big)dx ;
$$
la parenthèse tend vers $\rho_d(x)$ par continuité et est majorée par $\sup\rho_d$, et $\rho_c$ est
intégrable, donc la convergence dominée donne la limite. Enfin $N_\varepsilon\le\sum_{c<d}k_ck_d$,
donc $\mathbb EN_\varepsilon\le\big(\sum_{c<d}k_ck_d\big)P(N_\varepsilon\ge1)$, et
$N_\varepsilon\ge1$ entraîne $\delta<\varepsilon$. $\square$

Ainsi, pour deux blocs dont les spectres se chevauchent, $P(\delta<\varepsilon)$ est au moins
d'ordre $\varepsilon$, alors qu'à l'intérieur d'un bloc GOE $2\times2$ elle vaut
$1-e^{-\varepsilon^2/4}\le\varepsilon^2/4$ (proposition 5.13, espacement non normalisé), et pour un
bloc GOE $n\times n$ elle est d'ordre $\varepsilon^2$ (d'après la densité jointe ; non démontré
ici). Le test de répulsion des niveaux échoue pour le spectre fusionné, bien que chaque bloc soit un
GOE parfait.

*Nombreux secteurs (une dérivation, non une démonstration).* Pour des spectres stationnaires, la
probabilité $\mathcal G(s)$ de ne trouver aucun niveau dans un intervalle de longueur $s$ détermine
la densité des espacements, $p(s)=\mathcal G''(s)$ pour un espacement moyen unité (Mehta,
*op. cit.*). Des spectres indépendants se superposent en multipliant ces probabilités de lacune. Si
$C$ secteurs fournissent chacun une fraction $1/C$ des niveaux et ont chacun la probabilité de
lacune $\mathcal G_1$ dans leurs propres unités dépliées, alors $\mathcal G(s)=\mathcal G_1(s/C)^C$,
avec $\mathcal G_1(x)=1-x+\int_0^x(x-t)p_1(t)\,dt$. En dérivant deux fois en $0$,
$$
p(0)=1-\frac1C+\frac{p_1(0)}C=1-\frac1C ,
$$
l'analogue du théorème 5.15 pour des spectres dépliés, et $\mathcal G_1(s/C)^C\to e^{-s}$, la
probabilité de lacune d'un processus de Poisson : quand $C\to\infty$, le spectre fusionné ressemble
à des niveaux indépendants (N. Rosenzweig, C. E. Porter, Phys. Rev. 120 (1960), 1698–1714 ;
M. V. Berry, M. Robnik, J. Phys. A 17 (1984), 2413–2421 ; des versions rigoureuses de la limite de
Poisson sont des théorèmes limites pour les superpositions de processus ponctuels raréfiés, voir
B. Grigelionis, Theory Probab. Appl. 8 (1963)).

Dans la bibliothèque, les blocs sont séparés dès le départ : `eigenvalues()` renvoie un spectre par
secteur, et les mélanger exige un `torch.cat` explicite. Nous plaçons $C$ blocs GOE indépendants de
taille $64$ dans une même algèbre et comparons le rapport moyen par secteur à celui du spectre
fusionné :

```python
k, merged_r = 64, {}
for C in (1, 2, 4, 8, 16):
    alg_C = TypeIAlgebra([k] * C, [k] * C, complex_valued=False, precision='double')
    G = torch.randn(800 // C, C, k, k, dtype=torch.float64)
    lam = alg_C.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()    # one spectrum per sector
    per_sector = torch.stack([ratio(l[:, k // 4: 3 * k // 4]) for l in lam])
    merged = torch.cat(lam, dim=-1).sort(dim=-1).values                       # the spectrum on C^D
    r_m = ratio(merged[:, k * C // 4: 3 * k * C // 4])
    merged_r[C] = r_m.mean().item()
    print(f"C = {C:2d}: <r> per sector {per_sector.mean():.4f}, merged {merged_r[C]:.4f}"
          f" +- {r_m.std() / math.sqrt(len(r_m)):.4f}")
    assert abs(per_sector.mean() - 0.5307) < 0.01                    # every sector is GOE
assert merged_r[1] > merged_r[2] > merged_r[4] > merged_r[16]
assert abs(merged_r[16] - (2 * math.log(2) - 1)) < 0.01              # close to independent levels
```

Chaque secteur a $\langle r\rangle\approx0.53$ ; le spectre fusionné tombe à $0.42$ pour deux
secteurs et à $0.39$, la valeur des niveaux indépendants, pour $8$ et $16$.

Il en va de même dans un modèle physique. La chaîne de Heisenberg en champs aléatoires
(chapitre 9) conserve l'aimantation, de sorte que son hamiltonien appartient à
$\bigoplus_{N=0}^LM_{\binom LN}$ ; pour un désordre faible, les niveaux à l'intérieur de chaque
secteur suivent les statistiques du GOE (numériquement ; A. Pal, D. A. Huse, Phys. Rev. B 82
(2010), 174411). Ses blocs ne sont *pas* indépendants — ils dépendent tous des mêmes champs
aléatoires —, si bien que le théorème 5.15 ne s'applique pas tel quel ; l'expérience montre le même
effet.

```python
chain = SpinChain(10, complex_valued=False)                         # sectors N = 0, ..., 10
lam = chain.random_field_heisenberg(W=0.5, batch_size=20).eigenvalues()
k5 = chain.sector_dim(5)                                            # 252
r_sector = ratio(lam[5][:, k5 // 4: 3 * k5 // 4]).mean().item()
spectrum = torch.cat(lam, dim=-1).sort(dim=-1).values               # all 1024 levels mixed
r_all = ratio(spectrum[:, 256:768]).mean().item()
print(f"Heisenberg chain, L = 10: <r> in sector N = 5: {r_sector:.4f}, all sectors mixed: {r_all:.4f}")
assert r_sector > 0.49 and r_all < 0.42
```

Dans le secteur $N=5$, le rapport vaut $0.52$ (une valeur de taille finie légèrement inférieure à
$0.5307$) ; le spectre où tous les secteurs sont mélangés donne $0.39$.

> **Interprétation physique.** Des niveaux portant des nombres quantiques conservés différents
> n'interagissent pas, ils se croisent donc librement ; seuls les niveaux de même symétrie se
> repoussent. Un spectre qui mélange des secteurs de symétrie présente par conséquent des
> statistiques de type Poisson même lorsque la dynamique dans chaque secteur est pleinement
> chaotique, et une symétrie non résolue peut être prise à tort pour de l'intégrabilité ou de la
> localisation. La règle de l'art — résoudre toutes les symétries avant de calculer des
> statistiques de niveaux — revient à dire qu'il faut calculer dans l'algèbre
> $\bigoplus_cM_{k_c}$, et non dans $M_D$.

## 5.5 Points de comparaison

Unitaires de Haar : `scipy.stats.unitary_group` (la QR corrigée du théorème 5.3) et
`ortho_group`, les unitaires aléatoires de Haar de QuTiP, et toute QR de `numpy`/`torch` suivie de
la correction des phases. Calcul de Weingarten : les paquets symboliques IntU (Mathematica ;
Z. Puchała, J. A. Miszczak), RTNI (Mathematica et Python ; M. Fukuda, R. König, I. Nechita,
J. Phys. A 52 (2019), 425303) et le paquet Python Haarpy calculent symboliquement les moments de
tout ordre ; les formules démontrées dans ce chapitre en sont le cas d'ordre deux. Les statistiques
de niveaux se calculent d'ordinaire en quelques lignes de `numpy`, et les codes de diagonalisation
exacte comme QuSpin fournissent les spectres secteur par secteur. Ce qu'apporte la bibliothèque,
c'est le groupe unitaire d'une algèbre en tant qu'objet unique : des blocs de Haar indépendants de
toutes tailles tirés par des appels vectorisés, le twirl identifié à l'espérance conditionnelle
`channels.center_expectation`, et des spectres renvoyés par secteur, de sorte que mélanger les
secteurs n'est jamais le comportement par défaut.

## Exercices

1. Soit $U$ distribué selon la loi de Haar sur $U(\mathcal M)$ et
   $\operatorname{Tr}U=\sum_c\operatorname{Tr}U_c$. Montrer que $\mathbb E|\operatorname{Tr}U|^2=C$
   et que $\mathbb E|\operatorname{Tr}U|^4=\sum_cm(k_c)+2C(C-1)$ avec $m(1)=1$, $m(k)=2$ pour
   $k\ge2$ ; comparer avec $1$ et $2$ pour un unitaire de Haar de $U(D)$, $D\ge2$. Vérifier pour
   $\mathcal M=M_1\oplus M_2\oplus M_3$.
2. À l'aide du théorème 5.12, montrer que $\mathbb E|U_{11}|^{2n}=\frac{n!}{k(k+1)\cdots(k+n-1)}$
   pour $n\le k$. Calculer numériquement la matrice de Weingarten d'ordre $3$ pour $k=3$ et vérifier
   que $\mathbb E|U_{11}|^6=\frac1{10}$.
3. L'**ensemble orthogonal circulaire** (`measure='coe'`) est $S=W^TW$ avec $W$ distribué selon la
   loi de Haar sur $U(k)$. Montrer que $S$ est unitaire et symétrique, que sa loi est invariante
   par $S\mapsto V^TSV$ pour tout $V\in U(k)$, et que la densité moyenne de ses phases propres est
   uniforme. Vérifier numériquement que ses phases propres ont $\langle r\rangle\approx0.53$, comme
   le GOE.

<details>
<summary>Solutions</summary>

1. Posons $t_c=\operatorname{Tr}U_c$ ; les $t_c$ sont indépendants (corollaire 5.4). Multiplier
   $U_c$ par $e^{i\theta}$ ne change pas sa loi et multiplie $t_c$ par $e^{i\theta}$ ; un moment
   $\mathbb E[t_at_b\bar t_c\bar t_d]$ est donc nul sauf si chaque secteur apparaît aussi souvent
   dans $(a,b)$ que dans $(c,d)$, c'est-à-dire si $\{a,b\}=\{c,d\}$ en tant que multiensembles, et
   de même $\mathbb E[t_a\bar t_b]=0$ pour $a\neq b$. D'après la proposition 5.7,
   $\mathbb E[(U_c)_{ii}\overline{(U_c)_{jj}}]=\langle e_i,\mathbb E[U_ce_ie_j^*U_c^*]e_j\rangle=\delta_{ij}/k_c$,
   donc $\mathbb E|t_c|^2=1$ et $\mathbb E|\operatorname{Tr}U|^2=C$. Dans
   $\mathbb E|\sum_at_a|^4$, les termes avec $a=b$ donnent $\sum_a\mathbb E|t_a|^4$, ceux avec
   $a\neq b$ donnent $2\sum_{a\neq b}\mathbb E|t_a|^2\,\mathbb E|t_b|^2=2C(C-1)$ (deux choix de
   $(c,d)$). Pour $k=1$, $|t|=1$. Pour $k\ge2$,
   $\mathbb E|t|^4=\sum_{i,j,i',j'}\mathbb E[U_{ii}U_{jj}\overline{U_{i'i'}U_{j'j'}}]$ ; dans le
   théorème 5.10, les deux termes avec $\sigma=\pi$ imposent $\{i',j'\}=\{i,j\}$ d'une seule façon
   et donnent chacun $k^2\mathrm{Wg}(e)$ après sommation, les deux termes avec $\sigma\neq\pi$
   imposent $i=j=i'=j'$ et donnent chacun $k\,\mathrm{Wg}(\tau)$ : au total
   $2\big(\frac{k^2}{k^2-1}-\frac1{k^2-1}\big)=2$. Le cas $C=1$, $k=D$ donne $1$ et $2$ pour
   $U(D)$. Pour $M_1\oplus M_2\oplus M_3$ : $\mathbb E|\operatorname{Tr}U|^4=1+2+2+12=17$, contre
   $2$ pour un unitaire de Haar de $U(6)$.

```python
tr = U.trace                                                        # Tr U = sum_c Tr U_c, shape (B,)
mc(tr.abs() ** 2, 3.0, "E|Tr U|^2")
mc(tr.abs() ** 4, 17.0, "E|Tr U|^4")
```

2. Lorsque tous les indices valent $1$, chaque delta vaut $1$, donc
   $\mathbb E|U_{11}|^{2n}=\sum_{\sigma,\pi}\mathrm{Wg}(\sigma\pi^{-1})=\mathbf 1^TW\mathbf 1$, où
   $W=\Gamma^{-1}=(\mathrm{Wg}(\sigma^{-1}\pi))$ et $\Gamma=(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)})$.
   On a utilisé ici que $\mathrm{Wg}$ ne dépend que du type cyclique, et que
   $\sigma\pi^{-1}=\sigma(\pi^{-1}\sigma)\sigma^{-1}$ est conjuguée à
   $\pi^{-1}\sigma=(\sigma^{-1}\pi)^{-1}$, qui a le type cyclique de $\sigma^{-1}\pi$. Chaque somme
   de ligne de $\Gamma$ vaut $c_n(k)=\sum_{\rho\in S_n}k^{\#\mathrm{cycles}(\rho)}$, donc
   $\Gamma\mathbf 1=c_n(k)\mathbf 1$, $W\mathbf 1=\mathbf 1/c_n(k)$ et
   $\mathbf 1^TW\mathbf 1=n!/c_n(k)$. Enfin $c_n(x)=x(x+1)\cdots(x+n-1)$ par récurrence : une
   permutation de $\{1,\dots,n\}$ s'obtient à partir d'une permutation de $\{1,\dots,n-1\}$ soit en
   ajoutant $n$ comme point fixe (un cycle de plus, facteur $x$), soit en insérant $n$ après l'un
   des $n-1$ éléments dans son cycle ($n-1$ façons, même nombre de cycles) ; c'est une bijection,
   donc $c_n(x)=(x+n-1)c_{n-1}(x)$.

```python
import itertools

def n_cycles(p):                                                    # number of cycles of a permutation
    seen, count = set(), 0
    for i in range(len(p)):
        count += i not in seen
        while i not in seen:
            seen.add(i)
            i = p[i]
    return count

perms = list(itertools.permutations(range(3)))                      # S_3; s.index(j) = s^{-1}(j)
Gram = torch.tensor([[3.0 ** n_cycles([s.index(p[i]) for i in range(3)]) for p in perms] for s in perms],
                    dtype=torch.float64)                            # k^{#cycles(s^{-1} p)}, k = 3
Wg3 = torch.linalg.inv(Gram)                                        # the Weingarten matrix, n = k = 3
assert abs(Wg3.sum().item() - 6 / (3 * 4 * 5)) < 1e-12
mc(u3[:, 0, 0].abs() ** 6, 0.1, "E|U_11|^6")
```

3. $S^*S=W^*\bar WW^TW=1$ puisque $W^T$ est unitaire, et $S^T=W^TW=S$. Pour $V\in U(k)$, $WV$ a la
   loi de $W$ (invariance à droite, lemme 5.1), et $(WV)^T(WV)=V^TSV$. Avec $V=e^{i\theta/2}1$,
   $e^{i\theta}S$ a la loi de $S$, donc la loi de l'ensemble des phases propres est invariante par
   rotation, de même que la densité moyenne, qui est donc uniforme. En particulier, les espacements
   ne nécessitent aucun dépliement. Les statistiques locales du COE coïncident avec celles du GOE
   pour $k$ grand (F. J. Dyson, J. Math. Phys. 3 (1962), 140–156 ; Mehta, *op. cit.* ; énoncé sans
   démonstration).

```python
alg_s = TypeIAlgebra([32], [32], precision='double')
S = alg_s.random_unitary_operator(400, measure='coe').matrix[:, 0]
assert close(S, S.mT) and close(S @ S.mH, torch.eye(32, dtype=dt))
phases = torch.linalg.eigvals(S).angle().sort(dim=-1).values
s = torch.diff(torch.cat([phases, phases[:, :1] + 2 * math.pi], dim=-1), dim=-1)   # around the circle
r = torch.minimum(s, s.roll(-1, -1)) / torch.maximum(s, s.roll(-1, -1))
mc(r.mean(-1), 0.5307, "<r>, COE(32) eigenphases")
```

</details>

**Suite :** le [chapitre 6](06_quantum_channels.md) traite des canaux quantiques : représentations
de Kraus et de Stinespring, duaux relativement à une trace, et applications entre secteurs.
L'espérance conditionnelle $E$, rencontrée ici comme moyenne de Haar, en est le premier exemple.
