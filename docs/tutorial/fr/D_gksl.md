[English](../D_gksl.md) · [Español](../es/D_gksl.md) · **Français** · [Deutsch](../de/D_gksl.md) · [中文](../zh/D_gksl.md) · [日本語](../ja/D_gksl.md) · [Русский](../ru/D_gksl.md)

# Annexe D. Générateurs des semi-groupes dynamiques quantiques

Le chapitre 8 a démontré qu'un générateur de la forme de Gorini–Kossakowski–Sudarshan–Lindblad (GKSL) engendre un
semi-groupe dynamique quantique (théorème 8.6) et a utilisé la réciproque (théorème 8.7) et sa version pour
$\mathcal M=\bigoplus_cM_{k_c}$. Cette annexe démontre les deux en dimension finie, ainsi que la
liberté dans la forme GKSL et un critère d'unicité de l'état stationnaire. L'outil est la
matrice de Choi du chapitre 6 : un générateur n'est en général pas complètement positif, mais sa matrice de Choi est
positive sur l'orthogonal du vecteur maximalement intriqué. Le même type d'argument — positivité
d'une famille d'opérateurs, puis une décomposition — démontre dans la dernière section le théorème de correction
d'erreurs par algèbres d'opérateurs énoncé au chapitre 10 (théorème 10.7), avec un canal de récupération explicite.

**Vous apprendrez** pourquoi le générateur d'un semi-groupe complètement positif est conditionnellement complètement
positif, comment lire $H$ et les opérateurs de saut sur la matrice de Choi d'un générateur, en quel sens
ils sont uniques, pourquoi tout générateur sur une algèbre à secteurs a la forme utilisée à la section 8.3, quand
l'état stationnaire est unique, et quelles algèbres d'observables un canal de récupération peut protéger.

Dans toute cette annexe, $M_n=M_n(\mathbb C)$, les $E_{ij}=|i\rangle\langle j|$ sont ses unités matricielles, un
**superopérateur** est une application linéaire sur $M_n$ (ou sur $\mathcal M$), et
$$
\Omega=\sum_{i=1}^ne_i\otimes e_i\in\mathbb C^n\otimes\mathbb C^n,\qquad
J(\Psi)=\sum_{i,j=1}^nE_{ij}\otimes\Psi(E_{ij})=[\Psi(E_{ij})]_{i,j=1}^n
$$
sont le vecteur maximalement intriqué et la **matrice de Choi** d'un superopérateur $\Psi$ (théorème 6.2 ;
ici $\Psi$ n'est pas nécessairement positif). Pour des vecteurs $v,w$, $|v\rangle\langle w|$ est l'opérateur
$x\mapsto\langle w,x\rangle v$. Un superopérateur $\Psi$ **préserve l'hermiticité** si
$\Psi(A^*)=\Psi(A)^*$ pour tout $A$. Les semi-groupes dynamiques quantiques et leurs générateurs sont ceux de la
section 8.2 (lemme 8.5) : $T_t=e^{t\mathcal L}$ avec des canaux $T_t$.

## D.1 Matrices de Choi des briques élémentaires

**Lemme D.1.** Soit $\Psi$ un superopérateur sur $M_n$.

(a) $\Psi\mapsto J(\Psi)$ est une bijection linéaire de l'espace des superopérateurs sur $M_n$ sur $M_n(M_n)=M_{n^2}$.
(b) Tout $v\in\mathbb C^n\otimes\mathbb C^n$ s'écrit $v=(1\otimes K)\Omega$ pour exactement un $K\in M_n$, et
$\langle(1\otimes K)\Omega,(1\otimes L)\Omega\rangle=\operatorname{Tr}(K^*L)$. En particulier,
$\langle\Omega,(1\otimes K)\Omega\rangle=\operatorname{Tr}K$.
(c) Pour $K,L\in M_n$, le superopérateur $A\mapsto KAL^*$ a pour matrice de Choi
$|(1\otimes K)\Omega\rangle\langle(1\otimes L)\Omega|$. En particulier, $J(\mathrm{id})=|\Omega\rangle\langle\Omega|$.
(d) $\Psi$ préserve l'hermiticité si et seulement si $J(\Psi)$ est autoadjointe. Toute application CP préserve l'hermiticité.

*Démonstration.* (a) Le bloc $(i,j)$ de $J(\Psi)$ est $\Psi(E_{ij})$ et les $E_{ij}$ forment une base de $M_n$, donc
$J$ est linéaire et injective, et toute matrice par blocs $[X_{ij}]$ est la matrice de Choi du superopérateur
défini par $E_{ij}\mapsto X_{ij}$. (b) $(1\otimes K)\Omega=\sum_ie_i\otimes Ke_i$ a pour $i$-ième bloc $Ke_i$, la
$i$-ième colonne de $K$, qui détermine $K$, et tout $v$ s'obtient ainsi. Ensuite
$\langle\sum_ie_i\otimes Ke_i,\sum_je_j\otimes Le_j\rangle=\sum_i\langle Ke_i,Le_i\rangle=\operatorname{Tr}(K^*L)$ ;
pour la dernière assertion, on remplace $(K,L)$ par $(1,K)$. (c) Le bloc $(i,j)$ du membre de droite est
$Ke_i(Le_j)^*=KE_{ij}L^*$. (d) Le bloc $(i,j)$ de $J(\Psi)^*$ est $\Psi(E_{ji})^*$. Si $\Psi$
préserve l'hermiticité, $\Psi(E_{ji})^*=\Psi(E_{ji}^*)=\Psi(E_{ij})$. Réciproquement, si
$\Psi(E_{ji})^*=\Psi(E_{ij})$ pour tous $i,j$ et $A=\sum a_{ij}E_{ij}$, alors
$\Psi(A)^*=\sum\bar a_{ij}\Psi(E_{ij})^*=\sum\bar a_{ij}\Psi(E_{ji})=\Psi(A^*)$. Une application CP a une forme de Kraus
(théorème 6.2), et $(\sum_lK_lAK_l^*)^*=\sum_lK_lA^*K_l^*$. $\square$

**Lemme D.2 (décomposition).** Soient $q_1,\dots,q_C$ des vecteurs orthonormés de
$\mathbb C^n\otimes\mathbb C^n$, écrits $q_c=(1\otimes Q_c)\Omega$, soit $W=\operatorname{span}\{q_c\}$, et soit
$\Psi$ un superopérateur qui préserve l'hermiticité, avec $\langle v,J(\Psi)v\rangle\ge0$ pour tout $v\perp W$.
Posons $J=J(\Psi)$ et $\xi_c=Jq_c-\frac12\sum_d\langle q_d,Jq_c\rangle q_d$, et soit $K_c\in M_n$ défini par
$(1\otimes K_c)\Omega=\xi_c$. Alors il existe $L_1,\dots,L_r\in M_n$, où $r=\operatorname{rank}(P_W^\perp JP_W^\perp)$
et $P_W^\perp$ est le projecteur sur $W^\perp$, tels que
$$
\Psi(A)=\sum_{j=1}^rL_jAL_j^*+\sum_{c=1}^C\big(K_cAQ_c^*+Q_cAK_c^*\big)\qquad(A\in M_n),
$$
$\operatorname{Tr}(Q_c^*L_j)=0$ pour tous $c,j$, et $\operatorname{Tr}(L_j^*L_k)=0$ pour $j\neq k$.

*Démonstration.* $J$ est autoadjointe d'après le lemme D.1(d). Soient $P=\sum_c|q_c\rangle\langle q_c|$ le projecteur sur $W$
et $P^\perp=P_W^\perp=1-P$. En développant $J=(P+P^\perp)J(P+P^\perp)$ et en utilisant $PJP^\perp=PJ-PJP$ et
$P^\perp JP=JP-PJP$,
$$
J=P^\perp JP^\perp+\big(PJ+JP-PJP\big).
$$
Le premier terme est autoadjoint, et $\langle v,P^\perp JP^\perp v\rangle=\langle P^\perp v,JP^\perp v\rangle\ge0$
car $P^\perp v\perp W$ ; donc $P^\perp JP^\perp\ge0$. D'après le théorème 1.3,
$P^\perp JP^\perp=\sum_{j=1}^r|l_j\rangle\langle l_j|$, où les $l_j$ sont les vecteurs propres des valeurs propres non nulles
$\mu_j$, multipliés par $\sqrt{\mu_j}$ : ils sont deux à deux orthogonaux et appartiennent à l'image de
$P^\perp JP^\perp$, qui est orthogonale à $W$. Écrivons $l_j=(1\otimes L_j)\Omega$ ; d'après le lemme D.1(b),
$\operatorname{Tr}(Q_c^*L_j)=\langle q_c,l_j\rangle=0$ et $\operatorname{Tr}(L_j^*L_k)=\langle l_j,l_k\rangle=0$ pour $j\ne k$,
et d'après le lemme D.1(c), $P^\perp JP^\perp$ est la matrice de Choi de $A\mapsto\sum_jL_jAL_j^*$. Pour le second terme,
$JP=\sum_c|Jq_c\rangle\langle q_c|$ et $PJP=\sum_{c,d}\langle q_d,Jq_c\rangle|q_d\rangle\langle q_c|$, donc
$$
\sum_c|\xi_c\rangle\langle q_c|=JP-\tfrac12PJP,\qquad\sum_c|q_c\rangle\langle\xi_c|=PJ-\tfrac12PJP
$$
(la seconde est l'adjointe de la première ; $PJP$ est autoadjointe), et leur somme vaut $PJ+JP-PJP$. D'après le
lemme D.1(c), $|\xi_c\rangle\langle q_c|$ et $|q_c\rangle\langle\xi_c|$ sont les matrices de Choi de
$A\mapsto K_cAQ_c^*$ et $A\mapsto Q_cAK_c^*$. Ainsi les deux membres de l'identité annoncée ont pour matrice de Choi $J$,
et le lemme D.1(a) donne l'identité. $\square$

## D.2 Positivité complète conditionnelle

Un superopérateur $\Psi$ sur $M_n$ est **conditionnellement complètement positif** (CCP) si
$\langle v,J(\Psi)v\rangle\ge0$ pour tout $v\perp\Omega$. D'après le lemme D.1(b), ces $v$ sont les vecteurs
$(1\otimes K)\Omega$ avec $\operatorname{Tr}K=0$.

**Théorème D.3.** Pour un superopérateur $\mathcal L$ sur $M_n$, les assertions suivantes sont équivalentes :

1. $e^{t\mathcal L}$ est CP pour tout $t\ge0$ ;
2. $\mathcal L$ préserve l'hermiticité et est conditionnellement complètement positif ;
3. il existe une application CP $\Phi$ sur $M_n$ et $G\in M_n$ tels que $\mathcal L(A)=\Phi(A)+GA+AG^*$ pour tout $A$.

Dans 3, on peut prendre $\Phi(A)=\sum_{j=1}^rL_jAL_j^*$ avec $r\le n^2-1$, $\operatorname{Tr}L_j=0$ et
$\operatorname{Tr}(L_j^*L_k)=0$ pour $j\neq k$, et $G$ avec $\operatorname{Tr}G\in\mathbb R$.

*Démonstration.* (1)$\Rightarrow$(2). D'après le lemme 8.1, $(e^{t\mathcal L}-\mathrm{id})/t\to\mathcal L$ quand $t\to0$, et $J$ est
linéaire sur un espace de dimension finie, donc continue ; ainsi
$J(e^{t\mathcal L})=J(\mathrm{id})+t\,J\big((e^{t\mathcal L}-\mathrm{id})/t\big)$, où la dernière matrice de Choi tend vers
$J(\mathcal L)$. Pour $v\perp\Omega$, $\langle v,J(\mathrm{id})v\rangle=|\langle\Omega,v\rangle|^2=0$ (lemme D.1(c)), donc
pour $t>0$, comme $J(e^{t\mathcal L})\ge0$ (théorème 6.2),
$$
0\le\frac1t\langle v,J(e^{t\mathcal L})v\rangle=\big\langle v,J\big((e^{t\mathcal L}-\mathrm{id})/t\big)v\big\rangle
\;\xrightarrow[t\to0]{}\;\langle v,J(\mathcal L)v\rangle .
$$
Chaque $e^{t\mathcal L}$ est CP, donc préserve l'hermiticité (lemme D.1(d)), d'où
$\mathcal L(A^*)=\lim_t(e^{t\mathcal L}(A)^*-A^*)/t=\mathcal L(A)^*$.

(2)$\Rightarrow$(3). On applique le lemme D.2 avec $C=1$ et $q_1=\Omega/\sqrt n$, c'est-à-dire $Q_1=1/\sqrt n$ :
$\mathcal L(A)=\sum_jL_jAL_j^*+GA+AG^*$ avec $G=K_1/\sqrt n$, et $\Phi=\sum_jL_j\cdot L_j^*$ est CP
(théorème 6.2). De plus, $\operatorname{Tr}L_j=\sqrt n\operatorname{Tr}(Q_1^*L_j)=0$, $r\le\operatorname{rank}P_W^\perp=n^2-1$,
et d'après le lemme D.1(b),
$\operatorname{Tr}G=\langle\Omega,\xi_1\rangle/\sqrt n=\langle q_1,\xi_1\rangle=\frac12\langle q_1,J(\mathcal L)q_1\rangle\in\mathbb R$,
puisque $J(\mathcal L)$ est autoadjointe.

(3)$\Rightarrow$(1). Les étapes 1 à 4 de la démonstration du théorème 8.6 n'utilisent que le fait que le générateur est
$\mathcal L_0+\Phi$ avec $\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)$ pour une certaine matrice
$H_{\mathrm{eff}}$ et $\Phi$ CP. Avec $H_{\mathrm{eff}}=iG$, on a $\mathcal L_0(\rho)=G\rho+\rho G^*$, et ces
étapes montrent que $e^{t\mathcal L}=\sum_{m\ge0}S_m(t)$ est CP. $\square$

La positivité seule ne suffirait pas. Soient $n\ge2$, $\theta$ la transposition sur $M_n$ et $\mathcal L=\theta-\mathrm{id}$. Comme $\theta^2=\mathrm{id}$,
$$
e^{t\mathcal L}=e^{-t}\big(\cosh t\,\mathrm{id}+\sinh t\,\theta\big)=\tfrac{1+e^{-2t}}2\,\mathrm{id}+\tfrac{1-e^{-2t}}2\,\theta ,
$$
combinaison convexe d'applications positives qui préservent la trace (proposition 6.1). Mais $J(\theta)=[E_{ji}]=\sum_{i,j}E_{ij}\otimes E_{ji}$ est l'opérateur d'échange $e_k\otimes e_l\mapsto e_l\otimes e_k$, qui envoie donc
$v=e_1\otimes e_2-e_2\otimes e_1$, orthogonal à $\Omega$, sur $-v$ (comme dans la démonstration de la proposition 6.1) ; avec
$\|v\|^2=2$ et $\langle v,J(\mathrm{id})v\rangle=|\langle\Omega,v\rangle|^2=0$, on obtient
$\langle v,J(\mathcal L)v\rangle=-2$ : $\mathcal L$ n'est pas CCP, et de fait
$\langle v,J(e^{t\mathcal L})v\rangle=-(1-e^{-2t})<0$ pour tout $t>0$. Un semi-groupe d'applications positives qui préservent
la trace n'est pas nécessairement formé de canaux.

Nous vérifions le théorème D.3 sur un générateur GKSL à sauts arbitraires sur $M_3$, construit par `lindblad_superoperator`
(la matrice $S$ de $\mathcal L$ agissant sur $\operatorname{vec}\rho$ rangé ligne par ligne, de sorte que $J_{(i,a),(j,b)}=S_{(a,b),(i,j)}$ comme à la
section 6.1), et sur $\theta-\mathrm{id}$ :

```python
import math
import numpy as np
import scipy.linalg as sla
import torch
from torch_vn_algebra import TypeIAlgebra, channels, dynamics

torch.manual_seed(0)
rng = np.random.default_rng(0)
f64, c128 = torch.float64, torch.complex128

def choi_of(S, n):
    """J_{(i,a),(j,b)} = S_{(a,b),(i,j)}: the Choi matrix of the map with row-major matrix S (Section 6.1)."""
    return S.reshape(n, n, n, n).transpose(2, 0, 3, 1).reshape(n * n, n * n)

def superop(f, n):
    """Row-major matrix of a linear map f on n x n matrices."""
    E = np.eye(n * n).reshape(n * n, n, n)
    return np.stack([f(X).reshape(-1) for X in E], axis=1)

vec_of = lambda K: K.T.reshape(-1)                        # (1 (x) K) Omega = sum_i e_i (x) K e_i
op_of = lambda v, n: v.reshape(n, n).T                    # its inverse (Lemma D.1(b))
cplx = lambda *s: rng.standard_normal(s) + 1j * rng.standard_normal(s)

n = 3
M3 = TypeIAlgebra([n], [n], precision='double')
H = cplx(n, n)
H = (H + H.conj().T) / 2
Ls = [cplx(n, n) / 2 for _ in range(2)]                   # arbitrary jumps, not traceless
to_op = lambda X: M3.operator(torch.from_numpy(X).reshape(1, 1, n, n))
S = dynamics.lindblad_superoperator(M3, to_op(H), [to_op(L) for L in Ls])[0][0].numpy()
J = choi_of(S, n)
q = vec_of(np.eye(n)) / math.sqrt(n)                      # the unit vector Omega / sqrt(n)
Pperp = np.eye(n * n) - np.outer(q, q)
assert np.abs(J - J.conj().T).max() < 1e-12               # Hermiticity preserving (Lemma D.1(d))
print(f"min eigenvalue of J(L): {np.linalg.eigvalsh(J).min():.2f}; "
      f"on the complement of Omega: {np.linalg.eigvalsh(Pperp @ J @ Pperp).min():.1e}")
assert np.linalg.eigvalsh(J).min() < -0.1                 # L itself is not CP ...
assert np.linalg.eigvalsh(Pperp @ J @ Pperp).min() > -1e-12        # ... but CCP (Theorem D.3)
for t in (0.01, 0.1, 1.0):
    assert np.linalg.eigvalsh(choi_of(sla.expm(t * S), n)).min() > -1e-12

L_theta = superop(lambda X: X.T, n) - np.eye(n * n)     # L = transpose - id: Hermiticity preserving, not CCP
anti = vec_of(np.eye(n)[:, [1]] @ np.eye(n)[[0], :] - np.eye(n)[:, [0]] @ np.eye(n)[[1], :])   # e1(x)e2 - e2(x)e1
assert abs(anti @ q) < 1e-15 and abs(anti @ choi_of(L_theta, n) @ anti + 2) < 1e-12
X = cplx(n, n)
for t in (0.01, 0.1, 1.0):
    T_t = sla.expm(t * L_theta)
    b = (1 - math.exp(-2 * t)) / 2
    assert abs(anti @ choi_of(T_t, n) @ anti + (1 - math.exp(-2 * t))) < 1e-12        # <v, J(e^{tL}) v>
    assert np.abs(T_t - ((1 - b) * np.eye(n * n) + b * superop(lambda Y: Y.T, n))).max() < 1e-12
    assert abs(np.trace((T_t @ X.reshape(-1)).reshape(n, n)) - np.trace(X)) < 1e-12   # trace preserving
    assert np.linalg.eigvalsh(choi_of(T_t, n)).min() < -0.9 * b                     # but not CP
```

## D.3 Le théorème GKSL

**Théorème D.4 (Gorini–Kossakowski–Sudarshan ; Lindblad ; théorème 8.7).** Un superopérateur $\mathcal L$ sur
$M_n$ est le générateur d'un semi-groupe dynamique quantique sur $M_n$ si et seulement s'il existe $H=H^*$ et
$L_1,\dots,L_r$ dans $M_n$ tels que
$$
\mathcal L(\rho)=-i[H,\rho]+\sum_{j=1}^r\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
$$
On peut choisir $r\le n^2-1$, $\operatorname{Tr}H=0$, $\operatorname{Tr}L_j=0$ et $\operatorname{Tr}(L_j^*L_k)=0$ pour $j\neq k$.

*Démonstration.* Le sens « si » est le théorème 8.6. Réciproquement, soit $(T_t)$ un semi-groupe dynamique quantique de générateur
$\mathcal L$, de sorte que $T_t=e^{t\mathcal L}$ (lemme 8.5). Comme tout $T_t$ est CP, le théorème D.3 donne
$\mathcal L(\rho)=\sum_jL_j\rho L_j^*+G\rho+\rho G^*$ avec $r\le n^2-1$ opérateurs $L_j$ de trace nulle, deux à deux orthogonaux, et
$\operatorname{Tr}G\in\mathbb R$. Posons $H=\frac i2(G-G^*)$ et $R=\frac12(G+G^*)$, tous deux autoadjoints ; alors
$G=-iH+R$, $G^*=iH+R$ et $G\rho+\rho G^*=-i[H,\rho]+\{R,\rho\}$. De plus,
$\operatorname{Tr}H=\frac i2(\operatorname{Tr}G-\overline{\operatorname{Tr}G})=0$. Comme tout $T_t$ est TP,
$\operatorname{Tr}\mathcal L(\rho)=\frac d{dt}\operatorname{Tr}T_t(\rho)\big|_{t=0}=0$ (lemme 8.1), et par cyclicité
$$
0=\operatorname{Tr}\mathcal L(\rho)=\operatorname{Tr}\Big(\rho\,\Big(\sum_jL_j^*L_j+2R\Big)\Big)\qquad\text{for all }\rho .
$$
Avec $\rho=X=\sum_jL_j^*L_j+2R$, qui est autoadjoint, $\operatorname{Tr}(X^*X)=0$, donc $X=0$ et
$R=-\frac12\sum_jL_j^*L_j$, ce qui est la forme GKSL. $\square$

**Théorème D.5 (liberté dans la forme GKSL).** Soient $(H;L_1,\dots,L_r)$ et $(H';L'_1,\dots,L'_s)$ des données GKSL
($H,H'$ autoadjoints) de générateurs $\mathcal L$, $\mathcal L'$ comme dans le théorème D.4.

(a) $\mathcal L'=\mathcal L$ dans chacun des cas suivants :
(i) $L'_j=L_j$ et $H'=H+c1$, $c\in\mathbb R$ ;
(ii) $H'=H$ et, après avoir complété par des zéros la famille la plus courte, $L'_j=\sum_ku_{jk}L_k$ avec une matrice unitaire $u$ ;
(iii) $L'_j=L_j+a_j1$ avec $a_j\in\mathbb C$, et $H'=H-\frac i2\sum_j(\bar a_jL_j-a_jL_j^*)$, qui est autoadjoint.

(b) Si $\operatorname{Tr}H=\operatorname{Tr}H'=0$, si tous les $L_j,L'_j$ sont de trace nulle et si $\mathcal L=\mathcal L'$, alors $H=H'$ et,
après complétion par des zéros, $L'_j=\sum_ku_{jk}L_k$ avec une matrice unitaire $u$.

(c) Par conséquent, deux données GKSL ont le même générateur si et seulement si l'une s'obtient de l'autre par un nombre fini
d'opérations (i)–(iii). En particulier, l'opérateur $H$ de la forme normalisée du théorème D.4 et sa
partie dissipative $\Phi(\rho)=\sum_jL_j\rho L_j^*$ sont déterminés par $\mathcal L$.

*Démonstration.* (a) (i) $[c1,\rho]=0$. (ii) D'après le théorème 6.7, $\Phi'=\sum_jL'_j\cdot L_j'^*$ est égale à
$\Phi=\sum_jL_j\cdot L_j^*$ ; la duale ne dépend pas des opérateurs de Kraus (section 6.2), donc
$\sum_jL_j'^*L'_j=\Phi^*(1)=\sum_jL_j^*L_j$. (iii) Pour un seul saut $L$ et $a\in\mathbb C$, le développement donne
$$
(L+a)\rho(L+a)^*-\tfrac12\{(L+a)^*(L+a),\rho\}-\Big(L\rho L^*-\tfrac12\{L^*L,\rho\}\Big)=\tfrac12[B,\rho],
\qquad B=\bar aL-aL^* ,
$$
car les termes $\bar aL\rho+a\rho L^*$ du premier produit se combinent avec
$-\frac12(aL^*\rho+a\rho L^*+\bar aL\rho+\bar a\rho L)$ et les termes en $|a|^2$ se compensent. Comme $B^*=-B$, la
matrice $\frac i2B$ est autoadjointe, donc $H'=H-\frac i2\sum_jB_j$ est autoadjoint, et
$-i[H'-H,\rho]=-\frac12\sum_j[B_j,\rho]$ compense la variation du dissipateur.

(b) Soient $\Phi=\sum_jL_j\cdot L_j^*$ et $G=-iH-\frac12\Phi^*(1)$, de sorte que $\mathcal L(\rho)=\Phi(\rho)+G\rho+\rho G^*$.
D'après le lemme D.1(c),
$J(\mathcal L)=\sum_j|l_j\rangle\langle l_j|+|g\rangle\langle\Omega|+|\Omega\rangle\langle g|$ avec $l_j=(1\otimes L_j)\Omega$ et
$g=(1\otimes G)\Omega$, et $l_j\perp\Omega$ car $\operatorname{Tr}L_j=0$ (lemme D.1(b)). En notant $P^\perp$ le projecteur
sur $\Omega^\perp$, on a donc $P^\perp J(\mathcal L)P^\perp=\sum_j|l_j\rangle\langle l_j|=J(\Phi)$, et de même
$P^\perp J(\mathcal L')P^\perp=J(\Phi')$. Ainsi $J(\Phi)=J(\Phi')$, $\Phi=\Phi'$ (lemme D.1(a)), et le théorème 6.7 fournit la
matrice unitaire $u$. Alors $\Phi^*(1)=\Phi'^*(1)$, donc $[H-H',\rho]=0$ pour tout $\rho$. Une matrice $X$ qui commute avec toutes
les $E_{ij}$ est un multiple de $1$ (comparer les coefficients de $XE_{ij}$ et de $E_{ij}X$), donc $H-H'=\lambda1$, et
$\lambda=0$ en prenant les traces.

(c) D'après (a), les opérations (i)–(iii) ne changent pas le générateur, et l'inverse de chaque opération est une opération du
même type ((iii) avec $-a_j$, puisque $\bar a_jL'_j-a_jL_j'^*=\bar a_jL_j-a_jL_j^*$). Toute donnée se ramène à la
forme normalisée par (iii) avec $a_j=-\operatorname{Tr}L_j/n$, suivie de (i) avec $c=-\operatorname{Tr}H''/n$, où $H''$ est
l'hamiltonien après (iii) ; deux données normalisées de même générateur sont reliées par (ii), d'après (b). $\square$

La démonstration du théorème D.4 est un algorithme : $\Phi$ à partir de la compression de $J(\mathcal L)$ à $\Omega^\perp$, puis
$G$ à partir de $\xi_1$, puis $H$. Nous l'exécutons sur le générateur ci-dessus, reconstruisons le générateur avec la bibliothèque, et
comparons avec la forme normalisée des données d'origine calculée par le théorème D.5(iii) et (i) :

```python
def gksl_from_generator(S, n, tol=1e-10):
    """(H, [L_j]) with Tr H = Tr L_j = 0 from the matrix of a GKSL generator: the proof of Theorem D.4."""
    J = choi_of(S, n)
    q = vec_of(np.eye(n)) / math.sqrt(n)
    Pperp = np.eye(n * n) - np.outer(q, q)
    w, U = np.linalg.eigh(Pperp @ J @ Pperp)              # = J(Phi) >= 0
    jumps = [math.sqrt(x) * op_of(U[:, k], n) for k, x in enumerate(w) if x > tol]
    xi = J @ q - 0.5 * (q.conj() @ J @ q) * q
    G = op_of(xi, n) / math.sqrt(n)
    return 0.5j * (G - G.conj().T), jumps

H_rec, L_rec = gksl_from_generator(S, n)
S_rec = dynamics.lindblad_superoperator(M3, to_op(H_rec), [to_op(L) for L in L_rec])[0][0].numpy()
assert np.abs(S_rec - S).max() < 1e-12                    # the same generator
assert abs(np.trace(H_rec)) < 1e-12 and all(abs(np.trace(L)) < 1e-12 for L in L_rec)

a = [-np.trace(L) / n for L in Ls]                        # Theorem D.5(iii): make the jumps traceless
L_norm = [L + aj * np.eye(n) for L, aj in zip(Ls, a)]
H_norm = H - 0.5j * sum(np.conj(aj) * L - aj * L.conj().T for L, aj in zip(Ls, a))
H_norm = H_norm - np.trace(H_norm) / n * np.eye(n)        # Theorem D.5(i)
assert np.abs(H_norm - H_norm.conj().T).max() < 1e-12
assert np.abs(H_rec - H_norm).max() < 1e-12               # H is unique (Theorem D.5(b))
jump_choi = lambda Ks: sum(np.outer(vec_of(K), vec_of(K).conj()) for K in Ks)
assert np.abs(jump_choi(L_rec) - jump_choi(L_norm)).max() < 1e-12      # so is sum_j L_j . L_j^*
Wr, Wn = np.stack([vec_of(K) for K in L_rec], 1), np.stack([vec_of(K) for K in L_norm], 1)
u = np.linalg.lstsq(Wr, Wn, rcond=None)[0]                # L_norm_j = sum_k u_kj L_rec_k
assert np.abs(Wr @ u - Wn).max() < 1e-12 and np.abs(u.conj().T @ u - np.eye(2)).max() < 1e-12   # unitary
print("eigenvalues of the normalised H:", np.round(np.linalg.eigvalsh(H_rec), 4), " jumps:", len(L_rec))
```

**Théorème D.6 (unicité de l'état stationnaire).** Soit $\mathcal L$ le générateur GKSL du théorème D.4
de données $(H;L_1,\dots,L_r)$. Si $\mathcal L$ a un état stationnaire fidèle $\sigma$ ($\sigma>0$) et si les seules
matrices qui commutent avec $H$ et avec tous les $L_j$, $L_j^*$ sont les multiples de $1$, alors $\ker\mathcal L=\mathbb C\sigma$ :
$\sigma$ est le seul état stationnaire. De plus, $\ker\mathcal L^*=\mathbb C1$.

*Démonstration.* La duale pour l'accouplement par la trace, $\operatorname{Tr}(\mathcal L(\rho)A)=\operatorname{Tr}(\rho\,\mathcal L^*(A))$, est
$\mathcal L^*(A)=i[H,A]+\sum_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})$ : par cyclicité,
$\operatorname{Tr}(-i[H,\rho]A)=\operatorname{Tr}(\rho\,i[H,A])$, $\operatorname{Tr}(L\rho L^*A)=\operatorname{Tr}(\rho L^*AL)$ et
$\operatorname{Tr}(\{L^*L,\rho\}A)=\operatorname{Tr}(\rho\{L^*L,A\})$ (c'est l'exercice 1 du chapitre 8).
Un développement direct, dans lequel les termes hamiltoniens se compensent parce que $A\mapsto i[H,A]$ est une dérivation, donne
pour tout $X$
$$
\mathcal L^*(X^*X)-X^*\mathcal L^*(X)-\mathcal L^*(X^*)X=\sum_j[L_j,X]^*[L_j,X]:
$$
pour un seul saut, parmi les six termes d'anticommutateur, quatre se compensent deux à deux et les deux autres,
$\frac12X^*L^*LX+\frac12X^*L^*LX$, s'additionnent, de sorte que le membre de gauche vaut $L^*X^*XL-X^*L^*XL-L^*X^*LX+X^*L^*LX$, c'est-à-dire
$(XL-LX)^*(XL-LX)$. Soit $X\in\ker\mathcal L^*$. Alors $\mathcal L^*(X^*)=\mathcal L^*(X)^*=0$, car $\mathcal L^*$
préserve l'hermiticité, donc $\mathcal L^*(X^*X)=\sum_j[L_j,X]^*[L_j,X]=:Y\ge0$, et
$\operatorname{Tr}(\sigma Y)=\operatorname{Tr}(\mathcal L(\sigma)X^*X)=0$. Comme $\sigma>0$,
$\operatorname{Tr}(\sigma Y)=\operatorname{Tr}(\sigma^{1/2}Y\sigma^{1/2})$, avec $\sigma^{1/2}Y\sigma^{1/2}\ge0$, ne s'annule que si
$\sigma^{1/2}Y\sigma^{1/2}=0$, c'est-à-dire $Y=0$ ; en prenant les traces, $\sum_j\operatorname{Tr}([L_j,X]^*[L_j,X])=0$ avec des termes
positifs, donc $[L_j,X]=0$ pour tout $j$. Le même argument appliqué à $X^*\in\ker\mathcal L^*$
donne $[L_j,X^*]=0$, c'est-à-dire $[L_j^*,X]=0$. Alors $X$ commute avec $L_j^*L_j$, la partie dissipative de $\mathcal L^*(X)$
s'annule, et $0=\mathcal L^*(X)=i[H,X]$. Par hypothèse, $X\in\mathbb C1$ ; réciproquement, $\mathcal L^*(1)=0$. Donc
$\ker\mathcal L^*=\mathbb C1$. Dans les bases $(E_{ij})$ et $(E_{ji})$ de $M_n$, qui sont duales pour l'accouplement
$(\rho,A)\mapsto\operatorname{Tr}(\rho A)$, la matrice de $\mathcal L^*$ est la transposée de celle de $\mathcal L$ ; elles ont donc le
même rang, et $\dim\ker\mathcal L=\dim\ker\mathcal L^*=1$. Un état de $\mathbb C\sigma$ de trace $1$ est $\sigma$. $\square$

C'est le cas de dimension finie d'un critère de Frigerio [*Commun. Math. Phys.* 63, 269 (1978)],
utilisé à la section 8.5. Les deux hypothèses sont nécessaires. Dans l'exercice 4 du chapitre 8, l'état $\tfrac12\,1$ est fidèle
et stationnaire, mais $H$ et $\sigma_z$ commutent avec toutes les matrices diagonales, et les états stationnaires engendrent
l'espace de dimension deux des matrices diagonales. Inversement, sur $M_3$, soient $H=0$, $L_1=|e_1\rangle\langle e_3|$,
$L_2=|e_2\rangle\langle e_3|$. La $*$-algèbre engendrée par $L_1,L_2$ contient $L_1L_2^*=E_{12}$, $L_1^*L_1=E_{33}$ et donc
toutes les unités matricielles, de sorte que seuls les multiples de $1$ commutent avec les $L_j,L_j^*$ ; mais tout $\rho=P\rho P$ avec $P=E_{11}+E_{22}$
est stationnaire, puisque $L_j\rho=0=\rho L_j^*$ et $\{E_{33},\rho\}=0$, donc $\dim\ker\mathcal L\ge4$ ; d'après le théorème D.6, aucun état
stationnaire n'est fidèle. Pour le générateur aléatoire ci-dessus, le code vérifie les hypothèses et la conclusion
(la matrice de $\mathcal L^*$ sur les vecteurs de $A^T$ rangés ligne par ligne est $S^T$: $\operatorname{vec}(\mathcal L^*(A)^T)=S^T\operatorname{vec}(A^T)$), ainsi que le second exemple :

```python
def kernel_dim(S, tol=1e-9):
    return int((np.linalg.svd(S, compute_uv=False) < tol).sum())

def commutant_dim(ops, n):                    # dim {X : XA = AX for all A in ops}
    return kernel_dim(np.concatenate([np.kron(A, np.eye(n)) - np.kron(np.eye(n), A.T) for A in ops]))

w, V = np.linalg.eig(S)
sigma = V[:, np.argmin(abs(w))].reshape(n, n)             # row-major: vec(X)_{(a,b)} = X_ab
sigma = sigma / np.trace(sigma)
assert np.abs(S @ sigma.reshape(-1)).max() < 1e-12        # stationary
assert np.abs(sigma - sigma.conj().T).max() < 1e-10 and np.linalg.eigvalsh(sigma).min() > 1e-3   # faithful
assert commutant_dim([H] + Ls + [L.conj().T for L in Ls], n) == 1
assert kernel_dim(S) == 1 and kernel_dim(S.T) == 1        # Theorem D.6
assert np.abs(S.T @ np.eye(n).reshape(-1)).max() < 1e-12  # ker L^* = C 1

e3 = np.eye(3)
L_abs = [np.outer(e3[0], e3[2]), np.outer(e3[1], e3[2])]       # L_1 = |e_1><e_3|, L_2 = |e_2><e_3|
S_abs = superop(lambda R: sum(L @ R @ L.T - 0.5 * (L.T @ L @ R + R @ L.T @ L) for L in L_abs), 3)
assert commutant_dim(L_abs + [L.T for L in L_abs], 3) == 1
print("second example: dim ker L =", kernel_dim(S_abs))
assert kernel_dim(S_abs) == 4                             # no faithful stationary state, no uniqueness
```

## D.4 Algèbres à secteurs

Soit $\mathcal M=\bigoplus_{c=1}^CM_{k_c}\subseteq M_D$ l'algèbre des matrices diagonales par blocs, $D=\sum_ck_c$, de
projecteurs centraux minimaux $P_c$, et soit $E(X)=\sum_cP_cXP_c$ l'espérance conditionnelle de $M_D$ sur
$\mathcal M$. Comme dans la démonstration de la proposition 8.8, $E(B)=B$, $E(XB)=E(X)B$ et $E(BX)=BE(X)$ pour $B\in\mathcal M$,
et $E(X^*)=E(X)^*$. Pour une application CP $\Phi:\mathcal M\to\mathcal M$, $\Phi^*$ désigne sa duale brute (proposition 6.8).

**Théorème D.7.** Soit $\mathcal L$ le générateur d'un semi-groupe dynamique quantique $(T_t)$ sur $\mathcal M$. Il existe
une application CP $\Phi:\mathcal M\to\mathcal M$ et $H=H^*\in\mathcal M$ tels que
$$
\mathcal L(\rho)=\Phi(\rho)-\tfrac12\{\Phi^*(1),\rho\}-i[H,\rho]\qquad(\rho\in\mathcal M).
$$
D'après la proposition 6.4, $\Phi(\rho)_d=\sum_{c,l}K^{dc}_l\rho_cK^{dc*}_l$ avec des opérateurs de Kraus entre les secteurs ;
réciproquement, tout $\mathcal L$ de cette forme engendre un semi-groupe dynamique quantique sur $\mathcal M$ (proposition 8.8(1)).
Ainsi les générateurs de la section 8.3 (où $\Phi=\sum_jg_jL_j\cdot L_j^*+\sum_mg_m\Phi_m$) sont exactement les
générateurs des semi-groupes dynamiques quantiques sur $\mathcal M$.

*Démonstration.* *Étape 1 : une application CP sur $M_D$.* D'après la proposition 6.4,
$T_t(\rho)_d=\sum_{c,l}K^{dc}_l(t)\rho_cK^{dc}_l(t)^*$. Soit $N^{dc}_l(t)\in M_D$ la matrice ayant le bloc $K^{dc}_l(t)$ en
position $(d,c)$ et des zéros ailleurs. Pour $X\in M_D$, $N^{dc}_lXN^{dc*}_l$ est le bloc
$K^{dc}_lX_{cc}K^{dc*}_l$ en position $(d,d)$, où $X_{cc}$ est le $c$-ième bloc diagonal, c'est-à-dire le $c$-ième bloc
de $E(X)$. Par conséquent, $\hat T_t=T_t\circ E=\sum_{d,c,l}N^{dc}_l\cdot N^{dc*}_l$ est CP sur $M_D$ (théorème 6.2), avec
$\hat T_0=E$. D'après le lemme 8.5, $(\hat T_t-E)/t=\big((T_t-\mathrm{id})/t\big)\circ E\to\mathcal L\circ E$.

*Étape 2 : positivité hors des vecteurs des secteurs.* Soient $p_c=(1\otimes P_c)\Omega\in\mathbb C^D\otimes\mathbb C^D$ et
$W=\operatorname{span}\{p_c\}$. D'après le lemme D.1(c), $J(E)=\sum_c|p_c\rangle\langle p_c|$, donc pour $v\perp W$ et $t>0$,
$0\le\frac1t\langle v,J(\hat T_t)v\rangle=\langle v,J((\hat T_t-E)/t)v\rangle\to\langle v,J(\mathcal L\circ E)v\rangle$.
Tout $T_t$ préserve l'hermiticité (par sa forme de Kraus), donc $\mathcal L$ aussi (comme dans le théorème D.3), ainsi que $E$
et $\mathcal L\circ E$.

*Étape 3 : décomposition.* D'après le lemme D.1(b), $\langle p_c,p_d\rangle=\operatorname{Tr}(P_cP_d)=\delta_{cd}k_c$, donc les
$q_c=p_c/\sqrt{k_c}=(1\otimes Q_c)\Omega$, avec $Q_c=P_c/\sqrt{k_c}$, sont orthonormés. Le lemme D.2 appliqué à
$\Psi=\mathcal L\circ E$ fournit une application CP $\hat\Phi=\sum_j\hat L_j\cdot\hat L_j^*$ sur $M_D$ et des $K_c\in M_D$ tels que
$$
\mathcal L(E(A))=\hat\Phi(A)+\sum_c\big(K_cAP_c+P_cAK_c^*\big)/\sqrt{k_c}\qquad(A\in M_D).
$$
Pour $A\in\mathcal M$, $E(A)=A$ et $AP_c=P_cA$, donc la somme vaut $\tilde GA+A\tilde G^*$ avec
$\tilde G=\sum_cK_cP_c/\sqrt{k_c}$.

*Étape 4 : compression.* Pour $A\in\mathcal M$, $\mathcal L(A)\in\mathcal M$, donc
$\mathcal L(A)=E(\mathcal L(A))=E(\hat\Phi(A))+GA+AG^*$ avec $G=E(\tilde G)\in\mathcal M$. L'application
$\Phi=E\circ\hat\Phi$ sur $\mathcal M$ s'écrit
$\Phi(A)_d=\sum_{c,j}(P_d\hat L_jP_c)A_c(P_d\hat L_jP_c)^*$ : elle est CP d'après la proposition 6.4, et ses opérateurs de Kraus sont les
blocs des $\hat L_j$.

*Étape 5 : préservation de la trace.* Comme dans le théorème D.4, $\operatorname{Tr}\mathcal L(A)=0$ pour $A\in\mathcal M$, et
$\operatorname{Tr}\Phi(A)=\operatorname{Tr}(A\Phi^*(1))$, donc $\operatorname{Tr}(AX)=0$ pour tout $A\in\mathcal M$, où
$X=\Phi^*(1)+G+G^*\in\mathcal M$ est autoadjoint. Avec $A=X$, on obtient $X=0$ (la trace brute est fidèle, théorème 3.1).
Avec $H=\frac i2(G-G^*)\in\mathcal M$, on a $G=-iH-\frac12\Phi^*(1)$ et $G\rho+\rho G^*=-i[H,\rho]-\frac12\{\Phi^*(1),\rho\}$. $\square$

Les vecteurs $p_c$ remplacent $\Omega=\sum_cp_c$ parce que $\hat T_0=E$ et non $\mathrm{id}$ : $J(\mathcal L\circ E)$ n'est pas nécessairement
positive sur $\Omega^\perp$, et dans l'exemple ci-dessous elle ne l'est pas.

**Corollaire D.8 (représentation de Heisenberg).** Soit $(S_t)_{t\ge0}$ une famille continue d'applications CP unifères sur
$\mathcal M$ avec $S_0=\mathrm{id}$ et $S_{s+t}=S_sS_t$. Son générateur est $A\mapsto\Psi(A)+K^*A+AK$ avec
$\Psi$ CP sur $\mathcal M$, $K\in\mathcal M$ et $\Psi(1)+K+K^*=0$.

*Démonstration.* Soient $T_t=S_t^*$ (duales brutes). Elles sont CP (proposition 6.8) et TP (corollaire 6.9(b)),
$T_0=\mathrm{id}$, et $T_{s+t}=(S_tS_s)^*=S_s^*S_t^*=T_sT_t$, car
$\operatorname{Tr}(S_tS_s(x)y)=\operatorname{Tr}(xS_s^*S_t^*(y))$ et $S_tS_s=S_{s+t}$. Comme $\Psi\mapsto\Psi^*$ est linéaire sur un
espace de dimension finie, $t\mapsto T_t$ est continue et $(T_t-\mathrm{id})/t\to\mathcal L_S^*$, où $\mathcal L_S$ est le
générateur de $(S_t)$. Le théorème D.7 donne $\mathcal L_S^*(\rho)=\Phi(\rho)+G\rho+\rho G^*$ avec $G+G^*=-\Phi^*(1)$, et par
cyclicité, $\operatorname{Tr}((G\rho+\rho G^*)A)=\operatorname{Tr}(\rho(G^*A+AG))$, donc d'après le corollaire 6.9(c),
$\mathcal L_S(A)=\Phi^*(A)+G^*A+AG$. On prend $\Psi=\Phi^*$ (CP d'après la proposition 6.8) et $K=G$. $\square$

C'est le cas de dimension finie du théorème de Christensen et Evans cité après le théorème 8.7. Nous
vérifions le théorème D.7 sur $M_1\oplus M_2\subseteq M_3$ pour un générateur de `lindblad_rhs` avec un saut à l'intérieur des
secteurs et un saut aléatoire entre eux, en suivant les cinq étapes de la démonstration :

```python
alg = TypeIAlgebra([1, 2], [1, 2], precision='double')    # M = M_1 (+) M_2 in M_3
D, ks = 3, alg.k_factors
P = [np.diag([1.0, 0, 0]), np.diag([0, 1.0, 1.0])]
E = lambda X: sum(Pc @ X @ Pc for Pc in P)               # the conditional expectation onto M
def blocks(X):                                            # dense block-diagonal 3 x 3 -> library layout
    out = torch.zeros(1, 2, 2, 2, dtype=c128)
    out[0, 0, :1, :1] = torch.from_numpy(X[:1, :1])
    out[0, 1] = torch.from_numpy(X[1:, 1:])
    return out
dense = lambda Y: sla.block_diag(Y[0, 0, :1, :1].numpy(), Y[0, 1].numpy())
Hs = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
inside = alg.operator(torch.randn(1, 2, 2, 2, dtype=c128))       # a jump inside the sectors
between = channels.random_inter_sector_channel(alg, alg, kraus_rank=2)
rhs = dynamics.lindblad_rhs(Hs, [inside, between], rates=[0.7, 1.3], algebra=alg)
L_M = lambda X: dense(rhs(0.0, blocks(X)))               # the generator on M (X block diagonal)
J = choi_of(superop(lambda X: L_M(E(X)), D), D)          # J(L o E) on M_3
p = [vec_of(Pc) for Pc in P]
Q = np.stack([pc / np.linalg.norm(pc) for pc in p], 1)    # orthonormal q_c = p_c / sqrt(k_c)
Pw = np.eye(D * D) - Q @ Q.T                              # projection onto the complement of W
assert np.linalg.eigvalsh(Pw @ J @ Pw).min() > -1e-12     # Step 2
w, U = np.linalg.eigh(Pw @ J @ Pw)                        # Step 3 (Lemma D.2)
jumps = [math.sqrt(x) * op_of(U[:, k], D) for k, x in enumerate(w) if x > 1e-10]
xis = [J @ Q[:, c] - 0.5 * sum((Q[:, d] @ J @ Q[:, c]) * Q[:, d] for d in range(2)) for c in range(2)]
G = E(sum(op_of(xi, D) @ Pc / math.sqrt(k) for xi, Pc, k in zip(xis, P, ks)))   # Step 4
Phi = lambda X: E(sum(K @ X @ K.conj().T for K in jumps))
Phi_dual_one = E(sum(K.conj().T @ K for K in jumps))      # Phi^*(1)
H_M = 0.5j * (G - G.conj().T)
assert np.abs(G + G.conj().T + Phi_dual_one).max() < 1e-12      # Step 5
for _ in range(3):
    X = E(cplx(D, D))
    gksl = Phi(X) - 0.5 * (Phi_dual_one @ X + X @ Phi_dual_one) - 1j * (H_M @ X - X @ H_M)
    assert np.abs(L_M(X) - gksl).max() < 1e-12            # Theorem D.7
    dual = lambda A: E(sum(K.conj().T @ A @ K for K in jumps)) + G.conj().T @ A + A @ G   # Corollary D.8
    assert abs(np.trace(L_M(X) @ (Y := E(cplx(D, D)))) - np.trace(X @ dual(Y))) < 1e-12
Pomega = np.eye(D * D) - np.outer(sum(p), sum(p)) / D     # projection onto the complement of Omega
print(f"numerical observation: min eigenvalue of J(L o E) on the complement of Omega: "
      f"{np.linalg.eigvalsh(Pomega @ J @ Pomega).min():.2f}")
assert np.linalg.eigvalsh(Pomega @ J @ Pomega).min() < -0.1
```

## D.5 Correction d'erreurs par algèbres d'opérateurs

Nous reprenons le cadre de la section 10.3. $\Phi(\rho)=\sum_{i=1}^rE_i\rho E_i^*$ est un canal sur $M_D$, $\mathcal C=P\mathbb C^D$
un code, $B(\mathcal C)=\{X\in M_D:X=PXP\}$ l'ensemble des opérateurs sur $\mathcal C$, $\mathfrak A\subseteq B(\mathcal C)$ une $*$-algèbre
d'unité $P$ (une algèbre de von Neumann sur $\mathcal C$ d'après le théorème 2.2), et $\mathfrak A'$ son commutant dans
$B(\mathcal C)$. $\mathfrak A$ est **corrigeable** pour $\Phi$ s'il existe un canal $\mathcal R$ sur $M_D$ tel que
$P(\mathcal R\circ\Phi)^*(A)P=A$ pour tout $A\in\mathfrak A$. Si $\mathcal R$ a pour opérateurs de Kraus les $R_m$, alors
$\mathcal R\circ\Phi$ a pour opérateurs de Kraus les $R_mE_i$, et $\sum_mR_m^*R_m=1$.

**Théorème D.9 (Bény–Kempf–Kribs ; théorème 10.7).** $\mathfrak A$ est corrigeable pour $\Phi$ si et seulement si
$PE_i^*E_jP\in\mathfrak A'$ pour tous $i,j$. Dans ce cas, le canal $\mathcal R$ construit dans la démonstration corrige $\mathfrak A$.

*Démonstration.* *Condition nécessaire.* Supposons que $\mathcal R$ corrige $\mathfrak A$ et posons $W_{mi}=R_mE_iP$, de sorte que
$\sum_{m,i}W_{mi}^*XW_{mi}=P(\mathcal R\circ\Phi)^*(X)P$ pour tout $X\in M_D$. Ceci vaut $X$ pour $X\in\mathfrak A$, et $P$ pour
$X=1$, car $(\mathcal R\circ\Phi)^*$ est unifère (corollaire 6.9(a)). Soit $p\in\mathfrak A$ un projecteur ; alors $p\le P$
et $\sum W_{mi}^*(1-p)W_{mi}=P-p$. Pour $x=px$,
$\sum_{m,i}\|(1-p)W_{mi}x\|^2=\langle x,(P-p)x\rangle=0$, et pour $x=(P-p)x$,
$\sum_{m,i}\|pW_{mi}x\|^2=\langle x,px\rangle=0$. Par conséquent $(1-p)W_{mi}p=0$ et $pW_{mi}(P-p)=0$, et comme
$W_{mi}=W_{mi}P$,
$$
W_{mi}p=pW_{mi}p=pW_{mi}P=pW_{mi}.
$$
Tout $A\in\mathfrak A$ autoadjoint est une combinaison linéaire de projecteurs de $\mathfrak A$ : d'après le théorème 1.3 sur $\mathcal C$,
$A=\sum_k\lambda_kp_k$ avec des $\lambda_k$ distincts et des projecteurs orthogonaux $p_k$, $\sum_kp_k=P$, et
$p_k=\prod_{l\neq k}(A-\lambda_lP)/(\lambda_k-\lambda_l)\in\mathfrak A$ (appliquer les deux membres aux vecteurs propres dans $\mathcal C$ et à $\mathcal C^\perp$ ; le produit vide vaut $P$). Tout élément
de $\mathfrak A$ s'écrit $A_1+iA_2$ avec $A_1,A_2\in\mathfrak A$ autoadjoints. Ainsi $W_{mi}A=AW_{mi}$ pour tout $A\in\mathfrak A$ et,
en passant aux adjoints, $W_{mi}^*A=AW_{mi}^*$. Comme $\sum_mW_{mi}^*W_{mj}=PE_i^*(\sum_mR_m^*R_m)E_jP=PE_i^*E_jP$,
$$
A\,PE_i^*E_jP=\sum_mW_{mi}^*AW_{mj}=PE_i^*E_jP\,A\qquad(A\in\mathfrak A),
$$
c'est-à-dire $PE_i^*E_jP\in\mathfrak A'$.

*Condition suffisante.* D'après le théorème 2.6 appliqué à $\mathfrak A$ sur $\mathcal C$, il existe une base orthonormée $|c;a,b\rangle$
($c=1,\dots,C$, $a=1,\dots,k_c$, $b=1,\dots,m_c$) de $\mathcal C$ — les images des $e_a\otimes f_b$ du $c$-ième terme par
l'unitaire $V$ — dans laquelle $\mathfrak A$ est formée des opérateurs $\sum_{c,a,a',b}(A_c)_{aa'}|c;a,b\rangle\langle c;a',b|$,
$A_c\in M_{k_c}$, et $\mathfrak A'$ des opérateurs $\sum_{c,a,b,b'}(Y_c)_{bb'}|c;a,b\rangle\langle c;a,b'|$, $Y_c\in M_{m_c}$.
L'hypothèse fournit donc des matrices $\Lambda^c_{ij}\in M_{m_c}$ telles que
$$
\langle c;a,b|E_i^*E_j|d;a',b'\rangle=\delta_{cd}\,\delta_{aa'}\,(\Lambda^c_{ij})_{bb'} .
$$
*Étape 1.* Pour chaque $c$ et chaque couple $\alpha=(i,b)$, soient $F^c_\alpha:\mathbb C^{k_c}\to\mathbb C^D$,
$F^c_\alpha e_a=E_i|c;a,b\rangle$, et $\Lambda^c_{\alpha\beta}=(\Lambda^c_{ij})_{bb'}$ pour $\beta=(j,b')$. La formule ci-dessus s'écrit
$F^{c*}_\alpha F^d_\beta=\delta_{cd}\Lambda^c_{\alpha\beta}1_{k_c}$ : une condition de Knill–Laflamme pour les « erreurs » $F^c_\alpha$ sur le
code $\mathbb C^{k_c}$, avec en outre l'orthogonalité entre les secteurs.

*Étape 2.* Comme dans la démonstration du théorème 10.5, $\Lambda^c$ est autoadjointe, car
$\Lambda^c_{\beta\alpha}1=F^{c*}_\beta F^c_\alpha=(F^{c*}_\alpha F^c_\beta)^*=\bar\Lambda^c_{\alpha\beta}1$, et positive, car
$\sum\bar z_\alpha z_\beta\Lambda^c_{\alpha\beta}1=(\sum z_\alpha F^c_\alpha)^*(\sum z_\beta F^c_\beta)\ge0$. Écrivons
$\Lambda^c=U^c\operatorname{diag}(d^c_\kappa)U^{c*}$ avec $U^c$ unitaire et $d^c_\kappa\ge0$, et
$G^c_\kappa=\sum_\alpha u^c_{\alpha\kappa}F^c_\alpha$. Alors $G^{c*}_\kappa G^d_\lambda=\delta_{cd}\delta_{\kappa\lambda}d^c_\kappa1$, donc
$G^c_\kappa=0$ si $d^c_\kappa=0$, et sinon $V^c_\kappa=G^c_\kappa/\sqrt{d^c_\kappa}$ est une isométrie avec
$V^{c*}_\kappa V^d_\lambda=\delta_{cd}\delta_{\kappa\lambda}1$. En inversant, $F^c_\alpha=\sum_\kappa\bar u^c_{\alpha\kappa}G^c_\kappa
=\sum_{\kappa:d^c_\kappa>0}\bar u^c_{\alpha\kappa}\sqrt{d^c_\kappa}\,V^c_\kappa$, d'où
$$
V^{c'*}_\kappa F^d_\beta=\delta_{c'd}\,\bar u^d_{\beta\kappa}\sqrt{d^d_\kappa}\,1,\qquad
F^{c*}_\alpha V^{c'}_\kappa=\delta_{cc'}\,u^c_{\alpha\kappa}\sqrt{d^c_\kappa}\,1 .
$$
*Étape 3 : la récupération.* Les $V^c_\kappa V^{c*}_\kappa$ ($d^c_\kappa>0$) sont des projecteurs deux à deux orthogonaux ; soit $Q$ leur
somme. Soient $\iota_c:\mathbb C^{k_c}\to\mathbb C^D$, $\iota_ce_a=|c;a,1\rangle$, une isométrie, et
$$
\mathcal R(\sigma)=\sum_{c}\sum_{\kappa:d^c_\kappa>0}\iota_cV^{c*}_\kappa\,\sigma\,V^c_\kappa\iota_c^*+(1-Q)\sigma(1-Q).
$$
C'est un canal : $\sum_{c,\kappa}V^c_\kappa\iota_c^*\iota_cV^{c*}_\kappa+(1-Q)^2=Q+1-Q=1$. Il mesure le secteur $c$ et
le syndrome $\kappa$, annule l'erreur sur $\mathbb C^{k_c}$ et écrit le résultat dans la copie $b=1$.

*Étape 4 : correction.* $\mathcal R^*(X)=\sum_{c,\kappa}V^c_\kappa\iota_c^*X\iota_cV^{c*}_\kappa+(1-Q)X(1-Q)$, et pour
$A\in\mathfrak A$ de blocs $A_c$, $\iota_c^*A\iota_c=A_c$. Comme chaque $F^d_\beta$ est une combinaison des $V^d_\kappa$,
$(1-Q)F^d_\beta=0$. Avec l'étape 2, pour $\alpha=(i,b)$, $\beta=(j,b')$,
$$
\langle c;a,b|E_i^*\mathcal R^*(A)E_j|d;a',b'\rangle=\sum_{c',\kappa}\langle e_a,F^{c*}_\alpha V^{c'}_\kappa A_{c'}V^{c'*}_\kappa F^d_\beta e_{a'}\rangle
=\delta_{cd}\sum_\kappa u^c_{\alpha\kappa}d^c_\kappa\bar u^c_{\beta\kappa}(A_c)_{aa'}=\delta_{cd}(\Lambda^c_{ij})_{bb'}(A_c)_{aa'} .
$$
Prenons $j=i$ et sommons sur $i$ : d'après la formule de l'hypothèse avec $\sum_iE_i^*E_i=1$, $\sum_i\Lambda^c_{ii}=1_{m_c}$, donc
$$
\langle c;a,b|(\mathcal R\circ\Phi)^*(A)|d;a',b'\rangle=\delta_{cd}\,\delta_{bb'}(A_c)_{aa'}=\langle c;a,b|A|d;a',b'\rangle .
$$
$P(\mathcal R\circ\Phi)^*(A)P$ et $A$ appartiennent tous deux à $B(\mathcal C)$ et ont la même matrice dans la base $|c;a,b\rangle$ ; ils sont donc
égaux. $\square$

**Corollaire D.10.** (a) Pour $\mathfrak A=B(\mathcal C)$, $\mathfrak A'=\mathbb CP$, et le théorème D.9 est le théorème de
Knill–Laflamme 10.5 ($C=1$, $m_1=1$ ; l'étape 3 est la récupération du théorème 10.5). (b) Pour $\mathfrak A=\operatorname{span}\{P_1,\dots,P_C\}$ avec
des projecteurs orthogonaux $P_c\ne0$, $\sum_cP_c=P$, c'est l'équivalence de 2 et 3 dans la proposition 10.6. (c) (*Sous-systèmes
sans bruit.*) Avec les notations de la proposition 10.1(a), soient $\mathcal C=V(\mathbb C^{k_c}\otimes\mathbb C^{m_c})$ le $c$-ième
terme pour l'algèbre d'interaction $\mathcal A$ de $\Phi$, et $\mathfrak A=V(1_{k_c}\otimes M_{m_c})V^*$ sur $\mathcal C$, une
algèbre de la forme $M_{m}\otimes1_{k}$ après échange des deux facteurs. Alors $\mathfrak A$ est corrigée par
$\mathcal R=\mathrm{id}$, et $PE_i^*E_jP=V(E_{i,c}^*E_{j,c}\otimes1_{m_c})V^*\in\mathfrak A'$ ($P=VV^*$).

*Démonstration.* (a) $\mathfrak A'=\mathbb CP$ est le cas $C=1$, $m_1=1$ de la description de $\mathfrak A'$ ci-dessus, et
$PE_i^*E_jP\in\mathbb CP$ est la condition de Knill–Laflamme. La corrigeabilité de $B(\mathcal C)$ est la correction au sens de la
section 10.2 : si $\mathcal R(\Phi(\rho))=\rho$ pour tous les états $\rho$ à support dans $\mathcal C$, alors
$\operatorname{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(A)P)=\operatorname{Tr}(\rho A)$ pour chacun d'eux, et ces $\rho$ engendrent $B(\mathcal C)$ ; la
réciproque est montrée à la section 10.3. Avec $C=1$, $m_1=1$, les opérateurs de Kraus $\iota_1V^{1*}_\kappa$ de l'étape 3 sont ceux
du théorème 10.5, car $\iota_1\iota_1^*=P$. (b) $T\in B(\mathcal C)$ commute avec tous les $P_c$ si et seulement si
$P_cTP_d=0$ pour $c\ne d$, ce qui est la condition 2 de la proposition 10.6. (c) Tout $A\in\mathfrak A$, prolongé par $0$ aux
autres termes, appartient à $\mathcal A'=V(\bigoplus_c1_{k_c}\otimes M_{m_c})V^*$ (théorème 2.6), donc $\Phi^*(A)=A$ d'après la proposition
10.1(a), et $PAP=A$. La formule pour $PE_i^*E_jP$ découle de $V^*E_iV=\bigoplus_cE_{i,c}\otimes1_{m_c}$, et
$E_{i,c}^*E_{j,c}\otimes1_{m_c}$ commute avec $1_{k_c}\otimes M_{m_c}$. $\square$

D'après le théorème 2.6, toute $\mathfrak A$ est $*$-isomorphe à une algèbre $\bigoplus_cM_{k_c}$. Lorsque $C\ge2$ et que $k_c\ge2$ pour au moins un $c$, une
$\mathfrak A$ corrigeable protège l'étiquette de secteur $c$ (information classique) en même temps que l'information quantique de
chaque $M_{k_c}$ ; de tels codes sont appelés codes **hybrides**. Comme exemple, nous prenons $\mathfrak A=(M_2\otimes1_2)\oplus M_1$ sur un
code aléatoire de dimension cinq dans $\mathbb C^8$ et un bruit qui agit sur les facteurs de jauge $\mathbb C^{m_c}$ par
$\sum_iB_{c,i}^*B_{c,i}=1_{m_c}$, envoie les secteurs sur des sous-espaces orthogonaux, brouille $\mathcal C^\perp$, puis est mélangé par un
unitaire aléatoire (théorème 6.7). La condition de Knill–Laflamme est en défaut, la condition du théorème D.9 est vérifiée, et la
récupération de la démonstration restaure tout $A\in\mathfrak A$, mais pas un opérateur générique sur $\mathcal C$ :

```python
D, dims = 8, [(2, 2), (1, 1)]                             # A = (M_2 (x) 1_2) (+) (M_1 (x) 1_1), dim C = 5
U8 = np.linalg.qr(cplx(D, D))[0]
Vc = np.linalg.qr(cplx(D, 5))[0]                          # orthonormal basis |c; a, b> of the code
W = [Vc[:, :4], Vc[:, 4:]]                                # sector c: columns indexed by (a, b), row-major
Pcode = sum(Wc @ Wc.conj().T for Wc in W)
Y = [U8[:, :4], U8[:, 4:6]]                               # orthogonal images C^{k_c} (x) C^{s_c}, s_c = 2
r = 3
B = [np.linalg.qr(cplx(2 * r, m))[0].reshape(r, 2, m) for _, m in dims]   # sum_i B_i^* B_i = 1_{m_c}
E_code = [sum(Yc @ np.kron(np.eye(k), Bc[i]) @ Wc.conj().T for Yc, Bc, Wc, (k, m) in zip(Y, B, W, dims))
          for i in range(r)]
kraus = E_code + [np.linalg.qr(cplx(D, D))[0] @ (np.eye(D) - Pcode)]
mix = np.linalg.qr(cplx(r + 1, r + 1))[0]                 # Theorem 6.7: an equivalent family
kraus = [sum(mix[i, l] * K for l, K in enumerate(kraus)) for i in range(r + 1)]
assert np.abs(sum(K.conj().T @ K for K in kraus) - np.eye(D)).max() < 1e-12      # a channel

units = []                                                 # matrix units of A
for Wc, (k, m) in zip(W, dims):
    for a1 in range(k):
        for a2 in range(k):
            units.append(Wc @ np.kron(np.eye(k)[:, [a1]] @ np.eye(k)[[a2], :], np.eye(m)) @ Wc.conj().T)
comm = lambda X: max(np.abs(X @ A - A @ X).max() for A in units)
assert max(comm(Pcode @ Ki.conj().T @ Kj @ Pcode) for Ki in kraus for Kj in kraus) < 1e-12   # Theorem D.9
KL_dev = max(np.abs(Pcode @ Ki.conj().T @ Kj @ Pcode - np.trace(Pcode @ Ki.conj().T @ Kj) / 5 * Pcode).max()
             for Ki in kraus for Kj in kraus)
print(f"Knill-Laflamme condition violated by {KL_dev:.2f}")
assert KL_dev > 0.05

R_kraus, Qsum = [], np.zeros((D, D), complex)             # Steps 1-3 of the proof
for Wc, (k, m) in zip(W, dims):
    F = [Ki @ Wc[:, b::m] for Ki in kraus for b in range(m)]   # F_(i,b) e_a = E_i |c; a, b>
    Lam = np.array([[(Fa.conj().T @ Fb)[0, 0] for Fb in F] for Fa in F])
    assert max(np.abs(Fa.conj().T @ Fb - Lam[x, y] * np.eye(k)).max()
               for x, Fa in enumerate(F) for y, Fb in enumerate(F)) < 1e-12      # Step 1
    d, u = np.linalg.eigh(Lam)
    iota = Wc[:, 0::m]                                    # e_a -> |c; a, 1>
    for kk in range(len(d)):
        if d[kk] > 1e-12:
            Vk = sum(u[x, kk] * Fx for x, Fx in enumerate(F)) / math.sqrt(d[kk])
            R_kraus.append(iota @ Vk.conj().T)
            Qsum += Vk @ Vk.conj().T
R_kraus.append(np.eye(D) - Qsum)
assert np.abs(sum(R.conj().T @ R for R in R_kraus) - np.eye(D)).max() < 1e-10    # R is a channel
dual = lambda X: sum(Ki.conj().T @ R.conj().T @ X @ R @ Ki for R in R_kraus for Ki in kraus)
for _ in range(3):
    A = sum(complex(z) * Ue for z, Ue in zip(cplx(len(units)), units))
    assert np.abs(Pcode @ dual(A) @ Pcode - A).max() < 1e-10                       # A is corrected
B_gen = Pcode @ cplx(D, D) @ Pcode
print(f"a generic operator on the code: error {np.abs(Pcode @ dual(B_gen) @ Pcode - B_gen).max():.2f}")
assert np.abs(Pcode @ dual(B_gen) @ Pcode - B_gen).max() > 0.1
assert max(comm(R @ Ki @ Pcode) for R in R_kraus for Ki in kraus) < 1e-10   # proof of necessity: W_mi in A'
```

Les résultats de cette annexe sont dus à V. Gorini, A. Kossakowski et E. C. G. Sudarshan
[*J. Math. Phys.* 17, 821 (1976)] et G. Lindblad [*Commun. Math. Phys.* 48, 119 (1976)] (théorèmes D.3–D.5),
E. Christensen et D. E. Evans [*J. London Math. Soc.* 20, 358 (1979)] (corollaire D.8 dans des algèbres de von Neumann
générales), A. Frigerio (théorème D.6), et C. Bény, A. Kempf et D. W. Kribs [*Phys. Rev. Lett.* 98, 100502
(2007) ; *Phys. Rev. A* 76, 042303 (2007)] (théorème D.9) ; les démonstrations ci-dessus sont complètes pour les énoncés
de dimension finie utilisés dans le cours.
