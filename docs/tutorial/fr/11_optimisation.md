[English](../11_optimisation.md) · [Español](../es/11_optimisation.md) · **Français** · [Deutsch](../de/11_optimisation.md) · [中文](../zh/11_optimisation.md) · [日本語](../ja/11_optimisation.md) · [Русский](../ru/11_optimisation.md)

# 11. Optimisation et inégalités de traces

Une inégalité de traces est un énoncé sur les valeurs extrêmes d'une fonctionnelle d'opérateurs. La section 3.7 en
a présenté une qui est plausible mais fausse, ainsi qu'un test par échantillonnage qui ne détecte pas la violation.
Ce chapitre explique pourquoi l'échantillonnage échoue en grande dimension et le remplace par une optimisation sur
des ensembles de contraintes paramétrés de façon exacte ; lorsque les extrêmes sont connus sous forme close, nous
les démontrons, de sorte que l'optimiseur est confronté à un théorème.

**Vous apprendrez** comment la bibliothèque paramètre les unitaires et les spectres de contraste prescrit, et
pourquoi tout opérateur admissible est atteint ; l'inégalité de traces de von Neumann et les extrêmes exacts
qu'elle fournit ; une inégalité de traces avec son cas d'égalité ; pourquoi les échantillons de Haar restent loin des
extrêmes ; et ce que démontre une recherche multi-départs.

Dans tout ce chapitre, $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ avec $D=\sum_ck_c\ge2$, $\operatorname{Tr}$ est
la trace brute (section 3.1) et $\|A\|_1=\sum_c\|A_c\|_1$, la norme trace de $A$ vue comme matrice diagonale par blocs
dans $M_D(\mathbb C)$. Les valeurs singulières de $A\in M_k(\mathbb C)$ sont $s_1(A)\ge\dots\ge s_k(A)$ ;
$x^\downarrow$ et $x^\uparrow$ désignent les réarrangements décroissant et croissant d'un vecteur réel $x$. $\Delta$
est le contraste de Michelson (section 3.6), et pour $\delta\in[0,1)$
$$
\ell(\delta)=\frac{1-\delta}{1+\delta},\qquad 1-\ell(\delta)=\frac{2\delta}{1+\delta};
$$
un opérateur positif avec $\lambda_{\max}=1$ a pour contraste $\ell(\lambda_{\min})$ ; il a donc le contraste $\delta$
si et seulement si $\lambda_{\min}=\ell(\delta)$, puisque $\ell$ est une involution strictement décroissante de
$[0,1]$. Les exemples sont les trois fonctionnelles de l'article qui accompagne la bibliothèque (section *From
sampling to optimisation*) :

| | fonctionnelle | variables |
|---|---|---|
| Exp. 1 | $z=\lvert\operatorname{Tr}(XUY)\rvert-\operatorname{Tr}(XY)$ | $X,Y\ge0$, $U\in\mathcal M$ unitaire |
| Exp. 2 | $z=\operatorname{Tr}\lvert XY\rvert-\operatorname{Tr}(\lvert X\rvert\,\lvert Y\rvert)$ | $X=UX_0$, $Y=VY_0$, $X_0,Y_0\ge0$, $U,V$ unitaires |
| Exp. 3 | $z=\operatorname{Tr}(Y\lvert X\rvert Y)-\operatorname{Tr}\lvert YXY\rvert$ | $X=X^*$, $Y\ge0$ |

avec $|X|$, $|Y|$ normalisés par $\lambda_{\max}=1$ et leurs contrastes comme paramètres. Les **enveloppes** sont
$\overline z(\delta)=\sup\{z:\Delta(|X|)=\delta\}$ et $\underline z(\delta)=\inf\{z:\Delta(|X|)=\delta\}$, toutes les
autres variables étant libres.

## 11.1 Ensembles de contraintes et leurs paramétrages

**Définition.** Pour $\delta\in[0,1)$ et $\ell=\ell(\delta)$, la **couche de contraste** est
$$
\mathcal P_\delta=\{X\in\mathcal M:X\ge0,\ \lambda_{\max}(X)=1,\ \Delta(X)=\delta\}
=\{X=X^*:\operatorname{spec}X\subseteq[\ell,1],\ 1,\ell\in\operatorname{spec}X\}.
$$
Les **emplacements** de $\mathcal M$ sont les $D$ couples $(c,i)$, $i\le k_c$ ; pour $\lambda\in\mathbb R^D$,
$\operatorname{diag}(\lambda)=\bigoplus_c\operatorname{diag}(\lambda_{c,1},\dots,\lambda_{c,k_c})$.

`optimize.PositiveParam`, avec un contraste prescrit $\delta$, tire pour chaque départ deux emplacements distincts
$s_\pm=(c_\pm,i_\pm)$ et un unitaire $U_0\in\mathcal M$ selon la loi de Haar, et envoie $\theta\in\mathbb R^D$ et
$G=\bigoplus_cG_c$ sur
$$
X(\theta,G)=U\operatorname{diag}(\lambda)U^*,\qquad U=U_0\,e^{G-G^*},\qquad
\lambda_{s_+}=1,\quad\lambda_{s_-}=\ell,\quad\lambda_s=\ell+(1-\ell)\tfrac{1+\sin\theta_s}2\ \text{otherwise}.
$$
`UnitaryParam` est l'application $G\mapsto U_0e^{G-G^*}$ seule. Contrairement à une sigmoïde, $(1+\sin\theta)/2$
atteint $0$ et $1$ pour des $\theta$ finis, ce qui importe car les spectres extrémaux se trouvent souvent en $\ell$
et $1$. Un contraste libre dans $[a,b]$ est un paramètre de plus, $\delta=a+(b-a)(1+\sin\eta)/2$.

**Lemme 11.1.** Pour tout unitaire $U_0\in\mathcal M$, l'application $G\mapsto U_0e^{G-G^*}$ est surjective sur le
groupe unitaire de $\mathcal M$.

*Démonstration.* Soit $U$ unitaire. D'après le théorème 1.3 (bloc par bloc) et le corollaire 1.4(2),
$U_0^*U=\bigoplus_cW_c\operatorname{diag}(e^{i\varphi_{c,j}})_jW_c^*$ avec $W_c$ unitaire et $\varphi_{c,j}\in\mathbb R$.
Alors $S=\bigoplus_cW_c\operatorname{diag}(i\varphi_{c,j})_jW_c^*$ vérifie $S^*=-S$ et $e^S=U_0^*U$ (la série
$\sum_nS^n/n!$ se calcule terme à terme dans la base propre), et $G=S/2$ donne $G-G^*=S$. $\square$

**Proposition 11.2 (le paramétrage est surjectif).** Pour $\delta$, $s_\pm$ et $U_0$ fixés, l'image de
$(\theta,G)\mapsto X(\theta,G)$ est $\{X\in\mathcal P_\delta:1\in\operatorname{spec}X_{c_+},\ \ell\in\operatorname{spec}X_{c_-}\}$,
et $\mathcal P_\delta$ est la réunion de ces images sur tous les couples d'emplacements distincts.

*Démonstration.* L'image est contenue dans cet ensemble : $X$ a des valeurs propres $\lambda_s\in[\ell,1]$, avec
$\lambda_{s_\pm}$ dans les blocs $c_\pm$. Réciproquement, diagonalisons $X_c=W_c\operatorname{diag}(\mu_c)W_c^*$ et
permutons chaque $\mu_c$ (en remplaçant $W_c$ par $W_c\Pi_c$, où $\Pi_c$ est une matrice de permutation) de sorte
que $\mu_{s_+}=1$ et $\mu_{s_-}=\ell$. C'est possible : si $c_+=c_-$ et $\delta>0$, alors $1\ne\ell$ sont deux valeurs
propres du même bloc et occupent des positions différentes ; si $\delta=0$, toutes les valeurs propres valent $1$.
Si $\ell<1$, posons $\theta_s=\arcsin\big(2(\mu_s-\ell)/(1-\ell)-1\big)$ pour les autres emplacements, de sorte que
$\lambda=\mu$ (si $\ell=1$, tout $\theta$ convient). Le lemme 11.1 fournit $G$ avec $U_0e^{G-G^*}=\bigoplus_cW_c$.
*Réunion :* pour $X\in\mathcal P_\delta$, choisissons des blocs $c_\pm$ contenant les valeurs propres $1$, $\ell$ et
des emplacements distincts dans ces blocs ; si $c_+=c_-$ et $\delta>0$, le bloc a deux valeurs propres différentes,
donc $k_{c_+}\ge2$ ; si $\delta=0$, $X=1$ et deux emplacements quelconques conviennent. $\square$

Ainsi, des emplacements aléatoires permettent à une recherche multi-départs de couvrir toute la couche, et la
contrainte est satisfaite exactement à chaque itéré. `SelfAdjointParam` utilise $X=U\operatorname{diag}(\sigma\lambda)U^*$
avec des signes $\sigma\in\{\pm1\}^D$ tirés pour chaque départ (et modifiés uniquement par les changements de signe
de la section 11.5), de sorte que $|X|=U\operatorname{diag}(\lambda)U^*$, et (par la même démonstration) la réunion sur
les signes et les emplacements est $\{X=X^*:|X|\in\mathcal P_\delta\}$. Seul $S=G-G^*$ compte. L'application
$S\mapsto U_0e^S$ sur les $S$ antiautoadjoints a pour différentielle en $S=0$ l'application $S\mapsto U_0S$, une
bijection sur l'espace tangent $U_0\cdot\{S^*=-S\}$ du groupe unitaire ; c'est donc un difféomorphisme local au
voisinage de $S=0$ ; loin de $0$, elle n'est pas injective, et sa différentielle est singulière là où deux valeurs
propres de $S$ diffèrent d'un multiple non nul de $2\pi i$. C'est pourquoi `extremize` **recentre** après chaque
tour ($U_0\leftarrow U_0e^{G-G^*}$, $G\leftarrow0$). La carte de Cayley $S\mapsto U_0C(-S/2)$ (exercice 1) et la
rétraction polaire $S\mapsto$ facteur unitaire de $U_0(1+S)$ (théorème 1.8) ont la même différentielle en $0$.

Le code vérifie la contrainte, ainsi que la proposition 11.2 de manière constructive : les paramètres de la
démonstration reproduisent une cible.

```python
import math
import scipy.linalg
import torch
from torch_vn_algebra import TypeIAlgebra, optimize as opt

torch.manual_seed(0)
alg = TypeIAlgebra([3, 3], [3, 3], precision='double')          # M_3 (+) M_3, D = 6
ell = lambda delta: (1 - delta) / (1 + delta)                     # lambda_min for contrast delta
d = torch.tensor([0.0, 0.25, 0.5, 0.9], dtype=torch.float64)
x = opt.PositiveParam(alg, 4, delta=d).operator()                 # four starts, one contrast each
assert torch.allclose(x.michelson_contrast, d, atol=1e-12)        # the constraint holds exactly
assert torch.allclose(x.lambda_max, torch.ones(4, dtype=torch.float64), atol=1e-12)

P = opt.PositiveParam(alg, 1, delta=torch.tensor([0.5], dtype=torch.float64))
lo = ell(0.5)
mu = lo + (1 - lo) * torch.rand(1, 2, 3, dtype=torch.float64)    # target spectrum in [lo, 1] ...
(cp, ip), (cm, im) = P.spectrum.pin_max[0].tolist(), P.spectrum.pin_min[0].tolist()
mu[0, cp, ip], mu[0, cm, im] = 1.0, lo                            # ... with 1 and lo in the drawn slots
W = alg.random_unitary_operator(1).matrix                         # target eigenvectors
target = W @ torch.diag_embed(mu.to(W.dtype)) @ W.mH
with torch.no_grad():                                             # the parameters of the proof
    P.spectrum.theta.copy_(torch.asin((2 * (mu - lo) / (1 - lo) - 1).clamp(-1, 1)))
    V = P.unitary.base.mH @ W                                     # U0^* W = exp(S), G = S / 2
    S = torch.stack([torch.from_numpy(scipy.linalg.logm(V[0, c].numpy())) for c in range(2)])
    P.unitary.G.copy_(S[None] / 2)
err = (P.operator().matrix - target).abs().max().item()
print(f"slots of 1 and lo: {(cp, ip)}, {(cm, im)};  reconstruction error {err:.1e}")
assert err < 1e-10
```

> **Interprétation physique.** Le contraste de Michelson est la visibilité $(I_{\max}-I_{\min})/(I_{\max}+I_{\min})$
> de franges d'interférence. Pour $0\le X\le1$ avec $\lambda_{\max}=1$ — un filtre, un atténuateur, l'efficacité d'un
> détecteur —, il fixe la dynamique $\lambda_{\min}/\lambda_{\max}$, l'inverse d'un conditionnement.

## 11.2 L'inégalité de traces de von Neumann

Une matrice $(c_{ij})$ avec $c_{ij}\ge0$ dont toutes les sommes de lignes et de colonnes valent $1$ est dite
**bistochastique** ; pour un unitaire $G$, $(|G_{ij}|^2)$ et $(|G_{ji}|^2)$ sont bistochastiques, car les lignes et
les colonnes de $G$ sont des vecteurs unitaires.

**Lemme 11.3 (réarrangement).** Soient $p_1\ge\dots\ge p_k$, $q_1\ge\dots\ge q_k$ réels et $(c_{ij})$
bistochastique. Alors $\sum_ip_iq_{k+1-i}\le\sum_{i,j}p_iq_jc_{ij}\le\sum_ip_iq_i$.

*Démonstration.* *Majoration.* Pour $m\le k$, soit $r_j=\sum_{i\le m}c_{ij}\in[0,1]$, de sorte que $\sum_jr_j=m$.
Comme $q_j\ge q_m$ pour $j\le m$ et $q_j\le q_m$ pour $j>m$,
$$
\sum_{j\le m}q_j-\sum_jq_jr_j=\sum_{j\le m}q_j(1-r_j)-\sum_{j>m}q_jr_j\ge q_m\Big(m-\sum_jr_j\Big)=0 .
$$
Avec $x_i=\sum_jq_jc_{ij}$, cela signifie $\sum_{i\le m}x_i\le\sum_{j\le m}q_j$, avec égalité pour $m=k$. Par
sommation d'Abel et $p_m-p_{m+1}\ge0$,
$$
\sum_ip_ix_i=p_k\sum_{i\le k}x_i+\sum_{m<k}(p_m-p_{m+1})\sum_{i\le m}x_i
\le p_k\sum_{j\le k}q_j+\sum_{m<k}(p_m-p_{m+1})\sum_{j\le m}q_j=\sum_ip_iq_i .
$$
*Minoration.* Appliquons la majoration à $p$, à la suite décroissante $q'_j=-q_{k+1-j}$ et à la matrice
bistochastique $c'_{ij}=c_{i,k+1-j}$ :
$-\sum_{ij}p_iq_jc_{ij}=\sum_{ij}p_iq'_jc'_{ij}\le\sum_ip_iq'_i=-\sum_ip_iq_{k+1-i}$. $\square$

Toute $A\in M_k(\mathbb C)$ admet une **décomposition en valeurs singulières** $A=W\Sigma V^*$, avec $W,V$ unitaires
et $\Sigma=\operatorname{diag}(s_i(A))$ : écrivons $A=W'|A|$ (théorème 1.8), $|A|=V\Sigma V^*$ (théorème 1.3) et
$W=W'V$.

**Théorème 11.4 (inégalité de traces de von Neumann).** Pour $A,B\in M_k(\mathbb C)$,
$|\operatorname{Tr}(AB)|\le\sum_is_i(A)\,s_i(B)$.

*Démonstration.* Soient $A=W_A\Sigma_AV_A^*$, $B=W_B\Sigma_BV_B^*$, $a_i=s_i(A)$, $b_j=s_j(B)$, et $G=V_A^*W_B$,
$H=V_B^*W_A$. Par cyclicité, $\operatorname{Tr}(AB)=\operatorname{Tr}(\Sigma_AG\Sigma_BH)=\sum_{i,j}a_iG_{ij}b_jH_{ji}$,
donc, d'après $|g||h|\le\frac12(|g|^2+|h|^2)$ et $a_ib_j\ge0$,
$$
|\operatorname{Tr}(AB)|\le\sum_{i,j}a_ib_j\,d_{ij},\qquad d_{ij}=\tfrac12\big(|G_{ij}|^2+|H_{ji}|^2\big).
$$
$(d_{ij})$ est bistochastique, et le lemme 11.3 donne $\sum_{ij}a_ib_jd_{ij}\le\sum_ia_ib_i$. $\square$

**Corollaire 11.5.** (a) Pour $P,Q\in M_k(\mathbb C)$ autoadjointes de valeurs propres $p,q$ :
$\sum_ip^\downarrow_iq^\uparrow_i\le\operatorname{Tr}(PQ)\le\sum_ip^\downarrow_iq^\downarrow_i$.
(b) Pour $P,Q\ge0$ : $\operatorname{Tr}(PQ)\le\|PQ\|_1\le\sum_ip^\downarrow_iq^\downarrow_i$.

*Démonstration.* (a) Comme dans la démonstration du théorème 3.5, $\operatorname{Tr}(PQ)=\sum_{ij}p_iq_jc_{ij}$ avec la
matrice bistochastique $c_{ij}=|\langle u_i,v_j\rangle|^2$ pour des bases propres ordonnées de façon décroissante ;
appliquons le lemme 11.3 (sa minoration est $\sum_ip^\downarrow_iq^\downarrow_{k+1-i}=\sum_ip^\downarrow_iq^\uparrow_i$).
(b) La première inégalité est le théorème 1.9. D'après le théorème 1.9 et le théorème 11.4,
$\|PQ\|_1=\max_U|\operatorname{Tr}(P\,QU)|\le\sum_is_i(P)s_i(QU)$, et $s_i(QU)=q^\downarrow_i$ car
$|QU|=(U^*Q^2U)^{1/2}=U^*QU$ (théorème 1.6). $\square$

La violation de la proposition 3.9 est donc au plus $\sum_ix^\downarrow_i(y^\downarrow_i-y^\uparrow_i)$. Pour Exp. 2,
ce nombre est la réponse exacte :

**Proposition 11.6 (Exp. 2 à spectres fixés).** Dans Exp. 2, $|X|=X_0$, $|Y|=Y_0$ et
$z=\|X_0VY_0\|_1-\operatorname{Tr}(X_0Y_0)$. Si les spectres $x_c$ de $(X_0)_c$ et $y_c$ de $(Y_0)_c$ sont fixés et
que les bases propres et $V$ varient, alors
$$
\sup z=-\inf z=\sum_c\sum_ix^\downarrow_{c,i}\big(y^\downarrow_{c,i}-y^\uparrow_{c,i}\big),
$$
et les deux sont atteints.

*Démonstration.* $X^*X=X_0U^*UX_0=X_0^2$, donc $|X|=X_0$ (théorème 1.6) ; de même $|Y|=Y_0$, et $|UA|=|A|$ donne
$\operatorname{Tr}|XY|=\|X_0VY_0\|_1$. Toutes les variables agissent bloc par bloc et $z$ est une somme sur les blocs ;
on peut donc supposer $C=1$. Posons $Y'=VY_0V^*\ge0$ (de spectre $y$) ; comme $|AV|=V^*|A|V$,
$\|X_0VY_0\|_1=\|X_0Y'V\|_1=\|X_0Y'\|_1$. Le corollaire 11.5 appliqué aux couples $(X_0,Y')$ et $(X_0,Y_0)$, avec
$\|X_0Y'\|_1\ge\operatorname{Tr}(X_0Y')$, donne
$$
\sum_ix^\downarrow_iy^\uparrow_i-\sum_ix^\downarrow_iy^\downarrow_i\le\|X_0Y'\|_1-\operatorname{Tr}(X_0Y_0)\le\sum_ix^\downarrow_iy^\downarrow_i-\sum_ix^\downarrow_iy^\uparrow_i .
$$
Soit $J$ la permutation de retournement, $Je_i=e_{k+1-i}$ ; alors $J\operatorname{diag}(w)=\operatorname{diag}(w_k,\dots,w_1)J$,
et la norme trace de $\operatorname{diag}(w)J$ vaut $\sum_i|w_i|$. Avec $X_0=\operatorname{diag}(x^\downarrow)$ et $V=J$ :
pour $Y_0=\operatorname{diag}(y^\uparrow)$, $X_0VY_0=\operatorname{diag}(x^\downarrow_iy^\downarrow_i)J$ et
$\operatorname{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\uparrow_i$ (la borne supérieure) ; pour $Y_0=\operatorname{diag}(y^\downarrow)$,
$X_0VY_0=\operatorname{diag}(x^\downarrow_iy^\uparrow_i)J$ et $\operatorname{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\downarrow_i$
(la borne inférieure). $\square$

**Corollaire 11.7 (enveloppe d'Exp. 2).** Soient $\ell_X,\ell_Y\in(0,1]$, et supposons $k_c\ge2$ pour au moins un $c$. Sur tous les
$X_0,Y_0\ge0$ avec $\lambda_{\max}=1$, $\lambda_{\min}(X_0)=\ell_X$, $\lambda_{\min}(Y_0)=\ell_Y$, et tous les unitaires,
$$
\sup z=-\inf z=\sum_c\Big\lfloor\frac{k_c}2\Big\rfloor(1-\ell_X)(1-\ell_Y).
$$

*Démonstration.* Dans un bloc, écrivons $x_i=x^\downarrow_{c,i}$, $y_i=y^\downarrow_{c,i}$, $k=k_c$, de sorte que
$y^\uparrow_{c,i}=y_{k+1-i}$. En associant $i$ à $k+1-i$,
$x_i(y_i-y_{k+1-i})+x_{k+1-i}(y_{k+1-i}-y_i)=(x_i-x_{k+1-i})(y_i-y_{k+1-i})$, et pour $k$ impair l'indice du milieu
contribue $x_i(y_i-y_i)=0$. Par conséquent
$$
\sum_ix_i(y_i-y_{k+1-i})=\sum_{i\le\lfloor k/2\rfloor}(x_i-x_{k+1-i})(y_i-y_{k+1-i})\le\Big\lfloor\frac k2\Big\rfloor(1-\ell_X)(1-\ell_Y),
$$
puisque $x_i\in[\ell_X,1]$, $y_i\in[\ell_Y,1]$. L'égalité a lieu si, dans chaque bloc, les $\lfloor k_c/2\rfloor$ plus
grandes valeurs propres valent $1$ et les $\lfloor k_c/2\rfloor$ plus petites valent $\ell_X$ (resp. $\ell_Y$), celle
du milieu étant arbitraire ; ces spectres ont pour extrêmes globaux $1$ et $\ell_X$ ($\ell_Y$) parce que $k_c\ge2$ pour au moins un $c$.
La proposition 11.6 permet de conclure. $\square$

Ici $1-\ell_X=2\delta/(1+\delta)$ pour $\Delta(|X|)=\delta$ ; pour $k_c$ impair, la valeur propre du milieu est sans
importance. Le code vérifie le théorème 11.4, le corollaire 11.5 et les configurations extrémales de la
proposition 11.6 pour $k=3$, $C=2$.

```python
M5 = TypeIAlgebra([5], [5], precision='double')
A = M5.operator(torch.randn(2000, 1, 5, 5, dtype=torch.complex128))         # Ginibre, no structure
B = M5.operator(torch.randn(2000, 1, 5, 5, dtype=torch.complex128))
sA, sB = torch.linalg.svdvals(A.matrix[:, 0]), torch.linalg.svdvals(B.matrix[:, 0])   # decreasing
assert ((A @ B).trace.abs() <= (sA * sB).sum(-1) + 1e-12).all()             # Theorem 11.4
Ph, Qh = (M5.operator_from_eigenvalues(lambda k: torch.randn(2000, k), batch_size=2000, force_self_adjoint=True)
          for _ in range(2))
p, q = Ph.eigenvalues()[0], Qh.eigenvalues()[0]                             # increasing order
t = (Ph @ Qh).trace.real
assert (t <= (p * q).sum(-1) + 1e-12).all() and (t >= (p * q.flip(-1)).sum(-1) - 1e-12).all()   # Cor. 11.5(a)
Pp, Qp = Ph @ Ph, Qh @ Qh                                                   # positive, spectra p^2, q^2
nuc = (Pp @ Qp).trace_norm()
assert ((Pp @ Qp).trace.real <= nuc + 1e-12).all()                          # Cor. 11.5(b)
assert (nuc <= ((p ** 2).sort(-1)[0] * (q ** 2).sort(-1)[0]).sum(-1) + 1e-10).all()

alg3 = TypeIAlgebra([3, 3], [3, 3], precision='double')
x_dn = torch.tensor([1.0, 0.6, ell(0.5)], dtype=torch.float64)              # Delta = 0.5, middle value 0.6
y_dn = torch.tensor([1.0, 0.7, ell(0.8)], dtype=torch.float64)              # Delta = 0.8, middle value 0.7
J = torch.eye(3, dtype=torch.float64).flip(0)                               # reversal permutation
blocks = lambda m: alg3.from_blocks([m, m])                                 # the same block in both sectors

def z_exp2(X0, V, Y0):
    return ((X0 @ V @ Y0).trace_norm() - (X0 @ Y0).trace.real).item()

X0 = blocks(torch.diag(x_dn))
z_max = z_exp2(X0, blocks(J), blocks(torch.diag(y_dn.flip(0))))
z_min = z_exp2(X0, blocks(J), blocks(torch.diag(y_dn)))
exact = 2 * (x_dn * (y_dn - y_dn.flip(0))).sum().item()                     # Proposition 11.6
envelope = 2 * (1 - ell(0.5)) * (1 - ell(0.8))                              # Corollary 11.7, floor(3/2) = 1
print(f"Exp. 2, k = 3, C = 2: z_max = {z_max:.6f}, z_min = {z_min:.6f}, formula {envelope:.6f}")
assert abs(z_max - exact) < 1e-12 and abs(z_min + exact) < 1e-12 and abs(exact - envelope) < 1e-12
```

## 11.3 Une inégalité de traces avec cas d'égalité

Pour $X=X^*$, soit $X_\pm=f_\pm(X)$ avec $f_+(t)=\max(t,0)$, $f_-(t)=\max(-t,0)$. Alors $X_\pm\ge0$, $X=X_+-X_-$,
$|X|=X_++X_-$ et $X_+X_-=0$ (théorème 1.5(c)).

**Lemme 11.8.** Pour $A,B\ge0$ dans $M_n(\mathbb C)$, $\|A-B\|_1\le\operatorname{Tr}A+\operatorname{Tr}B$, avec égalité si et
seulement si $AB=0$.

*Démonstration.* Soient $H=A-B$ et $Q=g(H)$ avec $g=1$ sur $(0,\infty)$ et $g=0$ ailleurs. Comme $|t|=t(2g(t)-1)$,
on a $|H|=H(2Q-1)$ et
$$
\|H\|_1=\operatorname{Tr}\big(A(2Q-1)\big)-\operatorname{Tr}\big(B(2Q-1)\big)
=\operatorname{Tr}A+\operatorname{Tr}B-2\operatorname{Tr}\big((1-Q)A(1-Q)\big)-2\operatorname{Tr}(QBQ).
$$
Les traces soustraites sont $\ge0$, et $\operatorname{Tr}(RAR)=\|A^{1/2}R\|_2^2$ s'annule si et seulement si $AR=0$ ;
l'égalité a donc lieu si et seulement si $A(1-Q)=0$ et $BQ=0$, et alors $AB=AQB=A(BQ)^*=0$. Réciproquement, si
$AB=0$, alors $BA=0$, $(A-B)^2=A^2+B^2=(A+B)^2$, et le théorème 1.6 donne $|A-B|=A+B$. $\square$

**Proposition 11.9 (Exp. 3).** Pour $X=X^*$ et $Y\ge0$ dans $\mathcal M$,
$\operatorname{Tr}|YXY|\le\operatorname{Tr}(Y|X|Y)$, avec égalité si et seulement si $YX_+Y^2X_-Y=0$, en particulier si
$XY=YX$.

*Démonstration.* Dans $M_D(\mathbb C)$, $A=YX_+Y$ et $B=YX_-Y$ sont positifs (corollaire 1.7 avec $T=X_\pm^{1/2}Y$),
$YXY=A-B$ et $\operatorname{Tr}(Y|X|Y)=\operatorname{Tr}A+\operatorname{Tr}B$ ; appliquons le lemme 11.8. Si $XY=YX$, alors
$Y$ commute avec $X_\pm$ (théorème 1.5(d)), et $AB=Y^4X_+X_-=0$. $\square$

Ainsi $z\ge0$ dans Exp. 3, et $\underline z=0$ (prendre $Y=1$). La condition d'égalité est plus faible que la
commutation : pour $X=\operatorname{diag}(1,-1,0)$ et $Y$ de lignes $(1,0,a)$, $(0,1,0)$, $(a,0,1)$, $|a|<1$, on a
$AB=(Ye_1)\langle Ye_1,Ye_2\rangle(Ye_2)^*=0$, alors que $XY\ne YX$ pour $a\ne0$.

```python
def z_exp3(X, Y):
    return (Y @ X.abs() @ Y).trace.real - (Y @ X @ Y).trace_norm()

Xs = alg3.operator_from_eigenvalues(lambda k: torch.randn(500, k), batch_size=500, force_self_adjoint=True)
Yp = alg3.operator_from_eigenvalues(lambda k: torch.rand(500, k), batch_size=500, force_positive=True)
assert (z_exp3(Xs, Yp) >= -1e-12).all()                                     # Proposition 11.9

M3 = TypeIAlgebra([3], [3], precision='double')
Xe = M3.operator(torch.diag(torch.tensor([1.0, -1.0, 0.0]))[None, None])
Ye = M3.operator(torch.tensor([[1.0, 0, 0.5], [0, 1, 0], [0.5, 0, 1]])[None, None])
print(f"equality without commutation: z = {z_exp3(Xe, Ye).item():.1e}, "
      f"||XY - YX|| = {(Xe @ Ye - Ye @ Xe).operator_norm().item():.2f}")
assert abs(z_exp3(Xe, Ye).item()) < 1e-12 and (Xe @ Ye - Ye @ Xe).operator_norm().item() > 0.1
```

## 11.4 Échantillonnage contre optimisation

Quelle est la taille typique de $|\operatorname{Tr}(AU)|$, comparée à son maximum $\|A\|_1$ ? La bibliothèque tire $U$
comme le facteur $Q$ de la décomposition QR $Z=QR$, $R_{ii}>0$, d'une **matrice de Ginibre** $Z$ (coefficients
indépendants de densité $\pi^{-1}e^{-|z|^2}$ ; recette de Mezzadri), indépendamment dans chaque bloc.

**Lemme 11.10.** $Q$ est distribué selon la loi de Haar sur $U(k)$ ; en particulier, pour tout unitaire fixé
$V\in M_k(\mathbb C)$, $VQ$ a la même loi que $Q$. L'unitaire $\bigoplus_cQ_c$ à blocs indépendants est distribué
selon la loi de Haar sur le groupe unitaire de $\mathcal M$.

*Démonstration.* C'est le théorème 5.3 (l'invariance à gauche est l'étape d'équivariance de sa démonstration) et le
corollaire 5.4. $\square$

(Pour la concentration de la mesure et davantage sur ces groupes, voir E. Meckes, *The Random Matrix Theory of the
Classical Compact Groups*, Cambridge 2019.)

**Proposition 11.11.** Pour $U\in U(k)$ de loi de Haar et $A\in M_k(\mathbb C)$ : $\mathbb E\operatorname{Tr}(AU)=0$ et
$\mathbb E|\operatorname{Tr}(AU)|^2=\|A\|_2^2/k$.

*Démonstration.* $\mathbb EU_{ij}=0$ d'après le corollaire 5.5. La proposition 5.7 (avec $C=1$) pour
$X=|e_j\rangle\langle e_m|$ donne la formule de Weingarten d'ordre un,
$\mathbb E[U_{ij}\bar U_{lm}]=\mathbb E(UXU^*)_{il}=\delta_{il}\delta_{jm}/k$ ; pour $k\ge2$, elle découle aussi de la
formule d'ordre deux, le théorème 5.10, en posant $i_2=i'_2$, $j_2=j'_2$ et en sommant sur les deux indices, puisque
$\sum_{i_2,j_2}|U_{i_2j_2}|^2=k$. Avec $\operatorname{Tr}(AU)=\sum_{ij}A_{ji}U_{ij}$,
$\mathbb E|\operatorname{Tr}(AU)|^2=\sum_{ijlm}A_{ji}\bar A_{ml}\mathbb E[U_{ij}\bar U_{lm}]=\frac1k\sum_{ij}|A_{ji}|^2$. $\square$

**Corollaire 11.12.** Soient $X,Y\ge0$ dans $M_k(\mathbb C)$ à spectres dans $[a,b]$, $a>0$, et $U$ de loi de Haar.
Alors
$$
\frac{\big(\mathbb E|\operatorname{Tr}(XUY)|^2\big)^{1/2}}{\max_V|\operatorname{Tr}(XVY)|}=\frac{\|YX\|_2}{\sqrt k\,\|YX\|_1}\le\min\Big(\frac1{\sqrt k},\frac{b^2}{a^2k}\Big),
\qquad
\mathbb P\big(|\operatorname{Tr}(XUY)|\ge\operatorname{Tr}(XY)\big)\le\frac{b^4}{a^4k^2}.
$$

*Démonstration.* L'égalité résulte de la proposition 11.11 pour $A=YX$ et de la proposition 3.9.
$\sum_is_i^2\le(\sum_is_i)^2$ donne $1/\sqrt k$. Comme $s_1(A)\le\|Y\|\|X\|\le b^2$ et
$s_k(A)=\min_{\|v\|=1}\|YXv\|\ge a^2$, on a $\|A\|_2^2\le kb^4$ et $\|A\|_1\ge ka^2$. Enfin
$\operatorname{Tr}(XY)=\operatorname{Tr}(X^{1/2}YX^{1/2})\ge a\operatorname{Tr}X\ge ka^2$, et l'inégalité de Markov pour
$|\operatorname{Tr}(XUY)|^2$ majore la probabilité par $(kb^4/k)/(ka^2)^2$. $\square$

Ainsi, $N$ échantillons de Haar trouvent une violation de $|\operatorname{Tr}(XUY)|\le\operatorname{Tr}(XY)$ avec une
probabilité d'au plus $Nb^4/(a^4k^2)$, bien qu'il en existe une pour tout couple qui ne commute pas. La réalité est
plus radicale encore : $f(U)=|\operatorname{Tr}(AU)|$ est lipschitzienne de constante $L=\|A\|_2$ pour $\|U-V\|_2$, donc
$\mathbb P(f\ge\mathbb Ef+t)\le e^{-ckt^2/L^2}$ avec une constante absolue $c>0$ (concentration de la mesure, énoncée
sans démonstration ; Meckes, *loc. cit.*, chap. 5) : ici, un exposant de l'ordre de $k^2a^4/b^4$. Le code vérifie la
proposition 11.11 et affiche, pour des spectres uniformes dans $[0.2,1.2]$ ($b/a=6$), le rapport du
corollaire 11.12, $\sup z$ et le meilleur $z$ parmi $2000$ échantillons, tous deux rapportés à $\operatorname{Tr}(XY)$.

```python
def random_pair(alg_k):
    pos = lambda k: 0.2 + torch.rand(1, k)                                  # spectrum in [a, b] = [0.2, 1.2]
    return tuple(alg_k.operator_from_eigenvalues(pos, force_positive=True) for _ in range(2))

M6 = TypeIAlgebra([6], [6], precision='double')
Xf, Yf = random_pair(M6)
N = 20_000
t2 = (Yf @ Xf @ M6.random_unitary_operator(N)).trace.abs() ** 2
m2, err, exact2 = t2.mean().item(), t2.std().item() / math.sqrt(N), (Yf @ Xf).frobenius_norm().item() ** 2 / 6
print(f"E|Tr(AU)|^2 = {m2:.4f} +- {err:.4f}   (exact ||A||_2^2 / k = {exact2:.4f})")
assert abs(m2 - exact2) < 5 * err                                           # Proposition 11.11

print(" k   rms/max   sup z/TrXY   best of 2000 samples")
for k in (2, 4, 8, 16, 32):
    Mk = TypeIAlgebra([k], [k], precision='double')
    Xk, Yk = random_pair(Mk)
    trXY, nuc = (Xk @ Yk).trace.real.item(), (Yk @ Xk).trace_norm().item()
    ratio = (Yk @ Xk).frobenius_norm().item() / (math.sqrt(k) * nuc)
    best = (Xk @ Mk.random_unitary_operator(2000) @ Yk).trace.abs().max().item()
    print(f"{k:2d}   {ratio:.4f}    {(nuc - trXY) / trXY:.2e}     {(best - trXY) / trXY:+.3f}")
    assert ratio <= min(1 / math.sqrt(k), 36 / k) + 1e-12                   # Corollary 11.12
    assert nuc > trXY and best < nuc and (k < 4 or best < trXY)            # from k = 4 on: no violation found
```

Le rapport décroît comme $1/k$, et à partir de $k=4$ le meilleur échantillon n'atteint même pas $z=0$. Les extrêmes
sont des ensembles de mesure de Haar minuscule, qu'une méthode de gradient atteint directement. Pour le couple
`Xf, Yf`, les deux extrêmes en $U$ sont connus :

**Lemme 11.13.** Si $k\ge2$, alors $\min_U|\operatorname{Tr}(AU)|=0$ pour toute $A\in M_k(\mathbb C)$.

*Démonstration.* Soient $A=W\Sigma V^*$ et $P$ le décalage cyclique, $Pe_i=e_{(i\bmod k)+1}$, dont la diagonale est
nulle pour $k\ge2$. Pour $U=VPW^*$, $\operatorname{Tr}(AU)=\operatorname{Tr}(\Sigma P)=0$. $\square$

Ainsi $\inf_Uz=-\operatorname{Tr}(XY)$ et $\sup_Uz=\|YX\|_1-\operatorname{Tr}(XY)$. Un seul appel à `extremize` trouve les
deux, avec $12$ départs qui maximisent $z$ et $12$ qui maximisent $-z$. L'échantillonnage s'approche de l'infimum,
qui est voisin des valeurs typiques ($|\operatorname{Tr}(XUY)|$ est typiquement de l'ordre de $\|YX\|_2/\sqrt k$), mais
pas du supremum :

```python
trXY = (Xf @ Yf).trace.real
sup_exact, inf_exact = ((Yf @ Xf).trace_norm() - trXY).item(), -trXY.item()
U6 = opt.UnitaryParam(M6, 24)                                               # 24 starts
pm = torch.tensor([1.0, -1.0], dtype=torch.float64).repeat_interleave(12)   # 12 maximise z, 12 maximise -z
best = opt.extremize(lambda: pm * ((Xf @ U6.operator() @ Yf).trace.abs() - trXY), [U6],
                     maximize=True, rounds=6, steps=30).reshape(2, 12).max(1)[0]
best_up, best_down = best[0].item(), -best[1].item()
z_mc = (Xf @ M6.random_unitary_operator(N) @ Yf).trace.abs() - trXY
print(f"sup z: exact {sup_exact:.6e}, ascent {best_up:.6e}, best sample {z_mc.max().item():+.4f}")
print(f"inf z: exact {inf_exact:.6f}, descent {best_down:.6f}, worst sample {z_mc.min().item():+.6f}")
assert best_up <= sup_exact + 1e-12 and best_down >= inf_exact - 1e-12    # attained values: one-sided
assert sup_exact - best_up < 1e-3 * sup_exact and best_down - inf_exact < 1e-3 * abs(inf_exact)
assert z_mc.min().item() - inf_exact < 1e-2 * abs(inf_exact)               # sampling finds the typical side ...
assert z_mc.max().item() < 0                                               # ... but misses the violation
```

> **Interprétation physique.** La concentration de la mesure est la forme mathématique de la **typicité** en
> mécanique statistique : en grande dimension, presque tout état ou toute évolution donne la même valeur d'une
> grandeur macroscopique. Une inégalité optimale se joue sur des configurations finement ajustées — bases propres
> alignées, spectres à deux niveaux — qu'un unitaire aléatoire (un circuit aléatoire, une matrice de diffusion
> aléatoire) ne produit pour ainsi dire jamais.

## 11.5 Recherche multi-départs des enveloppes

`optimize.extremize(objective, params, maximize)` exécute Adam en `rounds` tours de `steps` pas (taux
d'apprentissage `lr * decay**round`), recentre les cartes unitaires après chaque tour et renvoie, **pour chaque
départ**, la meilleure valeur atteinte. Chaque élément du lot est un départ indépendant avec ses propres
emplacements, signes, $U_0$ et contraste. Les signes de `SelfAdjointParam` sont discrets : entre les tours,
`discrete_step` essaie huit changements de signe individuels par départ, l'un après l'autre, et conserve chaque
changement qui améliore l'objectif.

Comme chaque itéré est admissible (proposition 11.2), chaque valeur renvoyée est atteinte : une maximisation donne
une **borne inférieure** du supremum, une minimisation une borne supérieure de l'infimum, et rien de plus — une
recherche locale peut s'arrêter en un extremum local. L'autre sens exige une démonstration, comme le
corollaire 11.7 et la proposition 11.9, ou une relaxation convexe (section 11.6). Une enveloppe calculée est
confrontée (i) à des formes closes ; (ii) à des configurations construites à la main, telles que des sommes directes
de solutions plus petites, en dessous desquelles elle ne doit pas tomber ; (iii) à des échantillons structurés
(spectres à deux niveaux, bases de permutation) qui s'approchent des extrêmes en toute dimension, contrairement aux
échantillons de Haar (`scripts/bounds_checks.py`).

**Proposition 11.14 (enveloppe inférieure d'Exp. 1).** Supposons tous les $k_c\ge2$. Sur $X\in\mathcal P_\delta$,
$Y\ge0$ avec $\lambda_{\max}(Y)=1$ et les unitaires $U\in\mathcal M$, $\ \min z=-\big(D-2\delta/(1+\delta)\big)$.

*Démonstration.* $Y\le1$ donne $X^{1/2}YX^{1/2}\le X$, donc $\operatorname{Tr}(XY)\le\operatorname{Tr}X\le(D-1)+\ell(\delta)$,
puisqu'une valeur propre de $X$ vaut $\ell(\delta)$ et que les autres sont $\le1$ ; d'où
$z\ge-\operatorname{Tr}(XY)\ge-(D-2\delta/(1+\delta))$. L'égalité a lieu pour $X$ de valeurs propres $\ell(\delta)$ une
fois et $1$ sinon, $Y=1$, et $U=\bigoplus_cU_c$ avec $\operatorname{Tr}(X_cU_c)=0$ (lemme 11.13). $\square$

Le code vérifie la proposition 11.14 sur $M_2\oplus M_2$, le corollaire 11.7 sur $M_3\oplus M_3$ (les deux contrastes
prescrits, $k=3$ impair) et la proposition 11.9 ; suprema et infima partagent un même appel, comme ci-dessus. Le
mot-clé `is_self_adjoint` choisit l'algorithme de la norme trace (valeurs singulières, ou valeurs propres pour
$YXY$) sans test à chaque pas. La dernière vérification est de type (ii) : le $z$ d'Exp. 3 est additif sur les
sommes directes, et le contraste d'une somme directe est le plus grand contraste de ses termes (exercice 3), donc
$$
\overline z_{M_2\oplus M_2}(\delta)\ \ge\ \overline z_{M_2}(\delta)+\max_{\delta'\le\delta}\overline z_{M_2}(\delta').
$$

```python
torch.manual_seed(1)
S = 6                                                                       # starts per grid point
grid = torch.tensor([0.0, 0.5, 0.9], dtype=torch.float64)
alg22 = TypeIAlgebra([2, 2], [2, 2], precision='double')                    # D = 4

# Exp. 1, lower envelope along Delta(X); Delta(Y) free in [0, 0.999]
X1, Y1, U1 = (opt.PositiveParam(alg22, 3 * S, delta=grid.repeat_interleave(S)),
              opt.PositiveParam(alg22, 3 * S), opt.UnitaryParam(alg22, 3 * S))
def z1():
    x, y = X1.operator(), Y1.operator()
    return (x @ U1.operator() @ y).trace.abs() - (x @ y).trace.real
inf1 = opt.extremize(z1, [X1, Y1, U1], maximize=False, rounds=4, steps=30).reshape(3, S).min(1)[0]
exact1 = -(4 - 2 * grid / (1 + grid))                                       # Proposition 11.14
print("Exp. 1  inf:", [f"{v:.4f}" for v in inf1.tolist()], " exact:", [f"{v:.4f}" for v in exact1.tolist()])
assert (inf1 >= exact1 - 1e-12).all() and (inf1 - exact1).max() < 4e-3

# Exp. 2 on M_3 (+) M_3, Delta(X0) = 0.5, Delta(Y0) = 0.8: sup and inf in one call
sign = torch.tensor([1.0, -1.0], dtype=torch.float64).repeat_interleave(S)
X2, Y2, V2 = (opt.PositiveParam(alg3, 2 * S, delta=torch.full((2 * S,), 0.5, dtype=torch.float64)),
              opt.PositiveParam(alg3, 2 * S, delta=torch.full((2 * S,), 0.8, dtype=torch.float64)),
              opt.UnitaryParam(alg3, 2 * S))
def z2():
    x0, y0 = X2.operator(), Y2.operator()
    M = alg3.operator((x0 @ V2.operator() @ y0).matrix, is_self_adjoint=False)
    return sign * (M.trace_norm() - (x0 @ y0).trace.real)
best2 = opt.extremize(z2, [X2, Y2, V2], maximize=True, rounds=4, steps=30).reshape(2, S).max(1)[0]
sup2, inf2 = best2[0].item(), -best2[1].item()
print(f"Exp. 2  sup {sup2:.6f}, inf {inf2:.6f}, exact +-{envelope:.6f}")
assert sup2 <= envelope + 1e-12 and inf2 >= -envelope - 1e-12               # Corollary 11.7: one-sided ...
assert envelope - sup2 < 1e-3 * envelope and inf2 + envelope < 1e-3 * envelope   # ... and reached

# Exp. 3: maximise sign * z at the contrasts `deltas` of |X|; Delta(Y) free in [0, 0.999]
def exp3(alg_, deltas, sign):
    X3, Y3 = opt.SelfAdjointParam(alg_, len(deltas), delta=deltas), opt.PositiveParam(alg_, len(deltas))
    def z3():
        y = Y3.operator()
        M = alg_.operator((y @ X3.operator() @ y).matrix, is_self_adjoint=True)
        return sign * ((y @ X3.abs_operator() @ y).trace.real - M.trace_norm())
    return opt.extremize(z3, [X3, Y3], maximize=True, rounds=4, steps=30).reshape(2, S).max(1)[0]

best3 = exp3(alg22, torch.full((2 * S,), 0.5, dtype=torch.float64), sign)
sup_M22, inf3 = best3[0].item(), -best3[1].item()
s_M2 = exp3(TypeIAlgebra([2], [2], precision='double'),
            torch.tensor([0.0, 0.5], dtype=torch.float64).repeat_interleave(S), 1.0)
direct_sum = (s_M2[1] + s_M2.max()).item()                                 # delta = 0.5, best delta' <= 0.5
conj = [(1 - dl) * (2 * 0.999 / 1.999) ** 2 for dl in (0.0, 0.5)]          # unproven closed form, d = 0.999
print(f"Exp. 3  inf {inf3:.1e};  sup for M_2 at delta = 0, 0.5: {s_M2[0]:.5f}, {s_M2[1]:.5f}"
      f"  (closed form {conj[0]:.5f}, {conj[1]:.5f})")
print(f"        sup for M_2 (+) M_2 at 0.5: search {sup_M22:.5f}, direct sum {direct_sum:.5f}")
assert inf3 >= -1e-12 and inf3 < 1e-6                                       # Proposition 11.9
assert sup_M22 >= direct_sum - 5e-3 * direct_sum                           # (ii), up to 0.5 %
```

La recherche sur $M_2\oplus M_2$ reste en dessous de la somme directe d'environ $0.2\,\%$ : elle n'a pas entièrement
convergé, et la meilleure borne inférieure de $\overline z_{M_2\oplus M_2}(0.5)$ est le plus grand des deux nombres —
détecter cela est le but de la vérification (ii). Les valeurs pour $M_2$ coïncident à $10^{-4}$ près avec
$(1-\delta)\big(2d/(1+d)\big)^2=(1-\delta)(1-\ell(d))^2$, où $d=0.999$ est le plus grand contraste admis pour $Y$ :
une forme close que l'article qui accompagne la bibliothèque a observée numériquement (dans des algèbres réelles)
sans pouvoir la démontrer ; le code ne l'affirme pas. Pour des blocs plus grands, les enveloppes supérieures de
Exp. 1 et Exp. 3 ne sont pas connues sous forme close, à notre connaissance ; les courbes calculées sont alors des
bornes inférieures et des conjectures, rien de plus.

## 11.6 Points de comparaison

**Manopt**, **pymanopt** et **geoopt** fournissent des variétés matricielles (variétés de Stiefel, groupes orthogonaux
et unitaires) avec des rétractions telles que QR, polaire ou de Cayley, et des solveurs riemanniens (gradients
conjugués, régions de confiance, Adam riemannien) souvent plus efficaces qu'Adam dans une carte recentrée. **cvxpy**
avec un solveur SDP résout des problèmes convexes en variables opérateurs jusqu'à l'optimalité globale, avec un
certificat dual à la tolérance du solveur près, et les relaxations semi-définies donnent des bornes rigoureuses dans
le sens inaccessible à une recherche locale. **NCTSSOS**, **ncpol2sdpa** et **NCSOStools** bornent les valeurs propres
et les traces de polynômes non commutatifs en toutes dimensions à la fois (section 0.3). La bibliothèque apporte des
paramétrages exacts des contraintes spectrales et de contraste de $\bigoplus_cM_{k_c}$ (les extrêmes peuvent se
trouver dans n'importe quel secteur), des unitaires qui respectent les blocs, et un seul calcul par lots pour de
nombreux départs et valeurs de contrainte, sur CPU ou GPU. Elle produit des candidats extrémaux et des bornes
inférieures pour les suprema ; démonstrations ou relaxations doivent fournir le reste.

## Exercices

1. (**Transformation de Cayley.**) Pour $S^*=-S$, montrer que $C(S)=(1-S)(1+S)^{-1}$ est unitaire et que
   $S\mapsto C(S)$ est une bijection sur les unitaires sans valeur propre $-1$. Ainsi, contrairement à l'application
   exponentielle (lemme 11.1), une carte de Cayley $U_0C(S)$ manque $-U_0$.
2. (**Mesure de Haar réelle.**) Pour $O$ distribué selon la loi de Haar sur $O(k)$ (QR d'une matrice de Ginibre
   réelle) et $A$ réelle, montrer que $\mathbb E\operatorname{Tr}(AO)^2=\|A\|_2^2/k$, et le vérifier avec une algèbre
   réelle.
3. (**Sommes directes.**) Montrer que le $z$ d'Exp. 3 est additif sur les sommes directes et que
   $\Delta(A\oplus A')=\max(\Delta(A),\Delta(A'))$ pour $A,A'\ge0$ avec $\lambda_{\max}=1$. En déduire l'inégalité de
   somme directe de la section 11.5 (avec $\Delta(Y)$ libre dans un intervalle fixé $[0,b]$).

<details>
<summary>Solutions</summary>

1. $S$ a des valeurs propres imaginaires pures, donc $1\pm S$ sont inversibles ; ils commutent, et $(1+S)^*=1-S$.
   Par conséquent $C^*C=(1-S)^{-1}(1+S)(1-S)(1+S)^{-1}=1$, et $1+C=2(1+S)^{-1}$ est inversible, donc
   $-1\notin\operatorname{spec}C$. Réciproquement, pour un unitaire $U$ avec $1+U$ inversible, $S=(1-U)(1+U)^{-1}$ est
   antiautoadjoint : $S^*=(1+U^*)^{-1}(1-U^*)=(U+1)^{-1}U\,U^*(U-1)=-S$. Comme $C$ commute avec $S$, $C(1+S)=1-S$
   équivaut à $S(1+C)=1-C$, c'est-à-dire $S=(1-C)(1+C)^{-1}$ : les deux applications sont réciproques l'une de
   l'autre.

```python
Sk = torch.randn(4, 4, dtype=torch.complex128)
Sk = Sk - Sk.mH                                                   # skew-adjoint
I4 = torch.eye(4, dtype=torch.complex128)
Cay = (I4 - Sk) @ torch.linalg.inv(I4 + Sk)
assert torch.allclose(Cay.mH @ Cay, I4, atol=1e-12)
assert torch.allclose((I4 - Cay) @ torch.linalg.inv(I4 + Cay), Sk, atol=1e-10)
assert torch.linalg.eigvals(Cay).add(1).abs().min() > 1e-6        # -1 is never an eigenvalue
```

2. Le théorème 5.3 et sa démonstration valent pour une matrice gaussienne réelle (densité
   $\propto e^{-\|Z\|_2^2/2}$), le groupe orthogonal et $V$ orthogonale (la remarque qui suit le théorème 5.3), donc
   $VO$ a la même loi que $O$. La multiplication à gauche par $\operatorname{diag}(\pm1)$ donne
   $\mathbb E[O_{ij}O_{lm}]=0$ pour $i\ne l$, les permutations montrent que $\mathbb E[O_{ij}O_{im}]$ ne dépend pas de
   $i$, et $\sum_iO_{ij}O_{im}=\delta_{jm}$ ; donc $\mathbb E[O_{ij}O_{lm}]=\delta_{il}\delta_{jm}/k$ et
   $\mathbb E\operatorname{Tr}(AO)^2=\sum_{ij}A_{ji}^2/k$.

```python
Mr = TypeIAlgebra([5], [5], complex_valued=False, precision='double')
Ar = Mr.operator(torch.randn(1, 1, 5, 5, dtype=torch.float64))
tr2 = (Ar @ Mr.random_unitary_operator(20_000)).trace ** 2
assert abs(tr2.mean().item() - Ar.frobenius_norm().item() ** 2 / 5) < 5 * tr2.std().item() / math.sqrt(20_000)
```

3. Les traces et les normes trace d'opérateurs diagonaux par blocs sont des sommes sur les blocs, et
   $|X\oplus X'|=|X|\oplus|X'|$ ; $\lambda_{\min}(A\oplus A')=\min(\lambda_{\min}(A),\lambda_{\min}(A'))$, et $\Delta$
   décroît avec $\lambda_{\min}$. Si donc $(X,Y)$ est admissible en $\delta$ et $(X',Y')$ en $\delta'\le\delta$, alors
   $(X\oplus X',Y\oplus Y')$ est admissible en $\delta$, avec $z=z(X,Y)+z(X',Y')$ ; on passe au supremum sur les deux
   couples et sur $\delta'\le\delta$.

</details>

**Suite :** le [chapitre 12](12_equivariant_learning.md) traite de l'apprentissage automatique équivariant : les couches qui
commutent avec une symétrie sont des éléments d'un commutant (chapitre 2), construits à partir des mêmes opérateurs
différentiables et traités par lots.
