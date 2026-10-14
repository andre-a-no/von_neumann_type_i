[English](../03_traces_and_states.md) · [Español](../es/03_traces_and_states.md) · **Français** · [Deutsch](../de/03_traces_and_states.md) · [中文](../zh/03_traces_and_states.md) · [日本語](../ja/03_traces_and_states.md) · [Русский](../ru/03_traces_and_states.md)

# 3. Traces, états et entropie

Sur $M_n(\mathbb C)$, il n'existe pour l'essentiel qu'une seule trace. Sur
$\bigoplus_cM_{k_c}$, on dispose d'une famille à $C$ paramètres, et le choix importe :
densités, entropies et applications duales en dépendent toutes. Ce chapitre classifie les
traces, représente les états par des densités, démontre les inégalités d'entropie
fondamentales, et se termine par une inégalité de traces plausible qui est fausse.

**Vous apprendrez** les trois traces de la bibliothèque et comment passer de l'une à
l'autre, la partie classique d'un état, la décomposition de l'entropie de von Neumann
selon les secteurs, l'inégalité de Klein et ses conséquences, ainsi que le contraste de
Michelson.

Dans tout ce chapitre, $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ avec tous les
$k_c\ge1$, $D=\sum_ck_c$, les $P_c$ sont les projecteurs centraux minimaux, et
$A=\bigoplus_cA_c$.

## 3.1 Traces

Une **trace** sur $\mathcal M$ est une forme linéaire $\tau$ telle que $\tau(AB)=\tau(BA)$
pour tous $A,B$. Elle est **positive** si $\tau(A^*A)\ge0$ pour tout $A$, et **fidèle** si,
de plus, $\tau(A^*A)=0$ entraîne $A=0$.

**Théorème 3.1.** Toute trace sur $\mathcal M$ est de la forme
$\operatorname{tr}_w(A)=\sum_cw_c\operatorname{Tr}A_c$ avec $w\in\mathbb C^C$. Elle est positive
si et seulement si tous les $w_c\ge0$, et fidèle si et seulement si tous les $w_c>0$.

*Démonstration.* Soient $E^c_{ij}$ les unités matricielles du $c$-ième bloc. Pour une trace
$\tau$,
$\tau(E^c_{ij})=\tau(E^c_{i1}E^c_{1j})=\tau(E^c_{1j}E^c_{i1})=\delta_{ij}\tau(E^c_{11})$.
Posons $w_c=\tau(E^c_{11})$ ; par linéarité, $\tau(A)=\sum_c\sum_i(A_c)_{ii}w_c$.
Réciproquement, toute $\operatorname{tr}_w$ est une trace. Positivité :
$E^c_{11}=(E^c_{11})^*E^c_{11}$, donc la condition $w_c=\tau\big((E^c_{11})^*E^c_{11}\big)\ge0$ est
nécessaire, et $w_c>0$ est nécessaire pour la fidélité, puisque $E^c_{11}\neq0$. Réciproquement,
$\operatorname{tr}_w(A^*A)=\sum_cw_c\|A_c\|_2^2$, qui est $\ge0$ si tous les $w_c\ge0$, et ne
s'annule que pour $A=0$ si tous les $w_c>0$. $\square$

La bibliothèque implémente trois traces fidèles :

| nom | poids $w_c$ | valeur en $1$ | signification |
|---|---|---|---|
| `Tr_blunt` | $1$ | $D$ | trace de $\bigoplus_c\mathbb C^{k_c}$ (multiplicités $m_c=1$) |
| `Tr_norm` | $1/k_c$ | $C$ | chaque bloc normalisé |
| `tau_vN` | $1/(Ck_c)$ | $1$ | l'état tracial de poids égaux $1/C$ sur les secteurs |

D'après la section 2.3, la trace d'un espace sur lequel $\mathcal M$ agit avec les
multiplicités $m_c$ est $\operatorname{tr}_w$ avec $w_c=m_c$. Choisir une trace fidèle, c'est
donc choisir « combien compte chaque secteur ».

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, DensityMatrix, states, channels

torch.manual_seed(0)
alg = TypeIAlgebra([2, 3, 5], [2, 3, 5], precision='double')    # M_2 (+) M_3 (+) M_5, D = 10, C = 3
U = alg.random_unitary_operator(4)
A = alg.operator_from_eigenvalues(lambda d: torch.randn(4, d), batch_size=4, force_self_adjoint=True)
B = U @ A                                                       # a batch of four non-normal operators
one = alg.identity()

for f, value in ((lambda X: X.Tr_blunt(), 10), (lambda X: X.Tr_norm(), 3), (lambda X: X.tau_vN(), 1)):
    assert abs(f(one).item() - value) < 1e-12                   # normalisation
    assert (f(A @ B) - f(B @ A)).abs().max() < 1e-10            # trace property
assert (A.Tr_norm() - 3 * A.tau_vN()).abs().max() < 1e-12       # Tr_norm = C tau_vN
```

## 3.2 États et densités

Un **état** est une forme linéaire $\omega$ sur $\mathcal M$ positive
($\omega(A^*A)\ge0$) et normalisée ($\omega(1)=1$).

**Théorème 3.2.** Soit $\operatorname{tr}_w$ une trace fidèle. Pour toute forme linéaire
$\omega$, il existe un unique $\rho_w\in\mathcal M$ tel que
$\omega(A)=\operatorname{tr}_w(\rho_wA)$ pour tout $A$. La forme $\omega$ est positive si et
seulement si $\rho_w\ge0$, et c'est un état si et seulement si, de plus,
$\operatorname{tr}_w(\rho_w)=1$. Pour deux traces fidèles, les densités sont reliées bloc par
bloc par
$$
w_c\,(\rho_w)_c=w'_c\,(\rho_{w'})_c .
$$

*Démonstration.* $\langle X,Y\rangle_w=\operatorname{tr}_w(X^*Y)$ est un produit scalaire sur
$\mathcal M$ (théorème 3.1). D'après le théorème de représentation de Riesz dans l'espace
de dimension finie $\mathcal M$, il existe un unique $R$ tel que
$\omega(A)=\langle R,A\rangle_w=\operatorname{tr}_w(R^*A)$ ; on pose $\rho_w=R^*$. Si
$\rho_w\ge0$, alors $\omega(A^*A)=\sum_cw_c\operatorname{Tr}(A_c\,(\rho_w)_c\,A_c^*)\ge0$,
puisque $A_c(\rho_w)_cA_c^*\ge0$ (corollaire 1.7 avec $T=(\rho_w)_c^{1/2}A_c^*$), qu'une
matrice positive a une trace $\ge0$ (ses valeurs propres sont $\ge0$, corollaire 1.4(3)) et que
$w_c>0$. Réciproquement, si $\omega$ est positive, prenons $x\in\mathbb C^{k_c}$ et pour $A$
l'opérateur de rang un $|x\rangle\langle x|$ dans le bloc $c$ ($x\neq0$) ; alors $A=A^*A/\|x\|^2$ et
$0\le\omega(A)=w_c\langle x,(\rho_w)_cx\rangle$, donc $(\rho_w)_c\ge0$. La normalisation
s'écrit $\omega(1)=\operatorname{tr}_w(\rho_w)$. Enfin,
$\sum_cw_c\operatorname{Tr}((\rho_w)_cA_c)$ et $\sum_cw'_c\operatorname{Tr}((\rho_{w'})_cA_c)$
sont tous deux égaux à $\omega(A)$ pour tout $A$, donc $\bigoplus_cw_c(\rho_w)_c$ et
$\bigoplus_cw'_c(\rho_{w'})_c$ sont toutes deux des densités de $\omega$ relativement à la trace brute
($w_c=1$), et par unicité $w_c(\rho_w)_c=w'_c(\rho_{w'})_c$. $\square$

La bibliothèque stocke un état par sa densité $\rho=\rho_{\mathrm{blunt}}$ relativement à
la trace brute `Tr_blunt` (classe `DensityMatrix`) ; `density(trace)` et
`DensityMatrix.from_density` effectuent les conversions. La valeur de l'état ne dépend pas
de la convention :

```python
rho = states.random_density_matrix(alg, batch_size=4)          # Hilbert-Schmidt random states
for trace, tr in (('blunt', lambda X: X.Tr_blunt()), ('norm', lambda X: X.Tr_norm()),
                  ('tau_vN', lambda X: X.tau_vN())):
    d = alg.operator(rho.density(trace))
    assert (tr(d @ A) - (rho @ A).Tr_blunt()).abs().max() < 1e-10       # omega(A) for all traces
    assert (tr(d) - 1).abs().max() < 1e-10                               # tr_w(rho_w) = 1
    back = DensityMatrix.from_density(alg, rho.density(trace), trace=trace)
    assert (back.matrix - rho.matrix).abs().max() < 1e-12
```

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). $\omega(A)$ est la valeur moyenne de l'observable $A$. La
> densité $\rho$ relative à la trace brute est la matrice densité usuelle de la mécanique
> quantique sur $\bigoplus_c\mathbb C^{k_c}$. Les densités relatives à d'autres traces
> apparaissent lorsque la même algèbre agit avec des multiplicités (chapitre 2), par
> exemple sur la partie symétrique et la partie non symétrique d'un système de spins.

## 3.3 La partie classique d'un état

La restriction de $\omega$ au centre est une loi de probabilité sur les secteurs,
$p_c=\omega(P_c)=\operatorname{Tr}\rho_c$. Si $p_c>0$, l'**état conditionnel** dans le
secteur $c$ a pour densité $\sigma_c=\rho_c/p_c$ sur $\mathbb C^{k_c}$, et
$$
\rho=\bigoplus_cp_c\,\sigma_c
$$
(si $p_c=0$, alors $\rho_c=0$, puisque les valeurs propres de $\rho_c\ge0$ sont $\ge0$ et de
somme $\operatorname{Tr}\rho_c=0$, corollaire 1.4(3), donc toutes nulles, et que $\rho_c$ est
unitairement diagonalisable, théorème 1.3 ; le terme est omis).
Ainsi, tout état sur $\mathcal M$ est une combinaison convexe d'états des secteurs.
L'application $\rho\mapsto\sum_cP_c\rho P_c$ est l'identité sur $\mathcal M$, car tout élément de
$\mathcal M$ est diagonal par blocs.

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Un état
> sur $\mathcal M$ est un mélange classique d'états quantiques dans les secteurs. L'application
> $\rho\mapsto\sum_cP_c\rho P_c$ est la mesure de « quel secteur », et elle ne perturbe aucun
> état de $\mathcal M$. C'est là le contenu opérationnel d'une règle de supersélection.

```python
p = rho.sector_probabilities()                                   # shape (batch, C)
assert torch.allclose(p.sum(-1), torch.ones(4, dtype=p.dtype))
P = [alg.central([1.0 if c == d else 0.0 for d in range(alg.C)]) for c in range(alg.C)]
measure_sector = channels.measurement_channel(P)
assert (measure_sector(rho).matrix - rho.matrix).abs().max() < 1e-12   # no disturbance
```

## 3.4 Entropie

L'**entropie de von Neumann** d'un état de densité (brute) $\rho$ est
$S(\rho)=-\operatorname{Tr}\rho\log\rho=-\sum_i\lambda_i\log\lambda_i$, la somme portant sur
les valeurs propres de $\rho$, avec la convention $0\log0=0$ ; autrement dit,
$\rho\log\rho:=f(\rho)$ pour $f(t)=t\log t$, $f(0)=0$, et
$\operatorname{Tr}f(\rho)=\sum_if(\lambda_i)$ par définition de $f(\rho)=U\,\mathrm{diag}(f(\lambda_i))\,U^*$
(section 1.3).

**Théorème 3.3 (entropie d'un mélange de secteurs).** Avec $p_c$ et $\sigma_c$ comme
ci-dessus (les secteurs tels que $p_c=0$ étant omis),
$$
S(\rho)=H(p)+\sum_cp_c\,S(\sigma_c),\qquad H(p)=-\sum_cp_c\log p_c .
$$

*Démonstration.* Les valeurs propres non nulles de $\rho$ (avec multiplicité) sont les nombres $p_c\mu$, où
$\mu$ parcourt les valeurs propres non nulles des $\sigma_c$ (secteurs tels que $p_c>0$) ;
les valeurs propres nulles contribuent pour $0$ aux deux membres. Par conséquent
$-\sum_c\sum_{\mu>0}p_c\mu\log(p_c\mu)=-\sum_cp_c\log p_c\sum_{\mu>0}\mu-\sum_cp_c\sum_{\mu>0}\mu\log\mu$,
et $\sum_{\mu>0}\mu=\operatorname{Tr}\sigma_c=1$. $\square$

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). La formule sépare l'incertitude sur le système en deux parts :
> $H(p)$, l'incertitude classique sur la charge (le secteur où se trouve le système), et
> $\sum_cp_cS(\sigma_c)$, l'incertitude quantique moyenne à l'intérieur des secteurs. Sous une règle
> de supersélection, la première part est l'entropie d'une variable aléatoire classique : la
> charge.

L'entropie relativement à une autre trace est $S_w(\rho)=-\operatorname{tr}_w(\rho_w\log\rho_w)$,
avec la même convention $\rho_w\log\rho_w:=f(\rho_w)$. D'après le théorème 3.2, $(\rho_w)_c=\rho_c/w_c$, donc
$$
S_w(\rho)=-\sum_cw_c\operatorname{Tr}\Big(\frac{\rho_c}{w_c}\log\frac{\rho_c}{w_c}\Big)
=S(\rho)+\sum_cp_c\log w_c .
$$
L'entropie est une propriété de l'état *et* de la trace : un même état a des entropies
différentes relativement à `Tr_blunt` et à `tau_vN`.

```python
S = rho.von_neumann_entropy()
sigma = [rho.matrix[:, c, :k, :k] / p[:, c, None, None] for c, k in enumerate(alg.k_factors)]
S_sec = [-(w * torch.log(w)).sum(-1) for w in (torch.linalg.eigvalsh(s) for s in sigma)]
H = -(p * torch.log(p)).sum(-1)
assert torch.allclose(S, H + sum(p[:, c] * S_sec[c] for c in range(alg.C)))    # Theorem 3.3

w_tau = torch.tensor([1 / (alg.C * k) for k in alg.k_factors], dtype=torch.float64)
d_tau = alg.operator(rho.density('tau_vN'), is_positive=True)
S_tau = -(d_tau @ d_tau.log()).tau_vN().real
assert torch.allclose(S_tau, S + (p * torch.log(w_tau)).sum(-1))
```

## 3.5 Entropie relative et inégalité de Klein

Pour des états de densités $\rho,\sigma$, avec $\sigma>0$, l'**entropie relative** est
$S(\rho\|\sigma)=\operatorname{Tr}\rho\log\rho-\operatorname{Tr}\rho\log\sigma$, où
$\operatorname{Tr}\rho\log\rho:=\sum_i\lambda_i\log\lambda_i$, la somme portant sur les
valeurs propres de $\rho$, avec $0\log0=0$ (et de même dans chaque bloc). Elle ne dépend
pas de la trace : d'après le théorème 3.2, $(\rho_w)_c=\rho_c/w_c$ et $(\sigma_w)_c=\sigma_c/w_c$, donc
les valeurs propres de $(\rho_w)_c$ sont celles de $\rho_c$ divisées par $w_c$,
$\log(\sigma_c/w_c)=\log\sigma_c-(\log w_c)1$ (si $\sigma_c=U\,\mathrm{diag}(s_i)\,U^*$, alors
$\sigma_c/w_c=U\,\mathrm{diag}(s_i/w_c)\,U^*$, et théorème 1.5(a)), et
$\operatorname{tr}_w(\rho_w\log\rho_w)$ comme $\operatorname{tr}_w(\rho_w\log\sigma_w)$ diffèrent de
$\operatorname{Tr}\rho\log\rho$ et de $\operatorname{Tr}\rho\log\sigma$ de la même quantité
$-\sum_cp_c\log w_c$ (le logarithme de $\rho$ n'intervient que sur son support, par la convention
$0\log0=0$). La différence est donc inchangée.

**Lemme 3.4 (inégalité de la tangente).** Soit $f:I\to\mathbb R$ convexe sur un intervalle $I$
et dérivable en $b\in I$ (à droite ou à gauche si $b$ est une extrémité de $I$). Alors
$f(a)\ge f(b)+(a-b)f'(b)$ pour tout $a\in I$, avec inégalité stricte pour $a\neq b$ si $f$ est
strictement convexe, c'est-à-dire si $f((1-t)x+ty)<(1-t)f(x)+tf(y)$ pour $x\neq y$ dans $I$ et
$t\in(0,1)$. De plus, si $f$ est continue sur $I$ et deux fois dérivable avec $f''>0$ sur
l'intérieur de $I$, alors $f$ est strictement convexe.

*Démonstration.* Pour $t\in(0,1]$, la convexité donne $f(b+t(a-b))\le(1-t)f(b)+tf(a)$,
c'est-à-dire $\big(f(b+t(a-b))-f(b)\big)/t\le f(a)-f(b)$ ; en faisant tendre $t$ vers $0$,
on obtient l'inégalité. Si $f$ est strictement convexe et $a\neq b$, posons $m=(a+b)/2$ ;
alors $f(a)-f(b)>2\big(f(m)-f(b)\big)\ge2(m-b)f'(b)=(a-b)f'(b)$, la première étape par
stricte convexité en $m$, la seconde par l'inégalité qui vient d'être démontrée, appliquée
à $m$. Pour la dernière assertion, soient $x<y$ dans $I$, $t\in(0,1)$ et $z=(1-t)x+ty$. D'après
le théorème des accroissements finis, il existe $\xi_1\in(x,z)$ et $\xi_2\in(z,y)$ tels que
$f(z)-f(x)=t(y-x)f'(\xi_1)$ et $f(y)-f(z)=(1-t)(y-x)f'(\xi_2)$, et $f'(\xi_1)<f'(\xi_2)$ d'après
le théorème des accroissements finis appliqué à $f'$ sur $[\xi_1,\xi_2]$, puisque $f''>0$ sur cet
intervalle. En multipliant par $t(1-t)(y-x)>0$, on obtient
$(1-t)\big(f(z)-f(x)\big)<t\big(f(y)-f(z)\big)$, c'est-à-dire $f(z)<(1-t)f(x)+tf(y)$. $\square$

**Théorème 3.5 (inégalité de Klein).** Soit $f:I\to\mathbb R$ convexe sur un intervalle
$I$, et $X,Y$ autoadjoints de spectres contenus dans $I$, $f$ étant dérivable en chaque
valeur propre de $Y$ (à droite ou à gauche en une extrémité de $I$). Alors $\operatorname{Tr}[f(X)-f(Y)-(X-Y)f'(Y)]\ge0$. Si $f$ est
strictement convexe, l'égalité n'a lieu que pour $X=Y$.

*Démonstration.* Écrivons $X=\sum_ia_i|u_i\rangle\langle u_i|$ et
$Y=\sum_jb_j|v_j\rangle\langle v_j|$ dans des bases propres orthonormées, et posons
$c_{ij}=|\langle u_i,v_j\rangle|^2$, de sorte que $\sum_ic_{ij}=\sum_jc_{ij}=1$. Alors
$\operatorname{Tr}f(X)=\sum_{ij}c_{ij}f(a_i)$, $\operatorname{Tr}f(Y)=\sum_{ij}c_{ij}f(b_j)$ et
$\operatorname{Tr}Xf'(Y)=\sum_{ij}c_{ij}a_if'(b_j)$, $\operatorname{Tr}Yf'(Y)=\sum_{ij}c_{ij}b_jf'(b_j)$.
Le membre de gauche vaut donc $\sum_{ij}c_{ij}[f(a_i)-f(b_j)-(a_i-b_j)f'(b_j)]$, somme de
termes positifs d'après le lemme 3.4. Si $f$ est strictement convexe et que
la somme est nulle, alors $a_i=b_j$ dès que $c_{ij}>0$, donc
$Xv_j=\sum_ia_i\langle u_i,v_j\rangle u_i=b_jv_j$ pour tout $j$, et $X=Y$. $\square$

**Corollaire 3.6.** $S(\rho\|\sigma)\ge0$, avec égalité si et seulement si $\rho=\sigma$.
En particulier $S(\rho)\le\log D$, avec égalité si et seulement si $\rho=1/D$.

*Démonstration.* On applique le théorème 3.5 à $f(t)=t\log t$, $f(0)=0$, qui est
strictement convexe sur $I=[0,\infty)$ d'après le lemme 3.4 ($f$ est continue en $0$, car
$t\log t\to0$ quand $t\to0^+$, et $f''(t)=1/t>0$ sur $(0,\infty)$) et dérivable sur $(0,\infty)\supseteq\operatorname{spec}\sigma$ avec
$f'(t)=\log t+1$, et à $X=\rho$, $Y=\sigma$ : le membre de gauche vaut
$\operatorname{Tr}\rho\log\rho-\operatorname{Tr}\sigma\log\sigma-\operatorname{Tr}(\rho-\sigma)(\log\sigma+1)
=S(\rho\|\sigma)$, puisque $\operatorname{Tr}\rho=\operatorname{Tr}\sigma=1$. Pour la seconde
assertion, on prend $\sigma=1/D$ : $0\le S(\rho\|1/D)=-S(\rho)+\log D$. $\square$

**Théorème 3.7 (règle de la chaîne ; monotonie sous l'espérance conditionnelle sur le centre).** Soient
$\rho,\sigma$ des états avec $\sigma>0$, $p_c=\operatorname{Tr}\rho_c$,
$q_c=\operatorname{Tr}\sigma_c>0$, $\hat\sigma_c=\sigma_c/q_c$ et, pour $p_c>0$,
$\hat\rho_c=\rho_c/p_c$. Alors, les sommes portant sur les $c$ tels que $p_c>0$ (les autres
secteurs contribuent pour $0$ en vertu de $0\log0=0$),
$$
S(\rho\|\sigma)=\sum_cp_c\log\frac{p_c}{q_c}+\sum_cp_c\,S(\hat\rho_c\|\hat\sigma_c)
\ \ge\ \sum_cp_c\log\frac{p_c}{q_c}=S(E\rho\|E\sigma),
$$
où $E$ est l'espérance conditionnelle sur le centre (proposition 2.8).

*Démonstration.* Pour $p_c>0$,
$\operatorname{Tr}\rho_c\log\rho_c=p_c\log p_c+p_c\operatorname{Tr}\hat\rho_c\log\hat\rho_c$
(valeurs propres $p_c\mu$, comme dans le théorème 3.3) et
$\log\sigma_c=(\log q_c)1+\log\hat\sigma_c$ (on diagonalise $\hat\sigma_c=U\,\mathrm{diag}(s_i)\,U^*$, de
sorte que $\sigma_c=U\,\mathrm{diag}(q_cs_i)\,U^*$, et on utilise le théorème 1.5(a)), donc
$\operatorname{Tr}\rho_c\log\sigma_c=p_c\log q_c+p_c\operatorname{Tr}\hat\rho_c\log\hat\sigma_c$ et
$\operatorname{Tr}\rho_c\log\rho_c-\operatorname{Tr}\rho_c\log\sigma_c=p_c\log(p_c/q_c)+p_cS(\hat\rho_c\|\hat\sigma_c)$.
Si $p_c=0$, alors $\rho_c=0$ et le secteur contribue pour $0$.
On somme sur $c$ ; la seconde somme est $\ge0$ d'après le corollaire 3.6, appliqué dans
$M_{k_c}$ aux états $\hat\rho_c$ et $\hat\sigma_c>0$. Enfin
$E\rho=\bigoplus_c(p_c/k_c)1_{k_c}$ et $E\sigma=\bigoplus_c(q_c/k_c)1_{k_c}>0$ sont des états ; $E\rho$ est
la densité de l'état $\omega\circ E$ (où $\omega(A)=\operatorname{Tr}\rho A$), puisque
$\operatorname{Tr}(\rho E(A))=\sum_cp_c\operatorname{Tr}(A_c)/k_c=\operatorname{Tr}(E(\rho)A)$. Et
$S(E\rho\|E\sigma)=\sum_ck_c\frac{p_c}{k_c}\big(\log\frac{p_c}{k_c}-\log\frac{q_c}{k_c}\big)=\sum_cp_c\log\frac{p_c}{q_c}$,
car les valeurs propres de $E\rho$ dans le bloc $c$ sont $p_c/k_c$ ($k_c$ fois) et les secteurs
tels que $p_c=0$ contribuent pour $0$. $\square$

Ainsi, remplacer les deux états par leurs images par $E$ ne peut que diminuer l'entropie
relative ; ce qui subsiste est la divergence de Kullback–Leibler $\sum_cp_c\log(p_c/q_c)$ entre
les lois des secteurs.

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours).
> L'entropie relative mesure à quel point deux états peuvent être distingués. Oublier
> l'information quantique contenue dans les secteurs ne peut que rendre deux états plus
> difficiles à distinguer ; ce qui subsiste est la discernabilité classique des lois de la charge.

```python
sigma_state = states.random_density_matrix(alg, batch_size=4)
q = sigma_state.sector_probabilities()
D_rel = states.relative_entropy(rho, sigma_state)
kl = (p * torch.log(p / q)).sum(-1)
E = channels.center_expectation(alg)
assert (D_rel >= kl - 1e-12).all()                                       # Theorem 3.7
assert torch.allclose(states.relative_entropy(E(rho), E(sigma_state)), kl)
assert (rho.von_neumann_entropy() <= math.log(10) + 1e-12).all()         # Corollary 3.6
assert abs(states.maximally_mixed_state(alg).von_neumann_entropy().item() - math.log(10)) < 1e-12
```

## 3.6 Le contraste de Michelson

Pour $A\ge0$, $A\neq0$, le **contraste de Michelson** est
$\Delta(A)=(\lambda_{\max}-\lambda_{\min})/(\lambda_{\max}+\lambda_{\min})\in[0,1]$, les valeurs
propres extrêmes étant prises sur l'ensemble des secteurs (et $\Delta(0)=0$). Il mesure
l'écart entre $A$ et un multiple de l'identité : $\Delta(A)=0$ si et seulement si
$A=\lambda1$.

**Proposition 3.8.** Si $\Phi:\mathcal M\to\mathcal M$ est linéaire, positive et unifère
($\Phi(1)=1$), alors $\Delta(\Phi(A))\le\Delta(A)$ pour tout $A\ge0$.

*Démonstration.* $\lambda_{\min}1\le A\le\lambda_{\max}1$ (corollaire 1.4(3) pour
$A-\lambda_{\min}1$ et $\lambda_{\max}1-A$), et la positivité, la linéarité ainsi que
$\Phi(1)=1$ donnent $\lambda_{\min}1\le\Phi(A)\le\lambda_{\max}1$. Par conséquent $\Phi(A)$ est
autoadjoint (une matrice positive l'est, corollaire 1.4(3)), $\Phi(A)\ge0$ (car
$\lambda_{\min}\ge0$), et d'après le corollaire 1.4(3) $\lambda_{\min}(\Phi(A))\ge\lambda_{\min}$ et
$\lambda_{\max}(\Phi(A))\le\lambda_{\max}$. Si $\Phi(A)=0$, l'assertion est triviale. Sinon $A\neq0$,
$\lambda_{\max}$ et $\lambda_{\max}(\Phi(A))$ sont tous deux $>0$, et en posant
$r=\lambda_{\min}/\lambda_{\max}\in[0,1]$,
$\Delta=(1-r)/(1+r)$ est une fonction décroissante de $r$, et $r(\Phi(A))\ge r(A)$. $\square$

```python
A_pos = alg.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(4, d), batch_size=4, force_positive=True)
Phi = channels.random_mixed_unitary_channel(alg, 3, batch_size=4)       # unital
assert (Phi(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
assert (E(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
```

## 3.7 Une inégalité plausible mais fausse

Pour $X,Y$ positifs dans $M_k$ et $U$ unitaire, a-t-on
$|\operatorname{Tr}(XUY)|\le\operatorname{Tr}(XY)$ ? C'est vrai pour $U=1$, et cela ressemble à
un énoncé de type Cauchy–Schwarz. C'est pourtant faux dès que $X$ et $Y$ ne commutent pas.

**Proposition 3.9.** Soient $X,Y\ge0$ dans $M_k$. Alors $\max_U|\operatorname{Tr}(XUY)|=\|YX\|_1\ge\operatorname{Tr}(XY)$, le maximum
portant sur les unitaires $U\in M_k$, avec égalité si et seulement si $XY=YX$.

*Démonstration.* $\operatorname{Tr}(XUY)=\operatorname{Tr}(YX\,U)$, donc le maximum vaut
$\|YX\|_1$ d'après le théorème 1.9, et $\|YX\|_1\ge|\operatorname{Tr}(YX)|=\operatorname{Tr}(XY)$,
puisque $\operatorname{Tr}(XY)=\operatorname{Tr}(X^{1/2}YX^{1/2})\ge0$, car
$X^{1/2}YX^{1/2}=(Y^{1/2}X^{1/2})^*(Y^{1/2}X^{1/2})\ge0$ (corollaire 1.7). Si $XY=YX$, alors
$X^{1/2}$ commute avec $Y$ (théorème 1.5(d)), donc $YX=X^{1/2}YX^{1/2}\ge0$ et
$\|YX\|_1=\operatorname{Tr}(YX)$ (pour $P\ge0$, $\|P\|_1=\operatorname{Tr}P$, section 1.4).
Réciproquement, supposons $\|YX\|_1=\operatorname{Tr}(YX)$. Si $|YX|=0$, alors
$YX=0=(YX)^*=XY$. Sinon, écrivons $YX=W|YX|$ (théorème 1.8) et $R=|YX|^{1/2}\neq0$. Par
cyclicité, $\operatorname{Tr}(YX)=\operatorname{Tr}(RWR)=\langle R,WR\rangle_{HS}$, tandis que
$\|R\|_2^2=\|WR\|_2^2=\operatorname{Tr}|YX|=\|YX\|_1$. Donc $\langle R,WR\rangle_{HS}=\|R\|_2\|WR\|_2$ :
il y a égalité dans l'inégalité de Cauchy–Schwarz, donc $WR=\lambda R$ pour un scalaire
$\lambda$, et $\lambda\|R\|_2^2=\langle R,WR\rangle_{HS}=\|R\|_2^2$ donne $\lambda=1$. Ainsi
$W|YX|^{1/2}=|YX|^{1/2}$ et $YX=W|YX|=|YX|\ge0$. Une matrice positive est autoadjointe, donc
$YX=(YX)^*=XY$. $\square$

La bibliothèque trouve le point où le maximum est atteint grâce à la décomposition polaire :
d'après le théorème 1.9, pour $YX=W|YX|$, l'unitaire $U=W^*$ atteint $\|YX\|_1$. Une recherche
aléatoire fait moins bien en moyenne : pour un $U$ distribué selon la loi de Haar
(théorème 5.3 ; c'est ainsi que la bibliothèque tire $U$), la proposition 11.11 avec $A=YX$ donne
$\big(\mathbb E|\operatorname{Tr}(XUY)|^2\big)^{1/2}=\|YX\|_2/\sqrt k\le\|YX\|_1/\sqrt k=\max_V|\operatorname{Tr}(XVY)|/\sqrt k$
(car $\sum_is_i^2\le(\sum_is_i)^2$), et le corollaire 11.12 montre que, pour $X,Y$ de spectres
contenus dans $[a,b]$, $a>0$, un seul tirage atteint $\operatorname{Tr}(XY)$ avec une probabilité au
plus égale à $b^4/(a^4k^2)$ (ce qui n'est informatif que pour $k>b^2/a^2$ ; dans l'exécution
ci-dessous, $b/a=6$ et $k=6$, de sorte que cela ne dit rien). Le code calcule le supremum exact
pour un couple qui ne commute pas avec $k=6$ ainsi que le meilleur de $20000$ tirages de Haar,
tous deux diminués de $\operatorname{Tr}(XY)$. D'après la proposition 3.9, aucun tirage ne peut
dépasser le supremum ; l'écart entre le meilleur tirage et le supremum est une observation
numérique de cette exécution, pas un théorème.

```python
Mk = TypeIAlgebra([6], [6], precision='double')
X = Mk.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
Y = Mk.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
assert (X @ Y - Y @ X).matrix.abs().max().item() > 1e-6             # X and Y do not commute
YX = Y @ X
W = YX @ YX.abs().inverse()                                          # polar part of YX
gap = (X @ W.adjoint() @ Y).trace.abs() - (X @ Y).trace.real
assert abs(gap.item() - (YX.trace_norm() - (X @ Y).trace.real).item()) < 1e-10   # Theorem 1.9
assert gap.item() > 0                                  # Proposition 3.9: the inequality fails
samples = Mk.random_unitary_operator(20_000)
best = (X @ samples @ Y).trace.abs().max() - (X @ Y).trace.real
print(f"sup - Tr(XY) = {gap.item():.4f}; numerical observation: best of 20000 Haar samples "
      f"- Tr(XY) = {best.item():.4f}")
assert best.item() <= gap.item() + 1e-10               # Proposition 3.9: no sample beats the sup
```

## 3.8 Points de comparaison

QuTiP (`entropy_vn`, `entropy_relative`) et d'autres bibliothèques d'information quantique
(par exemple `toqito`) calculent les entropies et entropies relatives de matrices densité sur
un seul espace de Hilbert, avec la trace de cet espace ; sur $\bigoplus_c\mathbb C^{k_c}$, c'est
la trace brute d'ici (à la base du logarithme près, qui varie d'une bibliothèque à l'autre).
Ce qu'elles ne fournissent pas, c'est
la famille de traces du théorème 3.1 et les conversions du théorème 3.2, qui deviennent
nécessaires dès qu'une algèbre agit avec des multiplicités ou que les secteurs doivent être
pondérés différemment.

## Exercices

1. Montrer que la densité de l'état tracial $\tau_{\mathrm{vN}}$ (poids égaux sur les
   secteurs) relativement à $\tau_{\mathrm{vN}}$ est $1$, puis calculer sa densité brute et
   son entropie $S$.
2. Soit $\operatorname{tr}_w$ fidèle. Montrer que $S_w(\rho)\le\log\operatorname{tr}_w(1)$, avec
   égalité si et seulement si $\rho_w=1/\operatorname{tr}_w(1)$.
3. Montrer que le canal de mesure $\rho\mapsto\sum_cP_c\rho P_c$ sur $M_D$ (et non sur
   $\mathcal M$ !) perturbe bel et bien les états, et calculer $S(\rho)$ avant et après la
   mesure pour une superposition de deux secteurs.

<details>
<summary>Solutions</summary>

1. $\tau_{\mathrm{vN}}(A)=\tau_{\mathrm{vN}}(1\cdot A)$, donc $\rho_{\tau}=1$ par unicité.
   D'après le théorème 3.2, la densité brute est $\bigoplus_c1_{k_c}/(Ck_c)$, avec
   $p_c=1/C$ et $\sigma_c=1/k_c$, donc $S=\log C+\frac1C\sum_c\log k_c$ d'après le
   théorème 3.3.

```python
tau = states.tracial_state(alg)
assert (tau.density('tau_vN') - one.matrix).abs().max() < 1e-12
assert abs(tau.von_neumann_entropy().item() - (math.log(3) + sum(math.log(k) for k in (2, 3, 5)) / 3)) < 1e-12
```

2. Les sommes portent sur les $c$ tels que $p_c>0$. Alors
   $S_w(\rho)=S(\rho)+\sum_cp_c\log w_c=-\sum_cp_c\log p_c+\sum_cp_cS(\sigma_c)+\sum_cp_c\log w_c
   \le\sum_cp_c\log\frac{k_cw_c}{p_c}\le\log\sum_{p_c>0}k_cw_c\le\log\operatorname{tr}_w(1)$, en
   utilisant $S(\sigma_c)\le\log k_c$ (corollaire 3.6 dans $M_{k_c}$) et l'inégalité de Jensen finie
   $\sum_{c=1}^mp_c\log x_c\le\log\sum_{c=1}^mp_cx_c$ ($p_c>0$, $\sum_cp_c=1$, $x_c>0$), avec
   égalité seulement si tous les $x_c$ sont égaux, appliquée à $x_c=k_cw_c/p_c$. Preuve par
   récurrence sur $m$ : pour $m=1$, il n'y a rien à démontrer. Pour $m\ge2$, posons
   $q=1-p_m\in(0,1)$ et $\bar x=\sum_{c<m}(p_c/q)x_c$ ; par l'hypothèse de récurrence,
   $\sum_{c<m}(p_c/q)\log x_c\le\log\bar x$, avec égalité seulement si $x_1=\dots=x_{m-1}$, et par
   stricte concavité de $\log$ (lemme 3.4 pour $-\log$, dont la dérivée seconde est $1/t^2>0$)
   $q\log\bar x+p_m\log x_m\le\log(q\bar x+p_mx_m)$,
   avec égalité seulement si $\bar x=x_m$. En multipliant la première inégalité par $q$ et en
   additionnant, $\sum_cp_c\log x_c\le\log\sum_cp_cx_c$, et l'égalité impose
   $x_1=\dots=x_{m-1}=\bar x=x_m$. L'égalité partout impose que tous les $p_c$ soient $>0$ (car
   tous les $k_cw_c>0$), que
   $\sigma_c=1/k_c$ et que $k_cw_c/p_c$ soit constant, c'est-à-dire
   $p_c=k_cw_c/\operatorname{tr}_w(1)$, soit encore $(\rho_w)_c=\rho_c/w_c=1_{k_c}/\operatorname{tr}_w(1)$.
   Réciproquement, pour $\rho_w=1/\operatorname{tr}_w(1)$,
   $S_w(\rho)=-\operatorname{tr}_w(1)\cdot\frac1{\operatorname{tr}_w(1)}\log\frac1{\operatorname{tr}_w(1)}=\log\operatorname{tr}_w(1)$.

3. Dans $M_D$, une superposition pure $\psi=(\psi_1\oplus\psi_2)/\sqrt2$ de vecteurs
   unitaires de deux secteurs a pour entropie $S=0$ ; après la mesure, l'état est
   $\frac12|\psi_1\rangle\langle\psi_1|\oplus\frac12|\psi_2\rangle\langle\psi_2|$, d'entropie
   $S=\log2$. Restreint à $\mathcal M$, en revanche, l'état vectoriel associé à $\psi$
   coïncide avec le mélange :
   $\langle\psi,A\psi\rangle=\tfrac12\langle\psi_1,A_1\psi_1\rangle+\tfrac12\langle\psi_2,A_2\psi_2\rangle$
   pour tout $A\in\mathcal M$. Aucun élément de $\mathcal M$ ne détecte la cohérence
   (physiquement, c'est le contenu d'une règle de supersélection).

</details>

**Suite :** le [chapitre 4](04_symmetries.md) traite des symétries et des représentations :
comment l'action d'un groupe sur un système physique fait apparaître les secteurs, et
comment la bibliothèque construit les algèbres correspondantes pour les systèmes de spins
et les registres de qubits.
