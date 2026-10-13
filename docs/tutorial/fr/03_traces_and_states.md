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
$E^c_{11}=(E^c_{11})^*E^c_{11}$, donc la condition $w_c\ge0$ est nécessaire, et
$\operatorname{tr}_w(A^*A)=\sum_cw_c\|A_c\|_2^2$ montre qu'elle est suffisante, ainsi que le
critère de fidélité. $\square$

La bibliothèque implémente trois traces fidèles :

| nom | poids $w_c$ | valeur en $1$ | signification |
|---|---|---|---|
| `Tr_blunt` | $1$ | $D$ | trace de $\bigoplus_c\mathbb C^{k_c}$ (multiplicités $m_c=1$) |
| `Tr_norm` | $1/k_c$ | $C$ | chaque bloc normalisé |
| `tau_vN` | $1/(Ck_c)$ | $1$ | l'état tracial de poids égaux $1/C$ sur les secteurs |

D'après le chapitre 2, la trace d'un espace sur lequel $\mathcal M$ agit avec les
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
puisque $A_c(\rho_w)_cA_c^*\ge0$ (corollaire 1.7 avec $T=(\rho_w)_c^{1/2}A_c^*$).
Réciproquement, si $\omega$ est positive, prenons $x\in\mathbb C^{k_c}$ et pour $A$
l'opérateur de rang un $|x\rangle\langle x|$ dans le bloc $c$ ; alors $A=A^*A/\|x\|^2$ et
$0\le\omega(A)=w_c\langle x,(\rho_w)_cx\rangle$, donc $(\rho_w)_c\ge0$. La normalisation
s'écrit $\omega(1)=\operatorname{tr}_w(\rho_w)$. Enfin,
$\sum_cw_c\operatorname{Tr}((\rho_w)_cA_c)$ et $\sum_cw'_c\operatorname{Tr}((\rho_{w'})_cA_c)$
sont tous deux égaux à $\omega(A)$ pour tout $A$, et par unicité
$w_c(\rho_w)_c=w'_c(\rho_{w'})_c$. $\square$

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

> **Interprétation physique.** $\omega(A)$ est la valeur moyenne de l'observable $A$. La
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
(si $p_c=0$, alors $\rho_c=0$, puisque $\rho_c\ge0$ et $\operatorname{Tr}\rho_c=0$, et le
terme est omis).
Un état sur $\mathcal M$ est un mélange classique d'états quantiques dans les secteurs.
Mesurer « quel secteur » ne le perturbe pas : le canal de mesure
$\rho\mapsto\sum_cP_c\rho P_c$ est l'identité sur $\mathcal M$, car tout élément de
$\mathcal M$ est diagonal par blocs. C'est là le contenu opérationnel d'une règle de
supersélection.

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
les valeurs propres de $\rho$, avec la convention $0\log0=0$.

**Théorème 3.3 (entropie d'un mélange de secteurs).** Avec $p_c$ et $\sigma_c$ comme
ci-dessus (les secteurs tels que $p_c=0$ étant omis),
$$
S(\rho)=H(p)+\sum_cp_c\,S(\sigma_c),\qquad H(p)=-\sum_cp_c\log p_c .
$$

*Démonstration.* Les valeurs propres non nulles de $\rho$ sont les nombres $p_c\mu$, où
$\mu$ parcourt les valeurs propres non nulles des $\sigma_c$ (secteurs tels que $p_c>0$) ;
les valeurs propres nulles contribuent pour $0$ aux deux membres. Par conséquent
$-\sum_c\sum_{\mu>0}p_c\mu\log(p_c\mu)=-\sum_cp_c\log p_c\sum_{\mu>0}\mu-\sum_cp_c\sum_{\mu>0}\mu\log\mu$,
et $\sum_{\mu>0}\mu=\operatorname{Tr}\sigma_c=1$. $\square$

> **Interprétation physique.** La formule sépare l'incertitude sur le système en deux parts :
> $H(p)$, l'incertitude classique sur la charge (le secteur où se trouve le système), et
> $\sum_cp_cS(\sigma_c)$, l'incertitude quantique moyenne à l'intérieur des secteurs. Sous une règle
> de supersélection, la première part est l'entropie d'une variable aléatoire classique : la
> charge.

L'entropie relativement à une autre trace est $S_w(\rho)=-\operatorname{tr}_w(\rho_w\log\rho_w)$.
D'après le théorème 3.2, $(\rho_w)_c=\rho_c/w_c$, donc
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
les valeurs propres de $(\rho_w)_c$ sont celles de $\rho_c$ divisées par $w_c$, et
$\operatorname{tr}_w(\rho_w\log\rho_w)$ comme $\operatorname{tr}_w(\rho_w\log\sigma_w)$ diffèrent de
$\operatorname{Tr}\rho\log\rho$ et de $\operatorname{Tr}\rho\log\sigma$ de la même quantité
$-\sum_cp_c\log w_c$ (le logarithme de $\rho$ n'intervient que sur son support, par la convention
$0\log0=0$). La différence est donc inchangée.

**Lemme 3.4 (inégalité de la tangente).** Soit $f:I\to\mathbb R$ convexe sur un intervalle $I$
et dérivable en $b\in I$. Alors $f(a)\ge f(b)+(a-b)f'(b)$ pour tout $a\in I$, avec
inégalité stricte pour $a\neq b$ si $f$ est strictement convexe.

*Démonstration.* Pour $t\in(0,1]$, la convexité donne $f(b+t(a-b))\le(1-t)f(b)+tf(a)$,
c'est-à-dire $\big(f(b+t(a-b))-f(b)\big)/t\le f(a)-f(b)$ ; en faisant tendre $t$ vers $0$,
on obtient l'inégalité. Si $f$ est strictement convexe et $a\neq b$, posons $m=(a+b)/2$ ;
alors $f(a)-f(b)>2\big(f(m)-f(b)\big)\ge2(m-b)f'(b)=(a-b)f'(b)$, la première étape par
stricte convexité en $m$, la seconde par l'inégalité qui vient d'être démontrée, appliquée
à $m$. $\square$

**Théorème 3.5 (inégalité de Klein).** Soit $f:I\to\mathbb R$ convexe sur un intervalle
$I$, et $X,Y$ autoadjoints de spectres contenus dans $I$, $f$ étant dérivable en chaque
valeur propre de $Y$. Alors $\operatorname{Tr}[f(X)-f(Y)-(X-Y)f'(Y)]\ge0$. Si $f$ est
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
strictement convexe sur $I=[0,\infty)$ ($f''(t)=1/t>0$ sur $(0,\infty)$, et $f$ est
continue en $0$) et dérivable sur $(0,\infty)\supseteq\operatorname{spec}\sigma$ avec
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
$\log\sigma_c=(\log q_c)1+\log\hat\sigma_c$, donc
$\operatorname{Tr}\rho_c\log\sigma_c=p_c\log q_c+p_c\operatorname{Tr}\hat\rho_c\log\hat\sigma_c$ et
$\operatorname{Tr}\rho_c\log\rho_c-\operatorname{Tr}\rho_c\log\sigma_c=p_c\log(p_c/q_c)+p_cS(\hat\rho_c\|\hat\sigma_c)$.
Si $p_c=0$, alors $\rho_c=0$ et le secteur contribue pour $0$.
On somme sur $c$ ; la seconde somme est $\ge0$ d'après le corollaire 3.6 appliqué dans
chaque bloc. Enfin $E\rho=\bigoplus_c(p_c/k_c)1_{k_c}$, $E\sigma=\bigoplus_c(q_c/k_c)1_{k_c}$,
et $\operatorname{Tr}E\rho(\log E\rho-\log E\sigma)=\sum_ck_c\frac{p_c}{k_c}\log\frac{p_c}{q_c}$. $\square$

Ainsi, oublier l'information quantique contenue dans les secteurs ne peut que diminuer la
discernabilité de deux états ; ce qui subsiste est la divergence de Kullback–Leibler
classique entre les lois des secteurs.

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

*Démonstration.* $\lambda_{\min}1\le A\le\lambda_{\max}1$, et la positivité ainsi que
$\Phi(1)=1$ donnent $\lambda_{\min}1\le\Phi(A)\le\lambda_{\max}1$. Par conséquent
$\lambda_{\min}(\Phi(A))\ge\lambda_{\min}$ et $\lambda_{\max}(\Phi(A))\le\lambda_{\max}$. Si
$\Phi(A)=0$, l'assertion est triviale. Sinon, en posant $r=\lambda_{\min}/\lambda_{\max}$,
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
$YX=0=(YX)^*=XY$. Sinon, écrivons $YX=W|YX|$. D'après le théorème 1.9, l'égalité
$\operatorname{Tr}(YX\cdot1)=\|YX\|_1$ signifie qu'il y a égalité dans l'étape de
Cauchy–Schwarz pour $U=1$ : $|YX|^{1/2}$ et $W|YX|^{1/2}$ sont proportionnels avec un
facteur positif, qui vaut $1$ par comparaison des normes, donc $W|YX|^{1/2}=|YX|^{1/2}$ et
$YX=W|YX|=|YX|\ge0$. Une matrice positive est autoadjointe, donc $YX=(YX)^*=XY$. $\square$

La bibliothèque trouve le point où le maximum est atteint grâce à la décomposition polaire.
L'article qui accompagne la bibliothèque montre numériquement qu'un échantillonnage de
Monte-Carlo de $U$ n'y parvient pas, et ce d'autant moins que $k$ est grand — c'est l'un de
ses thèmes centraux ; le code ci-dessous l'illustre pour un couple. La raison en est que $U=1$ est déjà proche de
l'optimum, tandis que pour un $U$ tiré selon la loi de Haar en grande dimension,
$|\operatorname{Tr}(XUY)|$ est typiquement bien plus petit que $\operatorname{Tr}(XY)$ (un unitaire
aléatoire brouille les phases) : tous les échantillons restent loin en dessous du supremum exact.

```python
Mk = TypeIAlgebra([6], [6], precision='double')
X = Mk.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
Y = Mk.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
YX = Y @ X
W = YX @ YX.abs().inverse()                                          # polar part of YX
gap = (X @ W.adjoint() @ Y).trace.abs() - (X @ Y).trace.real
assert abs(gap.item() - (YX.trace_norm() - (X @ Y).trace.real).item()) < 1e-10
assert gap.item() > 0                                                # the inequality fails
samples = Mk.random_unitary_operator(20_000)
best = (X @ samples @ Y).trace.abs().max() - (X @ Y).trace.real
print(f"exact sup {gap.item():.4f}, best of 20000 Haar samples {best.item():.4f}")
assert best.item() < gap.item()
```

## 3.8 Points de comparaison

QuTiP (`entropy_vn`, `entropy_relative`) et les bibliothèques d'information quantique (par
exemple `qutip`, `toqito`) calculent les entropies et entropies relatives de matrices
densité sur un seul espace de Hilbert ; sur $\bigoplus_c\mathbb C^{k_c}$, elles coïncident
avec les conventions de la trace brute adoptées ici. Ce qu'elles ne fournissent pas, c'est
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
   utilisant $S(\sigma_c)\le\log k_c$ (corollaire 3.6) et l'inégalité de Jensen finie
   $\sum_cp_c\log x_c\le\log\sum_cp_cx_c$ pour la fonction concave $\log$ (qui découle de la
   définition de la concavité par récurrence sur le nombre de termes), avec
   $x_c=k_cw_c/p_c$ ; elle est stricte sauf si tous les $x_c$ sont égaux, puisque $\log$ est
   strictement concave. L'égalité partout impose que tous les $p_c$ soient $>0$, que
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
   pour tout $A\in\mathcal M$. Aucune observable de $\mathcal M$ ne détecte la cohérence :
   c'est le contenu d'une règle de supersélection.

</details>

**Suite :** le chapitre 4 (en préparation) traite des symétries et des représentations :
comment l'action d'un groupe sur un système physique fait apparaître les secteurs, et
comment la bibliothèque construit les algèbres correspondantes pour les systèmes de spins
et les registres de qubits.
