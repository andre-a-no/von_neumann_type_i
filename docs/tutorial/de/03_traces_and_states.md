[English](../03_traces_and_states.md) · [Español](../es/03_traces_and_states.md) · [Français](../fr/03_traces_and_states.md) · **Deutsch** · [中文](../zh/03_traces_and_states.md) · [日本語](../ja/03_traces_and_states.md) · [Русский](../ru/03_traces_and_states.md)

# 3. Spuren, Zustände und Entropie

Auf $M_n(\mathbb C)$ gibt es im Wesentlichen nur eine Spur. Auf $\bigoplus_cM_{k_c}$ gibt es
eine $C$-parametrige Familie, und die Wahl ist von Bedeutung: Dichten, Entropien und die
dualen Abbildungen hängen alle von ihr ab. Dieses Kapitel klassifiziert die Spuren,
stellt Zustände durch Dichten dar, beweist die grundlegenden Entropieungleichungen und
endet mit einer plausiblen Spurungleichung, die falsch ist.

**Sie lernen** die drei Spuren der Bibliothek und die Umrechnung zwischen ihnen kennen,
den klassischen Teil eines Zustands, die Zerlegung der von-Neumann-Entropie nach Sektoren,
die Kleinsche Ungleichung und ihre Folgerungen sowie den Michelson-Kontrast.

Im Folgenden sei $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ mit allen $k_c\ge1$,
$D=\sum_ck_c$, $P_c$ die minimalen zentralen Projektionen und $A=\bigoplus_cA_c$.

## 3.1 Spuren

Eine **Spur** auf $\mathcal M$ ist ein lineares Funktional $\tau$ mit $\tau(AB)=\tau(BA)$ für
alle $A,B$. Sie heißt **positiv**, wenn $\tau(A^*A)\ge0$ für alle $A$ gilt, und **treu**, wenn
außerdem aus $\tau(A^*A)=0$ stets $A=0$ folgt.

**Satz 3.1.** Jede Spur auf $\mathcal M$ hat die Form
$\operatorname{tr}_w(A)=\sum_cw_c\operatorname{Tr}A_c$ mit $w\in\mathbb C^C$. Sie ist genau dann
positiv, wenn alle $w_c\ge0$, und genau dann treu, wenn alle $w_c>0$ sind.

*Beweis.* Seien $E^c_{ij}$ die Matrixeinheiten des $c$-ten Blocks. Für eine Spur $\tau$ gilt
$\tau(E^c_{ij})=\tau(E^c_{i1}E^c_{1j})=\tau(E^c_{1j}E^c_{i1})=\delta_{ij}\tau(E^c_{11})$.
Man setze $w_c=\tau(E^c_{11})$; wegen der Linearität ist $\tau(A)=\sum_c\sum_i(A_c)_{ii}w_c$.
Umgekehrt ist jedes $\operatorname{tr}_w$ eine Spur. Positivität: $E^c_{11}=(E^c_{11})^*E^c_{11}$,
also ist $w_c\ge0$ notwendig, und $\operatorname{tr}_w(A^*A)=\sum_cw_c\|A_c\|_2^2$ zeigt, dass
es hinreichend ist, und liefert das Kriterium für Treue. $\square$

Die Bibliothek implementiert drei treue Spuren:

| Name | Gewichte $w_c$ | Wert auf $1$ | Bedeutung |
|---|---|---|---|
| `Tr_blunt` | $1$ | $D$ | Spur von $\bigoplus_c\mathbb C^{k_c}$ (Vielfachheiten $m_c=1$) |
| `Tr_norm` | $1/k_c$ | $C$ | jeder Block normiert |
| `tau_vN` | $1/(Ck_c)$ | $1$ | der Spurzustand mit gleichen Sektorgewichten $1/C$ |

Nach Kapitel 2 ist die Spur eines Raumes, auf dem $\mathcal M$ mit Vielfachheiten $m_c$
wirkt, gleich $\operatorname{tr}_w$ mit $w_c=m_c$. Eine treue Spur ist also die Festlegung,
„wie viel jeder Sektor zählt“.

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

## 3.2 Zustände und Dichten

Ein **Zustand** ist ein lineares Funktional $\omega$ auf $\mathcal M$, das positiv
($\omega(A^*A)\ge0$) und normiert ($\omega(1)=1$) ist.

**Satz 3.2.** Sei $\operatorname{tr}_w$ eine treue Spur. Zu jedem linearen Funktional
$\omega$ gibt es genau ein $\rho_w\in\mathcal M$ mit $\omega(A)=\operatorname{tr}_w(\rho_wA)$ für
alle $A$. Das Funktional $\omega$ ist genau dann positiv, wenn $\rho_w\ge0$, und genau dann
ein Zustand, wenn außerdem $\operatorname{tr}_w(\rho_w)=1$ gilt. Für zwei treue Spuren hängen
die Dichten blockweise zusammen durch
$$
w_c\,(\rho_w)_c=w'_c\,(\rho_{w'})_c .
$$

*Beweis.* $\langle X,Y\rangle_w=\operatorname{tr}_w(X^*Y)$ ist ein Skalarprodukt auf $\mathcal M$
(Satz 3.1). Nach dem Darstellungssatz von Riesz im endlichdimensionalen Raum
$\mathcal M$ gibt es genau ein $R$ mit $\omega(A)=\langle R,A\rangle_w=\operatorname{tr}_w(R^*A)$;
man setze $\rho_w=R^*$. Ist $\rho_w\ge0$, so gilt
$\omega(A^*A)=\sum_cw_c\operatorname{Tr}(A_c\,(\rho_w)_c\,A_c^*)\ge0$, da
$A_c(\rho_w)_cA_c^*\ge0$ (Korollar 1.7 mit $T=(\rho_w)_c^{1/2}A_c^*$). Ist umgekehrt
$\omega$ positiv, so nehme man $x\in\mathbb C^{k_c}$ und als $A$ den Rang-eins-Operator
$|x\rangle\langle x|$ im Block $c$; dann ist $A=A^*A/\|x\|^2$ und
$0\le\omega(A)=w_c\langle x,(\rho_w)_cx\rangle$, also $(\rho_w)_c\ge0$. Die Normierung lautet
$\omega(1)=\operatorname{tr}_w(\rho_w)$. Schließlich sind sowohl $\sum_cw_c\operatorname{Tr}((\rho_w)_cA_c)$ als auch
$\sum_cw'_c\operatorname{Tr}((\rho_{w'})_cA_c)$ für alle $A$ gleich $\omega(A)$, und wegen der
Eindeutigkeit gilt $w_c(\rho_w)_c=w'_c(\rho_{w'})_c$. $\square$

Die Bibliothek speichert einen Zustand durch seine Dichte $\rho=\rho_{\mathrm{blunt}}$ bezüglich
`Tr_blunt` (Klasse `DensityMatrix`); `density(trace)` und `DensityMatrix.from_density`
rechnen um. Der Wert des Zustands hängt nicht von der Konvention ab:

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

> **Physikalische Bedeutung.** $\omega(A)$ ist der Erwartungswert der Observablen $A$.
> Die Dichte $\rho$ bezüglich der ungewichteten Spur ist die gewöhnliche Dichtematrix der
> Quantenmechanik auf $\bigoplus_c\mathbb C^{k_c}$. Dichten bezüglich anderer Spuren treten
> auf, wenn dieselbe Algebra mit Vielfachheiten wirkt (Kapitel 2), etwa auf dem
> symmetrischen und dem nichtsymmetrischen Teil eines Spinsystems.

## 3.3 Der klassische Teil eines Zustands

Die Einschränkung von $\omega$ auf das Zentrum ist eine Wahrscheinlichkeitsverteilung auf den
Sektoren, $p_c=\omega(P_c)=\operatorname{Tr}\rho_c$. Ist $p_c>0$, so hat der **bedingte Zustand**
im Sektor $c$ die Dichte $\sigma_c=\rho_c/p_c$ auf $\mathbb C^{k_c}$, und es gilt
$$
\rho=\bigoplus_cp_c\,\sigma_c
$$
(ist $p_c=0$, so ist $\rho_c=0$, da $\rho_c\ge0$ und $\operatorname{Tr}\rho_c=0$, und der Term
entfällt).
Ein Zustand auf $\mathcal M$ ist eine klassische Mischung von Quantenzuständen in den
Sektoren. Die Messung von „welcher Sektor“ stört ihn nicht: Der Messkanal
$\rho\mapsto\sum_cP_c\rho P_c$ ist die Identität auf $\mathcal M$, da jedes Element von
$\mathcal M$ blockdiagonal ist. Das ist der operationale Gehalt einer Superauswahlregel.

```python
p = rho.sector_probabilities()                                   # shape (batch, C)
assert torch.allclose(p.sum(-1), torch.ones(4, dtype=p.dtype))
P = [alg.central([1.0 if c == d else 0.0 for d in range(alg.C)]) for c in range(alg.C)]
measure_sector = channels.measurement_channel(P)
assert (measure_sector(rho).matrix - rho.matrix).abs().max() < 1e-12   # no disturbance
```

## 3.4 Entropie

Die **von-Neumann-Entropie** eines Zustands mit (ungewichteter) Dichte $\rho$ ist
$S(\rho)=-\operatorname{Tr}\rho\log\rho=-\sum_i\lambda_i\log\lambda_i$, summiert über die
Eigenwerte von $\rho$, mit $0\log0=0$.

**Satz 3.3 (Entropie einer Mischung von Sektoren).** Mit $p_c$ und $\sigma_c$ wie oben
(Sektoren mit $p_c=0$ werden weggelassen) gilt
$$
S(\rho)=H(p)+\sum_cp_c\,S(\sigma_c),\qquad H(p)=-\sum_cp_c\log p_c .
$$

*Beweis.* Die von null verschiedenen Eigenwerte von $\rho$ sind die Zahlen $p_c\mu$ zu den
von null verschiedenen Eigenwerten $\mu$ der $\sigma_c$ (Sektoren mit $p_c>0$); Eigenwerte
null tragen zu beiden Seiten $0$ bei. Folglich ist
$-\sum_c\sum_{\mu>0}p_c\mu\log(p_c\mu)=-\sum_cp_c\log p_c\sum_{\mu>0}\mu-\sum_cp_c\sum_{\mu>0}\mu\log\mu$,
und $\sum_{\mu>0}\mu=\operatorname{Tr}\sigma_c=1$. $\square$

> **Physikalische Bedeutung.** Die Formel zerlegt die Unsicherheit über das System in zwei Teile:
> $H(p)$, die klassische Unsicherheit über die Ladung (in welchem Sektor sich das System befindet),
> und $\sum_cp_cS(\sigma_c)$, die mittlere Quantenunsicherheit innerhalb der Sektoren. Unter einer
> Superauswahlregel ist der erste Teil die Entropie einer klassischen Zufallsvariablen: der Ladung.

Die Entropie bezüglich einer anderen Spur ist $S_w(\rho)=-\operatorname{tr}_w(\rho_w\log\rho_w)$.
Nach Satz 3.2 ist $(\rho_w)_c=\rho_c/w_c$, also
$$
S_w(\rho)=-\sum_cw_c\operatorname{Tr}\Big(\frac{\rho_c}{w_c}\log\frac{\rho_c}{w_c}\Big)
=S(\rho)+\sum_cp_c\log w_c .
$$
Die Entropie ist eine Eigenschaft des Zustands *und* der Spur: Derselbe Zustand hat
bezüglich `Tr_blunt` und `tau_vN` verschiedene Entropien.

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

## 3.5 Relative Entropie und die Kleinsche Ungleichung

Für Zustände mit Dichten $\rho,\sigma$, $\sigma>0$, ist die **relative Entropie**
$S(\rho\|\sigma)=\operatorname{Tr}\rho\log\rho-\operatorname{Tr}\rho\log\sigma$, wobei
$\operatorname{Tr}\rho\log\rho:=\sum_i\lambda_i\log\lambda_i$, summiert über die Eigenwerte von
$\rho$, mit $0\log0=0$ (und ebenso in jedem Block). Sie hängt nicht von der Spur ab: Nach Satz 3.2
ist $(\rho_w)_c=\rho_c/w_c$ und $(\sigma_w)_c=\sigma_c/w_c$, die Eigenwerte von $(\rho_w)_c$ sind also
die von $\rho_c$, geteilt durch $w_c$, und sowohl $\operatorname{tr}_w(\rho_w\log\rho_w)$ als auch
$\operatorname{tr}_w(\rho_w\log\sigma_w)$ unterscheiden sich von $\operatorname{Tr}\rho\log\rho$ bzw.
$\operatorname{Tr}\rho\log\sigma$ um denselben Betrag $-\sum_cp_c\log w_c$ (der Logarithmus von $\rho$
wird nach der Konvention $0\log0=0$ nur auf dessen Träger gebraucht). Die Differenz bleibt
unverändert.

**Lemma 3.4 (Tangentenungleichung).** Sei $f:I\to\mathbb R$ konvex auf einem Intervall $I$ und
differenzierbar in $b\in I$. Dann gilt $f(a)\ge f(b)+(a-b)f'(b)$ für alle $a\in I$, und zwar
strikt für $a\neq b$, wenn $f$ strikt konvex ist.

*Beweis.* Für $t\in(0,1]$ liefert die Konvexität $f(b+t(a-b))\le(1-t)f(b)+tf(a)$, d. h.
$\big(f(b+t(a-b))-f(b)\big)/t\le f(a)-f(b)$; der Grenzübergang $t\to0$ ergibt die Ungleichung.
Ist $f$ strikt konvex und $a\neq b$, so sei $m=(a+b)/2$; dann gilt
$f(a)-f(b)>2\big(f(m)-f(b)\big)\ge2(m-b)f'(b)=(a-b)f'(b)$, der erste Schritt wegen strikter
Konvexität in $m$, der zweite nach der eben bewiesenen Ungleichung, angewandt auf $m$. $\square$

**Satz 3.5 (Kleinsche Ungleichung).** Sei $f:I\to\mathbb R$ konvex auf einem Intervall $I$,
und seien $X,Y$ selbstadjungiert mit Spektren in $I$, wobei $f$ in jedem Eigenwert von $Y$
differenzierbar ist. Dann gilt $\operatorname{Tr}[f(X)-f(Y)-(X-Y)f'(Y)]\ge0$. Ist $f$ strikt
konvex, so gilt Gleichheit nur für $X=Y$.

*Beweis.* Seien $X=\sum_ia_i|u_i\rangle\langle u_i|$ und $Y=\sum_jb_j|v_j\rangle\langle v_j|$ mit
Orthonormalbasen aus Eigenvektoren, und sei $c_{ij}=|\langle u_i,v_j\rangle|^2$, also
$\sum_ic_{ij}=\sum_jc_{ij}=1$. Dann ist $\operatorname{Tr}f(X)=\sum_{ij}c_{ij}f(a_i)$,
$\operatorname{Tr}f(Y)=\sum_{ij}c_{ij}f(b_j)$ und
$\operatorname{Tr}Xf'(Y)=\sum_{ij}c_{ij}a_if'(b_j)$, $\operatorname{Tr}Yf'(Y)=\sum_{ij}c_{ij}b_jf'(b_j)$.
Folglich ist die linke Seite gleich $\sum_{ij}c_{ij}[f(a_i)-f(b_j)-(a_i-b_j)f'(b_j)]$, nach
Lemma 3.4 eine Summe nichtnegativer Terme. Ist $f$ strikt konvex und die
Summe null, so gilt $a_i=b_j$, sobald $c_{ij}>0$, also
$Xv_j=\sum_ia_i\langle u_i,v_j\rangle u_i=b_jv_j$ für jedes $j$, und $X=Y$. $\square$

**Korollar 3.6.** $S(\rho\|\sigma)\ge0$, mit Gleichheit genau dann, wenn $\rho=\sigma$.
Insbesondere ist $S(\rho)\le\log D$, mit Gleichheit genau dann, wenn $\rho=1/D$.

*Beweis.* Man wende Satz 3.5 auf $f(t)=t\log t$, $f(0)=0$, an – eine Funktion, die auf
$I=[0,\infty)$ strikt konvex ist ($f''(t)=1/t>0$ auf $(0,\infty)$, und $f$ ist stetig in $0$)
und auf $(0,\infty)\supseteq\operatorname{spec}\sigma$ differenzierbar mit $f'(t)=\log t+1$ –,
mit $X=\rho$, $Y=\sigma$: Die linke Seite ist
$\operatorname{Tr}\rho\log\rho-\operatorname{Tr}\sigma\log\sigma-\operatorname{Tr}(\rho-\sigma)(\log\sigma+1)
=S(\rho\|\sigma)$, da $\operatorname{Tr}\rho=\operatorname{Tr}\sigma=1$. Für die zweite Behauptung
nehme man $\sigma=1/D$: $0\le S(\rho\|1/D)=-S(\rho)+\log D$. $\square$

**Satz 3.7 (Kettenregel; Monotonie unter der bedingten Erwartung auf das Zentrum).** Seien $\rho,\sigma$ Zustände mit
$\sigma>0$, $p_c=\operatorname{Tr}\rho_c$, $q_c=\operatorname{Tr}\sigma_c>0$, $\hat\sigma_c=\sigma_c/q_c$ und,
für $p_c>0$, $\hat\rho_c=\rho_c/p_c$. Dann gilt, mit Summen über die $c$ mit $p_c>0$ (die
übrigen Sektoren tragen wegen $0\log0=0$ den Wert $0$ bei),
$$
S(\rho\|\sigma)=\sum_cp_c\log\frac{p_c}{q_c}+\sum_cp_c\,S(\hat\rho_c\|\hat\sigma_c)
\ \ge\ \sum_cp_c\log\frac{p_c}{q_c}=S(E\rho\|E\sigma),
$$
wobei $E$ die bedingte Erwartung auf das Zentrum ist (Proposition 2.8).

*Beweis.* Für $p_c>0$ ist $\operatorname{Tr}\rho_c\log\rho_c=p_c\log p_c+p_c\operatorname{Tr}\hat\rho_c\log\hat\rho_c$
(Eigenwerte $p_c\mu$, wie in Satz 3.3) und $\log\sigma_c=(\log q_c)1+\log\hat\sigma_c$, also
$\operatorname{Tr}\rho_c\log\sigma_c=p_c\log q_c+p_c\operatorname{Tr}\hat\rho_c\log\hat\sigma_c$ und
$\operatorname{Tr}\rho_c\log\rho_c-\operatorname{Tr}\rho_c\log\sigma_c=p_c\log(p_c/q_c)+p_cS(\hat\rho_c\|\hat\sigma_c)$.
Ist $p_c=0$, so ist $\rho_c=0$, und der Sektor trägt $0$ bei.
Man summiere über $c$; die zweite Summe ist nach Korollar 3.6, angewandt in jedem Block,
$\ge0$. Schließlich ist
$E\rho=\bigoplus_c(p_c/k_c)1_{k_c}$, $E\sigma=\bigoplus_c(q_c/k_c)1_{k_c}$, und
$\operatorname{Tr}E\rho(\log E\rho-\log E\sigma)=\sum_ck_c\frac{p_c}{k_c}\log\frac{p_c}{q_c}$. $\square$

Das Vergessen der Quanteninformation innerhalb der Sektoren kann die Unterscheidbarkeit
zweier Zustände also nur verringern; was bleibt, ist die klassische
Kullback-Leibler-Divergenz der Sektorverteilungen.

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

## 3.6 Der Michelson-Kontrast

Für $A\ge0$, $A\neq0$, ist der **Michelson-Kontrast**
$\Delta(A)=(\lambda_{\max}-\lambda_{\min})/(\lambda_{\max}+\lambda_{\min})\in[0,1]$, wobei die
extremalen Eigenwerte über alle Sektoren genommen werden (und $\Delta(0)=0$). Er misst, wie
weit $A$ von einem Vielfachen der Identität entfernt ist: $\Delta(A)=0$ genau dann, wenn $A=\lambda1$.

**Proposition 3.8.** Ist $\Phi:\mathcal M\to\mathcal M$ linear, positiv und unital
($\Phi(1)=1$), so gilt $\Delta(\Phi(A))\le\Delta(A)$ für alle $A\ge0$.

*Beweis.* Es ist $\lambda_{\min}1\le A\le\lambda_{\max}1$, und Positivität und Unitalität
liefern $\lambda_{\min}1\le\Phi(A)\le\lambda_{\max}1$. Folglich ist
$\lambda_{\min}(\Phi(A))\ge\lambda_{\min}$ und $\lambda_{\max}(\Phi(A))\le\lambda_{\max}$. Ist
$\Phi(A)=0$, so ist die Behauptung trivial. Andernfalls ist mit $r=\lambda_{\min}/\lambda_{\max}$
die Größe $\Delta=(1-r)/(1+r)$ fallend in $r$, und $r(\Phi(A))\ge r(A)$. $\square$

```python
A_pos = alg.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(4, d), batch_size=4, force_positive=True)
Phi = channels.random_mixed_unitary_channel(alg, 3, batch_size=4)       # unital
assert (Phi(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
assert (E(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
```

## 3.7 Eine plausible Ungleichung, die falsch ist

Gilt für positive $X,Y$ in $M_k$ und ein unitäres $U$ die Ungleichung
$|\operatorname{Tr}(XUY)|\le\operatorname{Tr}(XY)$? Sie gilt für $U=1$ und sieht wie eine
Aussage vom Cauchy-Schwarz-Typ aus. Sie ist jedoch falsch, sobald $X$ und $Y$ nicht kommutieren.

**Proposition 3.9.** Seien $X,Y\ge0$ in $M_k$. Dann gilt $\max_U|\operatorname{Tr}(XUY)|=\|YX\|_1\ge\operatorname{Tr}(XY)$, wobei das Maximum
über die unitären $U\in M_k$ läuft, mit Gleichheit genau dann, wenn $XY=YX$.

*Beweis.* $\operatorname{Tr}(XUY)=\operatorname{Tr}(YX\,U)$, also ist das Maximum nach Satz 1.9
gleich $\|YX\|_1$, und $\|YX\|_1\ge|\operatorname{Tr}(YX)|=\operatorname{Tr}(XY)$, da
$\operatorname{Tr}(XY)=\operatorname{Tr}(X^{1/2}YX^{1/2})\ge0$, weil
$X^{1/2}YX^{1/2}=(Y^{1/2}X^{1/2})^*(Y^{1/2}X^{1/2})\ge0$ (Korollar 1.7). Ist $XY=YX$, so
kommutiert $X^{1/2}$ mit $Y$ (Satz 1.5(d)), also ist $YX=X^{1/2}YX^{1/2}\ge0$ und
$\|YX\|_1=\operatorname{Tr}(YX)$ (für $P\ge0$ ist $\|P\|_1=\operatorname{Tr}P$, Abschnitt 1.4).
Umgekehrt gelte $\|YX\|_1=\operatorname{Tr}(YX)$. Ist $|YX|=0$, so ist $YX=0=(YX)^*=XY$.
Andernfalls schreibe man $YX=W|YX|$. Nach Satz 1.9 bedeutet $\operatorname{Tr}(YX\cdot1)=\|YX\|_1$
Gleichheit im Cauchy-Schwarz-Schritt für $U=1$:
$|YX|^{1/2}$ und $W|YX|^{1/2}$ sind proportional mit einem nichtnegativen Faktor, der nach
Vergleich der Normen gleich $1$ ist; also ist $W|YX|^{1/2}=|YX|^{1/2}$ und $YX=W|YX|=|YX|\ge0$.
Eine positive Matrix ist selbstadjungiert, also $YX=(YX)^*=XY$. $\square$

Die Bibliothek findet die Maximalstelle mithilfe der Polarzerlegung. Die Begleitarbeit zeigt
numerisch, dass Monte-Carlo-Stichproben von $U$ dahinter zurückbleiben, und zwar umso
deutlicher, je größer $k$ wird – eines ihrer zentralen Themen; der folgende Code zeigt dies
für ein Paar. Der Grund: $U=1$ ist bereits nahezu optimal, während für ein Haar-zufälliges $U$
in hoher Dimension $|\operatorname{Tr}(XUY)|$ typischerweise viel kleiner als $\operatorname{Tr}(XY)$ ist
(eine zufällige unitäre Matrix verwürfelt die Phasen); daher bleiben alle Stichproben weit unter
dem exakten Supremum.

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

## 3.8 Womit man vergleichen sollte

QuTiP (`entropy_vn`, `entropy_relative`) und die Pakete für Quanteninformation (z. B.
`qutip`, `toqito`) berechnen Entropien und relative Entropien von Dichtematrizen auf einem
einzelnen Hilbertraum; auf $\bigoplus_c\mathbb C^{k_c}$ stimmen sie mit den hiesigen
Konventionen der ungewichteten Spur überein. Was sie nicht bereitstellen, ist die Familie
von Spuren aus Satz 3.1 und die Umrechnungen aus Satz 3.2, die notwendig werden, sobald
eine Algebra mit Vielfachheiten wirkt oder die Sektoren unterschiedlich gewichtet werden
müssen.

## Aufgaben

1. Man zeige, dass die Dichte des Spurzustands $\tau_{\mathrm{vN}}$ (gleiche Sektorgewichte)
   bezüglich $\tau_{\mathrm{vN}}$ gleich $1$ ist, und berechne seine ungewichtete Dichte und
   seine Entropie $S$.
2. Sei $\operatorname{tr}_w$ treu. Man zeige $S_w(\rho)\le\log\operatorname{tr}_w(1)$, mit Gleichheit
   genau dann, wenn $\rho_w=1/\operatorname{tr}_w(1)$.
3. Man zeige, dass der Messkanal $\rho\mapsto\sum_cP_c\rho P_c$ auf $M_D$ (nicht auf
   $\mathcal M$!) Zustände sehr wohl stört, und berechne $S(\rho)$ vorher und nachher für
   eine Superposition zweier Sektoren.

<details>
<summary>Lösungen</summary>

1. $\tau_{\mathrm{vN}}(A)=\tau_{\mathrm{vN}}(1\cdot A)$, also $\rho_{\tau}=1$ wegen der
   Eindeutigkeit. Nach Satz 3.2 ist die ungewichtete Dichte $\bigoplus_c1_{k_c}/(Ck_c)$, mit
   $p_c=1/C$ und $\sigma_c=1/k_c$, also $S=\log C+\frac1C\sum_c\log k_c$ nach Satz 3.3.

```python
tau = states.tracial_state(alg)
assert (tau.density('tau_vN') - one.matrix).abs().max() < 1e-12
assert abs(tau.von_neumann_entropy().item() - (math.log(3) + sum(math.log(k) for k in (2, 3, 5)) / 3)) < 1e-12
```

2. Summiert wird über die $c$ mit $p_c>0$. Dann gilt
   $S_w(\rho)=S(\rho)+\sum_cp_c\log w_c=-\sum_cp_c\log p_c+\sum_cp_cS(\sigma_c)+\sum_cp_c\log w_c
   \le\sum_cp_c\log\frac{k_cw_c}{p_c}\le\log\sum_{p_c>0}k_cw_c\le\log\operatorname{tr}_w(1)$, wobei
   $S(\sigma_c)\le\log k_c$ (Korollar 3.6) und die endliche Jensensche Ungleichung
   $\sum_cp_c\log x_c\le\log\sum_cp_cx_c$ für den konkaven $\log$ (die aus der Definition der
   Konkavität durch Induktion nach der Anzahl der Terme folgt) mit $x_c=k_cw_c/p_c$
   verwendet wurden; sie ist strikt, außer wenn alle $x_c$ gleich sind, da $\log$ strikt
   konkav ist. Gleichheit überall erzwingt, dass jedes $p_c>0$ ist, $\sigma_c=1/k_c$ und
   $k_cw_c/p_c$ konstant, d. h. $p_c=k_cw_c/\operatorname{tr}_w(1)$, d. h.
   $(\rho_w)_c=\rho_c/w_c=1_{k_c}/\operatorname{tr}_w(1)$.
   Umgekehrt gilt für $\rho_w=1/\operatorname{tr}_w(1)$
   $S_w(\rho)=-\operatorname{tr}_w(1)\cdot\frac1{\operatorname{tr}_w(1)}\log\frac1{\operatorname{tr}_w(1)}=\log\operatorname{tr}_w(1)$.

3. In $M_D$ hat eine reine Superposition $\psi=(\psi_1\oplus\psi_2)/\sqrt2$ von
   Einheitsvektoren in zwei Sektoren $S=0$; nach der Messung ist der Zustand
   $\frac12|\psi_1\rangle\langle\psi_1|\oplus\frac12|\psi_2\rangle\langle\psi_2|$ mit $S=\log2$.
   Eingeschränkt auf $\mathcal M$ stimmt der Vektorzustand von $\psi$ jedoch mit der Mischung
   überein:
   $\langle\psi,A\psi\rangle=\tfrac12\langle\psi_1,A_1\psi_1\rangle+\tfrac12\langle\psi_2,A_2\psi_2\rangle$
   für alle $A\in\mathcal M$. Keine Observable in $\mathcal M$ erkennt die Kohärenz; genau das
   ist der Gehalt einer Superauswahlregel.

</details>

**Weiter:** [Kapitel 4](04_symmetries.md) behandelt Symmetrien und Darstellungen: wie eine
Gruppenwirkung auf einem physikalischen System die Sektoren hervorbringt und wie die
Bibliothek die zugehörigen Algebren für Spinsysteme und Qubit-Register konstruiert.
