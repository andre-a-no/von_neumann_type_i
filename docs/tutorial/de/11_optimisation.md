[English](../11_optimisation.md) · [Español](../es/11_optimisation.md) · [Français](../fr/11_optimisation.md) · **Deutsch** · [中文](../zh/11_optimisation.md) · [日本語](../ja/11_optimisation.md) · [Русский](../ru/11_optimisation.md)

# 11. Optimierung und Spurungleichungen

Eine Spurungleichung ist eine Aussage über die Extremwerte eines Funktionals von Operatoren. Abschnitt 3.7 zeigte
eine, die plausibel und falsch ist, und einen Stichprobentest, der die Verletzung übersieht. Dieses Kapitel erklärt,
warum Stichproben in hoher Dimension versagen, und ersetzt sie durch Optimierung über exakt parametrisierte
Nebenbedingungsmengen; wo die Extrema in geschlossener Form bekannt sind, beweisen wir sie, sodass der Optimierer an
einem Satz geprüft wird.

**Sie lernen**, wie die Bibliothek unitäre Matrizen und Spektren mit vorgegebenem Kontrast parametrisiert und warum
jeder zulässige Operator erreicht wird; die von Neumannsche Spurungleichung und die exakten Extrema, die sie liefert;
eine Spurungleichung mit ihrer Gleichheitsbedingung; warum Haar-Stichproben weit von den Extrema entfernt bleiben;
und was eine Multistart-Suche beweist.

Im Folgenden sei $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ mit $D=\sum_ck_c\ge2$, $\operatorname{Tr}$ die
ungewichtete Spur (Abschnitt 3.1) und $\|A\|_1=\sum_c\|A_c\|_1$ die Spurnorm von $A$ als blockdiagonale Matrix in
$M_D(\mathbb C)$. Die Singulärwerte von $A\in M_k(\mathbb C)$ sind $s_1(A)\ge\dots\ge s_k(A)$; $x^\downarrow$ und
$x^\uparrow$ bezeichnen die fallende und die steigende Umordnung eines reellen Vektors $x$. $\Delta$ ist der
Michelson-Kontrast (Abschnitt 3.6), und für $\delta\in[0,1)$ sei
$$
\ell(\delta)=\frac{1-\delta}{1+\delta},\qquad 1-\ell(\delta)=\frac{2\delta}{1+\delta};
$$
ein positiver Operator mit $\lambda_{\max}=1$ hat den Kontrast $\ell(\lambda_{\min})$, also genau dann den Kontrast
$\delta$, wenn $\lambda_{\min}=\ell(\delta)$, da $\ell$ eine streng fallende Involution von $[0,1]$ ist. Die Beispiele
sind die drei Funktionale der Begleitarbeit (Abschnitt *From sampling to optimisation*):

| | Funktional | Variablen |
|---|---|---|
| Exp. 1 | $z=\lvert\operatorname{Tr}(XUY)\rvert-\operatorname{Tr}(XY)$ | $X,Y\ge0$, $U\in\mathcal M$ unitär |
| Exp. 2 | $z=\operatorname{Tr}\lvert XY\rvert-\operatorname{Tr}(\lvert X\rvert\,\lvert Y\rvert)$ | $X=UX_0$, $Y=VY_0$, $X_0,Y_0\ge0$, $U,V$ unitär |
| Exp. 3 | $z=\operatorname{Tr}(Y\lvert X\rvert Y)-\operatorname{Tr}\lvert YXY\rvert$ | $X=X^*$, $Y\ge0$ |

wobei $|X|$, $|Y|$ auf $\lambda_{\max}=1$ normiert sind und ihre Kontraste als Parameter dienen. Die
**Einhüllenden** sind $\overline z(\delta)=\sup\{z:\Delta(|X|)=\delta\}$ und
$\underline z(\delta)=\inf\{z:\Delta(|X|)=\delta\}$, wobei alle übrigen Variablen frei sind.

## 11.1 Nebenbedingungsmengen und ihre Parametrisierungen

**Definition.** Für $\delta\in[0,1)$ und $\ell=\ell(\delta)$ ist die **Kontrastschale**
$$
\mathcal P_\delta=\{X\in\mathcal M:X\ge0,\ \lambda_{\max}(X)=1,\ \Delta(X)=\delta\}
=\{X=X^*:\operatorname{spec}X\subseteq[\ell,1],\ 1,\ell\in\operatorname{spec}X\}.
$$
Die **Plätze** von $\mathcal M$ sind die $D$ Paare $(c,i)$, $i\le k_c$; für $\lambda\in\mathbb R^D$ sei
$\operatorname{diag}(\lambda)=\bigoplus_c\operatorname{diag}(\lambda_{c,1},\dots,\lambda_{c,k_c})$.

`optimize.PositiveParam` mit vorgegebenem Kontrast $\delta$ zieht pro Start zwei verschiedene Plätze
$s_\pm=(c_\pm,i_\pm)$ und eine Haar-verteilte unitäre Matrix $U_0\in\mathcal M$ und bildet $\theta\in\mathbb R^D$ und
$G=\bigoplus_cG_c$ ab auf
$$
X(\theta,G)=U\operatorname{diag}(\lambda)U^*,\qquad U=U_0\,e^{G-G^*},\qquad
\lambda_{s_+}=1,\quad\lambda_{s_-}=\ell,\quad\lambda_s=\ell+(1-\ell)\tfrac{1+\sin\theta_s}2\ \text{otherwise}.
$$
`UnitaryParam` ist die Abbildung $G\mapsto U_0e^{G-G^*}$ allein. Anders als eine Sigmoidfunktion erreicht
$(1+\sin\theta)/2$ die Werte $0$ und $1$ bei endlichem $\theta$; das ist wichtig, weil extremale Spektren oft bei
$\ell$ und $1$ liegen. Ein freier Kontrast in $[a,b]$ ist ein weiterer Parameter, $\delta=a+(b-a)(1+\sin\eta)/2$.

**Lemma 11.1.** Für jedes unitäre $U_0\in\mathcal M$ bildet $G\mapsto U_0e^{G-G^*}$ surjektiv auf die unitäre Gruppe
von $\mathcal M$ ab.

*Beweis.* Sei $U$ unitär. Nach Satz 1.3 (blockweise) und Korollar 1.4(2) ist
$U_0^*U=\bigoplus_cW_c\operatorname{diag}(e^{i\varphi_{c,j}})_jW_c^*$ mit unitären $W_c$ und $\varphi_{c,j}\in\mathbb R$.
Dann erfüllt $S=\bigoplus_cW_c\operatorname{diag}(i\varphi_{c,j})_jW_c^*$ die Gleichungen $S^*=-S$ und $e^S=U_0^*U$ (die
Reihe $\sum_nS^n/n!$ wird gliedweise in der Eigenbasis berechnet), und $G=S/2$ liefert $G-G^*=S$. $\square$

**Proposition 11.2 (die Parametrisierung ist surjektiv).** Für feste $\delta$, $s_\pm$ und $U_0$ ist das Bild von
$(\theta,G)\mapsto X(\theta,G)$ gleich $\{X\in\mathcal P_\delta:1\in\operatorname{spec}X_{c_+},\ \ell\in\operatorname{spec}X_{c_-}\}$,
und $\mathcal P_\delta$ ist die Vereinigung dieser Bilder über alle Paare verschiedener Plätze.

*Beweis.* Das Bild liegt in dieser Menge: $X$ hat Eigenwerte $\lambda_s\in[\ell,1]$, mit $\lambda_{s_\pm}$ in den
Blöcken $c_\pm$. Umgekehrt diagonalisiere man $X_c=W_c\operatorname{diag}(\mu_c)W_c^*$ und permutiere jedes $\mu_c$
(man ersetze $W_c$ durch $W_c\Pi_c$ mit einer Permutationsmatrix $\Pi_c$) so, dass $\mu_{s_+}=1$ und $\mu_{s_-}=\ell$.
Das ist möglich: Ist $c_+=c_-$ und $\delta>0$, so sind $1\ne\ell$ zwei Eigenwerte desselben Blocks und besetzen
verschiedene Positionen; ist $\delta=0$, so sind alle Eigenwerte gleich $1$. Ist $\ell<1$, so setze man
$\theta_s=\arcsin\big(2(\mu_s-\ell)/(1-\ell)-1\big)$ für die übrigen Plätze, sodass $\lambda=\mu$ (für $\ell=1$ ist
$\theta$ beliebig). Lemma 11.1 liefert $G$ mit $U_0e^{G-G^*}=\bigoplus_cW_c$. *Vereinigung:* Zu $X\in\mathcal P_\delta$
wähle man Blöcke $c_\pm$, die die Eigenwerte $1$ bzw. $\ell$ enthalten, und darin verschiedene Plätze; ist $c_+=c_-$
und $\delta>0$, so hat der Block zwei verschiedene Eigenwerte, also $k_{c_+}\ge2$; ist $\delta=0$, so ist $X=1$, und
zwei beliebige Plätze genügen. $\square$

Zufällige Plätze erlauben es einer Multistart-Suche also, die ganze Schale abzudecken, und die Nebenbedingung gilt in
jeder Iteration exakt. `SelfAdjointParam` verwendet $X=U\operatorname{diag}(\sigma\lambda)U^*$ mit pro Start gezogenen
Vorzeichen $\sigma\in\{\pm1\}^D$ (die nur durch die Vorzeichenwechsel aus Abschnitt 11.5 verändert werden), sodass
$|X|=U\operatorname{diag}(\lambda)U^*$, und (mit demselben Beweis) ist die Vereinigung über Vorzeichen und Plätze gleich
$\{X=X^*:|X|\in\mathcal P_\delta\}$. Es kommt nur auf $S=G-G^*$ an. Die Abbildung $S\mapsto U_0e^S$ auf
schiefadjungierten $S$ hat in $S=0$ das Differential $S\mapsto U_0S$, eine Bijektion auf den Tangentialraum
$U_0\cdot\{S^*=-S\}$ der unitären Gruppe; sie ist daher nahe $S=0$ ein lokaler Diffeomorphismus. Fern von $0$ ist sie
nicht injektiv, und ihr Differential ist singulär, wo sich zwei Eigenwerte von $S$ um ein von null verschiedenes
Vielfaches von $2\pi i$ unterscheiden. Deshalb **zentriert** `extremize` nach jeder Runde **neu**
($U_0\leftarrow U_0e^{G-G^*}$, $G\leftarrow0$). Die Cayley-Karte $S\mapsto U_0C(-S/2)$ (Aufgabe 1) und die polare
Retraktion $S\mapsto$ unitärer Faktor von $U_0(1+S)$ (Satz 1.8) haben in $0$ dasselbe Differential.

Der Code prüft die Nebenbedingung und, konstruktiv, Proposition 11.2: Die Parameter aus dem Beweis reproduzieren
eine Zielmatrix.

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

> **Physikalische Bedeutung.** Der Michelson-Kontrast ist die Sichtbarkeit $(I_{\max}-I_{\min})/(I_{\max}+I_{\min})$
> von Interferenzstreifen. Für $0\le X\le1$ mit $\lambda_{\max}=1$ – ein Filter, ein Abschwächer, eine
> Detektoreffizienz – legt er den Dynamikbereich $\lambda_{\min}/\lambda_{\max}$ fest, den Kehrwert einer Konditionszahl.

## 11.2 Die von Neumannsche Spurungleichung

Eine Matrix $(c_{ij})$ mit $c_{ij}\ge0$, deren Zeilen- und Spaltensummen alle gleich $1$ sind, heißt **doppelt
stochastisch**; für eine unitäre Matrix $G$ sind $(|G_{ij}|^2)$ und $(|G_{ji}|^2)$ doppelt stochastisch, da die Zeilen
und Spalten von $G$ Einheitsvektoren sind.

**Lemma 11.3 (Umordnung).** Seien $p_1\ge\dots\ge p_k$, $q_1\ge\dots\ge q_k$ reell und $(c_{ij})$ doppelt
stochastisch. Dann gilt $\sum_ip_iq_{k+1-i}\le\sum_{i,j}p_iq_jc_{ij}\le\sum_ip_iq_i$.

*Beweis.* *Obere Schranke.* Für $m\le k$ sei $r_j=\sum_{i\le m}c_{ij}\in[0,1]$, also $\sum_jr_j=m$. Da $q_j\ge q_m$
für $j\le m$ und $q_j\le q_m$ für $j>m$, gilt
$$
\sum_{j\le m}q_j-\sum_jq_jr_j=\sum_{j\le m}q_j(1-r_j)-\sum_{j>m}q_jr_j\ge q_m\Big(m-\sum_jr_j\Big)=0 .
$$
Mit $x_i=\sum_jq_jc_{ij}$ besagt dies $\sum_{i\le m}x_i\le\sum_{j\le m}q_j$, mit Gleichheit für $m=k$. Mit abelscher
partieller Summation und $p_m-p_{m+1}\ge0$ folgt
$$
\sum_ip_ix_i=p_k\sum_{i\le k}x_i+\sum_{m<k}(p_m-p_{m+1})\sum_{i\le m}x_i
\le p_k\sum_{j\le k}q_j+\sum_{m<k}(p_m-p_{m+1})\sum_{j\le m}q_j=\sum_ip_iq_i .
$$
*Untere Schranke.* Man wende die obere Schranke auf $p$, die fallende Folge $q'_j=-q_{k+1-j}$ und die doppelt
stochastische Matrix $c'_{ij}=c_{i,k+1-j}$ an:
$-\sum_{ij}p_iq_jc_{ij}=\sum_{ij}p_iq'_jc'_{ij}\le\sum_ip_iq'_i=-\sum_ip_iq_{k+1-i}$. $\square$

Jedes $A\in M_k(\mathbb C)$ besitzt eine **Singulärwertzerlegung** $A=W\Sigma V^*$ mit unitären $W,V$ und
$\Sigma=\operatorname{diag}(s_i(A))$: Man schreibe $A=W'|A|$ (Satz 1.8), $|A|=V\Sigma V^*$ (Satz 1.3) und setze
$W=W'V$.

**Satz 11.4 (von Neumannsche Spurungleichung).** Für $A,B\in M_k(\mathbb C)$ gilt
$|\operatorname{Tr}(AB)|\le\sum_is_i(A)\,s_i(B)$.

*Beweis.* Seien $A=W_A\Sigma_AV_A^*$, $B=W_B\Sigma_BV_B^*$, $a_i=s_i(A)$, $b_j=s_j(B)$ sowie $G=V_A^*W_B$,
$H=V_B^*W_A$. Wegen der Zyklizität ist
$\operatorname{Tr}(AB)=\operatorname{Tr}(\Sigma_AG\Sigma_BH)=\sum_{i,j}a_iG_{ij}b_jH_{ji}$, also wegen
$|g||h|\le\frac12(|g|^2+|h|^2)$ und $a_ib_j\ge0$
$$
|\operatorname{Tr}(AB)|\le\sum_{i,j}a_ib_j\,d_{ij},\qquad d_{ij}=\tfrac12\big(|G_{ij}|^2+|H_{ji}|^2\big).
$$
$(d_{ij})$ ist doppelt stochastisch, und Lemma 11.3 liefert $\sum_{ij}a_ib_jd_{ij}\le\sum_ia_ib_i$. $\square$

**Korollar 11.5.** (a) Für selbstadjungierte $P,Q\in M_k(\mathbb C)$ mit Eigenwerten $p,q$ gilt
$\sum_ip^\downarrow_iq^\uparrow_i\le\operatorname{Tr}(PQ)\le\sum_ip^\downarrow_iq^\downarrow_i$.
(b) Für $P,Q\ge0$ gilt $\operatorname{Tr}(PQ)\le\|PQ\|_1\le\sum_ip^\downarrow_iq^\downarrow_i$.

*Beweis.* (a) Wie im Beweis von Satz 3.5 ist $\operatorname{Tr}(PQ)=\sum_{ij}p_iq_jc_{ij}$ mit der doppelt
stochastischen Matrix $c_{ij}=|\langle u_i,v_j\rangle|^2$ für fallend geordnete Eigenbasen; man wende Lemma 11.3 an
(dessen untere Schranke ist $\sum_ip^\downarrow_iq^\downarrow_{k+1-i}=\sum_ip^\downarrow_iq^\uparrow_i$). (b) Die erste
Ungleichung ist Satz 1.9. Nach Satz 1.9 und Satz 11.4 ist
$\|PQ\|_1=\max_U|\operatorname{Tr}(P\,QU)|\le\sum_is_i(P)s_i(QU)$, und $s_i(QU)=q^\downarrow_i$, weil
$|QU|=(U^*Q^2U)^{1/2}=U^*QU$ (Satz 1.6). $\square$

Die Verletzung in Proposition 3.9 beträgt also höchstens $\sum_ix^\downarrow_i(y^\downarrow_i-y^\uparrow_i)$. Für Exp. 2
ist diese Zahl die exakte Antwort:

**Proposition 11.6 (Exp. 2 bei festen Spektren).** In Exp. 2 ist $|X|=X_0$, $|Y|=Y_0$ und
$z=\|X_0VY_0\|_1-\operatorname{Tr}(X_0Y_0)$. Sind die Spektren $x_c$ von $(X_0)_c$ und $y_c$ von $(Y_0)_c$ fest und
variieren die Eigenbasen und $V$, so gilt
$$
\sup z=-\inf z=\sum_c\sum_ix^\downarrow_{c,i}\big(y^\downarrow_{c,i}-y^\uparrow_{c,i}\big),
$$
und beide Werte werden angenommen.

*Beweis.* $X^*X=X_0U^*UX_0=X_0^2$, also $|X|=X_0$ (Satz 1.6); ebenso $|Y|=Y_0$, und $|UA|=|A|$ liefert
$\operatorname{Tr}|XY|=\|X_0VY_0\|_1$. Alle Variablen wirken blockweise, und $z$ ist eine Summe über die Blöcke; sei
also $C=1$. Man setze $Y'=VY_0V^*\ge0$ (Spektrum $y$); wegen $|AV|=V^*|A|V$ ist
$\|X_0VY_0\|_1=\|X_0Y'V\|_1=\|X_0Y'\|_1$. Korollar 11.5 für die Paare $(X_0,Y')$ und $(X_0,Y_0)$ liefert zusammen
mit $\|X_0Y'\|_1\ge\operatorname{Tr}(X_0Y')$
$$
\sum_ix^\downarrow_iy^\uparrow_i-\sum_ix^\downarrow_iy^\downarrow_i\le\|X_0Y'\|_1-\operatorname{Tr}(X_0Y_0)\le\sum_ix^\downarrow_iy^\downarrow_i-\sum_ix^\downarrow_iy^\uparrow_i .
$$
Sei $J$ die Umkehrpermutation, $Je_i=e_{k+1-i}$; dann ist $J\operatorname{diag}(w)=\operatorname{diag}(w_k,\dots,w_1)J$, und
die Spurnorm von $\operatorname{diag}(w)J$ ist $\sum_i|w_i|$. Mit $X_0=\operatorname{diag}(x^\downarrow)$ und $V=J$ gilt:
Für $Y_0=\operatorname{diag}(y^\uparrow)$ ist $X_0VY_0=\operatorname{diag}(x^\downarrow_iy^\downarrow_i)J$ und
$\operatorname{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\uparrow_i$ (die obere Schranke); für $Y_0=\operatorname{diag}(y^\downarrow)$
ist $X_0VY_0=\operatorname{diag}(x^\downarrow_iy^\uparrow_i)J$ und $\operatorname{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\downarrow_i$
(die untere Schranke). $\square$

**Korollar 11.7 (Einhüllende von Exp. 2).** Seien $\ell_X,\ell_Y\in(0,1]$ und $k_c\ge2$ für mindestens ein $c$. Über alle
$X_0,Y_0\ge0$ mit $\lambda_{\max}=1$, $\lambda_{\min}(X_0)=\ell_X$, $\lambda_{\min}(Y_0)=\ell_Y$ und alle unitären
Matrizen gilt
$$
\sup z=-\inf z=\sum_c\Big\lfloor\frac{k_c}2\Big\rfloor(1-\ell_X)(1-\ell_Y).
$$

*Beweis.* In einem Block schreibe man $x_i=x^\downarrow_{c,i}$, $y_i=y^\downarrow_{c,i}$, $k=k_c$, sodass
$y^\uparrow_{c,i}=y_{k+1-i}$. Paart man $i$ mit $k+1-i$, so ist
$x_i(y_i-y_{k+1-i})+x_{k+1-i}(y_{k+1-i}-y_i)=(x_i-x_{k+1-i})(y_i-y_{k+1-i})$, und für ungerades $k$ trägt der mittlere
Index $x_i(y_i-y_i)=0$ bei. Folglich gilt
$$
\sum_ix_i(y_i-y_{k+1-i})=\sum_{i\le\lfloor k/2\rfloor}(x_i-x_{k+1-i})(y_i-y_{k+1-i})\le\Big\lfloor\frac k2\Big\rfloor(1-\ell_X)(1-\ell_Y),
$$
da $x_i\in[\ell_X,1]$, $y_i\in[\ell_Y,1]$. Gleichheit gilt, wenn in jedem Block die $\lfloor k_c/2\rfloor$ größten
Eigenwerte gleich $1$ und die $\lfloor k_c/2\rfloor$ kleinsten gleich $\ell_X$ (bzw. $\ell_Y$) sind, während der
mittlere beliebig ist; diese Spektren haben die globalen Extrema $1$ und $\ell_X$ ($\ell_Y$), weil $k_c\ge2$ für mindestens ein $c$ gilt.
Proposition 11.6 schließt den Beweis ab. $\square$

Hier ist $1-\ell_X=2\delta/(1+\delta)$ für $\Delta(|X|)=\delta$; für ungerades $k_c$ spielt der mittlere Eigenwert keine
Rolle. Der Code prüft Satz 11.4, Korollar 11.5 und die extremalen Konfigurationen aus Proposition 11.6 für $k=3$,
$C=2$.

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

## 11.3 Eine Spurungleichung mit Gleichheitsbedingung

Für $X=X^*$ sei $X_\pm=f_\pm(X)$ mit $f_+(t)=\max(t,0)$, $f_-(t)=\max(-t,0)$. Dann gilt $X_\pm\ge0$, $X=X_+-X_-$,
$|X|=X_++X_-$ und $X_+X_-=0$ (Satz 1.5(c)).

**Lemma 11.8.** Für $A,B\ge0$ in $M_n(\mathbb C)$ gilt $\|A-B\|_1\le\operatorname{Tr}A+\operatorname{Tr}B$, mit Gleichheit
genau dann, wenn $AB=0$.

*Beweis.* Sei $H=A-B$ und $Q=g(H)$ mit $g=1$ auf $(0,\infty)$ und $g=0$ sonst. Wegen $|t|=t(2g(t)-1)$ ist
$|H|=H(2Q-1)$ und
$$
\|H\|_1=\operatorname{Tr}\big(A(2Q-1)\big)-\operatorname{Tr}\big(B(2Q-1)\big)
=\operatorname{Tr}A+\operatorname{Tr}B-2\operatorname{Tr}\big((1-Q)A(1-Q)\big)-2\operatorname{Tr}(QBQ).
$$
Die abgezogenen Spuren sind $\ge0$, und $\operatorname{Tr}(RAR)=\|A^{1/2}R\|_2^2$ verschwindet genau dann, wenn $AR=0$;
Gleichheit gilt also genau dann, wenn $A(1-Q)=0$ und $BQ=0$, und dann ist $AB=AQB=A(BQ)^*=0$. Ist umgekehrt
$AB=0$, so ist $BA=0$, $(A-B)^2=A^2+B^2=(A+B)^2$, und Satz 1.6 liefert $|A-B|=A+B$. $\square$

**Proposition 11.9 (Exp. 3).** Für $X=X^*$ und $Y\ge0$ in $\mathcal M$ gilt
$\operatorname{Tr}|YXY|\le\operatorname{Tr}(Y|X|Y)$, mit Gleichheit genau dann, wenn $YX_+Y^2X_-Y=0$, insbesondere
wenn $XY=YX$.

*Beweis.* In $M_D(\mathbb C)$ sind $A=YX_+Y$ und $B=YX_-Y$ positiv (Korollar 1.7 mit $T=X_\pm^{1/2}Y$), $YXY=A-B$
und $\operatorname{Tr}(Y|X|Y)=\operatorname{Tr}A+\operatorname{Tr}B$; man wende Lemma 11.8 an. Ist $XY=YX$, so kommutiert
$Y$ mit $X_\pm$ (Satz 1.5(d)), und $AB=Y^4X_+X_-=0$. $\square$

Also ist $z\ge0$ in Exp. 3, und $\underline z=0$ (man nehme $Y=1$). Die Gleichheitsbedingung ist schwächer als
Kommutativität: Für $X=\operatorname{diag}(1,-1,0)$ und $Y$ mit den Zeilen $(1,0,a)$, $(0,1,0)$, $(a,0,1)$, $|a|<1$, ist
$AB=(Ye_1)\langle Ye_1,Ye_2\rangle(Ye_2)^*=0$, während $XY\ne YX$ für $a\ne0$.

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

## 11.4 Stichproben oder Optimierung

Wie groß ist $|\operatorname{Tr}(AU)|$ typischerweise, verglichen mit seinem Maximum $\|A\|_1$? Die Bibliothek zieht $U$
als den Faktor $Q$ der QR-Zerlegung $Z=QR$, $R_{ii}>0$, einer **Ginibre-Matrix** $Z$ (unabhängige Einträge mit der
Dichte $\pi^{-1}e^{-|z|^2}$; Mezzadris Rezept), und zwar in jedem Block unabhängig.

**Lemma 11.10.** $Q$ ist Haar-verteilt auf $U(k)$; insbesondere hat $VQ$ für jede feste unitäre Matrix
$V\in M_k(\mathbb C)$ dieselbe Verteilung wie $Q$. Die unitäre Matrix $\bigoplus_cQ_c$ mit unabhängigen Blöcken ist
Haar-verteilt auf der unitären Gruppe von $\mathcal M$.

*Beweis.* Das ist Satz 5.3 (die Linksinvarianz ist der Äquivarianzschritt seines Beweises) zusammen mit
Korollar 5.4. $\square$

(Zur Maßkonzentration und mehr über diese Gruppen siehe E. Meckes, *The Random Matrix Theory of the Classical
Compact Groups*, Cambridge 2019.)

**Proposition 11.11.** Für Haar-verteiltes $U\in U(k)$ und $A\in M_k(\mathbb C)$ gilt $\mathbb E\operatorname{Tr}(AU)=0$
und $\mathbb E|\operatorname{Tr}(AU)|^2=\|A\|_2^2/k$.

*Beweis.* $\mathbb EU_{ij}=0$ nach Korollar 5.5. Proposition 5.7 (mit $C=1$) für $X=|e_j\rangle\langle e_m|$ liefert
die Weingarten-Formel erster Ordnung,
$\mathbb E[U_{ij}\bar U_{lm}]=\mathbb E(UXU^*)_{il}=\delta_{il}\delta_{jm}/k$; für $k\ge2$ folgt sie auch aus der Formel
zweiter Ordnung, Satz 5.10, indem man $i_2=i'_2$, $j_2=j'_2$ setzt und über beide summiert, da
$\sum_{i_2,j_2}|U_{i_2j_2}|^2=k$. Mit $\operatorname{Tr}(AU)=\sum_{ij}A_{ji}U_{ij}$ ist
$\mathbb E|\operatorname{Tr}(AU)|^2=\sum_{ijlm}A_{ji}\bar A_{ml}\mathbb E[U_{ij}\bar U_{lm}]=\frac1k\sum_{ij}|A_{ji}|^2$. $\square$

**Korollar 11.12.** Seien $X,Y\ge0$ in $M_k(\mathbb C)$ mit Spektren in $[a,b]$, $a>0$, und $U$ Haar-verteilt. Dann
gilt
$$
\frac{\big(\mathbb E|\operatorname{Tr}(XUY)|^2\big)^{1/2}}{\max_V|\operatorname{Tr}(XVY)|}=\frac{\|YX\|_2}{\sqrt k\,\|YX\|_1}\le\min\Big(\frac1{\sqrt k},\frac{b^2}{a^2k}\Big),
\qquad
\mathbb P\big(|\operatorname{Tr}(XUY)|\ge\operatorname{Tr}(XY)\big)\le\frac{b^4}{a^4k^2}.
$$

*Beweis.* Die Gleichung ist Proposition 11.11 für $A=YX$ zusammen mit Proposition 3.9.
$\sum_is_i^2\le(\sum_is_i)^2$ liefert $1/\sqrt k$. Da $s_1(A)\le\|Y\|\|X\|\le b^2$ und
$s_k(A)=\min_{\|v\|=1}\|YXv\|\ge a^2$, ist $\|A\|_2^2\le kb^4$ und $\|A\|_1\ge ka^2$. Schließlich ist
$\operatorname{Tr}(XY)=\operatorname{Tr}(X^{1/2}YX^{1/2})\ge a\operatorname{Tr}X\ge ka^2$, und die Markov-Ungleichung für
$|\operatorname{Tr}(XUY)|^2$ beschränkt die Wahrscheinlichkeit durch $(kb^4/k)/(ka^2)^2$. $\square$

$N$ Haar-Stichproben finden eine Verletzung von $|\operatorname{Tr}(XUY)|\le\operatorname{Tr}(XY)$ also höchstens mit
Wahrscheinlichkeit $Nb^4/(a^4k^2)$, obwohl es für jedes nichtkommutierende Paar eine gibt. Tatsächlich ist es noch
drastischer: $f(U)=|\operatorname{Tr}(AU)|$ ist Lipschitz-stetig mit der Konstanten $L=\|A\|_2$ bezüglich
$\|U-V\|_2$, also gilt $\mathbb P(f\ge\mathbb Ef+t)\le e^{-ckt^2/L^2}$ mit einer absoluten Konstanten $c>0$
(Maßkonzentration, ohne Beweis zitiert; Meckes, *loc. cit.*, Kap. 5): hier ein Exponent der Größenordnung
$k^2a^4/b^4$. Der Code prüft Proposition 11.11 und gibt für in $[0.2,1.2]$ gleichverteilte Spektren ($b/a=6$) den
Quotienten aus Korollar 11.12, $\sup z$ und das beste $z$ aus $2000$ Stichproben aus, beide relativ zu
$\operatorname{Tr}(XY)$.

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

Der Quotient fällt wie $1/k$, und ab $k=4$ erreicht die beste Stichprobe nicht einmal $z=0$. Die Extrema sind Mengen
von winzigem Haar-Maß, die ein Gradientenverfahren direkt ansteuert. Für das Paar `Xf, Yf` sind beide Extrema über
$U$ bekannt:

**Lemma 11.13.** Ist $k\ge2$, so gilt $\min_U|\operatorname{Tr}(AU)|=0$ für jedes $A\in M_k(\mathbb C)$.

*Beweis.* Sei $A=W\Sigma V^*$ und $P$ die zyklische Verschiebung, $Pe_i=e_{(i\bmod k)+1}$, deren Diagonale für $k\ge2$
verschwindet. Für $U=VPW^*$ ist $\operatorname{Tr}(AU)=\operatorname{Tr}(\Sigma P)=0$. $\square$

Also ist $\inf_Uz=-\operatorname{Tr}(XY)$ und $\sup_Uz=\|YX\|_1-\operatorname{Tr}(XY)$. Ein einziger Aufruf von `extremize`
findet beide, mit $12$ Starts, die $z$ maximieren, und $12$, die $-z$ maximieren. Stichproben kommen dem Infimum
nahe, das neben den typischen Werten liegt ($|\operatorname{Tr}(XUY)|$ ist typischerweise von der Größenordnung
$\|YX\|_2/\sqrt k$), nicht aber dem Supremum:

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

> **Physikalische Bedeutung.** Maßkonzentration ist die mathematische Gestalt der **Typizität** in der
> statistischen Mechanik: In hoher Dimension liefert fast jeder Zustand und fast jede Zeitentwicklung denselben Wert
> einer makroskopischen Größe. Über eine scharfe Ungleichung entscheiden fein abgestimmte Konfigurationen –
> ausgerichtete Eigenbasen, Zwei-Niveau-Spektren –, die eine zufällige unitäre Matrix (ein Zufallsschaltkreis, eine
> zufällige Streumatrix) praktisch nie hervorbringt.

## 11.5 Multistart-Suche nach Einhüllenden

`optimize.extremize(objective, params, maximize)` führt Adam in `rounds` Runden zu je `steps` Schritten aus
(Lernrate `lr * decay**round`), zentriert die unitären Karten nach jeder Runde neu und gibt **pro Start** den besten
erreichten Wert zurück. Jedes Batch-Element ist ein unabhängiger Start mit eigenen Plätzen, Vorzeichen, $U_0$ und
eigenem Kontrast. Die Vorzeichen von `SelfAdjointParam` sind diskret: Zwischen den Runden probiert `discrete_step`
pro Start acht einzelne Vorzeichenwechsel nacheinander aus und behält jeden Wechsel, der die Zielfunktion verbessert.

Da jede Iterierte zulässig ist (Proposition 11.2), wird jeder zurückgegebene Wert angenommen: Eine Maximierung
liefert eine **untere Schranke** für das Supremum, eine Minimierung eine obere Schranke für das Infimum, und nicht
mehr – eine lokale Suche kann in einem lokalen Extremum stehen bleiben. Die andere Richtung erfordert einen Beweis,
wie Korollar 11.7 und Proposition 11.9, oder eine konvexe Relaxierung (Abschnitt 11.6). Eine berechnete Einhüllende
wird geprüft (i) an geschlossenen Formeln; (ii) an von Hand gebauten Konfigurationen, etwa direkten Summen kleinerer
Lösungen, unter die sie nicht fallen darf; (iii) an strukturierten Stichproben (Zwei-Niveau-Spektren,
Permutationsbasen), die den Extrema in jeder Dimension nahekommen, anders als Haar-Stichproben
(`scripts/bounds_checks.py`).

**Proposition 11.14 (untere Einhüllende von Exp. 1).** Seien alle $k_c\ge2$. Über $X\in\mathcal P_\delta$, $Y\ge0$
mit $\lambda_{\max}(Y)=1$ und unitäre $U\in\mathcal M$ gilt $\ \min z=-\big(D-2\delta/(1+\delta)\big)$.

*Beweis.* Aus $Y\le1$ folgt $X^{1/2}YX^{1/2}\le X$, also $\operatorname{Tr}(XY)\le\operatorname{Tr}X\le(D-1)+\ell(\delta)$,
da ein Eigenwert von $X$ gleich $\ell(\delta)$ ist und die übrigen $\le1$ sind; daher
$z\ge-\operatorname{Tr}(XY)\ge-(D-2\delta/(1+\delta))$. Gleichheit gilt für $X$ mit dem Eigenwert $\ell(\delta)$ einmal
und $1$ sonst, $Y=1$ und $U=\bigoplus_cU_c$ mit $\operatorname{Tr}(X_cU_c)=0$ (Lemma 11.13). $\square$

Der Code prüft Proposition 11.14 auf $M_2\oplus M_2$, Korollar 11.7 auf $M_3\oplus M_3$ (beide Kontraste
vorgegeben, $k=3$ ungerade) und Proposition 11.9; Suprema und Infima teilen sich wie oben einen Aufruf. Das
Schlüsselwort `is_self_adjoint` wählt den Algorithmus für die Spurnorm (Singulärwerte oder, für $YXY$, Eigenwerte),
ohne in jedem Schritt zu testen. Die letzte Prüfung ist vom Typ (ii): $z$ aus Exp. 3 ist additiv über direkte
Summen, und der Kontrast einer direkten Summe ist der größte Kontrast ihrer Summanden (Aufgabe 3), also
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

Die Suche auf $M_2\oplus M_2$ bleibt um etwa $0.2\,\%$ hinter der direkten Summe zurück: Sie ist nicht vollständig
konvergiert, und die beste untere Schranke für $\overline z_{M_2\oplus M_2}(0.5)$ ist die größere Zahl – dies zu
erkennen ist der Zweck von Prüfung (ii). Die Werte für $M_2$ stimmen bis auf $10^{-4}$ mit
$(1-\delta)\big(2d/(1+d)\big)^2=(1-\delta)(1-\ell(d))^2$ überein, wobei $d=0.999$ der größte zugelassene Kontrast von
$Y$ ist: eine geschlossene Form, die die Begleitarbeit numerisch beobachtet (in reellen Algebren), aber nicht
bewiesen hat; der Code behauptet sie nicht. Für größere Blöcke sind die oberen Einhüllenden von Exp. 1 und Exp. 3
unseres Wissens nicht in geschlossener Form bekannt; berechnete Kurven sind dann untere Schranken und Vermutungen,
nicht mehr.

## 11.6 Womit man vergleichen sollte

**Manopt**, **pymanopt** und **geoopt** stellen Matrixmannigfaltigkeiten (Stiefel-Mannigfaltigkeiten, orthogonale
und unitäre Gruppen) mit Retraktionen wie QR, polar oder Cayley bereit sowie Riemannsche Löser (konjugierte
Gradienten, Trust-Region-Verfahren, Riemannsches Adam), die oft effizienter sind als Adam in einer neu zentrierten
Karte. **cvxpy** mit einem SDP-Löser löst konvexe Probleme in Operatorvariablen global optimal, mit einem dualen
Zertifikat bis auf die Lösertoleranz, und semidefinite Relaxierungen liefern rigorose Schranken in der Richtung, die
eine lokale Suche nicht erreichen kann. **NCTSSOS**, **ncpol2sdpa** und **NCSOStools** beschränken Eigenwerte und
Spuren nichtkommutativer Polynome in allen Dimensionen zugleich (Abschnitt 0.3). Die Bibliothek ergänzt exakte
Parametrisierungen der Spektral- und Kontrastnebenbedingungen von $\bigoplus_cM_{k_c}$ (die Extrema können in jedem
Sektor liegen), unitäre Matrizen, die die Blöcke respektieren, und eine einzige gebündelte Rechnung für viele Starts
und Werte der Nebenbedingung, auf CPU oder GPU. Sie liefert Kandidaten für Extremalstellen und untere Schranken für
Suprema; den Rest müssen Beweise oder Relaxierungen beisteuern.

## Aufgaben

1. (**Cayley-Transformation.**) Man zeige für $S^*=-S$, dass $C(S)=(1-S)(1+S)^{-1}$ unitär ist und dass
   $S\mapsto C(S)$ eine Bijektion auf die unitären Matrizen ohne Eigenwert $-1$ ist. Anders als die
   Exponentialabbildung (Lemma 11.1) verfehlt eine Cayley-Karte $U_0C(S)$ also $-U_0$.
2. (**Reelles Haar-Maß.**) Für Haar-verteiltes $O$ auf $O(k)$ (QR einer reellen Ginibre-Matrix) und reelles $A$ zeige
   man $\mathbb E\operatorname{Tr}(AO)^2=\|A\|_2^2/k$ und prüfe dies mit einer reellen Algebra.
3. (**Direkte Summen.**) Man zeige, dass $z$ aus Exp. 3 additiv über direkte Summen ist und dass
   $\Delta(A\oplus A')=\max(\Delta(A),\Delta(A'))$ für $A,A'\ge0$ mit $\lambda_{\max}=1$ gilt. Man folgere die
   Ungleichung für direkte Summen aus Abschnitt 11.5 (mit $\Delta(Y)$ frei in einem festen Intervall $[0,b]$).

<details>
<summary>Lösungen</summary>

1. $S$ hat rein imaginäre Eigenwerte, also sind $1\pm S$ invertierbar; sie kommutieren, und $(1+S)^*=1-S$. Daher
   ist $C^*C=(1-S)^{-1}(1+S)(1-S)(1+S)^{-1}=1$, und $1+C=2(1+S)^{-1}$ ist invertierbar, also
   $-1\notin\operatorname{spec}C$. Umgekehrt ist für eine unitäre Matrix $U$ mit invertierbarem $1+U$ die Matrix
   $S=(1-U)(1+U)^{-1}$ schiefadjungiert: $S^*=(1+U^*)^{-1}(1-U^*)=(U+1)^{-1}U\,U^*(U-1)=-S$. Da $C$ mit $S$
   kommutiert, ist $C(1+S)=1-S$ äquivalent zu $S(1+C)=1-C$, d. h. $S=(1-C)(1+C)^{-1}$: Die beiden Abbildungen sind
   zueinander invers.

```python
Sk = torch.randn(4, 4, dtype=torch.complex128)
Sk = Sk - Sk.mH                                                   # skew-adjoint
I4 = torch.eye(4, dtype=torch.complex128)
Cay = (I4 - Sk) @ torch.linalg.inv(I4 + Sk)
assert torch.allclose(Cay.mH @ Cay, I4, atol=1e-12)
assert torch.allclose((I4 - Cay) @ torch.linalg.inv(I4 + Cay), Sk, atol=1e-10)
assert torch.linalg.eigvals(Cay).add(1).abs().min() > 1e-6        # -1 is never an eigenvalue
```

2. Satz 5.3 und sein Beweis gelten für eine reelle Gaußsche Matrix (Dichte $\propto e^{-\|Z\|_2^2/2}$), die
   orthogonale Gruppe und orthogonales $V$ (die Bemerkung nach Satz 5.3), also hat $VO$ dieselbe Verteilung wie $O$.
   Linksmultiplikation mit $\operatorname{diag}(\pm1)$ liefert $\mathbb E[O_{ij}O_{lm}]=0$ für $i\ne l$, Permutationen
   zeigen, dass $\mathbb E[O_{ij}O_{im}]$ nicht von $i$ abhängt, und $\sum_iO_{ij}O_{im}=\delta_{jm}$; also ist
   $\mathbb E[O_{ij}O_{lm}]=\delta_{il}\delta_{jm}/k$ und $\mathbb E\operatorname{Tr}(AO)^2=\sum_{ij}A_{ji}^2/k$.

```python
Mr = TypeIAlgebra([5], [5], complex_valued=False, precision='double')
Ar = Mr.operator(torch.randn(1, 1, 5, 5, dtype=torch.float64))
tr2 = (Ar @ Mr.random_unitary_operator(20_000)).trace ** 2
assert abs(tr2.mean().item() - Ar.frobenius_norm().item() ** 2 / 5) < 5 * tr2.std().item() / math.sqrt(20_000)
```

3. Spuren und Spurnormen blockdiagonaler Operatoren sind Summen über die Blöcke, und $|X\oplus X'|=|X|\oplus|X'|$;
   $\lambda_{\min}(A\oplus A')=\min(\lambda_{\min}(A),\lambda_{\min}(A'))$, und $\Delta$ fällt mit wachsendem
   $\lambda_{\min}$. Ist also $(X,Y)$ zulässig bei $\delta$ und $(X',Y')$ bei $\delta'\le\delta$, so ist
   $(X\oplus X',Y\oplus Y')$ zulässig bei $\delta$, mit $z=z(X,Y)+z(X',Y')$; man gehe zum Supremum über beide Paare und
   über $\delta'\le\delta$ über.

</details>

**Weiter:** [Kapitel 12](12_equivariant_learning.md) behandelt äquivariantes maschinelles Lernen: Schichten, die mit einer
Symmetrie kommutieren, sind Elemente einer Kommutante (Kapitel 2), aufgebaut aus denselben gebündelten,
differenzierbaren Operatoren.
