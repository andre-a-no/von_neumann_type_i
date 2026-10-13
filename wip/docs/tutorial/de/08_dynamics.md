[English](../08_dynamics.md) · [Español](../es/08_dynamics.md) · [Français](../fr/08_dynamics.md) · **Deutsch** · [中文](../zh/08_dynamics.md) · [日本語](../ja/08_dynamics.md) · [Русский](../ru/08_dynamics.md)

# 8. Dynamik

Ein abgeschlossenes Quantensystem entwickelt sich gemäß einer einparametrigen Gruppe von Unitären,
ein offenes – ein System, das an eine nicht beobachtete Umgebung gekoppelt ist – gemäß einer
Halbgruppe von Kanälen. In einer Algebra mit Sektoren haben beide eine Eigenschaft, die auf einer
einzelnen Matrixalgebra unsichtbar bleibt: Ein Hamiltonoperator in $\mathcal M$ kann den Sektor
niemals ändern, die Ladungen sind also erhalten, während ein offenes System zwischen Sektoren
springen kann, etwa indem es ein Teilchen verliert. Dieses Kapitel beweist die grundlegenden
Tatsachen über beide Arten der Zeitentwicklung, stellt die Lindblad-Gleichung mit Sprüngen zwischen
Sektoren auf und zeigt, dass „kein Teilchen ging verloren“ eine Messung der Ladung ist, deren
postselektierter Zustand sich gemäß einem nicht-hermiteschen Hamiltonoperator entwickelt.

**Sie lernen**, warum die Sektorgewichte Erhaltungsgrößen sind, warum Generatoren in
Gorini-Kossakowski-Sudarshan-Lindblad-Form (GKSL-Form) Kanäle erzeugen, wie Sprünge zwischen
Sektoren in die Lindblad-Gleichung eingehen und wann die Sektorgewichte einer klassischen
Ratengleichung genügen, was die sprungfreie Entwicklung und die Postselektion bedeuten und wie man
stationäre Zustände und Relaxationsraten bestimmt.

Im Folgenden wirkt $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ auf
$\mathcal H=\bigoplus_c\mathbb C^{k_c}$, es ist $D=\sum_ck_c$, die $P_c$ sind die minimalen
zentralen Projektionen, Zustände werden durch ihre ungewichteten Dichten $\rho$ mit den
Sektorgewichten $p_c=\operatorname{Tr}\rho_c$ angegeben (Kapitel 3), $[A,B]=AB-BA$, $\{A,B\}=AB+BA$
und $\hbar=1$. Ein **Superoperator** ist eine lineare Abbildung auf $\mathcal M$ oder auf
$M_D(\mathbb C)$; alle Normen von Superoperatoren sind Operatornormen und daher submultiplikativ.
Für eine lineare Abbildung $B$ auf einem endlichdimensionalen Raum ist
$e^{tB}=\sum_{n\ge0}t^nB^n/n!$.

## 8.1 Abgeschlossene Systeme

Ein **Hamiltonoperator** ist ein selbstadjungiertes $H\in\mathcal M$. Ein Vektor entwickelt sich
gemäß der **Schrödinger-Gleichung** $i\psi'(t)=H\psi(t)$, eine Dichte gemäß der
**von-Neumann-Gleichung** $\rho'(t)=-i[H,\rho(t)]$. Beide sind lineare Gleichungen mit konstanten
Koeffizienten auf einem endlichdimensionalen Raum, und das folgende Lemma ist alles, was wir aus der
Theorie der Differentialgleichungen brauchen.

**Lemma 8.1.** Sei $B$ eine lineare Abbildung auf einem endlichdimensionalen normierten Raum $X$.
Dann gilt $\frac d{dt}e^{tB}=Be^{tB}=e^{tB}B$, $e^{(s+t)B}=e^{sB}e^{tB}$,
$\|e^{tB}\|\le e^{|t|\,\|B\|}$, und für jedes $x_0\in X$ ist $x(t)=e^{tB}x_0$ die einzige
differenzierbare Funktion $x$ auf einem Intervall $I\ni0$ mit $x'=Bx$ und $x(0)=x_0$.

*Beweis.* Die Reihe konvergiert absolut, mit der angegebenen Normschranke, und als Potenzreihe in
$t$ darf sie gliedweise differenziert werden; das liefert die erste Formel. Dann ist
$\frac d{dt}(e^{tB}e^{-tB})=Be^{tB}e^{-tB}-e^{tB}Be^{-tB}=0$, also $e^{tB}e^{-tB}=1$. Gilt $x'=Bx$,
so ist $\frac d{dt}(e^{-tB}x(t))=e^{-tB}(x'-Bx)=0$, also $e^{-tB}x(t)=x_0$ und $x(t)=e^{tB}x_0$.
Schließlich lösen $t\mapsto e^{(s+t)B}x_0$ und $t\mapsto e^{tB}e^{sB}x_0$ die Gleichung $y'=By$ mit
demselben Anfangswert. $\square$

**Proposition 8.2.** Sei $H=H^*\in\mathcal M$ und $U(t)=e^{-itH}$, definiert durch den
Funktionalkalkül aus Abschnitt 1.3 mit $f_t(\lambda)=e^{-it\lambda}$.

1. $U(t)=\bigoplus_ce^{-itH_c}$ ist ein unitäres Element von $\mathcal M$, $U(s+t)=U(s)U(t)$ und
   $U'(t)=-iHU(t)$; $U(t)$ stimmt mit der Exponentialreihe von $-itH$ überein.
2. $\psi(t)=U(t)\psi_0$ ist die einzige Lösung von $i\psi'=H\psi$ mit $\psi(0)=\psi_0$, und
   $\rho(t)=U(t)\rho_0U(t)^*$ ist die einzige Lösung in $\mathcal M$ von $\rho'=-i[H,\rho]$ mit
   $\rho(0)=\rho_0$. Ist $\rho_0$ ein Zustand, so auch jedes $\rho(t)$.
3. (*Heisenberg-Bild*) $\operatorname{Tr}(\rho(t)A)=\operatorname{Tr}(\rho_0A(t))$ mit
   $A(t)=U(t)^*AU(t)$, und $A'=i[H,A]$.

*Beweis.* (1) Man diagonalisiere jeden Block, $H_c=V_c\Lambda_cV_c^*$; dann diagonalisiert
$V=\bigoplus_cV_c$ die Matrix $H$, und $U(t)=Ve^{-it\Lambda}V^*$ ist blockdiagonal (nach Satz 1.5(a)
hängt das Ergebnis nicht von der diagonalisierenden unitären Matrix ab). Wegen $|f_t|^2=1$ und
$f_sf_t=f_{s+t}$ liefert Satz 1.5(c) die Unitarität und das Gruppengesetz. Eintragsweises
Differenzieren von $e^{-it\Lambda}$ ergibt $U'=V(-i\Lambda)e^{-it\Lambda}V^*=-iHU$, also ist $U(t)$
nach Lemma 8.1 die Reihe $e^{-itH}$ ($X=M_D$, $B$ die Linksmultiplikation mit $-iH$). (2)
$(U\psi_0)'=-iHU\psi_0$, und $(U\rho_0U^*)'=-iHU\rho_0U^*+U\rho_0U^*(iH)$, da $(U')^*=iU^*H$. Die
Eindeutigkeit ist Lemma 8.1 mit $B=-iH$ auf $\mathcal H$ bzw. $B=-i[H,\cdot\,]$ auf $\mathcal M$.
Positivität: $\langle x,U\rho_0U^*x\rangle=\langle U^*x,\rho_0U^*x\rangle\ge0$; die Spur wegen der
Zyklizität. (3) Zyklizität der Spur, und $(U^*AU)'=iU^*HAU-iU^*AHU=i[H,A(t)]$. $\square$

**Satz 8.3 (Ladungserhaltung).** Sei $t\mapsto H(t)$ eine Abbildung auf einem Intervall $I$ mit
Werten in den selbstadjungierten Elementen von $\mathcal M$.

1. Ist $\psi:I\to\mathcal H$ differenzierbar mit $i\psi'=H(t)\psi$, so ist jedes $\|P_c\psi(t)\|$
   konstant.
2. Ist $\rho:I\to\mathcal M$ differenzierbar mit $\rho'=-i[H(t),\rho]$, so ist
   $\operatorname{Tr}(Z\rho(t))$ für jedes zentrale $Z$ konstant; insbesondere ist jedes
   Sektorgewicht $p_c(t)=\operatorname{Tr}(P_c\rho(t))$ konstant.
3. Ist $H$ konstant und kommutiert $A\in\mathcal M$ mit $H$, so ist $\operatorname{Tr}(A\rho(t))$
   konstant.

*Beweis.* (1) Das Skalarprodukt ist im ersten Argument antilinear, also gilt
$\frac d{dt}\langle\psi,P_c\psi\rangle=\langle-iH\psi,P_c\psi\rangle+\langle\psi,-iP_cH\psi\rangle =i\langle\psi,HP_c\psi\rangle-i\langle\psi,P_cH\psi\rangle=0$,
da $P_c$ zentral ist. (2)
$\frac d{dt}\operatorname{Tr}(Z\rho)=-i\operatorname{Tr}(Z[H,\rho])=-i\operatorname{Tr}([Z,H]\rho)=0$
wegen der Zyklizität. (3) Dieselbe Rechnung mit $A$. $\square$

> **Physikalische Bedeutung.** Ein Hamiltonoperator in $\mathcal M$ kommutiert mit den Ladungen
> (Teilchenzahl, Magnetisierung, eine Symmetriequantenzahl), daher ist die Wahrscheinlichkeit jeder
> Ladung eine Erhaltungsgröße, selbst für getriebene, zeitabhängige Hamiltonoperatoren. Rechnerisch
> zerfällt die Zeitentwicklung in $C$ unabhängige Probleme: Eine Diagonalisierung kostet
> größenordnungsmäßig $\sum_ck_c^3$ Operationen statt $D^3$.

Die Bibliothek berechnet $U(t)$ mit einer Eigenzerlegung pro Blockgröße (`dynamics.propagator`), und
`dynamics.von_neumann` entwickelt eine Dichte in der Eigenbasis zu vielen Zeitpunkten auf einmal,
wobei $\rho(t)_{ij}=e^{-it(\lambda_i-\lambda_j)}\rho_{ij}$. Wir prüfen gegen `scipy.linalg.expm`:

```python
import math
import numpy as np
import scipy.integrate
import scipy.linalg as sla
import torch
from torch_vn_algebra import TypeIAlgebra, DensityMatrix, InterSectorChannel, channels, dynamics, states

torch.manual_seed(0)
f64, c128 = torch.float64, torch.complex128
alg = TypeIAlgebra([2, 3, 4], [2, 3, 4], precision='double')     # M_2 (+) M_3 (+) M_4, D = 9
H = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
t = 0.7
U = dynamics.propagator(H, t)                                      # exp(-itH), block by block
for c, k in enumerate(alg.k_factors):                              # independent reference
    ref = sla.expm(-1j * t * H.matrix[0, c, :k, :k].numpy())
    assert np.abs(U.matrix[0, c, :k, :k].numpy() - ref).max() < 1e-12
assert (U.adjoint() @ U - alg.identity()).matrix.abs().max() < 1e-12

rho0 = states.random_density_matrix(alg)
times = [0.0, 0.5, 1.0, 2.0]
traj = dynamics.von_neumann(rho0, H, times)                        # exact, all times at once
for s, rho_s in zip(times, traj):
    U_s = dynamics.propagator(H, s)
    assert (rho_s.matrix - (U_s @ rho0 @ U_s.adjoint()).matrix).abs().max() < 1e-12
    assert (rho_s.sector_probabilities() - rho0.sector_probabilities()).abs().max() < 1e-12   # Theorem 8.3
A = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
A_t = U.adjoint() @ A @ U                                          # Heisenberg picture
assert abs(rho0.expectation(A_t) - dynamics.von_neumann(rho0, H, [t])[0].expectation(A)) < 1e-12
```

Auf $\mathcal M$ definiert ein Vektor $\psi=\bigoplus_c\psi_c$ den Zustand der Mischung
$\bigoplus_c|\psi_c\rangle\langle\psi_c|$ (Kapitel 3, Aufgabe 3), daher müssen
`dynamics.schrodinger` (Vektoren im Blockformat `(batch, C, k_max)`) und `von_neumann` darauf
übereinstimmen. Für ein zeitabhängiges $H(t)$ verwendet die Bibliothek das klassische
Runge-Kutta-Verfahren mit fester Schrittweite (`dynamics.rk4`, `schrodinger_rk4`): Ein Schritt
multipliziert $\psi$ mit einem Polynom in Werten von $H(t)\in\mathcal M$, verschiebt also niemals
Amplitude zwischen Sektoren, ist aber nicht exakt unitär, und die Sektorgewichte driften in der
Größenordnung des Diskretisierungsfehlers. (Für einen einzelnen großen Sektor benötigt
`krylov.evolve` nur Matrix-Vektor-Produkte; siehe Kapitel 9.)

```python
rows = torch.tensor([[1, 1, 0, 0], [1, 1, 1, 0], [1, 1, 1, 1]], dtype=c128)   # entries of C^2 (+) C^3 (+) C^4
psi0 = torch.randn(1, alg.C, alg.k_max, dtype=c128) * rows
psi0 = psi0 / psi0.norm()
psi = dynamics.schrodinger(psi0, H, times)                         # shape (T, batch, C, k_max)
p_vec = (psi.abs() ** 2).sum(-1)                                   # ||P_c psi(t)||^2
assert (p_vec - p_vec[0]).abs().max() < 1e-12
mix = DensityMatrix(alg, matrix=psi0[..., :, None] * psi0[..., None, :].conj())   # the state of psi0 on M
outer = psi[..., :, None] * psi[..., None, :].conj()
assert all((o - r.matrix).abs().max() < 1e-12 for o, r in zip(outer, dynamics.von_neumann(mix, H, times)))

V = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
H_of_t = lambda s: H.matrix + math.cos(3 * s) * V.matrix          # driven, in M at every time
psi_rk = dynamics.schrodinger_rk4(psi0, H_of_t, times, substeps=50, algebra=alg)
drift = ((psi_rk.abs() ** 2).sum(-1) - p_vec[0]).abs().max().item()
print(f"RK4, time-dependent H: largest change of a sector weight {drift:.1e}")
assert drift < 1e-8 and (dynamics.schrodinger_rk4(psi0, H, times, substeps=50) - psi).abs().max() < 1e-7

```

Die Sektorgewichte ändern sich unter RK4 um $1.1\cdot10^{-9}$ und unter den exakten Lösern um
nichts, was über Rundungsfehler hinausgeht.

## 8.2 Quantendynamische Halbgruppen

CP- und TP-Abbildungen sowie Kanäle sind in Abschnitt 6.1 definiert ([Kapitel
6](06_quantum_channels.md)). Neben Satz 6.2 und Proposition 6.4 (eine Abbildung
$A\mapsto\sum_lK_lAK_l^*$, auch zwischen Sektoren, ist CP) benötigen wir drei
Abgeschlossenheitseigenschaften.

**Lemma 8.4.** (a) Summen, nichtnegative Vielfache und Kompositionen von CP-Abbildungen sind CP. Ist
$\Phi$ CP auf $M_D$ und bildet eine $*$-Unteralgebra $\mathcal N\subseteq M_D$ (etwa $\mathcal M$
oder $P_c\mathcal MP_c\cong M_{k_c}$ mit Einselement $P_c$) in sich ab, so ist $\Phi|_{\mathcal N}$
CP. (b) Konvergieren CP-Abbildungen $\Phi_m$ gegen $\Phi$ (d. h. $\Phi_m(A)\to\Phi(A)$ für alle
$A$), so ist $\Phi$ CP. (c) Ist $s\mapsto\Psi(s)$ eine stetige Abbildung von $[a,b]$ in die
CP-Abbildungen, so ist $\int_a^b\Psi(s)\,ds$ CP.

*Beweis.* (a) $\Phi\mapsto\Phi_n$ ist linear und $(\Phi\Psi)_n=\Phi_n\Psi_n$; Summen positiver
Matrizen sind positiv. Zur Einschränkung: $M_n(\mathcal N)\subseteq M_n(M_D)$, und die Positivität
ist in beiden dieselbe Bedingung (nach einer Permutation der Basis ist ein Element von
$M_n(\mathcal N)$ eine direkte Summe der Blöcke $[(A_{ij})_c]$ und von Nullen). (b) Die positiven
Matrizen bilden eine abgeschlossene Menge ($\langle x,Xx\rangle\ge0$ für jedes $x$ ist eine
abgeschlossene Bedingung), und $(\Phi_m)_n(X)\to\Phi_n(X)$. (c) Die Riemann-Summen sind nach (a) CP
und konvergieren gegen das Integral, da $\Psi$ gleichmäßig stetig ist; man wende (b) an. $\square$

**Definition.** Eine **quantendynamische Halbgruppe** auf $\mathcal M$ ist eine Familie
$(T_t)_{t\ge0}$ von Kanälen $\mathcal M\to\mathcal M$ mit $T_0=\mathrm{id}$, $T_{s+t}=T_sT_t$ und
stetigem $t\mapsto T_t$.

**Lemma 8.5.** Jede stetige Halbgruppe von Superoperatoren ($T_0=\mathrm{id}$, $T_{s+t}=T_sT_t$) ist
von der Form $T_t=e^{t\mathcal L}$ für genau einen Superoperator $\mathcal L$, ihren **Generator**,
und es gilt $\mathcal L=\lim_{t\to0}(T_t-\mathrm{id})/t$.

*Beweis.* Es gilt $V_h=\frac1h\int_0^hT_s\,ds\to\mathrm{id}$ für $h\to0$, also ist $V_h$ für ein
festes kleines $h$ invertierbar. Nun ist $T_tV_h=\frac1h\int_t^{t+h}T_s\,ds$ in $t$ differenzierbar
mit der Ableitung $\frac1h(T_{t+h}-T_t)=T_t\frac{T_h-\mathrm{id}}h$. Folglich ist $T_t$
differenzierbar und $T_t'=T_t\mathcal L$ mit $\mathcal L=\frac{T_h-\mathrm{id}}hV_h^{-1}$, und nach
Lemma 8.1 (mit $B$ der Komposition von rechts mit $\mathcal L$) ist $T_t=e^{t\mathcal L}$.
Eindeutigkeit: $\mathcal L=T'(0)$. $\square$

**Satz 8.6 (GKSL-Generatoren).** Seien $H=H^*$ und $L_1,\dots,L_J$ in $M_D$, und sei
$$
\mathcal L(\rho)=-i[H,\rho]+\sum_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
$$
Dann ist $(e^{t\mathcal L})_{t\ge0}$ eine quantendynamische Halbgruppe auf $M_D$. Genauer: Man setze
$H_{\mathrm{eff}}=H-\frac i2\sum_jL_j^*L_j$,
$\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)$ und
$\Phi(\rho)=\sum_jL_j\rho L_j^*$, sodass $\mathcal L=\mathcal L_0+\Phi$, und definiere
$$
S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},\qquad
S_n(t)=\int_0^tS_0(t-s)\,\Phi\,S_{n-1}(s)\,ds\quad(n\ge1).
$$
Dann ist jedes $S_n(t)$ CP, und $e^{t\mathcal L}=\sum_{n\ge0}S_n(t)$ (die **Dyson-Entwicklung**),
mit $\|S_n(t)\|\le e^{t\|\mathcal L_0\|}(t\|\Phi\|)^n/n!$.

*Beweis.* *Schritt 1.* $X(t)=e^{-itH_{\mathrm{eff}}}\rho e^{itH_{\mathrm{eff}}^*}$ erfüllt
$X'=\mathcal L_0(X)$ und $X(0)=\rho$, also ist $S_0(t)=e^{t\mathcal L_0}$ nach Lemma 8.1; es ist CP
mit dem einzigen Kraus-Operator $e^{-itH_{\mathrm{eff}}}$. *Schritt 2.* Man schreibe
$a=\|\mathcal L_0\|$, $b=\|\Phi\|$. Durch Induktion sind die $S_n$ stetig, und es gilt
$\|S_n(t)\|\le\int_0^te^{(t-s)a}\,b\,e^{sa}\frac{(bs)^{n-1}}{(n-1)!}ds=e^{ta}\frac{(bt)^n}{n!}$; die
Reihe konvergiert also gleichmäßig auf jedem $[0,T]$. *Schritt 3.* Da
$S_n(t)=S_0(t)\int_0^te^{-s\mathcal L_0}\Phi S_{n-1}(s)\,ds$, liefern die Produktregel und der
Hauptsatz der Differential- und Integralrechnung $S_n'=\mathcal L_0S_n+\Phi S_{n-1}$ ($n\ge1$) und
$S_0'=\mathcal L_0S_0$. Wegen der Schranken konvergiert die Reihe der Ableitungen gleichmäßig auf
$[0,T]$, also ist $Y(t)=\sum_nS_n(t)$ differenzierbar mit $Y'=\mathcal L_0Y+\Phi Y=\mathcal LY$,
$Y(0)=\mathrm{id}$, und nach Lemma 8.1 ist $Y(t)=e^{t\mathcal L}$. *Schritt 4.* $\Phi$ ist CP (Satz
6.2), also ist durch Induktion der Integrand $s\mapsto S_0(t-s)\Phi S_{n-1}(s)$ eine stetige Familie
von CP-Abbildungen und $S_n(t)$ CP (Lemma 8.4); ebenso die Partialsummen und ihr Grenzwert
$e^{t\mathcal L}$. *Schritt 5.* Wegen der Zyklizität ist
$\operatorname{Tr}\mathcal L(\rho)=\sum_j\big(\operatorname{Tr}(L_j\rho L_j^*)-\operatorname{Tr}(L_j^*L_j\rho)\big)=0$,
folglich
$\frac d{dt}\operatorname{Tr}e^{t\mathcal L}\rho=\operatorname{Tr}\mathcal L(e^{t\mathcal L}\rho)=0$:
Jedes $e^{t\mathcal L}$ ist TP. Halbgruppengesetz und Stetigkeit folgen aus Lemma 8.1. $\square$

**Satz 8.7 (Gorini-Kossakowski-Sudarshan; Lindblad; ohne Beweis).** Umgekehrt hat der Generator
jeder quantendynamischen Halbgruppe auf $M_D$ die Form aus Satz 8.6. [V. Gorini, A. Kossakowski, E.
C. G. Sudarshan, *J. Math. Phys.* 17, 821 (1976); G. Lindblad, *Commun. Math. Phys.* 48, 119
(1976).] Für eine normstetige Halbgruppe unitaler CP-Abbildungen auf einer allgemeinen
von-Neumann-Algebra $\mathcal M$ zeigten Christensen und Evans [*J. London Math. Soc.* 20, 358
(1979)], dass der Generator die Form $A\mapsto\Psi(A)+K^*A+AK$ mit einer CP-Abbildung $\Psi$ und
$K\in\mathcal M$ hat.

Auf $\mathcal M=\bigoplus_cM_{k_c}$ wende man dies auf die Halbgruppe $(T_t^*)$ im Heisenberg-Bild
an, deren Abbildungen CP und unital sind (Proposition 6.8, Korollar 6.9(b)) und deren Generator
$\mathcal L^*$ ist. Geht man wieder zu den ungewichteten dualen Abbildungen über, so erhält man
$\mathcal L(\rho)=\Phi(\rho)+G\rho+\rho G^*$ mit $\Phi=\Psi^*$ CP auf $\mathcal M$ (Kraus-Operatoren
zwischen Sektoren, Proposition 6.4) und $G=K$. Die Spurerhaltung,
$\operatorname{Tr}\big((\Phi^*(1)+G+G^*)\rho\big)=0$ für alle $\rho\in\mathcal M$, erzwingt
$G+G^*=-\Phi^*(1)$, d. h. $G=-iH-\frac12\Phi^*(1)$ mit $H=H^*\in\mathcal M$: genau die Generatoren
des nächsten Abschnitts.

> **Physikalische Bedeutung.** Ein System, das schwach mit einer großen Umgebung wechselwirkt, deren
> Korrelationen schnell abklingen (die Born-Markov-Näherung), verliert das Gedächtnis an seine
> Vergangenheit: Sein Zustand zur Zeit $s+t$ geht aus dem Zustand zur Zeit $s$ durch einen Kanal
> hervor, der nur von $t$ abhängt. Satz 8.7 besagt, dass eine solche Dynamik stets aus einem
> Hamiltonschen Anteil und „Sprüngen“ $L_j$ (Photonenemission, Teilchenverlust, Dephasierung)
> besteht, mit Raten, die durch die Normierung der $L_j$ festgelegt sind.

Die Bibliothek hat vier Lindblad-Routinen: `lindblad_rhs` gibt die rechte Seite als Funktion für
`rk4` zurück; `lindblad_evolve` integriert mit RK4; `lindblad_superoperator` gibt für jeden Sektor
die $k_c^2\times k_c^2$-Matrix von $\mathcal L$ auf dem zeilenweise gebildeten
$\operatorname{vec}\rho_c$ zurück; und `lindblad_channel` berechnet $e^{t\mathcal L}$ exakt, Sektor
für Sektor, als `Channel`. Wir prüfen Satz 8.6 für beliebige Sprünge auf $M_3$: Die Choi-Matrix von
$e^{t\mathcal L}$ ist positiv (Satz 6.2), die Abbildung ist TP, das Halbgruppengesetz gilt, RK4
stimmt mit der Exponentialfunktion überein, und der Sprunganteil
$e^{t\mathcal L}-S_0(t)=\sum_{n\ge1}S_n(t)$ ist selbst CP:

```python
M3 = TypeIAlgebra([3], [3], precision='double')
H3 = M3.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
jumps = [M3.operator(torch.randn(1, 1, 3, 3, dtype=c128) / 2) for _ in range(2)]     # arbitrary L_j
T = lambda s: dynamics.lindblad_channel(M3, H3, jumps, t=s)        # exp(sL), exact
S = lambda s: T(s).superoperator()[0][0]                           # its 9 x 9 matrix on row-major vec
choi = lambda X: X.reshape(3, 3, 3, 3).permute(2, 0, 3, 1).reshape(9, 9)   # J_{(i,a),(j,b)} = S_{(a,b),(i,j)}
assert (choi(S(0.8)) - T(0.8).choi()[0][0]).abs().max() < 1e-12
assert torch.linalg.eigvalsh(choi(S(0.8))).min() > -1e-12         # completely positive
assert T(0.8).trace_preservation_error().max() < 1e-12             # trace preserving
assert (S(0.3) @ S(0.5) - S(0.8)).abs().max() < 1e-12              # semigroup

rho3 = states.random_density_matrix(M3)
rk = dynamics.lindblad_evolve(rho3, H3, jumps, [0.0, 0.8], substeps=40)[-1]
assert (rk.matrix - T(0.8)(rho3).matrix).abs().max() < 1e-7        # RK4 against exp(tL)
assert abs(dynamics.lindblad_rhs(H3, jumps)(0.0, rho3.matrix)[0, 0].trace()) < 1e-12   # Tr L(rho) = 0

H_eff = H3.matrix - 0.5j * sum(L.adjoint().matrix @ L.matrix for L in jumps)
S0 = channels.Channel(M3, torch.linalg.matrix_exp(-0.8j * H_eff)[:, None]).superoperator()[0][0]
jump_part = torch.linalg.eigvalsh(choi(S(0.8) - S0))
print(f"smallest eigenvalue of the Choi matrix of the jump part: {jump_part.min().item():.1e}")
assert jump_part.min() > -1e-12                                    # sum_{n >= 1} S_n(t) is CP
```

## 8.3 Lindblad-Gleichungen in einer Algebra mit Sektoren

Auf $\mathcal M$ gibt es zwei Arten von Sprüngen. Ein Sprung **innerhalb der Sektoren** ist ein
Operator $L_j\in\mathcal M$. Ein Sprung **zwischen Sektoren** ist eine CP-Abbildung
$\Phi:\mathcal M\to\mathcal M$ mit Kraus-Operatoren $K^{dc}_l:\mathbb C^{k_c}\to\mathbb C^{k_d}$
(Proposition 6.4, ein `InterSectorChannel`): $\Phi(\rho)_d=\sum_{c,l}K^{dc}_l\rho_cK^{dc*}_l$. Nach
Proposition 6.8 erfüllt ihre ungewichtete duale Abbildung $\Phi^*(1)_c=\sum_{d,l}K^{dc*}_lK^{dc}_l$,
ein blockdiagonales Element von $\mathcal M$. Mit Raten $g_j,g_m>0$ lautet der Generator von
`lindblad_rhs`
$$
\mathcal L(\rho)=-i[H,\rho]+\sum_jg_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big)
+\sum_mg_m\Big(\Phi_m(\rho)-\tfrac12\{\Phi_m^*(1),\rho\}\Big).
$$

**Proposition 8.8.** Sei $J^{dc}_l\in M_D$ die Matrix mit dem Block $K^{dc}_l$ an der Position
$(d,c)$ und Nullen sonst.

1. $\mathcal L$ bildet $\mathcal M$ in $\mathcal M$ ab und ist die Einschränkung auf $\mathcal M$
   des GKSL-Generators $\hat{\mathcal L}$ auf $M_D$ mit den Sprüngen $\sqrt{g_j}L_j$ und
   $\sqrt{g_m}J^{dc}_{m,l}$ (alle $m,l,d,c$). Folglich ist $(e^{t\mathcal L})_{t\ge0}$ eine
   quantendynamische Halbgruppe auf $\mathcal M$.
2. Sei stattdessen $\tilde{\mathcal L}$ der GKSL-Generator auf $M_D$, in dem die Blöcke eines
   Kraus-Index zu einem einzigen Sprung summiert werden, $L_{m,l}=\sum_{d,c}J^{dc}_{m,l}$. Dann gilt
   $\mathcal L(\rho)=E\big(\tilde{\mathcal L}(\rho)\big)$ für $\rho\in\mathcal M$, wobei
   $E(X)=\sum_dP_dXP_d$ die bedingte Erwartung von $M_D$ auf $\mathcal M$ ist; im Allgemeinen ist
   $\tilde{\mathcal L}(\rho)\notin\mathcal M$.

*Beweis.* Da $\rho=\sum_cP_c\rho P_c$ und $J^{dc}P_{c'}=\delta_{cc'}J^{dc}$, ist
$J^{dc}\rho J^{dc*}$ der Block $K^{dc}\rho_cK^{dc*}$ an der Position $(d,d)$, und $J^{dc*}J^{dc}$
ist $K^{dc*}K^{dc}$ an der Position $(c,c)$. Summation über $d,l$ liefert $\Phi(\rho)$ und
$\Phi^*(1)$; das beweist die Identität in (1), und alle Terme liegen in $\mathcal M$. Also bilden
$\hat{\mathcal L}$ und alle seine Potenzen $\mathcal M$ in $\mathcal M$ ab, es gilt
$e^{t\hat{\mathcal L}}|_{\mathcal M}=e^{t\mathcal L}$, und dies ist nach Satz 8.6 und Lemma 8.4(a)
ein Kanal von $\mathcal M$. (2) Es ist $L\rho L^*=\sum_{d,d',c}J^{dc}\rho J^{d'c*}$ für
$L=\sum_{d,c}J^{dc}$, und $E$ behält die Terme mit $d=d'$. Für $B\in\mathcal M$ gilt
$P_dXBP_d=P_dXP_dB$, also $E(XB)=E(X)B$ und $E(BX)=BE(X)$; somit ist
$E\{L^*L,\rho\}=\{E(L^*L),\rho\}$, und $E(L^*L)=\sum_cP_cL^*LP_c=\sum_{c,d}J^{dc*}J^{dc}$, weil
$J^{dc*}J^{d'c}=0$ für $d\ne d'$. Schließlich ist $E(-i[H,\rho])=-i[H,\rho]$. $\square$

> **Physikalische Bedeutung.** Teil (2) besagt, dass die Lindblad-Gleichung auf $\mathcal M$ das
> ist, was von einer Gleichung auf $M_D$ übrig bleibt, wenn die von den Sprüngen erzeugten
> Kohärenzen zwischen Sektoren sofort gelöscht werden – der Gehalt einer Superauswahlregel. Beim
> Teilchenverlust wird nichts gelöscht: Der Vernichtungsoperator $a_j$ senkt die Teilchenzahl um
> genau eins, also sind $a_j\rho a_j^*$ und $a_j^*a_j=n_j$ blockdiagonal, und die Gleichung auf
> $\mathcal M$ ist die übliche.

**Proposition 8.9 (Sektorgewichte).** Sei $\rho(t)=e^{t\mathcal L}\rho_0$, und seien $p(t)$ seine
Sektorgewichte.

1. Gibt es keine Sprünge zwischen Sektoren, so ist jedes $p_c(t)$ konstant.
2. Gilt für jedes $m$ und jedes Paar $d\ne c$ die Beziehung
   $\sum_lK^{dc*}_{m,l}K^{dc}_{m,l}=r_{m,dc}1_{k_c}$ mit Zahlen $r_{m,dc}$, so ist $p'=Qp$ mit
   $Q_{dc}=\sum_mg_mr_{m,dc}$ für $d\ne c$ und $Q_{cc}=-\sum_{d\ne c}Q_{dc}$. Folglich hängt
   $p(t)=e^{tQ}p(0)$ weder von $H$ noch von den $L_j$ noch von den Zuständen innerhalb der Sektoren
   ab.

*Beweis.* $p_d'=\operatorname{Tr}(P_d\mathcal L(\rho))$. Der Hamiltonsche Term verschwindet wie in
Satz 8.3. Für $L_j\in\mathcal M$ ist
$\operatorname{Tr}(P_dL_j\rho L_j^*)-\operatorname{Tr}(P_dL_j^*L_j\rho) =\operatorname{Tr}(L_j^*P_dL_j\rho)-\operatorname{Tr}(L_j^*L_jP_d\rho)=0$,
da $P_d$ zentral ist. Für einen Sprung $\Phi$ zwischen Sektoren ist
$\operatorname{Tr}(P_d\Phi(\rho))=\sum_c\operatorname{Tr}\big(\sum_lK^{dc*}_lK^{dc}_l\rho_c\big)$
und
$\operatorname{Tr}(P_d\tfrac12\{\Phi^*(1),\rho\})=\sum_{d'}\operatorname{Tr}\big(\sum_lK^{d'd*}_lK^{d'd}_l\rho_d\big)$;
die Terme mit $c=d$ und $d'=d$ heben sich weg. Das beweist (1), und unter der Voraussetzung von (2)
sind die verbleibenden Terme $\sum_{c\ne d}r_{dc}p_c-\sum_{d'\ne d}r_{d'd}p_d$. Lemma 8.1 liefert
$p(t)=e^{tQ}p(0)$. $\square$

Unter der Voraussetzung von (2) ist für $d\ne c$ die Zahl $r_{m,dc}$ für jeden Zustand $\rho$ der
Eintrag $T_{dc}(\rho)$ der Übergangsmatrix von $\Phi_m$ (definiert vor Proposition 6.12). Für einen
einzigen Sprung mit der Rate $1$ ist also $Q=R-\operatorname{diag}(\sum_dR_{dc})_c$ mit $R$ =
`transition_matrix()` (die Diagonale von $R$ hebt sich weg).

**Beispiel: Teilchenverlust.** Zwei bosonische Moden mit höchstens zwei Teilchen; die Sektoren sind
$N=0,1,2$ mit den Basen $|00\rangle$; $|10\rangle,|01\rangle$; $|20\rangle,|11\rangle,|02\rangle$,
also $\mathcal M=M_1\oplus M_2\oplus M_3$. Der Bose-Hubbard-Dimer
$H=J(a_1^*a_2+a_2^*a_1)+\frac U2\sum_jn_j(n_j-1)$ kommutiert mit $N$, und jede Mode verliert
Teilchen mit der Rate $\gamma$: ein Sprung zwischen Sektoren mit den Kraus-Operatoren $a_1,a_2$,
eingeschränkt auf $N\to N-1$. Da $a_1^*a_1+a_2^*a_2=N$, ist Proposition 8.9(2) mit
$Q_{N-1,N}=\gamma N$ anwendbar, und ausgehend von $N=2$ gilt
$$
p_2(t)=q^2,\qquad p_1(t)=2q(1-q),\qquad p_0(t)=(1-q)^2,\qquad q=e^{-\gamma t},
$$
eine Binomialverteilung: Jedes Teilchen überlebt unabhängig mit der Wahrscheinlichkeit
$e^{-\gamma t}$, was auch immer das Tunneln (Hopping) und die Wechselwirkung bewirken (man prüfe
$p'=Qp$ durch Differenzieren).

```python
fock = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double', charges=[0, 1, 2])   # N = 0, 1, 2
r2 = math.sqrt(2)
a1 = {(0, 1): [[1.0, 0.0]], (1, 2): [[r2, 0.0, 0.0], [0.0, 1.0, 0.0]]}           # a_1 : N -> N - 1
a2 = {(0, 1): [[0.0, 1.0]], (1, 2): [[0.0, 1.0, 0.0], [0.0, 0.0, r2]]}           # a_2 : N -> N - 1
blocks = lambda a: {dc: torch.tensor(K, dtype=f64) for dc, K in a.items()}
loss = InterSectorChannel.from_blocks(fock, fock, {dc: [blocks(a1)[dc], blocks(a2)[dc]] for dc in a1})
J, U_int = 1.0, 0.8                                                # hopping and interaction
H_bh = fock.from_blocks([torch.zeros(1, 1, dtype=f64), torch.tensor([[0.0, J], [J, 0.0]], dtype=f64),
                         torch.tensor([[U_int, r2 * J, 0.0], [r2 * J, 0.0, r2 * J], [0.0, r2 * J, U_int]], dtype=f64)])

assert (loss.adjoint()(fock.identity()).matrix - fock.central([0.0, 1.0, 2.0]).matrix).abs().max() < 1e-12
R = loss.transition_matrix()[0]                                    # r_dc (columns: from sector c)
Q = R - torch.diag(R.sum(0))                                       # Proposition 8.9(2), for rate 1
print(Q)

gamma = 0.5
psi2 = torch.randn(3, dtype=c128)
psi2 = psi2 / psi2.norm()
m0 = torch.zeros(1, 3, 3, 3, dtype=c128)
m0[0, 2] = torch.outer(psi2, psi2.conj())                          # a pure state with N = 2
ts = [0.0, 0.5, 1.0, 2.0, 4.0]
traj = dynamics.lindblad_evolve(DensityMatrix(fock, matrix=m0), H_bh, [loss], ts, rates=[gamma], substeps=100)
for s, rho_s in zip(ts, traj):
    q = math.exp(-gamma * s)
    p = rho_s.sector_probabilities()[0]
    assert (p - torch.tensor([(1 - q) ** 2, 2 * q * (1 - q), q ** 2], dtype=f64)).abs().max() < 1e-8
    assert (p - torch.from_numpy(sla.expm(gamma * s * Q.numpy())[:, 2])).abs().max() < 1e-8
neg = min(w.min().item() for rho_s in traj for w in rho_s.eigenvalues())
print(f"smallest eigenvalue of rho(t) under RK4: {neg:.1e}")
assert neg > -1e-7                                                 # positive up to the RK4 error
```

RK4 ist nicht exakt vollständig positiv: Die Eigenwerte von $\rho(t)$ können in der Größenordnung
des Diskretisierungsfehlers negativ werden (hier $-2.5\cdot10^{-8}$, mit Schrittweiten bis $0.02$);
`lindblad_channel` ist exakt, behandelt aber nur Sprünge innerhalb der Sektoren. Schließlich prüfen
wir Proposition 8.8(2) für eine zufällige Abbildung zwischen den Sektoren, im Vergleich mit dem
GKSL-Generator auf $M_6$, der mit numpy aufgebaut wird:

```python
Phi = channels.random_inter_sector_channel(fock, fock, kraus_rank=2)
rho = states.random_density_matrix(fock)
on_M = dynamics.lindblad_rhs(None, [Phi], algebra=fock)(0.0, rho.matrix)[0].numpy()
ks, K = fock.k_factors, Phi.kraus[0].numpy()                       # K[l, d, c] = K_l^{dc}, zero padded
embed = lambda Y: sla.block_diag(*[Y[c, :k, :k] for c, k in enumerate(ks)])    # (+)_c Y_c in M_6
Ls = [np.block([[K[l, d, c, :kd, :kc] for c, kc in enumerate(ks)] for d, kd in enumerate(ks)])
      for l in range(2)]                                           # L_l = sum_{d,c} J_l^{dc}
X = embed(rho.matrix[0].numpy())
full = sum(L @ X @ L.conj().T - 0.5 * (L.conj().T @ L @ X + X @ L.conj().T @ L) for L in Ls)   # GKSL on M_6
E_full = full * embed(np.ones((3, 3, 3)))                          # E(X) = sum_d P_d X P_d
assert np.abs(E_full - embed(on_M)).max() < 1e-12                  # Proposition 8.8(2)
print(f"coherences between sectors created on M_6: {np.abs(full - E_full).max():.2f}")
```

## 8.4 Postselektion und die sprungfreie Entwicklung

Man fasse die Sprünge des Generators aus Abschnitt 8.3 zu
$\Phi_{\mathrm{tot}}(\rho)=\sum_jg_jL_j\rho L_j^*+\sum_mg_m\Phi_m(\rho)$ zusammen, mit
$A=\Phi_{\mathrm{tot}}^*(1)=\sum_jg_jL_j^*L_j+\sum_mg_m\Phi_m^*(1)\ge0$ in $\mathcal M$ und
$H_{\mathrm{eff}}=H-\frac i2A$. Nach Proposition 8.8(1) und Satz 8.6 (dessen $S_n$ den Raum
$\mathcal M$ in $\mathcal M$ abbilden, da $H_{\mathrm{eff}}\in\mathcal M$) gilt
$$
\rho(t)=\sum_{n\ge0}S_n(t)\rho,\qquad S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},
\qquad S_n(t)=\int_0^tS_0(t-s)\,\Phi_{\mathrm{tot}}\,S_{n-1}(s)\,ds .
$$

**Proposition 8.10.** Für einen Zustand $\rho$ sei $w_n(t)=\operatorname{Tr}S_n(t)\rho$. Dann ist
$w_n(t)\ge0$, $\sum_nw_n(t)=1$, und $w_0$ ist monoton nicht wachsend mit
$w_0'(t)=-\operatorname{Tr}\big(A\,S_0(t)\rho\big)$. Für einen reinen Zustand
$\rho=|\psi\rangle\langle\psi|$ ist $w_0(t)=\|e^{-itH_{\mathrm{eff}}}\psi\|^2$.

*Beweis.* $S_n(t)$ ist CP, also positiv, und $\sum_nw_n=\operatorname{Tr}\rho(t)=1$. Mit
$X=S_0(t)\rho\ge0$ ist
$w_0'=\operatorname{Tr}(-iH_{\mathrm{eff}}X+iXH_{\mathrm{eff}}^*)=-\operatorname{Tr}(AX)$, und
$\operatorname{Tr}(AX)=\operatorname{Tr}(X^{1/2}AX^{1/2})\ge0$. Die letzte Behauptung ist
$\operatorname{Tr}|\phi\rangle\langle\phi|=\|\phi\|^2$. $\square$

> **Physikalische Bedeutung.** Ein Detektor registriere jeden Sprung – jedes emittierte Photon,
> jedes verlorene Atom. Dann ist $w_n(t)$ die Wahrscheinlichkeit für genau $n$ Klicks bis zur Zeit
> $t$, $S_n(t)\rho/w_n(t)$ der Zustand unter der Bedingung von $n$ Klicks und $S_0(t)\rho/w_0(t)$
> der Zustand unter der Bedingung, dass es **keinen Klick** gab: Auch ein Nullergebnis ändert den
> Zustand, und zwar über das nicht-hermitesche $H_{\mathrm{eff}}$. Diese Deutung als Zählprozess [E.
> B. Davies, *Quantum Theory of Open Systems*, 1976; M. D. Srinivas und E. B. Davies, *Optica Acta*
> 28, 981 (1981)] ist die Grundlage der Methode der Quantentrajektorien
> (Monte-Carlo-Wellenfunktionen) [J. Dalibard, Y. Castin, K. Mølmer, *Phys. Rev. Lett.* 68, 580
> (1992)]; hier verwenden wir nur ihren mathematischen Teil.

Bei Verlusten ist das Zählen der Klicks dasselbe wie das Messen der Ladung:

**Proposition 8.11.** Die Sektoren seien mit $N=0,\dots,C-1$ bezeichnet, es gebe keine Sprünge
innerhalb der Sektoren, jeder Sprung senke die Ladung um eins ($K^{dc}_{m,l}=0$ außer für $d=c-1$),
und es sei $\rho=P_{N_0}\rho P_{N_0}$. Dann gilt $S_n(t)\rho=P_{N_0-n}\,\rho(t)\,P_{N_0-n}$ für
$0\le n\le N_0$ und $S_n(t)\rho=0$ für $n>N_0$. Insbesondere ist die Lüders-Aktualisierung von
$\rho(t)$ mit $P_{N_0}$ („kein Teilchen verloren“) der sprungfreie Zustand:
$$
\frac{P_{N_0}\rho(t)P_{N_0}}{p_{N_0}(t)}=\frac{e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*}}{w_0(t)},
\qquad p_{N_0}(t)=w_0(t).
$$

*Beweis.* $S_0(s)$ bildet Operatoren mit Träger im Sektor $N$ auf Operatoren mit Träger im Sektor
$N$ ab, da $H_{\mathrm{eff}}\in\mathcal M$; $\Phi_{\mathrm{tot}}$ bildet den Sektor $N\ge1$ auf den
Sektor $N-1$ ab und annulliert den Sektor $0$. Durch Induktion hat $S_n(s)\rho$ für alle $s$ seinen
Träger im Sektor $N_0-n$ und ist $0$ für $n>N_0$. Die Kompression von $\rho(t)=\sum_nS_n(t)\rho$ mit
$P_{N_0-n}$ behält genau den $n$-ten Term. $\square$

Die Funktion `states.lueders_update(rho, P)` der Bibliothek gibt
$P\rho P/\operatorname{Tr}(P\rho P)$ und $\operatorname{Tr}(\rho P)$ zurück. Nun mögen die beiden
Moden mit verschiedenen Raten $\gamma_1>\gamma_2$ Teilchen verlieren; dann ist
$A=\gamma_1n_1+\gamma_2n_2$ nicht zentral, Proposition 8.9(2) ist nicht anwendbar, und die
Postselektion auf $N=2$ liefert eine echt nicht-hermitesche Entwicklung. Wir prüfen Proposition 8.11
für $n=0$ und – indem wir das Integral für $S_1$ mit der Simpson-Regel berechnen – für $n=1$:

```python
g1, g2 = 0.6, 0.2                                                  # mode 1 leaks faster
loss1 = InterSectorChannel.from_blocks(fock, fock, {dc: [K] for dc, K in blocks(a1).items()})
loss2 = InterSectorChannel.from_blocks(fock, fock, {dc: [K] for dc, K in blocks(a2).items()})
m1 = torch.zeros(1, 3, 3, 3, dtype=c128)
m1[0, 2, 0, 0] = 1.0                                               # |20>: both particles in mode 1
ts = [0.0, 0.5, 1.0, 2.0]
traj = dynamics.lindblad_evolve(DensityMatrix(fock, matrix=m1), H_bh, [loss1, loss2], ts, rates=[g1, g2],
                                substeps=100)

n_1 = [np.diag(x) for x in ([0.0], [1.0, 0.0], [2.0, 1.0, 0.0])]   # n_1 in each sector; n_2 = N - n_1
H_eff = [H_bh.matrix[0, N, :N + 1, :N + 1].numpy() - 0.5j * (g1 * n + g2 * (N * np.eye(N + 1) - n))
         for N, n in enumerate(n_1)]

def S0(s, N, X):                                  # e^{-isH_eff} X e^{isH_eff^*} in sector N, s a time or an array
    E = torch.linalg.matrix_exp(torch.from_numpy(-1j * np.asarray(s)[..., None, None] * H_eff[N])).numpy()
    return E @ X @ E.conj().swapaxes(-1, -2)

rho_20 = m1[0, 2].numpy()
for s, rho_s in zip(ts, traj):
    w0 = np.trace(S0(s, 2, rho_20)).real                           # probability of no loss
    post, prob = states.lueders_update(rho_s, fock.central([0.0, 0.0, 1.0]))   # post-select N = 2
    assert abs(prob.item() - w0) < 1e-8                            # Proposition 8.11, n = 0
    assert np.abs(post.matrix[0, 2].numpy() - S0(s, 2, rho_20) / w0).max() < 1e-7
closed = dynamics.von_neumann(DensityMatrix(fock, matrix=m1), H_bh, [2.0])[0]   # no losses at all
dist = 0.5 * np.abs(np.linalg.eigvalsh(S0(2.0, 2, rho_20) / w0 - closed.matrix[0, 2].numpy())).sum()
print(f"t = 2: P(no loss) = {w0:.3f} (without hopping: {math.exp(-2 * g1 * 2):.3f}); "
      f"trace distance post-selected/closed system: {dist:.3f}")

K1, K2 = blocks(a1)[(1, 2)].numpy(), blocks(a2)[(1, 2)].numpy()   # S_1(2) rho = int_0^2 S_0(2-s) Phi S_0(s) rho ds
jump = lambda X: g1 * K1 @ X @ K1.T + g2 * K2 @ X @ K2.T
s_grid = np.linspace(0.0, 2.0, 201)
S1 = scipy.integrate.simpson(S0(2.0 - s_grid, 1, jump(S0(s_grid, 2, rho_20))), x=s_grid, axis=0)
assert np.abs(S1 - traj[-1].matrix[0, 1, :2, :2].numpy()).max() < 1e-8     # Proposition 8.11, n = 1
```

Bei $t=2$ beträgt die Wahrscheinlichkeit, kein Teilchen zu verlieren, $0.230$; ohne Tunneln wäre sie
$e^{-2\gamma_1t}=0.091$, hier hängen die Sektorgewichte also sehr wohl vom Hamiltonoperator ab. Der
postselektierte Zustand hat vom Zustand des abgeschlossenen Systems den Spurabstand $0.135$: Kein
Teilchen zu verlieren ist eine Information.

## 8.5 Stationäre Zustände und Relaxation

Ein Zustand $\sigma$ heißt **stationär**, wenn $e^{t\mathcal L}\sigma=\sigma$ für alle $t\ge0$ gilt,
äquivalent $\mathcal L(\sigma)=0$ (man differenziere bei $t=0$; umgekehrt verwende man die Reihe).

**Proposition 8.12.** Jede quantendynamische Halbgruppe auf $\mathcal M$ besitzt einen stationären
Zustand.

*Beweis.* Man fixiere einen Zustand $\rho$ und setze
$\sigma_\tau=\frac1\tau\int_0^\tau T_s\rho\,ds$. Dies ist ein Zustand (ein Grenzwert von
Konvexkombinationen von Zuständen; die Zustände bilden eine abgeschlossene Menge). Wegen der
Linearität ist
$T_t\sigma_\tau-\sigma_\tau=\frac1\tau\big(\int_\tau^{\tau+t}-\int_0^t\big)T_s\rho\,ds$, und dessen
Spurnorm ist höchstens $2t/\tau$, da $\|T_s\rho\|_1=1$. Die Zustände bilden eine beschränkte
abgeschlossene Menge in einem endlichdimensionalen Raum, also konvergiert eine Folge
$\sigma_{\tau_m}$, $\tau_m\to\infty$, gegen einen Zustand $\sigma$ (Bolzano-Weierstraß), und
$T_t\sigma-\sigma=\lim_m(T_t\sigma_{\tau_m}-\sigma_{\tau_m})=0$. $\square$

**Korollar 8.13.** Gibt es keine Sprünge zwischen Sektoren, so spannen die stationären Zustände
einen Raum der Dimension mindestens $C$ auf: Jeder Sektor trägt seinen eigenen stationären Zustand.

*Beweis.* $\mathcal L$ bildet $P_c\mathcal MP_c\cong M_{k_c}$ in sich ab ($H$ und alle $L_j$ sind
blockdiagonal), also ist die Einschränkung von $e^{t\mathcal L}$ eine quantendynamische Halbgruppe
auf $M_{k_c}$ (Lemma 8.4(a)) mit einem stationären Zustand $\sigma_c$ nach Proposition 8.12. Die
$\sigma_c$ haben disjunkte Träger. $\square$

**Proposition 8.14 (spontane Emission).** Auf $M_2$ mit der Basis $|g\rangle,|e\rangle$ seien
$H=\omega|e\rangle\langle e|$, ein Sprung $L=|g\rangle\langle e|$ und die Rate $\gamma>0$ gegeben.
Dann gilt
$$
\rho_{ee}(t)=e^{-\gamma t}\rho_{ee}(0),\qquad\rho_{ge}(t)=e^{(i\omega-\gamma/2)t}\rho_{ge}(0),
\qquad e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}},
$$
wobei $\mathcal U_t(\rho)=e^{-itH}\rho e^{itH}$ und $\mathcal A_p$ der Amplitudendämpfungskanal mit
Parameter $p$ ist (Abschnitt 6.4). Das Spektrum von $\mathcal L$ ist
$\{0,-\gamma,-\gamma/2\pm i\omega\}$, $|g\rangle\langle g|$ ist der einzige stationäre Zustand, und
es gilt
$\tfrac12\|\rho(t)-|g\rangle\langle g|\|_1=\big(\rho_{ee}(0)^2e^{-2\gamma t}+|\rho_{ge}(0)|^2e^{-\gamma t}\big)^{1/2}$.

*Beweis.* Hier ist $L^*L=|e\rangle\langle e|$ und $L\rho L^*=\rho_{ee}|g\rangle\langle g|$, also
eintragsweise $\rho_{ee}'=-\gamma\rho_{ee}$, $\rho_{gg}'=\gamma\rho_{ee}$ und
$\rho_{ge}'=(i\omega-\gamma/2)\rho_{ge}$, woraus die Lösung folgt. In der Basis
$(\rho_{gg},\rho_{ge},\rho_{eg},\rho_{ee})$ ist die Matrix von $\mathcal L$ eine obere
Dreiecksmatrix mit der Diagonale $(0,-\gamma/2+i\omega,-\gamma/2-i\omega,-\gamma)$.
$\mathcal L(\sigma)=0$ erzwingt $\sigma_{ee}=\sigma_{ge}=0$. Der Amplitudendämpfungskanal
multipliziert $\rho_{ee}$ mit $1-p$ und $\rho_{ge}$ mit $\sqrt{1-p}$ (seine Kraus-Operatoren sind
$|g\rangle\langle g|+\sqrt{1-p}\,|e\rangle\langle e|$ und $\sqrt p\,|g\rangle\langle e|$), und
$\mathcal U_t$ multipliziert $\rho_{ge}$ mit $e^{i\omega t}$; beide sind TP, was $\rho_{gg}$
festlegt. Schließlich ist $\rho(t)-|g\rangle\langle g|$ spurfrei und selbstadjungiert mit der
Diagonale $(-\rho_{ee}(t),\rho_{ee}(t))$, seine Eigenwerte sind also
$\pm(\rho_{ee}(t)^2+|\rho_{ge}(t)|^2)^{1/2}$. $\square$

Die Kohärenz zerfällt mit der halben Rate der Besetzung – die Beziehung $T_2=2T_1$ aus Abschnitt 6.4
ohne reine Dephasierung; die langsamste von null verschiedene Rate, die **Spektrallücke**
$\gamma/2$, bestimmt die Relaxation. Wir prüfen das Spektrum mit `lindblad_superoperator`, den Kanal
und den Spurabstand mit `lindblad_channel` sowie Korollar 8.13 für zufällige Sprünge innerhalb der
Sektoren von $M_2\oplus M_3$:

```python
q2 = TypeIAlgebra([2], [2], precision='double')                    # basis |g>, |e>
omega, gamma = 1.3, 0.4
H_at = q2.from_blocks([torch.tensor([[0.0, 0.0], [0.0, omega]], dtype=f64)])
sigma_minus = q2.from_blocks([torch.tensor([[0.0, 1.0], [0.0, 0.0]], dtype=f64)])   # |g><e|
gen = dynamics.lindblad_superoperator(q2, H_at, [sigma_minus], rates=[gamma])[0][0]
spec = sorted(torch.linalg.eigvals(gen).tolist(), key=lambda z: (z.real, z.imag))
print([f"{z.real:+.2f}{z.imag:+.2f}i" for z in spec])
assert np.allclose(spec, [-gamma, -gamma / 2 - 1j * omega, -gamma / 2 + 1j * omega, 0.0], atol=1e-12)

psi_at = torch.tensor([0.6, 0.8], dtype=c128)                      # rho_ee(0) = 0.64, rho_ge(0) = 0.48
rho_at = DensityMatrix(q2, matrix=torch.outer(psi_at, psi_at.conj())[None, None])
ground = DensityMatrix(q2, matrix=torch.tensor([[[[1.0, 0.0], [0.0, 0.0]]]], dtype=c128))
for s in (0.5, 2.0, 10.0):
    ch = dynamics.lindblad_channel(q2, H_at, [sigma_minus], t=s, rates=[gamma])
    UA = channels.unitary_channel(dynamics.propagator(H_at, s)).compose(
        channels.amplitude_damping_channel(q2, 1 - math.exp(-gamma * s)))
    assert (ch.superoperator()[0] - UA.superoperator()[0]).abs().max() < 1e-12
    r = ch(rho_at).matrix[0, 0]
    assert abs(r[1, 1].real - 0.64 * math.exp(-gamma * s)) < 1e-12
    assert abs(r[0, 1] - 0.48 * np.exp((1j * omega - gamma / 2) * s)) < 1e-12
    exact = math.sqrt(0.64 ** 2 * math.exp(-2 * gamma * s) + 0.48 ** 2 * math.exp(-gamma * s))
    assert abs(states.trace_distance(ch(rho_at), ground).item() - exact) < 1e-12

two = TypeIAlgebra([2, 3], [2, 3], precision='double')
H_two = two.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
L_two = [two.operator(torch.randn(1, 2, 3, 3, dtype=c128)) for _ in range(2)]
kernel = [int((torch.linalg.eigvals(G_c[0]).abs() < 1e-10).sum())
          for G_c in dynamics.lindblad_superoperator(two, H_two, L_two)]
print("stationary states per sector:", kernel)
assert all(dim >= 1 for dim in kernel)                             # Corollary 8.13
```

Für diese zufälligen Sprünge findet der Code genau einen stationären Zustand pro Sektor. Kriterien
für die Eindeutigkeit – etwa ein treuer stationärer Zustand zusammen mit
$\{H,L_j,L_j^*\}'=\mathbb C1$ [A. Frigerio, *Commun. Math. Phys.* 63, 269 (1978)] – werden hier
nicht bewiesen; Aufgabe 4 zeigt, dass die Eindeutigkeit selbst in einem Faktor verletzt sein kann.

## 8.6 Womit man vergleichen sollte

QuTiP (`sesolve`, `mesolve` für die Lindblad-Gleichung mit adaptiven ODE-Lösern, `mcsolve` für
Quantentrajektorien, `steadystate`, `liouvillian`), dynamiqs (JAX; gebatchte und differenzierbare
Schrödinger- und Lindblad-Löser auf GPUs), QuantumOptics.jl und SciPy (`scipy.linalg.expm` sowie
`scipy.sparse.linalg.expm_multiply` für $e^{t\mathcal L}\rho$ mit einem dünnbesetzten vektorisierten
Generator) lösen dieselben Gleichungen auf einem einzelnen Hilbertraum, und für die obigen Probleme
müssen sie dieselben Zahlen liefern. Auf $\mathcal M$ übergibt man ihnen blockdiagonale
$D\times D$-Matrizen; der vektorisierte Generator lebt dann auf einem $D^2$-dimensionalen Raum
(dünnbesetzt, aber mit allen Sektoren in einer gemeinsamen Rechnung gekoppelt), sofern der Benutzer
das Problem nicht von Hand in Sektoren zerlegt. Die Bibliothek hat die Blockstruktur eingebaut:
exakte Propagatoren und Lindblad-Kanäle Sektor für Sektor ($k_c^2\times k_c^2$-Superoperatoren),
Sprünge zwischen Sektoren als `InterSectorChannel`s mit der Kompression aus Proposition 8.8, Batches
von Hamiltonoperatoren und Sprungoperatoren (für Mittelungen über Unordnung) und automatische
Differentiation (die letzten beiden auch in dynamiqs). Sie hat keine adaptive Schrittweitensteuerung
(RK4 arbeitet mit fester Schrittweite und ist nur näherungsweise positiv), keinen Trajektorienlöser,
und `lindblad_channel` akzeptiert keine Sprünge zwischen Sektoren; dafür sind QuTiP oder dynamiqs
auf $M_D$ die Referenz.

## Aufgaben

1. (Heisenberg-Bild) Man zeige
   $\operatorname{Tr}(\mathcal L(\rho)A)=\operatorname{Tr}(\rho\,\mathcal L^*(A))$ für den Generator
   aus Abschnitt 8.3, mit
   $\mathcal L^*(A)=i[H,A]+\sum_jg_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})+\sum_mg_m(\Phi_m^*(A)-\frac12\{\Phi_m^*(1),A\})$,
   ferner $\mathcal L^*(1)=0$ sowie $\mathcal L^*(P_c)=0$, falls es keine Sprünge zwischen Sektoren
   gibt.
2. (Exzeptioneller Punkt) Im Sektor $N=1$ aus Abschnitt 8.4 ist
   $H_{\mathrm{eff}}=\begin{pmatrix}-i\gamma_1/2&J\\J&-i\gamma_2/2\end{pmatrix}$. Man bestimme seine
   Eigenwerte und zeige, dass er für $J=|\gamma_1-\gamma_2|/4$ nicht diagonalisierbar ist.
3. $N_0$ Bosonen in $m$ Moden mögen jeweils mit der Rate $\gamma$ verloren gehen (Sprünge
   $a_1,\dots,a_m$, eingeschränkt auf $N\to N-1$). Man zeige $p_N(t)=\binom{N_0}Nq^N(1-q)^{N_0-N}$
   mit $q=e^{-\gamma t}$ für jeden Hamiltonoperator in $\mathcal M$.
4. (Reine Dephasierung) Auf $M_2$ seien $H=\omega|e\rangle\langle e|$ und
   $L=\sigma_z=\operatorname{diag}(1,-1)$ mit der Rate $\gamma$. Man zeige, dass die Besetzungen
   konstant sind, dass $\rho_{ge}(t)=e^{(i\omega-2\gamma)t}\rho_{ge}(0)$ gilt und dass die
   stationären Zustände den zweidimensionalen Raum der Diagonalmatrizen aufspannen, obwohl $C=1$
   ist.

<details>
<summary>Lösungen</summary>

1. Wegen der Zyklizität gilt $\operatorname{Tr}(L\rho L^*A)=\operatorname{Tr}(\rho L^*AL)$,
   $\operatorname{Tr}(\{L^*L,\rho\}A)=\operatorname{Tr}(\rho\{L^*L,A\})$ und
   $\operatorname{Tr}(-i[H,\rho]A)=\operatorname{Tr}(\rho\,i[H,A])$; für Sprünge zwischen Sektoren
   ist $\operatorname{Tr}(\Phi(\rho)A)=\operatorname{Tr}(\rho\Phi^*(A))$ die Definition von $\Phi^*$
   (Abschnitt 6.3). Dann ist
   $\mathcal L^*(1)=\sum_jg_j(L_j^*L_j-L_j^*L_j)+\sum_mg_m(\Phi_m^*(1)-\Phi_m^*(1))=0$, was noch
   einmal die Spurerhaltung ist. Ohne Sprünge zwischen Sektoren kommutiert $P_c$ mit $H$ und allen
   $L_j$, also ist $L_j^*P_cL_j=L_j^*L_jP_c=\frac12\{L_j^*L_j,P_c\}$ und $\mathcal L^*(P_c)=0$; dann
   ist $p_c'=\operatorname{Tr}(\rho\mathcal L^*(P_c))=0$, ein zweiter Beweis von Proposition 8.9(1).

2. $H_{\mathrm{eff}}=-\frac i4(\gamma_1+\gamma_2)1+\begin{pmatrix}-i\delta&J\\J&i\delta\end{pmatrix}$
   mit $\delta=(\gamma_1-\gamma_2)/4$; die zweite Matrix hat die Spur $0$ und die Determinante
   $\delta^2-J^2$, die Eigenwerte sind also $-\frac i4(\gamma_1+\gamma_2)\pm\sqrt{J^2-\delta^2}$:
   zwei Zerfallsraten für $J<|\delta|$, eine Rate und zwei Frequenzen für $J>|\delta|$. Für
   $J=|\delta|>0$ erfüllt die zweite Matrix $Y$ die Gleichung $Y^2=(J^2-\delta^2)1=0$, aber $Y\ne0$;
   daher hat $H_{\mathrm{eff}}$ den einzigen Eigenwert $-\frac i4(\gamma_1+\gamma_2)$ und ist kein
   Vielfaches von $1$: Er ist nicht diagonalisierbar.

```python
lam = np.linalg.eigvals(H_eff[1])
delta = (g1 - g2) / 4
root = np.sqrt(J ** 2 - delta ** 2 + 0j)
assert np.allclose(sorted(lam, key=np.real), sorted(-0.25j * (g1 + g2) + np.array([root, -root]), key=np.real))
H_ep = np.array([[-0.5j * g1, delta], [delta, -0.5j * g2]])         # J = delta: the exceptional point
Y = H_ep + 0.25j * (g1 + g2) * np.eye(2)                           # H_ep minus its only eigenvalue
assert np.abs(Y @ Y).max() < 1e-15 and np.linalg.matrix_rank(Y) == 1   # one Jordan block of size 2
```

3. Auf dem Sektor $N$ gilt $\sum_ja_j^*a_j=N1$, und das einzige Sektorpaar ist $(N-1,N)$, also
   liefert Proposition 8.9(2) $p_N'=-\gamma Np_N+\gamma(N+1)p_{N+1}$. Mit $q'=-\gamma q$ ist die
   Ableitung des Binomialterms gleich
   $-\gamma Np_N+\gamma(N_0-N)\binom{N_0}Nq^{N+1}(1-q)^{N_0-N-1}$, und
   $(N+1)\binom{N_0}{N+1}=(N_0-N)\binom{N_0}N$ zeigt, dass der zweite Term gleich
   $\gamma(N+1)p_{N+1}$ ist. Bei $t=0$ ist dies die Punktmasse in $N_0$; die Eindeutigkeit ist Lemma
   8.1.

4. $L\rho L^*-\frac12\{L^*L,\rho\}=\sigma_z\rho\sigma_z-\rho$ multipliziert die
   Nichtdiagonaleinträge mit $-2$ und annulliert die Diagonale, und
   $-i[H,\rho]_{ge}=i\omega\rho_{ge}$. Also gilt $\mathcal L(\sigma)=0$ genau dann, wenn
   $\sigma_{ge}=0$ (für $\gamma>0$), und das Spektrum ist $\{0,0,-2\gamma\pm i\omega\}$. Da $H$ und
   $\sigma_z$ kommutieren und $X\mapsto\sigma_zX\sigma_z$ eine Involution ist, gilt
   $e^{t\mathcal L}(\rho)=\frac{1+e^{-2\gamma t}}2U_t\rho U_t^* +\frac{1-e^{-2\gamma t}}2U_t\sigma_z\rho\sigma_zU_t^*$
   mit $U_t=e^{-itH}$: ein unitaler Kanal, dessen Fixpunkte die mit seinen Kraus-Operatoren
   kommutierenden Matrizen sind (Satz 6.10), also die Diagonalmatrizen.

```python
deph = dynamics.lindblad_superoperator(q2, H_at, [q2.from_blocks([torch.diag(torch.tensor([1.0, -1.0], dtype=f64))])],
                                       rates=[gamma])[0][0]
ev = torch.linalg.eigvals(deph)
assert int((ev.abs() < 1e-12).sum()) == 2
assert np.allclose(sorted(ev.tolist(), key=lambda z: (z.real, z.imag))[:2],
                   [-2 * gamma - 1j * omega, -2 * gamma + 1j * omega], atol=1e-12)
```

</details>

**Weiter:** Kapitel 9 (in Vorbereitung) behandelt Spinketten: Sektoren von Magnetisierung und
Impuls, Lanczos- und Krylov-Zeitentwicklung für große Sektoren sowie Unordnung.
