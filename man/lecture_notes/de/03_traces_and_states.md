[English](../03_traces_and_states.md) · [Español](../es/03_traces_and_states.md) · [Français](../fr/03_traces_and_states.md) · **Deutsch** · [中文](../zh/03_traces_and_states.md) · [日本語](../ja/03_traces_and_states.md) · [Русский](../ru/03_traces_and_states.md)

# 3. Spuren, Zustände und Entropie

Auf $M_n(\mathbb C)$ gibt es im Wesentlichen nur eine Spur. Auf $\bigoplus_cM_{k_c}$ gibt es
eine $C$-parametrige Familie, und die Wahl ist von Bedeutung: Dichten, Entropien und die
dualen Abbildungen hängen alle von ihr ab. Dieses Kapitel klassifiziert die Spuren,
stellt Zustände durch Dichten dar, beweist die grundlegenden Entropieungleichungen und
endet mit einer plausiblen Spurungleichung, die falsch ist.

**Sie lernen** die Klassifikation der Spuren auf $\bigoplus_cM_{k_c}$ und den Zusammenhang zwischen den Dichten
eines Zustands bezüglich verschiedener Spuren kennen, den klassischen Teil eines Zustands, die Zerlegung der
von-Neumann-Entropie nach Sektoren, die Kleinsche Ungleichung und ihre Folgerungen sowie den Michelson-Kontrast.
Die numerischen Kontrollen dieses Kapitels stehen im
[Tutorial 3](../../tutorial/de/03_traces_and_states.md).

Im Folgenden sei $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ mit allen $k_c\ge1$,
$D=\sum_ck_c$, $P_c$ die minimalen zentralen Projektionen und $A=\bigoplus_cA_c$.

<a id="s-3-1"></a>

## 3.1 Spuren

Eine **Spur** auf $\mathcal M$ ist ein lineares Funktional $\tau$ mit $\tau(AB)=\tau(BA)$ für
alle $A,B$. Sie heißt **positiv**, wenn $\tau(A^*A)\ge0$ für alle $A$ gilt, und **treu**, wenn
außerdem aus $\tau(A^*A)=0$ stets $A=0$ folgt.

<a id="r-3-1"></a>

**Satz 3.1.** Jede Spur auf $\mathcal M$ hat die Form
$\mathrm{tr}_w(A)=\sum_cw_c\mathrm{Tr}\,A_c$ mit $w\in\mathbb C^C$. Sie ist genau dann
positiv, wenn alle $w_c\ge0$, und genau dann treu, wenn alle $w_c>0$ sind.

*Beweis.* Seien $E^c_{ij}$ die Matrixeinheiten des $c$-ten Blocks. Für eine Spur $\tau$ gilt
$\tau(E^c_{ij})=\tau(E^c_{i1}E^c_{1j})=\tau(E^c_{1j}E^c_{i1})=\delta_{ij}\tau(E^c_{11})$.
Man setze $w_c=\tau(E^c_{11})$; wegen der Linearität ist $\tau(A)=\sum_c\sum_i(A_c)_{ii}w_c$.
Umgekehrt ist jedes $\mathrm{tr}_w$ eine Spur. Positivität: $E^c_{11}=(E^c_{11})^*E^c_{11}$,
also ist $w_c=\tau\big((E^c_{11})^*E^c_{11}\big)\ge0$ notwendig, und $w_c>0$ ist für Treue notwendig,
da $E^c_{11}\neq0$. Umgekehrt ist $\mathrm{tr}_w(A^*A)=\sum_cw_c\|A_c\|_2^2$, was $\ge0$ ist, wenn
alle $w_c\ge0$, und nur für $A=0$ verschwindet, wenn alle $w_c>0$. $\square$

Die Bibliothek implementiert drei treue Spuren:

| Name | Gewichte $w_c$ | Wert auf $1$ | Bedeutung |
|---|---|---|---|
| `Tr_blunt` | $1$ | $D$ | Spur von $\bigoplus_c\mathbb C^{k_c}$ (Vielfachheiten $m_c=1$) |
| `Tr_norm` | $1/k_c$ | $C$ | jeder Block normiert |
| `tau_vN` | $1/(Ck_c)$ | $1$ | der Spurzustand mit gleichen Sektorgewichten $1/C$ |

Nach Abschnitt 2.3 ist die Spur eines Raumes, auf dem $\mathcal M$ mit Vielfachheiten $m_c$
wirkt, gleich $\mathrm{tr}_w$ mit $w_c=m_c$. Eine treue Spur ist also die Festlegung,
„wie viel jeder Sektor zählt“.

*In der Bibliothek:* Die Spuren `Tr_blunt`, `Tr_norm`, `tau_vN`, ihre Werte auf $1$ und die Spureigenschaft aus Satz 3.1 — [Tutorial, Abschnitt 3.1](../../tutorial/de/03_traces_and_states.md#s-3-1).

<a id="s-3-2"></a>

## 3.2 Zustände und Dichten

Ein **Zustand** ist ein lineares Funktional $\omega$ auf $\mathcal M$, das positiv
($\omega(A^*A)\ge0$) und normiert ($\omega(1)=1$) ist.

<a id="r-3-2"></a>

**Satz 3.2.** Sei $\mathrm{tr}_w$ eine treue Spur. Zu jedem linearen Funktional
$\omega$ gibt es genau ein $\rho_w\in\mathcal M$ mit $\omega(A)=\mathrm{tr}_w(\rho_wA)$ für
alle $A$. Das Funktional $\omega$ ist genau dann positiv, wenn $\rho_w\ge0$, und genau dann
ein Zustand, wenn außerdem $\mathrm{tr}_w(\rho_w)=1$ gilt. Für zwei treue Spuren hängen
die Dichten blockweise zusammen durch
$$
w_c\,(\rho_w)_c=w'_c\,(\rho_{w'})_c .
$$

*Beweis.* $\langle X,Y\rangle_w=\mathrm{tr}_w(X^*Y)$ ist ein Skalarprodukt auf $\mathcal M$
(Satz 3.1). Nach dem Darstellungssatz von Riesz im endlichdimensionalen Raum
$\mathcal M$ gibt es genau ein $R$ mit $\omega(A)=\langle R,A\rangle_w=\mathrm{tr}_w(R^*A)$;
man setze $\rho_w=R^*$. Ist $\rho_w\ge0$, so gilt
$\omega(A^*A)=\sum_cw_c\mathrm{Tr}(A_c\,(\rho_w)_c\,A_c^*)\ge0$, da
$A_c(\rho_w)_cA_c^*\ge0$ (Korollar 1.7 mit $T=(\rho_w)_c^{1/2}A_c^*$), eine positive Matrix
Spur $\ge0$ hat (ihre Eigenwerte sind $\ge0$, Korollar 1.4(3)) und $w_c>0$. Ist umgekehrt
$\omega$ positiv, so nehme man $x\in\mathbb C^{k_c}$ und als $A$ den Rang-eins-Operator
$|x\rangle\langle x|$ im Block $c$ ($x\neq0$); dann ist $A=A^*A/\|x\|^2$ und
$0\le\omega(A)=w_c\langle x,(\rho_w)_cx\rangle$, also $(\rho_w)_c\ge0$. Die Normierung lautet
$\omega(1)=\mathrm{tr}_w(\rho_w)$. Schließlich sind sowohl $\sum_cw_c\mathrm{Tr}((\rho_w)_cA_c)$ als auch
$\sum_cw'_c\mathrm{Tr}((\rho_{w'})_cA_c)$ für alle $A$ gleich $\omega(A)$, also sind
$\bigoplus_cw_c(\rho_w)_c$ und $\bigoplus_cw'_c(\rho_{w'})_c$ beide Dichten von $\omega$ bezüglich der
ungewichteten Spur ($w_c=1$), und wegen der Eindeutigkeit gilt $w_c(\rho_w)_c=w'_c(\rho_{w'})_c$. $\square$

*In der Bibliothek:* Ein Zustand ist eine `DensityMatrix`, gespeichert durch seine ungewichtete Dichte; die Umrechnungen aus Satz 3.2 sind `DensityMatrix.density` und `DensityMatrix.from_density` — [Tutorial, Abschnitt 3.2](../../tutorial/de/03_traces_and_states.md#s-3-2).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). $\omega(A)$ ist der Erwartungswert der Observablen $A$.
> Die Dichte $\rho$ bezüglich der ungewichteten Spur ist die gewöhnliche Dichtematrix der
> Quantenmechanik auf $\bigoplus_c\mathbb C^{k_c}$. Dichten bezüglich anderer Spuren treten
> auf, wenn dieselbe Algebra mit Vielfachheiten wirkt (Kapitel 2), etwa auf dem
> symmetrischen und dem nichtsymmetrischen Teil eines Spinsystems.

<a id="s-3-3"></a>

## 3.3 Der klassische Teil eines Zustands

Die Einschränkung von $\omega$ auf das Zentrum ist eine Wahrscheinlichkeitsverteilung auf den
Sektoren, $p_c=\omega(P_c)=\mathrm{Tr}\,\rho_c$. Ist $p_c>0$, so hat der **bedingte Zustand**
im Sektor $c$ die Dichte $\sigma_c=\rho_c/p_c$ auf $\mathbb C^{k_c}$, und es gilt
$$
\rho=\bigoplus_cp_c\,\sigma_c
$$
(ist $p_c=0$, so ist $\rho_c=0$, da die Eigenwerte von $\rho_c\ge0$ nichtnegativ sind und die
Summe $\mathrm{Tr}\,\rho_c=0$ haben, Korollar 1.4(3), also alle $0$ sind, und $\rho_c$ unitär
diagonalisierbar ist, Satz 1.3; der Term entfällt).
Also ist jeder Zustand auf $\mathcal M$ eine Konvexkombination von Zuständen der Sektoren. Die
Abbildung $\rho\mapsto\sum_cP_c\rho P_c$ ist die Identität auf $\mathcal M$, da jedes Element von
$\mathcal M$ blockdiagonal ist.

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Ein Zustand auf
> $\mathcal M$ ist eine klassische Mischung von Quantenzuständen in den Sektoren. Die Abbildung
> $\rho\mapsto\sum_cP_c\rho P_c$ ist die Messung von „welcher Sektor“, und sie stört keinen
> Zustand von $\mathcal M$. Das ist der operationale Gehalt einer Superauswahlregel.

*Rechnung:* Die Sektorwahrscheinlichkeiten $p_c$ (`DensityMatrix.sector_probabilities`) und die Messung $\rho\mapsto\sum_cP_c\rho P_c$ (`channels.measurement_channel`), die $\rho$ nicht stört — [Tutorial, Abschnitt 3.3](../../tutorial/de/03_traces_and_states.md#s-3-3).

<a id="s-3-4"></a>

## 3.4 Entropie

Die **von-Neumann-Entropie** eines Zustands mit (ungewichteter) Dichte $\rho$ ist
$S(\rho)=-\mathrm{Tr}\,\rho\log\rho=-\sum_i\lambda_i\log\lambda_i$, summiert über die
Eigenwerte von $\rho$, mit $0\log0=0$; d. h. $\rho\log\rho:=f(\rho)$ für $f(t)=t\log t$, $f(0)=0$, und
$\mathrm{Tr}\,f(\rho)=\sum_if(\lambda_i)$ nach der Definition von $f(\rho)=U\,\mathrm{diag}(f(\lambda_i))\,U^*$
(Abschnitt 1.3).

<a id="r-3-3"></a>

**Satz 3.3 (Entropie einer Mischung von Sektoren).** Mit $p_c$ und $\sigma_c$ wie oben
(Sektoren mit $p_c=0$ werden weggelassen) gilt
$$
S(\rho)=H(p)+\sum_cp_c\,S(\sigma_c),\qquad H(p)=-\sum_cp_c\log p_c .
$$

*Beweis.* Die von null verschiedenen Eigenwerte von $\rho$ (mit Vielfachheit) sind die Zahlen $p_c\mu$ zu den
von null verschiedenen Eigenwerten $\mu$ der $\sigma_c$ (Sektoren mit $p_c>0$); Eigenwerte
null tragen zu beiden Seiten $0$ bei. Folglich ist
$-\sum_c\sum_{\mu>0}p_c\mu\log(p_c\mu)=-\sum_cp_c\log p_c\sum_{\mu>0}\mu-\sum_cp_c\sum_{\mu>0}\mu\log\mu$,
und $\sum_{\mu>0}\mu=\mathrm{Tr}\,\sigma_c=1$. $\square$

*Rechnung:* Satz 3.3 für zufällige Zustände, mit `DensityMatrix.von_neumann_entropy` — [Tutorial, Abschnitt 3.4](../../tutorial/de/03_traces_and_states.md#s-3-4).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Die Formel zerlegt die Unsicherheit über das System in zwei Teile:
> $H(p)$, die klassische Unsicherheit über die Ladung (in welchem Sektor sich das System befindet),
> und $\sum_cp_cS(\sigma_c)$, die mittlere Quantenunsicherheit innerhalb der Sektoren. Unter einer
> Superauswahlregel ist der erste Teil die Entropie einer klassischen Zufallsvariablen: der Ladung.

Die Entropie bezüglich einer anderen Spur ist $S_w(\rho)=-\mathrm{tr}_w(\rho_w\log\rho_w)$,
mit derselben Konvention $\rho_w\log\rho_w:=f(\rho_w)$. Nach Satz 3.2 ist $(\rho_w)_c=\rho_c/w_c$, also
$$
S_w(\rho)=-\sum_cw_c\mathrm{Tr}\Big(\frac{\rho_c}{w_c}\log\frac{\rho_c}{w_c}\Big)
=S(\rho)+\sum_cp_c\log w_c .
$$
Die Entropie ist eine Eigenschaft des Zustands *und* der Spur: Derselbe Zustand hat
bezüglich `Tr_blunt` und `tau_vN` verschiedene Entropien.

*Rechnung:* Die Entropie bezüglich $\tau_{\mathrm{vN}}$ und die Formel für $S_w(\rho)$ — [Tutorial, Abschnitt 3.4](../../tutorial/de/03_traces_and_states.md#s-3-4).

<a id="s-3-5"></a>

## 3.5 Relative Entropie und die Kleinsche Ungleichung

Für Zustände mit Dichten $\rho,\sigma$, $\sigma>0$, ist die **relative Entropie**
$S(\rho\|\sigma)=\mathrm{Tr}\,\rho\log\rho-\mathrm{Tr}\,\rho\log\sigma$, wobei
$\mathrm{Tr}\,\rho\log\rho:=\sum_i\lambda_i\log\lambda_i$, summiert über die Eigenwerte von
$\rho$, mit $0\log0=0$ (und ebenso in jedem Block). Sie hängt nicht von der Spur ab: Nach Satz 3.2
ist $(\rho_w)_c=\rho_c/w_c$ und $(\sigma_w)_c=\sigma_c/w_c$, die Eigenwerte von $(\rho_w)_c$ sind also
die von $\rho_c$, geteilt durch $w_c$, es gilt $\log(\sigma_c/w_c)=\log\sigma_c-(\log w_c)1$ (ist
$\sigma_c=U\,\mathrm{diag}(s_i)\,U^*$, so ist $\sigma_c/w_c=U\,\mathrm{diag}(s_i/w_c)\,U^*$, und Satz 1.5(a)),
und sowohl $\mathrm{tr}_w(\rho_w\log\rho_w)$ als auch
$\mathrm{tr}_w(\rho_w\log\sigma_w)$ unterscheiden sich von $\mathrm{Tr}\,\rho\log\rho$ bzw.
$\mathrm{Tr}\,\rho\log\sigma$ um denselben Betrag $-\sum_cp_c\log w_c$ (der Logarithmus von $\rho$
wird nach der Konvention $0\log0=0$ nur auf dessen Träger gebraucht). Die Differenz bleibt
unverändert.

<a id="r-3-4"></a>

**Lemma 3.4 (Tangentenungleichung).** Sei $f:I\to\mathbb R$ konvex auf einem Intervall $I$ und
differenzierbar in $b\in I$ (einseitig, wenn $b$ ein Randpunkt von $I$ ist). Dann gilt
$f(a)\ge f(b)+(a-b)f'(b)$ für alle $a\in I$, und zwar strikt für $a\neq b$, wenn $f$ strikt
konvex ist, d. h. wenn $f((1-t)x+ty)<(1-t)f(x)+tf(y)$ für $x\neq y$ in $I$ und $t\in(0,1)$.
Ist $f$ außerdem stetig auf $I$ und im Inneren von $I$ zweimal differenzierbar mit $f''>0$, so
ist $f$ strikt konvex.

*Beweis.* Für $t\in(0,1]$ liefert die Konvexität $f(b+t(a-b))\le(1-t)f(b)+tf(a)$, d. h.
$\big(f(b+t(a-b))-f(b)\big)/t\le f(a)-f(b)$; der Grenzübergang $t\to0$ ergibt die Ungleichung.
Ist $f$ strikt konvex und $a\neq b$, so sei $m=(a+b)/2$; dann gilt
$f(a)-f(b)>2\big(f(m)-f(b)\big)\ge2(m-b)f'(b)=(a-b)f'(b)$, der erste Schritt wegen strikter
Konvexität in $m$, der zweite nach der eben bewiesenen Ungleichung, angewandt auf $m$. Für die
letzte Behauptung seien $x<y$ in $I$, $t\in(0,1)$ und $z=(1-t)x+ty$. Nach dem Mittelwertsatz gibt
es $\xi_1\in(x,z)$ und $\xi_2\in(z,y)$ mit $f(z)-f(x)=t(y-x)f'(\xi_1)$ und
$f(y)-f(z)=(1-t)(y-x)f'(\xi_2)$, und $f'(\xi_1)<f'(\xi_2)$ nach dem Mittelwertsatz für $f'$ auf
$[\xi_1,\xi_2]$, da dort $f''>0$. Multiplikation mit $t(1-t)(y-x)>0$ ergibt
$(1-t)\big(f(z)-f(x)\big)<t\big(f(y)-f(z)\big)$, d. h. $f(z)<(1-t)f(x)+tf(y)$. $\square$

<a id="r-3-5"></a>

**Satz 3.5 (Kleinsche Ungleichung).** Sei $f:I\to\mathbb R$ konvex auf einem Intervall $I$,
und seien $X,Y$ selbstadjungiert mit Spektren in $I$, wobei $f$ in jedem Eigenwert von $Y$
differenzierbar ist (einseitig in einem Randpunkt von $I$). Dann gilt $\mathrm{Tr}[f(X)-f(Y)-(X-Y)f'(Y)]\ge0$. Ist $f$ strikt
konvex, so gilt Gleichheit nur für $X=Y$.

*Beweis.* Seien $X=\sum_ia_i|u_i\rangle\langle u_i|$ und $Y=\sum_jb_j|v_j\rangle\langle v_j|$ mit
Orthonormalbasen aus Eigenvektoren, und sei $c_{ij}=|\langle u_i,v_j\rangle|^2$, also
$\sum_ic_{ij}=\sum_jc_{ij}=1$. Dann ist $\mathrm{Tr}\,f(X)=\sum_{ij}c_{ij}f(a_i)$,
$\mathrm{Tr}\,f(Y)=\sum_{ij}c_{ij}f(b_j)$ und
$\mathrm{Tr}\,Xf'(Y)=\sum_{ij}c_{ij}a_if'(b_j)$, $\mathrm{Tr}\,Yf'(Y)=\sum_{ij}c_{ij}b_jf'(b_j)$.
Folglich ist die linke Seite gleich $\sum_{ij}c_{ij}[f(a_i)-f(b_j)-(a_i-b_j)f'(b_j)]$, nach
Lemma 3.4 eine Summe nichtnegativer Terme. Ist $f$ strikt konvex und die
Summe null, so gilt $a_i=b_j$, sobald $c_{ij}>0$, also
$Xv_j=\sum_ia_i\langle u_i,v_j\rangle u_i=b_jv_j$ für jedes $j$, und $X=Y$. $\square$

<a id="r-3-6"></a>

**Korollar 3.6.** $S(\rho\|\sigma)\ge0$, mit Gleichheit genau dann, wenn $\rho=\sigma$.
Insbesondere ist $S(\rho)\le\log D$, mit Gleichheit genau dann, wenn $\rho=1/D$.

*Beweis.* Man wende Satz 3.5 auf $f(t)=t\log t$, $f(0)=0$, an – eine Funktion, die auf
$I=[0,\infty)$ nach Lemma 3.4 strikt konvex ist ($f$ ist stetig in $0$, da $t\log t\to0$ für
$t\to0^+$, und $f''(t)=1/t>0$ auf $(0,\infty)$)
und auf $(0,\infty)\supseteq\mathrm{spec}\,\sigma$ differenzierbar mit $f'(t)=\log t+1$ –,
mit $X=\rho$, $Y=\sigma$: Die linke Seite ist
$\mathrm{Tr}\,\rho\log\rho-\mathrm{Tr}\,\sigma\log\sigma-\mathrm{Tr}(\rho-\sigma)(\log\sigma+1)
=S(\rho\|\sigma)$, da $\mathrm{Tr}\,\rho=\mathrm{Tr}\,\sigma=1$. Für die zweite Behauptung
nehme man $\sigma=1/D$: $0\le S(\rho\|1/D)=-S(\rho)+\log D$. $\square$

*Rechnung:* $S(\rho)\le\log D$ und der maximal gemischte Zustand — [Tutorial, Abschnitt 3.5](../../tutorial/de/03_traces_and_states.md#s-3-5).

<a id="r-3-7"></a>

**Satz 3.7 (Kettenregel; Monotonie unter der bedingten Erwartung auf das Zentrum).** Seien $\rho,\sigma$ Zustände mit
$\sigma>0$, $p_c=\mathrm{Tr}\,\rho_c$, $q_c=\mathrm{Tr}\,\sigma_c>0$, $\hat\sigma_c=\sigma_c/q_c$ und,
für $p_c>0$, $\hat\rho_c=\rho_c/p_c$. Dann gilt, mit Summen über die $c$ mit $p_c>0$ (die
übrigen Sektoren tragen wegen $0\log0=0$ den Wert $0$ bei),
$$
S(\rho\|\sigma)=\sum_cp_c\log\frac{p_c}{q_c}+\sum_cp_c\,S(\hat\rho_c\|\hat\sigma_c)
\ \ge\ \sum_cp_c\log\frac{p_c}{q_c}=S(E\rho\|E\sigma),
$$
wobei $E$ die bedingte Erwartung auf das Zentrum ist (Proposition 2.8).

*Beweis.* Für $p_c>0$ ist $\mathrm{Tr}\,\rho_c\log\rho_c=p_c\log p_c+p_c\mathrm{Tr}\,\hat\rho_c\log\hat\rho_c$
(Eigenwerte $p_c\mu$, wie in Satz 3.3) und $\log\sigma_c=(\log q_c)1+\log\hat\sigma_c$ (man
diagonalisiere $\hat\sigma_c=U\,\mathrm{diag}(s_i)\,U^*$, also $\sigma_c=U\,\mathrm{diag}(q_cs_i)\,U^*$, und
verwende Satz 1.5(a)), also
$\mathrm{Tr}\,\rho_c\log\sigma_c=p_c\log q_c+p_c\mathrm{Tr}\,\hat\rho_c\log\hat\sigma_c$ und
$\mathrm{Tr}\,\rho_c\log\rho_c-\mathrm{Tr}\,\rho_c\log\sigma_c=p_c\log(p_c/q_c)+p_cS(\hat\rho_c\|\hat\sigma_c)$.
Ist $p_c=0$, so ist $\rho_c=0$, und der Sektor trägt $0$ bei.
Man summiere über $c$; die zweite Summe ist nach Korollar 3.6, angewandt in $M_{k_c}$ auf die
Zustände $\hat\rho_c$ und $\hat\sigma_c>0$, $\ge0$. Schließlich sind
$E\rho=\bigoplus_c(p_c/k_c)1_{k_c}$ und $E\sigma=\bigoplus_c(q_c/k_c)1_{k_c}>0$ Zustände; $E\rho$ ist die
Dichte des Zustands $\omega\circ E$ (mit $\omega(A)=\mathrm{Tr}\,\rho A$), da
$\mathrm{Tr}(\rho E(A))=\sum_cp_c\mathrm{Tr}(A_c)/k_c=\mathrm{Tr}(E(\rho)A)$. Und
$S(E\rho\|E\sigma)=\sum_ck_c\frac{p_c}{k_c}\big(\log\frac{p_c}{k_c}-\log\frac{q_c}{k_c}\big)=\sum_cp_c\log\frac{p_c}{q_c}$,
da die Eigenwerte von $E\rho$ im Block $c$ gleich $p_c/k_c$ sind ($k_c$-fach) und Sektoren mit
$p_c=0$ den Wert $0$ beitragen. $\square$

*Rechnung:* Satz 3.7 für zufällige Zustände, mit `states.relative_entropy` und `channels.center_expectation` — [Tutorial, Abschnitt 3.5](../../tutorial/de/03_traces_and_states.md#s-3-5).

Ersetzt man beide Zustände durch ihre Bilder unter $E$, so kann die relative Entropie also nur
abnehmen; was bleibt, ist die Kullback-Leibler-Divergenz $\sum_cp_c\log(p_c/q_c)$ der
Sektorverteilungen.

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Die relative
> Entropie misst, wie gut sich zwei Zustände unterscheiden lassen. Das Vergessen der
> Quanteninformation innerhalb der Sektoren kann zwei Zustände nur schwerer unterscheidbar
> machen; was bleibt, ist die klassische Unterscheidbarkeit der Ladungsverteilungen.

<a id="s-3-6"></a>

## 3.6 Der Michelson-Kontrast

Für $A\ge0$, $A\neq0$, ist der **Michelson-Kontrast**
$\Delta(A)=(\lambda_{\max}-\lambda_{\min})/(\lambda_{\max}+\lambda_{\min})\in[0,1]$, wobei die
extremalen Eigenwerte über alle Sektoren genommen werden (und $\Delta(0)=0$). Er misst, wie
weit $A$ von einem Vielfachen der Identität entfernt ist: $\Delta(A)=0$ genau dann, wenn $A=\lambda1$.

<a id="r-3-8"></a>

**Proposition 3.8.** Ist $\Phi:\mathcal M\to\mathcal M$ linear, positiv und unital
($\Phi(1)=1$), so gilt $\Delta(\Phi(A))\le\Delta(A)$ für alle $A\ge0$.

*Beweis.* Es ist $\lambda_{\min}1\le A\le\lambda_{\max}1$ (Korollar 1.4(3) für $A-\lambda_{\min}1$
und $\lambda_{\max}1-A$), und Positivität, Linearität und Unitalität liefern
$\lambda_{\min}1\le\Phi(A)\le\lambda_{\max}1$. Folglich ist $\Phi(A)$ selbstadjungiert (wie jede
positive Matrix, Korollar 1.4(3)), $\Phi(A)\ge0$ (da $\lambda_{\min}\ge0$), und nach Korollar 1.4(3)
ist $\lambda_{\min}(\Phi(A))\ge\lambda_{\min}$ und $\lambda_{\max}(\Phi(A))\le\lambda_{\max}$. Ist
$\Phi(A)=0$, so ist die Behauptung trivial. Andernfalls ist $A\neq0$, sowohl $\lambda_{\max}$ als
auch $\lambda_{\max}(\Phi(A))$ sind $>0$, und mit $r=\lambda_{\min}/\lambda_{\max}\in[0,1]$ ist
die Größe $\Delta=(1-r)/(1+r)$ fallend in $r$, und $r(\Phi(A))\ge r(A)$. $\square$

*Rechnung:* Proposition 3.8 für einen zufälligen gemischt-unitären Kanal und für $E$, mit `Operator.michelson_contrast` — [Tutorial, Abschnitt 3.6](../../tutorial/de/03_traces_and_states.md#s-3-6).

<a id="s-3-7"></a>

## 3.7 Eine plausible Ungleichung, die falsch ist

Gilt für positive $X,Y$ in $M_k$ und ein unitäres $U$ die Ungleichung
$|\mathrm{Tr}(XUY)|\le\mathrm{Tr}(XY)$? Sie gilt für $U=1$ und sieht wie eine
Aussage vom Cauchy-Schwarz-Typ aus. Sie ist jedoch falsch, sobald $X$ und $Y$ nicht kommutieren.

<a id="r-3-9"></a>

**Proposition 3.9.** Seien $X,Y\ge0$ in $M_k$. Dann gilt $\max_U|\mathrm{Tr}(XUY)|=\|YX\|_1\ge\mathrm{Tr}(XY)$, wobei das Maximum
über die unitären $U\in M_k$ läuft, mit Gleichheit genau dann, wenn $XY=YX$.

*Beweis.* $\mathrm{Tr}(XUY)=\mathrm{Tr}(YX\,U)$, also ist das Maximum nach Satz 1.9
gleich $\|YX\|_1$, und $\|YX\|_1\ge|\mathrm{Tr}(YX)|=\mathrm{Tr}(XY)$, da
$\mathrm{Tr}(XY)=\mathrm{Tr}(X^{1/2}YX^{1/2})\ge0$, weil
$X^{1/2}YX^{1/2}=(Y^{1/2}X^{1/2})^*(Y^{1/2}X^{1/2})\ge0$ (Korollar 1.7). Ist $XY=YX$, so
kommutiert $X^{1/2}$ mit $Y$ (Satz 1.5(d)), also ist $YX=X^{1/2}YX^{1/2}\ge0$ und
$\|YX\|_1=\mathrm{Tr}(YX)$ (für $P\ge0$ ist $\|P\|_1=\mathrm{Tr}\,P$, Abschnitt 1.4).
Umgekehrt gelte $\|YX\|_1=\mathrm{Tr}(YX)$. Ist $|YX|=0$, so ist $YX=0=(YX)^*=XY$.
Andernfalls schreibe man $YX=W|YX|$ (Satz 1.8) und $R=|YX|^{1/2}\neq0$. Wegen der Zyklizität ist
$\mathrm{Tr}(YX)=\mathrm{Tr}(RWR)=\langle R,WR\rangle_{HS}$, während
$\|R\|_2^2=\|WR\|_2^2=\mathrm{Tr}|YX|=\|YX\|_1$. Also ist $\langle R,WR\rangle_{HS}=\|R\|_2\|WR\|_2$:
In der Cauchy-Schwarz-Ungleichung gilt Gleichheit, folglich ist $WR=\lambda R$ für einen Skalar
$\lambda$, und $\lambda\|R\|_2^2=\langle R,WR\rangle_{HS}=\|R\|_2^2$ ergibt $\lambda=1$. Also ist
$W|YX|^{1/2}=|YX|^{1/2}$ und $YX=W|YX|=|YX|\ge0$. Eine positive Matrix ist selbstadjungiert, also
$YX=(YX)^*=XY$. $\square$

Nach Satz 1.9 erreicht für $YX=W|YX|$ das unitäre $U=W^*$ das Maximum $\|YX\|_1$. Eine zufällige unitäre Matrix
ist im Mittel schlechter: Für ein Haar-verteiltes $U$ (Satz 5.3) liefert
Proposition 11.11 mit $A=YX$
$\big(\mathbb E|\mathrm{Tr}(XUY)|^2\big)^{1/2}=\|YX\|_2/\sqrt k\le\|YX\|_1/\sqrt k=\max_V|\mathrm{Tr}(XVY)|/\sqrt k$
(da $\sum_is_i^2\le(\sum_is_i)^2$), und Korollar 11.12 zeigt, dass für $X,Y$ mit Spektren in
$[a,b]$, $a>0$, eine einzelne Stichprobe $\mathrm{Tr}(XY)$ höchstens mit Wahrscheinlichkeit
$b^4/(a^4k^2)$ erreicht (aussagekräftig nur für $k>b^2/a^2$).

*Rechnung:* Proposition 3.9 für ein nicht kommutierendes Paar in $M_6$, die Maximalstelle $W^*$ und eine Zufallssuche über Haar-verteilte unitäre Matrizen — [Tutorial, Abschnitt 3.7](../../tutorial/de/03_traces_and_states.md#s-3-7).

<a id="s-3-8"></a>

## 3.8 Womit man vergleichen sollte

*In der Bibliothek:* Vergleich mit anderer Software — [Tutorial, Abschnitt 3.8](../../tutorial/de/03_traces_and_states.md#s-3-8).

<a id="x-3"></a>

## Aufgaben

1. Man zeige, dass die Dichte des Spurzustands $\tau_{\mathrm{vN}}$ (gleiche Sektorgewichte)
   bezüglich $\tau_{\mathrm{vN}}$ gleich $1$ ist, und berechne seine ungewichtete Dichte und
   seine Entropie $S$.
2. Sei $\mathrm{tr}_w$ treu. Man zeige $S_w(\rho)\le\log\mathrm{tr}_w(1)$, mit Gleichheit
   genau dann, wenn $\rho_w=1/\mathrm{tr}_w(1)$.
3. Man zeige, dass der Messkanal $\rho\mapsto\sum_cP_c\rho P_c$ auf $M_D$ (nicht auf
   $\mathcal M$!) Zustände sehr wohl stört, und berechne $S(\rho)$ vorher und nachher für
   eine Superposition zweier Sektoren.

<details>
<summary>Lösungen</summary>

1. $\tau_{\mathrm{vN}}(A)=\tau_{\mathrm{vN}}(1\cdot A)$, also $\rho_{\tau}=1$ wegen der
   Eindeutigkeit. Nach Satz 3.2 ist die ungewichtete Dichte $\bigoplus_c1_{k_c}/(Ck_c)$, mit
   $p_c=1/C$ und $\sigma_c=1/k_c$, also $S=\log C+\frac1C\sum_c\log k_c$ nach Satz 3.3.

2. Summiert wird über die $c$ mit $p_c>0$. Dann gilt
   $S_w(\rho)=S(\rho)+\sum_cp_c\log w_c=-\sum_cp_c\log p_c+\sum_cp_cS(\sigma_c)+\sum_cp_c\log w_c
   \le\sum_cp_c\log\frac{k_cw_c}{p_c}\le\log\sum_{p_c>0}k_cw_c\le\log\mathrm{tr}_w(1)$, wobei
   $S(\sigma_c)\le\log k_c$ (Korollar 3.6 in $M_{k_c}$) und die endliche Jensensche Ungleichung
   $\sum_{c=1}^mp_c\log x_c\le\log\sum_{c=1}^mp_cx_c$ ($p_c>0$, $\sum_cp_c=1$, $x_c>0$), mit
   Gleichheit nur, wenn alle $x_c$ gleich sind, angewandt auf $x_c=k_cw_c/p_c$, verwendet wurden.
   Beweis durch Induktion nach $m$: Für $m=1$ ist nichts zu zeigen. Für $m\ge2$ setze man
   $q=1-p_m\in(0,1)$ und $\bar x=\sum_{c<m}(p_c/q)x_c$; nach Induktionsvoraussetzung ist
   $\sum_{c<m}(p_c/q)\log x_c\le\log\bar x$, mit Gleichheit nur, wenn $x_1=\dots=x_{m-1}$, und wegen
   der strikten Konkavität von $\log$ (Lemma 3.4 für $-\log$, dessen zweite Ableitung $1/t^2>0$ ist)
   gilt $q\log\bar x+p_m\log x_m\le\log(q\bar x+p_mx_m)$,
   mit Gleichheit nur, wenn $\bar x=x_m$. Multipliziert man die erste Ungleichung mit $q$ und
   addiert, so folgt $\sum_cp_c\log x_c\le\log\sum_cp_cx_c$, und Gleichheit erzwingt
   $x_1=\dots=x_{m-1}=\bar x=x_m$. Gleichheit überall erzwingt, dass jedes $p_c>0$ ist (da alle
   $k_cw_c>0$), $\sigma_c=1/k_c$ und
   $k_cw_c/p_c$ konstant, d. h. $p_c=k_cw_c/\mathrm{tr}_w(1)$, d. h.
   $(\rho_w)_c=\rho_c/w_c=1_{k_c}/\mathrm{tr}_w(1)$.
   Umgekehrt gilt für $\rho_w=1/\mathrm{tr}_w(1)$
   $S_w(\rho)=-\mathrm{tr}_w(1)\cdot\frac1{\mathrm{tr}_w(1)}\log\frac1{\mathrm{tr}_w(1)}=\log\mathrm{tr}_w(1)$.

3. In $M_D$ hat eine reine Superposition $\psi=(\psi_1\oplus\psi_2)/\sqrt2$ von
   Einheitsvektoren in zwei Sektoren $S=0$; nach der Messung ist der Zustand
   $\frac12|\psi_1\rangle\langle\psi_1|\oplus\frac12|\psi_2\rangle\langle\psi_2|$ mit $S=\log2$.
   Eingeschränkt auf $\mathcal M$ stimmt der Vektorzustand von $\psi$ jedoch mit der Mischung
   überein:
   $\langle\psi,A\psi\rangle=\tfrac12\langle\psi_1,A_1\psi_1\rangle+\tfrac12\langle\psi_2,A_2\psi_2\rangle$
   für alle $A\in\mathcal M$. Kein Element von $\mathcal M$ erkennt die Kohärenz (physikalisch
   ist genau das der Gehalt einer Superauswahlregel).

</details>

*Rechnung:* Aufgabe 1 numerisch — [Tutorial, Aufgaben](../../tutorial/de/03_traces_and_states.md#x-3).

**Weiter:** [Kapitel 4](04_symmetries.md) behandelt Symmetrien und Darstellungen: wie eine
Gruppenwirkung auf einem physikalischen System die Sektoren hervorbringt und wie die
Bibliothek die zugehörigen Algebren für Spinsysteme und Qubit-Register konstruiert.
