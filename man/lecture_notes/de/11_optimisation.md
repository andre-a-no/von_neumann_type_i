[English](../11_optimisation.md) · [Español](../es/11_optimisation.md) · [Français](../fr/11_optimisation.md) · **Deutsch** · [中文](../zh/11_optimisation.md) · [日本語](../ja/11_optimisation.md) · [Русский](../ru/11_optimisation.md)

# 11. Optimierung und Spurungleichungen

Eine Spurungleichung ist eine Aussage über die Extremwerte eines Funktionals von Operatoren. Abschnitt 3.7 zeigte
eine, die plausibel und falsch ist (Proposition 3.9); ein Stichprobenlauf ([Tutorial, Abschnitt 3.7](../../tutorial/de/03_traces_and_states.md#s-3-7)) hat die Verletzung numerisch nicht gefunden. Dieses Kapitel erklärt,
warum Stichproben in hoher Dimension versagen, und ersetzt sie durch Optimierung über exakt parametrisierte
Nebenbedingungsmengen; wo die Extrema in geschlossener Form bekannt sind, beweisen wir sie, sodass der Optimierer an
einem Satz geprüft wird.

**Sie lernen**, wie man unitäre Matrizen und Spektren mit vorgegebenem Kontrast parametrisiert und warum
jeder zulässige Operator erreicht wird; die von Neumannsche Spurungleichung und die exakten Extrema, die sie liefert;
eine Spurungleichung mit ihrer Gleichheitsbedingung; warum Haar-Stichproben im quadratischen Mittel weit von den
Extrema entfernt sind; und was eine Multistart-Suche beweist.
Die numerischen Kontrollen dieses Kapitels und der Optimierer der Bibliothek stehen im
[Tutorial 11](../../tutorial/de/11_optimisation.md).

Im Folgenden sei $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$ mit $`D=\sum_ck_c\ge2`$, $`\mathrm{Tr}`$ die
ungewichtete Spur (Abschnitt 3.1) und $`\|A\|_1=\sum_c\|A_c\|_1`$ die Spurnorm von $`A`$ als blockdiagonale Matrix in
$`M_D(\mathbb C)`$. Die Singulärwerte von $`A\in M_k(\mathbb C)`$ sind $`s_1(A)\ge\dots\ge s_k(A)`$; $`x^\downarrow`$ und
$`x^\uparrow`$ bezeichnen die fallende und die steigende Umordnung eines reellen Vektors $`x`$. $`\Delta`$ ist der
Michelson-Kontrast (Abschnitt 3.6), und für $`\delta\in[0,1)`$ sei
```math
\ell(\delta)=\frac{1-\delta}{1+\delta},\qquad 1-\ell(\delta)=\frac{2\delta}{1+\delta};
```
ein positiver Operator mit $`\lambda_{\max}=1`$ hat den Kontrast $`\ell(\lambda_{\min})`$, also genau dann den Kontrast
$`\delta`$, wenn $`\lambda_{\min}=\ell(\delta)`$, da $`\ell`$ eine streng fallende Involution von $`[0,1]`$ ist. Die Beispiele
sind die drei Funktionale der Begleitarbeit (Abschnitt *From sampling to optimisation*):

| | Funktional | Variablen |
|---|---|---|
| Exp. 1 | $`z=\lvert\mathrm{Tr}(XUY)\rvert-\mathrm{Tr}(XY)`$ | $`X,Y\ge0`$, $`U\in\mathcal M`$ unitär |
| Exp. 2 | $`z=\mathrm{Tr}\lvert XY\rvert-\mathrm{Tr}(\lvert X\rvert\,\lvert Y\rvert)`$ | $`X=UX_0`$, $`Y=VY_0`$, $`X_0,Y_0\ge0`$, $`U,V`$ unitär |
| Exp. 3 | $`z=\mathrm{Tr}(Y\lvert X\rvert Y)-\mathrm{Tr}\lvert YXY\rvert`$ | $`X=X^*`$, $`Y\ge0`$ |

wobei $`|X|`$, $`|Y|`$ auf $`\lambda_{\max}=1`$ normiert sind und ihre Kontraste als Parameter dienen. Die
**Einhüllenden** sind $`\overline z(\delta)=\sup\{z:\Delta(|X|)=\delta\}`$ und
$`\underline z(\delta)=\inf\{z:\Delta(|X|)=\delta\}`$, wobei alle übrigen Variablen frei sind.

<a id="s-11-1"></a>

## 11.1 Nebenbedingungsmengen und ihre Parametrisierungen

**Definition.** Für $`\delta\in[0,1)`$ und $`\ell=\ell(\delta)`$ ist die **Kontrastschale**
```math
\mathcal P_\delta=\{X\in\mathcal M:X\ge0,\ \lambda_{\max}(X)=1,\ \Delta(X)=\delta\}
=\{X=X^*:\mathrm{spec}\,X\subseteq[\ell,1],\ 1,\ell\in\mathrm{spec}\,X\}.
```
Die **Plätze** von $`\mathcal M`$ sind die $`D`$ Paare $`(c,i)`$, $`i\le k_c`$; für $`\lambda\in\mathbb R^D`$ sei
$`\mathrm{diag}(\lambda)=\bigoplus_c\mathrm{diag}(\lambda_{c,1},\dots,\lambda_{c,k_c})`$.

**Die Parametrisierung.** Für einen vorgegebenen Kontrast $`\delta`$ zieht eine Suche pro Start zwei verschiedene Plätze
$`s_\pm=(c_\pm,i_\pm)`$ ($`s_+`$ gleichverteilt unter den $`D`$ Plätzen, dann $`s_-`$ gleichverteilt unter den $`D-1`$ übrigen,
sodass jedes geordnete Paar verschiedener Plätze die Wahrscheinlichkeit $`1/(D(D-1))`$ hat) und eine Haar-verteilte
unitäre Matrix $`U_0\in\mathcal M`$ (Lemma 11.10) und bildet $`\theta\in\mathbb R^D`$ und
$`G=\bigoplus_cG_c`$ ab auf
```math
X(\theta,G)=U\mathrm{diag}(\lambda)U^*,\qquad U=U_0\,e^{G-G^*},\qquad
\lambda_{s_+}=1,\quad\lambda_{s_-}=\ell,\quad\lambda_s=\ell+(1-\ell)\tfrac{1+\sin\theta_s}2\ \text{otherwise}.
```
Unitäre Matrizen allein werden durch $`G\mapsto U_0e^{G-G^*}`$ parametrisiert. Anders als eine Sigmoidfunktion erreicht
$`(1+\sin\theta)/2`$ die Werte $`0`$ und $`1`$ bei endlichem $`\theta`$; das ist bequem, weil die extremalen Spektren dieses
Kapitels (Korollar 11.7, Proposition 11.14) bei $`\ell`$ und $`1`$ liegen. Ein freier Kontrast in $`[a,b]`$ ist ein weiterer
Parameter, $`\delta=a+(b-a)(1+\sin\eta)/2`$, etwa mit $`[a,b]=[0,0.999]`$. Ein Wert, der mit dieser Einschränkung angenommen wird, wird auch ohne sie angenommen und ist
daher weiterhin eine untere Schranke für ein uneingeschränktes Supremum (eine obere Schranke für ein
uneingeschränktes Infimum).

*In der Bibliothek:* Diese Parametrisierung ist `optimize.PositiveParam`, und $`G\mapsto U_0e^{G-G^*}`$ ist `optimize.UnitaryParam` — [Tutorial, Abschnitt 11.1](../../tutorial/de/11_optimisation.md#s-11-1).

<a id="r-11-1"></a>

**Lemma 11.1.** Für jedes unitäre $`U_0\in\mathcal M`$ bildet $`G\mapsto U_0e^{G-G^*}`$ surjektiv auf die unitäre Gruppe
von $`\mathcal M`$ ab.

*Beweis.* Sei $`U`$ unitär. Nach Satz 1.3 (blockweise) und Korollar 1.4(2) ist
$`U_0^*U=\bigoplus_cW_c\mathrm{diag}(e^{i\varphi_{c,j}})_jW_c^*`$ mit unitären $`W_c`$ und $`\varphi_{c,j}\in\mathbb R`$.
Dann erfüllt $`S=\bigoplus_cW_c\mathrm{diag}(i\varphi_{c,j})_jW_c^*`$ die Gleichungen $`S^*=-S`$ und $`e^S=U_0^*U`$ (die
Reihe $`\sum_nS^n/n!`$ wird gliedweise in der Eigenbasis berechnet), und $`G=S/2`$ liefert $`G-G^*=S`$. $`\square`$

<a id="r-11-2"></a>

**Proposition 11.2 (die Parametrisierung ist surjektiv).** Für feste $`\delta`$, $`s_\pm`$ und $`U_0`$ ist das Bild von
$`(\theta,G)\mapsto X(\theta,G)`$ gleich $`\{X\in\mathcal P_\delta:1\in\mathrm{spec}\,X_{c_+},\ \ell\in\mathrm{spec}\,X_{c_-}\}`$,
und $`\mathcal P_\delta`$ ist die Vereinigung dieser Bilder über alle Paare verschiedener Plätze.

*Beweis.* Das Bild liegt in dieser Menge: $`X`$ hat Eigenwerte $`\lambda_s\in[\ell,1]`$, mit $`\lambda_{s_\pm}`$ in den
Blöcken $`c_\pm`$. Umgekehrt diagonalisiere man $`X_c=W_c\mathrm{diag}(\mu_c)W_c^*`$ und permutiere jedes $`\mu_c`$
(man ersetze $`W_c`$ durch $`W_c\Pi_c`$ mit einer Permutationsmatrix $`\Pi_c`$) so, dass $`\mu_{s_+}=1`$ und $`\mu_{s_-}=\ell`$.
Das ist möglich: Ist $`c_+=c_-`$ und $`\delta>0`$, so sind $`1\ne\ell`$ zwei Eigenwerte desselben Blocks und besetzen
verschiedene Positionen; ist $`\delta=0`$, so sind alle Eigenwerte gleich $`1`$. Ist $`\ell<1`$, so setze man
$`\theta_s=\arcsin\big(2(\mu_s-\ell)/(1-\ell)-1\big)`$ für die übrigen Plätze, sodass $`\lambda=\mu`$ (für $`\ell=1`$ ist
$`\theta`$ beliebig). Lemma 11.1 liefert $`G`$ mit $`U_0e^{G-G^*}=\bigoplus_cW_c`$. *Vereinigung:* Zu $`X\in\mathcal P_\delta`$
wähle man Blöcke $`c_\pm`$, die die Eigenwerte $`1`$ bzw. $`\ell`$ enthalten, und darin verschiedene Plätze; ist $`c_+=c_-`$
und $`\delta>0`$, so hat der Block zwei verschiedene Eigenwerte, also $`k_{c_+}\ge2`$; ist $`\delta=0`$, so ist $`X=1`$, und
zwei beliebige Plätze genügen. $`\square`$

*Rechnung:* Proposition 11.2 konstruktiv: Die Parameter aus dem Beweis reproduzieren einen Zieloperator — [Tutorial, Abschnitt 11.1](../../tutorial/de/11_optimisation.md#s-11-1).

Da jedes Paar von Plätzen mit positiver Wahrscheinlichkeit gezogen wird, liegt jeder Punkt der Schale mit positiver
Wahrscheinlichkeit im Bild eines Starts, und (in exakter Arithmetik) gilt die Nebenbedingung in jeder Iteration.
Für selbstadjungiertes $`X`$ (Exp. 3) verwendet man $`X=U\mathrm{diag}(\sigma\lambda)U^*`$ mit pro Start gezogenen Vorzeichen
$`\sigma\in\{\pm1\}^D`$ (pro Block ist die Anzahl der Minuszeichen gleichverteilt auf $`\{0,\dots,k_c\}`$ und ihre
Positionen sind gleichverteilt zufällig, sodass jedes Vorzeichenmuster positive Wahrscheinlichkeit hat), die nur durch
die Vorzeichenwechsel aus Abschnitt 11.5 verändert werden. Dann ist $`|X|=U\mathrm{diag}(\lambda)U^*`$
(Satz 1.5(e), da $`\lambda>0`$), und die Vereinigung über Vorzeichen und Plätze ist $`\{X=X^*:|X|\in\mathcal P_\delta\}`$:
Für ein solches $`X`$ diagonalisiere man $`X_c=W_c\mathrm{diag}(\nu_c)W_c^*`$, sodass
$`|X|_c=W_c\mathrm{diag}(|\nu_c|)W_c^*`$ mit $`|\nu_{c,i}|\ge\ell>0`$, und wende den Beweis von Proposition 11.2 auf
$`\mu=|\nu|`$ mit $`\sigma_s=\nu_s/|\nu_s|`$ an. Es kommt nur auf $`S=G-G^*`$ an.

*In der Bibliothek:* `optimize.SelfAdjointParam` — [Tutorial, Abschnitt 11.1](../../tutorial/de/11_optimisation.md#s-11-1).

Die Karte $`S\mapsto U_0e^S`$ verhält sich nur nahe $`S=0`$ gutartig. Seien $`\mathfrak u=\{S\in\mathcal M:S^*=-S\}`$ und
$`\mathcal M_h=\{H\in\mathcal M:H^*=H\}`$, reelle Unterräume mit $`\mathcal M=\mathfrak u\oplus\mathcal M_h`$
($`X=\frac{X-X^*}2+\frac{X+X^*}2`$). *Nahe $`0`$ ist $`S\mapsto U_0e^S`$ ein Homöomorphismus einer offenen Kugel $`N\ni0`$ in
$`\mathfrak u`$ auf eine Umgebung von $`U_0`$ in der unitären Gruppe.* Denn $`\Psi(S,H)=U_0e^Se^H`$ ist eine $`C^1`$-Abbildung
$`\mathfrak u\times\mathcal M_h\to\mathcal M`$ (die Reihe $`\sum_{m,n\ge0}U_0S^mH^n/(m!\,n!)`$ und die Reihe der partiellen
Ableitungen ihrer Glieder nach den reellen Koordinaten konvergieren gleichmäßig auf beschränkten Mengen), und
$`e^{tS}e^{tH}=1+t(S+H)+O(t^2)`$ zeigt, dass ihr Differential in $`0`$ die Abbildung $`(S,H)\mapsto U_0(S+H)`$ ist, eine
Bijektion auf $`\mathcal M`$. Nach dem Satz über die Umkehrfunktion gibt es offene Kugeln $`N\ni0`$ in $`\mathfrak u`$ und
$`N_h\ni0`$ in $`\mathcal M_h`$, sodass $`\Psi`$ ein Homöomorphismus von $`N\times N_h`$ auf eine offene Menge
$`O\subseteq\mathcal M`$ ist. Ist $`\Psi(S,H)`$ unitär, so auch $`e^H=e^{-S}U_0^*\Psi(S,H)`$ (in einer Eigenbasis der
normalen Matrix $`S`$, deren Eigenwerte imaginär sind, sind $`e^{\pm S}`$ diagonal mit Einträgen vom Betrag $`1`$ und
zueinander invers); aber $`e^H`$ hat die positiven Eigenwerte $`e^{h_i}`$, wobei $`h_i`$ die Eigenwerte von $`H`$ sind, also
ist $`e^{h_i}=1`$ (Korollar 1.4(2)) und $`H=0`$. Folglich ist $`O\cap U(\mathcal M)=U_0e^N`$, wobei $`U(\mathcal M)`$ die
unitäre Gruppe von $`\mathcal M`$ bezeichnet: eine Umgebung von $`U_0`$ in $`U(\mathcal M)`$, und $`S\mapsto U_0e^S`$ ist eine
stetige Bijektion $`N\to U_0e^N`$, deren Umkehrung, die erste Komponente von $`\Psi^{-1}`$, stetig ist. Fern von $`0`$ ist
die Karte nicht injektiv: $`e^{2\pi iP}=1=e^0`$ für jede Projektion $`P\in\mathcal M`$. Deshalb **zentriert** eine Suche die Karte nach jeder Runde **neu** ($`U_0\leftarrow U_0e^{G-G^*}`$, $`G\leftarrow0`$). Zwei weitere Karten haben in $`0`$ dasselbe
Differential $`S\mapsto U_0S`$. Die Cayley-Karte $`S\mapsto U_0C(-S/2)`$ (Aufgabe 1):
$`(1+tS/2)(1-tS/2)^{-1}=(1+tS/2)(1+tS/2+O(t^2))=1+tS+O(t^2)`$. Die polare Retraktion, $`S\mapsto`$ der unitäre Faktor von
$`A=U_0(1+S)`$ (Satz 1.8): $`A^*A=(1-S)(1+S)=1+S^*S\ge1`$ ist invertierbar, also ist $`|A|`$ invertierbar und der unitäre
Faktor ist $`A|A|^{-1}`$; für $`tS`$ statt $`S`$ ist $`|A|^{-1}=(1+t^2S^*S)^{-1/2}=1+O(t^2)`$ (in einer Eigenbasis von
$`S^*S`$), also ist der Faktor $`U_0(1+tS)+O(t^2)`$.

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Der Michelson-Kontrast ist die
> Sichtbarkeit $`(I_{\max}-I_{\min})/(I_{\max}+I_{\min})`$ von Interferenzstreifen. Für $`0\le X\le1`$ mit
> $`\lambda_{\max}=1`$ – ein Filter, ein Abschwächer, eine Detektoreffizienz – legt er den Dynamikbereich
> $`\lambda_{\min}/\lambda_{\max}`$ fest, den Kehrwert einer Konditionszahl.

<a id="s-11-2"></a>

## 11.2 Die von Neumannsche Spurungleichung

Eine Matrix $`(c_{ij})`$ mit $`c_{ij}\ge0`$, deren Zeilen- und Spaltensummen alle gleich $`1`$ sind, heißt **doppelt
stochastisch**; für eine unitäre Matrix $`G`$ sind $`(|G_{ij}|^2)`$ und $`(|G_{ji}|^2)`$ doppelt stochastisch, da die Zeilen
und Spalten von $`G`$ Einheitsvektoren sind.

<a id="r-11-3"></a>

**Lemma 11.3 (Umordnung).** Seien $`p_1\ge\dots\ge p_k`$, $`q_1\ge\dots\ge q_k`$ reell und $`(c_{ij})`$ doppelt
stochastisch. Dann gilt $`\sum_ip_iq_{k+1-i}\le\sum_{i,j}p_iq_jc_{ij}\le\sum_ip_iq_i`$.

*Beweis.* *Obere Schranke.* Für $`m\le k`$ sei $`r_j=\sum_{i\le m}c_{ij}\in[0,1]`$, also $`\sum_jr_j=m`$. Da $`q_j\ge q_m`$
für $`j\le m`$ und $`q_j\le q_m`$ für $`j>m`$, gilt
```math
\sum_{j\le m}q_j-\sum_jq_jr_j=\sum_{j\le m}q_j(1-r_j)-\sum_{j>m}q_jr_j\ge q_m\Big(m-\sum_jr_j\Big)=0 .
```
Mit $`x_i=\sum_jq_jc_{ij}`$ besagt dies $`\sum_{i\le m}x_i\le\sum_{j\le m}q_j`$, mit Gleichheit für $`m=k`$. Mit abelscher
partieller Summation und $`p_m-p_{m+1}\ge0`$ folgt
```math
\sum_ip_ix_i=p_k\sum_{i\le k}x_i+\sum_{m<k}(p_m-p_{m+1})\sum_{i\le m}x_i
\le p_k\sum_{j\le k}q_j+\sum_{m<k}(p_m-p_{m+1})\sum_{j\le m}q_j=\sum_ip_iq_i .
```
*Untere Schranke.* Man wende die obere Schranke auf $`p`$, die fallende Folge $`q'_j=-q_{k+1-j}`$ und die doppelt
stochastische Matrix $`c'_{ij}=c_{i,k+1-j}`$ an:
$`-\sum_{ij}p_iq_jc_{ij}=\sum_{ij}p_iq'_jc'_{ij}\le\sum_ip_iq'_i=-\sum_ip_iq_{k+1-i}`$. $`\square`$

Jedes $`A\in M_k(\mathbb C)`$ besitzt eine **Singulärwertzerlegung** $`A=W\Sigma V^*`$ mit unitären $`W,V`$ und
$`\Sigma=\mathrm{diag}(s_i(A))`$: Man schreibe $`A=W'|A|`$ (Satz 1.8), $`|A|=V\Sigma V^*`$ (Satz 1.3) und setze
$`W=W'V`$.

<a id="r-11-4"></a>

**Satz 11.4 (von Neumannsche Spurungleichung).** Für $`A,B\in M_k(\mathbb C)`$ gilt
$`|\mathrm{Tr}(AB)|\le\sum_is_i(A)\,s_i(B)`$.

*Beweis.* Seien $`A=W_A\Sigma_AV_A^*`$, $`B=W_B\Sigma_BV_B^*`$, $`a_i=s_i(A)`$, $`b_j=s_j(B)`$ sowie $`G=V_A^*W_B`$,
$`H=V_B^*W_A`$. Wegen der Zyklizität ist
$`\mathrm{Tr}(AB)=\mathrm{Tr}(\Sigma_AG\Sigma_BH)=\sum_{i,j}a_iG_{ij}b_jH_{ji}`$, also wegen
$`|g||h|\le\frac12(|g|^2+|h|^2)`$ und $`a_ib_j\ge0`$
```math
|\mathrm{Tr}(AB)|\le\sum_{i,j}a_ib_j\,d_{ij},\qquad d_{ij}=\tfrac12\big(|G_{ij}|^2+|H_{ji}|^2\big).
```
$`(d_{ij})`$ ist doppelt stochastisch, und Lemma 11.3 liefert $`\sum_{ij}a_ib_jd_{ij}\le\sum_ia_ib_i`$. $`\square`$

<a id="r-11-5"></a>

**Korollar 11.5.** (a) Für selbstadjungierte $`P,Q\in M_k(\mathbb C)`$ mit Eigenwerten $`p,q`$ gilt
$`\sum_ip^\downarrow_iq^\uparrow_i\le\mathrm{Tr}(PQ)\le\sum_ip^\downarrow_iq^\downarrow_i`$.
(b) Für $`P,Q\ge0`$ gilt $`\mathrm{Tr}(PQ)\le\|PQ\|_1\le\sum_ip^\downarrow_iq^\downarrow_i`$.

*Beweis.* (a) Wie im Beweis von Satz 3.5 ist $`\mathrm{Tr}(PQ)=\sum_{ij}p_iq_jc_{ij}`$ mit der doppelt
stochastischen Matrix $`c_{ij}=|\langle u_i,v_j\rangle|^2`$ für fallend geordnete Eigenbasen; man wende Lemma 11.3 an
(dessen untere Schranke ist $`\sum_ip^\downarrow_iq^\downarrow_{k+1-i}=\sum_ip^\downarrow_iq^\uparrow_i`$). (b) Die erste
Ungleichung ist Satz 1.9. Nach Satz 1.9 und Satz 11.4 ist
$`\|PQ\|_1=\max_U|\mathrm{Tr}(P\,QU)|\le\sum_is_i(P)s_i(QU)`$, und $`s_i(QU)=q^\downarrow_i`$, weil
$`|QU|=(U^*Q^2U)^{1/2}=U^*QU`$ (Satz 1.6). $`\square`$

*Rechnung:* Satz 11.4 und Korollar 11.5 an Zufallsmatrizen — [Tutorial, Abschnitt 11.2](../../tutorial/de/11_optimisation.md#s-11-2).

Für $`X,Y\ge0`$ in $`M_k(\mathbb C)`$ mit Spektren $`x,y`$ beträgt die Verletzung $`\|YX\|_1-\mathrm{Tr}(XY)`$ in
Proposition 3.9 also höchstens $`\sum_ix^\downarrow_i(y^\downarrow_i-y^\uparrow_i)`$ (Korollar 11.5(b) für $`(Y,X)`$ und die
untere Schranke aus (a)). Für Exp. 2 ist diese Zahl die exakte Antwort:

<a id="r-11-6"></a>

**Proposition 11.6 (Exp. 2 bei festen Spektren).** In Exp. 2 ist $`|X|=X_0`$, $`|Y|=Y_0`$ und
$`z=\|X_0VY_0\|_1-\mathrm{Tr}(X_0Y_0)`$. Sind die Spektren $`x_c`$ von $`(X_0)_c`$ und $`y_c`$ von $`(Y_0)_c`$ fest und
variieren die Eigenbasen und $`V`$, so gilt
```math
\sup z=-\inf z=\sum_c\sum_ix^\downarrow_{c,i}\big(y^\downarrow_{c,i}-y^\uparrow_{c,i}\big),
```
und beide Werte werden angenommen.

*Beweis.* $`X^*X=X_0U^*UX_0=X_0^2`$, also $`|X|=X_0`$ (Satz 1.6); ebenso $`|Y|=Y_0`$, und $`|UA|=|A|`$ liefert
$`\mathrm{Tr}|XY|=\|X_0VY_0\|_1`$. Alle Variablen wirken blockweise, und $`z`$ ist eine Summe über die Blöcke; sei
also $`C=1`$. Man setze $`Y'=VY_0V^*\ge0`$ (Spektrum $`y`$); wegen $`|AV|=V^*|A|V`$ ist
$`\|X_0VY_0\|_1=\|X_0Y'V\|_1=\|X_0Y'\|_1`$. Korollar 11.5 für die Paare $`(X_0,Y')`$ und $`(X_0,Y_0)`$ liefert zusammen
mit $`\|X_0Y'\|_1\ge\mathrm{Tr}(X_0Y')`$
```math
\sum_ix^\downarrow_iy^\uparrow_i-\sum_ix^\downarrow_iy^\downarrow_i\le\|X_0Y'\|_1-\mathrm{Tr}(X_0Y_0)\le\sum_ix^\downarrow_iy^\downarrow_i-\sum_ix^\downarrow_iy^\uparrow_i .
```
Sei $`J`$ die Umkehrpermutation, $`Je_i=e_{k+1-i}`$; dann ist $`J\mathrm{diag}(w)=\mathrm{diag}(w_k,\dots,w_1)J`$, und
die Spurnorm von $`\mathrm{diag}(w)J`$ ist $`\sum_i|w_i|`$. Mit $`X_0=\mathrm{diag}(x^\downarrow)`$ und $`V=J`$ gilt:
Für $`Y_0=\mathrm{diag}(y^\uparrow)`$ ist $`X_0VY_0=\mathrm{diag}(x^\downarrow_iy^\downarrow_i)J`$ und
$`\mathrm{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\uparrow_i`$ (die obere Schranke); für $`Y_0=\mathrm{diag}(y^\downarrow)`$
ist $`X_0VY_0=\mathrm{diag}(x^\downarrow_iy^\uparrow_i)J`$ und $`\mathrm{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\downarrow_i`$
(die untere Schranke). $`\square`$

<a id="r-11-7"></a>

**Korollar 11.7 (Einhüllende von Exp. 2).** Seien $`\ell_X,\ell_Y\in(0,1]`$ und $`k_c\ge2`$ für mindestens ein $`c`$. Über alle
$`X_0,Y_0\ge0`$ mit $`\lambda_{\max}=1`$, $`\lambda_{\min}(X_0)=\ell_X`$, $`\lambda_{\min}(Y_0)=\ell_Y`$ und alle unitären
Matrizen gilt
```math
\sup z=-\inf z=\sum_c\Big\lfloor\frac{k_c}2\Big\rfloor(1-\ell_X)(1-\ell_Y).
```

*Beweis.* In einem Block schreibe man $`x_i=x^\downarrow_{c,i}`$, $`y_i=y^\downarrow_{c,i}`$, $`k=k_c`$, sodass
$`y^\uparrow_{c,i}=y_{k+1-i}`$. Paart man $`i`$ mit $`k+1-i`$, so ist
$`x_i(y_i-y_{k+1-i})+x_{k+1-i}(y_{k+1-i}-y_i)=(x_i-x_{k+1-i})(y_i-y_{k+1-i})`$, und für ungerades $`k`$ trägt der mittlere
Index $`x_i(y_i-y_i)=0`$ bei. Folglich gilt
```math
\sum_ix_i(y_i-y_{k+1-i})=\sum_{i\le\lfloor k/2\rfloor}(x_i-x_{k+1-i})(y_i-y_{k+1-i})\le\Big\lfloor\frac k2\Big\rfloor(1-\ell_X)(1-\ell_Y),
```
da $`x_i\in[\ell_X,1]`$, $`y_i\in[\ell_Y,1]`$ fallend sind, sodass beide Faktoren in $`[0,1-\ell_X]`$ bzw. $`[0,1-\ell_Y]`$
liegen. Gleichheit gilt, wenn in jedem Block die $`\lfloor k_c/2\rfloor`$ größten Eigenwerte gleich $`1`$ und die
$`\lfloor k_c/2\rfloor`$ kleinsten gleich $`\ell_X`$ (bzw. $`\ell_Y`$) sind, während der mittlere (bei ungeradem $`k_c`$)
beliebig in $`[\ell_X,1]`$ (bzw. $`[\ell_Y,1]`$) ist; diese Spektren haben die globalen Extrema $`1`$ und $`\ell_X`$ ($`\ell_Y`$), weil $`k_c\ge2`$ für mindestens ein $`c`$ gilt.
Proposition 11.6 schließt den Beweis ab. $`\square`$

Hier ist $`1-\ell_X=2\delta/(1+\delta)`$ für $`\Delta(|X|)=\delta`$; für ungerades $`k_c`$ spielt der mittlere Eigenwert keine
Rolle.

*Rechnung:* die extremalen Konfigurationen aus Proposition 11.6 und die Einhüllende aus Korollar 11.7 für $`k=3`$, $`C=2`$ — [Tutorial, Abschnitt 11.2](../../tutorial/de/11_optimisation.md#s-11-2); die von einer Multistart-Suche gefundene Einhüllende — [Tutorial, Abschnitt 11.5](../../tutorial/de/11_optimisation.md#s-11-5).

<a id="s-11-3"></a>

## 11.3 Eine Spurungleichung mit Gleichheitsbedingung

Für $`X=X^*`$ sei $`X_\pm=f_\pm(X)`$ mit $`f_+(t)=\max(t,0)`$, $`f_-(t)=\max(-t,0)`$. Dann gilt $`X_\pm\ge0`$, $`X=X_+-X_-`$,
$`|X|=X_++X_-`$ und $`X_+X_-=0`$ (Satz 1.5(c), und 1.5(e) für $`|X|=(X^2)^{1/2}=h(X)`$ mit
$`h(t)=\sqrt{t^2}=|t|=f_+(t)+f_-(t)`$; Korollar 1.4(3) für die Positivität).

<a id="r-11-8"></a>

**Lemma 11.8.** Für $`A,B\ge0`$ in $`M_n(\mathbb C)`$ gilt $`\|A-B\|_1\le\mathrm{Tr}\,A+\mathrm{Tr}\,B`$, mit Gleichheit
genau dann, wenn $`AB=0`$.

*Beweis.* Sei $`H=A-B`$ und $`Q=g(H)`$ mit $`g=1`$ auf $`(0,\infty)`$ und $`g=0`$ sonst; $`Q`$ ist eine Projektion, da
$`g=\bar g=g^2`$ (Satz 1.5(c)). Wegen $`|t|=t(2g(t)-1)`$ ist $`|H|=H(2Q-1)`$ (Satz 1.5(c),(e)), und mit
$`\mathrm{Tr}(AR)=\mathrm{Tr}(RAR)`$ für die Projektionen $`R=Q,1-Q`$ (Zyklizität) folgt
```math
\|H\|_1=\mathrm{Tr}\big(A(2Q-1)\big)-\mathrm{Tr}\big(B(2Q-1)\big)
=\mathrm{Tr}\,A+\mathrm{Tr}\,B-2\mathrm{Tr}\big((1-Q)A(1-Q)\big)-2\mathrm{Tr}(QBQ).
```
Die abgezogenen Spuren sind $`\ge0`$, und $`\mathrm{Tr}(RAR)=\|A^{1/2}R\|_2^2`$ verschwindet genau dann, wenn $`AR=0`$;
Gleichheit gilt also genau dann, wenn $`A(1-Q)=0`$ und $`BQ=0`$, und dann ist $`AB=AQB=A(BQ)^*=0`$. Ist umgekehrt
$`AB=0`$, so ist $`BA=0`$, $`(A-B)^2=A^2+B^2=(A+B)^2`$, und Satz 1.6 liefert $`|A-B|=A+B`$. $`\square`$

<a id="r-11-9"></a>

**Proposition 11.9 (Exp. 3).** Für $`X=X^*`$ und $`Y\ge0`$ in $`\mathcal M`$ gilt
$`\mathrm{Tr}|YXY|\le\mathrm{Tr}(Y|X|Y)`$, mit Gleichheit genau dann, wenn $`YX_+Y^2X_-Y=0`$, insbesondere
wenn $`XY=YX`$.

*Beweis.* In $`M_D(\mathbb C)`$ sind $`A=YX_+Y`$ und $`B=YX_-Y`$ positiv (Korollar 1.7 mit $`T=X_\pm^{1/2}Y`$), $`YXY=A-B`$
und $`\mathrm{Tr}(Y|X|Y)=\mathrm{Tr}\,A+\mathrm{Tr}\,B`$; man wende Lemma 11.8 an. Ist $`XY=YX`$, so kommutiert
$`Y`$ mit $`X_\pm`$ (Satz 1.5(d)), und $`AB=Y^4X_+X_-=0`$. $`\square`$

Also ist $`z\ge0`$ in Exp. 3, und $`\underline z=0`$ (man nehme $`Y=1`$). Die Gleichheitsbedingung ist schwächer als
Kommutativität: Für $`X=\mathrm{diag}(1,-1,0)`$ und $`Y`$ mit den Zeilen $`(1,0,a)`$, $`(0,1,0)`$, $`(a,0,1)`$, $`-1<a<1`$, ist
$`AB=(Ye_1)\langle Ye_1,Ye_2\rangle(Ye_2)^*=0`$, während $`XY\ne YX`$ für $`a\ne0`$.

*Rechnung:* Proposition 11.9 an zufälligen Paaren und das Beispiel für Gleichheit ohne Kommutativität — [Tutorial, Abschnitt 11.3](../../tutorial/de/11_optimisation.md#s-11-3).

<a id="s-11-4"></a>

## 11.4 Stichproben oder Optimierung

Wie groß ist $`|\mathrm{Tr}(AU)|`$ typischerweise, verglichen mit seinem Maximum $`\|A\|_1`$? Wir ziehen $`U`$
als den Faktor $`Q`$ der QR-Zerlegung $`Z=QR`$, $`R_{ii}>0`$, einer **Ginibre-Matrix** $`Z`$ (unabhängige Einträge mit der
Dichte $`\pi^{-1}e^{-|z|^2}`$; Mezzadris Rezept), und zwar in jedem Block unabhängig.

<a id="r-11-10"></a>

**Lemma 11.10.** $`Q`$ ist Haar-verteilt auf $`U(k)`$; insbesondere hat $`VQ`$ für jede feste unitäre Matrix
$`V\in M_k(\mathbb C)`$ dieselbe Verteilung wie $`Q`$. Die unitäre Matrix $`\bigoplus_cQ_c`$ mit unabhängigen Blöcken ist
Haar-verteilt auf der unitären Gruppe von $`\mathcal M`$.

*Beweis.* Das ist Satz 5.3 (die Linksinvarianz ist der Äquivarianzschritt seines Beweises; die Verteilung wird
durch das Eindeutigkeitslemma 5.1, dessen maßtheoretischer Schritt Korollar A.2 ist, als Haar-Maß identifiziert)
zusammen mit Korollar 5.4. $`\square`$

*In der Bibliothek:* Diese Stichprobenziehung ist `TypeIAlgebra.random_unitary_operator` — [Tutorial, Abschnitt 11.4](../../tutorial/de/11_optimisation.md#s-11-4).

<a id="r-11-11"></a>

**Proposition 11.11.** Für Haar-verteiltes $`U\in U(k)`$ und $`A\in M_k(\mathbb C)`$ gilt $`\mathbb E\mathrm{Tr}(AU)=0`$
und $`\mathbb E|\mathrm{Tr}(AU)|^2=\|A\|_2^2/k`$.

*Beweis.* $`\mathbb EU_{ij}=0`$ nach Korollar 5.5. Proposition 5.7 (mit $`C=1`$) für $`X=|e_j\rangle\langle e_m|`$ liefert
die Weingarten-Formel erster Ordnung,
$`\mathbb E[U_{ij}\bar U_{lm}]=\mathbb E(UXU^*)_{il}=\delta_{il}\delta_{jm}/k`$; für $`k\ge2`$ folgt sie auch aus der Formel
zweiter Ordnung, Satz 5.10, indem man $`i_2=i'_2`$, $`j_2=j'_2`$ setzt und über beide summiert, da
$`\sum_{i_2,j_2}|U_{i_2j_2}|^2=k`$. Mit $`\mathrm{Tr}(AU)=\sum_{ij}A_{ji}U_{ij}`$ ist
$`\mathbb E|\mathrm{Tr}(AU)|^2=\sum_{ijlm}A_{ji}\bar A_{ml}\mathbb E[U_{ij}\bar U_{lm}]=\frac1k\sum_{ij}|A_{ji}|^2`$. $`\square`$

<a id="r-11-12"></a>

**Korollar 11.12.** Seien $`X,Y\ge0`$ in $`M_k(\mathbb C)`$ mit Spektren in $`[a,b]`$, $`a>0`$, und $`U`$ Haar-verteilt. Dann
gilt
```math
\frac{\big(\mathbb E|\mathrm{Tr}(XUY)|^2\big)^{1/2}}{\max_V|\mathrm{Tr}(XVY)|}=\frac{\|YX\|_2}{\sqrt k\,\|YX\|_1}\le\min\Big(\frac1{\sqrt k},\frac{b^2}{a^2k}\Big),
\qquad
\mathbb P\big(|\mathrm{Tr}(XUY)|\ge\mathrm{Tr}(XY)\big)\le\frac{b^4}{a^4k^2}.
```

*Beweis.* Die Gleichung ist Proposition 11.11 für $`A=YX`$ (denn $`\mathrm{Tr}(XUY)=\mathrm{Tr}(AU)`$)
zusammen mit Proposition 3.9. Da $`\|A\|_2^2=\sum_is_i^2\le(\sum_is_i)^2=\|A\|_1^2`$ (Abschnitt 1.4), ist der Quotient
$`\le1/\sqrt k`$. Wegen $`\|Av\|=\||A|v\|`$ (Abschnitt 1.4) liefert Diagonalisieren von $`|A|`$
$`s_k(A)=\min_{\|v\|=1}\|Av\|`$ und $`s_1(A)=\max_{\|v\|=1}\|Av\|`$; und $`\|Xv\|^2=\langle v,X^2v\rangle\in[a^2,b^2]`$ für
Einheitsvektoren $`v`$ (Korollar 1.4(3) für $`X^2-a^2`$ und $`b^2-X^2`$), ebenso für $`Y`$. Also ist
$`a^2\le s_k(A)\le s_1(A)\le b^2`$, $`\|A\|_2^2\le kb^4`$ und $`\|A\|_1\ge ka^2`$. Schließlich ist
$`\mathrm{Tr}(XY)=\mathrm{Tr}(X^{1/2}YX^{1/2})\ge a\mathrm{Tr}\,X\ge ka^2`$, und die Markov-Ungleichung für
$`|\mathrm{Tr}(XUY)|^2`$ beschränkt die Wahrscheinlichkeit durch $`(kb^4/k)/(ka^2)^2`$. $`\square`$

$`N`$ Haar-Stichproben finden eine Verletzung von $`|\mathrm{Tr}(XUY)|\le\mathrm{Tr}(XY)`$ also höchstens mit
Wahrscheinlichkeit $`Nb^4/(a^4k^2)`$ (eine Vereinigung von $`N`$ Ereignissen), obwohl es für jedes nichtkommutierende
Paar eine gibt (Proposition 3.9). Die Schranke ist erst für $`k>\sqrt N\,b^2/a^2`$ aussagekräftig (für $`b/a=6`$ und
$`N=2000`$: $`k>1609`$).

*Rechnung:* Proposition 11.11 per Monte-Carlo und der Quotient aus Korollar 11.12 im Vergleich mit der besten von $`2000`$ Stichproben für $`k\le32`$ — [Tutorial, Abschnitt 11.4](../../tutorial/de/11_optimisation.md#s-11-4).

<a id="r-11-13"></a>

**Lemma 11.13.** Ist $`k\ge2`$, so gilt $`\min_U|\mathrm{Tr}(AU)|=0`$ für jedes $`A\in M_k(\mathbb C)`$.

*Beweis.* Sei $`A=W\Sigma V^*`$ und $`P`$ die zyklische Verschiebung, $`Pe_i=e_{(i\bmod k)+1}`$, deren Diagonale für $`k\ge2`$
verschwindet. Für $`U=VPW^*`$ ist $`\mathrm{Tr}(AU)=\mathrm{Tr}(\Sigma P)=0`$. $`\square`$

Also ist $`\inf_Uz=-\mathrm{Tr}(XY)`$ und $`\sup_Uz=\|YX\|_1-\mathrm{Tr}(XY)`$ (Proposition 3.9). Nach Proposition 11.11 ist das quadratische
Mittel von $`|\mathrm{Tr}(XUY)|`$ gleich $`\|YX\|_2/\sqrt k`$; mit den Schranken aus dem Beweis von Korollar 11.12
(Spektren in $`[a,b]`$) ist $`\mathbb E|\mathrm{Tr}(XUY)|^2\le b^4`$ und $`\mathrm{Tr}(XY)\ge ka^2`$, sodass die
Markov-Ungleichung $`\mathbb P\big(z\ge-(1-\varepsilon)\mathrm{Tr}(XY)\big)\le b^4/(\varepsilon^2a^4k^2)`$ für
$`\varepsilon>0`$ liefert: Für großes $`k`$ liegt $`z`$ mit hoher Wahrscheinlichkeit nahe bei seinem Infimum.

*Rechnung:* beide Extrema über $`U`$ durch Gradientenanstieg und -abstieg, im Vergleich mit $`20\,000`$ Haar-Stichproben — [Tutorial, Abschnitt 11.4](../../tutorial/de/11_optimisation.md#s-11-4).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Korollar 11.12 ist eine schwache,
> beweisbare Instanz dessen, was die statistische Mechanik **Typizität** nennt: In hoher Dimension erwartet man, dass
> fast jeder Zustand und fast jede Zeitentwicklung nahezu denselben Wert einer makroskopischen Größe liefert. Über
> eine scharfe Ungleichung entscheiden fein abgestimmte Konfigurationen – ausgerichtete Eigenbasen,
> Zwei-Niveau-Spektren –, die eine zufällige unitäre Matrix (ein Zufallsschaltkreis, eine zufällige Streumatrix)
> selten hervorbringt.

<a id="s-11-5"></a>

## 11.5 Multistart-Suche nach Einhüllenden

Eine **Multistart-Suche** führt ein Gradientenverfahren in Runden von vielen unabhängigen Starts aus, jeder mit
eigenen Plätzen, Vorzeichen, $`U_0`$ und eigenem Kontrast; sie zentriert die unitären Karten nach jeder Runde neu und gibt
pro Start den besten Wert der Zielfunktion an den Iterierten zurück. Die Vorzeichen der selbstadjungierten
Parametrisierung sind diskret; zwischen den Runden werden sie durch einzelne Vorzeichenwechsel an zufälligen Plätzen
verändert, und jeder Wechsel wird behalten, wenn er die Zielfunktion verbessert.

*In der Bibliothek:* `optimize.extremize` — [Tutorial, Abschnitt 11.5](../../tutorial/de/11_optimisation.md#s-11-5).

In exakter Arithmetik ist jede Iterierte zulässig (Proposition 11.2 und die Bemerkungen danach), also ist jeder
zurückgegebene Wert der Wert von $`z`$ an einem zulässigen Punkt; in Gleitkommaarithmetik ist er dieser Wert bis auf
Rundungsfehler. Daher liefert eine Maximierung eine **untere Schranke** für
das Supremum, eine Minimierung eine obere Schranke für das Infimum, und nicht mehr – eine lokale Suche kann in einem
lokalen Extremum stehen bleiben. Die andere Richtung erfordert einen Beweis, wie Korollar 11.7 und
Proposition 11.9, oder eine konvexe Relaxierung (Abschnitt 11.6). Eine berechnete Einhüllende wird geprüft (i) an
geschlossenen Formeln; (ii) an von Hand gebauten Konfigurationen, etwa direkten Summen kleinerer Lösungen, unter die
sie nicht fallen darf; (iii) an strukturierten Stichproben (Zwei-Niveau-Spektren, Permutationsbasen), deren Werte
weitere angenommene Werte sind.

<a id="r-11-14"></a>

**Proposition 11.14 (untere Einhüllende von Exp. 1).** Seien alle $`k_c\ge2`$. Über $`X\in\mathcal P_\delta`$, $`Y\ge0`$
mit $`\lambda_{\max}(Y)=1`$ und unitäre $`U\in\mathcal M`$ gilt $`\ \min z=-\big(D-2\delta/(1+\delta)\big)`$.

*Beweis.* Aus $`Y\le1`$ folgt $`X^{1/2}YX^{1/2}\le X`$, also $`\mathrm{Tr}(XY)\le\mathrm{Tr}\,X\le(D-1)+\ell(\delta)`$,
da ein Eigenwert von $`X`$ gleich $`\ell(\delta)`$ ist und die übrigen $`\le1`$ sind; daher
$`z\ge-\mathrm{Tr}(XY)\ge-(D-2\delta/(1+\delta))`$. Gleichheit gilt für $`X`$ mit dem Eigenwert $`\ell(\delta)`$ einmal
und $`1`$ sonst, $`Y=1`$ und $`U=\bigoplus_cU_c`$ mit $`\mathrm{Tr}(X_cU_c)=0`$ (Lemma 11.13). $`\square`$

Prüfungen vom Typ (ii) können direkte Summen verwenden: $`z`$ aus Exp. 3 ist additiv über direkte Summen, und der Kontrast einer direkten Summe ist der größte
Kontrast ihrer Summanden (Aufgabe 3), also gilt für die Einhüllenden von Exp. 3 mit $`\Delta(Y)`$ frei in einem festen Intervall $`[0,b]`$
```math
\overline z_{M_2\oplus M_2}(\delta)\ \ge\ \overline z_{M_2}(\delta)+\sup_{\delta'\le\delta}\overline z_{M_2}(\delta');
```
insbesondere ist die Summe zweier auf $`M_2`$ bei den Kontrasten $`\delta`$ und $`\delta'\le\delta`$ angenommener Werte ein
auf $`M_2\oplus M_2`$ beim Kontrast $`\delta`$ angenommener Wert.

*Rechnung:* Proposition 11.14 auf $`M_2\oplus M_2`$, Korollar 11.7 auf $`M_3\oplus M_3`$ und Proposition 11.9 per Multistart-Suche, mit einer Prüfung vom Typ (ii) — [Tutorial, Abschnitt 11.5](../../tutorial/de/11_optimisation.md#s-11-5).

Für Blöcke der Größe $`2`$ sind die oberen Einhüllenden von Exp. 1 und Exp. 3 exakt bekannt. In $`M_2(\mathbb C)`$
gehen die Eigenbasen nur über eine Zahl ein: Für Orthonormalbasen aus Eigenvektoren $`e_1,e_2`$ und $`f_1,f_2`$ hat die
doppelt stochastische Matrix $`(|\langle e_i,f_j\rangle|^2)`$ die Zeilen $`(\gamma,1-\gamma)`$ und $`(1-\gamma,\gamma)`$, $`\gamma=|\langle e_1,f_1\rangle|^2\in[0,1]`$, und
jedes $`\gamma`$ tritt schon für reelle Basen auf ($`f_1=\cos t\,e_1+\sin t\,e_2`$, $`f_2=-\sin t\,e_1+\cos t\,e_2`$, $`\gamma=\cos^2t`$).

<a id="r-11-15"></a>

**Proposition 11.15 (Exp. 1 und Exp. 3 auf $`M_2`$).** Sei $`Y\ge0`$ in $`M_2(\mathbb C)`$ mit Eigenwerten
$`y_1\ge y_2\ge0`$, Eigenvektoren $`f_1,f_2`$ und $`y_1>0`$.

(a) Sei $`X\ge0`$ mit Eigenwerten $`x_1\ge x_2\ge0`$, $`x_1>0`$, und Eigenvektoren $`e_1,e_2`$. Für jedes unitäre $`U`$ gilt
```math
\lvert\mathrm{Tr}(XUY)\rvert-\mathrm{Tr}(XY)\ \le\ \frac{(x_1-x_2)^2(y_1-y_2)^2}{4(x_1+x_2)(y_1+y_2)} .
```
Ist $`x_1>x_2`$ und $`y_1>y_2`$, so gilt Gleichheit genau dann, wenn $`|\mathrm{Tr}(XUY)|=\|YX\|_1`$ und
$`\gamma=\gamma_1:=\frac14+\frac{x_1y_2+x_2y_1}{2(x_1+x_2)(y_1+y_2)}`$; andernfalls ist die Schranke $`0`$, und sie wird genau dann
angenommen, wenn $`|\mathrm{Tr}(XUY)|=\|YX\|_1`$.

(b) Sei $`X=X^*`$, und $`|X|`$ habe die Eigenwerte $`x_1\ge x_2\ge0`$, $`x_1>0`$. Dann gilt
```math
\mathrm{Tr}(Y|X|Y)-\mathrm{Tr}|YXY|\ \le\ \frac{2x_1x_2}{x_1+x_2}\,(y_1-y_2)^2 .
```
Hat $`X`$ Eigenwerte beider Vorzeichen und ist $`y_1>y_2`$, so seien $`e_1`$, $`e_2`$ Eigenvektoren von $`X`$ zu seinen
Eigenwerten vom Betrag $`x_1`$ und $`x_2`$; dann gilt Gleichheit genau dann, wenn $`\gamma=\gamma_3:=\frac{x_1y_2+x_2y_1}{(x_1+x_2)(y_1+y_2)}`$.

Für alle Eigenwerte wie oben werden beide Schranken angenommen, und zwar von reellen $`X`$, $`Y`$ und $`U`$ (in (b) mit
$`X=\mathrm{diag}(x_1,-x_2)`$); insbesondere sind sie in $`M_2(\mathbb R)`$ dieselben.

*Beweis.* Wir schreiben $`z`$ für die linke Seite und setzen $`S=x_1y_2+x_2y_1`$, $`p=(x_1-x_2)(y_1-y_2)\ge0`$ und
$`q=(x_1+x_2)(y_1+y_2)>0`$, sodass $`q-2S=p`$.

(a) Nach Satz 1.9 ist $`|\mathrm{Tr}(XUY)|=|\mathrm{Tr}(YX\,U)|\le\|YX\|_1`$, mit Gleichheit für $`U=W^*`$, $`YX=W|YX|`$.
Wie im Beweis von Korollar 11.5 ist $`\mathrm{Tr}(PQ)=\sum_{i,j}a_ib_j|\langle e_i,f_j\rangle|^2`$ für $`P=\sum_ia_i|e_i\rangle\langle e_i|`$,
$`Q=\sum_jb_j|f_j\rangle\langle f_j|`$; daher
$`\mathrm{Tr}(XY)=(x_1y_2+x_2y_1)+\gamma(x_1-x_2)(y_1-y_2)=S+\gamma p`$ und
$`\mathrm{Tr}(X^2Y^2)=(x_1^2y_2^2+x_2^2y_1^2)+\gamma(x_1^2-x_2^2)(y_1^2-y_2^2)`$. Die Quadrate $`s_1^2,s_2^2`$ der
Singulärwerte von $`YX`$ sind die Eigenwerte von $`(YX)^*YX=XY^2X`$, also $`s_1^2+s_2^2=\mathrm{Tr}(X^2Y^2)`$ und
$`s_1s_2=|\det(YX)|=x_1x_2y_1y_2`$. Folglich
```math
\|YX\|_1^2=(s_1+s_2)^2=\mathrm{Tr}(X^2Y^2)+2x_1x_2y_1y_2=S^2+\gamma pq .
```
Für $`r\ge0`$ ist $`\sqrt r\le r/q+q/4`$, mit Gleichheit genau für $`\sqrt r=q/2`$ (man multipliziere $`(\sqrt r-q/2)^2\ge0`$ aus). Mit $`r=S^2+\gamma pq`$ folgt
```math
z\ \le\ \|YX\|_1-\mathrm{Tr}(XY)\ \le\ \frac{S^2+\gamma pq}q+\frac q4-S-\gamma p=\frac{(q/2-S)^2}q=\frac{p^2}{4q}.
```
Gleichheit gilt genau dann, wenn $`|\mathrm{Tr}(XUY)|=\|YX\|_1`$ und $`S^2+\gamma pq=q^2/4`$, d. h. $`\gamma pq=(q/2-S)(q/2+S)=\frac p2(q/2+S)`$.
Für $`p>0`$ heißt das $`\gamma=(q/2+S)/(2q)=\gamma_1`$, und das liegt in $`[\frac14,\frac12]`$, da $`0\le S\le q/2`$; für $`p=0`$ gilt es
für jedes $`\gamma`$.

(b) Ersetzt man $`X`$ durch $`-X`$, so ändert sich keine der beiden Seiten. Ist $`X\ge0`$, so ist $`YXY\ge0`$, $`|YXY|=YXY`$ und $`z=0`$;
ebenso für $`X\le0`$. Andernfalls ist $`x_2>0`$, und nach eventuellem Ersetzen von $`X`$ durch $`-X`$ gilt $`Xe_1=x_1e_1`$,
$`Xe_2=-x_2e_2`$. Sei $`u=\langle e_1,Y^2e_1\rangle=\gamma y_1^2+(1-\gamma)y_2^2`$ und $`v=\langle e_2,Y^2e_2\rangle=s-u`$, $`s=y_1^2+y_2^2`$.
Dann ist $`\mathrm{Tr}(Y|X|Y)=\mathrm{Tr}(|X|Y^2)=x_1u+x_2v=:A`$. Die selbstadjungierte Matrix $`YXY`$ hat die Determinante
$`\det X(\det Y)^2=-x_1x_2y_1^2y_2^2\le0`$, also Eigenwerte $`\mu_+\ge0\ge\mu_-`$, und
```math
\mathrm{Tr}|YXY|=\mu_+-\mu_-=\sqrt{(\mu_++\mu_-)^2-4\mu_+\mu_-}=\sqrt{L^2+K^2},\qquad L=\mathrm{Tr}(XY^2)=x_1u-x_2v,\quad K=2\sqrt{x_1x_2}\,y_1y_2 .
```
Sei $`\alpha=(x_1-x_2)/(x_1+x_2)`$ und $`\beta=2\sqrt{x_1x_2}/(x_1+x_2)>0`$, sodass $`\alpha^2+\beta^2=1`$. Wegen
$`A=(x_1-x_2)u+x_2s`$ und $`L=(x_1+x_2)u-x_2s`$ ist $`A-\alpha L=x_2s(1+\alpha)=2x_1x_2s/(x_1+x_2)`$ und
$`\beta K=4x_1x_2y_1y_2/(x_1+x_2)`$, also $`A-\alpha L-\beta K=\frac{2x_1x_2}{x_1+x_2}(y_1-y_2)^2=:m`$. Die
Cauchy–Schwarz-Ungleichung in $`\mathbb R^2`$ liefert $`\sqrt{L^2+K^2}\ge\alpha L+\beta K`$, also
```math
z=A-\sqrt{L^2+K^2}\ \le\ A-\alpha L-\beta K=m,
```
mit Gleichheit genau dann, wenn $`(L,K)`$ ein nichtnegatives Vielfaches von $`(\alpha,\beta)`$ ist, d. h. (da $`\beta>0`$, $`K\ge0`$)
genau dann, wenn $`\beta L=\alpha K`$, d. h. $`(x_1+x_2)u-x_2s=(x_1-x_2)y_1y_2`$. Das legt $`u`$ fest, und
$`u-y_2^2=(y_1-y_2)(x_1y_2+x_2y_1)/(x_1+x_2)`$; für $`y_1>y_2`$ liefert $`u=\gamma y_1^2+(1-\gamma)y_2^2`$ dann $`\gamma=\gamma_3\in[0,\frac12]`$.

*Reelle Maximierer.* Man nehme für $`e_i`$ die Standardbasis, für $`f_1,f_2`$ die gedrehte reelle Basis mit $`\cos^2t=\gamma_1`$ (bzw. $`\gamma_3`$)
und in (a) $`U=VW^{\mathrm T}`$ aus einer reellen Singulärwertzerlegung $`YX=W\Sigma V^{\mathrm T}`$ (durch Diagonalisieren
der reellen symmetrischen Matrix $`(YX)^{\mathrm T}YX`$), sodass $`\mathrm{Tr}(UYX)=\mathrm{Tr}\,\Sigma=\|YX\|_1`$. $`\square`$

<a id="r-11-16"></a>

**Korollar 11.16 (Einhüllende von Exp. 1 und Exp. 3 für Blöcke der Größe $`2`$).** Seien alle $`k_c=2`$ (also $`D=2C`$),
$`\delta\in[0,1)`$, $`\delta'\in[0,1]`$, $`\ell=\ell(\delta)`$, $`\ell'=\ell(\delta')`$ (mit $`\ell(1)=0`$), und $`Y`$
durchlaufe die positiven Operatoren mit $`\lambda_{\max}(Y)=1`$ und $`\Delta(Y)=\delta'`$, d. h. $`\lambda_{\min}(Y)=\ell'`$.

(a) Exp. 1, über $`X\in\mathcal P_\delta`$, solche $`Y`$ und unitäre $`U`$:
$`\ \sup z=C\,\frac{\delta^2\delta'^2}{(1+\delta)(1+\delta')}`$ und $`\inf z=-(2C-1+\ell\ell')`$.

(b) Exp. 3, über $`X=X^*`$ mit $`|X|\in\mathcal P_\delta`$ und solche $`Y`$:
$`\ \sup z=(C-\delta)\Big(\frac{2\delta'}{1+\delta'}\Big)^2`$ und $`\inf z=0`$.

Alle vier Werte werden angenommen. Folglich gilt mit dem anderen Kontrast frei in $`[0,d]`$, $`d\in(0,1)`$: Die obere
Einhüllende von Exp. 1 ist $`C\delta^2d^2/((1+\delta)(1+d))`$ entlang $`\Delta(X)=\delta`$ und (wegen der Symmetrie von (a))
$`C\delta'^2d^2/((1+\delta')(1+d))`$ entlang $`\Delta(Y)=\delta'`$; die von Exp. 3 ist $`(C-\delta)\big(2d/(1+d)\big)^2`$ entlang $`\Delta(|X|)=\delta`$ und
$`C\big(2\delta'/(1+\delta')\big)^2`$ entlang $`\Delta(Y)=\delta'`$, für jedes $`d`$. Für $`d\to1`$ wachsen die ersten beiden gegen
$`C\delta^2/(2(1+\delta))`$ und $`C-\delta`$, die Werte bei $`\Delta(Y)=1`$ (singuläres $`Y`$).

*Beweis.* Seien $`x_{c,1}\ge x_{c,2}`$ die Eigenwerte von $`X_c`$ (in (b): von $`|X|_c`$) und $`y_{c,1}\ge y_{c,2}`$ die von
$`Y_c`$; sie liegen in $`[\ell,1]`$ und $`[\ell',1]`$, und $`\ell`$ ist der kleinere in seinem Block. Es ist
$`1-\ell=2\delta/(1+\delta)`$, $`1+\ell=2/(1+\delta)`$, also $`(1-\ell)^2/(1+\ell)=2\delta^2/(1+\delta)`$ und
$`2\ell/(1+\ell)=1-\delta`$; ebenso für $`\ell'`$.

(a) *Supremum.* Nach der Dreiecksungleichung und Satz 1.9 in jedem Block ist
$`|\mathrm{Tr}(XUY)|\le\sum_c|\mathrm{Tr}(X_cU_cY_c)|\le\sum_c\|Y_cX_c\|_1`$, also nach Proposition 11.15(a) (ein Block mit
$`Y_c=0`$ trägt $`-\mathrm{Tr}(X_cY_c)=0`$ bei)
```math
z\le\sum_c\frac{x_{c,1}y_{c,1}}4\,g\Big(\frac{x_{c,2}}{x_{c,1}}\Big)\,g\Big(\frac{y_{c,2}}{y_{c,1}}\Big),\qquad g(r)=\frac{(1-r)^2}{1+r}.
```
$`g`$ fällt auf $`[0,1]`$, $`x_{c,2}/x_{c,1}\ge x_{c,2}\ge\ell`$, $`y_{c,2}/y_{c,1}\ge\ell'`$ und $`x_{c,1},y_{c,1}\le1`$, also ist
jeder Summand höchstens $`g(\ell)g(\ell')/4=\delta^2\delta'^2/((1+\delta)(1+\delta'))`$. Gleichheit: in jedem Block
$`X_c=\mathrm{diag}(1,\ell)`$, $`Y_c`$ mit Eigenwerten $`1,\ell'`$ bei $`\gamma=\gamma_1`$ und $`U_c=W_c^*`$, $`Y_cX_c=W_c|Y_cX_c|`$;
dann ist $`\mathrm{Tr}(X_cU_cY_c)=\|Y_cX_c\|_1\ge0`$ in jedem Block, und die Beträge addieren sich.
*Infimum.* $`z\ge-\mathrm{Tr}(XY)`$, mit Gleichheit für $`U=\bigoplus_cU_c`$, $`\mathrm{Tr}(Y_cX_cU_c)=0`$ (Lemma 11.13). Nach
Korollar 11.5(a) ist $`\mathrm{Tr}(XY)\le\sum_c(x_{c,1}y_{c,1}+x_{c,2}y_{c,2})`$. Liegen $`\ell`$ und $`\ell'`$ im selben
Block, so trägt dieser höchstens $`1+\ell\ell'`$ bei und jeder andere höchstens $`2`$; andernfalls tragen die beiden Blöcke
höchstens $`1+\ell`$ und $`1+\ell'`$ bei, und $`2C-2+\ell+\ell'=2C-1+\ell\ell'-(1-\ell)(1-\ell')`$. Also ist
$`\mathrm{Tr}(XY)\le2C-1+\ell\ell'`$, mit Gleichheit für $`X=\mathrm{diag}(1,\ell)\oplus1\oplus\dots\oplus1`$ und
$`Y=\mathrm{diag}(1,\ell')\oplus1\oplus\dots\oplus1`$.

(b) Spuren, Spurnormen und $`|\cdot|`$ wirken blockweise (Aufgabe 3), also $`z=\sum_cz_c`$, und nach Proposition 11.15(b) ist
$`z_c\le H_c(y_{c,1}-y_{c,2})^2\le H_c(1-\ell')^2`$ mit $`H_c=2x_{c,1}x_{c,2}/(x_{c,1}+x_{c,2})\le x_{c,1}\le1`$. In einem Block
$`c_-`$, in dem $`|X|`$ den Eigenwert $`\ell`$ hat, ist $`x_{c_-,2}=\ell`$ und $`H_{c_-}=2\ell/(1+\ell/x_{c_-,1})\le2\ell/(1+\ell)=1-\delta`$.
Daher ist $`z\le(C-1+1-\delta)(1-\ell')^2`$. Gleichheit: $`X_1=\mathrm{diag}(1,-\ell)`$, $`X_c=\mathrm{diag}(1,-1)`$ für $`c\ge2`$,
und jedes $`Y_c`$ mit Eigenwerten $`1,\ell'`$ beim $`\gamma=\gamma_3`$ seines Blocks (für $`c\ge2`$ ist $`\gamma_3=\frac12`$). Das Infimum ist
Proposition 11.9, mit Gleichheit für kommutierende $`X`$, $`Y`$.

*Freier Kontrast.* Das Supremum in (a) wächst in $`\delta'`$ (streng, falls $`\delta>0`$), das in (b) wächst streng in
$`\delta'`$ und fällt in $`\delta`$. Über den freien Kontrast in $`[0,d]`$ werden die Suprema daher bei $`d`$ angenommen
(entlang $`\Delta(Y)`$ in Exp. 3: bei $`\Delta(|X|)=0`$), und die Grenzwerte für $`d\to1`$ werden nur bei $`\Delta(Y)=1`$
erreicht (in Exp. 1 falls $`\delta>0`$). $`\square`$

Für $`C=1`$ und $`d=0.999`$ sind dies die geschlossenen Formen, die die Begleitarbeit zuerst numerisch beobachtet hat (in
reellen Algebren). Für jedes $`C`$ sind die oberen Einhüllenden die Werte direkter Summen von $`C`$ Lösungen auf $`M_2`$
($`C-\delta=(1-\delta)+(C-1)\cdot1`$ in Exp. 3): Auf $`M_2\oplus M_2`$ ist die obige Ungleichung für direkte Summen eine
Gleichheit. Für Blöcke der Größe $`k_c\ge3`$ kennen wir keine geschlossene Form dieser oberen Einhüllenden.

*Rechnung:* die Maximierer aus Proposition 11.15, die Einhüllenden aus Korollar 11.16 für $`C=1,2`$ und ihr Grenzwert
$`d\to1`$ sowie eine Multistart-Suche im Vergleich dazu — [Tutorial, Abschnitt 11.5](../../tutorial/de/11_optimisation.md#s-11-5).

<a id="s-11-6"></a>

## 11.6 Womit man vergleichen sollte

*In der Bibliothek:* der Vergleich mit anderer Software, darunter SDP-Löser, deren konvexe Relaxierungen ein Supremum nach oben beschränken — [Tutorial, Abschnitt 11.6](../../tutorial/de/11_optimisation.md#s-11-6).

<a id="x-11"></a>

## Aufgaben

1. (**Cayley-Transformation.**) Man zeige für $`S^*=-S`$, dass $`C(S)=(1-S)(1+S)^{-1}`$ unitär ist und dass
   $`S\mapsto C(S)`$ eine Bijektion auf die unitären Matrizen ohne Eigenwert $`-1`$ ist. Anders als die
   Exponentialabbildung (Lemma 11.1) verfehlt eine Cayley-Karte $`U_0C(S)`$ also $`-U_0`$.
2. (**Reelles Haar-Maß.**) Für Haar-verteiltes $`O`$ auf $`O(k)`$ (QR einer reellen Ginibre-Matrix) und reelles $`A`$ zeige
   man $`\mathbb E\mathrm{Tr}(AO)^2=\|A\|_2^2/k`$ und prüfe dies mit einer reellen Algebra.
3. (**Direkte Summen.**) Man zeige, dass $`z`$ aus Exp. 3 additiv über direkte Summen ist und dass
   $`\Delta(A\oplus A')=\max(\Delta(A),\Delta(A'))`$ für $`A,A'\ge0`$ mit $`\lambda_{\max}=1`$ gilt. Man folgere die
   Ungleichung für direkte Summen aus Abschnitt 11.5 (mit $`\Delta(Y)`$ frei in einem festen Intervall $`[0,b]`$).

<details>
<summary>Lösungen</summary>

1. $`S`$ ist normal und $`iS`$ selbstadjungiert, also hat $`S`$ rein imaginäre Eigenwerte (Korollar 1.4(2)), und $`1\pm S`$
   sind invertierbar; sie kommutieren, und $`(1+S)^*=1-S`$. Daher ist $`C^*C=(1-S)^{-1}(1+S)(1-S)(1+S)^{-1}=1`$, und
   $`1+C=2(1+S)^{-1}`$ ist invertierbar, also $`-1\notin\mathrm{spec}\,C`$. Umgekehrt ist für eine unitäre Matrix $`U`$
   mit invertierbarem $`1+U`$ die Matrix $`S=(1-U)(1+U)^{-1}`$ schiefadjungiert:
   $`S^*=(1+U^*)^{-1}(1-U^*)=(U+1)^{-1}U\,U^*(U-1)=(U+1)^{-1}(U-1)=-S`$, da $`U\pm1`$ kommutieren. Für kommutierende $`C`$
   und $`S`$ ist $`C(1+S)=1-S`$ äquivalent zu $`S(1+C)=1-C`$ (beides besagt $`C+CS+S=1`$). Nun kommutiert $`C(S)`$ mit $`S`$, und
   $`S=(1-U)(1+U)^{-1}`$ kommutiert mit $`U`$; also ist $`U=C(S)`$ genau dann, wenn $`S=(1-U)(1+U)^{-1}`$: Die beiden
   Abbildungen sind zueinander invers.

2. Nach Satz B.17 (auch Lemma C.2) ist der vorzeichenkorrigierte QR-Faktor $`O`$ einer reellen Gaußschen Matrix, den
   die Bibliothek für eine reelle Algebra zieht, Haar-verteilt auf $`O(k)`$; also hat $`VO`$ für $`V\in O(k)`$ dieselbe
   Verteilung wie $`O`$.
   Linksmultiplikation mit $`\mathrm{diag}(\pm1)`$ liefert $`\mathbb E[O_{ij}O_{lm}]=0`$ für $`i\ne l`$, Permutationen
   zeigen, dass $`\mathbb E[O_{ij}O_{im}]`$ nicht von $`i`$ abhängt, und $`\sum_iO_{ij}O_{im}=\delta_{jm}`$; also ist
   $`\mathbb E[O_{ij}O_{lm}]=\delta_{il}\delta_{jm}/k`$ und $`\mathbb E\mathrm{Tr}(AO)^2=\sum_{ij}A_{ji}^2/k`$.

3. Spuren und Spurnormen blockdiagonaler Operatoren sind Summen über die Blöcke, und $`|X\oplus X'|=|X|\oplus|X'|`$;
   $`\lambda_{\min}(A\oplus A')=\min(\lambda_{\min}(A),\lambda_{\min}(A'))`$, und $`\Delta`$ fällt mit wachsendem
   $`\lambda_{\min}`$. Ist also $`(X,Y)`$ zulässig bei $`\delta`$ und $`(X',Y')`$ bei $`\delta'\le\delta`$, so ist
   $`(X\oplus X',Y\oplus Y')`$ zulässig bei $`\delta`$, mit $`z=z(X,Y)+z(X',Y')`$; man gehe zum Supremum über beide Paare und
   über $`\delta'\le\delta`$ über.

</details>

*Rechnung:* die Aufgaben 1 und 2 numerisch — [Tutorial, Aufgaben](../../tutorial/de/11_optimisation.md#x-11).

**Weiter:** [Kapitel 12](12_equivariant_learning.md) behandelt äquivariantes maschinelles Lernen: Schichten, die mit einer
Symmetrie kommutieren, sind Elemente einer Kommutante (Kapitel 2), aufgebaut aus denselben gebündelten,
differenzierbaren Operatoren.
