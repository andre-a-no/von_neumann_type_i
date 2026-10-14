[English](../06_quantum_channels.md) · [Español](../es/06_quantum_channels.md) · [Français](../fr/06_quantum_channels.md) · **Deutsch** · [中文](../zh/06_quantum_channels.md) · [日本語](../ja/06_quantum_channels.md) · [Русский](../ru/06_quantum_channels.md)

# 6. Quantenkanäle

Ein Quantenkanal ist eine vollständig positive, spurerhaltende lineare Abbildung. Die Konjugation
mit einer unitären Matrix, die bedingte Erwartung auf das Zentrum, die durch eine Isometrie $`V`$
gegebene Abbildung $`\rho\mapsto\mathrm{Tr}_E(V\rho V^*)`$ und Abbildungen, die Gewicht von einem
Sektor in einen anderen verschieben, sind Beispiele; ihre physikalische Lesart (eine Zeitentwicklung,
eine Messung, deren Ergebnis vergessen wird, eine Umgebung, der Verlust eines Teilchens) wird nur in
den Kästen *Physikalische Bedeutung* besprochen. Dieses Kapitel beweist die Standardbeschreibungen solcher Abbildungen – Choi-Matrizen,
Kraus-Operatoren, die Stinespring-Dilatation – für Abbildungen zwischen endlichdimensionalen
von-Neumann-Algebren und untersucht dann, was an direkten Summen besonders ist: Abbildungen
zwischen Sektoren, duale Abbildungen bezüglich der Spuren aus Kapitel 3 und den klassischen Teil
eines Kanals. Die Kapitel 2 und 3 haben die Kanäle der Bibliothek bereits verwendet; hier folgt ihre
Theorie.

**Sie lernen**, warum Positivität nicht genügt, wie man zwischen Choi-Matrizen, Kraus-Operatoren,
Superoperatoren und Isometrien umrechnet, warum „spurerhaltend“ von der Spur abhängt, sobald eine
Abbildung den Sektor wechselt, welche Operatoren ein unitaler Kanal festlässt und was ein zufälliger
Kanal im Mittel mit einem reinen Zustand macht.

Die numerischen Überprüfungen dieses Kapitels stehen in [Tutorial 6](../../tutorial/de/06_quantum_channels.md).

Im Folgenden seien $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$ und
$`\mathcal N=\bigoplus_{d=1}^{C'}M_{k'_d}(\mathbb C)`$ mit minimalen zentralen Projektionen $`P_c`$ und
$`Q_d`$; $`A=\bigoplus_cA_c`$; $`M_{k'\times k}`$ ist der Raum der komplexen $`k'\times k`$-Matrizen;
$`E_{ij}=|i\rangle\langle j|`$ sind die Matrixeinheiten; $`\mathrm{Tr}`$ ist die ungewichtete Spur
$`\sum_c\mathrm{Tr}\,A_c`$ und $`\mathrm{tr}_w=\sum_cw_c\mathrm{Tr}\,A_c`$ eine treue Spur
(Satz 3.1). Blockmatrizen wirken auf $`\mathbb C^n\otimes\mathbb C^k`$, der Blockindex zuerst.

<a id="s-6-1"></a>

## 6.1 Vollständig positive Abbildungen und der Satz von Choi

Eine lineare Abbildung $`\Phi:\mathcal M\to\mathcal N`$ heißt **positiv**, wenn aus $`A\ge0`$ stets
$`\Phi(A)\ge0`$ folgt. Die $`n\times n`$-Matrizen $`[A_{ij}]`$ mit Einträgen in $`\mathcal M`$ bilden wieder
eine direkte Summe von Matrixalgebren, $`M_n(\mathcal M)\cong\bigoplus_cM_n(M_{k_c})`$ vermöge
$`[A_{ij}]\mapsto([(A_{ij})_c]_{ij})_c`$; daher bedeutet $`[A_{ij}]\ge0`$, dass jede Blockmatrix
$`[(A_{ij})_c]\in M_{nk_c}`$ positiv ist. $`\Phi`$ heißt **$`n`$-positiv**, wenn
$`\Phi_n:[A_{ij}]\mapsto[\Phi(A_{ij})]`$ positiv ist, und **vollständig positiv** (CP), wenn es für jedes
$`n`$ $`n`$-positiv ist. Es heißt **spurerhaltend** (TP), wenn $`\mathrm{Tr}\,\Phi(A)=\mathrm{Tr}\,A`$
gilt, **unital**, wenn $`\Phi(1)=1`$, und ein **Kanal**, wenn es CP und TP ist. Wir verwenden
wiederholt, dass aus $`X\ge0`$ für jede rechteckige Matrix $`T`$ auch $`TXT^*\ge0`$ folgt, da
$`\langle y,TXT^*y\rangle=\langle T^*y,XT^*y\rangle`$.

<a id="r-6-1"></a>

**Proposition 6.1.** Für $`k\ge2`$ ist die Transposition $`\theta(A)=A^T`$ auf $`M_k`$ positiv, aber nicht
$`2`$-positiv.

*Beweis.* Für $`A\ge0`$ ist $`\langle x,A^Tx\rangle=\sum_{ij}\bar x_iA_{ji}x_j=\langle\bar x,A\bar x\rangle\ge0`$.
Mit $`\Omega=e_1\otimes e_1+e_2\otimes e_2\in\mathbb C^2\otimes\mathbb C^k`$ ist
$`|\Omega\rangle\langle\Omega|=[E_{ij}]_{i,j=1}^2\ge0`$, und $`\theta_2`$ bildet dies auf
$`F=[E_{ji}]_{i,j=1}^2`$ ab, das $`e_1\otimes e_2-e_2\otimes e_1`$ auf sein Negatives abbildet. $`\square`$

*Rechnung:* Die Choi-Matrix der Transposition hat den Eigenwert $`-1`$ — [Tutorial, Abschnitt 6.1](../../tutorial/de/06_quantum_channels.md#s-6-1).

<a id="r-6-2"></a>

**Satz 6.2 (Choi).** Für eine lineare Abbildung $`\Phi:M_k\to M_{k'}`$ sind äquivalent:
(a) $`\Phi`$ ist CP; (b) $`\Phi`$ ist $`k`$-positiv; (c) die **Choi-Matrix**
```math
J(\Phi)=\sum_{i,j=1}^kE_{ij}\otimes\Phi(E_{ij})=[\Phi(E_{ij})]_{i,j=1}^k\in M_k(M_{k'})
```
ist positiv; (d) es gibt $`K_1,\dots,K_r\in M_{k'\times k}`$ (**Kraus-Operatoren**) mit
$`\Phi(A)=\sum_lK_lAK_l^*`$, und man kann $`r=\mathrm{rank}\,J(\Phi)\le kk'`$ wählen.

*Beweis.* (a)$`\Rightarrow`$(b) ist trivial. (b)$`\Rightarrow`$(c): Für
$`\Omega=\sum_ie_i\otimes e_i\in\mathbb C^k\otimes\mathbb C^k`$ ist
$`|\Omega\rangle\langle\Omega|=[E_{ij}]\ge0`$ und $`J(\Phi)=\Phi_k([E_{ij}])`$. (c)$`\Rightarrow`$(d): Nach
Satz 1.3 ist $`J(\Phi)=\sum_{l=1}^rv_lv_l^*`$ mit $`r=\mathrm{rank}\,J(\Phi)`$ (die Eigenvektoren zu den
von null verschiedenen Eigenwerten, skaliert mit deren Quadratwurzeln). Jedes
$`v\in\mathbb C^k\otimes\mathbb C^{k'}`$ lässt sich als $`v=\sum_ie_i\otimes Ke_i`$ mit genau einem
$`K\in M_{k'\times k}`$ schreiben ($`Ke_i`$ ist der $`i`$-te Block von $`v`$), und der $`(i,j)`$-Block von
$`vv^*`$ ist $`Ke_ie_j^*K^*=KE_{ij}K^*`$. Schreibt man jedes $`v_l`$ so mit $`K_l`$ und vergleicht die
Blöcke, erhält man $`\Phi(E_{ij})=\sum_lK_lE_{ij}K_l^*`$ und daraus (d) wegen der Linearität.
(d)$`\Rightarrow`$(a): $`\Phi_n([A_{ij}])=\sum_l(1_n\otimes K_l)[A_{ij}](1_n\otimes K_l)^*\ge0`$. $`\square`$

<a id="r-6-3"></a>

**Korollar 6.3.** Ist $`\Phi(A)=\sum_{l=1}^sL_lAL_l^*`$, so gilt $`s\ge\mathrm{rank}\,J(\Phi)`$. Die
minimale Anzahl von Kraus-Operatoren, der **Kraus-Rang**, ist also $`\mathrm{rank}\,J(\Phi)`$.

*Beweis.* Liest man die letzte Rechnung rückwärts, so ist $`J(\Phi)=\sum_{l=1}^sw_lw_l^*`$ mit
$`w_l=\sum_ie_i\otimes L_le_i`$, eine Summe von $`s`$ Matrizen vom Rang $`\le1`$. $`\square`$

*Rechnung:* Die Choi-Matrix und der Kraus-Rang eines zufälligen Kanals — [Tutorial, Abschnitt 6.1](../../tutorial/de/06_quantum_channels.md#s-6-1).

Eine Abbildung zwischen direkten Summen ist eine Matrix von Abbildungen zwischen Blöcken. Es setze
$`\iota_c:M_{k_c}\to\mathcal M`$ eine Matrix in den Block $`c`$, $`\pi_d:\mathcal N\to M_{k'_d}`$ greife den
Block $`d`$ heraus, und es sei $`\Phi_{dc}=\pi_d\circ\Phi\circ\iota_c:M_{k_c}\to M_{k'_d}`$. Wegen
$`A=\sum_c\iota_c(A_c)`$ ist $`\Phi(A)_d=\sum_c\Phi_{dc}(A_c)`$.

<a id="r-6-4"></a>

**Proposition 6.4.** $`\Phi:\mathcal M\to\mathcal N`$ ist genau dann $`n`$-positiv, wenn jedes
$`\Phi_{dc}`$ es ist. Folglich ist $`\Phi`$ genau dann CP, wenn es $`K^{dc}_l\in M_{k'_d\times k_c}`$,
$`l=1,\dots,r`$, gibt mit
```math
\Phi(A)_d=\sum_c\sum_{l=1}^rK^{dc}_l\,A_c\,K^{dc*}_l ,
```
und $`r=\max_{d,c}\mathrm{rank}\,J(\Phi_{dc})`$ genügt. $`\Phi`$ ist genau dann TP, wenn
$`\sum_{d,l}K^{dc*}_lK^{dc}_l=1_{k_c}`$ für jedes $`c`$ gilt.

*Beweis.* Sind alle $`\Phi_{dc}`$ $`n`$-positiv und ist $`[A_{ij}]\ge0`$, so ist jedes
$`X^c=[(A_{ij})_c]\ge0`$, und der Block $`d`$ von $`\Phi_n([A_{ij}])`$ ist $`\sum_c(\Phi_{dc})_n(X^c)\ge0`$.
Umgekehrt ist für $`[X_{ij}]\ge0`$ in $`M_n(M_{k_c})`$ auch $`[\iota_c(X_{ij})]\ge0`$ in $`M_n(\mathcal M)`$
(sein Block $`c`$ ist $`[X_{ij}]`$, die übrigen verschwinden), und der Block $`d`$ seines Bildes ist
$`[\Phi_{dc}(X_{ij})]`$. Die Kraus-Form folgt aus Satz 6.2 für jedes Paar $`(d,c)`$, aufgefüllt mit
Nulloperatoren auf ein gemeinsames $`r`$. Schließlich ist wegen der Zyklizität
$`\mathrm{Tr}\,\Phi(A)=\sum_c\mathrm{Tr}\big(A_c\sum_{d,l}K^{dc*}_lK^{dc}_l\big)`$, und dies ist
genau dann für alle $`A`$ gleich $`\sum_c\mathrm{Tr}\,A_c`$, wenn jede dieser Summen gleich $`1_{k_c}`$
ist (man setze Matrixeinheiten ein: $`\mathrm{Tr}(E_{ij}X)=X_{ji}`$). $`\square`$

Eine lineare Abbildung $`\Phi:M_k\to M_{k'}`$ erfüllt $`\Phi(X^*)=\Phi(X)^*`$ für alle $`X`$ genau dann, wenn
$`J(\Phi)^*=J(\Phi)`$: Der $`(i,j)`$-Block von $`J^*`$ ist $`\Phi(E_{ji})^*`$, und dieser ist wegen der Linearität
genau dann für alle $`i,j`$ gleich $`\Phi(E_{ij})`$, wenn $`\Phi`$ Adjungierte erhält.

*In der Bibliothek:* Sektorerhaltende Abbildungen sind ein `Channel`, die allgemeinen Abbildungen aus Proposition 6.4 ein `InterSectorChannel`; `Channel.choi()`, `Channel.superoperator()` und `Channel.from_superoperator` rechnen zwischen den Darstellungen um — [Tutorial, Abschnitt 6.1](../../tutorial/de/06_quantum_channels.md#s-6-1).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Ein System ist nie allein: Es kann mit einer unberührten Ancilla mit
> $`n`$ Niveaus korreliert sein. Ein physikalischer Prozess wirkt auf das Paar als
> $`\Phi_n=\mathrm{id}_n\otimes\Phi`$ und muss Zustände auf Zustände abbilden; also muss $`\Phi`$
> vollständig positiv sein. Die Transposition ist positiv, bildet aber den verschränkten Zustand
> $`|\Omega\rangle\langle\Omega|/2`$ auf einen Operator mit Eigenwert $`-1/2`$ ab. Die Choi-Matrix ist das
> Bild des maximal verschränkten Zustands: Ein einziger Zustand entscheidet über vollständige
> Positivität.

<a id="s-6-2"></a>

## 6.2 Stinespring-Dilatation und die Freiheit in den Kraus-Operatoren

Für $`\Phi(A)=\sum_{l=1}^rK_lAK_l^*`$ von $`M_k`$ nach $`M_{k'}`$ setze man $`\Phi^*(B)=\sum_lK_l^*BK_l`$; wegen
der Zyklizität ist $`\mathrm{Tr}(\Phi(A)B)=\mathrm{Tr}(A\Phi^*(B))`$, also ist $`\Phi^*`$ die
duale Abbildung zu $`\Phi`$ bezüglich der Spur (Abschnitt 6.3) und hängt nicht von der Wahl der $`K_l`$ ab
(ist $`\mathrm{Tr}(AX)=\mathrm{Tr}(AX')`$ für alle $`A`$, so ist $`X=X'`$: man nehme $`A=E_{ji}`$).
Auf $`\mathbb C^r\otimes\mathbb C^{k'}`$ (die **Umgebung** $`\mathbb C^r`$ zuerst) sei
$`\mathrm{Tr}_E\,X=\sum_l(e_l^*\otimes1)X(e_l\otimes1)`$.

<a id="r-6-5"></a>

**Satz 6.5 (Stinespring, endlichdimensional).** Sei
$`V=\sum_le_l\otimes K_l:\mathbb C^k\to\mathbb C^r\otimes\mathbb C^{k'}`$. Dann gilt
```math
\Phi(A)=\mathrm{Tr}_E(VAV^*),\qquad\Phi^*(B)=V^*(1_r\otimes B)V,
```
und $`\Phi`$ ist genau dann TP, wenn $`V`$ eine Isometrie ist. Umgekehrt definiert jedes $`V`$ durch die
erste Formel eine CP-Abbildung.

*Beweis.* $`(e_l^*\otimes1)V=K_l`$, also $`\mathrm{Tr}_E(VAV^*)=\sum_lK_lAK_l^*`$; umgekehrt ist
$`\mathrm{Tr}_E(VAV^*)`$ eine Kraus-Form mit $`K_l=(e_l^*\otimes1)V`$. Ferner ist
$`V^*(1\otimes B)V=\sum_{l,l'}\langle e_l,e_{l'}\rangle K_l^*BK_{l'}=\Phi^*(B)`$ und
$`V^*V=\sum_lK_l^*K_l`$, was genau dann $`1`$ ist, wenn $`\Phi`$ TP ist (Proposition 6.4). $`\square`$

Also ist das kleinste $`r`$, für das $`\Phi(A)=\mathrm{Tr}_E(VAV^*)`$ mit einem
$`V:\mathbb C^k\to\mathbb C^r\otimes\mathbb C^{k'}`$ gilt, der Kraus-Rang $`\mathrm{rank}\,J(\Phi)`$: Jedes
solche $`V`$ liefert eine Kraus-Form mit $`r`$ Operatoren, also ist $`r\ge\mathrm{rank}\,J(\Phi)`$
(Korollar 6.3), und $`r=\mathrm{rank}\,J(\Phi)`$ wird nach Satz 6.2(d) erreicht.

*Rechnung:* Satz 6.5 für einen zufälligen Kanal — [Tutorial, Abschnitt 6.2](../../tutorial/de/06_quantum_channels.md#s-6-2).

Die Formel $`\Phi^*(B)=V^*\pi(B)V`$ mit der Darstellung $`\pi(B)=1_r\otimes B`$ ist die Form, in der
W. F. Stinespring (Proc. AMS 6 (1955) 211) den Satz formuliert hat, in einem weit allgemeineren
Rahmen als dem hier bewiesenen endlichdimensionalen (historische Zuordnung; weiterführende Literatur:
V. Paulsen, *Completely Bounded Maps and Operator Algebras*, Cambridge 2002, Kap. 4). Für eine
Abbildung zwischen direkten Summen sei $`\iota'_d`$ die Inklusion des $`d`$-ten Summanden von
$`\bigoplus_d\mathbb C^{k'_d}`$ und
$`V_c=\sum_{d,l}e_l\otimes\iota'_dK^{dc}_l:\mathbb C^{k_c}\to\mathbb C^r\otimes\bigoplus_d\mathbb C^{k'_d}`$.
Dann ist $`(e_l^*\otimes1)V_c=\sum_d\iota'_dK^{dc}_l`$, also ist der $`(d,d)`$-Block von
$`\sum_c\mathrm{Tr}_E(V_cA_cV_c^*)=\sum_{c,l}\big(\sum_d\iota'_dK^{dc}_l\big)A_c\big(\sum_{d'}\iota'_{d'}K^{d'c}_l\big)^*`$
gleich $`\sum_{c,l}K^{dc}_lA_cK^{dc*}_l=\Phi(A)_d`$; und wegen $`\iota'^*_d\iota'_{d'}=\delta_{dd'}1`$ ist
$`V_c^*V_c=\sum_{d,l}K^{dc*}_lK^{dc}_l`$, also ist $`\Phi`$ genau dann TP, wenn jedes $`V_c`$ eine Isometrie
ist (Proposition 6.4). So ziehen `random_channel` und `random_inter_sector_channel` Kanäle
(Abschnitt 6.5).

*In der Bibliothek:* `channels.random_channel` und `channels.random_inter_sector_channel` — [Tutorial, Abschnitt 6.5](../../tutorial/de/06_quantum_channels.md#s-6-5).

<a id="r-6-6"></a>

**Lemma 6.6.** Sind $`X,Y\in M_{p\times n}`$ mit $`XX^*=YY^*`$, so ist $`Y=XW`$ mit einem unitären
$`W\in M_n`$.

*Beweis.* $`\|X^*z\|^2=\langle z,XX^*z\rangle=\|Y^*z\|^2`$ für alle $`z`$, also ist $`X^*z\mapsto Y^*z`$ eine
wohldefinierte lineare Isometrie vom Bild von $`X^*`$ auf das Bild von $`Y^*`$. Wie im Beweis von Satz 1.8
lässt sie sich zu einem unitären $`U`$ mit $`UX^*=Y^*`$ fortsetzen; man nehme $`W=U^*`$. $`\square`$

<a id="r-6-7"></a>

**Satz 6.7 (unitäre Freiheit).** Seien $`K_1,\dots,K_n`$ und $`L_1,\dots,L_n`$ in $`M_{k'\times k}`$ (die
kürzere Familie wird mit Nullen aufgefüllt). Dann gilt $`\sum_lK_lAK_l^*=\sum_jL_jAL_j^*`$ für alle $`A`$
genau dann, wenn $`L_j=\sum_lu_{jl}K_l`$ mit einem unitären $`u\in M_n`$.

*Beweis.* Ist $`L_j=\sum_lu_{jl}K_l`$, so gilt
$`\sum_jL_jAL_j^*=\sum_{l,l'}(u^*u)_{l'l}K_lAK_{l'}^*=\sum_lK_lAK_l^*`$. Umgekehrt haben gleiche
Abbildungen gleiche Choi-Matrizen, $`\sum_lv_lv_l^*=\sum_jw_jw_j^*`$ mit $`v_l=\sum_ie_i\otimes K_le_i`$
und $`w_j=\sum_ie_i\otimes L_je_i`$. Lemma 6.6, angewandt auf die Matrizen mit den Spalten $`v_l`$ und
$`w_j`$, liefert $`w_j=\sum_lW_{lj}v_l`$; da $`K\mapsto\sum_ie_i\otimes Ke_i`$ linear und injektiv ist, folgt
$`L_j=\sum_lW_{lj}K_l`$, und $`u=W^T`$ ist unitär. $`\square`$

*Rechnung:* Satz 6.7 mit $`n=3`$ — [Tutorial, Abschnitt 6.2](../../tutorial/de/06_quantum_channels.md#s-6-2).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Jeder Kanal ist eine isometrische Wechselwirkung des Systems mit einer
> Umgebung in einem festen Zustand, gefolgt vom Ignorieren der Umgebung. Eine Messung der Umgebung in
> der Basis $`e_l`$ liefert das Ergebnis $`l`$ mit Wahrscheinlichkeit $`\mathrm{Tr}\,K_l\rho K_l^*`$ und
> hinterlässt den Zustand $`K_l\rho K_l^*/\mathrm{Tr}(K_l\rho K_l^*)`$: Die Kraus-Operatoren sind
> die möglichen „Quantensprünge“. Eine andere Basis ergibt andere Sprünge, aber denselben Kanal
> (Satz 6.7), und die kleinste Umgebung hat als Dimension den Kraus-Rang, wie nach Satz 6.5 gezeigt.

<a id="s-6-3"></a>

## 6.3 Duale Abbildungen bezüglich einer Spur

**Definition.** Seien $`\mathrm{tr}_w`$, $`\mathrm{tr}_v`$ treue Spuren auf $`\mathcal M`$,
$`\mathcal N`$. Die **duale Abbildung** einer linearen Abbildung $`\Phi:\mathcal M\to\mathcal N`$
bezüglich $`(w,v)`$ ist die lineare Abbildung $`\Phi^\sharp:\mathcal N\to\mathcal M`$ mit
```math
\mathrm{tr}_v\big(\Phi(x)\,y\big)=\mathrm{tr}_w\big(x\,\Phi^\sharp(y)\big)\qquad(x\in\mathcal M,\ y\in\mathcal N).
```
Sie existiert und ist eindeutig: Für festes $`y`$ ist $`x\mapsto\mathrm{tr}_v(\Phi(x)y)`$ ein
lineares Funktional, nach Satz 3.2 also gleich $`\mathrm{tr}_w(\rho x)=\mathrm{tr}_w(x\rho)`$
für genau ein $`\rho=:\Phi^\sharp(y)`$, und die Eindeutigkeit macht $`\Phi^\sharp`$ linear. Für die
ungewichteten Spuren schreiben wir $`\Phi^*`$. Seien $`z_w=\sum_cw_cP_c`$ und $`z_v=\sum_dv_dQ_d`$ (zentral
und invertierbar).

<a id="r-6-8"></a>

**Proposition 6.8.** Hat $`\Phi`$ die Kraus-Operatoren $`K^{dc}_l`$ (Proposition 6.4), so gilt
```math
\Phi^\sharp(y)_c=\sum_d\frac{v_d}{w_c}\sum_lK^{dc*}_l\,y_d\,K^{dc}_l,\qquad\text{i.e.}\qquad
\Phi^\sharp(y)=z_w^{-1}\,\Phi^*(z_v\,y).
```
Insbesondere ist $`\Phi^\sharp`$ CP, mit den Kraus-Operatoren $`(v_d/w_c)^{1/2}K^{dc*}_l`$.

*Beweis.* Wegen der Zyklizität der Spur in jedem Block ist
```math
\mathrm{tr}_v(\Phi(x)y)=\sum_dv_d\sum_{c,l}\mathrm{Tr}\big(K^{dc}_lx_cK^{dc*}_ly_d\big)
=\sum_cw_c\mathrm{Tr}\Big(x_c\sum_d\frac{v_d}{w_c}\sum_lK^{dc*}_ly_dK^{dc}_l\Big),
```
und die Eindeutigkeit liefert die Formel. Für $`v=w=1`$ lautet sie $`\Phi^*(y)_c=\sum_{d,l}K^{dc*}_ly_dK^{dc}_l`$;
die Faktoren $`v_d`$ und $`1/w_c`$ sind die Multiplikation mit $`z_v`$ und $`z_w^{-1}`$. $`\square`$

$`\Phi`$ heißt **spurerhaltend bezüglich $`(w,v)`$**, wenn $`\mathrm{tr}_v\,\Phi(x)=\mathrm{tr}_w(x)`$
für alle $`x`$ gilt.

<a id="r-6-9"></a>

**Korollar 6.9.** (a) $`\Phi`$ ist genau dann TP bezüglich $`(w,v)`$, wenn $`\Phi^\sharp(1)=1`$, d. h.
$`\sum_{d,l}(v_d/w_c)K^{dc*}_lK^{dc}_l=1_{k_c}`$ für jedes $`c`$. (b) $`\Phi`$ ist genau dann unital, wenn
$`\Phi^\sharp`$ TP bezüglich $`(v,w)`$ ist. (c) Die duale Abbildung zu $`\Phi^\sharp`$ bezüglich $`(v,w)`$ ist
$`\Phi`$. (d) Ist $`\Phi:\mathcal M\to\mathcal M`$ sektorerhaltend und $`v=w`$, so ist $`\Phi^\sharp=\Phi^*`$
für jedes $`w`$, und $`\Phi`$ ist genau dann für eine Spur TP, wenn es für alle TP ist.

*Beweis.* (a) $`\mathrm{tr}_v\,\Phi(x)=\mathrm{tr}_w(x\Phi^\sharp(1))`$, und
$`\mathrm{tr}_w(x(\Phi^\sharp(1)-1))=0`$ gilt genau dann für alle $`x`$, wenn $`\Phi^\sharp(1)=1`$ (man
nehme $`x=(\Phi^\sharp(1)-1)^*`$; $`\mathrm{tr}_w`$ ist treu); dann verwende man Proposition 6.8.
(b) $`\mathrm{tr}_w\,\Phi^\sharp(y)=\mathrm{tr}_v(\Phi(1)y)`$ ist genau dann für alle $`y`$ gleich
$`\mathrm{tr}_v(y)`$, wenn $`\Phi(1)=1`$. (c) $`\mathrm{tr}_w(\Phi^\sharp(y)x)=\mathrm{tr}_v(y\Phi(x))`$
wegen der Spureigenschaft; man verwende die Eindeutigkeit. (d) Es kommen nur die Paare $`d=c`$ vor, für
die $`v_d/w_c=1`$ ist. $`\square`$

**Dichten im Heisenberg-Bild.** Sei $`\Psi:\mathcal N\to\mathcal M`$ linear, $`\omega`$ ein lineares
Funktional auf $`\mathcal M`$ mit der Dichte $`\rho`$ bezüglich $`\mathrm{tr}_w`$ (Satz 3.2) und
$`\Psi^\sharp:\mathcal M\to\mathcal N`$ die duale Abbildung zu $`\Psi`$ bezüglich $`(v,w)`$. Dann ist
$`\omega(\Psi(y))=\mathrm{tr}_w(\Psi(y)\rho)=\mathrm{tr}_v(y\Psi^\sharp(\rho))=\mathrm{tr}_v(\Psi^\sharp(\rho)y)`$,
also hat $`\omega\circ\Psi`$ die Dichte $`\Psi^\sharp(\rho)`$ bezüglich $`\mathrm{tr}_v`$; ist $`\Psi`$
unital, so ist $`\Psi^\sharp`$ nach Korollar 6.9(b) TP bezüglich $`(w,v)`$ (mit vertauschten Rollen von
$`\mathcal M`$ und $`\mathcal N`$).

Für Abbildungen, die den Sektor wechseln, gilt Korollar 6.9(d) nicht. Sei
$`\mathcal M=\mathcal N=M_1\oplus M_2`$ (Sektoren $`0`$ und $`1`$) und $`\Phi`$ die als **Teilchenverlust**
bezeichnete Abbildung (siehe den Kasten unten), $`\Phi(x)=(x_0+\mathrm{Tr}\,x_1)\oplus0`$, mit den Kraus-Operatoren
$`K^{00}=1`$, $`K^{01}_l=\langle l|`$. Sie ist TP für die ungewichteten Spuren. Für $`\mathrm{tr}_w`$ mit
$`w=v=(1,\tfrac12)`$ ist $`\Phi^\sharp(1)=1\oplus2\cdot1_2`$ und
$`\mathrm{tr}_v\,\Phi(x)=x_0+\mathrm{Tr}\,x_1\ne x_0+\tfrac12\mathrm{Tr}\,x_1=\mathrm{tr}_w(x)`$.

*In der Bibliothek:* `Channel.adjoint()` ist $`\Phi^*`$ und `InterSectorChannel.adjoint(trace=...)` ist $`\Phi^\sharp`$; Überprüfungen von Proposition 6.8, Korollar 6.9(a) und des Teilchenverlusts — [Tutorial, Abschnitt 6.3](../../tutorial/de/06_quantum_channels.md#s-6-3).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Das spurunabhängige Objekt ist die Heisenberg-Abbildung
> $`\Psi:\mathcal N\to\mathcal M`$ auf den Observablen: Ein Zustand $`\omega`$ vor dem Prozess wird danach
> zu $`\omega\circ\Psi`$, und $`\Psi`$ ist CP und unital. Dichten hängen von der Spur ab (Satz 3.2), und die
> Dichten transformieren sich mit der dualen Abbildung zu $`\Psi`$ bezüglich der gewählten Spuren (der
> Absatz nach Korollar 6.9). Wirken die Algebren mit Vielfachheiten, so hat die physikalische Spur die
> Gewichte $`w_c=m_c`$ (Abschnitt 2.3). Im Beispiel wird Sektor $`0`$ als „kein Teilchen“ und Sektor $`1`$
> als „ein Teilchen auf zwei Plätzen“ gelesen; die Abbildung entfernt dann das Teilchen. Sie erhält
> die Wahrscheinlichkeit bezüglich der ungewichteten Spur, die Zustände zählt, nicht aber bezüglich
> `Tr_norm`, das die Sektoren gleich gewichtet.

<a id="s-6-4"></a>

## 6.4 Dekohärenz

Die Standardkanäle der Bibliothek wirken innerhalb jedes Blocks ($`p,\gamma\in[0,1]`$):

| Kanal | $`\Phi(\rho)_c`$ | Kraus-Operatoren im Block $`c`$ |
|---|---|---|
| `dephasing_channel` | $`(1-p)\rho_c+p\mathrm{diag}\,\rho_c`$ | $`\sqrt{1-p}\,1`$, $`\sqrt p\,E_{ii}`$ |
| `depolarizing_channel` | $`(1-p)\rho_c+p\mathrm{Tr}(\rho_c)\,1/k_c`$ | $`\sqrt{1-p}\,1`$, $`\sqrt{p/k_c}\,E_{ij}`$ |
| `center_expectation` | $`\mathrm{Tr}(\rho_c)\,1/k_c`$ | $`\sqrt{1/k_c}\,E_{ij}`$ |
| `amplitude_damping_channel` | Niveaus $`i>0`$ zerfallen nach $`0`$ | $`E_{00}+\sqrt{1-\gamma}\sum_{i>0}E_{ii}`$, $`\sqrt\gamma\,E_{0i}`$ |
| `random_mixed_unitary_channel` | $`\sum_jp_jU_j\rho U_j^*`$ | $`\sqrt{p_j}\,U_j`$ |

Die Kraus-Formen prüft man durch direkte Rechnung, z. B.
$`\sum_{ij}\frac1kE_{ij}\rho E_{ji}=\frac1k\sum_{ij}\rho_{jj}E_{ii}`$. Alle fünf sind TP, und alle außer
der Amplitudendämpfung mit $`\gamma>0`$ sind unital (Proposition 6.4: man berechne $`\sum_lK_l^*K_l`$ und
$`\sum_lK_lK_l^*`$; für die Amplitudendämpfung in einem Block der Größe $`k`$ sind sie $`1`$ und
$`1+(k-1)\gamma E_{00}-\gamma\sum_{i>0}E_{ii}`$, und $`\sum_jp_j=1`$ für den gemischt-unitären Kanal).
Explizit bildet die Amplitudendämpfung $`\rho_{00}\mapsto\rho_{00}+\gamma\sum_{i>0}\rho_{ii}`$,
$`\rho_{ii}\mapsto(1-\gamma)\rho_{ii}`$ und $`\rho_{ij}\mapsto(1-\gamma)\rho_{ij}`$ für $`i\ne j`$, beide
positiv, sowie $`\rho_{0j}\mapsto\sqrt{1-\gamma}\,\rho_{0j}`$, $`\rho_{j0}\mapsto\sqrt{1-\gamma}\,\rho_{j0}`$
für $`j>0`$ ab. Die bedingte Erwartung $`E`$ auf das Zentrum (Proposition 2.8) ist die Depolarisierung mit $`p=1`$. Anders als die Messung des Sektors,
$`X\mapsto\sum_cP_cXP_c`$, die auf $`\mathcal M`$ die Identität ist (Abschnitt 3.3), behält sie nur die
$`p_c`$.

*In der Bibliothek:* Die fünf Kanäle der Tabelle sind Funktionen von `channels` — [Tutorial, Abschnitt 6.4](../../tutorial/de/06_quantum_channels.md#s-6-4).

<a id="r-6-10"></a>

**Satz 6.10 (Fixpunkte unitaler Kanäle).** Sei $`\Phi(X)=\sum_lK_lXK_l^*`$ mit $`K_l\in\mathcal M`$
spurerhaltend und unital. Dann gilt $`\Phi(X)=X`$ genau dann, wenn $`X`$ mit allen $`K_l`$ und $`K_l^*`$
kommutiert. Insbesondere bilden die Fixpunkte eine von-Neumann-Unteralgebra von $`\mathcal M`$.

*Beweis.* Kommutiert $`X`$ mit allen $`K_l`$, so ist $`\Phi(X)=X\sum_lK_lK_l^*=X`$. Umgekehrt ist
$`\Phi(X^*)=\Phi(X)^*`$, sodass die selbstadjungierten Teile $`(X+X^*)/2`$ und $`(X-X^*)/2i`$ eines
Fixpunkts Fixpunkte sind; sei also $`X=X^*`$. Mit $`Y_l=K_lX-XK_l`$ ist
```math
\sum_lY_lY_l^*=\sum_l\big(K_lX^2K_l^*-K_lXK_l^*X-XK_lXK_l^*+XK_lK_l^*X\big)
=\Phi(X^2)-\Phi(X)X-X\Phi(X)+X\Phi(1)X=\Phi(X^2)-X^2,
```
wobei $`Y_l^*=XK_l^*-K_l^*X`$, $`\Phi(X)=X`$ und $`\Phi(1)=1`$ benutzt wurden. Die Spur davon ist $`0`$, weil
$`\Phi`$ TP ist. Da $`\mathrm{Tr}\,Y_lY_l^*=\|Y_l\|_2^2\ge0`$, ist
jedes $`Y_l=0`$; durch Übergang zu den Adjungierten kommutiert $`X`$ auch mit $`K_l^*`$. Diese $`X`$ bilden
den Kommutanten einer selbstadjungierten Menge, geschnitten mit $`\mathcal M`$: eine $`*`$-Algebra, die
$`1`$ enthält, also eine von-Neumann-Algebra (Satz 2.2). $`\square`$

Also erhält die Dephasierung ($`p>0`$) genau die Diagonale jedes Blocks, und die Depolarisierung
($`p>0`$) und $`E`$ erhalten genau das Zentrum. Ohne Unitalität gilt der Satz nicht: Die
Amplitudendämpfung lässt $`E_{00}`$ fest, das nicht mit dem Sprung $`E_{01}`$ kommutiert. Für $`0<\gamma\le1`$
und $`k\ge2`$ haben ihre Fixpunkte in einem Block und der Kommutant ihrer Kraus-Operatoren beide die
Dimension $`1`$, aber es sind $`\mathbb CE_{00}`$ und $`\mathbb C1`$: Nach der expliziten Form oben hat ein
Fixpunkt $`\rho_{ij}=0`$ für $`i\ne j`$ und $`\rho_{ii}=0`$ für $`i>0`$ (diese Einträge werden mit
$`\sqrt{1-\gamma}\ne1`$ bzw. $`1-\gamma\ne1`$ multipliziert), und jedes Vielfache von $`E_{00}`$ ist fest;
und was mit allen $`E_{0i}`$, $`E_{i0}`$ ($`i>0`$) kommutiert, kommutiert mit ihren Produkten
$`E_{ij}=E_{i0}E_{0j}`$ und $`E_{00}=E_{01}E_{10}`$, also mit $`M_k`$, und ist daher ein Vielfaches von $`1`$
(Kapitel 2, Aufgabe 1).

*Rechnung:* Fixpunkte und Kommutanten der fünf Kanäle — [Tutorial, Abschnitt 6.4](../../tutorial/de/06_quantum_channels.md#s-6-4).

<a id="r-6-11"></a>

**Satz 6.11 (unitale Kanäle erhöhen die Entropie).** Sei $`\Phi:\mathcal M\to\mathcal N`$ positiv,
spurerhaltend und unital, und seien $`\rho`$ und $`\Phi(\rho)`$ Zustände mit den Eigenwerten $`\lambda_i`$
und $`\mu_j`$. Dann ist $`\mu=\Delta\lambda`$ mit einer doppelt stochastischen Matrix $`\Delta`$, also
$`\sum_jf(\mu_j)\ge\sum_if(\lambda_i)`$ für jede konkave Funktion $`f`$ auf $`[0,1]`$. Insbesondere gilt
$`S(\Phi(\rho))\ge S(\rho)`$ und $`\mathrm{Tr}\,\Phi(\rho)^2\le\mathrm{Tr}\,\rho^2`$.

*Beweis.* Man diagonalisiere jeden Block: $`\rho=\sum_i\lambda_i|e_i\rangle\langle e_i|`$ und
$`\Phi(\rho)=\sum_j\mu_j|f_j\rangle\langle f_j|`$ mit Orthonormalbasen von $`\bigoplus_c\mathbb C^{k_c}`$
und $`\bigoplus_d\mathbb C^{k'_d}`$, deren Vektoren jeweils in einem einzigen Summanden liegen, sodass
$`|e_i\rangle\langle e_i|\in\mathcal M`$. Dann ist $`\mu_j=\langle f_j,\Phi(\rho)f_j\rangle=\sum_i\Delta_{ji}\lambda_i`$
mit $`\Delta_{ji}=\langle f_j,\Phi(|e_i\rangle\langle e_i|)f_j\rangle\ge0`$,
$`\sum_j\Delta_{ji}=\mathrm{Tr}\,\Phi(|e_i\rangle\langle e_i|)=1`$ und
$`\sum_i\Delta_{ji}=\langle f_j,\Phi(1)f_j\rangle=1`$ (also ist $`\Delta`$ quadratisch: man summiere alle
Einträge). Nach der endlichen Jensenschen Ungleichung $`f(\sum_i\Delta_{ji}\lambda_i)\ge\sum_i\Delta_{ji}f(\lambda_i)`$
(durch Induktion nach der Anzahl der Terme aus der Definition der Konkavität, wie in der Lösung von
Aufgabe 2 in Kapitel 3) folgt
$`\sum_jf(\mu_j)\ge\sum_j\sum_i\Delta_{ji}f(\lambda_i)=\sum_if(\lambda_i)`$. Man wähle $`f(t)=-t\log t`$
($`f(0)=0`$) und $`f(t)=-t^2`$, die nach Lemma 3.4 auf $`[0,1]`$ konkav sind: $`t\log t`$ und $`t^2`$ sind
stetig auf $`[0,1]`$ mit den zweiten Ableitungen $`1/t>0`$ und $`2>0`$ auf $`(0,1)`$. $`\square`$

*Rechnung:* Entropie und Reinheit unter den unitalen Kanälen und der Entropieabfall unter Amplitudendämpfung — [Tutorial, Abschnitt 6.4](../../tutorial/de/06_quantum_channels.md#s-6-4).

Der klassische Teil einer Abbildung zwischen Sektoren ist ihre **Übergangsmatrix**
$`T_{dc}(\rho)=\mathrm{Tr}\,\Phi_{dc}(\rho_c)/\mathrm{Tr}\,\rho_c`$ ($`\rho_c\ne0`$); die Methode
`transition_matrix()` der Bibliothek ist $`T=T(1)`$, also für den maximal gemischten Zustand jedes
Eingangssektors.

<a id="r-6-12"></a>

**Proposition 6.12.** Sei $`\Phi:\mathcal M\to\mathcal N`$ positiv und TP, $`\rho`$ ein Zustand mit
$`p_c=\mathrm{Tr}\,\rho_c>0`$ und $`q_d=\mathrm{Tr}\,\Phi(\rho)_d`$. Dann ist $`T(\rho)`$
spaltenstochastisch und $`q=T(\rho)p`$. Mit den bedingten Erwartungen $`E_{\mathcal M}`$, $`E_{\mathcal N}`$
auf die Zentren (Proposition 2.8) gilt
```math
E_{\mathcal N}\Phi E_{\mathcal M}(\rho)=\bigoplus_d\frac{(Tp)_d}{k'_d}\,1_{k'_d},\qquad T=T(1).
```

*Beweis.* $`T_{dc}\ge0`$, da $`\Phi_{dc}`$ positiv ist, und
$`\sum_dT_{dc}=\mathrm{Tr}\,\Phi(\iota_c\rho_c)/\mathrm{Tr}\,\rho_c=1`$. Ferner ist
$`q_d=\sum_c\mathrm{Tr}\,\Phi_{dc}(\rho_c)=\sum_cT_{dc}(\rho)p_c`$. Schließlich ist
$`E_{\mathcal M}\rho=\bigoplus_cp_c1_{k_c}/k_c`$, also hat
$`\Phi(E_{\mathcal M}\rho)_d=\sum_cp_c\Phi_{dc}(1_{k_c}/k_c)`$ die Spur $`\sum_cT_{dc}p_c`$, und
$`E_{\mathcal N}`$ ersetzt den Block $`d`$ durch seine Spur mal $`1/k'_d`$. $`\square`$

*In der Bibliothek:* Die Übergangsmatrix ist `InterSectorChannel.transition_matrix()`; eine Überprüfung von Proposition 6.12 — [Tutorial, Abschnitt 6.4](../../tutorial/de/06_quantum_channels.md#s-6-4).

**Kohärenzen eines Qubits.** Auf $`M_2`$ sei $`D_p`$ die Dephasierung und $`A_\gamma`$ die
Amplitudendämpfung. Nach der Tabelle und der expliziten Form oben multipliziert $`D_p`$ die Einträge
$`\rho_{01}`$ und $`\rho_{10}`$ mit $`1-p`$ und lässt die Diagonale fest, während $`A_\gamma`$ die Einträge
$`\rho_{01}`$ und $`\rho_{10}`$ mit $`\sqrt{1-\gamma}`$ multipliziert und die Diagonale auf
$`(\rho_{00}+\gamma\rho_{11},(1-\gamma)\rho_{11})`$ abbildet, unabhängig von den Kohärenzen. Daher ist
$`D_pA_\gamma=A_\gamma D_p`$: Beide multiplizieren die Kohärenzen mit $`(1-p)\sqrt{1-\gamma}`$ und wirken auf
der Diagonale wie $`A_\gamma`$. Insbesondere multipliziert die Hintereinanderausführung für
$`T_1,T_\varphi>0`$, $`t\ge0`$, $`1-p=e^{-t/T_\varphi}`$ und $`1-\gamma=e^{-t/T_1}`$ die Kohärenzen mit
$`e^{-t/T_2}`$, wobei
```math
\frac1{T_2}=\frac1{2T_1}+\frac1{T_\varphi},
```
und diese Familien setzen sich additiv in $`t`$ zusammen: $`D_{p(t)}D_{p(s)}=D_{p(t+s)}`$ und
$`A_{\gamma(t)}A_{\gamma(s)}=A_{\gamma(t+s)}`$, da sich die Faktoren $`e^{-t/T}`$ multiplizieren und die Spur
erhalten bleibt.

*Rechnung:* Der Kohärenzfaktor und $`1/T_2=1/(2T_1)+1/T_\varphi`$ — [Tutorial, Abschnitt 6.4](../../tutorial/de/06_quantum_channels.md#s-6-4).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Für ein Qubit
> modelliert die Dephasierung die reine Dephasierung mit der Zeit $`T_\varphi`$ und die
> Amplitudendämpfung den Zerfall des angeregten Niveaus mit der Zeit $`T_1`$; die obige Rechnung liefert
> dann die Kohärenzzeit $`T_2`$ eines Qubits, das beidem ausgesetzt ist. Satz 6.10
> ist die algebraische Form dekohärenzfreier Information: Was mit allen Sprungoperatoren kommutiert,
> überlebt (Kapitel 10). Proposition 6.12 besagt, dass eine Abbildung zwischen Sektoren, durch das
> Zentrum betrachtet, eine klassische Markow-Matrix auf den Ladungen ist.

<a id="s-6-5"></a>

## 6.5 Zufällige Kanäle

Eine **Haar-zufällige Isometrie** $`V:\mathbb C^k\to\mathbb C^N`$ ist, wie in der Bibliothek
(`random_isometry`), der Faktor $`Q`$ von $`G=QR`$, wobei $`G`$ eine $`N\times k`$-Ginibre-Matrix ist
(unabhängige Einträge mit der Dichte $`\pi^{-1}e^{-|z|^2}`$), $`Q^*Q=1`$ und $`R`$ eine obere
Dreiecksmatrix mit positiver Diagonale (auf der CPU nimmt die Bibliothek eine LAPACK-Zerlegung
$`G=Q'R'`$ und multipliziert die $`i`$-te Spalte von $`Q'`$ mit $`R'_{ii}/|R'_{ii}|`$; das ergibt diese
Zerlegung, wie in Satz 5.3, sobald alle $`R'_{ii}\ne0`$ sind, d. h. sobald $`G`$ Rang $`k`$ hat, was nach
dem Beweis von Lemma 6.13 fast sicher gilt; für große Stapel auf einer GPU verwendet sie
Gram–Schmidt, das sie direkt liefert). `random_channel(alg, r)` zieht solche
$`V_c:\mathbb C^{k_c}\to\mathbb C^r\otimes\mathbb C^{k_c}`$ unabhängig in jedem Block und verwendet die
Kraus-Operatoren aus Satz 6.5.

*In der Bibliothek:* `random_isometry` und `channels.random_channel` — [Tutorial, Abschnitt 6.5](../../tutorial/de/06_quantum_channels.md#s-6-5).

<a id="r-6-13"></a>

**Lemma 6.13.** $`V`$ hat die Verteilung der ersten $`k`$ Spalten einer Haar-verteilten unitären Matrix
aus $`U(N)`$. Insbesondere hat $`UV`$ für jedes unitäre $`U\in M_N`$ dieselbe Verteilung wie $`V`$.

*Beweis.* Sei $`Z`$ eine $`N\times N`$-Ginibre-Matrix mit den ersten $`k`$ Spalten $`G`$ und
$`E\in M_{N\times k}`$ die ersten $`k`$ Spalten von $`1_N`$. Nach Satz 5.3 ist $`Z`$ fast sicher invertierbar
und $`Z=Q_ZR_Z`$ mit Haar-verteiltem $`Q_Z`$. Da $`R_Z`$ eine obere Dreiecksmatrix ist, ist
$`G=ZE=(Q_ZE)(E^*R_ZE)`$ eine Faktorisierung der verlangten Art, und zwar die einzige: Ist
$`G=QR=Q'R'`$, so gilt $`R^*R=G^*G=R'^*R'`$, also ist $`R'R^{-1}`$ unitär und eine obere Dreiecksmatrix mit
positiver Diagonale, folglich gleich $`1`$ (Beweis von Lemma 5.2). Also ist $`V=Q_ZE`$, und $`UV=(UQ_Z)E`$
hat die Verteilung von $`V`$, weil $`UQ_Z`$ die von $`Q_Z`$ hat. $`\square`$

<a id="r-6-14"></a>

**Lemma 6.14.** Sei $`u`$ ein zufälliger Einheitsvektor in $`\mathbb C^N`$, dessen Verteilung unter allen
unitären Matrizen invariant ist. Dann ist $`\mathbb E[u_a\bar u_b]=\delta_{ab}/N`$ und
```math
\mathbb E\big[u_a\bar u_b\,u_c\bar u_d\big]=\frac{\delta_{ab}\delta_{cd}+\delta_{ad}\delta_{cb}}{N(N+1)} .
```

*Beweis.* Die Invarianz unter $`\mathrm{diag}(e^{i\theta_1},\dots,e^{i\theta_N})`$ multipliziert
die Momente mit $`e^{i(\theta_a-\theta_b)}`$ bzw. $`e^{i(\theta_a-\theta_b+\theta_c-\theta_d)}`$; sie
verschwinden also, außer wenn $`a=b`$ bzw. $`\{a,c\}=\{b,d\}`$ als Multimengen. Wegen der Invarianz unter
Permutationsmatrizen ist $`\mathbb E|u_a|^2=1/N`$, und es bleiben nur $`\alpha=\mathbb E|u_1|^4`$ und
$`\beta=\mathbb E|u_1|^2|u_2|^2`$. Erstens liefert $`\sum_c|u_c|^2=1`$ die Gleichung
$`\alpha+(N-1)\beta=\mathbb E|u_1|^2=1/N`$. Zweitens ist $`(u_1+u_2)/\sqrt2`$ die erste Koordinate von
$`Uu`$ für ein unitäres $`U`$, also $`\mathbb E|u_1+u_2|^4=4\alpha`$. Mit $`s=|u_1|^2+|u_2|^2`$ und
$`z=u_1\bar u_2`$ ist $`|u_1+u_2|^4=s^2+2s(z+\bar z)+(z+\bar z)^2`$; die Invarianzen $`u_1\mapsto-u_1`$
und $`u_1\mapsto iu_1`$ liefern $`\mathbb E[s(z+\bar z)]=0`$ und $`\mathbb Ez^2=0`$, also
$`4\alpha=(2\alpha+2\beta)+2\beta`$, d. h. $`\alpha=2\beta`$. Daher ist $`\beta=1/(N(N+1))`$ und
$`\alpha=2\beta`$, und das ist die Formel in allen Fällen. (Hier wurde $`N\ge2`$ benutzt; für $`N=1`$ ist
$`|u_1|=1`$, und beide Seiten sind gleich $`1`$.) $`\square`$

<a id="r-6-15"></a>

**Satz 6.15 (Ausgangsreinheit eines zufälligen Kanals).** Sei $`\psi`$ ein reiner Zustand in einem
Sektor der Größe $`k`$ und $`\Phi`$ von `random_channel` mit Kraus-Rang $`r`$ gezogen. Dann gilt
```math
\mathbb E\mathrm{Tr}\,\Phi(\psi)^2=\frac{k+r}{kr+1},\qquad\mathbb E\,\Phi(\psi)=\frac1k1_k .
```

*Beweis.* Im Sektor von $`\psi`$ ist $`\Phi(\psi)=\mathrm{Tr}_E|u\rangle\langle u|`$ mit $`u=V\psi`$,
einem Einheitsvektor, dessen Verteilung nach Lemma 6.13 unitär invariant ist. Man indiziere
$`\mathbb C^r\otimes\mathbb C^k`$ durch Paare $`(e,x)`$. Dann ist $`\Phi(\psi)_{xy}=\sum_eu_{ex}\bar u_{ey}`$
mit dem Mittelwert $`r\delta_{xy}/(kr)`$, und
$`\mathrm{Tr}\,\Phi(\psi)^2=\sum_{x,y,e,f}u_{ex}\bar u_{ey}u_{fy}\bar u_{fx}`$. Nach Lemma 6.14 mit
$`N=kr`$ erzwingt der erste Term $`x=y`$ ($`kr^2`$ Indextupel) und der zweite $`e=f`$ ($`k^2r`$ Tupel), also
ist der Mittelwert $`(kr^2+k^2r)/(kr(kr+1))`$. $`\square`$

*Rechnung:* Satz 6.15 per Monte-Carlo für $`r=1,2,4`$ — [Tutorial, Abschnitt 6.5](../../tutorial/de/06_quantum_channels.md#s-6-5).

Für $`r=1`$ ist der Kanal unitär und die Reinheit gleich $`1`$; für $`r\to\infty`$ strebt sie gegen $`1/k`$,
die Reinheit des maximal gemischten Zustands. Der Beweis hat nur die unitäre Invarianz der
Verteilung von $`u`$ benutzt; er liefert daher auch Lubkins mittlere Reinheit des reduzierten Zustands
eines zufälligen reinen Zustands $`u=g/\|g\|`$ mit Gaußschem $`g`$ (dessen Verteilung nach
Satz A.20(b) unitär invariant ist), die Kapitel 7, Satz 7.11, mit Gaußschen Integralen beweist; die
Zahlen $`\alpha,\beta`$ aus Lemma 6.14 sind die aus Korollar 5.11 für die Einträge einer
Haar-verteilten unitären Matrix, und Korollar A.21 berechnet sie für $`u=g/\|g\|`$ erneut.

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Ein zufälliger Kanal modelliert generisches, unstrukturiertes Rauschen
> aus einer Umgebung der Dimension $`r`$. Schon $`r=4`$ drückt die Reinheit eines reinen Qutrit-Zustands
> im Mittel auf $`7/13\approx0.54`$ (Satz 6.15), verglichen mit $`1`$ ohne Rauschen und dem Minimum $`1/3`$
> (nach Cauchy–Schwarz ist $`\mathrm{Tr}\,\rho^2\ge(\mathrm{Tr}\,\rho)^2/k`$). In einem
> System mit Superauswahl wirkt das Rauschen in jedem Sektor getrennt, und maßgeblich ist die
> Dimension des Sektors, nicht die des ganzen Raumes.

<a id="s-6-6"></a>

## 6.6 Womit man vergleichen sollte

*In der Bibliothek:* Vergleich mit anderer Software — [Tutorial, Abschnitt 6.6](../../tutorial/de/06_quantum_channels.md#s-6-6).

<a id="x-6"></a>

## Aufgaben

1. (Werner–Holevo) Man zeige, dass $`\Phi(\rho)=\big(\mathrm{Tr}(\rho)1-\rho^T\big)/(k-1)`$ auf
   $`M_k`$, $`k\ge2`$, ein unitaler Kanal mit Kraus-Rang $`k(k-1)/2`$ ist, obwohl er aus der Transposition
   gebildet ist.
2. Man zeige, dass der Dephasierungskanal auf $`M_k`$ für $`0<p\le1`$ den Kraus-Rang $`k`$ hat, einen
   weniger als die Anzahl der Kraus-Operatoren, die die Bibliothek speichert.
3. Sei $`X`$ die zyklische Verschiebung auf $`\mathbb C^k`$ und $`Z=\mathrm{diag}(1,\omega,\dots,\omega^{k-1})`$,
   $`\omega=e^{2\pi i/k}`$. Man zeige $`\frac1{k^2}\sum_{a,b=0}^{k-1}X^aZ^b\rho Z^{-b}X^{-a}=\mathrm{Tr}(\rho)1/k`$:
   Die bedingte Erwartung auf das Zentrum von $`M_k`$ ist ein gemischt-unitärer Kanal.
4. Man zeige $`(\Psi\circ\Phi)^\sharp=\Phi^\sharp\circ\Psi^\sharp`$, wenn auf der mittleren Algebra
   dieselbe Spur verwendet wird.

*Rechnung:* Die vier Aufgaben numerisch — [Tutorial, Aufgaben](../../tutorial/de/06_quantum_channels.md#x-6).

<details>
<summary>Lösungen</summary>

1. Mit dem Swap-Operator $`F`$ ist $`J(\theta)=\sum_{ij}E_{ij}\otimes E_{ji}=F`$ und $`J(\rho\mapsto\mathrm{Tr}(\rho)1)=1\otimes1`$,
   also $`J(\Phi)=(1-F)/(k-1)=2P_-/(k-1)`$, wobei $`P_-=(1-F)/2`$ die Projektion auf die
   antisymmetrischen Vektoren ist ($`F=F^*`$, $`F^2=1`$), vom Rang $`\mathrm{Tr}\,P_-=(k^2-\mathrm{Tr}\,F)/2`$
   mit $`\mathrm{Tr}\,F=\sum_{ij}\langle e_i\otimes e_j,e_j\otimes e_i\rangle=k`$. Also ist $`J\ge0`$ mit Rang $`k(k-1)/2`$ (Satz 6.2,
   Korollar 6.3), und $`\mathrm{Tr}\,\Phi(\rho)=\mathrm{Tr}\,\rho`$, $`\Phi(1)=1`$ sieht man direkt.

2. $`J=(1-p)|\Omega\rangle\langle\Omega|+p\sum_iE_{ii}\otimes E_{ii}`$ verschwindet außerhalb des
   Aufspanns der $`e_i\otimes e_i`$, und auf diesem $`k`$-dimensionalen Raum ist es die Matrix
   $`(1-p)\mathbf 1\mathbf 1^T+p1`$ ($`\mathbf 1`$ der Vektor aus Einsen) mit den Eigenwerten $`p`$ und
   $`(1-p)k+p`$, die alle positiv sind. Also ist $`\mathrm{rank}\,J=k`$, während die Bibliothek
   $`1+k_{\max}`$ Operatoren speichert.

3. $`\mathrm{Tr}((X^aZ^b)^*X^{a'}Z^{b'})=k\,\delta_{aa'}\delta_{bb'}`$: $`Z^{-b}X^{a'-a}Z^{b'}`$ hat
   verschwindende Diagonale, außer wenn $`a=a'`$, und $`\mathrm{Tr}\,Z^{b'-b}=\sum_j\omega^{j(b'-b)}=k\delta_{bb'}`$.
   Also bilden die $`B_{ab}=X^aZ^b/\sqrt k`$ eine Orthonormalbasis von $`M_k`$ bezüglich
   $`\mathrm{Tr}(A^*B)`$, und je zwei Orthonormalbasen gehen durch eine unitäre Matrix
   auseinander hervor; nach Satz 6.7 ist $`\sum_mB_m\rho B_m^*`$ für alle Orthonormalbasen gleich. Für
   die Basis $`E_{ab}`$ ist es $`\sum_{ab}E_{ab}\rho E_{ba}=\mathrm{Tr}(\rho)1`$. Man teile durch $`k`$.

4. $`\mathrm{tr}(\Psi(\Phi(x))z)=\mathrm{tr}(\Phi(x)\Psi^\sharp(z))=\mathrm{tr}(x\,\Phi^\sharp(\Psi^\sharp(z)))`$
   für alle $`x`$ und $`z`$; die Eindeutigkeit der dualen Abbildung liefert die Behauptung.

</details>

**Weiter:** [Kapitel 7](07_tensor_products.md) behandelt Tensorprodukte und Verschränkung: die partielle
Spur – einen Kanal zwischen verschiedenen Algebren, dessen duale Abbildung die Einbettung
$`A\mapsto A\otimes1`$ ist –, die Page-Formel und die Fusion von Ladungen.
