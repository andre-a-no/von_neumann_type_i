[English](../D_gksl.md) · [Español](../es/D_gksl.md) · [Français](../fr/D_gksl.md) · **Deutsch** · [中文](../zh/D_gksl.md) · [日本語](../ja/D_gksl.md) · [Русский](../ru/D_gksl.md)

# Anhang D. Generatoren quantendynamischer Halbgruppen

Kapitel 8 hat bewiesen, dass ein Generator in Gorini–Kossakowski–Sudarshan–Lindblad-Form (GKSL-Form) eine
quantendynamische Halbgruppe erzeugt (Satz 8.6), und hat die Umkehrung (Satz 8.7) sowie ihre Version für
$`\mathcal M=\bigoplus_cM_{k_c}`$ verwendet. Dieser Anhang beweist beides in endlicher Dimension, zusammen mit der
Freiheit in der GKSL-Form und einem Kriterium für die Eindeutigkeit des stationären Zustands. Das Werkzeug ist die
Choi-Matrix aus Kapitel 6: Ein Generator ist im Allgemeinen nicht vollständig positiv, aber seine Choi-Matrix ist positiv auf
dem orthogonalen Komplement des maximal verschränkten Vektors. Dieselbe Art von Argument – Positivität
einer Familie von Operatoren, dann eine Zerlegung – beweist im letzten Abschnitt den in Kapitel 10 formulierten Satz
über Fehlerkorrektur mit Operatoralgebren (Satz 10.7), mit einem expliziten Wiederherstellungskanal.

**Sie lernen**, warum der Generator einer vollständig positiven Halbgruppe bedingt vollständig
positiv ist, wie man $`H`$ und die Sprungoperatoren an der Choi-Matrix eines Generators abliest, in welchem Sinn
sie eindeutig sind, warum jeder Generator auf einer Algebra mit Sektoren die in Abschnitt 8.3 verwendete Form hat, wann
der stationäre Zustand eindeutig ist und welche Observablenalgebren ein Wiederherstellungskanal schützen kann.

Im Folgenden ist $`M_n=M_n(\mathbb C)`$, die $`E_{ij}=|i\rangle\langle j|`$ sind seine Matrixeinheiten, ein
**Superoperator** ist eine lineare Abbildung auf $`M_n`$ (oder auf $`\mathcal M`$), und
```math
\Omega=\sum_{i=1}^ne_i\otimes e_i\in\mathbb C^n\otimes\mathbb C^n,\qquad
J(\Psi)=\sum_{i,j=1}^nE_{ij}\otimes\Psi(E_{ij})=[\Psi(E_{ij})]_{i,j=1}^n
```
sind der maximal verschränkte Vektor und die **Choi-Matrix** eines Superoperators $`\Psi`$ (Satz 6.2;
hier muss $`\Psi`$ nicht positiv sein). Für Vektoren $`v,w`$ ist $`|v\rangle\langle w|`$ der Operator
$`x\mapsto\langle w,x\rangle v`$. Ein Superoperator $`\Psi`$ heißt **hermitizitätserhaltend**, wenn
$`\Psi(A^*)=\Psi(A)^*`$ für alle $`A`$ gilt. Quantendynamische Halbgruppen und ihre Generatoren sind die aus
Abschnitt 8.2 (Lemma 8.5): $`T_t=e^{t\mathcal L}`$ mit Kanälen $`T_t`$.

*Rechnung:* Die numerischen Prüfungen dieses Anhangs — [Tutorial, Anhang D](../../tutorial/de/D_gksl.md).

<a id="s-D-1"></a>

## D.1 Choi-Matrizen der Bausteine

<a id="r-D-1"></a>

**Lemma D.1.** Sei $`\Psi`$ ein Superoperator auf $`M_n`$.

(a) $`\Psi\mapsto J(\Psi)`$ ist eine lineare Bijektion von den Superoperatoren auf $`M_n`$ auf $`M_n(M_n)=M_{n^2}`$.
(b) Jedes $`v\in\mathbb C^n\otimes\mathbb C^n`$ ist $`v=(1\otimes K)\Omega`$ für genau ein $`K\in M_n`$, und
$`\langle(1\otimes K)\Omega,(1\otimes L)\Omega\rangle=\mathrm{Tr}(K^*L)`$. Insbesondere ist
$`\langle\Omega,(1\otimes K)\Omega\rangle=\mathrm{Tr}\,K`$.
(c) Für $`K,L\in M_n`$ hat der Superoperator $`A\mapsto KAL^*`$ die Choi-Matrix
$`|(1\otimes K)\Omega\rangle\langle(1\otimes L)\Omega|`$. Insbesondere ist $`J(\mathrm{id})=|\Omega\rangle\langle\Omega|`$.
(d) $`\Psi`$ ist genau dann hermitizitätserhaltend, wenn $`J(\Psi)`$ selbstadjungiert ist. Jede CP-Abbildung ist hermitizitätserhaltend.

*Beweis.* (a) Der $`(i,j)`$-Block von $`J(\Psi)`$ ist $`\Psi(E_{ij})`$, und die $`E_{ij}`$ bilden eine Basis von $`M_n`$; also ist
$`J`$ linear und injektiv, und jede Blockmatrix $`[X_{ij}]`$ ist die Choi-Matrix des Superoperators
mit $`E_{ij}\mapsto X_{ij}`$. (b) $`(1\otimes K)\Omega=\sum_ie_i\otimes Ke_i`$ hat als $`i`$-ten Block $`Ke_i`$, die
$`i`$-te Spalte von $`K`$, die $`K`$ bestimmt, und jedes $`v`$ entsteht auf diese Weise. Ferner ist
$`\langle\sum_ie_i\otimes Ke_i,\sum_je_j\otimes Le_j\rangle=\sum_i\langle Ke_i,Le_i\rangle=\mathrm{Tr}(K^*L)`$;
für die letzte Aussage ersetze man $`(K,L)`$ durch $`(1,K)`$. (c) Der $`(i,j)`$-Block der rechten Seite ist
$`Ke_i(Le_j)^*=KE_{ij}L^*`$. (d) Der $`(i,j)`$-Block von $`J(\Psi)^*`$ ist $`\Psi(E_{ji})^*`$. Ist $`\Psi`$
hermitizitätserhaltend, so ist $`\Psi(E_{ji})^*=\Psi(E_{ji}^*)=\Psi(E_{ij})`$. Gilt umgekehrt
$`\Psi(E_{ji})^*=\Psi(E_{ij})`$ für alle $`i,j`$ und ist $`A=\sum a_{ij}E_{ij}`$, so ist
$`\Psi(A)^*=\sum\bar a_{ij}\Psi(E_{ij})^*=\sum\bar a_{ij}\Psi(E_{ji})=\Psi(A^*)`$. Eine CP-Abbildung hat eine Kraus-Form
(Satz 6.2), und $`(\sum_lK_lAK_l^*)^*=\sum_lK_lA^*K_l^*`$. $`\square`$

<a id="r-D-2"></a>

**Lemma D.2 (Zerlegung).** Seien $`q_1,\dots,q_C`$ orthonormale Vektoren in
$`\mathbb C^n\otimes\mathbb C^n`$, geschrieben als $`q_c=(1\otimes Q_c)\Omega`$, sei $`W=\mathrm{span}\{q_c\}`$, und sei
$`\Psi`$ ein hermitizitätserhaltender Superoperator mit $`\langle v,J(\Psi)v\rangle\ge0`$ für alle $`v\perp W`$.
Man setze $`J=J(\Psi)`$ und $`\xi_c=Jq_c-\frac12\sum_d\langle q_d,Jq_c\rangle q_d`$, und $`K_c\in M_n`$ sei definiert durch
$`(1\otimes K_c)\Omega=\xi_c`$. Dann gibt es $`L_1,\dots,L_r\in M_n`$, $`r=\mathrm{rank}(P_W^\perp JP_W^\perp)`$
mit der Projektion $`P_W^\perp`$ auf $`W^\perp`$, so dass
```math
\Psi(A)=\sum_{j=1}^rL_jAL_j^*+\sum_{c=1}^C\big(K_cAQ_c^*+Q_cAK_c^*\big)\qquad(A\in M_n),
```
$`\mathrm{Tr}(Q_c^*L_j)=0`$ für alle $`c,j`$ und $`\mathrm{Tr}(L_j^*L_k)=0`$ für $`j\neq k`$.

*Beweis.* $`J`$ ist nach Lemma D.1(d) selbstadjungiert. Sei $`P=\sum_c|q_c\rangle\langle q_c|`$ die Projektion auf $`W`$
und $`P^\perp=P_W^\perp=1-P`$. Ausmultiplizieren von $`J=(P+P^\perp)J(P+P^\perp)`$ mit $`PJP^\perp=PJ-PJP`$ und
$`P^\perp JP=JP-PJP`$ ergibt
```math
J=P^\perp JP^\perp+\big(PJ+JP-PJP\big).
```
Der erste Term ist selbstadjungiert, und $`\langle v,P^\perp JP^\perp v\rangle=\langle P^\perp v,JP^\perp v\rangle\ge0`$,
weil $`P^\perp v\perp W`$; also ist $`P^\perp JP^\perp\ge0`$. Nach Satz 1.3 ist
$`P^\perp JP^\perp=\sum_{j=1}^r|l_j\rangle\langle l_j|`$, wobei die $`l_j`$ die mit $`\sqrt{\mu_j}`$ skalierten Eigenvektoren zu den
von null verschiedenen Eigenwerten $`\mu_j`$ sind: Sie sind paarweise orthogonal und liegen im Bild von
$`P^\perp JP^\perp`$, das orthogonal zu $`W`$ ist. Man schreibe $`l_j=(1\otimes L_j)\Omega`$; nach Lemma D.1(b) ist
$`\mathrm{Tr}(Q_c^*L_j)=\langle q_c,l_j\rangle=0`$ und $`\mathrm{Tr}(L_j^*L_k)=\langle l_j,l_k\rangle=0`$ für $`j\ne k`$,
und nach Lemma D.1(c) ist $`P^\perp JP^\perp`$ die Choi-Matrix von $`A\mapsto\sum_jL_jAL_j^*`$. Für den zweiten Term gilt
$`JP=\sum_c|Jq_c\rangle\langle q_c|`$ und $`PJP=\sum_{c,d}\langle q_d,Jq_c\rangle|q_d\rangle\langle q_c|`$, also
```math
\sum_c|\xi_c\rangle\langle q_c|=JP-\tfrac12PJP,\qquad\sum_c|q_c\rangle\langle\xi_c|=PJ-\tfrac12PJP
```
(die zweite Gleichung ist die adjungierte der ersten; $`PJP`$ ist selbstadjungiert), und ihre Summe ist $`PJ+JP-PJP`$. Nach
Lemma D.1(c) sind $`|\xi_c\rangle\langle q_c|`$ und $`|q_c\rangle\langle\xi_c|`$ die Choi-Matrizen von
$`A\mapsto K_cAQ_c^*`$ und $`A\mapsto Q_cAK_c^*`$. Also haben beide Seiten der behaupteten Identität die Choi-Matrix $`J`$,
und Lemma D.1(a) liefert die Identität. $`\square`$

<a id="s-D-2"></a>

## D.2 Bedingte vollständige Positivität

Ein Superoperator $`\Psi`$ auf $`M_n`$ heißt **bedingt vollständig positiv** (CCP), wenn
$`\langle v,J(\Psi)v\rangle\ge0`$ für jedes $`v\perp\Omega`$ gilt. Nach Lemma D.1(b) sind diese $`v`$ die Vektoren
$`(1\otimes K)\Omega`$ mit $`\mathrm{Tr}\,K=0`$.

<a id="r-D-3"></a>

**Satz D.3.** Für einen Superoperator $`\mathcal L`$ auf $`M_n`$ sind äquivalent:

1. $`e^{t\mathcal L}`$ ist CP für jedes $`t\ge0`$;
2. $`\mathcal L`$ ist hermitizitätserhaltend und bedingt vollständig positiv;
3. es gibt eine CP-Abbildung $`\Phi`$ auf $`M_n`$ und $`G\in M_n`$ mit $`\mathcal L(A)=\Phi(A)+GA+AG^*`$ für alle $`A`$.

In 3 kann man $`\Phi(A)=\sum_{j=1}^rL_jAL_j^*`$ mit $`r\le n^2-1`$, $`\mathrm{Tr}\,L_j=0`$ und
$`\mathrm{Tr}(L_j^*L_k)=0`$ für $`j\neq k`$ sowie $`G`$ mit $`\mathrm{Tr}\,G\in\mathbb R`$ wählen.

*Beweis.* (1)$`\Rightarrow`$(2). Nach Lemma 8.1 gilt $`(e^{t\mathcal L}-\mathrm{id})/t\to\mathcal L`$ für $`t\to0`$, und $`J`$ ist
linear auf einem endlichdimensionalen Raum, also stetig; daher ist
$`J(e^{t\mathcal L})=J(\mathrm{id})+t\,J\big((e^{t\mathcal L}-\mathrm{id})/t\big)`$, wobei die letzte Choi-Matrix gegen
$`J(\mathcal L)`$ strebt. Für $`v\perp\Omega`$ ist $`\langle v,J(\mathrm{id})v\rangle=|\langle\Omega,v\rangle|^2=0`$ (Lemma D.1(c)), also
gilt für $`t>0`$, da $`J(e^{t\mathcal L})\ge0`$ (Satz 6.2),
```math
0\le\frac1t\langle v,J(e^{t\mathcal L})v\rangle=\big\langle v,J\big((e^{t\mathcal L}-\mathrm{id})/t\big)v\big\rangle
\;\xrightarrow[t\to0]{}\;\langle v,J(\mathcal L)v\rangle .
```
Jedes $`e^{t\mathcal L}`$ ist CP, also hermitizitätserhaltend (Lemma D.1(d)); daher ist
$`\mathcal L(A^*)=\lim_t(e^{t\mathcal L}(A)^*-A^*)/t=\mathcal L(A)^*`$.

(2)$`\Rightarrow`$(3). Man wende Lemma D.2 mit $`C=1`$ und $`q_1=\Omega/\sqrt n`$ an, d. h. $`Q_1=1/\sqrt n`$:
$`\mathcal L(A)=\sum_jL_jAL_j^*+GA+AG^*`$ mit $`G=K_1/\sqrt n`$, und $`\Phi=\sum_jL_j\cdot L_j^*`$ ist CP
(Satz 6.2). Ferner ist $`\mathrm{Tr}\,L_j=\sqrt n\mathrm{Tr}(Q_1^*L_j)=0`$, $`r\le\mathrm{rank}\,P_W^\perp=n^2-1`$,
und nach Lemma D.1(b) ist
$`\mathrm{Tr}\,G=\langle\Omega,\xi_1\rangle/\sqrt n=\langle q_1,\xi_1\rangle=\frac12\langle q_1,J(\mathcal L)q_1\rangle\in\mathbb R`$,
da $`J(\mathcal L)`$ selbstadjungiert ist.

(3)$`\Rightarrow`$(1). Die Schritte 1–4 des Beweises von Satz 8.6 verwenden nur, dass der Generator die Form
$`\mathcal L_0+\Phi`$ hat, mit $`\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)`$ für eine Matrix
$`H_{\mathrm{eff}}`$ und CP-Abbildung $`\Phi`$. Mit $`H_{\mathrm{eff}}=iG`$ ist $`\mathcal L_0(\rho)=G\rho+\rho G^*`$, und diese
Schritte zeigen, dass $`e^{t\mathcal L}=\sum_{m\ge0}S_m(t)`$ CP ist. $`\square`$

*Rechnung:* Satz D.3 an einem GKSL-Generator auf $`M_3`$ und das Beispiel $`\theta-\mathrm{id}`$ weiter unten — [Tutorial, Abschnitt D.2](../../tutorial/de/D_gksl.md#s-D-2). *In der Bibliothek:* Die Matrix eines GKSL-Generators, `dynamics.lindblad_superoperator` — [Tutorial, Abschnitt 8.2](../../tutorial/de/08_dynamics.md#s-8-2).

Positivität allein würde nicht genügen. Sei $`n\ge2`$, $`\theta`$ die Transposition auf $`M_n`$ und $`\mathcal L=\theta-\mathrm{id}`$. Wegen $`\theta^2=\mathrm{id}`$ ist
```math
e^{t\mathcal L}=e^{-t}\big(\cosh t\,\mathrm{id}+\sinh t\,\theta\big)=\tfrac{1+e^{-2t}}2\,\mathrm{id}+\tfrac{1-e^{-2t}}2\,\theta ,
```
eine Konvexkombination positiver spurerhaltender Abbildungen (Proposition 6.1). Aber $`J(\theta)=[E_{ji}]=\sum_{i,j}E_{ij}\otimes E_{ji}`$ ist der Swap-Operator $`e_k\otimes e_l\mapsto e_l\otimes e_k`$ und bildet daher
$`v=e_1\otimes e_2-e_2\otimes e_1`$, das orthogonal zu $`\Omega`$ ist, auf $`-v`$ ab (wie im Beweis von Proposition 6.1); mit
$`\|v\|^2=2`$ und $`\langle v,J(\mathrm{id})v\rangle=|\langle\Omega,v\rangle|^2=0`$ folgt
$`\langle v,J(\mathcal L)v\rangle=-2`$: $`\mathcal L`$ ist nicht CCP, und tatsächlich ist
$`\langle v,J(e^{t\mathcal L})v\rangle=-(1-e^{-2t})<0`$ für jedes $`t>0`$. Eine Halbgruppe positiver spurerhaltender Abbildungen
besteht nicht notwendig aus Kanälen.

<a id="s-D-3"></a>

## D.3 Der GKSL-Satz

<a id="r-D-4"></a>

**Satz D.4 (Gorini–Kossakowski–Sudarshan; Lindblad; Satz 8.7).** Ein Superoperator $`\mathcal L`$ auf
$`M_n`$ ist genau dann der Generator einer quantendynamischen Halbgruppe auf $`M_n`$, wenn es $`H=H^*`$ und
$`L_1,\dots,L_r`$ in $`M_n`$ gibt mit
```math
\mathcal L(\rho)=-i[H,\rho]+\sum_{j=1}^r\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
```
Man kann $`r\le n^2-1`$, $`\mathrm{Tr}\,H=0`$, $`\mathrm{Tr}\,L_j=0`$ und $`\mathrm{Tr}(L_j^*L_k)=0`$ für $`j\neq k`$ wählen.

*Beweis.* „Wenn“ ist Satz 8.6. Sei umgekehrt $`(T_t)`$ eine quantendynamische Halbgruppe mit Generator
$`\mathcal L`$, also $`T_t=e^{t\mathcal L}`$ (Lemma 8.5). Da jedes $`T_t`$ CP ist, liefert Satz D.3
$`\mathcal L(\rho)=\sum_jL_j\rho L_j^*+G\rho+\rho G^*`$ mit $`r\le n^2-1`$ spurfreien, paarweise orthogonalen $`L_j`$ und
$`\mathrm{Tr}\,G\in\mathbb R`$. Man setze $`H=\frac i2(G-G^*)`$ und $`R=\frac12(G+G^*)`$, beide selbstadjungiert; dann ist
$`G=-iH+R`$, $`G^*=iH+R`$ und $`G\rho+\rho G^*=-i[H,\rho]+\{R,\rho\}`$. Außerdem ist
$`\mathrm{Tr}\,H=\frac i2(\mathrm{Tr}\,G-\overline{\mathrm{Tr}\,G})=0`$. Da jedes $`T_t`$ TP ist, gilt
$`\mathrm{Tr}\,\mathcal L(\rho)=\frac d{dt}\mathrm{Tr}\,T_t(\rho)\big|_{t=0}=0`$ (Lemma 8.1), und wegen der Zyklizität
```math
0=\mathrm{Tr}\,\mathcal L(\rho)=\mathrm{Tr}\Big(\rho\,\Big(\sum_jL_j^*L_j+2R\Big)\Big)\qquad\text{for all }\rho .
```
Mit $`\rho=X=\sum_jL_j^*L_j+2R`$, das selbstadjungiert ist, folgt $`\mathrm{Tr}(X^*X)=0`$, also $`X=0`$ und
$`R=-\frac12\sum_jL_j^*L_j`$; das ist die GKSL-Form. $`\square`$

*Rechnung:* Der Beweis von Satz D.4, als Algorithmus ausgeführt, der $`H`$ und die Sprünge aus $`\mathcal L`$ zurückgewinnt — [Tutorial, Abschnitt D.3](../../tutorial/de/D_gksl.md#s-D-3).

<a id="r-D-5"></a>

**Satz D.5 (Freiheit in der GKSL-Form).** Seien $`(H;L_1,\dots,L_r)`$ und $`(H';L'_1,\dots,L'_s)`$ GKSL-Daten
($`H,H'`$ selbstadjungiert) mit Generatoren $`\mathcal L`$, $`\mathcal L'`$ wie in Satz D.4.

(a) $`\mathcal L'=\mathcal L`$ in jedem der folgenden Fälle:
(i) $`L'_j=L_j`$ und $`H'=H+c1`$, $`c\in\mathbb R`$;
(ii) $`H'=H`$ und, nach Auffüllen der kürzeren Familie mit Nullen, $`L'_j=\sum_ku_{jk}L_k`$ mit einer unitären Matrix $`u`$;
(iii) $`L'_j=L_j+a_j1`$ mit $`a_j\in\mathbb C`$ und $`H'=H-\frac i2\sum_j(\bar a_jL_j-a_jL_j^*)`$, das selbstadjungiert ist.

(b) Ist $`\mathrm{Tr}\,H=\mathrm{Tr}\,H'=0`$, sind alle $`L_j,L'_j`$ spurfrei und ist $`\mathcal L=\mathcal L'`$, so ist $`H=H'`$ und,
nach Auffüllen, $`L'_j=\sum_ku_{jk}L_k`$ mit einer unitären Matrix $`u`$.

(c) Folglich haben zwei GKSL-Daten genau dann denselben Generator, wenn die einen aus den anderen durch endlich viele
Schritte (i)–(iii) hervorgehen. Insbesondere sind der Operator $`H`$ der normierten Form aus Satz D.4 und ihr
dissipativer Anteil $`\Phi(\rho)=\sum_jL_j\rho L_j^*`$ durch $`\mathcal L`$ bestimmt.

*Beweis.* (a) (i) $`[c1,\rho]=0`$. (ii) Nach Satz 6.7 ist $`\Phi'=\sum_jL'_j\cdot L_j'^*`$ gleich
$`\Phi=\sum_jL_j\cdot L_j^*`$; die duale Abbildung hängt nicht von den Kraus-Operatoren ab (Abschnitt 6.2), also ist
$`\sum_jL_j'^*L'_j=\Phi^*(1)=\sum_jL_j^*L_j`$. (iii) Für einen einzelnen Sprung $`L`$ und $`a\in\mathbb C`$ ergibt Ausmultiplizieren
```math
(L+a)\rho(L+a)^*-\tfrac12\{(L+a)^*(L+a),\rho\}-\Big(L\rho L^*-\tfrac12\{L^*L,\rho\}\Big)=\tfrac12[B,\rho],
\qquad B=\bar aL-aL^* ,
```
denn die Terme $`\bar aL\rho+a\rho L^*`$ aus dem ersten Produkt verbinden sich mit
$`-\frac12(aL^*\rho+a\rho L^*+\bar aL\rho+\bar a\rho L)`$, und die Terme mit $`|a|^2`$ heben sich weg. Wegen $`B^*=-B`$ ist die
Matrix $`\frac i2B`$ selbstadjungiert, also ist $`H'=H-\frac i2\sum_jB_j`$ selbstadjungiert, und
$`-i[H'-H,\rho]=-\frac12\sum_j[B_j,\rho]`$ hebt die Änderung des Dissipators auf.

(b) Sei $`\Phi=\sum_jL_j\cdot L_j^*`$ und $`G=-iH-\frac12\Phi^*(1)`$, so dass $`\mathcal L(\rho)=\Phi(\rho)+G\rho+\rho G^*`$.
Nach Lemma D.1(c) ist
$`J(\mathcal L)=\sum_j|l_j\rangle\langle l_j|+|g\rangle\langle\Omega|+|\Omega\rangle\langle g|`$ mit $`l_j=(1\otimes L_j)\Omega`$ und
$`g=(1\otimes G)\Omega`$, und $`l_j\perp\Omega`$, weil $`\mathrm{Tr}\,L_j=0`$ (Lemma D.1(b)). Mit der Projektion $`P^\perp`$
auf $`\Omega^\perp`$ ist daher $`P^\perp J(\mathcal L)P^\perp=\sum_j|l_j\rangle\langle l_j|=J(\Phi)`$, und ebenso
$`P^\perp J(\mathcal L')P^\perp=J(\Phi')`$. Also ist $`J(\Phi)=J(\Phi')`$, $`\Phi=\Phi'`$ (Lemma D.1(a)), und Satz 6.7 liefert die
unitäre Matrix $`u`$. Dann ist $`\Phi^*(1)=\Phi'^*(1)`$, also $`[H-H',\rho]=0`$ für alle $`\rho`$. Eine Matrix $`X`$, die mit allen
$`E_{ij}`$ kommutiert, ist ein Vielfaches von $`1`$ (man vergleiche die Einträge von $`XE_{ij}`$ und $`E_{ij}X`$), also ist $`H-H'=\lambda1`$, und
$`\lambda=0`$ durch Spurbildung.

(c) Nach (a) ändern die Schritte (i)–(iii) den Generator nicht, und die Umkehrung jedes Schritts ist ein Schritt
derselben Art ((iii) mit $`-a_j`$, da $`\bar a_jL'_j-a_jL_j'^*=\bar a_jL_j-a_jL_j^*`$). Beliebige Daten werden durch (iii) mit
$`a_j=-\mathrm{Tr}\,L_j/n`$ und anschließend (i) mit $`c=-\mathrm{Tr}\,H''/n`$, wobei $`H''`$ der
Hamiltonoperator nach (iii) ist, in normierte Form gebracht; zwei normierte Daten mit demselben Generator hängen nach (b) über (ii) zusammen. $`\square`$

*Rechnung:* Satz D.5: die normierte Form durch (iii) und (i) sowie die Eindeutigkeit in (b) — [Tutorial, Abschnitt D.3](../../tutorial/de/D_gksl.md#s-D-3).

<a id="r-D-6"></a>

**Satz D.6 (Eindeutigkeit des stationären Zustands).** Sei $`\mathcal L`$ der GKSL-Generator aus Satz D.4
mit den Daten $`(H;L_1,\dots,L_r)`$. Hat $`\mathcal L`$ einen treuen stationären Zustand $`\sigma`$ ($`\sigma>0`$) und sind die einzigen
Matrizen, die mit $`H`$ und allen $`L_j`$, $`L_j^*`$ kommutieren, die Vielfachen von $`1`$, so ist $`\ker\mathcal L=\mathbb C\sigma`$:
$`\sigma`$ ist der einzige stationäre Zustand. Außerdem ist $`\ker\mathcal L^*=\mathbb C1`$.

*Beweis.* Die duale Abbildung bezüglich der Spurpaarung, $`\mathrm{Tr}(\mathcal L(\rho)A)=\mathrm{Tr}(\rho\,\mathcal L^*(A))`$, ist
$`\mathcal L^*(A)=i[H,A]+\sum_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})`$: Wegen der Zyklizität ist
$`\mathrm{Tr}(-i[H,\rho]A)=\mathrm{Tr}(\rho\,i[H,A])`$, $`\mathrm{Tr}(L\rho L^*A)=\mathrm{Tr}(\rho L^*AL)`$ und
$`\mathrm{Tr}(\{L^*L,\rho\}A)=\mathrm{Tr}(\rho\{L^*L,A\})`$ (das ist Aufgabe 1 aus Kapitel 8).
Direktes Ausmultiplizieren, bei dem sich die Hamiltonschen Terme wegheben, weil $`A\mapsto i[H,A]`$ eine Derivation ist, ergibt
für jedes $`X`$
```math
\mathcal L^*(X^*X)-X^*\mathcal L^*(X)-\mathcal L^*(X^*)X=\sum_j[L_j,X]^*[L_j,X]:
```
Für einen einzelnen Sprung heben sich von den sechs Antikommutatortermen vier paarweise weg, und die beiden anderen,
$`\frac12X^*L^*LX+\frac12X^*L^*LX`$, addieren sich, so dass die linke Seite $`L^*X^*XL-X^*L^*XL-L^*X^*LX+X^*L^*LX`$ ist, also
$`(XL-LX)^*(XL-LX)`$. Sei $`X\in\ker\mathcal L^*`$. Dann ist $`\mathcal L^*(X^*)=\mathcal L^*(X)^*=0`$, da $`\mathcal L^*`$
hermitizitätserhaltend ist, also $`\mathcal L^*(X^*X)=\sum_j[L_j,X]^*[L_j,X]=:Y\ge0`$, und
$`\mathrm{Tr}(\sigma Y)=\mathrm{Tr}(\mathcal L(\sigma)X^*X)=0`$. Wegen $`\sigma>0`$ verschwindet
$`\mathrm{Tr}(\sigma Y)=\mathrm{Tr}(\sigma^{1/2}Y\sigma^{1/2})`$ mit $`\sigma^{1/2}Y\sigma^{1/2}\ge0`$ nur, wenn
$`\sigma^{1/2}Y\sigma^{1/2}=0`$ ist, d. h. $`Y=0`$; Spurbildung ergibt $`\sum_j\mathrm{Tr}([L_j,X]^*[L_j,X])=0`$ mit nichtnegativen
Summanden, also $`[L_j,X]=0`$ für alle $`j`$. Dasselbe Argument für $`X^*\in\ker\mathcal L^*`$
liefert $`[L_j,X^*]=0`$, d. h. $`[L_j^*,X]=0`$. Dann kommutiert $`X`$ mit $`L_j^*L_j`$, der dissipative Anteil von $`\mathcal L^*(X)`$
verschwindet, und $`0=\mathcal L^*(X)=i[H,X]`$. Nach Voraussetzung ist $`X\in\mathbb C1`$; umgekehrt ist $`\mathcal L^*(1)=0`$. Also ist
$`\ker\mathcal L^*=\mathbb C1`$. In den Basen $`(E_{ij})`$ und $`(E_{ji})`$ von $`M_n`$, die bezüglich der Paarung
$`(\rho,A)\mapsto\mathrm{Tr}(\rho A)`$ dual zueinander sind, ist die Matrix von $`\mathcal L^*`$ die Transponierte derjenigen von $`\mathcal L`$; also haben beide
denselben Rang, und $`\dim\ker\mathcal L=\dim\ker\mathcal L^*=1`$. Ein Zustand in $`\mathbb C\sigma`$ mit Spur $`1`$ ist $`\sigma`$. $`\square`$

Dies ist der endlichdimensionale Fall eines Kriteriums von Frigerio [*Commun. Math. Phys.* 63, 269 (1978)],
das in Abschnitt 8.5 verwendet wird. Beide Voraussetzungen sind nötig. In Aufgabe 4 aus Kapitel 8 ist der Zustand $`\tfrac12\,1`$ treu
und stationär, aber $`H`$ und $`\sigma_z`$ kommutieren mit allen Diagonalmatrizen, und die stationären Zustände spannen den
zweidimensionalen Raum der Diagonalmatrizen auf. Umgekehrt seien auf $`M_3`$ $`H=0`$, $`L_1=|e_1\rangle\langle e_3|`$,
$`L_2=|e_2\rangle\langle e_3|`$. Die von $`L_1,L_2`$ erzeugte $`*`$-Algebra enthält $`L_1L_2^*=E_{12}`$, $`L_1^*L_1=E_{33}`$ und damit
alle Matrixeinheiten, also kommutieren nur die Vielfachen von $`1`$ mit $`L_j,L_j^*`$; aber jedes $`\rho=P\rho P`$ mit $`P=E_{11}+E_{22}`$
ist stationär, da $`L_j\rho=0=\rho L_j^*`$ und $`\{E_{33},\rho\}=0`$, also ist $`\dim\ker\mathcal L\ge4`$; nach Satz D.6 ist kein stationärer
Zustand treu.

*Rechnung:* Satz D.6 für den Generator aus Abschnitt D.2 und das zweite Beispiel — [Tutorial, Abschnitt D.3](../../tutorial/de/D_gksl.md#s-D-3).

<a id="s-D-4"></a>

## D.4 Algebren mit Sektoren

Sei $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}\subseteq M_D`$ die Algebra der blockdiagonalen Matrizen, $`D=\sum_ck_c`$, mit
minimalen zentralen Projektionen $`P_c`$, und sei $`E(X)=\sum_cP_cXP_c`$ die bedingte Erwartung von $`M_D`$ auf
$`\mathcal M`$. Wie im Beweis von Proposition 8.8 gilt $`E(B)=B`$, $`E(XB)=E(X)B`$ und $`E(BX)=BE(X)`$ für $`B\in\mathcal M`$
sowie $`E(X^*)=E(X)^*`$. Für eine CP-Abbildung $`\Phi:\mathcal M\to\mathcal M`$ ist $`\Phi^*`$ ihre ungewichtete duale Abbildung (Proposition 6.8).

<a id="r-D-7"></a>

**Satz D.7.** Sei $`\mathcal L`$ der Generator einer quantendynamischen Halbgruppe $`(T_t)`$ auf $`\mathcal M`$. Dann gibt
es eine CP-Abbildung $`\Phi:\mathcal M\to\mathcal M`$ und $`H=H^*\in\mathcal M`$ mit
```math
\mathcal L(\rho)=\Phi(\rho)-\tfrac12\{\Phi^*(1),\rho\}-i[H,\rho]\qquad(\rho\in\mathcal M).
```
Nach Proposition 6.4 ist $`\Phi(\rho)_d=\sum_{c,l}K^{dc}_l\rho_cK^{dc*}_l`$ mit Kraus-Operatoren zwischen den Sektoren;
umgekehrt erzeugt jedes solche $`\mathcal L`$ eine quantendynamische Halbgruppe auf $`\mathcal M`$ (Proposition 8.8(1)).
Also sind die Generatoren aus Abschnitt 8.3 (mit $`\Phi=\sum_jg_jL_j\cdot L_j^*+\sum_mg_m\Phi_m`$) genau die
Generatoren quantendynamischer Halbgruppen auf $`\mathcal M`$.

*Beweis.* *Schritt 1: eine CP-Abbildung auf $`M_D`$.* Nach Proposition 6.4 ist
$`T_t(\rho)_d=\sum_{c,l}K^{dc}_l(t)\rho_cK^{dc}_l(t)^*`$. Sei $`N^{dc}_l(t)\in M_D`$ die Matrix mit dem Block $`K^{dc}_l(t)`$ an der
Position $`(d,c)`$ und Nullen sonst. Für $`X\in M_D`$ ist $`N^{dc}_lXN^{dc*}_l`$ der Block
$`K^{dc}_lX_{cc}K^{dc*}_l`$ an der Position $`(d,d)`$, wobei $`X_{cc}`$ der $`c`$-te Diagonalblock ist, d. h. der $`c`$-te Block
von $`E(X)`$. Daher ist $`\hat T_t=T_t\circ E=\sum_{d,c,l}N^{dc}_l\cdot N^{dc*}_l`$ CP auf $`M_D`$ (Satz 6.2), mit
$`\hat T_0=E`$. Nach Lemma 8.5 gilt $`(\hat T_t-E)/t=\big((T_t-\mathrm{id})/t\big)\circ E\to\mathcal L\circ E`$.

*Schritt 2: Positivität abseits der Sektorvektoren.* Sei $`p_c=(1\otimes P_c)\Omega\in\mathbb C^D\otimes\mathbb C^D`$ und
$`W=\mathrm{span}\{p_c\}`$. Nach Lemma D.1(c) ist $`J(E)=\sum_c|p_c\rangle\langle p_c|`$, also gilt für $`v\perp W`$ und $`t>0`$
$`0\le\frac1t\langle v,J(\hat T_t)v\rangle=\langle v,J((\hat T_t-E)/t)v\rangle\to\langle v,J(\mathcal L\circ E)v\rangle`$.
Jedes $`T_t`$ ist hermitizitätserhaltend (wegen seiner Kraus-Form), also auch $`\mathcal L`$ (wie in Satz D.3), und ebenso $`E`$
und $`\mathcal L\circ E`$.

*Schritt 3: Zerlegung.* Nach Lemma D.1(b) ist $`\langle p_c,p_d\rangle=\mathrm{Tr}(P_cP_d)=\delta_{cd}k_c`$, also sind die
$`q_c=p_c/\sqrt{k_c}=(1\otimes Q_c)\Omega`$ mit $`Q_c=P_c/\sqrt{k_c}`$ orthonormal. Lemma D.2 für
$`\Psi=\mathcal L\circ E`$ liefert eine CP-Abbildung $`\hat\Phi=\sum_j\hat L_j\cdot\hat L_j^*`$ auf $`M_D`$ und $`K_c\in M_D`$ mit
```math
\mathcal L(E(A))=\hat\Phi(A)+\sum_c\big(K_cAP_c+P_cAK_c^*\big)/\sqrt{k_c}\qquad(A\in M_D).
```
Für $`A\in\mathcal M`$ ist $`E(A)=A`$ und $`AP_c=P_cA`$, also ist die Summe gleich $`\tilde GA+A\tilde G^*`$ mit
$`\tilde G=\sum_cK_cP_c/\sqrt{k_c}`$.

*Schritt 4: Kompression.* Für $`A\in\mathcal M`$ ist $`\mathcal L(A)\in\mathcal M`$, also
$`\mathcal L(A)=E(\mathcal L(A))=E(\hat\Phi(A))+GA+AG^*`$ mit $`G=E(\tilde G)\in\mathcal M`$. Die Abbildung
$`\Phi=E\circ\hat\Phi`$ auf $`\mathcal M`$ ist
$`\Phi(A)_d=\sum_{c,j}(P_d\hat L_jP_c)A_c(P_d\hat L_jP_c)^*`$: Sie ist CP nach Proposition 6.4, mit den Blöcken der $`\hat L_j`$ als
Kraus-Operatoren.

*Schritt 5: Spurerhaltung.* Wie in Satz D.4 ist $`\mathrm{Tr}\,\mathcal L(A)=0`$ für $`A\in\mathcal M`$, und
$`\mathrm{Tr}\,\Phi(A)=\mathrm{Tr}(A\Phi^*(1))`$, also $`\mathrm{Tr}(AX)=0`$ für alle $`A\in\mathcal M`$, wobei
$`X=\Phi^*(1)+G+G^*\in\mathcal M`$ selbstadjungiert ist. Mit $`A=X`$ folgt $`X=0`$ (die ungewichtete Spur ist treu, Satz 3.1).
Mit $`H=\frac i2(G-G^*)\in\mathcal M`$ ist $`G=-iH-\frac12\Phi^*(1)`$ und $`G\rho+\rho G^*=-i[H,\rho]-\frac12\{\Phi^*(1),\rho\}`$. $`\square`$

Die Vektoren $`p_c`$ treten an die Stelle von $`\Omega=\sum_cp_c`$, weil $`\hat T_0=E`$ ist und nicht $`\mathrm{id}`$: $`J(\mathcal L\circ E)`$ muss
auf $`\Omega^\perp`$ nicht positiv sein, und im folgenden Beispiel ist sie es nicht.

*Rechnung:* Satz D.7 auf $`M_1\oplus M_2\subseteq M_3`$, entlang der fünf Schritte des Beweises; dort ist $`J(\mathcal L\circ E)`$ auf $`\Omega^\perp`$ nicht positiv — [Tutorial, Abschnitt D.4](../../tutorial/de/D_gksl.md#s-D-4). *In der Bibliothek:* Generatoren auf einer Algebra mit Sektoren, `dynamics.lindblad_rhs` mit `channels.random_inter_sector_channel` — [Tutorial, Abschnitt 8.3](../../tutorial/de/08_dynamics.md#s-8-3).

<a id="r-D-8"></a>

**Korollar D.8 (Heisenberg-Bild).** Sei $`(S_t)_{t\ge0}`$ eine stetige Familie unitaler CP-Abbildungen auf
$`\mathcal M`$ mit $`S_0=\mathrm{id}`$ und $`S_{s+t}=S_sS_t`$. Ihr Generator ist $`A\mapsto\Psi(A)+K^*A+AK`$ mit
einer CP-Abbildung $`\Psi`$ auf $`\mathcal M`$, $`K\in\mathcal M`$ und $`\Psi(1)+K+K^*=0`$.

*Beweis.* Sei $`T_t=S_t^*`$ (ungewichtete duale Abbildungen). Sie sind CP (Proposition 6.8) und TP (Korollar 6.9(b)),
$`T_0=\mathrm{id}`$, und $`T_{s+t}=(S_tS_s)^*=S_s^*S_t^*=T_sT_t`$, weil
$`\mathrm{Tr}(S_tS_s(x)y)=\mathrm{Tr}(xS_s^*S_t^*(y))`$ und $`S_tS_s=S_{s+t}`$. Da $`\Psi\mapsto\Psi^*`$ linear auf einem
endlichdimensionalen Raum ist, ist $`t\mapsto T_t`$ stetig und $`(T_t-\mathrm{id})/t\to\mathcal L_S^*`$, wobei $`\mathcal L_S`$ der
Generator von $`(S_t)`$ ist. Satz D.7 liefert $`\mathcal L_S^*(\rho)=\Phi(\rho)+G\rho+\rho G^*`$ mit $`G+G^*=-\Phi^*(1)`$, und wegen der
Zyklizität ist $`\mathrm{Tr}((G\rho+\rho G^*)A)=\mathrm{Tr}(\rho(G^*A+AG))`$, also nach Korollar 6.9(c)
$`\mathcal L_S(A)=\Phi^*(A)+G^*A+AG`$. Man nehme $`\Psi=\Phi^*`$ (CP nach Proposition 6.8) und $`K=G`$. $`\square`$

Dies ist der endlichdimensionale Fall des nach Satz 8.7 zitierten Satzes von Christensen und Evans.

*Rechnung:* Der Generator aus Korollar D.8 als Dual von $`\mathcal L`$ — [Tutorial, Abschnitt D.4](../../tutorial/de/D_gksl.md#s-D-4).

<a id="s-D-5"></a>

## D.5 Fehlerkorrektur mit Operatoralgebren

Wir verwenden den Rahmen aus Abschnitt 10.3. $`\Phi(\rho)=\sum_{i=1}^rE_i\rho E_i^*`$ ist ein Kanal auf $`M_D`$, $`\mathcal C=P\mathbb C^D`$
ein Code, $`B(\mathcal C)=\{X\in M_D:X=PXP\}`$ die Operatoren auf $`\mathcal C`$, $`\mathfrak A\subseteq B(\mathcal C)`$ eine $`*`$-Algebra
mit Einselement $`P`$ (nach Satz 2.2 eine von-Neumann-Algebra auf $`\mathcal C`$) und $`\mathfrak A'`$ ihr Kommutant in
$`B(\mathcal C)`$. $`\mathfrak A`$ heißt **korrigierbar** für $`\Phi`$, wenn es einen Kanal $`\mathcal R`$ auf $`M_D`$ mit
$`P(\mathcal R\circ\Phi)^*(A)P=A`$ für alle $`A\in\mathfrak A`$ gibt. Hat $`\mathcal R`$ die Kraus-Operatoren $`R_m`$, so hat
$`\mathcal R\circ\Phi`$ die Kraus-Operatoren $`R_mE_i`$, und $`\sum_mR_m^*R_m=1`$.

<a id="r-D-9"></a>

**Satz D.9 (Bény–Kempf–Kribs; Satz 10.7).** $`\mathfrak A`$ ist genau dann korrigierbar für $`\Phi`$, wenn
$`PE_i^*E_jP\in\mathfrak A'`$ für alle $`i,j`$. In diesem Fall korrigiert der im Beweis konstruierte Kanal $`\mathcal R`$ die Algebra $`\mathfrak A`$.

*Beweis.* *Notwendigkeit.* $`\mathcal R`$ korrigiere $`\mathfrak A`$; man setze $`W_{mi}=R_mE_iP`$, so dass
$`\sum_{m,i}W_{mi}^*XW_{mi}=P(\mathcal R\circ\Phi)^*(X)P`$ für alle $`X\in M_D`$. Dies ist gleich $`X`$ für $`X\in\mathfrak A`$ und gleich $`P`$ für
$`X=1`$, weil $`(\mathcal R\circ\Phi)^*`$ unital ist (Korollar 6.9(a)). Sei $`p\in\mathfrak A`$ eine Projektion; dann ist $`p\le P`$
und $`\sum W_{mi}^*(1-p)W_{mi}=P-p`$. Für $`x=px`$ ist
$`\sum_{m,i}\|(1-p)W_{mi}x\|^2=\langle x,(P-p)x\rangle=0`$, und für $`x=(P-p)x`$ ist
$`\sum_{m,i}\|pW_{mi}x\|^2=\langle x,px\rangle=0`$. Daher ist $`(1-p)W_{mi}p=0`$ und $`pW_{mi}(P-p)=0`$, und wegen
$`W_{mi}=W_{mi}P`$ gilt
```math
W_{mi}p=pW_{mi}p=pW_{mi}P=pW_{mi}.
```
Jedes selbstadjungierte $`A\in\mathfrak A`$ ist eine Linearkombination von Projektionen in $`\mathfrak A`$: Nach Satz 1.3 auf $`\mathcal C`$ ist
$`A=\sum_k\lambda_kp_k`$ mit verschiedenen $`\lambda_k`$ und orthogonalen Projektionen $`p_k`$, $`\sum_kp_k=P`$, und
$`p_k=\prod_{l\neq k}(A-\lambda_lP)/(\lambda_k-\lambda_l)\in\mathfrak A`$ (man wende beide Seiten auf die Eigenvektoren in $`\mathcal C`$ und auf $`\mathcal C^\perp`$ an; das leere Produkt ist $`P`$). Jedes Element
von $`\mathfrak A`$ ist $`A_1+iA_2`$ mit selbstadjungierten $`A_1,A_2\in\mathfrak A`$. Also ist $`W_{mi}A=AW_{mi}`$ für alle $`A\in\mathfrak A`$ und,
durch Übergang zu Adjungierten, $`W_{mi}^*A=AW_{mi}^*`$. Da $`\sum_mW_{mi}^*W_{mj}=PE_i^*(\sum_mR_m^*R_m)E_jP=PE_i^*E_jP`$, folgt
```math
A\,PE_i^*E_jP=\sum_mW_{mi}^*AW_{mj}=PE_i^*E_jP\,A\qquad(A\in\mathfrak A),
```
d. h. $`PE_i^*E_jP\in\mathfrak A'`$.

*Hinlänglichkeit.* Nach Satz 2.6 für $`\mathfrak A`$ auf $`\mathcal C`$ gibt es eine Orthonormalbasis $`|c;a,b\rangle`$
($`c=1,\dots,C`$, $`a=1,\dots,k_c`$, $`b=1,\dots,m_c`$) von $`\mathcal C`$ – die Bilder der $`e_a\otimes f_b`$ im $`c`$-ten Summanden unter
dem unitären $`V`$ –, in der $`\mathfrak A`$ aus den Operatoren $`\sum_{c,a,a',b}(A_c)_{aa'}|c;a,b\rangle\langle c;a',b|`$,
$`A_c\in M_{k_c}`$, und $`\mathfrak A'`$ aus den Operatoren $`\sum_{c,a,b,b'}(Y_c)_{bb'}|c;a,b\rangle\langle c;a,b'|`$, $`Y_c\in M_{m_c}`$, besteht.
Die Voraussetzung liefert daher Matrizen $`\Lambda^c_{ij}\in M_{m_c}`$ mit
```math
\langle c;a,b|E_i^*E_j|d;a',b'\rangle=\delta_{cd}\,\delta_{aa'}\,(\Lambda^c_{ij})_{bb'} .
```
*Schritt 1.* Für jedes $`c`$ und jedes Paar $`\alpha=(i,b)`$ sei $`F^c_\alpha:\mathbb C^{k_c}\to\mathbb C^D`$,
$`F^c_\alpha e_a=E_i|c;a,b\rangle`$, und $`\Lambda^c_{\alpha\beta}=(\Lambda^c_{ij})_{bb'}`$ für $`\beta=(j,b')`$. Die obige Gleichung besagt
$`F^{c*}_\alpha F^d_\beta=\delta_{cd}\Lambda^c_{\alpha\beta}1_{k_c}`$: eine Knill–Laflamme-Bedingung für die „Fehler“ $`F^c_\alpha`$ auf dem
Code $`\mathbb C^{k_c}`$, zusätzlich mit Orthogonalität zwischen den Sektoren.

*Schritt 2.* Wie im Beweis von Satz 10.5 ist $`\Lambda^c`$ selbstadjungiert, da
$`\Lambda^c_{\beta\alpha}1=F^{c*}_\beta F^c_\alpha=(F^{c*}_\alpha F^c_\beta)^*=\bar\Lambda^c_{\alpha\beta}1`$, und positiv, da
$`\sum\bar z_\alpha z_\beta\Lambda^c_{\alpha\beta}1=(\sum z_\alpha F^c_\alpha)^*(\sum z_\beta F^c_\beta)\ge0`$. Man schreibe
$`\Lambda^c=U^c\mathrm{diag}(d^c_\kappa)U^{c*}`$ mit unitärem $`U^c`$ und $`d^c_\kappa\ge0`$, und
$`G^c_\kappa=\sum_\alpha u^c_{\alpha\kappa}F^c_\alpha`$. Dann ist $`G^{c*}_\kappa G^d_\lambda=\delta_{cd}\delta_{\kappa\lambda}d^c_\kappa1`$, also
$`G^c_\kappa=0`$ für $`d^c_\kappa=0`$, und andernfalls ist $`V^c_\kappa=G^c_\kappa/\sqrt{d^c_\kappa}`$ eine Isometrie mit
$`V^{c*}_\kappa V^d_\lambda=\delta_{cd}\delta_{\kappa\lambda}1`$. Umkehrung ergibt $`F^c_\alpha=\sum_\kappa\bar u^c_{\alpha\kappa}G^c_\kappa =\sum_{\kappa:d^c_\kappa>0}\bar u^c_{\alpha\kappa}\sqrt{d^c_\kappa}\,V^c_\kappa`$, also
```math
V^{c'*}_\kappa F^d_\beta=\delta_{c'd}\,\bar u^d_{\beta\kappa}\sqrt{d^d_\kappa}\,1,\qquad
F^{c*}_\alpha V^{c'}_\kappa=\delta_{cc'}\,u^c_{\alpha\kappa}\sqrt{d^c_\kappa}\,1 .
```
*Schritt 3: die Wiederherstellung.* Die $`V^c_\kappa V^{c*}_\kappa`$ ($`d^c_\kappa>0`$) sind paarweise orthogonale Projektionen; sei $`Q`$ ihre
Summe. Sei $`\iota_c:\mathbb C^{k_c}\to\mathbb C^D`$, $`\iota_ce_a=|c;a,1\rangle`$, eine Isometrie, und
```math
\mathcal R(\sigma)=\sum_{c}\sum_{\kappa:d^c_\kappa>0}\iota_cV^{c*}_\kappa\,\sigma\,V^c_\kappa\iota_c^*+(1-Q)\sigma(1-Q).
```
Dies ist ein Kanal: $`\sum_{c,\kappa}V^c_\kappa\iota_c^*\iota_cV^{c*}_\kappa+(1-Q)^2=Q+1-Q=1`$. Er misst den Sektor $`c`$ und
das Syndrom $`\kappa`$, macht den Fehler auf $`\mathbb C^{k_c}`$ rückgängig und schreibt das Ergebnis in die Kopie $`b=1`$.

*Schritt 4: Korrektheit.* $`\mathcal R^*(X)=\sum_{c,\kappa}V^c_\kappa\iota_c^*X\iota_cV^{c*}_\kappa+(1-Q)X(1-Q)`$, und für
$`A\in\mathfrak A`$ mit Blöcken $`A_c`$ ist $`\iota_c^*A\iota_c=A_c`$. Da jedes $`F^d_\beta`$ eine Kombination der $`V^d_\kappa`$ ist,
gilt $`(1-Q)F^d_\beta=0`$. Mit Schritt 2 folgt für $`\alpha=(i,b)`$, $`\beta=(j,b')`$
```math
\langle c;a,b|E_i^*\mathcal R^*(A)E_j|d;a',b'\rangle=\sum_{c',\kappa}\langle e_a,F^{c*}_\alpha V^{c'}_\kappa A_{c'}V^{c'*}_\kappa F^d_\beta e_{a'}\rangle
=\delta_{cd}\sum_\kappa u^c_{\alpha\kappa}d^c_\kappa\bar u^c_{\beta\kappa}(A_c)_{aa'}=\delta_{cd}(\Lambda^c_{ij})_{bb'}(A_c)_{aa'} .
```
Man setze $`j=i`$ und summiere über $`i`$: Nach der Gleichung der Voraussetzung mit $`\sum_iE_i^*E_i=1`$ ist $`\sum_i\Lambda^c_{ii}=1_{m_c}`$, also
```math
\langle c;a,b|(\mathcal R\circ\Phi)^*(A)|d;a',b'\rangle=\delta_{cd}\,\delta_{bb'}(A_c)_{aa'}=\langle c;a,b|A|d;a',b'\rangle .
```
Sowohl $`P(\mathcal R\circ\Phi)^*(A)P`$ als auch $`A`$ liegen in $`B(\mathcal C)`$ und haben in der Basis $`|c;a,b\rangle`$ dieselbe Matrix, also sind sie
gleich. $`\square`$

*Rechnung:* Satz D.9 und die Wiederherstellung aus seinem Beweis an einem hybriden Code in $`\mathbb C^8`$ — [Tutorial, Abschnitt D.5](../../tutorial/de/D_gksl.md#s-D-5).

<a id="r-D-10"></a>

**Korollar D.10.** (a) Für $`\mathfrak A=B(\mathcal C)`$ ist $`\mathfrak A'=\mathbb CP`$, und Satz D.9 ist der Satz
10.5 von Knill–Laflamme ($`C=1`$, $`m_1=1`$; Schritt 3 ist die Wiederherstellung aus Satz 10.5). (b) Für $`\mathfrak A=\mathrm{span}\{P_1,\dots,P_C\}`$ mit
orthogonalen Projektionen $`P_c\ne0`$, $`\sum_cP_c=P`$, ist er die Äquivalenz von 2 und 3 in Proposition 10.6. (c) (*Rauschfreie
Teilsysteme.*) In der Notation von Proposition 10.1(a) sei $`\mathcal C=V(\mathbb C^{k_c}\otimes\mathbb C^{m_c})`$ der $`c`$-te
Summand zur Wechselwirkungsalgebra $`\mathcal A`$ von $`\Phi`$ und $`\mathfrak A=V(1_{k_c}\otimes M_{m_c})V^*`$ auf $`\mathcal C`$, eine
Algebra der Form $`M_{m}\otimes1_{k}`$ nach Vertauschen der beiden Faktoren. Dann wird $`\mathfrak A`$ durch
$`\mathcal R=\mathrm{id}`$ korrigiert, und $`PE_i^*E_jP=V(E_{i,c}^*E_{j,c}\otimes1_{m_c})V^*\in\mathfrak A'`$ ($`P=VV^*`$).

*Beweis.* (a) $`\mathfrak A'=\mathbb CP`$ ist der Fall $`C=1`$, $`m_1=1`$ der obigen Beschreibung von $`\mathfrak A'`$, und
$`PE_i^*E_jP\in\mathbb CP`$ ist die Knill–Laflamme-Bedingung. Korrigierbarkeit von $`B(\mathcal C)`$ ist Korrektur im Sinne von
Abschnitt 10.2: Ist $`\mathcal R(\Phi(\rho))=\rho`$ für alle Zustände $`\rho`$ mit Träger in $`\mathcal C`$, so ist
$`\mathrm{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(A)P)=\mathrm{Tr}(\rho A)`$ für alle diese, und diese $`\rho`$ spannen $`B(\mathcal C)`$ auf; die
Umkehrung wird in Abschnitt 10.3 gezeigt. Mit $`C=1`$, $`m_1=1`$ sind die Kraus-Operatoren $`\iota_1V^{1*}_\kappa`$ aus Schritt 3 die
aus Satz 10.5, weil $`\iota_1\iota_1^*=P`$. (b) $`T\in B(\mathcal C)`$ kommutiert genau dann mit allen $`P_c`$, wenn
$`P_cTP_d=0`$ für $`c\ne d`$, und das ist Bedingung 2 aus Proposition 10.6. (c) Jedes $`A\in\mathfrak A`$, durch $`0`$ auf die
anderen Summanden fortgesetzt, liegt in $`\mathcal A'=V(\bigoplus_c1_{k_c}\otimes M_{m_c})V^*`$ (Satz 2.6), also ist $`\Phi^*(A)=A`$ nach Proposition
10.1(a), und $`PAP=A`$. Die Formel für $`PE_i^*E_jP`$ folgt aus $`V^*E_iV=\bigoplus_cE_{i,c}\otimes1_{m_c}`$, und
$`E_{i,c}^*E_{j,c}\otimes1_{m_c}`$ kommutiert mit $`1_{k_c}\otimes M_{m_c}`$. $`\square`$

Nach Satz 2.6 ist jede $`\mathfrak A`$ $`*`$-isomorph zu einer Algebra $`\bigoplus_cM_{k_c}`$. Ist $`C\ge2`$ und $`k_c\ge2`$ für mindestens ein $`c`$, so
schützt eine korrigierbare $`\mathfrak A`$ die Sektorbezeichnung $`c`$ (klassische Information) zusammen mit der Quanteninformation
jedes $`M_{k_c}`$; solche Codes heißen **hybride** Codes.

Die Ergebnisse dieses Anhangs gehen zurück auf V. Gorini, A. Kossakowski und E. C. G. Sudarshan
[*J. Math. Phys.* 17, 821 (1976)] und G. Lindblad [*Commun. Math. Phys.* 48, 119 (1976)] (Sätze D.3–D.5),
E. Christensen und D. E. Evans [*J. London Math. Soc.* 20, 358 (1979)] (Korollar D.8 in allgemeinen von-Neumann-Algebren),
A. Frigerio (Satz D.6) sowie C. Bény, A. Kempf und D. W. Kribs [*Phys. Rev. Lett.* 98, 100502
(2007); *Phys. Rev. A* 76, 042303 (2007)] (Satz D.9); die obigen Beweise sind für die im Kurs verwendeten
endlichdimensionalen Aussagen vollständig.
