[English](../07_tensor_products.md) · [Español](../es/07_tensor_products.md) · [Français](../fr/07_tensor_products.md) · **Deutsch** · [中文](../zh/07_tensor_products.md) · [日本語](../ja/07_tensor_products.md) · [Русский](../ru/07_tensor_products.md)

# 7. Tensorprodukte und Verschränkung

Dieses Kapitel untersucht das Tensorprodukt zweier Algebren mit Sektoren (physikalisch die
Algebra zweier Systeme zusammen; siehe die Kästen *Physikalische Bedeutung*). Mit Sektoren gibt es
zwei Produkte: das einfache, in dem jeder Teil seine eigene Ladung behält, und das *fusionierte*,
in dem nur die Gesamtladung erhalten ist. Wir beweisen die Struktur beider, definieren die
partielle Spur als duale Abbildung einer Einbettung und untersuchen die Verschränkung.

**Sie lernen,** wie sich Sektoren multiplizieren, warum die
partielle Spur ein Kanal und, normiert, eine bedingte Erwartung von endlichem Index ist, wie der
reduzierte Zustand eines ladungserhaltenden Zustands aussieht, wie die Verschränkung in einen
Teilchenzahlanteil und einen Konfigurationsanteil zerfällt, und wie groß die mittlere Reinheit
und die mittlere Entropie zufälliger Zustände sind. Die numerischen Überprüfungen dieses Kapitels
stehen in [Tutorial 7](../../tutorial/de/07_tensor_products.md).

Im Folgenden wirke $`\mathcal M_1=\bigoplus_{c=1}^{C_1}M_{k_c}(\mathbb C)`$ auf
$`H_1=\bigoplus_c\mathbb C^{k_c}`$ mit minimalen zentralen Projektionen $`P_c`$, und
$`\mathcal M_2=\bigoplus_{d=1}^{C_2}M_{l_d}(\mathbb C)`$ wirke auf $`H_2=\bigoplus_d\mathbb C^{l_d}`$
mit minimalen zentralen Projektionen $`Q_d`$. $`\mathrm{Tr}`$ ist die Spur des Raumes, auf dem ein
Operator wirkt; auf $`\mathcal M_i`$ ist sie die ungewichtete Spur aus Kapitel 3. Für $`X`$ auf
$`\mathbb C^k\otimes\mathbb C^l`$ schreiben wir $`X_{(a,b),(a',b')}=\langle e_a\otimes e_b,X\,e_{a'}\otimes e_{b'}\rangle`$
und identifizieren $`\mathbb C^k\otimes\mathbb C^l`$ mit $`\mathbb C^{kl}`$ durch $`e_a\otimes e_b\mapsto e_{(a-1)l+b}`$
(die Kronecker-Konvention), sodass $`(A\otimes B)_{(a,b),(a',b')}=A_{aa'}B_{bb'}`$ gilt.

<a id="s-7-1"></a>

## 7.1 Das Tensorprodukt zweier Algebren

Für $`*`$-Algebren $`\mathcal M_1`$ auf $`H_1`$ und $`\mathcal M_2`$ auf $`H_2`$ ist das **Tensorprodukt**
$`\mathcal M_1\otimes\mathcal M_2`$ die lineare Hülle der Operatoren $`A\otimes B`$ auf
$`H_1\otimes H_2`$ mit $`A\in\mathcal M_1`$, $`B\in\mathcal M_2`$. Wegen
$`(A\otimes B)(A'\otimes B')=AA'\otimes BB'`$, $`(A\otimes B)^*=A^*\otimes B^*`$ und $`1\otimes1=1`$ ist es
eine $`*`$-Algebra und damit eine von-Neumann-Algebra (Satz 2.2).

<a id="r-7-1"></a>

**Satz 7.1.** (a) $`H_1\otimes H_2=\bigoplus_{c,d}\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$ ist die
orthogonale Summe der Bilder der Projektionen $`P_c\otimes Q_d`$, und
```math
\mathcal M_1\otimes\mathcal M_2=\bigoplus_{c,d}M_{k_c}\otimes M_{l_d}\cong\bigoplus_{c,d}M_{k_cl_d}(\mathbb C):
```
es besteht aus allen Operatoren $`X=\sum_{c,d}(P_c\otimes Q_d)X(P_c\otimes Q_d)`$ mit beliebigen Blöcken
auf den $`\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$. (b) Seine minimalen zentralen Projektionen sind die
$`C_1C_2`$ Projektionen $`P_c\otimes Q_d`$, und es gilt
$`\dim(\mathcal M_1\otimes\mathcal M_2)=\dim\mathcal M_1\dim\mathcal M_2`$.
(c) Für Spuren $`\mathrm{tr}_w`$ auf $`\mathcal M_1`$ und $`\mathrm{tr}_v`$ auf $`\mathcal M_2`$ (Satz 3.1) gilt
$`\mathrm{tr}_w(A)\mathrm{tr}_v(B)=\mathrm{tr}_{w\otimes v}(A\otimes B)`$ mit
$`(w\otimes v)_{cd}=w_cv_d`$. Insbesondere sind `Tr_blunt`, `Tr_norm` und `tau_vN` des Produkts die
Produkte derselben Spuren der Faktoren.

*Beweis.* (a) Die $`P_c\otimes Q_d`$ sind paarweise orthogonale Projektionen mit Summe $`1\otimes1`$,
und das Bild von $`P_c\otimes Q_d`$ ist $`\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$. Da $`P_c`$ und $`Q_d`$ zentral
sind, gilt $`(P_c\otimes Q_d)(A\otimes B)=AP_c\otimes BQ_d=(A\otimes B)(P_c\otimes Q_d)`$, also ist jedes
Element von $`\mathcal M_1\otimes\mathcal M_2`$ blockdiagonal. Umgekehrt ist mit den Matrixeinheiten
$`E^c_{aa'}\in\mathcal M_1`$ und $`F^d_{bb'}\in\mathcal M_2`$ der Blöcke $`c`$ und $`d`$
$`E^c_{aa'}\otimes F^d_{bb'}=|e_a\otimes e_b\rangle\langle e_{a'}\otimes e_{b'}|`$ innerhalb von
$`\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$; dies sind die Matrixeinheiten von $`M_{k_cl_d}(\mathbb C)`$, also
liegt jeder im Block $`(c,d)`$ getragene Operator in $`\mathcal M_1\otimes\mathcal M_2`$. (b) Das Zentrum
einer direkten Summe voller Matrixalgebren wird von den Projektionen auf die Summanden aufgespannt
(Kapitel 2, Aufgabe 1), und $`\sum_{c,d}k_c^2l_d^2=\sum_ck_c^2\sum_dl_d^2`$. (c)
$`\mathrm{Tr}(X\otimes Y)=\sum_{a,b}X_{aa}Y_{bb}=\mathrm{Tr}\,X\mathrm{Tr}\,Y`$, also
$`\mathrm{tr}_{w\otimes v}(A\otimes B)=\sum_{c,d}w_cv_d\mathrm{Tr}\,A_c\mathrm{Tr}\,B_d=\mathrm{tr}_w(A)\mathrm{tr}_v(B)`$.
Die Gewichte $`1`$, $`1/k_c`$, $`1/(C_1k_c)`$ der drei Spuren (Abschnitt 3.1) multiplizieren sich zu $`1`$,
$`1/(k_cl_d)`$, $`1/(C_1C_2k_cl_d)`$, den Gewichten derselben Spuren auf dem Produkt. $`\square`$

*In der Bibliothek:* Das Produkt ist `tensor_product` und $`A\otimes B`$ ist `kron`; eine Überprüfung von Satz 7.1 an $`(\mathbb C\oplus M_2)\otimes(M_2\oplus M_3)`$ — [Tutorial, Abschnitt 7.1](../../tutorial/de/07_tensor_products.md#s-7-1).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Gehorcht jedes von zwei Systemen seiner eigenen
> Superauswahlregel – zwei Kästen, die keine Teilchen austauschen können –, so hat das Gesamtsystem einen Sektor für
> jedes *Paar* von Ladungen, und keine Observable verbindet verschiedene Paare. Ein gemeinsamer
> Sektor zählt mit dem Produkt der Gewichte (Satz 7.1(c)).

<a id="s-7-2"></a>

## 7.2 Die partielle Spur

**Definition.** Sei $`\rho\in\mathcal M_1\otimes\mathcal M_2`$. Die **partielle Spur**
$`\mathrm{Tr}_2\,\rho\in\mathcal M_1`$ ist das Element mit
```math
\mathrm{Tr}\big((\mathrm{Tr}_2\,\rho)\,A\big)=\mathrm{Tr}\big(\rho\,(A\otimes1)\big)\qquad\text{for all }A\in\mathcal M_1 .
```
Sie existiert und ist eindeutig nach Satz 3.2, angewandt auf das lineare Funktional
$`A\mapsto\mathrm{Tr}(\rho(A\otimes1))`$ und die treue Spur $`\mathrm{Tr}`$ von $`\mathcal M_1`$: Die
partielle Spur ist die duale Abbildung der Einbettung $`A\mapsto A\otimes1`$ bezüglich der
ungewichteten Spuren. $`\mathrm{Tr}_1`$ wird ebenso definiert; für einen Zustand $`\rho`$ sind
$`\mathrm{Tr}_2\,\rho`$ und $`\mathrm{Tr}_1\,\rho`$ die **reduzierten Zustände**. Für die vollen Algebren
$`\mathcal M_i=B(H_i)`$ ist dies die **partielle Spur des Hilbertraums**
$`\mathrm{Tr}_{H_2}:B(H_1\otimes H_2)\to B(H_1)`$.

<a id="r-7-2"></a>

**Proposition 7.2.** (a) Für $`\rho=\bigoplus_{c,d}\rho_{cd}`$ gilt
```math
(\mathrm{Tr}_2\,\rho)_c=\sum_d\mathrm{tr}_{l_d}\,\rho_{cd},\qquad
(\mathrm{tr}_l\,X)_{aa'}=\sum_{b=1}^lX_{(a,b),(a',b)}\quad(X\text{ on }\mathbb C^k\otimes\mathbb C^l),
```
und $`\mathrm{Tr}_2\,\rho=\mathrm{Tr}_{H_2}\,\rho`$, wobei $`(\mathrm{Tr}_{H_2}\,X)_{ii'}=\sum_jX_{(i,j),(i',j)}`$.
(b) $`\mathrm{Tr}_2(A\otimes B)=(\mathrm{Tr}\,B)\,A`$.
(c) $`\mathrm{Tr}_2\big((A\otimes1)X(A'\otimes1)\big)=A\,(\mathrm{Tr}_2\,X)\,A'`$ für $`A,A'\in\mathcal M_1`$.
(d) $`\mathrm{Tr}(\mathrm{Tr}_2\,X)=\mathrm{Tr}\,X`$, und $`\mathrm{Tr}_2`$ bildet Zustände auf Zustände ab.

*Beweis.* (a) Mit $`A=\bigoplus_cA_c`$ ist $`A\otimes1=\bigoplus_{c,d}A_c\otimes1_{l_d}`$ und
```math
\mathrm{Tr}\big(\rho(A\otimes1)\big)=\sum_{c,d}\sum_{a,a',b}(\rho_{cd})_{(a,b),(a',b)}(A_c)_{a'a}
=\sum_c\mathrm{Tr}\Big(\Big(\sum_d\mathrm{tr}_{l_d}\,\rho_{cd}\Big)A_c\Big),
```
also erfüllt das Element auf der rechten Seite die Definition und ist gleich $`\mathrm{Tr}_2\,\rho`$.
Dieselbe Rechnung für $`B(H_i)`$ liefert die Formel für $`\mathrm{Tr}_{H_2}`$; für blockdiagonales
$`\rho`$ verschwinden die Einträge $`\rho_{(i,j),(i',j)}`$ mit $`i,i'`$ in verschiedenen Blöcken, also hat
$`\mathrm{Tr}_{H_2}\,\rho`$ die Blöcke aus (a). (b) $`\mathrm{tr}_l(X\otimes Y)=(\mathrm{Tr}\,Y)X`$. (c) Für
$`B\in\mathcal M_1`$ ist
$`\mathrm{Tr}\big(\mathrm{Tr}_2((A\otimes1)X(A'\otimes1))B\big)=\mathrm{Tr}\big(X(A'BA\otimes1)\big) =\mathrm{Tr}\big(A(\mathrm{Tr}_2\,X)A'B\big)`$; man verwende die Eindeutigkeit. (d) Man nehme $`A=1`$. Ist
$`\rho`$ ein Zustand, so ist $`A\mapsto\mathrm{Tr}(\rho(A\otimes1))`$ positiv, da
$`A^*A\otimes1=(A\otimes1)^*(A\otimes1)`$, und hat bei $`A=1`$ den Wert $`1`$; also ist seine Dichte
$`\mathrm{Tr}_2\,\rho`$ ein Zustand (Satz 3.2). $`\square`$

*In der Bibliothek:* $`\mathrm{Tr}_2`$ ist `partial_trace`; eine Überprüfung von Proposition 7.2(a), (b) gegen `numpy.einsum` — [Tutorial, Abschnitt 7.2](../../tutorial/de/07_tensor_products.md#s-7-2).

Kanäle sind das Thema von Kapitel 6; hier brauchen wir nur diese Definition: Eine lineare
Abbildung $`\Phi`$ heißt **vollständig positiv**, wenn $`\mathrm{id}_n\otimes\Phi`$ für jedes $`n`$ auf
$`M_n(\mathbb C)\otimes(\cdot)`$ positiv ist, und ein **Kanal**, wenn sie außerdem spurerhaltend ist.

<a id="r-7-3"></a>

**Satz 7.3.** Es laufe $`j`$ über die Standardbasis von $`H_2`$, und es sei
$`K_j=1\otimes\langle j|:H_1\otimes H_2\to H_1`$, $`K_j(x\otimes y)=\langle j,y\rangle x`$. Dann gilt
```math
\mathrm{Tr}_{H_2}\,X=\sum_jK_jXK_j^*,\qquad\sum_jK_j^*K_j=1,\qquad\sum_jK_j^*AK_j=A\otimes1 .
```
Folglich ist $`\mathrm{Tr}_2`$ ein Kanal von $`\mathcal M_1\otimes\mathcal M_2`$ nach $`\mathcal M_1`$, und
seine duale Abbildung bezüglich der ungewichteten Spuren ist $`A\mapsto A\otimes1`$.

*Beweis.* $`K_j^*x=x\otimes e_j`$, also $`\langle e_a,K_jXK_j^*e_{a'}\rangle=X_{(a,j),(a',j)}`$; Summation
über $`j`$ liefert Proposition 7.2(a). Ferner ist
$`\sum_jK_j^*K_j(x\otimes y)=\sum_j\langle j,y\rangle x\otimes e_j=x\otimes y`$
und $`\sum_jK_j^*AK_j=\sum_jA\otimes|j\rangle\langle j|=A\otimes1`$. Für $`Y\ge0`$ ist
$`(\mathrm{id}_n\otimes\mathrm{Tr}_{H_2})(Y)=\sum_j(1_n\otimes K_j)Y(1_n\otimes K_j)^*\ge0`$, denn
$`\langle x,TYT^*x\rangle=\langle T^*x,YT^*x\rangle\ge0`$; dies bleibt auf der Unteralgebra
$`\mathcal M_1\otimes\mathcal M_2`$ gültig. Schließlich ist
$`\mathrm{Tr}\,\sum_jK_jXK_j^*A=\mathrm{Tr}\,X\sum_jK_j^*AK_j`$, was für $`A=1`$ die Spurerhaltung und im
Allgemeinen die Dualität ist. $`\square`$

*In der Bibliothek:* Der Kanal aus Satz 7.3 ist `composite.partial_trace_channel`, mit seiner dualen Abbildung — [Tutorial, Abschnitt 7.2](../../tutorial/de/07_tensor_products.md#s-7-2).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Der reduzierte Zustand $`\mathrm{Tr}_2\,\rho`$
> ist alles, was ein Beobachter erfährt, der nur das erste System messen kann – die Definition besagt, dass er jeden
> Erwartungswert einer Observablen $`A\otimes1`$ reproduziert. Gemeinsame Zustände mit gleichen
> reduzierten Zuständen unterscheiden sich in ihren Korrelationen.

**Die partielle Spur als bedingte Erwartung.** Sei $`\mathcal M_2=M_m(\mathbb C)`$, sodass die Sektoren
von $`\mathcal M_1\otimes M_m`$ die von $`\mathcal M_1`$ sind (Satz 7.1), sei $`\iota(A)=A\otimes1_m`$ und
$`\mathcal E=\frac1m\mathrm{Tr}_2:\mathcal M_1\otimes M_m\to\mathcal M_1`$, und seien $`e_{ab}`$ die
Matrixeinheiten von $`M_m`$.

<a id="r-7-4"></a>

**Proposition 7.4.** (a) $`\iota\circ\mathcal E`$ ist eine bedingte Erwartung von
$`\mathcal M_1\otimes M_m`$ auf $`\mathcal M_1\otimes1`$ (Abschnitt 2.6), und für jede treue Spur
$`\mathrm{tr}_w`$ von $`\mathcal M_1\otimes M_m`$ ist sie die einzige, die $`\mathrm{tr}_w`$ erhält. (b) Die
$`m^2`$ Elemente $`u_{ab}=\sqrt m\,(1\otimes e_{ab})`$ erfüllen
$`X=\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*X))`$ für alle $`X`$, und es gilt
$`\sum_{a,b}u_{ab}u_{ab}^*=m^2\,1`$. (c) Jede Familie $`(v_k)`$ mit
$`X=\sum_kv_k\,\iota(\mathcal E(v_k^*X))`$ für alle $`X`$ erfüllt $`\sum_kv_kv_k^*=m^2\,1`$.

*Beweis.* (a) $`\mathcal E`$ ist positiv (Satz 7.3), ebenso $`\iota`$
($`T^*T\otimes1=(T\otimes1)^*(T\otimes1)`$); nach Proposition 7.2(b) gilt $`\mathcal E(A\otimes1)=A`$, und die
Bimoduleigenschaft ist Proposition 7.2(c). Wir schreiben $`\mathrm{tr}_w`$ auch für die Spur mit
denselben Gewichten $`w_c`$ auf $`\mathcal M_1`$; dann ist $`\mathrm{tr}_w(Y\otimes1)=m\mathrm{tr}_w(Y)`$. Die
Definition von $`\mathrm{Tr}_2`$, angewandt auf $`AP_c`$, mit $`w_c`$ multipliziert und über $`c`$ summiert,
ergibt $`\mathrm{tr}_w(X(A\otimes1))=\mathrm{tr}_w((\mathrm{Tr}_2\,X)A)`$, und $`A=1`$ ergibt
$`\mathrm{tr}_w(\iota\mathcal E(X))=\mathrm{tr}_w(X)`$. Ist $`F=\iota\circ F_0`$ eine weitere bedingte Erwartung
auf $`\mathcal M_1\otimes1`$, die $`\mathrm{tr}_w`$ erhält, so gilt wegen
$`F(X(A\otimes1))=F(X)(A\otimes1)`$
$`\mathrm{tr}_w(X(A\otimes1))=\mathrm{tr}_w(F(X)(A\otimes1))=m\mathrm{tr}_w(F_0(X)A)`$, also
$`\mathrm{tr}_w((\mathrm{Tr}_2\,X-mF_0(X))A)=0`$ für alle $`A\in\mathcal M_1`$, und wegen der Treue ist
$`F_0=\mathcal E`$ (man nehme $`A=(\mathrm{Tr}_2\,X-mF_0(X))^*`$). (b) Man schreibe
$`X=\sum_{a,b}X_{ab}\otimes e_{ab}`$ mit $`X_{ab}\in\mathcal M_1`$ (Satz 7.1(a)). Nach Proposition 7.2(b) ist
$`\mathcal E(u_{ab}^*X)=\frac{\sqrt m}m\sum_{b'}\mathrm{Tr}_2(X_{ab'}\otimes e_{bb'})=X_{ab}/\sqrt m`$, also
$`u_{ab}\,\iota(\mathcal E(u_{ab}^*X))=X_{ab}\otimes e_{ab}`$; und
$`\sum_{a,b}u_{ab}u_{ab}^*=m\sum_{a,b}1\otimes e_{aa}=m^21`$. (c) Nach Proposition 7.2(a) ist
$`\mathrm{Tr}_2(Y^*)=(\mathrm{Tr}_2\,Y)^*`$, also ist die Adjungierte der Entwicklung von $`Y=u_{ab}`$ gleich
$`u_{ab}^*=\sum_k\iota(\mathcal E(u_{ab}^*v_k))\,v_k^*`$, und nach (b) gilt
$`\sum_{a,b}u_{ab}u_{ab}^*=\sum_k\big(\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*v_k))\big)v_k^*=\sum_kv_kv_k^*`$. $`\square`$

Eine Familie wie in (c) ist eine *Quasibasis* von $`\mathcal E`$, und $`\sum_kv_kv_k^*`$, das nach (c)
nicht von ihr abhängt und nach (b) existiert, ist der **Index** $`\mathrm{Ind}\,\mathcal E=m^2`$ von
$`\mathcal E`$ (eine Definition; der Name folgt Y. Watatani, Mem. Amer. Math. Soc. 83, no. 424, 1990).
Die **Inklusionsmatrix** einer Inklusion von $`\bigoplus_cM_{k_c}`$ in $`\bigoplus_sM_{n_s}`$ hat die
Einträge $`\Lambda_{sc}`$, die Anzahlen der Kopien des Blocks $`c`$ der Unteralgebra im Block $`s`$ der
größeren (Kapitel 9 verwendet sie für Türme von Algebren). Für einen Faktor,
$`M_k\otimes1\subset M_k\otimes M_m`$: Ordnet man die Basis $`e_a\otimes e_b`$ von $`\mathbb C^{km}`$ zuerst nach
$`b`$, so ist $`A\otimes1_m=\mathrm{diag}(A,\dots,A)`$ mit $`m`$ Kopien, also $`\Lambda=(m)`$, und
$`m^2=\|\Lambda\|^2=\dim(M_k\otimes M_m)/\dim M_k`$. (Weiterführende Literatur, hier nicht verwendet:
V. F. R. Jones, Invent. Math. 72, 1, 1983, zum Index von Subfaktoren, und F. M. Goodman, P. de la
Harpe und V. F. R. Jones, *Coxeter graphs and towers of algebras*, Springer, 1989, zu
Inklusionsmatrizen.) Für allgemeines $`\mathcal M_1`$ ist die
Inklusion die direkte Summe der $`C_1`$ Inklusionen $`M_{k_c}\otimes1\subset M_{k_c}\otimes M_m`$,
$`\Lambda=m1_{C_1}`$.

*In der Bibliothek:* `inclusions.from_tensor_product`, mit der bedingten Erwartung, der Quasibasis und dem Index aus Proposition 7.4 — [Tutorial, Abschnitt 7.2](../../tutorial/de/07_tensor_products.md#s-7-2).

<a id="s-7-3"></a>

## 7.3 Fusion von Ladungen

Nun trage der Sektor $`c`$ von $`\mathcal M_1`$ eine Ladung $`q_c`$, der Sektor $`d`$ von $`\mathcal M_2`$
eine Ladung $`r_d`$, und eine **Fusionsregel** $`f`$ liefere die Gesamtladung $`f(q_c,r_d)`$ – die
Addition für Teilchenzahlen, die Addition modulo $`2`$ für die Fermionenparität oder irgendeine
andere Abbildung. Man setze
```math
H_q=\bigoplus_{f(q_c,r_d)=q}\mathbb C^{k_c}\otimes\mathbb C^{l_d},\qquad
\Pi_q=\sum_{f(q_c,r_d)=q}P_c\otimes Q_d,\qquad n_q=\sum_{f(q_c,r_d)=q}k_cl_d .
```
Das **fusionierte Produkt** $`\mathcal M_1\otimes_f\mathcal M_2=\{\Pi_q\}'`$ besteht aus den Operatoren
auf $`H_1\otimes H_2`$, die die Gesamtladung erhalten.

<a id="r-7-5"></a>

**Proposition 7.5.** (a) $`\mathcal M_1\otimes_f\mathcal M_2=\bigoplus_qB(H_q)\cong\bigoplus_qM_{n_q}(\mathbb C)`$
mit den minimalen zentralen Projektionen $`\Pi_q`$. (b) $`\mathcal M_1\otimes\mathcal M_2\subseteq\mathcal M_1\otimes_f\mathcal M_2`$,
mit Gleichheit genau dann, wenn jedes $`H_q`$ nur ein einziges Paar $`(c,d)`$ enthält. (c) Sind die
Ladungen reelle Zahlen und ist $`f`$ die Addition, so gilt $`\mathcal M_1\otimes_f\mathcal M_2=\{N\}'`$ mit
$`N=N_1\otimes1+1\otimes N_2`$, $`N_1=\sum_cq_cP_c`$, $`N_2=\sum_dr_dQ_d`$.

*Beweis.* (a) Jedes Paar liegt in genau einem $`H_q`$, also sind die $`\Pi_q`$ paarweise orthogonale
Projektionen mit Summe $`1`$, und $`X`$ kommutiert genau dann mit allen, wenn $`X=\sum_q\Pi_qX\Pi_q`$
gilt, mit beliebigen Blöcken in $`B(H_q)`$; das Zentrum ergibt sich wie in Satz 7.1(b). (b) Jedes
Element von $`\mathcal M_1\otimes\mathcal M_2`$ kommutiert mit jedem $`P_c\otimes Q_d`$, also auch mit
deren Summen $`\Pi_q`$. Ist jedes $`H_q`$ ein einzelnes $`\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$, so sind
beide Algebren gleich $`\bigoplus_{c,d}B(\mathbb C^{k_c}\otimes\mathbb C^{l_d})`$; enthält $`H_q`$ zwei
Paare, so enthält $`B(H_q)`$ Operatoren, die den Block des einen Paares in den des anderen abbilden
und nicht blockdiagonal sind. (c) $`N=\sum_{c,d}(q_c+r_d)P_c\otimes Q_d=\sum_qq\,\Pi_q`$, also sind die
von null verschiedenen $`\Pi_q`$ die Spektralprojektionen von $`N`$. Sie sind Polynome in $`N`$ (Beweis
von Satz 1.5(a)), und $`N`$ ist eine Linearkombination von ihnen; daher kommutiert $`X`$ genau dann mit
$`N`$, wenn es mit jedem $`\Pi_q`$ kommutiert. $`\square`$

**Beispiel: Hard-Core-Bosonen.** Auf $`L`$ Plätzen, $`(\mathbb C^2)^{\otimes L}`$
mit der Besetzungszahl $`N=\sum_in^{(i)}`$, $`n=\mathrm{diag}(0,1)`$, wird der Eigenraum $`N=q`$ von den
$`\binom Lq`$ Produktbasisvektoren mit $`q`$ besetzten Plätzen aufgespannt, also bilden die
teilchenzahlerhaltenden Operatoren $`\{N\}'`$ die Algebra $`\bigoplus_qM_{\binom Lq}`$ (Proposition 4.10(2)):
für zwei Plätze $`M_1\oplus M_2\oplus M_1`$ (Ladungen $`0,1,2`$), für drei Plätze
$`M_1\oplus M_3\oplus M_3\oplus M_1`$. Nach Proposition 7.5(c) ist das fusionierte Produkt $`\{N\}'`$ für
die Gesamtbesetzungszahl von fünf Plätzen, bis auf die Reihenfolge der Basis, und tatsächlich ist
$`n_q=\sum_c\binom2c\binom3{q-c}=\binom5q`$ (Vandermondesche Identität: Eine $`q`$-elementige Teilmenge
von fünf Plätzen besteht aus $`c`$ der ersten beiden und $`q-c`$ der letzten drei).

*In der Bibliothek:* Das fusionierte Produkt ist `fused_tensor_product`; die Hard-Core-Bosonen und Proposition 7.5 — [Tutorial, Abschnitt 7.3](../../tutorial/de/07_tensor_products.md#s-7-3).

Die fusionierte Algebra hat die Dimension $`\sum_q\binom5q^2=\binom{10}5=252`$ (wieder Vandermonde),
das einfache Produkt $`(1+4+1)(1+9+9+1)=120`$ (Satz 7.1(b)). Nach Satz 7.1(a) und Proposition 7.5(b)
besteht das einfache Produkt aus den Elementen des fusionierten, die bezüglich der Paare $`(c,d)`$
blockdiagonal sind; die übrigen $`132`$ Dimensionen werden von den Operatoren aufgespannt, die einen
Paarblock $`(c,d)`$ in einen anderen Paarblock $`(c',d')`$ desselben $`H_q`$ abbilden, $`c'\ne c`$ (also
$`d'\ne d`$): Sie ändern die Ladung jedes Teils bei erhaltener Gesamtladung, d. h. hier verschieben
sie Teilchen zwischen den beiden Teilen.

Ein Zustand von $`\mathcal M_1\otimes_f\mathcal M_2`$ kann Kohärenzen zwischen den Blöcken
verschiedener Paare haben, sodass Proposition 7.2 nicht anwendbar ist. Sei $`E_1(X)=\sum_cP_cXP_c`$
auf $`B(H_1)`$ (die bedingte Erwartung aus Abschnitt 2.6); sie bildet auf
$`\mathcal M_1`$ ab, und es gilt $`\mathrm{Tr}(E_1(X)A)=\mathrm{Tr}(XA)`$ für $`A\in\mathcal M_1`$, da
$`A=\sum_cP_cAP_c`$.

<a id="r-7-6"></a>

**Satz 7.6.** Sei $`\rho`$ die Dichte eines Zustands $`\omega`$ auf $`\mathcal M_1\otimes_f\mathcal M_2`$. Die
Dichte seiner Einschränkung $`A\mapsto\omega(A\otimes1)`$ auf $`\mathcal M_1`$ ist
```math
\rho_1=E_1\big(\mathrm{Tr}_{H_2}\,\rho\big).
```
Der Block $`(c,c')`$ von $`\mathrm{Tr}_{H_2}\,\rho`$ ist $`\sum_d\mathrm{tr}_{l_d}\,\rho_{(c,d),(c',d)}`$, summiert
über die $`d`$ mit $`f(q_c,r_d)=f(q_{c'},r_d)`$, wobei $`\rho_{(c,d),(c',d)}`$ der Block von $`\rho`$ vom
Paar $`(c',d)`$ zum Paar $`(c,d)`$ ist. Folglich wirkt $`E_1`$ genau dann auf allen $`\mathrm{Tr}_{H_2}\,\rho`$
trivial, wenn es keine $`c\neq c'`$ und $`d`$ mit $`f(q_c,r_d)=f(q_{c'},r_d)`$ gibt; zum Beispiel, wenn $`f`$
die Addition ist und die $`q_c`$ paarweise verschieden sind.

*Beweis.* $`A\otimes1\in\mathcal M_1\otimes_f\mathcal M_2`$ (Proposition 7.5(b)) und
$`\omega(A\otimes1)=\mathrm{Tr}(\rho(A\otimes1))=\mathrm{Tr}\big((\mathrm{Tr}_{H_2}\,\rho)A\big) =\mathrm{Tr}\big(E_1(\mathrm{Tr}_{H_2}\,\rho)A\big)`$. Da $`E_1(\mathrm{Tr}_{H_2}\,\rho)\in\mathcal M_1`$ ist, ist es
nach der Eindeutigkeit in Satz 3.2 die Dichte. Die Blockformel ist Proposition 7.2(a) für
$`\mathrm{Tr}_{H_2}`$, mit $`i`$ im Block $`c`$, $`i'`$ im Block $`c'`$ und $`j`$ im Block $`d`$; der Block von $`\rho`$
zwischen $`(c',d)`$ und $`(c,d)`$ verschwindet, außer wenn beide Paare im selben $`H_q`$ liegen. Gibt es
$`c\neq c'`$ und $`d`$ mit $`f(q_c,r_d)=f(q_{c'},r_d)=q`$, so nehme man Einheitsvektoren
$`x\in\mathbb C^{k_c}`$, $`x'\in\mathbb C^{k_{c'}}`$, $`y\in\mathbb C^{l_d}`$ und
$`\psi=(x+x')\otimes y/\sqrt2\in H_q`$. Dann ist $`|\psi\rangle\langle\psi|\in B(H_q)`$, und
$`\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|=\frac12|x+x'\rangle\langle x+x'|`$ hat den von null verschiedenen
Block $`\frac12|x\rangle\langle x'|`$, den $`E_1`$ entfernt. $`\square`$

Nach der Blockformel von Satz 7.6 ist $`(\rho_1)_c=\sum_d\mathrm{tr}_{l_d}\,\rho_{(c,d),(c,d)}`$.

*In der Bibliothek:* `partial_trace` und `composite.partial_trace_channel` berechnen $`\rho_1`$ für jede Fusionsregel; eine Überprüfung von Satz 7.6 für die Addition und für die Paritätsregel — [Tutorial, Abschnitt 7.3](../../tutorial/de/07_tensor_products.md#s-7-3).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). In der BCS-Beschreibung eines
> Supraleiters ist nur die Parität der Elektronenzahl erhalten, sodass der globale Zustand die Teilchenzahlen $`0`$ und $`2`$ eines Gebiets
> überlagern kann. Ein Beobachter, dessen Observablen die Teilchenzahl des Gebiets erhalten, kann
> diese Kohärenz nicht sehen: Für ihn ist der reduzierte Zustand der blockdiagonal gemachte Zustand
> aus Satz 7.6 (der mathematische Inhalt dieses Kastens), eine Mischung von Teilchenzahlen.

<a id="s-7-4"></a>

## 7.4 Verschränkung reiner Zustände

<a id="r-7-7"></a>

**Lemma 7.7 (Schmidt-Zerlegung).** Jeder Einheitsvektor $`\psi\in\mathbb C^a\otimes\mathbb C^b`$ lässt
sich als $`\psi=\sum_{i=1}^rs_i\,u_i\otimes v_i`$ schreiben, mit $`s_i>0`$, $`\sum_is_i^2=1`$ und
orthonormalen $`u_i\in\mathbb C^a`$, $`v_i\in\mathbb C^b`$. Dann ist
$`\mathrm{Tr}_2|\psi\rangle\langle\psi|=\sum_is_i^2|u_i\rangle\langle u_i|`$
und $`\mathrm{Tr}_1|\psi\rangle\langle\psi|=\sum_is_i^2|v_i\rangle\langle v_i|`$: Die beiden reduzierten
Zustände haben dieselben von null verschiedenen Eigenwerte $`s_i^2`$.

*Beweis.* Sei $`\Psi`$ die $`a\times b`$-Matrix $`\Psi_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,\psi\rangle`$.
Nach dem Spektralsatz ist $`\Psi\Psi^*=\sum_{i=1}^a\lambda_i|u_i\rangle\langle u_i|`$ mit einer
Orthonormalbasis $`u_i`$ und $`\lambda_i=\|\Psi^*u_i\|^2\ge0`$, so angeordnet, dass $`\lambda_i>0`$ genau
für $`i\le r`$ gilt. Man setze $`s_i=\sqrt{\lambda_i}`$ und $`w_i=\Psi^*u_i/s_i`$ für $`i\le r`$; dann ist
$`\langle w_i,w_j\rangle=\langle u_i,\Psi\Psi^*u_j\rangle/(s_is_j)=\delta_{ij}`$, und $`\Psi^*u_i=0`$ für $`i>r`$.
Also ist $`\Psi=\sum_i|u_i\rangle\langle u_i|\Psi=\sum_{i\le r}|u_i\rangle\langle\Psi^*u_i|=\sum_{i\le r}s_i|u_i\rangle\langle w_i|`$,
d. h. $`\Psi_{\alpha\beta}=\sum_is_i(u_i)_\alpha\overline{(w_i)_\beta}`$: Das ist die Behauptung mit
$`v_i=\bar w_i`$, und $`\sum_is_i^2=\mathrm{Tr}\,\Psi\Psi^*=\|\psi\|^2=1`$. Schließlich ist
$`|\psi\rangle\langle\psi|=\sum_{i,j}s_is_j|u_i\rangle\langle u_j|\otimes|v_i\rangle\langle v_j|`$, und
Proposition 7.2(b) liefert die reduzierten Zustände. $`\square`$

Die **Verschränkungsentropie** von $`\psi`$ ist
$`E(\psi)=S(\mathrm{Tr}_2|\psi\rangle\langle\psi|)=S(\mathrm{Tr}_1|\psi\rangle\langle\psi|)=-\sum_is_i^2\log s_i^2`$.
Sie verschwindet genau dann, wenn $`r=1`$ (jeder Term $`-s_i^2\log s_i^2`$ ist $`\ge0`$ und nur für
$`s_i^2=1`$ gleich $`0`$), d. h. genau dann, wenn $`\psi=u\otimes v`$ (dann hat $`\Psi=uv^T`$ Rang $`1`$, und
$`r=\mathrm{rank}\,\Psi\Psi^*=\mathrm{rank}\,\Psi`$), und es gilt $`E(\psi)\le\log\min(a,b)`$, mit Gleichheit
genau dann, wenn der reduzierte Zustand auf dem kleineren Faktor maximal gemischt ist
(Korollar 3.6, angewandt auf beide reduzierten Zustände).

Mit Ladungen nehmen wir einen Einheitsvektor $`\psi`$ mit bestimmter Gesamtladung, d. h. in einem
Sektor $`H_q`$ des fusionierten Produkts (eine Überlagerung verschiedener $`H_q`$ liefert auf
$`\mathcal M_1\otimes_f\mathcal M_2`$ denselben Zustand wie eine Mischung, Kapitel 3, Aufgabe 3); wir
schreiben $`\psi=\sum_{(c,d)}\psi_{cd}`$ mit $`\psi_{cd}\in\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$.

<a id="r-7-8"></a>

**Proposition 7.8 (symmetrieaufgelöste Verschränkung).** Sei $`f`$ die Addition, und seien die $`q_c`$
paarweise verschieden und die $`r_d`$ paarweise verschieden, sodass es zu jedem $`c`$ höchstens ein
$`d=d(c)`$ mit $`q_c+r_d=q`$ gibt. Man setze $`\psi_c=\psi_{c\,d(c)}`$ (bzw. $`0`$), $`p_c=\|\psi_c\|^2`$ und,
für $`p_c>0`$, $`\hat\psi_c=\psi_c/\sqrt{p_c}`$. Dann liegt $`\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|`$ in
$`\mathcal M_1`$ (es kommutiert mit $`N_1`$) und ist gleich
$`\bigoplus_cp_c\mathrm{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|`$; der reduzierte Zustand auf $`\mathcal M_2`$
hat bei $`d(c)`$ die Sektorwahrscheinlichkeit $`p_c`$; und
```math
E(\psi)=H(p)+\sum_cp_c\,E(\hat\psi_c),\qquad H(p)=-\sum_cp_c\log p_c .
```

*Beweis.* Nach der Blockformel aus Satz 7.6 ist der Block $`(c,c')`$ von
$`\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|`$ die Summe der $`\mathrm{tr}_{l_d}|\psi_{cd}\rangle\langle\psi_{c'd}|`$
über die $`d`$ mit $`q_c+r_d=q=q_{c'}+r_d`$, was $`c=c'`$ und $`d=d(c)`$ erzwingt. Der reduzierte Zustand
ist also $`\bigoplus_c\mathrm{Tr}_2|\psi_c\rangle\langle\psi_c|`$ mit den Blockspuren $`p_c`$; Vertauschen der
Rollen der Faktoren liefert die Aussage für $`\mathcal M_2`$. Nach Lemma 7.7 ist $`E(\psi)`$ die
Entropie dieses Zustands, und Satz 3.3 mit $`\sigma_c=\mathrm{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|`$,
$`S(\sigma_c)=E(\hat\psi_c)`$, liefert die Formel. $`\square`$

$`H(p)`$ ist die **Teilchenzahlverschränkung** und $`\sum_cp_cE(\hat\psi_c)`$ die
**Konfigurationsverschränkung**.

*Rechnung:* Lemma 7.7 und Proposition 7.8 für einen zufälligen Zustand von zwei Teilchen auf fünf Plätzen — [Tutorial, Abschnitt 7.4](../../tutorial/de/07_tensor_products.md#s-7-4).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Parteien, die an die Superauswahlregel der
> Teilchenzahl gebunden sind, können nur den Konfigurationsanteil als Ressource nutzen: die *Verschränkung von Teilchen*
> nach H. M. Wiseman und J. A. Vaccaro (Phys. Rev. Lett. 91, 097902, 2003). $`H(p)`$ wird in kalten
> Atomen als *Teilchenzahlentropie* (number entropy) gemessen (A. Lukin et al., Science 364, 256,
> 2019).

<a id="s-7-5"></a>

## 7.5 Zufällige Zustände und die Pagesche Formel

Sei $`g\in\mathbb C^N`$ mit unabhängigen Koordinaten, deren Real- und Imaginärteile unabhängig
$`N(0,\tfrac12)`$-verteilt sind (ein komplexer Standard-Gauß-Vektor, Proposition A.16(a)), sodass
$`\mathbb E|g_i|^2=1`$; fast sicher ist $`g\neq0`$ (Satz A.20). Ein **Haar-zufälliger reiner Zustand** ist
$`\psi=g/\|g\|`$. Seine Verteilung ist unitär invariant: $`Ug`$ hat für jedes unitäre $`U`$ die Verteilung
von $`g`$ (Proposition A.16(a)), und $`Ug/\|Ug\|=U\psi`$. Nach Proposition A.22 ist sie die einzige unitär
invariante Verteilung auf der Einheitssphäre, und, wie in deren Beweis gezeigt, ist sie die
Verteilung der ersten Spalte einer Haar-verteilten unitären Matrix (Satz 5.3).

<a id="r-7-9"></a>

**Lemma 7.9.** $`\|g\|`$ und $`\psi=g/\|g\|`$ sind unabhängig.

*Beweis.* (Dies ist ein Teil von Satz A.20(b).) Seien $`h\ge0`$ und $`\phi\ge0`$ Borel-Funktionen auf
$`(0,\infty)`$ bzw. auf der Einheitssphäre $`S^{2N-1}`$ von $`\mathbb C^N\cong\mathbb R^{2N}`$. Da $`g`$ die Dichte
$`\pi^{-N}e^{-|x|^2}`$ hat (Proposition A.16(a)), liefern Polarkoordinaten (Satz A.18 mit $`n=2N`$ und dem
dort definierten Maß $`\mu`$)
```math
\mathbb E\big[h(\|g\|)\phi(\psi)\big]=\pi^{-N}\!\int_{\mathbb R^{2N}}\!h(|x|)e^{-|x|^2}\phi\Big(\frac x{|x|}\Big)dx
=2N\pi^{-N}\int_0^\infty h(t)e^{-t^2}t^{2N-1}dt\int_{S^{2N-1}}\phi\,d\mu=:\alpha(h)\beta(\phi).
```
Setzt man $`\phi=1`$, dann $`h=1`$, dann beides, so erhält man $`\mathbb Eh(\|g\|)=\alpha(h)\beta(1)`$,
$`\mathbb E\phi(\psi)=\alpha(1)\beta(\phi)`$ und $`\alpha(1)\beta(1)=1`$, also
$`\mathbb E[h(\|g\|)\phi(\psi)]=\mathbb Eh(\|g\|)\,\mathbb E\phi(\psi)`$. $`\square`$

<a id="r-7-10"></a>

**Lemma 7.10.** $`\mathbb E\|g\|^4=N(N+1)`$. Ist $`G`$ eine $`a\times b`$-Matrix mit unabhängigen, wie die
$`g_i`$ verteilten Einträgen (eine *Ginibre*-Matrix), so ist $`\mathbb E\mathrm{Tr}(GG^*)^2=ab(a+b)`$.

*Beweis.* Mit $`g_i=x+iy`$ und unabhängigen $`N(0,\frac12)`$-verteilten $`x,y`$ ist $`\mathbb Ex^4=3(\mathbb Ex^2)^2=\frac34`$
(Satz A.24: drei Paarungen), also
$`\mathbb E|g_i|^4=\mathbb Ex^4+2\mathbb Ex^2\mathbb Ey^2+\mathbb Ey^4=\tfrac34+\tfrac12+\tfrac34=2`$ und
$`\mathbb E\|g\|^4=\sum_i\mathbb E|g_i|^4+\sum_{i\ne j}\mathbb E|g_i|^2\mathbb E|g_j|^2=2N+N(N-1)`$. Weiter
ist $`\mathrm{Tr}(GG^*)^2=\sum_{i,i',j,j'}G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}`$. Die Einträge und ihre
Konjugierten sind komplexe Linearkombinationen der Koordinaten eines einzigen zentrierten
Gauß-Vektors (der Real- und Imaginärteile aller Einträge), mit
$`\mathbb E[G_{ij}G_{kl}]=\mathbb E[\bar G_{ij}\bar G_{kl}]=0`$ und $`\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}`$
(Proposition A.16(a)). Nach der Isserlis-Wick-Formel (Satz A.24, komplexe Form, mit den drei
Paarungen der vier Faktoren) ist
$`\mathbb E[G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}]=\delta_{ii'}\delta_{i'i}+0+\delta_{jj'}\delta_{jj'}=\delta_{ii'}+\delta_{jj'}`$.
Summation über $`i,i'\le a`$ und $`j,j'\le b`$ ergibt $`ab^2+a^2b`$. $`\square`$

<a id="r-7-11"></a>

**Satz 7.11 (E. Lubkin, J. Math. Phys. 19, 1028, 1978).** Für ein Haar-zufälliges
$`\psi\in\mathbb C^a\otimes\mathbb C^b`$ und $`\rho_A=\mathrm{Tr}_2|\psi\rangle\langle\psi|`$ gilt
```math
\mathbb E\mathrm{Tr}\,\rho_A^2=\frac{a+b}{ab+1}.
```

*Beweis.* Sei $`G_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,g\rangle`$, eine Ginibre-Matrix; die Matrix
von $`\psi`$ ist $`\Psi=G/\|g\|`$, und nach Proposition 7.2(a) ist
$`(\rho_A)_{\alpha\alpha'}=\sum_\beta\Psi_{\alpha\beta}\bar\Psi_{\alpha'\beta}`$, d. h.
$`\rho_A=\Psi\Psi^*=GG^*/\|g\|^2`$. Folglich ist
$`\mathrm{Tr}(GG^*)^2=\|g\|^4\mathrm{Tr}\,\rho_A^2`$, wobei $`\mathrm{Tr}\,\rho_A^2\in[0,1]`$ eine Funktion von
$`\psi`$ ist. Nach Lemma 7.9 ist $`\mathbb E\mathrm{Tr}(GG^*)^2=\mathbb E\|g\|^4\,\mathbb E\mathrm{Tr}\,\rho_A^2`$, und nach
Lemma 7.10 ist $`\mathbb E\mathrm{Tr}\,\rho_A^2=ab(a+b)/(ab(ab+1))`$. $`\square`$

Kapitel 6 erhält denselben Mittelwert aus den vierten Momenten eines Einheitsvektors mit unitär
invarianter Verteilung (Lemma 6.14), als Ausgangsreinheit eines zufälligen Kanals (Satz 6.15 und
die Bemerkung danach).

<a id="r-7-12"></a>

**Korollar 7.12.** $`\log\frac{ab+1}{a+b}\le\mathbb E\,S(\rho_A)\le\log\min(a,b)`$.

*Beweis.* Nach Korollar 3.6 ist $`S(\rho_A)\le\log a`$, und nach Lemma 7.7 ist
$`S(\rho_A)=S(\mathrm{Tr}_1|\psi\rangle\langle\psi|)\le\log b`$. Die endliche Jensensche Ungleichung für den
konkaven $`\log`$ (Lösung von Aufgabe 2 in Kapitel 3) mit den Gewichten $`\lambda_i>0`$, den von null
verschiedenen Eigenwerten von $`\rho_A`$, liefert $`\sum_i\lambda_i\log\lambda_i\le\log\sum_i\lambda_i^2`$,
d. h. $`S(\rho_A)\ge-\log X`$ mit $`X=\mathrm{Tr}\,\rho_A^2\in(0,1]`$. Sei $`m=\mathbb EX=\frac{a+b}{ab+1}`$
(Satz 7.11). Lemma 3.4 für die konvexe Funktion $`-\log`$ auf $`(0,\infty)`$ im Punkt $`m`$ liefert
$`-\log X\ge-\log m-(X-m)/m`$; beide Seiten sind beschränkt ($`X\ge1/a`$ nach Cauchy-Schwarz), und
Bilden der Erwartungswerte ergibt $`\mathbb E[-\log X]\ge-\log m`$. $`\square`$

Für festes $`a`$ und $`b\to\infty`$ treffen sich die Schranken: Die mittlere Entropie eines kleinen
Teils eines zufälligen reinen Zustands strebt gegen ihr Maximum $`\log a`$. Der exakte Mittelwert ist
die Pagesche Formel, die in Anhang C bewiesen wird.

<a id="r-7-13"></a>

**Satz 7.13 (Pagesche Formel).** Für $`a\le b`$ gilt
```math
\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}.
```
Außerdem ist $`\log a-\frac{a-1}{ab+1}-\frac{a-1}{2b}\le\mathbb E\,S(\rho_A)\le\log a-\frac{a-1}{2b}`$, also
$`\big|\mathbb E\,S(\rho_A)-\log a+\frac a{2b}\big|\le\frac1{2b}`$; für zwei Qubits ist $`\mathbb E\,S(\rho_A)=\frac13`$.

*Beweis.* Dies ist Satz C.17 zusammen mit Korollar C.18. Der Beweis in Anhang C geht von
$`\rho_A=GG^*/\mathrm{Tr}\,GG^*`$ aus (Beweis von Satz 7.11), berechnet die gemeinsame Dichte der
Eigenwerte der Wishart-Matrix $`GG^*`$ (Satz C.13) und von $`\rho_A`$ (Satz C.14) und integriert
$`x\log x`$ gegen die Laguerre-Niveaudichte (Lemmata C.15 und C.16). $`\square`$

*Rechnung:* Die Sätze 7.11 und 7.13 und Korollar 7.12 an $`10^4`$ zufälligen reinen Zuständen von $`\mathbb C^3\otimes\mathbb C^6`$ — [Tutorial, Abschnitt 7.5](../../tutorial/de/07_tensor_products.md#s-7-5).

Die Formel wurde von D. N. Page vermutet (Phys. Rev. Lett. 71, 1291, 1993) und von S. K. Foong und
S. Kanno (Phys. Rev. Lett. 72, 1148, 1994), J. Sánchez-Ruiz (Phys. Rev. E 52, 5653, 1995) und
S. Sen (Phys. Rev. Lett. 77, 1, 1996) bewiesen.

Für $`a>b`$ vertauschen $`a`$ und $`b`$ ihre Rollen. Denn die unitäre Abbildung
$`e_\alpha\otimes e_\beta\mapsto e_\beta\otimes e_\alpha`$ von $`\mathbb C^a\otimes\mathbb C^b`$ auf $`\mathbb C^b\otimes\mathbb C^a`$
permutiert die Koordinaten von $`g`$, bildet also einen Haar-zufälligen reinen Zustand $`\psi`$ auf
einen Haar-zufälligen reinen Zustand $`\psi'`$ von $`\mathbb C^b\otimes\mathbb C^a`$ ab; dessen Matrix ist
$`\Psi^T`$, also ist $`\mathrm{Tr}_2|\psi'\rangle\langle\psi'|=\Psi^T\bar\Psi=\mathrm{Tr}_1|\psi\rangle\langle\psi|`$
(Proposition 7.2(a)), dessen Entropie nach Lemma 7.7 gleich $`S(\rho_A)`$ ist. Mit $`a'=\min(a,b)`$ und
$`b'=\max(a,b)`$ besagt Satz 7.13 daher, dass die mittlere Verschränkungsentropie ihr Maximum
$`\log a'`$ um mindestens $`\frac{a'-1}{2b'}`$ und höchstens $`\frac{a'-1}{2b'}+\frac{a'-1}{a'b'+1}`$ unterschreitet,
was weniger als ein Nat ist: Für $`b'=1`$ ist es $`0`$, und sonst ist $`\frac{a'-1}{2b'}\le\frac12`$ und
$`\frac{a'-1}{a'b'+1}<\frac1{b'}\le\frac12`$. Für festes $`ab`$ liegt der Mittelwert also, als Funktion von
$`a`$, weniger als ein Nat unter $`\log\min(a,b)=\min(\log a,\log(ab)-\log a)`$.

Die Verteilung von $`\rho_A`$ ist das **induzierte Maß** (der Name stammt von K. Życzkowski und
H.-J. Sommers, J. Phys. A 34, 7111, 2001; Satz C.14 gibt die Dichte seiner Eigenwerte an).
`states.random_density_matrix(alg, rank=r)` bildet $`G_cG_c^*`$ mit unabhängigen
$`k_c\times r`$-Ginibre-Blöcken $`G_c`$ und teilt durch $`\sum_c\mathrm{Tr}\,G_cG_c^*`$. Für einen einzigen
Block ist dies die obige Verteilung von $`\rho_A`$ mit $`b=r`$ (Beweis von Satz 7.11). Für mehrere
Blöcke ist es die Verteilung der Dichte der Einschränkung des Vektorzustands eines
Haar-zufälligen reinen Zustands $`\psi=g/\|g\|`$ von $`H_1\otimes\mathbb C^r`$ auf $`\mathcal M_1`$. Denn sei
$`G_{i\beta}=\langle e_i\otimes e_\beta,g\rangle`$, eine $`D\times r`$-Ginibre-Matrix ($`D=\sum_ck_c`$), und
seien $`G_c`$ ihre Zeilen im Block $`c`$: unabhängige $`k_c\times r`$-Ginibre-Matrizen. Man nehme
$`\mathcal M_2=M_r`$ (ein Sektor) und eine konstante Fusionsregel; dann gibt es eine einzige
Gesamtladung, $`\Pi=1`$, und $`\mathcal M_1\otimes_f\mathcal M_2=B(H_1\otimes\mathbb C^r)`$. Nach Satz 7.6 hat die
Einschränkung die Dichte $`E_1(\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|)`$, und
$`\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|=GG^*/\|g\|^2`$ (Proposition 7.2(a), wie im Beweis von
Satz 7.11), dessen Diagonalblöcke $`G_cG_c^*/\|g\|^2`$ sind, mit $`\|g\|^2=\sum_c\mathrm{Tr}\,G_cG_c^*`$.

*In der Bibliothek:* `states.random_density_matrix(alg, rank=r)` — [Tutorial, Abschnitt 7.5](../../tutorial/de/07_tensor_products.md#s-7-5).

<a id="r-7-14"></a>

**Proposition 7.14.** Für $`\rho`$, gezogen mit `random_density_matrix(alg, rank=r)` auf
$`\bigoplus_cM_{k_c}`$, und $`D=\sum_ck_c`$ gilt
```math
\mathbb E\mathrm{Tr}\,\rho^2=\frac{\sum_ck_c(k_c+r)}{D(Dr+1)} .
```

*Beweis.* $`\mathrm{Tr}\,\rho^2=\sum_c\mathrm{Tr}(G_cG_c^*)^2/\|g\|^4`$, wobei $`g\in\mathbb C^{Dr}`$ die
Gesamtheit aller Einträge ist, ist eine Funktion von $`g/\|g\|`$ mit Werten in $`[0,1]`$. Wie in
Satz 7.11 liefert Lemma 7.9 $`\sum_c\mathbb E\mathrm{Tr}(G_cG_c^*)^2=\mathbb E\|g\|^4\,\mathbb E\mathrm{Tr}\,\rho^2`$, und
Lemma 7.10 liefert $`\mathbb E\mathrm{Tr}(G_cG_c^*)^2=k_cr(k_c+r)`$ und $`\mathbb E\|g\|^4=Dr(Dr+1)`$, also
$`\mathbb E\mathrm{Tr}\,\rho^2=\sum_ck_cr(k_c+r)/(Dr(Dr+1))`$. $`\square`$

*Rechnung:* Proposition 7.14 für $`M_1\oplus M_2\oplus M_3`$, $`r=2`$ — [Tutorial, Abschnitt 7.5](../../tutorial/de/07_tensor_products.md#s-7-5).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Die mittlere Verschränkungsentropie
> eines zufälligen reinen Zustands, aufgetragen gegen die Größe des Teils bei fester
> Gesamtdimension (die nach Satz 7.13 beschriebene Kurve), heißt *Page-Kurve*. In der Physik
> beschreibt man mit ihr, wie Information ein Schwarzes Loch in der Hawking-Strahlung verlässt und
> warum die Teile eines abgeschlossenen Systems thermisch aussehen; das sind Deutungen, keine
> Ergebnisse dieses Kurses.

<a id="s-7-6"></a>

## 7.6 Womit man vergleichen sollte

*In der Bibliothek:* Vergleich mit anderer Software — [Tutorial, Abschnitt 7.6](../../tutorial/de/07_tensor_products.md#s-7-6).

<a id="x-7"></a>

## Aufgaben

1. Man berechne die duale Abbildung von $`\mathrm{Tr}_2`$ bezüglich $`\tau_{\mathrm{vN}}`$ auf
   $`\mathcal M_1\otimes\mathcal M_2`$ und auf $`\mathcal M_1`$.
2. (Subadditivität.) Sei $`\rho`$ ein Zustand auf $`\mathcal M_1\otimes\mathcal M_2`$ mit invertierbaren
   reduzierten Zuständen $`\rho_1,\rho_2`$. Man zeige $`S(\rho)\le S(\rho_1)+S(\rho_2)`$, mit Gleichheit
   genau dann, wenn $`\rho=\rho_1\otimes\rho_2`$.
3. Sei `site` die Algebra $`\mathbb C\oplus\mathbb C`$ mit den Ladungen $`0,1`$ (ein leerer oder besetzter
   Platz), mit sich selbst durch Addition fusioniert. Man berechne $`E(\psi)`$ und seine beiden Anteile
   für $`\psi=(|01\rangle+|10\rangle)/\sqrt2`$, ein auf zwei Plätze verteiltes Teilchen, und
   beschreibe die Einschränkung des Vektorzustands von $`(|00\rangle+|11\rangle)/\sqrt2`$ auf die
   fusionierte Algebra.

*Rechnung:* Die drei Aufgaben numerisch — [Tutorial, Aufgaben](../../tutorial/de/07_tensor_products.md#x-7).

<details>
<summary>Lösungen</summary>

1. Die Gewichte sind $`1/(C_1k_c)`$ auf $`\mathcal M_1`$ und $`1/(C_1C_2k_cl_d)`$ auf dem Produkt. Für
   $`A\in\mathcal M_1`$ gilt nach Proposition 7.2(a)
   $`\tau_{\mathrm{vN}}((\mathrm{Tr}_2\,X)A)=\sum_c\frac1{C_1k_c}\sum_d\mathrm{Tr}(X_{cd}(A_c\otimes1))    =\sum_{c,d}\frac1{C_1C_2k_cl_d}\mathrm{Tr}\big(X_{cd}\,C_2l_d(A_c\otimes1)\big)`$, also ist die duale
   Abbildung $`A\mapsto\bigoplus_{c,d}C_2l_d\,A_c\otimes1_{l_d}=A\otimes Z`$ mit dem zentralen Element
   $`Z=\bigoplus_dC_2l_d1_{l_d}`$ von $`\mathcal M_2`$. Sie ist nur für $`Z=1`$ unital, d. h. für
   $`\mathcal M_2=\mathbb C`$, und im Beispiel nicht ($`Z=4\cdot1_2\oplus6\cdot1_3`$). Für `Tr_norm` liefert
   dieselbe Rechnung $`A\otimes\bigoplus_dl_d1_{l_d}`$, unital nur, wenn alle $`l_d=1`$ sind. Die
   ungewichteten Spuren machen $`\mathrm{Tr}_2`$ stets zur dualen Abbildung von $`A\mapsto A\otimes1`$;
   allgemein leisten das $`\mathrm{tr}_w`$ auf $`\mathcal M_1`$ und die Gewichte $`w_{cd}=w_c`$ auf dem
   Produkt, nach derselben Rechnung.

2. $`\rho_1\otimes\rho_2`$ ist ein Zustand (Satz 7.1(c)), invertierbar, mit Eigenvektoren
   $`u_i\otimes v_j`$ und Eigenwerten $`\lambda_i\mu_j`$, also gilt
   $`\log(\rho_1\otimes\rho_2)=\log\rho_1\otimes1+1\otimes\log\rho_2`$ (Satz 1.5(a)). Nach der Definition
   der partiellen Spur ist $`\mathrm{Tr}\,\rho(\log\rho_1\otimes1)=\mathrm{Tr}\,\rho_1\log\rho_1`$, und ebenso
   für den zweiten Faktor. Folglich ist $`S(\rho\|\rho_1\otimes\rho_2)=S(\rho_1)+S(\rho_2)-S(\rho)`$, und
   dies ist $`\ge0`$, mit Gleichheit genau dann, wenn $`\rho=\rho_1\otimes\rho_2`$ (Korollar 3.6). Diese
   Differenz ist die Transinformation (wechselseitige Information).

3. Die Sektoren sind $`N=0,1,2`$ mit den Größen $`1,2,1`$. $`\psi`$ liegt in $`N=1`$, mit
   $`\psi_0=|01\rangle/\sqrt2`$, $`\psi_1=|10\rangle/\sqrt2`$, also ist $`p=(\frac12,\frac12)`$, und da die
   Blöcke eindimensional sind, ist $`E(\hat\psi_c)=0`$: Nach Proposition 7.8 ist
   $`E(\psi)=H(p)=\log2`$ vollständig Teilchenzahlverschränkung, und der
   Konfigurationsanteil ist $`0`$. Der Vektor $`(|00\rangle+|11\rangle)/\sqrt2`$ überlagert $`N=0`$ und
   $`N=2`$; auf der blockdiagonalen fusionierten Algebra ist sein Vektorzustand die Mischung
   $`\frac12|00\rangle\langle00|+\frac12|11\rangle\langle11|`$ von Produktzuständen (Kapitel 3,
   Aufgabe 3).

</details>

**Weiter:** [Kapitel 8](08_dynamics.md) behandelt die Dynamik: die Schrödinger- und die
Lindblad-Gleichung in einer Algebra mit Sektoren, Sprünge zwischen Sektoren und Postselektion.
