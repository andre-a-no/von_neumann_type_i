[English](../04_symmetries.md) · [Español](../es/04_symmetries.md) · [Français](../fr/04_symmetries.md) · **Deutsch** · [中文](../zh/04_symmetries.md) · [日本語](../ja/04_symmetries.md) · [Русский](../ru/04_symmetries.md)

# 4. Symmetrien und Darstellungen

Kapitel 2 hat gezeigt, dass eine Symmetriegruppe eine von-Neumann-Algebra hervorbringt, den Spann der
Symmetrieoperatoren, sowie deren Kommutante, die symmetrischen Observablen (Proposition 2.7), und hat ein
Beispiel von Hand zerlegt. Dieses Kapitel liefert die zugrunde liegende Darstellungstheorie, macht daraus einen
Algorithmus, den die Bibliothek auf jede Liste von Matrizen anwendet, und wendet sie auf die beiden Symmetrien von
Spinsystemen an: $`U(1)`$ (eine Teilchenzahl oder eine Magnetisierung) und $`SU(2)`$ (die Drehungen der Spins).

**Sie lernen** die vollständige Reduzibilität und das Schursche Lemma kennen, wie die Zerlegung in irreduzible
Darstellungen die Kommutante und die erzeugte Algebra bestimmt, wie Charaktere Vielfachheiten zählen, warum ein
zufälliges Element der Kommutante die Zerlegung numerisch findet, was Ladungssektoren sind, die irreduziblen
Spindarstellungen und wie $`L`$ Spins $`\frac12`$ zerfallen, warum $`SU(2)`$-invariante Hamiltonoperatoren viel
kleinere Blöcke haben als $`U(1)`$-invariante, und wie die komplexe Konjugation die Sektoren einer von reellen Matrizen
erzeugten Algebra in reellen, komplexen und quaternionischen Typ einteilt. Die numerischen Kontrollen dieses Kapitels stehen im
[Tutorial 4](../../tutorial/de/04_symmetries.md).

Im Folgenden ist $`S'`$ die Kommutante von $`S\subseteq M_n(\mathbb C)`$ (Kapitel 2). Ist $`S^*=S`$, so ist $`S''`$ die
**von $`S`$ erzeugte von-Neumann-Algebra**: eine $`*`$-Algebra, die $`S`$ enthält, und jede $`*`$-Algebra
$`\mathcal B\supseteq S`$ erfüllt $`\mathcal B=\mathcal B''\supseteq S''`$ (Satz 2.2). $`\mathcal A'_h=\{K\in\mathcal A':K=K^*\}`$ ist
ein reeller Vektorraum. $`\mathbb C^k\otimes\mathbb C^m`$ trägt die Kronecker-Reihenfolge aus Kapitel 7
($`e_a\otimes f_j`$ ist der Basisvektor $`am+j`$, gezählt ab $`0`$), sodass $`1_k\otimes x`$ blockdiagonal mit $`k`$ Kopien
von $`x`$ ist. $`|G|`$ ist die Ordnung einer endlichen Gruppe, und $`\hbar=1`$.

<a id="s-4-1"></a>

## 4.1 Darstellungen und das Schursche Lemma

**Definition.** Sei $`X`$ eine Menge mit einer Involution $`x\mapsto x^*`$ ($`x^{**}=x`$). Eine
**Darstellung** von $`X`$ auf $`\mathbb C^k`$ ist eine Abbildung $`\rho:X\to M_k(\mathbb C)`$ mit
$`\rho(x^*)=\rho(x)^*`$. Ein Unterraum $`K`$ heißt **invariant**, wenn $`\rho(x)K\subseteq K`$ für alle $`x`$ gilt; $`\rho`$
heißt **irreduzibel**, wenn $`k\ge1`$ ist und $`0,\mathbb C^k`$ die einzigen invarianten Unterräume sind. Ein
**Verflechtungsoperator** von $`\rho`$ in eine Darstellung $`\rho'`$ auf $`\mathbb C^{k'}`$ ist eine
$`k'\times k`$-Matrix $`T`$ mit $`T\rho(x)=\rho'(x)T`$ für alle $`x`$; sie bilden einen Vektorraum
$`\mathrm{Hom}(\rho,\rho')`$, und es ist $`\mathrm{Hom}(\rho,\rho)=\rho(X)'`$. $`\rho`$ und $`\rho'`$ heißen
**äquivalent**, $`\rho\cong\rho'`$, wenn es einen invertierbaren Verflechtungsoperator gibt.

Drei Fälle kommen vor: eine Gruppe $`G`$ mit $`g^*=g^{-1}`$, wobei wir verlangen, dass $`\rho`$ ein Homomorphismus
$`G\to U(k)`$ ist (eine **unitäre Darstellung**; dann ist $`\rho(g^{-1})=\rho(g)^{-1}=\rho(g)^*`$); eine $`*`$-Algebra mit
ihrer Adjungierten, wobei $`\rho`$ ein unitaler $`*`$-Homomorphismus ist; und Bezeichnungen von Erzeugern mit $`x^*=x`$,
etwa $`\{x,y,z\}`$, dargestellt durch Tupel hermitescher Matrizen wie die Spinmatrizen (Abschnitt 4.5). Ist $`K`$
invariant und bilden die Spalten von $`Q\in M_{k\times d}`$ eine Orthonormalbasis von $`K`$, so gilt
$`\rho(x)Q=QQ^*\rho(x)Q`$, also ist $`\rho_K(x)=Q^*\rho(x)Q`$ wieder eine Darstellung (ein Homomorphismus, wenn $`\rho`$
einer ist), die **Einschränkung** von $`\rho`$ auf $`K`$.

<a id="r-4-1"></a>

**Lemma 4.1 (vollständige Reduzibilität).** Sei $`\rho`$ eine Darstellung von $`X`$ auf $`\mathbb C^n`$.
(a) Ist $`K`$ invariant, so auch $`K^\perp`$. (b) $`\mathbb C^n`$ ist eine orthogonale direkte Summe invarianter
Unterräume, auf denen $`\rho`$ jeweils irreduzibel ist.

*Beweis.* (a) Für $`y\in K^\perp`$ und $`z\in K`$ ist $`\langle z,\rho(x)y\rangle=\langle\rho(x^*)z,y\rangle=0`$, da
$`\rho(x^*)z\in K`$. (b) Induktion nach $`n`$. Ein von null verschiedener invarianter Unterraum $`K_1`$ minimaler
Dimension ist irreduzibel; $`K_1^\perp`$ ist nach (a) invariant, und die invarianten Unterräume der Einschränkung auf
$`K_1^\perp`$ sind invariante Unterräume von $`\rho`$, sodass die Induktionsvoraussetzung auf sie anwendbar ist. $`\square`$

Für eine endliche Gruppe kostet Unitarität nichts: Sind die $`\rho(g)`$ nur invertierbar, so erfüllt
$`P=\sum_g\rho(g)^*\rho(g)`$ die Ungleichung $`P\ge\rho(e)^*\rho(e)=1`$, also ist $`P>0`$, und
$`\rho(h)^*P\rho(h)=\sum_g\rho(gh)^*\rho(gh)=P`$ (man substituiere $`g\to gh^{-1}`$); daher ist
$`(P^{1/2}\rho(h)P^{-1/2})^*(P^{1/2}\rho(h)P^{-1/2})=P^{-1/2}\rho(h)^*P\rho(h)P^{-1/2}=1`$, und
$`h\mapsto P^{1/2}\rho(h)P^{-1/2}`$ ist eine zu $`\rho`$ äquivalente unitäre Darstellung (Weylscher Trick).

<a id="r-4-2"></a>

**Lemma 4.2 (Schursches Lemma).** Seien $`\rho,\rho'`$ irreduzible Darstellungen von $`X`$ auf
$`\mathbb C^k,\mathbb C^{k'}`$ und $`T\in\mathrm{Hom}(\rho,\rho')`$.
(a) $`T=0`$, oder $`T`$ ist invertierbar. (b) $`\mathrm{Hom}(\rho,\rho)=\mathbb C1`$.
(c) $`T^*\in\mathrm{Hom}(\rho',\rho)`$, und ist $`T\neq0`$, so gilt $`T=\lambda U`$ mit $`\lambda>0`$ und unitärem $`U`$.
Folglich ist $`\dim\mathrm{Hom}(\rho,\rho')`$ gleich $`1`$, falls $`\rho\cong\rho'`$, und $`0`$ sonst, und äquivalente
irreduzible Darstellungen sind unitär äquivalent.

*Beweis.* (a) und (b) wurden für Gruppen nach Proposition 2.7 bewiesen, und der Beweis verwendet nur, dass
$`\ker T`$ invariant unter $`\rho`$ ist ($`T\rho(x)v=\rho'(x)Tv=0`$, falls $`Tv=0`$) und das Bild von $`T`$ invariant unter
$`\rho'`$: also $`\ker T\in\{0,\mathbb C^k\}`$ und $`\mathrm{ran}\,T\in\{0,\mathbb C^{k'}\}`$; und für (b) ist
$`T-\lambda1`$ mit einem Eigenwert $`\lambda`$ von $`T`$ ein nicht invertierbarer Verflechtungsoperator, also $`0`$.
(c) Die Adjungierte von $`T\rho(x^*)=\rho'(x^*)T`$ ist $`\rho(x)T^*=T^*\rho'(x)`$. Dann ist
$`T^*T\in\mathrm{Hom}(\rho,\rho)`$, also $`T^*T=\lambda^21`$ mit $`\lambda>0`$, falls $`T\neq0`$ ($`T^*T\ge0`$ ist von null
verschieden), und $`U=T/\lambda`$ ist nach (a) quadratisch mit $`U^*U=1`$. Sind schließlich
$`T_1,T_2\in\mathrm{Hom}(\rho,\rho')`$ von null verschieden, so ist
$`T_1^{-1}T_2\in\mathrm{Hom}(\rho,\rho)=\mathbb C1`$. $`\square`$

Vertauschen alle $`\rho(x)`$ miteinander (eine abelsche Gruppe), so ist ein irreduzibles $`\rho`$ eindimensional: Jedes
$`\rho(x)`$ liegt in $`\mathrm{Hom}(\rho,\rho)=\mathbb C1`$, also ist jeder Unterraum invariant.

<a id="s-4-2"></a>

## 4.2 Die isotypische Zerlegung und die Kommutante

<a id="r-4-3"></a>

**Satz 4.3 (isotypische Zerlegung).** Sei $`\rho`$ eine Darstellung von $`X`$ auf $`\mathbb C^n`$.
Es gibt irreduzible, paarweise inäquivalente Darstellungen $`\rho_1,\dots,\rho_C`$ von $`X`$ auf
$`\mathbb C^{k_c}`$, ganze Zahlen $`m_c\ge1`$ mit $`\sum_ck_cm_c=n`$ und eine unitäre Abbildung
$`V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n`$ mit
$`V^*\rho(x)V=\bigoplus_c\rho_c(x)\otimes1_{m_c}`$ für alle $`x`$. Für jede solche Zerlegung gilt
```math
V^*\rho(X)'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),\qquad
V^*\rho(X)''V=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},
```
und $`m_c=\dim\mathrm{Hom}(\rho_c,\rho)`$; für eine irreduzible Darstellung $`\sigma`$ von $`X`$ ist
$`\mathrm{Hom}(\sigma,\rho)\neq0`$ genau dann, wenn $`\sigma\cong\rho_c`$ für genau ein $`c`$ gilt.

*Beweis.* *Existenz.* Lemma 4.1 liefert paarweise orthogonale invariante Unterräume $`W_1,\dots,W_M`$ mit Summe
$`\mathbb C^n`$ und irreduziblen Einschränkungen $`\sigma_i(x)=Q_i^*\rho(x)Q_i`$. Man fasse die $`W_i`$ zu Klassen
äquivalenter Einschränkungen zusammen und wähle in jeder Klasse $`c`$ eine Referenz $`\rho_c=\sigma_{i_0}`$. Für jedes
andere $`i`$ der Klasse liefert Lemma 4.2(c) ein unitäres $`U_i`$ mit $`\sigma_i(x)=U_i\rho_c(x)U_i^*`$; ersetzt man $`Q_i`$
durch $`Q_iU_i`$, so wird $`Q_i^*\rho(x)Q_i=\rho_c(x)`$. Man nummeriere die Mitglieder der Klasse $`c`$ mit
$`j=0,\dots,m_c-1`$, schreibe $`Q_{c,j}`$ für ihre Basen und setze $`V(e_a\otimes f_j)=Q_{c,j}e_a`$ auf dem $`c`$-ten
Summanden. Da die $`W_i`$ invariant und orthogonal sind, gilt $`Q_{c,j}^*\rho(x)Q_{d,l}=\delta_{cd}\delta_{jl}\rho_c(x)`$,
und das ist die behauptete Form.

*Kommutante.* Sei nun $`V`$ irgendeine Zerlegung wie in der Aussage, und sei $`Q_{c,j}`$ die $`n\times k_c`$-Isometrie
$`e_a\mapsto V(e_a\otimes f_j)`$ (im $`c`$-ten Summanden); für das soeben konstruierte $`V`$ sind dies die obigen $`Q_{c,j}`$. Die
Identität $`V^*\rho(x)V=\bigoplus_c\rho_c(x)\otimes1_{m_c}`$ besagt $`\rho(x)Q_{d,l}=Q_{d,l}\rho_d(x)`$, und ihre Adjungierte,
angewandt auf $`x^*`$, liefert $`Q_{c,j}^*\rho(x)=\rho_c(x)Q_{c,j}^*`$. Da $`V`$ unitär ist, vertauscht $`Y\in M_n`$ genau dann mit allen $`\rho(x)`$,
wenn jeder Block $`Y_{cj,dl}=Q_{c,j}^*YQ_{d,l}`$ die Gleichung $`Y_{cj,dl}\rho_d(x)=\rho_c(x)Y_{cj,dl}`$ erfüllt, d. h. in
$`\mathrm{Hom}(\rho_d,\rho_c)`$ liegt, und dieser Raum ist $`0`$ für $`c\neq d`$ und $`\mathbb C1`$ für $`c=d`$ (Lemma 4.2).
Also ist $`Y\in\rho(X)'`$ genau dann, wenn $`Y_{cj,dl}=\delta_{cd}(y_c)_{jl}1_{k_c}`$, d. h.
$`V^*YV=\bigoplus_c1_{k_c}\otimes y_c`$.

*Bikommutante.* Vertauscht $`Z`$ mit $`\rho(X)'`$, so vertauscht es mit den Projektionen auf die Summanden (man nehme
$`y_c=1`$ und $`y_d=0`$ für $`d\neq c`$), also ist $`V^*ZV=\bigoplus_cZ_c`$, und $`Z_c`$ vertauscht mit
$`1_{k_c}\otimes M_{m_c}`$; nach Lemma 2.5, Schritt 5 (mit vertauschten Tensorfaktoren), ist
$`Z_c\in M_{k_c}\otimes1_{m_c}`$. Umgekehrt vertauschen diese Operatoren mit $`\rho(X)'`$.

*Vielfachheiten.* Für $`T\in\mathrm{Hom}(\sigma,\rho)`$ mit irreduziblem $`\sigma`$ zeigt dasselbe Argument
$`Q_{c,j}^*T\in\mathrm{Hom}(\sigma,\rho_c)`$, und dieser Raum ist $`\mathbb Ct_{c}`$ mit invertierbarem $`t_c`$, falls
$`\sigma\cong\rho_c`$, und $`0`$ sonst; ferner ist $`T=\sum_{c,j}Q_{c,j}Q_{c,j}^*T`$. Da die $`\rho_c`$ paarweise inäquivalent
sind, ist $`\sigma`$ zu höchstens einer von ihnen äquivalent, etwa zu $`\rho_c`$ (andernfalls ist $`T=0`$), und dann ist
$`T=\sum_ja_jQ_{c,j}t_c`$; umgekehrt liegt jedes $`Q_{c,j}t_c`$ in $`\mathrm{Hom}(\sigma,\rho)`$, und diese $`m_c`$
Abbildungen haben orthogonale Bilder, also ist $`\dim\mathrm{Hom}(\sigma,\rho)=m_c`$. $`\square`$

*In der Bibliothek:* Eine Zerlegung wie in Satz 4.3 berechnen `decompose.from_group` und `decompose.from_generators` — [Tutorial, Abschnitt 4.3](../../tutorial/de/04_symmetries.md#s-4-3).

Für eine Gruppe ist $`\rho(G)''`$ der Spann von $`\rho(G)`$ (Abschnitt 2.4), und Satz 4.3 ist Proposition 2.7, in der
anderen Richtung gelesen: Dort brachte der Struktursatz die irreduziblen Darstellungen hervor, hier bringen sie die
Struktur hervor; in beiden Fällen ist $`\dim\rho(X)'=\sum_cm_c^2`$ und $`\dim\rho(X)''=\sum_ck_c^2`$. Das Bild von
$`P_c=\sum_jQ_{c,j}Q_{c,j}^*`$, einer minimalen zentralen Projektion von $`\rho(X)''`$, ist die **isotypische
Komponente** von $`\rho_c`$.
Nach Satz 4.3 ist jedes hermitesche $`H\in\rho(X)'`$ von der Form $`V(\bigoplus_c1_{k_c}\otimes h_c)V^*`$ mit hermiteschen
$`h_c\in M_{m_c}`$, und $`1_{k_c}\otimes h_c`$ ist blockdiagonal mit $`k_c`$ Kopien von $`h_c`$: Das Spektrum von $`H`$ ist die
Vereinigung der Spektren der $`m_c\times m_c`$-Blöcke $`h_c`$, wobei jeder Eigenwert von $`h_c`$ $`k_c`$-mal wiederholt wird.

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Ist $`\rho(X)`$ eine Symmetrie und
> $`H\in\rho(X)'`$ der Hamiltonoperator, so sind die eben beschriebenen $`k_c`$-fach entarteten Niveaus die
> **Multipletts** der Symmetrie, und die Blöcke, die man diagonalisieren muss, haben die Größen $`m_c`$ der
> Vielfachheiten, nicht die Größen der Multipletts.

Für eine endliche Gruppe lassen sich die Vielfachheiten berechnen, ohne irgendetwas zu zerlegen.

<a id="r-4-4"></a>

**Satz 4.4 (Schursche Orthogonalitätsrelationen).** Seien $`\rho,\rho'`$ irreduzible unitäre Darstellungen einer
endlichen Gruppe $`G`$ auf $`\mathbb C^k,\mathbb C^{k'}`$. Dann gilt
```math
\frac1{|G|}\sum_{g\in G}\rho(g)_{ij}\,\overline{\rho'(g)_{i'j'}}=
\begin{cases}\delta_{ii'}\delta_{jj'}/k&\text{if }\rho=\rho',\\0&\text{if }\rho\not\cong\rho'.\end{cases}
```

*Beweis.* Für eine $`k\times k'`$-Matrix $`Y`$ sei $`T=\frac1{|G|}\sum_g\rho(g)Y\rho'(g)^*`$. Für $`h\in G`$ liefern die
Substitution $`g\to hg`$ und $`\rho'(hg)^*=\rho'(g)^*\rho'(h)^{-1}`$ die Gleichung $`\rho(h)T=T\rho'(h)`$, also ist
$`T\in\mathrm{Hom}(\rho',\rho)`$. Ist $`\rho\not\cong\rho'`$, so ist $`T=0`$; ist $`\rho=\rho'`$, so ist $`T=\lambda1`$ mit
$`k\lambda=\mathrm{Tr}\,T=\mathrm{Tr}\,Y`$ (Lemma 4.2). Für $`Y=|j\rangle\langle j'|`$ ist der Eintrag $`T_{ii'}`$ die
linke Seite, und $`\mathrm{Tr}\,Y=\delta_{jj'}`$. $`\square`$

Der **Charakter** einer Darstellung ist $`\chi_\rho(g)=\mathrm{Tr}\,\rho(g)`$; äquivalente Darstellungen haben gleiche
Charaktere, da $`\mathrm{Tr}(T\rho(g)T^{-1})=\mathrm{Tr}\,\rho(g)`$. Wir schreiben
$`\langle\varphi,\psi\rangle_G=\frac1{|G|}\sum_g\overline{\varphi(g)}\psi(g)`$.

<a id="r-4-5"></a>

**Korollar 4.5 (Charaktere).** Sei $`G`$ endlich.
(a) Für irreduzible unitäre $`\rho,\rho'`$ ist $`\langle\chi_\rho,\chi_{\rho'}\rangle_G`$ gleich $`1`$, falls $`\rho\cong\rho'`$, und $`0`$ sonst.
(b) Für eine wie in Satz 4.3 zerlegte unitäre Darstellung $`\pi`$ und ein irreduzibles $`\rho`$ ist die Vielfachheit
von $`\rho`$ in $`\pi`$
```math
m_\rho=\frac1{|G|}\sum_{g\in G}\overline{\chi_\rho(g)}\,\chi_\pi(g)
```
($`m_\rho=m_c`$, falls $`\rho\cong\rho_c`$, und $`m_\rho=0`$, falls $`\rho`$ zu keinem $`\rho_c`$ äquivalent ist).
(c) $`\dim\pi(G)'=\sum_cm_c^2=\frac1{|G|}\sum_{g\in G}|\chi_\pi(g)|^2`$. Insbesondere ist $`\pi`$ genau dann irreduzibel,
wenn diese Zahl $`1`$ ist.

*Beweis.* (a) Man ersetze $`\rho'`$ durch eine äquivalente Darstellung (gleicher Charakter), summiere Satz 4.4 über
$`i=j`$, $`i'=j'`$ und konjugiere: Das Ergebnis ist $`\sum_{i,i'}\delta_{ii'}/k=1`$, falls $`\rho=\rho'`$, und $`0`$, falls
$`\rho\not\cong\rho'`$. (b) Die Spur ist unitär invariant und $`\mathrm{Tr}(a\otimes1_m)=m\mathrm{Tr}\,a`$, also
$`\chi_\pi=\sum_cm_c\chi_{\rho_c}`$; man wende (a) an. (c) Mit derselben Entwicklung und (a) ist
$`\langle\chi_\pi,\chi_\pi\rangle_G=\sum_cm_c^2`$, und das ist $`\dim\pi(G)'`$ nach Satz 4.3. $`\square`$

*Rechnung:* Korollar 4.5 für $`S_3`$, das drei Qubits permutiert — [Tutorial, Abschnitt 4.2](../../tutorial/de/04_symmetries.md#s-4-2); Satz 4.4 für die Standarddarstellung von $`S_3`$ — [Tutorial, Abschnitt 4.3](../../tutorial/de/04_symmetries.md#s-4-3).

Ein zweiter Beweis von (c) zeigt, wie die Bibliothek Kommutanten abtastet. Der **Twirl**
$`P(Y)=\frac1{|G|}\sum_g\pi(g)Y\pi(g)^*`$ ist idempotent mit Bild $`\pi(G)'`$ (Kapitel 2, Aufgabe 3) und selbstadjungiert
bezüglich $`\langle Y,Z\rangle_{HS}=\mathrm{Tr}(Y^*Z)`$ (man substituiere $`g\to g^{-1}`$), also die orthogonale
Projektion auf $`\pi(G)'`$, und $`\dim\pi(G)'`$ ist ihre Spur. Die Abbildung $`Y\mapsto AYB`$ hat die Spur
$`\mathrm{Tr}\,A\mathrm{Tr}\,B`$ (ihre Diagonaleinträge sind $`A_{ii}B_{jj}`$), was
$`\frac1{|G|}\sum_g|\chi_\pi(g)|^2`$ ergibt.

*In der Bibliothek:* Der Twirl ist `decompose.twirl` — [Tutorial, Abschnitt 4.3](../../tutorial/de/04_symmetries.md#s-4-3); $`\dim\pi(G)'`$ ist `decompose.commutant_dimension` — [Tutorial, Abschnitt 4.2](../../tutorial/de/04_symmetries.md#s-4-2).

Wir rechnen das Beispiel aus Abschnitt 2.4 mit Charakteren nach. Eine Permutationsmatrix hat als Spur die Anzahl
der Basisvektoren, die sie festhält, und $`\pi(\sigma)`$ hält $`|b_1b_2b_3\rangle`$ genau dann fest, wenn $`b`$ auf den
Zyklen von $`\sigma`$ konstant ist (Lemma B.1(c)), also ist $`\chi_\pi(\sigma)=2^{\#\text{cycles}(\sigma)}`$. Drei unitäre
Darstellungen von $`S_3`$ liegen auf der Hand. Sei $`P_\sigma e_i=e_{\sigma(i)}`$ die Permutationsdarstellung auf $`\mathbb C^3`$
(ein Homomorphismus in die Permutationsmatrizen). Die **triviale** Darstellung ist $`\sigma\mapsto1`$. Das **Signum** ist
$`\mathrm{sgn}\,\sigma=\det P_\sigma`$, ein Homomorphismus nach $`\{\pm1\}\subseteq U(1)`$, da die Determinante
multiplikativ ist und die Determinante einer Permutationsmatrix $`\pm1`$ ist. Die Gerade $`\mathbb C(1,1,1)`$ ist invariant,
also auch $`E=(1,1,1)^\perp`$ (Lemma 4.1(a)); die **Standarddarstellung** ist die Einschränkung von $`P`$ auf $`E`$. In einer an
$`\mathbb C(1,1,1)\oplus E`$ angepassten Orthonormalbasis ist $`P_\sigma=1\oplus\mathrm{std}(\sigma)`$, also
$`\chi_{\rm std}(\sigma)=\mathrm{Tr}\,P_\sigma-1`$, die Anzahl der Fixpunkte minus eins. Die sechs Elemente von $`S_3`$
sind die Identität, drei Transpositionen und zwei $`3`$-Zyklen, mit $`3,2,1`$ Zyklen und $`3,1,0`$ Fixpunkten; eine
Transposition vertauscht zwei Spalten von $`1_3`$, ihr Signum ist also $`-1`$, und ein $`3`$-Zyklus ist ein Produkt zweier
Transpositionen ($`(1\,2\,3)=(1\,3)(1\,2)`$), sein Signum ist also $`+1`$. Das ergibt die Tabelle

| | Identität ($`1`$) | Transpositionen ($`3`$) | $`3`$-Zyklen ($`2`$) |
|---|---|---|---|
| $`\chi_{\rm triv}`$ | $`1`$ | $`1`$ | $`1`$ |
| $`\mathrm{sgn}`$ | $`1`$ | $`-1`$ | $`1`$ |
| $`\chi_{\rm std}`$ | $`2`$ | $`0`$ | $`-1`$ |
| $`\chi_\pi`$ | $`8`$ | $`4`$ | $`2`$ |

(also $`\mathrm{sgn}\,\sigma=(-1)^{3-\#\text{cycles}(\sigma)}`$; alle Charaktere sind reell). Mit klassenweise
summiertem $`\langle\varphi,\psi\rangle_G`$:
$`\langle\chi_{\rm triv},\chi_{\rm triv}\rangle_G=\frac{1+3+2}6=1`$, $`\langle\mathrm{sgn},\mathrm{sgn}\,\rangle_G=\frac{1+3+2}6=1`$ und
$`\langle\chi_{\rm std},\chi_{\rm std}\rangle_G=\frac{4+0+2}6=1`$, also sind die drei Darstellungen nach Korollar 4.5(c)
irreduzibel, und sie sind paarweise inäquivalent, da ihre Charaktere verschieden sind. Korollar 4.5(b) liefert die
Vielfachheiten in $`\pi`$:
```math
m_{\rm triv}=\tfrac{8+3\cdot4+2\cdot2}6=4,\qquad m_{\rm sgn}=\tfrac{8-3\cdot4+2\cdot2}6=0,\qquad
m_{\rm std}=\tfrac{2\cdot8+0-2\cdot2}6=2 .
```
Also sind die triviale und die Standarddarstellung äquivalent zu Konstituenten $`\rho_c`$ von $`\pi`$ (Satz 4.3) mit
$`m_c=4`$ und $`m_c=2`$; sie machen $`4\cdot1+2\cdot2=8=\sum_ck_cm_c`$ Dimensionen aus, also gibt es keinen weiteren
Konstituenten: $`k=(1,2)`$, $`m=(4,2)`$, und $`\dim\pi(S_3)'=4^2+2^2=\frac16(64+3\cdot16+2\cdot4)=20`$ (Korollar 4.5(c)).

Das ist die Zerlegung aus Abschnitt 2.4, nun gewonnen aus einer Tabelle von sechs Zahlen pro Darstellung.

*Rechnung:* Diese Tabelle, aus $`\pi`$ neu berechnet — [Tutorial, Abschnitt 4.2](../../tutorial/de/04_symmetries.md#s-4-2).

<a id="s-4-3"></a>

## 4.3 Von den Erzeugern zur Algebra: ein randomisierter Algorithmus

Satz 4.3 beschreibt die Zerlegung, berechnet sie aber nicht: Sein Beweis benötigt die irreduziblen invarianten
Unterräume. Murota, Kanno, Kojima und Kojima (*Japan J. Indust. Appl. Math.* 27 (2010) 125–160; siehe auch Maehara
und Murota, *ibid.* 263–293) haben beobachtet, dass die Eigenräume eines *zufälligen Elements der Kommutante* sie
finden. Sei $`\mathcal A`$ die von gegebenen Matrizen $`S_1,\dots,S_J`$ und ihren Adjungierten erzeugte
von-Neumann-Algebra, mit einer Zerlegung $`V^*\mathcal AV=\bigoplus_cM_{k_c}\otimes1_{m_c}`$,
$`V^*\mathcal A'V=\bigoplus_c1_{k_c}\otimes M_{m_c}`$ (Satz 2.6); nur die Analyse verwendet $`V`$.

<a id="r-4-6"></a>

**Lemma 4.6.** Sei $`K\in\mathcal A'_h`$ mit $`V^*KV=\bigoplus_c1_{k_c}\otimes\kappa_c`$. Jeder Eigenraum von $`K`$ ist
$`\mathcal A`$-invariant. Sind die $`N=\sum_cm_c`$ Eigenwerte aller $`\kappa_c`$ zusammen paarweise verschieden, so sind
die Eigenräume von $`K`$ die $`N`$ Unterräume $`V(\mathbb C^{k_c}\otimes e)`$ (im $`c`$-ten Summanden, mit einem
Eigenvektor $`e`$ von $`\kappa_c`$), und $`\mathcal A`$ wirkt auf jedem von ihnen irreduzibel, als $`M_{k_c}`$.

*Beweis.* Ist $`Kv=\lambda v`$ und $`a\in\mathcal A`$, so ist $`Kav=aKv=\lambda av`$. Man schreibe
$`\kappa_c=\sum_j\mu_{cj}|e_{cj}\rangle\langle e_{cj}|`$ mit orthonormalen $`e_{cj}`$ (Satz 1.3). Dann ist
$`V^*KV=\sum_{c,j}\mu_{cj}R_{cj}`$, wobei $`R_{cj}`$ die Projektion auf $`\mathbb C^{k_c}\otimes e_{cj}`$ ist; die $`R_{cj}`$
sind orthogonal mit Summe $`1`$ und die $`\mu_{cj}`$ verschieden, also sind dies die Spektralprojektionen. Auf
$`\mathbb C^{k_c}\otimes e_{cj}`$ wirkt das Element $`V^*aV=\bigoplus_da_d\otimes1`$ als $`a_c`$, und die $`a_c`$ füllen
$`M_{k_c}`$ aus, das außer $`0`$ und $`\mathbb C^{k_c}`$ keine invarianten Unterräume hat. $`\square`$

<a id="r-4-7"></a>

**Proposition 4.7 (Generizität).** Seien $`K`$ und $`Y`$ unabhängige zufällige Elemente von $`\mathcal A'_h`$ und
$`\mathcal A'`$, deren Verteilungen Dichten bezüglich des Lebesgue-Maßes auf diesen reellen Vektorräumen haben
(Lemma A.3). Mit Wahrscheinlichkeit eins gilt: (a) Die Voraussetzung von Lemma 4.6 ist erfüllt; (b)
$`\Pi_iY\Pi_j\neq0`$ für je zwei verschiedene Spektralprojektionen $`\Pi_i,\Pi_j`$ von $`K`$, deren Bilder im selben
Summanden liegen.

*Beweis.* (a) Die hermitesche Matrix $`\kappa=\bigoplus_c\kappa_c\in M_N`$ hängt reell-linear von $`K`$ ab. Die Gramsche
Determinante $`\delta(K)=\det[\mathrm{Tr}\,\kappa^{i+j}]_{i,j=0}^{N-1}`$ von $`1,\kappa,\dots,\kappa^{N-1}`$ bezüglich
$`\langle\cdot,\cdot\rangle_{HS}`$ ist ein reelles Polynom in den Koordinaten von $`K`$ und genau dann von null
verschieden, wenn diese Potenzen linear unabhängig sind. Da $`p(\kappa)=U\mathrm{diag}(p(\lambda_l))U^*`$ für
$`\kappa=U\mathrm{diag}(\lambda)U^*`$ gilt, annulliert ein von null verschiedenes Polynom vom Grad $`<N`$ die Matrix
$`\kappa`$ genau dann, wenn es auf $`\mathrm{spec}\,\kappa`$ verschwindet, und das ist genau dann möglich, wenn
$`\kappa`$ weniger als $`N`$ verschiedene Eigenwerte hat (man nehme $`\prod_\mu(t-\mu)`$ über diese). Also ist
$`\delta(K)\neq0`$ die Voraussetzung von Lemma 4.6, und $`\delta\not\equiv0`$: Es ist von null verschieden in
$`V(\bigoplus_c1\otimes\mathrm{diag}(\text{verschiedene reelle Zahlen}))V^*\in\mathcal A'_h`$. Daher ist
$`\mathbb P(\delta(K)=0)=0`$ nach Korollar A.7(b).

(b) Die Menge $`O=\{\delta\neq0\}\subseteq\mathcal A'_h`$ ist offen ($`\delta`$ ist stetig). Für $`K\in O`$ hat $`K`$ genau $`N`$
verschiedene Eigenwerte $`\mu_1(K)<\dots<\mu_N(K)`$ (die von $`\kappa`$, mit Vielfachheiten $`k_c`$), also sind seine
Spektralprojektionen $`\Pi_1(K),\dots,\Pi_N(K)`$, nach aufsteigenden Eigenwerten nummeriert, stetig auf $`O`$
(Proposition A.13). Man halte $`l\neq l'`$ fest, und sei $`\mathcal N_{ll'}`$ die Menge der $`(K,Y)\in O\times\mathcal A'`$ mit
$`P_c\Pi_l(K)=\Pi_l(K)`$ und $`P_c\Pi_{l'}(K)=\Pi_{l'}(K)`$ für ein $`c`$ ($`P_c`$ die zentralen Projektionen) und
$`\Pi_l(K)Y\Pi_{l'}(K)=0`$; sie ist durch stetige Gleichungen auf der Borel-Menge $`O\times\mathcal A'`$ definiert, also eine
Borel-Menge. Man halte ein $`K\in O`$ fest, für das die beiden Bilder im Summanden $`c`$ liegen. Nach Lemma 4.6
projizieren $`V^*\Pi_lV`$ und $`V^*\Pi_{l'}V`$ auf $`\mathbb C^{k_c}\otimes e`$ und $`\mathbb C^{k_c}\otimes e'`$ mit orthonormalen
Eigenvektoren $`e,e'`$ von $`\kappa_c`$, und für $`V^*YV=\bigoplus_d1\otimes y_d`$ ist $`\Pi_lY\Pi_{l'}=0`$ genau dann, wenn
$`\langle e,y_ce'\rangle=0`$. Die Abbildung $`Y\mapsto\langle e,y_ce'\rangle`$ ist reell-linear $`\mathcal A'\to\mathbb C`$ und
nicht null (ihr Wert an $`V(0\oplus\dots\oplus1\otimes|e\rangle\langle e'|\oplus\dots\oplus0)V^*`$ ist $`1`$), also ist ihr Kern
ein echter Unterraum, eine Nullmenge (Korollar A.7(c)). Also ist jeder Schnitt $`\{Y:(K,Y)\in\mathcal N_{ll'}\}`$ eine
Nullmenge (er ist leer, wenn die Bilder in verschiedenen Summanden liegen), und da $`Y`$ von $`K`$ unabhängig ist, gilt
$`\mathbb P((K,Y)\in\mathcal N_{ll'})=0`$ nach Korollar A.7(d). Das Ereignis aus (b) scheitert nur auf
$`\{\delta(K)=0\}\cup\bigcup_{l\neq l'}\{(K,Y)\in\mathcal N_{ll'}\}`$, einer endlichen Vereinigung von Ereignissen der
Wahrscheinlichkeit $`0`$. $`\square`$

<a id="r-4-8"></a>

**Lemma 4.8 (Kompressionen; Ausrichten der Kopien).** (a) Sei $`W`$ ein $`\mathcal A`$-invarianter Unterraum mit
Orthonormalbasis $`Q\in M_{n\times d}`$ und $`\sigma(a)=Q^*aQ`$ die Einschränkung. Dann ist $`\sigma(\mathcal A)`$ eine
$`*`$-Algebra, $`\{Q^*YQ:Y\in\mathcal A'\}=\sigma(\mathcal A)'`$ und $`\{Q^*KQ:K\in\mathcal A'_h\}=\sigma(\mathcal A)'_h`$. $`\mathcal A`$
wirkt genau dann irreduzibel auf $`W`$, wenn $`Q^*YQ\in\mathbb C1`$ für jedes $`Y\in\mathcal A'`$. Ist das nicht der Fall und
ist $`X`$ ein zufälliges Element von $`\mathcal A'`$ oder von $`\mathcal A'_h`$ mit Dichte, so ist mit Wahrscheinlichkeit eins
$`Q^*XQ\notin\mathbb C1`$.

(b) Seien $`W_i,W_j`$ invariante Unterräume, auf denen $`\mathcal A`$ irreduzibel
wirkt, mit Orthonormalbasen $`Q_i,Q_j`$ und Einschränkungen $`\sigma_i(a)=Q_i^*aQ_i`$, und sei $`Y\in\mathcal A'`$. Dann ist
$`T=Q_i^*YQ_j\in\mathrm{Hom}(\sigma_j,\sigma_i)`$. Ist $`T\neq0`$, so ist $`T=\lambda U`$ mit $`\lambda>0`$ und unitärem
$`U`$, und $`Q_j'=Q_jU^*`$ erfüllt $`Q_j'^*aQ_j'=\sigma_i(a)`$ für alle $`a\in\mathcal A`$. Liegen $`W_i,W_j`$ in verschiedenen
Summanden ($`W_i\subseteq P_c\mathbb C^n`$, $`W_j\subseteq P_d\mathbb C^n`$, $`c\neq d`$), so ist $`T=0`$.

*Beweis.* (a) Wie im Beweis von Satz 4.3 ($`W`$ ist invariant unter $`a`$ und $`a^*`$) gilt $`aQ=Q\sigma(a)`$ und
$`Q^*a=\sigma(a)Q^*`$ für $`a\in\mathcal A`$. Daher ist $`\sigma(ab)=Q^*aQ\sigma(b)=\sigma(a)\sigma(b)`$, $`\sigma(a^*)=\sigma(a)^*`$ und
$`\sigma(1)=Q^*Q=1_d`$, also ist $`\sigma`$ ein unitaler $`*`$-Homomorphismus und $`\sigma(\mathcal A)`$ eine $`*`$-Algebra; und für
$`Y\in\mathcal A'`$ gilt $`Q^*YQ\sigma(a)=Q^*YaQ=Q^*aYQ=\sigma(a)Q^*YQ`$. Umgekehrt sei für $`T\in\sigma(\mathcal A)'`$
$`Y=QTQ^*`$; dann ist $`aY=Q\sigma(a)TQ^*=QT\sigma(a)Q^*=Ya`$, also $`Y\in\mathcal A'`$, mit $`Q^*YQ=T`$, und $`Y`$ ist hermitesch,
wenn $`T`$ es ist. Ist $`\sigma`$ irreduzibel, so ist $`\sigma(\mathcal A)'=\mathrm{Hom}(\sigma,\sigma)=\mathbb C1`$
(Lemma 4.2(b)). Ist $`E\subseteq\mathbb C^d`$ ein $`\sigma`$-invarianter Unterraum außer $`0,\mathbb C^d`$, so auch $`E^\perp`$
(Lemma 4.1(a)), und die Orthogonalprojektion auf $`E`$ ist ein nicht skalares Element von $`\sigma(\mathcal A)'_h`$. In
diesem Fall sind $`\mathbb C1\cap\sigma(\mathcal A)'`$ und $`\mathbb R1\cap\sigma(\mathcal A)'_h`$ echte reelle Unterräume, ihre
Urbilder unter den surjektiven reell-linearen Abbildungen $`X\mapsto Q^*XQ`$ sind echte Unterräume von $`\mathcal A'`$ und
$`\mathcal A'_h`$, und diese sind Nullmengen (Korollar A.7(c)).

(b) Wie in (a) gilt $`aQ_j=Q_j\sigma_j(a)`$ und $`Q_i^*a=\sigma_i(a)Q_i^*`$, also
$`T\sigma_j(a)=Q_i^*YaQ_j=Q_i^*aYQ_j=\sigma_i(a)T`$. Lemma 4.2 liefert $`T=\lambda U`$ und
$`U\sigma_j(a)U^*=\sigma_i(a)`$, und das ist $`Q_j'^*aQ_j'`$. Schließlich ist $`P_c\in\mathcal A\cap\mathcal A'`$, $`Q_i=P_cQ_i`$
und $`P_cQ_j=0`$, also $`T=Q_i^*P_cYQ_j=Q_i^*YP_cQ_j=0`$. $`\square`$

*Rechnung:* Lemma 4.6 und Lemma 4.8(b) für ein zufälliges Element von $`\pi(S_3)'`$ — [Tutorial, Abschnitt 4.3](../../tutorial/de/04_symmetries.md#s-4-3).

Der Algorithmus (`decompose.wedderburn`, mit den Schnittstellen `from_group` und `from_generators`) lautet nun:

1. *Die Kommutante abtasten.* Für alle Elemente einer endlichen Gruppe nimmt man $`Y=P(Z)`$ und $`K=P(Z'+Z'^*)`$ mit dem
   Twirl $`P`$ und unabhängigen Ginibre-Matrizen $`Z,Z'`$. $`P`$ ist eine lineare Projektion von $`M_n`$ auf $`\mathcal A'`$ mit
   $`P(Z)^*=P(Z^*)`$, also sind $`Z\mapsto P(Z)`$ und $`Z'\mapsto P(Z'+Z'^*)`$ reell-lineare Surjektionen auf $`\mathcal A'`$ und
   auf $`\mathcal A'_h`$ (für $`H\in\mathcal A'_h`$ ist $`P(\frac H2+\frac H2)=H`$); eine Ginibre-Matrix ist ein nicht
   ausgearteter zentrierter Gauß-Vektor in $`M_n`$ (Proposition A.16), also sind $`Y`$ und $`K`$ nicht ausgeartete zentrierte
   Gauß-Vektoren in $`\mathcal A'`$ und $`\mathcal A'_h`$ und haben Dichten (Satz A.15). Für Erzeuger ist
   $`\mathcal A'=\{X:[S_j,X]=[S_j^*,X]=0\}`$ der Kern einer linearen Abbildung auf $`M_n`$ (die Bibliothek verwendet die
   hermiteschen Teile $`\frac12(S_j+S_j^*)`$, $`\frac1{2i}(S_j-S_j^*)`$, die denselben komplexen Raum aufspannen wie
   $`S_j,S_j^*`$ und daher dieselbe Kommutante haben); mit einer Basis $`(B_r)_{r\le D}`$ davon ist $`Y=\sum_rc_rB_r`$ und
   $`K=Y'+Y'^*`$ mit einer unabhängigen Kopie $`Y'`$ von $`Y`$ und unabhängigen standard-komplex-Gaußschen $`c_r`$: Bilder
   eines nicht ausgearteten Gauß-Vektors in $`\mathbb C^D`$ unter der reell-linearen Bijektion $`c\mapsto\sum_rc_rB_r`$ auf
   $`\mathcal A'`$ und der reell-linearen Surjektion $`c\mapsto X+X^*`$ ($`X=\sum_rc_rB_r`$) auf $`\mathcal A'_h`$, also wieder mit
   Dichten. Jedes zufällige Element ist eine neue Ziehung, unabhängig von den früheren; die Bibliothek normiert sie
   auf Norm eins, was keinen Eigenraum und kein Verschwinden eines Blocks ändert.
2. *Aufspalten.* Die Eigenräume von $`K`$ sind irreduzibel (Lemma 4.6, Proposition 4.7(a)); ihre Dimensionen sind
   die $`k_c`$.
3. *Gruppieren und ausrichten.* Die Blöcke $`Q_i^*YQ_j`$ zwischen Eigenräumen gleicher Dimension sind genau innerhalb
   eines Summanden von null verschieden (Proposition 4.7(b), Lemma 4.8(b)); die Klassen der Relation „durch einen von
   null verschiedenen Block verbunden“ ergeben die Sektoren, ihre Größen die $`m_c`$, und die Polarfaktoren richten die
   Kopien an einer Referenzkopie pro Klasse aus (Lemma 4.8(b); eine an einer ausgerichteten Kopie ausgerichtete Kopie
   ist an der Referenz ausgerichtet).
4. *Zusammensetzen:* $`V`$ wird in Kronecker-Reihenfolge aufgebaut, sodass $`V^*aV=\bigoplus_ca_c\otimes1_{m_c}`$.

In exakter Arithmetik, mit exakt getrennten Eigenräumen von $`K`$, liefert der Algorithmus also mit Wahrscheinlichkeit
eins eine Zerlegung wie in Satz 4.3: Nach Lemma 4.6 enthält der Summand $`c`$ genau $`m_c`$ Eigenräume der Dimension
$`k_c`$, und nach dem Ausrichten ist die Konstruktion der Existenzteil des Beweises von Satz 4.3. In Gleitkommaarithmetik
behandelt die Bibliothek Eigenwerte, die näher als eine Toleranz beieinanderliegen, als gleich, sodass ein berechneter
„Eigenraum“ $`W`$ eine Summe mehrerer Eigenräume von $`K`$ sein kann; er ist trotzdem $`\mathcal A`$-invariant (Lemma 4.6),
und die Bibliothek komprimiert ein neues zufälliges $`X\in\mathcal A'_h`$ auf ihn. Nach Lemma 4.8(a) ist $`Q^*XQ`$ ein
Vielfaches von $`1`$, wann immer $`\mathcal A`$ irreduzibel auf $`W`$ wirkt, und andernfalls mit Wahrscheinlichkeit eins
nicht; dann zerlegen seine Eigenräume, die unter $`\sigma(\mathcal A)`$ invariant sind (Lemma 4.6, angewandt auf
$`\sigma(\mathcal A)`$), $`W`$ in kleinere invariante Unterräume, und die Bibliothek wiederholt dies, bis jede Kompression
skalar ist. Was auch immer in Gleitkommaarithmetik geschieht, das Ergebnis wird am Ende zertifiziert:

<a id="r-4-9"></a>

**Proposition 4.9 (Zertifikat).** Sei $`\mathcal A`$ von $`S_1,\dots,S_J`$ und ihren Adjungierten erzeugt und
$`\mathcal B=V(\bigoplus_cM_{k_c}\otimes1_{m_c})V^*`$ mit unitärem $`V`$. Liegen alle $`S_j`$ in $`\mathcal B`$ und ist
$`\dim\mathcal A'=\sum_cm_c^2`$, so ist $`\mathcal A=\mathcal B`$.

*Beweis.* $`\mathcal B`$ ist eine $`*`$-Algebra, die alle $`S_j,S_j^*`$ enthält, also ist $`\mathcal A\subseteq\mathcal B`$ und
$`\mathcal B'\subseteq\mathcal A'`$. Satz 4.3, angewandt auf die Darstellung von $`\mathcal B`$ auf $`\mathbb C^n`$ (die
Blockabbildungen $`b\mapsto b_c`$ sind irreduzibel und inäquivalent, da die $`c`$-te zentrale Projektion von einer auf
$`1`$ und von den anderen auf $`0`$ abgebildet wird), liefert $`\dim\mathcal B'=\sum_cm_c^2=\dim\mathcal A'`$. Also ist
$`\mathcal B'=\mathcal A'`$ und $`\mathcal B=\mathcal B''=\mathcal A''=\mathcal A`$. $`\square`$

*In der Bibliothek:* Der Algorithmus, `decompose.wedderburn` (Schnittstellen `from_group`, `from_generators`), das Zertifikat aus Proposition 4.9, `Decomposition.verify()`, und ihr Aufwand — [Tutorial, Abschnitt 4.3](../../tutorial/de/04_symmetries.md#s-4-3).

<a id="s-4-4"></a>

## 4.4 $`U(1)`$ und Ladungssektoren

<a id="r-4-10"></a>

**Proposition 4.10.** Sei $`Q`$ eine hermitesche Matrix auf $`\mathbb C^n`$ mit ganzzahligen Eigenwerten $`q`$ (den
**Ladungen**), Spektralprojektionen $`P_q`$ und $`d_q=\mathrm{rank}\,P_q`$, und sei $`U_\theta=e^{i\theta Q}`$; wir
schreiben $`U(1)'=\{U_\theta:\theta\in\mathbb R\}'`$ und $`U(1)''`$ für deren Kommutante.

1. $`\theta\mapsto U_\theta`$ ist eine unitäre Darstellung von $`U(1)=\mathbb R/2\pi\mathbb Z`$, und jede stetige unitäre
   Darstellung von $`U(1)`$ auf $`\mathbb C^n`$ hat diese Form.
2. Für $`X\in M_n`$ sind äquivalent: (i) $`XU_\theta=U_\theta X`$ für alle $`\theta`$; (ii) $`[X,Q]=0`$;
   (iii) $`X=\sum_qP_qXP_q`$. Daher ist $`U(1)'\cong\bigoplus_qM_{d_q}(\mathbb C)`$ (die blockdiagonalen Operatoren) und
   $`U(1)''=\mathrm{span}\{P_q\}`$: Alle $`k_c`$ sind $`1`$, und die Vielfachheiten sind die $`d_q`$.
3. (*Erhaltene Ladung.*) Für hermitesches $`H`$ gilt $`e^{itH}Qe^{-itH}=Q`$ für alle $`t`$ genau dann, wenn $`[H,Q]=0`$,
   d. h. genau dann, wenn $`H\in U(1)'`$. Dann ist jedes Sektorgewicht $`\mathrm{Tr}(P_q\rho(t))`$ von
   $`\rho(t)=e^{-itH}\rho e^{itH}`$ konstant.

*Beweis.* (1) $`U_\theta=\sum_qe^{iq\theta}P_q`$ ist unitär, $`U_\theta U_{\theta'}=U_{\theta+\theta'}`$ (Satz 1.5(c)) und
$`U_{2\pi}=1`$. Umgekehrt ist für eine stetige Darstellung $`(U_\theta)_{\theta\ge0}`$ eine stetige Halbgruppe von
Matrizen, also $`U_\theta=e^{\theta B}`$ für $`\theta\ge0`$ nach Lemma 8.5 (bewiesen in Kapitel 8; der Beweis benutzt
nur die endliche Dimension), und für $`\theta<0`$ ist $`U_\theta=U_{-\theta}^{-1}=e^{\theta B}`$. Ableiten von
$`U_\theta^*U_\theta=1`$ bei $`0`$ (Lemma 8.1) ergibt $`B^*+B=0`$, also $`B=iQ`$ mit hermiteschem $`Q`$. Die Exponentialreihe
$`e^{i\theta Q}`$ stimmt mit dem Funktionalkalkül überein (Proposition 8.2(1) mit $`H=-Q`$), also erzwingt
$`U_{2\pi}=\sum_qe^{2\pi iq}P_q=1`$ die Gleichung $`e^{2\pi iq}=1`$, $`q\in\mathbb Z`$, für jeden Eigenwert $`q`$ von $`Q`$. (2) (i)$`\Rightarrow`$(ii): Ableiten bei $`\theta=0`$.
(ii)$`\Rightarrow`$(iii): $`X`$ vertauscht mit $`P_q=f(Q)`$ (Satz 1.5(d)), also $`P_qXP_{q'}=XP_qP_{q'}=0`$ für
$`q\neq q'`$. (iii)$`\Rightarrow`$(i): $`XU_\theta=\sum_qe^{iq\theta}P_qXP_q=U_\theta X`$. Die blockdiagonalen Operatoren
bilden $`\bigoplus_qM_{d_q}\otimes1_1`$ (in einer an die Bilder der $`P_q`$ angepassten Orthonormalbasis), deren Kommutante
$`\bigoplus_q1_{d_q}\otimes M_1=\mathrm{span}\{P_q\}`$ ist: Ein Operator, der mit ihnen vertauscht, vertauscht mit
jedem $`P_q`$, ist also blockdiagonal, und jeder Block vertauscht mit $`M_{d_q}`$, ist also ein Skalar (Lemma 2.5,
Schritt 5, mit $`m=1`$). (3) $`\frac d{dt}e^{itH}Qe^{-itH}=ie^{itH}[H,Q]e^{-itH}`$ (Proposition 8.2) verschwindet genau dann für
alle $`t`$, wenn $`[H,Q]=0`$. Dann vertauscht $`e^{itH}`$ mit $`P_q`$ (Satz 1.5(d), zweimal), und
$`\mathrm{Tr}(P_q\rho(t))=\mathrm{Tr}(e^{itH}P_qe^{-itH}\rho)=\mathrm{Tr}(P_q\rho)`$. $`\square`$

Die irreduziblen Darstellungen von $`U(1)`$ sind also die Charaktere $`\theta\mapsto e^{iq\theta}`$, $`q\in\mathbb Z`$
(eindimensional, wie für jede abelsche Gruppe), und die Vielfachheit der Ladung $`q`$ ist $`d_q`$. Ladungen addieren
sich unter Tensorprodukten: Sind $`u_k`$ und $`w_l`$ Orthonormalbasen aus Eigenvektoren von $`Q_1`$ und $`Q_2`$ mit
Eigenwerten $`\alpha_k`$ und $`\beta_l`$, so besteht die Orthonormalbasis $`u_k\otimes w_l`$ aus Eigenvektoren von
$`Q_1\otimes1+1\otimes Q_2`$ (Eigenwert $`\alpha_k+\beta_l`$) und von $`e^{i\theta Q_1}\otimes e^{i\theta Q_2}`$ (Eigenwert
$`e^{i\theta\alpha_k}e^{i\theta\beta_l}`$), also nach Definition des Funktionalkalküls
```math
e^{i\theta(Q_1\otimes1+1\otimes Q_2)}=e^{i\theta Q_1}\otimes e^{i\theta Q_2},
```
und dieser Operator wirkt auf dem Bild von $`P_{q_1}\otimes P_{q_2}`$ mit Ladung $`q_1+q_2`$. Das ist die Fusionsregel $`(q_1,q_2)\mapsto q_1+q_2`$ der
fusionierten Tensorprodukte aus Kapitel 7.

Für $`n`$ Qubits ist die Ladung der **Teilchenzahloperator** $`\hat N=\sum_jn_j`$, $`n_j=|1\rangle\langle1|`$ auf Qubit
$`j`$: $`\hat N|b\rangle=|b|\,|b\rangle`$ mit dem Hamming-Gewicht $`|b|`$. Seine Sektoren haben die Dimension
$`\binom nw`$, also bilden die Operatoren, die das Hamming-Gewicht erhalten,
```math
U(1)'\cong\bigoplus_{w=0}^nM_{\binom nw}(\mathbb C),
```
und das ist das fusionierte Produkt von $`n`$ Kopien von $`M_1\oplus M_1`$ mit den
Ladungen $`0,1`$ (Bit $`1`$ = Spin aufwärts, $`\hat N=S^z+\frac n2`$). Jeder
XXZ-Hamiltonoperator $`\sum_{\langle ij\rangle}J(s^x_is^x_j+s^y_is^y_j+\Delta s^z_is^z_j)`$ erhält $`\hat N`$: $`s^z_is^z_j`$ ist diagonal in der
Basis $`|b\rangle`$, und der Hüpfterm $`s^x_is^x_j+s^y_is^y_j=\frac12(s^+_is^-_j+s^-_is^+_j)`$ (man setze $`s^\pm=s^x\pm is^y`$
ein; die Kreuzterme heben sich weg) bildet $`|b\rangle`$ auf $`0`$ oder auf die Konfiguration ab, in der ein aufwärts
gerichteter Spin von einem der Plätze $`i,j`$ auf den anderen gewandert ist, bildet also jeden Eigenraum von $`\hat N`$
in sich ab.

*Rechnung:* Proposition 4.10 für vier Spins: die Ladungssektoren, allein aus $`\hat N`$ gefunden, und ein XXZ-Ring, sektorweise diagonalisiert — [Tutorial, Abschnitt 4.4](../../tutorial/de/04_symmetries.md#s-4-4). *In der Bibliothek:* $`U(1)'`$ für $`n`$ Qubits ist `SpinChain(n).algebra`, und der XXZ-Hamiltonoperator ist `SpinChain.xxz` — [Tutorial, Abschnitt 4.4](../../tutorial/de/04_symmetries.md#s-4-4).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). $`U(1)`$-Symmetrie ist die Erhaltung einer Teilchenzahl, einer Ladung oder einer
> Magnetisierung. Proposition 4.10(3) ist das Noether-Theorem in diesem Rahmen: Die Ladung ist genau dann
> erhalten, wenn der Hamiltonoperator symmetrisch ist, und dann verlässt die Dynamik nie einen Sektor (Satz 8.3).
> Kann zudem nichts Kohärenzen zwischen Sektoren messen, so sind die Sektoren Superauswahlsektoren (Kapitel 3).

<a id="s-4-5"></a>

## 4.5 $`SU(2)`$: Spin

**Definition.** Eine **Spindarstellung** auf $`\mathbb C^d`$ ist ein Tripel hermitescher Matrizen $`S_x,S_y,S_z`$ mit
```math
[S_x,S_y]=iS_z,\qquad[S_y,S_z]=iS_x,\qquad[S_z,S_x]=iS_y,
```
d. h. eine Darstellung der Bezeichnungen $`\{x,y,z\}`$ (Abschnitt 4.1) mit diesen Relationen. Das Grundbeispiel ist
$`s_a=\frac12\sigma_a`$ auf $`\mathbb C^2`$ (Pauli-Matrizen). Sind $`(S_a)`$ und $`(T_a)`$ Spindarstellungen, so auch
$`S_a\otimes1+1\otimes T_a`$ (die Kreuzterme vertauschen); daher ist der **Gesamtspin** $`S_a=\sum_is_a^{(i)}`$ von $`L`$
Spins $`\frac12`$ (`decompose.total_spin(L)`) eine solche. Die Drehungen $`e^{-i\theta\,\mathbf n\cdot\mathbf S}`$
($`\theta\in\mathbb R`$, $`\mathbf n\in\mathbb R^3`$ ein Einheitsvektor, Funktionalkalkül) vertauschen genau dann für alle
$`\theta,\mathbf n`$ mit $`X`$, wenn $`X`$ mit den drei $`S_a`$ vertauscht (Ableiten bei $`\theta=0`$, Proposition 8.2(1);
umgekehrt Satz 1.5(d)), also ist die Kommutante „der $`SU(2)`$-Wirkung“ gleich $`\{S_x,S_y,S_z\}'`$, und nach Satz 2.2
erzeugen die Drehungen und die $`S_a`$ dieselbe Algebra (beide Mengen sind unter Adjunktion abgeschlossen). Für $`L`$
Spins $`\frac12`$ sind die Drehungen die $`V^{\otimes L}`$, $`V\in SU(2)`$. In der Tat ist $`\mathbf n\cdot\mathbf S=\sum_ih^{(i)}`$
mit $`h=\mathbf n\cdot\mathbf s`$, und wie bei den Ladungen in Abschnitt 4.4 (Induktion nach $`L`$) ist
$`e^{-i\theta\,\mathbf n\cdot\mathbf S}=(e^{-i\theta h})^{\otimes L}`$ mit $`e^{-i\theta h}\in SU(2)`$ (seine Eigenwerte sind
$`e^{\mp i\theta/2}`$, da $`h`$ spurfrei ist mit $`h^2=\frac14`$). Umgekehrt ist ein $`V\in SU(2)`$ normal mit Eigenwerten vom
Betrag $`1`$ und Produkt $`1`$, also $`V=W\mathrm{diag}(e^{i\varphi},e^{-i\varphi})W^*`$ mit unitärem $`W`$, und $`V=e^{iA}`$
für die spurfreie hermitesche Matrix $`A=W\mathrm{diag}(\varphi,-\varphi)W^*`$. Eine spurfreie hermitesche
$`2\times2`$-Matrix ist eine reelle Linearkombination $`A=2\,\mathbf c\cdot\mathbf s`$ der Pauli-Matrizen, also ist
$`V=e^{-i\theta\,\mathbf n\cdot\mathbf s}`$ mit $`\theta=2|\mathbf c|`$, $`\mathbf n=-\mathbf c/|\mathbf c|`$ (und $`\theta=0`$,
falls $`\mathbf c=0`$). Satz B.6(a),(c) und Korollar B.7 (mit $`d=2`$) zeigen direkt, dass die $`V^{\otimes L}`$,
$`V\in SU(2)`$, die von den $`S_a`$ erzeugte Algebra aufspannen.

Man setze $`S_\pm=S_x\pm iS_y`$, sodass $`S_-=S_+^*`$, und $`\mathbf S^2=S_x^2+S_y^2+S_z^2`$. Die Relationen ergeben
```math
[S_z,S_\pm]=\pm S_\pm,\qquad[S_+,S_-]=2S_z,\qquad\mathbf S^2=S_-S_++S_z^2+S_z,
```
und $`\mathbf S^2`$ vertauscht mit jedem $`S_a`$ (z. B.
$`[\mathbf S^2,S_z]=S_x[S_x,S_z]+[S_x,S_z]S_x+S_y[S_y,S_z]+[S_y,S_z]S_y=-i\{S_x,S_y\}+i\{S_y,S_x\}=0`$).

<a id="r-4-11"></a>

**Satz 4.11 (irreduzible Spindarstellungen).** Für jedes $`s\in\{0,\frac12,1,\frac32,\dots\}`$ bilden die Matrizen
von `decompose.spin_matrices(s)` – auf $`\mathbb C^{2s+1}`$ mit der Basis $`|m\rangle`$, $`m=s,s-1,\dots,-s`$,
```math
S_z|m\rangle=m|m\rangle,\qquad S_+|m\rangle=\sqrt{s(s+1)-m(m+1)}\;|m+1\rangle,\qquad S_-=S_+^*,
```
eine irreduzible Spindarstellung (**Spin $`s`$**) mit $`\mathbf S^2=s(s+1)1`$. Jede irreduzible Spindarstellung ist zu
genau einer von ihnen unitär äquivalent.

*Beweis.* *Spin $`s`$ ist eine Spindarstellung.* $`[S_z,S_+]=S_+`$ ist klar. $`S_+S_-|m\rangle=(s(s+1)-(m-1)m)|m\rangle`$
und $`S_-S_+|m\rangle=(s(s+1)-m(m+1))|m\rangle`$, also $`[S_+,S_-]=2S_z`$; diese beiden Relationen und ihre Adjungierten
sind äquivalent zu den drei definierenden. $`\mathbf S^2|m\rangle=(s(s+1)-m(m+1)+m^2+m)|m\rangle=s(s+1)|m\rangle`$.
*Irreduzibel:* Ein invariantes $`W\neq0`$ ist invariant unter $`S_z`$, dessen Eigenwerte einfach sind, enthält also ein
$`|m\rangle`$; die Koeffizienten von $`S_\pm`$ verschwinden nur bei $`m=\pm s`$, also liefert die Anwendung von $`S_\pm`$
alle $`|m'\rangle`$.

*Klassifikation.* Sei $`(S_a)`$ irreduzibel auf $`\mathbb C^d`$, $`\lambda`$ der größte Eigenwert von $`S_z`$ und $`v_0`$ ein
normierter Eigenvektor. Es ist $`S_zS_+v_0=S_+(S_z+1)v_0=(\lambda+1)S_+v_0`$, also $`S_+v_0=0`$ wegen der Maximalität. Sei
$`v_j=S_-^jv_0`$; dann ist $`S_zv_j=(\lambda-j)v_j`$, und per Induktion mit $`S_+S_-=S_-S_++2S_z`$ gilt
```math
S_+v_j=j(2\lambda-j+1)\,v_{j-1}
```
($`S_+v_j=S_-S_+v_{j-1}+2S_zv_{j-1}=[(j-1)(2\lambda-j+2)+2(\lambda-j+1)]v_{j-1}`$). Von null verschiedene $`v_j`$ sind
Eigenvektoren von $`S_z`$ zu verschiedenen Eigenwerten, also gibt es einen letzten von null verschiedenen, $`v_r`$. Der
Spann von $`v_0,\dots,v_r`$ ist invariant unter $`S_z,S_\pm`$, also gleich $`\mathbb C^d`$, und $`d=r+1`$. Aus
$`0=S_+v_{r+1}=(r+1)(2\lambda-r)v_r`$ folgt $`\lambda=r/2=:s`$. Schließlich ist
$`\|v_j\|^2=\langle v_{j-1},S_+v_j\rangle=j(2s-j+1)\|v_{j-1}\|^2`$, sodass die Einheitsvektoren
$`|s-j\rangle=v_j/\|v_j\|`$ die Gleichung $`S_+|m\rangle=\sqrt{(s-m)(s+m+1)}\,|m+1\rangle`$ mit $`m=s-j`$ erfüllen, also die
obige Formel. Verschiedene $`s`$ haben verschiedene Dimensionen. $`\square`$

*In der Bibliothek:* Spin $`s`$ ist `decompose.spin_matrices(s)`, an Satz 4.11 geprüft in [Tutorial, Abschnitt 4.6](../../tutorial/de/04_symmetries.md#s-4-6); der Gesamtspin von $`L`$ Spins $`\frac12`$ ist `decompose.total_spin(L)` — [Tutorial, Abschnitt 4.5](../../tutorial/de/04_symmetries.md#s-4-5).

<a id="r-4-12"></a>

**Korollar 4.12.** Jede Spindarstellung auf $`\mathbb C^n`$ ist unitär äquivalent zu
$`\bigoplus_s(\text{Spin }s)\otimes1_{m_s}`$, und mit $`d_M=\dim\ker(S_z-M)`$ gilt
```math
m_s=d_s-d_{s+1}.
```

*Beweis.* Satz 4.3 und Satz 4.11. Jede Kopie von Spin $`s'`$ enthält genau einen Eigenvektor von $`S_z`$ zum Eigenwert
$`M`$, falls $`|M|\le s'`$ und $`s'-M\in\mathbb Z`$, und sonst keinen; also ist
$`d_M=\sum_{s'\ge|M|,\,s'-M\in\mathbb Z}m_{s'}`$, und die Differenz zweier aufeinanderfolgender Werte ist $`m_s`$. $`\square`$

<a id="r-4-13"></a>

**Proposition 4.13 ($`L`$ Spins $`\frac12`$).** Der Gesamtspin auf $`(\mathbb C^2)^{\otimes L}`$ enthält Spin $`S`$ genau
dann, wenn $`S\in\{\frac L2,\frac L2-1,\dots\}`$, $`S\ge0`$, und zwar mit der Vielfachheit
```math
m_S=\binom L{\frac L2-S}-\binom L{\frac L2-S-1}\qquad\Big(\tbinom L{-1}=0\Big),
```
sodass die von $`S_x,S_y,S_z`$ erzeugte Algebra $`\bigoplus_SM_{2S+1}\otimes1_{m_S}`$ ist und ihre Kommutante
$`\bigoplus_S1_{2S+1}\otimes M_{m_S}`$.

*Beweis.* $`S_z=\sum_is_z^{(i)}`$ ist diagonal in der Produktbasis, mit Eigenwert $`u-\frac L2`$ auf einem Produkt aus
$`u`$ Spins aufwärts und $`L-u`$ Spins abwärts, also $`d_M=\binom L{L/2+M}=\binom L{L/2-M}`$, falls
$`\frac L2+M\in\{0,\dots,L\}`$, und $`0`$ sonst. Man wende Korollar 4.12 und Satz 4.3 an. Die Werte sind positiv, weil
$`\binom L{j+1}/\binom Lj=\frac{L-j}{j+1}>1`$ für $`j<\frac{L-1}2`$. $`\square`$

Für $`L=2`$: $`\frac12\otimes\frac12=1\oplus0`$ (**Clebsch-Gordan**): das Triplett $`S_-^j|{\uparrow\uparrow}\rangle`$,
$`j=0,1,2`$, symmetrisch unter der Vertauschung $`F`$ der beiden Spins, und das Singulett
$`\frac1{\sqrt2}(|{\uparrow\downarrow}\rangle-|{\downarrow\uparrow}\rangle)`$, antisymmetrisch und von $`S_+`$ und $`S_z`$
annulliert. Für $`L=4`$: $`m_2=1`$, $`m_1=4-1=3`$, $`m_0=6-4=2`$, also
```math
\{S_x,S_y,S_z\}''\cong M_5\oplus M_3\oplus M_1\ \text{ with multiplicities }1,3,2,\qquad
\{S_x,S_y,S_z\}'\cong M_1\oplus M_3\oplus M_2,
```
mit den Dimensionen $`35`$ und $`14`$ – den Zahlen aus Kapitel 10.

*Rechnung:* Die Zerlegung für $`L=4`$ allein aus den drei Erzeugern, mit Proposition 4.13 und $`\mathbf S^2=S(S+1)`$ auf Spin $`S`$ — [Tutorial, Abschnitt 4.5](../../tutorial/de/04_symmetries.md#s-4-5).

Welche Operatoren vertauschen mit allen Drehungen? Die Platzpermutationen tun es: Konjugation mit $`\pi(\sigma)`$
bildet $`s^{(i)}_a`$ auf $`s^{(\sigma(i))}_a`$ ab (Lemma B.1(b)) und hält daher $`\sum_is^{(i)}_a`$ fest. Umgekehrt tut es
nichts anderes:

<a id="r-4-14"></a>

**Satz 4.14 (Schur-Weyl-Dualität).** Sei $`\pi`$ die Darstellung von $`S_L`$, die die Faktoren von
$`(\mathbb C^2)^{\otimes L}`$ permutiert. (a) $`\mathrm{span}\,\pi(S_L)\subseteq\{S_x,S_y,S_z\}'`$. (b) Es gilt
Gleichheit; äquivalent dazu (Satz 2.2) ist die vom Gesamtspin erzeugte Algebra gleich $`\pi(S_L)'`$.

*Beweis.* Teil (a) wurde soeben bewiesen. Teil (b) ist Korollar B.7, der Fall $`d=2`$ von Satz B.6, der die
Dualität für $`(\mathbb C^d)^{\otimes L}`$ und die diagonale Wirkung von $`U(d)`$ für alle $`d`$ und $`L`$ beweist. $`\square`$

(Zum Weiterlesen: W. Fulton und J. Harris, *Representation Theory*, Springer 1991, Lecture 6; R. Goodman und
N. Wallach, *Symmetry, Representations, and Invariants*, Springer 2009, Abschnitt 9.1.)

*Rechnung:* Satz 4.14 für $`L\le5`$, die Dimensionen $`2,5,14,42`$ von $`\{S_a\}'`$ (Korollar B.10) und $`\pi(S_4)''`$, dessen Paare $`(k_c,m_c)`$ die des Gesamtspins mit vertauschten Rollen sind (Korollar B.9) — [Tutorial, Abschnitt 4.5](../../tutorial/de/04_symmetries.md#s-4-5).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Ein Spin-$`S`$-Multiplett besteht für jeden drehinvarianten Hamiltonoperator aus
> $`2S+1`$ Zuständen gleicher Energie, die sich durch $`S^z=M`$ unterscheiden. Zwei Spins $`\frac12`$ bilden ein Triplett
> und ein Singulett; vier Spins bilden ein Quintett, drei Tripletts und zwei Singuletts. Die Vielfachheitsräume
> der beiden Singuletts und der drei Tripletts sind das dekohärenzfreie Qubit und das rauschfreie Qutrit aus
> Kapitel 10.

<a id="s-4-6"></a>

## 4.6 Spinsysteme

Das **Heisenberg-Modell** auf $`L`$ Spins $`\frac12`$ ist $`H=\sum_{i<j}J_{ij}\,\mathbf s^{(i)}\cdot\mathbf s^{(j)}`$ mit
reellen Kopplungen $`J_{ij}`$.

<a id="r-4-15"></a>

**Proposition 4.15.** (a) (*Dirac-Identität*) $`\mathbf s^{(i)}\cdot\mathbf s^{(j)}=\frac12F_{ij}-\frac14`$, wobei $`F_{ij}`$
die Spins $`i`$ und $`j`$ vertauscht. Daher liegt jeder Heisenberg-Hamiltonoperator in
$`\mathrm{span}\,\pi(S_L)\subseteq\{S_a\}'`$, und $`V^*HV=\bigoplus_S1_{2S+1}\otimes H_S`$ mit Blöcken $`H_S\in M_{m_S}`$
in der Zerlegung aus Proposition 4.13: Jeder Eigenwert von $`H_S`$ ist ein $`(2S+1)`$-fach entartetes Niveau von $`H`$.
(b) Der Raum der **Vektoren höchsten Gewichts** $`\ker S_+\cap\ker(S_z-S)`$ hat die Dimension $`m_S`$, ist invariant
unter $`H`$, und $`H`$ wirkt darauf als $`H_S`$ (bis auf einen unitären Basiswechsel).

*Beweis.* (a) Auf zwei Spins ist
$`\mathbf s^{(1)}\cdot\mathbf s^{(2)}=\frac12(\mathbf S^2-\mathbf s^{(1)2}-\mathbf s^{(2)2})=\frac12\mathbf S^2-\frac34`$
(Satz 4.11 mit $`s=\frac12`$). Die Triplettvektoren $`S_-^j|{\uparrow\uparrow}\rangle`$ sind $`|{\uparrow\uparrow}\rangle`$,
$`|{\downarrow\uparrow}\rangle+|{\uparrow\downarrow}\rangle`$, $`2|{\downarrow\downarrow}\rangle`$, und zusammen mit dem Singulett bilden sie
eine Orthogonalbasis von $`\mathbb C^2\otimes\mathbb C^2`$. Auf dem Triplett ist $`\mathbf S^2=2`$ und $`F=1`$:
$`\mathbf S^2=S_-S_++S_z^2+S_z`$ ist $`0+1+1=2`$ auf $`|{\uparrow\uparrow}\rangle`$, und $`\mathbf S^2`$ und $`F`$ vertauschen mit
$`S_-`$, während $`F`$ den Vektor $`|{\uparrow\uparrow}\rangle`$ festhält. Auf dem Singulett ist $`\mathbf S^2=0`$ (es wird von
$`S_+`$ und $`S_z`$ annulliert) und $`F=-1`$. Beide Seiten stimmen auf diesen vier Vektoren überein. Für Plätze $`i,j`$ von
$`L`$ bilde man das Tensorprodukt mit der Identität. Der Rest ist Satz 4.14(a) und Satz 4.3. (b) Wir beweisen mehr:
Für jede wie in Korollar 4.12 zerlegte Spindarstellung und jedes $`\mu\in\frac12\mathbb Z`$ ist
$`\ker S_+\cap\ker(S_z-\mu)`$ (das Bild unter der Unitären aus Korollar 4.12 von) $`|\mu,\mu\rangle\otimes\mathbb C^{m_\mu}`$,
falls $`\mu\ge0`$, und $`0`$, falls $`\mu<0`$. In der Tat ist in dieser Basis ein Vektor mit $`S_z=\mu`$ von der Form
$`\sum_{S'}|S',\mu\rangle\otimes w_{S'}`$, summiert über die $`S'\ge|\mu|`$ mit $`S'-\mu\in\mathbb Z`$, und nach Satz 4.11 bildet
$`S_+`$ ihn auf $`\sum_{S'}\sqrt{(S'-\mu)(S'+\mu+1)}\,|S',\mu+1\rangle\otimes w_{S'}`$ ab, wobei $`S'+\mu+1\ge1`$; er wird also
genau dann annulliert, wenn $`w_{S'}=0`$ für alle $`S'\neq\mu`$. Für das Heisenbergsche $`H`$ wirkt $`1\otimes H_S`$ auf
$`|S,S\rangle\otimes\mathbb C^{m_S}`$ als $`H_S`$. $`\square`$

*Rechnung:* Proposition 4.15 an einem Ring aus vier Spins und an einem Ring aus acht Spins — [Tutorial, Abschnitt 4.6](../../tutorial/de/04_symmetries.md#s-4-6).

Teil (b) beschreibt, wie das Modul `su2` (Spinketten, Kapitel 9) $`H_S`$ ohne numerische Zerlegung aufbaut: Es
konstruiert eine Orthonormalbasis der Vektoren höchsten Gewichts, indem es Platz für Platz hinzufügt, mit den
Clebsch-Gordan-Koeffizienten von $`j\otimes\frac12`$. Wir leiten sie durch Absteigen vom höchsten Gewicht her. Nach
Satz 4.11 und $`S_-=S_+^*`$ erfüllt Spin $`j`$
```math
S_\pm|j,m\rangle=\sqrt{(j\mp m)(j\pm m+1)}\;|j,m\pm1\rangle ,
```
und wir setzen $`|j,m\rangle=0`$ für $`|m|>j`$ (verträglich mit der Formel, deren Koeffizient bei $`m=\pm j`$ verschwindet).

<a id="r-4-16"></a>

**Satz 4.16 (Clebsch-Gordan-Koeffizienten von $`j\otimes\frac12`$).** Sei $`j\in\{\frac12,1,\frac32,\dots\}`$, sei
$`S_a=S^{(j)}_a\otimes1+1\otimes s_a`$ auf $`\mathbb C^{2j+1}\otimes\mathbb C^2`$ (Spin $`j`$ aus Satz 4.11 mit Basis $`|j,m\rangle`$,
gekoppelt an einen Spin $`\frac12`$ mit $`|{\uparrow}\rangle=|\frac12,\frac12\rangle`$, $`|{\downarrow}\rangle=|\frac12,-\frac12\rangle`$), und
sei $`\alpha_M=\sqrt{(j+M+\frac12)/(2j+1)}`$, $`\beta_M=\sqrt{(j-M+\frac12)/(2j+1)}`$. Für $`J=j\pm\frac12`$ und
$`M=J,J-1,\dots,-J`$ sei
```math
\begin{aligned}
|j+\tfrac12,M\rangle\!\rangle&=\alpha_M\,|j,M-\tfrac12\rangle\otimes|{\uparrow}\rangle+\beta_M\,|j,M+\tfrac12\rangle\otimes|{\downarrow}\rangle,\\
|j-\tfrac12,M\rangle\!\rangle&=-\beta_M\,|j,M-\tfrac12\rangle\otimes|{\uparrow}\rangle+\alpha_M\,|j,M+\tfrac12\rangle\otimes|{\downarrow}\rangle .
\end{aligned}
```
(a) Diese $`2(2j+1)`$ Vektoren bilden eine Orthonormalbasis von $`\mathbb C^{2j+1}\otimes\mathbb C^2`$.
(b) $`|j+\frac12,j+\frac12\rangle\!\rangle=|j,j\rangle\otimes|{\uparrow}\rangle`$, und $`|j-\frac12,j-\frac12\rangle\!\rangle`$ ist bis auf einen
Faktor vom Betrag eins der einzige Einheitsvektor $`x`$ mit $`S_zx=(j-\frac12)x`$ und $`S_+x=0`$; der Faktor wird durch einen
positiven Koeffizienten von $`|j,j\rangle\otimes|{\downarrow}\rangle`$ festgelegt (die Condon-Shortley-Konvention).
(c) $`S_z|J,M\rangle\!\rangle=M|J,M\rangle\!\rangle`$ und $`S_-|J,M\rangle\!\rangle=\sqrt{(J+M)(J-M+1)}\,|J,M-1\rangle\!\rangle`$: Die Zustände
$`|J,M\rangle\!\rangle`$ entstehen aus $`|J,J\rangle\!\rangle`$ durch Absteigen und Normieren.
(d) Daher ist $`S_a|J,M\rangle\!\rangle=\sum_{M'}(S^{(J)}_a)_{M'M}|J,M'\rangle\!\rangle`$ für $`a=x,y,z`$: Die Unitäre
$`|J,M\rangle\mapsto|J,M\rangle\!\rangle`$ überführt Spin $`(j+\frac12)\oplus`$ Spin $`(j-\frac12)`$ in $`S_a`$, also ist
$`j\otimes\frac12\cong(j+\frac12)\oplus(j-\frac12)`$, und die Clebsch-Gordan-Koeffizienten $`\langle j,m;\frac12,\sigma|J,M\rangle`$
(die Koordinaten von $`|J,M\rangle\!\rangle`$ in der Produktbasis) sind $`\pm\alpha_M,\beta_M`$ wie angegeben.

*Beweis.* (a) Ist ein Basisvektor in der Formel $`0`$, so ist auch sein Koeffizient $`0`$: $`|j,M-\frac12\rangle=0`$ nur für
$`M=-j-\frac12`$, wo $`\alpha_M=0`$, und $`|j,M+\frac12\rangle=0`$ nur für $`M=j+\frac12`$, wo $`\beta_M=0`$. Also hat jeder Vektor die
Norm $`\alpha_M^2+\beta_M^2=1`$. Vektoren mit verschiedenem $`M`$ enthalten verschiedene Produktbasisvektoren, und für
gleiches $`M`$ ist das Skalarprodukt $`-\alpha_M\beta_M+\beta_M\alpha_M=0`$. Es sind $`(2j+2)+2j=2(2j+1)`$ Vektoren.
(b) Die erste Aussage ist $`\alpha_{j+1/2}=1`$, $`\beta_{j+1/2}=0`$. Der Eigenraum $`S_z=j-\frac12`$ wird von
$`|j,j-1\rangle\otimes|{\uparrow}\rangle`$ und $`|j,j\rangle\otimes|{\downarrow}\rangle`$ aufgespannt, und mit $`s_+|{\uparrow}\rangle=0`$,
$`s_+|{\downarrow}\rangle=|{\uparrow}\rangle`$ gilt
```math
S_+\big(p\,|j,j-1\rangle\otimes|{\uparrow}\rangle+q\,|j,j\rangle\otimes|{\downarrow}\rangle\big)=(p\sqrt{2j}+q)\,|j,j\rangle\otimes|{\uparrow}\rangle ,
```
was genau dann verschwindet, wenn $`q=-\sqrt{2j}\,p`$; mit $`|p|^2+|q|^2=1`$ und $`q>0`$ (also $`p<0`$) ist das
$`p=-\beta_{j-1/2}=-1/\sqrt{2j+1}`$, $`q=\alpha_{j-1/2}=\sqrt{2j/(2j+1)}`$. (c) $`S_z`$: Beide Terme haben
$`S_z=(M\mp\frac12)\pm\frac12=M`$. Für $`S_-`$ sei
$`x_M=p\,|j,M-\frac12\rangle\otimes|{\uparrow}\rangle+q\,|j,M+\frac12\rangle\otimes|{\downarrow}\rangle`$; wegen
$`s_-|{\uparrow}\rangle=|{\downarrow}\rangle`$ und $`s_-|{\downarrow}\rangle=0`$ ist
```math
S_-x_M=p\sqrt{(j+M-\tfrac12)(j-M+\tfrac32)}\;|j,M-\tfrac32\rangle\otimes|{\uparrow}\rangle
+\Big(p+q\sqrt{(j+M+\tfrac12)(j-M+\tfrac12)}\Big)|j,M-\tfrac12\rangle\otimes|{\downarrow}\rangle .
```
Wir verwenden $`\beta_M\sqrt{(j+M+\frac12)(j-M+\frac12)}=(j-M+\frac12)\alpha_M`$ und
$`\alpha_M\sqrt{(j+M+\frac12)(j-M+\frac12)}=(j+M+\frac12)\beta_M`$.
Für $`J=j+\frac12`$ ist $`(p,q)=(\alpha_M,\beta_M)`$ und $`c=\sqrt{(J+M)(J-M+1)}=\sqrt{(j+M+\frac12)(j-M+\frac32)}`$: Der erste
Koeffizient ist $`\sqrt{(j+M+\frac12)(j+M-\frac12)(j-M+\frac32)/(2j+1)}=c\,\alpha_{M-1}`$, der zweite
$`\alpha_M(j-M+\frac32)=c\,\beta_{M-1}`$.
Für $`J=j-\frac12`$ ist $`(p,q)=(-\beta_M,\alpha_M)`$ und $`c=\sqrt{(j+M-\frac12)(j-M+\frac12)}`$: Der erste Koeffizient ist
$`-\sqrt{(j-M+\frac12)(j+M-\frac12)(j-M+\frac32)/(2j+1)}=-c\,\beta_{M-1}`$, der zweite $`\beta_M(j+M-\frac12)=c\,\alpha_{M-1}`$.
In beiden Fällen ist $`S_-|J,M\rangle\!\rangle=c\,|J,M-1\rangle\!\rangle`$ für $`M>-J`$. Für $`M=-J`$ verschwinden beide Koeffizienten
($`p=\alpha_M=0`$, falls $`J=j+\frac12`$; $`j+M-\frac12=0`$, falls $`J=j-\frac12`$), also $`S_-|J,-J\rangle\!\rangle=0`$. (d) Nach (a) und
(c) sind in der Orthonormalbasis $`|J,M\rangle\!\rangle`$ die Matrizen von $`S_z`$ und $`S_-`$ die von Spin
$`(j+\frac12)\oplus`$ Spin $`(j-\frac12)`$ (Satz 4.11 und die obige Formel); dasselbe gilt für die Matrix von $`S_+=S_-^*`$,
also für die von $`S_x=\frac12(S_++S_-)`$ und $`S_y=\frac1{2i}(S_+-S_-)`$. $`\square`$

*Rechnung:* Satz 4.16 für $`j\le2`$, an dichten Matrizen und gegen `su2.clebsch_gordan` — [Tutorial, Abschnitt 4.6](../../tutorial/de/04_symmetries.md#s-4-6).

Aufgabe 2 erweitert die Zerlegung (aber nicht die Koeffizienten) auf $`j_1\otimes j_2`$. Für die Bibliothek zählen nur
die beiden Vektoren höchsten Gewichts aus Satz 4.16, und diese funktionieren in jeder Spindarstellung, ob irreduzibel
oder nicht:

<a id="r-4-17"></a>

**Korollar 4.17 (Vektoren höchsten Gewichts Platz für Platz).** Sei $`(T_a)`$ eine Spindarstellung auf $`E`$, sei
$`\mathcal H_\mu(E)=\ker T_+\cap\ker(T_z-\mu)`$ mit einer Orthonormalbasis $`B_\mu`$ ($`B_\mu=\emptyset`$, falls
$`\mathcal H_\mu(E)=0`$), und sei $`S_a=T_a\otimes1+1\otimes s_a`$ auf $`E\otimes\mathbb C^2`$. Für jedes $`S\in\{0,\frac12,1,\dots\}`$
bilden die Vektoren
```math
v\otimes|{\uparrow}\rangle\quad(v\in B_{S-1/2}),\qquad
\sqrt{\tfrac{2S+1}{2S+2}}\;v\otimes|{\downarrow}\rangle-\tfrac1{\sqrt{(2S+1)(2S+2)}}\;T_-v\otimes|{\uparrow}\rangle\quad(v\in B_{S+1/2})
```
eine Orthonormalbasis von $`\ker S_+\cap\ker(S_z-S)`$ (die erste Familie ist für $`S=0`$ leer). Ausgehend von einem Spin
($`E=\mathbb C^2`$, $`B_{1/2}=\{|{\uparrow}\rangle\}`$, alle anderen $`B_\mu`$ leer) und durch Hinzufügen von $`L-1`$ Spins erhält man
so Orthonormalbasen der Räume der Vektoren höchsten Gewichts von $`L`$ Spins $`\frac12`$, der Dimensionen $`m_S`$
(Proposition 4.15(b)).

*Beweis.* Sei $`a=\sqrt{(2S+1)/(2S+2)}`$, $`b=-1/\sqrt{(2S+1)(2S+2)}`$,
$`u(v)=a\,v\otimes|{\downarrow}\rangle+b\,T_-v\otimes|{\uparrow}\rangle`$. Für $`v\in\mathcal H_{S+1/2}(E)`$ ist
$`T_+T_-v=(T_-T_++2T_z)v=(2S+1)v`$, und $`a+(2S+1)b=0`$. *Höchstes Gewicht:* $`S_z`$ ist $`S`$ auf beiden Familien
($`T_z=S\mp\frac12`$ auf $`v`$ und $`S-\frac12`$ auf $`T_-v`$); $`S_+(v\otimes|{\uparrow}\rangle)=T_+v\otimes|{\uparrow}\rangle=0`$, und
$`S_+u(v)=a\,T_+v\otimes|{\downarrow}\rangle+a\,v\otimes|{\uparrow}\rangle+b\,T_+T_-v\otimes|{\uparrow}\rangle=(a+(2S+1)b)\,v\otimes|{\uparrow}\rangle=0`$.
*Orthonormal:* $`\langle v\otimes|{\uparrow}\rangle,u(v')\rangle=b\langle v,T_-v'\rangle=b\langle T_+v,v'\rangle=0`$ für
$`v\in B_{S-1/2}`$, und
$`\langle u(v),u(v')\rangle=a^2\langle v,v'\rangle+b^2\langle v,T_+T_-v'\rangle=(a^2+(2S+1)b^2)\langle v,v'\rangle=\langle v,v'\rangle`$.
*Aufspannend:* Sei $`S_zw=Sw`$ und $`S_+w=0`$, und man schreibe $`w=x\otimes|{\uparrow}\rangle+y\otimes|{\downarrow}\rangle`$. Dann ist
$`T_zx=(S-\frac12)x`$, $`T_zy=(S+\frac12)y`$ und $`0=S_+w=(T_+x+y)\otimes|{\uparrow}\rangle+T_+y\otimes|{\downarrow}\rangle`$, also
$`T_+y=0`$, d. h. $`y\in\mathcal H_{S+1/2}(E)`$, und $`T_+x=-y`$. Dann erfüllt $`x'=x+\frac1{2S+1}T_-y`$ die Gleichung
$`T_+x'=-y+y=0`$, also $`x'\in\mathcal H_{S-1/2}(E)`$, und das ist $`0`$ für $`S=0`$ (Beweis von Proposition 4.15(b)). Da
$`b/a=-1/(2S+1)`$, liegt $`w=x'\otimes|{\uparrow}\rangle+a^{-1}u(y)`$, mit linear fortgesetztem $`u`$, im Spann der beiden
Familien. Für $`L`$ Spins wende man dies mit $`E=(\mathbb C^2)^{\otimes(n-1)}`$ und dessen Gesamtspin $`T_a`$ an,
$`n=2,\dots,L`$; auf einem Spin ist $`\mathcal H_{1/2}=\mathbb C|{\uparrow}\rangle`$ und alle anderen $`\mathcal H_\mu=0`$ (der einzige
Spin ist $`\frac12`$). $`\square`$

Die zweite Familie besteht aus den Vektoren $`|S,S\rangle\!\rangle`$ aus Satz 4.16 mit $`j=S+\frac12`$, geschrieben innerhalb der
von den $`T_-^kv`$ aufgespannten Kopie von Spin $`j`$ (Beweis von Satz 4.11, $`|j,j-1\rangle=T_-v/\sqrt{2j}`$), und die erste
Familie ist $`|j+\frac12,j+\frac12\rangle\!\rangle`$ mit $`j=S-\frac12`$. Bezeichnet man jeden Basisvektor durch seine Folge von
Zwischenspins $`S_1=\frac12,S_2,\dots,S_L=S`$, $`|S_n-S_{n-1}|=\frac12`$, so erhält man die Yamanouchi-Basis.

*In der Bibliothek:* `su2.SU2Chain` baut die Yamanouchi-Basis mit diesen Koeffizienten auf, innerhalb des $`U(1)`$-Sektors $`\hat N=\frac L2+S`$ von `SpinChain`; Korollar 4.17 wird für acht Spins geprüft, mit den $`U(1)`$- und $`SU(2)`$-Blöcken des Heisenberg-Rings und ihrem Aufwand — [Tutorial, Abschnitt 4.6](../../tutorial/de/04_symmetries.md#s-4-6).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Das Heisenberg-Modell beschreibt die
> Austauschwechselwirkung lokalisierter Elektronenspins in magnetischen Isolatoren; $`J_{ij}>0`$ begünstigt
> antiparallele Nachbarn. Sein Spektrum besteht aus Multipletts, und ein Magnetfeld $`h\,S^z`$ spaltet jedes Multiplett
> in $`2S+1`$ äquidistante Niveaus auf (Aufgabe 4). Satz 9.9 (bipartite Kopplungen; für gerade antiferromagnetische
> Ketten und Ringe ein Singulett-Grundzustand) geht auf E. Lieb und D. Mattis zurück (*J. Math. Phys.* 3 (1962) 749).

<a id="s-4-7"></a>

## 4.7 Reelle Algebren: reeller, komplexer und quaternionischer Typ

Symmetrien sind oft durch reelle Matrizen gegeben: die Drehungen des $`\mathbb R^3`$, Permutationen, Hamiltonoperatoren,
die in einer Basis reell sind. Die komplexe Algebra, die sie erzeugen, wird durch Satz 2.6 zerlegt, aber die *reellen*
Matrizen in ihr bilden eine reelle Algebra, deren einfache Bausteine von dreierlei Art sind, $`M_k(\mathbb R)`$,
$`M_k(\mathbb C)`$ und $`M_k(\mathbb H)`$. Wir leiten sie aus der komplexen Zerlegung und einer weiteren Angabe her: was
die komplexe Konjugation mit den Sektoren macht.

Für $`a\in M_n(\mathbb C)`$ sei $`\bar a`$ die eintragsweise komplex Konjugierte (ebenso für Vektoren, und
$`\bar S=\{\bar s:s\in S\}`$ für Mengen). Die Konjugation ist antilinear und bijektiv, $`\overline{ab}=\bar a\bar b`$,
$`\overline{a^*}=\bar a^*`$ und $`\bar{\bar a}=a`$; es ist $`\bar a=a`$ genau dann, wenn $`a`$ reell ist, und dann gilt
$`a^*=a^T`$. Für eine Menge $`\mathcal B`$ von Matrizen (oder von Vektoren) sei
$`\mathcal B_{\mathbb R}=\mathcal B\cap M_n(\mathbb R)`$ (bzw. $`\mathcal B\cap\mathbb R^n`$). Das Skalarprodukt ist
$`\langle u,v\rangle=\sum_i\bar u_iv_i`$, also
$`\langle\bar u,\bar v\rangle=\overline{\langle u,v\rangle}=\langle v,u\rangle`$.

<a id="r-4-18"></a>

**Lemma 4.18 (Realteile).** (a) Sei $`\mathcal B`$ ein komplexer Unterraum von $`M_n(\mathbb C)`$ oder von $`\mathbb C^n`$
mit $`\bar{\mathcal B}=\mathcal B`$. Dann gilt $`\mathcal B=\mathcal B_{\mathbb R}\oplus i\mathcal B_{\mathbb R}`$ und
$`\dim_{\mathbb R}\mathcal B_{\mathbb R}=\dim_{\mathbb C}\mathcal B`$.
(b) Sei $`\mathcal A`$ eine $`*`$-Algebra auf $`\mathbb C^n`$ mit $`\bar{\mathcal A}=\mathcal A`$. Dann ist
$`\mathcal A_{\mathbb R}`$ eine reelle Algebra mit Einselement, die unter Transposition abgeschlossen ist; es gilt
$`\overline{\mathcal A'}=\mathcal A'`$; und $`(\mathcal A')_{\mathbb R}`$ ist die Menge der reellen $`X`$ mit $`Xa=aX`$ für
alle $`a\in\mathcal A_{\mathbb R}`$.
(c) Erfüllt $`S\subseteq M_n(\mathbb C)`$ die Bedingungen $`S^*=S`$ und $`\bar S=S`$ – etwa
$`S=\{S_1,\dots,S_J,S_1^T,\dots,S_J^T\}`$ mit reellen $`S_j`$ oder eine Gruppe reeller orthogonaler Matrizen –, so gilt
$`\overline{S''}=S''`$.

*Beweis.* (a) Für $`b\in\mathcal B`$ sind $`\mathrm{Re}\,b=\frac12(b+\bar b)`$ und
$`\mathrm{Im}\,b=\frac1{2i}(b-\bar b)`$ reelle Elemente von $`\mathcal B`$ mit
$`b=\mathrm{Re}\,b+i\mathrm{Im}\,b`$, und ein Element von $`\mathcal B_{\mathbb R}\cap i\mathcal B_{\mathbb R}`$
ist reell und rein imaginär, also $`0`$. Eine Basis $`b_1,\dots,b_r`$ des reellen Raums $`\mathcal B_{\mathbb R}`$ spannt
daher $`\mathcal B`$ über $`\mathbb C`$ auf, und sie ist linear unabhängig über $`\mathbb C`$: Aus
$`\sum_j(\alpha_j+i\beta_j)b_j=0`$ mit reellen $`\alpha_j,\beta_j`$ folgt durch Real- und Imaginärteil
$`\sum_j\alpha_jb_j=\sum_j\beta_jb_j=0`$. (b) $`\mathcal A_{\mathbb R}`$ enthält $`1`$ und ist unter Produkten
abgeschlossen, ebenso unter $`a\mapsto a^T=a^*`$. $`Y`$ vertauscht genau dann mit allen $`t\in T`$, wenn $`\bar Y`$ mit allen
$`\bar t`$ vertauscht, also gilt $`\overline{T'}=\bar T'`$ für jede Menge $`T`$, und
$`\overline{\mathcal A'}=\bar{\mathcal A}'=\mathcal A'`$. Ein reelles $`X`$, das mit $`\mathcal A_{\mathbb R}`$ vertauscht,
vertauscht mit $`\mathcal A_{\mathbb R}+i\mathcal A_{\mathbb R}=\mathcal A`$ (nach (a)). (c) Man wende
$`\overline{T'}=\bar T'`$ zweimal an: $`\overline{S''}=\bar S''=S''`$. $`\square`$

**Feste Bezeichnungen.** Im Rest dieses Abschnitts ist $`\mathcal A`$ eine $`*`$-Algebra auf $`\mathbb C^n`$ mit
$`\bar{\mathcal A}=\mathcal A`$, und $`V`$ ist eine Zerlegung wie in Satz 2.6. $`V_c`$ ist der $`n\times k_cm_c`$-Spaltenblock
von $`V`$, der zum $`c`$-ten Summanden gehört, sodass $`V_c^*V_d=\delta_{cd}1`$ und $`V_cV_c^*=P_c`$, die minimalen zentralen
Projektionen, und
```math
aV_c=V_c\big(\pi_c(a)\otimes1_{m_c}\big)\quad(a\in\mathcal A),\qquad XV_c=V_c\big(1_{k_c}\otimes x_c\big)\quad(X\in\mathcal A'),
```
mit $`\pi_c(a)\in M_{k_c}`$ und $`x_c\in M_{m_c}`$. Jedes $`\pi_c`$ ist ein unitaler $`*`$-Homomorphismus von $`\mathcal A`$ auf
$`M_{k_c}`$: Für $`x\in M_{k_c}`$ hat das Element $`V_c(x\otimes1)V_c^*`$ von $`\mathcal A`$ die Komponenten $`\pi_c=x`$ und
$`\pi_e=0`$ für $`e\neq c`$. Schließlich sei $`\Omega_2=\begin{pmatrix}0&-1\\1&0\end{pmatrix}`$ und
$`\Omega_{2r}=1_r\otimes\Omega_2`$ (blockdiagonal mit $`r`$ Kopien von $`\Omega_2`$).

Eine Abbildung $`\theta:E\to E`$ auf einem Unterraum $`E\subseteq\mathbb C^k`$ heißt **antiunitär**, wenn
$`\theta(\alpha u+\beta v)=\bar\alpha\theta u+\bar\beta\theta v`$ und
$`\langle\theta u,\theta v\rangle=\langle v,u\rangle`$ gilt. Für eine unitäre Matrix $`J\in U(k)`$ ist $`\theta_Jv=J\bar v`$
antiunitär auf $`\mathbb C^k`$ ($`\langle J\bar u,J\bar v\rangle=\langle\bar u,\bar v\rangle=\langle v,u\rangle`$), mit
$`\theta_J^2v=J\bar Jv`$.

<a id="r-4-19"></a>

**Lemma 4.19 (antiunitäre Abbildungen).** Sei $`\theta`$ antiunitär auf $`E\neq0`$ mit
$`\theta^2=\varepsilon\,\mathrm{id}`$, $`\varepsilon\in\{1,-1\}`$.
(a) Ist $`\varepsilon=1`$, so hat $`E`$ eine Orthonormalbasis aus Vektoren $`w`$ mit $`\theta w=w`$.
(b) Ist $`\varepsilon=-1`$, so gilt $`\langle u,\theta u\rangle=0`$ für alle $`u\in E`$; für jeden Einheitsvektor $`u`$ ist
$`\mathrm{span}\{u,\theta u\}`$ eine $`\theta`$-invariante Ebene; und $`E`$ hat eine Orthonormalbasis
$`w_1,\theta w_1,\dots,w_r,\theta w_r`$. Insbesondere ist $`\dim E=2r`$ gerade.
(c) Sei $`J\in U(k)`$ mit $`J\bar J=\varepsilon1`$. Ist $`\varepsilon=1`$, so gibt es eine unitäre Matrix $`W`$ mit
$`J\bar W=W`$, d. h. $`J=WW^T`$. Ist $`\varepsilon=-1`$, so ist $`k`$ gerade, und es gibt eine unitäre Matrix $`W`$ mit
$`J\bar W=W\Omega_k`$, d. h. $`J=W\Omega_kW^T`$.

*Beweis.* (a) Sei $`F=\{v\in E:\theta v=v\}`$, ein reeller Unterraum. Jedes $`v\in E`$ ist $`v_1+iv_2`$ mit
$`v_1=\frac12(v+\theta v)`$ und $`v_2=\frac1{2i}(v-\theta v)`$, und beide liegen in $`F`$
($`\theta v_2=\frac1{-2i}(\theta v-v)=v_2`$). Für $`u,v\in F`$ ist
$`\langle u,v\rangle=\langle\theta u,\theta v\rangle=\langle v,u\rangle`$ reell, also ist $`\langle\cdot,\cdot\rangle`$
ein reelles Skalarprodukt auf $`F`$; sei $`w_1,\dots,w_s`$ eine Orthonormalbasis von $`F`$ (Gram–Schmidt). Die $`w_j`$ sind
orthonormal in $`E`$, also linear unabhängig über $`\mathbb C`$, und sie spannen $`E`$ über $`\mathbb C`$ auf, da $`E=F+iF`$.
(b) Wendet man $`\langle\theta x,\theta y\rangle=\langle y,x\rangle`$ auf $`x=\theta u`$, $`y=u`$ an, so folgt
$`\langle u,\theta u\rangle=\langle\theta^2u,\theta u\rangle=-\langle u,\theta u\rangle`$. Für einen Einheitsvektor $`u`$
ist also $`\theta u`$ ein zu $`u`$ orthogonaler Einheitsvektor, und $`E_1=\mathrm{span}\{u,\theta u\}`$ ist invariant:
$`\theta(\alpha u+\beta\theta u)=\bar\alpha\theta u-\bar\beta u`$. Auch das orthogonale Komplement von $`E_1`$ in $`E`$ ist
invariant: Für $`x\perp E_1`$ und $`y\in E_1`$ ist $`y=-\theta(\theta y)`$ und
$`\langle\theta x,y\rangle=-\langle\theta x,\theta(\theta y)\rangle=-\langle\theta y,x\rangle=0`$, da $`\theta y\in E_1`$.
Induktion nach $`\dim E`$ liefert die Basis. (c) $`\theta_J`$ ist antiunitär auf $`\mathbb C^k`$ mit
$`\theta_J^2=J\bar J=\varepsilon`$. Sei $`W`$ die unitäre Matrix, deren Spalten die Basis aus (a) bzw. die Basis
$`w_1,\theta_Jw_1,w_2,\theta_Jw_2,\dots`$ aus (b) sind. Die Spalten von $`J\bar W`$ sind die Bilder der Spalten von $`W`$
unter $`\theta_J`$. Im Fall (a) besagt das $`J\bar W=W`$. Im Fall (b) sind die Bilder $`\theta_Jw_j`$ und
$`\theta_J^2w_j=-w_j`$, also die Spalten $`2j-1`$ und $`2j`$ von $`W\Omega_k`$ (da $`\Omega_2e_1=e_2`$, $`\Omega_2e_2=-e_1`$).
Wegen $`\bar W^{-1}=\overline{W^*}=W^T`$ folgt $`J=WW^T`$ bzw. $`J=W\Omega_kW^T`$. $`\square`$

<a id="r-4-20"></a>

**Satz 4.20 (Konjugation der Sektoren).** (a) Es gibt eine Permutation $`\sigma`$ der Sektoren mit
$`\sigma\circ\sigma=\mathrm{id}`$ und unitäre Matrizen $`J_c\in U(k_c)`$, $`M_c\in U(m_c)`$, sodass
$`\bar P_c=P_{\sigma(c)}`$, $`k_{\sigma(c)}=k_c`$, $`m_{\sigma(c)}=m_c`$ und für jedes $`c`$ und alle $`a\in\mathcal A`$
Folgendes gilt:
```math
\overline{V_c}=V_{\sigma(c)}\,(J_c\otimes M_c),\qquad\pi_{\sigma(c)}(\bar a)=J_c\,\overline{\pi_c(a)}\,J_c^*.
```
$`J_c\otimes M_c`$ ist eindeutig, und $`J_c,M_c`$ sind eindeutig bis auf $`(J_c,M_c)\mapsto(\lambda J_c,\bar\lambda M_c)`$
mit $`|\lambda|=1`$.
(b) Ist $`\sigma(c)=c`$, so gibt es ein $`\varepsilon_c\in\{1,-1\}`$ mit $`J_c\bar J_c=\varepsilon_c1_{k_c}`$ und
$`M_c\bar M_c=\varepsilon_c1_{m_c}`$. Ist $`\varepsilon_c=-1`$, so sind $`k_c`$ und $`m_c`$ gerade.
(c) $`\sigma`$ und die $`\varepsilon_c`$ hängen nicht von der Wahl von $`V`$ ab.
(d) Für $`y\in M_{m_c}`$ gilt
$`\overline{V_c(1\otimes y)V_c^*}=V_{\sigma(c)}\big(1\otimes M_c\bar yM_c^*\big)V_{\sigma(c)}^*`$. Hat also
$`X\in\mathcal A'`$ die Blöcke $`x_c`$, so hat $`\bar X\in\mathcal A'`$ die Blöcke $`(\bar X)_{\sigma(c)}=M_c\bar x_cM_c^*`$:
Auf der Kommutante spielt $`M_c`$ die Rolle, die $`J_c`$ auf $`\mathcal A`$ spielt.

*Beweis.* (a) Die Konjugation bildet $`\mathcal A`$ auf $`\mathcal A`$ ab und erhält Summen, Produkte und Adjungierte;
also bildet sie das Zentrum $`Z(\mathcal A)`$ auf sich ab und die Menge der von null verschiedenen Projektionen von
$`Z(\mathcal A)`$ bijektiv auf sich, unter Erhaltung der Ordnung $`p\le q`$ ($`pq=p`$). Daher permutiert sie die minimalen
unter ihnen, und das sind die $`P_c`$ (Lemma 2.3): $`\bar P_c=P_{\sigma(c)}`$, und $`\sigma\circ\sigma=\mathrm{id}`$, weil
$`\bar{\bar P}_c=P_c`$. Sei $`d=\sigma(c)`$. Das Bild von $`\overline{V_c}`$ ist das Konjugierte des Bildes von $`P_c`$, also
das Bild von $`\bar P_c=P_d`$; somit $`\overline{V_c}=P_d\overline{V_c}=V_dT`$ mit $`T=V_d^*\overline{V_c}`$, und
$`T^*T=\overline{V_c}^*P_d\overline{V_c}=\overline{V_c^*V_c}=1`$. Wegen
$`k_cm_c=\mathrm{rank}\,P_c=\mathrm{rank}\,\bar P_c=k_dm_d`$ ist die Isometrie $`T`$ quadratisch, also unitär.
Konjugiert man $`aV_c=V_c(\pi_c(a)\otimes1)`$ und benutzt $`\bar aV_d=V_d(\pi_d(\bar a)\otimes1)`$, so erhält man
$`V_d(\pi_d(\bar a)\otimes1)T=V_dT(\overline{\pi_c(a)}\otimes1)`$; Multiplikation mit $`V_d^*`$ ergibt
```math
\big(\pi_d(\bar a)\otimes1_{m_d}\big)\,T=T\,\big(\overline{\pi_c(a)}\otimes1_{m_c}\big)\qquad(a\in\mathcal A).
```
Man schreibe $`T=\sum_{j,l}T_{jl}\otimes|j\rangle\langle l|`$ mit $`k_d\times k_c`$-Blöcken $`T_{jl}`$. Die Gleichung besagt
$`\pi_d(\bar a)T_{jl}=T_{jl}\overline{\pi_c(a)}`$: Jedes $`T_{jl}`$ ist ein Verflechtungsoperator von
$`\rho(a)=\overline{\pi_c(a)}`$ in $`\rho'(a)=\pi_d(\bar a)`$. Das sind Darstellungen von $`\mathcal A`$ mit seiner
Adjungierten im Sinne von Abschnitt 4.1 ($`\rho(a^*)=\rho(a)^*`$, $`\rho'(a^*)=\rho'(a)^*`$), und sie sind irreduzibel, da
ihre Bilder ganz $`M_{k_c}`$ bzw. ganz $`M_{k_d}`$ sind. Wegen $`T\neq0`$ ist ein $`T_{jl}\neq0`$, also $`k_c=k_d`$ und
$`\mathrm{Hom}(\rho,\rho')=\mathbb CJ_c`$ mit einer unitären Matrix $`J_c`$ (Lemma 4.2). Daher ist
$`T_{jl}=\mu_{jl}J_c`$, d. h. $`T=J_c\otimes M_c`$ mit $`M_c=(\mu_{jl})`$; es gilt $`1=T^*T=1\otimes M_c^*M_c`$ und
$`m_c=k_cm_c/k_c=m_d`$, also ist $`M_c`$ unitär, und die Gleichung wird zur Formel für $`\pi_d(\bar a)`$. $`T`$ ist durch $`V`$
bestimmt; ist auch $`T=J'\otimes M'`$ mit unitärem $`J'`$, so zeigen die Blöcke $`T_{jl}=M'_{jl}J'`$, dass
$`J'\in\mathrm{Hom}(\rho,\rho')`$, also $`J'=\lambda J_c`$ mit $`|\lambda|=1`$ und dann $`M'=\bar\lambda M_c`$.

(b) Sei $`\sigma(c)=c`$ und $`J=J_c`$, $`M=M_c`$. Konjugieren von $`\overline{V_c}=V_c(J\otimes M)`$ ergibt
$`V_c=\overline{V_c}(\bar J\otimes\bar M)=V_c(J\bar J\otimes M\bar M)`$, also $`J\bar J\otimes M\bar M=1`$. Zweimalige
Anwendung der Formel aus (a) liefert mit $`\overline{J^*}=J^T`$
```math
\pi_c(a)=\pi_c(\bar{\bar a})=J\,\overline{J\,\overline{\pi_c(a)}\,J^*}\,J^*=(J\bar J)\,\pi_c(a)\,(J\bar J)^*,
```
also vertauscht die unitäre Matrix $`J\bar J`$ mit $`\pi_c(\mathcal A)=M_{k_c}`$ und ist ein Skalar $`\lambda1`$ (Lemma 2.5,
Schritt 5, mit $`m=1`$), $`|\lambda|=1`$. Dann ist $`\bar JJ=\overline{J\bar J}=\bar\lambda1`$ und
$`\lambda J=(J\bar J)J=J(\bar JJ)=\bar\lambda J`$, also $`\lambda=\bar\lambda\in\{1,-1\}`$; wir nennen es $`\varepsilon_c`$.
Es ändert sich nicht unter $`J\mapsto\lambda J`$, und $`\varepsilon_c1\otimes M\bar M=1`$ ergibt $`M\bar M=\varepsilon_c1`$.
Ist $`\varepsilon_c=-1`$, so ist $`(-1)^{k_c}=\det(J\bar J)=\det J\cdot\overline{\det J}=|\det J|^2>0`$, also ist $`k_c`$
gerade, und ebenso $`m_c`$.

(c) Die $`P_c`$ hängen nicht von $`V`$ ab (Satz 2.6), also auch $`\sigma`$ nicht. Sei $`V'`$ eine weitere Zerlegung, so
nummeriert, dass $`V'_cV_c'^*=P_c`$ (Satz 2.6), mit Homomorphismen $`\pi'_c`$ und Daten $`J'_c,M'_c`$, und sei
$`\sigma(c)=c`$. Dann ist $`V'_c=P_cV'_c=V_cY`$ mit $`Y=V_c^*V'_c`$, einer quadratischen Isometrie, also unitär, und
$`(\pi_c(a)\otimes1)Y=Y(\pi'_c(a)\otimes1)`$ für $`a\in\mathcal A`$. Das Argument aus (a) mit den irreduziblen
Darstellungen $`\pi'_c,\pi_c`$ anstelle von $`\rho,\rho'`$ ergibt $`Y=u\otimes w`$ mit unitären $`u,w`$. Folglich ist
$`\overline{V'_c}=\overline{V_c}(\bar u\otimes\bar w)=V_c(J_c\bar u\otimes M_c\bar w)=V'_c(u^*J_c\bar u\otimes w^*M_c\bar w)`$,
also $`J'_c=\lambda u^*J_c\bar u`$ mit $`|\lambda|=1`$ nach der Eindeutigkeit in (a), und wegen
$`\bar uu^T=\overline{uu^*}=1`$ gilt $`J'_c\bar J'_c=u^*J_c\bar uu^T\bar J_cu=u^*J_c\bar J_cu=\varepsilon_c1`$.

(d) Nach (a) ist
$`\overline{V_c(1\otimes y)V_c^*}=V_{\sigma(c)}(J_c\otimes M_c)(1\otimes\bar y)(J_c\otimes M_c)^*V_{\sigma(c)}^*`$, und
$`(J_c\otimes M_c)(1\otimes\bar y)(J_c\otimes M_c)^*=J_cJ_c^*\otimes M_c\bar yM_c^*=1\otimes M_c\bar yM_c^*`$. Man
summiere über $`c`$ für $`X=\sum_cV_c(1\otimes x_c)V_c^*`$. $`\square`$

Ein Sektor $`c`$ ist vom **reellen Typ**, wenn $`\sigma(c)=c`$ und $`\varepsilon_c=1`$, vom **quaternionischen Typ**, wenn
$`\sigma(c)=c`$ und $`\varepsilon_c=-1`$, und vom **komplexen Typ**, wenn $`\sigma(c)\neq c`$; nach Satz 4.20(c) ist das
eine Eigenschaft von $`\mathcal A`$ und $`P_c`$. Ein Sektor mit $`\sigma(c)=c`$ und ungeradem $`k_c`$ oder ungeradem $`m_c`$ ist
nach Satz 4.20(b) vom reellen Typ.

*In der Bibliothek:* $`\sigma`$, $`J_c`$, $`M_c`$, $`\varepsilon_c`$ und die Typen der Sektoren berechnet `Decomposition.real_structure()` — [Tutorial, Abschnitt 4.7](../../tutorial/de/04_symmetries.md#s-4-7).

<a id="r-4-21"></a>

**Satz 4.21 (die reelle Algebra).** $`\Gamma`$ enthalte jeden Sektor $`c`$ mit $`\sigma(c)=c`$ und einen Sektor aus jedem
Paar $`\{c,\sigma(c)\}`$ mit $`\sigma(c)\neq c`$. Für $`\sigma(c)=c`$ sei $`R_c=\{x\in M_{k_c}(\mathbb C):xJ_c=J_c\bar x\}`$.
Dann ist $`a\mapsto(\pi_c(a))_{c\in\Gamma}`$ ein Isomorphismus reeller Algebren, der die Transponierte in die
Adjungierte überführt,
```math
\mathcal A_{\mathbb R}\;\cong\bigoplus_{c\in\Gamma,\ \sigma(c)=c}R_c\;\oplus\bigoplus_{c\in\Gamma,\ \sigma(c)\neq c}M_{k_c}(\mathbb C),
```
wobei die letzten Summanden als reelle Algebren aufgefasst werden (**komplexer Typ**). Ferner:
(a) (*reeller Typ*) Ist $`\varepsilon_c=1`$, so schreibe man $`J_c=W_JW_J^T`$ und $`M_c=W_MW_M^T`$ mit $`J_c\bar W_J=W_J`$,
$`M_c\bar W_M=W_M`$ (Lemma 4.19(c)). Dann bildet $`x\mapsto W_J^*xW_J`$ die Menge $`R_c`$ isomorph auf $`M_{k_c}(\mathbb R)`$
ab. Die Spalten von $`O_c=V_c(W_J\otimes W_M)`$ sind reell und orthonormal, sie spannen das Bild von $`P_c`$ über
$`\mathbb C`$ auf, und $`aO_c=O_c(\alpha_c(a)\otimes1_{m_c})`$ mit der reellen Matrix $`\alpha_c(a)=W_J^*\pi_c(a)W_J`$ für
$`a\in\mathcal A_{\mathbb R}`$.
(b) (*quaternionischer Typ*) Ist $`\varepsilon_c=-1`$, so schreibe man $`J_c=W_J\Omega_{k_c}W_J^T`$ mit
$`J_c\bar W_J=W_J\Omega_{k_c}`$ (Lemma 4.19(c)). Dann ist $`x\mapsto\mathbb M^{-1}(W_J^*xW_J)`$ ein Isomorphismus reeller
Algebren $`R_c\cong M_{k_c/2}(\mathbb H)`$, der die Adjungierte in die quaternionische Adjungierte
$`(q^*)_{ij}=\bar q_{ji}`$ überführt; dabei ersetzt $`\mathbb M:M_r(\mathbb H)\to M_{2r}(\mathbb C)`$ jeden Eintrag $`q`$
durch den $`2\times2`$-Block $`M(q)`$ aus Lemma E.15(2), und sein Bild ist
$`\{y\in M_{2r}(\mathbb C):y\Omega_{2r}=\Omega_{2r}\bar y\}`$.
Insbesondere ist $`\dim_{\mathbb R}\mathcal A_{\mathbb R}=\sum_ck_c^2=\dim_{\mathbb C}\mathcal A`$.

*Beweis.* Für reelles $`a`$ ist $`\bar a=a`$, und Satz 4.20(a) lautet $`\pi_{\sigma(c)}(a)=J_c\overline{\pi_c(a)}J_c^*`$.
Für $`\sigma(c)=c`$ besagt das $`\pi_c(a)\in R_c`$. Also bildet $`\Phi(a)=(\pi_c(a))_{c\in\Gamma}`$ die Algebra
$`\mathcal A_{\mathbb R}`$ in die rechte Seite ab; es ist ein Homomorphismus reeller Algebren, und
$`\Phi(a^T)=\Phi(a^*)=\Phi(a)^*`$. Ist $`\Phi(a)=0`$, so liefert dieselbe Formel $`\pi_{\sigma(c)}(a)=0`$ für $`c\in\Gamma`$,
also sind alle $`\pi_e(a)=0`$ und $`a=0`$. *Surjektivität:* Seien $`c\in\Gamma`$, $`x\in R_c`$, falls $`\sigma(c)=c`$, und sonst
$`x\in M_{k_c}`$, und $`a=V_c(x\otimes1)V_c^*\in\mathcal A`$. Nach Satz 4.20(a) für den Sektor $`\sigma(e)`$ ist
$`\pi_e(\bar a)=J_{\sigma(e)}\overline{\pi_{\sigma(e)}(a)}J_{\sigma(e)}^*`$, und das ist $`0`$, außer wenn $`\sigma(e)=c`$.
Ist $`\sigma(c)=c`$, so ist $`\pi_c(\bar a)=J_c\bar xJ_c^*=x`$ (Definition von $`R_c`$), also hat die reelle Matrix
$`\mathrm{Re}\,a=\frac12(a+\bar a)\in\mathcal A_{\mathbb R}`$ den Wert $`\Phi(\mathrm{Re}\,a)=x`$ an der Stelle
$`c`$ und $`0`$ sonst. Ist $`\sigma(c)\neq c`$, so ist $`\pi_c(\bar a)=0`$, und $`2\mathrm{Re}\,a`$ hat $`\pi_c=x`$ und
$`\pi_e=0`$ für $`e\notin\{c,\sigma(c)\}`$; wegen $`\sigma(c)\notin\Gamma`$ ist $`\Phi(2\mathrm{Re}\,a)`$ gleich $`x`$ an
der Stelle $`c`$ und $`0`$ sonst. Summen solcher Elemente ergeben die ganze rechte Seite.

(a) Man multipliziere $`xW_JW_J^T=W_JW_J^T\bar x`$ von links mit $`W_J^*`$ und von rechts mit $`\bar W_J`$; wegen
$`W_J^T\bar W_J=\overline{W_J^*W_J}=1`$ wird daraus $`W_J^*xW_J=W_J^T\bar x\bar W_J=\overline{W_J^*xW_J}`$, und die
Schritte lassen sich umkehren. Also ist $`R_c=W_JM_{k_c}(\mathbb R)W_J^*`$, und $`x\mapsto W_J^*xW_J`$ ist ein mit der
Adjungierten verträglicher Isomorphismus reeller Algebren. Weiter ist
$`\overline{O_c}=\overline{V_c}(\bar W_J\otimes\bar W_M)=V_c(J_c\bar W_J\otimes M_c\bar W_M)=V_c(W_J\otimes W_M)=O_c`$,
$`O_c^*O_c=1`$ und $`aO_c=V_c(\pi_c(a)W_J\otimes W_M)=O_c(W_J^*\pi_c(a)W_J\otimes1)`$, mit reellem $`W_J^*\pi_c(a)W_J`$ für
$`a\in\mathcal A_{\mathbb R}`$.

(b) Dieselbe Rechnung mit $`J_c=W_J\Omega W_J^T`$ ($`\Omega=\Omega_{k_c}`$) zeigt: $`xJ_c=J_c\bar x`$ gilt genau dann, wenn
$`y=W_J^*xW_J`$ die Gleichung $`y\Omega=\Omega\bar y`$ erfüllt. Man schreibe $`y`$ in $`2\times2`$-Blöcken $`y_{ij}`$; da
$`\Omega`$ blockdiagonal ist, gilt $`y\Omega=\Omega\bar y`$ genau dann, wenn $`y_{ij}\Omega_2=\Omega_2\bar y_{ij}`$ für alle
$`i,j`$, und für $`y_{ij}=\begin{pmatrix}\alpha&\beta\\\gamma&\delta\end{pmatrix}`$ lautet das
$`\begin{pmatrix}\beta&-\alpha\\\delta&-\gamma\end{pmatrix}=\begin{pmatrix}-\bar\gamma&-\bar\delta\\\bar\alpha&\bar\beta\end{pmatrix}`$,
d. h. $`\delta=\bar\alpha`$ und $`\beta=-\bar\gamma`$:
$`y_{ij}=\begin{pmatrix}\alpha&-\bar\gamma\\\gamma&\bar\alpha\end{pmatrix}=M(\alpha,\bar\gamma)`$ mit
$`M(a,b)=\begin{pmatrix}a&-b\\\bar b&\bar a\end{pmatrix}`$. Also ist das Bild von $`\mathbb M`$ wie angegeben. Nach Lemma
E.15(2) ist $`M`$ ein injektiver reell-linearer Homomorphismus $`\mathbb H\to M_2(\mathbb C)`$ mit $`M(\bar q)=M(q)^*`$; da
ein Matrixprodukt nur Summen und Produkte von Einträgen verwendet, ist $`\mathbb M`$ ein injektiver Homomorphismus
reeller Algebren mit $`\mathbb M(q^*)=\mathbb M(q)^*`$ (wie im Beweis von Satz E.24(1)). Daher ist
$`x\mapsto\mathbb M^{-1}(W_J^*xW_J)`$ der behauptete Isomorphismus. *Dimensionen:* $`\dim_{\mathbb R}R_c=k_c^2`$ in beiden
Fällen ($`\dim_{\mathbb R}M_r(\mathbb H)=4r^2`$ mit $`2r=k_c`$), und ein Paar trägt
$`\dim_{\mathbb R}M_{k_c}(\mathbb C)=2k_c^2=k_c^2+k_{\sigma(c)}^2`$ bei. $`\square`$

Dasselbe gilt für die Kommutante, mit vertauschten Faktoren: $`(\mathcal A')_{\mathbb R}`$ ist die direkte Summe von
$`M_{m_c}(\mathbb R)`$, $`M_{m_c/2}(\mathbb H)`$ und $`M_{m_c}(\mathbb C)`$ über dieselben Sektoren, mit demselben Typ. Denn
$`\mathcal A'`$ ist eine $`*`$-Algebra mit $`\overline{\mathcal A'}=\mathcal A'`$ (Lemma 4.18(b)), und der obige Beweis hat
nur drei Tatsachen benutzt, die für $`\mathcal A'`$ mit $`X\mapsto x_c`$ anstelle von $`\pi_c`$ und $`M_c`$ anstelle von $`J_c`$
gelten: $`X\mapsto x_c`$ ist ein unitaler $`*`$-Homomorphismus auf $`M_{m_c}`$ ($`V_c(1\otimes y)V_c^*\in\mathcal A'`$ hat den
Block $`y`$ im Sektor $`c`$ und $`0`$ sonst), $`(\bar X)_{\sigma(c)}=M_c\bar x_cM_c^*`$ (Satz 4.20(d)), und
$`M_c\bar M_c=\varepsilon_c1`$ für $`\sigma(c)=c`$ (Satz 4.20(b)).

*In der Bibliothek:* Eine an Satz 4.21 angepasste reelle Orthogonalbasis, `Decomposition.real_basis()` — [Tutorial, Abschnitt 4.7](../../tutorial/de/04_symmetries.md#s-4-7).

Drei Beispiele. (i) Die Erzeuger $`(L_a)_{bc}=-\epsilon_{abc}`$ der Drehungen
des $`\mathbb R^3`$, die auf $`\mathbb R^3\otimes\mathbb R^2`$ als $`L_a\otimes1_2`$ wirken: Es gibt einen Sektor, $`k=3`$, $`m=2`$ ($`iL_x,iL_y,iL_z`$ ist eine Spindarstellung auf $`\mathbb C^3`$, in der $`iL_z`$ die
einfachen Eigenwerte $`1,0,-1`$ hat, also nach Korollar 4.12 Spin $`1`$); er ist vom reellen Typ, weil $`k`$ ungerade ist, und $`\mathcal A_{\mathbb R}\cong M_3(\mathbb R)`$.
(ii) Der Drehungserzeuger $`\Omega_2`$ des $`\mathbb R^2`$, der auf $`\mathbb R^2\otimes\mathbb R^2`$ als
$`\Omega_2\otimes1_2`$ wirkt: Wegen $`\Omega_2^2=-1`$ und $`\Omega_2^*=-\Omega_2`$ ist die Algebra
$`\mathrm{span}\{1,\Omega_2\otimes1\}`$; die zueinander konjugierten Eigenvektoren $`(1,\mp i)/\sqrt2`$ von
$`\Omega_2`$ zu $`\pm i`$ ergeben zwei Sektoren mit $`k=1`$, $`m=2`$, die von der Konjugation vertauscht werden: ein Paar vom
komplexen Typ, und $`\mathcal A_{\mathbb R}=\{\alpha+\beta\,\Omega_2\otimes1\}\cong\mathbb C`$. (iii) Die
Linksmultiplikation $`L_q`$ mit einer Quaternion $`q=(a,b)\in\mathbb H`$ (Abschnitt E.6) auf $`\mathbb H=\mathbb R^4`$, mit
den Koordinaten $`(\mathrm{Re}\,a,\mathrm{Im}\,a,\mathrm{Re}\,b,\mathrm{Im}\,b)`$. $`q\mapsto L_q`$ ist
ein injektiver Homomorphismus reeller Algebren ($`\mathbb H`$ ist assoziativ, Lemma E.15(2)), und $`L_{\bar q}=L_q^T`$,
weil $`\langle qx,y\rangle=\langle x,\bar qy\rangle`$ (Lemma E.14(2)) für die Norm $`N(a,b)=|a|^2+|b|^2`$, die euklidische
Norm des $`\mathbb R^4`$. Die von $`L_I,L_J,L_K`$ ($`I=(i,0)`$, $`J=(0,1)`$, $`K=(0,i)`$) und ihren Transponierten erzeugte
reelle Algebra ist also $`\{L_q:q\in\mathbb H\}\cong\mathbb H`$; die komplexe Algebra $`\mathcal A`$ ist der komplexe
Spann der Produkte der Erzeuger, $`\{L_q\}+i\{L_q\}`$, also $`\mathcal A_{\mathbb R}=\{L_q\}\cong\mathbb H`$ (Lemma
4.18(a)). Das legt die Zerlegung bereits fest:
$`\sum_ck_c^2=\dim_{\mathbb C}\mathcal A=\dim_{\mathbb R}\mathcal A_{\mathbb R}=4`$ (Satz 4.21), und nicht alle $`k_c`$
sind $`1`$, denn sonst wäre $`\mathcal A`$ und damit $`\mathcal A_{\mathbb R}\cong\mathbb H`$ kommutativ; also gibt es einen
einzigen Sektor, mit $`k=2`$ und $`\sigma(c)=c`$, und er ist nicht vom reellen Typ, da Satz 4.21(a) sonst
$`\mathcal A_{\mathbb R}\cong M_2(\mathbb R)`$ ergäbe, das Nullteiler hat, während $`\mathbb H`$ keine hat (Lemma
E.15(3)): Er ist vom quaternionischen Typ, mit $`m=4/2=2`$ auf $`\mathbb R^4`$. (Die Quaternionen der Norm $`1`$ wirken hier
als $`SU(2)`$: $`M(q)`$ ist unitär mit Determinante $`|a|^2+|b|^2=1`$, und jedes
$`U=\begin{pmatrix}\alpha&\beta\\\gamma&\delta\end{pmatrix}\in SU(2)`$ ist von dieser Form, denn
$`U^*=U^{-1}=\begin{pmatrix}\delta&-\beta\\-\gamma&\alpha\end{pmatrix}`$ ergibt $`\delta=\bar\alpha`$,
$`\beta=-\bar\gamma`$, also $`U=M(\alpha,\bar\gamma)`$.) Auf $`\mathbb H^2=\mathbb R^8`$, wo $`q`$ als $`1_2\otimes L_q`$ wirkt,
ergibt dasselbe Argument einen Sektor mit $`k=2`$, $`m=4`$, vom quaternionischen Typ,
$`\mathcal A_{\mathbb R}\cong M_1(\mathbb H)=\mathbb H`$.

*Rechnung:* Die drei Beispiele, mit `Decomposition.real_structure()` und `Decomposition.real_basis()` — [Tutorial, Abschnitt 4.7](../../tutorial/de/04_symmetries.md#s-4-7).

In der Darstellungstheorie einer Gruppe durch reelle Matrizen ist das Folgende die reelle Form des Schurschen Lemmas;
es zeigt, dass sich die drei Typen an der reellen Darstellung selbst ablesen lassen.

<a id="r-4-22"></a>

**Satz 4.22 (reelles Schursches Lemma).** Sei $`U\subseteq\mathbb R^n`$ ein von null verschiedener Unterraum mit
$`aU\subseteq U`$ für alle $`a\in\mathcal A_{\mathbb R}`$, der **irreduzibel** ist: $`0`$ und $`U`$ sind seine einzigen
Unterräume mit dieser Eigenschaft. Sei $`D_U`$ die Menge der reell-linearen Abbildungen $`T:U\to U`$ mit $`Tau=aTu`$ für
alle $`a\in\mathcal A_{\mathbb R}`$, $`u\in U`$.
(a) $`D_U`$ ist eine assoziative reelle Algebra mit Einselement und ohne Nullteiler, also isomorph zu $`\mathbb R`$,
$`\mathbb C`$ oder $`\mathbb H`$ (Satz E.19).
(b) Für genau eine Bahn $`\mathcal O=\{c,\sigma(c)\}`$ von $`\sigma`$ liegt $`U`$ im Bild von
$`P_{\mathcal O}=\sum_{e\in\mathcal O}P_e`$, und $`\dim_{\mathbb R}D_U=1,2,4`$ und $`\dim_{\mathbb R}U=k_c,2k_c,2k_c`$, je
nachdem, ob $`c`$ vom reellen, komplexen oder quaternionischen Typ ist. Also ist
$`D_U\cong\mathbb R,\mathbb C,\mathbb H`$.
(c) Zu jeder Bahn $`\mathcal O`$ gibt es ein solches $`U`$ im Bild von $`P_{\mathcal O}`$. $`\mathbb R^n`$ selbst ist genau
dann irreduzibel, wenn $`\sigma`$ eine einzige Bahn $`\{c,\sigma(c)\}`$ hat und $`m_c=1`$ (reeller und komplexer Typ) bzw.
$`m_c=2`$ (quaternionischer Typ) ist.

*Beweis.* (a) Für $`T\in D_U`$ sind $`\ker T`$ und $`TU`$ invariante Unterräume, also ist $`T=0`$ oder $`T`$ bijektiv, und dann
ist $`T^{-1}\in D_U`$. $`D_U`$ ist unter der (assoziativen) Komposition abgeschlossen, enthält $`\mathrm{id}_U\neq0`$, und
aus $`ST=0`$ mit $`S\neq0`$ folgt $`T=S^{-1}ST=0`$.

(b) *Invariante Unterräume und reelle Projektionen.* Ist $`U'\subseteq\mathbb R^n`$ invariant, so auch sein orthogonales
Komplement in $`\mathbb R^n`$ (für $`y\perp U'`$, $`u\in U'`$ und $`a\in\mathcal A_{\mathbb R}`$ ist
$`u\cdot ay=a^Tu\cdot y=0`$, da $`a^T\in\mathcal A_{\mathbb R}`$); also vertauscht die Orthogonalprojektion $`e`$ auf $`U'`$,
ein reelles symmetrisches Idempotent, mit $`\mathcal A_{\mathbb R}`$, und nach Lemma 4.18(b) ist
$`e\in(\mathcal A')_{\mathbb R}`$. Umgekehrt ist $`e\mathbb R^n`$ für jede Projektion $`e\in(\mathcal A')_{\mathbb R}`$
invariant. Also ist $`U'\mapsto e`$ eine inklusionserhaltende Bijektion von den invarianten Unterräumen auf die Menge
$`\mathcal P`$ der Projektionen in $`(\mathcal A')_{\mathbb R}`$, und $`U`$ ist genau dann irreduzibel, wenn seine
Projektion $`e`$ ein minimales von null verschiedenes Element von $`\mathcal P`$ ist. Nach Lemma 4.18(a), angewandt auf
das Bild $`e\mathbb C^n`$ (dessen reelle Vektoren $`e\mathbb R^n`$ bilden), ist $`\dim_{\mathbb R}U=\mathrm{rank}\,e`$.

*$`D_U\cong e(\mathcal A')_{\mathbb R}e`$.* Für $`T\in D_U`$ sei $`\hat T`$ die reelle $`n\times n`$-Matrix $`x\mapsto T(ex)`$.
Für $`a\in\mathcal A_{\mathbb R}`$ ist $`\hat Tax=T(eax)=T(aex)=aT(ex)=a\hat Tx`$, also
$`\hat T\in(\mathcal A')_{\mathbb R}`$ (Lemma 4.18(b)) und $`\hat T=e\hat Te`$. Umgekehrt liegt für
$`X\in(\mathcal A')_{\mathbb R}`$ die Einschränkung von $`eXe`$ auf $`U`$ in $`D_U`$. Diese Abbildungen sind zueinander invers
und multiplikativ. Außerdem ist $`e\mathcal A'e`$ unter Konjugation abgeschlossen, und seine reellen Elemente sind die
$`e(\mathcal A')_{\mathbb R}e`$ (ist $`eXe`$ reell, so ist $`eXe=e\mathrm{Re}(X)e`$ mit
$`\mathrm{Re}\,X\in(\mathcal A')_{\mathbb R}`$), also $`\dim_{\mathbb R}D_U=\dim_{\mathbb C}e\mathcal A'e`$ nach Lemma
4.18(a).

*Projektionen in der Zerlegung.* Eine Projektion $`e\in\mathcal A'`$ ist $`e=\sum_cV_c(1\otimes f_c)V_c^*`$ mit
Projektionen $`f_c\in M_{m_c}`$ (Satz 2.6). Dann ist $`\mathrm{rank}\,e=\sum_ck_c\mathrm{rank}\,f_c`$,
$`e\mathcal A'e=\sum_cV_c(1\otimes f_cM_{m_c}f_c)V_c^*`$ hat die komplexe Dimension $`\sum_c(\mathrm{rank}\,f_c)^2`$,
und nach Satz 4.20(d) gilt $`\bar e=e`$ genau dann, wenn $`f_{\sigma(c)}=M_c\bar f_cM_c^*`$ für alle $`c`$.

*Minimale Elemente.* Sei $`e`$ minimal in $`\mathcal P`$. Für jede Bahn $`\mathcal O`$ von $`\sigma`$ ist $`P_{\mathcal O}`$
eine reelle zentrale Projektion ($`\bar P_{\mathcal O}=P_{\mathcal O}`$ nach Satz 4.20(a)) in
$`Z(\mathcal A)\subseteq\mathcal A'`$, die mit $`e`$ vertauscht, also $`eP_{\mathcal O}\in\mathcal P`$ und
$`eP_{\mathcal O}\le e`$; wegen $`e=\sum_{\mathcal O}eP_{\mathcal O}`$ hat genau eine Bahn $`\mathcal O`$ die Eigenschaft
$`eP_{\mathcal O}=e`$, d. h. $`U=e\mathbb R^n`$ liegt im Bild von $`P_{\mathcal O}`$, und $`f_c=0`$ für $`c\notin\mathcal O`$.
Ist $`\mathcal O=\{c,d\}`$ mit $`d\neq c`$, so hat $`f_d=M_c\bar f_cM_c^*`$ den Rang von $`f_c`$; für eine
Rang-eins-Projektion $`g\le f_c`$ ist $`g'=Y+\bar Y`$ mit $`Y=V_c(1\otimes g)V_c^*`$ eine Summe zweier orthogonaler
Projektionen von $`\mathcal A'`$ (in den Sektoren $`c`$ und $`d`$), reell, und $`g'\le e`$ (da $`Y\le e`$ und
$`\bar Y\le\bar e=e`$); also $`g'=e`$, $`\mathrm{rank}\,f_c=\mathrm{rank}\,f_d=1`$,
$`\dim_{\mathbb C}e\mathcal A'e=2`$ und $`\mathrm{rank}\,e=2k_c`$. Ist $`\mathcal O=\{c\}`$, so seien $`f=f_c`$, $`M=M_c`$
und $`\theta v=M\bar v`$, antiunitär auf $`\mathbb C^{m_c}`$ mit $`\theta^2=M\bar M=\varepsilon_c`$ (Satz 4.20(b)). Aus
$`fM=M\bar f`$ folgt $`\theta(fv)=M\bar f\bar v=fM\bar v=f\theta v`$, also ist das Bild $`F`$ von $`f`$ $`\theta`$-invariant.
Ist $`F_1\subseteq F`$ $`\theta`$-invariant mit Projektion $`g`$, so ist $`M\bar gM^*`$ die Projektion auf
$`M\bar F_1=\theta F_1=F_1`$ ($`\theta`$ ist injektiv), also $`g=M\bar gM^*`$, und $`V_c(1\otimes g)V_c^*\in\mathcal P`$ liegt
unterhalb von $`e`$. Wegen der Minimalität hat $`F`$ keine $`\theta`$-invarianten Unterräume außer $`0`$ und $`F`$. Ist
$`\varepsilon_c=1`$, so spannt ein von $`\theta`$ festgehaltener Einheitsvektor in $`F`$ (Lemma 4.19(a)) einen solchen auf,
also $`\mathrm{rank}\,f=1`$, $`\dim_{\mathbb C}e\mathcal A'e=1`$ und $`\mathrm{rank}\,e=k_c`$. Ist
$`\varepsilon_c=-1`$, so ist die Ebene $`\mathrm{span}\{u,\theta u\}`$ eines Einheitsvektors $`u\in F`$ ein solcher
(Lemma 4.19(b)), also $`\mathrm{rank}\,f=2`$, $`\dim_{\mathbb C}e\mathcal A'e=4`$ und $`\mathrm{rank}\,e=2k_c`$. Da
$`\mathbb R`$, $`\mathbb C`$, $`\mathbb H`$ die Dimensionen $`1,2,4`$ haben, liefert (a)
$`D_U\cong\mathbb R,\mathbb C,\mathbb H`$.

(c) $`P_{\mathcal O}\in\mathcal P`$ ist von null verschieden, und ein Element von $`\mathcal P\setminus\{0\}`$ unterhalb
davon mit minimalem Rang ist minimal. $`\mathbb R^n`$ ist genau dann irreduzibel, wenn $`e=1`$ minimal ist; seine Blöcke
sind $`f_c=1_{m_c}`$, also verlangt das nach dem Beweis von (b) eine einzige Bahn und $`m_c=1,1,2`$ in den drei Fällen.
Umgekehrt hat unter diesen Bedingungen ein von null verschiedenes $`e'\in\mathcal P`$ Blöcke $`f'_c\neq0`$ für die
Sektoren der Bahn ($`f'_{\sigma(c)}=M_c\bar f'_cM_c^*`$), also $`f'_c=1`$, falls $`m_c=1`$; im quaternionischen Fall ist das
Bild von $`f'_c`$ $`\theta`$-invariant, also von gerader Dimension (Lemma 4.19(b)), also $`f'_c=1_2`$. Somit ist $`e'=1`$.
$`\square`$

Für eine Gruppe $`G`$, die durch reelle orthogonale Matrizen wirkt, ist $`\mathcal A=\mathrm{span}\,\pi(G)`$ unter
Konjugation abgeschlossen (Lemma 4.18(c)), $`\mathcal A_{\mathbb R}`$ ist der reelle Spann von $`\pi(G)`$ (Lemma 4.18(a):
$`\mathcal A_{\mathbb R}`$ besteht aus den Realteilen der Elemente von $`\mathcal A`$), und $`D_U`$ ist die Kommutante einer
irreduziblen reellen Darstellung von $`G`$: Sie ist $`\mathbb R`$, $`\mathbb C`$ oder $`\mathbb H`$. Nach dem Beweis von (b)
trägt die Komplexifizierung $`U+iU=e\mathbb C^n=V_c(\mathbb C^{k_c}\otimes F)`$ (bei einem Paar zuzüglich ihrer
Konjugierten) in den drei Fällen $`\pi_c`$, $`\pi_c\oplus\pi_{\sigma(c)}`$ bzw. $`\pi_c\oplus\pi_c`$: Sie ist irreduzibel,
die Summe zweier inäquivalenter irreduzibler Darstellungen oder die Summe zweier äquivalenter.

*Rechnung:* Satz 4.22 an den irreduziblen reellen Moduln $`\mathbb R^3`$, $`\mathbb R^2`$, $`\mathbb R^4`$ der drei Beispiele — [Tutorial, Abschnitt 4.7](../../tutorial/de/04_symmetries.md#s-4-7).

<a id="r-4-23"></a>

**Proposition 4.23 (reelle Elemente der Kommutante).** Sei $`X\in\mathcal A'`$ reell, mit Blöcken $`x_c`$. Dann gilt
$`x_{\sigma(c)}=M_c\bar x_cM_c^*`$ für jedes $`c`$. Ist außerdem $`X=X^T`$, so gilt:
(a) Für $`c`$ vom reellen Typ ist $`y_c=W_M^*x_cW_M`$ ($`W_M`$ wie in Satz 4.21(a)) eine reelle symmetrische Matrix, und
$`XO_c=O_c(1_{k_c}\otimes y_c)`$ mit den reellen orthonormalen Spalten $`O_c`$ aus Satz 4.21(a).
(b) Für $`c`$ vom quaternionischen Typ hat jeder Eigenwert von $`x_c`$ gerade Vielfachheit, sodass jeder Eigenwert von $`X`$
auf dem Bild von $`P_c`$ eine durch $`2k_c`$ teilbare Vielfachheit hat (**Kramers-Entartung**).
(c) Für $`c`$ vom komplexen Typ haben $`x_c`$ und $`x_{\sigma(c)}`$ dieselben Eigenwerte mit denselben Vielfachheiten.
Ohne $`X=X^T`$ ist $`y_c`$ in (a) immer noch reell.

*Beweis.* Es ist $`X=\bar X`$, also liefert Satz 4.20(d) $`x_{\sigma(c)}=M_c\bar x_cM_c^*`$. Ist $`X=X^T`$, so ist $`X^*=X`$,
und jedes $`x_c`$ ist hermitesch. (a) Mit $`M_c=W_MW_M^T`$, $`M_c^*=\bar W_MW_M^*`$ und $`W_M^T\bar W_M=1`$ gilt
$`\bar y_c=W_M^T\bar x_c\bar W_M=W_M^TM_c^*x_cM_c\bar W_M=W_M^*x_cW_M=y_c`$, und $`y_c`$ ist hermitesch, also reell
symmetrisch. Ferner ist $`XO_c=XV_c(W_J\otimes W_M)=V_c(W_J\otimes x_cW_M)=O_c(1\otimes y_c)`$. (b) Sei
$`\theta v=M_c\bar v`$, antiunitär mit $`\theta^2=-1`$ (Satz 4.20(b)). Aus $`x_cM_c=M_c\bar x_c`$ folgt
$`x_c\theta v=M_c\bar x_c\bar v=\theta(x_cv)`$; für einen Eigenwert $`\lambda\in\mathbb R`$ von $`x_c`$ erfüllt der
Eigenraum $`E_\lambda`$ also $`\theta E_\lambda\subseteq E_\lambda`$ ($`x_c\theta v=\theta(\lambda v)=\lambda\theta v`$),
und $`\dim E_\lambda`$ ist nach Lemma 4.19(b) gerade. Auf dem Bild von $`P_c`$ wirkt $`X`$ als
$`V_c(1_{k_c}\otimes x_c)V_c^*`$, mit den Eigenwerten von $`x_c`$, jeder $`k_c`$-mal wiederholt. (c) $`\bar x_c`$ hat die
Eigenwerte $`\bar\lambda=\lambda`$ von $`x_c`$, mit den konjugierten Eigenvektoren, und $`x_{\sigma(c)}`$ ist unitär
äquivalent zu $`\bar x_c`$. $`\square`$

Im quaternionischen Fall erfüllt $`y_c=W_M^*x_cW_M`$ mit $`M_c=W_M\Omega_{m_c}W_M^T`$ wie in Lemma 4.19(c) die Gleichung
$`y_c\Omega=\Omega\bar y_c`$ (die Rechnung aus Satz 4.21(b) mit $`M_c`$ anstelle von $`J_c`$), also $`y_c=\mathbb M(h)`$ für
eine quaternionische Matrix $`h`$; ist $`X=X^T`$, so gilt $`h^*=h`$, da $`\mathbb M(h^*)=y_c^*=y_c`$ und $`\mathbb M`$ injektiv
ist.

*Rechnung:* Proposition 4.23 für ein zufälliges reelles symmetrisches Element der Kommutante in jedem der drei Beispiele — [Tutorial, Abschnitt 4.7](../../tutorial/de/04_symmetries.md#s-4-7).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Vertauscht ein Hamiltonoperator mit
> einer Symmetrie und ist er reell in einer Basis, in der die Symmetrie reell ist – die Wirkung einer
> Zeitumkehrsymmetrie, die dort als komplexe Konjugation wirkt –, so sagt Proposition 4.23, wie seine Blöcke aussehen:
> reell symmetrisch in Sektoren vom reellen Typ, komplex hermitesch mit einem Partner gleichen Spektrums in Paaren vom
> komplexen Typ und quaternionisch hermitesch mit zweifach entarteten Niveaus in Sektoren vom quaternionischen Typ,
> wie bei der Zeitumkehr einer ungeraden Anzahl von Spins $`\frac12`$, deren Quadrat $`-1`$ ist (Satz von Kramers). Dysons
> *dreifacher Weg* (*threefold way*; F. J. Dyson, *J. Math. Phys.* 3 (1962) 1199) modelliert die Blöcke eines
> generischen Hamiltonoperators jeder Art durch Zufallsmatrizen der entsprechenden Klasse – das GOE und das GUE aus
> Abschnitt 5.3 sowie das Gaußsche symplektische Ensemble (GSE), das dieser Kurs nicht definiert – und sagt danach die
> Niveaustatistik jedes Blocks voraus.

<a id="s-4-8"></a>

## 4.8 Womit man vergleichen sollte

*In der Bibliothek:* Der Vergleich mit anderer Software — [Tutorial, Abschnitt 4.8](../../tutorial/de/04_symmetries.md#s-4-8).

<a id="x-4"></a>

## Aufgaben

1. Man zeige, dass für die Darstellung von $`S_L`$ auf $`(\mathbb C^d)^{\otimes L}`$
   $`\chi_\pi(\sigma)=d^{\#\text{cycles}(\sigma)}`$ und $`\dim\pi(S_L)'=\binom{d^2+L-1}L`$ gilt. Man prüfe $`20`$
   (Abschnitt 2.4) und $`35`$ ($`=\dim\{S_a\}''`$ für vier Spins).
2. (*Clebsch-Gordan.*) Man zeige $`j_1\otimes j_2\cong\bigoplus_jj`$, $`j=|j_1-j_2|,|j_1-j_2|+1,\dots,j_1+j_2`$, für die
   Spindarstellung $`S^{(j_1)}_a\otimes1+1\otimes S^{(j_2)}_a`$ – die nichtabelsche Fusionsregel von $`SU(2)`$, die Kapitel 7 erwähnt.
3. Sei $`\theta/2\pi`$ irrational. Man zeige, dass $`\{U_\theta\}'=U(1)'`$ für ein einzelnes $`U_\theta=e^{i\theta Q}`$ gilt,
   aber $`\{e^{i\pi\hat N}\}'\cong M_{2^{n-1}}\oplus M_{2^{n-1}}`$ für $`n`$ Qubits.
4. (*Zeeman-Aufspaltung.*) Sei $`H`$ ein Heisenberg-Hamiltonoperator auf $`L\ge1`$ Spins und $`H_h=H+hS_z`$. Man zeige,
   dass $`H_h`$ mit $`S_z`$ und $`\mathbf S^2`$ vertauscht, für $`h\neq0`$ aber nicht mit $`S_x`$, und dass seine Eigenwerte
   $`E+hM`$, $`M=-S,\dots,S`$, für die Eigenwerte $`E`$ von $`H_S`$ sind.
5. (*Die Quaternionengruppe.*) Sei $`Q_8=\{\pm1,\pm I,\pm J,\pm K\}\subseteq\mathbb H`$ mit $`I=(i,0)`$, $`J=(0,1)`$,
   $`K=(0,i)`$ (Abschnitt E.6), und sei $`\pi`$ ihre reguläre Darstellung auf $`\mathbb R^8`$, $`\pi(g)e_h=e_{gh}`$. Man
   zeige, dass $`\mathcal A=\mathrm{span}\,\pi(Q_8)`$ vier Sektoren vom reellen Typ mit $`k=m=1`$ und einen vom
   quaternionischen Typ mit $`k=m=2`$ hat, sodass $`\mathcal A_{\mathbb R}\cong\mathbb R^4\oplus\mathbb H`$.

<details>
<summary>Lösungen</summary>

1. Nach Lemma B.1(c) ist $`\pi(\sigma)`$ eine Permutationsmatrix, die $`|b_1\dots b_L\rangle`$ genau dann festhält, wenn $`b`$
   auf den Zyklen von $`\sigma`$ konstant ist, und es gibt $`d^{\#\text{cycles}(\sigma)}`$ solche $`b`$. Nach Korollar 4.5(c)
   ist $`\dim\pi(S_L)'=\frac1{L!}\sum_\sigma(d^2)^{\#\text{cycles}(\sigma)}=\frac{x(x+1)\cdots(x+L-1)}{L!}`$ mit $`x=d^2`$
   (die Identität $`\sum_\sigma x^{\#\text{cycles}(\sigma)}=x(x+1)\cdots(x+L-1)`$ wird in der Lösung von Kapitel 5,
   Aufgabe 2, bewiesen), und das ist $`\binom{d^2+L-1}L`$ (Proposition B.3 beweist dieselbe Dimension durch Abzählen von
   Multimengen): $`\binom63=20`$ und $`\binom74=35`$. Dass $`35`$ auch $`\dim\{S_a\}''`$ ist, ist Satz 4.14(b) für $`L=4`$:
   $`\{S_a\}''=\pi(S_4)'`$.

2. $`S_a=S^{(j_1)}_a\otimes1+1\otimes S^{(j_2)}_a`$ ist eine Spindarstellung (Abschnitt 4.5), und die Unitäre
   $`x\otimes y\mapsto y\otimes x`$ verflicht sie mit der für $`(j_2,j_1)`$; also sei $`j_1\ge j_2`$. $`S_z`$ hat auf
   $`|m_1\rangle\otimes|m_2\rangle`$ den Eigenwert $`m_1+m_2`$, also liegen seine Eigenwerte in $`j_1+j_2+\mathbb Z`$, und für
   solche $`M`$ ist $`d_M`$ die Anzahl der $`m_2\in\{-j_2,\dots,j_2\}`$ mit $`|M-m_2|\le j_1`$ (dann ist automatisch
   $`m_1=M-m_2\in j_1+\mathbb Z`$); $`d_M=0`$ für $`M\notin j_1+j_2+\mathbb Z`$. Für $`M\ge0`$ lautet die Bedingung
   $`M-j_1\le m_2\le M+j_1`$, und $`M+j_1\ge j_2`$, also läuft $`m_2`$ von $`\max(-j_2,M-j_1)`$ bis $`j_2`$: $`d_M=2j_2+1`$ für
   $`0\le M\le j_1-j_2`$, $`d_M=j_1+j_2-M+1`$ für $`j_1-j_2\le M\le j_1+j_2`$ und $`d_M=0`$ für $`M>j_1+j_2`$. Nach Korollar 4.12
   ist die Vielfachheit von Spin $`j\ge0`$ gleich $`m_j=d_j-d_{j+1}`$. Sie ist $`0`$, falls $`j\notin j_1+j_2+\mathbb Z`$ oder
   $`j>j_1+j_2`$; $`m_{j_1+j_2}=1-0=1`$; für $`j_1-j_2\le j<j_1+j_2`$ (mit $`j\in j_1+j_2+\mathbb Z`$, also $`j+1\le j_1+j_2`$) ist
   $`m_j=(j_1+j_2-j+1)-(j_1+j_2-j)=1`$; und für $`0\le j<j_1-j_2`$ (also $`j+1\le j_1-j_2`$) ist $`m_j=(2j_2+1)-(2j_2+1)=0`$.
   Daher ist $`j_1\otimes j_2\cong\bigoplus_{j=j_1-j_2}^{j_1+j_2}j`$, jeder Spin einmal. Für $`j_2=\frac12`$ ist das
   Satz 4.16(d), der zusätzlich die Koeffizienten liefert.

3. $`U_\theta=\sum_qe^{iq\theta}P_q`$, und $`e^{iq\theta}=e^{iq'\theta}`$ gilt genau dann, wenn $`(q-q')\theta\in2\pi\mathbb Z`$,
   was für irrationales $`\theta/2\pi`$ $`q=q'`$ erzwingt. Also sind die Eigenräume von $`U_\theta`$ die Bilder der $`P_q`$,
   und ein Operator vertauscht genau dann mit der normalen Matrix $`U_\theta`$, wenn er ihre Eigenräume erhält
   (Satz 1.5(d) und Lagrange-Interpolation wie in Abschnitt 2.2), d. h. genau dann, wenn $`X=\sum_qP_qXP_q`$; das ist
   nach Proposition 4.10 $`U(1)'`$. Das Argument verwendet nur, dass die $`e^{iq\theta}`$ für die auftretenden Ladungen
   $`q`$ verschieden sind.
   Für $`\theta=\pi`$ sieht der Eigenwert $`(-1)^w`$ nur die Parität des Hamming-Gewichts: zwei Eigenräume der Dimension
   $`2^{n-1}`$ (für $`n\ge1`$ die Konfigurationen geraden und ungeraden Gewichts; das Umklappen des ersten Bits ist eine
   Bijektion zwischen ihnen), und die Kommutante ist $`M_{2^{n-1}}\oplus M_{2^{n-1}}`$.

4. $`S_z`$ und $`\mathbf S^2`$ vertauschen mit $`H`$ (Proposition 4.15(a)) und mit $`S_z`$, also mit $`H_h`$; und
   $`[H_h,S_x]=h[S_z,S_x]=ihS_y\neq0`$. In der Zerlegung aus Proposition 4.13 ist
   $`V^*S_zV=\bigoplus_SS^{(S)}_z\otimes1_{m_S}`$ und $`V^*HV=\bigoplus_S1_{2S+1}\otimes H_S`$ (Proposition 4.15), also auf
   dem Spin-$`S`$-Summanden $`V^*H_hV=hS^{(S)}_z\otimes1+1\otimes H_S`$, eine Summe vertauschender Operatoren, die durch
   die Produkte von Eigenvektoren diagonalisiert wird, mit den Eigenwerten $`hM+E`$.

5. Nach dem Beweis von Satz E.19 gilt $`I^2=J^2=K^2=-1`$ und $`IJ=-JI=K`$, $`JK=-KJ=I`$, $`KI=-IK=J`$, also ist $`Q_8`$ eine
   Gruppe, und $`h\mapsto gh`$ ist eine Bijektion von $`Q_8`$: Die $`\pi(g)`$ sind Permutationsmatrizen, eine unitäre
   Darstellung durch reelle Matrizen, und $`\mathcal A=\pi(Q_8)''`$ ist unter Konjugation abgeschlossen (Lemma 4.18(c)).
   $`\chi_\pi(g)`$ ist die Anzahl der $`h`$ mit $`gh=h`$, also $`8`$ für $`g=1`$ und sonst $`0`$; nach Korollar 4.5(b) hat daher
   eine irreduzible unitäre Darstellung $`\rho`$ die Vielfachheit
   $`m_\rho=\frac18\overline{\chi_\rho(1)}\cdot8=\dim\rho`$. *Eindimensional:* Bis auf Vorzeichen multiplizieren sich
   $`I,J,K`$ wie die nichtneutralen Elemente von $`\mathbb Z_2\times\mathbb Z_2`$, also ist für $`s,t\in\{\pm1\}`$ die
   Abbildung $`\chi(\pm1)=1`$, $`\chi(\pm I)=s`$, $`\chi(\pm J)=t`$, $`\chi(\pm K)=st`$ ein Homomorphismus; das sind vier
   verschiedene Charaktere, also vier inäquivalente irreduzible Darstellungen, jede mit $`m=1`$. *Zweidimensional:*
   $`\rho(g)=M(g)`$ (Lemma E.15(2)) ist eine unitäre Darstellung ($`M(g)^*M(g)=M(\bar gg)=1`$) mit $`\chi_\rho(\pm1)=\pm2`$
   und sonst $`\chi_\rho=0`$ ($`\mathrm{Tr}\,M(a,b)=2\mathrm{Re}\,a`$), also
   $`\langle\chi_\rho,\chi_\rho\rangle_{Q_8}=\frac{4+4}8=1`$, und $`\rho`$ ist irreduzibel (Korollar 4.5(c)), mit $`m=2`$.
   Für eine Zerlegung $`V`$ wie in Satz 2.6 sind die $`\pi_c\circ\pi`$ irreduzibel und paarweise inäquivalent (Proposition
   2.7), also ist $`V`$ eine Zerlegung wie in Satz 4.3, und jede der fünf obigen Darstellungen ist wegen $`m_\rho\ge1`$ zu
   genau einem $`\pi_c\circ\pi`$ äquivalent, verschiedene zu verschiedenen $`c`$ (Korollar 4.5(b)). Wegen
   $`4\cdot1\cdot1+2\cdot2=8=\sum_ck_cm_c`$ gibt es keine weiteren Sektoren. Nach Satz 4.20(c) dürfen wir eine Zerlegung
   mit $`\pi_c(\pi(g))=\rho_c(g)`$ für diese Darstellungen verwenden (Lemma 4.2(c) liefert unitäre $`u`$ mit
   $`\rho_c=u\,\pi_c(\pi(\cdot))\,u^*`$; man ersetze $`V_c`$ durch $`V_c(u^*\otimes1)`$). Da $`\pi(g)`$ reell ist, liefert
   Satz 4.20(a) $`\pi_{\sigma(c)}(\pi(g))=J_c\overline{\rho_c(g)}J_c^*`$. Für einen Charakter $`\chi`$ ist das $`\chi(g)`$,
   also $`\pi_{\sigma(c)}\circ\pi\cong\pi_c\circ\pi`$, $`\sigma(c)=c`$ (Proposition 2.7: Verschiedene Sektoren tragen
   inäquivalente Darstellungen), und $`J_c\bar J_c=|J_c|^2=1`$: reeller Typ. Für $`\rho=M`$ ist $`\sigma(c)=c`$ (der einzige
   Sektor mit $`k=2`$), und $`\overline{M(g)}=\Omega_2^*M(g)\Omega_2`$ (nach dem Beweis von Satz 4.21(b) gilt
   $`M(q)\Omega_2=\Omega_2\overline{M(q)}`$), also $`M(g)=(J_c\Omega_2^*)M(g)(J_c\Omega_2^*)^*`$. Die Matrizen $`M(1)=1`$,
   $`M(I)=\mathrm{diag}(i,-i)`$, $`M(J)=\Omega_2`$, $`M(K)=\begin{pmatrix}0&-i\\-i&0\end{pmatrix}`$ sind linear
   unabhängig, spannen also $`M_2(\mathbb C)`$ auf; daher ist $`J_c\Omega_2^*=\lambda1`$, $`J_c=\lambda\Omega_2`$ und
   $`J_c\bar J_c=|\lambda|^2\Omega_2^2=-1`$: quaternionischer Typ. Nach Satz 4.21 ist
   $`\mathcal A_{\mathbb R}\cong\mathbb R^4\oplus M_1(\mathbb H)`$, von der Dimension $`8=|Q_8|`$.

</details>

*Rechnung:* die Aufgaben 1–5 numerisch — [Tutorial, Aufgaben](../../tutorial/de/04_symmetries.md#x-4).

**Weiter:** [Kapitel 5](05_random_matrices.md) zieht Zufallsmatrizen in einer Algebra: Haar-verteilte unitäre
Matrizen von $`\prod_cU(k_c)`$, den Weingarten-Kalkül – dessen Fall zweiter Ordnung auf der Kommutante von
$`\{V\otimes V\}`$ beruht, dem kleinsten Fall der Schur-Weyl-Dualität – und die Niveaustatistik Sektor für Sektor.
