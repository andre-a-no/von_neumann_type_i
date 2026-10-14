[English](../04_symmetries.md) · [Español](../es/04_symmetries.md) · [Français](../fr/04_symmetries.md) · **Deutsch** · [中文](../zh/04_symmetries.md) · [日本語](../ja/04_symmetries.md) · [Русский](../ru/04_symmetries.md)

# 4. Symmetrien und Darstellungen

Kapitel 2 hat gezeigt, dass eine Symmetriegruppe eine von-Neumann-Algebra hervorbringt, den Spann der
Symmetrieoperatoren, sowie deren Kommutante, die symmetrischen Observablen (Proposition 2.7), und hat ein
Beispiel von Hand zerlegt. Dieses Kapitel liefert die zugrunde liegende Darstellungstheorie, macht daraus einen
Algorithmus, den die Bibliothek auf jede Liste von Matrizen anwendet, und wendet sie auf die beiden Symmetrien von
Spinsystemen an: $U(1)$ (eine Teilchenzahl oder eine Magnetisierung) und $SU(2)$ (die Drehungen der Spins).

**Sie lernen** die vollständige Reduzibilität und das Schursche Lemma kennen, wie die Zerlegung in irreduzible
Darstellungen die Kommutante und die erzeugte Algebra bestimmt, wie Charaktere Vielfachheiten zählen, warum ein
zufälliges Element der Kommutante die Zerlegung numerisch findet, was Ladungssektoren sind, die irreduziblen
Spindarstellungen und wie $L$ Spins $\frac12$ zerfallen, und warum $SU(2)$-invariante Hamiltonoperatoren viel
kleinere Blöcke haben als $U(1)$-invariante.

Im Folgenden ist $S'$ die Kommutante von $S\subseteq M_n(\mathbb C)$ (Kapitel 2). Ist $S^*=S$, so ist $S''$ die
**von $S$ erzeugte von-Neumann-Algebra**: eine $*$-Algebra, die $S$ enthält, und jede $*$-Algebra
$\mathcal B\supseteq S$ erfüllt $\mathcal B=\mathcal B''\supseteq S''$ (Satz 2.2). $\mathcal A'_h=\{K\in\mathcal A':K=K^*\}$ ist
ein reeller Vektorraum. $\mathbb C^k\otimes\mathbb C^m$ trägt die Kronecker-Reihenfolge aus Kapitel 7
($e_a\otimes f_j$ ist der Basisvektor $am+j$, gezählt ab $0$), sodass $1_k\otimes x$ blockdiagonal mit $k$ Kopien
von $x$ ist. $|G|$ ist die Ordnung einer endlichen Gruppe, und $\hbar=1$.

## 4.1 Darstellungen und das Schursche Lemma

**Definition.** Sei $X$ eine Menge mit einer Involution $x\mapsto x^*$ ($x^{**}=x$). Eine
**Darstellung** von $X$ auf $\mathbb C^k$ ist eine Abbildung $\rho:X\to M_k(\mathbb C)$ mit
$\rho(x^*)=\rho(x)^*$. Ein Unterraum $K$ heißt **invariant**, wenn $\rho(x)K\subseteq K$ für alle $x$ gilt; $\rho$
heißt **irreduzibel**, wenn $k\ge1$ ist und $0,\mathbb C^k$ die einzigen invarianten Unterräume sind. Ein
**Verflechtungsoperator** von $\rho$ in eine Darstellung $\rho'$ auf $\mathbb C^{k'}$ ist eine
$k'\times k$-Matrix $T$ mit $T\rho(x)=\rho'(x)T$ für alle $x$; sie bilden einen Vektorraum
$\operatorname{Hom}(\rho,\rho')$, und es ist $\operatorname{Hom}(\rho,\rho)=\rho(X)'$. $\rho$ und $\rho'$ heißen
**äquivalent**, $\rho\cong\rho'$, wenn es einen invertierbaren Verflechtungsoperator gibt.

Drei Fälle kommen vor: eine Gruppe $G$ mit $g^*=g^{-1}$, wobei wir verlangen, dass $\rho$ ein Homomorphismus
$G\to U(k)$ ist (eine **unitäre Darstellung**; dann ist $\rho(g^{-1})=\rho(g)^{-1}=\rho(g)^*$); eine $*$-Algebra mit
ihrer Adjungierten, wobei $\rho$ ein unitaler $*$-Homomorphismus ist; und Bezeichnungen von Erzeugern mit $x^*=x$,
etwa $\{x,y,z\}$, dargestellt durch Tupel hermitescher Matrizen wie die Spinmatrizen (Abschnitt 4.5). Ist $K$
invariant und bilden die Spalten von $Q\in M_{k\times d}$ eine Orthonormalbasis von $K$, so gilt
$\rho(x)Q=QQ^*\rho(x)Q$, also ist $\rho_K(x)=Q^*\rho(x)Q$ wieder eine Darstellung (ein Homomorphismus, wenn $\rho$
einer ist), die **Einschränkung** von $\rho$ auf $K$.

**Lemma 4.1 (vollständige Reduzibilität).** Sei $\rho$ eine Darstellung von $X$ auf $\mathbb C^n$.
(a) Ist $K$ invariant, so auch $K^\perp$. (b) $\mathbb C^n$ ist eine orthogonale direkte Summe invarianter
Unterräume, auf denen $\rho$ jeweils irreduzibel ist.

*Beweis.* (a) Für $y\in K^\perp$ und $z\in K$ ist $\langle z,\rho(x)y\rangle=\langle\rho(x^*)z,y\rangle=0$, da
$\rho(x^*)z\in K$. (b) Induktion nach $n$. Ein von null verschiedener invarianter Unterraum $K_1$ minimaler
Dimension ist irreduzibel; $K_1^\perp$ ist nach (a) invariant, und die invarianten Unterräume der Einschränkung auf
$K_1^\perp$ sind invariante Unterräume von $\rho$, sodass die Induktionsvoraussetzung auf sie anwendbar ist. $\square$

Für eine endliche Gruppe kostet Unitarität nichts: Sind die $\rho(g)$ nur invertierbar, so erfüllt
$P=\sum_g\rho(g)^*\rho(g)$ die Ungleichung $P\ge\rho(e)^*\rho(e)=1$, also ist $P>0$, und
$\rho(h)^*P\rho(h)=\sum_g\rho(gh)^*\rho(gh)=P$ (man substituiere $g\to gh^{-1}$); daher ist
$(P^{1/2}\rho(h)P^{-1/2})^*(P^{1/2}\rho(h)P^{-1/2})=P^{-1/2}\rho(h)^*P\rho(h)P^{-1/2}=1$, und
$h\mapsto P^{1/2}\rho(h)P^{-1/2}$ ist eine zu $\rho$ äquivalente unitäre Darstellung (Weylscher Trick).

**Lemma 4.2 (Schursches Lemma).** Seien $\rho,\rho'$ irreduzible Darstellungen von $X$ auf
$\mathbb C^k,\mathbb C^{k'}$ und $T\in\operatorname{Hom}(\rho,\rho')$.
(a) $T=0$, oder $T$ ist invertierbar. (b) $\operatorname{Hom}(\rho,\rho)=\mathbb C1$.
(c) $T^*\in\operatorname{Hom}(\rho',\rho)$, und ist $T\neq0$, so gilt $T=\lambda U$ mit $\lambda>0$ und unitärem $U$.
Folglich ist $\dim\operatorname{Hom}(\rho,\rho')$ gleich $1$, falls $\rho\cong\rho'$, und $0$ sonst, und äquivalente
irreduzible Darstellungen sind unitär äquivalent.

*Beweis.* (a) und (b) wurden für Gruppen nach Proposition 2.7 bewiesen, und der Beweis verwendet nur, dass
$\ker T$ invariant unter $\rho$ ist ($T\rho(x)v=\rho'(x)Tv=0$, falls $Tv=0$) und das Bild von $T$ invariant unter
$\rho'$: also $\ker T\in\{0,\mathbb C^k\}$ und $\operatorname{ran}T\in\{0,\mathbb C^{k'}\}$; und für (b) ist
$T-\lambda1$ mit einem Eigenwert $\lambda$ von $T$ ein nicht invertierbarer Verflechtungsoperator, also $0$.
(c) Die Adjungierte von $T\rho(x^*)=\rho'(x^*)T$ ist $\rho(x)T^*=T^*\rho'(x)$. Dann ist
$T^*T\in\operatorname{Hom}(\rho,\rho)$, also $T^*T=\lambda^21$ mit $\lambda>0$, falls $T\neq0$ ($T^*T\ge0$ ist von null
verschieden), und $U=T/\lambda$ ist nach (a) quadratisch mit $U^*U=1$. Sind schließlich
$T_1,T_2\in\operatorname{Hom}(\rho,\rho')$ von null verschieden, so ist
$T_1^{-1}T_2\in\operatorname{Hom}(\rho,\rho)=\mathbb C1$. $\square$

Vertauschen alle $\rho(x)$ miteinander (eine abelsche Gruppe), so ist ein irreduzibles $\rho$ eindimensional: Jedes
$\rho(x)$ liegt in $\operatorname{Hom}(\rho,\rho)=\mathbb C1$, also ist jeder Unterraum invariant.

## 4.2 Die isotypische Zerlegung und die Kommutante

**Satz 4.3 (isotypische Zerlegung).** Sei $\rho$ eine Darstellung von $X$ auf $\mathbb C^n$.
Es gibt irreduzible, paarweise inäquivalente Darstellungen $\rho_1,\dots,\rho_C$ von $X$ auf
$\mathbb C^{k_c}$, ganze Zahlen $m_c\ge1$ mit $\sum_ck_cm_c=n$ und eine unitäre Abbildung
$V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n$ mit
$V^*\rho(x)V=\bigoplus_c\rho_c(x)\otimes1_{m_c}$ für alle $x$. Für jede solche Zerlegung gilt
$$
V^*\rho(X)'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),\qquad
V^*\rho(X)''V=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},
$$
und $m_c=\dim\operatorname{Hom}(\rho_c,\rho)$; für eine irreduzible Darstellung $\sigma$ von $X$ ist
$\operatorname{Hom}(\sigma,\rho)\neq0$ genau dann, wenn $\sigma\cong\rho_c$ für genau ein $c$ gilt.

*Beweis.* *Existenz.* Lemma 4.1 liefert paarweise orthogonale invariante Unterräume $W_1,\dots,W_M$ mit Summe
$\mathbb C^n$ und irreduziblen Einschränkungen $\sigma_i(x)=Q_i^*\rho(x)Q_i$. Man fasse die $W_i$ zu Klassen
äquivalenter Einschränkungen zusammen und wähle in jeder Klasse $c$ eine Referenz $\rho_c=\sigma_{i_0}$. Für jedes
andere $i$ der Klasse liefert Lemma 4.2(c) ein unitäres $U_i$ mit $\sigma_i(x)=U_i\rho_c(x)U_i^*$; ersetzt man $Q_i$
durch $Q_iU_i$, so wird $Q_i^*\rho(x)Q_i=\rho_c(x)$. Man nummeriere die Mitglieder der Klasse $c$ mit
$j=0,\dots,m_c-1$, schreibe $Q_{c,j}$ für ihre Basen und setze $V(e_a\otimes f_j)=Q_{c,j}e_a$ auf dem $c$-ten
Summanden. Da die $W_i$ invariant und orthogonal sind, gilt $Q_{c,j}^*\rho(x)Q_{d,l}=\delta_{cd}\delta_{jl}\rho_c(x)$,
und das ist die behauptete Form.

*Kommutante.* Sei nun $V$ irgendeine Zerlegung wie in der Aussage, und sei $Q_{c,j}$ die $n\times k_c$-Isometrie
$e_a\mapsto V(e_a\otimes f_j)$ (im $c$-ten Summanden); für das soeben konstruierte $V$ sind dies die obigen $Q_{c,j}$. Die
Identität $V^*\rho(x)V=\bigoplus_c\rho_c(x)\otimes1_{m_c}$ besagt $\rho(x)Q_{d,l}=Q_{d,l}\rho_d(x)$, und ihre Adjungierte,
angewandt auf $x^*$, liefert $Q_{c,j}^*\rho(x)=\rho_c(x)Q_{c,j}^*$. Da $V$ unitär ist, vertauscht $Y\in M_n$ genau dann mit allen $\rho(x)$,
wenn jeder Block $Y_{cj,dl}=Q_{c,j}^*YQ_{d,l}$ die Gleichung $Y_{cj,dl}\rho_d(x)=\rho_c(x)Y_{cj,dl}$ erfüllt, d. h. in
$\operatorname{Hom}(\rho_d,\rho_c)$ liegt, und dieser Raum ist $0$ für $c\neq d$ und $\mathbb C1$ für $c=d$ (Lemma 4.2).
Also ist $Y\in\rho(X)'$ genau dann, wenn $Y_{cj,dl}=\delta_{cd}(y_c)_{jl}1_{k_c}$, d. h.
$V^*YV=\bigoplus_c1_{k_c}\otimes y_c$.

*Bikommutante.* Vertauscht $Z$ mit $\rho(X)'$, so vertauscht es mit den Projektionen auf die Summanden (man nehme
$y_c=1$ und $y_d=0$ für $d\neq c$), also ist $V^*ZV=\bigoplus_cZ_c$, und $Z_c$ vertauscht mit
$1_{k_c}\otimes M_{m_c}$; nach Lemma 2.5, Schritt 5 (mit vertauschten Tensorfaktoren), ist
$Z_c\in M_{k_c}\otimes1_{m_c}$. Umgekehrt vertauschen diese Operatoren mit $\rho(X)'$.

*Vielfachheiten.* Für $T\in\operatorname{Hom}(\sigma,\rho)$ mit irreduziblem $\sigma$ zeigt dasselbe Argument
$Q_{c,j}^*T\in\operatorname{Hom}(\sigma,\rho_c)$, und dieser Raum ist $\mathbb Ct_{c}$ mit invertierbarem $t_c$, falls
$\sigma\cong\rho_c$, und $0$ sonst; ferner ist $T=\sum_{c,j}Q_{c,j}Q_{c,j}^*T$. Da die $\rho_c$ paarweise inäquivalent
sind, ist $\sigma$ zu höchstens einer von ihnen äquivalent, etwa zu $\rho_c$ (andernfalls ist $T=0$), und dann ist
$T=\sum_ja_jQ_{c,j}t_c$; umgekehrt liegt jedes $Q_{c,j}t_c$ in $\operatorname{Hom}(\sigma,\rho)$, und diese $m_c$
Abbildungen haben orthogonale Bilder, also ist $\dim\operatorname{Hom}(\sigma,\rho)=m_c$. $\square$

Für eine Gruppe ist $\rho(G)''$ der Spann von $\rho(G)$ (Abschnitt 2.4), und Satz 4.3 ist Proposition 2.7, in der
anderen Richtung gelesen: Dort brachte der Struktursatz die irreduziblen Darstellungen hervor, hier bringen sie die
Struktur hervor; in beiden Fällen ist $\dim\rho(X)'=\sum_cm_c^2$ und $\dim\rho(X)''=\sum_ck_c^2$. Das Bild von
$P_c=\sum_jQ_{c,j}Q_{c,j}^*$, einer minimalen zentralen Projektion von $\rho(X)''$, ist die **isotypische
Komponente** von $\rho_c$.
Nach Satz 4.3 ist jedes hermitesche $H\in\rho(X)'$ von der Form $V(\bigoplus_c1_{k_c}\otimes h_c)V^*$ mit hermiteschen
$h_c\in M_{m_c}$, und $1_{k_c}\otimes h_c$ ist blockdiagonal mit $k_c$ Kopien von $h_c$: Das Spektrum von $H$ ist die
Vereinigung der Spektren der $m_c\times m_c$-Blöcke $h_c$, wobei jeder Eigenwert von $h_c$ $k_c$-mal wiederholt wird.

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Ist $\rho(X)$ eine Symmetrie und
> $H\in\rho(X)'$ der Hamiltonoperator, so sind die eben beschriebenen $k_c$-fach entarteten Niveaus die
> **Multipletts** der Symmetrie, und die Blöcke, die man diagonalisieren muss, haben die Größen $m_c$ der
> Vielfachheiten, nicht die Größen der Multipletts.

Für eine endliche Gruppe lassen sich die Vielfachheiten berechnen, ohne irgendetwas zu zerlegen.

**Satz 4.4 (Schursche Orthogonalitätsrelationen).** Seien $\rho,\rho'$ irreduzible unitäre Darstellungen einer
endlichen Gruppe $G$ auf $\mathbb C^k,\mathbb C^{k'}$. Dann gilt
$$
\frac1{|G|}\sum_{g\in G}\rho(g)_{ij}\,\overline{\rho'(g)_{i'j'}}=
\begin{cases}\delta_{ii'}\delta_{jj'}/k&\text{if }\rho=\rho',\\0&\text{if }\rho\not\cong\rho'.\end{cases}
$$

*Beweis.* Für eine $k\times k'$-Matrix $Y$ sei $T=\frac1{|G|}\sum_g\rho(g)Y\rho'(g)^*$. Für $h\in G$ liefern die
Substitution $g\to hg$ und $\rho'(hg)^*=\rho'(g)^*\rho'(h)^{-1}$ die Gleichung $\rho(h)T=T\rho'(h)$, also ist
$T\in\operatorname{Hom}(\rho',\rho)$. Ist $\rho\not\cong\rho'$, so ist $T=0$; ist $\rho=\rho'$, so ist $T=\lambda1$ mit
$k\lambda=\operatorname{Tr}T=\operatorname{Tr}Y$ (Lemma 4.2). Für $Y=|j\rangle\langle j'|$ ist der Eintrag $T_{ii'}$ die
linke Seite, und $\operatorname{Tr}Y=\delta_{jj'}$. $\square$

Der **Charakter** einer Darstellung ist $\chi_\rho(g)=\operatorname{Tr}\rho(g)$; äquivalente Darstellungen haben gleiche
Charaktere, da $\operatorname{Tr}(T\rho(g)T^{-1})=\operatorname{Tr}\rho(g)$. Wir schreiben
$\langle\varphi,\psi\rangle_G=\frac1{|G|}\sum_g\overline{\varphi(g)}\psi(g)$.

**Korollar 4.5 (Charaktere).** Sei $G$ endlich.
(a) Für irreduzible unitäre $\rho,\rho'$ ist $\langle\chi_\rho,\chi_{\rho'}\rangle_G$ gleich $1$, falls $\rho\cong\rho'$, und $0$ sonst.
(b) Für eine wie in Satz 4.3 zerlegte unitäre Darstellung $\pi$ und ein irreduzibles $\rho$ ist die Vielfachheit
von $\rho$ in $\pi$
$$
m_\rho=\frac1{|G|}\sum_{g\in G}\overline{\chi_\rho(g)}\,\chi_\pi(g)
$$
($m_\rho=m_c$, falls $\rho\cong\rho_c$, und $m_\rho=0$, falls $\rho$ zu keinem $\rho_c$ äquivalent ist).
(c) $\dim\pi(G)'=\sum_cm_c^2=\frac1{|G|}\sum_{g\in G}|\chi_\pi(g)|^2$. Insbesondere ist $\pi$ genau dann irreduzibel,
wenn diese Zahl $1$ ist.

*Beweis.* (a) Man ersetze $\rho'$ durch eine äquivalente Darstellung (gleicher Charakter), summiere Satz 4.4 über
$i=j$, $i'=j'$ und konjugiere: Das Ergebnis ist $\sum_{i,i'}\delta_{ii'}/k=1$, falls $\rho=\rho'$, und $0$, falls
$\rho\not\cong\rho'$. (b) Die Spur ist unitär invariant und $\operatorname{Tr}(a\otimes1_m)=m\operatorname{Tr}a$, also
$\chi_\pi=\sum_cm_c\chi_{\rho_c}$; man wende (a) an. (c) Mit derselben Entwicklung und (a) ist
$\langle\chi_\pi,\chi_\pi\rangle_G=\sum_cm_c^2$, und das ist $\dim\pi(G)'$ nach Satz 4.3. $\square$

Ein zweiter Beweis von (c) zeigt, wie die Bibliothek Kommutanten abtastet. Der **Twirl**
$P(Y)=\frac1{|G|}\sum_g\pi(g)Y\pi(g)^*$ ist idempotent mit Bild $\pi(G)'$ (Kapitel 2, Aufgabe 3) und selbstadjungiert
bezüglich $\langle Y,Z\rangle_{HS}=\operatorname{Tr}(Y^*Z)$ (man substituiere $g\to g^{-1}$), also die orthogonale
Projektion auf $\pi(G)'$, und $\dim\pi(G)'$ ist ihre Spur. Die Abbildung $Y\mapsto AYB$ hat die Spur
$\operatorname{Tr}A\operatorname{Tr}B$ (ihre Diagonaleinträge sind $A_{ii}B_{jj}$), was
$\frac1{|G|}\sum_g|\chi_\pi(g)|^2$ ergibt.

Wir rechnen das Beispiel aus Abschnitt 2.4 mit Charakteren nach. Eine Permutationsmatrix hat als Spur die Anzahl
der Basisvektoren, die sie festhält, und $\pi(\sigma)$ hält $|b_1b_2b_3\rangle$ genau dann fest, wenn $b$ auf den
Zyklen von $\sigma$ konstant ist (Lemma B.1(c)), also ist $\chi_\pi(\sigma)=2^{\#\text{cycles}(\sigma)}$. Drei unitäre
Darstellungen von $S_3$ liegen auf der Hand. Sei $P_\sigma e_i=e_{\sigma(i)}$ die Permutationsdarstellung auf $\mathbb C^3$
(ein Homomorphismus in die Permutationsmatrizen). Die **triviale** Darstellung ist $\sigma\mapsto1$. Das **Signum** ist
$\operatorname{sgn}\sigma=\det P_\sigma$, ein Homomorphismus nach $\{\pm1\}\subseteq U(1)$, da die Determinante
multiplikativ ist und die Determinante einer Permutationsmatrix $\pm1$ ist. Die Gerade $\mathbb C(1,1,1)$ ist invariant,
also auch $E=(1,1,1)^\perp$ (Lemma 4.1(a)); die **Standarddarstellung** ist die Einschränkung von $P$ auf $E$. In einer an
$\mathbb C(1,1,1)\oplus E$ angepassten Orthonormalbasis ist $P_\sigma=1\oplus\operatorname{std}(\sigma)$, also
$\chi_{\rm std}(\sigma)=\operatorname{Tr}P_\sigma-1$, die Anzahl der Fixpunkte minus eins. Die sechs Elemente von $S_3$
sind die Identität, drei Transpositionen und zwei $3$-Zyklen, mit $3,2,1$ Zyklen und $3,1,0$ Fixpunkten; eine
Transposition vertauscht zwei Spalten von $1_3$, ihr Signum ist also $-1$, und ein $3$-Zyklus ist ein Produkt zweier
Transpositionen ($(1\,2\,3)=(1\,3)(1\,2)$), sein Signum ist also $+1$. Das ergibt die Tabelle

| | Identität ($1$) | Transpositionen ($3$) | $3$-Zyklen ($2$) |
|---|---|---|---|
| $\chi_{\rm triv}$ | $1$ | $1$ | $1$ |
| $\operatorname{sgn}$ | $1$ | $-1$ | $1$ |
| $\chi_{\rm std}$ | $2$ | $0$ | $-1$ |
| $\chi_\pi$ | $8$ | $4$ | $2$ |

(also $\operatorname{sgn}\sigma=(-1)^{3-\#\text{cycles}(\sigma)}$; alle Charaktere sind reell). Mit klassenweise
summiertem $\langle\varphi,\psi\rangle_G$:
$\langle\chi_{\rm triv},\chi_{\rm triv}\rangle_G=\frac{1+3+2}6=1$, $\langle\operatorname{sgn},\operatorname{sgn}\rangle_G=\frac{1+3+2}6=1$ und
$\langle\chi_{\rm std},\chi_{\rm std}\rangle_G=\frac{4+0+2}6=1$, also sind die drei Darstellungen nach Korollar 4.5(c)
irreduzibel, und sie sind paarweise inäquivalent, da ihre Charaktere verschieden sind. Korollar 4.5(b) liefert die
Vielfachheiten in $\pi$:
$$
m_{\rm triv}=\tfrac{8+3\cdot4+2\cdot2}6=4,\qquad m_{\rm sgn}=\tfrac{8-3\cdot4+2\cdot2}6=0,\qquad
m_{\rm std}=\tfrac{2\cdot8+0-2\cdot2}6=2 .
$$
Also sind die triviale und die Standarddarstellung äquivalent zu Konstituenten $\rho_c$ von $\pi$ (Satz 4.3) mit
$m_c=4$ und $m_c=2$; sie machen $4\cdot1+2\cdot2=8=\sum_ck_cm_c$ Dimensionen aus, also gibt es keinen weiteren
Konstituenten: $k=(1,2)$, $m=(4,2)$, und $\dim\pi(S_3)'=4^2+2^2=\frac16(64+3\cdot16+2\cdot4)=20$ (Korollar 4.5(c)). Der
Code berechnet die Tabelle neu: Er gewinnt $P_\sigma$ als Einschränkung von $\pi(\sigma)$ auf den invarianten Spann von
$|100\rangle,|010\rangle,|001\rangle$ (auf dem $\pi(\sigma)$ die Positionen der einzigen $1$ permutiert, Lemma B.1(a)), und
$\#\text{cycles}(\sigma)=\dim\ker(P_\sigma-1)$ (die auf den Zyklen konstanten Vektoren).

```python
import math
import torch
import torch_vn_algebra.decompose as dc

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    return (torch.as_tensor(X) - torch.as_tensor(Y)).abs().max().item() < tol

G = dc.permutation_representation(3, 2)                  # (6, 8, 8): pi(sigma) permutes three qubits
chi_pi = G.diagonal(dim1=-2, dim2=-1).sum(-1).real       # characters Tr pi(sigma)
P3 = G[:, [4, 2, 1]][:, :, [4, 2, 1]].real              # on |100>, |010>, |001>: S_3 permuting C^3
cycles = 3 - torch.linalg.matrix_rank(P3 - torch.eye(3, dtype=torch.float64)).double()   # dim ker(P - 1)
assert close(chi_pi, 2.0 ** cycles)
irreps = {'trivial': torch.ones(6, dtype=torch.float64), 'sign': (-1.0) ** (3 - cycles),
          'standard': P3.diagonal(dim1=-2, dim2=-1).sum(-1) - 1}                 # fixed points - 1
chis = torch.stack(list(irreps.values()))
assert close(chis @ chis.T / 6, torch.eye(3, dtype=torch.float64))    # Corollary 4.5(a): all irreducible
mult = chis @ chi_pi / 6                                               # Corollary 4.5(b)
print(dict(zip(irreps, mult.tolist())), (chi_pi ** 2).mean().item())
assert close(mult, torch.tensor([4.0, 0.0, 2.0], dtype=torch.float64))
assert dc.commutant_dimension(G) == 20                                # Corollary 4.5(c): M_4 (+) M_2
```

Das ist die Zerlegung aus Abschnitt 2.4, nun gewonnen aus einer Tabelle von sechs Zahlen pro Darstellung.

## 4.3 Von den Erzeugern zur Algebra: ein randomisierter Algorithmus

Satz 4.3 beschreibt die Zerlegung, berechnet sie aber nicht: Sein Beweis benötigt die irreduziblen invarianten
Unterräume. Murota, Kanno, Kojima und Kojima (*Japan J. Indust. Appl. Math.* 27 (2010) 125–160; siehe auch Maehara
und Murota, *ibid.* 263–293) haben beobachtet, dass die Eigenräume eines *zufälligen Elements der Kommutante* sie
finden. Sei $\mathcal A$ die von gegebenen Matrizen $S_1,\dots,S_J$ und ihren Adjungierten erzeugte
von-Neumann-Algebra, mit einer Zerlegung $V^*\mathcal AV=\bigoplus_cM_{k_c}\otimes1_{m_c}$,
$V^*\mathcal A'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$ (Satz 2.6); nur die Analyse verwendet $V$.

**Lemma 4.6.** Sei $K\in\mathcal A'_h$ mit $V^*KV=\bigoplus_c1_{k_c}\otimes\kappa_c$. Jeder Eigenraum von $K$ ist
$\mathcal A$-invariant. Sind die $N=\sum_cm_c$ Eigenwerte aller $\kappa_c$ zusammen paarweise verschieden, so sind
die Eigenräume von $K$ die $N$ Unterräume $V(\mathbb C^{k_c}\otimes e)$ (im $c$-ten Summanden, mit einem
Eigenvektor $e$ von $\kappa_c$), und $\mathcal A$ wirkt auf jedem von ihnen irreduzibel, als $M_{k_c}$.

*Beweis.* Ist $Kv=\lambda v$ und $a\in\mathcal A$, so ist $Kav=aKv=\lambda av$. Man schreibe
$\kappa_c=\sum_j\mu_{cj}|e_{cj}\rangle\langle e_{cj}|$ mit orthonormalen $e_{cj}$ (Satz 1.3). Dann ist
$V^*KV=\sum_{c,j}\mu_{cj}R_{cj}$, wobei $R_{cj}$ die Projektion auf $\mathbb C^{k_c}\otimes e_{cj}$ ist; die $R_{cj}$
sind orthogonal mit Summe $1$ und die $\mu_{cj}$ verschieden, also sind dies die Spektralprojektionen. Auf
$\mathbb C^{k_c}\otimes e_{cj}$ wirkt das Element $V^*aV=\bigoplus_da_d\otimes1$ als $a_c$, und die $a_c$ füllen
$M_{k_c}$ aus, das außer $0$ und $\mathbb C^{k_c}$ keine invarianten Unterräume hat. $\square$

**Proposition 4.7 (Generizität).** Seien $K$ und $Y$ unabhängige zufällige Elemente von $\mathcal A'_h$ und
$\mathcal A'$, deren Verteilungen Dichten bezüglich des Lebesgue-Maßes auf diesen reellen Vektorräumen haben
(Lemma A.3). Mit Wahrscheinlichkeit eins gilt: (a) Die Voraussetzung von Lemma 4.6 ist erfüllt; (b)
$\Pi_iY\Pi_j\neq0$ für je zwei verschiedene Spektralprojektionen $\Pi_i,\Pi_j$ von $K$, deren Bilder im selben
Summanden liegen.

*Beweis.* (a) Die hermitesche Matrix $\kappa=\bigoplus_c\kappa_c\in M_N$ hängt reell-linear von $K$ ab. Die Gramsche
Determinante $\delta(K)=\det[\operatorname{Tr}\kappa^{i+j}]_{i,j=0}^{N-1}$ von $1,\kappa,\dots,\kappa^{N-1}$ bezüglich
$\langle\cdot,\cdot\rangle_{HS}$ ist ein reelles Polynom in den Koordinaten von $K$ und genau dann von null
verschieden, wenn diese Potenzen linear unabhängig sind. Da $p(\kappa)=U\operatorname{diag}(p(\lambda_l))U^*$ für
$\kappa=U\operatorname{diag}(\lambda)U^*$ gilt, annulliert ein von null verschiedenes Polynom vom Grad $<N$ die Matrix
$\kappa$ genau dann, wenn es auf $\operatorname{spec}\kappa$ verschwindet, und das ist genau dann möglich, wenn
$\kappa$ weniger als $N$ verschiedene Eigenwerte hat (man nehme $\prod_\mu(t-\mu)$ über diese). Also ist
$\delta(K)\neq0$ die Voraussetzung von Lemma 4.6, und $\delta\not\equiv0$: Es ist von null verschieden in
$V(\bigoplus_c1\otimes\operatorname{diag}(\text{verschiedene reelle Zahlen}))V^*\in\mathcal A'_h$. Daher ist
$\mathbb P(\delta(K)=0)=0$ nach Korollar A.7(b).

(b) Die Menge $O=\{\delta\neq0\}\subseteq\mathcal A'_h$ ist offen ($\delta$ ist stetig). Für $K\in O$ hat $K$ genau $N$
verschiedene Eigenwerte $\mu_1(K)<\dots<\mu_N(K)$ (die von $\kappa$, mit Vielfachheiten $k_c$), also sind seine
Spektralprojektionen $\Pi_1(K),\dots,\Pi_N(K)$, nach aufsteigenden Eigenwerten nummeriert, stetig auf $O$
(Proposition A.13). Man halte $l\neq l'$ fest, und sei $\mathcal N_{ll'}$ die Menge der $(K,Y)\in O\times\mathcal A'$ mit
$P_c\Pi_l(K)=\Pi_l(K)$ und $P_c\Pi_{l'}(K)=\Pi_{l'}(K)$ für ein $c$ ($P_c$ die zentralen Projektionen) und
$\Pi_l(K)Y\Pi_{l'}(K)=0$; sie ist durch stetige Gleichungen auf der Borel-Menge $O\times\mathcal A'$ definiert, also eine
Borel-Menge. Man halte ein $K\in O$ fest, für das die beiden Bilder im Summanden $c$ liegen. Nach Lemma 4.6
projizieren $V^*\Pi_lV$ und $V^*\Pi_{l'}V$ auf $\mathbb C^{k_c}\otimes e$ und $\mathbb C^{k_c}\otimes e'$ mit orthonormalen
Eigenvektoren $e,e'$ von $\kappa_c$, und für $V^*YV=\bigoplus_d1\otimes y_d$ ist $\Pi_lY\Pi_{l'}=0$ genau dann, wenn
$\langle e,y_ce'\rangle=0$. Die Abbildung $Y\mapsto\langle e,y_ce'\rangle$ ist reell-linear $\mathcal A'\to\mathbb C$ und
nicht null (ihr Wert an $V(0\oplus\dots\oplus1\otimes|e\rangle\langle e'|\oplus\dots\oplus0)V^*$ ist $1$), also ist ihr Kern
ein echter Unterraum, eine Nullmenge (Korollar A.7(c)). Also ist jeder Schnitt $\{Y:(K,Y)\in\mathcal N_{ll'}\}$ eine
Nullmenge (er ist leer, wenn die Bilder in verschiedenen Summanden liegen), und da $Y$ von $K$ unabhängig ist, gilt
$\mathbb P((K,Y)\in\mathcal N_{ll'})=0$ nach Korollar A.7(d). Das Ereignis aus (b) scheitert nur auf
$\{\delta(K)=0\}\cup\bigcup_{l\neq l'}\{(K,Y)\in\mathcal N_{ll'}\}$, einer endlichen Vereinigung von Ereignissen der
Wahrscheinlichkeit $0$. $\square$

**Lemma 4.8 (Kompressionen; Ausrichten der Kopien).** (a) Sei $W$ ein $\mathcal A$-invarianter Unterraum mit
Orthonormalbasis $Q\in M_{n\times d}$ und $\sigma(a)=Q^*aQ$ die Einschränkung. Dann ist $\sigma(\mathcal A)$ eine
$*$-Algebra, $\{Q^*YQ:Y\in\mathcal A'\}=\sigma(\mathcal A)'$ und $\{Q^*KQ:K\in\mathcal A'_h\}=\sigma(\mathcal A)'_h$. $\mathcal A$
wirkt genau dann irreduzibel auf $W$, wenn $Q^*YQ\in\mathbb C1$ für jedes $Y\in\mathcal A'$. Ist das nicht der Fall und
ist $X$ ein zufälliges Element von $\mathcal A'$ oder von $\mathcal A'_h$ mit Dichte, so ist mit Wahrscheinlichkeit eins
$Q^*XQ\notin\mathbb C1$.

(b) Seien $W_i,W_j$ invariante Unterräume, auf denen $\mathcal A$ irreduzibel
wirkt, mit Orthonormalbasen $Q_i,Q_j$ und Einschränkungen $\sigma_i(a)=Q_i^*aQ_i$, und sei $Y\in\mathcal A'$. Dann ist
$T=Q_i^*YQ_j\in\operatorname{Hom}(\sigma_j,\sigma_i)$. Ist $T\neq0$, so ist $T=\lambda U$ mit $\lambda>0$ und unitärem
$U$, und $Q_j'=Q_jU^*$ erfüllt $Q_j'^*aQ_j'=\sigma_i(a)$ für alle $a\in\mathcal A$. Liegen $W_i,W_j$ in verschiedenen
Summanden ($W_i\subseteq P_c\mathbb C^n$, $W_j\subseteq P_d\mathbb C^n$, $c\neq d$), so ist $T=0$.

*Beweis.* (a) Wie im Beweis von Satz 4.3 ($W$ ist invariant unter $a$ und $a^*$) gilt $aQ=Q\sigma(a)$ und
$Q^*a=\sigma(a)Q^*$ für $a\in\mathcal A$. Daher ist $\sigma(ab)=Q^*aQ\sigma(b)=\sigma(a)\sigma(b)$, $\sigma(a^*)=\sigma(a)^*$ und
$\sigma(1)=Q^*Q=1_d$, also ist $\sigma$ ein unitaler $*$-Homomorphismus und $\sigma(\mathcal A)$ eine $*$-Algebra; und für
$Y\in\mathcal A'$ gilt $Q^*YQ\sigma(a)=Q^*YaQ=Q^*aYQ=\sigma(a)Q^*YQ$. Umgekehrt sei für $T\in\sigma(\mathcal A)'$
$Y=QTQ^*$; dann ist $aY=Q\sigma(a)TQ^*=QT\sigma(a)Q^*=Ya$, also $Y\in\mathcal A'$, mit $Q^*YQ=T$, und $Y$ ist hermitesch,
wenn $T$ es ist. Ist $\sigma$ irreduzibel, so ist $\sigma(\mathcal A)'=\operatorname{Hom}(\sigma,\sigma)=\mathbb C1$
(Lemma 4.2(b)). Ist $E\subseteq\mathbb C^d$ ein $\sigma$-invarianter Unterraum außer $0,\mathbb C^d$, so auch $E^\perp$
(Lemma 4.1(a)), und die Orthogonalprojektion auf $E$ ist ein nicht skalares Element von $\sigma(\mathcal A)'_h$. In
diesem Fall sind $\mathbb C1\cap\sigma(\mathcal A)'$ und $\mathbb R1\cap\sigma(\mathcal A)'_h$ echte reelle Unterräume, ihre
Urbilder unter den surjektiven reell-linearen Abbildungen $X\mapsto Q^*XQ$ sind echte Unterräume von $\mathcal A'$ und
$\mathcal A'_h$, und diese sind Nullmengen (Korollar A.7(c)).

(b) Wie in (a) gilt $aQ_j=Q_j\sigma_j(a)$ und $Q_i^*a=\sigma_i(a)Q_i^*$, also
$T\sigma_j(a)=Q_i^*YaQ_j=Q_i^*aYQ_j=\sigma_i(a)T$. Lemma 4.2 liefert $T=\lambda U$ und
$U\sigma_j(a)U^*=\sigma_i(a)$, und das ist $Q_j'^*aQ_j'$. Schließlich ist $P_c\in\mathcal A\cap\mathcal A'$, $Q_i=P_cQ_i$
und $P_cQ_j=0$, also $T=Q_i^*P_cYQ_j=Q_i^*YP_cQ_j=0$. $\square$

Der Algorithmus (`decompose.wedderburn`, mit den Schnittstellen `from_group` und `from_generators`) lautet nun:

1. *Die Kommutante abtasten.* Für alle Elemente einer endlichen Gruppe nimmt man $Y=P(Z)$ und $K=P(Z'+Z'^*)$ mit dem
   Twirl $P$ und unabhängigen Ginibre-Matrizen $Z,Z'$. $P$ ist eine lineare Projektion von $M_n$ auf $\mathcal A'$ mit
   $P(Z)^*=P(Z^*)$, also sind $Z\mapsto P(Z)$ und $Z'\mapsto P(Z'+Z'^*)$ reell-lineare Surjektionen auf $\mathcal A'$ und
   auf $\mathcal A'_h$ (für $H\in\mathcal A'_h$ ist $P(\frac H2+\frac H2)=H$); eine Ginibre-Matrix ist ein nicht
   ausgearteter zentrierter Gauß-Vektor in $M_n$ (Proposition A.16), also sind $Y$ und $K$ nicht ausgeartete zentrierte
   Gauß-Vektoren in $\mathcal A'$ und $\mathcal A'_h$ und haben Dichten (Satz A.15). Für Erzeuger ist
   $\mathcal A'=\{X:[S_j,X]=[S_j^*,X]=0\}$ der Kern einer linearen Abbildung auf $M_n$ (die Bibliothek verwendet die
   hermiteschen Teile $\frac12(S_j+S_j^*)$, $\frac1{2i}(S_j-S_j^*)$, die denselben komplexen Raum aufspannen wie
   $S_j,S_j^*$ und daher dieselbe Kommutante haben); mit einer Basis $(B_r)_{r\le D}$ davon ist $Y=\sum_rc_rB_r$ und
   $K=Y'+Y'^*$ mit einer unabhängigen Kopie $Y'$ von $Y$ und unabhängigen standard-komplex-Gaußschen $c_r$: Bilder
   eines nicht ausgearteten Gauß-Vektors in $\mathbb C^D$ unter der reell-linearen Bijektion $c\mapsto\sum_rc_rB_r$ auf
   $\mathcal A'$ und der reell-linearen Surjektion $c\mapsto X+X^*$ ($X=\sum_rc_rB_r$) auf $\mathcal A'_h$, also wieder mit
   Dichten. Jedes zufällige Element ist eine neue Ziehung, unabhängig von den früheren; die Bibliothek normiert sie
   auf Norm eins, was keinen Eigenraum und kein Verschwinden eines Blocks ändert.
2. *Aufspalten.* Die Eigenräume von $K$ sind irreduzibel (Lemma 4.6, Proposition 4.7(a)); ihre Dimensionen sind
   die $k_c$.
3. *Gruppieren und ausrichten.* Die Blöcke $Q_i^*YQ_j$ zwischen Eigenräumen gleicher Dimension sind genau innerhalb
   eines Summanden von null verschieden (Proposition 4.7(b), Lemma 4.8(b)); die Klassen der Relation „durch einen von
   null verschiedenen Block verbunden“ ergeben die Sektoren, ihre Größen die $m_c$, und die Polarfaktoren richten die
   Kopien an einer Referenzkopie pro Klasse aus (Lemma 4.8(b); eine an einer ausgerichteten Kopie ausgerichtete Kopie
   ist an der Referenz ausgerichtet).
4. *Zusammensetzen:* $V$ wird in Kronecker-Reihenfolge aufgebaut, sodass $V^*aV=\bigoplus_ca_c\otimes1_{m_c}$.

In exakter Arithmetik, mit exakt getrennten Eigenräumen von $K$, liefert der Algorithmus also mit Wahrscheinlichkeit
eins eine Zerlegung wie in Satz 4.3: Nach Lemma 4.6 enthält der Summand $c$ genau $m_c$ Eigenräume der Dimension
$k_c$, und nach dem Ausrichten ist die Konstruktion der Existenzteil des Beweises von Satz 4.3. In Gleitkommaarithmetik
behandelt die Bibliothek Eigenwerte, die näher als eine Toleranz beieinanderliegen, als gleich, sodass ein berechneter
„Eigenraum“ $W$ eine Summe mehrerer Eigenräume von $K$ sein kann; er ist trotzdem $\mathcal A$-invariant (Lemma 4.6),
und die Bibliothek komprimiert ein neues zufälliges $X\in\mathcal A'_h$ auf ihn. Nach Lemma 4.8(a) ist $Q^*XQ$ ein
Vielfaches von $1$, wann immer $\mathcal A$ irreduzibel auf $W$ wirkt, und andernfalls mit Wahrscheinlichkeit eins
nicht; dann zerlegen seine Eigenräume, die unter $\sigma(\mathcal A)$ invariant sind (Lemma 4.6, angewandt auf
$\sigma(\mathcal A)$), $W$ in kleinere invariante Unterräume, und die Bibliothek wiederholt dies, bis jede Kompression
skalar ist. Was auch immer in Gleitkommaarithmetik geschieht, das Ergebnis wird am Ende zertifiziert:

**Proposition 4.9 (Zertifikat).** Sei $\mathcal A$ von $S_1,\dots,S_J$ und ihren Adjungierten erzeugt und
$\mathcal B=V(\bigoplus_cM_{k_c}\otimes1_{m_c})V^*$ mit unitärem $V$. Liegen alle $S_j$ in $\mathcal B$ und ist
$\dim\mathcal A'=\sum_cm_c^2$, so ist $\mathcal A=\mathcal B$.

*Beweis.* $\mathcal B$ ist eine $*$-Algebra, die alle $S_j,S_j^*$ enthält, also ist $\mathcal A\subseteq\mathcal B$ und
$\mathcal B'\subseteq\mathcal A'$. Satz 4.3, angewandt auf die Darstellung von $\mathcal B$ auf $\mathbb C^n$ (die
Blockabbildungen $b\mapsto b_c$ sind irreduzibel und inäquivalent, da die $c$-te zentrale Projektion von einer auf
$1$ und von den anderen auf $0$ abgebildet wird), liefert $\dim\mathcal B'=\sum_cm_c^2=\dim\mathcal A'$. Also ist
$\mathcal B'=\mathcal A'$ und $\mathcal B=\mathcal B''=\mathcal A''=\mathcal A$. $\square$

`Decomposition.verify()` meldet die relativen Abstände der Erzeuger (und neuer Elemente der Kommutante) von
$\mathcal B$ (bzw. $\mathcal B'$) und vergleicht $\sum_cm_c^2$ mit $\dim\mathcal A'$ (für eine Gruppe aus den
Charakteren, Korollar 4.5(c); für Erzeuger aus dem Kern). Wir wenden es auf $S_3$ und drei Qubits an und schauen
dann hinein: die Eigenräume eines zufälligen $K$ und die Verflechtungsoperatoren eines zufälligen $Y$.

```python
D = dc.from_group(G)                                     # the twirl samples pi(S_3)'
print(D)
assert D.k_factors == (1, 2) and D.multiplicities == (4, 2)            # Section 2.4: k = (1, 2), m = (4, 2)
assert D.verify()['ok']
rho = D.to_algebra(G).matrix                             # (6, 2, 2, 2): the blocks rho_c(sigma), padded
assert close(rho[:, 0, 0, 0], torch.ones(6))                           # sector 1: trivial representation
assert close(rho[:, 1].diagonal(dim1=-2, dim2=-1).sum(-1), chis[2])    # sector 2: the standard one
schur = torch.einsum('gij,gkl->ijkl', rho[:, 1], rho[:, 1].conj()) / 6
assert close(schur.reshape(4, 4), torch.eye(4) / 2)                    # Theorem 4.4 with k = 2
assert close(D.from_algebra(D.to_algebra(G)), G)

Z = torch.randn(8, 8, dtype=dt)
K = dc.twirl(G, Z + Z.mH)                                # a random element of pi(S_3)'_h
w, Q = torch.linalg.eigh(K)
cuts = [0] + [i + 1 for i in range(7) if w[i + 1] - w[i] > 1e-8] + [8]
spaces = [Q[:, a:b] for a, b in zip(cuts[:-1], cuts[1:])]
assert sorted(W.shape[1] for W in spaces) == [1, 1, 1, 1, 2, 2]      # Lemma 4.6: m_c pieces of size k_c
assert all(close(G @ W, W @ (W.mH @ G @ W)) for W in spaces)           # invariant subspaces
Y = dc.twirl(G, torch.randn(8, 8, dtype=dt))
Wa, Wb = [W for W in spaces if W.shape[1] == 2]
s = torch.linalg.svdvals(Wa.mH @ Y @ Wb)                               # Lemma 4.8(b): a multiple of a unitary
assert s[1] > 1e-3 and abs(s[0] - s[1]) < 1e-10
assert all(close(W.mH @ Y @ Wa, 0) for W in spaces if W.shape[1] == 1)   # different summands: zero
```

Die Kommutante bestimmt die Kosten: Eine Twirl-Stichprobe kostet $|G|$ Produkte von $n\times n$-Matrizen. Für Erzeuger
diagonalisiert die Bibliothek zunächst ein hermitesches Element $R=U\operatorname{diag}(\lambda)U^*$ von $\mathcal A$
(eine zufällige reelle Linearkombination der hermiteschen Teile der Erzeuger plus ein Vielfaches eines
symmetrisierten Produkts zweier solcher Kombinationen). Jedes $X\in\mathcal A'$ vertauscht mit $R$, also mit dessen
Spektralprojektionen (Satz 1.5(d)), sodass $U^*XU$ blockdiagonal bezüglich der Eigenräume von $R$ der Dimensionen
$d_\lambda$ ist; dies wird vorab erzwungen, und das Kernproblem hat $N_R=\sum_\lambda d_\lambda^2$ statt $n^2$ Unbekannte.
Das Zusammenfassen mehrerer Eigenwerte vergrößert nur die Blöcke, sodass das reduzierte Problem dieselbe
Lösungsmenge $\mathcal A'$ hat. In der Zerlegung ist $V^*RV=\bigoplus_cr_c\otimes1_{m_c}$ mit hermiteschen
$r_c\in M_{k_c}$; hat $\lambda$ in $r_c$ die Vielfachheit $\mu_c(\lambda)$, so ist $d_\lambda=\sum_c\mu_c(\lambda)m_c$, und da
$\mu^2\ge\mu$ für ganze Zahlen $\mu\ge0$ gilt,
$$
N_R=\sum_\lambda\Big(\sum_c\mu_c(\lambda)m_c\Big)^2\ge\sum_{\lambda,c}\mu_c(\lambda)^2m_c^2\ge\sum_{\lambda,c}\mu_c(\lambda)m_c^2=\sum_ck_cm_c^2,
$$
mit Gleichheit genau dann, wenn die Eigenwerte aller $r_c$ einfach und paarweise verschieden sind. Für den
Gesamtspin von $L=3,4,5$ Spins findet der folgende Code Gleichheit (numerische Beobachtung). Selbst das Minimum ist
groß: Für acht Spins $\frac12$ ($n=256$) ist $\sum_ck_cm_c^2=4900$, und das reduzierte System, eine
$(n^2+N_R)\times N_R$-Matrix, hat $3.5\cdot10^8$ komplexe Einträge. Für Spinketten verwendet Abschnitt 4.6 daher das
Modul `su2`, das die Antwort im Voraus kennt.

```python
for n in (3, 4, 5):
    Dn = dc.from_generators(dc.total_spin(n))
    _, info = dc.commutant_basis(dc.total_spin(n), return_info=True)
    N_min = sum(k * m * m for k, m in zip(Dn.k_factors, Dn.multiplicities))
    assert info['unknowns'] >= N_min                                  # N_R >= sum_c k_c m_c^2
    print(n, info['unknowns'], N_min)                                # observed: equality
```

## 4.4 $U(1)$ und Ladungssektoren

**Proposition 4.10.** Sei $Q$ eine hermitesche Matrix auf $\mathbb C^n$ mit ganzzahligen Eigenwerten $q$ (den
**Ladungen**), Spektralprojektionen $P_q$ und $d_q=\operatorname{rank}P_q$, und sei $U_\theta=e^{i\theta Q}$; wir
schreiben $U(1)'=\{U_\theta:\theta\in\mathbb R\}'$ und $U(1)''$ für deren Kommutante.

1. $\theta\mapsto U_\theta$ ist eine unitäre Darstellung von $U(1)=\mathbb R/2\pi\mathbb Z$, und jede stetige unitäre
   Darstellung von $U(1)$ auf $\mathbb C^n$ hat diese Form.
2. Für $X\in M_n$ sind äquivalent: (i) $XU_\theta=U_\theta X$ für alle $\theta$; (ii) $[X,Q]=0$;
   (iii) $X=\sum_qP_qXP_q$. Daher ist $U(1)'\cong\bigoplus_qM_{d_q}(\mathbb C)$ (die blockdiagonalen Operatoren) und
   $U(1)''=\operatorname{span}\{P_q\}$: Alle $k_c$ sind $1$, und die Vielfachheiten sind die $d_q$.
3. (*Erhaltene Ladung.*) Für hermitesches $H$ gilt $e^{itH}Qe^{-itH}=Q$ für alle $t$ genau dann, wenn $[H,Q]=0$,
   d. h. genau dann, wenn $H\in U(1)'$. Dann ist jedes Sektorgewicht $\operatorname{Tr}(P_q\rho(t))$ von
   $\rho(t)=e^{-itH}\rho e^{itH}$ konstant.

*Beweis.* (1) $U_\theta=\sum_qe^{iq\theta}P_q$ ist unitär, $U_\theta U_{\theta'}=U_{\theta+\theta'}$ (Satz 1.5(c)) und
$U_{2\pi}=1$. Umgekehrt ist für eine stetige Darstellung $(U_\theta)_{\theta\ge0}$ eine stetige Halbgruppe von
Matrizen, also $U_\theta=e^{\theta B}$ für $\theta\ge0$ nach Lemma 8.5 (bewiesen in Kapitel 8; der Beweis benutzt
nur die endliche Dimension), und für $\theta<0$ ist $U_\theta=U_{-\theta}^{-1}=e^{\theta B}$. Ableiten von
$U_\theta^*U_\theta=1$ bei $0$ (Lemma 8.1) ergibt $B^*+B=0$, also $B=iQ$ mit hermiteschem $Q$. Die Exponentialreihe
$e^{i\theta Q}$ stimmt mit dem Funktionalkalkül überein (Proposition 8.2(1) mit $H=-Q$), also erzwingt
$U_{2\pi}=\sum_qe^{2\pi iq}P_q=1$ die Gleichung $e^{2\pi iq}=1$, $q\in\mathbb Z$, für jeden Eigenwert $q$ von $Q$. (2) (i)$\Rightarrow$(ii): Ableiten bei $\theta=0$.
(ii)$\Rightarrow$(iii): $X$ vertauscht mit $P_q=f(Q)$ (Satz 1.5(d)), also $P_qXP_{q'}=XP_qP_{q'}=0$ für
$q\neq q'$. (iii)$\Rightarrow$(i): $XU_\theta=\sum_qe^{iq\theta}P_qXP_q=U_\theta X$. Die blockdiagonalen Operatoren
bilden $\bigoplus_qM_{d_q}\otimes1_1$ (in einer an die Bilder der $P_q$ angepassten Orthonormalbasis), deren Kommutante
$\bigoplus_q1_{d_q}\otimes M_1=\operatorname{span}\{P_q\}$ ist: Ein Operator, der mit ihnen vertauscht, vertauscht mit
jedem $P_q$, ist also blockdiagonal, und jeder Block vertauscht mit $M_{d_q}$, ist also ein Skalar (Lemma 2.5,
Schritt 5, mit $m=1$). (3) $\frac d{dt}e^{itH}Qe^{-itH}=ie^{itH}[H,Q]e^{-itH}$ (Proposition 8.2) verschwindet genau dann für
alle $t$, wenn $[H,Q]=0$. Dann vertauscht $e^{itH}$ mit $P_q$ (Satz 1.5(d), zweimal), und
$\operatorname{Tr}(P_q\rho(t))=\operatorname{Tr}(e^{itH}P_qe^{-itH}\rho)=\operatorname{Tr}(P_q\rho)$. $\square$

Die irreduziblen Darstellungen von $U(1)$ sind also die Charaktere $\theta\mapsto e^{iq\theta}$, $q\in\mathbb Z$
(eindimensional, wie für jede abelsche Gruppe), und die Vielfachheit der Ladung $q$ ist $d_q$. Ladungen addieren
sich unter Tensorprodukten: Sind $u_k$ und $w_l$ Orthonormalbasen aus Eigenvektoren von $Q_1$ und $Q_2$ mit
Eigenwerten $\alpha_k$ und $\beta_l$, so besteht die Orthonormalbasis $u_k\otimes w_l$ aus Eigenvektoren von
$Q_1\otimes1+1\otimes Q_2$ (Eigenwert $\alpha_k+\beta_l$) und von $e^{i\theta Q_1}\otimes e^{i\theta Q_2}$ (Eigenwert
$e^{i\theta\alpha_k}e^{i\theta\beta_l}$), also nach Definition des Funktionalkalküls
$$
e^{i\theta(Q_1\otimes1+1\otimes Q_2)}=e^{i\theta Q_1}\otimes e^{i\theta Q_2},
$$
und dieser Operator wirkt auf dem Bild von $P_{q_1}\otimes P_{q_2}$ mit Ladung $q_1+q_2$. Das ist die Fusionsregel $(q_1,q_2)\mapsto q_1+q_2$ der
fusionierten Tensorprodukte aus Kapitel 7.

Für $n$ Qubits ist die Ladung der **Teilchenzahloperator** $\hat N=\sum_jn_j$, $n_j=|1\rangle\langle1|$ auf Qubit
$j$: $\hat N|b\rangle=|b|\,|b\rangle$ mit dem Hamming-Gewicht $|b|$. Seine Sektoren haben die Dimension
$\binom nw$, also bilden die Operatoren, die das Hamming-Gewicht erhalten,
$$
U(1)'\cong\bigoplus_{w=0}^nM_{\binom nw}(\mathbb C),
$$
und das ist `chains.SpinChain(n).algebra`, das fusionierte Produkt von $n$ Kopien von $M_1\oplus M_1$ mit den
Ladungen $0,1$ (Bit $1$ = Spin aufwärts, $\hat N=S^z+\frac n2$; in `decompose.total_spin` ist Spin aufwärts dagegen
der *erste* Basisvektor von $\mathbb C^2$, was die Basis umbenennt, nicht aber die Sektoren). Jeder
XXZ-Hamiltonoperator $\sum_{\langle ij\rangle}J(s^x_is^x_j+s^y_is^y_j+\Delta s^z_is^z_j)$ erhält $\hat N$: $s^z_is^z_j$ ist diagonal in der
Basis $|b\rangle$, und der Hüpfterm $s^x_is^x_j+s^y_is^y_j=\frac12(s^+_is^-_j+s^-_is^+_j)$ (man setze $s^\pm=s^x\pm is^y$
ein; die Kreuzterme heben sich weg) bildet $|b\rangle$ auf $0$ oder auf die Konfiguration ab, in der ein aufwärts
gerichteter Spin von einem der Plätze $i,j$ auf den anderen gewandert ist, bildet also jeden Eigenraum von $\hat N$
in sich ab.

```python
import functools
from torch_vn_algebra import SpinChain

L = 4
Sx, Sy, Sz = dc.total_spin(L)                              # (16, 16) each, site 0 most significant
N_op = Sz + L / 2 * torch.eye(2 ** L, dtype=dt)            # number of up spins
DU = dc.from_generators(N_op)                              # the algebra generated by U(1)
print(DU)
chain = SpinChain(L, boundary='periodic')                  # U(1)' = (+)_w M_binom(L, w)
assert set(DU.k_factors) == {1} and sorted(DU.multiplicities) == sorted(chain.algebra.k_factors)
assert sorted(chain.algebra.k_factors) == sorted(math.comb(L, w) for w in range(L + 1))

sx, sy, sz = dc.spin_matrices(0.5)
def two_site(a, b, i, j, n=L):                             # a on site i, b on site j of n sites
    ops = [a if k == i else b if k == j else torch.eye(2, dtype=dt) for k in range(n)]
    return functools.reduce(torch.kron, ops)

H = sum(two_site(sx, sx, i, j) + two_site(sy, sy, i, j) + 0.5 * two_site(sz, sz, i, j)
        for i, j in chain.bonds)                           # XXZ ring, Delta = 0.5
assert close(H @ N_op, N_op @ H)                           # [H, N] = 0 ...
assert close(DU.from_commutant(DU.to_commutant(H)), H)     # ... so H is block diagonal (Proposition 4.10)
Hc = chain.xxz(J=1.0, Delta=0.5)                           # the same H, built sector by sector
E_sec = torch.cat([torch.linalg.eigvalsh(Hc.matrix[0, N, :k, :k]) for N, k in enumerate(chain.algebra.k_factors)])
assert close(torch.sort(E_sec).values, torch.linalg.eigvalsh(H))
```

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). $U(1)$-Symmetrie ist die Erhaltung einer Teilchenzahl, einer Ladung oder einer
> Magnetisierung. Proposition 4.10(3) ist das Noether-Theorem in diesem Rahmen: Die Ladung ist genau dann
> erhalten, wenn der Hamiltonoperator symmetrisch ist, und dann verlässt die Dynamik nie einen Sektor (Satz 8.3).
> Kann zudem nichts Kohärenzen zwischen Sektoren messen, so sind die Sektoren Superauswahlsektoren (Kapitel 3).

## 4.5 $SU(2)$: Spin

**Definition.** Eine **Spindarstellung** auf $\mathbb C^d$ ist ein Tripel hermitescher Matrizen $S_x,S_y,S_z$ mit
$$
[S_x,S_y]=iS_z,\qquad[S_y,S_z]=iS_x,\qquad[S_z,S_x]=iS_y,
$$
d. h. eine Darstellung der Bezeichnungen $\{x,y,z\}$ (Abschnitt 4.1) mit diesen Relationen. Das Grundbeispiel ist
$s_a=\frac12\sigma_a$ auf $\mathbb C^2$ (Pauli-Matrizen). Sind $(S_a)$ und $(T_a)$ Spindarstellungen, so auch
$S_a\otimes1+1\otimes T_a$ (die Kreuzterme vertauschen); daher ist der **Gesamtspin** $S_a=\sum_is_a^{(i)}$ von $L$
Spins $\frac12$ (`decompose.total_spin(L)`) eine solche. Die Drehungen $e^{-i\theta\,\mathbf n\cdot\mathbf S}$
($\theta\in\mathbb R$, $\mathbf n\in\mathbb R^3$ ein Einheitsvektor, Funktionalkalkül) vertauschen genau dann für alle
$\theta,\mathbf n$ mit $X$, wenn $X$ mit den drei $S_a$ vertauscht (Ableiten bei $\theta=0$, Proposition 8.2(1);
umgekehrt Satz 1.5(d)), also ist die Kommutante „der $SU(2)$-Wirkung“ gleich $\{S_x,S_y,S_z\}'$, und nach Satz 2.2
erzeugen die Drehungen und die $S_a$ dieselbe Algebra (beide Mengen sind unter Adjunktion abgeschlossen). Für $L$
Spins $\frac12$ sind die Drehungen die $V^{\otimes L}$, $V\in SU(2)$. In der Tat ist $\mathbf n\cdot\mathbf S=\sum_ih^{(i)}$
mit $h=\mathbf n\cdot\mathbf s$, und wie bei den Ladungen in Abschnitt 4.4 (Induktion nach $L$) ist
$e^{-i\theta\,\mathbf n\cdot\mathbf S}=(e^{-i\theta h})^{\otimes L}$ mit $e^{-i\theta h}\in SU(2)$ (seine Eigenwerte sind
$e^{\mp i\theta/2}$, da $h$ spurfrei ist mit $h^2=\frac14$). Umgekehrt ist ein $V\in SU(2)$ normal mit Eigenwerten vom
Betrag $1$ und Produkt $1$, also $V=W\operatorname{diag}(e^{i\varphi},e^{-i\varphi})W^*$ mit unitärem $W$, und $V=e^{iA}$
für die spurfreie hermitesche Matrix $A=W\operatorname{diag}(\varphi,-\varphi)W^*$. Eine spurfreie hermitesche
$2\times2$-Matrix ist eine reelle Linearkombination $A=2\,\mathbf c\cdot\mathbf s$ der Pauli-Matrizen, also ist
$V=e^{-i\theta\,\mathbf n\cdot\mathbf s}$ mit $\theta=2|\mathbf c|$, $\mathbf n=-\mathbf c/|\mathbf c|$ (und $\theta=0$,
falls $\mathbf c=0$). Satz B.6(a),(c) und Korollar B.7 (mit $d=2$) zeigen direkt, dass die $V^{\otimes L}$,
$V\in SU(2)$, die von den $S_a$ erzeugte Algebra aufspannen.

Man setze $S_\pm=S_x\pm iS_y$, sodass $S_-=S_+^*$, und $\mathbf S^2=S_x^2+S_y^2+S_z^2$. Die Relationen ergeben
$$
[S_z,S_\pm]=\pm S_\pm,\qquad[S_+,S_-]=2S_z,\qquad\mathbf S^2=S_-S_++S_z^2+S_z,
$$
und $\mathbf S^2$ vertauscht mit jedem $S_a$ (z. B.
$[\mathbf S^2,S_z]=S_x[S_x,S_z]+[S_x,S_z]S_x+S_y[S_y,S_z]+[S_y,S_z]S_y=-i\{S_x,S_y\}+i\{S_y,S_x\}=0$).

**Satz 4.11 (irreduzible Spindarstellungen).** Für jedes $s\in\{0,\frac12,1,\frac32,\dots\}$ bilden die Matrizen
von `decompose.spin_matrices(s)` – auf $\mathbb C^{2s+1}$ mit der Basis $|m\rangle$, $m=s,s-1,\dots,-s$,
$$
S_z|m\rangle=m|m\rangle,\qquad S_+|m\rangle=\sqrt{s(s+1)-m(m+1)}\;|m+1\rangle,\qquad S_-=S_+^*,
$$
eine irreduzible Spindarstellung (**Spin $s$**) mit $\mathbf S^2=s(s+1)1$. Jede irreduzible Spindarstellung ist zu
genau einer von ihnen unitär äquivalent.

*Beweis.* *Spin $s$ ist eine Spindarstellung.* $[S_z,S_+]=S_+$ ist klar. $S_+S_-|m\rangle=(s(s+1)-(m-1)m)|m\rangle$
und $S_-S_+|m\rangle=(s(s+1)-m(m+1))|m\rangle$, also $[S_+,S_-]=2S_z$; diese beiden Relationen und ihre Adjungierten
sind äquivalent zu den drei definierenden. $\mathbf S^2|m\rangle=(s(s+1)-m(m+1)+m^2+m)|m\rangle=s(s+1)|m\rangle$.
*Irreduzibel:* Ein invariantes $W\neq0$ ist invariant unter $S_z$, dessen Eigenwerte einfach sind, enthält also ein
$|m\rangle$; die Koeffizienten von $S_\pm$ verschwinden nur bei $m=\pm s$, also liefert die Anwendung von $S_\pm$
alle $|m'\rangle$.

*Klassifikation.* Sei $(S_a)$ irreduzibel auf $\mathbb C^d$, $\lambda$ der größte Eigenwert von $S_z$ und $v_0$ ein
normierter Eigenvektor. Es ist $S_zS_+v_0=S_+(S_z+1)v_0=(\lambda+1)S_+v_0$, also $S_+v_0=0$ wegen der Maximalität. Sei
$v_j=S_-^jv_0$; dann ist $S_zv_j=(\lambda-j)v_j$, und per Induktion mit $S_+S_-=S_-S_++2S_z$ gilt
$$
S_+v_j=j(2\lambda-j+1)\,v_{j-1}
$$
($S_+v_j=S_-S_+v_{j-1}+2S_zv_{j-1}=[(j-1)(2\lambda-j+2)+2(\lambda-j+1)]v_{j-1}$). Von null verschiedene $v_j$ sind
Eigenvektoren von $S_z$ zu verschiedenen Eigenwerten, also gibt es einen letzten von null verschiedenen, $v_r$. Der
Spann von $v_0,\dots,v_r$ ist invariant unter $S_z,S_\pm$, also gleich $\mathbb C^d$, und $d=r+1$. Aus
$0=S_+v_{r+1}=(r+1)(2\lambda-r)v_r$ folgt $\lambda=r/2=:s$. Schließlich ist
$\|v_j\|^2=\langle v_{j-1},S_+v_j\rangle=j(2s-j+1)\|v_{j-1}\|^2$, sodass die Einheitsvektoren
$|s-j\rangle=v_j/\|v_j\|$ die Gleichung $S_+|m\rangle=\sqrt{(s-m)(s+m+1)}\,|m+1\rangle$ mit $m=s-j$ erfüllen, also die
obige Formel. Verschiedene $s$ haben verschiedene Dimensionen. $\square$

**Korollar 4.12.** Jede Spindarstellung auf $\mathbb C^n$ ist unitär äquivalent zu
$\bigoplus_s(\text{Spin }s)\otimes1_{m_s}$, und mit $d_M=\dim\ker(S_z-M)$ gilt
$$
m_s=d_s-d_{s+1}.
$$

*Beweis.* Satz 4.3 und Satz 4.11. Jede Kopie von Spin $s'$ enthält genau einen Eigenvektor von $S_z$ zum Eigenwert
$M$, falls $|M|\le s'$ und $s'-M\in\mathbb Z$, und sonst keinen; also ist
$d_M=\sum_{s'\ge|M|,\,s'-M\in\mathbb Z}m_{s'}$, und die Differenz zweier aufeinanderfolgender Werte ist $m_s$. $\square$

**Proposition 4.13 ($L$ Spins $\frac12$).** Der Gesamtspin auf $(\mathbb C^2)^{\otimes L}$ enthält Spin $S$ genau
dann, wenn $S\in\{\frac L2,\frac L2-1,\dots\}$, $S\ge0$, und zwar mit der Vielfachheit
$$
m_S=\binom L{\frac L2-S}-\binom L{\frac L2-S-1}\qquad\Big(\tbinom L{-1}=0\Big),
$$
sodass die von $S_x,S_y,S_z$ erzeugte Algebra $\bigoplus_SM_{2S+1}\otimes1_{m_S}$ ist und ihre Kommutante
$\bigoplus_S1_{2S+1}\otimes M_{m_S}$.

*Beweis.* $S_z=\sum_is_z^{(i)}$ ist diagonal in der Produktbasis, mit Eigenwert $u-\frac L2$ auf einem Produkt aus
$u$ Spins aufwärts und $L-u$ Spins abwärts, also $d_M=\binom L{L/2+M}=\binom L{L/2-M}$, falls
$\frac L2+M\in\{0,\dots,L\}$, und $0$ sonst. Man wende Korollar 4.12 und Satz 4.3 an. Die Werte sind positiv, weil
$\binom L{j+1}/\binom Lj=\frac{L-j}{j+1}>1$ für $j<\frac{L-1}2$. $\square$

Für $L=2$: $\frac12\otimes\frac12=1\oplus0$ (**Clebsch-Gordan**): das Triplett $S_-^j|{\uparrow\uparrow}\rangle$,
$j=0,1,2$, symmetrisch unter der Vertauschung $F$ der beiden Spins, und das Singulett
$\frac1{\sqrt2}(|{\uparrow\downarrow}\rangle-|{\downarrow\uparrow}\rangle)$, antisymmetrisch und von $S_+$ und $S_z$
annulliert. Für $L=4$: $m_2=1$, $m_1=4-1=3$, $m_0=6-4=2$, also
$$
\{S_x,S_y,S_z\}''\cong M_5\oplus M_3\oplus M_1\ \text{ with multiplicities }1,3,2,\qquad
\{S_x,S_y,S_z\}'\cong M_1\oplus M_3\oplus M_2,
$$
mit den Dimensionen $35$ und $14$ – den Zahlen aus Kapitel 10. Die Bibliothek findet dasselbe allein aus den drei
Erzeugern (Sektoren sortiert nach $(k_c,m_c)$):

```python
S4 = dc.total_spin(4)
assert close(S4[0] @ S4[1] - S4[1] @ S4[0], 1j * S4[2])                # [S_x, S_y] = i S_z
D4 = dc.from_generators(S4)
print(D4)
assert D4.k_factors == (1, 3, 5) and D4.multiplicities == (2, 3, 1)     # M_1, M_3, M_5 with m = 2, 3, 1
assert D4.verify()['ok'] and D4.dims == (35, 14)
m_S = {S: math.comb(4, 2 - S) - (math.comb(4, 1 - S) if S < 2 else 0) for S in (0, 1, 2)}
assert [m_S[(k - 1) // 2] for k in D4.k_factors] == list(D4.multiplicities)    # Proposition 4.13
casimir = D4.to_algebra((S4 @ S4).sum(0))                                  # S^2 = S(S+1) on spin S
for c, k in enumerate(D4.k_factors):
    assert close(casimir.matrix[0, c, :k, :k], (k - 1) / 2 * (k + 1) / 2 * torch.eye(k))
```

Welche Operatoren vertauschen mit allen Drehungen? Die Platzpermutationen tun es: Konjugation mit $\pi(\sigma)$
bildet $s^{(i)}_a$ auf $s^{(\sigma(i))}_a$ ab (Lemma B.1(b)) und hält daher $\sum_is^{(i)}_a$ fest. Umgekehrt tut es
nichts anderes:

**Satz 4.14 (Schur-Weyl-Dualität).** Sei $\pi$ die Darstellung von $S_L$, die die Faktoren von
$(\mathbb C^2)^{\otimes L}$ permutiert. (a) $\operatorname{span}\pi(S_L)\subseteq\{S_x,S_y,S_z\}'$. (b) Es gilt
Gleichheit; äquivalent dazu (Satz 2.2) ist die vom Gesamtspin erzeugte Algebra gleich $\pi(S_L)'$.

*Beweis.* Teil (a) wurde soeben bewiesen. Teil (b) ist Korollar B.7, der Fall $d=2$ von Satz B.6, der die
Dualität für $(\mathbb C^d)^{\otimes L}$ und die diagonale Wirkung von $U(d)$ für alle $d$ und $L$ beweist. $\square$

(Zum Weiterlesen: W. Fulton und J. Harris, *Representation Theory*, Springer 1991, Lecture 6; R. Goodman und
N. Wallach, *Symmetry, Representations, and Invariants*, Springer 2009, Abschnitt 9.1.) Für $L\le5$ vergleicht der
Code die beiden Räume numerisch: Wegen (a) genügt es, dass sie dieselbe Dimension haben.

```python
for n in range(2, 6):
    comm = dc.commutant_basis(dc.total_spin(n)).reshape(-1, 4 ** n)          # a basis of {S_a}'
    perm = dc.permutation_representation(n, 2).reshape(-1, 4 ** n)          # all n! permutations
    dims = [torch.linalg.matrix_rank(A).item() for A in (comm, perm, torch.cat([comm, perm]))]
    print(n, dims)
    assert dims[0] == dims[1] == dims[2]                                     # span pi(S_n) = {S_a}'
Dp = dc.from_group(dc.permutation_representation(4, 2))                   # pi(S_4)'' = {S_a}'
print(Dp)
assert sorted(zip(Dp.multiplicities, Dp.k_factors)) == sorted(zip(D4.k_factors, D4.multiplicities))
```

Die ausgegebenen Dimensionen sind $2,5,14,42$, die Catalan-Zahlen $\binom{2L}L/(L+1)$; Korollar B.10 beweist
$\dim\{S_a\}'=\sum_Sm_S^2=\binom{2L}L/(L+1)$ für jedes $L$. Die letzte Zeile zeigt die
Dualität in Aktion: $\pi(S_4)''$ hat Blöcke $k=(1,2,3)$ mit Vielfachheiten $(5,1,3)$ – die Paare $(k_c,m_c)$ des
Gesamtspins mit vertauschten Rollen. Allgemein (Korollar B.9) trägt der Vielfachheitsraum von Spin $S$ eine
irreduzible Darstellung von $S_L$ der Dimension $m_S$, und diese sind paarweise inäquivalent; hier sind die
Dimensionen $1,3,2$ für $S=2,1,0$.

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Ein Spin-$S$-Multiplett besteht für jeden drehinvarianten Hamiltonoperator aus
> $2S+1$ Zuständen gleicher Energie, die sich durch $S^z=M$ unterscheiden. Zwei Spins $\frac12$ bilden ein Triplett
> und ein Singulett; vier Spins bilden ein Quintett, drei Tripletts und zwei Singuletts. Die Vielfachheitsräume
> der beiden Singuletts und der drei Tripletts sind das dekohärenzfreie Qubit und das rauschfreie Qutrit aus
> Kapitel 10.

## 4.6 Spinsysteme

Das **Heisenberg-Modell** auf $L$ Spins $\frac12$ ist $H=\sum_{i<j}J_{ij}\,\mathbf s^{(i)}\cdot\mathbf s^{(j)}$ mit
reellen Kopplungen $J_{ij}$.

**Proposition 4.15.** (a) (*Dirac-Identität*) $\mathbf s^{(i)}\cdot\mathbf s^{(j)}=\frac12F_{ij}-\frac14$, wobei $F_{ij}$
die Spins $i$ und $j$ vertauscht. Daher liegt jeder Heisenberg-Hamiltonoperator in
$\operatorname{span}\pi(S_L)\subseteq\{S_a\}'$, und $V^*HV=\bigoplus_S1_{2S+1}\otimes H_S$ mit Blöcken $H_S\in M_{m_S}$
in der Zerlegung aus Proposition 4.13: Jeder Eigenwert von $H_S$ ist ein $(2S+1)$-fach entartetes Niveau von $H$.
(b) Der Raum der **Vektoren höchsten Gewichts** $\ker S_+\cap\ker(S_z-S)$ hat die Dimension $m_S$, ist invariant
unter $H$, und $H$ wirkt darauf als $H_S$ (bis auf einen unitären Basiswechsel).

*Beweis.* (a) Auf zwei Spins ist
$\mathbf s^{(1)}\cdot\mathbf s^{(2)}=\frac12(\mathbf S^2-\mathbf s^{(1)2}-\mathbf s^{(2)2})=\frac12\mathbf S^2-\frac34$
(Satz 4.11 mit $s=\frac12$). Die Triplettvektoren $S_-^j|{\uparrow\uparrow}\rangle$ sind $|{\uparrow\uparrow}\rangle$,
$|{\downarrow\uparrow}\rangle+|{\uparrow\downarrow}\rangle$, $2|{\downarrow\downarrow}\rangle$, und zusammen mit dem Singulett bilden sie
eine Orthogonalbasis von $\mathbb C^2\otimes\mathbb C^2$. Auf dem Triplett ist $\mathbf S^2=2$ und $F=1$:
$\mathbf S^2=S_-S_++S_z^2+S_z$ ist $0+1+1=2$ auf $|{\uparrow\uparrow}\rangle$, und $\mathbf S^2$ und $F$ vertauschen mit
$S_-$, während $F$ den Vektor $|{\uparrow\uparrow}\rangle$ festhält. Auf dem Singulett ist $\mathbf S^2=0$ (es wird von
$S_+$ und $S_z$ annulliert) und $F=-1$. Beide Seiten stimmen auf diesen vier Vektoren überein. Für Plätze $i,j$ von
$L$ bilde man das Tensorprodukt mit der Identität. Der Rest ist Satz 4.14(a) und Satz 4.3. (b) Wir beweisen mehr:
Für jede wie in Korollar 4.12 zerlegte Spindarstellung und jedes $\mu\in\frac12\mathbb Z$ ist
$\ker S_+\cap\ker(S_z-\mu)$ (das Bild unter der Unitären aus Korollar 4.12 von) $|\mu,\mu\rangle\otimes\mathbb C^{m_\mu}$,
falls $\mu\ge0$, und $0$, falls $\mu<0$. In der Tat ist in dieser Basis ein Vektor mit $S_z=\mu$ von der Form
$\sum_{S'}|S',\mu\rangle\otimes w_{S'}$, summiert über die $S'\ge|\mu|$ mit $S'-\mu\in\mathbb Z$, und nach Satz 4.11 bildet
$S_+$ ihn auf $\sum_{S'}\sqrt{(S'-\mu)(S'+\mu+1)}\,|S',\mu+1\rangle\otimes w_{S'}$ ab, wobei $S'+\mu+1\ge1$; er wird also
genau dann annulliert, wenn $w_{S'}=0$ für alle $S'\neq\mu$. Für das Heisenbergsche $H$ wirkt $1\otimes H_S$ auf
$|S,S\rangle\otimes\mathbb C^{m_S}$ als $H_S$. $\square$

Teil (b) beschreibt, wie das Modul `su2` (Spinketten, Kapitel 9) $H_S$ ohne numerische Zerlegung aufbaut: Es
konstruiert eine Orthonormalbasis der Vektoren höchsten Gewichts, indem es Platz für Platz hinzufügt, mit den
Clebsch-Gordan-Koeffizienten von $j\otimes\frac12$. Wir leiten sie durch Absteigen vom höchsten Gewicht her. Nach
Satz 4.11 und $S_-=S_+^*$ erfüllt Spin $j$
$$
S_\pm|j,m\rangle=\sqrt{(j\mp m)(j\pm m+1)}\;|j,m\pm1\rangle ,
$$
und wir setzen $|j,m\rangle=0$ für $|m|>j$ (verträglich mit der Formel, deren Koeffizient bei $m=\pm j$ verschwindet).

**Satz 4.16 (Clebsch-Gordan-Koeffizienten von $j\otimes\frac12$).** Sei $j\in\{\frac12,1,\frac32,\dots\}$, sei
$S_a=S^{(j)}_a\otimes1+1\otimes s_a$ auf $\mathbb C^{2j+1}\otimes\mathbb C^2$ (Spin $j$ aus Satz 4.11 mit Basis $|j,m\rangle$,
gekoppelt an einen Spin $\frac12$ mit $|{\uparrow}\rangle=|\frac12,\frac12\rangle$, $|{\downarrow}\rangle=|\frac12,-\frac12\rangle$), und
sei $\alpha_M=\sqrt{(j+M+\frac12)/(2j+1)}$, $\beta_M=\sqrt{(j-M+\frac12)/(2j+1)}$. Für $J=j\pm\frac12$ und
$M=J,J-1,\dots,-J$ sei
$$
\begin{aligned}
|j+\tfrac12,M\rangle\!\rangle&=\alpha_M\,|j,M-\tfrac12\rangle\otimes|{\uparrow}\rangle+\beta_M\,|j,M+\tfrac12\rangle\otimes|{\downarrow}\rangle,\\
|j-\tfrac12,M\rangle\!\rangle&=-\beta_M\,|j,M-\tfrac12\rangle\otimes|{\uparrow}\rangle+\alpha_M\,|j,M+\tfrac12\rangle\otimes|{\downarrow}\rangle .
\end{aligned}
$$
(a) Diese $2(2j+1)$ Vektoren bilden eine Orthonormalbasis von $\mathbb C^{2j+1}\otimes\mathbb C^2$.
(b) $|j+\frac12,j+\frac12\rangle\!\rangle=|j,j\rangle\otimes|{\uparrow}\rangle$, und $|j-\frac12,j-\frac12\rangle\!\rangle$ ist bis auf einen
Faktor vom Betrag eins der einzige Einheitsvektor $x$ mit $S_zx=(j-\frac12)x$ und $S_+x=0$; der Faktor wird durch einen
positiven Koeffizienten von $|j,j\rangle\otimes|{\downarrow}\rangle$ festgelegt (die Condon-Shortley-Konvention).
(c) $S_z|J,M\rangle\!\rangle=M|J,M\rangle\!\rangle$ und $S_-|J,M\rangle\!\rangle=\sqrt{(J+M)(J-M+1)}\,|J,M-1\rangle\!\rangle$: Die Zustände
$|J,M\rangle\!\rangle$ entstehen aus $|J,J\rangle\!\rangle$ durch Absteigen und Normieren.
(d) Daher ist $S_a|J,M\rangle\!\rangle=\sum_{M'}(S^{(J)}_a)_{M'M}|J,M'\rangle\!\rangle$ für $a=x,y,z$: Die Unitäre
$|J,M\rangle\mapsto|J,M\rangle\!\rangle$ überführt Spin $(j+\frac12)\oplus$ Spin $(j-\frac12)$ in $S_a$, also ist
$j\otimes\frac12\cong(j+\frac12)\oplus(j-\frac12)$, und die Clebsch-Gordan-Koeffizienten $\langle j,m;\frac12,\sigma|J,M\rangle$
(die Koordinaten von $|J,M\rangle\!\rangle$ in der Produktbasis) sind $\pm\alpha_M,\beta_M$ wie angegeben.

*Beweis.* (a) Ist ein Basisvektor in der Formel $0$, so ist auch sein Koeffizient $0$: $|j,M-\frac12\rangle=0$ nur für
$M=-j-\frac12$, wo $\alpha_M=0$, und $|j,M+\frac12\rangle=0$ nur für $M=j+\frac12$, wo $\beta_M=0$. Also hat jeder Vektor die
Norm $\alpha_M^2+\beta_M^2=1$. Vektoren mit verschiedenem $M$ enthalten verschiedene Produktbasisvektoren, und für
gleiches $M$ ist das Skalarprodukt $-\alpha_M\beta_M+\beta_M\alpha_M=0$. Es sind $(2j+2)+2j=2(2j+1)$ Vektoren.
(b) Die erste Aussage ist $\alpha_{j+1/2}=1$, $\beta_{j+1/2}=0$. Der Eigenraum $S_z=j-\frac12$ wird von
$|j,j-1\rangle\otimes|{\uparrow}\rangle$ und $|j,j\rangle\otimes|{\downarrow}\rangle$ aufgespannt, und mit $s_+|{\uparrow}\rangle=0$,
$s_+|{\downarrow}\rangle=|{\uparrow}\rangle$ gilt
$$
S_+\big(p\,|j,j-1\rangle\otimes|{\uparrow}\rangle+q\,|j,j\rangle\otimes|{\downarrow}\rangle\big)=(p\sqrt{2j}+q)\,|j,j\rangle\otimes|{\uparrow}\rangle ,
$$
was genau dann verschwindet, wenn $q=-\sqrt{2j}\,p$; mit $|p|^2+|q|^2=1$ und $q>0$ (also $p<0$) ist das
$p=-\beta_{j-1/2}=-1/\sqrt{2j+1}$, $q=\alpha_{j-1/2}=\sqrt{2j/(2j+1)}$. (c) $S_z$: Beide Terme haben
$S_z=(M\mp\frac12)\pm\frac12=M$. Für $S_-$ sei
$x_M=p\,|j,M-\frac12\rangle\otimes|{\uparrow}\rangle+q\,|j,M+\frac12\rangle\otimes|{\downarrow}\rangle$; wegen
$s_-|{\uparrow}\rangle=|{\downarrow}\rangle$ und $s_-|{\downarrow}\rangle=0$ ist
$$
S_-x_M=p\sqrt{(j+M-\tfrac12)(j-M+\tfrac32)}\;|j,M-\tfrac32\rangle\otimes|{\uparrow}\rangle
+\Big(p+q\sqrt{(j+M+\tfrac12)(j-M+\tfrac12)}\Big)|j,M-\tfrac12\rangle\otimes|{\downarrow}\rangle .
$$
Wir verwenden $\beta_M\sqrt{(j+M+\frac12)(j-M+\frac12)}=(j-M+\frac12)\alpha_M$ und
$\alpha_M\sqrt{(j+M+\frac12)(j-M+\frac12)}=(j+M+\frac12)\beta_M$.
Für $J=j+\frac12$ ist $(p,q)=(\alpha_M,\beta_M)$ und $c=\sqrt{(J+M)(J-M+1)}=\sqrt{(j+M+\frac12)(j-M+\frac32)}$: Der erste
Koeffizient ist $\sqrt{(j+M+\frac12)(j+M-\frac12)(j-M+\frac32)/(2j+1)}=c\,\alpha_{M-1}$, der zweite
$\alpha_M(j-M+\frac32)=c\,\beta_{M-1}$.
Für $J=j-\frac12$ ist $(p,q)=(-\beta_M,\alpha_M)$ und $c=\sqrt{(j+M-\frac12)(j-M+\frac12)}$: Der erste Koeffizient ist
$-\sqrt{(j-M+\frac12)(j+M-\frac12)(j-M+\frac32)/(2j+1)}=-c\,\beta_{M-1}$, der zweite $\beta_M(j+M-\frac12)=c\,\alpha_{M-1}$.
In beiden Fällen ist $S_-|J,M\rangle\!\rangle=c\,|J,M-1\rangle\!\rangle$ für $M>-J$. Für $M=-J$ verschwinden beide Koeffizienten
($p=\alpha_M=0$, falls $J=j+\frac12$; $j+M-\frac12=0$, falls $J=j-\frac12$), also $S_-|J,-J\rangle\!\rangle=0$. (d) Nach (a) und
(c) sind in der Orthonormalbasis $|J,M\rangle\!\rangle$ die Matrizen von $S_z$ und $S_-$ die von Spin
$(j+\frac12)\oplus$ Spin $(j-\frac12)$ (Satz 4.11 und die obige Formel); dasselbe gilt für die Matrix von $S_+=S_-^*$,
also für die von $S_x=\frac12(S_++S_-)$ und $S_y=\frac1{2i}(S_+-S_-)$. $\square$

Aufgabe 2 erweitert die Zerlegung (aber nicht die Koeffizienten) auf $j_1\otimes j_2$. Für die Bibliothek zählen nur
die beiden Vektoren höchsten Gewichts aus Satz 4.16, und diese funktionieren in jeder Spindarstellung, ob irreduzibel
oder nicht:

**Korollar 4.17 (Vektoren höchsten Gewichts Platz für Platz).** Sei $(T_a)$ eine Spindarstellung auf $E$, sei
$\mathcal H_\mu(E)=\ker T_+\cap\ker(T_z-\mu)$ mit einer Orthonormalbasis $B_\mu$ ($B_\mu=\emptyset$, falls
$\mathcal H_\mu(E)=0$), und sei $S_a=T_a\otimes1+1\otimes s_a$ auf $E\otimes\mathbb C^2$. Für jedes $S\in\{0,\frac12,1,\dots\}$
bilden die Vektoren
$$
v\otimes|{\uparrow}\rangle\quad(v\in B_{S-1/2}),\qquad
\sqrt{\tfrac{2S+1}{2S+2}}\;v\otimes|{\downarrow}\rangle-\tfrac1{\sqrt{(2S+1)(2S+2)}}\;T_-v\otimes|{\uparrow}\rangle\quad(v\in B_{S+1/2})
$$
eine Orthonormalbasis von $\ker S_+\cap\ker(S_z-S)$ (die erste Familie ist für $S=0$ leer). Ausgehend von einem Spin
($E=\mathbb C^2$, $B_{1/2}=\{|{\uparrow}\rangle\}$, alle anderen $B_\mu$ leer) und durch Hinzufügen von $L-1$ Spins erhält man
so Orthonormalbasen der Räume der Vektoren höchsten Gewichts von $L$ Spins $\frac12$, der Dimensionen $m_S$
(Proposition 4.15(b)).

*Beweis.* Sei $a=\sqrt{(2S+1)/(2S+2)}$, $b=-1/\sqrt{(2S+1)(2S+2)}$,
$u(v)=a\,v\otimes|{\downarrow}\rangle+b\,T_-v\otimes|{\uparrow}\rangle$. Für $v\in\mathcal H_{S+1/2}(E)$ ist
$T_+T_-v=(T_-T_++2T_z)v=(2S+1)v$, und $a+(2S+1)b=0$. *Höchstes Gewicht:* $S_z$ ist $S$ auf beiden Familien
($T_z=S\mp\frac12$ auf $v$ und $S-\frac12$ auf $T_-v$); $S_+(v\otimes|{\uparrow}\rangle)=T_+v\otimes|{\uparrow}\rangle=0$, und
$S_+u(v)=a\,T_+v\otimes|{\downarrow}\rangle+a\,v\otimes|{\uparrow}\rangle+b\,T_+T_-v\otimes|{\uparrow}\rangle=(a+(2S+1)b)\,v\otimes|{\uparrow}\rangle=0$.
*Orthonormal:* $\langle v\otimes|{\uparrow}\rangle,u(v')\rangle=b\langle v,T_-v'\rangle=b\langle T_+v,v'\rangle=0$ für
$v\in B_{S-1/2}$, und
$\langle u(v),u(v')\rangle=a^2\langle v,v'\rangle+b^2\langle v,T_+T_-v'\rangle=(a^2+(2S+1)b^2)\langle v,v'\rangle=\langle v,v'\rangle$.
*Aufspannend:* Sei $S_zw=Sw$ und $S_+w=0$, und man schreibe $w=x\otimes|{\uparrow}\rangle+y\otimes|{\downarrow}\rangle$. Dann ist
$T_zx=(S-\frac12)x$, $T_zy=(S+\frac12)y$ und $0=S_+w=(T_+x+y)\otimes|{\uparrow}\rangle+T_+y\otimes|{\downarrow}\rangle$, also
$T_+y=0$, d. h. $y\in\mathcal H_{S+1/2}(E)$, und $T_+x=-y$. Dann erfüllt $x'=x+\frac1{2S+1}T_-y$ die Gleichung
$T_+x'=-y+y=0$, also $x'\in\mathcal H_{S-1/2}(E)$, und das ist $0$ für $S=0$ (Beweis von Proposition 4.15(b)). Da
$b/a=-1/(2S+1)$, liegt $w=x'\otimes|{\uparrow}\rangle+a^{-1}u(y)$, mit linear fortgesetztem $u$, im Spann der beiden
Familien. Für $L$ Spins wende man dies mit $E=(\mathbb C^2)^{\otimes(n-1)}$ und dessen Gesamtspin $T_a$ an,
$n=2,\dots,L$; auf einem Spin ist $\mathcal H_{1/2}=\mathbb C|{\uparrow}\rangle$ und alle anderen $\mathcal H_\mu=0$ (der einzige
Spin ist $\frac12$). $\square$

Die zweite Familie besteht aus den Vektoren $|S,S\rangle\!\rangle$ aus Satz 4.16 mit $j=S+\frac12$, geschrieben innerhalb der
von den $T_-^kv$ aufgespannten Kopie von Spin $j$ (Beweis von Satz 4.11, $|j,j-1\rangle=T_-v/\sqrt{2j}$), und die erste
Familie ist $|j+\frac12,j+\frac12\rangle\!\rangle$ mit $j=S-\frac12$. Bezeichnet man jeden Basisvektor durch seine Folge von
Zwischenspins $S_1=\frac12,S_2,\dots,S_L=S$, $|S_n-S_{n-1}|=\frac12$, so erhält man die Yamanouchi-Basis, die
`su2.SU2Chain` aufbaut (ihr Quellcode verwendet genau diese Koeffizienten; der folgende Code prüft nur die Aussagen von
Korollar 4.17), innerhalb des $U(1)$-Sektors $\hat N=\frac L2+S$ von `SpinChain` (dessen Anordnung der Plätze eine
Umbenennung der Tensorfaktoren ist, die nach Lemma B.1(b) mit dem Gesamtspin vertauscht). Der Code prüft Satz 4.16
an dichten Matrizen und gegen `su2.clebsch_gordan` der Bibliothek, dann die Aussagen von Korollar 4.17 für die Basen
von acht Spins. Außerdem prüfen wir Proposition 4.15(a) an einem Ring aus vier Spins mit `D4` und vergleichen die
Blockgrößen für acht Spins.

```python
from torch_vn_algebra import su2

def ket(j, m):                                            # |j, m> of Theorem 4.11 (zero if |m| > j)
    v = torch.zeros(int(2 * j + 1), dtype=dt)
    if abs(m) <= j:
        v[int(j - m)] = 1
    return v

up, dn = ket(0.5, 0.5), ket(0.5, -0.5)
for j in (0.5, 1.0, 1.5, 2.0):                            # Theorem 4.16
    Sj = dc.spin_matrices(j)
    sp_j = [math.sqrt(j * (j + 1) - (j - k) * (j - k + 1)) for k in range(1, int(2 * j) + 1)]
    assert close(Sj[0] + 1j * Sj[1], torch.diag(torch.tensor(sp_j, dtype=dt), 1))     # S_+ of Theorem 4.11
    assert close(Sj[2], torch.diag(torch.tensor([j - k for k in range(int(2 * j) + 1)], dtype=dt)))
    Sc = torch.stack([torch.kron(a, torch.eye(2, dtype=dt)) + torch.kron(torch.eye(len(a), dtype=dt), b)
                      for a, b in zip(Sj, dc.spin_matrices(0.5))])
    Sp, Sm = Sc[0] + 1j * Sc[1], Sc[0] - 1j * Sc[1]
    cols = []
    for J in (j + 0.5, j - 0.5):
        for M in [J - k for k in range(int(2 * J + 1))]:
            al, be = math.sqrt((j + M + 0.5) / (2 * j + 1)), math.sqrt((j - M + 0.5) / (2 * j + 1))
            p, q = (al, be) if J > j else (-be, al)
            x = p * torch.kron(ket(j, M - 0.5), up) + q * torch.kron(ket(j, M + 0.5), dn)
            assert close(Sc[2] @ x, M * x)
            if M == J:
                assert close(Sp @ x, 0)                               # highest weight
            else:
                assert close(Sm @ cols[-1], math.sqrt((J + M + 1) * (J - M)) * x)   # lowering, (c)
            for s2, e in ((0.5, up), (-0.5, dn)):                     # the library's Racah formula agrees
                cg = su2.clebsch_gordan(j, M - s2, 0.5, s2, J, M) if abs(M - s2) <= j else 0.0
                assert abs(cg - torch.vdot(torch.kron(ket(j, M - s2), e), x).real.item()) < 1e-12
            cols.append(x)
    U = torch.stack(cols, 1)
    assert U.shape == (2 * len(Sj[0]), 2 * len(Sj[0])) and close(U.mH @ U, torch.eye(U.shape[1]))   # (a)

def swap(i, j, n):
    p = list(range(n))
    p[i], p[j] = j, i
    return dc.permutation_representation(n, 2, [p])[0]

ring = [(i, (i + 1) % 4) for i in range(4)]
H4 = sum(swap(i, j, 4) / 2 - torch.eye(16, dtype=dt) / 4 for i, j in ring)        # Dirac's identity
assert close(H4, sum(two_site(a, a, i, j) for a in (sx, sy, sz) for i, j in ring))
h = D4.to_commutant(H4)                                   # blocks H_S of sizes m_S = 2, 3, 1
assert close(D4.from_commutant(h), H4)                    # H4 lies in {S_a}'
levels = torch.cat([torch.linalg.eigvalsh(h.matrix[0, c, :m, :m]).repeat_interleave(k)
                    for c, (k, m) in enumerate(zip(D4.k_factors, D4.multiplicities))])
assert close(torch.sort(levels).values, torch.linalg.eigvalsh(H4))   # level of H_S repeated 2S+1 times

u1, sym = SpinChain(8, boundary='periodic'), su2.SU2Chain(8, boundary='periodic')
for S, m in zip(sym.spins, sym.multiplicities):           # Corollary 4.17: orthonormal highest-weight bases
    B, N = sym.hw_basis(S), int(4 + S)                     # columns in the U(1) sector N = L/2 + S
    assert B.shape == (math.comb(8, N), m) and close(B.T @ B, torch.eye(m, dtype=B.dtype))
    assert N == 8 or close(sym.raise_(B, N), 0)            # S_+ B = 0 (S_z = S holds in sector N)
H_u1, H_su2 = u1.xxz(J=1.0, Delta=1.0), sym.heisenberg(J=1.0)
E_u1 = torch.sort(torch.cat([torch.linalg.eigvalsh(H_u1.matrix[0, N, :k, :k])
                             for N, k in enumerate(u1.algebra.k_factors)])).values
E_su2, S_lev, _ = sym.full_spectrum(H_su2, return_quantum_numbers=True)   # level of block S, 2S+1 times
E_su2, S_lev = E_su2[0], S_lev[0]
assert close(E_u1, E_su2)
assert S_lev[0] == 0 and (S_lev[1:4] == 1).all() and E_su2[4] - E_su2[3] > 1e-3     # singlet, then triplet
assert E_su2[1] - E_su2[0] > 1e-3 and close(E_su2[1:4], E_su2[1].expand(3))       # simple singlet (Theorem 9.9)
assert close(sym.restrict(H_u1).matrix, H_su2.matrix)      # Proposition 4.15(b)
flops = lambda ks: sum(k ** 3 for k in ks)
print(u1.algebra.k_factors, sym.multiplicities, flops(u1.algebra.k_factors), flops(sym.multiplicities))
print(f"ground state energy {E_su2[0].item():.6f}, first excited {E_su2[1].item():.6f} (x3)")
```

Für den Ring aus acht Spins haben die $U(1)$-Blöcke die Größen $\binom8w=1,8,28,56,70,56,28,8,1$ (Abschnitt 4.4), die
$SU(2)$-Blöcke die Größen $m_S=14,28,20,7,1$ für $S=0,\dots,4$ (Proposition 4.13). Der Aufwand der Diagonalisierung,
$\sum k^3$, sinkt von $739162$ auf $33040$, um einen Faktor $22$; die Dimension der Algebra, $\sum k^2$, von
$\sum_w\binom8w\binom8{8-w}=\binom{16}8=12870$ (man wähle $8$ von $16$ Elementen, $w$ davon unter den ersten $8$) auf
$\sum_Sm_S^2=1430=\binom{16}8/9$ (Korollar B.10). Nach Satz 9.9 (der Ring hat gerade Länge und $J>0$) ist der
Grundzustand ein einzelnes Singulett, und $E_0(0)<E_0(1)$. *Numerische Beobachtung* (berechnet vom obigen Code): Die
Grundzustandsenergie ist $-3.651093$, und das erste angeregte Niveau ist ein Triplett bei $-3.128419$.

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Das Heisenberg-Modell beschreibt die
> Austauschwechselwirkung lokalisierter Elektronenspins in magnetischen Isolatoren; $J_{ij}>0$ begünstigt
> antiparallele Nachbarn. Sein Spektrum besteht aus Multipletts, und ein Magnetfeld $h\,S^z$ spaltet jedes Multiplett
> in $2S+1$ äquidistante Niveaus auf (Aufgabe 4). Satz 9.9 (bipartite Kopplungen; für gerade antiferromagnetische
> Ketten und Ringe ein Singulett-Grundzustand) geht auf E. Lieb und D. Mattis zurück (*J. Math. Phys.* 3 (1962) 749).

## 4.7 Womit man vergleichen sollte

Für endliche Gruppen berechnen GAP und SageMath Charaktertafeln exakt, und Korollar 4.5 liefert dann die
Vielfachheiten; RepLAB (MATLAB; D. Rosset, F. Montealegre-Mora, J.-D. Bancal) zerlegt Darstellungen endlicher und
kompakter Gruppen numerisch, ebenfalls mit zufälligen Elementen der Kommutante, und SymbolicWedderburn.jl (Julia)
berechnet Wedderburn-Zerlegungen zur Symmetriereduktion semidefiniter Programme. In der Physik bringen QuSpin und
andere Codes für exakte Diagonalisierung $U(1)$- und Gittersymmetrien in Blockdiagonalform, ITensor und TeNPy
führen abelsche Ladungen in Tensornetzwerken mit, und TensorKit.jl behandelt $SU(2)$ mit Fusionsbäumen. Das Modul
`decompose` der Bibliothek ist ein dichtes Verfahren für kleine Matrizen ($n$ bis zu einigen Hundert); es akzeptiert
jede Liste von Matrizen (Gruppenelemente oder Erzeuger), zertifiziert sein Ergebnis (Proposition 4.9) und gibt die
Algebra und ihre Kommutante als `TypeIAlgebra`-Objekte mit `to_algebra`, `from_algebra`, `to_commutant`,
`from_commutant` zurück, sodass die Kapitel 3 und 5–8 unmittelbar anwendbar sind. Für Spin-$\frac12$-Ketten liefert
`su2.SU2Chain` die $SU(2)$-Blöcke direkt (Korollar 4.17); es akzeptiert bis zu $L=18$ Spins, eine Grenze, die durch
den Speicherbedarf der Basen gesetzt ist.

## Aufgaben

1. Man zeige, dass für die Darstellung von $S_L$ auf $(\mathbb C^d)^{\otimes L}$
   $\chi_\pi(\sigma)=d^{\#\text{cycles}(\sigma)}$ und $\dim\pi(S_L)'=\binom{d^2+L-1}L$ gilt. Man prüfe $20$
   (Abschnitt 2.4) und $35$ ($=\dim\{S_a\}''$ für vier Spins).
2. (*Clebsch-Gordan.*) Man zeige $j_1\otimes j_2\cong\bigoplus_jj$, $j=|j_1-j_2|,|j_1-j_2|+1,\dots,j_1+j_2$, für die
   Spindarstellung $S^{(j_1)}_a\otimes1+1\otimes S^{(j_2)}_a$ – die nichtabelsche Fusionsregel von $SU(2)$, die
   Kapitel 7 erwähnt – und prüfe $1\otimes\frac12$ und $1\otimes1$ mit `from_generators`.
3. Sei $\theta/2\pi$ irrational. Man zeige, dass $\{U_\theta\}'=U(1)'$ für ein einzelnes $U_\theta=e^{i\theta Q}$ gilt,
   aber $\{e^{i\pi\hat N}\}'\cong M_{2^{n-1}}\oplus M_{2^{n-1}}$ für $n$ Qubits.
4. (*Zeeman-Aufspaltung.*) Sei $H$ ein Heisenberg-Hamiltonoperator auf $L\ge1$ Spins und $H_h=H+hS_z$. Man zeige,
   dass $H_h$ mit $S_z$ und $\mathbf S^2$ vertauscht, für $h\neq0$ aber nicht mit $S_x$, und dass seine Eigenwerte
   $E+hM$, $M=-S,\dots,S$, für die Eigenwerte $E$ von $H_S$ sind. Man prüfe das für den Ring aus acht Spins.

<details>
<summary>Lösungen</summary>

1. Nach Lemma B.1(c) ist $\pi(\sigma)$ eine Permutationsmatrix, die $|b_1\dots b_L\rangle$ genau dann festhält, wenn $b$
   auf den Zyklen von $\sigma$ konstant ist, und es gibt $d^{\#\text{cycles}(\sigma)}$ solche $b$. Nach Korollar 4.5(c)
   ist $\dim\pi(S_L)'=\frac1{L!}\sum_\sigma(d^2)^{\#\text{cycles}(\sigma)}=\frac{x(x+1)\cdots(x+L-1)}{L!}$ mit $x=d^2$
   (die Identität $\sum_\sigma x^{\#\text{cycles}(\sigma)}=x(x+1)\cdots(x+L-1)$ wird in der Lösung von Kapitel 5,
   Aufgabe 2, bewiesen), und das ist $\binom{d^2+L-1}L$ (Proposition B.3 beweist dieselbe Dimension durch Abzählen von
   Multimengen): $\binom63=20$ und $\binom74=35$. Dass $35$ auch $\dim\{S_a\}''$ ist, ist Satz 4.14(b) für $L=4$:
   $\{S_a\}''=\pi(S_4)'$.

```python
for n, value in ((3, 20), (4, 35)):
    assert dc.commutant_dimension(dc.permutation_representation(n, 2)) == math.comb(n + 3, n) == value
```

2. $S_a=S^{(j_1)}_a\otimes1+1\otimes S^{(j_2)}_a$ ist eine Spindarstellung (Abschnitt 4.5), und die Unitäre
   $x\otimes y\mapsto y\otimes x$ verflicht sie mit der für $(j_2,j_1)$; also sei $j_1\ge j_2$. $S_z$ hat auf
   $|m_1\rangle\otimes|m_2\rangle$ den Eigenwert $m_1+m_2$, also liegen seine Eigenwerte in $j_1+j_2+\mathbb Z$, und für
   solche $M$ ist $d_M$ die Anzahl der $m_2\in\{-j_2,\dots,j_2\}$ mit $|M-m_2|\le j_1$ (dann ist automatisch
   $m_1=M-m_2\in j_1+\mathbb Z$); $d_M=0$ für $M\notin j_1+j_2+\mathbb Z$. Für $M\ge0$ lautet die Bedingung
   $M-j_1\le m_2\le M+j_1$, und $M+j_1\ge j_2$, also läuft $m_2$ von $\max(-j_2,M-j_1)$ bis $j_2$: $d_M=2j_2+1$ für
   $0\le M\le j_1-j_2$, $d_M=j_1+j_2-M+1$ für $j_1-j_2\le M\le j_1+j_2$ und $d_M=0$ für $M>j_1+j_2$. Nach Korollar 4.12
   ist die Vielfachheit von Spin $j\ge0$ gleich $m_j=d_j-d_{j+1}$. Sie ist $0$, falls $j\notin j_1+j_2+\mathbb Z$ oder
   $j>j_1+j_2$; $m_{j_1+j_2}=1-0=1$; für $j_1-j_2\le j<j_1+j_2$ (mit $j\in j_1+j_2+\mathbb Z$, also $j+1\le j_1+j_2$) ist
   $m_j=(j_1+j_2-j+1)-(j_1+j_2-j)=1$; und für $0\le j<j_1-j_2$ (also $j+1\le j_1-j_2$) ist $m_j=(2j_2+1)-(2j_2+1)=0$.
   Daher ist $j_1\otimes j_2\cong\bigoplus_{j=j_1-j_2}^{j_1+j_2}j$, jeder Spin einmal. Für $j_2=\frac12$ ist das
   Satz 4.16(d), der zusätzlich die Koeffizienten liefert.

```python
def coupled(a, b):                                        # S^(a) (x) 1 + 1 (x) S^(b)
    Sa, Sb = dc.spin_matrices(a), dc.spin_matrices(b)
    ea, eb = torch.eye(Sa.shape[-1], dtype=dt), torch.eye(Sb.shape[-1], dtype=dt)
    return torch.stack([torch.kron(x, eb) + torch.kron(ea, y) for x, y in zip(Sa, Sb)])
for a, b, ks in ((1.0, 0.5, (2, 4)), (1.0, 1.0, (1, 3, 5))):  # 1 (x) 1/2 = 1/2 (+) 3/2, 1 (x) 1 = 0 (+) 1 (+) 2
    D_cg = dc.from_generators(coupled(a, b))
    assert D_cg.k_factors == ks and set(D_cg.multiplicities) == {1}
```

3. $U_\theta=\sum_qe^{iq\theta}P_q$, und $e^{iq\theta}=e^{iq'\theta}$ gilt genau dann, wenn $(q-q')\theta\in2\pi\mathbb Z$,
   was für irrationales $\theta/2\pi$ $q=q'$ erzwingt. Also sind die Eigenräume von $U_\theta$ die Bilder der $P_q$,
   und ein Operator vertauscht genau dann mit der normalen Matrix $U_\theta$, wenn er ihre Eigenräume erhält
   (Satz 1.5(d) und Lagrange-Interpolation wie in Abschnitt 2.2), d. h. genau dann, wenn $X=\sum_qP_qXP_q$; das ist
   nach Proposition 4.10 $U(1)'$. Das Argument verwendet nur, dass die $e^{iq\theta}$ für die auftretenden Ladungen
   $q$ verschieden sind; der Code nimmt $\theta=1$ und die Ladungen $0,\dots,4$ von vier Qubits, für die
   $0<|q-q'|\le4<2\pi$ genügt.
   Für $\theta=\pi$ sieht der Eigenwert $(-1)^w$ nur die Parität des Hamming-Gewichts: zwei Eigenräume der Dimension
   $2^{n-1}$ (für $n\ge1$ die Konfigurationen geraden und ungeraden Gewichts; das Umklappen des ersten Bits ist eine
   Bijektion zwischen ihnen), und die Kommutante ist $M_{2^{n-1}}\oplus M_{2^{n-1}}$.

```python
for theta, m in ((1.0, (1, 1, 4, 4, 6)), (math.pi, (8, 8))):   # theta = 1: e^{iq}, q = 0..4, distinct as 4 < 2 pi
    D_theta = dc.from_generators(torch.linalg.matrix_exp(1j * theta * N_op))
    assert D_theta.multiplicities == m and set(D_theta.k_factors) == {1}
```

4. $S_z$ und $\mathbf S^2$ vertauschen mit $H$ (Proposition 4.15(a)) und mit $S_z$, also mit $H_h$; und
   $[H_h,S_x]=h[S_z,S_x]=ihS_y\neq0$. In der Zerlegung aus Proposition 4.13 ist
   $V^*S_zV=\bigoplus_SS^{(S)}_z\otimes1_{m_S}$ und $V^*HV=\bigoplus_S1_{2S+1}\otimes H_S$ (Proposition 4.15), also auf
   dem Spin-$S$-Summanden $V^*H_hV=hS^{(S)}_z\otimes1+1\otimes H_S$, eine Summe vertauschender Operatoren, die durch
   die Produkte von Eigenvektoren diagonalisiert wird, mit den Eigenwerten $hM+E$.

```python
assert not close(H4 @ Sx - Sx @ H4 + 0.3 * (Sz @ Sx - Sx @ Sz), 0)          # [H4 + 0.3 S_z, S_x] != 0
E_field = sym.full_spectrum(sym.heisenberg(J=1.0, h=0.3))[0]                # E + h M on every multiplet
H_field = u1.xxz(J=1.0, Delta=1.0, h=0.3)
E_u1_field = torch.cat([torch.linalg.eigvalsh(H_field.matrix[0, N, :k, :k])
                        for N, k in enumerate(u1.algebra.k_factors)])
assert close(torch.sort(E_u1_field).values, E_field)
```

</details>

**Weiter:** [Kapitel 5](05_random_matrices.md) zieht Zufallsmatrizen in einer Algebra: Haar-verteilte unitäre
Matrizen von $\prod_cU(k_c)$, den Weingarten-Kalkül – dessen Fall zweiter Ordnung auf der Kommutante von
$\{V\otimes V\}$ beruht, dem kleinsten Fall der Schur-Weyl-Dualität – und die Niveaustatistik Sektor für Sektor.
