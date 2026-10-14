[English](../12_equivariant_learning.md) · [Español](../es/12_equivariant_learning.md) · [Français](../fr/12_equivariant_learning.md) · **Deutsch** · [中文](../zh/12_equivariant_learning.md) · [日本語](../ja/12_equivariant_learning.md) · [Русский](../ru/12_equivariant_learning.md)

# 12. Äquivariantes maschinelles Lernen

Eine Schicht eines neuronalen Netzes heißt **äquivariant**, wenn man dasselbe Ergebnis erhält, ob man die Eingabe mit
einer Symmetrie transformiert und dann die Schicht anwendet oder die Schicht anwendet und dann die Ausgabe
transformiert: Ein Umnummerieren der Knoten eines Graphen nummeriert seine Merkmale auf dieselbe Weise um. Für
lineare Schichten ist das eine Vertauschungsrelation, also bilden die äquivarianten Schichten eine Kommutante, die
nach Kapitel 2 bis auf einen Basiswechsel die Algebra $\bigoplus_cM_{m_c}(\mathbb C)$ ist. Dieses Kapitel zählt ihre
Parameter, baut und trainiert sie mit `torch_vn_algebra.decompose` und nutzt die Algebrenstruktur für unitäre,
positive und spektrale Schichten, die nach Konstruktion äquivariant sind.

**Sie lernen**, warum die äquivarianten linearen Abbildungen zwischen zwei Darstellungen Blockmatrizen sind, deren
Größen durch die Vielfachheiten gegeben sind; wie man sie über Bahnen oder über Charaktere abzählt (DeepSets und die
Bell-Zahlen von Graphennetzen); warum die Mittelung über die Gruppe die Orthogonalprojektion auf die äquivarianten
Abbildungen ist; wie Funktionalkalkül und unitäre Parametrisierungen blockweise auf der Kommutante wirken; und wie
man invariante Vektoren, invariante quadratische Formen und invariante Polynome jedes Grades (Molien-Formel) abzählt.

Die Schichten werden in [Tutorial 12](../../tutorial/de/12_equivariant_learning.md) mit der Bibliothek konstruiert und trainiert; dort werden auch die Ergebnisse dieses Kapitels numerisch überprüft.

Im Folgenden ist $G$ eine endliche Gruppe, und eine **Darstellung** von $G$ auf $\mathbb C^n$ ist ein Homomorphismus
$\pi:G\to M_n(\mathbb C)$ mit unitären Werten (das ist keine Einschränkung: Invertierbare $\pi(g)$ sind unitär
bezüglich des Skalarprodukts $\sum_g\langle\pi(g)x,\pi(g)y\rangle$, Abschnitt 4.1). Ihr **Charakter** ist
$\chi_\pi(g)=\mathrm{Tr}\,\pi(g)$. Irreduzible und äquivalente Darstellungen sowie die **Vielfachheit**
$m_\rho(\pi)$ einer irreduziblen $\rho$ in $\pi$ ($0$, falls $\rho$ nicht auftritt) sind wie in Abschnitt 2.4 erklärt,
wo bewiesen wird, dass $m_\rho(\pi)$ nicht von der Zerlegung abhängt, und wo das Lemma von Schur bewiesen wird;
[Kapitel 4](04_symmetries.md) behandelt beides ausführlich (Abschnitte 4.1–4.2). Auf Rechteckmatrizen ist
$\langle S,T\rangle=\mathrm{Tr}(S^*T)$ das **Hilbert–Schmidt**-Skalarprodukt.

<a id="s-12-1"></a>

## 12.1 Äquivariante lineare Abbildungen

**Definition.** Seien $\pi_1,\pi_2$ Darstellungen von $G$ auf $\mathbb C^{n_1}$, $\mathbb C^{n_2}$. Eine Matrix
$T\in M_{n_2\times n_1}(\mathbb C)$ heißt **äquivariant** (ein **Verflechtungsoperator**), wenn $T\pi_1(g)=\pi_2(g)T$
für alle $g\in G$ gilt. Diese Matrizen bilden einen Vektorraum $\mathrm{Hom}_G(\pi_1,\pi_2)$; für
$\pi_1=\pi_2=\pi$ ist er die Kommutante, $\mathrm{Hom}_G(\pi,\pi)=\pi(G)'$.

Die im maschinellen Lernen häufigsten Darstellungen permutieren Koordinaten. Die Gruppe $G$ wirke auf einer endlichen
Menge $\Omega$; die **Permutationsdarstellung** auf $\mathbb C^\Omega$ ist $\pi(g)e_x=e_{gx}$. Beispiele: $S_n$, die die
$n$ Punkte einer Menge permutiert ($\Omega=[n]=\{1,\dots,n\}$), und $S_n$, die die Knoten eines Graphen umnummeriert
und dabei auf $\Omega=[n]^2$ wirkt, d. h. auf die Adjazenzmatrix durch $X\mapsto P_\sigma XP_\sigma^T$.

<a id="r-12-1"></a>

**Lemma 12.1 (Bahnen).** Seien $\pi_i$ die Permutationsdarstellungen von $G$ auf $\mathbb C^{\Omega_i}$. Dann gilt
$T\in\mathrm{Hom}_G(\pi_1,\pi_2)$ genau dann, wenn $T_{gy,gx}=T_{yx}$ für alle $g$, $y\in\Omega_2$,
$x\in\Omega_1$. Die Indikatormatrizen der Bahnen von $G$ auf $\Omega_2\times\Omega_1$ bilden eine Basis, und
$\dim\mathrm{Hom}_G(\pi_1,\pi_2)$ ist die Anzahl dieser Bahnen.

*Beweis.* $\pi_2(g)T\pi_1(g)^*e_x=\pi_2(g)Te_{g^{-1}x}=\sum_yT_{y,g^{-1}x}e_{gy}$, also ist der $(gy,x)$-Eintrag von
$\pi_2(g)T\pi_1(g)^*$ gleich $T_{y,g^{-1}x}$. Daher ist $T$ genau dann äquivariant, d. h. $\pi_2(g)T\pi_1(g)^*=T$ für
alle $g$, wenn $T_{gy,x}=T_{y,g^{-1}x}$ für alle $g,y,x$ gilt, und das ist (mit $gx$ statt $x$) die angegebene
Bedingung: $T$ ist auf den Bahnen konstant. Die Indikatormatrizen der Bahnen haben disjunkte Träger, sind also linear
unabhängig, und jede auf den Bahnen konstante Funktion ist eine Linearkombination von ihnen. $\square$

*Rechnung:* Lemma 12.1 für $S_4$ auf Knoten- und Kantenmerkmalen, durch Abzählen der Bahnen — [Tutorial, Abschnitt 12.1](../../tutorial/de/12_equivariant_learning.md#s-12-1).

<a id="r-12-2"></a>

**Korollar 12.2.** (a) (*DeepSets.*) Für $S_n$ auf $\mathbb C^n$, $n\ge2$, ist die Kommutante
$\mathrm{span}\{1,J\}$, wobei $J$ die Matrix ist, deren Einträge alle gleich eins sind, und
$\mathbb C^n=\mathbb C\mathbf 1\oplus\mathbf 1^\perp$ ist eine Summe zweier inäquivalenter irreduzibler
Darstellungen, der trivialen und der **Standarddarstellung** der Dimension $n-1$. (b) (*Graphennetze.*) Die
äquivarianten Abbildungen von $\mathbb C^{[n]^p}$ nach $\mathbb C^{[n]^q}$ (Tensoren der Stufe $p$ und $q$ auf $n$
Knoten) werden durch die Partitionen einer $(p+q)$-elementigen Menge in höchstens $n$ Blöcke abgezählt; für
$n\ge p+q$ ist ihre Anzahl die **Bell-Zahl** $B_{p+q}$, die Anzahl aller Partitionen einer $(p+q)$-elementigen Menge
($B_2=2$, $B_3=5$, $B_4=15$).

*Beweis.* (b) Zwei Tupel $x,y\in[n]^{p+q}$ liegen genau dann in derselben Bahn, wenn sie dasselbe
**Gleichheitsmuster** haben, d. h. dieselbe Partition der Positionen $\{1,\dots,p+q\}$ durch
$i\sim j\Leftrightarrow x_i=x_j$. Das Muster ist offensichtlich invariant. Haben umgekehrt $x,y$ dasselbe Muster, so
ist $x_i\mapsto y_i$ eine wohldefinierte Injektion von $\{x_i\}$ in $[n]$, die sich zu einer Permutation $\sigma$ mit
$\sigma x=y$ fortsetzt. Jede Partition mit höchstens $n$ Blöcken ist ein Muster (man gebe den Blöcken verschiedene
Werte). Lemma 12.1 schließt den Beweis ab. Für $n\ge p+q$ ist die Schranke an die Anzahl der Blöcke leer. Zählt man
nach der Anzahl der Blöcke: $B_2=1+1$; $B_3=1+3+1$ (eine Partition einer $3$-elementigen Menge in zwei Blöcke ist
durch ihren einelementigen Block festgelegt); $B_4=1+7+6+1$ (zwei Blöcke: $4$ mit den Größen $3+1$ und $3$ mit den
Größen $2+2$; drei Blöcke: die $\binom42=6$ Wahlen des Paares). (a) Für $p=q=1$ gibt es zwei Bahnen auf $[n]^2$, die Diagonale und ihr
Komplement, mit den Indikatoren $1$ und $J-1$. Also ist $\dim\pi(S_n)'=2=\sum_cm_c^2$ (Satz 2.6, Proposition 2.7),
was zwei Sektoren mit $m_c=1$ erzwingt, d. h. zwei inäquivalente irreduzible Summanden, und $\pi(S_n)'$ wird von
deren Projektionen $P_1,P_2$ aufgespannt. Die Projektion auf die invariante Gerade $\mathbb C\mathbf 1$ liegt in
$\pi(S_n)'$ (Lemma 2.1), ist also $P_1$ oder $P_2$: Die Gerade, die die triviale Darstellung trägt, ist ein Sektor,
und $\mathbf 1^\perp$ ist der andere. $\square$

Die DeepSets-Schicht $x\mapsto\lambda x+\gamma\,(\mathbf 1^Tx)\mathbf 1$ (Zaheer et al. 2017) ist (a); die Basen von
Maron et al. (2019) für Graphennetze sind (b).

*Rechnung:* Korollar 12.2 für $n=4$ ($B_2=2$, $B_3=5$, $B_4=15$) und die Zerlegung der Knoten- und Kantenmerkmale von $S_4$ mit `dc.from_group` (Kantenmerkmale: Sektoren der Dimensionen $k=(1,2,3,3)$ mit den Vielfachheiten $m=(2,1,1,3)$) — [Tutorial, Abschnitt 12.1](../../tutorial/de/12_equivariant_learning.md#s-12-1).

Für zwei Darstellungen reduziert sich alles auf eine einzige Kommutante:

<a id="r-12-3"></a>

**Lemma 12.3 (direkte Summe).** Sei $\pi=\pi_1\oplus\pi_2$ auf $\mathbb C^{n_1}\oplus\mathbb C^{n_2}$, mit den
isometrischen Inklusionen $Q_i$ und den Projektionen $P_i=Q_iQ_i^*$. Dann gilt $P_1,P_2\in\pi(G)'$,
$\mathrm{Hom}_G(\pi_1,\pi_2)=Q_2^*\,\pi(G)'\,Q_1$, und $T\mapsto Q_2TQ_1^*$ ist eine lineare Bijektion von
$\mathrm{Hom}_G(\pi_1,\pi_2)$ auf $P_2\,\pi(G)'P_1$.

*Beweis.* $\pi(g)Q_i=Q_i\pi_i(g)$ und $Q_i^*\pi(g)=\pi_i(g)Q_i^*$, weil $\pi(g)$ blockdiagonal ist; daher gilt
$P_i\pi(g)=Q_i\pi_i(g)Q_i^*=\pi(g)P_i$. Ist $X\in\pi(G)'$, so gilt
$Q_2^*XQ_1\pi_1(g)=Q_2^*X\pi(g)Q_1=Q_2^*\pi(g)XQ_1=\pi_2(g)Q_2^*XQ_1$. Ist $T$ äquivariant, so gilt
$\pi(g)Q_2TQ_1^*=Q_2\pi_2(g)TQ_1^*=Q_2T\pi_1(g)Q_1^*=Q_2TQ_1^*\pi(g)$, also $Q_2TQ_1^*\in\pi(G)'$; diese Matrix ist
gleich $P_2(Q_2TQ_1^*)P_1$, und $Q_2^*(Q_2TQ_1^*)Q_1=T$. Umgekehrt liegt jedes $Y=P_2XP_1$ mit $X\in\pi(G)'$ in
$\pi(G)'$, und $Y=Q_2(Q_2^*YQ_1)Q_1^*$. $\square$

*Rechnung:* Die Ecken $Q_2^*\,\pi(G)'\,Q_1$ für Knoten- und Kantenmerkmale — [Tutorial, Abschnitt 12.2](../../tutorial/de/12_equivariant_learning.md#s-12-2).

<a id="r-12-4"></a>

**Satz 12.4 (äquivariante Abbildungen sind Blockmatrizen).** Seien $\rho_1,\dots,\rho_C$ die paarweise inäquivalenten
irreduziblen Darstellungen, die in $\pi_1\oplus\pi_2$ auftreten, $k_c=\dim\rho_c$, $a_c=m_{\rho_c}(\pi_1)$,
$b_c=m_{\rho_c}(\pi_2)$. Es gibt unitäre Abbildungen
$V_1:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{a_c}\to\mathbb C^{n_1}$ und
$V_2:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{b_c}\to\mathbb C^{n_2}$ mit
$V_i^*\pi_i(g)V_i=\bigoplus_c\rho_c(g)\otimes1$, sodass
$$
\mathrm{Hom}_G(\pi_1,\pi_2)=\Big\{V_2\Big(\bigoplus_{c=1}^C1_{k_c}\otimes y_c\Big)V_1^*:\ y_c\in M_{b_c\times a_c}(\mathbb C)\Big\},
$$
und $(y_c)_c\mapsto T$ ist injektiv. Folglich ist $\mathrm{Hom}_G(\pi_1,\pi_2)\cong\bigoplus_cM_{b_c\times a_c}(\mathbb C)$,
$$
\dim\mathrm{Hom}_G(\pi_1,\pi_2)=\sum_\rho m_\rho(\pi_1)\,m_\rho(\pi_2),
$$
wobei die Summe über alle irreduziblen Darstellungen bis auf Äquivalenz läuft; für $\pi_1=\pi_2$ ist das
$\pi(G)'\cong\bigoplus_cM_{m_c}$ (Proposition 2.7), und $\mathrm{Hom}_G(\pi_1,\pi_2)=0$ gilt genau dann, wenn
keine irreduzible Darstellung in beiden auftritt.

*Beweis.* Man wende Satz 2.6 und Proposition 2.7 auf $\pi=\pi_1\oplus\pi_2$ an: Es gibt eine unitäre Matrix $V$ mit
$V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$ und $V^*\pi(G)'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$. Nach der
Eindeutigkeit der Vielfachheiten (Abschnitt 2.4) ist jede in $\pi$ auftretende irreduzible Darstellung zu genau
einem $\rho_c$ äquivalent, sodass die $\rho_c$ Repräsentanten der Klassen aus der Aussage sind. Beliebige andere
Repräsentanten sind $\rho_c'=w_c\rho_cw_c^*$ mit unitären $w_c$ (Lemma 4.2(c)), und ersetzt man $V$ durch
$V\bigoplus_c(w_c^*\otimes1)$, so geht die erste Normalform in $\bigoplus_c\rho'_c(g)\otimes1_{m_c}$ über und die
zweite bleibt unverändert, da $w_c^*\otimes1$ mit $1\otimes M_{m_c}$ vertauscht; wir dürfen also die $\rho_c$ aus
Proposition 2.7 verwenden. Da
$X\mapsto V^*XV$ ein $*$-Isomorphismus ist und $P_1\in\pi(G)'$ (Lemma 12.3), gilt $V^*P_1V=\bigoplus_c1\otimes p_c$ mit
Projektionen $p_c$. Man wähle unitäre $u_c$ mit $u_c^*p_cu_c=1_{a_c}\oplus0_{b_c}$ (Satz 1.3), wobei nun
$a_c=\mathrm{rank}\,p_c$, $b_c=m_c-a_c$, und ersetze $V$ durch $V\bigoplus_c(1\otimes u_c)$; das ändert keine der
beiden Normalformen, da $1\otimes u_c$ mit $\rho_c(g)\otimes1$ vertauscht. Seien
$\iota_c:\mathbb C^{a_c}\to\mathbb C^{m_c}$ und $\kappa_c:\mathbb C^{b_c}\to\mathbb C^{m_c}$ die Inklusionen der ersten
$a_c$ bzw. der letzten $b_c$ Koordinaten, und
$$
V_1=Q_1^*V\bigoplus_c(1\otimes\iota_c),\qquad V_2=Q_2^*V\bigoplus_c(1\otimes\kappa_c).
$$
$W_1=V\bigoplus_c(1\otimes\iota_c)$ ist eine Isometrie mit Bild $\mathrm{range}\,P_1=\mathrm{range}\,Q_1$, also
$V_1^*V_1=W_1^*P_1W_1=1$, und $V_1$ ist unitär, weil $\sum_ck_ca_c=\mathrm{rank}\,P_1=n_1$. Außerdem gilt
$\pi_1(g)V_1=Q_1^*\pi(g)W_1=Q_1^*V\big(\bigoplus_c\rho_c(g)\otimes1_{m_c}\big)\bigoplus_c(1\otimes\iota_c)
=V_1\bigoplus_c\rho_c(g)\otimes1_{a_c}$. Also ist $\pi_1\cong\bigoplus_c\rho_c\otimes1_{a_c}$ mit inäquivalenten
irreduziblen $\rho_c$, und $a_c=m_{\rho_c}(\pi_1)$ wegen der Eindeutigkeit der Vielfachheiten (Abschnitt 2.4); ebenso
für $V_2$ und $b_c$. Schließlich ist für $X=V\big(\bigoplus_c1\otimes x_c\big)V^*\in\pi(G)'$
$P_2XP_1=V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*$, wobei $y_c=\kappa_c^*x_c\iota_c$ die linke untere
$b_c\times a_c$-Ecke von $x_c$ ist, eine beliebige Matrix. Nach Lemma 12.3 besteht
$\mathrm{Hom}_G(\pi_1,\pi_2)=Q_2^*P_2\,\pi(G)'P_1Q_1$ aus den Matrizen
$Q_2^*V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*Q_1=V_2\big(\bigoplus_c1\otimes y_c\big)V_1^*$, und
$(y_c)\mapsto T$ ist injektiv, da $V_1,V_2$ unitär sind. Eine irreduzible $\rho$, die in $\pi_1\oplus\pi_2$ nicht
auftritt, hat $m_\rho(\pi_1)m_\rho(\pi_2)=0$. $\square$

*In der Bibliothek:* Der Beweis von Satz 12.4 als Algorithmus, die Schichten `EquivariantLinear` und `EquivariantMap`, aufgebaut auf `dc.from_group`, `D.commutant_algebra()` und `D.from_commutant`; von Knoten- zu Kantenmerkmalen von $S_4$ ist $\dim\mathrm{Hom}_G=1\cdot2+1\cdot3=5=B_3$ — [Tutorial, Abschnitt 12.2](../../tutorial/de/12_equivariant_learning.md#s-12-2).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Satz 12.4 ist der
> linear-algebraische Kern der **Auswahlregeln** und der skalare Fall des Wigner–Eckart-Theorems: Ein Operator, der mit der Symmetrie vertauscht, kann nur Zustände desselben
> Symmetrietyps verbinden, und innerhalb eines Typs wirkt er als $1_{k_c}\otimes y_c$, d. h. durch ein „reduziertes
> Matrixelement“ $y_c$, das blind für die Position innerhalb des Multipletts ist. (Operatoren, die sich selbst nach
> einer nichttrivialen Darstellung transformieren, etwa Vektoroperatoren, benötigen zusätzlich die
> Clebsch–Gordan-Koeffizienten.) Ein symmetrischer Hamiltonoperator ist der Fall $\pi_1=\pi_2$; seine Entartungen
> $k_c$ sind die Dimensionen der irreduziblen Darstellungen.

<a id="s-12-2"></a>

## 12.2 Eine äquivariante Schicht mit der Bibliothek

*In der Bibliothek:* Äquivariante Schichten, parametrisiert durch die Blöcke aus Satz 12.4 — [Tutorial, Abschnitt 12.2](../../tutorial/de/12_equivariant_learning.md#s-12-2).

Die zweite Art, eine Schicht äquivariant zu machen, ist die Mittelung über die Gruppe. Der **Twirl** ist
$$
\mathcal P(W)=\frac1{|G|}\sum_{g\in G}\pi_2(g)\,W\,\pi_1(g)^*,\qquad W\in M_{n_2\times n_1}(\mathbb C).
$$

<a id="r-12-5"></a>

**Satz 12.5 (der Twirl ist die Orthogonalprojektion).** $\mathcal P$ ist die bezüglich des Hilbert–Schmidt-Skalarprodukts
orthogonale Projektion von $M_{n_2\times n_1}(\mathbb C)$ auf $\mathrm{Hom}_G(\pi_1,\pi_2)$, und
$$
\dim\mathrm{Hom}_G(\pi_1,\pi_2)=\frac1{|G|}\sum_{g\in G}\overline{\chi_{\pi_1}(g)}\,\chi_{\pi_2}(g).
$$

*Beweis.* Für $h\in G$ ist $\pi_2(h)\mathcal P(W)\pi_1(h)^*=\frac1{|G|}\sum_g\pi_2(hg)W\pi_1(hg)^*=\mathcal P(W)$, da
$g\mapsto hg$ eine Bijektion von $G$ ist; also ist $\mathcal P(W)$ äquivariant. Ist $W$ äquivariant, so ist jeder
Summand gleich $W$, also $\mathcal P(W)=W$. Daher gilt $\mathcal P^2=\mathcal P$ mit Bild
$\mathrm{Hom}_G(\pi_1,\pi_2)$. Die Abbildung ist selbstadjungiert:
$\langle S,\mathcal P(W)\rangle=\frac1{|G|}\sum_g\mathrm{Tr}\big((\pi_2(g)^*S\pi_1(g))^*W\big)$
$=\langle\mathcal P(S),W\rangle$, wobei wir $g\mapsto g^{-1}$ substituiert und $\pi_i(g^{-1})=\pi_i(g)^*$ benutzt haben.
Ein selbstadjungiertes Idempotent ist die Orthogonalprojektion auf sein Bild. Seine Spur ist sein Rang. Die lineare
Abbildung $W\mapsto AWB$ auf $M_{n_2\times n_1}$ hat die Spur $\mathrm{Tr}\,A\mathrm{Tr}\,B$ (der Koeffizient
der Matrixeinheit $E_{ij}$ in $AE_{ij}B$ ist $A_{ii}B_{jj}$); mit $A=\pi_2(g)$, $B=\pi_1(g)^*$ ergibt das die
Formel. $\square$

Zusammen mit Satz 12.4 folgt
$\sum_\rho m_\rho(\pi_1)m_\rho(\pi_2)=|G|^{-1}\sum_g\overline{\chi_{\pi_1}(g)}\chi_{\pi_2}(g)$; für $\pi_1=\pi_2$ ist das
die Charakterformel für $\dim\pi(G)'$ aus Abschnitt 4.2. Die Parametrisierung
$W\mapsto\mathcal P(W)$ erreicht jede äquivariante Abbildung, $\mathcal P(W)$ ist die zu $W$ nächstgelegene
äquivariante Abbildung, und da $\mathcal P$ selbstadjungiert ist, ist der Gradient nach $W$ das Bild unter
$\mathcal P$ des Gradienten nach $T=\mathcal P(W)$. Aber $n_1n_2-\dim\mathrm{Hom}_G$ ihrer Parameter werden an den
Kern verschwendet.

*Rechnung:* Satz 12.5 auf den Kantenmerkmalen und die Dimensionsformel von Knoten- zu Kantenmerkmalen — [Tutorial, Abschnitt 12.2](../../tutorial/de/12_equivariant_learning.md#s-12-2). *In der Bibliothek:* Der Twirl `dc.twirl`, dieselbe Projektion aus der Zerlegung berechnet, `D.conditional_expectation(W, 'commutant')`, und die Charakterformel `dc.commutant_dimension` — [Tutorial, Abschnitt 12.2](../../tutorial/de/12_equivariant_learning.md#s-12-2).

Die bisherigen Schichten sind linear. Ein Netz wechselt sie mit Nichtlinearitäten ab, und auch diese müssen
äquivariant sein. Für eine Permutationsdarstellung ist eine koordinatenweise angewandte Funktion äquivariant, da
$(\pi(g)x)_{gy}=x_y$; für andere Darstellungen ist sie es nicht, und man verwendet Nichtlinearitäten, die auf Normen
wirken oder auf den Koeffizienten bezüglich einer Basis, die $G$ permutiert (Aufgabe 3).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Der Twirl ist die bedingte Erwartung
> auf die Kommutante (Kapitel 2, Aufgabe 3): Mittelt man eine Kopplung über alle Symmetrieoperationen, so bleibt genau
> ihr symmetrischer Anteil übrig. Einen Hamiltonoperator oder eine Wechselwirkung durch Summation über die
> Symmetriegruppe zu „symmetrisieren“ ist dieselbe Projektion.

<a id="s-12-3"></a>

## 12.3 Eingeschränkte und spektrale Schichten

Die Gewichte $w=(w_c)_c$ einer äquivarianten Schicht $x\mapsto Tx$, $T\in\pi(G)'$, leben in der $*$-Algebra
$\bigoplus_cM_{m_c}$ (Satz 12.4), also lässt sich jede Konstruktion der Kapitel 1 und 11 auf sie anwenden. Man schreibe
$\Phi(w)=V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*$ für $w\in\bigoplus_cM_{m_c}$.

<a id="r-12-6"></a>

**Proposition 12.6.** (a) $\Phi$ ist ein $*$-Isomorphismus von $\bigoplus_cM_{m_c}$ auf $\pi(G)'$, und $\Phi(w)$ ist
genau dann selbstadjungiert, unitär bzw. positiv, wenn jedes $w_c$ es ist. (b) Für normales $w$ und eine Funktion $f$
auf $\bigcup_c\mathrm{spec}\,w_c$ gilt $f(\Phi(w))=\Phi\big((f(w_c))_c\big)$; insbesondere ist $f(X)\in\pi(G)'$ für
jedes normale $X\in\pi(G)'$. (c) Die äquivarianten unitären Matrizen sind die $\Phi(u_0e^{s-s^*})$, wobei $s$ ganz
$\bigoplus_cM_{m_c}$ durchläuft, für jedes feste unitäre $u_0\in\bigoplus_cM_{m_c}$; insbesondere sind es die
Exponentiale $e^S$ der schiefadjungierten $S\in\pi(G)'$.

*Beweis.* (a) $\Phi$ ist linear, unital, multiplikativ und mit der Involution verträglich, weil
$(1\otimes x)(1\otimes y)=1\otimes xy$ und $(1\otimes x)^*=1\otimes x^*$; sie ist injektiv und nach Satz 2.6 surjektiv
auf $\pi(G)'$. Daher gilt $\Phi(w)^*=\Phi(w)$ genau dann, wenn $w^*=w$, und $\Phi(w)^*\Phi(w)=1=\Phi(w)\Phi(w)^*$ genau
dann, wenn $w^*w=1=ww^*$. Sind alle $w_c=t_c^*t_c$, so ist $\Phi(w)=\Phi(t)^*\Phi(t)\ge0$ (Korollar 1.7). Ist
umgekehrt $\Phi(w)\ge0$, so gilt für $v\in\mathbb C^{m_c}$ und den Vektor $\xi=V(e_1\otimes v)$ im Sektor $c$
$0\le\langle\xi,\Phi(w)\xi\rangle=\langle v,w_cv\rangle$.
(b) Sei $w_c=u_c\Lambda_cu_c^*$ (Satz 1.3). Dann ist $\Phi(w)=U\big(\bigoplus_c1\otimes\Lambda_c\big)U^*$ mit der
unitären Matrix $U=V\bigoplus_c(1\otimes u_c)$, eine Diagonalisierung, in der die Eigenwerte von $w_c$ je $k_c$-mal
wiederholt auftreten; nach Definition und Satz 1.5(a) gilt
$f(\Phi(w))=U\big(\bigoplus_c1\otimes f(\Lambda_c)\big)U^*=\Phi\big((f(w_c))_c\big)$.
(c) Nach (a) sind die äquivarianten unitären Matrizen die $\Phi(u)$ mit unitärem $u$ in $\bigoplus_cM_{m_c}$, und
diese sind nach Lemma 11.1 die $u_0e^{s-s^*}$. Mit $u_0=1$ liefert (b) $\Phi(e^{s-s^*})=e^{\Phi(s-s^*)}$, und
$\Phi(s-s^*)$ durchläuft die schiefadjungierten Elemente von $\pi(G)'$. $\square$

*In der Bibliothek:* Unitäre und positive äquivariante Schichten mit `optimize.UnitaryParam` und `PositiveParam` auf `D.commutant_algebra()`, spektrale Schichten mit `apply_function`; das Training einer äquivarianten unitären Schicht und (b) für $e^{-tL}$ mit seiner Ableitung nach $t$ — [Tutorial, Abschnitt 12.3](../../tutorial/de/12_equivariant_learning.md#s-12-3).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Ein symmetrischer Hamiltonoperator
> $H\in\pi(G)'$ erzeugt eine symmetrische Zeitentwicklung $e^{-itH}\in\pi(G)'$ (Kapitel 8), die nie Sektoren mischt
> und nur auf den Vielfachheitsräumen wirkt. Eine äquivariante unitäre Schicht ist eine solche Entwicklung, gelernt
> statt vorgegeben.

<a id="s-12-4"></a>

## 12.4 Invarianten

Ein Vektor $x$ heißt **invariant**, wenn $\pi(g)x=x$ für alle $g$ gilt; eine Funktion $q$ auf $\mathbb C^n$ heißt
invariant, wenn $q(\pi(g)x)=q(x)$. Invariante Merkmale sind das, was ein invariantes Netz (etwa ein
Graphenklassifikator) ausgeben kann.

<a id="r-12-7"></a>

**Proposition 12.7.** (a) Die invarianten Vektoren bilden das Bild der Orthogonalprojektion
$P_0=|G|^{-1}\sum_g\pi(g)$; ihre Dimension ist $m_{\mathrm{triv}}(\pi)=|G|^{-1}\sum_g\chi_\pi(g)$, und wenn die triviale
Darstellung auftritt, als $\rho_{c_0}$, dann ist $P_0$ die minimale zentrale Projektion $P_{c_0}$ von $\pi(G)''$.
(b) Eine hermitesche Form $q_A(x)=\langle x,Ax\rangle$, $A=A^*$, ist genau dann invariant, wenn $A\in\pi(G)'$; die
invarianten hermiteschen Formen bilden einen reellen Vektorraum der Dimension $\sum_cm_c^2$.
(c) Die homogenen invarianten Polynome vom Grad $2$ in den Koordinaten $x_1,\dots,x_n$ bilden einen Raum der Dimension
$$
\frac1{|G|}\sum_{g\in G}\frac{\chi_\pi(g)^2+\chi_\pi(g^2)}2 .
$$
Für eine Permutationsdarstellung ist das die Anzahl der Bahnen von $G$ auf den Multimengen $\{x,y\}$ aus zwei
Elementen von $\Omega$.

*Beweis.* (a) Die invarianten Vektoren sind $\mathrm{Hom}_G(\mathrm{triv},\pi)$, mit der trivialen Darstellung
$\mathrm{triv}$ auf $\mathbb C=M_{1\times1}$; für sie ist der Twirl $\mathcal P(x)=P_0x$, also liefert Satz 12.5 die
Projektion und die Dimension, und Satz 12.4 liefert $m_{\mathrm{triv}}(\pi)$. In der Basis aus Satz 2.6 schreibe man
einen invarianten Vektor als $\xi=\sum_c\xi_c$ mit
$\xi_c=\sum_jv_{c,j}\otimes e_j\in\mathbb C^{k_c}\otimes\mathbb C^{m_c}$. Da $\pi(g)$ die Sektoren erhält, gilt
$\rho_c(g)v_{c,j}=v_{c,j}$ für alle $g$. Ein von null verschiedenes $v_{c,j}$ spannt eine invariante Gerade auf, die
die triviale Darstellung trägt; dann ist $\rho_c$ als irreduzible Darstellung trivial und $c=c_0$. Daher sind die
invarianten Vektoren die Vektoren des Sektors $c_0$, auf dem $\pi(g)$ als $1$ wirkt: das Bild von $P_{c_0}$, das
folglich $P_0$ ist.
(b) $q_A(\pi(g)x)=\langle x,\pi(g)^*A\pi(g)x\rangle$, und eine Matrix $D$ mit $\langle x,Dx\rangle=0$ für alle $x$
verschwindet (Polarisation, Beweis von Korollar 1.4). Also ist $q_A$ genau dann invariant, wenn $\pi(g)^*A\pi(g)=A$
für alle $g$, d. h. $A\in\pi(G)'$. Nach Proposition 12.6(a) sind die selbstadjungierten Elemente von $\pi(G)'$ die
$\Phi(w)$ mit hermiteschen $w_c$, ein reeller Raum der Dimension $\sum_cm_c^2$.
(c) Ein quadratisches Polynom $p(x)=x^TBx$ ist durch seine symmetrische Matrix $B$ bestimmt. Man setze
$(g\cdot p)(x)=p(\pi(g)^{-1}x)$; dann ist $g\cdot(h\cdot p)=(gh)\cdot p$, und $g\cdot p$ hat die Matrix
$\sigma(g)(B)=A^TBA$, $A=\pi(g)^{-1}$. Das Mittel $M=|G|^{-1}\sum_g\sigma(g)$ erfüllt $\sigma(h)M=M$ und $MB=B$ für
invariantes $B$, ist also (wie im Beweis von Satz 12.5) idempotent mit den invarianten Polynomen als Bild, und deren
Dimension ist $\mathrm{Tr}\,M=|G|^{-1}\sum_g\mathrm{Tr}\,\sigma(g)$ (ein Idempotent ist $1\oplus0$ in einer an
Bild $\oplus$ Kern angepassten Basis, also ist seine Spur sein Rang). $A$ ist unitär, also $A=R\Lambda R^{-1}$
mit $\Lambda=\mathrm{diag}(\lambda)$ (Satz 1.3). Die symmetrischen Matrizen
$B_{ij}=(R^{-1})^T(e_ie_j^T+e_je_i^T)R^{-1}$, $i\le j$, bilden eine Basis, und
$A^TB_{ij}A=(R^{-1})^T\Lambda(e_ie_j^T+e_je_i^T)\Lambda R^{-1}=\lambda_i\lambda_jB_{ij}$. Also ist
$\mathrm{Tr}\,\sigma(g)=\sum_{i\le j}\lambda_i\lambda_j=\frac12\big((\mathrm{Tr}\,A)^2+\mathrm{Tr}\,A^2\big)
=\frac12\big(\chi_\pi(g^{-1})^2+\chi_\pi(g^{-2})\big)$, und Summation über $g^{-1}$ statt $g$ ergibt die Formel. Für
eine Permutationsdarstellung permutiert $G$ die Monome $x_ax_b$, also ist ein Polynom genau dann invariant, wenn seine
Koeffizienten auf den Bahnen der Monome konstant sind. $\square$

*Rechnung:* (a), (b) und die beiden Abzählungen aus (c) für die Kantenmerkmale von $S_4$ — [Tutorial, Abschnitt 12.4](../../tutorial/de/12_equivariant_learning.md#s-12-4).

Für Graphen werden die invarianten Kantenmerkmale von den Indikatoren der Diagonale und ihres Komplements aufgespannt
($B_2=2$), und die invarianten linearen Ausleseabbildungen (*readouts*) sind die Summe der Diagonaleinträge und die
der Nichtdiagonaleinträge. Auf den Knotenmerkmalen $\mathbb C^4$ zählt (c) die Bahnen von $S_4$ auf den Multimengen
aus zwei Knoten: Es gibt zwei.

Auf $\mathbb C^4$ sind es $\sum_ix_i^2$ und $\sum_{i<j}x_ix_j$. Der Beweis von (c) funktioniert in jedem Grad:

<a id="r-12-8"></a>

**Proposition 12.8 (Molien-Formel).** Sei $\mathbb C[x]_d$ der Raum der homogenen Polynome vom Grad $d$ in
$x_1,\dots,x_n$ und $\mathbb C[x]^G_d$ der Unterraum der invarianten, $p(\pi(g)x)=p(x)$ für alle $g$. Für $|t|<1$ gilt
$$
\sum_{d\ge0}\dim\mathbb C[x]^G_d\,t^d=\frac1{|G|}\sum_{g\in G}\frac1{\det\big(1-t\pi(g)\big)} .
$$
Für eine Permutationsdarstellung ist $\dim\mathbb C[x]^G_d$ die Anzahl der Bahnen von $G$ auf den Monomen vom Grad $d$
(den Multimengen aus $d$ Elementen von $\Omega$). Der Koeffizient von $t^2$ ist Proposition 12.7(c).

*Beweis.* Für eine invertierbare Matrix $M$ sei $\phi_M(p)(x)=p(Mx)$ (Substitution). $\phi_M$ ist linear, bildet
$\mathbb C[x]_d$ in sich ab, und $\phi_N\circ\phi_M=\phi_{MN}$, also ist $\phi_M$ bijektiv mit Inverser
$\phi_{M^{-1}}$. Daher erfüllt $\sigma(g)=\phi_{\pi(g)^{-1}}$ die Beziehung $\sigma(g)\sigma(h)=\sigma(gh)$, und $p$ ist
genau dann invariant, wenn $\sigma(g)p=p$ für alle $g$. Wie im Beweis von Proposition 12.7(c) ist
$M_d=|G|^{-1}\sum_g\sigma(g)$ auf $\mathbb C[x]_d$ idempotent mit Bild $\mathbb C[x]^G_d$, also
$\dim\mathbb C[x]^G_d=\mathrm{Tr}\,M_d=|G|^{-1}\sum_g\mathrm{Tr}\,\sigma_d(g)$, mit $\sigma_d(g)$ der Einschränkung
von $\sigma(g)$ auf $\mathbb C[x]_d$. Man fixiere $g$ und schreibe die unitäre Matrix $A=\pi(g)^{-1}$ als
$A=R\Lambda R^{-1}$, $\Lambda=\mathrm{diag}(\lambda_1,\dots,\lambda_n)$, $|\lambda_i|=1$ (Satz 1.3). Mit den
Linearformen $\ell_i(x)=(R^{-1}x)_i$ bilden die Produkte $\ell^\alpha=\prod_i\ell_i^{\alpha_i}=\phi_{R^{-1}}(x^\alpha)$,
$|\alpha|=d$, eine Basis von $\mathbb C[x]_d$ (das Bild der Monombasis unter der Bijektion $\phi_{R^{-1}}$), und wegen
$R^{-1}A=\Lambda R^{-1}$ gilt
$$
(\sigma(g)\ell^\alpha)(x)=\prod_i\big((R^{-1}Ax)_i\big)^{\alpha_i}=\prod_i\big(\lambda_i\ell_i(x)\big)^{\alpha_i}=\lambda^\alpha\ell^\alpha(x).
$$
Also ist $\mathrm{Tr}\,\sigma_d(g)=\sum_{|\alpha|=d}\lambda^\alpha$. Für $|t|<1$ darf das Produkt der $n$ absolut
konvergenten Reihen $(1-t\lambda_i)^{-1}=\sum_{a\ge0}t^a\lambda_i^a$ ausmultipliziert und nach Graden umgeordnet
werden:
$$
\frac1{\det(1-tA)}=\prod_{i=1}^n\frac1{1-t\lambda_i}=\sum_\alpha t^{|\alpha|}\lambda^\alpha=\sum_{d\ge0}t^d\mathrm{Tr}\,\sigma_d(g).
$$
Mittelung über $g$ und Ersetzen von $g$ durch $g^{-1}$ (sodass $A=\pi(g)$) liefert die Formel; beide Seiten
konvergieren, da $\dim\mathbb C[x]^G_d\le\dim\mathbb C[x]_d=\binom{n+d-1}d$. Für eine Permutationsdarstellung ist
$(\pi(g)^{-1}x)_y=x_{gy}$, also bildet $\sigma(g)$ das Monom $\prod_jx_{y_j}$ auf das Monom $\prod_jx_{gy_j}$ ab: $G$
permutiert die Monome vom Grad $d$, ein Polynom ist genau dann invariant, wenn seine Koeffizienten auf den Bahnen
konstant sind, und die Bahnsummen bilden eine Basis von $\mathbb C[x]^G_d$. $\square$

(Das ist der Satz von Molien aus dem Jahr 1897; zur Vertiefung siehe B. Sturmfels, *Algorithms in Invariant
Theory*, Springer 2008, Abschnitt 2.2.) Für $S_4$ auf $\mathbb C^4$ sind die Anzahlen $1,1,2,3,5$ für $d\le4$ die
Partitionen von $d$ in höchstens vier Teile (die Bahnen von Multimengen sind durch ihre sortierten Vielfachheiten
festgelegt).

*Rechnung:* Die Molien-Reihe bis $t^4$ im Vergleich mit den Bahnen der Monome, für Knoten- und Kantenmerkmale — [Tutorial, Abschnitt 12.4](../../tutorial/de/12_equivariant_learning.md#s-12-4).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). In der Landau-Theorie der
> Phasenübergänge ist die freie Energie ein invariantes Polynom in einem reellen Ordnungsparameter $x$, auf den $G$
> durch reelle orthogonale Matrizen wirkt. Ihr quadratischer Anteil ist $x^TBx$ mit $B$ in der Kommutante (da
> $\pi(g)^T=\pi(g)^{-1}$), also nach Proposition 12.6(a) blockdiagonal bezüglich der Sektoren. Die Landau-Theorie
> postuliert dann, dass die symmetrische Phase in dem Sektor instabil wird, in dem zuerst ein Eigenwert von $B$ das
> Vorzeichen wechselt, und dass generisch die irreduzible Darstellung dieses Sektors die Symmetrie der geordneten
> Phase bestimmt.

<a id="s-12-5"></a>

## 12.5 Womit man vergleichen sollte

*In der Bibliothek:* Vergleich mit anderer Software (escnn, e3nn, EMLP) und was die Bibliothek hinzufügt — [Tutorial, Abschnitt 12.5](../../tutorial/de/12_equivariant_learning.md#s-12-5).

<a id="x-12"></a>

## Aufgaben

1. (**DeepSets mit Kanälen.**) $S_n$ ($n\ge2$) wirkt auf $\mathbb C^n\otimes\mathbb C^d$ durch
   $P_\sigma\otimes1_d$. Man zeige, dass die äquivarianten Abbildungen von $\mathbb C^n\otimes\mathbb C^{d_1}$ nach
   $\mathbb C^n\otimes\mathbb C^{d_2}$ genau die $1\otimes A+J\otimes B$ mit $A,B\in M_{d_2\times d_1}$ sind.
2. (**Reelle Schichten.**) Die Darstellungen $\pi_1,\pi_2$ mögen reelle Matrizen haben (z. B.
   Permutationsdarstellungen). Man zeige, dass die reellen äquivarianten Matrizen einen reellen Vektorraum der
   Dimension $\dim_{\mathbb C}\mathrm{Hom}_G(\pi_1,\pi_2)$ bilden und dass $T\mapsto\mathrm{Re}\,T$ den Raum
   $\mathrm{Hom}_G(\pi_1,\pi_2)$ auf ihn abbildet.
3. (**Nichtlinearitäten.**) (a) Man zeige, dass die koordinatenweise Anwendung einer Funktion
   $f:\mathbb C\to\mathbb C$ mit jeder Permutationsdarstellung vertauscht. (b) Man zeige, dass sie für
   $f(t)=\max(t,0)$ nicht mit der Drehung von $\mathbb R^2$ um $\pi/2$ vertauscht. (c) Man zeige, dass
   $x\mapsto\varphi(\|x\|)x$ mit jeder unitären Darstellung vertauscht.
4. (**Wie viele Stichproben legen eine Schicht fest?**) Seien $T=\Phi(w)\in\pi(G)'$ und
   $x_1,\dots,x_N\in\mathbb C^n$. Man schreibe $V^*x_i$ im Sektor $c$ als $\sum_{a,j}(X_{i,c})_{aj}e_a\otimes e_j$ mit
   $X_{i,c}\in M_{k_c\times m_c}$, und sei $S_c\in M_{Nk_c\times m_c}$ die Matrix, die durch Übereinanderstapeln von
   $X_{1,c},\dots,X_{N,c}$ entsteht. Man zeige, dass die äquivarianten $T$ mit $Tx_1=\dots=Tx_N=0$ einen Raum der
   Dimension $\sum_cm_c(m_c-\mathrm{rank}\,S_c)$ bilden; für unabhängige standardnormalverteilte komplexe
   Gauß-Vektoren $x_i$ ist dies fast sicher $\sum_cm_c\max(0,m_c-Nk_c)$. Man folgere, dass ein äquivariantes $T$
   genau dann (unter den äquivarianten Abbildungen) durch $Tx_1,\dots,Tx_N$ festgelegt ist, wenn diese Dimension $0$
   ist, und dass für die Kantendarstellung zwei Gaußsche Stichproben eine äquivariante Schicht fast sicher festlegen
   und eine einzige Stichprobe dies nie tut, während eine beliebige $16\times16$-Matrix nie durch weniger als
   sechzehn Stichproben festgelegt ist.

*Rechnung:* Aufgaben 1, 3 und 4 numerisch — [Tutorial, Aufgaben](../../tutorial/de/12_equivariant_learning.md#x-12).

<details>
<summary>Lösungen</summary>

1. $P_\sigma J=J=JP_\sigma$, also sind diese Abbildungen äquivariant. Sie sind linear unabhängig: Aus
   $1\otimes A+J\otimes B=0$ folgt mit den Nichtdiagonalblöcken $B=0$ und dann $A=0$. Nach Korollar 12.2(a) und
   Satz 2.6 ist $P_\sigma\otimes1_d\cong(\mathrm{triv}\otimes1_d)\oplus(\mathrm{std}\otimes1_d)$, also haben beide
   irreduziblen Darstellungen die Vielfachheit $d$, und Satz 12.4 liefert $\dim\mathrm{Hom}_G=2d_1d_2$, die
   Dimension des Raums der $1\otimes A+J\otimes B$.

2. Ist $T$ äquivariant, so auch $\bar T$, da $\overline{\pi_i(g)}=\pi_i(g)$. Daher sind
   $\mathrm{Re}\,T=(T+\bar T)/2$ und $\mathrm{Im}\,T=(T-\bar T)/2i$ reell und äquivariant, also ist die
   Abbildung $T\mapsto\mathrm{Re}\,T$ surjektiv auf die reellen äquivarianten Matrizen (auf ihnen ist sie die
   Identität), und $\mathrm{Hom}_G(\pi_1,\pi_2)=H_{\mathbb R}\oplus iH_{\mathbb R}$, wobei $H_{\mathbb R}$ die
   reellen äquivarianten Matrizen bezeichnet. Eine reelle Basis von $H_{\mathbb R}$ ist dann eine komplexe Basis von
   $\mathrm{Hom}_G$: Sie spannt den Raum auf, und eine komplexe Relation $\sum_j(\alpha_j+i\beta_j)T_j=0$ mit
   reellen $T_j$ zerfällt in Real- und Imaginärteil. In Lemma 12.1 bilden die Bahnindikatoren eine solche reelle
   Basis. (Die Blöcke $y_c$ aus Satz 12.4 sind selbst für reelle Daten komplex, da das von der Bibliothek
   berechnete $V$ nicht reell sein muss; $\mathrm{Re}\,T$ durch komplexe $y_c$ zu parametrisieren verwendet
   doppelt so viele reelle Parameter wie nötig.)

3. (a) $(\pi(g)f(x))_{gy}=f(x)_y=f(x_y)=f((\pi(g)x)_{gy})$. (b) Für $R(x_1,x_2)=(-x_2,x_1)$ und $x=(1,-1)$ ist
   $f(Rx)=f(1,1)=(1,1)$, aber $Rf(x)=R(1,0)=(0,1)$. (c) $\|\pi(g)x\|=\|x\|$, also
   $\varphi(\|\pi(g)x\|)\pi(g)x=\pi(g)\varphi(\|x\|)x$.

4. $(1\otimes w_c)(e_a\otimes e_j)=\sum_le_a\otimes(w_c)_{lj}e_l$, also hat $V^*Tx_i$ im Sektor $c$ die
   Koeffizientenmatrix $X_{i,c}w_c^T$. Daher gilt $Tx_i=0$ für alle $i$ genau dann, wenn $S_cw_c^T=0$ für alle $c$,
   d. h. genau dann, wenn jede Spalte von $w_c^T$ in $\ker S_c$ liegt, einem Raum der Dimension
   $m_c-\mathrm{rank}\,S_c$: Die $w_c$ mit dieser Eigenschaft bilden einen Raum der Dimension
   $m_c(m_c-\mathrm{rank}\,S_c)$, und $\Phi$ ist injektiv. Für Gaußsche $x_i$ sind die Koordinaten von $V^*x_i$
   diejenigen von $x_i$ in der Orthonormalbasis der Spalten von $V$, wieder unabhängige standardnormalverteilte
   komplexe Zufallsvariablen (Proposition A.16(a)), und unabhängig in $i$; also ist $S_c$ eine
   $Nk_c\times m_c$-Ginibre-Matrix, die eine Dichte hat (Proposition A.16(b)). Ihr Rang ist höchstens
   $\min(Nk_c,m_c)$, mit Gleichheit, außer wenn die linke obere quadratische Teilmatrix dieser Größe singulär ist;
   $|\det|^2$ dieser Teilmatrix ist ein reelles Polynom in den Real- und Imaginärteilen der Einträge, das nicht
   identisch verschwindet (es ist $1$, wenn die Teilmatrix die Einheitsmatrix ist), also ist der Rang nach
   Korollar A.7(b) fast sicher $\min(Nk_c,m_c)$. Für die Kantendarstellung, $(k,m)=(1,2),(2,1),(3,1),(3,3)$ (die
   berechnete Zerlegung), trägt nur der Sektor $(1,2)$ bei: $2\cdot(2-1)=2$ für $N=1$ und $0$ für $N=2$. Eine
   beliebige $16\times16$-Matrix ist der Fall der trivialen Gruppe (ein Sektor, $k=1$, $m=16$): Die Dimension ist
   $16(16-\mathrm{rank}\,S)\ge16(16-N)>0$ für jede Wahl von $N<16$ Stichproben, da $\mathrm{rank}\,S\le N$,
   und fast sicher $0$ für $N=16$ Gaußsche Stichproben. Da $T\mapsto(Tx_1,\dots,Tx_N)$ linear ist, stimmen zwei
   äquivariante Abbildungen genau dann auf den Stichproben überein, wenn ihre Differenz im eben berechneten Raum
   liegt, was das Kriterium liefert. Für $N=1$ hat der Raum für *jedes* $x_1$ mindestens die Dimension $2$, da
   $\mathrm{rank}\,S_c\le Nk_c=1<m_c=2$ im Sektor $(1,2)$.

</details>

**Wie es weitergeht.** Das letzte Kapitel, [Kapitel 13](13_jordan_algebras.md), fragt, warum die Quantenmechanik *-Algebren und nicht Jordan-Algebren verwendet. [Kapitel 4](04_symmetries.md) entwickelt die hier
verwendete Darstellungstheorie: das Lemma von Schur, Charaktere und den randomisierten Algorithmus hinter
`dc.from_group` (Abschnitte 4.1–4.3). [Kapitel 9](09_spin_chains.md) wendet dieselben Sektoren auf Spinketten an, und das
Haar-Maß auf $\prod_cU(m_c)$ aus [Kapitel 5](05_random_matrices.md) liefert die natürliche zufällige Initialisierung
einer äquivarianten unitären Schicht. Zur Mathematik siehe J.-P. Serre, *Linear Representations of Finite Groups*
(Springer 1977) und W. Fulton und J. Harris, *Representation Theory* (Springer 1991); zum äquivarianten Lernen
M. Bronstein, J. Bruna, T. Cohen und P. Veličković, *Geometric Deep Learning* (2021), und M. Weiler et al.,
*Equivariant and Coordinate Independent Convolutional Networks* (2023); zum Zerlegungsalgorithmus die Literaturangaben
in `torch_vn_algebra/decompose.py`.
