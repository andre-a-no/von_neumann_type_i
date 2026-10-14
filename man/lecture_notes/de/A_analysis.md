[English](../A_analysis.md) · [Español](../es/A_analysis.md) · [Français](../fr/A_analysis.md) · **Deutsch** · [中文](../zh/A_analysis.md) · [日本語](../ja/A_analysis.md) · [Русский](../ru/A_analysis.md)

# Anhang A. Tatsachen aus Analysis und Wahrscheinlichkeitstheorie

Die Kapitel verwenden einige Tatsachen aus Analysis und Wahrscheinlichkeitstheorie, die über die Voraussetzungen des
Kurses hinausgehen, aber keine eigenständigen Sätze sind: Nullstellenmengen von Polynomen sind Nullmengen, die
geordneten Eigenwerte einer hermiteschen Matrix hängen stetig von ihr ab, nicht ausgeartete lineare Bilder Gaußscher
Vektoren haben Dichten, die Richtung eines Gauß-Vektors ist gleichverteilt auf der Sphäre und unabhängig von seiner Länge,
Gaußsche Momente sind Summen über Paarungen, und die Abstände gleichverteilter Punkte sind normierte exponentialverteilte
Variablen. Dieser Anhang beweist sie, damit die Kapitel sie mit ihrer Nummer zitieren können.

Die numerischen Prüfungen dieses Anhangs stehen in [Tutorial, Anhang A](../../tutorial/de/A_analysis.md).

**Verwendete Voraussetzungen.** Lineare Algebra; reelle Analysis (Grenzwerte, Stetigkeit, Kompaktheit und der Satz von
Bolzano–Weierstraß im $\mathbb R^n$, Potenzreihen in einer Variablen, das Riemann-Integral einer stetigen Funktion, die
Transformationsformel); Maß- und Integrationstheorie (Lebesgue-Maß, Produktmaße und Fubini–Tonelli, monotone und majorisierte
Konvergenz, Bildmaße, Dichten, Unabhängigkeit); elementare Wahrscheinlichkeitstheorie (die Verteilung $N(0,s^2)$ mit Dichte
$(2\pi s^2)^{-1/2}e^{-y^2/2s^2}$ und momenterzeugender Funktion $\mathbb Ee^{aY}=e^{a^2s^2/2}$). Für eine invertierbare
lineare Abbildung $L$ des $\mathbb R^m$ und eine Borelmenge $E$ gilt $\lambda_m(LE)=|\det L|\,\lambda_m(E)$, wobei
$\lambda_m$ das Lebesgue-Maß ist; das ist der lineare Fall der Transformationsformel.

Durchweg bedeutet „messbar“ Borel-messbar, $|x|$ ist die euklidische Norm von $x\in\mathbb R^n$ oder
$\mathbb C^n$, $\|A\|=\sup_{|x|=1}|Ax|$ die Operatornorm und $\|A\|_2=(\mathrm{Tr}\,A^*A)^{1/2}$ die
Frobeniusnorm.

<a id="s-A-1"></a>

## A.1 Maße und Integrale

<a id="r-A-1"></a>

**Lemma A.1 (Dynkin).** Seien $\mu,\nu$ endliche Maße auf einem Messraum $(\Omega,\mathcal F)$ mit
$\mu(\Omega)=\nu(\Omega)$, und sei $\mathcal P\subseteq\mathcal F$ abgeschlossen unter endlichen Durchschnitten und erzeuge
$\mathcal F$ als $\sigma$-Algebra. Gilt $\mu(P)=\nu(P)$ für alle $P\in\mathcal P$, so ist $\mu=\nu$.

*Beweis.* Man nennt $\mathcal D\subseteq\mathcal F$ ein *Dynkin-System*, wenn $\Omega\in\mathcal D$, $B\setminus A\in\mathcal D$
für alle $A\subseteq B$ aus $\mathcal D$ und $\bigcup_kA_k\in\mathcal D$ für jede aufsteigende Folge
$A_k\in\mathcal D$ gilt. Die Familie $\{A\in\mathcal F:\mu(A)=\nu(A)\}$ ist ein Dynkin-System, da
$\mu(B\setminus A)=\mu(B)-\mu(A)$ für endliche Maße gilt und Maße entlang aufsteigender Vereinigungen stetig sind. Sei $d$
der Durchschnitt aller Dynkin-Systeme, die $\mathcal P$ enthalten; er ist selbst ein Dynkin-System.
Die Familie $\mathcal D_1=\{A\in d:A\cap P\in d\text{ für alle }P\in\mathcal P\}$ ist ein Dynkin-System (für
$A\subseteq A'$ ist $(A'\setminus A)\cap P=(A'\cap P)\setminus(A\cap P)$) und enthält $\mathcal P$, also ist
$\mathcal D_1=d$. Dann ist $\mathcal D_2=\{A\in d:A\cap B\in d\text{ für alle }B\in d\}$ ein Dynkin-System, das
$\mathcal P$ enthält (weil $\mathcal D_1=d$), also $\mathcal D_2=d$: $d$ ist abgeschlossen unter endlichen Durchschnitten.
Folglich ist $d$ eine $\sigma$-Algebra: $\Omega\setminus A\in d$, $A\cup B=\Omega\setminus\big((\Omega\setminus A)\cap(\Omega\setminus B)\big)\in d$,
und eine abzählbare Vereinigung ist die aufsteigende Vereinigung endlicher Vereinigungen. Also ist
$\mathcal F=\sigma(\mathcal P)\subseteq d\subseteq\{\mu=\nu\}$. $\square$

<a id="r-A-2"></a>

**Korollar A.2.** Zwei endliche Borel-Maße $\mu,\nu$ auf einem metrischen Raum $M$ mit $\int f\,d\mu=\int f\,d\nu$
für jede beschränkte stetige Funktion $f:M\to\mathbb R$ sind gleich.

*Beweis.* Mit $f=1$ ist $\mu(M)=\nu(M)$. Für eine nichtleere abgeschlossene Menge $K$ (für $K=\emptyset$ ist nichts zu zeigen) nehmen die stetigen Funktionen
$f_k(x)=\max(0,1-k\mathrm{dist}(x,K))$ Werte in $[0,1]$ an und fallen gegen die Indikatorfunktion von $K$
(für $x\notin K$ ist $\mathrm{dist}(x,K)>0$). Nach dem Satz von der majorisierten Konvergenz ist $\mu(K)=\nu(K)$. Die abgeschlossenen Mengen sind
abgeschlossen unter endlichen Durchschnitten und erzeugen die Borelmengen; man wende Lemma A.1 an. $\square$

Lemma 5.1 verwendet genau dieses Argument.

<a id="r-A-3"></a>

**Lemma A.3 (Lebesgue-Maß auf einem Vektorraum).** Sei $W$ ein reeller Vektorraum der Dimension $m$. Für eine
lineare Bijektion $\varphi:\mathbb R^m\to W$ setze man $\lambda_\varphi(E)=\lambda_m(\varphi^{-1}E)$ für Borelmengen $E\subseteq W$.

(a) Für eine weitere lineare Bijektion $\psi$ ist $\lambda_\psi=|\det(\psi^{-1}\varphi)|\,\lambda_\varphi$. Die Nullmengen von
$W$ und die Frage, ob ein Zufallsvektor in $W$ eine Dichte hat, hängen also nicht von der Wahl von $\varphi$ ab.

(b) Für jede lineare Bijektion $L:W\to W$ ist $\lambda_\varphi(LE)=|\det L|\,\lambda_\varphi(E)$. Ist $L$ orthogonal
bezüglich eines Skalarprodukts auf $W$, so ist $|\det L|=1$, und $L$ erhält $\lambda_\varphi$.

(c) Hat ein Zufallsvektor $X$ in $W$ die Dichte $f$ bezüglich $\lambda_\varphi$, so hat $LX$ die Dichte
$y\mapsto f(L^{-1}y)/|\det L|$.

*Beweis.* (a) $\lambda_\psi(E)=\lambda_m\big((\psi^{-1}\varphi)\varphi^{-1}E\big)=|\det(\psi^{-1}\varphi)|\,\lambda_m(\varphi^{-1}E)$.
(b) $\lambda_\varphi(LE)=\lambda_m\big((\varphi^{-1}L\varphi)\varphi^{-1}E\big)=|\det(\varphi^{-1}L\varphi)|\,\lambda_\varphi(E)$, und
$\det(\varphi^{-1}L\varphi)=\det L$ hängt nicht von der Basis ab. In einer Orthonormalbasis hat ein orthogonales $L$ eine
Matrix $O$ mit $O^TO=1$, also $(\det L)^2=1$. (c) Nach (b) ist das Bildmaß von $\lambda_\varphi$ unter $L$ gleich
$|\det L|^{-1}\lambda_\varphi$, also gilt für Borelmengen $E$
$\mathbb P(LX\in E)=\int1_E(Lx)f(x)\,d\lambda_\varphi(x)=|\det L|^{-1}\int1_E(y)f(L^{-1}y)\,d\lambda_\varphi(y)$. $\square$

*Konventionen.* „Lebesgue-Maß“ auf einem reellen Vektorraum bedeutet ein beliebiges $\lambda_\varphi$; nach (a) spielt die Wahl
für Nullmengen und Dichten nie eine Rolle. Auf $\mathbb C^n$ und auf $M_{a\times b}(\mathbb C)$ verwenden wir Real- und Imaginärteile
der Einträge als Koordinaten; sie sind orthonormal für das reelle Skalarprodukt
$\mathrm{Re}\,\mathrm{Tr}(X^*Y)$, bezüglich dessen jede Abbildung $X\mapsto VXW$ mit unitären $V,W$ orthogonal ist
(sie erhält $\mathrm{Tr}(X^*Y)$). Auf den hermiteschen Matrizen $\mathrm{Herm}_n$ verwenden wir die Koordinaten
$h_{ii}$ und $\mathrm{Re}\,h_{ij},\mathrm{Im}\,h_{ij}$ ($i<j$); sie sind nicht orthonormal für das reelle Skalarprodukt
$\mathrm{Tr}(XY)$, aber nach (b) erhalten die Abbildungen $H\mapsto UHU^*$ ($U$ unitär), die bezüglich dieses Skalarprodukts
orthogonal sind, dennoch das Lebesgue-Maß. Ebenso verwenden wir auf den reellen symmetrischen Matrizen $\mathrm{Sym}_n(\mathbb R)$ die
Koordinaten $h_{ii}$ und $h_{ij}$ ($i<j$), und die Abbildungen $H\mapsto OHO^T$ ($O\in O(n)$), orthogonal für $\mathrm{Tr}(XY)$,
erhalten das Lebesgue-Maß.

<a id="r-A-4"></a>

**Lemma A.4 (Integrale vektorwertiger Funktionen).** Sei $V$ ein endlichdimensionaler reeller oder komplexer Vektorraum
und $F:[a,b]\to V$ stetig. Man definiere $\int_a^bF(t)\,dt$ koordinatenweise in einer Basis (Real- und Imaginärteil
getrennt). Dann gilt:

(a) für Zerlegungen $a=t_0<\dots<t_N=b$ mit Zwischenstellen $\tau_i\in[t_{i-1},t_i]$ konvergieren die Riemann-Summen
$\sum_iF(\tau_i)(t_i-t_{i-1})$ gegen $\int_a^bF$, wenn $\max_i(t_i-t_{i-1})\to0$; insbesondere hängt das Integral
nicht von der Basis ab;

(b) $L\int_a^bF=\int_a^bLF$ für jede ($\mathbb R$-)lineare Abbildung $L:V\to V'$ in einen endlichdimensionalen Raum $V'$;

(c) $\|\int_a^bF\|\le\int_a^b\|F(t)\|\,dt$ für jede Norm auf $V$;

(d) liegen alle $F(t)$ in einer abgeschlossenen konvexen Menge $K$, so ist $\frac1{b-a}\int_a^bF\in K$; ist $K$ zudem ein Kegel,
so ist $\int_a^bF\in K$;

(e) $t\mapsto\int_a^tF$ ist differenzierbar mit Ableitung $F(t)$; $\int_a^bF=\int_a^cF+\int_c^bF$, und
$\int_a^bF(s+c)\,ds=\int_{a+c}^{b+c}F(s)\,ds$.

*Beweis.* Jede Koordinate von $F$ ist eine stetige reelle Funktion, also gelten (a) und (e) koordinatenweise; (b) folgt
aus (a), da $L$ Riemann-Summen auf Riemann-Summen abbildet und stetig ist. Eine Norm ist stetig
($|\|x\|-\|y\||\le\|x-y\|\le\sum_k|x_k-y_k|\,\|e_k\|$), also ist $t\mapsto\|F(t)\|$ stetig, und (c) folgt aus
(a) und der Dreiecksungleichung für Riemann-Summen. Für (d) nehme man die Zerlegungen in $N$ gleiche Teile: Die
durch $b-a$ geteilten Riemann-Summen sind Konvexkombinationen $\frac1N\sum_iF(\tau_i)\in K$, und $K$ ist abgeschlossen. $\square$

Das rechtfertigt Lemma 8.4(c) (Riemann-Summen von CP-Abbildungen) und die Cesàro-Mittel
$\sigma_\tau=\frac1\tau\int_0^\tau T_s\rho\,ds$ aus Proposition 8.12: $\sigma_\tau$ ist nach (d) ein Zustand, und
$\|T_t\sigma_\tau-\sigma_\tau\|_1\le2t/\tau$ nach (b), (c) und (e).

<a id="s-A-2"></a>

## A.2 Nullstellenmengen analytischer Funktionen

Für einen Multiindex $\alpha\in\mathbb N^n$ schreiben wir $|\alpha|=\alpha_1+\dots+\alpha_n$ und $y^\alpha=y_1^{\alpha_1}\cdots y_n^{\alpha_n}$.
Für $a\in\mathbb R^n$ und $0<\rho\le\infty$ sei $Q_\rho(a)=\{x\in\mathbb R^n:\max_i|x_i-a_i|<\rho\}$ der offene
Würfel ($Q_\infty(a)=\mathbb R^n$).

**Definition.** Eine Funktion $f:Q_\rho(a)\to\mathbb R$ ist eine **Potenzreihe auf $Q_\rho(a)$**, wenn
$f(x)=\sum_\alpha c_\alpha(x-a)^\alpha$ mit $\sum_\alpha|c_\alpha|r^{|\alpha|}<\infty$ für alle $0<r<\rho$.

Polynome sind Potenzreihen auf $\mathbb R^n$ (nur endlich viele $c_\alpha\neq0$); in mehreren Variablen wendet der Kurs
Satz A.6 nur auf Polynome an. (Eine Funktion, die auf einer zusammenhängenden offenen Menge nur lokal eine Potenzreihe ist, würde
zusätzlich den Identitätssatz in mehreren Variablen erfordern, den der Kurs nicht verwendet.) Auf dem abgeschlossenen Würfel vom Radius
$r<\rho$ sind die Terme durch $|c_\alpha|r^{|\alpha|}$ beschränkt, also konvergiert die Reihe dort absolut und gleichmäßig
(Weierstraßsches Majorantenkriterium): Ihre Summe hängt nicht von der Reihenfolge der Terme ab, und $f$ ist stetig.

<a id="r-A-5"></a>

**Lemma A.5 (eine Variable).** Sei $g(t)=\sum_kc_k(t-a)^k$ eine Potenzreihe auf $I=(a-\rho,a+\rho)$.

(a) Für $b\in I$ und $|t-b|<\rho_b=\rho-|b-a|$ ist $g(t)=\sum_jd_j(t-b)^j$ mit
$d_j=\sum_{k\ge j}\binom kjc_k(b-a)^{k-j}$, einer Potenzreihe auf $(b-\rho_b,b+\rho_b)$.

(b) Ist ein $d_j\neq0$, so hat $g$ in einer punktierten Umgebung von $b$ keine Nullstelle.

(c) Häufen sich die Nullstellen von $g$ in einem Punkt von $I$, so ist $g=0$ auf $I$. Ist $g\not\equiv0$ auf $I$, so ist ihre
Nullstellenmenge abzählbar, also eine Lebesgue-Nullmenge. Ist $g\equiv0$ auf $I$, so sind alle $c_k=0$.

*Beweis.* (a) Für $0\le s<\rho_b$ ist
$\sum_k\sum_{j\le k}|c_k|\binom kjs^j|b-a|^{k-j}=\sum_k|c_k|(s+|b-a|)^k<\infty$, also konvergiert die Doppelreihe, die man
durch binomisches Ausmultiplizieren von $(t-a)^k=\big((t-b)+(b-a)\big)^k$ erhält, für $|t-b|\le s$ absolut und darf in
beliebiger Reihenfolge summiert werden; summiert man zuerst über $k$, erhält man die $d_j$, mit $\sum_j|d_j|s^j<\infty$.
(b) Sei $j_0$ das kleinste $j$ mit $d_j\ne0$. Dann ist $g(t)=(t-b)^{j_0}h(t)$ mit
$h(t)=\sum_{j\ge j_0}d_j(t-b)^{j-j_0}$, das nahe $b$ stetig ist (Majorantenkriterium) mit $h(b)=d_{j_0}\neq0$; also ist
$h\neq0$ nahe $b$. (c) Sei $A$ die Menge der $b\in I$, in deren Nähe $g$ identisch verschwindet; $A$ ist offen. Häufen sich die
Nullstellen in $b\in I$, so verschwinden nach (b) alle $d_j$, und nach (a) ist $b\in A$. Ist $b\in I$ Grenzwert von
Punkten aus $A$, so häufen sich die Nullstellen in $b$, also $b\in A$: $A$ ist abgeschlossen in $I$. Da $I$ zusammenhängend ist, ist $A$
leer oder gleich $I$. Ist $g\not\equiv0$, so ist $A=\emptyset$, also ist kein Punkt von $I$ Häufungspunkt von Nullstellen; jedes
kompakte Intervall $J\subset I$ enthält dann nur endlich viele Nullstellen (eine unendliche Teilmenge von $J$ hat nach
Bolzano–Weierstraß einen Häufungspunkt in $J$), und $I$ ist abzählbare Vereinigung solcher $J$. Ist $g\equiv0$, so zeigt
(b) an der Stelle $b=a$ (wo $d_j=c_j$), dass alle $c_k=0$ sind. $\square$

<a id="r-A-6"></a>

**Satz A.6 (Nullstellenmengen sind Nullmengen).** Sei $f$ eine Potenzreihe auf einem Würfel $Q\subseteq\mathbb R^n$. Verschwindet $f$
auf einer Teilmenge von $Q$ mit positivem Lebesgue-Maß, so ist $f=0$ auf $Q$. Äquivalent: Ist $f\not\equiv0$ auf $Q$,
so ist $\{x\in Q:f(x)=0\}$ eine Lebesgue-Nullmenge.

*Beweis.* Induktion über $n$; der Fall $n=1$ ist Lemma A.5(c). Sei $n\ge2$; man schreibe $x=(x',t)$ mit
$x'\in\mathbb R^{n-1}$, $Q=Q'\times I$ und $a=(a',a_n)$ und fasse die Terme nach Potenzen von $t-a_n$ zusammen:
$$
f(x',t)=\sum_kh_k(x')(t-a_n)^k,\qquad h_k(x')=\sum_{\alpha'}c_{(\alpha',k)}(x'-a')^{\alpha'} .
$$
Das ist zulässig: Für $x'\in Q'$ und $\max_i|x'_i-a'_i|\le r<\rho$ ist
$\sum_k\sum_{\alpha'}|c_{(\alpha',k)}||x'-a'|^{\alpha'}r^k\le\sum_\alpha|c_\alpha|r^{|\alpha|}<\infty$ (mit
$|y|^{\alpha'}=\prod_i|y_i|^{\alpha'_i}$). Dieselbe Abschätzung zeigt, dass jedes $h_k$ eine Potenzreihe auf $Q'$ ist und dass für
festes $x'\in Q'$ die Funktion $g_{x'}(t)=f(x',t)$ eine Potenzreihe auf $I$ mit Koeffizienten $h_k(x')$ ist.

Sei $Z=\{x\in Q:f(x)=0\}$, abgeschlossen in $Q$, da $f$ stetig ist, und es gelte $\lambda_n(Z)>0$. Nach dem Satz von
Tonelli ist $\lambda_n(Z)=\int_{Q'}\lambda_1(Z_{x'})\,dx'$ mit $Z_{x'}=\{t\in I:f(x',t)=0\}$, also hat die Menge
$F=\{x'\in Q':\lambda_1(Z_{x'})>0\}$ positives $(n-1)$-dimensionales Maß. Für $x'\in F$ liefert Lemma A.5(c)
$g_{x'}\equiv0$ auf $I$ und damit $h_k(x')=0$ für alle $k$. Also verschwindet jedes $h_k$ auf $F$, und nach
Induktionsvoraussetzung ist $h_k=0$ auf $Q'$. Folglich ist $f(x',t)=\sum_kh_k(x')(t-a_n)^k=0$ auf $Q$. $\square$

<a id="r-A-7"></a>

**Korollar A.7.** (a) Ein reelles Polynom $p$ auf $\mathbb R^n$, das nicht identisch null ist, verschwindet nur auf einer
Lebesgue-Nullmenge. Ein Polynom mit einem von null verschiedenen Koeffizienten ist nicht identisch null.

(b) Sei $W$ ein endlichdimensionaler reeller Vektorraum und $p:W\to\mathbb R$ ein Polynom in den Koordinaten
irgendeiner Basis (dann ist es eines in jeder Basis), nicht identisch null. Hat ein Zufallsvektor $X$ in $W$ eine Dichte,
so ist $\mathbb P(p(X)=0)=0$.

(c) Ein echter linearer Unterraum $U\subsetneq W$ ist eine Nullmenge, und ein Zufallsvektor mit Dichte liegt fast sicher außerhalb
von $U$. Ist insbesondere $E\neq0$ ein Unterraum von $\mathbb C^n$ und hat $v$ eine Dichte auf $\mathbb C^n$
(z. B. ein Gaußscher Startvektor), so hat $v$ fast sicher eine von null verschiedene Komponente in $E$.

(d) Seien $X$ und $Y$ unabhängige Zufallsvektoren in endlichdimensionalen reellen Räumen $W_1$ und $W_2$, $Y$ mit
Dichte, und $N\subseteq W_1\times W_2$ eine Borelmenge, deren Schnitte $N_x=\{y:(x,y)\in N\}$ für alle $x$ Nullmengen
sind. Dann ist $\mathbb P((X,Y)\in N)=0$.

*Beweis.* (a) Die erste Aussage ist Satz A.6 mit $Q=\mathbb R^n$. Für die zweite Induktion über $n$: Man schreibe
$p(x',t)=\sum_kp_k(x')t^k$; ein $p_k$ hat einen von null verschiedenen Koeffizienten, also ist nach Induktionsvoraussetzung
$p_k(x'_0)\neq0$ für ein $x'_0$ (für $n=1$ sind die $p_k$ Konstanten, und es gibt kein $x'$); dann ist $t\mapsto p(x'_0,t)$ ein von null verschiedenes Polynom in einer Variablen, das nur endlich viele Nullstellen hat.
(b) Ein Basiswechsel ist linear, also bleibt $p$ ein Polynom, und $\{p=0\}$ ist nach (a) und Lemma A.3(a) eine Nullmenge;
dann ist $\mathbb P(p(X)=0)=\int_{\{p=0\}}f\,d\lambda=0$. (c) Man wähle ein von null verschiedenes lineares Funktional $\ell$, das auf $U$
verschwindet, und wende (b) auf $p=\ell$ an. Für das Beispiel ist $v\mapsto P_Ev$ ($P_E$ die Orthogonalprojektion) eine von null
verschiedene $\mathbb R$-lineare Abbildung, und $\{P_Ev=0\}$ ist ihr Kern, ein echter Unterraum. (d) Die Verteilung von $(X,Y)$ ist das Produkt der Verteilungen, also ist nach
dem Satz von Tonelli $\mathbb P((X,Y)\in N)=\int\mathbb P(Y\in N_x)\,d\mathbb P_X(x)=0$. $\square$

*Rechnung:* $\det Z\neq0$ fast sicher für eine Gaußsche Matrix $Z$ — [Tutorial, Abschnitt A.2](../../tutorial/de/A_analysis.md#s-A-2).

Der Kurs wendet Korollar A.7 an auf $\det Z$ für eine Ginibre-Matrix $Z$ ($|\det Z|^2$ ist ein Polynom in den Real- und
Imaginärteilen der Einträge, gleich $1$ bei $Z=1$; Satz 5.3, Abschnitt 1.6), auf die Gram-Determinante aus
Proposition 4.7(a), auf den Rang Gaußscher Matrizen (Kapitel 12, Aufgabe 4) und auf Startvektoren des
Lanczos-Verfahrens (Kapitel 9). Teil (d) ist der Schritt „da $Y$ von $K$ unabhängig ist, liefert der Satz von Tonelli
Wahrscheinlichkeit $0$“ in Proposition 4.7(b).

<a id="s-A-3"></a>

## A.3 Eigenwerte hermitescher Matrizen

Für $A\in\mathrm{Herm}_n$ seien $\lambda_1(A)\le\dots\le\lambda_n(A)$ seine Eigenwerte mit Vielfachheit und
$v_1,\dots,v_n$ eine Orthonormalbasis aus Eigenvektoren (Satz 1.3, Korollar 1.4). Für einen Einheitsvektor $x$ ist
$\langle x,Ax\rangle=\sum_i\lambda_i|\langle v_i,x\rangle|^2$ mit $\sum_i|\langle v_i,x\rangle|^2=1$, also
$\lambda_1(A)\le\langle x,Ax\rangle\le\lambda_n(A)$; und $\|A\|=\max_i|\lambda_i(A)|$ (Korollar 1.4).

<a id="r-A-8"></a>

**Satz A.8 (Courant–Fischer).** Für $A\in\mathrm{Herm}_n$ und $1\le j\le n$ gilt
$$
\lambda_j(A)=\min_{\dim U=j}\ \max_{x\in U,\,|x|=1}\langle x,Ax\rangle=\max_{\dim U=n-j+1}\ \min_{x\in U,\,|x|=1}\langle x,Ax\rangle ,
$$
über die Unterräume $U\subseteq\mathbb C^n$ der angegebenen Dimension. (Die inneren Maxima und Minima existieren, weil die
Einheitssphäre von $U$ kompakt ist; die äußeren werden angenommen.)

*Beweis.* Sei $\dim U=j$ und $W=\mathrm{span}(v_j,\dots,v_n)$. Wegen $\dim U+\dim W=n+1>n$ gibt es einen Einheitsvektor
$x\in U\cap W$, und $\langle x,Ax\rangle=\sum_{i\ge j}\lambda_i|\langle v_i,x\rangle|^2\ge\lambda_j$. Also ist das Maximum über $U$
mindestens $\lambda_j$. Für $U=\mathrm{span}(v_1,\dots,v_j)$ ist $\langle x,Ax\rangle=\sum_{i\le j}\lambda_i|\langle v_i,x\rangle|^2\le\lambda_j$,
mit Gleichheit bei $x=v_j$; also wird das Minimum über $U$ angenommen und ist gleich $\lambda_j$. Für die zweite Formel
wende man die erste auf $-A$ an, dessen Eigenwerte $\lambda_j(-A)=-\lambda_{n+1-j}(A)$ sind:
$-\lambda_j(A)=\lambda_{n+1-j}(-A)=\min_{\dim U=n+1-j}\max_x(-\langle x,Ax\rangle)=-\max_{\dim U=n+1-j}\min_x\langle x,Ax\rangle$. $\square$

*Rechnung:* Satz A.8 für zufällige Unterräume — [Tutorial, Abschnitt A.3](../../tutorial/de/A_analysis.md#s-A-3).

<a id="r-A-9"></a>

**Korollar A.9 (Kompressionen).** Sei $Q\in M_{n\times m}(\mathbb C)$ mit $Q^*Q=1_m$ und $A\in\mathrm{Herm}_n$. Dann ist
$\lambda_j(Q^*AQ)=\min\max\langle x,Ax\rangle$ über die $j$-dimensionalen $U\subseteq\mathrm{ran}\,Q$ und die Einheitsvektoren $x\in U$, und
$$
\lambda_j(A)\le\lambda_j(Q^*AQ)\le\lambda_{j+n-m}(A)\qquad(1\le j\le m).
$$

*Beweis.* $U'\mapsto QU'$ ist eine Bijektion von den $j$-dimensionalen Unterräumen des $\mathbb C^m$ auf die von
$\mathrm{ran}\,Q$, $|Qy|=|y|$ und $\langle y,Q^*AQy\rangle=\langle Qy,AQy\rangle$; daher liefert Satz A.8 für $Q^*AQ$ die
Formel. Minimiert man über die größere Familie aller $j$-dimensionalen $U\subseteq\mathbb C^n$, erhält man $\lambda_j(A)$, das
somit kleiner ist. Ebenso läuft die Max-Min-Formel für $\lambda_j(Q^*AQ)$ über die $(m-j+1)$-dimensionalen
$U\subseteq\mathrm{ran}\,Q$, und über alle $U\subseteq\mathbb C^n$ dieser Dimension liefert sie die größere Zahl
$\lambda_{n-(m-j+1)+1}(A)=\lambda_{j+n-m}(A)$. $\square$

*Rechnung:* Die Verschränkung der Eigenwerte für eine zufällige Kompression — [Tutorial, Abschnitt A.3](../../tutorial/de/A_analysis.md#s-A-3).

<a id="r-A-10"></a>

**Satz A.10 (Weyl).** Für $A,E\in\mathrm{Herm}_n$ und alle $j$ gilt
$$
\lambda_j(A)+\lambda_1(E)\le\lambda_j(A+E)\le\lambda_j(A)+\lambda_n(E),
\qquad\text{hence}\qquad
|\lambda_j(A+E)-\lambda_j(A)|\le\|E\|\le\|E\|_2 .
$$

*Beweis.* Für Einheitsvektoren $x$ ist $\langle x,(A+E)x\rangle\le\langle x,Ax\rangle+\lambda_n(E)$. Man nehme das Maximum über die
Einheitsvektoren eines $j$-dimensionalen $U$ und dann das Minimum über $U$ (Satz A.8): $\lambda_j(A+E)\le\lambda_j(A)+\lambda_n(E)$. Angewandt auf
$A+E$ und $-E$ ergibt das $\lambda_j(A)\le\lambda_j(A+E)+\lambda_n(-E)=\lambda_j(A+E)-\lambda_1(E)$. Folglich ist
$|\lambda_j(A+E)-\lambda_j(A)|\le\max(|\lambda_1(E)|,|\lambda_n(E)|)=\|E\|$, und $\|E\|^2=\max_i\lambda_i(E)^2\le\sum_i\lambda_i(E)^2=\|E\|_2^2$. $\square$

*Rechnung:* Satz A.10 für zufällige hermitesche Matrizen — [Tutorial, Abschnitt A.3](../../tutorial/de/A_analysis.md#s-A-3).

<a id="r-A-11"></a>

**Korollar A.11.** Jedes $\lambda_j:\mathrm{Herm}_n\,\to\mathbb R$ ist $1$-Lipschitz-stetig bezüglich $\|\cdot\|$, also stetig und
Borel-messbar. Für eine zufällige hermitesche Matrix $H$ ist $(\lambda_1(H),\dots,\lambda_n(H))$ ein Zufallsvektor, ebenso seine Abstände und
Abstandsverhältnisse (auf der Borelmenge, auf der ihre Nenner von null verschieden sind); die Eigenwertvektoren unabhängiger zufälliger hermitescher Matrizen sind unabhängig. Die Matrizen mit
$n$ verschiedenen Eigenwerten bilden eine offene Menge.

*Beweis.* Satz A.10 mit $E=B-A$. Funktionen unabhängiger Zufallsvariablen sind unabhängig. Ist
$\lambda_j(A)<\lambda_{j+1}(A)$ für alle $j$, so bleiben die strikten Ungleichungen für $\|B-A\|<\frac12\min_j(\lambda_{j+1}(A)-\lambda_j(A))$ erhalten. $\square$

<a id="r-A-12"></a>

**Lemma A.12 (Lagrange-Formel).** Sei $A\in\mathrm{Herm}_n$ mit den verschiedenen Eigenwerten $\mu_1<\dots<\mu_r$ und den
Spektralprojektionen $P_1,\dots,P_r$ (Orthogonalprojektionen auf die Eigenräume, also $A=\sum_i\mu_iP_i$,
$P_iP_l=\delta_{il}P_i$, $\sum_iP_i=1$). Dann ist $p(A)=\sum_ip(\mu_i)P_i$ für jedes Polynom $p$, und
$$
P_i=\prod_{l\neq i}\frac{A-\mu_l}{\mu_i-\mu_l} .
$$

*Beweis.* $A^k=\sum_i\mu_i^kP_i$ durch Induktion über $k$ (mit $P_iP_l=\delta_{il}P_i$), und $p(A)=\sum_ip(\mu_i)P_i$ wegen der Linearität.
Das Polynom $\ell_i(t)=\prod_{l\neq i}(t-\mu_l)/(\mu_i-\mu_l)$ erfüllt $\ell_i(\mu_l)=\delta_{il}$, also ist $\ell_i(A)=P_i$. $\square$

*Rechnung:* Die Lagrange-Formel für Eigenwerte der Vielfachheiten $2,3,1$ — [Tutorial, Abschnitt A.3](../../tutorial/de/A_analysis.md#s-A-3).

<a id="r-A-13"></a>

**Proposition A.13 (Stetigkeit der Eigenwerte und Spektralprojektionen).** Sei $S$ eine Teilmenge eines endlichdimensionalen
reellen Vektorraums, $A:S\to\mathrm{Herm}_n$ stetig, und $A(x)$ habe für jedes $x\in S$ genau $r$ verschiedene Eigenwerte
$\mu_1(x)<\dots<\mu_r(x)$. Dann sind ihre Vielfachheiten auf $S$ lokal konstant, und die $\mu_i$ sowie
die Spektralprojektionen $P_i(x)$ von $A(x)$ sind stetig auf $S$. Insbesondere sind auf der offenen Menge der hermiteschen
Matrizen mit einfachem Spektrum die Eigenwerte und die Spektralprojektionen vom Rang eins stetig.

*Beweis.* Ist $r=1$, so ist $A(x)=\mu_1(x)1$ mit $\mu_1=\mathrm{Tr}\,A/n$ und $P_1=1$. Sei $r\ge2$, man fixiere $x_0\in S$ und setze
$d=\min_i(\mu_{i+1}(x_0)-\mu_i(x_0))>0$. Sei $J_i$ die Menge der Indizes $j$ mit $\lambda_j(A(x_0))=\mu_i(x_0)$, ein Block
aufeinanderfolgender ganzer Zahlen. Wegen der Stetigkeit von $A$ gibt es eine Umgebung $N$ von $x_0$ in $S$ mit $\|A(x)-A(x_0)\|<d/2$ für
$x\in N$; nach Satz A.10 liegt $\lambda_j(A(x))$ für $j\in J_i$ im Intervall $I_i=(\mu_i(x_0)-d/2,\mu_i(x_0)+d/2)$, und
diese $r$ Intervalle sind disjunkt. Also hat $A(x)$ in jedem $I_i$ mindestens einen der verschiedenen Eigenwerte und, da es genau $r$ hat,
genau einen: $\lambda_j(A(x))$ ist für alle $j\in J_i$ derselbe. Daher sind die Vielfachheiten auf $N$ die Konstanten $|J_i|$, und
$\mu_i(x)=\lambda_{\max J_i}(A(x))$ ist stetig (Korollar A.11). Nach Lemma A.12 ist
$P_i(x)=\prod_{l\neq i}(A(x)-\mu_l(x))/(\mu_i(x)-\mu_l(x))$, stetig auf $N$, weil die Nenner von $0$ weg beschränkt bleiben.
Die letzte Aussage ist der Fall $r=n$, mit $S$ der offenen Menge aus Korollar A.11. $\square$

In Proposition 4.7 ist $S$ die Menge, auf der die Gram-Determinante $\delta$ von null verschieden ist und auf der $K$ genau $N$ verschiedene
Eigenwerte hat; also ist $K\mapsto\Pi_i(K)$ dort stetig, und das Ereignis aus Proposition 4.7(b) ist eine Borelmenge. Satz 5.15
verwendet Korollar A.11, und Satz 9.14 verwendet Korollar A.9.

<a id="s-A-4"></a>

## A.4 Gauß-Vektoren

Ein **Standard-Gauß-Vektor** im $\mathbb R^m$ ist $g=(g_1,\dots,g_m)$ mit unabhängigen $N(0,1)$-verteilten Koordinaten. Seine Verteilung
ist das Produkt der Verteilungen der Koordinaten, also hat er nach dem Satz von Tonelli die Dichte
$\varphi_m(x)=(2\pi)^{-m/2}e^{-|x|^2/2}$.

<a id="r-A-14"></a>

**Proposition A.14 (orthogonale Invarianz).** (a) $Og$ hat die Verteilung von $g$ für jedes orthogonale $O\in O(m)$.

(b) Sei $W$ ein reeller Skalarproduktraum der Dimension $m$. Für eine Orthonormalbasis $(w_i)$ hängt die Verteilung von
$\sum_ig_iw_i$ nicht von der Basis ab; einen solchen Vektor nennen wir einen **Standard-Gauß-Vektor von $W$**. Seine Koordinaten
in jeder Orthonormalbasis sind unabhängig und $N(0,1)$-verteilt, seine Verteilung ist invariant unter den orthogonalen Abbildungen von $W$, und ist
$W=W_1\oplus W_2$ eine orthogonale Zerlegung, so sind seine Komponenten unabhängige Standard-Gauß-Vektoren von $W_1$ und $W_2$.

(c) Für $v\in\mathbb R^m$ hat $\langle v,g\rangle$ die Verteilung $N(0,|v|^2)$; insbesondere ist $\sum_ia_ig_i\sim N(0,\sum_ia_i^2)$.

*Beweis.* (a) Nach Lemma A.3(c) hat $Og$ die Dichte $\varphi_m(O^{-1}y)/|\det O|=\varphi_m(y)$. (b) Ist $w'_i=\sum_jO_{ji}w_j$
eine weitere Orthonormalbasis, so ist $O$ orthogonal und $\sum_ig_iw'_i=\sum_j(Og)_jw_j$, was nach (a) die Verteilung von
$\sum_jg_jw_j$ hat. Die Koordinaten von $X=\sum_ig_iw_i$ in einer weiteren Orthonormalbasis $(w'_i)$ sind $O^Tg$, nach (a) wieder
standard. Eine orthogonale Abbildung $L$ bildet $(w_i)$ auf eine Orthonormalbasis ab, und $LX=\sum_ig_iLw_i$. Für die letzte
Aussage wähle man $(w_i)$ angepasst an $W_1\oplus W_2$. (c) Für $v=0$ sind beide Seiten das Punktmaß in $0$. Für $v\neq0$
ergänze man $v/|v|$ zu einer Orthonormalbasis; nach (b) ist die Koordinate $\langle v/|v|,g\rangle$ $N(0,1)$-verteilt, also hat
$\langle v,g\rangle=|v|\,\langle v/|v|,g\rangle$ die Dichte $(2\pi|v|^2)^{-1/2}e^{-y^2/2|v|^2}$ (Lemma A.3(c) mit $m=1$). $\square$

**Definition.** Ein **zentrierter Gauß-Vektor** in einem endlichdimensionalen reellen Vektorraum $W$ ist $X=Ag$ für eine lineare
Abbildung $A:\mathbb R^M\to W$ und einen Standard-Gauß-Vektor $g$ des $\mathbb R^M$. In Koordinaten ($W=\mathbb R^m$) ist seine
**Kovarianz** $C=\mathbb E[XX^T]=AA^T$. Er heißt **nicht ausgeartet**, wenn $A$ surjektiv ist, äquivalent, wenn $C$ positiv definit ist
($\langle x,Cx\rangle=|A^Tx|^2$).

Jeder Zufallsvektor des Kurses, der linear aus unabhängigen Gaußschen Einträgen aufgebaut ist – komplexe Gauß-Vektoren,
Ginibre-Matrizen, GOE- und GUE-Matrizen, ihre Bilder unter Twirls –, ist ein zentrierter Gauß-Vektor.

<a id="r-A-15"></a>

**Satz A.15 (lineare Bilder und Dichten).** Sei $X=Ag$ ein zentrierter Gauß-Vektor in $W$.

(a) Für jede lineare Abbildung $L:W\to W'$ ist $LX$ ein zentrierter Gauß-Vektor; er ist nicht ausgeartet, wenn $X$ es ist und $L$ surjektiv ist.

(b) Ist $X$ nicht ausgeartet, so hat er eine Dichte bezüglich des Lebesgue-Maßes auf $W$; in Koordinaten
$W=\mathbb R^m$ ist
$$
f_X(x)=(2\pi)^{-m/2}(\det C)^{-1/2}\,e^{-\langle x,C^{-1}x\rangle/2}.
$$
Ist $X$ ausgeartet, so liegt er im echten Unterraum $\mathrm{ran}\,A$, einer Nullmenge (Korollar A.7(c)), und hat keine Dichte.

(c) Die Verteilung von $X$ hängt nur von seiner Kovarianz ab.

*Beweis.* (a) $LX=(LA)g$, und $LA$ ist surjektiv, wenn $A$ und $L$ es sind. (b) Sei $K=\ker A$, also
$\mathbb R^M=K^\perp\oplus K$. Nach Proposition A.14(b) ist $g=g_1+g_2$ mit einem Standard-Gauß-Vektor $g_1$ von
$K^\perp$, und $X=Ag_1=Bg_1$ mit $B=A|_{K^\perp}$, das injektiv und, da $A$ surjektiv ist, bijektiv
$K^\perp\to W$ ist; also ist $\dim K^\perp=m$. In einer Orthonormalbasis von $K^\perp$, ergänzt durch eine von $K$, ist die Matrix von $A$ gleich
$(B\ 0)$, also $C=AA^T=BB^T$. Die Koordinaten von $g_1$ haben die Dichte $\varphi_m$, und nach Lemma A.3(c) hat $X$
die Dichte $\varphi_m(B^{-1}x)/|\det B|$, mit $|B^{-1}x|^2=\langle x,(BB^T)^{-1}x\rangle$ und $|\det B|=(\det C)^{1/2}$. Wäre
$X$ ausgeartet und hätte eine Dichte $f$, so wäre $1=\mathbb P(X\in\mathrm{ran}\,A)=\int_{\mathrm{ran}\,A}f\,d\lambda=0$.
(c) $\ker C=\ker A^T$ (aus $AA^Tx=0$ folgt $|A^Tx|^2=0$), also enthält $\mathrm{ran}\,C=(\ker C)^\perp=(\ker A^T)^\perp=\mathrm{ran}\,A$
den Vektor $X$. Sei $d=\mathrm{rank}\,C$ (ist $d=0$, so ist $X=0$); man wähle eine Orthonormalbasis $p_1,\dots,p_d$ von
$\mathrm{ran}\,C$ (sie hängt nur von $C$ ab) und sei $P$ die $d\times m$-Matrix mit den Zeilen $p_k^T$, sodass $PP^T=1_d$ und
$P^TP$ die Orthogonalprojektion auf $\mathrm{ran}\,C$ ist. Dann ist $PX=(PA)g$ mit $PA$ surjektiv auf $\mathbb R^d$, ein nicht ausgearteter
zentrierter Gauß-Vektor mit Kovarianz $PCP^T$; nach (b) hängt seine Verteilung nur von $C$ ab, und damit auch die Verteilung von
$X=P^T(PX)$. $\square$

*Rechnung:* Die Dichte $f_X$ im Vergleich mit einer Monte-Carlo-Schätzung — [Tutorial, Abschnitt A.4](../../tutorial/de/A_analysis.md#s-A-4).

So ist etwa in Kapitel 4 der Twirl $P$ eine lineare Projektion von $M_n(\mathbb C)$ auf $\mathcal A'$, also sind
$Z\mapsto P(Z)$ und $Z\mapsto P(Z+Z^*)$ $\mathbb R$-lineare Surjektionen auf $\mathcal A'$ bzw. auf $\mathcal A'_h$ (für
$H\in\mathcal A'_h$ ist $P(\frac H2+\frac H2)=H$); für eine Ginibre-Matrix $Z$ (nach Proposition A.16 nicht ausgeartet) haben $Y=P(Z)$
und $K=P(Z+Z^*)$ nach Satz A.15 Dichten auf diesen Räumen. Eine Gaußsche Linearkombination $\sum_k\xi_kB_k$ einer Basis
$(B_k)$ von $\mathcal A'$ ist das Bild des Standard-Gauß-Vektors $\xi$ unter einer linearen Bijektion.

<a id="r-A-16"></a>

**Proposition A.16 (komplexe Gauß-Vektoren, Ginibre-Matrizen).** (a) Sei $z\in\mathbb C^n$ mit unabhängigen
Koordinaten, deren Real- und Imaginärteile unabhängig und $N(0,\frac12)$-verteilt sind (ein **komplexer Standard-Gauß-Vektor**).
Dann ist $\sqrt2z$ der Standard-Gauß-Vektor von $\mathbb C^n\cong\mathbb R^{2n}$ mit dem Skalarprodukt
$\mathrm{Re}\langle u,w\rangle$; $z$ hat die Dichte $\pi^{-n}e^{-|z|^2}$; $Uz$ hat für jedes
$U\in U(n)$ die Verteilung von $z$, also sind die Koordinaten $\langle w_i,z\rangle$ in jeder Orthonormalbasis wieder unabhängige
komplexe Standard-Gauß-Variablen; und $\mathbb E[z_i\bar z_j]=\delta_{ij}$, $\mathbb E[z_iz_j]=0$.

(b) Eine $a\times b$-Ginibre-Matrix $G$ (unabhängige komplexe Standard-Gauß-Einträge) hat die Dichte
$\pi^{-ab}e^{-\|G\|_2^2}$ auf $M_{a\times b}(\mathbb C)$, und $VGW$ hat für alle $V\in U(a)$, $W\in U(b)$ die Verteilung von $G$.

(c) Eine reelle $a\times b$-Matrix $G$ mit unabhängigen $N(0,1)$-Einträgen hat die Dichte $(2\pi)^{-ab/2}e^{-\|G\|_2^2/2}$, und
$OGO'$ hat für alle $O\in O(a)$, $O'\in O(b)$ die Verteilung von $G$.

*Beweis.* (a) Die $2n$ Zahlen $\sqrt2\mathrm{Re}\,z_i,\sqrt2\mathrm{Im}\,z_i$ sind unabhängig und $N(0,1)$-verteilt, und sie sind die
Koordinaten von $\sqrt2z$ in der Orthonormalbasis $(e_i,ie_i)$. Nach Lemma A.3(c) mit $L=2^{-1/2}$ hat $z$ die Dichte
$(2\pi)^{-n}e^{-|\sqrt2z|^2/2}2^{n}=\pi^{-n}e^{-|z|^2}$. Ein unitäres $U$ ist $\mathbb R$-linear und erhält
$\mathrm{Re}\langle u,w\rangle$, ist also orthogonal; man wende Proposition A.14(b) an. Die Koordinaten in der Basis $(w_i)$ sind
$W^*z$, wobei $W$ die unitäre Matrix mit den Spalten $w_i$ ist. Schließlich ist mit $z_j=x_j+iy_j$
$\mathbb E[z_i\bar z_j]=\mathbb E[x_ix_j+y_iy_j]+i\,\mathbb E[y_ix_j-x_iy_j]=\delta_{ij}$ und $\mathbb E[z_iz_j]=\mathbb E[x_ix_j-y_iy_j]+i\,\mathbb E[x_iy_j+y_ix_j]=0$.
(b) $M_{a\times b}(\mathbb C)\cong\mathbb C^{ab}$ mit $|G|=\|G\|_2$, und $X\mapsto VXW$ ist unitär bezüglich
$\mathrm{Tr}(X^*Y)$ (Konventionen nach Lemma A.3); man wende (a) an. (c) Ebenso mit $\mathbb R^{ab}$ und Proposition A.14. $\square$

<a id="r-A-17"></a>

**Proposition A.17 (GOE und GUE).** Sei $H=(G+G^T)/2$ mit einer reellen $n\times n$-Matrix $G$ mit unabhängigen $N(0,1)$-Einträgen
(GOE) oder $H=(G+G^*)/\sqrt2$ mit einer Ginibre-Matrix $G$ (GUE). Dann sind die Einträge $h_{ij}$, $i\le j$, unabhängig,
$h_{ii}\sim N(0,1)$, und für $i<j$ ist $h_{ij}\sim N(0,\frac12)$ (GOE) bzw. sind $\mathrm{Re}\,h_{ij},\mathrm{Im}\,h_{ij}$
unabhängig und $N(0,\frac12)$-verteilt (GUE). Die Dichte von $H$ in den nach Lemma A.3 festgelegten Koordinaten ist
$$
(2\pi)^{-n/2}\pi^{-n(n-1)/4}e^{-\mathrm{Tr}\,H^2/2}\ \ \text{(GOE)},\qquad
(2\pi)^{-n/2}\pi^{-n(n-1)/2}e^{-\mathrm{Tr}\,H^2/2}\ \ \text{(GUE)},
$$
und $OHO^T$ ($O\in O(n)$) bzw. $UHU^*$ ($U\in U(n)$) hat die Verteilung von $H$.

*Beweis.* GOE: $G$ ist der Standard-Gauß-Vektor von $M_n(\mathbb R)$ mit $\mathrm{Tr}(X^TY)$, und die Matrizen
$e_{ii}$, $(e_{ij}\pm e_{ji})/\sqrt2$ ($i<j$) bilden eine Orthonormalbasis. Nach Proposition A.14(b) sind die Koordinaten
$G_{ii}$, $(G_{ij}+G_{ji})/\sqrt2$, $(G_{ij}-G_{ji})/\sqrt2$ unabhängig und $N(0,1)$-verteilt. Nun ist $h_{ii}=G_{ii}$ und
$h_{ij}=2^{-1/2}\cdot(G_{ij}+G_{ji})/\sqrt2$. GUE: $\sqrt2G$ ist der Standard-Gauß-Vektor von $M_n(\mathbb C)$ mit
$\mathrm{Re}\,\mathrm{Tr}(X^*Y)$ (Proposition A.16), mit der Orthonormalbasis $e_{ii}$, $ie_{ii}$,
$(e_{ij}\pm e_{ji})/\sqrt2$, $i(e_{ij}\pm e_{ji})/\sqrt2$; die Koordinaten von $\sqrt2G$ sind $\sqrt2\mathrm{Re}\,G_{ii}$,
$\sqrt2\mathrm{Im}\,G_{ii}$, $\mathrm{Re}(G_{ij}\pm G_{ji})$, $\mathrm{Im}(G_{ij}\pm G_{ji})$, unabhängig und $N(0,1)$-verteilt. Nun ist
$h_{ii}=\sqrt2\mathrm{Re}\,G_{ii}$, $\mathrm{Re}\,h_{ij}=2^{-1/2}\mathrm{Re}(G_{ij}+G_{ji})$ und
$\mathrm{Im}\,h_{ij}=2^{-1/2}\mathrm{Im}(G_{ij}-G_{ji})$. In beiden Fällen ist die Dichte das Produkt der
eindimensionalen Dichten, $(2\pi)^{-1/2}e^{-h_{ii}^2/2}$ und $\pi^{-1/2}e^{-t^2}$ für jede reelle Koordinate $t$ der Varianz
$\frac12$, und $\mathrm{Tr}\,H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2$. Die Abbildung $H\mapsto UHU^*$ ist linear, invertierbar und
orthogonal bezüglich $\mathrm{Tr}(XY)$, also hat $UHU^*$ nach Lemma A.3(b),(c) die Dichte $f(U^*HU)=f(H)$, da $f$ nur von
$\mathrm{Tr}\,H^2$ abhängt; ebenso für $O$. $\square$

*Rechnung:* $\mathbb E\mathrm{Tr}\,H^2$ für GOE und GUE sowie die Invarianz des GOE — [Tutorial, Abschnitt A.4](../../tutorial/de/A_analysis.md#s-A-4).

<a id="s-A-5"></a>

## A.5 Die Richtung eines Gauß-Vektors

Sei $S^{n-1}=\{x\in\mathbb R^n:|x|=1\}$, $V_n=\lambda_n(\{|x|\le1\})$, und für eine Borelmenge $B\subseteq S^{n-1}$ sei
$C(B)=\{x\in\mathbb R^n:0<|x|\le1,\ x/|x|\in B\}$, eine Borelmenge ($x\mapsto(|x|,x/|x|)$ ist stetig auf $\mathbb R^n\setminus0$).

<a id="r-A-18"></a>

**Satz A.18 (Polarkoordinaten).** Für Borelmengen $A\subseteq(0,\infty)$ und $B\subseteq S^{n-1}$ gilt
$$
\lambda_n\big(\{x\neq0:|x|\in A,\ x/|x|\in B\}\big)=n\,\lambda_n(C(B))\int_Ar^{n-1}\,dr .
$$
Folglich gilt für Borel-messbare $h:(0,\infty)\to[0,\infty]$ und $\phi:S^{n-1}\to[0,\infty]$, mit $\mu(B)=\lambda_n(C(B))$,
$$
\int_{\mathbb R^n}h(|x|)\,\phi(x/|x|)\,dx=n\int_0^\infty h(r)r^{n-1}\,dr\int_{S^{n-1}}\phi\,d\mu .
$$

*Beweis.* Man fixiere $B$ und seien $\nu_1(A)$ und $\nu_2(A)$ die beiden Seiten, Borel-Maße auf $(0,\infty)$. Die Streckung
$x\mapsto rx$ hat die Determinante $r^n$ und bildet $C(B)$ auf $\{0<|x|\le r,\ x/|x|\in B\}$ ab, da $rx/|rx|=x/|x|$; also ist
$\nu_1((0,r])=r^n\lambda_n(C(B))=\nu_2((0,r])$. Für $R>0$ sind die Einschränkungen von $\nu_1,\nu_2$ auf $(0,R]$ endlich
(höchstens $\lambda_n(\{|x|\le R\})$) und stimmen auf den Intervallen $(0,r]$, $r\le R$, überein, die abgeschlossen unter Durchschnitten sind,
$(0,R]$ enthalten und die Borelmengen von $(0,R]$ erzeugen (sie liefern $(s,r]=(0,r]\setminus(0,s]$, und jede relativ offene
Teilmenge von $(0,R]$ ist abzählbare Vereinigung solcher Intervalle). Nach Lemma A.1 ist $\nu_1=\nu_2$ auf $(0,R]$; man lasse $R\to\infty$. Die
Integralformel gilt für $h=1_A$, $\phi=1_B$ (die eben bewiesene Identität; $\{0\}$ ist eine Nullmenge), also wegen der Linearität für nichtnegative
einfache $h,\phi$, und für allgemeine $h,\phi$ nach dem Satz von der monotonen Konvergenz entlang einfacher Funktionen $h_k\uparrow h$, $\phi_k\uparrow\phi$. $\square$

<a id="r-A-19"></a>

**Korollar A.19 (radiale Integrale).** Für Borel-messbares $h\ge0$ auf $(0,\infty)$ ist
$\int_{\mathbb R^n}h(|x|)\,dx=s_n\int_0^\infty h(r)r^{n-1}\,dr$ mit $s_n=nV_n=(2\pi)^{n/2}/I_n$,
$I_n=\int_0^\infty r^{n-1}e^{-r^2/2}dr$. Dabei ist $I_1=\sqrt{\pi/2}$, $I_2=1$, $I_{n+2}=nI_n$, also $s_1=2$, $s_2=2\pi$, $s_3=4\pi$,
$s_4=2\pi^2$; allgemein $s_n=2\pi^{n/2}/\Gamma(n/2)$.

*Beweis.* Die Formel ist Satz A.18 mit $\phi=1$, $\mu(S^{n-1})=V_n$. Für $h(r)=e^{-r^2/2}$ ist die linke Seite
$\big(\int_{\mathbb R}e^{-t^2/2}dt\big)^n=(2\pi)^{n/2}$ (Fubini), was $s_n$ ergibt. $I_1=\frac12\sqrt{2\pi}$, $I_2=1$, und
partielle Integration liefert $I_{n+2}=\big[-r^ne^{-r^2/2}\big]_0^\infty+n\int_0^\infty r^{n-1}e^{-r^2/2}dr=nI_n$; die Substitution
$r=\sqrt{2u}$ liefert $I_n=2^{n/2-1}\Gamma(n/2)$. $\square$

*Rechnung:* $V_3=4\pi/3$ per Monte-Carlo — [Tutorial, Abschnitt A.5](../../tutorial/de/A_analysis.md#s-A-5).

Für $n=2$ und $n=3$ sind das die Polar- und Kugelkoordinaten aus Proposition 5.13.

**Definition.** Die **Gleichverteilung** auf $S^{n-1}$ ist $\sigma_{n-1}(B)=\lambda_n(C(B))/V_n$. Auf der Einheitssphäre von
$\mathbb C^N\cong\mathbb R^{2N}$ ist sie $\sigma_{2N-1}$.

<a id="r-A-20"></a>

**Satz A.20 (Polarzerlegung radialer Verteilungen).** Sei $X$ ein Zufallsvektor im $\mathbb R^n$ mit einer Dichte der
Form $f(|x|)$. Dann ist $X\neq0$ fast sicher, $|X|$ hat die Dichte $s_nr^{n-1}f(r)$ auf $(0,\infty)$, $X/|X|$ hat die Verteilung
$\sigma_{n-1}$, und $|X|$ und $X/|X|$ sind unabhängig. Außerdem ist $\sigma_{n-1}(OB)=\sigma_{n-1}(B)$ für jedes $O\in O(n)$.
Insbesondere gilt:

(a) für einen Standard-Gauß-Vektor $g$ des $\mathbb R^n$ ist $g/|g|\sim\sigma_{n-1}$, unabhängig von $|g|$;

(b) für einen komplexen Standard-Gauß-Vektor $z\in\mathbb C^N$ ist $\psi=z/|z|\sim\sigma_{2N-1}$, eine unter $U(N)$ invariante Verteilung,
unabhängig von $|z|$; und $|z|^2$ hat die Dichte $t^{N-1}e^{-t}/(N-1)!$ auf $(0,\infty)$, also ist $\mathbb E|z|^4=N(N+1)$.

*Beweis.* Nach Satz A.18 mit $h=f\cdot1_A$, $\phi=1_B$ ist
$\mathbb P(|X|\in A,\ X/|X|\in B)=n\,\lambda_n(C(B))\int_Af(r)r^{n-1}dr$. Mit $A=(0,\infty)$ und $B=S^{n-1}$ ist das
$1=nV_n\int_0^\infty f(r)r^{n-1}dr$ (und $\mathbb P(X=0)=0$). Die Wahl $A=(0,\infty)$ liefert $\mathbb P(X/|X|\in B)=\sigma_{n-1}(B)$,
die Wahl $B=S^{n-1}$ liefert $\mathbb P(|X|\in A)=s_n\int_Af(r)r^{n-1}dr$, und der allgemeine Fall ist das Produkt der beiden.
Für $O\in O(n)$ ist $OC(B)=C(OB)$ und $\lambda_n(OE)=\lambda_n(E)$ (Lemma A.3(b)). (a) $\varphi_n(x)=(2\pi)^{-n/2}e^{-|x|^2/2}$ ist radial.
(b) $z=g/\sqrt2$ mit einem Standard-Gauß-Vektor $g$ im $\mathbb R^{2N}$ (Proposition A.16), also $z/|z|=g/|g|$ und $|z|=|g|/\sqrt2$; unitäre Matrizen sind
orthogonal auf $\mathbb R^{2N}$. Die Dichte von $|g|$ ist $s_{2N}(2\pi)^{-N}r^{2N-1}e^{-r^2/2}$ mit $s_{2N}=2\pi^N/(N-1)!$, und die
Substitution $t=r^2/2$ liefert die Dichte von $|z|^2$; dann ist $\mathbb E|z|^4=\int_0^\infty t^{N+1}e^{-t}dt/(N-1)!=(N+1)N$. $\square$

Lemma 7.9 ist Teil (b). Für $N=1$ ist $\sigma_1$ die normierte Bogenlänge, d. h. das Bildmaß von
$\frac1{2\pi}\lambda_1$ auf $[0,2\pi)$ unter $\theta\mapsto e^{i\theta}$: Für eine Borelmenge $B\subseteq S^1$ sei
$\Theta(B)=\{\theta\in(0,2\pi):e^{i\theta}\in B\}$. Die Abbildung $(r,\theta)\mapsto re^{i\theta}$ ist eine $C^1$-Bijektion der offenen Menge
$(0,1)\times(0,2\pi)$ auf die offene Einheitskreisscheibe ohne die Strecke $[0,1)$, mit Jacobi-Determinante $r$, und $C(B)$ unterscheidet sich vom
Bild von $(0,1)\times\Theta(B)$ um eine Teilmenge dieser Strecke und des Einheitskreises, also um eine Nullmenge. Nach der
Transformationsformel ist $\lambda_2(C(B))=\int_{\Theta(B)}\int_0^1r\,dr\,d\theta=\frac12\lambda_1(\Theta(B))$, und mit $V_2=\pi$ folgt
$\sigma_1(B)=\lambda_1(\Theta(B))/2\pi$ (der Punkt $\theta=0$ ist eine Nullmenge).

*Rechnung:* Satz A.20(b) in $\mathbb C^3$ und die gleichverteilte Phase von $z_1$ — [Tutorial, Abschnitt A.5](../../tutorial/de/A_analysis.md#s-A-5).

<a id="r-A-21"></a>

**Korollar A.21 (Momente eines gleichverteilten Einheitsvektors).** Ist $u\sim\sigma_{2N-1}$ auf der Einheitssphäre von $\mathbb C^N$, so ist
$\mathbb E[u_a\bar u_b]=\delta_{ab}/N$, $\mathbb E|u_1|^4=\frac2{N(N+1)}$ und $\mathbb E|u_1|^2|u_2|^2=\frac1{N(N+1)}$ ($N\ge2$). Ist
$u\sim\sigma_{n-1}$ auf $S^{n-1}\subset\mathbb R^n$, so ist $\mathbb Eu_1^2=\frac1n$, $\mathbb Eu_1^4=\frac3{n(n+2)}$ und
$\mathbb Eu_1^2u_2^2=\frac1{n(n+2)}$ ($n\ge2$).

*Beweis.* Man nehme $u=z/|z|$ wie in Satz A.20(b). Da $z_a\bar z_b=|z|^2u_a\bar u_b$ und $|z|$, $u$ unabhängig sind, ist
$\delta_{ab}=\mathbb E[z_a\bar z_b]=\mathbb E|z|^2\,\mathbb E[u_a\bar u_b]=N\,\mathbb E[u_a\bar u_b]$ ($\mathbb E|z|^2=N$ nach
Proposition A.16(a)). Ebenso ist $\mathbb E|z_1|^4=N(N+1)\mathbb E|u_1|^4$ und
$\mathbb E|z_1|^2|z_2|^2=N(N+1)\mathbb E|u_1|^2|u_2|^2$, wobei $\mathbb E|z_1|^2|z_2|^2=1$ wegen der Unabhängigkeit und
$\mathbb E|z_1|^4=\mathbb E(x^2+y^2)^2=\frac34+2\cdot\frac14+\frac34=2$ für unabhängige $x,y\sim N(0,\frac12)$. (Für
$Y\sim N(0,s^2)$ liefert partielle Integration $\mathbb EY^4=3s^2\,\mathbb EY^2=3s^4$, da
$y^3\cdot ye^{-y^2/2s^2}=-s^2y^3\frac d{dy}e^{-y^2/2s^2}$.) Im reellen Fall ist $\mathbb E|g|^2=n$,
$\mathbb E|g|^4=s_n(2\pi)^{-n/2}I_{n+4}=I_{n+4}/I_n=n(n+2)$ (Korollar A.19), $\mathbb Eg_1^4=3$, $\mathbb Eg_1^2g_2^2=1$, und
$|g|^2$, $|g|^4$ sind unabhängig von $u=g/|g|$ (Satz A.20(a)). $\square$

*Rechnung:* Die Momente eines gleichverteilten Einheitsvektors — [Tutorial, Abschnitt A.5](../../tutorial/de/A_analysis.md#s-A-5).

Diese Werte stimmen mit Lemma 6.14 überein, das sie allein aus der unitären Invarianz herleitet. Umgekehrt bestimmt die Invarianz die Verteilung:

<a id="r-A-22"></a>

**Proposition A.22 (Eindeutigkeit der invarianten Verteilung).** Sei $u$ ein zufälliger Einheitsvektor in $\mathbb C^N$, dessen Verteilung
unter allen unitären Matrizen invariant ist. Dann ist $u\sim\sigma_{2N-1}$. Für $N=1$: Jedes rotationsinvariante Borel-Wahrscheinlichkeitsmaß
auf dem Kreis ist die normierte Bogenlänge.

*Beweis.* Sei $Q$ Haar-verteilt auf $U(N)$ (Satz 5.3) und unabhängig von $u$. Für einen festen Einheitsvektor $v$ wähle man
$W\in U(N)$ mit $We_1=v$. Nach Lemma 5.1 ist $QW$ Haar-verteilt, also hat $Qv=(QW)e_1$ die Verteilung von $Qe_1$. In
Satz 5.3 ist $Q$ der Gram-Schmidt-Faktor einer Ginibre-Matrix (Lemma 5.2), also ist $Qe_1=z/|z|$ für ihre erste Spalte $z$, die nach
Satz A.20(b) die Verteilung $\sigma_{2N-1}$ hat. Bedingt man auf $u$ (Fubini, Unabhängigkeit), so ist $\mathbb P(Qu\in B)=\sigma_{2N-1}(B)$
für jede Borelmenge $B$. Bedingt man stattdessen auf $Q$ und verwendet, dass $Uu$ für jedes feste unitäre $U$ die Verteilung von $u$ hat,
so ist $\mathbb P(Qu\in B)=\mathbb P(u\in B)$. Ein Wahrscheinlichkeitsmaß auf dem Kreis ist die Verteilung der identischen Abbildung, und
$U(1)$ wirkt durch Drehungen. $\square$

<a id="s-A-6"></a>

## A.6 Gaußsche Momente: die Isserlis-Wick-Formel

<a id="r-A-23"></a>

**Lemma A.23.** Sei $Y$ ein zentrierter Gauß-Vektor im $\mathbb R^n$ mit Kovarianz $C=(c_{ij})$. Dann ist
$c_{ij}=\mathbb E[Y_iY_j]$, $\mathbb Ee^{\langle t,Y\rangle}=e^{\langle t,Ct\rangle/2}$ und $\mathbb Ee^{|\langle t,Y\rangle|}<\infty$ für alle
$t\in\mathbb R^n$, und alle Momente $\mathbb E|Y_i|^k$ sind endlich.

*Beweis.* Man schreibe $Y=Ag$. Dann ist $\mathbb E[Y_iY_j]=\sum_{k,l}A_{ik}A_{jl}\mathbb E[g_kg_l]=(AA^T)_{ij}$. Mit $a=A^Tt$ ist
$\langle t,Y\rangle=\sum_ka_kg_k$, und wegen der Unabhängigkeit ist
$\mathbb Ee^{\langle t,Y\rangle}=\prod_k\mathbb Ee^{a_kg_k}=\prod_ke^{a_k^2/2}=e^{|A^Tt|^2/2}=e^{\langle t,Ct\rangle/2}$. Schließlich ist
$e^{|s|}\le e^s+e^{-s}$ und $|Y_i|^k\le k!\,e^{|Y_i|}$. $\square$

Eine **Paarung** von $\{1,\dots,2p\}$ ist eine Zerlegung in $p$ zweielementige Blöcke; es gibt $(2p-1)!!=(2p)!/(2^pp!)$ davon
(der Partner von $1$ kann auf $2p-1$ Arten gewählt werden, danach werden die übrigen $2p-2$ Elemente gepaart; Induktion über $p$).

<a id="r-A-24"></a>

**Satz A.24 (Isserlis–Wick).** Seien $Y_1,\dots,Y_m$ Koordinaten eines zentrierten Gauß-Vektors (Wiederholungen
erlaubt). Ist $m$ ungerade, so ist $\mathbb E[Y_1\cdots Y_m]=0$; ist $m=2p$, so gilt
$$
\mathbb E[Y_1\cdots Y_{2p}]=\sum_\pi\prod_{\{a,b\}\in\pi}\mathbb E[Y_aY_b],
$$
mit der Summe über alle Paarungen $\pi$ von $\{1,\dots,2p\}$. Dieselbe Formel gilt für komplexe Zufallsvariablen
$Z_1,\dots,Z_m$, die komplexe Linearkombinationen der Koordinaten eines einzigen zentrierten Gauß-Vektors sind (zum Beispiel die
Einträge einer Ginibre- oder GUE-Matrix und ihre komplex Konjugierten), mit $\mathbb E[Z_aZ_b]$ (ohne Konjugation) in den Paarfaktoren.

*Beweis.* Der Vektor $(Y_1,\dots,Y_m)$ ist das Bild des gegebenen Gauß-Vektors unter einer linearen Abbildung, also ein zentrierter
Gauß-Vektor (Satz A.15(a)); sei $C$ seine Kovarianz. Man fixiere $t\in\mathbb R^m$ und setze $S=\langle t,Y\rangle$. Für $s\in\mathbb R$ sind die
Partialsummen von $\sum_k(sS)^k/k!$ durch $e^{|sS|}$ beschränkt, das nach Lemma A.23 integrierbar ist; nach dem Satz von der majorisierten
Konvergenz und Lemma A.23 ist also
$$
\sum_{k\ge0}\frac{s^k}{k!}\mathbb E[S^k]=\mathbb Ee^{sS}=e^{s^2q/2}=\sum_{p\ge0}\frac{s^{2p}q^p}{2^pp!},\qquad q=\langle t,Ct\rangle .
$$
Beide Seiten sind Potenzreihen auf $\mathbb R$ im Sinne von Abschnitt A.2: Rechts ist das klar, und links ist
$\sum_k\frac{\rho^k}{k!}|\mathbb E[S^k]|\le\mathbb E\sum_k\frac{\rho^k|S|^k}{k!}=\mathbb Ee^{\rho|S|}<\infty$ für jedes $\rho>0$
(monotone Konvergenz und Lemma A.23). Also stimmen ihre Koeffizienten überein (Lemma A.5(c), angewandt auf die Differenz). Folglich ist $\mathbb E\langle t,Y\rangle^m/m!$ gleich $0$ für ungerades $m$ und gleich $\langle t,Ct\rangle^p/(2^pp!)$ für $m=2p$, für jedes
$t\in\mathbb R^m$. Beide Seiten sind Polynome in $t$, also stimmen ihre Koeffizienten überein (Korollar A.7(a), angewandt auf die Differenz).
In $\mathbb E\langle t,Y\rangle^m=\sum_{i_1,\dots,i_m}t_{i_1}\cdots t_{i_m}\mathbb E[Y_{i_1}\cdots Y_{i_m}]$ ($|Y_{i_1}\cdots Y_{i_m}|\le\sum_k|Y_{i_k}|^m$, A.23) stammt das Monom $t_1t_2\cdots t_m$
von den $m!$ Folgen $(i_1,\dots,i_m)$, die Permutationen von $(1,\dots,m)$ sind, also ist sein Koeffizient in
$\mathbb E\langle t,Y\rangle^m/m!$ gleich $\mathbb E[Y_1\cdots Y_m]$. Für ungerades $m$ ist er daher $0$. Für $m=2p$ ist
$q^p=\sum c_{i_1j_1}\cdots c_{i_pj_p}\,t_{i_1}t_{j_1}\cdots t_{i_p}t_{j_p}$ über alle Folgen von $p$ geordneten Paaren; das Monom
$t_1\cdots t_{2p}$ entsteht genau dann, wenn die ungeordneten Paare $\{i_k,j_k\}$ eine Paarung $\pi$ bilden, und jede Paarung entsteht aus
$p!\,2^p$ Folgen (Reihenfolge der Paare, Reihenfolge innerhalb jedes Paares), die jeweils $\prod_{\{a,b\}\in\pi}c_{ab}$ beitragen, da $C$
symmetrisch ist. Division durch $2^pp!$ liefert die Formel, da $c_{ab}=\mathbb E[Y_aY_b]$.
Für die komplexe Aussage schreibe man $Z_a=\sum_i\alpha_{ai}Y_i$. Beide Seiten sind multilinear in $(Z_1,\dots,Z_m)$ – die linke wegen der
Linearität von $\mathbb E$, die rechte, weil jedes $Z_a$ in genau einem Faktor jedes Terms vorkommt –, und sie stimmen überein, wenn jedes
$Z_a$ eines der $Y_i$ ist; man multipliziere aus. $\square$

Zum Beispiel ist $\mathbb E|z|^4=\mathbb E[z\bar zz\bar z]=1\cdot1+0+1\cdot1=2$ für eine komplexe Standard-Gauß-Variable $z$. Proposition 1.10
und Lemma 7.10 verwenden den Fall $m=4$; Satz A.24 beweist ihn, ohne unter dem Erwartungswert zu differenzieren.

*Rechnung:* Satz A.24 für Momente der Ordnung $3$, $4$ und $6$ sowie $\mathbb E|z|^4=2$ — [Tutorial, Abschnitt A.6](../../tutorial/de/A_analysis.md#s-A-6).

<a id="s-A-7"></a>

## A.7 Exponentialverteilte Variablen und gleichverteilte Abstände

Eine Zufallsvariable $e$ heißt **exponentialverteilt** ($e\sim\mathrm{Exp}(1)$), wenn sie die Dichte $e^{-y}$ auf $(0,\infty)$ hat;
$e/\beta$ ist exponentialverteilt mit Rate $\beta$. Sei $\Delta_n=\{x\in\mathbb R^n:0<x_1<\dots<x_n<1\}$.

<a id="r-A-25"></a>

**Satz A.25.** Seien $e_0,\dots,e_n$ ($n\ge1$) unabhängige $\mathrm{Exp}(1)$-verteilte Variablen, $T=e_0+\dots+e_n$ und
$X_k=(e_0+\dots+e_{k-1})/T$ für $k=1,\dots,n$. Dann hat $(X_1,\dots,X_n)$ die Dichte $n!$ auf $\Delta_n$, $T$ hat die Dichte
$t^ne^{-t}/n!$ auf $(0,\infty)$, und $(X_1,\dots,X_n)$ und $T$ sind unabhängig. Insbesondere ist
$(e_0,\dots,e_n)/T$ unabhängig von $T$.

*Beweis.* Die Abbildung $\Phi(x,t)=t\,(x_1,\ x_2-x_1,\ \dots,\ x_n-x_{n-1},\ 1-x_n)$ ist eine Bijektion der offenen Menge
$\Delta_n\times(0,\infty)$ auf $(0,\infty)^{n+1}$ mit der Umkehrung $y\mapsto\big(((y_0+\dots+y_{k-1})/\sum_iy_i)_{k=1}^n,\ \sum_iy_i\big)$;
beide sind $C^1$. In der Jacobi-Matrix von $\Phi$ hat die Spalte von $\partial/\partial t$ die Eintragssumme $1$, und jede Spalte
$\partial/\partial x_j$ (Einträge $t$ und $-t$) hat die Eintragssumme $0$. Ersetzt man die letzte Zeile durch die Summe aller Zeilen (was
die Determinante nicht ändert), erhält man die Zeile $(0,\dots,0,1)$, und der verbleibende $n\times n$-Block ist eine untere Dreiecksmatrix mit
Diagonale $t$; also ist $|\det D\Phi|=t^n$. Der Vektor $e=(e_0,\dots,e_n)$ hat die Dichte $e^{-\sum_iy_i}$, und
$\sum_i\Phi(x,t)_i=t$. Nach der Transformationsformel gilt für Borelmengen $E\subseteq\Delta_n\times(0,\infty)$
$$
\mathbb P\big((X,T)\in E\big)=\mathbb P\big(e\in\Phi(E)\big)=\int_{\Phi(E)}e^{-\sum_iy_i}dy=\int_Ee^{-t}t^n\,dx\,dt .
$$
Die Dichte $e^{-t}t^n=n!\cdot\frac{t^ne^{-t}}{n!}$ ist ein Produkt, und $\int_0^\infty t^ne^{-t}dt=n!$, also ist
$\lambda_n(\Delta_n)=1/n!$; nach Fubini ist $\mathbb P(X\in A,T\in A')=\mathbb P(X\in A)\mathbb P(T\in A')$. Schließlich ist $(e_0,\dots,e_n)/T$ gleich
$(X_1,X_2-X_1,\dots,1-X_n)$, eine Funktion von $X$. $\square$

*Rechnung:* Die Unabhängigkeit von $(e_0,\dots,e_n)/T$ und $T$ — [Tutorial, Abschnitt A.7](../../tutorial/de/A_analysis.md#s-A-7).

<a id="r-A-26"></a>

**Korollar A.26 (gleichverteilte Abstände).** Seien $U_1,\dots,U_n$ unabhängig und gleichverteilt auf $(0,1)$, mit den Ordnungsstatistiken
$U_{(1)}<\dots<U_{(n)}$ und den Abständen $s_0=U_{(1)}$, $s_k=U_{(k+1)}-U_{(k)}$ ($1\le k\le n-1$), $s_n=1-U_{(n)}$. Dann gilt:

(a) Bindungen haben Wahrscheinlichkeit $0$, und $(U_{(1)},\dots,U_{(n)})$ hat die Dichte $n!$ auf $\Delta_n$;

(b) $(s_0,\dots,s_n)$ hat die Verteilung von $(e_0,\dots,e_n)/T$, mit $e_i$, $T$ wie in Satz A.25;

(c) $\mathbb Es_i=\frac1{n+1}$, $\mathbb Es_i^2=\frac2{(n+1)(n+2)}$ und $\mathbb E[s_is_j]=\frac1{(n+1)(n+2)}$ für $i\neq j$;

(d) für $i\neq j$ hat $s_i/s_j$ die Verteilung von $e_i/e_j$, mit $\mathbb P(e_i/e_j\le t)=\frac t{1+t}$ ($t>0$); also hat für
$0\le i\le n-1$ das Abstandsverhältnis $\min(s_i,s_{i+1})/\max(s_i,s_{i+1})$ die in Proposition 5.14(a) bestimmte Verteilung von
$\min(e_i,e_{i+1})/\max(e_i,e_{i+1})$, $\mathbb P(r\le t)=\frac{2t}{1+t}$ auf $[0,1]$, mit Mittelwert $2\ln2-1$.

Für $n$ unabhängige gleichverteilte Punkte auf einem Intervall $[a,b]$ multipliziert die affine Abbildung $u\mapsto a+(b-a)u$ alle Abstände
mit $b-a$: Die durch $b-a$ geteilten Abstände erfüllen (b) und (c), und die Verhältnisse in (d) bleiben unverändert.

*Beweis.* (a) Die Menge, auf der zwei Koordinaten übereinstimmen, ist eine endliche Vereinigung von Hyperebenen, also eine Nullmenge (Korollar A.7(c)). Außerhalb
davon gibt es genau eine Permutation $\pi$ mit $U_{\pi(1)}<\dots<U_{\pi(n)}$; für Borelmengen $E\subseteq\Delta_n$ ist das Ereignis
$\{(U_{(1)},\dots,U_{(n)})\in E\}$ also bis auf eine Nullmenge die disjunkte Vereinigung über $\pi$ der Ereignisse
$\{(U_{\pi(1)},\dots,U_{\pi(n)})\in E\}$. Jedes davon hat die Wahrscheinlichkeit $\lambda_n(E)$, weil $(U_1,\dots,U_n)$ die Dichte
$1_{(0,1)^n}$ hat und das Permutieren der Koordinaten sowohl $\lambda_n$ als auch $(0,1)^n$ erhält. Die Wahrscheinlichkeit ist also $n!\,\lambda_n(E)$.
(b) Nach (a) und Satz A.25 haben die Vektoren $(U_{(k)})_k$ und $(X_k)_k$ dieselbe Verteilung, und die Abstände sind dieselbe lineare
Funktion von beiden. (c) Sei $w=(e_0,\dots,e_n)/T$, unabhängig von $T$ (Satz A.25), mit $\mathbb ET=n+1$ und
$\mathbb ET^2=\int_0^\infty t^{n+2}e^{-t}dt/n!=(n+1)(n+2)$. Aus $e_i=Tw_i$ folgt $1=\mathbb Ee_i=\mathbb ET\,\mathbb Ew_i$,
$2=\mathbb Ee_i^2=\mathbb ET^2\,\mathbb Ew_i^2$ und $1=\mathbb E[e_ie_j]=\mathbb ET^2\,\mathbb E[w_iw_j]$. (d) $s_i/s_j=w_i/w_j=e_i/e_j$, und nach
Tonelli ist $\mathbb P(e_i\le te_j)=\int_0^\infty e^{-y}(1-e^{-ty})\,dy=1-\frac1{1+t}$. Das Abstandsverhältnis ist eine Funktion von
$s_i/s_{i+1}$ (nämlich $\min(x,1/x)$ an der Stelle $x=s_i/s_{i+1}$, fast sicher definiert, da alle $e_k>0$), hat also die Verteilung derselben
Funktion von $e_i/e_{i+1}$, und das ist das Verhältnis $r$ aus Proposition 5.14(a) für die unabhängigen $\mathrm{Exp}(1)$-verteilten
Variablen $e_i,e_{i+1}$. $\square$

*Rechnung:* Korollar A.26(b)–(d) für $n=5$ — [Tutorial, Abschnitt A.7](../../tutorial/de/A_analysis.md#s-A-7).

Man betrachte insbesondere die „Poisson-Niveaus“ $\lambda_k=e_0+\dots+e_{k-1}$, $k=1,\dots,n+1$ (unabhängige exponentialverteilte
Abstände $e_0,\dots,e_n$), und $n$ sortierte gleichverteilte Punkte mit ihren $n+1$ Abständen $s_0,\dots,s_n$. Der Vektor der
Verhältnisse aufeinanderfolgender Abstände $r_i=\min(d_i,d_{i+1})/\max(d_i,d_{i+1})$, $0\le i<n$, hat für
$d=(e_0,\dots,e_n)$ und für $d=(s_0,\dots,s_n)$ dieselbe gemeinsame Verteilung, weil sich die Verhältnisse nicht ändern, wenn alle $d_i$ durch $T$ geteilt werden, und
$(e_0,\dots,e_n)/T$ die Verteilung von $(s_0,\dots,s_n)$ hat (Teil (b)); und nach Satz A.25 haben die ersten $n$ Poisson-Niveaus, geteilt durch das $(n+1)$-te, genau die Verteilung von $n$ sortierten
gleichverteilten Punkten.

**Zurück zum** [Inhaltsverzeichnis](README.md).
