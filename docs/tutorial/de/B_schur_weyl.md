[English](../B_schur_weyl.md) · [Español](../es/B_schur_weyl.md) · [Français](../fr/B_schur_weyl.md) · **Deutsch** · [中文](../zh/B_schur_weyl.md) · [日本語](../ja/B_schur_weyl.md) · [Русский](../ru/B_schur_weyl.md)

# Anhang B. Schur-Weyl-Dualität und der Weingarten-Kalkül

Die unitäre Gruppe $U(d)$ wirkt auf $(\mathbb C^d)^{\otimes n}$ durch $V\mapsto V^{\otimes n}$, die symmetrische Gruppe
$S_n$ durch Permutation der Tensorfaktoren. Jede der beiden Operatorfamilien ist die Kommutante der anderen
(**Schur-Weyl-Dualität**). Dieser Anhang beweist die Dualität vollständig und leitet daraus die Weingarten-Formel für
die Momente einer Haar-verteilten unitären Matrix beliebiger Ordnung her. Nebenbei beweist er die Tatsachen, die die
Kapitel 4, 5, 9 und 11 von hier übernehmen: Satz 4.14(b) für alle $L$, die Irreduzibilität der Darstellungen von
$S_L$ auf den Vielfachheitsräumen des Gesamtspins, die Catalan-Zahlen aus Abschnitt 4.5, Satz 5.12, das
Wigner-Eckart-Theorem für $SU(2)$, das in Abschnitt 9.3 verwendet wird, und das Haar-Maß der orthogonalen Gruppe, das
in Kapitel 11 verwendet wird.

Der Weg ist elementar: Eine Polarisationsidentität identifiziert die Kommutante der Permutationen mit dem Spann aller
$A^{\otimes n}$, $A\in M_d(\mathbb C)$; ein Argument mit der Cayley-Transformation zeigt, dass die unitären
$V^{\otimes n}$ denselben Raum aufspannen; der Bikommutantensatz (Satz 2.2) macht daraus die Dualität. Die
Weingarten-Formel folgt dann, weil das Haar-Mittel die Orthogonalprojektion auf die Kommutante ist.

Durchweg sei $d,n\ge1$ und $\mathcal H_n=(\mathbb C^d)^{\otimes n}$ mit der Orthonormalbasis
$e_{\mathbf i}=e_{i_1}\otimes\dots\otimes e_{i_n}$, wobei ein Multiindex $\mathbf i=(i_1,\dots,i_n)\in[d]^n$ als Abbildung
$[n]\to[d]$, $m\mapsto i_m$ gelesen wird ($[d]=\{1,\dots,d\}$). $\operatorname{End}\mathcal H_n=M_{d^n}(\mathbb C)$ trägt das
Hilbert–Schmidt-Skalarprodukt $\langle X,Y\rangle=\operatorname{Tr}(X^*Y)$. Für $X\in M_d$ und $1\le m\le n$ sei
$X^{(m)}=1\otimes\dots\otimes X\otimes\dots\otimes1$ ($X$ an Position $m$) und $d\Pi(X)=\sum_{m=1}^nX^{(m)}$. Wir schreiben
$$
\mathcal T_n=\operatorname{span}\{A^{\otimes n}:A\in M_d(\mathbb C)\},\qquad
\mathcal U_n=\operatorname{span}\{V^{\otimes n}:V\in U(d)\}.
$$
Produkte von Kronecker-Produkten werden faktorweise berechnet, $(A_1\otimes\dots\otimes A_n)(B_1\otimes\dots\otimes B_n)=A_1B_1\otimes\dots\otimes A_nB_n$,
und $(A_1\otimes\dots\otimes A_n)^*=A_1^*\otimes\dots\otimes A_n^*$. Permutationen werden als Abbildungen verkettet,
$(\sigma\tau)(m)=\sigma(\tau(m))$, und $\#\text{cycles}(\sigma)$ zählt die Zyklen von $\sigma$, Fixpunkte eingeschlossen.

## B.1 Permutationen der Tensorfaktoren

Für $\sigma\in S_n$ sei $\pi(\sigma)e_{\mathbf i}=e_{\mathbf i\circ\sigma^{-1}}$, d. h. der Basisvektor an Stelle $m$ wandert an
Stelle $\sigma(m)$. Das ist die Konvention von `decompose.permutation_representation` und der Kapitel 2 und 4.

**Lemma B.1.** (a) $\pi$ ist eine unitäre Darstellung von $S_n$, und
$\pi(\sigma)(x_1\otimes\dots\otimes x_n)=x_{\sigma^{-1}(1)}\otimes\dots\otimes x_{\sigma^{-1}(n)}$ für alle $x_m\in\mathbb C^d$.
(b) $\pi(\sigma)(A_1\otimes\dots\otimes A_n)\pi(\sigma)^*=A_{\sigma^{-1}(1)}\otimes\dots\otimes A_{\sigma^{-1}(n)}$ für alle
$A_m\in M_d$. Insbesondere vertauscht $\pi(\sigma)$ mit jedem $A^{\otimes n}$ und jedem $d\Pi(X)$, und
$\pi(\sigma)X^{(m)}\pi(\sigma)^*=X^{(\sigma(m))}$.
(c) $\operatorname{Tr}\pi(\sigma)=d^{\#\text{cycles}(\sigma)}$ und $\langle\pi(\sigma),\pi(\tau)\rangle=d^{\#\text{cycles}(\sigma^{-1}\tau)}$.
Die Anzahl der Zyklen ist invariant unter Konjugation und Inversion, und zwei Permutationen sind genau dann konjugiert,
wenn sie denselben Zyklentyp haben.
(d) Die $n!$ Operatoren $\pi(\sigma)$ sind genau dann linear unabhängig, wenn $n\le d$.

*Beweis.* (a) $\mathbf i\mapsto\mathbf i\circ\sigma^{-1}$ ist eine Bijektion von $[d]^n$, also permutiert $\pi(\sigma)$ eine
Orthonormalbasis und ist unitär. $\pi(\sigma)\pi(\tau)e_{\mathbf i}=\pi(\sigma)e_{\mathbf i\circ\tau^{-1}}=e_{\mathbf i\circ\tau^{-1}\circ\sigma^{-1}}=\pi(\sigma\tau)e_{\mathbf i}$.
Mit $x_m=\sum_ax_m(a)e_a$ ist $x_1\otimes\dots\otimes x_n=\sum_{\mathbf i}\prod_mx_m(i_m)\,e_{\mathbf i}$; Anwenden von $\pi(\sigma)$ und
die Substitution $\mathbf j=\mathbf i\circ\sigma^{-1}$ ergeben $\sum_{\mathbf j}\prod_mx_m(j_{\sigma(m)})\,e_{\mathbf j}=\sum_{\mathbf j}\prod_{p}x_{\sigma^{-1}(p)}(j_p)\,e_{\mathbf j}$,
also den behaupteten Produktvektor.
(b) Nach (a), auf einem Produktvektor: $\pi(\sigma)^*=\pi(\sigma^{-1})$ setzt $x_{\sigma(m)}$ an Stelle $m$, der Operator
wendet dort $A_m$ an, und $\pi(\sigma)$ bringt Stelle $m$ an Stelle $\sigma(m)$; an Stelle $p$ steht also am Ende
$A_{\sigma^{-1}(p)}x_p$. Produktvektoren spannen $\mathcal H_n$ auf. Für $X^{(m)}$ wandert der Faktor $X$ von Stelle $m$ an
Stelle $\sigma(m)$, also bleibt $d\Pi(X)$ fest.
(c) $\pi(\sigma)$ ist eine Permutationsmatrix, und $e_{\mathbf i}$ bleibt genau dann fest, wenn $i_{\sigma^{-1}(p)}=i_p$ für
alle $p$, d. h. genau dann, wenn $\mathbf i$ auf den Zyklen von $\sigma$ konstant ist: Es gibt $d^{\#\text{cycles}(\sigma)}$
solche $\mathbf i$. Weiter ist
$\langle\pi(\sigma),\pi(\tau)\rangle=\operatorname{Tr}(\pi(\sigma^{-1})\pi(\tau))=\operatorname{Tr}\pi(\sigma^{-1}\tau)$. Die Zyklen von
$\kappa\rho\kappa^{-1}$ sind die Bilder der Zyklen von $\rho$ unter $\kappa$, und $\rho^{-1}$ durchläuft dieselben Zyklen
rückwärts. Haben umgekehrt $\rho,\rho'$ denselben Zyklentyp, so schreibe man ihre Zyklen, nach Länge geordnet,
untereinander und lasse $\kappa$ jeden Eintrag der ersten Zeile auf den Eintrag darunter abbilden; dann ist
$\kappa\rho\kappa^{-1}=\rho'$.
(d) Sei $n\le d$ und $\mathbf u=(1,2,\dots,n)$. Die Vektoren $\pi(\sigma)e_{\mathbf u}=e_{\mathbf u\circ\sigma^{-1}}$ sind für
verschiedene $\sigma$ verschiedene Basisvektoren, also liefert $\sum_\sigma c_\sigma\pi(\sigma)=0$, angewandt auf
$e_{\mathbf u}$, alle $c_\sigma=0$. Sei $n>d$ und $\mathcal A=\sum_\sigma\operatorname{sgn}(\sigma)\pi(\sigma)$. Jedes
$\mathbf i\in[d]^n$ hat $i_a=i_b$ für gewisse $a\neq b$; für die Transposition $t=(a\,b)$ gilt $\pi(t)e_{\mathbf i}=e_{\mathbf i}$
und $\mathcal A\pi(t)=\sum_\sigma\operatorname{sgn}(\sigma)\pi(\sigma t)=-\mathcal A$. Daher ist
$\mathcal Ae_{\mathbf i}=\mathcal A\pi(t)e_{\mathbf i}=-\mathcal Ae_{\mathbf i}=0$ für alle $\mathbf i$: $\mathcal A=0$ ist eine
nichttriviale lineare Relation. $\square$

## B.2 Polarisation: die Kommutante der Permutationen

**Lemma B.2 (Polarisation).** Für $A_1,\dots,A_n\in M_d(\mathbb C)$ gilt mit $|\varepsilon|=\sum_m\varepsilon_m$
$$
\sum_{\sigma\in S_n}A_{\sigma(1)}\otimes\dots\otimes A_{\sigma(n)}=
\sum_{\varepsilon\in\{0,1\}^n}(-1)^{n-|\varepsilon|}\Big(\sum_{m=1}^n\varepsilon_mA_m\Big)^{\otimes n}.
$$

*Beweis.* Wegen der Multilinearität ist $\big(\sum_m\varepsilon_mA_m\big)^{\otimes n}=\sum_f\varepsilon_{f(1)}\cdots\varepsilon_{f(n)}\,A_{f(1)}\otimes\dots\otimes A_{f(n)}$,
summiert über alle Abbildungen $f:[n]\to[n]$. Da $\varepsilon_m^r=\varepsilon_m$ für $r\ge1$, ist der Koeffizient des
Terms $f$ auf der rechten Seite
$$
\sum_{\varepsilon}\prod_{m=1}^n(-1)^{1-\varepsilon_m}\prod_{m\in f([n])}\varepsilon_m
=\prod_{m\notin f([n])}\Big(\sum_{\varepsilon_m=0,1}(-1)^{1-\varepsilon_m}\Big)\prod_{m\in f([n])}\Big(\sum_{\varepsilon_m=0,1}(-1)^{1-\varepsilon_m}\varepsilon_m\Big).
$$
Faktoren der ersten Art sind $-1+1=0$, die der zweiten Art $1$. Der Koeffizient ist also $1$, wenn $f$ surjektiv,
d. h. eine Permutation, ist, und sonst $0$. $\square$

Der **Symmetrisierer** ist $\mathcal S(X)=\frac1{n!}\sum_\sigma\pi(\sigma)X\pi(\sigma)^*$ auf $\operatorname{End}\mathcal H_n$.

**Proposition B.3.** $\mathcal S$ ist die Orthogonalprojektion von $\operatorname{End}\mathcal H_n$ auf $\pi(S_n)'$, und
$$
\pi(S_n)'=\mathcal T_n,\qquad\dim\pi(S_n)'=\binom{d^2+n-1}n .
$$

*Beweis.* $\mathcal S$ ist der Twirl der endlichen Gruppe $S_n$ in der Darstellung $\pi$, also die Orthogonalprojektion
auf $\pi(S_n)'$ (Satz 12.5 mit $\pi_1=\pi_2=\pi$; siehe auch Abschnitt 4.2). Nach Lemma B.1(b) ist
$\mathcal T_n\subseteq\pi(S_n)'$. Umgekehrt wird $\operatorname{End}\mathcal H_n$ von den Produkten
$A_1\otimes\dots\otimes A_n$ von Matrixeinheiten aufgespannt, und nach Lemma B.1(b) und der Substitution
$\sigma\to\sigma^{-1}$ ist $\mathcal S(A_1\otimes\dots\otimes A_n)=\frac1{n!}\sum_\sigma A_{\sigma(1)}\otimes\dots\otimes A_{\sigma(n)}$,
was nach Lemma B.2 in $\mathcal T_n$ liegt. Also ist $\pi(S_n)'=\mathcal S(\operatorname{End}\mathcal H_n)\subseteq\mathcal T_n$.

*Dimension.* Sei $P=[d]\times[d]$ und $E_\alpha=E_{\alpha_1}\otimes\dots\otimes E_{\alpha_n}$ für $\alpha\in P^n$, mit den
Matrixeinheiten $E_{(a,b)}=|e_a\rangle\langle e_b|$: eine Basis von $\operatorname{End}\mathcal H_n$. Nach Lemma B.1(b) ist
$\mathcal S(E_\alpha)=\frac1{n!}\sum_\sigma E_{\alpha\circ\sigma^{-1}}$, und das ist $|O_\alpha|^{-1}\sum_{\beta\in O_\alpha}E_\beta$ für die
Bahn $O_\alpha=\{\alpha\circ\sigma\}$ (jedes Element der Bahn tritt für dieselbe Anzahl $n!/|O_\alpha|$ von $\sigma$ auf).
Summen über verschiedene Bahnen enthalten disjunkte Mengen von Basiselementen, also hat das Bild von $\mathcal S$ eine
durch die Bahnen indizierte Basis. Zwei Tupel liegen genau dann in derselben Bahn, wenn jedes Element von $P$ in
beiden gleich oft vorkommt: Die Bahnen sind die Multimengen der Größe $n$ aus den $N=d^2$ Elementen von $P$. Nummeriert
man $P$ als $1,\dots,N$, so ist eine solche Multimenge eine Folge $1\le c_1\le\dots\le c_n\le N$, und
$c_m\mapsto c_m+m-1$ ist eine Bijektion auf die streng wachsenden Folgen in $[N+n-1]$, d. h. die $n$-elementigen
Teilmengen; davon gibt es $\binom{N+n-1}n$. $\square$

## B.3 Schur-Weyl-Dualität

**Lemma B.4 (Polynome, die auf der unitären Gruppe verschwinden).** Sei $p:M_d(\mathbb C)\to\mathbb C$ ein Polynom in
den $d^2$ Einträgen, homogen vom Grad $n\ge1$. Ist $p(V)=0$ für alle $V\in U(d)$, so ist $p=0$.

*Beweis.* *Schritt 1 (Cayley-Transformation).* Für $A\in M_d$ sei $q(A)=p\big((1+iA)\operatorname{adj}(1-iA)\big)$, wobei
$\operatorname{adj}$ die adjunkte Matrix bezeichnet; $q$ ist ein Polynom in den Einträgen von $A$. Sei $H$ hermitesch.
Seine Eigenwerte sind reell (Korollar 1.4), also ist $1-iH$ invertierbar, $\operatorname{adj}(1-iH)=\det(1-iH)(1-iH)^{-1}$
(Cramersche Regel), und wegen der Homogenität ist $q(H)=\det(1-iH)^n\,p(C)$ mit $C=(1+iH)(1-iH)^{-1}$. Die Matrizen
$1\pm iH$ und ihre Inversen sind Funktionen von $H$ und vertauschen (Satz 1.5), und $(1+iH)^*=1-iH$, also
$C^*C=(1+iH)^{-1}(1-iH)(1+iH)(1-iH)^{-1}=1$. Daher ist $q(H)=0$.

*Schritt 2 (von hermiteschen zu allen Matrizen).* Man wähle eine Basis $B_1,\dots,B_{d^2}$ des reellen Vektorraums der
hermiteschen Matrizen. Sie ist auch eine Basis von $M_d(\mathbb C)$ über $\mathbb C$, da jedes $A$ gleich
$\frac{A+A^*}2+i\frac{A-A^*}{2i}$ ist und die Dimensionen übereinstimmen. $r(t)=q(\sum_kt_kB_k)$ ist ein Polynom in
$t\in\mathbb C^{d^2}$, das auf $\mathbb R^{d^2}$ verschwindet. Ein solches Polynom hat lauter Koeffizienten $0$: In einer
Variablen hat es unendlich viele Nullstellen; für $N$ Variablen schreibe man $r=\sum_jr_j(t')t_N^j$ – für jedes reelle
$t'$ verschwindet das Polynom in $t_N$ auf $\mathbb R$, also verschwindet jedes $r_j$ auf $\mathbb R^{N-1}$, und Induktion
greift. Daher ist $q(A)=r(t)=0$ für jedes $A=\sum_kt_kB_k\in M_d$.

*Schritt 3 (Umkehrung der Transformation).* Sei $B\in M_d$ mit $1+B$ invertierbar und $A=-i(1+B)^{-1}(B-1)$. Dann ist
$1-iA=1-(1+B)^{-1}(B-1)=2(1+B)^{-1}$ und $1+iA=(1+B)^{-1}2B$; da $B$ mit $1+B$ vertauscht, ist
$(1+iA)(1-iA)^{-1}=(1+B)^{-1}B(1+B)=B$. Also ist $0=q(A)=\det(1-iA)^np(B)$ mit $\det(1-iA)\neq0$: $p(B)=0$.

*Schritt 4.* Für jedes $B\neq0$ und $0<t<1/\|B\|$ ist $-1$ kein Eigenwert von $tB$ (aus $tBv=-v$ folgt
$\|v\|=t\|Bv\|\le t\|B\|\|v\|<\|v\|$, außer wenn $v=0$), also ist $t^np(B)=p(tB)=0$. Schließlich ist $p(0)=0$ wegen
der Homogenität. $\square$

**Lemma B.5 (Newton-Identitäten).** Seien $x_1,\dots,x_n$ paarweise vertauschende Elemente einer Algebra,
$x_S=\prod_{m\in S}x_m$, $e_k=\sum_{|S|=k}x_S$ ($e_0=1$) und $p_j=\sum_mx_m^j$. Dann gilt
$k\,e_k=\sum_{j=1}^k(-1)^{j-1}e_{k-j}\,p_j$ für $1\le k\le n$.

*Beweis.* Für $1\le j\le k$ sei $a_j=\sum x_Sx_m^j$, summiert über die Paare $(S,m)$ mit $|S|=k-j$ und $m\notin S$, und
$a_{k+1}=0$. In $e_{k-j}p_j=\sum_{|S|=k-j}\sum_mx_Sx_m^j$ liefern die Terme mit $m\notin S$ den Wert $a_j$, und die Terme
mit $m\in S$, geschrieben mit $S'=S\setminus\{m\}$, sind $x_{S'}x_m^{j+1}$ und liefern $a_{j+1}$. Die alternierende
Summe teleskopiert also zu $a_1$, und $a_1=\sum_{|T|=k}|T|\,x_T=k\,e_k$ (jede $k$-Menge $T$ entsteht aus den $k$
Wahlen von $m\in T$). $\square$

**Satz B.6 (Schur-Weyl-Dualität).** Sei $d,n\ge1$.
(a) $\pi(S_n)'=\mathcal T_n=\mathcal U_n=\operatorname{span}\{V^{\otimes n}:V\in SU(d)\}$. Dieser Raum ist eine $*$-Algebra, also
die von den $V^{\otimes n}$, $V\in U(d)$, erzeugte Algebra.
(b) Für jede Menge $\mathcal X\subseteq M_d$, deren Spann $M_d$ ist, ist die von $\{d\Pi(X):X\in\mathcal X\}$ erzeugte
unitale Algebra gleich $\mathcal T_n$.
(c) $\{V^{\otimes n}:V\in U(d)\}'=\{V^{\otimes n}:V\in SU(d)\}'=\{d\Pi(X):X\in M_d\}'=\operatorname{span}\pi(S_n)$.

*Beweis.* (a) $\pi(S_n)'=\mathcal T_n$ ist Proposition B.3, und $\mathcal U_n\subseteq\mathcal T_n$. Angenommen,
$\mathcal U_n\neq\mathcal T_n$. Dann ist ein $A_0^{\otimes n}\notin\mathcal U_n$; sei $Y\neq0$ die zu $\mathcal U_n$ orthogonale
Komponente von $A_0^{\otimes n}$. Die Funktion $p(A)=\langle Y,A^{\otimes n}\rangle$ ist ein homogenes Polynom vom Grad $n$
in den Einträgen von $A$ (jeder Eintrag von $A^{\otimes n}$ ist ein Produkt von $n$ Einträgen von $A$), sie verschwindet
auf $U(d)$, und $p(A_0)=\|Y\|^2\neq0$ – ein Widerspruch zu Lemma B.4.
Für $V\in U(d)$ wähle man $\theta$ mit $\det V=e^{i\theta}$; dann ist $W=e^{-i\theta/d}V\in SU(d)$ und
$V^{\otimes n}=e^{in\theta/d}W^{\otimes n}$, also stimmen die beiden Spanne überein. $\mathcal U_n$ enthält $1=1^{\otimes n}$,
$V^{\otimes n}W^{\otimes n}=(VW)^{\otimes n}$ und $(V^{\otimes n})^*=(V^*)^{\otimes n}$, ist also eine $*$-Algebra, und jede
Algebra, die die $V^{\otimes n}$ enthält, enthält ihren Spann.

(b) Sei $\mathcal A$ diese Algebra. $d\Pi$ ist linear, also enthält $\mathcal A$ die $d\Pi(X)$ für alle $X\in M_d$. Nach
Lemma B.1(b) liegt jedes $d\Pi(X)$ in der unitalen Algebra $\pi(S_n)'=\mathcal T_n$, also ist $\mathcal A\subseteq\mathcal T_n$.
Umgekehrt fixiere man $A\in M_d$ und wende Lemma B.5 auf die paarweise vertauschenden $x_m=A^{(m)}$ an: Dann ist
$e_n=A^{\otimes n}$ und $p_j=d\Pi(A^j)\in\mathcal A$, und per Induktion nach $k$ liegt jedes
$e_k=\frac1k\sum_j(-1)^{j-1}e_{k-j}p_j$ in $\mathcal A$. Also ist $A^{\otimes n}\in\mathcal A$ und $\mathcal T_n\subseteq\mathcal A$.

(c) Die Kommutante einer Menge ist gleich der Kommutante der von ihr erzeugten Algebra. Nach (a) und (b) erzeugen die
drei Mengen $\mathcal T_n=\pi(S_n)'$, deren Kommutante $\pi(S_n)''$ ist. Der Spann von $\pi(S_n)$ ist eine $*$-Algebra
($\pi(\sigma)\pi(\tau)=\pi(\sigma\tau)$, $\pi(\sigma)^*=\pi(\sigma^{-1})$, $\pi(e)=1$), also ist
$\pi(S_n)''=\operatorname{span}\pi(S_n)$ nach Satz 2.2. $\square$

**Korollar B.7 (Satz 4.14(b)).** Für den Gesamtspin $S_a=\sum_is_a^{(i)}$ von $L$ Spins $\frac12$ ist die von
$S_x,S_y,S_z$ erzeugte Algebra gleich $\{S_a\}''=\pi(S_L)'$, und $\{S_a\}'=\operatorname{span}\pi(S_L)$.

*Beweis.* $M_2=\operatorname{span}\{1,s_x,s_y,s_z\}$ (die Pauli-Matrizen und $1$ bilden eine Basis), $d\Pi(1)=L\cdot1$
und $d\Pi(s_a)=S_a$. Nach Satz B.6(b) ist die von $S_x,S_y,S_z$ erzeugte unitale Algebra $\mathcal A$ gleich
$\mathcal T_L=\pi(S_L)'$. Da die $S_a$ hermitesch sind, ist $\mathcal A$ eine $*$-Algebra (die Adjungierte eines Produkts
von Erzeugern ist das Produkt in umgekehrter Reihenfolge), also $\mathcal A=\mathcal A''\supseteq\{S_a\}''$ nach Satz 2.2,
und $\{S_a\}''\supseteq\mathcal A$, weil $\{S_a\}''$ eine unitale Algebra ist, die die $S_a$ enthält. Schließlich ist
$\{S_a\}'=\{S_a\}'''=\pi(S_L)''=\operatorname{span}\pi(S_L)$ (Abschnitt 2.1 und Satz B.6(c)). $\square$

Nach Satz B.6(a) für $d=2$ wird dieselbe Algebra von den Drehungen $V^{\otimes L}$, $V\in SU(2)$, aufgespannt; dafür
braucht man die in Abschnitt 4.5 erwähnte exponentielle Parametrisierung von $SU(2)$ nicht.

Der Code vergleicht für kleine $d$ und $n$ vier Dimensionen: die von $\operatorname{span}\pi(S_n)$, die der Kommutante
zweier zufälliger $V^{\otimes n}$ (die $\{V^{\otimes n}\}'=\operatorname{span}\pi(S_n)$ enthält und in allen ausgegebenen
Fällen damit übereinstimmt) sowie die der Spanne vieler zufälliger $A^{\otimes n}$ und $V^{\otimes n}$ mit
$\binom{d^2+n-1}n$. Außerdem prüft er, dass `permutation_representation` aus der Bibliothek der Konvention von
Lemma B.1(a) folgt, sowie Lemma B.2 und die Newton-Rekursion aus dem Beweis von Satz B.6(b).

```python
import functools
import itertools
import math
import torch
import torch_vn_algebra.decompose as dc

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    return (torch.as_tensor(X) - torch.as_tensor(Y)).abs().max().item() < tol

def haar(d, B):                                        # Theorem 5.3: QR of Ginibre matrices, phases fixed
    Q, R = torch.linalg.qr(torch.randn(B, d, d, dtype=dt))
    r = R.diagonal(dim1=-2, dim2=-1)
    return (Q * (r / r.abs()).unsqueeze(-2)).contiguous()

def power(A, n):                                       # A^{(x) n} for a batch A of shape (B, d, d)
    out = A
    for _ in range(n - 1):
        out = torch.einsum('bij,bkl->bikjl', out, A).reshape(A.shape[0], out.shape[1] * A.shape[1], -1)
    return out

def rank(X, tol=1e-9):                                 # dimension of the span of a batch of matrices
    return torch.linalg.matrix_rank(X.reshape(X.shape[0], -1), atol=tol).item()

def commutant_dim(ops, tol=1e-9):                      # dim {X : XA = AX for all A in ops}
    eye = torch.eye(ops.shape[-1], dtype=dt)
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    return int((torch.linalg.svdvals(L) < tol).sum())

for d, n in ((2, 2), (2, 3), (2, 4), (3, 2), (3, 3)):
    P = dc.permutation_representation(n, d)            # the n! operators pi(sigma)
    N_sym = math.comb(d * d + n - 1, n)                # Proposition B.3
    span_A = rank(power(torch.randn(N_sym + 5, d, d, dtype=dt), n))
    span_V = rank(power(haar(d, N_sym + 5), n))
    comm_V = commutant_dim(power(haar(d, 2), n))       # commutant of two random V^{(x) n}
    print(d, n, rank(P), comm_V, span_A, span_V, N_sym)
    assert span_A == span_V == N_sym == commutant_dim(P)                     # Theorem B.6(a)
    assert comm_V == rank(P)              # numerical observation: two V suffice here (B.6(c) gives >=)
    assert rank(P) == (math.factorial(n) if n <= d else math.comb(2 * n, n) // (n + 1))   # B.1(d), B.10

A = torch.randn(3, 2, 2, dtype=dt)                     # Lemma B.2 for n = 3
lhs = sum(functools.reduce(torch.kron, [A[s[0]], A[s[1]], A[s[2]]]) for s in itertools.permutations(range(3)))
rhs = sum((-1) ** (3 - sum(eps)) * functools.reduce(torch.kron, [sum(e * a for e, a in zip(eps, A))] * 3)
          for eps in itertools.product((0, 1), repeat=3))
assert close(lhs, rhs)

xs = torch.randn(3, 2, dtype=dt)                       # Lemma B.1(a): pi(s) x_1 (x) x_2 (x) x_3 = (x) x_{s^-1(p)}
for s in itertools.permutations(range(3)):
    inv = [s.index(m) for m in range(3)]
    Ps = dc.permutation_representation(3, 2, [list(s)])[0]
    assert close(Ps @ functools.reduce(torch.kron, list(xs)), functools.reduce(torch.kron, [xs[inv[p]] for p in range(3)]))

def d_pi(X, n):                                        # sum_m 1 (x) .. (x) X (x) .. (x) 1
    eye = torch.eye(X.shape[0], dtype=dt)
    return sum(functools.reduce(torch.kron, [X if l == m else eye for l in range(n)]) for m in range(n))

n, A = 4, torch.randn(2, 2, dtype=dt)
e = [torch.eye(2 ** n, dtype=dt)]
for k in range(1, n + 1):                              # Lemma B.5 with x_m = A^(m), p_j = dPi(A^j)
    e.append(sum((-1) ** (j - 1) * e[k - j] @ d_pi(torch.linalg.matrix_power(A, j), n)
                 for j in range(1, k + 1)) / k)
assert close(e[n], functools.reduce(torch.kron, [A] * n), 1e-8)
```

## B.4 Die Zerlegung von $(\mathbb C^d)^{\otimes n}$

**Satz B.8.** Es gibt $C\ge1$, ganze Zahlen $k_c,m_c\ge1$ mit $\sum_ck_cm_c=d^n$ und eine unitäre Abbildung
$W:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathcal H_n$ mit
$$
W^*\pi(\sigma)W=\bigoplus_{c=1}^C\rho_c(\sigma)\otimes1_{m_c},\qquad W^*V^{\otimes n}W=\bigoplus_{c=1}^C1_{k_c}\otimes\tau_c(V)
\qquad(\sigma\in S_n,\ V\in U(d)),
$$
wobei die $\rho_c$ paarweise inäquivalente irreduzible unitäre Darstellungen von $S_n$ auf $\mathbb C^{k_c}$ sind und
die $\tau_c$ paarweise inäquivalente irreduzible unitäre Darstellungen von $U(d)$ auf $\mathbb C^{m_c}$. Ferner ist
$W^*\operatorname{span}\pi(S_n)W=\bigoplus_cM_{k_c}\otimes1_{m_c}$, $W^*\mathcal U_nW=\bigoplus_c1_{k_c}\otimes M_{m_c}$,
$$
\sum_ck_c^2=\dim\operatorname{span}\pi(S_n),\qquad\sum_cm_c^2=\binom{d^2+n-1}n,
$$
und $\sum_ck_c^2=n!$, falls $n\le d$.

*Beweis.* Proposition 2.7 für die endliche Gruppe $S_n$ und Satz 2.6 für $\mathcal M=\operatorname{span}\pi(S_n)=\pi(S_n)''$
liefern $W$, die $\rho_c$ mit den angegebenen Eigenschaften und $W^*\pi(S_n)'W=\bigoplus_c1_{k_c}\otimes M_{m_c}$. Nach
Satz B.6(a) ist $\pi(S_n)'=\mathcal U_n$. Also ist $W^*V^{\otimes n}W=\bigoplus_c1\otimes\tau_c(V)$ mit gewissen
$\tau_c(V)\in M_{m_c}$, und $V\mapsto\tau_c(V)$ ist eine unitäre Darstellung, weil $V\mapsto V^{\otimes n}$ eine ist. Da die
$V^{\otimes n}$ den Raum $\mathcal U_n$ aufspannen, ist die lineare Abbildung
$\sum_la_lV_l^{\otimes n}\mapsto(\sum_la_l\tau_c(V_l))_c$ surjektiv auf $\bigoplus_cM_{m_c}$; Irreduzibilität und
Inäquivalenz der $\tau_c$ folgen wörtlich wie im Beweis von Proposition 2.7. Die Dimensionen sind Satz 2.6,
Proposition B.3 und Lemma B.1(d). $\square$

Für $d=2$ werden die Sektoren durch den Gesamtspin identifiziert:

**Korollar B.9 ($L$ Spins $\frac12$).** Sei $V$ unitär mit $V^*S_aV=\bigoplus_SS^{(S)}_a\otimes1_{m_S}$ für den Gesamtspin
von $L$ Spins $\frac12$ (Korollar 4.12, Proposition 4.13). Dann ist
$V^*\pi(\sigma)V=\bigoplus_S1_{2S+1}\otimes\rho_S(\sigma)$, wobei die $\rho_S$ paarweise inäquivalente irreduzible unitäre
Darstellungen von $S_L$ der Dimension
$$
m_S=\binom Lj-\binom L{j-1},\qquad j=\tfrac L2-S\in\{0,1,\dots,\lfloor L/2\rfloor\}.
$$
sind. Insbesondere ist $\dim\operatorname{span}\pi(S_L)=\sum_Sm_S^2$ und $\dim\pi(S_L)'=\sum_S(2S+1)^2=\binom{L+3}3$.

*Beweis.* Nach Satz 4.3 ist $V^*\{S_a\}'V=\bigoplus_S1_{2S+1}\otimes M_{m_S}$, und $\{S_a\}'=\operatorname{span}\pi(S_L)$
nach Korollar B.7. Also ist $V^*\pi(\sigma)V=\bigoplus_S1\otimes\rho_S(\sigma)$ mit einer unitären Darstellung $\rho_S$, und
die Abbildung $\sum_\sigma a_\sigma\pi(\sigma)\mapsto(\sum_\sigma a_\sigma\rho_S(\sigma))_S$ ist surjektiv auf $\bigoplus_SM_{m_S}$;
Irreduzibilität und Inäquivalenz folgen wie im Beweis von Proposition 2.7. Die Werte von $m_S$ sind Proposition 4.13.
Nach Satz 4.3 ist $\dim\{S_a\}'=\sum_Sm_S^2$ und $\dim\{S_a\}''=\sum_S(2S+1)^2$; nach Korollar B.7 ist
$\{S_a\}'=\operatorname{span}\pi(S_L)$ und $\{S_a\}''=\pi(S_L)'$, dessen Dimension nach Proposition B.3 mit $d=2$ gleich
$\binom{4+L-1}L=\binom{L+3}3$ ist. $\square$

**Korollar B.10 (Catalan-Zahlen).** Für $L\ge1$ ist $\displaystyle\sum_Sm_S^2=\binom{2L}L\frac1{L+1}$; das ist
$\dim\{S_a\}'=\dim\operatorname{span}\pi(S_L)$ auf $(\mathbb C^2)^{\otimes L}$.

*Beweis.* Sei $N=2^L$ und $\operatorname{vec}:M_N\to\mathbb C^N\otimes\mathbb C^N$ die lineare Bijektion
$|x\rangle\langle y|\mapsto x\otimes\bar y$. Dann ist $\operatorname{vec}(AX-XA)=(A\otimes1-1\otimes A^T)\operatorname{vec}X$, da
$A|x\rangle\langle y|B=|Ax\rangle\langle B^*y|\mapsto Ax\otimes B^T\bar y$. Also wird $\{S_a\}'$ auf den gemeinsamen Kern der
$S_a\otimes1-1\otimes S_a^T$, $a=x,y,z$, abgebildet, wobei $S_a^T=\sum_i(s_a^T)^{(i)}$. Aus $s_x^T=s_x$, $s_y^T=-s_y$,
$s_z^T=s_z$ und $\sigma_ys_x\sigma_y=-s_x$, $\sigma_ys_y\sigma_y=s_y$, $\sigma_ys_z\sigma_y=-s_z$ folgt
$\sigma_y(-s_a^T)\sigma_y=s_a$, also liefert die unitäre Matrix $Y=1\otimes\sigma_y^{\otimes L}$
$Y(S_a\otimes1-1\otimes S_a^T)Y^*=S_a\otimes1+1\otimes S_a$, den Gesamtspin von $2L$ Spins $\frac12$ auf
$(\mathbb C^2)^{\otimes2L}=\mathbb C^N\otimes\mathbb C^N$. In der Zerlegung aus Korollar 4.12 für $2L$ Spins wird ein Vektor,
der von allen drei Komponenten annulliert wird, von $\mathbf S^2$ annulliert, das auf den Summanden mit Spin
$S\ge\frac12$ gleich $S(S+1)\ne0$ ist; umgekehrt verschwinden die Komponenten auf dem Summanden mit Spin $0$. Der
gemeinsame Kern ist also dieser Summand, der Dimension
$m_0(2L)=\binom{2L}L-\binom{2L}{L-1}=\binom{2L}L\big(1-\frac L{L+1}\big)$ (Proposition 4.13). Die erste Gleichheit der
Dimensionen ist Satz 4.3, die zweite Korollar B.7. $\square$

Der Code zerlegt $\pi(S_3)$ auf $(\mathbb C^3)^{\otimes3}$ und prüft die Dimensionen aus Satz B.8; dann zerlegt er
$\pi(S_L)$ auf $(\mathbb C^2)^{\otimes L}$ für $L=5,6$ und vergleicht die Paare $(k_c,m_c)$ mit den Paaren $(m_S,2S+1)$
aus Korollar B.9.

```python
D33 = dc.from_group(dc.permutation_representation(3, 3))          # S_3 on (C^3)^{(x)3}
print(D33)
assert D33.verify()['ok']
assert sum(k * k for k in D33.k_factors) == 6 == math.factorial(3)            # Theorem B.8, n <= d
assert sum(m * m for m in D33.multiplicities) == math.comb(11, 3)              # 165

for L in (5, 6):
    DL = dc.from_group(dc.permutation_representation(L, 2))
    pairs = sorted((math.comb(L, j) - (math.comb(L, j - 1) if j else 0), L - 2 * j + 1)
                   for j in range(L // 2 + 1))                                 # (m_S, 2S + 1), S = L/2 - j
    print(L, DL.k_factors, DL.multiplicities)
    assert sorted(zip(DL.k_factors, DL.multiplicities)) == pairs               # Corollary B.9
    assert sum(k * k for k in DL.k_factors) == math.comb(2 * L, L) // (L + 1)  # Corollary B.10
    assert sum(m * m for m in DL.multiplicities) == math.comb(L + 3, 3)
```

## B.5 Der Weingarten-Kalkül

**Lemma B.11 (der Haar-Twirl ist die Orthogonalprojektion auf die Kommutante).** Sei $U$ Haar-verteilt auf $U(d)$ und
$\Phi(X)=\mathbb E[U^{\otimes n}X(U^{\otimes n})^*]$ für $X\in\operatorname{End}\mathcal H_n$. Dann ist $\Phi$ die
Orthogonalprojektion von $\operatorname{End}\mathcal H_n$ auf $\{V^{\otimes n}:V\in U(d)\}'=\operatorname{span}\pi(S_n)$.

*Beweis.* Die Einträge von $U^{\otimes n}X(U^{\otimes n})^*$ sind Polynome in den Einträgen von $U$ und $\bar U$, die
durch $1$ beschränkt sind, also ist $\Phi$ eine wohldefinierte lineare Abbildung. Für $V\in U(d)$ hat $VU$ die
Verteilung von $U$, also ist $V^{\otimes n}\Phi(X)(V^{\otimes n})^*=\mathbb E[(VU)^{\otimes n}X((VU)^{\otimes n})^*]=\Phi(X)$:
Das Bild liegt in der Kommutante. Liegt $Y$ in der Kommutante, so ist $U^{\otimes n}Y(U^{\otimes n})^*=Y$ und
$\Phi(Y)=Y$. Also ist $\Phi$ idempotent mit der Kommutante als Bild. Nach Lemma 5.1 hat $U^*=U^{-1}$ die Verteilung von
$U$, daher ist
$\langle Y,\Phi(X)\rangle=\mathbb E\operatorname{Tr}\big(((U^*)^{\otimes n}YU^{\otimes n})^*X\big)=\langle\Phi(Y),X\rangle$: $\Phi$ ist
selbstadjungiert. Eine selbstadjungierte idempotente Abbildung ist die Orthogonalprojektion auf ihr Bild
($\langle X-\Phi(X),\Phi(Z)\rangle=\langle\Phi(X-\Phi(X)),Z\rangle=0$). Die Kommutante ist $\operatorname{span}\pi(S_n)$ nach
Satz B.6(c). $\square$

Um diese Projektion mit Hilfe der aufspannenden Familie $\pi(\sigma)$ zu schreiben, die für $n>d$ linear abhängig ist,
verwenden wir die **Pseudoinverse** einer hermiteschen Matrix $G$: $G^+=f(G)$ mit $f(\lambda)=1/\lambda$ für
$\lambda\neq0$ und $f(0)=0$ (Funktionalkalkül, Satz 1.5). Sie ist hermitesch, $G^+=G^{-1}$, falls $G$ invertierbar ist,
$GG^+=\mathbf 1_{\neq0}(G)$ ist die Projektion auf $(\ker G)^\perp$, $G^+v=v/\lambda$, falls $Gv=\lambda v$ mit
$\lambda\neq0$, und $(R^*GR)^+=R^*G^+R$ für unitäres $R$ (ist $G=Q\Lambda Q^*$, so ist $R^*GR=(R^*Q)\Lambda(R^*Q)^*$).

**Lemma B.12.** Seien $v_1,\dots,v_N$ Vektoren in einem endlichdimensionalen Innenproduktraum $\mathcal K$,
$G_{ab}=\langle v_a,v_b\rangle$ ihre Gram-Matrix und $W=G^+$. Dann ist die Orthogonalprojektion auf
$\operatorname{span}\{v_a\}$ gleich $x\mapsto\sum_{a,b}v_a\,W_{ab}\,\langle v_b,x\rangle$.

*Beweis.* Sei $T:\mathbb C^N\to\mathcal K$, $Tc=\sum_ac_av_a$. Dann ist $T^*x=(\langle v_b,x\rangle)_b$, $G=T^*T$, und die
fragliche Abbildung ist $Q=TG^+T^*$. Wegen $\langle c,Gc\rangle=\|Tc\|^2$ ist $\ker G=\ker T$, also $TG^+G=T(1-\Pi)=T$
für die Projektion $\Pi$ auf $\ker T$. Daher ist $QT=T$: $Q$ ist die Identität auf $\operatorname{ran}T$. $Q$ verschwindet
auf $(\operatorname{ran}T)^\perp=\ker T^*$, und $Q^*=Q$. Also ist $Q$ die Orthogonalprojektion auf $\operatorname{ran}T$. $\square$

**Satz B.13 (Weingarten-Formel beliebiger Ordnung).** Sei $n,d\ge1$, $U$ Haar-verteilt auf $U(d)$, $G$ die
$n!\times n!$-Matrix $G_{\sigma\tau}=d^{\#\text{cycles}(\sigma^{-1}\tau)}$ ($\sigma,\tau\in S_n$) und $W=G^+$.
(a) $W_{\sigma\tau}=\mathrm{Wg}(\sigma^{-1}\tau)$ für die Funktion $\mathrm{Wg}(\rho)=\mathrm{Wg}(\rho,d)=W_{e\rho}$ auf $S_n$, die
nur vom Zyklentyp von $\rho$ abhängt; insbesondere ist $\mathrm{Wg}(\rho^{-1})=\mathrm{Wg}(\rho)$.
(b) $G$ ist genau dann invertierbar, wenn $n\le d$, und dann ist $W=G^{-1}$.
(c) Für $X\in\operatorname{End}\mathcal H_n$ ist
$\mathbb E[U^{\otimes n}X(U^{\otimes n})^*]=\sum_{\sigma,\tau}\mathrm{Wg}(\sigma^{-1}\tau)\operatorname{Tr}(\pi(\tau)^*X)\,\pi(\sigma)$.
(d) Für alle $\mathbf i,\mathbf j,\mathbf i',\mathbf j'\in[d]^n$ gilt
$$
\mathbb E\big[U_{i_1j_1}\cdots U_{i_nj_n}\overline{U_{i'_1j'_1}\cdots U_{i'_nj'_n}}\big]
=\sum_{\sigma,\tau\in S_n}\prod_{m=1}^n\delta_{i_mi'_{\sigma(m)}}\delta_{j_mj'_{\tau(m)}}\;\mathrm{Wg}(\sigma\tau^{-1}),
$$
d. h. die Summe läuft über die Paare mit $\mathbf i=\mathbf i'\circ\sigma$ und $\mathbf j=\mathbf j'\circ\tau$.
(e) Momente mit verschiedenen Anzahlen $p\neq q$ von Faktoren $U$ und $\bar U$ verschwinden.

Das ist Satz 5.12 (mit $k=d$), nun bewiesen.

*Beweis.* (a) Man schreibe $g(\rho)=d^{\#\text{cycles}(\rho)}$, sodass $G_{\sigma\tau}=g(\sigma^{-1}\tau)$, und seien
$\Lambda_\kappa,R_\kappa$ die Permutationsmatrizen auf $\mathbb C^{S_n}$ mit $\Lambda_\kappa e_\sigma=e_{\kappa\sigma}$,
$R_\kappa e_\sigma=e_{\sigma\kappa}$. Dann ist $(\Lambda_\kappa^*G\Lambda_\kappa)_{\sigma\tau}=g(\sigma^{-1}\kappa^{-1}\kappa\tau)=G_{\sigma\tau}$
und $(R_\kappa^*GR_\kappa)_{\sigma\tau}=g(\kappa^{-1}\sigma^{-1}\tau\kappa)=G_{\sigma\tau}$ (Lemma B.1(c)). Nach den Eigenschaften der
Pseudoinversen hat $W$ dieselben beiden Invarianzen: $W_{\kappa\sigma,\kappa\tau}=W_{\sigma\tau}$ und
$W_{\sigma\kappa,\tau\kappa}=W_{\sigma\tau}$. Die erste mit $\kappa=\sigma^{-1}$ liefert $W_{\sigma\tau}=\mathrm{Wg}(\sigma^{-1}\tau)$;
die zweite liefert $\mathrm{Wg}(\kappa^{-1}\rho\kappa)=\mathrm{Wg}(\rho)$, also ist $\mathrm{Wg}$ konstant auf den
Konjugationsklassen, und diese sind die Zyklentypen (Lemma B.1(c)); $\rho^{-1}$ hat den Zyklentyp von $\rho$.
(b) $G$ ist die Gram-Matrix der $\pi(\sigma)$ (Lemma B.1(c)), und mit $T$ wie im Beweis von Lemma B.12 ist $G=T^*T$
genau dann invertierbar, wenn $T$ injektiv ist, d. h. wenn die $\pi(\sigma)$ linear unabhängig sind: nach Lemma B.1(d)
genau dann, wenn $n\le d$.
(c) Lemma B.11 und Lemma B.12 mit $\mathcal K=\operatorname{End}\mathcal H_n$ und $v_\sigma=\pi(\sigma)$.
(d) Die Einträge eines Kronecker-Produkts sind $\langle e_{\mathbf i},U^{\otimes n}e_{\mathbf j}\rangle=\prod_mU_{i_mj_m}$, also ist
das Moment gleich $\langle e_{\mathbf i},\Phi(X)e_{\mathbf i'}\rangle$ für $X=|e_{\mathbf j}\rangle\langle e_{\mathbf j'}|$. Dabei ist
$\operatorname{Tr}(\pi(\tau)^*X)=\langle\pi(\tau)e_{\mathbf j'},e_{\mathbf j}\rangle=\delta_{\mathbf j,\mathbf j'\circ\tau^{-1}}$ und
$\langle e_{\mathbf i},\pi(\sigma)e_{\mathbf i'}\rangle=\delta_{\mathbf i,\mathbf i'\circ\sigma^{-1}}$. Nach (c) ist das Moment
$\sum_{\sigma,\tau}\delta_{\mathbf i,\mathbf i'\circ\sigma^{-1}}\delta_{\mathbf j,\mathbf j'\circ\tau^{-1}}\mathrm{Wg}(\sigma^{-1}\tau)$; die
Substitution $\sigma\to\sigma^{-1}$, $\tau\to\tau^{-1}$ liefert die Formel.
(e) Für reelles $\theta$ ist $e^{i\theta}1\in U(d)$, also hat $e^{i\theta}U$ die Verteilung von $U$, und das Moment ist
gleich $e^{i(p-q)\theta}$ mal sich selbst. Man wähle $\theta$ mit $e^{i(p-q)\theta}\neq1$. $\square$

Für $n=2\le d$ ist $G=\begin{pmatrix}d^2&d\\d&d^2\end{pmatrix}$, und $G^{-1}$ liefert $\mathrm{Wg}(e)=\frac1{d^2-1}$,
$\mathrm{Wg}(\tau)=-\frac1{d(d^2-1)}$: Satz 5.10. Für $d=1$ ist jedes $G_{\sigma\tau}=1$, $W=\frac1{(n!)^2}G$, und das Moment
der Phase $U=U_{11}$ ist $\sum_{\sigma,\tau}\frac1{(n!)^2}=1=|U_{11}|^{2n}$, wie es sein muss.

**Korollar B.14.** Für alle $n,d\ge1$ ist $\mathbb E|U_{11}|^{2n}=\dfrac{n!}{d(d+1)\cdots(d+n-1)}$.

*Beweis.* Mit allen Indizes gleich $1$ in Satz B.13(d) ist das Moment $\sum_{\sigma,\tau}\mathrm{Wg}(\sigma\tau^{-1})=\sum_{\sigma,\tau}W_{\sigma^{-1},\tau^{-1}}=\mathbf 1^TW\mathbf 1$
(Teil (a)). Jede Zeilensumme von $G$ ist $c_n(d)=\sum_{\rho\in S_n}d^{\#\text{cycles}(\rho)}\ge d^n>0$, also ist
$G\mathbf 1=c_n(d)\mathbf 1$, $W\mathbf 1=\mathbf 1/c_n(d)$, und das Moment ist $n!/c_n(d)$. Schließlich ist
$c_n(x)=x(x+1)\cdots(x+n-1)$ (Kapitel 5, Aufgabe 2, Lösung). $\square$

Das erweitert Kapitel 5, Aufgabe 2, auf $n>d$, wo $G$ singulär ist. Der Code berechnet $\mathrm{Wg}$ aus der
Pseudoinversen, prüft (a), (b) und Satz 5.10 und vergleicht den gesamten Momententensor der Ordnung $n=3$ für $d=2$ –
einen Fall mit singulärem $G$ – sowie den Twirl eines zufälligen $X$ mit Monte-Carlo-Mitteln über Haar-verteilte
unitäre Matrizen.

```python
def n_cycles(p):                                        # number of cycles of a permutation in one-line form
    seen, count = set(), 0
    for i in range(len(p)):
        count += i not in seen
        while i not in seen:
            seen.add(i)
            i = p[i]
    return count

def compose(s, t):                                     # (s t)(m) = s(t(m))
    return tuple(s[m] for m in t)

def inverse(s):
    return tuple(sorted(range(len(s)), key=lambda m: s[m]))

def weingarten(n, d):
    perms = list(itertools.permutations(range(n)))
    G = torch.tensor([[float(d) ** n_cycles(compose(inverse(s), t)) for t in perms] for s in perms],
                     dtype=torch.float64)
    return perms, G, torch.linalg.pinv(G, hermitian=True)

for n, d in ((2, 2), (3, 2), (3, 3), (4, 2), (4, 3)):
    perms, G, W = weingarten(n, d)
    pos = {p: a for a, p in enumerate(perms)}
    Wg = {p: W[0, pos[p]].item() for p in perms}        # perms[0] is the identity
    assert all(abs(W[a, b] - Wg[compose(inverse(s), t)]) < 1e-12
               for a, s in enumerate(perms) for b, t in enumerate(perms))       # Theorem B.13(a)
    assert all(abs(Wg[p] - Wg[compose(compose(inverse(r), p), r)]) < 1e-12 for p in perms for r in perms)
    assert (torch.linalg.matrix_rank(G).item() == math.factorial(n)) == (n <= d)  # Theorem B.13(b)
for d in (2, 3, 5):                                     # Theorem 5.10
    _, _, W = weingarten(2, d)
    assert close(W[0], torch.tensor([1 / (d * d - 1), -1 / (d * (d * d - 1))], dtype=torch.float64), 1e-12)

n, d, B = 3, 2, 100_000
perms, G, W = weingarten(n, d)
P = dc.permutation_representation(n, d, [list(p) for p in perms])        # pi(sigma), same order as W
Phi = lambda X: torch.einsum('st,tij,ij,skl->kl', W.to(dt), P.conj(), X, P)   # Theorem B.13(c)
U3 = power(haar(d, B), n)                                                 # B samples of U^{(x)3}
moments = torch.einsum('bx,by->xy', U3.reshape(B, -1), U3.reshape(B, -1).conj()) / B
second = torch.einsum('bx,by->xy', U3.reshape(B, -1).abs() ** 2, U3.reshape(B, -1).abs() ** 2) / B
err = (second - moments.abs() ** 2).clamp(min=0).sqrt() / math.sqrt(B)   # standard errors
N = d ** n
exact = torch.stack([Phi(torch.outer(torch.eye(N, dtype=dt)[j], torch.eye(N, dtype=dt)[jp]))
                     for j in range(N) for jp in range(N)])               # [(j, j'), i, i']
exact = exact.reshape(N, N, N, N).permute(2, 0, 3, 1).reshape(N * N, N * N)   # [(i, j), (i', j')]
z = ((moments - exact).abs() / (err + 1e-15)).max().item()
print(f"{N ** 4} moments of order 3 for d = 2: largest deviation {z:.2f} standard errors")
assert ((moments - exact).abs() <= 5 * err + 1e-12).all()                 # Theorem B.13(d), singular G
assert abs(exact[0, 0].item() - 1 / 4) < 1e-12                            # Corollary B.14: 3!/(2*3*4)

X = torch.randn(N, N, dtype=dt)
T = Phi(X)
coef = torch.linalg.lstsq(P.reshape(len(perms), -1).T, T.reshape(-1, 1)).solution
assert close(Phi(T), T) and close(torch.einsum('s,sij->ij', coef[:, 0], P), T)   # a projection into span pi(S_n)
Vn = power(haar(d, 1), n)[0]
assert close(Vn @ T, T @ Vn)                                             # = {V^{(x)n}}'
assert close(Phi(P[3]), P[3])
Y = U3 @ X @ U3.mH                                                        # samples of U X U*
assert ((Y.mean(0) - T).abs() <= 5 * Y.std(0) / math.sqrt(B) + 1e-12).all()   # Lemma B.11
```

## B.6 Tensoroperatoren und das Wigner-Eckart-Theorem für $SU(2)$

Abschnitt 9.3 baut Spinkorrelationen aus dem Wigner-Eckart-Theorem auf. Wir beweisen es für die Spindarstellungen aus
Abschnitt 4.5. Spin $k$ ist die Darstellung $S^{(k)}_a$ aus Satz 4.11 auf $\mathbb C^{2k+1}$ mit der Basis $|k,q\rangle$,
$q=k,\dots,-k$; zur Erinnerung: $S^{(k)}_\pm|k,q\rangle=\sqrt{k(k+1)-q(q\pm1)}\,|k,q\pm1\rangle$ (für $S_-=S_+^*$ folgt das
aus der Formel für $S_+$).

**Definition.** Sei $(S_a)$ eine Spindarstellung auf $\mathbb C^N$ und $k\in\{0,\frac12,1,\dots\}$. Ein **Tensoroperator
vom Rang $k$** ist eine Familie $(T_q)_{q=-k,\dots,k}$ von $N\times N$-Matrizen mit
$$
[S_z,T_q]=q\,T_q,\qquad[S_\pm,T_q]=\sqrt{k(k+1)-q(q\pm1)}\;T_{q\pm1}
$$
(der Koeffizient verschwindet für $q=\pm k$, wo $T_{\pm(k+1)}$ nicht definiert ist). Ein **Vektoroperator** ist ein
Tripel $(T_x,T_y,T_z)$ mit $[S_a,T_b]=i\sum_c\varepsilon_{abc}T_c$, d. h. $[S_x,T_y]=iT_z$, $[S_x,T_z]=-iT_y$ und zyklisch,
sowie $[S_a,T_a]=0$.

**Lemma B.15.** (a) Ist $(T_a)$ ein Vektoroperator, so bilden seine **sphärischen Komponenten**
$T_{\pm1}=\mp\frac1{\sqrt2}(T_x\pm iT_y)$, $T_0=T_z$ einen Tensoroperator vom Rang $1$, und für zwei Vektoroperatoren gilt
$\sum_aT_aT'_a=\sum_q(-1)^qT_qT'_{-q}$.
(b) Für $L$ Spins $\frac12$ und den Gesamtspin $S_a$ sind jeder Platzspin $(s^{(i)}_a)_a$ und $(S_a)_a$ selbst
Vektoroperatoren.
(c) $(T_q)$ ist genau dann ein Tensoroperator vom Rang $k$, wenn die Abbildung
$\Psi:\mathbb C^N\otimes\mathbb C^{2k+1}\to\mathbb C^N$, $\Psi(v\otimes|k,q\rangle)=T_qv$, die Gleichung
$\Psi\,(S_a\otimes1+1\otimes S^{(k)}_a)=S_a\Psi$ für $a=x,y,z$ erfüllt.

*Beweis.* (a) Aus den Relationen folgen $[S_z,T_x\pm iT_y]=\pm(T_x\pm iT_y)$, $[S_\pm,T_z]=\mp(T_x\pm iT_y)$,
$[S_\pm,T_x\mp iT_y]=\pm2T_z$ und $[S_\pm,T_x\pm iT_y]=0$ (z. B. $[S_x+iS_y,T_x-iT_y]=-i[S_x,T_y]+i[S_y,T_x]=T_z+T_z$). In
sphärischen Komponenten: $[S_z,T_q]=qT_q$, $[S_\pm,T_0]=\sqrt2\,T_{\pm1}$, $[S_\pm,T_{\mp1}]=\sqrt2\,T_0$,
$[S_\pm,T_{\pm1}]=0$; das sind die definierenden Relationen für $k=1$ ($\sqrt{2-q(q\pm1)}$ ist $\sqrt2$ für $q=0$ und
$q=\mp1$ und $0$ für $q=\pm1$). Für das Produkt gilt
$T_{+1}T'_{-1}+T_{-1}T'_{+1}=-\frac12\big((T_x+iT_y)(T'_x-iT'_y)+(T_x-iT_y)(T'_x+iT'_y)\big)=-(T_xT'_x+T_yT'_y)$.
(b) $[S_a,s^{(i)}_b]=[s^{(i)}_a,s^{(i)}_b]$ (andere Plätze vertauschen), und die Spin-$\frac12$-Matrizen erfüllen die
Spinrelationen; für $S_a$ summiere man über $i$. (c)
$\Psi(S_a\otimes1+1\otimes S^{(k)}_a)(v\otimes|k,q\rangle)=T_qS_av+\sum_{q'}(S^{(k)}_a)_{q'q}T_{q'}v$, also besagt die
Verflechtungsrelation für $S_a$, dass $[S_a,T_q]=\sum_{q'}(S^{(k)}_a)_{q'q}T_{q'}$ für alle $q$. Für $S_z,S_\pm$ sind das die
definierenden Relationen, und die Relationen für $S_x,S_y,S_z$ und für $S_z,S_\pm$ sind Linearkombinationen
voneinander. $\square$

**Satz B.16 (Wigner-Eckart).** Sei $(S_a)$ eine Spindarstellung auf $\mathbb C^N$, $V$ unitär mit
$V^*S_aV=\bigoplus_SS^{(S)}_a\otimes1_{m_S}$ (Korollar 4.12), $|S,M;b\rangle=V(|S,M\rangle\otimes f_b)$ die zugehörige
Orthonormalbasis ($f_b$ die Standardbasis von $\mathbb C^{m_S}$) und $(T_q)$ ein Tensoroperator vom Rang $k$. Die Spins
$S,S'$ mögen vorkommen ($m_S,m_{S'}\ge1$).
(a) Ist $S'\notin\{|S-k|,|S-k|+1,\dots,S+k\}$, so ist $\langle S',M';a|T_q|S,M;b\rangle=0$ für alle $M',M,q,a,b$.
(b) Andernfalls bilden die Verflechtungsoperatoren von Spin $S\otimes$ Spin $k$ nach Spin $S'$ einen eindimensionalen
Raum; für jedes solche $C=C^{S'}_{S,k}:\mathbb C^{2S+1}\otimes\mathbb C^{2k+1}\to\mathbb C^{2S'+1}$ mit $C\neq0$ gibt es eine
eindeutige $m_{S'}\times m_S$-Matrix $R=\langle S'\|T\|S\rangle$, das **reduzierte Matrixelement**, mit
$$
\langle S',M';a|\,T_q\,|S,M;b\rangle=\langle S',M'|\,C\,(|S,M\rangle\otimes|k,q\rangle)\rangle\;R_{ab}\qquad\text{for all }M',M,q,a,b .
$$
(c) Das Element verschwindet, außer wenn $M'=M+q$. Im Fall (b) gibt es ein $q$ mit
$\langle S',S'|C(|S,S'-q\rangle\otimes|k,q\rangle)\rangle\neq0$, sodass $R$ und damit der ganze Block von $T$ zwischen den
Multipletts $S$ und $S'$ durch die Matrixelemente $\langle S',S';a|T_q|S,S'-q;b\rangle$ zwischen den Vektoren höchsten
Gewichts von $S'$ und je einem Vektor jeder Kopie von $S$ bestimmt ist.

*Beweis.* Man fixiere $a\le m_{S'}$, $b\le m_S$, und sei $J:\mathbb C^{2S+1}\to\mathbb C^N$, $Jx=V(x\otimes f_b)$ (im
Summanden $S$), und $K:\mathbb C^N\to\mathbb C^{2S'+1}$ die Adjungierte von $x\mapsto V(x\otimes f_a)$ (im Summanden $S'$).
Wir schreiben $\alpha\in\{x,y,z\}$ für den Komponentenindex. Aus der Gestalt von $V^*S_\alpha V$ folgt
$S_\alpha J=JS^{(S)}_\alpha$ und, durch Adjungieren derselben Relation für $S'$ (die Spinmatrizen sind hermitesch),
$KS_\alpha=S^{(S')}_\alpha K$. Sei $t=K\Psi(J\otimes1)$ mit $\Psi$ aus Lemma B.15(c); seine Matrixelemente sind
$\langle S',M'|t(|S,M\rangle\otimes|k,q\rangle)\rangle=\langle S',M';a|T_q|S,M;b\rangle$, und
$$
t\,(S^{(S)}_\alpha\otimes1+1\otimes S^{(k)}_\alpha)=K\Psi(S_\alpha\otimes1+1\otimes S^{(k)}_\alpha)(J\otimes1)=KS_\alpha\Psi(J\otimes1)=S^{(S')}_\alpha\,t,
$$
also ist $t$ ein Verflechtungsoperator von $\rho=$ Spin $S\otimes$ Spin $k$ (der Spindarstellung
$S^{(S)}_\alpha\otimes1+1\otimes S^{(k)}_\alpha$, Abschnitt 4.5) nach $\sigma=$ Spin $S'$. Adjungieren ist eine lineare
Bijektion zwischen $\operatorname{Hom}(\rho,\sigma)$ und $\operatorname{Hom}(\sigma,\rho)$ (die Adjungierte von
$t\rho(\alpha)=\sigma(\alpha)t$ ist $\rho(\alpha)t^*=t^*\sigma(\alpha)$, da alle beteiligten Matrizen hermitesch sind). Nach
Korollar 4.12 und Kapitel 4, Aufgabe 2 (in deren Lösung allein aus Korollar 4.12 bewiesen), ist $\rho$ unitär
äquivalent zu $\bigoplus_j(\text{Spin }j)$, mit je einer Kopie jedes $j\in\{|S-k|,|S-k|+1,\dots,S+k\}$. Diese Spins sind
irreduzibel und paarweise inäquivalent (Satz 4.11), also ist das eine Zerlegung wie in Satz 4.3, und deren letzte
Aussage liefert $\dim\operatorname{Hom}(\sigma,\rho)=1$, falls $S'$ eines dieser $j$ ist, und sonst
$\operatorname{Hom}(\sigma,\rho)=0$ (Spin $S'$ ist dann nach Satz 4.11 zu keinem von ihnen äquivalent). Das liefert
(a), und im Fall (b) ist $t=R_{ab}C$ mit einer eindeutigen Zahl $R_{ab}$.
(c) Die Verflechtungsrelation für $S_z$ liefert $M'\langle S',M'|t(|S,M\rangle\otimes|k,q\rangle)\rangle=(M+q)\langle S',M'|t(|S,M\rangle\otimes|k,q\rangle)\rangle$.
Das Bild von $C\neq0$ ist ein von null verschiedener invarianter Unterraum des irreduziblen Spins $S'$, also der ganze
Raum, sodass $\langle S',S'|C(|S,M\rangle\otimes|k,q\rangle)\rangle\neq0$ für einen Basisvektor gilt, und dann ist $M=S'-q$
nach derselben Relation für den Verflechtungsoperator $C$. Man dividiere durch diese Zahl. $\square$

Nach Lemma B.15(a) ist der Block zwischen den Multipletts $S$ und $S'$ eines Skalarprodukts wie
$\mathbf s^{(i)}\cdot\mathbf s^{(j)}=\sum_q(-1)^q s^{(i)}_qs^{(j)}_{-q}$ gleich
$\sum_{S''}\sum_q(-1)^q(\text{Block }S''\to S'\text{ von }s^{(i)}_q)(\text{Block }S\to S''\text{ von }s^{(j)}_{-q})$, und jeder
Faktor ist nach Satz B.16 durch ein reduziertes Matrixelement und eine Clebsch-Gordan-Abbildung gegeben. Das ist das
Prinzip hinter `su2.spin_correlation`; die Korrektheit seiner Implementierung prüft Kapitel 9 numerisch gegen die
Methode `'direct'`. Der Code prüft Lemma B.15 für die Platzspins von sechs Spins $\frac12$ und Satz B.16 in einer
basisfreien Form: In der von der Bibliothek berechneten Zerlegung (deren Sektorbasen nicht die $|S,M\rangle$ sein
müssen) hat der Block von $(T_q)_q$ zwischen zwei Multipletts, als Matrix mit Zeilen $(M',M,q)$ und Spalten $(a,b)$,
höchstens Rang eins, und er ist außerhalb des Bereichs aus (a) null.

```python
L = 6
S6 = dc.total_spin(L)
Sp, Sm = S6[0] + 1j * S6[1], S6[0] - 1j * S6[1]
D6 = dc.from_generators(S6)
print(D6)
off = [0]
for k, m in zip(D6.k_factors, D6.multiplicities):
    off.append(off[-1] + k * m)
sx, sy, sz = dc.spin_matrices(0.5)
site = lambda x, i: functools.reduce(torch.kron, [x if l == i else torch.eye(2, dtype=dt) for l in range(L)])
spherical = lambda X, Y, Z: {1: -(X + 1j * Y) / math.sqrt(2), 0: Z, -1: (X - 1j * Y) / math.sqrt(2)}
for i in (0, 3):
    T = spherical(site(sx, i), site(sy, i), site(sz, i))
    zero = torch.zeros_like(T[0])
    for q in (-1, 0, 1):                                                   # Lemma B.15(a), (b)
        assert close(S6[2] @ T[q] - T[q] @ S6[2], q * T[q])
        assert close(Sp @ T[q] - T[q] @ Sp, math.sqrt(2 - q * (q + 1)) * T.get(q + 1, zero))
        assert close(Sm @ T[q] - T[q] @ Sm, math.sqrt(2 - q * (q - 1)) * T.get(q - 1, zero))
    for c, (kc, mc) in enumerate(zip(D6.k_factors, D6.multiplicities)):
        for e, (ke, me) in enumerate(zip(D6.k_factors, D6.multiplicities)):
            Vc, Ve = D6.V[:, off[c]:off[c + 1]], D6.V[:, off[e]:off[e + 1]]
            block = torch.stack([(Ve.mH @ T[q] @ Vc).reshape(ke, me, kc, mc) for q in (-1, 0, 1)])
            sv = torch.linalg.svdvals(block.permute(0, 1, 3, 2, 4).reshape(3 * ke * kc, me * mc))
            S, S2 = (kc - 1) / 2, (ke - 1) / 2
            if abs(S - S2) <= 1 and S + S2 >= 1:
                assert sv.numel() == 1 or sv[1] < 1e-10                        # Theorem B.16(b): rank one
            else:
                assert sv[0] < 1e-10                                           # Theorem B.16(a)
X, Y, Z = (site(s, 0) for s in (sx, sy, sz))
X2, Y2, Z2 = (site(s, 1) for s in (sx, sy, sz))
T, T2 = spherical(X, Y, Z), spherical(X2, Y2, Z2)
assert close(X @ X2 + Y @ Y2 + Z @ Z2, sum((-1) ** q * T[q] @ T2[-q] for q in (-1, 0, 1)))   # Lemma B.15(a)
```

## B.7 Im Kurs verwendete Haar-Maße

Der Kurs verwendet die Haar-Maße der folgenden Gruppen: $U(k)$ (Satz 5.3), $U(1)$ (der Fall $k=1$; Kapitel 5,
Aufgabe 3), $U(\mathcal M)\cong\prod_cU(k_c)$ (Korollar 5.4), $SU(k)$ (Abschnitt 5.1, nach Korollar 5.4), die orthogonale
Gruppe $O(k)$ (Abschnitt 5.1, Kapitel 11, Aufgabe 2, und reelle Algebren) sowie die Diagonalgruppen $T$ aus
Lemma C.2. Die Eindeutigkeit ist für alle Lemma 5.1, da jede eine abgeschlossene Untergruppe eines $U(n)$ ist; die
Existenz ist in jedem Fall eine explizite Konstruktion (für $T$ in Lemma C.2), und die Konstruktion für $O(k)$ ist die
folgende. Ein allgemeiner Existenzsatz für das Haar-Maß einer kompakten Gruppe wird nirgends im Kurs verwendet.

**Satz B.17 (Haar-Maß der orthogonalen Gruppe).** Sei $Z$ eine reelle $k\times k$-Matrix mit unabhängigen
$N(0,1)$-Einträgen. Dann ist $Z$ mit Wahrscheinlichkeit eins invertierbar, $Z=OR$ mit $O$ orthogonal und $R$ reell,
oberer Dreiecksmatrix mit positiver Diagonale, und $O$ ist Haar-verteilt auf $O(k)$. Ist $Z=O'R'$ irgendeine reelle
QR-Zerlegung ($O'$ orthogonal, $R'$ reelle obere Dreiecksmatrix), so ist $O'\operatorname{diag}(\operatorname{sgn}R'_{ii})=O$.
Für Haar-verteiltes $O$ und $X\in M_k(\mathbb C)$ ist $\mathbb E[O]=0$ und $\mathbb E[OXO^T]=\frac{\operatorname{Tr}X}k1$;
insbesondere ist $\mathbb E[O_{ij}O_{lm}]=\delta_{il}\delta_{jm}/k$.

*Beweis.* $\det Z$ ist ein Polynom in den $k^2$ reellen Einträgen, nicht identisch null ($\det1=1$), und $Z$ hat eine
Dichte, also ist $\det Z\neq0$ fast sicher (Korollar A.7(b), wie im Beweis von Satz 5.3). Gram-Schmidt, angewandt auf
reelle Spalten, liefert reelle $q_j$ und reelle Koeffizienten, d. h. eine reelle Zerlegung $Z=OR$ mit $R_{jj}>0$;
wegen der Eindeutigkeit in Lemma 5.2 ist sie die Zerlegung aus Lemma 5.2, also hängt $O$ stetig von $Z$ ab. Die
Dichte $(2\pi)^{-k^2/2}e^{-\|Z\|_2^2/2}$ von $Z$ ist invariant unter $Z\mapsto VZ$ für $V\in O(k)$, eine lineare
Abbildung von $\mathbb R^{k^2}$, die $\|Z\|_2^2=\operatorname{Tr}(Z^TZ)$ erhält und daher orthogonal mit Jacobi-Determinante
$1$ ist; also hat $VZ$ die Verteilung von $Z$ (Proposition A.14(a)). Da $VZ=(VO)R$, liefert die Eindeutigkeit
$O(VZ)=VO(Z)$, also ist die Verteilung von $O$ linksinvariant auf der kompakten Gruppe $O(k)$, und Lemma 5.1 greift.
Die Vorzeichenkorrektur ist das Ende des Beweises von Satz 5.3, mit $\operatorname{sgn}R'_{ii}=R'_{ii}/|R'_{ii}|$
($R'_{ii}\neq0$, da $Z$ invertierbar ist). Da $-1\in O(k)$, hat $-O$ die Verteilung von $O$, also ist
$\mathbb E[O]=-\mathbb E[O]=0$ (die Einträge sind durch $1$ beschränkt). Schließlich erfüllt $T=\mathbb E[OXO^T]$ die
Gleichung $VTV^T=T$ für $V\in O(k)$, da $VO$ die Verteilung von $O$ hat. Mit $V$ der Diagonalmatrix mit $-1$ an
Position $i$ und $1$ sonst folgt $T_{ij}=-T_{ij}$ für $j\neq i$; mit Permutationsmatrizen sind alle $T_{ii}$ gleich
(wie in Lemma 5.6). Also ist $T=\lambda1$ und $k\lambda=\operatorname{Tr}T=\operatorname{Tr}X$. Für $X=|e_j\rangle\langle e_m|$
ist $(OXO^T)_{il}=O_{ij}O_{lm}$. $\square$

Der Code prüft per Monte-Carlo die zweiten Momente, $\mathbb E\operatorname{Tr}O=0$ und
$\mathbb E(\operatorname{Tr}O)^2=\sum_{i,l}\mathbb E[O_{ii}O_{ll}]=1$. Außerdem hält er eine *numerische Beobachtung* über die
QR-Routine von `torch.linalg.qr` fest (eine Aussage über diese Routine, kein Satz): Ohne die Vorzeichenkorrektur liegt
$\mathbb E\operatorname{Tr}Q'$ für ihren Faktor $Q'$ viele Standardfehler von $0$ entfernt, er ist also nicht
Haar-verteilt.

```python
k, B = 3, 100_000
Zr = torch.randn(B, k, k, dtype=torch.float64)
Q, R = torch.linalg.qr(Zr)                              # LAPACK: diag(R) of either sign
O = Q * torch.sign(R.diagonal(dim1=-2, dim2=-1)).unsqueeze(-2)     # Theorem B.17
assert close(O.mT @ O, torch.eye(k, dtype=torch.float64).expand(B, k, k))
assert (R.diagonal(dim1=-2, dim2=-1) * torch.sign(R.diagonal(dim1=-2, dim2=-1)) > 0).all()
second = torch.einsum('bij,blm->ijlm', O, O) / B                   # E O_ij O_lm
expected = torch.einsum('il,jm->ijlm', torch.eye(k, dtype=torch.float64), torch.eye(k, dtype=torch.float64)) / k
err = (torch.einsum('bij,blm->ijlm', O ** 2, O ** 2) / B - second ** 2).clamp(min=0).sqrt() / math.sqrt(B)
assert ((second - expected).abs() <= 5 * err).all()
tr, tr_raw = O.diagonal(dim1=-2, dim2=-1).sum(-1), Q.diagonal(dim1=-2, dim2=-1).sum(-1)
print(f"E Tr O = {tr.mean():.4f}, E (Tr O)^2 = {(tr ** 2).mean():.4f};  without the sign correction "
      f"E Tr Q = {tr_raw.mean():.4f}")
assert abs(tr.mean()) < 5 * tr.std() / math.sqrt(B)                # E Tr O = 0
assert abs((tr ** 2).mean() - 1) < 5 * (tr ** 2).std() / math.sqrt(B)   # E (Tr O)^2 = sum_i 1/k = 1
assert abs(tr_raw.mean()) > 20 * tr_raw.std() / math.sqrt(B)       # LAPACK's signs: not Haar
```

**Zurück zu:** [Kapitel 4](04_symmetries.md) (Satz 4.14), [Kapitel 5](05_random_matrices.md) (Satz 5.12),
[Kapitel 9](09_spin_chains.md) (Korrelationen).
