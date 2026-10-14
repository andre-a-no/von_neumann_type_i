[English](../13_jordan_algebras.md) · [Español](../es/13_jordan_algebras.md) · [Français](../fr/13_jordan_algebras.md) · **Deutsch** · [中文](../zh/13_jordan_algebras.md) · [日本語](../ja/13_jordan_algebras.md) · [Русский](../ru/13_jordan_algebras.md)

# 13. Jordan-Algebren: warum die Quantenmechanik *-Algebren verwendet

Das Produkt $xy$ zweier Observablen ist keine Observable, es sei denn, $x$ und $y$ kommutieren, das symmetrisierte
Produkt $x\circ y=\frac12(xy+yx)$ dagegen immer. 1933 schlug P. Jordan vor, die Observablen mit diesem Produkt als
Grundstruktur der Quantenmechanik zu nehmen, und Jordan, von Neumann und Wigner klassifizierten die
endlichdimensionalen Algebren, die dabei auftreten. Dieses Kapitel entwickelt die endlichdimensionale Theorie, mit
vollständigen Beweisen hier oder, für die vier langen (die Klassifikation, den Satz von Koecher–Vinberg, den Satz von
Kadison und die Jordan-Automorphismen von $\bigoplus_cM_{k_c}$), in den Anhängen E und F, und fragt, was der
Jordan-Struktur entgeht. Die Antwort ist das Vorzeichen des Kommutators: Die Jordan-Algebra eines Systems bestimmt den
Kommutator nur bis auf das Vorzeichen (Abschnitt 13.5), zusammengesetzte Systeme brauchen aber das Vorzeichen: Es
entscheidet, welche gemeinsamen Operatoren positiv sind, und der daraus gebildete PPT-Test weist Verschränkung nach
(Proposition 13.12).

**Sie lernen** das Jordan-Produkt und die quadratische Darstellung $P(y)x=yxy$ kennen, die im Gegensatz zu $\circ$
die Positivität erhält; Spinfaktoren, ihre Spektraltheorie und den Lorentz-Kegel; Jordan-Unteralgebren hermitescher
Matrizen, die keine selbstadjungierten Teile von $*$-Algebren sind; die Klassifikation von Jordan–von Neumann–Wigner;
warum die Transposition ein Jordan-Automorphismus und kein $*$-Automorphismus ist und was sie auf einem
zusammengesetzten System bewirkt; und eine Spurungleichung aus Kapitel 11, in Jordan-Sprache geschrieben und in einem
Spinfaktor bewiesen.

Die numerischen Überprüfungen dieses Kapitels stehen in [Tutorial 13](../../tutorial/de/13_jordan_algebras.md).

Im Folgenden ist $\mathrm{Herm}_n(\mathbb K)$ der reelle Vektorraum der hermiteschen $n\times n$-Matrizen über
$\mathbb K=\mathbb R$ oder $\mathbb C$ (Abschnitt 13.4 verwendet auch die Quaternionen $\mathbb H$ und die Oktonionen
$\mathbb O$, konstruiert in Abschnitt E.6, und $\mathrm{Herm}_n$ über ihnen, Abschnitt E.7), $M_{sa}$ der
selbstadjungierte Teil von $\mathcal M=\bigoplus_cM_{k_c}(\mathbb C)$,
$X,Y,Z$ sind die Pauli-Matrizen und $\sigma=(X,Y,Z)$, also $u\cdot\sigma=u_1X+u_2Y+u_3Z$ für $u\in\mathbb R^3$.
$[a,b]=ab-ba$, und $|u|$ ist die euklidische Norm von $u\in\mathbb R^m$.

<a id="s-13-1"></a>

## 13.1 Observablen und das Jordan-Produkt

**Definition.** Eine **Jordan-Algebra** ist ein reeller Vektorraum $J$ mit einem bilinearen, kommutativen Produkt
$x\circ y$, das die **Jordan-Identität** $(x^2\circ y)\circ x=x^2\circ(y\circ x)$ erfüllt, wobei $x^2=x\circ x$. Sie
heißt **speziell**, wenn sie zu einem Unterraum einer assoziativen Algebra isomorph ist, der unter $\frac12(xy+yx)$
abgeschlossen ist. Auf $J$ sei $L(y)x=y\circ x$; die **quadratische Darstellung** ist $P(y)=2L(y)^2-L(y^2)$.

<a id="r-13-1"></a>

**Lemma 13.1.** Sei $\mathcal A$ eine assoziative reelle oder komplexe Algebra und $x\circ y=\frac12(xy+yx)$. Dann
erfüllt $\circ$ die Jordan-Identität, und die Jordan-Potenzen $x^{\circ1}=x$, $x^{\circ(k+1)}=x\circ x^{\circ k}$
stimmen mit den assoziativen Potenzen $x^k$ überein. Insbesondere sind $M_{sa}$, $\mathrm{Herm}_n(\mathbb C)$
und $\mathrm{Herm}_n(\mathbb R)$ spezielle Jordan-Algebren.

*Beweis.* Ausmultiplizieren ergibt $4(x^2\circ y)\circ x=x^2yx+yx^3+x^3y+xyx^2$ und
$4x^2\circ(y\circ x)=x^2yx+x^3y+yx^3+xyx^2$, dieselben vier Terme. Durch Induktion ist
$x\circ x^k=\frac12(x^{k+1}+x^{k+1})=x^{k+1}$. Sind $x,y$ hermitesch, so ist $(xy+yx)^*=yx+xy$, also sind die drei Räume
unter $\circ$ abgeschlossen; für reelle symmetrische $x,y$ ist das Produkt reell. $\square$

*Rechnung:* Die Jordan-Identität und die Potenzassoziativität in $M_2\oplus M_3$ — [Tutorial, Abschnitt 13.1](../../tutorial/de/13_jordan_algebras.md#s-13-1).

<a id="r-13-2"></a>

**Proposition 13.2.** Seien $x,y,z\in M_{sa}$.

1. $P(y)x=yxy$, und es gilt die **Fundamentalformel** $P(P(y)x)=P(y)P(x)P(y)$.
2. $\mathrm{Tr}\,P(y)x=\mathrm{Tr}\,y^2x$.
3. Ist $x\ge0$, so ist $P(y)x\ge0$. Ist $y$ invertierbar, so ist $P(y)$ eine Bijektion des positiven Kegels auf sich,
   mit Inverser $P(y^{-1})$.
4. Das Jordan-Produkt positiver Elemente muss nicht positiv sein: Für $p=\begin{pmatrix}1&0\\0&0\end{pmatrix}$ und
   $q=\frac12\begin{pmatrix}1&1\\1&1\end{pmatrix}$ ist $\det(p\circ q)=-\frac1{16}$.

*Beweis.* (1) $2y\circ(y\circ x)=\frac12(y^2x+2yxy+xy^2)$ und $y^2\circ x=\frac12(y^2x+xy^2)$; man subtrahiere. Dann
ist $P(P(y)x)z=yxy\,z\,yxy=P(y)P(x)P(y)z$ wegen der Assoziativität. (2) Zyklizität der Spur. (3)
$yxy=(x^{1/2}y)^*(x^{1/2}y)\ge0$ (Korollar 1.7), und $P(y)P(y^{-1})=P(y^{-1})P(y)=\mathrm{id}$. (4)
$p\circ q=\frac14\begin{pmatrix}2&1\\1&0\end{pmatrix}$. $\square$

*Rechnung:* Proposition 13.2 in $M_2\oplus M_3$ und das Beispiel (4) — [Tutorial, Abschnitt 13.1](../../tutorial/de/13_jordan_algebras.md#s-13-1). *In der Bibliothek:* `jordan.jordan_product`, `jordan.jordan_power` und `jordan.quadratic_representation`, gebündelt über die Sektoren einer `TypeIAlgebra` — [Tutorial, Abschnitt 13.1](../../tutorial/de/13_jordan_algebras.md#s-13-1).

Schreibt man $P(y)$ als $2L(y)^2-L(y^2)$, so ist jede Identität in (1) allein mit $\circ$ ausgedrückt. Der Beweis von
(1) benutzt nur die Assoziativität des Matrixprodukts; in jeder speziellen Jordan-Algebra $J\subseteq\mathcal A$
($\mathcal A$ assoziativ) liefert er daher $P(y)x=yxy$ und die Fundamentalformel $P(P(y)x)=P(y)P(x)P(y)$ auf $J$. (Das
deckt alle Algebren ab, in denen der Kurs sie verwendet; für nicht spezielle Algebren wird die Fundamentalformel nicht
benutzt.)

Unter den Produkten $\alpha xy+\beta yx$ ($\alpha,\beta\in\mathbb R$) auf $\mathrm{Herm}_n(\mathbb C)$, $n\ge2$,
ist $\frac12(xy+yx)$ das einzige, das in $x,y$ symmetrisch ist und für $y=x$ gleich $x^2$ ist: Symmetrie bedeutet
$(\alpha-\beta)(xy-yx)=0$ für alle $x,y$, und $XZ\ne ZX$ (im linken oberen $2\times2$-Block) erzwingt $\alpha=\beta$;
dann liefert $(\alpha+\beta)1=1$ für $x=y=1$ die Werte $\alpha=\beta=\frac12$.

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). $x\circ y$ ist das symmetrisierte
> Produkt von Observablen (die einzige solche Kombination, wie oben gezeigt), und es ist wieder eine Observable. Die quadratische Darstellung ist die **Lüders-Regel**: Misst man eine Ja/Nein-Frage $p$ und behält das
> Ergebnis „ja“, so geht eine Dichte $\rho$ in $P(p)\rho/\mathrm{Tr}\,P(p)\rho=p\rho p/\mathrm{Tr}(p\rho)$
> über. Proposition 13.2(4) besagt, dass $\circ$ selbst für positive Elemente keine solche operationelle Bedeutung hat.

<a id="s-13-2"></a>

## 13.2 Spinfaktoren und der Lorentz-Kegel

**Definition.** Für $m\ge1$ hat der **Spinfaktor** $V_m=\mathbb R\oplus\mathbb R^m$ die Elemente $x=(a,u)$, das
Produkt $(a,u)\circ(b,v)=(ab+u\cdot v,\;av+bu)$ und die Eins $1=(1,0)$. Seine **Spur** ist $\mathrm{tr}\,x=2a$.

Das Produkt ist bilinear und kommutativ; die Jordan-Identität folgt aus einer Matrixdarstellung.

<a id="r-13-3"></a>

**Lemma 13.3 (antikommutierende Symmetrien).** Seien $p=\lfloor m/2\rfloor$ und $D=2^p$ ($D=2$ für $m=1$). Auf
$(\mathbb C^2)^{\otimes p}$ setze man $s_{2j-1}=Z^{\otimes(j-1)}\otimes X\otimes1$, $s_{2j}=Z^{\otimes(j-1)}\otimes Y\otimes1$
für $j\le p$ und $s_m=Z^{\otimes p}$, falls $m$ ungerade ist ($s_1=Z$ für $m=1$). Dann gilt $s_i=s_i^*$ und
$s_is_j+s_js_i=2\delta_{ij}1$. Die Abbildung $\phi(a,u)=a1+\sum_iu_is_i$ ist injektiv und erfüllt
$\phi(x\circ y)=\phi(x)\circ\phi(y)$; folglich ist $V_m$ eine spezielle Jordan-Algebra. Für $m\ge2$ erzeugen die
$s_i$ die $*$-Algebra $M_D(\mathbb C)$.

*Beweis.* Pauli-Matrizen quadrieren zu $1$ und antikommutieren paarweise; dasselbe gilt für $Z$ mit $X$ und $Y$. Zwei
verschiedene Ketten $s_i,s_j$ unterscheiden sich an genau einem Platz in antikommutierender Weise: für
$s_{2j-1},s_{2j}$ am Platz $j$ ($X$ gegen $Y$); für Ketten mit Indizes $j<j'$ am Platz $j$ ($X$ oder $Y$ gegen $Z$),
während sie an den übrigen Plätzen gleiche Faktoren oder einen Faktor $1$ tragen; dasselbe gilt für
$s_m=Z^{\otimes p}$. Also ist $s_is_j=-s_js_i$ für $i\ne j$ und $s_i^2=1$. Somit ist
$(u\cdot s)\circ(v\cdot s)=\sum_{ij}u_iv_j\,s_i\circ s_j=(u\cdot v)1$, was $\phi(x)\circ\phi(y)=\phi(x\circ y)$ ergibt.
Jedes $s_i$ ist spurfrei (ein Tensorprodukt mit einem spurfreien Faktor), und
$\mathrm{Tr}(s_is_j)=\mathrm{Tr}(s_i\circ s_j)=D\delta_{ij}$, also sind $1,s_1,\dots,s_m$ orthogonal
bezüglich $\mathrm{Tr}(ab)$, und $\phi$ ist injektiv. Die Jordan-Identität gilt in $\mathrm{Herm}_D$
(Lemma 13.1) und wird durch $\phi$ zurücktransportiert. Schließlich ist
$s_{2j-1}s_{2j}=1\otimes\cdots\otimes XY\otimes\cdots\otimes1=iZ_j$ ($Z$ am Platz $j$), also enthält die erzeugte
$*$-Algebra $Z_j$, dann $X_j=Z_1\cdots Z_{j-1}s_{2j-1}$ und $Y_j$, also alle Tensorprodukte von Pauli-Matrizen, eine
Basis von $M_D(\mathbb C)$. $\square$

*In der Bibliothek:* Die $s_i$ sind `jordan.anticommuting_symmetries(m)`, und $V_m$ mit der Abbildung $\phi$ ist die Klasse `SpinFactor` (`to_matrix`) — [Tutorial, Abschnitt 13.2](../../tutorial/de/13_jordan_algebras.md#s-13-2).

Bis auf Vorzeichen sind $s_{2j-1}$ und $s_{2j}$ die Majorana-Operatoren $c_j+c_j^*$ und $i(c_j-c_j^*)$ der
Jordan-Wigner-Fermionen $c_j$ aus Abschnitt 9.4, wobei Platz $j-1$ aus Kapitel 9 der $j$-te Tensorfaktor ist: Dort
ist $s^a=\frac12\sigma_a$, $n=s^z+\frac12$ und $s^\pm=s^x\pm is^y$, also $(-1)^n=1-2n=-Z$ und $s^-=\frac12(X-iY)$,
folglich $c_j=(-1)^{j-1}Z^{\otimes(j-1)}\otimes s^-\otimes1$. Da $s^-+s^+=X$ und $i(s^--s^+)=i(-iY)=Y$, gilt
$$
s_{2j-1}=(-1)^{j-1}(c_j+c_j^*),\qquad s_{2j}=(-1)^{j-1}\,i(c_j-c_j^*).
$$

*Rechnung:* Diese Identitäten für $p=3$ — [Tutorial, Abschnitt 13.2](../../tutorial/de/13_jordan_algebras.md#s-13-2).

<a id="r-13-4"></a>

**Satz 13.4 (Spektraltheorie von $V_m$).** Seien $x=(a,u)\in V_m$ und $\lambda_\pm=a\pm|u|$.

1. Die Idempotente von $V_m$ sind $0$, $1$ und $c_w=\frac12(1,w)$ mit $|w|=1$; und $c_w\circ c_{w'}=0$ genau dann,
   wenn $w'=-w$.
2. Ist $u\ne0$, so ist $x=\lambda_+e_++\lambda_-e_-$ mit $e_\pm=c_{\pm u/|u|}$; ist $u=0$, so ist $x=a1$. Ist
   $x=\mu_1c_w+\mu_2c_{-w}$, so gilt $\{\mu_1,\mu_2\}=\{\lambda_+,\lambda_-\}$. Wir nennen $\lambda_\pm$ die
   **Eigenwerte** von $x$.
3. Für $f:\{\lambda_+,\lambda_-\}\to\mathbb R$ setze man $f(x)=f(\lambda_+)e_++f(\lambda_-)e_-$ ($f(a)1$ für $u=0$).
   Dann ist $x^{\circ k}=\lambda_+^ke_++\lambda_-^ke_-$, $\phi(f(x))=f(\phi(x))$, die Eigenwerte von $\phi(x)$ sind
   $\lambda_\pm$, jeweils mit Vielfachheit $D/2$, und $\mathrm{Tr}\,\phi(x)=\frac D2\mathrm{tr}\,x$. Außerdem
   ist $\mathrm{tr}\,x=\lambda_++\lambda_-$ und $\det x:=\lambda_+\lambda_-=a^2-|u|^2$.
4. Der Kegel der Quadrate $Q=\{x^2:x\in V_m\}$ ist der **Lorentz-Kegel** $\{(a,u):a\ge|u|\}=\{x:\lambda_-\ge0\}=\{x:\phi(x)\ge0\}$.
   Er ist **selbstdual**: $y\in Q$ genau dann, wenn $\mathrm{tr}(x\circ y)\ge0$ für alle $x\in Q$.
5. $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$; folglich ist $P(y)Q\subseteq Q$ für jedes $y\in V_m$.

*Beweis.* (1) $(a,u)^2=(a^2+|u|^2,2au)$ ist genau dann gleich $(a,u)$, wenn $u=0$, $a\in\{0,1\}$, oder $a=\frac12$,
$|u|=\frac12$. $c_w\circ c_{w'}=\frac14(1+w\cdot w',\,w+w')$ verschwindet genau dann, wenn $w'=-w$. (2)
$\lambda_+e_++\lambda_-e_-=(a,0)+|u|(0,u/|u|)=x$, und $\mu_1c_w+\mu_2c_{-w}=(\frac{\mu_1+\mu_2}2,\frac{\mu_1-\mu_2}2w)$,
also $\mu_1+\mu_2=2a$ und $|\mu_1-\mu_2|=2|u|$. (3) Nach (1) ist $(\alpha e_++\beta e_-)\circ(\gamma e_++\delta e_-)=\alpha\gamma
e_++\beta\delta e_-$, was durch Induktion die Potenzen liefert. Mit $w=u/|u|$ ist $(w\cdot s)^2=1$, also sind
$\phi(e_\pm)=\frac12(1\pm w\cdot s)$ Orthogonalprojektionen mit Summe $1$ und Produkt $0$, jede mit Spur $D/2$; somit
ist $\phi(x)=\lambda_+\phi(e_+)+\lambda_-\phi(e_-)$ eine Spektralzerlegung, und
$f(\phi(x))=f(\lambda_+)\phi(e_+)+f(\lambda_-)\phi(e_-)=\phi(f(x))$ (Satz 1.5). (4) $x^2$ hat die Eigenwerte
$\lambda_\pm^2\ge0$; ist umgekehrt $\lambda_-\ge0$, so ist $x=(\sqrt x)^2$. Es ist $\lambda_-\ge0$ genau dann, wenn
$a\ge|u|$, und nach (3) genau dann, wenn $\phi(x)\ge0$. Zur Selbstdualität: $\mathrm{tr}(x\circ y)=2(ab+u\cdot v)$.
Ist $a\ge|u|$ und $b\ge|v|$, so ist $ab\ge|u||v|\ge-u\cdot v$. Ist $b<|v|$, so liefert $x=(|v|,-v)\in Q$ den Wert
$2(|v|b-|v|^2)<0$ (für $v=0$, $b<0$ nehme man $x=1$). (5) $P(y)x=2y\circ(y\circ x)-y^2\circ x$, und $\phi$ erhält
$\circ$, also ist $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$ nach Proposition 13.2(1); man wende (4) und Proposition 13.2(3)
an. $\square$

*Rechnung:* Satz 13.4 in $V_5$: die Spektralzerlegung, die Eigenwerte von $\phi(x)$, $\phi(e^x)=e^{\phi(x)}$ und der Kegel der Quadrate — [Tutorial, Abschnitt 13.2](../../tutorial/de/13_jordan_algebras.md#s-13-2). *In der Bibliothek:* `SpinFactor.spectral_decomposition`, `exp`, `sqrt`, `in_cone` — [Tutorial, Abschnitt 13.2](../../tutorial/de/13_jordan_algebras.md#s-13-2).

<a id="r-13-5"></a>

**Proposition 13.5.** (a) Für $m=3$ ist $\phi(a,u)=a1+u\cdot\sigma$ ein Jordan-Isomorphismus
$V_3\cong\mathrm{Herm}_2(\mathbb C)$. (b) Die Jordan-Automorphismen von $V_m$ sind die Abbildungen
$(a,u)\mapsto(a,Ou)$, $O\in O(m)$. (c) Für $m\ge2$ hat $\phi(V_m)$ die Dimension $m+1$, und der selbstadjungierte Teil
der von ihm erzeugten $*$-Algebra hat die Dimension $D^2=4^{\lfloor m/2\rfloor}$; gleich sind sie nur für $m=3$.

*Beweis.* (a) Für $m=3$ ist $(s_1,s_2,s_3)=(X,Y,Z)$, und $\phi$ ist injektiv zwischen Räumen der Dimension $4$. (b)
Die Abbildungen $(a,u)\mapsto(a,Ou)$ erhalten $u\cdot v$, also $\circ$. Sei umgekehrt $\alpha$ ein
Jordan-Automorphismus. Dann ist $\alpha(1)$ eine Eins von $V_m$, also $\alpha(1)=1$. Nach der Formel für $x^2$ gilt
$x^2=1$ genau dann, wenn $x=\pm1$ oder $x=(0,w)$ mit $|w|=1$. Für $|w|=1$ ist $\alpha(0,w)^2=1$ und
$\alpha(0,w)\ne\pm1=\alpha(\pm1)$, also $\alpha(0,w)=(0,w')$ mit $|w'|=1$. Schreibt man $\alpha(0,u)=(\ell(u),Ou)$ mit
linearen $\ell,O$, so besagt dies $\ell=0$ und $|Ou|=|u|$, d. h. $O\in O(m)$. (c) Lemma 13.3;
$4^{\lfloor m/2\rfloor}>m+1$ für $m=2$ und $m\ge4$. $\square$

Für $m=1$ erzeugt $s_1=Z$ die Diagonalmatrizen, der Dimension $2=m+1$, und $V_1\cong\mathbb R\oplus\mathbb R$ ist
nicht einfach.

*Rechnung:* (b) für ein zufälliges $O\in O(5)$ und (c) für $m\le6$ mit `jordan.compare_with_star_algebra` — [Tutorial, Abschnitt 13.2](../../tutorial/de/13_jordan_algebras.md#s-13-2).

**Zustände von $V_m$.** Ein **Zustand** von $V_m$ ist ein lineares Funktional $\omega$ mit $\omega(Q)\ge0$ und
$\omega(1)=1$. Die Zustände sind genau die Abbildungen $x\mapsto\mathrm{tr}(\rho\circ x)$ mit
$\rho=\frac12(1,r)$, $|r|\le1$; der Zustandsraum ist also eine Kugel der Dimension $m$. *Beweis.*
$\langle x,y\rangle=\mathrm{tr}(x\circ y)=2(ab+u\cdot v)$ ist ein Skalarprodukt auf $V_m$, also ist jedes lineare
$\omega$ von der Form $\mathrm{tr}(\rho\circ\cdot)$ für genau ein $\rho=(b,v)$ (endlichdimensionale
Riesz-Darstellung: in einer Orthonormalbasis ist $\rho=\sum_k\omega(e_k)e_k$). Nach der Selbstdualität in Satz
13.4(4) gilt $\omega(Q)\ge0$ genau dann, wenn $\rho\in Q$, d. h. $b\ge|v|$; und $\omega(1)=\mathrm{tr}\,\rho=2b=1$.
Also ist $\rho=\frac12(1,r)$ mit $r=2v$, $|r|\le1$; umgekehrt liefert jedes solche $\rho$ einen Zustand. $\square$

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Für $m=3$ ist der Zustandsraum von
> $V_m$ die Bloch-Kugel eines Qubits (Proposition 13.5(a)). Der Lorentz-Kegel ist der Zukunftslichtkegel der
> speziellen Relativitätstheorie in $1+m$ Dimensionen, und
> Proposition 13.5(b) besagt, dass die Jordan-Automorphismen die orthogonalen Abbildungen des Raumes sind, die die
> Zeitachse festhalten.

Die Kegel $Q$ sind die Nebenbedingungsmengen der **Kegelprogrammierung zweiter Ordnung** (SOCP): lineare
Optimierung über Produkten von Lorentz-Kegeln. Die Größen aus Satz 13.4 (Eigenwerte, $\det$ und $P(y)$) sind
diejenigen, in denen Innere-Punkte-Verfahren für SOCP üblicherweise formuliert werden; diese Algorithmen werden in
diesem Kurs nicht behandelt (weiterführende Literatur: F. Alizadeh, D. Goldfarb, *Second-order cone programming*, Math.
Program. 95 (2003) 3–51).

<a id="s-13-3"></a>

## 13.3 Jordan-Unteralgebren, die keine *-Unteralgebren sind

Eine Jordan-Unteralgebra $J\subseteq\mathrm{Herm}_n(\mathbb C)$ erzeugt eine komplexe $*$-Algebra $\mathcal A$
mit $J\subseteq\mathcal A_{sa}$. Proposition 13.5(c) lieferte Beispiele mit $J\ne\mathcal A_{sa}$; hier sind zwei
weitere.

<a id="r-13-6"></a>

**Proposition 13.6.** (a) $\mathrm{Herm}_n(\mathbb R)$ ist eine Jordan-Unteralgebra von
$\mathrm{Herm}_n(\mathbb C)$ der Dimension $n(n+1)/2$; die von ihr erzeugte $*$-Algebra ist $M_n(\mathbb C)$,
deren selbstadjungierter Teil die Dimension $n^2$ hat. (b) Zu $x,y\in\mathrm{Herm}_2(\mathbb C)$ gibt es eine
unitäre Matrix $U$, sodass $UxU^*$ und $UyU^*$ reell symmetrisch sind. Die von $1,x,y$ erzeugte Jordan-Algebra ist
$\mathrm{span}\{1,x,y\}$; sind $x,y,1$ linear unabhängig, so ist sie $U^*\mathrm{Herm}_2(\mathbb R)U\cong V_2$,
der Dimension $3$, während die von $1,x,y$ erzeugte $*$-Algebra $M_2(\mathbb C)$ ist, der Dimension $4$.

*Beweis.* (a) Lemma 13.1. Für $i\ne j$ ist $E_{ij}=E_{ii}(E_{ij}+E_{ji})$ ein Produkt reeller symmetrischer Matrizen,
also enthält die erzeugte Algebra alle Matrixeinheiten. (b) Man schreibe $x=a1+u\cdot\sigma$, $y=b1+v\cdot\sigma$ und
wähle einen Einheitsvektor $n\perp u,v$. Die Matrizen $n\cdot\sigma$ und $Y$ sind hermitesch, spurfrei und quadrieren
zu $1$, also haben beide die Eigenwerte $1,-1$, und $U(n\cdot\sigma)U^*=Y$ für ein unitäres $U$ (Satz 1.3). Konjugation
mit $U$ erhält $\mathrm{Tr}$ und $\circ$, bildet also $u\cdot\sigma$ auf $(Ru)\cdot\sigma$ ab, mit linearem $R$
und $Ru\cdot Rv=\frac12\mathrm{Tr}\big((Ru\cdot\sigma)(Rv\cdot\sigma)\big)=u\cdot v$: $R\in O(3)$, und $Rn=e_2$.
Daher ist $Ru,Rv\perp e_2$, und $UxU^*=a1+(Ru)_1X+(Ru)_3Z$ ist reell, da $X,Z$ es sind; ebenso $UyU^*$. Ferner ist
$x\circ y=(ab+u\cdot v)1+av\cdot\sigma+bu\cdot\sigma\in\mathrm{span}\{1,x,y\}$, und ebenso $x^2,y^2$, also ist der
Spann abgeschlossen. Hat er die Dimension $3$, so ist sein Konjugiertes unter $U$ ganz
$\mathrm{Herm}_2(\mathbb R)=\mathrm{span}\{1,X,Z\}$, das vermöge $(a,u)\mapsto a1+u_1X+u_2Z$ zu $V_2$
isomorph ist (der Beweis von Lemma 13.3 benutzt nur, dass die $s_i$ antikommutierende hermitesche Symmetrien sind, und
$X,Z$ sind solche) und nach (a) $M_2(\mathbb C)$ erzeugt. $\square$

*Rechnung:* Proposition 13.6(a) für $n\le5$ und (b) für zufällige hermitesche $2\times2$-Matrizen, mit `jordan.compare_with_star_algebra` — [Tutorial, Abschnitt 13.3](../../tutorial/de/13_jordan_algebras.md#s-13-3).

Drei generische hermitesche $2\times2$-Matrizen erzeugen ganz $\mathrm{Herm}_2(\mathbb C)$. In einem Block $M_2$ von
$\mathcal M$ verraten zwei Observablen die komplexe Struktur also nie; eine dritte kann es.

<a id="s-13-4"></a>

## 13.4 Die Klassifikation

Eine endlichdimensionale Jordan-Algebra heißt **formal reell**, wenn aus $\sum_ix_i^2=0$ folgt, dass alle $x_i=0$ sind,
und **euklidisch**, wenn sie ein Einselement $1$ und ein Skalarprodukt mit
$\langle x\circ y,z\rangle=\langle y,x\circ z\rangle$ besitzt. Die beiden Begriffe stimmen überein (Satz E.5(4); die
Richtung „euklidisch $\Rightarrow$ formal reell“ ist eine Zeile:
$0=\sum_i\langle x_i^2,1\rangle=\sum_i\langle x_i,x_i\rangle$). $M_{sa}$ mit $\langle x,y\rangle=\mathrm{Tr}(xy)$
ist euklidisch: $\mathrm{Tr}(x^2)=\sum_{ij}|x_{ij}|^2$, und wegen der Zyklizität ist
$\mathrm{Tr}((x\circ y)z)=\frac12\mathrm{Tr}(xyz+yxz)=\mathrm{Tr}(y(x\circ z))$. Ebenso $V_m$ mit
$\langle x,y\rangle=\mathrm{tr}(x\circ y)=2(ab+u\cdot v)$: Für $z=(c,w)$ ist
$\mathrm{tr}((x\circ y)\circ z)=2(abc+c\,u\cdot v+a\,v\cdot w+b\,u\cdot w)$ symmetrisch in $x,y,z$.

<a id="r-13-7"></a>

**Satz 13.7 (Jordan–von Neumann–Wigner).** Jede endlichdimensionale formal reelle Jordan-Algebra ist eine direkte
Summe einfacher Ideale, und jede einfache ist zu genau einer der folgenden isomorph: $\mathbb R$;
$\mathrm{Herm}_r(\mathbb K)$ für $\mathbb K=\mathbb R,\mathbb C,\mathbb H$ und $r\ge3$; den Spinfaktoren $V_m$,
$m\ge2$; der **Albert-Algebra** $\mathrm{Herm}_3(\mathbb O)$ der hermiteschen $3\times3$-Matrizen über den
Oktonionen (Dimension $27$). Außerdem gilt $\mathrm{Herm}_1(\mathbb K)\cong\mathbb R$ und
$\mathrm{Herm}_2(\mathbb K)\cong V_{\dim\mathbb K+1}$ für $\mathbb K=\mathbb R,\mathbb C,\mathbb H,\mathbb O$.

*Beweis.* Anhang E beweist ihn vollständig: Die Zerlegung ist Satz E.9, die Liste und „genau eine“ sind Satz
E.27(1),(2), und die Isomorphismen in Rang eins und zwei sind Satz E.24(3). $\square$

(P. Jordan, J. von Neumann, E. Wigner, *On an algebraic generalization of the quantum mechanical formalism*, Ann. of
Math. 35 (1934) 29–64; der Beweis in Anhang E folgt [FK94, Kapitel II–V].)

Alle außer der Albert-Algebra sind speziell: $\mathbb R=\mathrm{Herm}_1(\mathbb R)$ und
$\mathrm{Herm}_r(\mathbb K)$, $\mathbb K=\mathbb R,\mathbb C$, nach Lemma 13.1, $\mathrm{Herm}_r(\mathbb H)$,
weil das Schreiben von Quaternionen als komplexe $2\times2$-Matrizen es in $\mathrm{Herm}_{2r}(\mathbb C)$
einbettet (Satz E.24(1)), und $V_m$ nach Lemma 13.3. Die Albert-Algebra ist eine Jordan-Algebra (Satz E.24(2)) und ist
**exzeptionell**: Sie ist zu keiner Jordan-Unteralgebra einer assoziativen reellen oder komplexen Algebra isomorph
(Satz E.26, ein Argument von der Art, wie es A. A. Albert verwendete, *On a certain algebra of quantum mechanics*, Ann.
of Math. 35 (1934) 65–73: Eine Jordan-Einbettung in eine assoziative Algebra würde das Oktonionenprodukt assoziativ
machen, was es nicht ist, Lemma E.15(2)). Für $r\ge4$ erfüllt $\mathrm{Herm}_r(\mathbb O)$ die Jordan-Identität
nicht (Korollar E.25), und eine einfache Algebra mit einem Jordan-Rahmen aus $r\ge4$ Elementen ist
$\mathrm{Herm}_r(D)$ mit einer assoziativen Kompositionsalgebra $D\ne\mathbb O$ (Satz E.22), weshalb die
Oktonionen nur für $r\le3$ auftreten. Daher lässt sich $\mathrm{Herm}_3(\mathbb O)$ nicht treu durch reelle
oder komplexe Matrizen (oder innerhalb irgendeiner assoziativen Algebra) mit dem Produkt $\frac12(xy+yx)$ darstellen;
eine Rechnung in ihr muss statt eines Matrixprodukts das Oktonionenprodukt (oder seine Strukturkonstanten) verwenden.

Sei $\Omega$ ein nichtleerer offener konvexer Kegel in einem euklidischen Raum $V$ (einem endlichdimensionalen reellen
Skalarproduktraum $\ne\{0\}$). Sein **dualer Kegel** ist
$\Omega^*=\{y:\langle x,y\rangle>0\text{ für alle }x\in\bar\Omega\setminus\{0\}\}$; $\Omega$ heißt **selbstdual**,
wenn $\Omega^*=\Omega$, **homogen**, wenn die linearen Abbildungen $g$ mit $g\Omega=\Omega$ transitiv auf $\Omega$
wirken, und **symmetrisch**, wenn beides gilt.

<a id="r-13-8"></a>

**Satz 13.8 (Koecher–Vinberg).** Ein nichtleerer offener konvexer Kegel $\Omega$ in einem euklidischen Raum $V$ ist
genau dann symmetrisch, wenn es ein Produkt $\circ$ auf $V$ gibt, sodass $(V,\circ,\langle\cdot,\cdot\rangle)$ eine
euklidische Jordan-Algebra und $\Omega$ das Innere ihres Kegels der Quadrate ist.

*Beweis.* Satz F.12 („wenn“: Sätze F.5(2) und F.6(2)). $\square$

(M. Koecher, 1957; E. B. Vinberg, 1960; der Beweis in Anhang F folgt [FK94, Kapitel I–III].) Satz 13.4(4) beweist die
abgeschlossene Form der Selbstdualität für den Lorentz-Kegel ($y\in Q$ genau dann, wenn $\mathrm{tr}(x\circ y)\ge0$
für alle $x\in Q$), und Aufgabe 1 seine Homogenität. Der Kegel der Quadrate einer direkten Summe euklidischer
Jordan-Algebren (mit der orthogonalen Summe der Skalarprodukte) ist das Produkt ihrer Kegel der Quadrate, da
$(x_1,\dots,x_k)^2=(x_1^2,\dots,x_k^2)$. Nach der Richtung „wenn“ sind daher die Inneren der Kegel der semidefiniten
Programmierung (Produkte von $\mathrm{Herm}_n(\mathbb R)^+$ oder $\mathrm{Herm}_n(\mathbb C)^+$, nach Satz
1.6 die Kegel der Quadrate; für reelles $x\ge0$ ist die Quadratwurzel reell, da ihre komplex Konjugierte ebenfalls eine
positive Quadratwurzel von $x$ ist), der SOCP (Produkte von Lorentz-Kegeln) und der linearen Programmierung
($\mathbb R_+^n$, der Kegel der Quadrate von $\mathbb R^n=\mathbb R\oplus\dots\oplus\mathbb R$) symmetrische Kegel.

[FK94]: J. Faraut, A. Korányi, *Analysis on Symmetric Cones*, Oxford University Press, 1994.

<a id="s-13-5"></a>

## 13.5 Jordan-Abbildungen und *-Abbildungen

Ein **Jordan-$*$-Homomorphismus** von von-Neumann-Algebren ist eine komplex-lineare Abbildung $\Phi$ mit
$\Phi(x^*)=\Phi(x)^*$ und $\Phi(x\circ y)=\Phi(x)\circ\Phi(y)$; auf den selbstadjungierten Teilen ist er ein
Jordan-Homomorphismus, und jeder Jordan-Homomorphismus selbstadjungierter Teile setzt sich durch
$\Phi(x+iy)=\Phi(x)+i\Phi(y)$ ($x,y$ selbstadjungiert) zu einem solchen fort. (Dieses $\Phi$ ist komplex-linear und
erfüllt $\Phi(x-iy)=\Phi(x+iy)^*$; und da $\circ$ komplex-bilinear ist und Paare selbstadjungierter Elemente auf
selbstadjungierte Elemente abbildet, wird $(x+iy)\circ(x'+iy')=(x\circ x'-y\circ y')+i(x\circ y'+y\circ x')$ auf
$\Phi(x+iy)\circ\Phi(x'+iy')$ abgebildet.)

<a id="r-13-9"></a>

**Proposition 13.9 (die Transposition).** Sei $T(A)=\bigoplus_cA_c^T$ auf $\mathcal M$. (a) $T$ ist ein unitaler
Jordan-$*$-Automorphismus mit $T(AB)=T(B)T(A)$; er erhält Spektren, also sind $T$ und $T^{-1}=T$ positiv. (b) Auf
$M_k$, $k\ge2$, ist die Choi-Matrix $\sum_{ij}E_{ij}\otimes T(E_{ij})$ der Swap $F(e_a\otimes e_b)=e_b\otimes e_a$,
mit den Eigenwerten $+1$ (Vielfachheit $k(k+1)/2$) und $-1$ (Vielfachheit $k(k-1)/2$). Also ist $T$ nicht vollständig
positiv.

*Beweis.* (a) $(A^T)^*=(A^*)^T$, $(AB)^T=B^TA^T$, also $(A\circ B)^T=B^T\circ A^T=A^T\circ B^T$;
$\det(A^T-\lambda)=\det(A-\lambda)$. (b) $\sum_{ij}E_{ij}e_a\otimes E_{ji}e_b=e_b\otimes e_a$. $F^2=1$, und $F$ ist $+1$
auf den symmetrischen und $-1$ auf den antisymmetrischen Tensoren, die zusammen den Raum aufspannen. Nach Satz 6.2
(und Proposition 6.1) ist $T$ nicht CP. $\square$

*Rechnung:* Proposition 13.9 mit der Choi-Matrix der Transposition auf $M_3$ — [Tutorial, Abschnitt 13.5](../../tutorial/de/13_jordan_algebras.md#s-13-5). *In der Bibliothek:* Die Transposition $T$ ist `jordan.transpose_map` — [Tutorial, Abschnitt 13.5](../../tutorial/de/13_jordan_algebras.md#s-13-5).

Ein Jordan-Automorphismus muss also kein $*$-Automorphismus sein. Bis auf dieses Phänomen bestimmt die Positivität die
Jordan-Struktur:

<a id="r-13-10"></a>

**Satz 13.10 (Kadison, endliche Dimension).** Seien $\mathcal M=\bigoplus_cM_{k_c}(\mathbb C)$ und
$\mathcal N=\bigoplus_dM_{k'_d}(\mathbb C)$, und sei $\Phi:\mathcal M\to\mathcal N$ eine lineare Bijektion mit
$\Phi(1)=1$, sodass $\Phi$ und $\Phi^{-1}$ positiv sind. Dann ist $\Phi$ ein Jordan-$*$-Isomorphismus. Umgekehrt ist
jeder Jordan-$*$-Isomorphismus $\mathcal M\to\mathcal N$ unital, und er und seine Inverse sind positiv.

*Beweis.* Satz F.15. (Die Idee: $\Phi$ bildet die konvexe Menge $\{0\le a\le1\}$ auf die entsprechende Menge von
$\mathcal N$ ab, also Extrempunkte, und das sind die Projektionen, auf Projektionen; es erhält $p+q\le1$, was für
Projektionen $pq=0$ bedeutet; also bildet es Spektralzerlegungen auf Spektralzerlegungen ab, und
$\Phi(x^2)=\Phi(x)^2$.) $\square$

R. V. Kadison bewies den Satz für alle unitalen C\*-Algebren (*A generalized Schwarz inequality and algebraic
invariants for operator algebras*, Ann. of Math. 56 (1952) 494–503); in diesem Kurs wird nur der endlichdimensionale
Fall bewiesen und verwendet. Für die Faktoren $M_n$ sind die Jordan-Automorphismen genau bekannt: Jeder
Jordan-$*$-Automorphismus von $M_n(\mathbb C)$ ist $A\mapsto UAU^*$ oder $A\mapsto UA^TU^*$ mit unitärem $U$ (Satz
F.17, für alle $n$; das Ergebnis stammt von N. Jacobson, C. E. Rickart, Trans. Amer. Math. Soc. 69 (1950) 479–502,
und I. N. Herstein, Trans. Amer. Math. Soc. 81 (1956) 331–341). Die ersten Abbildungen sind multiplikativ, die zweiten
antimultiplikativ; ein Jordan-$*$-Automorphismus von $M_n$ erhält also entweder alle Kommutatoren oder kehrt alle um:
$\Phi([x,y])=\pm[\Phi(x),\Phi(y)]$. Auf $\bigoplus_cM_{k_c}$ permutiert ein Jordan-$*$-Isomorphismus Blöcke gleicher
Größe und wirkt auf jedem Block auf eine der beiden Arten (Satz F.18); er ist unital und in beiden Richtungen positiv
(Satz 13.10), und nach Korollar F.19 ist er genau dann vollständig positiv, wenn er keinen Block der Größe $\ge2$
transponiert, also genau dann, wenn er ein $*$-Isomorphismus ist. Hier ist ein direkter Beweis des Falls $n=2$ in der
Sprache der Spinfaktoren.

<a id="r-13-11"></a>

**Proposition 13.11.** Jeder Jordan-Automorphismus von $\mathrm{Herm}_2(\mathbb C)$ ist $A\mapsto UAU^*$ oder
$A\mapsto UA^TU^*$ mit $U\in U(2)$.

*Beweis.* Nach Proposition 13.5(a,b) wirkt ein Jordan-Automorphismus als $a1+u\cdot\sigma\mapsto a1+(Ou)\cdot\sigma$
mit $O\in O(3)$. Man schreibe $R_U\in O(3)$ für die Matrix von $A\mapsto UAU^*$ in diesem Sinne; $R_{UV}=R_UR_V$. (i)
$\det R_U=1$: $U=W\mathrm{diag}(e^{i\theta_1},e^{i\theta_2})W^*$ ist durch den Weg
$U_t=W\mathrm{diag}(e^{it\theta_1},e^{it\theta_2})W^*$ mit $1$ verbunden, und $\det R_{U_t}\in\{\pm1\}$ hängt
stetig von $t$ ab. (ii) Jedes $R\in SO(3)$ ist ein $R_U$: Wie in Proposition 13.6 gibt es $U_1$ mit
$U_1ZU_1^*=(Re_3)\cdot\sigma$, d. h. $R_{U_1}e_3=Re_3$. Dann hält $R_{U_1}^{-1}R\in SO(3)$ den Vektor $e_3$ fest, ist
also eine Drehung um einen Winkel $\beta$ um $e_3$. Eine direkte Rechnung ergibt für
$U_2=\mathrm{diag}(e^{-i\beta/2},e^{i\beta/2})$, dass $U_2XU_2^*=\cos\beta\,X+\sin\beta\,Y$ und
$U_2YU_2^*=-\sin\beta\,X+\cos\beta\,Y$, also $R=R_{U_1U_2}$. (iii) Die Transposition hält $X,Z$ fest und bildet
$Y\mapsto-Y$ ab: $O_T=\mathrm{diag}(1,-1,1)$, $\det O_T=-1$. Daher ist $O=R_U$, falls $\det O=1$, und $OO_T=R_U$,
d. h. der Automorphismus ist $A\mapsto UA^TU^*$, falls $\det O=-1$. $\square$

*Rechnung:* $O_T=\mathrm{diag}(1,-1,1)$ — [Tutorial, Abschnitt 13.5](../../tutorial/de/13_jordan_algebras.md#s-13-5).

**Zusammengesetzte Systeme.** Die Jordan-Struktur eines Faktors kann $T$ nicht von der Identität unterscheiden. Bei
zwei Systemen wird diese Mehrdeutigkeit sichtbar:

<a id="r-13-12"></a>

**Proposition 13.12.** Seien $x,x'\in\mathrm{Herm}_n(\mathbb C)$, $y,y'\in\mathrm{Herm}_m(\mathbb C)$. In
$\mathrm{Herm}_{nm}(\mathbb C)$, wirkend auf $\mathbb C^n\otimes\mathbb C^m$, gilt
$$
(x\otimes y)\circ(x'\otimes y')=(x\circ x')\otimes(y\circ y')+\tfrac14[x,x']\otimes[y,y'].
$$
Seien $\theta=\mathrm{id}\otimes T$ (die **partielle Transposition**) und
$A\circ_\theta B=\theta(\theta(A)\circ\theta(B))$. Dann ist $\circ_\theta$ ein Jordan-Produkt auf
$\mathrm{Herm}_{nm}(\mathbb C)$, über $\theta$ isomorph zu $\circ$, das mit $\circ$ auf allen Paaren lokaler
Observablen $x\otimes1$, $1\otimes y$ übereinstimmt, und
$(x\otimes y)\circ_\theta(x'\otimes y')=(x\circ x')\otimes(y\circ y')-\frac14[x,x']\otimes[y,y']$. Sein Kegel der
Quadrate ist $\theta(\mathrm{Herm}_{nm}^+)$. Für $n=m=2$ und $\Omega=e_1\otimes e_1+e_2\otimes e_2$ liegt der
Zustand $\omega=\frac12|\Omega\rangle\langle\Omega|\ge0$ nicht in diesem Kegel: $\theta(\omega)=\frac12F$ hat den
Eigenwert $-\frac12$.

*Beweis.* Mit $j=x\circ x'$, $c=\frac12[x,x']$, $k=y\circ y'$, $d=\frac12[y,y']$ ist $xx'=j+c$, $x'x=j-c$, und ebenso
für $y$, also
$\frac12(xx'\otimes yy'+x'x\otimes y'y)=\frac12\big((j+c)\otimes(k+d)+(j-c)\otimes(k-d)\big)=j\otimes k+c\otimes d$.
$\theta$ ist eine reell-lineare Involution von $\mathrm{Herm}_{nm}$: Sie bildet $x\otimes y$ auf $x\otimes y^T$
ab, und die $x\otimes y$ spannen $\mathrm{Herm}_{nm}$ über $\mathbb R$ auf (die Produkte
$h_\alpha\otimes g_\beta$ reeller Basen von $\mathrm{Herm}_n,\mathrm{Herm}_m$ bilden eine komplexe Basis
von $M_{nm}$, und eine hermitesche Matrix hat darin reelle Koeffizienten). Also ist $\circ_\theta$ ein Jordan-Produkt
und $\theta$ ein Isomorphismus. Da $(y^T\circ y'^T)^T=y\circ y'$ und $[y^T,y'^T]^T=-[y,y']$, ergibt die Anwendung von
$\theta$ auf die Formel für $(x\otimes y^T)\circ(x'\otimes y'^T)$ die Formel für $\circ_\theta$. Auf lokalen
Observablen verschwindet der Kommutatorterm ($[1,1]=0$, $[x,1]=0$). Die Quadrate sind
$A\circ_\theta A=\theta(\theta(A)^2)$, und $\omega$ liegt genau dann in $\theta(\mathrm{Herm}_{nm}^+)$, wenn
$\theta(\omega)\ge0$, da $\theta^2=\mathrm{id}$. Schließlich ist
$\theta(\omega)=\frac12\sum_{ij}E_{ij}\otimes E_{ji}=\frac12F$ (Proposition 13.9(b); das ist die Rechnung aus dem
Beweis von Proposition 6.1). $\square$

*Rechnung:* Die Produktformel, $\circ_\theta$ und der Eigenwert $-\frac12$ von $\theta(\omega)$ — [Tutorial, Abschnitt 13.5](../../tutorial/de/13_jordan_algebras.md#s-13-5).

Die beiden Produkte sind gleich gute Jordan-Algebren, sie stimmen in allem überein, was sich an jeweils einem
Teilsystem messen lässt, und sie unterscheiden sich darin, welche gemeinsamen Operatoren positiv sind, also darin,
was die Zustände des zusammengesetzten Systems sind. Was sie unterscheidet, ist das Vorzeichen des Kommutatorterms,
und dieses Vorzeichen wird durch die $*$-Algebra $M_n\otimes M_m$ festgelegt und nicht durch die Jordan-Algebren der
Faktoren. ($\circ_\theta$ ist das Jordan-Produkt des assoziativen Produkts $A\star B=\theta(\theta(A)\theta(B))$, für
das $(x\otimes y)\star(x'\otimes y')=xx'\otimes y'y$ gilt: die Algebra $M_n\otimes M_m^{\rm op}$.) Der Kegel
$\theta(\mathrm{Herm}_{nm}^+)$ ist es, der das **PPT-Kriterium** liefert (A. Peres, Phys. Rev. Lett. 77 (1996)
1413). Ein Zustand $\rho$ von $M_n\otimes M_m$ heißt **separabel**, wenn $\rho=\sum_kp_ka_k\otimes b_k$ (eine endliche
Summe) mit $p_k\ge0$, $a_k,b_k\ge0$, und andernfalls **verschränkt**. Ein separabler Zustand liegt in beiden Kegeln, da
$\theta(\rho)=\sum_kp_ka_k\otimes b_k^T$ mit $b_k^T\ge0$ (Proposition 13.9(a)) und
$a\otimes b=(a^{1/2}\otimes b^{1/2})^2\ge0$. Also ist $\omega$, das nicht in $\theta(\mathrm{Herm}_4^+)$ liegt,
verschränkt. (Weiterführende Literatur zu Zusammensetzungen von Jordan-Algebren, hier nicht verwendet: H.
Hanche-Olsen, *On the structure and tensor products of JC-algebras*, Canad. J. Math. 35 (1983) 1059–1074; H. Barnum,
M. Graydon, A. Wilce, *Composites and categories of Euclidean Jordan algebras*, Quantum 4 (2020) 359.)

**Dynamik.** Man könnte hoffen, dass die Zeitentwicklung das assoziative Produkt erzwingt. Das tut sie nicht:

<a id="r-13-13"></a>

**Proposition 13.13.** Seien $H\in\mathrm{Herm}_n(\mathbb C)$ und $\delta(z)=i[H,z]$. (a) $\delta$ bildet
$\mathrm{Herm}_n$ in sich ab und ist eine Derivation von $\circ$; $e^{t\delta}(z)=e^{itH}ze^{-itH}$ ist ein
Jordan-Automorphismus. (b) $[L(x),L(y)]z=\frac14[[x,y],z]$ für $x,y,z\in\mathrm{Herm}_n$. (c)
$\delta=4\sum_k[L(x_k),L(y_k)]$ für gewisse $x_k,y_k\in\mathrm{Herm}_n$: Die Hamiltonsche Dynamik wird allein
vom Jordan-Produkt erzeugt.

*Beweis.* (a) $(i[H,z])^*=-i[z,H]=i[H,z]$; $[H,\cdot]$ ist eine Derivation des assoziativen Produkts, also auch von
$\circ$. $z(t)=e^{itH}ze^{-itH}$ erfüllt $z'=i[H,z]$ und $z(0)=z$ (Proposition 8.2(3) mit $\mathcal M=M_n$), also
ist $z(t)=e^{t\delta}(z)$ nach der Eindeutigkeit in Lemma 8.1; Konjugation mit dem Unitären $e^{itH}$ erhält Produkte
und Adjungierte, also $\circ$ und $\mathrm{Herm}_n$.
(b) $4x\circ(y\circ z)=xyz+xzy+yzx+zyx$; man subtrahiere denselben Ausdruck mit vertauschten $x,y$. (c) Nach
(b) genügt es, $iH-\frac{\mathrm{Tr}\,H}n i1$ (das dieselben Kommutatoren hat wie $iH$) als Summe von $[x_k,y_k]$
zu schreiben. Für $i\ne j$ und $a=E_{ij}+E_{ji}$, $b=i(E_{ij}-E_{ji})$ gilt $[a,E_{ii}]=E_{ji}-E_{ij}$,
$[b,E_{ii}]=-i(E_{ij}+E_{ji})$ und $[a,b]=2i(E_{jj}-E_{ii})$; diese spannen über $\mathbb R$ die spurfreien
antihermiteschen Matrizen auf. $\square$

*Rechnung:* (a) und (b) für zufällige $3\times3$-Matrizen — [Tutorial, Abschnitt 13.5](../../tutorial/de/13_jordan_algebras.md#s-13-5).

Was die Jordan-Algebra nicht enthält, ist die **Zuordnung** $H\mapsto i[H,\cdot]$ zwischen einer Observablen und der
Dynamik, die sie erzeugt. In $\mathrm{Herm}_n(\mathbb R)$ etwa stammen die Derivationen
$[L(x),L(y)]=\frac14[[x,y],\cdot]$ von reellen antisymmetrischen $[x,y]$, die nicht $i$ mal eine reelle symmetrische
Matrix sind (Aufgabe 4). Die Korrespondenz ist durch $\circ$ nicht bestimmt: Auf $\mathrm{Herm}_n(\mathbb C)$
schicken sowohl $H\mapsto i[H,\cdot]$ als auch $H\mapsto-i[H,\cdot]$ jede Observable auf eine Derivation von $\circ$
(Proposition 13.13(a)), und die Transposition, ein Jordan-Automorphismus, vertauscht sie, da
$T(i[H,T(z)])=-i[H^T,z]$. Die zweite ist $i$ mal der Kommutator der entgegengesetzten Algebra
($a\cdot_{\rm op}b=ba$): $-i[H,z]=i(H\cdot_{\rm op}z-z\cdot_{\rm op}H)$ – die beiden Vorzeichen aus Proposition 13.12.
(Weiterführende Literatur zu solchen „dynamischen Korrespondenzen“, hier nicht verwendet: E. M. Alfsen, F. W. Shultz,
*Geometry of State Spaces of Operator Algebras*, Birkhäuser, 2003.)

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Für ein spinloses Teilchen ist die
> Transposition (komplexe Konjugation in einer reellen Basis) die **Zeitumkehr**. Sie ist eine Symmetrie der Observablen eines einzelnen Systems, aber sie ist antiunitär,
> und auf nur ein Teilsystem angewandt ist sie keine physikalische Operation: Proposition 13.12 zeigt, dass sie den
> Bell-Zustand in einen Operator mit einem negativen Eigenwert verwandelt. Die Quantenmechanik verwendet komplexe
> $*$-Algebren, weil zusammengesetzte Systeme den Kommutator mit seinem Vorzeichen brauchen, nicht nur das
> symmetrisierte Produkt.

<a id="s-13-6"></a>

## 13.6 Zurück zu den Spurungleichungen

Proposition 11.9 (Exp. 3) besagt $\mathrm{Tr}|YXY|\le\mathrm{Tr}(Y|X|Y)$ für $X=X^*$ und $Y\ge0$ (hier
sind $X,Y$ Operatoren aus $\mathcal M$, keine Pauli-Matrizen). Ihr Beweis verwendet $Y\ge0$ nirgends: $YX_\pm Y=T^*T$
mit $T=X_\pm^{1/2}Y$ braucht nur $Y=Y^*$, und Lemma 11.8 wird auf $YX_+Y$ und $YX_-Y$ angewandt. Sie gilt also für
alle $Y=Y^*$ und lautet in Jordan-Sprache
$$
\mathrm{tr}|P(y)x|\le\mathrm{tr}\,P(y)|x|,
$$
was in jeder euklidischen Jordan-Algebra sinnvoll ist, sobald eine Spur festgelegt ist. Wir beweisen sie in $V_m$ mit einem Argument, das nur drei
Eigenschaften verwendet.

<a id="r-13-14"></a>

**Proposition 13.14 (Exp. 3 in einem Spinfaktor).** Für $x,y\in V_m$ gilt
$\mathrm{tr}|P(y)x|\le\mathrm{tr}\,P(y)|x|$, mit Gleichheit, falls $x\in Q$ oder $-x\in Q$.

*Beweis.* Seien $x_\pm=f_\pm(x)$ mit $f_+(t)=\max(t,0)$, $f_-(t)=\max(-t,0)$ (Satz 13.4(3)): $x_\pm\in Q$,
$x=x_+-x_-$ und $|x|=x_++x_-$. Man setze $A=P(y)x_+$ und $B=P(y)x_-$; sie liegen nach Satz 13.4(5) in $Q$, und
$P(y)x=A-B$. Sei $P(y)x=\mu_1c_1+\mu_2c_2$ eine Spektralzerlegung (Satz 13.4(2); $c_1=c_w$, $c_2=c_{-w}$, also
$c_1+c_2=1$ und $c_1\circ c_2=0$), $c$ die Summe der $c_i$ mit $\mu_i>0$ und
$s=2c-1=\sum_i\mathrm{sgn}'(\mu_i)c_i$ mit $\mathrm{sgn}'(t)=1$ für $t>0$ und $-1$ sonst. Dann ist
$s\circ P(y)x=\sum_i|\mu_i|c_i=|P(y)x|$ (nach Satz 13.4(2,3) hängt das für $\mu_1=\mu_2$ nicht von der Wahl von $w$
ab), und $c,1-c\in Q$ (Idempotente sind Quadrate). Nach der leichten Hälfte der Selbstdualität (Satz 13.4(4):
$\mathrm{tr}(a\circ b)\ge0$ für $a,b\in Q$), angewandt auf $\mathrm{tr}(A\circ(1-c))$ und
$\mathrm{tr}(B\circ c)$, gilt
$$
\mathrm{tr}|P(y)x|=\mathrm{tr}(A\circ c)-\mathrm{tr}(A\circ(1-c))-\mathrm{tr}(B\circ c)+\mathrm{tr}(B\circ(1-c))
\le\mathrm{tr}\,A+\mathrm{tr}\,B=\mathrm{tr}\,P(y)|x|.
$$
Ist $x\in Q$, so ist $x_-=0$, $B=0$ und $|P(y)x|=A$, also sind beide Seiten gleich $\mathrm{tr}\,A$; ebenso für
$-x\in Q$. $\square$

*Rechnung:* Proposition 13.14 für $m=2,4,7$ und über $\phi$ die Ungleichung aus Proposition 11.9 für $Y=Y^*$ — [Tutorial, Abschnitt 13.6](../../tutorial/de/13_jordan_algebras.md#s-13-6).

Der Beweis verwendet nur: (i) $P(y)Q\subseteq Q$ für jedes $y$; (ii) den Spektralsatz: Jedes $z$ ist
$\sum_i\mu_ic_i$ mit orthogonalen Idempotenten $c_i$ der Summe $1$, und $f(z)=\sum_if(\mu_i)c_i$ (sodass $x_\pm$ und
die Idempotente Quadrate sind); (iii) $\mathrm{tr}$ ist linear, $\mathrm{tr}(a\circ1)=\mathrm{tr}\,a$,
und $\mathrm{tr}(a\circ b)\ge0$ für $a,b\in Q$. Alle drei gelten in jeder euklidischen Jordan-Algebra $V$ mit Eins
$e$ für die Linearform $\tau(x)=\langle x,e\rangle$, für die $\tau(a\circ b)=\langle b,a\circ e\rangle=\langle a,b\rangle$
gilt: (ii) ist Satz F.3, (i) ist Satz F.6(4), und $\langle a,b\rangle\ge0$ für $a,b\in Q$ ist Satz F.5(1). Korollar
F.7 führt das Argument aus: $\tau(|P(y)x|)\le\tau(P(y)|x|)$ in jeder euklidischen Jordan-Algebra, mit Gleichheit, falls
$x\in Q$ oder $-x\in Q$. Auf $M_{sa}$ mit $\langle x,y\rangle=\mathrm{Tr}(xy)$ ist $\tau=\mathrm{Tr}$ und
$P(y)x=yxy$; die Idempotente von $(M_{sa},\circ)$ sind die Projektionen ($c\circ c=c^2$), und für Projektionen gilt
$c\circ c'=0$ genau dann, wenn $cc'=0$ (Lemma F.14(2)); nach der Eindeutigkeit in Satz F.3 ist seine Zerlegung von $x$
also die Spektralzerlegung in Eigenprojektionen, und $|x|$ ist der Matrix-Absolutbetrag. Damit ist Proposition 11.9
für alle $Y=Y^*$ erneut bewiesen. Auf $V_m$ mit $\langle x,y\rangle=\mathrm{tr}(x\circ y)$,
$\tau=\mathrm{tr}$, ist die Zerlegung aus Satz F.3 die aus Satz 13.4(2) (mit $x=a1$, falls $u=0$), und das ist
Proposition 13.14. In einer einfachen euklidischen Jordan-Algebra $J$ vom Rang $r$ erfüllt die Spur
$\mathrm{tr}\,x=\frac r{\dim J}\mathrm{Tr}\,L(x)$ aus Korollar E.28 die Gleichung
$\mathrm{tr}(\sum_i\mu_ic_i)=\sum_i\mu_i$ für jeden Jordan-Rahmen $(c_i)$ (Satz E.7), und
$\langle x,y\rangle=\mathrm{tr}(x\circ y)$ ist ein assoziatives Skalarprodukt; für dieses ist
$\tau=\mathrm{tr}$, also gilt die Ungleichung für diese Spur. Das trifft auf $\mathrm{Herm}_r(\mathbb H)$
und die Albert-Algebra zu, die einfache euklidische Jordan-Algebren sind (Satz E.24). Für $r\ge3$ sind sie für kein
$\mathcal M=\bigoplus_cM_{k_c}$ zu $M_{sa}$ Jordan-isomorph: $M_{sa}=\bigoplus_c\mathrm{Herm}_{k_c}(\mathbb C)$ ist
eine direkte Summe von Idealen, also nur bei einem einzigen Block einfach, und $\mathrm{Herm}_k(\mathbb C)$
($\cong\mathbb R$ für $k=1$, $\cong V_3$ für $k=2$) ist nach der Eindeutigkeit in Satz 13.7 weder zu
$\mathrm{Herm}_r(\mathbb H)$, $r\ge3$, noch zu $\mathrm{Herm}_3(\mathbb O)$ isomorph. In $V_m$ folgt die
Ungleichung auch aus Proposition 11.9 (oben auf $Y=Y^*$ erweitert) über $\phi$, da $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$
und $\phi(|x|)=|\phi(x)|$.

<a id="s-13-7"></a>

## 13.7 Womit man vergleichen sollte

*In der Bibliothek:* Vergleich mit anderer Software (SageMath, Magma, SOCP-Löser) und was die Bibliothek hinzufügt — [Tutorial, Abschnitt 13.7](../../tutorial/de/13_jordan_algebras.md#s-13-7).

<a id="x-13"></a>

## Aufgaben

1. (**Homogenität des Lorentz-Kegels.**) Sei $y\in V_m$ mit $\det y\ne0$. Man zeige, dass $y^{-1}=(a,-u)/\det y$ die
   Gleichung $y\circ y^{-1}=1$ erfüllt, dass $P(y)P(y^{-1})=\mathrm{id}$, dass $P(y)$ das Innere von $Q$ auf sich
   abbildet und dass $P(x^{1/2})1=x$ für $x$ im Inneren. Man folgere, dass das Innere von $Q$ homogen ist.
2. (**Isotrope Zustände.**) Für $p\in[0,1]$ sei $\rho_p=p\,\omega+(1-p)\frac14$ auf $\mathbb C^2\otimes\mathbb C^2$
   (Proposition 13.12). Man berechne das Spektrum von $\theta(\rho_p)$ und zeige, dass $\rho_p$ für $p>\frac13$
   verschränkt ist.
3. (**Produkte von Randpunkten.**) Seien $|w|=|w'|=1$ in $\mathbb R^m$. Man zeige, dass $(1,w)\circ(1,w')\in Q$ genau
   dann gilt, wenn $w'=w$ oder $w'=-w$: Das Jordan-Produkt zweier Randstrahlen des Lorentz-Kegels liegt fast nie im
   Kegel.
4. (**Reelle Quantenmechanik.**) Man zeige, dass für $x,y\in\mathrm{Herm}_n(\mathbb R)$ die Derivation
   $[L(x),L(y)]$ von $\mathrm{Herm}_n(\mathbb R)$ gleich $z\mapsto[K,z]$ mit reellem antisymmetrischem $K$ ist
   und dass für $n\ge2$ eine solche Abbildung nie von der Form $z\mapsto i[H,z]$ mit
   $H\in\mathrm{Herm}_n(\mathbb R)$ ist, es sei denn, sie ist $0$.

*Rechnung:* Aufgaben 1–3 numerisch — [Tutorial, Aufgaben](../../tutorial/de/13_jordan_algebras.md#x-13).

<details>
<summary>Lösungen</summary>

1. $(a,u)\circ(a,-u)=(a^2-|u|^2,0)$, was $y\circ y^{-1}=1$ ergibt. $\phi(y)$ und $\phi(y^{-1})$ kommutieren (beide
   liegen in $\mathrm{span}\{1,u\cdot s\}$), also ist $\phi(y)\phi(y^{-1})=\phi(y)\circ\phi(y^{-1})=1$, und nach
   Satz 13.4(5) ist $\phi(P(y)P(y^{-1})x)=\phi(x)$; $\phi$ ist injektiv. $P(y)$ bildet $Q$ in $Q$ ab, ebenso seine
   Inverse $P(y^{-1})$, also ist $P(y)Q=Q$; als linearer Homöomorphismus bildet es das Innere auf das Innere ab.
   Schließlich ist $P(y)1=2y\circ y-y^2=y^2$, also $P(x^{1/2})1=x$, und $P(x^{1/2})$ ist invertierbar, wenn $x$ im
   Inneren liegt ($\det x^{1/2}=\sqrt{\lambda_+\lambda_-}>0$; das Innere von $Q$ ist $\{a>|u|\}=\{\lambda_->0\}$, da
   $a-|u|$ stetig ist und $(a-\varepsilon,u)\notin Q$ für $a=|u|$). Je zwei innere Punkte werden durch
   $P(x'^{1/2})P(x^{1/2})^{-1}$ ineinander überführt.

2. $\theta(\rho_p)=\frac p2F+\frac{1-p}4$ nach Proposition 13.12. $F$ hat die Eigenwerte $1$ (dreifach) und $-1$
   (einfach), also besteht das Spektrum aus $\frac{1+p}4$ (dreifach) und $\frac{1-3p}4$, das genau dann negativ ist,
   wenn $p>\frac13$; dann ist $\rho_p$ nach dem PPT-Kriterium nicht separabel.

3. $(1,w)\circ(1,w')=(1+t,\,w+w')$ mit $t=w\cdot w'\in[-1,1]$ und $|w+w'|=\sqrt{2+2t}$. Es liegt genau dann in $Q$,
   wenn $1+t\ge\sqrt{2(1+t)}$, d. h. genau dann, wenn $t=-1$ (das Produkt ist $0$) oder $1+t\ge2$, d. h. $t=1$.

4. Nach Proposition 13.13(b) ist $[L(x),L(y)]z=[K,z]$ mit $K=\frac14[x,y]$, und $K^T=\frac14[y^T,x^T]=-K$ ist reell.
   Für reelle $H$ und $z$ ist $i[H,z]$ rein imaginär, während $[K,z]$ reell ist; sie stimmen also nur überein, wenn
   beide verschwinden. In der reellen Quantenmechanik sind die Erzeugenden der Dynamik, $K$, also keine Observablen.

</details>

**Wie es weitergeht.** Dies ist das letzte Kapitel des Kurses; die Anhänge E und F vervollständigen seine Beweise. Als
weiterführende Literatur ist die Standardreferenz für die endlichdimensionale Theorie [FK94]; K. McCrimmon, *A Taste of Jordan Algebras* (Springer 2004), behandelt die Algebra und die
Geschichte, H. Hanche-Olsen und E. Størmer, *Jordan Operator Algebras* (Pitman 1984), sowie Alfsen–Shultz (2003) die
unendlichdimensionale Theorie und ihre Beziehung zu C\*-Algebren. Zur Rolle zusammengesetzter Systeme bei
Rekonstruktionen der Quantentheorie aus operationellen Axiomen siehe die oben zitierten Arbeiten von Barnum und Wilce.
Zu Spurungleichungen über Kapitel 11 hinaus: E. Carlen, *Trace inequalities and quantum entropy* (Contemp. Math. 529,
2010). Die quadratische Darstellung und der Lorentz-Kegel tauchen in der konischen Optimierung wieder auf:
Alizadeh–Goldfarb (2003) und die dort besprochenen Jordan-algebraischen Innere-Punkte-Verfahren.
