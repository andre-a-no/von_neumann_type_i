[English](../C_eigenvalue_densities.md) · [Español](../es/C_eigenvalue_densities.md) · [Français](../fr/C_eigenvalue_densities.md) · **Deutsch** · [中文](../zh/C_eigenvalue_densities.md) · [日本語](../ja/C_eigenvalue_densities.md) · [Русский](../ru/C_eigenvalue_densities.md)

# Anhang C. Eigenwertdichten und die Pagesche Formel

Kapitel 5 verwendet die gemeinsame Dichte der Eigenwerte der Gaußschen Ensembles, und Kapitel 7 verwendet
die Pagesche Formel für die mittlere Verschränkungsentropie eines zufälligen reinen Zustands. Dieser Anhang beweist
beides. Der Weg führt über die Weylsche Integrationsformel (die Jacobi-Determinante der „Diagonalisierung“), die
Gaußschen Ensembles, die Einpunktfunktion einer determinantiellen Dichte mittels Orthogonalpolynomen, das
Wishart-Laguerre-Ensemble und das induzierte Maß auf reduzierten Zuständen und schließlich zu einer exakten
Rechnung mit Laguerre-Polynomen.

**Sie lernen**, warum sich die Eigenwerte einer invarianten Zufallsmatrix wie
$`\prod_{i<j}|\lambda_i-\lambda_j|^\beta`$ abstoßen, warum die Eigenvektoren dann Haar-verteilt und unabhängig
von den Eigenwerten sind, wie sich eine Determinante orthonormaler Polynome ausintegrieren lässt und wie
$`\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}`$ aus einer Handvoll Integralen folgt.

Neben den Voraussetzungen des Kurses (Maßtheorie, die Transformationsformel im
$`\mathbb R^N`$ für injektive $`C^1`$-Abbildungen mit invertierbarer Ableitung, der Satz über die Umkehrfunktion und der Satz über
implizite Funktionen) verwenden wir aus Kapitel 5 Lemma 5.1 (Eindeutigkeit und Rechtsinvarianz des Haar-Maßes),
Lemma 5.2 und Satz 5.3 (Existenz des Haar-Maßes auf $`U(n)`$ durch Gaußsche QR-Zerlegung), aus Anhang B
Satz B.17 (dasselbe für $`O(n)`$), aus Anhang A Proposition A.16(a) (unitäre Invarianz komplexer
Gaußscher Vektoren) und Satz A.20(b) (die Verteilung von $`|z|^2`$) sowie aus Kapitel 7 die
Darstellung $`\rho_A=GG^*/\mathrm{Tr}\,GG^*`$ aus dem Beweis von Satz 7.11.

**Bezeichnungen.** $`\beta\in\{1,2\}`$ und entsprechend $`\mathbb F=\mathbb R`$ oder $`\mathbb C`$. $`\mathcal H`$ ist der
Raum der reellen symmetrischen ($`\beta=1`$) bzw. komplexen hermiteschen ($`\beta=2`$) $`n\times n`$-Matrizen mit dem
Lebesgue-Maß $`dH`$ in den Koordinaten $`h_{ii}`$, $`\mathrm{Re}\,h_{ij}`$ und (für $`\beta=2`$)
$`\mathrm{Im}\,h_{ij}`$, $`i<j`$; seine Dimension ist $`n+\beta n(n-1)/2`$. $`K=O(n)`$ bzw. $`U(n)`$, $`T`$ ist die
Gruppe der Diagonalelemente von $`K`$ (Vorzeichen bzw. Phasen), $`\mathfrak k=\{X\in M_n(\mathbb F):X^*=-X\}`$,
$`\mathfrak t`$ sein Diagonalanteil und $`\mathfrak m`$ sein Nichtdiagonalanteil (Diagonale null), mit den
Koordinaten $`\mathrm{Re}\,X_{ij}`$, $`\mathrm{Im}\,X_{ij}`$, $`i<j`$. Für $`\lambda\in\mathbb R^n`$ sei
$`\Lambda=\mathrm{diag}(\lambda)`$ und $`\Delta(\lambda)=\prod_{i<j}(\lambda_j-\lambda_i)`$; $`C=\{\lambda:\lambda_1<\dots<\lambda_n\}`$,
auf dem $`\Delta>0`$ gilt. Für $`H\in\mathcal H`$ ist $`\lambda(H)`$ der Vektor seiner Eigenwerte in aufsteigender
Reihenfolge. $`\|X\|_2^2=\mathrm{Tr}\,X^*X`$. Für $`V\in K`$ ist die Abbildung $`H\mapsto VHV^*`$ linear auf $`\mathcal H`$
und erhält $`\|\cdot\|_2`$; eine lineare Isometrie hat $`|\det|=1`$, sie **erhält also $`dH`$**.

*Rechnung:* Die numerischen Prüfungen dieses Anhangs — [Tutorial, Anhang C](../../tutorial/de/C_eigenvalue_densities.md).

<a id="s-C-1"></a>

## C.1 Einfache Spektren und Eigenrahmen

<a id="r-C-1"></a>

**Lemma C.1.** (a) Die Menge $`\Omega\subseteq\mathcal H`$ der Matrizen mit $`n`$ verschiedenen Eigenwerten ist
offen, und $`\mathcal H\setminus\Omega`$ ist eine Lebesgue-Nullmenge. (b) Auf $`\Omega`$ ist die Abbildung $`H\mapsto\lambda(H)`$
$`C^\infty`$, und die Spektralprojektionen $`P_j(H)=\prod_{k\neq j}\frac{H-\lambda_k(H)}{\lambda_j(H)-\lambda_k(H)}`$ sind stetig.

*Beweis.* (a) Sei $`A_{ij}=\lambda_j^{i-1}`$ mit $`\lambda=\lambda(H)`$. Dann ist $`\det A=\Delta(\lambda)`$ (Vandermonde) und
$`(AA^T)_{ij}=\sum_k\lambda_k^{i+j-2}=\mathrm{Tr}\,H^{i+j-2}`$, also ist
$`\mathrm{disc}(H)=\det\big(\mathrm{Tr}\,H^{i+j-2}\big)_{i,j=1}^n=\Delta(\lambda(H))^2`$ ein Polynom
in den Koordinaten von $`H`$, ungleich null in $`\mathrm{diag}(1,\dots,n)`$, und $`\Omega=\{\mathrm{disc}\,\neq0\}`$ ist offen.
Die Nullstellenmenge eines vom Nullpolynom verschiedenen Polynoms $`p`$ auf $`\mathbb R^N`$ ist eine Nullmenge, per Induktion nach $`N`$: Für $`N=1`$ ist sie
endlich; andernfalls schreibe man $`p=\sum_jp_j(x')x_N^j`$ mit einem $`p_{j_0}\neq0`$; für $`x'`$ außerhalb der Nullstellenmenge von
$`p_{j_0}`$ (nach Induktion eine Nullmenge) hat $`x_N\mapsto p(x',x_N)`$ endlich viele Nullstellen, und der Satz von Tonelli, angewandt
auf die Indikatorfunktion der Nullstellenmenge, liefert das Maß $`0`$.
(b) $`\chi(H,t)=\det(t-H)=\prod_i(t-\lambda_i)`$ ist ein reelles Polynom in $`(H,t)`$. In $`H_0\in\Omega`$ ist
$`\partial_t\chi(H_0,\lambda_j(H_0))=\prod_{k\neq j}(\lambda_j-\lambda_k)\neq0`$, also liefert der Satz über implizite Funktionen
$`C^\infty`$-Funktionen $`r_j`$ nahe $`H_0`$ mit $`\chi(H,r_j(H))=0`$ und $`r_j(H_0)=\lambda_j(H_0)`$. Aus Stetigkeitsgründen gilt
$`r_1<\dots<r_n`$ nahe $`H_0`$: Das sind $`n`$ verschiedene Nullstellen eines Polynoms vom Grad $`n`$, also alle, der Reihe
nach, und $`r_j=\lambda_j`$. Ist $`H=\sum_k\lambda_kQ_k`$ die Spektralzerlegung, so liefert die Lagrange-Interpolation
$`P_j(H)=Q_j`$; die Stetigkeit folgt aus der von $`\lambda`$. $`\square`$

Für $`H\in\Omega`$ ist ein **Eigenrahmen** ein $`V\in K`$ mit $`H=V\Lambda V^*`$, $`\Lambda=\mathrm{diag}\,\lambda(H)`$. Sind $`V,V'`$
Eigenrahmen, so vertauscht $`M=V'^*V`$ mit $`\Lambda`$, also ist $`(\lambda_i-\lambda_j)M_{ij}=0`$ und $`M\in T`$; umgekehrt ist $`VD`$
für jedes $`D\in T`$ ein Eigenrahmen. Einen **Borelschen Eigenrahmen** $`F(H)`$ erhält man, indem man als $`j`$-te Spalte
$`P_j(H)e_m/\|P_j(H)e_m\|`$ nimmt, wobei $`m`$ der kleinste Index mit $`P_j(H)e_m\neq0`$ ist: Er ist stetig auf jeder der
Borelmengen, auf denen diese Indizes fest sind, seine Spalten sind orthonormale Eigenvektoren, und er ist reell, wenn $`H`$ es ist.

<a id="r-C-2"></a>

**Lemma C.2 (Haar-Maße).** $`U(n)`$, $`O(n)`$ und $`T`$ besitzen Haar-Maße, und $`\int_Kf(VW)\,dV=\int_Kf\,dV`$
für alle $`W\in K`$ und alle beschränkten Borel-Funktionen $`f`$.

*Beweis.* $`U(n)`$: Satz 5.3. $`O(n)`$: Satz B.17 (die reelle Gaußsche QR-Konstruktion: Eine reelle Gaußsche Matrix ist
fast sicher invertierbar, und der orthogonale Faktor ihrer QR-Zerlegung ist Haar-verteilt). $`T`$: Das Gleichmaß auf $`\{\pm1\}^n`$ bzw. das Produkt der normierten Bogenlängen auf den
Kreisen ist invariant unter Multiplikation mit einem festen Element. Die Rechtsinvarianz ist Lemma 5.1. $`\square`$

<a id="r-C-3"></a>

**Lemma C.3 (Eigenrahmen invarianter Ensembles).** Sei $`H`$ ein zufälliges Element von $`\mathcal H`$ mit
$`P(H\in\Omega)=1`$, dessen Verteilung für jedes $`V\in K`$ invariant unter $`H\mapsto VHV^*`$ ist. Sei $`D`$ gleichverteilt auf $`T`$ und
unabhängig von $`H`$, und $`\tilde U=F(H)D`$. Dann ist $`\tilde U`$ ein Eigenrahmen von $`H`$, Haar-verteilt auf
$`K`$ und unabhängig von $`\lambda(H)`$.

*Beweis.* Sei $`V\in K`$ fest. Sowohl $`F(VHV^*)`$ als auch $`VF(H)`$ sind Eigenrahmen von $`VHV^*`$, also ist
$`F(VHV^*)=VF(H)E(H)`$ mit $`E(H)=(VF(H))^*F(VHV^*)\in T`$, einer Borel-Funktion von $`H`$. Für beschränkte Borel-Funktionen
$`\varphi`$ auf $`\mathbb R^n`$ und $`\psi`$ auf $`K`$ gilt, da $`VHV^*`$ die Verteilung von $`H`$ hat und $`\lambda(VHV^*)=\lambda(H)`$,
```math
\mathbb E\big[\varphi(\lambda(H))\psi(\tilde U)\big]=\mathbb E\big[\varphi(\lambda(H))\psi(VF(H)E(H)D)\big]
=\mathbb E\big[\varphi(\lambda(H))\psi(V\tilde U)\big],
```
wobei der zweite Schritt zuerst über $`D`$ integriert (Fubini): Für festes $`H`$ ist $`E(H)D`$ gleichverteilt auf $`T`$. Sei nun $`V`$
Haar-verteilt und von allem anderen unabhängig; nach Fubini und Lemma C.2 ist
$`\mathbb E[\varphi(\lambda(H))\psi(\tilde U)]=\mathbb E\big[\varphi(\lambda(H))\int_K\psi(V\tilde U)\,dV\big] =\mathbb E[\varphi(\lambda(H))]\int_K\psi\,dV`$. Mit $`\varphi=1`$ ist $`\tilde U`$ Haar-verteilt; die Produktform ist die
Unabhängigkeit. $`\square`$

<a id="s-C-2"></a>

## C.2 Die Weylsche Integrationsformel

Für $`\varepsilon>0`$ sei $`\mathfrak m_\varepsilon=\{X\in\mathfrak m:\|X\|_2<\varepsilon\}`$ und
```math
\Psi(X,\lambda)=e^X\Lambda e^{-X},\qquad (X,\lambda)\in\mathfrak m_\varepsilon\times C .
```
Für $`X\in\mathfrak k`$ ist $`(e^X)^*=e^{-X}=(e^X)^{-1}`$, also $`e^X\in K`$.

<a id="r-C-4"></a>

**Lemma C.4 (eine Karte).** Es gibt $`\varepsilon>0`$ und eine stetige Funktion $`j:\overline{\mathfrak m_\varepsilon}\to(0,\infty)`$
mit $`j(0)=1`$, sodass gilt: (a) $`|\det D\Psi(X,\lambda)|=\Delta(\lambda)^\beta j(X)`$; (b) $`\Psi`$ ist injektiv; (c) für
$`H\in\Omega`$ ist $`H\in\Psi(\mathfrak m_\varepsilon\times C)`$ genau dann, wenn $`F(H)\in S=\{e^XD:X\in\mathfrak m_\varepsilon,D\in T\}`$, und
$`S`$ ist eine Borelmenge mit $`ST=S`$.

*Beweis.* (a) $`\Psi`$ ist $`C^\infty`$ (Potenzreihen). Für $`Y\in\mathfrak m`$ setze $`A_X(Y)=e^{-X}\frac{d}{dt}\big|_0e^{X+tY}`$.
Ableiten von $`W(t)^*W(t)=1`$ für $`W(t)=e^{X+tY}\in K`$ ergibt $`A_X(Y)\in\mathfrak k`$, und mit
$`\frac{d}{dt}W^{-1}=-W^{-1}W'W^{-1}`$ gilt für $`g=e^X`$
```math
D\Psi(X,\lambda)(Y,\mu)=g\big([A_X(Y),\Lambda]+\mathrm{diag}(\mu)\big)g^{-1}.
```
Diagonalmatrizen vertauschen mit $`\Lambda`$, also ist $`[A_X(Y),\Lambda]=[B_XY,\Lambda]`$ mit $`B_X=\pi\circ A_X:\mathfrak m\to\mathfrak m`$, wobei
$`\pi`$ die Projektion von $`\mathfrak k`$ auf $`\mathfrak m`$ längs $`\mathfrak t`$ ist. Somit ist $`D\Psi=\mathrm{Ad}_g\,\circ L_\lambda\circ(B_X\times1)`$
mit $`L_\lambda(Z,\mu)=[Z,\Lambda]+\mathrm{diag}(\mu)`$. Wegen $`([Z,\Lambda])_{ij}=Z_{ij}(\lambda_j-\lambda_i)`$ multipliziert $`L_\lambda`$ in den Koordinaten
$`\mathrm{Re}\,Z_{ij}`$ (und $`\mathrm{Im}\,Z_{ij}`$) mit $`\lambda_j-\lambda_i`$ und bildet $`\mu`$ auf die Diagonale ab:
$`|\det L_\lambda|=\Delta(\lambda)^\beta`$. $`\mathrm{Ad}_g`$ erhält $`dH`$. Setze $`j(X)=|\det B_X|`$: Sie ist stetig, und
$`B_0=1`$, weil $`\frac d{dt}e^{tY}|_0=Y`$; man wähle $`\varepsilon`$ mit $`j>0`$ auf $`\overline{\mathfrak m_\varepsilon}`$.
(b) Ist $`\Psi(X,\lambda)=\Psi(X',\lambda')=H`$, so ist $`\lambda=\lambda'=\lambda(H)`$, und $`e^X,e^{X'}`$ sind Eigenrahmen von $`H\in\Omega`$,
also $`e^X=e^{X'}D`$ mit $`D\in T`$. Sei $`G(X,t)=e^Xe^t`$ auf $`\mathfrak m\times\mathfrak t`$ (für $`\beta=1`$ ist $`\mathfrak t=\{0\}`$); $`DG(0,0)(Y,s)=Y+s`$ ist injektiv,
und da $`Y`$ außerhalb der Diagonale liegt und $`s`$ diagonal ist, gilt $`\|Y+s\|_2^2=\|Y\|_2^2+\|s\|_2^2`$, also $`\|DG(0,0)v\|=\|v\|`$ für die Produktnorm; setze $`m=1`$. Man wähle $`\delta>0`$ mit $`\|DG(z)-DG(0)\|\le m/2`$ für $`\|z\|<\delta`$; für $`u,v`$ in dieser
Kugel liefert $`G(u)-G(v)=\int_0^1DG(v+s(u-v))(u-v)\,ds`$ die Abschätzung $`\|G(u)-G(v)\|\ge\frac m2\|u-v\|`$, also ist $`G`$ dort injektiv.
($`\delta`$ hängt nur von $`n,\beta`$ ab; $`\varepsilon`$ wird danach gewählt.)
Nun ist $`\|D-1\|_2=\|e^{X'}(D-1)\|_2=\|e^X-e^{X'}\|_2\le\|e^X-1\|_2+\|e^{X'}-1\|_2\le2(e^\varepsilon-1)`$ nach der Potenzreihe und
$`\|AB\|_2\le\|A\|_2\|B\|_2`$. Für $`\beta=1`$ ist $`D=\mathrm{diag}(\pm1)`$, und ein Eintrag $`-1`$ ergäbe
$`\|D-1\|_2\ge2`$; für $`\varepsilon<\log2`$ ist also $`D=1=e^t`$ mit $`t=0`$. Für $`\beta=2`$ ist $`D=e^t`$ mit $`t=i\mathrm{diag}(\theta)`$,
$`\theta_j\in(-\pi,\pi]`$, und $`|\theta_j|\le\frac\pi2|e^{i\theta_j}-1|`$ (da $`|e^{i\theta}-1|=2|\sin\frac\theta2|`$ und $`\sin x\ge2x/\pi`$ auf $`[0,\pi/2]`$), also
$`\|t\|_2\le\frac\pi2\|D-1\|_2\le\pi(e^\varepsilon-1)`$. Für kleines $`\varepsilon`$ ($`\varepsilon^2+\pi^2(e^\varepsilon-1)^2<\delta^2`$) liegen $`(X,0)`$ und $`(X',t)`$ in der $`\delta`$-Kugel, und
$`G(X,0)=e^X=e^{X'}D=G(X',t)`$, also $`X=X'`$.
(c) Ist $`H=\Psi(X,\lambda')`$, so ist $`e^X`$ ein Eigenrahmen, also $`F(H)\in e^XT\subseteq S`$. Ist $`F(H)=e^XD`$, so ist
$`H=F(H)\Lambda F(H)^*=\Psi(X,\lambda(H))`$. $`S=G(\mathfrak m_\varepsilon\times\mathfrak t)`$ (bzw. $`e^{\mathfrak m_\varepsilon}T`$, $`T`$ endlich) ist ein
stetiges Bild einer $`\sigma`$-kompakten Menge, also $`\sigma`$-kompakt und Borelsch; $`ST=S`$ ist klar. $`\square`$

*Rechnung:* Die Jacobi-Determinante von $`\Psi`$ für $`n=2,3`$ und $`\beta=1,2`$ — [Tutorial, Abschnitt C.2](../../tutorial/de/C_eigenvalue_densities.md#s-C-2).

<a id="r-C-5"></a>

**Satz C.5 (Weylsche Integrationsformel).** Es gibt eine Konstante $`c_{n,\beta}>0`$, sodass für jede Borel-Funktion
$`f:\mathcal H\to[0,\infty]`$ gilt
```math
\int_{\mathcal H}f(H)\,dH=c_{n,\beta}\int_C\Delta(\lambda)^\beta\Big(\int_Kf(V\Lambda V^*)\,dV\Big)d\lambda .
```
Äquivalent: $`\Phi(V,\lambda)=V\Lambda V^*`$ bildet das Maß $`c_{n,\beta}\,dV\otimes\Delta(\lambda)^\beta d\lambda`$ auf $`K\times C`$
auf das Lebesgue-Maß auf $`\mathcal H`$ ab. Der Faktor $`\Delta(\lambda)^\beta`$ ist die Jacobi-Determinante der Karte $`\Psi`$ in $`X=0`$
(Lemma C.4(a)); die Konstante $`c_{n,\beta}`$ spiegelt die Normierung des Haar-Maßes wider und wird im Folgenden nie
explizit benötigt, da jede aus der Formel gewonnene Dichte durch ihre Gesamtmasse normiert wird.

*Beweis.* *Referenzensemble.* $`H_0`$ habe die Dichte $`p_0(H)=Z_0^{-1}e^{-\mathrm{Tr}\,H^2/2}`$; wegen
$`\mathrm{Tr}\,H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2`$ ist $`Z_0`$ ein endliches Produkt Gaußscher Integrale. $`p_0`$ und $`dH`$
sind invariant unter $`H\mapsto VHV^*`$, also auch die Verteilung von $`H_0`$; $`P(H_0\in\Omega)=1`$ nach Lemma C.1. Nach Lemma C.3
ist $`\tilde U`$ Haar-verteilt, unabhängig von $`\lambda=\lambda(H_0)`$, und $`H_0=\tilde U\Lambda\tilde U^*`$.
*Die Verteilung von $`\lambda`$.* Sei $`A\subseteq C`$ eine Borelmenge. Nach Lemma C.4(c) und $`ST=S`$ ist
$`\{H_0\in\Psi(\mathfrak m_\varepsilon\times A)\}=\{\lambda\in A,\ \tilde U\in S\}`$ bis auf ein Nullereignis, dessen Wahrscheinlichkeit
$`P(\lambda\in A)\,h`$ mit $`h=\int_K1_S\,dV`$ ist. Andererseits ist $`\Psi`$ injektiv und $`C^1`$ mit nichtverschwindender Jacobi-Determinante auf der
offenen Menge $`\mathfrak m_\varepsilon\times C`$, und $`p_0(\Psi(X,\lambda))=p_0(\Lambda)`$; die Transformationsformel und Lemma C.4(a) liefern
```math
P\big(H_0\in\Psi(\mathfrak m_\varepsilon\times A)\big)=\int_{\mathfrak m_\varepsilon}\!\int_Ap_0(\Lambda)\Delta(\lambda)^\beta j(X)\,d\lambda\,dX
=\kappa\int_Ap_0(\Lambda)\Delta(\lambda)^\beta\,d\lambda,\qquad \kappa=\int_{\mathfrak m_\varepsilon}j>0 .
```
Mit $`A=C`$ ist $`h=\kappa\int_Cp_0\Delta^\beta>0`$. Also hat $`\lambda`$ die Dichte $`c\,p_0(\Lambda)\Delta(\lambda)^\beta`$ auf $`C`$, $`c=\kappa/h`$.
*Allgemeines $`f`$.* Da $`p_0>0`$, ist $`\int f\,dH=\mathbb E\big[(f/p_0)(\tilde U\Lambda\tilde U^*)\big]`$. Die gemeinsame Verteilung von $`(\tilde U,\lambda)`$ ist
$`dV\otimes c\,p_0(\Lambda)\Delta^\beta d\lambda`$, und $`p_0(V\Lambda V^*)=p_0(\Lambda)`$; der Satz von Tonelli liefert die Formel mit $`c_{n,\beta}=c`$.
Die Konstante ist durch die Formel selbst bestimmt (man nehme $`f=p_0`$). $`\square`$

<a id="r-C-6"></a>

**Korollar C.6 (invariante Dichten).** $`H`$ habe eine Dichte $`p`$ auf $`\mathcal H`$ mit $`p(VHV^*)=p(H)`$ für alle
$`V\in K`$. Dann hat $`\lambda(H)`$ die Dichte $`c_{n,\beta}\,p(\Lambda)\Delta(\lambda)^\beta`$ auf $`C`$; der Eigenrahmen $`\tilde U`$ aus Lemma C.3
ist Haar-verteilt und unabhängig von $`\lambda(H)`$; und für symmetrische Borel-Funktionen $`g\ge0`$ auf $`\mathbb R^n`$ gilt
```math
\mathbb E\,g(\lambda(H))=\frac{c_{n,\beta}}{n!}\int_{\mathbb R^n}g(x)\,p(\mathrm{diag}\,x)\,|\Delta(x)|^\beta\,dx .
```

*Beweis.* $`P(H\in\Omega)=1`$ nach Lemma C.1. Für eine Borelmenge $`A\subseteq C`$ wende man Satz C.5 auf $`f=1_A(\lambda(\cdot))\,p`$ an, wofür
$`f(V\Lambda V^*)=1_A(\lambda)p(\Lambda)`$ gilt. Die Verteilung von $`H`$ ist invariant, also ist Lemma C.3 anwendbar. Bis auf die Nullmenge, auf der zwei
Koordinaten übereinstimmen, ist $`\mathbb R^n`$ die disjunkte Vereinigung der $`n!`$ Bilder von $`C`$ unter Koordinatenpermutationen, und
$`g`$, $`|\Delta|`$ und $`x\mapsto p(\mathrm{diag}\,x)`$ sind permutationsinvariant (eine Permutationsmatrix liegt in $`K`$). $`\square`$

<a id="s-C-3"></a>

## C.3 Die Gaußschen Ensembles

<a id="r-C-7"></a>

**Satz C.7 (GOE und GUE).** Sei $`H=(G+G^T)/2`$ mit reellem $`G`$ mit unabhängigen $`N(0,1)`$-Einträgen ($`\beta=1`$, GOE) oder
$`H=(G+G^*)/\sqrt2`$ mit einer Ginibre-Matrix $`G`$ ($`\beta=2`$, GUE), wie in Abschnitt 5.3. Dann hat $`H`$ die Dichte
$`(2\pi)^{-n/2}\pi^{-\beta n(n-1)/4}e^{-\mathrm{Tr}\,H^2/2}`$ auf $`\mathcal H`$, seine Eigenwerte haben die Dichte
```math
Z_{n,\beta}^{-1}\,\Delta(\lambda)^\beta\,e^{-\sum_i\lambda_i^2/2}\quad\text{on }C,\qquad Z_{n,\beta}=\int_C\Delta^\beta e^{-|\lambda|^2/2}\,d\lambda\in(0,\infty),
```
und der Eigenrahmen $`\tilde U`$ ist Haar-verteilt auf $`K`$ und unabhängig von den Eigenwerten.

*Beweis.* GOE: $`h_{ii}=G_{ii}\sim N(0,1)`$ und $`h_{ij}=(G_{ij}+G_{ji})/2\sim N(0,\frac12)`$ für $`i<j`$. GUE: $`h_{ii}=\sqrt2\mathrm{Re}\,G_{ii}\sim N(0,1)`$,
und $`\mathrm{Re}\,h_{ij}=(\mathrm{Re}\,G_{ij}+\mathrm{Re}\,G_{ji})/\sqrt2`$ sowie $`\mathrm{Im}\,h_{ij}=(\mathrm{Im}\,G_{ij}-\mathrm{Im}\,G_{ji})/\sqrt2`$ sind $`N(0,\frac12)`$-verteilt.
Verschiedene Koordinaten sind Funktionen disjunkter Familien unabhängiger Variablen, also unabhängig, und die
Dichte ist $`\prod_i(2\pi)^{-1/2}e^{-h_{ii}^2/2}\prod_{i<j}\pi^{-\beta/2}e^{-|h_{ij}|^2}`$, also die angegebene, weil
$`\mathrm{Tr}\,H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2`$. Sie ist invariant, und Korollar C.6 ist anwendbar; die Konstante ist durch die
Normierung festgelegt. $`Z_{n,\beta}<\infty`$, weil $`|\lambda_j-\lambda_i|\le(1+|\lambda_i|)(1+|\lambda_j|)`$ und daher
$`|\Delta(\lambda)|\le\prod_i(1+|\lambda_i|)^{n-1}`$. $`\square`$

<a id="r-C-8"></a>

**Korollar C.8 (Niveaudichte).** Für das GOE oder GUE der Größe $`n`$ sei $`q(x)=(n!Z_{n,\beta})^{-1}|\Delta(x)|^\beta e^{-|x|^2/2}`$ auf
$`\mathbb R^n`$ und $`\rho(x)=n\int_{\mathbb R^{n-1}}q(x,y)\,dy`$. Dann gilt $`\mathbb E\sum_i\varphi(\lambda_i(H))=\int\varphi\rho`$ für beschränkte Borel-Funktionen
$`\varphi`$, $`\int\rho=n`$, und $`\rho`$ ist stetig mit $`0\le\rho(x)\le C_n(1+|x|)^{\beta(n-1)}e^{-x^2/2}`$. Insbesondere ist die
Voraussetzung von Satz 5.15 für unabhängige GOE- oder GUE-Blöcke erfüllt.

*Beweis.* Nach Korollar C.6 mit $`g(x)=\sum_i\varphi(x_i)`$ (für $`\varphi\ge0`$; im Allgemeinen wende man es auf $`\varphi^\pm=\max(\pm\varphi,0)`$ an) und der Symmetrie von $`q`$ ist
$`\mathbb E\sum_i\varphi(\lambda_i)=n\int\varphi(x_1)q(x)\,dx=\int\varphi\rho`$ (Fubini). Nach der Schranke für $`|\Delta|`$ ist
$`q(x,y)\le(n!Z_{n,\beta})^{-1}(1+|x|)^{\beta(n-1)}e^{-x^2/2}\prod_{i}(1+|y_i|)^{\beta(n-1)}e^{-y_i^2/2}`$, was nach Integration über $`y`$ die Schranke
liefert und, da $`\sup_x(1+|x|)^{\beta(n-1)}e^{-x^2/2}<\infty`$, eine von $`x`$ unabhängige integrierbare Majorante: Die Stetigkeit
folgt aus dem Satz von der majorisierten Konvergenz. $`\square`$

*Rechnung:* Satz C.7 und Korollar C.8 für das GOE mit $`n=3`$ — [Tutorial, Abschnitt C.3](../../tutorial/de/C_eigenvalue_densities.md#s-C-3).

<a id="r-C-9"></a>

**Proposition C.9 (kleinster Abstand).** Für das GOE oder GUE der Größe $`n\ge2`$ sei $`\delta=\min_i(\lambda_{i+1}-\lambda_i)`$. Es gibt
Konstanten $`0<c\le C`$ (abhängig von $`n,\beta`$) mit $`c\,\varepsilon^{\beta+1}\le P(\delta<\varepsilon)\le C\varepsilon^{\beta+1}`$ für $`0<\varepsilon\le1`$.

*Beweis.* *Obere Schranke.* $`P(\delta<\varepsilon)\le\sum_iP(\lambda_{i+1}-\lambda_i<\varepsilon)`$. Auf $`\{0<u=x_{i+1}-x_i<\varepsilon\}`$ schätze man den Faktor des
Paares $`(i,i+1)`$ durch $`u^\beta`$ ab, die übrigen wie oben, sowie $`1+|x_{i+1}|\le2(1+|x_i|)`$ und $`e^{-x_{i+1}^2/2}\le1`$; mit der Substitution
$`x_{i+1}=x_i+u`$ folgt
```math
P(\lambda_{i+1}-\lambda_i<\varepsilon)\le\frac{2^{\beta(n-1)}}{Z_{n,\beta}}\int_0^\varepsilon u^\beta du\int(1+|x|)^{2\beta(n-1)}e^{-\frac{x^2}2}dx
\Big(\int(1+|x|)^{\beta(n-1)}e^{-\frac{x^2}2}dx\Big)^{n-2}.
```
*Untere Schranke.* Sei $`R_\varepsilon=\{x:0<x_1<1,\ 0<x_2-x_1<\varepsilon,\ 3k<x_k<3k+1\ (3\le k\le n)\}\subseteq C`$. Auf $`R_\varepsilon`$ ist
$`\delta\le x_2-x_1`$, jedes Paar außer $`(1,2)`$ hat Abstand $`\ge1`$, und $`|x|^2\le M`$ mit einem von $`\varepsilon\le1`$ unabhängigen $`M`$; also ist
$`P(\delta<\varepsilon)\ge P(\lambda\in R_\varepsilon)\ge Z_{n,\beta}^{-1}e^{-M/2}\int_0^\varepsilon u^\beta du`$. $`\square`$

Für das GOE ist das die am Ende von Abschnitt 5.4 zitierte Ordnung $`\varepsilon^2`$. Für Haar-verteiltes $`U\in O(3)`$ ist $`\mathbb E|U_{11}|^2=\frac13`$: Linksmultiplikation mit einer Permutationsmatrix permutiert die
Einträge der ersten Spalte und erhält das Haar-Maß, also sind die drei $`\mathbb E|U_{i1}|^2`$ gleich, und ihre Summe ist $`1`$.

<a id="s-C-4"></a>

## C.4 Orthogonalpolynome und die Einpunktfunktion

<a id="r-C-10"></a>

**Lemma C.10 (Einpunktfunktion).** Sei $`I\subseteq\mathbb R`$ ein Intervall und $`w\ge0`$ eine Borel-Funktion auf $`I`$, die auf einer
Menge positiven Maßes positiv ist, mit $`\int_I|x|^mw<\infty`$ für alle $`m`$. Seien $`p_0,\dots,p_{n-1}`$ reelle Polynome, $`\deg p_k=k`$,
mit Leitkoeffizient $`\kappa_k\neq0`$, orthonormal in $`L^2(w)`$ (sie existieren: Man wende Gram-Schmidt auf $`1,x,\dots,x^{n-1}`$ an, die
in $`L^2(w)`$ linear unabhängig sind, da ein vom Nullpolynom verschiedenes Polynom nur endlich viele Nullstellen hat; das Vorzeichen jedes $`p_k`$ kann beliebig gewählt werden). Sei $`q=Z_n^{-1}\Delta(x)^2\prod_iw(x_i)`$ auf
$`I^n`$ mit $`Z_n=\int_{I^n}\Delta^2\prod_iw(x_i)\,dx`$. Dann ist $`Z_n=n!\prod_k\kappa_k^{-2}`$ und
```math
n\int_{I^{n-1}}q(x,y)\,dy=K_n(x)\,w(x),\qquad K_n(x)=\sum_{k=0}^{n-1}p_k(x)^2=v(x)^TM^{-1}v(x),
```
wobei $`v(x)=(1,x,\dots,x^{n-1})`$ und $`M_{ij}=\int_Ix^{i+j}w`$, $`0\le i,j<n`$.

*Beweis.* $`\Delta(x)=\det(x_j^{i-1})_{i,j}`$. Addiert man zu jeder Zeile Vielfache der darüberliegenden Zeilen, so ersetzt man $`x_j^{i-1}`$ durch
$`p_{i-1}(x_j)/\kappa_{i-1}`$, also ist $`\Delta(x)=(\prod_k\kappa_k)^{-1}\det(p_{i-1}(x_j))`$ und
```math
\Delta(x)^2=\Big(\prod_k\kappa_k\Big)^{-2}\sum_{\sigma,\tau\in S_n}\mathrm{sgn}\,\sigma\mathrm{sgn}\,\tau\prod_{j=1}^np_{\sigma(j)-1}(x_j)\,p_{\tau(j)-1}(x_j).
```
Integriert man über $`x_2,\dots,x_n`$ bezüglich $`\prod_{j\ge2}w(x_j)`$, so liefert der Faktor $`j`$ den Wert $`\delta_{\sigma(j)\tau(j)}`$; also bleibt nur $`\sigma=\tau`$
übrig, und das Ergebnis ist $`(\prod_k\kappa_k)^{-2}\sum_\sigma p_{\sigma(1)-1}(x_1)^2=(\prod_k\kappa_k)^{-2}(n-1)!\,K_n(x_1)`$.
Integration bezüglich $`w(x_1)`$ ergibt $`Z_n=(\prod_k\kappa_k)^{-2}n!`$, und Division ergibt die Formel. Ist schließlich
$`p_k(x)=\sum_iC_{ki}x^i`$, so lautet die Orthonormalität $`CMC^T=1`$ mit invertierbarem $`C`$, also ist $`M^{-1}=C^TC`$ und
$`v^TM^{-1}v=|Cv|^2=\sum_kp_k^2`$. $`\square`$

Nach Korollar C.6 haben die ungeordneten Eigenwerte des GUE die Dichte $`q`$ aus Lemma C.10 mit $`w(x)=e^{-x^2/2}`$ auf
$`\mathbb R`$, also ist die Niveaudichte des GUE gleich $`K_n(x)e^{-x^2/2}`$.

*Rechnung:* Die Niveaudichte des GUE für $`n=3`$ — [Tutorial, Abschnitt C.4](../../tutorial/de/C_eigenvalue_densities.md#s-C-4).

<a id="s-C-5"></a>

## C.5 Wishart-Matrizen und das induzierte Maß

Sei $`a\le b`$, $`G`$ eine $`a\times b`$-Ginibre-Matrix (unabhängige Einträge mit unabhängigen $`N(0,\frac12)`$-verteilten Real- und Imaginärteilen,
Dichte $`\pi^{-ab}e^{-\mathrm{Tr}\,GG^*}`$) und $`W=GG^*`$. Sei $`T_+`$ die Menge der unteren Dreiecksmatrizen der Größe $`a\times a`$ mit positiver
Diagonale, eine offene Teilmenge von $`\mathbb R^{a^2}`$ in den Koordinaten $`L_{ii}`$, $`\mathrm{Re}\,L_{ij}`$, $`\mathrm{Im}\,L_{ij}`$ ($`i>j`$), und
$`P_a`$ die Menge der positiv definiten Matrizen, offen in den hermiteschen $`a\times a`$-Matrizen.

<a id="r-C-11"></a>

**Lemma C.11 (Bartlett-Zerlegung).** Fast sicher sind die Zeilen $`r_1,\dots,r_a`$ von $`G`$ linear unabhängig, und
Gram-Schmidt liefert $`G=LV`$ mit $`L\in T_+`$ und $`V`$ mit orthonormalen Zeilen, also $`W=LL^*`$. Die Einträge $`L_{ij}`$, $`i\ge j`$, sind
unabhängig, $`L_{ij}`$ ($`i>j`$) ist standard-komplex-Gaußsch, und $`L_{ii}^2`$ hat die Dichte $`x^{b-i}e^{-x}/(b-i)!`$ auf $`(0,\infty)`$.
Folglich hat $`L`$ die Dichte $`c_L\prod_iL_{ii}^{2(b-i)+1}e^{-\mathrm{Tr}\,LL^*}`$ auf $`T_+`$.

*Beweis.* Sei $`\mathcal F_{i-1}=\sigma(r_1,\dots,r_{i-1})`$, und es seien (Induktion) $`r_1,\dots,r_{i-1}`$ f.s. linear unabhängig mit den
Gram-Schmidt-Vektoren $`v_1,\dots,v_{i-1}`$. Sei $`u_1,\dots,u_b`$ die Orthonormalbasis, die Gram-Schmidt aus
$`v_1,\dots,v_{i-1},e_1,\dots,e_b`$ liefert (wobei abhängige Vektoren weggelassen werden); $`u_k=v_k`$ für $`k<i`$, und die $`u_k`$ sind Borel-Funktionen von
$`r_1,\dots,r_{i-1}`$. Setze $`c_k=\langle u_k,r_i\rangle`$. Da $`r_i`$ von $`\mathcal F_{i-1}`$ unabhängig ist und seine Dichte $`\pi^{-b}e^{-\|r\|^2}`$
unter unitären Abbildungen invariant ist (Proposition A.16(a)), ist der Vektor $`c`$ bedingt auf $`\mathcal F_{i-1}`$ standard-komplex-Gaußsch (Fubini). Nun ist
$`L_{ij}=c_j`$ für $`j<i`$ und $`r_i-\sum_{j<i}c_jv_j=\sum_{k\ge i}c_ku_k`$, also $`L_{ii}^2=\sum_{k=i}^b|c_k|^2>0`$ f.s. Die bedingte Verteilung
der Zeile $`i`$ von $`L`$ gegeben $`\mathcal F_{i-1}`$ ist somit eine feste Verteilung, also ist Zeile $`i`$ von den früheren Zeilen unabhängig, mit unabhängigen Koordinaten
$`c_1,\dots,c_{i-1}`$ und $`L_{ii}`$; per Induktion sind alle Zeilen von $`L`$ unabhängig. Nach Satz A.20(b) mit $`N=m=b-i+1`$
hat $`L_{ii}^2=\|(c_i,\dots,c_b)\|^2`$ die Dichte $`s^{m-1}e^{-s}/(m-1)!`$. Die Dichte von $`L_{ii}=\sqrt s`$ ist dann
$`2t^{2(b-i)+1}e^{-t^2}/(b-i)!`$, und $`\mathrm{Tr}\,LL^*=\sum_{i\ge j}|L_{ij}|^2`$. $`\square`$

*Rechnung:* Lemma C.11 über den Cholesky-Faktor von $`W`$ — [Tutorial, Abschnitt C.5](../../tutorial/de/C_eigenvalue_densities.md#s-C-5).

<a id="r-C-12"></a>

**Lemma C.12 (Cholesky-Abbildung).** $`\chi(L)=LL^*`$ ist eine $`C^\infty`$-Bijektion von $`T_+`$ auf $`P_a`$ mit
$`|\det D\chi(L)|=2^a\prod_{i=1}^aL_{ii}^{2(a-i)+1}`$.

*Beweis.* *Surjektiv:* Für $`W\in P_a`$ wende man Lemma 5.2 auf $`W^{1/2}=QR`$ an; dann ist $`W=R^*R`$ und $`L=R^*\in T_+`$. *Injektiv:* Ist
$`LL^*=MM^*`$, so ist $`M^{-1}L=M^*(L^*)^{-1}`$ eine untere und obere Dreiecksmatrix, also diagonal mit positiver Diagonale, und
unitär, also gleich $`1`$. *Jacobi-Determinante:* Für $`i\ge j`$ ist $`W_{ij}=L_{ij}L_{jj}+\sum_{k<j}L_{ik}\overline{L_{jk}}`$. Man ordne die Koordinaten von $`L`$ und des
unteren Dreiecks von $`W`$ spaltenweise und innerhalb der Spalte $`j`$ die Diagonale zuerst. Dann hängt $`W_{ij}`$ nur von
$`L_{ij}`$ und von früheren Koordinaten ab, also ist $`D\chi`$ eine Dreiecksmatrix mit den Diagonaleinträgen $`2L_{jj}`$ ($`\partial W_{jj}/\partial L_{jj}`$) und,
für jedes $`i>j`$, zweimal $`L_{jj}`$ (Real- und Imaginärteil von $`W_{ij}`$). Das Produkt ist $`\prod_j2L_{jj}\cdot L_{jj}^{2(a-j)}`$. $`\square`$

<a id="r-C-13"></a>

**Satz C.13 (Wishart- und Laguerre-Ensemble).** $`W=GG^*`$ hat die Dichte $`c_W\det(W)^{b-a}e^{-\mathrm{Tr}\,W}`$ auf $`P_a`$
(null außerhalb), und seine Eigenwerte haben die Dichte
```math
Z_{a,b}^{-1}\,\Delta(\lambda)^2\prod_{i=1}^a\lambda_i^{\,b-a}e^{-\lambda_i}\quad\text{on }C_+=\{0<\lambda_1<\dots<\lambda_a\}.
```

*Beweis.* Für Borel-Funktionen $`g\ge0`$ liefern Lemma C.11 und die Variablentransformation $`W=\chi(L)`$ (Lemma C.12)
$`\mathbb E\,g(W)=\int_{P_a}g(W)\,c_L\prod_iL_{ii}^{2(b-i)+1}e^{-\mathrm{Tr}\,W}\big/\big(2^a\prod_iL_{ii}^{2(a-i)+1}\big)dW`$ mit
$`L=\chi^{-1}(W)`$, und $`\prod_iL_{ii}^{2(b-a)}=\det(W)^{b-a}`$. Diese Dichte ist invariant unter $`W\mapsto VWV^*`$, $`V\in U(a)`$, und
$`\Lambda\in P_a`$ genau dann, wenn $`\lambda_1>0`$; man wende Korollar C.6 mit $`\beta=2`$ an. $`\square`$

<a id="r-C-14"></a>

**Satz C.14 (induziertes Maß).** Sei $`\psi`$ ein Haar-zufälliger reiner Zustand von $`\mathbb C^a\otimes\mathbb C^b`$ (Abschnitt 7.5),
$`\rho_A=\mathrm{Tr}_2|\psi\rangle\langle\psi|=W/s`$ mit $`s=\mathrm{Tr}\,W`$ (Beweis von Satz 7.11), und seien
$`\mu_1<\dots<\mu_a`$ die Eigenwerte von $`\rho_A`$. Dann sind $`s`$ und $`\mu`$ unabhängig, $`s`$ hat die Dichte
$`s^{ab-1}e^{-s}/(ab-1)!`$, und $`(\mu_1,\dots,\mu_{a-1})`$ hat die Dichte
```math
Z'^{-1}_{a,b}\,\Delta(\mu)^2\prod_{i=1}^a\mu_i^{\,b-a}\quad\text{on }D_a=\{0<\mu_1<\dots<\mu_{a-1}<\mu_a\},\qquad\mu_a=1-\sum_{i<a}\mu_i ,
```
die Dichte „$`\propto\delta(1-\sum_i\mu_i)\,\Delta(\mu)^2\prod_i\mu_i^{b-a}`$“ der physikalischen Literatur.

*Beweis.* $`\Theta(s,\mu_1,\dots,\mu_{a-1})=s\,(\mu_1,\dots,\mu_{a-1},1-\sum_{i<a}\mu_i)`$ ist eine $`C^\infty`$-Bijektion von $`(0,\infty)\times D_a`$ auf
$`C_+`$ mit der Umkehrabbildung $`\lambda\mapsto(\sum\lambda_i,\lambda_1/\sum\lambda_i,\dots)`$. Addiert man alle Zeilen von $`D\Theta`$ zur letzten, so wird diese zu
$`(1,0,\dots,0)`$, also ist $`|\det D\Theta|=s^{a-1}`$. Da $`\Delta`$ homogen vom Grad $`a(a-1)/2`$ ist, ist die Dichte von $`(s,\mu)`$
```math
Z_{a,b}^{-1}\,s^{a(a-1)+a(b-a)+a-1}e^{-s}\cdot\Delta(\mu)^2\prod_i\mu_i^{b-a}=Z_{a,b}^{-1}\,s^{ab-1}e^{-s}\cdot\Delta(\mu)^2\prod_i\mu_i^{b-a},
```
ein Produkt, was die Unabhängigkeit beweist und nach Normierung beide Faktoren identifiziert. $`\square`$

*Rechnung:* Satz C.14 für $`a=2`$ — [Tutorial, Abschnitt C.5](../../tutorial/de/C_eigenvalue_densities.md#s-C-5). *In der Bibliothek:* Das induzierte Maß, `states.random_density_matrix` mit `rank=b` — [Tutorial, Abschnitt C.5](../../tutorial/de/C_eigenvalue_densities.md#s-C-5).

<a id="s-C-6"></a>

## C.6 Laguerre-Polynome und die Pagesche Formel

<a id="r-C-15"></a>

**Lemma C.15 (Laguerre-Polynome).** Für ganze Zahlen $`\alpha,k\ge0`$ sei
$`L_k^{(\alpha)}(x)=\sum_{j=0}^k(-1)^j\binom{k+\alpha}{k-j}\frac{x^j}{j!}`$ und $`L_{-1}^{(\alpha)}=0`$. Dann gilt:

(a) $`L_k^{(\alpha)}(x)\,x^\alpha e^{-x}=\frac1{k!}\frac{d^k}{dx^k}\big(x^{k+\alpha}e^{-x}\big)`$; $`L_k^{(\alpha)}`$ hat den Grad $`k`$ und den Leitkoeffizienten $`(-1)^k/k!`$.

(b) $`\frac{d}{dx}L_k^{(\alpha)}=-L_{k-1}^{(\alpha+1)}`$, also $`\frac{d^i}{dx^i}L_k^{(\alpha)}=(-1)^iL_{k-i}^{(\alpha+i)}`$ für $`i\le k`$; und $`L_k^{(\alpha)}=L_k^{(\alpha+1)}-L_{k-1}^{(\alpha+1)}`$.

(c) $`\int_0^\infty x^mL_k^{(\alpha)}x^\alpha e^{-x}dx=0`$ für $`m<k`$, und $`\int_0^\infty L_k^{(\alpha)}L_l^{(\alpha)}x^\alpha e^{-x}dx=\delta_{kl}\frac{(k+\alpha)!}{k!}`$.

Folglich haben die Eigenwerte von $`W`$ nach Lemma C.10 und Satz C.13 die Niveaudichte
($`\mathbb E\sum_i\varphi(\lambda_i(W))=\int_0^\infty\varphi\,\rho_{a,b}`$ für Borel-Funktionen $`\varphi\ge0`$ und für jedes $`\varphi`$ mit $`\int|\varphi|\rho_{a,b}<\infty`$)
```math
\rho_{a,b}(x)=\sum_{k=0}^{a-1}\frac{k!}{(k+\alpha)!}L_k^{(\alpha)}(x)^2\,x^\alpha e^{-x},\qquad\alpha=b-a .
```

*Beweis.* (a) Nach Leibniz ist $`\frac{d^k}{dx^k}(x^{k+\alpha}e^{-x})=\sum_j\binom kj(-1)^j\frac{(k+\alpha)!}{(\alpha+j)!}x^{\alpha+j}e^{-x}`$ (mit $`j`$ Ableitungen
auf $`e^{-x}`$), und $`\frac1{k!}\binom kj\frac{(k+\alpha)!}{(\alpha+j)!}=\binom{k+\alpha}{k-j}\frac1{j!}`$. (b) Gliedweises Differenzieren ergibt
$`-\sum_{i=0}^{k-1}(-1)^i\binom{k+\alpha}{k-1-i}\frac{x^i}{i!}=-L_{k-1}^{(\alpha+1)}`$; die zweite Identität ist die Pascalsche Regel
$`\binom{k+\alpha+1}{k-j}-\binom{k+\alpha}{k-1-j}=\binom{k+\alpha}{k-j}`$. (c) Man integriere (a) $`k`$-mal partiell gegen $`x^m`$: Die
Randterme verschwinden, weil $`\frac{d^r}{dx^r}(x^{k+\alpha}e^{-x})`$ den Faktor $`x^{k+\alpha-r}`$, $`r<k`$, hat und bei $`\infty`$ abfällt. Das
Ergebnis ist $`\frac{(-1)^k}{k!}\int\frac{d^kx^m}{dx^k}x^{k+\alpha}e^{-x}dx`$, also $`0`$ für $`m<k`$ und $`(-1)^k(k+\alpha)!`$ für $`m=k`$. Für $`l<k`$
ist $`L_l^{(\alpha)}`$ eine Linearkombination von $`x^m`$, $`m<k`$; für $`l=k`$ zählt nur der Leitkoeffizient $`(-1)^k/k!`$. Also sind die
$`p_k=\big(k!/(k+\alpha)!\big)^{1/2}L_k^{(\alpha)}`$ orthonormal für $`w(x)=x^\alpha e^{-x}`$ auf $`(0,\infty)`$, und die ungeordneten
Eigenwerte von $`W`$ haben nach Korollar C.6 die Dichte $`q`$ aus Lemma C.10 (mit $`n=a`$), sodass wie im Beweis von Korollar C.8
$`\mathbb E\sum_i\varphi(\lambda_i)=a\int\varphi(x_1)q(x)\,dx=\int\varphi\,K_aw=\int\varphi\,\rho_{a,b}`$ für Borel-Funktionen $`\varphi\ge0`$ gilt, und wegen Linearität für
$`\varphi`$ mit $`\int|\varphi|\rho_{a,b}<\infty`$. $`\square`$

*Rechnung:* Die Niveaudichte $`\rho_{a,b}`$ im Vergleich mit einem Histogramm — [Tutorial, Abschnitt C.6](../../tutorial/de/C_eigenvalue_densities.md#s-C-6).

<a id="r-C-16"></a>

**Lemma C.16 (logarithmische Integrale).** Sei $`I_m=\int_0^\infty x^m\log x\,e^{-x}dx`$ und $`H_m=\sum_{i=1}^m\frac1i`$. Dann ist
$`I_m=m!\,(I_0+H_m)`$. Für ganze Zahlen $`\alpha\ge0`$ und $`0\le k\le m`$ gilt
```math
A^{(\alpha)}_{k,m}:=\int_0^\infty\log x\,L_k^{(\alpha)}(x)L_m^{(\alpha)}(x)\,x^\alpha e^{-x}dx=
\begin{cases}I_{k+\alpha}/k!, & m=k,\\ -\dfrac{(k+\alpha)!}{k!\,(m-k)}, & m>k.\end{cases}
```
($`I_0`$ hebt sich in Satz C.17 weg, sein Wert wird also nicht benötigt.)

*Beweis.* Für $`m\ge1`$ liefert partielle Integration mit $`(-e^{-x})'=e^{-x}`$ die Beziehung $`I_m=mI_{m-1}+(m-1)!`$, d. h.
$`I_m/m!=I_{m-1}/(m-1)!+1/m`$. Für $`A_{k,m}`$ schreibe man $`L_m^{(\alpha)}x^\alpha e^{-x}=\frac1{m!}g^{(m)}`$ mit $`g=x^{m+\alpha}e^{-x}`$
(Lemma C.15(a)) und $`f=\log x\,L_k^{(\alpha)}`$ und integriere $`m`$-mal partiell. Die Randterme
$`f^{(i)}g^{(m-1-i)}`$ verschwinden bei $`\infty`$ und bei $`0`$, wo $`f^{(i)}=O(|\log x|+x^{-i})`$ und $`g^{(m-1-i)}=O(x^{\alpha+1+i})`$. Also ist
$`A_{k,m}=\frac{(-1)^m}{m!}\int_0^\infty f^{(m)}g\,dx`$ (der Integrand ist $`O(x^\alpha(1+|\log x|))`$ bei $`0`$) mit
$`f^{(m)}=\sum_{i\le k}\binom mi(\log x)^{(m-i)}\frac{d^i}{dx^i}L_k^{(\alpha)}`$. Für $`m-i\ge1`$ ist $`(\log x)^{(m-i)}=(-1)^{m-i-1}(m-i-1)!\,x^{i-m}`$,
und nach Lemma C.15(b,c) gilt
```math
\int_0^\infty x^{i-m}\frac{d^iL_k^{(\alpha)}}{dx^i}\,x^{m+\alpha}e^{-x}dx=(-1)^i\int_0^\infty L_{k-i}^{(\alpha+i)}\,x^{\alpha+i}e^{-x}dx=(-1)^i(\alpha+i)!\,\delta_{ik}.
```
Ist $`m>k`$, so trägt nur $`i=k`$ bei: $`A_{k,m}=\frac{(-1)^m}{m!}\binom mk(-1)^{m-k-1}(m-k-1)!(-1)^k(\alpha+k)!=-\frac{(k+\alpha)!}{k!(m-k)}`$.
Ist $`m=k`$, so verschwinden die Terme $`i<k`$, und der Term $`i=k`$ ist $`\log x\cdot(-1)^k`$, was $`A_{k,k}=\frac{(-1)^k}{k!}(-1)^kI_{k+\alpha}`$ ergibt. $`\square`$

*Rechnung:* Die Integrale $`A^{(\alpha)}_{k,m}`$ durch Quadratur — [Tutorial, Abschnitt C.6](../../tutorial/de/C_eigenvalue_densities.md#s-C-6).

<a id="r-C-17"></a>

**Satz C.17 (Pagesche Formel).** Für einen Haar-zufälligen reinen Zustand von $`\mathbb C^a\otimes\mathbb C^b`$ mit $`a\le b`$ gilt
```math
\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}.
```

*Beweis.* *Reduktion.* In der Notation von Satz C.14 ist $`\lambda_i=s\mu_i`$ und $`\sum_i\mu_i=1`$, also
$`\sum_i\lambda_i\log\lambda_i=s\log s-s\,S(\rho_A)`$. Alle Terme sind integrierbar ($`|x\log x|\le1+x^2`$, und $`\mathbb Es^2<\infty`$,
$`0\le S(\rho_A)\le a/e`$, da $`0\le-x\log x\le1/e`$ auf $`[0,1]`$), nach Lemma C.15 ist $`\mathbb E\sum_i\lambda_i\log\lambda_i=\int x\log x\,\rho_{a,b}(x)\,dx`$, und $`s`$ ist unabhängig von
$`S(\rho_A)`$, also $`\mathbb E\sum_i\lambda_i\log\lambda_i=\mathbb E[s\log s]-\mathbb E[s]\,\mathbb E\,S(\rho_A)`$. Aus der Dichte von $`s`$ folgt
$`\mathbb Es=ab`$ und $`\mathbb E[s\log s]=I_{ab}/(ab-1)!=ab\,(I_0+H_{ab})`$ (Lemma C.16). Daher ist
```math
\mathbb E\,S(\rho_A)=I_0+H_{ab}-\frac1{ab}\int_0^\infty x\log x\,\rho_{a,b}(x)\,dx .
```
*Das Integral.* Nach Lemma C.15(b) ist $`L_k^{(\alpha)}=L_k^{(\alpha+1)}-L_{k-1}^{(\alpha+1)}`$, also mit $`x\cdot x^\alpha=x^{\alpha+1}`$ und Lemma C.16
für $`\alpha+1`$
```math
\frac{k!}{(k+\alpha)!}\int_0^\infty x\log x\,L_k^{(\alpha)}(x)^2x^\alpha e^{-x}dx=\frac{k!}{(k+\alpha)!}\Big(A^{(\alpha+1)}_{k,k}-2A^{(\alpha+1)}_{k-1,k}+A^{(\alpha+1)}_{k-1,k-1}\Big)
```
```math
=(k+\alpha+1)(I_0+H_{k+\alpha+1})+2k+k\,(I_0+H_{k+\alpha}),
```
auch für $`k=0`$. Da $`\sum_{k=0}^{a-1}(2k+\alpha+1)=a(a-1)+a(b-a+1)=ab`$, folgt
```math
\int_0^\infty x\log x\,\rho_{a,b}\,dx=ab\,I_0+\Sigma,\qquad \Sigma=\sum_{k=0}^{a-1}\big[(k+\alpha+1)H_{k+\alpha+1}+kH_{k+\alpha}+2k\big].
```
*Die Summe.* Setze $`m=\alpha+1=b-a+1`$ und $`j=k+m\in\{m,\dots,b\}`$. Wegen $`(j-m)H_{j-1}=(j-m)H_j-1+\frac mj`$ ist
$`\Sigma=\sum_{j=m}^b(2j-m)H_j-a+m(H_b-H_{m-1})+a(a-1)`$. Mit $`F(j)=j(j+1-m)`$ gilt $`F(j)-F(j-1)=2j-m`$, $`F(m-1)=0`$,
$`F(b)=ab`$, und Abelsche partielle Summation liefert
```math
\sum_{j=m}^b(2j-m)H_j=abH_b-\sum_{j=m}^{b-1}\frac{F(j)}{j+1}=abH_b-\sum_{j=m}^{b-1}\Big(j-m+\frac m{j+1}\Big)=abH_b-\frac{(a-1)(a-2)}2-m(H_b-H_m).
```
Wegen $`m(H_m-H_{m-1})=1`$ ist $`\Sigma=abH_b-\frac{(a-1)(a-2)}2+1-a+a(a-1)=abH_b+\frac{a(a-1)}2`$. Folglich ist
$`\mathbb E\,S(\rho_A)=H_{ab}-H_b-\frac{a-1}{2b}`$. $`\square`$

*Rechnung:* Die Identität für $`\Sigma`$ in exakter Arithmetik und die Pagesche Formel durch Quadratur ($`a=2`$) und per Monte-Carlo — [Tutorial, Abschnitt C.6](../../tutorial/de/C_eigenvalue_densities.md#s-C-6). *In der Bibliothek:* Die Entropie $`S(\rho_A)`$ ist `DensityMatrix.von_neumann_entropy` — [Tutorial, Abschnitt 3.4](../../tutorial/de/03_traces_and_states.md#s-3-4).

<a id="r-C-18"></a>

**Korollar C.18.** Für $`a\le b`$ gilt
$`\log a-\frac{a-1}{ab+1}-\frac{a-1}{2b}\le\mathbb E\,S(\rho_A)\le\log a-\frac{a-1}{2b}`$; insbesondere ist
$`\big|\mathbb E\,S(\rho_A)-\log a+\frac a{2b}\big|\le\frac1{2b}`$. Für zwei Qubits ist $`\mathbb E\,S(\rho_A)=\frac13`$.

*Beweis.* Da $`1/x`$ fällt, ist $`\int_k^{k+1}\frac{dx}x\le\frac1k\le\int_{k-1}^k\frac{dx}x`$; Summation über $`b<k\le ab`$ ergibt
$`\log\frac{ab+1}{b+1}\le H_{ab}-H_b\le\log a`$, und $`\log\frac{ab+1}{b+1}=\log a-\log\big(1+\frac{a-1}{ab+1}\big)\ge\log a-\frac{a-1}{ab+1}`$.
Schließlich ist $`\frac{a-1}{ab+1}\le\frac1b`$, und $`a=b=2`$ ergibt $`\frac13+\frac14-\frac14`$. $`\square`$

*Rechnung:* Die beiden Schranken für $`a<12`$, $`b<40`$ — [Tutorial, Abschnitt C.6](../../tutorial/de/C_eigenvalue_densities.md#s-C-6).

Die Formel wurde von D. N. Page vermutet (Phys. Rev. Lett. 71, 1291, 1993) und von S. K. Foong und S. Kanno
(Phys. Rev. Lett. 72, 1148, 1994), J. Sánchez-Ruiz (Phys. Rev. E 52, 5653, 1995) und S. Sen (Phys. Rev. Lett. 77, 1, 1996) bewiesen;
die obige Rechnung folgt Sens Reduktion auf das Laguerre-Ensemble, summiert aber die Einpunktdichte Term für Term,
statt die Christoffel-Darboux-Formel zu verwenden. Das induzierte Maß ist das von K. Życzkowski und H.-J. Sommers
(J. Phys. A 34, 7111, 2001); die Zerlegung aus Lemma C.11 geht auf M. S. Bartlett (1933) zurück, Satz C.5 auf H. Weyl.

**Zurück:** [Kapitel 5](05_random_matrices.md) verwendet Satz C.7, Korollar C.8 und Proposition C.9;
[Kapitel 7](07_tensor_products.md) verwendet die Sätze C.14 und C.17 und Korollar C.18.
