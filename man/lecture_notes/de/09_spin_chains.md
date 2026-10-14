[English](../09_spin_chains.md) · [Español](../es/09_spin_chains.md) · [Français](../fr/09_spin_chains.md) · **Deutsch** · [中文](../zh/09_spin_chains.md) · [日本語](../ja/09_spin_chains.md) · [Русский](../ru/09_spin_chains.md)

# 9. Spinketten

Eine Kette aus $`L`$ Spins $`\frac12`$ lebt auf $`(\mathbb C^2)^{\otimes L}`$, und dichte lineare Algebra in Dimension
$`2^L`$ wird jenseits von $`L\approx14`$ teuer (eine dichte Matrix hat $`4^L`$ Einträge). Symmetrien verschieben diese Grenze: Die Magnetisierung zerlegt den Raum in $`L+1`$
Sektoren (Kapitel 4), Translationen zerlegen jeden Sektor erneut, und $`SU(2)`$ verkleinert die Blöcke noch einmal.
Dieses Kapitel konstruiert diese Sektoren für die XXZ-Kette, beweist, was sie über Grundzustände aussagen, prüft
die Bibliothek an freien Fermionen und am Bethe-Ansatz und geht über dichte Matrizen (Krylov-Verfahren) und
saubere Ketten (Unordnung) hinaus.

**Sie lernen**, warum die Algebra einer Kette $`\bigoplus_NM_{\binom LN}`$ ist und warum der Turm dieser Algebren
das Pascalsche Dreieck als Bratteli-Diagramm und den Index $`4\cos^2\frac\pi{2L+4}`$ hat, wie man Impulssektoren
abzählt, die Marshallsche Vorzeichenregel und ihre Folgen für Impuls und Spin von Grundzuständen (Lieb–Mattis),
die Jordan-Wigner-Lösung der XX-Kette, die Bethe-Gleichungen für zwei Magnonen, warum Lanczos monoton konvergiert
und wie man Niveaustatistik und die Verschränkung freier Fermionen misst. Die numerischen Kontrollen dieses
Kapitels stehen im [Tutorial 9](../../tutorial/de/09_spin_chains.md).

Im ganzen Kapitel sind die Plätze $`0,\dots,L-1`$, die Basis von $`(\mathbb C^2)^{\otimes L}`$ besteht aus den
Konfigurationen $`|b\rangle=|b_0\dots b_{L-1}\rangle`$ mit $`b_i=1`$ für Spin nach oben (Konvention der Bibliothek: Bit
$`i`$ eines Labels ist Platz $`i`$), $`s^a_i`$ ($`a=x,y,z`$) sind die Spinmatrizen $`\frac12\sigma_a`$ auf Platz $`i`$,
$`s^\pm_i=s^x_i\pm is^y_i`$, $`n_i=s^z_i+\frac12`$, $`S^a=\sum_is^a_i`$, und $`N=S^z+\frac L2`$ ist die Anzahl der Spins nach
oben. Eine Kette ist **offen** (Bindungen $`(i,i+1)`$, $`i<L-1`$) oder **periodisch** (zusätzlich die Bindung
$`(L-1,0)`$; dann ist $`L\ge3`$, und Indizes werden modulo $`L`$ genommen). $`J>0`$ (antiferromagnetisch), $`\hbar=1`$, und
$`\log`$ ist der natürliche Logarithmus.

<a id="s-9-1"></a>

## 9.1 Die XXZ-Kette und ihre $`U(1)`$-Sektoren

Die **XXZ-Kette** mit Anisotropie $`\Delta\in\mathbb R`$ und Feldern $`h_i\in\mathbb R`$ ist
```math
H=J\sum_{(i,j)}\Big[\tfrac12\big(s^+_is^-_j+s^-_is^+_j\big)+\Delta\,s^z_is^z_j\Big]+\sum_ih_is^z_i ,
```
summiert über die Bindungen; $`\frac12(s^+_is^-_j+s^-_is^+_j)=s^x_is^x_j+s^y_is^y_j`$ bewegt einen Spin nach oben entlang
der Bindung. Für $`\Delta=1`$ ist dies die Heisenberg-Kette $`J\sum\mathbf s_i\cdot\mathbf s_j`$, für $`\Delta=0`$ die
XX-Kette.

<a id="r-9-1"></a>

**Proposition 9.1.** $`[H,S^z]=0`$. Folglich ist $`H\in U(1)'\cong\bigoplus_{N=0}^LM_{\binom LN}(\mathbb C)`$, wobei der
Block $`N`$ auf dem Spann der Konfigurationen mit $`N`$ Spins nach oben wirkt.

*Beweis.* Die Diagonalterme kommutieren mit dem diagonalen $`S^z`$. Für $`i\neq j`$ kommutiert $`s^+_is^-_j`$ mit $`s^z_l`$,
$`l\ne i,j`$, und $`[s^z_i+s^z_j,s^+_is^-_j]=[s^z_i,s^+_i]s^-_j+s^+_i[s^z_j,s^-_j]=s^+_is^-_j-s^+_is^-_j=0`$ wegen
$`[s^z,s^\pm]=\pm s^\pm`$ (Abschnitt 4.5); für $`s^-_is^+_j`$ geht man zu den Adjungierten über. Folglich
kommutiert $`H`$ mit $`N=S^z+\frac L2`$, das die ganzzahligen Eigenwerte $`0,\dots,L`$ hat, mit Eigenräumen, die von den
Konfigurationen mit $`N`$ Spins nach oben aufgespannt werden und die Dimension $`\binom LN`$ haben; der Rest ist
Proposition 4.10(2) mit $`Q=N`$. $`\square`$

*Rechnung:* Proposition 9.1 und die Blöcke $`H_N`$, die `SpinChain.xxz` aus den Konfigurationen jedes Sektors (`chain.labels[N]`) baut, im Vergleich mit der vollen $`2^L\times2^L`$-Matrix — [Tutorial, Abschnitt 9.1](../../tutorial/de/09_spin_chains.md#s-9-1).

**Der Turm der Ketten.** Sei $`\mathcal A_L=\bigoplus_NM_{\binom LN}`$, wirkend auf $`L`$ Plätzen. Das Hinzufügen eines
Platzes, $`x\mapsto x\otimes1_2`$, bildet $`\mathcal A_L`$ in $`\mathcal A_{L+1}`$ ab, denn $`x\otimes1`$ kommutiert mit
$`S^z\otimes1+1\otimes s^z`$. Das ergibt einen Turm $`\mathcal A_0\subset\mathcal A_1\subset\mathcal A_2\subset\cdots`$.

<a id="r-9-2"></a>

**Proposition 9.2.** Die Inklusion $`\mathcal A_L\subset\mathcal A_{L+1}`$ hat die $`(L+2)\times(L+1)`$-Inklusionsmatrix
$`\Lambda_L`$ mit $`\lambda_{N,N}=\lambda_{N,N-1}=1`$ und allen übrigen Einträgen $`0`$; die Blockgrößen gehorchen der
Pascalschen Regel $`\binom{L+1}N=\binom LN+\binom L{N-1}`$, das Bratteli-Diagramm des Turms ist das Pascalsche Dreieck,
und der bipartite Graph von $`\Lambda_L`$ ist ein Pfad mit $`2L+3`$ Ecken.

*Beweis.* Die Konfigurationen mit $`N`$ Spins nach oben auf $`L+1`$ Plätzen sind diejenigen mit $`N`$ Spins nach oben auf
den ersten $`L`$ Plätzen und dem neuen Spin nach unten sowie diejenigen mit $`N-1`$ und dem neuen Spin nach oben. Auf
diesen beiden Unterräumen wirkt $`x\otimes1`$ als $`x_N`$ bzw. $`x_{N-1}`$, also enthält Block $`N`$ von $`\mathcal A_{L+1}`$
die Blöcke $`N`$ und $`N-1`$ von $`\mathcal A_L`$ je einmal; die Dimensionen addieren sich gemäß der Pascalschen Regel. Im
bipartiten Graphen ist der Sektor $`N`$ von $`\mathcal A_{L+1}`$ mit $`N-1`$ und $`N`$ von $`\mathcal A_L`$ verbunden, was den
Pfad $`0'-0-1'-1-\cdots-L'-L-(L+1)'`$ ergibt (Striche für $`\mathcal A_{L+1}`$). $`\square`$

Ein Pfad von der Spitze des Bratteli-Diagramms zum Sektor $`N`$ auf Stufe $`L`$ ist eine Konfiguration $`b`$; davon gibt
es $`\binom LN`$, die Größe des Blocks. Für eine zusammenhängende Inklusion mit Inklusionsmatrix $`\Lambda`$
**definieren** wir den **Index** als $`\|\Lambda\|^2`$ (Operatornorm); für $`M_k\otimes1\subset M_k\otimes M_m`$ ist die
Inklusionsmatrix die $`1\times1`$-Matrix $`(m)`$, also ist der Index das $`m^2`$ aus Abschnitt 7.2. `Inclusion.index`
berechnet ihn. (Der Name ist der Indextheorie von V. Jones entlehnt, *Invent. Math.* 72 (1983) 1; in diesem Kurs ist
der Index einer Inklusion endlichdimensionaler Algebren diese Definition und nichts weiter.)

*In der Bibliothek:* Der Index ist `Inclusion.index` — [Tutorial, Abschnitt 9.1](../../tutorial/de/09_spin_chains.md#s-9-1).

<a id="r-9-3"></a>

**Lemma 9.3.** Die Adjazenzmatrix $`A_n`$ des Pfades mit den Ecken $`1,\dots,n`$ hat die Eigenwerte
$`2\cos\frac{\pi k}{n+1}`$, $`k=1,\dots,n`$, mit Eigenvektoren $`v_k(j)=\sin\frac{\pi jk}{n+1}`$; folglich ist
$`\|A_n\|=2\cos\frac\pi{n+1}`$.

*Beweis.* Man setze $`v_k(0)=v_k(n+1)=0`$. Dann ist $`(A_nv_k)(j)=v_k(j-1)+v_k(j+1)=2\cos\frac{\pi k}{n+1}\,v_k(j)`$ wegen
$`\sin(a-b)+\sin(a+b)=2\sin a\cos b`$, und $`v_k\ne0`$, da $`v_k(1)=\sin\frac{\pi k}{n+1}\neq0`$. Die $`n`$ Eigenwerte sind
verschieden ($`\cos`$ ist auf $`(0,\pi)`$ streng monoton fallend), also sind es alle, und die Norm der symmetrischen
Matrix $`A_n`$ ist der größte Betrag, angenommen bei $`k=1`$ und $`k=n`$. $`\square`$

<a id="r-9-4"></a>

**Korollar 9.4.** $`[\mathcal A_{L+1}:\mathcal A_L]=\|\Lambda_L\|^2=4\cos^2\frac\pi{2L+4}`$, wachsend gegen $`4`$.

*Beweis.* Ordnet man die Ecken entlang des Pfades aus Proposition 9.2, so wird die Adjazenzmatrix
$`B=\begin{pmatrix}0&\Lambda_L\\\Lambda_L^T&0\end{pmatrix}`$ zu $`A_{2L+3}`$ (eine Ähnlichkeit mittels Permutation). Nach
Lemma 1.1 ist $`\|B\|^2=\|B^2\|`$ und $`\|\Lambda_L\Lambda_L^T\|=\|\Lambda_L^T\Lambda_L\|=\|\Lambda_L\|^2`$, und
$`B^2=\Lambda_L\Lambda_L^T\oplus\Lambda_L^T\Lambda_L`$, also $`\|B\|=\|\Lambda_L\|`$; Lemma 9.3 mit $`n=2L+3`$ liefert den Wert. $`\square`$

*Rechnung:* Proposition 9.2, Lemma 9.3 und Korollar 9.4 für $`L\le8`$, mit `inclusions.chain_tower` und `inclusions.bratteli_diagram` — [Tutorial, Abschnitt 9.1](../../tutorial/de/09_spin_chains.md#s-9-1).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). $`S^z`$ ist erhalten, sobald die Wechselwirkung invariant unter Drehungen um die
> $`z`$-Achse ist; $`N`$ zählt die Spins nach oben, oder die Teilchen nach der Jordan-Wigner-Abbildung aus
> Abschnitt 9.4. Im Turm wird jeder Ladungssektor der längeren Kette von zwei Sektoren der kürzeren gespeist.

<a id="s-9-2"></a>

## 9.2 Translationen, Impuls und die Marshallsche Vorzeichenregel

Auf einer periodischen Kette sei $`T`$ die Translation $`T|b_0b_1\dots b_{L-1}\rangle=|b_{L-1}b_0\dots b_{L-2}\rangle`$, die
Platz $`i`$ auf Platz $`i+1`$ verschiebt. $`T`$ ist eine Permutationsmatrix, also unitär, $`T^L=1`$ und
$`Ts^a_iT^*=s^a_{i+1}`$. Folglich ist $`TS^zT^*=S^z`$ und $`THT^*=H`$ für ein homogenes Feld $`h_i=h`$: $`T`$ bildet Bindungen
auf Bindungen ab. Seine Eigenwerte haben daher die Form $`e^{ik}`$ mit **Impulsen** $`k=2\pi m/L`$, $`m=0,\dots,L-1`$.

<a id="r-9-5"></a>

**Proposition 9.5 (Impulssektoren).** Sei $`N`$ fest. Für eine Konfiguration $`a`$ mit $`N`$ Spins nach oben sei $`R_a`$
die Größe ihrer Bahn $`\{T^ja\}`$. (a) $`R_a`$ teilt $`L`$; für $`k=2\pi m/L`$ mit $`mR_a\equiv0\pmod L`$ ist der Vektor
$`|a,k\rangle=\sum_{j=0}^{R_a-1}e^{-ikj}T^j|a\rangle`$ von null verschieden mit $`T|a,k\rangle=e^{ik}|a,k\rangle`$, und diese $`R_a`$
Vektoren bilden eine Orthogonalbasis des Spanns der Bahn. (b) Der Eigenwert $`e^{2\pi im/L}`$ von $`T`$ auf dem Sektor
$`N`$ hat die Vielfachheit $`d_{N,m}=\#\{\text{Bahnen }:mR_a\equiv0\pmod L\}`$, und jeder Operator, der mit $`S^z`$ und $`T`$
kommutiert – etwa die translationsinvariante XXZ-Kette –, ist blockdiagonal, liegt also in
$`\bigoplus_N\bigoplus_mM_{d_{N,m}}`$.

*Beweis.* (a) Der Stabilisator von $`a`$ in der zyklischen Gruppe $`\{T^j\}\cong\mathbb Z_L`$ ist eine Untergruppe, und
$`R_a`$ ist ihr Index, ein Teiler von $`L`$; die Vektoren $`T^j|a\rangle`$, $`0\le j<R_a`$, sind orthonormal, und
$`T^{R_a}|a\rangle=|a\rangle`$. Ist $`e^{ikR_a}=1`$, so gilt $`T|a,k\rangle=e^{ik}\sum_{j=1}^{R_a}e^{-ikj}T^j|a\rangle=e^{ik}|a,k\rangle`$,
da der Term $`j=R_a`$ gleich dem Term $`j=0`$ ist; der Koeffizient von $`|a\rangle`$ ist $`1`$. Die zulässigen
$`m\in\{0,\dots,L-1\}`$ sind die Vielfachen von $`L/R_a`$, insgesamt $`R_a`$ Stück, mit verschiedenen Eigenwerten des
unitären $`T`$; also sind die Vektoren orthogonal und spannen den $`R_a`$-dimensionalen Spann der Bahn auf. (b) Der
Sektor $`N`$ ist die orthogonale Summe der Spanne seiner Bahnen, jeder davon invariant unter $`T`$, und nach (a) tritt
der Eigenwert $`e^{2\pi im/L}`$ im Spann der Bahn von $`a`$ einmal auf, wenn $`mR_a\equiv0\pmod L`$, und sonst gar nicht.
Ein $`X`$, das mit $`S^z`$ und $`T`$ kommutiert, kommutiert
mit deren Spektralprojektionen $`P_N`$, $`Q_m`$ (Satz 1.5(d)), also ist $`X=\sum_{N,m}P_NQ_mXP_NQ_m`$. $`\square`$

*Rechnung:* Proposition 9.5(b) für $`L=12`$ und die Impulsblöcke, die `SpinChain.xxz_momentum` baut — [Tutorial, Abschnitt 9.2](../../tutorial/de/09_spin_chains.md#s-9-2).

<a id="r-9-6"></a>

**Lemma 9.6 (Perron-Frobenius, symmetrischer Fall).** Sei $`M`$ eine reelle symmetrische Matrix mit $`M_{ij}\ge0`$ für
$`i\ne j`$, deren Graph ($`i\sim j`$ genau dann, wenn $`M_{ij}>0`$) zusammenhängend ist. Dann ist ihr größter Eigenwert
$`\lambda`$ einfach, mit einem Eigenvektor, dessen Einträge alle $`>0`$ sind.

*Beweis.* $`\lambda=\max_{\|x\|=1}\langle x,Mx\rangle`$, und das Maximum wird genau auf den Einheitsvektoren des
$`\lambda`$-Eigenraums angenommen (Entwicklung in einer Orthonormalbasis aus Eigenvektoren). Für ein solches $`x`$ erfüllt
der Vektor $`|x|`$ der Beträge $`\langle|x|,M|x|\rangle\ge\mathrm{Re}\langle x,Mx\rangle=\lambda`$, denn die Diagonalterme
stimmen überein und $`M_{ij}|x_i||x_j|\ge\mathrm{Re}\,M_{ij}\bar x_ix_j`$ für $`i\ne j`$; also ist $`|x|`$ ein Eigenvektor.
Ist $`|x|_i=0`$, so erzwingt $`0=(\lambda-M_{ii})|x|_i=\sum_{j\ne i}M_{ij}|x|_j`$, dass $`|x|_j=0`$ für alle Nachbarn $`j`$ von
$`i`$ gilt, und wegen des Zusammenhangs ist $`x=0`$. Also hat kein von null verschiedener Eigenvektor einen
Nulleintrag; ein Eigenraum der Dimension $`\ge2`$ enthält aber einen von null verschiedenen Vektor mit $`x_1=0`$. Daher
ist $`\lambda`$ einfach, und $`|x|>0`$ ist ein Eigenvektor. $`\square`$

Wir beweisen die Vorzeichenregel für jede bipartite Anordnung der Kopplungen; die Kette ist ein Spezialfall.

<a id="r-9-7"></a>

**Satz 9.7 (Marshallsche Vorzeichenregel).** Die Plätze $`\{0,\dots,L-1\}`$ seien in zwei Mengen $`A\sqcup B`$ zerlegt,
und es sei
```math
H=\sum_{i<j}\Big[\tfrac{J_{ij}}2\big(s^+_is^-_j+s^-_is^+_j\big)+K_{ij}\,s^z_is^z_j\Big]+\sum_ih_is^z_i
```
mit reellen $`K_{ij},h_i`$ und reellen $`J_{ij}`$, sodass $`J_{ij}\ge0`$ gilt, wenn $`i,j`$ in verschiedenen Mengen liegen,
$`J_{ij}\le0`$, wenn sie in derselben Menge liegen, und der Graph $`G`$ auf den Plätzen mit einer Kante $`i\sim j`$ genau
dann, wenn $`J_{ij}\ne0`$, zusammenhängend ist. Sei $`N_B(b)`$ die Anzahl der Spins nach oben von $`b`$ auf $`B`$. In jedem
Sektor $`N`$ ist der kleinste Eigenwert von $`H_N`$ einfach, und sein Eigenvektor kann so gewählt werden, dass
$`(-1)^{N_B(b)}\psi(b)>0`$ für jede Konfiguration $`b`$ mit $`N`$ Spins nach oben gilt.

Insbesondere gilt dies für die XXZ-Kette mit $`J>0`$ und beliebigen $`\Delta,h_i`$, wenn die Kette offen oder periodisch
mit geradem $`L`$ ist, mit $`A`$ den geraden und $`B`$ den ungeraden Plätzen: Dann ist $`J_{ij}=J`$, $`K_{ij}=J\Delta`$ auf den
Bindungen und $`0`$ sonst, jede Bindung verbindet einen geraden mit einem ungeraden Platz (beim Ring, weil $`L`$ gerade
ist), und $`G`$ enthält den Pfad $`0-1-\cdots-(L-1)`$.

*Beweis.* Sei $`V|b\rangle=(-1)^{N_B(b)}|b\rangle`$; $`V`$ ist diagonal, $`V=V^*=V^{-1}`$. Der Operator $`s^+_is^-_j`$ bewegt einen
Spin nach oben von $`j`$ nach $`i`$; er ändert $`N_B`$ um $`\pm1`$, wenn $`i,j`$ in verschiedenen Mengen liegen, und lässt es
sonst unverändert, also ist entsprechend $`Vs^+_is^-_jV=\mp s^+_is^-_j`$, und ebenso für $`s^-_is^+_j`$. Zwei verschiedene
Konfigurationen $`b,b'`$ sind genau dann durch einen solchen Term verbunden, wenn $`b'`$ aus $`b`$ durch Vertauschen der
entgegengesetzten Spins auf zwei Plätzen $`i,j`$ entsteht, und dann sind $`i,j`$ durch $`b,b'`$ bestimmt (die Plätze, an
denen sie sich unterscheiden). Daher hat $`M=-VH_NV`$ die Nichtdiagonaleinträge $`J_{ij}/2`$ (verschiedene Mengen) bzw.
$`-J_{ij}/2`$ (dieselbe Menge), also in beiden Fällen $`|J_{ij}|/2\ge0`$, positiv genau dann, wenn $`i\sim j`$ in $`G`$, und
$`0`$ zwischen Konfigurationen, die sich anders unterscheiden. Der Graph von $`M`$ ist zusammenhängend. Denn man
schreibe $`(i\,j)`$ für die Vertauschung der Spins auf den Plätzen $`i,j`$; für verschiedene $`a,b,c`$ ist
$`(a\,c)=(a\,b)(b\,c)(a\,b)`$, also ist durch Induktion entlang eines Pfades in $`G`$ von $`i`$ nach $`j`$ jedes $`(i\,j)`$ eine
Komposition von Vertauschungen entlang Kanten von $`G`$, und jede solche Vertauschung fixiert entweder eine
Konfiguration (gleiche Spins) oder läuft entlang einer Kante des Graphen von $`M`$. Haben $`b\ne b'`$ je $`N`$ Spins nach
oben, so gibt es Plätze $`i`$, oben in $`b`$ und unten in $`b'`$, und $`j`$, unten in $`b`$ und oben in $`b'`$; $`(i\,j)b`$
unterscheidet sich von $`b'`$ an zwei Plätzen weniger als $`b`$, also führen endlich viele solche Schritte von $`b`$ zu
$`b'`$. Nach Lemma 9.6 ist der größte Eigenwert von $`M`$ – das Negative des kleinsten von $`H_N`$ – einfach, mit einem
Eigenvektor $`\varphi>0`$; man nehme $`\psi=V\varphi`$. $`\square`$

<a id="r-9-8"></a>

**Korollar 9.8 (Impuls des Grundzustands).** Auf einer periodischen Kette mit geradem $`L`$, $`J>0`$ und homogenem Feld
hat der Grundzustand des Sektors $`N`$ den Impuls $`k=0`$, wenn $`N`$ gerade ist, und $`k=\pi`$, wenn $`N`$ ungerade ist.

*Beweis.* $`T`$ kommutiert mit $`H_N`$, und der Grundzustand $`\psi=V\varphi`$ ist einfach, also gilt $`T\psi=e^{ik}\psi`$.
Nun ist $`(T\psi)(b)=\psi(T^{-1}b)`$, und $`T^{-1}b`$ hat auf Platz $`i`$ den Spin $`b_{i+1}`$; da $`L`$ gerade ist, ist $`i+1`$
für ungerades $`i`$ gerade, also $`N_B(T^{-1}b)=N-N_B(b)`$ und
$`(VT\psi)(b)=(-1)^{N_B(b)}(-1)^{N-N_B(b)}\varphi(T^{-1}b)=(-1)^N\varphi(T^{-1}b)`$. Außerdem ist $`VT\psi=e^{ik}\varphi`$.
Sowohl $`\varphi`$ als auch $`\varphi\circ T^{-1}`$ sind positiv, also $`e^{ik}=(-1)^N`$. $`\square`$

*Rechnung:* Satz 9.7 für offene und periodische Ketten mit Zufallsfeldern und Korollar 9.8 für $`L=12`$ — [Tutorial, Abschnitt 9.2](../../tutorial/de/09_spin_chains.md#s-9-2).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). In einem Antiferromagneten legt
> das Gitter die Vorzeichen des Grundzustands fest: Einen Spin über eine Bindung zu bewegen, ändert das Vorzeichen.
> Physiker nutzen das, um Quanten-Monte-Carlo für bipartite Antiferromagneten ohne Vorzeichenproblem zu rechnen;
> bei frustrierten Ketten, in denen eine antiferromagnetische Kopplung zwei Plätze derselben Menge verbindet, ist
> die Voraussetzung von Satz 9.7 verletzt.

<a id="s-9-3"></a>

## 9.3 $`SU(2)`$ bei $`\Delta=1`$

Für $`\Delta=1`$ und $`h=0`$ kommutiert $`H=J\sum\mathbf s_i\cdot\mathbf s_j`$ mit $`S^x,S^y,S^z`$ (Proposition 4.15), also
ist $`H\cong\bigoplus_S1_{2S+1}\otimes H_S`$ mit Blöcken $`H_S\in M_{m_S}`$, $`m_S=\binom L{L/2-S}-\binom L{L/2-S-1}`$
(Proposition 4.13), und das Modul `su2` baut $`H_S`$ auf den Vektoren höchsten Gewichts (Proposition 4.15(b)).
Dasselbe gilt für jedes $`H=\sum_{i<j}J_{ij}\,\mathbf s_i\cdot\mathbf s_j`$ mit reellen $`J_{ij}`$: Nach Proposition 4.15(a) liegt
jedes $`\mathbf s_i\cdot\mathbf s_j`$ in $`\{S_a\}'`$. Sei $`E_0(S)`$ der kleinste Eigenwert von $`H_S`$, für die auftretenden Spins $`S`$.

*In der Bibliothek:* Die Blöcke $`H_S`$ baut `SU2Chain` — [Tutorial, Abschnitt 9.3](../../tutorial/de/09_spin_chains.md#s-9-3).

<a id="r-9-9"></a>

**Satz 9.9 (Lieb-Mattis).** Sei $`H=\sum_{i<j}J_{ij}\,\mathbf s_i\cdot\mathbf s_j`$ mit $`A,B,J_{ij}`$ wie in Satz 9.7 (also
$`K_{ij}=J_{ij}`$, $`h_i=0`$), $`n_A=|A|\ge n_B=|B|`$ und $`S_*=\frac{n_A-n_B}2`$.
(a) Für $`S_*\le S\le\frac L2`$ (mit ganzzahligem $`S-\frac L2`$) hat der Grundzustand des Sektors $`S^z=S`$ den Gesamtspin
$`S`$; für $`0\le M\le S_*`$ (mit ganzzahligem $`M-\frac L2`$) hat der Grundzustand des Sektors $`S^z=M`$ den Gesamtspin $`S_*`$.
(b) Folglich ist $`E_0(S)`$ für $`S\ge S_*`$ ein einfacher Eigenwert von $`H_S`$, $`E_0(S_*)<E_0(S_*+1)<\dots<E_0(\frac L2)`$,
und $`E_0(S)>E_0(S_*)`$ für $`S<S_*`$: Der Grundraum von $`H`$ ist ein einziges Multiplett, vom Spin $`S_*`$.
Für die Heisenberg-Kette mit $`J>0`$, offen oder periodisch mit geradem $`L`$, sind $`A`$ und $`B`$ die geraden und
ungeraden Plätze, $`S_*=0`$ für gerades $`L`$ und $`S_*=\frac12`$ für ungerades $`L`$ (offene Kette), und
$`E_0(S_*)<E_0(S_*+1)<\dots<E_0(\frac L2)`$.

*Beweis.* (a) Sei $`M\ge0`$ mit $`M-\frac L2\in\mathbb Z`$, $`s=\max(M,S_*)`$, und sei $`\psi`$ der Grundzustand des Sektors
$`S^z=M`$, einfach nach Satz 9.7: $`\psi=V\varphi`$ mit $`\varphi>0`$. $`\mathbf S^2`$ kommutiert mit $`H`$ und mit $`S^z`$, erhält also
den eindimensionalen Grundraum des Sektors: $`\psi`$ hat einen Gesamtspin $`S'`$. Wir konstruieren einen Vektor $`\chi`$
des Sektors mit $`\mathbf S^2\chi=s(s+1)\chi`$ und $`\langle\psi,\chi\rangle\ne0`$; da Eigenvektoren von $`\mathbf S^2`$ zu
verschiedenen Eigenwerten orthogonal sind, folgt $`S'=s`$. Man setze $`p=\frac L2-s`$; dann ist $`0\le p\le n_B`$, denn
$`s\ge S_*`$ bedeutet $`p\le\frac{n_A+n_B}2-\frac{n_A-n_B}2`$. Man wähle verschiedene $`a_1,\dots,a_p\in A`$ und verschiedene
$`b_1,\dots,b_p\in B`$, sei $`R`$ die Menge der übrigen $`2s`$ Plätze, und sei
```math
\chi_0=\bigotimes_{l=1}^p\tfrac1{\sqrt2}\big(|{\uparrow}_{a_l}{\downarrow}_{b_l}\rangle-|{\downarrow}_{a_l}{\uparrow}_{b_l}\rangle\big)\otimes|{\uparrow\cdots\uparrow}\rangle_R ,
\qquad \chi=(S^-)^{s-M}\chi_0 .
```
Jedes Singulett wird von $`s^\pm_{a_l}+s^\pm_{b_l}`$ annihiliert und die Spins nach oben auf $`R`$ von $`s^+_i`$, also ist
$`S^+\chi_0=0`$; mit $`S^z\chi_0=s\chi_0`$ ergibt das $`\mathbf S^2\chi_0=(S^-S^++(S^z)^2+S^z)\chi_0=s(s+1)\chi_0`$
(Abschnitt 4.5), und $`\mathbf S^2\chi=s(s+1)\chi`$, da $`\mathbf S^2`$ mit $`S^-`$ kommutiert. Schreibt man
$`S^-=\sum_l(s^-_{a_l}+s^-_{b_l})+S^-_R`$, so kommutieren die Summanden, und die ersten annihilieren die Singuletts;
also ist $`\chi`$ das Produkt der Singuletts mit $`(S^-_R)^{s-M}|{\uparrow\cdots\uparrow}\rangle_R=(s-M)!\sum_D|D\rangle`$, der
Summe über die Konfigurationen von $`R`$ mit genau $`s-M`$ Spins nach unten (diese existieren, da $`s-M\le2s=|R|`$);
somit ist $`\chi\ne0`$, $`S^z\chi=M\chi`$, und die Koeffizienten von $`\chi`$ sind reell. Sein Vorzeichen: Ein Singulett
trägt $`+`$ bei, wenn sein Spin nach oben auf $`A`$ sitzt, und $`-`$, wenn er auf $`B`$ sitzt, und der $`R`$-Anteil ist
positiv. Ist $`M=s`$, so zeigen alle Plätze von $`R`$ nach oben; ist $`M<s`$, so ist $`s=S_*`$, $`p=n_B`$ und $`R\subseteq A`$. In
beiden Fällen ist $`N_B(b)`$ die Anzahl der Singuletts mit dem Spin nach oben auf $`B`$ plus eine Konstante, $`|R\cap B|`$
(falls $`M=s`$) bzw. $`0`$, also hat $`(-1)^{N_B(b)}\chi(b)`$ für alle $`b`$ mit $`\chi(b)\ne0`$ dasselbe Vorzeichen. Daher ist
$`\langle\psi,\chi\rangle=\langle\varphi,V\chi\rangle=\pm\sum_b\varphi(b)|\chi(b)|\ne0`$.

(b) Nach Proposition 4.15(a), Korollar 4.12 und Satz 4.11 trifft ein Multiplett vom Spin $`S'`$ den Sektor $`S^z=M`$ in
einem Vektor, wenn $`S'\ge|M|`$, und sonst gar nicht; also ist das Spektrum von $`H`$ auf dem Sektor $`S^z=M\ge0`$ die
Vereinigung der Spektren der $`H_{S'}`$, $`S'\ge M`$, mit Vielfachheiten. Nach (a) ist sein kleinster Eigenwert einfach
und gehört zum Spin $`s=\max(M,S_*)`$: Also ist $`E_0(s)`$ einfach in $`H_s`$, und $`E_0(s)<E_0(S')`$ für alle $`S'\ge M`$,
$`S'\ne s`$. Für $`M=S\ge S_*`$ ergibt das die Einfachheit und $`E_0(S)<E_0(S')`$ für $`S'>S`$; für $`M=S<S_*`$ ergibt es
$`E_0(S_*)<E_0(S)`$. Also ist $`E_0(S_*)`$ der kleinste Eigenwert von $`H`$, angenommen nur in $`H_{S_*}`$, wo er einfach ist.
Bei der Kette übertrifft die Anzahl der geraden Plätze die der ungeraden um $`0`$ ($`L`$ gerade) oder $`1`$ ($`L`$
ungerade). $`\square`$

Satz 9.7 geht auf W. Marshall (*Proc. Roy. Soc. A* 232 (1955) 48) zurück und Satz 9.9 auf E. Lieb und D. Mattis
(*J. Math. Phys.* 3 (1962) 749). Die **Singulett-Triplett-Lücke** $`E_0(1)-E_0(0)`$ des geraden Heisenberg-Rings ist
nach Satz 9.9 positiv.

*Rechnung:* Satz 9.9 an den $`U(1)`$-Sektoren für $`L=12`$, die Singulett-Triplett-Lücke sowie die Sätze 9.7 und 9.9 auf einem bipartiten Graphen, der keine Kette ist — [Tutorial, Abschnitt 9.3](../../tutorial/de/09_spin_chains.md#s-9-3).

**Korrelationen.** Auf jedem geraden Heisenberg-Ring wechseln die Grundzustandskorrelationen
$`\langle\mathbf s_0\cdot\mathbf s_r\rangle`$ das Vorzeichen; die Sätze 9.7 und 9.9 erzwingen dies. Sei $`\psi`$ der reelle
normierte Grundzustand, ein Singulett (Satz 9.9) und daher der Grundzustand des Sektors $`S^z=0`$. Aus
$`\mathbf S^2\psi=0`$ folgt $`\|S^\pm\psi\|^2=\langle\psi,S^\mp S^\pm\psi\rangle=\langle\psi,(\mathbf S^2-(S^z)^2\mp S^z)\psi\rangle=0`$, also
$`S^\pm\psi=0`$ und $`\langle\psi,[S^+,[S^-,X]]\psi\rangle=0`$ für jedes $`X`$ (jeder der vier Terme hat $`S^\pm`$ neben $`\psi`$). Für
$`X=s^z_0s^z_r`$, $`r\ne0`$, liefern die Relationen $`[s^\pm,s^z]=\mp s^\pm`$, $`[s^+,s^-]=2s^z`$ die Identitäten
$`[S^-,X]=s^-_0s^z_r+s^z_0s^-_r`$ und $`[S^+,[S^-,X]]=4s^z_0s^z_r-s^+_0s^-_r-s^-_0s^+_r`$. Daher ist
$`\langle s^z_0s^z_r\rangle=\frac14\langle s^+_0s^-_r+s^-_0s^+_r\rangle`$ und
```math
\langle\mathbf s_0\cdot\mathbf s_r\rangle=\langle s^z_0s^z_r\rangle+\tfrac12\langle s^+_0s^-_r+s^-_0s^+_r\rangle
=\tfrac34\sum_{b_0\ne b_r}\psi(b)\psi(b') ,
```
wobei $`b'`$ aus $`b`$ durch Vertauschen der Spins auf den Plätzen $`0`$ und $`r`$ entsteht. Nach Satz 9.7 hat
$`\psi(b)\psi(b')`$ das Vorzeichen von $`(-1)^{N_B(b)+N_B(b')}`$, und $`N_B(b')-N_B(b)`$ ist für ungerades $`r`$ gleich $`\pm1`$ und
für gerades $`r`$ gleich $`0`$; die Summe ist nicht leer, also ist $`\langle\mathbf s_0\cdot\mathbf s_r\rangle`$ für $`0<r<L`$ von
null verschieden mit dem Vorzeichen $`(-1)^r`$.

*Rechnung:* Die Korrelationen $`\langle\mathbf s_0\cdot\mathbf s_r\rangle`$ für $`L=12`$, mit `SU2Chain.spin_correlation`, das den Satz von Wigner-Eckart (Satz B.16) verwendet — [Tutorial, Abschnitt 9.3](../../tutorial/de/09_spin_chains.md#s-9-3).

<a id="s-9-4"></a>

## 9.4 Exakte Prüfungen: freie Fermionen und der Bethe-Ansatz

**Jordan-Wigner.** Auf einer offenen Kette setze man $`c_j=\big(\prod_{l<j}(-1)^{n_l}\big)s^-_j`$.

<a id="r-9-10"></a>

**Satz 9.10 (die XX-Kette ist frei).** (a) $`\{c_i,c_j^*\}=\delta_{ij}`$, $`\{c_i,c_j\}=0`$ und $`c_j^*c_j=n_j`$.
(b) Die offene XX-Kette ist $`H=\frac J2\sum_{j<L-1}(s^+_js^-_{j+1}+s^-_js^+_{j+1})=\sum_{i,j}h_{ij}c_i^*c_j`$ mit
$`h=\frac J2A_L`$. (c) Seien $`u_k`$ orthonormale reelle Eigenvektoren von $`h`$ zu den Eigenwerten
$`\varepsilon_k=J\cos\frac{\pi k}{L+1}`$ (Lemma 9.3) und $`b_k=\sum_ju_k(j)c_j`$. Für $`K=\{k_1<\dots<k_N\}`$ bilden die
Vektoren $`b_K^*|0\rangle=b_{k_1}^*\cdots b_{k_N}^*|0\rangle`$ ($`|0\rangle`$: alle Spins nach unten) eine orthonormale
Eigenbasis des Sektors $`N`$, mit Energien $`\sum_{k\in K}\varepsilon_k`$.

*Beweis.* (a) Operatoren auf verschiedenen Plätzen kommutieren, und auf einem Platz antikommutiert $`(-1)^n`$ mit
$`s^\pm`$. Für $`i<j`$ treffen $`c_i`$ und $`c_j^*`$ (oder $`c_j`$) auf Platz $`i`$ als $`s^-_i`$ und $`(-1)^{n_i}`$ aufeinander –
eine Antivertauschung – und kommutieren ansonsten, also antikommutieren sie. Für $`i=j`$ quadrieren sich die Strings
zu $`1`$, und $`c_jc_j^*+c_j^*c_j=s^-_js^+_j+s^+_js^-_j=1`$, $`c_j^*c_j=s^+_js^-_j=n_j`$, $`c_j^2=(s^-_j)^2=0`$; der Fall $`i>j`$
folgt durch Vertauschen der Rollen von $`i,j`$ (und Übergang zu den Adjungierten). (b)
$`c_j^*c_{j+1}=s^+_j(-1)^{n_j}s^-_{j+1}=s^+_js^-_{j+1}`$, da $`s^+_j\neq0`$ nur auf $`n_j=0`$; dann geht man zu den
Adjungierten über. (c) $`\{b_k,b_l^*\}=\sum_ju_k(j)u_l(j)=\delta_{kl}`$ und $`\{b_k,b_l\}=0`$, und
$`H=\sum_k\varepsilon_kb_k^*b_k`$ nach der Spektralzerlegung von $`h`$. Aus den CAR folgt $`[b_l^*b_l,b_k^*]=\delta_{lk}b_k^*`$,
also $`[H,b_k^*]=\varepsilon_kb_k^*`$; mit $`H|0\rangle=0`$ ergibt das $`Hb_K^*|0\rangle=(\sum_K\varepsilon_k)b_K^*|0\rangle`$. Da
$`b_k|0\rangle=0`$ (jedes $`c_j`$ senkt einen Spin ab), liefert das Durchschieben der Vernichter nach rechts mit den CAR
$`\langle b_{K'}^*0,b_K^*0\rangle=\delta_{KK'}`$ durch Induktion über $`N`$. Jedes $`c_j^*`$ erhöht $`N`$ um eins, also liegt
$`b_K^*|0\rangle`$ im Sektor $`|K|`$: $`\binom LN`$ orthonormale Vektoren in einem Raum dieser Dimension. $`\square`$

*Rechnung:* Alle $`1024`$ Niveaus der offenen XX-Kette mit $`L=10`$, im Vergleich mit Satz 9.10 — [Tutorial, Abschnitt 9.4](../../tutorial/de/09_spin_chains.md#s-9-4).

**Magnonen.** Auf dem Ring mit $`h=0`$ hat der vollständig polarisierte Zustand $`|0\rangle`$ die Energie
$`E_F=J\Delta L/4`$. Ein Spin nach oben bei $`x`$ bricht zwei Bindungen auf:
$`H|x\rangle=(E_F-J\Delta)|x\rangle+\frac J2(|x-1\rangle+|x+1\rangle)`$ ($`L\ge3`$).

<a id="r-9-11"></a>

**Proposition 9.11 (ein Magnon).** Im Sektor $`N=1`$ des Rings erfüllt $`|k\rangle=\sum_xe^{-ikx}|x\rangle`$, $`k=2\pi m/L`$,
die Gleichungen $`T|k\rangle=e^{ik}|k\rangle`$ und $`H|k\rangle=\big(E_F+J(\cos k-\Delta)\big)|k\rangle`$.

*Beweis.* $`\sum_xe^{-ikx}(|x-1\rangle+|x+1\rangle)=\sum_y(e^{-ik(y+1)}+e^{-ik(y-1)})|y\rangle=2\cos k\,|k\rangle`$ unter
Verwendung von $`e^{ikL}=1`$; dieselbe Umindizierung liefert $`T|k\rangle=e^{ik}|k\rangle`$. $`\square`$

Zwei Magnonen streuen aneinander. Wir schreiben $`|x_1,x_2\rangle`$, $`x_1<x_2`$, für die Spins nach oben bei $`x_1,x_2`$; ein
Paar von Nachbarn bricht zwei Bindungen auf, ein getrenntes Paar vier.

<a id="r-9-12"></a>

**Proposition 9.12 (zwei Magnonen; Bethe-Gleichungen).** Seien $`k_1,k_2\in\mathbb C`$, $`K=k_1+k_2`$ und
$`(A_{12},A_{21})\neq0`$ mit (i) $`A_{12}(1+e^{iK}-2\Delta e^{ik_2})+A_{21}(1+e^{iK}-2\Delta e^{ik_1})=0`$ und
(ii) $`A_{21}=e^{ik_2L}A_{12}`$, $`A_{12}=e^{ik_1L}A_{21}`$. Seien
$`\psi(x_1,x_2)=A_{12}e^{i(k_1x_1+k_2x_2)}+A_{21}e^{i(k_2x_1+k_1x_2)}`$ auf $`\mathbb Z^2`$ und
$`\Psi=\sum_{0\le x_1<x_2<L}\psi(x_1,x_2)|x_1,x_2\rangle`$. Ist $`\Psi\neq0`$, so gilt
$`H\Psi=\big(E_F+J(\cos k_1+\cos k_2-2\Delta)\big)\Psi`$. Verschwinden die Nenner nicht, so besagen (i) und (ii)
```math
e^{ik_1L}=-\frac{1+e^{iK}-2\Delta e^{ik_1}}{1+e^{iK}-2\Delta e^{ik_2}},\qquad
e^{ik_2L}=-\frac{1+e^{iK}-2\Delta e^{ik_2}}{1+e^{iK}-2\Delta e^{ik_1}} .
```

*Beweis.* Sei $`E=E_F+J(\cos k_1+\cos k_2-2\Delta)`$. (1) Für alle $`(x_1,x_2)\in\mathbb Z^2`$ liefert jede ebene Welle Werte
$`\psi(x_1\pm1,x_2)`$ und $`\psi(x_1,x_2\pm1)`$ mit der Summe $`2(\cos k_1+\cos k_2)\psi(x_1,x_2)`$. (2) Durch Einsetzen ist
$`\psi(x,x)+\psi(x+1,x+1)-2\Delta\psi(x,x+1)`$ gleich $`e^{iKx}`$ mal der linken Seite von (i), verschwindet also für alle
$`x`$. (3) Nach (ii) ist
$`\psi(x_2,x_1+L)=A_{12}e^{ik_2L}e^{i(k_2x_1+k_1x_2)}+A_{21}e^{ik_1L}e^{i(k_1x_1+k_2x_2)}=\psi(x_1,x_2)`$, und zweimalige
Anwendung ergibt $`\psi(x_1+L,x_2+L)=\psi(x_1,x_2)`$. Die Repräsentanten einer Konfiguration $`\{y_1<y_2\}`$ in
$`D=\{x_1<x_2<x_1+L\}`$ sind $`(y_1+nL,y_2+nL)`$ und $`(y_2+nL,y_1+(n+1)L)`$, also ist $`\psi`$ eine Funktion der
Konfigurationen. Sei $`(x_1,x_2)\in D`$. Ist $`2\le x_2-x_1\le L-2`$, so bleiben alle vier Sprünge in $`D`$, und der
Diagonalterm ist $`E_F-2J\Delta`$, also $`(H\Psi)(x_1,x_2)=(E_F-2J\Delta)\psi+\frac J2\cdot2(\cos k_1+\cos k_2)\psi=E\psi`$
nach (1). Ist $`x_2=x_1+1`$, so sind die Sprünge $`(x_1+1,x_2)`$, $`(x_1,x_2-1)`$ blockiert, und der Diagonalterm ist
$`E_F-J\Delta`$; nach (1) und (2) ist
$`(H\Psi)(x_1,x_2)=(E_F-J\Delta)\psi+\frac J2\big(2(\cos k_1+\cos k_2)\psi-\psi(x_1,x_1)-\psi(x_1+1,x_1+1)\big)=E\psi`$.
Ist $`x_2=x_1+L-1`$, so verwende man den Repräsentanten $`(x_2,x_1+L)`$, der vom vorigen Typ ist ($`L\ge3`$). $`\square`$

Proposition 9.12 sagt nur, dass Lösungen Eigenvektoren liefern; sie sagt nicht, dass alle Eigenvektoren so
entstehen. Wir schreiben $`k_{1,2}=K/2\pm q`$ und $`c=\cos\frac K2`$. Dann ist
$`1+e^{iK}-2\Delta e^{ik_{1,2}}=2e^{iK/2}(c-\Delta e^{\pm iq})`$, also lösen $`A_{12}=c-\Delta e^{iq}`$,
$`A_{21}=-(c-\Delta e^{-iq})`$ die Gleichung (i). Nach (ii) und $`(A_{12},A_{21})\ne0`$ sind beide Koeffizienten von null
verschieden, also ergibt Multiplikation der beiden Gleichungen (ii) $`e^{iKL}=1`$,
$`K=2\pi n/L`$, und dann ist (ii) äquivalent zu $`F(q)=(-1)^ne^{iqL}(c-\Delta e^{-iq})+c-\Delta e^{iq}=0`$, einer Gleichung
in einer Variablen. Wir suchen ihre Nullstellen $`q\in(0,\pi)`$, wo $`e^{-iqL/2}F(q)`$ reell oder rein imaginär ist,
weil $`c-\Delta e^{\pm iq}`$ zueinander komplex konjugiert sind (Streuzustände), und auf den Geraden $`q=i\kappa`$,
$`\pi+i\kappa`$, $`\kappa>0`$, wo $`F`$ reell ist (gebundene Zustände).

Für gerades $`L`$ und $`\Delta\ne0`$ sind die Bethe-Vektoren nicht vollständig. Ist $`\Delta\ne0`$, so erfüllt $`\Phi=\sum_x(-1)^x|x,x+1\rangle`$ die
Gleichungen $`T\Phi=-\Phi`$ und $`H\Phi=(E_F-J\Delta)\Phi`$ (die Sprünge nach $`|x-1,x+1\rangle`$ und $`|x,x+2\rangle`$ heben sich
paarweise auf). Ein Bethe-Vektor erfüllt $`T\Psi=e^{-iK}\Psi`$, da $`\psi(x_1-1,x_2-1)=e^{-iK}\psi(x_1,x_2)`$. Ist
$`e^{-iK}\ne-1`$, so ist $`\Psi\perp\Phi`$, weil beide Eigenvektoren des unitären $`T`$ zu verschiedenen Eigenwerten sind;
ist $`e^{iK}=-1`$, so gilt $`\cos k_1+\cos k_2=2\cos\frac K2\cos\frac{k_1-k_2}2=0`$, die Energie von $`\Psi`$ ist
$`E_F-2J\Delta\ne E_F-J\Delta`$, und $`\Psi\perp\Phi`$ als Eigenvektoren von $`H`$ zu verschiedenen Eigenwerten. Also ist $`\Phi`$
zu allen Bethe-Vektoren aus Proposition 9.12 orthogonal, die daher den Sektor $`N=2`$ eines geraden Rings nie
aufspannen. Die Zwei-Magnonen-Rechnung ist der Fall $`N=2`$
des Ansatzes von H. Bethe (*Z. Phys.* 71 (1931) 205); keinen anderen Fall verwenden wir.

*Rechnung:* Die Nullstellen von $`F`$ für $`L=9`$, $`\Delta=0.6`$, die alle $`36`$ Zustände des Sektors $`N=2`$ liefern, und der Vektor $`\Phi`$ für $`L=8`$ — [Tutorial, Abschnitt 9.4](../../tutorial/de/09_spin_chains.md#s-9-4).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Ein Magnon ist eine Spinwelle:
> ein umgeklappter Spin, delokalisiert über den Ring, mit Impuls $`k`$. Zwei Magnonen streuen mit einer
> Phasenverschiebung (reelles $`q`$) oder binden aneinander (komplexes $`q`$).

<a id="s-9-5"></a>

## 9.5 Lanczos- und Krylov-Verfahren

Der Sektor $`N=L/2`$ hat die Dimension $`12870`$ für $`L=16`$ und $`184756`$ für $`L=20`$, aber $`H`$ ist dünn besetzt ($`O(L)`$
von null verschiedene Einträge pro Zeile, `xxz_sparse`), und Krylov-Verfahren brauchen nur Produkte $`Hv`$. Sei $`H`$
hermitesch von der Größe $`n`$ und $`v\ne0`$. Der **Krylov-Raum** ist
$`\mathcal K_m=\mathrm{span}\{v,Hv,\dots,H^{m-1}v\}`$.

<a id="r-9-13"></a>

**Proposition 9.13 (Lanczos).** Sei $`\dim\mathcal K_m=m`$. Sei $`q_1=v/\|v\|`$, und für $`j<m`$ sei $`q_{j+1}`$ die normierte
Komponente von $`Hq_j`$ orthogonal zu $`\mathcal K_j`$. Dann ist $`q_1,\dots,q_m`$ eine Orthonormalbasis von $`\mathcal K_m`$;
mit $`Q_m=(q_1\cdots q_m)`$ ist $`T_m=Q_m^*HQ_m`$ reell symmetrisch und tridiagonal mit Diagonale
$`\alpha_j=\langle q_j,Hq_j\rangle`$ und Nebendiagonale $`\beta_j=\langle q_{j+1},Hq_j\rangle>0`$; und $`HQ_m=Q_mT_m+r_me_m^T`$, wobei
$`r_m\perp\mathcal K_m`$ genau dann verschwindet, wenn $`\mathcal K_{m+1}=\mathcal K_m`$ (andernfalls $`r_m=\beta_mq_{m+1}`$,
$`q_{m+1}=r_m/\|r_m\|`$).

*Beweis.* Es gelte $`\mathrm{span}\{q_1,\dots,q_j\}=\mathcal K_j`$ (wahr für $`j=1`$). Dann ist
$`H\mathcal K_j\subseteq\mathcal K_j+\mathbb CHq_j`$, da $`Hq_i\in\mathcal K_{i+1}`$; also
$`\mathcal K_{j+1}=\mathcal K_j+\mathbb CH^jv\subseteq\mathcal K_j+\mathbb CHq_j\subseteq\mathcal K_{j+1}`$ und
$`\mathcal K_{j+1}=\mathcal K_j+\mathbb CHq_j`$. Für $`j<m`$ wächst die Dimension, also ist die zu $`\mathcal K_j`$ orthogonale
Komponente von $`Hq_j`$ von null verschieden, und die Induktion läuft weiter; für $`j=m`$ verschwindet diese
Komponente, $`r_m`$, genau dann, wenn $`\mathcal K_{m+1}=\mathcal K_m`$. Für $`i>j+1`$ ist $`q_i\perp\mathcal K_{j+1}\ni Hq_j`$, also
$`(T_m)_{ij}=0`$, und $`(T_m)_{ji}=0`$, da $`T_m`$ hermitesch ist. $`\alpha_j`$ ist reell, und $`\beta_j=\langle q_{j+1},Hq_j\rangle`$
ist die Norm der Komponente, $`>0`$. Also ist $`Hq_j=\beta_{j-1}q_{j-1}+\alpha_jq_j+\beta_jq_{j+1}`$ für $`j<m`$ und
$`Hq_m=\beta_{m-1}q_{m-1}+\alpha_mq_m+r_m`$. $`\square`$

Die Eigenwerte $`\theta^{(m)}_1\le\dots\le\theta^{(m)}_m`$ von $`T_m`$ sind die **Ritz-Werte**.

<a id="r-9-14"></a>

**Satz 9.14 (monotone Konvergenz).** Seien $`\lambda_1\le\dots\le\lambda_n`$ die Eigenwerte von $`H`$ und
$`\dim\mathcal K_m=m`$. Dann gilt $`\theta^{(m)}_j\ge\lambda_j`$, und $`\theta^{(m)}_j\ge\theta^{(m+1)}_j`$, falls
$`\dim\mathcal K_{m+1}=m+1`$ ($`j\le m`$). Hat $`v`$ eine von null verschiedene Komponente im Eigenraum von $`\lambda_1`$ und
ist $`\mathcal K_{m+1}=\mathcal K_m`$, so ist $`\theta^{(m)}_1=\lambda_1`$.

*Beweis.* Nach Korollar A.9 (Courant–Fischer für Kompressionen, Satz A.8), angewandt auf $`Q=Q_m`$, ist
$`\theta^{(m)}_j=\min_{U}\max_{x\in U,\|x\|=1}\langle x,Hx\rangle`$ über die $`j`$-dimensionalen $`U\subseteq\mathcal K_m=\mathrm{ran}\,Q_m`$.
Dieselbe Formel über $`U\subseteq\mathcal K_{m+1}`$ liefert $`\theta^{(m+1)}_j`$ und über $`U\subseteq\mathbb C^n`$ liefert
$`\lambda_j`$ (Satz A.8); Minimieren über
mehr Unterräume ergibt einen kleineren Wert. Ist $`\mathcal K_{m+1}=\mathcal K_m`$, so ist $`\mathcal K_m`$ $`H`$-invariant und
enthält $`p(H)v`$ für jedes Polynom $`p`$, insbesondere $`P_1v\ne0`$ ($`P_1`$ die Spektralprojektion zu $`\lambda_1`$, nach
Lagrange-Interpolation ein Polynom in $`H`$, wie in Satz 1.5); also ist $`\theta^{(m)}_1\le\langle x,Hx\rangle=\lambda_1`$ bei
$`x=P_1v/\|P_1v\|\in\mathcal K_m`$. $`\square`$

Ein kleines Residuum garantiert einen nahen Eigenwert: Für hermitesches $`H`$, eine Zahl $`E`$ und $`\|\psi\|=1`$ ist
$`\|H\psi-E\psi\|^2=\sum_i|\lambda_i-E|^2|\langle u_i,\psi\rangle|^2\ge\min_i|\lambda_i-E|^2`$ ($`u_i`$ eine Orthonormalbasis aus
Eigenvektoren). Es garantiert nicht, dass dieser Eigenwert der kleinste ist. Für $`L=4`$ ist die Grundzustandsenergie
des Heisenberg-Rings exakt: Mit $`\mathbf S_A=\mathbf s_0+\mathbf s_2`$, $`\mathbf S_B=\mathbf s_1+\mathbf s_3`$ ist der
Ring $`H=J\,\mathbf S_A\cdot\mathbf S_B=\frac J2(\mathbf S^2-\mathbf S_A^2-\mathbf S_B^2)\ge\frac J2(0-2-2)=-2J`$, da $`\mathbf S^2\ge0`$ und
$`\mathbf S_A^2,\mathbf S_B^2`$ die Eigenwerte $`0,2`$ haben (Proposition 4.13 mit $`L=2`$); Gleichheit gilt auf einem Vektor vom
Spin $`0`$ im Produkt der beiden Tripletts, der nach Kapitel 4, Aufgabe 2 existiert ($`1\otimes1\cong2\oplus1\oplus0`$).

*Rechnung:* Satz 9.14 für die offene XX-Kette, mit `krylov.lanczos` auf dem dünn besetzten Sektor aus `SpinChain.xxz_sparse`, und Grundzustände des Heisenberg-Rings bis $`L=16`$ mit `krylov.ground_state` und seinem Residuum — [Tutorial, Abschnitt 9.5](../../tutorial/de/09_spin_chains.md#s-9-5).

**Zeitentwicklung.** Für $`e^{-i\tau H}v`$ projiziert man auf $`\mathcal K_m`$: $`w_\tau=\|v\|\,Q_me^{-i\tau T_m}e_1`$.

<a id="r-9-15"></a>

**Proposition 9.15.** Sei $`\dim\mathcal K_m=m`$. Für jedes $`\tau`$ gilt $`\|w_\tau\|=\|v\|`$ und
$`\langle w_\tau,Hw_\tau\rangle=\langle v,Hv\rangle`$. Ist $`\mathcal K_m`$ $`H`$-invariant, so ist $`w_\tau=e^{-i\tau H}v`$.

*Beweis.* $`Q_m`$ ist eine Isometrie und $`e^{-i\tau T_m}`$ unitär ($`T_m`$ ist reell symmetrisch), das liefert die Norm. Da
$`Q_m^*HQ_m=T_m`$ mit $`e^{-i\tau T_m}`$ kommutiert, ist
$`\langle w_\tau,Hw_\tau\rangle=\|v\|^2\langle e^{-i\tau T_m}e_1,T_me^{-i\tau T_m}e_1\rangle=\|v\|^2\langle e_1,T_me_1\rangle=\langle v,Hv\rangle`$. Ist
$`\mathcal K_m`$ invariant, so gilt $`HQ_m=Q_mT_m`$ (Proposition 9.13, $`r_m=0`$), also $`H^jQ_m=Q_mT_m^j`$ und, durch
Summation der Exponentialreihe (Proposition 8.2(1)), $`e^{-i\tau H}Q_m=Q_me^{-i\tau T_m}`$; man wende dies auf $`e_1`$ an, mit $`Q_me_1=v/\|v\|`$. $`\square`$

*Rechnung:* Proposition 9.15 für `krylov.evolve`, das schrittweise in Krylov-Räumen rechnet, bei einem Quench aus dem Néel-Zustand, im Vergleich mit dem dichten Propagator — [Tutorial, Abschnitt 9.5](../../tutorial/de/09_spin_chains.md#s-9-5).

<a id="s-9-6"></a>

## 9.6 Unordnung: Niveaustatistik und Verschränkung

Die **Heisenberg-Kette im Zufallsfeld** hat unabhängige, auf $`[-W,W]`$ gleichverteilte Felder $`h_i`$. Sie erhält
$`S^z`$ (Proposition 9.1), daher wird die Niveaustatistik pro Sektor $`N`$ berechnet (Satz 5.15 und Abschnitt 5.4: Das
Mischen von Sektoren täuscht unabhängige Niveaus vor). Unterscheiden sich zwei Felder, so ist
$`THT^*-H=\sum_i(h_{i-1}-h_i)s^z_i\ne0`$, Translationen sind also keine Symmetrie mehr.

*Rechnung:* Das mittlere Verhältnis $`\langle r\rangle`$ der Kette im Zufallsfeld, im Vergleich mit GOE-Matrizen und mit unabhängigen Niveaus (numerische Beobachtungen) — [Tutorial, Abschnitt 9.6](../../tutorial/de/09_spin_chains.md#s-9-6).

**Verschränkung freier Fermionen.** Für die Eigenzustände der offenen XX-Kette ist die Verschränkungsentropie der
ersten $`n`$ Plätze durch eine $`n\times n`$-Matrix bestimmt.

<a id="r-9-16"></a>

**Satz 9.16 (Peschels Formel).** Sei $`\Phi=b_{k_1}^*\cdots b_{k_N}^*|0\rangle`$ wie in Satz 9.10,
$`C_{ij}=\langle\Phi,c_i^*c_j\Phi\rangle`$ und $`C_A=(C_{ij})_{i,j<n}`$. Dann gilt $`0\le C_A\le1`$, und die
Verschränkungsentropie der Plätze $`A=\{0,\dots,n-1\}`$ ist, wobei $`\nu`$ die Eigenwerte von $`C_A`$ mit Vielfachheit
durchläuft,
```math
S_A=\sum_\nu \eta(\nu),\qquad \eta(\nu)=-\nu\log\nu-(1-\nu)\log(1-\nu)\quad(0\log0=0).
```

*Beweis.* Für einen reellen Vektor $`x`$ setze man $`b(x)=\sum_jx_jc_j`$; dann ist $`b(x)^*b(y)^*=-b(y)^*b(x)^*`$ und
$`\{b(x),b(y)^*\}=\langle x,y\rangle`$. Sei $`U\subseteq\mathbb R^L`$ der Spann der $`u_{k_a}`$, und seien $`P_A`$, $`P_B=1-P_A`$ die
Koordinatenprojektionen auf $`\{0,\dots,n-1\}`$ und dessen Komplement. Diagonalisieren der Form $`\langle P_Ax,P_Ay\rangle`$
auf $`U`$ liefert eine Orthonormalbasis $`\phi_1,\dots,\phi_N`$ von $`U`$ mit $`\langle P_A\phi_a,P_A\phi_b\rangle=\nu_a\delta_{ab}`$,
wobei $`\nu_a=\|P_A\phi_a\|^2\in[0,1]`$, also $`\langle P_B\phi_a,P_B\phi_b\rangle=(1-\nu_a)\delta_{ab}`$. Da $`(x_1,\dots,x_N)\mapsto b(x_1)^*\cdots b(x_N)^*|0\rangle`$
multilinear und alternierend ist, gilt $`\Phi=\pm\prod_a\big(\sqrt{\nu_a}f_a^*+\sqrt{1-\nu_a}\,g_a^*\big)|0\rangle`$ mit
$`f_a=b(P_A\phi_a/\sqrt{\nu_a})`$, $`g_a=b(P_B\phi_a/\sqrt{1-\nu_a})`$ (Terme mit Koeffizient null weggelassen).
Ausmultiplizieren ergibt $`\Phi=\sum_I\pm\prod_{a\in I}\sqrt{\nu_a}\prod_{a\notin I}\sqrt{1-\nu_a}\,F_I^*G_{I^c}^*|0\rangle`$,
summiert über die Teilmengen $`I`$, mit den geordneten Produkten $`F_I^*`$, $`G_{I^c}^*`$. Für $`j\ge n`$ enthält der String
von $`c_j^*`$ ganz $`A`$: $`c_j^*=(-1)^{N_A}\otimes\tilde c_j^*`$ mit Jordan-Wigner-Operatoren $`\tilde c_j`$ der Plätze $`\ge n`$,
während für $`i<n`$ gilt $`c_i^*=\hat c_i^*\otimes1`$. Da $`(-1)^{N_A}`$ auf $`|0_A\rangle`$ gleich $`1`$ ist, folgt
$`F_I^*G_{I^c}^*|0\rangle=\xi_I\otimes\chi_{I^c}`$ mit $`\xi_I=\hat F_I^*|0_A\rangle`$, $`\chi_{I^c}=\tilde G_{I^c}^*|0_B\rangle`$, und
diese bilden nach der CAR-Rechnung aus Satz 9.10 Orthonormalsysteme. Das ist eine Schmidt-Zerlegung (Lemma 7.7) mit
Gewichten $`p_I=\prod_{a\in I}\nu_a\prod_{a\notin I}(1-\nu_a)`$ (die Vorzeichen gehen in $`\xi_I`$ ein), der Verteilung
unabhängiger Bernoulli-Variablen, also, wegen $`\log p_I=\sum_{a\in I}\log\nu_a+\sum_{a\notin I}\log(1-\nu_a)`$ und
$`\sum_{I\ni a}p_I=\nu_a`$, $`S_A=-\sum_Ip_I\log p_I=\sum_a\eta(\nu_a)`$. Schließlich ist $`C_{ij}=\sum_a\phi_a(i)\phi_a(j)`$ (man entwickle $`c_i^*c_j`$
nach den $`b(\phi_a)`$), also hat $`C_A=\sum_a(P_A\phi_a)(P_A\phi_a)^T`$ die Eigenwerte $`\nu_a`$ auf den orthogonalen Vektoren
$`P_A\phi_a`$ und sonst $`0`$, und $`\eta(0)=0`$. $`\square`$

Satz 9.16 stammt von I. Peschel (*J. Phys. A* 36 (2003) L205).

*Rechnung:* Satz 9.16 für die Grundzustände der offenen XX-Kette mit $`L=12`$, im Vergleich mit `SpinChain.entanglement_entropy` — [Tutorial, Abschnitt 9.6](../../tutorial/de/09_spin_chains.md#s-9-6).

<a id="s-9-7"></a>

## 9.7 Womit man vergleichen sollte

*Rechnung:* Der Vergleich mit anderer Software (QuSpin, Codes für exakte Diagonalisierung, Tensornetzwerke) — [Tutorial, Abschnitt 9.7](../../tutorial/de/09_spin_chains.md#s-9-7).

<a id="x-9"></a>

## Aufgaben

1. (Impulssektoren von Hand) Man bestimme für $`L=4`$ und $`N=2`$ die Bahnen der sechs Konfigurationen unter $`T`$ und
   die Dimensionen $`d_{2,m}`$, $`m=0,\dots,3`$.
2. (Summenregeln) Sei $`\psi`$ der Grundzustand des Heisenberg-Rings mit geradem $`L`$. Man zeige
   $`E_0=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle`$ und $`\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle=0`$.
3. (Zu viel Symmetrie) Man zeige, dass im Sektor $`N=L/2`$ der Vektor mit lauter Koeffizienten $`1`$ ein Vielfaches
   von $`(S^-)^{L/2}|{\uparrow\cdots\uparrow}\rangle`$ und ein Eigenvektor des Heisenberg-Rings zur Energie $`JL/4`$ ist und
   dass jeder Krylov-Raum $`\mathcal K_m(H,v)`$ in jedem $`H`$-invarianten Unterraum liegt, der $`v`$ enthält. Was liefert
   `krylov.ground_state`, wenn man von diesem Vektor aus startet?

<details>
<summary>Lösungen</summary>

1. $`0011\to1001\to1100\to0110\to0011`$ ist eine Bahn mit $`R=4`$, und $`0101\to1010\to0101`$ eine mit $`R=2`$. Nach
   Proposition 9.5 trägt die erste zu jedem $`m`$ bei, die zweite zu den $`m`$ mit $`2m\equiv0\pmod4`$, d. h. $`m=0,2`$:
   $`d_{2,m}=2,1,2,1`$.

2. Der Grundzustand ist ein einfaches Singulett (Satz 9.9), und $`T`$ kommutiert mit $`H`$, also $`T\psi=e^{ik}\psi`$. Da
   $`\mathbf s_i\cdot\mathbf s_{i+1}=T^i(\mathbf s_0\cdot\mathbf s_1)T^{*i}`$, ist $`\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle`$ für alle
   $`L`$ Bindungen gleich, und $`E_0=J\sum_i\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle`$. Da
   $`\mathbf S^2\psi=0`$, ist
   $`0=\langle\mathbf S^2\rangle=\sum_{i,j}\langle\mathbf s_i\cdot\mathbf s_j\rangle=L\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle`$, wiederum
   wegen der Translationsinvarianz.

3. $`S^-`$ senkt einen der Spins nach oben ab; $`L/2`$-malige Anwendung auf den vollständig polarisierten Zustand
   erreicht jede Konfiguration mit $`L/2`$ Spins nach oben in $`(L/2)!`$ Reihenfolgen, alle mit Koeffizient $`1`$.
   $`|{\uparrow\cdots\uparrow}\rangle`$ hat alle $`L`$ Bindungen parallel, Energie $`JL/4`$, und $`S^-`$ kommutiert mit $`H`$,
   also ist der Vektor ein Eigenvektor zur selben Energie. Ist $`W`$ invariant und $`v\in W`$, so ist $`H^jv\in W`$ für
   alle $`j`$. Hier ist $`\mathcal K_m=\mathbb Cv`$ für alle $`m`$, und `ground_state` liefert $`E=JL/4`$, das obere Ende des
   Spektrums (jedes $`\mathbf s_i\cdot\mathbf s_j=\frac12F_{ij}-\frac14\le\frac14`$ nach Proposition 4.15(a), da $`F_{ij}`$ eine
   symmetrische Permutationsmatrix mit Eigenwerten $`\pm1`$ ist). Ein zufälliger Startvektor vermeidet diese Falle: Ein
   gaußsches $`v`$ hat eine Dichte, also hat es fast sicher eine von null verschiedene Komponente $`P_1v`$ im Eigenraum
   von $`\lambda_1`$ (Korollar A.7(c); für ein reelles $`v`$ wende man es in $`W=\mathbb R^n`$ auf den echten Unterraum
   $`\{v:P_1v=0\}`$ an, echt, weil die reelle Matrix $`P_1\ne0`$ eine von null verschiedene Spalte hat), und dann hat nach
   Satz 9.14 ein $`H`$-invarianter Krylov-Raum $`\mathcal K_m(H,v)`$ den kleinsten Ritz-Wert $`\lambda_1`$, während hier der
   invariante Raum $`\mathcal K_1=\mathbb Cv`$ den Ritz-Wert $`JL/4`$ hat.

</details>

*Rechnung:* Aufgaben 1–3 numerisch — [Tutorial, Aufgaben](../../tutorial/de/09_spin_chains.md#x-9).

**Weiter:** [Kapitel 10](10_quantum_computing.md) behandelt Quantencomputing mit Sektoren: dekohärenzfreie
Unterräume und rauschfreie Teilsysteme, Fehlerkorrektur in Operatoralgebren und Mixer, die ein erhaltenes
Hamming-Gewicht respektieren – die Algebra $`\bigoplus_wM_{\binom nw}`$ dieses Kapitels in einer neuen Rolle.
