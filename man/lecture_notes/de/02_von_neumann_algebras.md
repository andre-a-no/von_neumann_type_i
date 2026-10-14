[English](../02_von_neumann_algebras.md) · [Español](../es/02_von_neumann_algebras.md) · [Français](../fr/02_von_neumann_algebras.md) · **Deutsch** · [中文](../zh/02_von_neumann_algebras.md) · [日本語](../ja/02_von_neumann_algebras.md) · [Русский](../ru/02_von_neumann_algebras.md)

# 2. Endlichdimensionale von-Neumann-Algebren

Eine von-Neumann-Algebra ist eine Operatoralgebra, die unter der Bildung der „zweifachen
Kommutante“ abgeschlossen ist. In endlicher Dimension ist dies eine rein algebraische
Bedingung, und der Hauptsatz dieses Kapitels besagt, dass jede solche Algebra bis auf
einen unitären Basiswechsel eine direkte Summe voller Matrixalgebren ist, von denen jede
mit einer gewissen Vielfachheit auftritt. Auf dieser Struktur ist die Bibliothek aufgebaut.

*In der Bibliothek:* Die Algebra $`\bigoplus_cM_{k_c}(\mathbb C)`$ ist die Klasse `TypeIAlgebra` — [Tutorial, Abschnitt 2.5](../../tutorial/de/02_von_neumann_algebras.md#s-2-5).

**Sie lernen** den Bikommutantensatz, den Struktursatz mit vollständigem Beweis, wie eine
Symmetriegruppe eine solche Algebra erzeugt, was das Zentrum ist und was eine bedingte
Erwartung ist. Die numerischen Kontrollen dieses Kapitels stehen im
[Tutorial 2](../../tutorial/de/02_von_neumann_algebras.md).

Im Folgenden bedeutet *$`*`$-Algebra auf $`\mathbb C^n`$* einen komplexen Untervektorraum
$`\mathcal M\subseteq M_n(\mathbb C)`$, der unter Produkten und Adjungierten abgeschlossen ist
und die Identität $`1`$ enthält.

<a id="s-2-1"></a>

## 2.1 Kommutanten und der Bikommutantensatz

Für eine Menge $`S\subseteq M_n(\mathbb C)`$ ist die **Kommutante**
$`S'=\{X\in M_n(\mathbb C):XA=AX\text{ für alle }A\in S\}`$. Sie ist stets eine unitale Algebra;
ist $`S`$ unter Adjungierten abgeschlossen ($`S^*=S`$), so ist $`S'`$ eine $`*`$-Algebra, denn aus
$`XA=AX`$ für alle $`A\in S`$ folgt $`X^*A^*=A^*X^*`$ für alle $`A\in S`$. Jedes Element von $`S`$
kommutiert mit jedem Element von $`S'`$, also gilt $`S\subseteq S''`$; auf $`S'`$ angewandt ergibt
dies $`S'\subseteq S'''`$, und Übergang zu den Kommutanten in $`S\subseteq S''`$ (der Inklusionen
umkehrt) ergibt $`S'''\subseteq S'`$. Folglich ist $`S'=S'''`$.

*Rechnung:* Kommutanten als Kerne, numerisch — [Tutorial, Abschnitt 2.4](../../tutorial/de/02_von_neumann_algebras.md#s-2-4). *In der Bibliothek:* `decompose.commutant_basis`, `decompose.commutant_dimension` — [Tutorial, Abschnitt 4.2](../../tutorial/de/04_symmetries.md#s-4-2).

Wirkt eine Algebra auf einem Unterraum $`E`$ von $`\mathbb C^n`$, so identifizieren wir $`E`$
mittels einer Orthonormalbasis mit $`\mathbb C^{\dim E}`$, sodass alle für $`\mathbb C^n`$
formulierten Ergebnisse darauf anwendbar sind.

**Definition.** Eine **von-Neumann-Algebra** auf $`\mathbb C^n`$ ist eine $`*`$-Algebra
$`\mathcal M`$ mit $`\mathcal M=\mathcal M''`$. Ihr **Zentrum** ist
$`Z(\mathcal M)=\mathcal M\cap\mathcal M'`$. Ein **Faktor** ist eine von-Neumann-Algebra mit
$`Z(\mathcal M)=\mathbb C1`$.

Die Definition ist für Operatoren auf unendlichdimensionalen Hilberträumen gemacht, die dieser Kurs
nicht behandelt. In endlicher Dimension ist die Bedingung $`\mathcal M=\mathcal M''`$ automatisch
erfüllt:

<a id="r-2-1"></a>

**Lemma 2.1.** Sei $`\mathcal M`$ eine $`*`$-Algebra auf $`\mathbb C^n`$ und $`K\subseteq\mathbb C^n`$
ein Unterraum mit $`\mathcal MK\subseteq K`$. Dann liegt die Orthogonalprojektion $`P`$ auf $`K`$
in $`\mathcal M'`$.

*Beweis.* Für $`A\in\mathcal M`$ gilt $`AK\subseteq K`$ und $`A^*K\subseteq K`$ (da $`A^*\in\mathcal M`$).
Aus der zweiten Inklusion folgt $`AK^\perp\subseteq K^\perp`$: Für $`y\in K^\perp`$, $`x\in K`$ ist
$`\langle x,Ay\rangle=\langle A^*x,y\rangle=0`$. Also bildet $`A`$ den Raum $`K`$ in $`K`$ und
$`K^\perp`$ in $`K^\perp`$ ab, d. h. $`AP=PAP=PA`$. $`\square`$

<a id="r-2-2"></a>

**Satz 2.2 (Bikommutantensatz).** Jede $`*`$-Algebra $`\mathcal M`$ auf
$`\mathbb C^n`$ erfüllt $`\mathcal M=\mathcal M''`$.

*Beweis.* Es gilt $`\mathcal M\subseteq\mathcal M''`$; sei $`T\in\mathcal M''`$.

*Schritt 1: ein Vektor.* Für $`x\in\mathbb C^n`$ ist der Unterraum $`K=\mathcal Mx`$ invariant
unter $`\mathcal M`$ und enthält $`x=1x`$. Nach Lemma 2.1 liegt seine Projektion $`P`$ in
$`\mathcal M'`$, also $`TP=PT`$ und $`Tx=TPx=PTx\in K`$: Es gibt ein $`A\in\mathcal M`$ mit $`Tx=Ax`$.

*Schritt 2: alle Vektoren zugleich.* Auf $`\mathbb C^n\otimes\mathbb C^n=(\mathbb C^n)^{\oplus n}`$
betrachte man die $`*`$-Algebra $`\tilde{\mathcal M}=\{A\oplus\dots\oplus A:A\in\mathcal M\}`$. Eine
Blockmatrix $`[X_{ij}]_{i,j=1}^n`$ kommutiert genau dann mit jedem $`A\oplus\dots\oplus A`$, wenn
$`X_{ij}A=AX_{ij}`$ für alle $`i,j`$ gilt; also ist $`\tilde{\mathcal M}'=\{[X_{ij}]:X_{ij}\in\mathcal M'\}`$.
Der Operator $`\tilde T=T\oplus\dots\oplus T`$ kommutiert mit jedem solchen $`[X_{ij}]`$, denn
$`TX_{ij}=X_{ij}T`$. Folglich ist $`\tilde T\in\tilde{\mathcal M}''`$. Man wende Schritt 1 auf
$`\tilde{\mathcal M}`$, $`\tilde T`$ und den Vektor $`\xi=(e_1,\dots,e_n)`$ an: Es gibt ein
$`A\in\mathcal M`$ mit $`(Te_1,\dots,Te_n)=(Ae_1,\dots,Ae_n)`$, d. h. $`T=A\in\mathcal M`$. $`\square`$

*Rechnung:* Satz 2.2 an einem Beispiel — [Tutorial, Abschnitt 2.4](../../tutorial/de/02_von_neumann_algebras.md#s-2-4).

In endlicher Dimension sind „von-Neumann-Algebra“ und „unitale $`*`$-Algebra von Matrizen“
also dasselbe. Der Bikommutantensatz ist dennoch der richtige Gesichtspunkt: Physikalisch
ist eine Observablenalgebra sehr oft als Kommutante *gegeben* – als die Operatoren, die mit
einer Symmetrie vertauschen (Abschnitt 2.4).

<a id="s-2-2"></a>

## 2.2 Das Zentrum

<a id="r-2-3"></a>

**Lemma 2.3.** Sei $`\mathcal Z`$ eine kommutative $`*`$-Algebra auf $`\mathbb C^n`$. Es gibt
von null verschiedene, paarweise orthogonale Projektionen $`P_1,\dots,P_C\in\mathcal Z`$ mit
$`\sum_cP_c=1`$, sodass $`\mathcal Z=\mathrm{span}\{P_1,\dots,P_C\}`$. Sie sind die
minimalen von null verschiedenen Projektionen von $`\mathcal Z`$ und bis auf die Reihenfolge
eindeutig.

*Beweis.* Jedes $`z\in\mathcal Z`$ ist normal ($`z^*\in\mathcal Z`$ kommutiert mit $`z`$). Es mögen
$`z_1,\dots,z_r`$ den Raum $`\mathcal Z`$ aufspannen. Wir zerlegen $`\mathbb C^n`$ in gemeinsame
Eigenräume. Die Eigenräume von $`z_1`$ sind paarweise orthogonal und spannen $`\mathbb C^n`$ auf
(Satz 1.3); sie sind invariant unter $`z_2,\dots,z_r`$ und deren Adjungierten (die mit $`z_1`$
kommutieren), also ist jedes $`z_j`$, eingeschränkt auf einen Eigenraum von $`z_1`$, wieder
normal, und wir können jeden Eigenraum von $`z_1`$ in Eigenräume von $`z_2`$ zerlegen, und so
fort. Das Ergebnis ist eine orthogonale Zerlegung $`\mathbb C^n=E_1\oplus\dots\oplus E_C`$ in
von null verschiedene Unterräume, auf denen jedes $`z\in\mathcal Z`$ als Skalar $`\chi_c(z)`$
wirkt, wobei wir Unterräume mit derselben Funktion $`\chi_c:\mathcal Z\to\mathbb C`$
zusammengefasst haben, sodass $`\chi_c\neq\chi_d`$ für $`c\neq d`$. Sei $`P_c`$ die Projektion auf
$`E_c`$. Dann gilt $`z=\sum_c\chi_c(z)P_c`$ für jedes $`z\in\mathcal Z`$.

Es bleibt $`P_c\in\mathcal Z`$ zu zeigen. Zu $`d\neq c`$ wähle man $`z_{cd}\in\mathcal Z`$ mit
$`\chi_c(z_{cd})\neq\chi_d(z_{cd})`$. Das Produkt
$`\prod_{d\neq c}\big(z_{cd}-\chi_d(z_{cd})1\big)/\big(\chi_c(z_{cd})-\chi_d(z_{cd})\big)`$
liegt in $`\mathcal Z`$, wirkt als $`1`$ auf $`E_c`$ und als $`0`$ auf jedem $`E_d`$, $`d\neq c`$; es ist
also $`P_c`$ (für $`C=1`$ ist $`P_1=1`$: das leere Produkt ist $`1`$). Minimalität und Eindeutigkeit: Für eine Projektion $`\sum_c a_cP_c`$ gilt
$`a_c\in\{0,1\}`$, also sind die von null verschiedenen Projektionen von $`\mathcal Z`$ genau die
Summen über nichtleere Mengen der $`P_c`$. $`\square`$

Angewandt auf $`\mathcal Z=Z(\mathcal M)`$ heißen die Projektionen $`P_c`$ die **minimalen
zentralen Projektionen** von $`\mathcal M`$. Jedes $`A\in\mathcal M`$ zerfällt als
$`A=\sum_cAP_c`$ mit $`AP_c=P_cAP_c`$, also ist
$`\mathcal M=\bigoplus_c\mathcal MP_c`$ (eine orthogonale direkte Summe von Algebren, die auf
$`E_c=P_c\mathbb C^n`$ wirken).

*In der Bibliothek:* Die zentralen Elemente $`\sum_cz_cP_c`$, darunter die minimalen zentralen Projektionen, erzeugt `TypeIAlgebra.central` — [Tutorial, Abschnitt 2.5](../../tutorial/de/02_von_neumann_algebras.md#s-2-5).

<a id="r-2-4"></a>

**Lemma 2.4.** Jedes $`\mathcal M_c=\mathcal MP_c`$ ist, als Algebra auf $`E_c`$ mit Einselement
$`P_c`$, ein Faktor.

*Beweis.* Sei $`y\in\mathcal M_c`$ mit $`\mathcal M_c`$ vertauschbar. Für $`A\in\mathcal M`$ gilt
$`Ay=AP_cy=(AP_c)y=y(AP_c)=yP_cA=yA`$, da $`P_c`$ zentral ist und $`y=yP_c=P_cy`$.
Also ist $`y\in Z(\mathcal M)`$, $`y=\sum_d\chi_d(y)P_d`$, und $`y=yP_c=\chi_c(y)P_c`$. $`\square`$

<a id="s-2-3"></a>

## 2.3 Der Struktursatz

<a id="r-2-5"></a>

**Lemma 2.5 (Faktoren).** Sei $`\mathcal N`$ ein Faktor auf einem Raum $`E`$. Es gibt ganze Zahlen
$`k,m\ge1`$ und eine unitäre Abbildung $`V:\mathbb C^k\otimes\mathbb C^m\to E`$ mit
$`V^*\mathcal NV=M_k(\mathbb C)\otimes1_m`$. Ferner gilt $`V^*\mathcal N'V=1_k\otimes M_m(\mathbb C)`$.

*Beweis.* *Schritt 1: eine minimale Projektion.* $`\mathcal N`$ enthält von null verschiedene
Projektionen (z. B. $`1`$). Sei $`e\in\mathcal N`$ eine von null verschiedene Projektion minimalen
Rangs. Wir behaupten $`e\mathcal Ne=\mathbb Ce`$. Sei $`x=x^*\in e\mathcal Ne`$. Jede
Spektralprojektion $`q`$ von $`x`$ zu einem von null verschiedenen Eigenwert $`\lambda`$ ist
$`q=p(x)`$, wobei $`p`$ das Lagrange-Polynom mit $`p(\lambda)=1`$ und $`p(\mu)=0`$ für
$`\mu\in(\mathrm{spec}\,x\cup\{0\})\setminus\{\lambda\}`$ ist (Satz 1.5(a),(b): $`p`$ stimmt auf
$`\mathrm{spec}\,x`$ mit der Indikatorfunktion von $`\{\lambda\}`$ überein). Wegen $`p(0)=0`$ ist
$`p(t)=t\,r(t)`$ mit einem Polynom $`r`$, also $`q=x\,r(x)=r(x)\,x\in\mathcal N`$, und aus $`x=ex=xe`$ folgt
$`q=eq=qe`$; also liegt das Bild von $`q`$ im Bild von $`e`$. Ist $`q\neq0`$, so liefert die Minimalität
$`\mathrm{rank}\,q\ge\mathrm{rank}\,e`$, also stimmen die beiden Bilder überein und $`q=e`$. Da die
Spektralprojektionen zu verschiedenen Eigenwerten orthogonal sind, hat $`x`$ höchstens einen von null
verschiedenen Eigenwert $`\lambda`$, mit Spektralprojektion $`e`$, also $`x=\lambda e`$ (und $`x=0=0\cdot e`$, falls
es keinen gibt). Ein allgemeines Element von $`e\mathcal Ne`$ ist $`x_1+ix_2`$ mit
$`x_{1,2}=x_{1,2}^*\in e\mathcal Ne`$.

*Schritt 2: $`e`$ erzeugt $`E`$.* Sei $`F=\mathrm{span}\,\mathcal Ne E`$. Dieser Raum ist
invariant unter $`\mathcal N`$, also liegt seine Projektion $`p_F`$ in $`\mathcal N'`$ (Lemma 2.1).
Er ist auch invariant unter der $`*`$-Algebra $`\mathcal N'`$: Für $`b\in\mathcal N'`$ ist
$`bAe\xi=Aeb\xi\in F`$; nach Lemma 2.1, angewandt auf $`\mathcal N'`$, gilt also
$`p_F\in\mathcal N''=\mathcal N`$ (Satz 2.2). Folglich ist $`p_F\in Z(\mathcal N)=\mathbb C1`$,
und $`p_F\neq0`$ (da $`eE\subseteq F`$), also $`F=E`$.

*Schritt 3: Matrixeinheiten.* Ein $`v\in\mathcal N`$ mit $`v^*v=e`$ heiße *$`e`$-Isometrie*; dann ist
$`vv^*`$ eine Projektion in $`\mathcal N`$. (Zunächst ist $`v=ve`$: Für jedes $`x`$ gilt
$`\|v(1-e)x\|^2=\langle(1-e)x,v^*v(1-e)x\rangle=\langle(1-e)x,e(1-e)x\rangle=0`$. Daher ist
$`vv^*v=ve=v`$ und $`(vv^*)^2=vev^*=vv^*`$. Für die $`v_i`$ unten liefert dies $`e_iv_i=v_i`$ und
$`v_i^*e_i=v_i^*`$.) Man wähle
$`e`$-Isometrien $`v_1=e,v_2,\dots,v_k`$ mit paarweise orthogonalen Bildern $`e_i=v_iv_i^*`$, $`k`$
maximal (die Bilder sind orthogonale Unterräume derselben Dimension $`\mathrm{rank}\,e`$,
also $`k\le\dim E/\mathrm{rank}\,e`$). Wir behaupten $`\sum_ie_i=1`$. Andernfalls ist
$`g=1-\sum_ie_i\neq0`$ eine Projektion in $`\mathcal N`$. Nach Schritt 2 ist $`gAe\neq0`$ für ein
$`A\in\mathcal N`$ (wäre $`gAe=0`$ für alle $`A`$, so wäre $`gAe\xi=0`$ für alle $`A,\xi`$, also
verschwände $`g`$ auf $`F=E`$). Sei $`y=gAe`$. Dann ist $`y^*y\in e\mathcal Ne`$ positiv und von null
verschieden, also $`y^*y=\lambda e`$ mit $`\lambda>0`$, und $`v=y/\sqrt\lambda`$ ist eine
$`e`$-Isometrie mit Bild in $`gE`$, orthogonal zu allen $`e_i`$ – im Widerspruch zur Maximalität.

Für $`i\neq j`$ ist $`v_i^*v_j=v_i^*e_ie_jv_j=0`$. Man setze $`E_{ij}=v_iv_j^*`$. Dann gilt
$`E_{ij}E_{kl}=v_iv_j^*v_kv_l^*=\delta_{jk}v_iev_l^*=\delta_{jk}E_{il}`$,
$`E_{ij}^*=E_{ji}`$ und $`\sum_iE_{ii}=1`$. Für $`A\in\mathcal N`$ ist
$`A=\sum_{ij}e_iAe_j=\sum_{ij}v_i(v_i^*Av_j)v_j^*`$ und $`v_i^*Av_j\in e\mathcal Ne=\mathbb Ce`$,
also $`A=\sum_{ij}a_{ij}E_{ij}`$ mit $`a_{ij}\in\mathbb C`$. Somit ist $`\mathcal N=\mathrm{span}\{E_{ij}\}`$.

*Schritt 4: die unitäre Abbildung.* Sei $`m=\mathrm{rank}\,e`$ und $`f_1,\dots,f_m`$ eine
Orthonormalbasis von $`eE`$. Die Vektoren $`v_if_r`$ sind orthonormal:
$`\langle v_if_r,v_jf_s\rangle=\langle f_r,v_i^*v_jf_s\rangle=\delta_{ij}\langle f_r,ef_s\rangle=\delta_{ij}\delta_{rs}`$.
Sie spannen $`E`$ auf, da $`E=\bigoplus_ie_iE`$ und $`e_iE=v_ieE`$. Man definiere
$`V(|i\rangle\otimes f_r)=v_if_r`$, wobei $`|1\rangle,\dots,|k\rangle`$ die Standardbasis von
$`\mathbb C^k`$ ist und wir $`\mathbb C^m`$ mittels der $`f_r`$ mit $`eE`$ identifizieren.
Dann gilt $`E_{ij}v_lf_r=\delta_{jl}v_if_r`$, d. h. $`V^*E_{ij}V=|i\rangle\langle j|\otimes1_m`$, und
$`V^*\mathcal NV=M_k(\mathbb C)\otimes1_m`$.

*Schritt 5: die Kommutante.* Man schreibe $`X\in M_k\otimes M_m`$ als Blockmatrix
$`X=\sum_{ij}|i\rangle\langle j|\otimes X_{ij}`$. Für feste $`a,b`$ ist
$`(|a\rangle\langle b|\otimes1)X=\sum_j|a\rangle\langle j|\otimes X_{bj}`$ und
$`X(|a\rangle\langle b|\otimes1)=\sum_i|i\rangle\langle b|\otimes X_{ia}`$. Sind beide gleich,
so liefert der Vergleich der Blöcke an der Position $`(a,j)`$ mit $`j\neq b`$ die Gleichung
$`X_{bj}=0`$ und an der Position $`(a,b)`$ die Gleichung $`X_{bb}=X_{aa}`$. Da $`a,b`$ beliebig sind,
ist $`X=1_k\otimes Y`$ mit $`Y=X_{11}`$; umgekehrt kommutiert jedes $`1_k\otimes Y`$ mit
$`M_k\otimes1_m`$. $`\square`$

<a id="r-2-6"></a>

**Satz 2.6 (Struktur endlichdimensionaler von-Neumann-Algebren).** Sei
$`\mathcal M`$ eine von-Neumann-Algebra auf $`\mathbb C^n`$. Es gibt $`C\ge1`$, ganze Zahlen
$`k_c,m_c\ge1`$ mit $`\sum_ck_cm_c=n`$ und eine unitäre Abbildung
$`V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n`$, sodass
```math
V^*\mathcal MV=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},\qquad
V^*\mathcal M'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),
```
und $`V^*Z(\mathcal M)V=\mathrm{span}\{V^*P_1V,\dots,V^*P_CV\}`$, wobei $`P_1,\dots,P_C`$ die
minimalen zentralen Projektionen von $`\mathcal M`$ sind (Lemma 2.3) und $`V^*P_cV`$ die Projektion
auf den $`c`$-ten Summanden ist. Die Paare $`(k_c,m_c)`$ sind bis auf die Reihenfolge eindeutig: Für
*jede* unitäre Abbildung $`V`$ mit der ersten Eigenschaft bildet $`V(\cdot)V^*`$ die Projektionen auf die
Summanden auf die $`P_c`$ ab (siehe den Beweis), also sind $`k_c^2=\dim\mathcal MP_c`$ und
$`m_c=\mathrm{rank}\,P_c/k_c`$ durch $`\mathcal M`$ bestimmt. Insbesondere ist $`\mathcal M\cong\bigoplus_cM_{k_c}(\mathbb C)`$ als $`*`$-Algebren, und
$`\dim\mathcal M=\sum_ck_c^2`$, $`\dim\mathcal M'=\sum_cm_c^2`$.

*Beweis.* Nach Lemma 2.3 und der anschließenden Bemerkung ist $`\mathcal M=\bigoplus_c\mathcal M_c`$,
wobei $`\mathcal M_c`$ auf $`E_c`$ wirkt, und jedes $`\mathcal M_c`$ ist ein Faktor (Lemma 2.4).
Lemma 2.5 liefert unitäre Abbildungen $`V_c`$ mit $`V_c^*\mathcal M_cV_c=M_{k_c}\otimes1_{m_c}`$;
man setze $`V=\bigoplus_cV_c`$. Zur Kommutante: Jedes $`X\in\mathcal M'`$ kommutiert mit jedem
$`P_c`$, also ist $`X=\sum_cP_cXP_c`$, und $`P_cXP_c`$ kommutiert mit $`\mathcal M_c`$; umgekehrt
kommutiert jede solche Summe mit $`\mathcal M`$. Lemma 2.5 bestimmt $`\mathcal M_c'`$ auf $`E_c`$, und Konjugation
mit $`V`$ liefert die Formel für $`V^*\mathcal M'V`$.

*Eindeutigkeit.* Sei $`V`$ eine beliebige unitäre Abbildung mit $`V^*\mathcal MV=\mathcal B:=\bigoplus_{c=1}^{C'}M_{k_c}\otimes1_{m_c}`$,
und sei $`Q_c`$ die Projektion auf den $`c`$-ten Summanden. Ein Element $`\bigoplus_cX_c`$ von $`\mathcal B`$ ($`X_c\in M_{k_c}\otimes1_{m_c}`$)
ist genau dann zentral, wenn jedes $`X_c`$ mit $`M_{k_c}\otimes1_{m_c}`$ kommutiert, d. h. (Lemma 2.5, Schritt 5) $`X_c\in(M_{k_c}\otimes1_{m_c})\cap(1_{k_c}\otimes M_{m_c})=\mathbb C1`$
(ein Element $`a\otimes1=1\otimes b`$ hat lauter Diagonalblöcke gleich $`b`$ und lauter Nebendiagonalblöcke $`0`$, also ist $`a`$ ein Skalar).
Daher ist $`Z(\mathcal B)=\mathrm{span}\{Q_1,\dots,Q_{C'}\}`$, und $`Z(\mathcal M)=VZ(\mathcal B)V^*`$ wird von den von null verschiedenen,
paarweise orthogonalen Projektionen $`VQ_cV^*`$ mit Summe $`1`$ aufgespannt. Jede von ihnen ist eine von null verschiedene Projektion von $`Z(\mathcal M)`$, also (Lemma 2.3)
die Summe über eine nichtleere Menge $`S_c`$ der $`P_d`$; die $`S_c`$ sind disjunkt (Orthogonalität) und überdecken alle $`d`$ (Summe $`1`$), und wegen
$`C'=\dim Z(\mathcal M)=C`$ ist jedes $`S_c`$ einelementig. Nach Umnummerierung ist also $`VQ_cV^*=P_c`$, woraus
$`\dim\mathcal MP_c=\dim(M_{k_c}\otimes1_{m_c})=k_c^2`$ und $`\mathrm{rank}\,P_c=k_cm_c`$ folgen. $`\square`$

*In der Bibliothek:* Die Algebra $`\bigoplus_cM_{k_c}(\mathbb C)`$ ist eine `TypeIAlgebra` — [Tutorial, Abschnitt 2.5](../../tutorial/de/02_von_neumann_algebras.md#s-2-5); das unitäre $`V`$ und die Paare $`(k_c,m_c)`$ berechnet `decompose.wedderburn` (mit den Schnittstellen `from_group`, `from_generators`) — [Tutorial, Abschnitt 4.3](../../tutorial/de/04_symmetries.md#s-4-3).

Die Zahlen $`k_c`$ beschreiben die Algebra selbst. Die **Vielfachheiten** $`m_c`$ beschreiben,
wie sie in $`M_n(\mathbb C)`$ liegt: Sie sind für die algebraische Struktur unsichtbar,
bestimmen aber die Einschränkung der Spur von $`\mathbb C^n`$ auf $`\mathcal M`$,
```math
\mathrm{Tr}_{\mathbb C^n}(A)=\sum_cm_c\mathrm{Tr}(A_c),\qquad A=V\Big(\bigoplus_cA_c\otimes1_{m_c}\Big)V^*.
```
Deshalb arbeitet die Bibliothek mit einer Familie von Spuren (Kapitel 3). Ebenso hat für selbstadjungiertes
$`A`$ die Einschränkung von $`A`$ auf $`V(\mathbb C^{k_c}\otimes\mathbb C^{m_c})`$ die Eigenwerte des Blocks $`A_c`$, jeden mit
dem $`m_c`$-Fachen seiner Vielfachheit in $`A_c`$, da $`A_c\otimes1_{m_c}`$ aus $`m_c`$ Kopien von $`A_c`$ besteht.

*In der Bibliothek:* Die Spuren `Tr_blunt`, `Tr_norm`, `tau_vN` — [Tutorial, Abschnitt 3.1](../../tutorial/de/03_traces_and_states.md#s-3-1).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Die $`P_c`$ sind **Superauswahlsektoren** oder
> **Symmetriesektoren**: Keine Observable in $`\mathcal M`$ verbindet verschiedene Sektoren,
> daher ist eine Messung von „welcher Sektor“ mit jeder anderen Messung verträglich. Das
> Zentrum ist der *klassische* Teil des Systems. Die Vielfachheit $`m_c`$ ist eine Entartung jeder Observablen in $`\mathcal M`$.

<a id="s-2-4"></a>

## 2.4 Beispiel: eine Symmetriegruppe und ihre Kommutante

Eine endliche Gruppe $`G`$ wirke auf $`\mathbb C^n`$ durch eine unitäre Darstellung $`\pi`$. Die
Operatoren, die mit der Symmetrie vertauschen, bilden $`\pi(G)'`$, eine von-Neumann-Algebra;
deren Kommutante $`\pi(G)''`$ ist der Spann der $`\pi(g)`$ (Satz 2.2, angewandt auf diesen
Spann, der eine $`*`$-Algebra ist: Er enthält $`\pi(e)=1`$, die Produkte $`\pi(g)\pi(h)=\pi(gh)`$
und die Adjungierten $`\pi(g)^*=\pi(g^{-1})`$).

<a id="r-2-7"></a>

**Proposition 2.7.** Man wende Satz 2.6 auf $`\mathcal M=\pi(G)''`$ an. Dann gilt
$`V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}`$, wobei jedes $`\rho_c`$ eine irreduzible
unitäre Darstellung von $`G`$ auf $`\mathbb C^{k_c}`$ ist und $`\rho_c,\rho_d`$ für $`c\neq d`$
inäquivalent sind. Folglich gilt
```math
\pi\cong\bigoplus_c\rho_c\otimes1_{m_c},\qquad
\pi(G)''\cong\bigoplus_cM_{\dim\rho_c},\qquad
\pi(G)'\cong\bigoplus_cM_{m_c}:
```
Die Blöcke der symmetrischen Operatoren $`\pi(G)'`$ haben als Größen die Vielfachheiten der
irreduziblen Darstellungen, und die Rollen von $`k_c`$ und $`m_c`$ werden zwischen einer
Algebra und ihrer Kommutante vertauscht.

*Beweis.* Es ist $`\pi(g)\in\mathcal M`$, also $`V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}`$
mit gewissen $`\rho_c(g)\in M_{k_c}`$, und $`g\mapsto\rho_c(g)`$ ist eine unitäre Darstellung,
weil $`\pi`$ eine ist. Da die $`\pi(g)`$ den Raum $`\mathcal M`$ aufspannen, spannen die
$`\rho_c(g)`$ den Raum $`M_{k_c}`$ auf, und stärker noch ist die Abbildung
$`\sum_ga_g\pi(g)\mapsto\big(\sum_ga_g\rho_c(g)\big)_c`$ surjektiv auf $`\bigoplus_cM_{k_c}`$.
*Irreduzibel:* Ein unter allen $`\rho_c(g)`$ invarianter Unterraum von $`\mathbb C^{k_c}`$ ist
invariant unter deren Spann $`M_{k_c}`$, also $`0`$ oder $`\mathbb C^{k_c}`$. *Inäquivalent:* Ist
$`\rho_d(g)=T\rho_c(g)T^{-1}`$ für alle $`g`$ und ein invertierbares $`T`$, $`c\neq d`$, so erfüllt
jedes Element $`(x_{c'})_{c'}`$ des Bildes $`x_d=Tx_cT^{-1}`$, was das Element mit
$`x_c=1`$, $`x_d=0`$ ausschließt – ein Widerspruch zur Surjektivität. Die letzten beiden
Isomorphien sind Satz 2.6. $`\square`$

*In der Bibliothek:* Die Zerlegung einer Darstellung einer endlichen Gruppe, `decompose.from_group` — [Tutorial, Abschnitt 4.3](../../tutorial/de/04_symmetries.md#s-4-3).

Die Zerlegung einer Darstellung in irreduzible ist bis auf Äquivalenz eindeutig, sodass die
$`m_c`$ *die* Vielfachheiten in jeder Zerlegung sind. *Beweis.* Seien $`\rho,\rho'`$ irreduzibel und
$`T`$ ein Verflechtungsoperator, $`T\rho(g)=\rho'(g)T`$ für alle $`g`$. Kern und Bild von $`T`$ sind
invariant, also ist $`T=0`$ oder $`T`$ invertierbar (Lemma von Schur); ist $`\rho'=\rho`$, so gilt
dasselbe für $`T-\lambda1`$ mit einem Eigenwert $`\lambda`$ von $`T`$, also $`T=\lambda1`$. Daher bilden
die Verflechtungsoperatoren von $`\rho`$ nach $`\rho'`$ einen Raum der Dimension $`1`$, falls
$`\rho\cong\rho'`$ (für einen invertierbaren Verflechtungsoperator $`S`$ identifiziert $`T\mapsto S^{-1}T`$ sie mit
den Verflechtungsoperatoren von $`\rho`$ nach sich selbst, und diese sind $`\mathbb C1`$), und $`0`$ sonst (ein von null
verschiedener wäre invertierbar, also eine Äquivalenz). Ist $`\pi\cong\bigoplus_j\sigma_j`$ mit irreduziblen $`\sigma_j`$
(Wiederholungen erlaubt), so ist ein Verflechtungsoperator von einer irreduziblen $`\rho`$ nach
$`\bigoplus_j\sigma_j`$ ein Tupel $`(T_j)_j`$ von Verflechtungsoperatoren von $`\rho`$ nach $`\sigma_j`$, und die Komposition
mit der Äquivalenz ändert keine Dimensionen; also bilden die Verflechtungsoperatoren von $`\rho`$ nach
$`\pi`$ einen Raum der Dimension $`\#\{j:\sigma_j\cong\rho\}`$, eine Zahl, die nur von $`\pi`$
und $`\rho`$ abhängt. Für die Zerlegung aus Proposition 2.7, deren $`\rho_c`$ paarweise inäquivalent
sind, ist sie $`m_c`$ für $`\rho\cong\rho_c`$ und $`0`$, wenn $`\rho`$ zu keinem $`\rho_c`$ äquivalent ist.
$`\square`$

Wir nehmen $`G=S_3`$, das die drei Faktoren von $`(\mathbb C^2)^{\otimes3}`$ permutiert, also
drei Qubits oder drei Spins $`\tfrac12`$: $`\pi(\sigma)|x\rangle=|\sigma\cdot x\rangle`$ für $`x\in\{0,1\}^3`$, wobei
$`(\sigma\cdot x)_{\sigma(i)}=x_i`$. Sei $`F=\{v:\pi(g)v=v\text{ für alle }g\}`$ der **symmetrische Unterraum**. Die
Zerlegung lautet: $`F`$ (Dimension $`4`$) trägt die triviale Darstellung mit Vielfachheit $`4`$; sein
orthogonales Komplement (Dimension $`4`$) trägt die zweidimensionale **Standarddarstellung** (die
Permutationsdarstellung von $`S_3`$ auf den Vektoren von $`\mathbb C^3`$ mit Koordinatensumme $`0`$) mit
Vielfachheit $`2`$; die Signumsdarstellung tritt nicht auf. Also gilt $`\pi(S_3)''\cong M_1\oplus M_2`$
(Dimension $`5`$) und $`\pi(S_3)'\cong M_4\oplus M_2`$ (Dimension $`20`$).

*Beweis.* Man schreibe $`V^*\pi(S_3)''V=\bigoplus_cM_{k_c}\otimes1_{m_c}`$ (Satz 2.6), sodass
$`\sum_ck_cm_c=8`$ und $`\sum_cm_c^2=\dim\pi(S_3)'`$.

*Schritt 1: $`\dim\pi(S_3)'=20`$.* Wegen $`\pi(g)^*|g\cdot x\rangle=|x\rangle`$ hat die Matrix
$`\pi(g)X\pi(g)^*`$ den Eintrag $`X_{xy}`$ an der Position $`(g\cdot x,g\cdot y)`$. Daher ist $`X\in\pi(S_3)'`$, d. h.
$`\pi(g)X\pi(g)^*=X`$ für alle $`g`$, genau dann, wenn $`(x,y)\mapsto X_{xy}`$ konstant auf den Bahnen von $`S_3`$ ist, das
auf Paaren $`(x,y)`$ durch $`g\cdot(x,y)=(g\cdot x,g\cdot y)`$ wirkt; die Indikatormatrizen der Bahnen bilden eine Basis von
$`\pi(S_3)'`$, dessen Dimension also die Anzahl der Bahnen ist. Ein Paar $`(x,y)`$ ist ein Tripel
$`\big((x_1,y_1),(x_2,y_2),(x_3,y_3)\big)`$ von Elementen der vierelementigen Menge $`\{0,1\}^2`$, und $`S_3`$
permutiert die drei Komponenten; eine Bahn ist also eine Multimenge der Größe $`3`$ aus vier Elementen: $`4`$ der Form
$`\{a,a,a\}`$, $`4\cdot3=12`$ der Form $`\{a,a,b\}`$ mit $`b\neq a`$ und $`4`$ der Form $`\{a,b,c\}`$ mit
paarweise verschiedenen Elementen, insgesamt $`20`$.

*Schritt 2: $`k=(1,2)`$, $`m=(4,2)`$.* Wegen $`k_c\ge1`$ ist $`\sum_cm_c\le\sum_ck_cm_c=8`$. Sei $`M=\max_cm_c`$.
Dann gilt $`20=\sum_cm_c^2\le M\sum_cm_c\le8M`$ und $`M^2\le20`$, also $`M\in\{3,4\}`$. Ist $`M=4`$, so haben die übrigen
$`m_c`$ die Quadratsumme $`4`$, sind also $`(2)`$ oder $`(1,1,1,1)`$. Ist $`M=3`$, so haben die übrigen $`m_c`$ die
Quadratsumme $`11`$ und eine Summe von höchstens $`5`$; wären alle $`\le2`$, so wäre ihre Quadratsumme höchstens
$`2\cdot5=10`$, also ist eines von ihnen gleich $`3`$, und die restlichen haben die Quadratsumme $`2`$: Sie sind $`(1,1)`$.
Folglich ist $`m`$ bis auf die Reihenfolge $`(4,2)`$, $`(4,1,1,1,1)`$ oder $`(3,3,1,1)`$. In den letzten beiden Fällen
erzwingt $`\sum_cm_c=8=\sum_ck_cm_c`$, dass alle $`k_c=1`$ sind, und dann wäre $`\pi(S_3)''\cong\bigoplus_cM_1`$
kommutativ, während die Transpositionen $`(12)`$ und $`(23)`$ nicht vertauschen:
$`\pi((12))\pi((23))|001\rangle=|100\rangle`$, aber $`\pi((23))\pi((12))|001\rangle=|010\rangle`$. Also ist
$`m=(4,2)`$ (wobei die Summanden so nummeriert sind, dass $`m_1=4`$), und $`4k_1+2k_2=8`$ mit $`k_1,k_2\ge1`$ ergibt $`k=(1,2)`$. Insbesondere ist
$`\dim\pi(S_3)''=\sum_ck_c^2=5`$.

*Schritt 3: die triviale und die Signumsdarstellung.* $`\rho_2`$ hält keinen von null verschiedenen Vektor $`u`$ fest: Die Gerade
$`\mathbb Cu`$ wäre ein invarianter Unterraum von $`\mathbb C^2`$, im Widerspruch zur Irreduzibilität
(Proposition 2.7). Man schreibe $`V^*v=v_1\oplus v_2`$ mit $`v_1\in\mathbb C^1\otimes\mathbb C^4`$ und
$`v_2=\sum_{r=1}^2u_r\otimes f_r\in\mathbb C^2\otimes\mathbb C^2`$ ($`f_1,f_2`$ die Standardbasis). Dann ist
$`v\in F`$ genau dann, wenn $`\rho_1(g)v_1=v_1`$ und $`\rho_2(g)u_r=u_r`$ für alle $`g`$ und $`r`$; also hat jedes $`v\in F`$
die Komponente $`v_2=0`$. Der von null verschiedene Vektor $`|000\rangle`$ liegt in $`F`$, also hält der Skalar $`\rho_1(g)`$ einen von null verschiedenen Vektor fest:
$`\rho_1`$ ist trivial, $`F=V(\mathbb C^1\otimes\mathbb C^4)`$ hat die Dimension $`4`$, und sein orthogonales
Komplement ist $`V(\mathbb C^2\otimes\mathbb C^2)`$. Die Signumsdarstellung ist weder zu
$`\rho_1`$ äquivalent (da $`\mathrm{sgn}(12)=-1`$) noch zu $`\rho_2`$ (die Dimensionen sind verschieden); nach der oben gezeigten
Eindeutigkeit der Vielfachheiten tritt sie also $`0`$-mal auf.

*Schritt 4: $`\rho_2`$ ist die Standarddarstellung.* Die Vektoren $`e_1=|100\rangle`$, $`e_2=|010\rangle`$,
$`e_3=|001\rangle`$ werden von den $`\pi(\sigma)`$ so permutiert, wie $`\sigma`$ die Positionen permutiert; also ist
$`W=\{\sum_ix_ie_i:\sum_ix_i=0\}`$ invariant, und $`\pi`$ wirkt darauf durch die Standarddarstellung.
$`W\perp F`$: Für $`v\in F`$ und $`w=\sum_ix_ie_i\in W`$ ist $`\langle v,w\rangle=\langle\pi(g)v,w\rangle=\langle v,\pi(g)^{-1}w\rangle`$
für jedes $`g`$, also $`\langle v,w\rangle=\langle v,\frac16\sum_g\pi(g)w\rangle`$, und
$`\frac16\sum_g\pi(g)w=\frac13\big(\sum_ix_i\big)(e_1+e_2+e_3)=0`$ (jedes $`e_i`$ wird von
genau zwei der sechs Permutationen auf jedes $`e_j`$ abgebildet). Daher ist $`V^*W\subseteq\mathbb C^2\otimes\mathbb C^2`$. Die Abbildungen
$`T_r:\mathbb C^2\otimes\mathbb C^2\to\mathbb C^2`$, $`u\otimes f\mapsto\langle f_r,f\rangle u`$, erfüllen
$`T_r(\rho_2(g)\otimes1)=\rho_2(g)T_r`$. Man wähle $`r`$ mit $`T_rV^*|_W\neq0`$ (verschwänden beide, so wäre $`V^*W`$
gleich $`0`$). Dann ist $`T_rV^*|_W`$ ein von null verschiedener Verflechtungsoperator von $`\pi|_W`$ nach $`\rho_2`$; sein Bild ist ein von null
verschiedener invarianter Unterraum von $`\mathbb C^2`$, also $`\mathbb C^2`$, und wegen $`\dim W=2`$ ist er bijektiv. Also ist $`\rho_2`$
äquivalent zur Standarddarstellung (die nebenbei irreduzibel ist). $`\square`$

*Rechnung:* Die Dimensionen $`20`$ und $`5`$ sowie Satz 2.2 an diesem Beispiel — [Tutorial, Abschnitt 2.4](../../tutorial/de/02_von_neumann_algebras.md#s-2-4). *In der Bibliothek:* `decompose.permutation_representation`, `decompose.commutant_dimension` — [Tutorial, Abschnitt 4.2](../../tutorial/de/04_symmetries.md#s-4-2).

Die Vielfachheiten sind im Spektrum sichtbar. Jedes selbstadjungierte $`K\in\pi(S_3)'`$ hat die Form
$`V^*KV=(1_1\otimes\kappa_1)\oplus(1_2\otimes\kappa_2)`$ mit selbstadjungierten $`\kappa_1\in M_4`$, $`\kappa_2\in M_2`$
(Satz 2.6). Auf $`F`$ sind seine Eigenwerte die von $`\kappa_1`$; auf dem Komplement von $`F`$ sind es die
von $`\kappa_2`$, jeder doppelt, da $`1_2\otimes\kappa_2`$ aus zwei Kopien von $`\kappa_2`$ besteht. Also hat *jeder* Eigenwert von
$`K`$ auf dem Komplement von $`F`$ gerade Vielfachheit. Für ein zufälliges $`K`$ gilt mehr. Die Mittelung über die
Gruppe, $`E(X)=\frac1{|G|}\sum_g\pi(g)X\pi(g)^*`$, bildet $`M_8(\mathbb C)`$ auf die Kommutante ab und ist dort
die Identität (Aufgabe 3); für eine Ginibre-Matrix $`Z`$ (unabhängige komplexe Standard-Gauß-Einträge)
hat $`K=E(Z+Z^*)`$ daher eine Dichte auf dem reellen Vektorraum $`\pi(S_3)'_h`$ der selbstadjungierten Elemente
(Satz A.15 und die anschließende Bemerkung). Nach Proposition 4.7(a) (ein Vorgriff; ihr Beweis verwendet
Satz 2.6 und Anhang A, nicht diesen Absatz) sind mit Wahrscheinlichkeit eins die sechs
Eigenwerte von $`\kappa_1`$ und $`\kappa_2`$ paarweise verschieden: $`K`$ hat dann $`4`$ einfache Eigenwerte auf $`F`$ und
$`2`$ zweifach entartete auf dessen Komplement. Schließlich ist $`P_F=\frac1{|G|}\sum_g\pi(g)`$ die Orthogonalprojektion
auf $`F`$: $`\pi(h)P_F=P_F`$ für alle $`h`$, also liegt ihr Bild in $`F`$; $`P_Fv=v`$ für $`v\in F`$; und
$`P_F^*=P_F`$ durch die Substitution $`g\to g^{-1}`$.

*Rechnung:* Die zweifach entarteten Eigenwerte eines zufälligen $`K`$ — [Tutorial, Abschnitt 2.4](../../tutorial/de/02_von_neumann_algebras.md#s-2-4).

<a id="s-2-5"></a>

## 2.5 Die Algebra in der Bibliothek

*In der Bibliothek:* Wie `TypeIAlgebra` die Algebra $`\bigoplus_cM_{k_c}(\mathbb C)`$, ihre Operatoren und ihr Zentrum speichert — [Tutorial, Abschnitt 2.5](../../tutorial/de/02_von_neumann_algebras.md#s-2-5).

<a id="s-2-6"></a>

## 2.6 Bedingte Erwartungen

Eine **bedingte Erwartung** von $`\mathcal M`$ auf eine von-Neumann-Unteralgebra
$`\mathcal N\subseteq\mathcal M`$ (eine in $`\mathcal M`$ enthaltene $`*`$-Algebra mit demselben Einselement) ist eine
lineare Abbildung $`E:\mathcal M\to\mathcal N`$, die positiv ist ($`A\ge0\Rightarrow E(A)\ge0`$), $`E(B)=B`$ für $`B\in\mathcal N`$ erfüllt und eine
$`\mathcal N`$-Bimodulabbildung ist: $`E(B_1AB_2)=B_1E(A)B_2`$ für $`B_1,B_2\in\mathcal N`$.

<a id="r-2-8"></a>

**Proposition 2.8.** Für $`\mathcal M=\bigoplus_cM_{k_c}`$ auf $`\bigoplus_c\mathbb C^{k_c}`$, mit $`P_c`$ der Projektion auf
den $`c`$-ten Block, und $`\mathcal N=Z(\mathcal M)=\mathrm{span}\{P_1,\dots,P_C\}`$ (Aufgabe 1) ist
```math
E(A)=\sum_c\frac{\mathrm{Tr}\,A_c}{k_c}\,P_c
```
eine bedingte Erwartung. Sie erhält jede Spur der Form
$`\sum_cw_c\mathrm{Tr}\,A_c`$.

*Beweis.* Die Linearität ist klar; aus $`A\ge0`$ folgt $`\mathrm{Tr}\,A_c=\sum_i\langle e_i,Ae_i\rangle\ge0`$
(Summe über die Standardbasisvektoren $`e_i`$ des $`c`$-ten Blocks), also $`E(A)\ge0`$. Für $`B=\sum_cb_cP_c`$ ist $`E(B)=\sum_cb_ck_cP_c/k_c=B`$. Für $`B_i=\sum_cb^{(i)}_cP_c`$ ist
$`\mathrm{Tr}(B_1AB_2)_c=b^{(1)}_cb^{(2)}_c\mathrm{Tr}\,A_c`$, woraus die
Bimoduleigenschaft folgt. Schließlich ist $`\mathrm{Tr}(E(A)_c)=\mathrm{Tr}\,A_c`$ für jedes $`c`$. $`\square`$

*In der Bibliothek:* $`E`$ ist `channels.center_expectation` — [Tutorial, Abschnitt 2.6](../../tutorial/de/02_von_neumann_algebras.md#s-2-6).

Die Gruppenmittelung aus Abschnitt 2.4 ist eine bedingte Erwartung auf die Kommutante
$`\pi(G)'`$ innerhalb von $`M_n(\mathbb C)`$, und sie erhält $`\mathrm{Tr}`$ (Aufgabe 3). Ebenso ist auf der
vollen Matrixalgebra $`M_D(\mathbb C)`$, $`D=\sum_ck_c`$, die $`\mathcal M=\bigoplus_cM_{k_c}`$ als ihre
Blockdiagonalmatrizen enthält, die Abbildung $`E_{\mathcal M}(X)=\sum_cP_cXP_c`$ eine bedingte Erwartung auf
$`\mathcal M`$, die $`\mathrm{Tr}`$ erhält. *Beweis.* $`P_cXP_c`$ verschwindet außerhalb des $`c`$-ten Diagonalblocks,
also ist $`E_{\mathcal M}(X)\in\mathcal M`$. Für $`X\ge0`$ ist $`\langle v,P_cXP_cv\rangle=\langle P_cv,XP_cv\rangle\ge0`$. Ein Element
$`B\in\mathcal M`$ ist blockdiagonal, also $`BP_c=P_cB`$ und $`B=\sum_cP_cBP_c`$; daher ist $`E_{\mathcal M}(B)=B`$ und
$`E_{\mathcal M}(B_1XB_2)=\sum_cB_1P_cXP_cB_2=B_1E_{\mathcal M}(X)B_2`$ für $`B_1,B_2\in\mathcal M`$. Schließlich ist
$`\sum_c\mathrm{Tr}(P_cXP_c)=\sum_c\mathrm{Tr}(P_cX)=\mathrm{Tr}\,X`$. $`\square`$

*In der Bibliothek:* Die Gruppenmittelung ist `decompose.twirl` — [Tutorial, Abschnitt 4.3](../../tutorial/de/04_symmetries.md#s-4-3).

Für einen **Zustand** $`\omega`$ auf $`\mathcal M`$ (ein lineares Funktional mit $`\omega(A^*A)\ge0`$ für alle $`A`$ und $`\omega(1)=1`$;
Zustände werden in Kapitel 3 untersucht) liefern Proposition 2.8 und die Linearität
$`\omega(E(A))=\sum_c\omega(P_c)\mathrm{Tr}(A_c)/k_c`$, und $`\omega\circ E=\omega`$ auf dem Zentrum, da dort
$`E(z)=z`$ gilt. (Bedingte Erwartungen sind im Allgemeinen nicht eindeutig: Auf $`\mathbb C1`$ liefert jeder
Zustand $`\varphi`$ eine, nämlich $`A\mapsto\varphi(A)1`$; auf $`M_2(\mathbb C)`$ liefern die Zustände $`A\mapsto A_{11}`$ und
$`A\mapsto A_{22}`$ zwei verschiedene.) All dies sind Beispiele derselben Idee:
die Projektion auf diejenigen Operatoren, die eine gewisse Struktur „nicht sehen“.

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Die Abbildung $`E_{\mathcal M}`$ löscht die
> Kohärenzen zwischen Superauswahlsektoren, und genau das geschieht, wenn man die Ladung (den
> Sektor) misst, ohne das Ergebnis zu notieren. Die bedingte Erwartung $`E`$ auf das Zentrum geht
> weiter: Sie behält nur die klassische Information, das Gewicht $`\omega(P_c)`$ jedes Sektors, und ersetzt
> den Zustand innerhalb jedes Sektors durch den maximal gemischten (die normierte Spur $`\mathrm{Tr}(A_c)/k_c`$).

<a id="s-2-7"></a>

## 2.7 Womit man vergleichen sollte

*Rechnung:* Der Vergleich mit anderer Software (GAP, SageMath, QuSpin, QuTiP) — [Tutorial, Abschnitt 2.7](../../tutorial/de/02_von_neumann_algebras.md#s-2-7).

<a id="x-2"></a>

## Aufgaben

1. Man zeige, dass das Zentrum von $`M_k(\mathbb C)`$ gleich $`\mathbb C1`$ ist (man verwende die
   Matrixeinheiten $`|i\rangle\langle j|`$), und dass $`Z(\bigoplus_cM_{k_c})`$ die Dimension $`C`$ hat.
2. Sei $`S`$ die zyklische Verschiebung auf $`\mathbb C^3`$, $`Se_i=e_{i+1\bmod3}`$. Man bestimme
   $`\{S\}'`$ und $`\{S\}''`$.
3. Für eine unitäre Darstellung $`\pi`$ einer endlichen Gruppe zeige man, dass
   $`E(X)=\frac1{|G|}\sum_g\pi(g)X\pi(g)^*`$ eine bedingte Erwartung von
   $`M_n(\mathbb C)`$ auf $`\pi(G)'`$ ist, die $`\mathrm{Tr}`$ erhält.

<details>
<summary>Lösungen</summary>

1. Kommutiert $`X`$ mit allen $`|i\rangle\langle j|`$, so gilt
   $`X|i\rangle\langle j|=|i\rangle\langle j|X`$; wendet man beide Seiten auf $`e_j`$ an und
   vergleicht, so erhält man $`Xe_i=X_{jj}e_i`$ für alle $`i,j`$, also $`X=X_{11}1`$. In einer
   direkten Summe ist ein Element $`\bigoplus X_c`$ genau dann zentral, wenn jedes $`X_c`$ es
   ist; das Zentrum ist also $`\{\bigoplus_cz_c1_{k_c}\}`$ und hat die Dimension $`C`$.

2. $`S`$ ist unitär mit den Eigenwerten $`1,\omega,\omega^2`$, $`\omega=e^{2\pi i/3}`$, jeder einfach: Die
   Vektoren $`u_j=\sum_i\omega^{-ij}e_i`$ ($`j=0,1,2`$) erfüllen $`Su_j=\omega^ju_j`$ und sind orthogonal.
   Ein mit $`S`$ vertauschbarer Operator erhält dessen Eigenräume, ist also in der
   Eigenbasis diagonal, und jede Diagonalmatrix kommutiert mit $`S`$: $`\{S\}'\cong\mathbb C^3`$. Sie enthält
   $`1,S,S^2`$, die in der Eigenbasis $`\mathrm{diag}(\lambda^0)`$, $`\mathrm{diag}(\lambda)`$,
   $`\mathrm{diag}(\lambda^2)`$ für $`\lambda=(1,\omega,\omega^2)`$ sind, linear unabhängig (Vandermonde: Die Eigenwerte
   sind paarweise verschieden); also ist $`\{S\}'=\mathrm{span}\{1,S,S^2\}`$, die zirkulanten Matrizen. Dann ist auch
   $`\{S\}''=\{S\}'`$, denn in der Eigenbasis von $`S`$ ist eine Matrix, die mit allen
   Diagonalmatrizen kommutiert, diagonal (man vergleiche die Einträge von $`XD`$ und $`DX`$
   für ein $`D`$ mit paarweise verschiedenen Diagonaleinträgen).

3. $`E`$ ist linear und positiv (eine Summe von Abbildungen $`X\mapsto UXU^*`$). Für $`h\in G`$ ist
   $`\pi(h)E(X)\pi(h)^*=\frac1{|G|}\sum_g\pi(hg)X\pi(hg)^*=E(X)`$, also $`E(X)\in\pi(G)'`$. Ist
   $`X\in\pi(G)'`$, so ist jeder Summand gleich $`X`$, also $`E(X)=X`$. Für $`B_1,B_2\in\pi(G)'`$ ist
   $`\pi(g)B_1XB_2\pi(g)^*=B_1\pi(g)X\pi(g)^*B_2`$, woraus die Bimoduleigenschaft folgt.
   Schließlich ist $`\mathrm{Tr}(\pi(g)X\pi(g)^*)=\mathrm{Tr}\,X`$ für jedes $`g`$, also
   $`\mathrm{Tr}\,E(X)=\mathrm{Tr}\,X`$.

*Rechnung:* Aufgabe 2 numerisch — [Tutorial, Aufgaben](../../tutorial/de/02_von_neumann_algebras.md#x-2).

</details>

**Weiter:** [Kapitel 3](03_traces_and_states.md) klassifiziert die Spuren auf
$`\bigoplus_cM_{k_c}`$, führt Zustände und Dichten ein und beweist die später verwendeten
Entropie- und Spurungleichungen.
