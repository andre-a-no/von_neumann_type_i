[English](../02_von_neumann_algebras.md) · [Español](../es/02_von_neumann_algebras.md) · [Français](../fr/02_von_neumann_algebras.md) · **Deutsch** · [中文](../zh/02_von_neumann_algebras.md) · [日本語](../ja/02_von_neumann_algebras.md) · [Русский](../ru/02_von_neumann_algebras.md)

# 2. Endlichdimensionale von-Neumann-Algebren

Eine von-Neumann-Algebra ist eine Operatoralgebra, die unter der Bildung der „zweifachen
Kommutante“ abgeschlossen ist. In endlicher Dimension ist dies eine rein algebraische
Bedingung, und der Hauptsatz dieses Kapitels besagt, dass jede solche Algebra bis auf
einen unitären Basiswechsel eine direkte Summe voller Matrixalgebren ist, von denen jede
mit einer gewissen Vielfachheit auftritt. Auf dieser Struktur ist die Bibliothek aufgebaut.

**Sie lernen** den Bikommutantensatz, den Struktursatz mit vollständigem Beweis, wie eine
Symmetriegruppe eine solche Algebra erzeugt, was das Zentrum ist und wie die Bibliothek
all dies darstellt.

Im Folgenden bedeutet *$*$-Algebra auf $\mathbb C^n$* einen komplexen Untervektorraum
$\mathcal M\subseteq M_n(\mathbb C)$, der unter Produkten und Adjungierten abgeschlossen ist
und die Identität $1$ enthält.

## 2.1 Kommutanten und der Bikommutantensatz

Für eine Menge $S\subseteq M_n(\mathbb C)$ ist die **Kommutante**
$S'=\{X\in M_n(\mathbb C):XA=AX\text{ für alle }A\in S\}$. Sie ist stets eine unitale Algebra;
ist $S$ unter Adjungierten abgeschlossen ($S^*=S$), so ist $S'$ eine $*$-Algebra, denn aus
$XA=AX$ für alle $A\in S$ folgt $X^*A^*=A^*X^*$ für alle $A\in S$. Jedes Element von $S$
kommutiert mit jedem Element von $S'$, also gilt $S\subseteq S''$; auf $S'$ angewandt ergibt
dies $S'\subseteq S'''$, und Übergang zu den Kommutanten in $S\subseteq S''$ (der Inklusionen
umkehrt) ergibt $S'''\subseteq S'$. Folglich ist $S'=S'''$.

Wirkt eine Algebra auf einem Unterraum $E$ von $\mathbb C^n$, so identifizieren wir $E$
mittels einer Orthonormalbasis mit $\mathbb C^{\dim E}$, sodass alle für $\mathbb C^n$
formulierten Ergebnisse darauf anwendbar sind.

**Definition.** Eine **von-Neumann-Algebra** auf $\mathbb C^n$ ist eine $*$-Algebra
$\mathcal M$ mit $\mathcal M=\mathcal M''$. Ihr **Zentrum** ist
$Z(\mathcal M)=\mathcal M\cap\mathcal M'$. Ein **Faktor** ist eine von-Neumann-Algebra mit
$Z(\mathcal M)=\mathbb C1$.

Die Definition ist für Operatoren auf unendlichdimensionalen Hilberträumen gemacht, die dieser Kurs
nicht behandelt. In endlicher Dimension ist die Bedingung $\mathcal M=\mathcal M''$ automatisch
erfüllt:

**Lemma 2.1.** Sei $\mathcal M$ eine $*$-Algebra auf $\mathbb C^n$ und $K\subseteq\mathbb C^n$
ein Unterraum mit $\mathcal MK\subseteq K$. Dann liegt die Orthogonalprojektion $P$ auf $K$
in $\mathcal M'$.

*Beweis.* Für $A\in\mathcal M$ gilt $AK\subseteq K$ und $A^*K\subseteq K$ (da $A^*\in\mathcal M$).
Aus der zweiten Inklusion folgt $AK^\perp\subseteq K^\perp$: Für $y\in K^\perp$, $x\in K$ ist
$\langle x,Ay\rangle=\langle A^*x,y\rangle=0$. Also bildet $A$ den Raum $K$ in $K$ und
$K^\perp$ in $K^\perp$ ab, d. h. $AP=PAP=PA$. $\square$

**Satz 2.2 (Bikommutantensatz).** Jede $*$-Algebra $\mathcal M$ auf
$\mathbb C^n$ erfüllt $\mathcal M=\mathcal M''$.

*Beweis.* Es gilt $\mathcal M\subseteq\mathcal M''$; sei $T\in\mathcal M''$.

*Schritt 1: ein Vektor.* Für $x\in\mathbb C^n$ ist der Unterraum $K=\mathcal Mx$ invariant
unter $\mathcal M$ und enthält $x=1x$. Nach Lemma 2.1 liegt seine Projektion $P$ in
$\mathcal M'$, also $TP=PT$ und $Tx=TPx=PTx\in K$: Es gibt ein $A\in\mathcal M$ mit $Tx=Ax$.

*Schritt 2: alle Vektoren zugleich.* Auf $\mathbb C^n\otimes\mathbb C^n=(\mathbb C^n)^{\oplus n}$
betrachte man die $*$-Algebra $\tilde{\mathcal M}=\{A\oplus\dots\oplus A:A\in\mathcal M\}$. Eine
Blockmatrix $[X_{ij}]_{i,j=1}^n$ kommutiert genau dann mit jedem $A\oplus\dots\oplus A$, wenn
$X_{ij}A=AX_{ij}$ für alle $i,j$ gilt; also ist $\tilde{\mathcal M}'=\{[X_{ij}]:X_{ij}\in\mathcal M'\}$.
Der Operator $\tilde T=T\oplus\dots\oplus T$ kommutiert mit jedem solchen $[X_{ij}]$, denn
$TX_{ij}=X_{ij}T$. Folglich ist $\tilde T\in\tilde{\mathcal M}''$. Man wende Schritt 1 auf
$\tilde{\mathcal M}$, $\tilde T$ und den Vektor $\xi=(e_1,\dots,e_n)$ an: Es gibt ein
$A\in\mathcal M$ mit $(Te_1,\dots,Te_n)=(Ae_1,\dots,Ae_n)$, d. h. $T=A\in\mathcal M$. $\square$

In endlicher Dimension sind „von-Neumann-Algebra“ und „unitale $*$-Algebra von Matrizen“
also dasselbe. Der Bikommutantensatz ist dennoch der richtige Gesichtspunkt: Physikalisch
ist eine Observablenalgebra sehr oft als Kommutante *gegeben* – als die Operatoren, die mit
einer Symmetrie vertauschen (Abschnitt 2.4).

## 2.2 Das Zentrum

**Lemma 2.3.** Sei $\mathcal Z$ eine kommutative $*$-Algebra auf $\mathbb C^n$. Es gibt
von null verschiedene, paarweise orthogonale Projektionen $P_1,\dots,P_C\in\mathcal Z$ mit
$\sum_cP_c=1$, sodass $\mathcal Z=\operatorname{span}\{P_1,\dots,P_C\}$. Sie sind die
minimalen von null verschiedenen Projektionen von $\mathcal Z$ und bis auf die Reihenfolge
eindeutig.

*Beweis.* Jedes $z\in\mathcal Z$ ist normal ($z^*\in\mathcal Z$ kommutiert mit $z$). Es mögen
$z_1,\dots,z_r$ den Raum $\mathcal Z$ aufspannen. Wir zerlegen $\mathbb C^n$ in gemeinsame
Eigenräume. Die Eigenräume von $z_1$ sind paarweise orthogonal und spannen $\mathbb C^n$ auf
(Satz 1.3); sie sind invariant unter $z_2,\dots,z_r$ und deren Adjungierten (die mit $z_1$
kommutieren), also ist jedes $z_j$, eingeschränkt auf einen Eigenraum von $z_1$, wieder
normal, und wir können jeden Eigenraum von $z_1$ in Eigenräume von $z_2$ zerlegen, und so
fort. Das Ergebnis ist eine orthogonale Zerlegung $\mathbb C^n=E_1\oplus\dots\oplus E_C$ in
von null verschiedene Unterräume, auf denen jedes $z\in\mathcal Z$ als Skalar $\chi_c(z)$
wirkt, wobei wir Unterräume mit derselben Funktion $\chi_c:\mathcal Z\to\mathbb C$
zusammengefasst haben, sodass $\chi_c\neq\chi_d$ für $c\neq d$. Sei $P_c$ die Projektion auf
$E_c$. Dann gilt $z=\sum_c\chi_c(z)P_c$ für jedes $z\in\mathcal Z$.

Es bleibt $P_c\in\mathcal Z$ zu zeigen. Zu $d\neq c$ wähle man $z_{cd}\in\mathcal Z$ mit
$\chi_c(z_{cd})\neq\chi_d(z_{cd})$. Das Produkt
$\prod_{d\neq c}\big(z_{cd}-\chi_d(z_{cd})1\big)/\big(\chi_c(z_{cd})-\chi_d(z_{cd})\big)$
liegt in $\mathcal Z$, wirkt als $1$ auf $E_c$ und als $0$ auf jedem $E_d$, $d\neq c$; es ist
also $P_c$ (für $C=1$ ist $P_1=1$: das leere Produkt ist $1$). Minimalität und Eindeutigkeit: Für eine Projektion $\sum_c a_cP_c$ gilt
$a_c\in\{0,1\}$, also sind die von null verschiedenen Projektionen von $\mathcal Z$ genau die
Summen über nichtleere Mengen der $P_c$. $\square$

Angewandt auf $\mathcal Z=Z(\mathcal M)$ heißen die Projektionen $P_c$ die **minimalen
zentralen Projektionen** von $\mathcal M$. Jedes $A\in\mathcal M$ zerfällt als
$A=\sum_cAP_c$ mit $AP_c=P_cAP_c$, also ist
$\mathcal M=\bigoplus_c\mathcal MP_c$ (eine orthogonale direkte Summe von Algebren, die auf
$E_c=P_c\mathbb C^n$ wirken).

**Lemma 2.4.** Jedes $\mathcal M_c=\mathcal MP_c$ ist, als Algebra auf $E_c$ mit Einselement
$P_c$, ein Faktor.

*Beweis.* Sei $y\in\mathcal M_c$ mit $\mathcal M_c$ vertauschbar. Für $A\in\mathcal M$ gilt
$Ay=AP_cy=(AP_c)y=y(AP_c)=yP_cA=yA$, da $P_c$ zentral ist und $y=yP_c=P_cy$.
Also ist $y\in Z(\mathcal M)$, $y=\sum_d\chi_d(y)P_d$, und $y=yP_c=\chi_c(y)P_c$. $\square$

## 2.3 Der Struktursatz

**Lemma 2.5 (Faktoren).** Sei $\mathcal N$ ein Faktor auf einem Raum $E$. Es gibt ganze Zahlen
$k,m\ge1$ und eine unitäre Abbildung $V:\mathbb C^k\otimes\mathbb C^m\to E$ mit
$V^*\mathcal NV=M_k(\mathbb C)\otimes1_m$. Ferner gilt $V^*\mathcal N'V=1_k\otimes M_m(\mathbb C)$.

*Beweis.* *Schritt 1: eine minimale Projektion.* $\mathcal N$ enthält von null verschiedene
Projektionen (z. B. $1$). Sei $e\in\mathcal N$ eine von null verschiedene Projektion minimalen
Rangs. Wir behaupten $e\mathcal Ne=\mathbb Ce$. Sei $x=x^*\in e\mathcal Ne$. Jede
Spektralprojektion $q$ von $x$ zu einem von null verschiedenen Eigenwert $\lambda$ ist
$q=p(x)$, wobei $p$ das Lagrange-Polynom mit $p(\lambda)=1$ und $p(\mu)=0$ für
$\mu\in(\operatorname{spec}x\cup\{0\})\setminus\{\lambda\}$ ist (Satz 1.5(a),(b): $p$ stimmt auf
$\operatorname{spec}x$ mit der Indikatorfunktion von $\{\lambda\}$ überein). Wegen $p(0)=0$ ist
$p(t)=t\,r(t)$ mit einem Polynom $r$, also $q=x\,r(x)=r(x)\,x\in\mathcal N$, und aus $x=ex=xe$ folgt
$q=eq=qe$; also liegt das Bild von $q$ im Bild von $e$. Ist $q\neq0$, so liefert die Minimalität
$\operatorname{rank}q\ge\operatorname{rank}e$, also stimmen die beiden Bilder überein und $q=e$. Da die
Spektralprojektionen zu verschiedenen Eigenwerten orthogonal sind, hat $x$ höchstens einen von null
verschiedenen Eigenwert $\lambda$, mit Spektralprojektion $e$, also $x=\lambda e$ (und $x=0=0\cdot e$, falls
es keinen gibt). Ein allgemeines Element von $e\mathcal Ne$ ist $x_1+ix_2$ mit
$x_{1,2}=x_{1,2}^*\in e\mathcal Ne$.

*Schritt 2: $e$ erzeugt $E$.* Sei $F=\operatorname{span}\mathcal Ne E$. Dieser Raum ist
invariant unter $\mathcal N$, also liegt seine Projektion $p_F$ in $\mathcal N'$ (Lemma 2.1).
Er ist auch invariant unter der $*$-Algebra $\mathcal N'$: Für $b\in\mathcal N'$ ist
$bAe\xi=Aeb\xi\in F$; nach Lemma 2.1, angewandt auf $\mathcal N'$, gilt also
$p_F\in\mathcal N''=\mathcal N$ (Satz 2.2). Folglich ist $p_F\in Z(\mathcal N)=\mathbb C1$,
und $p_F\neq0$ (da $eE\subseteq F$), also $F=E$.

*Schritt 3: Matrixeinheiten.* Ein $v\in\mathcal N$ mit $v^*v=e$ heiße *$e$-Isometrie*; dann ist
$vv^*$ eine Projektion in $\mathcal N$. (Zunächst ist $v=ve$: Für jedes $x$ gilt
$\|v(1-e)x\|^2=\langle(1-e)x,v^*v(1-e)x\rangle=\langle(1-e)x,e(1-e)x\rangle=0$. Daher ist
$vv^*v=ve=v$ und $(vv^*)^2=vev^*=vv^*$. Für die $v_i$ unten liefert dies $e_iv_i=v_i$ und
$v_i^*e_i=v_i^*$.) Man wähle
$e$-Isometrien $v_1=e,v_2,\dots,v_k$ mit paarweise orthogonalen Bildern $e_i=v_iv_i^*$, $k$
maximal (die Bilder sind orthogonale Unterräume derselben Dimension $\operatorname{rank}e$,
also $k\le\dim E/\operatorname{rank}e$). Wir behaupten $\sum_ie_i=1$. Andernfalls ist
$g=1-\sum_ie_i\neq0$ eine Projektion in $\mathcal N$. Nach Schritt 2 ist $gAe\neq0$ für ein
$A\in\mathcal N$ (wäre $gAe=0$ für alle $A$, so wäre $gAe\xi=0$ für alle $A,\xi$, also
verschwände $g$ auf $F=E$). Sei $y=gAe$. Dann ist $y^*y\in e\mathcal Ne$ positiv und von null
verschieden, also $y^*y=\lambda e$ mit $\lambda>0$, und $v=y/\sqrt\lambda$ ist eine
$e$-Isometrie mit Bild in $gE$, orthogonal zu allen $e_i$ – im Widerspruch zur Maximalität.

Für $i\neq j$ ist $v_i^*v_j=v_i^*e_ie_jv_j=0$. Man setze $E_{ij}=v_iv_j^*$. Dann gilt
$E_{ij}E_{kl}=v_iv_j^*v_kv_l^*=\delta_{jk}v_iev_l^*=\delta_{jk}E_{il}$,
$E_{ij}^*=E_{ji}$ und $\sum_iE_{ii}=1$. Für $A\in\mathcal N$ ist
$A=\sum_{ij}e_iAe_j=\sum_{ij}v_i(v_i^*Av_j)v_j^*$ und $v_i^*Av_j\in e\mathcal Ne=\mathbb Ce$,
also $A=\sum_{ij}a_{ij}E_{ij}$ mit $a_{ij}\in\mathbb C$. Somit ist $\mathcal N=\operatorname{span}\{E_{ij}\}$.

*Schritt 4: die unitäre Abbildung.* Sei $m=\operatorname{rank}e$ und $f_1,\dots,f_m$ eine
Orthonormalbasis von $eE$. Die Vektoren $v_if_r$ sind orthonormal:
$\langle v_if_r,v_jf_s\rangle=\langle f_r,v_i^*v_jf_s\rangle=\delta_{ij}\langle f_r,ef_s\rangle=\delta_{ij}\delta_{rs}$.
Sie spannen $E$ auf, da $E=\bigoplus_ie_iE$ und $e_iE=v_ieE$. Man definiere
$V(|i\rangle\otimes f_r)=v_if_r$, wobei $|1\rangle,\dots,|k\rangle$ die Standardbasis von
$\mathbb C^k$ ist und wir $\mathbb C^m$ mittels der $f_r$ mit $eE$ identifizieren.
Dann gilt $E_{ij}v_lf_r=\delta_{jl}v_if_r$, d. h. $V^*E_{ij}V=|i\rangle\langle j|\otimes1_m$, und
$V^*\mathcal NV=M_k(\mathbb C)\otimes1_m$.

*Schritt 5: die Kommutante.* Man schreibe $X\in M_k\otimes M_m$ als Blockmatrix
$X=\sum_{ij}|i\rangle\langle j|\otimes X_{ij}$. Für feste $a,b$ ist
$(|a\rangle\langle b|\otimes1)X=\sum_j|a\rangle\langle j|\otimes X_{bj}$ und
$X(|a\rangle\langle b|\otimes1)=\sum_i|i\rangle\langle b|\otimes X_{ia}$. Sind beide gleich,
so liefert der Vergleich der Blöcke an der Position $(a,j)$ mit $j\neq b$ die Gleichung
$X_{bj}=0$ und an der Position $(a,b)$ die Gleichung $X_{bb}=X_{aa}$. Da $a,b$ beliebig sind,
ist $X=1_k\otimes Y$ mit $Y=X_{11}$; umgekehrt kommutiert jedes $1_k\otimes Y$ mit
$M_k\otimes1_m$. $\square$

**Satz 2.6 (Struktur endlichdimensionaler von-Neumann-Algebren).** Sei
$\mathcal M$ eine von-Neumann-Algebra auf $\mathbb C^n$. Es gibt $C\ge1$, ganze Zahlen
$k_c,m_c\ge1$ mit $\sum_ck_cm_c=n$ und eine unitäre Abbildung
$V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n$, sodass
$$
V^*\mathcal MV=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},\qquad
V^*\mathcal M'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),
$$
und $V^*Z(\mathcal M)V=\operatorname{span}\{V^*P_1V,\dots,V^*P_CV\}$, wobei $P_1,\dots,P_C$ die
minimalen zentralen Projektionen von $\mathcal M$ sind (Lemma 2.3) und $V^*P_cV$ die Projektion
auf den $c$-ten Summanden ist. Die Paare $(k_c,m_c)$ sind bis auf die Reihenfolge eindeutig: Für
*jede* unitäre Abbildung $V$ mit der ersten Eigenschaft bildet $V(\cdot)V^*$ die Projektionen auf die
Summanden auf die $P_c$ ab (siehe den Beweis), also sind $k_c^2=\dim\mathcal MP_c$ und
$m_c=\operatorname{rank}P_c/k_c$ durch $\mathcal M$ bestimmt. Insbesondere ist $\mathcal M\cong\bigoplus_cM_{k_c}(\mathbb C)$ als $*$-Algebren, und
$\dim\mathcal M=\sum_ck_c^2$, $\dim\mathcal M'=\sum_cm_c^2$.

*Beweis.* Nach Lemma 2.3 und der anschließenden Bemerkung ist $\mathcal M=\bigoplus_c\mathcal M_c$,
wobei $\mathcal M_c$ auf $E_c$ wirkt, und jedes $\mathcal M_c$ ist ein Faktor (Lemma 2.4).
Lemma 2.5 liefert unitäre Abbildungen $V_c$ mit $V_c^*\mathcal M_cV_c=M_{k_c}\otimes1_{m_c}$;
man setze $V=\bigoplus_cV_c$. Zur Kommutante: Jedes $X\in\mathcal M'$ kommutiert mit jedem
$P_c$, also ist $X=\sum_cP_cXP_c$, und $P_cXP_c$ kommutiert mit $\mathcal M_c$; umgekehrt
kommutiert jede solche Summe mit $\mathcal M$. Lemma 2.5 bestimmt $\mathcal M_c'$ auf $E_c$, und Konjugation
mit $V$ liefert die Formel für $V^*\mathcal M'V$.

*Eindeutigkeit.* Sei $V$ eine beliebige unitäre Abbildung mit $V^*\mathcal MV=\mathcal B:=\bigoplus_{c=1}^{C'}M_{k_c}\otimes1_{m_c}$,
und sei $Q_c$ die Projektion auf den $c$-ten Summanden. Ein Element $\bigoplus_cX_c$ von $\mathcal B$ ($X_c\in M_{k_c}\otimes1_{m_c}$)
ist genau dann zentral, wenn jedes $X_c$ mit $M_{k_c}\otimes1_{m_c}$ kommutiert, d. h. (Lemma 2.5, Schritt 5) $X_c\in(M_{k_c}\otimes1_{m_c})\cap(1_{k_c}\otimes M_{m_c})=\mathbb C1$
(ein Element $a\otimes1=1\otimes b$ hat lauter Diagonalblöcke gleich $b$ und lauter Nebendiagonalblöcke $0$, also ist $a$ ein Skalar).
Daher ist $Z(\mathcal B)=\operatorname{span}\{Q_1,\dots,Q_{C'}\}$, und $Z(\mathcal M)=VZ(\mathcal B)V^*$ wird von den von null verschiedenen,
paarweise orthogonalen Projektionen $VQ_cV^*$ mit Summe $1$ aufgespannt. Jede von ihnen ist eine von null verschiedene Projektion von $Z(\mathcal M)$, also (Lemma 2.3)
die Summe über eine nichtleere Menge $S_c$ der $P_d$; die $S_c$ sind disjunkt (Orthogonalität) und überdecken alle $d$ (Summe $1$), und wegen
$C'=\dim Z(\mathcal M)=C$ ist jedes $S_c$ einelementig. Nach Umnummerierung ist also $VQ_cV^*=P_c$, woraus
$\dim\mathcal MP_c=\dim(M_{k_c}\otimes1_{m_c})=k_c^2$ und $\operatorname{rank}P_c=k_cm_c$ folgen. $\square$

Die Zahlen $k_c$ beschreiben die Algebra selbst. Die **Vielfachheiten** $m_c$ beschreiben,
wie sie in $M_n(\mathbb C)$ liegt: Sie sind für die algebraische Struktur unsichtbar,
bestimmen aber die Einschränkung der Spur von $\mathbb C^n$ auf $\mathcal M$,
$$
\operatorname{Tr}_{\mathbb C^n}(A)=\sum_cm_c\operatorname{Tr}(A_c),\qquad A=V\Big(\bigoplus_cA_c\otimes1_{m_c}\Big)V^*.
$$
Deshalb arbeitet die Bibliothek mit einer Familie von Spuren (Kapitel 3). Ebenso hat für selbstadjungiertes
$A$ die Einschränkung von $A$ auf $V(\mathbb C^{k_c}\otimes\mathbb C^{m_c})$ die Eigenwerte des Blocks $A_c$, jeden mit
dem $m_c$-Fachen seiner Vielfachheit in $A_c$, da $A_c\otimes1_{m_c}$ aus $m_c$ Kopien von $A_c$ besteht.

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Die $P_c$ sind **Superauswahlsektoren** oder
> **Symmetriesektoren**: Keine Observable in $\mathcal M$ verbindet verschiedene Sektoren,
> daher ist eine Messung von „welcher Sektor“ mit jeder anderen Messung verträglich. Das
> Zentrum ist der *klassische* Teil des Systems. Die Vielfachheit $m_c$ ist eine Entartung jeder Observablen in $\mathcal M$.

## 2.4 Beispiel: eine Symmetriegruppe und ihre Kommutante

Eine endliche Gruppe $G$ wirke auf $\mathbb C^n$ durch eine unitäre Darstellung $\pi$. Die
Operatoren, die mit der Symmetrie vertauschen, bilden $\pi(G)'$, eine von-Neumann-Algebra;
deren Kommutante $\pi(G)''$ ist der Spann der $\pi(g)$ (Satz 2.2, angewandt auf diesen
Spann, der eine $*$-Algebra ist: Er enthält $\pi(e)=1$, die Produkte $\pi(g)\pi(h)=\pi(gh)$
und die Adjungierten $\pi(g)^*=\pi(g^{-1})$).

**Proposition 2.7.** Man wende Satz 2.6 auf $\mathcal M=\pi(G)''$ an. Dann gilt
$V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$, wobei jedes $\rho_c$ eine irreduzible
unitäre Darstellung von $G$ auf $\mathbb C^{k_c}$ ist und $\rho_c,\rho_d$ für $c\neq d$
inäquivalent sind. Folglich gilt
$$
\pi\cong\bigoplus_c\rho_c\otimes1_{m_c},\qquad
\pi(G)''\cong\bigoplus_cM_{\dim\rho_c},\qquad
\pi(G)'\cong\bigoplus_cM_{m_c}:
$$
Die Blöcke der symmetrischen Operatoren $\pi(G)'$ haben als Größen die Vielfachheiten der
irreduziblen Darstellungen, und die Rollen von $k_c$ und $m_c$ werden zwischen einer
Algebra und ihrer Kommutante vertauscht.

*Beweis.* Es ist $\pi(g)\in\mathcal M$, also $V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$
mit gewissen $\rho_c(g)\in M_{k_c}$, und $g\mapsto\rho_c(g)$ ist eine unitäre Darstellung,
weil $\pi$ eine ist. Da die $\pi(g)$ den Raum $\mathcal M$ aufspannen, spannen die
$\rho_c(g)$ den Raum $M_{k_c}$ auf, und stärker noch ist die Abbildung
$\sum_ga_g\pi(g)\mapsto\big(\sum_ga_g\rho_c(g)\big)_c$ surjektiv auf $\bigoplus_cM_{k_c}$.
*Irreduzibel:* Ein unter allen $\rho_c(g)$ invarianter Unterraum von $\mathbb C^{k_c}$ ist
invariant unter deren Spann $M_{k_c}$, also $0$ oder $\mathbb C^{k_c}$. *Inäquivalent:* Ist
$\rho_d(g)=T\rho_c(g)T^{-1}$ für alle $g$ und ein invertierbares $T$, $c\neq d$, so erfüllt
jedes Element $(x_{c'})_{c'}$ des Bildes $x_d=Tx_cT^{-1}$, was das Element mit
$x_c=1$, $x_d=0$ ausschließt – ein Widerspruch zur Surjektivität. Die letzten beiden
Isomorphien sind Satz 2.6. $\square$

Die Zerlegung einer Darstellung in irreduzible ist bis auf Äquivalenz eindeutig, sodass die
$m_c$ *die* Vielfachheiten in jeder Zerlegung sind. *Beweis.* Seien $\rho,\rho'$ irreduzibel und
$T$ ein Verflechtungsoperator, $T\rho(g)=\rho'(g)T$ für alle $g$. Kern und Bild von $T$ sind
invariant, also ist $T=0$ oder $T$ invertierbar (Lemma von Schur); ist $\rho'=\rho$, so gilt
dasselbe für $T-\lambda1$ mit einem Eigenwert $\lambda$ von $T$, also $T=\lambda1$. Daher bilden
die Verflechtungsoperatoren von $\rho$ nach $\rho'$ einen Raum der Dimension $1$, falls
$\rho\cong\rho'$ (für einen invertierbaren Verflechtungsoperator $S$ identifiziert $T\mapsto S^{-1}T$ sie mit
den Verflechtungsoperatoren von $\rho$ nach sich selbst, und diese sind $\mathbb C1$), und $0$ sonst (ein von null
verschiedener wäre invertierbar, also eine Äquivalenz). Ist $\pi\cong\bigoplus_j\sigma_j$ mit irreduziblen $\sigma_j$
(Wiederholungen erlaubt), so ist ein Verflechtungsoperator von einer irreduziblen $\rho$ nach
$\bigoplus_j\sigma_j$ ein Tupel $(T_j)_j$ von Verflechtungsoperatoren von $\rho$ nach $\sigma_j$, und die Komposition
mit der Äquivalenz ändert keine Dimensionen; also bilden die Verflechtungsoperatoren von $\rho$ nach
$\pi$ einen Raum der Dimension $\#\{j:\sigma_j\cong\rho\}$, eine Zahl, die nur von $\pi$
und $\rho$ abhängt. Für die Zerlegung aus Proposition 2.7, deren $\rho_c$ paarweise inäquivalent
sind, ist sie $m_c$ für $\rho\cong\rho_c$ und $0$, wenn $\rho$ zu keinem $\rho_c$ äquivalent ist.
$\square$

Wir nehmen $G=S_3$, das die drei Faktoren von $(\mathbb C^2)^{\otimes3}$ permutiert, also
drei Qubits oder drei Spins $\tfrac12$. Die Zerlegung lautet: Der symmetrische Unterraum
(Gesamtspin $\tfrac32$, Dimension $4$) trägt die triviale Darstellung mit Vielfachheit $4$; sein
Komplement (Dimension $4$) trägt die zweidimensionale Standarddarstellung mit Vielfachheit
$2$; die Signumsdarstellung tritt nicht auf. Also gilt $\pi(S_3)''\cong M_1\oplus M_2$
(Dimension $5$) und $\pi(S_3)'\cong M_4\oplus M_2$ (Dimension $20$).

*Beweis.* Man schreibe $\pi(S_3)''\cong\bigoplus_cM_{k_c}\otimes1_{m_c}$ (Satz 2.6). Dann gilt
$\sum_ck_cm_c=8$, $\sum_ck_c^2=\dim\pi(S_3)''=5$ und $\sum_cm_c^2=\dim\pi(S_3)'=20$; beide
Dimensionen berechnet der Code unten mit exakter linearer Algebra (die Kommutante einer Menge $S$
ist der Kern der linearen Abbildung $X\mapsto(AX-XA)_{A\in S}$). Die Zahl $20$ lässt sich auch von
Hand abzählen: $X$ vertauscht genau dann mit allen $\pi(g)$, wenn seine Einträge $X_{xy}$,
$x,y\in\{0,1\}^3$, auf den Bahnen von $S_3$ auf den Paaren $(x,y)$ konstant sind, d. h. auf
Tripeln von Elementen von $\{0,1\}^2$, und es gibt $\binom63=20$ solche Bahnen (Multimengen der
Größe $3$ aus $4$ Elementen). Aus $\sum_ck_c^2=5$ folgt $k=(1,2)$ oder $k=(1,1,1,1,1)$. Im ersten
Fall ergeben $m_1+2m_2=8$ und $m_1^2+m_2^2=20$ die Lösung $m=(4,2)$. Im zweiten Fall müssen fünf
positive ganze Zahlen mit Summe $8$ und Quadratsumme $20$ gleich $m=(4,1,1,1,1)$ sein; dann wäre
aber $\pi(S_3)''$ kommutativ, während die Transpositionen $(12)$ und $(23)$ nicht vertauschen:
$\pi((12))\pi((23))|001\rangle=|100\rangle$, aber $\pi((23))\pi((12))|001\rangle=|010\rangle$.
Also ist $k=(1,2)$ und $m=(4,2)$. Die eindimensionale Darstellung $\rho_1$ ist die triviale oder
die Signumsdarstellung. Sie ist nicht das Signum: Jeder Basisvektor $|b_1b_2b_3\rangle$ hat zwei
gleiche Bits und wird daher von einer Transposition $t$ festgehalten, und
$A=\sum_g\operatorname{sgn}(g)\pi(g)$, das $A\pi(t)=-A$ erfüllt, liefert
$A|b_1b_2b_3\rangle=A\pi(t)|b_1b_2b_3\rangle=-A|b_1b_2b_3\rangle$; also ist $A=0$, und die
Signumsdarstellung tritt nicht auf. Somit ist $\rho_1$ trivial, ihre $4$ Kopien spannen den
symmetrischen Unterraum auf, und $\rho_2$ ist eine irreduzible zweidimensionale Darstellung
(Proposition 2.7): die Standarddarstellung. $\square$

```python
import itertools
import torch

torch.manual_seed(0)
dt = torch.complex128

def permutation_operator(perm, n_qubits=3):
    """Unitary on (C^2)^{(x)3} that permutes the tensor factors."""
    dim = 2 ** n_qubits
    P = torch.zeros(dim, dim, dtype=dt)
    for bits in itertools.product([0, 1], repeat=n_qubits):
        src = int(''.join(map(str, bits)), 2)
        dst = int(''.join(str(bits[perm[i]]) for i in range(n_qubits)), 2)
        P[dst, src] = 1
    return P

group = [permutation_operator(p) for p in itertools.permutations(range(3))]

def commutant_basis(ops, tol=1e-10):
    """Orthonormal basis (Hilbert-Schmidt) of {X : X A = A X for all A in ops}."""
    d = ops[0].shape[0]
    eye = torch.eye(d, dtype=dt)
    # vec(AX - XA) = (A (x) 1 - 1 (x) A^T) vec(X) for row-major vec
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    _, s, Vh = torch.linalg.svd(L)
    null = Vh[(s > tol).sum():].conj().resolve_conj()   # rows span the null space
    return [v.reshape(d, d) for v in null]

def span_dim(ops, tol=1e-10):
    return int((torch.linalg.svdvals(torch.stack([A.reshape(-1) for A in ops])) > tol).sum())

comm = commutant_basis(group)
print(len(comm), span_dim(group), len(commutant_basis(comm)))
assert len(comm) == 20                                   # pi(S3)'  = M_4 (+) M_2
assert span_dim(group) == 5                              # pi(S3)'' = M_1 (+) M_2
assert len(commutant_basis(comm)) == 5                   # double commutant theorem
```

Die Vielfachheiten sind im Spektrum sichtbar. Die Mittelung über die Gruppe,
$X\mapsto\frac1{|G|}\sum_g\pi(g)X\pi(g)^*$, bildet jeden Operator in die Kommutante ab
(Aufgabe 3). Ein zufälliges selbstadjungiertes Element von $\pi(S_3)'$ hat daher (fast
sicher) $4$ einfache Eigenwerte auf dem symmetrischen Unterraum und $2$ zweifach entartete
auf dessen Komplement, weil der Block von $\pi(S_3)'$ dort $1_2\otimes M_2$ ist; der Code
überprüft die entarteten Paare.

```python
def group_average(X):
    return sum(g @ X @ g.conj().T for g in group) / len(group)

H0 = torch.randn(8, 8, dtype=dt)
Hsym = group_average(H0 + H0.conj().T)                  # a random symmetric Hamiltonian
assert all(torch.allclose(g @ Hsym, Hsym @ g) for g in group)

P_sym = sum(group) / len(group)                         # projection onto the symmetric subspace
w, Q = torch.linalg.eigh(P_sym)
Q_rest = Q[:, w.real < 0.5]                             # orthonormal basis of the complement (dim 4)
e_rest = torch.linalg.eigvalsh(Q_rest.conj().T @ Hsym @ Q_rest)
print(e_rest)
assert Q_rest.shape[1] == 4
assert torch.allclose(e_rest[0], e_rest[1]) and torch.allclose(e_rest[2], e_rest[3])
```

## 2.5 Die Algebra in der Bibliothek

Die Bibliothek speichert die abstrakte Algebra $\bigoplus_cM_{k_c}(\mathbb C)$ aus Satz 2.6
auf $\bigoplus_c\mathbb C^{k_c}$, ohne die Vielfachheiten: Ein Operator ist die Liste seiner
Blöcke $A_c$, aufgefüllt auf eine gemeinsame Größe. Die Kommutante aus dem Beispiel ist
`TypeIAlgebra([4, 2], [4, 2])`.

```python
from torch_vn_algebra import TypeIAlgebra, channels

alg = TypeIAlgebra([4, 2], [4, 2], precision='double')    # M_4 (+) M_2
print(alg.k_factors, alg.C)

A = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
print(A.matrix.shape)              # (1, 2, 4, 4): block 1 occupies the top-left 2 x 2 corner
assert A.matrix[0, 1, 2:, :].abs().max() == 0 and A.matrix[0, 1, :, 2:].abs().max() == 0
```

Das Zentrum besteht aus den Elementen $\sum_cz_cP_c$; `alg.central` konstruiert sie, und
sie kommutieren mit allem:

```python
def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

P0, P1 = alg.central([1.0, 0.0]), alg.central([0.0, 1.0])   # minimal central projections
Z = alg.central([2.0, -1.0])
assert P0.is_projection and P1.is_projection
assert close(P0 + P1, alg.identity())
assert close(Z @ A, A @ Z)
assert close(A @ P0 + A @ P1, A)                         # A = sum_c A P_c
```

## 2.6 Bedingte Erwartungen

Eine **bedingte Erwartung** von $\mathcal M$ auf eine von-Neumann-Unteralgebra
$\mathcal N\subseteq\mathcal M$ ist eine lineare Abbildung $E:\mathcal M\to\mathcal N$, die
positiv ist ($A\ge0\Rightarrow E(A)\ge0$), $E(B)=B$ für $B\in\mathcal N$ erfüllt und eine
$\mathcal N$-Bimodulabbildung ist: $E(B_1AB_2)=B_1E(A)B_2$ für $B_1,B_2\in\mathcal N$.

**Proposition 2.8.** Für $\mathcal M=\bigoplus_cM_{k_c}$ und $\mathcal N=Z(\mathcal M)$ ist
$$
E(A)=\sum_c\frac{\operatorname{Tr}A_c}{k_c}\,P_c
$$
eine bedingte Erwartung. Sie erhält jede Spur der Form
$\sum_cw_c\operatorname{Tr}A_c$.

*Beweis.* Die Linearität ist klar; aus $A\ge0$ folgt $\operatorname{Tr}A_c\ge0$, also
$E(A)\ge0$. Für $B=\sum_cb_cP_c$ ist $E(B)=\sum_cb_ck_cP_c/k_c=B$. Für $B_i=\sum_cb^{(i)}_cP_c$ ist
$\operatorname{Tr}(B_1AB_2)_c=b^{(1)}_cb^{(2)}_c\operatorname{Tr}A_c$, woraus die
Bimoduleigenschaft folgt. Schließlich ist $\operatorname{Tr}(E(A)_c)=\operatorname{Tr}A_c$ für jedes $c$. $\square$

In der Bibliothek ist diese Abbildung der Kanal `channels.center_expectation(alg)`. Wendet
man ihn auf Operatoren an, erhält man den klassischen (über die Sektoren gemittelten) Teil:

```python
E = channels.center_expectation(alg)
EA = E(A)
assert close(E(EA), EA)                                  # idempotent
assert close(E(Z @ A @ P1), Z @ EA @ P1)                 # bimodule property
A2 = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
assert close(EA @ A2, A2 @ EA)                           # E(A) commutes with a fresh element
for c in range(alg.C):                                   # Tr E(A)_c = Tr A_c
    assert abs(EA.matrix[0, c].diagonal().sum() - A.matrix[0, c].diagonal().sum()) < 1e-10
```

Die Gruppenmittelung aus Abschnitt 2.4 ist eine bedingte Erwartung auf die Kommutante
$\pi(G)'$ innerhalb von $M_n(\mathbb C)$, und sie erhält $\operatorname{Tr}$ (Aufgabe 3).
(Bedingte Erwartungen sind im Allgemeinen nicht eindeutig: Auf $\mathbb C1$ liefert jeder
Zustand $\varphi$ eine, nämlich $A\mapsto\varphi(A)1$.) Beides sind Beispiele derselben Idee:
die Projektion auf diejenigen Operatoren, die eine gewisse Struktur „nicht sehen“.

> **Physikalische Bedeutung.** Auf der vollen Matrixalgebra $M_D(\mathbb C)$, $D=\sum_ck_c$, ist
> die Abbildung $X\mapsto\sum_cP_cXP_c$ die bedingte Erwartung auf $\mathcal M$: Sie löscht die
> Kohärenzen zwischen Superauswahlsektoren, und genau das geschieht, wenn man die Ladung (den
> Sektor) misst, ohne das Ergebnis zu notieren. Die bedingte Erwartung $E$ auf das Zentrum geht
> weiter: Sie behält nur die klassische Information, das Gewicht $p_c$ jedes Sektors, und ersetzt
> den Zustand innerhalb jedes Sektors durch den maximal gemischten. Ein Zustand $\omega$ und der
> Zustand $\omega\circ E$ stimmen auf jeder Observablen überein, die nur fragt: „Welcher Sektor?“

## 2.7 Womit man vergleichen sollte

Die Beispiele aus Abschnitt 2.4 lassen sich mit jedem Paket für lineare Algebra
überprüfen; Software für Darstellungstheorie (GAP, SageMath) liefert die Zerlegung
$\pi\cong\bigoplus\rho_\lambda\otimes1_{m_\lambda}$ symbolisch. Programme zur exakten
Diagonalisierung in der Physik (QuSpin, die Symmetriemodule von QuTiP) nutzen dieselbe
Blockstruktur, um Hamilton-Operatoren Sektor für Sektor zu diagonalisieren. Die Bibliothek
nimmt die Struktur als Ausgangspunkt: Sind die Blöcke einmal bekannt, so bleiben alle
weiteren Rechnungen – Funktionen, Spuren, Zustände, Kanäle, Dynamik – innerhalb der Algebra.

## Aufgaben

1. Man zeige, dass das Zentrum von $M_k(\mathbb C)$ gleich $\mathbb C1$ ist (man verwende die
   Matrixeinheiten $|i\rangle\langle j|$), und dass $Z(\bigoplus_cM_{k_c})$ die Dimension $C$ hat.
2. Sei $S$ die zyklische Verschiebung auf $\mathbb C^3$, $Se_i=e_{i+1\bmod3}$. Man bestimme
   $\{S\}'$ und $\{S\}''$ und überprüfe die Dimensionen numerisch.
3. Für eine unitäre Darstellung $\pi$ einer endlichen Gruppe zeige man, dass
   $E(X)=\frac1{|G|}\sum_g\pi(g)X\pi(g)^*$ eine bedingte Erwartung von
   $M_n(\mathbb C)$ auf $\pi(G)'$ ist, die $\operatorname{Tr}$ erhält.

<details>
<summary>Lösungen</summary>

1. Kommutiert $X$ mit allen $|i\rangle\langle j|$, so gilt
   $X|i\rangle\langle j|=|i\rangle\langle j|X$; wendet man beide Seiten auf $e_j$ an und
   vergleicht, so erhält man $Xe_i=X_{jj}e_i$ für alle $i,j$, also $X=X_{11}1$. In einer
   direkten Summe ist ein Element $\bigoplus X_c$ genau dann zentral, wenn jedes $X_c$ es
   ist; das Zentrum ist also $\{\bigoplus_cz_c1_{k_c}\}$ und hat die Dimension $C$.

2. $S$ ist unitär mit den Eigenwerten $1,\omega,\omega^2$, $\omega=e^{2\pi i/3}$, jeder einfach.
   Ein mit $S$ vertauschbarer Operator erhält dessen Eigenräume, ist also in der
   Eigenbasis diagonal: $\{S\}'\cong\mathbb C^3$, die zirkulanten Matrizen. Dann ist auch
   $\{S\}''\cong\mathbb C^3$, denn in der Eigenbasis von $S$ ist eine Matrix, die mit allen
   Diagonalmatrizen kommutiert, diagonal (man vergleiche die Einträge von $XD$ und $DX$
   für ein $D$ mit paarweise verschiedenen Diagonaleinträgen); sie wird von $1,S,S^2$
   aufgespannt (Vandermonde: Die Eigenwerte sind paarweise verschieden).

```python
S = torch.roll(torch.eye(3, dtype=dt), 1, dims=0)
assert len(commutant_basis([S])) == 3
assert len(commutant_basis(commutant_basis([S]))) == 3
assert span_dim([torch.eye(3, dtype=dt), S, S @ S]) == 3
```

3. $E$ ist linear und positiv (eine Summe von Abbildungen $X\mapsto UXU^*$). Für $h\in G$ ist
   $\pi(h)E(X)\pi(h)^*=\frac1{|G|}\sum_g\pi(hg)X\pi(hg)^*=E(X)$, also $E(X)\in\pi(G)'$. Ist
   $X\in\pi(G)'$, so ist jeder Summand gleich $X$, also $E(X)=X$. Für $B_1,B_2\in\pi(G)'$ ist
   $\pi(g)B_1XB_2\pi(g)^*=B_1\pi(g)X\pi(g)^*B_2$, woraus die Bimoduleigenschaft folgt.
   Schließlich ist $\operatorname{Tr}(\pi(g)X\pi(g)^*)=\operatorname{Tr}X$ für jedes $g$, also
   $\operatorname{Tr}E(X)=\operatorname{Tr}X$.

</details>

**Weiter:** [Kapitel 3](03_traces_and_states.md) klassifiziert die Spuren auf
$\bigoplus_cM_{k_c}$, führt Zustände und Dichten ein und beweist die später verwendeten
Entropie- und Spurungleichungen.
