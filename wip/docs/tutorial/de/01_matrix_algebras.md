[English](../01_matrix_algebras.md) · [Español](../es/01_matrix_algebras.md) · [Français](../fr/01_matrix_algebras.md) · **Deutsch** · [中文](../zh/01_matrix_algebras.md) · [日本語](../ja/01_matrix_algebras.md) · [Русский](../ru/01_matrix_algebras.md)

# 1. Matrixalgebren und Funktionalkalkül

Dieses Kapitel legt die Bezeichnungen fest und beweist die Tatsachen über
$M_n(\mathbb C)$, die alle späteren Kapitel verwenden: den Spektralsatz, den
Funktionalkalkül, die Polarzerlegung und die durch Singulärwerte definierten Normen. Auf
jedes Ergebnis folgen der zugehörige Bibliothekscode und eine numerische Überprüfung.

**Sie lernen**, wie man Operatoren mit vorgegebenem Spektrum konstruiert, wie man eine
Funktion auf einen Operator anwendet, wie Stapel (Batches) von Operatoren gespeichert
werden und wie man eine Rechnung an einem exakten Ergebnis überprüft.

**Konventionen.** Skalarprodukte $\langle x,y\rangle=\sum_i\bar x_iy_i$ sind linear im
zweiten Argument. $A^*$ ist die konjugiert Transponierte, also $\langle x,Ay\rangle=\langle A^*x,y\rangle$.
$e_1,\dots,e_n$ ist die Standardbasis, $1$ die Einheitsmatrix, $\operatorname{Tr}$ die
gewöhnliche Spur.

## 1.1 Die Algebra $M_n(\mathbb C)$

$M_n(\mathbb C)$ ist die Algebra der komplexen $n\times n$-Matrizen mit der Involution
$A\mapsto A^*$ und der Operatornorm $\|A\|=\sup_{\|x\|=1}\|Ax\|$. Ein Element $A$ heißt

| Bezeichnung | Definition |
|---|---|
| normal | $A^*A=AA^*$ |
| selbstadjungiert | $A=A^*$ |
| positiv, $A\ge0$ | $\langle x,Ax\rangle\ge0$ für alle $x$ |
| unitär | $A^*A=AA^*=1$ |
| Projektion | $A=A^*=A^2$ |

Wir schreiben $A\le B$, falls $B-A\ge0$, und $A>0$, falls $A\ge0$ und $A$ invertierbar ist.

**Lemma 1.1 (C\*-Identität).** $\|A^*\|=\|A\|$ und $\|A^*A\|=\|A\|^2$.

*Beweis.* $\|Ax\|=\sup_{\|y\|=1}|\langle y,Ax\rangle|=\sup_{\|y\|=1}|\langle A^*y,x\rangle|$,
also $\|A\|=\sup_{\|x\|=\|y\|=1}|\langle A^*y,x\rangle|=\|A^*\|$. Ferner ist
$\|Ax\|^2=\langle x,A^*Ax\rangle\le\|A^*A\|$ für $\|x\|=1$, folglich $\|A\|^2\le\|A^*A\|\le\|A^*\|\,\|A\|=\|A\|^2$. $\square$

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). In der
> Quantenmechanik wird ein System mit $n$ unterscheidbaren Konfigurationen durch $\mathbb C^n$
> beschrieben. Selbstadjungierte Matrizen sind *Observablen* (Energie, Spin, Ort auf einem
> Gitter); ihre Eigenwerte sind die möglichen Messergebnisse. Projektionen sind Ja/Nein-Fragen. Unitäre Matrizen sind
> Symmetrien und Zeitentwicklungen, $e^{-itH}$ zu einem Hamilton-Operator $H$. Positive
> Matrizen mit Spur eins sind *Zustände* (Dichtematrizen, Kapitel 3).

## 1.2 Der Spektralsatz

**Lemma 1.2 (Schur).** Zu jedem $A\in M_n(\mathbb C)$ gibt es ein unitäres $U$, sodass
$U^*AU$ eine obere Dreiecksmatrix ist.

*Beweis.* Induktion nach $n$; der Fall $n=1$ ist trivial. Das charakteristische Polynom
von $A$ hat eine Nullstelle $\lambda$ (Fundamentalsatz der Algebra), also gilt
$Av=\lambda v$ für einen Einheitsvektor $v$. Man ergänze $v$ zu einer Orthonormalbasis
$v,w_2,\dots,w_n$ und sei $V$ die unitäre Matrix mit diesen Spalten. Dann ist
$V^*AV=\begin{pmatrix}\lambda&*\\0&A'\end{pmatrix}$ mit $A'\in M_{n-1}(\mathbb C)$.
Nach Induktionsvoraussetzung ist $U'^*A'U'$ für ein unitäres $U'$ eine obere
Dreiecksmatrix, und $U=V\,(1\oplus U')$ leistet das Gewünschte. $\square$

**Satz 1.3 (Spektralsatz).** $A$ ist genau dann normal, wenn
$A=U\,\mathrm{diag}(\lambda_1,\dots,\lambda_n)\,U^*$ mit einem unitären $U$ gilt. Die
$\lambda_i$ sind die Eigenwerte von $A$, mit Vielfachheit gezählt.

*Beweis.* Ist $A=U\Lambda U^*$ mit diagonalem $\Lambda$, so gilt
$A^*A=U\bar\Lambda\Lambda U^*=U\Lambda\bar\Lambda U^*=AA^*$. Sei umgekehrt $A$
normal und $T=U^*AU$ eine obere Dreiecksmatrix (Lemma 1.2); $T$ ist normal. Vergleich
der $(1,1)$-Einträge: $(T^*T)_{11}=|t_{11}|^2$ und $(TT^*)_{11}=\sum_j|t_{1j}|^2$, also
$t_{1j}=0$ für $j>1$. Damit ist $T=t_{11}\oplus T'$ mit einer normalen oberen
Dreiecksmatrix $T'$, und Induktion zeigt, dass $T$ diagonal ist. Die Diagonale von
$U^*AU$ besteht aus den Nullstellen von $\det(A-\lambda)=\det(T-\lambda)$. $\square$

**Korollar 1.4.** Sei $A$ normal mit Eigenwerten $\lambda_i$. Dann gilt:

1. $\|A\|=\max_i|\lambda_i|$;
2. $A$ ist genau dann selbstadjungiert, wenn alle $\lambda_i\in\mathbb R$; genau dann unitär, wenn alle $|\lambda_i|=1$;
   genau dann eine Projektion, wenn alle $\lambda_i\in\{0,1\}$;
3. $A\ge0$ genau dann, wenn $A$ selbstadjungiert ist und alle $\lambda_i\ge0$ sind. Außerdem ist jede
   positive Matrix selbstadjungiert, also normal; daher charakterisiert (3) die Positivität unter
   *allen* Matrizen, und der folgende Funktionalkalkül ist auf jedes $A\ge0$ anwendbar.

*Beweis.* Mit $y=U^*x$ ist $\|Ax\|^2=\sum_i|\lambda_i|^2|y_i|^2$ und $\|y\|=\|x\|$,
woraus (1) folgt, mit Gleichheit für $x=Ue_i$, $|\lambda_i|$ maximal. In (2) ist $A^*=U\bar\Lambda U^*$,
$A^*A=U|\Lambda|^2U^*$ und $A^2=U\Lambda^2U^*$, und die Bedingungen gelten genau dann, wenn
sie eintragsweise für $\Lambda$ gelten. Zu (3) zunächst: Jedes positive $A$ ist
selbstadjungiert, denn
$\langle x,(A-A^*)x\rangle=\langle x,Ax\rangle-\overline{\langle x,Ax\rangle}=0$ für alle $x$,
und eine Matrix $D$ mit $\langle x,Dx\rangle=0$ für alle $x$ verschwindet (Polarisation:
$4\langle x,Dy\rangle=\sum_{k=0}^3 i^{-k}\langle x+i^ky,D(x+i^ky)\rangle$). Ferner ist
$\langle x,Ax\rangle=\sum_i\lambda_i|y_i|^2$, und dies ist genau dann für alle $x$ $\ge0$, wenn alle $\lambda_i\ge0$. $\square$

## 1.3 Der Funktionalkalkül

Sei $A=U\Lambda U^*$ normal und $f$ eine komplexwertige Funktion auf $\operatorname{spec}A$. Man definiert
$$
f(A)=U\,\mathrm{diag}\big(f(\lambda_1),\dots,f(\lambda_n)\big)\,U^*.
$$

**Satz 1.5.** (a) $f(A)$ hängt nicht von der Wahl von $U$ ab. (b) Für ein
Polynom $p$ ist $p(A)$ das gewöhnliche Polynom in $A$. (c) $f\mapsto f(A)$ ist
linear und multiplikativ, $\bar f(A)=f(A)^*$, und $\operatorname{spec}f(A)=f(\operatorname{spec}A)$.
(d) Kommutiert $B$ mit $A$, so kommutiert $B$ mit $f(A)$. (e) Ist $g$ auf
$f(\operatorname{spec}A)$ definiert, so gilt $g(f(A))=(g\circ f)(A)$.

*Beweis.* (b) $A^k=U\Lambda^kU^*$, also $p(A)=Up(\Lambda)U^*$ wegen der Linearität. (a) Das
Spektrum ist endlich, daher liefert Lagrange-Interpolation ein Polynom $p$ mit $p=f$ auf
$\operatorname{spec}A$: Sind $\mu_1,\dots,\mu_r$ die verschiedenen Eigenwerte, so nehme man
$p(t)=\sum_if(\mu_i)\prod_{j\neq i}(t-\mu_j)/(\mu_i-\mu_j)$. Dann ist $U\,f(\Lambda)\,U^*=U\,p(\Lambda)\,U^*=p(A)$ für *jedes*
unitäre $U$, das $A$ diagonalisiert, und $p(A)$ ist ohne Bezug auf $U$ definiert. (c) gilt
für die Diagonalmatrizen $f(\Lambda)$ und bleibt unter Konjugation mit $U$ erhalten; das
Spektrum von $Uf(\Lambda)U^*$ ist die Diagonale von $f(\Lambda)$. (d) $B$ kommutiert mit jedem
Polynom in $A$, insbesondere mit $p(A)=f(A)$. (e) $f(A)=Uf(\Lambda)U^*$ ist normal und wird
durch dasselbe $U$ diagonalisiert, also gilt nach (a) $g(f(A))=U\,g(f(\Lambda))\,U^*=(g\circ f)(A)$. $\square$

**Satz 1.6 (Quadratwurzel).** Zu $A\ge0$ gibt es genau ein $B\ge0$ mit
$B^2=A$, nämlich das durch den Funktionalkalkül definierte $B=\sqrt A$.

*Beweis.* $\sqrt A\ge0$ und $(\sqrt A)^2=A$ nach Satz 1.5(c) und Korollar 1.4(3).
Sei $B\ge0$ mit $B^2=A$, $B=V\,\mathrm{diag}(\mu_j)\,V^*$, $\mu_j\ge0$. Dann ist
$A=V\,\mathrm{diag}(\mu_j^2)\,V^*$, also $\operatorname{spec}A=\{\mu_j^2\}$. Sei $p$ ein
Polynom mit $p(t)=\sqrt t$ auf $\operatorname{spec}A$. Dann gilt
$p(A)=\sqrt A$ und $p(A)=p(B^2)=V\,\mathrm{diag}(p(\mu_j^2))\,V^*=V\,\mathrm{diag}(\mu_j)\,V^*=B$. $\square$

**Korollar 1.7.** $A\ge0$ genau dann, wenn $A=T^*T$ für ein $T$ gilt. (Ist $A\ge0$, so nehme man $T=\sqrt A$;
umgekehrt ist $\langle x,T^*Tx\rangle=\|Tx\|^2\ge0$.)

## 1.4 Polarzerlegung und Singulärwerte

Für jedes $A$ ist die Matrix $A^*A$ positiv (Korollar 1.7); ihre Quadratwurzel
$|A|=(A^*A)^{1/2}$ heißt der **Betrag** von $A$.

**Satz 1.8 (Polarzerlegung).** Jedes $A\in M_n(\mathbb C)$ lässt sich als
$A=W|A|$ mit unitärem $W$ schreiben.

*Beweis.* Für jedes $x$ ist $\||A|x\|^2=\langle x,|A|^2x\rangle=\langle x,A^*Ax\rangle=\|Ax\|^2$.
Daher ist $W_0:|A|x\mapsto Ax$ eine wohldefinierte lineare Isometrie vom Bild
$R_1$ von $|A|$ auf das Bild $R_2$ von $A$: Aus $|A|x=|A|x'$ folgt
$\|A(x-x')\|=\||A|(x-x')\|=0$. Insbesondere ist $\dim R_1=\dim R_2$, also haben auch die
orthogonalen Komplemente gleiche Dimension; sei $W_1:R_1^\perp\to R_2^\perp$ eine beliebige
Isometrie zwischen ihnen (man bilde eine Orthonormalbasis auf eine Orthonormalbasis ab). Sei
$W=W_0\oplus W_1$, d. h. $W(u_1+u_2)=W_0u_1+W_1u_2$ für $u_1\in R_1$, $u_2\in R_1^\perp$. Wegen
$W_0u_1\in R_2\perp R_2^\perp\ni W_1u_2$ liefert der Satz des Pythagoras $\|Wu\|^2=\|u_1\|^2+\|u_2\|^2=\|u\|^2$
für alle $u$. Also ist $\langle u,(W^*W-1)u\rangle=0$ für alle $u$ und damit $W^*W=1$ durch Polarisation
(Beweis von Korollar 1.4(3)); folglich ist $W$ invertierbar mit $W^{-1}=W^*$, also auch $WW^*=1$: $W$ ist
unitär, und $W|A|x=W_0|A|x=Ax$. $\square$

Die Eigenwerte $s_1\ge\dots\ge s_n\ge0$ von $|A|$ heißen die **Singulärwerte** von
$A$. Durchgehend werden drei Normen verwendet:
$$
\|A\|=s_1,\qquad \|A\|_2=\Big(\sum_is_i^2\Big)^{1/2}=(\operatorname{Tr}A^*A)^{1/2},
\qquad\|A\|_1=\sum_is_i=\operatorname{Tr}|A|.
$$
Die erste Identität gilt, weil $\|Ax\|=\||A|x\|$ für jedes $x$ ist, also
$\|A\|=\||A|\|=s_1$ nach Korollar 1.4(1); die zweite, weil
$\operatorname{Tr}A^*A=\operatorname{Tr}|A|^2=\sum_is_i^2$.

**Satz 1.9 (Dualität der Spurnorm).**
$\|A\|_1=\max_{U\text{ unitär}}|\operatorname{Tr}(AU)|$, und das Maximum wird bei
$U=W^*$ angenommen, wobei $A=W|A|$. Insbesondere ist $|\operatorname{Tr}A|\le\|A\|_1$.

*Beweis.* $\operatorname{Tr}(AW^*)=\operatorname{Tr}(W|A|W^*)=\operatorname{Tr}|A|$. Für jedes unitäre $U$
gilt wegen der Zyklizität und der Cauchy-Schwarz-Ungleichung für $\langle X,Y\rangle_{HS}=\operatorname{Tr}(X^*Y)$
$$
|\operatorname{Tr}(AU)|=|\operatorname{Tr}(|A|^{1/2}\,UW\,|A|^{1/2})|
=|\langle|A|^{1/2},\,UW|A|^{1/2}\rangle_{HS}|
\le\||A|^{1/2}\|_2\,\|UW|A|^{1/2}\|_2=\operatorname{Tr}|A|,
$$
da $\|UW X\|_2=\|X\|_2$ für unitäres $UW$ und $\||A|^{1/2}\|_2^2=\operatorname{Tr}|A|$. $\square$

Insbesondere sind die drei obigen Größen Normen. $\|\cdot\|$ ist eine Operatornorm;
$\|\cdot\|_2$ ist die Norm des Skalarprodukts $\langle X,Y\rangle_{HS}$; und für
$\|\cdot\|_1$ liefert Satz 1.9
$\|A+B\|_1=\max_U|\operatorname{Tr}(AU)+\operatorname{Tr}(BU)|\le\|A\|_1+\|B\|_1$ und
$\|cA\|_1=\max_U|c\operatorname{Tr}(AU)|=|c|\,\|A\|_1$, während
aus $\|A\|_1=0$ zwingend $|A|=0$ und $A=W|A|=0$ folgt. Für $P\ge0$ ist $|P|=(P^2)^{1/2}=P$ nach Satz 1.6,
also $\|P\|_1=\operatorname{Tr}P$.

## 1.5 Operatoren in der Bibliothek

Die Bibliothek arbeitet mit Algebren $\bigoplus_cM_{k_c}(\mathbb C)$ (Kapitel 2). Ein einzelner
Block ist der Fall $C=1$: `TypeIAlgebra([n], [n])` ist $M_n(\mathbb C)$. (Die beiden Listen
sind die nominellen und die tatsächlichen Blockgrößen; vorerst stimmen sie überein.)

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra

torch.manual_seed(0)
M3 = TypeIAlgebra([3], [3], precision='double')       # the algebra M_3(C), complex128

# a positive operator with eigenvalues 1, 2, 4 and Haar-random eigenvectors
A = M3.operator_from_eigenvalues(lambda d: torch.tensor([[1.0, 2.0, 4.0]]),
                                 force_positive=True)
print(A.matrix.shape)          # (batch, C, k_max, k_max) = (1, 1, 3, 3)
print(A.eigenvalues()[0])      # one tensor per block, eigenvalues in ascending order
```

Jeder Operator wird als ein Tensor der Gestalt `(batch, C, k_max, k_max)` gespeichert: ein
Stapel von Operatoren, jeder eine Liste von `C` Blöcken. Hier ist `batch = C = 1`. Die
Konstruktion ist genau Satz 1.3, rückwärts gelesen: $A=U\Lambda U^*$ mit Haar-verteiltem
$U$ (die Bibliothek erzeugt $U$ mit der QR-Konstruktion aus Satz 5.3). Das Flag
`force_positive=True` wird an den vorgegebenen Eigenwerten überprüft (die Bibliothek meldet einen
Fehler, wenn einer von ihnen kleiner als $-10^{-12}$ ist); nach Korollar 1.4(3) ist das Ergebnis
dann positiv (exakt, wenn die vorgegebenen Eigenwerte $\ge0$ sind, wie hier; im Allgemeinen bis auf
diese Toleranz), und die Bibliothek vermerkt dies, statt es erneut zu testen. Aus demselben Grund
gibt `A.eigenvalues()` oben das vorgegebene Spektrum sortiert zurück, ohne etwas zu berechnen; das
ist keine Überprüfung. Die Überprüfung von $\operatorname{spec}f(A)=f(\operatorname{spec}A)$ weiter
unten rechnet dagegen tatsächlich: `apply_function` diagonalisiert die Matrix von $A$, und die
Eigenwerte der resultierenden Matrix $e^A$ werden von Grund auf neu berechnet.

Der Funktionalkalkül für selbstadjungierte Operatoren ist `apply_function` (die Methode
diagonalisiert jeden Block mit `eigh`). Die Funktion $f$ muss elementweise wirken: Sie erhält die
Eigenwerte aller Sektoren einer Größe $k$ als einen Tensor der Form `(batch, n·k)`, wobei $n$ die
Zahl dieser Sektoren ist, und gibt einen Tensor derselben Form zurück. Die gebräuchlichen
Funktionen haben eigene Methoden. Wir überprüfen die Sätze 1.5 und 1.6:

```python
def close(X, Y, tol=1e-10):
    """Maximal entry of |X - Y| for two operators (or tensors) is below tol."""
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

R = A.sqrt()
assert close(R @ R, A)                                  # Theorem 1.6
assert close(A.log().expm(), A)                         # exp(log A) = A, Theorem 1.5(e)
assert close(A.apply_function(lambda w: w ** 3), A @ A @ A)   # Theorem 1.5(b)
assert close(A.power(-1), A.inverse())                  # f(t) = 1/t
spec = A.apply_function(torch.exp).eigenvalues()[0]     # spec f(A) = f(spec A)
assert torch.allclose(spec, torch.exp(torch.tensor([[1.0, 2.0, 4.0]], dtype=spec.dtype)))
```

`R @ R` berechnet noch nichts: Produkte, Summen und Funktionen bauen eine Vorschrift auf,
die erst ausgewertet wird, wenn `.matrix` oder eine Zahl wie etwa eine Spur angefordert wird.

## 1.6 Normen und die Dualität der Spurnorm

Eine allgemeine (nicht normale) Matrix wird über `alg.operator(tensor)` eingegeben. Wir
überprüfen die drei Normen an direkt mit PyTorch berechneten Singulärwerten, die
Polarzerlegung (Satz 1.8) und die Dualität (Satz 1.9):

```python
G = torch.randn(1, 1, 3, 3, dtype=torch.complex128)    # Ginibre matrix: no symmetry at all
X = M3.operator(G)
s = torch.linalg.svdvals(G[0, 0])                       # singular values, descending

assert abs(X.operator_norm().item() - s[0].item()) < 1e-12
assert abs(X.frobenius_norm().item() - s.norm().item()) < 1e-12
assert abs(X.trace_norm().item() - s.sum().item()) < 1e-12

absX = X.abs()
W = X @ absX.inverse()                                  # |X| is invertible almost surely, see below
assert close(W.adjoint() @ W, M3.identity())            # W is unitary
assert abs((X @ W.adjoint()).trace.item() - s.sum().item()) < 1e-10    # Tr(X W*) = ||X||_1
assert abs(X.trace.item()) <= X.trace_norm().item()     # |Tr X| <= ||X||_1
```

`torch.randn` mit komplexem dtype erzeugt unabhängige Einträge mit unabhängigen Real- und
Imaginärteilen $\sim N(0,\tfrac12)$; `G` ist also eine Ginibre-Matrix, und diese hat eine Dichte
(Proposition A.16(b)). $|X|$ ist fast sicher invertierbar:
$(\det|X|)^2=\det(|X|^2)=\det(X^*X)=|\det X|^2$, und $|\det X|^2$ ist ein Polynom in den Real- und
Imaginärteilen der Einträge, das bei $X=1$ den Wert $1$ hat; also verschwindet es mit
Wahrscheinlichkeit null (Korollar A.7(b); dies ist auch der erste Schritt von Satz 5.3).

## 1.7 Stapel: viele Operatoren auf einmal

Die Stapelachse ist der Grund, weshalb es die Bibliothek gibt: Zufallsmatrixtheorie,
Monte-Carlo-Verfahren und Optimierung mit mehreren Startpunkten benötigen Tausende von
Operatoren, und ein einziger gestapelter Aufruf ersetzt eine Python-Schleife über sie.
Als Beispiel überprüfen wir ein exaktes Moment des **Gaußschen unitären Ensembles** (GUE).

Die Einträge von $G$ seien unabhängig mit unabhängigen Real- und Imaginärteilen
$\sim N(0,\tfrac12)$, also $\mathbb E|G_{ij}|^2=1$, und sei $H=(G+G^*)/\sqrt2$.

**Proposition 1.10.** $\mathbb E\,\operatorname{Tr}H^4=2n^3+n$. Folglich gilt für
$X=H/\sqrt n$: $\ \mathbb E\,\tfrac1n\operatorname{Tr}X^4=2+\tfrac1{n^2}$.

*Beweis.* *Formel von Isserlis (Wick)* (Satz A.24, Fall $m=4$). Die Real- und Imaginärteile
der Einträge von $G$ sind unabhängige zentrierte Gaußsche Zufallsvariablen, also die Koordinaten
eines einzigen zentrierten Gauß-Vektors, und die Einträge von $H$ und von $\bar H$ sind komplexe
Linearkombinationen davon. Satz A.24 liefert dann für je vier solche Zufallsvariablen
$Z_1,\dots,Z_4$ (Wiederholungen erlaubt)
$$
\mathbb E[Z_1Z_2Z_3Z_4]=\mathbb E[Z_1Z_2]\,\mathbb E[Z_3Z_4]+\mathbb E[Z_1Z_3]\,\mathbb E[Z_2Z_4]+\mathbb E[Z_1Z_4]\,\mathbb E[Z_2Z_3],
$$
die Summe über die drei Paarungen von $\{1,2,3,4\}$ der Produkte der Paarerwartungswerte (ohne
komplexe Konjugation). Der Beweis in Anhang A vergleicht die Potenzreihe der momenterzeugenden
Funktion mit den Gaußschen Momenten und differenziert nicht unter dem Erwartungswert.

*Das Moment.* Die Paarerwartungswerte sind
$\mathbb E[H_{ij}H_{kl}]=\tfrac12\mathbb E[(G_{ij}+\bar G_{ji})(G_{kl}+\bar G_{lk})]=\delta_{il}\delta_{jk}$,
denn $\mathbb E[G_{ij}G_{kl}]=0$ und $\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}$ (Proposition A.16(a);
für $(i,j)\neq(k,l)$ sind die Einträge unabhängig und zentriert, und für $G_{ij}=x+iy$ ist
$\mathbb E[G_{ij}^2]=\mathbb E[x^2-y^2]+2i\,\mathbb E[xy]=0$ und $\mathbb E|G_{ij}|^2=\tfrac12+\tfrac12=1$).
In $\operatorname{Tr}H^4=\sum_{ijkl}H_{ij}H_{jk}H_{kl}H_{li}$ liefern die Paarungen:
$\{12\}\{34\}$: $\delta_{ik}\delta_{ki}$, mit Summe $n^3$ (freie Indizes $i=k$, $j$, $l$);
$\{14\}\{23\}$: $\delta_{jl}\delta_{jl}$, mit Summe $n^3$;
$\{13\}\{24\}$: $\delta_{il}\delta_{jk}\delta_{ji}\delta_{kl}$, was $i=j=k=l$ erzwingt, mit Summe $n$.
Division durch $n\cdot n^2$ ergibt die zweite Aussage. $\square$

Die beiden *nichtkreuzenden* Paarungen $\{12\}\{34\}$ und $\{14\}\{23\}$ (die Paare lassen sich als
sich nicht schneidende Bögen über den Punkten $1,2,3,4$ zeichnen) liefern den führenden Term $2$, für
jedes $n$; die kreuzende Paarung $\{13\}\{24\}$ liefert die Korrektur $1/n^2$.

```python
n, B = 4, 200_000
G = torch.complex(torch.randn(B, 1, n, n), torch.randn(B, 1, n, n)).to(torch.complex128) / math.sqrt(2)
M4 = TypeIAlgebra([n], [n], precision='double')
X = M4.operator((G + G.conj().transpose(-2, -1)) / math.sqrt(2 * n), is_self_adjoint=True)

X2 = X @ X
m4 = (X2 @ X2).trace.real / n                           # shape (B,): one value per sample
mean, err = m4.mean().item(), m4.std().item() / math.sqrt(B)
print(f"E tr X^4 / n = {mean:.4f} +- {err:.4f}   (exact {2 + 1 / n**2:.4f})")
assert abs(mean - (2 + 1 / n**2)) < 5 * err
```

Eine Monte-Carlo-Überprüfung ist ein statistischer Test: Der Schätzwert muss mit dem
exakten Wert bis auf wenige Standardfehler übereinstimmen. Die Zufallsvariable
$\tfrac1n\operatorname{Tr}X^4$ ist ein Polynom in Gaußschen Zufallsvariablen und hat daher endliche
Varianz $\sigma^2$ (Lemma A.23: alle Gaußschen Momente sind endlich; Produkte behandelt man mit
$|y_1\cdots y_r|\le\sum_i|y_i|^r$), und nach der Tschebyschow-Ungleichung weicht der Mittelwert von $B$
unabhängigen Stichproben mit Wahrscheinlichkeit höchstens $1/25$ um mehr als $5\sigma/\sqrt B$ vom
exakten Wert ab. Der Code verwendet die Stichprobenstandardabweichung anstelle von $\sigma$, und der
feste Seed macht die Überprüfung reproduzierbar. Nach diesem Muster verläuft in diesem Kurs jede
Überprüfung an bewiesenen Ergebnissen: eine exakte Aussage, ein numerischer Schätzwert und eine
explizite Toleranz.

## 1.8 Womit man vergleichen sollte

Für einzelne Matrizen sind `numpy.linalg` und `scipy.linalg` (`sqrtm`, `expm`, `logm`,
`polar`) die Referenzimplementierungen, und die obigen Überprüfungen lassen sich mit ihnen
wiederholen. Die Bibliothek fügt drei Dinge hinzu, die später wichtig werden: Stapelverarbeitung
(Abschnitt 1.7), Operatoren, die eine Blockstruktur tragen (Kapitel 2), und automatisches
Differenzieren (PyTorchs Autograd) durch ihre Operationen hindurch.

## Aufgaben

1. Man überprüfe $(fg)(A)=f(A)g(A)$ für $f=\sin$, $g=\exp$ an einem zufälligen selbstadjungierten $A$.
2. Man zeige mithilfe von Satz 1.9, dass $\|A\|_1=\|A^*\|_1$ und $\|AB\|_1\le\|A\|\,\|B\|_1$.
3. Das **geometrische Mittel** von $A,B>0$ ist $A\#B=A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}$.
   Man beweise, dass $A\#B$ die eindeutige positive Lösung $X$ von $XA^{-1}X=B$ ist, folgere
   $A\#B=B\#A$ und überprüfe beides numerisch.

<details>
<summary>Lösungen</summary>

1. Nach Satz 1.5(c).

```python
Hs = M3.operator_from_eigenvalues(lambda d: torch.randn(1, d), force_self_adjoint=True)
assert close(Hs.apply_function(lambda w: torch.sin(w) * torch.exp(w)),
             Hs.apply_function(torch.sin) @ Hs.apply_function(torch.exp))
```

2. $|\operatorname{Tr}(A^*U)|=|\overline{\operatorname{Tr}(U^*A)}|=|\operatorname{Tr}(AU^*)|$, und
   $U\mapsto U^*$ ist eine Bijektion der unitären Gruppe, also stimmen die Maxima in Satz 1.9
   überein. Für die zweite Behauptung sei $B=V|B|$. Für unitäres $U$ ist
   $|\operatorname{Tr}(ABU)|=|\operatorname{Tr}(|B|^{1/2}UAV|B|^{1/2})|
   \le\||B|^{1/2}\|_2\,\|UAV|B|^{1/2}\|_2\le\|A\|\operatorname{Tr}|B|$, wobei
   $\|YX\|_2\le\|Y\|\,\|X\|_2$ (man wende $\|Yx\|\le\|Y\|\,\|x\|$ auf jede Spalte von $X$ an) und
   $\|UAV\|=\|A\|$ für unitäre $U,V$ verwendet wurden (sie erhalten die Länge von Vektoren).

3. Man setze $T=A^{-1/2}XA^{-1/2}$. Dann gilt $XA^{-1}X=B$ genau dann, wenn $T^2=A^{-1/2}BA^{-1/2}$, und
   $X\ge0$ genau dann, wenn $T\ge0$, da $\langle x,Xx\rangle=\langle A^{1/2}x,TA^{1/2}x\rangle$ und
   $A^{1/2}$ invertierbar ist (hier sind $A^{\pm1/2}$ durch den Funktionalkalkül mit $t^{\pm1/2}$ auf
   $\operatorname{spec}A\subset(0,\infty)$ gegeben, also $A^{-1/2}A^{1/2}=1$ nach Satz 1.5(c)). Außerdem
   ist $A^{-1/2}BA^{-1/2}\ge0$, da $\langle x,A^{-1/2}BA^{-1/2}x\rangle=\langle A^{-1/2}x,BA^{-1/2}x\rangle\ge0$.
   Nach Satz 1.6 ist $T=(A^{-1/2}BA^{-1/2})^{1/2}$ die eindeutige
   positive Lösung, d. h. $X=A\#B$. Jede Lösung $X$ von $XA^{-1}X=B$ ist invertierbar,
   denn $\det(X)^2\det(A)^{-1}=\det B\neq0$. Daher ist $XA^{-1}X=B$ äquivalent zu
   $XB^{-1}X=A$ (aus $A^{-1}=X^{-1}BX^{-1}$ durch Invertieren beider Seiten), und die
   Eindeutigkeit liefert $A\#B=B\#A$.

```python
def geometric_mean(A, B):
    a, ai = A.sqrt(), A.power(-0.5)
    return a @ (ai @ B @ ai).sqrt() @ a

P = M3.operator_from_eigenvalues(lambda d: 0.5 + torch.rand(1, d), force_positive=True)
Q = M3.operator_from_eigenvalues(lambda d: 0.5 + torch.rand(1, d), force_positive=True)
GM = geometric_mean(P, Q)
assert close(GM, geometric_mean(Q, P), 1e-9)
assert close(GM @ P.inverse() @ GM, Q, 1e-9)
```

</details>

**Weiter:** [Kapitel 2](02_von_neumann_algebras.md) ersetzt $M_n$ durch eine direkte Summe von
Matrixalgebren und beweist, dass jede endlichdimensionale von-Neumann-Algebra diese
Gestalt hat.
