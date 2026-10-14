[English](../01_matrix_algebras.md) · [Español](../es/01_matrix_algebras.md) · [Français](../fr/01_matrix_algebras.md) · **Deutsch** · [中文](../zh/01_matrix_algebras.md) · [日本語](../ja/01_matrix_algebras.md) · [Русский](../ru/01_matrix_algebras.md)

# Tutorial 1. Matrixalgebren und Funktionalkalkül

Dieses Tutorial begleitet [Kapitel 1](../../lecture_notes/de/01_matrix_algebras.md) der Vorlesungsnotizen, das den Spektralsatz, den
Funktionalkalkül, die Polarzerlegung und die Dualität der Spurnorm in $M_n(\mathbb C)$ beweist. Hier wird jedes dieser
Ergebnisse mit der Bibliothek überprüft.

**Sie lernen**, wie man Operatoren mit vorgegebenem Spektrum konstruiert, wie man eine
Funktion auf einen Operator anwendet, wie Stapel (Batches) von Operatoren gespeichert
werden und wie man eine Rechnung an einem exakten Ergebnis überprüft.

<a id="s-1-5"></a>

## 1.5 Operatoren in der Bibliothek

*Vorlesungsnotizen:* [Abschnitt 1.2](../../lecture_notes/de/01_matrix_algebras.md#s-1-2) und [Abschnitt 1.3](../../lecture_notes/de/01_matrix_algebras.md#s-1-3): [Satz 1.3](../../lecture_notes/de/01_matrix_algebras.md#r-1-3), [Korollar 1.4](../../lecture_notes/de/01_matrix_algebras.md#r-1-4), [Satz 1.5](../../lecture_notes/de/01_matrix_algebras.md#r-1-5), [Satz 1.6](../../lecture_notes/de/01_matrix_algebras.md#r-1-6).

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
ist keine Überprüfung. Die Überprüfung von $\mathrm{spec}\,f(A)=f(\mathrm{spec}\,A)$ weiter
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

<a id="s-1-6"></a>

## 1.6 Normen und die Dualität der Spurnorm

*Vorlesungsnotizen:* [Abschnitt 1.4](../../lecture_notes/de/01_matrix_algebras.md#s-1-4): [Satz 1.8](../../lecture_notes/de/01_matrix_algebras.md#r-1-8), die drei Normen, [Satz 1.9](../../lecture_notes/de/01_matrix_algebras.md#r-1-9).

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

<a id="s-1-7"></a>

## 1.7 Stapel: viele Operatoren auf einmal

*Vorlesungsnotizen:* [Abschnitt 1.7](../../lecture_notes/de/01_matrix_algebras.md#s-1-7), [Proposition 1.10](../../lecture_notes/de/01_matrix_algebras.md#r-1-10).

Die Stapelachse ist der Grund, weshalb es die Bibliothek gibt: Zufallsmatrixtheorie,
Monte-Carlo-Verfahren und Optimierung mit mehreren Startpunkten benötigen Tausende von
Operatoren, und ein einziger gestapelter Aufruf ersetzt eine Python-Schleife über sie.
Als Beispiel überprüfen wir ein exaktes Moment des **Gaußschen unitären Ensembles** (GUE).

Nach [Vorlesungsnotizen, Proposition 1.10](../../lecture_notes/de/01_matrix_algebras.md#r-1-10) gilt: Hat $G$ unabhängige Einträge mit unabhängigen
Real- und Imaginärteilen $\sim N(0,\tfrac12)$ und ist $X=(G+G^*)/\sqrt{2n}$, so ist $\mathbb E\,\tfrac1n\mathrm{Tr}\,X^4=2+\tfrac1{n^2}$.

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
$\tfrac1n\mathrm{Tr}\,X^4$ ist ein Polynom in Gaußschen Zufallsvariablen und hat daher endliche
Varianz $\sigma^2$ (Lemma A.23: alle Gaußschen Momente sind endlich; Produkte behandelt man mit
$|y_1\cdots y_r|\le\sum_i|y_i|^r$), und nach der Tschebyschow-Ungleichung weicht der Mittelwert von $B$
unabhängigen Stichproben mit Wahrscheinlichkeit höchstens $1/25$ um mehr als $5\sigma/\sqrt B$ vom
exakten Wert ab. Der Code verwendet die Stichprobenstandardabweichung anstelle von $\sigma$, und der
feste Seed macht die Überprüfung reproduzierbar. Nach diesem Muster verläuft in diesem Kurs jede
Überprüfung an bewiesenen Ergebnissen: eine exakte Aussage, ein numerischer Schätzwert und eine
explizite Toleranz.

<a id="s-1-8"></a>

## 1.8 Womit man vergleichen sollte

Für einzelne Matrizen sind `numpy.linalg` und `scipy.linalg` (`sqrtm`, `expm`, `logm`,
`polar`) die Referenzimplementierungen, und die obigen Überprüfungen lassen sich mit ihnen
wiederholen. Die Bibliothek fügt drei Dinge hinzu, die später wichtig werden: Stapelverarbeitung
(Abschnitt 1.7), Operatoren, die eine Blockstruktur tragen (Kapitel 2), und automatisches
Differenzieren (PyTorchs Autograd) durch ihre Operationen hindurch.

<a id="x-1"></a>

## Aufgaben

[Aufgaben 1 und 3](../../lecture_notes/de/01_matrix_algebras.md#x-1) der Vorlesungsnotizen, numerisch:
- Man überprüfe $(fg)(A)=f(A)g(A)$ für $f=\sin$, $g=\exp$ an einem zufälligen selbstadjungierten $A$ (Aufgabe 1).
- Das **geometrische Mittel** von $A,B>0$ ist $A\#B=A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}$, die eindeutige
  positive Lösung $X$ von $XA^{-1}X=B$ (Aufgabe 3). Man überprüfe diese Gleichung und $A\#B=B\#A$ numerisch.

<details>
<summary>Lösungen</summary>

```python
Hs = M3.operator_from_eigenvalues(lambda d: torch.randn(1, d), force_self_adjoint=True)
assert close(Hs.apply_function(lambda w: torch.sin(w) * torch.exp(w)),
             Hs.apply_function(torch.sin) @ Hs.apply_function(torch.exp))
```

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

**Weiter:** [Tutorial 2](02_von_neumann_algebras.md): Kommutanten, das Zentrum und die Algebra $\bigoplus_cM_{k_c}(\mathbb C)$ in der Bibliothek.
