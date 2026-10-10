[English](../00_overview.md) · [Español](../es/00_overview.md) · [Français](../fr/00_overview.md) · **Deutsch** · [中文](../zh/00_overview.md) · [日本語](../ja/00_overview.md) · [Русский](../ru/00_overview.md)

# 0. Worum es in diesem Kurs geht

Von-Neumann-Algebren werden fast ausschließlich mit Papier und Bleistift untersucht. In
endlicher Dimension ist jede von ihnen eine direkte Summe voller Matrixalgebren (Kapitel 2),
sodass sich im Prinzip alles berechnen lässt – in der Praxis wird eine Algebra jedoch meist zu
einer einzigen großen Matrix verflacht, und ihre Struktur (die Sektoren, das Zentrum, die
vielen Spuren) verschwindet aus der Rechnung. Dieser Kurs und die Bibliothek
`torch_vn_algebra` gehen den umgekehrten Weg: **Die Algebra selbst ist der Datentyp**, und
Rechnungen in ihr werden effizient ausgeführt.

## 0.1 Die Beobachtung

Eine endlichdimensionale von-Neumann-Algebra ist bis auf einen unitären Basiswechsel
$$
\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c}
$$
(Satz 2.6). Als abstrakte Algebra ist sie $\bigoplus_cM_{k_c}(\mathbb C)$: Ein Element ist eine
Liste von $C$ quadratischen Blöcken. Das ist genau die Gestalt eines Tensors mit einer
*Sektorachse*: Ein Element von $\mathcal M$ ist ein Array der Gestalt $(C,k,k)$ (die Blöcke
aufgefüllt auf $k=\max_ck_c$), und $B$ Elemente – Monte-Carlo-Stichproben, Startpunkte einer
Optimierung, Zeitschritte, Realisierungen der Unordnung – bilden ein Array der Gestalt
$(B,C,k,k)$. Moderne Array-Bibliotheken (PyTorch auf CPUs und GPUs) sind gerade für solche
„Stapel kleiner Matrizen“ gebaut, weil neuronale Netze sie benötigen.

Unter dieser Identifikation sind die Operationen der Algebra gestapelte Array-Operationen,
exakt und ohne Näherung:

| in der Algebra $\mathcal M$ | auf dem Array der Gestalt $(B,C,k,k)$ |
|---|---|
| Produkt $AB$, Summe, Adjungierte | gestapeltes Matrixprodukt, Summe, konjugierte Transposition über $(B,C)$ |
| Funktionalkalkül $f(A)$ (Kapitel 1) | gestapelte Eigenzerlegung der Blöcke |
| eine Spur $\sum_cw_c\operatorname{Tr}A_c$ (Kapitel 3) | Diagonalsumme pro Block, gewichtete Summe über die Sektorachse |
| das Zentrum, zentrale Projektionen (Kapitel 2) | eine Zahl pro Sektor |
| bedingte Erwartung auf das Zentrum | Blockspur geteilt durch $k_c$, zurück auf die Blöcke verteilt |
| Haar-Maß auf der unitären Gruppe $\prod_cU(k_c)$ | gestapelte QR-Zerlegung Gaußscher Blöcke (Kapitel 5) |
| Zustände, Kanäle, Dynamik (Kapitel 3, 6, 8) | dieselben Operationen, hintereinander ausgeführt |

Zwei Konsequenzen machen dies zu mehr als bloßer Buchführung.

- **Aufwand.** Ein Produkt oder eine Spektralzerlegung in $\mathcal M$ kostet $\sum_ck_c^3$
  Operationen, gegenüber $D^3$ für die Einbettung in $M_D(\mathbb C)$, $D=\sum_ck_c$; für $C$
  gleich große Blöcke ist das ein Faktor $C^2$. Produkte, Summen und Spuren verarbeiten alle
  $B\cdot C$ Blöcke in einem einzigen Aufruf; Zerlegungen (Eigenwerte, Singulärwerte) laufen als ein
  Aufruf pro vorkommender Blockgröße – ein einziger Aufruf, wenn alle Sektoren gleich groß sind.
  Das lastet eine GPU aus.
- **Korrektheit durch Konstruktion.** Alles, was berechnet wird, bleibt in $\mathcal M$:
  Zufällige unitäre Elemente sind Haar-verteilt auf $\prod_cU(k_c)$ (einer Nullmenge in $U(D)$,
  die man durch Stichproben aus $U(D)$ nie erreicht), Gradientenschritte verlassen die Algebra
  nicht, und die Spur, die Dichten und die bedingten Erwartungen beziehen sich auf die Algebra
  und nicht auf eine zufällige Einbettung.

## 0.2 Eine erste Rechnung

Wir vergleichen dieselbe Rechnung in $\mathcal M=\bigoplus_{c=1}^{64}M_{16}(\mathbb C)$ und in
dessen Einbettung in $M_{1024}(\mathbb C)$: die Quadratwurzel eines zufälligen positiven Elements
und deren Spur.

```python
import time
import torch
from torch_vn_algebra import TypeIAlgebra

torch.manual_seed(0)
C, k = 64, 16
alg = TypeIAlgebra([k] * C, [k] * C, precision='double')       # M = (+)_{c=1}^{64} M_16(C)

A = alg.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(1, d), force_positive=True)

t0 = time.perf_counter()
root_blocks = A.sqrt().matrix                                    # (1, 64, 16, 16): 64 small eigenproblems
t_blocks = time.perf_counter() - t0

big = torch.block_diag(*A.matrix[0])                             # the same element in M_1024(C)
t0 = time.perf_counter()
w, V = torch.linalg.eigh(big)                                    # one large eigenproblem
root_big = (V * w.sqrt()) @ V.conj().T
t_big = time.perf_counter() - t0

assert torch.allclose(torch.block_diag(*root_blocks[0]), root_big, atol=1e-10)   # same result
assert abs(A.sqrt().Tr_blunt().item() - root_big.trace().item()) < 1e-9
print(f"in the algebra: {t_blocks * 1e3:.1f} ms, as one 1024 x 1024 matrix: {t_big * 1e3:.1f} ms")
```

Die Ergebnisse stimmen bis auf Rundungsfehler überein. Die Version in der Algebra ist schneller –
bei dieser Größe auf einer typischen CPU etwa um das Zehnfache (etwa 16 ms gegenüber 200 ms auf der
4-Kern-Maschine, auf der dieser Kurs getestet wurde; alle 64 Blöcke gehen in einem einzigen Aufruf
an den Eigenwertlöser) –, allerdings weit weniger, als die Zahl der
Operationen $C^2=4096$ nahelegt, weil ein einziger Aufruf eines großen Eigenwertlösers die
Hardware besser nutzt als viele kleine. Der Abstand wächst mit der Zahl der Sektoren, und die
Speicherersparnis beträgt einen Faktor $C$. Bei einem Stapel von $B$ Elementen bleibt die Version in der Algebra ein einziger
Aufruf; die Version mit der großen Matrix wird zu $B$ großen Eigenwertproblemen.

## 0.3 Was es gibt und was neu ist

Wir haben nach Software gesucht, die numerisch *innerhalb* endlichdimensionaler
von-Neumann-Algebren rechnet. Was wir gefunden haben, fällt in vier Gruppen, alle nützlich und
alle mit unterschiedlichem Zweck:

- **Allgemeine Matrixsoftware** (NumPy/SciPy, QuTiP und andere Quanten-Toolboxen) arbeitet mit
  jeweils einer Matrix. Eine etwa vorhandene Blockstruktur muss von Hand verwaltet werden.
- **Bibliotheken für symmetrische Tensornetzwerke** speichern Tensoren als blockdünnbesetzte
  Arrays über Ladungssektoren einer abelschen Symmetrie, einige mit einem PyTorch-Backend für
  GPUs und automatisches Differenzieren (z. B. YASTN [1]; ferner TeNPy, ITensor, TensorKit).
  Technisch sind sie die nächsten Verwandten. Ihr Gegenstand ist jedoch ein zu kontrahierendes
  Tensornetzwerk, keine Operatoralgebra: Sie stellen weder die Familie der Spuren, Zustände und
  Dichten bereit noch das Zentrum und bedingte Erwartungen, Haar-Maße auf $\prod_cU(k_c)$ oder
  Kanäle zwischen Sektoren.
- **Numerische Blockdiagonalisierung** von Matrix-$*$-Algebren [2, 3, 4] findet die Zerlegung
  $\mathcal M\cong\bigoplus_cM_{k_c}\otimes1_{m_c}$ einer gegebenen Algebra – typischerweise, um
  semidefinite Programme mit Symmetrie zu verkleinern. Das ist komplementär: Jene Algorithmen
  finden die Blöcke; diese Bibliothek rechnet in ihnen. Exakte symbolische Zerlegungen von
  Gruppenalgebren gibt es in der Computeralgebra (z. B. das GAP-Paket Wedderga).
- **Nichtkommutative polynomiale Optimierung** (NCSOStools, ncpol2sdpa) schätzt Eigenwerte und
  Spuren von Polynomen in nichtkommutierenden Variablen über alle Matrixgrößen zugleich durch
  semidefinite Relaxierungen ab [5]. Damit wird eine andere
  Frage (eine für alle Algebren gültige Schranke) mit anderen Mitteln (konvexe Relaxierung)
  beantwortet.

Unseres Wissens gibt es keine allgemeine Bibliothek, in der eine endlichdimensionale
von-Neumann-Algebra der grundlegende Datentyp für gestapelte, differenzierbare numerische
Rechnungen ist. Das Neue ist daher keine neue Mathematik, sondern **effizientes Rechnen**: die
Identifikation aus §0.1, so implementiert, dass die Struktur der Algebra erhalten bleibt und
jede Operation als gestapelte Kernel läuft (einer pro vorkommender Blockgröße bei Zerlegungen). Wenn Sie verwandte Software kennen,
freuen wir uns über einen Hinweis.

## 0.4 Wie man diesen Kurs liest

| Kapitel | Inhalt |
|---|---|
| 1–3 | die Mathematik: Matrixalgebren, der Struktursatz, Spuren, Zustände und Entropie – mit vollständigen Beweisen |
| 4–8 | woher diese Algebren kommen und was man in ihnen berechnet: Symmetrien, Zufallsmatrizen, Kanäle, Tensorprodukte, Dynamik |
| 9–12 | Anwendungen: Spinketten, Quantencomputing, Optimierung, äquivariantes maschinelles Lernen |

Jedes Kapitel formuliert Ergebnisse mit Beweisen, zeigt den Bibliothekscode, überprüft den Code
an einem exakten Ergebnis und nennt die Standardwerkzeuge zum Vergleich.

## Literatur

1. M. M. Rams, G. Wójtowicz, A. Sinha, J. Hasik, *YASTN: Yet another symmetric tensor networks; A
   Python library for Abelian symmetric tensor network calculations*, SciPost Phys. Codebases 52
   (2025).
2. K. Murota, Y. Kanno, M. Kojima, S. Kojima, *A numerical algorithm for block-diagonal
   decomposition of matrix $*$-algebras*, Technical Report METR 2007-52, University of Tokyo (2007).
3. T. Maehara, K. Murota, *A numerical algorithm for block-diagonal decomposition of matrix
   $*$-algebras, Part II: general algorithm*, Optimization Online (2009).
4. E. de Klerk, C. Dobre, D. V. Pasechnik, *Numerical block diagonalization of matrix
   $*$-algebras with application to semidefinite programming*, Math. Program. 129 (2011) 91–111.
5. I. Klep, J. Povh, *Semidefinite programming and sums of hermitian squares of noncommutative
   polynomials*, J. Pure Appl. Algebra 214 (2010) 740–749.

**Weiter:** [Kapitel 1](01_matrix_algebras.md).
