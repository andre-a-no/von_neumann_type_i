[English](../00_overview.md) · [Español](../es/00_overview.md) · [Français](../fr/00_overview.md) · **Deutsch** · [中文](../zh/00_overview.md) · [日本語](../ja/00_overview.md) · [Русский](../ru/00_overview.md)

# Tutorial 0. Worum es in diesem Kurs geht

Dieses Tutorial begleitet [Kapitel 0](../../lecture_notes/de/00_overview.md) der Vorlesungsnotizen, das ein Element einer
von-Neumann-Algebra $`\bigoplus_cM_{k_c}(\mathbb C)`$ mit einem Array der Gestalt $`(C,k,k)`$ identifiziert und beweist,
dass ihre Operationen gestapelte Array-Operationen sind. Hier sehen wir, was das in der Bibliothek
`torch_vn_algebra` ergibt: was es kostet, eine erste Rechnung und wie sich die Bibliothek zu vorhandener Software verhält.

Von-Neumann-Algebren werden meist mit Papier und Bleistift untersucht. In endlicher Dimension
ist jede von ihnen $`*`$-isomorph zu einer direkten Summe voller Matrixalgebren (Satz 2.6), sodass
sich im Prinzip alles berechnen lässt – wird eine solche Algebra jedoch numerisch behandelt, so
wird sie oft zu einer einzigen großen Matrix verflacht, und ihre Struktur (die Sektoren, das
Zentrum, die vielen Spuren) verschwindet aus der Rechnung. Dieser Kurs und die Bibliothek
`torch_vn_algebra` gehen den umgekehrten Weg: **Die Algebra selbst ist der Datentyp**, und
Rechnungen in ihr werden effizient ausgeführt.

<a id="s-0-1"></a>

## 0.1 Die Beobachtung

*Vorlesungsnotizen:* [Abschnitt 0.1](../../lecture_notes/de/00_overview.md#s-0-1), der aus [Satz 2.6](../../lecture_notes/de/02_von_neumann_algebras.md#r-2-6),
[Satz 3.1](../../lecture_notes/de/03_traces_and_states.md#r-3-1), [Satz 5.3](../../lecture_notes/de/05_random_matrices.md#r-5-3) und
[Korollar 5.4](../../lecture_notes/de/05_random_matrices.md#r-5-4) beweist, dass die ersten sechs Zeilen der folgenden Tabelle Identitäten sind:

| in der Algebra $`\mathcal M`$ | auf dem Array der Gestalt $`(B,C,k,k)`$ |
|---|---|
| Produkt $`AB`$, Summe, Adjungierte | gestapeltes Matrixprodukt, Summe, konjugierte Transposition über $`(B,C)`$ |
| Funktionalkalkül $`f(A)`$ (Kapitel 1) | gestapelte Eigenzerlegung der Blöcke |
| eine Spur $`\sum_cw_c\mathrm{Tr}\,A_c`$ (Kapitel 3) | Diagonalsumme pro Block, gewichtete Summe über die Sektorachse |
| das Zentrum, zentrale Projektionen (Kapitel 2) | eine Zahl pro Sektor |
| bedingte Erwartung auf das Zentrum | Blockspur geteilt durch $`k_c`$, zurück auf die Blöcke verteilt |
| Haar-Maß auf der unitären Gruppe $`\prod_cU(k_c)`$ | gestapelte QR-Zerlegung Gaußscher Blöcke, mit der Phasenkorrektur aus Satz 5.3 |
| Zustände, Kanäle, Dynamik (Kapitel 3, 6, 8) | aus den obigen Operationen aufgebaut |

Zwei Konsequenzen machen dies zu mehr als bloßer Buchführung.

- **Aufwand.** Eine Spektralzerlegung in $`\mathcal M`$ wird blockweise berechnet (Zeile 2 der
  Tabelle). Die üblichen dichten Eigenwertlöser (LAPACK, das PyTorch auf einer CPU aufruft)
  benötigen eine Anzahl von Gleitkommaoperationen, die wie die dritte Potenz der Matrixgröße
  wächst; das ist eine Beschreibung dieser Algorithmen, kein Ergebnis dieses Kurses. In dem
  Modell, in dem ein $`k\times k`$-Eigenwertproblem eine Konstante mal $`k^3`$ kostet, kostet eine
  Spektralzerlegung in $`\mathcal M`$ $`\sum_ck_c^3`$ gegenüber $`D^3`$ für die Einbettung in
  $`M_D(\mathbb C)`$, $`D=\sum_ck_c`$; für $`C`$ gleich große Blöcke der Größe $`k`$ ist das Verhältnis
  $`(Ck)^3/(Ck^3)=C^2`$. Produkte laufen auf dem aufgefüllten Tensor; das Schulbuchprodukt zweier
  $`k\times k`$-Matrizen braucht $`k^3`$ Multiplikationen ($`k^2`$ Einträge, jeder eine Summe von $`k`$
  Produkten), also braucht ein Produkt in $`\mathcal M`$ $`C\,k_{\max}^3`$, was mit $`\sum_ck_c^3`$
  übereinstimmt, wenn alle Blöcke gleich groß sind. Ein Element wird als $`C\,k_{\max}^2`$ Zahlen
  gespeichert gegenüber $`D^2`$; für gleich große Blöcke ist das Verhältnis $`(Ck)^2/(Ck^2)=C`$. Was
  der Code tut: Produkte, Summen und Spuren sind einzelne Tensoroperationen auf dem ganzen Array
  $`(B,C,k_{\max},k_{\max})`$; hermitesche Eigenzerlegungen, Singulärwerte und QR-Zerlegungen werden
  auf den Stapel aller $`B\cdot n_k`$ Blöcke jeder vorkommenden Größe $`k`$ ($`n_k`$ die Zahl der Sektoren
  dieser Größe) in einem gestapelten Aufruf angewandt – ein einziger Aufruf, wenn alle Sektoren
  gleich groß sind. Auf einer CPU mit $`t\ge2`$ PyTorch-Threads teilt die Bibliothek einen solchen
  Aufruf in $`t`$ Teilstücke auf, die in parallelen Threads laufen, sofern der Stapel mindestens
  $`2t`$ Blöcke hat und (Zahl der Blöcke)$`\,\cdot k^3\ge2^{18}`$ gilt (Modul `cost`, Funktion `_run`); die
  Ergebnisse sind dieselben LAPACK-Aufrufe auf denselben Matrizen. Ziel ist, eine GPU oder alle
  Kerne einer CPU auszulasten.
- **Korrektheit durch Konstruktion.** Jedes mit Nullen aufgefüllte Array der Gestalt
  $`(C,k_{\max},k_{\max})`$ ist ein Element von $`\mathcal M`$; was aus den Blöcken berechnet wird,
  bleibt also in $`\mathcal M`$. Die zufälligen unitären Elemente der Bibliothek werden mit der
  Konstruktion aus Satz 5.3 unabhängig in jedem Block gezogen (`random_unitary_operator` ruft den
  blockweisen Sampler `random_unitary` einmal pro Blockgröße auf, und Abschnitt 5.1 prüft
  `random_unitary` an dieser Konstruktion), was nach Korollar 5.4 das Haar-Maß von
  $`\prod_cU(k_c)`$ ergibt (für ideale Gaußsche Zufallszahlen); ein Haar-verteiltes unitäres Element
  von $`U(D)`$ liegt für $`C\ge2`$ mit Wahrscheinlichkeit null in dieser Untergruppe (Abschnitt 5.1,
  nach Korollar 5.4). Ein Gradientenschritt verändert die Blöcke und kann die Algebra daher nicht
  verlassen. Und die Spuren sind die Funktionale $`\sum_cw_c\mathrm{Tr}\,A_c`$ aus Satz 3.1 mit
  explizit gewählten Gewichten (und die Dichten und bedingten Erwartungen werden bezüglich dieser
  gebildet), nicht die Spur irgendeines $`\mathbb C^n`$, auf dem $`\mathcal M`$ zufällig wirkt, die
  $`w_c=m_c`$ festlegen würde (Abschnitt 2.3).

<a id="s-0-2"></a>

## 0.2 Eine erste Rechnung

*Vorlesungsnotizen:* [Abschnitt 0.1](../../lecture_notes/de/00_overview.md#s-0-1), Zeile 2 der Tabelle: $`f(A)=\bigoplus_cf(A_c)`$, eine Eigenzerlegung pro Block.

Wir vergleichen dieselbe Rechnung in $`\mathcal M=\bigoplus_{c=1}^{64}M_{16}(\mathbb C)`$ und in
dessen Einbettung in $`M_{1024}(\mathbb C)`$: die Quadratwurzel eines zufälligen positiven Elements
und ihre Spur. `TypeIAlgebra(n_factors, k_factors)` nimmt zwei Listen von Blockgrößen: die
nominellen Größen $`n_c`$ und die Größen $`k_c\le n_c`$ der tatsächlich verwendeten Blöcke
(Abschnitt 1.5). Jede Algebra, die der Kurs mit `TypeIAlgebra` konstruiert, hat zwei gleiche
Listen.

```python
import time
import torch
from torch_vn_algebra import TypeIAlgebra

torch.manual_seed(0)
C, k = 64, 16
alg = TypeIAlgebra([k] * C, [k] * C, precision='double')       # M = (+)_{c=1}^{64} M_16(C)

A = alg.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(1, d), force_positive=True)
blocks = A.matrix                                                # build A before timing

t0 = time.perf_counter()
root_blocks = A.sqrt().matrix                                    # (1, 64, 16, 16): 64 small eigenproblems
t_blocks = time.perf_counter() - t0

big = torch.block_diag(*blocks[0])                               # the same element in M_1024(C)
t0 = time.perf_counter()
w, V = torch.linalg.eigh(big)                                    # one large eigenproblem
root_big = (V * w.sqrt()) @ V.conj().T
t_big = time.perf_counter() - t0

assert torch.allclose(torch.block_diag(*root_blocks[0]), root_big, atol=1e-10)   # row 2: f(A) = (+)_c f(A_c)
assert abs(A.sqrt().Tr_blunt().item() - root_big.trace().item()) < 1e-9
assert big.numel() == C * blocks.numel()                          # storage: a factor C = 64
print(f"in the algebra: {t_blocks * 1e3:.1f} ms, as one 1024 x 1024 matrix: {t_big * 1e3:.1f} ms")
```

Die Ergebnisse stimmen bis auf Rundungsfehler überein (Zeile 2 der Tabelle), und das
Speicherverhältnis ist $`C`$, wie oben gezählt. Die 64 Blöcke gehen als ein gestapelter Aufruf auf
einem Stapel der Gestalt $`(64,16,16)`$ an den Eigenwertlöser; mit $`2\le t\le32`$ PyTorch-Threads
teilt die Bibliothek ihn in $`t`$ Teilstücke auf ($`64\cdot16^3=2^{18}`$ erreicht die Schwelle des
Moduls `cost` genau), mit einem Thread nicht. *Numerische Beobachtung* (die Laufzeiten hängen
von der Hardware ab; der Code gibt sie aus): Auf zwei 4-Kern-CPUs brauchte die Version in der
Algebra 3–10 ms und die Version mit der großen Matrix 400–760 ms, ein Faktor zwischen einigen
Zehn und einigen Hundert (auf der zweiten Maschine mit 4 Threads und auch mit einem). Der
gemessene Faktor liegt weit unter dem Verhältnis $`C^2=4096`$ des kubischen Kostenmodells: Eine
Operationszählung ist ein Modell der Laufzeit, nicht die Laufzeit, und auf diesen Maschinen
liefen die kleinen Eigenwertprobleme mit einer viel geringeren Rate an Operationen pro Sekunde
als das große. Im Kostenmodell wächst das Verhältnis $`C^2`$ mit der Zahl der Sektoren. Bei einem
Stapel von $`B`$ Elementen bleibt die Version in der Algebra ein einziger gestapelter Aufruf; die
Version mit der großen Matrix muss $`B`$ Eigenwertprobleme der Größe $`1024`$ lösen.

<a id="s-0-3"></a>

## 0.3 Was es gibt und was neu ist

Wir haben nach Software gesucht, die numerisch *innerhalb* endlichdimensionaler
von-Neumann-Algebren rechnet. Was wir gefunden haben, fällt in vier Gruppen, alle nützlich und
alle mit unterschiedlichem Zweck:

- **Allgemeine Matrixsoftware** (NumPy/SciPy, QuTiP und andere Quanten-Toolboxen) arbeitet mit
  jeweils einer Matrix. Eine etwa vorhandene Blockstruktur muss von Hand verwaltet werden.
- **Bibliotheken für symmetrische Tensornetzwerke** speichern Tensoren als blockdünnbesetzte
  Arrays über den Ladungssektoren einer Symmetrie (abelsch in YASTN [1], TeNPy und ITensor;
  TensorKit behandelt auch nichtabelsche); YASTN hat ein PyTorch-Backend für GPUs und
  automatisches Differenzieren. Technisch sind sie die nächsten Verwandten. Ihr Gegenstand ist
  jedoch ein zu kontrahierendes Tensornetzwerk, keine Operatoralgebra: Soweit wir ihrer
  Dokumentation entnehmen konnten, stellen sie weder die Familie der Spuren, Zustände und Dichten
  bereit noch das Zentrum und bedingte Erwartungen, Haar-Maße auf $`\prod_cU(k_c)`$ oder Kanäle
  zwischen Sektoren.
- **Numerische Blockdiagonalisierung** von Matrix-$`*`$-Algebren [2, 3, 4] findet die Zerlegung
  $`\mathcal M\cong\bigoplus_cM_{k_c}\otimes1_{m_c}`$ einer gegebenen Algebra – typischerweise, um
  semidefinite Programme mit Symmetrie zu verkleinern. Das ist komplementär: Jene Algorithmen
  finden die Blöcke; diese Bibliothek rechnet in ihnen. Ihr Eingabeschritt (Modul `decompose`,
  Abschnitt 4.3) verwendet die Idee aus [2, 3], dass die Eigenräume eines zufälligen Elements der
  Kommutante die Blöcke finden; Abschnitt 4.3 beweist, dass dies in exakter Arithmetik mit
  Wahrscheinlichkeit eins gelingt (Lemma 4.6, Proposition 4.7, Lemma 4.8), und Proposition 4.9
  ist das Kriterium, das `Decomposition.verify()` numerisch prüft. Exakte symbolische Zerlegungen
  von Gruppenalgebren gibt es in der Computeralgebra (z. B. das GAP-Paket Wedderga).
- **Nichtkommutative polynomiale Optimierung** (NCSOStools, ncpol2sdpa) schätzt Eigenwerte und
  Spuren von Polynomen in nichtkommutierenden Variablen über alle Matrixgrößen zugleich durch
  semidefinite Relaxierungen ab [5]. Damit wird eine andere Frage (eine für alle Matrixgrößen
  gültige Schranke) mit anderen Mitteln (konvexe Relaxierung) beantwortet.

*Vorlesungsnotizen:* der randomisierte Zerlegungsalgorithmus und sein Kriterium — [Abschnitt 4.3](../../lecture_notes/de/04_symmetries.md#s-4-3), [Proposition 4.7](../../lecture_notes/de/04_symmetries.md#r-4-7), [Proposition 4.9](../../lecture_notes/de/04_symmetries.md#r-4-9).

Diese Übersicht beschreibt die Software, wie wir sie vorgefunden haben; sie ist kein Satz.
Unseres Wissens gibt es keine allgemeine Bibliothek, in der eine endlichdimensionale
von-Neumann-Algebra der grundlegende Datentyp für gestapelte, differenzierbare numerische
Rechnungen ist. Das Neue ist daher keine neue Mathematik, sondern **effizientes Rechnen**: die
Identifikation aus §0.1, so implementiert, dass die Struktur der Algebra erhalten bleibt und die
Operationen als gestapelte Kernel laufen (ein Aufruf pro vorkommender Blockgröße bei
Zerlegungen, auf einer Mehrkern-CPU oberhalb der Größenschwelle aus §0.1 in ein Teilstück pro
Thread aufgeteilt). Wenn Sie verwandte Software kennen, freuen wir uns über einen Hinweis.

<a id="s-0-4"></a>

## 0.4 Wie man diesen Kurs liest

| Kapitel | Inhalt |
|---|---|
| 1–3 | die Mathematik: Matrixalgebren, der Struktursatz, Spuren, Zustände und Entropie – mit vollständigen Beweisen |
| 4–8 | woher diese Algebren kommen und was man in ihnen berechnet: Symmetrien, Zufallsmatrizen, Kanäle, Tensorprodukte, Dynamik |
| 9–12 | Anwendungen: Spinketten, Quantencomputing, Optimierung, äquivariantes maschinelles Lernen |
| 13 | Grundlagen: Jordan-Algebren und warum die Quantenmechanik *-Algebren verwendet |
| Anhänge | vollständige Beweise der Hilfsresultate, die die Kapitel mit Nummer zitieren (z. B. Satz A.24) |

Jedes Tutorial begleitet das Kapitel der [Vorlesungsnotizen](../../lecture_notes/de/00_overview.md) mit derselben Nummer: Es zeigt
den Bibliothekscode, überprüft den Code numerisch an dort bewiesenen Ergebnissen und nennt die Standardwerkzeuge
zum Vergleich. Jeder Abschnitt beginnt mit einer Zeile *Vorlesungsnotizen:*, die auf den Abschnitt und die von ihm
verwendeten Ergebnisse verweist. Die Codeblöcke eines Tutorials sind dafür gedacht, der Reihe nach in einer
einzigen Python-Sitzung ausgeführt zu werden. Zahlen, die die Vorlesungsnotizen nicht beweisen (etwa Laufzeiten),
erscheinen nur als *numerische Beobachtungen* zum eigenen Code des Tutorials.

## Literatur

1. M. M. Rams, G. Wójtowicz, A. Sinha, J. Hasik, *YASTN: Yet another symmetric tensor networks; A
   Python library for Abelian symmetric tensor network calculations*, SciPost Phys. Codebases 52
   (2025).
2. K. Murota, Y. Kanno, M. Kojima, S. Kojima, *A numerical algorithm for block-diagonal
   decomposition of matrix $`*`$-algebras with application to semidefinite programming*, Japan J.
   Indust. Appl. Math. 27 (2010) 125–160 (Preprint: Technical Report METR 2007-52, University of
   Tokyo, 2007).
3. T. Maehara, K. Murota, *A numerical algorithm for block-diagonal decomposition of matrix
   $`*`$-algebras with general irreducible components*, Japan J. Indust. Appl. Math. 27 (2010)
   263–293.
4. E. de Klerk, C. Dobre, D. V. Pasechnik, *Numerical block diagonalization of matrix
   $`*`$-algebras with application to semidefinite programming*, Math. Program. 129 (2011) 91–111.
5. I. Klep, J. Povh, *Semidefinite programming and sums of hermitian squares of noncommutative
   polynomials*, J. Pure Appl. Algebra 214 (2010) 740–749.

**Weiter:** [Tutorial 1](01_matrix_algebras.md): Operatoren, Funktionen von Operatoren und Stapel in der Bibliothek.
