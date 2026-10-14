[English](../README.md) · [Español](../es/README.md) · [Français](../fr/README.md) · **Deutsch** · [中文](../zh/README.md) · [日本語](../ja/README.md) · [Русский](../ru/README.md)

# Endlichdimensionale von-Neumann-Algebren: Vorlesungsnotizen

Dieser Kurs behandelt die Theorie der endlichdimensionalen von-Neumann-Algebren – der
Algebren $\bigoplus_cM_{k_c}(\mathbb C)$ – zusammen mit ihren Anwendungen in der
Quantenphysik, der Quanteninformation, der Theorie der Zufallsmatrizen und dem
maschinellen Lernen.

Er richtet sich zugleich an zwei Arten von Lesern:

- **Mathematiker** finden einen in sich geschlossenen Text: Jeder Satz, jede Formel, jede
  Konstante und jede Ungleichung, die formuliert wird, wird im Kurs aus den unten
  aufgeführten Voraussetzungen bewiesen — im Kapitel selbst oder in einem der Anhänge A–F,
  die mit Nummer zitiert werden. Zahlen, die nur berechnet werden (etwa die Niveaustatistik
  für großes $n$), sind als *numerische Beobachtung* gekennzeichnet und werden nie als
  Prämissen oder Testziele verwendet. Literaturangaben dienen der historischen Zuordnung
  oder als weiterführende Lektüre, nie als Ersatz für einen Beweis.
- **Physiker und Ingenieure** finden die physikalische Bedeutung jedes Begriffs
  in Kästen mit der Überschrift *Physikalische Bedeutung*. Diese Kästen sind physikalische Deutung, keine Sätze des
  Kurses: Nichts, was anderswo bewiesen wird, hängt von ihnen ab. Lauffähigen Code, den sie
  anpassen können, enthält das [Tutorial](../../tutorial/de/README.md).

Die Kapitel 1–13 sind gleich aufgebaut: Definitionen und Ergebnisse mit vollständigen Beweisen, durchgerechnete
Beispiele und Aufgaben mit Lösungen. Nach jedem Ergebnis, das die Bibliothek `torch_vn_algebra` berechnet oder
überprüft, verweist eine Zeile *Rechnung:* oder *In der Bibliothek:* auf die Stelle im [Tutorial](../../tutorial/de/README.md),
an der dies geschieht.

## Inhalt

| | Kapitel | Hauptergebnisse |
|---|---|---|
| 0 | [Worum es in diesem Kurs geht](00_overview.md) | die Algebra als Datentyp: Blockstruktur = gestapelte Tensoren; Aufwand; was es gibt und was neu ist |
| 1 | [Matrixalgebren und Funktionalkalkül](01_matrix_algebras.md) | Spektralsatz, Funktionalkalkül, Polarzerlegung, Dualität der Spurnorm, ein GUE-Moment |
| 2 | [Endlichdimensionale von-Neumann-Algebren](02_von_neumann_algebras.md) | Bikommutantensatz, Struktursatz, Zentrum, Symmetriegruppen und Kommutanten, bedingte Erwartungen |
| 3 | [Spuren, Zustände und Entropie](03_traces_and_states.md) | Klassifikation der Spuren, Dichten, Entropie der Sektoren, Kleinsche Ungleichung, Monotonie unter der bedingten Erwartung auf das Zentrum, Michelson-Kontrast |
| 4 | [Symmetrien und Darstellungen](04_symmetries.md) | Schursches Lemma, Ladungssektoren, $U(1)$ und $SU(2)$, Spinsysteme, reeller, komplexer und quaternionischer Typ |
| 5 | [Zufallsmatrizen in einer Algebra](05_random_matrices.md) | Haar-Maß auf $\prod_cU(k_c)$, Weingarten-Kalkül, Niveaustatistik und warum Sektoren nicht vermischt werden dürfen |
| 6 | [Quantenkanäle](06_quantum_channels.md) | Kraus und Stinespring, duale Abbildungen bezüglich einer Spur, Abbildungen zwischen Sektoren, Dekohärenz |
| 7 | [Tensorprodukte und Verschränkung](07_tensor_products.md) | partielle Spur, Pagesche Formel, Fusion von Ladungen |
| 8 | [Dynamik](08_dynamics.md) | Schrödinger- und Lindblad-Gleichung, Sprünge zwischen Sektoren, Postselektion |
| 9 | [Spinketten](09_spin_chains.md) | Sektoren von Magnetisierung und Impuls, Lanczos, Überprüfungen mit dem Bethe-Ansatz, Unordnung |
| 10 | [Quantencomputing](10_quantum_computing.md) | dekohärenzfreie Unterräume, Fehlerkorrektur mit Operatoralgebren, Mixer mit Nebenbedingungen |
| 11 | [Optimierung und Spurungleichungen](11_optimisation.md) | Parametrisierungen mit Nebenbedingungen, Stichprobenziehung versus Optimierung |
| 12 | [Äquivariantes maschinelles Lernen](12_equivariant_learning.md) | äquivariante Schichten als Kommutanten |
| 13 | [Jordan-Algebren: warum die Quantenmechanik *-Algebren verwendet](13_jordan_algebras.md) | Jordan-Produkt und quadratische Darstellung, Spinfaktoren und Lorentz-Kegel, die Klassifikation von Jordan–von Neumann–Wigner, Transposition und partielle Transposition, eine Spurungleichung in Jordan-Form |
| A | [Tatsachen aus Analysis und Wahrscheinlichkeitstheorie](A_analysis.md) | Potenzreihen und Nullmengen, Gaußsche Vektoren, sphärische Maße, Eigenwertungleichungen, Isserlis–Wick |
| B | [Schur-Weyl-Dualität und der Weingarten-Kalkül](B_schur_weyl.md) | Schur-Weyl-Dualität, Weingarten-Funktionen, Satz von Wigner–Eckart, Haar-Maß auf $O(k)$ |
| C | [Eigenwertdichten und die Pagesche Formel](C_eigenvalue_densities.md) | Weylsche Integrationsformel, GOE/GUE, Wishart- und induzierte Dichten, Pagesche Formel |
| D | [Generatoren quantendynamischer Halbgruppen](D_gksl.md) | der GKSL-Satz, Generatoren auf direkten Summen, Satz von Frigerio, Bény–Kempf–Kribs |
| E | [Klassifikation der euklidischen Jordan-Algebren](E_jordan_classification.md) | Peirce-Zerlegung, Sätze von Hurwitz und Frobenius, Koordinatisierung, die Albert-Algebra, der Satz von Jordan–von Neumann–Wigner |
| F | [Symmetrische Kegel, Ordnungsisomorphismen und Jordan-Abbildungen](F_cones_and_order.md) | Kegel der Quadrate, der Satz von Koecher–Vinberg, der Satz von Kadison, Jordan-Automorphismen von $\bigoplus_cM_{k_c}$ |

<a id="prerequisites"></a>

## Voraussetzungen

Der Kurs setzt nur das Folgende voraus und beweist alles Übrige.

- Lineare Algebra über $\mathbb R$ und $\mathbb C$ (endlichdimensional) und der
  Fundamentalsatz der Algebra.
- Standardstoff der reellen Analysis und der Differentialrechnung in mehreren Variablen:
  Grenzwerte, Stetigkeit, Kompaktheit im $\mathbb R^n$, Differentiation, der
  Mittelwertsatz, die Taylorsche Formel, der Satz über die Umkehrfunktion und der Satz über
  implizite Funktionen.
- Standardstoff der Maß- und Integrationstheorie: Lebesgue-Maß, Produktmaße und
  Fubini–Tonelli, majorisierte und monotone Konvergenz, Variablentransformation im
  $\mathbb R^n$, Bildmaße, Dichten, Unabhängigkeit.
- Elementare Wahrscheinlichkeitstheorie: Erwartungswert, Varianz, Unabhängigkeit, die
  Normalverteilung und ihre momenterzeugende Funktion, die Markow- und die
  Tschebyschow-Ungleichung, das Gesetz der großen Zahlen.

Physik und Funktionalanalysis werden nicht vorausgesetzt.

## Weiterführende Literatur

M. Takesaki, *Theory of Operator Algebras I* (Springer, 2002); V. Jones,
*von Neumann algebras*, Vorlesungsskript; D. Petz, *Quantum Information Theory and Quantum
Statistics* (Springer, 2008); J. Watrous, *The Theory of Quantum Information* (Cambridge,
2018); M. L. Mehta, *Random Matrices* (Elsevier, 2004).

Die Bibliothek `torch_vn_algebra` und ihre Verwendung im Kontext dieser Notizen sind Gegenstand des
[Tutorials](../../tutorial/de/README.md). Die Begleitarbeit zur Bibliothek ist
[`man/article/main.pdf`](https://github.com/andre-a-no/von_neumann_type_i/blob/main/man/article/main.pdf).
