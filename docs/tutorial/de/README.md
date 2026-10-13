[English](../README.md) · [Español](../es/README.md) · [Français](../fr/README.md) · **Deutsch** · [中文](../zh/README.md) · [日本語](../ja/README.md) · [Русский](../ru/README.md)

# Endlichdimensionale von-Neumann-Algebren: ein Kurs mit `torch_vn_algebra`

Dieser Kurs behandelt die Theorie der endlichdimensionalen von-Neumann-Algebren – der
Algebren $\bigoplus_cM_{k_c}(\mathbb C)$ – zusammen mit ihren Anwendungen in der
Quantenphysik, der Quanteninformation, der Theorie der Zufallsmatrizen und dem
maschinellen Lernen, und zeigt bei jedem Schritt, wie man mit ihnen in `torch_vn_algebra`
rechnet.

Er richtet sich zugleich an zwei Arten von Lesern:

- **Mathematiker** finden vollständige Beweise aller verwendeten
  Ergebnisse. Ohne Beweis zitiert werden nur Standardtatsachen der linearen Algebra, der
  Fundamentalsatz der Algebra und elementare Tatsachen aus Analysis und
  Wahrscheinlichkeitstheorie (Grenzwerte, Konvexität, die momenterzeugende Funktion der
  Normalverteilung); wird ein Ergebnis formuliert, bevor es bewiesen ist, oder nur
  angewandt, so wird dies ausdrücklich gesagt.
- **Physiker und Ingenieure** finden die physikalische Bedeutung jedes Begriffs
  in Kästen mit der Überschrift *Physikalische Bedeutung* sowie lauffähigen Code,
  den sie anpassen können.

Die Kapitel 1–13 sind gleich aufgebaut: Theorie mit Beweisen; der Bibliothekscode zu jedem
Ergebnis; eine Überprüfung des Codes an einem exakten Ergebnis (einem Satz, einer
geschlossenen Formel oder einer unabhängigen Rechnung); ein Abschnitt *Womit man
vergleichen sollte*, der die Standardwerkzeuge für dieselbe Aufgabe nennt; und Aufgaben
mit Lösungen.

Der gesamte Code des Kurses wird von der Testsuite (`tests/test_tutorial.py`) ausgeführt;
jedes Beispiel funktioniert also mit der aktuellen Version der Bibliothek.

## Inhalt

| | Kapitel | Hauptergebnisse |
|---|---|---|
| 0 | [Worum es in diesem Kurs geht](00_overview.md) | die Algebra als Datentyp: Blockstruktur = gestapelte Tensoren; Aufwand; was es gibt und was neu ist |
| 1 | [Matrixalgebren und Funktionalkalkül](01_matrix_algebras.md) | Spektralsatz, Funktionalkalkül, Polarzerlegung, Dualität der Spurnorm, ein GUE-Moment |
| 2 | [Endlichdimensionale von-Neumann-Algebren](02_von_neumann_algebras.md) | Bikommutantensatz, Struktursatz, Zentrum, Symmetriegruppen und Kommutanten, bedingte Erwartungen |
| 3 | [Spuren, Zustände und Entropie](03_traces_and_states.md) | Klassifikation der Spuren, Dichten, Entropie der Sektoren, Kleinsche Ungleichung, Monotonie unter der bedingten Erwartung auf das Zentrum, Michelson-Kontrast |
| 4 | [Symmetrien und Darstellungen](04_symmetries.md) | Schursches Lemma, Ladungssektoren, $U(1)$ und $SU(2)$, Spinsysteme |
| 5 | [Zufallsmatrizen in einer Algebra](05_random_matrices.md) | Haar-Maß auf $\prod_cU(k_c)$, Weingarten-Kalkül, Niveaustatistik und warum Sektoren nicht vermischt werden dürfen |
| 6 | [Quantenkanäle](06_quantum_channels.md) | Kraus und Stinespring, duale Abbildungen bezüglich einer Spur, Abbildungen zwischen Sektoren, Dekohärenz |
| 7 | [Tensorprodukte und Verschränkung](07_tensor_products.md) | partielle Spur, Pagesche Formel, Fusion von Ladungen |
| 8 | [Dynamik](08_dynamics.md) | Schrödinger- und Lindblad-Gleichung, Sprünge zwischen Sektoren, Postselektion |
| 9 | [Spinketten](09_spin_chains.md) | Sektoren von Magnetisierung und Impuls, Lanczos, Überprüfungen mit dem Bethe-Ansatz, Unordnung |
| 10 | [Quantencomputing](10_quantum_computing.md) | dekohärenzfreie Unterräume, Fehlerkorrektur mit Operatoralgebren, Mixer mit Nebenbedingungen |
| 11 | [Optimierung und Spurungleichungen](11_optimisation.md) | Parametrisierungen mit Nebenbedingungen, Stichprobenziehung versus Optimierung |
| 12 | [Äquivariantes maschinelles Lernen](12_equivariant_learning.md) | äquivariante Schichten als Kommutanten |
| 13 | [Jordan-Algebren: warum die Quantenmechanik *-Algebren verwendet](13_jordan_algebras.md) | Jordan-Produkt und quadratische Darstellung, Spinfaktoren und Lorentz-Kegel, die Klassifikation von Jordan–von Neumann–Wigner, Transposition und partielle Transposition, eine Spurungleichung in Jordan-Form |

## Voraussetzungen

Lineare Algebra über $\mathbb C$ (Eigenwerte, Skalarprodukte, unitäre Matrizen),
Grundlagen der Wahrscheinlichkeitstheorie und Python. Physik und Funktionalanalysis
werden nicht vorausgesetzt.

## Wie man den Code ausführt

```bash
pip install -e .                # from the root of the repository
python -m pytest tests/test_tutorial.py -q
```

Alle Codeblöcke eines Kapitels können der Reihe nach in eine einzige Python-Sitzung oder
ein Jupyter-Notebook eingefügt werden. Die Blöcke eines Kapitels teilen sich ihre Variablen.

## Weiterführende Literatur

M. Takesaki, *Theory of Operator Algebras I* (Springer, 2002); V. Jones,
*von Neumann algebras*, Vorlesungsskript; D. Petz, *Quantum Information Theory and Quantum
Statistics* (Springer, 2008); J. Watrous, *The Theory of Quantum Information* (Cambridge,
2018); M. L. Mehta, *Random Matrices* (Elsevier, 2004). Den physikalischen Hintergrund in
mathematischer Sprache bietet [Physik für Mathematiker](../../physics_for_mathematicians.de.md).
Die Begleitarbeit zur Bibliothek ist [`paper/main.pdf`](https://github.com/andre-a-no/von_neumann_type_i/blob/main/paper/main.pdf).
