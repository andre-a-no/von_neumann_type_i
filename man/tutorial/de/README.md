[English](../README.md) · [Español](../es/README.md) · [Français](../fr/README.md) · **Deutsch** · [中文](../zh/README.md) · [日本語](../ja/README.md) · [Русский](../ru/README.md)

# Tutorial: torch_vn_algebra im Kontext der Vorlesungsnotizen

Dieses Tutorial zeigt, wie man in der Bibliothek `torch_vn_algebra` mit den endlichdimensionalen von-Neumann-Algebren
$`\bigoplus_cM_{k_c}(\mathbb C)`$ rechnet, im Kontext der [Vorlesungsnotizen](../../lecture_notes/de/README.md). Jedes Tutorial begleitet
das Kapitel der Vorlesungsnotizen mit derselben Nummer: Es zeigt den Bibliothekscode zu den dort bewiesenen Ergebnissen,
überprüft den Code an einem exakten Ergebnis (einem Satz, einer geschlossenen Formel oder einer unabhängigen Rechnung)
und nennt in einem Abschnitt *Womit man vergleichen sollte* die Standardwerkzeuge für dieselbe Aufgabe. Jeder Abschnitt
beginnt mit einer Zeile *Vorlesungsnotizen:*, die auf den Abschnitt und die von ihm verwendeten Ergebnisse verweist.
Zahlen, die nur berechnet werden (etwa die Niveaustatistik für großes $`n`$), sind als *numerische Beobachtung*
gekennzeichnet und werden nie als Prämissen oder Testziele verwendet. Neben den
[Voraussetzungen](../../lecture_notes/de/README.md#prerequisites) der Vorlesungsnotizen setzt das Tutorial Python voraus.

Der gesamte Code des Tutorials wird von der Testsuite (`tests/test_man.py`) ausgeführt;
jedes Beispiel funktioniert also mit der aktuellen Version der Bibliothek.

## Wie man den Code ausführt

```bash
pip install -e .                # from the root of the repository
python -m pytest tests/test_man.py -q
```

Alle Codeblöcke eines Kapitels können der Reihe nach in eine einzige Python-Sitzung oder
ein Jupyter-Notebook eingefügt werden. Die Blöcke eines Kapitels teilen sich ihre Variablen.

## Inhalt

| | Kapitel | Hauptergebnisse |
|---|---|---|
| 0 | [Worum es in diesem Kurs geht](00_overview.md) | die Algebra als Datentyp: Blockstruktur = gestapelte Tensoren; Aufwand; was es gibt und was neu ist |
| 1 | [Matrixalgebren und Funktionalkalkül](01_matrix_algebras.md) | Spektralsatz, Funktionalkalkül, Polarzerlegung, Dualität der Spurnorm, ein GUE-Moment |
| 2 | [Endlichdimensionale von-Neumann-Algebren](02_von_neumann_algebras.md) | Bikommutantensatz, Struktursatz, Zentrum, Symmetriegruppen und Kommutanten, bedingte Erwartungen |
| 3 | [Spuren, Zustände und Entropie](03_traces_and_states.md) | Klassifikation der Spuren, Dichten, Entropie der Sektoren, Kleinsche Ungleichung, Monotonie unter der bedingten Erwartung auf das Zentrum, Michelson-Kontrast |
| 4 | [Symmetrien und Darstellungen](04_symmetries.md) | Schursches Lemma, Ladungssektoren, $`U(1)`$ und $`SU(2)`$, Spinsysteme, reeller, komplexer und quaternionischer Typ |
| 5 | [Zufallsmatrizen in einer Algebra](05_random_matrices.md) | Haar-Maß auf $`\prod_cU(k_c)`$, Weingarten-Kalkül, Niveaustatistik und warum Sektoren nicht vermischt werden dürfen |
| 6 | [Quantenkanäle](06_quantum_channels.md) | Kraus und Stinespring, duale Abbildungen bezüglich einer Spur, Abbildungen zwischen Sektoren, Dekohärenz |
| 7 | [Tensorprodukte und Verschränkung](07_tensor_products.md) | partielle Spur, Pagesche Formel, Fusion von Ladungen |
| 8 | [Dynamik](08_dynamics.md) | Schrödinger- und Lindblad-Gleichung, Sprünge zwischen Sektoren, Postselektion |
| 9 | [Spinketten](09_spin_chains.md) | Sektoren von Magnetisierung und Impuls, Lanczos, Überprüfungen mit dem Bethe-Ansatz, Unordnung |
| 10 | [Quantencomputing](10_quantum_computing.md) | dekohärenzfreie Unterräume, Fehlerkorrektur mit Operatoralgebren, Mixer mit Nebenbedingungen |
| 11 | [Optimierung und Spurungleichungen](11_optimisation.md) | Parametrisierungen mit Nebenbedingungen, Stichprobenziehung versus Optimierung |
| 12 | [Äquivariantes maschinelles Lernen](12_equivariant_learning.md) | äquivariante Schichten als Kommutanten |
| 13 | [Jordan-Algebren: warum die Quantenmechanik *-Algebren verwendet](13_jordan_algebras.md) | Jordan-Produkt und quadratische Darstellung, Spinfaktoren und Lorentz-Kegel, die Klassifikation von Jordan–von Neumann–Wigner, Transposition und partielle Transposition, eine Spurungleichung in Jordan-Form |
| A | [Tatsachen aus Analysis und Wahrscheinlichkeitstheorie](A_analysis.md) | Potenzreihen und Nullmengen, Gaußsche Vektoren, sphärische Maße, Eigenwertungleichungen, Isserlis–Wick |
| B | [Schur-Weyl-Dualität und der Weingarten-Kalkül](B_schur_weyl.md) | Schur-Weyl-Dualität, Weingarten-Funktionen, Satz von Wigner–Eckart, Haar-Maß auf $`O(k)`$ |
| C | [Eigenwertdichten und die Pagesche Formel](C_eigenvalue_densities.md) | Weylsche Integrationsformel, GOE/GUE, Wishart- und induzierte Dichten, Pagesche Formel |
| D | [Generatoren quantendynamischer Halbgruppen](D_gksl.md) | der GKSL-Satz, Generatoren auf direkten Summen, Satz von Frigerio, Bény–Kempf–Kribs |
| E | [Klassifikation der euklidischen Jordan-Algebren](E_jordan_classification.md) | Peirce-Zerlegung, Sätze von Hurwitz und Frobenius, Koordinatisierung, die Albert-Algebra, der Satz von Jordan–von Neumann–Wigner |
| F | [Symmetrische Kegel, Ordnungsisomorphismen und Jordan-Abbildungen](F_cones_and_order.md) | Kegel der Quadrate, der Satz von Koecher–Vinberg, der Satz von Kadison, Jordan-Automorphismen von $`\bigoplus_cM_{k_c}`$ |

Die Mathematik – Definitionen, Sätze und vollständige Beweise – steht in den [Vorlesungsnotizen](../../lecture_notes/de/README.md),
zusammen mit den Voraussetzungen, weiterführender Literatur und der Begleitarbeit zur Bibliothek.
