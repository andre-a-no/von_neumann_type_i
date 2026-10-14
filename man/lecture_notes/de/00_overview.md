[English](../00_overview.md) · [Español](../es/00_overview.md) · [Français](../fr/00_overview.md) · **Deutsch** · [中文](../zh/00_overview.md) · [日本語](../ja/00_overview.md) · [Русский](../ru/00_overview.md)

# 0. Worum es in diesem Kurs geht

Von-Neumann-Algebren werden meist mit Papier und Bleistift untersucht. In endlicher Dimension
ist jede von ihnen $`*`$-isomorph zu einer direkten Summe voller Matrixalgebren (Satz 2.6), sodass
sich im Prinzip alles berechnen lässt – wird eine solche Algebra jedoch numerisch behandelt, so
wird sie oft zu einer einzigen großen Matrix verflacht, und ihre Struktur (die Sektoren, das
Zentrum, die vielen Spuren) verschwindet aus der Rechnung. Dieser Kurs und die Bibliothek
`torch_vn_algebra` gehen den umgekehrten Weg: **Die Algebra selbst ist der Datentyp**, und
Rechnungen in ihr werden effizient ausgeführt.

<a id="s-0-1"></a>

## 0.1 Die Beobachtung

Eine von-Neumann-Algebra $`\mathcal M`$ auf $`\mathbb C^n`$ ist bis auf einen unitären Basiswechsel
```math
\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c}
```
(Satz 2.6). Als abstrakte Algebra ist sie $`\bigoplus_cM_{k_c}(\mathbb C)`$: Ein Element ist eine
Liste von $`C`$ quadratischen Blöcken. Das ist genau die Gestalt eines Tensors mit einer
*Sektorachse*: Ein Element von $`\mathcal M`$ ist ein Array der Gestalt $`(C,k,k)`$ (die Blöcke mit
Nullen aufgefüllt auf $`k=\max_ck_c`$; so speichert es die Bibliothek, Abschnitt 2.5), und $`B`$
Elemente – Monte-Carlo-Stichproben, Startpunkte einer Optimierung, Zeitschritte, Realisierungen
der Unordnung – bilden ein Array der Gestalt $`(B,C,k,k)`$. Moderne Array-Bibliotheken (PyTorch auf
CPUs und GPUs) stellen gestapelte Routinen – Matrixprodukte, Eigenzerlegungen, QR-Zerlegungen –
gerade für solche „Stapel kleiner Matrizen“ bereit, die auch neuronale Netze verwenden.

*In der Bibliothek:* Ein Element von $`\mathcal M`$ ist ein Operator einer `TypeIAlgebra`, gespeichert als ein Tensor der Gestalt $`(B,C,k_{\max},k_{\max})`$ — [Tutorial, Abschnitt 1.5](../../tutorial/de/01_matrix_algebras.md#s-1-5) und [Tutorial, Abschnitt 2.5](../../tutorial/de/02_von_neumann_algebras.md#s-2-5).

Unter dieser Identifikation sind die Operationen der Algebra gestapelte Array-Operationen. Die
Entsprechung ist exakt (jede der ersten sechs Zeilen der Tabelle ist eine Identität, die unten
aus Ergebnissen des Kurses bewiesen wird); nur die numerischen Routinen, die die Operationen
ausführen, arbeiten in Gleitkommaarithmetik:

| in der Algebra $`\mathcal M`$ | auf dem Array der Gestalt $`(B,C,k,k)`$ |
|---|---|
| Produkt $`AB`$, Summe, Adjungierte | gestapeltes Matrixprodukt, Summe, konjugierte Transposition über $`(B,C)`$ |
| Funktionalkalkül $`f(A)`$ (Kapitel 1) | gestapelte Eigenzerlegung der Blöcke |
| eine Spur $`\sum_cw_c\mathrm{Tr}\,A_c`$ (Kapitel 3) | Diagonalsumme pro Block, gewichtete Summe über die Sektorachse |
| das Zentrum, zentrale Projektionen (Kapitel 2) | eine Zahl pro Sektor |
| bedingte Erwartung auf das Zentrum | Blockspur geteilt durch $`k_c`$, zurück auf die Blöcke verteilt |
| Haar-Maß auf der unitären Gruppe $`\prod_cU(k_c)`$ | gestapelte QR-Zerlegung Gaußscher Blöcke, mit der Phasenkorrektur aus Satz 5.3 |
| Zustände, Kanäle, Dynamik (Kapitel 3, 6, 8) | aus den obigen Operationen aufgebaut |

*Beweise der ersten sechs Zeilen.* Nach Satz 2.6 gibt es ein Unitäres $`V`$ mit
$`\mathcal M=V\big(\bigoplus_cM_{k_c}\otimes1_{m_c}\big)V^*`$, und
$`\varphi\big(V(\bigoplus_cA_c\otimes1_{m_c})V^*\big)=\bigoplus_cA_c`$ ist ein $`*`$-Isomorphismus von $`\mathcal M`$ auf
$`\bigoplus_cM_{k_c}(\mathbb C)`$ (wohldefiniert, da $`A_c\otimes1_{m_c}`$ das $`A_c`$ bestimmt; multiplikativ, da
$`(A_c\otimes1)(B_c\otimes1)=A_cB_c\otimes1`$). Die Bibliothek speichert $`\varphi(A)`$ (Abschnitt 2.5).
(1) In $`\bigoplus_cM_{k_c}`$ werden Summen, Produkte und Adjungierte blockweise gebildet.
(2) Seien $`A\in\mathcal M`$ normal, $`\varphi(A)=\bigoplus_cA_c`$ und $`f`$ eine Funktion auf $`\mathrm{spec}\,A`$. Jedes $`A_c`$ ist
normal, da $`\varphi(A)^*\varphi(A)=\varphi(A^*A)=\varphi(AA^*)=\varphi(A)\varphi(A)^*`$ blockweise gilt. Ist
$`A_cx=\lambda x`$ mit $`x\neq0`$, so ist $`V(0\oplus\dots\oplus x\otimes e_1\oplus\dots\oplus0)`$ ($`x\otimes e_1`$ im
$`c`$-ten Summanden) ein Eigenvektor von $`A`$ zu $`\lambda`$, also $`\mathrm{spec}\,A_c\subseteq\mathrm{spec}\,A`$. Sei $`p`$ ein
Polynom mit $`p=f`$ auf $`\mathrm{spec}\,A`$ (Lagrange-Interpolation, wie im Beweis von Satz 1.5). Für jedes normale
$`X`$ mit $`\mathrm{spec}\,X\subseteq\mathrm{spec}\,A`$ und $`X=U\Lambda U^*`$ gilt $`f(X)=Uf(\Lambda)U^*=Up(\Lambda)U^*=p(X)`$ nach
Satz 1.5(b); also $`f(A)=p(A)`$, $`f(A_c)=p(A_c)`$, und da $`\varphi`$ ein unitaler Homomorphismus ist,
```math
\varphi(f(A))=p(\varphi(A))=\bigoplus_cp(A_c)=\bigoplus_cf(A_c),\qquad f(A_c)=U_cf(\Lambda_c)U_c^*
```
für jede Diagonalisierung $`A_c=U_c\Lambda_cU_c^*`$ (Satz 1.3): eine Eigenzerlegung pro Block.
(3) Die Spuren von $`\bigoplus_cM_{k_c}`$ sind die Funktionale $`\sum_cw_c\mathrm{Tr}\,A_c`$ (Satz 3.1).
(4) $`Z(\mathcal M)`$ wird von den minimalen zentralen Projektionen $`P_c`$ aufgespannt, und $`\varphi(P_c)`$ ist $`1_{k_c}`$ im Block
$`c`$ und $`0`$ sonst (Satz 2.6), also $`\varphi(Z(\mathcal M))=\{\bigoplus_cz_c1_{k_c}:z\in\mathbb C^C\}`$.
(5) Proposition 2.8. (6) Korollar 5.4 (unabhängige Blöcke, jeder Haar-verteilt auf $`U(k_c)`$) und Satz 5.3
(jeder Block aus der QR-Zerlegung einer Gaußschen Matrix). Die letzte Zeile ist nur ein Verweis: Zustände
werden durch Dichten dargestellt (Satz 3.2), Kanäle zwischen Sektoren durch Kraus-Operatoren zwischen Blöcken
(Proposition 6.4) und die geschlossene Dynamik durch $`e^{-itH}=\bigoplus_ce^{-itH_c}`$ (Proposition 8.2).

*In der Bibliothek:* Produkte und Funktionen wie `sqrt` und `apply_function` — [Tutorial, Abschnitt 1.5](../../tutorial/de/01_matrix_algebras.md#s-1-5); die Spuren, z. B. `Tr_blunt` — [Tutorial, Abschnitt 3.1](../../tutorial/de/03_traces_and_states.md#s-3-1); das Zentrum, `TypeIAlgebra.central` — [Tutorial, Abschnitt 2.5](../../tutorial/de/02_von_neumann_algebras.md#s-2-5); die bedingte Erwartung auf das Zentrum, `channels.center_expectation` — [Tutorial, Abschnitt 2.6](../../tutorial/de/02_von_neumann_algebras.md#s-2-6); Haar-verteilte unitäre Elemente, `random_unitary_operator` — [Tutorial, Abschnitt 5.1](../../tutorial/de/05_random_matrices.md#s-5-1).

*Rechnung:* Zeile 2 der Tabelle auf $`\bigoplus_{c=1}^{64}M_{16}(\mathbb C)`$ — [Tutorial, Abschnitt 0.2](../../tutorial/de/00_overview.md#s-0-2); was das blockweise Rechnen kostet — [Tutorial, Abschnitt 0.1](../../tutorial/de/00_overview.md#s-0-1).

<a id="s-0-2"></a>

## 0.2 Eine erste Rechnung

*Rechnung:* die Quadratwurzel und die Spur eines zufälligen positiven Elements von $`\bigoplus_{c=1}^{64}M_{16}(\mathbb C)`$, blockweise und in $`M_{1024}(\mathbb C)`$ berechnet — [Tutorial, Abschnitt 0.2](../../tutorial/de/00_overview.md#s-0-2).

<a id="s-0-3"></a>

## 0.3 Was es gibt und was neu ist

*In der Bibliothek:* verwandte Software und was die Bibliothek ihr hinzufügt — [Tutorial, Abschnitt 0.3](../../tutorial/de/00_overview.md#s-0-3).

<a id="s-0-4"></a>

## 0.4 Wie man diesen Kurs liest

| Kapitel | Inhalt |
|---|---|
| 1–3 | die Mathematik: Matrixalgebren, der Struktursatz, Spuren, Zustände und Entropie – mit vollständigen Beweisen |
| 4–8 | woher diese Algebren kommen und was man in ihnen berechnet: Symmetrien, Zufallsmatrizen, Kanäle, Tensorprodukte, Dynamik |
| 9–12 | Anwendungen: Spinketten, Quantencomputing, Optimierung, äquivariantes maschinelles Lernen |
| 13 | Grundlagen: Jordan-Algebren und warum die Quantenmechanik *-Algebren verwendet |
| Anhänge | vollständige Beweise der Hilfsresultate, die die Kapitel mit Nummer zitieren (z. B. Satz A.24) |

Jedes der Kapitel 1–13 formuliert Ergebnisse mit vollständigen Beweisen. Die Voraussetzungen (lineare Algebra,
reelle Analysis, Maßtheorie, elementare Wahrscheinlichkeitsrechnung) sind im [README](../../lecture_notes/de/README.md#prerequisites)
aufgeführt; alles Weitere wird in den Vorlesungsnotizen bewiesen, und Ergebnisse anderer Kapitel und der Anhänge
werden mit Nummer zitiert. Die Kästen mit der Überschrift *Physikalische Bedeutung* sind physikalische Deutung,
keine Sätze der Vorlesungsnotizen.

Jedes Kapitel hat ein Gegenstück im [Tutorial](../../tutorial/de/00_overview.md), das den Bibliothekscode zeigt, ihn numerisch an den hier bewiesenen Ergebnissen überprüft und die Standardwerkzeuge zum Vergleich nennt. Die Zeilen *Rechnung:* und *In der Bibliothek:* nach einem Ergebnis verweisen auf die Stelle im Tutorial, an der es berechnet oder verwendet wird.

**Weiter:** [Kapitel 1](01_matrix_algebras.md).
