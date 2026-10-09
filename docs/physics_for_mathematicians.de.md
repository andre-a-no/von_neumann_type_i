[English](physics_for_mathematicians.md) · [Español](physics_for_mathematicians.es.md) · [Français](physics_for_mathematicians.fr.md) · **Deutsch** · [中文](physics_for_mathematicians.zh.md) · [日本語](physics_for_mathematicians.ja.md) · [Русский](physics_for_mathematicians.ru.md)

# Physik für Mathematiker: ein Wörterbuch zu `torch_vn_algebra`

*Übersetzung der englischen Fassung; bei Abweichungen ist die englische Fassung maßgeblich.*

Dieser Leitfaden richtet sich an Leserinnen und Leser aus den Operatoralgebren, die wissen möchten, was die
physikalischen Modelle in der Bibliothek und im Artikel bedeuten. Alles wird in der Sprache der Algebra
M = ⊕_c M_{k_c}(C) formuliert, die auf H = ⊕_c C^{k_c} wirkt. Zu jedem physikalischen Begriff geben wir an,
was er in M ist und wo er sich in der Bibliothek findet.

---

## 0. Der Grundgedanke

Eine physikalische Symmetrie ist eine unitäre Darstellung π einer Gruppe G auf H. Alles, was mit der
Symmetrie verträglich ist (der Hamiltonoperator, die Zeitentwicklung, die Observablen unter einer
Superauswahlregel), liegt in der **Kommutante** π(G)'. Nach dem Bikommutantensatz und dem Schurschen Lemma gilt

    H ≅ ⊕_λ V_λ ⊗ C^{m_λ},     π(G)' ≅ ⊕_λ M_{m_λ}(C),

also ist π(G)' eine endlichdimensionale Von-Neumann-Algebra vom Typ I. Ihre Sektoren sind die irreduziblen
Darstellungen λ, die Blockgrößen sind die Multiplizitäten m_λ, und das Zentrum besteht aus Funktionen der
„Ladungen“. Wer als Physiker ein System mit Symmetrie modelliert, rechnet also in einer Algebra der Form (1),
meist ohne sie so zu nennen. Für eine einzelne Erhaltungsgröße Q (einen selbstadjungierten Operator) sind die
Sektoren die Eigenwerte von Q und die Blöcke ihre Eigenräume.

Damit ist auch die Frage „Warum nicht einfach eine große Matrix?“ beantwortet. Die Einbettung M ⊂ M_N
(N = Σk_c) ist exakt, aber dabei gehen verloren:
- **Rechenaufwand**: Σk_c³ statt N³; bei C gleich großen Blöcken ist das ein Faktor C² in der Laufzeit;
- **Maß**: Die unitäre Gruppe von M ist Π U(k_c), eine Nullmenge in U(N). Ein Haar-zufälliges Unitäres
  in U(N) liegt nicht in M, und ein Gradientenschritt in M_N führt aus M heraus;
- **Spuren**: Auf M_N ist die Spur bis auf einen Faktor eindeutig, auf M bilden die treuen Spuren die
  C-parametrige Familie Σ w_c Tr_c. „Spurerhaltend“ und duale Abbildungen hängen von den Gewichten ab;
- **das Zentrum**: Die Sektorprojektionen bilden den klassischen Teil des Systems (Abschnitt 2).
  Bedingte Wahrscheinlichkeiten und Postselektion sind Operationen mit dem Zentrum; arbeitet man mit
  einer einzigen Matrix, muss man sie aus Blockindizes rekonstruieren.

Der Fall C = 1 wird ebenfalls unterstützt: Er ist gebatchte lineare Algebra in M_n.

---

## 1. Zustände, Observablen, Messungen

| Physik | Mathematik | in der Bibliothek |
|---|---|---|
| Observable (Energie, Spin, Teilchenzahl) | selbstadjungiertes A ∈ M | `Operator` |
| mögliche Messergebnisse | Spektrum von A | `eigh` |
| (gemischter) Zustand | normales positives Funktional ω, ω(1) = 1; Dichte ρ ≥ 0, Tr ρ = 1 | `DensityMatrix` |
| reiner Zustand, „Wellenfunktion“ | Vektor ψ ∈ H, ρ = \|ψ⟩⟨ψ\| (eine minimale Projektion) | `vector_in_sector`, `basis_state` |
| Erwartungswert | ω(A) = Tr(ρA) | `expectation` |
| Wahrscheinlichkeit des Ergebnisses λ | ω(P_λ), P_λ die Spektralprojektion | `apply_function` |
| Zustand nach einer Messung (Lüders-Regel) | ρ ↦ PρP / Tr(Pρ) | `lueders_update`, `condition_on` |
| Von-Neumann-Entropie | S(ρ) = −Tr ρ log ρ | `entropy` |
| Temperatur, thermischer Zustand | ρ = e^{−βH}/Tr e^{−βH} (Gibbs-Zustand = KMS-Zustand) | `gibbs_state` |
| unendliche Temperatur | Spurzustand τ / τ(1) | `tracial_state` |

**Physikalische Bedeutung der drei Spuren.**
- Tr_blunt (Gewichte 1) ist die physikalische Spur auf H. Normiert ist sie der Zustand bei unendlicher
  Temperatur: Alle N Basiszustände sind gleich wahrscheinlich, und Sektor c hat die Wahrscheinlichkeit k_c/N.
- Tr_norm (Gewichte 1/k_c): der Gleichverteilungszustand (mikrokanonische Zustand) innerhalb jedes einzelnen
  Sektors.
- τ_vN (Gewichte 1/(C k_c)): die gleichgewichtete Mischung der mikrokanonischen Zustände der Sektoren. Jeder
  *Sektor* (nicht jeder Basiszustand) hat die Wahrscheinlichkeit 1/C.

Die Einschränkung eines beliebigen Zustands auf das Zentrum ist eine klassische Wahrscheinlichkeitsverteilung
auf den Sektoren (`sector_probabilities`). Das ist es, was Physiker „die Wahrscheinlichkeit der Ladung q“ nennen.

---

## 2. Superauswahl, Dekohärenz, klassisch und quantenmechanisch

- **Superauswahlregel** (Wick–Wightman–Wigner): Superpositionen von Zuständen mit verschiedener Ladung
  (elektrische Ladung, Fermionenparität) sind nicht beobachtbar. Mathematisch bilden die zulässigen
  Observablen die Kommutante der Ladung, also die Algebra M und nicht ganz B(H).
- **Das Zentrum von M** besteht aus Größen, die sich störungsfrei messen lassen und mit allem kommutieren:
  den „klassischen“ Variablen. Sind alle k_c = 1, ist das System rein klassisch (eine kommutative Algebra,
  d. h. ein Wahrscheinlichkeitsraum auf C Punkten). Ist C = 1, ist es rein quantenmechanisch. Der allgemeine
  Fall ist ein Hybrid.
- **Dekohärenz / Einselektion** (Zurek): Die Wechselwirkung mit einer Umgebung unterdrückt die
  Nebendiagonalblöcke ρ_{cd}, c ≠ d. Das ist die bedingte Erwartung E: B(H) → M, ρ ↦ Σ P_c ρ P_c, und die
  effektive Observablenalgebra schrumpft auf M. In der Bibliothek liegen die Operatoren bereits in M; die
  Rechnungen beginnen also dort, wo diese bedingte Erwartung bereits angewandt ist.

---

## 3. Dynamik abgeschlossener Systeme

| Physik | Mathematik | in der Bibliothek |
|---|---|---|
| Hamiltonoperator | selbstadjungiertes H ∈ M | — |
| Schrödinger-Gleichung i dψ/dt = Hψ | ψ(t) = e^{−itH}ψ(0) | `dynamics.schrodinger`, `propagator` |
| Von-Neumann-Gleichung dρ/dt = −i[H, ρ] | ρ(t) = U ρ U*, U = e^{−itH} | `von_neumann` |
| Heisenberg-Bild | A(t) = U* A U (ein Automorphismus von M) | `solve_operator_ode` |
| Erhaltungsgröße | Q mit [H, Q] = 0 | Sektoren |

Ist H ∈ M, so ist auch e^{−itH} ∈ M, und die Dynamik verlässt die Sektoren nie. Jeder Sektor entwickelt sich
unabhängig, und die Exponentialfunktion wird blockweise berechnet.

---

## 4. Offene Systeme und Kanäle

- Ein **Quantenkanal** ist eine vollständig positive, spurerhaltende Abbildung (Schrödinger-Bild). Ihre
  Adjungierte (Heisenberg-Bild) ist vollständig positiv und unital. Kraus-Darstellung:
  Φ(ρ) = Σ_i K_i ρ K_i*. In der Bibliothek: `Channel`.
- Ein **Kanal zwischen Sektoren** (`InterSectorChannel`) beschreibt einen Prozess, der die Ladung ändert:
  Teilchenverlust, Zerfall, das Heraustunneln eines Elektrons. Seine Kraus-Blöcke sind K^{dc}: C^{k_c} → C^{m_d}.
  Spurerhaltung und duale Abbildung hängen von den Spurgewichten auf Ein- und Ausgangsseite ab; daher die
  Faktoren w_out/w_in.
- Die **Lindblad-Gleichung** (GKSL) erzeugt eine Halbgruppe von Kanälen e^{tL}:
  dρ/dt = −i[H, ρ] + Σ_j γ_j (L_j ρ L_j* − ½{L_j* L_j, ρ}).
  Die L_j heißen „Sprungoperatoren“; jeder beschreibt eine Art, wie die Umgebung auf das System einwirkt.
  - **T1** (Energierelaxation, Amplitudendämpfung): L = σ⁻ = |0⟩⟨1|, ein Qubit zerfällt vom angeregten in
    den Grundzustand. Die Zahl der Anregungen sinkt um eins, es handelt sich also um einen Sprung **zwischen**
    Sektoren. Beispiel: `site_amplitude_damping`, `lowering`.
  - **T2** (Dephasierung): L = σ^z. Dieser Operator ist diagonal und ändert den Sektor nicht; er zerstört nur
    Phasen. Das ist ein Sprung **innerhalb** eines Sektors.
  - T1 und T2 sind die Standard-Kohärenzzeiten, die für jeden Quantencomputer angegeben werden.

---

## 5. Spinketten

- **Spin ½ an einem Gitterplatz**: H_site = C², Basis |↑⟩, |↓⟩ (oder |1⟩, |0⟩). Pauli-Matrizen σ^x, σ^y, σ^z;
  S^a = σ^a/2; S^± = S^x ± iS^y sind die Auf- und Absteigeoperatoren.
- **Eine Kette aus L Spins**: H = (C²)^{⊗L}, dim = 2^L. Dieses exponentielle Wachstum ist der Grund, warum
  Physiker auf Symmetrien angewiesen sind.
- **XXZ-Hamiltonoperator**: H = Σ_⟨ij⟩ [J/2 (S_i⁺S_j⁻ + S_i⁻S_j⁺) + JΔ S_i^z S_j^z] + Σ_i h_i S_i^z.
  Der erste Term vertauscht benachbarte ↑↓ (ein „Teilchenhüpfen“), der zweite ist ihre Wechselwirkung, der
  dritte ein äußeres Feld. Δ = 1 ist das Heisenberg-Modell (ein Magnet), Δ = 0 ist äquivalent zu freien
  Fermionen (über die Jordan–Wigner-Transformation), was exakte Vergleichswerte liefert.
- **Erhaltene Ladung**: S^z_tot bzw. N = Anzahl der ↑ (das Hamming-Gewicht eines Bitstrings). Jeder Term von
  H erhält sie, also gilt H ∈ ⊕_{N=0}^{L} M_{C(L,N)}: Die Blöcke haben Binomialgrößen.
- **Fusion von Sektoren** (`fused_tensor_product`): Unter dem Tensorprodukt addieren sich die Ladungen. Das ist
  die Fusionsregel von U(1): (⊕_q M_{n_q}) ⊗ (⊕_r M_{m_r}) ⊂ ⊕_s M_{Σ_{q+r=s} n_q m_r}. `SpinChain` baut seine
  Basis auf diese Weise auf, Gitterplatz für Gitterplatz.
- **Impuls** (`momentum_algebra`): Ein Ring besitzt zusätzlich die Translationssymmetrie T (die Gruppe Z_L).
  Gemeinsame Eigenräume (N, k), wobei e^{ik} ein Eigenwert von T ist, liefern eine feinere Zerlegung mit
  Blöcken, die etwa um den Faktor L kleiner sind.
- Der **Grundzustand** ist der Eigenvektor von H zum kleinsten Eigenwert, d. h. er beschreibt das Verhalten
  am absoluten Nullpunkt. Für große Blöcke bestimmt man ihn mit dem Lanczos-Verfahren (`krylov`).
- Was im Artikel überprüft wird und warum das Standardtests sind:
  - Grundzustandsenergien des Heisenberg-Rings (bekannte exakte Werte, Bethe-Ansatz);
  - die **Marshall-Vorzeichenregel**: Für das Heisenberg-Modell auf einem bipartiten Gitter hat der
    Grundzustand die Vorzeichen (−1)^{number of ↑ on one sublattice}. Das ist ein Satz und damit ein sauberer Test;
  - **SU(2)-Verschachtelung**: Das Heisenberg-Modell besitzt die volle SU(2)-Symmetrie, daher ist für N < L/2
    das Spektrum von Sektor N in dem von Sektor N+1 enthalten (Multipletts);
  - **Schmelzen einer Domänenwand**: |↑…↑↓…↓⟩ zerfließt bei Δ = 0 mit einem bekannten Profil;
  - **Peschel-Formel**: die Verschränkungsentropie freier Fermionen aus der Korrelationsmatrix.
- **Unordnung und Chaos** (`random_field_heisenberg`, `level_spacing_ratio`): Zufallsfelder
  h_i ∈ [−W, W]. Für kleines W ist das System „chaotisch“, und die Niveaus stoßen sich ab wie im GOE,
  ⟨r⟩ ≈ 0.53. Für großes W lokalisiert es (Vielteilchenlokalisierung, MBL), und die Niveaus sind unabhängig
  wie bei einem Poisson-Prozess, ⟨r⟩ ≈ 0.386. Niveaustatistik ist nur **innerhalb eines Sektors** sinnvoll:
  Das Vermischen von Sektoren erzeugt künstlich Poisson-Statistik. Das ist ein gutes physikalisches Argument
  dafür, in der Algebra mit ihren Sektoren zu arbeiten.

---

## 6. Verschränkung und Teilsysteme

- Ein zusammengesetztes System A+B: H_A ⊗ H_B, Algebra M_A ⊗ M_B.
- Der **reduzierte Zustand** ρ_A = Tr_B ρ ist die Einschränkung des Funktionals ω auf die Unteralgebra M_A ⊗ 1.
  Auf der Ebene der Dichten ist er eine bedingte Erwartung (die partielle Spur). In der Bibliothek:
  `partial_trace`, `reduced_state`.
- Die **Verschränkungsentropie** S(ρ_A) eines reinen ρ misst die Quantenkorrelation zwischen A und B.
- **Page-Formel**: die exakte mittlere Entropie eines Teilsystems eines zufälligen reinen Zustands. Mit ihr
  wird geprüft, dass der Sampler Haar-verteilt ist.
- **Symmetrieaufgelöste Verschränkung**: Ist die Ladung erhalten, so ist ρ_A blockdiagonal, und seine Entropie
  zerfällt in einen „klassischen“ Anteil (die Entropie der Verteilung auf dem Zentrum) und einen
  „quantenmechanischen“ Anteil (die mittlere Entropie der Blöcke). Das ist eine direkte Anwendung des Zentrums.

---

## 7. Quantencomputing

- **Qubit** = C², n Qubits = (C²)^{⊗n}; ein **Gatter** ist ein unitärer Operator; ein **Schaltkreis** ist ein
  Produkt von Gattern; eine **Messung** in der Rechenbasis verwendet die Spektralprojektionen von σ^z.
- **QAOA** (quantum approximate / alternating operator ansatz): Man wechselt ab zwischen e^{−iγ H_C} (H_C ist
  diagonal und kodiert die Zielfunktion) und e^{−iβ H_M} (dem „Mixer“). Für Probleme mit der Nebenbedingung
  „genau k Einsen“ (wähle k von n Objekten) verwendet man den **XY-Mixer** Σ(X_iX_j + Y_iY_j) =
  2Σ(σ⁺σ⁻ + h.c.). Er erhält das Hamming-Gewicht, sodass der gesamte Schaltkreis in ⊕_N M_{C(n,N)} liegt und
  die zulässigen Lösungen genau den Sektor N = k bilden.
- **Rauschen**: T1 führt aus dem zulässigen Sektor heraus (eine Eins geht verloren), T2 nicht.
- **Postselektion**: Man misst das Hamming-Gewicht und verwirft fehlerhafte Durchläufe. Das ist die
  Lüders-Regel mit einer **zentralen** Projektion P_k: ρ ↦ P_kρP_k / Tr(P_kρ), mit Erfolgswahrscheinlichkeit
  ω(P_k). Genau das macht `ex9_qaoa_xy_mixer.py`. Bei reinem T1-Rauschen liefert die Postselektion exakt den
  idealen Zustand (F = 1), weil jeder T1-Sprung den Sektor verlässt.
- Dekohärenzfreie Unterräume und operatoralgebraische Fehlerkorrektur: Geschützte Information wird in einem
  Block (oder einem Tensorfaktor eines Blocks) einer Algebra kodiert, auf die das Rauschen nicht einwirkt.
  Auch das ist die Sprache der Typ-I-Algebren (Knill–Laflamme–Viola, Bény–Kempf–Kribs).

---

## 8. Zufallsmatrizen, reell und komplex

Dieser Abschnitt hängt direkt mit dem Parameter `complex_valued` zusammen.
- Eine **antiunitäre Symmetrie** ist die Zeitumkehr T. Gilt T² = 1 (spinlose Teilchen, kein Magnetfeld), so ist
  H in einer geeigneten Basis **reell**: GOE, orthogonale Ensembles, COE. Reelle Arithmetik genügt; sie halbiert
  den Speicherbedarf und ist schneller (siehe die Laufzeittabelle im Artikel).
- Ohne eine solche Symmetrie (Magnetfeld, komplexe Phasen, Impuls k ≠ 0, π) ist H komplex: GUE, CUE.
  Komplexe Zahlen sind dann unverzichtbar.
- Gilt T² = −1 (halbzahliger Spin mit Spin-Bahn-Kopplung), liegt der quaternionische Fall vor: GSE, CSE.
- Das XXZ-Modell mit reellen Feldern ist reell. Impulssektoren mit k ∉ {0, π} sind komplex. Die Dynamik
  e^{−itH} ist stets komplex, auch für reelles H. Das ist die Tabelle „was wann genügt“ im Artikel.
- Doppelte Genauigkeit wird benötigt, wenn es auf kleine Unterschiede ankommt: entartete Spektren, Identitäten,
  die bis 1e−12 geprüft werden, lange Zeitentwicklungen, Lanczos mit vielen Schritten. Einfache Genauigkeit
  (und TF32 auf Tensor Cores) genügt für Statistiken über viele Stichproben.

---

## 9. Maschinelles Lernen

Eine äquivariante lineare Schicht ist eine Abbildung W: V → V, die mit einer Gruppendarstellung vertauscht
(Translationen bei Faltungen, Permutationen bei Graph-Netzen, Drehungen bei Molekülen). Nach dem Schurschen
Lemma ist der Raum solcher W gleich ⊕_λ M_{m_λ} mit den Multiplizitäten m_λ; eine äquivariante Schicht wird
also durch ein Element einer Algebra der Form (1) parametrisiert. Die Bibliothek stellt für solche
Parametrisierungen zufällige Initialisierung mit dem richtigen Maß, spektrale Nebenbedingungen und Gradienten
bereit.

---

## 10. Einordnung: was die Bibliothek ist und was nicht

- Die Bibliothek liefert **keine** neue Physik. Sie ist ein validiertes Werkzeug für Rechnungen in Algebren mit
  Sektoren. Die physikalischen Modelle im Artikel sind bekannte Systeme mit bekannten Ergebnissen, die zur
  Überprüfung des Werkzeugs dienen.
- Für ein einzelnes großes System (ein einzelner Grundzustand bei L = 30, Tensornetzwerke, Simulation großer
  Schaltkreise) sind die spezialisierten Pakete (QuSpin, ITensor/TeNPy, Qiskit) die bessere Wahl.
- Ihre Nische sind viele mittelgroße Operatoren mit expliziter Sektorstruktur: Unordnungsensembles, Statistiken
  über zufällige Zustände und Kanäle, extremale Operatoren unter spektralen Nebenbedingungen, Kanäle zwischen
  Sektoren mit konsistenten Spuren, und das zusammen mit automatischem Differenzieren und GPUs.
- Die ursprüngliche mathematische Fragestellung (Spurungleichungen, Michelson-Kontrast) ist ein Beispiel für ein
  Extremalproblem, bei dem Monte-Carlo-Sampling versagt und Optimierung unter Nebenbedingungen funktioniert.

---

## Weiterführende Literatur

- M. A. Nielsen, I. L. Chuang, *Quantum Computation and Quantum Information* — Kap. 2 (Postulate),
  Kap. 8 (Kanäle, T1/T2).
- M. M. Wilde, *Quantum Information Theory* — Kanäle und Zustände, nah an den Operatoralgebren.
- A. W. Sandvik, "Computational studies of quantum spin systems", AIP Conf. Proc. 1297 (2010) —
  Sektoren und exakte Diagonalisierung in der Praxis.
- H.-P. Breuer, F. Petruccione, *The Theory of Open Quantum Systems* — die Lindblad-Gleichung.
- F. Haake, *Quantum Signatures of Chaos* — Zufallsmatrizen, Zeitumkehr, GOE/GUE/GSE.
- S. D. Bartlett, T. Rudolph, R. W. Spekkens, Rev. Mod. Phys. 79, 555 (2007) — Superauswahl und
  Bezugssysteme, gut lesbar für Mathematiker.
