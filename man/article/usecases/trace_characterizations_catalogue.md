# Inequalities that characterise the trace: a catalogue

Compiled for the use case "quantitative failure of trace characterisations" (`trace_characterizations.tex`,
`scripts/trace_characterizations.py`).

**How this list was made.** The only source was web search. We could not open mathnet.ru, Springer, ScienceDirect,
Cambridge Core or arXiv full texts from this machine. Every entry below rests on an **abstract, a publisher or
database listing, or a secondary description** (another paper's summary). None comes from the full text. The author
should check every statement against the paper. *Later update:* the full texts of Abed–Nikolaeva–Novikov (B9,
arXiv v9) and Moslehian–Bikchentaev (B2) were read; the entries that use them are marked **F**.

The column "seen" says what we saw:

- **A**: the abstract on the publisher or database page, as shown by the search engine;
- **L**: a bibliographic listing only;
- **S**: a secondary description, i.e. how another paper or listing summarises the result;
- **F**: the full text.

"Form used" is the exact inequality evaluated in `scripts/trace_characterizations.py`. Where the source does not fix
the form, we say which reading we chose.

Throughout, φ is a positive linear functional. On M_n it has the form φ = Tr(W ·) with a density W ≥ 0. "Trace" means
φ(XY) = φ(YX); on M_n this is equivalent to φ = c·Tr.

## A. Inequalities evaluated numerically

| # | inequality (form used) | operators | algebra | statement of the characterisation | reference | seen |
|---|---|---|---|---|---|---|
| 1 | \|φ(X)\| ≤ φ(\|X\|) | arbitrary X | unital C*-algebra (general) | a positive functional is a trace iff the inequality holds for all X. Gardner's Thm. 1, as quoted in Abed–Nikolaeva–Novikov (B9), p. 20: "The finite traces on a C\*-algebra A are precisely those (positive) linear functionals φ on A which satisfy \|φ(x)\| ≤ φ(\|x\|) for all x ∈ A" | L. T. Gardner, *An inequality characterizes the trace*, Canad. J. Math. 31 (1979), no. 6, 1322–1328, doi:10.4153/CJM-1979-109-9 | A, S |
| 2 | φ(\|X+Y\|) ≤ φ(\|X\|) + φ(\|Y\|) | self-adjoint X, Y (we optimise over self-adjoint ones) | von Neumann algebra (normal semifinite weights) | a normal semifinite weight satisfying the triangle inequality for all self-adjoint a₁, a₂ is a trace. Gardner (1979) states the same for maps on C*-algebras with the triangle inequality for all x, y (abstract only) | A. I. Stolyarov, O. E. Tikhonov, A. N. Sherstnev, *Characterization of normal traces on von Neumann algebras by inequalities for the modulus*, Mat. Zametki 72 (2002), no. 3, 448–454; Math. Notes 72 (2002), no. 3–4, 411–416, doi:10.1023/A:1020559623287; and Gardner (1979) | A |
| 3 | φ(A + A\*) ≤ 2φ(\|A\*\|) | arbitrary A | von Neumann algebra | stated as a corollary: tracial functionals are characterised by this inequality | A. M. Bikchentaev, S. A. Abed, *Projections and traces on von Neumann algebras*, Lobachevskii J. Math. 40 (2019), 1260–1267, doi:10.1134/S1995080219090051 | S |
| 4 | \|φ(A²)\| ≤ φ(A\*A) (also \|φ(Re A²)\| ≤ φ(A\*A)) | arbitrary A | von Neumann algebra, positive normal φ | φ is tracial ⟺ \|φ(Re(A²))\| ≤ φ(\|A\|²) ∀A ⟺ \|φ(A²)\| ≤ φ(\|A\|²) ∀A | A. M. Bikchentaev, M. S. Moslehian, *Trace inequalities and characterizations of tracial functionals in operator algebras*, Positivity (2026), doi:10.1007/s11117-026-01180-0 (volume and pages not seen) | A |
| 5 | φ(\|P−R\|) ≤ φ(\|P−Q\|) + φ(\|Q−R\|) | projections | von Neumann algebra | ρ_φ(P,Q) = φ(\|P−Q\|) (or d_φ(P,Q) = φ(P∨Q − P∧Q)) satisfies the triangle inequality iff φ is tracial. That ρ_τ is a metric needs a faithful trace | A. M. Bikchentaev, *Metrics on projections of the von Neumann algebra associated with tracial functionals*, Sib. Math. J. 60 (2019), 952–956, doi:10.1134/S003744661906003X | A |
| 6 | φ(A ♮ B) ≤ (φ(A) φ(B))^{1/2}, with A ♮ B = (A⁻¹#B)^{1/2} A (A⁻¹#B)^{1/2} | positive definite A, B | M_n | φ is a positive multiple of Tr iff the inequality holds for all positive definite A, B. The proof uses "nearly parallel pure states". The paper also characterises the trace by φ(A♮B) ≤ φ((A+B)/2) | A. M. Bikchentaev, T. H. Dinh, A. V. Le, M. S. Moslehian, *Spectral geometric mean and trace characterizations*, Positivity 30 (2026), Art. 46, doi:10.1007/s11117-026-01202-x; arXiv:2605.18888 | A |
| 7 | φ((B^{1/2}AB^{1/2})^{rq}) ≤ φ((B^{r/2}A^rB^{r/2})^q). We use **r = 2, q = 1**: φ((B^{1/2}AB^{1/2})²) ≤ φ(BA²B) | positive A, B | C*-algebra | "either the Peierls–Bogoliubov or the Araki–Lieb–Thirring inequality characterises tracial functionals among positive functionals"; the exponents for which this is proved were not seen. Moslehian–Bikchentaev (B2), p. 1, lists Peierls–Bogoliubov and Araki–Lieb–Thirring among inequalities that characterise traces on M_n *when restricted to projections* (citing Math. Notes 89 (2011), this paper and the 2024 book B10), so the published hypothesis may be on projections only | A. M. Bikchentaev, *The Peierls–Bogoliubov inequality in C\*-algebras and characterization of tracial functionals*, Lobachevskii J. Math. 32 (2011), no. 3, 175–179, doi:10.1134/S1995080211030061 | A |
| 8 | φ(\|AB\*\|) ≤ ½φ(\|A\|²) + ½φ(\|B\|²) (Young, p = q = 2, in the form of Ando's trace inequality) | arbitrary A, B | M_n | Young's inequality for φ, for conjugate exponents, characterises the trace. A search snippet of the paper gives the form φ(\|AB\|) ≤ φ(A^p)/p + φ(B^q)/q for **positive semidefinite** A, B; our form extends it to arbitrary A, B (modelled on Ando's trace form), and V = h − 1 (Prop. tcB(a)) is the envelope of this extension. For positive A, B an independent L-BFGS search (25 random starts, n = 2, 3) gives V ≈ κ − 1 (0.1546 at δ = 0.5, 1.294 at δ = 0.9), a quadratic failure; observation, not proof | A. M. Bikchentaev, O. E. Tikhonov, *Characterization of the trace by Young's inequality*, J. Inequal. Pure Appl. Math. (JIPAM) 6 (2005), no. 2, Art. 49. See also K. Cho, T. Sano, *Young's inequality and trace*, Linear Algebra Appl. 431 (2009), no. 8, 1218–1222 (cited as a characterisation by B2 and B9; not seen) | S |
| 9 | φ(PQ + QP) ≤ 2φ((QPQ)^p). We use **p = 1**: φ(PQ+QP) ≤ 2φ(QPQ) | projections | von Neumann algebra, positive normal φ | φ is tracial iff the inequality holds for all projections P, Q and some p = p(P,Q) ∈ (0, 1]. The same abstract lists φ(PQP) ≤ φ(P)^{1/p} φ(Q)^{1/q}, 1/p + 1/q = 1, p ≠ 2, as equivalent; we did not use it because the role of p ≠ 2 is unclear to us from the snippet | Bikchentaev and Abed (2019), as in #3 | S |
| 10 | 0 ≤ A ≤ B ⟹ φ(A^p) ≤ φ(B^p). We use **p = 2** | positive A ≤ B | M_n | if the implication holds for all A ≤ B (fixed p > 1), then φ is a non-negative multiple of Tr | A. M. Bikchentaev, O. E. Tikhonov, *Characterization of the trace by monotonicity inequalities*, Linear Algebra Appl. 422 (2007), no. 1, 274–278, doi:10.1016/j.laa.2006.10.005 | A |
| 11 | (det A)^{1/n} ≤ φ(A)/n with φ(1) = n (an AM–GM analogue) | positive A | M_n | if "some analog of the classical inequalities for the determinant and trace (or permanent and trace)" holds for φ with φ(1) = n, then φ = Tr. **Which determinant inequality the paper uses was not seen**; AM–GM is our choice. This item characterises the *canonical* trace, not traces in general | A. M. Bikchentaev, *Inequalities for determinants and characterization of the trace*, Sib. Math. J. 61 (2020), no. 2, 248–254, doi:10.1134/S0037446620020068 | A |

## A′. Inequalities evaluated numerically (second round: Hoa, Tikhonov, Novikov, Stolyarov, Petz–Zemánek and others)

These were added at the author's request. Items 12–16 are evaluated in the paper. The column "seen" uses the same
codes as above.

| # | inequality (form used) | operators | algebra | statement of the characterisation | reference | seen |
|---|---|---|---|---|---|---|
| 12 | φ((A+B)^{1/2}) ≤ φ(A^{1/2}) + φ(B^{1/2}). We use **f(t) = t^{1/2}** | positive A, B | von Neumann algebra, positive normal φ | the abstract asks under which assumptions on φ and on a Borel function f with f(0) = 0 the subadditivity inequality φ(f(A+B)) ≤ φ(f(A)) + φ(f(B)) holds for all positive A, B. The visible abstract breaks off there, but later papers cite it as a characterisation of tracial functionals. We use the concave choice f = t^{1/2}; for the trace this is Rotfel'd's inequality | O. E. Tikhonov, *Subadditivity inequalities in von Neumann algebras and characterization of tracial functionals*, Positivity 9 (2005), 259–264, doi:10.1007/s11117-005-2711-1 | A (truncated) |
| 13 | φ(A+B−\|A−B\|) ≤ 2φ(f(A)^{1/2} g(B) f(A)^{1/2}) with g(t) = t/f(t). We use **f = g = t^{1/2}**: φ(A+B−\|A−B\|) ≤ 2φ(A^{1/4}B^{1/2}A^{1/4}) | positive A, B | C*-algebra, positive linear φ | the generalised Powers–Størmer inequality holds for f operator monotone on [0, ∞), and it "characterizes the tracial functionals on C*-algebras". The exact form shown here comes from a slide summary, not from the paper | D. T. Hoa, H. Osaka, H. M. Toan, *On generalized Powers–Størmer's inequality*, Linear Algebra Appl. 438 (2013), no. 1, 242–249, doi:10.1016/j.laa.2012.07.053; arXiv:1204.6665 | A, S |
| 14 | \|φ(PQ)\| ≤ φ(P)^{1/p} φ(Q)^{1/q}. We use **p = 3**, q = 3/2 | projections | M_n; C*-algebras | "Hölder, Cauchy–Schwarz(–Bunyakovskii), Golden–Thompson … inequalities, when they hold only for projections, characterize tracial functionals". The exact projection form and exponents were not seen. For p = 2 the inequality is the Cauchy–Schwarz inequality and holds for every positive φ, which explains the restriction p ≠ 2 in #9 (Bikchentaev–Abed) | A. M. Bikchentaev, *Commutation of projections and trace characterization on von Neumann algebras. II*, Math. Notes 89 (2011), no. 3–4, 461–471, credited in Moslehian–Bikchentaev, J. Math. Anal. Appl. 552 (2025), p. 1 (full text read: "Hölder, Cauchy–Schwarz–Buniakowski, Golden–Thompson, Peierls–Bogoliubov, Araki–Lieb–Thirring, when restricted to projections, characterize the tracial functionals among all positive functionals ϕ on Mn"); form and exponents not given there. Compare Bikchentaev–Abed (2019), #9 | S (via F) |
| 15 | φ(e^{P+Q}) ≤ φ(e^{P/2} e^Q e^{P/2}) (symmetric Golden–Thompson) | projections | M_n | same source as #14: Golden–Thompson restricted to projections characterises tracial functionals. **Form not seen.** We do **not** attribute it to Petz–Zemánek (see B11) | as #14 (Math. Notes 89 (2011), credited in JMAA 552 (2025), p. 1) | S (via F) |
| 16 | φ(e^A) exp(φ(e^{A/2}Be^{A/2})/φ(e^A)) ≤ φ(e^{A+B}) (symmetric Peierls–Bogoliubov) | self-adjoint A, B | C*-algebra | "either the Peierls–Bogoliubov or the Araki–Lieb–Thirring inequality characterises tracial functionals among positive functionals"; it also answers a question of J. Zemánek. The form of φ(e^A B) used in the paper was not seen; we use the symmetric (real) version. As for #7, the published hypothesis may be on projections only (B2, p. 1) | Bikchentaev (2011), as #7 | A |

## B. Further characterisations (not evaluated: the form was not visible, or the statement is not an inequality)

| # | statement (as far as seen) | algebra | reference | seen |
|---|---|---|---|---|
| B1 | "Characterizations of the trace": the content was not visible | M_n | D. Petz, J. Zemánek, *Characterizations of the trace*, Linear Algebra Appl. 111 (1988), 43–52, doi:10.1016/0024-3795(88)90050-X | L |
| B2 | Thm. 2.6: if φ(A) = φ(\|A\|) ⇒ A ≥ 0 for every nonsingular A ∈ M_n, then φ = k·Tr, k > 0; Prop. 2.1 and 2.3: the same implication for dissipative and for normal A, for faithful φ. Thm. 3.3: for φ(1) = n, φ = Tr ⟺ (φ(A^{1/2}BA^{1/2}))^{1/2} ≤ φ(A+B)/2 for all A, B ≥ 0 (the root is on the value; the inequality is due to Yang, JMAA 133 (1988) 573–574). Like #11, this characterises the canonical trace | M_n | M. S. Moslehian, A. M. Bikchentaev, *Characterizations of the canonical trace on full matrix algebras*, J. Math. Anal. Appl. 552 (2025), no. 2, 129764 | F |
| B3 | characterisations of Tr on M_n by a commutator inequality for Hermitian matrices and by a generalised AM–GM inequality for all matrices and conjugate exponents (forms not seen) | M_n | M. S. Moslehian, *More characterizations of the matrix trace*, Rend. Circ. Mat. Palermo (2) 74 (2025), Art. 130, doi:10.1007/s12215-025-01253-1 | A |
| B4 | inequalities for the spread λ_max − λ_min characterise Tr among positive φ with φ(1) = n, for n ≥ 3; three further characterisations hold for n ≥ 2 | M_n | M. Dalloul, *Characterization of tracial functionals by inequalities for the spread of a matrix*, Lobachevskii J. Math. 45 (2024), 2478–2484, doi:10.1134/S1995080224602777 | A |
| B5 | an inequality characterising tracial functionals among positive normal functionals on a von Neumann algebra, "strengthening Gardner". One secondary snippet gives it as φ(A) ≤ φ(\|A + iB\|) (A, B self-adjoint). If so, it is #3 with Z = A + iB, and its envelope is κ − 1 | von Neumann algebra | H. Alhasan, K. Fawwaz, *Characterization of tracial functionals on von Neumann algebras*, Lobachevskii J. Math. 42 (2021), 2273–2279, doi:10.1134/S1995080221100024 | A, S |
| B6 | operator inequalities that are necessary and sufficient for two projections to commute, used to characterise traces among positive normal functionals (part III: arbitrary weights, and positive functionals on C*-algebras) | von Neumann algebra | A. M. Bikchentaev, *Commutativity of projections and characterization of traces on von Neumann algebras*, Sib. Math. J. 51 (2010), 971–977, doi:10.1007/s11202-010-0096-2; *… II*, Math. Notes 89 (2011), 461–471, doi:10.1134/S0001434611030175; *… III*, Int. J. Theor. Phys. 54 (2015), doi:10.1007/s10773-015-2639-6 | A |
| B7 | a new inequality for a trace on a unital C*-algebra that characterises traces among positive functionals, and a commutativity criterion | unital C*-algebra | A. M. Bikchentaev, *Inequality for a trace on a unital C\*-algebra*, Math. Notes 99 (2016), 487–491, doi:10.1134/S0001434616030214 | A |
| B8 | tracial states among all states on a Rickart C*-algebra characterised via polar decomposition; a triangle inequality for Hermitian elements and traces | Rickart C*-algebra | A. M. Bikchentaev, *Trace inequalities for Rickart C\*-algebras*, Positivity 25 (2021), 1943–1957, doi:10.1007/s11117-021-00852-3 | A |
| B9 | *related, quantitative* (full text read, arXiv v9): the Michelson contrast Δ of positive operators and the factor-wise contrast Δ′ (Sect. 4). Sect. 5.1, Fig. 1 is a Monte Carlo scatter (real, n = 2, 3, ‖X‖ = ‖Y‖ = 1) of the unnormalised Gardner violation \|Tr(XY)\| − Tr(X\|Y\|) against Δ(X), where X is the density of φ = Tr(X·), and against Δ(Y). Sect. 7 conjectures contrast-dependent upper bounds and lists monotonicity, subadditivity, Young and others as open. No envelope, bound or closed form; only Gardner (item 1) | countably decomposable W*-algebras | S. A. Abed, I. A. Nikolaeva, A. A. Novikov, *Generalisation of Michelson contrast for operators and its properties*, Lobachevskii J. Math. 45 (2024), no. 8, 3835–3848, doi:10.1134/S1995080224603400; arXiv:2010.10130 | F |
| B10 | survey of trace inequalities, including characterisations of traces | matrices, operators | A. M. Bikchentaev, F. Kittaneh, M. S. Moslehian, Y. Seo, *Trace Inequalities: For Matrices and Hilbert Space Operators*, Forum for Interdisciplinary Mathematics, Springer, 2024, doi:10.1007/978-981-97-6520-1 | L |

### B′. Second round (not evaluated)

| # | statement (as far as seen) | algebra | reference | seen |
|---|---|---|---|---|
| B11 | Petz–Zemánek: "a number of equivalent conditions which characterize the trace among the linear functionals on the matrix algebra; some of these results are extended to more general operator algebras". **The conditions themselves were not visible.** No source confirms that they include Golden–Thompson. Bikchentaev (2011) answers "a question of J. Zemánek" about the Peierls–Bogoliubov inequality | M_n, operator algebras | Petz, Zemánek (1988), B1 | A |
| B12 | weighted monotonicity: Tr(w(A) f(A)) ≤ Tr(w(A) f(B)) for Hermitian A ≤ B and a weight w ≥ 0; the paper gives conditions that characterise the trace through this inequality | M_n | D. T. Hoa, O. E. Tikhonov, *Weighted trace inequalities of monotonicity*, Lobachevskii J. Math. 26 (2007), 63–67 | A |
| B13 | operator version, with the symmetric form τ(w(A)^{1/2} f(A) w(A)^{1/2}). It characterises traces among all weights on a von Neumann algebra, and among positive functionals on a C*-algebra, by operator inequalities. The pairs (w, f) were not seen, so the item was not evaluated: with w(t) = t and an operator monotone f, φ_W(A^{1/2}(f(B)−f(A))A^{1/2}) ≥ 0 holds for every positive φ, so no violation is possible for those choices | von Neumann / C*-algebras | D. T. Hoa, O. E. Tikhonov, *Weighted monotonicity inequalities for traces on operator algebras*, Math. Notes 88 (2010), 177–182, doi:10.1134/S0001434610070175 | A, S |
| B14 | a function is operator monotone (operator convex) if the corresponding inequalities hold for some normal state on B(H), dim H = ∞; operator convex functions with respect to a von Neumann algebra, by the type of its summands. This is a characterisation of *functions*, not of traces | B(H), von Neumann algebras | D. T. Hoa, O. E. Tikhonov, *To the theory of operator monotone and operator convex functions*, Izv. Vyssh. Uchebn. Zaved. Mat. 2010, no. 3, 9–14; Russian Math. (Iz. VUZ) 54:3 (2010), 7–11, doi:10.3103/S1066369X10030023 | A |
| B15 | characterisation of operator monotone functions by Powers–Størmer type inequalities (functions, not traces) | matrices | D. T. Hoa, H. Osaka, J. Tomiyama, *Characterization of operator monotone functions by Powers–Størmer type inequalities*, Linear Multilinear Algebra 63 (2015), no. 8, 1577–1589 (author list from a citing paper); see also arXiv:1207.5201 | S |
| B16 | survey chapter on matrix inequalities and characterisations of operator monotone functions | matrices | D. T. Hoa, H. Osaka, O. E. Tikhonov, *Matrix Inequalities and Characterizations of Operator Monotone Functions*, Springer book chapter (2022), doi:10.1007/978-3-031-02104-6_4 | L |
| B17 | a list of inequalities that characterise **central elements** of von Neumann algebras and C*-algebras: a positive a is central iff certain inequalities hold for all positive (normal) functionals. Key lemma: φ(sas) < φ(a) for a symmetry s produces two functionals separating a. This is the dual problem (fixed element, all functionals), not evaluated here | von Neumann / C*-algebras | A. A. Novikov, O. E. Tikhonov, *Characterization of central elements of operator algebras by inequalities*, Lobachevskii J. Math. 36 (2015), 208–210, doi:10.1134/S1995080215020146; arXiv:1502.01267 | A |
| B18 | related: measures on orthoideals and L₁-spaces associated with positive operators; L₁-space for a positive operator affiliated with a von Neumann algebra (no trace characterisation seen) | von Neumann algebras | A. A. Novikov, O. E. Tikhonov, Lobachevskii J. Math. 37 (2016), no. 4, 497–499; A. A. Novikov, *L₁-space for a positive operator affiliated with von Neumann algebra*, Positivity (2017), doi:10.1007/s11117-016-0422-4 | L |
| B19 | centrality characterised by local monotonicity of functions such as the exponential on C*-algebras (dual problem, as B17) | C*-algebras | D. Virosztek, *Connections between centrality and local monotonicity of certain functions on C\*-algebras*, J. Math. Anal. Appl. 453 (2017), 221–226; arXiv:1608.05409 | L |
| B20 | Hölder-type trace inequalities for Rényi means, and in the limit for the Log-Euclidean mean of positive definite matrices, characterise the canonical trace among positive linear functionals; also the Fiedler–Pták spectral geometric mean (forms not seen) | M_n | L. Molnár, *Rényi means, Hölder inequality and trace characterizations*, Positivity 30 (2026), Art. 70, doi:10.1007/s11117-026-01224-5 | A |
| B21 | commutativity of projections expressed by operator inequalities, applied to trace characterisation among positive normal functionals (part I of B6) | von Neumann algebras | A. M. Bikchentaev, *Commutativity of projectors and trace characterization on von Neumann algebras. I*, Russian Math. (Iz. VUZ) 53 (2009), no. 12, doi:10.3103/S1066369X0912010X (volume and issue inferred from the DOI; check) | L |
| B22 | a survey by Petz of trace inequalities, which lists Petz–Zemánek (1988) among trace characterisations | matrices | D. Petz, *A survey of certain trace inequalities*, Banach Center Publ. 30 (1994) | L |
| B23 | traces and subadditive measures on projections of JBW- and von Neumann algebras. Related, but not a characterisation of traces among functionals by an operator inequality | JBW/von Neumann | L. J. Bunce, J. Hamhalter, Proc. Amer. Math. Soc. (1995) (title from a citing listing) | L |

**On A. I. Stolyarov.** Beyond Stolyarov–Tikhonov–Sherstnev (2002), #2, the search found no further paper by Stolyarov
on trace characterisations.

**Conflicting attribution.** One search summary credits the |A+iB| characterisation (B5) to
"Bikchentaev and Tikhonov, 2021". The Springer listing names H. Alhasan and K. Fawwaz. The author should check.

## C. What we did not find

We found no paper that gives **quantitative** information of the following kind: how large the violation of a
trace-characterising inequality can be for a fixed non-tracial φ, as a function of a spectral quantity of its density.
Examples would be the Michelson contrast, or the condition number a/b of W. The closest work is B9. Its Fig. 1
plots the Gardner violation against the contrast of the density (Monte Carlo, unnormalised, n = 2, 3), and its
Sect. 7 conjectures contrast-dependent upper bounds, but it computes no supremum. A survey of the characterisations
themselves is A. M. Bikchentaev, M. S. Moslehian, *Characterizations of tracial functionals on matrix and operator
algebras*, Bull. Belg. Math. Soc. Simon Stevin 32 (2025), no. 3, 343–361 (not seen).

Our search was limited to abstracts and listings, so this statement is "not found", not "does not exist". The proofs
in the papers above may contain such estimates implicitly. For instance, the proof via "nearly parallel pure states"
in #6 must produce some violation size.

We found no paper that characterises traces by Golden–Thompson or Hölder for *general* operators. The search
summaries give only the projection versions, #14 and #15. We also found no trace characterisation by Kosaki, by
Hiai, or by Cho–Ito–Kittaneh, and none by a commutator inequality |φ(XY−YX)| ≤ …, beyond the commutator inequality
for Hermitian matrices in B3, whose form we did not see. Hiai–Petz ("The Golden–Thompson trace inequality
is complemented") appears in Petz's survey, but as a trace inequality, not as a characterisation.

## D. Proof details moved out of the paper subsection (for length)

**Proposition tcB(c)** (item 9, n = 2). Ranks 0 and 2 give no violation.

1. Set-up. Take rank-one P = uu\*, Q = vv\* and write u = c v + σ v_⊥. Normalise a + b = 2, and put β = v\*Wv,
   ρ = \|v\*Wv_⊥\| and γ = 2 − β.
2. The ratio (L − R)/N equals r/(1 + r), where, with c = cos θ,
   r = ρ sin 2θ / ((3β + γ)/2 + ((β − γ)/2) cos 2θ).
3. The maximum over θ is r = ρ/(2√β).
4. Use ρ² = β(2 − β) − (1 − δ²). The maximum over β is at β = √(1 − δ²), where r² = (1 − √(1 − δ²))/2 = sin²(ψ/2).

**Proposition tcC, item 7** (ALT, r = 2, q = 1).

1. Set-up. Work in the basis f₁ = e₊, f₂ = e₋. Take B = η f₁f₁\* + f₂f₂\* and A = uu\* with u = f₁ + sη f₂.
2. Limit. As η → 0, the ratio tends to W₁₁/(W₁₁ + 2sW₁₂ + s²W₂₂).
3. Maximum. At s = −W₁₂/W₂₂ this equals W₁₁W₂₂/(W₁₁W₂₂ − W₁₂²) = κ².

**Proposition tcC, item 10** (monotonicity, p = 2).

1. Set-up. Take A = uu\* and B = A + t vv\*, with u and v at angles α and π − α to e₁. Write y = −cos 2α and normalise
   a + b = 2.
2. Then Re (u\*v)(v\*Wu) = −y(δ − y) and u\*Wu = v\*Wv = 1 − δy.
3. With the optimal t = y(δ − y)/(1 − δy), the value is y²(δ − y)²/((1 − 2δy + y²)(1 − y²)).
4. Maximising over y gives m(δ).

**Item 14, leading order** (Hölder for projections).

1. Set-up. Take W = (1 + δD)/2 with D = diag(1, −1), and real unit vectors u, v at angles m ± ε/2. Put P = uu\* and
   Q = vv\*.
2. Expansion. φ(PQ) = (cos²ε + δ cos ε cos 2m)/2. With weights 1/p and 1/q,
   φ(P)^{1/p} φ(Q)^{1/q} = (1 + δ cos 2m + δ(1/q − 1/p) ε sin 2m)/2 + O(δ²ε²).
3. Hence L − R = (−ε² − δ(1 − 2/p) ε sin 2m)/2 + O(δ³) for ε = O(δ).
4. Maximise at ε = −δ(1 − 2/p) sin 2m / 2 and m = π/4. The result is V ≥ (1 − 2/p)² δ²/4 + O(δ³).
5. The search, the reduced SciPy model and the tests give equality to leading order, δ²/36 for p = 3.

## E. Full-text checks (2026-10-16)

- Bikchentaev, Math. Notes 89 (2011), no. 3–4, 461–471, doi:10.1134/S0001434611030175, Theorem 2: for φ ∈ M₊* and
  projections p, q the trace is characterised by (i) |φ(pq)| ≤ φ(p)^{1/α} φ(q)^{1/β}, α ≠ 2; (ii) φ(|pq|) ≤
  φ(p)^{1/α} φ(q)^{1/β}, α ≠ 2; (iii) φ(|pq|^{1/2}) ≤ φ(p)^{1/2} φ(q)^{1/2}; (iv) φ(e^{p+q}) ≤ φ(e^{p/2} e^q e^{p/2});
  (v) φ(e^{qpq}) ≤ φ(e^p); (vi) φ(Re pq) ≥ 0; (vii) φ(Re pq) ≤ φ(|pq|); (viii) φ(Im pq) ≥ 0; (ix) φ(f(p+q)) ≤
  φ(f(p)) + φ(f(q)) for concave f, f(0) = 0, f′(1) ≠ f(1). The remark after the theorem: α = 2 in (i) does not
  characterise traces. Items 14 and 15 are (i) with α = 3 and (iv).
- Bikchentaev, Sib. Math. J. 51 (2010), no. 6, 971–977, doi:10.1007/s11202-010-0096-2 (Russian original Sib. Mat.
  Zh. 51:6, 1228–1236), Theorem 3.4; announced without proofs in Russian Math. 53 (2009), no. 12, 68–71,
  doi:10.3103/S1066369X0912010X, Theorem 5: for φ ∈ M₊* the trace is characterised by φ(s_r(pqp)) = φ(s_r(qpq)),
  φ(p + q) ≤ φ(e + p ∧ q), φ(pqp) ≤ φ(q), φ(|p − q|) ≤ φ(p + q) (p, q projections). Not yet in the table.
- Bikchentaev, Lobachevskii J. Math. 32 (2011), no. 3, 175–179, doi:10.1134/S1995080211030061: Lemma 1 / Theorem 1
  (Peierls–Bogoliubov, φ(e^H) exp(φ(e^{H/2}Ke^{H/2})/φ(e^H)) ≤ φ(e^{H+K}), positive H, K) and Lemma 2 / Theorem 2
  (Araki–Lieb–Thirring, φ((A^{1/2}BA^{1/2})^{rp}) ≤ φ((A^{r/2}B^rA^{r/2})^p) for some r > 1, p > 0, positive A, B).
  Item 7 is Lemma 2 with r = 2, p = 1; item 16 is Lemma 1 (invariant under H ↦ H + c, K ↦ K + c).
- Tikhonov, Positivity 9 (2005), no. 2, 259–264, doi:10.1007/s11117-005-2711-1, Corollary 1: for a nonlinear
  concave f on R₊ with f(0) = 0, φ ∈ M₊* is tracial iff φ(f(A+B)) ≤ φ(f(A)) + φ(f(B)) for all positive A, B.
  Item 12 (f = t^{1/2}) is this corollary. Corollary 2: φ(AB + BA) ≥ 0 for all positive A, B.
- Bikchentaev, Sib. Math. J. 61 (2020), no. 2, 248–254, doi:10.1134/S0037446620020068, Theorem 3: for a positive
  φ on M_n with φ(I) = n, φ = tr iff (ii) det(P(exp A)) ≥ exp φ(A) for all pinchings P and A ≥ 0; (iii) det(A)^{1/n}
  ≤ φ(A)/n for A ≥ 0; (iv) per(A) ≤ φ(Aⁿ)/n for nonnegative A; (v) det(I + εA) = 1 + εφ(A) + o(ε); for faithful φ
  also (vi) det(exp A) ≤ exp φ(A) and (vii) φ(A^p)^{1/p} ≤ φ(A^q)^{1/q}, 0 < q < p. Theorem 4: Schur-type
  inequalities, e.g. Σ λ_t(A)² ≤ φ(A²), tr(A)² ≤ φ(√A)… for A ≥ 0. Item 11 is Theorem 3(iii).
