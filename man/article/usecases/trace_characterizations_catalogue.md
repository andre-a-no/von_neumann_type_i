# Inequalities that characterise the trace: a catalogue

Compiled for the use case "quantitative failure of trace characterisations" (`trace_characterizations.tex`,
`scripts/trace_characterizations.py`).

**How this list was made.** The only source was web search. We could not open mathnet.ru, Springer, ScienceDirect,
Cambridge Core or arXiv full texts from this machine. Every entry below rests on an **abstract, a publisher or
database listing, or a secondary description** (another paper's summary). None comes from the full text. The author
should check every statement against the paper.

The column "seen" says what we saw:

- **A**: the abstract on the publisher or database page, as shown by the search engine;
- **L**: a bibliographic listing only;
- **S**: a secondary description, i.e. how another paper or listing summarises the result.

"Form used" is the exact inequality evaluated in `scripts/trace_characterizations.py`. Where the source does not fix
the form, we say which reading we chose.

Throughout, φ is a positive linear functional. On M_n it has the form φ = Tr(W ·) with a density W ≥ 0. "Trace" means
φ(XY) = φ(YX); on M_n this is equivalent to φ = c·Tr.

## A. Inequalities evaluated numerically

| # | inequality (form used) | operators | algebra | statement of the characterisation | reference | seen |
|---|---|---|---|---|---|---|
| 1 | \|φ(X)\| ≤ φ(\|X\|) | arbitrary X | unital C*-algebra (general) | a positive functional is a trace iff the inequality holds for all X | L. T. Gardner, *An inequality characterizes the trace*, Canad. J. Math. 31 (1979), no. 6, 1322–1328, doi:10.4153/CJM-1979-109-9 | A |
| 2 | φ(\|X+Y\|) ≤ φ(\|X\|) + φ(\|Y\|) | self-adjoint X, Y (we optimise over self-adjoint ones) | von Neumann algebra (normal semifinite weights) | a normal semifinite weight satisfying the triangle inequality for all self-adjoint a₁, a₂ is a trace. Gardner (1979) states the same for maps on C*-algebras with the triangle inequality for all x, y (abstract only) | A. I. Stolyarov, O. E. Tikhonov, A. N. Sherstnev, *Characterization of normal traces on von Neumann algebras by inequalities for the modulus*, Mat. Zametki 72 (2002), no. 3, 448–454; Math. Notes 72 (2002), 411–416, doi:10.1023/A:1020559623287; and Gardner (1979) | A |
| 3 | φ(A + A\*) ≤ 2φ(\|A\*\|) | arbitrary A | von Neumann algebra | stated as a corollary: tracial functionals are characterised by this inequality | A. M. Bikchentaev, S. A. Abed, *Projections and traces on von Neumann algebras*, Lobachevskii J. Math. 40 (2019), 1260–1267, doi:10.1134/S1995080219090051 | S |
| 4 | \|φ(A²)\| ≤ φ(A\*A) (also \|φ(Re A²)\| ≤ φ(A\*A)) | arbitrary A | von Neumann algebra, positive normal φ | φ is tracial ⟺ \|φ(Re(A²))\| ≤ φ(\|A\|²) ∀A ⟺ \|φ(A²)\| ≤ φ(\|A\|²) ∀A | A. M. Bikchentaev, M. S. Moslehian, *Trace inequalities and characterizations of tracial functionals in operator algebras*, Positivity (2026), doi:10.1007/s11117-026-01180-0 (volume and pages not seen) | A |
| 5 | φ(\|P−R\|) ≤ φ(\|P−Q\|) + φ(\|Q−R\|) | projections | von Neumann algebra | ρ_φ(P,Q) = φ(\|P−Q\|) (or d_φ(P,Q) = φ(P∨Q − P∧Q)) satisfies the triangle inequality iff φ is tracial. That ρ_τ is a metric needs a faithful trace | A. M. Bikchentaev, *Metrics on projections of the von Neumann algebra associated with tracial functionals*, Sib. Math. J. 60 (2019), 952–956, doi:10.1134/S003744661906003X | A |
| 6 | φ(A ♮ B) ≤ (φ(A) φ(B))^{1/2}, with A ♮ B = (A⁻¹#B)^{1/2} A (A⁻¹#B)^{1/2} | positive definite A, B | M_n | φ is a positive multiple of Tr iff the inequality holds for all positive definite A, B. The proof uses "nearly parallel pure states". The paper also characterises the trace by φ(A♮B) ≤ φ((A+B)/2) | A. M. Bikchentaev, T. H. Dinh, A. V. Le, M. S. Moslehian, *Spectral geometric mean and trace characterizations*, Positivity 30 (2026), Art. 46, doi:10.1007/s11117-026-01202-x; arXiv:2605.18888 | A |
| 7 | φ((B^{1/2}AB^{1/2})^{rq}) ≤ φ((B^{r/2}A^rB^{r/2})^q). We use **r = 2, q = 1**: φ((B^{1/2}AB^{1/2})²) ≤ φ(BA²B) | positive A, B | C*-algebra | "either the Peierls–Bogoliubov or the Araki–Lieb–Thirring inequality characterises tracial functionals among positive functionals"; the exponents for which this is proved were not seen | A. M. Bikchentaev, *The Peierls–Bogoliubov inequality in C\*-algebras and characterization of tracial functionals*, Lobachevskii J. Math. 32 (2011), no. 3, 175–179, doi:10.1134/S1995080211030061 | A |
| 8 | φ(\|AB\*\|) ≤ ½φ(\|A\|²) + ½φ(\|B\|²) (Young, p = q = 2, in the form of Ando's trace inequality) | arbitrary A, B | M_n | Young's inequality for φ, for conjugate exponents, characterises the trace. **The exact form (which of \|AB\*\|, \|A\*B\|, which normalisation) was not seen**; our form is a guess modelled on Ando's trace form | A. M. Bikchentaev, O. E. Tikhonov, *Characterization of the trace by Young's inequality*, J. Inequal. Pure Appl. Math. (JIPAM) 6 (2005), no. 2, Art. 49 | S |
| 9 | φ(PQ + QP) ≤ 2φ((QPQ)^p). We use **p = 1**: φ(PQ+QP) ≤ 2φ(QPQ) | projections | von Neumann algebra, positive normal φ | φ is tracial iff the inequality holds for all projections P, Q and some p = p(P,Q) ∈ (0, 1]. The same abstract lists φ(PQP) ≤ φ(P)^{1/p} φ(Q)^{1/q}, 1/p + 1/q = 1, p ≠ 2, as equivalent; we did not use it because the role of p ≠ 2 is unclear to us from the snippet | Bikchentaev and Abed (2019), as in #3 | S |
| 10 | 0 ≤ A ≤ B ⟹ φ(A^p) ≤ φ(B^p). We use **p = 2** | positive A ≤ B | M_n | if the implication holds for all A ≤ B (fixed p > 1), then φ is a non-negative multiple of Tr | A. M. Bikchentaev, O. E. Tikhonov, *Characterization of the trace by monotonicity inequalities*, Linear Algebra Appl. 422 (2007), no. 1, 274–278, doi:10.1016/j.laa.2006.10.005 | A |
| 11 | (det A)^{1/n} ≤ φ(A)/n with φ(1) = n (an AM–GM analogue) | positive A | M_n | if "some analog of the classical inequalities for the determinant and trace (or permanent and trace)" holds for φ with φ(1) = n, then φ = Tr. **Which determinant inequality the paper uses was not seen**; AM–GM is our choice. This item characterises the *canonical* trace, not traces in general | A. M. Bikchentaev, *Inequalities for determinants and characterization of the trace*, Sib. Math. J. 61 (2020), 248–254, doi:10.1134/S0037446620020068 | A |

## B. Further characterisations (not evaluated: the form was not visible, or the statement is not an inequality)

| # | statement (as far as seen) | algebra | reference | seen |
|---|---|---|---|---|
| 12 | "Characterizations of the trace": the content was not visible | M_n | D. Petz, J. Zemánek, *Characterizations of the trace*, Linear Algebra Appl. 111 (1988), 43–52, doi:10.1016/0024-3795(88)90050-X | L |
| 13 | a positive functional φ on M_n is a positive multiple of Tr iff φ(A) = φ(\|A\|) implies A ≥ 0. For φ(1) = n, a Yang-type inequality φ(A^{1/2}BA^{1/2})^{1/2} ≤ φ(A+B)/2 (A, B ≥ 0) singles out Tr; the placement of the square root was not clear to us | M_n | M. S. Moslehian, A. M. Bikchentaev, *Characterizations of the canonical trace on full matrix algebras*, J. Math. Anal. Appl. 552 (2025), no. 2, 129764 | A |
| 14 | characterisations of Tr on M_n by a commutator inequality for Hermitian matrices and by a generalised AM–GM inequality for all matrices and conjugate exponents (forms not seen) | M_n | M. S. Moslehian, *More characterizations of the matrix trace*, Rend. Circ. Mat. Palermo (2) 74 (2025), Art. 130, doi:10.1007/s12215-025-01253-1 | A |
| 15 | inequalities for the spread λ_max − λ_min characterise Tr among positive φ with φ(1) = n, for n ≥ 3; three further characterisations hold for n ≥ 2 | M_n | M. Dalloul, *Characterization of tracial functionals by inequalities for the spread of a matrix*, Lobachevskii J. Math. 45 (2024), 2478–2484, doi:10.1134/S1995080224602777 | A |
| 16 | an inequality characterising tracial functionals among positive normal functionals on a von Neumann algebra, "strengthening Gardner". One secondary snippet gives it as φ(A) ≤ φ(\|A + iB\|) (A, B self-adjoint). If so, it is #3 with Z = A + iB, and its envelope is κ − 1 | von Neumann algebra | H. Alhasan, K. Fawwaz, *Characterization of tracial functionals on von Neumann algebras*, Lobachevskii J. Math. 42 (2021), 2273–2279, doi:10.1134/S1995080221100024 | A, S |
| 17 | operator inequalities that are necessary and sufficient for two projections to commute, used to characterise traces among positive normal functionals (part III: arbitrary weights, and positive functionals on C*-algebras) | von Neumann algebra | A. M. Bikchentaev, *Commutativity of projections and characterization of traces on von Neumann algebras*, Sib. Math. J. 51 (2010), 971–977, doi:10.1007/s11202-010-0096-2; *… II*, Math. Notes 89 (2011), 461–471, doi:10.1134/S0001434611030175; *… III*, Int. J. Theor. Phys. 54 (2015), doi:10.1007/s10773-015-2639-6 | A |
| 18 | a new inequality for a trace on a unital C*-algebra that characterises traces among positive functionals, and a commutativity criterion | unital C*-algebra | A. M. Bikchentaev, *Inequality for a trace on a unital C\*-algebra*, Math. Notes 99 (2016), 487–491, doi:10.1134/S0001434616030214 | A |
| 19 | tracial states among all states on a Rickart C*-algebra characterised via polar decomposition; a triangle inequality for Hermitian elements and traces | Rickart C*-algebra | A. M. Bikchentaev, *Trace inequalities for Rickart C\*-algebras*, Positivity 25 (2021), 1943–1957, doi:10.1007/s11117-021-00852-3 | A |
| 20 | *related, quantitative*: the Michelson contrast of positive operators in W*-algebras, and "how the inequalities characterizing traces interplay with the Michelson contrasts of operator **variables**". This concerns the contrast of the variables, not of the density of φ | countably decomposable W*-algebras | S. A. Abed, I. A. Nikolaeva, A. A. Novikov, *Generalisation of Michelson contrast for operators and its properties*, Lobachevskii J. Math. 45 (2024), 3835–3848, doi:10.1134/S1995080224603400; arXiv:2010.10130 | A |
| 21 | survey of trace inequalities, including characterisations of traces | matrices, operators | A. M. Bikchentaev, F. Kittaneh, M. S. Moslehian, Y. Seo, *Trace Inequalities: For Matrices and Hilbert Space Operators*, Forum for Interdisciplinary Mathematics, Springer, 2024, doi:10.1007/978-981-97-6520-1 | L |

## C. What we did not find

We found no paper that gives **quantitative** information of the following kind: how large the violation of a
trace-characterising inequality can be for a fixed non-tracial φ, as a function of a spectral quantity of its density.
Examples would be the Michelson contrast, or the condition number a/b of W. The closest work is #20, which asks
about the contrast of the operator variables.

Our search was limited to abstracts and listings, so this statement is "not found", not "does not exist". The proofs
in the papers above may contain such estimates implicitly. For instance, the proof via "nearly parallel pure states"
in #6 must produce some violation size.

We found no reference for a characterisation by a Golden–Thompson-type inequality φ(e^{A+B}) ≤ φ(e^A e^B), nor by
Hölder's inequality for φ, so neither is listed.
