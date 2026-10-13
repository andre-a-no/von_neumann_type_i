# Inventory of unproved statements in the English course (docs/tutorial)

Scope: `README.md`, `00_overview.md` … `13_jordan_algebras.md` (English originals only), read in full.
Criterion: the README promise (l.12–16) "Nothing is quoted without proof except standard linear algebra,
the fundamental theorem of algebra and elementary facts of analysis and probability (limits, convexity,
the Gaussian moment generating function)". Everything else stated or used without a complete proof
in the course is listed. Library claims count when they are mathematical.

Action codes: **P** = prove in place (route + length); **A** = prove in an appendix (length);
**N** = replace by a numerical observation computed by the chapter's own code, clearly labelled
(only for numbers/asymptotics not used as premises); **R** = remove or rephrase so that nothing is
claimed (inessential). **PHYS** = physical/empirical claim in a *Physical meaning* box or in prose:
not provable mathematics; action R, or keep under an explicit README disclaimer (see plan).
Lengths are lines of Markdown in the current style (~110 characters per line).
Severity: **[premise]** = used in a proof or as the target of an `assert`; **[stated]** = stated
but not used; **[README]** = listed among the "main results" of the README.

---

## README.md

| ID | line | statement | what is missing | action |
|---|---|---|---|---|
| R.1 | 12–16 | allowed-without-proof list | The course actually uses Lebesgue measure, Fubini–Tonelli, dominated convergence, change of variables (incl. polar/spherical coordinates, Jacobian of orthogonal maps), Dynkin π–λ, Borel measures on metric spaces, term-by-term differentiation of power series, Riemann integrals of continuous matrix functions, inverse function theorem (Ch.11) — none on the list. | Amend list (2–3 lines) to "standard measure and integration theory (Lebesgue measure, Fubini–Tonelli, dominated convergence, change of variables) and elementary calculus (power series, Riemann integral)", and collect the derived lemmas in Appendix A. Proving measure theory itself is not realistic (>400 lines). |
| R.2 | 39 | Ch.7 main result "Page's formula" | Theorem 7.13 is stated without proof. | Change table entry (1 line) unless Page is proved (see 7.6). |
| R.3 | 45 | Ch.13 main result "the Jordan–von Neumann–Wigner classification" | Theorem 13.7 is stated without proof. | Change table entry (1 line) unless proved (see 13.7). |
| R.4 | 17–18 | *Physical meaning* boxes | They contain empirical physics claims (MBL, BCS, CFT, Page curve, …) that cannot be proved. | Add one sentence: boxes are interpretation, cite physics literature, and are never used in proofs or checks; then remove from them every *mathematical* claim that is not proved (listed below). |

## Chapter 0 — `00_overview.md`

| ID | line | statement | what is missing | action |
|---|---|---|---|---|
| 0.1 | 118–120 | $\prod_cU(k_c)$ "a null set inside $U(D)$, never reached by sampling $U(D)$" | Proved only in §5.1 (l.110–112); no pointer. | P: add "(Section 5.1)". 0–1 line. |
| 0.2 | 110–116 | spectral decomposition "costs $\sum_ck_c^3$ operations against $D^3$" | Operation count of dense eigensolvers is an algorithmic fact, not proved. | R/relabel: "dense eigensolvers need $O(k^3)$ operations for a $k\times k$ matrix". 1 line. Low priority. |

## Chapter 1 — `01_matrix_algebras.md`

| ID | line | statement | what is missing | action |
|---|---|---|---|---|
| 1.1 | 236 | code comment: $\lvert X\rvert$ invertible "almost surely, for Ginibre" | Forward fact (Theorem 5.3: $\det Z\neq0$ a.s.). | P: pointer "(Theorem 5.3)". 0 lines. |
| 1.2 | 257–262 | differentiation under the expectation "by dominated convergence" (Isserlis) | DCT is measure theory, not on the README list (the MGF is). | Covered by R.1/Appendix A (quote DCT). 0 lines here. |
| 1.3 | 278–279 | "the Catalan number 2 of the semicircle law" | Semicircle law never proved. | R: "the two planar pairings give the limit $2$ (a Catalan number)". 0 lines. |

## Chapter 2 — `02_von_neumann_algebras.md`

| ID | line | statement | what is missing | action |
|---|---|---|---|---|
| 2.1 | 248–249 | symmetric subspace = "total spin $\tfrac32$" | Proved only in Ch.4 (Prop. 4.13). | P: pointer. 0 lines. |
| 2.2 | 256–258 | $\sum_ck_c^2=\dim\pi(S_3)''=5$ "computed by the code below by exact linear algebra" — **[premise]** of the proof of the decomposition | A floating-point SVD rank is a numerical computation, not a proof. The hand count gives only $\dim\pi(S_3)'=20$. | P (3–5 lines): case analysis from $\sum k_cm_c=8$, $\sum m_c^2=20$ alone: $m\in\{(4,2),(4,1,1,1,1),(3,3,1,1),(2,2,2,2,2)\}$; all but $(4,2)$ force all $k_c=1$ or are infeasible, contradicting non-commutativity; then $k=(1,2)$. (Or forward pointer to the exact character count of §4.2.) |
| 2.3 | 271–273 | the 2-dim irreducible $\rho_2$ "is the standard representation" | Uniqueness of the 2-dim irrep of $S_3$ not shown here. | P: say "an irreducible 2-dimensional representation; §4.2 identifies it with the standard one by characters". 1 line. |
| 2.4 | 316–318 | a random element of $\pi(S_3)'$ "has (almost surely) 4 simple eigenvalues … and 2 doubly degenerate ones" | Genericity (Prop. 4.7, forward) + the twirled GUE having a density on $\mathcal A'_h$ (see 4.4). | P: pointer to Prop. 4.7 and Appendix A.2. 1 line. |
| 2.5 | 411–414 (box) | $X\mapsto\sum_cP_cXP_c$ "is the conditional expectation onto $\mathcal M$" | Not verified (positivity, bimodule). | P: 2 lines (as in Prop. 8.8 proof: $E(XB)=E(X)B$, $E(B)=B$, $P_cXP_c\ge0$). |

## Chapter 3 — `03_traces_and_states.md`

| ID | line | statement | what is missing | action |
|---|---|---|---|---|
| 3.1 | 307–309 | "The companion paper shows numerically that Monte Carlo sampling falls short …, increasingly so as $k$ grows" | Numerical claim from outside the course. | R or N (the code below already shows one pair). 0 lines. |
| 3.2 | 310–312 | for Haar $U$ in high dimension $\lvert\mathrm{Tr}(XUY)\rvert$ "is typically much smaller than $\mathrm{Tr}(XY)$" | Unproved here; it *is* proved later (Prop. 11.11, Cor. 11.12). | P: forward pointer to Cor. 11.12. 1 line. |

Otherwise Ch.3 is complete (Klein, Theorem 3.7, Prop. 3.9 equality case use only Cauchy–Schwarz equality, standard).

## Chapter 4 — `04_symmetries.md`

| ID | line | statement | what is missing | action |
|---|---|---|---|---|
| 4.1 | 230–233, 250 | zero set of a non-zero real polynomial is Lebesgue-null (proved "by Tonelli") | Uses Lebesgue measure + Tonelli (not on the list). Same lemma re-proved in 5.3 (l.75–79) and cited in Ch.12 Ex.4. | A: Appendix A.1 (8–10 lines) once; quote Tonelli via R.1. |
| 4.2 | 244–245 | on the set where (a) holds "the spectral projections are polynomials in $K$ whose coefficients are continuous functions of the eigenvalues" — **[premise]** for measurability in Tonelli | Continuity of the ordered eigenvalues of a Hermitian matrix in its entries is never proved (also used in 5.15). | A: Appendix A.3 — Courant–Fischer (6 lines) ⇒ Weyl $\lvert\lambda_i(A)-\lambda_i(B)\rvert\le\lVert A-B\rVert$ (4 lines) ⇒ continuity of Lagrange-polynomial spectral projections on the simple-spectrum set (3 lines). ~13 lines. |
| 4.3 | 265–267 | "images of Gaussian laws under linear surjections onto $\mathcal A'$ and $\mathcal A'_h$, hence with densities" — **[premise]** of Prop. 4.7 in the algorithm | Push-forward of a non-degenerate Gaussian under a linear surjection has a density: needs change of variables. | A: Appendix A.2 (split $\mathbb R^N=\ker L^\perp\oplus\ker L$, Fubini, linear change of variables). ~8 lines. |
| 4.4 | 275–277 | floating point: "close eigenvalues of different copies are kept together and split by the compression of a fresh random element of the commutant (module docstring)" | Library claim (decompose docstring: "compression $Q^*YQ$ … is scalar iff the restriction is irreducible") not proved. Output correctness is covered by the certificate (Prop. 4.9), so the claim is not needed. | R ("the result is certified, Prop. 4.9"), or P ~10 lines: for an $\mathcal A$-invariant $W$, $Q^*\mathcal A'Q$ is the commutant of $\mathcal A\vert_W$ (structure theorem restricted to $W$) + genericity as in Prop. 4.7. |
| 4.5 | 320–322 | library reduces the generator problem "with a random element of $\mathcal A$ to $N=\sum_ck_cm_c^2$ unknowns" | Genericity for a random element of $\mathcal A$ (not $\mathcal A'$) not proved. | P: 5 lines, mirror of Prop. 4.7(a) with roles of $k_c,m_c$ exchanged. |
| 4.6 | 416–417 | "every such $V\in SU(2)$ is some $e^{-i\theta\,\mathbf n\cdot\mathbf s}$", and $e^{-i\theta\mathbf n\cdot\mathbf S}=(e^{-i\theta\mathbf n\cdot\mathbf s})^{\otimes L}$ | Neither proved; the second needs $e^{X+Y}=e^Xe^Y$ for commuting $X,Y$ and $e^{A\otimes1}=e^A\otimes1$ (not proved anywhere; Lemma 8.1 only gives $e^{(s+t)B}$). | P: 4 lines (spectral theorem: $V=W\mathrm{diag}(e^{i\varphi},e^{-i\varphi})W^*$, $\log$ traceless Hermitian $=$ real combination of Paulis; binomial series for commuting exponentials). |
| 4.7 | 503–512 | **Theorem 4.14(b) Schur–Weyl duality "stated without proof for general $L$"**; for $L\le5$ "we check it" (floating-point ranks) | Full proof absent; the $L\le5$ check is numerical. **[premise]** for l.528–529 (multiplicity spaces carry irreducible $S_L$-representations), Ex.1 solution l.637–638 ("35 is also $\dim\{S_a\}''$"), and the general Weingarten Theorem 5.12. | A: Appendix B (shared with 5.8). Route: (i) polarisation: permutation-invariant elements of $M_d^{\otimes L}$ are spanned by $X^{\otimes L}$; (ii) $\{X^{\otimes L}\}$ and $\{\sum_iX^{(i)}\}$ (hence $\{S_a\}$, or $U(d)^{\otimes L}$) generate the same algebra (exponential/derivative + polynomial-in-$t$ argument); (iii) double commutant theorem. ~35–45 lines. |
| 4.8 | 526 | printed dimensions "$2,5,14,42$, the Catalan numbers $\binom{2L}L/(L+1)$" | General identity $\sum_Sm_S^2=C_L$ unproved; only $L\le5$ printed. | N: "for $L=2,\dots,5$ these are the Catalan numbers" (0 lines). Optional P 8 lines: $\dim\{S_a\}'=$ multiplicity of spin 0 in $2L$ spins $=m_0(2L)=C_L$. |
| 4.9 | 556–558 | `su2` builds an orthonormal basis of highest-weight vectors "by Clebsch–Gordan coupling one site at a time" | Library correctness; explicit CG coefficients for $j\otimes\frac12$ are never derived (Ex.2 gives only the decomposition). Checked numerically. | P ~12 lines (explicit $j\otimes\frac12$ coefficients from $S_-$ and orthogonality) or N ("implementation detail, verified by the code"). |
| 4.10 | 591–594 | block sizes, cost factor 22, ground energy $-3.651093$, gap level $-3.128419$ | Numerical outputs presented as facts. | N: label "computed by the code above". 0 lines. |
| 4.11 | 598–601 (box) | Lieb–Mattis theorem for bipartite lattices "stated without proof"; "the ring of eight spins above is an example" | Unproved (general lattices). For chains, Theorem 9.9 proves it. | R: replace by forward pointer to Theorem 9.9. (A alternative: generalise Thm 9.7/9.9 to connected bipartite graphs ~40 lines.) |
| 4.12 | 596–598 (box) | Heisenberg model describes exchange in magnetic insulators | PHYS. | PHYS. |

## Chapter 5 — `05_random_matrices.md`

| ID | line | statement | what is missing | action |
|---|---|---|---|---|
| 5.1 | 38–46 | Lemma 5.1 (uniqueness of Haar) uses Fubini, dominated convergence, Dynkin π–λ, Borel measures on metric spaces | Measure theory beyond the list. | Quote via R.1 / Appendix A (list of quoted facts). 0 lines here. |
| 5.2 | 52–54 | "Every compact group has a Haar measure (Folland; stated here without proof)" | Unproved; **not used** (all groups used — $U(k)$, $U(\mathcal M)$, $SU(k)$, $O(k)$, $U(1)$ — get Haar measure by construction). | R. 0 lines. (A would be ~80–120 lines via Kakutani/Riesz or averaging; not worth it.) |
| 5.3 | 75–79 | zero set of a polynomial is null "(… conclude by Fubini)" | Duplicate of 4.1. | A: Appendix A.1; replace by pointer. |
| 5.4 | 81–84 | orthogonal map of $\mathbb R^{2k^2}$ has "Jacobian 1" and preserves Lebesgue measure | Change of variables. | Quote via R.1 / Appendix A. |
| 5.5 | 93 | "The same proof with real Gaussian matrices gives the Haar measure of $O(k)$" — **[premise]** for Ch.11 Ex.2 and for real algebras | Not written out (sign correction, Lemma 5.1 for $O(k)$). | P: state as Theorem 5.3′ with a 3-line proof. |
| 5.6 | 94–97, 177 | LAPACK's sign convention "introduces a bias"; $\mathbb E\,\mathrm{Re\,Tr}\,Q\approx-0.84$ for $k=2$ | Software behaviour + a printed number. | N: label as numerical observation of the code. 0 lines. |
| 5.7 | 128–129 | `'cse'` is Dyson's circular symplectic ensemble | Library claim; CSE never defined or checked. | R, or P 5 lines (define $S=W^DW$ with quaternion dual, check unitarity/self-duality). |
| 5.8 | 306–321 | **Theorem 5.12 general Weingarten formula "stated without proof"** (incl. class-function property, Gram inverse for $n\le k$, pseudo-inverse for $n>k$, vanishing of unbalanced moments) — **[premise]** of Ex.2 and of the check $\mathbb E\lvert U_{11}\rvert^6=1/10$ (l.665) | No proof; rests on Schur–Weyl for $U(k)$ on $(\mathbb C^k)^{\otimes n}$. | A: Appendix B (after 4.7): twirl lies in span of permutations (5 lines), linear independence of permutation operators iff $n\le k$ (10 lines), Gram/pseudo-inverse (12 lines), class function (3 lines), phase argument for unbalanced moments (2 lines). ~40 lines on top of Schur–Weyl. |
| 5.9 | 378 | density of GOE/GUE "proportional to $e^{-\mathrm{Tr}H^2/2}$" | Not computed. | P: 2 lines ($\mathrm{Tr}H^2=\sum_ih_{ii}^2+2\sum_{i<j}\lvert h_{ij}\rvert^2$). |
| 5.10 | 391–396 | polar and spherical coordinates in Prop. 5.13 | Change of variables. | Quote via R.1. |
| 5.11 | 402–407 | joint eigenvalue density $\propto\prod_{i<j}\lvert\lambda_i-\lambda_j\rvert^\beta e^{-\sum\lambda_i^2/2}$ for every $n$; Gaudin–Mehta limit "(stated without proof)" | Unproved. The joint density is a **[premise]** at l.486 ("For the GOE and the GUE this [bounded continuous level density] follows from the joint density quoted above") and l.516. | Gaudin–Mehta: R. Joint density: R (make Theorem 5.15 explicitly conditional on the density hypothesis, delete l.486 sentence, label the GOE experiment numerical), or A ~100–150 lines (Weyl integration formula: $H=U\Lambda U^*$, Jacobian $\prod\lvert\Delta\rvert^\beta$, needs manifold change of variables). Recommend R. |
| 5.12 | 409–415 (box) | Wigner's hypothesis, von Neumann–Wigner level repulsion, BGS and Berry–Tabor conjectures | PHYS (already called conjectures). | PHYS. |
| 5.13 | 439–442 | uniform levels: spacings $\sim(e_1,\dots,e_{n+1})/\sum e_j$ (Feller, "stated without proof") — **[premise]** of the Monte-Carlo assert l.468–469 (uniform levels vs $2\ln2-1$) | Unproved. | Best: R + change the code to Poisson levels (cumulative sums of i.i.d. exponentials), to which Prop. 5.14 applies directly (2 code lines, 2 text lines). Alternative P ~12 lines (order statistics of uniforms via density $n!$ on the simplex + change of variables). |
| 5.14 | 443–445 | $\langle r\rangle_{GOE}=0.5307$, $\langle r\rangle_{GUE}=0.5996$ (literature numerics), $3\times3$ surmise values $4-2\sqrt3$, $2\sqrt3/\pi-\tfrac12$ "stated without proof" — 0.5307 is an **assert target** at l.475 (5σ test for GOE(64)), l.550, l.682 (COE), Ch.9 l.578 | Literature constants used as test targets; $3\times3$ values unproved. | N: compute the reference value with the course's own code (e.g. GOE $n=200$, enough samples; print "numerical estimate $\langle r\rangle_{GOE}\approx0.531$") and compare against it with a stated tolerance; drop GUE value or compute it likewise. $3\times3$ values: R (A would need the $n=3$ joint density + a 2-D integral, ~50 lines). ~10 code + 5 text lines. |
| 5.15 | 501 | "The ordered eigenvalues are continuous functions of the matrix" — **[premise]** of Theorem 5.15 | Unproved (same as 4.2). | A: Appendix A.3. |
| 5.16 | 502–511 | Fubini/conditioning, dominated convergence in Theorem 5.15 | Measure theory. | Quote via R.1. |
| 5.17 | 514–517 | for an $n\times n$ GOE block $P(\delta<\varepsilon)$ is of order $\varepsilon^2$ "(from the joint density; not proved here)" | Unproved. | R (keep the proven $2\times2$ statement), or N. |
| 5.18 | 519–532 | "Many sectors (a derivation, not a proof)": $p(s)=\mathcal G''(s)$ (Mehta), superposition multiplies gap probabilities, $p(0)=1-1/C$, Poisson limit (Rosenzweig–Porter, Berry–Robnik, Grigelionis) | Explicitly unproved heuristic. | R (the experiment that follows already shows the effect numerically: label it N). Rigorous version needs stationary point-process theory (>100 lines). |
| 5.19 | 555–556, 575–576 | $\langle r\rangle\approx0.53$; merged $0.42$, $0.39$; sector $N=5$: $0.52$ "a finite-size value slightly below 0.5307" | Numerical; comparison with the literature constant. | N: label as output of the code; compare with the course-computed constant (5.14). |
| 5.20 | 558–562 | "for weak disorder the levels inside each sector have GOE statistics (numerically; Pal–Huse)" | Literature numerics. | N: "the code finds 0.52". |
| 5.21 | 578–583 (box) | levels with different quantum numbers cross freely; unresolved symmetry mistaken for integrability | Mathematical half = Theorem 5.15; rest PHYS. | PHYS. |
| 5.22 | 669–673, 682 (Ex.3) | "local statistics of the COE coincide with those of the GOE for large $k$ (Dyson; Mehta; stated without proof)" — **assert** `mc(r, 0.5307)` for COE(32) | Unproved; used as 5σ test target. | N: compare with the course-computed GOE constant with a stated tolerance, labelled numerical; R the Dyson sentence. |
| 5.23 | 669–671 (Ex.3) | rotation-invariant mean eigenphase density "is therefore uniform" | Needs uniqueness of rotation-invariant probability on the circle. | P: 1 line, Lemma 5.1 for $U(1)$ (Haar exists by Theorem 5.3 with $k=1$). |

## Chapter 6 — `06_quantum_channels.md`

| ID | line | statement | what is missing | action |
|---|---|---|---|---|
| 6.1 | 188–192 | general Stinespring theorem for unital C\*-algebras "(stated without proof)" | Unproved; not used. | R. |
| 6.2 | 460–464 (box) | $1/T_2=1/(2T_1)+1/T_\varphi$ | Follows from Prop. 8.14 + Ch.8 Ex.4 once the two generators are added (they commute on coherences), not shown. | P 3 lines, or PHYS. |

Otherwise Ch.6 is complete (Choi, Kraus, Stinespring in finite dimensions, unitary freedom, duals, Theorems 6.10, 6.11, 6.15).

## Chapter 7 — `07_tensor_products.md`

| ID | line | statement | what is missing | action |
|---|---|---|---|---|
| 7.1 | 237–244 | $m^2=\lVert\Lambda\rVert^2=\dim(M_k\otimes M_m)/\dim M_k$ "is its Jones index $[\mathcal M:\mathcal N]$" | Identification with Jones' index (coupling constant, subfactors) unproved. | R: present $\lVert\Lambda\rVert^2$ as a definition with a historical pointer, no claim. 0 lines. |
| 7.2 | 378–381 (box) | BCS: only parity conserved; superposition of 0 and 2 particles | PHYS. | PHYS (or P 2 lines if the BCS state $\prod(u+vc^*c^*)\lvert0\rangle$ is defined). |
| 7.3 | 462–465 (box) | Wiseman–Vaccaro resource claim; number entropy measured by Lukin et al. | PHYS. | PHYS. |
| 7.4 | 471–472 | uniform law is "the only invariant law on the sphere, a fact we do not use" | Unproved, unused. | R. |
| 7.5 | 471, 476–483 | unitary invariance of Lebesgue measure; polar coordinates in $\mathbb R^{2N}$ with surface measure (Lemma 7.9) | Change of variables / polar coordinates. | Quote via R.1 (Appendix A list). |
| 7.6 | 520–530 | **Theorem 7.13 Page's formula "stated without proof"**; asymptotic $\approx\log a-\frac a{2b}$; — **assert** l.570 (5σ test of mean entropy against Page's value) **[README]** | Unproved. | Either A (~200–300 lines: joint eigenvalue density of the induced/Wishart–Laguerre ensemble via Weyl integration formula ~100, plus Sen's or Sánchez-Ruiz' Laguerre-polynomial computation ~120; the asymptotic then 5 lines), or N: keep the proven bounds (Cor. 7.12) and Lubkin's purity, replace Page's value by a numerical observation ("the code finds $\mathbb ES\approx\dots$; Page conjectured, and Foong–Kanno, Sánchez-Ruiz and Sen proved, a closed formula; not proved here" → under strict rule: drop the formula), and change the assert to the proven bounds. Recommend N/R + README update. |
| 7.7 | 532–537 | induced measure with several blocks = restriction of a Haar pure state of $H_1\otimes\mathbb C^r$ (sketch "apply Theorem 7.6") | Only sketched. | P: 4 lines. |
| 7.8 | 579–582 (box) | Page curve, black-hole information, thermalisation | PHYS. | PHYS. |

## Chapter 8 — `08_dynamics.md`

| ID | line | statement | what is missing | action |
|---|---|---|---|---|
| 8.1 | 212–217 | **Theorem 8.7 (GKSL converse) "stated without proof"**; Christensen–Evans for general von Neumann algebras | Unproved. | A: Appendix D, finite-dimensional proof ~45 lines on $M_D$ (Choi matrix of $\mathcal L$; $T_t$ CP ⇒ $\mathcal L$ conditionally completely positive, i.e. $P_\Omega^\perp J(\mathcal L)P_\Omega^\perp\ge0$; read off $\Phi$ and $K$; trace preservation fixes $K$). Or R. |
| 8.2 | 219–224 | derivation, *from Christensen–Evans*, that every generator on $\mathcal M=\bigoplus M_{k_c}$ has the form of §8.3 | Rests on 8.1's unproved theorem. | A: +15 lines in Appendix D (extend $T_t$ to $M_D$ via the conditional expectation $E(X)=\sum P_cXP_c$; $T_t\circ E$ is a CP semigroup on $M_D$ with $\hat T_0=E$, not a QDS; compress the GKSL form of its generator). Or R together with 8.1. |
| 8.3 | 226–230 (box) | Born–Markov approximation leads to semigroups | PHYS. | PHYS. |
| 8.4 | 408–414 (box) | counting-process interpretation (Davies, Srinivas–Davies, MCWF) | PHYS (text says only the math is used). | PHYS. |
| 8.5 | 560–563 | uniqueness criteria (Frigerio) "are not proved here" | Unproved, unused. | R (keep "the code finds one stationary state per sector" as N). |
| 8.6 | 147, 368–369, 475–477 | $1.1\cdot10^{-9}$, $\approx-2\cdot10^{-8}$, $0.230$, $0.091$, $0.135$ | Numerical outputs. | N: already outputs of code; label. 0 lines. |

## Chapter 9 — `09_spin_chains.md`

| ID | line | statement | what is missing | action |
|---|---|---|---|---|
| 9.1 | 88–91 | $\lVert\Lambda\rVert^2$ is "the index of the Markov trace (Jones; Goodman–de la Harpe–Jones)" | Identification unproved. | R: definition only (as 7.1). |
| 9.2 | 109–110 | Jones: $4\cos^2\frac\pi n$ are exactly the subfactor index values below 4 | Unproved (deep). | R. |
| 9.3 | 193–194 | "exchanges of neighbours … generate all permutations" | Standard but unproved (bubble sort). | P: 1–2 lines. |
| 9.4 | 219–221 (box) | QMC sign problem | PHYS. | PHYS. |
| 9.5 | 247–250 | general Lieb–Mattis theorem "(stated without proof)" | Unproved; not used. | R (A ~40 lines possible, as 4.11). |
| 9.6 | 250–252, 269–270, 274–275 | $L(E_0(1)-E_0(0))\to\pi^2/2$ with log corrections (des Cloizeaux–Pearson, Affleck et al.; "stated without proof") — **assert** `gap[12] < math.pi**2/2` | Unproved asymptotic (Bethe ansatz + CFT) used as an upper bound in an assert. | N/R: keep the computed values and the proven positivity; drop $\pi^2/2$ from the assert and the text (or "literature value, not proved"; not allowed under strict rule). 2 code + 2 text lines. Proof is out of reach. |
| 9.7 | 277–281 | **Wigner–Eckart theorem "(stated without proof)"**; default method of `spin_correlation` | Unproved; library correctness rests on it (checked numerically against `'direct'`). | A: Appendix G ~45–55 lines (spherical components of a vector operator, CG coefficients for $1\otimes S$, Schur via Theorem 4.3 on $\mathrm{Hom}$ spaces). Or N/R: present only the `'direct'` method (proved by Prop. 4.15(b), since $\mathbf s_i\cdot\mathbf s_j\in\{S_a\}'$) and call Wigner–Eckart an implementation detail verified by the code. |
| 9.8 | 309–310 | correlation decay $(-1)^r\sqrt{\log r}/r$ "(stated without proof)" | Unproved. | R. |
| 9.9 | 380, 424–430 | completeness of the 2-magnon Bethe ansatz "checked here, not proved"; "With $\Delta=0.6$ and $4\le L\le13$ the search above finds every state for odd $L$, and every state but $\Phi$ for even $L$" | The code checks only $L=9$ (and $\Phi$ for $L=8$); the claim for $4\le L\le13$ is not computed in the chapter. | N: add a loop over $L$ (≈5 code lines) or restrict the sentence to $L=8,9$. |
| 9.10 | 430–431 | $N$-magnon ansatz: $N!$ terms, $N$ coupled Bethe equations | Unproved pointer. | R or keep as reference without claim. |
| 9.11 | 433–437 (box), 509, 513–515 | factorised scattering; ground-state energy per site $\frac14-\log2$ (Hulthén, "stated without proof"); finite-$L$ values "approach the value … from below" | Unproved constant and asymptotic; printed by the code. | R/N: print only the computed energies; drop the infinite-chain value or label it explicitly as a literature value not used. |
| 9.12 | 466 | Courant–Fischer "(standard linear algebra)" | Allowed by README wording, but it is also the basis of Appendix A.3; one proof serves both. | A: Appendix A.3 (6 lines). Optional. |
| 9.13 | 474–476 | Kaniel–Paige–Saad convergence bounds "(stated without proof)" | Unproved, unused. | R. |
| 9.14 | 495–500 | reference energies "tabulated to nine digits in `scripts/chains.py` (reference values of the companion paper)" — $L=14,16$ values are **assert targets** (l.510) | Computed outside the chapter. | N: label as numerical values produced by repository scripts (dense/Lanczos) — acceptable as numerical; or compute $L=14$ in-chapter. 1 line. |
| 9.15 | 528–530 | `krylov.evolve` a-posteriori error estimate (Hochbruck–Lubich; "stated without proof") | Library claim. | R (or N: the code compares with the dense propagator for $L=12$). |
| 9.16 | 559–564 | "generically nothing else" is conserved; MBL: "for strong disorder the eigenstates localise … and the levels become independent, $\langle r\rangle=2\log2-1$" — **asserts** l.578 against 0.5307 and 0.3863 | Physics claims; literature constant as target. | PHYS for MBL; R "generically nothing else"; N for the asserts (compare with the course-computed GOE constant and with the proven Poisson value, labelled numerical). |
| 9.17 | 582–588 | "at $W=0$ the chain is integrable"; crossover near $W\approx3.7$ (Luitz et al.); debated (Šuntajs et al.) | Unproved / literature. | PHYS or R. |
| 9.18 | 636–638 | entanglement grows like $\frac16\log L$, "the CFT law for an open chain (Peschel; Calabrese–Cardy; stated without proof)" | Unproved asymptotic. | R, or N (fit with the chapter's Peschel code for several $L$, labelled numerical; ~6 code lines). |
| 9.19 | 689–690 | "a random start vector avoids this with probability one" | Needs Appendix A.1 (Gaussian vector has non-zero component in the ground eigenspace a.s.). | P: 1 line + pointer. |

## Chapter 10 — `10_quantum_computing.md`

| ID | line | statement | what is missing | action |
|---|---|---|---|---|
| 10.1 | 145–148 (box) | trapped-ion DFS experiment (Kielpinski et al.) | PHYS. | PHYS. |
| 10.2 | 307–316 | **Theorem 10.7 (Bény–Kempf–Kribs) "stated without proof"**; consequence: operator (subsystem) codes contain noiseless subsystems | Unproved. Only the two extreme cases (Thm. 10.5, Prop. 10.6) are proved. | A: Appendix E ~55–70 lines (finite dimensions: write $\mathfrak A=\bigoplus_cM_{k_c}\otimes1_{m_c}$; necessity as in KL with the algebra; sufficiency = Prop. 10.6-type sector measurement + KL recovery on each $\mathbb C^{k_c}$ factor with the subsystem tensor structure). Or R (keep 10.5, 10.6). |
| 10.3 | 347 (box) | "hybrid codes correspond to $\bigoplus_cM_{k_c}$" | Remark depending on 10.7. | R or keep after 10.2 is proved. |
| 10.4 | 487–490 (box) | $T_1$ a leading error of superconducting qubits | PHYS. | PHYS. |

Otherwise complete (Thm. 10.2, Cor. 10.3, Prop. 10.4, KL theorem, Prop. 10.6, 10.8–10.10).

## Chapter 11 — `11_optimisation.md`

| ID | line | statement | what is missing | action |
|---|---|---|---|---|
| 11.1 | 82–86 | $S\mapsto U_0e^S$ "is a local diffeomorphism near $S=0$"; "its differential is singular where two eigenvalues of $S$ differ by a nonzero multiple of $2\pi i$" | Inverse function theorem (not on the list) and the $d\exp$ formula, both unproved. Motivation only. | R (keep "re-centring keeps the chart near $S=0$"), or P ~10 lines ($d\exp_S$ in the eigenbasis of $S$ multiplies entries by $(e^{s_i}-e^{s_j})/(s_i-s_j)$). |
| 11.2 | 86–88 | Cayley chart and polar retraction "have the same differential at 0" | Unproved. | P: 3 lines, or R. |
| 11.3 | 335–338 | concentration of measure on $U(k)$, $\mathbb P(f\ge\mathbb Ef+t)\le e^{-ckt^2/L^2}$ (Meckes; "stated without proof") | Unproved. | R (the proven Markov bound of Cor. 11.12 suffices). A would need log-Sobolev/Bakry–Émery on $SU(k)$: ~120–200 lines. |
| 11.4 | 366 | "The ratio falls like $1/k$" | Only the upper bound $b^2/(a^2k)$ is proved. | N: "is below the bound $b^2/(a^2k)$ and, numerically, falls like $1/k$". 0 lines. |
| 11.5 | 411–414 | structured samples "come close to the extremes in every dimension, unlike Haar samples (`scripts/bounds_checks.py`)" | Unproved, outside the course. | R or N. |
| 11.6 | 395–398 (box) | typicality in statistical mechanics | PHYS. | PHYS. |
| 11.7 | 487–493 | closed form $(1-\delta)(2d/(1+d))^2$ observed numerically, "could not prove; the code does not assert it" | Already labelled numerical. | N (compliant). 0 lines. |

## Chapter 12 — `12_equivariant_learning.md`

| ID | line | statement | what is missing | action |
|---|---|---|---|---|
| 12.1 | 336–338 | spectral layers get "gradients by the Daleckii–Krein formula (`apply_function`)" | Library/autograd claim; differentiability of $f(A)$ and the DK formula unproved (the code checks one derivative). | R ("gradients by automatic differentiation; checked below"), or A ~25 lines (DK formula for $C^1$ $f$: polynomials first, then uniform approximation on a neighbourhood of the spectrum). |
| 12.2 | 454–456 | Molien's formula "(stated without proof)" | Unproved; Prop. 12.7(c) proves the $t^2$ coefficient only. | P ~12 lines: same proof as 12.7(c) with $\mathrm{Sym}^d$: $\mathrm{Tr}\,\sigma_d(g)=h_d(\lambda)$, $\sum_dh_d(\lambda)t^d=\prod_i(1-t\lambda_i)^{-1}=\det(1-tA)^{-1}$. |
| 12.3 | 458–462 (box) | Landau theory: "generically the irreducible representation of that sector determines the symmetry of the ordered phase" | PHYS ("generically" unproved). | PHYS. |
| 12.4 | 472–475 | Kondor–Trivedi: equivariant layers on homogeneous spaces of compact groups are generalised convolutions | A theorem quoted in "What to compare with". | R (keep as citation without the statement) — 0 lines. |
| 12.5 | 501–502, 553–557 (Ex.4) | Gaussian $x_i$ ⇒ $V^*x_i$ Gaussian; $\mathrm{rank}\,S_c=\min(Nk_c,m_c)$ a.s. | Uses Gaussian invariance + polynomial null set (proved in §4.3 only modulo Tonelli/Lebesgue). | A: pointer to Appendix A.1–A.2. 0 lines. |

## Chapter 13 — `13_jordan_algebras.md`

| ID | line | statement | what is missing | action |
|---|---|---|---|---|
| 13.1 | 55–57 | "the fundamental formula holds in every Jordan algebra (McCrimmon), which we do not reproduce" | Unproved; not used. | R (A = Macdonald's theorem, 150+ lines). |
| 13.2 | 95–97 (box) | $\frac12(xy+yx)$ is "the only combination $\alpha xy+\beta yx$ symmetric in $x,y$ and reducing to $x^2$" | Trivial but unproved. | P: 1 line. |
| 13.3 | 127–128 | "Up to signs, $s_{2j-1},s_{2j}$ are the Majorana operators $c_j+c_j^*$, $i(c_j-c_j^*)$" | Unchecked (and basis conventions differ: Ch.9 strings use $(-1)^{n_l}$, spin-up = bit 1). | P: 3 lines, or R. |
| 13.4 | 196–198 (box) | states of $V_m$ are $x\mapsto\mathrm{tr}(\rho\circ x)$, $\rho=\frac12(1,r)$, $\lvert r\rvert\le1$ | Riesz + self-duality, not written out. | P: 3 lines. |
| 13.5 | 202–204, 286–289 | interior-point methods for SOCP "use exactly Theorem 13.4"; "interior-point methods can be written once for all Euclidean Jordan algebras" | Algorithmic/literature claims. | R or rephrase as citation. |
| 13.6 | 257–259 | formally real ⇔ Euclidean in finite dimensions [FK94] | Unproved (only "Euclidean ⇒ formally real" is 3 lines). | R (define only "Euclidean"; show $M_{sa}$, $V_m$ Euclidean as now), or A ~40–60 lines (spectral theory of formally real Jordan algebras, positivity of the trace form). |
| 13.7 | 262–269 | **Theorem 13.7 Jordan–von Neumann–Wigner classification "stated without proof"**: decomposition into simples, the list, uniqueness, rank-2 isomorphisms $\mathrm{Herm}_2(\mathbb K)\cong V_{\dim\mathbb K+1}$ for $\mathbb K=\mathbb R,\mathbb C,\mathbb H,\mathbb O$ **[README]** | Unproved. Honest size of a complete proof: spectral theorem for EJAs (~100), Peirce decomposition w.r.t. a Jordan frame (~80), simple ⇒ Peirce-connected (~40), rank 2 ⇒ spin factor (~40), coordinatisation for rank ≥3 (Jacobson) with composition-algebra structure (~150–200), Hurwitz theorem on composition algebras incl. Frobenius for associative ones (~100), associativity for $r\ge4$ / alternativity for $r=3$ (~50), $\mathrm{Herm}_3(\mathbb O)$ satisfies the Jordan identity (~60–100), uniqueness by dimension/rank (~20): **~650–900 lines**. | Recommend R: move to a clearly labelled "Outlook (not proved, not used)" or delete, and update README R.3. Cheap partial P/A: rank-2 isomorphisms for $\mathbb R,\mathbb C,\mathbb H$ (~15 lines; needs definition of $\mathbb H$), $\mathbb O$ (+10 with octonion multiplication table). |
| 13.8 | 271–272 | $\mathrm{Herm}_r(\mathbb H)$ embeds into $\mathrm{Herm}_{2r}(\mathbb C)$ (special) | Unproved; $\mathbb H$ never defined. | P ~8 lines (quaternions as $2\times2$ complex matrices; Hermitian ↦ Hermitian; $\circ$ preserved). |
| 13.9 | 273–280 | Albert algebra exceptional (Albert 1934); Glennie identities; $\mathrm{Herm}_r(\mathbb O)$ Jordan only for $r\le3$; "cannot be represented faithfully … inside any associative algebra" | Unproved (deep; Glennie's identity verification alone is long). | R (A ≥ 200 lines). |
| 13.10 | 282–289 | **Theorem 13.8 Koecher–Vinberg "stated without proof"**; "So the cones of SDP, SOCP and LP are all of this form" | Unproved. Only the easy direction for specific cones is needed. | R the theorem (A ~250 lines: EJA theory + Vinberg's characteristic function); P ~15 lines: $\mathrm{Herm}_n(\mathbb R)^+$, $\mathrm{Herm}_n(\mathbb C)^+$, $\mathbb R_+^n$ are self-dual and homogeneous (as Ex.1: $P(x^{1/2})1=x$, trace pairing). |
| 13.11 | 312–317 | **Theorem 13.10 Kadison "stated without proof"** (unital C\*-algebras), and the isometry theorem | Unproved. | A/P ~30–35 lines in finite dimensions (all that the course needs): a unital order isomorphism maps $\{0\le a\le1\}$ onto itself, hence extreme points (= projections, 8 lines) to projections; $p+q\le1$ ⇒ $pq=0$ (3 lines), so orthogonal projections go to orthogonal projections; spectral decomposition ⇒ $\Phi(x^2)=\Phi(x)^2$ (4 lines); complexify (3 lines). R the C\*-algebra statement and the isometry theorem. |
| 13.12 | 318–322 | every Jordan $*$-automorphism of $M_n(\mathbb C)$ is $UAU^*$ or $UA^TU^*$ (Jacobson–Rickart, Herstein) | Proved only for $n=2$ (Prop. 13.11). | R (state only $n=2$), or A ~60–80 lines (Jordan automorphisms preserve rank-one projections and $\mathrm{Tr}(pq)$ ⇒ Wigner's theorem, Bargmann-style proof ~45 lines). |
| 13.13 | 376–378 | Hanche-Olsen; Barnum–Graydon–Wilce "for results on which composites exist" | Pointer to results. | R or keep as reading list without claims. |
| 13.14 | 396–399 | Alfsen–Shultz: JBW-algebra is the self-adjoint part of a vN algebra iff it has a dynamical correspondence; state-space characterisation "(stated without proof)" | Unproved. | R. |
| 13.15 | 447–451 (box) | transpose = time reversal for a spinless particle, antiunitary; "quantum mechanics uses complex $*$-algebras because …" | PHYS/interpretive. | PHYS. |
| 13.16 | 483–494 | in every EJA: Jordan-frame spectral theorem, $\mathrm{tr}z=\sum\mu_i$, associativity of the trace form, $P(y)Q\subseteq Q$ "[FK94]; we do not prove these facts"; conclusion "So the inequality holds in every Euclidean Jordan algebra — including $\mathrm{Herm}_r(\mathbb H)$ and the Albert algebra" | Conclusion rests on unproved facts. | P (0–2 lines): state the proven conditional "in every Jordan algebra with properties (i)–(iii)"; R the unconditional claim. (A for the facts: ~120 lines of EJA spectral theory.) |
| 13.17 | 566–569 (Ex.2) | Horodecki: for $2\times2$ systems PPT ⇒ separable "we do not prove it" | Unproved. | R (A ≥ 150 lines: Størmer–Woronowicz decomposability of positive maps $M_2\to M_2$). |

---

## Totals

Counted items: **119 entries** in 15 files (README 4, Ch.0 2, Ch.1 3, Ch.2 5, Ch.3 2, Ch.4 12, Ch.5 23, Ch.6 2, Ch.7 8, Ch.8 6, Ch.9 19, Ch.10 4, Ch.11 7, Ch.12 5, Ch.13 17).

| category | count | main items |
|---|---|---|
| Major theorems stated without proof | 12 | Schur–Weyl (4.14b), general Weingarten (5.12), GOE/GUE joint density + Gaudin–Mehta, Page (7.13), GKSL converse (8.7) + derivation for $\mathcal M$, Wigner–Eckart, Bény–Kempf–Kribs (10.7), JvNW (13.7), Koecher–Vinberg (13.8), Kadison (13.10), Jacobson–Rickart/Herstein, Lieb–Mattis (general) |
| Measure-theory / analysis beyond the README list | ~12 uses | Fubini–Tonelli, DCT, π–λ, change of variables, polar coordinates, Lebesgue invariance, inverse function theorem, eigenvalue continuity, Gaussian push-forward densities |
| Literature constants / asymptotics as assert targets | 7 asserts | 0.5307 (×4: l.5-475, 5-550, 5-682, 9-578), Page value (7-570), $\pi^2/2$ (9-270), RING table (9-510) |
| Other stated-without-proof remarks (unused) | ~25 | Haar existence, general Stinespring, Feller, Dyson COE, Jones index theorem, KPS, Hochbruck–Lubich, concentration of measure, Frigerio, Molien, Daleckii–Krein, Hulthén, CFT laws, des Cloizeaux–Pearson, Albert/Glennie, Alfsen–Shultz, Horodecki, fundamental formula, … |
| Short gaps fixable in place (P) | ~25 | 2.2 (numerical rank used as premise), 2.5, 4.5, 4.6, 5.5, 5.9, 5.23, 7.7, 9.3, 9.19, 11.2, 12.2, 13.2–13.4, 13.8, … |
| Physical/empirical claims (PHYS) | ~20 | MBL, BCS, CFT, Page curve/black holes, Born–Markov, BGS/Berry–Tabor, QMC sign problem, $T_1$/$T_2$, trapped ions, typicality, Landau theory, time reversal, … |

Note on the "decompose" algorithm (user's example): in exact arithmetic the claim "a random element of the
commutant finds the decomposition with probability 1" **is** proved (Lemma 4.6, Prop. 4.7, Lemma 4.8, Prop. 4.9),
modulo Appendix A (4.1–4.3: polynomial null sets, Gaussian push-forwards, eigenvalue continuity). Unproved are only the
floating-point refinement (4.4) and the generic size $N$ of the reduced generator problem (4.5).

## Recommended plan

**Step 0 — README (≈6 lines).** Extend the allowed list to standard measure and integration theory and
elementary calculus (R.1); add the *Physical meaning* disclaimer (R.4); fix the Ch.7 and Ch.13 rows (R.2, R.3).

**Appendices to create (core, recommended):**

| appendix | content | lines |
|---|---|---|
| **A. Measure-theoretic and analytic facts used** | A.0 list of quoted theorems (Lebesgue measure, Fubini–Tonelli, DCT, change of variables incl. polar coordinates, π–λ); A.1 null sets of polynomials (4.1, 5.3, 9.19, 12.5); A.2 Gaussian vectors: orthogonal invariance, push-forward under linear surjections has a density (4.3, 2.4); A.3 Courant–Fischer ⇒ Weyl's inequality ⇒ continuity of ordered eigenvalues and of spectral projections on the simple-spectrum set (4.2, 5.15, 9.12) | 60–80 |
| **B. Schur–Weyl duality and the Weingarten calculus** | polarisation; commutant of $\{V^{\otimes n}\}$ = span of permutations; Theorem 4.14(b) as corollary; Theorem 5.12 complete (Gram matrix, independence for $n\le k$, pseudo-inverse, class function, unbalanced moments); Catalan identity optional | 100–130 |
| **D. The GKSL theorem in finite dimensions** | converse on $M_D$ via conditional complete positivity of the Choi matrix; extension to $\bigoplus M_{k_c}$ | 60–70 |
| **F. Jordan algebras: what is proved** | Kadison's theorem in finite dimensions (13.11); $\mathrm{Herm}_r(\mathbb H)\subset\mathrm{Herm}_{2r}(\mathbb C)$ (13.8); rank-2 isomorphisms for $\mathbb R,\mathbb C,\mathbb H$; self-dual homogeneous cones for SDP/LP (13.10 easy part) | 70–90 |

**Optional appendices (only if the author wants to keep the statements):**

| appendix | content | lines |
|---|---|---|
| C. Eigenvalue densities and Page's formula | Weyl integration formula for GOE/GUE/Laguerre; level density bounded continuous; Page via Laguerre polynomials | 220–330 |
| E. Operator-algebra error correction | Bény–Kempf–Kribs in finite dimensions | 55–70 |
| G. Wigner–Eckart for SU(2) vector operators | spherical components, CG for $1\otimes S$, theorem; justifies `su2.spin_correlation` | 45–55 |
| H. Further Jordan theory | EJA spectral theorem + trace associativity (13.6, 13.16) ~120; Jacobson–Rickart/Herstein via Wigner's theorem ~70; full JvNW classification incl. Peirce decomposition, coordinatisation, Hurwitz ~650–900; Koecher–Vinberg ~250; Albert/Glennie ≥200 | 1300–1550 |
| — | general Lieb–Mattis ~40; concentration of measure ~150; Daleckii–Krein ~25; existence of Haar measure on compact groups ~100 | ~315 |

**In-place work.** P items ≈ 25 × 1–15 lines ≈ **110–140 lines**. N items (course-computed GOE/GUE $\langle r\rangle$
constant used by all four asserts, COE comparison, MBL asserts, $\pi^2/2$ assert removal, Page assert → proven bounds,
Poisson levels instead of uniform levels, Bethe loop over $L$, labels on printed numbers) ≈ **40–60 lines** of code and
text. R items: mostly deletions (net −80 to −120 lines), notably Theorem 13.7/13.8/Alfsen–Shultz/Albert paragraph,
general Stinespring, Haar existence, Gaudin–Mehta, the "many sectors" derivation, concentration of measure, CFT and
Hulthén values, KPS, Hochbruck–Lubich, Frigerio, Horodecki, fundamental formula, Jones' index theorem.

**Estimated total added lines.**
- Recommended self-contained version (README fix + Appendices A, B, D, F + P + N, with R for the deep theorems):
  **≈ 450–560 lines added**, ≈ 100 lines removed. (Without Appendix D, i.e. R for the GKSL converse: ≈ 390–490.)
- Keeping every current statement with full proofs (adding C, E, G, H and the misc. items): **≈ 2 800–3 300 lines**
  added — dominated by the Jordan–von Neumann–Wigner classification (650–900), Koecher–Vinberg (~250) and
  Page/eigenvalue densities (220–330). Not recommended for a course whose stated scope is
  $\bigoplus_cM_{k_c}(\mathbb C)$.

Priority order (premises and asserts first): 2.2, 4.2–4.3 + A, 5.13–5.14/5.22/9.16 (constants as assert targets),
7.6 (Page assert), 9.6 ($\pi^2/2$ assert), 4.7/5.8 + B, 5.11 (level density premise), 8.1–8.2 + D, 13.11 + F,
README; then the stated-only remarks (R) and the PHYS disclaimer.

Note: the six translations (`es`, `fr`, `de`, `zh`, `ja`, `ru`) will need the same edits.
