# `torch_vn_algebra.certify`: plan and status

This file is the plan of the module; there is no separate plan document. Status as of 2026-10-31.

Goal: replace "numerical search gives ≈ …" in the papers by "Lemma (computer-assisted)" with a reference to a
certificate. For an envelope V(δ) = sup F the module returns a rigorous two-sided interval [L, U] ∋ V(δ) and a
certificate file that a small independent program checks.

## Why

Before this module the two halves were separate.

- Lower bounds: `optimize.extremize` finds points, `scripts/trace_characterizations.py` recomputes them in mpmath.
  These are lower bounds only.
- Upper bounds: `scripts/certify_envelopes.py` has hand-derived formulas and a branch and bound with hand constants,
  in float64 with a tolerance, so rounding is not controlled.

As a result 9 of the 25 envelopes of the trace preprint are not marked as exact: items 6, 7, 10, 12, 13, 14, 15, 16
and 25.

## Design and status

| part | method | status |
|---|---|---|
| lower bound L: exact witness | the numerical point is replaced by an exact rational configuration of the original problem: orthogonal/unitary matrices by Cayley, U = (1 − K)(1 + K)⁻¹ (real) or (1 − iK)(1 + iK)⁻¹ (complex, K rational Hermitian), positive matrices U diag(qᵢ) U*; the violation is checked without rounding | real case done (`witness.py`, item 14); complex case to do |
| lower bound L: Arb recomputation | for problems where no exact rational witness exists (matrix functions with irrational values), the found point is recomputed in Arb ball arithmetic and L is the lower end of the ball | to do |
| upper bound U in fixed dimension | branch and bound on boxes with rational endpoints, Arb bounds on every box: the smallest of the natural interval extension, the mean-value form and the second-order Taylor form (Hessian enclosure by forward AD); monotonicity test; smear splitting | done (`engine.py`, `ad.py`). The original design had two bounds; the Taylor form was added because item 14 has a ridge just below L along m = m′, where the mean-value form alone needs 5 times more boxes |
| matrix functions | for \|X\|, X^{1/2}, e^X, A♮B: value at the box centre plus the radius times a Lipschitz constant from Kittaneh-type inequalities, so no curvature is derived by hand | to do |
| reductions | applied before the search: unitary invariance, removal of phases by the triangle inequality, compression to a subspace, reduction to two-level spectra, the two-projection theorem (Halmos) | written by hand per problem, with proofs in the docstring and references in the certificate (item 14). A reusable library of reductions applied by the engine: to do |
| all finite von Neumann algebras | polynomial problems (items 7, 10, 14, 25): trace-polynomial SDP hierarchy of Klep, Magron and Volčič; the numerical solution becomes a proof by rational rounding of the Gram matrices (Peyrl and Parrilo). If the bound meets a matrix witness, the constant is sharp for all finite algebras, II₁ included | to do |
| certificate | gzipped JSON: problem, rational δ, reductions, witness, preorder search tree or rational Gram matrices, software versions | done for box certificates, format `torch_vn_algebra-certificate/1`; Gram-matrix part to do |
| verifier | `torch_vn_algebra/certify/verify.py`: python-flint and the standard library only, a few hundred lines (280 now), does not import the package; replays the tree with exact rational boxes, checks the witness exactly; the reductions are lemmas of the text and are not checked | done for box certificates; Gram-matrix check to do |
| claims of the papers | every numerical statement of the papers checked against its certificate: `certificates/claims.json` and `tests/test_certified_claims.py` (the name `tests/test_article_claims.py` is taken by the library paper) | done for item 14 (3 table values and all search values of the preprint); the other items to do |

## Papers

The trace preprint is split into four papers. A single submission to LAA is no longer planned.

| paper | certify items | what certify provides |
|---|---|---|
| Б (B) | 15 | a short computer-assisted lemma on M₂ |
| Г (D) | 6, 7, 10, 12, 13, 16, 25 | the whole paper rests on the module: box certificates in fixed dimension, the SDP layer for 7, 10, 25, the analytic lemma near V = 0 for 6 and 25 |
| А (A), В (C) | none assigned | — |

## Item 14 (done)

Hölder's inequality for projections, |φ(PQ)| ≤ φ(P)^{1/3} φ(Q)^{2/3}. The question of arbitrary ranks is closed: the
two-projection theorem, compression and Hölder's inequality reduce every n and all ranks of P and Q to one 2×2 block,
V_n = V_2 for n ≥ 2. With the phases removed and half-angle coordinates this is a problem in three variables on
[0, 1]² × [(1 − δ)/(1 + δ), 1] (proof in the module docstring of `problems.py`).

- 13 certificates on the δ grid of the preprint (`certificates/`, 172 KB), relative width 10⁻⁶; all pass the verifier
  (`VN_VERIFY_CERTIFICATES=1`, 3.5 min).
- Example: 2.79878567081·10⁻⁴ ≤ V(0.1) ≤ 2.79878846960·10⁻⁴. The earlier prototype gave
  [2.79877975·10⁻⁴, 2.79878570·10⁻⁴]; the two agree (intersection [2.79878567·10⁻⁴, 2.79878570·10⁻⁴]).
- The search values of the preprint are all ≤ U, so they are correct lower bounds, but all are below L, by up to 17 %
  at δ = 0.01 (n = 4). The three table values (2.8·10⁻⁴, 0.0086, 0.0711) agree with the certified intervals.
- Cost: from 5 s and 5.6·10³ boxes (δ = 0.95) to 4 min and 3.6·10⁵ boxes (δ = 0.01), roughly ∝ δ⁻¹. For small δ the
  reduced function has a ridge of height ≈ 0.97 L along the diagonal m = m′, where F ≡ 0.

## Next

1. Item 15 (paper Б): the short computer lemma on M₂. Needs the Arb recomputation for e^X and the Lipschitz bound
   for the matrix exponential.
2. Paper Г: items 12, 13, 16 (box certificates with Lipschitz bounds for matrix functions), 7, 10, 25 (SDP layer), 6.
   Each needs its reductions written out with references.
3. SDP layer: trace-polynomial hierarchy, Peyrl–Parrilo rounding, Gram-matrix check in the verifier. Item 14 as
   the first test case: the box certificates give the value it has to reach.
4. Items with V ≈ 0 at small δ (6 and 25): branch and bound cannot prove F ≤ 0 near the equality set; a separate
   analytic lemma is needed there, and the boxes cover the rest of the domain.
5. Complex Cayley witnesses (Gaussian rationals) and a reusable library of reductions applied by the engine.
6. The trace preprint (and paper Г): quote the certified values of item 14 in the table and the figure, and
   replace "numerical search gives ≈ …" by computer-assisted lemmas with certificate references; one claims entry
   per number.
7. `scripts/certify_envelopes.py`: re-run its float64 branch and bound (block model, k = 4, 6) through the engine,
   so the upper bounds of the library paper are certified too.
8. Optional speed-up for small δ: a reduction that removes the direction along the ridge of item 14 analytically.
