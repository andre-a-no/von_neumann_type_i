# `torch_vn_algebra.certify`: plan and status

Goal: replace "numerical search gives ≈ …" in the papers by "Lemma (computer-assisted)" with a reference to a
certificate. For an envelope V(δ) = sup F the module returns a rigorous interval [L, U] ∋ V(δ) and a certificate file
that a small independent program checks.

## Why

Before this module the two halves were separate. Lower bounds came from `optimize.extremize` and the mpmath
recomputation in `scripts/trace_characterizations.py`. Upper bounds came from `scripts/certify_envelopes.py`
(hand-derived formulas, branch and bound with hand constants, float64 with a tolerance, so rounding was not
controlled). As a result 9 of the 25 envelopes of the trace preprint were not marked as exact: items 6, 7, 10, 12, 13,
14, 15, 16 and 25.

## Design

| part | method | status |
|---|---|---|
| lower bound L | witness of the original problem in exact rationals: orthogonal/unitary matrices by Cayley, U = (1 − K)(1 + K)⁻¹ (real) or (1 − iK)(1 + iK)⁻¹ (complex, K Hermitian), positive matrices U diag(qᵢ) U*; violation checked without rounding | real case done (`witness.py`); complex case to do |
| upper bound U in fixed dimension | branch and bound, Arb ball bounds on each box: the smallest of the natural extension, the mean-value form and the second-order Taylor form (Hessian enclosure by forward AD), monotonicity test, smear splitting. The plan had two bounds; the Taylor form was added because item 14 has a ridge just below L along m = m′, where the mean-value form needs 5 times more boxes | done (`engine.py`, `ad.py`) |
| matrix functions | Lipschitz bounds for \|X\|, X^{1/2}, e^X, A♮B (Kittaneh-type inequalities): value at the centre plus radius times the Lipschitz constant, no hand-derived curvature | to do |
| reductions | unitary invariance, phases by the triangle inequality, compression, two-level spectra, the two-projection theorem (Halmos) | listed per problem with references; used for item 14 |
| all finite von Neumann algebras | trace-polynomial SDP hierarchy (Klep, Magron, Volčič), rational rounding of the Gram matrices (Peyrl, Parrilo); if it meets a matrix witness, the constant is sharp for all finite algebras, II₁ included | to do (items 7, 10, 14, 25) |
| certificate | gzipped JSON: problem, rational δ, reductions, witness, preorder search tree (or rational Gram matrices), versions | done, format `torch_vn_algebra-certificate/1` |
| verifier | `torch_vn_algebra/certify/verify.py`, python-flint and the standard library only, does not import the package | done for box certificates; Gram-matrix check to do |
| claims | `certificates/claims.json`, `tests/test_certified_claims.py` (the existing `tests/test_article_claims.py` covers the library paper) | done for item 14 |

## Papers

The trace preprint is split into four papers; the certify work feeds two of them.

| paper | items | what certify provides |
|---|---|---|
| Б (B) | 15 | a short computer-assisted lemma on M₂ |
| Г (D) | 6, 7, 10, 12, 13, 16, 25 | the whole paper rests on the module: box certificates in fixed dimension, the SDP layer for 7, 10, 25, the analytic lemma near V = 0 for 6 and 25 |

The earlier plan of a single submission to LAA is dropped. The question of arbitrary ranks for item 14 is closed:
the reduction below holds for every n and all ranks of P and Q.

## Item 14 (done)

Hölder's inequality for projections, α = 3. Reduced to three variables on [0, 1]² × [(1 − δ)/(1 + δ), 1] for every n
and all ranks (two-projection theorem, compression, Hölder, phases; module docstring of `problems.py`). 13
certificates on the δ grid of the preprint, relative width 10⁻⁶, from 5 s (δ = 0.95) to 4 min (δ = 0.01, 3.6·10⁵
boxes). Example: 2.79878567081·10⁻⁴ ≤ V(0.1) ≤ 2.79878846960·10⁻⁴; the earlier prototype gave
[2.79877975·10⁻⁴, 2.79878570·10⁻⁴], consistent with it (intersection [2.79878567·10⁻⁴, 2.79878570·10⁻⁴]). The search values of the preprint are all below
L, by up to 17 % at δ = 0.01 (n = 4), so the table should quote the certified values.

The cost grows like δ⁻¹: for small δ the reduced function has a ridge of height ≈ 0.97 L along the diagonal m = m′,
where F ≡ 0. A reduction that removes this direction analytically would make small δ cheap.

## Next

1. Port the remaining items to `problems.py`: 15 (short computer lemma on M₂ for paper B), then 6, 7, 10, 12, 13, 16,
   25 (paper D). Each needs its reductions written out with references.
2. Lipschitz box bounds for matrix functions (items 6, 12, 13, 15, 16).
3. Complex Cayley witnesses (Gaussian rationals).
4. SDP layer with Peyrl–Parrilo rounding and the Gram-matrix check in the verifier (items 7, 10, 14, 25).
5. Items with V ≈ 0 at small δ (6 and 25): branch and bound cannot prove F ≤ 0 near the equality set; a separate
   analytic lemma near that set is needed, the boxes cover the rest.
