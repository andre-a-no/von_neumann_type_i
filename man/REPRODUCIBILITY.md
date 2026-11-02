# Numerical results and source revisions

The papers describe version 1.2.0 and the corrections listed under **Unreleased** in
[CHANGELOG.md](https://github.com/andre-a-no/von_neumann_type_i/blob/main/CHANGELOG.md). Their stored numerical tables and research figures are historical results.
Rebuilding the PDFs updates the text and captions; it does not repeat the numerical experiments.

The numerical corrections have regression tests, including analytic and dense references in
`tests/test_article_contracts.py`, `tests/test_review_regressions.py` and `tests/test_followup_regressions.py`.
Passing these tests does not establish that every published table or benchmark is reproduced by the corrected code.
Analytic theorems and their proofs are separate from these numerical validation claims.

## Stored runs

The following provenance comes from the `environment` and `mode` fields of the committed JSON files in
`man/article/generated/`. These results have **not been regenerated after the numerical corrections**.

| Study | JSON | Recorded revision | Mode |
| --- | --- | --- | --- |
| Connes–Størmer | `connes_stormer.json` | `e8d5873` | full |
| Boundary time crystal | `permutational.json` | `c33eb35+dirty` | check |
| Symmetry-resolved entanglement | `symmetry_resolved.json` | `b7568f6+dirty` | full |
| Trace characterisations | `trace_characterizations.json` | `38fd587+dirty` | full |
| Rényi thresholds | `renyi_threshold.json` | `3197dfc` | full |
| Channel learning | `channel_learning.json` | `497ba6e` | one full-size configuration |
| Residual networks | `residual_markov.json` | `39d970d+dirty` | one full-size configuration |

`+dirty` reports source changes that were not committed. The revision string does not preserve those changes;
checking out that commit alone is insufficient to reconstruct the run. Some older GPU records, including
`benchmark.json` and `numerics.json`, contain environment information but no source revision. This limits their
provenance too. The absence of a new run is an evidence gap, not a demonstrated error in the stored numbers.

## High-precision evaluation and certification

The trace, Connes–Størmer and Rényi scripts independently evaluate explicit points with ordinary mpmath arithmetic
at 40–50 significant digits. They re-normalise unitaries, partitions or Kraus operators at that precision and compare
with working-precision values. These are high-precision numerical checks, without an interval enclosure or a proved
bound on the final rounding error. For compatibility, existing helper names such as `certify` and JSON keys such as
`certified` retain their historical spelling; they mean **high-precision recomputation**, not a rigorous certificate.

The exact objective at any exactly admissible point is a lower bound for its supremum. Its floating-point
approximation is not automatically a rigorous numerical lower bound. Such a certificate would require a guaranteed
error enclosure and use its lower endpoint. Local optimisation supplies no global upper bound; absence of a found
violation does not prove an inequality. Analytic bounds and separately proved branch-and-bound guarantees are not
reclassified by this terminology change.

## Releasing updated numerical results

1. Commit the source and record a clean revision, dependency versions, device, precision, seeds and run parameters.
   Keep validation outputs separate from the historical data until the comparison is reviewed.
2. Run the regression suite. Repeat the affected studies: spectral calculus and trainable channels affect channel
   learning and Rényi searches; operator flags, caches and norms affect trace searches; Krylov normalisation affects
   chain calculations. Repeat GPU benchmarks on a GPU before quoting performance for the corrected code.
3. Use the documented script modes and preserve the original CPU/GPU configuration when comparing with a published
   run. In particular, the stored BTC results use **check** mode, while symmetry-resolved entanglement uses **full**.
   `bash scripts/run_all.sh full cuda` requests a different BTC configuration and does not reproduce that CPU run.
4. Save the new JSON, compare numerical differences and statistical or optimisation tolerances with the original
   results, and update the generated tables and figures together. Report any changes to scientific conclusions.
5. Rebuild all four paper PDFs and inspect the numerical status statements. Remove the historical-run qualification
   only for results actually regenerated and checked on the identified revision.

The entry points and reproduction commands remain in each paper and in `scripts/run_all.sh`. A quick or check run
can validate execution and selected reference cases; it must not be presented as repetition of a full-size study.
