# Contributing to torch_vn_algebra

Thank you for your interest. Bug reports, questions, examples and pull requests are welcome.

## Reporting a problem

Open an issue at <https://github.com/andre-a-no/von_neumann_type_i/issues> with

- a minimal script that reproduces the problem (algebra sizes, precision, device),
- what you expected (ideally a dense NumPy/SciPy reference value) and what you got,
- the versions printed by `python -c "import torch, torch_vn_algebra as t; print(t.__version__, torch.__version__)"`.

Wrong numerical results have the highest priority; please say so in the title.

## Development set-up

```bash
git clone https://github.com/andre-a-no/von_neumann_type_i
cd von_neumann_type_i
pip install -e ".[test,scripts]"
pytest -q                      # about a minute on a CPU; GPU tests run when CUDA is available
```

## Pull requests

1. One topic per pull request, with a test that fails before the change and passes after it.
2. Numerical code is checked against an independent reference: a dense matrix built with NumPy/SciPy, an exact
   formula, or a brute-force computation (see `tests/test_regressions.py` for the pattern).
3. Keep the design rules of the library:
   - operators are tensors `(batch, C, k_max, k_max)` with zero padding; decompositions run as one batched call
     per distinct block size (`Operator._groups`, `TypeIAlgebra._size_groups`), never in a Python loop over sectors;
   - no host synchronisations (`.item()`, `bool(tensor)`, `.tolist()`) in computational paths
     (`tests/test_sync_free.py` counts them);
   - property flags set by the library must be mathematically guaranteed;
   - single and double precision, real and complex algebras, CPU and GPU.
4. Run `pytest -q` before pushing; CI runs the test suite on Python 3.9–3.12.
5. Changes visible to users go into `CHANGELOG.md` under "Unreleased".

## Documentation

All documentation lives in `man/`: the lecture notes (`man/lecture_notes`, the theory with complete proofs, no
code), the tutorial (`man/tutorial`, the library in the context of the theory) and the paper (`man/article`). The
lecture notes and the tutorial have one file per chapter with the same name and link to each other through explicit
anchors (`<a id="r-4-20"></a>` before Theorem 4.20, `<a id="s-4-7"></a>` before Section 4.7). Both exist in seven
languages; the English text is authoritative and `tests/test_man.py` runs every code block of the tutorial and checks
that code, displayed formulas, headings, anchors and links agree across translations. A change to an English chapter
needs the same change in every translation (or an issue asking for it). Formulas use GitHub's verbatim math syntax:
inline ``$`\mathrm{Tr}\,A`$``, displays in a ```` ```math ```` block; plain `$...$` and `$$...$$` lose `\,`, `\{` and `^*` to
Markdown on GitHub. Write `\mathrm{Tr}`, not `\operatorname{Tr}`.

## Code of conduct

Participation in this project is governed by the [Code of Conduct](CODE_OF_CONDUCT.md).
