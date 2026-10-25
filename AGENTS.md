# Instructions for AI coding agents

This file is for agents (Claude Code, Codex, Copilot, Cursor and others) that work *in this repository*. To *use* the
library, read [`torch_vn_algebra/guide.md`](torch_vn_algebra/guide.md) (`python -m torch_vn_algebra.guide`). The
rules for people are in [`CONTRIBUTING.md`](CONTRIBUTING.md) and apply to agents as well. This file adds what an agent
tends to get wrong.

## Layout

| path | contents |
|---|---|
| `torch_vn_algebra/` | the library: one module per topic. `algebra.py` holds `TypeIAlgebra` and its nested `Operator` |
| `tests/` | pytest suite: about 5000 tests, a few minutes on a CPU. `test_man.py` runs every code block of the tutorial and `test_guide.py` those of the guide |
| `man/lecture_notes/`, `man/tutorial/` | theory (no code) and tutorial (executed code), one file per chapter, in seven languages |
| `man/article/` | the library paper. Its numbers come from `generated/*.tex`, which the scripts write |
| `man/trace_envelopes/`, `man/connes_stormer/`, `man/btc_sre/` | the other preprints, each with its own generated numbers |
| `scripts/` | experiments that write the papers' tables and figures (`run_all.sh`), plus the arXiv packaging (`make_arxiv*.sh`, `flatten_arxiv.py`) |
| `docs_site/hooks.py`, `mkdocs.yml` | the documentation website. The API pages are generated from the docstrings (ast, nothing imported) |
| `llms.txt` | index for language models, published at the site root together with `llms-full.txt` |

## Commands

```bash
pip install -e ".[test,scripts]"
pytest -q                                   # everything
pytest -q tests/test_projections.py -x      # one module while iterating
python -m torch_vn_algebra.guide --api states   # exact signatures of a module
mkdocs build --strict                       # the website (needs docs_site/requirements.txt)
bash scripts/make_arxiv.sh                  # the arXiv package of the library paper (compiles it twice)
```

## Rules

1. **Compatibility.** CI runs Python 3.8 with torch 1.13, numpy 1.21 and scipy 1.7, and Python 3.8–3.13 with the
   latest versions. So:
   - no `match`, no `X | Y` type unions, no `ast.unparse` in the library;
   - an f-string must not reuse its own quote inside `{...}` (allowed only from Python 3.12);
   - no `torch.func` without a fallback;
   - no complex parameters in `torch.optim.LBFGS` (torch 1.13).
2. **Layout invariants.**
   - Operators are `(batch, C, k_max, k_max)` with zero padding. Decompositions run as one batched call per block
     size, never as a Python loop over sectors.
   - No `.item()`, `bool(tensor)` or `.tolist()` in computational paths: `tests/test_sync_free.py` counts host
     synchronisations.
   - Set a property flag (`is_positive`, `is_projection`, …) only when it is guaranteed mathematically.
3. **Every numerical change gets an independent check:** a dense NumPy/SciPy reference, an exact formula, or brute
   force. Add the test to the module's file in `tests/`.
4. **Numbers in the papers are never typed by hand.**
   - Change the script, rerun it, and commit the regenerated `generated/*.tex` and `.json`.
   - The quick mode of `run_all.sh` writes to `.quick_output/`, never into the paper.
   - Only full runs (`bash scripts/run_all.sh full cuda`) write `man/article/generated/`.
5. **Documentation.**
   - An English chapter of `man/` changes together with its six translations.
   - Inline math is ``$`...`$``, displays use a ```` ```math ```` block, and `\mathrm{Tr}` replaces `\operatorname`.
   - A new public module goes into:
     - the module table of `README.md`;
     - `API_MODULES` in `docs_site/hooks.py`;
     - `MODULES` in `torch_vn_algebra/guide.py`;
     - `llms.txt`;
     - `CHANGELOG.md` under "Unreleased".
   - `tests/test_guide.py` checks the guide, `guide.py` and `llms.txt`.
6. **Claims.**
   - Say "new" only after checking the literature.
   - Cite the source of every known result.
   - Keep statements about the papers' results exactly as strong as the proof or the computation behind them.
7. **arXiv packages are flat:**
   - one self-contained `.tex` file (`scripts/flatten_arxiv.py` inlines the generated numbers and the sections);
   - the figures next to it;
   - no subdirectories.

   Check every package by compiling it in an empty directory (`make_arxiv*.sh` does this).
8. **Releases** follow `RELEASING.md`. The tag `vX.Y.Z` must equal `__version__` in `torch_vn_algebra/__init__.py`:
   `release.yml` checks this before publishing to PyPI.
