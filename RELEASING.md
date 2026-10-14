# Releasing torch_vn_algebra

## One-time set-up (repository owner, about 10 minutes)

1. **PyPI trusted publishing** (no token in the repository):
   - create an account on <https://pypi.org> (with two-factor authentication);
   - *Your account → Publishing → Add a new pending publisher*: PyPI project name `torch_vn_algebra`,
     owner `andre-a-no`, repository `von_neumann_type_i`, workflow `release.yml`, environment `pypi`;
   - on GitHub: *Settings → Environments → New environment* `pypi` (optionally restricted to tags `v*`).
2. **Zenodo DOI for every release**: log in at <https://zenodo.org> with GitHub, *GitHub* tab → switch on
   `andre-a-no/von_neumann_type_i`. Zenodo reads `.zenodo.json` and archives every published GitHub release.

## Every release

1. All tests pass (`pytest -q`); for a release that changes numerical code, the GPU check too
   (`python scripts/gpu_check.py`).
2. Set the version in `torch_vn_algebra/__init__.py` (`__version__`) and in `CITATION.cff` (`version`);
   move the "Unreleased" entries of `CHANGELOG.md` under the new version with today's date.
3. Commit, merge into `main`.
4. On GitHub: *Releases → Draft a new release*, tag `vX.Y.Z` on `main`, title `vX.Y.Z`, notes = the CHANGELOG
   section → *Publish release*.
5. The `release` workflow builds the sdist and wheel, checks that the tag equals the version and uploads to PyPI;
   Zenodo archives the release and mints a DOI. Add the DOI badge to `README.md` and the DOI to `CITATION.cff`
   (`doi:` field) after the first release.

## Local check before tagging

```bash
pip install build twine
python -m build && twine check --strict dist/*
```
