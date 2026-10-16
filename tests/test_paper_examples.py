"""The code listings of the article (man/article/listings) must run."""
import pathlib
import runpy

import pytest

LISTINGS = pathlib.Path(__file__).parent.parent / 'man' / 'article' / 'listings'
# connes_stormer.py, permutational.py and symmetry_resolved.py are run by the tests of their use cases
EXAMPLES = sorted(LISTINGS.glob('ex*.py')) + [LISTINGS / f for f in ('inclusions.py', 'symmetry.py',
                                                                      'trace_characterizations.py')]


@pytest.mark.parametrize('path', EXAMPLES, ids=[p.stem for p in EXAMPLES])
def test_paper_example_runs(path, capsys):
    runpy.run_path(str(path), run_name='__main__')
    assert capsys.readouterr().out
