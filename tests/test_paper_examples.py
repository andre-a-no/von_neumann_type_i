"""The code listings of the paper (examples/paper) must run."""
import pathlib
import runpy

import pytest

EXAMPLES = sorted((pathlib.Path(__file__).parent.parent / 'examples' / 'paper').glob('ex*.py'))


@pytest.mark.parametrize('path', EXAMPLES, ids=[p.stem for p in EXAMPLES])
def test_paper_example_runs(path, capsys):
    runpy.run_path(str(path), run_name='__main__')
    assert capsys.readouterr().out
