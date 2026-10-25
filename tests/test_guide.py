"""
The guide for AI coding assistants (torch_vn_algebra/guide.md, printed by python -m torch_vn_algebra.guide): its code
blocks run in order in one namespace, every name it mentions exists, the API index covers the public API, and
llms.txt lists every module.
"""
import re
from pathlib import Path

import pytest

import torch_vn_algebra
from torch_vn_algebra import guide

ROOT = Path(__file__).resolve().parents[1]


def test_guide_code_runs():
    namespace = {'__name__': 'guide'}
    blocks = re.findall(r'^```python\n(.*?)^```', guide.text(), flags=re.S | re.M)
    assert len(blocks) >= 5
    for i, code in enumerate(blocks):
        try:
            exec(compile(code, f'guide.md[block {i + 1}]', 'exec'), namespace)
        except Exception as e:
            raise AssertionError(f'guide.md, code block {i + 1} failed: {e!r}') from e


def test_api_index_covers_the_public_api():
    text = guide.api()
    for name in torch_vn_algebra.__all__:
        obj = getattr(torch_vn_algebra, name)
        if callable(obj) and name != 'Operator':
            assert re.search(rf'^(class )?{name}\(', text, flags=re.M), name
    assert 'Operator.apply_function(' in text and 'class Operator' in text


def test_modules_are_listed_everywhere():
    modules = {p.stem for p in (ROOT / 'torch_vn_algebra').glob('*.py')} - {'__init__', 'guide'}
    assert set(guide.MODULES) == modules
    hooks = (ROOT / 'docs_site' / 'hooks.py').read_text(encoding='utf-8')
    api_modules = re.search(r'^API_MODULES = (\[.*?\])', hooks, flags=re.S | re.M).group(1)
    assert guide.MODULES == eval(api_modules)                    # the website documents the same modules
    llms = (ROOT / 'llms.txt').read_text(encoding='utf-8')
    for m in modules:
        assert f'torch_vn_algebra.{m}' in llms, m


def test_names_in_the_table_exist():
    """Every `module.name` and `Class.method` in backticks of the guide is an attribute of the library."""
    text = guide.text()
    for mod, name in set(re.findall(r'`(?:torch_vn_algebra\.)?([a-z_]+)\.([A-Za-z_]+)', text)):
        if mod in guide.MODULES:
            module = __import__(f'torch_vn_algebra.{mod}', fromlist=['_'])
            assert hasattr(module, name), f'{mod}.{name}'


def test_cli(capsys):
    guide.main([])
    assert 'guide for AI coding assistants' in capsys.readouterr().out
    guide.main(['--api', 'states'])
    out = capsys.readouterr().out
    assert 'random_density_matrix(' in out and 'torch_vn_algebra.channels' not in out
