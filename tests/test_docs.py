"""Consistency of docs/physics_for_mathematicians*.md: the translations must follow the English original."""
import re
from pathlib import Path

import pytest

import torch_vn_algebra as vn

DOCS = Path(__file__).resolve().parents[1] / 'docs'
LANGS = ['', 'es', 'fr', 'de', 'zh', 'ja', 'ru']          # order of the language bar
FILES = [DOCS / (f'physics_for_mathematicians{"." + l if l else ""}.md') for l in LANGS]


def _code(text):
    return sorted(set(re.findall(r'`([^`]+)`', text)))


def _sections(text):
    return [line for line in text.splitlines() if line.startswith('## ')]


@pytest.fixture(scope='module')
def english():
    return FILES[0].read_text(encoding='utf-8')


@pytest.mark.parametrize('path', FILES[1:], ids=LANGS[1:])
def test_translation_matches_english(path, english):
    text = path.read_text(encoding='utf-8')
    assert len(_sections(text)) == len(_sections(english))
    assert _code(text) == _code(english), "code identifiers must not be translated, added or dropped"
    assert text.count('\n|') == english.count('\n|'), "tables must have the same rows"


@pytest.mark.parametrize('path', FILES, ids=LANGS[0:1] + LANGS[1:])
def test_language_bar(path):
    bar = path.read_text(encoding='utf-8').splitlines()[0]
    labels = re.findall(r'\[([^\]]+)\]|\*\*([^*]+)\*\*', bar)
    labels = [a or b for a, b in labels]
    assert labels == ['English', 'Español', 'Français', 'Deutsch', '中文', '日本語', 'Русский']
    for target in re.findall(r'\]\(([^)]+)\)', bar):
        assert (DOCS / target).exists()


def test_identifiers_exist(english):
    names = set()
    for mod in (vn, vn.SpinChain, vn.Operator, vn.DensityMatrix, vn.states, vn.channels, vn.dynamics,
                vn.krylov, vn.composite, vn.TypeIAlgebra):
        names |= set(dir(mod))
    for ident in _code(english):
        if ident.endswith('.py') or ident == 'torch_vn_algebra':
            continue
        head, _, tail = ident.partition('.')
        assert (tail or head) in names, ident


def test_every_reference_of_the_paper_is_cited():
    tex = (Path(__file__).resolve().parents[1] / 'paper' / 'v2' / 'main.tex').read_text(encoding='utf-8')
    keys = re.findall(r'\\bibitem\{([^}]+)\}', tex)
    cited = {k.strip() for group in re.findall(r'\\cite\{([^}]+)\}', tex) for k in group.split(',')}
    assert keys and not [k for k in keys if k not in cited]
