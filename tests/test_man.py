"""
The documentation in man/: the lecture notes (theory, no code) and the tutorial (the library in the context of the
theory; every ```python block is executed, in order, in one namespace per file), each in seven languages, with
cross-links between them.
"""
import re
from pathlib import Path

import pytest

MAN = Path(__file__).resolve().parents[1] / 'man'
SIDES = {'lecture_notes': MAN / 'lecture_notes', 'tutorial': MAN / 'tutorial'}
LANGS = ['es', 'fr', 'de', 'zh', 'ja', 'ru']
LABELS = ['English', 'Español', 'Français', 'Deutsch', '中文', '日本語', 'Русский']


def _chapters(d):
    return sorted(d.glob('[0-9][0-9]_*.md')) + sorted(d.glob('[A-Z]_*.md'))   # chapters, then appendices


def _files(d):
    return [d / 'README.md'] + _chapters(d)


ENGLISH = [(side, p) for side, d in SIDES.items() for p in _files(d)]
ALL = ENGLISH + [(side, p.parent / l / p.name) for side, p in ENGLISH for l in LANGS]


def _blocks(text):
    return re.findall(r'^```python\n(.*?)^```', text, flags=re.S | re.M)


def _headings(text):
    return [line.split(' ', 1)[0] for line in text.splitlines() if line.startswith('#')]


def _anchors(text):
    return re.findall(r'<a id="([^"]+)"></a>', text)


def _id(item):
    side, p = item
    return f'{side}/{p.parent.name + "/" if p.parent.name in LANGS else ""}{p.name}'


@pytest.mark.parametrize('path', _chapters(SIDES['tutorial']), ids=lambda p: p.stem)
def test_tutorial_runs(path, capsys):
    namespace = {'__name__': f'tutorial_{path.stem}'}
    blocks = _blocks(path.read_text(encoding='utf-8'))
    assert blocks, 'a tutorial chapter without code'
    for i, code in enumerate(blocks):
        try:
            exec(compile(code, f'{path.name}[block {i + 1}]', 'exec'), namespace)
        except Exception as e:                               # name the failing block
            raise AssertionError(f'{path.name}, code block {i + 1} failed: {e!r}') from e


def test_lecture_notes_and_tutorial_are_twins():
    names = lambda d: sorted(p.name for p in _chapters(d))
    assert names(SIDES['lecture_notes']) == names(SIDES['tutorial'])


@pytest.mark.parametrize('item', ALL, ids=_id)
def test_no_code_in_lecture_notes_and_no_operatorname(item):
    side, path = item
    text = path.read_text(encoding='utf-8')
    assert '\\operatorname' not in text, 'use \\mathrm{...}: \\operatorname is not rendered by every Markdown engine'
    if side == 'lecture_notes':
        assert '```' not in text, 'the lecture notes contain no code; code goes to the tutorial'


@pytest.mark.parametrize('side', SIDES)
def test_contents_lists_every_chapter(side):
    index = (SIDES[side] / 'README.md').read_text(encoding='utf-8')
    for path in _chapters(SIDES[side]):
        assert f'({path.name})' in index, path.name


@pytest.mark.parametrize('lang', LANGS)
@pytest.mark.parametrize('english', ENGLISH, ids=_id)
def test_translation_follows_english(lang, english):
    """Translations keep the code byte for byte, the sections, tables, displays and anchors."""
    side, path = english
    tr_path = path.parent / lang / path.name
    assert tr_path.exists(), f'missing translation {lang}/{path.name}'
    en, tr = path.read_text(encoding='utf-8'), tr_path.read_text(encoding='utf-8')
    assert _blocks(tr) == _blocks(en), 'code blocks must be identical to the English original'
    assert _headings(tr) == _headings(en), 'same sections'
    assert _anchors(tr) == _anchors(en), 'same anchors'
    assert tr.count('\n|') == en.count('\n|'), 'same table rows'
    assert tr.count('$$') == en.count('$$'), 'same displayed formulas'
    assert tr.count('<details>') == en.count('<details>')


@pytest.mark.parametrize('item', ALL, ids=_id)
def test_language_bar(item):
    _, path = item
    bar = path.read_text(encoding='utf-8').splitlines()[0]
    labels = [a or b for a, b in re.findall(r'\[([^\]]+)\]|\*\*([^*]+)\*\*', bar)]
    assert labels == LABELS
    for target in re.findall(r'\]\(([^)]+)\)', bar):
        assert (path.parent / target).resolve().exists(), target


@pytest.mark.parametrize('item', ALL, ids=_id)
def test_links_resolve(item):
    """Every relative link points to an existing file, and every #fragment to an explicit anchor there."""
    _, path = item
    text = path.read_text(encoding='utf-8')
    text = re.sub(r'^```.*?^```', '', text, flags=re.S | re.M)
    bad = []
    for target in re.findall(r'\]\(([^)\s]+)\)', text):
        if re.match(r'[a-z]+:', target):
            continue
        file, _, frag = target.partition('#')
        dest = (path.parent / file).resolve() if file else path
        if not dest.exists():
            bad.append(target)
        elif frag and dest.suffix == '.md' and f'<a id="{frag}"></a>' not in dest.read_text(encoding='utf-8'):
            bad.append(target)
    assert not bad, f'broken links: {bad}'
