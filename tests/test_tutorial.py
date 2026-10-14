"""Every ```python block of the tutorial is executed, in order, in one namespace per chapter."""
import re
from pathlib import Path

import pytest

TUTORIAL = Path(__file__).resolve().parents[1] / 'docs' / 'tutorial'
CHAPTERS = sorted(TUTORIAL.glob('[0-9][0-9]_*.md')) + sorted(TUTORIAL.glob('[A-Z]_*.md'))   # chapters, then appendices


def _blocks(text):
    return re.findall(r'^```python\n(.*?)^```', text, flags=re.S | re.M)


@pytest.mark.parametrize('path', CHAPTERS, ids=[p.stem for p in CHAPTERS])
def test_chapter_runs(path, capsys):
    namespace = {'__name__': f'tutorial_{path.stem}'}
    blocks = _blocks(path.read_text(encoding='utf-8'))
    assert blocks, 'a chapter without code'
    for i, code in enumerate(blocks):
        try:
            exec(compile(code, f'{path.name}[block {i + 1}]', 'exec'), namespace)
        except Exception as e:                               # name the failing block
            raise AssertionError(f'{path.name}, code block {i + 1} failed: {e!r}') from e


def test_contents_lists_every_chapter():
    index = (TUTORIAL / 'README.md').read_text(encoding='utf-8')
    for path in CHAPTERS:
        assert f'({path.name})' in index, path.name


LANGS = ['es', 'fr', 'de', 'zh', 'ja', 'ru']
LABELS = ['English', 'Español', 'Français', 'Deutsch', '中文', '日本語', 'Русский']
FILES = [TUTORIAL / 'README.md'] + CHAPTERS


def _headings(text):
    return [line.split(' ', 1)[0] for line in text.splitlines() if line.startswith('#')]


@pytest.mark.parametrize('lang', LANGS)
@pytest.mark.parametrize('english', FILES, ids=[p.name for p in FILES])
def test_translation_follows_english(lang, english):
    """Translations keep the code byte for byte, the section structure, the tables and the displays."""
    if not (TUTORIAL / lang).is_dir():
        pytest.skip(f'no {lang} translation yet')
    path = TUTORIAL / lang / english.name
    assert path.exists(), f'missing translation {lang}/{english.name}'
    en, tr = english.read_text(encoding='utf-8'), path.read_text(encoding='utf-8')
    assert _blocks(tr) == _blocks(en), 'code blocks must be identical to the English original'
    assert _headings(tr) == _headings(en), 'same sections'
    assert tr.count('\n|') == en.count('\n|'), 'same table rows'
    assert tr.count('$$') == en.count('$$'), 'same displayed formulas'
    assert tr.count('<details>') == en.count('<details>')


@pytest.mark.parametrize('path', FILES + [TUTORIAL / l / p.name for l in LANGS for p in FILES],
                         ids=lambda p: f'{p.parent.name}/{p.name}')
def test_language_bar(path):
    if not path.exists():
        pytest.skip('translation missing (reported by test_translation_follows_english)')
    if path.parent == TUTORIAL:                                  # English: links only to existing translations
        pytest.skip('checked once the translations exist') if not all((TUTORIAL / l).is_dir() for l in LANGS) else None
    bar = path.read_text(encoding='utf-8').splitlines()[0]
    labels = [a or b for a, b in re.findall(r'\[([^\]]+)\]|\*\*([^*]+)\*\*', bar)]
    assert labels == LABELS
    for target in re.findall(r'\]\(([^)]+)\)', bar):
        assert (path.parent / target).resolve().exists(), target
