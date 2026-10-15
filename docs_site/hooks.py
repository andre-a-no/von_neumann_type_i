"""
MkDocs hooks for the documentation website (mkdocs.yml).

The Markdown files of man/ are written for GitHub (inline math $`...`$, displays in ```math blocks, <details> blocks,
links between the files). This hook makes them a website without changing them:

* the home page is README.md, with its links rewritten (files under man/ become pages, other files GitHub links);
* the navigation lists the lecture notes and the tutorial in English and the six translations;
* $`...`$ becomes \\(...\\) for pymdownx.arithmatex (```math blocks are a superfences custom fence in mkdocs.yml),
  outside code blocks and code spans, so MathJax renders exactly the formulas that GitHub renders;
* <details> blocks get markdown="1" (md_in_html), as GitHub renders Markdown inside them;
* the API reference is generated from the source (ast; nothing is imported), one page per module: signatures and
  docstrings verbatim (the docstrings are plain text with ASCII formulas, which Markdown would mangle);
* every page gets a meta description (its first paragraph unless it sets one), keywords, Open Graph tags, its
  language in <html lang>, and hreflang links to its translations; the home page gets a schema.org JSON-LD record;
* robots.txt (with the sitemap that MkDocs writes) and llms.txt are added at the site root.
"""
import ast
import html
import inspect
import json
import posixpath
import re
from pathlib import Path

from mkdocs.structure.files import File

ROOT = Path(__file__).resolve().parents[1]
REPO = 'https://github.com/andre-a-no/von_neumann_type_i'
BLOB = REPO + '/blob/main/'
LANGS = ['en', 'es', 'fr', 'de', 'zh', 'ja', 'ru']
NATIVE = {'en': 'English', 'es': 'Español', 'fr': 'Français', 'de': 'Deutsch', 'zh': '中文', 'ja': '日本語',
          'ru': 'Русский'}
SIDES = ['lecture_notes', 'tutorial']
SIDE_NAME = {'lecture_notes': {'en': 'Lecture notes', 'es': 'Apuntes de clase', 'fr': 'Notes de cours',
                               'de': 'Vorlesungsnotizen', 'zh': '讲义', 'ja': '講義ノート', 'ru': 'Конспект лекций'},
             'tutorial': {'en': 'Tutorial', 'es': 'Tutorial', 'fr': 'Tutoriel', 'de': 'Tutorial', 'zh': '教程',
                          'ja': 'チュートリアル', 'ru': 'Руководство'}}
HOME_TITLE = 'Block-diagonal operators, symmetry sectors and quantum channels in PyTorch'
API_MODULES = ['algebra', 'functional', 'hilbert_space', 'states', 'channels', 'composite', 'dynamics', 'chains',
               'krylov', 'su2', 'decompose', 'inclusions', 'divergences', 'jordan', 'optimize', 'cost']

_files = {}


# ---------------------------------------------------------------------------------------------------- navigation

def _chapters(d: Path):
    return sorted(p.name for p in d.glob('[0-9][0-9]_*.md')) + sorted(p.name for p in d.glob('[A-Z]_*.md'))


def _side_nav(side, lang):
    d = ROOT / 'man' / side / ('' if lang == 'en' else lang)
    prefix = f'{side}/' + ('' if lang == 'en' else f'{lang}/')
    return [prefix + 'README.md'] + [prefix + name for name in _chapters(d)]


def on_config(config):
    nav = [{'Home': 'index.md'}]
    nav += [{SIDE_NAME[s]['en']: _side_nav(s, 'en')} for s in SIDES]
    nav.append({'API reference': ['api/index.md'] + [f'api/{m}.md' for m in API_MODULES]})
    nav.append({'Paper (PDF)': 'article/main.pdf'})
    nav.append({'Translations': [{NATIVE[l]: [{SIDE_NAME[s][l]: _side_nav(s, l)} for s in SIDES]}
                                 for l in LANGS[1:]]})
    config['nav'] = nav
    return config


# ------------------------------------------------------------------------------------------------- generated files

def _readme_link(m):
    text, target = m.group(1), m.group(2)
    if re.match(r'[a-z]+:|#', target):
        return m.group(0)
    path, _, frag = target.partition('#')
    if path.startswith('man/') and (path.endswith('.md') or path == 'man/article/main.pdf'):
        new = path[len('man/'):]
    else:
        new = BLOB + path.rstrip('/')
    return f'[{text}]({new}{"#" + frag if frag else ""})'


def _home_page():
    text = (ROOT / 'README.md').read_text(encoding='utf-8')
    # links outside code blocks; [text](target) with text that may itself be an image
    parts = re.split(r'(^```.*?^```)', text, flags=re.S | re.M)
    for i in range(0, len(parts), 2):
        parts[i] = re.sub(r'\[((?:!\[[^\]]*\]\([^)]*\))|[^\]]*)\]\(([^)\s]+)\)', _readme_link, parts[i])
    meta = ('---\n'
            f'title: "{HOME_TITLE}"\n'
            'description: "torch_vn_algebra is a PyTorch library for block-diagonal matrices with symmetry sectors '
            '(Type I von Neumann algebras): random unitaries, quantum channels, Lindblad dynamics, spin chains, '
            'on CPU and GPU."\n'
            '---\n\n')
    return meta + ''.join(parts)


def _api_index():
    esc = lambda t: re.sub(r'([\\`*_{}\[\]|])', r'\\\1', html.escape(t, quote=False))
    rows = '\n'.join(f'| [`torch_vn_algebra.{m}`]({m}.md) | {esc(_module_summary(m))} |' for m in API_MODULES)
    return ('---\ndescription: "API reference of torch_vn_algebra, generated from the docstrings: algebras, '
            'operators, states, channels, dynamics, spin chains, decompositions, divergences."\n---\n\n'
            '# API reference\n\n'
            'Generated from the docstrings of the package. The package exports the main classes at the top level '
            '(`from torch_vn_algebra import TypeIAlgebra, SpinChain, SU2Chain, ...`); the modules are:\n\n'
            '| module | contents |\n|---|---|\n' + rows + '\n\n'
            'The [tutorial](../tutorial/README.md) shows every module at work, with tested code.\n')


def _module_summary(m):
    tree = ast.parse((ROOT / 'torch_vn_algebra' / f'{m}.py').read_text(encoding='utf-8'))
    doc = ast.get_docstring(tree)
    if not doc:                                   # no module docstring: the first public class or function
        doc = next((ast.get_docstring(n) for n in _public(tree) if ast.get_docstring(n)), '')
    first = doc.strip().split('\n\n')[0]
    first = re.sub(r':math:`([^`]*)`', r'\1', ' '.join(first.split()))
    first = re.sub(r'``([^`]*)``', r'\1', first)
    return first


def _public(node, names=None):
    for n in node.body:
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if (names is not None and n.name in names) or (names is None and not n.name.startswith('_')):
                yield n


def _signature(n):
    if isinstance(n, ast.ClassDef):
        init = next((b for b in n.body if isinstance(b, ast.FunctionDef) and b.name == '__init__'), None)
        args = ast.unparse(init.args) if init else ''
        args = re.sub(r'^self(, )?', '', args)
        return f'class {n.name}({args})'
    args = ast.unparse(n.args)
    ret = f' -> {ast.unparse(n.returns)}' if n.returns else ''
    deco = ''.join(f'@{ast.unparse(d)}\n' for d in n.decorator_list if ast.unparse(d) in
                   ('property', 'staticmethod', 'classmethod', 'cached_property', 'functools.cached_property'))
    return f'{deco}def {n.name}({args}){ret}'


def _block(text, lang='text'):
    return f'````{lang}\n{text}\n````\n\n'


def _api_entries(node, prefix, level):
    out = []
    for n in _public(node):
        name = f'{prefix}{n.name}'
        out.append(f'{"#" * level} {name}\n\n' + _block(_signature(n), 'python'))
        doc = ast.get_docstring(n)
        if doc:
            out.append(_block(doc))
        if isinstance(n, ast.ClassDef) and level < 4:
            out += _api_entries(n, name + '.', level + 1)
    return out


def _api_page(m):
    tree = ast.parse((ROOT / 'torch_vn_algebra' / f'{m}.py').read_text(encoding='utf-8'))
    names = None
    for n in tree.body:
        if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == '__all__' for t in n.targets):
            names = set(ast.literal_eval(n.value))
    summary = _module_summary(m)
    desc = f'torch_vn_algebra.{m}: ' + summary.replace('"', "'")
    desc = desc if len(desc) <= 160 else desc[:157].rsplit(' ', 1)[0] + '…'
    out = [f'---\ntitle: "torch_vn_algebra.{m}"\ndescription: "{desc}"\n---\n\n'
           f'# torch_vn_algebra.{m}\n\n'
           f'Source: [`torch_vn_algebra/{m}.py`]({BLOB}torch_vn_algebra/{m}.py). '
           'Signatures and docstrings as in the source.\n\n']
    doc = ast.get_docstring(tree)
    if doc:
        out.append('## Module documentation\n\n' + _block(doc))
    top = [n for n in _public(tree, names)]
    holder = ast.Module(body=top, type_ignores=[])
    out += _api_entries(holder, '', 2)
    return ''.join(out)


def _robots(config):
    return f'User-agent: *\nAllow: /\n\nSitemap: {config["site_url"]}sitemap.xml\n'


def on_files(files, config):
    def add(uri, content):
        files.append(File.generated(config, uri, content=content))

    add('index.md', _home_page())
    add('api/index.md', _api_index())
    for m in API_MODULES:
        add(f'api/{m}.md', _api_page(m))
    add('robots.txt', _robots(config))
    llms = ROOT / 'llms.txt'
    if llms.exists():
        add('llms.txt', llms.read_text(encoding='utf-8'))
    _files.clear()
    _files.update({f.src_uri: f for f in files})
    return files


# ------------------------------------------------------------------------------------------------------ markdown

_INLINE = re.compile(r'(?P<math>\$`(?P<tex>[^`]+)`\$)|(?P<code>(?P<ticks>`+)[^`].*?(?P=ticks))', re.S)
_FENCE = re.compile(r'^(\s*)(`{3,}|~{3,})(.*)$')


def _convert_inline(segment, table_row=False):
    def sub(m):
        if m.group('code'):
            return m.group(0)
        tex = m.group('tex')
        if table_row:
            tex = tex.replace('\\|', '\\Vert ').replace('|', '\\vert ')
        return '\\(' + tex + '\\)'
    return _INLINE.sub(sub, segment)


def _convert(markdown):
    """$`...`$ -> \\(...\\) outside fenced blocks and code spans; <details> -> <details markdown="1">."""
    out, prose, fence = [], [], None

    def flush():
        if prose:
            text = '\n'.join(prose)
            # inline code and math may not span table rows; convert table rows separately
            lines = text.split('\n')
            if any(l.lstrip().startswith('|') for l in lines):
                text = '\n'.join(_convert_inline(l, table_row=l.lstrip().startswith('|')) for l in lines)
            else:
                text = _convert_inline(text)
            out.append(text.replace('<details>', '<details markdown="1">'))
            prose.clear()

    for line in markdown.split('\n'):
        m = _FENCE.match(line)
        if fence is None and m:
            flush()
            fence = m.group(2)[0] * len(m.group(2))
            out.append(line)
        elif fence is not None:
            out.append(line)
            if m and m.group(2).startswith(fence) and not m.group(3).strip():
                fence = None
        else:
            prose.append(line)
    flush()
    return '\n'.join(out)


def _plain(text):
    text = re.sub(r'\$`([^`]+)`\$', lambda m: m.group(1) if '\\' not in m.group(1) else '', text)
    text = re.sub(r'!?\[([^\]]*)\]\([^)]*\)', r'\1', text)
    text = re.sub(r'<[^>]+>', '', text)
    text = text.replace('`', '')
    text = re.sub(r'\*\*(.+?)\*\*', r'\1', text)
    text = re.sub(r'(?<![\w*])\*(?![\s-])(.+?)(?<!\s)\*(?![\w*])', r'\1', text)   # emphasis, not '*-algebra'
    text = re.sub(r'\(\s*\)', '', text)
    text = ' '.join(text.split())
    cjk = '\u3000-\u30ff\u3400-\u9fff\uff00-\uffef'
    return re.sub(f'(?<=[{cjk}]) (?=[{cjk}])', '', text)


def _description(markdown):
    body = re.sub(r'^```.*?^```', '', markdown, flags=re.S | re.M)
    for para in re.split(r'\n\s*\n', body):
        p = para.strip()
        if not p or p[0] in '#|<[*-!>' or p.startswith('**English**'):
            continue
        text = _plain(p)
        if len(text) < 60:
            continue
        if len(text) > 160:
            cut = text[:157]
            cut = cut.rsplit(' ', 1)[0] if ' ' in cut[100:] else cut
            text = cut.rstrip(' ,;:—–-') + '…'
        return text
    return None


def on_page_markdown(markdown, page, config, files):
    if 'title' not in page.meta:            # the H1 (MkDocs looks only at the first line, the language bar)
        h1 = re.search(r'^# (.+)$', re.sub(r'^```.*?^```', '', markdown, flags=re.S | re.M), flags=re.M)
        if h1:
            page.meta['title'] = _plain(h1.group(1))
    if 'description' not in page.meta:
        desc = _description(markdown)
        if desc:
            page.meta['description'] = desc
    return _convert(markdown)


def _strip_math(title):
    return re.sub(r'\\\((.*?)\\\)', r'\1', title) if title else title


def on_page_content(content, page, config, files):
    def walk(items):
        for item in items:
            item.title = _strip_math(item.title)
            walk(item.children)
    walk(page.toc)
    return content


# ---------------------------------------------------------------------------------------------------------- head

def _lang_and_key(src_uri):
    parts = src_uri.split('/')
    if len(parts) == 3 and parts[0] in SIDES and parts[1] in LANGS:
        return parts[1], (parts[0], parts[2])
    if len(parts) == 2 and parts[0] in SIDES:
        return 'en', (parts[0], parts[1])
    return 'en', None


def on_post_page(output, page, config):
    src = page.file.src_uri
    lang, key = _lang_and_key(src)
    site = config['site_url']
    head = []
    if key:
        side, name = key
        for l in LANGS:
            f = _files.get(f'{side}/{name}' if l == 'en' else f'{side}/{l}/{name}')
            if f:
                head.append(f'<link rel="alternate" hreflang="{l}" href="{site}{f.url}">')
                if l == 'en':
                    head.append(f'<link rel="alternate" hreflang="x-default" href="{site}{f.url}">')
    keywords = ' '.join(config['extra'].get('keywords', '').split())
    if keywords:
        head.append(f'<meta name="keywords" content="{html.escape(keywords)}">')
    title = (page.title + ' - ' if page.title else '') + config['site_name']
    if key and key[1] != 'README.md' and key[0] == 'lecture_notes':
        title = f'{page.title} - {SIDE_NAME[key[0]][lang]} - {config["site_name"]}'
        output = re.sub(r'<title>.*?</title>', f'<title>{html.escape(title)}</title>', output, count=1, flags=re.S)
    desc = page.meta.get('description') or config['site_description']
    url = page.canonical_url or site
    head += [f'<meta property="og:type" content="{"website" if page.is_homepage else "article"}">',
             f'<meta property="og:site_name" content="{html.escape(config["site_name"])}">',
             f'<meta property="og:title" content="{html.escape(_strip_math(title))}">',
             f'<meta property="og:description" content="{html.escape(desc)}">',
             f'<meta property="og:url" content="{url}">',
             f'<meta property="og:locale" content="{lang}">',
             '<meta name="twitter:card" content="summary">']
    if page.is_homepage:
        record = {
            '@context': 'https://schema.org', '@type': 'SoftwareSourceCode',
            'name': 'torch_vn_algebra', 'description': desc, 'url': site,
            'codeRepository': REPO, 'programmingLanguage': 'Python', 'runtimePlatform': 'PyTorch',
            'license': 'https://opensource.org/licenses/MIT',
            'keywords': [k.strip() for k in keywords.split(',') if k.strip()],
            'author': [{'@type': 'Person', 'name': 'Irina Nikolaeva'}, {'@type': 'Person', 'name': 'Andrej Novikov'}],
            'citation': config['extra'].get('paper'),
        }
        head.append('<script type="application/ld+json">' + json.dumps(record, ensure_ascii=False) + '</script>')
    output = output.replace('</head>', '\n'.join(head) + '\n</head>', 1)
    if lang != 'en':
        output = re.sub(r'<html lang="[^"]*"', f'<html lang="{lang}"', output, count=1)
    return output
