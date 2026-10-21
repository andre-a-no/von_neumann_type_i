"""
Flatten man/article/vn_algebra_pytorch.tex into one self-contained .tex file for arXiv.

The paper reads its numbers from generated/*.tex and its sections from sec_*.tex and usecases/*.tex, and falls
back to placeholders when a generated file is missing. On arXiv a missing file would then not stop the build but
silently change the paper, so the upload contains a single .tex file with everything inlined and the fallbacks
resolved: every generated file is inserted in place of the \\gen line, every \\input of a section is replaced by
its text, \\IfFileExists{generated/...}{A}{B} becomes A, \\G is redefined to take its first branch, and \\fig
includes the figure unconditionally (a missing figure is then a hard error). The figures go next to the .tex
file, without subdirectories (arXiv has lost files in subdirectories of earlier uploads).

    python scripts/flatten_arxiv.py SRC_DIR OUT_FILE
"""
import re
import sys
from pathlib import Path


def group(s: str, i: int) -> int:
    """Index just after the balanced {...} group that starts at s[i] == '{'."""
    assert s[i] == '{', s[i:i + 40]
    depth = 0
    j = i
    while j < len(s):
        c = s[j]
        if c == '\\':                                  # a control symbol: \{, \} and \\ are not braces
            j += 2
            continue
        if c == '{':
            depth += 1
        elif c == '}':
            depth -= 1
            if depth == 0:
                return j + 1
        j += 1
    raise ValueError('unbalanced group at ' + s[i:i + 60])


def args(s: str, i: int, n: int):
    """The n brace groups starting at s[i] (whitespace between them skipped): their contents and the end index."""
    out = []
    for _ in range(n):
        while s[i] in ' \t\n':
            i += 1
        j = group(s, i)
        out.append(s[i + 1:j - 1])
        i = j
    return out, i


def replace_macro(s: str, name: str, n: int, fn) -> str:
    """Replace every use of the n-argument macro \\name (not its definition) by fn(arguments)."""
    pat = re.compile(r'\\' + name + r'(?![A-Za-z@])')
    out, pos = [], 0
    for m in pat.finditer(s):
        if m.start() < pos:
            continue
        before = s[max(0, m.start() - 20):m.start()]
        if re.search(r'\\(newcommand|renewcommand|providecommand)\{?$', before):
            continue
        a, end = args(s, m.end(), n)
        r = fn(*a)
        if r is None:
            continue
        out.append(s[pos:m.start()])
        out.append(r)
        pos = end
    out.append(s[pos:])
    return ''.join(out)


def keep_space(a: str) -> str:
    """A branch ending in a control word swallows the space that followed the macro call; {} keeps it."""
    return a + '{}' if re.search(r'\\[A-Za-z@]+$', a) else a


def strip_comment_lines(t: str) -> str:
    return '\n'.join(line for line in t.splitlines() if not line.lstrip().startswith('%'))


def main(src: Path, out: Path):
    gen = src / 'generated'

    def generated(name):
        p = gen / (name + '.tex')
        if not p.exists():
            raise SystemExit(f'missing {p}')
        return p

    def section(path):
        p = src / path
        if not p.exists():
            raise SystemExit(f'missing {p}')
        return expand(p.read_text())

    def ifexists(path, a, b):
        if path.startswith('generated/'):
            generated(path[len('generated/'):-len('.tex')])
            return keep_space(a)
        if path.startswith('figures/'):
            return None                                        # \fig keeps its own test until \fig is redefined
        return a if (src / path).exists() else b

    def expand(t):
        t = replace_macro(t, 'IfFileExists', 3, ifexists)
        t = re.sub(r'\\input\{generated/[A-Za-z_]+\.tex\}', r'\\relax', t)    # already inlined via \gen
        t = replace_macro(t, 'input', 1, lambda p: '%\n' + section(p) + '\n' if not p.startswith('generated/') else None)
        return t

    t = (src / 'vn_algebra_pytorch.tex').read_text()
    # numbers: inline every \gen{...} of the preamble
    m = re.search(r'^(\\gen\{[A-Za-z_]+\})+\s*$', t, re.M)
    names = re.findall(r'\\gen\{([A-Za-z_]+)\}', m.group(0))
    block = ''.join(f'% ---- generated/{n}.tex\n{strip_comment_lines(generated(n).read_text())}\n' for n in names)
    t = t[:m.start()] + block + t[m.end():]
    t = re.sub(r'\\newcommand\{\\gen\}\[1\]\{.*\}\n', '', t)
    # fallbacks resolved: \G takes its first branch, \fig includes the figure
    t = re.sub(r'\\newcommand\{\\G\}\[3\]\{.*\}\n', r'\\newcommand{\\G}[3]{#2}\n', t)
    t = re.sub(r'\\newcommand\{\\fig\}\[2\]\{.*\}\n', r'\\newcommand{\\fig}[2]{\\includegraphics[width=#2]{#1}}\n', t)
    t = expand(t)
    t = re.sub(r'\\graphicspath\{\{figures/\}\}\n', '', t)    # the figures sit next to the .tex file
    for f in re.findall(r'\\G\{([A-Za-z_]+)\}', t):      # \G now takes its first branch; its file must exist
        generated(f)
    out.write_text(t)


if __name__ == '__main__':
    main(Path(sys.argv[1]), Path(sys.argv[2]))
