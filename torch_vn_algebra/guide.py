"""
A guide to the library for AI coding assistants (and for people in a hurry).

    python -m torch_vn_algebra.guide              # conventions, pitfalls and tested recipes (guide.md)
    python -m torch_vn_algebra.guide --api        # every public class, method and function with its signature
    python -m torch_vn_algebra.guide --api states # the same for one module

The guide is the file guide.md next to this module; its code blocks are run by the test suite. The API index is
built at run time from the installed version (inspect; nothing is read from the documentation), so it always matches
the code that will be imported.
"""
import importlib
import inspect
import sys
from pathlib import Path

MODULES = ['algebra', 'functional', 'hilbert_space', 'states', 'channels', 'composite', 'dynamics', 'chains',
           'krylov', 'su2', 'decompose', 'inclusions', 'divergences', 'jordan', 'optimize', 'cost', 'measurement',
           'projections', 'direct_integral', 'anyons', 'learn', 'interop', 'permutational']


def text() -> str:
    """The guide: data layout, conventions, pitfalls, a map from tasks to functions, and tested recipes."""
    return (Path(__file__).with_name('guide.md')).read_text(encoding='utf-8')


def _summary(obj) -> str:
    doc = inspect.getdoc(obj) or ''
    para = doc.split('\n\n', 1)[0]
    return ' '.join(para.split())


def _signature(obj) -> str:
    try:
        return str(inspect.signature(obj))
    except (TypeError, ValueError):
        return '(...)'


def _public(module):
    names = getattr(module, '__all__', None)
    if names is None:
        names = [n for n, v in vars(module).items() if not n.startswith('_') and
                 (inspect.isclass(v) or inspect.isfunction(v)) and getattr(v, '__module__', None) == module.__name__]
    return [(n, getattr(module, n)) for n in names if inspect.isclass(getattr(module, n, None)) or
            inspect.isfunction(getattr(module, n, None))]


def _members(cls):
    out = []
    for n, v in vars(cls).items():
        if n.startswith('_') and n != '__call__':
            continue
        if isinstance(v, (staticmethod, classmethod)):
            out.append((n, v.__func__, '@' + type(v).__name__ + ' '))
        elif isinstance(v, property):
            out.append((n, v.fget, 'property '))
        elif inspect.isfunction(v):
            out.append((n, v, ''))
    return out


def api(module: str = None) -> str:
    """Every public class (with its methods and properties) and function of the library, one line each."""
    lines = []
    for m in ([module] if module else MODULES):
        mod = importlib.import_module(f'torch_vn_algebra.{m}')
        lines.append(f'## torch_vn_algebra.{m}: {_summary(mod)}')
        for name, obj in _public(mod):
            if getattr(obj, '__module__', '').split('.')[-1] != m:
                continue                                         # re-exported from another module
            if inspect.isclass(obj):
                lines.append(f'class {name}{_signature(obj)}: {_summary(obj)}')
                for mn, f, kind in _members(obj):
                    sig = '' if kind == 'property ' else _signature(f)
                    lines.append(f'    {kind}{name}.{mn}{sig}: {_summary(f)}')
                if name == 'TypeIAlgebra':                       # its nested Operator class
                    op = obj.Operator
                    lines.append(f'class Operator (TypeIAlgebra.Operator): {_summary(op)}')
                    for mn, f, kind in _members(op):
                        sig = '' if kind == 'property ' else _signature(f)
                        lines.append(f'    {kind}Operator.{mn}{sig}: {_summary(f)}')
            else:
                lines.append(f'{name}{_signature(obj)}: {_summary(obj)}')
        lines.append('')
    return '\n'.join(lines)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == '--api':
        print(api(argv[1] if len(argv) > 1 else None))
    else:
        print(text())


if __name__ == '__main__':
    main()
