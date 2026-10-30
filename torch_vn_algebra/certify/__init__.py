r"""
Certified envelopes: a numerical envelope V(delta) = sup F becomes a rigorous interval [L, U] containing V(delta),
with a certificate file that a separate small program checks (:mod:`.verify`, which depends on python-flint only).

    from torch_vn_algebra import certify
    cert = certify.certify("holder_projections", delta="1/10")      # item 14 of the trace catalogue
    cert.lower, cert.upper                                          # exact rationals, L <= V(delta) <= U
    cert.save("certificates/holder_projections_alpha3_delta_1_10.json.gz")

    $ python -m torch_vn_algebra.certify.verify certificates/holder_projections_alpha3_delta_1_10.json.gz

Lower bound L: the numerical maximiser is replaced by an exact rational configuration of the original problem
(Cayley transforms, rational spectra; :mod:`.witness`), and the violation is checked in exact rational arithmetic.
Upper bound U: branch and bound on the reduced variables (:mod:`.engine`), with Arb bounds on every box: the smaller
of the natural interval extension and the mean-value form, plus a monotonicity test. The reductions (two-projection
theorem, compression, Hölder, phases, ...) are listed in the certificate with references; they are proved in the text,
not by the verifier.

Requires python-flint (``pip install torch_vn_algebra[certify]``). Not imported by ``import torch_vn_algebra``.

Not implemented yet: the SDP layer (trace-polynomial hierarchies of Klep, Magron and Volčič with rational rounding of
Peyrl and Parrilo) for bounds over all finite von Neumann algebras, and Lipschitz box bounds for matrix functions
(|X|, X^{1/2}, e^X, the spectral geometric mean) via Kittaneh-type inequalities; see ROADMAP_certify.md.
"""
import datetime
import gzip
import json
import platform
import time

from flint import fmpq, ctx, __version__ as _flint_version

from .problems import PROBLEMS, Problem, HolderProjections, _round_up
from .engine import upper_bound
from .witness import parse_rational, fmt_rational

FORMAT = "torch_vn_algebra-certificate/1"
__all__ = ["certify", "Certificate", "PROBLEMS", "Problem", "HolderProjections", "load", "FORMAT"]


class Certificate:
    """A certified interval [lower, upper] for V(delta), with everything the verifier needs."""

    def __init__(self, data):
        self.data = data

    @property
    def lower(self):
        return parse_rational(self.data["claim"]["lower"])

    @property
    def upper(self):
        return parse_rational(self.data["claim"]["upper"])

    @property
    def delta(self):
        return parse_rational(self.data["problem"]["delta"])

    def __repr__(self):
        return (f"Certificate({self.data['problem']['id']}, delta={self.data['problem']['delta']}: "
                f"[{float(self.lower):.12g}, {float(self.upper):.12g}], {self.data['upper']['leaves']} boxes)")

    def save(self, path):
        raw = json.dumps(self.data, indent=1).encode()
        if str(path).endswith(".gz"):
            with gzip.GzipFile(path, "wb", mtime=0) as f:
                f.write(raw)
        else:
            with open(path, "wb") as f:
                f.write(raw)
        return path


def load(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as f:
        return Certificate(json.loads(f.read().decode()))


def _problem(problem, params):
    if isinstance(problem, Problem):
        return problem
    return PROBLEMS[problem](**params)


def certify(problem, delta, rel_tol=1e-6, prec=128, starts=40, max_leaves=2_000_000, **params):
    """
    Certified interval for V(delta) of `problem` (a key of PROBLEMS or a Problem instance; `params` go to its
    constructor). `delta` is converted to an exact rational (pass a string such as "1/10" to avoid float input).
    The target width is U - L = rel_tol * L (rel_tol > 0). Returns a :class:`Certificate`.
    """
    P = _problem(problem, params)
    delta = parse_rational(delta) if isinstance(delta, str) else fmpq(*_as_ratio(delta))
    if not (0 < delta < 1):
        raise ValueError("delta must lie in (0, 1)")
    t0 = time.time()
    old = ctx.prec
    ctx.prec = prec
    try:
        point, fval = P.search(delta, starts=starts)
        wit, L = P.witness(point, delta)
        U = _round_up(float(L) * (1 + rel_tol))
        F = lambda x: P.objective(x, delta)  # noqa: E731
        bb = upper_bound(F, P.domain(delta), U, prec=prec, max_leaves=max_leaves)
    finally:
        ctx.prec = old
    data = dict(
        format=FORMAT,
        problem=dict(P.describe(), delta=fmt_rational(delta)),
        claim=dict(lower=fmt_rational(L), upper=fmt_rational(U), lower_float=float(L), upper_float=float(U),
                   search_value=fval),
        witness=wit,
        upper=dict(method="branch and bound, Arb balls: min(natural extension, mean-value form), monotonicity test",
                   precision_bits=prec, domain=[[fmt_rational(lo), fmt_rational(hi)] for lo, hi in P.domain(delta)],
                   leaves=bb.leaves, max_depth=bb.max_depth, tree=" ".join(bb.tokens)),
        software=dict(python_flint=_flint_version, python=platform.python_version(),
                      torch_vn_algebra=_package_version()),
        created=datetime.date.today().isoformat(), seconds=round(time.time() - t0, 2),
    )
    return Certificate(data)


def _as_ratio(x):
    from fractions import Fraction
    f = Fraction(x).limit_denominator(10 ** 6) if isinstance(x, float) else Fraction(x)
    return f.numerator, f.denominator


def _package_version():
    try:
        from .. import __version__
        return __version__
    except Exception:  # pragma: no cover
        return None
