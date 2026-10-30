"""
Numerical statements of the papers against their certificates (certificates/claims.json).

* 'rounded' claims: the printed number must be what both ends of the certified interval [L, U] round to;
* 'lower' claims (search values): every printed value must be <= U.

The certificates themselves are checked by the independent verifier: the δ = 9/10 certificate always (a few
seconds), all of them with VN_VERIFY_CERTIFICATES=1 (several minutes; δ = 1/10 is also verified in test_certify.py).
"""
import gzip
import importlib.util
import json
import os
import re
from decimal import Decimal, ROUND_HALF_EVEN
from fractions import Fraction

import pytest

pytest.importorskip("flint")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLAIMS = json.load(open(os.path.join(ROOT, "certificates", "claims.json")))["claims"]
ALL = bool(os.environ.get("VN_VERIFY_CERTIFICATES"))


def _cert(rel):
    with gzip.open(os.path.join(ROOT, rel), "rb") as f:
        return json.loads(f.read().decode())


def _interval(c):
    return Fraction(c["claim"]["lower"]), Fraction(c["claim"]["upper"])


def _printed(s):
    """'2.8\\times10^{-4}' or '0.0086' -> (Decimal value, quantum of the last printed digit)."""
    m = re.fullmatch(r"\s*([0-9.]+)\s*(?:\\times\s*10\^\{(-?\d+)\})?\s*", s)
    if not m:
        raise ValueError(f"cannot parse printed number {s!r}")
    mant, exp = Decimal(m.group(1)), int(m.group(2) or 0)
    return mant.scaleb(exp), Decimal(1).scaleb(mant.as_tuple().exponent + exp)


def _macro(src, name):
    text = open(os.path.join(ROOT, src)).read()
    m = re.search(r"\\newcommand\{\\" + name + r"\}\{(.*)\}\s*$", text, re.M)
    assert m, f"macro {name} not found in {src}"
    return m.group(1)


ROUNDED = [c for c in CLAIMS if c["kind"] == "rounded"]
LOWER = [c for c in CLAIMS if c["kind"] == "lower"]


@pytest.mark.parametrize("claim", ROUNDED, ids=[c["id"] for c in ROUNDED])
def test_rounded_claim_matches_the_certified_interval(claim):
    printed = _macro(claim["source"], claim["macro"])
    assert printed == claim["printed"], "the paper changed: update certificates/claims.json"
    value, q = _printed(printed)
    L, U = _interval(_cert(claim["certificate"]))
    for end in (L, U):
        d = (Decimal(end.numerator) / Decimal(end.denominator)).quantize(q, rounding=ROUND_HALF_EVEN)
        assert d == value, f"{claim['id']}: {float(end)} does not round to {printed}"


@pytest.mark.parametrize("claim", LOWER, ids=[c["id"] for c in LOWER])
def test_search_values_are_below_the_certified_upper_bound(claim):
    data = json.load(open(os.path.join(ROOT, claim["source"])))
    grid = data
    for k in claim["grid"]:
        grid = grid[k]
    node = data
    for k in claim["path"][:-2]:                     # ... "dims", "*" (every n), "certified"
        node = node[k]
    checked = 0
    assert claim["path"][-2] == "*"
    for dim, rec in node.items():
        for d, v in zip(grid, rec[claim["path"][-1]]):
            if d == 0:
                continue
            fr = Fraction(d).limit_denominator(1000)
            L, U = _interval(_cert(claim["certificate"].format(p=fr.numerator, q=fr.denominator)))
            assert v <= U, f"{claim['id']}: n = {dim}, delta = {fr}: {v} > U = {float(U)}"
            checked += 1
    assert checked > 0


def _verifier():
    path = os.path.join(ROOT, "torch_vn_algebra", "certify", "verify.py")
    spec = importlib.util.spec_from_file_location("tvna_verify_standalone_claims", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CERTS = sorted(f for f in os.listdir(os.path.join(ROOT, "certificates")) if f.endswith(".json.gz"))


@pytest.mark.parametrize("name", CERTS)
def test_certificate_verifies(name):
    if not ALL and name != "holder_projections_alpha3_delta_9_10.json.gz":
        pytest.skip("set VN_VERIFY_CERTIFICATES=1 to verify every certificate")
    L, U = _verifier().verify(os.path.join(ROOT, "certificates", name), log=lambda *a: None)
    assert 0 < L <= U and (U - L) <= 1e-6 * L * 1.01


def test_summary_matches_the_certificates():
    summ = json.load(open(os.path.join(ROOT, "certificates", "holder_projections_summary.json")))
    for d, r in summ.items():
        c = _cert(r["file"])
        assert c["problem"]["delta"] == d
        assert (c["claim"]["lower"], c["claim"]["upper"]) == (r["lower"], r["upper"])
