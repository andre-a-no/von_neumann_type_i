"""
torch_vn_algebra.certify: exactness of the witnesses, soundness of the Arb box bounds, the independent verifier,
and rejection of tampered certificates.
"""
import gzip
import importlib.util
import json
import os
import random

import numpy as np
import pytest

flint = pytest.importorskip("flint")
from flint import fmpq, ctx  # noqa: E402

from torch_vn_algebra import certify as C  # noqa: E402
from torch_vn_algebra.certify import engine, witness, ad  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CERT = os.path.join(ROOT, "certificates", "holder_projections_alpha3_delta_1_10.json.gz")
VERIFY = os.path.join(ROOT, "torch_vn_algebra", "certify", "verify.py")


def _verifier():
    spec = importlib.util.spec_from_file_location("tvna_verify_standalone", VERIFY)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _cert():
    with gzip.open(CERT, "rb") as f:
        return json.loads(f.read().decode())


@pytest.fixture(autouse=True)
def _prec():
    old = ctx.prec
    ctx.prec = 128
    yield
    ctx.prec = old


# ------------------------------------------------------------------------------------------------- witnesses
def test_cayley_is_exactly_orthogonal():
    rng = random.Random(0)
    for n in (2, 3, 5):
        K = [[fmpq(0)] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                K[i][j] = fmpq(rng.randint(-9, 9), rng.randint(1, 9))
                K[j][i] = -K[i][j]
        U = witness.cayley_orthogonal(K)
        assert U.transpose() * U == flint.fmpq_mat(n, n, [int(i == j) for i in range(n) for j in range(n)])
        P = witness.rational_psd(U, [1] + [0] * (n - 1))
        assert P * P == P and P.transpose() == P


def test_cayley_rejects_non_skew():
    with pytest.raises(ValueError):
        witness.cayley_orthogonal([[1, 0], [0, 0]])


def test_half_angle_rotation():
    m = fmpq(1, 3)
    U = witness.rotation_from_half_angle(m)
    assert (U[0, 0], U[1, 0]) == (fmpq(4, 5), fmpq(3, 5))


# ------------------------------------------------------------------------------------------------ AD and bounds
def _F14(delta="1/10"):
    P = C.HolderProjections(3)
    d = witness.parse_rational(delta)
    return P, d, (lambda x: P.objective(x, d))


def test_jet_matches_finite_differences():
    P, d, F = _F14()
    x0 = [0.31, 0.37, 0.85]
    J = F([ad.Jet.variable(witness.rational(v), i, 3, order=2) for i, v in enumerate(x0)])
    h = 1e-5
    for i in range(3):
        e = np.eye(3)[i] * h
        fd = (P._float(*(x0 + e)) - P._float(*(x0 - e))) / (2 * h)
        assert abs(float(J.grad[i].mid()) - fd) < 1e-8
        for j in range(3):
            ej = np.eye(3)[j] * h
            fdd = (P._float(*(x0 + e + ej)) - P._float(*(x0 + e - ej)) - P._float(*(x0 - e + ej))
                   + P._float(*(x0 - e - ej))) / (4 * h * h)
            assert abs(float(J.hess[i][j].mid()) - fdd) < 1e-4


@pytest.mark.parametrize("seed", range(5))
def test_box_bounds_enclose_sampled_values(seed):
    """Soundness: every one of the three bounds is >= F at every sampled point of the box."""
    P, d, F = _F14()
    rng = np.random.default_rng(seed)
    for _ in range(4):
        w = 2.0 ** -rng.integers(4, 9)
        lo = [rng.uniform(0, 1 - w), rng.uniform(0, 1 - w), rng.uniform(float(P.a_of(d)), 1 - w)]
        box = [(witness.rational(v, 2 ** 20), witness.rational(v, 2 ** 20) + witness.rational(w, 2 ** 20))
               for v in lo]
        bs, _ = engine.bounds(F, box)
        pts = rng.uniform([float(a) for a, _ in box], [float(b) for _, b in box], size=(400, 3))
        top = max(P._float(*p) for p in pts)
        for bnd in bs:
            assert top <= float(bnd.upper().mid()) + 1e-15


def test_loose_target_tree_replays_in_the_verifier():
    """A tree from the engine (splits, monotonicity steps, leaves) is accepted by the independent replay."""
    P, d, F = _F14()
    bb = engine.BranchAndBound(F, P.domain(d), fmpq(1, 100)).run()
    assert 1 < bb.leaves < 2000 and any(t.startswith("S") for t in bb.tokens)
    v = _verifier()
    prob = v.Holder({"alpha": 3}, d)
    assert v.replay(prob.F, prob.domain(), bb.tokens, fmpq(1, 100)) == bb.leaves
    with pytest.raises(AssertionError):
        v.replay(prob.F, prob.domain(), bb.tokens[:-1], fmpq(1, 100))      # unfinished tree


def test_too_small_target_is_refused():
    P, d, F = _F14()
    with pytest.raises(RuntimeError):
        engine.BranchAndBound(F, P.domain(d), fmpq(1, 10 ** 4), max_leaves=200).run()


# --------------------------------------------------------------------------------------- witness of item 14
def test_item14_witness_is_exact_and_beats_the_search_of_the_preprint():
    P = C.HolderProjections(3)
    pt, fval = P.search(fmpq(1, 10), starts=10)
    data, L = P.witness(pt, fmpq(1, 10))
    assert abs(float(L) - fval) < 1e-12
    assert float(L) > 2.7864788517748744e-04                   # trace preprint, certified search value at 0.1
    v = _verifier()
    v.Holder({"alpha": 3}, fmpq(1, 10)).check_witness(data, L)
    with pytest.raises(AssertionError):
        v.Holder({"alpha": 3}, fmpq(1, 10)).check_witness(data, L + fmpq(1, 10 ** 9))


def test_item14_asymptotics_match_the_preprint():
    """V(delta) = (1 - 2/alpha)^2 delta^2 / 4 + O(delta^3) (trace preprint); the witness lower bound follows it."""
    P = C.HolderProjections(3)
    for d in ("1/100", "1/50"):
        pt, _ = P.search(witness.parse_rational(d), starts=8)
        _, L = P.witness(pt, witness.parse_rational(d))
        dd = float(witness.parse_rational(d))
        assert abs(float(L) / (dd * dd / 36) - 1) < 3 * dd


# -------------------------------------------------------------------------------------------------- verifier
def test_verifier_does_not_import_the_library():
    src = open(VERIFY).read()
    assert "torch_vn_algebra" not in "".join(line for line in src.splitlines()
                                              if line.startswith(("import", "from")))
    assert "import torch" not in src


def test_stored_certificate_verifies():
    L, U = _verifier().verify(CERT, log=lambda *a: None)
    assert 2.798785e-4 < L <= U < 2.798789e-4


def _tampered(tmp_path, edit):
    c = _cert()
    edit(c)
    p = tmp_path / "t.json"
    p.write_text(json.dumps(c))
    return p


@pytest.mark.parametrize("zeroed", [("P",), ("Q",), ("P", "Q")], ids=["P", "Q", "both"])
def test_zero_projection_cannot_certify_a_tight_lower_bound(tmp_path, zeroed):
    """A zero denominator makes the powered witness check vacuous, even with L = U."""
    def edit(c):
        for name in zeroed:
            c["witness"][name] = [["0" for _ in row] for row in c["witness"][name]]
        c["claim"]["lower"] = c["claim"]["upper"]

    with pytest.raises(AssertionError, match="phi\\(P\\) and phi\\(Q\\) must be positive"):
        _verifier().verify(_tampered(tmp_path, edit), log=lambda *a: None)


@pytest.mark.parametrize("name", sorted(f for f in os.listdir(os.path.join(ROOT, "certificates"))
                                       if f.endswith(".json.gz")))
def test_published_witness_is_admissible(name):
    """Check the exact lower-bound contract for every published certificate, without replaying its tree."""
    with gzip.open(os.path.join(ROOT, "certificates", name), "rt") as f:
        c = json.load(f)
    v = _verifier()
    prob = v.Holder(c["problem"]["params"], v.Q(c["problem"]["delta"]))
    prob.check_witness(c["witness"], v.Q(c["claim"]["lower"]))


@pytest.mark.parametrize("edit", [
    lambda c: c["witness"]["W"][0].__setitem__(0, "56/100"),                 # contrast of W changed (trace 1.01)
    lambda c: c["claim"].__setitem__("lower", "3/10000"),                     # L above the witness value
    lambda c: c["claim"].__setitem__("upper", c["claim"]["lower"] + "0"),     # U < L
    lambda c: c["upper"].__setitem__("domain", [["0", "1"], ["0", "1"], ["9/10", "1"]]),   # domain too small
    lambda c: c["upper"].__setitem__("tree", "L"),                            # root box is not below U
    lambda c: c["problem"].__setitem__("id", "unknown"),
    lambda c: c["problem"].__setitem__("delta", "1/9"),
], ids=["W", "L", "U<L", "domain", "tree", "id", "delta"])
def test_tampered_certificates_fail(tmp_path, edit):
    with pytest.raises(AssertionError):
        _verifier().verify(_tampered(tmp_path, edit), log=lambda *a: None)


def test_small_certificate_round_trip(tmp_path):
    """certify -> save -> verify for a large contrast (fast), through the public API."""
    c = C.certify("holder_projections", "9/10", rel_tol=1e-3, starts=8)
    assert c.lower <= c.upper
    p = c.save(str(tmp_path / "c.json"))
    L, U = _verifier().verify(p, log=lambda *a: None)
    assert (L, U) == (float(c.lower), float(c.upper))
    assert C.load(p).data == c.data


# --------------------------------------------------------------------- public contrast and relative-target contract
@pytest.mark.parametrize("contrast, expected", [
    (1e-8, fmpq(*float(1e-8).as_integer_ratio())),
    (0.99999999, fmpq(*float(0.99999999).as_integer_ratio())),
    (0.123456789, fmpq(*float(0.123456789).as_integer_ratio())),
    (0.1, fmpq(*float(0.1).as_integer_ratio())),
    (fmpq(1, 10), fmpq(1, 10)),
], ids=["small-float", "near-one-float", "general-float", "decimal-float", "native-fmpq"])
def test_requested_contrast_reaches_search_without_approximation(monkeypatch, contrast, expected):
    """A bounded-denominator replacement changes the optimization problem, even before its search starts."""
    class SearchReached(Exception):
        pass

    def search(self, delta, starts):
        assert delta == expected
        raise SearchReached

    monkeypatch.setattr(C.HolderProjections, "search", search)
    with pytest.raises(SearchReached):
        C.certify("holder_projections", contrast)


def test_zero_envelope_alpha_two_stops_before_search(monkeypatch):
    """Cauchy--Schwarz gives V=0; a relative positive-width search cannot resolve that case."""
    def search(*args, **kwargs):
        pytest.fail("the exactly zero envelope must not enter numerical search")

    monkeypatch.setattr(C.HolderProjections, "search", search)
    with pytest.raises(ValueError, match="alpha=2 has zero envelope"):
        C.certify("holder_projections", "1/10", alpha=2, max_leaves=1)


@pytest.mark.parametrize("lower", [fmpq(0), -fmpq(1, 10 ** 20)], ids=["zero", "negative"])
def test_nonpositive_witness_bound_cannot_start_relative_cover(monkeypatch, lower):
    """An unresolved search must fail explicitly, not send U<=0 to an unbounded leftmost proof branch."""
    monkeypatch.setattr(C.HolderProjections, "search", lambda self, delta, starts: ([0, 0, 1], 0.0))
    monkeypatch.setattr(C.HolderProjections, "witness", lambda self, point, delta: (_cert()["witness"], lower))

    def cover(*args, **kwargs):
        pytest.fail("a nonpositive lower bound must not enter branch and bound")

    monkeypatch.setattr(C, "upper_bound", cover)
    ctx.prec = 96
    with pytest.raises(ValueError, match="strictly positive witness lower bound"):
        C.certify("holder_projections", "1/10", prec=128)
    assert ctx.prec == 96


@pytest.mark.parametrize("tolerance", [0, -1e-6, float("nan"), float("inf")])
def test_invalid_relative_tolerance_stops_before_search(monkeypatch, tolerance):
    def search(*args, **kwargs):
        pytest.fail("invalid relative tolerance must be rejected before search")

    monkeypatch.setattr(C.HolderProjections, "search", search)
    with pytest.raises(ValueError, match="rel_tol must be a finite positive number"):
        C.certify("holder_projections", "1/10", rel_tol=tolerance)


def test_float_contrast_exact_relative_width_and_native_delta_round_trip(tmp_path, monkeypatch):
    """The complete producer/replay path uses the supplied binary contrast and an exact rational target."""
    from fractions import Fraction

    delta = 0.9
    c = C.certify("holder_projections", delta, rel_tol=Fraction(1, 1000), starts=8)
    assert c.delta == fmpq(*delta.as_integer_ratio())
    assert c.upper / c.lower == fmpq(1001, 1000)
    p = c.save(tmp_path / "binary-contrast.json.gz")
    assert _verifier().verify(p, log=lambda *args: None) == (float(c.lower), float(c.upper))

    class SearchReached(Exception):
        pass

    def search(self, requested, starts):
        assert requested == c.delta
        raise SearchReached

    monkeypatch.setattr(C.HolderProjections, "search", search)
    with pytest.raises(SearchReached):
        C.certify("holder_projections", C.load(p).delta)


def test_relative_target_does_not_underflow_through_binary_float():
    """A constant positive objective has a known exact supremum even below the binary64 range."""
    tiny = fmpq(1, 10 ** 400)

    class ConstantEnvelope(C.Problem):
        id = "constant-envelope-test"
        statement = "supremum of the constant 10^-400 on [0,1]"
        variables = ("x",)
        reductions = []

        def domain(self, delta):
            return [(fmpq(0), fmpq(1))]

        def objective(self, x, delta):
            return x[0] * 0 + tiny

        def search(self, delta, starts=40, seed=0):
            return [0], 0.0

        def witness(self, point, delta):
            return {"point": "0"}, tiny

    c = C.certify(ConstantEnvelope(), "1/2", rel_tol=fmpq(1, 1000))
    assert c.lower == tiny
    assert c.upper == fmpq(1001, 1000 * 10 ** 400)
    assert c.data["upper"]["leaves"] == 1


def test_positive_but_impossible_upper_target_obeys_cover_budget(monkeypatch):
    """An optimizer can miss a peak: an unclosed leftmost path must still consume the leaf budget."""
    class MissedPeak(C.Problem):
        id = "missed-peak-budget-test"
        statement = "supremum of 1-x^2 on [0,1] is 1"
        variables = ("x",)
        reductions = []

        def domain(self, delta):
            return [(fmpq(0), fmpq(1))]

        def objective(self, x, delta):
            return 1 - x[0] ** 2

        def search(self, delta, starts=40, seed=0):
            return [fmpq(3, 4)], 7 / 16

        def witness(self, point, delta):
            return {"point": "3/4"}, fmpq(7, 16)

    original = engine.BranchAndBound.closed
    calls = 0

    def bounded(self, box):
        nonlocal calls
        calls += 1
        if calls > 32:
            pytest.fail("the leaf budget did not limit an unclosed proof branch")
        return original(self, box)

    monkeypatch.setattr(engine.BranchAndBound, "closed", bounded)
    with pytest.raises(RuntimeError, match="more than 8 leaves"):
        C.certify(MissedPeak(), "1/2", rel_tol=fmpq(1, 100), max_leaves=8)
    assert calls <= 8


@pytest.mark.parametrize("budget", [0, -1, 1.5])
def test_cover_requires_positive_integer_leaf_budget(budget):
    with pytest.raises(ValueError, match="max_leaves must be a positive integer"):
        engine.BranchAndBound(lambda x: x[0], [(fmpq(0), fmpq(1))], fmpq(2), max_leaves=budget)


def test_single_proved_root_fits_one_leaf_budget():
    bb = engine.BranchAndBound(lambda x: x[0], [(fmpq(0), fmpq(1))], fmpq(2), max_leaves=1).run()
    assert bb.leaves == bb.potential_leaves == 1
    assert bb.tokens == ["L"]
