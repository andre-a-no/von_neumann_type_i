#!/usr/bin/env python3
r"""
Independent verifier of torch_vn_algebra envelope certificates.

    python verify.py CERT.json[.gz] [CERT2 ...]        (or: python -m torch_vn_algebra.certify.verify ...)

Depends on python-flint and the standard library only, and imports nothing from torch_vn_algebra: copy this file
next to a certificate and run it. It re-derives everything it checks:

1. Problem: the id and parameters must be known to this file (VERIFIERS below), and delta must be a rational in (0, 1).
2. Lower bound L: the witness matrices are checked in exact rational arithmetic in the ORIGINAL problem (P, Q
   symmetric idempotents; W positive diagonal with Tr W = 1 and Michelson contrast exactly delta), and the violation
   |phi(PQ)|^alpha >= (1 + L)^alpha phi(P) phi(Q)^(alpha - 1) is checked exactly. Hence V(delta) >= L.
3. Upper bound U: the reduced domain is rebuilt from delta; the search tree is replayed from the root box with exact
   rational midpoints (coverage by construction); on every leaf, one of three Arb bounds (natural extension,
   mean-value form, second-order Taylor form) must be certainly < U; every monotonicity step must have a certain
   sign of the partial derivative on its box. Hence sup F < U on the reduced domain.
4. The reductions from the original to the reduced problem are listed in the certificate with references. They are
   mathematical lemmas and are NOT checked here: the certificate is a computer proof of the reduced statement.

Exit status 0 if every certificate passes.
"""
import gzip
import json
import sys

from flint import arb, fmpq, ctx

FORMAT = "torch_vn_algebra-certificate/1"


# ---------------------------------------------------------------------------------------------------------- numbers
def Q(s):
    p, _, q = str(s).partition("/")
    return fmpq(int(p), int(q) if q else 1)


def B(x):
    """Ball containing the exact rational x."""
    x = fmpq(x)
    return arb(x.p) / arb(x.q)


# ----------------------------------------------------------------------- second-order forward AD on Arb balls
class T:
    """(value, gradient, Hessian) enclosures; written independently of torch_vn_algebra.certify.ad."""
    __slots__ = ("v", "g", "h")

    def __init__(s, v, g, h):
        s.v, s.g, s.h = v, g, h

    @staticmethod
    def var(ball, i, n):
        return T(ball, [arb(int(k == i)) for k in range(n)], [[arb(0)] * n for _ in range(n)])

    def lift(s, o):
        if isinstance(o, T):
            return o
        n = len(s.g)
        return T(B(o) if not isinstance(o, arb) else o, [arb(0)] * n, [[arb(0)] * n for _ in range(n)])

    def __add__(s, o):
        o = s.lift(o)
        return T(s.v + o.v, [a + b for a, b in zip(s.g, o.g)], [[a + b for a, b in zip(x, y)] for x, y in zip(s.h, o.h)])
    __radd__ = __add__

    def __neg__(s):
        return T(-s.v, [-a for a in s.g], [[-a for a in x] for x in s.h])

    def __sub__(s, o):
        return s + (-s.lift(o))

    def __rsub__(s, o):
        return s.lift(o) + (-s)

    def __mul__(s, o):
        o = s.lift(o)
        n = len(s.g)
        return T(s.v * o.v, [s.g[i] * o.v + s.v * o.g[i] for i in range(n)],
                 [[s.h[i][j] * o.v + s.g[i] * o.g[j] + s.g[j] * o.g[i] + s.v * o.h[i][j] for j in range(n)]
                  for i in range(n)])
    __rmul__ = __mul__

    def chain(s, f0, f1, f2):
        n = len(s.g)
        return T(f0, [f1 * a for a in s.g], [[f2 * s.g[i] * s.g[j] + f1 * s.h[i][j] for j in range(n)] for i in range(n)])

    def power(s, r):
        """s ** r, r rational; the argument must be positive (otherwise Arb gives an indeterminate ball)."""
        rb = B(r)
        return s.chain(s.v ** r, rb * s.v ** (r - 1), rb * (rb - 1) * s.v ** (r - 2))

    def inv(s):
        w = 1 / s.v
        return s.chain(w, -w * w, 2 * w * w * w)

    def __truediv__(s, o):
        return s * s.lift(o).inv()

    def __rtruediv__(s, o):
        return s.lift(o) * s.inv()


# --------------------------------------------------------------------------------------------- problem: item 14
class Holder:
    """|phi(PQ)| <= phi(P)^(1/alpha) phi(Q)^(1-1/alpha) on projections, phi = Tr(W .), contrast(W) = delta."""
    nvar = 3

    def __init__(self, params, delta):
        self.alpha = int(params["alpha"])
        if self.alpha != params["alpha"] or self.alpha < 2:
            raise ValueError("alpha must be an integer >= 2")
        self.delta = delta

    def domain(self):
        d = self.delta
        return [(fmpq(0), fmpq(1)), (fmpq(0), fmpq(1)), ((1 - d) / (1 + d), fmpq(1))]

    def F(self, x):
        """Reduced objective: u = (cos s, sin s), v = (cos t, sin t), W = diag(1, b); cos, sin via m = tan(s/2)."""
        m, mm, b = x
        cu, su = (1 - m * m) / (1 + m * m), 2 * m / (1 + m * m)
        cv, sv = (1 - mm * mm) / (1 + mm * mm), 2 * mm / (1 + mm * mm)
        uv = cu * cv + su * sv                     # <u, v>
        vwu = cu * cv + b * (su * sv)              # <v, W u>
        pu = cu * cu + b * (su * su)               # <u, W u>
        pv = cv * cv + b * (sv * sv)               # <v, W v>
        a = fmpq(1, self.alpha)
        den = pu.power(a) * pv.power(1 - a) if isinstance(pu, T) else pu ** a * pv ** (1 - a)
        return uv * vwu / den - 1

    def check_witness(self, w, L):
        """Exact check in M_n; returns a description. Raises on failure."""
        def mat(rows):
            return [[Q(v) for v in r] for r in rows]
        W, P, Qm = mat(w["W"]), mat(w["P"]), mat(w["Q"])
        n = len(W)
        if not all(len(M) == n and all(len(r) == n for r in M) for M in (W, P, Qm)):
            raise AssertionError("witness: shapes")

        def mul(X, Y):
            return [[sum((X[i][k] * Y[k][j] for k in range(n)), fmpq(0)) for j in range(n)] for i in range(n)]

        def tr(X):
            return sum((X[i][i] for i in range(n)), fmpq(0))

        for name, M in (("P", P), ("Q", Qm)):
            if mul(M, M) != M or any(M[i][j] != M[j][i] for i in range(n) for j in range(n)):
                raise AssertionError(f"witness: {name} is not an orthogonal projection")
        if any(W[i][j] != 0 for i in range(n) for j in range(n) if i != j):
            raise AssertionError("witness: W must be diagonal")
        ev = [W[i][i] for i in range(n)]
        if min(ev) <= 0 or tr(W) != 1:
            raise AssertionError("witness: W must be positive definite with trace 1")
        if (max(ev) - min(ev)) / (max(ev) + min(ev)) != self.delta:
            raise AssertionError("witness: contrast of W differs from delta")
        x, p, q = tr(mul(W, mul(P, Qm))), tr(mul(W, P)), tr(mul(W, Qm))
        A = self.alpha
        if not ((1 + L) ** A * p * q ** (A - 1) <= x ** A):              # |x| = x up to sign: x^alpha, alpha odd
            if not ((1 + L) ** A * p * q ** (A - 1) <= abs(x) ** A):
                raise AssertionError("witness: the violation is smaller than L")
        return f"exact: |phi(PQ)|^{A} >= (1 + L)^{A} phi(P) phi(Q)^{A - 1} in M_{n}"


VERIFIERS = {"holder_projections": Holder}


# ------------------------------------------------------------------------------------------------- box bounds
def balls(box):
    return [B(lo) if lo == hi else B(lo).union(B(hi)) for lo, hi in box]


def leaf_ok(F, box, Ub):
    n = len(box)
    Fb = F([T.var(b, i, n) for i, b in enumerate(balls(box))])
    if Fb.v < Ub:
        return True
    Fc = F([T.var(B((lo + hi) / 2), i, n) for i, (lo, hi) in enumerate(box)])
    r = [B((hi - lo) / 2) for lo, hi in box]
    act = [i for i in range(n) if box[i][0] != box[i][1]]
    mv = Fc.v + sum((abs(Fb.g[i]) * r[i] for i in act), arb(0))
    if mv < Ub:
        return True
    t2 = Fc.v + sum((abs(Fc.g[i]) * r[i] for i in act), arb(0)) \
        + sum((abs(Fb.h[i][j]) * r[i] * r[j] for i in act for j in act), arb(0)) / 2
    return bool(t2 < Ub)


def grad_sign(F, box, k):
    n = len(box)
    g = F([T.var(b, i, n) for i, b in enumerate(balls(box))]).g[k]
    return "-" if g < 0 else "+" if g > 0 else None


def replay(F, root, tokens, U):
    """Replay the preorder tree; every leaf must be certainly < U. Returns the number of leaves."""
    Ub = B(U)
    stack, it, leaves = [list(root)], iter(tokens), 0
    for tok in it:
        if not stack:
            raise AssertionError("tree: more tokens than boxes")
        box = stack.pop()
        if tok == "L":
            if not leaf_ok(F, box, Ub):
                raise AssertionError(f"tree: leaf {[(str(a), str(b)) for a, b in box]} is not certainly below U")
            leaves += 1
        elif tok[0] == "M":
            k, sign = int(tok[1:-1]), tok[-1]
            lo, hi = box[k]
            if lo == hi or grad_sign(F, box, k) != sign:
                raise AssertionError(f"tree: monotonicity step {tok} not justified on {box}")
            face = list(box)
            face[k] = (lo, lo) if sign == "-" else (hi, hi)
            stack.append(face)
        elif tok[0] == "S":
            k = int(tok[1:])
            lo, hi = box[k]
            if lo == hi:
                raise AssertionError("tree: split of a degenerate side")
            mid = (lo + hi) / 2
            left, right = list(box), list(box)
            left[k], right[k] = (lo, mid), (mid, hi)
            stack.append(right)
            stack.append(left)
        else:
            raise AssertionError(f"tree: unknown token {tok!r}")
    if stack:
        raise AssertionError("tree: unfinished (boxes left without a decision)")
    return leaves


# ------------------------------------------------------------------------------------------------------- main
def verify(path, log=print):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as f:
        c = json.loads(f.read().decode())
    if c.get("format") != FORMAT:
        raise AssertionError(f"unknown format {c.get('format')!r}")
    pr = c["problem"]
    delta = Q(pr["delta"])
    if not (0 < delta < 1):
        raise AssertionError("delta must lie in (0, 1)")
    if pr["id"] not in VERIFIERS:
        raise AssertionError(f"no verifier for problem {pr['id']!r}")
    prob = VERIFIERS[pr["id"]](pr["params"], delta)
    L, U = Q(c["claim"]["lower"]), Q(c["claim"]["upper"])
    if not L <= U:
        raise AssertionError("claim: L > U")
    how = prob.check_witness(c["witness"], L)
    up = c["upper"]
    ctx.prec = int(up["precision_bits"])
    root = prob.domain()
    if [[Q(a), Q(b)] for a, b in up["domain"]] != [[a, b] for a, b in root]:
        raise AssertionError("upper: domain differs from the one implied by delta")
    leaves = replay(prob.F, root, up["tree"].split(), U)
    if leaves != up["leaves"]:
        raise AssertionError("upper: leaf count differs")
    log(f"OK {path}: {pr['id']} {pr['params']} delta = {pr['delta']}: "
        f"{float(L):.12g} <= V <= {float(U):.12g}  [L {how}; U: {leaves} boxes at {ctx.prec} bits]")
    log("   reductions assumed (proved in the text): " + ", ".join(r["id"] for r in pr["reductions"]))
    return float(L), float(U)


def main(argv):
    if not argv:
        print(__doc__)
        return 2
    bad = 0
    for p in argv:
        try:
            verify(p)
        except Exception as e:      # report and continue with the next certificate
            bad += 1
            print(f"FAIL {p}: {type(e).__name__}: {e}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
