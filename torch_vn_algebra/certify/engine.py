"""
Branch and bound in Arb: a rigorous upper bound U >= sup F on a box with rational endpoints.

Each box gets the smallest of three upper bounds, all rigorous (c the centre, r_i the half-widths; the derivative
enclosures come from forward-mode AD on the balls that contain the box, :mod:`.ad`):

* natural interval extension: F evaluated on the box;
* mean-value form: F(c) + sum_i |dF/dx_i (box)| r_i;
* second-order Taylor form: F(c) + sum_i |dF/dx_i (c)| r_i + 1/2 sum_ij |d2F/dx_i dx_j (box)| r_i r_j. This one is
  needed where F has a ridge just below U (item 14 has one along the diagonal m = m', where F = 0 identically): the
  mean-value form loses a factor |grad F| r there, the Taylor form only (r^2).

A box is closed when one of them is certainly below the target U. Otherwise, if some partial derivative has a fixed
sign on the box, the box is replaced by the face where F is larger (monotonicity test: sup over the box equals sup
over that face); else it is bisected along its relatively widest side. The search tree is recorded in preorder:

    'L'    leaf: bound < U on this box
    'Mk-'  dF/dx_k < 0 on the box: continue on the face x_k = lo        (one child)
    'Mk+'  dF/dx_k > 0 on the box: continue on the face x_k = hi        (one child)
    'Sk'   bisect at the exact midpoint of side k                       (two children: lower half, upper half)

The side is chosen by the smear rule (largest |dF/dx_k| times width); the verifier does not need to know the rule.

so the verifier can rebuild every box exactly from the root box, and coverage holds by construction.
"""
from flint import arb, fmpq, ctx

from .ad import Jet


def _ball(q):
    return arb(q.p) / arb(q.q)


def box_balls(box):
    return [_ball(lo) if lo == hi else _ball(lo).union(_ball(hi)) for lo, hi in box]


def bounds(F, box, order=2):
    """([natural, mean_value(, taylor2)], gradient on the box): balls whose upper ends bound F on the box."""
    n = len(box)
    Fb = F([Jet.variable(b, i, n, order) for i, b in enumerate(box_balls(box))])
    Fc = F([Jet.variable(_ball((lo + hi) / 2), i, n, 1) for i, (lo, hi) in enumerate(box)])
    r = [_ball((hi - lo) / 2) if lo != hi else None for lo, hi in box]
    act = [i for i in range(n) if r[i] is not None]
    mv = Fc.val
    for i in act:
        mv = mv + abs(Fb.grad[i]) * r[i]
    out = [Fb.val, mv]
    if order >= 2:
        t2 = Fc.val
        for i in act:
            t2 = t2 + abs(Fc.grad[i]) * r[i]
        quad = arb(0)
        for i in act:
            for j in act:
                quad = quad + abs(Fb.hess[i][j]) * r[i] * r[j]
        out.append(t2 + quad / 2)
    return out, Fb.grad


class BranchAndBound:
    def __init__(self, F, root, U, max_leaves=2_000_000, monotone=True, order=2, split="smear"):
        self.F, self.root, self.U, self.Ub = F, list(root), fmpq(U), _ball(fmpq(U))
        self.max_leaves, self.monotone, self.order, self.split = max_leaves, monotone, order, split
        self.tokens, self.leaves, self.max_depth = [], 0, 0
        self.widths = [hi - lo for lo, hi in root]

    def split_side(self, box, grad):
        """The side to bisect: largest |dF/dx_k (box)| * width_k (smear rule), relatively widest side as fallback."""
        act = [i for i in range(len(box)) if box[i][0] != box[i][1]]
        if not act:
            return None
        if self.split == "smear":
            sm = [(abs(grad[i]).upper() * _ball(box[i][1] - box[i][0]), i) for i in act]
            if all(v.is_finite() for v, _ in sm):
                return max(sm, key=lambda t: float(t[0].mid()))[1]
        return max(act, key=lambda i: (box[i][1] - box[i][0]) / self.widths[i])

    def closed(self, box):
        bs, grad = bounds(self.F, box, self.order)
        return any(b < self.Ub for b in bs), grad

    def run(self):
        stack = [(self.root, 0)]
        while stack:
            box, depth = stack.pop()
            self.max_depth = max(self.max_depth, depth)
            ok, grad = self.closed(box)
            if ok:
                self.tokens.append("L")
                self.leaves += 1
                if self.leaves > self.max_leaves:
                    raise RuntimeError(f"branch and bound: more than {self.max_leaves} leaves; raise U or max_leaves")
                continue
            if self.monotone:
                k = next((i for i, g in enumerate(grad) if box[i][0] != box[i][1] and (g < 0 or g > 0)), None)
                if k is not None:
                    lo, hi = box[k]
                    sign = "-" if grad[k] < 0 else "+"
                    face = list(box)
                    face[k] = (lo, lo) if sign == "-" else (hi, hi)
                    self.tokens.append(f"M{k}{sign}")
                    stack.append((face, depth + 1))
                    continue
            k = self.split_side(box, grad)
            if k is None:
                raise RuntimeError(f"branch and bound: a point box {box} is not below U = {self.U}; U < sup F?")
            lo, hi = box[k]
            mid = (lo + hi) / 2
            left, right = list(box), list(box)
            left[k], right[k] = (lo, mid), (mid, hi)
            self.tokens.append(f"S{k}")
            stack.append((right, depth + 1))       # pushed first, popped second: preorder lower, upper
            stack.append((left, depth + 1))
        return self


def upper_bound(F, root, U, prec=128, **kw):
    """Run the branch and bound with target U at `prec` bits; returns the finished BranchAndBound."""
    old = ctx.prec
    ctx.prec = prec
    try:
        return BranchAndBound(F, root, U, **kw).run()
    finally:
        ctx.prec = old
