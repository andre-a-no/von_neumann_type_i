r"""
Envelope problems V(delta) = sup F that :func:`certify` can enclose.

A problem fixes, for a rational contrast delta,

* ``domain(delta)``: the box of reduced variables (rational endpoints) on which the branch and bound runs, after the
  reductions listed in ``reductions`` (each with a proof reference; the verifier does not check them, they are the
  human-readable part of the certificate);
* ``objective(x, delta)``: F on that box, written with + - * / and rational powers only, so that it can be evaluated
  on Arb balls and on :class:`.ad.Dual` numbers;
* ``search(delta)``: a float maximiser (any heuristic; it only has to be good);
* ``witness(point, delta)``: an exact rational configuration of the ORIGINAL problem near the point, and the exact
  rational lower bound it proves.

Item 14 of the trace catalogue (Hölder's inequality for projections)
--------------------------------------------------------------------
For a positive functional phi = Tr(W .) on M_n with Michelson contrast Delta(W) = delta and projections P, Q,

    F(P, Q) = |phi(PQ)| / (phi(P)^{1/alpha} phi(Q)^{1 - 1/alpha}) - 1,      V(delta) = sup F       (alpha = 3).

Reduction (valid for every n and all ranks of P, Q).

R1 (Halmos). C^n is an orthogonal sum of subspaces H_i of dimension 1 or 2 invariant under P and Q; on each
2-dimensional H_i, P and Q are rank one. Since P, Q are block diagonal, phi(PQ) = sum_i Tr(W_i P_i Q_i),
phi(P) = sum_i Tr(W_i P_i), phi(Q) = sum_i Tr(W_i Q_i), with W_i the compression of W to H_i.
R2 (compression). The spectrum of W_i lies in [lambda_min(W), lambda_max(W)], so Delta(W_i) <= delta.
R3 (Hölder). With x_i = Tr(W_i P_i Q_i), p_i = Tr(W_i P_i), q_i = Tr(W_i Q_i):
|sum x_i| <= sum |x_i| <= (1 + V_2) sum p_i^{1/alpha} q_i^{1 - 1/alpha} <= (1 + V_2) p^{1/alpha} q^{1 - 1/alpha},
where V_2 is the supremum over one block (2x2, contrast <= delta, rank-one P, Q). One-dimensional blocks give
x_i = p_i = q_i or x_i = 0, i.e. no violation. Hence V_n <= V_2, and V_n >= V_2 by embedding (n >= 2).
R4 (unitary invariance). On a block, diagonalise W_i = diag(1, b) (up to scale, F is homogeneous of degree 0 in W),
b = lambda_min/lambda_max in [a, 1], a = (1 - delta)/(1 + delta).
R5 (phases). P = uu^*, Q = vv^*, u = (cos s, e^{i alpha} sin s), v = (cos t, e^{i beta} sin t), s, t in [0, pi/2].
diag(1, e^{-i alpha}) commutes with W and removes alpha. Then |u^*v| = |cos s cos t + e^{i beta} sin s sin t| and
|v^*Wu| = |cos s cos t + b e^{-i beta} sin s sin t| are both maximal at beta = 0 (triangle inequality, all terms
>= 0), while phi(P), phi(Q) do not depend on beta.
R6 (rational parametrisation). cos s = (1 - m^2)/(1 + m^2), sin s = 2m/(1 + m^2) with m = tan(s/2) in [0, 1]; the
same for t with m'. The reduced problem is

    F(m, m', b) = (c c' + s s')(c c' + b s s') / ((c^2 + b s^2)^{1/alpha} (c'^2 + b s'^2)^{1 - 1/alpha}) - 1

on [0, 1] x [0, 1] x [a, 1]: three variables, rational except for the two real powers of quantities >= b >= a > 0.
"""
import math

import numpy as np
from flint import fmpq, arb

from .witness import rational, rotation_from_half_angle, rational_psd, mat_to_json, fmt_rational


class Problem:
    id = None
    variables = ()

    def __init__(self, **params):
        self.params = params

    def domain(self, delta):
        raise NotImplementedError

    def objective(self, x, delta):
        raise NotImplementedError

    def search(self, delta, starts=40, seed=0):
        raise NotImplementedError

    def witness(self, point, delta):
        raise NotImplementedError

    def describe(self):
        return dict(id=self.id, params={k: fmt_rational(v) if isinstance(v, fmpq) else v
                                        for k, v in self.params.items()},
                    statement=self.statement, variables=list(self.variables), reductions=self.reductions)


def _cs(m):
    """(cos, sin) of 2 arctan(m), as rational functions of m (works on fmpq, arb and Dual)."""
    d = 1 + m * m
    return (1 - m * m) / d, 2 * m / d


class HolderProjections(Problem):
    """Item 14: |phi(PQ)| <= phi(P)^{1/alpha} phi(Q)^{1-1/alpha} for projections, alpha an integer >= 3."""
    id = "holder_projections"
    variables = ("m", "m'", "b")
    statement = ("V(delta) = sup |phi(PQ)| / (phi(P)^(1/alpha) phi(Q)^(1-1/alpha)) - 1 over projections P, Q in M_n "
                 "(every n >= 2, all ranks) and densities W of phi = Tr(W .) with Michelson contrast delta")
    reductions = [
        dict(id="halmos", text="two-projection theorem: orthogonal sum of invariant blocks of dimension 1 or 2, "
                               "rank-one P, Q on 2-dimensional blocks", ref="Halmos, Trans. AMS 144 (1969); module "
                                                                           "docstring R1"),
        dict(id="compression", text="compressions of W to the blocks have contrast <= delta", ref="R2"),
        dict(id="holder", text="the block values combine by the triangle and Hölder inequalities, so V_n = V_2 "
                               "for every n >= 2", ref="R3"),
        dict(id="unitary_invariance", text="W = diag(1, b), b in [(1-delta)/(1+delta), 1]", ref="R4"),
        dict(id="phases", text="phases removed by a diagonal unitary and the triangle inequality: real u, v in "
                               "the first quadrant", ref="R5"),
        dict(id="half_angle", text="cos, sin by the half-angle tangent m in [0, 1]", ref="R6"),
    ]

    def __init__(self, alpha=3):
        if int(alpha) != alpha or alpha < 2:
            raise ValueError("alpha must be an integer >= 2 (the exact witness check raises to the power alpha)")
        super().__init__(alpha=int(alpha))
        self.alpha = int(alpha)

    @staticmethod
    def a_of(delta):
        delta = fmpq(delta)
        return (1 - delta) / (1 + delta)

    def domain(self, delta):
        return [(fmpq(0), fmpq(1)), (fmpq(0), fmpq(1)), (self.a_of(delta), fmpq(1))]

    def objective(self, x, delta):
        m1, m2, b = x
        c1, s1 = _cs(m1)
        c2, s2 = _cs(m2)
        cc, ss = c1 * c2, s1 * s2
        p = c1 * c1 + b * s1 * s1
        q = c2 * c2 + b * s2 * s2
        r = fmpq(1, self.alpha)
        return (cc + ss) * (cc + b * ss) / (p ** r * q ** (1 - r)) - 1

    def _float(self, m1, m2, b):
        c1, s1 = math.cos(2 * math.atan(m1)), math.sin(2 * math.atan(m1))
        c2, s2 = math.cos(2 * math.atan(m2)), math.sin(2 * math.atan(m2))
        p, q = c1 * c1 + b * s1 * s1, c2 * c2 + b * s2 * s2
        r = 1 / self.alpha
        return (c1 * c2 + s1 * s2) * (c1 * c2 + b * s1 * s2) / (p ** r * q ** (1 - r)) - 1

    def search(self, delta, starts=40, seed=0):
        from scipy.optimize import minimize
        b = float(self.a_of(delta))
        rng = np.random.default_rng(seed)
        best = (-np.inf, None)
        for _ in range(starts):
            res = minimize(lambda z: -self._float(*np.clip(z, 0, 1), b), rng.uniform(0, 1, 2), method="L-BFGS-B",
                           bounds=[(0, 1), (0, 1)], options=dict(ftol=1e-16, gtol=1e-14))
            if -res.fun > best[0]:
                best = (-float(res.fun), np.clip(res.x, 0, 1))
        return [float(best[1][0]), float(best[1][1]), b], best[0]

    def witness(self, point, delta, max_denominator=10 ** 9):
        """Exact witness in M_2: P = uu^T, Q = vv^T with u, v columns of Cayley rotations, W = diag((1+d)/2, (1-d)/2).

        Returns (json-able dict, L) where L is a rational with (1 + L)^alpha phi(P) phi(Q)^(alpha-1) <= |phi(PQ)|^alpha
        checked exactly, i.e. V(delta) >= F(P, Q) >= L."""
        delta = fmpq(delta)
        m1, m2 = rational(point[0], max_denominator), rational(point[1], max_denominator)
        U1, U2 = rotation_from_half_angle(m1), rotation_from_half_angle(m2)
        P = rational_psd(U1, [1, 0])
        Q = rational_psd(U2, [1, 0])
        W = rational_psd(rotation_from_half_angle(0), [(1 + delta) / 2, (1 - delta) / 2])
        tr = lambda M: sum((M[i, i] for i in range(M.nrows())), fmpq(0))  # noqa: E731
        x, p, q = tr(W * P * Q), tr(W * P), tr(W * Q)
        A = self.alpha
        ratio = x ** A / (p * q ** (A - 1))                            # (1 + F)^alpha, exact
        lam = (arb(ratio.p) / arb(ratio.q)) ** fmpq(1, A) - 1          # F, as a ball
        L = _round_down(lam)
        assert (1 + L) ** A * p * q ** (A - 1) <= x ** A, "witness: exact check failed"
        data = dict(n=2, W=mat_to_json(W), P=mat_to_json(P), Q=mat_to_json(Q),
                    construction=f"P = uu^T, Q = vv^T, u, v first columns of the Cayley rotations of m = "
                                 f"{fmt_rational(m1)}, {fmt_rational(m2)}; W = diag((1+delta)/2, (1-delta)/2)",
                    point=[fmt_rational(m1), fmt_rational(m2), fmt_rational(self.a_of(delta))])
        return data, L


def _round_down(ball, digits=15):
    """A rational <= every point of the ball, with about `digits` significant digits."""
    lo = ball.lower()
    f = float(lo.mid())                                   # nearest double to the lower end, then step down
    e = max(-300, int(math.floor(math.log10(abs(f)))) - digits) if f != 0 else -300
    k = math.floor(f / 10.0 ** e) - 2
    L = fmpq(k, 1) * fmpq(10) ** e if e >= 0 else fmpq(k, 10 ** (-e))
    if not (arb(L.p) / arb(L.q) < ball):
        raise ArithmeticError("round_down: could not separate the bound from the ball")
    return L


def _round_up(x, digits=15):
    """A rational >= x (x a float or ball), with about `digits` significant digits."""
    f = float(x.upper().mid()) if isinstance(x, arb) else float(x)
    e = int(math.floor(math.log10(abs(f)))) - digits if f != 0 else -300
    k = math.ceil(f / 10.0 ** e) + 2
    return fmpq(k, 1) * fmpq(10) ** e if e >= 0 else fmpq(k, 10 ** (-e))


PROBLEMS = {"holder_projections": HolderProjections}
