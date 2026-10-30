"""
Forward-mode automatic differentiation over Arb balls (python-flint), to first or second order.

A :class:`Jet` carries a ball ``val``, the balls ``grad[i]`` of the partial derivatives and, at order 2, the balls
``hess[i][j]`` of the second derivatives. Evaluated at balls that contain a box, each component encloses the range of
the corresponding derivative on the box. This is what the box bounds of :mod:`.engine` need:

    mean-value form       f(box) subset f(c) + sum_i grad_i(box) [x_i - c_i]
    second-order Taylor   f(box) subset f(c) + sum_i grad_i(c) [x_i - c_i] + 1/2 sum_ij hess_ij(box) [x_i - c_i][x_j - c_j]

Every operation is rigorous because Arb arithmetic is: the result ball contains every value of the exact operation
on points of the argument balls. Only + - * / and rational powers are implemented (enough for the problems of
:mod:`.problems`); a power x ** r with non-integer r needs x > 0 on the ball, otherwise Arb returns an indeterminate
ball and the box is simply not closed.
"""
from flint import arb, fmpq


def _ball(x):
    if isinstance(x, arb):
        return x
    if isinstance(x, fmpq):
        return arb(x.p) / arb(x.q)
    if isinstance(x, int):
        return arb(x)
    raise TypeError(f"not an exact scalar: {type(x).__name__} (floats are not allowed in certified code)")


class Jet:
    """f, its gradient and (order 2) its Hessian, as Arb balls, with respect to ``n`` variables."""

    __slots__ = ("val", "grad", "hess")

    def __init__(self, val, grad, hess=None):
        self.val, self.grad, self.hess = val, grad, hess

    @classmethod
    def variable(cls, ball, index, n, order=1):
        grad = [arb(1) if i == index else arb(0) for i in range(n)]
        hess = [[arb(0)] * n for _ in range(n)] if order >= 2 else None
        return cls(_ball(ball), grad, hess)

    def _const(self, x):
        n = len(self.grad)
        hess = [[arb(0)] * n for _ in range(n)] if self.hess is not None else None
        return Jet(_ball(x), [arb(0)] * n, hess)

    def _lift(self, other):
        return other if isinstance(other, Jet) else self._const(other)

    def _unary(self, f0, f1, f2):
        """f(self) for a scalar function with value f0, f' = f1, f'' = f2 at self.val (chain rule)."""
        g = [f1 * a for a in self.grad]
        h = None
        if self.hess is not None:
            n = len(g)
            h = [[f2 * self.grad[i] * self.grad[j] + f1 * self.hess[i][j] for j in range(n)] for i in range(n)]
        return Jet(f0, g, h)

    def __add__(self, other):
        o = self._lift(other)
        h = None if self.hess is None else [[a + b for a, b in zip(r, s)] for r, s in zip(self.hess, o.hess)]
        return Jet(self.val + o.val, [a + b for a, b in zip(self.grad, o.grad)], h)

    __radd__ = __add__

    def __neg__(self):
        h = None if self.hess is None else [[-a for a in r] for r in self.hess]
        return Jet(-self.val, [-a for a in self.grad], h)

    def __sub__(self, other):
        return self + (-self._lift(other))

    def __rsub__(self, other):
        return self._lift(other) - self

    def __mul__(self, other):
        if not isinstance(other, Jet):
            c = _ball(other)
            h = None if self.hess is None else [[c * a for a in r] for r in self.hess]
            return Jet(c * self.val, [c * a for a in self.grad], h)
        a, b = self, other
        g = [x * b.val + a.val * y for x, y in zip(a.grad, b.grad)]
        h = None
        if a.hess is not None:
            n = len(g)
            h = [[a.hess[i][j] * b.val + a.grad[i] * b.grad[j] + a.grad[j] * b.grad[i] + a.val * b.hess[i][j]
                  for j in range(n)] for i in range(n)]
        return Jet(a.val * b.val, g, h)

    __rmul__ = __mul__

    def reciprocal(self):
        v = self.val
        r = 1 / v
        return self._unary(r, -r * r, 2 * r * r * r)

    def __truediv__(self, other):
        if not isinstance(other, Jet):
            return self * (1 / _ball(other))
        return self * other.reciprocal()

    def __rtruediv__(self, other):
        return self._lift(other) * self.reciprocal()

    def __pow__(self, r):
        """x ** r for an exact rational (or integer) exponent r; x > 0 is required unless r is a natural number."""
        if isinstance(r, int) and r >= 0:
            out = self._const(1)
            for _ in range(r):
                out = out * self
            return out
        r = fmpq(r)
        rb, v = _ball(r), self.val
        return self._unary(v ** r, rb * v ** (r - 1), rb * (rb - 1) * v ** (r - 2))


Dual = Jet          # first-order use


def value(x):
    """The ball of a Jet or of a plain ball."""
    return x.val if isinstance(x, Jet) else _ball(x)
