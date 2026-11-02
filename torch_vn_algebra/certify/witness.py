"""
Exact rational witnesses.

A numerically found extremal point is replaced by a nearby point with rational coordinates at which the
constraints hold exactly, so that the violation can be checked without rounding:

* orthogonal (real unitary) matrices by the Cayley transform U = (1 - K)(1 + K)^{-1} of a rational skew-symmetric K
  (1 + K is invertible since the spectrum of K is imaginary; U^T U = 1 exactly);
* positive matrices as U diag(q_i) U^T with rational q_i >= 0 (projections: q_i in {0, 1});
* rational numbers from floats by continued fractions with a bounded denominator.

The complex case (U = (1 - iK)(1 + iK)^{-1}, K Hermitian) needs Gaussian rationals; the problems implemented so far
reduce to real configurations, so only the real case is provided.
"""
from fractions import Fraction

from flint import fmpq, fmpq_mat


def rational(x, max_denominator=10 ** 12):
    """The best rational approximation of the float x with denominator <= max_denominator, as fmpq."""
    f = Fraction(x).limit_denominator(max_denominator)
    return fmpq(f.numerator, f.denominator)


def parse_rational(s):
    """fmpq from 'p/q' or 'p'."""
    if isinstance(s, fmpq):
        return s
    p, _, q = str(s).partition("/")
    return fmpq(int(p), int(q) if q else 1)


def fmt_rational(x):
    x = fmpq(x)
    return str(x.p) if x.q == 1 else f"{x.p}/{x.q}"


def cayley_orthogonal(K):
    """U = (1 - K)(1 + K)^{-1} for a rational skew-symmetric K (fmpq_mat or nested lists); U^T U = 1 exactly."""
    K = K if isinstance(K, fmpq_mat) else fmpq_mat([[fmpq(v) for v in row] for row in K])
    n = K.nrows()
    if K.transpose() != -K:
        raise ValueError("cayley_orthogonal: K must be skew-symmetric")
    one = fmpq_mat(n, n, [1 if i == j else 0 for i in range(n) for j in range(n)])
    return (one - K) * (one + K).inv()


def rational_psd(U, q):
    """U diag(q) U^T with rational q (exactly positive semidefinite when q >= 0)."""
    n = U.nrows()
    D = fmpq_mat(n, n, [fmpq(q[i]) if i == j else 0 for i in range(n) for j in range(n)])
    return U * D * U.transpose()


def rotation_from_half_angle(m):
    """The Cayley transform of [[0, m], [-m, 0]]: [[c, -s], [s, c]], c = (1 - m^2)/(1 + m^2), s = 2m/(1 + m^2)."""
    m = fmpq(m)
    return cayley_orthogonal([[0, m], [-m, 0]])


def mat_to_json(M):
    return [[fmt_rational(M[i, j]) for j in range(M.ncols())] for i in range(M.nrows())]
