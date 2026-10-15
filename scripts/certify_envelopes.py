#!/usr/bin/env python3
r"""
Certified upper bounds for the upper envelopes of Exp. 1 and Exp. 3 (paper, Section 6), C = 1, k = 2, 3, 4 (and 6).

The multistart search of scripts/bounds_search.py only gives LOWER bounds for

    sup z,  Exp 1:  z = |Tr(X U Y)| - Tr(X Y)        X, Y >= 0,  U unitary
            Exp 3:  z = Tr(Y |X| Y) - Tr|Y X Y|      X = X^*,    Y >= 0

with the normalisation lambda_max = 1 and one contrast prescribed, the other free in [0, D] (D = 0.999 as in the
paper). This script proves upper bounds. Most of them are exact (Theorems 1 and 4 below; for k = 2 they give the
closed forms of Section 6 with independent proofs); the case left open by the theorems (Exp. 3, k even >= 4,
envelope along Delta(|X|)) is reduced to a model in 2m - 1 = k - 1 angles and certified by a rigorous branch and
bound.
Notation: a = (1 - delta)/(1 + delta) is lambda_min of the operator with prescribed contrast delta, b = (1 - D)/(1 + D)
the smallest admissible lambda_min of the free one, n the dimension. All statements hold in the real and in the
complex algebra (every reduction ends in real configurations), so the real envelopes of the paper have the same
suprema.

    s1(a, b) = (1 - a)^2 (1 - b)^2 / (4 (1 + a)(1 + b)) = delta^2 D^2 / ((1 + delta)(1 + D))
    s3(a, b) = 2 a (1 - b)^2 / (1 + a)                  = (1 - delta) (2 D / (1 + D))^2,     s0(b) = s3(1, b) = (1 - b)^2

Lemma 0 (two projections). For projections P, Q on a finite-dimensional Hilbert space H (real or complex), H is an
orthogonal sum of subspaces of dimension 1 or 2 invariant under P and Q, and on each 2-dimensional one P and Q are
rank-one projections. Proof: PQP restricted to ran P is self-adjoint; take an orthonormal eigenbasis v_i,
PQP v_i = mu_i v_i, and put W_i = span{v_i, Q v_i}. P Q v_i = PQP v_i = mu_i v_i, so W_i is invariant under P and Q.
For i != j, <v_i, Q v_j> = <v_i, PQP v_j> = 0 and <Q v_i, Q v_j> = <v_i, PQP v_j> = 0, so the W_i are orthogonal.
Their sum contains ran P; its orthogonal complement is invariant under Q and lies in ker P: diagonalise Q there. If
dim W_i = 2, P|W_i = v_i v_i^* and Q|W_i (Q W_i = span{Q v_i}) have rank one. In a 2-dimensional block with basis
(f, f'), P = f f^* and Q = q q^*, q = cos(theta) f + e^{i psi} sin(theta) f'; the unitary diag(1, e^{-i psi}) commutes
with P and removes the phase, and the reflection f' -> -f' maps theta to -theta: WLOG theta in [0, pi/2].

Theorem 1 (Exp. 1, C sectors of size k, Delta(X) = delta, Delta(Y) in [0, D] free):
        sup z = C floor(k/2) s1(a, b) = C floor(k/2) delta^2 D^2 / ((1 + delta)(1 + D)),
attained by direct sums of the k = 2 solution. By symmetry the same along Delta(Y) (swap delta and D), and with both
contrasts prescribed sup z = C floor(k/2) s1(a_X, a_Y).
Proof. (i) |Tr(X U Y)| = |Tr(U Y X)| <= ||Y X||_1 = ||X Y||_1, with equality for the (block-diagonal; real if X, Y
are real) polar unitary of (Y X)^*; so sup z = sup F, F(X, Y) = ||X Y||_1 - Tr(X Y). (ii) Fix the eigenbasis of X;
x -> F is convex (a norm of a linear map minus a linear map). The admissible spectra {x in [a, 1]^n : some x_i = 1,
some x_j = a} form a union of faces of the cube, and a convex function on a box attains its maximum at a vertex: so
X = a + (1 - a) P with a projection P. The same for Y: y in [b, 1]^n with some y_i = 1 (Delta(Y) <= D), so
Y = b + (1 - b) Q (or Y = 1, where F = 0). (iii) X, Y and hence F split over the blocks of Lemma 0 (in every sector):
||.||_1 and Tr are additive over orthogonal invariant subspaces. On a 1-dimensional block X Y = x y >= 0
contributes 0. On a 2-dimensional block X = diag(1, a), Y = b + (1 - b) q q^*; with c = cos^2(theta),
alpha = 1 + a b, beta = a + b one has Tr(XY) = c alpha + (1 - c) beta and, for a 2x2 matrix M,
||M||_1^2 = ||M||_F^2 + 2 |det M| = Tr(X^2 Y^2) + 2 a b = c alpha^2 + (1 - c) beta^2. So the block contributes
g(c) = sqrt(c alpha^2 + (1 - c) beta^2) - (c alpha + (1 - c) beta), a concave function of c (square root of an affine
function minus an affine function). g'(c) = 0 gives sqrt(...) = (alpha + beta)/2, c* = (alpha + 3 beta) /
(4 (alpha + beta)) in [0, 1], g(c*) = (alpha - beta)^2 / (4 (alpha + beta)) = s1(a, b) (alpha - beta =
(1 - a)(1 - b), alpha + beta = (1 + a)(1 + b)). (iv) A sector of size k has at most floor(k/2) two-dimensional
blocks, and P, Q with floor(k/2) blocks at c = c* attain the bound. QED. (So the direct sums of the paper, check
(ii), are optimal for Exp. 1, for every k and C.)

Lemma 2 (Exp. 3). Let X = X^* be invertible, P the projection onto its positive part, S = 2P - 1, and Y > 0.
 (a) z = min_{-1 <= T <= 1} [Tr(|X| Y^2) - Tr(T Y X Y)] = 2 min_{0 <= E <= 1} Tr(E A + (1 - E) B),
     A = Y X_+ Y, B = Y X_- Y   (||H||_1 = max_{-1<=T<=1} Tr(T H) for H = H^*; E = (1 - T)/2).
 (b) z is non-decreasing in |X| among operators commuting with P: if |X| <= |X'|, both commute with P and
     X' = S |X'|, then z(X, Y) <= z(X', Y). Proof: z = f(A, B) with f(A, B) = Tr A + Tr B - ||A - B||_1, and for
     A' = A + D, D >= 0, ||A' - B||_1 <= ||A - B||_1 + Tr D, so f(A', B) >= f(A, B); A = Y |X|^{1/2} P |X|^{1/2} Y
     is non-decreasing in |X| (P commutes with |X|), likewise B.
 (c) For fixed X and a fixed eigenbasis of Y, z is convex in the eigenvalues of Y^2. Proof: N = |X|^{1/2} Y is
     invertible and E -> F = N E N^* maps {0 <= E <= 1} onto {0 <= F <= G}, G = N N^* = |X|^{1/2} Y^2 |X|^{1/2};
     Tr(E A) = Tr(P F), Tr((1 - E) B) = Tr(P^perp (G - F)). Hence z = 2 min_{0 <= F <= G} Tr(P F + P^perp (G - F)),
     the partial minimum of a linear function over the convex cone {(F, G): 0 <= F <= G}: convex in G, and G is
     affine in Y^2.
 (d) z(-X, Y) = z(X, Y); z = 0 if X >= 0 or if X, Y commute (paper, Proposition in Section 6).
Reduction. Let |X| have contrast delta (eigenvalues in [a, 1], both ends attained) and let e be a unit eigenvector
of |X| for the eigenvalue a, chosen in ran P or ran P^perp (|X| commutes with P) and in one sector; by (d) WLOG
e in ran P. Then |X| <= |X'| := 1 - (1 - a) e e^*, and X' = S - (1 - a) e e^* has the same contrast (n >= 2), so by
(b) z(X, Y) <= z(X', Y). By (c) and the vertex argument of Theorem 1 (ii) applied to y^2 in [b^2, 1]^n,
Y = b + (1 - b) Q (Y = 1 gives z = 0). Apply Lemma 0 to (P, Q) in every sector.

Lemma 3 (k = 2). For X = diag(a, -1), Y = b + (1 - b) q q^*, sup over theta of z = s3(a, b), attained at
c* = cos^2(theta*) = (1 + a b)/((1 + a)(1 + b)). Proof. With K = Y^2 = b^2 + (1 - b^2) q q^*, det(Y X Y) < 0, so
||Y X Y||_1 = sqrt(Tr(YXY)^2 - 4 det(YXY)) and z = T - sqrt(T^2 - w), T = a K_11 + K_22, w = 4 a K_12^2 =
4 a (1 - b^2)^2 c (1 - c). The root r = T - sqrt(T^2 - w) of r (2T - r) = w lies in [0, T], where r -> r (2T - r) is
increasing. With R = s3(a, b): if R >= T then r <= T <= R; otherwise R (2T - R) - w = 4 a (1 - b^2)^2 (c - c*)^2 >= 0
(a polynomial identity in c: the quadratic has a double root; checked symbolically in the tests), so r <= R, with
equality at c = c*. QED.

Theorem 4 (Exp. 3, C sectors of size k). (i) Along Delta(Y) = delta_Y, Delta(|X|) in [0, D] free:
sup z = C floor(k/2) s0(a_Y) = C floor(k/2) (2 delta_Y / (1 + delta_Y))^2. (ii) Along Delta(|X|) = delta in (0, 1),
Delta(Y) in [0, D] free: k odd: sup z = C floor(k/2) s0(b), independent of delta; k = 2: sup z = (C - 1) s0(b) +
s3(a, b) (for C = 1 the paper's closed form); k even:
    (C k/2 - 1) s0(b) + s3(a, b) <= sup z <= (C - 1)(k/2) s0(b) + max((k/2 - 1) s0(b), Z_{k/2}(a, b)),
where Z_m is the supremum of the m-block model below (the lower bound is the direct sum; the script certifies
Z_2 and Z_3 numerically, i.e. k = 4, 6).
Proof. By the reduction, z <= z(X', Y) with X' = S - (1 - a) e e^*, Y = b + (1 - b) Q; also z(X', Y) <= z(S, Y) by
(b) (|X'| <= 1), and z(S, Y) is the sum over the blocks of Lemma 0 of k = 2 values at delta = 0, each <= s0(b) by
Lemma 3, so a sector with m two-dimensional blocks contributes <= m s0(b), m <= floor(k/2). (i): along Delta(Y) the
free contrast of |X| allows |X| = 1, so (b) gives z <= z(S, Y) <= C floor(k/2) s0(a_Y) directly (now
Y = a_Y + (1 - a_Y) Q by (c)), attained by floor(k/2) blocks per sector at theta = pi/4 and |X| = 1. (ii) The
sectors without e contribute <= floor(k/2) s0(b) each. k odd: so does the sector of e, and the bound is attained with
blocks at theta = pi/4, |X| = 1 and the eigenvalue a on a remaining 1-dimensional block (which contributes 0). k even:
if the sector of e has m < k/2 blocks it contributes <= (k/2 - 1) s0(b); if m = k/2 its blocks fill it,
ran P = span{f_1, ..., f_m} and e = sum_j c_j f_j: this is the m-block model
    X = S - (1 - a) e e^*,  Y = b + (1 - b) sum_j q_j q_j^*,  S = diag(1, -1, 1, -1, ...),
    q_j = cos(theta_j) f_j + sin(theta_j) f_j',  e = sum_j c_j f_j,
with supremum Z_m (Z_1 = s3(a, b) by Lemma 3, which gives k = 2). Multiplying block j by a phase or a sign (this
commutes with S and Q) makes c_j >= 0, and Lemma 0 makes theta_j in [0, pi/2]: Z_m is the supremum over
theta in [0, pi/2]^m and e in the positive orthant, e = e(phi) in hyperspherical coordinates
e_l = sin(phi_1)...sin(phi_l) cos(phi_{l+1}), e_{m-1} = sin(phi_1)...sin(phi_{m-1}), phi in [0, pi/2]^{m-1}. For m = 2
the exchange of the blocks maps (theta_1, theta_2, phi) to (theta_2, theta_1, pi/2 - phi), so phi in [0, pi/4]
suffices. The direct sum is e = f_1, theta_1 = theta*, theta_j = pi/4 (j >= 2): value s3(a, b) + (m - 1) s0(b). QED.

Corollary 5 (all contrasts at once). sup z along Delta(|X|) is non-increasing in delta: for delta <= delta' (a >= a')
the reduction maps a configuration at delta' to X' = S - (1 - a') e e^*, and S - (1 - a) e e^* has contrast delta and
larger |X|, so by Lemma 2 (b) at least the same z. Hence for delta in [delta_i, delta_{i+1}],
dsum(delta) <= sup z(delta) <= UB(delta_i), and max_i UB(delta_i) - dsum(delta_{i+1}) bounds the gap to the direct
sums dsum (also non-increasing) on all of [0, delta_N] (certify_all_contrasts).

Certificate for Z_m (branch and bound). For a box B with centre c and half-widths r, let T_c be the sign of
Y X Y at c (Y X Y is invertible: an orthogonal symmetry; numerically renormalised to ||T_c|| <= 1, see below). By
Lemma 2 (a), z(p) <= Phi(p) := Tr(|X| Y^2) - Tr(T_c Y X Y) for every p in B, and Phi is smooth. Along p(t) = c + t h,
Taylor's formula with integral remainder gives Phi(c + h) <= Phi(c) + grad Phi(c) . h + M / 2, M = sup_t |Phi''(t)|.
Derivatives of the rank-one projections: Q_j' = h_j R_j and Q_j'' = h_j^2 R_j' with R_j, R_j'/2 symmetries of the
block, so ||Y'|| <= (1 - b) alpha, ||Y'||_1 <= 2 (1 - b) alpha, ||Y''||_1 <= 4 (1 - b) alpha^2,
alpha = sum_j |h_theta_j|. For rho = e e^*: d_i e has norm prod_{q<i} sin(phi_q), d_i d_i e = -(the part of e from
index i - 1 on) the same norm, and d_i d_j e (i < j) has norm |cos(phi_i)| prod_{q<j, q != i} sin(phi_q); on a box with
phi_q <= phi_q^max <= pi/2 all three are <= w_min(i,j), w_i = prod_{q<i} sin(phi_q^max). Hence
|e'| <= gamma_1 = sum_i w_i |h_phi_i|, |e''| <= gamma_2 = sum_{i,j} w_min(i,j) |h_phi_i| |h_phi_j|, and as e' is
orthogonal to e, ||rho'|| = |e'| <= gamma_1, ||rho'||_1 = 2 |e'| <= 2 gamma_1, ||rho''||_1 <= 2 |e''| + 2 |e'|^2 <=
2 gamma_2 + 2 gamma_1^2; X' = |X|' = -(1 - a) rho'. Bounding every term of Phi'' with |Tr(A B)| <= ||A||_1 ||B|| and
||X||, ||Y||, |||X|||, ||T|| <= 1 (Tr(|X|'' Y^2) <= (1-a) ||rho''||_1; 2 Tr(|X|' (Y'Y + YY')) <=
8 (1-a)(1-b) alpha gamma_1; Tr(|X| (Y''Y + 2 Y'^2 + Y Y'')) <= 8 (1-b) alpha^2 + 4 (1-b)^2 alpha^2; the same three
bounds for Tr(T Y X Y), with X in place of |X|):
    M <= (16 (1 - b) + 8 (1 - b)^2) alpha^2 + 4 (1 - a)(gamma_2 + gamma_1^2) + 16 (1 - a)(1 - b) alpha gamma_1,
(curvature_bound; the weights w_i matter for m >= 3: phi_2 is irrelevant where phi_1 = 0) and, with |h_i| <= r_i,
    sup_B z <= UB(B) = Phi(c) + sum_i |d_i Phi(c)| r_i + M/2 + EPS.
The search keeps the best value z(c) found (a lower bound L for Z_m: each centre is admissible), splits every box with
UB(B) > L + tol in half along the coordinate with the largest share of UB(B) - Phi(c), and accepts the others; the
accepted boxes cover the domain, so L <= Z_m <= max UB <= L + tol.
Rounding (EPS). Phi(c) and grad Phi(c) are computed in double precision from n x n matrices (n = 2m <= 6) with entries
of modulus <= 2 and absolute errors eta <= 4u (u = 2^-53; at most four operations on libm sines and cosines, which
are accurate to 1 ulp), by products of at most four factors and a trace. By the standard error bound for
floating-point matrix products, |fl(M1 ... Mj) - M1 ... Mj| <= ((j - 1) gamma_n + j eta) |M1| ... |Mj| entrywise
(gamma_n = n u/(1 - n u)), so each computed trace is within n * n^3 * 2^4 * (3 gamma_6 + 16 u) < 1e-10 of the exact
one; UB(B) contains at most 2 + 4 (2m - 1) <= 22 such traces (each gradient term is multiplied by r_i <= pi/4 < 1).
||T_c|| <= 1 is enforced by dividing T_c by sqrt(1 + ||T_c^T T_c - 1||_F + 1e-12) (||T||^2 = ||T^T T|| <=
1 + ||T^T T - 1||_F; the slack 1e-12 covers the rounding of this computation). EPS = 1e-8 is therefore a safe margin,
below the tolerances used (>= 1e-7).

Remark (convex relaxations are not enough). z is concave in |X| (Lemma 2 (a): a minimum of functions affine in
(X_+, X_-)), so averaging e e^* over the sign changes of the blocks gives the valid bound relaxation_bound; it exceeds
the direct sum by up to s0/3 (k = 4, delta -> 1), because it forgets that e e^* has rank one. The certificate
therefore uses branch and bound.

    python scripts/certify_envelopes.py                         # table for k = 2, 3, 4 (about 1 s without --search)
    python scripts/certify_envelopes.py --search --continuum 1000
    python scripts/certify_envelopes.py --k6 --tol 1e-5         # also k = 6 (minutes)
"""
import argparse
import math
import time

import numpy as np

DMAX = 0.999                    # upper end of the free contrast (as in bounds_search.py / optimize)
EPS = 1e-8                      # rounding margin of the branch and bound (docstring)
BOX_LIMIT = 10 ** 9               # maximal number of boxes evaluated


def lo(delta):
    """lambda_min = (1 - delta)/(1 + delta) for lambda_max = 1 and Michelson contrast delta."""
    return (1 - np.asarray(delta, dtype=float)) / (1 + np.asarray(delta, dtype=float))


def s1(delta, D=DMAX):
    """k = 2 upper envelope of Exp. 1 (Theorem 1): delta^2 D^2 / ((1 + delta)(1 + D))."""
    delta = np.asarray(delta, dtype=float)
    return delta ** 2 * D ** 2 / ((1 + delta) * (1 + D))


def s3(delta, D=DMAX):
    """k = 2 upper envelope of Exp. 3 along Delta(|X|) (Lemma 3): (1 - delta) (2D/(1 + D))^2."""
    return (1 - np.asarray(delta, dtype=float)) * (2 * D / (1 + D)) ** 2


def s0(D=DMAX):
    """k = 2 value of Exp. 3 at delta = 0, the supremum of one block with |X| = 1: (2D/(1 + D))^2 = (1 - b)^2."""
    return (2 * D / (1 + D)) ** 2


def theta_star(delta, D=DMAX):
    """Optimal block angle of Lemma 3: cos^2 theta* = (1 + a b)/((1 + a)(1 + b))."""
    a, b = lo(delta), lo(D)
    return np.arccos(np.sqrt((1 + a * b) / ((1 + a) * (1 + b))))


def direct_sum(exp, k, delta, D=DMAX, C=1):
    """Value attained by the direct sums of the paper, check (ii), extended to odd k by a 1-dimensional summand
    (C sectors of size k, envelope along the first contrast): Exp 1 C floor(k/2) s1; Exp 3 s3 + (C k/2 - 1) s0 for k
    even and C floor(k/2) s0 (blocks at delta' = 0, the eigenvalue lambda_min on a 1-dimensional summand) for k odd."""
    m = k // 2
    if exp == 1:
        return C * m * s1(delta, D)
    if k % 2 == 0:
        return s3(delta, D) + (C * m - 1) * s0(D)
    return C * m * s0(D) + 0 * np.asarray(delta, dtype=float)


# ---------------------------------------------------------------------------------------------------------------
# the m-block model of Exp. 3 and its branch and bound
# ---------------------------------------------------------------------------------------------------------------

def unit_vector(phi):
    """e(phi) in R^m (m = phi.shape[1] + 1), hyperspherical coordinates, and its partial derivatives (N, m-1, m):
    e_l = sin(phi_1) ... sin(phi_l) cos(phi_{l+1}), e_{m-1} = sin(phi_1) ... sin(phi_{m-1})."""
    phi = np.asarray(phi, dtype=float)
    N, m1 = phi.shape
    m = m1 + 1
    c, s = np.cos(phi), np.sin(phi)

    def vec(cc, ss):
        e = np.ones((N, m))
        for l in range(m1):
            e[:, l] *= cc[:, l]
            e[:, l + 1:] *= ss[:, l:l + 1]
        return e
    e = vec(c, s)
    de = np.zeros((N, m1, m))
    for i in range(m1):                       # d/dphi_i: replace (cos, sin) of phi_i by (-sin, cos)
        cc, ss = c.copy(), s.copy()
        cc[:, i], ss[:, i] = -s[:, i], c[:, i]
        d = vec(cc, ss)
        d[:, :i] = 0.0
        de[:, i] = d
    return e, de


def block_model(theta, phi, a, b):
    """
    Matrices of the m-block model (Theorem 4), batched: theta (N, m) block angles, phi (N, m - 1) the
    hyperspherical angles of e = sum_j e_j f_j (unit_vector). Basis (f_1, f_1', f_2, f_2', ...).
    Returns X, |X|, Y and the derivatives dY/dtheta_j (N, m, n, n), dX/dphi_i (N, m - 1, n, n).
    """
    theta = np.asarray(theta, dtype=float)
    N, m = theta.shape
    n = 2 * m
    c2, s2 = np.cos(2 * theta), np.sin(2 * theta)
    Q = np.zeros((N, n, n))
    dQ = np.zeros((N, m, n, n))
    for j in range(m):
        i = 2 * j
        Q[:, i, i], Q[:, i + 1, i + 1] = (1 + c2[:, j]) / 2, (1 - c2[:, j]) / 2
        Q[:, i, i + 1] = Q[:, i + 1, i] = s2[:, j] / 2
        dQ[:, j, i, i], dQ[:, j, i + 1, i + 1] = -s2[:, j], s2[:, j]
        dQ[:, j, i, i + 1] = dQ[:, j, i + 1, i] = c2[:, j]
    Y = b * np.eye(n) + (1 - b) * Q
    dY = (1 - b) * dQ
    em, dem = unit_vector(np.zeros((N, 0)) if m == 1 else phi)
    e, de = np.zeros((N, n)), np.zeros((N, m - 1, n))
    e[:, 0::2], de[:, :, 0::2] = em, dem
    S = np.diag(np.tile([1.0, -1.0], m))
    rho = e[:, :, None] * e[:, None, :]
    X = S - (1 - a) * rho
    absX = np.eye(n) - (1 - a) * rho
    dX = -(1 - a) * (de[:, :, :, None] * e[:, None, None, :] + e[:, None, :, None] * de[:, :, None, :])
    return X, absX, Y, dY, dX


def z3(X, Y):
    """Exp. 3 functional Tr(Y|X|Y) - Tr|YXY| for symmetric X, Y (batched, NumPy, independent of the library)."""
    w, V = np.linalg.eigh(X)
    absX = (V * np.abs(w)[..., None, :]) @ np.swapaxes(V, -1, -2)
    M = Y @ X @ Y
    return np.einsum('...ii->...', Y @ absX @ Y) - np.abs(np.linalg.eigvalsh((M + np.swapaxes(M, -1, -2)) / 2)).sum(-1)


def _tr(A):
    return np.einsum('...ii->...', A)


def curvature_bound(r_theta, r_phi, phi_hi, a, b):
    """
    Bound M on |Phi''| along every segment c + t h, t in [0, 1], |h_i| <= r_i, inside a box whose phi-coordinates
    lie in [0, phi_hi] (subset of [0, pi/2]) (docstring):
        M = (16 (1 - b) + 8 (1 - b)^2) alpha^2 + 4 (1 - a)(gamma_2 + gamma_1^2) + 16 (1 - a)(1 - b) alpha gamma_1,
    alpha = sum_j r_theta_j, gamma_1 = sum_i w_i r_phi_i, gamma_2 = sum_{i,j} w_min(i,j) r_phi_i r_phi_j,
    w_i = prod_{q < i} sin(phi_hi_q). Returns M (N,) and a split of M/2 over the coordinates (N, 2m - 1).
    """
    N, m = r_theta.shape
    w = np.ones_like(r_phi)
    if m > 2:
        w[:, 1:] = np.cumprod(np.sin(np.minimum(phi_hi[:, :-1], np.pi / 2)), axis=1)
    wr = w * r_phi
    alpha, g1 = r_theta.sum(1), wr.sum(1)
    tail = np.cumsum(r_phi[:, ::-1], axis=1)[:, ::-1]     # sum_{j >= i} r_j
    g2 = (wr * (2 * tail - r_phi)).sum(1)                  # sum_i w_i r_i (r_i + 2 sum_{j > i} r_j)
    ct = 16 * (1 - b) + 8 * (1 - b) ** 2
    M = ct * alpha ** 2 + 4 * (1 - a) * (g2 + g1 ** 2) + 16 * (1 - a) * (1 - b) * alpha * g1
    share = np.concatenate([(ct * alpha / 2 + 8 * (1 - a) * (1 - b) * g1)[:, None] * r_theta,
                            (2 * (1 - a) * (2 * tail - r_phi + g1[:, None]) + 8 * (1 - a) * (1 - b) * alpha[:, None])
                            * wr], axis=1)
    return M, share


def box_bounds(lo_, hi_, a, b):
    """
    For boxes [lo_, hi_] (N, 2m - 1) of the m-block model (coordinates theta_1..theta_m, phi_1..phi_{m-1}): the
    value z at the centres (attained, a lower bound for Z_m), the certified upper bound UB over every box
    (docstring), and the per-coordinate contributions to UB used to choose the splitting direction.
    """
    c, r = (lo_ + hi_) / 2, (hi_ - lo_) / 2
    m = (c.shape[1] + 1) // 2
    X, absX, Y, dY, dX = block_model(c[:, :m], c[:, m:], a, b)
    M = Y @ X @ Y
    M = (M + np.swapaxes(M, -1, -2)) / 2
    w, V = np.linalg.eigh(M)
    Y2 = Y @ Y
    zc = _tr(absX @ Y2) - np.abs(w).sum(-1)
    T = (V * np.sign(w)[:, None, :]) @ np.swapaxes(V, -1, -2)
    defect = np.linalg.norm(np.swapaxes(T, -1, -2) @ T - np.eye(T.shape[-1]), axis=(-2, -1))
    T = T / np.sqrt(1 + defect + 1e-12)[:, None, None]
    phi_c = _tr(absX @ Y2) - _tr(T @ M)                    # Phi(c) >= z(c) (T renormalised)
    grads = []
    for j in range(m):
        Yj = dY[:, j]
        grads.append(_tr(absX @ (Yj @ Y + Y @ Yj)) - _tr(T @ (Yj @ X @ Y + Y @ X @ Yj)))
    for i in range(m - 1):
        Xi = dX[:, i]
        grads.append(_tr(Xi @ Y2) - _tr(T @ (Y @ Xi @ Y)))
    G = np.abs(np.stack(grads, axis=1))
    Mb, share = curvature_bound(r[:, :m], r[:, m:], hi_[:, m:], a, b)
    ub = phi_c + (G * r).sum(1) + Mb / 2 + EPS
    contrib = G * r + share                                # share of every coordinate in UB - Phi(c)
    return zc, ub, contrib


def certify_blocks(m, delta, D=DMAX, tol=1e-6, max_boxes=BOX_LIMIT, chunk=100_000):
    """
    Branch and bound for Z_m (m = 1: the k = 2 envelope of Exp. 3 along Delta(|X|); m = 2, 3: the open cases
    k = 4, 6).
    Returns dict(lower, upper, boxes, seconds): lower is the best attained value found, upper the certified bound;
    lower <= Z_m <= upper <= lower + tol (+ EPS).
    """
    t0 = time.time()
    a, b = float(lo(delta)), float(lo(D))
    # domain (docstring): theta_j in [0, pi/2], e in the positive orthant (phi_i in [0, pi/2]); for m = 2 the
    # exchange of the blocks halves it to phi in [0, pi/4]
    lo_ = np.zeros((1, 2 * m - 1))
    hi_ = np.full((1, 2 * m - 1), np.pi / 2)
    if m == 2:
        hi_[0, 2] = np.pi / 4
    start = np.concatenate([[theta_star(delta, D)], np.full(m - 1, np.pi / 4), np.zeros(m - 1)])[None]
    lower = float(box_bounds(start, start, a, b)[0][0])          # the direct sum: an attained value
    upper, boxes = -np.inf, 0
    stack = [(lo_, hi_)]                                  # depth first in chunks: bounded memory
    while stack:
        lo_, hi_ = stack.pop()
        if len(lo_) > chunk:
            stack.append((lo_[chunk:], hi_[chunk:]))
            lo_, hi_ = lo_[:chunk], hi_[:chunk]
        boxes += len(lo_)
        zc, ub, contrib = box_bounds(lo_, hi_, a, b)
        lower = max(lower, float(zc.max()))
        split = ub > lower + tol
        if (~split).any():
            upper = max(upper, float(ub[~split].max()))
        lo_, hi_, contrib = lo_[split], hi_[split], contrib[split]
        if not len(lo_):
            continue
        if boxes > max_boxes:
            raise RuntimeError(f"more than {max_boxes} boxes; increase tol")
        d = np.argmax(contrib, axis=1)
        idx = np.arange(len(lo_))
        mid = (lo_[idx, d] + hi_[idx, d]) / 2
        lo2, hi1 = lo_.copy(), hi_.copy()
        hi1[idx, d] = mid
        lo2[idx, d] = mid
        stack.append((np.concatenate([lo_, lo2]), np.concatenate([hi1, hi_])))
    # every accepted box satisfies ub <= lower(at its time) + tol <= lower + tol
    return dict(lower=lower, upper=max(upper, lower), boxes=boxes, seconds=time.time() - t0)


# ---------------------------------------------------------------------------------------------------------------
# certified envelopes
# ---------------------------------------------------------------------------------------------------------------

def certified_sup(exp, k, delta, D=DMAX, tol=1e-6, C=1):
    """
    Certified sup z for C sectors of size k, envelope along the prescribed contrast delta of X (|X| for Exp. 3),
    the other contrast free in [0, D]. Returns dict(lower, upper, method): lower is attained (a configuration with
    this value exists, up to rounding), upper is a proven bound (Theorems 1 and 4, or the branch and bound).
    """
    m = k // 2
    if exp == 1:
        v = float(C * m * s1(delta, D))
        return dict(lower=v, upper=v, method='Theorem 1 (exact)')
    if exp != 3:
        raise ValueError("Exp. 2 is known exactly (paper); only Exp. 1 and 3 here")
    rest = float((C - 1) * m * s0(D))                    # the sectors without e (Theorem 4 (ii))
    if delta == 0 or k % 2:
        v = float(C * m * s0(D))
        return dict(lower=v, upper=v, method='Theorem 4 (exact)')
    if k == 2:
        v = rest + float(s3(delta, D))
        return dict(lower=v, upper=v, method='Lemma 3 / Theorem 4 (exact)')
    if k in (4, 6):
        r = certify_blocks(m, delta, D, tol)
        return dict(lower=rest + r['lower'], upper=rest + max(r['upper'], float((m - 1) * s0(D))), boxes=r['boxes'],
                    seconds=r['seconds'], method='Theorem 4 + branch and bound')
    raise ValueError("certified here for even k <= 6 (k odd: any k)")


def certify_all_contrasts(n=1000, D=DMAX, tol=1e-6, dmax=DMAX):
    """
    Exp. 3, k = 4, C = 1: a certificate for EVERY contrast delta in [0, dmax]. sup z is non-increasing in delta
    (for delta <= delta', a >= a', and |X''| = 1 - (1 - a) e e^* >= 1 - (1 - a') e e^*, so Lemma 2 (b) maps every
    configuration at delta' to one at delta with at least the same z). Hence on [delta_i, delta_{i+1}],
    dsum(delta) <= sup z(delta) <= UB(delta_i), and the gap to the direct sums is at most
    max_i UB(delta_i) - dsum(delta_{i+1}) (the direct sums dsum are non-increasing as well).
    Returns dict(gap, seconds, boxes).
    """
    t0 = time.time()
    grid = np.linspace(0.0, dmax, n + 1)
    ub = np.array([certified_sup(3, 4, d, D, tol)['upper'] for d in grid])
    ds = direct_sum(3, 4, grid, D)
    return dict(gap=float(np.max(ub[:-1] - ds[1:])), pointwise_gap=float(np.max(ub - ds)), seconds=time.time() - t0)


def relaxation_bound(m, delta, D=DMAX):
    """
    A valid but loose bound for the m-block model, from concavity alone: z is concave in |X| (Lemma 2 (a): a
    minimum of functions affine in (X_+, X_-)); averaging e e^* over the sign changes of the blocks (which commute
    with S and Q) gives block-diagonal |X| with eigenvalues a_j = 1 - (1 - a) c_j^2 on f_j, so by Lemma 3
    z <= (1 - b)^2 sum_j 2 a_j / (1 + a_j) <= m (1 - b)^2 h(1 - (1 - a)/m), h(x) = 2x/(1 + x) concave (Jensen).
    It exceeds the direct sum by up to s0/3 (m = 2, delta -> 1): the 'pure state' constraint on e e^* matters, which
    is why the certificate uses branch and bound rather than a convex relaxation.
    """
    a = lo(delta)
    x = 1 - (1 - a) / m
    return m * s0(D) * 2 * x / (1 + x)


# ---------------------------------------------------------------------------------------------------------------
# multistart search (the method of bounds_search.py), for comparison
# ---------------------------------------------------------------------------------------------------------------

def search_sup(exp, k, deltas, starts=16, rounds=4, steps=40, seed=0, D=DMAX):
    """Upper envelope by the multistart gradient search of torch_vn_algebra.optimize (real algebra, C = 1, as in
    scripts/bounds_search.py): a lower bound for sup z at every contrast in `deltas`."""
    import torch
    from torch_vn_algebra import TypeIAlgebra, optimize as opt
    torch.manual_seed(seed)
    alg = TypeIAlgebra([k], [k], complex_valued=False, precision='double', device='cpu')
    G = len(deltas)
    d = torch.tensor(np.asarray(deltas, dtype=float), dtype=torch.float64).repeat_interleave(starts)
    if exp == 3:
        X = opt.SelfAdjointParam(alg, G * starts, delta=d)
    else:
        X = opt.PositiveParam(alg, G * starts, delta=d)
    Y = opt.PositiveParam(alg, G * starts, delta_range=(0.0, D))
    if exp == 1:
        U = opt.UnitaryParam(alg, G * starts)

        def z():
            x, y = X.operator(), Y.operator()
            return (x @ U.operator() @ y).trace.abs() - (x @ y).trace.real
        params = [X, Y, U]
    else:
        def z():
            x, y = X.operator(), Y.operator()
            M = y @ x @ y
            M._is_self_adjoint = True
            return (y @ X.abs_operator() @ y).trace.real - M.trace_norm()
        params = [X, Y]
    vals = opt.extremize(z, params, maximize=True, rounds=rounds, steps=steps).reshape(G, starts)
    vals = torch.nan_to_num(vals, nan=-math.inf)
    return vals.max(dim=1)[0].numpy()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--grid', type=float, nargs='*', default=[0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.999])
    p.add_argument('--tol', type=float, default=1e-7, help='branch-and-bound gap (absolute)')
    p.add_argument('--search', action='store_true', help='also run the multistart search for comparison')
    p.add_argument('--k6', action='store_true', help='also k = 6 (Exp. 3: branch and bound in 5 angles, minutes)')
    p.add_argument('--continuum', type=int, default=0, metavar='N',
                   help='Exp. 3, k = 4: certify all contrasts in [0, 0.999] with N + 1 grid points')
    args = p.parse_args()
    print(f"C = 1; envelopes along the prescribed contrast (Delta(X), resp. Delta(|X|)); the other contrast free in "
          f"[0, {DMAX}]")
    print(f"{'exp':>3s} {'k':>2s} {'delta':>6s} {'direct sum':>12s} {'certified UB':>13s} {'UB-dsum':>9s} "
          f"{'relaxation':>11s} {'search':>10s} {'boxes':>7s} {'time[s]':>7s}  method")
    t_all = time.time()
    for exp in (1, 3):
        for k in (2, 3, 4) + ((6,) if args.k6 else ()):
            found = search_sup(exp, k, args.grid) if args.search else [float('nan')] * len(args.grid)
            for delta, sv in zip(args.grid, found):
                r = certified_sup(exp, k, delta, tol=args.tol)
                ds = float(direct_sum(exp, k, delta))
                rel = float(max((k // 2 - 1) * s0(), relaxation_bound(k // 2, delta))) \
                    if exp == 3 and k % 2 == 0 and k > 2 else float('nan')
                print(f"{exp:3d} {k:2d} {delta:6.3f} {ds:12.8f} {r['upper']:13.8f} {r['upper'] - ds:9.1e} "
                      f"{rel:11.6f} {sv:10.6f} {r.get('boxes', 0):7d} {r.get('seconds', 0.0):7.2f}  {r['method']}",
                      flush=True)
    if args.continuum:
        c = certify_all_contrasts(args.continuum, tol=args.tol)
        print(f"Exp 3, k = 4, all delta in [0, {DMAX}] ({args.continuum + 1} grid points): sup z - direct sum <= "
              f"{c['gap']:.2e} (at the grid points <= {c['pointwise_gap']:.2e}); {c['seconds']:.1f} s")
    print(f"total {time.time() - t_all:.1f} s")


if __name__ == '__main__':
    main()
