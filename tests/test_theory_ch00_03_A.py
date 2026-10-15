"""
Theory tests for the lecture notes, Chapters 0-3 and Appendix A (man/lecture_notes/{00,01,02,03,A}_*.md), through
the library: every theorem, proposition, lemma, example and exercise with a computable consequence is checked on
exact small cases, random and adversarial inputs, and (for distributional statements) by Monte Carlo with a z-score
bound. Deterministic (fixed seeds), CPU, double precision.

Run with: OMP_NUM_THREADS=2 python -m pytest -q tests/test_theory_ch00_03_A.py
"""
import itertools
import math

import numpy as np
import pytest
import torch

import torch_vn_algebra as tv
from torch_vn_algebra import TypeIAlgebra, DensityMatrix, channels, states, decompose

DT = torch.complex128
RDT = torch.float64

# the cost module's run-time estimate is a heuristic (it can overestimate small batched QR calls); not under test here
pytestmark = pytest.mark.filterwarnings('ignore::torch_vn_algebra.cost.CostWarning')

# sector layouts: one block, k = 1, C = 1, equal blocks, unequal blocks with padding
LAYOUTS = [[1], [3], [1, 1, 1], [2, 3, 5], [4, 4], [1, 4, 2]]
LAYOUT_IDS = ['-'.join(map(str, k)) for k in LAYOUTS]


def alg_of(ks, complex_valued=True):
    return TypeIAlgebra(list(ks), list(ks), precision='double', complex_valued=complex_valued)


def seed(s):
    torch.manual_seed(s)
    return np.random.default_rng(s)


def blocks(op, b=0):
    """The active blocks A_c of sample b of an operator (list of tensors)."""
    m = op.matrix if hasattr(op, 'matrix') else op
    return [m[b, c, :k, :k] for c, k in enumerate(op.algebra.k_factors)]


def dense(op, mult=None, b=0):
    """The element (+)_c A_c (x) 1_{m_c} of M_n (multiplicities m_c, default 1)."""
    ks = op.algebra.k_factors
    mult = mult or [1] * len(ks)
    return torch.block_diag(*[torch.kron(A, torch.eye(m, dtype=A.dtype)) for A, m in zip(blocks(op, b), mult)])


def maxdiff(X, Y):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item()


def ginibre(alg, B=1):
    """Operator with independent standard complex Gaussian entries in every block (non-normal)."""
    return alg.operator(torch.randn(B, alg.C, alg.k_max, alg.k_max, dtype=DT))


def hermitian(alg, B=1):
    G = torch.randn(B, alg.C, alg.k_max, alg.k_max, dtype=DT)
    return alg.operator((G + G.mH) / 2, is_self_adjoint=True)


def positive(alg, B=1, shift=0.0):
    G = torch.randn(B, alg.C, alg.k_max, alg.k_max, dtype=DT)
    eye = alg.identity(B).matrix
    return alg.operator(G @ G.mH + shift * eye, is_self_adjoint=True, is_positive=True)


def central_projection(alg, c, B=1):
    return alg.central([1.0 if d == c else 0.0 for d in range(alg.C)], batch_size=B)


def zscore(samples, exact):
    samples = np.asarray(samples, dtype=float).ravel()
    sd = samples.std(ddof=1)
    if sd == 0:                                    # a constant sample: exact agreement or an infinite z-score
        return 0.0 if abs(samples.mean() - exact) < 1e-12 else math.inf
    return (samples.mean() - exact) / (sd / math.sqrt(samples.size))


# =====================================================================================================================
# Chapter 0: the table of Section 0.1 (each row an identity)
# =====================================================================================================================

@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_ch0_row1_products_sums_adjoints_are_blockwise(ks):
    """Section 0.1, row 1: product, sum and adjoint of M are taken block by block (phi is a *-homomorphism)."""
    seed(1)
    alg = alg_of(ks)
    A, B = ginibre(alg, 3), ginibre(alg, 3)
    for b in range(3):
        assert maxdiff(dense(A @ B, b=b), dense(A, b=b) @ dense(B, b=b)) < 1e-12
        assert maxdiff(dense(A + B, b=b), dense(A, b=b) + dense(B, b=b)) < 1e-12
        assert maxdiff(dense(A.adjoint(), b=b), dense(A, b=b).mH) < 1e-15
        assert maxdiff(dense(A * (2 - 3j), b=b), (2 - 3j) * dense(A, b=b)) < 1e-14
    assert maxdiff((A @ B).matrix * ~alg._const('block_mask', torch.bool), torch.zeros(1)) == 0   # padding stays 0


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
@pytest.mark.parametrize('mult', ['ones', 'many'])
def test_ch0_row2_functional_calculus_is_blockwise(ks, mult):
    """Section 0.1, row 2: phi(f(A)) = (+)_c f(A_c), also for the algebra acting with multiplicities m_c > 1."""
    seed(2)
    alg = alg_of(ks)
    m = [1] * len(ks) if mult == 'ones' else [c % 3 + 1 for c in range(len(ks))]
    A = positive(alg, 2, shift=0.1)
    for f_lib, f_np in ((lambda X: X.sqrt(), np.sqrt), (lambda X: X.log(), np.log),
                        (lambda X: X.power(-0.7), lambda w: w ** -0.7),
                        (lambda X: X.apply_function(lambda w: torch.sin(3 * w)), lambda w: np.sin(3 * w))):
        F = f_lib(A)
        for b in range(2):
            big = dense(A, m, b).numpy()
            w, V = np.linalg.eigh(big)                               # one large eigenproblem in M_n
            f_big = (V * f_np(w)) @ V.conj().T
            assert np.abs(dense(F, m, b).numpy() - f_big).max() < 1e-10


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_ch0_row2_spectrum_of_blocks_in_spectrum(ks):
    """Section 0.1, proof of row 2: spec A_c is contained in spec A, and spec A is the union of the spec A_c."""
    seed(3)
    alg = alg_of(ks)
    H = hermitian(alg)
    w_blocks = H.eigenvalues()
    w_all = torch.linalg.eigvalsh(dense(H))
    union = torch.sort(torch.cat([w[0] for w in w_blocks])).values
    assert maxdiff(union, w_all) < 1e-12


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_ch0_row3_traces_are_weighted_block_traces(ks):
    """Section 0.1, row 3: the traces sum_c w_c Tr A_c are diagonal sums per block, weighted over the sector axis."""
    seed(4)
    alg = alg_of(ks)
    A = ginibre(alg, 2)
    C = len(ks)
    for b in range(2):
        tr = [torch.trace(X) for X in blocks(A, b)]
        assert abs(A.Tr_blunt()[b] - sum(tr)) < 1e-12
        assert abs(A.Tr_norm()[b] - sum(t / k for t, k in zip(tr, ks))) < 1e-12
        assert abs(A.tau_vN()[b] - sum(t / (C * k) for t, k in zip(tr, ks))) < 1e-12
        assert abs(A.trace[b] - torch.trace(dense(A, b=b))) < 1e-12


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_ch0_row4_centre_is_one_number_per_sector(ks):
    """Section 0.1, row 4: phi(Z(M)) = {(+)_c z_c 1_{k_c}}; central elements commute with M, dim Z(M) = C."""
    seed(5)
    alg = alg_of(ks)
    z = torch.randn(len(ks), dtype=DT)
    Z = alg.central([complex(x) for x in z])
    A = ginibre(alg, 4)
    assert maxdiff(Z @ A, A @ Z) < 1e-12
    for c, X in enumerate(blocks(Z)):
        assert maxdiff(X, z[c] * torch.eye(ks[c], dtype=DT)) < 1e-15
    # the commutant of M inside M (the centre) has dimension C: null space of X -> ([E_ij^c, X])
    basis = [alg.from_blocks([torch.eye(k, dtype=DT)[i][:, None] * torch.eye(k, dtype=DT)[j][None, :] if d == c
                              else torch.zeros(k, k, dtype=DT) for d, k in enumerate(ks)])
             for c, kc in enumerate(ks) for i in range(kc) for j in range(kc)]
    gens = torch.stack([dense(E) for E in basis])
    D = sum(ks)
    L = torch.cat([torch.kron(g, torch.eye(D, dtype=DT)) - torch.kron(torch.eye(D, dtype=DT), g.T.contiguous()) for g in gens])
    s = torch.linalg.svdvals(L)
    dim_comm_MD = int((s < 1e-10).sum())                 # commutant of M in M_D: (+)_c C 1_{k_c}, dimension C
    assert dim_comm_MD == len(ks)


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_ch0_row5_centre_expectation_is_normalised_block_trace(ks):
    """Section 0.1, row 5: the conditional expectation onto the centre is Tr A_c / k_c, broadcast back."""
    seed(6)
    alg = alg_of(ks)
    A = ginibre(alg, 3)
    EA = channels.center_expectation(alg)(A)
    for b in range(3):
        for X, Y, k in zip(blocks(A, b), blocks(EA, b), ks):
            assert maxdiff(Y, torch.trace(X) / k * torch.eye(k, dtype=DT)) < 1e-12


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_ch0_row6_random_unitary_operator_blocks_unitary_and_haar_moments(ks):
    """Section 0.1, row 6 (Corollary 5.4, Theorem 5.3): independent Haar blocks: E|U_ij|^2 = 1/k_c, E U_11 = 0."""
    seed(7)
    alg = alg_of(ks)
    B = 4000
    U = alg.random_unitary_operator(B)
    assert maxdiff(U.adjoint() @ U, alg.identity(B)) < 1e-12
    assert maxdiff(U @ U.adjoint(), alg.identity(B)) < 1e-12
    for c, k in enumerate(ks):
        Uc = U.matrix[:, c, :k, :k]
        assert abs(zscore((Uc[:, 0, 0].abs() ** 2).numpy(), 1 / k)) < 5 or k == 1
        # phase correction of Theorem 5.3: E U_11 = 0 (an uncorrected QR factor has a biased diagonal)
        assert abs(zscore(Uc[:, 0, 0].real.numpy(), 0.0)) < 5
        assert abs(zscore(Uc[:, 0, 0].imag.numpy(), 0.0)) < 5
    if len(ks) >= 2:                               # independent blocks: E[Tr U_0 conj Tr U_1] = 0
        t0 = U.matrix[:, 0].diagonal(dim1=-2, dim2=-1).sum(-1)
        t1 = U.matrix[:, 1].diagonal(dim1=-2, dim2=-1).sum(-1)
        assert abs(zscore((t0 * t1.conj()).real.numpy(), 0.0)) < 5
        assert abs(zscore((t0.abs() ** 2 * t1.abs() ** 2).numpy(), 1.0)) < 5     # E|Tr U_0|^2 |Tr U_1|^2 = 1 * 1


@pytest.mark.parametrize('k', [1, 2, 3, 5])
def test_ch0_row6_haar_trace_moments(k):
    """Section 0.1, row 6 (Theorem 5.3): Haar on U(k): E Tr U^j = 0 and E|Tr U^j|^2 = min(j, k) (Diaconis-Shahshahani)."""
    seed(8)
    alg = alg_of([k])
    U = alg.random_unitary_operator(20000)
    P = U
    for j in (1, 2, 3):
        t = P.trace
        assert abs(zscore(t.real.numpy(), 0.0)) < 5 and abs(zscore(t.imag.numpy(), 0.0)) < 5
        assert abs(zscore((t.abs() ** 2).numpy(), min(j, k))) < 5
        P = P @ U


# =====================================================================================================================
# Chapter 1: matrix algebras and functional calculus
# =====================================================================================================================

@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_lemma_1_1_c_star_identity(ks):
    """Lemma 1.1: ||A^*|| = ||A|| and ||A^*A|| = ||A||^2 (operator norm), on non-normal operators."""
    seed(11)
    alg = alg_of(ks)
    A = ginibre(alg, 5)
    n = A.operator_norm()
    assert (A.adjoint().operator_norm() - n).abs().max() < 1e-12
    assert ((A.adjoint() @ A).operator_norm() - n ** 2).abs().max() < 1e-11 * n.max() ** 2
    # the operator norm of the algebra is the max over the blocks
    for b in range(5):
        assert abs(n[b] - torch.linalg.matrix_norm(dense(A, b=b), ord=2)) < 1e-12


def test_lemma_1_2_and_theorem_1_3_schur_and_spectral():
    """Lemma 1.2 / Theorem 1.3: normal iff unitarily diagonalisable; a non-normal matrix is only triangularisable."""
    seed(12)
    alg = alg_of([2, 3, 5])
    lam = lambda d: torch.randn(3, d, dtype=DT)                    # complex spectrum
    N = alg.operator_from_eigenvalues(lam, batch_size=3)
    assert N.is_normal and not N.is_self_adjoint
    M = N.matrix
    assert maxdiff(M @ M.mH, M.mH @ M) < 1e-12                     # U diag U^* is normal
    G = ginibre(alg)
    assert not G.is_normal                                          # a Ginibre matrix is not normal (a.s.)
    import scipy.linalg as sla
    T, Z = sla.schur(dense(G).numpy(), output='complex')            # Lemma 1.2: U^* A U upper triangular
    assert np.abs(np.tril(T, -1)).max() < 1e-12 and np.abs(Z @ T @ Z.conj().T - dense(G).numpy()).max() < 1e-12
    # Theorem 1.3, "the lambda_i are the eigenvalues with multiplicity"
    w = torch.linalg.eigvals(dense(N, b=1))
    lam_all = torch.cat([N._eigenvalues[c][1] for c in range(3)])
    assert maxdiff(torch.sort(w.real).values, torch.sort(lam_all.real).values) < 1e-12


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_theorem_1_3_operator_from_eigenvalues_has_prescribed_spectrum(ks):
    """Theorem 1.3 read backwards: U diag(lambda) U^* with Haar U has spectrum lambda (recomputed from the matrix)."""
    seed(13)
    alg = alg_of(ks)
    A = alg.operator_from_eigenvalues(lambda d: torch.randn(4, d), batch_size=4, force_self_adjoint=True)
    B = alg.operator(A.matrix.clone())                                # forget the prescribed spectrum
    for c, k in enumerate(ks):
        prescribed = torch.sort(A._eigenvalues[c].real, dim=-1).values
        assert maxdiff(B.eigenvalues()[c], prescribed) < 1e-12
        assert maxdiff(A.eigenvalues()[c], prescribed) == 0


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_corollary_1_4(ks):
    """Corollary 1.4: ||A|| = max|lambda| for normal A; unitary iff |lambda|=1; projection iff lambda in {0,1}; A>=0 iff lambda>=0."""
    seed(14)
    alg = alg_of(ks)
    lam = [torch.randn(2, k, dtype=DT) for k in ks]
    it = iter(lam)
    N = alg.operator_from_eigenvalues(lambda d: next(it), batch_size=2)
    Nn = alg.operator(N.matrix.clone())                              # norm via the SVD path
    expected = torch.stack([torch.cat([l[b].abs() for l in lam]).max() for b in range(2)])
    assert maxdiff(Nn.operator_norm(), expected) < 1e-12
    # unitary iff |lambda| = 1
    Uph = alg.operator_from_eigenvalues(lambda d: torch.exp(2j * math.pi * torch.rand(2, d, dtype=RDT)), batch_size=2)
    assert maxdiff(Uph.adjoint() @ Uph, alg.identity(2)) < 1e-12
    # projection iff lambda in {0, 1}
    P = alg.operator_from_eigenvalues(lambda d: (torch.rand(2, d, dtype=RDT) < 0.5).double(), batch_size=2)
    assert P.is_projection and maxdiff(P @ P, P) < 1e-12 and maxdiff(P.adjoint(), P) < 1e-15
    # positivity: decided from the computed spectrum for an untagged matrix
    H = alg.operator(hermitian(alg, 1).matrix)
    lmin = H.lambda_min.item()
    assert alg.operator((H - alg.identity() * (lmin - 0.01)).matrix).is_positive
    assert not alg.operator((H - alg.identity() * (lmin + 0.01)).matrix).is_positive
    assert not alg.operator(ginibre(alg).matrix).is_positive or sum(ks) == 0


def test_corollary_1_4_positive_implies_self_adjoint_polarisation():
    """Corollary 1.4(3), proof: 4<x,Dy> = sum_k i^-k <x + i^k y, D(x + i^k y)> (polarisation), so <x,Dx>=0 for all x forces D=0."""
    seed(15)
    D = torch.randn(6, 6, dtype=DT)
    x, y = torch.randn(6, dtype=DT), torch.randn(6, dtype=DT)
    q = lambda v: torch.vdot(v, D @ v)
    rhs = sum((1j) ** (-k) * q(x + (1j) ** k * y) for k in range(4))
    assert abs(4 * torch.vdot(x, D @ y) - rhs) < 1e-12
    # a non-Hermitian matrix has a non-real quadratic form somewhere
    alg = alg_of([6])
    assert not alg.operator(D.reshape(1, 1, 6, 6)).is_positive


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_theorem_1_5_degenerate_spectrum_lagrange_polynomial(ks):
    """Theorem 1.5(a),(b): f(A) equals the Lagrange interpolation polynomial p(A), also for degenerate spectra."""
    seed(16)
    alg = alg_of(ks)
    levels = torch.tensor([-1.0, 0.5, 2.0], dtype=RDT)
    A = alg.operator_from_eigenvalues(lambda d: levels[torch.arange(d) % 3], force_self_adjoint=True)
    f = lambda w: torch.exp(w) * torch.cos(w)
    F = A.apply_function(f)
    I = alg.identity()
    p = alg.zero()
    for i, mu in enumerate(levels.tolist()):
        term = I * float(f(torch.tensor(mu, dtype=RDT)))
        for j, nu in enumerate(levels.tolist()):
            if j != i:
                term = term @ ((A - I * nu) * (1.0 / (mu - nu)))
        p = p + term
    assert maxdiff(F, p) < 1e-12
    # (b) polynomials: the functional calculus of t^3 - 2t + 1 is the polynomial in A
    poly = A.apply_function(lambda w: w ** 3 - 2 * w + 1)
    assert maxdiff(poly, A @ A @ A - A * 2.0 + I) < 1e-12


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_theorem_1_5_c_homomorphism_conjugation_spectral_mapping(ks):
    """Theorem 1.5(c): f -> f(A) linear and multiplicative, conj(f)(A) = f(A)^*, spec f(A) = f(spec A)."""
    seed(17)
    alg = alg_of(ks)
    A = alg.operator_from_eigenvalues(lambda d: torch.randn(3, d), batch_size=3, force_self_adjoint=True)
    f = lambda w: torch.exp(1j * w) + w ** 2
    g = lambda w: torch.sin(w) - 2j
    Fa, Ga = A.apply_function(f), A.apply_function(g)
    assert maxdiff(A.apply_function(lambda w: f(w) * g(w)), Fa @ Ga) < 1e-12           # multiplicative
    assert maxdiff(A.apply_function(lambda w: 2 * f(w) - 3j * g(w)), Fa * 2 - Ga * 3j) < 1e-12
    assert maxdiff(A.apply_function(lambda w: torch.conj(f(w))), Fa.adjoint()) < 1e-12
    S = A.apply_function(torch.sinh)
    for c in range(len(ks)):
        lam = torch.sort(A._eigenvalues[c].real, dim=-1).values
        assert maxdiff(alg.operator(S.matrix.clone()).eigenvalues()[c],
                       torch.sort(torch.sinh(lam), dim=-1).values) < 1e-12


def test_theorem_1_5_d_commutant_of_degenerate_operator():
    """Theorem 1.5(d): B commuting with A (not a function of A: it mixes a degenerate eigenspace) commutes with f(A)."""
    seed(18)
    alg = alg_of([5])
    U = alg.random_unitary(5).reshape(1, 1, 5, 5)
    Lam = torch.diag(torch.tensor([1.0, 1.0, 1.0, 3.0, -2.0], dtype=DT))
    A = alg.operator(U @ Lam @ U.mH, is_self_adjoint=True)
    blk = torch.zeros(5, 5, dtype=DT)
    blk[:3, :3] = torch.randn(3, 3, dtype=DT)                         # arbitrary on the 3-dim eigenspace
    blk[3, 3], blk[4, 4] = 2.0, -1.0 + 1j
    Bm = alg.operator(U @ blk @ U.mH)
    assert maxdiff(A @ Bm, Bm @ A) < 1e-12
    for f in (torch.exp, lambda w: w ** 4 - w, torch.atan):
        F = A.apply_function(f)
        assert maxdiff(F @ Bm, Bm @ F) < 1e-11


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_theorem_1_5_e_composition(ks):
    """Theorem 1.5(e): g(f(A)) = (g o f)(A): exp(log A) = A, (A^p)^q = A^{pq}, sqrt(A^2) = A, log(exp H) = H."""
    seed(19)
    alg = alg_of(ks)
    A = positive(alg, 2, shift=0.2)
    assert maxdiff(A.log().expm(), A) < 1e-10
    assert maxdiff(A.power(0.3).power(2.5), A.power(0.75)) < 1e-10
    assert maxdiff((A @ A).sqrt(), A) < 1e-9
    H = hermitian(alg, 2)
    assert maxdiff(H.expm().log(), H) < 1e-10
    assert maxdiff(H.apply_function(torch.exp), H.expm()) < 1e-10     # expm (Pade) = functional calculus


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_theorem_1_6_unique_positive_square_root(ks):
    """Theorem 1.6: sqrt(A) >= 0, sqrt(A)^2 = A, and every B >= 0 with B^2 = A is sqrt(A)."""
    seed(20)
    alg = alg_of(ks)
    B = positive(alg, 3)                                              # some positive B (with small eigenvalues)
    A = alg.operator((B @ B).matrix, is_self_adjoint=True)
    R = A.sqrt()
    assert maxdiff(R, B) < 1e-8 * max(1.0, B.operator_norm().max().item())
    assert maxdiff(R @ R, A) < 1e-10 * A.operator_norm().max().item()
    assert (alg.operator(R.matrix.clone()).lambda_min >= -1e-10).all()
    # a non-positive square root of A exists (e.g. -B, or B with a sign flipped on an eigenvector) but is not sqrt(A)
    assert maxdiff((B * -1.0) @ (B * -1.0), A) < 1e-10 * A.operator_norm().max().item()


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_corollary_1_7_TstarT_positive(ks):
    """Corollary 1.7: A >= 0 iff A = T^*T: T^*T is positive, and sqrt(A)^* sqrt(A) = A."""
    seed(21)
    alg = alg_of(ks)
    T = ginibre(alg, 3)
    A = alg.operator((T.adjoint() @ T).matrix)
    assert A.is_self_adjoint and A.is_positive
    R = A.sqrt()
    assert maxdiff(R.adjoint() @ R, A) < 1e-10 * A.operator_norm().max().item()


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_theorem_1_8_polar_decomposition(ks):
    """Theorem 1.8: A = W|A| with W unitary, |A| = (A^*A)^{1/2} (Operator.abs and Operator.polar)."""
    seed(22)
    alg = alg_of(ks)
    A = ginibre(alg, 3)
    absA = A.abs()
    assert maxdiff(absA @ absA, A.adjoint() @ A) < 1e-10
    assert maxdiff(absA, (A.adjoint() @ A).sqrt()) < 1e-10
    W, P = A.polar()
    assert maxdiff(P, absA) < 1e-10
    assert maxdiff(W.adjoint() @ W, alg.identity(3)) < 1e-10 and maxdiff(W @ W.adjoint(), alg.identity(3)) < 1e-10
    assert maxdiff(W @ P, A) < 1e-10
    # || |A| x || = || A x ||
    x = torch.randn(3, alg.C, alg.k_max, 1, dtype=DT) * alg._const('rows')[None, :, :, None]
    assert maxdiff(torch.linalg.vector_norm(absA.matrix @ x, dim=(-3, -2, -1)),
                   torch.linalg.vector_norm(A.matrix @ x, dim=(-3, -2, -1))) < 1e-10


def test_theorem_1_8_polar_singular():
    """Theorem 1.8 for a singular A (rank 2 in M_4): still A = W|A| (the library's W is a partial isometry there)."""
    seed(23)
    alg = alg_of([4])
    X = torch.randn(1, 1, 4, 2, dtype=DT) @ torch.randn(1, 1, 2, 4, dtype=DT)
    A = alg.operator(X)
    W, P = A.polar()
    assert maxdiff(W @ P, A) < 1e-10
    assert maxdiff(P @ P, A.adjoint() @ A) < 1e-10
    assert abs(torch.linalg.matrix_rank(P.matrix[0, 0]).item() - 2) == 0


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_section_1_4_three_norms(ks):
    """Section 1.4: ||A|| = s_1, ||A||_2 = (sum s_i^2)^{1/2} = (Tr A^*A)^{1/2}, ||A||_1 = sum s_i = Tr|A| (padding ignored)."""
    seed(24)
    alg = alg_of(ks)
    A = ginibre(alg, 4)
    for b in range(4):
        s = torch.linalg.svdvals(dense(A, b=b))
        assert abs(A.operator_norm()[b] - s[0]) < 1e-12
        assert abs(A.frobenius_norm()[b] - s.norm()) < 1e-12
        assert abs(A.trace_norm()[b] - s.sum()) < 1e-11
    assert maxdiff(A.frobenius_norm() ** 2, (A.adjoint() @ A).Tr_blunt().real) < 1e-11
    assert maxdiff(A.trace_norm(), A.abs().Tr_blunt().real) < 1e-10
    # norms: triangle inequality, homogeneity, ||P||_1 = Tr P for P >= 0
    B = ginibre(alg, 4)
    for nrm in ('operator_norm', 'frobenius_norm', 'trace_norm'):
        lhs = getattr(A + B, nrm)()
        assert (lhs <= getattr(A, nrm)() + getattr(B, nrm)() + 1e-12).all()
        assert maxdiff(getattr(A * (3 - 4j), nrm)(), 5 * getattr(A, nrm)()) < 1e-10
    P = positive(alg, 4)
    assert maxdiff(P.trace_norm(), P.Tr_blunt().real) < 1e-10
    # ordering s_1 <= ||A||_2 <= ||A||_1 <= sqrt(D) ||A||_2
    D = sum(ks)
    assert (A.operator_norm() <= A.frobenius_norm() + 1e-12).all()
    assert (A.frobenius_norm() <= A.trace_norm() + 1e-12).all()
    assert (A.trace_norm() <= math.sqrt(D) * A.frobenius_norm() + 1e-12).all()


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_theorem_1_9_trace_norm_duality(ks):
    """Theorem 1.9: ||A||_1 = max_U |Tr(AU)|, attained at U = W^* (A = W|A|); |Tr A| <= ||A||_1."""
    seed(25)
    alg = alg_of(ks)
    A = ginibre(alg, 2)
    W, _ = A.polar()
    t1 = A.trace_norm()
    assert maxdiff((A @ W.adjoint()).trace, t1.to(DT)) < 1e-10
    assert (A.trace.abs() <= t1 + 1e-12).all()
    U = alg.random_unitary_operator(2000)
    A_rep = alg.operator(A.matrix[:1].expand(2000, -1, -1, -1))
    vals = (A_rep @ U).trace.abs()
    assert (vals <= t1[0] + 1e-10).all()
    assert vals.max() > 0.5 * t1[0] or sum(ks) > 6                  # random unitaries do reach a fair fraction


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_exercise_1_2_trace_norm_adjoint_and_hoelder(ks):
    """Chapter 1, Exercise 2: ||A^*||_1 = ||A||_1 and ||AB||_1 <= ||A|| ||B||_1."""
    seed(26)
    alg = alg_of(ks)
    A, B = ginibre(alg, 5), ginibre(alg, 5)
    assert maxdiff(A.adjoint().trace_norm(), A.trace_norm()) < 1e-11
    assert ((A @ B).trace_norm() <= A.operator_norm() * B.trace_norm() + 1e-10).all()
    assert ((B @ A).trace_norm() <= A.operator_norm() * B.trace_norm() + 1e-10).all()
    # sharp case: B = |B| a rank-one projection aligned with the top singular vector of A
    if ks == [3]:
        U, s, Vh = torch.linalg.svd(A.matrix[0, 0])
        v = Vh[0].conj()
        Bm = alg.operator(torch.outer(v, v.conj()).reshape(1, 1, 3, 3))
        assert abs((A @ Bm).trace_norm()[0] - s[0]) < 1e-12


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_exercise_1_1_functional_calculus_multiplicative(ks):
    """Chapter 1, Exercise 1: (fg)(A) = f(A) g(A) (here f = sin, g = exp, on a random self-adjoint A)."""
    seed(27)
    alg = alg_of(ks)
    H = hermitian(alg, 3)
    assert maxdiff(H.apply_function(lambda w: torch.sin(w) * torch.exp(w)),
                   H.apply_function(torch.sin) @ H.apply_function(torch.exp)) < 1e-11


def _geometric_mean(A, B):
    a, ai = A.sqrt(), A.power(-0.5)
    return a @ (ai @ B @ ai).sqrt() @ a


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_exercise_1_3_geometric_mean(ks):
    """Chapter 1, Exercise 3: A#B >= 0 solves X A^-1 X = B, A#B = B#A, A#A = A, A#B = sqrt(AB) if AB = BA."""
    seed(28)
    alg = alg_of(ks)
    A, B = positive(alg, 2, shift=0.3), positive(alg, 2, shift=0.3)
    GM = _geometric_mean(A, B)
    assert maxdiff(GM @ A.inverse() @ GM, B) < 1e-9
    assert maxdiff(GM, _geometric_mean(B, A)) < 1e-9
    assert (alg.operator(GM.matrix, is_self_adjoint=True).lambda_min > 0).all()
    assert maxdiff(_geometric_mean(A, A), A) < 1e-9
    # commuting pair: functions of the same A
    A2 = A.apply_function(lambda w: w ** 2 + 1)
    assert maxdiff(_geometric_mean(A, A2), (A @ A2).sqrt()) < 1e-9
    # uniqueness among positive solutions: Riccati X A^-1 X = B has the negative solution -A#B as well
    Xn = GM * -1.0
    assert maxdiff(Xn @ A.inverse() @ Xn, B) < 1e-9


@pytest.mark.parametrize('n', [1, 2, 3, 5])
def test_proposition_1_10_gue_fourth_moment(n):
    """Proposition 1.10: E (1/n) Tr X^4 = 2 + 1/n^2 for X = (G + G^*)/sqrt(2n) (GUE); also E (1/n) Tr X^2 = 1."""
    seed(30 + n)
    B = 40000
    alg = alg_of([n])
    G = torch.randn(B, 1, n, n, dtype=DT)                          # Re, Im independent N(0, 1/2)
    X = alg.operator((G + G.mH) / math.sqrt(2 * n), is_self_adjoint=True)
    X2 = X @ X
    m4 = ((X2 @ X2).trace.real / n).numpy()
    m2 = (X2.trace.real / n).numpy()
    assert abs(zscore(m4, 2 + 1 / n ** 2)) < 5
    assert abs(zscore(m2, 1.0)) < 5


def test_proposition_1_10_pair_covariances():
    """Proposition 1.10, proof: E[H_ij H_kl] = delta_il delta_jk, E[G_ij G_kl] = 0, E[G_ij conj G_kl] = delta_ik delta_jl."""
    seed(36)
    n, B = 3, 100000
    G = torch.randn(B, n, n, dtype=DT)
    H = (G + G.mH) / math.sqrt(2)
    for (i, j, k, l) in [(0, 1, 1, 0), (0, 1, 0, 1), (0, 0, 0, 0), (0, 2, 1, 0), (1, 2, 2, 1), (2, 2, 1, 1)]:
        exact = float(i == l and j == k)
        v = (H[:, i, j] * H[:, k, l])
        assert abs(zscore(v.real.numpy(), exact)) < 5 and abs(zscore(v.imag.numpy(), 0.0)) < 5
    v = G[:, 0, 1] ** 2
    assert abs(zscore(v.real.numpy(), 0.0)) < 5 and abs(zscore(v.imag.numpy(), 0.0)) < 5
    assert abs(zscore((G[:, 0, 1].abs() ** 2).numpy(), 1.0)) < 5


def _gue_moment_by_pairings(p, n):
    """E Tr H^{2p} for E[H_ij H_kl] = delta_il delta_jk, by the Isserlis-Wick formula: sum over pairings pi of
    n^{#cycles of gamma o pi}, gamma the cyclic shift."""
    def pairings(idx):
        if not idx:
            yield []
            return
        for t in range(1, len(idx)):
            for rest in pairings(idx[1:t] + idx[t + 1:]):
                yield [(idx[0], idx[t])] + rest
    m = 2 * p
    total = 0
    for pi in pairings(list(range(m))):
        partner = {}
        for a, b in pi:
            partner[a], partner[b] = b, a
        perm = [(partner[a] + 1) % m for a in range(m)]              # gamma o pi
        seen, cycles = set(), 0
        for a in range(m):
            if a not in seen:
                cycles += 1
                while a not in seen:
                    seen.add(a)
                    a = perm[a]
        total += n ** cycles
    return total


@pytest.mark.parametrize('n', [1, 2, 3, 4, 7])
def test_proposition_1_10_exact_by_wick_enumeration(n):
    """Proposition 1.10 / Theorem A.24: the pairings give E Tr H^4 = 2n^3 + n (and E Tr H^6 = 5n^4 + 10n^2)."""
    assert _gue_moment_by_pairings(1, n) == n ** 2
    assert _gue_moment_by_pairings(2, n) == 2 * n ** 3 + n
    assert _gue_moment_by_pairings(3, n) == 5 * n ** 4 + 10 * n ** 2


# =====================================================================================================================
# Chapter 2: finite-dimensional von Neumann algebras
# =====================================================================================================================

def _commutant_dim(mats, tol=1e-9):
    """dim {X : XA = AX for all A in mats} by the null space of X -> (AX - XA) (row-major vec)."""
    d = mats[0].shape[0]
    eye = torch.eye(d, dtype=DT)
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in mats])
    s = torch.linalg.svdvals(L)
    return int((s <= tol * max(1.0, s[0].item())).sum()) + max(0, d * d - L.shape[0])


def _random_vn_algebra(ks, ms, n_gens=3, s=0):
    """Generators V((+)_c a_c (x) 1_{m_c})V^* of a von Neumann algebra with known (k_c, m_c), and V."""
    g = torch.Generator().manual_seed(s)
    n = sum(k * m for k, m in zip(ks, ms))
    V, _ = torch.linalg.qr(torch.randn(n, n, dtype=DT, generator=g))
    gens = []
    for _ in range(n_gens):
        blk = torch.block_diag(*[torch.kron(torch.randn(k, k, dtype=DT, generator=g), torch.eye(m, dtype=DT))
                                 for k, m in zip(ks, ms)])
        gens.append(V @ blk @ V.mH)
    return gens, V


STRUCTURES = [((1,), (3,)), ((3,), (1,)), ((2,), (2,)), ((1, 2), (2, 1)), ((2, 1, 3), (1, 2, 1)), ((1, 1), (1, 2))]


@pytest.mark.parametrize('ks,ms', STRUCTURES, ids=[f"k{k}m{m}" for k, m in STRUCTURES])
def test_theorem_2_2_and_2_6_dimensions_and_double_commutant(ks, ms):
    """Theorems 2.2 and 2.6: dim M = sum k_c^2, dim M' = sum m_c^2, M'' = M (wedderburn recovers (k_c, m_c))."""
    seed(200)
    gens, V = _random_vn_algebra(ks, ms, s=sum(ks) + 10 * sum(ms))
    comm = decompose.commutant_basis(gens)
    assert comm.shape[0] == sum(m * m for m in ms)
    assert _commutant_dim(gens) == sum(m * m for m in ms)
    double = decompose.commutant_basis(comm)
    assert double.shape[0] == sum(k * k for k in ks)                    # dim M'' = dim M
    # every element of M lies in M'' (S subset S''), so M'' = M by dimensions
    proj = torch.einsum('aij,bij->ab', double.conj(), torch.stack(gens))   # coordinates in an orthonormal basis of M''
    resid = torch.stack(gens) - torch.einsum('ab,aij->bij', proj, double)
    assert resid.abs().max() < 1e-9
    dec = decompose.wedderburn(gens)
    assert sorted(zip(dec.k_factors, dec.multiplicities)) == sorted(zip(ks, ms))
    assert dec.dims == (sum(k * k for k in ks), sum(m * m for m in ms))
    # V^* M V block form, and the minimal central projections have rank k_c m_c and sum to 1
    Pc = dec.central_projections()
    assert maxdiff(Pc.sum(0), torch.eye(dec.n, dtype=DT)) < 1e-9
    assert sorted(int(round(torch.trace(P).real.item())) for P in Pc) == sorted(k * m for k, m in zip(ks, ms))
    for P in Pc:
        assert maxdiff(P @ P, P) < 1e-9 and all(maxdiff(P @ g, g @ P) < 1e-9 for g in gens)


@pytest.mark.parametrize('ks,ms', STRUCTURES, ids=[f"k{k}m{m}" for k, m in STRUCTURES])
def test_section_2_3_trace_and_eigenvalue_multiplicities(ks, ms):
    """Section 2.3: Tr_{C^n}(A) = sum_c m_c Tr A_c, and every eigenvalue of A_c appears m_c times in A."""
    seed(41)
    alg = alg_of(list(ks))
    H = hermitian(alg)
    big = dense(H, list(ms))
    assert abs(torch.trace(big) - sum(m * torch.trace(X) for m, X in zip(ms, blocks(H)))) < 1e-12
    w_big = torch.linalg.eigvalsh(big)
    w_exp = torch.sort(torch.cat([w[0].repeat_interleave(m) for w, m in zip(H.eigenvalues(), ms)])).values
    assert maxdiff(w_big, w_exp) < 1e-12
    # the trace of C^n restricted to M is tr_w with w_c = m_c: a DensityMatrix-independent check through Tr_blunt
    W = alg.central([float(m) for m in ms])
    assert abs((W @ H).Tr_blunt()[0] - torch.trace(big)) < 1e-12


def test_section_2_1_triple_commutant():
    """Section 2.1: S' = S''' for a self-adjoint set S (dimensions of the commutants)."""
    seed(42)
    A = torch.randn(4, 4, dtype=DT)
    S = [A + A.mH, torch.diag(torch.tensor([1.0, 1.0, 2.0, 2.0], dtype=DT))]
    c1 = decompose.commutant_basis(S)
    c2 = decompose.commutant_basis(c1)
    c3 = decompose.commutant_basis(c2)
    assert c1.shape[0] == c3.shape[0]
    # S' subset S''': elements of S' commute with S''
    for X in c1:
        for Y in c2:
            assert maxdiff(X @ Y, Y @ X) < 1e-9


def test_lemma_2_1_invariant_subspace_projection_in_commutant():
    """Lemma 2.1: the projection onto an M-invariant subspace lies in M' (here the subspace M x)."""
    seed(201)
    gens, V = _random_vn_algebra((2, 1), (2, 1), s=5)
    x = torch.zeros(5, dtype=DT)
    x[0] = 1.0
    # K = M x: span of words applied to x; M is spanned by products of the generators (dim 5)
    basis = decompose.commutant_basis(decompose.commutant_basis(gens))
    K = torch.stack([b @ x for b in basis], dim=1)
    Q, R = torch.linalg.qr(K)
    r = int((R.diagonal().abs() > 1e-10).sum())
    U, s, _ = torch.linalg.svd(K)
    Q = U[:, : int((s > 1e-10).sum())]
    P = Q @ Q.mH
    for g in gens:
        assert maxdiff(g @ P, P @ g) < 1e-10
    assert Q.shape[1] == r


def test_lemma_2_3_commutative_algebra_minimal_projections():
    """Lemma 2.3: a commutative *-algebra is the span of C orthogonal projections (all k_c = 1); wedderburn finds them."""
    g = torch.Generator().manual_seed(7)
    n = 6
    V, _ = torch.linalg.qr(torch.randn(n, n, dtype=DT, generator=g))
    labels = torch.tensor([0, 0, 1, 2, 2, 2])                      # three joint eigenspaces, dims 2, 1, 3
    gens = [V @ torch.diag(torch.randn(3, dtype=DT, generator=g)[labels]) @ V.mH for _ in range(2)]
    dec = decompose.wedderburn(gens)
    assert dec.k_factors == (1, 1, 1) and sorted(dec.multiplicities) == [1, 2, 3]
    Pc = dec.central_projections()
    for c, P in enumerate(Pc):
        for d, Q in enumerate(Pc):
            assert maxdiff(P @ Q, P if c == d else torch.zeros_like(P)) < 1e-9


@pytest.mark.parametrize('k,m', [(1, 1), (2, 3), (3, 2), (4, 1)])
def test_lemma_2_5_factors(k, m):
    """Lemma 2.5: a factor is V(M_k (x) 1_m)V^*, with commutant V(1_k (x) M_m)V^*; its centre is C1."""
    seed(202)
    gens, V = _random_vn_algebra((k,), (m,), s=k * 10 + m)
    dec = decompose.wedderburn(gens)
    assert dec.k_factors == (k,) and dec.multiplicities == (m,)
    comm = decompose.commutant_basis(gens)
    assert comm.shape[0] == m * m
    # centre M cap M' = C 1
    alg_basis = decompose.commutant_basis(comm)
    centre = decompose.commutant_basis(torch.cat([alg_basis, comm]))
    assert centre.shape[0] == 1


def test_section_2_4_s3_example():
    """Section 2.4: S_3 on (C^2)^{(x)3}: dim pi' = 20, dim pi'' = 5, k = (1, 2), m = (4, 2), trivial on F (dim 4)."""
    seed(203)
    G = decompose.permutation_representation(3, 2)
    assert G.shape == (6, 8, 8)
    assert decompose.commutant_dimension(G) == 20
    assert _commutant_dim(list(G)) == 20
    span = torch.linalg.svdvals(G.reshape(6, 64))
    assert int((span > 1e-10).sum()) == 5
    dec = decompose.from_group(G)
    assert sorted(zip(dec.k_factors, dec.multiplicities)) == [(1, 4), (2, 2)]
    comm = decompose.commutant_basis(list(G))
    assert comm.shape[0] == 20 and decompose.commutant_basis(comm).shape[0] == 5      # Theorem 2.2: pi'' = span pi(g)
    # Step 1: the number of orbits of S_3 on pairs (x, y), i.e. multisets of size 3 from 4 elements
    pairs = itertools.product(itertools.product([0, 1], repeat=3), repeat=2)
    orbits = {tuple(sorted(zip(x, y))) for x, y in pairs}
    assert len(orbits) == 20
    # Step 3: F = fixed vectors has dimension 4, P_F = |G|^-1 sum pi(g) is its orthogonal projection
    PF = G.mean(0)
    assert maxdiff(PF @ PF, PF) < 1e-12 and maxdiff(PF.mH, PF) < 1e-12
    assert abs(torch.trace(PF).real - 4) < 1e-12
    # Step 4: rho_2 is the standard representation: character of pi = 4 chi_triv + 2 chi_std
    chi = torch.diagonal(G, dim1=-2, dim2=-1).sum(-1).real
    perms = list(itertools.permutations(range(3)))
    fixed = torch.tensor([sum(p[i] == i for i in range(3)) for p in perms], dtype=RDT)
    chi_std = fixed - 1                                               # permutation character minus trivial
    assert maxdiff(chi, 4 + 2 * chi_std) < 1e-12
    sign = torch.tensor([np.linalg.det(np.eye(3)[list(p)]) for p in perms], dtype=RDT)
    assert abs((chi * sign).mean()) < 1e-12                           # the sign representation occurs 0 times
    # Step 2: the transpositions (12), (23) do not commute
    assert maxdiff(G[1] @ G[2], G[2] @ G[1]) > 0.5


def test_section_2_4_s3_spectrum_of_random_symmetric_hamiltonian():
    """Section 2.4: K = E(Z + Z^*) has 4 simple eigenvalues on F and 2 doubly degenerate ones on F-perp (a.s.)."""
    seed(43)
    G = decompose.permutation_representation(3, 2)
    for _ in range(5):
        Z = torch.randn(8, 8, dtype=DT)
        K = decompose.twirl(G, Z + Z.mH)
        assert all(maxdiff(g @ K, K @ g) < 1e-12 for g in G)
        PF = G.mean(0)
        w, Q = torch.linalg.eigh(PF)
        on_F = torch.linalg.eigvalsh(Q[:, w > 0.5].mH @ K @ Q[:, w > 0.5])
        rest = torch.linalg.eigvalsh(Q[:, w < 0.5].mH @ K @ Q[:, w < 0.5])
        assert abs(rest[0] - rest[1]) < 1e-12 and abs(rest[2] - rest[3]) < 1e-12
        six = torch.cat([on_F, rest[::2]])
        gaps = (six[:, None] - six[None, :]).abs() + torch.eye(6, dtype=RDT)
        assert gaps.min() > 1e-6


def test_proposition_2_7_irreducible_inequivalent_blocks():
    """Proposition 2.7: V^* pi(g) V = (+)_c rho_c(g) (x) 1_{m_c}; the rho_c are irreducible and inequivalent."""
    seed(204)
    G = decompose.permutation_representation(3, 2)
    dec = decompose.from_group(G)
    A = dec.to_algebra(G)                                            # blocks rho_c(g)
    assert maxdiff(dec.from_algebra(A), G) < 1e-10
    for c, k in enumerate(dec.k_factors):
        rho = A.matrix[:, c, :k, :k]
        # unitary representation
        assert maxdiff(rho.mH @ rho, torch.eye(k, dtype=DT).expand(6, k, k)) < 1e-10
        # irreducible: |G|^-1 sum |chi|^2 = 1
        chi = torch.diagonal(rho, dim1=-2, dim2=-1).sum(-1)
        assert abs((chi.abs() ** 2).mean() - 1) < 1e-10
        # the rho_c(g) span M_k
        assert int((torch.linalg.svdvals(rho.reshape(6, -1)) > 1e-10).sum()) == k * k
    chis = [torch.diagonal(A.matrix[:, c, :k, :k], dim1=-2, dim2=-1).sum(-1) for c, k in enumerate(dec.k_factors)]
    assert abs((chis[0] * chis[1].conj()).mean()) < 1e-10           # inequivalent: orthogonal characters


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_proposition_2_8_conditional_expectation_onto_centre(ks):
    """Proposition 2.8: E(A) = sum_c Tr A_c / k_c P_c is positive, unital, the identity on Z(M), a Z(M)-bimodule map, trace preserving."""
    seed(44)
    alg = alg_of(ks)
    E = channels.center_expectation(alg)
    A = ginibre(alg, 3)
    z1 = torch.randn(len(ks), dtype=DT)
    z2 = torch.randn(len(ks), dtype=DT)
    Z1 = alg.central([complex(x) for x in z1])
    Z2 = alg.central([complex(x) for x in z2])
    assert maxdiff(E(Z1), Z1) < 1e-12
    assert maxdiff(E(Z1 @ A @ Z2), Z1 @ E(A) @ Z2) < 1e-12
    assert maxdiff(E(E(A)), E(A)) < 1e-12
    assert maxdiff(E(alg.identity()), alg.identity()) < 1e-12
    P = positive(alg, 3)
    assert (alg.operator(E(P).matrix, is_self_adjoint=True).lambda_min >= -1e-12).all()
    w = torch.rand(len(ks), dtype=RDT) + 0.1                          # every trace tr_w
    W = alg.central([float(x) for x in w])
    assert maxdiff((W @ E(A)).Tr_blunt(), (W @ A).Tr_blunt()) < 1e-12
    for tr in ('Tr_blunt', 'Tr_norm', 'tau_vN'):
        assert maxdiff(getattr(E(A), tr)(), getattr(A, tr)()) < 1e-12
    # E(A) is central: it commutes with every element
    B = ginibre(alg, 3)
    assert maxdiff(E(A) @ B, B @ E(A)) < 1e-12
    # omega(E(A)) = sum_c omega(P_c) Tr(A_c) / k_c for a state omega
    rho = states.random_density_matrix(alg, 3)
    p = rho.sector_probabilities()
    trA = torch.stack([A.matrix[:, c, :k, :k].diagonal(dim1=-2, dim2=-1).sum(-1) / k for c, k in enumerate(ks)], -1)
    assert maxdiff(rho.expectation(E(A)), (p * trA).sum(-1)) < 1e-12


@pytest.mark.parametrize('ks', [[2, 3], [1, 4, 2], [3, 3]], ids=['2-3', '1-4-2', '3-3'])
def test_section_2_6_block_diagonal_expectation_on_MD(ks):
    """Section 2.6: E_M(X) = sum_c P_c X P_c is a Tr-preserving conditional expectation of M_D onto (+)_c M_{k_c}."""
    seed(45)
    D = sum(ks)
    MD = alg_of([D])
    Ps = []
    off = 0
    for k in ks:
        P = torch.zeros(1, 1, D, D, dtype=DT)
        P[0, 0, off:off + k, off:off + k] = torch.eye(k, dtype=DT)
        Ps.append(MD.operator(P, is_self_adjoint=True, is_projection=True, is_positive=True))
        off += k
    EM = channels.measurement_channel(Ps)
    X = ginibre(MD, 2)
    EX = EM(X)
    mask = torch.block_diag(*[torch.ones(k, k) for k in ks]).bool()
    assert EX.matrix[:, 0][:, ~mask].abs().max() == 0                     # lands in the block-diagonal algebra
    assert maxdiff(EX.matrix[:, 0][:, mask], X.matrix[:, 0][:, mask]) == 0
    assert maxdiff(EM(EX), EX) < 1e-14                                   # identity on M
    B1 = MD.operator(EM(ginibre(MD, 2)).matrix)
    B2 = MD.operator(EM(ginibre(MD, 2)).matrix)
    assert maxdiff(EM(B1 @ X @ B2), B1 @ EX @ B2) < 1e-12                # M-bimodule
    assert maxdiff(EX.trace, X.trace) < 1e-12                            # Tr preserving
    P = positive(MD, 2)
    assert (MD.operator(EM(P).matrix, is_self_adjoint=True).lambda_min >= -1e-12).all()


def test_section_2_6_conditional_expectations_not_unique():
    """Section 2.6: on M_2, A -> A_11 1 and A -> A_22 1 are two different conditional expectations onto C1."""
    seed(46)
    A = torch.randn(2, 2, dtype=DT)
    E1, E2 = A[0, 0] * torch.eye(2, dtype=DT), A[1, 1] * torch.eye(2, dtype=DT)
    assert maxdiff(E1, E2) > 1e-3
    # positivity and bimodule property over C1 hold trivially; unitality: E(1) = 1 for both; both differ from Tr/2
    alg = alg_of([2])
    E = channels.center_expectation(alg)
    EA = E(alg.operator(A.reshape(1, 1, 2, 2))).matrix[0, 0]
    assert maxdiff(EA, (A[0, 0] + A[1, 1]) / 2 * torch.eye(2, dtype=DT)) < 1e-14


def test_exercise_2_1_centre_of_matrix_algebra():
    """Chapter 2, Exercise 1: Z(M_k) = C1, and Z((+)_c M_{k_c}) has dimension C."""
    seed(205)
    for k in (1, 2, 4):
        units = [torch.outer(torch.eye(k, dtype=DT)[i], torch.eye(k, dtype=DT)[j]) for i in range(k) for j in range(k)]
        assert _commutant_dim(units) == 1
    ks = [2, 1, 3]
    D = sum(ks)
    units = []
    off = 0
    for k in ks:
        for i in range(k):
            for j in range(k):
                E = torch.zeros(D, D, dtype=DT)
                E[off + i, off + j] = 1
                units.append(E)
        off += k
    # centre = M cap M': elements of span(units) commuting with all units
    comm = decompose.commutant_basis(units)
    assert comm.shape[0] == len(ks)                                  # M' = (+)_c C 1_{k_c} here (m_c = 1)


def test_exercise_2_2_cyclic_shift():
    """Chapter 2, Exercise 2: {S}' = {S}'' = span{1, S, S^2} (circulants), dimension 3."""
    seed(206)
    S = torch.roll(torch.eye(3, dtype=DT), 1, dims=0)
    c1 = decompose.commutant_basis([S])
    assert c1.shape[0] == 3
    assert decompose.commutant_basis(c1).shape[0] == 3
    # S u_j = omega^j u_j for u_j = sum_i omega^{-ij} e_i
    w = np.exp(2j * np.pi / 3)
    for j in range(3):
        u = torch.tensor([w ** (-i * j) for i in range(3)], dtype=DT)
        assert maxdiff(S @ u, w ** j * u) < 1e-12
    # every circulant commutes with S, and the commutant basis lies in span{1, S, S^2}
    span = torch.stack([torch.eye(3, dtype=DT), S, S @ S]).reshape(3, 9)
    coeff = torch.linalg.lstsq(span.T, c1.reshape(3, 9).T, driver='gelsd').solution
    assert maxdiff(span.T @ coeff, c1.reshape(3, 9).T) < 1e-12


def test_exercise_2_3_twirl_is_conditional_expectation():
    """Chapter 2, Exercise 3: the group average is a Tr-preserving conditional expectation of M_n onto pi(G)'."""
    seed(47)
    G = decompose.permutation_representation(3, 2)
    X = torch.randn(3, 8, 8, dtype=DT)
    EX = decompose.twirl(G, X)
    assert all(maxdiff(g @ EX, EX @ g) < 1e-12 for g in G)
    assert maxdiff(decompose.twirl(G, EX), EX) < 1e-12
    assert maxdiff(torch.diagonal(EX, dim1=-2, dim2=-1).sum(-1), torch.diagonal(X, dim1=-2, dim2=-1).sum(-1)) < 1e-12
    B1 = decompose.twirl(G, torch.randn(8, 8, dtype=DT))
    B2 = decompose.twirl(G, torch.randn(8, 8, dtype=DT))
    assert maxdiff(decompose.twirl(G, B1 @ X @ B2), B1 @ EX @ B2) < 1e-12
    P = X @ X.mH
    assert torch.linalg.eigvalsh(decompose.twirl(G, P)).min() > -1e-12
    # the Hilbert-Schmidt orthogonal projection: <EX, Y - EY> = 0
    Y = torch.randn(8, 8, dtype=DT)
    assert abs(torch.vdot(EX[0].reshape(-1), (Y - decompose.twirl(G, Y)).reshape(-1))) < 1e-12


# =====================================================================================================================
# Chapter 3: traces, states and entropy
# =====================================================================================================================

@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_theorem_3_1_traces_classification(ks):
    """Theorem 3.1: tr_w(E^c_ij) = delta_ij w_c; the library traces have w = 1, 1/k_c, 1/(C k_c); trace property; faithful."""
    seed(51)
    alg = alg_of(ks)
    C, D = len(ks), sum(ks)
    weights = {'Tr_blunt': [1.0] * C, 'Tr_norm': [1 / k for k in ks], 'tau_vN': [1 / (C * k) for k in ks]}
    for name, w in weights.items():
        for c, k in enumerate(ks):
            for i, j in [(0, 0), (k - 1, k - 1), (0, k - 1)]:
                blk = [torch.zeros(kk, kk, dtype=DT) for kk in ks]
                blk[c][i, j] = 1
                val = getattr(alg.from_blocks(blk), name)()[0]
                assert abs(val - (w[c] if i == j else 0)) < 1e-15
    one = alg.identity()
    assert abs(one.Tr_blunt()[0] - D) < 1e-12 and abs(one.Tr_norm()[0] - C) < 1e-12 and abs(one.tau_vN()[0] - 1) < 1e-12
    A, B = ginibre(alg, 4), ginibre(alg, 4)
    for name in weights:
        assert maxdiff(getattr(A @ B, name)(), getattr(B @ A, name)()) < 1e-11
        # positivity and faithfulness: tr_w(A^*A) = sum_c w_c ||A_c||_2^2 > 0
        val = getattr(A.adjoint() @ A, name)()
        expect = sum(weights[name][c] * (A.matrix[:, c, :k, :k].abs() ** 2).sum((-2, -1)) for c, k in enumerate(ks))
        assert maxdiff(val.real, expect) < 1e-11 and (val.real > 0).all() and val.imag.abs().max() < 1e-12
    # tau_vN is the tracial state with sector weights 1/C
    for c in range(C):
        assert abs(central_projection(alg, c).tau_vN()[0] - 1 / C) < 1e-15


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_theorem_3_2_densities_for_all_traces(ks):
    """Theorem 3.2: omega(A) = tr_w(rho_w A) for each faithful trace, rho_w >= 0, tr_w(rho_w) = 1, w_c rho_w = w'_c rho_w'."""
    seed(52)
    alg = alg_of(ks)
    C = len(ks)
    rho = states.random_density_matrix(alg, 3)
    A = ginibre(alg, 3)
    omega = rho.expectation(A)
    w = {'blunt': [1.0] * C, 'norm': [1 / k for k in ks], 'tau_vN': [1 / (C * k) for k in ks]}
    tr = {'blunt': 'Tr_blunt', 'norm': 'Tr_norm', 'tau_vN': 'tau_vN'}
    dens = {}
    for name in w:
        d = rho.density(name)
        dens[name] = d
        dop = alg.operator(d, is_self_adjoint=True)
        assert maxdiff(getattr(dop @ A, tr[name])(), omega) < 1e-12
        assert maxdiff(getattr(dop, tr[name])(), torch.ones(3, dtype=DT)) < 1e-12
        assert (dop.lambda_min >= -1e-14).all()
        back = DensityMatrix.from_density(alg, d, trace=name)
        assert maxdiff(back, rho) < 1e-14
    for a, b in itertools.combinations(w, 2):
        for c, k in enumerate(ks):
            assert maxdiff(w[a][c] * dens[a][:, c, :k, :k], w[b][c] * dens[b][:, c, :k, :k]) < 1e-14
    # uniqueness: the density of a linear functional (omega(A) = Tr(X A)) is recovered by the matrix units
    X = ginibre(alg)
    vals = []
    for c, k in enumerate(ks):
        blkd = torch.zeros(k, k, dtype=DT)
        for i in range(k):
            for j in range(k):
                blk = [torch.zeros(kk, kk, dtype=DT) for kk in ks]
                blk[c][j, i] = 1                                         # E_ji: Tr(X E_ji) = X_ij
                blkd[i, j] = (X @ alg.from_blocks(blk)).Tr_blunt()[0]
        vals.append(blkd)
    assert all(maxdiff(v, X.matrix[0, c, :k, :k]) < 1e-14 for c, (v, k) in enumerate(zip(vals, ks)))


def test_theorem_3_2_positive_functional_iff_positive_density():
    """Theorem 3.2: omega positive iff rho_w >= 0: a density with a negative eigenvalue gives omega(|x><x|) < 0."""
    seed(53)
    alg = alg_of([3, 2])
    d = alg.operator_from_eigenvalues(lambda k: torch.tensor([-0.1] + [0.5] * (k - 1)), force_self_adjoint=True)
    with pytest.raises(ValueError):
        DensityMatrix(alg, matrix=(d * (1 / d.Tr_blunt().real.item())).matrix)
    # the vector state of the eigenvector with negative eigenvalue
    w, V = torch.linalg.eigh(d.matrix[0, 0, :3, :3])
    x = V[:, 0]
    proj = alg.from_blocks([torch.outer(x, x.conj()), torch.zeros(2, 2, dtype=DT)])
    assert (d @ proj).Tr_blunt().real.item() < 0


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_section_3_3_classical_part(ks):
    """Section 3.3: p_c = omega(P_c) = Tr rho_c is a distribution; rho = (+)_c p_c sigma_c; the sector measurement does not disturb rho."""
    seed(54)
    alg = alg_of(ks)
    rho = states.random_density_matrix(alg, 4)
    p = rho.sector_probabilities()
    assert maxdiff(p.sum(-1), torch.ones(4, dtype=RDT)) < 1e-12 and (p >= 0).all()
    for c in range(len(ks)):
        assert maxdiff(rho.expectation(central_projection(alg, c)).real, p[:, c]) < 1e-14
    P = [central_projection(alg, c) for c in range(len(ks))]
    assert maxdiff(channels.measurement_channel(P)(rho), rho) < 1e-14


def _entropy_np(w):
    w = w[w > 1e-15]
    return float(-(w * np.log(w)).sum())


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_theorem_3_3_entropy_of_mixture_of_sectors(ks):
    """Theorem 3.3: S(rho) = H(p) + sum_c p_c S(sigma_c), also with an empty sector (p_c = 0) and pure sector states."""
    seed(55)
    alg = alg_of(ks)
    for case in ('random', 'empty_sector', 'pure_sectors'):
        if case == 'random':
            rho = states.random_density_matrix(alg, 2)
        elif case == 'empty_sector':
            if len(ks) == 1:
                continue
            G = [torch.randn(2, k, k, dtype=DT) for k in ks]
            G[0] = torch.zeros_like(G[0])
            mats = [g @ g.mH for g in G]
            m = alg.from_blocks(mats).matrix
            rho = DensityMatrix(alg, matrix=m / m.diagonal(dim1=-2, dim2=-1).sum((-2, -1))[:, None, None, None])
        else:
            vs = [torch.randn(2, k, 1, dtype=DT) for k in ks]
            mats = [v @ v.mH for v in vs]
            m = alg.from_blocks(mats).matrix
            rho = DensityMatrix(alg, matrix=m / m.diagonal(dim1=-2, dim2=-1).sum((-2, -1))[:, None, None, None])
        S = rho.von_neumann_entropy()
        p = rho.sector_probabilities()
        for b in range(2):
            H = _entropy_np(p[b].numpy())
            Ssec = 0.0
            for c, k in enumerate(ks):
                if p[b, c] > 0:
                    sig = rho.matrix[b, c, :k, :k] / p[b, c]
                    Ssec += p[b, c].item() * _entropy_np(torch.linalg.eigvalsh(sig).numpy())
            assert abs(S[b].item() - (H + Ssec)) < 1e-10
            if case == 'pure_sectors':
                assert abs(S[b].item() - H) < 1e-10                    # S(sigma_c) = 0


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_section_3_4_entropy_relative_to_trace(ks):
    """Section 3.4: S_w(rho) = -tr_w(rho_w log rho_w) = S(rho) + sum_c p_c log w_c (Tr_norm and tau_vN)."""
    seed(56)
    alg = alg_of(ks)
    C = len(ks)
    rho = states.random_density_matrix(alg, 3)
    S, p = rho.von_neumann_entropy(), rho.sector_probabilities()
    for name, tr, w in (('norm', 'Tr_norm', [1 / k for k in ks]), ('tau_vN', 'tau_vN', [1 / (C * k) for k in ks]),
                        ('blunt', 'Tr_blunt', [1.0] * C)):
        d = alg.operator(rho.density(name), is_self_adjoint=True, is_positive=True)
        Sw = -getattr(d @ d.log(), tr)().real
        assert maxdiff(Sw, S + (p * torch.log(torch.tensor(w, dtype=RDT))).sum(-1)) < 1e-10


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_section_3_5_relative_entropy_independent_of_trace(ks):
    """Section 3.5: tr_w(rho_w log rho_w) - tr_w(rho_w log sigma_w) does not depend on the faithful trace."""
    seed(57)
    alg = alg_of(ks)
    rho, sig = states.random_density_matrix(alg, 3), states.random_density_matrix(alg, 3)
    D0 = states.relative_entropy(rho, sig)
    for name, tr in (('blunt', 'Tr_blunt'), ('norm', 'Tr_norm'), ('tau_vN', 'tau_vN')):
        r = alg.operator(rho.density(name), is_self_adjoint=True, is_positive=True)
        s = alg.operator(sig.density(name), is_self_adjoint=True, is_positive=True)
        Dw = getattr(r @ r.log(), tr)().real - getattr(r @ s.log(), tr)().real
        assert maxdiff(Dw, D0) < 1e-10


def test_lemma_3_4_tangent_inequality_and_strict_convexity():
    """Lemma 3.4: f(a) >= f(b) + (a - b) f'(b) for convex f, strictly for a != b if f'' > 0 (incl. one-sided at 0)."""
    rng = seed(58)
    cases = [(lambda t: t * np.log(t, where=t > 0, out=np.zeros_like(t)), lambda t: np.log(t) + 1, (0.0, 5.0)),
             (np.exp, np.exp, (-3.0, 3.0)), (lambda t: t ** 4, lambda t: 4 * t ** 3, (-2.0, 2.0)),
             (lambda t: -np.log(t), lambda t: -1 / t, (0.01, 4.0))]
    for f, fp, (lo, hi) in cases:
        a = rng.uniform(lo, hi, 20000)
        b = rng.uniform(max(lo, 1e-3), hi, 20000)
        gap = f(a) - f(b) - (a - b) * fp(b)
        assert (gap >= -1e-12).all() and (gap[np.abs(a - b) > 1e-3] > 0).all()
    # f(t) = t log t at the endpoint a = 0: f(0) - f(b) + b f'(b) = b > 0
    b = rng.uniform(0.01, 3, 100)
    assert np.allclose(0 - b * np.log(b) + b * (np.log(b) + 1), b)


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_theorem_3_5_klein_inequality(ks):
    """Theorem 3.5: Tr[f(X) - f(Y) - (X - Y) f'(Y)] >= 0 for convex f; equality for X = Y (f strictly convex)."""
    seed(59)
    alg = alg_of(ks)
    X, Y = hermitian(alg, 6), hermitian(alg, 6)
    for f, fp in ((torch.exp, torch.exp), (lambda w: w ** 2, lambda w: 2 * w), (lambda w: w ** 4, lambda w: 4 * w ** 3),
                  (lambda w: torch.abs(w) ** 3, lambda w: 3 * w * torch.abs(w))):
        val = (X.apply_function(f) - Y.apply_function(f) - (X - Y) @ Y.apply_function(fp)).Tr_blunt()
        assert (val.real >= -1e-10).all() and val.imag.abs().max() < 1e-10
        assert (val.real > 1e-8).all()                                   # X != Y: strict (strictly convex f)
        val0 = (X.apply_function(f) - X.apply_function(f) - (X - X) @ X.apply_function(fp)).Tr_blunt()
        assert val0.abs().max() < 1e-12
    # f(t) = t^2: Tr[X^2 - Y^2 - 2(X - Y)Y] = ||X - Y||_2^2 exactly
    val = (X @ X - Y @ Y - (X - Y) @ Y * 2.0).Tr_blunt().real
    assert maxdiff(val, (X - Y).frobenius_norm() ** 2) < 1e-10


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_corollary_3_6_relative_entropy_positive_and_max_entropy(ks):
    """Corollary 3.6: S(rho||sigma) >= 0 with equality iff rho = sigma; S(rho) <= log D, with equality iff rho = 1/D."""
    seed(60)
    alg = alg_of(ks)
    D = sum(ks)
    rho, sig = states.random_density_matrix(alg, 8), states.random_density_matrix(alg, 8)
    assert (states.relative_entropy(rho, sig) > 0).all() or D == 1      # D = 1: the only state is 1
    assert states.relative_entropy(rho, rho).abs().max() < 1e-10
    assert (rho.von_neumann_entropy() <= math.log(D) + 1e-12).all()
    mm = states.maximally_mixed_state(alg)
    assert abs(mm.von_neumann_entropy().item() - math.log(D)) < 1e-12
    # S(rho||1/D) = log D - S(rho)
    mm8 = states.maximally_mixed_state(alg, 8)
    assert maxdiff(states.relative_entropy(rho, mm8), math.log(D) - rho.von_neumann_entropy()) < 1e-10
    # pure states have entropy 0 (rank-one random states)
    pure = states.random_density_matrix(alg, 3, rank=1)
    if len(ks) == 1:
        assert pure.von_neumann_entropy().abs().max() < 1e-10
    # near-maximal states: entropy strictly below log D
    eps_state = mm8.mix(rho, 1e-3)
    assert (eps_state.von_neumann_entropy() < math.log(D)).all() or D == 1


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_theorem_3_7_chain_rule_and_monotonicity(ks):
    """Theorem 3.7: S(rho||sigma) = sum p log(p/q) + sum p_c S(rho_c^ || sigma_c^) >= S(E rho || E sigma) = KL(p||q)."""
    seed(61)
    alg = alg_of(ks)
    rho, sig = states.random_density_matrix(alg, 4), states.random_density_matrix(alg, 4)
    p, q = rho.sector_probabilities(), sig.sector_probabilities()
    Drel = states.relative_entropy(rho, sig)
    kl = (p * torch.log(p / q)).sum(-1)
    E = channels.center_expectation(alg)
    assert maxdiff(states.relative_entropy(E(rho), E(sig)), kl) < 1e-10
    assert (Drel >= kl - 1e-12).all()
    # the chain rule as an equality, with the conditional states in M_{k_c}
    cond = torch.zeros(4, dtype=RDT)
    for c, k in enumerate(ks):
        a = alg_of([k])
        rc = DensityMatrix(a, matrix=(rho.matrix[:, c, :k, :k] / p[:, c, None, None]).unsqueeze(1))
        sc = DensityMatrix(a, matrix=(sig.matrix[:, c, :k, :k] / q[:, c, None, None]).unsqueeze(1))
        cond = cond + p[:, c] * states.relative_entropy(rc, sc)
    assert maxdiff(Drel, kl + cond) < 1e-10
    # E rho is the density of omega o E: Tr(rho E(A)) = Tr(E(rho) A)
    A = ginibre(alg, 4)
    assert maxdiff(rho.expectation(E(A)), E(rho).Tr_blunt() * 0 + (E(rho) @ A).Tr_blunt()) < 1e-12
    # E rho = (+)_c (p_c / k_c) 1
    for c, k in enumerate(ks):
        assert maxdiff(E(rho).matrix[:, c, :k, :k], (p[:, c] / k)[:, None, None] * torch.eye(k, dtype=DT)) < 1e-14


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_proposition_3_8_michelson_contrast_monotone(ks):
    """Proposition 3.8: Delta(Phi(A)) <= Delta(A) for positive unital Phi (mixed unitary, E, dephasing, depolarising)."""
    seed(62)
    alg = alg_of(ks)
    A = alg.operator_from_eigenvalues(lambda d: 0.05 + torch.rand(6, d), batch_size=6, force_positive=True)
    d0 = A.michelson_contrast
    assert ((d0 >= 0) & (d0 <= 1)).all()
    maps = [channels.random_mixed_unitary_channel(alg, 3, batch_size=6), channels.center_expectation(alg),
            channels.dephasing_channel(alg, 0.4), channels.depolarizing_channel(alg, 0.7),
            channels.unitary_channel(alg.random_unitary_operator(6))]
    for Phi in maps:
        assert maxdiff(Phi(alg.identity(6)), alg.identity(6)) < 1e-12              # unital
        out = alg.operator(Phi(A).matrix, is_self_adjoint=True)
        out._is_positive = True
        assert (out.michelson_contrast <= d0 + 1e-12).all()
    # unitary conjugation preserves the spectrum, hence Delta
    U = channels.unitary_channel(alg.random_unitary_operator(6))
    out = alg.operator(U(A).matrix, is_self_adjoint=True, is_positive=True)
    assert maxdiff(out.michelson_contrast, d0) < 1e-12


def test_section_3_6_michelson_contrast_exact_values():
    """Section 3.6: Delta = (l_max - l_min)/(l_max + l_min) over all sectors; Delta = 0 iff A = lambda 1; Delta(0) = 0; Delta(E A) exact."""
    alg = alg_of([2, 3])
    A = alg.from_blocks([torch.diag(torch.tensor([1.0, 3.0], dtype=DT)), torch.diag(torch.tensor([2.0, 0.5, 2.5], dtype=DT))])
    A = alg.operator(A.matrix, is_self_adjoint=True, is_positive=True)
    assert abs(A.michelson_contrast.item() - (3 - 0.5) / (3 + 0.5)) < 1e-15
    assert alg.identity().michelson_contrast.item() == 0
    assert (alg.identity() * 7.0).michelson_contrast.item() == 0
    assert alg.zero().michelson_contrast.item() == 0
    # scale invariance and the conditional expectation onto the centre: blocks 2 and 5/3
    assert abs((A * 1e-8).michelson_contrast.item() - A.michelson_contrast.item()) < 1e-12
    EA = alg.operator(channels.center_expectation(alg)(A).matrix, is_self_adjoint=True, is_positive=True)
    assert abs(EA.michelson_contrast.item() - (2 - 5 / 3) / (2 + 5 / 3)) < 1e-14
    # a positive A that is not a multiple of 1 but has equal extremes in every sector except across sectors
    B = alg.central([1.0, 2.0])
    assert abs(B.michelson_contrast.item() - 1 / 3) < 1e-15
    # rank-deficient positive A: Delta = 1
    Pp = alg.operator_from_eigenvalues(lambda d: torch.tensor([0.0] + [1.0] * (d - 1)), force_positive=True)
    assert abs(Pp.michelson_contrast.item() - 1) < 1e-15


@pytest.mark.parametrize('k', [2, 3, 6])
def test_proposition_3_9_false_inequality(k):
    """Proposition 3.9: max_U |Tr(XUY)| = ||YX||_1 >= Tr(XY), attained at W^* (YX = W|YX|), equality iff XY = YX."""
    seed(63 + k)
    alg = alg_of([k])
    X = alg.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
    Y = alg.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
    YX = Y @ X
    W, _ = YX.polar()
    sup = (X @ W.adjoint() @ Y).trace.abs()
    t1 = YX.trace_norm()
    tXY = (X @ Y).trace
    assert abs(tXY.imag).max() < 1e-12 and (tXY.real > 0).all()
    assert maxdiff(sup, t1) < 1e-12
    assert (t1 - tXY.real > 1e-8).all()                                  # non-commuting: strict
    U = alg.random_unitary_operator(5000)
    Xr = alg.operator(X.matrix.expand(5000, -1, -1, -1))
    Yr = alg.operator(Y.matrix.expand(5000, -1, -1, -1))
    vals = (Xr @ U @ Yr).trace.abs()
    assert (vals <= t1 + 1e-10).all()
    # Proposition 11.11: E|Tr(XUY)|^2 = ||YX||_2^2 / k
    assert abs(zscore((vals ** 2).numpy(), (YX.frobenius_norm() ** 2 / k).item())) < 5
    # commuting pair: functions of the same operator, equality
    Yc = X.apply_function(lambda w: 1 + w ** 2)
    Yc._is_positive = True
    assert abs((Yc @ X).trace_norm() - (X @ Yc).trace.real).max() < 1e-12


def test_corollary_11_12_probability_bound_is_vacuous_in_tutorial_regime():
    """Section 3.7: the bound b^4/(a^4 k^2) of Corollary 11.12 is >= 1 for b/a = 6, k = 6 (says nothing)."""
    a, b, k = 0.2, 1.2, 6
    assert b ** 4 / (a ** 4 * k ** 2) >= 1
    assert k <= b ** 2 / a ** 2                                         # informative only for k > b^2 / a^2


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_exercise_3_1_tracial_state(ks):
    """Chapter 3, Exercise 1: the density of tau_vN w.r.t. tau_vN is 1; S = log C + (1/C) sum_c log k_c."""
    alg = alg_of(ks)
    C = len(ks)
    tau = states.tracial_state(alg)
    assert maxdiff(tau.density('tau_vN'), alg.identity().matrix) < 1e-14
    for c, k in enumerate(ks):
        assert maxdiff(tau.matrix[0, c, :k, :k], torch.eye(k, dtype=DT) / (C * k)) < 1e-15
    assert abs(tau.von_neumann_entropy().item() - (math.log(C) + sum(math.log(k) for k in ks) / C)) < 1e-12
    # tau_vN(A) = Tr(rho_tau A)
    seed(64)
    A = ginibre(alg, 2)
    assert maxdiff(tau.expectation(A), A.tau_vN()) < 1e-14


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_exercise_3_2_maximal_entropy_relative_to_trace(ks):
    """Chapter 3, Exercise 2: S_w(rho) <= log tr_w(1), with equality iff rho_w = 1/tr_w(1) (Tr_norm, tau_vN, blunt)."""
    seed(65)
    alg = alg_of(ks)
    C, D = len(ks), sum(ks)
    rho = states.random_density_matrix(alg, 6)
    S, p = rho.von_neumann_entropy(), rho.sector_probabilities()
    for w, trw1, maximiser in (([1.0] * C, D, states.maximally_mixed_state(alg)),
                               ([1 / k for k in ks], C, states.tracial_state(alg)),
                               ([1 / (C * k) for k in ks], 1.0, states.tracial_state(alg))):
        lw = torch.log(torch.tensor(w, dtype=RDT))
        Sw = S + (p * lw).sum(-1)
        assert (Sw <= math.log(trw1) + 1e-12).all()
        Sm = maximiser.von_neumann_entropy() + (maximiser.sector_probabilities() * lw).sum(-1)
        assert abs(Sm.item() - math.log(trw1)) < 1e-12
    # the finite Jensen inequality, with equality only for equal x_c
    rng = np.random.default_rng(65)
    pp = rng.dirichlet(np.ones(5), 1000)
    x = rng.uniform(0.1, 3, (1000, 5))
    assert ((pp * np.log(x)).sum(1) <= np.log((pp * x).sum(1)) + 1e-14).all()
    assert np.allclose((pp * np.log(np.full_like(x, 2.0))).sum(1), np.log(2.0))


def test_exercise_3_3_measurement_on_MD_disturbs_superposition():
    """Chapter 3, Exercise 3: on M_D the sector measurement maps a superposition (S = 0) to a mixture (S = log 2)."""
    seed(66)
    ks = [2, 3]
    D = sum(ks)
    MD = alg_of([D])
    psi1 = torch.randn(2, dtype=DT)
    psi2 = torch.randn(3, dtype=DT)
    psi = torch.cat([psi1 / psi1.norm(), psi2 / psi2.norm()]) / math.sqrt(2)
    rho = DensityMatrix(MD, matrix=torch.outer(psi, psi.conj()).reshape(1, 1, D, D))
    assert abs(rho.von_neumann_entropy().item()) < 1e-10
    Ps = []
    for off, k in ((0, 2), (2, 3)):
        P = torch.zeros(1, 1, D, D, dtype=DT)
        P[0, 0, off:off + k, off:off + k] = torch.eye(k, dtype=DT)
        Ps.append(MD.operator(P, is_self_adjoint=True, is_projection=True))
    after = channels.measurement_channel(Ps)(rho)
    assert abs(after.von_neumann_entropy().item() - math.log(2)) < 1e-10
    # restricted to M = M_2 (+) M_3 the vector state is the mixture: same expectations of block-diagonal A
    alg = alg_of(ks)
    A = ginibre(alg)
    Ad = MD.operator(dense(A).reshape(1, 1, D, D))
    assert abs(rho.expectation(Ad) - after.expectation(Ad)).max() < 1e-12


# =====================================================================================================================
# Appendix A: facts from analysis and probability
# =====================================================================================================================

def test_lemma_A_3_unitary_conjugation_is_orthogonal_and_preserves_volume():
    """Lemma A.3 / conventions: X -> V X W (V, W unitary) preserves Re Tr(X^*Y); on Herm_n in the coordinates (h_ii, Re h_ij, Im h_ij) it has |det| = 1."""
    seed(70)
    alg = alg_of([3])
    V, W = alg.random_unitary(3), alg.random_unitary(3)
    X, Y = torch.randn(3, 3, dtype=DT), torch.randn(3, 3, dtype=DT)
    ip = lambda a, b: torch.trace(a.mH @ b).real
    assert abs(ip(V @ X @ W, V @ Y @ W) - ip(X, Y)) < 1e-12
    n = 3
    def coords(H):
        out = [H[i, i].real for i in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                out += [H[i, j].real, H[i, j].imag]
        return torch.stack(out)
    basis = []
    for idx in range(n * n):
        e = torch.zeros(n * n, dtype=RDT)
        e[idx] = 1
        H = torch.zeros(n, n, dtype=DT)
        t = 0
        for i in range(n):
            H[i, i] = e[t]
            t += 1
        for i in range(n):
            for j in range(i + 1, n):
                H[i, j] = e[t] + 1j * e[t + 1]
                H[j, i] = e[t] - 1j * e[t + 1]
                t += 2
        basis.append(coords(V @ H @ V.mH))
    L = torch.stack(basis, dim=1)
    assert abs(abs(torch.linalg.det(L).item()) - 1) < 1e-12


def test_lemma_A_4_riemann_sums_and_cesaro_mean():
    """Lemma A.4 / Proposition 8.12: the Cesaro mean of a unitary orbit is a state, and ||T_t s - s||_1 <= 2t/tau."""
    seed(71)
    alg = alg_of([2, 3])
    H = hermitian(alg)
    rho = states.random_density_matrix(alg)
    tau, N = 3.0, 600
    def orbit(t):
        U = (H * (-1j * t)).expm()
        return (U @ rho @ U.adjoint()).matrix
    ts = (torch.arange(N, dtype=RDT) + 0.5) * tau / N
    sig = sum(orbit(t.item()) for t in ts) / N                         # midpoint Riemann sum / tau
    sop = alg.operator(sig, is_self_adjoint=True)
    assert abs(sop.Tr_blunt().real.item() - 1) < 1e-12 and sop.lambda_min.item() > -1e-12
    for t in (0.1, 0.5, 1.0):
        U = (H * (-1j * t)).expm()
        shifted = alg.operator((U @ sop @ U.adjoint()).matrix)
        assert (shifted - sop).trace_norm().item() <= 2 * t / tau + 1e-9


def test_lemma_A_5_and_corollary_A_7_zero_sets():
    """Lemma A.5 / Corollary A.7: det Z = 0 has probability 0; P(|det Z| < eps) decreases like eps for Gaussian Z."""
    seed(72)
    alg = alg_of([3])
    Z = torch.randn(200000, 3, 3, dtype=RDT)
    d = torch.linalg.det(Z).abs()
    assert d.min() > 0
    f = [(d < e).double().mean().item() for e in (1e-1, 1e-2, 1e-3)]
    assert f[1] < 0.2 * f[0] and f[2] < 0.2 * f[1]
    # Haar unitaries from Ginibre matrices never fail (det Z != 0 a.s.): all unitary
    U = alg.random_unitary(3, batch_size=20000)
    assert (U.mH @ U - torch.eye(3, dtype=DT)).abs().max() < 1e-12
    # Corollary A.7(c): a Gaussian vector has a non-zero component in a given subspace almost surely
    g = torch.randn(100000, 5, dtype=DT)
    assert (g[:, :2].abs().sum(-1) > 0).all()


@pytest.mark.parametrize('ks', [[6], [2, 3, 5], [1, 4]], ids=['6', '2-3-5', '1-4'])
def test_theorem_A_8_courant_fischer_and_A_10_weyl(ks):
    """Theorems A.8 and A.10: max over j-dim U of <x,Ax> >= lambda_j; lambda_j(A) + lambda_1(E) <= lambda_j(A+E) <= lambda_j(A) + lambda_n(E)."""
    seed(73)
    alg = alg_of(ks)
    A, E = hermitian(alg, 50), hermitian(alg, 50) * 0.3
    la, lae, le = A.eigenvalues(), (A + E).eigenvalues(), E.eigenvalues()
    for c in range(len(ks)):
        assert (lae[c] >= la[c] + le[c][:, :1] - 1e-12).all() and (lae[c] <= la[c] + le[c][:, -1:] + 1e-12).all()
        nE = torch.maximum(le[c][:, :1].abs(), le[c][:, -1:].abs())
        assert ((lae[c] - la[c]).abs() <= nE + 1e-12).all()
    # Weyl bound in the algebra: |lambda_j(A+E) - lambda_j(A)| <= ||E|| <= ||E||_2 for the merged spectrum
    merged = lambda X: torch.sort(torch.cat(X.eigenvalues(), dim=-1), dim=-1).values
    assert ((merged(A + E) - merged(A)).abs().amax(-1) <= E.operator_norm() + 1e-12).all()
    assert (E.operator_norm() <= E.frobenius_norm() + 1e-12).all()
    n = ks[0]
    a0 = A.matrix[0, 0, :n, :n]
    lam = torch.linalg.eigvalsh(a0)
    for j in range(1, n + 1):
        Q = alg.random_isometry(n, j, batch_size=200)
        top = torch.linalg.eigvalsh(Q.mH @ a0 @ Q)[:, -1]
        assert (top >= lam[j - 1] - 1e-12).all()
        bottom = torch.linalg.eigvalsh(Q.mH @ a0 @ Q)[:, 0]
        assert (bottom <= lam[n - j] + 1e-12).all()                     # max-min form
    # the optimal subspaces attain the bound
    w, V = torch.linalg.eigh(a0)
    for j in range(1, n + 1):
        assert abs(torch.linalg.eigvalsh(V[:, :j].mH @ a0 @ V[:, :j])[-1] - lam[j - 1]) < 1e-12


@pytest.mark.parametrize('n,m', [(6, 4), (5, 1), (4, 4), (7, 3)])
def test_corollary_A_9_interlacing(n, m):
    """Corollary A.9: lambda_j(A) <= lambda_j(Q^*AQ) <= lambda_{j+n-m}(A) for an isometry Q (n x m)."""
    seed(74)
    alg = alg_of([n])
    A = hermitian(alg, 100).matrix[:, 0]
    Q = alg.random_isometry(n, m, batch_size=100)
    th = torch.linalg.eigvalsh(Q.mH @ A @ Q)
    lam = torch.linalg.eigvalsh(A)
    assert (th >= lam[:, :m] - 1e-12).all() and (th <= lam[:, n - m:] + 1e-12).all()


def test_corollary_A_11_lipschitz_and_simple_spectrum_open():
    """Corollary A.11: each lambda_j is 1-Lipschitz for the operator norm; simple spectra persist under ||B - A|| < gap/2."""
    seed(75)
    alg = alg_of([5])
    A, B = hermitian(alg, 200), hermitian(alg, 200)
    la, lb = A.eigenvalues()[0], B.eigenvalues()[0]
    assert ((la - lb).abs().amax(-1) <= (A - B).operator_norm() + 1e-12).all()
    gap = (la[:, 1:] - la[:, :-1]).amin(-1)
    E = hermitian(alg, 200)
    E = E * (0.49 * gap / E.operator_norm())
    lp = (A + E).eigenvalues()[0]
    assert ((lp[:, 1:] - lp[:, :-1]) > 0).all()


@pytest.mark.parametrize('ks', [[6], [3, 6], [1, 6, 2]], ids=['6', '3-6', '1-6-2'])
def test_lemma_A_12_lagrange_spectral_projections(ks):
    """Lemma A.12: the spectral projection P_i = prod_{l != i} (A - mu_l)/(mu_i - mu_l) equals 1_{mu_i}(A) (functional calculus)."""
    seed(76)
    alg = alg_of(ks)
    mu = torch.tensor([-1.0, 0.5, 2.0])
    A = alg.operator_from_eigenvalues(lambda d: mu[torch.arange(d) % 3], force_self_adjoint=True)
    I = alg.identity()
    total = alg.zero()
    for i in range(3):
        P = I
        for l in range(3):
            if l != i:
                P = P @ ((A - I * mu[l].item()) * (1.0 / (mu[i] - mu[l]).item()))
        ind = A.apply_function(lambda w, m=mu[i].item(): (torch.abs(w - m) < 1e-9).to(w.dtype))
        assert maxdiff(P, ind) < 1e-12
        assert maxdiff(P @ P, P) < 1e-12
        rank = sum((torch.arange(k) % 3 == i).sum().item() for k in ks)
        assert abs(P.Tr_blunt().real.item() - rank) < 1e-12
        total = total + P * mu[i].item()
    assert maxdiff(total, A) < 1e-12                                      # A = sum_i mu_i P_i


def test_proposition_A_13_continuity_of_spectral_projections():
    """Proposition A.13: on simple spectra, eigenvalues and rank-one spectral projections depend continuously on A."""
    seed(77)
    alg = alg_of([4])
    A = alg.operator_from_eigenvalues(lambda d: torch.tensor([-1.0, 0.0, 1.5, 3.0]), force_self_adjoint=True)
    A = alg.operator(A.matrix, is_self_adjoint=True)
    E = hermitian(alg)
    def proj(X):
        w, V = torch.linalg.eigh(X.matrix[0, 0])
        return torch.einsum('ij,kj->jik', V, V.conj())
    P0 = proj(A)
    errs = []
    for t in (1e-2, 1e-3, 1e-4):
        errs.append((proj(alg.operator((A + E * t).matrix, is_self_adjoint=True)) - P0).abs().max().item())
    assert errs[0] > errs[1] > errs[2] and errs[2] < 1e-3
    assert 5 < errs[0] / errs[1] < 20 and 5 < errs[1] / errs[2] < 20     # Lipschitz-like (first order)


def test_proposition_A_14_and_A_16_gaussian_vectors():
    """Propositions A.14, A.16: torch.randn complex is standard complex Gaussian: E z_i conj z_j = delta, E z_i z_j = 0, Uz ~ z."""
    seed(78)
    B, n = 200000, 3
    z = torch.randn(B, n, dtype=DT)
    for i in range(n):
        for j in range(n):
            v = z[:, i] * z[:, j].conj()
            assert abs(zscore(v.real.numpy(), float(i == j))) < 5 and abs(zscore(v.imag.numpy(), 0.0)) < 5
            u = z[:, i] * z[:, j]
            assert abs(zscore(u.real.numpy(), 0.0)) < 5 and abs(zscore(u.imag.numpy(), 0.0)) < 5
    assert abs(zscore((z.real[:, 0] ** 2).numpy(), 0.5)) < 5
    # unitary invariance: the coordinates of Uz have the moments of z (E|.|^4 = 2)
    alg = alg_of([n])
    U = alg.random_unitary(n)
    w = z @ U.T
    assert abs(zscore((w[:, 0].abs() ** 4).numpy(), 2.0)) < 5
    assert abs(zscore((w[:, 0] * w[:, 1].conj()).real.numpy(), 0.0)) < 5
    # Proposition A.14(c): <v, g> ~ N(0, |v|^2) for a real standard Gaussian vector
    g = torch.randn(B, 4, dtype=RDT)
    v = torch.tensor([1.0, -2.0, 0.5, 3.0], dtype=RDT)
    s = g @ v
    assert abs(zscore((s ** 2).numpy(), (v ** 2).sum().item())) < 5
    assert abs(zscore((s ** 4).numpy(), 3 * (v ** 2).sum().item() ** 2)) < 5


def test_theorem_A_15_gaussian_density():
    """Theorem A.15(b): a non-degenerate centred Gaussian X = Ag has density (2 pi)^{-m/2} det(C)^{-1/2} e^{-<x, C^-1 x>/2}."""
    rng = seed(79)
    Amat = rng.standard_normal((2, 4))
    C = Amat @ Amat.T
    X = rng.standard_normal((400000, 4)) @ Amat.T
    h = 1 / 200
    grid = (np.arange(200) + 0.5) * h
    xx = np.stack(np.meshgrid(grid, grid, indexing='ij'), -1).reshape(-1, 2)
    dens = np.exp(-0.5 * np.einsum('ni,ij,nj->n', xx, np.linalg.inv(C), xx)) / (2 * math.pi * math.sqrt(np.linalg.det(C)))
    inside = ((X >= 0) & (X <= 1)).all(axis=1)
    assert abs(zscore(inside, dens.sum() * h * h)) < 5
    # (c) the law depends only on C: B = A O with O orthogonal gives the same covariance and the same probability
    O, _ = np.linalg.qr(rng.standard_normal((4, 4)))
    X2 = rng.standard_normal((400000, 4)) @ (Amat @ O).T
    assert abs(zscore(((X2 >= 0) & (X2 <= 1)).all(axis=1), dens.sum() * h * h)) < 5


@pytest.mark.parametrize('n', [1, 2, 4])
def test_proposition_A_17_goe_gue(n):
    """Proposition A.17: GOE/GUE entries: h_ii ~ N(0,1), h_ij ~ N(0,1/2) (resp. Re, Im each), E Tr H^2 = n(n+1)/2 resp. n^2; invariance."""
    seed(80 + n)
    B = 100000
    G = torch.randn(B, n, n, dtype=RDT)
    goe = (G + G.mT) / 2
    Gc = torch.randn(B, n, n, dtype=DT)
    gue = (Gc + Gc.mH) / math.sqrt(2)
    assert abs(zscore(torch.einsum('bij,bji->b', goe, goe).numpy(), n * (n + 1) / 2)) < 5
    assert abs(zscore(torch.einsum('bij,bji->b', gue, gue).real.numpy(), n * n)) < 5
    assert abs(zscore((goe[:, 0, 0] ** 2).numpy(), 1.0)) < 5 and abs(zscore((gue[:, 0, 0].real ** 2).numpy(), 1.0)) < 5
    if n >= 2:
        assert abs(zscore((goe[:, 0, 1] ** 2).numpy(), 0.5)) < 5
        assert abs(zscore((gue[:, 0, 1].real ** 2).numpy(), 0.5)) < 5
        assert abs(zscore((gue[:, 0, 1].imag ** 2).numpy(), 0.5)) < 5
        assert abs(zscore((gue[:, 0, 1].real * gue[:, 0, 1].imag).numpy(), 0.0)) < 5
        assert abs(zscore((goe[:, 0, 0] * goe[:, 0, 1]).numpy(), 0.0)) < 5
        # invariance: U H U^* has the moments of H (fourth moment of an off-diagonal entry: E|h_12|^4 = 2 for the GUE)
        alg = alg_of([n])
        U = alg.random_unitary(n)
        rot = U @ gue @ U.mH
        assert abs(zscore((rot[:, 0, 1].abs() ** 4).numpy(), 2.0)) < 5
        assert abs(zscore((rot[:, 0, 0].real ** 4).numpy(), 3.0)) < 5
        O = alg_of([n], complex_valued=False).random_unitary(n)
        rot = O @ goe @ O.T
        assert abs(zscore((rot[:, 0, 1] ** 4).numpy(), 3 * 0.25)) < 5
    # the normalising constants of the densities: integrate the GUE density for n = 1, 2 against e^{-Tr H^2/2}
    if n == 2:
        cst_goe = (2 * math.pi) ** (-n / 2) * math.pi ** (-n * (n - 1) / 4)
        cst_gue = (2 * math.pi) ** (-n / 2) * math.pi ** (-n * (n - 1) / 2)
        # product of 1-d Gaussian integrals: (2 pi)^{n/2} for h_ii, pi^{1/2} per off-diagonal real coordinate
        assert abs(cst_goe * (2 * math.pi) ** (n / 2) * math.pi ** (n * (n - 1) / 4) - 1) < 1e-14
        assert abs(cst_gue * (2 * math.pi) ** (n / 2) * math.pi ** (n * (n - 1) / 2) - 1) < 1e-14


def test_corollary_A_19_sphere_areas_and_ball_volumes():
    """Corollary A.19: s_n = n V_n = (2 pi)^{n/2}/I_n = 2 pi^{n/2}/Gamma(n/2), I_{n+2} = n I_n; V_3 = 4 pi/3 by Monte Carlo."""
    from scipy import integrate
    I = lambda n: integrate.quad(lambda r: r ** (n - 1) * math.exp(-r * r / 2), 0, math.inf)[0]
    assert abs(I(1) - math.sqrt(math.pi / 2)) < 1e-10 and abs(I(2) - 1) < 1e-10
    for n in range(1, 9):
        assert abs(I(n + 2) - n * I(n)) < 1e-9 * max(1, I(n + 2))
        assert abs(I(n) - 2 ** (n / 2 - 1) * math.gamma(n / 2)) < 1e-9 * max(1, I(n))
        s = (2 * math.pi) ** (n / 2) / I(n)
        assert abs(s - 2 * math.pi ** (n / 2) / math.gamma(n / 2)) < 1e-9 * s
    for n, s in ((1, 2), (2, 2 * math.pi), (3, 4 * math.pi), (4, 2 * math.pi ** 2)):
        assert abs(2 * math.pi ** (n / 2) / math.gamma(n / 2) - s) < 1e-12
    rng = seed(81)
    x = rng.uniform(-1, 1, (400000, 3))
    assert abs(zscore((x ** 2).sum(1) <= 1, (4 * math.pi / 3) / 8)) < 5
    x = rng.uniform(-1, 1, (400000, 4))
    assert abs(zscore((x ** 2).sum(1) <= 1, (2 * math.pi ** 2 / 4) / 16)) < 5      # V_4 = s_4 / 4


@pytest.mark.parametrize('N', [1, 2, 3, 5])
def test_theorem_A_20_and_corollary_A_21_complex_direction(N):
    """Theorem A.20(b), Corollary A.21: |z|^2 ~ Gamma(N), E|z|^4 = N(N+1), z/|z| uniform and independent of |z|; moments of u."""
    seed(90 + N)
    B = 200000
    z = torch.randn(B, N, dtype=DT)
    R2 = (z.abs() ** 2).sum(-1)
    u = z / R2.sqrt()[:, None]
    t = 1.5
    p_gamma = 1 - math.exp(-t) * sum(t ** k / math.factorial(k) for k in range(N))
    assert abs(zscore((R2 <= t).double().numpy(), p_gamma)) < 5
    assert abs(zscore((R2 ** 2).numpy(), N * (N + 1))) < 5
    assert abs(zscore((u[:, 0].abs() ** 2).numpy(), 1 / N)) < 5 or N == 1
    assert abs(zscore((u[:, 0].abs() ** 4).numpy(), 2 / (N * (N + 1)))) < 5 or N == 1
    if N >= 2:
        assert abs(zscore(((u[:, 0] * u[:, 1]).abs() ** 2).numpy(), 1 / (N * (N + 1)))) < 5
        assert abs(zscore(((R2 <= t).double() * u[:, 0].abs() ** 4).numpy(), p_gamma * 2 / (N * (N + 1)))) < 5
        assert abs(zscore((u[:, 0] * u[:, 1].conj()).real.numpy(), 0.0)) < 5
    # the phase of z_1 is uniform (sigma_1 = normalised arc length)
    ph = torch.angle(z[:, 0])
    for lo, hi in ((0.0, 1.0), (-3.0, -2.0), (2.0, math.pi)):
        assert abs(zscore(((ph >= lo) & (ph < hi)).double().numpy(), (hi - lo) / (2 * math.pi))) < 5


@pytest.mark.parametrize('n', [2, 3, 6])
def test_corollary_A_21_real_sphere(n):
    """Corollary A.21 (real case): E u_1^2 = 1/n, E u_1^4 = 3/(n(n+2)), E u_1^2 u_2^2 = 1/(n(n+2)); E|g|^4 = n(n+2)."""
    seed(100 + n)
    g = torch.randn(200000, n, dtype=RDT)
    r2 = (g ** 2).sum(-1)
    u = g / r2.sqrt()[:, None]
    assert abs(zscore((u[:, 0] ** 2).numpy(), 1 / n)) < 5
    assert abs(zscore((u[:, 0] ** 4).numpy(), 3 / (n * (n + 2)))) < 5
    assert abs(zscore((u[:, 0] ** 2 * u[:, 1] ** 2).numpy(), 1 / (n * (n + 2)))) < 5
    assert abs(zscore((r2 ** 2).numpy(), n * (n + 2))) < 5


@pytest.mark.parametrize('N', [2, 3, 5])
def test_proposition_A_22_haar_columns_uniform(N):
    """Proposition A.22 / Corollary A.21: Qv (Q Haar, v a fixed unit vector) is uniform on the sphere: E|u_1|^4 = 2/(N(N+1))."""
    seed(110 + N)
    alg = alg_of([N])
    Q = alg.random_unitary(N, batch_size=100000)
    v = torch.randn(N, dtype=DT)
    v = v / v.norm()
    u = Q @ v
    assert maxdiff(u.norm(dim=-1), torch.ones(100000, dtype=RDT)) < 1e-12
    assert abs(zscore((u[:, 0].abs() ** 4).numpy(), 2 / (N * (N + 1)))) < 5
    assert abs(zscore(((u[:, 0] * u[:, 1]).abs() ** 2).numpy(), 1 / (N * (N + 1)))) < 5
    assert abs(zscore((u[:, 0] * u[:, 1].conj()).real.numpy(), 0.0)) < 5
    # rows of a Haar unitary as well (left and right invariance)
    r = Q[:, 0, :]
    assert abs(zscore((r[:, -1].abs() ** 4).numpy(), 2 / (N * (N + 1)))) < 5


def test_lemma_A_23_moment_generating_function():
    """Lemma A.23: E[Y_i Y_j] = c_ij and E e^{<t,Y>} = e^{<t,Ct>/2} for a centred Gaussian Y = Ag."""
    rng = seed(120)
    Amat = rng.standard_normal((3, 5)) / 2
    C = Amat @ Amat.T
    Y = rng.standard_normal((400000, 5)) @ Amat.T
    for i in range(3):
        for j in range(3):
            assert abs(zscore(Y[:, i] * Y[:, j], C[i, j])) < 5
    t = np.array([0.3, -0.2, 0.25])
    assert abs(zscore(np.exp(Y @ t), math.exp(t @ C @ t / 2))) < 5


def _pairings(idx):
    if not idx:
        yield []
        return
    for i in range(1, len(idx)):
        for rest in _pairings(idx[1:i] + idx[i + 1:]):
            yield [(idx[0], idx[i])] + rest


def _wick(cov, idx):
    return sum(math.prod(cov[idx[a], idx[b]] for a, b in p) for p in _pairings(list(range(len(idx)))))


def test_theorem_A_24_isserlis_wick():
    """Theorem A.24: E[Y_1...Y_2p] = sum over (2p-1)!! pairings of products of covariances; odd moments vanish; complex form."""
    assert [len(list(_pairings(list(range(2 * p))))) for p in range(1, 6)] == [1, 3, 15, 105, 945]
    assert all(len(list(_pairings(list(range(2 * p))))) == math.factorial(2 * p) // (2 ** p * math.factorial(p))
               for p in range(1, 6))
    rng = seed(121)
    Amat = rng.standard_normal((3, 3)) / math.sqrt(3)
    cov = Amat @ Amat.T
    Y = rng.standard_normal((1000000, 3)) @ Amat.T
    for idx in ((0, 1, 2, 0), (0, 0, 0, 0), (0, 1, 1, 2, 2, 0), (1, 1), (2, 2, 2, 2, 2, 2)):
        assert abs(zscore(np.prod(Y[:, list(idx)], axis=1), _wick(cov, idx))) < 5
    for idx in ((0, 1, 2), (0, 0, 0), (0, 1, 1, 2, 2)):
        assert abs(zscore(np.prod(Y[:, list(idx)], axis=1), 0.0)) < 5
    # complex form: (z, conj z) = M (x, y), E[Z_a Z_b] = [[0, 1], [1, 0]]; E|z|^4 = 2, E|z|^6 = 6, E|z|^8 = 24
    M = np.array([[1, 1j], [1, -1j]])
    cz = M @ M.T / 2
    assert np.abs(cz - np.array([[0, 1], [1, 0]])).max() < 1e-15
    for p in (1, 2, 3, 4):
        assert abs(_wick(cz, (0, 1) * p) - math.factorial(p)) < 1e-12
    z = torch.randn(400000, dtype=DT)
    for p in (2, 3):
        assert abs(zscore((z.abs() ** (2 * p)).numpy(), math.factorial(p))) < 5


def test_section_A_6_gue_sixth_moment_monte_carlo():
    """Theorem A.24 applied as in Proposition 1.10: E Tr H^6 = 5n^4 + 10n^2 for the GUE (Monte Carlo through the library)."""
    seed(122)
    n, B = 3, 100000
    alg = alg_of([n])
    G = torch.randn(B, 1, n, n, dtype=DT)
    H = alg.operator((G + G.mH) / math.sqrt(2), is_self_adjoint=True)
    H3 = H @ H @ H
    m6 = (H3 @ H3).trace.real.numpy()
    assert abs(zscore(m6, 5 * n ** 4 + 10 * n ** 2)) < 5
    # the same through the functional calculus (eigenvalues)
    m6e = (torch.cat(H.eigenvalues(), -1) ** 6).sum(-1).numpy()
    assert np.abs(m6 - m6e).max() < 1e-8 * np.abs(m6).max()


def test_theorem_A_25_normalised_exponentials():
    """Theorem A.25: (e_0..e_n)/T is independent of T ~ Gamma(n+1); the partial sums / T have density n! on the simplex."""
    rng = seed(123)
    n, B = 4, 400000
    e = rng.exponential(1.0, (B, n + 1))
    T = e.sum(1)
    w = e / T[:, None]
    assert abs(zscore(w[:, 0] * T ** 2, (n + 2))) < 5                # E w_0 E T^2 = (n+1)(n+2)/(n+1)
    assert abs(zscore((T <= 4.0).astype(float), 1 - math.exp(-4) * sum(4 ** k / math.factorial(k) for k in range(n + 1)))) < 5
    assert abs(zscore(((w[:, 0] < 0.1) & (T > 5)).astype(float), (w[:, 0] < 0.1).mean() * (T > 5).mean())) < 5
    # lambda_n(Delta_n) = 1/n!: a uniform point of [0,1]^n is sorted with probability 1/n!
    x = rng.uniform(0, 1, (B, n))
    assert abs(zscore((np.diff(x, axis=1) > 0).all(1), 1 / math.factorial(n))) < 5


@pytest.mark.parametrize('n', [1, 2, 5])
def test_corollary_A_26_uniform_spacings(n):
    """Corollary A.26: spacings of n uniform points: E s = 1/(n+1), E s^2 = 2/((n+1)(n+2)), E s_i s_j = 1/((n+1)(n+2)), ratio laws."""
    rng = seed(130 + n)
    B = 200000
    lv = np.sort(rng.uniform(0, 1, (B, n)), axis=1)
    s = np.diff(np.concatenate([np.zeros((B, 1)), lv, np.ones((B, 1))], axis=1), axis=1)
    for i in (0, n):
        assert abs(zscore(s[:, i], 1 / (n + 1))) < 5
        assert abs(zscore(s[:, i] ** 2, 2 / ((n + 1) * (n + 2)))) < 5
    assert abs(zscore(s[:, 0] * s[:, n], 1 / ((n + 1) * (n + 2)))) < 5
    for t in (0.3, 1.0, 4.0):
        assert abs(zscore(s[:, 0] / s[:, 1] <= t, t / (1 + t))) < 5
    r = np.minimum(s[:, 0], s[:, 1]) / np.maximum(s[:, 0], s[:, 1])
    for t in (0.25, 0.5, 0.9):
        assert abs(zscore(r <= t, 2 * t / (1 + t))) < 5
    assert abs(zscore(r, 2 * math.log(2) - 1)) < 5
    # affine invariance: points on [a, b]: spacings scale by b - a, ratios unchanged
    a, b = -2.0, 5.0
    lv2 = a + (b - a) * lv
    s2 = np.diff(np.concatenate([np.full((B, 1), a), lv2, np.full((B, 1), b)], axis=1), axis=1)
    assert np.abs(s2 / (b - a) - s).max() < 1e-12


def test_section_A_7_poisson_levels_ratio_law():
    """Section A.7: the consecutive spacing ratios of Poisson levels and of uniform spacings have the same joint law."""
    rng = seed(140)
    n, B = 5, 200000
    e = rng.exponential(1.0, (B, n + 1))
    lv = np.sort(rng.uniform(0, 1, (B, n)), axis=1)
    s = np.diff(np.concatenate([np.zeros((B, 1)), lv, np.ones((B, 1))], axis=1), axis=1)
    ratio = lambda d: np.minimum(d[:, :-1], d[:, 1:]) / np.maximum(d[:, :-1], d[:, 1:])
    re, rs = ratio(e), ratio(s)
    # joint statistic of two consecutive ratios and the marginal law 2t/(1+t)
    for stat, x in (('prod', lambda r: r[:, 0] * r[:, 1]), ('min', lambda r: np.minimum(r[:, 1], r[:, 2]))):
        a, b = x(re), x(rs)
        z = (a.mean() - b.mean()) / math.sqrt(a.var() / B + b.var() / B)
        assert abs(z) < 5
    assert abs(zscore(re[:, 2] <= 0.5, 2 / 3)) < 5
    # the first n Poisson levels divided by the (n+1)-st are n sorted uniform points
    lev = np.cumsum(e, axis=1)
    x = lev[:, :n] / lev[:, n:]
    for k in range(n):
        assert abs(zscore(x[:, k], (k + 1) / (n + 1))) < 5           # E U_(k) = k/(n+1)


# =====================================================================================================================
# Further cases: normal (non-self-adjoint) operators, non-normal operators, real algebras, the decomposition maps
# =====================================================================================================================

@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_theorem_1_5_normal_operators_degenerate_and_near_degenerate(ks):
    """Theorem 1.5 for normal A with complex (degenerate, nearly degenerate) spectrum: z^2, conj(z), exp(z) give A^2, A^*, e^A."""
    seed(150)
    alg = alg_of(ks)
    spectra = [torch.exp(1j * torch.tensor([0.0, 0.0, 1.5, 1.5, 3.0], dtype=RDT)),
               torch.tensor([1 + 1j, 1 + 1j + 1e-10, -2.0, 0.5j, 3.0], dtype=DT)]
    for lam in spectra:
        A = alg.operator_from_eigenvalues(lambda d: lam[:d])
        A2 = alg.operator(A.matrix.clone())                              # untagged: normality is detected
        for X in (A, A2):
            assert X.is_normal
            assert maxdiff(X.apply_function(lambda z: z ** 2), X @ X) < 1e-12
            assert maxdiff(X.apply_function(torch.conj), X.adjoint()) < 1e-12
            assert maxdiff(X.apply_function(torch.exp), X.expm()) < 1e-12
            # spec f(A) = f(spec A) for the unitary case (|lambda| = 1): f(A) = A^3 is unitary
            if lam.abs().sub(1).abs().max() < 1e-15:
                F = X.apply_function(lambda z: z ** 3)
                assert maxdiff(F.adjoint() @ F, alg.identity()) < 1e-12


@pytest.mark.parametrize('ks', LAYOUTS, ids=LAYOUT_IDS)
def test_theorem_1_5_b_polynomials_of_non_normal_operators(ks):
    """Theorem 1.5(b) beyond normal operators: the holomorphic calculus of a polynomial is the polynomial in A (non-normal A)."""
    seed(151)
    alg = alg_of(ks)
    A = ginibre(alg, 2) * 0.5
    I = alg.identity(2)
    P = A.holomorphic_function(lambda z: z ** 3 - 2 * z + 1)
    assert maxdiff(P, A @ A @ A - A * 2.0 + I) < 1e-9
    assert maxdiff(A.holomorphic_function(torch.exp), A.expm()) < 1e-9


@pytest.mark.parametrize('ks', [[3], [2, 3, 5], [1, 4, 2]], ids=['3', '2-3-5', '1-4-2'])
def test_chapter_1_3_real_algebras(ks):
    """Chapters 1 and 3 over the reals (O(k) blocks): norms, polar decomposition, traces, Klein, densities."""
    seed(152)
    alg = alg_of(ks, complex_valued=False)
    G = alg.operator(torch.randn(3, alg.C, alg.k_max, alg.k_max, dtype=RDT))
    for b in range(3):
        s = torch.linalg.svdvals(dense(G, b=b))
        assert abs(G.operator_norm()[b] - s[0]) < 1e-12 and abs(G.trace_norm()[b] - s.sum()) < 1e-11
    W, P = G.polar()
    assert maxdiff(W @ P, G) < 1e-10 and maxdiff(W.adjoint() @ W, alg.identity(3)) < 1e-10
    U = alg.random_unitary_operator(3)
    assert U.matrix.dtype == RDT and maxdiff(U.adjoint() @ U, alg.identity(3)) < 1e-12
    H1, H2 = alg.operator((G + G.adjoint()).matrix, is_self_adjoint=True), hermitian(alg_of(ks), 3)
    X = alg.operator(H2.matrix.real.contiguous(), is_self_adjoint=True)
    val = (H1.expm() - X.expm() - (H1 - X) @ X.expm()).Tr_blunt()
    assert (val >= -1e-10).all()
    rho, sig = states.random_density_matrix(alg, 3), states.random_density_matrix(alg, 3)
    assert (states.relative_entropy(rho, sig) > 0).all()
    for name, tr in (('norm', 'Tr_norm'), ('tau_vN', 'tau_vN')):
        d = alg.operator(rho.density(name), is_self_adjoint=True)
        assert maxdiff(getattr(d @ H1, tr)(), rho.expectation(H1)) < 1e-12


@pytest.mark.parametrize('ks,ms', STRUCTURES, ids=[f"k{k}m{m}" for k, m in STRUCTURES])
def test_theorem_2_6_decomposition_maps(ks, ms):
    """Theorem 2.6: V^* M V = (+)_c M_{k_c} (x) 1_{m_c}, V^* M' V = (+)_c 1_{k_c} (x) M_{m_c}, V^* P_c V = projection onto summand c."""
    seed(153)
    gens, _ = _random_vn_algebra(ks, ms, s=7 + sum(ks))
    dec = decompose.wedderburn(gens)
    A = dec.to_algebra(torch.stack(gens))
    assert maxdiff(dec.from_algebra(A), torch.stack(gens)) < 1e-9          # gens lie in V((+) M_k (x) 1_m)V^*
    # random elements of the commutant: V (+)_c (1 (x) x_c) V^* commutes with M
    calg = dec.commutant_algebra()
    X = dec.from_commutant(ginibre(calg, 3))
    for g in gens:
        assert (g @ X - X @ g).abs().max() < 1e-9
    # the minimal central projections are the central elements e_c of the abstract algebra
    Pc = dec.central_projections()
    for c in range(dec.C):
        e = [1.0 if d == c else 0.0 for d in range(dec.C)]
        assert maxdiff(dec.to_algebra(Pc[c]), dec.algebra.central(e)) < 1e-9
        assert maxdiff(dec.to_commutant(Pc[c]), calg.central(e)) < 1e-9
    # the product of the abstract algebra is the product in M_n (phi is a *-isomorphism)
    a, b = ginibre(dec.algebra), ginibre(dec.algebra)
    assert maxdiff(dec.from_algebra(a @ b), dec.from_algebra(a) @ dec.from_algebra(b)) < 1e-9
    assert maxdiff(dec.from_algebra(a.adjoint()), dec.from_algebra(a).mH) < 1e-12
    # Section 2.3: Tr_{C^n} = sum_c m_c Tr a_c
    W = dec.algebra.central([float(m) for m in dec.multiplicities])
    assert abs(torch.trace(dec.from_algebra(a)[0]) - (W @ a).Tr_blunt()[0]) < 1e-9


def test_theorem_2_6_uniqueness_of_pairs_under_change_of_basis():
    """Theorem 2.6, uniqueness: the pairs (k_c, m_c) do not depend on the unitary V (two bases of the same algebra)."""
    seed(207)
    gens, V = _random_vn_algebra((2, 1, 3), (1, 2, 1), s=9)
    W = alg_of([7]).random_unitary(7)
    gens2 = [W @ g @ W.mH for g in gens]
    d1, d2 = decompose.wedderburn(gens), decompose.wedderburn(gens2, seed=5)
    assert sorted(zip(d1.k_factors, d1.multiplicities)) == sorted(zip(d2.k_factors, d2.multiplicities))
    for P in d1.central_projections():                              # k_c^2 = dim M P_c, m_c = rank P_c / k_c
        r = round(torch.trace(P).real.item())
        assert r in [k * m for k, m in zip(d1.k_factors, d1.multiplicities)]
