"""
Theory tests for Chapter 5 (random matrices in an algebra) and Appendix C (eigenvalue densities and Page's
formula) of the lecture notes, checked through the library: Haar measure by corrected QR, U(M) and SU(k), the
twirl, the Weingarten calculus of order two and three, the circular ensembles, the Wigner surmise, the spacing
ratio, the absence of repulsion between sectors, the Weyl chart, the Gaussian and Laguerre level densities, the
Bartlett decomposition, the induced measure, Page's formula and the purity of random states.

Monte Carlo checks compare a sample mean with an exact value in units of its standard error (|z| <= 5, fixed
seeds); exact identities are checked to machine precision.
"""
import itertools
import math
from fractions import Fraction

import numpy as np
import pytest
import torch
from scipy import integrate, special

from torch_vn_algebra import TypeIAlgebra, SpinChain, channels, states, cost
from torch_vn_algebra.decompose import commutant_basis, wedderburn

pytestmark = pytest.mark.filterwarnings("ignore", category=cost.CostWarning)

CD, RD = torch.complex128, torch.float64
Z_MAX = 5.0
ratio = SpinChain.level_spacing_ratio


# ----------------------------------------------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------------------------------------------

def zscore(samples, exact):
    s = torch.as_tensor(samples, dtype=RD).flatten()
    se = s.std().item() / math.sqrt(s.numel())
    return (s.mean().item() - exact) / se if se > 0 else (0.0 if abs(s.mean().item() - exact) < 1e-12 else math.inf)


def assert_mc(samples, exact, name="", z=Z_MAX):
    """|mean - exact| <= z standard errors (plus 1e-12 for quantities that are constant up to rounding)."""
    s = torch.as_tensor(samples, dtype=RD).flatten()
    mean, se = s.mean().item(), s.std().item() / math.sqrt(s.numel())
    assert abs(mean - exact) <= z * se + 1e-12, f"{name}: z = {zscore(s, exact):.2f} (mean {mean:.6g}, exact {exact:.6g})"


def assert_mc_complex(samples, exact, name="", z=Z_MAX):
    assert_mc(samples.real, complex(exact).real, name + " (Re)", z)
    assert_mc(samples.imag, complex(exact).imag, name + " (Im)", z)


def alg_c(ks, ns=None):
    return TypeIAlgebra(list(ns or ks), list(ks), precision='double')


def alg_r(ks, ns=None):
    return TypeIAlgebra(list(ns or ks), list(ks), complex_valued=False, precision='double')


def haar(k, B, seed):
    torch.manual_seed(seed)
    return alg_c([k]).random_unitary(k, batch_size=B)


def swap(k):
    F = torch.zeros(k * k, k * k, dtype=CD)
    for i in range(k):
        for j in range(k):
            F[j * k + i, i * k + j] = 1
    return F


def n_cycles(p):
    seen, count = set(), 0
    for i in range(len(p)):
        if i not in seen:
            count += 1
            while i not in seen:
                seen.add(i)
                i = p[i]
    return count


def compose(a, b):                       # (a o b)(i) = a(b(i))
    return tuple(a[b[i]] for i in range(len(b)))


def inverse(a):
    inv = [0] * len(a)
    for i, ai in enumerate(a):
        inv[ai] = i
    return tuple(inv)


def weingarten(n, k):
    """Wg on S_n (dict perm -> value) from the pseudo-inverse of the Gram matrix k^{#cycles(s^-1 t)} (Thm 5.12)."""
    perms = list(itertools.permutations(range(n)))
    G = np.array([[float(k) ** n_cycles(compose(inverse(s), t)) for t in perms] for s in perms])
    W = np.linalg.pinv(G)
    e = perms.index(tuple(range(n)))
    return perms, G, W, {p: W[e, j] for j, p in enumerate(perms)}


def wg_moment(i, j, ip, jp, k):
    """E[U_{i1j1}..U_{injn} conj(U_{i'1j'1}..U_{i'nj'n})] by the Weingarten formula of Theorem 5.12."""
    n = len(i)
    perms, _, _, wg = weingarten(n, k)
    tot = 0.0
    for s in perms:
        if any(i[m] != ip[s[m]] for m in range(n)):
            continue
        for p in perms:
            if any(j[m] != jp[p[m]] for m in range(n)):
                continue
            tot += wg[compose(s, inverse(p))]
    return tot


def goe(B, n, seed):
    torch.manual_seed(seed)
    G = torch.randn(B, n, n, dtype=RD)
    return (G + G.mT) / 2


def gue(B, n, seed):
    torch.manual_seed(seed)
    G = torch.randn(B, n, n, dtype=CD)
    return (G + G.mH) / math.sqrt(2)


def vandermonde(lam):
    n = lam.shape[-1]
    out = torch.ones(lam.shape[:-1], dtype=lam.dtype)
    for a in range(n):
        for b in range(a + 1, n):
            out = out * (lam[..., b] - lam[..., a])
    return out


def laguerre(k, al, x):
    return sum((-1) ** j * math.comb(k + al, k - j) * x ** j / math.factorial(j) for j in range(k + 1))


def rho_laguerre(x, a, b):
    al = b - a
    return sum(math.factorial(k) / math.factorial(k + al) * laguerre(k, al, x) ** 2 for k in range(a)) * x ** al * np.exp(-x)


def gue_density(n):
    """GUE level density K_n(x) e^{-x^2/2} (Lemma C.10, Hermite weight), as a numpy function."""
    def rho(x):
        x = np.asarray(x, dtype=float)
        K = np.zeros_like(x)
        for m in range(n):
            He = special.eval_hermitenorm(m, x)
            K += He ** 2 / (math.sqrt(2 * math.pi) * math.factorial(m))
        return K * np.exp(-x ** 2 / 2)
    return rho


def harmonic(m):
    return sum(Fraction(1, i) for i in range(1, m + 1))


def page(a, b):
    return harmonic(a * b) - harmonic(b) - Fraction(a - 1, 2 * b)


# ================================================================================================================
# 5.1 Haar measure
# ================================================================================================================

@pytest.mark.parametrize("k", [1, 2, 3, 5])
def test_theorem_5_3_sampler_is_phase_corrected_qr(k):
    """Theorem 5.3: random_unitary is Q' diag(R'_ii/|R'_ii|) of the QR of the same Ginibre matrix, and Z = QR with diag(R) > 0."""
    alg = alg_c([k])
    torch.manual_seed(11)
    Z = torch.randn(64, k, k, dtype=CD)
    torch.manual_seed(11)
    U = alg.random_unitary(k, batch_size=64)
    Q, R = torch.linalg.qr(Z)
    d = R.diagonal(dim1=-2, dim2=-1)
    assert torch.allclose(U, Q * (d / d.abs()).unsqueeze(-2), atol=1e-13)
    assert torch.allclose(U.mH @ U, torch.eye(k, dtype=CD).expand(64, k, k), atol=1e-13)
    Rp = U.mH @ Z                                                   # Lemma 5.2: the unique R
    assert torch.allclose(Rp, Rp.triu(), atol=1e-12)
    diag = Rp.diagonal(dim1=-2, dim2=-1)
    assert diag.imag.abs().max() < 1e-12 and diag.real.min() > 0


@pytest.mark.parametrize("k", [1, 2, 3])
def test_theorem_5_3_real_case_sign_corrected_qr(k):
    """Theorem 5.3 (real case, Theorem B.17): in a real algebra random_unitary is O' diag(sgn R'_ii), orthogonal with E O = 0, E O_ij^2 = 1/k."""
    alg = alg_r([k])
    torch.manual_seed(12)
    Z = torch.randn(50_000, k, k, dtype=RD)
    torch.manual_seed(12)
    O = alg.random_unitary(k, batch_size=50_000)
    Q, R = torch.linalg.qr(Z)
    assert torch.allclose(O, Q * R.diagonal(dim1=-2, dim2=-1).sign().unsqueeze(-2), atol=1e-13)
    assert torch.allclose(O.mT @ O, torch.eye(k, dtype=RD).expand_as(O), atol=1e-12)
    for i in range(k):
        for j in range(k):
            assert_mc(O[:, i, j], 0.0, f"E O_{i}{j}")
            assert_mc(O[:, i, j] ** 2, 1 / k, f"E O_{i}{j}^2")


def test_theorem_5_3_uncorrected_qr_is_not_haar():
    """Remark after Theorem 5.3: without the phase correction the LAPACK Q' has E Re Tr Q' far from 0, so it is not Haar."""
    torch.manual_seed(13)
    Q, _ = torch.linalg.qr(torch.randn(50_000, 2, 2, dtype=CD))
    tr = Q.diagonal(dim1=-2, dim2=-1).sum(-1).real
    assert abs(zscore(tr, 0.0)) > 20


@pytest.mark.parametrize("k", [1, 2, 3, 4])
def test_corollary_5_5_first_and_second_moments_all_entries(k):
    """Corollary 5.5: E U_ij = 0 and E|U_ij|^2 = 1/k for every entry of a Haar unitary."""
    U = haar(k, 60_000, 21 + k)
    for i in range(k):
        for j in range(k):
            assert_mc_complex(U[:, i, j], 0.0, f"E U_{i}{j}")
            assert_mc(U[:, i, j].abs() ** 2, 1 / k, f"E|U_{i}{j}|^2")


@pytest.mark.parametrize("k", [2, 3])
def test_lemma_5_1_left_and_right_invariance(k):
    """Lemma 5.1 / Theorem 5.3: VU, UV and U^* have the law of U (checked on E|(.)_11|^4 = 2/(k(k+1)) and E Re (.)_12 (.)_21-bar = 0)."""
    U = haar(k, 100_000, 30 + k)
    torch.manual_seed(99)
    V = alg_c([k]).random_unitary(k)                     # a fixed unitary
    for X, name in ((V @ U, "VU"), (U @ V, "UV"), (U.mH, "U*")):
        assert_mc(X[:, 0, 0].abs() ** 4, 2 / (k * (k + 1)), f"{name}: E|X_11|^4")
        assert_mc_complex(X[:, 0, 1] * X[:, 1, 0].conj(), 0.0, f"{name}: E X_12 conj X_21")
        assert_mc(X[:, 0, 0].abs() ** 2 * X[:, 1, 1].abs() ** 2, wg_moment((0, 1), (0, 1), (0, 1), (0, 1), k),
                  f"{name}: E|X_11|^2|X_22|^2")


@pytest.mark.parametrize("ks,ns", [([1], [1]), ([3], [3]), ([1, 2, 3], [1, 2, 3]), ([2, 2], [2, 2]),
                                   ([4, 1], [4, 1]), ([2, 3], [3, 5])])
def test_corollary_5_4_random_unitary_operator_is_in_U_M(ks, ns):
    """Corollary 5.4: random_unitary_operator returns unitaries of M = (+)_c M_{k_c} (zero padding, U*U = UU* = 1)."""
    alg = alg_c(ks, ns)
    torch.manual_seed(41)
    U = alg.random_unitary_operator(32)
    one = alg.identity(32).matrix
    assert torch.allclose((U.adjoint() @ U).matrix, one, atol=1e-12)
    assert torch.allclose((U @ U.adjoint()).matrix, one, atol=1e-12)
    for c, k in enumerate(ks):
        if k < max(ks):                                           # zero padding outside the k_c x k_c block
            assert U.matrix[:, c, k:, :].abs().max().item() == 0
            assert U.matrix[:, c, :, k:].abs().max().item() == 0


def test_corollary_5_4_blocks_independent_and_haar():
    """Corollary 5.4: the blocks of a Haar unitary of U(M_1 + M_2 + M_3) are Haar on U(k_c) and independent."""
    alg = alg_c([1, 2, 3])
    torch.manual_seed(42)
    U = alg.random_unitary_operator(100_000).matrix
    u = [U[:, c, :k, :k] for c, k in enumerate(alg.k_factors)]
    for c, k in enumerate(alg.k_factors):
        assert_mc(u[c][:, 0, 0].abs() ** 2, 1 / k, f"E|(U_{c})_11|^2")
        if k >= 2:
            assert_mc(u[c][:, 0, 0].abs() ** 4, 2 / (k * (k + 1)), f"E|(U_{c})_11|^4")
    assert_mc(u[1][:, 0, 0].abs() ** 2 * u[2][:, 0, 0].abs() ** 2, 1 / 6, "E|U_1 11|^2 |U_2 11|^2")
    assert_mc(u[0][:, 0, 0].abs() ** 2 * u[2][:, 1, 2].abs() ** 2, 1 / 3, "E|U_0 11|^2 |U_2 23|^2")
    assert_mc_complex(u[1][:, 0, 0] * u[2][:, 0, 0].conj(), 0.0, "E U_1 11 conj U_2 11")
    assert_mc_complex(u[0][:, 0, 0] * u[1][:, 1, 1].conj(), 0.0, "E U_0 11 conj U_1 22")


def test_section_5_1_haar_of_U_D_leaves_M_almost_surely():
    """Section 5.1: a Haar unitary of U(D) has all first-column entries nonzero a.s., so it is not in M when C >= 2."""
    U = haar(5, 20_000, 43)
    assert U[:, :, 0].abs().min().item() > 0


@pytest.mark.parametrize("k", [2, 3])
def test_section_5_1_haar_su(k):
    """Section 5.1: measure='haar_su' gives det V = 1, E V_ij = 0, E|V_ij|^2 = 1/k and E[V_11...V_kk] = 1/k! (as E det V = 1)."""
    alg = alg_c([k])
    torch.manual_seed(50 + k)
    V = alg.random_unitary_operator(100_000, measure='haar_su').matrix[:, 0]
    assert torch.allclose(torch.linalg.det(V), torch.ones(V.shape[0], dtype=CD), atol=1e-12)
    for i in range(k):
        assert_mc_complex(V[:, i, (i + 1) % k], 0.0, "E V_ij")
        assert_mc(V[:, i, 0].abs() ** 2, 1 / k, "E|V_i1|^2")
    prod = torch.ones(V.shape[0], dtype=CD)
    for i in range(k):
        prod = prod * V[:, i, i]
    assert_mc_complex(prod, 1 / math.factorial(k), "E prod V_ii")
    # left invariance under a fixed W in SU(k): same moment for W V
    torch.manual_seed(7)
    W = alg.random_unitary(k, measure='haar_su')
    WV = W @ V
    prod = torch.ones(V.shape[0], dtype=CD)
    for i in range(k):
        prod = prod * WV[:, i, i]
    assert_mc_complex(prod, 1 / math.factorial(k), "E prod (WV)_ii")


@pytest.mark.parametrize("k", [2, 3])
def test_section_5_1_haar_su_real_is_so(k):
    """Section 5.1 (real case): measure='haar_su' of a real algebra gives SO(k): orthogonal, det = +1, E O_11^2 = 1/k."""
    alg = alg_r([k])
    torch.manual_seed(52)
    O = alg.random_unitary(k, measure='haar_su', batch_size=50_000)
    assert torch.allclose(torch.linalg.det(O), torch.ones(O.shape[0], dtype=RD), atol=1e-12)
    assert_mc(O[:, 0, 0] ** 2, 1 / k, "E O_11^2")
    assert_mc(O[:, 0, 1], 0.0, "E O_12")


# ================================================================================================================
# 5.2 Moments and the Weingarten calculus
# ================================================================================================================

@pytest.mark.parametrize("k", [1, 2, 3, 4])
def test_lemma_5_6_commutant_of_diagonal_unitaries_and_permutations(k):
    """Lemma 5.6: the commutant of the diagonal unitaries and the permutation matrices is C 1."""
    gens = [torch.diag(torch.exp(1j * torch.arange(1, k + 1, dtype=RD))).to(CD)]
    gens += [torch.eye(k, dtype=CD)[list(p)] for p in itertools.permutations(range(k))]
    basis = commutant_basis(gens)
    assert basis.shape[0] == 1
    b = basis[0] / basis[0][0, 0]
    assert torch.allclose(b, torch.eye(k, dtype=CD), atol=1e-10)


@pytest.mark.parametrize("ks", [[1], [3], [1, 2, 3], [2, 4]])
def test_proposition_5_7_twirl_is_center_expectation(ks):
    """Proposition 5.7: E[U X U*] over U Haar on U(M) is E(X) = (+)_c (Tr X_c / k_c) 1_{k_c}, which is channels.center_expectation."""
    alg = alg_c(ks)
    torch.manual_seed(70)
    X = alg.from_blocks([torch.randn(k, k, dtype=CD) for k in ks])
    EX = channels.center_expectation(alg)(X).matrix[0]
    exact = torch.zeros_like(EX)
    for c, k in enumerate(ks):
        exact[c, :k, :k] = torch.eye(k, dtype=CD) * X.matrix[0, c, :k, :k].trace() / k
    assert torch.allclose(EX, exact, atol=1e-12)
    U = alg.random_unitary_operator(60_000)
    T = (U @ X @ U.adjoint()).matrix
    for c, k in enumerate(ks):
        for i in range(k):
            for j in range(k):
                assert_mc_complex(T[:, c, i, j], exact[c, i, j].item(), f"twirl entry {c},{i},{j}")


def test_proposition_5_7_haar_on_U_D_loses_sector_weights():
    """Remarks after Proposition 5.7: Haar on U(D) gives Tr X / D times 1; the pi/4 rotation maps P_1 out of M (entry 1/2 at (1, k_1 + 1))."""
    D, k1 = 4, 2
    torch.manual_seed(71)
    X = torch.diag(torch.tensor([3.0, -1.0, 0.5, 2.0], dtype=RD)).to(CD) + 0.3 * torch.randn(D, D, dtype=CD)
    U = haar(D, 60_000, 72)
    T = U @ X @ U.mH
    for i in range(D):
        for j in range(D):
            assert_mc_complex(T[:, i, j], (X.trace() / D).item() if i == j else 0.0, f"entry {i},{j}")
    P1 = torch.diag(torch.tensor([1.0] * k1 + [0.0] * (D - k1), dtype=RD)).to(CD)
    R = torch.eye(D, dtype=CD)
    c = s = 1 / math.sqrt(2)
    R[0, 0], R[0, k1], R[k1, 0], R[k1, k1] = c, -s, s, c
    assert abs((R @ P1 @ R.mH)[0, k1].item() - 0.5) < 1e-14


@pytest.mark.parametrize("k", [1, 2, 3, 4])
def test_lemma_5_8_commutant_of_V_tensor_V(k):
    """Lemma 5.8: the commutant of {V (x) V : V in U(k)} is span{1, F} (dimension 2 for k >= 2, 1 for k = 1); F commutes with V (x) V."""
    torch.manual_seed(80 + k)
    Vs = alg_c([k]).random_unitary(k, batch_size=3)
    VV = [torch.kron(V, V) for V in Vs]
    F = swap(k)
    assert all(torch.allclose(F @ W, W @ F, atol=1e-13) for W in VV)
    basis = commutant_basis(VV)
    assert basis.shape[0] == (2 if k >= 2 else 1)
    # 1 and F lie in the span of the basis
    flat = basis.reshape(basis.shape[0], -1)
    for A in (torch.eye(k * k, dtype=CD), F):
        coef = torch.linalg.lstsq(flat.T, A.reshape(-1, 1)).solution
        assert torch.allclose(flat.T @ coef, A.reshape(-1, 1), atol=1e-10)


@pytest.mark.parametrize("k", [2, 3, 4])
def test_corollary_5_9_sym_and_alt_are_inequivalent_irreducibles(k):
    """Corollary 5.9: V -> V (x) V decomposes into two inequivalent irreducibles of dimensions k(k+1)/2 and k(k-1)/2, multiplicity 1."""
    torch.manual_seed(90 + k)
    Vs = alg_c([k]).random_unitary(k, batch_size=3)
    dec = wedderburn([torch.kron(V, V) for V in Vs])
    assert sorted(zip(dec.k_factors, dec.multiplicities)) == sorted([(k * (k + 1) // 2, 1), (k * (k - 1) // 2, 1)])
    F = swap(k)
    for sgn, dim in ((1, k * (k + 1) // 2), (-1, k * (k - 1) // 2)):
        P = (torch.eye(k * k, dtype=CD) + sgn * F) / 2
        assert torch.allclose(P @ P, P, atol=1e-14) and torch.allclose(P, P.mH)
        assert abs(P.trace().real.item() - dim) < 1e-12


@pytest.mark.parametrize("k", [2, 3])
def test_theorem_5_10_twirl_of_order_two(k):
    """Theorem 5.10: E[(U(x)U) X (U(x)U)*] = alpha 1 + beta F with the stated alpha, beta (random X, Monte Carlo)."""
    torch.manual_seed(100 + k)
    X = torch.randn(k * k, k * k, dtype=CD)
    F = swap(k)
    trX, trFX = X.trace(), (F @ X).trace()
    alpha = (trX - trFX / k) / (k ** 2 - 1)
    beta = (trFX - trX / k) / (k ** 2 - 1)
    target = alpha * torch.eye(k * k, dtype=CD) + beta * F
    U = haar(k, 40_000 if k == 2 else 15_000, 101 + k)
    UU = torch.einsum('bij,bkl->bikjl', U, U).reshape(-1, k * k, k * k)
    T = UU @ X @ UU.mH
    for a in range(k * k):
        for b in range(k * k):
            assert_mc_complex(T[:, a, b], target[a, b].item(), f"entry {a},{b}")
    # traces are preserved exactly, as in the proof
    assert torch.allclose(T.diagonal(dim1=-2, dim2=-1).sum(-1), trX.expand(T.shape[0]), atol=1e-11)
    assert torch.allclose((F @ T).diagonal(dim1=-2, dim2=-1).sum(-1), trFX.expand(T.shape[0]), atol=1e-11)


@pytest.mark.parametrize("k", [2, 3, 5])
def test_theorem_5_10_weingarten_values_are_inverse_gram(k):
    """Theorem 5.10: (Wg(s^-1 t)) is the inverse of the Gram matrix [[k^2, k], [k, k^2]], Wg(e) = 1/(k^2-1), Wg(tau) = -1/(k(k^2-1))."""
    _, G, W, wg = weingarten(2, k)
    assert np.allclose(G, [[k * k, k], [k, k * k]])
    assert np.allclose(W @ G, np.eye(2), atol=1e-12)
    assert abs(wg[(0, 1)] - 1 / (k * k - 1)) < 1e-12 and abs(wg[(1, 0)] + 1 / (k * (k * k - 1))) < 1e-12


@pytest.mark.parametrize("k", [2, 3, 4])
def test_corollary_5_11_three_moments(k):
    """Corollary 5.11: E|U11|^4 = 2/(k(k+1)), E|U11|^2|U12|^2 = 1/(k(k+1)), E U11 U22 conj(U12 U21) = -1/((k-1)k(k+1))."""
    U = haar(k, 150_000, 110 + k)
    assert_mc(U[:, 0, 0].abs() ** 4, 2 / (k * (k + 1)), "E|U11|^4")
    assert_mc(U[:, 0, 0].abs() ** 2 * U[:, 0, 1].abs() ** 2, 1 / (k * (k + 1)), "E|U11|^2|U12|^2")
    assert_mc_complex(U[:, 0, 0] * U[:, 1, 1] * (U[:, 0, 1] * U[:, 1, 0]).conj(), -1 / ((k - 1) * k * (k + 1)),
                      "E U11 U22 conj(U12 U21)")
    # the same numbers from the general Weingarten formula
    assert abs(wg_moment((0, 0), (0, 0), (0, 0), (0, 0), k) - 2 / (k * (k + 1))) < 1e-12
    assert abs(wg_moment((0, 0), (0, 1), (0, 0), (0, 1), k) - 1 / (k * (k + 1))) < 1e-12
    assert abs(wg_moment((0, 1), (0, 1), (0, 1), (1, 0), k) + 1 / ((k - 1) * k * (k + 1))) < 1e-12


@pytest.mark.parametrize("n,k", [(2, 2), (3, 2), (3, 3), (3, 4)])
def test_theorem_5_12_weingarten_every_order_monte_carlo(n, k):
    """Theorem 5.12: random moments of order n (pseudo-inverse when n > k) agree with Monte Carlo; Wg is a class function."""
    perms, G, W, wg = weingarten(n, k)
    assert (np.linalg.matrix_rank(G) == len(perms)) == (n <= k)
    for s in perms:                                       # class function
        for t in perms:
            assert abs(wg[compose(compose(t, s), inverse(t))] - wg[s]) < 1e-12
    U = haar(k, 120_000, 120 + 10 * n + k)
    rng = np.random.default_rng(n * 10 + k)
    for _ in range(4):
        i = tuple(rng.integers(0, k, n))
        j = tuple(rng.integers(0, k, n))
        sig, pi = rng.permutation(n), rng.permutation(n)
        ip = tuple(i[sig[m]] for m in range(n))           # balanced index sets give nonzero moments
        jp = tuple(j[pi[m]] for m in range(n))
        prod = torch.ones(U.shape[0], dtype=CD)
        for m in range(n):
            prod = prod * U[:, i[m], j[m]] * U[:, ip[m], jp[m]].conj()
        assert_mc_complex(prod, wg_moment(i, j, ip, jp, k), f"moment {i},{j},{ip},{jp}")


@pytest.mark.parametrize("k", [1, 2, 3])
def test_theorem_5_12_unbalanced_moments_vanish(k):
    """Theorem 5.12: moments with different numbers of U and conj(U) vanish (E U11, E U11^2 conj U11, E U11 U22)."""
    U = haar(k, 60_000, 130 + k)
    assert_mc_complex(U[:, 0, 0], 0.0, "E U11")
    assert_mc_complex(U[:, 0, 0] ** 2 * U[:, 0, 0].conj(), 0.0, "E U11^2 conj U11")
    assert_mc_complex(U[:, 0, 0] * U[:, -1, -1], 0.0, "E U11 Ukk")


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_exercise_5_2_moments_of_one_entry_exact(n):
    """Exercise 5.2: E|U11|^{2n} = 1^T W 1 = n!/(k(k+1)...(k+n-1)) (also n > k by pseudo-inverse), and c_n(x) = x(x+1)...(x+n-1)."""
    for k in range(1, 6):
        perms, G, W, _ = weingarten(n, k)
        assert abs(W.sum() - math.factorial(n) / math.prod(range(k, k + n))) < 1e-10
        assert sum(k ** n_cycles(p) for p in perms) == math.prod(range(k, k + n))


@pytest.mark.parametrize("k,n", [(3, 3), (2, 3), (4, 2)])
def test_exercise_5_2_monte_carlo(k, n):
    """Exercise 5.2 (and Corollary B.14 for n > k): E|U11|^{2n} = n!/(k(k+1)...(k+n-1)) by Monte Carlo."""
    U = haar(k, 150_000, 140 + k + n)
    assert_mc(U[:, 0, 0].abs() ** (2 * n), math.factorial(n) / math.prod(range(k, k + n)), "E|U11|^2n")


@pytest.mark.parametrize("k,n,exact", [(1, 1, 1), (1, 2, 1), (2, 1, 1), (2, 2, 2), (3, 2, 2), (2, 3, 5), (3, 3, 6)])
def test_theorem_5_12_moments_of_the_trace(k, n, exact):
    """Theorem 5.12: E|Tr U|^{2n} = n! (W G)_{ee} (= n! for n <= k; 1 for k = 1; 5 for k = 2, n = 3), exact and by Monte Carlo."""
    perms, G, W, _ = weingarten(n, k)
    e = perms.index(tuple(range(n)))
    assert abs(math.factorial(n) * (W @ G)[e, e] - exact) < 1e-10
    U = haar(k, 120_000, 150 + 10 * k + n)
    assert_mc(U.diagonal(dim1=-2, dim2=-1).sum(-1).abs() ** (2 * n), float(exact), "E|Tr U|^2n")


@pytest.mark.parametrize("k", [2, 3])
def test_theorem_5_12_commutant_of_tensor_cube_is_span_of_permutations(k):
    """Theorem 5.12 / B.6(c): the commutant of {V^(x)3} has dimension rank G = dim span pi(S_3) (5 for k = 2, 6 for k >= 3)."""
    torch.manual_seed(160 + k)
    Vs = alg_c([k]).random_unitary(k, batch_size=3)
    basis = commutant_basis([torch.kron(torch.kron(V, V), V) for V in Vs])
    _, G, _, _ = weingarten(3, k)
    assert basis.shape[0] == np.linalg.matrix_rank(G) == (5 if k == 2 else 6)


@pytest.mark.parametrize("ks,expected4", [([1], 1), ([2], 2), ([1, 1], 6), ([3, 3], 8), ([1, 2, 3], 17)])
def test_exercise_5_1_moments_of_trace_over_U_M(ks, expected4):
    """Exercise 5.1: E|Tr U|^2 = C and E|Tr U|^4 = sum_c m(k_c) + 2C(C-1) for U Haar on U(M)."""
    C = len(ks)
    assert sum(1 if k == 1 else 2 for k in ks) + 2 * C * (C - 1) == expected4
    alg = alg_c(ks)
    torch.manual_seed(170 + sum(ks))
    tr = alg.random_unitary_operator(100_000).trace
    assert_mc(tr.abs() ** 2, float(C), "E|Tr U|^2")
    assert_mc(tr.abs() ** 4, float(expected4), "E|Tr U|^4")


# ================================================================================================================
# Circular ensembles (Exercise 5.3 and the library's CSE)
# ================================================================================================================

@pytest.mark.parametrize("k", [1, 2, 3, 6])
def test_exercise_5_3_coe_structure(k):
    """Exercise 5.3: the library's COE is W^T W with W Haar (same random numbers); it is unitary and symmetric."""
    alg = alg_c([k])
    torch.manual_seed(180)
    S = alg.random_unitary_operator(200, measure='coe').matrix[:, 0]
    torch.manual_seed(180)
    W = alg.random_unitary(k, batch_size=200)
    assert torch.allclose(S, W.mT @ W, atol=1e-13)
    assert torch.allclose(S, S.mT, atol=1e-13)
    assert torch.allclose(S @ S.mH, torch.eye(k, dtype=CD).expand_as(S), atol=1e-12)


@pytest.mark.parametrize("k", [1, 2, 3, 5])
def test_exercise_5_3_coe_moments(k):
    """Exercise 5.3: E Tr S^m = 0 (m = 1, 2, 3, -1); and Theorem 5.10 gives E|Tr S|^2 = 2k/(k+1) for the COE."""
    alg = alg_c([k])
    torch.manual_seed(190 + k)
    S = alg.random_unitary(k, measure='coe', batch_size=60_000)
    for m in (1, 2, 3, -1):
        trm = torch.linalg.matrix_power(S, m).diagonal(dim1=-2, dim2=-1).sum(-1)
        assert_mc_complex(trm, 0.0, f"E Tr S^{m}")
    assert_mc(S.diagonal(dim1=-2, dim2=-1).sum(-1).abs() ** 2, 2 * k / (k + 1), "E|Tr S|^2")


def test_exercise_5_3_coe_invariance_and_uniform_phases():
    """Exercise 5.3: V^T S V has the law of S (moments), and the mean eigenphase density is uniform (E N_B = k|B|/2pi on 4 arcs)."""
    k = 4
    alg = alg_c([k])
    torch.manual_seed(200)
    S = alg.random_unitary(k, measure='coe', batch_size=50_000)
    torch.manual_seed(201)
    V = alg.random_unitary(k) @ torch.diag(torch.exp(1j * torch.tensor([0.3, 1.0, 2.0, -0.4], dtype=RD))).to(CD)
    S2 = V.mT @ S @ V
    for X, name in ((S, "S"), (S2, "V^T S V")):
        assert_mc(X[:, 0, 0].abs() ** 2, 2 / (k + 1), f"{name}: E|S_11|^2")
        assert_mc(X.diagonal(dim1=-2, dim2=-1).sum(-1).abs() ** 2, 2 * k / (k + 1), f"{name}: E|Tr S|^2")
    phases = torch.linalg.eigvals(S).angle()                      # in (-pi, pi]
    for lo in (-math.pi, -math.pi / 2, 0.0, math.pi / 2):
        cnt = ((phases > lo) & (phases <= lo + math.pi / 2)).sum(-1).double()
        assert_mc(cnt, k / 4, f"E N(({lo:.2f}, {lo + math.pi / 2:.2f}])")


def test_exercise_5_3_coe_entry_second_moment_exact():
    """Exercise 5.3 with Theorem 5.10: E|S_11|^2 = E|sum_j W_j1^2|^2 = 2/(k+1) for S = W^T W (k = 1, ..., 5), exact."""
    for k in range(1, 6):
        idx = range(k)
        val = sum(wg_moment((j, j), (0, 0), (l, l), (0, 0), k) for j in idx for l in idx)
        assert abs(val - 2 / (k + 1)) < 1e-12


@pytest.mark.parametrize("n", [2, 4, 6])
def test_library_cse_self_dual_and_kramers(n):
    """Circular symplectic ensemble (library, Dyson): S = W^R W is unitary, self-dual (J S^T J^T = S), with doubly degenerate eigenvalues."""
    alg = alg_c([n])
    torch.manual_seed(210 + n)
    S = alg.random_unitary(n, measure='cse', batch_size=200)
    h = n // 2
    J = torch.zeros(n, n, dtype=CD)
    J[:h, h:] = torch.eye(h, dtype=CD)
    J[h:, :h] = -torch.eye(h, dtype=CD)
    assert torch.allclose(S @ S.mH, torch.eye(n, dtype=CD).expand_as(S), atol=1e-12)
    assert torch.allclose(J @ S.mT @ J.mT, S, atol=1e-12)
    torch.manual_seed(210 + n)
    W = alg.random_unitary(n, batch_size=200)
    assert torch.allclose(S, J @ W.mT @ J.mT @ W, atol=1e-12)
    ph = torch.linalg.eigvals(S).angle().sort(dim=-1).values
    assert (ph[:, 0::2] - ph[:, 1::2]).abs().max() < 1e-6


def test_library_cse_mean_density_uniform():
    """Circular symplectic ensemble: e^{i a} S has the law of S, so E Tr S^m = 0 (m = 1, 2) and the phases are uniform on average."""
    n = 4
    alg = alg_c([n])
    torch.manual_seed(220)
    S = alg.random_unitary(n, measure='cse', batch_size=60_000)
    for m in (1, 2):
        assert_mc_complex(torch.linalg.matrix_power(S, m).diagonal(dim1=-2, dim2=-1).sum(-1), 0.0, f"E Tr S^{m}")
    phases = torch.linalg.eigvals(S).angle()
    cnt = ((phases > 0) & (phases <= 1.0)).sum(-1).double()
    assert_mc(cnt, n / (2 * math.pi), "E #phases in (0, 1]")


# ================================================================================================================
# 5.3 Level statistics
# ================================================================================================================

def test_section_5_3_goe_gue_entry_laws():
    """Section 5.3: GOE/GUE diagonal entries N(0,1), off-diagonal N(0,1/2) (GUE: Re, Im independent N(0,1/2)), uncorrelated."""
    n = 4
    H1, H2 = goe(100_000, n, 230), gue(100_000, n, 231)
    for H in (H1, H2):
        assert_mc(H[:, 1, 1].real ** 2, 1.0, "Var h_ii")
        assert_mc(H[:, 0, 2].real ** 2, 0.5, "Var Re h_ij")
        assert_mc(H[:, 0, 2].real * H[:, 1, 1].real, 0.0, "Cov")
    assert_mc(H2[:, 0, 2].imag ** 2, 0.5, "Var Im h_ij")
    assert_mc(H2[:, 0, 2].real * H2[:, 0, 2].imag, 0.0, "Cov Re Im")


@pytest.mark.parametrize("beta", [1, 2])
def test_proposition_5_13_wigner_surmise(beta):
    """Proposition 5.13: for 2x2 GOE/GUE, E s = sqrt(pi) resp. 4/sqrt(pi), and S = s/Es has the stated distribution functions."""
    if beta == 1:
        alg = alg_r([2])
        torch.manual_seed(240)
        G = torch.randn(300_000, 1, 2, 2, dtype=RD)
        H = alg.operator((G + G.mT) / 2, is_self_adjoint=True)
        mean_s = math.sqrt(math.pi)
        cdf = lambda t: 1 - math.exp(-math.pi * t ** 2 / 4)
    else:
        alg = alg_c([2])
        torch.manual_seed(241)
        G = torch.randn(300_000, 1, 2, 2, dtype=CD)
        H = alg.operator((G + G.mH) / math.sqrt(2), is_self_adjoint=True)
        mean_s = 4 / math.sqrt(math.pi)
        cdf = lambda t: math.erf(2 * t / math.sqrt(math.pi)) - 4 * t / math.pi * math.exp(-4 * t ** 2 / math.pi)
    lam = H.eigenvalues()[0]
    s = lam[:, 1] - lam[:, 0]
    assert_mc(s, mean_s, "E s")
    S = s / mean_s
    for t in (0.1, 0.3, 0.6, 1.0, 1.5, 2.2):
        assert_mc((S <= t).double(), cdf(t), f"P(S <= {t})")
    # the densities integrate to 1 and have mean 1; their CDFs are the stated ones (quadrature)
    p = (lambda x: math.pi / 2 * x * math.exp(-math.pi * x ** 2 / 4)) if beta == 1 else \
        (lambda x: 32 / math.pi ** 2 * x ** 2 * math.exp(-4 * x ** 2 / math.pi))
    assert abs(integrate.quad(p, 0, np.inf)[0] - 1) < 1e-10
    assert abs(integrate.quad(lambda x: x * p(x), 0, np.inf)[0] - 1) < 1e-10
    for t in (0.2, 1.0, 3.0):
        assert abs(integrate.quad(p, 0, t)[0] - cdf(t)) < 1e-10


@pytest.mark.parametrize("beta", [1, 2])
def test_proposition_5_13_small_spacing_exponent(beta):
    """Proposition 5.13: P(S <= eps) = c eps^{beta+1} (1 + O(eps^2)) for the exact distribution functions."""
    cdf = (lambda t: -math.expm1(-math.pi * t ** 2 / 4)) if beta == 1 else \
        (lambda t: math.erf(2 * t / math.sqrt(math.pi)) - 4 * t / math.pi * math.exp(-4 * t ** 2 / math.pi))
    c = math.pi / 4 if beta == 1 else 32 / (3 * math.pi ** 2)
    for eps in (1e-2, 1e-3):
        assert abs(cdf(eps) / (c * eps ** (beta + 1)) - 1) < 2 * eps


def test_proposition_5_14a_independent_levels():
    """Proposition 5.14(a): for independent exponential spacings P(r <= t) = 2t/(1+t) and E r = 2 ln 2 - 1; also for uniform levels (Cor. A.26)."""
    torch.manual_seed(250)
    e = -torch.log(torch.rand(400_000, 2, dtype=RD))
    levels = torch.cat([torch.zeros(400_000, 1, dtype=RD), e.cumsum(-1)], dim=-1)   # 3 levels, 2 spacings
    r = ratio(levels)
    assert_mc(r, 2 * math.log(2) - 1, "E r")
    for t in (0.1, 0.5, 0.9):
        assert_mc((r <= t).double(), 2 * t / (1 + t), f"P(r <= {t})")
    assert abs(integrate.quad(lambda t: 2 * t / (1 + t) ** 2, 0, 1)[0] - (2 * math.log(2) - 1)) < 1e-12
    u = torch.rand(20_000, 30, dtype=RD)
    assert_mc(ratio(u), 2 * math.log(2) - 1, "uniform levels")


@pytest.mark.parametrize("beta", [1, 2])
def test_proposition_5_14b_three_by_three_exact(beta):
    """Proposition 5.14(b): I_0 = 4/27 resp. 2pi/(81 sqrt 3), E r = I_1/I_0 = 4 - 2 sqrt 3 resp. 2 sqrt 3/pi - 1/2, and f(1/x) = x^2 f(x)."""
    f = lambda x: (x + x * x) ** beta / (1 + x + x * x) ** (1 + 1.5 * beta)
    I0 = integrate.quad(f, 0, 1, epsabs=1e-14, epsrel=1e-13)[0]
    I1 = integrate.quad(lambda x: x * f(x), 0, 1, epsabs=1e-14, epsrel=1e-13)[0]
    I0_exact = 4 / 27 if beta == 1 else 2 * math.pi / (81 * math.sqrt(3))
    r_exact = 4 - 2 * math.sqrt(3) if beta == 1 else 2 * math.sqrt(3) / math.pi - 0.5
    assert abs(I0 - I0_exact) < 1e-13 and abs(I1 / I0 - r_exact) < 1e-12
    assert abs(integrate.quad(f, 0, np.inf, epsrel=1e-13)[0] - 2 * I0) < 1e-11
    for x in (0.3, 1.7, 5.0):
        assert abs(f(1 / x) - x * x * f(x)) < 1e-14
    # the recursion 3a L_{a+1} = 3^{1-a} - 1 + (4a - 2) L_a
    L = lambda a: integrate.quad(lambda x: (1 + x + x * x) ** (-a), 0, 1, epsabs=1e-14)[0]
    for a in (0.5, 1.0, 1.5, 2.0, 3.0):
        assert abs(3 * a * L(a + 1) - (3 ** (1 - a) - 1 + (4 * a - 2) * L(a))) < 1e-13
    assert abs(L(1.5) - 2 / 3 * (math.sqrt(3) - 1)) < 1e-14
    assert abs(L(2.5) - (18 * math.sqrt(3) - 22) / 27) < 1e-14
    assert abs(L(1.0) - math.pi / (3 * math.sqrt(3))) < 1e-14


@pytest.mark.parametrize("beta", [1, 2])
def test_proposition_5_14b_three_by_three_monte_carlo(beta):
    """Proposition 5.14(b): Monte Carlo of the 3x3 GOE/GUE spacing ratio: E r and P(r <= t) from the density f_beta/I_0."""
    H = goe(200_000, 3, 260) if beta == 1 else gue(200_000, 3, 261)
    alg = alg_r([3]) if beta == 1 else alg_c([3])
    r = ratio(alg.operator(H.unsqueeze(1), is_self_adjoint=True).eigenvalues()[0])
    f = lambda x: (x + x * x) ** beta / (1 + x + x * x) ** (1 + 1.5 * beta)
    I0 = integrate.quad(f, 0, 1)[0]
    assert_mc(r, 4 - 2 * math.sqrt(3) if beta == 1 else 2 * math.sqrt(3) / math.pi - 0.5, "E r")
    for t in (0.2, 0.5, 0.8):
        assert_mc((r <= t).double(), integrate.quad(f, 0, t)[0] / I0, f"P(r <= {t})")
    # x = s2/s1 on (0, inf) has density f / (2 I0)
    lam = torch.linalg.eigvalsh(H)
    x = (lam[:, 2] - lam[:, 1]) / (lam[:, 1] - lam[:, 0])
    for t in (0.5, 2.0):
        assert_mc((x <= t).double(), integrate.quad(f, 0, t)[0] / (2 * I0), f"P(x <= {t})")


def test_proposition_5_14b_change_of_variables():
    """Proof of Proposition 5.14(b): the map (m, s1, s2) -> lambda has determinant 1, Delta = s1 s2 (s1+s2), |lambda|^2 = 3m^2 + 2/3 (s1^2+s1 s2+s2^2)."""
    M = np.array([[1, -2 / 3, -1 / 3], [1, 1 / 3, -1 / 3], [1, 1 / 3, 2 / 3]])
    assert abs(np.linalg.det(M) - 1) < 1e-14
    rng = np.random.default_rng(0)
    for _ in range(5):
        m, s1, s2 = rng.normal(), rng.random(), rng.random()
        lam = M @ np.array([m, s1, s2])
        assert abs(lam[1] - lam[0] - s1) < 1e-14 and abs(lam[2] - lam[1] - s2) < 1e-14
        assert abs((lam[1] - lam[0]) * (lam[2] - lam[0]) * (lam[2] - lam[1]) - s1 * s2 * (s1 + s2)) < 1e-14
        assert abs(lam @ lam - (3 * m * m + 2 / 3 * (s1 * s1 + s1 * s2 + s2 * s2))) < 1e-13


def test_section_5_3_spacing_ratio_affine_invariance():
    """Definition of the spacing ratio: r_i is unchanged by lambda -> a lambda + b, a > 0 (and level_spacing_ratio sorts)."""
    torch.manual_seed(270)
    lam = torch.randn(50, 20, dtype=RD)
    r = ratio(lam)
    assert torch.allclose(ratio(3.7 * lam - 2.1), r, atol=1e-13)
    assert torch.allclose(ratio(lam[:, torch.randperm(20)]), r, atol=1e-13)


# ================================================================================================================
# 5.4 Why sectors must not be mixed
# ================================================================================================================

def _gue_blocks(ks, B, seed):
    alg = alg_c(ks)
    torch.manual_seed(seed)
    G = torch.randn(B, len(ks), alg.k_max, alg.k_max, dtype=CD)
    return alg.operator((G + G.mH) / math.sqrt(2), is_self_adjoint=True).eigenvalues()


def _pair_integral(rc, rd, eps, grid):
    """int int 1{|x-y| < eps} rho_c(x) rho_d(y) dx dy by quadrature on a fine grid."""
    h = grid[1] - grid[0]
    Fd = np.concatenate([[0.0], integrate.cumulative_trapezoid(rd(grid), grid)])
    inner = np.interp(grid + eps, grid, Fd) - np.interp(grid - eps, grid, Fd)
    return integrate.trapezoid(rc(grid) * inner, grid)


def test_theorem_5_15_no_repulsion_between_sectors():
    """Theorem 5.15: for independent GUE blocks of sizes 1, 2, 3, E N_eps = sum_{c<d} int int 1{|x-y|<eps} rho_c rho_d, and P(delta < eps) >= E N_eps / sum k_c k_d."""
    ks = [1, 2, 3]
    lam = _gue_blocks(ks, 120_000, 280)
    for (c, k) in enumerate(ks):                              # each block is a GUE(k_c): check E Tr H_c^2 = k_c^2
        assert_mc((lam[c] ** 2).sum(-1), float(k * k), f"E Tr H_{c}^2")
    grid = np.linspace(-9, 9, 36_001)
    rho = [gue_density(k) for k in ks]
    for k, r in zip(ks, rho):
        assert abs(integrate.trapezoid(r(grid), grid) - k) < 1e-10
    merged = torch.cat(lam, dim=-1).sort(dim=-1).values
    delta = torch.diff(merged, dim=-1).min(-1).values
    for eps in (0.05, 0.2, 0.5):
        N = sum(((lam[c][:, :, None] - lam[d][:, None, :]).abs() < eps).sum((-2, -1))
                for c in range(3) for d in range(c + 1, 3)).double()
        exact = sum(_pair_integral(rho[c], rho[d], eps, grid) for c in range(3) for d in range(c + 1, 3))
        assert_mc(N, exact, f"E N_{eps}")
        assert (delta < eps).double().mean() >= N.mean() / sum(ks[c] * ks[d] for c in range(3) for d in range(c + 1, 3))
    # the limit E N_eps / (2 eps) -> sum_{c<d} int rho_c rho_d (deterministic)
    limit = sum(integrate.trapezoid(rho[c](grid) * rho[d](grid), grid) for c in range(3) for d in range(c + 1, 3))
    eps = 1e-3
    small = sum(_pair_integral(rho[c], rho[d], eps, grid) for c in range(3) for d in range(c + 1, 3))
    assert abs(small / (2 * eps) / limit - 1) < 1e-4


def test_theorem_5_15_merged_spectrum_has_linear_small_spacing():
    """Theorem 5.15 vs Proposition C.9: two GOE(2) blocks merged: P(delta < eps) ~ eps (ratio 2 on doubling eps), one GOE(2) block: 1 - e^{-eps^2/4}."""
    alg = alg_r([2, 2])
    torch.manual_seed(290)
    G = torch.randn(300_000, 2, 2, 2, dtype=RD)
    lam = alg.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()
    merged = torch.cat(lam, dim=-1).sort(dim=-1).values
    delta = torch.diff(merged, dim=-1).min(-1).values
    p1, p2 = (delta < 0.02).double().mean().item(), (delta < 0.04).double().mean().item()
    assert 1.75 < p2 / p1 < 2.3                                  # exponent 1, not 2
    s = lam[0][:, 1] - lam[0][:, 0]
    for eps in (0.2, 0.5):
        assert_mc((s < eps).double(), -math.expm1(-eps ** 2 / 4), f"P(s < {eps}) in one block")


def test_section_5_4_merged_spectrum_is_union_of_sector_spectra():
    """Section 5.4: the spectrum of H = (+)_c H_c on C^D is the union of the block spectra (eigenvalues() per sector)."""
    ks = [1, 3, 2]
    alg = alg_c(ks)
    torch.manual_seed(300)
    G = torch.randn(10, 3, 3, 3, dtype=CD)
    H = alg.operator(G + G.mH, is_self_adjoint=True)
    merged = torch.cat(H.eigenvalues(), dim=-1).sort(dim=-1).values
    full = torch.zeros(10, 6, 6, dtype=CD)
    off = 0
    for c, k in enumerate(ks):
        full[:, off:off + k, off:off + k] = H.matrix[:, c, :k, :k]
        off += k
    assert torch.allclose(torch.linalg.eigvalsh(full), merged, atol=1e-12)


# ================================================================================================================
# Appendix C.1 - C.2: simple spectra, eigenframes, the Weyl chart
# ================================================================================================================

@pytest.mark.parametrize("n", [2, 3, 5])
@pytest.mark.parametrize("beta", [1, 2])
def test_lemma_C_1_discriminant_and_spectral_projections(n, beta):
    """Lemma C.1: det(Tr H^{i+j-2}) = Delta(lambda)^2, and the Lagrange polynomials P_j(H) are the eigenprojections."""
    H = (goe(5, n, 310) if beta == 1 else gue(5, n, 311)).to(CD)
    lam, V = torch.linalg.eigh(H)
    pw = [torch.linalg.matrix_power(H, m).diagonal(dim1=-2, dim2=-1).sum(-1).real for m in range(2 * n - 1)]
    M = torch.stack([torch.stack([pw[i + j] for j in range(n)], -1) for i in range(n)], -2)
    D2 = vandermonde(lam) ** 2
    assert torch.allclose(torch.linalg.det(M), D2, rtol=1e-8)
    for j in range(n):
        P = torch.eye(n, dtype=CD).expand_as(H).clone()
        for m in range(n):
            if m != j:
                P = P @ (H - lam[:, m, None, None] * torch.eye(n, dtype=CD)) / (lam[:, j] - lam[:, m])[:, None, None]
        Q = V[:, :, j:j + 1] @ V[:, :, j:j + 1].mH
        assert torch.allclose(P, Q, atol=1e-9)


def _weyl_psi(p, n, beta):
    iu = torch.triu_indices(n, n, 1)
    d = beta * n * (n - 1) // 2
    c = p[:d]
    X = torch.zeros(n, n, dtype=CD)
    q = iu.shape[1]
    X[iu[0], iu[1]] = (c[:q] + 1j * c[q:]) if beta == 2 else c[:q].to(CD)
    X = X - X.mH
    E = torch.linalg.matrix_exp(X)
    H = E @ torch.diag(p[d:]).to(CD) @ E.mH
    off = H[iu[0], iu[1]]
    return torch.cat([H.diagonal().real, off.real] + ([off.imag] if beta == 2 else []))


def _weyl_jac(p, n, beta, h=1e-6):
    J = torch.stack([(_weyl_psi(p + h * e, n, beta) - _weyl_psi(p - h * e, n, beta)) / (2 * h)
                     for e in torch.eye(len(p), dtype=RD)], dim=1)
    return torch.linalg.det(J).abs().item()


@pytest.mark.parametrize("n,beta", [(2, 1), (2, 2), (3, 1), (3, 2), (4, 2)])
def test_lemma_C_4_jacobian_of_the_weyl_chart(n, beta):
    """Lemma C.4(a): |det D Psi(X, lambda)| = Delta(lambda)^beta j(X) with j(0) = 1 and j(X) independent of lambda."""
    torch.manual_seed(320 + 10 * n + beta)
    d = beta * n * (n - 1) // 2
    X = 0.3 * torch.randn(d, dtype=RD)
    js = []
    for _ in range(3):
        lam = torch.randn(n, dtype=RD).sort().values
        D = vandermonde(lam).item() ** beta
        assert abs(_weyl_jac(torch.cat([torch.zeros(d, dtype=RD), lam]), n, beta) / D - 1) < 1e-6
        js.append(_weyl_jac(torch.cat([X, lam]), n, beta) / D)
    assert max(js) - min(js) < 1e-6 and min(js) > 0


@pytest.mark.parametrize("beta", [1, 2])
def test_lemma_C_3_theorem_C_7_eigenframe_haar_and_independent(beta):
    """Lemma C.3 / Theorem C.7: the eigenframe of GOE(3)/GUE(3) is Haar on O(3)/U(3) and independent of the eigenvalues."""
    n = 3
    H = goe(150_000, n, 330) if beta == 1 else gue(150_000, n, 331)
    lam, V = torch.linalg.eigh(H)
    u11 = V[:, 0, 0].abs() ** 2
    assert_mc(u11, 1 / 3, "E|U11|^2")
    u4 = 3 / (n * (n + 2)) if beta == 1 else 2 / (n * (n + 1))     # first column uniform on the sphere
    assert_mc(u11 ** 2, u4, "E|U11|^4")
    for j in range(n):
        assert_mc((u11 - 1 / 3) * lam[:, j], 0.0, f"Cov(|U11|^2, lambda_{j})")
        assert_mc((V[:, 1, j].abs() ** 2 - 1 / 3) * (lam[:, 2] - lam[:, 0]), 0.0, "Cov(|U2j|^2, spread)")


# ================================================================================================================
# Appendix C.3 - C.4: Gaussian ensembles, level densities, smallest spacing, orthogonal polynomials
# ================================================================================================================

def test_theorem_C_7_density_normalisation():
    """Theorem C.7: the density (2 pi)^{-n/2} pi^{-beta n(n-1)/4} e^{-Tr H^2/2} integrates to 1 and matches the sampled E Tr H^2 = n + beta n(n-1)/2."""
    for beta in (1, 2):
        for n in (1, 2, 3, 4):
            Z = math.sqrt(2 * math.pi) ** n * math.sqrt(math.pi) ** (beta * n * (n - 1) // 2)
            assert abs(Z * (2 * math.pi) ** (-n / 2) * math.pi ** (-beta * n * (n - 1) / 4) - 1) < 1e-12
    for beta in (1, 2):
        n = 4
        H = goe(100_000, n, 340) if beta == 1 else gue(100_000, n, 341)
        trH2 = (H.abs() ** 2).sum((-2, -1))
        assert_mc(trH2, n + beta * n * (n - 1) / 2, "E Tr H^2 = sum h_ii^2 + 2 sum |h_ij|^2")


def test_theorem_C_7_joint_density_moments():
    """Theorem C.7 (via Corollary C.6): E g(lambda) = int g |Delta|^beta e^{-|x|^2/2} / int |Delta|^beta e^{-|x|^2/2} for GUE(2), GUE(3) (Gauss-Hermite, exact)."""
    xg, wg = np.polynomial.hermite_e.hermegauss(20)                # exact for polynomial x e^{-x^2/2}
    for n in (2, 3):
        grids = np.meshgrid(*([xg] * n), indexing='ij')
        w = np.prod(np.meshgrid(*([wg] * n), indexing='ij'), axis=0)
        D2 = np.ones_like(grids[0])
        for a in range(n):
            for b in range(a + 1, n):
                D2 = D2 * (grids[b] - grids[a]) ** 2
        Zt = (w * D2).sum()
        H = gue(200_000, n, 350 + n)
        lam = torch.linalg.eigvalsh(H)
        for g_np, g_t, name in ((lambda x: sum(xi ** 4 for xi in x), lambda l: (l ** 4).sum(-1), "sum x^4"),
                                (lambda x: np.prod(x, axis=0) ** 2, lambda l: l.prod(-1) ** 2, "prod x^2"),
                                (lambda x: sum(xi ** 2 for xi in x) ** 2, lambda l: (l ** 2).sum(-1) ** 2, "(sum x^2)^2")):
            exact = (w * D2 * g_np(np.array(grids))).sum() / Zt
            assert_mc(g_t(lam), exact, f"GUE({n}) E {name}")
        # E Tr H^4 = 2 n^3 + n (Proposition 1.10) from the joint density
        assert abs((w * D2 * sum(g ** 4 for g in grids)).sum() / Zt - (2 * n ** 3 + n)) < 1e-9


@pytest.mark.parametrize("n", [2, 3])
def test_corollary_C_8_goe_level_density(n):
    """Corollary C.8: the GOE(n) level density rho = n int q(x, y) dy (numerical marginal of the joint density) against bin counts."""
    h = 0.05
    xs = np.arange(-6, 6 + h / 2, h)
    if n == 2:
        X, Y = np.meshgrid(xs, xs, indexing='ij')
        marg = (np.abs(X - Y) * np.exp(-(X ** 2 + Y ** 2) / 2)).sum(1) * h
    else:
        Y, W = np.meshgrid(xs, xs, indexing='ij')
        base = np.abs(W - Y) * np.exp(-(Y ** 2 + W ** 2) / 2)
        marg = np.array([(np.abs((Y - x) * (W - x)) * base).sum() * h * h * math.exp(-x * x / 2) for x in xs])
    rho = n * marg / (marg.sum() * h)
    # closed form for n = 2: rho(x) = e^{-x^2/2} [x (2 Phi(x) - 1) + 2 phi(x)] / sqrt 2
    if n == 2:
        closed = np.exp(-xs ** 2 / 2) * (xs * special.erf(xs / math.sqrt(2)) + 2 * np.exp(-xs ** 2 / 2) / math.sqrt(2 * math.pi)) / math.sqrt(2)
        assert np.abs(rho - closed).max() < 2e-3
    lam = torch.linalg.eigvalsh(goe(200_000, n, 360 + n))
    for lo in (-3.0, -1.5, 0.0, 1.0, 2.5):
        hi = lo + 1.0
        sel = (xs >= lo - 1e-9) & (xs <= hi + 1e-9)
        count = ((lam >= lo) & (lam < hi)).sum(-1).double()
        assert_mc(count, integrate.trapezoid(rho[sel], xs[sel]), f"GOE({n}) bin [{lo}, {hi})")


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_lemma_C_10_gue_one_point_function(n):
    """Lemma C.10: K_n = v^T M^{-1} v = sum p_k^2 (Hermite), int K_n w = n, Z_n = n! prod kappa_k^{-2} (Gauss-Hermite), and GUE bin counts."""
    x = np.linspace(-4, 4, 9)
    M = np.array([[integrate.quad(lambda t: t ** (i + j) * math.exp(-t * t / 2), -np.inf, np.inf)[0]
                   for j in range(n)] for i in range(n)])
    v = np.stack([x ** i for i in range(n)], -1)
    K1 = np.einsum('xi,ij,xj->x', v, np.linalg.inv(M), v)
    K2 = gue_density(n)(x) * np.exp(x ** 2 / 2)
    assert np.allclose(K1, K2, rtol=1e-9)
    assert abs(integrate.quad(gue_density(n), -np.inf, np.inf)[0] - n) < 1e-9
    xg, wg = np.polynomial.hermite_e.hermegauss(12)
    grids = np.meshgrid(*([xg] * n), indexing='ij')
    w = np.prod(np.meshgrid(*([wg] * n), indexing='ij'), axis=0)
    D2 = np.ones_like(grids[0])
    for a in range(n):
        for b in range(a + 1, n):
            D2 = D2 * (grids[b] - grids[a]) ** 2
    Zn = math.factorial(n) * math.prod(math.sqrt(2 * math.pi) * math.factorial(k) for k in range(n))   # kappa_k^{-2}
    assert abs((w * D2).sum() / Zn - 1) < 1e-10
    lam = torch.linalg.eigvalsh(gue(150_000, n, 370 + n))
    for lo in (-3.0, -1.0, 0.5, 2.0):
        exact = integrate.quad(gue_density(n), lo, lo + 1.0)[0]
        assert_mc(((lam >= lo) & (lam < lo + 1.0)).sum(-1).double(), exact, f"GUE({n}) bin {lo}")


def test_lemma_C_10_laguerre_normalisation():
    """Lemma C.10 with w = x^alpha e^{-x}: Z_n = n! prod_k k! (k+alpha)! (Gauss-Laguerre, exact for polynomials)."""
    for al in (0, 1, 3):
        xg, wg = special.roots_genlaguerre(14, al)
        for n in (1, 2, 3):
            grids = np.meshgrid(*([xg] * n), indexing='ij')
            w = np.prod(np.meshgrid(*([wg] * n), indexing='ij'), axis=0)
            D2 = np.ones_like(grids[0])
            for a in range(n):
                for b in range(a + 1, n):
                    D2 = D2 * (grids[b] - grids[a]) ** 2
            Zn = math.factorial(n) * math.prod(math.factorial(k) * math.factorial(k + al) for k in range(n))
            assert abs((w * D2).sum() / Zn - 1) < 1e-9


@pytest.mark.parametrize("beta", [1, 2])
def test_proposition_C_9_smallest_spacing_exponent(beta):
    """Proposition C.9: P(delta < eps) ~ eps^{beta+1} for GOE(3)/GUE(3) (doubling eps multiplies it by about 2^{beta+1})."""
    H = goe(300_000, 3, 380) if beta == 1 else gue(300_000, 3, 381)
    delta = torch.diff(torch.linalg.eigvalsh(H), dim=-1).min(-1).values
    e1, e2 = (0.1, 0.2) if beta == 1 else (0.2, 0.4)
    p1, p2 = (delta < e1).double().mean().item(), (delta < e2).double().mean().item()
    expo = math.log(p2 / p1) / math.log(2)
    assert abs(expo - (beta + 1)) < 0.3


def test_proposition_C_9_explicit_bounds():
    """Proof of Proposition C.9: c eps^{beta+1} <= P(delta < eps) <= C eps^{beta+1} with the constants of the proof (GOE(3), GUE(3), n = 2)."""
    for beta in (1, 2):
        for n in (2, 3):
            # Z_{n,beta} = (1/n!) int |Delta|^beta e^{-|x|^2/2} by Monte Carlo over iid N(0,1) (exact enough for loose bounds)
            torch.manual_seed(390 + n + beta)
            g = torch.randn(200_000, n, dtype=RD).sort(-1).values
            Z = (2 * math.pi) ** (n / 2) * (vandermonde(g).abs() ** beta).mean().item() / math.factorial(n)
            A = lambda p: integrate.quad(lambda x: (1 + abs(x)) ** (p * beta * (n - 1)) * math.exp(-x * x / 2), -np.inf, np.inf)[0]
            C = (n - 1) * 2 ** (beta * (n - 1)) / Z / (beta + 1) * A(2) * A(1) ** (n - 2)
            Mx = 1 + 4 + sum((3 * k + 1) ** 2 for k in range(3, n + 1))
            c = math.exp(-Mx / 2) / Z / (beta + 1)
            H = goe(100_000, n, 395) if beta == 1 else gue(100_000, n, 396)
            delta = torch.diff(torch.linalg.eigvalsh(H), dim=-1).min(-1).values
            for eps in (0.25, 0.5, 1.0):
                p = (delta < eps).double().mean().item()
                assert c * eps ** (beta + 1) <= p <= C * eps ** (beta + 1)


# ================================================================================================================
# Appendix C.5: Bartlett, Cholesky, Wishart, the induced measure
# ================================================================================================================

def test_lemma_C_11_bartlett_decomposition():
    """Lemma C.11: the Cholesky factor of W = G G* has L_ii^2 ~ Gamma(b - i + 1), L_ij standard complex Gaussian, all independent."""
    a, b, B = 3, 5, 150_000
    torch.manual_seed(400)
    G = torch.randn(B, a, b, dtype=CD)
    L = torch.linalg.cholesky(G @ G.mH)
    for i in range(a):                                    # 0-based i: shape b - i
        m = b - i
        x = L[:, i, i].real ** 2
        assert_mc(x, m, f"E L_{i}{i}^2")
        assert_mc(x ** 2, m * (m + 1), f"E L_{i}{i}^4")
        for t in (m / 2, m, 2 * m):
            assert_mc((x <= t).double(), special.gammainc(m, t), f"P(L_{i}{i}^2 <= {t})")
    for (i, j) in ((1, 0), (2, 0), (2, 1)):
        assert_mc(L[:, i, j].abs() ** 2, 1.0, f"E|L_{i}{j}|^2")
        assert_mc_complex(L[:, i, j] ** 2, 0.0, f"E L_{i}{j}^2")
        assert_mc(L[:, i, j].abs() ** 4, 2.0, f"E|L_{i}{j}|^4")
    assert_mc(L[:, 2, 0].abs() ** 2 * L[:, 1, 1].real ** 2, 1.0 * (b - 1), "independence")
    assert_mc(L[:, 2, 1].abs() ** 2 * L[:, 2, 2].real ** 2, 1.0 * (b - 2), "independence (same row)")
    # E det W = prod (b - i + 1) = b!/(b-a)!
    assert_mc(torch.linalg.det(G @ G.mH).real, math.factorial(b) / math.factorial(b - a), "E det W")


@pytest.mark.parametrize("a", [1, 2, 3])
def test_lemma_C_12_cholesky_jacobian(a):
    """Lemma C.12: |det D chi(L)| = 2^a prod_i L_ii^{2(a-i)+1} for chi(L) = L L* (real coordinates, finite differences)."""
    torch.manual_seed(410 + a)
    il = torch.tril_indices(a, a, -1)
    q = il.shape[1]

    def unpack(p):
        L = torch.zeros(a, a, dtype=CD)
        L[range(a), range(a)] = p[:a].to(CD)
        L[il[0], il[1]] = p[a:a + q] + 1j * p[a + q:]
        return L

    def chi(p):
        L = unpack(p)
        W = L @ L.mH
        off = W[il[0], il[1]]
        return torch.cat([W.diagonal().real, off.real, off.imag])

    p = torch.cat([0.5 + torch.rand(a, dtype=RD), torch.randn(2 * q, dtype=RD)])
    h = 1e-6
    J = torch.stack([(chi(p + h * e) - chi(p - h * e)) / (2 * h) for e in torch.eye(len(p), dtype=RD)], dim=1)
    diag = p[:a]
    exact = 2 ** a * math.prod(diag[i].item() ** (2 * (a - 1 - i) + 1) for i in range(a))   # 0-based i
    assert abs(torch.linalg.det(J).abs().item() / exact - 1) < 1e-7


def test_theorem_C_13_wishart_level_density():
    """Theorem C.13 + Lemma C.15: the eigenvalues of W = G G* (a = 3, b = 5) have the level density rho_{a,b} (bins and moments)."""
    a, b = 3, 5
    torch.manual_seed(420)
    G = torch.randn(120_000, a, b, dtype=CD)
    lam = torch.linalg.eigvalsh(G @ G.mH)
    for lo, hi in ((0, 1), (1, 3), (3, 6), (6, 10), (10, 25)):
        exact = integrate.quad(lambda x: rho_laguerre(x, a, b), lo, hi)[0]
        assert_mc(((lam >= lo) & (lam < hi)).sum(-1).double(), exact, f"bin [{lo}, {hi})")
    for p, exact in ((1, a * b), (2, a * b * (a + b))):
        assert abs(integrate.quad(lambda x: x ** p * rho_laguerre(x, a, b), 0, np.inf)[0] - exact) < 1e-8
        assert_mc((lam ** p).sum(-1), float(exact), f"E Tr W^{p}")
    assert abs(integrate.quad(lambda x: rho_laguerre(x, a, b), 0, np.inf)[0] - a) < 1e-10


@pytest.mark.parametrize("a,b", [(1, 3), (2, 2), (2, 5), (3, 4)])
def test_theorem_C_14_library_sampler_is_induced_measure(a, b):
    """Theorem C.14: states.random_density_matrix(rank=b) is G G*/Tr G G* with G an a x b Ginibre matrix (same random numbers)."""
    alg = alg_c([a])
    torch.manual_seed(430)
    rho = states.random_density_matrix(alg, 50, rank=b)
    torch.manual_seed(430)
    G = torch.randn(50, 1, a, b, dtype=CD)
    W = G @ G.mH
    exact = W / W.diagonal(dim1=-2, dim2=-1).sum(-1).real[..., None, None]
    assert torch.allclose(rho.matrix, exact, atol=1e-13)


def test_theorem_C_14_trace_independent_of_spectrum():
    """Theorem C.14: s = Tr W ~ Gamma(ab), independent of the eigenvalues mu of rho_A; for a = 2, mu_1 has density prop. to (1-2t)^2 (t(1-t))^{b-2}."""
    a, b = 2, 3
    torch.manual_seed(440)
    G = torch.randn(200_000, a, b, dtype=CD)
    W = G @ G.mH
    s = W.diagonal(dim1=-2, dim2=-1).sum(-1).real
    mu = torch.linalg.eigvalsh(W / s[:, None, None])
    assert_mc(s, a * b, "E s")
    assert_mc((s - a * b) ** 2, a * b, "Var s")
    for t in (4.0, 6.0, 9.0):
        assert_mc((s <= t).double(), special.gammainc(a * b, t), f"P(s <= {t})")
    m1 = mu[:, 0]
    assert_mc((s - a * b) * (m1 - m1.mean()), 0.0, "Cov(s, mu_1)")
    assert_mc((s - a * b) * ((mu ** 2).sum(-1) - (a + b) / (a * b + 1)), 0.0, "Cov(s, purity)")
    dens = lambda t: (1 - 2 * t) ** 2 * (t * (1 - t)) ** (b - 2)
    Zd = integrate.quad(dens, 0, 0.5)[0]
    for lo, hi in ((0, 0.1), (0.1, 0.25), (0.25, 0.4)):
        assert_mc(((m1 >= lo) & (m1 < hi)).double(), integrate.quad(dens, lo, hi)[0] / Zd, f"P({lo} <= mu_1 < {hi})")
    # the Jacobian of Theta(s, mu) = s (mu_1, ..., mu_{a-1}, 1 - sum mu) is s^{a-1}
    for aa in (2, 3, 4):
        sv, m = 1.7, np.full(aa - 1, 0.1)
        J = np.zeros((aa, aa))
        J[:, 0] = np.append(m, 1 - m.sum())
        for i in range(aa - 1):
            J[i, 1 + i] = sv
            J[aa - 1, 1 + i] = -sv
        assert abs(abs(np.linalg.det(J)) - sv ** (aa - 1)) < 1e-12


# ================================================================================================================
# Appendix C.6: Laguerre polynomials, logarithmic integrals, Page's formula, purity
# ================================================================================================================

def test_lemma_C_15_laguerre_polynomials():
    """Lemma C.15: the explicit L_k^(alpha) equal scipy's, L' = -L_{k-1}^(alpha+1), L^(a) = L^(a+1)_k - L^(a+1)_{k-1}, orthogonality with norm (k+alpha)!/k!."""
    x = np.linspace(0.05, 12, 41)
    for al in range(4):
        for k in range(6):
            assert np.allclose(laguerre(k, al, x), special.eval_genlaguerre(k, al, x), rtol=1e-11, atol=1e-11)
            h = 1e-5
            der = (laguerre(k, al, x + h) - laguerre(k, al, x - h)) / (2 * h)
            prev = laguerre(k - 1, al + 1, x) if k >= 1 else 0 * x
            assert np.allclose(der, -prev, atol=1e-6 * (1 + np.abs(prev)))
            assert np.allclose(laguerre(k, al, x), laguerre(k, al + 1, x) - prev, atol=1e-9)
        xg, wg = special.roots_genlaguerre(12, al)
        for k in range(6):
            for l in range(6):
                v = (wg * laguerre(k, al, xg) * laguerre(l, al, xg)).sum()
                ex = math.factorial(k + al) / math.factorial(k) if k == l else 0.0
                assert abs(v - ex) < 1e-8 * max(1.0, ex)


def test_lemma_C_16_logarithmic_integrals():
    """Lemma C.16: I_m = m! (I_0 + H_m) with I_0 = -gamma, and A_{k,m}^(alpha) = I_{k+alpha}/k! (m = k), -(k+alpha)!/(k! (m-k)) (m > k)."""
    def logint(g):
        return integrate.quad(lambda t: math.log(t) * g(t), 0, 1, limit=200)[0] + \
            integrate.quad(lambda t: math.log(t) * g(t), 1, np.inf, limit=200)[0]
    I0 = logint(lambda t: math.exp(-t))
    assert abs(I0 + np.euler_gamma) < 1e-10
    Im = lambda m: math.factorial(m) * (I0 + float(harmonic(m)))
    for m in range(7):
        assert abs(logint(lambda t: t ** m * math.exp(-t)) - Im(m)) < 1e-8 * math.factorial(m)
    for al in range(3):
        for k in range(4):
            for m in range(k, 5):
                A = logint(lambda t: laguerre(k, al, t) * laguerre(m, al, t) * t ** al * math.exp(-t))
                ex = Im(k + al) / math.factorial(k) if m == k else \
                    -math.factorial(k + al) / (math.factorial(k) * (m - k))
                assert abs(A - ex) < 1e-8 * max(1.0, abs(ex))


def test_theorem_C_17_sigma_identity_and_log_integral():
    """Proof of Theorem C.17: Sigma = ab H_b + a(a-1)/2 exactly, and int x log x rho_{a,b} = ab I_0 + Sigma by quadrature."""
    for a in range(1, 6):
        for b in range(a, 9):
            al = b - a
            Sig = sum((k + al + 1) * harmonic(k + al + 1) + k * harmonic(k + al) + 2 * k for k in range(a))
            assert Sig == a * b * harmonic(b) + Fraction(a * (a - 1), 2)
            assert sum(2 * k + al + 1 for k in range(a)) == a * b
    for a, b in ((1, 1), (2, 2), (2, 5), (3, 4)):
        al = b - a
        Sig = float(sum((k + al + 1) * harmonic(k + al + 1) + k * harmonic(k + al) + 2 * k for k in range(a)))
        J = integrate.quad(lambda t: t * math.log(t) * rho_laguerre(t, a, b), 0, 1)[0] + \
            integrate.quad(lambda t: t * math.log(t) * rho_laguerre(t, a, b), 1, np.inf, limit=200)[0]
        assert abs(J - (-a * b * np.euler_gamma + Sig)) < 1e-8
        # Page's formula from the reduction E S = I_0 + H_ab - J/(ab)
        assert abs(-np.euler_gamma + float(harmonic(a * b)) - J / (a * b) - float(page(a, b))) < 1e-9


def test_theorem_C_17_page_formula_a2_by_quadrature():
    """Theorem C.17 for a = 2: integrating S(rho_A) against the density of Theorem C.14 gives Page's value (1/3, 9/20, 107/210, ...)."""
    assert [page(2, b) for b in (2, 3, 4)] == [Fraction(1, 3), Fraction(9, 20), Fraction(107, 210)]
    for b in range(2, 8):
        dens = lambda t: (1 - 2 * t) ** 2 * (t * (1 - t)) ** (b - 2)
        S = lambda t: -(t * math.log(t) + (1 - t) * math.log(1 - t))
        val = integrate.quad(lambda t: dens(t) * S(t), 0, 0.5, epsabs=1e-14)[0] / integrate.quad(dens, 0, 0.5)[0]
        assert abs(val - float(page(2, b))) < 1e-10


@pytest.mark.parametrize("a,b", [(1, 4), (2, 2), (2, 3), (3, 3), (3, 5), (4, 4)])
def test_theorem_C_17_page_formula_monte_carlo(a, b):
    """Theorem C.17: E S(rho_A) = sum_{k=b+1}^{ab} 1/k - (a-1)/(2b) for the library's induced measure (von_neumann_entropy)."""
    torch.manual_seed(450 + 10 * a + b)
    rho = states.random_density_matrix(alg_c([a]), 60_000, rank=b)
    assert_mc(rho.von_neumann_entropy(), float(page(a, b)), f"E S, a={a}, b={b}")


@pytest.mark.parametrize("a,b", [(1, 3), (2, 2), (2, 3), (3, 4), (4, 4)])
def test_purity_of_random_states(a, b):
    """Theorem C.14 + Lemma C.15 (Theorem 7.11): E Tr rho_A^2 = E Tr W^2 / E s^2 = ab(a+b)/(ab(ab+1)) = (a+b)/(ab+1)."""
    m2 = integrate.quad(lambda x: x * x * rho_laguerre(x, a, b), 0, np.inf)[0]
    assert abs(m2 / (a * b * (a * b + 1)) - (a + b) / (a * b + 1)) < 1e-10
    torch.manual_seed(460 + 10 * a + b)
    rho = states.random_density_matrix(alg_c([a]), 60_000, rank=b)
    assert_mc(rho.purity(), (a + b) / (a * b + 1), f"E purity, a={a}, b={b}")


def test_corollary_C_18_bounds():
    """Corollary C.18: log a - (a-1)/(ab+1) - (a-1)/(2b) <= E S <= log a - (a-1)/(2b), |E S - log a + a/(2b)| <= 1/(2b); two qubits 1/3."""
    assert page(2, 2) == Fraction(1, 3)
    for a in range(1, 12):
        for b in range(a, 40):
            p = float(page(a, b))
            assert math.log(a) - (a - 1) / (a * b + 1) - (a - 1) / (2 * b) - 1e-12 <= p <= math.log(a) - (a - 1) / (2 * b) + 1e-12
            assert abs(p - math.log(a) + a / (2 * b)) <= 1 / (2 * b) + 1e-12
            assert math.log((a * b + 1) / (b + 1)) - 1e-12 <= float(harmonic(a * b) - harmonic(b)) <= math.log(a) + 1e-12
