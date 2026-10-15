"""
Stability of the numerical Wedderburn decomposition (lecture notes, Section 4.9): Lemma 4.24 (perturbation of singular
values, spectral subspaces, orthonormal bases and polar factors), Theorem 4.25 (the commutant at a tolerance),
Proposition 4.26 (a-posteriori certificate, ``Decomposition.certificate``), Proposition 4.27 (random elements of the
commutant) and Theorem 4.28 (stability of the decomposition), checked against dense computations that do not use the
library, and the noise-scaled thresholds of ``decompose.wedderburn``.
"""
import math
import warnings

import pytest
import torch

import torch_vn_algebra.decompose as dc

C128 = torch.complex128


# ----------------------------------------------------------------------------
# independent helpers (dense linear algebra, no library internals)
# ----------------------------------------------------------------------------

def gen(seed):
    return torch.Generator().manual_seed(seed)


def cgauss(g, *shape):
    return torch.randn(*shape, dtype=C128, generator=g)


def haar(n, g):
    Q, R = torch.linalg.qr(cgauss(g, n, n))
    return Q * (torch.diagonal(R) / torch.diagonal(R).abs())


def hidden(ks, ms, g, J=2):
    """J random Hermitian generators of V((+)_c M_{k_c} (x) 1_{m_c})V^*, unit Frobenius norm, and V."""
    n = sum(k * m for k, m in zip(ks, ms))
    V = haar(n, g)
    Hs = []
    for _ in range(J):
        blocks = []
        for k, m in zip(ks, ms):
            a = cgauss(g, k, k)
            blocks.append(torch.kron(a + a.mH, torch.eye(m, dtype=C128)))
        Hs.append(V @ torch.block_diag(*blocks) @ V.mH)
    H = torch.stack(Hs)
    return H / torch.linalg.matrix_norm(H)[:, None, None], V


def comm_matrix(Hs):
    """Matrix of X -> ([H_j, X])_j on row-major vec(X)."""
    n = Hs.shape[-1]
    eye = torch.eye(n, dtype=C128)
    return torch.cat([torch.kron(H, eye) - torch.kron(eye, H.T.contiguous()) for H in Hs])


def svals_up(M):
    return torch.linalg.svdvals(M).flip(0)


def null_projector(L, d):
    """Projection (on vec) onto the right singular vectors of the d smallest singular values, and their rows."""
    Vh = torch.linalg.svd(L)[2]
    F = Vh[-d:].conj()
    return F.T @ F.conj(), F


def opnorm(X):
    return torch.linalg.matrix_norm(X, ord=2).item()


def fro(X):
    return torch.linalg.matrix_norm(X).item()


def noisy(H, eta, g):
    """H + E_j with Hermitian E_j of operator norm eta / sqrt(J): (sum_j ||E_j||^2)^{1/2} = eta."""
    E = cgauss(g, *H.shape)
    E = E + E.mH
    E = E / torch.linalg.matrix_norm(E, ord=2)[:, None, None] * eta / math.sqrt(H.shape[0])
    return H + E


def unit(S):
    return S / torch.linalg.matrix_norm(S)[:, None, None]


def project_onto_B(D, H):
    """E_B(H) for B = V((+)_c M_{k_c} (x) 1_{m_c})V^*, computed directly from V (average of the m_c copies)."""
    M = D.V.mH @ H @ D.V
    out = torch.zeros_like(M)
    o = 0
    for k, m in zip(D.k_factors, D.multiplicities):
        blk = M[..., o:o + k * m, o:o + k * m].reshape(*M.shape[:-2], k, m, k, m)
        a = torch.einsum('...ijkj->...ik', blk) / m
        out[..., o:o + k * m, o:o + k * m] = torch.kron(a, torch.eye(m, dtype=C128)) if a.dim() == 2 else \
            torch.stack([torch.kron(x, torch.eye(m, dtype=C128)) for x in a])
        o += k * m
    return D.V @ out @ D.V.mH


def quiet(f, *a, **k):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', dc.DecompositionWarning)
        return f(*a, **k)


# ----------------------------------------------------------------------------
# Lemma 4.24
# ----------------------------------------------------------------------------

def test_lemma_4_24_ab_singular_values_and_spectral_subspaces():
    """Lemma 4.24(a),(b): Weyl for singular values and the 2 Delta / gap bounds for spectral projections."""
    g = gen(0)
    for trial in range(20):
        p, q, d = 7, 9, 3
        T = cgauss(g, q, p)
        U, s, Vh = torch.linalg.svd(T, full_matrices=False)
        s[-d:] = 0                                                      # kernel of dimension d
        T = U @ torch.diag(s).to(C128) @ Vh
        gamma = svals_up(T)[d].item()
        F = cgauss(g, q, p)
        Delta = gamma * 10 ** (-0.5 - 2 * torch.rand(1, generator=g).item())
        Tt = T + Delta * F / opnorm(F)
        assert (svals_up(Tt) - svals_up(T)).abs().max() <= Delta * (1 + 1e-10)
        if Delta < gamma / 2:
            P, _ = null_projector(T, d)
            Pt, _ = null_projector(Tt, d)
            assert opnorm(Pt - P) <= 2 * Delta / gamma
        # (b): an eigenvalue mu of multiplicity 2 of a Hermitian matrix
        Q = haar(6, g)
        lam = torch.tensor([0.0, 0.0, 1.0, -1.3, 2.0, 0.7], dtype=torch.float64)
        A = Q @ torch.diag(lam).to(C128) @ Q.mH
        P = Q[:, :2] @ Q[:, :2].mH
        gap = 0.7
        Eh = cgauss(g, 6, 6)
        Eh = Eh + Eh.mH
        Dl = gap * 10 ** (-0.4 - 2 * torch.rand(1, generator=g).item())
        At = A + Dl * Eh / opnorm(Eh)
        w, Qt = torch.linalg.eigh(At)
        inside = (w.abs() <= Dl + 1e-12)
        assert int(inside.sum()) == 2 and (w[~inside].abs() >= gap - Dl - 1e-12).all()
        Pt = Qt[:, inside] @ Qt[:, inside].mH
        assert opnorm(Pt - P) <= 2 * Dl / gap


def test_lemma_4_24_cd_bases_and_polar_factors():
    """Lemma 4.24(c),(d): the nearest orthonormal basis of ran P and the perturbation of polar factors."""
    g = gen(1)
    for trial in range(30):
        n, k = 7, 3
        Q0 = haar(n, g)[:, :k]
        P = Q0 @ Q0.mH
        X = cgauss(g, n, n)
        X = X + X.mH
        Pt_full = torch.linalg.matrix_exp(1j * 0.1 * torch.rand(1, generator=g).item() * X)
        Qt = Pt_full @ Q0                                            # a rotated subspace
        Pt = Qt @ Qt.mH
        eps = opnorm(Pt - P)
        if eps > 0.5:
            continue
        S = Qt.mH @ P @ Qt
        w, U = torch.linalg.eigh(S)
        assert w.min() >= 1 - eps ** 2 - 1e-12
        Q = P @ Qt @ (U @ torch.diag(w ** -0.5).to(C128) @ U.mH)
        assert torch.allclose(Q.mH @ Q, torch.eye(k, dtype=C128), atol=1e-12)
        assert torch.allclose(P @ Q, Q, atol=1e-12)
        assert opnorm(Q - Qt) <= 2 * eps + 1e-12
        # (d)
        Uk = haar(k, g)
        lam = 0.5 + torch.rand(1, generator=g).item()
        Fk = cgauss(g, k, k)
        zeta = lam * torch.rand(1, generator=g).item() * 0.99
        M = lam * Uk + zeta * Fk / opnorm(Fk)
        A_, _, Bh = torch.linalg.svd(M)
        pol = A_ @ Bh
        w2, V2 = torch.linalg.eigh(M.mH @ M)
        pol2 = M @ V2 @ torch.diag(w2 ** -0.5).to(C128) @ V2.mH            # M (M^* M)^{-1/2}
        assert torch.allclose(pol, pol2, atol=1e-10)
        assert opnorm(pol - Uk) <= 2 * zeta / lam + 1e-12


# ----------------------------------------------------------------------------
# Theorem 4.25
# ----------------------------------------------------------------------------

def test_theorem_4_25_bounds_on_random_perturbations():
    """Theorem 4.25(a),(b): singular values within 2 eta; in the window dim d, *-closed, within 4 eta / gamma."""
    g = gen(2)
    cases = [((1, 2), (2, 1)), ((2, 1), (1, 2)), ((1, 1, 2), (2, 1, 1)), ((2,), (2,)), ((1, 3), (3, 1))]
    for trial in range(25):
        ks, ms = cases[trial % len(cases)]
        H, _ = hidden(ks, ms, g)
        n, d = H.shape[-1], sum(m * m for m in ms)
        L = comm_matrix(H)
        s = svals_up(L)
        gamma = s[d].item()
        assert s[d - 1] < 1e-12
        eta = gamma * 10 ** (-0.7 - 4 * torch.rand(1, generator=g).item())
        Ht = noisy(H, eta, g)
        eta_true = math.sqrt(sum(opnorm(x) ** 2 for x in Ht - H))
        assert abs(eta_true - eta) < 1e-12 * max(1, eta)
        st = svals_up(comm_matrix(Ht))
        assert (st - s).abs().max() <= 2 * eta * (1 + 1e-9)
        if eta < gamma / 4:
            theta = math.sqrt(2 * eta * (gamma - 2 * eta))           # inside (2 eta, gamma - 2 eta)
            assert int((st <= theta).sum()) == d
            P, _ = null_projector(L, d)
            Pt, F = null_projector(comm_matrix(Ht), d)
            assert opnorm(Pt - P) <= 4 * eta / gamma
            Xs = F.reshape(d, n, n).mH.reshape(d, -1)                 # adjoints of the basis of C_theta
            assert torch.allclose(Pt @ Xs.T, Xs.T, atol=1e-9)


def test_theorem_4_25_example_attains_bound():
    """Theorem 4.25(c): the example diag(1,1,0) + eta(E_12 + E_21) attains 2 eta; dimensions 3, 5, 9."""
    for eta in (1e-6, 1e-3, 0.05):
        H = torch.diag(torch.tensor([1.0, 1.0, 0.0], dtype=C128))
        Ht = H.clone()
        Ht[0, 1] = Ht[1, 0] = eta
        st = svals_up(comm_matrix(Ht[None]))
        expected = torch.tensor([0, 0, 0, 2 * eta, 2 * eta, 1 - eta, 1 - eta, 1 + eta, 1 + eta], dtype=torch.float64)
        assert torch.allclose(st, expected, atol=1e-13)
        s = svals_up(comm_matrix(H[None]))
        assert torch.allclose(s, torch.tensor([0.0] * 5 + [1.0] * 4, dtype=torch.float64), atol=1e-14)
        # the library: below 2 eta the symmetry is lost, inside the window it is found, above the gap all is "symmetric"
        scale = fro(Ht)
        assert quiet(dc.commutant_basis, Ht, tol=0.5 * eta / scale).shape[0] == 3
        assert quiet(dc.commutant_basis, Ht, tol=math.sqrt(2 * eta * (1 - 2 * eta)) / scale).shape[0] == 5
        assert quiet(dc.commutant_basis, Ht, tol=(1 + 3 * eta) / scale).shape[0] == 9


@pytest.mark.parametrize('name', ['spin3', 'perm3', 'hidden'])
def test_theorem_4_25_library_recovers_structure_in_the_window(name):
    """Theorems 4.25, 4.28 in the library: tol in the window recovers (k_c, m_c), certified, eta = 1e-9..1e-3."""
    g = gen(3)
    if name == 'spin3':
        S, km = dc.total_spin(3), ((2, 4), (2, 1))
    elif name == 'perm3':
        S, km = dc.permutation_representation(3, 2, 'adjacent').to(C128), ((1, 2), (4, 2))
    else:
        S, _ = hidden((1, 2, 2), (2, 1, 2), g, J=3)
        km = ((1, 2, 2), (2, 1, 2))
    H = unit(S)
    d = sum(m * m for m in km[1])
    gamma = svals_up(comm_matrix(H))[d].item()
    for eta in (1e-9, 1e-7, 1e-5, 1e-3):
        St = noisy(S, eta * float(torch.linalg.matrix_norm(S).min()), g)
        Ht = unit(St)
        eta_n = math.sqrt(sum(fro(x) ** 2 for x in Ht - H))
        st = svals_up(comm_matrix(Ht))
        assert st[d - 1] <= 2 * eta_n + 1e-15 and st[d] >= gamma - 2 * eta_n
        tol = math.sqrt(max(st[d - 1].item(), 1e-16) * st[d].item()) / max(1.0, st[-1].item())
        D = quiet(dc.from_generators, St, tol=tol)
        assert (D.k_factors, D.multiplicities) == km
        c = D.certificate(warn=False)
        assert c.certified and c.trusted and c.ok
        assert c.residuals['backward'] <= 10 * eta_n / gamma + 1e-12   # Theorem 4.28(c): V within O(eta / gamma)
        assert abs(c.gap - gamma) <= 6 * eta_n + 1e-9              # Lemma 4.24(a): ||H_B - H|| <= eta + r <= 3 eta


def test_noise_scaled_thresholds_regression():
    """Regression: noise-scaled thresholds min(1e3 nu, sqrt nu) (were 1e3 nu ~ 1) recover noisy structure."""
    g = gen(4)
    S = dc.total_spin(3)
    E = cgauss(g, *S.shape)
    St = S + 1e-3 * (E + E.mH) / torch.linalg.matrix_norm(E + E.mH, ord=2)[:, None, None]
    D = quiet(dc.from_generators, St, tol=1e-2)
    assert (D.k_factors, D.multiplicities) == ((2, 4), (2, 1))
    th = D.diagnostics['thresholds']
    ci = D.diagnostics['commutant_info']
    nu = max(ci['sigma_null_max'], dc._eps(C128) * max(1.0, ci['sigma_max'])) / ci['sigma_kept_min']
    assert th['small'] == pytest.approx(max(math.sqrt(dc._eps(C128)), min(1e3 * nu, math.sqrt(nu))))
    assert th['small'] < 0.1
    # clean input: unchanged thresholds (sqrt(eps))
    D0 = dc.from_generators(S)
    assert D0.diagnostics['thresholds']['small'] == pytest.approx(math.sqrt(dc._eps(C128)))


# ----------------------------------------------------------------------------
# Proposition 4.26 and Decomposition.certificate
# ----------------------------------------------------------------------------

def test_proposition_4_26_gap_formula_and_generation():
    """Proposition 4.26(a)-(c): block formula for gamma_B, B generated by E(H~_j), certified window."""
    g = gen(5)
    for ks, ms in [((1, 2), (2, 1)), ((2, 3), (1, 2)), ((2, 2), (1, 1)), ((1, 1, 2), (3, 1, 1))]:
        S, _ = hidden(ks, ms, g, J=2)
        St = noisy(S, 1e-5, g)
        D = quiet(dc.from_generators, St, tol=1e-3)
        c = D.certificate(warn=False)
        d = D.commutant_dim
        Ht = unit(St)
        Hb = project_onto_B(D, Ht)
        sb = svals_up(comm_matrix(Hb))
        assert sb[d - 1] < 1e-10 and abs(sb[d].item() - c.gap) < 1e-9
        r = math.sqrt(sum(fro(x) ** 2 for x in Ht - Hb))
        assert abs(r - c.residuals['backward']) < 1e-12
        st = svals_up(comm_matrix(Ht))
        for t in torch.linspace(0.02, 0.98, 7).tolist():
            theta = math.exp((1 - t) * math.log(c.window[0]) + t * math.log(c.window[1]))
            assert int((st <= theta).sum()) == d
        # no model with more symmetry within gamma_B / 2 - r: the exact generators are farther away or have dim <= d
        assert D.commutant_dim == sum(m * m for m in ms)


def test_certificate_values_against_independent_computation():
    """Decomposition.certificate(): values against direct computations, tol_window = window / scale."""
    S = dc.total_spin(3)
    g = gen(6)
    St = S + 1e-5 * (lambda E: E + E.mH)(cgauss(g, *S.shape))
    D = quiet(dc.from_generators, St, tol=1e-3)
    c = D.certificate(warn=False)
    eye = torch.eye(D.n, dtype=C128)
    assert c.residuals['unitarity'] == pytest.approx(opnorm(D.V.mH @ D.V - eye), abs=1e-15)
    assert c.residuals['unitarity'] < 1e-13
    res = D.verify()
    assert c.residuals['algebra'] == res['algebra'] and c.ok == res['ok'] and c.commutant_dim == (5, 5)
    Ht = unit(St)
    r = math.sqrt(sum(fro(x) ** 2 for x in Ht - project_onto_B(D, Ht)))
    assert c.residuals['backward'] == pytest.approx(r, rel=1e-8)
    ci = D.diagnostics['commutant_info']
    assert c.rank['threshold'] == ci['threshold'] and c.rank['sigma_null_max'] == ci['sigma_null_max']
    assert c.rank['gap'] == pytest.approx(ci['sigma_kept_min'] / ci['sigma_null_max'])
    scale = max(1.0, ci['sigma_max'])
    assert c.tol_window == pytest.approx((c.window[0] / scale, c.window[1] / scale))
    lo, hi = c.tol_window
    dim = lambda tol: quiet(dc.commutant_basis, St, tol=tol).shape[0]
    assert all(dim(lo ** (1 - t) * hi ** t) == 5 for t in (0.1, 0.3, 0.5, 0.7, 0.9))
    st = svals_up(comm_matrix(Ht))                                          # full map: Proposition 4.26(c)
    assert st[4] < c.window[0] and st[5] > c.window[1]
    assert c.links['tree_min'] > c.links['threshold'] > c.links['cross_max']
    assert 'Certificate(' in repr(c)


def test_certificate_flags_untrusted_cut():
    """Decomposition.certificate(): warning and trusted=False below 2r or for foreign generators."""
    S = dc.total_spin(3)
    g = gen(7)
    E = cgauss(g, *S.shape)
    St = S + 1e-4 * (E + E.mH)
    D = quiet(dc.from_generators, St, tol=1e-2)
    r, scale = D.certificate(warn=False).residuals['backward'], max(1.0, D.diagnostics['commutant_info']['sigma_max'])
    D_mid = quiet(dc.from_generators, St, tol=1.5 * r / scale)          # above the noise, below 2 r
    assert (D_mid.k_factors, D_mid.multiplicities) == ((2, 4), (2, 1))
    with pytest.warns(dc.DecompositionWarning, match='not certified'):
        c = D_mid.certificate()
    assert c.certified and not c.trusted and c.rank['threshold'] < c.window[0]
    with warnings.catch_warnings():
        warnings.simplefilter('error', dc.DecompositionWarning)
        c0 = dc.from_generators(S).certificate()
        cg = dc.from_group(dc.permutation_representation(3, 2)).certificate()
    assert c0.trusted and c0.residuals['backward'] < 1e-12
    assert cg.trusted and cg.rank is None and cg.tol_window is None               # twirl: no rank cut
    # a decomposition that does not decompose the generators (sigma_z on the first spin breaks SU(2)) is not certified
    D_wrong = dc.from_generators(dc.total_spin(2))
    with pytest.warns(dc.DecompositionWarning):
        sz1 = torch.diag(torch.tensor([1.0, 1, -1, -1], dtype=C128))
        cw = D_wrong.certificate(generators=dc.total_spin(2) + 0.3 * sz1)
    assert not cw.trusted


def test_certificate_given_commutant_and_backward_compatibility():
    """Decomposition.certificate() with a given commutant; verify() unchanged (backward compatibility)."""
    S, _ = hidden((1, 2), (2, 2), gen(8))
    C = dc.commutant_basis(S)
    D = dc.from_generators(S, commutant=C)
    c = D.certificate(warn=False)
    assert c.rank is None and c.trusted and c.residuals['backward'] < 1e-12
    assert set(D.verify()) == {'unitarity', 'algebra', 'commutant', 'commutant_dim', 'ok', 'tol'}
    assert 'Certificate' in dc.__all__


# ----------------------------------------------------------------------------
# Proposition 4.27
# ----------------------------------------------------------------------------

def _gaussian_blocks(D, nsamp, g):
    """Blocks kappa_c of standard Gaussian vectors of A'_h and z_c of standard complex Gaussian vectors of A'."""
    B = D.basis('commutant')
    Z = torch.einsum('br,rij->bij', cgauss(g, nsamp, B.shape[0]), B)
    Z2 = torch.einsum('br,rij->bij', cgauss(g, nsamp, B.shape[0]), B)
    G = (Z2 + Z2.mH) / math.sqrt(2)
    return D.to_commutant(G).matrix, D.to_commutant(Z).matrix, G, Z


@pytest.mark.parametrize('ks,ms', [((1, 2), (2, 1)), ((1, 2, 3), (3, 2, 1)), ((2, 1), (2, 2))])
def test_proposition_4_27_cross_sector_gap_bound(ks, ms):
    """Proposition 4.27(a),(b),(d) by Monte Carlo: block laws, cross-sector gap bound, relative bound."""
    g = gen(9)
    S, _ = hidden(ks, ms, g)
    D = dc.from_generators(S)
    order = sorted(range(len(ks)), key=lambda c: (ks[c], ms[c]))
    assert (D.k_factors, D.multiplicities) == (tuple(ks[c] for c in order), tuple(ms[c] for c in order))
    ks, ms = D.k_factors, D.multiplicities
    nsamp = 40000
    kap, _, G, _ = _gaussian_blocks(D, nsamp, g)
    d = D.commutant_dim
    ev = [torch.linalg.eigvalsh(kap[:, c, :m, :m]) for c, m in enumerate(ms)]
    for c, (k, m) in enumerate(zip(ks, ms)):
        assert (kap[:, c, :m, :m].abs() ** 2).sum((-2, -1)).mean().item() == pytest.approx(m * m / k, rel=0.03)
    assert (G.abs() ** 2).sum((-2, -1)).mean().item() == pytest.approx(d, rel=0.03)
    cross = torch.full((nsamp,), float('inf'), dtype=torch.float64)
    within = torch.full((nsamp,), float('inf'), dtype=torch.float64)
    Ssum = 0.0
    for c in range(len(ks)):
        if ms[c] > 1:
            within = torch.minimum(within, (ev[c][:, 1:] - ev[c][:, :-1]).min(1).values)
        for e in range(c + 1, len(ks)):
            cross = torch.minimum(cross, (ev[c][:, :, None] - ev[e][:, None, :]).abs().flatten(1).min(1).values)
            Ssum += ms[c] * ms[e] / math.sqrt(1 / (ks[c] * ms[c]) + 1 / (ks[e] * ms[e]))
    norm = torch.linalg.matrix_norm(G)
    for x in (0.003, 0.01, 0.03):
        assert (cross < x).double().mean().item() <= math.sqrt(2 / math.pi) * x * Ssum
    for eps in (1e-3, 3e-3):
        p = (torch.minimum(cross, within) < eps * norm).double().mean().item()
        assert p <= 4 * math.sqrt(d / math.pi) * eps * Ssum                     # within-sector term is O(eps^3)
        assert p > 0.05 * 4 * math.sqrt(d / math.pi) * eps * Ssum               # and the order eps is attained


def test_proposition_4_27_within_sector_and_intertwiners():
    """Proposition 4.27(c),(e) by Monte Carlo: 2 x 2 spacing law, Exp(1) law of intertwiner blocks."""
    g = gen(10)
    S, _ = hidden((3, 1), (2, 1), g)
    D = dc.from_generators(S)
    assert D.k_factors == (1, 3) and D.multiplicities == (1, 2)
    nsamp = 40000
    kap, zb, _, _ = _gaussian_blocks(D, nsamp, g)
    k = 3
    w, U = torch.linalg.eigh(kap[:, 1, :2, :2])
    spacing = w[:, 1] - w[:, 0]
    for x in (0.3, 0.5, 0.8):
        y = k ** 0.5 * x                                                    # spacing of the GUE block H_c
        exact = math.erf(y / 2) - y / math.sqrt(math.pi) * math.exp(-y * y / 4)   # Proposition 5.13
        assert exact <= y ** 3 / (6 * math.sqrt(math.pi))
        assert (spacing < x).double().mean().item() == pytest.approx(exact, rel=0.08)
    # ||Pi Z Pi'||_2^2 = k |<e, z_c e'>|^2 for the eigenvectors e, e' of kappa_c (Lemma 4.6)
    t = k * (torch.einsum('bi,bij,bj->b', U[:, :, 0].conj(), zb[:, 1, :2, :2], U[:, :, 1]).abs() ** 2)
    for x in (0.5, 1.0, 2.0):
        assert (t < x).double().mean().item() == pytest.approx(1 - math.exp(-x), abs=0.01)
    assert t.mean().item() == pytest.approx(1.0, rel=0.03)


# ----------------------------------------------------------------------------
# Theorem 4.28 (idealised algorithm, dense)
# ----------------------------------------------------------------------------

def _polar(M):
    A, _, Bh = torch.linalg.svd(M)
    return A @ Bh


def test_theorem_4_28_idealised_algorithm():
    """Theorem 4.28(a)-(c): clusters, eps_1, zeta and the N^{1/2}(2 eps_1 + 3^h xi) bound for V."""
    g = gen(11)
    cases = [((1, 2), (3, 2)), ((2, 1), (2, 2)), ((1, 1, 2), (2, 1, 2))]
    checked_V = 0
    for trial in range(18):
        ks, ms = cases[trial % 3]
        H, _ = hidden(ks, ms, g)
        n, d, N = H.shape[-1], sum(m * m for m in ms), sum(ms)
        L = comm_matrix(H)
        gamma = svals_up(L)[d].item()
        Pi, F = null_projector(L, d)
        eta = gamma * 10 ** (-3 - 3 * torch.rand(1, generator=g).item())
        Pit, Ft = null_projector(comm_matrix(noisy(H, eta, g)), d)
        eps0 = opnorm(Pit - Pi)
        assert eps0 <= 4 * eta / gamma
        basis_t = Ft.reshape(d, n, n)                                      # orthonormal basis of the computed commutant
        Yp = torch.einsum('r,rij->ij', cgauss(g, d), basis_t)
        Kt = Yp + Yp.mH
        Yt = torch.einsum('r,rij->ij', cgauss(g, d), basis_t)
        K = (Pi @ Kt.reshape(-1)).reshape(n, n)
        Y = (Pi @ Yt.reshape(-1)).reshape(n, n)
        assert torch.allclose(K, K.mH, atol=1e-12)
        w, Q = torch.linalg.eigh(K)
        cuts = [0] + [i + 1 for i in range(n - 1) if w[i + 1] - w[i] > 1e-8] + [n]
        groups = list(zip(cuts[:-1], cuts[1:]))
        assert len(groups) == N
        delta = min((w[b0] - w[a1 - 1]).item() for (a0, a1), (b0, b1) in zip(groups[:-1], groups[1:]))
        rho = eps0 * fro(Kt)
        eps1 = 2 * rho / delta
        if not 4 * rho < delta:
            continue
        P = [Q[:, a:b] @ Q[:, a:b].mH for a, b in groups]
        wt, Qt = torch.linalg.eigh(Kt)
        cuts_t = [0] + [i + 1 for i in range(n - 1) if wt[i + 1] - wt[i] > 2 * rho] + [n]
        assert cuts_t == cuts                                              # (a): tau = 2 rho gives the N groups
        Qtl = [Qt[:, a:b] for a, b in groups]
        Pt = [q @ q.mH for q in Qtl]
        assert max(opnorm(Pt[l] - P[l]) for l in range(N)) <= eps1
        size = [b - a for a, b in groups]
        same = {(l, lp): opnorm(P[l] @ Y @ P[lp]) > 1e-10 for l in range(N) for lp in range(N) if l != lp}
        zeta = (eps0 + 2 * eps1) * fro(Yt)
        beta = min([fro(P[l] @ Y @ P[lp]) for (l, lp), s in same.items() if s] or [float('inf')])
        for (l, lp), s in same.items():
            if size[l] == size[lp]:
                assert abs(fro(Pt[l] @ Yt @ Pt[lp]) - fro(P[l] @ Y @ P[lp])) <= zeta
            if not s:
                assert fro(P[l] @ Y @ P[lp]) < 1e-10                       # Lemma 4.8(b)
        if not 2 * zeta < beta:
            continue
        # (b): links at theta' = zeta are exactly the same-sector pairs
        for (l, lp), s in same.items():
            if size[l] == size[lp]:
                assert (fro(Pt[l] @ Yt @ Pt[lp]) > zeta) == s
        # (c): star alignment of every sector to its first copy (a tree of depth 1), exact and computed
        kmax, h = max(ks), 1
        xi = kmax ** 0.5 * (eps0 + 4 * eps1) * fro(Yt) / beta
        if 3 ** h * xi > 1:
            continue
        Qex = []
        for l in range(N):                                                 # Lemma 4.24(c)
            S_ = Qtl[l].mH @ P[l] @ Qtl[l]
            ws, Us = torch.linalg.eigh(S_)
            Qex.append(P[l] @ Qtl[l] @ Us @ torch.diag(ws ** -0.5).to(C128) @ Us.mH)
        classes = {}
        for l in range(N):
            root = next((r for r in classes if same.get((r, l), False)), None)
            classes.setdefault(l if root is None else root, []).append(l)
        cols_t, cols_e = [], []
        for root, mem in classes.items():
            for l in mem:
                if l == root:
                    Wt = We = torch.eye(size[l], dtype=C128)
                else:
                    Wt = _polar(Qtl[root].mH @ Yt @ Qtl[l]).mH
                    We = _polar(Qex[root].mH @ Y @ Qex[l]).mH
                cols_t.append(Qtl[l] @ Wt)
                cols_e.append(Qex[l] @ We)
            # the exact aligned copies carry the same matrices: V^* H V is block diagonal with equal blocks
            ref = Qex[root]
            for l, col in zip(mem, cols_e[-len(mem):]):
                for Hj in H:
                    assert torch.allclose(col.mH @ Hj @ col, ref.mH @ Hj @ ref, atol=1e-9)
        Vt, Vh = torch.cat(cols_t, 1), torch.cat(cols_e, 1)
        assert torch.allclose(Vh.mH @ Vh, torch.eye(n, dtype=C128), atol=1e-10)
        assert opnorm(Vt - Vh) <= math.sqrt(N) * (2 * eps1 + 3 ** h * xi)
        checked_V += 1
    assert checked_V >= 6


def test_theorem_4_28_density_comparison_map():
    """Theorem 4.28(d): (1 - eps0^2)|x|^2 <= |Pi x|^2 <= |x|^2 on the computed commutant."""
    g = gen(12)
    H, _ = hidden((1, 2), (2, 2), g)
    n, d = H.shape[-1], 8
    L = comm_matrix(H)
    gamma = svals_up(L)[d].item()
    Pi, F = null_projector(L, d)
    _, Ft = null_projector(comm_matrix(noisy(H, 0.05 * gamma, g)), d)
    eps0 = opnorm(Ft.T @ Ft.conj() - Pi)
    B = F.conj() @ Ft.T                                                     # matrix of Pi: ran Pi~ -> ran Pi
    sv = torch.linalg.svdvals(B)
    assert sv.max() <= 1 + 1e-12 and sv.min() ** 2 >= 1 - eps0 ** 2 - 1e-12
    assert torch.prod(sv).item() ** 2 >= (1 - eps0 ** 2) ** d - 1e-12       # real determinant = |complex det|^2
