"""
Renyi divergences, Petz recovery maps and L^p norms (torch_vn_algebra.divergences), checked against independent
dense references: block-diagonal NumPy matrices with scipy.linalg.fractional_matrix_power / numpy.linalg.eigh, the
classical formulas for commuting states, generalised eigenvalues for D_max, and known theorems (Araki-Lieb-Thirring,
data processing, Petz's equality condition, the recoverability bound of Junge-Renner-Sutter-Wilde-Winter).
"""
import math
import warnings

import numpy as np
import pytest
import scipy.linalg as sl
import torch

from torch_vn_algebra import TypeIAlgebra, DensityMatrix, states as st, channels as ch, divergences as dv
from torch_vn_algebra.composite import tensor_product, kron, partial_trace_channel
from torch_vn_algebra.states import _trace_weights

KINDS = ('petz', 'sandwiched')
DIV = {'petz': dv.petz_renyi_divergence, 'sandwiched': dv.sandwiched_renyi_divergence}


# ----------------------------------------------------------------------------------------------------------------
# dense references
# ----------------------------------------------------------------------------------------------------------------
def _alg(ks=(3, 0, 2, 4), complex_=True, precision='double'):
    return TypeIAlgebra([max(k, 1) for k in ks], list(ks), complex_valued=complex_, precision=precision)


def _dense(alg, m):
    """Block-diagonal numpy matrix of one sample (empty sectors skipped)."""
    return sl.block_diag(*[m[c, :k, :k].detach().numpy() for c, k in enumerate(alg.k_factors) if k])


def _blocks_np(alg, m):
    return [m[c, :k, :k].detach().numpy() for c, k in enumerate(alg.k_factors) if k]


def _pow(A, p, rtol=None):
    """A^p on the support (eigenvalues above rtol * lambda_max, rtol = n eps by default), 0 elsewhere."""
    A = 0.5 * (A + A.conj().T)
    w, V = np.linalg.eigh(A)
    thr = (A.shape[0] * np.finfo(float).eps if rtol is None else rtol) * max(w.max(), 0)
    big = w > thr
    f = np.zeros_like(w)
    f[big] = w[big] ** p
    return (V * f) @ V.conj().T


def _fpow(A, p):
    """A^p of a positive definite A by scipy (Schur-Pade; independent of eigh)."""
    return np.asarray(sl.fractional_matrix_power(A, p))


def _petz_ref(r, s, a, power=_pow):
    return np.log(np.trace(power(r, a) @ power(s, 1 - a)).real / np.trace(r).real) / (a - 1)


def _sand_ref(r, s, a, power=_pow):
    """~D_alpha with the nonzero spectrum of S r S (S = s^g) from the r x r matrix G^* S^2 G, r = G G^* with
    G = V_r diag(w_r)^{1/2} of rank r (rank-revealing: no rounding-level eigenvalues enter)."""
    S = power(s, (1 - a) / (2 * a))
    w, V = np.linalg.eigh(0.5 * (r + r.conj().T))
    keep = w > r.shape[0] * np.finfo(float).eps * w.max()
    G = V[:, keep] * np.sqrt(w[keep])
    lam = np.linalg.eigvalsh(G.conj().T @ S @ S @ G)
    return np.log(np.sum(np.clip(lam, 0, None) ** a) / np.trace(r).real) / (a - 1)


REF = {'petz': _petz_ref, 'sandwiched': _sand_ref}


def _umegaki_ref(r, s):
    wr = np.linalg.eigvalsh(r)
    wr = wr[wr > 1e-14]
    L = sl.logm(s)
    return float(np.sum(wr * np.log(wr)) - np.trace(r @ L).real)


def _dmax_ref(r, s):
    """log of the largest generalised eigenvalue r x = lambda s x (s positive definite)."""
    return float(np.log(sl.eigh(r, s, eigvals_only=True).max()))


def _random_states(alg, B, seed, rank=None):
    torch.manual_seed(seed)
    return st.random_density_matrix(alg, B, rank=rank), st.random_density_matrix(alg, B)


def _positive(alg, m):
    return alg.operator(m, is_self_adjoint=True, is_positive=True)


# ----------------------------------------------------------------------------------------------------------------
# definitions
# ----------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('complex_', [True, False])
@pytest.mark.parametrize('kind', KINDS)
@pytest.mark.parametrize('alpha', [0.3, 0.7, 1.5, 2.5])
def test_dense_reference_full_rank(complex_, kind, alpha):
    """Unequal sectors (with an empty one), batch 3, against scipy's fractional matrix powers."""
    alg = _alg(complex_=complex_)
    rho, sigma = _random_states(alg, 3, 0)
    D = DIV[kind](rho, sigma, alpha)
    for b in range(3):
        r, s = _dense(alg, rho.matrix[b]), _dense(alg, sigma.matrix[b])
        assert abs(D[b].item() - REF[kind](r, s, alpha, _fpow)) < 1e-10
    assert D.dtype == torch.float64 and D.shape == (3,)


@pytest.mark.parametrize('complex_', [True, False])
@pytest.mark.parametrize('kind', KINDS)
def test_rank_deficient_states(complex_, kind):
    alg = _alg(complex_=complex_)
    torch.manual_seed(1)
    rho_low = st.random_density_matrix(alg, 3, rank=1)            # rank 1 in every block
    sigma = st.random_density_matrix(alg, 3)
    tau = st.random_density_matrix(alg, 3, rank=1)
    sigma_low = rho_low.mix(tau, 0.5)                              # supp rho_low < supp sigma_low (singular)
    for alpha in (0.3, 0.7, 1.5, 2.5):
        # rank-deficient rho, faithful sigma: finite for every alpha
        D = DIV[kind](rho_low, sigma, alpha)
        # both singular with supp rho in supp sigma: finite (inverse powers on the support)
        D2 = DIV[kind](rho_low, sigma_low, alpha)
        # faithful rho, singular sigma: +inf for alpha > 1, finite for alpha < 1
        D3 = DIV[kind](sigma, sigma_low, alpha)
        for b in range(3):
            r, s, s2 = (_dense(alg, x.matrix[b]) for x in (rho_low, sigma, sigma_low))
            assert abs(D[b].item() - REF[kind](r, s, alpha)) < 1e-10
            assert abs(D2[b].item() - REF[kind](r, s2, alpha)) < 1e-9
            if alpha > 1:
                assert D3[b].item() == math.inf
            else:
                assert abs(D3[b].item() - REF[kind](s, s2, alpha)) < 1e-9


def test_support_conventions():
    """Orthogonal supports give +inf for every alpha; supp rho not in supp sigma gives +inf for alpha > 1, for the
    Umegaki relative entropy and for D_max; D_0 = -log Tr(Pi_rho sigma)."""
    alg = _alg((2, 3))
    z = torch.zeros(1, 2, 3, 3, dtype=alg.dtype)
    a, b = z.clone(), z.clone()
    a[0, 0, :2, :2] = torch.eye(2) / 2                            # rho in sector 0
    b[0, 1, :3, :3] = torch.eye(3) / 3                            # sigma in sector 1
    rho, sigma = DensityMatrix(alg, matrix=a), DensityMatrix(alg, matrix=b)
    for kind in KINDS:
        for alpha in (0.3, 0.8, 1.5):
            assert DIV[kind](rho, sigma, alpha).item() == math.inf
    assert dv.umegaki_relative_entropy(rho, sigma).item() == math.inf
    assert dv.max_relative_entropy(rho, sigma).item() == math.inf
    assert dv.min_relative_entropy(rho, sigma).item() == math.inf
    mix = rho.mix(sigma, 0.5)
    for kind in KINDS:
        assert DIV[kind](mix, sigma, 1.5).item() == math.inf
        assert abs(DIV[kind](mix, sigma, 0.5).item() - (-2 * math.log(0.5 ** 0.5) if kind == 'petz'
                                                         else -math.log(0.5))) < 1e-12
    assert abs(dv.min_relative_entropy(sigma, mix).item() - math.log(2)) < 1e-12        # -log Tr(Pi_sigma mix)
    assert abs(dv.umegaki_relative_entropy(sigma, mix).item() - math.log(2)) < 1e-12
    assert abs(dv.max_relative_entropy(sigma, mix).item() - math.log(2)) < 1e-12


def test_commuting_states_are_classical():
    """For commuting rho, sigma both families are the classical Renyi divergence of the eigenvalue distributions
    (with zero probabilities in p); also Umegaki = KL, D_max = log max p/q, D_0 = -log q(supp p)."""
    torch.manual_seed(2)
    alg = _alg((3, 2, 4))
    B = 4
    N = sum(alg.k_factors)
    p = torch.rand(B, N, dtype=torch.float64)
    p[:, ::3] = 0.0                                                 # some zero probabilities
    q = torch.rand(B, N, dtype=torch.float64) + 0.05
    p, q = p / p.sum(-1, keepdim=True), q / q.sum(-1, keepdim=True)

    def diag(v):
        m = torch.zeros(B, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype)
        o = 0
        for c, k in enumerate(alg.k_factors):
            m[:, c, :k, :k] = torch.diag_embed(v[:, o:o + k]).to(alg.dtype)
            o += k
        return m
    U = alg.random_unitary_operator(B).matrix
    rot = lambda m: U @ m @ U.conj().transpose(-2, -1)
    rho, sigma = DensityMatrix(alg, matrix=rot(diag(p))), DensityMatrix(alg, matrix=rot(diag(q)))
    pn, qn = p.numpy(), q.numpy()
    m = pn > 0
    for alpha in (0.2, 0.5, 0.9, 1.3, 2.0, 4.0):
        ref = np.array([np.log(np.sum(pn[b, m[b]] ** alpha * qn[b, m[b]] ** (1 - alpha))) / (alpha - 1)
                        for b in range(B)])
        for kind in KINDS:
            assert np.allclose(DIV[kind](rho, sigma, alpha).numpy(), ref, atol=1e-10)
    kl = np.array([np.sum(pn[b, m[b]] * np.log(pn[b, m[b]] / qn[b, m[b]])) for b in range(B)])
    assert np.allclose(dv.umegaki_relative_entropy(rho, sigma).numpy(), kl, atol=1e-10)
    assert np.allclose(dv.max_relative_entropy(rho, sigma).numpy(), np.log((pn / qn).max(-1)), atol=1e-10)
    assert np.allclose(dv.min_relative_entropy(rho, sigma).numpy(),
                       -np.log((qn * m).sum(-1)), atol=1e-10)


@pytest.mark.parametrize('complex_', [True, False])
def test_limits(complex_):
    """alpha -> 1: Umegaki (error O(|alpha - 1|), from both sides); alpha -> inf (sandwiched): D_max from the
    generalised eigenvalues; alpha = 1/2 (sandwiched): -log F; alpha -> 0 (Petz): D_0."""
    alg = _alg(complex_=complex_)
    rho, sigma = _random_states(alg, 3, 3)
    D1 = dv.umegaki_relative_entropy(rho, sigma)
    assert torch.allclose(D1, st.relative_entropy(rho, sigma), atol=1e-12)
    for b in range(3):
        assert abs(D1[b].item() - _umegaki_ref(_dense(alg, rho.matrix[b]), _dense(alg, sigma.matrix[b]))) < 1e-10
    d = 1e-4
    for kind in KINDS:
        lo, hi = DIV[kind](rho, sigma, 1 - d), DIV[kind](rho, sigma, 1 + d)
        assert ((D1 - lo).abs() < 10 * d).all() and ((hi - D1).abs() < 10 * d).all()
        assert (lo <= D1 + 1e-12).all() and (D1 <= hi + 1e-12).all()
        assert torch.equal(dv.renyi_divergence(rho, sigma, 1.0, kind), D1)
    dmax = dv.max_relative_entropy(rho, sigma)
    for b in range(3):
        assert abs(dmax[b].item() - _dmax_ref(_dense(alg, rho.matrix[b]), _dense(alg, sigma.matrix[b]))) < 1e-10
    big = dv.sandwiched_renyi_divergence(rho, sigma, 1e6)
    assert ((dmax - big).abs() < 1e-4).all() and (big <= dmax + 1e-12).all()
    assert torch.equal(dv.renyi_divergence(rho, sigma, math.inf, 'sandwiched'), dmax)
    assert torch.allclose(dv.sandwiched_renyi_divergence(rho, sigma, 0.5), -torch.log(st.fidelity(rho, sigma)),
                          atol=1e-12)
    low = st.random_density_matrix(alg, 3, rank=1)
    d0 = dv.min_relative_entropy(low, sigma)
    assert torch.equal(dv.renyi_divergence(low, sigma, 0.0, 'petz'), d0)
    assert ((dv.petz_renyi_divergence(low, sigma, 1e-7) - d0).abs() < 1e-5).all()
    for b in range(3):
        r, s = _dense(alg, low.matrix[b]), _dense(alg, sigma.matrix[b])
        assert abs(d0[b].item() + np.log(np.trace(_pow(r, 0) @ s).real)) < 1e-10


def test_ordering_monotonicity_positivity():
    """~D_alpha <= D_alpha (Araki-Lieb-Thirring); both nondecreasing in alpha (Petz on (0, inf), sandwiched on
    [1/2, inf)); D >= 0 with equality iff rho = sigma."""
    alg = _alg((3, 2, 4))
    rho, sigma = _random_states(alg, 16, 4)
    rank1 = st.random_density_matrix(alg, 16, rank=1)
    alphas = [0.1, 0.25, 0.5, 0.75, 0.95, 1.05, 1.5, 2.0, 3.0, 6.0]
    for r in (rho, rank1):
        P = torch.stack([dv.petz_renyi_divergence(r, sigma, a) for a in alphas])
        S = torch.stack([dv.sandwiched_renyi_divergence(r, sigma, a) for a in alphas])
        assert (S <= P + 1e-12).all()
        assert (P[1:] >= P[:-1] - 1e-12).all()
        assert (S[3:] >= S[2:-1] - 1e-12).all()                         # alpha >= 1/2
        assert (P > 1e-6).all() and (S[2:] > 1e-6).all()
    for kind in KINDS:
        for a in (0.3, 0.5, 2.0):
            assert DIV[kind](rho, rho, a).abs().max() < 1e-12
            assert DIV[kind](rank1, rank1, a).abs().max() < 1e-12
    assert dv.umegaki_relative_entropy(rho, rho).abs().max() < 1e-12
    assert dv.max_relative_entropy(rho, rho).abs().max() < 1e-12


def test_per_sample_alpha_and_broadcasting():
    alg = _alg()
    rho, sigma = _random_states(alg, 5, 5)
    alphas = torch.tensor([0.3, 0.6, 1.4, 2.0, 7.0], dtype=torch.float64)
    for kind in KINDS:
        D = DIV[kind](rho, sigma, alphas)
        ref = torch.stack([DIV[kind](rho, sigma, float(a))[i] for i, a in enumerate(alphas)])
        assert torch.allclose(D, ref, atol=1e-13, rtol=0)
        # one state against a batch, one alpha per sample
        one = DensityMatrix(alg, matrix=rho.matrix[:1], validate=False)
        D1 = DIV[kind](one, sigma, alphas)
        ref1 = torch.stack([DIV[kind](one, sigma, float(a))[i] for i, a in enumerate(alphas)])
        assert torch.allclose(D1, ref1, atol=1e-13, rtol=0)
        with pytest.raises(ValueError):
            DIV[kind](rho, sigma, torch.ones(3, dtype=torch.float64))   # batch 3 against 5


def test_unnormalised_inputs_scale():
    """For positive rho, sigma: D(a rho||b sigma) = D(rho||sigma) + log(a/b) (normalisation by Tr rho)."""
    alg = _alg()
    rho, sigma = _random_states(alg, 2, 6)
    a, b = 3.0, 0.25
    R, S = _positive(alg, rho.matrix * a), _positive(alg, sigma.matrix * b)
    for kind in KINDS:
        for alpha in (0.4, 2.0):
            assert torch.allclose(DIV[kind](R, S, alpha), DIV[kind](rho, sigma, alpha) + math.log(a / b), atol=1e-12)
    assert torch.allclose(dv.umegaki_relative_entropy(R, S), dv.umegaki_relative_entropy(rho, sigma)
                          + math.log(a / b), atol=1e-12)
    assert torch.allclose(dv.max_relative_entropy(R, S), dv.max_relative_entropy(rho, sigma) + math.log(a / b),
                          atol=1e-12)


def test_single_precision():
    alg = _alg(precision='single')
    torch.manual_seed(7)
    rho, sigma = st.random_density_matrix(alg, 3), st.random_density_matrix(alg, 3)
    alg64 = alg.like(precision='double')
    r64, s64 = rho.cast(alg64), sigma.cast(alg64)
    for kind in KINDS:
        for alpha in (0.4, 1.7, 50.0):
            D, D64 = DIV[kind](rho, sigma, alpha), DIV[kind](r64, s64, alpha)
            assert D.dtype == torch.float32 and torch.isfinite(D).all()
            assert torch.allclose(D.double(), D64, rtol=2e-4, atol=2e-4)
    assert torch.allclose(dv.sandwiched_renyi_divergence(rho, sigma, 1e5).double(),
                          dv.max_relative_entropy(r64, s64), atol=1e-3)


# ----------------------------------------------------------------------------------------------------------------
# data processing
# ----------------------------------------------------------------------------------------------------------------
def _channels(alg, B):
    torch.manual_seed(8)
    out_alg = TypeIAlgebra([2, 5], [2, 5], precision='double')
    return {'random Stinespring': ch.random_channel(alg, 2, B),
            'random inter-sector': ch.random_inter_sector_channel(alg, out_alg, 2, B),
            'centre expectation': ch.center_expectation(alg),
            'dephasing': ch.dephasing_channel(alg, 0.7),
            'amplitude damping': ch.amplitude_damping_channel(alg, 0.4)}


@pytest.mark.parametrize('name', ['random Stinespring', 'random inter-sector', 'centre expectation', 'dephasing',
                                  'amplitude damping'])
def test_data_processing(name):
    """Sandwiched DPI for alpha >= 1/2 (Frank-Lieb, Beigi), Petz DPI for alpha in (0, 2] (Petz), Umegaki, D_max,
    on 64 random pairs of states (rank-deficient and faithful)."""
    alg = _alg((3, 2, 4))
    B = 64
    torch.manual_seed(9)
    rho = st.random_density_matrix(alg, B)
    rho = DensityMatrix(alg, matrix=torch.cat([rho.matrix[:32], st.random_density_matrix(alg, 32, rank=1).matrix]),
                        validate=False)
    sigma = st.random_density_matrix(alg, B)
    Phi = _channels(alg, B)[name]
    pr, ps = Phi(rho), Phi(sigma)
    tol = 1e-10
    for a in (0.5, 0.75, 1.5, 3.0, 10.0):
        assert (dv.sandwiched_renyi_divergence(pr, ps, a) <= dv.sandwiched_renyi_divergence(rho, sigma, a) + tol).all()
    for a in (0.1, 0.5, 0.9, 1.5, 2.0):
        assert (dv.petz_renyi_divergence(pr, ps, a) <= dv.petz_renyi_divergence(rho, sigma, a) + tol).all()
    assert (dv.umegaki_relative_entropy(pr, ps) <= dv.umegaki_relative_entropy(rho, sigma) + tol).all()
    assert (dv.max_relative_entropy(pr, ps) <= dv.max_relative_entropy(rho, sigma) + tol).all()


# ----------------------------------------------------------------------------------------------------------------
# trace conventions and sector decomposition
# ----------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('trace', ['norm', 'tau_vN'])
def test_trace_convention_invariance(trace):
    """The divergences computed with tau_w from the densities d = rho_c / w_c equal those of the states: checked
    against the weighted formulas evaluated independently in numpy, and through the trace argument."""
    alg = _alg((3, 2, 4))
    rho, sigma = _random_states(alg, 3, 10)
    w = _trace_weights(alg, trace).numpy()
    dr, ds = rho.density(trace), sigma.density(trace)
    Dr, Ds = _positive(alg, dr), _positive(alg, ds)
    for a in (0.4, 0.5, 1.6, 3.0):
        for b in range(3):
            br, bs = _blocks_np(alg, dr[b]), _blocks_np(alg, ds[b])
            qp = sum(wc * np.trace(_fpow(x, a) @ _fpow(y, 1 - a)).real for wc, x, y in zip(w, br, bs))
            g = (1 - a) / (2 * a)
            qs = sum(wc * np.sum(np.linalg.eigvalsh(_fpow(y, g) @ x @ _fpow(y, g)) ** a)
                     for wc, x, y in zip(w, br, bs))
            assert abs(np.log(qp) / (a - 1) - dv.petz_renyi_divergence(rho, sigma, a)[b].item()) < 1e-10
            assert abs(np.log(qs) / (a - 1) - dv.sandwiched_renyi_divergence(rho, sigma, a)[b].item()) < 1e-10
        for kind in KINDS:
            assert torch.allclose(DIV[kind](Dr, Ds, a, trace=trace), DIV[kind](rho, sigma, a), atol=1e-12)
    for f in (dv.umegaki_relative_entropy, dv.max_relative_entropy, dv.min_relative_entropy):
        assert torch.allclose(f(Dr, Ds, trace=trace), f(rho, sigma), atol=1e-12)
    with pytest.raises(ValueError):
        dv.petz_renyi_divergence(rho, sigma, 0.5, trace=trace)        # a DensityMatrix is a Tr_blunt density


@pytest.mark.parametrize('kind', ['petz', 'sandwiched', 'umegaki'])
def test_sector_contributions(kind):
    """Q = sum_c Q_c and Q_c = p_c^alpha q_c^{1-alpha} exp((alpha - 1) D(rho_c/p_c||sigma_c/q_c)) (Renyi),
    D = D(p||q) + sum_c p_c D(rho_c/p_c||sigma_c/q_c) (Umegaki); conditional divergences against numpy; NaN for
    sectors without weight (an empty sector and a sector where rho vanishes)."""
    alg = _alg((3, 0, 2, 4))
    rho, sigma = _random_states(alg, 3, 11)
    m = rho.matrix.clone()
    m[:, 2] = 0                                                      # rho has no weight in sector 2
    rho = DensityMatrix(alg, matrix=m / torch.diagonal(m, dim1=-2, dim2=-1).sum((-2, -1))[:, None, None, None])
    for a in ((0.4, 1.8) if kind != 'umegaki' else (None,)):
        sc = dv.sector_contributions(rho, sigma, a, kind)
        assert torch.isnan(sc.conditional[:, 1]).all() and torch.isnan(sc.conditional[:, 2]).all()
        ok = [0, 3]
        p, q, cond = sc.p[:, ok], sc.q[:, ok], sc.conditional[:, ok]
        if kind == 'umegaki':
            D = dv.umegaki_relative_entropy(rho, sigma)
            assert torch.allclose(sc.terms.sum(-1), D, atol=1e-12)
            kl = (p * torch.log(p / sc.q[:, ok])).sum(-1)
            assert torch.allclose(kl + (p * cond).sum(-1), D, atol=1e-12)
        else:
            D = DIV[kind](rho, sigma, a)
            assert torch.allclose(torch.log(sc.terms.sum(-1)) / (a - 1), D, atol=1e-12)
            assert torch.allclose(sc.terms[:, ok], p ** a * q ** (1 - a) * torch.exp((a - 1) * cond), atol=1e-12)
            assert torch.allclose(torch.log((p ** a * q ** (1 - a) * torch.exp((a - 1) * cond)).sum(-1)) / (a - 1), D,
                                  atol=1e-12)
        for b in range(3):
            for j, c in enumerate(ok):
                k = alg.k_factors[c]
                r = rho.matrix[b, c, :k, :k].numpy() / p[b, j].item()
                s = sigma.matrix[b, c, :k, :k].numpy() / q[b, j].item()
                ref = _umegaki_ref(r, s) if kind == 'umegaki' else REF[kind](r, s, a)
                assert abs(cond[b, j].item() - ref) < 1e-10


# ----------------------------------------------------------------------------------------------------------------
# Petz recovery map
# ----------------------------------------------------------------------------------------------------------------
def _choi_of_map(R, alg_out, b=0):
    """Choi matrices sum_ij |i><j| (x) R(E_ij) of a map M_out -> M_in, for every pair of sectors (d, c), sample b."""
    out = []
    for d, m in enumerate(alg_out.k_factors):
        cols = []
        for i in range(m):
            for j in range(m):
                E = torch.zeros(1, alg_out.C, alg_out.k_max, alg_out.k_max, dtype=alg_out.dtype)
                E[0, d, i, j] = 1
                cols.append(R(alg_out.operator(E)).matrix[b])
        for c, k in enumerate(R.alg_in.k_factors):
            J = np.zeros((m * k, m * k), dtype=complex)
            for i in range(m):
                for j in range(m):
                    J[i * k:(i + 1) * k, j * k:(j + 1) * k] = cols[i * m + j][c, :k, :k].numpy()
            out.append(J)
    return out


def _dense_petz(alg_in, alg_out, sigma, kraus, X, b=0):
    """R(X) for a sector-preserving channel from dense block-diagonal matrices."""
    s = _dense(alg_in, sigma[b])
    Ks = [_dense(alg_in, K) for K in kraus[min(b, kraus.shape[0] - 1)]]
    phis = sum(K @ s @ K.conj().T for K in Ks)
    Pm = _fpow(phis, -0.5)
    Y = Pm @ _dense(alg_out, X[b]) @ Pm
    return _fpow(s, 0.5) @ sum(K.conj().T @ Y @ K for K in Ks) @ _fpow(s, 0.5)


@pytest.mark.parametrize('complex_', [True, False])
def test_petz_map_dense_and_identities(complex_):
    alg = _alg((3, 2, 4), complex_=complex_)
    torch.manual_seed(12)
    B = 3
    sigma, rho = st.random_density_matrix(alg, B), st.random_density_matrix(alg, B)
    Phi = ch.random_channel(alg, 2, B)
    R = dv.PetzRecovery(sigma, Phi)
    out = R(Phi(rho))
    for b in range(B):
        assert np.abs(_dense(alg, out.matrix[b]) - _dense_petz(alg, alg, sigma.matrix, Phi.kraus, Phi(rho).matrix, b)
                      ).max() < 1e-10
    assert (R(Phi(sigma)).matrix - sigma.matrix).abs().max() < 1e-10            # R(Phi(sigma)) = sigma
    # trace preserving (Phi(sigma) faithful) and completely positive
    X = alg.operator(torch.randn(B, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype))
    assert torch.allclose(R(X).trace, X.trace, atol=1e-10)
    for J in _choi_of_map(R, alg, b=1):
        assert np.linalg.eigvalsh(0.5 * (J + J.conj().T)).min() > -1e-10
    assert out._is_positive and out._is_self_adjoint


def test_petz_map_inter_sector_channel():
    alg_in = _alg((3, 2), complex_=True)
    alg_out = TypeIAlgebra([2, 3, 1], [2, 3, 1], precision='double')
    torch.manual_seed(13)
    B = 2
    sigma = st.random_density_matrix(alg_in, B)
    Phi = ch.random_inter_sector_channel(alg_in, alg_out, 2, B)
    R = dv.PetzRecovery(sigma, Phi)
    assert (R(Phi(sigma)).matrix - sigma.matrix).abs().max() < 1e-10
    X = st.random_density_matrix(alg_out, B)
    assert torch.allclose(R(X).trace.real, torch.ones(B, dtype=torch.float64), atol=1e-10)
    # dense reference with the inter-sector Kraus blocks
    for b in range(B):
        s = _dense(alg_in, sigma.matrix[b])
        K = Phi.kraus[b]                                                  # (r, D, C, m, k)
        offs_in = np.cumsum([0] + alg_in.k_factors)
        offs_out = np.cumsum([0] + alg_out.k_factors)
        Ks = []
        for i in range(K.shape[0]):
            Kd = np.zeros((offs_out[-1], offs_in[-1]), dtype=complex)
            for d, m in enumerate(alg_out.k_factors):
                for c, k in enumerate(alg_in.k_factors):
                    Kd[offs_out[d]:offs_out[d] + m, offs_in[c]:offs_in[c] + k] = K[i, d, c, :m, :k].numpy()
            Ks.append(Kd)
        phis = sum(Kd @ s @ Kd.conj().T for Kd in Ks)
        assert np.allclose(_dense(alg_out, Phi(sigma).matrix[b]), phis, atol=1e-12)
        Pm = _fpow(phis, -0.5)
        Y = Pm @ _dense(alg_out, X.matrix[b]) @ Pm
        ref = _fpow(s, 0.5) @ sum(Kd.conj().T @ Y @ Kd for Kd in Ks) @ _fpow(s, 0.5)
        assert np.abs(_dense(alg_in, R(X).matrix[b]) - ref).max() < 1e-10
    for J in _choi_of_map(R, alg_out):
        assert np.linalg.eigvalsh(0.5 * (J + J.conj().T)).min() > -1e-10


def test_petz_map_of_centre_expectation():
    """For sigma in the centre and Phi = E_Z (self-adjoint for Tr_blunt), R_{sigma,E} = E on M."""
    alg = _alg((3, 2, 4))
    torch.manual_seed(14)
    E = ch.center_expectation(alg)
    sigma = E(st.random_density_matrix(alg, 2))
    R = dv.PetzRecovery(sigma, E)
    X = alg.operator(torch.randn(2, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype))
    assert (R(X).matrix - E(X).matrix).abs().max() < 1e-12


def test_petz_equality_condition():
    """Petz's theorem: D(rho||sigma) = D(Phi rho||Phi sigma) iff R_{sigma,Phi}(Phi(rho)) = rho. Sufficient instances:
    (i) rho = rho_A (x) tau, sigma = sigma_A (x) tau under the partial trace; (ii) equal conditional states in every
    sector under the conditional expectation onto the centre. Generic instances: strict inequality and imperfect
    recovery."""
    torch.manual_seed(15)
    a1, a2 = _alg((2, 1)), _alg((3,))
    a12 = tensor_product(a1, a2)
    B = 3
    tau = st.random_density_matrix(a2, B)
    rA, sA = st.random_density_matrix(a1, B), st.random_density_matrix(a1, B)
    rho, sigma = kron(rA, tau, a12), kron(sA, tau, a12)
    Tr2 = partial_trace_channel(a12, keep=1)
    R = dv.PetzRecovery(sigma, Tr2)
    for f in (dv.umegaki_relative_entropy, lambda r, s: dv.petz_renyi_divergence(r, s, 0.6),
              lambda r, s: dv.sandwiched_renyi_divergence(r, s, 1.7)):
        assert torch.allclose(f(Tr2(rho), Tr2(sigma)), f(rho, sigma), atol=1e-12)
    assert (R(Tr2(rho)).matrix - rho.matrix).abs().max() < 1e-12
    # generic: strict
    rho_g = st.random_density_matrix(a12, B)
    gap = dv.umegaki_relative_entropy(rho_g, sigma) - dv.umegaki_relative_entropy(Tr2(rho_g), Tr2(sigma))
    assert (gap > 1e-3).all()
    assert ((R(Tr2(rho_g)).matrix - rho_g.matrix).abs().amax(dim=(1, 2, 3)) > 1e-3).all()
    # (ii) centre: rho_c = p_c omega_c, sigma_c = q_c omega_c
    alg = _alg((3, 2, 4))
    omega = st.random_density_matrix(alg, B)
    pw = torch.rand(B, alg.C, dtype=torch.float64)
    qw = torch.rand(B, alg.C, dtype=torch.float64)
    cond = omega.matrix / torch.diagonal(omega.matrix, dim1=-2, dim2=-1).sum(-1)[:, :, None, None]
    rho = DensityMatrix(alg, matrix=cond * (pw / pw.sum(-1, keepdim=True))[:, :, None, None])
    sigma = DensityMatrix(alg, matrix=cond * (qw / qw.sum(-1, keepdim=True))[:, :, None, None])
    E = ch.center_expectation(alg)
    assert torch.allclose(dv.umegaki_relative_entropy(E(rho), E(sigma)), dv.umegaki_relative_entropy(rho, sigma),
                          atol=1e-12)
    assert (dv.PetzRecovery(sigma, E)(E(rho)).matrix - rho.matrix).abs().max() < 1e-12


def test_rotated_and_universal_recovery():
    """R^t(Phi(sigma)) = sigma for every t; the quadrature of the universal map is converged; the recoverability
    bound of Junge, Renner, Sutter, Wilde, Winter (2018, Thm. 2.1),
    D(rho||sigma) - D(Phi rho||Phi sigma) >= -int beta_0(t) log F(rho, R^{t/2}(Phi rho)) dt
                                         >= -log F(rho, R_univ(Phi rho)),
    on random instances (F = ||sqrt(rho) sqrt(omega)||_1^2)."""
    alg = _alg((3, 2))
    torch.manual_seed(16)
    B = 6
    rho, sigma = st.random_density_matrix(alg, B), st.random_density_matrix(alg, B)
    Phi = ch.random_channel(alg, 2, B)
    R = dv.PetzRecovery(sigma, Phi)
    for t in (0.0, 0.7, -2.3):
        assert (R.rotated(Phi(sigma), t).matrix - sigma.matrix).abs().max() < 1e-10
    out = R.rotated(Phi(rho), 0.9)
    assert torch.allclose(out.trace.real, torch.ones(B, dtype=torch.float64), atol=1e-10)
    U = dv.universal_recovery(sigma, Phi, Phi(rho))
    U2 = dv.universal_recovery(sigma, Phi, Phi(rho), step=0.0625, cutoff=14.0)
    assert (U.matrix - U2.matrix).abs().max() < 1e-10
    assert (dv.universal_recovery(sigma, Phi, Phi(sigma)).matrix - sigma.matrix).abs().max() < 1e-10
    gap = dv.umegaki_relative_entropy(rho, sigma) - dv.umegaki_relative_entropy(Phi(rho), Phi(sigma))
    bound_univ = -torch.log(st.fidelity(rho, U))
    ts, ws = dv.universal_recovery_nodes()
    bound_int = -sum(w * torch.log(st.fidelity(rho, R.rotated(Phi(rho), t / 2))) for t, w in zip(ts, ws))
    assert (gap >= bound_int - 1e-10).all() and (bound_int >= bound_univ - 1e-10).all()
    assert (bound_univ > 0).all()                                     # the bound is not trivial here
    with pytest.raises(ValueError):
        dv.PetzRecovery(st.random_density_matrix(_alg((2,), complex_=False)),
                        ch.random_channel(_alg((2,), complex_=False), 2), t=0.5)


# ----------------------------------------------------------------------------------------------------------------
# L^p norms
# ----------------------------------------------------------------------------------------------------------------
def _lp_ref(alg, m, p, w):
    s = [np.linalg.svd(x, compute_uv=False) for x in _blocks_np(alg, m)]
    if p == math.inf:
        return max(x.max() for x, wc in zip(s, w) if wc > 0)
    return sum(wc * np.sum(x ** p) for x, wc in zip(s, w)) ** (1 / p)


@pytest.mark.parametrize('trace', ['blunt', 'norm', 'tau_vN'])
@pytest.mark.parametrize('complex_', [True, False])
def test_lp_norms(trace, complex_):
    alg = _alg((3, 0, 2, 4), complex_=complex_)
    torch.manual_seed(17)
    B = 3
    A = alg.operator(torch.randn(B, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype))
    H = alg.operator(A.matrix + A.matrix.conj().transpose(-2, -1), is_self_adjoint=True)
    w = _trace_weights(alg, trace).numpy()
    act = [wc for wc, k in zip(w, alg.k_factors) if k]
    for X in (A, H):
        for p in (1.0, 1.5, 2.0, 3.0, 7.0, math.inf):
            n = dv.lp_norm(X, p, trace)
            for b in range(B):
                assert abs(n[b].item() - _lp_ref(alg, X.matrix[b], p, act)) < 1e-10 * max(1, n[b].item())
        assert torch.allclose(dv.lp_norm(X, math.inf, trace), X.operator_norm(), atol=1e-12)
        fro = (torch.tensor(w) * X.matrix.abs().pow(2).sum((-2, -1))).sum(-1).sqrt()
        assert torch.allclose(dv.lp_norm(X, 2, trace), fro, atol=1e-12)
    if trace == 'blunt':
        assert torch.allclose(dv.lp_norm(A, 1), A.trace_norm(), atol=1e-12)
        assert torch.allclose(dv.lp_norm(H, 1), H.trace_norm(), atol=1e-12)


def test_lp_holder_duality_commutative_and_weights():
    alg = _alg((3, 2, 4))
    torch.manual_seed(18)
    B = 8
    A = alg.operator(torch.randn(B, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype))
    Bop = alg.operator(torch.randn(B, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype))
    for trace in ('blunt', 'norm', 'tau_vN'):
        w = _trace_weights(alg, trace)
        pair = lambda X, Y: (w * torch.einsum('bcij,bcji->bc', X.matrix, Y.matrix)).sum(-1)
        for p in (1.0, 1.25, 2.0, 3.0, math.inf):
            q = math.inf if p == 1 else 1.0 if p == math.inf else p / (p - 1)
            assert (pair(A, Bop).abs() <= dv.lp_norm(A, p, trace) * dv.lp_norm(Bop, q, trace) + 1e-12).all()
            if p == math.inf:
                continue
            # the dual element B* = |A|^{p-1} W^* / ||A||_p^{p-1} attains the supremum (A = W |A|)
            n = dv.lp_norm(A, p, trace)
            Bs = torch.zeros_like(A.matrix)
            for c, k in enumerate(alg.k_factors):
                U, s, Vh = torch.linalg.svd(A.matrix[:, c, :k, :k])
                Bs[:, c, :k, :k] = (Vh.mH * (s ** (p - 1)).to(Vh.dtype).unsqueeze(-2)) @ U.mH
            Bs = alg.operator(Bs / (n ** (p - 1))[:, None, None, None].to(Bs.dtype))
            assert torch.allclose(pair(A, Bs).real, n, atol=1e-10)
            assert torch.allclose(dv.lp_norm(Bs, q, trace), torch.ones(B, dtype=torch.float64), atol=1e-10)
    # commutative algebra: weighted l^p norms
    comm = _alg((1, 1, 1, 1, 1))
    v = torch.randn(B, comm.C, dtype=torch.complex128)
    Z = comm.operator(v[:, :, None, None])
    wts = torch.rand(B, comm.C, dtype=torch.float64)
    for p in (1.0, 2.5, math.inf):
        ref = v.abs().amax(-1) if p == math.inf else (wts * v.abs() ** p).sum(-1) ** (1 / p)
        assert torch.allclose(dv.lp_norm(Z, p, weights=wts), ref, atol=1e-12)
    assert torch.allclose(dv.lp_norm(Z, 3.0, 'tau_vN'), ((v.abs() ** 3).sum(-1) / 5) ** (1 / 3), atol=1e-12)
    with pytest.raises(ValueError):
        dv.lp_norm(Z, 0.5)


def test_kosaki_norms():
    """||1||_{p,rho} = 1; ||x||_{p,rho} <= ||x||_inf and nondecreasing in p (symmetric embedding); left / right /
    symmetric against numpy; independent of the trace; Beigi's formula ~D_alpha = alpha/(alpha-1) log
    ||sigma^{-1/2} rho sigma^{-1/2}||_{alpha,sigma}."""
    alg = _alg((3, 2, 4))
    torch.manual_seed(19)
    B = 4
    rho, sigma = _random_states(alg, B, 19)
    x = alg.operator(torch.randn(B, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype))
    one = alg.identity(B)
    for p in (1.0, 2.0, 3.5):
        for eta in (0.0, 0.5, 1.0):
            assert torch.allclose(dv.kosaki_norm(one, rho, p, eta), torch.ones(B, dtype=torch.float64), atol=1e-12)
            n = dv.kosaki_norm(x, rho, p, eta)
            for b in range(B):
                r, X = _dense(alg, rho.matrix[b]), _dense(alg, x.matrix[b])
                Y = _fpow(r, eta / p) @ X @ _fpow(r, (1 - eta) / p)
                ref = np.sum(np.linalg.svd(Y, compute_uv=False) ** p) ** (1 / p)
                assert abs(n[b].item() - ref) < 1e-10
    ps = [1.0, 1.5, 2.0, 3.0, 6.0, math.inf]
    ns = torch.stack([dv.kosaki_norm(x, rho, p) for p in ps])
    assert (ns[1:] >= ns[:-1] - 1e-12).all()
    assert torch.allclose(ns[-1], x.operator_norm(), atol=1e-12)
    for trace in ('norm', 'tau_vN'):
        d = _positive(alg, rho.density(trace))
        assert torch.allclose(dv.kosaki_norm(x, d, 2.5, 0.5, trace=trace), dv.kosaki_norm(x, rho, 2.5), atol=1e-12)
    inv = sigma.power(-0.5)
    Gamma = alg.operator((inv @ rho @ inv).matrix)
    for a in (1.5, 2.0, 4.0):
        assert torch.allclose(a / (a - 1) * torch.log(dv.kosaki_norm(Gamma, sigma, a)),
                              dv.sandwiched_renyi_divergence(rho, sigma, a), atol=1e-10)


# ----------------------------------------------------------------------------------------------------------------
# gradients, synchronisation, errors
# ----------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('complex_', [True, False])
@pytest.mark.parametrize('kind', KINDS)
def test_gradcheck(complex_, kind):
    """Gradients in rho and sigma (Daleckii-Krein through the spectral calculus) against finite differences, double
    precision, non-degenerate spectra; also for D_max, Umegaki, the L^p and Kosaki norms and the Petz map."""
    alg = _alg((2, 3), complex_=complex_)
    torch.manual_seed(20)
    sigma = st.random_density_matrix(alg, 2)
    G0 = torch.randn(2, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype) * alg._const('block_mask')

    def state(G):
        m = G @ G.conj().transpose(-2, -1) + 0.05 * alg._const('eye')
        m = m / torch.diagonal(m, dim1=-2, dim2=-1).sum((-2, -1)).real[:, None, None, None].to(m.dtype)
        return DensityMatrix(alg, matrix=m, validate=False)
    G = G0.clone().requires_grad_(True)
    for a in (0.6, 1.7):
        assert torch.autograd.gradcheck(lambda g: DIV[kind](state(g), sigma, a), (G,), atol=1e-7)
        assert torch.autograd.gradcheck(lambda g: DIV[kind](sigma, state(g), a), (G,), atol=1e-7)
    if kind == 'petz':
        assert torch.autograd.gradcheck(lambda g: dv.umegaki_relative_entropy(state(g), sigma), (G,), atol=1e-7)
        assert torch.autograd.gradcheck(lambda g: dv.max_relative_entropy(state(g), sigma), (G,), atol=1e-7)
        X = alg.operator(torch.randn(2, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype))
        assert torch.autograd.gradcheck(lambda g: dv.lp_norm(alg.operator(g), 3.0, 'norm'), (G,), atol=1e-7)
        assert torch.autograd.gradcheck(lambda g: dv.kosaki_norm(X, state(g), 2.5), (G,), atol=1e-7)
        Phi = ch.random_channel(alg, 2, 2)
        assert torch.autograd.gradcheck(
            lambda g: dv.PetzRecovery(state(g), Phi)(Phi(sigma)).matrix.real.sum((-3, -2, -1)), (G,), atol=1e-7)


_NAMES = ('item', 'tolist', '__bool__', '__float__', '__int__', '__index__', 'cpu', 'numpy')


def _count_host_reads(fn):
    calls = []
    saved = {n: getattr(torch.Tensor, n) for n in _NAMES}

    def wrap(name):
        f = saved[name]

        def g(self, *a, **k):
            calls.append(name)
            return f(self, *a, **k)
        return g
    for n in _NAMES:
        setattr(torch.Tensor, n, wrap(n))
    try:
        fn()
    finally:
        for n, f in saved.items():
            setattr(torch.Tensor, n, f)
    return calls


def test_no_host_synchronisation():
    """The computational paths read no values back to the host (states with known positivity)."""
    alg = _alg((6, 6, 4, 6))
    torch.manual_seed(21)
    B = 8
    rho, sigma = st.random_density_matrix(alg, B), st.random_density_matrix(alg, B)
    Phi = ch.random_channel(alg, 3, B)
    Phi(rho)                                                         # one-off trace-preservation check
    A = torch.linspace(0.3, 3.0, B, dtype=torch.float64)
    x = alg.operator(torch.randn(B, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype))
    cases = {
        'petz': lambda: dv.petz_renyi_divergence(rho, sigma, 1.5),
        'sandwiched, alpha per sample': lambda: dv.sandwiched_renyi_divergence(rho, sigma, A),
        'umegaki, D_max, D_0': lambda: (dv.umegaki_relative_entropy(rho, sigma), dv.max_relative_entropy(rho, sigma),
                                        dv.min_relative_entropy(rho, sigma)),
        'sector contributions': lambda: dv.sector_contributions(rho, sigma, 0.7, 'sandwiched'),
        'Petz map': lambda: dv.PetzRecovery(sigma, Phi)(Phi(rho)).matrix,
        'L^p and Kosaki norms': lambda: (dv.lp_norm(x, 3.0, 'norm'), dv.kosaki_norm(x, rho, 2.0)),
    }
    for name, fn in cases.items():
        fn()
        calls = _count_host_reads(fn)
        assert not calls, f"{name}: {len(calls)} host reads ({sorted(set(calls))})"


def test_errors():
    alg = _alg((2, 3))
    rho, sigma = _random_states(alg, 1, 22)
    for bad in (0.0, 1.0, -0.5, math.inf, math.nan):
        with pytest.raises(ValueError):
            dv.petz_renyi_divergence(rho, sigma, bad)
    with pytest.raises(ValueError):
        dv.sandwiched_renyi_divergence(rho, st.random_density_matrix(_alg((3, 2)), 1), 0.5)
    with pytest.raises(TypeError):
        dv.petz_renyi_divergence(rho.matrix, sigma, 0.5)
    with pytest.raises(RuntimeError):
        dv.petz_renyi_divergence(alg.operator(-rho.matrix, is_self_adjoint=True), sigma, 0.5)
    with pytest.raises(ValueError):
        dv.kosaki_norm(alg.identity(), rho, 2.0, eta=1.5)
    with pytest.raises(ValueError):
        dv.renyi_divergence(rho, sigma, 0.5, kind='other')
