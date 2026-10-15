"""
Time steppers of the Lindblad equation (Section 8.7 of the lecture notes): the RK4 step is the polynomial T_4(hL)
and is not positive (Proposition 8.15); the splitting of Theorem 8.16 (lindblad_evolve(method='kraus')) and the exact
step exp(hL) (method='exact') are channels for every step size, including jumps between sectors. References are built
independently with numpy/scipy on M_D (Proposition 8.8(1)). CPU, double precision, fixed seeds.
"""
import math

import numpy as np
import pytest
import scipy.linalg as sla
import torch

from torch_vn_algebra import TypeIAlgebra, DensityMatrix, InterSectorChannel, channels, dynamics, states

C128, F64 = torch.complex128, torch.float64


def T4(z):
    return 1 + z + z ** 2 / 2 + z ** 3 / 6 + z ** 4 / 24


def qubit(omega, gamma):
    q = TypeIAlgebra([2], [2], precision='double')                 # basis |g>, |e>
    H = q.from_blocks([torch.tensor([[0.0, 0.0], [0.0, omega]], dtype=F64)])
    sm = q.from_blocks([torch.tensor([[0.0, 1.0], [0.0, 0.0]], dtype=F64)])   # |g><e|
    return q, H, [sm], [gamma]


def fock_loss():
    """Two bosonic modes, N <= 2: M_1 (+) M_2 (+) M_3, Bose-Hubbard dimer, loss of either mode (Section 8.3)."""
    fock = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double', charges=[0, 1, 2])
    r2 = math.sqrt(2)
    a1 = {(0, 1): [[1.0, 0.0]], (1, 2): [[r2, 0.0, 0.0], [0.0, 1.0, 0.0]]}
    a2 = {(0, 1): [[0.0, 1.0]], (1, 2): [[0.0, 1.0, 0.0], [0.0, 0.0, r2]]}
    blocks = lambda a: {dc: torch.tensor(K, dtype=F64) for dc, K in a.items()}
    loss = InterSectorChannel.from_blocks(fock, fock, {dc: [blocks(a1)[dc], blocks(a2)[dc]] for dc in a1})
    H = fock.from_blocks([torch.zeros(1, 1, dtype=F64), torch.tensor([[0.0, 1.0], [1.0, 0.0]], dtype=F64),
                          torch.tensor([[0.8, r2, 0.0], [r2, 0.0, r2], [0.0, r2, 0.8]], dtype=F64)])
    return fock, H, loss


def mixed_problem(seed=0, batch=1):
    """H, a jump inside the sectors and a random jump between the sectors on M_1 (+) M_2 (+) M_3."""
    torch.manual_seed(seed)
    alg = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')
    H = alg.operator_from_eigenvalues(lambda d: torch.randn(batch, d, dtype=F64), batch_size=batch,
                                      force_self_adjoint=True)
    L = alg.operator(torch.randn(batch, 3, 3, 3, dtype=C128) / 2)
    Phi = channels.random_inter_sector_channel(alg, alg, kraus_rank=2)
    return alg, H, [L, Phi], [1.0, 0.7]


def embed(alg, X):
    """(+)_c X_c as a D x D numpy matrix."""
    return sla.block_diag(*[X[c, :k, :k] for c, k in enumerate(alg.k_factors)])


def full_jumps(alg, jumps, rates, b=0):
    """The jumps on M_D of Proposition 8.8(1): sqrt(g) L_j and sqrt(g) J^{dc}_l, one block each."""
    ks, off = alg.k_factors, np.cumsum([0] + list(alg.k_factors))
    out = []
    for g, L in zip(rates, jumps):
        if isinstance(L, InterSectorChannel):
            K = L.kraus[min(b, L.kraus.shape[0] - 1)].numpy()
            for l in range(K.shape[0]):
                for d, kd in enumerate(ks):
                    for c, kc in enumerate(ks):
                        J = np.zeros((off[-1], off[-1]), dtype=complex)
                        J[off[d]:off[d] + kd, off[c]:off[c] + kc] = K[l, d, c, :kd, :kc]
                        out.append(math.sqrt(g) * J)
        else:
            out.append(math.sqrt(g) * embed(alg, L.matrix[min(b, L.matrix.shape[0] - 1)].numpy()))
    return out


def full_generator(Hd, Ls):
    """Matrix of the GKSL generator on M_D acting on row-major vec."""
    D = Hd.shape[0]
    I = np.eye(D)
    G = -1j * (np.kron(Hd, I) - np.kron(I, Hd.T))
    for L in Ls:
        LL = L.conj().T @ L
        G = G + np.kron(L, L.conj()) - 0.5 * (np.kron(LL, I) + np.kron(I, LL.T))
    return G


def step_superoperator(alg, H, jumps, rates, h, method):
    """The step map rho -> R_h(rho) of lindblad_evolve(substeps=1) as a dict {(d, c): array (k_d^2, k_c^2)} on
    row-major vec, reconstructed by linearity from its values on Hermitian matrices."""
    ks = alg.k_factors
    basis = []
    for c, k in enumerate(ks):
        for i in range(k):
            for j in range(k):
                E = torch.zeros(alg.C, alg.k_max, alg.k_max, dtype=C128)
                E[c, i, j] = 1.0
                basis.append(E)
    E = torch.stack(basis)
    X1 = E + E.conj().transpose(-2, -1)
    X2 = 1j * (E - E.conj().transpose(-2, -1))
    out = []
    for X in (X1, X2):
        out.append(dynamics.lindblad_evolve(alg.operator(X), H, jumps, [0.0, h], rates=rates, substeps=1,
                                            method=method)[-1].matrix)
    R = (0.5 * (out[0] - 1j * out[1])).numpy()                      # R_h(E) by linearity
    S, n = {}, 0
    cols = {}
    for c, k in enumerate(ks):
        cols[c] = range(n, n + k * k)
        n += k * k
    for d, kd in enumerate(ks):
        for c, kc in enumerate(ks):
            S[d, c] = np.stack([R[m, d, :kd, :kd].reshape(-1) for m in cols[c]], axis=1)
    return S


def choi_blocks(alg, S):
    """Choi matrices J(R_dc) = [R_dc(E_ij)] of the sector components (Proposition 6.4: CP iff all are positive)."""
    ks = alg.k_factors
    out = []
    for (d, c), Sdc in S.items():
        kd, kc = ks[d], ks[c]
        out.append(Sdc.reshape(kd, kd, kc, kc).transpose(2, 0, 3, 1).reshape(kc * kd, kc * kd))
    return out


# ----------------------------------------------------------------------------------------------------------------
# Proposition 8.15: RK4 is not positive
# ----------------------------------------------------------------------------------------------------------------
def test_rk4_step_is_T4_of_the_generator():
    """Section 8.7: one RK4 step of rho' = L(rho) is the polynomial T_4(hL) of the generator."""
    q, H, J, g = qubit(1.3, 0.4)
    G = dynamics.lindblad_superoperator(q, H, J, rates=g)[0][0].numpy()
    for h in (0.1, 1.0, 3.0):
        S = step_superoperator(q, H, J, g, h, 'rk4')[0, 0]
        P = sum(np.linalg.matrix_power(h * G, n) / math.factorial(n) for n in range(5))
        assert np.abs(S - P).max() < 1e-13


def test_rk4_qubit_formula_and_criterion():
    """Proposition 8.15(1): R_h(rho) = [[rho_gg + (1-tau) rho_ee, a rho_ge], [conj(a) rho_eg, tau rho_ee]], and R_h is
    CP iff tau <= 1 and |a|^2 <= tau."""
    rng = np.random.default_rng(1)
    for omega, gamma, h in [(1.0, 1.0, 1.0), (2.0, 1.0, 1.0), (0.0, 1.0, 2.0), (0.0, 1.0, 3.0), (3.0, 0.5, 0.2)]:
        q, H, J, g = qubit(omega, gamma)
        tau, a = T4(-h * gamma), T4(h * (1j * omega - gamma / 2))
        S = step_superoperator(q, H, J, g, h, 'rk4')[0, 0]
        X = rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2))
        Y = (S @ X.reshape(-1)).reshape(2, 2)
        ref = np.array([[X[0, 0] + (1 - tau) * X[1, 1], a * X[0, 1]], [np.conj(a) * X[1, 0], tau * X[1, 1]]])
        assert np.abs(Y - ref).max() < 1e-12
        cp = np.linalg.eigvalsh(choi_blocks(q, {(0, 0): S})[0]).min() > -1e-14
        assert cp == (tau <= 1 and abs(a) ** 2 <= tau)


def test_rk4_positivity_polynomial():
    """Proposition 8.15(2),(3): x_* is the real root of x^3 - 4x^2 + 12x - 24; tau - |a|^2 = h^5 (c5 + c6 h + c7 h^2
    + c8 h^3); for omega = gamma it is negative for every h > 0, for omega = 0 positive for h gamma <= x_*."""
    roots = np.roots([1, -4, 12, -24])
    xs = roots[np.abs(roots.imag) < 1e-12].real
    assert len(xs) == 1 and abs(xs[0] - 2.7852935634) < 1e-9
    x_star = xs[0]
    assert T4(-x_star + 1e-9) < 1 < T4(-x_star - 1e-9)
    rng = np.random.default_rng(2)
    for _ in range(20):
        omega, gamma, h = rng.uniform(0, 3), rng.uniform(0.1, 3), rng.uniform(0.01, 2)
        u = gamma ** 2 + 4 * omega ** 2
        c = [-gamma * u * (4 * omega ** 2 - 3 * gamma ** 2) / 384, u ** 2 * (4 * omega ** 2 - 5 * gamma ** 2) / 4608,
             gamma * u ** 3 / 9216, -u ** 4 / 147456]
        lhs = T4(-h * gamma) - abs(T4(h * (1j * omega - gamma / 2))) ** 2
        rhs = h ** 5 * sum(ci * h ** i for i, ci in enumerate(c))
        assert abs(lhs - rhs) < 1e-12 * max(1.0, abs(lhs))
    x = np.geomspace(0.1, 100, 400)                     # x = h gamma (below 0.1 the difference is below round-off)
    assert (T4(-x) - np.abs(T4(x * (1j - 0.5))) ** 2 < 0).all()                  # omega = gamma
    x = np.linspace(0.1, x_star, 200)
    assert (T4(-x) - T4(-x / 2) ** 2 > 0).all()                                   # omega = 0
    q0 = lambda x: 1 / 128 - 5 * x / 4608 + x ** 2 / 9216 - x ** 3 / 147456      # omega = 0: tau - |a|^2 = x^5 q0(x)
    assert q0(3.0) > 0 and 1 / 128 - 15 / 4608 - 27 / 147456 > 0
    # y^3 - 16 y^2 + 32 y + 384 > 0 for y > 0 (the case omega = gamma of the proof)
    y = np.linspace(0, 100, 100001)
    assert (y ** 3 - 16 * y ** 2 + 32 * y + 384 > 0).all()


def test_rk4_counterexample_negative_eigenvalue():
    """Proposition 8.15(3), omega = gamma = 1, h = 1: one RK4 step of lindblad_evolve maps a pure state of the qubit,
    and the maximally entangled state of the qubit and an idle ancilla, to matrices with a negative eigenvalue; the
    splitting and the exact step do not."""
    q, H, J, g = qubit(1.0, 1.0)
    th = 0.2746
    psi = torch.tensor([math.cos(th / 2), math.sin(th / 2)], dtype=C128)
    rho = DensityMatrix(q, matrix=torch.outer(psi, psi.conj())[None, None])
    out = {m: dynamics.lindblad_evolve(rho, H, J, [0.0, 1.0], rates=g, substeps=1, method=m)[-1]
           for m in dynamics.LINDBLAD_METHODS}
    lam = torch.linalg.eigvalsh(out['rk4'].matrix[0, 0])
    tau, a = T4(-1.0), T4(1j - 0.5)
    c, s = math.cos(th / 2), math.sin(th / 2)
    ref = np.linalg.eigvalsh(np.array([[c * c + (1 - tau) * s * s, a * c * s], [np.conj(a) * c * s, tau * s * s]]))
    assert np.abs(lam.numpy() - ref).max() < 1e-14
    assert lam.min() < -8e-5
    for m in ('kraus', 'exact'):
        assert torch.linalg.eigvalsh(out[m].matrix[0, 0]).min() > 0
    # with an idle ancilla: M_4 = M_2 (x) M_2, generator L (x) id, initial state the maximally entangled state
    q4 = TypeIAlgebra([4], [4], precision='double')
    I2 = torch.eye(2, dtype=F64)
    H4 = q4.from_blocks([torch.kron(H.matrix[0, 0], I2)])
    J4 = [q4.from_blocks([torch.kron(J[0].matrix[0, 0], I2)])]
    omega_vec = torch.tensor([1.0, 0.0, 0.0, 1.0], dtype=C128) / math.sqrt(2)
    rho4 = DensityMatrix(q4, matrix=torch.outer(omega_vec, omega_vec)[None, None])
    choi = dynamics.lindblad_evolve(rho4, H4, J4, [0.0, 1.0], rates=g, substeps=1)[-1]
    lam4 = torch.linalg.eigvalsh(choi.matrix[0, 0])
    ref4 = 0.5 * ((1 + tau) - math.sqrt((1 - tau) ** 2 + 4 * abs(a) ** 2)) / 2
    assert abs(lam4.min().item() - ref4) < 1e-14 and ref4 < -3e-3
    kr = dynamics.lindblad_evolve(rho4, H4, J4, [0.0, 1.0], rates=g, substeps=1, method='kraus')[-1]
    assert torch.linalg.eigvalsh(kr.matrix[0, 0]).min() > -1e-15


def test_rk4_large_step_negative_population():
    """Proposition 8.15(2): for h gamma > x_* one RK4 step maps |e><e| to diag(1 - tau, tau) with 1 - tau < 0."""
    q, H, J, g = qubit(0.0, 1.0)
    rho = DensityMatrix(q, matrix=torch.tensor([[[[0.0, 0.0], [0.0, 1.0]]]], dtype=C128))
    for h in (2.7, 2.9, 3.0):
        p = dynamics.lindblad_evolve(rho, None, J, [0.0, h], rates=g, substeps=1)[-1].matrix[0, 0].diagonal().real
        assert abs(p[0].item() - (1 - T4(-h))) < 1e-14
        assert (p[0].item() < 0) == (h > 2.7852935634)


# ----------------------------------------------------------------------------------------------------------------
# Theorem 8.16: the splitting is a channel for every step size; the exact step
# ----------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("method", ['kraus', 'exact'])
@pytest.mark.parametrize("h", [1e-3, 0.1, 1.0, 10.0, 1e3])
def test_step_is_a_channel_for_every_h(method, h):
    """Theorem 8.16(1) and Proposition 8.8(1): the steps of method='kraus' and 'exact' are CP (all Choi matrices of
    the sector components positive, Proposition 6.4) and TP, with jumps inside and between the sectors."""
    for alg, H, J, g in [mixed_problem(3), fock_loss()[:2] + ([fock_loss()[2]], [0.5])]:
        S = step_superoperator(alg, H, J, g, h, method)
        for Jc in choi_blocks(alg, S):
            assert np.abs(Jc - Jc.conj().T).max() < 1e-12
            assert np.linalg.eigvalsh(Jc).min() > -1e-12 * max(1.0, np.abs(Jc).max())
        ks = alg.k_factors
        for c, kc in enumerate(ks):                                 # Tr sum_d R_dc(X) = Tr X
            tr = sum(np.trace(S[d, c].reshape(kd, kd, kc * kc), axis1=0, axis2=1) for d, kd in enumerate(ks))
            assert np.abs(tr - np.eye(kc).reshape(-1)).max() < 1e-12


def test_rk4_step_not_cp_with_sectors():
    """Proposition 8.15 in an algebra with sectors: for the particle-loss example the RK4 step with h = 3 is not CP
    (a Choi matrix of a sector component has a negative eigenvalue), the splitting is."""
    fock, H, loss = fock_loss()
    S = step_superoperator(fock, H, [loss], [1.0], 3.0, 'rk4')
    assert min(np.linalg.eigvalsh(Jc).min() for Jc in choi_blocks(fock, S)) < -0.1


def test_kraus_step_matches_numpy_reference():
    """Theorem 8.16: one step of method='kraus' equals Psi_h(M^{-1/2} rho M^{-1/2}) computed with numpy on M_D from
    the jumps of Proposition 8.8(1), E = exp(-i h/2 H_eff), Psi_h = S o (id + h Phi + h^2/2 Phi^2) o S, M = Psi_h^*(1)."""
    alg, H, J, g = mixed_problem(4)
    rho = states.random_density_matrix(alg)
    Hd = embed(alg, H.matrix[0].numpy())
    Ls = full_jumps(alg, J, g)
    A = sum(L.conj().T @ L for L in Ls)
    for h in (0.05, 0.7, 5.0):
        E = sla.expm(-0.5j * h * (Hd - 0.5j * A))
        phi = lambda X: sum(L @ X @ L.conj().T for L in Ls)
        phid = lambda Z: sum(L.conj().T @ Z @ L for L in Ls)
        Psi = lambda X: E @ (lambda Y: Y + h * phi(Y) + h * h / 2 * phi(phi(Y)))(E @ X @ E.conj().T) @ E.conj().T
        W = E.conj().T @ E
        M = E.conj().T @ (W + h * phid(W) + h * h / 2 * phid(phid(W))) @ E
        w, V = np.linalg.eigh(M)
        assert w.min() > 0
        N = (V / np.sqrt(w)) @ V.conj().T
        ref = Psi(N @ embed(alg, rho.matrix[0].numpy()) @ N)
        out = dynamics.lindblad_evolve(rho, H, J, [0.0, h], rates=g, substeps=1, method='kraus')[-1]
        assert np.abs(embed(alg, out.matrix[0].numpy()) - ref).max() < 1e-12


def test_exact_matches_expm_on_MD():
    """Proposition 8.8(1): method='exact' equals exp(tL) of the GKSL generator on M_D with the jumps sqrt(g) J^{dc}_l,
    computed with scipy, at every time of the grid."""
    alg, H, J, g = mixed_problem(5)
    rho = states.random_density_matrix(alg)
    Hd = embed(alg, H.matrix[0].numpy())
    G = full_generator(Hd, full_jumps(alg, J, g))
    ts = [0.0, 0.3, 1.0, 2.5]
    traj = dynamics.lindblad_evolve(rho, H, J, ts, rates=g, method='exact')
    X0 = embed(alg, rho.matrix[0].numpy()).reshape(-1)
    for t, r in zip(ts, traj):
        ref = (sla.expm(t * G) @ X0).reshape(Hd.shape)
        assert np.abs(embed(alg, r.matrix[0].numpy()) - ref).max() < 1e-12
        assert np.abs(ref * (1 - embed(alg, np.ones((alg.C, alg.k_max, alg.k_max))))).max() < 1e-12   # stays in M


def test_kraus_without_jumps_is_the_unitary_evolution():
    """Theorem 8.16 with Phi = 0: E_h = exp(-ihH/2), M_h = 1, and every step is exactly U(h) . U(h)^*."""
    torch.manual_seed(6)
    alg = TypeIAlgebra([2, 3], [2, 3], precision='double')
    H = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=F64), force_self_adjoint=True)
    rho = states.random_density_matrix(alg)
    ts = [0.0, 0.4, 2.0]
    for m in ('kraus', 'exact'):
        traj = dynamics.lindblad_evolve(rho, H, [], ts, substeps=3, method=m)
        for a, b in zip(traj, dynamics.von_neumann(rho, H, ts)):
            assert (a.matrix - b.matrix).abs().max() < 1e-12


# ----------------------------------------------------------------------------------------------------------------
# order of accuracy
# ----------------------------------------------------------------------------------------------------------------
def _errors(method, subs, alg, H, J, g, rho, t=1.0):
    ref = dynamics.lindblad_evolve(rho, H, J, [0.0, t], rates=g, method='exact')[-1].matrix
    return [(dynamics.lindblad_evolve(rho, H, J, [0.0, t], rates=g, substeps=n, method=method)[-1].matrix - ref)
            .abs().max().item() for n in subs]


def test_kraus_is_second_order_and_rk4_fourth():
    """Theorem 8.16(3): the global error of the splitting is O(h^2) (error ratio ~4 when h is halved), against the
    exact solution; RK4 has ratio ~16 (fourth order)."""
    alg, H, J, g = mixed_problem(7)
    rho = states.random_density_matrix(alg)
    e = _errors('kraus', [8, 16, 32, 64], alg, H, J, g, rho)
    ratios = [a / b for a, b in zip(e[:-1], e[1:])]
    assert all(3.5 < r < 4.5 for r in ratios), ratios
    e4 = _errors('rk4', [4, 8, 16], alg, H, J, g, rho)
    ratios4 = [a / b for a, b in zip(e4[:-1], e4[1:])]
    assert all(12 < r < 22 for r in ratios4), ratios4


def test_kraus_local_error_is_cubic():
    """Theorem 8.16(2): the local error ||Psi~_h - exp(hL)|| is O(h^3): ratio ~8 when h is halved."""
    alg, H, J, g = mixed_problem(8)
    errs = []
    for h in (0.08, 0.04, 0.02, 0.01):
        Sk = step_superoperator(alg, H, J, g, h, 'kraus')
        Se = step_superoperator(alg, H, J, g, h, 'exact')
        errs.append(max(np.abs(Sk[key] - Se[key]).max() for key in Sk))
    ratios = [a / b for a, b in zip(errs[:-1], errs[1:])]
    assert all(7.0 < r < 9.0 for r in ratios), ratios


def test_kraus_global_bound_long_run():
    """Theorem 8.16(3) and (1): over a long run with large steps the splitting keeps a state (positive, trace one)
    and its trace distance to the exact solution is O(h^2)."""
    fock, H, loss = fock_loss()
    torch.manual_seed(9)
    rho = states.random_density_matrix(fock)
    ts = [0.0, 5.0, 20.0]
    ex = dynamics.lindblad_evolve(rho, H, [loss], ts, rates=[0.5], method='exact')
    for n in (1, 4, 16):
        kr = dynamics.lindblad_evolve(rho, H, [loss], ts, rates=[0.5], substeps=n, method='kraus')
        for a, b in zip(kr, ex):
            assert min(w.min().item() for w in a.eigenvalues()) > -1e-14
            assert abs(a.matrix.diagonal(dim1=-2, dim2=-1).sum().item() - 1) < 1e-13
            assert states.trace_distance(a, b).item() < 0.5 * 20.0 * (5.0 / n) ** 2


# ----------------------------------------------------------------------------------------------------------------
# library behaviour
# ----------------------------------------------------------------------------------------------------------------
def test_default_method_is_rk4_and_validation():
    """lindblad_evolve: the default is method='rk4' (unchanged results); an unknown method or decreasing times for
    'kraus'/'exact' raise."""
    alg, H, J, g = mixed_problem(10)
    rho = states.random_density_matrix(alg)
    a = dynamics.lindblad_evolve(rho, H, J, [0.0, 0.5, 1.0], rates=g, substeps=5)
    b = dynamics.lindblad_evolve(rho, H, J, [0.0, 0.5, 1.0], rates=g, substeps=5, method='rk4')
    c = dynamics.rk4(dynamics.lindblad_rhs(H, J, g), rho.matrix.to(C128), [0.0, 0.5, 1.0], substeps=5)
    for x, y, z in zip(a, b, c):
        assert torch.equal(x.matrix, y.matrix)
        assert torch.equal(x.matrix, 0.5 * (z + z.conj().transpose(-2, -1)))
    with pytest.raises(ValueError):
        dynamics.lindblad_evolve(rho, H, J, [0.0, 1.0], rates=g, method='euler')
    for m in ('kraus', 'exact'):
        with pytest.raises(ValueError):
            dynamics.lindblad_evolve(rho, H, J, [1.0, 0.0], rates=g, method=m)


def test_batched_and_real_and_time_dependent():
    """method='kraus'/'exact': a batch of Hamiltonians equals the runs one by one; a real algebra with real jumps
    stays real; a time-dependent H (generator frozen at the step midpoint) converges to RK4 at second order."""
    alg, H, J, g = mixed_problem(11, batch=3)
    rho = states.random_density_matrix(alg)
    for m in ('kraus', 'exact'):
        allb = dynamics.lindblad_evolve(rho, H, J, [0.0, 0.7], rates=g, substeps=4, method=m)[-1].matrix
        for b in range(3):
            Hb = alg.operator(H.matrix[b:b + 1])
            Lb = alg.operator(J[0].matrix[b:b + 1])
            one = dynamics.lindblad_evolve(rho, Hb, [Lb, J[1]], [0.0, 0.7], rates=g, substeps=4, method=m)[-1]
            assert (allb[b] - one.matrix[0]).abs().max() < 1e-12
    real = TypeIAlgebra([2, 3], [2, 3], precision='double', complex_valued=False)
    torch.manual_seed(12)
    Lr = real.operator(torch.randn(1, 2, 3, 3, dtype=F64))
    rr = states.random_density_matrix(real)
    outs = {m: dynamics.lindblad_evolve(rr, None, [Lr], [0.0, 1.0], substeps=50, method=m)[-1].matrix
            for m in dynamics.LINDBLAD_METHODS}
    assert all(not torch.is_complex(x) for x in outs.values())
    assert (outs['kraus'] - outs['exact']).abs().max() < 1e-4 and (outs['rk4'] - outs['exact']).abs().max() < 1e-6
    alg, H, J, g = mixed_problem(13)
    V = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=F64), force_self_adjoint=True)
    Ht = lambda s: H.matrix + math.sin(2 * s) * V.matrix
    rho = states.random_density_matrix(alg)
    ref = dynamics.lindblad_evolve(rho, Ht, J, [0.0, 1.0], rates=g, substeps=400)[-1].matrix
    for m in ('kraus', 'exact'):
        e = [(dynamics.lindblad_evolve(rho, Ht, J, [0.0, 1.0], rates=g, substeps=n, method=m)[-1].matrix - ref)
             .abs().max().item() for n in (8, 16, 32)]
        assert all(3.3 < a / b < 4.7 for a, b in zip(e[:-1], e[1:])), e
