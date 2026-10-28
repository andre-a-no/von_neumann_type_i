"""
The mathematical statements and exact values of the article (man/article/vn_algebra_pytorch.tex), checked through the public API
of torch_vn_algebra: CPU, double precision, fixed seeds. Every test names the statement it checks by the section of
the article ("Sec. 2", ...) or by its numbered object ("Proposition 1") or by the program it refers to ("ex3_channels.py"). Distributional claims are
checked with z-scores; exact claims to round-off.
"""
import math
from fractions import Fraction
from math import comb

import numpy as np
import pytest
import torch

import torch_vn_algebra as tv
from torch_vn_algebra import (TypeIAlgebra, DensityMatrix, InterSectorChannel, SpinChain, SU2Chain, channels, states,
                              dynamics, krylov, composite, inclusions as inc, decompose as dc, divergences as dv,
                              optimize as opt, su2)

pytestmark = [pytest.mark.filterwarnings("ignore::torch_vn_algebra.cost.CostWarning"),
              pytest.mark.filterwarnings("ignore:sector sizes are very uneven")]

DT = torch.complex128
LAYOUTS = [[1], [3], [2, 3, 5], [1, 1], [4, 1, 2]]
TRACES = ('blunt', 'norm', 'tau_vN')


@pytest.fixture(autouse=True)
def _seed():
    torch.manual_seed(1234)
    np.random.seed(1234)


# ----------------------------------------------------------------------------------------------------------------
# helpers (dense references built from the blocks)
# ----------------------------------------------------------------------------------------------------------------

def alg_of(ks, complex_valued=True, ns=None):
    return TypeIAlgebra(list(ns or ks), list(ks), complex_valued=complex_valued, precision='double')


def from_blocks(alg, blocks):
    """Operator with batch B from a list of (B, k_c, k_c) tensors."""
    B = blocks[0].shape[0]
    M = torch.zeros(B, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype)
    for c, (k, X) in enumerate(zip(alg.k_factors, blocks)):
        M[:, c, :k, :k] = X
    return alg.operator(M)


def blocks_of(op):
    return [op.matrix[:, c, :k, :k] for c, k in enumerate(op.algebra.k_factors)]


def rand_blocks(alg, B, kind='general'):
    out = []
    for k in alg.k_factors:
        G = torch.randn(B, k, k, dtype=alg.dtype)
        if kind == 'herm':
            G = G + G.mH
        elif kind == 'pos':
            G = G @ G.mH + 0.05 * torch.eye(k, dtype=alg.dtype)
        out.append(G)
    return out


def rand_op(alg, B, kind='general'):
    return from_blocks(alg, rand_blocks(alg, B, kind))


def weights(alg, trace):
    k = torch.tensor(alg.k_factors, dtype=torch.float64)
    C = sum(1 for x in alg.k_factors if x > 0)
    return {'blunt': torch.ones_like(k), 'norm': 1 / k, 'tau_vN': 1 / (C * k)}[trace]


def tr_w(op, trace):
    return {'blunt': op.Tr_blunt, 'norm': op.Tr_norm, 'tau_vN': op.tau_vN}[trace]()


def ref_trace(op, w):
    return sum(w[c] * torch.diagonal(X, dim1=-2, dim2=-1).sum(-1) for c, X in enumerate(blocks_of(op)))


def mpow(A, p):
    """Power of a Hermitian positive semidefinite matrix (numpy), 0^p = 0."""
    w, V = np.linalg.eigh(A)
    w = np.clip(w, 0, None)
    wp = np.where(w > 1e-14, np.abs(w) ** p, 0.0)
    return (V * wp) @ V.conj().T


def rand_state_np(n, rng, floor=0.0):
    G = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    r = G @ G.conj().T
    r = r / np.trace(r).real
    return (1 - floor) * r + floor * np.eye(n) / n


def state_op(alg, rho_np):
    return DensityMatrix(alg, matrix=torch.as_tensor(rho_np, dtype=alg.dtype)[None, None])


# ================================================================================================================
# Section 2: mathematical setting
# ================================================================================================================

@pytest.mark.parametrize('ks', LAYOUTS, ids=str)
def test_sec2_three_traces_have_the_stated_weights(ks):
    """Sec. 2: Tr_blunt, Tr_norm, tau_vN are tr_w with w_c = 1, 1/k_c, 1/(C k_c); values on 1 are sum k_c, C, 1."""
    alg = alg_of(ks)
    A = rand_op(alg, 6)
    for t in TRACES:
        assert torch.allclose(tr_w(A, t), ref_trace(A, weights(alg, t)).to(DT), atol=1e-12)
    one = alg.identity()
    vals = [tr_w(one, t).real.item() for t in TRACES]
    assert vals == pytest.approx([sum(ks), len(ks), 1.0], abs=1e-12)
    assert torch.allclose(A.Tr_norm(), len(ks) * A.tau_vN(), atol=1e-12)          # Tr_norm = C tau_vN


def test_sec2_corner_of_nominal_factors_and_empty_sectors():
    """Sec. 2: with n_c >= k_c the traces refer to the corner (sum k_c, C, 1); empty sectors carry no weight."""
    alg = TypeIAlgebra([4, 5], [2, 3], precision='double')
    one = alg.identity()
    assert [tr_w(one, t).real.item() for t in TRACES] == pytest.approx([5, 2, 1], abs=1e-12)
    alg0 = TypeIAlgebra([2, 1, 3], [2, 0, 3], precision='double')                 # one empty sector
    one = alg0.identity()
    assert [tr_w(one, t).real.item() for t in TRACES] == pytest.approx([5, 2, 1], abs=1e-12)


def test_sec5_trace_functionals_table():
    """Sec. 5, Table 'traces': tracial, unitarily invariant and faithful with tr(A*A)/||A||_F^2 >= min_c w_c."""
    alg = alg_of([2, 3, 5])
    A, Bo = rand_op(alg, 50), rand_op(alg, 50)
    U = alg.random_unitary_operator(50)
    for t in TRACES:
        assert (tr_w(A @ Bo, t) - tr_w(Bo @ A, t)).abs().max() < 1e-11
        assert (tr_w(U @ A @ U.adjoint(), t) - tr_w(A, t)).abs().max() < 1e-11
        ratio = tr_w(A.adjoint() @ A, t).real / A.frobenius_norm() ** 2
        assert ratio.min() >= weights(alg, t).min() - 1e-12 and ratio.max() <= weights(alg, t).max() + 1e-12
        # the bound is sharp: an operator supported on the sector of smallest weight attains it
        c = int(torch.argmin(weights(alg, t)))
        bl = [torch.zeros(1, k, k, dtype=DT) for k in alg.k_factors]
        bl[c] = torch.randn(1, alg.k_factors[c], alg.k_factors[c], dtype=DT)
        X = from_blocks(alg, bl)
        assert (tr_w(X.adjoint() @ X, t).real / X.frobenius_norm() ** 2).item() == pytest.approx(
            weights(alg, t).min().item(), rel=1e-12)


def test_sec2_michelson_contrast():
    """Sec. 2: Delta(A) = (l_max - l_min)/(l_max + l_min) over all sectors, in [0, 1]; Delta(0) = 0; Delta = 0 iff A = c 1."""
    alg = alg_of([2, 3])
    A = rand_op(alg, 20, 'pos')
    ev = torch.cat([e for e in A.eigenvalues()], dim=-1).real
    ref = (ev.max(-1)[0] - ev.min(-1)[0]) / (ev.max(-1)[0] + ev.min(-1)[0])
    assert torch.allclose(A.michelson_contrast, ref, atol=1e-12)
    assert ((A.michelson_contrast >= 0) & (A.michelson_contrast <= 1)).all()
    assert alg.zero().michelson_contrast.abs().max() == 0
    assert (2.5 * alg.identity()).michelson_contrast.abs().max() < 1e-14
    # scalar in every block but different scalars: not a multiple of 1, contrast (3 - 1)/(3 + 1)
    X = alg.central([1.0, 3.0])
    assert X.michelson_contrast.item() == pytest.approx(0.5, abs=1e-14)
    # a projection: lambda_min = 0, contrast 1
    P = alg.central([1.0, 0.0])
    assert P.michelson_contrast.item() == pytest.approx(1.0, abs=1e-14)


@pytest.mark.parametrize('ks', [[2, 3], [1, 4, 2]], ids=str)
def test_sec2_densities_for_the_three_traces(ks):
    """Sec. 2: the density of omega w.r.t. tr_w is rho_c / w_c; omega(A) and p_c = Tr rho_c do not depend on the trace."""
    alg = alg_of(ks)
    rho = states.random_density_matrix(alg, batch_size=4)
    A = rand_op(alg, 4, 'herm')
    p = torch.stack([torch.diagonal(X, dim1=-2, dim2=-1).sum(-1).real for X in blocks_of(rho)], dim=-1)
    assert torch.allclose(rho.sector_probabilities(), p, atol=1e-12)
    assert torch.allclose(p.sum(-1), torch.ones(4, dtype=torch.float64), atol=1e-12)
    for t in TRACES:
        d = rho.density(t)
        w = weights(alg, t)
        for c, k in enumerate(ks):
            assert torch.allclose(d[:, c, :k, :k], rho.matrix[:, c, :k, :k] / w[c], atol=1e-12)
        assert torch.allclose(tr_w(alg.operator(d) @ A, t), rho.expectation(A).to(DT), atol=1e-12)
        back = DensityMatrix.from_density(alg, d, trace=t)
        assert (back.matrix - rho.matrix).abs().max() < 1e-12


def loss_channel(g, alg=None):
    """Particle loss on two fermionic modes (ex3_channels.py): Kraus sqrt(1 - gN), sqrt(g) a_1, sqrt(g) a_2."""
    alg = alg or alg_of([1, 2, 1])
    s, z = g ** 0.5, torch.zeros
    blocks = {
        (0, 0): [torch.ones(1, 1, dtype=DT)],
        (1, 1): [(1 - g) ** 0.5 * torch.eye(2, dtype=DT)],
        (2, 2): [(1 - 2 * g) ** 0.5 * torch.ones(1, 1, dtype=DT)],
        (0, 1): [z(1, 2, dtype=DT), s * torch.tensor([[1.0, 0.0]], dtype=DT), s * torch.tensor([[0.0, 1.0]], dtype=DT)],
        (1, 2): [z(2, 1, dtype=DT), s * torch.tensor([[0.0], [1.0]], dtype=DT), s * torch.tensor([[-1.0], [0.0]], dtype=DT)],
    }
    return alg, InterSectorChannel.from_blocks(alg, alg, blocks)


@pytest.mark.parametrize('g', [0.0, 0.2, 0.5])
def test_sec4_particle_loss_is_a_pure_death_process(g):
    """Sec. 4 (ex3_channels.py): the loss channel is trace preserving with T = [[1,g,0],[0,1-g,2g],[0,0,1-2g]] for every state."""
    alg, loss = loss_channel(g)
    assert loss.is_trace_preserving()
    T = torch.tensor([[1, g, 0], [0, 1 - g, 2 * g], [0, 0, 1 - 2 * g]], dtype=torch.float64)
    assert torch.allclose(loss.transition_matrix()[0].real, T, atol=1e-13)
    # sum_i K^* K is a multiple of 1_c in every block: T does not depend on the input state
    rho = states.random_density_matrix(alg, batch_size=5)
    assert torch.allclose(loss.transition_matrix(rho).real, T.expand(5, 3, 3), atol=1e-12)
    # and it is the classical part of the channel: p_out = T p_in
    p_out = loss(rho).sector_probabilities()
    assert torch.allclose(p_out, rho.sector_probabilities() @ T.T, atol=1e-12)


def test_sec2_transition_matrix_depends_on_the_state_in_general():
    """Sec. 2: T_dc = Tr Phi_dc(rho_c)/Tr rho_c depends on rho unless sum_i K^dc* K^dc is a multiple of 1_c."""
    a_in, a_out = alg_of([3, 2]), alg_of([2, 1, 2])
    ch = channels.random_inter_sector_channel(a_in, a_out, kraus_rank=2)
    rho = states.random_density_matrix(a_in, batch_size=6)
    T = ch.transition_matrix(rho).real
    assert torch.allclose(T.sum(-2), torch.ones(6, 2, dtype=torch.float64), atol=1e-12)     # columns sum to 1
    # T_dc from the definition: restrict rho to sector c, apply, read sector d
    for c in range(a_in.C):
        bl = [torch.zeros_like(X) for X in blocks_of(rho)]
        bl[c] = blocks_of(rho)[c]
        out = ch(from_blocks(a_in, bl))
        pc = torch.diagonal(bl[c], dim1=-2, dim2=-1).sum(-1).real
        for d in range(a_out.C):
            Od = blocks_of(out)[d]
            assert torch.allclose(T[:, d, c], torch.diagonal(Od, dim1=-2, dim2=-1).sum(-1).real / pc, atol=1e-12)
    assert T.std(0).max() > 1e-3                                                          # it varies with rho


@pytest.mark.parametrize('trace', TRACES)
def test_sec2_dual_of_a_map_between_sectors(trace):
    """Sec. 2: the dual w.r.t. tr_w satisfies tr_w(Phi(rho) A) = tr_w(rho Phi^{*,w}(A)) and depends on w across sectors."""
    a_in, a_out = alg_of([3, 2]), alg_of([2, 1, 2])
    ch = channels.random_inter_sector_channel(a_in, a_out, kraus_rank=2)
    rho, A = rand_op(a_in, 3, 'pos'), rand_op(a_out, 3)
    lhs = tr_w(ch(rho) @ A, trace)
    rhs = tr_w(rho @ ch.dual(trace)(A), trace)
    assert torch.allclose(lhs, rhs, atol=1e-11)
    # dual(w) = (w_out_d / w_in_c) * blunt dual sector by sector
    blunt = ch.dual('blunt')
    w_in, w_out = weights(a_in, trace), weights(a_out, trace)
    for d in range(a_out.C):
        bl = [torch.zeros_like(X) for X in blocks_of(A)]
        bl[d] = blocks_of(A)[d]
        Ad = from_blocks(a_out, bl)
        for c in range(a_in.C):
            assert torch.allclose(blocks_of(ch.dual(trace)(Ad))[c], w_out[d] / w_in[c] * blocks_of(blunt(Ad))[c],
                                  atol=1e-12)


def test_sec2_dual_of_a_sector_preserving_map_is_the_same_for_all_traces():
    """Sec. 2: for sector-preserving maps on one algebra the weights cancel and all three traces give the same dual."""
    alg = alg_of([2, 3])
    ch = channels.random_channel(alg, kraus_rank=3).to_inter_sector()
    A = rand_op(alg, 2)
    duals = [ch.dual(t)(A).matrix for t in TRACES]
    assert torch.allclose(duals[0], duals[1], atol=1e-12) and torch.allclose(duals[0], duals[2], atol=1e-12)


def test_sec4_unital_positive_maps_do_not_increase_the_contrast():
    """Sec. 4: l_min 1 <= A <= l_max 1 is kept by a unital positive map, so Delta(Phi(A)) <= Delta(A)."""
    alg = alg_of([2, 3, 4])
    X = rand_op(alg, 200, 'pos')
    mix = channels.random_mixed_unitary_channel(alg, 3, batch_size=200)
    Y = mix(X)
    assert (Y.michelson_contrast <= X.michelson_contrast + 1e-12).all()
    assert (Y.lambda_min >= X.lambda_min - 1e-12).all() and (Y.lambda_max <= X.lambda_max + 1e-12).all()
    _, loss = loss_channel(0.2)
    Z = rand_op(loss.algebra_in, 200, 'pos')
    assert (loss(Z).michelson_contrast > Z.michelson_contrast + 1e-9).any()             # non-unital: can increase


def test_sec2_fused_product_sizes_and_partial_trace_dual():
    """Sec. 2: n_q = sum_{q_c+q'_d=q} k_c m_d; Tr_2 lands in M_1 and its Tr_blunt-dual is A -> A (x) 1."""
    A1 = TypeIAlgebra([1, 2, 1], [1, 2, 1], precision='double', charges=[0, 1, 2])
    A2 = TypeIAlgebra([2, 3], [2, 3], precision='double', charges=[0, 1])
    F = composite.fused_tensor_product(A1, A2)
    sizes = {}
    for qc, kc in zip([0, 1, 2], [1, 2, 1]):
        for qd, md in zip([0, 1], [2, 3]):
            sizes[qc + qd] = sizes.get(qc + qd, 0) + kc * md
    assert list(F.k_factors) == [sizes[q] for q in sorted(sizes)]
    X = rand_op(F, 3)
    A = rand_op(A1, 3)
    one2 = A2.identity()
    lhs = (X @ composite.kron(A, one2, F)).Tr_blunt()
    rhs = (composite.partial_trace(X, keep=1) @ A).Tr_blunt()
    assert torch.allclose(lhs, rhs, atol=1e-11)


@pytest.mark.parametrize('L', [6, 8])
def test_sec4_saturation_field_of_the_heisenberg_ring(L):
    """Sec. 4 (ex7_spin_chain.py): E_min(N=1) - E(N=0) = h - 2J (one magnon at k=pi), so the ground state is N=0 iff h >= 2J."""
    ch = SpinChain(L, 'periodic')
    assert list(ch.algebra.k_factors) == [comb(L, N) for N in range(L + 1)]
    for h in (0.3, 1.9, 2.1, 3.0):
        w = [blk[0].real for blk in ch.heisenberg(J=1.0, h=h).eigenvalues()]
        assert (w[1].min() - w[0].min()).item() == pytest.approx(h - 2.0, abs=1e-11)
        gs_sector = min(range(L + 1), key=lambda N: w[N].min().item())
        assert (gs_sector == 0) == (h > 2.0)


def test_sec2_lindblad_with_jumps_in_M_conserves_sector_probabilities():
    """Sec. 2: the GKSL equation with H and L_j in M preserves every sector (p_c conserved: superselection)."""
    alg = alg_of([2, 3])
    H = rand_op(alg, 1, 'herm')
    jumps = [rand_op(alg, 1) for _ in range(2)]
    rho0 = states.random_density_matrix(alg)
    E = dynamics.lindblad_channel(alg, H, jumps, t=0.7, rates=[0.3, 0.5])
    rho = E(rho0)
    assert torch.allclose(rho.sector_probabilities(), rho0.sector_probabilities(), atol=1e-12)
    assert (rho.matrix - rho0.matrix).abs().max() > 1e-3


# ================================================================================================================
# Section 3: design (random ensembles, functional calculus, number field and precision)
# ================================================================================================================

def zscore(x, exact):
    x = x.double()
    return ((x.mean() - exact) / (x.std() / math.sqrt(x.numel()))).item()


@pytest.mark.parametrize('cplx', [True, False], ids=['U(n)', 'O(n)'])
@pytest.mark.parametrize('n', [2, 3, 5])
def test_sec3_haar_fourth_moments(cplx, n):
    """Sec. 3/5: E|U_11|^2 = 1/n, E|U_11|^4 = 2/(n(n+1)) on U(n) and 3/(n(n+2)) on O(n) (z-scores)."""
    alg = alg_of([n], complex_valued=cplx)
    U = alg.random_unitary(n, batch_size=40000)
    a = U[:, 0, 0].abs() ** 2
    m4 = 2 / (n * (n + 1)) if cplx else 3 / (n * (n + 2))
    assert abs(zscore(a, 1 / n)) < 5
    assert abs(zscore(a ** 2, m4)) < 5


def test_sec3_su_so_and_circular_ensembles():
    """Sec. 3: SU(n)/SO(n) have det 1; COE S = W^T W is symmetric unitary; CSE S = W^R W is self-dual unitary."""
    n, B = 4, 50
    a = alg_of([n])
    eye = torch.eye(n, dtype=DT)
    S = a.random_unitary(n, 'haar_su', B)
    assert (torch.linalg.det(S) - 1).abs().max() < 1e-12
    O = alg_of([n], complex_valued=False).random_unitary(n, 'haar_su', B)
    assert (torch.linalg.det(O) - 1).abs().max() < 1e-12
    coe = a.random_unitary(n, 'coe', B)
    assert (coe - coe.transpose(-2, -1)).abs().max() < 1e-12
    assert (coe @ coe.mH - eye).abs().max() < 1e-12
    cse = a.random_unitary(n, 'cse', B)
    J = torch.zeros(n, n, dtype=DT)
    J[:n // 2, n // 2:] = torch.eye(n // 2)
    J[n // 2:, :n // 2] = -torch.eye(n // 2)
    dual = J @ cse.transpose(-2, -1) @ J.T
    assert (dual - cse).abs().max() < 1e-12
    assert (cse @ cse.mH - eye).abs().max() < 1e-12


def test_sec3_functional_calculus_levels():
    """Sec. 3: f(A) = V f(Lambda) V^* (self-adjoint), the normal calculus for unitaries, |A| from A^*A, and the contour integral."""
    alg = alg_of([3, 4])
    H = rand_op(alg, 4, 'herm')
    E = H.apply_function(torch.exp)
    for X, Y in zip(blocks_of(E), blocks_of(H)):
        assert torch.allclose(X, torch.linalg.matrix_exp(Y), atol=1e-10, rtol=1e-12)
    P = rand_op(alg, 4, 'pos')
    R = P.sqrt()
    assert (R @ R - P).frobenius_norm().max() < 1e-11
    U = alg.random_unitary_operator(4)
    U3 = U.apply_function(lambda z: z ** 3)
    assert ((U3 - U @ U @ U).frobenius_norm()).max() < 1e-10
    A = rand_op(alg, 4)
    absA = A.abs()
    for X, Y in zip(blocks_of(absA), blocks_of(A)):
        w, V = torch.linalg.eigh(Y.mH @ Y)
        assert torch.allclose(X, (V * w.clamp(min=0).sqrt().unsqueeze(-2)) @ V.mH, atol=1e-10)
        assert torch.allclose(torch.diagonal(X, dim1=-2, dim2=-1).sum(-1).real, torch.linalg.svdvals(Y).sum(-1),
                              atol=1e-10)
    W, Pp = A.polar()
    assert (W @ Pp - A).frobenius_norm().max() < 1e-10
    assert (W.adjoint() @ W - alg.identity()).frobenius_norm().max() < 1e-10
    Ex = A.holomorphic_function(torch.exp)
    for X, Y in zip(blocks_of(Ex), blocks_of(A)):
        assert torch.allclose(X, torch.linalg.matrix_exp(Y), atol=1e-9, rtol=1e-9)


def test_sec3_daleckii_krein_derivative_at_a_degenerate_spectrum():
    """Sec. 3: dF = V (L o V^* dA V) V^* with L_ij = divided differences (f' on ties), finite at degenerate spectra."""
    alg = alg_of([4])
    V = alg.random_unitary(4, batch_size=1)
    lam = torch.tensor([1.0, 1.0, 2.0, 2.0], dtype=torch.float64)            # doubly degenerate twice
    A = (V * lam.to(DT).unsqueeze(-2)) @ V.mH
    dA = torch.randn(1, 4, 4, dtype=DT)
    dA = dA + dA.mH
    Bm = torch.randn(1, 4, 4, dtype=DT)
    M = A[:, None].clone().requires_grad_(True)
    F = alg.operator(M, is_self_adjoint=True).apply_function(torch.exp)
    val = (F.matrix[:, 0] * Bm.transpose(-2, -1)).sum().real
    val.backward()
    autograd = (M.grad[:, 0].conj() * dA).sum().real.item()
    # Daleckii-Krein reference
    f = torch.exp(lam)
    Lm = torch.empty(4, 4, dtype=torch.float64)
    for i in range(4):
        for j in range(4):
            Lm[i, j] = f[i] if abs(lam[i] - lam[j]) < 1e-12 else (f[i] - f[j]) / (lam[i] - lam[j])
    dF = V @ (Lm.to(DT) * (V.mH @ dA @ V)) @ V.mH
    ref = (dF * Bm.transpose(-2, -1)).sum().real.item()
    assert math.isfinite(autograd)
    assert autograd == pytest.approx(ref, rel=1e-9, abs=1e-10)
    # and the derivative of the exponential itself (finite differences of matrix_exp)
    h = 1e-6
    fd = ((torch.linalg.matrix_exp(A + h * dA) - torch.linalg.matrix_exp(A - h * dA)) / (2 * h))
    assert (fd - dF).abs().max() < 1e-6


def test_sec3_storage_types_and_bytes():
    """Sec. 3: four storage types with 4, 8, 8, 16 bytes per entry; like() changes the numerics, cast() moves operators."""
    alg = TypeIAlgebra([2, 3], [2, 3], layout='padded')
    packed = TypeIAlgebra([2, 3], [2, 3], layout='packed')
    expect = {(False, 'single'): (torch.float32, 4), (False, 'double'): (torch.float64, 8),
              (True, 'single'): (torch.complex64, 8), (True, 'double'): (torch.complex128, 16)}
    for (cv, pr), (dt, nbytes) in expect.items():
        a2 = alg.like(complex_valued=cv, precision=pr)
        assert a2.dtype == dt and list(a2.k_factors) == [2, 3]
        assert a2.bytes_per_operator(7) == 7 * 2 * 3 * 3 * nbytes                      # C k_max^2 entries
        assert packed.like(complex_valued=cv, precision=pr).bytes_per_operator(7) == 7 * (4 + 9) * nbytes
    a64 = alg.like(complex_valued=True, precision='double')
    X = rand_op(alg.like(complex_valued=False, precision='double'), 2, 'herm')
    Y = X.cast(a64)
    assert Y.matrix.dtype == torch.complex128 and torch.allclose(Y.matrix.real, X.matrix)


def test_sec3_real_problems_have_real_extremisers():
    """Sec. 3: for real X, Y >= 0, sup_U |Tr(XUY)| = ||YX||_1 is attained by the real orthogonal polar factor of YX."""
    alg = alg_of([5], complex_valued=False)
    X, Y = rand_op(alg, 8, 'pos'), rand_op(alg, 8, 'pos')
    W, _ = (Y @ X).polar()                                     # YX = W |YX|, W orthogonal
    assert W.matrix.dtype == torch.float64
    assert (W.adjoint() @ W - alg.identity()).frobenius_norm().max() < 1e-10
    val = (X @ W.adjoint() @ Y).trace.abs()
    assert torch.allclose(val, (Y @ X).trace_norm(), rtol=1e-10)
    U = alg.random_unitary_operator(8)
    assert ((X @ U @ Y).trace.abs() <= (Y @ X).trace_norm() + 1e-10).all()


def test_sec3_real_algebra_refuses_unitary_dynamics():
    """Sec. 3: complex matrices are required for e^{-iHt}; the library refuses unitary dynamics in a real algebra."""
    alg = alg_of([2], complex_valued=False)
    with pytest.raises(ValueError):
        dynamics.von_neumann(alg.identity(), alg.identity(), [1.0])


def test_sec3_eigenvalues_without_eigenvectors():
    """Sec. 3: Operator.eigenvalues() gives the spectra of eigh() block by block."""
    alg = alg_of([2, 3, 5])
    H = rand_op(alg, 3, 'herm')
    for w, (w2, _) in zip(H.eigenvalues(), H.eigh()):
        assert torch.allclose(w.real, w2.real, atol=1e-12)


# ================================================================================================================
# Section 4: usage examples
# ================================================================================================================

def test_sec4_lueders_conditioning():
    """Sec. 4 (ex2_states_measurement.py): rho|P = P rho P / Tr(rho P), with probability Tr(rho P); conditioning on a block moves p into it."""
    alg = alg_of([2, 3])
    rho = states.random_density_matrix(alg, batch_size=5)
    P = from_blocks(alg, [torch.diag(torch.tensor([1.0, 0.0], dtype=DT))[None], torch.zeros(1, 3, 3, dtype=DT)])
    post, prob = rho.condition_on(P)
    ref_p = (rho @ P).Tr_blunt().real
    assert torch.allclose(prob, ref_p, atol=1e-12)
    assert torch.allclose(post.matrix, (P @ rho @ P).matrix / ref_p[:, None, None, None], atol=1e-12)
    assert torch.allclose(post.sector_probabilities(), torch.tensor([[1.0, 0.0]]).expand(5, 2).double(), atol=1e-12)


def test_sec4_spontaneous_emission_exact_solution():
    """Sec. 4 (ex4_dynamics.py): P_e(t) = P_e(0) e^{-gamma t}, |rho_eg(t)| = |rho_eg(0)| e^{-gamma t/2} (exp of the generator)."""
    alg = alg_of([2])
    omega, gamma = 1.0, 0.3
    H = from_blocks(alg, [torch.diag(torch.tensor([0.5 * omega, -0.5 * omega], dtype=DT))[None]])
    L = from_blocks(alg, [torch.tensor([[0.0, 0.0], [1.0, 0.0]], dtype=DT)[None]])
    psi = torch.tensor([0.6, 0.8], dtype=DT)
    rho0 = DensityMatrix(alg, matrix=torch.outer(psi, psi.conj())[None, None])
    for t in (0.5, 2.0, 5.0):
        m = dynamics.lindblad_channel(alg, H, [L], t=t, rates=[gamma])(rho0).matrix[0, 0]
        assert m[0, 0].real.item() == pytest.approx(0.36 * math.exp(-gamma * t), abs=1e-13)
        assert m[0, 1].abs().item() == pytest.approx(0.48 * math.exp(-gamma * t / 2), abs=1e-13)
    rk4 = dynamics.lindblad_evolve(rho0, H, [L], [0.0, 2.0], rates=[gamma], substeps=100)[-1].matrix[0, 0]
    assert rk4[0, 0].real.item() == pytest.approx(0.36 * math.exp(-gamma * 2.0), abs=1e-9)


def test_sec4_coherence_time_of_damping_and_dephasing():
    """Sec. 4: jumps S^- (rate g1) and S^z = Z/2 (rate g2) give T_1 = 1/g1, 1/T_phi = g2/2 and 1/T_2 = (g1 + g2)/2."""
    alg = alg_of([2])                                                   # basis |1> = up, |0> = down
    Sm = from_blocks(alg, [torch.tensor([[0.0, 0.0], [1.0, 0.0]], dtype=DT)[None]])
    Sz = from_blocks(alg, [torch.diag(torch.tensor([0.5, -0.5], dtype=DT))[None]])
    psi = torch.tensor([0.6, 0.8], dtype=DT)
    rho0 = DensityMatrix(alg, matrix=torch.outer(psi, psi.conj())[None, None])
    g1, g2, t = 0.4, 0.7, 1.3
    m = dynamics.lindblad_channel(alg, None, [Sm, Sz], t=t, rates=[g1, g2])(rho0).matrix[0, 0]
    assert m[0, 0].real.item() == pytest.approx(0.36 * math.exp(-g1 * t), abs=1e-13)
    assert m[0, 1].abs().item() == pytest.approx(0.48 * math.exp(-(g1 + g2) * t / 2), abs=1e-13)


def test_sec4_xy_mixer_post_selection():
    """Sec. 4 (ex9_qaoa_xy_mixer.py): damping only gives acceptance e^{-w g1 T} and fidelity 1; with dephasing the post-selected state does not depend on g1."""
    n, w, T = 4, 2, 1.0
    q = SpinChain(n, 'periodic')
    H = q.xxz(J=1.0, Delta=0.0)
    rho0 = q.basis_state('1100')
    ideal = dynamics.von_neumann(rho0, H, [T])[0]
    P = q.algebra.central([1.0 if N == w else 0.0 for N in range(n + 1)])
    damp = [q.lowering(i) for i in range(n)]
    deph = [q.sz(i) for i in range(n)]

    def run(g1, g2):
        rho = dynamics.lindblad_evolve(rho0, H, damp + deph, [0.0, T], rates=[g1] * n + [g2] * n, substeps=60)[-1]
        post, p = rho.condition_on(P)
        return post, p.item(), states.trace_of_product(post, ideal).real.item()

    _, p, F = run(0.1, 0.0)
    assert p == pytest.approx(math.exp(-w * 0.1 * T), abs=1e-8)                    # RK4 accuracy
    assert F == pytest.approx(1.0, abs=1e-8)
    postA, pA, FA = run(0.0, 0.2)
    postB, pB, FB = run(0.15, 0.2)
    assert pA == pytest.approx(1.0, abs=1e-12) and pB == pytest.approx(math.exp(-w * 0.15 * T), abs=1e-8)
    assert (postA.matrix - postB.matrix).abs().max() < 1e-8 and FA < 1 - 1e-3


def test_sec4_page_formula():
    """Sec. 4/5 (ex5_page_entropy.py): <S(rho_A)> = sum_{k=n+1}^{mn} 1/k - (m-1)/(2n) for Haar pure states on C^m (x) C^n (z-score)."""
    m, n, B = 2, 3, 20000
    A, Bs = alg_of([m]), alg_of([n])
    AB = composite.tensor_product(A, Bs)
    psi = states.random_density_matrix(AB, batch_size=B, rank=1)
    S = states.von_neumann_entropy(composite.partial_trace(psi, keep=1))
    page = sum(1.0 / k for k in range(n + 1, m * n + 1)) - (m - 1) / (2 * n)
    assert abs(zscore(S, page)) < 5


def test_sec4_exp1_supremum_is_the_trace_norm():
    """Sec. 4/6 (ex6_inequality_search.py): sup_U z = ||YX||_1 - Tr(XY) > 0 for non-commuting X, Y >= 0, attained at U = W^* (YX = W|YX|)."""
    alg = alg_of([4])
    X, Y = rand_op(alg, 6, 'pos'), rand_op(alg, 6, 'pos')
    sup = (Y @ X).trace_norm() - (X @ Y).trace.real
    assert (sup > 1e-6).all()
    W, _ = (Y @ X).polar()
    z = (X @ W.adjoint() @ Y).trace.abs() - (X @ Y).trace.real
    assert torch.allclose(z, sup, atol=1e-10)
    U = alg.random_unitary_operator(2000)
    zs = (X.matrix[:1] @ U.matrix @ Y.matrix[:1]).diagonal(dim1=-2, dim2=-1).sum((-1, -2)).abs() \
        - (X @ Y).trace.real[0]
    assert (zs <= sup[0] + 1e-10).all()
    # commuting X, Y: no violation
    Xc = from_blocks(alg, [torch.diag(torch.rand(1, 4, dtype=torch.float64) + 0.1).to(DT)])
    Yc = from_blocks(alg, [torch.diag(torch.rand(1, 4, dtype=torch.float64) + 0.1).to(DT)])
    assert ((Yc @ Xc).trace_norm() - (Xc @ Yc).trace.real).abs().item() < 1e-12


# ================================================================================================================
# Section 5: validation against exact results
# ================================================================================================================

def test_sec5_weingarten_value():
    """Sec. 5: E tau((A U B U^*)^2) = -1/(n^2 - 1) for traceless involutions A = B and Haar U in M_n (z-score)."""
    n, B = 4, 30000
    alg = alg_of([n])
    A = torch.diag(torch.tensor([1.0, 1.0, -1.0, -1.0], dtype=DT))
    U = alg.random_unitary(n, batch_size=B)
    X = A @ U @ A @ U.mH
    v = torch.diagonal(X @ X, dim1=-2, dim2=-1).sum(-1).real / n
    assert abs(zscore(v, -1 / (n ** 2 - 1))) < 5


@pytest.mark.parametrize('k,r', [(4, 2), (3, 3), (5, 1)])
def test_sec5_purity_of_induced_states(k, r):
    """Sec. 5: E Tr rho^2 = (k + r)/(k r + 1) for random states of M_k induced by an ancilla of dimension r (z-score)."""
    alg = alg_of([k])
    rho = states.random_density_matrix(alg, batch_size=20000, rank=r)
    pur = states.purity(rho).real
    exact = (k + r) / (k * r + 1)
    if r == 1:
        assert (pur - 1).abs().max() < 1e-12
    else:
        assert abs(zscore(pur, exact)) < 5


# ================================================================================================================
# Section 6: sampling versus optimisation
# ================================================================================================================

def exp3_z(X, Y):
    return (Y @ X.abs() @ Y).trace.real - (Y @ X @ Y).trace_norm()


@pytest.mark.parametrize('ks', [[1], [3], [2, 3]], ids=str)
def test_prop1_exp3_is_nonnegative(ks):
    """Proposition 1: Tr|YXY| <= Tr(Y|X|Y) for X = X^*, Y >= 0 (random, rank-deficient Y, scalar X)."""
    alg = alg_of(ks)
    X, Y = rand_op(alg, 200, 'herm'), rand_op(alg, 200, 'pos')
    assert (exp3_z(X, Y) >= -1e-11).all()
    Yr = from_blocks(alg, [G @ G.mH for G in [torch.randn(50, k, 1, dtype=DT) for k in ks]])   # rank one per block
    Xr = rand_op(alg, 50, 'herm')
    assert (exp3_z(Xr, Yr) >= -1e-11).all()
    assert exp3_z(-2.0 * alg.identity(), Y).abs().max() < 1e-10                  # X scalar: equality


def x_parts(Xb):
    w, V = torch.linalg.eigh(Xb)
    Xp = (V * w.clamp(min=0).to(DT).unsqueeze(-2)) @ V.mH
    Xm = (V * (-w).clamp(min=0).to(DT).unsqueeze(-2)) @ V.mH
    return Xp, Xm


def test_prop1_equality_condition():
    """Proposition 1: equality iff Y X_+ Y^2 X_- Y = 0 (commuting X, Y; also non-commuting with orthogonal supports)."""
    alg = alg_of([4])
    # commuting: Y a function of X
    Xb = blocks_of(rand_op(alg, 20, 'herm'))[0]
    Yb = torch.linalg.matrix_exp(0.3 * Xb)
    assert exp3_z(from_blocks(alg, [Xb]), from_blocks(alg, [Yb])).abs().max() < 1e-10
    # non-commuting but Y leaves ran X_+ and ran X_- invariant: Y = Y_1 (+) Y_2 on C^2 (+) C^2, X_+ on the first summand
    Gp, Gm = torch.randn(20, 2, 2, dtype=DT), torch.randn(20, 2, 2, dtype=DT)
    Xb = torch.zeros(20, 4, 4, dtype=DT)
    Xb[:, :2, :2] = Gp @ Gp.mH
    Xb[:, 2:, 2:] = -(Gm @ Gm.mH)
    Y1, Y2 = torch.randn(20, 2, 2, dtype=DT), torch.randn(20, 2, 2, dtype=DT)
    Yb = torch.zeros(20, 4, 4, dtype=DT)
    Yb[:, :2, :2] = Y1 @ Y1.mH
    Yb[:, 2:, 2:] = Y2 @ Y2.mH
    assert (Xb @ Yb - Yb @ Xb).abs().amax((-1, -2)).min() > 1e-3                    # X and Y do not commute
    assert exp3_z(from_blocks(alg, [Xb]), from_blocks(alg, [Yb])).abs().max() < 1e-10
    # generic: the condition fails and the inequality is strict
    X, Y = rand_op(alg, 50, 'herm'), rand_op(alg, 50, 'pos')
    Xp, Xm = x_parts(blocks_of(X)[0])
    Yb = blocks_of(Y)[0]
    cond = (Yb @ Xp @ Yb @ Yb @ Xm @ Yb).abs().amax((-1, -2))
    z = exp3_z(X, Y)
    assert ((cond > 1e-6) <= (z > 1e-12)).all() and (z > 1e-8).all()


def test_prop1_proof_lemma_trace_norm_of_a_difference():
    """Proposition 1 (proof): for A, B >= 0, ||A - B||_1 <= Tr A + Tr B with equality iff AB = 0."""
    alg = alg_of([4])
    A, Bp = rand_op(alg, 30, 'pos'), rand_op(alg, 30, 'pos')
    assert ((A - Bp).trace_norm() < A.trace.real + Bp.trace.real - 1e-8).all()
    V = alg.random_unitary(4, batch_size=30)
    a, b = torch.zeros(30, 4, dtype=DT), torch.zeros(30, 4, dtype=DT)
    a[:, :2], b[:, 2:] = torch.rand(30, 2) + 0.1, torch.rand(30, 2) + 0.1
    Ao = alg.operator(((V * a.unsqueeze(-2)) @ V.mH)[:, None])
    Bo = alg.operator(((V * b.unsqueeze(-2)) @ V.mH)[:, None])
    assert (Ao @ Bo).matrix.abs().max() < 1e-12
    assert torch.allclose((Ao - Bo).trace_norm(), Ao.trace.real + Bo.trace.real, atol=1e-12)


def spectral_op(alg, eig_blocks, unitaries=None):
    """Operator with prescribed eigenvalues per block (list of (B, k_c)) and eigenvectors (identity or given)."""
    bl = []
    for c, e in enumerate(eig_blocks):
        D = torch.diag_embed(e.to(alg.dtype))
        if unitaries is not None:
            U = unitaries[c]
            D = U @ D @ U.mH
        bl.append(D)
    return from_blocks(alg, bl)


@pytest.mark.parametrize('k,C,cplx', [(2, 1, True), (3, 2, False), (4, 3, True), (2, 4, False)])
def test_sec6_exp1_lower_envelope_exact(k, C, cplx):
    """Sec. 6: for Delta(X) = delta and all sectors of size k >= 2, inf z = -(kC - 2 delta/(1 + delta)) (Y = 1, cyclic U)."""
    alg = alg_of([k] * C, complex_valued=cplx)
    for delta in (0.0, 0.3, 0.9):
        lo = (1 - delta) / (1 + delta)
        eig = [torch.ones(1, k, dtype=torch.float64) for _ in range(C)]
        eig[0][0, -1] = lo                                              # kC - 1 eigenvalues 1, one lambda_min
        Vs = [alg.random_unitary(k, batch_size=1) for _ in range(C)]
        X = spectral_op(alg, eig, Vs)
        assert X.michelson_contrast.item() == pytest.approx(delta, abs=1e-12)
        shift = torch.roll(torch.eye(k, dtype=alg.dtype), 1, dims=0)
        U = from_blocks(alg, [Vs[c] @ shift.to(DT if cplx else torch.float64) @ Vs[c].mH for c in range(C)])
        Y = alg.identity()
        z = (X @ U @ Y).trace.abs() - (X @ Y).trace.real
        assert z.item() == pytest.approx(-(k * C - 2 * delta / (1 + delta)), abs=1e-12)
        # every admissible configuration lies above: random Y >= 0 with lambda_max(Y) = 1 and random U
        Yr = rand_op(alg, 300, 'pos')
        Yr = alg.operator(Yr.matrix / Yr.lambda_max[:, None, None, None])
        er = [lo + (1 - lo) * torch.rand(300, k, dtype=torch.float64) for _ in range(C)]
        er[0][:, 0], er[0][:, 1] = 1.0, lo                             # Delta(X) = delta exactly
        Xr = spectral_op(alg, er, [alg.random_unitary(k, batch_size=300) for _ in range(C)])
        Ur = alg.random_unitary_operator(300)
        zr = (Xr @ Ur @ Yr).trace.abs() - (Xr @ Yr).trace.real
        assert (zr >= -(k * C - 2 * delta / (1 + delta)) - 1e-10).all()


@pytest.mark.parametrize('ks', [[2], [3], [2, 3, 4], [1, 5]], ids=str)
def test_sec6_exp2_closed_form(ks):
    """Sec. 6: Exp. 2 has sup z = -inf z = sum_c sum_i x_ci^desc (y_ci^desc - y_ci^asc), attained; <= sum_c floor(k_c/2)(1-lx)(1-ly)."""
    alg = alg_of(ks)
    x = [torch.sort(torch.rand(1, k, dtype=torch.float64) + 0.2, descending=True)[0] for k in ks]
    y = [torch.sort(torch.rand(1, k, dtype=torch.float64) + 0.2, descending=True)[0] for k in ks]
    sup = sum((xc * (yc - yc.flip(-1))).sum() for xc, yc in zip(x, y)).item()
    rev = [torch.eye(k, dtype=DT).flip(0) for k in ks]
    X0 = spectral_op(alg, x)
    U = alg.random_unitary_operator(1)                                            # drops out of z
    Xg = U @ X0                                                                   # X = U X_0, |X| = X_0

    def z(Y0, V):
        Yg = V @ Y0
        return ((Xg @ Yg).trace_norm() - (Xg.abs() @ Yg.abs()).trace.real).item()

    assert (Xg.abs() - X0).matrix.abs().max() < 1e-10
    Y_asc = spectral_op(alg, [yc.flip(-1) for yc in y])                       # anti-aligned with X_0
    Y_desc = spectral_op(alg, y)                                                # aligned with X_0
    Vrev = from_blocks(alg, [r[None] for r in rev])
    # sup: Tr(X0 Y0) = sum x^desc y^asc and ||X0 V Y0||_1 = sum x^desc y^desc
    assert z(Y_asc, Vrev) == pytest.approx(sup, abs=1e-10)
    assert z(Y_desc, Vrev) == pytest.approx(-sup, abs=1e-10)
    # random configurations with the same spectra stay in [-sup, sup]
    for _ in range(20):
        Vs = [alg.random_unitary(k, batch_size=1) for k in ks]
        Yr = spectral_op(alg, y, Vs)
        zr = z(Yr, alg.random_unitary_operator(1))
        assert -sup - 1e-10 <= zr <= sup + 1e-10
    # pairing identity and the two-level bound
    lx = min(xc.min().item() for xc in x) / max(xc.max().item() for xc in x)
    ly = min(yc.min().item() for yc in y) / max(yc.max().item() for yc in y)
    xn = [xc / max(t.max() for t in x) for xc in x]
    yn = [yc / max(t.max() for t in y) for yc in y]
    pair = sum(sum((xc[0, i] - xc[0, -1 - i]) * (yc[0, i] - yc[0, -1 - i]) for i in range(xc.shape[1] // 2))
               for xc, yc in zip(xn, yn))
    full = sum((xc * (yc - yc.flip(-1))).sum() for xc, yc in zip(xn, yn))
    assert float(pair) == pytest.approx(float(full), abs=1e-12)
    assert float(full) <= sum(k // 2 for k in ks) * (1 - lx) * (1 - ly) + 1e-12


@pytest.mark.parametrize('k,C', [(2, 1), (4, 1), (2, 3), (16, 2)])
def test_sec6_exp2_two_level_envelope(k, C):
    """Sec. 6: with Delta(|X|) = delta, Delta(|Y|) = D, two-level spectra attain C floor(k/2) 2delta/(1+delta) (1 - lo(D))."""
    alg = alg_of([k] * C, complex_valued=False)
    delta, D = 0.6, 0.999
    lx, ly = (1 - delta) / (1 + delta), (1 - D) / (1 + D)
    xs = torch.tensor([1.0] * (k // 2) + [lx] * (k - k // 2), dtype=torch.float64)[None]
    ya = torch.tensor([ly] * (k - k // 2) + [1.0] * (k // 2), dtype=torch.float64)[None]       # ascending
    X0 = spectral_op(alg, [xs] * C)
    Y0 = spectral_op(alg, [ya] * C)
    Vrev = from_blocks(alg, [torch.eye(k, dtype=alg.dtype).flip(0)[None]] * C)
    z = ((X0 @ Vrev @ Y0).trace_norm() - (X0 @ Y0).trace.real).item()
    expected = C * (k // 2) * 2 * delta / (1 + delta) * (1 - 0.001 / 1.999)
    assert z == pytest.approx(expected, rel=1e-12)
    assert 1 - ly == pytest.approx(1 - 0.001 / 1.999, rel=1e-14)


def test_sec6_direct_sums_are_admissible():
    """Sec. 6 (ii): z is additive over orthogonal blocks (Exp. 1 after aligning signs with -U) and the contrast of a sum is the largest."""
    a2 = alg_of([2], complex_valued=False)
    a4 = alg_of([4], complex_valued=False)
    a22 = alg_of([2, 2], complex_valued=False)
    torch.manual_seed(7)
    lam = lambda d: torch.tensor([[1.0, (1 - d) / (1 + d)]], dtype=torch.float64)
    Xs = [spectral_op(a2, [lam(d)], [a2.random_unitary(2, batch_size=1)]) for d in (0.5, 0.3)]
    Ys = [spectral_op(a2, [lam(d)], [a2.random_unitary(2, batch_size=1)]) for d in (0.8, 0.2)]
    Us = [a2.random_unitary_operator(1) for _ in range(2)]

    def bd(ops, alg):                                   # block-diagonal embedding into M_4 or M_2 (+) M_2
        if alg is a4:
            M = torch.zeros(1, 4, 4, dtype=torch.float64)
            M[:, :2, :2], M[:, 2:, 2:] = blocks_of(ops[0])[0], blocks_of(ops[1])[0]
            return from_blocks(alg, [M])
        return from_blocks(alg, [blocks_of(o)[0] for o in ops])

    t1 = [(x @ u @ y).trace.item() for x, u, y in zip(Xs, Us, Ys)]
    Us_al = [u if t >= 0 else -1.0 * u for u, t in zip(Us, t1)]                  # align the signs of the traces
    z1 = [abs(t) - (x @ y).trace.item() for t, x, y in zip(t1, Xs, Ys)]
    z2 = [((x @ u @ y).trace_norm() - (x.abs() @ y.abs()).trace.real).item() for x, u, y in zip(Xs, Us, Ys)]
    z3 = [exp3_z(x - 0.75 * a2.identity(), y).item() for x, y in zip(Xs, Ys)]
    for alg in (a4, a22):
        X, Y, U = bd(Xs, alg), bd(Ys, alg), bd(Us_al, alg)
        assert X.michelson_contrast.item() == pytest.approx(0.5, abs=1e-12)          # the larger contrast
        assert Y.michelson_contrast.item() == pytest.approx(0.8, abs=1e-12)
        zz1 = ((X @ U @ Y).trace.abs() - (X @ Y).trace).item()
        assert zz1 == pytest.approx(sum(z1), abs=1e-12)
        zz2 = ((X @ bd(Us, alg) @ Y).trace_norm() - (X.abs() @ Y.abs()).trace.real).item()
        assert zz2 == pytest.approx(sum(z2), abs=1e-12)
        assert exp3_z(X - 0.75 * alg.identity(), Y).item() == pytest.approx(sum(z3), abs=1e-12)


def test_sec6_parametrisation_has_the_exact_contrast():
    """Sec. 6: the spectra are lambda_max = 1, lambda_min = (1-delta)/(1+delta), others lo + (1-lo)(1+sin theta)/2 (every iterate admissible)."""
    alg = alg_of([2, 3], complex_valued=False)
    d = torch.tensor([0.0, 0.25, 0.7, 0.999], dtype=torch.float64).repeat_interleave(5)
    X = opt.PositiveParam(alg, len(d), delta=d)
    with torch.no_grad():
        X.spectrum.theta.normal_(0, 3)
        Xop = X.operator()
    assert torch.allclose(Xop.michelson_contrast, d, atol=1e-12)
    assert torch.allclose(Xop.lambda_max, torch.ones_like(d), atol=1e-12)
    assert torch.allclose(Xop.lambda_min, (1 - d) / (1 + d), atol=1e-12)
    S = opt.SelfAdjointParam(alg, len(d), delta=d)
    absS = S.abs_operator()
    assert torch.allclose(absS.michelson_contrast, d, atol=1e-12)
    assert torch.allclose(S.operator().abs().matrix, absS.matrix, atol=1e-10)
    # (1 + sin theta)/2 reaches both ends at finite theta
    with torch.no_grad():
        X.spectrum.theta.fill_(math.pi / 2)
    ev = torch.cat(X.operator().eigenvalues(), dim=-1).real
    lo = ((1 - d) / (1 + d))[:, None]
    assert ((ev - 1).abs() < 1e-12).sum(-1).min() == 4                                # all but lambda_min at 1
    assert ((ev - lo).abs() < 1e-12).any(-1).all()


@pytest.mark.parametrize('exp', [1, 3])
def test_sec6_k2_closed_forms_of_the_upper_envelopes(exp):
    """Sec. 6: for k = 2, C = 1 the upper envelopes are d^2 D^2/((1+d)(1+D)) (Exp. 1) and (1-d)(2D/(1+D))^2 (Exp. 3)."""
    torch.manual_seed(3)
    alg = alg_of([2], complex_valued=False)
    D, starts = 0.999, 16
    deltas = torch.tensor([0.3, 0.6, 0.9], dtype=torch.float64)
    d = deltas.repeat_interleave(starts)
    X = opt.SelfAdjointParam(alg, len(d), delta=d) if exp == 3 else opt.PositiveParam(alg, len(d), delta=d)
    Y = opt.PositiveParam(alg, len(d), delta_range=(0.0, D))
    if exp == 1:
        U = opt.UnitaryParam(alg, len(d))
        obj = lambda: (X.operator() @ U.operator() @ Y.operator()).trace.abs() - (X.operator() @ Y.operator()).trace.real
        params = [X, Y, U]
    else:
        obj = lambda: (Y.operator() @ X.abs_operator() @ Y.operator()).trace.real \
            - (Y.operator() @ X.operator() @ Y.operator()).trace_norm()
        params = [X, Y]
    best = opt.extremize(obj, params, maximize=True, rounds=5, steps=40).reshape(3, starts).max(-1)[0]
    closed = deltas ** 2 * D ** 2 / ((1 + deltas) * (1 + D)) if exp == 1 else (1 - deltas) * (2 * D / (1 + D)) ** 2
    assert (best <= closed + 1e-10).all()                                  # never above the conjectured envelope
    assert torch.allclose(best, closed, rtol=2e-3)


# ================================================================================================================
# Section 6.1: Renyi divergences
# ================================================================================================================

def renyi_np(r, s, a, kind):
    if kind == 'petz':
        return math.log(np.trace(mpow(r, a) @ mpow(s, 1 - a)).real) / (a - 1)
    g = (1 - a) / (2 * a)
    return math.log(np.trace(mpow(mpow(s, g) @ r @ mpow(s, g), a)).real) / (a - 1)


@pytest.mark.parametrize('ks', [[3], [2, 2]], ids=str)
def test_sec61_divergence_definitions(ks):
    """Sec. 6.1: D_alpha, ~D_alpha as defined; ~D_alpha <= D_alpha; ~D_{1/2} = -log F with F = ||sqrt(rho) sqrt(sigma)||_1^2."""
    alg = alg_of(ks)
    rho = states.random_density_matrix(alg, batch_size=6)
    sig = states.random_density_matrix(alg, batch_size=6)
    for a in (0.3, 0.5, 0.8, 1.5, 2.0, 3.0):
        Dp = dv.petz_renyi_divergence(rho, sig, a)
        Ds = dv.sandwiched_renyi_divergence(rho, sig, a)
        assert (Ds <= Dp + 1e-10).all()
        for b in range(2):
            # dense reference on the block-diagonal matrices
            R = torch.block_diag(*[X[b] for X in blocks_of(rho)]).numpy()
            S = torch.block_diag(*[X[b] for X in blocks_of(sig)]).numpy()
            assert Dp[b].item() == pytest.approx(renyi_np(R, S, a, 'petz'), abs=1e-10)
            assert Ds[b].item() == pytest.approx(renyi_np(R, S, a, 'sand'), abs=1e-10)
    F = (rho.sqrt() @ sig.sqrt()).trace_norm() ** 2
    assert torch.allclose(dv.sandwiched_renyi_divergence(rho, sig, 0.5), -torch.log(F), atol=1e-10)


def test_sec61_limits_alpha_to_one_and_infinity():
    """Sec. 6.1: both families tend to the Umegaki relative entropy as alpha -> 1; ~D_alpha increases to D_max as alpha -> inf."""
    alg = alg_of([2, 3])
    rho = states.random_density_matrix(alg, batch_size=4)
    sig = states.random_density_matrix(alg, batch_size=4)
    D = dv.umegaki_relative_entropy(rho, sig)
    ref = torch.stack([torch.tensor(np.trace(R @ (logm_h(R) - logm_h(S))).real) for R, S in dense_pairs(rho, sig)])
    assert torch.allclose(D, ref, atol=1e-10)
    for eps in (1e-5, -1e-5):
        for f in (dv.petz_renyi_divergence, dv.sandwiched_renyi_divergence):
            assert torch.allclose(f(rho, sig, 1 + eps), D, atol=1e-4)
    alphas = [2.0, 5.0, 20.0, 100.0, 1e4]
    vals = torch.stack([dv.sandwiched_renyi_divergence(rho, sig, a) for a in alphas])
    assert (vals[1:] >= vals[:-1] - 1e-10).all()
    Dmax = dv.max_relative_entropy(rho, sig)
    assert (vals <= Dmax + 1e-10).all() and torch.allclose(vals[-1], Dmax, atol=2e-3)


def logm_h(A):
    w, V = np.linalg.eigh(A)
    return (V * np.log(w)) @ V.conj().T


def dense_pairs(rho, sig):
    for b in range(rho.matrix.shape[0]):
        yield (torch.block_diag(*[X[b] for X in blocks_of(rho)]).numpy(),
               torch.block_diag(*[X[b] for X in blocks_of(sig)]).numpy())


def test_sec61_trace_independence_and_sector_split():
    """Sec. 6.1: the divergences do not depend on the trace (central weights cancel) and split into classical and conditional parts."""
    alg = alg_of([2, 3])
    rho = states.random_density_matrix(alg, batch_size=3)
    sig = states.random_density_matrix(alg, batch_size=3)
    for a in (0.6, 2.5):
        for f in (dv.petz_renyi_divergence, dv.sandwiched_renyi_divergence):
            base = f(rho, sig, a)
            for t in ('norm', 'tau_vN'):
                v = f(alg.operator(rho.density(t)), alg.operator(sig.density(t)), a, trace=t)
                assert torch.allclose(v, base, atol=1e-11)
            # exp((a-1) D) = sum_c p_c^a q_c^{1-a} exp((a-1) D(rho_c/p_c || sigma_c/q_c))
            p, q = rho.sector_probabilities(), sig.sector_probabilities()
            tot = 0
            for c, k in enumerate(alg.k_factors):
                ac = alg_of([k])
                rc = DensityMatrix(ac, matrix=(blocks_of(rho)[c] / p[:, c, None, None])[:, None])
                sc = DensityMatrix(ac, matrix=(blocks_of(sig)[c] / q[:, c, None, None])[:, None])
                tot = tot + p[:, c] ** a * q[:, c] ** (1 - a) * torch.exp((a - 1) * f(rc, sc, a))
            assert torch.allclose(torch.exp((a - 1) * base), tot, atol=1e-11)


def stinespring_channel(alg, d_env, gen):
    """A channel of M_d with a d_env-dimensional environment, from V = Z (Z^*Z)^{-1/2}."""
    d = alg.k_factors[0]
    Z = torch.randn(d * d_env, d, dtype=DT, generator=gen)
    w, Vv = torch.linalg.eigh(Z.mH @ Z)
    Vi = Z @ (Vv * w.rsqrt().to(DT).unsqueeze(-2)) @ Vv.mH
    K = Vi.reshape(d_env, d, d)                                   # Kraus operators K_e = (<e| (x) 1) V
    return channels.kraus_channel([from_blocks(alg, [K[e][None]]) for e in range(d_env)])


def test_sec61_data_processing_inside_the_ranges():
    """Sec. 6.1: DPI holds for Petz alpha in (0,1) u (1,2] and sandwiched alpha >= 1/2 (random states, Stinespring channels, dephasing)."""
    gen = torch.Generator().manual_seed(11)
    alg = alg_of([3])
    rng = np.random.default_rng(5)
    chans = [stinespring_channel(alg, 3, gen) for _ in range(3)] + [channels.dephasing_channel(alg, 1.0)]
    for ch in chans:
        assert ch.is_trace_preserving()
        for _ in range(4):
            rho = state_op(alg, rand_state_np(3, rng, 0.05))
            sig = state_op(alg, rand_state_np(3, rng, 0.05))
            for a in (0.3, 0.7, 1.5, 2.0):
                v = dv.petz_renyi_divergence(ch(rho), ch(sig), a) - dv.petz_renyi_divergence(rho, sig, a)
                assert v.item() <= 1e-10
            for a in (0.5, 0.8, 1.5, 3.0, 10.0):
                v = dv.sandwiched_renyi_divergence(ch(rho), ch(sig), a) - dv.sandwiched_renyi_divergence(rho, sig, a)
                assert v.item() <= 1e-10


def test_sec61_sandwiched_dpi_fails_below_one_half():
    """Sec. 6.1: outside alpha >= 1/2 the sandwiched DPI fails: explicit qubit states and complete dephasing at alpha = 0.3."""
    alg = alg_of([2])

    def bloch(x, z):
        return np.array([[1 + z, x], [x, 1 - z]]) / 2

    u = np.array([3246.94, 28751.63]); u = 0.98 * u / np.linalg.norm(u)
    s = np.array([174960.71, -305199.04]); s = 0.98 * s / np.linalg.norm(s)
    R, S = bloch(*u), bloch(*s)
    rho, sig = state_op(alg, R), state_op(alg, S)
    deph = channels.dephasing_channel(alg, 1.0)
    for a, violated in ((0.3, True), (0.5, False), (0.7, False)):
        v = (dv.sandwiched_renyi_divergence(deph(rho), deph(sig), a) - dv.sandwiched_renyi_divergence(rho, sig, a)).item()
        ref = renyi_np(np.diag(np.diag(R)), np.diag(np.diag(S)), a, 'sand') - renyi_np(R, S, a, 'sand')
        assert v == pytest.approx(ref, abs=1e-10)
        assert (v > 0.1) if violated else (v <= 1e-12)


def test_sec61_petz_dpi_fails_above_two():
    """Sec. 6.1: outside alpha <= 2 the Petz DPI fails: explicit two-qubit states and the partial trace at alpha = 2.5, 3."""
    eps = 1e-3
    psi = np.array([-0.38, 0.91, -0.06, 0.15])
    psi /= np.linalg.norm(psi)
    phi = np.array([0.9, 0.37, 0.21, 0.09])
    phi -= (phi @ psi) * psi
    phi /= np.linalg.norm(phi)
    R = (1 - eps) * np.outer(psi, psi) + eps * np.eye(4) / 4
    S = (1 - eps) * (0.64 * np.outer(phi, phi) + 0.36 * np.outer(psi, psi)) + eps * np.eye(4) / 4
    A = alg_of([2])
    AB = composite.tensor_product(A, A)
    rho, sig = state_op(AB, R), state_op(AB, S)
    ptr = lambda X: np.trace(X.reshape(2, 2, 2, 2), axis1=1, axis2=3)
    for a, violated in ((1.5, False), (2.0, False), (2.5, True), (3.0, True)):
        v = (dv.petz_renyi_divergence(composite.partial_trace(rho, keep=1), composite.partial_trace(sig, keep=1), a)
             - dv.petz_renyi_divergence(rho, sig, a)).item()
        ref = renyi_np(ptr(R), ptr(S), a, 'petz') - renyi_np(R, S, a, 'petz')
        assert v == pytest.approx(ref, abs=1e-9)
        assert (v > 0.5) if violated else (v < 0)


def test_sec61_petz_recovery_and_universal_recovery_bound():
    """Sec. 6.1: R_{sigma,Phi}(Phi(sigma)) = sigma, and D(rho||sigma) - D(Phi rho||Phi sigma) >= -log F(rho, R_univ(Phi(rho)))."""
    gen = torch.Generator().manual_seed(2)
    alg = alg_of([3])
    rng = np.random.default_rng(9)
    ch = stinespring_channel(alg, 2, gen)
    for _ in range(3):
        rho = state_op(alg, rand_state_np(3, rng, 0.05))
        sig = state_op(alg, rand_state_np(3, rng, 0.05))
        rec = dv.petz_recovery(sig, ch, ch(sig))
        assert (rec.matrix - sig.matrix).abs().max() < 1e-10
        gap = dv.umegaki_relative_entropy(rho, sig) - dv.umegaki_relative_entropy(ch(rho), ch(sig))
        R = dv.universal_recovery(sig, ch, ch(rho))
        F = (rho.sqrt() @ R.sqrt()).trace_norm() ** 2
        assert gap.item() >= -math.log(F.item()) - 1e-9
    # a unitary channel is reversible: the Petz map recovers every state
    U = alg.random_unitary_operator(1)
    uch = channels.unitary_channel(U)
    rho = state_op(alg, rand_state_np(3, rng, 0.05))
    assert (dv.petz_recovery(sig, uch, uch(rho)).matrix - rho.matrix).abs().max() < 1e-10


# ================================================================================================================
# Section 7: symmetric tensor products and spin chains
# ================================================================================================================

def test_sec7_free_fermion_spectrum_of_the_open_xx_chain():
    """Sec. 7: open XX chain, sector N spectrum = all sums of N distinct eps_k = J cos(pi k/(L+1))."""
    L, J = 6, 1.3
    ch = SpinChain(L, 'open', complex_valued=False)
    eps = [J * math.cos(math.pi * k / (L + 1)) for k in range(1, L + 1)]
    from itertools import combinations
    for N, w in enumerate(ch.xxz(J=J, Delta=0.0).eigenvalues()):
        ref = sorted(sum(c) for c in combinations(eps, N))
        assert np.allclose(np.sort(w[0].real.numpy()), ref, atol=1e-12)


@pytest.mark.parametrize('L,E0', [(4, -2.0), (6, -2.802775638), (8, -3.651093409), (10, -4.515446354)])
def test_sec7_heisenberg_ring_energies(L, E0):
    """Sec. 7, Table 'ring': ground-state energies of H = sum S_i.S_{i+1} (nine digits), and SU(2) nesting of the sectors."""
    ch = SpinChain(L, 'periodic', complex_valued=False)
    ws = [w[0].real.numpy() for w in ch.heisenberg().eigenvalues()]
    assert ws[L // 2].min() == pytest.approx(E0, abs=5e-10)
    assert min(w.min() for w in ws) == pytest.approx(E0, abs=5e-10)
    for N in range(L // 2):                                     # every level of sector N reappears in N + 1
        dist = np.abs(ws[N][:, None] - ws[N + 1][None, :]).min(-1)
        assert dist.max() < 1e-10


def single_particle(L):
    h = np.zeros((L, L))
    for i in range(L - 1):
        h[i, i + 1] = h[i + 1, i] = 0.5                       # J/2 (S+S- + h.c.) with J = 1
    return h


def test_sec7_domain_wall_melting_matches_free_fermions():
    """Sec. 7: <n_i(t)> = [e^{iht} C(0) e^{-iht}]_ii for the domain wall in the open XX chain (sector N = L/2)."""
    L, t = 8, 2.5
    ch = SpinChain(L, 'open')
    N = L // 2
    psi0 = ch.vector_in_sector('1' * N + '0' * N).to(DT)
    H = ch.xxz(J=1.0, Delta=0.0, sector=N)
    w, V = H.eigh()[0]
    U = (V[0] * torch.exp(-1j * w[0] * t).unsqueeze(-2)) @ V[0].mH
    psit = U @ psi0
    n_lib = (psit.abs() ** 2) @ ch.bits(N).to(torch.float64)
    h = single_particle(L)
    from scipy.linalg import expm
    C0 = np.diag([1.0] * N + [0.0] * N)
    Ct = expm(1j * h * t) @ C0 @ expm(-1j * h * t)
    assert np.allclose(n_lib.numpy(), np.diag(Ct).real, atol=1e-11)


def test_sec7_peschel_entanglement_formula():
    """Sec. 7: block entropies of the open XX ground state equal -sum[nu ln nu + (1-nu) ln(1-nu)] over the restricted correlation matrix."""
    L = 8
    N = L // 2
    ch = SpinChain(L, 'open')
    w, V = ch.xxz(J=1.0, Delta=0.0, sector=N).eigh()[0]
    psi = V[0, :, 0]
    rho = DensityMatrix(ch.sector_algebra(N), matrix=torch.outer(psi, psi.conj())[None, None])
    e, phi = np.linalg.eigh(single_particle(L))
    occ = phi[:, :N]                                           # the N lowest single-particle levels
    Cm = occ @ occ.T
    for n in range(1, L):
        nu = np.clip(np.linalg.eigvalsh(Cm[:n, :n]), 1e-300, 1 - 1e-16)
        S = -np.sum(nu * np.log(nu) + (1 - nu) * np.log1p(-nu))
        assert ch.entanglement_entropy(rho, n, sector=N).item() == pytest.approx(S, abs=1e-10)


@pytest.mark.parametrize('L,m0', [(6, 3), (8, 0), (10, 5)])
def test_sec7_momentum_sectors_and_marshall(L, m0):
    """Sec. 7: momentum blocks of sector N = L/2 reproduce its spectrum; the ring ground state has k = 0 (L/2 even) or pi (L/2 odd)."""
    ring = SpinChain(L, 'periodic')
    N = L // 2
    dims = ring.momentum_dims(N)
    assert sum(dims) == comb(L, N) and len(dims) == L
    Hk = ring.xxz_momentum(J=1.0, Delta=1.0, sector=N)
    allk = np.sort(np.concatenate([w[0].real.numpy()[:dk] for w, dk in zip(Hk.eigenvalues(), dims) if dk > 0]))
    ref = np.sort(ring.heisenberg(sector=N).eigenvalues()[0][0].real.numpy())
    assert np.allclose(allk, ref, atol=1e-10)
    e0 = {m: w[0][:dk].real.min().item() for m, w, dk in zip(Hk.algebra.charges, Hk.eigenvalues(), dims) if dk > 0}
    assert min(e0, key=e0.get) == m0


def test_sec7_particle_loss_is_binomial():
    """Sec. 7: jumps S^-_i at rate g on every site give p_n(t) = binom(N0, n) e^{-n g t} (1 - e^{-g t})^{N0 - n}."""
    ch = SpinChain(4, 'periodic')
    g, t = 0.3, 1.0
    rho = dynamics.lindblad_evolve(ch.basis_state('1101'), ch.xxz(1.0, 0.7), [ch.lowering(i) for i in range(4)],
                                   [0.0, t], rates=[g] * 4, substeps=60)[-1]
    q = math.exp(-g * t)
    ref = [comb(3, n) * q ** n * (1 - q) ** (3 - n) for n in range(4)] + [0.0]
    assert np.allclose(rho.sector_probabilities()[0].numpy(), ref, atol=1e-9)


def test_sec7_lanczos_ground_states():
    """Sec. 7, Table 'krylov': Lanczos gives the free-fermion energy (sum of L/2 lowest cos(pi k/(L+1))) and the ring energy."""
    L = 12
    E0, psi, res = krylov.ground_state(SpinChain(L, 'open', complex_valued=False).xxz_sparse(J=1.0, Delta=0.0,
                                                                                             sector=L // 2))
    eps = np.sort(np.cos(np.pi * np.arange(1, L + 1) / (L + 1)))
    assert E0.item() == pytest.approx(eps[:L // 2].sum(), abs=1e-9)
    E0r, _, _ = krylov.ground_state(SpinChain(L, 'periodic', complex_valued=False).xxz_sparse(sector=L // 2))
    assert E0r.item() == pytest.approx(-5.387390917, abs=5e-10)


def test_sec7_krylov_evolution_matches_dense():
    """Sec. 7: Krylov e^{-itH} psi agrees with the dense exponential in a sector."""
    ch = SpinChain(8, 'open')
    psi0 = ch.vector_in_sector('11110000').to(DT)
    Hs = ch.xxz_sparse(J=1.0, Delta=0.5, sector=4)
    psit = krylov.evolve(Hs, psi0[None], [0.0, 1.7])
    w, V = ch.xxz(J=1.0, Delta=0.5, sector=4).eigh()[0]
    ref = (V[0] * torch.exp(-1j * w[0] * 1.7).unsqueeze(-2)) @ V[0].mH @ psi0
    assert (psit[-1, 0] - ref).abs().max() < 1e-8


# ================================================================================================================
# Section 8: inclusions, Wedderburn decomposition, SU(2)
# ================================================================================================================

def listing_inclusion():
    N = TypeIAlgebra([2, 1], [2, 1], precision='double')
    M = TypeIAlgebra([4, 3], [4, 3], precision='double')
    return N, M, inc.Inclusion(N, M, [[1, 2], [1, 1]], unitary=M.random_unitary_operator(3))


def test_sec8_listing_inclusion_index_and_expectation():
    """Sec. 8 (listings/inclusions.py): [M:N] = ||Lambda||^2 = (7 + sqrt 45)/2; E is unital, E o iota = id, bimodule, tau-preserving, iota's adjoint."""
    N, M, I = listing_inclusion()
    assert I.index() == pytest.approx((7 + math.sqrt(45)) / 2, rel=1e-13)
    for trace in ('blunt', 'norm', 'tau_vN', 'markov', torch.tensor([0.7, 1.9], dtype=torch.float64)):
        E = I.conditional_expectation(trace)
        w = I.trace_weights(trace)
        w = w[0] if w.dim() == 2 else w
        nu = I.induced_trace(trace)
        nu = nu[0] if nu.dim() == 2 else nu
        x, y, x2 = rand_op(N, 3), rand_op(M, 3), rand_op(N, 3)
        assert (E(M.identity()).matrix - N.identity().matrix).abs().max() < 1e-12
        assert (E(I.embed(x)).matrix - x.matrix).abs().max() < 1e-12
        assert (E(I.embed(x) @ y @ I.embed(x2)).matrix - (x @ E(y) @ x2).matrix).abs().max() < 1e-11
        assert torch.allclose(ref_trace(E(y), nu), ref_trace(y, w), atol=1e-11)            # tau|_N o E = tau
        lhs = ref_trace(I.embed(x).adjoint() @ y, w)
        rhs = ref_trace(x.adjoint() @ E(y), nu)
        assert torch.allclose(lhs, rhs, atol=1e-11)                                          # adjoint of iota
        Lam = torch.tensor([[1.0, 2.0], [1.0, 1.0]], dtype=torch.float64)
        assert torch.allclose(nu, Lam.T @ w, atol=1e-12)                                     # nu = Lambda^T w
        ind = E.index().matrix[0, :, 0, 0].real
        assert torch.allclose(ind, (Lam @ Lam.T @ w) / w, atol=1e-10)                         # Watatani formula
        assert torch.allclose(E.index(method='quasi_basis').matrix, E.index().matrix, atol=1e-10)
    Em = I.conditional_expectation('markov')
    assert torch.allclose(Em.index().matrix[0, :, 0, 0].real, torch.full((2,), I.index(), dtype=torch.float64),
                          atol=1e-10)


@pytest.mark.parametrize('make,index', [
    (lambda: inc.from_tensor_product(TypeIAlgebra([2], [2], precision='double'), 3), 9.0),     # (m/k)^2, m = 6, k = 2
    (lambda: inc.scalars(TypeIAlgebra([4], [4], precision='double')), 16.0),                    # n^2
    (lambda: inc.diagonal(5), 5.0),                                                             # n
])
def test_sec8_index_of_factors_scalars_and_masa(make, index):
    """Sec. 8: [M:N] = (m/k)^2 for M_k c M_m, n^2 for C c M_n and n for the diagonal C^n c M_n."""
    assert make().index() == pytest.approx(index, rel=1e-12)


def pp_brute_force(I, trace):
    """max{l : iota(E(y)) >= l y} over generic rank-one y in each sector d: 1/<xi|iota(E(y))^+|xi>, minimised over d."""
    E = I.conditional_expectation(trace)
    M = I.M
    best = math.inf
    for d, m in enumerate(M.k_factors):
        bl = [torch.zeros(1, k, k, dtype=DT) for k in M.k_factors]
        xi = torch.randn(m, dtype=DT)
        bl[d] = torch.outer(xi, xi.conj())[None]
        y = from_blocks(M, bl)
        G = blocks_of(I.embed(E(y)))[d][0]
        val = (xi.conj() @ torch.linalg.pinv(G, hermitian=True) @ xi).real.item()
        best = min(best, 1 / val)
    return best


def test_sec8_pimsner_popa_constants():
    """Sec. 8: lambda(E) = (max_d sum_c (nu_c/w_d) min(k_c, lambda_dc))^{-1}; for C c M_n, lambda(E) = 1/n and the complete constant 1/n^2."""
    N, M, I = listing_inclusion()
    Lam = np.array([[1, 2], [1, 1]])
    kN = [2, 1]
    for trace in ('blunt', 'markov', torch.tensor([0.3, 2.0], dtype=torch.float64)):
        w = I.trace_weights(trace)
        w = (w[0] if w.dim() == 2 else w).numpy()
        nu = Lam.T @ w
        formula = 1 / max(sum(nu[c] / w[d] * min(kN[c], Lam[d, c]) for c in range(2)) for d in range(2))
        lam = I.conditional_expectation(trace).pimsner_popa_constant().item()
        assert lam == pytest.approx(formula, rel=1e-12)
        assert pp_brute_force(I, trace) == pytest.approx(formula, rel=1e-9)
        # E(y) >= lambda y on random positive y
        E = I.conditional_expectation(trace)
        y = rand_op(M, 3, 'pos')                                  # batch of the unitary u of the inclusion
        gap = I.embed(E(y)) - lam * y
        assert min(e.real.min().item() for e in gap.eigenvalues()) > -1e-10
        complete = E.pimsner_popa_constant(complete=True).item()
        assert complete == pytest.approx(1 / E.index().matrix[0, :, 0, 0].real.max().item(), rel=1e-12)
    for n in (2, 3, 5):
        S = inc.scalars(TypeIAlgebra([n], [n], precision='double'))
        E = S.conditional_expectation('blunt')
        assert E.pimsner_popa_constant().item() == pytest.approx(1 / n, rel=1e-12)
        assert E.pimsner_popa_constant(complete=True).item() == pytest.approx(1 / n ** 2, rel=1e-12)


def test_sec8_basic_construction_and_jones_projection():
    """Sec. 8: M_1 has sizes Lambda^T m and inclusion matrix Lambda^T; e y e = E(y) e, e in N', Markov property tau_1(y e) = tau_1(y)/[M:N]."""
    N, M, I = listing_inclusion()
    J = I.basic_construction()
    Lam = np.array([[1, 2], [1, 1]])
    assert list(J.M.k_factors) == list(Lam.T @ np.array([4, 3]))
    assert np.array_equal(np.asarray(J.matrix), Lam.T)
    e = I.jones_projection('markov')
    E = I.conditional_expectation('markov')
    y = rand_op(M, 3)
    ey = e @ J.embed(y) @ e
    assert (ey - J.embed(I.embed(E(y))) @ e).matrix.abs().max() < 1e-10
    x = rand_op(N, 3)
    jx = J.embed(I.embed(x))
    assert (e @ jx - jx @ e).matrix.abs().max() < 1e-10
    assert (e @ e - e).matrix.abs().max() < 1e-10
    t1 = torch.as_tensor(J.markov_trace().t, dtype=torch.float64)
    lhs = ref_trace(J.embed(y) @ e, t1)
    rhs = ref_trace(J.embed(y), t1) / I.index()
    assert torch.allclose(lhs, rhs, atol=1e-11)
    R = I.relative_commutant()
    assert sum(k * k for k in R.N.k_factors) == int((Lam ** 2).sum())                # N' n M = (+) M_{lambda_dc}


@pytest.mark.parametrize('L', [0, 1, 2, 3, 4, 5])
def test_sec8_spin_chain_tower(L):
    """Sec. 8: A_L c A_{L+1} has the path incidence matrix and index 4 cos^2(pi/(2L+4)); E for the normalised trace is the partial trace."""
    T = inc.chain_tower(L + 1, L_min=L)[0]
    Lam = np.asarray(T.matrix)
    ref = np.zeros((L + 2, L + 1))
    for N in range(L + 1):
        ref[N, N] = ref[N + 1, N] = 1
    assert np.array_equal(Lam, ref)
    assert T.index() == pytest.approx(4 * math.cos(math.pi / (2 * L + 4)) ** 2, rel=1e-12)
    if L >= 1:
        ch, big = SpinChain(L), SpinChain(L + 1)
        E = T.conditional_expectation('blunt')
        nL = big.sz(L) + 0.5 * big.algebra.identity()                                # new spin up
        y = big.sz(0) @ nL
        assert (E(y).matrix - 0.5 * ch.sz(0).matrix).abs().max() < 1e-12
        assert (E(big.sz(0) @ big.sz(L)).matrix).abs().max() < 1e-12


@pytest.mark.parametrize('L', [3, 4, 5])
def test_sec8_schur_weyl_duality(L):
    """Sec. 8 (listings/symmetry.py): permutations of L qubits and total spin generate mutual commutants: k_c and m_c exchanged."""
    S = dc.from_group(dc.permutation_representation(L, 2))
    G = dc.from_generators(dc.total_spin(L))
    pairs_S = sorted(zip(S.k_factors, S.multiplicities))
    pairs_G = sorted(zip(G.multiplicities, G.k_factors))
    assert pairs_S == pairs_G
    ref = sorted((comb(L, j) - (comb(L, j - 1) if j else 0), L - 2 * j + 1) for j in range(L // 2 + 1))
    assert pairs_S == ref                                    # dims of S_L irreps (L-j, j) and 2S+1, S = L/2 - j
    assert sum(k * m for k, m in pairs_S) == 2 ** L
    assert S.commutant_dim == sum(m * m for _, m in pairs_S)
    if L == 4:
        assert (tuple(S.k_factors), tuple(S.multiplicities)) == ((1, 2, 3), (5, 1, 3))
        assert (tuple(G.k_factors), tuple(G.multiplicities)) == ((1, 3, 5), (2, 3, 1))


def test_sec8_real_complex_quaternionic_types():
    """Sec. 8: real_structure: Z_4 regular representation has two real and one complex-pair summand; the quaternions are of type H."""
    n = 4
    shift = torch.roll(torch.eye(n, dtype=DT), 1, dims=0)
    D = dc.from_group(torch.stack([torch.linalg.matrix_power(shift, j) for j in range(n)]))
    types = D.real_structure().types
    assert sorted(types) == ['C', 'C', 'R', 'R']
    # left multiplication by i and j on H = R^4: the real algebra H, complexified M_2(C) with multiplicity 2
    Li = torch.tensor([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]], dtype=DT)
    Lj = torch.tensor([[0, 0, -1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, -1, 0, 0]], dtype=DT)
    Q = dc.from_generators(torch.stack([Li, Lj]))
    assert tuple(Q.k_factors) == (2,) and tuple(Q.multiplicities) == (2,)
    assert list(Q.real_structure().types) == ['H']


def test_sec8_su2_multiplicities_and_cost():
    """Sec. 8: m_S = binom(L, L/2-S) - binom(L, L/2-S-1), sum (2S+1) m_S = 2^L; L = 14 gives 1001 and a cost factor of about 48."""
    for L in range(1, 15):
        S_vals = [Fraction(L - 2 * j, 2) for j in range(L // 2 + 1)]
        ms = {S: su2.multiplicity(L, S) for S in S_vals}
        for S, m in ms.items():
            a = int(Fraction(L, 2) - S)
            assert m == comb(L, a) - (comb(L, a - 1) if a >= 1 else 0)
        assert sum((2 * S + 1) * m for S, m in ms.items()) == 2 ** L
    m14 = [su2.multiplicity(14, S) for S in range(8)]
    assert m14 == [429, 1001, 1001, 637, 273, 77, 13, 1]
    ratio = sum(comb(14, N) ** 3 for N in range(15)) / sum(m ** 3 for m in m14)
    assert round(ratio) == 48
    m16 = [su2.multiplicity(16, S) for S in range(9)]
    assert max(m16) == 3640
    assert sum(m * comb(16, 8 + S) for S, m in enumerate(m16)) == pytest.approx(1.0e8, rel=0.01)


def test_sec8_su2_ring_lieb_mattis_and_majumdar_ghosh():
    """Sec. 8: SU2Chain reproduces the ring energy (S = 0), Lieb-Mattis ordering for open chains and E = -3L/8 (twice) at J2 = J/2."""
    ring = SU2Chain(8, 'periodic', complex_valued=False)
    gs = ring.ground_state(ring.heisenberg(J=1.0))
    assert gs.energy.item() == pytest.approx(-3.651093409, abs=5e-10) and gs.spin.item() == 0
    for L, S0 in ((5, 0.5), (6, 0.0)):
        ch = SU2Chain(L, 'open', complex_valued=False)
        H = ch.heisenberg(J=1.0)
        assert ch.ground_state(H).spin.item() == S0
        low = ch.lowest_energies(H)[0]
        assert (low[1:] > low[:-1] + 1e-10).all()
    mg = SU2Chain(8, 'periodic', complex_valued=False)
    H = mg.heisenberg(J=1.0, J2=0.5)
    w = H.eigenvalues()[mg.sector_index(0)][0].real
    assert w[0].item() == pytest.approx(-3.0, abs=1e-10) and w[1].item() == pytest.approx(-3.0, abs=1e-10)
    assert min(e[0].real.min().item() for e in H.eigenvalues()) == pytest.approx(-3.0, abs=1e-10)


def test_sec8_clebsch_gordan_condon_shortley():
    """Sec. 8: exact Clebsch-Gordan coefficients with Condon-Shortley phases (textbook values)."""
    h = Fraction(1, 2)
    cases = [((h, h, h, -h, 1, 0), 1 / math.sqrt(2)), ((h, h, h, -h, 0, 0), 1 / math.sqrt(2)),
             ((h, -h, h, h, 0, 0), -1 / math.sqrt(2)), ((1, 1, h, -h, h, h), math.sqrt(2 / 3)),
             ((1, 0, h, h, h, h), -math.sqrt(1 / 3)), ((1, 1, 1, -1, 0, 0), 1 / math.sqrt(3)),
             ((1, 0, 1, 0, 0, 0), -1 / math.sqrt(3)), ((1, 0, 1, 0, 1, 0), 0.0)]
    for args, val in cases:
        assert su2.clebsch_gordan(*args) == pytest.approx(val, abs=1e-14)
