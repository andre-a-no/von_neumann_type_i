"""The projection lattice, measures on it (Gleason, frame functions on 2x2 blocks, Kochen-Specker) and spectral
measures (torch_vn_algebra.projections, measurement.PVM)."""
import math

import pytest
import torch

from torch_vn_algebra import POVM, PVM, ProjectionMeasure, TypeIAlgebra, projections, spectral_measure, states


def alg_():
    return TypeIAlgebra([1, 3, 2], [1, 3, 2], precision='double')


def span_projection(alg, sector, vectors):
    """Projection onto the span of the given vectors (columns) in one sector, zero elsewhere."""
    V = torch.linalg.qr(torch.tensor(vectors, dtype=torch.complex128).T)[0]
    m = torch.zeros(1, alg.C, alg.k_max, alg.k_max, dtype=torch.complex128)
    k = alg.k_factors[sector]
    m[0, sector, :k, :k] = V @ V.conj().T
    return alg.operator(m, is_self_adjoint=True, is_projection=True)


def test_join_and_meet():
    alg = alg_()
    P = span_projection(alg, 1, [[1, 0, 0], [0, 1, 0]])
    Q = span_projection(alg, 1, [[0, 1, 0], [0, 0, 1]])
    assert torch.allclose(projections.meet(P, Q).matrix, span_projection(alg, 1, [[0, 1, 0]]).matrix, atol=1e-12)
    full = span_projection(alg, 1, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert torch.allclose(projections.join(P, Q).matrix, full.matrix, atol=1e-12)
    # generic position: ranks add as for the dimension function (modular law of the rank)
    torch.manual_seed(0)
    P = projections.random_projection(alg, [1, 2, 1], batch_size=16)
    Q = projections.random_projection(alg, [0, 2, 1], batch_size=16)
    r = projections.ranks
    lhs = r(projections.join(P, Q)) + r(projections.meet(P, Q))
    assert torch.equal(lhs, r(P) + r(Q))
    assert torch.equal(r(projections.join(P, Q)), torch.tensor([[1, 3, 2]]).expand(16, 3))
    assert torch.allclose(projections.complement(P).matrix + P.matrix, alg.identity().matrix.expand(16, -1, -1, -1))


def test_linear_measure_and_modular_defect():
    alg = alg_()
    torch.manual_seed(1)
    rho = states.random_density_matrix(alg, batch_size=1)
    mu = ProjectionMeasure.from_state(rho)
    P = projections.random_projection(alg, [1, 1, 1], batch_size=8)
    Q = projections.random_projection(alg, [0, 2, 1], batch_size=8)
    expected = (rho.matrix * P.matrix.transpose(-2, -1)).sum(dim=(-3, -2, -1)).real
    assert torch.allclose(mu(P), expected, atol=1e-12)
    assert mu.is_linear() and not mu.is_dimension_function()
    assert mu.modular_defect(P, Q).abs().max() > 1e-3                  # a non-tracial state breaks the modular law
    dim = ProjectionMeasure.dimension(alg, [1.0, 2.0, 0.5])
    assert dim.is_dimension_function()
    assert dim.modular_defect(P, Q).abs().max() < 1e-10
    w = torch.tensor([1.0, 2.0, 0.5], dtype=torch.float64)
    assert torch.allclose(dim(P), (projections.ranks(P).to(torch.float64) * w).sum(-1))


def frame(alg, g, w2=1.0):
    """Measure with weights (1, 1, w2) per sector, linear on sectors 0 and 1, frame function g on the 2x2 sector 2."""
    rho = alg.central([1.0, 1.0 / 3, w2 / 2])          # Tr rho_c = 1, 1, w2
    return ProjectionMeasure(rho, frames={2: g})


def test_frame_function_is_additive_but_not_linear():
    alg = alg_()
    mu = frame(alg, lambda n: 0.3 * n[..., 2] ** 3)
    torch.manual_seed(2)
    P = projections.random_projection(alg, [0, 1, 1], batch_size=32)
    Pc = projections.complement(P)
    assert torch.allclose(mu(P) + mu(Pc), mu.total().expand(32), atol=1e-12)      # finitely additive
    assert not mu.is_linear() and not mu.is_dimension_function()
    with pytest.raises(ValueError):
        mu.gleason_state()
    lin = frame(alg, lambda n: 0.2 * n[..., 0] - 0.1 * n[..., 2])                  # a linear frame function
    rho = lin.gleason_state()
    assert torch.allclose(lin(P), (rho.matrix * P.matrix.transpose(-2, -1)).sum(dim=(-3, -2, -1)).real, atol=1e-12)
    with pytest.raises(ValueError):
        frame(alg, lambda n: 0.3 * n[..., 2] ** 2)                                 # not odd
    with pytest.raises(ValueError):
        ProjectionMeasure(alg.identity(), frames={1: lambda n: n[..., 0]})          # sector 1 has size 3


def test_dispersion_free_measures_exist_only_on_small_blocks():
    alg = alg_()
    torch.manual_seed(3)
    P = projections.random_projection(alg, [1, 1, 1], batch_size=64)
    for sector in (0, 2):
        mu = ProjectionMeasure.dispersion_free(alg, sector)
        v = mu(P)
        assert bool((((v.abs() < 1e-12) | ((v - 1).abs() < 1e-12))).all())
        assert torch.allclose(v + mu(projections.complement(P)), torch.ones(64, dtype=torch.float64))
    with pytest.raises(ValueError):
        ProjectionMeasure.dispersion_free(alg, 1)                                   # k = 3: Kochen-Specker


def test_lattice_items_for_measures_on_m2_are_attained():
    """Items 19 and 23 of the trace paper on M_2: a measure with values in [1 - d, 1 + d] on rank-one projections
    (mu(1) = 2) reaches delta/(1 - delta) and delta exactly when it takes its extreme value at two distinct points,
    which no state does (Proposition prop:measures of the trace paper)."""
    alg = TypeIAlgebra([2], [2], precision='double')
    d = 0.3

    def g(n):                                           # -d on the cap n_z > 0.5, +d on the opposite cap, 0 between
        return d * (-(n[..., 2] > 0.5).to(n.dtype) + (n[..., 2] < -0.5).to(n.dtype))
    mu = ProjectionMeasure(alg.central([1.0]), frames={0: g})
    assert float(mu.total()) == pytest.approx(2.0)

    def proj(theta, phi):
        u = torch.tensor([math.cos(theta / 2), complex(math.cos(phi), math.sin(phi)) * math.sin(theta / 2)],
                         dtype=torch.complex128)
        return alg.operator((u[:, None] * u.conj()[None, :]).view(1, 1, 2, 2), is_self_adjoint=True, is_projection=True)
    P, Q = proj(0.3, 0.0), proj(0.5, 1.0)              # two distinct points of the cap
    PvQ = projections.join(P, Q)
    assert projections.ranks(projections.meet(P, Q)).tolist() == [[0]]
    item23 = (mu(PvQ) - mu(P) - mu(Q)) / (mu(P) + mu(Q))
    assert float(item23) == pytest.approx(d / (1 - d), abs=1e-12)
    Pm, Qm = proj(math.pi - 0.3, 0.0), proj(math.pi - 0.5, 1.0)
    item19 = (mu(Pm) + mu(Qm) - mu(alg.identity()) - mu(projections.meet(Pm, Qm))) / \
        (mu(alg.identity()) + mu(projections.meet(Pm, Qm)))
    assert float(item19) == pytest.approx(d, abs=1e-12)


def test_spectral_measure():
    alg = alg_()
    torch.manual_seed(4)
    U = alg.random_unitary_operator(1)
    # eigenvalue 1 in all three sectors, 2 in sectors 1 and 2, -1 in sector 1
    lam = [[1.0], [1.0, 2.0, -1.0], [1.0, 2.0]]
    m = torch.zeros(1, alg.C, alg.k_max, alg.k_max, dtype=torch.complex128)
    for c, l in enumerate(lam):
        k = len(l)
        u = U.matrix[0, c, :k, :k]
        m[0, c, :k, :k] = u @ torch.diag(torch.tensor(l, dtype=torch.complex128)) @ u.conj().T
    A = alg.operator(m, is_self_adjoint=True)
    E = spectral_measure(A)
    assert torch.allclose(E.values, torch.tensor([[-1.0, 1.0, 2.0]], dtype=torch.float64), atol=1e-10)
    assert E.is_projective()
    assert torch.allclose(E.observable().matrix, A.matrix, atol=1e-10)
    assert projections.ranks(E[1]).tolist() == [[1, 1, 1]]
    pos = E.projection(lambda v: v > 0)
    assert projections.ranks(pos).tolist() == [[1, 2, 2]]
    sq = E.pushforward(lambda v: v ** 2)                       # A^2 has eigenvalues 1 (from +-1) and 4
    assert torch.allclose(sq.values, torch.tensor([[1.0, 4.0]], dtype=torch.float64))
    assert torch.allclose(sq.observable().matrix, (A.matrix @ A.matrix), atol=1e-10)
    rho = states.random_density_matrix(alg, batch_size=3)
    assert torch.allclose(E.expectation(rho), states.expectation(rho, A).real, atol=1e-10)
    # one outcome per eigenvector, for a batch
    B = alg.operator(torch.stack([m[0], 2 * m[0]]), is_self_adjoint=True)
    F = spectral_measure(B, group=False)
    assert len(F) == 6 and F.values.shape == (2, 6)
    assert torch.allclose(F.observable().matrix, B.matrix, atol=1e-10)
    with pytest.raises(ValueError):
        spectral_measure(B)


def test_pvm_from_basis_and_checks():
    alg = alg_()
    vals = torch.arange(6, dtype=torch.float64)
    E = PVM.from_povm(POVM.basis(alg), vals)
    rho = states.random_density_matrix(alg, batch_size=2)
    assert torch.allclose(E.expectation(rho), (E.probabilities(rho) * vals).sum(-1))
    assert bool((E.variance(rho) >= 0).all())
    effects = POVM.basis(alg).effects
    with pytest.raises(ValueError):
        PVM(0.5 * torch.cat([effects, effects], dim=1), algebra=alg)              # a POVM, not projective


@pytest.mark.parametrize('n', [2, 3, 4])
def test_modular_law_halves(n):
    """Proposition prop:modular of the trace paper: for phi_W with contrast d in M_n, the half M_<= of the modular law
    is violated by at most d/(1 - d) and the half M_>= by at most d, in every dimension, and both are approached."""
    alg = TypeIAlgebra([n], [n], precision='double')
    gen = torch.manual_seed(5 + n)
    for d in (0.3, 0.7):
        ell = (1 - d) / (1 + d)
        ev = torch.cat([torch.tensor([1.0, ell], dtype=torch.float64),
                        ell + (1 - ell) * torch.rand(n - 2, dtype=torch.float64, generator=gen)])
        U = alg.random_unitary_operator(1).matrix[0, 0]
        W = U @ torch.diag(ev.to(torch.complex128)) @ U.conj().T
        mu = ProjectionMeasure.from_state(alg.operator(W.view(1, 1, n, n), is_self_adjoint=True))
        for r1, r2 in [(1, 1), (1, n - 1), (n - 1, n - 1)]:
            P = projections.random_projection(alg, [r1], 256)
            Q = projections.random_projection(alg, [r2], 256)
            J, Mt = projections.join(P, Q), projections.meet(P, Q)
            le = (mu(J) + mu(Mt) - mu(P) - mu(Q)) / (mu(P) + mu(Q))
            ge = (mu(P) + mu(Q) - mu(J) - mu(Mt)) / (mu(J) + mu(Mt))
            assert float(le.max()) <= d / (1 - d) + 1e-10
            assert float(ge.max()) <= d + 1e-10
        V = torch.linalg.eigh(W)[1]
        top, bot = V[:, -1], V[:, 0]

        def rank_one(x):
            return alg.operator((x[:, None] * x.conj()[None, :]).view(1, 1, n, n), is_self_adjoint=True,
                                is_projection=True)
        eps = 1e-3
        for e, f, half, bound in [(bot, top, 'le', d / (1 - d)), (top, bot, 'ge', d)]:
            P, Q = rank_one(e), rank_one(math.cos(eps) * e + math.sin(eps) * f)
            J, Mt = projections.join(P, Q, 1e-12), projections.meet(P, Q, 1e-12)
            if half == 'le':
                v = (mu(J) + mu(Mt) - mu(P) - mu(Q)) / (mu(P) + mu(Q))
            else:
                v = (mu(P) + mu(Q) - mu(J) - mu(Mt)) / (mu(J) + mu(Mt))
            assert bound - 1e-4 < float(v) <= bound + 1e-12


def test_spectral_measure_scale_and_support():
    alg = TypeIAlgebra([3], [3], precision='double')
    A = alg.operator(1e-6 * torch.diag(torch.tensor([1.0, 2.0, 3.0], dtype=torch.complex128)).view(1, 1, 3, 3),
                     is_self_adjoint=True)
    E = spectral_measure(A)
    assert len(E) == 3                                   # small eigenvalues are not merged by an absolute tolerance
    assert torch.allclose(E.observable().matrix, A.matrix, atol=1e-20)
    run = alg.operator(torch.diag(torch.tensor([0.0, 0.7e-5, 1.4e-5], dtype=torch.complex128)).view(1, 1, 3, 3),
                       is_self_adjoint=True)
    assert len(spectral_measure(run, tol=1e-5)) == 2     # no chaining of gaps below the tolerance
    small = alg.operator(1e-11 * alg.identity().matrix, is_self_adjoint=True)
    assert projections.ranks(projections.support(small)).tolist() == [[3]]       # relative tolerance
    with pytest.raises(ValueError):
        projections.support(alg.operator(-alg.identity().matrix, is_self_adjoint=True))
    other = TypeIAlgebra([1, 2], [1, 2], precision='double')
    with pytest.raises(ValueError):
        projections.join(projections.random_projection(alg, [1]), projections.random_projection(other, [0, 1]))
