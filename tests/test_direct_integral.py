"""Direct integrals over a measure space, discretised (torch_vn_algebra.direct_integral)."""
import math

import pytest
import torch
from scipy.special import ellipe

from torch_vn_algebra import DirectIntegral, projections


def ssh(t):
    """Bloch Hamiltonian of the SSH chain, H(k) = [[0, 1 + t e^{-ik}], [1 + t e^{ik}, 0]], for nodes (C, 1)."""
    def f(k):
        z = 1 + t * torch.exp(1j * k[:, 0])
        H = torch.zeros(k.shape[0], 2, 2, dtype=torch.complex128)
        H[:, 0, 1], H[:, 1, 0] = z.conj(), z
        return H
    return f


def ssh_energy(t):
    """Ground-state energy per cell, -int |1 + t e^{ik}| dk/2pi = -(2/pi)(1 + t) E(4t/(1 + t)^2)."""
    return -2 / math.pi * (1 + t) * ellipe(4 * t / (1 + t) ** 2)


def test_brillouin_zone_integral_converges_spectrally():
    t = 0.5
    errs = []
    for n in (8, 16, 32, 64):
        X = DirectIntegral.torus(n, dims=2)
        H = X.field(ssh(t), is_self_adjoint=True)
        e0 = X.trace(H.apply_function(lambda w: w.clamp(max=0)))
        errs.append(abs(float(e0.real) - ssh_energy(t)))
    assert errs[-1] < 1e-13 and errs[1] < 1e-2 * errs[0]           # periodic trapezoid: exponential convergence
    assert float(X.measure()) == pytest.approx(1.0)


def test_monte_carlo_and_gauss_rules():
    t = 0.5
    g = torch.Generator().manual_seed(0)
    X = DirectIntegral.monte_carlo(lambda n, gen: 2 * math.pi * torch.rand(n, generator=gen, dtype=torch.float64),
                                   20000, dims=2, generator=g)
    H = X.field(ssh(t), is_self_adjoint=True)
    e0 = float(X.trace(H.apply_function(lambda w: w.clamp(max=0))).real)
    assert abs(e0 - ssh_energy(t)) < 5 * 0.35 / math.sqrt(20000)    # |lambda_min| has standard deviation 0.345
    # Gauss-Legendre with a density: mu = 6 x (1 - x) dx on [0, 1] (mass 1), int x^2 dmu = 3/10
    Y = DirectIntegral.grid([(0, 1)], 6, dims=1, density=lambda x: 6 * x[:, 0] * (1 - x[:, 0]))
    assert float(Y.measure()) == pytest.approx(1.0, abs=1e-14)
    assert float(Y.integrate(Y.nodes[:, 0] ** 2)) == pytest.approx(0.3, abs=1e-14)


def test_varying_fibres_states_and_central_projections():
    # n(x) = 2 on [0, 1/2), 3 on [1/2, 1): int n dx = 5/2
    X = DirectIntegral.grid([(0, 1)], 10, dims=lambda x: 2 + (x[:, 0] >= 0.5).long(), rule='midpoint')
    assert X.dims == [2] * 5 + [3] * 5
    one = X.field(lambda x: torch.eye(2 if float(x[0]) < 0.5 else 3, dtype=torch.complex128))
    assert float(X.trace(one).real) == pytest.approx(2.5)
    assert float(X.normalized_trace(one).real) == pytest.approx(1.0)
    left = X.indicator(lambda x: x[:, 0] < 0.5)
    tau = X.state()                                                   # tracial state int tr . dx
    assert float(tau.expectation(left).real) == pytest.approx(0.5)
    rho = X.state(lambda x: torch.diag(torch.arange(1, 3 if float(x[0]) < 0.5 else 4,
                                                    dtype=torch.float64).to(torch.complex128)))
    assert float(rho.Tr_blunt().real) == pytest.approx(1.0)
    # phi(1_left) = int_left Tr rho / int Tr rho = (3/2) / (3/2 + 6/2)
    assert float(rho.expectation(left).real) == pytest.approx(1.5 / 4.5)
    s = X.scalar(lambda x: x[:, 0])
    assert float(X.normalized_trace(s).real) == pytest.approx(0.5)


def test_only_the_fibre_contrast_matters():
    """Theorem thm:integral of the trace paper for item 23: for a density field constant on cells, the violation on the
    direct integral is at most the largest fibre envelope delta_c/(1 - delta_c), however the cell measures differ."""
    torch.manual_seed(0)
    deltas = [0.2, 0.5, 0.3]
    X = DirectIntegral(torch.tensor([0.0, 1.0, 2.0]), torch.tensor([0.7, 0.01, 5.0]), dims=2)
    W = torch.stack([torch.diag(torch.tensor([1 + d, 1 - d], dtype=torch.complex128)) for d in deltas])
    Wop = X.field(lambda x: W, is_self_adjoint=True)

    def phi(A):
        return X.integrate(torch.einsum('cij,bcji->bc', W, A.matrix).real)
    P = projections.random_projection(X.algebra, [1, 1, 1], 4096)
    Q = projections.random_projection(X.algebra, [1, 1, 1], 4096)
    v = (phi(projections.join(P, Q)) - phi(P) - phi(Q)) / (phi(P) + phi(Q))
    bound = max(d / (1 - d) for d in deltas)
    assert float(v.max()) <= bound + 1e-12
    # approached by two close rank-one projections at the bottom eigenvector of the worst cell, zero elsewhere
    eps = 1e-4
    m = torch.zeros(2, 1, 3, 2, 2, dtype=torch.complex128)
    for j, a in enumerate((0.0, eps)):
        u = torch.tensor([math.sin(a), math.cos(a)], dtype=torch.complex128)
        m[j, 0, 1] = u[:, None] * u.conj()[None, :]
    P, Q = (X.algebra.operator(m[j], is_self_adjoint=True, is_projection=True) for j in range(2))
    v = (phi(projections.join(P, Q, 1e-12)) - phi(P) - phi(Q)) / (phi(P) + phi(Q))
    assert bound - 1e-6 < float(v) <= bound + 1e-12
    assert Wop.is_self_adjoint


def test_integer_values_and_predicates():
    X = DirectIntegral.grid([(0, 1)], 4, dims=1, rule='midpoint')
    assert float(X.integrate(torch.tensor([1, 1, 1, 1]))) == pytest.approx(1.0)        # no truncation of weights
    assert float(X.measure(lambda x: x[:, 0] > 0.5)) == pytest.approx(0.5)
    with pytest.raises(ValueError):
        X.measure(lambda x: (x[:, 0] > 0.5).long())
    with pytest.raises(ValueError):
        DirectIntegral.grid([(0, 1)], 4, dims=1, density=lambda x: x[:, 0] - 0.5)
    with pytest.raises(TypeError):
        DirectIntegral.torus(8, dims=1, rule='gauss')
