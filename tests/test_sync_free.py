"""
The computational paths read no values back to the host. On a GPU every .item(), bool(tensor), .tolist(), ... is a
synchronisation that stalls the device; per-sector synchronisations were the cause of the slow random generation
in the first full run. Here they are counted on the CPU by intercepting the Tensor methods that transfer values to
Python, for operators whose properties are known (validation of unknown properties may read values; that is
documented and tested separately).
"""
import contextlib
import warnings

import pytest
import torch

from torch_vn_algebra import TypeIAlgebra, states, channels, dynamics, optimize as opt, cost

_NAMES = ('item', 'tolist', '__bool__', '__float__', '__int__', '__index__', 'cpu', 'numpy')


@contextlib.contextmanager
def count_host_reads():
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
        yield calls
    finally:
        for n, f in saved.items():
            setattr(torch.Tensor, n, f)


def _setup(ks=(6, 6, 6, 4, 6), B=8):
    torch.manual_seed(0)
    warnings.simplefilter('ignore', cost.CostWarning)
    alg = TypeIAlgebra(list(ks), list(ks), precision='double')
    X = alg.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(B, d, dtype=torch.float64), batch_size=B,
                                      force_positive=True)
    M = X.matrix
    G = alg.operator(torch.randn(B, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype))
    rho = states.random_density_matrix(alg, B)
    Phi = channels.random_channel(alg, 3, B)
    return alg, B, M, G, rho, Phi


def _cases():
    alg, B, M, G, rho, Phi = _setup()
    H = lambda: alg.operator(M, is_self_adjoint=True)
    P = lambda: alg.operator(M, is_self_adjoint=True, is_positive=True)
    Xp, Yp, Up = opt.PositiveParam(alg, B), opt.PositiveParam(alg, B), opt.UnitaryParam(alg, B)

    def step():
        x, y = Xp.operator(), Yp.operator()
        ((x @ Up.operator() @ y).trace.abs() - (x @ y).trace.real).sum().backward()

    H0 = H()
    psi0 = torch.randn(B, alg.C, alg.k_max, dtype=alg.dtype)
    lind = dynamics.lindblad_rhs(H0, [alg.operator(0.2 * M)])
    return {
        'random unitary operator': lambda: alg.random_unitary_operator(B).matrix,
        'random state': lambda: states.random_density_matrix(alg, B).matrix,
        'random channel': lambda: channels.random_channel(alg, 3, B).kraus,
        'product, sum, trace': lambda: ((H() @ G + G) * 0.5).trace,
        'traces': lambda: (G.Tr_blunt(), G.Tr_norm(), G.tau_vN()),
        'f(A) of self-adjoint A': lambda: H().apply_function(torch.exp).matrix,
        'sqrt, log, power of positive A': lambda: (P().sqrt().matrix, P().log().matrix, P().power(0.3).matrix),
        'lambda_max, contrast': lambda: (P().lambda_max, P().michelson_contrast),
        '|A|, polar, inverse (general A)': lambda: (G.abs().matrix, G.polar()[0].matrix, G.inverse().matrix),
        'trace norm, operator norm': lambda: (G.trace_norm(), G.operator_norm(), H().trace_norm()),
        'expm': lambda: G.expm(0.1).matrix,
        'entropy, purity': lambda: (rho.von_neumann_entropy(), states.purity(rho)),
        'channel applied to a state': lambda: Phi(rho).matrix,
        'Choi matrices': lambda: Phi.choi(),
        'optimisation step': step,
        'schrodinger (exact)': lambda: dynamics.schrodinger(psi0, H0, [0.0, 0.5, 1.0]),
        'Lindblad right-hand side': lambda: lind(0.0, rho.matrix.to(alg.dtype)),
        'gibbs state': lambda: states.gibbs_state(H(), 0.5).matrix,
    }


CASES = _cases()


@pytest.mark.parametrize('name', list(CASES))
def test_no_host_reads(name):
    fn = CASES[name]
    fn()                                         # first call: one-off set-up (constants, trace-preservation check)
    with count_host_reads() as calls:
        fn()
    assert not calls, f"{name}: {len(calls)} host reads ({sorted(set(calls))})"


def test_validation_reads_once():
    """A property passed by the user is checked once, with one host read for all checks."""
    alg, B, M, *_ = _setup()
    with count_host_reads() as calls:
        A = TypeIAlgebra.Operator(alg, generator=lambda: M, is_self_adjoint=True, is_positive=True)
        A.matrix
    assert len(calls) == 1
