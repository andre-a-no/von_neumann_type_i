"""Residual channel networks stay in det > 0: the Pauli channel (XrX + YrY + ZrZ)/3, det = -1/27, is out of their reach
(scripts/residual_markov.py, Section sec:learn-markov of the paper), while one Stinespring layer fits it."""
import torch
from torch import nn

from torch_vn_algebra import DensityMatrix, TypeIAlgebra, learn

DT = torch.complex128
PAULIS = torch.stack([torch.eye(2, dtype=DT), torch.tensor([[0, 1], [1, 0]], dtype=DT),
                      torch.tensor([[0, -1j], [1j, 0]], dtype=DT), torch.tensor([[1, 0], [0, -1]], dtype=DT)])
TARGET = torch.diag(torch.tensor([1.0, -1 / 3, -1 / 3, -1 / 3], dtype=torch.float64))


def _probes(alg):
    kets = torch.tensor([[1, 0], [0, 1], [1, 1], [1, 1j]], dtype=DT)
    kets = kets / kets.norm(dim=1, keepdim=True)
    return DensityMatrix(alg, matrix=torch.einsum('bi,bj->bij', kets, kets.conj())[:, None], validate=False)


def _ptm(Y):
    """Pauli transfer matrix R_ab = Tr(s_a Phi(s_b))/2 from the images of |0>, |1>, |+>, |+i>."""
    Y = Y[:, 0]
    phi_i = Y[0] + Y[1]
    imgs = torch.stack([phi_i, 2 * Y[2] - phi_i, 2 * Y[3] - phi_i, Y[0] - Y[1]])
    return 0.5 * torch.einsum('aij,bji->ab', PAULIS, imgs).real


def _fit(net, X, steps, lr=0.05):
    opt = torch.optim.Adam(net.parameters(), lr=lr)
    for _ in range(steps):
        loss = ((_ptm(net(X).matrix) - TARGET) ** 2).sum()
        opt.zero_grad()
        loss.backward()
        opt.step()
    with torch.no_grad():
        return _ptm(net(X).matrix)


def test_residual_network_cannot_reach_the_negative_determinant_pauli_channel():
    torch.manual_seed(0)
    alg = TypeIAlgebra([2], [2], precision='double')
    X = _probes(alg)
    p = 0.2
    net = nn.Sequential(*[learn.Residual(learn.ChannelLayer(alg, kraus_rank=4), p=p) for _ in range(8)])
    with torch.no_grad():
        start = torch.linalg.matrix_norm(_ptm(net(X).matrix) - TARGET).item()
    R = _fit(net, X, 150)
    err = torch.linalg.matrix_norm(R - TARGET).item()
    assert err < start                                                     # the network did train
    # every block: eigenvalues 1 - p + p lambda, |lambda| <= 1, so |eigenvalue| >= 1 - 2p and det >= (1 - 2p)^4
    with torch.no_grad():
        for block in net:
            assert torch.linalg.det(_ptm(block(X).matrix)) >= (1 - 2 * p) ** 4 - 1e-12
    assert torch.linalg.det(R) > 0
    # nearest map with det >= 0 to diag(1, -1/3, -1/3, -1/3): distance 1/3 (Eckart-Young)
    assert err > 1 / 3
    assert torch.linalg.matrix_norm(R[1:, 1:] - TARGET[1:, 1:], ord=2) > 1 / 3


def test_one_stinespring_layer_fits_the_pauli_channel():
    torch.manual_seed(1)
    alg = TypeIAlgebra([2], [2], precision='double')
    R = _fit(learn.ChannelLayer(alg, kraus_rank=4), _probes(alg), 400)
    assert torch.linalg.matrix_norm(R - TARGET) < 1e-2                    # far below the residual floor 1/3
    assert torch.linalg.det(R) < -0.03                                     # target: -1/27


def test_damped_target_is_not_a_semigroup_on_any_branch():
    """Target (ii) of the study fails the test of Wolf, Eisert, Cubitt and Cirac on every branch of its logarithm."""
    import importlib.util
    import sys
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root / 'scripts'))
    spec = importlib.util.spec_from_file_location('residual_markov', root / 'scripts' / 'residual_markov.py')
    rm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rm)
    maps, _ = rm.target_maps(7)
    damped, gksl = rm.semigroup_test(maps['damped']), rm.semigroup_test(maps['gksl'])
    assert damped['all_branches_fail'] and not damped['degenerate']
    assert damped['best'] < -0.2
    assert gksl['branches'][0]['min_eig'] > -1e-10                 # the random GKSL target passes on k = 0
