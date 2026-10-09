"""Tests for sparse sector Hamiltonians and Krylov methods."""
import math
import warnings

import torch

from torch_vn_algebra import SpinChain, dynamics, krylov


def test_sparse_equals_dense():
    ch = SpinChain(8, 'periodic')
    h = torch.randn(3, 8, dtype=torch.float64)
    Hs = ch.xxz_sparse(1.0, 0.6, h, sector=4)
    Hd = ch.xxz(1.0, 0.6, h, sector=4).matrix[:, 0]
    assert torch.allclose(Hs.to_dense(), Hd)
    v = torch.randn(3, Hs.dim, dtype=torch.complex128)
    assert torch.allclose(Hs.matvec(v), (Hd @ v.unsqueeze(-1)).squeeze(-1))


def test_batched_lanczos_ground_state():
    ch = SpinChain(10, 'periodic')
    h = torch.randn(4, 10, dtype=torch.float64)
    E, psi, res = krylov.ground_state(ch.xxz_sparse(1.0, 1.0, h, sector=5))
    exact = torch.linalg.eigvalsh(ch.xxz(1.0, 1.0, h, sector=5).matrix[:, 0])[:, 0]
    assert torch.allclose(E, exact, atol=1e-10)
    assert torch.all(res < 1e-8)


def test_lanczos_heisenberg_ring_l16():
    ch = SpinChain(16, 'periodic', complex_valued=False)
    E, _, _ = krylov.ground_state(ch.xxz_sparse(sector=8))
    assert abs(E.item() + 7.142296361) < 1e-8


def test_krylov_evolution_matches_dense():
    ch = SpinChain(8)
    psi0 = ch.vector_in_sector('11110000')
    ts = [0.0, 0.7, 2.0, 4.0]
    dense = dynamics.schrodinger(psi0[None, None], ch.xxz(1.0, 0.5, sector=4), ts)[:, 0, 0]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        kry = krylov.evolve(ch.xxz_sparse(1.0, 0.5, sector=4), psi0[None], ts)[:, 0]
    assert torch.allclose(dense, kry, atol=1e-10)
    assert torch.allclose(torch.linalg.vector_norm(kry, dim=-1), torch.ones(4, dtype=torch.float64), atol=1e-12)
