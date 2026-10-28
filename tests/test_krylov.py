"""Tests for sparse sector Hamiltonians and Krylov methods."""
import math
import warnings

import pytest
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


def test_real_start_vector_keeps_complex_hamiltonian():
    # a real v0 must not drop Im H (momentum blocks are complex)
    from torch_vn_algebra.krylov import SparseSectorHamiltonian
    ch = SpinChain(8, 'periodic')
    Hm = ch.xxz_momentum(1.0, 1.0, 0.0, sector=4)
    c = Hm.algebra.charges.index(1)
    k = Hm.algebra.k_factors[c]
    B = Hm.matrix[0, c, :k, :k]
    H = SparseSectorHamiltonian(B.diagonal().real[None], (B - torch.diag(B.diagonal())).to_sparse(), torch.complex128)
    E, _, _ = krylov.ground_state(H, v0=torch.randn(1, k, dtype=torch.float64))
    assert abs(E.item() - torch.linalg.eigvalsh(B)[0].item()) < 1e-9
    x = torch.randn(1, k, dtype=torch.float64)
    assert torch.allclose(H.matvec(x), H.matvec(x.to(torch.complex128)))


def test_single_precision_converges_without_warning():
    import warnings
    ch = SpinChain(10, 'periodic', precision='single')
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        E, psi, _ = krylov.ground_state(ch.xxz_sparse(1.0, 1.0, sector=5))
    assert psi.dtype == torch.complex64 and abs(E.item() + 4.515446354) < 1e-4


# ----------------------------------------------------------------------------------------------------------------
# Matrix-free sector operators (protocol SectorOperator): functions, dense blocks, sparse matrices, combinations
# ----------------------------------------------------------------------------------------------------------------
def _herm(B, k, seed):
    g = torch.Generator().manual_seed(seed)
    X = torch.randn(B, k, k, dtype=torch.complex128, generator=g)
    return X + X.mH


def test_matrix_free_dense_and_sparse_operators_agree():
    from torch_vn_algebra import krylov, TypeIAlgebra
    A = _herm(2, 30, 0)
    ops = [krylov.MatrixFreeOperator(lambda V: (A @ V.unsqueeze(-1)).squeeze(-1), 30, 2),
           krylov.MatrixFreeOperator(lambda V: A @ V, 30, 2, batched_columns=True),
           krylov.DenseSectorOperator(A)]
    for H in ops:
        assert torch.allclose(H.to_dense(), A, atol=1e-12)
        E, psi, res = krylov.ground_state(H, generator=torch.Generator().manual_seed(1))
        assert torch.allclose(E, torch.linalg.eigvalsh(A)[:, 0], atol=1e-9) and (res < 1e-8).all()
        assert torch.allclose(H.expectation(psi).real, E, atol=1e-9)
    S = (A[0] * (A[0].abs() > 1.5)).to_sparse()                     # a shared sparse part plus batched diagonals
    d = torch.randn(3, 30, dtype=torch.float64)
    H = krylov.SparseSectorHamiltonian.from_sparse(S, diag=d)
    assert torch.allclose(H.to_dense(), torch.diag_embed(d).to(S.dtype) + S.to_dense(), atol=1e-12)
    i, j = S.coalesce().indices()
    H2 = krylov.SparseSectorHamiltonian.from_coo(i, j, S.coalesce().values(), 30, diag=d)
    assert torch.allclose(H2.to_dense(), H.to_dense())
    alg = TypeIAlgebra([3, 30], [3, 30], precision='double')
    op = alg.from_blocks([_herm(1, 3, 2)[0], A[0]])
    assert torch.allclose(krylov.from_operator(op, 1).to_dense()[0], A[0])


def test_combinations_with_per_sample_coefficients():
    from torch_vn_algebra import krylov
    A, B = _herm(1, 20, 3), _herm(1, 20, 4)
    HA, HB = krylov.DenseSectorOperator(A), krylov.DenseSectorOperator(B)
    c = torch.tensor([0.0, 0.5, 2.0], dtype=torch.float64)
    H = HA * 2.0 + HB * c - HA
    assert H.batch_size == 3 and isinstance(H, krylov.LinearCombination) and len(H.terms) == 3
    dense = A + c[:, None, None] * B
    assert torch.allclose(H.to_dense(), dense.to(H.dtype), atol=1e-12)
    E, _, _ = krylov.ground_state(H, generator=torch.Generator().manual_seed(0))
    assert torch.allclose(E, torch.linalg.eigvalsh(dense)[:, 0], atol=1e-9)
    with pytest.raises(ValueError):
        krylov.combine([(1.0, HA), (1.0, krylov.DenseSectorOperator(_herm(1, 5, 0)))])


def test_lowest_eigenvalues_by_deflation():
    from torch_vn_algebra import krylov, SpinChain
    A = _herm(2, 40, 5)
    E, psi = krylov.lowest_eigenvalues(krylov.DenseSectorOperator(A), levels=4,
                                       generator=torch.Generator().manual_seed(0))
    assert torch.allclose(E, torch.linalg.eigvalsh(A)[:, :4], atol=1e-8)
    ov = torch.einsum('bik,bjk->bij', psi.conj(), psi)
    assert torch.allclose(ov, torch.eye(4, dtype=ov.dtype).expand(2, -1, -1), atol=1e-6)
    ch = SpinChain(10, 'periodic')                                  # degenerate levels come with their multiplicity
    H = ch.xxz_sparse(Delta=1.0, sector=5)
    E, _ = krylov.lowest_eigenvalues(H, levels=5, generator=torch.Generator().manual_seed(1))
    assert torch.allclose(E[0], torch.linalg.eigvalsh(ch.xxz(Delta=1.0, sector=5).matrix[0, 0])[:5], atol=1e-8)
