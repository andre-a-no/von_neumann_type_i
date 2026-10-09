"""
Tensor products of Type I algebras.

    ((+)_c M_{k_c}) (x) ((+)_d M_{m_d}) = (+)_{(c, d)} M_{k_c m_d},

so the product is again of Type I, with one sector for every pair (c, d), stored in the order
s = c * D + d. Operators factorise blockwise, (A (x) B)_{(c,d)} = A_c (x) B_d, and the partial
traces are maps between sectors, (Tr_2 X)_c = sum_d Tr_{m_d} X_{(c,d)}, implemented as
InterSectorChannels.
"""
from typing import Tuple

import torch

from .algebra import TypeIAlgebra
from .channels import InterSectorChannel
from .states import DensityMatrix

Operator = TypeIAlgebra.Operator


def tensor_product(alg1: TypeIAlgebra, alg2: TypeIAlgebra) -> TypeIAlgebra:
    """The algebra alg1 (x) alg2 with sectors (c, d) -> c * alg2.C + d. Keeps a reference to the factors."""
    if alg1.hilbert.complex_valued != alg2.hilbert.complex_valued:
        raise ValueError("both factors must be real or both complex")
    n = [n1 * n2 for n1 in alg1.n_factors for n2 in alg2.n_factors]
    k = [k1 * k2 for k1 in alg1.k_factors for k2 in alg2.k_factors]
    alg = TypeIAlgebra(n, k, complex_valued=alg1.hilbert.complex_valued, device=alg1.hilbert.device)
    alg.tensor_factors = (alg1, alg2)
    return alg


def _factors(alg12: TypeIAlgebra) -> Tuple[TypeIAlgebra, TypeIAlgebra]:
    factors = getattr(alg12, 'tensor_factors', None)
    if factors is None:
        raise ValueError("algebra was not created by tensor_product")
    return factors


def sector_index(alg12: TypeIAlgebra, c: int, d: int) -> int:
    """Index of the sector (c, d) in the product algebra."""
    return c * _factors(alg12)[1].C + d


def kron(A: Operator, B: Operator, alg12: TypeIAlgebra) -> Operator:
    """A (x) B in alg12 = tensor_product(A.algebra, B.algebra) (batched, lazy). States give a state."""
    alg1, alg2 = _factors(alg12)
    assert A.algebra is alg1 and B.algebra is alg2

    def generator():
        a, b = A.matrix, B.matrix
        dtype = torch.promote_types(a.dtype, b.dtype)
        a, b = a.to(dtype), b.to(dtype)
        batch = max(a.shape[0], b.shape[0])
        out = torch.zeros(batch, alg12.C, alg12.k_max, alg12.k_max, dtype=dtype, device=a.device)
        for c, k_c in enumerate(alg1.k_factors):
            for d, m_d in enumerate(alg2.k_factors):
                if k_c * m_d == 0:
                    continue
                ac, bd = a[:, c, :k_c, :k_c], b[:, d, :m_d, :m_d]
                blk = (ac[:, :, None, :, None] * bd[:, None, :, None, :]).reshape(-1, k_c * m_d, k_c * m_d)
                out[:, c * alg2.C + d, :k_c * m_d, :k_c * m_d] = blk
        return out

    if isinstance(A, DensityMatrix) and isinstance(B, DensityMatrix):
        return DensityMatrix(alg12, generator=generator, validate=False)
    out = Operator(alg12, generator=generator)
    if A._is_self_adjoint and B._is_self_adjoint:
        out._is_self_adjoint = True
    if A._is_positive and B._is_positive:
        out._is_positive = True
    return out


def partial_trace_channel(alg12: TypeIAlgebra, keep: int = 1) -> InterSectorChannel:
    """
    Partial trace alg12 -> alg1 (keep=1, trace out the second factor) or -> alg2 (keep=2), as a
    trace-preserving map between sectors with Kraus operators 1 (x) <j| (resp. <j| (x) 1).
    Its blunt-trace adjoint is the embedding A -> A (x) 1 (resp. 1 (x) A).
    """
    alg1, alg2 = _factors(alg12)
    if keep not in (1, 2):
        raise ValueError("keep must be 1 or 2")
    out_alg = alg1 if keep == 1 else alg2
    traced = alg2 if keep == 1 else alg1
    r = max(traced.k_factors)
    K = torch.zeros(1, r, out_alg.C, alg12.C, out_alg.k_max, alg12.k_max,
                    dtype=alg12.hilbert.dtype, device=alg12.hilbert.device)
    for c, k_c in enumerate(alg1.k_factors):
        for d, m_d in enumerate(alg2.k_factors):
            s = c * alg2.C + d
            if keep == 1:
                for j in range(m_d):
                    for a in range(k_c):
                        K[0, j, c, s, a, a * m_d + j] = 1.0          # 1_{k_c} (x) <j|
            else:
                for j in range(k_c):
                    for b in range(m_d):
                        K[0, j, d, s, b, j * m_d + b] = 1.0          # <j| (x) 1_{m_d}
    return InterSectorChannel(alg12, out_alg, K)


def partial_trace(X: Operator, keep: int = 1) -> Operator:
    """Tr_2 X (keep=1) or Tr_1 X (keep=2) for X in a tensor-product algebra."""
    return partial_trace_channel(X.algebra, keep)(X)
