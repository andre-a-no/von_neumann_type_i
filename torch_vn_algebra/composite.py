"""
Tensor products of Type I algebras.

Plain product:

    ((+)_c M_{k_c}) (x) ((+)_d M_{m_d}) = (+)_{(c, d)} M_{k_c m_d},

again of Type I with one sector for every pair (c, d), stored in the order s = c * D + d.

Fused (symmetric) product: when the sectors carry conserved charges (e.g. particle numbers),
the operators of the composite system that conserve the *total* charge form the larger algebra

    (+)_q M_{n_q},     n_q = sum_{fuse(c, d) = q} k_c m_d,

in which all pairs with the same total charge q are merged into one sector (fused_tensor_product).
It contains the plain product: A (x) B sits block-diagonally inside every merged sector, but
operators such as hopping terms S+_i S-_j that move charge between the two factors only exist in
the fused algebra.

In both cases sector s is a list of pair blocks (c, d, offset): basis vector offset + a * m_d + b of
sector s is e_a (x) e_b with e_a in sector c of the first and e_b in sector d of the second factor.
Operators factorise blockwise, and the partial traces are maps between sectors,
implemented as InterSectorChannels.
"""
import operator
from typing import Callable, List, Tuple

import torch

from .algebra import TypeIAlgebra
from . import cost
from .channels import InterSectorChannel
from .states import DensityMatrix

Operator = TypeIAlgebra.Operator


def _make_product(alg1: TypeIAlgebra, alg2: TypeIAlgebra, groups: List[Tuple], charges: List) -> TypeIAlgebra:
    if alg1.hilbert.complex_valued != alg2.hilbert.complex_valued:
        raise ValueError("both factors must be real or both complex")
    blocks, n, k = [], [], []
    for pairs in groups:
        off, n_s, sec = 0, 0, []
        for c, d in pairs:
            sec.append((c, d, off))
            off += alg1.k_factors[c] * alg2.k_factors[d]
            n_s += alg1.n_factors[c] * alg2.n_factors[d]
        blocks.append(sec)
        k.append(off)
        n.append(max(n_s, 1))
    alg = TypeIAlgebra(n, k, complex_valued=alg1.hilbert.complex_valued, device=alg1.hilbert.device,
                       precision=alg1.hilbert.precision, charges=charges)
    alg.tensor_factors = (alg1, alg2)
    alg.product_blocks = blocks
    return alg


def tensor_product(alg1: TypeIAlgebra, alg2: TypeIAlgebra) -> TypeIAlgebra:
    """The algebra alg1 (x) alg2 with sectors (c, d) -> c * alg2.C + d and charges (q_c, q_d)."""
    groups = [[(c, d)] for c in range(alg1.C) for d in range(alg2.C)]
    charges = [(q1, q2) for q1 in alg1.charges for q2 in alg2.charges]
    return _make_product(alg1, alg2, groups, charges)


def fused_tensor_product(alg1: TypeIAlgebra, alg2: TypeIAlgebra,
                         fuse: Callable = operator.add) -> TypeIAlgebra:
    """
    Charge-conserving product: pairs (c, d) with equal fuse(q_c, q_d) are merged into one sector.
    Sectors are ordered by the fused charge (sorted when the charges are comparable), the pairs
    inside a sector by (c, d). Default fusion rule: addition (U(1) charges).
    """
    groups = {}
    for c, q1 in enumerate(alg1.charges):
        for d, q2 in enumerate(alg2.charges):
            groups.setdefault(fuse(q1, q2), []).append((c, d))
    try:
        keys = sorted(groups)
    except TypeError:
        keys = list(groups)
    return _make_product(alg1, alg2, [groups[q] for q in keys], keys)


def _factors(alg12: TypeIAlgebra) -> Tuple[TypeIAlgebra, TypeIAlgebra]:
    factors = getattr(alg12, 'tensor_factors', None)
    if factors is None:
        raise ValueError("algebra was not created by tensor_product or fused_tensor_product")
    return factors


def sector_index(alg12: TypeIAlgebra, c: int, d: int) -> int:
    """Index of the sector of alg12 that contains the pair (c, d)."""
    _factors(alg12)
    for s, pairs in enumerate(alg12.product_blocks):
        if any(p[0] == c and p[1] == d for p in pairs):
            return s
    raise KeyError((c, d))


def pair_offset(alg12: TypeIAlgebra, c: int, d: int) -> Tuple[int, int]:
    """(sector, offset) of the pair block (c, d) in alg12."""
    s = sector_index(alg12, c, d)
    return s, next(off for cc, dd, off in alg12.product_blocks[s] if cc == c and dd == d)


def kron(A: Operator, B: Operator, alg12: TypeIAlgebra) -> Operator:
    """A (x) B in alg12 = tensor_product / fused_tensor_product of (A.algebra, B.algebra) (batched,
    lazy). In a fused product A (x) B is block diagonal inside the merged sectors. States give a state."""
    alg1, alg2 = _factors(alg12)
    assert A.algebra is alg1 and B.algebra is alg2

    def generator():
        a, b = A.matrix, B.matrix
        dtype = torch.promote_types(a.dtype, b.dtype)
        a, b = a.to(dtype), b.to(dtype)
        batch = max(a.shape[0], b.shape[0])
        cost.check_memory(3 * cost.tensor_bytes((batch, alg12.C, alg12.k_max, alg12.k_max), dtype), a.device,
                          f"kron (batch {batch}, {alg12.C} sectors of size up to {alg12.k_max})")
        ci, ai, aj, di, bi, bj, mask = _kron_index(alg12)
        # one gather per factor: (A (x) B)[s, i, j] = A[c, a_i, a_j] B[d, b_i, b_j] inside a pair block, else 0
        return a[:, ci, ai, aj] * b[:, di, bi, bj] * mask.to(dtype)

    if isinstance(A, DensityMatrix) and isinstance(B, DensityMatrix):
        return DensityMatrix(alg12, generator=generator, validate=False)
    out = Operator(alg12, generator=generator)
    if A._is_self_adjoint and B._is_self_adjoint:
        out._is_self_adjoint = True
    if A._is_positive and B._is_positive:
        out._is_positive = True
    return out


def _kron_index(alg12: TypeIAlgebra):
    """Gather maps of kron, built once per product algebra: for every entry (s, i, j) of alg12 the sector and
    indices of the factors (c, a_i, a_j), (d, b_i, b_j), and a mask that is 1 inside the pair blocks."""
    key = ('kron_index',)
    idx = alg12._const_cache.get(key)
    if idx is not None:
        return idx
    alg1, alg2 = _factors(alg12)
    C, K = alg12.C, alg12.k_max
    pair = torch.full((C, K), -1, dtype=torch.long)                      # pair id of every row
    c_of, d_of, a_of, b_of = (torch.zeros(C, K, dtype=torch.long) for _ in range(4))
    pid = 0
    for s, pairs in enumerate(alg12.product_blocks):
        for c, d, off in pairs:
            k_c, m_d = alg1.k_factors[c], alg2.k_factors[d]
            if k_c * m_d == 0:
                continue
            r = torch.arange(k_c * m_d)
            pair[s, off:off + k_c * m_d] = pid
            c_of[s, off:off + k_c * m_d] = c
            d_of[s, off:off + k_c * m_d] = d
            a_of[s, off:off + k_c * m_d] = r // m_d
            b_of[s, off:off + k_c * m_d] = r % m_d
            pid += 1
    mask = (pair[:, :, None] == pair[:, None, :]) & (pair[:, :, None] >= 0)
    dev = alg12.hilbert.device
    idx = tuple(t.to(dev) for t in (c_of[:, :, None].expand(C, K, K), a_of[:, :, None].expand(C, K, K),
                                    a_of[:, None, :].expand(C, K, K), d_of[:, :, None].expand(C, K, K),
                                    b_of[:, :, None].expand(C, K, K), b_of[:, None, :].expand(C, K, K), mask))
    alg12._const_cache[key] = idx
    return idx


def _ptrace_index(alg12: TypeIAlgebra, keep: int):
    """Per group of pair blocks with equal sizes (k_c, m_d): (k, m, sector s, row indices (P, k m), output sector)
    for the direct partial trace; built once per product algebra."""
    key = ('ptrace_index', keep)
    groups = alg12._const_cache.get(key)
    if groups is not None:
        return groups
    alg1, alg2 = _factors(alg12)
    by_size = {}
    for s, pairs in enumerate(alg12.product_blocks):
        for c, d, off in pairs:
            k_c, m_d = alg1.k_factors[c], alg2.k_factors[d]
            if k_c * m_d:
                by_size.setdefault((k_c, m_d), []).append((s, off, c if keep == 1 else d))
    dev = alg12.hilbert.device
    groups = []
    for (k_c, m_d), entries in by_size.items():
        S = torch.tensor([e[0] for e in entries], dtype=torch.long, device=dev)
        rows = torch.tensor([e[1] for e in entries], dtype=torch.long, device=dev)[:, None] + \
            torch.arange(k_c * m_d, device=dev)[None, :]
        out = torch.tensor([e[2] for e in entries], dtype=torch.long, device=dev)
        groups.append((k_c, m_d, S, rows, out))
    alg12._const_cache[key] = groups
    return groups


def partial_trace_channel(alg12: TypeIAlgebra, keep: int = 1) -> InterSectorChannel:
    """
    Partial trace alg12 -> alg1 (keep=1, trace out the second factor) or -> alg2 (keep=2), as a
    trace-preserving map between sectors with Kraus operators 1 (x) <j| (resp. <j| (x) 1) on every pair
    block. Works for plain and fused products (for additive charges, blocks of a fused operator
    between different pairs never contribute). Its blunt-trace adjoint is the embedding A -> A (x) 1.
    """
    alg1, alg2 = _factors(alg12)
    if keep not in (1, 2):
        raise ValueError("keep must be 1 or 2")
    out_alg = alg1 if keep == 1 else alg2
    traced = alg2 if keep == 1 else alg1
    r_dim = max(traced.k_factors)
    # several pairs of one merged sector may map to the same output sector: separate Kraus slots
    slots = 1
    for pairs in alg12.product_blocks:
        counts = {}
        for c, d, _ in pairs:
            key = c if keep == 1 else d
            counts[key] = counts.get(key, 0) + 1
        slots = max(slots, max(counts.values(), default=1))
    K = torch.zeros(1, slots * r_dim, out_alg.C, alg12.C, out_alg.k_max, alg12.k_max,
                    dtype=alg12.hilbert.dtype, device=alg12.hilbert.device)
    dev = K.device
    for s, pairs in enumerate(alg12.product_blocks):
        used = {}
        for c, d, off in pairs:
            k_c, m_d = alg1.k_factors[c], alg2.k_factors[d]
            if k_c * m_d == 0:
                continue
            out_sec = c if keep == 1 else d
            slot = used.get(out_sec, 0)
            used[out_sec] = slot + 1
            if keep == 1:                                   # K_j = 1_{k_c} (x) <j|
                a = torch.arange(k_c, device=dev)[:, None]
                j = torch.arange(m_d, device=dev)[None, :]
                K[0, slot * r_dim + j, c, s, a, off + a * m_d + j] = 1.0
            else:                                           # K_j = <j| (x) 1_{m_d}
                j = torch.arange(k_c, device=dev)[:, None]
                b = torch.arange(m_d, device=dev)[None, :]
                K[0, slot * r_dim + j, d, s, b, off + j * m_d + b] = 1.0
    return InterSectorChannel(alg12, out_alg, K)


def partial_trace(X: Operator, keep: int = 1) -> Operator:
    """
    Tr_2 X (keep=1) or Tr_1 X (keep=2) for X in a tensor-product algebra, computed directly: every pair block
    (c, d) is reshaped to (k_c, m_d, k_c, m_d) and traced over one factor, one gather per group of equal block
    sizes (O(B sum k^2 m) operations; partial_trace_channel gives the same map as a channel, e.g. for its
    adjoint). States give states.
    """
    alg12 = X.algebra
    alg1, alg2 = _factors(alg12)
    if keep not in (1, 2):
        raise ValueError("keep must be 1 or 2")
    out_alg = alg1 if keep == 1 else alg2

    def generator():
        x = X.matrix
        B = x.shape[0]
        out = torch.zeros(B, out_alg.C, out_alg.k_max, out_alg.k_max, dtype=x.dtype, device=x.device)
        for k_c, m_d, S, rows, osec in _ptrace_index(alg12, keep):
            blk = x[:, S[:, None, None], rows[:, :, None], rows[:, None, :]]              # (B, P, km, km)
            blk = blk.reshape(B, len(S), k_c, m_d, k_c, m_d)
            red = blk.diagonal(dim1=3, dim2=5).sum(-1) if keep == 1 else blk.diagonal(dim1=2, dim2=4).sum(-1)
            n = red.shape[-1]
            pad = torch.zeros(B, len(S), out_alg.k_max, out_alg.k_max, dtype=x.dtype, device=x.device)
            pad[:, :, :n, :n] = red
            out.index_add_(1, osec, pad)
        return out

    if isinstance(X, DensityMatrix):
        return DensityMatrix(out_alg, generator=generator, validate=False)
    out = Operator(out_alg, generator=generator)
    out._is_self_adjoint = X._is_self_adjoint
    out._is_positive = X._is_positive
    return out
