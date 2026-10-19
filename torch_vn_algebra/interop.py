"""
Conversion between elements of M = (+)_c M_{k_c} and dense matrices, QuTiP objects and dynamiqs arrays.

An element of M becomes the block-diagonal D x D matrix, D = sum_c k_c, with the sectors in their order. A density
stays a density: the library stores states with respect to Tr_blunt, which is the ordinary trace of that matrix. A
map between algebras (InterSectorChannel) is a map into the algebra M_out, i.e. it ends with the pinching onto the
output sectors; its dense Kraus operators are therefore P_d K_i, one for every Kraus operator K_i and output sector d
(the sum over d of P_d K_i alone would create coherences between output sectors). A SectorJump, whose sectors each go
to a single sector, becomes the dense jump operator L itself.

When the algebra acts on a Hilbert space with a basis of its own, `layout` places the basis vectors there: either a
SpinChain (sector N, basis state with label x -> the computational basis of L qubits in QuTiP's convention: site 0
is the leftmost tensor factor and |0> is spin up, so S^z_i = sigmaz()/2 on site i and S^-_i = sigmam()), or a pair
(basis, dims) with basis a LongTensor of length D (position of every basis vector, sector after sector) and dims the
QuTiP dims of the target space.

    rho_q = interop.to_qutip(rho, layout=chain)            # Qobj (a list of Qobj for a batch)
    H_q = interop.to_qutip(H, layout=chain)
    c_ops = interop.jumps_to_qutip(jumps, layout=chain)    # SectorJump / InterSectorChannel / Operator
    H_d, L_d = interop.to_dynamiqs(H, layout=chain), interop.jumps_to_dynamiqs(jumps, layout=chain)
    rho = interop.from_qutip(rho_q, chain.algebra, layout=chain, state=True)
    D = interop.to_dense(A)                                 # torch (batch, D, D)

QuTiP and dynamiqs are optional; the functions that need them import them on use.
"""
from typing import List, Optional, Sequence, Tuple, Union

import torch

from .algebra import TypeIAlgebra
from .channels import InterSectorChannel, SectorJump
from .states import DensityMatrix

Operator = TypeIAlgebra.Operator


# ----------------------------------------------------------------------
# layouts
# ----------------------------------------------------------------------
def chain_basis(chain) -> torch.Tensor:
    """Positions of the basis states of a SpinChain (sector after sector) in the computational basis of L qubits,
    QuTiP's convention: index = sum_i (1 - b_i) 2^(L-1-i) for the spins b_i (1 = up) of the label."""
    L = chain.L
    lab = torch.cat(list(chain.labels))
    bits = (lab[:, None] >> torch.arange(L, device=lab.device)) & 1
    weights = 2 ** torch.arange(L - 1, -1, -1, device=lab.device)
    return ((1 - bits) * weights).sum(-1)


def _layout(algebra: TypeIAlgebra, layout) -> Tuple[torch.Tensor, int, Optional[list]]:
    """(positions (D,), dimension of the target space, QuTiP dims or None)."""
    D = sum(algebra.k_factors)
    if layout is None:
        return torch.arange(D), D, None
    if hasattr(layout, 'labels') and hasattr(layout, 'L'):             # a SpinChain
        if list(algebra.k_factors) != [len(x) for x in layout.labels]:
            raise ValueError("the SpinChain layout fits only the chain's full algebra (+)_N M_binom(L, N)")
        return chain_basis(layout).cpu(), 2 ** layout.L, [[2] * layout.L, [2] * layout.L]
    basis, dims = layout
    basis = torch.as_tensor(basis, dtype=torch.long).cpu()
    if basis.shape != (D,):
        raise ValueError(f"basis must have length D = {D}, got shape {tuple(basis.shape)}")
    if len(torch.unique(basis)) != D:
        raise ValueError("basis positions must be distinct")
    n = 1
    if dims is not None:
        for d in dims[0]:
            n *= d
    else:
        n = int(basis.max()) + 1
    if int(basis.max()) >= n:
        raise ValueError(f"basis position {int(basis.max())} outside a space of dimension {n}")
    return basis, n, dims


def _positions(algebra: TypeIAlgebra, basis: torch.Tensor) -> List[torch.Tensor]:
    out, start = [], 0
    for k in algebra.k_factors:
        out.append(basis[start:start + k])
        start += k
    return out


# ----------------------------------------------------------------------
# dense matrices
# ----------------------------------------------------------------------
def to_dense(x: Operator, layout=None) -> torch.Tensor:
    """(batch, n, n) matrix of an element of M (block diagonal for layout=None, n = D)."""
    alg = x.algebra
    basis, n, _ = _layout(alg, layout)
    M = x.matrix
    out = M.new_zeros(M.shape[0], n, n)
    for c, (k, pos) in enumerate(zip(alg.k_factors, _positions(alg, basis.to(M.device)))):
        if k:
            out[:, pos[:, None], pos[None, :]] = M[:, c, :k, :k]
    return out


def channel_to_dense(channel: Union[InterSectorChannel, SectorJump], layout_in=None, layout_out=None,
                     per_output_sector: Optional[bool] = None) -> torch.Tensor:
    """
    Kraus operators of a map between algebras as dense matrices, shape (batch, r', n_out, n_in). By default an
    InterSectorChannel gives r' = r * D operators P_d K_i (the dense form of the map into M_out, see the module
    docstring) and a SectorJump gives its operator L (r' = 1). per_output_sector=False returns the r operators
    sum_d P_d K_i, which agree with the channel on block-diagonal inputs only if every input sector goes to one
    output sector.
    """
    if per_output_sector is None:
        per_output_sector = not isinstance(channel, SectorJump)
    if isinstance(channel, SectorJump):
        channel = channel.to_inter_sector()
    a_in, a_out = channel.algebra_in, channel.algebra_out
    b_in, n_in, _ = _layout(a_in, layout_in)
    b_out, n_out, _ = _layout(a_out, layout_out)
    K = channel.kraus
    out = K.new_zeros(K.shape[0], K.shape[1], n_out, n_in)
    p_in, p_out = _positions(a_in, b_in.to(K.device)), _positions(a_out, b_out.to(K.device))
    for d, (m, po) in enumerate(zip(a_out.k_factors, p_out)):
        for c, (k, pi) in enumerate(zip(a_in.k_factors, p_in)):
            if m and k:
                out[:, :, po[:, None], pi[None, :]] = K[:, :, d, c, :m, :k]
    if not per_output_sector:
        return out
    pieces = []
    for m, po in zip(a_out.k_factors, p_out):                          # P_d K_i for every output sector d
        if m:
            P = torch.zeros(n_out, dtype=torch.bool, device=K.device)
            P[po] = True
            pieces.append(out * P[:, None])
    return torch.stack(pieces, 2).reshape(K.shape[0], -1, n_out, n_in)


def from_dense(matrix, algebra: TypeIAlgebra, layout=None, state: bool = False, atol: float = 1e-10):
    """
    The element of M with the dense matrix `matrix` ((n, n) or (batch, n, n); torch or NumPy). Raises if the matrix
    has weight outside the blocks of M (more than atol times its norm), i.e. if it is not in the algebra. state=True
    returns a DensityMatrix (the matrix is then the density with respect to the trace of the target space).
    """
    basis, n, _ = _layout(algebra, layout)
    A = torch.as_tensor(matrix)
    if A.dim() == 2:
        A = A.unsqueeze(0)
    if tuple(A.shape[1:]) != (n, n):
        raise ValueError(f"expected matrices of size {n} x {n}, got {tuple(A.shape[1:])}")
    A = A.to(device=algebra.hilbert.device, dtype=algebra.dtype)
    blocks = A.new_zeros(A.shape[0], algebra.C, algebra.k_max, algebra.k_max)
    inside = A.new_zeros(A.shape, dtype=torch.bool)
    for c, (k, pos) in enumerate(zip(algebra.k_factors, _positions(algebra, basis.to(A.device)))):
        if k:
            blocks[:, c, :k, :k] = A[:, pos[:, None], pos[None, :]]
            inside[:, pos[:, None], pos[None, :]] = True
    outside = (A * ~inside).abs().amax().item()
    if outside > atol * max(A.abs().amax().item(), 1.0):
        raise ValueError(f"the matrix is not in the algebra: entries of size {outside:.1e} outside its blocks")
    if state:
        return DensityMatrix(algebra, matrix=blocks)
    return algebra.operator(blocks)


# ----------------------------------------------------------------------
# QuTiP
# ----------------------------------------------------------------------
def _qutip():
    try:
        import qutip
    except ImportError as e:                                   # pragma: no cover - depends on the environment
        raise ImportError("this function needs QuTiP (pip install qutip)") from e
    return qutip


def _qobjs(dense: torch.Tensor, dims, batch_index):
    qutip = _qutip()
    arr = dense.detach().cpu().numpy()
    if batch_index is not None:
        return qutip.Qobj(arr[batch_index], dims=dims)
    qs = [qutip.Qobj(a, dims=dims) for a in arr]
    return qs[0] if len(qs) == 1 else qs


def to_qutip(x: Operator, layout=None, batch_index: Optional[int] = None):
    """QuTiP Qobj of an element of M (a list of Qobj for a batch, unless batch_index picks one sample)."""
    _, _, dims = _layout(x.algebra, layout)
    return _qobjs(to_dense(x, layout), dims, batch_index)


def jumps_to_qutip(jumps: Sequence, layout=None, rates: Optional[Sequence[float]] = None,
                   batch_index: int = 0) -> list:
    """
    Collapse operators for qutip.mesolve from the jump operators of dynamics.lindblad_evolve: an Operator (a jump
    inside the sectors), a SectorJump or an InterSectorChannel (each of its Kraus operators is one jump). With
    `rates`, the operators are multiplied by sqrt(rate) as lindblad_evolve does; a batch of coefficients (e.g.
    SectorJump.scaled) is resolved by batch_index.
    """
    out = []
    for j, J in enumerate(jumps):
        s = 1.0 if rates is None else float(rates[j]) ** 0.5
        if isinstance(J, (SectorJump, InterSectorChannel)):
            alg = J.algebra if isinstance(J, SectorJump) else J.algebra_in
            _, _, dims = _layout(alg, layout)
            K = channel_to_dense(J, layout, layout)
            b = min(batch_index, K.shape[0] - 1)
            out += [_qobjs(s * K[b, i][None], dims, 0) for i in range(K.shape[1])]
        else:
            q = to_qutip(J, layout, batch_index=min(batch_index, J.matrix.shape[0] - 1))
            out.append(s * q)
    return out


def from_qutip(q, algebra: TypeIAlgebra, layout=None, state: Optional[bool] = None, atol: float = 1e-10):
    """
    The element of M of a QuTiP operator, or the state of a QuTiP ket or density matrix. state=None treats kets as
    states and operators as operators; state=True turns an operator into a DensityMatrix. A list of Qobj gives a batch.
    """
    qs = q if isinstance(q, (list, tuple)) else [q]
    mats = []
    is_state = state
    for qi in qs:
        if qi.isket:
            v = torch.as_tensor(qi.full()).reshape(-1)
            mats.append(v[:, None] * v.conj()[None, :])
            is_state = True if state is None else state
        else:
            mats.append(torch.as_tensor(qi.full()))
    return from_dense(torch.stack(mats), algebra, layout, state=bool(is_state), atol=atol)


# ----------------------------------------------------------------------
# dynamiqs
# ----------------------------------------------------------------------
def jumps_to_dynamiqs(jumps: Sequence, layout=None, rates: Optional[Sequence[float]] = None,
                      batch_index: int = 0) -> list:
    """Jump operators for dynamiqs.mesolve, as jumps_to_qutip (one array per Kraus operator, sqrt(rate) folded in)."""
    out = []
    for j, J in enumerate(jumps):
        s = 1.0 if rates is None else float(rates[j]) ** 0.5
        a = to_dynamiqs(J, layout)
        if isinstance(J, (SectorJump, InterSectorChannel)):
            b = min(batch_index, a.shape[0] - 1)
            out += [s * a[b, i] for i in range(a.shape[1])]
        else:
            out.append(s * a[min(batch_index, a.shape[0] - 1)])
    return out


def to_dynamiqs(x, layout=None):
    """dynamiqs array (batch, n, n) of an element of M, or of the Kraus operators (batch, r, n, n) of a SectorJump
    or InterSectorChannel; the batch dimension is kept, so a batch of Hamiltonians or states runs as one dynamiqs
    batch."""
    try:
        import dynamiqs as dq
    except ImportError as e:                                   # pragma: no cover - depends on the environment
        raise ImportError("this function needs dynamiqs (pip install dynamiqs)") from e
    if isinstance(x, (SectorJump, InterSectorChannel)):
        alg = x.algebra if isinstance(x, SectorJump) else x.algebra_in
        dense = channel_to_dense(x, layout, layout)
    else:
        alg = x.algebra
        dense = to_dense(x, layout)
    _, _, dims = _layout(alg, layout)
    arr = dense.detach().cpu().numpy()
    return dq.asqarray(arr, dims=tuple(dims[0])) if dims is not None else dq.asqarray(arr)
