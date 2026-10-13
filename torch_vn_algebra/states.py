"""
States (density matrices) on M = (+)_c M_{k_c} and information-theoretic quantities.

A state is a normal positive normalised functional on M. It is represented by its density
rho = (+)_c rho_c with respect to the blunt trace: omega(A) = Tr_blunt(rho A), rho >= 0 and
sum_c Tr rho_c = 1, so p_c = Tr rho_c is the classical probability of sector c (the restriction
of the state to the centre of M). Densities with respect to Tr_norm or tau_vN differ by the
weights k_c and C k_c; see DensityMatrix.density and DensityMatrix.from_density.

All functions are batched over the leading dimension.
"""
from typing import Optional, Sequence, Tuple

import torch

from .algebra import TypeIAlgebra
from . import cost

Operator = TypeIAlgebra.Operator

TRACES = ('blunt', 'norm', 'tau_vN')


def _trace_weights(alg: TypeIAlgebra, trace: str) -> torch.Tensor:
    """w_c such that the given trace is sum_c w_c Tr(A_c), on the algebra's device (cached); empty sectors get
    the weight of a sector of size 1 (they never contribute)."""
    rdt = alg.hilbert.real_dtype
    if trace not in TRACES:
        raise ValueError(f"trace must be one of {TRACES}")
    key = ('trace_weights', trace, rdt)
    w = alg._const_cache.get(key)
    if w is None:
        k = alg._const('k', rdt).clamp(min=1)
        C = sum(1 for k_c in alg.k_factors if k_c)                 # empty sectors are not part of the algebra
        w = torch.ones_like(k) if trace == 'blunt' else 1.0 / k if trace == 'norm' else 1.0 / (C * k)
        alg._const_cache[key] = w
    return w


def trace_of_product(A: Operator, B: Operator) -> torch.Tensor:
    """Tr_blunt(A B) = sum_ij A_ij B_ji in O(k^2) per block, without forming A B."""
    a, b = A.matrix, B.matrix
    dtype = torch.promote_types(a.dtype, b.dtype)
    return (a.to(dtype) * b.to(dtype).transpose(-2, -1)).sum(dim=(-3, -2, -1))


class DensityMatrix(Operator):
    """
    State on M, stored as its density with respect to Tr_blunt (positive, unit blunt trace).

    Generic algebra operations (rho @ A, rho - sigma, ...) return plain Operators; operations
    that map states to states (trace-preserving channels, Hamiltonian and Lindblad evolution,
    mix, condition_on) return DensityMatrix again. With validate=True positivity and the
    normalisation are checked once, when the matrix is materialised.
    """

    def __init__(self, algebra: TypeIAlgebra, generator=None, matrix: Optional[torch.Tensor] = None,
                 validate: bool = True):
        super().__init__(algebra, generator=generator, matrix=matrix,
                         is_self_adjoint=True, is_positive=True)
        self._validate = validate
        if not validate:          # keep the flags, skip the checks
            self._is_self_adjoint_set = False
            self._is_positive_set = False
        elif matrix is not None:
            self._validate_properties()

    def _validate_properties(self):
        if not self._validate:                       # states produced by the library: nothing to check
            return
        super()._validate_properties()
        tr = torch.diagonal(self._matrix, dim1=-2, dim2=-1).sum(dim=(-2, -1)).real
        if bool(((tr - 1).abs() > 100 * self._tol_t()).any()):
            raise ValueError("DensityMatrix must have unit trace (Tr_blunt)")

    def __repr__(self) -> str:
        return "Density" + super().__repr__()

    # ----- conventions -----
    def density(self, trace: str = 'blunt') -> torch.Tensor:
        """Matrix d of the same state with respect to `trace`: omega(A) = trace(d A)."""
        w = _trace_weights(self.algebra, trace)
        return self.matrix / w[None, :, None, None].to(self.matrix.dtype)

    @classmethod
    def from_density(cls, algebra: TypeIAlgebra, d: torch.Tensor, trace: str = 'blunt',
                     normalize: bool = False) -> 'DensityMatrix':
        """State whose density with respect to `trace` is d (shape (batch, C, k_max, k_max))."""
        w = _trace_weights(algebra, trace).to(d.device)
        mat = d * w[None, :, None, None].to(d.dtype)
        if normalize:
            mat = _normalize(mat)
        return cls(algebra, matrix=mat.to(algebra.hilbert.device))

    # ----- functionals -----
    def expectation(self, A: Operator) -> torch.Tensor:
        """<A> = Tr(rho A), real for self-adjoint A."""
        return trace_of_product(self, A)

    def sector_probabilities(self) -> torch.Tensor:
        return sector_probabilities(self)

    def von_neumann_entropy(self, eps: float = 1e-12) -> torch.Tensor:
        return von_neumann_entropy(self, eps)

    def purity(self) -> torch.Tensor:
        return purity(self)

    def fidelity(self, other: 'DensityMatrix') -> torch.Tensor:
        return fidelity(self, other)

    def trace_distance(self, other: 'DensityMatrix') -> torch.Tensor:
        return trace_distance(self, other)

    def relative_entropy(self, other: 'DensityMatrix', eps: float = 1e-12) -> torch.Tensor:
        return relative_entropy(self, other, eps)

    # ----- state -> state -----
    def mix(self, other: 'DensityMatrix', p) -> 'DensityMatrix':
        """(1 - p) rho + p sigma for p in [0, 1]; p is a number or a tensor (batch,) with one weight per sample."""
        from .channels import _same_structure
        if not _same_structure(other.algebra, self.algebra):
            raise ValueError("mix: both states must live in the same algebra")
        if isinstance(p, (int, float)):                  # checked on the host
            if not 0 <= p <= 1:
                raise ValueError("mix: p must lie in [0, 1] (otherwise the result is not a state)")
            pt = p
        else:
            pt = torch.as_tensor(p).to(device=self.algebra.hilbert.device, dtype=self.algebra.hilbert.real_dtype)
            if bool(((pt < 0) | (pt > 1)).any()):        # one synchronisation (validation)
                raise ValueError("mix: p must lie in [0, 1] (otherwise the result is not a state)")
            if pt.dim() == 1:
                pt = pt.reshape(-1, 1, 1, 1)
            elif pt.dim() > 1:
                raise ValueError("mix: p must be a number or a tensor of shape (batch,)")
        return DensityMatrix(self.algebra, generator=lambda: (1 - pt) * self.matrix + pt * other.matrix,
                             validate=False)

    def condition_on(self, P: Operator, eps: float = 1e-12) -> Tuple['DensityMatrix', torch.Tensor]:
        """Lueders rule: (P rho P / Tr(P rho P), Tr(rho P))."""
        return lueders_update(self, P, eps)


def _normalize(mat: torch.Tensor) -> torch.Tensor:
    tr = torch.diagonal(mat, dim1=-2, dim2=-1).sum(dim=(-2, -1)).real
    mat = mat / tr[:, None, None, None].to(mat.dtype)
    return 0.5 * (mat + mat.conj().transpose(-2, -1))


def _density(alg: TypeIAlgebra, mat: torch.Tensor) -> DensityMatrix:
    return DensityMatrix(alg, matrix=_normalize(mat).to(alg.hilbert.device), validate=False)


def _ginibre(alg: TypeIAlgebra, batch_size: int, cols: Optional[int]) -> torch.Tensor:
    """
    Block c is a k_c x n_c Ginibre matrix, zero-padded, with n_c = cols if given; otherwise n_c = k_c for a
    complex and k_c + 1 for a real algebra: then G G^* / Tr has the flat (Hilbert-Schmidt) distribution on the
    states of the block (the Wishart density is proportional to det^((n - k - 1)/2) in the real and det^(n - k) in
    the complex case).
    """
    dev = alg.hilbert.device
    rows = alg._const('rows')                                             # (C, k_max): active rows
    if cols is None:
        extra = 0 if alg.hilbert.complex_valued else 1
        n_cols = alg.k_max + extra
        cmask = torch.arange(n_cols, device=dev) < (alg._const('k', torch.long) + extra)[:, None]   # (C, n_cols)
    else:
        n_cols = cols
        cmask = torch.ones(alg.C, cols, dtype=torch.bool, device=dev)
    G = torch.randn(batch_size, alg.C, alg.k_max, n_cols, dtype=alg.hilbert.dtype, device=dev)
    return G * (rows[:, :, None] & cmask[:, None, :])


def random_density_matrix(alg: TypeIAlgebra, batch_size: int = 1, rank: Optional[int] = None,
                          measure: str = 'hs') -> DensityMatrix:
    """
    Random state on M.

    measure='hs':    rho = G G^* / Tr(G G^*) with G = (+)_c G_c, G_c a k_c x rank Ginibre matrix
                     (by default k_c x k_c, and k_c x (k_c + 1) in a real algebra: the Hilbert-Schmidt
                     measure within each block; an explicit rank r gives the induced measure with an ancilla
                     of dimension r).
    measure='bures': rho ~ (1 + U) G G^* (1 + U)^* with U Haar unitary in M (Bures-type measure).

    The sector weights p_c = Tr rho_c are then random as well (proportional to the squared
    Frobenius norms of the blocks).
    """
    if rank is not None and rank < 1:
        raise ValueError("rank must be a positive integer (or None for the Hilbert-Schmidt measure)")
    r = rank or alg.k_max
    # G, rho and the temporaries of normalisation: about four batches of (k_max x max(k_max, r)) blocks
    cost.check_memory(4 * cost.tensor_bytes((batch_size, alg.C, alg.k_max, max(alg.k_max, r)), alg.hilbert.dtype),
                      alg.hilbert.device, f"random_density_matrix(batch={batch_size}, k_max={alg.k_max})")
    G = _ginibre(alg, batch_size, rank)                         # (B, C, k_max, rank or k_max)
    if measure == 'bures':
        U = alg.random_unitary_operator(batch_size).matrix
        G = G + U @ G                                            # (1 + U) G
    elif measure != 'hs':
        raise ValueError("measure must be 'hs' or 'bures'")
    return _density(alg, G @ G.conj().transpose(-2, -1))


def maximally_mixed_state(alg: TypeIAlgebra, batch_size: int = 1) -> DensityMatrix:
    """rho = 1 / sum_c k_c (the normalised trace of the full matrix algebra on H)."""
    return _density(alg, alg.identity(batch_size).matrix.clone())


def tracial_state(alg: TypeIAlgebra, batch_size: int = 1) -> DensityMatrix:
    """Density of tau_vN with respect to Tr_blunt: rho = (+)_c 1_c / (C k_c)."""
    return _density(alg, _tracial_matrix(alg, batch_size))


def _tracial_matrix(alg: TypeIAlgebra, batch_size: int) -> torch.Tensor:
    w = alg._const('w_norm', alg.hilbert.real_dtype).to(alg.dtype)
    return (alg._const('eye') * w[:, None, None]).expand(batch_size, -1, -1, -1).clone()


def gibbs_state(H: Operator, beta: float) -> DensityMatrix:
    """rho = exp(-beta H) / Z for self-adjoint H (eigenvalues shifted for numerical stability)."""
    if not H.is_self_adjoint:
        raise RuntimeError("gibbs_state requires a self-adjoint Hamiltonian")
    alg = H.algebra
    B = H.matrix.shape[0]
    dec = H._decompositions('hermitian')                     # one eigendecomposition per block size
    w_min = torch.stack([w.reshape(B, -1).amin(-1) for _, _, _, (w, _) in dec], dim=-1).amin(-1)     # (B,)
    mat = torch.zeros_like(H.matrix)
    for k, idx, _, (w, V) in dec:
        shift = w_min.repeat_interleave(len(idx))[:, None]
        e = torch.exp(-beta * (w - shift)).to(V.dtype)
        mat[:, alg._sel(idx), :k, :k] = ((V * e.unsqueeze(-2)) @ V.conj().transpose(-2, -1)).reshape(B, len(idx), k, k)
    return _density(alg, mat)


def partition_function(H: Operator, beta: float) -> torch.Tensor:
    """Z = Tr exp(-beta H) (blunt trace), computed from the spectrum (as exp(log Z), without overflow of the
    terms; Z itself may still overflow, see log_partition_function)."""
    return torch.exp(log_partition_function(H, beta))


def log_partition_function(H: Operator, beta: float) -> torch.Tensor:
    """log Z = log Tr exp(-beta H), by logsumexp over the spectrum."""
    w = torch.cat([w for w in H.eigenvalues() if w.numel()], dim=-1)
    return torch.logsumexp(-beta * w, dim=-1)


def sector_probabilities(rho: Operator) -> torch.Tensor:
    """p_c = Tr rho_c, shape (batch, C): the classical distribution over superselection sectors."""
    return torch.diagonal(rho.matrix, dim1=-2, dim2=-1).sum(dim=-1).real


def expectation(rho: Operator, A: Operator) -> torch.Tensor:
    """<A>_rho = Tr(rho A) (blunt trace)."""
    return trace_of_product(rho, A)


def born_probabilities(rho: Operator, effects: Sequence[Operator]) -> torch.Tensor:
    """Probabilities Tr(rho E_j) of a POVM {E_j}, shape (batch, len(effects))."""
    E = torch.stack(torch.broadcast_tensors(*[e.matrix for e in effects]), dim=1)       # (B, J, C, k, k)
    r = rho.matrix
    dtype = torch.promote_types(r.dtype, E.dtype)
    return (r.to(dtype).unsqueeze(1) * E.to(dtype).transpose(-2, -1)).sum(dim=(-3, -2, -1)).real


def lueders_update(rho: Operator, P: Operator, eps: float = 1e-12) -> Tuple[DensityMatrix, torch.Tensor]:
    """
    Conditional state after observing the outcome of projection P (Lueders rule):
    rho | P = P rho P / Tr(P rho P), together with the probability Tr(rho P).
    Samples with probability below eps are returned as zero operators.
    """
    post = (P @ rho @ P).matrix
    prob = torch.diagonal(post, dim1=-2, dim2=-1).sum(dim=(-2, -1)).real
    safe = torch.where(prob > eps, prob, torch.ones_like(prob))
    mat = torch.where(prob[:, None, None, None] > eps, post / safe[:, None, None, None].to(post.dtype),
                      torch.zeros_like(post))
    return DensityMatrix(rho.algebra, matrix=mat, validate=False), prob


def von_neumann_entropy(rho: Operator, eps: float = 1e-12) -> torch.Tensor:
    """S(rho) = -Tr rho log rho (natural logarithm)."""
    w = torch.cat([w for w in rho.eigenvalues() if w.numel()], dim=-1).clamp(min=0)
    return -(w * torch.log(w.clamp(min=eps))).sum(dim=-1)


def relative_entropy(rho: Operator, sigma: Operator, eps: float = 1e-12) -> torch.Tensor:
    """
    Umegaki relative entropy D(rho || sigma) = Tr rho (log rho - log sigma).
    Eigenvalues of sigma are clamped to eps, so a support mismatch gives a large finite value
    of order -log(eps) instead of +inf.
    """
    log_sigma = sigma.log(eps)
    return -von_neumann_entropy(rho, eps) - trace_of_product(rho, log_sigma).real


def fidelity(rho: Operator, sigma: Operator) -> torch.Tensor:
    """Uhlmann fidelity F = (Tr |sqrt(rho) sqrt(sigma)|)^2 (singular values of sqrt(rho) sqrt(sigma): accurate
    also for nearly pure states, unlike the eigenvalues of sqrt(rho) sigma sqrt(rho))."""
    M = rho.power(0.5) @ sigma.power(0.5)
    M._is_self_adjoint = False                     # not self-adjoint in general: singular values, no check
    return M.trace_norm() ** 2                     # batched over sectors of equal size


def trace_distance(rho: Operator, sigma: Operator) -> torch.Tensor:
    """T = (1/2) ||rho - sigma||_1."""
    return 0.5 * (rho - sigma).trace_norm()


def purity(rho: Operator) -> torch.Tensor:
    """Tr rho^2."""
    return trace_of_product(rho, rho).real
