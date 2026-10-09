"""
States (density matrices) on M = (+)_c M_{k_c} and information-theoretic quantities.

A state is a positive operator rho = (+)_c rho_c with Tr_blunt(rho) = sum_c Tr rho_c = 1, so
p_c = Tr rho_c is the classical probability of sector c (the restriction of rho to the
centre of M). All functions are batched over the leading dimension.
"""
from typing import Optional, Sequence, Tuple

import torch

from .algebra import TypeIAlgebra

Operator = TypeIAlgebra.Operator


def _density(alg: TypeIAlgebra, mat: torch.Tensor) -> Operator:
    tr = torch.diagonal(mat, dim1=-2, dim2=-1).sum(dim=(-2, -1)).real
    mat = mat / tr[:, None, None, None].to(mat.dtype)
    mat = 0.5 * (mat + mat.conj().transpose(-2, -1))
    return alg.operator(mat, is_self_adjoint=True, is_positive=True)


def _ginibre(alg: TypeIAlgebra, batch_size: int, cols: int) -> torch.Tensor:
    G = torch.randn(batch_size, alg.C, alg.k_max, cols, dtype=alg.hilbert.dtype, device=alg.hilbert.device)
    for c, k_c in enumerate(alg.k_factors):
        G[:, c, k_c:, :] = 0
    return G


def random_density_matrix(alg: TypeIAlgebra, batch_size: int = 1, rank: Optional[int] = None,
                          measure: str = 'hs') -> Operator:
    """
    Random state on M.

    measure='hs':    rho = G G^* / Tr(G G^*) with G = (+)_c G_c, G_c a k_c x rank Ginibre matrix
                     (rank = k_max by default: Hilbert-Schmidt measure within each block).
    measure='bures': rho ~ (1 + U) G G^* (1 + U)^* with U Haar unitary in M (Bures-type measure).

    The sector weights p_c = Tr rho_c are then random as well (proportional to the squared
    Frobenius norms of the blocks).
    """
    G = _ginibre(alg, batch_size, rank or alg.k_max)            # (B, C, k_max, rank)
    if measure == 'bures':
        U = alg.random_unitary_operator(batch_size).matrix
        eye = alg.identity(batch_size).matrix.to(U.dtype)
        G = (eye + U) @ G
    elif measure != 'hs':
        raise ValueError("measure must be 'hs' or 'bures'")
    return _density(alg, G @ G.conj().transpose(-2, -1))


def maximally_mixed_state(alg: TypeIAlgebra, batch_size: int = 1) -> Operator:
    """rho = 1 / sum_c k_c (the normalised trace of the full matrix algebra on H)."""
    return _density(alg, alg.identity(batch_size).matrix.clone())


def tracial_state(alg: TypeIAlgebra, batch_size: int = 1) -> Operator:
    """Density of tau_vN with respect to Tr_blunt: rho = (+)_c 1_c / (C k_c)."""
    return _density(alg, _tracial_matrix(alg, batch_size))


def _tracial_matrix(alg: TypeIAlgebra, batch_size: int) -> torch.Tensor:
    mat = torch.zeros(batch_size, alg.C, alg.k_max, alg.k_max, dtype=alg.hilbert.dtype, device=alg.hilbert.device)
    for c, k_c in enumerate(alg.k_factors):
        if k_c > 0:
            mat[:, c, :k_c, :k_c] = torch.eye(k_c, dtype=mat.dtype, device=mat.device) / k_c
    return mat


def gibbs_state(H: Operator, beta: float) -> Operator:
    """rho = exp(-beta H) / Z for self-adjoint H (eigenvalues shifted for numerical stability)."""
    w_min = H.lambda_min
    blocks = H.eigh()
    alg = H.algebra
    mat = torch.zeros_like(H.matrix)
    for c, (w, V) in enumerate(blocks):
        k_c = alg.k_factors[c]
        if k_c == 0:
            continue
        e = torch.exp(-beta * (w - w_min.unsqueeze(-1))).to(V.dtype)
        mat[:, c, :k_c, :k_c] = (V * e.unsqueeze(-2)) @ V.conj().transpose(-2, -1)
    return _density(alg, mat)


def partition_function(H: Operator, beta: float) -> torch.Tensor:
    """Z = Tr exp(-beta H) (blunt trace), computed from the spectrum."""
    w = torch.cat([wv[0] for wv in H.eigh() if wv[0].numel()], dim=-1)
    return torch.exp(-beta * w).sum(dim=-1)


def sector_probabilities(rho: Operator) -> torch.Tensor:
    """p_c = Tr rho_c, shape (batch, C): the classical distribution over superselection sectors."""
    return torch.diagonal(rho.matrix, dim1=-2, dim2=-1).sum(dim=-1).real


def expectation(rho: Operator, A: Operator) -> torch.Tensor:
    """<A>_rho = Tr(rho A) (blunt trace)."""
    return (rho @ A).trace


def born_probabilities(rho: Operator, effects: Sequence[Operator]) -> torch.Tensor:
    """Probabilities Tr(rho E_j) of a POVM {E_j}, shape (batch, len(effects))."""
    return torch.stack([expectation(rho, E).real for E in effects], dim=-1)


def lueders_update(rho: Operator, P: Operator, eps: float = 1e-12) -> Tuple[Operator, torch.Tensor]:
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
    return rho.algebra.operator(mat, is_self_adjoint=True, is_positive=True), prob


def von_neumann_entropy(rho: Operator, eps: float = 1e-12) -> torch.Tensor:
    """S(rho) = -Tr rho log rho (natural logarithm)."""
    w = torch.cat([wv[0] for wv in rho.eigh() if wv[0].numel()], dim=-1).clamp(min=0)
    return -(w * torch.log(w.clamp(min=eps))).sum(dim=-1)


def relative_entropy(rho: Operator, sigma: Operator, eps: float = 1e-12) -> torch.Tensor:
    """
    Umegaki relative entropy D(rho || sigma) = Tr rho (log rho - log sigma).
    Eigenvalues of sigma are clamped to eps, so a support mismatch gives a large finite value
    of order -log(eps) instead of +inf.
    """
    log_sigma = sigma.log(eps)
    return -von_neumann_entropy(rho, eps) - (rho @ log_sigma).trace.real


def fidelity(rho: Operator, sigma: Operator) -> torch.Tensor:
    """Uhlmann fidelity F = (Tr |sqrt(rho) sqrt(sigma)|)^2."""
    s = 0.0
    a, b = rho.power(0.5).matrix, sigma.power(0.5).matrix
    for c, k_c in enumerate(rho.algebra.k_factors):
        if k_c > 0:
            s = s + torch.linalg.svdvals(a[:, c, :k_c, :k_c] @ b[:, c, :k_c, :k_c]).sum(dim=-1)
    return s ** 2


def trace_distance(rho: Operator, sigma: Operator) -> torch.Tensor:
    """T = (1/2) ||rho - sigma||_1."""
    return 0.5 * (rho - sigma).trace_norm()


def purity(rho: Operator) -> torch.Tensor:
    """Tr rho^2."""
    return (rho @ rho).trace.real
