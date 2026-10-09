"""
Quantum channels on a Type I algebra M = (+)_c M_{k_c} that preserve the sectors.

A channel is stored by its Kraus operators, one tensor of shape
(batch, r, C, k_max, k_max): r Kraus operators K_i = (+)_c K_{i,c}, block c in the
top-left k_c x k_c corner as for operators. In the Schroedinger picture

    Phi(rho) = sum_i K_i rho K_i^*,

trace preserving iff sum_i K_i^* K_i = 1 and unital iff sum_i K_i K_i^* = 1.
The Heisenberg-picture (dual) channel has Kraus operators K_i^*. Because every block is
mapped to itself, the duality Tr(Phi(rho) A) = Tr(rho Phi^*(A)) holds blockwise and hence for
each of the traces Tr_blunt, Tr_norm and tau_vN.

The Kraus batch dimension is 1 (one channel for all samples) or equal to the batch of the
operators it acts on.
"""
from typing import List, Optional, Sequence

import torch

from .algebra import TypeIAlgebra

Operator = TypeIAlgebra.Operator


class Channel:
    """Sector-preserving completely positive map given by Kraus operators."""

    def __init__(self, algebra: TypeIAlgebra, kraus: torch.Tensor):
        if kraus.dim() == 4:                   # (r, C, k, k): a single channel
            kraus = kraus.unsqueeze(0)
        expected = (algebra.C, algebra.k_max, algebra.k_max)
        if kraus.dim() != 5 or tuple(kraus.shape[2:]) != expected:
            raise ValueError(f"Kraus tensor must have shape (batch, r, {algebra.C}, {algebra.k_max}, "
                             f"{algebra.k_max}), got {tuple(kraus.shape)}")
        self.algebra = algebra
        self.kraus = kraus.to(device=algebra.hilbert.device)

    # ------------------------------------------------------------------
    @property
    def batch_size(self) -> int:
        return self.kraus.shape[0]

    @property
    def kraus_rank(self) -> int:
        return self.kraus.shape[1]

    def __repr__(self) -> str:
        return f"Channel(batch={self.batch_size}, kraus_rank={self.kraus_rank}, C={self.algebra.C}, k_max={self.algebra.k_max})"

    # ------------------------------------------------------------------
    def apply(self, op: Operator) -> Operator:
        """Schroedinger picture: rho -> sum_i K_i rho K_i^*."""
        assert op.algebra is self.algebra
        K = self.kraus

        def generator():
            rho = op.matrix.unsqueeze(1)
            dtype = torch.promote_types(rho.dtype, K.dtype)
            Kd = K.to(dtype)
            return (Kd @ rho.to(dtype) @ Kd.conj().transpose(-2, -1)).sum(dim=1)

        out = Operator(self.algebra, generator=generator)
        if op._is_self_adjoint:
            out._is_self_adjoint = True
        if op._is_positive:
            out._is_positive = True
        return out

    __call__ = apply

    def adjoint(self) -> 'Channel':
        """Heisenberg-picture channel A -> sum_i K_i^* A K_i."""
        return Channel(self.algebra, self.kraus.conj().transpose(-2, -1))

    dual = adjoint

    def compose(self, other: 'Channel') -> 'Channel':
        """self o other: first `other`, then `self` (Kraus operators K_i L_j)."""
        assert other.algebra is self.algebra
        K = self.kraus.unsqueeze(2)            # (B, r1, 1, C, k, k)
        L = other.kraus.unsqueeze(1)           # (B, 1, r2, C, k, k)
        KL = K @ L
        B = KL.shape[0]
        return Channel(self.algebra, KL.reshape(B, -1, *KL.shape[3:]))

    __matmul__ = compose

    def mix(self, other: 'Channel', p: float) -> 'Channel':
        """Convex combination (1 - p) self + p other."""
        assert other.algebra is self.algebra
        B = max(self.batch_size, other.batch_size)
        a = (1 - p) ** 0.5 * self.kraus.expand(B, *self.kraus.shape[1:])
        b = p ** 0.5 * other.kraus.expand(B, *other.kraus.shape[1:]).to(a.dtype)
        return Channel(self.algebra, torch.cat([a, b], dim=1))

    # ------------------------------------------------------------------
    def _active_identity(self) -> torch.Tensor:
        alg = self.algebra
        eye = torch.zeros(alg.C, alg.k_max, alg.k_max, dtype=self.kraus.dtype, device=self.kraus.device)
        for c, k_c in enumerate(alg.k_factors):
            eye[c, :k_c, :k_c] = torch.eye(k_c, dtype=self.kraus.dtype, device=self.kraus.device)
        return eye

    def _deviation(self, gram: torch.Tensor) -> torch.Tensor:
        return (gram - self._active_identity()).abs().amax(dim=(-3, -2, -1))

    def trace_preservation_error(self) -> torch.Tensor:
        """max |sum_i K_i^* K_i - 1| per channel in the batch."""
        K = self.kraus
        return self._deviation((K.conj().transpose(-2, -1) @ K).sum(dim=1))

    def unitality_error(self) -> torch.Tensor:
        """max |sum_i K_i K_i^* - 1| per channel in the batch."""
        K = self.kraus
        return self._deviation((K @ K.conj().transpose(-2, -1)).sum(dim=1))

    def _tol(self) -> float:
        return 1e-5 if self.kraus.dtype in (torch.float32, torch.complex64) else 1e-10

    def is_trace_preserving(self, tol: Optional[float] = None) -> bool:
        return bool(torch.all(self.trace_preservation_error() <= (tol or self._tol())))

    def is_unital(self, tol: Optional[float] = None) -> bool:
        return bool(torch.all(self.unitality_error() <= (tol or self._tol())))

    # ------------------------------------------------------------------
    def choi(self) -> List[torch.Tensor]:
        """
        Choi matrices per channel c, J_c = sum_ij |i><j| (x) Phi(|i><j|), shape (batch, k_c^2, k_c^2).
        Positive semidefinite for a completely positive map; Tr_out J = 1 iff trace preserving.
        """
        out = []
        for c, k_c in enumerate(self.algebra.k_factors):
            K = self.kraus[:, :, c, :k_c, :k_c]                       # (B, r, a, i)
            v = K.transpose(-2, -1).reshape(K.shape[0], K.shape[1], k_c * k_c)  # index (i, a)
            out.append(torch.einsum('brx,bry->bxy', v, v.conj()))
        return out

    def superoperator(self) -> List[torch.Tensor]:
        """Matrix of Phi on row-major vec(rho), per channel: S_c = sum_i K_i (x) conj(K_i), (batch, k_c^2, k_c^2)."""
        out = []
        for c, k_c in enumerate(self.algebra.k_factors):
            K = self.kraus[:, :, c, :k_c, :k_c]
            S = torch.einsum('nrai,nrbj->nabij', K, K.conj()).reshape(K.shape[0], k_c * k_c, k_c * k_c)
            out.append(S)
        return out

    @staticmethod
    def from_superoperator(algebra: TypeIAlgebra, blocks: Sequence[torch.Tensor], tol: float = 0.0) -> 'Channel':
        """
        Channel from per-channel superoperators S_c of shape (batch, k_c^2, k_c^2) acting on row-major
        vec(rho). Kraus operators are obtained from the eigendecomposition of the Choi matrix; S_c must
        describe a completely positive map (eigenvalues of the Choi matrix below -tol raise).
        """
        assert len(blocks) == algebra.C
        k_max, C = algebra.k_max, algebra.C
        batch = max(S.shape[0] if S.dim() == 3 else 1 for S in blocks)
        r = k_max * k_max
        kraus = None
        for c, (S, k_c) in enumerate(zip(blocks, algebra.k_factors)):
            if k_c == 0:
                continue
            S = S if S.dim() == 3 else S.unsqueeze(0)
            S = S.expand(batch, *S.shape[1:])
            # J[(i,a),(j,b)] = S[(a,b),(i,j)]
            J = S.reshape(batch, k_c, k_c, k_c, k_c).permute(0, 3, 1, 4, 2).reshape(batch, k_c * k_c, k_c * k_c)
            J = 0.5 * (J + J.conj().transpose(-2, -1))
            w, V = torch.linalg.eigh(J)
            scale = max(1.0, w.abs().max().item())
            if torch.any(w < -(tol or 1e-4) * scale):
                raise ValueError(f"superoperator of channel {c} is not completely positive "
                                 f"(Choi eigenvalue {w.min().item():.3g})")
            w = torch.clamp(w, min=0.0)
            if kraus is None:
                kraus = torch.zeros(batch, r, C, k_max, k_max, dtype=algebra.hilbert.dtype, device=V.device)
            # K_m[a, i] = sqrt(w_m) V[(i, a), m]
            Km = (V * w.sqrt().unsqueeze(-2).to(V.dtype)).transpose(-2, -1)
            Km = Km.reshape(batch, k_c * k_c, k_c, k_c).transpose(-2, -1)
            if torch.is_complex(Km) and not torch.is_complex(kraus):
                if Km.imag.abs().max() > 1e-5 * max(1.0, Km.abs().max().item()):
                    raise ValueError("complex Kraus operators require a complex-valued algebra")
                Km = Km.real
            kraus[:, :k_c * k_c, c, :k_c, :k_c] = Km.to(kraus.dtype)
        return Channel(algebra, kraus)


# ======================================================================
# Constructors
# ======================================================================
def _empty_kraus(alg: TypeIAlgebra, r: int, batch: int = 1, dtype=None) -> torch.Tensor:
    return torch.zeros(batch, r, alg.C, alg.k_max, alg.k_max,
                       dtype=dtype or alg.hilbert.dtype, device=alg.hilbert.device)


def identity_channel(alg: TypeIAlgebra) -> Channel:
    K = _empty_kraus(alg, 1)
    for c, k_c in enumerate(alg.k_factors):
        K[0, 0, c, :k_c, :k_c] = torch.eye(k_c, dtype=K.dtype, device=K.device)
    return Channel(alg, K)


def unitary_channel(U: Operator) -> Channel:
    """rho -> U rho U^* (U unitary in the algebra, possibly batched)."""
    return Channel(U.algebra, U.matrix.unsqueeze(1))


def kraus_channel(ops: Sequence[Operator]) -> Channel:
    """Channel with the given Kraus operators (all from the same algebra)."""
    alg = ops[0].algebra
    mats = [op.matrix for op in ops]
    B = max(m.shape[0] for m in mats)
    return Channel(alg, torch.stack([m.expand(B, *m.shape[1:]) for m in mats], dim=1))


def dephasing_channel(alg: TypeIAlgebra, p: float) -> Channel:
    """rho_c -> (1 - p) rho_c + p diag(rho_c) in the standard basis of each block (unital)."""
    K = _empty_kraus(alg, 1 + alg.k_max)
    for c, k_c in enumerate(alg.k_factors):
        K[0, 0, c, :k_c, :k_c] = (1 - p) ** 0.5 * torch.eye(k_c, dtype=K.dtype, device=K.device)
        for i in range(k_c):
            K[0, 1 + i, c, i, i] = p ** 0.5
    return Channel(alg, K)


def depolarizing_channel(alg: TypeIAlgebra, p: float) -> Channel:
    """
    rho_c -> (1 - p) rho_c + p Tr(rho_c) 1/k_c in every block (unital, keeps the weight of each
    sector). p = 1 is the trace-preserving conditional expectation onto the centre of M.
    """
    km = alg.k_max
    K = _empty_kraus(alg, 1 + km * km)
    for c, k_c in enumerate(alg.k_factors):
        if k_c == 0:
            continue
        K[0, 0, c, :k_c, :k_c] = (1 - p) ** 0.5 * torch.eye(k_c, dtype=K.dtype, device=K.device)
        for i in range(k_c):
            for j in range(k_c):
                K[0, 1 + i * km + j, c, i, j] = (p / k_c) ** 0.5
    return Channel(alg, K)


def center_expectation(alg: TypeIAlgebra) -> Channel:
    """Conditional expectation onto the centre Z(M): A -> sum_c (Tr A_c / k_c) P_c."""
    return depolarizing_channel(alg, 1.0)


def amplitude_damping_channel(alg: TypeIAlgebra, gamma: float) -> Channel:
    """
    Decay of every level i > 0 to the level 0 of its block with probability gamma:
    K_0 = |0><0| + sqrt(1 - gamma) sum_{i>0} |i><i|,  K_i = sqrt(gamma) |0><i|.
    Trace preserving, not unital for gamma > 0.
    """
    K = _empty_kraus(alg, alg.k_max)
    for c, k_c in enumerate(alg.k_factors):
        if k_c == 0:
            continue
        K[0, 0, c, 0, 0] = 1.0
        for i in range(1, k_c):
            K[0, 0, c, i, i] = (1 - gamma) ** 0.5
            K[0, i, c, 0, i] = gamma ** 0.5
    return Channel(alg, K)


def measurement_channel(projections: Sequence[Operator]) -> Channel:
    """Non-selective Lueders measurement rho -> sum_j P_j rho P_j for a resolution of the identity."""
    return kraus_channel(projections)


def random_channel(alg: TypeIAlgebra, kraus_rank: int, batch_size: int = 1) -> Channel:
    """
    Random trace-preserving channel from a Haar-random Stinespring isometry
    V_c: C^{k_c} -> C^{k_c} (x) C^r in every block, K_{i,c} = (<i| (x) 1) V_c.
    """
    K = _empty_kraus(alg, kraus_rank, batch_size)
    for c, k_c in enumerate(alg.k_factors):
        if k_c == 0:
            continue
        V = alg.random_unitary(k_c * kraus_rank, batch_size=batch_size)[..., :k_c]   # (B, r k, k)
        K[:, :, c, :k_c, :k_c] = V.reshape(batch_size, kraus_rank, k_c, k_c)
    return Channel(alg, K)


def random_mixed_unitary_channel(alg: TypeIAlgebra, n_unitaries: int, batch_size: int = 1,
                                 weights: Optional[torch.Tensor] = None) -> Channel:
    """rho -> sum_j p_j U_j rho U_j^* with Haar unitaries in M and Dirichlet(1) (or given) weights; unital."""
    if weights is None:
        weights = torch.distributions.Dirichlet(torch.ones(n_unitaries)).sample((batch_size,))
    weights = weights.to(alg.hilbert.device)
    Us = torch.stack([alg.random_unitary_operator(batch_size).matrix for _ in range(n_unitaries)], dim=1)
    return Channel(alg, Us * weights.sqrt()[:, :, None, None, None].to(Us.dtype))
