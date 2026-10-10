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
import warnings
from typing import List, Optional, Sequence

import torch

from .algebra import TypeIAlgebra
from . import cost
from .states import DensityMatrix, _trace_weights

Operator = TypeIAlgebra.Operator


def _same_structure(a: TypeIAlgebra, b: TypeIAlgebra) -> bool:
    return a is b or (a.k_factors == b.k_factors and a.n_factors == b.n_factors and a.charges == b.charges
                      and a.dtype == b.dtype and a.hilbert.device == b.hilbert.device)


def _check_input(op, alg: TypeIAlgebra, what: str) -> None:
    if not isinstance(op, Operator):
        raise TypeError(f"{what}: expected an Operator or DensityMatrix of {alg}, got {type(op).__name__} "
                        f"(wrap tensors with alg.operator(...))")
    if not _same_structure(op.algebra, alg):
        raise ValueError(f"{what}: the operator belongs to {op.algebra} ({op.algebra.dtype}), the channel to "
                         f"{alg} ({alg.dtype}); convert with op.cast(alg)")


class Channel:
    """Sector-preserving completely positive map given by Kraus operators."""

    def __init__(self, algebra: TypeIAlgebra, kraus: torch.Tensor):
        if kraus.dim() == 4:                   # (r, C, k, k): a single channel
            kraus = kraus.unsqueeze(0)
        expected = (algebra.C, algebra.k_max, algebra.k_max)
        if kraus.dim() != 5 or tuple(kraus.shape[2:]) != expected:
            raise ValueError(f"Kraus tensor must have shape (batch, r, {algebra.C}, {algebra.k_max}, "
                             f"{algebra.k_max}), got {tuple(kraus.shape)}")
        kraus = algebra._to_field(kraus, "Kraus operators")
        self.algebra = algebra
        self.kraus = kraus.to(device=algebra.hilbert.device)
        self._trace_preserving = None

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
        """
        Schroedinger picture: rho -> sum_i K_i rho K_i^* (lazy).
        A DensityMatrix is mapped to a DensityMatrix if the channel is trace preserving;
        otherwise the result is a plain Operator and a warning is issued.
        """
        _check_input(op, self.algebra, "Channel.apply")
        K = self.kraus

        def generator():
            rho = op.matrix.unsqueeze(1)
            dtype = torch.promote_types(rho.dtype, K.dtype)
            B = max(rho.shape[0], K.shape[0])
            cost.check_memory(3 * cost.tensor_bytes((B, *K.shape[1:]), dtype), K.device,
                              f"Channel.apply (batch {B}, Kraus rank {K.shape[1]})")
            Kd = K.to(dtype)
            return (Kd @ rho.to(dtype) @ Kd.conj().transpose(-2, -1)).sum(dim=1)

        if isinstance(op, DensityMatrix):
            if self._trace_preserving is None:
                self._trace_preserving = self.is_trace_preserving()
            if self._trace_preserving:
                return DensityMatrix(self.algebra, generator=generator, validate=False)
            warnings.warn("channel is not trace preserving: the image of a DensityMatrix is "
                          "returned as an unnormalised Operator", stacklevel=2)
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
        if isinstance(other, InterSectorChannel):
            return self.to_inter_sector().compose(other)
        if not isinstance(other, Channel) or not _same_structure(other.algebra, self.algebra):
            raise ValueError("compose: both channels must act on the same algebra")
        K = self.kraus.unsqueeze(2)            # (B, r1, 1, C, k, k)
        L = other.kraus.unsqueeze(1)           # (B, 1, r2, C, k, k)
        KL = K @ L
        B = KL.shape[0]
        return Channel(self.algebra, KL.reshape(B, -1, *KL.shape[3:]))

    __matmul__ = compose

    def to_inter_sector(self) -> 'InterSectorChannel':
        """The same map as an InterSectorChannel with algebra_in = algebra_out."""
        K = self.kraus
        B, r, C, k, _ = K.shape
        G = torch.zeros(B, r, C, C, k, k, dtype=K.dtype, device=K.device)
        idx = torch.arange(C, device=K.device)
        G[:, :, idx, idx] = K
        return InterSectorChannel(self.algebra, self.algebra, G)

    def mix(self, other: 'Channel', p: float) -> 'Channel':
        """Convex combination (1 - p) self + p other."""
        if not _same_structure(other.algebra, self.algebra):
            raise ValueError("mix: both channels must act on the same algebra")
        if isinstance(p, torch.Tensor) and p.numel() != 1:
            raise ValueError("Channel.mix: p must be a single number (mix channels of a batch one by one)")
        p = float(p)
        if not 0 <= p <= 1:
            raise ValueError("mix: p must lie in [0, 1]")
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
        return bool(torch.all(self.trace_preservation_error() <= (self._tol() if tol is None else tol)))

    def is_unital(self, tol: Optional[float] = None) -> bool:
        return bool(torch.all(self.unitality_error() <= (self._tol() if tol is None else tol)))

    # ------------------------------------------------------------------
    def choi(self) -> List[torch.Tensor]:
        """
        Choi matrices per channel c, J_c = sum_ij |i><j| (x) Phi(|i><j|), shape (batch, k_c^2, k_c^2).
        Positive semidefinite for a completely positive map; Tr_out J = 1 iff trace preserving.
        """
        out = []
        cost.check_memory(sum(cost.tensor_bytes((self.batch_size, k ** 2, k ** 2), self.kraus.dtype)
                              for k in self.algebra.k_factors), self.kraus.device, "Channel.choi")
        for c, k_c in enumerate(self.algebra.k_factors):
            K = self.kraus[:, :, c, :k_c, :k_c]                       # (B, r, a, i)
            v = K.transpose(-2, -1).reshape(K.shape[0], K.shape[1], k_c * k_c)  # index (i, a)
            out.append(torch.einsum('brx,bry->bxy', v, v.conj()))
        return out

    def superoperator(self) -> List[torch.Tensor]:
        """Matrix of Phi on row-major vec(rho), per channel: S_c = sum_i K_i (x) conj(K_i), (batch, k_c^2, k_c^2)."""
        out = []
        cost.check_memory(sum(cost.tensor_bytes((self.batch_size, k ** 2, k ** 2), self.kraus.dtype)
                              for k in self.algebra.k_factors), self.kraus.device, "Channel.superoperator")
        for c, k_c in enumerate(self.algebra.k_factors):
            K = self.kraus[:, :, c, :k_c, :k_c]
            S = torch.einsum('nrai,nrbj->nabij', K, K.conj()).reshape(K.shape[0], k_c * k_c, k_c * k_c)
            out.append(S)
        return out

    @staticmethod
    def from_superoperator(algebra: TypeIAlgebra, blocks: Sequence[torch.Tensor], tol: Optional[float] = None) -> 'Channel':
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
            w, V = cost.batched_call(torch.linalg.eigh, J, f"Choi eigendecomposition (sector {c})")
            scale = max(1.0, w.abs().max().item())
            if torch.any(w < -(1e-4 if tol is None else tol) * scale):
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
    rdt = alg.hilbert.real_dtype
    if weights is None:                       # in the algebra's precision, so that the channel is TP to round-off
        weights = torch.distributions.Dirichlet(torch.ones(n_unitaries, dtype=rdt)).sample((batch_size,))
    weights = weights.to(device=alg.hilbert.device, dtype=rdt)
    weights = weights / weights.sum(dim=-1, keepdim=True)
    Us = torch.stack([alg.random_unitary_operator(batch_size).matrix for _ in range(n_unitaries)], dim=1)
    return Channel(alg, Us * weights.sqrt()[:, :, None, None, None].to(Us.dtype))


# ======================================================================
# Maps between sectors / between different algebras
# ======================================================================
class InterSectorChannel:
    """
    Completely positive map Phi: M_in = (+)_c M_{k_c} -> M_out = (+)_d M_{m_d} between (possibly
    different) Type I algebras, in Schroedinger form

        Phi(rho)_d = sum_c sum_i K_i^{dc} rho_c (K_i^{dc})^*,      K_i^{dc}: C^{k_c} -> C^{m_d}.

    Every normal CP map between such algebras has this form: Phi = sum_{d,c} Phi_{dc} with
    Phi_{dc}: M_{k_c} -> M_{m_d} CP. Kraus operators are stored as one tensor of shape
    (batch, r, D, C, m_max, k_max); blocks of rank < r are zero padded. The classical part of
    Phi is the matrix of sector transitions T[d, c] = Tr Phi_{dc}(rho_c) / Tr rho_c.

    Trace preserving iff sum_{d,i} K_i^{dc*} K_i^{dc} = 1_c for every input sector c.
    """

    def __init__(self, algebra_in: TypeIAlgebra, algebra_out: TypeIAlgebra, kraus: torch.Tensor):
        if kraus.dim() == 5:
            kraus = kraus.unsqueeze(0)
        expected = (algebra_out.C, algebra_in.C, algebra_out.k_max, algebra_in.k_max)
        if kraus.dim() != 6 or tuple(kraus.shape[2:]) != expected:
            raise ValueError(f"Kraus tensor must have shape (batch, r, {expected[0]}, {expected[1]}, "
                             f"{expected[2]}, {expected[3]}), got {tuple(kraus.shape)}")
        kraus = algebra_out._to_field(kraus, "Kraus operators")
        self.algebra_in = algebra_in
        self.algebra_out = algebra_out
        self.kraus = kraus.to(device=algebra_out.hilbert.device)
        self._trace_preserving = None

    def __repr__(self) -> str:
        return (f"InterSectorChannel(batch={self.kraus.shape[0]}, kraus_rank={self.kraus.shape[1]}, "
                f"C_in={self.algebra_in.C}, C_out={self.algebra_out.C})")

    # ------------------------------------------------------------------
    @staticmethod
    def from_blocks(algebra_in: TypeIAlgebra, algebra_out: TypeIAlgebra, blocks: dict,
                    batch_size: int = 1) -> 'InterSectorChannel':
        """
        Build from {(d, c): [K_1, K_2, ...]} with K_i of shape (m_d, k_c) or (batch, m_d, k_c);
        missing pairs (d, c) are zero maps.
        """
        r = max(len(v) for v in blocks.values())
        dtype = algebra_out.hilbert.dtype
        K = torch.zeros(batch_size, r, algebra_out.C, algebra_in.C, algebra_out.k_max, algebra_in.k_max,
                        dtype=dtype, device=algebra_out.hilbert.device)
        for (d, c), ops in blocks.items():
            m_d, k_c = algebra_out.k_factors[d], algebra_in.k_factors[c]
            for i, op in enumerate(ops):
                op = torch.as_tensor(op)
                assert op.shape[-2:] == (m_d, k_c), f"block ({d}, {c}) must be {m_d} x {k_c}"
                K[:, i, d, c, :m_d, :k_c] = algebra_out._to_field(op, f"Kraus block ({d}, {c})").to(dtype)
        return InterSectorChannel(algebra_in, algebra_out, K)

    # ------------------------------------------------------------------
    def apply(self, op: Operator) -> Operator:
        """Schroedinger picture (lazy). Trace-preserving maps send a DensityMatrix to a DensityMatrix."""
        _check_input(op, self.algebra_in, "InterSectorChannel.apply")
        K = self.kraus

        def generator():
            rho = op.matrix
            dtype = torch.promote_types(rho.dtype, K.dtype)
            B, r, D, C, m, _ = K.shape
            B = max(B, rho.shape[0])
            cost.check_memory(3 * cost.tensor_bytes((B, r, D, C, m, max(m, rho.shape[-1])), dtype), K.device,
                              f"InterSectorChannel.apply (batch {B}, Kraus rank {r}, {D}x{C} sector pairs)")
            Kd = K.to(dtype)
            rho = rho.to(dtype)[:, None, None]                             # (B, 1, 1, C, k, k)
            out = Kd @ rho @ Kd.conj().transpose(-2, -1)                   # (B, r, D, C, m, m)
            return out.sum(dim=(1, 3))                                      # (B, D, m, m)

        if isinstance(op, DensityMatrix):
            if self._trace_preserving is None:
                self._trace_preserving = self.is_trace_preserving()
            if self._trace_preserving:
                return DensityMatrix(self.algebra_out, generator=generator, validate=False)
            warnings.warn("map is not trace preserving: the image of a DensityMatrix is returned "
                          "as an unnormalised Operator", stacklevel=2)
        out = Operator(self.algebra_out, generator=generator)
        if op._is_self_adjoint:
            out._is_self_adjoint = True
        if op._is_positive:
            out._is_positive = True
        return out

    __call__ = apply

    def adjoint(self, trace: str = 'blunt') -> 'InterSectorChannel':
        """
        Dual map Phi^*: M_out -> M_in with tr_out(Phi(rho) A) = tr_in(rho Phi^*(A)) for the chosen
        trace ('blunt', 'norm' or 'tau_vN' on both algebras). With weights w_c of the trace,
        Phi^*(A)_c = sum_{d,i} (w^out_d / w^in_c) K_i^{dc*} A_d K_i^{dc}.
        """
        Kd = self.kraus.conj().transpose(-2, -1).transpose(2, 3)          # (B, r, C, D, k, m)
        if trace != 'blunt':
            w_in = _trace_weights(self.algebra_in, trace)
            w_out = _trace_weights(self.algebra_out, trace)
            scale = (w_out[None, :] / w_in[:, None]).sqrt()               # (C, D)
            Kd = Kd * scale[None, None, :, :, None, None].to(device=Kd.device, dtype=Kd.dtype)
        return InterSectorChannel(self.algebra_out, self.algebra_in, Kd)

    dual = adjoint

    def compose(self, other: 'InterSectorChannel') -> 'InterSectorChannel':
        """self o other. Kraus operators K_j^{ed} L_i^{dc}, indexed by (j, i, d)."""
        if isinstance(other, Channel):
            other = other.to_inter_sector()
        if not _same_structure(other.algebra_out, self.algebra_in):
            raise ValueError(f"compose: the output algebra of the first map ({other.algebra_out}) is not the input "
                             f"algebra of the second ({self.algebra_in})")
        K = self.kraus.unsqueeze(2).unsqueeze(5)          # (B, r2, 1, E, D, 1, m_e, m_d)
        L = other.kraus.unsqueeze(1).unsqueeze(3)         # (B, 1, r1, 1, D, C, m_d, k_c)
        KL = K @ L                                        # (B, r2, r1, E, D, C, m_e, k_c)
        B, r2, r1, E, D, C = KL.shape[:6]
        KL = KL.permute(0, 1, 2, 4, 3, 5, 6, 7).reshape(B, r2 * r1 * D, E, C, *KL.shape[-2:])
        return InterSectorChannel(other.algebra_in, self.algebra_out, KL)

    __matmul__ = compose

    # ------------------------------------------------------------------
    def _identity_in(self) -> torch.Tensor:
        alg = self.algebra_in
        eye = torch.zeros(alg.C, alg.k_max, alg.k_max, dtype=self.kraus.dtype, device=self.kraus.device)
        for c, k_c in enumerate(alg.k_factors):
            eye[c, :k_c, :k_c] = torch.eye(k_c, dtype=eye.dtype, device=eye.device)
        return eye

    def trace_preservation_error(self) -> torch.Tensor:
        K = self.kraus
        gram = (K.conj().transpose(-2, -1) @ K).sum(dim=(1, 2))           # (B, C, k, k)
        return (gram - self._identity_in()).abs().amax(dim=(-3, -2, -1))

    def is_trace_preserving(self, tol: Optional[float] = None) -> bool:
        tol = (1e-5 if self.kraus.dtype in (torch.float32, torch.complex64) else 1e-10) if tol is None else tol
        return bool(torch.all(self.trace_preservation_error() <= tol))

    def is_unital(self, tol: Optional[float] = None) -> bool:
        """Phi(1) = 1 (Schroedinger picture, blunt trace)."""
        tol = (1e-5 if self.kraus.dtype in (torch.float32, torch.complex64) else 1e-10) if tol is None else tol
        alg = self.algebra_in
        one = alg.operator(self._identity_in().unsqueeze(0).expand(self.kraus.shape[0], -1, -1, -1).clone())
        out = self.apply(one).matrix
        eye = torch.zeros_like(out)
        for d, m_d in enumerate(self.algebra_out.k_factors):
            eye[:, d, :m_d, :m_d] = torch.eye(m_d, dtype=out.dtype, device=out.device)
        return bool((out - eye).abs().max() <= tol)

    def transition_matrix(self, rho: Optional[Operator] = None) -> torch.Tensor:
        """
        Sector transition probabilities T[b, d, c] = Tr Phi_{dc}(rho_c) / Tr rho_c, for the given state
        or, by default, for the maximally mixed state of every input sector. Columns sum to 1 for a
        trace-preserving map.
        """
        alg = self.algebra_in
        if rho is None:
            mat = self._identity_in().unsqueeze(0)
            for c, k_c in enumerate(alg.k_factors):
                if k_c:
                    mat[:, c] = mat[:, c] / k_c
        else:
            mat = rho.matrix
            tr = torch.diagonal(mat, dim1=-2, dim2=-1).sum(-1).real
            mat = mat / torch.where(tr > 0, tr, torch.ones_like(tr))[:, :, None, None].to(mat.dtype)
        K = self.kraus
        dtype = torch.promote_types(mat.dtype, K.dtype)
        out = K.to(dtype) @ mat.to(dtype)[:, None, None] @ K.to(dtype).conj().transpose(-2, -1)
        return torch.diagonal(out, dim1=-2, dim2=-1).sum(-1).real.sum(dim=1)   # (B, D, C)


def random_inter_sector_channel(algebra_in: TypeIAlgebra, algebra_out: TypeIAlgebra, kraus_rank: int = 1,
                                batch_size: int = 1) -> InterSectorChannel:
    """
    Random trace-preserving map from a Haar-random Stinespring isometry per input sector,
    V_c: C^{k_c} -> ((+)_d C^{m_d}) (x) C^r, split into the blocks K_i^{dc}.
    """
    if algebra_in.hilbert.complex_valued != algebra_out.hilbert.complex_valued:
        raise ValueError("both algebras must be real or both complex")
    m = algebra_out.k_factors
    n_out = sum(m)
    K = torch.zeros(batch_size, kraus_rank, algebra_out.C, algebra_in.C, algebra_out.k_max, algebra_in.k_max,
                    dtype=algebra_out.hilbert.dtype, device=algebra_out.hilbert.device)
    offsets = [sum(m[:d]) for d in range(len(m))]
    for c, k_c in enumerate(algebra_in.k_factors):
        if k_c == 0:
            continue
        if n_out * kraus_rank < k_c:
            raise ValueError("output too small for an isometry: increase kraus_rank")
        V = algebra_out.random_unitary(n_out * kraus_rank, batch_size=batch_size)[..., :k_c]
        V = V.reshape(batch_size, kraus_rank, n_out, k_c)
        for d, m_d in enumerate(m):
            K[:, :, d, c, :m_d, :k_c] = V[:, :, offsets[d]:offsets[d] + m_d, :]
    return InterSectorChannel(algebra_in, algebra_out, K)
