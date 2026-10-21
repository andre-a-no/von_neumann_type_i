"""
Measurements on M = (+)_c M_{k_c}: POVMs, quantum instruments and the classical Fisher information.

A POVM is a family of effects E_j in M with E_j >= 0 and sum_j E_j = 1; the outcome j of a state rho has the Born
probability p_j = Tr_blunt(rho E_j). The effects are stored as one tensor of shape (batch, J, C, k_max, k_max) (batch
1: one POVM for all samples), zero padded as operators. An instrument {Phi_j} refines a POVM by the state after the
measurement: Phi_j(rho) = sum_i K_{ji} rho K_{ji}^* with Kraus operators in M (sector preserving), sum_j Phi_j trace
preserving, E_j = sum_i K_{ji}^* K_{ji}. As a single channel it is the map from M into (+)_j M (J copies of the sectors
of M, charges (j, c)) that writes the outcome into the centre of the output algebra: the classical register is the
sector label j, its distribution the sector probabilities of the output state.

    povm = POVM.basis(chain.algebra)                     # computational basis of a SpinChain, outcomes (c, a)
    p = povm.probabilities(rho)                          # (batch, J): Tr(rho E_j)
    counts = povm.sample(rho, shots=1000, generator=g)   # (batch, J) integer counts
    nll = -povm.log_likelihood(sigma, counts)            # differentiable in sigma (maximum-likelihood tomography)
    shadow = POVM.combine([POVM.from_unitary(U) for U in Us])    # randomised measurement, outcomes (s, a)
    inst = Instrument.from_povm(povm)                    # Kraus operators sqrt(E_j); Instrument.lueders for projections
    post, p = inst.post_states(rho)                      # conditional states (one DensityMatrix per outcome)
    out = inst.as_channel()(rho)                         # state on (+)_j M, the outcome in the centre
    F = fisher_information(lambda t: povm.probabilities(state(t)), [theta])   # (batch, P, P)

Everything is differentiable in the state, the effects and the Kraus operators. The effects are dense: a POVM costs
J * C * k_max^2 entries per sample, so POVM.basis is meant for algebras of moderate size.
"""
from typing import Callable, List, Optional, Sequence, Tuple

import torch

from .algebra import TypeIAlgebra, _as_tensor
from . import cost
from .channels import Channel, InterSectorChannel, _check_input, _same_structure
from .states import DensityMatrix

Operator = TypeIAlgebra.Operator


def _default_tol(dtype: torch.dtype) -> float:
    return 1e-5 if dtype in (torch.float32, torch.complex64) else 1e-10


def _pad_mask(alg: TypeIAlgebra, t: torch.Tensor) -> torch.Tensor:
    """Zero padding of a tensor (..., C, k_max, k_max) of block-diagonal matrices."""
    if any(k_c < alg.k_max for k_c in alg.k_factors):
        t = t * alg._const('block_mask', torch.bool)
    return t


def _stack_operators(ops: Sequence[Operator], what: str) -> Tuple[TypeIAlgebra, torch.Tensor]:
    """(algebra, tensor (batch, n, C, k_max, k_max)) of operators of one algebra, broadcast over the batch."""
    if not ops:
        raise ValueError(f"{what}: need at least one operator")
    alg = ops[0].algebra
    for op in ops:
        _check_input(op, alg, what)
    mats = [op.matrix for op in ops]
    B = max(m.shape[0] for m in mats)
    dtype = mats[0].dtype
    for m in mats:
        dtype = torch.promote_types(dtype, m.dtype)
    return alg, torch.stack([m.to(dtype).expand(B, *m.shape[1:]) for m in mats], dim=1)


def _born(rho: torch.Tensor, E: torch.Tensor) -> torch.Tensor:
    """Tr(rho E_j) for rho (B, C, k, k) and effects (B', J, C, k, k): (max(B, B'), J), real part."""
    dtype = torch.promote_types(rho.dtype, E.dtype)
    return (rho.to(dtype).unsqueeze(1) * E.to(dtype).transpose(-2, -1)).sum(dim=(-3, -2, -1)).real


def _check_batch(b1: int, b2: int, what: str) -> None:
    if b1 != b2 and 1 not in (b1, b2):
        raise ValueError(f"{what}: batch sizes {b1} and {b2} do not broadcast")


# ======================================================================
# POVM
# ======================================================================
class POVM:
    """
    Positive operator-valued measure {E_j}, j = 0..J-1, in a Type I algebra.

    effects: a sequence of Operators of one algebra, or a tensor (J, C, k_max, k_max) / (batch, J, C, k_max, k_max)
    together with `algebra`. check=True verifies once (one synchronisation) that every E_j is self-adjoint and positive
    and that sum_j E_j = 1, to `tol` (default 1e-10 in double, 1e-5 in single precision). `outcomes` labels the
    outcomes (default 0..J-1; the constructors use (c, a), the sector charges or (setting, outcome)).
    """

    def __init__(self, effects, algebra: Optional[TypeIAlgebra] = None, check: bool = True,
                 outcomes: Optional[Sequence] = None, tol: Optional[float] = None):
        if isinstance(effects, torch.Tensor):
            if algebra is None:
                raise ValueError("POVM: a tensor of effects needs the algebra (POVM(E, algebra=alg))")
            E = effects.unsqueeze(0) if effects.dim() == 4 else effects
            expected = (algebra.C, algebra.k_max, algebra.k_max)
            if E.dim() != 5 or tuple(E.shape[2:]) != expected:
                raise ValueError(f"effects must have shape (batch, J, {expected[0]}, {expected[1]}, {expected[2]}), "
                                 f"got {tuple(effects.shape)}")
            E = _pad_mask(algebra, algebra._to_field(E, "POVM effects"))
        else:
            alg, E = _stack_operators(list(effects), "POVM")
            if algebra is not None and not _same_structure(alg, algebra):
                raise ValueError(f"POVM: the effects belong to {alg}, not to {algebra}")
            algebra = alg
        self.algebra = algebra
        self.effects = E
        self.outcomes = list(range(E.shape[1])) if outcomes is None else list(outcomes)
        if len(self.outcomes) != E.shape[1]:
            raise ValueError(f"POVM: {len(self.outcomes)} outcome labels for {E.shape[1]} effects")
        if check:
            self._check(tol)

    def _check(self, tol: Optional[float]) -> None:
        E = self.effects.detach()
        tol = _default_tol(E.dtype) if tol is None else tol
        Eh = E.conj().transpose(-2, -1)
        skew = cost.abs_max(E - Eh)
        w = torch.linalg.eigvalsh(0.5 * (E + Eh))                          # padding: eigenvalue 0
        total = (E.sum(dim=1) - self.algebra._const('eye', E.dtype)).abs().amax()
        bad = torch.stack([skew > tol, w.amin() < -tol, total > tol]).tolist()   # one synchronisation
        if bad[0]:
            raise ValueError("POVM: the effects must be self-adjoint")
        if bad[1]:
            raise ValueError("POVM: the effects must be positive")
        if bad[2]:
            raise ValueError("POVM: the effects must sum to the identity")

    # ------------------------------------------------------------------
    def __len__(self) -> int:
        return self.effects.shape[1]

    @property
    def batch_size(self) -> int:
        return self.effects.shape[0]

    def __getitem__(self, j: int) -> Operator:
        """The effect E_j as an Operator."""
        return self.algebra.operator(self.effects[:, j], is_self_adjoint=True, is_positive=True)

    def __repr__(self) -> str:
        return f"POVM(batch={self.batch_size}, outcomes={len(self)}, C={self.algebra.C}, k_max={self.algebra.k_max})"

    # ------------------------------------------------------------------
    def probabilities(self, rho: Operator) -> torch.Tensor:
        """Born probabilities p_j = Tr(rho E_j), shape (batch, J) (as states.born_probabilities)."""
        _check_input(rho, self.algebra, "POVM.probabilities")
        _check_batch(rho.matrix.shape[0], self.batch_size, "POVM.probabilities")
        return _born(rho.matrix, self.effects)

    def sample(self, rho: Operator, shots: int, generator: Optional[torch.Generator] = None) -> torch.Tensor:
        """
        Outcome counts of `shots` independent measurements of every sample, shape (batch, J), int64 (multinomial;
        reproducible with a torch.Generator on the algebra's device). Rounding-level negative probabilities count as 0.
        """
        if int(shots) != shots or shots < 0:
            raise ValueError("sample: shots must be a non-negative integer")
        shots = int(shots)
        with torch.no_grad():
            p = self.probabilities(rho).clamp(min=0)
        B, J = p.shape
        counts = torch.zeros(B, J, dtype=torch.long, device=p.device)
        chunk = max(1, 2 ** 22 // max(B, 1))          # outcome indices of at most ~4M draws in memory at a time
        done = 0
        while done < shots:
            n = min(chunk, shots - done)
            idx = torch.multinomial(p, n, replacement=True, generator=generator)        # (B, n)
            counts.scatter_add_(1, idx, torch.ones_like(idx))
            done += n
        return counts

    def log_likelihood(self, rho: Operator, counts, eps: float = 1e-12) -> torch.Tensor:
        """sum_j n_j log max(p_j, eps) for counts n_j of shape (J,) or (batch, J): shape (batch,) (the log-likelihood
        up to the multinomial coefficient, which does not depend on the state)."""
        p = self.probabilities(rho)
        n = _as_tensor(counts).to(device=p.device, dtype=p.dtype)
        if n.shape[-1] != p.shape[-1]:
            raise ValueError(f"log_likelihood: counts must have {p.shape[-1]} outcomes, got shape {tuple(n.shape)}")
        return (n * torch.log(p.clamp(min=eps))).sum(dim=-1)

    def is_projective(self, tol: Optional[float] = None) -> bool:
        """E_j^2 = E_j for every outcome."""
        E = self.effects
        tol = _default_tol(E.dtype) if tol is None else tol
        return bool(cost.abs_max(E @ E - E) <= tol)

    # ------------------------------------------------------------------
    @staticmethod
    def sectors(algebra: TypeIAlgebra) -> 'POVM':
        """The central projections P_c (measurement of the charge): p_c = Tr rho_c, outcomes = algebra.charges."""
        eye = algebra._const('eye')
        idx = torch.arange(algebra.C, device=eye.device)
        E = torch.zeros(1, algebra.C, *eye.shape, dtype=eye.dtype, device=eye.device)
        E[0, idx, idx] = eye
        return POVM(E, algebra, check=False, outcomes=algebra.charges)

    @staticmethod
    def basis(algebra: TypeIAlgebra) -> 'POVM':
        """Rank-one projections onto the basis vectors of all sectors, sector after sector (outcomes (c, a)); for a
        SpinChain the computational basis, outcome (N, a) being the configuration chain.labels[N][a]."""
        pairs = [(c, a) for c, k_c in enumerate(algebra.k_factors) for a in range(k_c)]
        J = len(pairs)
        cost.check_memory(cost.tensor_bytes((J, algebra.C, algebra.k_max, algebra.k_max), algebra.dtype),
                          algebra.hilbert.device, f"POVM.basis ({J} effects)")
        dev = algebra.hilbert.device
        c = torch.tensor([p[0] for p in pairs], device=dev)
        a = torch.tensor([p[1] for p in pairs], device=dev)
        E = torch.zeros(1, J, algebra.C, algebra.k_max, algebra.k_max, dtype=algebra.dtype, device=dev)
        E[0, torch.arange(J, device=dev), c, a, a] = 1
        return POVM(E, algebra, check=False, outcomes=pairs)

    @staticmethod
    def from_unitary(U: Operator, base: Optional['POVM'] = None) -> 'POVM':
        """The POVM {U E_j U^*} of `base` (default: POVM.basis) rotated by a unitary U of the algebra (possibly a
        batch): measuring rho in the rotated basis is measuring U^* rho U in the base, e.g. a tomography setting."""
        alg = U.algebra
        base = POVM.basis(alg) if base is None else base
        _check_input(U, base.algebra, "POVM.from_unitary")
        u = U.matrix.unsqueeze(1)                                          # (B, 1, C, k, k)
        E = u @ base.effects.to(u.dtype) @ u.conj().transpose(-2, -1)
        return POVM(E, base.algebra, check=False, outcomes=base.outcomes)

    @staticmethod
    def combine(povms: Sequence['POVM'], weights=None) -> 'POVM':
        """
        Randomised measurement: setting s is chosen with probability w_s (default uniform; non-negative weights are
        normalised) and then measured. Effects w_s E^s_a, outcomes (s, outcome of setting s).
        """
        povms = list(povms)
        if not povms:
            raise ValueError("combine: need at least one POVM")
        alg = povms[0].algebra
        for P in povms:
            if not _same_structure(P.algebra, alg):
                raise ValueError("combine: all POVMs must belong to the same algebra")
        S = len(povms)
        w = torch.full((S,), 1.0 / S, dtype=alg.hilbert.real_dtype) if weights is None else _as_tensor(weights)
        w = w.to(device=alg.hilbert.device, dtype=alg.hilbert.real_dtype)
        if w.shape != (S,):
            raise ValueError(f"combine: need one weight per POVM ({S}), got shape {tuple(w.shape)}")
        if bool((w.detach() < 0).any()) or not bool(w.detach().sum() > 0):
            raise ValueError("combine: the weights must be non-negative with a positive sum")
        w = w / w.sum()
        B = max(P.batch_size for P in povms)
        E = torch.cat([w[s].to(P.effects.dtype) * P.effects.expand(B, *P.effects.shape[1:])
                       for s, P in enumerate(povms)], dim=1)
        return POVM(E, alg, check=False, outcomes=[(s, o) for s, P in enumerate(povms) for o in P.outcomes])


# ======================================================================
# PVM
# ======================================================================
class PVM(POVM):
    """
    Projection-valued measure: a POVM whose effects E_j are projections (hence mutually orthogonal, as they sum to 1),
    with a real value lambda_j per outcome: the spectral measure of the self-adjoint A = sum_j lambda_j E_j
    (projections.spectral_measure builds it from A). values: shape (J,) or (batch, J), default 0..J-1.

        E = PVM.from_povm(POVM.basis(alg), values)        # a measurement in the computational basis with values
        E.observable(), E.expectation(rho), E.variance(rho)
        E.projection(lambda v: v > 0)                     # spectral projection E(Delta) for Delta = (0, oo)
        E.pushforward(torch.abs)                          # spectral measure of f(A), equal values merged
    """

    def __init__(self, projections, values=None, algebra: Optional[TypeIAlgebra] = None, check: bool = True,
                 outcomes: Optional[Sequence] = None, tol: Optional[float] = None):
        super().__init__(projections, algebra, check=check, outcomes=outcomes, tol=tol)
        J = len(self)
        real = self.algebra.hilbert.real_dtype
        v = torch.arange(J, dtype=real) if values is None else _as_tensor(values)
        if torch.is_complex(v):
            raise ValueError("PVM: the values must be real")
        v = v.to(device=self.effects.device, dtype=real)
        v = v.unsqueeze(0) if v.dim() == 1 else v
        if v.dim() != 2 or v.shape[-1] != J:
            raise ValueError(f"PVM: need one value per outcome ({J}), got shape {tuple(v.shape)}")
        _check_batch(v.shape[0], self.batch_size, "PVM")
        self.values = v
        if check and not self.is_projective(tol):
            raise ValueError("PVM: the effects must be projections")

    @staticmethod
    def from_povm(povm: POVM, values=None, tol: Optional[float] = None) -> 'PVM':
        """A projective POVM with values attached (checked to be projective)."""
        return PVM(povm.effects, values, povm.algebra, check=True, outcomes=povm.outcomes, tol=tol)

    def __repr__(self) -> str:
        return f"PVM(batch={self.batch_size}, outcomes={len(self)}, C={self.algebra.C}, k_max={self.algebra.k_max})"

    def observable(self) -> Operator:
        """A = sum_j lambda_j E_j, self-adjoint."""
        v = self.values.to(self.effects.dtype)[..., None, None, None]
        return self.algebra.operator((v * self.effects).sum(dim=1), is_self_adjoint=True)

    def expectation(self, rho: Operator) -> torch.Tensor:
        """sum_j lambda_j Tr(rho E_j), shape (batch,)."""
        return (self.probabilities(rho) * self.values).sum(-1)

    def variance(self, rho: Operator) -> torch.Tensor:
        """sum_j p_j (lambda_j - mean)^2 for the normalised Born distribution, shape (batch,)."""
        p = self.probabilities(rho)
        p = p / p.sum(-1, keepdim=True)
        m = (p * self.values).sum(-1, keepdim=True)
        return (p * (self.values - m) ** 2).sum(-1)

    def projection(self, predicate: Callable[[torch.Tensor], torch.Tensor]) -> Operator:
        """The spectral projection E(Delta) = sum of E_j over the outcomes with predicate(lambda_j) true (predicate
        maps the values tensor (batch, J) to a boolean tensor of the same shape)."""
        mask = predicate(self.values)
        if mask.shape != self.values.shape or mask.dtype != torch.bool:
            raise ValueError("PVM.projection: the predicate must return a boolean tensor shaped like the values")
        m = (mask.to(self.effects.dtype)[..., None, None, None] * self.effects).sum(dim=1)
        return self.algebra.operator(m, is_self_adjoint=True, is_positive=True, is_projection=True)

    def pushforward(self, f: Callable[[torch.Tensor], torch.Tensor], tol: Optional[float] = None) -> 'PVM':
        """The spectral measure of f(A): values f(lambda_j); for batch 1, outcomes whose values lie within tol of the
        smallest of their group are merged (their projections added); default tol sqrt(eps) * max |f|."""
        v = f(self.values)
        if self.values.shape[0] != 1:
            return PVM(self.effects, v, self.algebra, check=False, outcomes=self.outcomes)
        vals = v[0].tolist()
        if tol is None:                                  # as projections.spectral_measure: relative to max |value|
            tol = _default_tol(self.effects.dtype) ** 0.5 * (max(abs(x) for x in vals) or 1.0)
        order = torch.argsort(v[0]).tolist()
        groups = []
        for j in order:                                  # compare with the first member: no chaining of small gaps
            if groups and vals[j] - vals[groups[-1][0]] <= tol:
                groups[-1].append(j)
            else:
                groups.append([j])
        E = torch.stack([self.effects[:, g].sum(dim=1) for g in groups], dim=1)
        new = torch.stack([v[0, g].mean() for g in groups]).unsqueeze(0)
        return PVM(E, new, self.algebra, check=False)


# ======================================================================
# Instrument
# ======================================================================
class Instrument:
    """
    Quantum instrument {Phi_j}: Phi_j(rho) = sum_i K_{ji} rho K_{ji}^*, sector preserving, sum_j Phi_j trace preserving.

    kraus: a sequence with one entry per outcome, each a Channel (its Kraus operators), an Operator (one Kraus
    operator) or a sequence of Operators; or a tensor (J, r, C, k_max, k_max) / (batch, J, r, C, k_max, k_max) together
    with `algebra` (outcomes with fewer Kraus operators are padded with zeros). check=True verifies
    sum_{j,i} K_{ji}^* K_{ji} = 1 once.
    """

    def __init__(self, kraus, algebra: Optional[TypeIAlgebra] = None, check: bool = True,
                 outcomes: Optional[Sequence] = None, tol: Optional[float] = None):
        if isinstance(kraus, torch.Tensor):
            if algebra is None:
                raise ValueError("Instrument: a tensor of Kraus operators needs the algebra")
            K = kraus.unsqueeze(0) if kraus.dim() == 5 else kraus
            expected = (algebra.C, algebra.k_max, algebra.k_max)
            if K.dim() != 6 or tuple(K.shape[3:]) != expected:
                raise ValueError(f"Kraus tensor must have shape (batch, J, r, {expected[0]}, {expected[1]}, "
                                 f"{expected[2]}), got {tuple(kraus.shape)}")
            K = _pad_mask(algebra, algebra._to_field(K, "Kraus operators"))
        else:
            algebra, K = self._from_sequence(list(kraus), algebra)
        self.algebra = algebra
        self.kraus = K
        self.outcomes = list(range(K.shape[1])) if outcomes is None else list(outcomes)
        if len(self.outcomes) != K.shape[1]:
            raise ValueError(f"Instrument: {len(self.outcomes)} outcome labels for {K.shape[1]} outcomes")
        self._trace_preserving = None
        self._register = None
        if check:
            tol = _default_tol(K.dtype) if tol is None else tol
            err = (self.average_channel().trace_preservation_error()).amax()
            if bool(err > tol):
                raise ValueError(f"Instrument: sum_j Phi_j is not trace preserving (deviation {float(err):.1e})")
            self._trace_preserving = True

    @staticmethod
    def _from_sequence(items: list, algebra: Optional[TypeIAlgebra]) -> Tuple[TypeIAlgebra, torch.Tensor]:
        if not items:
            raise ValueError("Instrument: need at least one outcome")
        per = []
        for item in items:
            if isinstance(item, Channel):
                per.append((item.algebra, item.kraus))
            elif isinstance(item, Operator):
                per.append(_stack_operators([item], "Instrument"))
            else:
                per.append(_stack_operators(list(item), "Instrument"))
        alg = algebra or per[0][0]
        for a, _ in per:
            if not _same_structure(a, alg):
                raise ValueError("Instrument: all Kraus operators must belong to the same algebra")
        B = max(t.shape[0] for _, t in per)
        r = max(t.shape[1] for _, t in per)
        dtype = per[0][1].dtype
        for _, t in per:
            dtype = torch.promote_types(dtype, t.dtype)
        pad = lambda t: torch.cat([t, t.new_zeros(t.shape[0], r - t.shape[1], *t.shape[2:])], dim=1) \
            if t.shape[1] < r else t
        K = torch.stack([pad(t.to(dtype)).expand(B, -1, -1, -1, -1) for _, t in per], dim=1)
        return alg, K

    # ------------------------------------------------------------------
    def __len__(self) -> int:
        return self.kraus.shape[1]

    @property
    def batch_size(self) -> int:
        return self.kraus.shape[0]

    def __repr__(self) -> str:
        return (f"Instrument(batch={self.batch_size}, outcomes={len(self)}, kraus_rank={self.kraus.shape[2]}, "
                f"C={self.algebra.C}, k_max={self.algebra.k_max})")

    def povm(self) -> POVM:
        """The POVM of the instrument, E_j = sum_i K_{ji}^* K_{ji}."""
        K = self.kraus
        E = (K.conj().transpose(-2, -1) @ K).sum(dim=2)
        return POVM(E, self.algebra, check=False, outcomes=self.outcomes)

    def probabilities(self, rho: Operator) -> torch.Tensor:
        """p_j = Tr Phi_j(rho), shape (batch, J)."""
        return self.povm().probabilities(rho)

    def _branches(self, rho: Operator) -> torch.Tensor:
        """Phi_j(rho) for all outcomes, shape (batch, J, C, k_max, k_max)."""
        _check_input(rho, self.algebra, "Instrument")
        r, K = rho.matrix, self.kraus
        _check_batch(r.shape[0], K.shape[0], "Instrument")
        dtype = torch.promote_types(r.dtype, K.dtype)
        B = max(r.shape[0], K.shape[0])
        cost.check_memory(2 * cost.tensor_bytes((B, *K.shape[1:]), dtype), K.device,
                          f"Instrument (batch {B}, {K.shape[1]} outcomes, Kraus rank {K.shape[2]})")
        Kd = K.to(dtype)
        X = Kd @ r.to(dtype)[:, None, None] @ Kd.conj().transpose(-2, -1)   # (B, J, r, C, k, k)
        return X.sum(dim=2)

    def post_states(self, rho: Operator, eps: float = 1e-12) -> Tuple[List[DensityMatrix], torch.Tensor]:
        """
        Conditional states Phi_j(rho) / p_j (one DensityMatrix of the batch per outcome) and the probabilities p_j,
        shape (batch, J). Samples with p_j <= eps get the zero operator, as in states.lueders_update.
        """
        X = self._branches(rho)
        p = torch.diagonal(X, dim1=-2, dim2=-1).sum(dim=(-2, -1)).real                # (B, J)
        ok = (p > eps)[:, :, None, None, None]
        safe = torch.where(p > eps, p, torch.ones_like(p))
        mat = torch.where(ok, X / safe[:, :, None, None, None].to(X.dtype), torch.zeros_like(X))
        return [DensityMatrix(self.algebra, matrix=mat[:, j], validate=False) for j in range(len(self))], p

    def average_channel(self) -> Channel:
        """The non-selective measurement sum_j Phi_j as a Channel (Kraus operators K_{ji} for all j, i)."""
        K = self.kraus
        return Channel(self.algebra, K.reshape(K.shape[0], -1, *K.shape[3:]), trace_preserving=self._trace_preserving)

    @property
    def register_algebra(self) -> TypeIAlgebra:
        """(+)_j M: J copies of the sectors of M, sector (j, c) at index j * C + c with charge (j, charge_c)."""
        if self._register is None:
            alg, J = self.algebra, len(self)
            self._register = TypeIAlgebra(alg.n_factors * J, alg.k_factors * J, complex_valued=alg.complex_valued,
                                          device=alg.hilbert.device, precision=alg.precision,
                                          charges=[(j, q) for j in range(J) for q in alg.charges])
        return self._register

    def as_channel(self) -> InterSectorChannel:
        """
        The instrument as one channel M -> (+)_j M, rho -> (+)_j Phi_j(rho): the outcome is recorded in the centre of
        the output algebra (sector (j, c) holds Phi_j(rho)_c). Its sector probabilities, summed over c, are the outcome
        probabilities; the pinching onto the register forgets nothing.
        """
        K, alg = self.kraus, self.algebra
        B, J, r, C, k = K.shape[:5]
        G = K.new_zeros(B, r, J, C, C, k, k)
        idx = torch.arange(C, device=K.device)
        G[:, :, :, idx, idx] = K.transpose(1, 2)                           # block (j, c) <- c
        return InterSectorChannel(alg, self.register_algebra, G.reshape(B, r, J * C, C, k, k),
                                  trace_preserving=self._trace_preserving)

    # ------------------------------------------------------------------
    @staticmethod
    def lueders(povm: POVM, tol: Optional[float] = None) -> 'Instrument':
        """Lueders instrument rho -> P_j rho P_j of a projective measurement (effects P_j = P_j^2)."""
        if not povm.is_projective(tol):
            raise ValueError("Instrument.lueders: the effects must be projections (use Instrument.from_povm)")
        return Instrument(povm.effects.unsqueeze(2), povm.algebra, check=False, outcomes=povm.outcomes)._tp()

    @staticmethod
    def from_povm(povm: POVM) -> 'Instrument':
        """The instrument with Kraus operators sqrt(E_j) (equal to the Lueders instrument for projections)."""
        E = povm.effects
        B, J = E.shape[:2]
        alg = povm.algebra
        op = alg.operator(E.reshape(B * J, *E.shape[2:]), is_self_adjoint=True, is_positive=True)
        K = op.sqrt().matrix.reshape(B, J, 1, *E.shape[2:])
        return Instrument(K, alg, check=False, outcomes=povm.outcomes)._tp()

    def _tp(self) -> 'Instrument':
        self._trace_preserving = True
        return self


# ======================================================================
# Fisher information
# ======================================================================
def fisher_information(prob_fn, params: Sequence[torch.Tensor], per_sample: bool = False, eps: float = 1e-12,
                       create_graph: bool = False) -> torch.Tensor:
    """
    Classical Fisher information F = J^T diag(1/p) J of the outcome distribution p = prob_fn(*params), shape
    (batch, J), with respect to the real parameters: all entries of `params`, flattened in order, complex ones as
    (real, imaginary) pairs. Returns (batch, P, P). Outcomes with p_j <= eps do not contribute.

    per_sample=True: every parameter has a leading batch dimension and sample b of prob_fn depends only on the
    parameters of sample b (one parameter point per sample); F[b] is the Fisher information of these, P being the
    number of real parameters per sample. Otherwise the parameters are shared by all samples (e.g. a batch of states
    that depend on the same theta). create_graph=True makes F differentiable with respect to tensors that prob_fn
    captures (e.g. the parameters of the measurement, to optimise it).
    """
    params = list(params)
    cplx = [torch.is_complex(t) for t in params]
    reals = tuple(torch.view_as_real(t).detach() if c else t.detach() for t, c in zip(params, cplx))
    store = {}

    def f(*xs):
        ps = [torch.complex(x[..., 0], x[..., 1]) if c else x for x, c in zip(xs, cplx)]
        p = prob_fn(*ps)
        p = p.unsqueeze(0) if p.dim() == 1 else p
        store['p'] = p.detach()
        return p.sum(dim=0) if per_sample else p

    jac = torch.autograd.functional.jacobian(f, reals, create_graph=create_graph)
    p = store['p']
    B, J = p.shape
    if per_sample:                    # d (sum_b p_b) / d x_b = d p_b / d x_b: (J, B, ...) -> (B, J, P_b)
        Jm = torch.cat([g.reshape(J, B, -1).transpose(0, 1) for g in jac], dim=-1)
    else:                             # (B, J, ...) -> (B, J, P)
        Jm = torch.cat([g.reshape(B, J, -1) for g in jac], dim=-1)
    w = torch.where(p > eps, 1.0 / torch.where(p > eps, p, torch.ones_like(p)), torch.zeros_like(p))
    return Jm.transpose(-2, -1) @ (w.unsqueeze(-1).to(Jm.dtype) * Jm)
