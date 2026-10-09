"""
Batched gradient-based optimisation over operators of a Type I algebra.

Monte Carlo sampling explores typical operators; extremal values of a functional (the
sharp constants in an inequality) sit on sets of tiny measure that random sampling almost never
hits in higher dimensions. This module parametrises the usual constraint sets differentiably,
so that sup / inf of a functional can be searched for with many independent starts at once
(one start per batch element):

    UnitaryParam            U = U0 exp(G - G^*)  (Haar-random U0 per start, chart re-centred)
    SpectrumParam           spectra in [lo, 1] with lambda_max = 1 and lambda_min = lo, i.e. a
                            prescribed (or bounded) Michelson contrast
    PositiveParam           X = U diag(lambda) U^*  with such a spectrum
    SelfAdjointParam        X = U diag(s * lambda) U^*  with fixed random signs s per start

    best = optimize.extremize(objective, params, maximize=True)

The objective maps the current operators to one value per start, shape (batch,). The search is
local: the result is a lower bound for a supremum (an upper bound for an infimum) and should be
checked with more starts / steps.
"""
from typing import Callable, List, Optional, Sequence

import torch

from .algebra import TypeIAlgebra
from . import cost

Operator = TypeIAlgebra.Operator


class Param:
    """Base class: a differentiable parametrisation of a batch of operators."""

    def parameters(self) -> List[torch.Tensor]:
        raise NotImplementedError

    def recenter(self):
        """Move the reference point to the current point (default: nothing)."""


class UnitaryParam(Param):
    """Block-diagonal unitaries U = U0 exp(G - G^*), one per start; U0 Haar-random by default."""

    def __init__(self, alg: TypeIAlgebra, batch_size: int, base: Optional[torch.Tensor] = None):
        self.alg = alg
        self.base = alg.random_unitary_operator(batch_size).matrix if base is None else base
        self.G = torch.zeros_like(self.base, requires_grad=True)

    def parameters(self):
        return [self.G]

    def operator(self) -> Operator:
        skew = self.G - self.G.conj().transpose(-2, -1)
        return self.alg.operator(self.base) @ self.alg.operator(skew).expm()

    def recenter(self):
        self.base = self.operator().matrix.detach()
        with torch.no_grad():
            self.G.zero_()


def _contrast_to_lo(delta: torch.Tensor) -> torch.Tensor:
    return (1 - delta) / (1 + delta)


class SpectrumParam(Param):
    """
    Spectra (per block, padded) in [lo, 1] with the global maximum 1 and minimum lo, so that the
    Michelson contrast equals delta = (1 - lo) / (1 + lo).

    delta: tensor (batch,) of prescribed contrasts, or None to optimise the contrast as well
    within [delta_range[0], delta_range[1]]. The maximum and the minimum are pinned to the first
    two eigenvalues of sector `pin_sector` (which needs k >= 2); for algebras whose sectors all have
    the same size this loses no generality.
    """

    def __init__(self, alg: TypeIAlgebra, batch_size: int, delta: Optional[torch.Tensor] = None,
                 delta_range=(0.0, 0.999), pin_sector: Optional[int] = None):
        if pin_sector is None:                     # default: the first sector with at least two dimensions
            pin_sector = next((c for c, k in enumerate(alg.k_factors) if k >= 2), 0)
        if alg.k_factors[pin_sector] < 2:
            raise ValueError("the pinned sector needs at least two dimensions")
        self.alg = alg
        rdt, dev = alg.hilbert.real_dtype, alg.hilbert.device
        self.theta = torch.randn(batch_size, alg.C, alg.k_max, dtype=rdt, device=dev).requires_grad_(True)
        self.fixed_delta = None if delta is None else torch.as_tensor(delta, dtype=rdt, device=dev).expand(batch_size)
        self.delta_range = delta_range
        self.eta = None if delta is not None else torch.randn(batch_size, dtype=rdt, device=dev).requires_grad_(True)
        mask = torch.zeros(alg.C, alg.k_max, dtype=torch.bool, device=dev)
        for c, k_c in enumerate(alg.k_factors):
            mask[c, :k_c] = True
        self.mask = mask
        self.pin = pin_sector

    def parameters(self):
        return [self.theta] + ([self.eta] if self.eta is not None else [])

    def delta(self) -> torch.Tensor:
        if self.fixed_delta is not None:
            return self.fixed_delta
        a, b = self.delta_range
        return a + (b - a) * torch.sigmoid(self.eta)

    def eigenvalues(self) -> torch.Tensor:
        """(batch, C, k_max), zero on padding."""
        lo = _contrast_to_lo(self.delta())[:, None, None]
        lam = lo + (1 - lo) * torch.sigmoid(self.theta)
        pin = torch.zeros_like(lam, dtype=torch.bool)
        pin[:, self.pin, :2] = True
        pinned = torch.zeros_like(lam)
        pinned[:, self.pin, 0] = 1.0
        pinned[:, self.pin, 1] = lo[:, 0, 0]
        lam = torch.where(pin, pinned, lam)
        return lam * self.mask


class PositiveParam(Param):
    """Positive operators X = U diag(lambda) U^* with a prescribed or bounded Michelson contrast."""

    def __init__(self, alg: TypeIAlgebra, batch_size: int, delta: Optional[torch.Tensor] = None,
                 delta_range=(0.0, 0.999), pin_sector: Optional[int] = None):
        self.alg = alg
        self.spectrum = SpectrumParam(alg, batch_size, delta, delta_range, pin_sector)
        self.unitary = UnitaryParam(alg, batch_size)

    def parameters(self):
        return self.spectrum.parameters() + self.unitary.parameters()

    def recenter(self):
        self.unitary.recenter()

    def _signed(self, lam):
        return lam

    def _build(self, lam: torch.Tensor) -> Operator:
        U = self.unitary.operator().matrix
        lam = lam.to(U.dtype)
        return self.alg.operator((U * lam.unsqueeze(-2)) @ U.conj().transpose(-2, -1))

    def operator(self) -> Operator:
        return self._build(self._signed(self.spectrum.eigenvalues()))

    def abs_operator(self) -> Operator:
        """|X| = U diag(lambda) U^*, exact (no SVD, so gradients stay finite)."""
        return self._build(self.spectrum.eigenvalues())

    def contrast(self) -> torch.Tensor:
        return self.spectrum.delta()


class SelfAdjointParam(PositiveParam):
    """
    Self-adjoint X = U diag(s * lambda) U^* where |X| has a prescribed or bounded contrast and the
    signs s = +-1 are drawn at random per start and kept fixed (discrete variables; the multistart
    covers the sign patterns).
    """

    def __init__(self, alg: TypeIAlgebra, batch_size: int, delta: Optional[torch.Tensor] = None,
                 delta_range=(0.0, 0.999), pin_sector: Optional[int] = None):
        super().__init__(alg, batch_size, delta, delta_range, pin_sector)
        dev = alg.hilbert.device
        self.signs = (2 * torch.randint(0, 2, (batch_size, alg.C, alg.k_max), device=dev) - 1).to(alg.hilbert.real_dtype)

    def _signed(self, lam):
        return lam * self.signs


def extremize(objective: Callable[[], torch.Tensor], params: Sequence[Param], maximize: bool = True,
              rounds: int = 10, steps: int = 50, lr: float = 0.05, decay: float = 0.6,
              progress: bool = False) -> torch.Tensor:
    """
    Gradient ascent (maximize=True) or descent on objective() -> tensor (batch,), every batch
    element being an independent start. Adam with learning rate lr * decay**round; unitary charts
    are re-centred after every round. Returns the objective at the final point (no gradient).
    """
    sign = -1.0 if maximize else 1.0
    timer = cost.StepTimer(rounds * steps, "optimize.extremize", progress)
    for rnd in range(rounds):
        tensors = [t for p in params for t in p.parameters()]
        opt = torch.optim.Adam(tensors, lr=lr * decay ** rnd)
        for _ in range(steps):
            loss = sign * objective().sum()
            opt.zero_grad()
            loss.backward()
            opt.step()
            timer.step(tensors[0].device)
        for p in params:
            p.recenter()
    timer.close()
    with torch.no_grad():
        return objective().detach()
