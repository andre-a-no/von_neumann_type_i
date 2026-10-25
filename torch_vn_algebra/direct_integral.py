"""
Continuous centres: direct integrals of matrix algebras, discretised.

A Type I von Neumann algebra with finite-dimensional fibres is a direct integral M = int^+_X M_{n(x)} dmu(x) over a
measure space (X, mu): its centre is L^oo(X, mu), its elements are bounded measurable fields x -> A(x) in M_{n(x)},
and its trace is Tr A = int Tr A(x) dmu(x). TypeIAlgebra is the case of a finite set X with point masses. A
DirectIntegral represents (X, mu) by nodes x_i in X with weights w_i > 0 -- a quadrature rule, Monte Carlo samples,
or the cells of a partition with their measures -- and works in the algebra (+)_i M_{n(x_i)}, one sector per node:

    field x -> A(x)            ->  operator (+)_i A(x_i)                            (DirectIntegral.field)
    int Tr A(x) dmu(x)         ->  sum_i w_i Tr A(x_i)                              (trace)
    int g(x) dmu(x)            ->  sum_i w_i g(x_i)                                 (integrate)
    state int Tr(rho(x) .)     ->  DensityMatrix (+)_i w_i rho(x_i), normalised      (state)
    central projection 1_S     ->  (+)_i [x_i in S] 1                               (indicator)
    ess sup of a function      ->  max over the nodes                               (ess_sup)

Everything that acts fibrewise (products, functional calculus, spectra, norms, channels that preserve the centre) is
exact at the nodes; only integrals carry the error of the rule. For fields that are constant on the cells of a
partition (simple fields) the cells are the sectors and the rule is exact. All nodes of one fibre dimension go through
one batched kernel, so thousands of nodes cost about as much as thousands of samples.

    X = DirectIntegral.torus(64, dims=2)                      # Brillouin zone [0, 2pi), measure dk/2pi, midpoint rule
    H = X.field(bloch_hamiltonian, is_self_adjoint=True)      # (+)_i H(k_i), one 2x2 sector per k-point
    e0 = X.trace(H.apply_function(lambda w: w.clamp(max=0)))  # int Tr min(H(k), 0) dk/2pi: ground energy per cell
    X = DirectIntegral.grid([(0, 1), (0, 1)], 32, dims=3, rule='gauss', density=lambda x: ...)   # mu = rho dx
    X = DirectIntegral.monte_carlo(lambda n, g: torch.rand(n, 2, generator=g), 4096, dims=2)    # mu(X) = 1
"""
import math
from typing import Callable, Optional, Sequence, Union

import torch

from .algebra import TypeIAlgebra
from .channels import _check_input
from .states import DensityMatrix

Operator = TypeIAlgebra.Operator

__all__ = ['DirectIntegral']


def _rule_1d(a: float, b: float, n: int, rule: str):
    """Nodes and weights of a one-dimensional rule on [a, b] (float64)."""
    if n < 1:
        raise ValueError("a rule needs at least one node")
    if rule == 'gauss':
        import numpy as np
        x, w = np.polynomial.legendre.leggauss(n)
        x = torch.tensor(x, dtype=torch.float64)
        w = torch.tensor(w, dtype=torch.float64)
        return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w
    if rule == 'midpoint':
        h = (b - a) / n
        return a + h * (torch.arange(n, dtype=torch.float64) + 0.5), torch.full((n,), h, dtype=torch.float64)
    if rule == 'periodic':                  # trapezoid on a period: nodes a + j h, spectrally accurate for smooth f
        h = (b - a) / n
        return a + h * torch.arange(n, dtype=torch.float64), torch.full((n,), h, dtype=torch.float64)
    raise ValueError(f"unknown rule {rule!r} (use 'gauss', 'midpoint' or 'periodic')")


class DirectIntegral:
    """
    A measure space (X, mu) given by nodes x_i (tensor (C, d), or (C,) for d = 1) with weights w_i > 0, and fibre
    dimensions n_i = n(x_i) (one int for all nodes, or one per node). The algebra (+)_i M_{n_i} is `self.algebra`
    (precision and device as for TypeIAlgebra; the sector charges are the node indices).
    """

    def __init__(self, nodes, weights, dims: Union[int, Sequence[int]], precision: str = 'double',
                 device: Optional[torch.device] = None):
        x = torch.as_tensor(nodes, dtype=torch.float64)
        x = x.unsqueeze(-1) if x.dim() == 1 else x
        w = torch.as_tensor(weights, dtype=torch.float64).reshape(-1)
        if x.dim() != 2 or x.shape[0] != w.shape[0]:
            raise ValueError(f"DirectIntegral: need nodes of shape (C, d) and C weights, got {tuple(x.shape)} and "
                             f"{tuple(w.shape)}")
        if bool((w <= 0).any()):
            raise ValueError("DirectIntegral: the weights must be positive (drop nodes of measure zero)")
        C = x.shape[0]
        dims = [int(dims)] * C if isinstance(dims, int) else [int(n) for n in dims]
        if len(dims) != C or min(dims) < 1:
            raise ValueError(f"DirectIntegral: need one fibre dimension >= 1 per node ({C}), got {dims[:5]}...")
        self.algebra = TypeIAlgebra(dims, dims, precision=precision, device=device)
        dev = self.algebra.hilbert.device
        real = self.algebra.hilbert.real_dtype
        self.nodes = x.to(device=dev, dtype=real)
        self.weights = w.to(device=dev, dtype=real)
        self.dims = dims

    # ------------------------------------------------------------------ constructors
    @staticmethod
    def grid(box: Sequence, n: Union[int, Sequence[int]], dims: Union[int, Callable], rule: str = 'gauss',
             density: Optional[Callable[[torch.Tensor], torch.Tensor]] = None, **kw) -> 'DirectIntegral':
        """
        Product rule on the box prod_j [a_j, b_j] (box = [(a_1, b_1), ...]) with n_j nodes per direction: rule 'gauss'
        (Gauss-Legendre), 'midpoint' or 'periodic' (trapezoid on a period). mu = density(x) dx (Lebesgue measure by
        default; nodes where the density vanishes are dropped). dims: an int, or a function of the nodes (C, d)
        returning the fibre dimension of each node.
        """
        box = [tuple(map(float, ab)) for ab in box]
        ns = [int(n)] * len(box) if isinstance(n, int) else [int(m) for m in n]
        rules = [_rule_1d(a, b, m, rule) for (a, b), m in zip(box, ns)]
        mesh = torch.meshgrid(*[r[0] for r in rules], indexing='ij')
        wmesh = torch.meshgrid(*[r[1] for r in rules], indexing='ij')
        x = torch.stack([m.reshape(-1) for m in mesh], dim=-1)
        w = torch.ones(x.shape[0], dtype=torch.float64)
        for m in wmesh:
            w = w * m.reshape(-1)
        if density is not None:
            rho = torch.as_tensor(density(x), dtype=torch.float64).reshape(-1)
            if bool((rho < 0).any()):
                raise ValueError("grid: the density must be non-negative")
            w = w * rho
            keep = w > 0
            x, w = x[keep], w[keep]
        d = dims(x) if callable(dims) else dims
        d = d if isinstance(d, int) else [int(v) for v in torch.as_tensor(d).reshape(-1)]
        return DirectIntegral(x, w, d, **kw)

    @staticmethod
    def torus(n: Union[int, Sequence[int]], dims: Union[int, Callable], d: int = 1, **kw) -> 'DirectIntegral':
        """The torus [0, 2 pi)^d with the normalised Haar measure dk/(2 pi)^d (a Brillouin zone), n nodes per
        direction (the momenta 2 pi j / n of a periodic chain of n cells), periodic trapezoid rule. Note: `dims` is
        the fibre dimension n(x), `d` the dimension of the torus: torus(64, dims=2) is a circle with 2x2 fibres."""
        if 'rule' in kw:
            raise TypeError("torus: the rule is the periodic trapezoid rule; use grid for other rules")
        ns = [int(n)] * d if isinstance(n, int) else list(n)
        X = DirectIntegral.grid([(0.0, 2 * math.pi)] * len(ns), ns, dims, rule='periodic', **kw)
        X.weights = X.weights / (2 * math.pi) ** len(ns)
        return X

    @staticmethod
    def monte_carlo(sampler: Callable[[int, Optional[torch.Generator]], torch.Tensor], samples: int,
                    dims: Union[int, Callable], mass: float = 1.0, generator: Optional[torch.Generator] = None,
                    **kw) -> 'DirectIntegral':
        """Nodes sampled from the probability measure mu/mass (sampler(n, generator) -> (n, d) or (n,)), weights
        mass / samples: integrals are Monte Carlo estimates, with error of order samples^(-1/2)."""
        x = torch.as_tensor(sampler(samples, generator), dtype=torch.float64)
        x = x.unsqueeze(-1) if x.dim() == 1 else x
        w = torch.full((samples,), float(mass) / samples, dtype=torch.float64)
        d = dims(x) if callable(dims) else dims
        d = d if isinstance(d, int) else [int(v) for v in torch.as_tensor(d).reshape(-1)]
        return DirectIntegral(x, w, d, **kw)

    # ------------------------------------------------------------------ fields
    @property
    def C(self) -> int:
        return self.algebra.C

    def __repr__(self) -> str:
        return (f"DirectIntegral(nodes={self.C}, d={self.nodes.shape[-1]}, fibres={sorted(set(self.dims))}, "
                f"mass={float(self.weights.sum()):.6g})")

    def field(self, f: Callable[[torch.Tensor], torch.Tensor], **properties) -> Operator:
        """
        The operator (+)_i A(x_i) of a field of matrices. If all fibres have the same dimension n, f receives all nodes
        (C, d) and returns (C, n, n) or (batch, C, n, n); otherwise f is called per node with a tensor (d,) and returns
        (n_i, n_i) or (batch, n_i, n_i). Keyword arguments tag properties (is_self_adjoint=True, ...).
        """
        alg = self.algebra
        dt = alg.dtype
        if len(set(self.dims)) == 1:
            m = torch.as_tensor(f(self.nodes)).to(device=self.nodes.device, dtype=dt)
            m = m.unsqueeze(0) if m.dim() == 3 else m
            n = self.dims[0]
            if tuple(m.shape[1:]) != (self.C, n, n):
                raise ValueError(f"field: f must return ({self.C}, {n}, {n}) or (batch, {self.C}, {n}, {n}), "
                                 f"got {tuple(m.shape)}")
            return alg.operator(m, **properties)
        blocks = []
        for i in range(self.C):
            b = torch.as_tensor(f(self.nodes[i])).to(device=self.nodes.device, dtype=dt)
            b = b.unsqueeze(0) if b.dim() == 2 else b
            if tuple(b.shape[1:]) != (self.dims[i], self.dims[i]):
                raise ValueError(f"field: node {i} has fibre dimension {self.dims[i]}, f returned {tuple(b.shape)}")
            blocks.append(b)
        B = max(b.shape[0] for b in blocks)
        m = torch.zeros(B, self.C, alg.k_max, alg.k_max, dtype=dt, device=self.nodes.device)
        for i, b in enumerate(blocks):
            m[:, i, :self.dims[i], :self.dims[i]] = b
        return alg.operator(m, **properties)

    def scalar(self, g: Callable[[torch.Tensor], torch.Tensor]) -> Operator:
        """The central element (+)_i g(x_i) 1 of a scalar function (g maps (C, d) to (C,) or (batch, C))."""
        v = torch.as_tensor(g(self.nodes)).to(self.nodes.device)
        v = v.unsqueeze(0) if v.dim() == 1 else v
        dt = self.algebra.dtype
        eye = self.algebra._const('eye', dt)
        return self.algebra.operator(v.to(dt)[..., None, None] * eye, is_self_adjoint=not torch.is_complex(v))

    def indicator(self, predicate: Callable[[torch.Tensor], torch.Tensor]) -> Operator:
        """The central projection 1_S of S = {x : predicate(x)} (predicate maps (C, d) to a boolean (C,))."""
        mask = torch.as_tensor(predicate(self.nodes)).to(self.nodes.device)
        if mask.dtype != torch.bool or mask.shape != (self.C,):
            raise ValueError(f"indicator: the predicate must return a boolean tensor of shape ({self.C},)")
        dt = self.algebra.dtype
        m = mask.to(dt)[None, :, None, None] * self.algebra._const('eye', dt)
        return self.algebra.operator(m, is_self_adjoint=True, is_positive=True, is_projection=True)

    # ------------------------------------------------------------------ integrals
    def integrate(self, values: torch.Tensor) -> torch.Tensor:
        """sum_i w_i g(x_i) for values g(x_i) of shape (C,) or (batch, C): shape () or (batch,)."""
        v = torch.as_tensor(values).to(self.nodes.device)
        if v.shape[-1] != self.C:
            raise ValueError(f"integrate: need one value per node ({self.C}), got shape {tuple(v.shape)}")
        dt = torch.promote_types(v.dtype, self.weights.dtype)        # integer values must not truncate the weights
        return (v.to(dt) * self.weights.to(dt)).sum(-1)

    def measure(self, predicate: Optional[Callable[[torch.Tensor], torch.Tensor]] = None) -> torch.Tensor:
        """mu(S) for S = {x : predicate(x)} (mu(X) without a predicate)."""
        if predicate is None:
            return self.weights.sum()
        mask = torch.as_tensor(predicate(self.nodes)).to(self.nodes.device)
        if mask.dtype != torch.bool or mask.shape != (self.C,):
            raise ValueError(f"measure: the predicate must return a boolean tensor of shape ({self.C},)")
        return self.weights[mask].sum()

    def fibre_traces(self, A: Operator) -> torch.Tensor:
        """Tr A(x_i), shape (batch, C)."""
        _check_input(A, self.algebra, "DirectIntegral")
        return torch.diagonal(A.matrix, dim1=-2, dim2=-1).sum(-1)

    def trace(self, A: Operator) -> torch.Tensor:
        """The trace int Tr A(x) dmu(x) ~ sum_i w_i Tr A(x_i), shape (batch,)."""
        return self.integrate(self.fibre_traces(A))

    def normalized_trace(self, A: Operator) -> torch.Tensor:
        """int tr A(x) dmu(x) with the normalised fibre trace tr = Tr / n(x), shape (batch,) (a tracial state if
        mu(X) = 1)."""
        n = torch.tensor(self.dims, dtype=self.weights.dtype, device=self.nodes.device)
        return self.integrate(self.fibre_traces(A) / n)

    def state(self, density: Optional[Callable[[torch.Tensor], torch.Tensor]] = None) -> DensityMatrix:
        """
        The normal state phi(A) = int Tr(rho(x) A(x)) dmu(x) / int Tr rho dmu of a field of positive matrices rho
        (called like `field`), as a DensityMatrix of the algebra: blocks w_i rho(x_i), normalised. Without a density:
        the tracial state int tr A dmu / mu(X).
        """
        if density is None:
            n = torch.tensor(self.dims, dtype=self.weights.dtype, device=self.nodes.device)
            scal = (self.weights / n).to(self.algebra.dtype)
            m = scal[None, :, None, None] * self.algebra._const('eye', self.algebra.dtype)
        else:
            m = self.field(density).matrix * self.weights.to(self.algebra.dtype)[None, :, None, None]
        tot = torch.diagonal(m, dim1=-2, dim2=-1).sum((-2, -1)).real
        return DensityMatrix(self.algebra, matrix=m / tot.to(m.dtype)[:, None, None, None], validate=False)

    def ess_sup(self, values: torch.Tensor) -> torch.Tensor:
        """The essential supremum of a function given by its node values (C,) or (batch, C): the maximum over the
        nodes (all of which have positive weight)."""
        return torch.as_tensor(values).amax(dim=-1)
