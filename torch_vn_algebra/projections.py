"""
The projection lattice of M = (+)_c M_{k_c}, measures on it, and spectral measures.

Lattice. P v Q and P ^ Q are the projections onto ran P + ran Q and ran P n ran Q. For projections P, Q the operator
P + Q has spectrum in [0, 2], P v Q is its support and P ^ Q its spectral projection for the eigenvalue 2 (since
(P + Q)x = 2x forces Px = Qx = x). Both are computed blockwise by functional calculus with an indicator function,
to a tolerance `tol` on the eigenvalues of P + Q: when P and Q are close (principal angle theta), P + Q has the
eigenvalues 1 +- cos(theta), so ^ and v are only defined to angles above about sqrt(2 tol).

Measures. A measure on the projections of M is a function mu >= 0 with mu(P + Q) = mu(P) + mu(Q) whenever PQ = 0.
Since P = sum_c P_c, mu = sum_c mu_c with mu_c a measure on the projections of M_{k_c}. By Gleason's theorem
mu_c(P) = Tr(rho_c P) for a positive rho_c when k_c >= 3, and trivially when k_c = 1. On a 2x2 block mu_c is any frame
function: mu_c(1) = w_c and, for rank-one P = (1 + n.sigma)/2 with Bloch vector n, mu_c(P) = w_c/2 + g_c(n) with g_c
odd and |g_c| <= w_c/2, linear (a state) only if g_c is. ProjectionMeasure stores a positive operator rho (the linear
part, mu(1) = Tr rho) and, for sectors of size 2, optional odd functions g_c that replace the linear part there.
Consequences: a measure is a state on every algebra without 2x2 blocks; it is unitarily invariant (a weighted
dimension function, i.e. a trace restricted to projections) iff rho is central and every g_c vanishes; and a
{0, 1}-valued (dispersion-free) measure exists iff some sector has k_c <= 2 (Kochen-Specker for k_c >= 3).

Spectral measures. spectral_measure(A) returns the projection-valued measure of a self-adjoint A, a PVM: a POVM of
mutually orthogonal projections E_j with real values lambda_j, A = sum_j lambda_j E_j (measurement.PVM).

    P = projections.random_projection(alg, ranks=[0, 1, 2], batch_size=8)    # rank r_c in sector c, Haar random
    PvQ, PQ = projections.join(P, Q), projections.meet(P, Q)
    mu = projections.ProjectionMeasure.from_state(rho)                      # mu(P) = Tr(rho P)
    mu = projections.ProjectionMeasure(alg.identity(), frames={1: g})       # non-linear on the 2x2 sector 1
    mu(P), mu.modular_defect(P, Q)                                          # mu(PvQ) + mu(P^Q) - mu(P) - mu(Q)
    E = projections.spectral_measure(A)                                     # PVM: E.values, E[j], E.observable()
"""
from typing import Callable, Dict, Optional, Sequence

import torch

from . import cost
from .algebra import TypeIAlgebra
from .channels import _check_input
from .measurement import PVM, _default_tol

Operator = TypeIAlgebra.Operator

__all__ = ['complement', 'support', 'join', 'meet', 'ranks', 'random_projection', 'bloch_vector', 'ProjectionMeasure',
           'spectral_measure']


def _require_projection(P: Operator, what: str) -> None:
    if not P.is_projection:
        raise ValueError(f"{what}: expected a projection (P = P* = P^2)")


def _indicator(A: Operator, threshold: float) -> Operator:
    """Spectral projection of the self-adjoint A onto the eigenvalues > threshold."""
    out = A.apply_function(lambda w: (w > threshold).to(w.dtype))
    return A.algebra.operator(out.matrix, is_self_adjoint=True, is_positive=True, is_projection=True)


# ======================================================================
# lattice
# ======================================================================
def complement(P: Operator) -> Operator:
    """The orthocomplement 1 - P."""
    _require_projection(P, "complement")
    alg = P.algebra
    m = alg._const('eye', P.matrix.dtype) - P.matrix
    return alg.operator(m, is_self_adjoint=True, is_positive=True, is_projection=True)


def support(A: Operator, tol: Optional[float] = None) -> Operator:
    """The support (range projection) of a positive operator: its spectral projection onto the eigenvalues
    > tol * ||A|| (relative tolerance, per sample). Raises ValueError for an operator with eigenvalues below
    -tol * ||A|| (not positive)."""
    _check_input(A, A.algebra, "support")
    if not A.is_self_adjoint:
        raise ValueError("support: expected a self-adjoint (positive) operator")
    tol = _default_tol(A.matrix.dtype) if tol is None else tol
    norm = A.operator_norm()
    if bool((A.lambda_min < -tol * norm).any()):
        raise ValueError("support: expected a positive operator")
    B = A.matrix.shape[0]
    thr = (tol * norm).reshape(B, 1)
    out = A.apply_function(lambda w: (w > thr.to(w.dtype)).to(w.dtype))
    return A.algebra.operator(out.matrix, is_self_adjoint=True, is_positive=True, is_projection=True)


def join(P: Operator, Q: Operator, tol: Optional[float] = None) -> Operator:
    """P v Q, the projection onto ran P + ran Q: the support of P + Q (eigenvalues > tol)."""
    _check_input(Q, P.algebra, "join")
    _require_projection(P, "join")
    _require_projection(Q, "join")
    tol = _default_tol(P.matrix.dtype) if tol is None else tol
    S = P.algebra.operator(P.matrix + Q.matrix, is_self_adjoint=True, is_positive=True)
    return _indicator(S, tol)


def meet(P: Operator, Q: Operator, tol: Optional[float] = None) -> Operator:
    """P ^ Q, the projection onto ran P n ran Q: the spectral projection of P + Q for the eigenvalues > 2 - tol."""
    _check_input(Q, P.algebra, "meet")
    _require_projection(P, "meet")
    _require_projection(Q, "meet")
    tol = _default_tol(P.matrix.dtype) if tol is None else tol
    S = P.algebra.operator(P.matrix + Q.matrix, is_self_adjoint=True, is_positive=True)
    return _indicator(S, 2 - tol)


def ranks(P: Operator) -> torch.Tensor:
    """Rank of every block of a projection, shape (batch, C), int64 (the rounded block traces)."""
    _require_projection(P, "ranks")
    t = torch.diagonal(P.matrix, dim1=-2, dim2=-1).sum(-1).real
    return torch.round(t).to(torch.long)


def random_projection(alg: TypeIAlgebra, ranks: Sequence[int], batch_size: int = 1) -> Operator:
    """Haar random projections with rank ranks[c] in sector c: U diag(1, .., 1, 0, .., 0) U*, U Haar in every block."""
    ranks = list(ranks)
    if len(ranks) != alg.C:
        raise ValueError(f"random_projection: need one rank per sector ({alg.C}), got {len(ranks)}")
    if any(not 0 <= r <= k for r, k in zip(ranks, alg.k_factors)):
        raise ValueError(f"random_projection: ranks {ranks} must lie between 0 and the block sizes {alg.k_factors}")
    U = alg.random_unitary_operator(batch_size).matrix
    d = torch.zeros(alg.C, alg.k_max, dtype=U.dtype, device=U.device)
    for c, r in enumerate(ranks):
        d[c, :r] = 1
    m = (U * d.unsqueeze(-2)) @ U.conj().transpose(-2, -1)
    return alg.operator(m, is_self_adjoint=True, is_positive=True, is_projection=True)


def _paulis(dtype, device) -> torch.Tensor:
    m = torch.tensor([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], dtype=torch.complex128)
    return m.to(device=device, dtype=dtype if dtype.is_complex else torch.complex128)


def bloch_vector(P: Operator, sector: int) -> torch.Tensor:
    """Bloch vector n_i = Tr(P_c sigma_i) of the 2x2 block c of P, shape (batch, 3) (a unit vector if P_c has rank
    one)."""
    alg = P.algebra
    if alg.k_factors[sector] != 2:
        raise ValueError(f"bloch_vector: sector {sector} has size {alg.k_factors[sector]}, not 2")
    blk = P.matrix[:, sector, :2, :2]
    s = _paulis(blk.dtype, blk.device)
    return torch.einsum('bij,sji->bs', blk.to(s.dtype), s).real


# ======================================================================
# measures on the projection lattice
# ======================================================================
class ProjectionMeasure:
    """
    A finitely additive measure on the projections of M (see the module docstring).

    density: a positive Operator rho (batch 1 or B); mu(P) = Tr(rho P) on every sector without a frame function.
    frames: {sector: g} for sectors of size 2, g mapping Bloch vectors (..., 3) to real values (...,), odd and with
    |g| <= w_c/2, w_c = Tr rho_c; on such a sector mu(P_c) = 0, w_c/2 + g(n) or w_c for P_c of rank 0, 1, 2.
    check=True tests positivity of rho and, on `samples` random unit vectors, oddness and the bound of every g.
    """

    def __init__(self, density: Operator, frames: Optional[Dict[int, Callable[[torch.Tensor], torch.Tensor]]] = None,
                 check: bool = True, samples: int = 256, tol: Optional[float] = None):
        alg = density.algebra
        self.algebra = alg
        self.density = density
        self.frames = dict(frames or {})
        for c in self.frames:
            if not 0 <= c < alg.C or alg.k_factors[c] != 2:
                raise ValueError(f"ProjectionMeasure: frame functions only on sectors of size 2; sector {c} has size "
                                 f"{alg.k_factors[c] if 0 <= c < alg.C else 'none'}")
        if check:
            self._check(samples, tol)

    def _check(self, samples: int, tol: Optional[float]) -> None:
        rho = self.density
        tol = _default_tol(rho.matrix.dtype) if tol is None else tol
        if not rho.is_self_adjoint:
            raise ValueError("ProjectionMeasure: the density must be self-adjoint")
        if min(float(w.min()) for w, k in zip(rho.eigenvalues(), self.algebra.k_factors) if k) < -tol:
            raise ValueError("ProjectionMeasure: the density must be positive")
        if not self.frames:
            return
        g = torch.Generator().manual_seed(0)
        n = torch.randn(samples, 3, generator=g, dtype=torch.float64)
        n = (n / n.norm(dim=-1, keepdim=True)).to(device=rho.algebra.hilbert.device,
                                                  dtype=rho.algebra.hilbert.real_dtype)
        w = self.weights()
        for c, f in self.frames.items():
            a, b = f(n), f(-n)
            if bool(((a + b).abs() > tol).any()):
                raise ValueError(f"ProjectionMeasure: the frame function of sector {c} must be odd, g(-n) = -g(n)")
            if bool((a.abs() > 0.5 * w[:, c].abs().amin() + tol).any()):
                raise ValueError(f"ProjectionMeasure: the frame function of sector {c} must satisfy |g| <= w_c/2")

    # ------------------------------------------------------------------
    @staticmethod
    def from_state(rho: Operator) -> 'ProjectionMeasure':
        """The linear measure mu(P) = Tr(rho P) of a state (or any positive operator)."""
        return ProjectionMeasure(rho, check=False)

    @staticmethod
    def dimension(algebra: TypeIAlgebra, weights: Optional[Sequence[float]] = None) -> 'ProjectionMeasure':
        """The weighted dimension function mu(P) = sum_c w_c rank(P_c) (a trace on projections; w_c = 1 by default)."""
        weights = [1.0] * algebra.C if weights is None else list(weights)
        return ProjectionMeasure(algebra.central(weights), check=False)

    @staticmethod
    def dispersion_free(algebra: TypeIAlgebra, sector: int, axis=(0.0, 0.0, 1.0)) -> 'ProjectionMeasure':
        """
        A {0, 1}-valued measure concentrated on `sector`: mu(1_c) = 1 and, on a 2x2 block, mu(P) = 1 for rank-one P
        with Bloch vector in the closed half-sphere of `axis` (ties broken by a fixed second and third axis). It exists
        only for k_c <= 2: for k_c >= 3 there is none (Kochen-Specker), and a ValueError is raised.
        """
        k = algebra.k_factors[sector]
        if k >= 3:
            raise ValueError(f"dispersion_free: sector {sector} has size {k} >= 3; by the Kochen-Specker theorem its "
                             "projections carry no {0, 1}-valued measure")
        w = [0.0] * algebra.C
        w[sector] = 1.0 / k                                 # mu(1_c) = Tr rho_c = 1
        rho = algebra.central(w)
        if k < 2:
            return ProjectionMeasure(rho, check=False)
        a = torch.tensor(axis, dtype=torch.float64)
        a = a / a.norm()
        b = torch.tensor([1.0, 0.0, 0.0], dtype=torch.float64) if abs(float(a[0])) < 0.9 else \
            torch.tensor([0.0, 1.0, 0.0], dtype=torch.float64)
        b = b - (b @ a) * a
        b = b / b.norm()
        e = torch.linalg.cross(a, b)
        frame = torch.stack([a, b, e])                                             # (3, 3), rows a, b, e

        def g(n):
            u = torch.movedim(n.to(torch.float64) @ frame.to(n.device).T, -1, 0)    # coordinates (3, ...)
            eps = 1e-12
            s = torch.where(u[0].abs() > eps, u[0].sign(), torch.where(u[1].abs() > eps, u[1].sign(), u[2].sign()))
            return (0.5 * s).to(n.dtype)
        return ProjectionMeasure(rho, frames={sector: g}, check=False)

    # ------------------------------------------------------------------
    def weights(self) -> torch.Tensor:
        """mu(1_c) = Tr rho_c, shape (batch, C), real."""
        return torch.diagonal(self.density.matrix, dim1=-2, dim2=-1).sum(-1).real

    def total(self) -> torch.Tensor:
        """mu(1), shape (batch,)."""
        return self.weights().sum(-1)

    def __call__(self, P: Operator) -> torch.Tensor:
        """mu(P) for a projection P (possibly a batch), shape (batch,), real."""
        _check_input(P, self.algebra, "ProjectionMeasure")
        _require_projection(P, "ProjectionMeasure")
        rho = self.density.matrix
        prod = rho.to(P.matrix.dtype) * P.matrix.transpose(-2, -1)
        per = prod.sum(dim=(-2, -1)).real                                     # (B, C): Tr(rho_c P_c)
        if self.frames:
            r = ranks(P)
            w = self.weights().to(per.dtype)
            per = per.clone()
            for c, g in self.frames.items():
                val = 0.5 * w[:, c] + g(bloch_vector(P, c)).to(per.dtype)
                per[:, c] = torch.where(r[:, c] == 0, torch.zeros_like(val),
                                        torch.where(r[:, c] == 2, w[:, c].expand_as(val), val))
        return per.sum(-1)

    def modular_defect(self, P: Operator, Q: Operator, tol: Optional[float] = None) -> torch.Tensor:
        """mu(P v Q) + mu(P ^ Q) - mu(P) - mu(Q): zero for every P, Q iff mu is a weighted dimension function."""
        return self(join(P, Q, tol)) + self(meet(P, Q, tol)) - self(P) - self(Q)

    def is_linear(self, samples: int = 256, tol: Optional[float] = None) -> bool:
        """True if mu is a state (Tr(rho' P) for some positive rho'): every frame function is linear on `samples`
        random Bloch vectors (Gleason: the other sectors are linear anyway)."""
        try:
            self.gleason_state(samples, tol)
            return True
        except ValueError:
            return False

    def gleason_state(self, samples: int = 256, tol: Optional[float] = None) -> Operator:
        """
        The positive operator rho' with mu(P) = Tr(rho' P) for all P (Gleason). On sectors with a frame function g,
        g(n) = a.n is fitted by least squares on `samples` random Bloch vectors and rho'_c = w_c/2 + a.sigma; a
        ValueError is raised if the fit misses by more than tol (then mu is not a state).
        """
        rho = self.density
        if not self.frames:
            return rho
        tol = _default_tol(rho.matrix.dtype) if tol is None else tol
        gen = torch.Generator().manual_seed(1)
        n = torch.randn(samples, 3, generator=gen, dtype=torch.float64)
        n = n / n.norm(dim=-1, keepdim=True)
        m = rho.matrix.clone()
        for c, g in self.frames.items():
            y = g(n.to(device=m.device, dtype=rho.algebra.hilbert.real_dtype)).to(device='cpu', dtype=torch.float64)
            a = torch.linalg.lstsq(n, y.unsqueeze(-1)).solution.squeeze(-1).to(m.device)    # (3,)
            if float((n @ a - y).abs().max()) > tol:
                raise ValueError(f"gleason_state: the frame function of sector {c} is not linear; mu is not a state")
            s = _paulis(m.dtype, m.device)
            w = self.weights()[:, c].to(m.dtype)
            blk = 0.5 * w.view(-1, 1, 1) * torch.eye(2, dtype=m.dtype, device=m.device) + \
                torch.einsum('s,sij->ij', a.to(m.dtype), s)
            m[:, c, :2, :2] = blk
        return self.algebra.operator(m, is_self_adjoint=True, is_positive=True)

    def is_dimension_function(self, samples: int = 256, tol: Optional[float] = None) -> bool:
        """True if mu is unitarily invariant, mu(P) = sum_c (w_c/k_c) rank(P_c): a state with central rho (Gleason fit
        on the sectors with a frame function), i.e. a trace restricted to the projections."""
        try:
            rho = self.gleason_state(samples, tol)
        except ValueError:
            return False
        tol = _default_tol(rho.matrix.dtype) if tol is None else tol
        m = rho.matrix
        k = torch.tensor(self.algebra.k_factors, device=m.device).clamp(min=1)
        mean = torch.diagonal(m, dim1=-2, dim2=-1).sum(-1) / k
        dev = m - mean[..., None, None] * self.algebra._const('eye', m.dtype)
        return bool(dev.abs().amax() <= tol)

    def __repr__(self) -> str:
        return f"ProjectionMeasure(C={self.algebra.C}, frames={sorted(self.frames)})"


# ======================================================================
# spectral measures
# ======================================================================
def spectral_measure(A: Operator, tol: Optional[float] = None, group: bool = True) -> PVM:
    """
    The spectral measure of a self-adjoint A as a PVM. group=True (batch 1): one outcome per distinct eigenvalue
    (eigenvalues within tol of the smallest of their group are merged, across sectors; default tol
    sqrt(eps) * ||A||), E_j the spectral projection, values ascending.
    group=False (any batch): one rank-one outcome per eigenvector, sector after sector, values (batch, sum_c k_c).
    """
    _check_input(A, A.algebra, "spectral_measure")
    if not A.is_self_adjoint:
        raise ValueError("spectral_measure: A must be self-adjoint")
    alg = A.algebra
    B = A.matrix.shape[0]
    dec = A.eigh()
    dtype = A.matrix.dtype
    dev = A.matrix.device
    J = sum(alg.k_factors)
    cost.check_memory(cost.tensor_bytes((B, J, alg.C, alg.k_max, alg.k_max), dtype), dev,
                      f"spectral_measure ({J} outcomes of {alg.C} sectors up to {alg.k_max})")
    if not group:
        E, vals, outs = [], [], []
        for c, (w, V) in enumerate(dec):
            for a in range(alg.k_factors[c]):
                e = torch.zeros(B, alg.C, alg.k_max, alg.k_max, dtype=dtype, device=dev)
                v = V[:, :, a]
                e[:, c, :alg.k_factors[c], :alg.k_factors[c]] = v.unsqueeze(-1) * v.conj().unsqueeze(-2)
                E.append(e)
                vals.append(w[:, a])
                outs.append((c, a))
        return PVM(torch.stack(E, dim=1), torch.stack(vals, dim=-1), alg, check=False, outcomes=outs)
    if B != 1:
        raise ValueError("spectral_measure: group=True needs a single operator (batch 1); use group=False for batches")
    items = sorted(((float(w[0, a]), c, a) for c, (w, V) in enumerate(dec) for a in range(alg.k_factors[c])))
    scale = max(abs(x[0]) for x in items) or 1.0
    tol = (_default_tol(dtype) ** 0.5) * scale if tol is None else tol    # relative to ||A||
    groups = []
    for lam, c, a in items:                              # compare with the first member: no chaining of small gaps
        if groups and lam - groups[-1][0][0] <= tol:
            groups[-1].append((lam, c, a))
        else:
            groups.append([(lam, c, a)])
    E = torch.zeros(1, len(groups), alg.C, alg.k_max, alg.k_max, dtype=dtype, device=dev)
    vals = []
    for j, grp in enumerate(groups):
        for lam, c, a in grp:
            k_c = alg.k_factors[c]
            v = dec[c][1][0, :, a]
            E[0, j, c, :k_c, :k_c] += v.unsqueeze(-1) * v.conj().unsqueeze(-2)
        vals.append(torch.stack([dec[c][0][0, a] for _, c, a in grp]).mean())     # differentiable in A
    values = torch.stack(vals).unsqueeze(0).to(alg.hilbert.real_dtype)
    return PVM(E, values, alg, check=False)
