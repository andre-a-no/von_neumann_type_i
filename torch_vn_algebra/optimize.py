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
    StinespringParam        channels: an isometry V = Z (Z^* Z)^{-1/2} per input sector (Haar-random at the start),
                            Kraus operators its blocks; within one algebra (Channel) or between algebras with a
                            prescribed set of allowed sector transitions (InterSectorChannel)
    LindbladParam           generators: a self-adjoint H and free jump operators (inside the sectors and, with
                            allowed transitions, between them), i.e. the Kossakowski matrix as L^* L; for
                            dynamics.lindblad_evolve

    best = optimize.extremize(objective, params, maximize=True)     # best value per start
    res = optimize.fit(loss, params)                                # minimisation; the parameters end at the fit

Initialisation (`init`, for a fraction `init_fraction` of the starts; the default 'haar' is unchanged): the extremes
of many trace functionals sit at aligned eigenbases, a set of Haar measure zero, and in a real algebra a unitary chart
U0 exp(G - G^T) never leaves the connected component (det = +-1 per sector) of U0. `structured_unitary` provides
'diagonal' (random signs on the diagonal), 'pairs' (2x2 rotations of consecutive basis vectors at random angles) and
'oriented' (Haar, with det U0_c = +1 in every sector of a real algebra) starting points; `spectrum_init` gives the
structured starts two-level spectra. With X on 'diagonal' and Y on 'pairs', the relative eigenbasis of X and Y is
a direct sum of 2x2 rotations, the structure of the maximisers of Exp. 1 and Exp. 3 (paper, Section 6), in every
dimension.

The objective maps the current operators to one value per start, shape (batch,). The search is
local: the result is a lower bound for a supremum (an upper bound for an infimum) and should be
checked with more starts / steps.
"""
from typing import Callable, List, Optional, Sequence

import torch

from .algebra import TypeIAlgebra
from .channels import Channel, InterSectorChannel
from . import cost

Operator = TypeIAlgebra.Operator


INITS = ('haar', 'diagonal', 'pairs', 'rotations', 'oriented')


def structured_unitary(alg: TypeIAlgebra, batch_size: int, kind: str = 'pairs') -> torch.Tensor:
    """
    Block-diagonal unitaries (batch, C, k_max, k_max) (zero on the padding) as starting points of a search:
      'haar'      Haar-random in every sector (random_unitary_operator);
      'diagonal'  diag(+-1) with independent random signs;
      'pairs'     R diag(+-1): R rotates the consecutive basis vectors (e_0, e_1), (e_2, e_3), ... of every sector
                  by independent angles uniform in [0, pi) (an odd sector keeps its last vector), with random signs;
      'rotations' R alone (det = +1 in every sector);
      'oriented'  Haar-random, and in a real algebra with det = +1 in every sector (first column negated where
                  det = -1); in a complex algebra the same as 'haar'.
    """
    if kind not in INITS:
        raise ValueError(f"unknown init {kind!r}; choose from {INITS}")
    base = alg.random_unitary_operator(batch_size).matrix
    if kind == 'haar':
        return base
    if kind == 'oriented':
        if alg.complex_valued:
            return base
        k = torch.tensor(alg.k_factors, device=base.device)
        pad = (torch.arange(alg.k_max, device=base.device) >= k[:, None]).to(base.dtype)       # (C, k_max)
        det = torch.linalg.det(base + torch.diag_embed(pad))                                   # padding as 1
        col = torch.ones(batch_size, alg.C, alg.k_max, dtype=base.dtype, device=base.device)
        col[..., 0] = torch.where(det < 0, -1.0, 1.0).to(base.dtype)                          # negate column 0
        return base * col[..., None, :]
    rdt, dev, B, C, n = alg.hilbert.real_dtype, base.device, batch_size, alg.C, alg.k_max
    k = torch.tensor(alg.k_factors, device=dev)
    act = (torch.arange(n, device=dev) < k[:, None]).to(rdt)                                     # (C, n)
    signs = torch.where(torch.rand(B, C, n, device=dev) < 0.5, -1.0, 1.0).to(rdt) * act
    if kind == 'diagonal':
        return torch.diag_embed(signs).to(base.dtype)
    R = torch.diag_embed(act.expand(B, C, n)).clone()
    t = torch.pi * torch.rand(B, C, n // 2, dtype=rdt, device=dev)
    i = torch.arange(n // 2, device=dev)
    full = (2 * i + 1 < k[:, None]).to(rdt)                                                     # (C, n//2): pair active
    c, s = torch.cos(t) * full + (1 - full), torch.sin(t) * full
    R[..., 2 * i, 2 * i], R[..., 2 * i + 1, 2 * i + 1] = c * act[:, 2 * i], c * act[:, 2 * i + 1]
    R[..., 2 * i, 2 * i + 1], R[..., 2 * i + 1, 2 * i] = -s, s
    if kind == 'rotations':
        return R.to(base.dtype)
    return (R * signs[..., None, :]).to(base.dtype)                                             # R diag(+-1)


def _initial_base(alg: TypeIAlgebra, batch_size: int, init: str, init_fraction: float):
    """Haar bases with some starts replaced by structured_unitary(init): a random fraction `init_fraction` (float) of
    them, or those selected by a boolean tensor (batch,) (e.g. the `structured` mask of another parameter, so that X
    and Y start structured at the same starts). Returns (base, mask)."""
    base = alg.random_unitary_operator(batch_size).matrix
    mask = torch.zeros(batch_size, dtype=torch.bool, device=base.device)
    if init != 'haar':
        if torch.is_tensor(init_fraction):
            mask = init_fraction.to(device=base.device, dtype=torch.bool).expand(batch_size).clone()
        elif init_fraction > 0:
            mask = torch.rand(batch_size, device=base.device) < init_fraction
        base = torch.where(mask[:, None, None, None], structured_unitary(alg, batch_size, init), base)
    return base, mask


class Param:
    """Base class: a differentiable parametrisation of a batch of operators."""

    def parameters(self) -> List[torch.Tensor]:
        raise NotImplementedError

    def recenter(self):
        """Move the reference point to the current point (default: nothing)."""


class UnitaryParam(Param):
    """Block-diagonal unitaries U = U0 exp(G - G^*), one per start; U0 Haar-random by default, or a structured start
    (`init`, see structured_unitary) for a random fraction `init_fraction` of the starts."""

    def __init__(self, alg: TypeIAlgebra, batch_size: int, base: Optional[torch.Tensor] = None, init: str = 'haar',
                 init_fraction=1.0):
        self.alg = alg
        if init not in INITS:
            raise ValueError(f"unknown init {init!r}; choose from {INITS}")
        self.base = _initial_base(alg, batch_size, init, init_fraction)[0] if base is None else base
        self.G = torch.zeros_like(self.base, requires_grad=True)
        self._cache = (None, None)
        self.G.register_hook(self._invalidate)       # a backward pass frees the graph of the cached U

    def _invalidate(self, grad):
        self._cache = (None, None)
        return None

    def parameters(self):
        return [self.G]

    def matrix(self) -> torch.Tensor:
        """U as a tensor; computed once per parameter state (an objective that uses U several times, e.g. X and
        |X| of a PositiveParam, pays one matrix exponential)."""
        key = (self.G._version, id(self.base), torch.is_grad_enabled())
        if self._cache[0] != key:
            skew = self.G - self.G.conj().transpose(-2, -1)
            self._cache = (key, self.base @ self.alg.operator(skew).expm().matrix)
        return self._cache[1]

    def operator(self) -> Operator:
        return self.alg.operator(self.matrix(), is_normal=None)

    def recenter(self):
        with torch.no_grad():
            self.base = self.matrix().detach()
            self.G.zero_()
        self._cache = (None, None)


def _contrast_to_lo(delta: torch.Tensor) -> torch.Tensor:
    return (1 - delta) / (1 + delta)


class SpectrumParam(Param):
    """
    Spectra (per block, padded) in [lo, 1] with the global maximum 1 and minimum lo, so that the
    Michelson contrast equals delta = (1 - lo) / (1 + lo).

    delta: tensor (batch,) of prescribed contrasts, or None to optimise the contrast as well
    within [delta_range[0], delta_range[1]]. The maximum 1 and the minimum lo are pinned to two
    eigenvalue slots; the other eigenvalues lie in [lo, 1] (through (1 + sin theta) / 2). By default
    the two slots are drawn at random per start among all active slots (the extremes may lie in the
    same sector or in different ones), so that a multistart search covers every placement and the
    union over placements is dense in the constraint set. With `pin_sector` both extremes sit in that
    sector (which needs k >= 2) for every start; for C >= 2 this is a restriction.
    """

    def __init__(self, alg: TypeIAlgebra, batch_size: int, delta: Optional[torch.Tensor] = None,
                 delta_range=(0.0, 0.999), pin_sector: Optional[int] = None):
        self.alg = alg
        rdt, dev = alg.hilbert.real_dtype, alg.hilbert.device
        slots = [(c, i) for c, k_c in enumerate(alg.k_factors) for i in range(k_c)]   # active (sector, index)
        if pin_sector is not None:
            if alg.k_factors[pin_sector] < 2:
                raise ValueError("the pinned sector needs at least two dimensions")
            i_max = torch.full((batch_size,), slots.index((pin_sector, 0)), dtype=torch.long)
            i_min = i_max + 1
        else:
            if len(slots) < 2:
                raise ValueError("a prescribed contrast needs at least two dimensions")
            i_max = torch.randint(len(slots), (batch_size,))                         # two distinct slots
            i_min = (i_max + 1 + torch.randint(len(slots) - 1, (batch_size,))) % len(slots)
        sl = torch.tensor(slots, dtype=torch.long)
        self.pin_max = sl[i_max].to(dev)          # (batch, 2): sector and index of lambda_max
        self.pin_min = sl[i_min].to(dev)          # (batch, 2): sector and index of lambda_min
        self.theta = torch.randn(batch_size, alg.C, alg.k_max, dtype=rdt, device=dev).requires_grad_(True)
        self.fixed_delta = None if delta is None else torch.as_tensor(delta, dtype=rdt, device=dev).expand(batch_size)
        self.delta_range = delta_range
        self.eta = None if delta is not None else torch.randn(batch_size, dtype=rdt, device=dev).requires_grad_(True)
        self.mask = alg._const('rows')
        b = torch.arange(batch_size, device=dev)
        is_max = torch.zeros(batch_size, alg.C, alg.k_max, dtype=torch.bool, device=dev)
        is_min = torch.zeros_like(is_max)
        is_max[b, self.pin_max[:, 0], self.pin_max[:, 1]] = True
        is_min[b, self.pin_min[:, 0], self.pin_min[:, 1]] = True
        self._is_max, self._is_min = is_max, is_min               # constant: built once

    def parameters(self):
        return [self.theta] + ([self.eta] if self.eta is not None else [])

    def delta(self) -> torch.Tensor:
        if self.fixed_delta is not None:
            return self.fixed_delta
        a, b = self.delta_range
        return a + (b - a) * 0.5 * (1 + torch.sin(self.eta))      # reaches both ends (e.g. Y = 1 at contrast 0)

    def eigenvalues(self) -> torch.Tensor:
        """(batch, C, k_max), zero on padding."""
        lo = _contrast_to_lo(self.delta())[:, None, None]
        # (1 + sin)/2 reaches the ends of [lo, 1] at finite theta (a sigmoid only asymptotically, which
        # slows the search down: extremal spectra often have all eigenvalues at lo or 1)
        lam = lo + (1 - lo) * 0.5 * (1 + torch.sin(self.theta))
        lam = torch.where(self._is_max, torch.ones_like(lam), torch.where(self._is_min, lo.expand_as(lam), lam))
        return lam * self.mask


class PositiveParam(Param):
    """Positive operators X = U diag(lambda) U^* with a prescribed or bounded Michelson contrast. The eigenbasis
    U = U0 exp(G - G^*) starts at a Haar-random U0, or at structured_unitary(init) for a random fraction
    `init_fraction` of the starts (`self.structured` marks them)."""

    def __init__(self, alg: TypeIAlgebra, batch_size: int, delta: Optional[torch.Tensor] = None,
                 delta_range=(0.0, 0.999), pin_sector: Optional[int] = None, init: str = 'haar',
                 init_fraction=1.0, spectrum_init: Optional[str] = None):
        self.alg = alg
        self.spectrum = SpectrumParam(alg, batch_size, delta, delta_range, pin_sector)
        if init not in INITS:
            raise ValueError(f"unknown init {init!r}; choose from {INITS}")
        base, self.structured = _initial_base(alg, batch_size, init, init_fraction)
        self.unitary = UnitaryParam(alg, batch_size, base=base)
        if spectrum_init is not None:
            self._init_spectrum(spectrum_init)

    def _init_spectrum(self, kind: str, t: float = 1.5):
        """Two-level starting spectra for the structured starts: 'pairs' puts one eigenvalue near lambda_max and one
        near lambda_min in every pair of consecutive slots (random order), 'high' all near lambda_max; a free
        contrast starts near the upper end of its range. (sin(+-t) = +-0.9975: near the ends, where the gradient of
        (1 + sin)/2 does not vanish.)"""
        if kind not in ('pairs', 'high'):
            raise ValueError("spectrum_init must be 'pairs', 'high' or None")
        sp, m = self.spectrum, self.structured
        B, C, n = sp.theta.shape
        with torch.no_grad():
            if kind == 'high':
                th = torch.full_like(sp.theta, t)
            else:
                first = torch.rand(B, C, (n + 1) // 2, device=sp.theta.device) < 0.5
                hi = torch.stack([first, ~first], dim=-1).reshape(B, C, -1)[..., :n]
                # consistent with the pinned slots: lambda_max on a high slot and lambda_min on a low one
                b = torch.arange(B, device=hi.device)
                k = torch.tensor(self.alg.k_factors, device=hi.device)
                for pin, val in ((sp.pin_max, True), (sp.pin_min, False)):
                    c, i = pin.unbind(-1)
                    j = i ^ 1                                                   # partner slot in the pair
                    ok = j < k[c]
                    hi[b[ok], c[ok], j[ok]] = not val
                for pin, val in ((sp.pin_max, True), (sp.pin_min, False)):
                    hi[b, pin[:, 0], pin[:, 1]] = val
                th = torch.where(hi, torch.full_like(sp.theta, t), torch.full_like(sp.theta, -t))
            sp.theta.copy_(torch.where(m[:, None, None], th, sp.theta))
            if sp.eta is not None:
                sp.eta.copy_(torch.where(m, torch.full_like(sp.eta, t), sp.eta))

    def parameters(self):
        return self.spectrum.parameters() + self.unitary.parameters()

    def recenter(self):
        self.unitary.recenter()

    def _signed(self, lam):
        return lam

    def _build(self, lam: torch.Tensor) -> Operator:
        U = self.unitary.matrix()
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
    signs s = +-1 are drawn at random per start (discrete variables; the multistart covers the sign patterns, and
    discrete_step flips single signs between rounds). For the starts with a structured eigenbasis (init != 'haar'),
    the signs alternate within the consecutive pairs of slots (one + and one - per pair, in random order), the sign
    pattern of the maximisers of Exp. 3 (paper, Section 6).
    """

    def __init__(self, alg: TypeIAlgebra, batch_size: int, delta: Optional[torch.Tensor] = None,
                 delta_range=(0.0, 0.999), pin_sector: Optional[int] = None, init: str = 'haar',
                 init_fraction=1.0, spectrum_init: Optional[str] = None):
        super().__init__(alg, batch_size, delta, delta_range, pin_sector, init, init_fraction, spectrum_init)
        dev = alg.hilbert.device
        # U absorbs permutations within a sector, so the sign pattern matters only through the number n_c of
        # negative eigenvalues in every sector. Independent signs would almost never give a pattern that is
        # the same in all C sectors; n_c is therefore drawn uniformly in {0, ..., k_c}, for half of the starts
        # independently per sector and for the other half as one fraction common to all sectors.
        B, C = batch_size, alg.C
        k = torch.tensor(alg.k_factors, device=dev, dtype=torch.float64)
        frac = torch.rand(B, C, device=dev, dtype=torch.float64)
        common = torch.rand(B, device=dev) < 0.5
        frac = torch.where(common[:, None], frac[:, :1].expand(B, C), frac)
        n_neg = torch.floor(frac * (k + 1)).clamp(max=k)                              # (B, C)
        keys = torch.rand(B, C, alg.k_max, device=dev, dtype=torch.float64)
        keys = keys + (torch.arange(alg.k_max, device=dev) >= k[:, None]).to(keys.dtype)   # padding ranks last
        rank = torch.argsort(torch.argsort(keys, dim=-1), dim=-1)                   # random order within a sector
        neg = rank < n_neg[..., None]
        if self.structured.any():                                       # alternating within consecutive pairs
            first = torch.rand(B, C, (alg.k_max + 1) // 2, device=dev) < 0.5
            alt = torch.stack([first, ~first], dim=-1).reshape(B, C, -1)[..., :alg.k_max]
            odd = torch.arange(alg.k_max, device=dev) == (k[:, None] - 1).long()          # last slot of odd k
            odd = odd & (k[:, None].long() % 2 == 1)
            alt = torch.where(odd, torch.rand(B, C, alg.k_max, device=dev) < 0.5, alt)
            neg = torch.where(self.structured[:, None, None], alt, neg)
        self.signs = torch.where(neg, -1.0, 1.0).to(alg.hilbert.real_dtype)

    def _signed(self, lam):
        return lam * self.signs

    def discrete_step(self, evaluate: Callable[[], torch.Tensor], maximize: bool, trials: int = 8):
        """Local search over the signs: flip the sign of one random active eigenvalue per start and keep the
        flip where the objective improves (between the gradient rounds of extremize)."""
        alg, B = self.alg, self.signs.shape[0]
        slots = alg._const_cache.get('active_slots')
        if slots is None:
            slots = alg._const_cache['active_slots'] = torch.tensor(
                [(c, i) for c, k_c in enumerate(alg.k_factors) for i in range(k_c)], device=self.signs.device)
        b = torch.arange(B, device=self.signs.device)
        with torch.no_grad():
            best = evaluate()
            for _ in range(trials):
                c, i = slots[torch.randint(len(slots), (B,), device=slots.device)].unbind(-1)
                self.signs[b, c, i] *= -1
                v = evaluate()
                better = v > best if maximize else v < best
                # undo where it did not help (no data-dependent indexing: no synchronisation)
                self.signs[b, c, i] *= torch.where(better, 1.0, -1.0).to(self.signs.dtype)
                best = torch.where(better, v, best)


class LindbladParam(Param):
    """
    GKSL generators of `alg` (one per start) for dynamics.lindblad_evolve: a self-adjoint Hamiltonian
    H = (X + X^*)/2 and free jump operators, n_jumps inside the sectors (Operators) and n_inter between sectors
    (InterSectorChannels of alg into itself, nonzero only for the allowed transitions[d, c], e.g. N -> N - 1).
    Free jump operators L_j are a Cholesky-type factor of the Kossakowski matrix, so every parameter value is a
    valid generator; the rates are absorbed into the L_j, and unitary mixing of the jumps (L_j -> sum_l u_jl L_l)
    is a gauge symmetry. The jumps start at the scale `init_scale` (weak dissipation).

        P = optimize.LindbladParam(alg, starts, n_jumps=2, n_inter=1, transitions=lowering)
        rho = dynamics.lindblad_evolve(rho0, P.hamiltonian(), P.jumps(), times, method='exact')
    """

    def __init__(self, alg: TypeIAlgebra, batch_size: int, n_jumps: int = 1, n_inter: int = 0, transitions=None,
                 hamiltonian: bool = True, init_scale: float = 0.1, generator: Optional[torch.Generator] = None):
        self.alg = alg
        dev = alg.hilbert.device
        C, k = alg.C, alg.k_max
        mask = alg._const('block_mask', torch.bool).to(dev)

        def free(shape, scale):
            Z = scale * torch.randn(shape, dtype=alg.dtype, device=dev, generator=generator)
            return Z.requires_grad_(True)
        self.X = free((batch_size, C, k, k), 1.0) if hamiltonian else None
        self.L = free((batch_size, n_jumps, C, k, k), init_scale) if n_jumps else None
        self._mask = mask
        self.inter_mask = None
        self.K = None
        if n_inter:
            allowed = (torch.ones(C, C, dtype=torch.bool, device=dev) if transitions is None
                       else torch.as_tensor(transitions, dtype=torch.bool, device=dev).reshape(C, C))
            rows = alg._const('rows').to(dev)                                       # (C, k)
            self.inter_mask = allowed[:, :, None, None] & rows[:, None, :, None] & rows[None, :, None, :]
            self.K = free((batch_size, n_inter, C, C, k, k), init_scale)

    def parameters(self):
        return [t for t in (self.X, self.L, self.K) if t is not None]

    def hamiltonian(self) -> Optional[Operator]:
        if self.X is None:
            return None
        X = self.X * self._mask
        return self.alg.operator(0.5 * (X + X.conj().transpose(-2, -1)), is_self_adjoint=True)

    def jumps(self) -> list:
        out = []
        if self.L is not None:
            L = self.L * self._mask
            out += [self.alg.operator(L[:, j]) for j in range(L.shape[1])]
        if self.K is not None:
            K = self.K * self.inter_mask
            out += [InterSectorChannel(self.alg, self.alg, K[:, j:j + 1], trace_preserving=False)
                    for j in range(K.shape[1])]
        return out


class StinespringParam(Param):
    """
    Trace-preserving CP maps from the algebra `alg` (one per start): for every input sector c an isometry
    V_c = Z_c (Z_c^* Z_c)^{-1/2} from C^{k_c} into the environment C^r times the output sectors, and the Kraus
    operators are the blocks of V_c. Every value of the free complex matrices Z is a channel (no constraint, no
    projection); Gaussian Z gives Haar-random isometries, i.e. the distribution of channels.random_channel and
    channels.random_inter_sector_channel. Gradients go through Z^* Z by the Daleckii-Krein formula. The Kraus gauge
    (V -> (W (x) 1) V with a unitary W of the environment) is a symmetry of every objective of the channel.

    alg_out=None: sector-preserving channels of `alg` (channel() is a Channel, Kraus rank r). Otherwise
    InterSectorChannels into alg_out; `transitions` (an (D, C) boolean matrix, default all True) fixes which output
    sectors d an input sector c may reach, e.g. a charge-lowering map. The isometry needs r * sum_{d allowed} m_d
    >= k_c for every input sector c.
    """

    def __init__(self, alg: TypeIAlgebra, batch_size: int, kraus_rank: int = 1, alg_out: Optional[TypeIAlgebra] = None,
                 transitions=None, generator: Optional[torch.Generator] = None):
        self.alg, self.alg_out, self.r = alg, alg_out, kraus_rank
        dev = alg.hilbert.device
        C, k = alg.C, alg.k_max
        if alg_out is None:
            if transitions is not None:
                raise ValueError("transitions need alg_out (a map between algebras)")
            D, m = C, k
            allowed = torch.eye(C, dtype=torch.bool, device=dev)
        else:
            D, m = alg_out.C, alg_out.k_max
            allowed = (torch.ones(D, C, dtype=torch.bool, device=dev) if transitions is None
                       else torch.as_tensor(transitions, dtype=torch.bool, device=dev).reshape(D, C))
        out_rows = (alg_out or alg)._const('rows').to(dev)                        # (D, m) active rows of the outputs
        width = (kraus_rank * (allowed.to(torch.long) * out_rows.sum(-1)[:, None]).sum(0)).tolist()
        for c, k_c in enumerate(alg.k_factors):
            if k_c and width[c] < k_c:
                raise ValueError(f"input sector {c} (dimension {k_c}) cannot be mapped isometrically into "
                                 f"{width[c]} dimensions: raise kraus_rank or allow more transitions")
        # mask of Z, shape (C, r, D, m, k): allowed (d, c), active rows of d, active columns of c
        in_cols = alg._const('rows').to(dev)                                       # (C, k)
        self._mask = (allowed.t()[:, None, :, None, None] & out_rows[None, None, :, :, None]
                      & in_cols[:, None, None, None, :]).expand(C, kraus_rank, D, m, k)
        shape = (batch_size, C, kraus_rank, D, m, k)
        Z = torch.randn(shape, dtype=alg.dtype, device=dev, generator=generator)
        self.Z = (Z * self._mask).requires_grad_(True)
        self._D, self._m = D, m

    def parameters(self):
        return [self.Z]

    def kraus(self) -> torch.Tensor:
        """Kraus tensor (batch, r, D, C, m, k) (between algebras) or (batch, r, C, k, k) (within one algebra)."""
        B, C, r, D, m, k = self.Z.shape
        Z = (self.Z * self._mask).reshape(B, C, r * D * m, k)
        gram = self.alg.operator(Z.conj().transpose(-2, -1) @ Z, is_self_adjoint=True, is_positive=True)
        V = Z @ gram.power(-0.5).matrix                                            # (B, C, r D m, k), V_c^* V_c = 1
        K = V.reshape(B, C, r, D, m, k).permute(0, 2, 3, 1, 4, 5)                  # (B, r, D, C, m, k)
        if self.alg_out is None:
            idx = torch.arange(C, device=K.device)
            K = K[:, :, idx, idx]                                                  # (B, r, C, k, k): d = c
        return K

    def channel(self):
        """The current channels (trace preserving by construction)."""
        if self.alg_out is None:
            return Channel(self.alg, self.kraus(), trace_preserving=True)
        return InterSectorChannel(self.alg, self.alg_out, self.kraus(), trace_preserving=True)


def extremize(objective: Callable[[], torch.Tensor], params: Sequence[Param], maximize: bool = True,
              rounds: int = 10, steps: int = 50, lr: float = 0.05, decay: float = 0.6,
              progress: bool = False) -> torch.Tensor:
    """
    Gradient ascent (maximize=True) or descent on objective() -> tensor (batch,), every batch
    element being an independent start. Adam with learning rate lr * decay**round; unitary charts
    are re-centred after every round. Returns, per start, the best objective value attained at any
    iterate (every iterate is an admissible point; no gradient).
    """
    sign = -1.0 if maximize else 1.0
    timer = cost.StepTimer(rounds * steps, "optimize.extremize", progress)
    best = None

    def keep(v):
        nonlocal best
        v = torch.nan_to_num(v.detach(), nan=-float('inf') if maximize else float('inf'))
        best = v if best is None else (torch.maximum(best, v) if maximize else torch.minimum(best, v))

    for rnd in range(rounds):
        tensors = [t for p in params for t in p.parameters()]
        opt = torch.optim.Adam(tensors, lr=lr * decay ** rnd)
        for _ in range(steps):
            v = objective()
            keep(v)
            loss = sign * v.sum()
            opt.zero_grad()
            loss.backward()
            opt.step()
            timer.step(tensors[0].device)
        for p in params:
            p.recenter()
        if rnd < rounds - 1:
            for p in params:
                if hasattr(p, 'discrete_step'):
                    p.discrete_step(objective, maximize)
    timer.close()
    with torch.no_grad():
        keep(objective())
    # a start whose objective was never finite has no attained value: NaN rather than -inf / +inf
    return torch.where(torch.isinf(best), torch.full_like(best, float('nan')), best)


class FitResult:
    """Outcome of fit: `loss` per start at the final point, `best` the index of the smallest, `history` the mean
    loss over the starts after every round. The parameters themselves are left at the final point."""

    def __init__(self, loss: torch.Tensor, history: List[float]):
        self.loss, self.history = loss, history
        finite = torch.nan_to_num(loss, nan=float('inf'))
        self.best = int(finite.argmin()) if finite.numel() else 0

    def __repr__(self) -> str:
        return f"FitResult(best loss {float(self.loss[self.best]):.3e} at start {self.best}, {self.loss.numel()} starts)"


def fit(loss: Callable[[], torch.Tensor], params: Sequence[Param], rounds: int = 10, steps: int = 50,
        lr: float = 0.05, decay: float = 0.6, progress: bool = False) -> FitResult:
    """
    Minimise loss() -> tensor (batch,) over the parameters, every batch element an independent start (the same
    schedule as extremize: Adam with learning rate lr * decay**round, unitary charts re-centred after every round).
    Unlike extremize, which reports the best value seen, fit leaves the parameters at the final point, so that the
    fitted operators or channels can be read from them (e.g. StinespringParam.channel()), and returns the final loss
    per start. A local search: compare several starts, and convex problems with an exact solver.
    """
    timer = cost.StepTimer(rounds * steps, "optimize.fit", progress)
    history = []
    for rnd in range(rounds):
        tensors = [t for p in params for t in p.parameters()]
        opt = torch.optim.Adam(tensors, lr=lr * decay ** rnd)
        for _ in range(steps):
            v = loss()
            opt.zero_grad()
            torch.nan_to_num(v, nan=0.0).sum().backward()
            opt.step()
            timer.step(tensors[0].device)
        for p in params:
            p.recenter()
        history.append(float(v.detach().mean()))
    timer.close()
    with torch.no_grad():
        final = loss().detach()
    return FitResult(final, history)
