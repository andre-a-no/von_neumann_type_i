"""
Quantum channels as trainable networks: layers that are channels of finite-dimensional von Neumann algebras, as
torch.nn.Modules, composable with torch.nn.Sequential and trained with any torch.optim optimiser.

Every layer maps a DensityMatrix (a batch of states) to a DensityMatrix, and every value of its parameters is a
physical map, so training never leaves the set of channels:

    ChannelLayer      a channel from Stinespring isometries (optimize.StinespringParam): within one algebra, or
                      between algebras with a prescribed set of allowed sector transitions (width = Kraus rank)
    UnitaryLayer      rho -> U rho U^*, U = exp(-i H) with a free self-adjoint H in the algebra
    LindbladLayer     rho -> exp(t L)(rho) for a free GKSL generator (optimize.LindbladParam): continuous depth
    Residual          rho -> (1 - p) rho + p Phi(rho): the skip connection that keeps complete positivity; a stack
                      of such blocks with small p is a time step of a jump process with generator Phi - id
    PostSelect        rho -> P rho P / Tr(P rho P) for a central projection P: the only nonlinear layer, at the cost
                      of the acceptance Tr(P rho P), which it records (use log(acceptance) as a regulariser)
    FeedForward       measure an instrument, record the outcome j in the centre, apply the branch j:
                      rho -> sum_j Phi'_j(Phi_j(rho)); classical information is the only thing that may be copied
                      forward, so this is the quantum form of a skip that fans out

    net = torch.nn.Sequential(learn.ChannelLayer(alg, kraus_rank=2), learn.Residual(learn.LindbladLayer(alg), p=0.1))
    out = net(rho)                                             # DensityMatrix
    loss = learn.infidelity(out, target).mean()
    loss.backward(); opt.step()

Losses: infidelity, trace_distance, relative_entropy (states), negative_log_likelihood (measurement.POVM data).
count_parameters gives the number of real parameters of a network (sector structure shows up here: a channel that
respects the sectors has far fewer than a dense one).

Channel, unitary, Lindblad, post-selection and feed-forward layers support Module.to(device=..., dtype=...) for
single/double real or complex arithmetic. Their algebra metadata follows the registered tensors; the original
shared algebra is unchanged. Convert inputs with rho.cast(layer.algebra) after a transfer. Assigned checkpoints
also update this metadata. FeedForward checkpoints include the fixed instrument and its outcome labels.
PostSelect returns a conditional state only for acceptance above eps; other samples are zero operators.
"""
import inspect
import math
from typing import Optional, Sequence, Union

import torch
from torch import nn

from .algebra import TypeIAlgebra
from .channels import Channel, InterSectorChannel
from .states import DensityMatrix
from . import dynamics, measurement, optimize, states

Operator = TypeIAlgebra.Operator
_APPLY_HAS_RECURSE = 'recurse' in inspect.signature(nn.Module._apply).parameters


def _as_state(alg: TypeIAlgebra, mat: torch.Tensor) -> DensityMatrix:
    return DensityMatrix(alg, matrix=mat, validate=False)


def _sync_loaded_layer(module, incompatible_keys):
    # assign=True can change dtype/device as well as replace a Parameter or buffer.
    module._sync_algebra()


def _algebra_like_tensor(alg: TypeIAlgebra, tensor: torch.Tensor) -> TypeIAlgebra:
    if tensor.dtype not in (torch.float32, torch.float64, torch.complex64, torch.complex128):
        raise TypeError("algebra layers support float32, float64, complex64 and complex128")
    if alg.dtype == tensor.dtype and torch.device(alg.hilbert.device) == tensor.device:
        return alg
    return alg.like(complex_valued=torch.is_complex(tensor),
                    precision='double' if tensor.dtype in (torch.float64, torch.complex128) else 'single',
                    device=tensor.device)


class _AlgebraLayer(nn.Module):
    """Keep non-Module algebra metadata in step with the registered tensors."""

    def __init__(self, algebra: TypeIAlgebra):
        super().__init__()
        self.algebra = algebra
        # Also gives a parameter-free Lindblad identity layer a dtype/device to transfer.
        self.register_buffer('_algebra_anchor', torch.empty(0, dtype=algebra.dtype, device=algebra.hilbert.device),
                             persistent=False)
        self.register_load_state_dict_post_hook(_sync_loaded_layer)

    def _apply(self, fn, recurse=True):
        # Older supported PyTorch versions do not expose the recurse argument.
        super()._apply(fn, **({'recurse': recurse} if _APPLY_HAS_RECURSE else {}))
        self._sync_algebra()
        return self

    def _sync_algebra(self):
        tensors = list(self.parameters(recurse=False))
        tensors += [t for name, t in self._buffers.items() if name != '_algebra_anchor' and t is not None]
        tensor = tensors[0] if tensors else self._algebra_anchor
        old = self.algebra
        self.algebra = _algebra_like_tensor(old, tensor)
        if hasattr(self, 'algebra_out'):
            self.algebra_out = (self.algebra if self.algebra_out is old else
                                _algebra_like_tensor(self.algebra_out, tensor))
        self._algebra_anchor = tensor.new_empty(0)
        self._sync_parameterisation(tensor.device)

    def _sync_parameterisation(self, device):
        pass


# ----------------------------------------------------------------------
# layers
# ----------------------------------------------------------------------
class ChannelLayer(_AlgebraLayer):
    """
    A trainable channel: Stinespring isometries V_c = Z_c (Z_c^* Z_c)^{-1/2}, one per input sector, Kraus operators
    their blocks (optimize.StinespringParam with one start). algebra_out=None: sector-preserving channels of
    `algebra`; otherwise maps into algebra_out with the allowed sector transitions `transitions` ((D, C) booleans).
    Gaussian initial Z gives a Haar-random channel; init_scale < 1 is not used here, see Residual for near-identity
    blocks. The polar factors are evaluated directly from Z (without a Gram matrix); deficient rank is completed
    to an isometry, while gradients require numerical full column rank.
    """

    def __init__(self, algebra: TypeIAlgebra, kraus_rank: int = 1, algebra_out: Optional[TypeIAlgebra] = None,
                 transitions=None, generator: Optional[torch.Generator] = None):
        super().__init__(algebra)
        self._param = optimize.StinespringParam(algebra, 1, kraus_rank, alg_out=algebra_out, transitions=transitions,
                                                generator=generator)
        self.Z = nn.Parameter(self._param.Z.detach().clone())
        self._param.Z = self.Z
        self.algebra, self.algebra_out = algebra, algebra_out or algebra

    def _sync_parameterisation(self, device):
        self._param.alg = self.algebra
        if self._param.alg_out is not None:
            self._param.alg_out = self.algebra_out
        self._param.Z = self.Z
        self._param._mask = self._param._mask.to(device=device)
        self._param._polar_groups = [(k, sectors.to(device=device), rows.to(device=device))
                                     for k, sectors, rows in self._param._polar_groups]

    def channel(self) -> Union[Channel, InterSectorChannel]:
        # Loading with assign=True can replace the registered Parameter.
        self._param.Z = self.Z
        return self._param.channel()

    def forward(self, rho: DensityMatrix) -> DensityMatrix:
        return self.channel()(rho)

    def extra_repr(self) -> str:
        return f"kraus_rank={self._param.r}, sectors {self.algebra.k_factors} -> {self.algebra_out.k_factors}"


class UnitaryLayer(_AlgebraLayer):
    """rho -> U rho U^* with U = exp(-i H), H = (X + X^*)/2 free in the algebra (a unitary of prod_c U(k_c))."""

    def __init__(self, algebra: TypeIAlgebra, init_scale: float = 1.0, generator: Optional[torch.Generator] = None):
        super().__init__(algebra)
        X = init_scale * torch.randn(1, algebra.C, algebra.k_max, algebra.k_max, dtype=algebra.dtype,
                                     device=algebra.hilbert.device, generator=generator)
        self.X = nn.Parameter(X)

    def unitary(self) -> Operator:
        X = self.X * self.algebra._const('block_mask', torch.bool)
        H = self.algebra.operator(0.5 * (X + X.conj().transpose(-2, -1)), is_self_adjoint=True)
        return H.expm(-1j)

    def forward(self, rho: DensityMatrix) -> DensityMatrix:
        U = self.unitary().matrix
        return _as_state(self.algebra, U @ rho.matrix @ U.conj().transpose(-2, -1))


class LindbladLayer(_AlgebraLayer):
    """
    rho -> exp(time L)(rho) for a trainable GKSL generator L (optimize.LindbladParam): a self-adjoint Hamiltonian,
    n_jumps free jump operators inside the sectors and n_inter between sectors (allowed transitions[d, c]). The
    evolution is computed by dynamics.lindblad_evolve (method 'exact' by default, or the CPTP splitting 'kraus');
    every value of the parameters gives a channel, and the layer is a block of continuous depth.
    Other integrators, including 'rk4', are rejected: a Runge-Kutta step need not be a positive map.
    """

    def __init__(self, algebra: TypeIAlgebra, time: float = 1.0, n_jumps: int = 1, n_inter: int = 0, transitions=None,
                 hamiltonian: bool = True, init_scale: float = 0.1, method: str = 'exact', substeps: int = 20,
                 generator: Optional[torch.Generator] = None):
        super().__init__(algebra)
        if method not in ('exact', 'kraus'):
            raise ValueError("LindbladLayer: method must be 'exact' or 'kraus' to preserve complete positivity")
        self.algebra, self.time, self.method, self.substeps = algebra, float(time), method, substeps
        self._param = optimize.LindbladParam(algebra, 1, n_jumps=n_jumps, n_inter=n_inter, transitions=transitions,
                                             hamiltonian=hamiltonian, init_scale=init_scale, generator=generator)
        for name in ('X', 'L', 'K'):
            t = getattr(self._param, name)
            if t is not None:
                p = nn.Parameter(t.detach().clone())
                setattr(self._param, name, p)
                self.register_parameter(name, p)

    def _sync_parameterisation(self, device):
        self._param.alg = self.algebra
        self._param._mask = self._param._mask.to(device=device)
        if self._param.inter_mask is not None:
            self._param.inter_mask = self._param.inter_mask.to(device=device)
        for name in ('X', 'L', 'K'):
            if getattr(self._param, name) is not None:
                setattr(self._param, name, getattr(self, name))

    def generator(self):
        """(H, jumps) for dynamics.lindblad_evolve."""
        for name in ('X', 'L', 'K'):
            if getattr(self._param, name) is not None:
                setattr(self._param, name, getattr(self, name))
        return self._param.hamiltonian(), self._param.jumps()

    def forward(self, rho: DensityMatrix) -> DensityMatrix:
        H, jumps = self.generator()
        return dynamics.lindblad_evolve(rho, H, jumps, [0.0, self.time], substeps=self.substeps,
                                        method=self.method)[-1]


class Residual(nn.Module):
    """
    rho -> (1 - p) rho + p Phi(rho) for a layer Phi: a convex combination of channels, hence a channel (the plain
    rho + Phi(rho) of a ResNet is not). With learn_p the weight is p = sigmoid(logit), trained with the rest. A stack
    of n blocks with weight p = t/n and a fixed Phi tends to exp(t (Phi - id)), the semigroup of the GKSL generator
    Phi - id whose jump operators are the Kraus operators of Phi: residual channel networks with small blocks express
    Markovian (divisible) evolutions.
    """

    def __init__(self, layer: nn.Module, p: float = 0.1, learn_p: bool = False):
        super().__init__()
        if not 0.0 <= p <= 1.0:
            raise ValueError("Residual: p must lie in [0, 1]")
        self.layer = layer
        logit = math.log(p / (1 - p)) if 0.0 < p < 1.0 else (math.inf if p == 1.0 else -math.inf)
        t = torch.tensor(logit, dtype=torch.float64)
        if learn_p:
            self.logit = nn.Parameter(t)
        else:
            self.register_buffer('logit', t)

    @property
    def p(self) -> torch.Tensor:
        return torch.sigmoid(self.logit)

    def forward(self, rho: DensityMatrix) -> DensityMatrix:
        p = self.p.to(rho.matrix.real.dtype)
        out = self.layer(rho)
        return _as_state(rho.algebra, (1 - p) * rho.matrix + p * out.matrix.to(rho.matrix.dtype))


class PostSelect(_AlgebraLayer):
    """
    rho -> P rho P / Tr(P rho P) for a projection P of the algebra, given as an Operator or as the list of its weights
    on the sectors (a central projection, e.g. a fixed charge). The nonlinearity of the network comes from this
    renormalisation; its price is the acceptance Tr(P rho P), stored in self.acceptance (shape (batch,)) after every
    forward pass, e.g. for a penalty -log(acceptance).
    """

    def __init__(self, algebra: TypeIAlgebra, projection: Union[Operator, Sequence[float]], eps: float = 1e-12):
        super().__init__(algebra)
        self.algebra, self.eps = algebra, eps
        P = projection if isinstance(projection, Operator) else algebra.central(list(projection))
        self.register_buffer('P', P.matrix.detach().clone())
        self.acceptance = None

    def forward(self, rho: DensityMatrix) -> DensityMatrix:
        P = self.algebra.operator(self.P.to(rho.matrix.dtype))
        post, p = states.lueders_update(rho, P, self.eps)
        self.acceptance = p
        return post


class FeedForward(_AlgebraLayer):
    """
    Measure the instrument {Phi_j} (measurement.Instrument), keep the outcome j as classical information (in the
    centre) and apply the branch layer j: rho -> sum_j Phi'_j(Phi_j(rho)), a channel. Branches are modules (e.g.
    ChannelLayer, nn.Identity()). Quantum states cannot be copied forward (no-cloning, no-broadcasting); outcomes in
    the centre can, so this is the form a fan-out skip connection takes for channels.
    The instrument's Kraus tensor is a registered buffer; outcome labels are extra state, so state_dict
    restores the fixed measurement as well as branch weights. The constructor must have matching sector,
    outcome and Kraus dimensions. Legacy checkpoints without instrument data require reconstructing the
    original instrument explicitly and loading with strict=False.
    """

    def __init__(self, instrument, branches: Sequence[nn.Module]):
        super().__init__(instrument.algebra)
        if len(branches) != len(instrument):
            raise ValueError(f"FeedForward: {len(instrument)} outcomes but {len(branches)} branches")
        self.register_buffer('instrument_kraus', instrument.kraus.clone())
        self._instrument_outcomes = list(instrument.outcomes)
        self._instrument_trace_preserving = instrument._trace_preserving
        self.branches = nn.ModuleList(branches)

    @property
    def instrument(self):
        """Instrument backed by the current checkpointed Kraus tensor and algebra."""
        inst = measurement.Instrument(self.instrument_kraus, self.algebra, check=False,
                                      outcomes=self._instrument_outcomes)
        inst._trace_preserving = self._instrument_trace_preserving
        return inst

    def get_extra_state(self):
        return dict(version=1, outcomes=list(self._instrument_outcomes),
                    trace_preserving=self._instrument_trace_preserving)

    def set_extra_state(self, state):
        if state.get('version') != 1 or len(state['outcomes']) != len(self.branches):
            raise ValueError("FeedForward: incompatible instrument metadata in checkpoint")
        self._instrument_outcomes = list(state['outcomes'])
        self._instrument_trace_preserving = state['trace_preserving']

    def forward(self, rho: DensityMatrix) -> DensityMatrix:
        X = self.instrument._branches(rho)                                  # (B, J, C, k, k): Phi_j(rho), unnormalised
        out = 0
        for j, branch in enumerate(self.branches):
            out = out + branch(_as_state(rho.algebra, X[:, j])).matrix
        return _as_state(rho.algebra, out)


# ----------------------------------------------------------------------
# losses and bookkeeping
# ----------------------------------------------------------------------
def infidelity(rho: Operator, target: Operator) -> torch.Tensor:
    """1 - F(rho, target), F the Uhlmann fidelity (squared), per sample."""
    return 1 - states.fidelity(rho, target)


def trace_distance(rho: Operator, sigma: Operator) -> torch.Tensor:
    """(1/2) || rho - sigma ||_1 per sample."""
    return states.trace_distance(rho, sigma)


def relative_entropy(rho: Operator, sigma: Operator, eps: float = 1e-12) -> torch.Tensor:
    """D(rho || sigma) per sample."""
    return states.relative_entropy(rho, sigma, eps)


def negative_log_likelihood(povm, rho: Operator, counts, eps: float = 1e-12) -> torch.Tensor:
    """-sum_j n_j log Tr(rho E_j) per sample: the loss of maximum-likelihood tomography from measurement counts."""
    return -povm.log_likelihood(rho, counts, eps)


def count_parameters(module: nn.Module, effective: bool = True) -> int:
    """
    Number of real parameters of a network (complex entries count twice). With effective=True only entries that can
    influence the map are counted: the active blocks, the allowed transitions of channel and generator layers, and a
    Hamiltonian as a self-adjoint element (sum_c k_c^2 real numbers). The gauge freedom of Kraus operators is not
    subtracted.
    """
    total = 0
    for mod in module.modules():
        for name, p in mod.named_parameters(recurse=False):
            if not effective:
                total += 2 * p.numel() if p.is_complex() else p.numel()
                continue
            if isinstance(mod, ChannelLayer):
                total += 2 * int(mod._param._mask.sum())
            elif isinstance(mod, (UnitaryLayer, LindbladLayer)) and name in ('X', 'L', 'K'):
                block = int(sum(k * k for k in mod.algebra.k_factors))
                if name == 'X':
                    total += block
                elif name == 'L':
                    total += 2 * p.shape[1] * block
                else:
                    total += 2 * p.shape[1] * int(mod._param.inter_mask.sum())
            else:
                total += 2 * p.numel() if p.is_complex() else p.numel()
    return total
