"""
Time evolution on M = (+)_c M_{k_c} for Hamiltonians and Lindblad generators in M.

Because H and the jump operators are block diagonal, every sector evolves independently
(superselection) and all routines work block by block and batched over samples.

Exact solvers (time-independent generators):
    propagator, schrodinger, von_neumann        spectral decomposition of H
    lindblad_channel                            exp of the Lindblad superoperator (k_c^2 x k_c^2)
Fixed-step RK4 solvers (any generator, also time dependent):
    rk4, schrodinger_rk4, lindblad_evolve, solve_operator_ode
"""
from typing import Callable, List, Optional, Sequence, Union

import torch

from .algebra import TypeIAlgebra
from .channels import Channel, InterSectorChannel
from . import cost
from .states import DensityMatrix

Operator = TypeIAlgebra.Operator
OperatorLike = Union[Operator, torch.Tensor]


def _mat(x: OperatorLike) -> torch.Tensor:
    return x.matrix if isinstance(x, Operator) else x


def _require_complex(alg: TypeIAlgebra):
    if not alg.hilbert.complex_valued:
        raise ValueError("unitary dynamics needs a complex-valued algebra (complex_valued=True)")


def _times(times, device) -> torch.Tensor:
    return torch.as_tensor(times, dtype=torch.float64, device=device).reshape(-1)


# ----------------------------------------------------------------------
# Closed systems, exact
# ----------------------------------------------------------------------
def propagator(H: Operator, t: float) -> Operator:
    """U(t) = exp(-i t H) for self-adjoint H."""
    _require_complex(H.algebra)
    U = H.apply_function(lambda w: torch.exp(-1j * t * w))
    U._is_normal = True
    U._is_invertible = True
    return U


def schrodinger(psi0: torch.Tensor, H: Operator, times) -> torch.Tensor:
    """
    psi(t) = exp(-i t H) psi0 for state vectors psi0 of shape (batch, C, k_max) (block layout,
    global norm 1). Returns a tensor of shape (len(times), batch, C, k_max).
    """
    _require_complex(H.algebra)
    if psi0.dim() == 4:
        psi0 = psi0.squeeze(-1)
    ts = _times(times, psi0.device)
    out = torch.zeros(len(ts), *psi0.shape, dtype=H.matrix.dtype, device=psi0.device)
    for c, (w, V) in enumerate(H.eigh()):
        k_c = H.algebra.k_factors[c]
        if k_c == 0:
            continue
        coeff = (V.conj().transpose(-2, -1) @ psi0[:, c, :k_c, None].to(V.dtype)).squeeze(-1)  # (B, k)
        phase = torch.exp(-1j * ts[:, None, None] * w[None].double()).to(V.dtype)          # (T, B, k)
        out[:, :, c, :k_c] = (V[None] @ (phase * coeff[None]).unsqueeze(-1)).squeeze(-1)
    return out


def von_neumann(rho0: Operator, H: Operator, times) -> List[Operator]:
    """rho(t) = U(t) rho0 U(t)^* for every t in `times` (exact)."""
    _require_complex(H.algebra)
    alg = H.algebra
    blocks = H.eigh()
    rho = rho0.matrix.to(H.matrix.dtype)
    out = []
    for t in _times(times, rho.device).tolist():
        mat = torch.zeros_like(rho)
        for c, (w, V) in enumerate(blocks):
            k_c = alg.k_factors[c]
            if k_c == 0:
                continue
            U = (V * torch.exp(-1j * t * w.double()).to(V.dtype).unsqueeze(-2)) @ V.conj().transpose(-2, -1)
            mat[:, c, :k_c, :k_c] = U @ rho[:, c, :k_c, :k_c] @ U.conj().transpose(-2, -1)
        if isinstance(rho0, DensityMatrix):
            out.append(DensityMatrix(alg, matrix=mat, validate=False))
        else:
            out.append(alg.operator(mat, is_self_adjoint=rho0._is_self_adjoint, is_positive=rho0._is_positive))
    return out


# ----------------------------------------------------------------------
# Generic RK4
# ----------------------------------------------------------------------
def rk4(f: Callable[[float, torch.Tensor], torch.Tensor], y0: torch.Tensor, times,
        substeps: int = 10, progress: bool = False) -> torch.Tensor:
    """
    Classical 4th-order Runge-Kutta for dy/dt = f(t, y) on tensors. `times` is an increasing grid;
    each interval is split into `substeps` steps. Returns y at every grid point, shape (T, *y0.shape).
    After the first step a CostWarning reports the expected run time if it is long;
    progress=True shows a progress bar (tqdm if installed).
    """
    ts = _times(times, y0.device).tolist()
    y = y0
    out = [y0]
    cost.check_memory(len(ts) * y0.numel() * y0.element_size(), y0.device, f"rk4 trajectory ({len(ts)} time points)")
    timer = cost.StepTimer(max(len(ts) - 1, 0) * substeps, "rk4", progress)
    for t0, t1 in zip(ts[:-1], ts[1:]):
        h = (t1 - t0) / substeps
        t = t0
        for _ in range(substeps):
            k1 = f(t, y)
            k2 = f(t + h / 2, y + h / 2 * k1)
            k3 = f(t + h / 2, y + h / 2 * k2)
            k4 = f(t + h, y + h * k3)
            y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            t += h
            timer.step(y.device)
        out.append(y)
    timer.close()
    return torch.stack(out)


def solve_operator_ode(f: Callable[[float, torch.Tensor], torch.Tensor], X0: Operator, times,
                       substeps: int = 10, progress: bool = False) -> List[Operator]:
    """
    dX/dt = f(t, X) for X in the algebra, f acting on matrices of shape (batch, C, k_max, k_max)
    and returning block-diagonal matrices of the same shape (e.g. Heisenberg, Lyapunov or
    Riccati-type right-hand sides). Returns X at every time of the grid.
    """
    ys = rk4(f, X0.matrix, times, substeps, progress)
    return [X0.algebra.operator(y) for y in ys]


def schrodinger_rk4(psi0: torch.Tensor, H_of_t: Callable[[float], OperatorLike], times,
                    substeps: int = 10, progress: bool = False) -> torch.Tensor:
    """i dpsi/dt = H(t) psi for a (possibly time-dependent) Hamiltonian; shape (T, batch, C, k_max)."""
    if psi0.dim() == 4:
        psi0 = psi0.squeeze(-1)

    def f(t, psi):
        H = _mat(H_of_t(t)).to(psi.dtype)
        return -1j * (H @ psi.unsqueeze(-1)).squeeze(-1)
    return rk4(f, psi0.to(torch.promote_types(psi0.dtype, torch.complex64)), times, substeps, progress)


# ----------------------------------------------------------------------
# Open systems: Lindblad master equation
# ----------------------------------------------------------------------
def _rates(jumps, rates):
    if rates is None:
        return [1.0] * len(jumps)
    assert len(rates) == len(jumps)
    return list(rates)


def lindblad_rhs(H: Optional[OperatorLike], jumps: Sequence[OperatorLike] = (),
                 rates: Optional[Sequence[float]] = None) -> Callable[[float, torch.Tensor], torch.Tensor]:
    """
    Right-hand side of the GKSL equation
        drho/dt = -i[H, rho] + sum_j g_j (L_j rho L_j^* - 1/2 {L_j^* L_j, rho})
    as a function (t, rho_matrix) -> drho/dt. H and L_j are constant, or callables of t.

    A jump that moves weight between sectors (particle loss N -> N-1, ...) is passed as an
    InterSectorChannel of the algebra into itself: each of its Kraus operators K_i is one jump operator
    (blocks K_i^{dc}), contributing sum_i (K_i rho K_i^* - 1/2 {K_i^* K_i, rho}); K_i^* K_i is block
    diagonal, so the generator still maps M into M and conserves the total trace.
    """
    rates = _rates(jumps, rates)
    sector_jumps = {}
    for idx, L in enumerate(jumps):
        if isinstance(L, InterSectorChannel):
            if L.algebra_in.k_factors != L.algebra_out.k_factors:
                raise ValueError("a jump between sectors must map the algebra into itself")
            K = L.kraus                                                    # (B, r, D, C, m, k)
            sector_jumps[idx] = (K, (K.conj().transpose(-2, -1) @ K).sum(dim=(1, 2)))   # (B, C, k, k)

    def get(x, t):
        return _mat(x(t) if callable(x) else x)

    def f(t, rho):
        out = torch.zeros_like(rho)
        if H is not None:
            Hm = get(H, t).to(rho.dtype)
            out = out - 1j * (Hm @ rho - rho @ Hm)
        for idx, (g, L) in enumerate(zip(rates, jumps)):
            if idx in sector_jumps:
                K, KdK = sector_jumps[idx]
                Kd = K.to(rho.dtype)
                gain = (Kd @ rho[:, None, None] @ Kd.conj().transpose(-2, -1)).sum(dim=(1, 3))
                KdK = KdK.to(rho.dtype)
                out = out + g * (gain - 0.5 * (KdK @ rho + rho @ KdK))
                continue
            Lm = get(L, t).to(rho.dtype)
            Ld = Lm.conj().transpose(-2, -1)
            LdL = Ld @ Lm
            out = out + g * (Lm @ rho @ Ld - 0.5 * (LdL @ rho + rho @ LdL))
        return out
    return f


def lindblad_evolve(rho0: Operator, H: Optional[OperatorLike], jumps: Sequence[OperatorLike] = (),
                    times=(0.0, 1.0), rates: Optional[Sequence[float]] = None,
                    substeps: int = 20, progress: bool = False) -> List[Operator]:
    """Integrate the Lindblad equation with RK4; returns rho(t) for every t in `times`."""
    rdt = rho0.matrix.dtype
    dtype = torch.promote_types(rdt, torch.complex64) if rho0.algebra.hilbert.complex_valued or H is not None else rdt
    if H is not None:
        _require_complex(rho0.algebra)
    ys = rk4(lindblad_rhs(H, jumps, rates), rho0.matrix.to(dtype), times, substeps, progress)
    alg = rho0.algebra
    if isinstance(rho0, DensityMatrix):
        return [DensityMatrix(alg, matrix=0.5 * (y + y.conj().transpose(-2, -1)), validate=False) for y in ys]
    return [alg.operator(0.5 * (y + y.conj().transpose(-2, -1)), is_self_adjoint=True) for y in ys]


def _kron(A: torch.Tensor, B: torch.Tensor) -> torch.Tensor:
    k = A.shape[-1]
    batch = torch.broadcast_shapes(A.shape[:-2], B.shape[:-2])
    return (A[..., :, None, :, None] * B[..., None, :, None, :]).reshape(*batch, k * k, k * k)


def lindblad_superoperator(alg: TypeIAlgebra, H: Optional[OperatorLike], jumps: Sequence[OperatorLike] = (),
                           rates: Optional[Sequence[float]] = None) -> List[torch.Tensor]:
    """Per-channel generator matrices acting on row-major vec(rho), shape (batch, k_c^2, k_c^2)."""
    rates = _rates(jumps, rates)
    mats = ([_mat(H)] if H is not None else []) + [_mat(L) for L in jumps]
    if not mats:
        raise ValueError("need a Hamiltonian or at least one jump operator")
    dtype = torch.complex128 if H is not None else torch.promote_types(mats[0].dtype, torch.float64)
    batch = max(m.shape[0] for m in mats)
    # peak: the generator, one Kronecker product per term and the exponential / Choi / eigenvector
    # workspaces of lindblad_channel, all of size k_c^4 per sample
    cost.check_memory((len(jumps) + 10) * max(cost.tensor_bytes((batch, k * k, k * k), dtype) for k in alg.k_factors),
                      mats[0].device, f"Lindblad superoperator / exp(tL) (batch {batch}, k_c^2 x k_c^2 per sector)")
    out = []
    for c, k_c in enumerate(alg.k_factors):
        eye = torch.eye(k_c, dtype=dtype, device=mats[0].device)
        G = torch.zeros(batch, k_c * k_c, k_c * k_c, dtype=dtype, device=mats[0].device)
        if H is not None:
            Hc = _mat(H)[:, c, :k_c, :k_c].to(dtype)
            G = G - 1j * (_kron(Hc, eye) - _kron(eye, Hc.transpose(-2, -1)))
        for g, L in zip(rates, jumps):
            Lc = _mat(L)[:, c, :k_c, :k_c].to(dtype)
            LdL = Lc.conj().transpose(-2, -1) @ Lc
            G = G + g * (_kron(Lc, Lc.conj()) - 0.5 * _kron(LdL, eye) - 0.5 * _kron(eye, LdL.transpose(-2, -1)))
        out.append(G)
    return out


def lindblad_channel(alg: TypeIAlgebra, H: Optional[OperatorLike], jumps: Sequence[OperatorLike] = (),
                     t: float = 1.0, rates: Optional[Sequence[float]] = None) -> Channel:
    """
    The channel exp(t L) of a time-independent Lindblad generator, computed exactly in double
    precision from the k_c^2 x k_c^2 superoperators (practical for k_c up to a few tens).
    """
    if H is not None:
        _require_complex(alg)
    S = [cost.batched_call(torch.linalg.matrix_exp, t * G, f"exp(tL) (sector {c}, size {G.shape[-1]})", kind='expm')
         for c, G in enumerate(lindblad_superoperator(alg, H, jumps, rates))]
    return Channel.from_superoperator(alg, S)
