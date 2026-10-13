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

import warnings

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
    psi0 = psi0.to(H.matrix.device)                       # a state built on the CPU may meet H on a GPU
    B = max(psi0.shape[0], H.matrix.shape[0])             # one initial state for a batch of Hamiltonians, or vice versa
    psi0 = psi0.expand(B, *psi0.shape[1:])
    ts = _times(times, psi0.device)
    if psi0.dtype in (torch.float64, torch.complex128) and H.matrix.dtype in (torch.float32, torch.complex64):
        warnings.warn("schrodinger: psi0 is in double but H in single precision; the evolution is computed in single "
                      "precision (create H in a double-precision algebra)", UserWarning, stacklevel=2)
    alg = H.algebra
    Bh = H.matrix.shape[0]
    out = torch.zeros(len(ts), *psi0.shape, dtype=H.matrix.dtype, device=psi0.device)
    for k, idx, _, (w, V) in H._decompositions('hermitian'):      # one decomposition per block size
        n = len(idx)
        sel = alg._sel(idx)
        V = V.reshape(Bh, n, k, k)
        coeff = (V.conj().transpose(-2, -1) @ psi0[:, sel, :k, None].to(V.dtype)).squeeze(-1)   # (B, n, k)
        phase = torch.exp(-1j * ts[:, None, None, None] * w.reshape(Bh, n, k)[None].double()).to(V.dtype)
        out[:, :, sel, :k] = (V[None] @ (phase * coeff[None]).unsqueeze(-1)).squeeze(-1)     # (T, B, n, k)
    return out


def von_neumann(rho0: Operator, H: Operator, times) -> List[Operator]:
    """rho(t) = U(t) rho0 U(t)^* for every t in `times` (exact)."""
    _require_complex(H.algebra)
    alg = H.algebra
    rho = rho0.matrix.to(H.matrix.dtype)
    Bh = H.matrix.shape[0]
    B = max(rho.shape[0], Bh)                             # one initial state for a batch of Hamiltonians, or vice versa
    rho = rho.expand(B, *rho.shape[1:])
    ts = _times(times, rho.device)
    T = len(ts)
    cost.check_memory(2 * T * rho.numel() * rho.element_size(), rho.device, f"von_neumann ({T} time points)")
    traj = torch.zeros(T, *rho.shape, dtype=rho.dtype, device=rho.device)
    for k, idx, _, (w, V) in H._decompositions('hermitian'):      # one decomposition per block size
        n = len(idx)
        sel = alg._sel(idx)
        V = V.reshape(Bh, n, k, k)
        Vh = V.conj().transpose(-2, -1)
        R = Vh @ rho[:, sel, :k, :k] @ V                                  # rho in the eigenbasis, (B, n, k, k)
        wd = w.reshape(Bh, n, k).double()
        # all times at once: rho(t)_ij = exp(-i t (w_i - w_j)) rho_ij in the eigenbasis (phases in double)
        phase = torch.exp(-1j * ts[:, None, None, None, None] * (wd[..., :, None] - wd[..., None, :])[None])
        traj[:, :, sel, :k, :k] = V[None] @ (phase.to(R.dtype) * R[None]) @ Vh[None]
    out = []
    for mat in traj:
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
    ts = _times(times, 'cpu').tolist()                     # the grid stays on the host (no device round trip)
    y = y0
    cost.check_memory(len(ts) * y0.numel() * y0.element_size(), y0.device, f"rk4 trajectory ({len(ts)} time points)")
    out = y0.new_empty((len(ts), *y0.shape))
    out[0] = y0
    timer = cost.StepTimer(max(len(ts) - 1, 0) * substeps, "rk4", progress)
    for i, (t0, t1) in enumerate(zip(ts[:-1], ts[1:])):
        h = (t1 - t0) / substeps
        for s_ in range(substeps):
            t = t0 + s_ * h                                 # no accumulated rounding of t
            k1 = f(t, y)
            k2 = f(t + h / 2, y + h / 2 * k1)
            k3 = f(t + h / 2, y + h / 2 * k2)
            k4 = f(t + h, y + h * k3)
            y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            timer.step(y.device)
        out[i + 1] = y
    timer.close()
    return out


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
    const = None if callable(H_of_t) and not isinstance(H_of_t, Operator) else _mat(H_of_t)
    if const is not None:
        psi0 = psi0.to(const.device)
    y0 = psi0.to(torch.promote_types(psi0.dtype, torch.complex64))
    Hc = None if const is None else (-1j * const.to(y0.dtype))      # hoisted: constant H is cast once

    def f(t, psi):
        A = Hc if Hc is not None else -1j * _mat(H_of_t(t)).to(psi.dtype)
        return (A @ psi.unsqueeze(-1)).squeeze(-1)
    return rk4(f, y0, times, substeps, progress)


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
    InterSectorChannel Phi of the algebra into itself; it contributes Phi(rho) - 1/2 {Phi^*(1), rho}, where
    Phi^*(1)_c = sum_{d,i} K_i^{dc*} K_i^{dc} is block diagonal. This is the GKSL generator of the jumps
    K_i^{dc} compressed to M (the conditional expectation onto M applied to the generator on B(H)); it maps M
    into M and conserves the total trace. Constant H and jumps are combined once (H_eff = H - i/2 sum g L^* L);
    callables of t are evaluated at every step.
    """
    rates = _rates(jumps, rates)
    for L in jumps:
        if isinstance(L, Channel):
            raise TypeError("a sector-preserving Channel is not a jump operator; pass its Kraus operators, or "
                            "channel.to_inter_sector() for a jump between sectors")
    inter = [(g, L) for g, L in zip(rates, jumps) if isinstance(L, InterSectorChannel)]
    for _, L in inter:
        if L.algebra_in.k_factors != L.algebra_out.k_factors:
            raise ValueError("a jump between sectors must map the algebra into itself")
    plain = [(g, L) for g, L in zip(rates, jumps) if not isinstance(L, InterSectorChannel)]
    timedep = (callable(H) and not isinstance(H, Operator)) or any(callable(L) and not isinstance(L, Operator)
                                                                     for _, L in plain)

    def get(x, t):
        return _mat(x(t) if callable(x) and not isinstance(x, Operator) else x)

    def parts(t, dtype, device):
        """(H_eff, stacked jumps (B, J, C, k, k) with their rates, inter-sector Kraus tensors): with
        A = sum_j g_j L_j^* L_j, drho/dt = -i (H_eff rho - rho H_eff^*) + sum_j g_j L_j rho L_j^* + (between
        sectors) with H_eff = H - (i/2) A."""
        Heff = None
        if H is not None:
            Heff = get(H, t).to(device=device, dtype=dtype)
        A = None
        Ls, gs = None, None
        if plain:
            Ls = torch.stack(torch.broadcast_tensors(*[get(L, t).to(device=device, dtype=dtype) for _, L in plain]),
                             dim=1)                                           # (B, J, C, k, k)
            gs = torch.tensor([g for g, _ in plain], dtype=Ls.real.dtype, device=device)
            A = torch.einsum('j,bjcai,bjcak->bcik', gs.to(dtype), Ls.conj(), Ls)
        Ks = []
        for g, L in inter:
            K = L.kraus.to(device=device, dtype=dtype)                        # (B, r, D, C, m, k)
            KdK = (K.conj().transpose(-2, -1) @ K).sum(dim=(1, 2))            # (B, C, k, k), block diagonal
            A = g * KdK if A is None else A + g * KdK
            Ks.append((g, K))
        if A is not None:
            Heff = -0.5j * A if Heff is None else Heff - 0.5j * A
        return Heff, Ls, gs, Ks

    cache = {}

    def f(t, rho):
        key = (rho.dtype, rho.device)
        if timedep or key not in cache:
            cache[key] = parts(t, rho.dtype, rho.device)
        Heff, Ls, gs, Ks = cache[key]
        out = torch.zeros_like(rho) if Heff is None else -1j * (Heff @ rho - rho @ Heff.conj().transpose(-2, -1))
        if Ls is not None:                     # sum_j g_j L_j rho L_j^*: one product over the merged (j, column) index
            X = (Ls @ rho.unsqueeze(1)) * gs[:, None, None, None].to(rho.dtype)          # (B, J, C, k, k)
            Bx, J, C, k = X.shape[0], X.shape[1], X.shape[2], X.shape[-1]
            Xc = X.permute(0, 2, 3, 1, 4).reshape(Bx, C, k, J * k)
            Lc = Ls.permute(0, 2, 3, 1, 4).reshape(Ls.shape[0], C, k, J * k)
            out = out + Xc @ Lc.conj().transpose(-2, -1)
        for g, K in Ks:                        # gain between sectors: sum_{i, c} K^{dc} rho_c K^{dc*}
            B_, r, D, C, m, k = K.shape
            X = K @ rho[:, None, None]                                          # (B, r, D, C, m, k)
            Bx = X.shape[0]
            Xc = X.permute(0, 2, 4, 1, 3, 5).reshape(Bx, D, m, r * C * k)
            Kc = K.permute(0, 2, 4, 1, 3, 5).reshape(B_, D, m, r * C * k)
            out = out + g * (Xc @ Kc.conj().transpose(-2, -1))
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
    # one initial state for a batch of Hamiltonians / jump operators (disorder averages), or vice versa
    t0 = float(_times(times, 'cpu')[0])

    def _batch(X):                                     # batch size of a (possibly time-dependent) H or jump
        if isinstance(X, InterSectorChannel):
            return X.kraus.shape[0]
        X = _mat(X(t0)) if callable(X) and not isinstance(X, Operator) else _mat(X)
        return X.shape[0] if isinstance(X, torch.Tensor) else 1
    batches = [rho0.matrix.shape[0]] + [_batch(X) for X in ([H] if H is not None else []) + list(jumps)]
    y0 = rho0.matrix.to(dtype)
    y0 = y0.expand(max(batches), *y0.shape[1:]).clone()
    ys = rk4(lindblad_rhs(H, jumps, rates), y0, times, substeps, progress)
    alg = rho0.algebra
    ys = 0.5 * (ys + ys.conj().transpose(-2, -1))          # once for the whole trajectory
    if isinstance(rho0, DensityMatrix):
        return [DensityMatrix(alg, matrix=y, validate=False) for y in ys]
    return [alg.operator(y, is_self_adjoint=True) for y in ys]


def _kron(A: torch.Tensor, B: torch.Tensor) -> torch.Tensor:
    k = A.shape[-1]
    batch = torch.broadcast_shapes(A.shape[:-2], B.shape[:-2])
    return (A[..., :, None, :, None] * B[..., None, :, None, :]).reshape(*batch, k * k, k * k)


def _lindblad_generators(alg: TypeIAlgebra, H: Optional[OperatorLike], jumps: Sequence[OperatorLike] = (),
                         rates: Optional[Sequence[float]] = None):
    """[(k, idx, G)] with G (batch, n, k^2, k^2) the generators of the n sectors of block size k (row-major vec)."""
    rates = _rates(jumps, rates)
    for L in jumps:
        if type(L).__name__ == 'InterSectorChannel':
            raise ValueError("lindblad_superoperator / lindblad_channel act sector by sector and do not accept jumps "
                             "between sectors (InterSectorChannel); use lindblad_evolve for those")
        if not isinstance(_mat(L), torch.Tensor):
            raise TypeError(f"jump operators must be Operators or tensors of shape (batch, C, k_max, k_max), "
                            f"got {type(L).__name__}")
    mats = ([_mat(H)] if H is not None else []) + [_mat(L) for L in jumps]
    if not mats:
        raise ValueError("need a Hamiltonian or at least one jump operator")
    dtype = torch.float64
    for m in mats:
        dtype = torch.promote_types(dtype, m.dtype)
    if H is not None:
        dtype = torch.complex128
    dev = mats[0].device
    batch = max(m.shape[0] for m in mats)
    # peak: the generator, a few Kronecker products and the exponential / Choi / eigenvector workspaces of
    # lindblad_channel, all of size k_c^4 per sector and sample
    cost.check_memory(12 * sum(cost.tensor_bytes((batch, k * k, k * k), dtype) for k in alg.k_factors),
                      dev, f"Lindblad superoperator / exp(tL) (batch {batch}, k_c^2 x k_c^2 per sector)")
    out = []
    for k, idx in alg._size_groups():
        sel = alg._sel(idx)
        n = len(idx)
        eye = torch.eye(k, dtype=dtype, device=dev)
        G = torch.zeros(batch, n, k * k, k * k, dtype=dtype, device=dev)
        A = None
        if jumps:
            Ls = torch.stack(torch.broadcast_tensors(*[_mat(L)[:, sel, :k, :k].to(device=dev, dtype=dtype)
                                                       for L in jumps]), dim=1)      # (B, J, n, k, k)
            g = torch.tensor(rates, dtype=dtype, device=dev)
            # sum_j g_j L_j (x) conj(L_j) in one contraction, and A = sum_j g_j L_j^* L_j
            G = G + torch.einsum('j,bjnai,bjnck->bnacik', g, Ls, Ls.conj()).reshape(Ls.shape[0], n, k * k, k * k)
            A = torch.einsum('j,bjnai,bjnak->bnik', g, Ls.conj(), Ls)
        if H is not None:
            Hc = _mat(H)[:, sel, :k, :k].to(device=dev, dtype=dtype)
            G = G - 1j * (_kron(Hc, eye) - _kron(eye, Hc.transpose(-2, -1)))
        if A is not None:
            G = G - 0.5 * (_kron(A, eye) + _kron(eye, A.transpose(-2, -1)))
        out.append((k, idx, G))
    return out


def lindblad_superoperator(alg: TypeIAlgebra, H: Optional[OperatorLike], jumps: Sequence[OperatorLike] = (),
                           rates: Optional[Sequence[float]] = None) -> List[torch.Tensor]:
    """Per-channel generator matrices acting on row-major vec(rho), shape (batch, k_c^2, k_c^2)."""
    groups = _lindblad_generators(alg, H, jumps, rates)
    dtype, dev = groups[0][2].dtype, groups[0][2].device
    out = [torch.zeros(groups[0][2].shape[0], 0, 0, dtype=dtype, device=dev)] * alg.C
    for k, idx, G in groups:
        for j, c in enumerate(idx):
            out[c] = G[:, j]
    return out


def lindblad_channel(alg: TypeIAlgebra, H: Optional[OperatorLike], jumps: Sequence[OperatorLike] = (),
                     t: float = 1.0, rates: Optional[Sequence[float]] = None) -> Channel:
    """
    The channel exp(t L) of a time-independent Lindblad generator, computed exactly in double
    precision from the k_c^2 x k_c^2 superoperators (practical for k_c up to a few tens); one matrix exponential
    and one Choi eigendecomposition per block size.
    """
    if H is not None:
        _require_complex(alg)
    groups = _lindblad_generators(alg, H, jumps, rates)
    S = [None] * alg.C
    for k, idx, G in groups:
        E = cost.batched_call(torch.linalg.matrix_exp, (t * G).reshape(-1, k * k, k * k),
                              f"exp(tL) (blocks of size {k * k})", kind='expm').reshape(G.shape)
        for j, c in enumerate(idx):
            S[c] = E[:, j]
    for c, k_c in enumerate(alg.k_factors):
        if S[c] is None:
            S[c] = torch.zeros(groups[0][2].shape[0], 0, 0, dtype=groups[0][2].dtype, device=groups[0][2].device)
    ch = Channel.from_superoperator(alg, S)
    ch._trace_preserving = True                   # exp(tL) of a GKSL generator is trace preserving
    return ch
