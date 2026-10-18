"""
Time evolution on M = (+)_c M_{k_c} for Hamiltonians and Lindblad generators in M.

Because H and the jump operators are block diagonal, every sector evolves independently
(superselection) and all routines work block by block and batched over samples.

Exact solvers (time-independent generators):
    propagator, schrodinger, von_neumann        spectral decomposition of H
    lindblad_channel                            exp of the Lindblad superoperator (k_c^2 x k_c^2)
Fixed-step RK4 solvers (any generator, also time dependent):
    rk4, schrodinger_rk4, lindblad_evolve, solve_operator_ode
lindblad_evolve(method='kraus' | 'exact'): time steps that are channels (CPTP) for every step size
"""
from typing import Callable, List, Optional, Sequence, Union

import warnings

import torch

from .algebra import TypeIAlgebra
from .channels import Channel, InterSectorChannel, SectorJump
from . import cost
from .states import DensityMatrix

Operator = TypeIAlgebra.Operator
OperatorLike = Union[Operator, torch.Tensor]


def _mat(x: OperatorLike, alg: Optional[TypeIAlgebra] = None) -> torch.Tensor:
    """Matrix of an operator; a raw tensor is masked to the blocks of `alg` (its padding may hold garbage)."""
    if isinstance(x, Operator):
        return x.matrix
    if alg is not None and min(alg.k_factors) != alg.k_max:
        return x * alg._const('block_mask').to(x.device)
    return x


def _algebra_of(*xs) -> Optional[TypeIAlgebra]:
    for x in xs:
        if isinstance(x, Operator):
            return x.algebra
    return None


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
    if torch.is_grad_enabled() and H.matrix.requires_grad:        # differentiable in H: propagators (Daleckii-Krein)
        return torch.stack([(propagator(H, float(t)).matrix @ psi0.to(H.matrix.dtype).unsqueeze(-1)).squeeze(-1)
                            for t in ts])
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
    differentiable = torch.is_grad_enabled() and H.matrix.requires_grad
    if differentiable:                                            # differentiable in H: propagators (Daleckii-Krein)
        Us = [propagator(H, float(t)).matrix for t in ts]
        traj = torch.stack([U @ rho @ U.conj().transpose(-2, -1) for U in Us])
    else:
        traj = torch.zeros(T, *rho.shape, dtype=rho.dtype, device=rho.device)
    for k, idx, _, (w, V) in ([] if differentiable else H._decompositions('hermitian')):
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
                    substeps: int = 10, progress: bool = False,
                    algebra: Optional[TypeIAlgebra] = None) -> torch.Tensor:
    """i dpsi/dt = H(t) psi for a (possibly time-dependent) Hamiltonian; shape (T, batch, C, k_max).
    Raw tensors (psi0, or H given as a tensor) are restricted to the blocks of `algebra` (taken from H when
    it is an Operator); without an algebra their padding must be zero."""
    if psi0.dim() == 4:
        psi0 = psi0.squeeze(-1)
    if algebra is None:
        algebra = _algebra_of(H_of_t)
    const = None if callable(H_of_t) and not isinstance(H_of_t, Operator) else _mat(H_of_t, algebra)
    if const is not None:
        psi0 = psi0.to(const.device)
    if algebra is not None and min(algebra.k_factors) != algebra.k_max:
        psi0 = psi0 * algebra._const('rows').to(psi0.device)        # entries in the padding are not part of psi
    y0 = psi0.to(torch.promote_types(psi0.dtype, torch.complex64))
    Hc = None if const is None else (-1j * const.to(y0.dtype))      # hoisted: constant H is cast once

    def f(t, psi):
        A = Hc if Hc is not None else -1j * _mat(H_of_t(t), algebra).to(psi.dtype)
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


def _lindblad_terms(H, jumps, rates, algebra):
    """Validate the data of a GKSL generator (lindblad_rhs) and split the jumps: returns
    (plain [(g, L)], inter [(g, InterSectorChannel)], timedep, algebra)."""
    rates = _rates(jumps, rates)
    for L in jumps:
        if isinstance(L, Channel):
            raise TypeError("a sector-preserving Channel is not a jump operator; pass its Kraus operators, or "
                            "channel.to_inter_sector() for a jump between sectors")
    inter = [(g, L) for g, L in zip(rates, jumps) if isinstance(L, (InterSectorChannel, SectorJump))]
    for _, L in inter:
        if isinstance(L, InterSectorChannel) and L.algebra_in.k_factors != L.algebra_out.k_factors:
            raise ValueError("a jump between sectors must map the algebra into itself")
    plain = [(g, L) for g, L in zip(rates, jumps) if not isinstance(L, (InterSectorChannel, SectorJump))]
    timedep = (callable(H) and not isinstance(H, Operator)) or any(callable(L) and not isinstance(L, Operator)
                                                                     for _, L in plain)
    if algebra is None:
        algebra = _algebra_of(H, *[L for _, L in plain])
    return plain, inter, timedep, algebra


def _lindblad_parts(H, plain, inter, algebra, t, dtype, device):
    """(H, A, stacked jumps (B, J, C, k, k) with their rates, inter-sector Kraus tensors) at time t: with
    A = sum_j g_j L_j^* L_j, drho/dt = -i (H rho - rho H^*) - 1/2 (A rho + rho A) + sum_j g_j L_j rho L_j^*
    + (between sectors). A is kept apart from H so that a real algebra with real jumps stays real."""
    def get(x):
        return _mat(x(t) if callable(x) and not isinstance(x, Operator) else x, algebra)

    Hm = None
    if H is not None:
        Hm = get(H).to(device=device, dtype=dtype)
    A = None
    Ls, gs = None, None
    if plain:
        Ls = torch.stack(torch.broadcast_tensors(*[get(L).to(device=device, dtype=dtype) for _, L in plain]),
                         dim=1)                                           # (B, J, C, k, k)
        gs = torch.tensor([g for g, _ in plain], dtype=Ls.real.dtype, device=device)
        A = torch.einsum('j,bjcai,bjcak->bcik', gs.to(dtype), Ls.conj(), Ls)
    Ks = []
    for g, L in inter:
        if isinstance(L, SectorJump):                                     # sparse: L^*L diagonal, never formed densely
            LdL = L.LdL().to(device=device, dtype=dtype)
            A = g * LdL if A is None else A + g * LdL
            Ks.append((g, L))
            continue
        K = L.kraus.to(device=device, dtype=dtype)                        # (B, r, D, C, m, k)
        KdK = (K.conj().transpose(-2, -1) @ K).sum(dim=(1, 2))            # (B, C, k, k), block diagonal
        A = g * KdK if A is None else A + g * KdK
        Ks.append((g, K))
    return Hm, A, Ls, gs, Ks


def _add_gain(out: torch.Tensor, rho: torch.Tensor, Ls, gs, Ks) -> torch.Tensor:
    """out + Phi(rho), Phi(rho) = sum_j g_j L_j rho L_j^* + sum_m g_m Phi_m(rho) (the jump part of the generator)."""
    if Ls is not None:                     # sum_j g_j L_j rho L_j^*: one product over the merged (j, column) index
        X = (Ls @ rho.unsqueeze(1)) * gs[:, None, None, None].to(rho.dtype)          # (B, J, C, k, k)
        Bx, J, C, k = X.shape[0], X.shape[1], X.shape[2], X.shape[-1]
        Xc = X.permute(0, 2, 3, 1, 4).reshape(Bx, C, k, J * k)
        Lc = Ls.permute(0, 2, 3, 1, 4).reshape(Ls.shape[0], C, k, J * k)
        out = out + Xc @ Lc.conj().transpose(-2, -1)
    for g, K in Ks:                        # gain between sectors: sum_{i, c} K^{dc} rho_c K^{dc*}
        if isinstance(K, SectorJump):      # sparse jump: a weighted permutation of the entries, O(sum_c k_c^2)
            out = out + g * K.sandwich(rho)
            continue
        B_, r, D, C, m, k = K.shape
        X = K @ rho[:, None, None]                                          # (B, r, D, C, m, k)
        Bx = X.shape[0]
        Xc = X.permute(0, 2, 4, 1, 3, 5).reshape(Bx, D, m, r * C * k)
        Kc = K.permute(0, 2, 4, 1, 3, 5).reshape(B_, D, m, r * C * k)
        out = out + g * (Xc @ Kc.conj().transpose(-2, -1))
    return out


def _gain_dual(Z: torch.Tensor, Ls, gs, Ks) -> torch.Tensor:
    """Phi^*(Z) = sum_j g_j L_j^* Z L_j + sum_m g_m Phi_m^*(Z), the dual of the jump part for the blunt trace
    (Phi_m^*(Z)_c = sum_{d,i} K_i^{dc*} Z_d K_i^{dc})."""
    out = torch.zeros_like(Z)
    if Ls is not None:
        X = Ls.conj().transpose(-2, -1) @ Z.unsqueeze(1) @ Ls                    # (B, J, C, k, k)
        out = out + (X * gs[:, None, None, None].to(Z.dtype)).sum(dim=1)
    for g, K in Ks:
        if isinstance(K, SectorJump):
            out = out + g * K.dual_sandwich(Z)
            continue
        X = K.conj().transpose(-2, -1) @ Z[:, None, :, None] @ K            # (B, r, D, C, k, k)
        out = out + g * X.sum(dim=(1, 2))
    return out


def lindblad_rhs(H: Optional[OperatorLike], jumps: Sequence[OperatorLike] = (),
                 rates: Optional[Sequence[float]] = None,
                 algebra: Optional[TypeIAlgebra] = None) -> Callable[[float, torch.Tensor], torch.Tensor]:
    """
    Right-hand side of the GKSL equation
        drho/dt = -i[H, rho] + sum_j g_j (L_j rho L_j^* - 1/2 {L_j^* L_j, rho})
    as a function (t, rho_matrix) -> drho/dt. H and L_j are constant, or callables of t.

    A jump that moves weight between sectors (particle loss N -> N-1, ...) is passed as an
    InterSectorChannel Phi of the algebra into itself (or, for spin flips and other jumps that map basis vectors to
    multiples of basis vectors, as a sparse SectorJump); it contributes Phi(rho) - 1/2 {Phi^*(1), rho}, where
    Phi^*(1)_c = sum_{d,i} K_i^{dc*} K_i^{dc} is block diagonal. This is the GKSL generator of the jumps
    K_i^{dc} compressed to M (the conditional expectation onto M applied to the generator on B(H)); it maps M
    into M and conserves the total trace. For constant H and jumps, A = sum g L^* L and the stacked jumps are
    prepared once; callables of t are evaluated at every step. Raw tensors are restricted to the blocks of
    `algebra` (taken from H or the jumps when they are Operators); without an algebra their padding must be zero.
    """
    plain, inter, timedep, algebra = _lindblad_terms(H, jumps, rates, algebra)
    cache = {}

    def f(t, rho):
        key = (rho.dtype, rho.device)
        if timedep or key not in cache:
            cache[key] = _lindblad_parts(H, plain, inter, algebra, t, rho.dtype, rho.device)
        Hm, A, Ls, gs, Ks = cache[key]
        out = torch.zeros_like(rho) if Hm is None else -1j * (Hm @ rho - rho @ Hm.conj().transpose(-2, -1))
        if A is not None:
            out = out - 0.5 * (A @ rho + rho @ A)
        return _add_gain(out, rho, Ls, gs, Ks)
    return f


def _kraus_step(parts, h: float, alg: TypeIAlgebra, dtype, device) -> Callable[[torch.Tensor], torch.Tensor]:
    """
    One step of the completely positive, trace-preserving splitting of Theorem 8.16 for the frozen generator
    `parts`: with E = exp(-i h/2 H_eff), H_eff = H - i/2 A, S(X) = E X E^*, and Phi the jump part,
        Psi_h = S o (id + h Phi + h^2/2 Phi^2) o S,    M = Psi_h^*(1),    rho -> Psi_h(M^{-1/2} rho M^{-1/2}).
    Every factor is CP and the last one makes the step exactly trace preserving. In floating point the trace is
    preserved up to about eps * cond(M); M is ill-conditioned when h times the spread of the decay rates (the
    eigenvalues of A) is large, and then the step is taken as two steps of size h/2 (recursively), each a channel.
    """
    Hm, A, Ls, gs, Ks = parts
    mask = alg._const('block_mask', torch.bool).to(device)
    eye = alg._const('eye', dtype).to(device)
    G = None                                                          # -i h/2 H_eff = -i h/2 H - h/4 A
    if Hm is not None:
        G = (-0.5j * h) * Hm
    if A is not None:
        G = -0.25 * h * A if G is None else G - 0.25 * h * A
    if G is None:
        E = eye[None]
    else:
        k = G.shape[-1]
        E = cost.batched_call(torch.linalg.matrix_exp, G.reshape(-1, k, k), "exp(-i h/2 H_eff)",
                              kind='expm').reshape(G.shape) * mask      # the padding of exp(G) is the identity
    jumps = Ls is not None or len(Ks) > 0

    def phi2(Y):                                                      # Y + h Phi(Y) + h^2/2 Phi(Phi(Y))
        if not jumps:
            return Y
        Y1 = _add_gain(torch.zeros_like(Y), Y, Ls, gs, Ks)
        Y2 = _add_gain(torch.zeros_like(Y1), Y1, Ls, gs, Ks)
        return Y + h * Y1 + (0.5 * h * h) * Y2

    def phi2_dual(W):
        if not jumps:
            return W
        W1 = _gain_dual(W, Ls, gs, Ks)
        W2 = _gain_dual(W1, Ls, gs, Ks)
        return W + h * W1 + (0.5 * h * h) * W2

    Eh = E.conj().transpose(-2, -1)
    M = Eh @ phi2_dual(Eh @ E) @ E                                    # Psi_h^*(1) >= (E^2)^* E^2 > 0
    M = 0.5 * (M + M.conj().transpose(-2, -1))
    pad = torch.eye(M.shape[-1], dtype=dtype, device=device) - eye   # identity on the padding, for eigh
    w, V = torch.linalg.eigh(M + pad)
    if h > 0 and (w[..., -1] / w[..., 0].clamp(min=torch.finfo(w.dtype).tiny)).max().item() \
            > torch.finfo(w.dtype).eps ** -0.25:
        half = _kraus_step(parts, 0.5 * h, alg, dtype, device)
        return lambda rho: half(half(rho))
    N = (V * w.clamp(min=torch.finfo(w.dtype).tiny).rsqrt().to(dtype).unsqueeze(-2)) @ V.conj().transpose(-2, -1)
    N = N * mask                                                      # M^{-1/2}, block diagonal

    def step(rho):
        Y = E @ (N @ rho @ N) @ Eh
        return E @ phi2(Y) @ Eh
    return step


def _generator_matrix(parts, alg: TypeIAlgebra, dtype, device) -> torch.Tensor:
    """Matrix of the generator on (+)_c M_{k_max} (zero padded blocks, row-major vec of every block, the sectors one
    after the other), shape (batch, C k_max^2, C k_max^2); jumps between sectors give off-diagonal sector blocks."""
    Hm, A, Ls, gs, Ks = parts
    C, k = alg.C, alg.k_max
    eye = torch.eye(k, dtype=dtype, device=device)
    blocks = torch.zeros(1, C, k * k, k * k, dtype=dtype, device=device)
    if Hm is not None:
        blocks = blocks - 1j * (_kron(Hm, eye) - _kron(eye, Hm.transpose(-2, -1)))
    if A is not None:
        blocks = blocks - 0.5 * (_kron(A, eye) + _kron(eye, A.transpose(-2, -1)))
    if Ls is not None:                                                # sum_j g_j L_j (x) conj(L_j)
        blocks = blocks + torch.einsum('j,bjcai,bjcel->bcaeil', gs.to(dtype), Ls, Ls.conj()).reshape(
            Ls.shape[0], C, k * k, k * k)
    G = torch.einsum('bcij,cd->bcidj', blocks, torch.eye(C, dtype=dtype, device=device))
    for g, K in Ks:                                                   # sum_i K_i^{dc} (x) conj(K_i^{dc}) in block (d, c)
        if isinstance(K, SectorJump):
            K = K.to_inter_sector().kraus.to(device=device, dtype=dtype)
        G = G + g * torch.einsum('brdcai,brdcel->bdaecil', K, K.conj()).reshape(K.shape[0], C, k * k, C, k * k)
    return G.reshape(G.shape[0], C * k * k, C * k * k)


def _exact_step(G: torch.Tensor, h: float, alg: TypeIAlgebra) -> Callable[[torch.Tensor], torch.Tensor]:
    """rho -> exp(h L) rho with the generator matrix G of _generator_matrix."""
    S = cost.batched_call(torch.linalg.matrix_exp, h * G, f"exp(hL) (size {G.shape[-1]})", kind='expm')
    mask = alg._const('block_mask', torch.bool).to(G.device)

    def step(rho):
        v = rho.reshape(rho.shape[0], -1, 1).to(S.dtype)
        return (S @ v).reshape(max(S.shape[0], rho.shape[0]), *rho.shape[1:]) * mask
    return step


LINDBLAD_METHODS = ('rk4', 'kraus', 'exact')


def lindblad_evolve(rho0: Operator, H: Optional[OperatorLike], jumps: Sequence[OperatorLike] = (),
                    times=(0.0, 1.0), rates: Optional[Sequence[float]] = None,
                    substeps: int = 20, progress: bool = False, method: str = 'rk4') -> List[Operator]:
    """
    Integrate the Lindblad equation of lindblad_rhs (jumps inside and between sectors); returns rho(t) for every t
    in `times`. Every interval of `times` is split into `substeps` steps of size h.

    method='rk4' (default): classical Runge-Kutta, fourth order, but a step is not a positive map in general, even
        for small h (Proposition 8.15 of the lecture notes): rho(t) can acquire small negative eigenvalues.
    method='kraus': the splitting of Theorem 8.16, exp(-i h/2 H_eff) . exp(-i h/2 H_eff)^* around
        id + h Phi + h^2/2 Phi^2 (Phi the jump part), renormalised to be exactly trace preserving. Every step is a
        channel (completely positive and trace preserving) for every h; second order (error O(t h^2)).
    method='exact': exp(h L) from the matrix of the generator on (+)_c M_{k_c} (size C k_max^2, so for small
        sectors; jumps between sectors included). For a constant generator the result is exact and one exponential
        is computed per interval length (`substeps` is not used).
    For time-dependent H or jumps, 'kraus' and 'exact' freeze the generator at the midpoint of each step; every
    step is still a channel. `times` must be nondecreasing for these two methods.
    """
    if method not in LINDBLAD_METHODS:
        raise ValueError(f"method must be one of {LINDBLAD_METHODS}, got {method!r}")
    rdt = rho0.matrix.dtype
    dtype = torch.promote_types(rdt, torch.complex64) if rho0.algebra.hilbert.complex_valued or H is not None else rdt
    if H is not None:
        _require_complex(rho0.algebra)
    # one initial state for a batch of Hamiltonians / jump operators (disorder averages), or vice versa
    t0 = float(_times(times, 'cpu')[0])

    def _batch(X):                                     # batch size of a (possibly time-dependent) H or jump
        if isinstance(X, InterSectorChannel):
            return X.kraus.shape[0]
        if isinstance(X, SectorJump):
            return X.val.shape[0]
        X = _mat(X(t0)) if callable(X) and not isinstance(X, Operator) else _mat(X)
        return X.shape[0] if isinstance(X, torch.Tensor) else 1
    batches = [rho0.matrix.shape[0]] + [_batch(X) for X in ([H] if H is not None else []) + list(jumps)]
    y0 = rho0.matrix.to(dtype)
    y0 = y0.expand(max(batches), *y0.shape[1:]).clone()
    alg = rho0.algebra
    if method == 'rk4':
        ys = rk4(lindblad_rhs(H, jumps, rates, algebra=alg), y0, times, substeps, progress)
    else:
        ys = _lindblad_splitting(alg, H, jumps, rates, y0, times, substeps, progress, method)
    ys = 0.5 * (ys + ys.conj().transpose(-2, -1))          # once for the whole trajectory
    if isinstance(rho0, DensityMatrix):
        return [DensityMatrix(alg, matrix=y, validate=False) for y in ys]
    return [alg.operator(y, is_self_adjoint=True) for y in ys]


def _lindblad_splitting(alg, H, jumps, rates, y0, times, substeps, progress, method) -> torch.Tensor:
    """The 'kraus' and 'exact' time steppers of lindblad_evolve: rho at every point of `times`."""
    plain, inter, timedep, algebra = _lindblad_terms(H, jumps, rates, alg)
    ts = _times(times, 'cpu').tolist()
    if any(b < a for a, b in zip(ts[:-1], ts[1:])):
        raise ValueError(f"lindblad_evolve(method={method!r}) needs nondecreasing times")
    dtype, dev = y0.dtype, y0.device
    cost.check_memory(len(ts) * y0.numel() * y0.element_size(), dev, f"Lindblad trajectory ({len(ts)} time points)")
    n_sub = 1 if (method == 'exact' and not timedep) else substeps
    if method == 'exact':
        n = alg.C * alg.k_max ** 2
        cost.check_memory(8 * cost.tensor_bytes((y0.shape[0], n, n), dtype), dev,
                          f"Lindblad generator on (+)_c M_k_max ({n} x {n})")
    out = y0.new_empty((len(ts), *y0.shape))
    out[0] = y0
    y = y0
    const_parts, const_G, steps = None, None, {}
    timer = cost.StepTimer(max(len(ts) - 1, 0) * n_sub, f"lindblad_evolve ({method})", progress)
    for i, (ta, tb) in enumerate(zip(ts[:-1], ts[1:])):
        h = (tb - ta) / n_sub
        for s_ in range(n_sub):
            if timedep or h not in steps:
                if timedep:
                    parts = _lindblad_parts(H, plain, inter, algebra, ta + (s_ + 0.5) * h, dtype, dev)   # midpoint
                else:
                    if const_parts is None:
                        const_parts = _lindblad_parts(H, plain, inter, algebra, ta, dtype, dev)
                    parts = const_parts
                if method == 'kraus':
                    step = _kraus_step(parts, h, alg, dtype, dev)
                else:
                    G = const_G if (const_G is not None and not timedep) else _generator_matrix(parts, alg, dtype, dev)
                    if not timedep:
                        const_G = G
                    step = _exact_step(G, h, alg)
                if not timedep:
                    steps[h] = step
            else:
                step = steps[h]
            y = step(y)
            timer.step(y.device)
        out[i + 1] = y
    timer.close()
    return out


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
