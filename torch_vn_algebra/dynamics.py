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
from . import cost, functional
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


def _packed(alg: Optional[TypeIAlgebra]) -> bool:
    return alg is not None and alg._packed


def _slabs_of(x: OperatorLike, alg: Optional[TypeIAlgebra]) -> List[torch.Tensor]:
    """An operator (or a raw padded tensor, masked to the blocks) as the slabs of `alg`: [padded tensor] or, in the
    packed layout, the group tensors (TypeIAlgebra._slab_groups). Without an algebra: [padded tensor]."""
    if alg is None:
        return [_mat(x)]
    if isinstance(x, Operator):
        return list(x.parts) if alg._packed else [x.matrix]
    x = _mat(x, alg)
    return alg._unpad(x) if alg._packed else [x]


def _y_of(slabs: List[torch.Tensor], alg: TypeIAlgebra) -> torch.Tensor:
    """Slabs as the one tensor the integrators work on: the padded tensor, or the flattened blocks when packed."""
    return alg._flatten(slabs) if alg._packed else slabs[0]


def _slabs_of_y(y: torch.Tensor, alg: TypeIAlgebra) -> List[torch.Tensor]:
    return alg._unflatten(y) if alg._packed else [y]


def _apply_vec(op: Operator, psi: torch.Tensor) -> torch.Tensor:
    """op psi for state vectors psi of shape (batch, C, k_max) (block layout), in either storage layout."""
    alg = op.algebra
    if not alg._packed:
        m = op.matrix
        return (m @ psi.to(m.dtype).unsqueeze(-1)).squeeze(-1)
    dtype = torch.promote_types(op._dtype(), psi.dtype)
    out = torch.zeros(max(op.batch_size, psi.shape[0]), *psi.shape[1:], dtype=dtype, device=psi.device)
    for (k, idx), p in zip(alg.size_groups, op.parts):
        sel = alg._sel(idx)
        out[:, sel, :k] = (p.to(dtype) @ psi[:, sel, :k, None].to(dtype)).squeeze(-1)
    return out


def _batch_of(x) -> int:
    if isinstance(x, Operator):
        return x.batch_size
    return x.shape[0] if isinstance(x, torch.Tensor) else 1


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
    psi0 = psi0.to(H._device())                           # a state built on the CPU may meet H on a GPU
    B = max(psi0.shape[0], H.batch_size)                  # one initial state for a batch of Hamiltonians, or vice versa
    psi0 = psi0.expand(B, *psi0.shape[1:])
    ts = _times(times, psi0.device)
    hdt = H._dtype()
    if psi0.dtype in (torch.float64, torch.complex128) and hdt in (torch.float32, torch.complex64):
        warnings.warn("schrodinger: psi0 is in double but H in single precision; the evolution is computed in single "
                      "precision (create H in a double-precision algebra)", UserWarning, stacklevel=2)
    alg = H.algebra
    Bh = H.batch_size
    if torch.is_grad_enabled() and H._requires_grad():           # differentiable in H: propagators (Daleckii-Krein)
        return torch.stack([_apply_vec(propagator(H, float(t)), psi0.to(hdt)) for t in ts])
    out = torch.zeros(len(ts), *psi0.shape, dtype=hdt, device=psi0.device)
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
    hdt = H._dtype()
    Bh = H.batch_size
    B = max(rho0.batch_size, Bh)                          # one initial state for a batch of Hamiltonians, or vice versa
    slabs = [s.to(hdt).expand(B, *s.shape[1:]) for s in _slabs_of(rho0, alg)]
    ts = _times(times, slabs[0].device)
    T = len(ts)
    nbytes = sum(s.numel() for s in slabs) * slabs[0].element_size()
    cost.check_memory(2 * T * nbytes, slabs[0].device, f"von_neumann ({T} time points)")
    differentiable = torch.is_grad_enabled() and H._requires_grad()
    if differentiable:                                            # differentiable in H: propagators (Daleckii-Krein)
        R0 = Operator(alg, matrix=None if alg._packed else slabs[0], parts=slabs if alg._packed else None)
        trajs = []
        for t in ts:
            U = propagator(H, float(t))
            trajs.append(_slabs_of((U @ R0 @ U.adjoint()), alg))
    else:
        trajs = [[torch.zeros_like(s) for s in slabs] for _ in range(T)] if not alg._packed else [[] for _ in range(T)]
    for g, (k, idx, _, (w, V)) in enumerate([] if differentiable else H._decompositions('hermitian')):
        n = len(idx)
        sel = alg._sel(idx)
        V = V.reshape(Bh, n, k, k)
        Vh = V.conj().transpose(-2, -1)
        rho_g = slabs[g] if alg._packed else slabs[0][:, sel, :k, :k]
        R = Vh @ rho_g @ V                                                # rho in the eigenbasis, (B, n, k, k)
        wd = w.reshape(Bh, n, k).double()
        # all times at once: rho(t)_ij = exp(-i t (w_i - w_j)) rho_ij in the eigenbasis (phases in double)
        phase = torch.exp(-1j * ts[:, None, None, None, None] * (wd[..., :, None] - wd[..., None, :])[None])
        Rt = V[None] @ (phase.to(R.dtype) * R[None]) @ Vh[None]           # (T, B, n, k, k)
        for i in range(T):
            if alg._packed:
                trajs[i].append(Rt[i])
            else:
                trajs[i][0][:, sel, :k, :k] = Rt[i]
    out = []
    for sl in trajs:
        data = dict(parts=sl) if alg._packed else dict(matrix=sl[0])
        if isinstance(rho0, DensityMatrix):
            out.append(DensityMatrix(alg, validate=False, **data))
        else:
            out.append(Operator(alg, is_self_adjoint=rho0._is_self_adjoint, is_positive=rho0._is_positive, **data))
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
    is_const = not (callable(H_of_t) and not isinstance(H_of_t, Operator))
    packed_op = isinstance(H_of_t, Operator) and H_of_t.algebra._packed
    const = None if not is_const or packed_op else _mat(H_of_t, algebra)
    if is_const:
        psi0 = psi0.to(H_of_t._device() if isinstance(H_of_t, Operator) else const.device)
    if algebra is not None and min(algebra.k_factors) != algebra.k_max:
        psi0 = psi0 * algebra._const('rows').to(psi0.device)        # entries in the padding are not part of psi
    y0 = psi0.to(torch.promote_types(psi0.dtype, torch.complex64))
    Hc = None if const is None else (-1j * const.to(y0.dtype))      # hoisted: constant H is cast once

    def f(t, psi):
        if Hc is not None:
            return (Hc @ psi.unsqueeze(-1)).squeeze(-1)
        Ht = H_of_t if packed_op else H_of_t(t)
        if isinstance(Ht, Operator) and Ht.algebra._packed:          # blocks of each size, no padding
            return -1j * _apply_vec(Ht, psi).to(psi.dtype)
        return ((-1j * _mat(Ht, algebra).to(psi.dtype)) @ psi.unsqueeze(-1)).squeeze(-1)
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
    """(H, A, stacked jumps with their rates, jumps between sectors) at time t, as slabs of `algebra` (lists with one
    tensor per slab: (B, n, k, k) for H and A, (B, J, n, k, k) for the stacked jumps); with
    A = sum_j g_j L_j^* L_j, drho/dt = -i (H rho - rho H^*) - 1/2 (A rho + rho A) + sum_j g_j L_j rho L_j^*
    + (between sectors). A is kept apart from H so that a real algebra with real jumps stays real."""
    def get(x):
        x = x(t) if callable(x) and not isinstance(x, Operator) else x
        return [s.to(device=device, dtype=dtype) for s in _slabs_of(x, algebra)]

    Hm = get(H) if H is not None else None
    A = None
    Ls, gs = None, None
    if plain:
        per = [get(L) for _, L in plain]                                  # per jump, per slab
        Ls = [torch.stack(torch.broadcast_tensors(*[p_[s] for p_ in per]), dim=1) for s in range(len(per[0]))]
        gs = torch.tensor([g for g, _ in plain], dtype=Ls[0].real.dtype, device=device)
        A = [torch.einsum('j,bjcai,bjcak->bcik', gs.to(dtype), L.conj(), L) for L in Ls]
    Ks = []
    for g, L in inter:
        if isinstance(L, SectorJump):                                     # sparse: L^*L diagonal, never formed densely
            LdL = _slabs_of(_jump_operator(L.LdL(), L.algebra), algebra)
            LdL = [x.to(device=device, dtype=dtype) for x in LdL]
            A = [g * x for x in LdL] if A is None else [a + g * x for a, x in zip(A, LdL)]
            Ks.append((g, L))
            continue
        K = L.kraus.to(device=device, dtype=dtype)                        # (B, r, D, C, m, k)
        KdK = _slabs_of((K.conj().transpose(-2, -1) @ K).sum(dim=(1, 2)), algebra)   # block diagonal
        A = [g * x for x in KdK] if A is None else [a + g * x for a, x in zip(A, KdK)]
        Ks.append((g, K))
    return Hm, A, Ls, gs, Ks


def _jump_operator(data, alg: TypeIAlgebra) -> Operator:
    """Stored data of `alg` (a padded tensor or a list of group tensors) as an Operator (no copy)."""
    return Operator(alg, parts=data) if isinstance(data, list) else Operator(alg, matrix=data)


def _sector_jump(K: SectorJump, slabs, alg: TypeIAlgebra, dual: bool = False):
    """K rho K^* (or K^* Z K) on slabs of `alg`; the jump may come from an algebra of the other layout."""
    f = K.dual_sandwich if dual else K.sandwich
    if alg is not None and K.algebra._packed == alg._packed:
        return alg._slabs(f(alg._unslab(slabs)))
    pad = alg._pad(slabs) if _packed(alg) else slabs[0]                  # through the padded tensor
    res = f(pad)                                                         # (a packed jump pads and unpads itself)
    return alg._unpad(res) if _packed(alg) else [res]


def _kraus_between(K: torch.Tensor, slabs, alg: TypeIAlgebra, dual: bool = False):
    """sum_{i,c} K^{dc} rho_c K^{dc*} (or the dual) for padded Kraus tensors (B, r, D, C, m, k), on slabs."""
    rho = alg._pad(slabs) if _packed(alg) else slabs[0]
    if dual:
        X = (K.conj().transpose(-2, -1) @ rho[:, None, :, None] @ K).sum(dim=(1, 2))   # (B, C, k, k)
    else:
        B_, r, D, C, m, k = K.shape
        X = K @ rho[:, None, None]                                          # (B, r, D, C, m, k)
        Bx = X.shape[0]
        Xc = X.permute(0, 2, 4, 1, 3, 5).reshape(Bx, D, m, r * C * k)
        Kc = K.permute(0, 2, 4, 1, 3, 5).reshape(B_, D, m, r * C * k)
        X = Xc @ Kc.conj().transpose(-2, -1)
    return alg._unpad(X) if _packed(alg) else [X]


def _add_gain(out, rho, Ls, gs, Ks, alg: TypeIAlgebra):
    """out + Phi(rho) on slabs, Phi(rho) = sum_j g_j L_j rho L_j^* + sum_m g_m Phi_m(rho) (the jump part of the
    generator)."""
    out = list(out)
    if Ls is not None:                     # sum_j g_j L_j rho L_j^*: one product over the merged (j, column) index
        for s, (L, r) in enumerate(zip(Ls, rho)):
            X = (L @ r.unsqueeze(1)) * gs[:, None, None, None].to(r.dtype)              # (B, J, n, k, k)
            Bx, J, C, k = X.shape[0], X.shape[1], X.shape[2], X.shape[-1]
            Xc = X.permute(0, 2, 3, 1, 4).reshape(Bx, C, k, J * k)
            Lc = L.permute(0, 2, 3, 1, 4).reshape(L.shape[0], C, k, J * k)
            out[s] = out[s] + Xc @ Lc.conj().transpose(-2, -1)
    for g, K in Ks:                        # gain between sectors: sum_{i, c} K^{dc} rho_c K^{dc*}
        if isinstance(K, SectorJump):      # sparse jump: a weighted permutation of the entries, O(sum_c k_c^2)
            Y = _sector_jump(K, rho, alg)
        else:
            Y = _kraus_between(K, rho, alg)
        out = [o + g * y for o, y in zip(out, Y)]
    return out


def _gain_dual(Z, Ls, gs, Ks, alg: TypeIAlgebra):
    """Phi^*(Z) = sum_j g_j L_j^* Z L_j + sum_m g_m Phi_m^*(Z) on slabs, the dual of the jump part for the blunt
    trace (Phi_m^*(Z)_c = sum_{d,i} K_i^{dc*} Z_d K_i^{dc})."""
    out = [torch.zeros_like(z) for z in Z]
    if Ls is not None:
        for s, (L, z) in enumerate(zip(Ls, Z)):
            X = L.conj().transpose(-2, -1) @ z.unsqueeze(1) @ L                  # (B, J, n, k, k)
            out[s] = out[s] + (X * gs[:, None, None, None].to(z.dtype)).sum(dim=1)
    for g, K in Ks:
        Y = _sector_jump(K, Z, alg, dual=True) if isinstance(K, SectorJump) else _kraus_between(K, Z, alg, dual=True)
        out = [o + g * y for o, y in zip(out, Y)]
    return out


def _rhs_slabs(H, jumps, rates, algebra):
    """The GKSL right-hand side (t, slabs) -> slabs of lindblad_rhs, in the layout of `algebra`."""
    plain, inter, timedep, algebra = _lindblad_terms(H, jumps, rates, algebra)
    cache = {}

    def f(t, rho):
        key = (rho[0].dtype, rho[0].device)
        if timedep or key not in cache:
            cache[key] = _lindblad_parts(H, plain, inter, algebra, t, rho[0].dtype, rho[0].device)
        Hm, A, Ls, gs, Ks = cache[key]
        if Hm is None:
            out = [torch.zeros_like(r) for r in rho]
        else:
            out = [-1j * (h @ r - r @ h.conj().transpose(-2, -1)) for h, r in zip(Hm, rho)]
        if A is not None:
            out = [o - 0.5 * (a @ r + r @ a) for o, a, r in zip(out, A, rho)]
        return _add_gain(out, rho, Ls, gs, Ks, algebra)
    return f, algebra


def lindblad_rhs(H: Optional[OperatorLike], jumps: Sequence[OperatorLike] = (),
                 rates: Optional[Sequence[float]] = None,
                 algebra: Optional[TypeIAlgebra] = None) -> Callable[[float, torch.Tensor], torch.Tensor]:
    """
    Right-hand side of the GKSL equation
        drho/dt = -i[H, rho] + sum_j g_j (L_j rho L_j^* - 1/2 {L_j^* L_j, rho})
    as a function (t, rho_matrix) -> drho/dt on padded tensors (batch, C, k_max, k_max). H and L_j are constant, or
    callables of t.

    A jump that moves weight between sectors (particle loss N -> N-1, ...) is passed as an
    InterSectorChannel Phi of the algebra into itself (or, for spin flips and other jumps that map basis vectors to
    multiples of basis vectors, as a sparse SectorJump); it contributes Phi(rho) - 1/2 {Phi^*(1), rho}, where
    Phi^*(1)_c = sum_{d,i} K_i^{dc*} K_i^{dc} is block diagonal. This is the GKSL generator of the jumps
    K_i^{dc} compressed to M (the conditional expectation onto M applied to the generator on B(H)); it maps M
    into M and conserves the total trace. For constant H and jumps, A = sum g L^* L and the stacked jumps are
    prepared once; callables of t are evaluated at every step. Raw tensors are restricted to the blocks of
    `algebra` (taken from H or the jumps when they are Operators); without an algebra their padding must be zero.
    In the packed layout the computation runs on the blocks of each size (lindblad_evolve integrates them directly).
    """
    g, algebra = _rhs_slabs(H, jumps, rates, algebra)
    if not _packed(algebra):
        return lambda t, rho: g(t, [rho])[0]
    return lambda t, rho: algebra._pad(g(t, algebra._unpad(rho)))


def _kraus_step(parts, h: float, alg: TypeIAlgebra, dtype, device) -> Callable:
    """
    One step of the completely positive, trace-preserving splitting of Theorem 8.16 for the frozen generator
    `parts`: with E = exp(-i h/2 H_eff), H_eff = H - i/2 A, S(X) = E X E^*, and Phi the jump part,
        Psi_h = S o (id + h Phi + h^2/2 Phi^2) o S,    M = Psi_h^*(1),    rho -> Psi_h(M^{-1/2} rho M^{-1/2}).
    Every factor is CP and the last one makes the step exactly trace preserving. In floating point the trace is
    preserved up to about eps * cond(M); M is ill-conditioned when h times the spread of the decay rates (the
    eigenvalues of A) is large, and then the step is taken as two steps of size h/2 (recursively), each a channel.
    The step maps slabs to slabs.
    """
    Hm, A, Ls, gs, Ks = parts
    eyes = [e.to(device=device, dtype=dtype) for e in alg._slab_eyes(dtype)]
    masks = [None] * len(eyes) if alg._packed else [alg._const('block_mask', torch.bool).to(device)]
    Es = []
    for s, (eye, mask) in enumerate(zip(eyes, masks)):
        G = None                                                      # -i h/2 H_eff = -i h/2 H - h/4 A
        if Hm is not None:
            G = (-0.5j * h) * Hm[s]
        if A is not None:
            G = -0.25 * h * A[s] if G is None else G - 0.25 * h * A[s]
        if G is None:
            Es.append(eye)
            continue
        k = G.shape[-1]
        E = cost.batched_call(torch.linalg.matrix_exp, G.reshape(-1, k, k), "exp(-i h/2 H_eff)",
                              kind='expm').reshape(G.shape)
        Es.append(E if mask is None else E * mask)                    # the padding of exp(G) is the identity
    jumps = Ls is not None or len(Ks) > 0

    def phi2(Y):                                                      # Y + h Phi(Y) + h^2/2 Phi(Phi(Y))
        if not jumps:
            return Y
        Y1 = _add_gain([torch.zeros_like(y) for y in Y], Y, Ls, gs, Ks, alg)
        Y2 = _add_gain([torch.zeros_like(y) for y in Y1], Y1, Ls, gs, Ks, alg)
        return [y + h * y1 + (0.5 * h * h) * y2 for y, y1, y2 in zip(Y, Y1, Y2)]

    def phi2_dual(W):
        if not jumps:
            return W
        W1 = _gain_dual(W, Ls, gs, Ks, alg)
        W2 = _gain_dual(W1, Ls, gs, Ks, alg)
        return [w + h * w1 + (0.5 * h * h) * w2 for w, w1, w2 in zip(W, W1, W2)]

    Ehs = [E.conj().transpose(-2, -1) for E in Es]
    D = phi2_dual([Eh @ E for E, Eh in zip(Es, Ehs)])
    Ns, ratios = [], []
    decs = []
    for s, (E, Eh, eye) in enumerate(zip(Es, Ehs, eyes)):
        M = Eh @ D[s] @ E                                             # Psi_h^*(1) >= (E^2)^* E^2 > 0
        M = 0.5 * (M + M.conj().transpose(-2, -1))
        pad = torch.eye(M.shape[-1], dtype=dtype, device=device) - eye   # identity on the padding, for eigh
        with torch.no_grad():
            w, V = torch.linalg.eigh((M + pad).detach())
        decs.append((M + pad, w, V))
        # condition of M, measured against 1 = Psi_0^*(1) as well: a block whose M has underflowed towards 0 must be
        # split too (in the padded layout the eigenvalues 1 of the padding play this part)
        top = w[..., -1] if not alg._packed else w[..., -1].clamp(min=1.0)
        ratios.append((top / w[..., 0].clamp(min=torch.finfo(w.dtype).tiny)).max())
    if h > 0 and torch.stack(ratios).max().item() > torch.finfo(decs[0][1].dtype).eps ** -0.25:
        half = _kraus_step(parts, 0.5 * h, alg, dtype, device)
        return lambda rho: half(half(rho))
    tiny = torch.finfo(decs[0][1].dtype).tiny                         # M^{-1/2}, block diagonal; Daleckii-Krein
    for (Mp, w, V), mask in zip(decs, masks):
        N = functional.spectral_apply(Mp, lambda x: x.clamp(min=tiny).rsqrt(),   # backward: finite gradients at
                                      decomposition=(w, V)).to(dtype)           # the degenerate spectra of M
        Ns.append(N if mask is None else N * mask)

    def step(rho):
        Y = [E @ (N @ r @ N) @ Eh for E, N, Eh, r in zip(Es, Ns, Ehs, rho)]
        return [E @ y @ Eh for E, y, Eh in zip(Es, phi2(Y), Ehs)]
    return step


def _packed_index(alg: TypeIAlgebra) -> torch.Tensor:
    """Positions of the flattened blocks (TypeIAlgebra._flatten, packed) inside the flattened padded tensor."""
    K = alg.k_max
    out = []
    for k, idx in alg.size_groups:
        a = torch.arange(k, device=alg.hilbert.device)
        rc = (a[:, None] * K + a[None, :]).reshape(-1)
        out += [c * K * K + rc for c in idx]
    return torch.cat(out)


def _generator_matrix(parts, alg: TypeIAlgebra, dtype, device) -> torch.Tensor:
    """Matrix of the generator on the flattened blocks (row-major vec of every block, the sectors one after the
    other: (+)_c M_{k_max} zero padded, or the packed blocks without padding), shape (batch, n, n); jumps between
    sectors give off-diagonal sector blocks."""
    Hm, A, Ls, gs, Ks = parts
    if alg._packed:                                                   # built on the padded blocks, then restricted
        Hm = None if Hm is None else alg._pad(Hm)
        A = None if A is None else alg._pad(A)
        Ls = None if Ls is None else torch.stack([alg._pad([L[:, j] for L in Ls]) for j in range(Ls[0].shape[1])],
                                                 dim=1)
    else:
        Hm = None if Hm is None else Hm[0]
        A = None if A is None else A[0]
        Ls = None if Ls is None else Ls[0]
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
    G = G.reshape(G.shape[0], C * k * k, C * k * k)
    if alg._packed:
        I = _packed_index(alg).to(device)
        G = G[:, I][:, :, I]
    return G


def _exact_step(G: torch.Tensor, h: float, alg: TypeIAlgebra) -> Callable[[torch.Tensor], torch.Tensor]:
    """y -> exp(h L) y with the generator matrix G of _generator_matrix (y: the padded tensor, or the flattened
    blocks when packed)."""
    S = cost.batched_call(torch.linalg.matrix_exp, h * G, f"exp(hL) (size {G.shape[-1]})", kind='expm')
    if alg._packed:
        return lambda y: (S @ y.to(S.dtype).unsqueeze(-1)).squeeze(-1)
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
    step is still a channel. `times` must be nondecreasing for these two methods. In the packed layout every method
    works on the blocks of each size (the state is never padded).
    """
    if method not in LINDBLAD_METHODS:
        raise ValueError(f"method must be one of {LINDBLAD_METHODS}, got {method!r}")
    rdt = rho0._dtype()
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
        return _batch_of(X(t0) if callable(X) and not isinstance(X, Operator) else X)
    batches = [rho0.batch_size] + [_batch(X) for X in ([H] if H is not None else []) + list(jumps)]
    alg = rho0.algebra
    y0 = _y_of([s.to(dtype) for s in _slabs_of(rho0, alg)], alg)
    y0 = y0.expand(max(batches), *y0.shape[1:]).clone()
    if method == 'rk4':
        g, _ = _rhs_slabs(H, jumps, rates, alg)
        ys = rk4(lambda t, y: _y_of(g(t, _slabs_of_y(y, alg)), alg), y0, times, substeps, progress)
    else:
        ys = _lindblad_splitting(alg, H, jumps, rates, y0, times, substeps, progress, method)
    out = []
    for y in ys:
        sl = [0.5 * (s + s.conj().transpose(-2, -1)) for s in _slabs_of_y(y, alg)]
        data = dict(parts=sl) if alg._packed else dict(matrix=sl[0])
        out.append(DensityMatrix(alg, validate=False, **data) if isinstance(rho0, DensityMatrix)
                   else Operator(alg, is_self_adjoint=True, **data))
    return out


def _lindblad_splitting(alg, H, jumps, rates, y0, times, substeps, progress, method) -> torch.Tensor:
    """The 'kraus' and 'exact' time steppers of lindblad_evolve: y (the padded tensor, or the flattened blocks when
    packed) at every point of `times`."""
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
                    kstep = _kraus_step(parts, h, alg, dtype, dev)
                    step = (lambda ks: lambda y_: _y_of(ks(_slabs_of_y(y_, alg)), alg))(kstep)
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
        if not isinstance(L, (Operator, torch.Tensor)):
            raise TypeError(f"jump operators must be Operators or tensors of shape (batch, C, k_max, k_max), "
                            f"got {type(L).__name__}")
    mats = ([H] if H is not None else []) + list(jumps)
    if not mats:
        raise ValueError("need a Hamiltonian or at least one jump operator")
    dtype = torch.float64
    for m in mats:
        dtype = torch.promote_types(dtype, m._dtype() if isinstance(m, Operator) else m.dtype)
    if H is not None:
        dtype = torch.complex128
    dev = mats[0]._device() if isinstance(mats[0], Operator) else mats[0].device
    batch = max(_batch_of(m) for m in mats)

    def group(x, g, k, sel):                          # the blocks of size k of x (no padded tensor for Operators)
        if isinstance(x, Operator):
            return x.parts[g] if alg._packed else x.matrix[:, sel, :k, :k]
        return x[:, sel, :k, :k]
    # peak: the generator, a few Kronecker products and the exponential / Choi / eigenvector workspaces of
    # lindblad_channel, all of size k_c^4 per sector and sample
    cost.check_memory(12 * sum(cost.tensor_bytes((batch, k * k, k * k), dtype) for k in alg.k_factors),
                      dev, f"Lindblad superoperator / exp(tL) (batch {batch}, k_c^2 x k_c^2 per sector)")
    out = []
    for gi, (k, idx) in enumerate(alg._size_groups()):
        sel = alg._sel(idx)
        n = len(idx)
        eye = torch.eye(k, dtype=dtype, device=dev)
        G = torch.zeros(batch, n, k * k, k * k, dtype=dtype, device=dev)
        A = None
        if jumps:
            Ls = torch.stack(torch.broadcast_tensors(*[group(L, gi, k, sel).to(device=dev, dtype=dtype)
                                                       for L in jumps]), dim=1)      # (B, J, n, k, k)
            g = torch.tensor(rates, dtype=dtype, device=dev)
            # sum_j g_j L_j (x) conj(L_j) in one contraction, and A = sum_j g_j L_j^* L_j
            G = G + torch.einsum('j,bjnai,bjnck->bnacik', g, Ls, Ls.conj()).reshape(Ls.shape[0], n, k * k, k * k)
            A = torch.einsum('j,bjnai,bjnak->bnik', g, Ls.conj(), Ls)
        if H is not None:
            Hc = group(H, gi, k, sel).to(device=dev, dtype=dtype)
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
