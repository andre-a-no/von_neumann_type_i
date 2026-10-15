"""
Permutation-invariant ensembles of N two-level systems: the Dicke-basis algebra, collective operators, and local
(single-site) dissipators as maps between spin sectors.

Mathematics
-----------
By Schur-Weyl duality the symmetric group S_N, permuting the factors of H = (C^2)^{(x) N}, and the collective
SU(2) generate each other's commutants, and

    H  =  (+)_j  V_j (x) C^{m_j},      dim V_j = 2j + 1,   j = N/2, N/2 - 1, ..., >= 0,

with m_j = binom(N, N/2 - j) - binom(N, N/2 - j - 1) (`su2.multiplicity`). The permutation-invariant operators are

    M  =  pi(S_N)'  =  (+)_j  M_{2j+1}(C) (x) 1_{m_j},

a Type I algebra with one block per total spin j; the trace of H restricts to Tr_H(A) = sum_j m_j tr(A_j). The
collective operators J_a = sum_i s_a^(i) (s = sigma / 2) act as the spin-j matrices in block j (standard basis
|j, m>, m = j, j-1, ..., -j, Condon-Shortley phases). A permutation-invariant state is rho = (+)_j sigma_j (x)
1_{m_j} / m_j; the block sigma_j = tr_{C^{m_j}}(block j of rho) is its density with respect to the blunt trace of
the library (sum_j tr sigma_j = 1), and Tr_H(rho A) = sum_j tr(sigma_j A_j). This is the "Dicke basis" density of
PIQS (Shammah et al. 2018).

Local dissipators. For a single-site operator L (sigma_-, sigma_+ or sigma_z) the map Phi(X) = sum_i L_i X L_i^*
commutes with the permutations and so maps M into M; as a normal completely positive map between Type I algebras
it has the form Phi(sigma)_{j'} = sum_{j, r} K_r^{j'j} sigma_j K_r^{j'j *} on blunt densities (see
`channels.InterSectorChannel`), with j' - j in {-1, 0, 1}. Restricting to the last site and branching N -> N - 1
(the multiplicity space of j for N sites is the direct sum of those of j1 = j +- 1/2 for N - 1 sites, and L_N
preserves j1 and the label in it; the other sites contribute equally by symmetry),

    K_{j1}^{j'j}  =  sqrt( N m^{(N-1)}_{j1} / m^{(N)}_j )  < j1, 1/2; j' | L_N | j1, 1/2; j >,

the last factor being the (2j'+1) x (2j+1) matrix of L acting on the second factor of V_{j1} (x) C^2 between the
coupled multiplets j and j' (Clebsch-Gordan coefficients of j1 (x) 1/2, Condon-Shortley phases), for j1 in
{j +- 1/2} cap {j' +- 1/2}. Each K has a single non-zero diagonal (L changes m by s = -1, +1, 0). The no-jump part
is fixed by Phi^*(1) = sum_i L_i^* L_i, which is diagonal: N/2 + J_z, N/2 - J_z and N for sigma_-, sigma_+,
sigma_z. Closed forms of the resulting coefficients are in Chase and Geremia, PRA 78, 052101 (2008), and Shammah
et al., PRA 98, 063815 (2018); the generator here agrees with the one of qutip.piqs (tests/test_permutational.py)
and with the exact Lindblad generator on (C^2)^{(x) N}.

Two forms of the same generator
-------------------------------
* `PermutationalSystem.local_jump(kind)` gives Phi as an `InterSectorChannel` (Kraus rank 2: the two branches
  j1), so the generic `dynamics.lindblad_evolve` / `lindblad_rhs` handle collective and local dissipation together
  (jumps inside sectors as Operators, between sectors as InterSectorChannel). Its dense Kraus tensor has
  2 C^2 k_max^2 entries (C = number of sectors, k_max = N + 1): practical for N up to a few tens.
* `PermutationalLindbladian` is the same generator for H = omega_x J_x + omega_y J_y + omega_z J_z, collective
  D[J_-], D[J_+], D[J_z] and local D[sigma_-], D[sigma_+], D[sigma_z], with every coefficient batched (one
  parameter set per batch element). All its operators are banded in the basis |j, m>, so one evaluation costs
  O(B C k_max^2) = O(B N^3) operations instead of O(B N^4) for dense block products, and no inter-sector tensor
  is stored. `rhs` plugs into `dynamics.rk4`; `evolve` integrates and returns DensityMatrix objects. Without
  omega_x and omega_z the generator is real and runs in a real algebra (complex_valued=False), at half the memory.

Truncation: PermutationalSystem(N, j_min=...) keeps the sectors j >= j_min only (e.g. j_min = N/2 for purely
collective dynamics, a single block of size N + 1); see the class docstring for the error bound.

    ps = PermutationalSystem(N=40)
    L = PermutationalLindbladian(ps, omega_x=2.0, collective_emission=1 / 20, emission=[0.0, 0.05, 0.1])
    times = torch.linspace(0, 20, 201, dtype=torch.float64)
    rho = L.evolve(ps.dicke_state(20, -20), times)          # batch of 3 trajectories, rho[t] a DensityMatrix
    jz = torch.stack([r.expectation(ps.jz()) for r in rho])  # (201, 3)

References: R. H. Dicke, Phys. Rev. 93, 99 (1954); B. A. Chase and J. M. Geremia, PRA 78, 052101 (2008);
B. Q. Baragiola, B. A. Chase and J. Geremia, PRA 81, 032104 (2010); N. Shammah, S. Ahmed, N. Lambert,
S. De Liberato and F. Nori, PRA 98, 063815 (2018).
"""
import math
from fractions import Fraction
from typing import Dict, List, Sequence, Tuple, Union

import torch

from .algebra import TypeIAlgebra
from .channels import InterSectorChannel
from .states import DensityMatrix
from .su2 import multiplicity
from . import cost

Operator = TypeIAlgebra.Operator
Number = Union[float, Sequence[float], torch.Tensor]

LOCAL_KINDS = ('-', '+', 'z')
_SHIFT = {'-': -1, '+': 1, 'z': 0}               # L changes m by this amount


def _branch_coefficients(two_j1: int):
    """(u, d) for the coupled multiplets J = j1 +- 1/2 of V_{j1} (x) C^2 (Condon-Shortley):
    |J M> = u(J, M) |j1, M - 1/2>|up> + d(J, M) |j1, M + 1/2>|down>. Returns a function (two_J, two_M) -> (u, d)."""
    den = two_j1 + 1

    def coeff(two_J, two_M):
        a = math.sqrt(max(two_j1 + two_M + 1, 0) / (2 * den))      # sqrt((j1 + M + 1/2) / (2 j1 + 1))
        b = math.sqrt(max(two_j1 - two_M + 1, 0) / (2 * den))      # sqrt((j1 - M + 1/2) / (2 j1 + 1))
        if two_J == two_j1 + 1:
            return a, b
        return -b, a
    return coeff


def _site_element(kind: str, two_j1: int, two_Jout: int, two_Jin: int, two_M: int) -> float:
    """<j1, 1/2; J_out, M + s | L | j1, 1/2; J_in, M> for L = sigma_-, sigma_+, sigma_z on the second factor."""
    c = _branch_coefficients(two_j1)
    u_in, d_in = c(two_Jin, two_M)
    s = _SHIFT[kind]
    u_out, d_out = c(two_Jout, two_M + 2 * s)
    if kind == '-':                       # |up> -> |down>: the u-component of the input, d-component of the output
        return d_out * u_in
    if kind == '+':
        return u_out * d_in
    return u_out * u_in - d_out * d_in


class PermutationalSystem:
    """
    N two-level systems with permutation symmetry: the algebra M = (+)_j M_{2j+1} (sectors ordered by increasing
    j, charges 2j, basis m = j, j-1, ..., -j in every block), its multiplicities m_j (the trace weights of Tr_H),
    collective spin operators, Dicke and coherent spin states, and the local jump maps between sectors.

    j_min (default: the smallest total spin) keeps only the sectors j >= j_min, i.e. the algebra P M P with P the
    central projection onto them. Local jumps into a removed sector are then dropped from the gain but kept in the
    loss, so the truncated generator is the compression of the full one: it generates a trace-decreasing CP
    semigroup with 0 <= P rho_trunc(t) P <= P rho(t) P, and for every A in P M the error is at most
    ||A|| (2 (1 - Tr rho_trunc(t))), the lost weight being observable during the run.
    """

    def __init__(self, N: int, precision: str = 'double', device=None, complex_valued: bool = True,
                 j_min=None):
        if N < 1:
            raise ValueError("N must be at least 1")
        self.N = N
        t_min = N % 2 if j_min is None else int(round(2 * float(j_min)))
        if t_min < N % 2 or t_min > N or (N - t_min) % 2:
            raise ValueError(f"j_min = {j_min} is not a total spin of {N} spins 1/2")
        self.two_j = list(range(t_min, N + 1, 2))                    # 2j, increasing
        self.spins = [Fraction(t, 2) for t in self.two_j]
        self.multiplicities = [multiplicity(N, Fraction(t, 2)) for t in self.two_j]
        k = [t + 1 for t in self.two_j]
        self.algebra = TypeIAlgebra(k, k, precision=precision, device=device, charges=list(self.two_j),
                                    complex_valued=complex_valued)
        self.C, self.k_max = self.algebra.C, self.algebra.k_max
        self.device = self.algebra.hilbert.device
        self.dtype = self.algebra.hilbert.dtype
        self.real_dtype = self.algebra.hilbert.real_dtype
        self._kraus_cache = {}

    def __repr__(self) -> str:
        return f"PermutationalSystem(N={self.N}, sectors j={self.two_j[0] / 2}..{self.N / 2}, k_max={self.k_max})"

    # ------------------------------------------------------------------ structure
    def sector_index(self, j) -> int:
        t = int(round(2 * float(j)))
        if t not in self.two_j or abs(2 * float(j) - t) > 1e-9:
            raise ValueError(f"j = {j} is not a total spin of {self.N} spins 1/2")
        return self.two_j.index(t)

    def trace_weights(self) -> torch.Tensor:
        """m_j as floats (C,): Tr_H(A) = sum_j m_j tr(A_j) for A in M (large N: float64 overflows beyond N ~ 1000)."""
        return torch.tensor([float(m) for m in self.multiplicities], dtype=torch.float64)

    def m_values(self) -> torch.Tensor:
        """(C, k_max) magnetic quantum numbers m = j - a of the basis vectors a = 0..2j (0 in the padding)."""
        a = torch.arange(self.k_max, dtype=torch.float64)
        j = torch.tensor([t / 2 for t in self.two_j], dtype=torch.float64)
        m = j[:, None] - a[None, :]
        return torch.where(a[None, :] <= 2 * j[:, None], m, torch.zeros_like(m))

    def raising_coefficients(self) -> torch.Tensor:
        """p (C, k_max): <a-1| J_+ |a> = sqrt(j(j+1) - m(m+1)), m = j - a, for 1 <= a <= 2j (else 0)."""
        a = torch.arange(self.k_max, dtype=torch.float64)
        j = torch.tensor([t / 2 for t in self.two_j], dtype=torch.float64)[:, None]
        m = j - a[None, :]
        p = (j * (j + 1) - m * (m + 1)).clamp(min=0).sqrt()
        return torch.where((a[None, :] >= 1) & (a[None, :] <= 2 * j), p, torch.zeros_like(p))

    # ------------------------------------------------------------------ collective operators
    def _op(self, mat: torch.Tensor, **flags) -> Operator:
        return self.algebra.operator(mat.to(device=self.device, dtype=self.dtype).unsqueeze(0), **flags)

    def jz(self) -> Operator:
        return self._op(torch.diag_embed(self.m_values()), is_self_adjoint=True)

    def jplus(self) -> Operator:
        p = self.raising_coefficients()
        return self._op(torch.diag_embed(p[:, 1:], offset=1))

    def jminus(self) -> Operator:
        p = self.raising_coefficients()
        return self._op(torch.diag_embed(p[:, 1:], offset=-1))

    def jx(self) -> Operator:
        p = self.raising_coefficients()
        return self._op(0.5 * (torch.diag_embed(p[:, 1:], offset=1) + torch.diag_embed(p[:, 1:], offset=-1)),
                        is_self_adjoint=True)

    def jy(self) -> Operator:
        if not self.algebra.hilbert.complex_valued:
            raise ValueError("J_y is not a real matrix: use a complex algebra (complex_valued=True)")
        p = self.raising_coefficients().to(torch.complex128)
        return self._op(-0.5j * (torch.diag_embed(p[:, 1:], offset=1) - torch.diag_embed(p[:, 1:], offset=-1)),
                        is_self_adjoint=True)

    def j2(self) -> Operator:
        """Casimir J^2 = (+)_j j(j+1) 1."""
        return self.algebra.central([t / 2 * (t / 2 + 1) for t in self.two_j])

    # ------------------------------------------------------------------ states
    def dicke_state(self, j, m) -> DensityMatrix:
        """|j, m><j, m| (x) 1_{m_j} / m_j: the permutation-invariant state of the Dicke vector in block j."""
        c = self.sector_index(j)
        a = int(round(float(j) - float(m)))
        if not 0 <= a <= self.two_j[c]:
            raise ValueError(f"m = {m} is not in -j..j for j = {j}")
        mat = torch.zeros(1, self.C, self.k_max, self.k_max, dtype=self.dtype, device=self.device)
        mat[0, c, a, a] = 1
        return DensityMatrix(self.algebra, matrix=mat, validate=False)

    def coherent_state(self, theta: Number, phi: Number = 0.0) -> DensityMatrix:
        """
        Coherent spin state (cos(theta/2)|up> + e^{i phi} sin(theta/2)|down>)^{(x) N} (theta = 0: all up), in the
        block j = N/2; theta, phi may be tensors (batch,). <J> = N/2 (sin theta cos phi, sin theta sin phi, cos theta).
        """
        th = torch.as_tensor(theta, dtype=torch.float64).reshape(-1)
        ph = torch.as_tensor(phi, dtype=torch.float64).reshape(-1)
        th, ph = torch.broadcast_tensors(th, ph)
        N = self.N
        a = torch.arange(N + 1, dtype=torch.float64)                  # a = N/2 - m = number of down spins
        logb = torch.tensor([math.lgamma(N + 1) - math.lgamma(k + 1) - math.lgamma(N - k + 1) for k in range(N + 1)],
                            dtype=torch.float64)
        c, s = torch.cos(th / 2)[:, None], torch.sin(th / 2)[:, None]
        amp = torch.exp(0.5 * logb)[None] * torch.where(a[None] < N, c.abs() ** (N - a)[None], torch.ones_like(c)) \
            * torch.where(a[None] > 0, s.abs() ** a[None], torch.ones_like(s))
        amp = amp * torch.sign(c) ** (N - a)[None] * torch.sign(s) ** a[None]
        psi = amp.to(torch.complex128) * torch.exp(1j * ph[:, None] * a[None])
        block = psi[:, :, None] * psi[:, None, :].conj()
        if not self.algebra.hilbert.complex_valued:
            if bool((block.imag.abs() > 1e-12).any()):
                raise ValueError("a real algebra holds coherent states with phi = 0 or pi only")
            block = block.real
        mat = torch.zeros(th.shape[0], self.C, self.k_max, self.k_max, dtype=self.dtype, device=self.device)
        mat[:, -1, :N + 1, :N + 1] = block.to(device=self.device, dtype=self.dtype)
        return DensityMatrix(self.algebra, matrix=mat, validate=False)

    def sector_probabilities(self, rho: Operator) -> torch.Tensor:
        """p_j = tr sigma_j (batch, C): the distribution of the total spin j."""
        return torch.diagonal(rho.matrix, dim1=-2, dim2=-1).sum(-1).real

    # ------------------------------------------------------------------ local jumps
    def local_kraus_coefficients(self, kind: str) -> Dict[Tuple[int, int], torch.Tensor]:
        """
        {(dj, branch): c} with c of shape (C, k_max): the single non-zero diagonal of the Kraus operators of
        Phi(X) = sum_i L_i X L_i^*: K[c + dj][a + dj - s, a] = c[c, a] (s the shift of m, a = j - m), where
        branch 0 / 1 is j1 = j - 1/2 / j + 1/2. Entries of invalid transitions are 0; with j_min, the transitions
        out of the lowest kept sector (c + dj = -1) are included (they enter the loss, not the gain).
        """
        if kind not in LOCAL_KINDS:
            raise ValueError(f"kind must be one of {LOCAL_KINDS}")
        if kind in self._kraus_cache:
            return self._kraus_cache[kind]
        N, s = self.N, _SHIFT[kind]
        out = {}
        for dj in (-1, 0, 1):
            for br in (0, 1):
                coef = torch.zeros(self.C, self.k_max, dtype=torch.float64)
                for c, tj in enumerate(self.two_j):
                    tj1 = tj - 1 if br == 0 else tj + 1
                    tjo = tj + 2 * dj
                    if tj1 < 0 or tj1 > N - 1 or tjo < 0 or abs(tjo - tj1) != 1:
                        continue
                    w = math.sqrt(N * multiplicity(N - 1, Fraction(tj1, 2)) / multiplicity(N, Fraction(tj, 2)))
                    for a in range(tj + 1):
                        tm = tj - 2 * a
                        tmo = tm + 2 * s
                        if abs(tmo) > tjo:
                            continue
                        coef[c, a] = w * _site_element(kind, tj1, tjo, tj, tm)
                if coef.abs().max() > 0:
                    out[(dj, br)] = coef
        self._kraus_cache[kind] = out
        return out

    def local_loss_diagonal(self, kind: str) -> torch.Tensor:
        """(C, k_max) diagonal of Phi^*(1) = sum K^* K = sum_i L_i^* L_i, from the Kraus coefficients
        (equal to N/2 + m, N/2 - m and N for sigma_-, sigma_+, sigma_z, also with j_min)."""
        out = torch.zeros(self.C, self.k_max, dtype=torch.float64)
        for coef in self.local_kraus_coefficients(kind).values():
            out = out + coef ** 2
        return out

    def local_jump(self, kind: str, rate: Number = 1.0) -> InterSectorChannel:
        """
        Phi(sigma) = rate sum_i L_i sigma L_i^* on blunt densities, L = sigma_- ('-'), sigma_+ ('+') or sigma_z
        ('z'), as an InterSectorChannel of the algebra into itself with Kraus rank 2; `rate` may be a tensor
        (batch,). As a jump of dynamics.lindblad_evolve it gives the local dissipator rate sum_i D[L_i].
        Memory: 2 C^2 k_max^2 complex numbers per batch element.
        """
        r = torch.as_tensor(rate, dtype=torch.float64).reshape(-1)
        B, C, k = r.shape[0], self.C, self.k_max
        cost.check_memory(2 * B * cost.tensor_bytes((2, C, C, k, k), self.dtype), self.device,
                          f"local jump of N = {self.N} spins as a dense InterSectorChannel")
        K = torch.zeros(B, 2, C, C, k, k, dtype=self.dtype, device=self.device)
        s = _SHIFT[kind]
        for (dj, br), coef in self.local_kraus_coefficients(kind).items():
            for c in range(C):
                d = c + dj
                if not 0 <= d < C:
                    continue
                for a in range(self.two_j[c] + 1):
                    if coef[c, a] != 0:
                        K[:, br, d, c, a + dj - s, a] = coef[c, a]
        K = K * r.sqrt().to(device=self.device, dtype=self.dtype)[:, None, None, None, None, None]
        return InterSectorChannel(self.algebra, self.algebra, K)


def _rate(x: Number) -> torch.Tensor:
    return torch.as_tensor(x, dtype=torch.float64).reshape(-1)


def _shifted(dc: int, da: int, C: int, k: int):
    """Index pairs (source, destination) for X[:, c, a, b] -> out[:, c + dc, a + da, b + da] (clipped)."""
    cs, cd = (slice(0, C - dc), slice(dc, C)) if dc >= 0 else (slice(-dc, C), slice(0, C + dc))
    a_s, a_d = (slice(0, k - da), slice(da, k)) if da >= 0 else (slice(-da, k), slice(0, k + da))
    return (slice(None), cs, a_s, a_s), (slice(None), cd, a_d, a_d)


class PermutationalLindbladian:
    """
    The generator

        d sigma/dt = -i[omega_x J_x + omega_y J_y + omega_z J_z, sigma]
                     + collective_emission D[J_-] + collective_pumping D[J_+] + collective_dephasing D[J_z]
                     + emission sum_i D[sigma_-^i] + pumping sum_i D[sigma_+^i] + dephasing sum_i D[sigma_z^i]

    on blunt densities of `system.algebra`, with D[L] X = L X L^* - 1/2 {L^* L, X}. Every coefficient is a number or
    a tensor (batch,); they are broadcast to a common batch size B, and one call integrates B parameter sets.
    In the notation of PIQS (qutip.piqs.Dicke): emission = gamma_E, pumping = gamma_P, dephasing = gamma_D / 4
    (PIQS uses J_z^i = sigma_z^i / 2), and the collective rates coincide.

    With omega_x = omega_z = 0 every term maps real matrices to real matrices (J_y is i times a real matrix, the
    Kraus operators are real), so a real algebra (PermutationalSystem(..., complex_valued=False)) can be used; the
    drive about y is the drive about x rotated by exp(-i pi/2 J_z), which commutes with all dissipators.

    Implementation: every operator above is banded in the basis |j, m> (J_+-, sigma_+- shift m by one, the local
    jumps change j by at most one), so the generator is a sum of nine shifted elementwise products
    W_{dc,da} * sigma -> out[c + dc, a + da, b + da] and the tridiagonal Hamiltonian, O(B C k_max^2) operations per
    evaluation. `rhs` assumes a self-adjoint argument (states and every Runge-Kutta stage are).
    """

    PARAMETERS = ('omega_x', 'omega_y', 'omega_z', 'collective_emission', 'collective_pumping',
                  'collective_dephasing', 'emission', 'pumping', 'dephasing')

    def __init__(self, system: PermutationalSystem, omega_x: Number = 0.0, omega_y: Number = 0.0,
                 omega_z: Number = 0.0, collective_emission: Number = 0.0, collective_pumping: Number = 0.0,
                 collective_dephasing: Number = 0.0, emission: Number = 0.0, pumping: Number = 0.0,
                 dephasing: Number = 0.0):
        self.system = system
        given = (omega_x, omega_y, omega_z, collective_emission, collective_pumping, collective_dephasing,
                 emission, pumping, dephasing)
        vals = torch.broadcast_tensors(*[_rate(v) for v in given])
        self.params = dict(zip(self.PARAMETERS, vals))
        self.batch_size = vals[0].shape[0]
        for name in self.PARAMETERS[3:]:
            if bool((self.params[name] < 0).any()):
                raise ValueError(f"{name} must be non-negative")
        self.real = not bool((self.params['omega_x'] != 0).any() or (self.params['omega_z'] != 0).any())
        if not system.algebra.hilbert.complex_valued and not self.real:
            raise ValueError("omega_x J_x and omega_z J_z need a complex algebra (drive about y in a real one)")
        self._prepared = {}

    def _prepare(self, dtype, device):
        key = (dtype, str(device))
        if key in self._prepared:
            return self._prepared[key]
        ps, P = self.system, self.params
        rdt = torch.float32 if dtype in (torch.float32, torch.complex64) else torch.float64
        C, k = ps.C, ps.k_max

        def col(x):                                                    # (B,) -> (B, 1, 1, 1)
            return x.to(device=device, dtype=rdt)[:, None, None, None]

        def nz(name):
            return bool((P[name] != 0).any())

        m = ps.m_values().to(device=device, dtype=rdt)                  # (C, k)
        p = ps.raising_coefficients().to(device=device, dtype=rdt)       # (C, k): <a-1|J_+|a>
        p_next = torch.zeros_like(p)
        p_next[:, :-1] = p[:, 1:]                                       # <a|J_+|a+1>
        mask = ps.algebra._const('block_mask', torch.bool).to(device)
        # diagonal of A = sum_j g_j L_j^* L_j (B, C, k)
        A = (col(P['collective_emission']) * (p_next ** 2)[None, :, :, None]
             + col(P['collective_pumping']) * (p ** 2)[None, :, :, None])[..., 0]
        terms = {}                                                      # (dc, da) -> W (B, C, k, k), real

        def add(key_, W):
            terms[key_] = W if key_ not in terms else terms[key_] + W

        if nz('collective_emission'):                                   # J_- s J_+: [a, b] <- p[a] p[b] s[a-1, b-1]
            add((0, 1), col(P['collective_emission']) * (p_next[:, :, None] * p_next[:, None, :])[None])
        if nz('collective_pumping'):                                    # J_+ s J_-: [a, b] <- p[a+1] p[b+1] s[a+1, b+1]
            add((0, -1), col(P['collective_pumping']) * (p[:, :, None] * p[:, None, :])[None])
        for kind, name in (('-', 'emission'), ('+', 'pumping'), ('z', 'dephasing')):
            if not nz(name):
                continue
            g = col(P[name])
            A = A + g[..., 0] * ps.local_loss_diagonal(kind).to(device=device, dtype=rdt)[None]
            for (dj, br), coef in ps.local_kraus_coefficients(kind).items():
                cf = coef.to(device=device, dtype=rdt)
                add((dj, dj - _SHIFT[kind]), g * (cf[:, :, None] * cf[:, None, :])[None])
        dm = m[:, :, None] - m[:, None, :]
        E = -0.5 * (A[..., :, None] + A[..., None, :])                   # (B, C, k, k)
        if nz('collective_dephasing'):
            E = E - 0.5 * col(P['collective_dephasing']) * (dm ** 2)[None]
        if (0, 0) in terms:                                             # merge the unshifted gain into E
            E = E + terms.pop((0, 0))
        E = E * mask
        if nz('omega_z'):
            E = E.to(dtype) - 1j * col(P['omega_z']).to(dtype) * dm[None].to(dtype)
        E = E.to(dtype)
        # Z = -i H_xy s with H_xy = omega_x J_x + omega_y J_y tridiagonal:
        # (H s)[a, b] = hu[a] s[a+1, b] + hl[a] s[a-1, b], hu = (w_x - i w_y)/2 p[a+1], hl = (w_x + i w_y)/2 p[a]
        ham = nz('omega_x') or nz('omega_y')
        zu = zl = None
        if ham:
            wx = P['omega_x'].to(device=device, dtype=torch.float64)
            wy = P['omega_y'].to(device=device, dtype=torch.float64)
            cu = (-0.5j * (wx - 1j * wy))[:, None, None, None]
            cl = (-0.5j * (wx + 1j * wy))[:, None, None, None]
            pn, pp = p_next.double()[None, :, :, None], p.double()[None, :, :, None]
            zu, zl = cu * pn, cl * pp
            if dtype in (torch.float32, torch.float64):
                zu, zl = zu.real, zl.real
            zu, zl = zu.to(dtype)[:, :, :-1], zl.to(dtype)[:, :, 1:]
        prep = dict(E=E, terms=[(_shifted(dc, da, C, k), W.to(dtype)) for (dc, da), W in sorted(terms.items())],
                    zu=zu, zl=zl)
        self._prepared[key] = prep
        return prep

    def rhs(self, t: float, rho: torch.Tensor) -> torch.Tensor:
        """d sigma/dt for self-adjoint blunt densities of shape (B or 1, C, k_max, k_max); time independent."""
        P = self._prepare(rho.dtype, rho.device)
        B = max(self.batch_size, rho.shape[0])
        out = P['E'] * rho
        if out.shape[0] != B:
            out = out.expand(B, *out.shape[1:]).clone()
        for (src, dst), W in P['terms']:
            Wd = W[src]
            out[dst].addcmul_(Wd, rho[src])
        if P['zu'] is not None:
            Z = torch.zeros_like(out)
            Z[:, :, :-1, :].addcmul_(P['zu'], rho[:, :, 1:, :])
            Z[:, :, 1:, :].addcmul_(P['zl'], rho[:, :, :-1, :])
            out += Z
            out += Z.transpose(-2, -1).conj() if Z.is_complex() else Z.transpose(-2, -1)
        return out

    def evolve(self, rho0: Operator, times, substeps: int = 10, progress: bool = False) -> List[DensityMatrix]:
        """Classical RK4 (dynamics.rk4) on the grid `times`, `substeps` steps per interval; rho(t) for every t.
        RK4 is stable for h |lambda| < 2.8 with lambda the eigenvalues of the generator; spectral_radius_bound gives
        an upper bound. From pure initial states the far tails of the blocks underflow to denormal numbers, which
        are slow on CPUs: torch.set_flush_denormal(True) avoids this."""
        from .dynamics import rk4
        y0 = rho0.matrix.to(self.system.dtype)
        y0 = y0.expand(max(self.batch_size, y0.shape[0]), *y0.shape[1:]).clone()
        ys = rk4(self.rhs, y0, times, substeps, progress)
        ys = 0.5 * (ys + ys.conj().transpose(-2, -1))
        return [DensityMatrix(self.system.algebra, matrix=y, validate=False) for y in ys]

    def spectral_radius_bound(self) -> torch.Tensor:
        """(B,) an upper bound for the moduli of the eigenvalues of the generator (Gershgorin-type, crude)."""
        P, N = self.params, self.system.N
        w = (P['omega_x'] ** 2 + P['omega_y'] ** 2).sqrt() * N + P['omega_z'].abs() * N
        coll = (P['collective_emission'] + P['collective_pumping']) * (N / 2) * (N / 2 + 1) * 2 \
            + P['collective_dephasing'] * N * N / 2
        loc = 2 * N * (P['emission'] + P['pumping']) + 4 * N * P['dephasing']
        return w + coll + loc

    def dynamics_arguments(self):
        """(H, jumps, rates) of the same generator for dynamics.lindblad_evolve / lindblad_rhs: H an Operator,
        collective jumps as Operators, local jumps as InterSectorChannel; the rates are folded into the jumps (one
        set per batch element) and returned as ones."""
        ps, P = self.system, self.params
        alg = ps.algebra
        H = None
        if not self.real or bool((P['omega_y'] != 0).any()):
            if not alg.hilbert.complex_valued:
                raise ValueError("dynamics_arguments with a Hamiltonian needs a complex algebra")
            Hm = sum(P[n].to(ps.dtype)[:, None, None, None] * op.matrix
                     for n, op in (('omega_x', ps.jx()), ('omega_y', ps.jy()), ('omega_z', ps.jz())))
            H = alg.operator(Hm.to(device=ps.device), is_self_adjoint=True)
        jumps = []
        for name, op in (('collective_emission', ps.jminus()), ('collective_pumping', ps.jplus()),
                         ('collective_dephasing', ps.jz())):
            g = P[name]
            if bool((g != 0).any()):
                jumps.append(alg.operator(g.sqrt()[:, None, None, None].to(device=ps.device, dtype=ps.dtype)
                                          * op.matrix))
        for kind, name in (('-', 'emission'), ('+', 'pumping'), ('z', 'dephasing')):
            g = P[name]
            if bool((g != 0).any()):
                jumps.append(ps.local_jump(kind, g))
        return H, jumps, [1.0] * len(jumps)
