**English** · [Español](es/08_dynamics.md) · [Français](fr/08_dynamics.md) · [Deutsch](de/08_dynamics.md) · [中文](zh/08_dynamics.md) · [日本語](ja/08_dynamics.md) · [Русский](ru/08_dynamics.md)

# Tutorial 8. Dynamics

This tutorial accompanies [Chapter 8](../lecture_notes/08_dynamics.md) of the lecture notes, which proves that a Hamiltonian
in $`\mathcal M`$ conserves the sector weights ([Theorem 8.3](../lecture_notes/08_dynamics.md#r-8-3)), that generators of GKSL form produce
quantum dynamical semigroups ([Theorem 8.6](../lecture_notes/08_dynamics.md#r-8-6)), and treats jumps between sectors, post-selection and
stationary states. We use the exact and the Runge–Kutta solvers of `dynamics` for the Schrödinger, von Neumann
and Lindblad equations, and check the results against SciPy and against GKSL generators on $`M_D`$.

<a id="s-8-1"></a>

## 8.1 Closed systems

*Lecture notes:* [Section 8.1](../lecture_notes/08_dynamics.md#s-8-1); [Proposition 8.2](../lecture_notes/08_dynamics.md#r-8-2) (the propagator $`U(t)=e^{-itH}`$), [Theorem 8.3](../lecture_notes/08_dynamics.md#r-8-3) (conservation of charge).

By Proposition 8.2(1) the evolution splits into $`C`$ independent problems of sizes $`k_c`$, one per block.
The library computes $`U(t)`$ with one eigendecomposition per block size (`dynamics.propagator`), and
`dynamics.von_neumann` evolves a density to many times at once in the eigenbasis, where
$`\rho(t)_{ij}=e^{-it(\lambda_i-\lambda_j)}\rho_{ij}`$. We check against `scipy.linalg.expm`:

```python
import math
import numpy as np
import scipy.integrate
import scipy.linalg as sla
import torch
from torch_vn_algebra import TypeIAlgebra, DensityMatrix, InterSectorChannel, channels, dynamics, states

torch.manual_seed(0)
f64, c128 = torch.float64, torch.complex128
alg = TypeIAlgebra([2, 3, 4], [2, 3, 4], precision='double')     # M_2 (+) M_3 (+) M_4, D = 9
H = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
t = 0.7
U = dynamics.propagator(H, t)                                      # exp(-itH), block by block
for c, k in enumerate(alg.k_factors):                              # independent reference
    ref = sla.expm(-1j * t * H.matrix[0, c, :k, :k].numpy())
    assert np.abs(U.matrix[0, c, :k, :k].numpy() - ref).max() < 1e-12
assert (U.adjoint() @ U - alg.identity()).matrix.abs().max() < 1e-12

rho0 = states.random_density_matrix(alg)
times = [0.0, 0.5, 1.0, 2.0]
traj = dynamics.von_neumann(rho0, H, times)                        # exact, all times at once
for s, rho_s in zip(times, traj):
    U_s = dynamics.propagator(H, s)
    assert (rho_s.matrix - (U_s @ rho0 @ U_s.adjoint()).matrix).abs().max() < 1e-12
    assert (rho_s.sector_probabilities() - rho0.sector_probabilities()).abs().max() < 1e-12   # Theorem 8.3
A = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
A_t = U.adjoint() @ A @ U                                          # Heisenberg picture
assert abs(rho0.expectation(A_t) - dynamics.von_neumann(rho0, H, [t])[0].expectation(A)) < 1e-12
```

A unit vector $`\psi=\bigoplus_c\psi_c`$ defines on $`\mathcal M`$ the state $`\bigoplus_c|\psi_c\rangle\langle\psi_c|`$, which
evolves by $`U(t)`$ ([lecture notes, Section 8.1](../lecture_notes/08_dynamics.md#s-8-1)), so `dynamics.schrodinger` (vectors in the block layout
`(batch, C, k_max)`) and `von_neumann` must agree on it. For a time-dependent $`H(t)`$
the library uses fixed-step classical Runge–Kutta (`dynamics.rk4`, `schrodinger_rk4`): a step multiplies
$`\psi`$ by a polynomial in values of $`H(t)\in\mathcal M`$, so it never moves amplitude between sectors,
but it is not unitary: for constant $`H`$ a step of size $`h`$ multiplies an eigenvector of $`H`$ with eigenvalue
$`\lambda`$ by $`T(-ih\lambda)`$, $`T(z)=\sum_{k=0}^4z^k/k!`$, and $`|T(ix)|^2=1-\frac{x^6}{72}+\frac{x^8}{576}<1`$ for
$`0<|x|<2\sqrt2`$ (expand). So RK4 does not conserve the sector weights exactly: for constant $`H`$, a sector
component with a nonzero coefficient on an eigenvector with $`0<|h\lambda|<2\sqrt2`$ loses weight in every step.
How much they change in a run is a numerical observation (below).
(For one large sector, `krylov.evolve` needs only matrix–vector products; see [Tutorial 9, Section 9.5](../tutorial/09_spin_chains.md#s-9-5).)

```python
rows = torch.tensor([[1, 1, 0, 0], [1, 1, 1, 0], [1, 1, 1, 1]], dtype=c128)   # entries of C^2 (+) C^3 (+) C^4
psi0 = torch.randn(1, alg.C, alg.k_max, dtype=c128) * rows
psi0 = psi0 / psi0.norm()
psi = dynamics.schrodinger(psi0, H, times)                         # shape (T, batch, C, k_max)
p_vec = (psi.abs() ** 2).sum(-1)                                   # ||P_c psi(t)||^2
assert (p_vec - p_vec[0]).abs().max() < 1e-12
mix = DensityMatrix(alg, matrix=psi0[..., :, None] * psi0[..., None, :].conj())   # the state of psi0 on M
outer = psi[..., :, None] * psi[..., None, :].conj()
assert all((o - r.matrix).abs().max() < 1e-12 for o, r in zip(outer, dynamics.von_neumann(mix, H, times)))

V = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
H_of_t = lambda s: H.matrix + math.cos(3 * s) * V.matrix          # driven, in M at every time
psi_rk = dynamics.schrodinger_rk4(psi0, H_of_t, times, substeps=50, algebra=alg)
drift = ((psi_rk.abs() ** 2).sum(-1) - p_vec[0]).abs().max().item()
print(f"RK4, time-dependent H: largest change of a sector weight {drift:.1e}")
assert drift < 1e-8 and (dynamics.schrodinger_rk4(psi0, H, times, substeps=50) - psi).abs().max() < 1e-7

```

Numerical observation: in this run the sector weights change by $`1.1\cdot10^{-9}`$ under RK4 (the printed
value) and by less than the tolerance $`10^{-12}`$ of the asserts under the exact solvers, for which
Theorem 8.3 gives exact conservation.

<a id="s-8-2"></a>

## 8.2 Quantum dynamical semigroups

*Lecture notes:* [Section 8.2](../lecture_notes/08_dynamics.md#s-8-2); [Theorem 8.6](../lecture_notes/08_dynamics.md#r-8-6) (GKSL generators and the Dyson expansion), [Theorem 8.7](../lecture_notes/08_dynamics.md#r-8-7).

The library has four Lindblad routines: `lindblad_rhs` returns the right-hand side as a function for
`rk4`; `lindblad_evolve` integrates with RK4; `lindblad_superoperator` returns the
$`k_c^2\times k_c^2`$ matrix of $`\mathcal L`$ on row-major $`\mathrm{vec}\,\rho_c`$ for every sector; and
`lindblad_channel` computes $`e^{t\mathcal L}`$ exactly, sector by sector, as a `Channel`. We check
Theorem 8.6 for arbitrary jumps on $`M_3`$: the Choi matrix of $`e^{t\mathcal L}`$ is positive (Theorem
6.2), the map is TP, the semigroup law holds, RK4 agrees with the exponential, and the jump part
$`e^{t\mathcal L}-S_0(t)=\sum_{n\ge1}S_n(t)`$ is itself CP:

```python
M3 = TypeIAlgebra([3], [3], precision='double')
H3 = M3.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
jumps = [M3.operator(torch.randn(1, 1, 3, 3, dtype=c128) / 2) for _ in range(2)]     # arbitrary L_j
T = lambda s: dynamics.lindblad_channel(M3, H3, jumps, t=s)        # exp(sL), exact
S = lambda s: T(s).superoperator()[0][0]                           # its 9 x 9 matrix on row-major vec
choi = lambda X: X.reshape(3, 3, 3, 3).permute(2, 0, 3, 1).reshape(9, 9)   # J_{(i,a),(j,b)} = S_{(a,b),(i,j)}
assert (choi(S(0.8)) - T(0.8).choi()[0][0]).abs().max() < 1e-12
assert torch.linalg.eigvalsh(choi(S(0.8))).min() > -1e-12         # completely positive
assert T(0.8).trace_preservation_error().max() < 1e-12             # trace preserving
assert (S(0.3) @ S(0.5) - S(0.8)).abs().max() < 1e-12              # semigroup

rho3 = states.random_density_matrix(M3)
rk = dynamics.lindblad_evolve(rho3, H3, jumps, [0.0, 0.8], substeps=40)[-1]
assert (rk.matrix - T(0.8)(rho3).matrix).abs().max() < 1e-7        # RK4 against exp(tL)
assert abs(dynamics.lindblad_rhs(H3, jumps)(0.0, rho3.matrix)[0, 0].trace()) < 1e-12   # Tr L(rho) = 0

H_eff = H3.matrix - 0.5j * sum(L.adjoint().matrix @ L.matrix for L in jumps)
S0 = channels.Channel(M3, torch.linalg.matrix_exp(-0.8j * H_eff)[:, None]).superoperator()[0][0]
jump_part = torch.linalg.eigvalsh(choi(S(0.8) - S0))
print(f"smallest eigenvalue of the Choi matrix of the jump part: {jump_part.min().item():.1e}")
assert jump_part.min() > -1e-12                                    # sum_{n >= 1} S_n(t) is CP
```

<a id="s-8-3"></a>

## 8.3 Lindblad equations in an algebra with sectors

*Lecture notes:* [Section 8.3](../lecture_notes/08_dynamics.md#s-8-3): the generator $`\mathcal L`$ with jumps inside and between the sectors, which `dynamics.lindblad_rhs` implements; [Proposition 8.8](../lecture_notes/08_dynamics.md#r-8-8), [Proposition 8.9](../lecture_notes/08_dynamics.md#r-8-9) (sector weights).

**Example: particle loss** ([lecture notes, Section 8.3](../lecture_notes/08_dynamics.md#s-8-3)). Two bosonic modes with at most two
particles, $`\mathcal M=M_1\oplus M_2\oplus M_3`$ (sectors $`N=0,1,2`$), the Bose–Hubbard dimer
$`H=J(a_1^*a_2+a_2^*a_1)+\frac U2\sum_jn_j(n_j-1)`$, and each mode loses particles at rate $`\gamma`$: one
`InterSectorChannel` with the Kraus operators $`a_1,a_2`$, restricted to $`N\to N-1`$. By
[Proposition 8.9(2)](../lecture_notes/08_dynamics.md#r-8-9) the sector weights obey $`p'=Qp`$ with $`Q=R-\mathrm{diag}(\sum_dR_{dc})_c`$ for rate $`1`$,
$`R`$ = `transition_matrix()`, and from $`N=2`$ they follow the binomial law $`p_2=q^2`$, $`p_1=2q(1-q)`$, $`p_0=(1-q)^2`$,
$`q=e^{-\gamma t}`$, whatever the hopping and the interaction do. The code checks
$`\Phi^*(1)=0\oplus1\oplus2`$, i.e. $`a_1^*a_1+a_2^*a_2=N`$ on the blocks used, and $`p(t)`$ against both formulas.

```python
fock = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double', charges=[0, 1, 2])   # N = 0, 1, 2
r2 = math.sqrt(2)
a1 = {(0, 1): [[1.0, 0.0]], (1, 2): [[r2, 0.0, 0.0], [0.0, 1.0, 0.0]]}           # a_1 : N -> N - 1
a2 = {(0, 1): [[0.0, 1.0]], (1, 2): [[0.0, 1.0, 0.0], [0.0, 0.0, r2]]}           # a_2 : N -> N - 1
blocks = lambda a: {dc: torch.tensor(K, dtype=f64) for dc, K in a.items()}
loss = InterSectorChannel.from_blocks(fock, fock, {dc: [blocks(a1)[dc], blocks(a2)[dc]] for dc in a1})
J, U_int = 1.0, 0.8                                                # hopping and interaction
H_bh = fock.from_blocks([torch.zeros(1, 1, dtype=f64), torch.tensor([[0.0, J], [J, 0.0]], dtype=f64),
                         torch.tensor([[U_int, r2 * J, 0.0], [r2 * J, 0.0, r2 * J], [0.0, r2 * J, U_int]], dtype=f64)])

assert (loss.adjoint()(fock.identity()).matrix - fock.central([0.0, 1.0, 2.0]).matrix).abs().max() < 1e-12
R = loss.transition_matrix()[0]                                    # r_dc (columns: from sector c)
Q = R - torch.diag(R.sum(0))                                       # Proposition 8.9(2), for rate 1
print(Q)

gamma = 0.5
psi2 = torch.randn(3, dtype=c128)
psi2 = psi2 / psi2.norm()
m0 = torch.zeros(1, 3, 3, 3, dtype=c128)
m0[0, 2] = torch.outer(psi2, psi2.conj())                          # a pure state with N = 2
ts = [0.0, 0.5, 1.0, 2.0, 4.0]
traj = dynamics.lindblad_evolve(DensityMatrix(fock, matrix=m0), H_bh, [loss], ts, rates=[gamma], substeps=100)
for s, rho_s in zip(ts, traj):
    q = math.exp(-gamma * s)
    p = rho_s.sector_probabilities()[0]
    assert (p - torch.tensor([(1 - q) ** 2, 2 * q * (1 - q), q ** 2], dtype=f64)).abs().max() < 1e-8
    assert (p - torch.from_numpy(sla.expm(gamma * s * Q.numpy())[:, 2])).abs().max() < 1e-8
neg = min(w.min().item() for rho_s in traj for w in rho_s.eigenvalues())
print(f"smallest eigenvalue of rho(t) under RK4: {neg:.1e}")
assert neg > -1e-7                                                 # positive up to the RK4 error

sm = torch.tensor([[0.0, 1.0], [0.0, 0.0]], dtype=c128)           # |g><e|; H = 0, rate 3, one step h = 1
decay = lambda s, X: 3.0 * (sm @ X @ sm.mH - 0.5 * (sm.mH @ sm @ X + X @ sm.mH @ sm))
step = dynamics.rk4(decay, torch.diag(torch.tensor([0.0, 1.0], dtype=c128)), [0.0, 1.0], substeps=1)[-1]
assert (step - torch.diag(torch.tensor([-3 / 8, 11 / 8], dtype=c128))).abs().max() < 1e-14   # not positive
```

RK4 is not guaranteed to preserve positivity: for $`\rho'=\mathcal L\rho`$ a step of size $`h`$ is
$`\rho\mapsto T(h\mathcal L)\rho`$ with the polynomial $`T`$ of Section 8.1 (insert the four stages), which is
trace preserving (as $`\mathrm{Tr}\,\mathcal L(X)=0`$ for all $`X`$) but need not be positive. For example, on $`M_2`$ with
basis $`|g\rangle,|e\rangle`$, $`H=0`$ and one jump $`|g\rangle\langle e|`$ with rate $`\gamma`$, $`\mathcal L`$ maps
$`\mathrm{diag}(a,b)`$ to $`\gamma b\mathrm{diag}(1,-1)`$, so one step with $`h\gamma=3`$ maps $`\mathrm{diag}(0,1)`$ to
$`\mathrm{diag}(1-T(-3),T(-3))=\mathrm{diag}(-\frac38,\frac{11}8)`$. Numerical observation: in this run the smallest eigenvalue of
$`\rho(t)`$ is $`\approx-2\cdot10^{-8}`$ (the printed value, with steps up to $`0.02`$). `lindblad_channel` is exact
but treats only jumps inside the sectors. Finally we check Proposition 8.8(2) for a random map between the
sectors and Proposition 8.8(3) for the particle loss, against GKSL generators on $`M_6`$ built with numpy:

```python
Phi = channels.random_inter_sector_channel(fock, fock, kraus_rank=2)
rho = states.random_density_matrix(fock)
on_M = dynamics.lindblad_rhs(None, [Phi], algebra=fock)(0.0, rho.matrix)[0].numpy()
ks, K = fock.k_factors, Phi.kraus[0].numpy()                       # K[l, d, c] = K_l^{dc}, zero padded
embed = lambda Y: sla.block_diag(*[Y[c, :k, :k] for c, k in enumerate(ks)])    # (+)_c Y_c in M_6
Ls = [np.block([[K[l, d, c, :kd, :kc] for c, kc in enumerate(ks)] for d, kd in enumerate(ks)])
      for l in range(2)]                                           # L_l = sum_{d,c} J_l^{dc}
X = embed(rho.matrix[0].numpy())
full = sum(L @ X @ L.conj().T - 0.5 * (L.conj().T @ L @ X + X @ L.conj().T @ L) for L in Ls)   # GKSL on M_6
E_full = full * embed(np.ones((3, 3, 3)))                          # E(X) = sum_d P_d X P_d
assert np.abs(E_full - embed(on_M)).max() < 1e-12                  # Proposition 8.8(2)
print(f"coherences between sectors created on M_6: {np.abs(full - E_full).max():.2f}")

zero = lambda d, c: np.zeros((ks[d], ks[c]))
L_loss = [np.block([[blocks(a)[(d, c)].numpy() if (d, c) in a else zero(d, c) for c in range(3)] for d in range(3)])
          for a in (a1, a2)]                                       # a_1, a_2 on M_6: one block per block row and column
full_loss = sum(L @ X @ L.T - 0.5 * (L.T @ L @ X + X @ L.T @ L) for L in L_loss)
loss_on_M = dynamics.lindblad_rhs(None, [loss], algebra=fock)(0.0, rho.matrix)[0].numpy()
assert np.abs(full_loss - embed(loss_on_M)).max() < 1e-12         # Proposition 8.8(3): no compression needed
```

<a id="s-8-4"></a>

## 8.4 Post-selection and the no-jump evolution

*Lecture notes:* [Section 8.4](../lecture_notes/08_dynamics.md#s-8-4); [Proposition 8.10](../lecture_notes/08_dynamics.md#r-8-10) (the no-jump probability $`w_0`$), [Proposition 8.11](../lecture_notes/08_dynamics.md#r-8-11) (post-selection on "no particle lost").

The library's `states.lueders_update(rho, P)` returns $`P\rho P/\mathrm{Tr}(P\rho P)`$ and
$`\mathrm{Tr}(\rho P)`$. Let the two modes now leak at different rates $`\gamma_1>\gamma_2`$; then
$`A=\gamma_1n_1+\gamma_2n_2`$ is not central (in the sector $`N=2`$ it is $`\mathrm{diag}(2\gamma_1,\gamma_1+\gamma_2,2\gamma_2)`$ in the
basis $`|20\rangle,|11\rangle,|02\rangle`$), Proposition 8.9(2) does not apply, and the post-selected state on $`N=2`$ is
the normalised evolution by $`H_{\mathrm{eff}}`$ whose anti-Hermitian part $`-\frac i2A`$ is not a multiple of $`1`$ in that sector. We check Proposition 8.11 for $`n=0`$ and, by computing
the integral for $`S_1`$ with Simpson's rule, for $`n=1`$:

```python
g1, g2 = 0.6, 0.2                                                  # mode 1 leaks faster
loss1 = InterSectorChannel.from_blocks(fock, fock, {dc: [K] for dc, K in blocks(a1).items()})
loss2 = InterSectorChannel.from_blocks(fock, fock, {dc: [K] for dc, K in blocks(a2).items()})
m1 = torch.zeros(1, 3, 3, 3, dtype=c128)
m1[0, 2, 0, 0] = 1.0                                               # |20>: both particles in mode 1
ts = [0.0, 0.5, 1.0, 2.0]
traj = dynamics.lindblad_evolve(DensityMatrix(fock, matrix=m1), H_bh, [loss1, loss2], ts, rates=[g1, g2],
                                substeps=100)

n_1 = [np.diag(x) for x in ([0.0], [1.0, 0.0], [2.0, 1.0, 0.0])]   # n_1 in each sector; n_2 = N - n_1
H_eff = [H_bh.matrix[0, N, :N + 1, :N + 1].numpy() - 0.5j * (g1 * n + g2 * (N * np.eye(N + 1) - n))
         for N, n in enumerate(n_1)]

def S0(s, N, X):                                  # e^{-isH_eff} X e^{isH_eff^*} in sector N, s a time or an array
    E = torch.linalg.matrix_exp(torch.from_numpy(-1j * np.asarray(s)[..., None, None] * H_eff[N])).numpy()
    return E @ X @ E.conj().swapaxes(-1, -2)

rho_20 = m1[0, 2].numpy()
for s, rho_s in zip(ts, traj):
    w0 = np.trace(S0(s, 2, rho_20)).real                           # probability of no loss
    post, prob = states.lueders_update(rho_s, fock.central([0.0, 0.0, 1.0]))   # post-select N = 2
    assert abs(prob.item() - w0) < 1e-8                            # Proposition 8.11, n = 0
    assert np.abs(post.matrix[0, 2].numpy() - S0(s, 2, rho_20) / w0).max() < 1e-7
closed = dynamics.von_neumann(DensityMatrix(fock, matrix=m1), H_bh, [2.0])[0]   # no losses at all
dist = 0.5 * np.abs(np.linalg.eigvalsh(S0(2.0, 2, rho_20) / w0 - closed.matrix[0, 2].numpy())).sum()
print(f"t = 2: P(no loss) = {w0:.3f} (without hopping: {math.exp(-2 * g1 * 2):.3f}); "
      f"trace distance post-selected/closed system: {dist:.3f}")

K1, K2 = blocks(a1)[(1, 2)].numpy(), blocks(a2)[(1, 2)].numpy()   # S_1(2) rho = int_0^2 S_0(2-s) Phi S_0(s) rho ds
jump = lambda X: g1 * K1 @ X @ K1.T + g2 * K2 @ X @ K2.T
s_grid = np.linspace(0.0, 2.0, 201)
S1 = scipy.integrate.simpson(S0(2.0 - s_grid, 1, jump(S0(s_grid, 2, rho_20))), x=s_grid, axis=0)
assert np.abs(S1 - traj[-1].matrix[0, 1, :2, :2].numpy()).max() < 1e-8     # Proposition 8.11, n = 1
```

Numerical observation: at $`t=2`$ the weight $`p_2(2)=w_0(2)`$ of the sector $`N=2`$ ("no loss") is $`0.230`$ (printed
above). Without hopping ($`J=0`$) it would be exactly $`e^{-2\gamma_1t}=e^{-2.4}=0.091`$, since then $`|20\rangle`$ is an eigenvector of
$`H_{\mathrm{eff}}`$ with eigenvalue $`U-i\gamma_1`$ and $`w_0(t)=|e^{-it(U-i\gamma_1)}|^2`$ (Proposition 8.10); so
here the sector weights do depend on the Hamiltonian. The post-selected state is at trace distance
$`0.135`$ (printed above) from the state of the closed system: conditioning on $`N=2`$ changes the state, beyond the unitary evolution.

<a id="s-8-5"></a>

## 8.5 Stationary states and relaxation

*Lecture notes:* [Section 8.5](../lecture_notes/08_dynamics.md#s-8-5); [Corollary 8.13](../lecture_notes/08_dynamics.md#r-8-13), [Proposition 8.14](../lecture_notes/08_dynamics.md#r-8-14) (spontaneous emission).

By [Proposition 8.14](../lecture_notes/08_dynamics.md#r-8-14), spontaneous emission on $`M_2`$ ($`H=\omega|e\rangle\langle e|`$, one jump
$`|g\rangle\langle e|`$ with rate $`\gamma`$) has $`e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}}`$, the spectrum
$`\{0,-\gamma,-\gamma/2\pm i\omega\}`$ and the spectral gap $`\gamma/2`$. We check the spectrum
with `lindblad_superoperator`, the channel and the trace distance with `lindblad_channel`, and
[Corollary 8.13](../lecture_notes/08_dynamics.md#r-8-13) for random jumps inside the sectors of $`M_2\oplus M_3`$:

```python
q2 = TypeIAlgebra([2], [2], precision='double')                    # basis |g>, |e>
omega, gamma = 1.3, 0.4
H_at = q2.from_blocks([torch.tensor([[0.0, 0.0], [0.0, omega]], dtype=f64)])
sigma_minus = q2.from_blocks([torch.tensor([[0.0, 1.0], [0.0, 0.0]], dtype=f64)])   # |g><e|
gen = dynamics.lindblad_superoperator(q2, H_at, [sigma_minus], rates=[gamma])[0][0]
spec = sorted(torch.linalg.eigvals(gen).tolist(), key=lambda z: (z.real, z.imag))
print([f"{z.real:+.2f}{z.imag:+.2f}i" for z in spec])
assert np.allclose(spec, [-gamma, -gamma / 2 - 1j * omega, -gamma / 2 + 1j * omega, 0.0], atol=1e-12)

psi_at = torch.tensor([0.6, 0.8], dtype=c128)                      # rho_ee(0) = 0.64, rho_ge(0) = 0.48
rho_at = DensityMatrix(q2, matrix=torch.outer(psi_at, psi_at.conj())[None, None])
ground = DensityMatrix(q2, matrix=torch.tensor([[[[1.0, 0.0], [0.0, 0.0]]]], dtype=c128))
for s in (0.5, 2.0, 10.0):
    ch = dynamics.lindblad_channel(q2, H_at, [sigma_minus], t=s, rates=[gamma])
    UA = channels.unitary_channel(dynamics.propagator(H_at, s)).compose(
        channels.amplitude_damping_channel(q2, 1 - math.exp(-gamma * s)))
    assert (ch.superoperator()[0] - UA.superoperator()[0]).abs().max() < 1e-12
    r = ch(rho_at).matrix[0, 0]
    assert abs(r[1, 1].real - 0.64 * math.exp(-gamma * s)) < 1e-12
    assert abs(r[0, 1] - 0.48 * np.exp((1j * omega - gamma / 2) * s)) < 1e-12
    exact = math.sqrt(0.64 ** 2 * math.exp(-2 * gamma * s) + 0.48 ** 2 * math.exp(-gamma * s))
    assert abs(states.trace_distance(ch(rho_at), ground).item() - exact) < 1e-12

two = TypeIAlgebra([2, 3], [2, 3], precision='double')
H_two = two.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
L_two = [two.operator(torch.randn(1, 2, 3, 3, dtype=c128)) for _ in range(2)]
kernel = [int((torch.linalg.eigvals(G_c[0]).abs() < 1e-10).sum())
          for G_c in dynamics.lindblad_superoperator(two, H_two, L_two)]
print("stationary states per sector:", kernel)
assert all(dim >= 1 for dim in kernel)                             # Corollary 8.13
```

Numerical observation: for these random jumps the code finds exactly one stationary state per sector.
By Theorem D.6, applied sector by sector ([lecture notes, Section 8.5](../lecture_notes/08_dynamics.md#s-8-5)), this holds if every sector
has a faithful stationary state and the only matrices commuting with $`H_c`$ and all $`L_{j,c}`$, $`L_{j,c}^*`$ are
the multiples of $`1`$; the code checks these two hypotheses (in floating point, for the rates $`g_j=1`$ used
above) and then asserts the conclusion.

```python
for c, k in enumerate(two.k_factors):                              # Theorem D.6 in every sector
    G_c = dynamics.lindblad_superoperator(two, H_two, L_two)[c][0].numpy()   # on row-major vec(rho_c)
    w_c, V_c = np.linalg.eig(G_c)
    sigma_c = V_c[:, np.argmin(abs(w_c))].reshape(k, k)
    sigma_c = sigma_c / np.trace(sigma_c)
    assert np.abs(G_c @ sigma_c.reshape(-1)).max() < 1e-12            # stationary
    assert np.abs(sigma_c - sigma_c.conj().T).max() < 1e-10 and np.linalg.eigvalsh(sigma_c).min() > 1e-3   # faithful
    ops = [H_two.matrix[0, c, :k, :k].numpy()] + [L.matrix[0, c, :k, :k].numpy() for L in L_two]
    ops += [A.conj().T for A in ops[1:]]
    comm = np.concatenate([np.kron(A, np.eye(k)) - np.kron(np.eye(k), A.T) for A in ops])   # X -> [A, X]
    assert int((np.linalg.svd(comm, compute_uv=False) < 1e-9).sum()) == 1   # commutant = C 1
assert kernel == [1, 1]                                            # the conclusion of Theorem D.6
```

<a id="s-8-6"></a>

## 8.6 What to compare with

QuTiP (`sesolve`, `mesolve` for the Lindblad equation with adaptive ODE solvers, `mcsolve` for
quantum trajectories, `steadystate`, `liouvillian`), dynamiqs (JAX; batched and differentiable
Schrödinger and Lindblad solvers on GPUs), QuantumOptics.jl, and SciPy (`scipy.linalg.expm`, and
`scipy.sparse.linalg.expm_multiply` for $`e^{t\mathcal L}\rho`$ with a sparse vectorised generator) solve
the same equations on one Hilbert space, and for the problems above they must give the same numbers.
On $`\mathcal M`$ one passes them $`D\times D`$ block-diagonal matrices; the vectorised generator then
lives on a $`D^2`$-dimensional space (sparse, but with the sectors coupled in one solve) unless the user
splits the problem into sectors by hand. The library builds the block structure in: exact propagators
and Lindblad channels sector by sector ($`k_c^2\times k_c^2`$ superoperators), jumps between sectors as
`InterSectorChannel`s with the compression of Proposition 8.8, batches of Hamiltonians and jump
operators (for disorder averages) and automatic differentiation (the last two also in dynamiqs). It
has no adaptive step-size control (RK4 is fixed step and only approximately positive), no trajectory
solver, and `lindblad_channel` does not accept jumps between sectors; for those, QuTiP or dynamiqs on
$`M_D`$ are the reference. Besides RK4, `lindblad_evolve` has two time steppers whose steps are channels for
every step size (Section 8.7): the second-order splitting `method='kraus'` and `method='exact'`, which also
treats jumps between sectors, through the matrix of $`\mathcal L`$ on $`\mathcal M`$ (for small sectors).

<a id="s-8-7"></a>

## 8.7 Time steps that are channels

*Lecture notes:* [Section 8.7](../lecture_notes/08_dynamics.md#s-8-7); [Proposition 8.15](../lecture_notes/08_dynamics.md#r-8-15) (RK4 is not positive), [Theorem 8.16](../lecture_notes/08_dynamics.md#r-8-16) (a splitting whose steps are channels).

`lindblad_evolve` has three time steppers, chosen by `method`: `'rk4'` (the default, unchanged) is the step
$`R_h=T_4(h\mathcal L)`$; `'kraus'` is the splitting $`\tilde\Psi_h`$ of Theorem 8.16, built from
$`E_h=e^{-ihH_{\mathrm{eff}}/2}`$, the jump part $`\Phi_{\mathrm{tot}}`$ (jumps inside and between the sectors) and the correction
$`M_h^{-1/2}`$; `'exact'` is $`e^{h\mathcal L}`$, computed from the matrix of $`\mathcal L`$ on $`\bigoplus_cM_{k_c}`$ (padded to
$`C k_{\max}^2`$ coordinates; for a constant generator one exponential per interval length, and `substeps` is not used).
We first check Proposition 8.15 on the two-level system with $`\omega=\gamma=1`$ and $`h=1`$: one RK4 step is
$`T_4(h\mathcal L)`$ with the matrix of `lindblad_superoperator`, $`\tau-|a|^2<0`$, and the pure state
$`\cos\frac\theta2|g\rangle+\sin\frac\theta2|e\rangle`$, $`\theta=0.2746`$, is mapped to a matrix with a negative eigenvalue;
the other two methods give states.

```python
T4 = lambda z: 1 + z + z ** 2 / 2 + z ** 3 / 6 + z ** 4 / 24
om, gm, h = 1.0, 1.0, 1.0                                          # omega = gamma: Proposition 8.15(3)
qb = TypeIAlgebra([2], [2], precision='double')                    # basis |g>, |e>
H_qb = qb.from_blocks([torch.tensor([[0.0, 0.0], [0.0, om]], dtype=f64)])
L_qb = [qb.from_blocks([torch.tensor([[0.0, 1.0], [0.0, 0.0]], dtype=f64)])]   # |g><e|
step = lambda rho, method: dynamics.lindblad_evolve(rho, H_qb, L_qb, [0.0, h], rates=[gm], substeps=1,
                                                    method=method)[-1]
gen = dynamics.lindblad_superoperator(qb, H_qb, L_qb, rates=[gm])[0][0].numpy()
R_h = sum(np.linalg.matrix_power(h * gen, n) / math.factorial(n) for n in range(5))   # T_4(hL)
X = states.random_density_matrix(qb)
assert np.abs(step(X, 'rk4').matrix[0, 0].numpy().reshape(-1) - R_h @ X.matrix[0, 0].numpy().reshape(-1)).max() < 1e-14

tau, a = T4(-h * gm), T4(h * (1j * om - gm / 2))
u = gm ** 2 + 4 * om ** 2
c = [-gm * u * (4 * om ** 2 - 3 * gm ** 2) / 384, u ** 2 * (4 * om ** 2 - 5 * gm ** 2) / 4608, gm * u ** 3 / 9216,
     -u ** 4 / 147456]
assert abs(tau - abs(a) ** 2 - h ** 5 * sum(ci * h ** i for i, ci in enumerate(c))) < 1e-15   # Proposition 8.15(3)
print(f"tau - |a|^2 = {tau - abs(a) ** 2:.4f}")

th = 0.2746
psi_th = torch.tensor([math.cos(th / 2), math.sin(th / 2)], dtype=c128)
pure = DensityMatrix(qb, matrix=torch.outer(psi_th, psi_th.conj())[None, None])
low = {m: torch.linalg.eigvalsh(step(pure, m).matrix[0, 0]).min().item() for m in dynamics.LINDBLAD_METHODS}
print({m: f"{v:.2e}" for m, v in low.items()})
assert low['rk4'] < -8e-5 and low['kraus'] > 0 and low['exact'] > 0
x_star = max(r.real for r in np.roots([1, -4, 12, -24]) if abs(r.imag) < 1e-12)
print(f"x_* = {x_star:.4f}")
```

Numerical observation: here $`\tau-|a|^2=-0.0091`$ and the smallest eigenvalue after one RK4 step is $`-8.5\cdot10^{-5}`$ (the
printed values); $`\theta=0.2746`$ is, to four digits, the angle that makes it smallest. Proposition 8.15 is about one qubit;
in an algebra with sectors the same defect appears in the Choi matrices of the sector components $`R_{dc}`$, which are all
positive iff the step is CP (Proposition 6.4). The function `step_map` below reconstructs a step from `lindblad_evolve` on
Hermitian matrices (the output is made self-adjoint, so the step is evaluated on $`E_{ij}+E_{ji}`$ and $`i(E_{ij}-E_{ji})`$
and extended by linearity). For the particle loss of Section 8.3 we check Theorem 8.16(1) and Proposition 8.8(1): the steps
of `'kraus'` and `'exact'` are CP and TP for step sizes from $`0.1`$ to $`1000`$, while the RK4 step is not CP already
for $`h=0.1`$.

```python
def step_map(alg, H, jumps, rates, h, method):
    """{(d, c): matrix of the sector component R_dc of one step, on row-major vec}"""
    ks_, out = alg.k_factors, {}
    for c, kc in enumerate(ks_):
        E = torch.zeros(kc * kc, alg.C, alg.k_max, alg.k_max, dtype=c128)
        for m in range(kc * kc):
            E[m, c, m // kc, m % kc] = 1.0
        Y = [dynamics.lindblad_evolve(alg.operator(Z), H, jumps, [0.0, h], rates=rates, substeps=1,
                                      method=method)[-1].matrix for Z in (E + E.mH, 1j * (E - E.mH))]
        RE = 0.5 * (Y[0] - 1j * Y[1])                                  # R(E_ij) by linearity
        for d, kd in enumerate(ks_):
            out[d, c] = RE[:, d, :kd, :kd].reshape(kc * kc, kd * kd).T.numpy()
    return out

def choi_min(alg, S):                                              # smallest eigenvalue of all J(R_dc)
    ks_ = alg.k_factors
    return min(np.linalg.eigvalsh(S[d, c].reshape(ks_[d], ks_[d], ks_[c], ks_[c]).transpose(2, 0, 3, 1)
                                  .reshape(ks_[c] * ks_[d], -1)).min() for d, c in S)

def tp_error(alg, S):                                              # Tr sum_d R_dc(X) = Tr X
    ks_ = alg.k_factors
    return max(np.abs(sum(np.trace(S[d, c].reshape(kd, kd, -1), axis1=0, axis2=1) for d, kd in enumerate(ks_))
                      - np.eye(kc).reshape(-1)).max() for c, kc in enumerate(ks_))

for h in (0.1, 3.0, 1000.0):
    for m in dynamics.LINDBLAD_METHODS:
        S = step_map(fock, H_bh, [loss], [0.5], h, m)
        print(f"h = {h:6}, {m:5}: smallest Choi eigenvalue {choi_min(fock, S):+.1e}, trace error {tp_error(fock, S):.0e}")
        if m != 'rk4':
            assert choi_min(fock, S) > -1e-12 and tp_error(fock, S) < 1e-12        # a channel
assert choi_min(fock, step_map(fock, H_bh, [loss], [0.5], 0.1, 'rk4')) < -5e-5
```

Numerical observation: for $`h=0.1`$ the smallest eigenvalue of a Choi matrix of the RK4 step is $`-9.7\cdot10^{-5}`$;
for $`h=3`$, where some eigenvalues $`z`$ of $`h\mathcal L`$ have $`|T_4(z)|>1`$ (RK4 is unstable), it is $`-9.6\cdot10^{2}`$
(printed above); the other two methods stay channels up to rounding. For $`h=0.1`$ all eigenvalues satisfy
$`|T_4(z)|\le1`$: the failure of positivity is not an instability. The renormalisation
by $`M_h^{-1/2}`$ is computed from the eigenvalues of $`M_h`$, whose condition number can grow exponentially with $`h`$
when the decay rates in a sector differ; when it exceeds $`\varepsilon^{-1/4}`$ ($`\varepsilon`$ the machine precision) the
library takes the step as two steps of size $`h/2`$, recursively, so that the computed step stays trace preserving up to
rounding.

Finally the accuracy. `'exact'` agrees with $`e^{t\hat{\mathcal L}}`$ of the generator $`\hat{\mathcal L}`$ on $`M_6`$ with
the jumps $`\sqrt\gamma a_1,\sqrt\gamma a_2`$ (Proposition 8.8(1) and (3), computed with SciPy from the matrices `L_loss` of
Section 8.3), and the trace distance of the other two methods to it falls like $`h^2`$ (Theorem 8.16(3)) and like $`h^4`$, respectively:

```python
rho_f = states.random_density_matrix(fock)
ref = dynamics.lindblad_evolve(rho_f, H_bh, [loss], [0.0, 2.0], rates=[0.5], method='exact')[-1]
H6, I6 = embed(H_bh.matrix[0].numpy()), np.eye(6)
G6 = -1j * (np.kron(H6, I6) - np.kron(I6, H6.T))
for L6 in L_loss:
    L6 = math.sqrt(0.5) * L6
    G6 = G6 + np.kron(L6, L6.conj()) - 0.5 * (np.kron(L6.T @ L6, I6) + np.kron(I6, (L6.T @ L6).T))
exact6 = sla.expm(2.0 * G6) @ embed(rho_f.matrix[0].numpy()).reshape(-1)
assert np.abs(exact6 - embed(ref.matrix[0].numpy()).reshape(-1)).max() < 1e-12

for m, subs in (('kraus', (4, 8, 16, 32)), ('rk4', (4, 8, 16, 32))):
    err = [states.trace_distance(dynamics.lindblad_evolve(rho_f, H_bh, [loss], [0.0, 2.0], rates=[0.5], substeps=n,
                                                          method=m)[-1], ref).item() for n in subs]
    ratios = [e0 / e1 for e0, e1 in zip(err[:-1], err[1:])]
    print(m, [f"{e:.1e}" for e in err], "ratios", [f"{r:.1f}" for r in ratios])
    assert all((3.5 < r < 4.5) if m == 'kraus' else (12 < r < 20) for r in ratios)
```

Numerical observation: halving $`h`$ divides the error of `'kraus'` by about $`4`$ and that of RK4 by about $`16`$ (the
printed ratios): the splitting buys positivity for every step size with a lower order. For a time-dependent $`H`$ or
time-dependent jumps, `'kraus'` and `'exact'` freeze the generator at the midpoint of each step, so every step is still a
channel by Theorem 8.16(1); the order of this frozen-midpoint scheme is not proved in the lecture notes.

<a id="x-8"></a>

## Exercises

[Exercises 2 and 4](../lecture_notes/08_dynamics.md#x-8) of the lecture notes: the eigenvalues
$`-\frac i4(\gamma_1+\gamma_2)\pm\sqrt{J^2-\delta^2}`$, $`\delta=(\gamma_1-\gamma_2)/4`$, of $`H_{\mathrm{eff}}`$ in the sector $`N=1`$ of
Section 8.4, with a single Jordan block of size $`2`$ at $`J=|\delta|`$; and the spectrum $`\{0,0,-2\gamma\pm i\omega\}`$ of
pure dephasing on $`M_2`$ ($`L=\sigma_z`$). Check both numerically.

<details>
<summary>Solutions</summary>

```python
lam = np.linalg.eigvals(H_eff[1])
delta = (g1 - g2) / 4
root = np.sqrt(J ** 2 - delta ** 2 + 0j)
assert np.allclose(sorted(lam, key=np.real), sorted(-0.25j * (g1 + g2) + np.array([root, -root]), key=np.real))
H_ep = np.array([[-0.5j * g1, delta], [delta, -0.5j * g2]])         # J = delta: the exceptional point
Y = H_ep + 0.25j * (g1 + g2) * np.eye(2)                           # H_ep minus its only eigenvalue
assert np.abs(Y @ Y).max() < 1e-15 and np.linalg.matrix_rank(Y) == 1   # one Jordan block of size 2
```

```python
deph = dynamics.lindblad_superoperator(q2, H_at, [q2.from_blocks([torch.diag(torch.tensor([1.0, -1.0], dtype=f64))])],
                                       rates=[gamma])[0][0]
ev = torch.linalg.eigvals(deph)
assert int((ev.abs() < 1e-12).sum()) == 2
assert np.allclose(sorted(ev.tolist(), key=lambda z: (z.real, z.imag))[:2],
                   [-2 * gamma - 1j * omega, -2 * gamma + 1j * omega], atol=1e-12)
```

</details>

**Next:** [Tutorial 9. Spin chains](09_spin_chains.md).
