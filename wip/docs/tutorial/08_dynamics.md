**English** · [Español](es/08_dynamics.md) · [Français](fr/08_dynamics.md) · [Deutsch](de/08_dynamics.md) · [中文](zh/08_dynamics.md) · [日本語](ja/08_dynamics.md) · [Русский](ru/08_dynamics.md)

# 8. Dynamics

A closed quantum system evolves by a one-parameter group of unitaries, an open one — a system
coupled to an environment that is not observed — by a semigroup of channels. In an algebra with
sectors both have a feature that is invisible on a single matrix algebra: a Hamiltonian in
$\mathcal M$ can never change the sector, so the charges are conserved, while an open system can
jump between sectors, for instance by losing a particle. This chapter proves the basic facts about
both kinds of evolution, writes the Lindblad equation with jumps between sectors, and shows that
"no particle was lost" is a measurement of the charge whose post-selected state follows a
non-Hermitian evolution.

**You will learn** why sector weights are constants of motion, why the generators of
Gorini–Kossakowski–Sudarshan–Lindblad (GKSL) form produce channels, how jumps between sectors enter
the Lindblad equation and when the sector weights obey a classical rate equation, what the no-jump
evolution and post-selection mean, and how to find stationary states and relaxation rates.

Throughout, $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ acts on
$\mathcal H=\bigoplus_c\mathbb C^{k_c}$, $D=\sum_ck_c$, $P_c$ are the minimal central projections,
states are given by their blunt densities $\rho$ with sector weights $p_c=\operatorname{Tr}\rho_c$
(Chapter 3), $[A,B]=AB-BA$, $\{A,B\}=AB+BA$, and $\hbar=1$. A **superoperator** is a linear map
on $\mathcal M$ or on $M_D(\mathbb C)$; all norms of superoperators are operator norms, hence
submultiplicative. For a linear map $B$ on a finite-dimensional space, $e^{tB}=\sum_{n\ge0}t^nB^n/n!$.

## 8.1 Closed systems

A **Hamiltonian** is a self-adjoint $H\in\mathcal M$. A vector evolves by the **Schrödinger
equation** $i\psi'(t)=H\psi(t)$, a density by the **von Neumann equation** $\rho'(t)=-i[H,\rho(t)]$.
Both are linear equations with constant coefficients on a finite-dimensional space, and the
following lemma is all the theory of differential equations we need.

**Lemma 8.1.** Let $B$ be a linear map on a finite-dimensional normed space $X$. Then
$\frac d{dt}e^{tB}=Be^{tB}=e^{tB}B$, $e^{(s+t)B}=e^{sB}e^{tB}$, $\|e^{tB}\|\le e^{|t|\,\|B\|}$, and for
every $x_0\in X$ the map $x(t)=e^{tB}x_0$ is the only differentiable $x$ on an interval
$I\ni0$ with $x'=Bx$ and $x(0)=x_0$.

*Proof.* The series converges absolutely, with the norm bound, and as a power series in $t$ it can
be differentiated term by term, which gives the first formula. Then
$\frac d{dt}(e^{tB}e^{-tB})=Be^{tB}e^{-tB}-e^{tB}Be^{-tB}=0$, so $e^{tB}e^{-tB}=1$. If $x'=Bx$, then
$\frac d{dt}(e^{-tB}x(t))=e^{-tB}(x'-Bx)=0$, so $e^{-tB}x(t)=x_0$ and $x(t)=e^{tB}x_0$. Finally
$t\mapsto e^{(s+t)B}x_0$ and $t\mapsto e^{tB}e^{sB}x_0$ solve $y'=By$ with the same initial value. $\square$

**Proposition 8.2.** Let $H=H^*\in\mathcal M$ and $U(t)=e^{-itH}$, defined by the functional
calculus of Section 1.3 with $f_t(\lambda)=e^{-it\lambda}$.

1. $U(t)=\bigoplus_ce^{-itH_c}$ is a unitary in $\mathcal M$, $U(s+t)=U(s)U(t)$ and $U'(t)=-iHU(t)$;
   $U(t)$ equals the exponential series of $-itH$.
2. $\psi(t)=U(t)\psi_0$ is the only solution of $i\psi'=H\psi$ with $\psi(0)=\psi_0$, and
   $\rho(t)=U(t)\rho_0U(t)^*$ is the only solution in $\mathcal M$ of $\rho'=-i[H,\rho]$ with
   $\rho(0)=\rho_0$. If $\rho_0$ is a state, so is every $\rho(t)$.
3. (*Heisenberg picture*) $\operatorname{Tr}(\rho(t)A)=\operatorname{Tr}(\rho_0A(t))$ with
   $A(t)=U(t)^*AU(t)$, and $A'=i[H,A]$.

*Proof.* (1) Diagonalise every block, $H_c=V_c\Lambda_cV_c^*$; then $V=\bigoplus_cV_c$ diagonalises
$H$, and $U(t)=Ve^{-it\Lambda}V^*$ is block diagonal (by Theorem 1.5(a) the result does not depend on
the diagonalising unitary). Since $|f_t|^2=1$ and $f_sf_t=f_{s+t}$, Theorem 1.5(c) gives unitarity and
the group law. Differentiating $e^{-it\Lambda}$ entrywise gives $U'=V(-i\Lambda)e^{-it\Lambda}V^*=-iHU$, so
$U(t)$ is the series $e^{-itH}$ by Lemma 8.1 ($X=M_D$, $B$ the left multiplication by $-iH$).
(2) $(U\psi_0)'=-iHU\psi_0$, and $(U\rho_0U^*)'=-iHU\rho_0U^*+U\rho_0U^*(iH)$ because $(U')^*=iU^*H$.
Uniqueness is Lemma 8.1 with $B=-iH$ on $\mathcal H$, resp. $B=-i[H,\cdot\,]$ on $\mathcal M$. Positivity:
$\langle x,U\rho_0U^*x\rangle=\langle U^*x,\rho_0U^*x\rangle\ge0$; the trace by cyclicity. (3) Cyclicity of
the trace, and $(U^*AU)'=iU^*HAU-iU^*AHU=i[H,A(t)]$. $\square$

**Theorem 8.3 (conservation of charge).** Let $t\mapsto H(t)$ take values in the self-adjoint
elements of $\mathcal M$, on an interval $I$.

1. If $\psi:I\to\mathcal H$ is differentiable with $i\psi'=H(t)\psi$, then every $\|P_c\psi(t)\|$ is constant.
2. If $\rho:I\to\mathcal M$ is differentiable with $\rho'=-i[H(t),\rho]$, then $\operatorname{Tr}(Z\rho(t))$ is
   constant for every central $Z$; in particular every sector weight $p_c(t)=\operatorname{Tr}(P_c\rho(t))$ is.
3. If $H$ is constant and $A\in\mathcal M$ commutes with $H$, then $\operatorname{Tr}(A\rho(t))$ is constant.

*Proof.* (1) The inner product is antilinear in the first argument, so
$\frac d{dt}\langle\psi,P_c\psi\rangle=\langle-iH\psi,P_c\psi\rangle+\langle\psi,-iP_cH\psi\rangle
=i\langle\psi,HP_c\psi\rangle-i\langle\psi,P_cH\psi\rangle=0$, as $P_c$ is central. (2)
$\frac d{dt}\operatorname{Tr}(Z\rho)=-i\operatorname{Tr}(Z[H,\rho])=-i\operatorname{Tr}([Z,H]\rho)=0$ by
cyclicity. (3) The same computation with $A$. $\square$

> **Physical meaning.** A Hamiltonian in $\mathcal M$ commutes with the charges (particle number,
> magnetisation, a symmetry quantum number), so the probability of each charge is a constant of
> motion, even for driven, time-dependent Hamiltonians. Computationally the evolution splits into
> $C$ independent problems: one diagonalisation costs $\sum_ck_c^3$ operations instead of $D^3$.

The library computes $U(t)$ with one eigendecomposition per block size (`dynamics.propagator`), and
`dynamics.von_neumann` evolves a density to many times at once in the eigenbasis, where
$\rho(t)_{ij}=e^{-it(\lambda_i-\lambda_j)}\rho_{ij}$. We check against `scipy.linalg.expm`:

```python
import math
import numpy as np
import scipy.integrate
import scipy.linalg as sla
import torch
from torch_vn_algebra import TypeIAlgebra, DensityMatrix, InterSectorChannel, channels, dynamics, krylov, states

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

A vector $\psi=\bigoplus_c\psi_c$ defines on $\mathcal M$ the state
$A\mapsto\sum_c\langle\psi_c,A_c\psi_c\rangle$, the mixture $\bigoplus_c|\psi_c\rangle\langle\psi_c|$
(Chapter 3, Exercise 3), so `dynamics.schrodinger` (vectors in the block layout
`(batch, C, k_max)`) and `von_neumann` must agree on it. For a time-dependent $H(t)$ the library uses
the classical fixed-step Runge–Kutta method (`dynamics.rk4`, `schrodinger_rk4`). One step multiplies
$\psi$ by a polynomial in values of $H(t)\in\mathcal M$, so it never moves amplitude between sectors;
it is not exactly unitary, and the sector weights drift at the level of the discretisation error.
For one large sector, `krylov.evolve` needs only matrix–vector products (Lanczos; Chapter 9, in
preparation):

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

n, hop = 60, torch.diag(torch.ones(59, dtype=f64), 1)              # one larger sector: a disordered chain
chain = krylov.SparseSectorHamiltonian(torch.randn(1, n, dtype=f64), (hop + hop.T).to_sparse(), c128)
v0 = torch.zeros(1, n, dtype=c128)
v0[0, n // 2] = 1.0                                                # a particle on the middle site
psi_kr = krylov.evolve(chain, v0, [0.0, 2.0], m=20, dt=1.0)        # matrix-vector products only
chain_H = TypeIAlgebra([n], [n], precision='double').operator(chain.to_dense()[:, None])
assert (psi_kr[:, 0] - dynamics.schrodinger(v0[:, None], chain_H, [0.0, 2.0])[:, 0, 0]).abs().max() < 1e-10
```

The sector weights change by $1.1\cdot10^{-9}$ under RK4 and by nothing beyond rounding under the
exact solvers.

## 8.2 Quantum dynamical semigroups

Completely positive (CP) and trace-preserving (TP) maps, and channels, are defined in Section 6.1
([Chapter 6](06_quantum_channels.md)).
Besides Theorem 6.2 and Proposition 6.4 (a map $A\mapsto\sum_lK_lAK_l^*$, also between sectors, is
CP) we need three closure properties.

**Lemma 8.4.** (a) Sums, nonnegative multiples and compositions of CP maps are CP. If
$\Phi$ is CP on $M_D$ and maps a $*$-subalgebra $\mathcal N\subseteq M_D$ (such as $\mathcal M$, or
$P_c\mathcal MP_c\cong M_{k_c}$ with unit $P_c$) into itself, then $\Phi|_{\mathcal N}$ is CP.
(b) If CP maps $\Phi_m$ converge to $\Phi$ (that is, $\Phi_m(A)\to\Phi(A)$ for all $A$), then $\Phi$ is CP.
(c) If $s\mapsto\Psi(s)$ is a continuous map from $[a,b]$ into the CP maps, then $\int_a^b\Psi(s)\,ds$ is CP.

*Proof.* (a) $\Phi\mapsto\Phi_n$ is linear and $(\Phi\Psi)_n=\Phi_n\Psi_n$; sums of positive matrices are
positive. For the restriction: $M_n(\mathcal N)\subseteq M_n(M_D)$, and positivity is the same
condition in both (positivity of a matrix). (b) The positive matrices form a closed set ($\langle x,Xx\rangle\ge0$ for each
$x$ is a closed condition), and $(\Phi_m)_n(X)\to\Phi_n(X)$. (c) The Riemann sums are CP by (a) and
converge to the integral, since $\Psi$ is uniformly continuous; apply (b). $\square$

**Definition.** A **quantum dynamical semigroup** on $\mathcal M$ is a family $(T_t)_{t\ge0}$ of
channels $\mathcal M\to\mathcal M$ with $T_0=\mathrm{id}$, $T_{s+t}=T_sT_t$ and $t\mapsto T_t$ continuous.

**Lemma 8.5.** Every continuous semigroup of superoperators ($T_0=\mathrm{id}$, $T_{s+t}=T_sT_t$) is $T_t=e^{t\mathcal L}$
for exactly one superoperator $\mathcal L$, its **generator**, and $\mathcal L=\lim_{t\to0}(T_t-\mathrm{id})/t$.

*Proof.* $V_h=\frac1h\int_0^hT_s\,ds\to\mathrm{id}$ as $h\to0$, so $V_h$ is invertible for some fixed
small $h$. Now $T_tV_h=\frac1h\int_t^{t+h}T_s\,ds$ is differentiable in $t$ with derivative
$\frac1h(T_{t+h}-T_t)=T_t\frac{T_h-\mathrm{id}}h$. Hence $T_t$ is differentiable and $T_t'=T_t\mathcal L$
with $\mathcal L=\frac{T_h-\mathrm{id}}hV_h^{-1}$, and $T_t=e^{t\mathcal L}$ by Lemma 8.1 (with $B$ the
right composition by $\mathcal L$). Uniqueness: $\mathcal L=T'(0)$. $\square$

**Theorem 8.6 (GKSL generators).** Let $H=H^*$ and $L_1,\dots,L_J$ be in $M_D$, and
$$
\mathcal L(\rho)=-i[H,\rho]+\sum_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
$$
Then $(e^{t\mathcal L})_{t\ge0}$ is a quantum dynamical semigroup on $M_D$. More precisely, put
$H_{\mathrm{eff}}=H-\frac i2\sum_jL_j^*L_j$, $\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)$
and $\Phi(\rho)=\sum_jL_j\rho L_j^*$, so that $\mathcal L=\mathcal L_0+\Phi$, and define
$$
S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},\qquad
S_n(t)=\int_0^tS_0(t-s)\,\Phi\,S_{n-1}(s)\,ds\quad(n\ge1).
$$
Then every $S_n(t)$ is CP, and $e^{t\mathcal L}=\sum_{n\ge0}S_n(t)$ (the **Dyson expansion**), with
$\|S_n(t)\|\le e^{t\|\mathcal L_0\|}(t\|\Phi\|)^n/n!$.

*Proof.* *Step 1.* $X(t)=e^{-itH_{\mathrm{eff}}}\rho e^{itH_{\mathrm{eff}}^*}$ satisfies $X'=\mathcal L_0(X)$
and $X(0)=\rho$, so $S_0(t)=e^{t\mathcal L_0}$ by Lemma 8.1; it is CP with the single Kraus operator
$e^{-itH_{\mathrm{eff}}}$. *Step 2.* Write $a=\|\mathcal L_0\|$, $b=\|\Phi\|$. By induction the $S_n$
are continuous and
$\|S_n(t)\|\le\int_0^te^{(t-s)a}\,b\,e^{sa}\frac{(bs)^{n-1}}{(n-1)!}ds=e^{ta}\frac{(bt)^n}{n!}$, so the
series converges uniformly on every $[0,T]$. *Step 3.* Since
$S_n(t)=S_0(t)\int_0^te^{-s\mathcal L_0}\Phi S_{n-1}(s)\,ds$, the product rule and the fundamental
theorem of calculus give $S_n'=\mathcal L_0S_n+\Phi S_{n-1}$ ($n\ge1$) and $S_0'=\mathcal L_0S_0$. By the
bounds the series of derivatives converges uniformly on $[0,T]$, so $Y(t)=\sum_nS_n(t)$ is
differentiable with $Y'=\mathcal L_0Y+\Phi Y=\mathcal LY$, $Y(0)=\mathrm{id}$, and $Y(t)=e^{t\mathcal L}$ by
Lemma 8.1. *Step 4.* $\Phi$ is CP (Theorem 6.2), so by induction the integrand
$s\mapsto S_0(t-s)\Phi S_{n-1}(s)$ is a continuous family of CP maps and $S_n(t)$ is CP (Lemma 8.4);
so are the partial sums and their limit $e^{t\mathcal L}$. *Step 5.* By cyclicity,
$\operatorname{Tr}\mathcal L(\rho)=\sum_j\big(\operatorname{Tr}(L_j\rho L_j^*)-\operatorname{Tr}(L_j^*L_j\rho)\big)=0$,
hence $\frac d{dt}\operatorname{Tr}e^{t\mathcal L}\rho=\operatorname{Tr}\mathcal L(e^{t\mathcal L}\rho)=0$: every
$e^{t\mathcal L}$ is TP. The semigroup law and continuity are Lemma 8.1. $\square$

The converse is deep, and we state it without proof.

**Theorem 8.7 (Gorini–Kossakowski–Sudarshan; Lindblad).** The generator of every quantum
dynamical semigroup on $M_D$ has the form of Theorem 8.6. [V. Gorini, A. Kossakowski,
E. C. G. Sudarshan, *J. Math. Phys.* 17, 821 (1976); G. Lindblad, *Commun. Math. Phys.* 48, 119
(1976).] For a general von Neumann algebra, Christensen and Evans [*J. London Math. Soc.* 20, 358
(1979)] showed that the generator in the Heisenberg picture is $A\mapsto\Psi(A)+K^*A+AK$ with
$\Psi$ CP and $K\in\mathcal M$.

On $\mathcal M=\bigoplus_cM_{k_c}$ the Schrödinger form of this is $\mathcal L(\rho)=\Phi(\rho)+G\rho+\rho G^*$
with $\Phi$ CP on $\mathcal M$ (with Kraus operators between sectors, Proposition 6.4) and $G\in\mathcal M$.
Trace preservation, $\operatorname{Tr}\big((\Phi^*(1)+G+G^*)\rho\big)=0$ for all $\rho\in\mathcal M$, forces
$G+G^*=-\Phi^*(1)$, i.e. $G=-iH-\frac12\Phi^*(1)$ with $H=H^*\in\mathcal M$: exactly the generators of
the next section.

> **Physical meaning.** A system that interacts weakly with a large environment whose correlations
> decay quickly (the Born–Markov approximation) loses memory of its past: its state at time $s+t$ is
> obtained from the state at time $s$ by a channel that depends only on $t$. Theorem 8.7 says that
> such a dynamics is always a Hamiltonian part plus "jumps" $L_j$ (photon emission, particle loss,
> dephasing) at rates fixed by the normalisation of $L_j$.

The library has four Lindblad routines: `lindblad_rhs` returns the right-hand side as a function for
`rk4`; `lindblad_evolve` integrates with RK4; `lindblad_superoperator` returns the
$k_c^2\times k_c^2$ matrix of $\mathcal L$ on row-major $\operatorname{vec}\rho_c$ for every sector; and
`lindblad_channel` computes $e^{t\mathcal L}$ exactly, sector by sector, as a `Channel`. We check
Theorem 8.6 for arbitrary jumps on $M_3$: the Choi matrix of $e^{t\mathcal L}$ is positive (Theorem
6.2), the map is TP, the semigroup law holds, RK4 agrees with the exponential, and the jump part
$e^{t\mathcal L}-S_0(t)=\sum_{n\ge1}S_n(t)$ is itself CP:

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

## 8.3 Lindblad equations in an algebra with sectors

On $\mathcal M$ there are two kinds of jumps. A jump **inside the sectors** is an operator
$L_j\in\mathcal M$. A jump **between sectors** is a CP map $\Phi:\mathcal M\to\mathcal M$ with Kraus
operators $K^{dc}_l:\mathbb C^{k_c}\to\mathbb C^{k_d}$ (Proposition 6.4, an `InterSectorChannel`):
$\Phi(\rho)_d=\sum_{c,l}K^{dc}_l\rho_cK^{dc*}_l$. By Proposition 6.8 its blunt dual satisfies
$\Phi^*(1)_c=\sum_{d,l}K^{dc*}_lK^{dc}_l$, a block-diagonal element of $\mathcal M$. With rates
$g_j,g_m>0$ the generator of `lindblad_rhs` is
$$
\mathcal L(\rho)=-i[H,\rho]+\sum_jg_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big)
+\sum_mg_m\Big(\Phi_m(\rho)-\tfrac12\{\Phi_m^*(1),\rho\}\Big).
$$

**Proposition 8.8.** Let $J^{dc}_l\in M_D$ be the matrix with the block $K^{dc}_l$ in position $(d,c)$
and zeros elsewhere.

1. $\mathcal L$ maps $\mathcal M$ into $\mathcal M$ and is the restriction to $\mathcal M$ of the GKSL generator
   $\hat{\mathcal L}$ on $M_D$ with jumps $\sqrt{g_j}L_j$ and $\sqrt{g_m}J^{dc}_{m,l}$ (all $m,l,d,c$). Hence
   $(e^{t\mathcal L})_{t\ge0}$ is a quantum dynamical semigroup on $\mathcal M$.
2. Let instead $\tilde{\mathcal L}$ be the GKSL generator on $M_D$ in which the blocks of one Kraus
   index are summed into one jump, $L_{m,l}=\sum_{d,c}J^{dc}_{m,l}$. Then
   $\mathcal L(\rho)=E\big(\tilde{\mathcal L}(\rho)\big)$ for $\rho\in\mathcal M$, where $E(X)=\sum_dP_dXP_d$ is
   the conditional expectation of $M_D$ onto $\mathcal M$; in general $\tilde{\mathcal L}(\rho)\notin\mathcal M$.

*Proof.* Since $\rho=\sum_cP_c\rho P_c$ and $J^{dc}P_{c'}=\delta_{cc'}J^{dc}$,
$J^{dc}\rho J^{dc*}$ is the block $K^{dc}\rho_cK^{dc*}$ in position $(d,d)$, and $J^{dc*}J^{dc}$ is
$K^{dc*}K^{dc}$ in position $(c,c)$. Summing over $d,l$ gives $\Phi(\rho)$ and $\Phi^*(1)$, which proves the
identity in (1), and all terms lie in $\mathcal M$. So $\hat{\mathcal L}$ and all its powers map $\mathcal M$
into $\mathcal M$, $e^{t\hat{\mathcal L}}|_{\mathcal M}=e^{t\mathcal L}$, and this is a channel of $\mathcal M$ by
Theorem 8.6 and Lemma 8.4(a). (2) $L\rho L^*=\sum_{d,d',c}J^{dc}\rho J^{d'c*}$ for $L=\sum_{d,c}J^{dc}$, and
$E$ keeps the terms $d=d'$. For $B\in\mathcal M$, $P_dXBP_d=P_dXP_dB$, so $E(XB)=E(X)B$ and
$E(BX)=BE(X)$; thus $E\{L^*L,\rho\}=\{E(L^*L),\rho\}$, and $E(L^*L)=\sum_cP_cL^*LP_c=\sum_{c,d}J^{dc*}J^{dc}$
because $J^{dc*}J^{d'c}=0$ for $d\ne d'$. Finally $E(-i[H,\rho])=-i[H,\rho]$. $\square$

> **Physical meaning.** Part (2) says that the Lindblad equation on $\mathcal M$ is what remains of an
> equation on $M_D$ when the coherences between sectors that the jumps create are erased at once —
> the content of a superselection rule. For particle loss nothing is erased: the annihilation
> operator $a_j$ lowers the particle number by exactly one, so $a_j\rho a_j^*$ and $a_j^*a_j=n_j$ are
> block diagonal, and the equation on $\mathcal M$ is the usual one.

**Proposition 8.9 (sector weights).** Let $\rho(t)=e^{t\mathcal L}\rho_0$ and $p(t)$ its sector weights.

1. If there are no jumps between sectors, every $p_c(t)$ is constant.
2. If for every $m$ and every pair $d\ne c$, $\sum_lK^{dc*}_{m,l}K^{dc}_{m,l}=r_{m,dc}1_{k_c}$ with numbers
   $r_{m,dc}$, then $p'=Qp$ with $Q_{dc}=\sum_mg_mr_{m,dc}$ for $d\ne c$ and $Q_{cc}=-\sum_{d\ne c}Q_{dc}$.
   Hence $p(t)=e^{tQ}p(0)$ does not depend on $H$, on the $L_j$ or on the states inside the sectors.

*Proof.* $p_d'=\operatorname{Tr}(P_d\mathcal L(\rho))$. The Hamiltonian term vanishes as in Theorem 8.3.
For $L_j\in\mathcal M$, $\operatorname{Tr}(P_dL_j\rho L_j^*)-\operatorname{Tr}(P_dL_j^*L_j\rho)
=\operatorname{Tr}(L_j^*P_dL_j\rho)-\operatorname{Tr}(L_j^*L_jP_d\rho)=0$, as $P_d$ is central. For a jump
$\Phi$ between sectors, $\operatorname{Tr}(P_d\Phi(\rho))=\sum_c\operatorname{Tr}\big(\sum_lK^{dc*}_lK^{dc}_l\rho_c\big)$
and $\operatorname{Tr}(P_d\tfrac12\{\Phi^*(1),\rho\})=\sum_{d'}\operatorname{Tr}\big(\sum_lK^{d'd*}_lK^{d'd}_l\rho_d\big)$;
the terms $c=d$ and $d'=d$ cancel. This proves (1), and under the hypothesis of (2) the remaining terms
are $\sum_{c\ne d}r_{dc}p_c-\sum_{d'\ne d}r_{d'd}p_d$. Lemma 8.1 gives $p(t)=e^{tQ}p(0)$. $\square$

Under the hypothesis of (2), $r_{m,dc}$ is the transition matrix $T_{dc}$ of $\Phi_m$ (Proposition 6.12)
for every state; this is what `transition_matrix()` returns.

**Example: particle loss.** Two bosonic modes with at most two particles; the sectors are
$N=0,1,2$ with bases $|00\rangle$; $|10\rangle,|01\rangle$; $|20\rangle,|11\rangle,|02\rangle$, so
$\mathcal M=M_1\oplus M_2\oplus M_3$. The Bose–Hubbard dimer
$H=J(a_1^*a_2+a_2^*a_1)+\frac U2\sum_jn_j(n_j-1)$ commutes with $N$, and each mode loses particles at
rate $\gamma$: one jump between sectors with the Kraus operators $a_1,a_2$, restricted to $N\to N-1$.
Since $a_1^*a_1+a_2^*a_2=N$, Proposition 8.9(2) applies with $Q_{N-1,N}=\gamma N$, and from $N=2$
$$
p_2(t)=q^2,\qquad p_1(t)=2q(1-q),\qquad p_0(t)=(1-q)^2,\qquad q=e^{-\gamma t},
$$
a binomial law: each particle survives independently with probability $e^{-\gamma t}$, whatever the
hopping and the interaction do (check $p'=Qp$ by differentiation).

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
    assert min(w.min().item() for w in rho_s.eigenvalues()) > -1e-7          # positive up to the RK4 error
print(traj[-1].sector_probabilities())
```

RK4 is not exactly completely positive: the eigenvalues of $\rho(t)$ can be negative at the level of
the discretisation error ($10^{-7}$ here, with steps $0.02$); `lindblad_channel` is exact but treats
only jumps inside the sectors. Finally we check Proposition 8.8(2) for a random map between the
sectors, against the GKSL generator on $M_6$ built with numpy:

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
```

## 8.4 Post-selection and the no-jump evolution

Collect the jumps of the generator of Section 8.3 into $\Phi_{\mathrm{tot}}(\rho)=\sum_jg_jL_j\rho L_j^*+\sum_mg_m\Phi_m(\rho)$,
with $A=\Phi_{\mathrm{tot}}^*(1)=\sum_jg_jL_j^*L_j+\sum_mg_m\Phi_m^*(1)\ge0$ in $\mathcal M$, and
$H_{\mathrm{eff}}=H-\frac i2A$. By Proposition 8.8(1) and Theorem 8.6 (whose $S_n$ map $\mathcal M$ into
$\mathcal M$, because $H_{\mathrm{eff}}\in\mathcal M$),
$$
\rho(t)=\sum_{n\ge0}S_n(t)\rho,\qquad S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},
\qquad S_n(t)=\int_0^tS_0(t-s)\,\Phi_{\mathrm{tot}}\,S_{n-1}(s)\,ds .
$$

**Proposition 8.10.** For a state $\rho$ let $w_n(t)=\operatorname{Tr}S_n(t)\rho$. Then $w_n(t)\ge0$,
$\sum_nw_n(t)=1$, and $w_0$ is nonincreasing with $w_0'(t)=-\operatorname{Tr}\big(A\,S_0(t)\rho\big)$. For a
pure state $\rho=|\psi\rangle\langle\psi|$, $w_0(t)=\|e^{-itH_{\mathrm{eff}}}\psi\|^2$.

*Proof.* $S_n(t)$ is CP, hence positive, and $\sum_nw_n=\operatorname{Tr}\rho(t)=1$. With
$X=S_0(t)\rho\ge0$, $w_0'=\operatorname{Tr}(-iH_{\mathrm{eff}}X+iXH_{\mathrm{eff}}^*)=-\operatorname{Tr}(AX)$,
and $\operatorname{Tr}(AX)=\operatorname{Tr}(X^{1/2}AX^{1/2})\ge0$. The last claim is
$\operatorname{Tr}|\phi\rangle\langle\phi|=\|\phi\|^2$. $\square$

> **Physical meaning.** Let a detector record every jump — every emitted photon, every lost atom.
> Then $w_n(t)$ is the probability of exactly $n$ clicks up to time $t$, $S_n(t)\rho/w_n(t)$ the state
> given $n$ clicks, and $S_0(t)\rho/w_0(t)$ the state given **no click**: a null result also changes the
> state, through the non-Hermitian $H_{\mathrm{eff}}$. This counting-process interpretation [E. B. Davies,
> *Quantum Theory of Open Systems*, 1976; M. D. Srinivas and E. B. Davies, *Optica Acta* 28, 981 (1981)]
> is the basis of the quantum-trajectory (Monte Carlo wave-function) method [J. Dalibard, Y. Castin,
> K. Mølmer, *Phys. Rev. Lett.* 68, 580 (1992)]; here we only use its mathematical part.

For losses, counting the clicks is the same as measuring the charge:

**Proposition 8.11.** Let the sectors be labelled $N=0,\dots,C-1$, let there be no jumps inside the
sectors, let every jump lower the charge by one ($K^{dc}_{m,l}=0$ unless $d=c-1$), and let
$\rho=P_{N_0}\rho P_{N_0}$. Then $S_n(t)\rho=P_{N_0-n}\,\rho(t)\,P_{N_0-n}$ for $0\le n\le N_0$ and
$S_n(t)\rho=0$ for $n>N_0$. In particular the Lüders update of $\rho(t)$ by $P_{N_0}$ ("no particle lost")
is the no-jump state:
$$
\frac{P_{N_0}\rho(t)P_{N_0}}{p_{N_0}(t)}=\frac{e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*}}{w_0(t)},
\qquad p_{N_0}(t)=w_0(t).
$$

*Proof.* $S_0(s)$ maps operators supported in sector $N$ to operators supported in sector $N$, as
$H_{\mathrm{eff}}\in\mathcal M$; $\Phi_{\mathrm{tot}}$ maps sector $N\ge1$ to sector $N-1$ and sector $0$ to $0$.
By induction $S_n(s)\rho$ is supported in sector $N_0-n$ for all $s$, and is $0$ for $n>N_0$. Compressing
$\rho(t)=\sum_nS_n(t)\rho$ with $P_{N_0-n}$ keeps exactly the $n$-th term. $\square$

The library's `states.lueders_update(rho, P)` returns $P\rho P/\operatorname{Tr}(P\rho P)$ and
$\operatorname{Tr}(\rho P)$. Let the two modes now leak at different rates $\gamma_1>\gamma_2$; then
$A=\gamma_1n_1+\gamma_2n_2$ is not central, Proposition 8.9(2) does not apply, and post-selection on
$N=2$ gives a genuinely non-Hermitian evolution. We check Proposition 8.11 for $n=0$ and, by computing
the integral for $S_1$ with Simpson's rule, for $n=1$:

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

At $t=2$ the probability of no loss is $0.230$; without hopping it would be $e^{-2\gamma_1t}=0.091$, so
here the sector weights do depend on the Hamiltonian. The post-selected state is at trace distance
$0.135$ from the state of the closed system: not losing a particle is information.

## 8.5 Stationary states and relaxation

A state $\sigma$ is **stationary** if $e^{t\mathcal L}\sigma=\sigma$ for all $t\ge0$, equivalently
$\mathcal L(\sigma)=0$ (differentiate at $t=0$; conversely use the series).

**Proposition 8.12.** Every quantum dynamical semigroup on $\mathcal M$ has a stationary state.

*Proof.* Fix a state $\rho$ and let $\sigma_\tau=\frac1\tau\int_0^\tau T_s\rho\,ds$. It is a state (a limit of
convex combinations of states; the states form a closed set). By linearity,
$T_t\sigma_\tau-\sigma_\tau=\frac1\tau\big(\int_\tau^{\tau+t}-\int_0^t\big)T_s\rho\,ds$, whose trace norm is
at most $2t/\tau$, as $\|T_s\rho\|_1=1$. The states form a bounded closed set in a finite-dimensional
space, so some sequence $\sigma_{\tau_m}$, $\tau_m\to\infty$, converges to a state $\sigma$
(Bolzano–Weierstrass), and $T_t\sigma-\sigma=\lim_m(T_t\sigma_{\tau_m}-\sigma_{\tau_m})=0$. $\square$

**Corollary 8.13.** If there are no jumps between sectors, the stationary states span a space of
dimension at least $C$: every sector carries its own stationary state.

*Proof.* $\mathcal L$ maps $P_c\mathcal MP_c\cong M_{k_c}$ into itself (all of $H$, $L_j$ are block diagonal),
so the restriction of $e^{t\mathcal L}$ is a quantum dynamical semigroup on $M_{k_c}$ (Lemma 8.4(a)) with
a stationary state $\sigma_c$ by Proposition 8.12. The $\sigma_c$ have disjoint supports. $\square$

**Proposition 8.14 (spontaneous emission).** On $M_2$ with basis $|g\rangle,|e\rangle$ let
$H=\omega|e\rangle\langle e|$, one jump $L=|g\rangle\langle e|$ and rate $\gamma>0$. Then
$$
\rho_{ee}(t)=e^{-\gamma t}\rho_{ee}(0),\qquad\rho_{ge}(t)=e^{(i\omega-\gamma/2)t}\rho_{ge}(0),
\qquad e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}},
$$
where $\mathcal U_t(\rho)=e^{-itH}\rho e^{itH}$ and $\mathcal A_p$ is the amplitude damping channel with
parameter $p$ (Section 6.4). The spectrum of $\mathcal L$ is $\{0,-\gamma,-\gamma/2\pm i\omega\}$,
$|g\rangle\langle g|$ is the only stationary state, and
$\tfrac12\|\rho(t)-|g\rangle\langle g|\|_1=\big(\rho_{ee}(0)^2e^{-2\gamma t}+|\rho_{ge}(0)|^2e^{-\gamma t}\big)^{1/2}$.

*Proof.* Here $L^*L=|e\rangle\langle e|$ and $L\rho L^*=\rho_{ee}|g\rangle\langle g|$, so entrywise
$\rho_{ee}'=-\gamma\rho_{ee}$, $\rho_{gg}'=\gamma\rho_{ee}$ and $\rho_{ge}'=(i\omega-\gamma/2)\rho_{ge}$, which
gives the solution. In the basis $(\rho_{gg},\rho_{ge},\rho_{eg},\rho_{ee})$ the matrix of $\mathcal L$ is upper
triangular with diagonal $(0,-\gamma/2+i\omega,-\gamma/2-i\omega,-\gamma)$. $\mathcal L(\sigma)=0$ forces
$\sigma_{ee}=\sigma_{ge}=0$. The amplitude damping channel multiplies $\rho_{ee}$ by $1-p$ and $\rho_{ge}$ by
$\sqrt{1-p}$ (its Kraus operators $|g\rangle\langle g|+\sqrt{1-p}\,|e\rangle\langle e|$ and
$\sqrt p\,|g\rangle\langle e|$), and $\mathcal U_t$ multiplies $\rho_{ge}$ by $e^{i\omega t}$; both are TP, which
fixes $\rho_{gg}$. Finally $\rho(t)-|g\rangle\langle g|$ is traceless and self-adjoint with diagonal
$(-\rho_{ee}(t),\rho_{ee}(t))$, so its eigenvalues are $\pm(\rho_{ee}(t)^2+|\rho_{ge}(t)|^2)^{1/2}$. $\square$

The coherence decays at half the rate of the population, the relation $T_2=2T_1$ of Section 6.4; the
slowest nonzero rate, the **spectral gap** $\gamma/2$, controls the relaxation. We check the spectrum
with `lindblad_superoperator`, the channel and the trace distance with `lindblad_channel`, and
Corollary 8.13 for random jumps inside the sectors of $M_2\oplus M_3$:

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

For these random jumps the code finds exactly one stationary state per sector. Criteria for
uniqueness — for instance, a faithful stationary state together with $\{H,L_j,L_j^*\}'=\mathbb C1$
[A. Frigerio, *Commun. Math. Phys.* 63, 269 (1978)] — are not proved here; Exercise 4 shows that
uniqueness can fail even in a factor.

## 8.6 What to compare with

QuTiP (`sesolve`, `mesolve` for the Lindblad equation with adaptive ODE solvers, `mcsolve` for
quantum trajectories, `steadystate`, `liouvillian`), dynamiqs (JAX; batched and differentiable
Schrödinger and Lindblad solvers on GPUs), QuantumOptics.jl, and SciPy (`scipy.linalg.expm`, and
`scipy.sparse.linalg.expm_multiply` for $e^{t\mathcal L}\rho$ with a sparse vectorised generator) solve
the same equations on one Hilbert space, and for the problems above they must give the same numbers.
On $\mathcal M$ they can be used with $D\times D$ block-diagonal matrices, at the price of superoperators
of size $D^2\times D^2$ instead of $k_c^2\times k_c^2$ per sector. The library adds the block structure:
exact propagators and Lindblad channels sector by sector, batches of Hamiltonians and jump operators
(for disorder averages), jumps between sectors as `InterSectorChannel`s with the compression of
Proposition 8.8, and automatic differentiation. It has no adaptive step-size control (RK4 is fixed
step and only approximately positive), no trajectory solver, and `lindblad_channel` does not accept
jumps between sectors; for those, QuTiP or dynamiqs on $M_D$ are the reference.

## Exercises

1. (Heisenberg picture) Show that $\operatorname{Tr}(\mathcal L(\rho)A)=\operatorname{Tr}(\rho\,\mathcal L^*(A))$ for
   the generator of Section 8.3, with
   $\mathcal L^*(A)=i[H,A]+\sum_jg_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})+\sum_mg_m(\Phi_m^*(A)-\frac12\{\Phi_m^*(1),A\})$,
   that $\mathcal L^*(1)=0$, and that $\mathcal L^*(P_c)=0$ if there are no jumps between sectors.
2. (Exceptional point) In the sector $N=1$ of Section 8.4,
   $H_{\mathrm{eff}}=\begin{pmatrix}-i\gamma_1/2&J\\J&-i\gamma_2/2\end{pmatrix}$. Find its eigenvalues and show that
   for $J=|\gamma_1-\gamma_2|/4$ it is not diagonalisable.
3. Let $N_0$ bosons in $m$ modes each be lost at rate $\gamma$ (jumps $a_1,\dots,a_m$ restricted to
   $N\to N-1$). Show that $p_N(t)=\binom{N_0}Nq^N(1-q)^{N_0-N}$ with $q=e^{-\gamma t}$, for every Hamiltonian
   in $\mathcal M$.
4. (Pure dephasing) On $M_2$ let $H=\omega|e\rangle\langle e|$ and $L=\sigma_z=\operatorname{diag}(1,-1)$ with rate
   $\gamma$. Show that the populations are constant, $\rho_{ge}(t)=e^{(i\omega-2\gamma)t}\rho_{ge}(0)$, and that
   the stationary states span the two-dimensional space of diagonal matrices although $C=1$.

<details>
<summary>Solutions</summary>

1. By cyclicity $\operatorname{Tr}(L\rho L^*A)=\operatorname{Tr}(\rho L^*AL)$,
   $\operatorname{Tr}(\{L^*L,\rho\}A)=\operatorname{Tr}(\rho\{L^*L,A\})$ and
   $\operatorname{Tr}(-i[H,\rho]A)=\operatorname{Tr}(\rho\,i[H,A])$; for jumps between sectors,
   $\operatorname{Tr}(\Phi(\rho)A)=\operatorname{Tr}(\rho\Phi^*(A))$ is the definition of $\Phi^*$
   (Proposition 6.8). Then $\mathcal L^*(1)=\sum_jg_j(L_j^*L_j-L_j^*L_j)+\sum_mg_m(\Phi_m^*(1)-\Phi_m^*(1))=0$,
   which is trace preservation once more. Without jumps between sectors, $P_c$ commutes with $H$ and
   all $L_j$, so $L_j^*P_cL_j=L_j^*L_jP_c=\frac12\{L_j^*L_j,P_c\}$ and $\mathcal L^*(P_c)=0$; then
   $p_c'=\operatorname{Tr}(\rho\mathcal L^*(P_c))=0$, a second proof of Proposition 8.9(1).

2. $H_{\mathrm{eff}}=-\frac i4(\gamma_1+\gamma_2)1+\begin{pmatrix}-i\delta&J\\J&i\delta\end{pmatrix}$ with
   $\delta=(\gamma_1-\gamma_2)/4$; the second matrix has trace $0$ and determinant $\delta^2-J^2$, so the
   eigenvalues are $-\frac i4(\gamma_1+\gamma_2)\pm\sqrt{J^2-\delta^2}$: two decay rates for $J<|\delta|$,
   one rate and two frequencies for $J>|\delta|$. For $J=|\delta|>0$ the second matrix $Y$ satisfies
   $Y^2=(J^2-\delta^2)1=0$ but $Y\ne0$, so $H_{\mathrm{eff}}$ has the single eigenvalue
   $-\frac i4(\gamma_1+\gamma_2)$ and is not a multiple of $1$: it is not diagonalisable.

```python
lam = np.linalg.eigvals(H_eff[1])
delta = (g1 - g2) / 4
root = np.sqrt(J ** 2 - delta ** 2 + 0j)
assert np.allclose(sorted(lam, key=np.real), sorted(-0.25j * (g1 + g2) + np.array([root, -root]), key=np.real))
Y = np.array([[-1j * delta, delta], [delta, 1j * delta]])          # J = delta: the exceptional point
assert np.abs(Y @ Y).max() < 1e-15 and np.abs(Y).max() > 0
```

3. On the sector $N$, $\sum_ja_j^*a_j=N1$, and the only pair of sectors is $(N-1,N)$, so Proposition
   8.9(2) gives $p_N'=-\gamma Np_N+\gamma(N+1)p_{N+1}$. With $q'=-\gamma q$, the derivative of the binomial
   term is $-\gamma Np_N+\gamma(N_0-N)\binom{N_0}Nq^{N+1}(1-q)^{N_0-N-1}$, and
   $(N+1)\binom{N_0}{N+1}=(N_0-N)\binom{N_0}N$ shows that the second term is $\gamma(N+1)p_{N+1}$. At $t=0$
   it is the point mass at $N_0$; uniqueness is Lemma 8.1.

4. $L\rho L^*-\frac12\{L^*L,\rho\}=\sigma_z\rho\sigma_z-\rho$ multiplies the off-diagonal entries by $-2$ and
   kills the diagonal, and $-i[H,\rho]_{ge}=i\omega\rho_{ge}$. So $\mathcal L(\sigma)=0$ iff $\sigma_{ge}=0$ (for
   $\gamma>0$), and the spectrum is $\{0,0,-2\gamma\pm i\omega\}$. Since $H$ and $\sigma_z$ commute and
   $(\sigma_z\cdot\sigma_z)^2=\mathrm{id}$, $e^{t\mathcal L}(\rho)=\frac{1+e^{-2\gamma t}}2U_t\rho U_t^*
   +\frac{1-e^{-2\gamma t}}2U_t\sigma_z\rho\sigma_zU_t^*$ with $U_t=e^{-itH}$: a unital channel whose fixed
   points are the matrices commuting with its Kraus operators (Theorem 6.10), the diagonal ones.

```python
deph = dynamics.lindblad_superoperator(q2, H_at, [q2.from_blocks([torch.diag(torch.tensor([1.0, -1.0], dtype=f64))])],
                                       rates=[gamma])[0][0]
ev = torch.linalg.eigvals(deph)
assert int((ev.abs() < 1e-12).sum()) == 2
assert np.allclose(sorted(ev.tolist(), key=lambda z: (z.real, z.imag))[:2],
                   [-2 * gamma - 1j * omega, -2 * gamma + 1j * omega], atol=1e-12)
```

</details>

**Next:** Chapter 9 (in preparation) treats spin chains: sectors of magnetisation and momentum,
Lanczos and the Krylov time evolution used above at scale, and disorder.
