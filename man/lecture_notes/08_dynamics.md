**English** · [Español](es/08_dynamics.md) · [Français](fr/08_dynamics.md) · [Deutsch](de/08_dynamics.md) · [中文](zh/08_dynamics.md) · [日本語](ja/08_dynamics.md) · [Русский](ru/08_dynamics.md)

# 8. Dynamics

In quantum mechanics a closed system is modelled by a one-parameter group of unitaries, an open one — a
system coupled to an environment that is not observed — by a semigroup of channels. In an algebra with
sectors both have a feature that is invisible on a single matrix algebra: a Hamiltonian in
$`\mathcal M`$ can never change the sector, so the charges are conserved, while an open system can
jump between sectors, for instance by losing a particle. This chapter proves the basic facts about
both kinds of evolution, writes the Lindblad equation with jumps between sectors, and shows that "no
particle was lost" is a measurement of the charge whose post-selected state evolves by a non-Hermitian
Hamiltonian.

**You will learn** why sector weights are constants of motion, why the generators of
Gorini–Kossakowski–Sudarshan–Lindblad (GKSL) form produce channels, how jumps between sectors enter
the Lindblad equation and when the sector weights obey a classical rate equation, what the no-jump
evolution and post-selection mean, and how to find stationary states and relaxation rates. The numerical
checks of this chapter are in [Tutorial 8](../tutorial/08_dynamics.md).

Throughout, $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$ acts on
$`\mathcal H=\bigoplus_c\mathbb C^{k_c}`$, $`D=\sum_ck_c`$, $`P_c`$ are the minimal central projections,
states are given by their blunt densities $`\rho`$ with sector weights $`p_c=\mathrm{Tr}\,\rho_c`$
(Chapter 3), $`[A,B]=AB-BA`$, $`\{A,B\}=AB+BA`$, and $`\hbar=1`$. A **superoperator** is a linear map
on $`\mathcal M`$ or on $`M_D(\mathbb C)`$; all norms of superoperators are operator norms, hence
submultiplicative. For a linear map $`B`$ on a finite-dimensional space, $`e^{tB}=\sum_{n\ge0}t^nB^n/n!`$.

<a id="s-8-1"></a>

## 8.1 Closed systems

A **Hamiltonian** is a self-adjoint $`H\in\mathcal M`$. A vector evolves by the **Schrödinger
equation** $`i\psi'(t)=H\psi(t)`$, a density by the **von Neumann equation** $`\rho'(t)=-i[H,\rho(t)]`$.
Both are linear equations with constant coefficients on a finite-dimensional space, and the
following lemma is all the theory of differential equations we need.

<a id="r-8-1"></a>

**Lemma 8.1.** Let $`B`$ be a linear map on a finite-dimensional normed space $`X`$. Then
$`\frac d{dt}e^{tB}=Be^{tB}=e^{tB}B`$, $`e^{(s+t)B}=e^{sB}e^{tB}`$, $`\|e^{tB}\|\le e^{|t|\,\|B\|}`$, and for
every $`x_0\in X`$ the map $`x(t)=e^{tB}x_0`$ is the only differentiable $`x`$ on an interval
$`I\ni0`$ with $`x'=Bx`$ and $`x(0)=x_0`$.

*Proof.* The series converges absolutely, with the norm bound, and as a power series in $`t`$ it can
be differentiated term by term, which gives the first formula. Then
$`\frac d{dt}(e^{tB}e^{-tB})=Be^{tB}e^{-tB}-e^{tB}Be^{-tB}=0`$, so $`e^{tB}e^{-tB}=1`$. If $`x'=Bx`$, then
$`\frac d{dt}(e^{-tB}x(t))=e^{-tB}(x'-Bx)=0`$, so $`e^{-tB}x(t)=x_0`$ and $`x(t)=e^{tB}x_0`$. Finally
$`t\mapsto e^{(s+t)B}x_0`$ and $`t\mapsto e^{tB}e^{sB}x_0`$ solve $`y'=By`$ with the same initial value. $`\square`$

<a id="r-8-2"></a>

**Proposition 8.2.** Let $`H=H^*\in\mathcal M`$ and $`U(t)=e^{-itH}`$, defined by the functional
calculus of Section 1.3 with $`f_t(\lambda)=e^{-it\lambda}`$.

1. $`U(t)=\bigoplus_ce^{-itH_c}`$ is a unitary in $`\mathcal M`$, $`U(s+t)=U(s)U(t)`$ and $`U'(t)=-iHU(t)`$;
   $`U(t)`$ equals the exponential series of $`-itH`$.
2. $`\psi(t)=U(t)\psi_0`$ is the only solution of $`i\psi'=H\psi`$ with $`\psi(0)=\psi_0`$, and
   $`\rho(t)=U(t)\rho_0U(t)^*`$ is the only solution in $`\mathcal M`$ of $`\rho'=-i[H,\rho]`$ with
   $`\rho(0)=\rho_0`$. If $`\rho_0`$ is a state, so is every $`\rho(t)`$.
3. (*Heisenberg picture*) $`\mathrm{Tr}(\rho(t)A)=\mathrm{Tr}(\rho_0A(t))`$ with
   $`A(t)=U(t)^*AU(t)`$, and $`A'=i[H,A]`$.

*Proof.* (1) Diagonalise every block, $`H_c=V_c\Lambda_cV_c^*`$; then $`V=\bigoplus_cV_c`$ diagonalises
$`H`$, and $`U(t)=Ve^{-it\Lambda}V^*`$ is block diagonal (by Theorem 1.5(a) the result does not depend on
the diagonalising unitary). Since $`|f_t|^2=1`$ and $`f_sf_t=f_{s+t}`$, Theorem 1.5(c) gives unitarity and
the group law. Differentiating $`e^{-it\Lambda}`$ entrywise gives $`U'=V(-i\Lambda)e^{-it\Lambda}V^*=-iHU`$, so
$`U(t)`$ is the series $`e^{-itH}`$ by Lemma 8.1 ($`X=M_D`$, $`B`$ the left multiplication by $`-iH`$).
(2) $`(U\psi_0)'=-iHU\psi_0`$, and $`(U\rho_0U^*)'=-iHU\rho_0U^*+U\rho_0U^*(iH)`$ because $`(U')^*=iU^*H`$.
Uniqueness is Lemma 8.1 with $`B=-iH`$ on $`\mathcal H`$, resp. $`B=-i[H,\cdot\,]`$ on $`\mathcal M`$. Positivity:
$`\langle x,U\rho_0U^*x\rangle=\langle U^*x,\rho_0U^*x\rangle\ge0`$; the trace by cyclicity. (3) Cyclicity of
the trace, and $`(U^*AU)'=iU^*HAU-iU^*AHU=i[H,A(t)]`$. $`\square`$

*In the library:* $`U(t)=e^{-itH}`$, block by block, is `dynamics.propagator`; the solutions of Proposition 8.2(2) are `dynamics.schrodinger` and `dynamics.von_neumann` — [tutorial, Section 8.1](../tutorial/08_dynamics.md#s-8-1).

<a id="r-8-3"></a>

**Theorem 8.3 (conservation of charge).** Let $`t\mapsto H(t)`$ take values in the self-adjoint
elements of $`\mathcal M`$, on an interval $`I`$.

1. If $`\psi:I\to\mathcal H`$ is differentiable with $`i\psi'=H(t)\psi`$, then every $`\|P_c\psi(t)\|`$ is constant.
2. If $`\rho:I\to\mathcal M`$ is differentiable with $`\rho'=-i[H(t),\rho]`$, then $`\mathrm{Tr}(Z\rho(t))`$ is
   constant for every central $`Z`$; in particular every sector weight $`p_c(t)=\mathrm{Tr}(P_c\rho(t))`$ is.
3. If $`H`$ is constant and $`A\in\mathcal M`$ commutes with $`H`$, then $`\mathrm{Tr}(A\rho(t))`$ is constant.

*Proof.* (1) The inner product is antilinear in the first argument, so
$`\frac d{dt}\langle\psi,P_c\psi\rangle=\langle-iH\psi,P_c\psi\rangle+\langle\psi,-iP_cH\psi\rangle =i\langle\psi,HP_c\psi\rangle-i\langle\psi,P_cH\psi\rangle=0`$, as $`P_c`$ is central. (2)
$`\frac d{dt}\mathrm{Tr}(Z\rho)=-i\mathrm{Tr}(Z[H,\rho])=-i\mathrm{Tr}([Z,H]\rho)=0`$ by
cyclicity. (3) The same computation with $`A`$. $`\square`$

> **Physical meaning** (physical interpretation, not a theorem of this course). A Hamiltonian in
> $`\mathcal M`$ commutes with the charges (particle number, magnetisation, a symmetry quantum number),
> so the probability of each charge is a constant of motion, even for driven, time-dependent Hamiltonians.

On $`\mathcal M`$ a unit vector $`\psi=\bigoplus_c\psi_c`$ defines the state of the mixture
$`\bigoplus_c|\psi_c\rangle\langle\psi_c|`$: every $`A\in\mathcal M`$ is block diagonal, so
$`\langle\psi,A\psi\rangle=\sum_c\langle\psi_c,A_c\psi_c\rangle=\mathrm{Tr}\big(A\bigoplus_c|\psi_c\rangle\langle\psi_c|\big)`$
(as in Chapter 3, Exercise 3). If $`\psi(t)`$ solves $`i\psi'=H\psi`$, then
$`\bigoplus_c|\psi_c(t)\rangle\langle\psi_c(t)|=U(t)\big(\bigoplus_c|\psi_c(0)\rangle\langle\psi_c(0)|\big)U(t)^*`$, as $`U(t)`$ is block
diagonal.

*Computation:* Theorem 8.3 and the agreement of the evolutions of $`\psi`$ and of $`\bigoplus_c|\psi_c\rangle\langle\psi_c|`$, with the exact solvers and with the Runge–Kutta solver `dynamics.schrodinger_rk4` for a driven Hamiltonian — [tutorial, Section 8.1](../tutorial/08_dynamics.md#s-8-1).

<a id="s-8-2"></a>

## 8.2 Quantum dynamical semigroups

CP and TP maps and channels are defined in Section 6.1 ([Chapter 6](06_quantum_channels.md)). Besides
Theorem 6.2 and Proposition 6.4 (a map $`A\mapsto\sum_lK_lAK_l^*`$, also between sectors, is CP) we need
three closure properties.

<a id="r-8-4"></a>

**Lemma 8.4.** (a) Sums, nonnegative multiples and compositions of CP maps are CP. Let $`\Phi`$ be CP on
$`\mathcal A=M_D`$ or on $`\mathcal A=\mathcal M`$, and let $`\mathcal N\subseteq\mathcal A`$ be $`\mathcal M`$ (if
$`\mathcal A=M_D`$) or one of the corners $`P_c\mathcal MP_c\cong M_{k_c}`$ (with unit $`P_c`$). If $`\Phi`$ maps
$`\mathcal N`$ into itself, then $`\Phi|_{\mathcal N}`$ is CP.
(b) If CP maps $`\Phi_m`$ converge to $`\Phi`$ (that is, $`\Phi_m(A)\to\Phi(A)`$ for all $`A`$), then $`\Phi`$ is CP.
(c) If $`s\mapsto\Psi(s)`$ is a continuous map from $`[a,b]`$ into the CP maps, then $`\int_a^b\Psi(s)\,ds`$ is CP.

*Proof.* (a) $`\Phi\mapsto\Phi_n`$ is linear and $`(\Phi\Psi)_n=\Phi_n\Psi_n`$; sums of positive matrices are
positive. For the restriction: $`M_n(\mathcal N)\subseteq M_n(\mathcal A)`$, and positivity is the same
condition in both. Indeed, by Section 6.1 an element $`[A_{ij}]`$ of $`M_n(\mathcal M)`$ is positive iff every
block matrix $`[(A_{ij})_c]`$ is; for $`\mathcal A=M_D`$ this is positivity of $`[A_{ij}]\in M_{nD}`$, because after a
permutation of the basis $`[A_{ij}]`$ is the direct sum of these blocks (all $`A_{ij}`$ are block diagonal);
and an element of $`M_n(P_c\mathcal MP_c)`$ has the single nonzero block $`[(A_{ij})_c]`$, whose positivity is
positivity in $`M_n(M_{k_c})`$. (b) The positive matrices form a closed set ($`\langle x,Xx\rangle\ge0`$ for each
$`x`$ is a closed condition), and $`(\Phi_m)_n(X)\to\Phi_n(X)`$. (c) The Riemann sums are CP by (a) and
converge to the integral (Lemma A.4(a)); apply (b). $`\square`$

**Definition.** A **quantum dynamical semigroup** on $`\mathcal M`$ is a family $`(T_t)_{t\ge0}`$ of
channels $`\mathcal M\to\mathcal M`$ with $`T_0=\mathrm{id}`$, $`T_{s+t}=T_sT_t`$ and $`t\mapsto T_t`$ continuous.

<a id="r-8-5"></a>

**Lemma 8.5.** Every continuous semigroup $`(T_t)_{t\ge0}`$ of superoperators ($`T_0=\mathrm{id}`$, $`T_{s+t}=T_sT_t`$) is
$`T_t=e^{t\mathcal L}`$ for exactly one superoperator $`\mathcal L`$, its **generator**, and
$`\mathcal L=\lim_{t\to0}(T_t-\mathrm{id})/t`$. The same holds, with the same proof, for continuous semigroups of
linear maps on any finite-dimensional space (for example of matrices).

*Proof.* (Integrals of continuous superoperator-valued functions are those of Lemma A.4.)
$`\|V_h-\mathrm{id}\|\le\max_{0\le s\le h}\|T_s-\mathrm{id}\|\to0`$ for $`V_h=\frac1h\int_0^hT_s\,ds`$ (Lemma A.4(c)), so
$`\det V_h\neq0`$ for some fixed small $`h`$, by continuity of the determinant. Now
$`T_tV_h=\frac1h\int_0^hT_{t+s}\,ds=\frac1h\int_t^{t+h}T_s\,ds`$ (Lemma A.4(b),(e)) is differentiable in $`t`$ with derivative
$`\frac1h(T_{t+h}-T_t)=T_t\frac{T_h-\mathrm{id}}h`$. Hence $`T_t`$ is differentiable and $`T_t'=T_t\mathcal L`$
with $`\mathcal L=\frac{T_h-\mathrm{id}}hV_h^{-1}`$, and $`T_t=e^{t\mathcal L}`$ by Lemma 8.1 (with $`B`$ the
right composition by $`\mathcal L`$). Uniqueness: $`\mathcal L=T'(0)`$. $`\square`$

<a id="r-8-6"></a>

**Theorem 8.6 (GKSL generators).** Let $`H=H^*`$ and $`L_1,\dots,L_J`$ be in $`M_D`$, and
```math
\mathcal L(\rho)=-i[H,\rho]+\sum_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
```
Then $`(e^{t\mathcal L})_{t\ge0}`$ is a quantum dynamical semigroup on $`M_D`$. More precisely, put
$`H_{\mathrm{eff}}=H-\frac i2\sum_jL_j^*L_j`$, $`\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)`$
and $`\Phi(\rho)=\sum_jL_j\rho L_j^*`$, so that $`\mathcal L=\mathcal L_0+\Phi`$, and define
```math
S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},\qquad
S_n(t)=\int_0^tS_0(t-s)\,\Phi\,S_{n-1}(s)\,ds\quad(n\ge1).
```
Then every $`S_n(t)`$ is CP, and $`e^{t\mathcal L}=\sum_{n\ge0}S_n(t)`$ (the **Dyson expansion**), with
$`\|S_n(t)\|\le e^{t\|\mathcal L_0\|}(t\|\Phi\|)^n/n!`$.

*Proof.* *Step 1.* $`X(t)=e^{-itH_{\mathrm{eff}}}\rho e^{itH_{\mathrm{eff}}^*}`$ satisfies $`X'=\mathcal L_0(X)`$
and $`X(0)=\rho`$, so $`S_0(t)=e^{t\mathcal L_0}`$ by Lemma 8.1; it is CP with the single Kraus operator
$`e^{-itH_{\mathrm{eff}}}`$. *Step 2.* Write $`a=\|\mathcal L_0\|`$, $`b=\|\Phi\|`$. By induction the $`S_n`$
are continuous and
$`\|S_n(t)\|\le\int_0^te^{(t-s)a}\,b\,e^{sa}\frac{(bs)^{n-1}}{(n-1)!}ds=e^{ta}\frac{(bt)^n}{n!}`$, so the
series converges uniformly on every $`[0,T]`$. *Step 3.* Since
$`S_n(t)=S_0(t)\int_0^te^{-s\mathcal L_0}\Phi S_{n-1}(s)\,ds`$, the product rule and the fundamental
theorem of calculus (Lemma A.4(e)) give $`S_n'=\mathcal L_0S_n+\Phi S_{n-1}`$ ($`n\ge1`$) and $`S_0'=\mathcal L_0S_0`$. By the
bounds the series of derivatives converges uniformly on $`[0,T]`$, so $`Y(t)=\sum_nS_n(t)`$ is
differentiable with $`Y'=\mathcal L_0Y+\Phi Y=\mathcal LY`$, $`Y(0)=\mathrm{id}`$, and $`Y(t)=e^{t\mathcal L}`$ by
Lemma 8.1. *Step 4.* $`\Phi`$ is CP (Theorem 6.2), so by induction the integrand
$`s\mapsto S_0(t-s)\Phi S_{n-1}(s)`$ is a continuous family of CP maps and $`S_n(t)`$ is CP (Lemma 8.4);
so are the partial sums and their limit $`e^{t\mathcal L}`$. *Step 5.* By cyclicity,
$`\mathrm{Tr}\,\mathcal L(\rho)=\sum_j\big(\mathrm{Tr}(L_j\rho L_j^*)-\mathrm{Tr}(L_j^*L_j\rho)\big)=0`$,
hence $`\frac d{dt}\mathrm{Tr}\,e^{t\mathcal L}\rho=\mathrm{Tr}\,\mathcal L(e^{t\mathcal L}\rho)=0`$: every
$`e^{t\mathcal L}`$ is TP. The semigroup law and continuity are Lemma 8.1. $`\square`$

*Computation:* Theorem 8.6 for arbitrary jumps on $`M_3`$, including the complete positivity of the jump part $`\sum_{n\ge1}S_n(t)`$ — [tutorial, Section 8.2](../tutorial/08_dynamics.md#s-8-2). *In the library:* $`e^{t\mathcal L}`$ is `dynamics.lindblad_channel`, its matrix $`\mathcal L`$ is `dynamics.lindblad_superoperator` — [tutorial, Section 8.2](../tutorial/08_dynamics.md#s-8-2).

<a id="r-8-7"></a>

**Theorem 8.7 (Gorini–Kossakowski–Sudarshan; Lindblad).** (1) Conversely, the generator of every
quantum dynamical semigroup on $`M_D`$ has the form of Theorem 8.6, with $`J\le D^2-1`$ jumps. One can choose
$`\mathrm{Tr}\,H=0`$, $`\mathrm{Tr}\,L_j=0`$ and $`\mathrm{Tr}(L_j^*L_k)=0`$ for $`j\neq k`$; with these
normalisations $`H`$ and $`\Phi=\sum_jL_j\cdot L_j^*`$ are determined by $`\mathcal L`$, and the $`L_j`$ up to
$`L'_j=\sum_ku_{jk}L_k`$ with a unitary $`u`$ (after padding the shorter family with zeros).

(2) The generator of every quantum dynamical semigroup on $`\mathcal M=\bigoplus_cM_{k_c}`$ is
```math
\mathcal L(\rho)=\Phi(\rho)-\tfrac12\{\Phi^*(1),\rho\}-i[H,\rho]
```
with a CP map $`\Phi:\mathcal M\to\mathcal M`$ (Kraus operators between the sectors, Proposition 6.4) and
$`H=H^*\in\mathcal M`$: the form of Section 8.3, with one jump $`\Phi_1=\Phi`$ between sectors and rate $`1`$.

*Proof.* (1) is Theorem D.4 with the uniqueness statements of Theorem D.5(b),(c), and (2) is Theorem D.7,
both proved in [Appendix D](D_gksl.md). $`\square`$

The key step (Theorem D.3) is that the Choi matrix of the generator of a semigroup of CP maps is positive on
the orthogonal complement of the maximally entangled vector; the semigroup property alone, with positive
maps $`T_t`$, would not suffice (the example after Theorem D.3). For (2) the maps $`T_t\circ E`$ with the
conditional expectation $`E(X)=\sum_cP_cXP_c`$ are CP on $`M_D`$; the Choi matrix of the derivative
$`\mathcal L\circ E`$ at $`t=0`$ is positive on the orthogonal complement of the vectors $`(1\otimes P_c)\Omega`$, which
gives a decomposition of $`\mathcal L\circ E`$ (Lemma D.2), and this is compressed back to $`\mathcal M`$ by $`E`$.
In the Heisenberg picture (Corollary D.8) the generator of a continuous semigroup of unital CP maps on $`\mathcal M`$ is $`A\mapsto\Psi(A)+K^*A+AK`$ with $`\Psi`$ CP and $`K\in\mathcal M`$. Historically, (1) is
due to V. Gorini, A. Kossakowski, E. C. G. Sudarshan [*J. Math. Phys.* 17, 821 (1976)] and G. Lindblad
[*Commun. Math. Phys.* 48, 119 (1976)], and the Heisenberg-picture form, for norm-continuous semigroups on
general von Neumann algebras, to E. Christensen and D. E. Evans [*J. London Math. Soc.* 20, 358 (1979)].

> **Physical meaning** (physical interpretation, not a theorem of this course). A system that interacts
> weakly with a large environment whose correlations decay quickly (the Born–Markov approximation) is
> expected to lose memory of its past: its state at time $`s+t`$ is then obtained from the state at time $`s`$
> by a channel that depends only on $`t`$. *If* this is so, Theorem 8.7 says that the dynamics is a
> Hamiltonian part plus "jumps" $`L_j`$ (photon emission, particle loss, dephasing) at rates fixed by the
> normalisation of $`L_j`$.

<a id="s-8-3"></a>

## 8.3 Lindblad equations in an algebra with sectors

On $`\mathcal M`$ there are two kinds of jumps. A jump **inside the sectors** is an operator
$`L_j\in\mathcal M`$. A jump **between sectors** is a CP map $`\Phi:\mathcal M\to\mathcal M`$ with Kraus
operators $`K^{dc}_l:\mathbb C^{k_c}\to\mathbb C^{k_d}`$ (Proposition 6.4, an `InterSectorChannel`):
$`\Phi(\rho)_d=\sum_{c,l}K^{dc}_l\rho_cK^{dc*}_l`$. By Proposition 6.8 its blunt dual satisfies
$`\Phi^*(1)_c=\sum_{d,l}K^{dc*}_lK^{dc}_l`$, a block-diagonal element of $`\mathcal M`$. With rates
$`g_j,g_m>0`$ the generator of `lindblad_rhs` is
```math
\mathcal L(\rho)=-i[H,\rho]+\sum_jg_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big)
+\sum_mg_m\Big(\Phi_m(\rho)-\tfrac12\{\Phi_m^*(1),\rho\}\Big).
```
By Proposition 8.8(1) below it generates a quantum dynamical semigroup on $`\mathcal M`$; conversely, by
Theorem 8.7(2) (Theorem D.7), every quantum dynamical semigroup on $`\mathcal M`$ has a generator of this form.

*In the library:* a jump between sectors is an `InterSectorChannel`, and $`\mathcal L(\rho)`$ is `dynamics.lindblad_rhs` — [tutorial, Section 8.3](../tutorial/08_dynamics.md#s-8-3).

<a id="r-8-8"></a>

**Proposition 8.8.** Let $`J^{dc}_l\in M_D`$ be the matrix with the block $`K^{dc}_l`$ in position $`(d,c)`$
and zeros elsewhere.

1. $`\mathcal L`$ maps $`\mathcal M`$ into $`\mathcal M`$ and is the restriction to $`\mathcal M`$ of the GKSL generator
   $`\hat{\mathcal L}`$ on $`M_D`$ with jumps $`\sqrt{g_j}L_j`$ and $`\sqrt{g_m}J^{dc}_{m,l}`$ (all $`m,l,d,c`$). Hence
   $`(e^{t\mathcal L})_{t\ge0}`$ is a quantum dynamical semigroup on $`\mathcal M`$.
2. Let instead $`\tilde{\mathcal L}`$ be the GKSL generator on $`M_D`$ in which the blocks of one Kraus
   index are summed into one jump, $`L_{m,l}=\sum_{d,c}J^{dc}_{m,l}`$ (with the same Hamiltonian, the same
   jumps $`\sqrt{g_j}L_j`$ and the jumps $`\sqrt{g_m}L_{m,l}`$). Then
   $`\mathcal L(\rho)=E\big(\tilde{\mathcal L}(\rho)\big)`$ for $`\rho\in\mathcal M`$, where $`E(X)=\sum_dP_dXP_d`$ is
   the conditional expectation of $`M_D`$ onto $`\mathcal M`$; in general $`\tilde{\mathcal L}(\rho)\notin\mathcal M`$.
3. If every $`L_{m,l}`$ has at most one nonzero block in each block row and in each block column (for
   instance if it maps sector $`c`$ into sector $`c-1`$ only, as the annihilation operators of the example below),
   then $`\tilde{\mathcal L}(\rho)=\mathcal L(\rho)`$ for $`\rho\in\mathcal M`$.

*Proof.* Since $`\rho=\sum_cP_c\rho P_c`$ and $`J^{dc}P_{c'}=\delta_{cc'}J^{dc}`$,
$`J^{dc}\rho J^{dc*}`$ is the block $`K^{dc}\rho_cK^{dc*}`$ in position $`(d,d)`$, and $`J^{dc*}J^{dc}`$ is
$`K^{dc*}K^{dc}`$ in position $`(c,c)`$. Summing over $`d,l`$ gives $`\Phi(\rho)`$ and $`\Phi^*(1)`$, which proves the
identity in (1), and all terms lie in $`\mathcal M`$. So $`\hat{\mathcal L}`$ and all its powers map $`\mathcal M`$
into $`\mathcal M`$, $`e^{t\hat{\mathcal L}}|_{\mathcal M}=e^{t\mathcal L}`$, and this is a channel of $`\mathcal M`$ by
Theorem 8.6 and Lemma 8.4(a). (2) $`L\rho L^*=\sum_{d,d',c}J^{dc}\rho J^{d'c*}`$ for $`L=\sum_{d,c}J^{dc}`$, and
$`E`$ keeps the terms $`d=d'`$. For $`B\in\mathcal M`$, $`P_dXBP_d=P_dXP_dB`$, so $`E(XB)=E(X)B`$ and
$`E(BX)=BE(X)`$; thus $`E\{L^*L,\rho\}=\{E(L^*L),\rho\}`$, and $`E(L^*L)=\sum_cP_cL^*LP_c=\sum_{c,d}J^{dc*}J^{dc}`$
because $`J^{dc*}J^{d'c}=0`$ for $`d\ne d'`$. Finally $`E(-i[H,\rho])=-i[H,\rho]`$, and the terms of the jumps
$`\sqrt{g_j}L_j\in\mathcal M`$ lie in $`\mathcal M`$. (3) Let $`L=\sum_{(d,c)\in S}J^{dc}`$, where the set $`S`$ of pairs
contains at most one pair with a given $`d`$ and at most one with a given $`c`$. For $`\rho\in\mathcal M`$,
$`J^{dc}\rho J^{d'c'*}=J^{dc}P_c\rho P_{c'}J^{d'c'*}`$ vanishes unless $`c=c'`$, and then $`d=d'`$; and
$`J^{dc*}J^{d'c'}`$ vanishes unless $`d=d'`$, and then $`c=c'`$. So $`L\rho L^*`$ and $`L^*L`$ are block diagonal, all
terms of $`\tilde{\mathcal L}(\rho)`$ lie in $`\mathcal M`$, and $`\tilde{\mathcal L}(\rho)=E(\tilde{\mathcal L}(\rho))=\mathcal L(\rho)`$
by (2). $`\square`$

> **Physical meaning** (physical interpretation, not a theorem of this course). Part (2) says that the
> Lindblad equation on $`\mathcal M`$ is what remains of an equation on $`M_D`$ when the coherences between
> sectors that the jumps create are erased at once — the content of a superselection rule. For particle
> loss nothing is erased (part (3)): the usual Lindblad equation with the annihilation operators as jumps
> is already an equation on $`\mathcal M`$.

*Computation:* Proposition 8.8(2) for a random map between the sectors and Proposition 8.8(3) for particle loss, against GKSL generators on $`M_6`$ — [tutorial, Section 8.3](../tutorial/08_dynamics.md#s-8-3).

<a id="r-8-9"></a>

**Proposition 8.9 (sector weights).** Let $`\rho(t)=e^{t\mathcal L}\rho_0`$ and $`p(t)`$ its sector weights.

1. If there are no jumps between sectors ($`K^{dc}_{m,l}=0`$ for all $`m,l`$ and $`d\ne c`$), every $`p_c(t)`$ is constant.
2. If for every $`m`$ and every pair $`d\ne c`$, $`\sum_lK^{dc*}_{m,l}K^{dc}_{m,l}=r_{m,dc}1_{k_c}`$ with numbers
   $`r_{m,dc}`$, then $`p'=Qp`$ with $`Q_{dc}=\sum_mg_mr_{m,dc}`$ for $`d\ne c`$ and $`Q_{cc}=-\sum_{d\ne c}Q_{dc}`$.
   Hence $`p(t)=e^{tQ}p(0)`$ does not depend on $`H`$, on the $`L_j`$ or on the states inside the sectors.

*Proof.* $`p_d'=\mathrm{Tr}(P_d\mathcal L(\rho))`$ by Lemma 8.1. The Hamiltonian term vanishes as in Theorem 8.3.
For $`L_j\in\mathcal M`$, $`\mathrm{Tr}(P_dL_j\rho L_j^*)-\mathrm{Tr}(P_dL_j^*L_j\rho) =\mathrm{Tr}(L_j^*P_dL_j\rho)-\mathrm{Tr}(L_j^*L_jP_d\rho)=0`$, as $`P_d`$ is central. For a jump
$`\Phi`$ between sectors, $`\mathrm{Tr}(P_d\Phi(\rho))=\sum_c\mathrm{Tr}\big(\sum_lK^{dc*}_lK^{dc}_l\rho_c\big)`$
and $`\mathrm{Tr}(P_d\tfrac12\{\Phi^*(1),\rho\})=\sum_{d'}\mathrm{Tr}\big(\sum_lK^{d'd*}_lK^{d'd}_l\rho_d\big)`$;
the terms $`c=d`$ and $`d'=d`$ cancel. This proves (1), and under the hypothesis of (2) the remaining terms
are $`\sum_{c\ne d}r_{dc}p_c-\sum_{d'\ne d}r_{d'd}p_d`$. Lemma 8.1 gives $`p(t)=e^{tQ}p(0)`$. $`\square`$

Under the hypothesis of (2), for $`d\ne c`$, $`r_{m,dc}`$ is the entry $`T_{dc}(\rho)`$ of the transition
matrix of $`\Phi_m`$ (defined before Proposition 6.12) for every state $`\rho`$ with $`p_c>0`$, since
$`\mathrm{Tr}\,\Phi_{m,dc}(\rho_c)=\mathrm{Tr}\big(\sum_lK^{dc*}_{m,l}K^{dc}_{m,l}\rho_c\big)=r_{m,dc}p_c`$. So, for one jump with rate
$`1`$, $`Q=R-\mathrm{diag}(\sum_dR_{dc})_c`$ with $`R`$ = `transition_matrix()` (the diagonal of $`R`$ cancels).

*In the library:* $`R`$ is `InterSectorChannel.transition_matrix` — [tutorial, Section 8.3](../tutorial/08_dynamics.md#s-8-3).

**Example: particle loss.** Two bosonic modes with at most two particles; the sectors are
$`N=0,1,2`$ with bases $`|00\rangle`$; $`|10\rangle,|01\rangle`$; $`|20\rangle,|11\rangle,|02\rangle`$, so
$`\mathcal M=M_1\oplus M_2\oplus M_3`$. The Bose–Hubbard dimer
$`H=J(a_1^*a_2+a_2^*a_1)+\frac U2\sum_jn_j(n_j-1)`$ commutes with $`N`$, and each mode loses particles at
rate $`\gamma`$: one jump between sectors with the Kraus operators $`a_1,a_2`$, restricted to $`N\to N-1`$.
Since $`a_1^*a_1+a_2^*a_2=N`$, Proposition 8.9(2) applies with $`Q_{N-1,N}=\gamma N`$, and from $`N=2`$
```math
p_2(t)=q^2,\qquad p_1(t)=2q(1-q),\qquad p_0(t)=(1-q)^2,\qquad q=e^{-\gamma t},
```
the binomial law with parameters $`2`$ and $`q`$ (the law of the number of survivors if each particle survives
independently with probability $`q`$), whatever the hopping and the interaction do: both sides solve
$`p'=Qp`$ (differentiate) with $`p(0)=(0,0,1)`$, so they agree by Lemma 8.1.

*Computation:* the example, with $`\Phi^*(1)=0\oplus1\oplus2`$ and $`p(t)`$ checked against both formulas — [tutorial, Section 8.3](../tutorial/08_dynamics.md#s-8-3).

<a id="s-8-4"></a>

## 8.4 Post-selection and the no-jump evolution

Collect the jumps of the generator of Section 8.3 into $`\Phi_{\mathrm{tot}}(\rho)=\sum_jg_jL_j\rho L_j^*+\sum_mg_m\Phi_m(\rho)`$,
with $`A=\Phi_{\mathrm{tot}}^*(1)=\sum_jg_jL_j^*L_j+\sum_mg_m\Phi_m^*(1)\ge0`$ in $`\mathcal M`$, and
$`H_{\mathrm{eff}}=H-\frac i2A`$. Apply Theorem 8.6 to the generator $`\hat{\mathcal L}`$ on $`M_D`$ of Proposition 8.8(1):
by the proof of Proposition 8.8 its jump part maps $`\mathcal M`$ into $`\mathcal M`$, where it equals $`\Phi_{\mathrm{tot}}`$, and
its sum $`\sum\hat L^*\hat L`$ over the jumps is $`A`$, so its $`H_{\mathrm{eff}}`$ is $`H-\frac i2A\in\mathcal M`$. Hence its $`S_n(t)`$ map
$`\mathcal M`$ into $`\mathcal M`$ (by induction on $`n`$), and for $`\rho\in\mathcal M`$
```math
\rho(t)=\sum_{n\ge0}S_n(t)\rho,\qquad S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},
\qquad S_n(t)=\int_0^tS_0(t-s)\,\Phi_{\mathrm{tot}}\,S_{n-1}(s)\,ds .
```

<a id="r-8-10"></a>

**Proposition 8.10.** For a state $`\rho`$ let $`w_n(t)=\mathrm{Tr}\,S_n(t)\rho`$. Then $`w_n(t)\ge0`$,
$`\sum_nw_n(t)=1`$, and $`w_0`$ is nonincreasing with $`w_0'(t)=-\mathrm{Tr}\big(A\,S_0(t)\rho\big)`$. For a
pure state $`\rho=|\psi\rangle\langle\psi|`$, $`w_0(t)=\|e^{-itH_{\mathrm{eff}}}\psi\|^2`$.

*Proof.* $`S_n(t)`$ is CP, hence positive, and $`\sum_nw_n=\mathrm{Tr}\,\rho(t)=1`$. With
$`X=S_0(t)\rho\ge0`$, $`w_0'=\mathrm{Tr}(-iH_{\mathrm{eff}}X+iXH_{\mathrm{eff}}^*)=-\mathrm{Tr}(AX)`$,
and $`\mathrm{Tr}(AX)=\mathrm{Tr}(X^{1/2}AX^{1/2})\ge0`$. The last claim is
$`\mathrm{Tr}|\phi\rangle\langle\phi|=\|\phi\|^2`$. $`\square`$

> **Physical meaning** (physical interpretation, not a theorem of this course). Let a detector record
> every jump — every emitted photon, every lost atom. Then $`w_n(t)`$ is the probability of exactly $`n`$ clicks
> up to time $`t`$, $`S_n(t)\rho/w_n(t)`$ the state
> given $`n`$ clicks, and $`S_0(t)\rho/w_0(t)`$ the state given **no click**: a null result also changes the
> state, through the non-Hermitian $`H_{\mathrm{eff}}`$. This counting-process interpretation [E. B. Davies,
> *Quantum Theory of Open Systems*, 1976; M. D. Srinivas and E. B. Davies, *Optica Acta* 28, 981 (1981)]
> is the basis of the quantum-trajectory (Monte Carlo wave-function) method [J. Dalibard, Y. Castin,
> K. Mølmer, *Phys. Rev. Lett.* 68, 580 (1992)]; here we only use its mathematical part.

For losses the Dyson terms are the sector components of $`\rho(t)`$ (in the language of the box: counting the
clicks is the same as measuring the charge). The **Lüders update** of a state $`\rho`$ by a projection $`P`$ with
$`\mathrm{Tr}(P\rho)>0`$ is the state $`P\rho P/\mathrm{Tr}(P\rho P)`$ (note $`\mathrm{Tr}(P\rho P)=\mathrm{Tr}(P\rho)`$).

*In the library:* the Lüders update is `states.lueders_update` — [tutorial, Section 8.4](../tutorial/08_dynamics.md#s-8-4).

<a id="r-8-11"></a>

**Proposition 8.11.** Let the sectors be labelled $`N=0,\dots,C-1`$, let there be no jumps inside the
sectors, let every jump lower the charge by one ($`K^{dc}_{m,l}=0`$ unless $`d=c-1`$), and let
$`\rho=P_{N_0}\rho P_{N_0}`$. Then $`S_n(t)\rho=P_{N_0-n}\,\rho(t)\,P_{N_0-n}`$ for $`0\le n\le N_0`$ and
$`S_n(t)\rho=0`$ for $`n>N_0`$. In particular $`p_{N_0}(t)=w_0(t)>0`$, and the Lüders update of $`\rho(t)`$ by $`P_{N_0}`$
("no particle lost") is the no-jump state:
```math
\frac{P_{N_0}\rho(t)P_{N_0}}{p_{N_0}(t)}=\frac{e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*}}{w_0(t)},
\qquad p_{N_0}(t)=w_0(t).
```

*Proof.* $`S_0(s)`$ maps operators supported in sector $`N`$ to operators supported in sector $`N`$, as
$`H_{\mathrm{eff}}\in\mathcal M`$; $`\Phi_{\mathrm{tot}}`$ maps sector $`N\ge1`$ to sector $`N-1`$ and annihilates sector $`0`$.
By induction $`S_n(s)\rho`$ is supported in sector $`N_0-n`$ for all $`s`$, and is $`0`$ for $`n>N_0`$. Compressing
$`\rho(t)=\sum_nS_n(t)\rho`$ with $`P_{N_0-n}`$ keeps exactly the $`n`$-th term. Finally $`w_0(t)>0`$, because
$`e^{-itH_{\mathrm{eff}}}`$ is invertible (Lemma 8.1), so $`S_0(t)\rho`$ is a nonzero positive matrix. $`\square`$

*Computation:* Proposition 8.11 for $`n=0`$ and $`n=1`$, for two modes that leak at different rates — [tutorial, Section 8.4](../tutorial/08_dynamics.md#s-8-4).

<a id="s-8-5"></a>

## 8.5 Stationary states and relaxation

A state $`\sigma`$ is **stationary** if $`e^{t\mathcal L}\sigma=\sigma`$ for all $`t\ge0`$, equivalently
$`\mathcal L(\sigma)=0`$ (differentiate at $`t=0`$; conversely use the series).

<a id="r-8-12"></a>

**Proposition 8.12.** Every quantum dynamical semigroup on $`\mathcal M`$ has a stationary state.

*Proof.* Fix a state $`\rho`$ and let $`\sigma_\tau=\frac1\tau\int_0^\tau T_s\rho\,ds`$. It is a state by Lemma A.4(d)
(the states form a closed convex set). By Lemma A.4(b),(e),
$`T_t\sigma_\tau-\sigma_\tau=\frac1\tau\big(\int_\tau^{\tau+t}-\int_0^t\big)T_s\rho\,ds`$, whose trace norm is
at most $`2t/\tau`$ by Lemma A.4(c), as $`\|T_s\rho\|_1=1`$. The states form a bounded closed set in a finite-dimensional
space, so some sequence $`\sigma_{\tau_m}`$, $`\tau_m\to\infty`$, converges to a state $`\sigma`$
(Bolzano–Weierstrass), and $`T_t\sigma-\sigma=\lim_m(T_t\sigma_{\tau_m}-\sigma_{\tau_m})=0`$. $`\square`$

<a id="r-8-13"></a>

**Corollary 8.13.** If there are no jumps between sectors (as in Proposition 8.9(1)), the stationary states span a space of
dimension at least $`C`$: every sector carries its own stationary state.

*Proof.* $`\mathcal L`$ maps $`P_c\mathcal MP_c\cong M_{k_c}`$ into itself ($`H`$ and the $`L_j`$ lie in $`\mathcal M`$, and
$`\Phi_m(X)=\sum_lK^{cc}_{m,l}X_cK^{cc*}_{m,l}`$ and $`\Phi_m^*(1)`$ are block diagonal), hence so does every $`e^{t\mathcal L}`$.
Its restriction is CP (Lemma 8.4(a)) and TP, so it is a quantum dynamical semigroup on $`M_{k_c}`$, with
a stationary state $`\sigma_c`$ by Proposition 8.12; as an element of $`\mathcal M`$, $`\sigma_c`$ is a stationary state of
$`e^{t\mathcal L}`$. The $`\sigma_c`$ lie in the different blocks, hence are linearly independent. $`\square`$

Theorem D.6 (the finite-dimensional case of a criterion of A. Frigerio) gives a sufficient condition for
uniqueness: if a GKSL generator on $`M_k`$ with data $`(H;L_1,\dots,L_r)`$ has a faithful stationary state $`\sigma>0`$
and the only matrices commuting with $`H`$ and all $`L_j`$, $`L_j^*`$ are the multiples of $`1`$, then $`\sigma`$ is the only
stationary state. Without jumps between sectors the restriction of $`\mathcal L`$ to $`P_c\mathcal MP_c\cong M_{k_c}`$
is the GKSL generator with data $`(H_c;\sqrt{g_1}L_{1,c},\sqrt{g_2}L_{2,c},\dots)`$, so the theorem applies sector by
sector. Exercise 4 shows that uniqueness can fail even in a factor: there the second hypothesis fails.

*Computation:* Corollary 8.13 and the hypotheses of Theorem D.6 for random jumps inside the sectors of $`M_2\oplus M_3`$ — [tutorial, Section 8.5](../tutorial/08_dynamics.md#s-8-5).

<a id="r-8-14"></a>

**Proposition 8.14 (spontaneous emission).** On $`M_2`$ with basis $`|g\rangle,|e\rangle`$ let
$`H=\omega|e\rangle\langle e|`$, one jump $`L=|g\rangle\langle e|`$ and rate $`\gamma>0`$. Then
```math
\rho_{ee}(t)=e^{-\gamma t}\rho_{ee}(0),\qquad\rho_{ge}(t)=e^{(i\omega-\gamma/2)t}\rho_{ge}(0),
\qquad e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}},
```
where $`\mathcal U_t(\rho)=e^{-itH}\rho e^{itH}`$ and $`\mathcal A_p`$ is the amplitude damping channel with
parameter $`p`$ (Section 6.4). The spectrum of $`\mathcal L`$ is $`\{0,-\gamma,-\gamma/2\pm i\omega\}`$,
$`|g\rangle\langle g|`$ is the only stationary state, and
$`\tfrac12\|\rho(t)-|g\rangle\langle g|\|_1=\big(\rho_{ee}(0)^2e^{-2\gamma t}+|\rho_{ge}(0)|^2e^{-\gamma t}\big)^{1/2}`$.

*Proof.* Here $`L^*L=|e\rangle\langle e|`$ and $`L\rho L^*=\rho_{ee}|g\rangle\langle g|`$, so entrywise
$`\rho_{ee}'=-\gamma\rho_{ee}`$, $`\rho_{gg}'=\gamma\rho_{ee}`$ and $`\rho_{ge}'=(i\omega-\gamma/2)\rho_{ge}`$, which
gives the solution. In the basis $`(\rho_{gg},\rho_{ge},\rho_{eg},\rho_{ee})`$ the matrix of $`\mathcal L`$ is upper
triangular with diagonal $`(0,-\gamma/2+i\omega,-\gamma/2-i\omega,-\gamma)`$. $`\mathcal L(\sigma)=0`$ forces
$`\sigma_{ee}=\sigma_{ge}=0`$. The amplitude damping channel multiplies $`\rho_{ee}`$ by $`1-p`$ and $`\rho_{ge}`$ by
$`\sqrt{1-p}`$ (its Kraus operators $`|g\rangle\langle g|+\sqrt{1-p}\,|e\rangle\langle e|`$ and
$`\sqrt p\,|g\rangle\langle e|`$), and $`\mathcal U_t`$ multiplies $`\rho_{ge}`$ by $`e^{i\omega t}`$; both are TP, which
fixes $`\rho_{gg}`$. Finally $`\rho(t)-|g\rangle\langle g|`$ is traceless and self-adjoint with diagonal
$`(-\rho_{ee}(t),\rho_{ee}(t))`$, so its eigenvalues are $`\pm(\rho_{ee}(t)^2+|\rho_{ge}(t)|^2)^{1/2}`$. $`\square`$

The coherence decays at half the rate of the population: with $`T_1=1/\gamma`$ and $`T_2=2/\gamma`$ this is the relation
$`\frac1{T_2}=\frac1{2T_1}+\frac1{T_\varphi}`$ of Section 6.4 without the dephasing term. The slowest nonzero decay rate,
the **spectral gap** $`\min\{-\mathrm{Re}\,\mu:\mu\ne0\text{ an eigenvalue of }\mathcal L\}=\gamma/2`$, is the rate of
relaxation here: by the last formula the
distance to $`|g\rangle\langle g|`$ is at most $`e^{-\gamma t/2}`$ (as $`|\rho_{ge}|^2\le\rho_{gg}\rho_{ee}`$ by positivity,
$`\rho_{ee}^2+|\rho_{ge}|^2\le\rho_{ee}\le1`$), and $`e^{\gamma t/2}`$ times it tends to $`|\rho_{ge}(0)|`$ as $`t\to\infty`$.

*Computation:* the spectrum of $`\mathcal L`$, the factorisation $`e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}}`$ and the trace distance of Proposition 8.14, with `dynamics.lindblad_superoperator` and `dynamics.lindblad_channel` — [tutorial, Section 8.5](../tutorial/08_dynamics.md#s-8-5).

<a id="s-8-6"></a>

## 8.6 What to compare with

*Computation:* the comparison with other software (QuTiP, dynamiqs, QuantumOptics.jl, SciPy) — [tutorial, Section 8.6](../tutorial/08_dynamics.md#s-8-6).

<a id="s-8-7"></a>

## 8.7 Time steps that are channels

A numerical solver replaces $`e^{t\mathcal L}`$ by $`n`$ steps of a map $`R_h`$, $`h=t/n`$. Every $`e^{h\mathcal L}`$ is a
channel of $`\mathcal M`$ (Proposition 8.8(1)); if every step is a channel too, the computed densities are states,
whatever $`h`$ is. For $`\rho'=\mathcal L(\rho)`$ the classical Runge–Kutta method (RK4) is the step
$`\rho\mapsto\rho+\frac h6(k_1+2k_2+2k_3+k_4)`$ with $`k_1=\mathcal L(\rho)`$, $`k_2=\mathcal L(\rho+\frac h2k_1)`$,
$`k_3=\mathcal L(\rho+\frac h2k_2)`$ and $`k_4=\mathcal L(\rho+hk_3)`$. Inserting the stages one into the next,
$`k_2=(\mathcal L+\frac h2\mathcal L^2)\rho`$, $`k_3=(\mathcal L+\frac h2\mathcal L^2+\frac{h^2}4\mathcal L^3)\rho`$ and
$`k_4=(\mathcal L+h\mathcal L^2+\frac{h^2}2\mathcal L^3+\frac{h^3}4\mathcal L^4)\rho`$, so
```math
R_h=T_4(h\mathcal L),\qquad T_4(z)=1+z+\frac{z^2}2+\frac{z^3}6+\frac{z^4}{24}.
```
$`R_h`$ is TP, because $`\mathrm{Tr}\,\mathcal L(X)=0`$ for every $`X`$ (Step 5 of the proof of Theorem 8.6, with
Proposition 8.8(1)), and it maps self-adjoint elements to self-adjoint ones, because $`\mathcal L(X^*)=\mathcal L(X)^*`$.
It differs from $`e^{h\mathcal L}`$ by the terms of order $`\ge5`$ of the exponential series, but it need not be positive.

<a id="r-8-15"></a>

**Proposition 8.15 (RK4 is not positive).** On $`M_2`$ with basis $`|g\rangle,|e\rangle`$ let $`\mathcal L`$ be the
generator of Proposition 8.14 ($`H=\omega|e\rangle\langle e|`$ with $`\omega\in\mathbb R`$, one jump $`|g\rangle\langle e|`$
with rate $`\gamma>0`$), $`\lambda=i\omega-\gamma/2`$, $`h>0`$, $`\tau=T_4(-h\gamma)`$ and $`a=T_4(h\lambda)`$.

1. $`R_h(\rho)=\begin{pmatrix}\rho_{gg}+(1-\tau)\rho_{ee}&a\rho_{ge}\\\bar a\rho_{eg}&\tau\rho_{ee}\end{pmatrix}`$, and the
   following are equivalent: $`R_h`$ is positive; $`R_h`$ is CP; $`\tau\le1`$ and $`|a|^2\le\tau`$.
2. $`\tau\le1`$ iff $`h\gamma\le x_*`$, where $`x_*=2.7853\ldots`$ is the only real root of $`x^3-4x^2+12x-24`$. For
   $`h\gamma>x_*`$ the state $`|e\rangle\langle e|`$ is mapped to $`\mathrm{diag}(1-\tau,\tau)`$, whose first entry is negative.
3. With $`u=\gamma^2+4\omega^2`$,
   ```math
   \tau-|a|^2=h^5\Big(-\frac{\gamma u(4\omega^2-3\gamma^2)}{384}+\frac{u^2(4\omega^2-5\gamma^2)}{4608}\,h
   +\frac{\gamma u^3}{9216}\,h^2-\frac{u^4}{147456}\,h^3\Big).
   ```
   Hence, if $`2\omega>\sqrt3\,\gamma`$, there is $`h_0>0`$ such that $`R_h`$ is not positive for $`0<h<h_0`$; if
   $`\omega=\gamma`$, $`R_h`$ is not positive for any $`h>0`$; if $`\omega=0`$, $`R_h`$ is CP iff $`h\gamma\le x_*`$.

*Proof.* (1) By the proof of Proposition 8.14, $`\mathcal L(|g\rangle\langle g|)=0`$,
$`\mathcal L(|e\rangle\langle e|)=\gamma(|g\rangle\langle g|-|e\rangle\langle e|)`$, $`\mathcal L(|g\rangle\langle e|)=\lambda|g\rangle\langle e|`$
and $`\mathcal L(|e\rangle\langle g|)=\bar\lambda|e\rangle\langle g|`$. By induction
$`\mathcal L^n(|e\rangle\langle e|)=(-\gamma)^n(|e\rangle\langle e|-|g\rangle\langle g|)`$ for $`n\ge1`$, so
$`R_h(|e\rangle\langle e|)=|e\rangle\langle e|+(\tau-1)(|e\rangle\langle e|-|g\rangle\langle g|)`$, while $`R_h`$ fixes
$`|g\rangle\langle g|`$ and multiplies $`|g\rangle\langle e|`$ by $`T_4(h\lambda)=a`$ and $`|e\rangle\langle g|`$ by
$`T_4(h\bar\lambda)=\bar a`$ ($`T_4`$ has real coefficients). This is the formula. If $`\tau\le1`$ and $`|a|^2\le\tau`$, the
Kraus operators $`|g\rangle\langle g|+\bar a|e\rangle\langle e|`$, $`\sqrt{1-\tau}\,|g\rangle\langle e|`$ and
$`\sqrt{\tau-|a|^2}\,|e\rangle\langle e|`$ give $`R_h`$ (compare the four entries), so $`R_h`$ is CP (Theorem 6.2), in
particular positive. Conversely let $`R_h`$ be positive. Then $`R_h(|e\rangle\langle e|)\ge0`$ gives $`1-\tau\ge0`$. For
$`\psi=c|g\rangle+s|e\rangle`$ with $`c,s>0`$, $`c^2+s^2=1`$, the matrix $`R_h(|\psi\rangle\langle\psi|)\ge0`$ has the
determinant $`s^2\big(c^2(\tau-|a|^2)+\tau(1-\tau)s^2\big)\ge0`$; letting $`s\to0`$ gives $`\tau-|a|^2\ge0`$.

(2) With $`x=h\gamma`$, $`\tau-1=\frac x{24}q(x)`$, $`q(x)=x^3-4x^2+12x-24`$. Since $`q'(x)=3x^2-8x+12>0`$ (its
discriminant is $`64-144<0`$), $`q`$ is strictly increasing, and $`q(2)=-8<0<3=q(3)`$; so $`q`$ has exactly one real root
$`x_*\in(2,3)`$, and $`\tau\le1`$ iff $`x\le x_*`$.

(3) $`\tau-|a|^2=T_4(-h\gamma)-T_4(h\lambda)T_4(h\bar\lambda)`$ is a polynomial in $`h`$ of degree $`8`$. The coefficient of
$`h^n`$ in $`T_4(h\lambda)T_4(h\bar\lambda)`$ is $`\sum\lambda^j\bar\lambda^l/(j!\,l!)`$ over $`j+l=n`$, $`0\le j,l\le4`$. For
$`n\le4`$ all such pairs occur, and the sum is $`(\lambda+\bar\lambda)^n/n!=(-\gamma)^n/n!`$ by the binomial theorem, the
coefficient of $`h^n`$ in $`T_4(-h\gamma)`$; so the coefficients of $`h^0,\dots,h^4`$ in $`\tau-|a|^2`$ vanish. For
$`n=5,\dots,8`$, in terms of $`\sigma=\lambda+\bar\lambda=-\gamma`$ and $`\pi=\lambda\bar\lambda=u/4`$, the sums are
$`\frac{\pi(\lambda^3+\bar\lambda^3)}{24}+\frac{\pi^2\sigma}{12}=\frac{\pi\sigma^3-\pi^2\sigma}{24}`$,
$`\frac{\pi^2(\lambda^2+\bar\lambda^2)}{48}+\frac{\pi^3}{36}=\frac{\pi^2\sigma^2}{48}-\frac{\pi^3}{72}`$, $`\frac{\pi^3\sigma}{144}`$ and
$`\frac{\pi^4}{576}`$ (with $`\lambda^2+\bar\lambda^2=\sigma^2-2\pi`$, $`\lambda^3+\bar\lambda^3=\sigma^3-3\pi\sigma`$), and
$`T_4(-h\gamma)`$ has no such terms; substituting $`\sigma`$ and $`\pi`$ gives the formula. If $`2\omega>\sqrt3\gamma`$, the
coefficient of $`h^5`$ is negative, so $`\tau-|a|^2<0`$ for all small $`h>0`$, and $`R_h`$ is not positive by (1). If
$`\omega=\gamma`$, then $`u=5\gamma^2`$, and with $`y=5h\gamma`$ the formula becomes
$`\tau-|a|^2=-\frac{5(h\gamma)^5}{147456}\,p(y)`$, $`p(y)=y^3-16y^2+32y+384`$. Now $`p(y)>0`$ for $`y>0`$:
$`p(y)=y(y-8)^2+(384-32y)`$ is a sum of two nonnegative terms, not both zero, for $`y\le12`$; for $`12\le y\le16`$,
$`y(y-8)^2\ge192`$ (it increases for $`y\ge8`$) and $`384-32y\ge-128`$; and $`p(y)=y^2(y-16)+32y+384>0`$ for $`y\ge16`$. So
$`\tau<|a|^2`$ for every $`h>0`$. If $`\omega=0`$, then $`u=\gamma^2`$ and, with $`x=h\gamma`$,
$`\tau-|a|^2=x^5\big(\frac1{128}-\frac{5x}{4608}+\frac{x^2}{9216}-\frac{x^3}{147456}\big)\ge
x^5\big(\frac1{128}-\frac{15}{4608}-\frac{27}{147456}\big)>0`$ for $`0<x\le3`$; with (1) and (2), $`R_h`$ is CP iff $`x\le x_*`$. $`\square`$

So for a weakly damped two-level system ($`2\omega>\sqrt3\,\gamma`$) the RK4 step is not positive however small the step
is, although neither the damping alone ($`\omega=0`$, by (3)) nor the rotation alone has this defect for small steps: for
$`\gamma=0`$ the proof of (1) applies with $`\tau=1`$, and $`|T_4(ix)|^2=1-\frac{x^6}{72}+\frac{x^8}{576}\le1`$ for $`x^2\le8`$
(expand). The violation is small, $`\tau-|a|^2=O(h^5)`$, but it is not a rounding effect, and it is not cured by
the trace preservation of $`R_h`$. Stability theory of RK4 (the region where $`|T_4(z)|\le1`$) does not see it either: here
$`h\lambda`$ and $`-h\gamma`$ lie in that region for small $`h`$.

*Computation:* Proposition 8.15 with `lindblad_evolve(..., method='rk4')`: the step as $`T_4(h\mathcal L)`$, a pure state
mapped to a matrix with a negative eigenvalue, and a Choi matrix with a negative eigenvalue for the particle loss of
Section 8.3 — [tutorial, Section 8.7](../tutorial/08_dynamics.md#s-8-7).

The following scheme is a channel for every step size. It uses the decomposition $`\mathcal L=\mathcal L_0+\Phi_{\mathrm{tot}}`$
of Section 8.4, $`\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)`$ with $`H_{\mathrm{eff}}=H-\frac i2A`$ and
$`A=\Phi_{\mathrm{tot}}^*(1)`$, and the dual $`\Theta^*`$ of a superoperator $`\Theta`$ on $`\mathcal M`$ for the trace,
$`\mathrm{Tr}(\Theta(X)Z)=\mathrm{Tr}(X\Theta^*(Z))`$. For superoperators $`\Theta`$ on $`\mathcal M`$ we write
$`\|\Theta\|_1=\max\{\|\Theta(X)\|_1:X\in\mathcal M,\ \|X\|_1\le1\}`$, with the trace norm of Section 1.4; this norm is
submultiplicative.

<a id="r-8-16"></a>

**Theorem 8.16 (a splitting whose steps are channels).** For $`h\ge0`$ let $`E_h=e^{-ihH_{\mathrm{eff}}/2}`$,
$`S_h(\rho)=E_h\rho E_h^*`$ and
```math
\Psi_h=S_h\circ\Big(\mathrm{id}+h\Phi_{\mathrm{tot}}+\frac{h^2}2\Phi_{\mathrm{tot}}^2\Big)\circ S_h,\qquad M_h=\Psi_h^*(1),\qquad
\tilde\Psi_h(\rho)=\Psi_h\big(M_h^{-1/2}\rho M_h^{-1/2}\big).
```

1. $`M_h\in\mathcal M`$ and $`M_h\ge e^{ihH_{\mathrm{eff}}^*}e^{-ihH_{\mathrm{eff}}}`$, which is positive definite; $`\tilde\Psi_h`$
   is a channel of $`\mathcal M`$ for every $`h\ge0`$.
2. (*Local error.*) There is a constant $`C`$, depending on $`\mathcal L`$, with $`\|\tilde\Psi_h-e^{h\mathcal L}\|_1\le Ch^3`$
   for $`0\le h\le1`$.
3. (*Global error.*) For every state $`\rho`$, every $`t>0`$ and every integer $`n\ge t`$, with $`h=t/n`$,
   $`\|\tilde\Psi_h^n(\rho)-e^{t\mathcal L}\rho\|_1\le Cth^2`$.

*Proof.* *Step 0.* If $`\Theta`$ is a positive TP map on $`\mathcal M`$ and $`Y=Y^*\in\mathcal M`$, then
$`\|\Theta(Y)\|_1\le\|Y\|_1`$. Indeed, write $`Y=Y_+-Y_-`$ with $`Y_\pm\ge0`$ the parts of the spectral decomposition with
positive and with negative eigenvalues. The singular values of $`Y`$ are the absolute values of its eigenvalues
($`|Y|=(Y^2)^{1/2}`$), so $`\|Y\|_1=\mathrm{Tr}\,Y_++\mathrm{Tr}\,Y_-`$, and by the triangle inequality and
$`\|P\|_1=\mathrm{Tr}\,P`$ for $`P\ge0`$ (Section 1.4),
$`\|\Theta(Y)\|_1\le\mathrm{Tr}\,\Theta(Y_+)+\mathrm{Tr}\,\Theta(Y_-)=\|Y\|_1`$. An arbitrary $`X\in\mathcal M`$ is
$`Y_1+iY_2`$ with the self-adjoint $`Y_1=\frac12(X+X^*)`$, $`Y_2=\frac1{2i}(X-X^*)`$ of trace norm $`\le\|X\|_1`$
(as $`\|X^*\|_1=\|X\|_1`$, Chapter 1, Exercise 2); hence $`\|\Theta\|_1\le2`$.

*Step 1: (1).* $`H_{\mathrm{eff}}\in\mathcal M`$ (Section 8.4), so $`E_h\in\mathcal M`$ (its exponential series has all terms
in the subspace $`\mathcal M`$). $`\Phi_{\mathrm{tot}}`$ maps $`\mathcal M`$ into $`\mathcal M`$ and is CP, with the Kraus operators
$`\sqrt{g_j}L_j`$ and $`\sqrt{g_m}K^{dc}_{m,l}`$ (Proposition 6.4). So $`\mathrm{id}+h\Phi_{\mathrm{tot}}+\frac{h^2}2\Phi_{\mathrm{tot}}^2`$
is CP (Lemma 8.4(a)), $`S_h`$ is CP with the Kraus operator $`E_h`$, and $`\Psi_h`$ is CP. By cyclicity of the trace,
$`S_h^*(Z)=E_h^*ZE_h`$, and $`\Psi_h^*=S_h^*\circ(\mathrm{id}+h\Phi_{\mathrm{tot}}^*+\frac{h^2}2(\Phi_{\mathrm{tot}}^*)^2)\circ S_h^*`$,
where $`\Phi_{\mathrm{tot}}^*(Z)=\sum_jg_jL_j^*ZL_j+\sum_mg_m\Phi_m^*(Z)`$ (Exercise 1) is positive (Proposition 6.8: the
Kraus operators of $`\Phi_m^*`$ are the $`K^{dc*}_{m,l}`$). With $`W=S_h^*(1)=E_h^*E_h\ge0`$, the middle factor maps $`W`$ to
$`W`$ plus a positive element, and $`S_h^*`$ preserves the order, so
$`M_h\ge S_h^*(W)=(E_h^2)^*E_h^2=e^{ihH_{\mathrm{eff}}^*}e^{-ihH_{\mathrm{eff}}}`$ ($`E_h^2=e^{-ihH_{\mathrm{eff}}}`$ by Lemma 8.1). This
matrix is positive definite: $`\langle x,(E_h^2)^*E_h^2x\rangle=\|E_h^2x\|^2>0`$ for $`x\ne0`$, because $`E_h^2`$ has the
inverse $`e^{ihH_{\mathrm{eff}}}`$ (Lemma 8.1). So all eigenvalues of the self-adjoint $`M_h\in\mathcal M`$ are positive, and
$`N_h=M_h^{-1/2}\in\mathcal M`$ (functional calculus, Section 1.3) satisfies $`N_hM_hN_h=1`$. The map
$`\rho\mapsto N_h\rho N_h`$ is CP, hence so is $`\tilde\Psi_h`$ (Lemma 8.4(a)), and
$`\mathrm{Tr}\,\tilde\Psi_h(\rho)=\mathrm{Tr}(N_h\rho N_hM_h)=\mathrm{Tr}(\rho N_hM_hN_h)=\mathrm{Tr}\,\rho`$.

*Step 2: Taylor coefficients.* Put $`\alpha=\frac12\mathcal L_0`$ and $`\beta=\Phi_{\mathrm{tot}}`$, superoperators on
$`\mathcal M`$. As in Step 1 of the proof of Theorem 8.6, $`X(s)=E_s\rho E_s^*`$ solves $`X'=\alpha(X)`$ (differentiate
$`E_s=e^{-isH_{\mathrm{eff}}/2}`$ by Lemma 8.1), so $`S_h=e^{h\alpha}`$ by Lemma 8.1, and $`\Psi_h=e^{h\alpha}(1+h\beta+\frac{h^2}2\beta^2)e^{h\alpha}`$, while $`e^{h\mathcal L}=e^{h(2\alpha+\beta)}`$.
Multiplying out the absolutely convergent series, $`\Psi_h=\sum_{n\ge0}h^nP_n`$ and $`e^{h\mathcal L}=\sum_{n\ge0}h^nQ_n`$ with
$`\sum_n\|P_n\|_1\le e^{2\|\alpha\|_1}(1+\|\beta\|_1+\frac12\|\beta\|_1^2)`$ and $`\sum_n\|Q_n\|_1\le e^{\|2\alpha+\beta\|_1}`$. Here
$`P_0=Q_0=\mathrm{id}`$, $`P_1=2\alpha+\beta=Q_1`$, and
$`P_2=\frac{\alpha^2}2+\alpha\beta+\alpha^2+\beta\alpha+\frac{\beta^2}2+\frac{\alpha^2}2=\frac12(2\alpha+\beta)^2=Q_2`$. Hence, for
$`0\le h\le1`$, $`\|\Psi_h-e^{h\mathcal L}\|_1\le C_1h^3`$ with $`C_1=\sum_{n\ge3}(\|P_n\|_1+\|Q_n\|_1)`$, and
$`\|\Psi_h\|_1\le K=\sum_n\|P_n\|_1`$.

*Step 3: (2).* Let $`x`$ be a unit eigenvector of $`M_h`$ with eigenvalue $`\mu`$. As $`M_h`$ is block diagonal, $`x`$ can be
chosen in one sector, so $`X=|x\rangle\langle x|\in\mathcal M`$, $`\|X\|_1=1`$. Since $`e^{h\mathcal L}`$ is TP,
$`\mu-1=\mathrm{Tr}(XM_h)-\mathrm{Tr}\,X=\mathrm{Tr}\,\Psi_h(X)-\mathrm{Tr}\,e^{h\mathcal L}(X)`$, and $`|\mathrm{Tr}\,Y|\le\|Y\|_1`$
(Section 1.4) gives $`|\mu-1|\le C_1h^3`$. Let $`h_1=\min\{1,(2C_1)^{-1/3}\}`$. For $`h\le h_1`$ all $`\mu\ge\frac12`$, and $`N_h`$ has
the eigenvalues $`\mu^{-1/2}\le\sqrt2`$ with $`|\mu^{-1/2}-1|=\frac{|1-\mu|}{\sqrt\mu(1+\sqrt\mu)}\le|\mu-1|`$, as
$`\sqrt\mu(1+\sqrt\mu)\ge2^{-1/2}(1+2^{-1/2})>1`$. The operator norm of a self-adjoint matrix is its largest absolute
eigenvalue, so $`\|N_h-1\|\le C_1h^3`$ and $`\|N_h\|\le\sqrt2`$. Now
```math
\tilde\Psi_h(X)-e^{h\mathcal L}(X)=\Psi_h\big((N_h-1)XN_h+X(N_h-1)\big)+\big(\Psi_h-e^{h\mathcal L}\big)(X),
```
and $`\|AB\|_1\le\|A\|\,\|B\|_1`$, $`\|BA\|_1\le\|A\|\,\|B\|_1`$ (Chapter 1, Exercise 2; the second by taking adjoints) give
$`\|\tilde\Psi_h-e^{h\mathcal L}\|_1\le\big((1+\sqrt2)K+1\big)C_1h^3`$ for $`h\le h_1`$. For $`h_1<h\le1`$ both maps are positive
and TP, so $`\|\tilde\Psi_h-e^{h\mathcal L}\|_1\le4\le4h^3/h_1^3`$ by Step 0. Take
$`C=\max\{((1+\sqrt2)K+1)C_1,\,4h_1^{-3}\}`$.

*Step 4: (3).* The sum
$`\tilde\Psi_h^n-e^{nh\mathcal L}=\sum_{k=0}^{n-1}\tilde\Psi_h^{\,n-1-k}\big(\tilde\Psi_h-e^{h\mathcal L}\big)e^{kh\mathcal L}`$
telescopes. Applied to $`\rho`$, the $`k`$-th term is $`\tilde\Psi_h^{\,n-1-k}(Y_k)`$ with
$`Y_k=(\tilde\Psi_h-e^{h\mathcal L})(e^{kh\mathcal L}\rho)`$, which is self-adjoint (CP maps preserve self-adjointness, by the Kraus
form) with $`\|Y_k\|_1\le Ch^3`$ by (2), as the state $`e^{kh\mathcal L}\rho`$ has trace norm $`1`$. By (1) and Step 0,
$`\|\tilde\Psi_h^{\,n-1-k}(Y_k)\|_1\le Ch^3`$, and the $`n`$ terms add up to at most $`nCh^3=Cth^2`$. $`\square`$

Step 4 is where the steps being channels pays off: the errors made in the single steps are never amplified, whatever $`n`$
is, and no stability condition on $`h`$ is needed. Expanding the definition, the Kraus operators of $`\tilde\Psi_h`$ are
$`E_h^2N_h`$, $`\sqrt h\,E_hJE_hN_h`$ and $`\frac h{\sqrt2}E_hJ'JE_hN_h`$, where $`J,J'`$ run over the Kraus operators of
$`\Phi_{\mathrm{tot}}`$: half a step of the no-jump evolution of Section 8.4, at most two jumps, the other half step, and the
correction $`N_h`$ of the probabilities, with $`\|N_h-1\|\le C_1h^3`$ by Step 3.

*In the library:* $`\tilde\Psi_h`$ is `lindblad_evolve(..., method='kraus')`; `method='exact'` computes $`e^{h\mathcal L}`$
itself from the matrix of $`\mathcal L`$ on $`\mathcal M`$ (of size $`\sum_ck_c^2`$, jumps between sectors included), whose steps
are channels by Proposition 8.8(1); the default `method='rk4'` is $`R_h`$ — [tutorial, Section 8.7](../tutorial/08_dynamics.md#s-8-7).

<a id="x-8"></a>

## Exercises

1. (Heisenberg picture) Show that $`\mathrm{Tr}(\mathcal L(\rho)A)=\mathrm{Tr}(\rho\,\mathcal L^*(A))`$ for
   the generator of Section 8.3, with
   $`\mathcal L^*(A)=i[H,A]+\sum_jg_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})+\sum_mg_m(\Phi_m^*(A)-\frac12\{\Phi_m^*(1),A\})`$,
   that $`\mathcal L^*(1)=0`$, and that $`\mathcal L^*(P_c)=0`$ if there are no jumps between sectors.
2. (Exceptional point) In the sector $`N=1`$ of Section 8.4,
   $`H_{\mathrm{eff}}=\begin{pmatrix}-i\gamma_1/2&J\\J&-i\gamma_2/2\end{pmatrix}`$. Find its eigenvalues and show that
   for $`J=|\gamma_1-\gamma_2|/4>0`$ it is not diagonalisable.
3. Let $`N_0`$ bosons in $`m`$ modes each be lost at rate $`\gamma`$ (jumps $`a_1,\dots,a_m`$ restricted to
   $`N\to N-1`$). Show that $`p_N(t)=\binom{N_0}Nq^N(1-q)^{N_0-N}`$ with $`q=e^{-\gamma t}`$, for every Hamiltonian
   in $`\mathcal M`$.
4. (Pure dephasing) On $`M_2`$ let $`H=\omega|e\rangle\langle e|`$ and $`L=\sigma_z=\mathrm{diag}(1,-1)`$ with rate
   $`\gamma`$. Show that the populations are constant, $`\rho_{ge}(t)=e^{(i\omega-2\gamma)t}\rho_{ge}(0)`$, and that
   the stationary states span the two-dimensional space of diagonal matrices although $`C=1`$.

<details>
<summary>Solutions</summary>

1. By cyclicity $`\mathrm{Tr}(L\rho L^*A)=\mathrm{Tr}(\rho L^*AL)`$,
   $`\mathrm{Tr}(\{L^*L,\rho\}A)=\mathrm{Tr}(\rho\{L^*L,A\})`$ and
   $`\mathrm{Tr}(-i[H,\rho]A)=\mathrm{Tr}(\rho\,i[H,A])`$; for jumps between sectors,
   $`\mathrm{Tr}(\Phi(\rho)A)=\mathrm{Tr}(\rho\Phi^*(A))`$ is the definition of $`\Phi^*`$
   (Section 6.3). Then $`\mathcal L^*(1)=\sum_jg_j(L_j^*L_j-L_j^*L_j)+\sum_mg_m(\Phi_m^*(1)-\Phi_m^*(1))=0`$,
   which is trace preservation once more. Without jumps between sectors, $`P_c`$ commutes with $`H`$ and
   all $`L_j`$, so $`L_j^*P_cL_j=L_j^*L_jP_c=\frac12\{L_j^*L_j,P_c\}`$, every $`\Phi_m`$ has only the
   blocks $`K^{cc}_{m,l}`$, so $`\Phi_m^*(P_c)=\Phi_m^*(1)P_c=\frac12\{\Phi_m^*(1),P_c\}`$, and $`\mathcal L^*(P_c)=0`$; then
   $`p_c'=\mathrm{Tr}(\rho\mathcal L^*(P_c))=0`$, a second proof of Proposition 8.9(1).

2. $`H_{\mathrm{eff}}=-\frac i4(\gamma_1+\gamma_2)1+\begin{pmatrix}-i\delta&J\\J&i\delta\end{pmatrix}`$ with
   $`\delta=(\gamma_1-\gamma_2)/4`$; the second matrix has trace $`0`$ and determinant $`\delta^2-J^2`$, so the
   eigenvalues are $`-\frac i4(\gamma_1+\gamma_2)\pm\sqrt{J^2-\delta^2}`$: two decay rates for $`J<|\delta|`$,
   one rate and two frequencies for $`J>|\delta|`$. For $`J=|\delta|>0`$ the second matrix $`Y`$ satisfies
   $`Y^2=(J^2-\delta^2)1=0`$ but $`Y\ne0`$, so $`H_{\mathrm{eff}}`$ has the single eigenvalue
   $`-\frac i4(\gamma_1+\gamma_2)`$ and is not a multiple of $`1`$: it is not diagonalisable.

3. On the sector $`N`$, $`\sum_ja_j^*a_j=N1`$, and the only pair of sectors is $`(N-1,N)`$, so Proposition
   8.9(2) gives $`p_N'=-\gamma Np_N+\gamma(N+1)p_{N+1}`$. With $`q'=-\gamma q`$, the derivative of the binomial
   term is $`-\gamma Np_N+\gamma(N_0-N)\binom{N_0}Nq^{N+1}(1-q)^{N_0-N-1}`$, and
   $`(N+1)\binom{N_0}{N+1}=(N_0-N)\binom{N_0}N`$ shows that the second term is $`\gamma(N+1)p_{N+1}`$. At $`t=0`$
   it is the point mass at $`N_0`$; uniqueness is Lemma 8.1.

4. $`L\rho L^*-\frac12\{L^*L,\rho\}=\sigma_z\rho\sigma_z-\rho`$ multiplies the off-diagonal entries by $`-2`$ and
   kills the diagonal, and $`-i[H,\rho]_{ge}=i\omega\rho_{ge}`$. So $`\mathcal L(\sigma)=0`$ iff $`\sigma_{ge}=0`$ (for
   $`\gamma>0`$), and the spectrum is $`\{0,0,-2\gamma\pm i\omega\}`$. With $`U_t=e^{-itH}`$,
   $`e^{t\mathcal L}(\rho)=\frac{1+e^{-2\gamma t}}2U_t\rho U_t^*+\frac{1-e^{-2\gamma t}}2U_t\sigma_z\rho\sigma_zU_t^*`$: since
   $`\sigma_z`$ commutes with $`U_t`$ and $`\sigma_z^2=1`$, differentiation shows that the right side $`R(t)`$ satisfies
   $`R'=-i[H,R]+\gamma(\sigma_zR\sigma_z-R)=\mathcal L(R)`$ and $`R(0)=\rho`$, so Lemma 8.1 applies. For $`t>0`$ this is a unital channel whose fixed
   points are the matrices commuting with its Kraus operators (Theorem 6.10), the diagonal ones.

</details>

*Computation:* Exercises 2 and 4 numerically — [tutorial, Exercises](../tutorial/08_dynamics.md#x-8).

**Next:** [Chapter 9](09_spin_chains.md) treats spin chains: sectors of magnetisation and momentum,
Lanczos and Krylov time evolution for large sectors, and disorder.
