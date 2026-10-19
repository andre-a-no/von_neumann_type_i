# Quantum channels as neural networks: a roadmap

*Design proposal for the next versions of `torch_vn_algebra` (current release: 1.1.0; 1.2.0 in preparation, see
Section 11). Status: discussion draft. Nothing described under "proposed" exists in the library yet, except where a
status column or Section 11 says otherwise; every code block marked **NOT IMPLEMENTED** is a sketch of an API, not
working code. References are listed at the end with full bibliographic data; all of them were checked, and items we
could not check were left out.*

**Contents.** 0 Scope · 1 Dictionary channel → network · 2 Parametrisations · 3 Expressivity and stability ·
4 Training problems as variational problems · 5 Training as a dynamics of the channel · 6 Architectures and skip
connections · 7 Where nonlinearity comes from · 8 Tomography as training · 9 Encoders and bottlenecks · 10 Proposed
API · 11 Plan, risks, measures of success · 12 Positioning · References.

## 0. Scope and conventions

The object of study is the **quantum channel**: a completely positive map between finite-dimensional von Neumann
algebras

```math
\Phi:\ \mathcal M=\bigoplus_{c=1}^{C}M_{k_c}\ \longrightarrow\ \mathcal N=\bigoplus_{d=1}^{D}M_{m_d},
\qquad \Phi(\rho)_d=\sum_{c}\sum_i K_i^{dc}\,\rho_c\,K_i^{dc*},\quad K_i^{dc}:\mathbb C^{k_c}\to\mathbb C^{m_d},
```

in the Schrödinger picture (on densities), trace preserving (TP) when $`\sum_{d,i}K_i^{dc*}K_i^{dc}=1_{k_c}`$ for every
$`c`$, unital in the Heisenberg picture. The input and output may have different sector structures. This covers the
library's `Channel` (sector preserving, Kraus tensor `(B, r, C, k, k)`), `InterSectorChannel` (Kraus tensor
`(B, r, D, C, m_max, k_max)`), `SectorJump` (sparse jumps), Lindblad semigroups (`dynamics.lindblad_channel`,
`lindblad_evolve`), conditional expectations (`inclusions.ConditionalExpectation`) and, as shown in Section 1,
instruments.

Neural-network machinery is the **tool**: layers, width and depth, parametrisations, backpropagation, training on data,
regularisation and the dynamics of training are used to represent, parametrise, learn and optimise channels. The
proposal does *not* turn the library into a general deep-learning framework.

**Marking.** From Section 5 on, claims are marked **[E]** established (with references), **[P]** proposal (our
reading, construction or planned experiment; "elementary" when a short proof is given or sketched), or **[N]** possibly
new. Every [N] carries a note on the literature check done so far; "[N?]" means the check is incomplete, and no [N]
item may be claimed as new in a publication before a full check.

## 1. Dictionary: channel → network

| channel notion | network notion | in the library now | missing |
|---|---|---|---|
| Stinespring isometry $`V:\mathcal H_{\rm in}\to\mathcal H_{\rm out}\otimes\mathbb C^r`$, then $`\operatorname{Tr}_{\mathbb C^r}`$ | linear (isometric) layer followed by a contraction / pooling over an environment index (Stinespring 1955) | `channels.random_channel` (Haar isometry); `composite.partial_trace_channel` | `optimize.StinespringParam` (sector preserving and inter-sector; unreleased, 1.2.0); an `nn.Module` layer |
| Kraus rank / environment dimension $`r`$ | width | `Channel.kraus_rank` | rank control (truncation, growth) during training |
| composition $`\Phi_L\circ\dots\circ\Phi_1`$ | depth | `Channel.compose`, `InterSectorChannel.compose` | a sequential container with algebra checks |
| $`e^{t\mathcal L}`$ of a GKSL generator | continuous depth (neural ODE, Chen et al. 2018; E 2017) | `lindblad_channel`, `lindblad_evolve`, differentiable | `optimize.LindbladParam` (unreleased, 1.2.0) |
| CPTP splitting step $`\Psi_h(M^{-1/2}\,\cdot\,M^{-1/2})`$ (`dynamics._kraus_step`) | a second-order residual block that stays inside the constraint set (Section 6) | `lindblad_evolve(method='kraus')` | exposure as a layer |
| conditional expectation $`E:\mathcal M\to\mathcal N\subset\mathcal M`$ | pooling onto a subalgebra (coarse-graining) | `inclusions.ConditionalExpectation`, `channels.center_expectation` | — (usable as fixed layers) |
| sector (block) structure | structured sparsity and weight sharing (strong/weak symmetry of Lindbladians: Buča and Prosen 2012) | the storage format itself; `decompose` turns a symmetry into an algebra | sparse lists of allowed pairs $`(d,c)`$ |
| Choi matrix $`J(\Phi)`$ | the weight tensor | `Channel.choi`, `superoperator`, `from_superoperator` | Choi-based losses |
| TP / unitality | normalisation constraints (as for stochastic matrices) | `trace_preservation_error`, `unitality_error` | constraint-preserving parametrisations |
| instrument, post-selection on a central projection | the source of nonlinearity (normalisation); instruments: Davies and Lewis (1970) | `measurement_channel`, `lueders_update`, `DensityMatrix.condition_on`, `born_probabilities` | `measurement.Instrument` (outcome in the centre; in progress, 1.2.0) |

Four entries need comment.

**Width and universality of one layer.** By Choi's theorem every CP map $`M_k\to M_m`$ has a Kraus form with at most
$`km`$ operators (the rank of its Choi matrix). Hence one Stinespring layer with $`r=\max_{d,c}k_cm_d`$ represents
every CP map $`\mathcal M\to\mathcal N`$, because the library's Kraus tensor shares $`r`$ over all pairs $`(d,c)`$.
Width bounds compose: the Kraus rank of $`\Phi_2\circ\Phi_1`$ is at most $`r_1r_2`$ (Kraus operators
$`K^{(2)}_jK^{(1)}_i`$). Depth therefore adds no expressive power for unrestricted layers. What depth does is build
structured families, such as products of local channels, sector-preserving steps or Trotter steps of a semigroup, out
of small layers.

**Sector structure as weight sharing.** A CPTP map $`M_k\to M_m`$ has $`k^2m^2-k^2`$ real parameters (a positive Choi
matrix minus the TP constraint). Channels between algebras have $`\sum_{d,c}k_c^2m_d^2-\sum_c k_c^2`$, and channels
that preserve every sector have $`\sum_c(k_c^4-k_c^2)`$. For $`n`$ qubits with U(1) sectors $`k_N=\binom nN`$ the
counts are as follows (computed exactly):

| $`n`$ | dense $`4^{2n}-4^n`$ | sector to sector, all pairs | sector preserving |
|---|---|---|---|
| 6 | 16 773 120 | 852 852 | 262 920 |
| 8 | 4 294 901 760 | 165 624 030 | 44 903 628 |

Restricting further to the pairs a physical noise model allows, such as $`N\to N`$ and $`N\to N\pm1`$ for local
emission and absorption, reduces the count again. This is the same mechanism as weight sharing in equivariant networks
(Cohen and Welling 2016; Kondor and Trivedi 2018; Finzi et al. 2021). Here the commutant is computed by `decompose`,
and the result is an algebra in which products, adjoints and functional calculus are available.

**Normalisation: the classical case is a softmax.** For commutative $`\mathcal M=\mathbb C^n`$ and
$`\mathcal N=\mathbb C^m`$ (all $`k_c=m_d=1`$) a channel is a column-stochastic matrix $`T_{dc}=\sum_i|K_i^{dc}|^2`$. The
Kraus normalisation $`K\mapsto KM^{-1/2}`$ with $`M=\sum_iK_i^*K_i`$, used by `_kraus_step`, becomes column
normalisation $`T_{dc}=w_{dc}/\sum_{d'}w_{d'c}`$. With $`w=e^{s}`$ this is exactly a softmax over the output index. On
the stacked Kraus column $`Z=[K_1;\dots;K_r]`$ the same normalisation is the polar retraction
$`V=Z(Z^*Z)^{-1/2}`$ onto the Stiefel manifold. So "softmax", "Kraus normalisation" and "polar retraction" are a
single construction at three levels of commutativity.

**Where nonlinearity enters.** Every channel is linear in $`\rho`$. Measurement followed by classical feed-forward,
$`\rho\mapsto\sum_j\Psi_j(\Phi_j(\rho))`$ for an instrument $`(\Phi_j)`$, is still a channel and still linear. A
genuine nonlinearity enters only in three ways:

1. renormalisation after post-selection, $`\rho\mapsto P\rho P/\operatorname{Tr}(P\rho)`$;
2. nonlinear functions of the output statistics in the loss;
3. several copies, as in $`\rho\mapsto\Phi(\rho\otimes\rho)`$.

The algebraic formulation is clean. An instrument with $`J`$ outcomes is a single channel
$`\mathcal M\to\mathcal N\otimes\mathbb C^J=\bigoplus_{j}\mathcal N`$. The outcome is recorded in the centre of the
output algebra, and post-selection is conditioning on a central projection. The existing `InterSectorChannel` can
already represent this, with $`J\cdot D`$ output sectors. Because the target of channel learning is itself linear, no
nonlinearity is needed to *represent* channels. It is needed only for tasks that are not channels, such as
classifying states. Section 7 treats the sources of nonlinearity in detail, with a fourth one (mean-field limits)
and their costs.

## 2. Parametrisations in which every parameter is a valid channel

| family | parametrisation | gradient path | status |
|---|---|---|---|
| CPTP, Kraus rank $`r`$ | free $`Z_c\in\mathbb C^{(\sum_d m_d r)\times k_c}`$ per input sector; $`V_c=Z_c(Z_c^*Z_c)^{-1/2}`$ (polar), or QR, or Cayley updates (Wen and Yin 2013); Riemannian view in Absil et al. (2008) | inverse square root by `functional.spectral_apply`: Daleckii–Krein backward, finite at degenerate spectra | `optimize.StinespringParam`, sector preserving and inter-sector with allowed transitions (unreleased, 1.2.0) |
| CPTP, Kraus form | $`K_i\mapsto K_iM^{-1/2}`$ (identical to the polar retraction above) | same | inside `dynamics._kraus_step` only |
| unitary channels in $`\mathcal M`$ | $`U=U_0e^{G-G^*}`$ (cf. Lezcano-Casado and Martínez-Rubio 2019; Arjovsky et al. 2016) | `expm` | `optimize.UnitaryParam` |
| GKSL generators | $`H=(A+A^*)/2`$; Kossakowski matrix $`K=LL^*\succeq0`$ in a jump basis $`F_\alpha`$ adapted to the sectors (Gorini, Kossakowski and Sudarshan 1976; Lindblad 1976) | differentiable `lindblad_evolve` (`'kraus'`: every step CPTP; `'exact'`: $`e^{h\mathcal L}`$) | `optimize.LindbladParam`: free jumps as a Cholesky-type factor, in-sector and inter-sector (unreleased, 1.2.0) |
| states | $`\rho=(1-\varepsilon)GG^*/\operatorname{Tr}GG^*+\varepsilon/d`$, or Gibbs $`e^{-\beta H}/Z`$ | `states.gibbs_state` | `StateParam` in scripts only |
| unital CPTP | alternate the TP and unital normalisations (an operator version of Sinkhorn scaling) | — | proposal; convergence not analysed |

The Stinespring and Lindblad parametrisations have **gauge symmetries**. The environment unitary
$`V\mapsto(1\otimes W)V`$, $`W\in U(r)`$, does not change $`\Phi`$, and unitary mixing of jump operators does not
change $`\mathcal L`$. The polar parametrisation is also invariant under $`Z\mapsto cZ`$ for $`c>0`$. These
symmetries make the Hessian singular, which matters for natural-gradient methods. They also give conservation laws
(Section 5), which serve as cheap correctness checks of a training run.

**Missing pieces** (after the 1.2.0 work listed in Section 11, Stage 1):
- an `nn.Module` wrapper (Section 10);
- moving `StateParam` from `scripts/renyi_threshold.py` into the library, and rewriting the script on
  `optimize.StinespringParam`;
- the unital (Sinkhorn-type) parametrisation.

## 3. Expressivity and stability, phrased the network way

**Universality and depth.** One unrestricted layer is universal (Section 1). Continuous depth with a constant
generator is not:

- Wolf, Eisert, Cubitt and Cirac (2008) give computable criteria for whether a channel is of Markovian form.
- Wolf and Cirac (2008) show that some channels are indivisible, that is, they admit no non-trivial factorisation
  into channels.

A Lindblad-parametrised model therefore carries an inductive bias, Markovianity. The gap between the best Lindblad
fit and the best Stinespring fit to the same data is a quantitative witness of non-Markovianity. In the other
direction, Kliesch et al. (2011) show that the evolution of local, possibly time-dependent Liouvillians is approximated
by circuits of local channels whose size is polynomial in time and system size, with explicit Trotter error bounds.
This is a "deep composition approximates continuous depth" theorem with a known depth scaling. Section 6.2 sharpens
this for residual architectures: deep networks of near-identity blocks express exactly what Wolf and Cirac call the
infinitesimal divisible channels.

**Conditioning: Michelson contrast and rapidity.** For $`A>0`$ let $`\Delta(A)=(\lambda_{\max}-\lambda_{\min})/(\lambda_{\max}+\lambda_{\min})`$
(the library's `michelson_contrast`) and $`\eta(A)=\operatorname{artanh}\Delta(A)=\tfrac12\log\operatorname{cond}(A)`$.
The following statements are proved, either elementarily or in the cited sources.

1. *Products.* For invertible $`A,B`$, $`\operatorname{cond}(AB)\le\operatorname{cond}(A)\operatorname{cond}(B)`$.
   This follows from $`\sigma_{\max}(AB)\le\sigma_{\max}(A)\sigma_{\max}(B)`$ and
   $`\sigma_{\min}(AB)\ge\sigma_{\min}(A)\sigma_{\min}(B)`$. For $`A,B>0`$ it reads
   $`\eta(A^{1/2}BA^{1/2})\le\eta(A)+\eta(B)`$: rapidities add (companion paper on trace characterisations,
   proposition "Multiplicative composites"). Equivalently, this is the triangle inequality for Hilbert's projective
   metric $`d_H`$, with $`\eta(A)=\tfrac12d_H(A,1)`$.
2. *Positive maps.* Let $`\Phi`$ be positive with $`\Phi(1)`$ invertible, and let $`A>0`$. Then
   $`\lambda_{\min}(A)\Phi(1)\le\Phi(A)\le\lambda_{\max}(A)\Phi(1)`$, so
   $`\eta(\Phi(A))\le\eta(A)+\eta(\Phi(1))`$. Unital maps cannot increase the rapidity. By induction, for a deep
   composition
   ```math
   \eta\big(\Phi_L\circ\cdots\circ\Phi_1(A)\big)\ \le\ \eta(A)+\sum_{l=1}^{L}\eta\big(\Phi_l(1)\big).
   ```
   For TP layers acting on faithful states, $`\eta(\Phi_l(1))`$ measures how non-unital layer $`l`$ is.
3. *Why this matters for backpropagation.* Losses built from $`\log\rho`$ (Umegaki, Petz–Rényi) are differentiated by
   the Daleckii–Krein formula. The divided differences $`(\log a-\log b)/(a-b)`$ are bounded by
   $`1/\min(a,b)`$, so in Frobenius norm the derivative of $`\log`$ at $`\rho`$ has norm at most
   $`1/\lambda_{\min}(\rho)`$. For a state on $`\mathbb C^n`$ we have $`\lambda_{\max}\ge1/n`$, hence
   $`1/\lambda_{\min}\le n\,e^{2\eta(\rho)}`$. Together with (2), the regulariser
   $`\sum_l\eta(\Phi_l(1))`$ controls the gradient amplification of a log-based loss after $`L`$ layers. The bound is
   elementary. We have not found it used as a regulariser, and it must be written out with full proofs before any claim
   of novelty is made.
4. *Forgetting.* Positive maps contract $`d_H`$. If the image of the cone has finite projective diameter $`\delta`$, the
   contraction coefficient is $`\tanh(\delta/4)`$ (Birkhoff–Hopf; for quantum channels see Reeb, Kastoryano and Wolf
   2011). Deep compositions of strictly positive layers therefore forget their input geometrically. This is the
   projective-metric counterpart of noise-induced barren plateaus (Wang et al. 2021; McClean et al. 2018 for the
   noiseless case). The relation is an analogy, not a theorem about gradient variances.

## 4. Training problems as variational problems over channels

The first question for each problem is whether it is convex in the Choi matrix.

*Convex problems: semidefinite programming is the better tool.* Examples:

- least-squares or maximum-likelihood process tomography over *all* channels;
- optimal recovery for entanglement fidelity, which is linear in the Choi matrix of the recovery (Fletcher, Shor and
  Win 2007);
- diamond-norm and discrimination problems (Watrous 2018).

Restricting the channel to sectors is a linear constraint and keeps such problems convex. A solver such as CVXPY
(Diamond and Boyd 2016) returns the optimum together with a dual certificate. Gradient training is useful here only
as a fast batched approximation, or as the ground truth against which training dynamics are measured (Section 5).

*Non-convex problems: gradient training is the tool.* Examples:

- fixed Kraus rank, fixed depth or Markovian ($`e^{\mathcal L}`$) families;
- minimal output entropy, which minimises a concave function, so the optimum sits at extreme points (on additivity see
  Hastings 2009);
- maximal violation of data processing over $`(\rho,\sigma,\Phi)`$;
- control problems.

**Certification** follows `scripts/renyi_threshold.py` and `scripts/connes_stormer.py`. The final iterate is
re-orthonormalised in mpmath, so that it is an exactly TP channel, and the objective is recomputed at 50 significant
digits. A supremum then has a certified lower bound, because the point is feasible and its value is recomputed far
above the error of the computation. This is a high-precision recomputation, not interval arithmetic. Upper bounds need
duality or relaxations, not training.

**Target applications and first experiments.**

- **(a) Channels and Lindbladians from data, in sectors.** Process tomography as supervised learning, with
  $`(\rho_{\rm in},\rho_{\rm out})`$ pairs or Born probabilities as data. Related work: Ahmed, Quijandría and Kockum
  (2023) train Kraus operators on the Stiefel manifold; Torlai et al. (2023) use tensor networks; Lindblad tomography is
  done experimentally by Samach et al. (2022), and Bairey et al. (2020) learn Lindbladians from steady states.
  *First experiment:* synthetic data from a 6-qubit register with local damping, excitation and dephasing. Fit (i) a
  dense Stinespring model, (ii) a sector model restricted to $`N\to N,N\pm1`$, and (iii) a Lindblad model. Measure
  sample complexity against the parameter counts of Section 1, and use the Lindblad–Stinespring gap as the
  non-Markovianity witness of Section 3.
- **(b) Noise-aware control with post-selection.** This extends the paper's use case on noisy circuits with a conserved
  charge: an XY mixer, local emission, excitation and dephasing, and post-selection on the Hamming weight, with the
  verification idea of Bonet-Monroig et al. (2018) and McArdle, Yuan and Benjamin (2019). The idea is GRAPE (Khaneja
  et al. 2005) through the open dynamics and the Lüders conditioning. The parameters are piecewise-constant,
  sector-preserving controls. The objective is to maximise $`F_{\rm post}=F/p`$ subject to $`p\ge p_{\min}`$.
  `lindblad_evolve(method='kraus')` keeps every step a channel, and `condition_on` is differentiable.
  *First experiment:* $`n=6,\dots,10`$. Compare optimised and constant mixers, and check that the optimised schedules
  respect the structure the proposition in the paper predicts: post-selection removes pure loss to first order, and
  dephasing remains.
- **(c) Optimal recovery.** Compare three recoveries: Petz and universal recovery (`divergences.PetzRecovery`,
  `universal_recovery`; Petz 1986; Junge et al. 2018), the SDP optimum (Fletcher et al. 2007), and trained Stinespring
  recoveries of Kraus rank $`r=1,2,\dots`$. Barnum and Knill (2002) show the Petz-type recovery is near-optimal.
  *First experiment:* measure how small $`r`$ can be while staying within $`10^{-6}`$ of the SDP value.
- **(d) Extremal problems of operator algebra.** The data-processing thresholds of the Rényi divergences (sandwiched,
  $`\alpha\ge1/2`$: Frank and Lieb 2013, Müller-Lennert et al. 2013) were already rediscovered by certified training,
  in the paper's Section 6.1. Also in this class: minimal output entropy, and the Connes–Størmer variational problems
  of the companion paper. *First experiment:* reproduce the Rényi thresholds with the library classes of Section 10
  (`optimize.StinespringParam` in 1.2.0) instead of the script.

## 5. Training is itself a dynamics of the channel

A training run is a trajectory $`t\mapsto\Phi_t`$ in the manifold of channels. For the Stinespring parametrisation it
is a flow on a Stiefel manifold, modulo the $`U(r)`$ gauge, mapped to channels. Each item below is marked
**[E]** established, with references, **[P]** proposal, or **[N]** would be new as far as we know, to be checked.

**(a) First order: gradient flows and their metric.**

- **[E]** Gradient descent discretises $`\dot\theta=-\nabla\ell(\theta)`$. The induced flow on channels is
  $`\dot\Phi=-JJ^*\nabla_\Phi\ell`$ with $`J=D\Phi(\theta)`$. This depends on the parametrisation.
- **[E]** The natural gradient $`\dot\theta=-G^{+}\nabla\ell`$, with $`G=J^*gJ`$ the pull-back of a metric $`g`$, is
  parametrisation-invariant (Amari 1998). The quantum natural gradient uses the real part of the quantum geometric
  tensor (Fubini–Study) for parametrised pure states (Stokes et al. 2020). Monotone metrics on states are
  classified by operator monotone functions (Petz 1996). Bures (SLD), Kubo–Mori and the Wigner–Yanase metrics are
  examples, and all contract under channels.
- **[P]** Natural gradient for channels: use $`g`$ = a monotone metric on the normalised Choi state
  $`J(\Phi)/\dim\mathcal H_{\rm in}`$, or the embedded Stiefel metric on $`V`$. For divergence losses the Kubo–Mori
  metric is the Hessian of the Umegaki relative entropy, the analogue of Amari's Fisher metric for KL. Because of the
  gauge, $`G`$ is singular, so the flow uses a damped pseudo-inverse. This is affordable only for small problems,
  through `torch.func` Jacobians.
- **[E]** Symmetry gives conservation laws for gradient flows (Kunin et al. 2021, deep networks). In our setting the
  polar parametrisation is scale invariant, so $`\mathrm{Re}\langle\nabla\ell(Z),Z\rangle=0`$. Gradient flow then
  conserves $`\|Z\|_F^2`$ exactly. Plain gradient descent satisfies $`\|Z_{k+1}\|^2=\|Z_k\|^2+h^2\|\nabla\ell\|^2`$
  exactly, so the effective step shrinks. Adam breaks the identity. This gives an exact, free test of a training run.

**(b) Second order: least action.**

- **[E]** Momentum methods discretise damped second-order dynamics. The Euler–Lagrange equations of
  $`t^3\big(\tfrac12|\dot X|^2-f(X)\big)`$ give Nesterov's ODE $`\ddot X+\tfrac3t\dot X+\nabla f=0`$ (Su, Boyd and
  Candès 2016). Those of $`e^{\gamma t}\big(\tfrac12|\dot X|^2-f\big)`$ give the heavy ball
  $`\ddot X+\gamma\dot X+\nabla f=0`$. The Bregman Lagrangian of Wibisono, Wilson and Jordan (2016) generates the
  accelerated methods, including non-Euclidean ones, as Euler–Lagrange flows: optimisation by a least-action principle.
- **[P]** On channels: $`\mathcal L(V,\dot V,t)=e^{a(t)}\big[\tfrac12g_V(\dot V,\dot V)-e^{b(t)}\ell(V)\big]`$ on the
  Stiefel manifold (or the channel manifold of fixed rank). Its Euler–Lagrange flow is
  $`\nabla_t\dot V+\dot a\,\dot V+e^{b}\operatorname{grad}\ell=0`$: a geodesic flow with damping and a force. It is
  discretised by a retraction plus vector transport (Absil et al. 2008). Without the loss, the trajectories are
  geodesics of the chosen metric between channels.
- **[P, elementary]** Noether. Suppose a Lie group acts by isometries of $`g`$ and leaves $`\ell`$ invariant, with
  generator $`X`$ (skew-Hermitian). Then $`J_X=e^{a(t)}g(\dot V,X\cdot V)`$ is conserved along the Euler–Lagrange
  flow. Hence $`g(\dot V,XV)=J_X(0)\,e^{a(0)-a(t)}`$, which decays like $`e^{-\gamma t}`$ for the heavy ball and like
  $`t^{-3}`$ for Nesterov's ODE. Two groups are relevant:
  - the environment gauge $`U(r)`$, which is always present;
  - a physical symmetry $`G`$ acting by covariance, $`V\mapsto(U_g\otimes1)VU_g^*`$, whenever loss and data are
    $`G`$-invariant.
- **[E]** If loss and metric are $`G`$-invariant, the gradient is $`G`$-equivariant, and the set of covariant channels
  is invariant under the flow: training started on sector-structured channels stays there. By Palais' principle of
  symmetric criticality (1979, compact $`G`$), critical points of the restricted problem are critical points of the
  full one. This justifies training only the blocks. The global optimum of a non-convex problem may still break the
  symmetry. For convex problems, averaging an optimum over $`G`$ gives a covariant optimum.

**(c) The converse: dissipative physics is a gradient flow.**

- **[E]** The Fokker–Planck equation is the steepest descent of the free energy in the Wasserstein metric (Jordan,
  Kinderlehrer and Otto 1998). Carlen and Maas (2017) treat ergodic quantum Markov semigroups with detailed balance
  with respect to $`\sigma`$, on finite-dimensional algebras. For these, the evolution of states is the gradient flow
  of $`D(\cdot\|\sigma)`$ for a non-commutative 2-Wasserstein-type Riemannian metric, and they derive entropy decay
  inequalities from it. Mittnenzweig and Mielke (2017) show that Lindblad equations with detailed balance with respect
  to the thermal state are gradient systems driven by relative entropy, and they couple them to macroscopic systems in
  GENERIC form. Independently of detailed balance, $`D(\rho_t\|\sigma)`$ is non-increasing when $`\sigma`$ is
  stationary. The entropy production $`-\operatorname{Tr}\big(\mathcal L(\rho)(\log\rho-\log\sigma)\big)\ge0`$ is
  due to Spohn (1978).
- **[E]** $`D(\rho\|\sigma_\beta)=\beta\big(F(\rho)-F(\sigma_\beta)\big)`$ with $`F(\rho)=\operatorname{Tr}\rho H-\beta^{-1}S(\rho)`$.
  Imaginary-time evolution lowers the energy: for $`|\psi_\tau\rangle\propto e^{-\tau H}|\psi\rangle`$,
  $`\tfrac{d}{d\tau}\langle H\rangle=-2\operatorname{Var}_\psi H`$. McArdle et al. (2019) simulate it variationally
  through McLachlan's principle. Their update has the form of a natural-gradient step with a quantum-geometric metric,
  as in Stokes et al. (2020). We state this as an observation, not as a theorem of either paper.
- **[P]** "Learning = channel dynamics" and "dissipation = gradient flow of an entropy" meet:
  - A learning rule can be chosen to *be* a physical Markovian evolution. To learn $`\sigma`$, run a detailed-balance
    Lindbladian. Every `kraus` step is then a channel, and the flow is a gradient flow of $`D(\cdot\|\sigma)`$.
  - Conversely, a physical semigroup can be read as the training of its stationary state. The jump operators fix
    the metric.
  - The Hamiltonian part of a Lindbladian plays the role of the conservative part of the damped Hamiltonian flows in
    (b). This is an analogy to be made precise, not a theorem.

**(d) First experiments, possible with the current library.**

- **X1. Flows on a convex test problem.** Optimal recovery for a small code under amplitude damping, with the global
  optimum known from the SDP. Compare:
  1. gradient descent on $`Z`$;
  2. Riemannian gradient descent with polar retraction;
  3. natural gradient with the Bures metric on the Choi state;
  4. heavy ball;
  5. Nesterov ODE with retraction.

  Measure the gap to the SDP value against steps and time. Also measure three conservation laws:
  - $`\|Z\|^2`$, exact identity for gradient descent;
  - $`J_X`$ for the $`U(r)`$ gauge, predicted decay $`e^{-\gamma t}`$ or $`t^{-3}`$, started with non-zero velocity;
  - for a U(1)-covariant instance trained in the dense representation, the leakage out of the sector structure, which
    should stay at round-off.
- **X2. The non-convex counterpart.** The Rényi-threshold search of the paper's Section 6.1. Measure the number of
  starts and steps needed, under each flow, to find a certified violation near $`\alpha=1/2`$.
- **X3. Carlen–Maas in sectors.** Take $`n`$ sites, $`H=\omega\sum_iS^z_i`$, jumps $`S^\pm_i`$ with rates
  $`\gamma_\uparrow=\gamma_\downarrow e^{-\beta\omega}`$, so detailed balance holds with respect to the
  block-diagonal Gibbs state. The jumps move between U(1) sectors (`SpinChain.sector_jump`). Check:
  1. that $`D(\rho_t\|\sigma)`$ is monotone;
  2. that $`\tfrac{d}{dt}D`$ equals minus Spohn's production, with the derivative taken by autograd through the
     Daleckii–Krein backward of $`\log`$;
  3. that the production equals $`g_\rho(\operatorname{grad}D,\operatorname{grad}D)`$ for the Carlen–Maas metric,
     implemented for small $`n`$ from its definition;
  4. that natural-gradient training of a parametrised state toward $`\sigma`$ in that metric reproduces the Lindblad
     trajectory, up to discretisation error.

  As a control, add a drive that breaks detailed balance: (1) and (2) must survive, and (3) must fail.

## 6. Architectures and skip connections

Marking as in Section 0. Throughout, $`\|\cdot\|`$ on maps is the norm induced by the trace norm, in which every
channel has norm 1 and spectral radius 1.

### 6.1 Residual blocks

- **[E]** A residual block $`x_{l+1}=x_l+f_l(x_l)`$ (He et al. 2016) is a forward-Euler step of $`\dot x=f(x,t)`$. Deep
  residual networks are read as discretised ODEs and analysed for stability in that picture (E 2017; Haber and Ruthotto
  2017), and the limit of infinitely many small blocks is the neural ODE (Chen et al. 2018).
- **[P, elementary] The literal residual step leaves the channel set.** For a GKSL generator $`\mathcal L`$ the Euler
  step $`\mathrm{id}+h\mathcal L`$ is CP iff $`\mathcal L+h^{-1}\mathrm{id}`$ is CP. For a pure Hamiltonian generator
  $`\mathcal L=-i[H,\cdot\,]`$ on $`M_n`$ with $`H`$ not scalar this fails for every $`h>0`$. Write
  $`|\Omega\rangle=\sum_i|ii\rangle`$ and $`(H\otimes1)|\Omega\rangle=c|\Omega\rangle+|b\rangle`$ with
  $`b\perp\Omega`$. Then the Choi matrix is $`|\Omega\rangle\langle\Omega|-ih(|b\rangle\langle\Omega|-|\Omega\rangle\langle b|)`$.
  On $`\operatorname{span}\{\Omega,b\}`$ its determinant is $`-h^2\|\Omega\|^2\|b\|^2<0`$, so it has a negative eigenvalue
  $`-h^2\|H-\tfrac{\operatorname{tr}H}{n}1\|_2^2+O(h^3)`$. As for vectors in a simplex, $`x+f(x)`$ is admissible only
  when $`f`$ has a special form.
- **[E] The CP residual block is a mixture with the identity.** $`R_p=(1-p)\,\mathrm{id}+p\,\Phi`$, with $`\Phi`$ CPTP
  and $`0\le p\le1`$, is a channel, and $`R_p-\mathrm{id}=p(\Phi-\mathrm{id})`$. By Wolf and Cirac (2008, Lemma 1),
  $`e^{t(\Phi-\mathrm{id})}`$ is a CP semigroup for every channel $`\Phi`$. Its generator
  $`\Phi-\mathrm{id}=\sum_iK_i\cdot K_i^*-\tfrac12\{\sum_iK_i^*K_i,\cdot\}`$ is a Lindbladian whose jump operators are
  the Kraus operators of $`\Phi`$ and whose Hamiltonian is zero. Two readings follow, each elementary.
  - *Poisson jumps.* With $`p=\lambda h`$ and $`n=t/h`$, $`R_p^{\,n}\to e^{\lambda t(\Phi-\mathrm{id})}=\sum_{m\ge0}e^{-\lambda t}\tfrac{(\lambda t)^m}{m!}\Phi^m`$.
    One block is one time step of a Poisson process that applies $`\Phi`$ at its jumps. The limit is a network of
    Poisson-distributed depth.
  - *ResNet → Lindblad.* In the classical picture, residual networks lead to neural ODEs. For channels with
    mixture skips, they lead to Lindblad equations, and a stack of blocks $`R_{p_l}(\Phi_l)`$ is a time-dependent
    master equation with piecewise-constant generators $`(p_l/h)(\Phi_l-\mathrm{id})`$.
- **[P, elementary] Every Lindbladian is a limit of mixture generators.** Write
  $`\mathcal L(\rho)=\Psi(\rho)+K\rho+\rho K^*`$ with $`\Psi`$ CP and $`K=-iH-\tfrac12\Psi^*(1)`$. Then set
  $`\Phi_\kappa(\rho)=\kappa^{-1}\Psi(\rho)+(1+K/\kappa)\rho(1+K/\kappa)^*`$. It satisfies
  $`\Phi_\kappa^*(1)=1+K^*K/\kappa^2`$, and after the normalisation $`\tilde\Phi_\kappa=\Phi_\kappa(M^{-1/2}\cdot M^{-1/2})`$,
  $`M=\Phi_\kappa^*(1)`$, it is a channel with $`\kappa(\tilde\Phi_\kappa-\mathrm{id})=\mathcal L+O(1/\kappa)`$. The
  Hamiltonian part is reached only in the limit. A product of finitely many blocks $`R_{p}`$ with $`p<1`$ is a convex
  combination of words that contains $`\mathrm{id}`$ with positive weight. Such a product is therefore never a
  non-trivial unitary channel, because unitary channels are extreme points.
- **[E, library] The `kraus` step is a second-order residual block.** `dynamics._kraus_step` applies
  $`\Psi_h=S\circ(\mathrm{id}+h\Phi+\tfrac{h^2}2\Phi^2)\circ S`$ with $`S=E\cdot E^*`$, $`E=e^{-ih H_{\rm eff}/2}`$. It
  then normalises, $`\rho\mapsto\Psi_h(M^{-1/2}\rho M^{-1/2})`$. The jump part is a truncated exponential of
  $`h\Phi`$: the second-order version of the mixture block. The drift is applied by conjugation, and every factor is
  CP. This is the construction of the previous item at second order.
- **[P, elementary] Stochastic depth is the mixture skip on average.** Stochastic depth (Huang et al. 2016) keeps
  block $`l`$ with probability $`p_l`$ in training and replaces it by the identity otherwise. The average of
  $`x+b\,f(x)`$, $`b\sim\mathrm{Bernoulli}(p)`$, is $`(1-p)x+p\,(x+f(x))`$. For a channel block $`\Phi`$ applied with
  probability $`p`$, the average channel is exactly $`R_p`$. A stack of $`L`$ equal blocks gives the binomial mixture
  $`\sum_m\binom Lm p^m(1-p)^{L-m}\Phi^m=R_p^{\,L}`$, whose Poisson limit is the semigroup above.
- **Prior work and literature check.** Quantum residual networks exist. Wen et al. (2024) add an auxiliary qubit so
  that the unevolved state passes alongside the evolved one (a coherent combination). Khatun, Karim and Usman (2026,
  preprint) implement the residual connection as a deterministic mixture of the identity and a variational unitary:
  $`R_p`$ with unitary $`\Phi`$, motivated by avoiding post-selection. That $`\Phi-\mathrm{id}`$ generates a CP
  semigroup is Wolf and Cirac's Lemma 1. In two web searches (October 2026) we found neither the reading of
  mixture-skip networks as Lindblad and infinitesimal-divisible evolutions (Section 6.2) nor the stochastic-depth
  identity used for quantum channels. This is **[N?]**: not to be claimed before a full check.

### 6.2 Expressivity: residual channel networks and Markovianity

- **[E]** One block with $`p=1`$ is universal (Section 1). Expressivity is therefore a question only for *small*
  blocks, $`p\le\varepsilon`$ or more generally $`\|T_l-\mathrm{id}\|\le\varepsilon`$.
- **[E] Wolf and Cirac (2008).**
  - Some channels are *indivisible*: they admit no factorisation into non-unitary channels.
  - Definition 14: a channel is *infinitesimal divisible* if it lies in the closure $`\bar I`$ of the set $`I`$ of
    channels that, for every $`\varepsilon>0`$, are a finite product of channels $`T_i`$ with
    $`\|T_i-\mathrm{id}\|\le\varepsilon`$.
  - Proposition 15: infinitesimal divisible channels have $`\det T\ge0`$, the determinant of the superoperator.
  - Theorem 16: $`I=I'=J`$. Here $`I'`$ restricts the factors to Markovian channels $`e^{\mathcal L_i}`$, and $`J`$ is
    the set of elements of continuous completely positive evolutions. So every infinitesimal divisible channel is
    approximated by products of Markovian channels.
  - For qubits they give a complete characterisation and identify the solutions of time-dependent master equations
    with continuous generators. In higher dimensions they expect such a classification to be hard.
  - The neighbourhood of the identity contains channels of every type, including indivisible ones.
- **[E]** Wolf, Eisert, Cubitt and Cirac (2008) give a test for whether a channel is of the form $`e^{\mathcal L}`$.
  Cubitt, Eisert and Wolf (2012) show that deciding this is NP-hard in general.
- **[N?] Claim.** Let $`R_\varepsilon`$ be the set of finite products of mixture blocks $`(1-p_l)\mathrm{id}+p_l\Phi_l`$
  with $`p_l\le\varepsilon`$. Then
  ```math
  \bigcap_{\varepsilon>0}\overline{R_\varepsilon}\;=\;\bar I\quad\text{(Wolf and Cirac's infinitesimal divisible channels).}
  ```
  What is proved, and what is not:
  - "$`\supseteq`$" **[P, elementary]**. By Theorem 16 it suffices to approximate products of Markovian channels.
    Each $`e^{\mathcal L}`$ is a limit of $`R_p^{\,n}`$ with $`\Phi=\tilde\Phi_\kappa`$ from Section 6.1. Take any
    $`p\le\varepsilon`$ and $`n=\kappa/p`$ blocks. Since $`\|\tilde\Phi_\kappa-\mathrm{id}\|=O(1/\kappa)`$, the
    error of $`R_p^{\,n}`$ against $`e^{np(\tilde\Phi_\kappa-\mathrm{id})}`$ is $`O(p/\kappa)`$, and the generator
    error is $`O(1/\kappa)`$.
  - "$`\subseteq`$" needs Wolf and Cirac's error estimate in the proof of Theorem 16. Mixture blocks are
    $`2p`$-close to the identity, but our set takes closures for each $`\varepsilon`$, while their $`I`$ asks for
    exact products. To be checked.
  - The necessary condition is elementary and independent of that check. The eigenvalues of $`R_p`$ are
    $`1-p+p\lambda`$ with $`|\lambda|\le1`$. The real ones are $`\ge1-2p>0`$ for $`p<\tfrac12`$, and the complex ones
    come in conjugate pairs, so $`\det R_p>0`$. Every element of the left-hand side therefore has $`\det T\ge0`$.
  - *Literature check:* the equality for general near-identity factors is Wolf and Cirac's. The novelty, if any, is
    only the reduction to mixture blocks and the network reading. Two web searches found no source linking infinitesimal
    divisibility to residual or neural-ODE architectures.
- **[P] Consequences and a falsifiable prediction.**
  - Small-block residual networks carry the inductive bias "infinitesimal divisible", which is weaker than the
    semigroup bias $`e^{\mathcal L}`$ of a constant Lindblad model (Section 3). The three nested model classes,
    semigroup ⊂ small-block residual ⊂ unrestricted Stinespring, give two fit gaps. These gaps are witnesses of two
    different kinds of non-Markovianity.
  - Example: the Pauli channel $`\rho\mapsto\tfrac13(X\rho X+Y\rho Y+Z\rho Z)`$ has Pauli transfer matrix
    $`\operatorname{diag}(1,-\tfrac13,-\tfrac13,-\tfrac13)`$ and $`\det=-\tfrac1{27}<0`$. Its fit error by small-block
    residual networks must therefore stay above a positive constant, whatever the depth. This is checked in Stage 3
    (Section 11).

### 6.3 Copying skips, routing and coherent skips

- **[E] Copying is forbidden.** DenseNet (Huang et al. 2017) and U-Net (Ronneberger, Fischer and Brox 2015) concatenate
  *copies* of earlier activations. For quantum states the copy map does not exist: no-cloning (Wootters and Zurek 1982)
  for pure states, and no-broadcasting for mixed states (Barnum et al. 1996). A family of states can be broadcast iff
  its members commute.
- **[P, elementary] Fan-out through the centre.** A broadcasting channel for *all* states of $`\mathcal M`$ exists iff
  $`\mathcal M`$ is commutative, by Barnum et al. applied to a non-commuting pair. The commutative part that is always
  available is the centre. The sector measurement $`\rho\mapsto\sum_cP_c\rho P_c`$, with central $`P_c`$, is the identity
  on states of $`\mathcal M`$. Its outcome, the sector label, can be copied as often as needed. So an instrument
  records its outcome in the centre of the output algebra (Section 1), and every later layer may be *classically
  controlled* by it: $`\rho\mapsto\sum_j\Psi_j(\Phi_j(\rho))`$ (feed-forward) is again a channel. Classical skip
  connections are allowed. Quantum ones are not.
- **[P] Mixture of experts = block-diagonal channel.** In a mixture of experts (Jacobs et al. 1991) a gate routes the
  input to one of several experts. With the sector projections $`P_k`$ (central) as the router and a channel
  $`\Phi_k`$ per sector, the network is $`\rho\mapsto\bigoplus_k\Phi_k(\rho_k)`$. More generally, on a control register,
  it is $`\sum_k\mathrm{Ad}(P_k)\otimes\Phi_k`$ with $`\mathrm{Ad}(P)(x)=PxP`$. This is exactly a sector-preserving
  `Channel`. The router is non-disturbing on $`\mathcal M`$ because $`P_k`$ is central, and the weight sharing of
  Section 1 is the expert structure. A soft router is an instrument on a control system followed by feed-forward.
- **[E] Coherent skips.**
  - Superposition of two channels by coherent control is possible, but the output depends on more than the two
    channels: it depends on their implementation. Chiribella and Kristjánsson (2019) describe this for superpositions
    of trajectories, and Abbott et al. (2020) for coherent control of channels.
  - For unitary blocks, a coherent residual $`(1+U)/2`$ is a linear combination of unitaries. It can be implemented
    with an ancilla and post-selection (Childs and Wiebe 2012), with success probability
    $`\|(1+U)\psi\|^2/4`$.
  - The renormalised output $`\rho\mapsto(1+U)\rho(1+U)^*/\operatorname{Tr}(\cdot)`$ is nonlinear in $`\rho`$
    (Section 7, item 1). Wen et al. (2024) realise residual connections coherently with auxiliary qubits; Khatun,
    Karim and Usman (2026) describe that construction as relying on post-selection, which theirs avoids.
  - **[P]** In the library a coherent skip is an `InterSectorChannel` into an algebra with an extra control sector,
    followed by conditioning on a central projection.

## 7. Where nonlinearity comes from

- **[E, elementary] Deep channel networks are deep linear networks.** Every channel is linear in $`\rho`$, so a
  composition $`\Phi_L\circ\dots\circ\Phi_1`$ is a product of linear maps with CPTP constraints. For unconstrained deep
  linear networks, the learning dynamics is nonlinear in the weights and shows long plateaus followed by rapid
  transitions (Saxe, McClelland and Ganguli 2014). **[P]** Similar plateaus are to be expected for deep channel
  networks; the Stage 3 experiments (Section 11) record loss curves to test this.
- **[E] Linear in the data, nonlinear in the parameters: kernel methods.** A model output
  $`f_\theta(\rho)=\operatorname{Tr}(O\,\Phi_\theta(\rho))=\operatorname{Tr}(\Phi_\theta^*(O)\rho)`$ is a linear
  functional of $`\rho`$. With a data encoding $`x\mapsto\rho(x)`$ it is a linear model in the feature space of density
  matrices, i.e. a kernel method with kernel $`\operatorname{Tr}\rho(x)\rho(x')`$ (Schuld 2021). No channel network
  computes, for example, the purity $`\operatorname{Tr}\rho^2`$ from one copy.
- **Sources of genuine nonlinearity in $`\rho`$**, with their costs:
  1. **Post-selection.** $`\rho\mapsto P\rho P/\operatorname{Tr}(P\rho)`$.
     - **[E]** The cost is the acceptance probability: for a chain of post-selected layers the acceptances
       multiply, $`p=\prod_lp_l`$ (conditional acceptances), and the expected number of repetitions is $`1/p`$.
     - **[E]** The quantum neuron of Cao, Guerreschi and Aspuru-Guzik (2017) obtains a sigmoid-like activation from
       repeat-until-success circuits. Aaronson (2005) shows that post-selection raises the power of polynomial-time
       quantum computation to PP, which indicates that unrestricted post-selection cannot be cheap.
     - **[P]** The natural regulariser is $`-\log p=\sum_l(-\log p_l)`$, additive over layers.
     - **[P]** In the algebra, post-selection on a central projection is `condition_on` / `Instrument.post_states`,
       and it is differentiable. A *post-selection-activated network* interleaves channel layers with central
       conditioning (Stage 3).
  2. **Several copies.** $`\rho^{\otimes k}`$ makes the output a polynomial of degree $`k`$ in $`\rho`$, e.g.
     $`\operatorname{Tr}\rho^2=\operatorname{Tr}(F\rho\otimes\rho)`$ with the swap $`F`$. Two-copy nonlinear
     transformations of qubits are constructed explicitly by Bechmann-Pasquinucci, Huttner and Gisin (1998). The cost
     is $`k`$ copies per evaluation. **[P]** The $`S_k`$-symmetric part of $`\mathcal M^{\otimes k}`$ is the setting of
     the library's `permutational` module (Schur–Weyl duality).
  3. **Feedback at the trajectory level.** **[E]** Under continuous measurement the conditional state obeys a
     nonlinear stochastic master equation, and feedback can act on it. Averaged over the records, the evolution is
     linear again (Wiseman and Milburn 2009). The nonlinearity is per trajectory, as for post-selection without
     discarding.
  4. **Mean-field limits.** **[E]** For $`N`$ subsystems with permutation-symmetric, mean-field-scaled linear
     dynamics, the one-body state obeys in the limit $`N\to\infty`$ a nonlinear evolution: the nonlinear quantum
     dynamical semigroups of Alicki and Messer (1983). **[P]** The library's `PermutationalSystem` and
     `PermutationalLindbladian` evolve the exact $`N`$-body dynamics in the Dicke-basis algebra
     $`\bigoplus_jM_{2j+1}`$, which has $`\sum_j(2j+1)^2=(N+1)(N+2)(N+3)/6`$ complex dimensions. The emergence of the
     nonlinear one-body dynamics can then be measured as $`N`$ grows ("mean-field emergence", Stage 3).

  Nonlinear functions of output statistics in the *loss* (Section 1) do not make the network itself nonlinear.

## 8. Tomography as training

Tomography is supervised learning of a state or channel from measurement counts. The loss is a negative
log-likelihood, and the parametrisations of Section 2 make every iterate physical.

- **[E] Neural-network and tensor-network tomography.**
  - Torlai et al. (2018) reconstruct many-body states from measurement data with a neural-network (restricted
    Boltzmann machine) ansatz.
  - Torlai et al. (2023) do process tomography with tensor networks and unsupervised learning.
  - Ahmed, Quijandría and Kockum (2023) learn Kraus operators by gradient descent on the Stiefel manifold.
- **[E] Classical shadows.**
  - Huang, Kueng and Preskill (2020) predict many properties of a state from few randomised measurements.
  - Zhao and Miyake (2024) adjust classical shadows to how a symmetry, e.g. global U(1), is corrupted by the device
    noise: *symmetry-adjusted classical shadows*, an error-mitigation method with rigorous sampling bounds under
    readout errors.
- **[E] Permutationally invariant tomography** (Tóth et al. 2010) reconstructs the permutationally invariant part of
  the density operator with an optimised measurement strategy, scalable to large registers. It was demonstrated on a
  photonic four-qubit symmetric Dicke state.
- **[P, elementary] Tomography in an algebra.** Let $`\mathcal M\subset B(\mathcal H)`$ be a unital $`*`$-subalgebra
  and $`E_{\mathcal M}`$ the trace-preserving conditional expectation onto it. Then measuring only observables of
  $`\mathcal M`$ identifies exactly $`E_{\mathcal M}(\rho)`$.
  - *Proof sketch.* The map $`\rho\mapsto(\operatorname{Tr}\rho A)_{A\in\mathcal M}`$ factors through
    $`E_{\mathcal M}`$, because $`\operatorname{Tr}\rho A=\operatorname{Tr}E_{\mathcal M}(\rho)A`$ for $`A\in\mathcal M`$.
    It is injective on $`\mathcal M`$: if $`\sigma\in\mathcal M`$ and $`\operatorname{Tr}\sigma A=0`$ for all
    $`A\in\mathcal M`$, then $`A=\sigma^*`$ gives $`\operatorname{Tr}\sigma\sigma^*=0`$. Two states with the same
    $`E_{\mathcal M}(\rho)`$ cannot be told apart, and for $`\rho\in\mathcal M`$ (a symmetric state) the
    reconstruction is complete.
  - *Instances.* Permutationally invariant tomography is the case $`\mathcal M=\pi(S_N)'`$. Observables that commute
    with the excitation number identify the pinching onto the U(1) sectors.
  - *Parameter count.* The model has $`\dim\mathcal M-1=\sum_ck_c^2-1`$ real parameters, against $`4^n-1`$. For the
    permutationally invariant algebra this is $`(N+1)(N+2)(N+3)/6-1`$.
  - *Library.* Maximum likelihood over states of $`\mathcal M`$ is `measurement.POVM.log_likelihood` with a state
    parameter (Stage 2).
  - *Literature check:* the identification statement is elementary and implicit in Tóth et al. We state it as [P]
    only for the uniform treatment of all symmetry algebras, not as new.
- **[P] Identifiability by Fisher information.** `measurement.fisher_information` (autograd) gives the classical
  Fisher matrix of the counts with respect to the parameters. Its rank equals the number of parameters minus the gauge
  dimension (Section 2) minus the directions the POVM does not see. A rank deficit beyond the gauge diagnoses an
  informationally incomplete measurement within $`\mathcal M`$.
- **[E] Gate set tomography and its gauge** (Nielsen et al. 2021). GST characterises preparations, measurements and
  gates self-consistently, without pre-calibration. All probabilities are invariant under
  $`\rho\mapsto S\rho`$, $`E\mapsto ES^{-1}`$, $`G_i\mapsto SG_iS^{-1}`$, a gauge freedom that must be fixed.
  - **[P]** This gauge changes the reported channels. The Kraus/Stinespring gauge $`V\mapsto(W\otimes1)V`$ (Section 2)
    acts inside the representation of one channel and changes nothing physical.
  - For gate sets that all lie in a common algebra, the relevant GST gauge group should reduce to superoperators that
    preserve the sector structure. To be worked out before any claim.
- **[E] Lindbladian learning.**
  - Bairey et al. (2020) reconstruct local Markovian dynamics from steady states with local measurements.
  - Samach et al. (2022) do Lindblad tomography of a superconducting processor.
  - Wang and Li (2024) treat learning as inference with a master-equation solver that preserves complete positivity
    with guaranteed accuracy, which matters when the measurement times are far apart. The library's
    `lindblad_evolve(method='kraus')` serves the same purpose.
- **[E] Markovianity from a tomographic snapshot.** Wolf, Eisert, Cubitt and Cirac (2008) decide whether a channel is
  compatible with a semigroup $`e^{\mathcal L}`$; in general this is NP-hard (Cubitt, Eisert and Wolf 2012). Luchnikov
  et al. (2020) learn non-Markovian dynamics by maximum likelihood, embedding them in a Markovian dynamics of the system
  and a finite effective reservoir.
- **[P] Residual-network fit error as a non-Markovianity witness.** Fit the three nested classes of Section 6.2 to the
  same data: semigroup, small-block residual network, and unrestricted Stinespring. The first gap tests time-independent
  Markovianity and can be compared with the Wolf–Eisert–Cubitt–Cirac criterion on small cases. The second tests
  infinitesimal divisibility, with $`\det T<0`$ as an exact certificate. The witness is one-sided, because a large gap
  may be an optimisation failure. Multistart (`optimize.fit`) and the SDP reference for the convex unrestricted fit
  control this.
- **[E] Learning theory.**
  - Huang et al. (2022) prove and demonstrate exponential advantages of learners with quantum memory for predicting
    properties, for principal component analysis of noisy states, and for learning approximate models of dynamics.
  - For Pauli channels: Flammia and Wallman (2020) give efficient estimation, and Chen, Zhou, Seif and Jiang (2022)
    give quantum advantages with entangled probes.
  - Caro (2024) learns processes and Hamiltonians via the Pauli transfer matrix.
  - **[P]** These sample-complexity results are the yardstick for the sector-restricted models: the question is how
    much the parameter reduction of Section 1 reduces the number of samples in practice.
- **[P / N?] A thermometry witness from trace-inequality violations.** Theorem A of the companion paper on trace
  characterisations states the following. For $`W>0`$ and $`X`$ in $`M_n`$ with $`\Delta(W)=\delta`$ and
  $`\Delta(|X|)=\delta'`$, $`\delta\delta'<1`$,
  $`|\varphi_W(X)|\le(1-\delta^2\delta'^2)^{-1/2}\varphi_W(|X|)`$, where $`\varphi_W=\operatorname{Tr}(W\,\cdot)`$; the
  constant is attained for even $`n`$. Let $`v=|\varphi_W(X)|/\varphi_W(|X|)-1>0`$ be the measured relative violation.
  Inverting the bound gives
  ```math
  \Delta(W)\ \ge\ \frac{1}{\delta'}\sqrt{1-(1+v)^{-2}}\ \ge\ \sqrt{1-(1+v)^{-2}}.
  ```
  For a Gibbs state $`\Delta(W)=\tanh\big(\beta(E_{\max}-E_{\min})/2\big)`$, so with $`\delta'=1`$
  ```math
  \beta\,(E_{\max}-E_{\min})\ \ge\ 2\operatorname{arccosh}(1+v).
  ```
  This uses $`\operatorname{artanh}\sqrt{1-c^{-2}}=\operatorname{arccosh}c`$; for general $`\delta'`$ the bound is
  $`2\operatorname{artanh}\big(\delta'^{-1}\sqrt{1-(1+v)^{-2}}\big)`$.
  - *Interpretation.* The observed violation certifies a minimum contrast (non-traciality) of the state; for Gibbs
    states, a lower bound on $`\beta\,\Delta E`$. Because the constant is attained, the witness is tight for even
    $`n`$.
  - *Measurement and cost.* The data are estimates of $`\operatorname{Tr}(WX)`$ (two observables, the real and
    imaginary parts) and of $`\operatorname{Tr}(W|X|)`$. Finite statistics give a confidence interval for $`v`$, hence
    a certified bound at a stated confidence level.
  - *Literature check:* two web searches (October 2026) found no use of trace-inequality violations as temperature or
    contrast witnesses. The companion paper reports that earlier work samples the unnormalised violation but proves no
    bound. To be checked further (thermometry, KMS-state and modular-theory literature).

## 9. Encoders and bottlenecks

- **[P] Autoencoders between algebras.** An autoencoder is an encoder channel $`\mathcal E:\mathcal M\to\mathcal L`$
  into a smaller *latent algebra* $`\mathcal L=\bigoplus_cM_{k_c}`$ and a decoder $`\mathcal D:\mathcal L\to\mathcal M`$,
  trained on the reconstruction loss of $`\mathcal D\circ\mathcal E`$ over the data. The type of $`\mathcal L`$ is a
  design choice:
  - the centre $`\mathbb C^C`$ holds *classical latents*, the sector labels;
  - the blocks $`M_{k_c}`$ hold *quantum latents*;
  - $`\mathcal L=\mathbb C^K`$ gives a measure-and-prepare (entanglement-breaking) bottleneck, $`\mathcal L=M_k`$ a
    fully quantum one, and the general case is *hybrid*.
  `optimize.StinespringParam` with `alg_out` already parametrises both $`\mathcal E`$ and $`\mathcal D`$.
- **[E] Compression.** Schumacher (1995) compresses a quantum source to its typical subspace. The quantum
  autoencoder of Romero, Olson and Aspuru-Guzik (2017) trains a unitary encoder so that a "trash" register can be
  reset, compressing a family of states into fewer qubits.
- **[E] Error correction is encoder–noise–decoder.**
  - The Knill–Laflamme conditions (1997) characterise when a code space can be recovered after the noise.
  - Operator algebra quantum error correction (Bény, Kempf and Kribs 2007, PRL and PRA) corrects *algebras*
    (observables) instead of subspaces. The correctable objects are algebras of the form $`\bigoplus_cM_{k_c}`$ (with
    multiplicities), the latent algebras above.
  - Hybrid memories: Kuperberg (2003) shows that $`N`$ copies of a hybrid memory $`A`$ embed into $`N(1+o(1))`$
    copies of $`B`$ iff $`\|\mathrm{shape}(A)\|_p\le\|\mathrm{shape}(B)\|_p`$ for all $`p\ge1`$. Here
    $`\mathrm{shape}=(k_c)_c`$.
  - Grassl, Lu and Zeng (2017) construct hybrid codes that transmit qudits and classical symbols together, with better
    parameters than codes built from stabiliser codes.
- **[E] Lossless autoencoding is sufficiency.** By Petz (1986, 1988), $`\mathcal E`$ can be inverted on a family
  $`\{\rho_\theta\}`$ by a channel iff it is sufficient for it. This holds iff equality holds in data processing for the
  relative entropy, and then the Petz map is a universal recovery. In the approximate case, Junge et al. (2018) bound
  the recovery error by the loss of relative entropy. The library has both (`divergences.PetzRecovery`,
  `universal_recovery`). They serve as reference decoders against which trained decoders are measured (Stage 3).
- **[P, elementary] The conditional expectation is the linear autoencoder (the analogue of PCA).** The conditional
  expectation $`E_{\mathcal N}`$ onto $`\mathcal N\subset\mathcal M`$ that preserves a faithful trace $`\tau`$ is the
  orthogonal projection onto $`\mathcal N`$ for the inner product $`\tau(x^*y)`$. It is unital and
  $`\tau`$-preserving. Linear autoencoders with squared loss learn the principal subspace (Baldi and Hornik 1989); here
  the code space is constrained to be a $`*`$-subalgebra. Two loss measures are available in the library:
  - The *Pimsner–Popa constant* $`\lambda(E)=\max\{l:E(y)\ge l\,y\ \forall y\ge0\}`$ (Pimsner and Popa 1986) is the
    worst-case retained fraction of any positive input; it is $`1/n`$ for $`\mathbb C\subset M_n`$
    (`ConditionalExpectation.pimsner_popa_constant`). With the index $`[\mathcal M:\mathcal N]`$ (`Inclusion.index`),
    it measures how much a subalgebra bottleneck discards.
  - *Link to Section 6.* Let $`\lambda_{\rm cb}`$ be the constant of the complete order, i.e. the largest $`l`$ with
    $`E-l\,\mathrm{id}`$ CP (`pimsner_popa_constant(complete=True)`). Then
    $`E=\lambda_{\rm cb}\,\mathrm{id}+(1-\lambda_{\rm cb})\Psi`$ with $`\Psi`$ CP, unital and trace preserving. A
    conditional expectation is a mixture skip block with identity weight $`\lambda_{\rm cb}`$.
- **[E] Data encoders (feature maps).**
  - Lloyd et al. (2020) train the embedding of classical data into quantum states.
  - Schuld (2021) shows that the encoding determines the model class, through the kernel.
  - Pérez-Salinas et al. (2020) obtain a universal classifier by re-uploading the data several times.
  - **[P]** In an algebra, a feature map $`x\mapsto\rho(x)\in S(\mathcal M)`$ produces block-sum kernels
    $`\sum_c\operatorname{tr}\rho_c(x)\rho_c(x')`$, so the sector structure restricts which kernels are reachable.
- **[N?] The type of the bottleneck.**
  - *Proposal.* At fixed algebra dimension $`D=\sum_ck_c^2`$, compare classical $`\mathbb C^D`$, quantum
    $`M_{\sqrt D}`$ and hybrid $`\bigoplus_cM_{k_c}`$ latents on the same data. Data: noisy sector-structured states
    of the XY register. Decoders: trained, Petz, and the SDP optimum for a fixed encoder, which is convex in the
    decoder's Choi matrix.
  - *Why $`D`$ is fixed.* At fixed Hilbert dimension $`\sum_ck_c`$ the comparison is trivial, because
    $`\bigoplus_cM_{k_c}\subset M_{\sum k_c}`$ and the quantum latent wins.
  - *Literature check:* Kuperberg (2003) settles the asymptotic, data-blind comparison of hybrid memories by $`p`$-norms
    of the shape. Hybrid codes (Grassl, Lu and Zeng 2017) fix code parameters. Hybrid classical–quantum autoencoders in
    quantum machine learning, e.g. Sakhnenko et al. (2022), insert a quantum circuit into the bottleneck of a classical
    autoencoder and measure it back to a classical latent. We found no single-shot,
    data-dependent comparison of latent-algebra types at fixed $`\sum_ck_c^2`$, after two web searches; incomplete.

## 10. Proposed API: `torch_vn_algebra.learn`

The module is named `learn` to avoid confusion with `torch.nn`. All classes are `torch.nn.Module`s whose parameters
are unconstrained tensors. The constraint is enforced in `forward`, by retraction or normalisation, so any
`torch.optim` optimiser applies. Operators keep the library layout `(B, C, k_max, k_max)`, and channels are returned
as `Channel` / `InterSectorChannel` objects.

*Placement (Section 11).* The parametrisations come first, as plain parameter objects in `optimize`:
`StinespringParam`, `LindbladParam` and `fit` in 1.2.0. POVMs, instruments and the Fisher information go into a new
`measurement` module, also in 1.2.0. The `learn` classes below wrap these as `nn.Module`s: `ChannelSequential` and
the skips of Section 6 in 1.4 (Stage 3), and `ChannelFlow` and `learn.certify` in Stage 4. `learn.Instrument` is a
layer view of `measurement.Instrument`.

1. **`StinespringChannel(alg_in, alg_out, kraus_rank, pairs=None, retraction='polar'|'qr'|'cayley', trace_preserving=True)`**.
   Parameters: one complex tensor `Z` of shape `(B, r, D, C, m_max, k_max)`, masked to `pairs`, the allowed
   `(d, c)`. Methods: `forward(rho) -> Operator`, `channel() -> InterSectorChannel` (TP by construction),
   `num_parameters()` (the exact counts of Section 1), `gauge_momentum(velocity)`.
2. **`LindbladGenerator(alg, H=None, jump_basis=..., kossakowski='cholesky', controls=None, n_slices=1)`**. Trainable
   $`H`$, Cholesky factor $`L`$ of $`K=LL^*`$, and piecewise-constant control amplitudes. Methods: `evolve(rho, times,
   method='kraus')`, `channel(t)` (exact, small blocks), `trotter(n) -> ChannelSequential`.
3. **`ChannelSequential(*layers)`**. Composition with algebra checks. Layers may be the two classes above, fixed
   `Channel`s, `ConditionalExpectation`s (pooling) or `Instrument`s.
4. **`Instrument(alg_in, alg_out, n_outcomes, ...)`** (a layer view of `measurement.Instrument`), with the
   constructor `Instrument.postselect(alg, central_projection)`.
   Stored as an `InterSectorChannel` into $`\bigoplus_j\mathcal N`$. `forward` returns `(states, probabilities)`;
   post-selection renormalises, differentiably.
5. **`ChannelFlow(module, kind='gradient'|'natural'|'heavy_ball'|'nesterov_ode', metric='embedded'|'bures'|'kubo_mori', step)`**.
   The training dynamics of Section 5 as an optimiser with retraction and vector transport. It logs energy, momentum
   maps and scale invariants.
6. **`ResidualChannel(block, p, kind='mixture'|'coherent')`** and **`CentreFeedForward(instrument, branches)`**
   (Section 6). `'mixture'` is $`(1-p)\mathrm{id}+p\,\Phi`$, CPTP for every $`p\in[0,1]`$. `'coherent'` adds a
   control sector and returns `(state, acceptance)` after conditioning. `CentreFeedForward` applies
   `branches[j]` to the post-measurement state of outcome `j`. Both compose in `ChannelSequential`.

Losses are plain functions in `learn.losses`: infidelity, trace distance, Umegaki and Rényi divergences (wrapping
`divergences`), entanglement and process fidelity from the Choi matrix, and the rapidity regulariser
$`\sum_l\eta(\Phi_l(1))`$. `learn.certify(channel, objective, digits=50)` generalises the mpmath recomputation of the
scripts.

```python
# NOT IMPLEMENTED: proposed API (Stage 2, 1.3; the same fit is possible in 1.2.0 with optimize.StinespringParam).
# Process tomography in U(1) sectors.
import torch
from math import comb
from torch_vn_algebra import TypeIAlgebra, learn

ks = [comb(6, N) for N in range(7)]                      # 6 qubits, sectors N = 0..6
M = TypeIAlgebra(ks, ks, precision='double')
pairs = [(c, c) for c in range(7)] + [(c - 1, c) for c in range(1, 7)]   # stay, or lose one excitation
model = learn.StinespringChannel(M, M, kraus_rank=4, pairs=pairs, retraction='polar')
opt = torch.optim.Adam(model.parameters(), lr=1e-2)
for rho_in, rho_out in loader:                           # DensityMatrix batches of M
    loss = learn.losses.infidelity(model(rho_in), rho_out).mean()
    opt.zero_grad(); loss.backward(); opt.step()
print(model.num_parameters(), model.channel().trace_preservation_error().max())
```

```python
# NOT IMPLEMENTED: proposed API (Stage 4). Noise-aware control with post-selection on the weight w.
gen = learn.LindbladGenerator(M, H=H_xy, controls=[H_1, H_2], n_slices=40,
                              jumps=noise_jumps, rates=noise_rates, trainable=('controls',))
post = learn.Instrument.postselect(M, M.central(onehot(w)))
opt = torch.optim.LBFGS(gen.parameters())
def closure():
    opt.zero_grad()
    rho_T = gen.evolve(rho0, times=[0.0, T], method='kraus')[-1]   # every step a channel
    rho_post, p = post(rho_T)
    loss = -learn.losses.fidelity(rho_post, psi_T) + lam * torch.relu(p_min - p)
    loss.mean().backward(); return loss.mean()
opt.step(closure)
```

```python
# NOT IMPLEMENTED: proposed API (Stage 4). Training dynamics with conservation checks (experiment X1).
model = learn.StinespringChannel(N_alg, M_alg, kraus_rank=2)
for kind in ('gradient', 'natural', 'heavy_ball', 'nesterov_ode'):
    flow = learn.ChannelFlow(model.reset(seed=0), kind=kind, metric='bures', step=1e-2)
    for k in range(2000):
        loss = -entanglement_fidelity(model.channel(), noise, code)
        flow.step(loss)
        log(kind, k, loss.item(), sdp_value - (-loss.item()),
            model.Z.norm()**2, flow.momentum_map('environment'))
```

## 11. Plan, risks, measures of success

The plan replaces the version table of earlier drafts. Stages 1–3 build the library parts in the order the
applications need them. Stage 4 keeps the experiments on training dynamics, and Stage 5 runs alongside as separate
notes. The status column records the state of the branch at the time of writing.

| stage | content | status | done when |
|---|---|---|---|
| **0** (now) | finish the two preprints in preparation; release **1.2.0** = `interop`, `SpinChain.bond`, the fix of the `kraus`-method gradients at degenerate spectra (all three in the changelog under "Unreleased"), **and Stage 1** | changelog entries written; Stage 1 partly committed | both preprints on arXiv; 1.2.0 tagged with every Stage 1 criterion met and the full test suite passing |
| **1** (1.2.0) parametrised channels and measurement | `optimize.StinespringParam`: sector-preserving `Channel` and `InterSectorChannel` with allowed transitions, polar retraction, Haar initialisation. `optimize.fit`: multistart minimisation that returns the final point. `measurement` module: `POVM`, sampled counts, log-likelihood, `Instrument` into an algebra with the outcome in the centre, Fisher information by autograd. Optional SDP reference via cvxpy (extra `sdp`). `optimize.LindbladParam`: Kossakowski matrix via a Cholesky-type factor, sector covariant | `StinespringParam`, `LindbladParam`, `fit` committed; `measurement` in progress; SDP extra not started | TP error $`\le10^{-12}`$ (double); `gradcheck` passes, including at degenerate spectra; parameter counts match the formulas of Section 1; Gaussian initialisation reproduces the distribution of `random_channel` (two-sample test); `fit` recovers a known channel from noise-free data to $`10^{-8}`$; on convex instances `fit` and the SDP reference agree to $`10^{-6}`$; Fisher information matches finite differences to $`10^{-6}`$; `renyi_threshold.py` rewritten on the library reproduces its certified table |
| **2** (1.3) tomography in the algebra (Section 8) | state MLE of $`E_{\mathcal M}(\rho)`$; reproduce Tóth et al. (2010) with the `permutational` module; U(1) sectors; covariant process tomography against the SDP reference (application (a) of Section 4); learning the noise rates of the XY register, including post-selection statistics; identifiability by the Fisher information | not started | MLE error $`\times`$ shots within 10 % of $`\operatorname{tr}F^{-1}`$ (Cramér–Rao) at large shot numbers; the Tóth et al. four-qubit Dicke example reproduced on simulated counts, and $`N\le20`$ run in the Dicke-basis algebra; covariant process tomography agrees with the SDP reference to $`10^{-6}`$ on noise-free data; sample complexity of sector models against the dense Stinespring baseline in the style of Ahmed et al.; XY noise rates recovered within three standard errors; Fisher rank = parameters − gauge dimension |
| **3** (1.4) architectures and encoders (Sections 6, 7, 9) | `ChannelSequential`, residual (mixture) and coherent skips, feed-forward through the centre; ResNet ↔ Markovianity experiment against the Wolf–Eisert–Cubitt–Cirac criterion; post-selection-activated networks; mean-field emergence; encoders: learned encoder–decoder pairs into sector algebras against Petz recovery and the SDP optimum (application (c) of Section 4); bottleneck-type experiment | not started | algebra mismatches rejected by `ChannelSequential`; residual fits of Markovian targets converge with depth, and on the $`\det<0`$ Pauli channel of Section 6.2 the error stays above a positive constant; agreement with the Wolf–Eisert–Cubitt–Cirac test on random qubit channels away from the boundary; a post-selection-activated network solves a task no single channel can (e.g. purity estimation), with the acceptance cost reported; one-body error against the mean-field equation measured as a function of $`N`$; trained decoders within $`10^{-6}`$ of the SDP value, with the smallest sufficient Kraus rank $`r`$ reported (application (c)) |
| **4** (1.4–2.0) training as a dynamics (Section 5) | `ChannelFlow`, `learn.certify`, the `nn.Module` layer of Section 10; experiments X1–X3; control with post-selection (application (b) of Section 4); a linear-map type for rectangular blocks $`\operatorname{Hom}_G=\bigoplus_\lambda M_{n_\lambda\times m_\lambda}`$ between algebras (not CP) | not started | X1: every flow reaches the SDP optimum to $`10^{-8}`$, and the predicted conservation and decay laws hold to discretisation order; X3: the Carlen–Maas identity checked to $`10^{-10}`$; (b): post-selected infidelity of optimised against constant schedules at $`n\le10`$ |
| **5** (parallel notes) | thermometry witness from trace-inequality violations (Section 8); contrast decay under dissipative dynamics via the Birkhoff–Hopf contraction (Section 3, item 4; Reeb, Kastoryano and Wolf 2011) | not started | witness: proof written out, equality at the even-$`n`$ extremisers of the companion paper checked to $`10^{-12}`$, no violation of the inverted bound in $`10^5`$ random Gibbs states; contrast decay: measured contraction of $`\eta`$ along `lindblad_evolve` compared with $`\tanh(\delta/4)`$ on spin chains |

**Cross-cutting.**

- A bucketed (unpadded) layout, which groups blocks by size, before any experiment with $`n\ge8`$ (see the padding
  risk below).
- A tutorial chapter "Learning channels" in `man/tutorial/`, built on the Stage 1–2 code.
- A literature check for every [N] item before it enters a preprint, recorded next to the item.

**Candidate preprints.**

1. *Tomography in symmetry algebras*: Stages 1–2 and Section 8. Maximum likelihood of $`E_{\mathcal M}(\rho)`$,
   permutational and U(1) instances, and identifiability.
2. *Residual channel networks and Markovianity*: Stage 3 and Section 6. Mixture skips, infinitesimal divisibility,
   and the fit-gap witnesses.
3. *A thermometry witness from trace-inequality violations*: Stage 5 and Section 8.

**Risks.**

- *Degenerate spectra.* First derivatives of $`f(\rho)`$ are finite (Daleckii–Krein), but derivatives through
  eigenvectors are not, and second derivatives are correct only at non-degenerate spectra. Natural-gradient metrics
  are affected.
- *Conditioning of $`M^{-1/2}`$ and $`\log\rho`$.* Monitor and regularise with the rapidity bound of Section 3. The
  `kraus` step already halves $`h`$ when $`M`$ is ill-conditioned.
- *Padding cost.* The inter-sector tensor costs $`r\,D\,C\,m_{\max}k_{\max}`$, against the useful
  $`r\sum_{(d,c)}m_dk_c`$. For chains with local jumps the waste is about $`C/3`$ before padding. Mitigations: sparse
  pair lists, grouping by block size (the bucketed layout above), `SectorJump`.
- *Scalability.* Choi matrices and superoperators are $`k^2\times k^2`$ per block, so blocks are limited to a few
  tens. The library does not compete with escnn (Cesa et al. 2022) or e3nn (Geiger and Smidt 2022) on large equivariant
  feature spaces, or with circuit simulators on many qubits.
- *Optimisation landscape.* Non-convexity, gauge-singular Hessians, and contraction with depth (Section 3, item 4).
  For convex problems, an SDP should be preferred.
- *Post-selection cost.* Acceptance probabilities multiply (Section 7). Every post-selected result is reported with
  its acceptance probability.
- *One-sided witnesses.* Fit gaps (Sections 6.2 and 8) can come from optimisation failure. They are reported with
  multistart statistics and, where the problem is convex, with the SDP value.
- *Over-claiming.* Most individual ingredients exist (Section 12). Claims of novelty need a literature check per item.

## 12. Positioning

- **Quantum neural networks.** Beer et al. (2020) build deep networks whose layers are CP maps: a unitary on input and
  output qudits, followed by a partial trace over the input. This is the Stinespring layer of Section 1, trained with
  fidelity losses. Quantum convolutional networks (Cong, Choi and Lukin 2019) pool by measurement and conditioned
  unitaries, which is an instrument followed by feed-forward. Equivariant QNNs (Nguyen et al. 2024) use the same
  commutant ideas as the library's `decompose`. Quantum residual networks exist with coherent (Wen et al. 2024) and
  mixture (Khatun, Karim and Usman 2026) skips (Section 6.1). The post-selection quantum neuron is due to Cao,
  Guerreschi and Aspuru-Guzik (2017).
- **Libraries.** PennyLane (Bergholm et al.) and TensorFlow Quantum (Broughton et al.) differentiate gate circuits,
  with noise as fixed channels. dynamiqs (JAX) offers batched, differentiable Lindblad simulation on GPUs.
- **Channel and Lindbladian learning.** See Ahmed et al. (2023), Torlai et al. (2023), Samach et al. (2022),
  Bairey et al. (2020), Wang and Li (2024), and for non-Markovian dynamics Luchnikov et al. (2020). Learning theory:
  Huang et al. (2022), Flammia and Wallman (2020), Chen et al. (2022), Caro (2024). Section 8 places these.
- **Quantum autoencoders and codes.** Romero, Olson and Aspuru-Guzik (2017); operator algebra error correction (Bény,
  Kempf and Kribs 2007); hybrid memories and codes (Kuperberg 2003; Grassl, Lu and Zeng 2017). Section 9 places these.
- **Neural-network states.** Carleo and Troyer (2017), and for open systems Nagy and Savona (2019) and Vicentini et
  al. (2019). This is a different meaning of the phrase: there, networks are ansätze for states.
- **C\*-algebra nets.** Hashimoto, Wang and Matsui (2022) and Hashimoto, Ikeda and Kadri (2024) put network
  parameters in a C\*-algebra: the dual direction.
- **SPD networks.** In Huang and Van Gool (2017), the BiMap layer $`X\mapsto WXW^{\mathsf T}`$ is a Kraus-rank-one CP
  map, and ReEig/LogEig are functional calculus.

**What this direction adds, honestly stated.** The ingredients above are largely known. The contribution is
integration with the finite-dimensional von Neumann algebra as the native type:

- channels between *different* sector structures;
- the centre as outcome register and post-selection;
- traces and conditional expectations as pooling;
- exact parameter counts from the sector structure;
- CPTP-by-construction steps;
- gradients that stay finite at degenerate spectra;
- batched multistart on GPUs;
- certified evaluation.

Sections 6–9 add a dictionary for architectures (mixture skips as Lindblad steps, the centre as the only fan-out,
mixtures of experts as sector-preserving channels), tomography restricted to an algebra, and bottlenecks of chosen
algebra type. Section 5 adds a testbed for the dynamics of training on channel manifolds, with falsifiable predictions: conservation
and decay laws, invariance of covariant sets, and the Carlen–Maas identity. Items marked [N] or "not found" must be
checked against the literature before any publication claims them.

## References

- Aaronson, S.: Quantum computing, postselection, and probabilistic polynomial-time. Proc. R. Soc. A 461(2063), 3473–3482 (2005). doi:10.1098/rspa.2005.1546; arXiv:quant-ph/0412187.
- Abbott, A. A., Wechs, J., Horsman, D., Mhalla, M., Branciard, C.: Communication through coherent control of quantum channels. Quantum 4, 333 (2020). doi:10.22331/q-2020-09-24-333; arXiv:1810.09826.
- Absil, P.-A., Mahony, R., Sepulchre, R.: *Optimization Algorithms on Matrix Manifolds*. Princeton University Press (2008). ISBN 978-0-691-13298-3, doi:10.1515/9781400830244.
- Ahmed, S., Quijandría, F., Kockum, A. F.: Gradient-descent quantum process tomography by learning Kraus operators. Phys. Rev. Lett. 130, 150402 (2023). doi:10.1103/PhysRevLett.130.150402; arXiv:2208.00812.
- Alicki, R., Messer, J.: Nonlinear quantum dynamical semigroups for many-body open systems. J. Stat. Phys. 32(2), 299–312 (1983). doi:10.1007/BF01012712.
- Amari, S.: Natural gradient works efficiently in learning. Neural Comput. 10(2), 251–276 (1998). doi:10.1162/089976698300017746.
- Arjovsky, M., Shah, A., Bengio, Y.: Unitary evolution recurrent neural networks. Proc. ICML 2016, PMLR 48, 1120–1128 (2016).
- Bairey, E., Guo, C., Poletti, D., Lindner, N. H., Arad, I.: Learning the dynamics of open quantum systems from their steady states. New J. Phys. 22, 032001 (2020). doi:10.1088/1367-2630/ab73cd.
- Baldi, P., Hornik, K.: Neural networks and principal component analysis: learning from examples without local minima. Neural Netw. 2(1), 53–58 (1989). doi:10.1016/0893-6080(89)90014-2.
- Barnum, H., Caves, C. M., Fuchs, C. A., Jozsa, R., Schumacher, B.: Noncommuting mixed states cannot be broadcast. Phys. Rev. Lett. 76, 2818–2821 (1996). doi:10.1103/PhysRevLett.76.2818; arXiv:quant-ph/9511010.
- Barnum, H., Knill, E.: Reversing quantum dynamics with near-optimal quantum and classical fidelity. J. Math. Phys. 43(5), 2097–2106 (2002). doi:10.1063/1.1459754; arXiv:quant-ph/0004088.
- Bechmann-Pasquinucci, H., Huttner, B., Gisin, N.: Non-linear quantum state transformation of spin-1/2. Phys. Lett. A 242(4–5), 198–204 (1998). doi:10.1016/S0375-9601(98)00189-3; arXiv:quant-ph/9708040.
- Beer, K., Bondarenko, D., Farrelly, T., Osborne, T. J., Salzmann, R., Scheiermann, D., Wolf, R.: Training deep quantum neural networks. Nat. Commun. 11, 808 (2020). doi:10.1038/s41467-020-14454-2.
- Bény, C., Kempf, A., Kribs, D. W.: Generalization of quantum error correction via the Heisenberg picture. Phys. Rev. Lett. 98, 100502 (2007). doi:10.1103/PhysRevLett.98.100502; arXiv:quant-ph/0608071.
- Bény, C., Kempf, A., Kribs, D. W.: Quantum error correction of observables. Phys. Rev. A 76, 042303 (2007). doi:10.1103/PhysRevA.76.042303; arXiv:0705.1574.
- Bergholm, V., et al.: PennyLane: Automatic differentiation of hybrid quantum-classical computations. arXiv:1811.04968 (2018, rev. 2022).
- Bonet-Monroig, X., Sagastizabal, R., Singh, M., O'Brien, T. E.: Low-cost error mitigation by symmetry verification. Phys. Rev. A 98, 062339 (2018).
- Broughton, M., et al.: TensorFlow Quantum: A software framework for quantum machine learning. arXiv:2003.02989 (2020, rev. 2021).
- Buča, B., Prosen, T.: A note on symmetry reductions of the Lindblad equation: transport in constrained open spin chains. New J. Phys. 14, 073007 (2012).
- Cao, Y., Guerreschi, G. G., Aspuru-Guzik, A.: Quantum neuron: an elementary building block for machine learning on quantum computers. arXiv:1711.11240 (2017), preprint.
- Carlen, E. A., Maas, J.: Gradient flow and entropy inequalities for quantum Markov semigroups with detailed balance. J. Funct. Anal. 273(5), 1810–1869 (2017). doi:10.1016/j.jfa.2017.05.003.
- Carleo, G., Troyer, M.: Solving the quantum many-body problem with artificial neural networks. Science 355(6325), 602–606 (2017). doi:10.1126/science.aag2302.
- Caro, M. C.: Learning quantum processes and Hamiltonians via the Pauli transfer matrix. ACM Trans. Quantum Comput. 5(2), Article 14, 1–53 (2024). doi:10.1145/3670418; arXiv:2212.04471.
- Cesa, G., Lang, L., Weiler, M.: A program to build E(N)-equivariant steerable CNNs. ICLR 2022 (library escnn).
- Chen, R. T. Q., Rubanova, Y., Bettencourt, J., Duvenaud, D.: Neural ordinary differential equations. Advances in Neural Information Processing Systems 31 (NeurIPS 2018). arXiv:1806.07366.
- Chen, S., Zhou, S., Seif, A., Jiang, L.: Quantum advantages for Pauli channel estimation. Phys. Rev. A 105, 032435 (2022). doi:10.1103/PhysRevA.105.032435; arXiv:2108.08488.
- Childs, A. M., Wiebe, N.: Hamiltonian simulation using linear combinations of unitary operations. Quantum Inf. Comput. 12(11–12), 901–924 (2012). arXiv:1202.5822.
- Chiribella, G., Kristjánsson, H.: Quantum Shannon theory with superpositions of trajectories. Proc. R. Soc. A 475, 20180903 (2019). doi:10.1098/rspa.2018.0903; arXiv:1812.05292.
- Choi, M.-D.: Completely positive linear maps on complex matrices. Linear Algebra Appl. 10, 285–290 (1975).
- Cohen, T., Welling, M.: Group equivariant convolutional networks. Proc. ICML 2016, PMLR 48, 2990–2999 (2016).
- Cong, I., Choi, S., Lukin, M. D.: Quantum convolutional neural networks. Nat. Phys. 15(12), 1273–1278 (2019). doi:10.1038/s41567-019-0648-8.
- Cubitt, T. S., Eisert, J., Wolf, M. M.: The complexity of relating quantum channels to master equations. Commun. Math. Phys. 310(2), 383–418 (2012). doi:10.1007/s00220-011-1402-y; arXiv:0908.2128.
- Davies, E. B., Lewis, J. T.: An operational approach to quantum probability. Commun. Math. Phys. 17, 239–260 (1970).
- Diamond, S., Boyd, S.: CVXPY: A Python-embedded modeling language for convex optimization. J. Mach. Learn. Res. 17(83), 1–5 (2016).
- E, W.: A proposal on machine learning via dynamical systems. Commun. Math. Stat. 5(1), 1–11 (2017). doi:10.1007/s40304-017-0103-z.
- Finzi, M., Welling, M., Wilson, A. G.: A practical method for constructing equivariant multilayer perceptrons for arbitrary matrix groups. Proc. ICML 2021, PMLR 139, 3318–3328 (2021).
- Flammia, S. T., Wallman, J. J.: Efficient estimation of Pauli channels. ACM Trans. Quantum Comput. 1(1), 1–32 (2020). doi:10.1145/3408039; arXiv:1907.12976.
- Fletcher, A. S., Shor, P. W., Win, M. Z.: Optimum quantum error recovery using semidefinite programming. Phys. Rev. A 75, 012338 (2007). doi:10.1103/PhysRevA.75.012338.
- Frank, R. L., Lieb, E. H.: Monotonicity of a relative Rényi entropy. J. Math. Phys. 54, 122201 (2013).
- Geiger, M., Smidt, T.: e3nn: Euclidean neural networks. arXiv:2207.09453 (2022).
- Gorini, V., Kossakowski, A., Sudarshan, E. C. G.: Completely positive dynamical semigroups of N-level systems. J. Math. Phys. 17(5), 821–825 (1976).
- Grassl, M., Lu, S., Zeng, B.: Codes for simultaneous transmission of quantum and classical information. Proc. IEEE Int. Symp. Inf. Theory (ISIT 2017), 1718–1722 (2017). doi:10.1109/ISIT.2017.8006823; arXiv:1701.06963.
- Haber, E., Ruthotto, L.: Stable architectures for deep neural networks. Inverse Probl. 34(1), 014004 (2017; journal volume dated 2018). doi:10.1088/1361-6420/aa9a90; arXiv:1705.03341.
- Hashimoto, Y., Ikeda, M., Kadri, H.: Position: C\*-algebraic machine learning — moving in a new direction. Proc. ICML 2024, PMLR 235, 17667–17679 (2024). arXiv:2402.02637.
- Hashimoto, Y., Wang, Z., Matsui, T.: C\*-algebra net: a new approach generalizing neural network parameters to C\*-algebra. Proc. ICML 2022, PMLR 162, 8523–8534 (2022).
- Hastings, M. B.: Superadditivity of communication capacity using entangled inputs. Nat. Phys. 5(4), 255–257 (2009). doi:10.1038/nphys1224.
- He, K., Zhang, X., Ren, S., Sun, J.: Deep residual learning for image recognition. Proc. IEEE CVPR 2016, 770–778 (2016). doi:10.1109/CVPR.2016.90; arXiv:1512.03385.
- Huang, G., Liu, Z., van der Maaten, L., Weinberger, K. Q.: Densely connected convolutional networks. Proc. IEEE CVPR 2017, 2261–2269 (2017). doi:10.1109/CVPR.2017.243; arXiv:1608.06993.
- Huang, G., Sun, Y., Liu, Z., Sedra, D., Weinberger, K. Q.: Deep networks with stochastic depth. In: Computer Vision – ECCV 2016, Lecture Notes in Computer Science 9908, 646–661. Springer (2016). doi:10.1007/978-3-319-46493-0_39; arXiv:1603.09382.
- Huang, H.-Y., Broughton, M., Cotler, J., Chen, S., Li, J., Mohseni, M., Neven, H., Babbush, R., Kueng, R., Preskill, J., McClean, J. R.: Quantum advantage in learning from experiments. Science 376(6598), 1182–1186 (2022). doi:10.1126/science.abn7293; arXiv:2112.00778.
- Huang, H.-Y., Kueng, R., Preskill, J.: Predicting many properties of a quantum system from very few measurements. Nat. Phys. 16, 1050–1057 (2020). doi:10.1038/s41567-020-0932-7; arXiv:2002.08953.
- Huang, Z., Van Gool, L.: A Riemannian network for SPD matrix learning. Proc. AAAI 31(1) (2017). doi:10.1609/aaai.v31i1.10866.
- Jacobs, R. A., Jordan, M. I., Nowlan, S. J., Hinton, G. E.: Adaptive mixtures of local experts. Neural Comput. 3(1), 79–87 (1991). doi:10.1162/neco.1991.3.1.79.
- Jordan, R., Kinderlehrer, D., Otto, F.: The variational formulation of the Fokker–Planck equation. SIAM J. Math. Anal. 29(1), 1–17 (1998).
- Junge, M., Renner, R., Sutter, D., Wilde, M. M., Winter, A.: Universal recovery maps and approximate sufficiency of quantum relative entropy. Ann. Henri Poincaré 19(10), 2955–2978 (2018).
- Khaneja, N., Reiss, T., Kehlet, C., Schulte-Herbrüggen, T., Glaser, S. J.: Optimal control of coupled spin dynamics: design of NMR pulse sequences by gradient ascent algorithms. J. Magn. Reson. 172, 296–305 (2005).
- Khatun, A., Karim, A., Usman, M.: A hardware efficient quantum residual neural network without post-selection. arXiv:2604.06866 (2026), preprint.
- Kliesch, M., Barthel, T., Gogolin, C., Kastoryano, M., Eisert, J.: Dissipative quantum Church–Turing theorem. Phys. Rev. Lett. 107, 120501 (2011); erratum ibid. 109, 119904 (2012).
- Knill, E., Laflamme, R.: Theory of quantum error-correcting codes. Phys. Rev. A 55(2), 900–911 (1997). doi:10.1103/PhysRevA.55.900; arXiv:quant-ph/9604034.
- Kondor, R., Trivedi, S.: On the generalization of equivariance and convolution in neural networks to the action of compact groups. Proc. ICML 2018, PMLR 80, 2747–2755 (2018).
- Kunin, D., Sagastuy-Brena, J., Ganguli, S., Yamins, D., Tanaka, H.: Neural mechanics: symmetry and broken conservation laws in deep learning dynamics. ICLR 2021. arXiv:2012.04728.
- Kuperberg, G.: The capacity of hybrid quantum memory. IEEE Trans. Inf. Theory 49(6), 1465–1473 (2003). doi:10.1109/TIT.2003.811917; arXiv:quant-ph/0203105.
- Lezcano-Casado, M., Martínez-Rubio, D.: Cheap orthogonal constraints in neural networks: a simple parametrization of the orthogonal and unitary group. Proc. ICML 2019, PMLR 97, 3794–3803 (2019).
- Lindblad, G.: On the generators of quantum dynamical semigroups. Commun. Math. Phys. 48, 119–130 (1976). doi:10.1007/BF01608499.
- Lloyd, S., Schuld, M., Ijaz, A., Izaac, J., Killoran, N.: Quantum embeddings for machine learning. arXiv:2001.03622 (2020), preprint.
- Luchnikov, I. A., Vintskevich, S. V., Grigoriev, D. A., Filippov, S. N.: Machine learning non-Markovian quantum dynamics. Phys. Rev. Lett. 124, 140502 (2020). doi:10.1103/PhysRevLett.124.140502; arXiv:1902.07019.
- McArdle, S., Jones, T., Endo, S., Li, Y., Benjamin, S. C., Yuan, X.: Variational ansatz-based quantum simulation of imaginary time evolution. npj Quantum Inf. 5, 75 (2019). doi:10.1038/s41534-019-0187-2.
- McArdle, S., Yuan, X., Benjamin, S.: Error-mitigated digital quantum simulation. Phys. Rev. Lett. 122, 180501 (2019).
- McClean, J. R., Boixo, S., Smelyanskiy, V. N., Babbush, R., Neven, H.: Barren plateaus in quantum neural network training landscapes. Nat. Commun. 9, 4812 (2018). doi:10.1038/s41467-018-07090-4.
- Mittnenzweig, M., Mielke, A.: An entropic gradient structure for Lindblad equations and couplings of quantum systems to macroscopic models. J. Stat. Phys. 167(2), 205–233 (2017). doi:10.1007/s10955-017-1756-4; arXiv:1609.05765.
- Müller-Lennert, M., Dupuis, F., Szehr, O., Fehr, S., Tomamichel, M.: On quantum Rényi entropies: a new generalization and some properties. J. Math. Phys. 54, 122203 (2013). doi:10.1063/1.4838856.
- Nagy, A., Savona, V.: Variational quantum Monte Carlo method with a neural-network ansatz for open quantum systems. Phys. Rev. Lett. 122, 250501 (2019).
- Nguyen, Q. T., Schatzki, L., Braccia, P., Ragone, M., Coles, P. J., Sauvage, F., Larocca, M., Cerezo, M.: Theory for equivariant quantum neural networks. PRX Quantum 5, 020328 (2024).
- Nielsen, E., Gamble, J. K., Rudinger, K., Scholten, T., Young, K., Blume-Kohout, R.: Gate set tomography. Quantum 5, 557 (2021). doi:10.22331/q-2021-10-05-557; arXiv:2009.07301.
- Palais, R. S.: The principle of symmetric criticality. Commun. Math. Phys. 69, 19–30 (1979).
- Pérez-Salinas, A., Cervera-Lierta, A., Gil-Fuster, E., Latorre, J. I.: Data re-uploading for a universal quantum classifier. Quantum 4, 226 (2020). doi:10.22331/q-2020-02-06-226; arXiv:1907.02085.
- Petz, D.: Sufficient subalgebras and the relative entropy of states of a von Neumann algebra. Commun. Math. Phys. 105(1), 123–131 (1986).
- Petz, D.: Sufficiency of channels over von Neumann algebras. Q. J. Math. 39(1), 97–108 (1988). doi:10.1093/qmath/39.1.97.
- Petz, D.: Monotone metrics on matrix spaces. Linear Algebra Appl. 244, 81–96 (1996).
- Pimsner, M., Popa, S.: Entropy and index for subfactors. Ann. Sci. Éc. Norm. Supér. (4) 19(1), 57–106 (1986). doi:10.24033/asens.1504.
- Reeb, D., Kastoryano, M. J., Wolf, M. M.: Hilbert's projective metric in quantum information theory. J. Math. Phys. 52(8), 082201 (2011).
- Romero, J., Olson, J. P., Aspuru-Guzik, A.: Quantum autoencoders for efficient compression of quantum data. Quantum Sci. Technol. 2, 045001 (2017). doi:10.1088/2058-9565/aa8072; arXiv:1612.02806.
- Ronneberger, O., Fischer, P., Brox, T.: U-Net: convolutional networks for biomedical image segmentation. In: MICCAI 2015, Lecture Notes in Computer Science 9351, 234–241. Springer (2015). doi:10.1007/978-3-319-24574-4_28; arXiv:1505.04597.
- Sakhnenko, A., O'Meara, C., Ghosh, K. J. B., Mendl, C. B., Cortiana, G., Bernabé-Moreno, J.: Hybrid classical-quantum autoencoder for anomaly detection. Quantum Mach. Intell. 4, 27 (2022). doi:10.1007/s42484-022-00075-z; arXiv:2112.08869.
- Samach, G. O., et al.: Lindblad tomography of a superconducting quantum processor. Phys. Rev. Applied 18, 064056 (2022). doi:10.1103/PhysRevApplied.18.064056.
- Saxe, A. M., McClelland, J. L., Ganguli, S.: Exact solutions to the nonlinear dynamics of learning in deep linear neural networks. ICLR 2014. arXiv:1312.6120.
- Schuld, M.: Supervised quantum machine learning models are kernel methods. arXiv:2101.11020 (2021), preprint.
- Schumacher, B.: Quantum coding. Phys. Rev. A 51(4), 2738–2747 (1995). doi:10.1103/PhysRevA.51.2738.
- Spohn, H.: Entropy production for quantum dynamical semigroups. J. Math. Phys. 19, 1227–1230 (1978). doi:10.1063/1.523789.
- Stinespring, W. F.: Positive functions on C\*-algebras. Proc. Amer. Math. Soc. 6, 211–216 (1955).
- Stokes, J., Izaac, J., Killoran, N., Carleo, G.: Quantum natural gradient. Quantum 4, 269 (2020). doi:10.22331/q-2020-05-25-269.
- Su, W., Boyd, S., Candès, E. J.: A differential equation for modeling Nesterov's accelerated gradient method: theory and insights. J. Mach. Learn. Res. 17(153), 1–43 (2016).
- Torlai, G., Mazzola, G., Carrasquilla, J., Troyer, M., Melko, R., Carleo, G.: Neural-network quantum state tomography. Nat. Phys. 14(5), 447–450 (2018). doi:10.1038/s41567-018-0048-5; arXiv:1703.05334.
- Torlai, G., Wood, C. J., Acharya, A., Carleo, G., Carrasquilla, J., Aolita, L.: Quantum process tomography with unsupervised learning and tensor networks. Nat. Commun. 14, 2858 (2023). doi:10.1038/s41467-023-38332-9; arXiv:2006.02424.
- Tóth, G., Wieczorek, W., Gross, D., Krischek, R., Schwemmer, C., Weinfurter, H.: Permutationally invariant quantum tomography. Phys. Rev. Lett. 105, 250403 (2010). doi:10.1103/PhysRevLett.105.250403; arXiv:1005.3313.
- Vicentini, F., Biella, A., Regnault, N., Ciuti, C.: Variational neural-network ansatz for steady states in open quantum systems. Phys. Rev. Lett. 122, 250503 (2019).
- Wang, K., Li, X.: Simulation-assisted learning of open quantum systems. Quantum 8, 1407 (2024). doi:10.22331/q-2024-07-11-1407; arXiv:2307.03858.
- Wang, S., Fontana, E., Cerezo, M., Sharma, K., Sone, A., Cincio, L., Coles, P. J.: Noise-induced barren plateaus in variational quantum algorithms. Nat. Commun. 12, 6961 (2021). doi:10.1038/s41467-021-27045-6.
- Watrous, J.: *The Theory of Quantum Information*. Cambridge University Press (2018). ISBN 978-1-107-18056-7.
- Wen, J., Huang, Z., Cai, D., Qian, L.: Enhancing the expressivity of quantum neural networks with residual connections. Commun. Phys. 7, 220 (2024). doi:10.1038/s42005-024-01719-1; arXiv:2401.15871.
- Wen, Z., Yin, W.: A feasible method for optimization with orthogonality constraints. Math. Program. 142, 397–434 (2013). doi:10.1007/s10107-012-0584-1.
- Wibisono, A., Wilson, A. C., Jordan, M. I.: A variational perspective on accelerated methods in optimization. Proc. Natl. Acad. Sci. USA 113(47), E7351–E7358 (2016). doi:10.1073/pnas.1614734113.
- Wiseman, H. M., Milburn, G. J.: *Quantum Measurement and Control*. Cambridge University Press (2009). doi:10.1017/CBO9780511813948.
- Wolf, M. M., Cirac, J. I.: Dividing quantum channels. Commun. Math. Phys. 279, 147–168 (2008). doi:10.1007/s00220-008-0411-y; arXiv:math-ph/0611057.
- Wolf, M. M., Eisert, J., Cubitt, T. S., Cirac, J. I.: Assessing non-Markovian quantum dynamics. Phys. Rev. Lett. 101, 150402 (2008).
- Wootters, W. K., Zurek, W. H.: A single quantum cannot be cloned. Nature 299, 802–803 (1982). doi:10.1038/299802a0.
- Zhao, A., Miyake, A.: Group-theoretic error mitigation enabled by classical shadows and symmetries. npj Quantum Inf. 10, 57 (2024). doi:10.1038/s41534-024-00854-5; arXiv:2310.03071.

Library-internal references: the companion paper on trace characterisations (`man/trace_envelopes/`), the
proposition "Multiplicative composites: rapidities add" and Theorem A (the two-contrast Gardner envelope, used in
Section 8); the library paper (`man/article/`), Section 6.1 (Rényi
thresholds) and the use case "noisy circuits with a conserved charge"; `scripts/renyi_threshold.py` (Stinespring and
state parametrisations, 50-digit certification); the modules `optimize`, `measurement`, `permutational` and
`inclusions` (Pimsner–Popa constant and index).
