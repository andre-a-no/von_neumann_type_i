# Quantum channels as neural networks: a roadmap

*Design proposal for the next major version of `torch_vn_algebra` (current: 1.1.0). Status: discussion draft.
Nothing described under "proposed" exists in the library yet; every code block marked **NOT IMPLEMENTED** is a
sketch of an API, not working code. References are listed at the end with full bibliographic data; all of them were
checked, and items we could not check were left out.*

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

## 1. Dictionary: channel → network

| channel notion | network notion | in the library now | missing |
|---|---|---|---|
| Stinespring isometry $`V:\mathcal H_{\rm in}\to\mathcal H_{\rm out}\otimes\mathbb C^r`$, then $`\operatorname{Tr}_{\mathbb C^r}`$ | linear (isometric) layer followed by a contraction / pooling over an environment index (Stinespring 1955) | `channels.random_channel` (Haar isometry); `composite.partial_trace_channel` | a trainable Stinespring layer (only `scripts/renyi_threshold.py:StinespringParam`, one block) |
| Kraus rank / environment dimension $`r`$ | width | `Channel.kraus_rank` | rank control (truncation, growth) during training |
| composition $`\Phi_L\circ\dots\circ\Phi_1`$ | depth | `Channel.compose`, `InterSectorChannel.compose` | a sequential container with algebra checks |
| $`e^{t\mathcal L}`$ of a GKSL generator | continuous depth (neural ODE, Chen et al. 2018; E 2017) | `lindblad_channel`, `lindblad_evolve`, differentiable | a trainable generator |
| CPTP splitting step $`\Psi_h(M^{-1/2}\,\cdot\,M^{-1/2})`$ (`dynamics._kraus_step`) | residual-type discretisation that stays inside the constraint set | `lindblad_evolve(method='kraus')` | exposure as a layer |
| conditional expectation $`E:\mathcal M\to\mathcal N\subset\mathcal M`$ | pooling onto a subalgebra (coarse-graining) | `inclusions.ConditionalExpectation`, `channels.center_expectation` | — (usable as fixed layers) |
| sector (block) structure | structured sparsity and weight sharing (strong/weak symmetry of Lindbladians: Buča and Prosen 2012) | the storage format itself; `decompose` turns a symmetry into an algebra | sparse lists of allowed pairs $`(d,c)`$ |
| Choi matrix $`J(\Phi)`$ | the weight tensor | `Channel.choi`, `superoperator`, `from_superoperator` | Choi-based losses |
| TP / unitality | normalisation constraints (as for stochastic matrices) | `trace_preservation_error`, `unitality_error` | constraint-preserving parametrisations |
| instrument, post-selection on a central projection | the source of nonlinearity (normalisation); instruments: Davies and Lewis (1970) | `measurement_channel`, `lueders_update`, `DensityMatrix.condition_on`, `born_probabilities` | an `Instrument` class that keeps the outcomes |

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
classifying states.

## 2. Parametrisations in which every parameter is a valid channel

| family | parametrisation | gradient path | status |
|---|---|---|---|
| CPTP, Kraus rank $`r`$ | free $`Z_c\in\mathbb C^{(\sum_d m_d r)\times k_c}`$ per input sector; $`V_c=Z_c(Z_c^*Z_c)^{-1/2}`$ (polar), or QR, or Cayley updates (Wen and Yin 2013); Riemannian view in Absil et al. (2008) | inverse square root by `functional.spectral_apply`: Daleckii–Krein backward, finite at degenerate spectra | single block in `scripts/renyi_threshold.py`; inter-sector version missing |
| CPTP, Kraus form | $`K_i\mapsto K_iM^{-1/2}`$ (identical to the polar retraction above) | same | inside `dynamics._kraus_step` only |
| unitary channels in $`\mathcal M`$ | $`U=U_0e^{G-G^*}`$ (cf. Lezcano-Casado and Martínez-Rubio 2019; Arjovsky et al. 2016) | `expm` | `optimize.UnitaryParam` |
| GKSL generators | $`H=(A+A^*)/2`$; Kossakowski matrix $`K=LL^*\succeq0`$ in a jump basis $`F_\alpha`$ adapted to the sectors (Gorini, Kossakowski and Sudarshan 1976; Lindblad 1976) | differentiable `lindblad_evolve` (`'kraus'`: every step CPTP; `'exact'`: $`e^{h\mathcal L}`$) | jumps are inputs; no trainable generator |
| states | $`\rho=(1-\varepsilon)GG^*/\operatorname{Tr}GG^*+\varepsilon/d`$, or Gibbs $`e^{-\beta H}/Z`$ | `states.gibbs_state` | `StateParam` in scripts only |
| unital CPTP | alternate the TP and unital normalisations (an operator version of Sinkhorn scaling) | — | proposal; convergence not analysed |

The Stinespring and Lindblad parametrisations have **gauge symmetries**. The environment unitary
$`V\mapsto(1\otimes W)V`$, $`W\in U(r)`$, does not change $`\Phi`$, and unitary mixing of jump operators does not
change $`\mathcal L`$. The polar parametrisation is also invariant under $`Z\mapsto cZ`$ for $`c>0`$. These
symmetries make the Hessian singular, which matters for natural-gradient methods. They also give conservation laws
(Section 5), which serve as cheap correctness checks of a training run.

**Missing pieces:**
- an `nn.Module` wrapper;
- an inter-sector Stinespring layer with a sparse list of allowed pairs;
- a trainable Kossakowski matrix;
- moving `StinespringParam` and `StateParam` from `scripts/renyi_threshold.py` into the library.

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
This is a "deep composition approximates continuous depth" theorem with a known depth scaling.

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
  of the companion paper. *First experiment:* reproduce the Rényi thresholds with the library classes of Section 6
  instead of the script.

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

## 6. Proposed API: `torch_vn_algebra.learn`

The module is named `learn` to avoid confusion with `torch.nn`. All classes are `torch.nn.Module`s whose parameters
are unconstrained tensors. The constraint is enforced in `forward`, by retraction or normalisation, so any
`torch.optim` optimiser applies. Operators keep the library layout `(B, C, k_max, k_max)`, and channels are returned
as `Channel` / `InterSectorChannel` objects.

1. **`StinespringChannel(alg_in, alg_out, kraus_rank, pairs=None, retraction='polar'|'qr'|'cayley', trace_preserving=True)`**.
   Parameters: one complex tensor `Z` of shape `(B, r, D, C, m_max, k_max)`, masked to `pairs`, the allowed
   `(d, c)`. Methods: `forward(rho) -> Operator`, `channel() -> InterSectorChannel` (TP by construction),
   `num_parameters()` (the exact counts of Section 1), `gauge_momentum(velocity)`.
2. **`LindbladGenerator(alg, H=None, jump_basis=..., kossakowski='cholesky', controls=None, n_slices=1)`**. Trainable
   $`H`$, Cholesky factor $`L`$ of $`K=LL^*`$, and piecewise-constant control amplitudes. Methods: `evolve(rho, times,
   method='kraus')`, `channel(t)` (exact, small blocks), `trotter(n) -> ChannelSequential`.
3. **`ChannelSequential(*layers)`**. Composition with algebra checks. Layers may be the two classes above, fixed
   `Channel`s, `ConditionalExpectation`s (pooling) or `Instrument`s.
4. **`Instrument(alg_in, alg_out, n_outcomes, ...)`**, with the constructor `Instrument.postselect(alg, central_projection)`.
   Stored as an `InterSectorChannel` into $`\bigoplus_j\mathcal N`$. `forward` returns `(states, probabilities)`;
   post-selection renormalises, differentiably.
5. **`ChannelFlow(module, kind='gradient'|'natural'|'heavy_ball'|'nesterov_ode', metric='embedded'|'bures'|'kubo_mori', step)`**.
   The training dynamics of Section 5 as an optimiser with retraction and vector transport. It logs energy, momentum
   maps and scale invariants.

Losses are plain functions in `learn.losses`: infidelity, trace distance, Umegaki and Rényi divergences (wrapping
`divergences`), entanglement and process fidelity from the Choi matrix, and the rapidity regulariser
$`\sum_l\eta(\Phi_l(1))`$. `learn.certify(channel, objective, digits=50)` generalises the mpmath recomputation of the
scripts.

```python
# NOT IMPLEMENTED: proposed API (target 1.2). Process tomography in U(1) sectors.
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
# NOT IMPLEMENTED: proposed API (target 1.3). Noise-aware control with post-selection on the weight w.
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
# NOT IMPLEMENTED: proposed API (target 1.3). Training dynamics with conservation checks (experiment X1).
model = learn.StinespringChannel(N_alg, M_alg, kraus_rank=2)
for kind in ('gradient', 'natural', 'heavy_ball', 'nesterov_ode'):
    flow = learn.ChannelFlow(model.reset(seed=0), kind=kind, metric='bures', step=1e-2)
    for k in range(2000):
        loss = -entanglement_fidelity(model.channel(), noise, code)
        flow.step(loss)
        log(kind, k, loss.item(), sdp_value - (-loss.item()),
            model.Z.norm()**2, flow.momentum_map('environment'))
```

## 7. Staged plan, risks, measures of success

| version | content | done when |
|---|---|---|
| **1.2** parametrised channels | `StinespringChannel` (inter-sector, polar/QR), `LindbladGenerator` (Cholesky Kossakowski), `learn.losses`; `StinespringParam` and `StateParam` moved from scripts | TP error $`\le10^{-12}`$ (double); `gradcheck` passes, including at degenerate spectra; `num_parameters()` matches the formulas of Section 1; `renyi_threshold.py` rewritten on the library reproduces its certified table |
| **1.3** composition and dynamics of training | `ChannelSequential`, `Instrument`, `ChannelFlow`, `learn.certify`; experiments X1–X3 | X1: every flow reaches the SDP optimum to $`10^{-8}`$, and the predicted conservation and decay laws hold to discretisation order; X3: the Carlen–Maas identity checked to $`10^{-10}`$ |
| **2.0** applications | sector-structured tomography and Lindbladian learning (a); control with post-selection (b); recovery (c); a linear-map type for rectangular blocks $`\operatorname{Hom}_G=\bigoplus_\lambda M_{n_\lambda\times m_\lambda}`$ between algebras (not CP) | (a): sample complexity against the dense Stinespring baseline in the style of Ahmed et al.; (b): post-selected infidelity of optimised versus constant schedules at $`n\le10`$; written up as a paper |

**Risks.**

- *Degenerate spectra.* First derivatives of $`f(\rho)`$ are finite (Daleckii–Krein), but derivatives through
  eigenvectors are not, and second derivatives are correct only at non-degenerate spectra. Natural-gradient metrics
  are affected.
- *Conditioning of $`M^{-1/2}`$ and $`\log\rho`$.* Monitor and regularise with the rapidity bound of Section 3. The
  `kraus` step already halves $`h`$ when $`M`$ is ill-conditioned.
- *Padding cost.* The inter-sector tensor costs $`r\,D\,C\,m_{\max}k_{\max}`$, against the useful
  $`r\sum_{(d,c)}m_dk_c`$. For chains with local jumps the waste is about $`C/3`$ before padding. Mitigations: sparse
  pair lists, grouping by block size, `SectorJump`.
- *Scalability.* Choi matrices and superoperators are $`k^2\times k^2`$ per block, so blocks are limited to a few
  tens. The library does not compete with escnn (Cesa et al. 2022) or e3nn (Geiger and Smidt 2022) on large equivariant
  feature spaces, or with circuit simulators on many qubits.
- *Optimisation landscape.* Non-convexity, gauge-singular Hessians, and contraction with depth (Section 3, item 4).
  For convex problems, an SDP should be preferred.
- *Over-claiming.* Most individual ingredients exist (Section 8). Claims of novelty need a literature check per item.

## 8. Positioning

- **Quantum neural networks.** Beer et al. (2020) build deep networks whose layers are CP maps: a unitary on input and
  output qudits, followed by a partial trace over the input. This is the Stinespring layer of Section 1, trained with
  fidelity losses. Quantum convolutional networks (Cong, Choi and Lukin 2019) pool by measurement and conditioned
  unitaries, which is an instrument followed by feed-forward. Equivariant QNNs (Nguyen et al. 2024) use the same
  commutant ideas as the library's `decompose`.
- **Libraries.** PennyLane (Bergholm et al.) and TensorFlow Quantum (Broughton et al.) differentiate gate circuits,
  with noise as fixed channels. dynamiqs (JAX) offers batched, differentiable Lindblad simulation on GPUs.
- **Channel and Lindbladian learning.** See Ahmed et al. (2023), Torlai et al. (2023), Samach et al. (2022) and
  Bairey et al. (2020).
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

Section 5 adds a testbed for the dynamics of training on channel manifolds, with falsifiable predictions: conservation
and decay laws, invariance of covariant sets, and the Carlen–Maas identity. Items marked [N] or "not found" must be
checked against the literature before any publication claims them.

## References

- Absil, P.-A., Mahony, R., Sepulchre, R.: *Optimization Algorithms on Matrix Manifolds*. Princeton University Press (2008). ISBN 978-0-691-13298-3, doi:10.1515/9781400830244.
- Ahmed, S., Quijandría, F., Kockum, A. F.: Gradient-descent quantum process tomography by learning Kraus operators. Phys. Rev. Lett. 130, 150402 (2023). doi:10.1103/PhysRevLett.130.150402; arXiv:2208.00812.
- Amari, S.: Natural gradient works efficiently in learning. Neural Comput. 10(2), 251–276 (1998). doi:10.1162/089976698300017746.
- Arjovsky, M., Shah, A., Bengio, Y.: Unitary evolution recurrent neural networks. Proc. ICML 2016, PMLR 48, 1120–1128 (2016).
- Bairey, E., Guo, C., Poletti, D., Lindner, N. H., Arad, I.: Learning the dynamics of open quantum systems from their steady states. New J. Phys. 22, 032001 (2020). doi:10.1088/1367-2630/ab73cd.
- Barnum, H., Knill, E.: Reversing quantum dynamics with near-optimal quantum and classical fidelity. J. Math. Phys. 43, 2097 (2002). arXiv:quant-ph/0004088.
- Beer, K., Bondarenko, D., Farrelly, T., Osborne, T. J., Salzmann, R., Scheiermann, D., Wolf, R.: Training deep quantum neural networks. Nat. Commun. 11, 808 (2020). doi:10.1038/s41467-020-14454-2.
- Bergholm, V., et al.: PennyLane: Automatic differentiation of hybrid quantum-classical computations. arXiv:1811.04968 (2018, rev. 2022).
- Bonet-Monroig, X., Sagastizabal, R., Singh, M., O'Brien, T. E.: Low-cost error mitigation by symmetry verification. Phys. Rev. A 98, 062339 (2018).
- Broughton, M., et al.: TensorFlow Quantum: A software framework for quantum machine learning. arXiv:2003.02989 (2020, rev. 2021).
- Buča, B., Prosen, T.: A note on symmetry reductions of the Lindblad equation: transport in constrained open spin chains. New J. Phys. 14, 073007 (2012).
- Carleo, G., Troyer, M.: Solving the quantum many-body problem with artificial neural networks. Science 355(6325), 602–606 (2017). doi:10.1126/science.aag2302.
- Carlen, E. A., Maas, J.: Gradient flow and entropy inequalities for quantum Markov semigroups with detailed balance. J. Funct. Anal. 273(5), 1810–1869 (2017). doi:10.1016/j.jfa.2017.05.003.
- Cesa, G., Lang, L., Weiler, M.: A program to build E(N)-equivariant steerable CNNs. ICLR 2022 (library escnn).
- Chen, R. T. Q., Rubanova, Y., Bettencourt, J., Duvenaud, D.: Neural ordinary differential equations. Advances in Neural Information Processing Systems 31 (NeurIPS 2018). arXiv:1806.07366.
- Choi, M.-D.: Completely positive linear maps on complex matrices. Linear Algebra Appl. 10, 285–290 (1975).
- Cohen, T., Welling, M.: Group equivariant convolutional networks. Proc. ICML 2016, PMLR 48, 2990–2999 (2016).
- Cong, I., Choi, S., Lukin, M. D.: Quantum convolutional neural networks. Nat. Phys. 15(12), 1273–1278 (2019). doi:10.1038/s41567-019-0648-8.
- Davies, E. B., Lewis, J. T.: An operational approach to quantum probability. Commun. Math. Phys. 17, 239–260 (1970).
- Diamond, S., Boyd, S.: CVXPY: A Python-embedded modeling language for convex optimization. J. Mach. Learn. Res. 17(83), 1–5 (2016).
- E, W.: A proposal on machine learning via dynamical systems. Commun. Math. Stat. 5(1), 1–11 (2017). doi:10.1007/s40304-017-0103-z.
- Finzi, M., Welling, M., Wilson, A. G.: A practical method for constructing equivariant multilayer perceptrons for arbitrary matrix groups. Proc. ICML 2021, PMLR 139, 3318–3328 (2021).
- Fletcher, A. S., Shor, P. W., Win, M. Z.: Optimum quantum error recovery using semidefinite programming. Phys. Rev. A 75, 012338 (2007). doi:10.1103/PhysRevA.75.012338.
- Frank, R. L., Lieb, E. H.: Monotonicity of a relative Rényi entropy. J. Math. Phys. 54, 122201 (2013).
- Geiger, M., Smidt, T.: e3nn: Euclidean neural networks. arXiv:2207.09453 (2022).
- Gorini, V., Kossakowski, A., Sudarshan, E. C. G.: Completely positive dynamical semigroups of N-level systems. J. Math. Phys. 17(5), 821–825 (1976).
- Hashimoto, Y., Wang, Z., Matsui, T.: C\*-algebra net: a new approach generalizing neural network parameters to C\*-algebra. Proc. ICML 2022, PMLR 162, 8523–8534 (2022).
- Hashimoto, Y., Ikeda, M., Kadri, H.: Position: C\*-algebraic machine learning — moving in a new direction. Proc. ICML 2024, PMLR 235, 17667–17679 (2024). arXiv:2402.02637.
- Hastings, M. B.: Superadditivity of communication capacity using entangled inputs. Nat. Phys. 5, 255 (2009). doi:10.1038/nphys1224.
- Huang, Z., Van Gool, L.: A Riemannian network for SPD matrix learning. Proc. AAAI 31(1) (2017). doi:10.1609/aaai.v31i1.10866.
- Jordan, R., Kinderlehrer, D., Otto, F.: The variational formulation of the Fokker–Planck equation. SIAM J. Math. Anal. 29(1), 1–17 (1998).
- Junge, M., Renner, R., Sutter, D., Wilde, M. M., Winter, A.: Universal recovery maps and approximate sufficiency of quantum relative entropy. Ann. Henri Poincaré 19(10), 2955–2978 (2018).
- Khaneja, N., Reiss, T., Kehlet, C., Schulte-Herbrüggen, T., Glaser, S. J.: Optimal control of coupled spin dynamics: design of NMR pulse sequences by gradient ascent algorithms. J. Magn. Reson. 172, 296–305 (2005).
- Kliesch, M., Barthel, T., Gogolin, C., Kastoryano, M., Eisert, J.: Dissipative quantum Church–Turing theorem. Phys. Rev. Lett. 107, 120501 (2011); erratum ibid. 109, 119904 (2012).
- Kondor, R., Trivedi, S.: On the generalization of equivariance and convolution in neural networks to the action of compact groups. Proc. ICML 2018, PMLR 80, 2747–2755 (2018).
- Kunin, D., Sagastuy-Brena, J., Ganguli, S., Yamins, D., Tanaka, H.: Neural mechanics: symmetry and broken conservation laws in deep learning dynamics. ICLR 2021. arXiv:2012.04728.
- Lezcano-Casado, M., Martínez-Rubio, D.: Cheap orthogonal constraints in neural networks: a simple parametrization of the orthogonal and unitary group. Proc. ICML 2019, PMLR 97, 3794–3803 (2019).
- Lindblad, G.: On the generators of quantum dynamical semigroups. Commun. Math. Phys. 48, 119–130 (1976). doi:10.1007/BF01608499.
- McArdle, S., Jones, T., Endo, S., Li, Y., Benjamin, S. C., Yuan, X.: Variational ansatz-based quantum simulation of imaginary time evolution. npj Quantum Inf. 5, 75 (2019). doi:10.1038/s41534-019-0187-2.
- McArdle, S., Yuan, X., Benjamin, S.: Error-mitigated digital quantum simulation. Phys. Rev. Lett. 122, 180501 (2019).
- McClean, J. R., Boixo, S., Smelyanskiy, V. N., Babbush, R., Neven, H.: Barren plateaus in quantum neural network training landscapes. Nat. Commun. 9, 4812 (2018). doi:10.1038/s41467-018-07090-4.
- Mittnenzweig, M., Mielke, A.: An entropic gradient structure for Lindblad equations and couplings of quantum systems to macroscopic models. J. Stat. Phys. 167(2), 205–233 (2017). doi:10.1007/s10955-017-1756-4; arXiv:1609.05765.
- Müller-Lennert, M., Dupuis, F., Szehr, O., Fehr, S., Tomamichel, M.: On quantum Rényi entropies: a new generalization and some properties. J. Math. Phys. 54, 122203 (2013). doi:10.1063/1.4838856.
- Nagy, A., Savona, V.: Variational quantum Monte Carlo method with a neural-network ansatz for open quantum systems. Phys. Rev. Lett. 122, 250501 (2019).
- Nguyen, Q. T., Schatzki, L., Braccia, P., Ragone, M., Coles, P. J., Sauvage, F., Larocca, M., Cerezo, M.: Theory for equivariant quantum neural networks. PRX Quantum 5, 020328 (2024).
- Palais, R. S.: The principle of symmetric criticality. Commun. Math. Phys. 69, 19–30 (1979).
- Petz, D.: Sufficient subalgebras and the relative entropy of states of a von Neumann algebra. Commun. Math. Phys. 105(1), 123–131 (1986).
- Petz, D.: Monotone metrics on matrix spaces. Linear Algebra Appl. 244, 81–96 (1996).
- Reeb, D., Kastoryano, M. J., Wolf, M. M.: Hilbert's projective metric in quantum information theory. J. Math. Phys. 52(8), 082201 (2011).
- Samach, G. O., et al.: Lindblad tomography of a superconducting quantum processor. Phys. Rev. Applied 18, 064056 (2022). doi:10.1103/PhysRevApplied.18.064056.
- Spohn, H.: Entropy production for quantum dynamical semigroups. J. Math. Phys. 19, 1227–1230 (1978). doi:10.1063/1.523789.
- Stinespring, W. F.: Positive functions on C\*-algebras. Proc. Amer. Math. Soc. 6, 211–216 (1955).
- Stokes, J., Izaac, J., Killoran, N., Carleo, G.: Quantum natural gradient. Quantum 4, 269 (2020). doi:10.22331/q-2020-05-25-269.
- Su, W., Boyd, S., Candès, E. J.: A differential equation for modeling Nesterov's accelerated gradient method: theory and insights. J. Mach. Learn. Res. 17(153), 1–43 (2016).
- Torlai, G., Wood, C. J., Acharya, A., Carleo, G., Carrasquilla, J., Aolita, L.: Quantum process tomography with unsupervised learning and tensor networks. Nat. Commun. 14 (2023). doi:10.1038/s41467-023-38332-9.
- Vicentini, F., Biella, A., Regnault, N., Ciuti, C.: Variational neural-network ansatz for steady states in open quantum systems. Phys. Rev. Lett. 122, 250503 (2019).
- Wang, S., Fontana, E., Cerezo, M., Sharma, K., Sone, A., Cincio, L., Coles, P. J.: Noise-induced barren plateaus in variational quantum algorithms. Nat. Commun. 12, 6961 (2021). doi:10.1038/s41467-021-27045-6.
- Watrous, J.: *The Theory of Quantum Information*. Cambridge University Press (2018). ISBN 978-1-107-18056-7.
- Wen, Z., Yin, W.: A feasible method for optimization with orthogonality constraints. Math. Program. 142, 397–434 (2013). doi:10.1007/s10107-012-0584-1.
- Wibisono, A., Wilson, A. C., Jordan, M. I.: A variational perspective on accelerated methods in optimization. Proc. Natl. Acad. Sci. USA 113(47), E7351–E7358 (2016). doi:10.1073/pnas.1614734113.
- Wolf, M. M., Cirac, J. I.: Dividing quantum channels. Commun. Math. Phys. 279, 147 (2008). doi:10.1007/s00220-008-0411-y.
- Wolf, M. M., Eisert, J., Cubitt, T. S., Cirac, J. I.: Assessing non-Markovian quantum dynamics. Phys. Rev. Lett. 101, 150402 (2008).

Library-internal references: the companion paper on trace characterisations (`man/trace_envelopes/`), the
proposition "Multiplicative composites: rapidities add"; the library paper (`man/article/`), Section 6.1 (Rényi
thresholds) and the use case "noisy circuits with a conserved charge"; `scripts/renyi_threshold.py` (Stinespring and
state parametrisations, 50-digit certification).
