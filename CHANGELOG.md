# Changelog

All notable changes to `torch_vn_algebra`. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versions follow [Semantic Versioning](https://semver.org/) from 1.0.0 on (the public API is everything exported by
`torch_vn_algebra/__init__.py` and documented in the module docstrings).

## [Unreleased]

### Added
- A guide for AI coding assistants, installed with the package:
  - `python -m torch_vn_algebra.guide` prints `torch_vn_algebra/guide.md`: the data layout, conventions, common
    mistakes, a task-to-function table and recipes run by the tests;
  - `python -m torch_vn_algebra.guide --api [module]` prints every public signature of the installed version.
- `AGENTS.md` and `CLAUDE.md`, instructions for agents that work in the repository.
- `skills/torch-vn-algebra/SKILL.md`, a Claude Code skill.
- The website serves `llms-full.txt` and a page with the guide. `llms.txt` now lists every module.

## [1.2.0] – 2026-10-25

Backward compatible: everything of 1.1.0 works unchanged. New: conversion to QuTiP and dynamiqs, channels and
Lindblad generators as trainable parameters, measurements and instruments, the projection lattice and measures on it.

### Added
- `anyons` module (`AnyonChain` exported): chains of SU(2)_k anyons in the fusion-path basis (paths on A_{k+1}),
  open (one sector per total charge) or periodic; Temperley–Lieb generators `generator(i)`, dense `hamiltonian` and
  `sparse_hamiltonian` (krylov.SparseSectorHamiltonian), `low_energies` (dense or restarted Lanczos with deflation),
  the tower TL_{L−1} ⊂ TL_L as an `Inclusion` in standard form (Jones index 4cos²(π/(k+2)), Markov trace = quantum
  dimensions), conformal data `central_charge`, `velocity`, `kac`. Tests `tests/test_anyons.py` (Temperley–Lieb
  relations, tower, Lanczos against dense, c = 7/10 of the golden chain, an Ising Kac ratio);
  `scripts/anyon_chain.py` (Table in the library paper).
- `direct_integral` module (`DirectIntegral` exported): finite type I algebras with a continuous centre,
  ∫^⊕_X M_{n(x)} dμ(x), discretised by nodes and weights (`grid`: Gauss–Legendre, midpoint or periodic product rules
  with a density; `torus`: a Brillouin zone with dk/2π; `monte_carlo`; or any partition), one sector per node:
  `field`, `scalar`, `indicator` (central projection of a measurable set), `integrate`, `measure`, `trace`,
  `normalized_trace`, `state` (normal states ∫ Tr(ρ(x) ·) dμ), `ess_sup`. Tests `tests/test_direct_integral.py`
  (SSH ground-state energy against the elliptic-integral formula; the direct-integral theorem of the trace paper);
  `scripts/direct_integral_demo.py`.
- `projections` module: the projection lattice (`join` = P ∨ Q, `meet` = P ∧ Q, `complement`, `support`, `ranks`,
  Haar `random_projection` with prescribed rank per sector, `bloch_vector` of a 2×2 block); `ProjectionMeasure`,
  finitely additive measures on the projections of ⊕M_k: linear (a state, `from_state`), weighted dimension functions
  (`dimension`), frame functions on 2×2 blocks (the only non-linear measures, by Gleason's theorem), dispersion-free
  measures (`dispersion_free`, which exist only on blocks of size ≤ 2, by Kochen–Specker), `modular_defect`,
  `gleason_state`, `is_linear`, `is_dimension_function`; `spectral_measure(A)`, the spectral measure of a self-adjoint
  operator (eigenvalues merged across sectors, or one outcome per eigenvector for a batch).
- `measurement.PVM` (exported): a projective POVM with real values per outcome; `observable`, `expectation`,
  `variance`, spectral projections `projection(predicate)`, `pushforward(f)` (the spectral measure of f(A)).
  Tests `tests/test_projections.py`, including the modular-law and measure propositions of the trace paper.
- `learn` module: quantum channels as trainable networks. `ChannelLayer` (Stinespring isometries, within one algebra
  or between algebras with allowed sector transitions), `UnitaryLayer`, `LindbladLayer` (continuous depth),
  `Residual` ((1−p)·id + p·Φ, the completely positive skip connection; n blocks with p = t/n tend to exp(t(Φ − id))),
  `PostSelect` (a projection of M, e.g. a set of sectors; records the acceptance), `FeedForward` (instrument outcome in the centre, one
  branch per outcome); losses and `count_parameters`. All are `torch.nn.Module`s; every parameter value is a channel.
- `measurement` module (`POVM` and `Instrument` exported at top level): POVMs with a one-time check of positivity
  and normalisation, Born probabilities, batched multinomial sampling of counts (reproducible with a generator) and
  log-likelihoods; constructors `POVM.sectors` (measurement of the charge), `POVM.basis` (rank-one
  projections, the computational basis of a `SpinChain`), `POVM.from_unitary` (rotated bases for tomography) and
  `POVM.combine` (randomised measurements, outcomes (setting, result)). `Instrument` (sector-preserving Kraus
  operators per outcome; `lueders`, `from_povm`): probabilities, conditional `post_states`, `average_channel` (a
  `Channel`) and `as_channel`, an `InterSectorChannel` from M into ⊕_j M (sectors (j, c)) whose centre records the
  outcome. `fisher_information(prob_fn, params)`: classical Fisher information J^T diag(1/p) J by autograd, shared or
  per-sample parameters, complex parameters as real and imaginary parts. Tests `tests/test_measurement.py` (agreement
  with `born_probabilities`, `sector_probabilities`, `lueders_update` and QuTiP's computational basis, sampling
  statistics, the analytic Fisher information of a qubit).
- `interop` module: conversion of operators, states, channels and Lindblad jump operators to and from dense
  matrices, `qutip.Qobj` and dynamiqs arrays (`to_dense`, `from_dense`, `channel_to_dense`, `to_qutip`, `from_qutip`,
  `jumps_to_qutip`, `to_dynamiqs`, `jumps_to_dynamiqs`). `layout=chain` places a `SpinChain` in the computational basis
  in QuTiP's convention (site 0 leftmost, |0> = spin up). Maps between algebras are converted with the pinching onto the
  output sectors (Kraus operators P_d K_i). `from_dense` rejects matrices outside the algebra. Tests check the
  conventions and that `lindblad_evolve` agrees with `qutip.mesolve` and `dynamiqs.mesolve`. Optional extra
  `interop = ["qutip>=5"]`.
- `optimize.StinespringParam`: channels as parameters. For every input sector an isometry V_c = Z_c (Z_c^* Z_c)^{-1/2}
  (Haar-random at the start), Kraus operators its blocks; every parameter value is a trace-preserving channel of one
  algebra (`Channel`) or between algebras with a prescribed set of allowed sector transitions (`InterSectorChannel`).
- `optimize.LindbladParam`: GKSL generators as parameters (a self-adjoint H and free jump operators inside the sectors
  and, for a prescribed set of allowed transitions, between them); a test learns a qubit generator from trajectories
  and predicts a new trajectory to 1e-4.
- `optimize.fit` / `FitResult`: batched multistart minimisation that leaves the parameters at the fitted point (the
  counterpart of `extremize` for learning); a test recovers a channel from 12 input/output pairs to 1e-5.
- `SpinChain.bond(i, j, Jxy, Jz)`: the two-site XXZ term on any pair of sites, for bond-resolved Hamiltonians and
  controls.

### Fixed
- `dynamics.lindblad_evolve(method='kraus')` returned NaN gradients when the normalisation operator of a step had a
  degenerate spectrum (for example M = 1 without jumps). Its inverse square root is now computed with the
  Daleckii–Krein backward pass of `functional.spectral_apply`, so gradients with respect to the Hamiltonian and the
  jump rates are finite and agree with finite differences in all three methods (new test).

## [1.1.0] – 2026-10-18

Backward compatible: everything of 1.0.0 works unchanged.

### Added
- `channels.SectorJump` (exported at top level): sparse jump operators that map basis vectors to multiples of basis
  vectors (spin flips S^±_i, S^z_i, particle loss). `dynamics.lindblad_evolve` (all three methods) and `lindblad_rhs`
  apply L ρ L^* as a weighted permutation of the entries, the dual L^* Z L as a gather and L^*L as a diagonal, in
  O(Σ_c k_c²) per application instead of O(C² k_max³) for a dense `InterSectorChannel`; `SectorJump.scaled(coef)` with a
  per-sample coefficient gives one jump per batch element (a sweep over noise rates in one evolution);
  `to_inter_sector()` converts to the dense form. Tests `tests/test_sector_jump.py`.
- `SpinChain.sector_jump(i, kind)` with `kind` in `'lower'`, `'raise'`, `'sz'`.
- Script `scripts/noisy_xy_mixer.py` and the main use case of the library paper: symmetry verification for a noisy
  XY mixer (Hamming-weight sectors, T1 at zero and finite temperature, dephasing), with exact statements checked to
  round-off.
- Script `scripts/architecture_bench.py`: the measured effect of the design decisions (known spectra, eigenvalues
  without eigenvectors, spectral cache, one call per block size, direct trace formulas, sparse jumps).
- Scripts for the companion paper on trace-characterising inequalities (`trace_contrast.py`, `trace_montecarlo.py`,
  `contrast_operations.py`, `trace_envelopes_figures.py`) with tests `tests/test_trace_contrast.py` and
  `tests/test_contrast_operations.py`.
- Module `permutational`: permutation-invariant ensembles of N two-level systems in ⊕_j M_{2j+1} ⊗ 1_{m_j}
  (`PermutationalSystem`, `PermutationalLindbladian`), with local dissipators as maps between spin sectors; tests against
  the full 2^N Lindblad generator and `qutip.piqs`. Script `scripts/permutational.py` and paper section "Use cases":
  lifetimes of a boundary time crystal under local emission and dephasing.
- Script `scripts/trace_characterizations.py` and paper use case: violation envelopes V(δ) of eleven inequalities
  that characterise the trace (Gardner, Bikchentaev and coauthors) for φ = Tr(W·) against the Michelson contrast δ of W;
  sharp values (Kantorovich constant, a/b, Specht ratio) proved where possible, violations certified in mpmath.
  Extended to sixteen inequalities: Tikhonov's subadditivity, a Powers–Størmer-type inequality (Hoa, Osaka, Toan),
  Hölder and Golden–Thompson inequalities for projections, and Peierls–Bogoliubov.
- Script `scripts/symmetry_resolved.py` and paper use case: symmetry-resolved entanglement of XXZ ground states
  (L = 24, eight anisotropies in one batched Lanczos run), with ρ_A assembled block by block in the fused product
  algebra; equipartition and its breaking, the Luttinger parameter from the charge statistics, free-fermion checks.
- Script `scripts/connes_stormer.py` and paper subsection "Use case: the Connes–Størmer entropy of two MASAs":
  H(Δ_n | UΔ_nU*) computed by batched optimisation over partitions of unity; closed form for n = 2 (one-dimensional
  maximisation, H > Choda's h(U) unless |u_11|² ∈ {0, 1/2, 1}); for n = 3, 4, 5 general partitions beat the best
  orthonormal basis found for some Haar-random U (values certified in mpmath); part of `scripts/run_all.sh`.
  Tests `tests/test_connes_stormer.py`.

### Changed
- The papers are split into four, each with its own file names and arXiv build script: the library paper
  (`man/article/vn_algebra_pytorch`, now a technical description of the architecture with the noisy-circuit use case),
  the trace-characterisation paper with S. A. Abed (`man/trace_envelopes/trace_inequalities_contrast`, which now also
  contains the envelopes of the three unnormalised trace functionals of Abed et al. and their proofs), the
  Connes–Størmer paper (`man/connes_stormer`) and the BTC/SRE paper (`man/btc_sre`).
- `scripts/trace_characterizations.py` covers 25 inequalities.

## [1.0.0] – 2026-10-15

First release, the version described in the paper (arXiv:2606.15882).

Requirements: Python ≥ 3.8, PyTorch ≥ 1.13, NumPy ≥ 1.21; the test suite runs in CI with these oldest versions
and with the latest ones on Python 3.8–3.13.

### Added
- `dynamics.lindblad_evolve(..., method='rk4' | 'kraus' | 'exact')` and `LINDBLAD_METHODS`: time steps that are
  quantum channels for every step size (course Theorem 8.16; RK4, the default, is not positive, Proposition 8.15).
- `krylov.krylov_error_bound`; `krylov.evolve(dt=None)` chooses adaptive steps from a proved error bound (course
  Theorem 9.18, Corollary 9.19) and `expm_multiply_step` returns that bound as `err`.
- `Decomposition.certificate()` and `Certificate`: a-posteriori certificate of a numerical Wedderburn decomposition
  (backward error, algebra gap, tolerance window; course Section 4.9).
- Theory tests (`tests/test_theory_*.py`, `tests/test_review_*.py`, `tests/test_article_claims.py`), single-precision
  and gradient tests, cross-checks against QuTiP, QuSpin, GAP and SymPy (`tests/test_crosscheck_external.py`).
- Documentation website (MkDocs, `mkdocs.yml`, `.github/workflows/pages.yml`).
- `TypeIAlgebra.operator(..., validate=True)` checks property flags passed by the user (they are trusted otherwise).
- `Operator.apply_function(f, per_block=True)` for functions that couple the eigenvalues of one block.
- `holomorphic_function(..., atol=...)`: absolute scale of the error estimate where f(X) may vanish.
- `schrodinger_rk4(..., algebra=...)` and `lindblad_rhs(..., algebra=...)` restrict raw tensors to the blocks.
- `states.log_partition_function`; release metadata (`CITATION.cff`, `.zenodo.json`), contribution guidelines.
- Module `inclusions`: unital inclusions N ⊂ M by inclusion matrices (`Inclusion`, with an optional batched unitary),
  the trace-preserving conditional expectation for any faithful trace (`ConditionalExpectation`), Markov trace and
  index ‖Λ‖², Watatani index and quasi-bases, Pimsner–Popa constants (order and complete order), basic construction,
  Jones projection and tower, relative commutant, composition and Bratteli diagrams, standard inclusions and the U(1)
  spin-chain tower A_L ⊂ A_{L+1} (index 4cos²(π/(2L+4))). Example `examples/example_inclusions.py`.
- Module `decompose`: numerical Wedderburn decomposition of the von Neumann algebra generated by matrices or by a
  finite group (`wedderburn`, `from_group`, `from_generators`, `Decomposition`), after Murota, Kanno, Kojima and
  Kojima (2010), with verification of every result; commutant bases, twirls, permutation representations and total
  spin. Example `examples/example_decompose.py` (Schur–Weyl duality).
  `Decomposition.real_structure()` and `real_basis()`: for an algebra closed under complex conjugation (real
  generators), the type of every simple summand of the real algebra — real M_k(R), complex M_k(C) (a pair of
  conjugate sectors) or quaternionic M_k(H) — from the sign of J conj(J) = ±1, and a real orthogonal basis adapted
  to it (Takagi factorisation); course Section 4.7.
- Module `su2`: `SU2Chain`, the SU(2)-invariant algebra ⊕_S M_{m_S} of spin-1/2 chains with highest-weight
  (Yamanouchi) bases in the U(1) sectors of `SpinChain`, Heisenberg Hamiltonians, full spectra with quantum numbers,
  ground states, exact Clebsch–Gordan coefficients, Wigner–Eckart reduced matrix elements and spin correlations.
  Example `examples/example_su2.py`.
- Module `divergences`: Petz and sandwiched Rényi divergences with exact support conventions and their limits
  (Umegaki, D_max, D_0), per-sector decomposition, Petz, rotated and universal recovery maps, noncommutative L^p and
  Kosaki norms. Example `examples/example_divergences.py`.
- Module `jordan`: Jordan product, quadratic representation, triple product and Jordan powers on operators; the
  transpose as a Jordan automorphism that is not completely positive; Jordan subalgebras of Herm_n generated by
  matrices, compared with the self-adjoint part of the generated *-algebra; spin factors V_m (spectral decomposition,
  functional calculus, Lorentz cone, representation by anticommuting symmetries). Example `examples/example_jordan.py`.
- Script `scripts/renyi_threshold.py` (paper, Section 6.1): locates the data-processing thresholds α = 1/2
  (sandwiched) and α = 2 (Petz) by batched multistart optimisation, with violations certified in mpmath; part of
  `scripts/run_all.sh`.
- Paper: section "Subalgebras, inclusions and symmetries" and subsection "Rediscovering the data-processing
  thresholds of the Rényi divergences".

### Fixed (independent correctness audit against dense NumPy/SciPy references)
- Gradients of functions of operators at (nearly) degenerate spectra (divided differences cancelled to 0; later the
  switching rule between quotient and midpoint was wrong for shifted spectra; now Proposition 1.13 of the course).
- `decompose`: the commutant ignored `tol`; noise-scaled thresholds saturated for noise >~ 1e-3.
- `divergences`: orthogonal supports gave finite values instead of +inf; NaN gradient of the Umegaki entropy at a
  zero eigenvalue.
- `states.gibbs_state`, `dynamics.schrodinger`, `dynamics.von_neumann` carried no gradient in H.
- `krylov.evolve`: steps far too long for large ||H|| (error 0.27), too short for a constant shift; float32 overflow
  of norms. `HilbertSpace.norm` overflowed in float32.
- `cost`: the timing probe overestimated ~100x (false `CostWarning`); Python lists became float32.
- `inverse()` and `power(p < 0)` of self-adjoint operators: the support threshold is per block again (small
  sectors next to large sectors of the same size were set to zero).
- `krylov.evolve` rounded the time grid to float32 (phase error of order E·3e-9).
- Lanczos stopped early when H carried a large constant shift (uniform field) and could return an excited state;
  the breakdown threshold is now tied to the machine epsilon. Zero initial vectors give zero, not NaN.
- Lindblad evolution in real algebras with jump operators failed (complex effective Hamiltonian).
- Garbage in the padding of raw tensors changed the results of `lindblad_rhs`/`lindblad_evolve`/`schrodinger_rk4`,
  `DensityMatrix(matrix=...)` and `DensityMatrix.from_density`.
- Property checks used an absolute tolerance floor; operators of small norm were misclassified (e.g. 1e-6 times a
  nilpotent as self-adjoint, -1e-6·1 as positive).
- `0 * A` kept stale property flags; `michelson_contrast` depended on the scale of A; `frobenius_norm` overflowed in
  single precision; `apply_function` of normal operators returned tensors attached to an internal graph.
- `holomorphic_function`: default contour for blocks that are scalar up to rounding.
- `gibbs_state` with negative beta (overflow); `random_mixed_unitary_channel` with invalid weights;
  `Channel.from_superoperator` tolerance in double precision; `random_density_matrix(measure='bures')` is rejected
  for real algebras (1 + O is singular for det O = -1).
- `Channel.from_superoperator` silently replaced a map that does not preserve Hermiticity by the Hermitian part of
  its Choi matrix (e.g. X ↦ iX became the zero map); it now raises. The complete-positivity tolerance is relative to
  each Choi matrix, not to the largest one of the same block size in the call.
- `HilbertSpace.norm/normalize` for flat vectors; `optimize.extremize` returns NaN (not ±inf) for starts without a
  finite value.

### Changed
- Minimum requirements lowered from Python 3.9 / PyTorch 2.0 to Python 3.8 / PyTorch 1.13 / NumPy 1.21
  (fallbacks for `torch.dtype.itemsize`, `check_invariants` of sparse tensors, overflow-safe complex scaling and
  version-independent gradients at the boundary of the positive cone). The workarounds are active only on the
  PyTorch versions that need them; on current versions the code paths and timings are unchanged. The wheel is pure
  Python and installs on Python 3.8; building from source needs Python ≥ 3.9 (setuptools ≥ 77, PEP 639 metadata).
- Decompositions run as one batched call per distinct block size, split into one chunk per thread on a multi-core
  CPU; no host synchronisations in the computational paths (GPU).
- `optimize`: spectra parametrised to reach the ends of [lambda_min, 1]; sign patterns of self-adjoint operators are
  searched; the best iterate of every start is returned.

## 0.x – development versions (not released)

- **0.5** (2026-10-10): functional calculus at four levels (self-adjoint, normal, holomorphic, general) with
  Daleckii–Krein gradients and a spectral cache.
- **0.4** (2026-10-09): sparse sector Hamiltonians, batched Lanczos and Krylov evolution, momentum sectors, Lindblad
  jumps between sectors, NumPy/SciPy baselines.
- **0.3** (2026-10-09): charge-conserving (fused) tensor products, U(1) spin chains.
- **0.2** (2026-10-09): channels, states, dynamics, maps between sectors, tensor products and partial traces, double
  precision, cost safeguards, constrained optimisation.
- **0.1** (2026-10-09): the code of the first version of the paper (arXiv:2606.15882v1), with its Haar sampler and
  spectrum bugs fixed.

[1.2.0]: https://github.com/andre-a-no/von_neumann_type_i/releases/tag/v1.2.0
[1.1.0]: https://github.com/andre-a-no/von_neumann_type_i/releases/tag/v1.1.0
[1.0.0]: https://github.com/andre-a-no/von_neumann_type_i/releases/tag/v1.0.0
