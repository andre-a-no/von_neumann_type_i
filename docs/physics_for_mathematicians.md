**English** · [Español](physics_for_mathematicians.es.md) · [Français](physics_for_mathematicians.fr.md) · [Deutsch](physics_for_mathematicians.de.md) · [中文](physics_for_mathematicians.zh.md) · [日本語](physics_for_mathematicians.ja.md) · [Русский](physics_for_mathematicians.ru.md)

# Physics for mathematicians: a dictionary for `torch_vn_algebra`

This guide is for readers who come from operator algebras and want to know what the physical models
in the library and the paper mean. Everything is stated in terms of the algebra
M = ⊕_c M_{k_c}(C) acting on H = ⊕_c C^{k_c}. For each physical notion we say what it is in M and
where it lives in the library.

---

## 0. The main idea

A physical symmetry is a unitary representation π of a group G on H. Everything compatible with the
symmetry (the Hamiltonian, the time evolution, the observables under a superselection rule) lies in the
**commutant** π(G)'. By the double commutant theorem and Schur's lemma,

    H ≅ ⊕_λ V_λ ⊗ C^{m_λ},     π(G)' ≅ ⊕_λ M_{m_λ}(C),

so π(G)' is a finite-dimensional Type I von Neumann algebra. Its sectors are the irreducible
representations λ, the block sizes are the multiplicities m_λ, and the centre consists of functions of
the "charges". A physicist who models a system with a symmetry therefore computes in an algebra of the
form (1), usually without calling it that. For a single conserved quantity Q (a self-adjoint operator)
the sectors are the eigenvalues of Q and the blocks are its eigenspaces.

This also answers the question "why not use one large matrix?". The embedding M ⊂ M_N (N = Σk_c) is
exact, but it loses
- **cost**: Σk_c³ instead of N³; for C equal blocks this is a factor C² in time;
- **measure**: the unitary group of M is Π U(k_c), a null set in U(N). A Haar-random unitary in U(N)
  is not in M, and a gradient step in M_N leaves M;
- **traces**: on M_N the trace is unique up to a factor, on M the faithful traces form the
  C-parameter family Σ w_c Tr_c. "Trace preserving" and dual maps depend on the weights;
- **the centre**: the sector projections are the classical part of the system (Section 2).
  Conditional probabilities and post-selection are operations with the centre, and with one matrix
  they have to be reconstructed from block indices.

The case C = 1 is supported as well: it is batched linear algebra in M_n.

---

## 1. States, observables, measurements

| physics | mathematics | in the library |
|---|---|---|
| observable (energy, spin, particle number) | self-adjoint A ∈ M | `Operator` |
| possible measurement outcomes | spectrum of A | `eigenvalues`, `eigh` |
| (mixed) state | normal positive functional ω, ω(1) = 1; density ρ ≥ 0, Tr ρ = 1 | `DensityMatrix` |
| pure state, "wave function" | vector ψ ∈ H, ρ = \|ψ⟩⟨ψ\| (a minimal projection) | `vector_in_sector`, `basis_state` |
| expectation value | ω(A) = Tr(ρA) | `expectation` |
| probability of outcome λ | ω(P_λ), P_λ the spectral projection | `apply_function` |
| state after a measurement (Lüders rule) | ρ ↦ PρP / Tr(Pρ) | `lueders_update`, `condition_on` |
| von Neumann entropy | S(ρ) = −Tr ρ log ρ | `entropy` |
| temperature, thermal state | ρ = e^{−βH}/Tr e^{−βH} (Gibbs state = KMS state) | `gibbs_state` |
| infinite temperature | tracial state τ / τ(1) | `tracial_state` |

**Physical meaning of the three traces.**
- Tr_blunt (weights 1) is the physical trace on H. Normalised, it is the infinite-temperature state:
  all N basis states are equally likely, and sector c has probability k_c/N.
- Tr_norm (weights 1/k_c): the uniform (microcanonical) state inside each sector separately.
- τ_vN (weights 1/(C k_c)): the uniform mixture of the microcanonical states of the sectors. Each
  *sector* (not each basis state) has probability 1/C.

The restriction of any state to the centre is a classical probability distribution on the sectors
(`sector_probabilities`). This is what physicists call "the probability of charge q".

---

## 2. Superselection, decoherence, classical and quantum

- **Superselection rule** (Wick–Wightman–Wigner): superpositions of states with different charge
  (electric charge, fermion parity) cannot be observed. Mathematically, the admissible observables form
  the commutant of the charge, that is the algebra M rather than all of B(H).
- **The centre of M** consists of quantities that can be measured without disturbance and commute with
  everything: the "classical" variables. If all k_c = 1, the system is purely classical (a commutative
  algebra, i.e. a probability space on C points). If C = 1, it is purely quantum. The general case is a
  hybrid.
- **Decoherence / einselection** (Zurek): interaction with an environment suppresses the off-diagonal
  blocks ρ_{cd}, c ≠ d. This is the conditional expectation E: B(H) → M, ρ ↦ Σ P_c ρ P_c, and the
  effective algebra of observables shrinks to M. In the library operators already lie in M, so the
  result of this conditional expectation is where computations start.

---

## 3. Closed-system dynamics

| physics | mathematics | in the library |
|---|---|---|
| Hamiltonian | self-adjoint H ∈ M | — |
| Schrödinger equation i dψ/dt = Hψ | ψ(t) = e^{−itH}ψ(0) | `dynamics.schrodinger`, `propagator` |
| von Neumann equation dρ/dt = −i[H, ρ] | ρ(t) = U ρ U*, U = e^{−itH} | `von_neumann` |
| Heisenberg picture | A(t) = U* A U (an automorphism of M) | `solve_operator_ode` |
| conserved quantity | Q with [H, Q] = 0 | sectors |

If H ∈ M then e^{−itH} ∈ M, and the dynamics never leaves the sectors. Each sector evolves
independently, and the exponential is computed block by block.

---

## 4. Open systems and channels

- A **quantum channel** is a completely positive trace-preserving map (Schrödinger picture). Its
  adjoint (Heisenberg picture) is completely positive and unital. Kraus form:
  Φ(ρ) = Σ_i K_i ρ K_i*. In the library: `Channel`.
- A **channel between sectors** (`InterSectorChannel`) describes a process that changes the charge:
  particle loss, decay, an electron tunnelling out. Its Kraus blocks are K^{dc}: C^{k_c} → C^{m_d}.
  Trace preservation and the dual map depend on the trace weights on input and output, hence the
  factors w_out/w_in.
- The **Lindblad equation** (GKSL) generates a semigroup of channels e^{tL}:
  dρ/dt = −i[H, ρ] + Σ_j γ_j (L_j ρ L_j* − ½{L_j* L_j, ρ}).
  The L_j are "jump operators"; each describes one way the environment acts on the system.
  - **T1** (energy relaxation, amplitude damping): L = σ⁻ = |0⟩⟨1|, a qubit decays from the excited
    to the ground state. The number of excitations drops by one, so this is a jump **between**
    sectors. Example: `site_amplitude_damping`, `lowering`.
  - **T2** (dephasing): L = σ^z. It is diagonal and does not change the sector; it only destroys
    phases. This is a jump **within** a sector.
  - T1 and T2 are the standard coherence times quoted for every quantum computer.

---

## 5. Spin chains

- **Spin ½ on a site**: H_site = C², basis |↑⟩, |↓⟩ (or |1⟩, |0⟩). Pauli matrices σ^x, σ^y, σ^z;
  S^a = σ^a/2; S^± = S^x ± iS^y are the raising and lowering operators.
- **A chain of L spins**: H = (C²)^{⊗L}, dim = 2^L. This exponential growth is why physicists need
  symmetries.
- **XXZ Hamiltonian**: H = Σ_⟨ij⟩ [J/2 (S_i⁺S_j⁻ + S_i⁻S_j⁺) + JΔ S_i^z S_j^z] + Σ_i h_i S_i^z.
  The first term swaps neighbouring ↑↓ (a "particle hop"), the second is their interaction, the third an
  external field. Δ = 1 is the Heisenberg model (a magnet), Δ = 0 is equivalent to free fermions (by the
  Jordan–Wigner transformation), which gives exact checks.
- **Conserved charge**: S^z_tot, or N = number of ↑ (the Hamming weight of a bit string). Every term of
  H conserves it, so H ∈ ⊕_{N=0}^{L} M_{C(L,N)}: the blocks have binomial sizes.
- **Sector fusion** (`fused_tensor_product`): under the tensor product charges add. This is the fusion
  rule of U(1): (⊕_q M_{n_q}) ⊗ (⊕_r M_{m_r}) ⊂ ⊕_s M_{Σ_{q+r=s} n_q m_r}. `SpinChain` builds its basis
  this way, site by site.
- **Momentum** (`momentum_algebra`): a ring also has the translation symmetry T (the group Z_L). Joint
  eigenspaces (N, k), with e^{ik} an eigenvalue of T, give a finer decomposition with blocks about L
  times smaller.
- The **ground state** is the eigenvector of H with the smallest eigenvalue, i.e. the behaviour at zero
  temperature. For large blocks it is found with the Lanczos method (`krylov`).
- What the paper checks, and why these are standard tests:
  - ground-state energies of the Heisenberg ring (known exact values, Bethe ansatz);
  - the **Marshall sign rule**: for the Heisenberg model on a bipartite lattice the ground state has
    signs (−1)^{number of ↑ on one sublattice}. It is a theorem, which makes it a clean test;
  - **SU(2) nesting**: the Heisenberg model has full SU(2) symmetry, so for N < L/2 the spectrum of
    sector N is contained in that of sector N+1 (multiplets);
  - **domain-wall melting**: |↑…↑↓…↓⟩ at Δ = 0 spreads with a known profile;
  - **Peschel's formula**: the entanglement entropy of free fermions from the correlation matrix.
- **Disorder and chaos** (`random_field_heisenberg`, `level_spacing_ratio`): random fields
  h_i ∈ [−W, W]. For small W the system is "chaotic" and levels repel as in the GOE, ⟨r⟩ ≈ 0.53. For
  large W it localises (many-body localisation, MBL) and levels are independent as for a Poisson
  process, ⟨r⟩ ≈ 0.386. Level statistics only make sense **inside one sector**: mixing sectors produces
  Poisson statistics artificially. This is a good physical argument for working in the algebra with its
  sectors.

---

## 6. Entanglement and subsystems

- A composite system A+B: H_A ⊗ H_B, algebra M_A ⊗ M_B.
- The **reduced state** ρ_A = Tr_B ρ is the restriction of the functional ω to the subalgebra M_A ⊗ 1.
  In terms of densities it is a conditional expectation (the partial trace). In the library:
  `partial_trace`, `reduced_state`.
- The **entanglement entropy** S(ρ_A) of a pure ρ measures the quantum correlation between A and B.
- **Page's formula**: the exact mean entropy of a subsystem of a random pure state. It tests that the
  sampler is Haar distributed.
- **Symmetry-resolved entanglement**: if the charge is conserved, ρ_A is block diagonal and its entropy
  splits into a "classical" part (the entropy of the distribution on the centre) and a "quantum" part
  (the average entropy of the blocks). This is a direct use of the centre.

---

## 7. Quantum computing

- **Qubit** = C², n qubits = (C²)^{⊗n}; a **gate** is a unitary; a **circuit** is a product of gates;
  a **measurement** in the computational basis uses the spectral projections of σ^z.
- **QAOA** (quantum approximate / alternating operator ansatz): alternate e^{−iγ H_C} (H_C is
  diagonal and encodes the objective) and e^{−iβ H_M} (the "mixer"). For problems with the constraint
  "exactly k ones" (choose k items out of n) one uses the **XY mixer** Σ(X_iX_j + Y_iY_j) =
  2Σ(σ⁺σ⁻ + h.c.). It conserves the Hamming weight, so the whole circuit lies in ⊕_N M_{C(n,N)} and
  the feasible solutions are exactly the sector N = k.
- **Noise**: T1 leaves the feasible sector (a one is lost), T2 does not.
- **Post-selection**: measure the Hamming weight and discard wrong runs. This is the Lüders rule with a
  **central** projection P_k: ρ ↦ P_kρP_k / Tr(P_kρ), with success probability ω(P_k). This is what
  `ex9_qaoa_xy_mixer.py` does. For pure T1 noise post-selection returns the ideal state exactly
  (F = 1), because every T1 jump leaves the sector.
- Decoherence-free subspaces and operator-algebra error correction: protected information is encoded in
  one block (or one tensor factor of a block) of an algebra that the noise does not touch. This, too,
  is the language of Type I algebras (Knill–Laflamme–Viola, Bény–Kempf–Kribs).

---

## 8. Random matrices, real and complex

This is directly related to the parameter `complex_valued`.
- An **antiunitary symmetry** is time reversal T. If T² = 1 (spinless particles, no magnetic field),
  H is **real** in a suitable basis: GOE, orthogonal ensembles, COE. Real arithmetic is enough; it
  halves the memory and is faster (see the timing table in the paper).
- Without such a symmetry (a magnetic field, complex phases, momentum k ≠ 0, π) H is complex: GUE, CUE.
  Complex numbers are required.
- If T² = −1 (half-integer spin with spin–orbit coupling), the case is quaternionic: GSE, CSE.
- The XXZ model with real fields is real. Momentum sectors with k ∉ {0, π} are complex. The dynamics
  e^{−itH} is always complex, even for real H. This is the "what suffices when" table of the paper.
- Double precision is needed when small differences matter: degenerate spectra, identities checked to
  1e−12, long time evolutions, Lanczos with many steps. Single precision (and TF32 on tensor cores) is
  enough for statistics over many samples.

---

## 9. Machine learning

An equivariant linear layer is a map W: V → V that commutes with a group representation (translations
for convolutions, permutations for graph networks, rotations for molecules). By Schur's lemma the space
of such W is ⊕_λ M_{m_λ} with multiplicities m_λ, so an equivariant layer is parametrised by an element
of an algebra of the form (1). The library provides random initialisation with the right measure,
spectral constraints and gradients for such parametrisations.

---

## 10. Scope: what the library is and is not

- The library does **not** produce new physics. It is a validated tool for computations in algebras
  with sectors. The physical models in the paper are known systems with known answers, used to check
  the tool.
- For one large system (a single ground state at L = 30, tensor networks, simulation of large circuits)
  the specialised packages (QuSpin, ITensor/TeNPy, Qiskit) are the better choice.
- Its niche is many moderately sized operators with explicit sector structure: disorder ensembles,
  statistics over random states and channels, extremal operators under spectral constraints, channels
  between sectors with consistent traces, together with automatic differentiation and GPUs.
- The original mathematical question (trace inequalities, Michelson contrast) is an example of an
  extremal problem where Monte Carlo sampling fails and constrained optimisation works.

---

## Further reading

- M. A. Nielsen, I. L. Chuang, *Quantum Computation and Quantum Information* — ch. 2 (postulates),
  ch. 8 (channels, T1/T2).
- M. M. Wilde, *Quantum Information Theory* — channels and states, close to operator algebras.
- A. W. Sandvik, "Computational studies of quantum spin systems", AIP Conf. Proc. 1297 (2010) —
  sectors and exact diagonalisation in practice.
- H.-P. Breuer, F. Petruccione, *The Theory of Open Quantum Systems* — the Lindblad equation.
- F. Haake, *Quantum Signatures of Chaos* — random matrices, time reversal, GOE/GUE/GSE.
- S. D. Bartlett, T. Rudolph, R. W. Spekkens, Rev. Mod. Phys. 79, 555 (2007) — superselection and
  reference frames, readable for mathematicians.
