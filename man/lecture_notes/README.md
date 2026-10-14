**English** · [Español](es/README.md) · [Français](fr/README.md) · [Deutsch](de/README.md) · [中文](zh/README.md) · [日本語](ja/README.md) · [Русский](ru/README.md)

# Finite-dimensional von Neumann algebras: lecture notes

This course teaches the theory of finite-dimensional von Neumann algebras — the algebras
$\bigoplus_cM_{k_c}(\mathbb C)$ — together with their uses in quantum physics, quantum information, random matrix
theory and machine learning.

It is written for two kinds of readers at once:

- **Mathematicians** get a self-contained text: every theorem, formula, constant and
  inequality that is stated is proved in the course — in the chapter itself or in one of the
  appendices A–F, cited by number — from the prerequisites listed below. Numbers that are
  only computed (for instance large-$n$ level statistics) are labelled *numerical
  observation* and are never used as premises or test targets. References are historical
  attribution or further reading, never a substitute for a proof.
- **Physicists and engineers** get the physical meaning of each notion in boxes marked
  *Physical meaning*. These boxes are physical
  interpretation, not theorems of the course: nothing proved elsewhere depends on them.
  Working code they can adapt is in the [tutorial](../tutorial/README.md).

Chapters 1–13 have the same structure: definitions and results with complete proofs, worked examples, and
exercises with solutions. After every result that the library `torch_vn_algebra` computes or checks, a line
*Computation:* or *In the library:* links to the place in the [tutorial](../tutorial/README.md) where this is done.

## Contents

| | chapter | main results |
|---|---|---|
| 0 | [What this course is about](00_overview.md) | the algebra as a data type: block structure = batched tensors; cost; what exists and what is new |
| 1 | [Matrix algebras and functional calculus](01_matrix_algebras.md) | spectral theorem, functional calculus, polar decomposition, trace-norm duality, a GUE moment |
| 2 | [Finite-dimensional von Neumann algebras](02_von_neumann_algebras.md) | double commutant theorem, structure theorem, centre, symmetry groups and commutants, conditional expectations |
| 3 | [Traces, states and entropy](03_traces_and_states.md) | classification of traces, densities, entropy of sectors, Klein's inequality, monotonicity under the conditional expectation onto the centre, Michelson contrast |
| 4 | [Symmetries and representations](04_symmetries.md) | Schur's lemma, charge sectors, $U(1)$ and $SU(2)$, spin systems, real, complex and quaternionic type |
| 5 | [Random matrices in an algebra](05_random_matrices.md) | Haar measure on $\prod_cU(k_c)$, Weingarten calculus, level statistics and why sectors must not be mixed |
| 6 | [Quantum channels](06_quantum_channels.md) | Kraus and Stinespring, duals with respect to a trace, maps between sectors, decoherence |
| 7 | [Tensor products and entanglement](07_tensor_products.md) | partial trace, Page's formula, fusion of charges |
| 8 | [Dynamics](08_dynamics.md) | Schrödinger and Lindblad equations, jumps between sectors, post-selection |
| 9 | [Spin chains](09_spin_chains.md) | sectors of magnetisation and momentum, Lanczos, Bethe ansatz checks, disorder |
| 10 | [Quantum computing](10_quantum_computing.md) | decoherence-free subspaces, operator-algebra error correction, constrained mixers |
| 11 | [Optimisation and trace inequalities](11_optimisation.md) | constrained parametrisations, sampling versus optimisation |
| 12 | [Equivariant machine learning](12_equivariant_learning.md) | equivariant layers as commutants |
| 13 | [Jordan algebras: why quantum mechanics uses *-algebras](13_jordan_algebras.md) | Jordan product and quadratic representation, spin factors and the Lorentz cone, the Jordan–von Neumann–Wigner classification, the transpose and the partial transpose, a trace inequality in Jordan form |
| A | [Facts from analysis and probability](A_analysis.md) | power series and null sets, Gaussian vectors, spherical measures, eigenvalue inequalities, Isserlis–Wick |
| B | [Schur–Weyl duality and the Weingarten calculus](B_schur_weyl.md) | Schur–Weyl duality, Weingarten functions, Wigner–Eckart theorem, Haar measure on $O(k)$ |
| C | [Eigenvalue densities and Page's formula](C_eigenvalue_densities.md) | Weyl integration formula, GOE/GUE, Wishart and induced densities, Page's formula |
| D | [Generators of quantum dynamical semigroups](D_gksl.md) | the GKSL theorem, generators on direct sums, Frigerio's theorem, Bény–Kempf–Kribs |
| E | [Classification of Euclidean Jordan algebras](E_jordan_classification.md) | Peirce decomposition, Hurwitz and Frobenius theorems, coordinatisation, the Albert algebra, the Jordan–von Neumann–Wigner theorem |
| F | [Symmetric cones, order isomorphisms and Jordan maps](F_cones_and_order.md) | cones of squares, the Koecher–Vinberg theorem, Kadison's theorem, Jordan automorphisms of $\bigoplus_cM_{k_c}$ |

<a id="prerequisites"></a>

## Prerequisites

The course assumes only the following, and proves everything else.

- Linear algebra over $\mathbb R$ and $\mathbb C$ (finite-dimensional), and the fundamental
  theorem of algebra.
- Standard real analysis and multivariable calculus: limits, continuity, compactness in
  $\mathbb R^n$, differentiation, the mean value theorem, Taylor's formula, the inverse and
  implicit function theorems.
- Standard measure and integration theory: Lebesgue measure, product measures and
  Fubini–Tonelli, dominated and monotone convergence, change of variables in $\mathbb R^n$,
  pushforward measures, densities, independence.
- Elementary probability: expectation, variance, independence, the Gaussian distribution and
  its moment generating function, the Markov and Chebyshev inequalities, the law of large
  numbers.

No physics and no functional analysis are assumed.

## Further reading

M. Takesaki, *Theory of Operator Algebras I* (Springer, 2002); V. Jones,
*von Neumann algebras*, lecture notes; D. Petz, *Quantum Information Theory and Quantum
Statistics* (Springer, 2008); J. Watrous, *The Theory of Quantum Information* (Cambridge,
2018); M. L. Mehta, *Random Matrices* (Elsevier, 2004).

The library `torch_vn_algebra` and its use in the context of these notes are the subject of the
[tutorial](../tutorial/README.md). The companion paper of the library is
[`man/article/main.pdf`](https://github.com/andre-a-no/von_neumann_type_i/blob/main/man/article/main.pdf).
