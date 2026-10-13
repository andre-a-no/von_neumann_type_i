**English** · [Español](es/README.md) · [Français](fr/README.md) · [Deutsch](de/README.md) · [中文](zh/README.md) · [日本語](ja/README.md) · [Русский](ru/README.md)

# Finite-dimensional von Neumann algebras: a course with `torch_vn_algebra`

This course teaches the theory of finite-dimensional von Neumann algebras — the algebras
$\bigoplus_cM_{k_c}(\mathbb C)$ — together with their uses in quantum physics, quantum
information, random matrix theory and machine learning, and shows at every step how to
compute with them in `torch_vn_algebra`.

It is written for two kinds of readers at once:

- **Mathematicians** get complete proofs of every result that is used. Nothing is quoted
  without proof except standard linear algebra, the fundamental theorem of algebra and
  elementary facts of analysis and probability (limits, convexity, the Gaussian moment
  generating function); where a result is stated before it is proved, or only applied, this
  is said explicitly.
- **Physicists and engineers** get the physical meaning of each notion in boxes marked
  *Physical meaning*, and working code they can adapt.

Chapters 1–12 have the same structure: theory with proofs; the library code for each
result; a check of the code against an exact result (a theorem, a closed formula or an
independent computation); a section *What to compare with* that names the standard tools
for the same task; and exercises with solutions.

All code in the course is executed by the test suite (`tests/test_tutorial.py`), so every
example works with the current version of the library.

## Contents

| | chapter | main results |
|---|---|---|
| 0 | [What this course is about](00_overview.md) | the algebra as a data type: block structure = batched tensors; cost; what exists and what is new |
| 1 | [Matrix algebras and functional calculus](01_matrix_algebras.md) | spectral theorem, functional calculus, polar decomposition, trace-norm duality, a GUE moment |
| 2 | [Finite-dimensional von Neumann algebras](02_von_neumann_algebras.md) | double commutant theorem, structure theorem, centre, symmetry groups and commutants, conditional expectations |
| 3 | [Traces, states and entropy](03_traces_and_states.md) | classification of traces, densities, entropy of sectors, Klein's inequality, monotonicity under the conditional expectation onto the centre, Michelson contrast |
| 4 | [Symmetries and representations](04_symmetries.md) | Schur's lemma, charge sectors, $U(1)$ and $SU(2)$, spin systems |
| 5 | [Random matrices in an algebra](05_random_matrices.md) | Haar measure on $\prod_cU(k_c)$, Weingarten calculus, level statistics and why sectors must not be mixed |
| 6 | [Quantum channels](06_quantum_channels.md) | Kraus and Stinespring, duals with respect to a trace, maps between sectors, decoherence |
| 7 | [Tensor products and entanglement](07_tensor_products.md) | partial trace, Page's formula, fusion of charges |
| 8 | [Dynamics](08_dynamics.md) | Schrödinger and Lindblad equations, jumps between sectors, post-selection |
| 9 | [Spin chains](09_spin_chains.md) | sectors of magnetisation and momentum, Lanczos, Bethe ansatz checks, disorder |
| 10 | [Quantum computing](10_quantum_computing.md) | decoherence-free subspaces, operator-algebra error correction, constrained mixers |
| 11 | [Optimisation and trace inequalities](11_optimisation.md) | constrained parametrisations, sampling versus optimisation |
| 12 | [Equivariant machine learning](12_equivariant_learning.md) | equivariant layers as commutants |

## Prerequisites

Linear algebra over $\mathbb C$ (eigenvalues, inner products, unitary matrices), basic
probability, and Python. No physics and no functional analysis are assumed.

## How to run the code

```bash
pip install -e .                # from the root of the repository
python -m pytest tests/test_tutorial.py -q
```

Every code block of a chapter can be pasted, in order, into one Python session or a
Jupyter notebook. The blocks of one chapter share their variables.

## Further reading

M. Takesaki, *Theory of Operator Algebras I* (Springer, 2002); V. Jones,
*von Neumann algebras*, lecture notes; D. Petz, *Quantum Information Theory and Quantum
Statistics* (Springer, 2008); J. Watrous, *The Theory of Quantum Information* (Cambridge,
2018); M. L. Mehta, *Random Matrices* (Elsevier, 2004). For the physical background in
mathematical language see [Physics for mathematicians](../physics_for_mathematicians.md). The
companion paper of the library is
[`paper/main.pdf`](https://github.com/andre-a-no/von_neumann_type_i/blob/main/paper/main.pdf).
