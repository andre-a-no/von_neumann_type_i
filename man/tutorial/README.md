**English** · [Español](es/README.md) · [Français](fr/README.md) · [Deutsch](de/README.md) · [中文](zh/README.md) · [日本語](ja/README.md) · [Русский](ru/README.md)

# Tutorial: torch_vn_algebra in the context of the lecture notes

This tutorial shows how to compute with the finite-dimensional von Neumann algebras $`\bigoplus_cM_{k_c}(\mathbb C)`$
in the library `torch_vn_algebra`, in the context of the [lecture notes](../lecture_notes/README.md). Each tutorial accompanies
the chapter of the lecture notes with the same number: it shows the library code for the results proved there, checks
the code against an exact result (a theorem, a closed formula or an independent computation), and names, in a section
*What to compare with*, the standard tools for the same task. Each section starts with a line *Lecture notes:* that
links to the section and the results it uses. Numbers that are only computed (for instance large-$`n`$ level
statistics) are labelled *numerical observation* and are never used as premises or test targets. Besides the
[prerequisites](../lecture_notes/README.md#prerequisites) of the lecture notes, the tutorial assumes Python.

All code in the tutorial is executed by the test suite (`tests/test_man.py`), so every
example works with the current version of the library.

## How to run the code

```bash
pip install -e .                # from the root of the repository
python -m pytest tests/test_man.py -q
```

Every code block of a chapter can be pasted, in order, into one Python session or a
Jupyter notebook. The blocks of one chapter share their variables.

## Contents

| | chapter | main results |
|---|---|---|
| 0 | [What this course is about](00_overview.md) | the algebra as a data type: block structure = batched tensors; cost; what exists and what is new |
| 1 | [Matrix algebras and functional calculus](01_matrix_algebras.md) | spectral theorem, functional calculus, polar decomposition, trace-norm duality, a GUE moment |
| 2 | [Finite-dimensional von Neumann algebras](02_von_neumann_algebras.md) | double commutant theorem, structure theorem, centre, symmetry groups and commutants, conditional expectations |
| 3 | [Traces, states and entropy](03_traces_and_states.md) | classification of traces, densities, entropy of sectors, Klein's inequality, monotonicity under the conditional expectation onto the centre, Michelson contrast |
| 4 | [Symmetries and representations](04_symmetries.md) | Schur's lemma, charge sectors, $`U(1)`$ and $`SU(2)`$, spin systems, real, complex and quaternionic type |
| 5 | [Random matrices in an algebra](05_random_matrices.md) | Haar measure on $`\prod_cU(k_c)`$, Weingarten calculus, level statistics and why sectors must not be mixed |
| 6 | [Quantum channels](06_quantum_channels.md) | Kraus and Stinespring, duals with respect to a trace, maps between sectors, decoherence |
| 7 | [Tensor products and entanglement](07_tensor_products.md) | partial trace, Page's formula, fusion of charges |
| 8 | [Dynamics](08_dynamics.md) | Schrödinger and Lindblad equations, jumps between sectors, post-selection |
| 9 | [Spin chains](09_spin_chains.md) | sectors of magnetisation and momentum, Lanczos, Bethe ansatz checks, disorder |
| 10 | [Quantum computing](10_quantum_computing.md) | decoherence-free subspaces, operator-algebra error correction, constrained mixers |
| 11 | [Optimisation and trace inequalities](11_optimisation.md) | constrained parametrisations, sampling versus optimisation |
| 12 | [Equivariant machine learning](12_equivariant_learning.md) | equivariant layers as commutants |
| 13 | [Jordan algebras: why quantum mechanics uses *-algebras](13_jordan_algebras.md) | Jordan product and quadratic representation, spin factors and the Lorentz cone, the Jordan–von Neumann–Wigner classification, the transpose and the partial transpose, a trace inequality in Jordan form |
| A | [Facts from analysis and probability](A_analysis.md) | power series and null sets, Gaussian vectors, spherical measures, eigenvalue inequalities, Isserlis–Wick |
| B | [Schur–Weyl duality and the Weingarten calculus](B_schur_weyl.md) | Schur–Weyl duality, Weingarten functions, Wigner–Eckart theorem, Haar measure on $`O(k)`$ |
| C | [Eigenvalue densities and Page's formula](C_eigenvalue_densities.md) | Weyl integration formula, GOE/GUE, Wishart and induced densities, Page's formula |
| D | [Generators of quantum dynamical semigroups](D_gksl.md) | the GKSL theorem, generators on direct sums, Frigerio's theorem, Bény–Kempf–Kribs |
| E | [Classification of Euclidean Jordan algebras](E_jordan_classification.md) | Peirce decomposition, Hurwitz and Frobenius theorems, coordinatisation, the Albert algebra, the Jordan–von Neumann–Wigner theorem |
| F | [Symmetric cones, order isomorphisms and Jordan maps](F_cones_and_order.md) | cones of squares, the Koecher–Vinberg theorem, Kadison's theorem, Jordan automorphisms of $`\bigoplus_cM_{k_c}`$ |

The mathematics — definitions, theorems and complete proofs — is in the [lecture notes](../lecture_notes/README.md),
together with the prerequisites, further reading and the companion paper of the library.
