**English** · [Español](es/00_overview.md) · [Français](fr/00_overview.md) · [Deutsch](de/00_overview.md) · [中文](zh/00_overview.md) · [日本語](ja/00_overview.md) · [Русский](ru/00_overview.md)

# Tutorial 0. What this course is about

This tutorial accompanies [Chapter 0](../lecture_notes/00_overview.md) of the lecture notes, which identifies an element of a
von Neumann algebra $`\bigoplus_cM_{k_c}(\mathbb C)`$ with an array of shape $`(C,k,k)`$ and proves that its operations
are batched array operations. Here we see what this gives in the library `torch_vn_algebra`: what it costs,
a first computation, and how the library relates to existing software.

Von Neumann algebras are mostly studied with pen and paper. In finite dimensions every one of
them is $`*`$-isomorphic to a direct sum of full matrix algebras (Theorem 2.6), so in principle
anything can be computed — but when such an algebra is handled numerically, it is often
flattened into one large matrix, and its structure (the sectors, the centre, the many traces)
disappears from the computation. This course and the library `torch_vn_algebra` take the
opposite approach: **the algebra itself is the data type**, and computations in it are carried
out efficiently.

<a id="s-0-1"></a>

## 0.1 The observation

*Lecture notes:* [Section 0.1](../lecture_notes/00_overview.md#s-0-1), which proves, from [Theorem 2.6](../lecture_notes/02_von_neumann_algebras.md#r-2-6),
[Theorem 3.1](../lecture_notes/03_traces_and_states.md#r-3-1), [Theorem 5.3](../lecture_notes/05_random_matrices.md#r-5-3) and
[Corollary 5.4](../lecture_notes/05_random_matrices.md#r-5-4), that the first six rows of the following table are identities:

| in the algebra $`\mathcal M`$ | on the array of shape $`(B,C,k,k)`$ |
|---|---|
| product $`AB`$, sum, adjoint | batched matrix product, sum, conjugate transpose over $`(B,C)`$ |
| functional calculus $`f(A)`$ (Chapter 1) | batched eigendecomposition of the blocks |
| a trace $`\sum_cw_c\mathrm{Tr}\,A_c`$ (Chapter 3) | diagonal sum per block, weighted sum over the sector axis |
| the centre, central projections (Chapter 2) | one number per sector |
| conditional expectation onto the centre | block trace divided by $`k_c`$, broadcast back |
| Haar measure on the unitary group $`\prod_cU(k_c)`$ | batched QR decomposition of Gaussian blocks, with the phase correction of Theorem 5.3 |
| states, channels, dynamics (Chapters 3, 6, 8) | built from the operations above |

Two consequences make this more than bookkeeping.

- **Cost.** A spectral decomposition in $`\mathcal M`$ is computed block by block (row 2 of the
  table). The standard dense eigensolvers (LAPACK, which PyTorch calls on a CPU) use a number of
  floating-point operations that grows like the cube of the matrix size; this is a description of
  those algorithms, not a result of this course. In the model in which a $`k\times k`$ eigenproblem
  costs a constant times $`k^3`$, a spectral decomposition in $`\mathcal M`$ costs $`\sum_ck_c^3`$ against $`D^3`$
  for the embedding into $`M_D(\mathbb C)`$, $`D=\sum_ck_c`$; for $`C`$ equal blocks of size $`k`$ the ratio is
  $`(Ck)^3/(Ck^3)=C^2`$. Products run on the padded tensor; the schoolbook product of two $`k\times k`$
  matrices uses $`k^3`$ multiplications ($`k^2`$ entries, each a sum of $`k`$ products), so a product in
  $`\mathcal M`$ uses $`C\,k_{\max}^3`$, which equals $`\sum_ck_c^3`$ when all blocks have the same size. An
  element is stored as $`C\,k_{\max}^2`$ numbers against $`D^2`$; for equal blocks the ratio is
  $`(Ck)^2/(Ck^2)=C`$. What the code does: products, sums and traces are single tensor operations on
  the whole $`(B,C,k_{\max},k_{\max})`$ array; Hermitian eigendecompositions, singular values and QR
  decompositions are applied to the stack of all $`B\cdot n_k`$ blocks of each distinct size $`k`$ ($`n_k`$
  the number of sectors of that size) in one batched call — a single call when all sectors have the
  same size. On a CPU with $`t\ge2`$ PyTorch threads the library splits such a call into $`t`$ chunks run
  in parallel threads, provided the stack has at least $`2t`$ blocks and (number of blocks)$`\,\cdot k^3\ge2^{18}`$
  (module `cost`, function `_run`); the results are the same LAPACK calls on the same matrices. The
  aim is to keep a GPU, or all cores of a CPU, busy.
- **Correctness by construction.** Every array of shape $`(C,k_{\max},k_{\max})`$ with zero padding is
  an element of $`\mathcal M`$, so whatever is computed from the blocks stays in $`\mathcal M`$. The
  library's random unitaries are drawn by the construction of Theorem 5.3 independently in every block
  (`random_unitary_operator` calls the per-block sampler `random_unitary` once per block size, and
  Section 5.1 checks `random_unitary` against this construction), which is the Haar measure of
  $`\prod_cU(k_c)`$ by Corollary 5.4 (for ideal Gaussian random numbers); a Haar unitary of $`U(D)`$ lies
  in this subgroup with probability zero when $`C\ge2`$ (Section 5.1, after Corollary 5.4). A gradient
  step changes the blocks and so cannot leave the algebra. And the traces are the functionals
  $`\sum_cw_c\mathrm{Tr}\,A_c`$ of Theorem 3.1 with explicitly chosen weights (and the densities and
  conditional expectations are taken with respect to them), not the trace of some $`\mathbb C^n`$ on which
  $`\mathcal M`$ happens to act, which would fix $`w_c=m_c`$ (Section 2.3).

<a id="s-0-2"></a>

## 0.2 A first computation

*Lecture notes:* [Section 0.1](../lecture_notes/00_overview.md#s-0-1), row 2 of the table: $`f(A)=\bigoplus_cf(A_c)`$, one eigendecomposition per block.

We compare the same computation done in $`\mathcal M=\bigoplus_{c=1}^{64}M_{16}(\mathbb C)`$ and in its
embedding into $`M_{1024}(\mathbb C)`$: the square root of a random positive element, and its trace.
`TypeIAlgebra(n_factors, k_factors)` takes two lists of block sizes: the nominal sizes $`n_c`$ and
the sizes $`k_c\le n_c`$ of the blocks actually used (Section 1.5). Every algebra that the course
constructs with `TypeIAlgebra` has two equal lists.

```python
import time
import torch
from torch_vn_algebra import TypeIAlgebra

torch.manual_seed(0)
C, k = 64, 16
alg = TypeIAlgebra([k] * C, [k] * C, precision='double')       # M = (+)_{c=1}^{64} M_16(C)

A = alg.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(1, d), force_positive=True)
blocks = A.matrix                                                # build A before timing

t0 = time.perf_counter()
root_blocks = A.sqrt().matrix                                    # (1, 64, 16, 16): 64 small eigenproblems
t_blocks = time.perf_counter() - t0

big = torch.block_diag(*blocks[0])                               # the same element in M_1024(C)
t0 = time.perf_counter()
w, V = torch.linalg.eigh(big)                                    # one large eigenproblem
root_big = (V * w.sqrt()) @ V.conj().T
t_big = time.perf_counter() - t0

assert torch.allclose(torch.block_diag(*root_blocks[0]), root_big, atol=1e-10)   # row 2: f(A) = (+)_c f(A_c)
assert abs(A.sqrt().Tr_blunt().item() - root_big.trace().item()) < 1e-9
assert big.numel() == C * blocks.numel()                          # storage: a factor C = 64
print(f"in the algebra: {t_blocks * 1e3:.1f} ms, as one 1024 x 1024 matrix: {t_big * 1e3:.1f} ms")
```

The results agree to round-off (row 2 of the table), and the storage ratio is $`C`$, as counted
above. The 64 blocks go to the eigensolver as one batched call on a stack of shape $`(64,16,16)`$;
with $`2\le t\le32`$ PyTorch threads the library splits it into $`t`$ chunks ($`64\cdot16^3=2^{18}`$ meets
the threshold of module `cost` exactly), with one thread it does not. *Numerical observation*
(timings depend on the hardware; the code prints them): on two 4-core CPUs the algebra version took
3–10 ms and the large-matrix version 400–760 ms, a factor between several tens and a few hundred
(on the second machine with 4 threads and also with one). The measured factor is far below the
ratio $`C^2=4096`$ of the cubic cost model: an operation count is a model of the running time, not the running time,
and on these machines the small eigenproblems ran at a much lower rate of operations per second
than the large one. In the cost model the ratio $`C^2`$ grows with the number of sectors. With a
batch of $`B`$ elements the algebra version stays one batched call; the large-matrix version has to
solve $`B`$ eigenproblems of size $`1024`$.

<a id="s-0-3"></a>

## 0.3 What exists, and what is new

We searched for software that computes numerically *inside* finite-dimensional von Neumann
algebras. What we found falls into four groups, all useful and all different in purpose:

- **General matrix software** (NumPy/SciPy, QuTiP and other quantum toolboxes) works with one
  matrix at a time. Block structure, if present, has to be maintained by hand.
- **Symmetric tensor-network libraries** store tensors as block-sparse arrays over the charge
  sectors of a symmetry (abelian in YASTN [1], TeNPy and ITensor; TensorKit also handles
  non-abelian ones); YASTN has a PyTorch backend for GPUs and automatic differentiation.
  Technically they are the closest relatives. Their object is a tensor network to be contracted,
  not an operator algebra: as far as we found in their documentation, they do not provide the
  family of traces, states and densities, the centre and conditional expectations, Haar measures
  on $`\prod_cU(k_c)`$, or channels between sectors.
- **Numerical block diagonalisation** of matrix $`*`$-algebras [2, 3, 4] finds the decomposition
  $`\mathcal M\cong\bigoplus_cM_{k_c}\otimes1_{m_c}`$ of a given algebra — typically to reduce the size
  of semidefinite programs with symmetry. This is complementary: those algorithms find the blocks;
  this library computes in them. Its input step (module `decompose`, Section 4.3) uses the idea of
  [2, 3] that the eigenspaces of a random element of the commutant find the blocks; Section 4.3
  proves that in exact arithmetic this succeeds with probability one (Lemma 4.6, Proposition 4.7,
  Lemma 4.8), and Proposition 4.9 is the criterion that `Decomposition.verify()` checks numerically.
  Exact symbolic decompositions of group algebras exist in computer algebra (e.g. the GAP package
  Wedderga).
- **Noncommutative polynomial optimisation** (NCSOStools, ncpol2sdpa) bounds eigenvalues and traces
  of polynomials in noncommuting variables, over all matrix sizes at once, by semidefinite
  relaxations [5]. This answers a different question (a bound valid for all matrix sizes) by
  different means (convex relaxation).

*Lecture notes:* the randomised decomposition algorithm and its criterion — [Section 4.3](../lecture_notes/04_symmetries.md#s-4-3), [Proposition 4.7](../lecture_notes/04_symmetries.md#r-4-7), [Proposition 4.9](../lecture_notes/04_symmetries.md#r-4-9).

This survey is a description of the software as we found it, not a theorem. To our knowledge
there is no general library in which a finite-dimensional von Neumann algebra is the basic data
type for batched, differentiable numerical computation. The novelty is therefore not new
mathematics but **efficient computation**: the identification of §0.1, implemented so that the
structure of the algebra is preserved and the operations run as batched kernels (one call per
distinct block size for decompositions, split into one chunk per thread on a multi-core CPU above
the size threshold of §0.1). If you know of related software, we would be glad to hear about it.

<a id="s-0-4"></a>

## 0.4 How to read this course

| chapters | content |
|---|---|
| 1–3 | the mathematics: matrix algebras, the structure theorem, traces, states and entropy — with complete proofs |
| 4–8 | where these algebras come from and what one computes in them: symmetries, random matrices, channels, tensor products, dynamics |
| 9–12 | applications: spin chains, quantum computing, optimisation, equivariant machine learning |
| 13 | foundations: Jordan algebras, and why quantum mechanics uses *-algebras |
| appendices | complete proofs of the auxiliary results that the chapters cite by number (e.g. Theorem A.24) |

Each tutorial accompanies the chapter of the [lecture notes](../lecture_notes/00_overview.md) with the same number: it shows the
library code, checks the code numerically against results proved there, and names the standard tools to compare with.
Each section starts with a line *Lecture notes:* that links to the section and the results it uses. The code blocks of
a tutorial are meant to be run in order, in one Python session. Numbers that the lecture notes do not prove (timings, for instance) appear
only as *numerical observations* of the tutorial's own code.

## References

1. M. M. Rams, G. Wójtowicz, A. Sinha, J. Hasik, *YASTN: Yet another symmetric tensor networks; A
   Python library for Abelian symmetric tensor network calculations*, SciPost Phys. Codebases 52
   (2025).
2. K. Murota, Y. Kanno, M. Kojima, S. Kojima, *A numerical algorithm for block-diagonal
   decomposition of matrix $`*`$-algebras with application to semidefinite programming*, Japan J.
   Indust. Appl. Math. 27 (2010) 125–160 (preprint: Technical Report METR 2007-52, University of
   Tokyo, 2007).
3. T. Maehara, K. Murota, *A numerical algorithm for block-diagonal decomposition of matrix
   $`*`$-algebras with general irreducible components*, Japan J. Indust. Appl. Math. 27 (2010)
   263–293.
4. E. de Klerk, C. Dobre, D. V. Pasechnik, *Numerical block diagonalization of matrix
   $`*`$-algebras with application to semidefinite programming*, Math. Program. 129 (2011) 91–111.
5. I. Klep, J. Povh, *Semidefinite programming and sums of hermitian squares of noncommutative
   polynomials*, J. Pure Appl. Algebra 214 (2010) 740–749.

**Next:** [Tutorial 1](01_matrix_algebras.md): operators, functions of operators and batches in the library.
