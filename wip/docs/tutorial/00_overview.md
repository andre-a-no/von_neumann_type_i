**English** · [Español](es/00_overview.md) · [Français](fr/00_overview.md) · [Deutsch](de/00_overview.md) · [中文](zh/00_overview.md) · [日本語](ja/00_overview.md) · [Русский](ru/00_overview.md)

# 0. What this course is about

Von Neumann algebras are studied almost exclusively with pen and paper. In finite dimensions
every one of them is a direct sum of full matrix algebras (Chapter 2), so in principle anything
can be computed — but in practice an algebra is usually flattened into one large matrix, and
its structure (the sectors, the centre, the many traces) disappears from the computation. This
course and the library `torch_vn_algebra` take the opposite approach: **the algebra itself is
the data type**, and computations in it are carried out efficiently.

## 0.1 The observation

A finite-dimensional von Neumann algebra is, up to a unitary change of basis,
$$
\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c}
$$
(Theorem 2.6). As an abstract algebra it is $\bigoplus_cM_{k_c}(\mathbb C)$: an element is a list
of $C$ square blocks. This is exactly the shape of a tensor with a *sector axis*: an element of
$\mathcal M$ is an array of shape $(C,k,k)$ (blocks padded to $k=\max_ck_c$), and $B$ elements —
Monte Carlo samples, optimisation starts, time steps, disorder realisations — form an array of
shape $(B,C,k,k)$. Modern array libraries (PyTorch on CPUs and GPUs) provide batched routines —
matrix products, eigendecompositions, QR decompositions — for precisely such "batches of small
matrices", which neural networks also use.

Under this identification the operations of the algebra are batched array operations. The
correspondence is exact (each row of the table is an identity proved in the course); only the
numerical routines that carry the operations out work in floating point:

| in the algebra $\mathcal M$ | on the array of shape $(B,C,k,k)$ |
|---|---|
| product $AB$, sum, adjoint | batched matrix product, sum, conjugate transpose over $(B,C)$ |
| functional calculus $f(A)$ (Chapter 1) | batched eigendecomposition of the blocks |
| a trace $\sum_cw_c\operatorname{Tr}A_c$ (Chapter 3) | diagonal sum per block, weighted sum over the sector axis |
| the centre, central projections (Chapter 2) | one number per sector |
| conditional expectation onto the centre | block trace divided by $k_c$, broadcast back |
| Haar measure on the unitary group $\prod_cU(k_c)$ | batched QR decomposition of Gaussian blocks (Chapter 5) |
| states, channels, dynamics (Chapters 3, 6, 8) | the same operations, composed |

The proofs: by Theorem 2.6, $\mathcal M\cong\bigoplus_cM_{k_c}(\mathbb C)$ as $*$-algebras, and in a direct sum
products, sums and adjoints are taken block by block. If $A=\bigoplus_cA_c$ is normal and
$A_c=U_c\Lambda_cU_c^*$ (Theorem 1.3), then $U=\bigoplus_cU_c$ diagonalises $A$, so
$f(A)=\bigoplus_cf(A_c)$ by Theorem 1.5(a). The traces are described by Theorem 3.1, the centre by
Theorem 2.6, the conditional expectation onto the centre by Proposition 2.8, and the Haar measure of
$\prod_cU(k_c)$ by Corollary 5.4 and Theorem 5.3.

Two consequences make this more than bookkeeping.

- **Cost.** A spectral decomposition in $\mathcal M$ is computed block by block (by the remark after
  the table). The standard dense eigensolvers (LAPACK, which PyTorch calls) use a number of
  floating-point operations that grows like the cube of the matrix size; this is a description of
  those algorithms, not a result of this course. In this cost model a spectral decomposition in
  $\mathcal M$ costs $\sum_ck_c^3$ against $D^3$ for the embedding into $M_D(\mathbb C)$, $D=\sum_ck_c$; for $C$
  equal blocks of size $k$ the ratio is $(Ck)^3/(Ck^3)=C^2$. Products run on the padded tensor; the
  schoolbook product of two $k\times k$ matrices uses $k^3$ multiplications, so a product in $\mathcal M$
  uses $C\,k_{\max}^3$, which equals $\sum_ck_c^3$ when all blocks have the same size. An element is
  stored as $C\,k_{\max}^2$ numbers against $D^2$; for equal blocks the ratio is $C$. As for the
  implementation: products, sums and traces process all $B\cdot C$ blocks in one call;
  decompositions (eigenvalues, singular values, QR) run as one batched call per distinct block
  size — a single call when all sectors have the same size — and on a multi-core CPU, for batches
  above a size threshold, the library splits each such call into one chunk per thread (module
  `cost`). The aim is to keep a GPU, or all cores of a CPU, busy.
- **Correctness by construction.** The data of an element are its blocks, so everything that is
  computed stays in $\mathcal M$: random unitaries are Haar distributed on $\prod_cU(k_c)$
  (Corollary 5.4), a subgroup that a Haar unitary of $U(D)$ hits with probability zero for $C\ge2$
  (Section 5.1, after Corollary 5.4); a gradient step changes the blocks and so cannot leave the
  algebra; and the trace, the densities and the conditional expectations refer to the algebra and
  not to an accidental embedding.

## 0.2 A first computation

We compare the same computation done in $\mathcal M=\bigoplus_{c=1}^{64}M_{16}(\mathbb C)$ and in its
embedding into $M_{1024}(\mathbb C)$: the square root of a random positive element, and its trace.
`TypeIAlgebra(n_factors, k_factors)` takes two lists of block sizes: the nominal sizes $n_c$ and
the sizes $k_c\le n_c$ of the blocks actually used; in this course the two lists are always equal
(Section 1.5).

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

assert torch.allclose(torch.block_diag(*root_blocks[0]), root_big, atol=1e-10)   # same result
assert abs(A.sqrt().Tr_blunt().item() - root_big.trace().item()) < 1e-9
assert big.numel() == C * blocks.numel()                          # storage: a factor C = 64
print(f"in the algebra: {t_blocks * 1e3:.1f} ms, as one 1024 x 1024 matrix: {t_big * 1e3:.1f} ms")
```

The results agree to round-off, and the storage ratio is $C$, as counted above. *Numerical
observation* (timings depend on the hardware; the code prints them): on the 4-core CPU used to
write this chapter the algebra version took about 10 ms and the large-matrix version between 400
and 700 ms, a factor of several tens. The 64 blocks go to the eigensolver as one batched call,
which the library splits into one chunk per thread. The measured factor is far below the ratio
$C^2=4096$ of the cubic cost model: an operation count is a model of the running time, not the
running time, and on this machine the small eigenproblems ran at a much lower rate of operations
per second than the large one. In the cost model the ratio $C^2$ grows with the number of
sectors. With a batch of $B$ elements the algebra version stays one batched call; the
large-matrix version has to solve $B$ eigenproblems of size $1024$.

## 0.3 What exists, and what is new

We searched for software that computes numerically *inside* finite-dimensional von Neumann
algebras. What we found falls into four groups, all useful and all different in purpose:

- **General matrix software** (NumPy/SciPy, QuTiP and other quantum toolboxes) works with one
  matrix at a time. Block structure, if present, has to be maintained by hand.
- **Symmetric tensor-network libraries** store tensors as block-sparse arrays over charge sectors
  of an abelian symmetry, some with a PyTorch backend for GPUs and automatic differentiation
  (e.g. YASTN [1]; also TeNPy, ITensor, TensorKit). Technically they are the closest relatives. Their
  object is a tensor network to be contracted, not an operator algebra: as far as we found in
  their documentation, they do not provide the family of traces, states and densities, the centre
  and conditional expectations, Haar measures on $\prod_cU(k_c)$, or channels between sectors.
- **Numerical block diagonalisation** of matrix $*$-algebras [2, 3, 4] finds the decomposition
  $\mathcal M\cong\bigoplus_cM_{k_c}\otimes1_{m_c}$ of a given algebra — typically to reduce the size
  of semidefinite programs with symmetry. This is complementary: those algorithms find the blocks;
  this library computes in them. Its input step (`decompose`, Chapter 4) uses the idea of [2, 3]
  that the eigenspaces of a random element of the commutant find the blocks; Chapter 4 proves that
  this succeeds with probability one in exact arithmetic. Exact symbolic decompositions of group
  algebras exist in computer algebra (e.g. the GAP package Wedderga).
- **Noncommutative polynomial optimisation** (NCSOStools, ncpol2sdpa) bounds eigenvalues and traces
  of polynomials in noncommuting variables, over all matrix sizes at once, by semidefinite
  relaxations [5]. This answers a different
  question (a bound valid for all algebras) by different means (convex relaxation).

To our knowledge there is no general library in which a finite-dimensional von Neumann algebra
is the basic data type for batched, differentiable numerical computation. The novelty is therefore
not new mathematics but **efficient computation**: the identification of §0.1, implemented so that
the structure of the algebra is preserved and the operations run as batched kernels (one per
distinct block size for decompositions, split into one chunk per thread on a multi-core CPU). If you
know of related software, we would be glad to hear about it.

## 0.4 How to read this course

| chapters | content |
|---|---|
| 1–3 | the mathematics: matrix algebras, the structure theorem, traces, states and entropy — with complete proofs |
| 4–8 | where these algebras come from and what one computes in them: symmetries, random matrices, channels, tensor products, dynamics |
| 9–12 | applications: spin chains, quantum computing, optimisation, equivariant machine learning |
| 13 | foundations: Jordan algebras, and why quantum mechanics uses *-algebras |
| appendices | complete proofs of the auxiliary results that the chapters cite by number (e.g. Theorem A.24) |

Each chapter states results with proofs, shows the library code, checks the code numerically
against results proved in the course, and names the standard tools to compare with. Results of
other chapters and of the appendices are cited by number. Numbers that the course does not prove
(timings, for instance) appear only as *numerical observations* of the chapter's own code, and
the boxes headed *Physical meaning* are physical interpretation, not theorems of this course.

## References

1. M. M. Rams, G. Wójtowicz, A. Sinha, J. Hasik, *YASTN: Yet another symmetric tensor networks; A
   Python library for Abelian symmetric tensor network calculations*, SciPost Phys. Codebases 52
   (2025).
2. K. Murota, Y. Kanno, M. Kojima, S. Kojima, *A numerical algorithm for block-diagonal
   decomposition of matrix $*$-algebras*, Technical Report METR 2007-52, University of Tokyo (2007).
3. T. Maehara, K. Murota, *A numerical algorithm for block-diagonal decomposition of matrix
   $*$-algebras, Part II: general algorithm*, Optimization Online (2009).
4. E. de Klerk, C. Dobre, D. V. Pasechnik, *Numerical block diagonalization of matrix
   $*$-algebras with application to semidefinite programming*, Math. Program. 129 (2011) 91–111.
5. I. Klep, J. Povh, *Semidefinite programming and sums of hermitian squares of noncommutative
   polynomials*, J. Pure Appl. Algebra 214 (2010) 740–749.

**Next:** [Chapter 1](01_matrix_algebras.md).
