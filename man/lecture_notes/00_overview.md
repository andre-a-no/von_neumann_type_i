**English** · [Español](es/00_overview.md) · [Français](fr/00_overview.md) · [Deutsch](de/00_overview.md) · [中文](zh/00_overview.md) · [日本語](ja/00_overview.md) · [Русский](ru/00_overview.md)

# 0. What this course is about

Von Neumann algebras are mostly studied with pen and paper. In finite dimensions every one of
them is $*$-isomorphic to a direct sum of full matrix algebras (Theorem 2.6), so in principle
anything can be computed — but when such an algebra is handled numerically, it is often
flattened into one large matrix, and its structure (the sectors, the centre, the many traces)
disappears from the computation. This course and the library `torch_vn_algebra` take the
opposite approach: **the algebra itself is the data type**, and computations in it are carried
out efficiently.

<a id="s-0-1"></a>

## 0.1 The observation

A von Neumann algebra $\mathcal M$ on $\mathbb C^n$ is, up to a unitary change of basis,
$$
\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c}
$$
(Theorem 2.6). As an abstract algebra it is $\bigoplus_cM_{k_c}(\mathbb C)$: an element is a list
of $C$ square blocks. This is exactly the shape of a tensor with a *sector axis*: an element of
$\mathcal M$ is an array of shape $(C,k,k)$ (blocks padded with zeros to $k=\max_ck_c$; this is how
the library stores it, Section 2.5), and $B$ elements —
Monte Carlo samples, optimisation starts, time steps, disorder realisations — form an array of
shape $(B,C,k,k)$. Modern array libraries (PyTorch on CPUs and GPUs) provide batched routines —
matrix products, eigendecompositions, QR decompositions — for precisely such "batches of small
matrices", which neural networks also use.

*In the library:* an element of $\mathcal M$ is an operator of a `TypeIAlgebra`, stored as one tensor of shape $(B,C,k_{\max},k_{\max})$ — [tutorial, Section 1.5](../tutorial/01_matrix_algebras.md#s-1-5) and [tutorial, Section 2.5](../tutorial/02_von_neumann_algebras.md#s-2-5).

Under this identification the operations of the algebra are batched array operations. The
correspondence is exact (each of the first six rows of the table is an identity, proved below
from results of the course); only the numerical routines that carry the operations out work in
floating point:

| in the algebra $\mathcal M$ | on the array of shape $(B,C,k,k)$ |
|---|---|
| product $AB$, sum, adjoint | batched matrix product, sum, conjugate transpose over $(B,C)$ |
| functional calculus $f(A)$ (Chapter 1) | batched eigendecomposition of the blocks |
| a trace $\sum_cw_c\mathrm{Tr}\,A_c$ (Chapter 3) | diagonal sum per block, weighted sum over the sector axis |
| the centre, central projections (Chapter 2) | one number per sector |
| conditional expectation onto the centre | block trace divided by $k_c$, broadcast back |
| Haar measure on the unitary group $\prod_cU(k_c)$ | batched QR decomposition of Gaussian blocks, with the phase correction of Theorem 5.3 |
| states, channels, dynamics (Chapters 3, 6, 8) | built from the operations above |

*Proofs of the first six rows.* By Theorem 2.6 there is a unitary $V$ with
$\mathcal M=V\big(\bigoplus_cM_{k_c}\otimes1_{m_c}\big)V^*$, and
$\varphi\big(V(\bigoplus_cA_c\otimes1_{m_c})V^*\big)=\bigoplus_cA_c$ is a $*$-isomorphism of $\mathcal M$ onto
$\bigoplus_cM_{k_c}(\mathbb C)$ (well defined, as $A_c\otimes1_{m_c}$ determines $A_c$; multiplicative, as
$(A_c\otimes1)(B_c\otimes1)=A_cB_c\otimes1$). The library stores $\varphi(A)$ (Section 2.5).
(1) In $\bigoplus_cM_{k_c}$ sums, products and adjoints are taken block by block.
(2) Let $A\in\mathcal M$ be normal, $\varphi(A)=\bigoplus_cA_c$, and $f$ a function on $\mathrm{spec}\,A$. Each $A_c$ is
normal, since $\varphi(A)^*\varphi(A)=\varphi(A^*A)=\varphi(AA^*)=\varphi(A)\varphi(A)^*$ holds block by block. If
$A_cx=\lambda x$ with $x\neq0$, then $V(0\oplus\dots\oplus x\otimes e_1\oplus\dots\oplus0)$ ($x\otimes e_1$ in the
$c$-th summand) is an eigenvector of $A$ for $\lambda$, so $\mathrm{spec}\,A_c\subseteq\mathrm{spec}\,A$. Let $p$ be a
polynomial with $p=f$ on $\mathrm{spec}\,A$ (Lagrange interpolation, as in the proof of Theorem 1.5). For any
normal $X$ with $\mathrm{spec}\,X\subseteq\mathrm{spec}\,A$ and $X=U\Lambda U^*$, $f(X)=Uf(\Lambda)U^*=Up(\Lambda)U^*=p(X)$ by
Theorem 1.5(b); so $f(A)=p(A)$, $f(A_c)=p(A_c)$, and as $\varphi$ is a unital homomorphism,
$$
\varphi(f(A))=p(\varphi(A))=\bigoplus_cp(A_c)=\bigoplus_cf(A_c),\qquad f(A_c)=U_cf(\Lambda_c)U_c^*
$$
for any diagonalisation $A_c=U_c\Lambda_cU_c^*$ (Theorem 1.3): one eigendecomposition per block.
(3) The traces of $\bigoplus_cM_{k_c}$ are the functionals $\sum_cw_c\mathrm{Tr}\,A_c$ (Theorem 3.1).
(4) $Z(\mathcal M)$ is spanned by the minimal central projections $P_c$, and $\varphi(P_c)$ is $1_{k_c}$ in block
$c$ and $0$ elsewhere (Theorem 2.6), so $\varphi(Z(\mathcal M))=\{\bigoplus_cz_c1_{k_c}:z\in\mathbb C^C\}$.
(5) Proposition 2.8. (6) Corollary 5.4 (independent blocks, each Haar on $U(k_c)$) and Theorem 5.3
(each block from the QR decomposition of a Gaussian matrix). The last row is only a pointer: states
are represented by densities (Theorem 3.2), channels between sectors by Kraus operators between blocks
(Proposition 6.4), and the closed dynamics by $e^{-itH}=\bigoplus_ce^{-itH_c}$ (Proposition 8.2).

*In the library:* products and functions such as `sqrt` and `apply_function` — [tutorial, Section 1.5](../tutorial/01_matrix_algebras.md#s-1-5); the traces, e.g. `Tr_blunt` — [tutorial, Section 3.1](../tutorial/03_traces_and_states.md#s-3-1); the centre, `TypeIAlgebra.central` — [tutorial, Section 2.5](../tutorial/02_von_neumann_algebras.md#s-2-5); the conditional expectation onto the centre, `channels.center_expectation` — [tutorial, Section 2.6](../tutorial/02_von_neumann_algebras.md#s-2-6); Haar unitaries, `random_unitary_operator` — [tutorial, Section 5.1](../tutorial/05_random_matrices.md#s-5-1).

*Computation:* row 2 of the table on $\bigoplus_{c=1}^{64}M_{16}(\mathbb C)$ — [tutorial, Section 0.2](../tutorial/00_overview.md#s-0-2); what computing block by block costs — [tutorial, Section 0.1](../tutorial/00_overview.md#s-0-1).

<a id="s-0-2"></a>

## 0.2 A first computation

*Computation:* the square root and the trace of a random positive element of $\bigoplus_{c=1}^{64}M_{16}(\mathbb C)$, computed block by block and in $M_{1024}(\mathbb C)$ — [tutorial, Section 0.2](../tutorial/00_overview.md#s-0-2).

<a id="s-0-3"></a>

## 0.3 What exists, and what is new

*In the library:* related software, and what the library adds to it — [tutorial, Section 0.3](../tutorial/00_overview.md#s-0-3).

<a id="s-0-4"></a>

## 0.4 How to read this course

| chapters | content |
|---|---|
| 1–3 | the mathematics: matrix algebras, the structure theorem, traces, states and entropy — with complete proofs |
| 4–8 | where these algebras come from and what one computes in them: symmetries, random matrices, channels, tensor products, dynamics |
| 9–12 | applications: spin chains, quantum computing, optimisation, equivariant machine learning |
| 13 | foundations: Jordan algebras, and why quantum mechanics uses *-algebras |
| appendices | complete proofs of the auxiliary results that the chapters cite by number (e.g. Theorem A.24) |

Each of Chapters 1–13 states results with complete proofs. The prerequisites (linear algebra, real analysis, measure theory, elementary probability) are
listed in the [README](../lecture_notes/README.md#prerequisites); everything else is proved in the lecture notes, and
results of other chapters and of the appendices are cited by number. The boxes headed *Physical meaning* are
physical interpretation, not theorems of the lecture notes.

Each chapter has a twin in the [tutorial](../tutorial/00_overview.md), which shows the library code, checks it numerically against the results proved here and names the standard tools to compare with. The lines *Computation:* and *In the library:* after a result link to the place in the tutorial where it is computed or used.

**Next:** [Chapter 1](01_matrix_algebras.md).
