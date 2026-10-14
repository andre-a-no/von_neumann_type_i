**English** · [Español](es/12_equivariant_learning.md) · [Français](fr/12_equivariant_learning.md) · [Deutsch](de/12_equivariant_learning.md) · [中文](zh/12_equivariant_learning.md) · [日本語](ja/12_equivariant_learning.md) · [Русский](ru/12_equivariant_learning.md)

# 12. Equivariant machine learning

A layer of a neural network is **equivariant** if transforming the input by a symmetry and then applying the
layer gives the same result as applying the layer and then transforming the output: relabelling the vertices of a
graph relabels its features in the same way. For linear layers this is a commutation relation, so the equivariant
layers form a commutant, which by Chapter 2 is, up to a change of basis, the algebra $`\bigoplus_cM_{m_c}(\mathbb C)`$.
This chapter counts their parameters, builds and trains them with `torch_vn_algebra.decompose`, and uses the
algebra structure for unitary, positive and spectral layers that are equivariant by construction.

**You will learn** why the equivariant linear maps between two representations are block matrices of
sizes given by multiplicities; how to count them by orbits or by characters (DeepSets and the Bell numbers of
graph networks); why averaging over the group is the orthogonal projection onto the equivariant maps; how
functional calculus and unitary parametrisations act blockwise on the commutant; and how invariant vectors,
invariant quadratic forms and invariant polynomials of every degree (Molien's formula) are counted.

The layers are built and trained with the library in [Tutorial 12](../tutorial/12_equivariant_learning.md), which also checks the results of this chapter numerically.

Throughout, $`G`$ is a finite group, and a **representation** of $`G`$ on $`\mathbb C^n`$ is a homomorphism
$`\pi:G\to M_n(\mathbb C)`$ with unitary values (no restriction: invertible $`\pi(g)`$ are unitary for the inner product
$`\sum_g\langle\pi(g)x,\pi(g)y\rangle`$, Section 4.1). Its **character** is $`\chi_\pi(g)=\mathrm{Tr}\,\pi(g)`$. Irreducible and
equivalent representations and the **multiplicity** $`m_\rho(\pi)`$ of an irreducible $`\rho`$ in $`\pi`$ ($`0`$ if $`\rho`$
does not occur) are as in Section 2.4, which proves that $`m_\rho(\pi)`$ does not depend on the decomposition and
proves Schur's lemma; [Chapter 4](04_symmetries.md) treats both at length (Sections 4.1–4.2). On rectangular matrices
$`\langle S,T\rangle=\mathrm{Tr}(S^*T)`$ is the **Hilbert–Schmidt** inner product.

<a id="s-12-1"></a>

## 12.1 Equivariant linear maps

**Definition.** Let $`\pi_1,\pi_2`$ be representations of $`G`$ on $`\mathbb C^{n_1}`$, $`\mathbb C^{n_2}`$. A matrix
$`T\in M_{n_2\times n_1}(\mathbb C)`$ is **equivariant** (an **intertwiner**) if $`T\pi_1(g)=\pi_2(g)T`$ for all
$`g\in G`$. They form a vector space $`\mathrm{Hom}_G(\pi_1,\pi_2)`$; for $`\pi_1=\pi_2=\pi`$ it is the
commutant, $`\mathrm{Hom}_G(\pi,\pi)=\pi(G)'`$.

The commonest representations in machine learning permute coordinates. Let $`G`$ act on a finite set $`\Omega`$;
the **permutation representation** on $`\mathbb C^\Omega`$ is $`\pi(g)e_x=e_{gx}`$. Examples: $`S_n`$ permuting the
$`n`$ points of a set ($`\Omega=[n]=\{1,\dots,n\}`$), and $`S_n`$ relabelling the vertices of a graph, acting on
$`\Omega=[n]^2`$, i.e. on the adjacency matrix by $`X\mapsto P_\sigma XP_\sigma^T`$.

<a id="r-12-1"></a>

**Lemma 12.1 (orbits).** Let $`\pi_i`$ be the permutation representations of $`G`$ on $`\mathbb C^{\Omega_i}`$. Then
$`T\in\mathrm{Hom}_G(\pi_1,\pi_2)`$ iff $`T_{gy,gx}=T_{yx}`$ for all $`g`$, $`y\in\Omega_2`$, $`x\in\Omega_1`$. The
indicator matrices of the orbits of $`G`$ on $`\Omega_2\times\Omega_1`$ are a basis, and
$`\dim\mathrm{Hom}_G(\pi_1,\pi_2)`$ is the number of these orbits.

*Proof.* $`\pi_2(g)T\pi_1(g)^*e_x=\pi_2(g)Te_{g^{-1}x}=\sum_yT_{y,g^{-1}x}e_{gy}`$, so the $`(gy,x)`$ entry of
$`\pi_2(g)T\pi_1(g)^*`$ is $`T_{y,g^{-1}x}`$. Hence $`T`$ is equivariant, i.e. $`\pi_2(g)T\pi_1(g)^*=T`$ for all $`g`$, iff
$`T_{gy,x}=T_{y,g^{-1}x}`$ for all $`g,y,x`$, which (replace $`x`$ by $`gx`$) is the stated condition: $`T`$ is constant
on orbits. The indicator matrices of the orbits have disjoint supports, so they are linearly independent, and
every function constant on orbits is a combination of them. $`\square`$

*Computation:* Lemma 12.1 for $`S_4`$ on vertex and edge features, by counting orbits — [tutorial, Section 12.1](../tutorial/12_equivariant_learning.md#s-12-1).

<a id="r-12-2"></a>

**Corollary 12.2.** (a) (*DeepSets.*) For $`S_n`$ on $`\mathbb C^n`$, $`n\ge2`$, the commutant is
$`\mathrm{span}\{1,J\}`$, $`J`$ the all-ones matrix, and $`\mathbb C^n=\mathbb C\mathbf 1\oplus\mathbf 1^\perp`$ is
a sum of two inequivalent irreducible representations, the trivial one and the **standard** one of dimension
$`n-1`$. (b) (*Graph networks.*) The equivariant maps from $`\mathbb C^{[n]^p}`$ to $`\mathbb C^{[n]^q}`$ (tensors
of order $`p`$ and $`q`$ on $`n`$ vertices) are counted by the partitions of a set of $`p+q`$ elements into at most $`n`$
blocks; for $`n\ge p+q`$ their number is the **Bell number** $`B_{p+q}`$, the number of all partitions of a set of
$`p+q`$ elements ($`B_2=2`$, $`B_3=5`$, $`B_4=15`$).

*Proof.* (b) Two tuples $`x,y\in[n]^{p+q}`$ lie in the same orbit iff they have the same **equality pattern**,
the partition of the positions $`\{1,\dots,p+q\}`$ by $`i\sim j\Leftrightarrow x_i=x_j`$. The pattern is clearly
invariant. Conversely, if $`x,y`$ have the same pattern, $`x_i\mapsto y_i`$ is a well-defined injection of
$`\{x_i\}`$ into $`[n]`$, which extends to a permutation $`\sigma`$ with $`\sigma x=y`$. Every partition with at most
$`n`$ blocks is a pattern (give the blocks distinct values). Lemma 12.1 concludes. For $`n\ge p+q`$ the bound on the number
of blocks is void. Counting by the number of blocks: $`B_2=1+1`$; $`B_3=1+3+1`$ (a two-block partition of a $`3`$-set is
fixed by its singleton); $`B_4=1+7+6+1`$ (two blocks: $`4`$ of sizes $`3+1`$ and $`3`$ of sizes $`2+2`$; three blocks: the
$`\binom42=6`$ choices of the pair). (a) For $`p=q=1`$ there are two
orbits on $`[n]^2`$, the diagonal and its complement, with indicators $`1`$ and $`J-1`$. So
$`\dim\pi(S_n)'=2=\sum_cm_c^2`$ (Theorem 2.6, Proposition 2.7), which forces two sectors with $`m_c=1`$, i.e. two
inequivalent irreducible summands, and $`\pi(S_n)'`$ is spanned by their projections $`P_1,P_2`$. The projection onto the
invariant line $`\mathbb C\mathbf 1`$ lies in $`\pi(S_n)'`$ (Lemma 2.1), so it is $`P_1`$ or $`P_2`$: the line, which
carries the trivial representation, is one sector, and $`\mathbf 1^\perp`$ is the other. $`\square`$

The DeepSets layer $`x\mapsto\lambda x+\gamma\,(\mathbf 1^Tx)\mathbf 1`$ (Zaheer et al. 2017) is (a); the bases of
Maron et al. (2019) for graph networks are (b).

*Computation:* Corollary 12.2 for $`n=4`$ ($`B_2=2`$, $`B_3=5`$, $`B_4=15`$), and the decomposition of the vertex and edge features of $`S_4`$ by `dc.from_group` (edge features: sectors of dimensions $`k=(1,2,3,3)`$ with multiplicities $`m=(2,1,1,3)`$) — [tutorial, Section 12.1](../tutorial/12_equivariant_learning.md#s-12-1).

For two representations, everything reduces to one commutant:

<a id="r-12-3"></a>

**Lemma 12.3 (direct sum).** Let $`\pi=\pi_1\oplus\pi_2`$ on $`\mathbb C^{n_1}\oplus\mathbb C^{n_2}`$, with the
isometric inclusions $`Q_i`$ and projections $`P_i=Q_iQ_i^*`$. Then $`P_1,P_2\in\pi(G)'`$,
$`\mathrm{Hom}_G(\pi_1,\pi_2)=Q_2^*\,\pi(G)'\,Q_1`$, and $`T\mapsto Q_2TQ_1^*`$ is a linear bijection of
$`\mathrm{Hom}_G(\pi_1,\pi_2)`$ onto $`P_2\,\pi(G)'P_1`$.

*Proof.* $`\pi(g)Q_i=Q_i\pi_i(g)`$ and $`Q_i^*\pi(g)=\pi_i(g)Q_i^*`$, because $`\pi(g)`$ is block diagonal; hence
$`P_i\pi(g)=Q_i\pi_i(g)Q_i^*=\pi(g)P_i`$. If $`X\in\pi(G)'`$, then
$`Q_2^*XQ_1\pi_1(g)=Q_2^*X\pi(g)Q_1=Q_2^*\pi(g)XQ_1=\pi_2(g)Q_2^*XQ_1`$. If $`T`$ is equivariant, then
$`\pi(g)Q_2TQ_1^*=Q_2\pi_2(g)TQ_1^*=Q_2T\pi_1(g)Q_1^*=Q_2TQ_1^*\pi(g)`$, so $`Q_2TQ_1^*\in\pi(G)'`$, it equals
$`P_2(Q_2TQ_1^*)P_1`$, and $`Q_2^*(Q_2TQ_1^*)Q_1=T`$. Conversely every $`Y=P_2XP_1`$ with $`X\in\pi(G)'`$ lies in
$`\pi(G)'`$ and $`Y=Q_2(Q_2^*YQ_1)Q_1^*`$. $`\square`$

*Computation:* the corners $`Q_2^*\,\pi(G)'\,Q_1`$ for vertex and edge features — [tutorial, Section 12.2](../tutorial/12_equivariant_learning.md#s-12-2).

<a id="r-12-4"></a>

**Theorem 12.4 (equivariant maps are block matrices).** Let $`\rho_1,\dots,\rho_C`$ be the pairwise inequivalent
irreducible representations occurring in $`\pi_1\oplus\pi_2`$, $`k_c=\dim\rho_c`$, $`a_c=m_{\rho_c}(\pi_1)`$,
$`b_c=m_{\rho_c}(\pi_2)`$. There are unitaries $`V_1:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{a_c}\to\mathbb C^{n_1}`$
and $`V_2:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{b_c}\to\mathbb C^{n_2}`$ with
$`V_i^*\pi_i(g)V_i=\bigoplus_c\rho_c(g)\otimes1`$, such that
```math
\mathrm{Hom}_G(\pi_1,\pi_2)=\Big\{V_2\Big(\bigoplus_{c=1}^C1_{k_c}\otimes y_c\Big)V_1^*:\ y_c\in M_{b_c\times a_c}(\mathbb C)\Big\},
```
and $`(y_c)_c\mapsto T`$ is injective. Hence $`\mathrm{Hom}_G(\pi_1,\pi_2)\cong\bigoplus_cM_{b_c\times a_c}(\mathbb C)`$,
```math
\dim\mathrm{Hom}_G(\pi_1,\pi_2)=\sum_\rho m_\rho(\pi_1)\,m_\rho(\pi_2),
```
the sum running over all irreducible representations up to equivalence; for $`\pi_1=\pi_2`$ this is
$`\pi(G)'\cong\bigoplus_cM_{m_c}`$ (Proposition 2.7), and $`\mathrm{Hom}_G(\pi_1,\pi_2)=0`$ iff no
irreducible representation occurs in both.

*Proof.* Apply Theorem 2.6 and Proposition 2.7 to $`\pi=\pi_1\oplus\pi_2`$: there is a unitary $`V`$ with
$`V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}`$ and $`V^*\pi(G)'V=\bigoplus_c1_{k_c}\otimes M_{m_c}`$. By the uniqueness
of multiplicities (Section 2.4) every irreducible representation occurring in $`\pi`$ is equivalent to exactly one
$`\rho_c`$, so the $`\rho_c`$ are representatives of the classes in the statement. Any other representatives are
$`\rho_c'=w_c\rho_cw_c^*`$ with unitaries $`w_c`$ (Lemma 4.2(c)), and replacing $`V`$ by $`V\bigoplus_c(w_c^*\otimes1)`$ turns
the first normal form into $`\bigoplus_c\rho'_c(g)\otimes1_{m_c}`$ and leaves the second unchanged, as $`w_c^*\otimes1`$
commutes with $`1\otimes M_{m_c}`$; so we may use the $`\rho_c`$ of Proposition 2.7. Since
$`X\mapsto V^*XV`$ is a $`*`$-isomorphism and $`P_1\in\pi(G)'`$ (Lemma 12.3), $`V^*P_1V=\bigoplus_c1\otimes p_c`$ with
projections $`p_c`$. Choose unitaries $`u_c`$ with $`u_c^*p_cu_c=1_{a_c}\oplus0_{b_c}`$ (Theorem 1.3), where now
$`a_c=\mathrm{rank}\,p_c`$, $`b_c=m_c-a_c`$, and replace $`V`$ by $`V\bigoplus_c(1\otimes u_c)`$; this changes neither
normal form, as $`1\otimes u_c`$ commutes with $`\rho_c(g)\otimes1`$. Let $`\iota_c:\mathbb C^{a_c}\to\mathbb C^{m_c}`$ and
$`\kappa_c:\mathbb C^{b_c}\to\mathbb C^{m_c}`$ be the inclusions of the first $`a_c`$ and the last $`b_c`$ coordinates, and
```math
V_1=Q_1^*V\bigoplus_c(1\otimes\iota_c),\qquad V_2=Q_2^*V\bigoplus_c(1\otimes\kappa_c).
```
$`W_1=V\bigoplus_c(1\otimes\iota_c)`$ is an isometry with range $`\mathrm{range}\,P_1=\mathrm{range}\,Q_1`$, so
$`V_1^*V_1=W_1^*P_1W_1=1`$, and $`V_1`$ is unitary because $`\sum_ck_ca_c=\mathrm{rank}\,P_1=n_1`$. Moreover
$`\pi_1(g)V_1=Q_1^*\pi(g)W_1=Q_1^*V\big(\bigoplus_c\rho_c(g)\otimes1_{m_c}\big)\bigoplus_c(1\otimes\iota_c) =V_1\bigoplus_c\rho_c(g)\otimes1_{a_c}`$. So $`\pi_1\cong\bigoplus_c\rho_c\otimes1_{a_c}`$ with inequivalent
irreducible $`\rho_c`$, and $`a_c=m_{\rho_c}(\pi_1)`$ by the uniqueness of multiplicities (Section 2.4); likewise for
$`V_2`$ and $`b_c`$. Finally, for $`X=V\big(\bigoplus_c1\otimes x_c\big)V^*\in\pi(G)'`$,
$`P_2XP_1=V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*`$, where $`y_c=\kappa_c^*x_c\iota_c`$ is the lower left
$`b_c\times a_c`$ corner of $`x_c`$, an arbitrary matrix. By Lemma 12.3,
$`\mathrm{Hom}_G(\pi_1,\pi_2)=Q_2^*P_2\,\pi(G)'P_1Q_1`$ consists of the matrices
$`Q_2^*V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*Q_1=V_2\big(\bigoplus_c1\otimes y_c\big)V_1^*`$, and
$`(y_c)\mapsto T`$ is injective since $`V_1,V_2`$ are unitary. An irreducible $`\rho`$ that does not occur in
$`\pi_1\oplus\pi_2`$ has $`m_\rho(\pi_1)m_\rho(\pi_2)=0`$. $`\square`$

*In the library:* the proof of Theorem 12.4 as an algorithm, the layers `EquivariantLinear` and `EquivariantMap` built on `dc.from_group`, `D.commutant_algebra()` and `D.from_commutant`; for vertex to edge features of $`S_4`$, $`\dim\mathrm{Hom}_G=1\cdot2+1\cdot3=5=B_3`$ — [tutorial, Section 12.2](../tutorial/12_equivariant_learning.md#s-12-2).

> **Physical meaning** (physical interpretation, not a theorem of this course). Theorem 12.4 is the linear-algebra
> core of **selection rules** and the scalar case of the Wigner–Eckart theorem: an operator that commutes with the
> symmetry can only connect states of the same symmetry type, and within a type it acts as $`1_{k_c}\otimes y_c`$,
> i.e. by a "reduced matrix element" $`y_c`$ that is blind to the position inside the multiplet. (Operators that
> themselves transform by a non-trivial representation, such as vectors, need in addition the Clebsch–Gordan
> coefficients.) A symmetric Hamiltonian is the case $`\pi_1=\pi_2`$; its degeneracies $`k_c`$ are the dimensions of
> the irreducible representations.

<a id="s-12-2"></a>

## 12.2 An equivariant layer with the library

*In the library:* equivariant layers parametrised by the blocks of Theorem 12.4 — [tutorial, Section 12.2](../tutorial/12_equivariant_learning.md#s-12-2).

The second way to make a layer equivariant is to average it over the group. The **twirl** is
```math
\mathcal P(W)=\frac1{|G|}\sum_{g\in G}\pi_2(g)\,W\,\pi_1(g)^*,\qquad W\in M_{n_2\times n_1}(\mathbb C).
```

<a id="r-12-5"></a>

**Theorem 12.5 (the twirl is the orthogonal projection).** $`\mathcal P`$ is the Hilbert–Schmidt orthogonal
projection of $`M_{n_2\times n_1}(\mathbb C)`$ onto $`\mathrm{Hom}_G(\pi_1,\pi_2)`$, and
```math
\dim\mathrm{Hom}_G(\pi_1,\pi_2)=\frac1{|G|}\sum_{g\in G}\overline{\chi_{\pi_1}(g)}\,\chi_{\pi_2}(g).
```

*Proof.* For $`h\in G`$, $`\pi_2(h)\mathcal P(W)\pi_1(h)^*=\frac1{|G|}\sum_g\pi_2(hg)W\pi_1(hg)^*=\mathcal P(W)`$, as
$`g\mapsto hg`$ is a bijection of $`G`$; so $`\mathcal P(W)`$ is equivariant. If $`W`$ is equivariant, every term equals
$`W`$, so $`\mathcal P(W)=W`$. Hence $`\mathcal P^2=\mathcal P`$ with range $`\mathrm{Hom}_G(\pi_1,\pi_2)`$. It is
self-adjoint: $`\langle S,\mathcal P(W)\rangle=\frac1{|G|}\sum_g\mathrm{Tr}\big((\pi_2(g)^*S\pi_1(g))^*W\big)`$
$`=\langle\mathcal P(S),W\rangle`$, where we substituted $`g\mapsto g^{-1}`$ and used $`\pi_i(g^{-1})=\pi_i(g)^*`$. A
self-adjoint idempotent is the orthogonal projection onto its range. Its trace is its rank. The linear map
$`W\mapsto AWB`$ on $`M_{n_2\times n_1}`$ has trace $`\mathrm{Tr}\,A\mathrm{Tr}\,B`$ (the coefficient of the
matrix unit $`E_{ij}`$ in $`AE_{ij}B`$ is $`A_{ii}B_{jj}`$); with $`A=\pi_2(g)`$, $`B=\pi_1(g)^*`$ this gives the formula. $`\square`$

With Theorem 12.4, $`\sum_\rho m_\rho(\pi_1)m_\rho(\pi_2)=|G|^{-1}\sum_g\overline{\chi_{\pi_1}(g)}\chi_{\pi_2}(g)`$;
for $`\pi_1=\pi_2`$ this is the character formula for $`\dim\pi(G)'`$ of Section 4.2. The parametrisation $`W\mapsto\mathcal P(W)`$ reaches every
equivariant map, $`\mathcal P(W)`$ is the equivariant map closest to $`W`$, and as $`\mathcal P`$ is self-adjoint, the
gradient with respect to $`W`$ is $`\mathcal P`$ of the gradient with respect to $`T=\mathcal P(W)`$. But
$`n_1n_2-\dim\mathrm{Hom}_G`$ of its parameters are wasted on the kernel.

*Computation:* Theorem 12.5 on the edge features and the dimension formula from vertex to edge features — [tutorial, Section 12.2](../tutorial/12_equivariant_learning.md#s-12-2). *In the library:* the twirl `dc.twirl`, the same projection computed from the decomposition, `D.conditional_expectation(W, 'commutant')`, and the character formula `dc.commutant_dimension` — [tutorial, Section 12.2](../tutorial/12_equivariant_learning.md#s-12-2).

The layers above are linear. A network alternates them with nonlinearities, and these must be equivariant too.
For a permutation representation, a function applied to each coordinate is equivariant, since
$`(\pi(g)x)_{gy}=x_y`$; for other representations it is not, and one uses nonlinearities that act on norms or on
the coefficients of a basis that $`G`$ permutes (Exercise 3).

> **Physical meaning** (physical interpretation, not a theorem of this course). The twirl is the conditional
> expectation onto the commutant (Exercise 3 of Chapter 2): averaging a coupling over all symmetry operations
> keeps exactly its symmetric part. "Symmetrising" a Hamiltonian or an interaction by summing over the symmetry
> group is the same projection.

<a id="s-12-3"></a>

## 12.3 Constrained and spectral layers

The weights $`w=(w_c)_c`$ of an equivariant layer $`x\mapsto Tx`$, $`T\in\pi(G)'`$, live in the $`*`$-algebra $`\bigoplus_cM_{m_c}`$
(Theorem 12.4), so every construction of Chapters 1 and 11 applies to them. Write $`\Phi(w)=V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*`$ for $`w\in\bigoplus_cM_{m_c}`$.

<a id="r-12-6"></a>

**Proposition 12.6.** (a) $`\Phi`$ is a $`*`$-isomorphism of $`\bigoplus_cM_{m_c}`$ onto $`\pi(G)'`$, and $`\Phi(w)`$ is
self-adjoint, unitary or positive iff every $`w_c`$ is. (b) For normal $`w`$ and a function $`f`$ on
$`\bigcup_c\mathrm{spec}\,w_c`$, $`f(\Phi(w))=\Phi\big((f(w_c))_c\big)`$; in particular $`f(X)\in\pi(G)'`$ for every
normal $`X\in\pi(G)'`$. (c) The equivariant unitaries are the $`\Phi(u_0e^{s-s^*})`$, $`s`$ ranging over
$`\bigoplus_cM_{m_c}`$, for any fixed unitary $`u_0\in\bigoplus_cM_{m_c}`$; in particular they are the exponentials $`e^S`$
of the skew-adjoint $`S\in\pi(G)'`$.

*Proof.* (a) $`\Phi`$ is linear, unital, multiplicative and $`*`$-preserving because $`(1\otimes x)(1\otimes y)=1\otimes xy`$
and $`(1\otimes x)^*=1\otimes x^*`$; it is injective, and onto $`\pi(G)'`$ by Theorem 2.6. Hence $`\Phi(w)^*=\Phi(w)`$ iff
$`w^*=w`$, and $`\Phi(w)^*\Phi(w)=1=\Phi(w)\Phi(w)^*`$ iff $`w^*w=1=ww^*`$. If all $`w_c=t_c^*t_c`$, then
$`\Phi(w)=\Phi(t)^*\Phi(t)\ge0`$ (Corollary 1.7). Conversely, if $`\Phi(w)\ge0`$, then for $`v\in\mathbb C^{m_c}`$ and
the vector $`\xi=V(e_1\otimes v)`$ in sector $`c`$, $`0\le\langle\xi,\Phi(w)\xi\rangle=\langle v,w_cv\rangle`$.
(b) Let $`w_c=u_c\Lambda_cu_c^*`$ (Theorem 1.3). Then $`\Phi(w)=U\big(\bigoplus_c1\otimes\Lambda_c\big)U^*`$ with the
unitary $`U=V\bigoplus_c(1\otimes u_c)`$, a diagonalisation with the eigenvalues of $`w_c`$ repeated $`k_c`$ times; by
definition and Theorem 1.5(a), $`f(\Phi(w))=U\big(\bigoplus_c1\otimes f(\Lambda_c)\big)U^*=\Phi\big((f(w_c))_c\big)`$.
(c) By (a), the equivariant unitaries are $`\Phi(u)`$ with $`u`$ unitary in $`\bigoplus_cM_{m_c}`$, and these are the
$`u_0e^{s-s^*}`$ by Lemma 11.1. With $`u_0=1`$, (b) gives $`\Phi(e^{s-s^*})=e^{\Phi(s-s^*)}`$, and $`\Phi(s-s^*)`$ runs over
the skew-adjoint elements of $`\pi(G)'`$. $`\square`$

*In the library:* unitary and positive equivariant layers with `optimize.UnitaryParam` and `PositiveParam` on `D.commutant_algebra()`, spectral layers with `apply_function`; training an equivariant unitary layer, and (b) for $`e^{-tL}`$ with its derivative in $`t`$ — [tutorial, Section 12.3](../tutorial/12_equivariant_learning.md#s-12-3).

> **Physical meaning** (physical interpretation, not a theorem of this course). A symmetric Hamiltonian
> $`H\in\pi(G)'`$ generates a symmetric time evolution $`e^{-itH}\in\pi(G)'`$ (Chapter 8), which never mixes sectors
> and acts on the multiplicity spaces only. An equivariant unitary layer is such an evolution, learned instead of
> given.

<a id="s-12-4"></a>

## 12.4 Invariants

A vector $`x`$ is **invariant** if $`\pi(g)x=x`$ for all $`g`$; a function $`q`$ on $`\mathbb C^n`$ is invariant if
$`q(\pi(g)x)=q(x)`$. Invariant features are what an invariant network (a graph classifier, say) can output.

<a id="r-12-7"></a>

**Proposition 12.7.** (a) The invariant vectors are the range of the orthogonal projection
$`P_0=|G|^{-1}\sum_g\pi(g)`$; their dimension is $`m_{\mathrm{triv}}(\pi)=|G|^{-1}\sum_g\chi_\pi(g)`$, and if the trivial
representation occurs, as $`\rho_{c_0}`$, then $`P_0`$ is the minimal central projection $`P_{c_0}`$ of $`\pi(G)''`$.
(b) A Hermitian form $`q_A(x)=\langle x,Ax\rangle`$, $`A=A^*`$, is invariant iff $`A\in\pi(G)'`$; the invariant Hermitian
forms are a real vector space of dimension $`\sum_cm_c^2`$.
(c) The homogeneous invariant polynomials of degree $`2`$ in the coordinates $`x_1,\dots,x_n`$ form a space of dimension
```math
\frac1{|G|}\sum_{g\in G}\frac{\chi_\pi(g)^2+\chi_\pi(g^2)}2 .
```
For a permutation representation it is the number of orbits of $`G`$ on the multisets $`\{x,y\}`$ of two elements of $`\Omega`$.

*Proof.* (a) Invariant vectors are $`\mathrm{Hom}_G(\mathrm{triv},\pi)`$, with $`\mathrm{triv}`$ the trivial
representation on $`\mathbb C=M_{1\times1}`$; for it the twirl is $`\mathcal P(x)=P_0x`$, so Theorem 12.5 gives the
projection and the dimension, and Theorem 12.4 gives $`m_{\mathrm{triv}}(\pi)`$. In the basis of Theorem 2.6, write an
invariant vector as $`\xi=\sum_c\xi_c`$ with $`\xi_c=\sum_jv_{c,j}\otimes e_j\in\mathbb C^{k_c}\otimes\mathbb C^{m_c}`$.
Since $`\pi(g)`$ preserves the sectors, $`\rho_c(g)v_{c,j}=v_{c,j}`$ for all $`g`$. A non-zero $`v_{c,j}`$ spans an invariant
line carrying the trivial representation, so then $`\rho_c`$, being irreducible, is trivial and $`c=c_0`$. Hence the
invariant vectors are the vectors of sector $`c_0`$, where $`\pi(g)`$ acts as $`1`$: the range of $`P_{c_0}`$, which is
therefore $`P_0`$.
(b) $`q_A(\pi(g)x)=\langle x,\pi(g)^*A\pi(g)x\rangle`$, and a matrix $`D`$ with $`\langle x,Dx\rangle=0`$ for all $`x`$ vanishes
(polarisation, proof of Corollary 1.4). So $`q_A`$ is invariant iff $`\pi(g)^*A\pi(g)=A`$ for all $`g`$, i.e. $`A\in\pi(G)'`$.
By Proposition 12.6(a) the self-adjoint elements of $`\pi(G)'`$ are the $`\Phi(w)`$ with Hermitian $`w_c`$, a real space of
dimension $`\sum_cm_c^2`$.
(c) A quadratic polynomial $`p(x)=x^TBx`$ is determined by its symmetric matrix $`B`$. Put $`(g\cdot p)(x)=p(\pi(g)^{-1}x)`$;
then $`g\cdot(h\cdot p)=(gh)\cdot p`$, and $`g\cdot p`$ has the matrix $`\sigma(g)(B)=A^TBA`$, $`A=\pi(g)^{-1}`$. The
average $`M=|G|^{-1}\sum_g\sigma(g)`$ satisfies $`\sigma(h)M=M`$ and $`MB=B`$ for invariant $`B`$, so (as in the proof of
Theorem 12.5) it is idempotent with range the invariant polynomials, and their number is
$`\mathrm{Tr}\,M=|G|^{-1}\sum_g\mathrm{Tr}\,\sigma(g)`$ (an idempotent is $`1\oplus0`$ in a basis adapted to range
$`\oplus`$ kernel, so its trace is its rank). $`A`$ is unitary, hence $`A=R\Lambda R^{-1}`$ with
$`\Lambda=\mathrm{diag}(\lambda)`$ (Theorem 1.3). The symmetric matrices
$`B_{ij}=(R^{-1})^T(e_ie_j^T+e_je_i^T)R^{-1}`$, $`i\le j`$, form a basis, and
$`A^TB_{ij}A=(R^{-1})^T\Lambda(e_ie_j^T+e_je_i^T)\Lambda R^{-1}=\lambda_i\lambda_jB_{ij}`$. So
$`\mathrm{Tr}\,\sigma(g)=\sum_{i\le j}\lambda_i\lambda_j=\frac12\big((\mathrm{Tr}\,A)^2+\mathrm{Tr}\,A^2\big) =\frac12\big(\chi_\pi(g^{-1})^2+\chi_\pi(g^{-2})\big)`$, and summing over $`g^{-1}`$ instead of $`g`$ gives the formula.
For a permutation representation, $`G`$ permutes the monomials $`x_ax_b`$, so a polynomial is invariant iff its
coefficients are constant on the orbits of monomials. $`\square`$

*Computation:* (a), (b) and the two counts of (c) for the edge features of $`S_4`$ — [tutorial, Section 12.4](../tutorial/12_equivariant_learning.md#s-12-4).

For graphs, the invariant edge features are spanned by the indicators of the diagonal and of its complement
($`B_2=2`$), and the invariant linear readouts are the sums of the diagonal and of the off-diagonal entries. On vertex features $`\mathbb C^4`$, (c) counts the orbits of $`S_4`$ on the multisets of two vertices: there are two.

On $`\mathbb C^4`$ they are $`\sum_ix_i^2`$ and $`\sum_{i<j}x_ix_j`$. The proof of (c) works in every degree:

<a id="r-12-8"></a>

**Proposition 12.8 (Molien's formula).** Let $`\mathbb C[x]_d`$ be the space of homogeneous polynomials of degree $`d`$ in
$`x_1,\dots,x_n`$ and $`\mathbb C[x]^G_d`$ the subspace of invariant ones, $`p(\pi(g)x)=p(x)`$ for all $`g`$. For $`|t|<1`$,
```math
\sum_{d\ge0}\dim\mathbb C[x]^G_d\,t^d=\frac1{|G|}\sum_{g\in G}\frac1{\det\big(1-t\pi(g)\big)} .
```
For a permutation representation, $`\dim\mathbb C[x]^G_d`$ is the number of orbits of $`G`$ on the monomials of degree $`d`$
(the multisets of $`d`$ elements of $`\Omega`$). The coefficient of $`t^2`$ is Proposition 12.7(c).

*Proof.* For an invertible matrix $`M`$ let $`\phi_M(p)(x)=p(Mx)`$ (substitution). $`\phi_M`$ is linear, maps $`\mathbb C[x]_d`$
into itself, and $`\phi_N\circ\phi_M=\phi_{MN}`$, so $`\phi_M`$ is bijective with inverse $`\phi_{M^{-1}}`$. Hence
$`\sigma(g)=\phi_{\pi(g)^{-1}}`$ satisfies $`\sigma(g)\sigma(h)=\sigma(gh)`$, and $`p`$ is invariant iff $`\sigma(g)p=p`$ for all
$`g`$. As in the proof of Proposition 12.7(c), $`M_d=|G|^{-1}\sum_g\sigma(g)`$ on $`\mathbb C[x]_d`$ is idempotent with range
$`\mathbb C[x]^G_d`$, so $`\dim\mathbb C[x]^G_d=\mathrm{Tr}\,M_d=|G|^{-1}\sum_g\mathrm{Tr}\,\sigma_d(g)`$, $`\sigma_d(g)`$ the
restriction of $`\sigma(g)`$ to $`\mathbb C[x]_d`$. Fix $`g`$ and write the unitary $`A=\pi(g)^{-1}`$ as $`A=R\Lambda R^{-1}`$,
$`\Lambda=\mathrm{diag}(\lambda_1,\dots,\lambda_n)`$, $`|\lambda_i|=1`$ (Theorem 1.3). With the linear forms
$`\ell_i(x)=(R^{-1}x)_i`$, the products $`\ell^\alpha=\prod_i\ell_i^{\alpha_i}=\phi_{R^{-1}}(x^\alpha)`$, $`|\alpha|=d`$, are a basis of
$`\mathbb C[x]_d`$ (the image of the monomial basis under the bijection $`\phi_{R^{-1}}`$), and since $`R^{-1}A=\Lambda R^{-1}`$,
```math
(\sigma(g)\ell^\alpha)(x)=\prod_i\big((R^{-1}Ax)_i\big)^{\alpha_i}=\prod_i\big(\lambda_i\ell_i(x)\big)^{\alpha_i}=\lambda^\alpha\ell^\alpha(x).
```
So $`\mathrm{Tr}\,\sigma_d(g)=\sum_{|\alpha|=d}\lambda^\alpha`$. For $`|t|<1`$ the product of the $`n`$ absolutely convergent
series $`(1-t\lambda_i)^{-1}=\sum_{a\ge0}t^a\lambda_i^a`$ may be expanded and regrouped by degree:
```math
\frac1{\det(1-tA)}=\prod_{i=1}^n\frac1{1-t\lambda_i}=\sum_\alpha t^{|\alpha|}\lambda^\alpha=\sum_{d\ge0}t^d\mathrm{Tr}\,\sigma_d(g).
```
Averaging over $`g`$ and replacing $`g`$ by $`g^{-1}`$ (so that $`A=\pi(g)`$) gives the formula; both sides converge, as
$`\dim\mathbb C[x]^G_d\le\dim\mathbb C[x]_d=\binom{n+d-1}d`$. For a permutation representation,
$`(\pi(g)^{-1}x)_y=x_{gy}`$, so $`\sigma(g)`$ maps the monomial $`\prod_jx_{y_j}`$ to the monomial $`\prod_jx_{gy_j}`$: $`G`$
permutes the monomials of degree $`d`$, a polynomial is invariant iff its coefficients are constant on the orbits, and
the orbit sums are a basis of $`\mathbb C[x]^G_d`$. $`\square`$

(This is Molien's theorem of 1897; see B. Sturmfels, *Algorithms in Invariant Theory*, Springer 2008, Section 2.2,
for further reading.) For $`S_4`$ on $`\mathbb C^4`$ the counts $`1,1,2,3,5`$ for $`d\le4`$ are the
partitions of $`d`$ into at most four parts (the orbits of multisets are determined by their multiplicities, sorted).

*Computation:* the Molien series up to $`t^4`$ against the orbits of monomials, for vertex and edge features — [tutorial, Section 12.4](../tutorial/12_equivariant_learning.md#s-12-4).

> **Physical meaning** (physical interpretation, not a theorem of this course). In Landau's theory of phase
> transitions the free energy is an invariant polynomial in a real order parameter $`x`$, on which $`G`$ acts by real
> orthogonal matrices. Its quadratic part is $`x^TBx`$ with $`B`$ in the commutant (as $`\pi(g)^T=\pi(g)^{-1}`$), so by
> Proposition 12.6(a) it is block diagonal over the sectors. Landau's theory then posits that the symmetric phase
> becomes unstable in the sector where an eigenvalue of $`B`$ first changes sign, and that, generically, the
> irreducible representation of that sector determines the symmetry of the ordered phase.

<a id="s-12-5"></a>

## 12.5 What to compare with

*In the library:* comparison with other software (escnn, e3nn, EMLP) and what the library adds — [tutorial, Section 12.5](../tutorial/12_equivariant_learning.md#s-12-5).

<a id="x-12"></a>

## Exercises

1. (**DeepSets with channels.**) $`S_n`$ ($`n\ge2`$) acts on $`\mathbb C^n\otimes\mathbb C^d`$ by $`P_\sigma\otimes1_d`$. Show
   that the equivariant maps from $`\mathbb C^n\otimes\mathbb C^{d_1}`$ to $`\mathbb C^n\otimes\mathbb C^{d_2}`$ are
   $`1\otimes A+J\otimes B`$ with $`A,B\in M_{d_2\times d_1}`$.
2. (**Real layers.**) Let $`\pi_1,\pi_2`$ have real matrices (e.g. permutation representations). Show that the real
   equivariant matrices form a real vector space of dimension $`\dim_{\mathbb C}\mathrm{Hom}_G(\pi_1,\pi_2)`$, and
   that $`T\mapsto\mathrm{Re}\,T`$ maps $`\mathrm{Hom}_G(\pi_1,\pi_2)`$ onto it.
3. (**Nonlinearities.**) (a) Show that applying a function $`f:\mathbb C\to\mathbb C`$ to every coordinate commutes
   with every permutation representation. (b) Show that it does not commute with the rotation by $`\pi/2`$ of
   $`\mathbb R^2`$ for $`f(t)=\max(t,0)`$. (c) Show that $`x\mapsto\varphi(\|x\|)x`$ commutes with every unitary
   representation.
4. (**How many samples determine a layer?**) Let $`T=\Phi(w)\in\pi(G)'`$ and $`x_1,\dots,x_N\in\mathbb C^n`$. Write
   $`V^*x_i`$ in sector $`c`$ as $`\sum_{a,j}(X_{i,c})_{aj}e_a\otimes e_j`$ with $`X_{i,c}\in M_{k_c\times m_c}`$, and let
   $`S_c\in M_{Nk_c\times m_c}`$ be the matrices $`X_{1,c},\dots,X_{N,c}`$ stacked. Show that the equivariant $`T`$ with
   $`Tx_1=\dots=Tx_N=0`$ form a space of dimension $`\sum_cm_c(m_c-\mathrm{rank}\,S_c)`$; for independent standard
   complex Gaussian vectors $`x_i`$ this is $`\sum_cm_c\max(0,m_c-Nk_c)`$ almost surely. Deduce that an equivariant $`T`$ is
   determined by $`Tx_1,\dots,Tx_N`$ (among equivariant maps) iff this dimension is $`0`$, and that for the edge
   representation two Gaussian samples determine an equivariant layer almost surely and one sample never does, while
   an arbitrary $`16\times16`$ matrix is never determined by fewer than sixteen samples.

*Computation:* Exercises 1, 3 and 4 numerically — [tutorial, Exercises](../tutorial/12_equivariant_learning.md#x-12).

<details>
<summary>Solutions</summary>

1. $`P_\sigma J=J=JP_\sigma`$, so these maps are equivariant. They are linearly independent: if
   $`1\otimes A+J\otimes B=0`$, the off-diagonal blocks give $`B=0`$ and then $`A=0`$. By Corollary 12.2(a) and
   Theorem 2.6, $`P_\sigma\otimes1_d\cong(\mathrm{triv}\otimes1_d)\oplus(\mathrm{std}\otimes1_d)`$, so both irreducible
   representations have multiplicity $`d`$, and Theorem 12.4 gives $`\dim\mathrm{Hom}_G=2d_1d_2`$, the dimension of
   the space of the $`1\otimes A+J\otimes B`$.

2. If $`T`$ is equivariant, so is $`\bar T`$, since $`\overline{\pi_i(g)}=\pi_i(g)`$. Hence
   $`\mathrm{Re}\,T=(T+\bar T)/2`$ and $`\mathrm{Im}\,T=(T-\bar T)/2i`$ are real and equivariant, so the map
   $`T\mapsto\mathrm{Re}\,T`$ is onto the real equivariant matrices (it is the identity on them), and
   $`\mathrm{Hom}_G(\pi_1,\pi_2)=H_{\mathbb R}\oplus iH_{\mathbb R}`$ with $`H_{\mathbb R}`$ the real equivariant
   matrices. A real basis of $`H_{\mathbb R}`$ is then a complex basis of $`\mathrm{Hom}_G`$: it spans, and a
   complex relation $`\sum_j(\alpha_j+i\beta_j)T_j=0`$ with real $`T_j`$ splits into real and imaginary parts. In
   Lemma 12.1 the orbit indicators are such a real basis. (The blocks $`y_c`$ of Theorem 12.4 are complex even for real
   data, as the $`V`$ computed by the library need not be real; parametrising $`\mathrm{Re}\,T`$ by complex $`y_c`$ uses
   twice as many real parameters as needed.)

3. (a) $`(\pi(g)f(x))_{gy}=f(x)_y=f(x_y)=f((\pi(g)x)_{gy})`$. (b) For $`R(x_1,x_2)=(-x_2,x_1)`$ and $`x=(1,-1)`$,
   $`f(Rx)=f(1,1)=(1,1)`$ but $`Rf(x)=R(1,0)=(0,1)`$. (c) $`\|\pi(g)x\|=\|x\|`$, so
   $`\varphi(\|\pi(g)x\|)\pi(g)x=\pi(g)\varphi(\|x\|)x`$.

4. $`(1\otimes w_c)(e_a\otimes e_j)=\sum_le_a\otimes(w_c)_{lj}e_l`$, so $`V^*Tx_i`$ has, in sector $`c`$, the coefficient
   matrix $`X_{i,c}w_c^T`$. Hence $`Tx_i=0`$ for all $`i`$ iff $`S_cw_c^T=0`$ for all $`c`$, i.e. iff every column of $`w_c^T`$ lies in
   $`\ker S_c`$, a space of dimension $`m_c-\mathrm{rank}\,S_c`$: the $`w_c`$ with this property form a space of
   dimension $`m_c(m_c-\mathrm{rank}\,S_c)`$, and $`\Phi`$ is injective. For Gaussian $`x_i`$, the coordinates of
   $`V^*x_i`$ are those of $`x_i`$ in the orthonormal basis of the columns of $`V`$, again independent standard complex
   Gaussians (Proposition A.16(a)), and independent over $`i`$; so $`S_c`$ is an $`Nk_c\times m_c`$ Ginibre matrix, which
   has a density (Proposition A.16(b)). Its rank is at most $`\min(Nk_c,m_c)`$, with equality unless the leading
   square submatrix of that size is singular; $`|\det|^2`$ of that submatrix is a real polynomial in the real and
   imaginary parts of the entries, not identically zero (it is $`1`$ when the submatrix is the identity), so by
   Corollary A.7(b) the rank is $`\min(Nk_c,m_c)`$ almost surely. For the edge representation,
   $`(k,m)=(1,2),(2,1),(3,1),(3,3)`$ (the computed decomposition), only the sector $`(1,2)`$ contributes:
   $`2\cdot(2-1)=2`$ for $`N=1`$, and $`0`$ for $`N=2`$. An arbitrary $`16\times16`$ matrix is the case of the trivial group
   (one sector, $`k=1`$, $`m=16`$): the dimension is $`16(16-\mathrm{rank}\,S)\ge16(16-N)>0`$ for every choice of
   $`N<16`$ samples, as $`\mathrm{rank}\,S\le N`$, and $`0`$ almost surely for $`N=16`$ Gaussian samples. Since
   $`T\mapsto(Tx_1,\dots,Tx_N)`$ is linear, two equivariant maps agree on the samples iff their difference lies in the
   space just computed, which gives the criterion. For $`N=1`$ the space has dimension at least $`2`$ for *every* $`x_1`$,
   because $`\mathrm{rank}\,S_c\le Nk_c=1<m_c=2`$ in the sector $`(1,2)`$.

</details>

**Where to go from here.** The last chapter, [Chapter 13](13_jordan_algebras.md), asks why quantum mechanics uses *-algebras rather than Jordan algebras. [Chapter 4](04_symmetries.md) develops the
representation theory used here: Schur's lemma, characters, and the randomised algorithm behind `dc.from_group`
(Sections 4.1–4.3). [Chapter 9](09_spin_chains.md) applies the same sectors to spin chains, and the Haar measure on
$`\prod_cU(m_c)`$ of [Chapter 5](05_random_matrices.md) gives the natural random initialisation of an equivariant
unitary layer. For the mathematics see J.-P. Serre, *Linear Representations of Finite Groups* (Springer 1977) and
W. Fulton and J. Harris, *Representation Theory* (Springer 1991); for equivariant learning, M. Bronstein, J. Bruna,
T. Cohen and P. Veličković, *Geometric Deep Learning* (2021), and M. Weiler et al., *Equivariant and Coordinate
Independent Convolutional Networks* (2023); for the decomposition algorithm, the references in
`torch_vn_algebra/decompose.py`.
