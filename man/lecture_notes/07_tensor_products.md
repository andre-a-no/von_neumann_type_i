**English** · [Español](es/07_tensor_products.md) · [Français](fr/07_tensor_products.md) · [Deutsch](de/07_tensor_products.md) · [中文](zh/07_tensor_products.md) · [日本語](ja/07_tensor_products.md) · [Русский](ru/07_tensor_products.md)

# 7. Tensor products and entanglement

This chapter studies the tensor product of two algebras with sectors (physically, the algebra of
two systems together; see the *Physical meaning* boxes). With sectors there are two products: the
plain one, in which each part keeps its own charge, and the *fused* one, in which only the total
charge is conserved. We prove the structure of
both, define the partial trace as the dual of an embedding, and study entanglement.

**You will learn** how sectors multiply, why the partial trace is
a channel and, normalised, a conditional expectation of finite index, what the reduced state of a charge-conserving state is, how entanglement splits into a
number part and a configurational part, and the average purity and entropy of random states. The numerical checks of this chapter are in
[Tutorial 7](../tutorial/07_tensor_products.md).

Throughout, $`\mathcal M_1=\bigoplus_{c=1}^{C_1}M_{k_c}(\mathbb C)`$ acts on
$`H_1=\bigoplus_c\mathbb C^{k_c}`$ with minimal central projections $`P_c`$, and
$`\mathcal M_2=\bigoplus_{d=1}^{C_2}M_{l_d}(\mathbb C)`$ acts on $`H_2=\bigoplus_d\mathbb C^{l_d}`$
with minimal central projections $`Q_d`$. $`\mathrm{Tr}`$ is the trace of the space on which an
operator acts; on $`\mathcal M_i`$ it is the blunt trace of Chapter 3. For $`X`$ on
$`\mathbb C^k\otimes\mathbb C^l`$ we write $`X_{(a,b),(a',b')}=\langle e_a\otimes e_b,X\,e_{a'}\otimes e_{b'}\rangle`$,
and we identify $`\mathbb C^k\otimes\mathbb C^l`$ with $`\mathbb C^{kl}`$ by $`e_a\otimes e_b\mapsto e_{(a-1)l+b}`$
(the Kronecker convention), so that $`(A\otimes B)_{(a,b),(a',b')}=A_{aa'}B_{bb'}`$.

<a id="s-7-1"></a>

## 7.1 The tensor product of two algebras

For $`*`$-algebras $`\mathcal M_1`$ on $`H_1`$ and $`\mathcal M_2`$ on $`H_2`$, the **tensor product**
$`\mathcal M_1\otimes\mathcal M_2`$ is the linear span of the operators $`A\otimes B`$ on
$`H_1\otimes H_2`$, $`A\in\mathcal M_1`$, $`B\in\mathcal M_2`$. Since
$`(A\otimes B)(A'\otimes B')=AA'\otimes BB'`$, $`(A\otimes B)^*=A^*\otimes B^*`$ and $`1\otimes1=1`$, it
is a $`*`$-algebra, hence a von Neumann algebra (Theorem 2.2).

<a id="r-7-1"></a>

**Theorem 7.1.** (a) $`H_1\otimes H_2=\bigoplus_{c,d}\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$, an
orthogonal sum of the ranges of the projections $`P_c\otimes Q_d`$, and
```math
\mathcal M_1\otimes\mathcal M_2=\bigoplus_{c,d}M_{k_c}\otimes M_{l_d}\cong\bigoplus_{c,d}M_{k_cl_d}(\mathbb C):
```
it consists of all operators $`X=\sum_{c,d}(P_c\otimes Q_d)X(P_c\otimes Q_d)`$, with arbitrary blocks
on the $`\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$. (b) Its minimal central projections are the $`C_1C_2`$
projections $`P_c\otimes Q_d`$, and $`\dim(\mathcal M_1\otimes\mathcal M_2)=\dim\mathcal M_1\dim\mathcal M_2`$.
(c) For traces $`\mathrm{tr}_w`$ on $`\mathcal M_1`$ and $`\mathrm{tr}_v`$ on $`\mathcal M_2`$ (Theorem 3.1),
$`\mathrm{tr}_w(A)\mathrm{tr}_v(B)=\mathrm{tr}_{w\otimes v}(A\otimes B)`$ with
$`(w\otimes v)_{cd}=w_cv_d`$. In particular `Tr_blunt`, `Tr_norm` and `tau_vN` of the product are the
products of the same traces of the factors.

*Proof.* (a) The $`P_c\otimes Q_d`$ are mutually orthogonal projections with sum $`1\otimes1`$, and the
range of $`P_c\otimes Q_d`$ is $`\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$. As $`P_c`$, $`Q_d`$ are central,
$`(P_c\otimes Q_d)(A\otimes B)=AP_c\otimes BQ_d=(A\otimes B)(P_c\otimes Q_d)`$, so every element of
$`\mathcal M_1\otimes\mathcal M_2`$ is block diagonal. Conversely, with the matrix units
$`E^c_{aa'}\in\mathcal M_1`$ and $`F^d_{bb'}\in\mathcal M_2`$ of the blocks $`c`$ and $`d`$,
$`E^c_{aa'}\otimes F^d_{bb'}=|e_a\otimes e_b\rangle\langle e_{a'}\otimes e_{b'}|`$ inside
$`\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$; these are the matrix units of $`M_{k_cl_d}(\mathbb C)`$, so every
operator supported in the block $`(c,d)`$ lies in $`\mathcal M_1\otimes\mathcal M_2`$. (b) The centre of a
direct sum of full matrix algebras is spanned by the projections onto the summands (Chapter 2,
Exercise 1), and $`\sum_{c,d}k_c^2l_d^2=\sum_ck_c^2\sum_dl_d^2`$. (c)
$`\mathrm{Tr}(X\otimes Y)=\sum_{a,b}X_{aa}Y_{bb}=\mathrm{Tr}\,X\mathrm{Tr}\,Y`$, so
$`\mathrm{tr}_{w\otimes v}(A\otimes B)=\sum_{c,d}w_cv_d\mathrm{Tr}\,A_c\mathrm{Tr}\,B_d=\mathrm{tr}_w(A)\mathrm{tr}_v(B)`$.
The weights $`1`$, $`1/k_c`$, $`1/(C_1k_c)`$ of the three traces (Section 3.1) multiply to $`1`$,
$`1/(k_cl_d)`$, $`1/(C_1C_2k_cl_d)`$, the weights of the same traces on the product. $`\square`$

*In the library:* the product is `tensor_product` and $`A\otimes B`$ is `kron`; a check of Theorem 7.1 on $`(\mathbb C\oplus M_2)\otimes(M_2\oplus M_3)`$ — [tutorial, Section 7.1](../tutorial/07_tensor_products.md#s-7-1).

> **Physical meaning** (physical interpretation, not a theorem of this course). If each of two systems obeys its own superselection rule — two boxes
> that cannot exchange particles — the joint system has a sector for every *pair* of charges, and
> no observable connects different pairs. A joint sector counts with the product of the weights
> (Theorem 7.1(c)).

<a id="s-7-2"></a>

## 7.2 The partial trace

**Definition.** Let $`\rho\in\mathcal M_1\otimes\mathcal M_2`$. The **partial trace**
$`\mathrm{Tr}_2\,\rho\in\mathcal M_1`$ is the element with
```math
\mathrm{Tr}\big((\mathrm{Tr}_2\,\rho)\,A\big)=\mathrm{Tr}\big(\rho\,(A\otimes1)\big)\qquad\text{for all }A\in\mathcal M_1 .
```
It exists and is unique by Theorem 3.2, applied to the linear functional
$`A\mapsto\mathrm{Tr}(\rho(A\otimes1))`$ and the faithful trace $`\mathrm{Tr}`$ of $`\mathcal M_1`$: the
partial trace is the dual, with respect to the blunt traces, of the embedding $`A\mapsto A\otimes1`$.
$`\mathrm{Tr}_1`$ is defined in the same way; for a state $`\rho`$, $`\mathrm{Tr}_2\,\rho`$ and
$`\mathrm{Tr}_1\,\rho`$ are the **reduced states**. For the full algebras $`\mathcal M_i=B(H_i)`$ this is the
**Hilbert-space partial trace** $`\mathrm{Tr}_{H_2}:B(H_1\otimes H_2)\to B(H_1)`$.

<a id="r-7-2"></a>

**Proposition 7.2.** (a) For $`\rho=\bigoplus_{c,d}\rho_{cd}`$,
```math
(\mathrm{Tr}_2\,\rho)_c=\sum_d\mathrm{tr}_{l_d}\,\rho_{cd},\qquad
(\mathrm{tr}_l\,X)_{aa'}=\sum_{b=1}^lX_{(a,b),(a',b)}\quad(X\text{ on }\mathbb C^k\otimes\mathbb C^l),
```
and $`\mathrm{Tr}_2\,\rho=\mathrm{Tr}_{H_2}\,\rho`$, where $`(\mathrm{Tr}_{H_2}\,X)_{ii'}=\sum_jX_{(i,j),(i',j)}`$.
(b) $`\mathrm{Tr}_2(A\otimes B)=(\mathrm{Tr}\,B)\,A`$.
(c) $`\mathrm{Tr}_2\big((A\otimes1)X(A'\otimes1)\big)=A\,(\mathrm{Tr}_2\,X)\,A'`$ for $`A,A'\in\mathcal M_1`$.
(d) $`\mathrm{Tr}(\mathrm{Tr}_2\,X)=\mathrm{Tr}\,X`$, and $`\mathrm{Tr}_2`$ maps states to states.

*Proof.* (a) With $`A=\bigoplus_cA_c`$, $`A\otimes1=\bigoplus_{c,d}A_c\otimes1_{l_d}`$ and
```math
\mathrm{Tr}\big(\rho(A\otimes1)\big)=\sum_{c,d}\sum_{a,a',b}(\rho_{cd})_{(a,b),(a',b)}(A_c)_{a'a}
=\sum_c\mathrm{Tr}\Big(\Big(\sum_d\mathrm{tr}_{l_d}\,\rho_{cd}\Big)A_c\Big),
```
so the element on the right satisfies the definition and is $`\mathrm{Tr}_2\,\rho`$. The same
computation for $`B(H_i)`$ gives the formula for $`\mathrm{Tr}_{H_2}`$; for block-diagonal $`\rho`$ the
entries $`\rho_{(i,j),(i',j)}`$ with $`i,i'`$ in different blocks vanish, so $`\mathrm{Tr}_{H_2}\,\rho`$ has
the blocks of (a). (b) $`\mathrm{tr}_l(X\otimes Y)=(\mathrm{Tr}\,Y)X`$. (c) For $`B\in\mathcal M_1`$,
$`\mathrm{Tr}\big(\mathrm{Tr}_2((A\otimes1)X(A'\otimes1))B\big)=\mathrm{Tr}\big(X(A'BA\otimes1)\big) =\mathrm{Tr}\big(A(\mathrm{Tr}_2\,X)A'B\big)`$; use uniqueness. (d) Take $`A=1`$. If $`\rho`$ is a state,
$`A\mapsto\mathrm{Tr}(\rho(A\otimes1))`$ is positive, as $`A^*A\otimes1=(A\otimes1)^*(A\otimes1)`$, and
is $`1`$ at $`A=1`$, so its density $`\mathrm{Tr}_2\,\rho`$ is a state (Theorem 3.2). $`\square`$

*In the library:* $`\mathrm{Tr}_2`$ is `partial_trace`; a check of Proposition 7.2(a), (b) against `numpy.einsum` — [tutorial, Section 7.2](../tutorial/07_tensor_products.md#s-7-2).

Channels are the subject of Chapter 6; here we need only this definition: a linear map
$`\Phi`$ is **completely positive** if $`\mathrm{id}_n\otimes\Phi`$ is positive on
$`M_n(\mathbb C)\otimes(\cdot)`$ for every $`n`$, and a **channel** if it also preserves the trace.

<a id="r-7-3"></a>

**Theorem 7.3.** Let $`j`$ run over the standard basis of $`H_2`$ and $`K_j=1\otimes\langle j|:H_1\otimes H_2\to H_1`$,
$`K_j(x\otimes y)=\langle j,y\rangle x`$. Then
```math
\mathrm{Tr}_{H_2}\,X=\sum_jK_jXK_j^*,\qquad\sum_jK_j^*K_j=1,\qquad\sum_jK_j^*AK_j=A\otimes1 .
```
Hence $`\mathrm{Tr}_2`$ is a channel from $`\mathcal M_1\otimes\mathcal M_2`$ to $`\mathcal M_1`$, and its
dual with respect to the blunt traces is $`A\mapsto A\otimes1`$.

*Proof.* $`K_j^*x=x\otimes e_j`$, so $`\langle e_a,K_jXK_j^*e_{a'}\rangle=X_{(a,j),(a',j)}`$; summing over $`j`$
gives Proposition 7.2(a). Further $`\sum_jK_j^*K_j(x\otimes y)=\sum_j\langle j,y\rangle x\otimes e_j=x\otimes y`$
and $`\sum_jK_j^*AK_j=\sum_jA\otimes|j\rangle\langle j|=A\otimes1`$. For $`Y\ge0`$,
$`(\mathrm{id}_n\otimes\mathrm{Tr}_{H_2})(Y)=\sum_j(1_n\otimes K_j)Y(1_n\otimes K_j)^*\ge0`$, because
$`\langle x,TYT^*x\rangle=\langle T^*x,YT^*x\rangle\ge0`$; this persists on the subalgebra
$`\mathcal M_1\otimes\mathcal M_2`$. Finally $`\mathrm{Tr}\,\sum_jK_jXK_j^*A=\mathrm{Tr}\,X\sum_jK_j^*AK_j`$,
which for $`A=1`$ is trace preservation and in general the duality. $`\square`$

*In the library:* the channel of Theorem 7.3 is `composite.partial_trace_channel`, with its dual — [tutorial, Section 7.2](../tutorial/07_tensor_products.md#s-7-2).

> **Physical meaning** (physical interpretation, not a theorem of this course). The reduced state $`\mathrm{Tr}_2\,\rho`$ is all that an observer who can
> only measure the first system learns — the definition says it reproduces every expectation of an
> observable $`A\otimes1`$. Joint states with equal reduced states differ in their correlations.

**The partial trace as a conditional expectation.** Let $`\mathcal M_2=M_m(\mathbb C)`$, so that the
sectors of $`\mathcal M_1\otimes M_m`$ are those of $`\mathcal M_1`$ (Theorem 7.1), let
$`\iota(A)=A\otimes1_m`$ and $`\mathcal E=\frac1m\mathrm{Tr}_2:\mathcal M_1\otimes M_m\to\mathcal M_1`$.
Let $`e_{ab}`$ be the matrix units of $`M_m`$.

<a id="r-7-4"></a>

**Proposition 7.4.** (a) $`\iota\circ\mathcal E`$ is a conditional expectation of $`\mathcal M_1\otimes M_m`$
onto $`\mathcal M_1\otimes1`$ (Section 2.6), and for every faithful trace $`\mathrm{tr}_w`$ of
$`\mathcal M_1\otimes M_m`$ it is the only one that preserves $`\mathrm{tr}_w`$. (b) The $`m^2`$ elements
$`u_{ab}=\sqrt m\,(1\otimes e_{ab})`$ satisfy $`X=\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*X))`$ for all $`X`$,
and $`\sum_{a,b}u_{ab}u_{ab}^*=m^2\,1`$. (c) Every family $`(v_k)`$ with
$`X=\sum_kv_k\,\iota(\mathcal E(v_k^*X))`$ for all $`X`$ has $`\sum_kv_kv_k^*=m^2\,1`$.

*Proof.* (a) $`\mathcal E`$ is positive (Theorem 7.3) and so is $`\iota`$
($`T^*T\otimes1=(T\otimes1)^*(T\otimes1)`$); $`\mathcal E(A\otimes1)=A`$ by Proposition 7.2(b), and the
bimodule property is Proposition 7.2(c). Write $`\mathrm{tr}_w`$ also for the trace with the same
weights $`w_c`$ on $`\mathcal M_1`$; then $`\mathrm{tr}_w(Y\otimes1)=m\mathrm{tr}_w(Y)`$. The definition of
$`\mathrm{Tr}_2`$ applied to $`AP_c`$, multiplied by $`w_c`$ and summed over $`c`$, gives
$`\mathrm{tr}_w(X(A\otimes1))=\mathrm{tr}_w((\mathrm{Tr}_2\,X)A)`$, and $`A=1`$ gives
$`\mathrm{tr}_w(\iota\mathcal E(X))=\mathrm{tr}_w(X)`$. If $`F=\iota\circ F_0`$ is another conditional
expectation onto $`\mathcal M_1\otimes1`$ preserving $`\mathrm{tr}_w`$, then, as $`F(X(A\otimes1))=F(X)(A\otimes1)`$,
$`\mathrm{tr}_w(X(A\otimes1))=\mathrm{tr}_w(F(X)(A\otimes1))=m\mathrm{tr}_w(F_0(X)A)`$, so
$`\mathrm{tr}_w((\mathrm{Tr}_2\,X-mF_0(X))A)=0`$ for all $`A\in\mathcal M_1`$, and $`F_0=\mathcal E`$ by faithfulness
(take $`A=(\mathrm{Tr}_2\,X-mF_0(X))^*`$). (b) Write $`X=\sum_{a,b}X_{ab}\otimes e_{ab}`$ with
$`X_{ab}\in\mathcal M_1`$ (Theorem 7.1(a)). By Proposition 7.2(b),
$`\mathcal E(u_{ab}^*X)=\frac{\sqrt m}m\sum_{b'}\mathrm{Tr}_2(X_{ab'}\otimes e_{bb'})=X_{ab}/\sqrt m`$, so
$`u_{ab}\,\iota(\mathcal E(u_{ab}^*X))=X_{ab}\otimes e_{ab}`$; and
$`\sum_{a,b}u_{ab}u_{ab}^*=m\sum_{a,b}1\otimes e_{aa}=m^21`$. (c) $`\mathrm{Tr}_2(Y^*)=(\mathrm{Tr}_2\,Y)^*`$ by
Proposition 7.2(a), so the adjoint of the expansion of $`Y=u_{ab}`$ is
$`u_{ab}^*=\sum_k\iota(\mathcal E(u_{ab}^*v_k))\,v_k^*`$, and by (b)
$`\sum_{a,b}u_{ab}u_{ab}^*=\sum_k\big(\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*v_k))\big)v_k^*=\sum_kv_kv_k^*`$. $`\square`$

A family as in (c) is a *quasi-basis* of $`\mathcal E`$, and $`\sum_kv_kv_k^*`$, which by (c) does not
depend on it and by (b) exists, is the **index** $`\mathrm{Ind}\,\mathcal E=m^2`$ of $`\mathcal E`$ (a
definition; the name follows Y. Watatani, Mem. Amer. Math. Soc. 83, no. 424, 1990). The
**inclusion matrix** of an inclusion of $`\bigoplus_cM_{k_c}`$ into $`\bigoplus_sM_{n_s}`$ has the entries
$`\Lambda_{sc}`$, the numbers of copies of the block $`c`$ of the subalgebra in the block $`s`$ of the larger one
(Chapter 9 uses it for towers of algebras). For a factor, $`M_k\otimes1\subset M_k\otimes M_m`$: ordering
the basis $`e_a\otimes e_b`$ of $`\mathbb C^{km}`$ by $`b`$ first, $`A\otimes1_m=\mathrm{diag}(A,\dots,A)`$ with $`m`$
copies, so $`\Lambda=(m)`$, and $`m^2=\|\Lambda\|^2=\dim(M_k\otimes M_m)/\dim M_k`$. (Further
reading, not used here: V. F. R. Jones, Invent. Math. 72, 1, 1983, on the index of subfactors, and
F. M. Goodman, P. de la Harpe and V. F. R. Jones, *Coxeter graphs and towers of algebras*, Springer,
1989, on inclusion matrices.) For general $`\mathcal M_1`$ the inclusion is the direct sum of the
$`C_1`$ inclusions $`M_{k_c}\otimes1\subset M_{k_c}\otimes M_m`$, $`\Lambda=m1_{C_1}`$.

*In the library:* `inclusions.from_tensor_product`, with the conditional expectation, the quasi-basis and the index of Proposition 7.4 — [tutorial, Section 7.2](../tutorial/07_tensor_products.md#s-7-2).

<a id="s-7-3"></a>

## 7.3 Fusion of charges

Suppose now that sector $`c`$ of $`\mathcal M_1`$ carries a charge $`q_c`$, sector $`d`$ of $`\mathcal M_2`$ a
charge $`r_d`$, and that a **fusion rule** $`f`$ gives the total charge $`f(q_c,r_d)`$ — addition for
particle numbers, addition modulo $`2`$ for fermion parity, or any other map. Put
```math
H_q=\bigoplus_{f(q_c,r_d)=q}\mathbb C^{k_c}\otimes\mathbb C^{l_d},\qquad
\Pi_q=\sum_{f(q_c,r_d)=q}P_c\otimes Q_d,\qquad n_q=\sum_{f(q_c,r_d)=q}k_cl_d .
```
The **fused product** $`\mathcal M_1\otimes_f\mathcal M_2=\{\Pi_q\}'`$ consists of the operators on
$`H_1\otimes H_2`$ that conserve the total charge.

<a id="r-7-5"></a>

**Proposition 7.5.** (a) $`\mathcal M_1\otimes_f\mathcal M_2=\bigoplus_qB(H_q)\cong\bigoplus_qM_{n_q}(\mathbb C)`$,
with minimal central projections $`\Pi_q`$. (b) $`\mathcal M_1\otimes\mathcal M_2\subseteq\mathcal M_1\otimes_f\mathcal M_2`$,
with equality iff every $`H_q`$ contains only one pair $`(c,d)`$. (c) If the charges are real numbers
and $`f`$ is addition, then $`\mathcal M_1\otimes_f\mathcal M_2=\{N\}'`$ with
$`N=N_1\otimes1+1\otimes N_2`$, $`N_1=\sum_cq_cP_c`$, $`N_2=\sum_dr_dQ_d`$.

*Proof.* (a) Every pair lies in exactly one $`H_q`$, so the $`\Pi_q`$ are mutually orthogonal
projections with sum $`1`$, and $`X`$ commutes with all of them iff $`X=\sum_q\Pi_qX\Pi_q`$, with
arbitrary blocks in $`B(H_q)`$; the centre is as in Theorem 7.1(b). (b) Every element of
$`\mathcal M_1\otimes\mathcal M_2`$ commutes with each $`P_c\otimes Q_d`$, hence with their sums $`\Pi_q`$.
If each $`H_q`$ is one $`\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$, both algebras are
$`\bigoplus_{c,d}B(\mathbb C^{k_c}\otimes\mathbb C^{l_d})`$; if $`H_q`$ contains two pairs, $`B(H_q)`$ contains
operators mapping one pair block into the other, which are not block diagonal. (c)
$`N=\sum_{c,d}(q_c+r_d)P_c\otimes Q_d=\sum_qq\,\Pi_q`$, so the nonzero $`\Pi_q`$ are the spectral
projections of $`N`$. They are polynomials in $`N`$ (proof of Theorem 1.5(a)) and $`N`$ is a combination
of them, so $`X`$ commutes with $`N`$ iff it commutes with every $`\Pi_q`$. $`\square`$

**Example: hard-core bosons.** On $`L`$ sites, $`(\mathbb C^2)^{\otimes L}`$ with the occupation number
$`N=\sum_in^{(i)}`$, $`n=\mathrm{diag}(0,1)`$, the eigenspace $`N=q`$ is spanned by the $`\binom Lq`$ product
basis vectors with $`q`$ occupied sites, so the number-conserving operators $`\{N\}'`$ form
$`\bigoplus_qM_{\binom Lq}`$ (Proposition 4.10(2)): for two sites $`M_1\oplus M_2\oplus M_1`$ (charges $`0,1,2`$), for
three sites $`M_1\oplus M_3\oplus M_3\oplus M_1`$. By Proposition 7.5(c) the fused product is $`\{N\}'`$ for
the total occupation number of five sites, up to the order of the basis, and indeed
$`n_q=\sum_c\binom2c\binom3{q-c}=\binom5q`$ (Vandermonde's identity: a $`q`$-element subset of five sites
consists of $`c`$ of the first two and $`q-c`$ of the last three).

*In the library:* the fused product is `fused_tensor_product`; the hard-core bosons and Proposition 7.5 — [tutorial, Section 7.3](../tutorial/07_tensor_products.md#s-7-3).

The fused algebra has dimension $`\sum_q\binom5q^2=\binom{10}5=252`$ (Vandermonde again), the plain
product $`(1+4+1)(1+9+9+1)=120`$ (Theorem 7.1(b)). By Theorem 7.1(a) and Proposition 7.5(b) the plain
product consists of the elements of the fused one that are block diagonal for the pairs $`(c,d)`$; the
remaining $`132`$ dimensions are spanned by the operators that map a pair block $`(c,d)`$ into a
different pair block $`(c',d')`$ of the same $`H_q`$, $`c'\ne c`$ (so $`d'\ne d`$): they change the charge
of each part while conserving the total, i.e. here they move particles between the two parts.

A state of $`\mathcal M_1\otimes_f\mathcal M_2`$ can have coherences between pair blocks, so Proposition
7.2 does not apply. Let $`E_1(X)=\sum_cP_cXP_c`$ on $`B(H_1)`$ (the conditional expectation of Section
2.6); it maps onto $`\mathcal M_1`$, and $`\mathrm{Tr}(E_1(X)A)=\mathrm{Tr}(XA)`$ for
$`A\in\mathcal M_1`$, because $`A=\sum_cP_cAP_c`$.

<a id="r-7-6"></a>

**Theorem 7.6.** Let $`\rho`$ be the density of a state $`\omega`$ on $`\mathcal M_1\otimes_f\mathcal M_2`$. The
density of its restriction $`A\mapsto\omega(A\otimes1)`$ to $`\mathcal M_1`$ is
```math
\rho_1=E_1\big(\mathrm{Tr}_{H_2}\,\rho\big).
```
The block $`(c,c')`$ of $`\mathrm{Tr}_{H_2}\,\rho`$ is $`\sum_d\mathrm{tr}_{l_d}\,\rho_{(c,d),(c',d)}`$, summed over
the $`d`$ with $`f(q_c,r_d)=f(q_{c'},r_d)`$, where $`\rho_{(c,d),(c',d)}`$ is the block of $`\rho`$ from the
pair $`(c',d)`$ to the pair $`(c,d)`$. Hence $`E_1`$ acts trivially on all $`\mathrm{Tr}_{H_2}\,\rho`$ iff there
are no $`c\neq c'`$ and $`d`$ with $`f(q_c,r_d)=f(q_{c'},r_d)`$; for instance if $`f`$ is addition and the
$`q_c`$ are pairwise distinct.

*Proof.* $`A\otimes1\in\mathcal M_1\otimes_f\mathcal M_2`$ (Proposition 7.5(b)) and
$`\omega(A\otimes1)=\mathrm{Tr}(\rho(A\otimes1))=\mathrm{Tr}\big((\mathrm{Tr}_{H_2}\,\rho)A\big) =\mathrm{Tr}\big(E_1(\mathrm{Tr}_{H_2}\,\rho)A\big)`$. As $`E_1(\mathrm{Tr}_{H_2}\,\rho)\in\mathcal M_1`$, it is
the density by the uniqueness in Theorem 3.2. The block formula is Proposition 7.2(a) for
$`\mathrm{Tr}_{H_2}`$, with $`i`$ in block $`c`$, $`i'`$ in block $`c'`$ and $`j`$ in block $`d`$; the block of
$`\rho`$ between $`(c',d)`$ and $`(c,d)`$ vanishes unless both pairs lie in the same $`H_q`$. If $`c\neq c'`$
and $`d`$ with $`f(q_c,r_d)=f(q_{c'},r_d)=q`$ exist, take unit vectors $`x\in\mathbb C^{k_c}`$,
$`x'\in\mathbb C^{k_{c'}}`$, $`y\in\mathbb C^{l_d}`$ and $`\psi=(x+x')\otimes y/\sqrt2\in H_q`$. Then
$`|\psi\rangle\langle\psi|\in B(H_q)`$, and $`\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|=\frac12|x+x'\rangle\langle x+x'|`$
has the nonzero block $`\frac12|x\rangle\langle x'|`$, which $`E_1`$ removes. $`\square`$

By the block formula of Theorem 7.6, $`(\rho_1)_c=\sum_d\mathrm{tr}_{l_d}\,\rho_{(c,d),(c,d)}`$.

*In the library:* `partial_trace` and `composite.partial_trace_channel` compute $`\rho_1`$ for any fusion rule; a check of Theorem 7.6 for addition and for the parity rule — [tutorial, Section 7.3](../tutorial/07_tensor_products.md#s-7-3).

> **Physical meaning** (physical interpretation, not a theorem of this course). In the BCS
> description of a superconductor only the parity of the electron
> number is conserved, so the global state may superpose the particle numbers $`0`$ and $`2`$ of a
> region. An observer whose observables conserve the particle number of the region cannot see
> this coherence: for them the reduced state is the pinched one of Theorem 7.6 (the mathematical
> content of this box), a mixture of particle numbers.

<a id="s-7-4"></a>

## 7.4 Entanglement of pure states

<a id="r-7-7"></a>

**Lemma 7.7 (Schmidt decomposition).** Every unit vector $`\psi\in\mathbb C^a\otimes\mathbb C^b`$ can be
written $`\psi=\sum_{i=1}^rs_i\,u_i\otimes v_i`$ with $`s_i>0`$, $`\sum_is_i^2=1`$ and orthonormal
$`u_i\in\mathbb C^a`$, $`v_i\in\mathbb C^b`$. Then $`\mathrm{Tr}_2|\psi\rangle\langle\psi|=\sum_is_i^2|u_i\rangle\langle u_i|`$
and $`\mathrm{Tr}_1|\psi\rangle\langle\psi|=\sum_is_i^2|v_i\rangle\langle v_i|`$: the two reduced states
have the same nonzero eigenvalues $`s_i^2`$.

*Proof.* Let $`\Psi`$ be the $`a\times b`$ matrix $`\Psi_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,\psi\rangle`$. By the
spectral theorem $`\Psi\Psi^*=\sum_{i=1}^a\lambda_i|u_i\rangle\langle u_i|`$ with an orthonormal basis
$`u_i`$ and $`\lambda_i=\|\Psi^*u_i\|^2\ge0`$, ordered so that $`\lambda_i>0`$ exactly for $`i\le r`$. Put
$`s_i=\sqrt{\lambda_i}`$ and $`w_i=\Psi^*u_i/s_i`$ for $`i\le r`$; then
$`\langle w_i,w_j\rangle=\langle u_i,\Psi\Psi^*u_j\rangle/(s_is_j)=\delta_{ij}`$, and $`\Psi^*u_i=0`$ for $`i>r`$.
Hence $`\Psi=\sum_i|u_i\rangle\langle u_i|\Psi=\sum_{i\le r}|u_i\rangle\langle\Psi^*u_i|=\sum_{i\le r}s_i|u_i\rangle\langle w_i|`$,
i.e. $`\Psi_{\alpha\beta}=\sum_is_i(u_i)_\alpha\overline{(w_i)_\beta}`$: the claim with $`v_i=\bar w_i`$, and
$`\sum_is_i^2=\mathrm{Tr}\,\Psi\Psi^*=\|\psi\|^2=1`$. Finally
$`|\psi\rangle\langle\psi|=\sum_{i,j}s_is_j|u_i\rangle\langle u_j|\otimes|v_i\rangle\langle v_j|`$, and
Proposition 7.2(b) gives the reduced states. $`\square`$

The **entanglement entropy** of $`\psi`$ is
$`E(\psi)=S(\mathrm{Tr}_2|\psi\rangle\langle\psi|)=S(\mathrm{Tr}_1|\psi\rangle\langle\psi|)=-\sum_is_i^2\log s_i^2`$.
It vanishes iff $`r=1`$ (each term $`-s_i^2\log s_i^2`$ is $`\ge0`$, and $`0`$ only for $`s_i^2=1`$), i.e. iff
$`\psi=u\otimes v`$ (then $`\Psi=uv^T`$ has rank $`1`$, and $`r=\mathrm{rank}\,\Psi\Psi^*=\mathrm{rank}\,\Psi`$), and
$`E(\psi)\le\log\min(a,b)`$ with equality iff the smaller reduced state is maximally mixed
(Corollary 3.6, applied to both reduced states).

With charges we take a unit vector $`\psi`$ with a definite total charge, i.e. in one sector $`H_q`$ of the
fused product (a superposition of different $`H_q`$ gives on $`\mathcal M_1\otimes_f\mathcal M_2`$ the same
state as a mixture, Chapter 3, Exercise 3); write $`\psi=\sum_{(c,d)}\psi_{cd}`$ with
$`\psi_{cd}\in\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$.

<a id="r-7-8"></a>

**Proposition 7.8 (symmetry-resolved entanglement).** Let $`f`$ be addition and let the $`q_c`$ be
pairwise distinct and the $`r_d`$ pairwise distinct, so that for each $`c`$ at most one $`d=d(c)`$ has
$`q_c+r_d=q`$. Put $`\psi_c=\psi_{c\,d(c)}`$ (or $`0`$), $`p_c=\|\psi_c\|^2`$ and, for $`p_c>0`$,
$`\hat\psi_c=\psi_c/\sqrt{p_c}`$. Then $`\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|`$ lies in $`\mathcal M_1`$ (it
commutes with $`N_1`$) and equals $`\bigoplus_cp_c\mathrm{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|`$; the
reduced state on $`\mathcal M_2`$ has the sector probability $`p_c`$ at $`d(c)`$; and
```math
E(\psi)=H(p)+\sum_cp_c\,E(\hat\psi_c),\qquad H(p)=-\sum_cp_c\log p_c .
```

*Proof.* By the block formula of Theorem 7.6, the block $`(c,c')`$ of $`\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|`$
is the sum of $`\mathrm{tr}_{l_d}|\psi_{cd}\rangle\langle\psi_{c'd}|`$ over the $`d`$ with
$`q_c+r_d=q=q_{c'}+r_d`$, which forces $`c=c'`$ and $`d=d(c)`$. So the reduced state is
$`\bigoplus_c\mathrm{Tr}_2|\psi_c\rangle\langle\psi_c|`$, with block traces $`p_c`$; exchanging the roles of
the factors gives the statement for $`\mathcal M_2`$. By Lemma 7.7, $`E(\psi)`$ is the entropy of this
state, and Theorem 3.3 with $`\sigma_c=\mathrm{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|`$, $`S(\sigma_c)=E(\hat\psi_c)`$,
gives the formula. $`\square`$

$`H(p)`$ is the **number entanglement** and $`\sum_cp_cE(\hat\psi_c)`$ the **configurational
entanglement**.

*Computation:* Lemma 7.7 and Proposition 7.8 for a random state of two particles on five sites — [tutorial, Section 7.4](../tutorial/07_tensor_products.md#s-7-4).

> **Physical meaning** (physical interpretation, not a theorem of this course). Parties bound by
> the particle-number superselection rule can use only the
> configurational part as a resource: the *entanglement of particles* of H. M. Wiseman and J. A.
> Vaccaro (Phys. Rev. Lett. 91, 097902, 2003). $`H(p)`$ is measured in cold atoms as the *number
> entropy* (A. Lukin et al., Science 364, 256, 2019).

<a id="s-7-5"></a>

## 7.5 Random states and Page's formula

Let $`g\in\mathbb C^N`$ have independent coordinates with independent $`N(0,\tfrac12)`$ real and imaginary
parts (a standard complex Gaussian vector, Proposition A.16(a)), so $`\mathbb E|g_i|^2=1`$; $`g\neq0`$
almost surely (Theorem A.20). A **Haar-random pure state** is $`\psi=g/\|g\|`$. Its law is unitarily
invariant: $`Ug`$ has the law of $`g`$ for every unitary $`U`$ (Proposition A.16(a)), and
$`Ug/\|Ug\|=U\psi`$. By Proposition A.22 it is the only unitarily invariant law on the unit sphere, and,
as shown in the proof there, it is the law of the first column of a Haar unitary (Theorem 5.3).

<a id="r-7-9"></a>

**Lemma 7.9.** $`\|g\|`$ and $`\psi=g/\|g\|`$ are independent.

*Proof.* (This is part of Theorem A.20(b).) Let $`h\ge0`$ and $`\phi\ge0`$ be Borel functions on
$`(0,\infty)`$ and on the unit sphere $`S^{2N-1}`$ of $`\mathbb C^N\cong\mathbb R^{2N}`$. As $`g`$ has the density
$`\pi^{-N}e^{-|x|^2}`$ (Proposition A.16(a)), polar coordinates (Theorem A.18 with $`n=2N`$ and the
measure $`\mu`$ defined there) give
```math
\mathbb E\big[h(\|g\|)\phi(\psi)\big]=\pi^{-N}\!\int_{\mathbb R^{2N}}\!h(|x|)e^{-|x|^2}\phi\Big(\frac x{|x|}\Big)dx
=2N\pi^{-N}\int_0^\infty h(t)e^{-t^2}t^{2N-1}dt\int_{S^{2N-1}}\phi\,d\mu=:\alpha(h)\beta(\phi).
```
Taking $`\phi=1`$, then $`h=1`$, then both, gives $`\mathbb Eh(\|g\|)=\alpha(h)\beta(1)`$,
$`\mathbb E\phi(\psi)=\alpha(1)\beta(\phi)`$ and $`\alpha(1)\beta(1)=1`$, hence
$`\mathbb E[h(\|g\|)\phi(\psi)]=\mathbb Eh(\|g\|)\,\mathbb E\phi(\psi)`$. $`\square`$

<a id="r-7-10"></a>

**Lemma 7.10.** $`\mathbb E\|g\|^4=N(N+1)`$. If $`G`$ is an $`a\times b`$ matrix with independent entries distributed
as the $`g_i`$ (a *Ginibre* matrix), then $`\mathbb E\mathrm{Tr}(GG^*)^2=ab(a+b)`$.

*Proof.* With $`g_i=x+iy`$ and $`x,y`$ independent $`N(0,\frac12)`$, $`\mathbb Ex^4=3(\mathbb Ex^2)^2=\frac34`$
(Theorem A.24: three pairings), so
$`\mathbb E|g_i|^4=\mathbb Ex^4+2\mathbb Ex^2\mathbb Ey^2+\mathbb Ey^4=\tfrac34+\tfrac12+\tfrac34=2`$ and
$`\mathbb E\|g\|^4=\sum_i\mathbb E|g_i|^4+\sum_{i\ne j}\mathbb E|g_i|^2\mathbb E|g_j|^2=2N+N(N-1)`$. Next,
$`\mathrm{Tr}(GG^*)^2=\sum_{i,i',j,j'}G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}`$. The entries and their
conjugates are complex linear combinations of the coordinates of one centred Gaussian vector (the
real and imaginary parts of all entries), with $`\mathbb E[G_{ij}G_{kl}]=\mathbb E[\bar G_{ij}\bar G_{kl}]=0`$ and
$`\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}`$ (Proposition A.16(a)). By the Isserlis–Wick formula
(Theorem A.24, complex form, with the three pairings of the four factors),
$`\mathbb E[G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}]=\delta_{ii'}\delta_{i'i}+0+\delta_{jj'}\delta_{jj'}=\delta_{ii'}+\delta_{jj'}`$.
Summing over $`i,i'\le a`$ and $`j,j'\le b`$ gives $`ab^2+a^2b`$. $`\square`$

<a id="r-7-11"></a>

**Theorem 7.11 (E. Lubkin, J. Math. Phys. 19, 1028, 1978).** For a Haar-random
$`\psi\in\mathbb C^a\otimes\mathbb C^b`$ and $`\rho_A=\mathrm{Tr}_2|\psi\rangle\langle\psi|`$,
```math
\mathbb E\mathrm{Tr}\,\rho_A^2=\frac{a+b}{ab+1}.
```

*Proof.* Let $`G_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,g\rangle`$, a Ginibre matrix; the matrix of $`\psi`$ is
$`\Psi=G/\|g\|`$, and by Proposition 7.2(a), $`(\rho_A)_{\alpha\alpha'}=\sum_\beta\Psi_{\alpha\beta}\bar\Psi_{\alpha'\beta}`$, i.e.
$`\rho_A=\Psi\Psi^*=GG^*/\|g\|^2`$. Hence
$`\mathrm{Tr}(GG^*)^2=\|g\|^4\mathrm{Tr}\,\rho_A^2`$, where $`\mathrm{Tr}\,\rho_A^2\in[0,1]`$ is a function of
$`\psi`$. By Lemma 7.9, $`\mathbb E\mathrm{Tr}(GG^*)^2=\mathbb E\|g\|^4\,\mathbb E\mathrm{Tr}\,\rho_A^2`$, and by Lemma 7.10
$`\mathbb E\mathrm{Tr}\,\rho_A^2=ab(a+b)/(ab(ab+1))`$. $`\square`$

Chapter 6 obtains the same mean from the fourth moments of a unitarily invariant unit vector
(Lemma 6.14), as the output purity of a random channel (Theorem 6.15 and the remark after it).

<a id="r-7-12"></a>

**Corollary 7.12.** $`\log\frac{ab+1}{a+b}\le\mathbb E\,S(\rho_A)\le\log\min(a,b)`$.

*Proof.* By Corollary 3.6, $`S(\rho_A)\le\log a`$, and by Lemma 7.7,
$`S(\rho_A)=S(\mathrm{Tr}_1|\psi\rangle\langle\psi|)\le\log b`$. The finite Jensen inequality for the concave
$`\log`$ (solution of Exercise 2 in Chapter 3) with the weights $`\lambda_i>0`$, the nonzero eigenvalues of
$`\rho_A`$, gives $`\sum_i\lambda_i\log\lambda_i\le\log\sum_i\lambda_i^2`$, i.e. $`S(\rho_A)\ge-\log X`$ with
$`X=\mathrm{Tr}\,\rho_A^2\in(0,1]`$. Let $`m=\mathbb EX=\frac{a+b}{ab+1}`$ (Theorem 7.11). Lemma 3.4 for the convex
function $`-\log`$ on $`(0,\infty)`$ at the point $`m`$ gives $`-\log X\ge-\log m-(X-m)/m`$; both sides are
bounded ($`X\ge1/a`$ by Cauchy–Schwarz), and taking expectations gives $`\mathbb E[-\log X]\ge-\log m`$. $`\square`$

For fixed $`a`$ and $`b\to\infty`$ the bounds meet: the mean entropy of a small part of a random pure
state tends to its maximum $`\log a`$. The exact mean is Page's formula, proved in Appendix C.

<a id="r-7-13"></a>

**Theorem 7.13 (Page's formula).** For $`a\le b`$,
```math
\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}.
```
Moreover $`\log a-\frac{a-1}{ab+1}-\frac{a-1}{2b}\le\mathbb E\,S(\rho_A)\le\log a-\frac{a-1}{2b}`$, so
$`\big|\mathbb E\,S(\rho_A)-\log a+\frac a{2b}\big|\le\frac1{2b}`$; for two qubits $`\mathbb E\,S(\rho_A)=\frac13`$.

*Proof.* This is Theorem C.17 together with Corollary C.18. The proof in Appendix C starts from
$`\rho_A=GG^*/\mathrm{Tr}\,GG^*`$ (proof of Theorem 7.11), computes the joint density of the eigenvalues
of the Wishart matrix $`GG^*`$ (Theorem C.13) and of $`\rho_A`$ (Theorem C.14), and integrates $`x\log x`$
against the Laguerre level density (Lemmas C.15 and C.16). $`\square`$

*Computation:* Theorems 7.11 and 7.13 and Corollary 7.12 on $`10^4`$ random pure states of $`\mathbb C^3\otimes\mathbb C^6`$ — [tutorial, Section 7.5](../tutorial/07_tensor_products.md#s-7-5).

The formula was conjectured by D. N. Page (Phys. Rev. Lett. 71, 1291, 1993) and proved by S. K.
Foong and S. Kanno (Phys. Rev. Lett. 72, 1148, 1994), J. Sánchez-Ruiz (Phys. Rev. E 52, 5653, 1995)
and S. Sen (Phys. Rev. Lett. 77, 1, 1996).

For $`a>b`$ the roles of $`a`$ and $`b`$ are exchanged. Indeed, the unitary $`e_\alpha\otimes e_\beta\mapsto e_\beta\otimes e_\alpha`$
of $`\mathbb C^a\otimes\mathbb C^b`$ onto $`\mathbb C^b\otimes\mathbb C^a`$ permutes the coordinates of $`g`$, so it maps a
Haar-random pure state $`\psi`$ to a Haar-random pure state $`\psi'`$ of $`\mathbb C^b\otimes\mathbb C^a`$; its matrix is
$`\Psi^T`$, so $`\mathrm{Tr}_2|\psi'\rangle\langle\psi'|=\Psi^T\bar\Psi=\mathrm{Tr}_1|\psi\rangle\langle\psi|`$ (Proposition 7.2(a)),
whose entropy is $`S(\rho_A)`$ by Lemma 7.7. Hence, with $`a'=\min(a,b)`$ and $`b'=\max(a,b)`$, Theorem 7.13 says
that the mean entanglement entropy falls short of its maximum $`\log a'`$ by at least $`\frac{a'-1}{2b'}`$ and
by at most $`\frac{a'-1}{2b'}+\frac{a'-1}{a'b'+1}`$, which is less than one nat: it is $`0`$ for $`b'=1`$, and
otherwise $`\frac{a'-1}{2b'}\le\frac12`$ and $`\frac{a'-1}{a'b'+1}<\frac1{b'}\le\frac12`$. So for fixed $`ab`$, as a
function of $`a`$, the mean is within one nat of $`\log\min(a,b)=\min(\log a,\log(ab)-\log a)`$.

The law of $`\rho_A`$ is the **induced measure** (the name is that of K. Życzkowski and H.-J. Sommers,
J. Phys. A 34, 7111, 2001; Theorem C.14 gives the density of its eigenvalues).
`states.random_density_matrix(alg, rank=r)` forms $`G_cG_c^*`$ with independent $`k_c\times r`$ Ginibre
blocks $`G_c`$ and divides by $`\sum_c\mathrm{Tr}\,G_cG_c^*`$. For one block this is the law of $`\rho_A`$
above with $`b=r`$ (proof of Theorem 7.11). For several blocks it is the law of the density of the
restriction to $`\mathcal M_1`$ of the vector state of a Haar-random pure state $`\psi=g/\|g\|`$ of
$`H_1\otimes\mathbb C^r`$. Indeed, let $`G_{i\beta}=\langle e_i\otimes e_\beta,g\rangle`$, a $`D\times r`$ Ginibre
matrix ($`D=\sum_ck_c`$), and $`G_c`$ its rows in block $`c`$: independent $`k_c\times r`$ Ginibre matrices.
Take $`\mathcal M_2=M_r`$ (one sector) and a constant fusion rule; then there is a single total charge,
$`\Pi=1`$, and $`\mathcal M_1\otimes_f\mathcal M_2=B(H_1\otimes\mathbb C^r)`$. By Theorem 7.6 the restriction has
the density $`E_1(\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|)`$, and $`\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|=GG^*/\|g\|^2`$
(Proposition 7.2(a), as in the proof of Theorem 7.11), whose diagonal blocks are $`G_cG_c^*/\|g\|^2`$,
with $`\|g\|^2=\sum_c\mathrm{Tr}\,G_cG_c^*`$.

*In the library:* `states.random_density_matrix(alg, rank=r)` — [tutorial, Section 7.5](../tutorial/07_tensor_products.md#s-7-5).

<a id="r-7-14"></a>

**Proposition 7.14.** For $`\rho`$ drawn by `random_density_matrix(alg, rank=r)` on $`\bigoplus_cM_{k_c}`$,
$`D=\sum_ck_c`$,
```math
\mathbb E\mathrm{Tr}\,\rho^2=\frac{\sum_ck_c(k_c+r)}{D(Dr+1)} .
```

*Proof.* $`\mathrm{Tr}\,\rho^2=\sum_c\mathrm{Tr}(G_cG_c^*)^2/\|g\|^4`$ with $`g\in\mathbb C^{Dr}`$ the collection of
all entries, a function of $`g/\|g\|`$ with values in $`[0,1]`$. As in Theorem 7.11, Lemma 7.9 gives
$`\sum_c\mathbb E\mathrm{Tr}(G_cG_c^*)^2=\mathbb E\|g\|^4\,\mathbb E\mathrm{Tr}\,\rho^2`$, and Lemma 7.10 gives
$`\mathbb E\mathrm{Tr}(G_cG_c^*)^2=k_cr(k_c+r)`$ and $`\mathbb E\|g\|^4=Dr(Dr+1)`$, so
$`\mathbb E\mathrm{Tr}\,\rho^2=\sum_ck_cr(k_c+r)/(Dr(Dr+1))`$. $`\square`$

*Computation:* Proposition 7.14 for $`M_1\oplus M_2\oplus M_3`$, $`r=2`$ — [tutorial, Section 7.5](../tutorial/07_tensor_products.md#s-7-5).

> **Physical meaning** (physical interpretation, not a theorem of this course). The mean
> entanglement entropy of a random pure state, plotted against the size of the part at fixed total
> dimension (the curve described after Theorem 7.13), is called the *Page curve*. In physics it
> is used to describe how information leaves a black hole in the Hawking radiation, and why the
> parts of a closed system look thermal; these are interpretations, not results of this course.

<a id="s-7-6"></a>

## 7.6 What to compare with

*In the library:* comparison with other software — [tutorial, Section 7.6](../tutorial/07_tensor_products.md#s-7-6).

<a id="x-7"></a>

## Exercises

1. Compute the dual of $`\mathrm{Tr}_2`$ with respect to $`\tau_{\mathrm{vN}}`$ on both
   $`\mathcal M_1\otimes\mathcal M_2`$ and $`\mathcal M_1`$.
2. (Subadditivity.) Let $`\rho`$ be a state on $`\mathcal M_1\otimes\mathcal M_2`$ with invertible reduced
   states $`\rho_1,\rho_2`$. Show that $`S(\rho)\le S(\rho_1)+S(\rho_2)`$, with equality iff
   $`\rho=\rho_1\otimes\rho_2`$.
3. Let `site` be $`\mathbb C\oplus\mathbb C`$ with charges $`0,1`$ (an empty or occupied site), fused with
   itself by addition. Compute $`E(\psi)`$ and its two parts for $`\psi=(|01\rangle+|10\rangle)/\sqrt2`$, one
   particle shared by two sites, and describe the restriction of the vector state of
   $`(|00\rangle+|11\rangle)/\sqrt2`$ to the fused algebra.

*Computation:* the three exercises numerically — [tutorial, Exercises](../tutorial/07_tensor_products.md#x-7).

<details>
<summary>Solutions</summary>

1. The weights are $`1/(C_1k_c)`$ on $`\mathcal M_1`$ and $`1/(C_1C_2k_cl_d)`$ on the product. For
   $`A\in\mathcal M_1`$, by Proposition 7.2(a),
   $`\tau_{\mathrm{vN}}((\mathrm{Tr}_2\,X)A)=\sum_c\frac1{C_1k_c}\sum_d\mathrm{Tr}(X_{cd}(A_c\otimes1))    =\sum_{c,d}\frac1{C_1C_2k_cl_d}\mathrm{Tr}\big(X_{cd}\,C_2l_d(A_c\otimes1)\big)`$, so the dual is
   $`A\mapsto\bigoplus_{c,d}C_2l_d\,A_c\otimes1_{l_d}=A\otimes Z`$ with the central element
   $`Z=\bigoplus_dC_2l_d1_{l_d}`$ of $`\mathcal M_2`$. It is unital only if $`Z=1`$, i.e. $`\mathcal M_2=\mathbb C`$,
   and not in the example ($`Z=4\cdot1_2\oplus6\cdot1_3`$). For `Tr_norm` the same computation gives
   $`A\otimes\bigoplus_dl_d1_{l_d}`$, unital only if all $`l_d=1`$. The blunt traces always make
   $`\mathrm{Tr}_2`$ the dual of $`A\mapsto A\otimes1`$; in general, $`\mathrm{tr}_w`$ on $`\mathcal M_1`$ and the
   weights $`w_{cd}=w_c`$ on the product do, by the same computation.

2. $`\rho_1\otimes\rho_2`$ is a state (Theorem 7.1(c)), invertible, with eigenvectors $`u_i\otimes v_j`$ and
   eigenvalues $`\lambda_i\mu_j`$, so $`\log(\rho_1\otimes\rho_2)=\log\rho_1\otimes1+1\otimes\log\rho_2`$
   (Theorem 1.5(a)). By the definition of the partial trace
   $`\mathrm{Tr}\,\rho(\log\rho_1\otimes1)=\mathrm{Tr}\,\rho_1\log\rho_1`$, and likewise for the second
   factor. Hence $`S(\rho\|\rho_1\otimes\rho_2)=S(\rho_1)+S(\rho_2)-S(\rho)`$, which is $`\ge0`$ with
   equality iff $`\rho=\rho_1\otimes\rho_2`$ (Corollary 3.6). This difference is the mutual information.

3. The sectors are $`N=0,1,2`$, of sizes $`1,2,1`$. $`\psi`$ lies in $`N=1`$, with $`\psi_0=|01\rangle/\sqrt2`$,
   $`\psi_1=|10\rangle/\sqrt2`$, so $`p=(\frac12,\frac12)`$ and, the blocks being one-dimensional,
   $`E(\hat\psi_c)=0`$: by Proposition 7.8, $`E(\psi)=H(p)=\log2`$ is all number entanglement, and the
   configurational part is $`0`$. The vector $`(|00\rangle+|11\rangle)/\sqrt2`$ superposes $`N=0`$ and
   $`N=2`$; on the block-diagonal fused algebra its vector state is the mixture
   $`\frac12|00\rangle\langle00|+\frac12|11\rangle\langle11|`$ of product states (Chapter 3, Exercise 3).

</details>

**Next:** [Chapter 8](08_dynamics.md) treats dynamics: the Schrödinger and Lindblad equations in an
algebra with sectors, jumps between sectors and post-selection.
