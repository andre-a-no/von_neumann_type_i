**English** · [Español](es/06_quantum_channels.md) · [Français](fr/06_quantum_channels.md) · [Deutsch](de/06_quantum_channels.md) · [中文](zh/06_quantum_channels.md) · [日本語](ja/06_quantum_channels.md) · [Русский](ru/06_quantum_channels.md)

# 6. Quantum channels

A quantum channel is a completely positive, trace-preserving linear map. Conjugation by a unitary,
the conditional expectation onto the centre, the map $`\rho\mapsto\mathrm{Tr}_E(V\rho V^*)`$ given by an
isometry $`V`$, and maps that move weight from one sector to another are examples; their physical
reading (an evolution, a measurement whose result is forgotten, an environment, the loss of a
particle) is discussed only in the *Physical meaning* boxes. This chapter
proves the standard descriptions of such maps — Choi matrices, Kraus operators, the Stinespring
dilation — for maps between finite-dimensional von Neumann algebras, and then studies what is
special about direct sums: maps between sectors, duals with respect to the traces of Chapter 3,
and the classical part of a channel. Chapters 2 and 3 already used the library's channels; here
is their theory.

**You will learn** why positivity is not enough, how to pass between Choi matrices, Kraus
operators, superoperators and isometries, why "trace preserving" depends on the trace once a map
changes the sector, which operators a unital channel leaves fixed, and what a random channel does
to a pure state on average.

The numerical checks of this chapter are in [Tutorial 6](../tutorial/06_quantum_channels.md).

Throughout, $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$ and
$`\mathcal N=\bigoplus_{d=1}^{C'}M_{k'_d}(\mathbb C)`$ with minimal central projections $`P_c`$ and
$`Q_d`$; $`A=\bigoplus_cA_c`$; $`M_{k'\times k}`$ is the space of complex $`k'\times k`$ matrices;
$`E_{ij}=|i\rangle\langle j|`$ are matrix units; $`\mathrm{Tr}`$ is the blunt trace
$`\sum_c\mathrm{Tr}\,A_c`$, and $`\mathrm{tr}_w=\sum_cw_c\mathrm{Tr}\,A_c`$ a faithful
trace (Theorem 3.1). Block matrices act on $`\mathbb C^n\otimes\mathbb C^k`$, block index first.

<a id="s-6-1"></a>

## 6.1 Completely positive maps and Choi's theorem

A linear map $`\Phi:\mathcal M\to\mathcal N`$ is **positive** if $`A\ge0`$ implies $`\Phi(A)\ge0`$.
The $`n\times n`$ matrices $`[A_{ij}]`$ with entries in $`\mathcal M`$ form again a direct sum of matrix
algebras, $`M_n(\mathcal M)\cong\bigoplus_cM_n(M_{k_c})`$ via $`[A_{ij}]\mapsto([(A_{ij})_c]_{ij})_c`$,
so $`[A_{ij}]\ge0`$ means that every block matrix $`[(A_{ij})_c]\in M_{nk_c}`$ is positive. $`\Phi`$ is
**$`n`$-positive** if $`\Phi_n:[A_{ij}]\mapsto[\Phi(A_{ij})]`$ is positive, and **completely positive**
(CP) if it is $`n`$-positive for every $`n`$. It is **trace preserving** (TP) if
$`\mathrm{Tr}\,\Phi(A)=\mathrm{Tr}\,A`$, **unital** if $`\Phi(1)=1`$, and a **channel** if it is
CP and TP. We use repeatedly that $`X\ge0`$ implies $`TXT^*\ge0`$ for every rectangular $`T`$, since
$`\langle y,TXT^*y\rangle=\langle T^*y,XT^*y\rangle`$.

<a id="r-6-1"></a>

**Proposition 6.1.** For $`k\ge2`$, the transpose $`\theta(A)=A^T`$ on $`M_k`$ is positive but not
$`2`$-positive.

*Proof.* For $`A\ge0`$, $`\langle x,A^Tx\rangle=\sum_{ij}\bar x_iA_{ji}x_j=\langle\bar x,A\bar x\rangle\ge0`$.
With $`\Omega=e_1\otimes e_1+e_2\otimes e_2\in\mathbb C^2\otimes\mathbb C^k`$,
$`|\Omega\rangle\langle\Omega|=[E_{ij}]_{i,j=1}^2\ge0`$, and $`\theta_2`$ maps it to
$`F=[E_{ji}]_{i,j=1}^2`$, which sends $`e_1\otimes e_2-e_2\otimes e_1`$ to its negative. $`\square`$

*Computation:* the Choi matrix of the transpose has the eigenvalue $`-1`$ — [tutorial, Section 6.1](../tutorial/06_quantum_channels.md#s-6-1).

<a id="r-6-2"></a>

**Theorem 6.2 (Choi).** For a linear map $`\Phi:M_k\to M_{k'}`$ the following are equivalent:
(a) $`\Phi`$ is CP; (b) $`\Phi`$ is $`k`$-positive; (c) the **Choi matrix**
```math
J(\Phi)=\sum_{i,j=1}^kE_{ij}\otimes\Phi(E_{ij})=[\Phi(E_{ij})]_{i,j=1}^k\in M_k(M_{k'})
```
is positive; (d) there are $`K_1,\dots,K_r\in M_{k'\times k}`$ (**Kraus operators**) with
$`\Phi(A)=\sum_lK_lAK_l^*`$, and one can take $`r=\mathrm{rank}\,J(\Phi)\le kk'`$.

*Proof.* (a)$`\Rightarrow`$(b) is trivial. (b)$`\Rightarrow`$(c): for
$`\Omega=\sum_ie_i\otimes e_i\in\mathbb C^k\otimes\mathbb C^k`$,
$`|\Omega\rangle\langle\Omega|=[E_{ij}]\ge0`$ and $`J(\Phi)=\Phi_k([E_{ij}])`$. (c)$`\Rightarrow`$(d): by
Theorem 1.3, $`J(\Phi)=\sum_{l=1}^rv_lv_l^*`$ with $`r=\mathrm{rank}\,J(\Phi)`$ (the eigenvectors of
the nonzero eigenvalues, scaled by their square roots). Every $`v\in\mathbb C^k\otimes\mathbb C^{k'}`$
is $`v=\sum_ie_i\otimes Ke_i`$ for exactly one $`K\in M_{k'\times k}`$ ($`Ke_i`$ is the $`i`$-th block of
$`v`$), and the $`(i,j)`$ block of $`vv^*`$ is $`Ke_ie_j^*K^*=KE_{ij}K^*`$. Writing $`v_l`$ in this way with
$`K_l`$ and comparing blocks gives $`\Phi(E_{ij})=\sum_lK_lE_{ij}K_l^*`$, hence (d) by linearity.
(d)$`\Rightarrow`$(a): $`\Phi_n([A_{ij}])=\sum_l(1_n\otimes K_l)[A_{ij}](1_n\otimes K_l)^*\ge0`$. $`\square`$

<a id="r-6-3"></a>

**Corollary 6.3.** If $`\Phi(A)=\sum_{l=1}^sL_lAL_l^*`$, then $`s\ge\mathrm{rank}\,J(\Phi)`$. So the
minimal number of Kraus operators, the **Kraus rank**, is $`\mathrm{rank}\,J(\Phi)`$.

*Proof.* Reading the last computation backwards, $`J(\Phi)=\sum_{l=1}^sw_lw_l^*`$ with
$`w_l=\sum_ie_i\otimes L_le_i`$, a sum of $`s`$ matrices of rank $`\le1`$. $`\square`$

*Computation:* the Choi matrix and the Kraus rank of a random channel — [tutorial, Section 6.1](../tutorial/06_quantum_channels.md#s-6-1).

A map between direct sums is a matrix of maps between blocks. Let $`\iota_c:M_{k_c}\to\mathcal M`$
put a matrix into block $`c`$, $`\pi_d:\mathcal N\to M_{k'_d}`$ take block $`d`$, and
$`\Phi_{dc}=\pi_d\circ\Phi\circ\iota_c:M_{k_c}\to M_{k'_d}`$. Since $`A=\sum_c\iota_c(A_c)`$,
$`\Phi(A)_d=\sum_c\Phi_{dc}(A_c)`$.

<a id="r-6-4"></a>

**Proposition 6.4.** $`\Phi:\mathcal M\to\mathcal N`$ is $`n`$-positive iff every $`\Phi_{dc}`$ is. Hence
$`\Phi`$ is CP iff there are $`K^{dc}_l\in M_{k'_d\times k_c}`$, $`l=1,\dots,r`$, with
```math
\Phi(A)_d=\sum_c\sum_{l=1}^rK^{dc}_l\,A_c\,K^{dc*}_l ,
```
and $`r=\max_{d,c}\mathrm{rank}\,J(\Phi_{dc})`$ suffices. $`\Phi`$ is TP iff
$`\sum_{d,l}K^{dc*}_lK^{dc}_l=1_{k_c}`$ for every $`c`$.

*Proof.* If all $`\Phi_{dc}`$ are $`n`$-positive and $`[A_{ij}]\ge0`$, then each $`X^c=[(A_{ij})_c]\ge0`$, and
block $`d`$ of $`\Phi_n([A_{ij}])`$ is $`\sum_c(\Phi_{dc})_n(X^c)\ge0`$. Conversely, for $`[X_{ij}]\ge0`$ in
$`M_n(M_{k_c})`$, $`[\iota_c(X_{ij})]\ge0`$ in $`M_n(\mathcal M)`$ (its block $`c`$ is $`[X_{ij}]`$, the others
vanish), and block $`d`$ of its image is $`[\Phi_{dc}(X_{ij})]`$. The Kraus form follows from Theorem 6.2
for each pair $`(d,c)`$, padded with zero operators to a common $`r`$. Finally, by cyclicity,
$`\mathrm{Tr}\,\Phi(A)=\sum_c\mathrm{Tr}\big(A_c\sum_{d,l}K^{dc*}_lK^{dc}_l\big)`$, which equals
$`\sum_c\mathrm{Tr}\,A_c`$ for all $`A`$ iff each sum is $`1_{k_c}`$ (take matrix units:
$`\mathrm{Tr}(E_{ij}X)=X_{ji}`$). $`\square`$

A linear map $`\Phi:M_k\to M_{k'}`$ satisfies $`\Phi(X^*)=\Phi(X)^*`$ for all $`X`$ iff $`J(\Phi)^*=J(\Phi)`$: the
$`(i,j)`$ block of $`J^*`$ is $`\Phi(E_{ji})^*`$, which equals $`\Phi(E_{ij})`$ for all $`i,j`$ iff $`\Phi`$ preserves
adjoints, by linearity.

*In the library:* sector-preserving maps are a `Channel`, the general maps of Proposition 6.4 an `InterSectorChannel`; `Channel.choi()`, `Channel.superoperator()` and `Channel.from_superoperator` pass between the representations — [tutorial, Section 6.1](../tutorial/06_quantum_channels.md#s-6-1).

> **Physical meaning** (physical interpretation, not a theorem of this course). A system is never alone: it may be correlated with an untouched ancilla
> with $`n`$ levels. A physical process acts as $`\Phi_n=\mathrm{id}_n\otimes\Phi`$ on the pair and must
> map states to states, so $`\Phi`$ must be completely positive. The transpose is positive, but it maps
> the entangled state $`|\Omega\rangle\langle\Omega|/2`$ to an operator with eigenvalue $`-1/2`$. The
> Choi matrix is the image of the maximally entangled state: one state decides complete positivity.

<a id="s-6-2"></a>

## 6.2 Stinespring dilation and the freedom in the Kraus operators

For $`\Phi(A)=\sum_{l=1}^rK_lAK_l^*`$ from $`M_k`$ to $`M_{k'}`$ put $`\Phi^*(B)=\sum_lK_l^*BK_l`$; by
cyclicity $`\mathrm{Tr}(\Phi(A)B)=\mathrm{Tr}(A\Phi^*(B))`$, so $`\Phi^*`$ is the dual of
$`\Phi`$ for the trace (Section 6.3) and does not depend on the choice of the $`K_l`$ (if
$`\mathrm{Tr}(AX)=\mathrm{Tr}(AX')`$ for all $`A`$, then $`X=X'`$: take $`A=E_{ji}`$). On
$`\mathbb C^r\otimes\mathbb C^{k'}`$ (the **environment** $`\mathbb C^r`$ first) let
$`\mathrm{Tr}_E\,X=\sum_l(e_l^*\otimes1)X(e_l\otimes1)`$.

<a id="r-6-5"></a>

**Theorem 6.5 (Stinespring, finite dimensions).** Let
$`V=\sum_le_l\otimes K_l:\mathbb C^k\to\mathbb C^r\otimes\mathbb C^{k'}`$. Then
```math
\Phi(A)=\mathrm{Tr}_E(VAV^*),\qquad\Phi^*(B)=V^*(1_r\otimes B)V,
```
and $`\Phi`$ is TP iff $`V`$ is an isometry. Conversely, every $`V`$ defines a CP map by the first formula.

*Proof.* $`(e_l^*\otimes1)V=K_l`$, so $`\mathrm{Tr}_E(VAV^*)=\sum_lK_lAK_l^*`$; conversely
$`\mathrm{Tr}_E(VAV^*)`$ is a Kraus form with $`K_l=(e_l^*\otimes1)V`$. Next
$`V^*(1\otimes B)V=\sum_{l,l'}\langle e_l,e_{l'}\rangle K_l^*BK_{l'}=\Phi^*(B)`$, and
$`V^*V=\sum_lK_l^*K_l`$, which is $`1`$ iff $`\Phi`$ is TP (Proposition 6.4). $`\square`$

So the smallest $`r`$ for which $`\Phi(A)=\mathrm{Tr}_E(VAV^*)`$ with some
$`V:\mathbb C^k\to\mathbb C^r\otimes\mathbb C^{k'}`$ is the Kraus rank $`\mathrm{rank}\,J(\Phi)`$: every
such $`V`$ gives a Kraus form with $`r`$ operators, so $`r\ge\mathrm{rank}\,J(\Phi)`$ (Corollary 6.3), and
$`r=\mathrm{rank}\,J(\Phi)`$ is attained by Theorem 6.2(d).

*Computation:* Theorem 6.5 for a random channel — [tutorial, Section 6.2](../tutorial/06_quantum_channels.md#s-6-2).

The formula $`\Phi^*(B)=V^*\pi(B)V`$ with the representation $`\pi(B)=1_r\otimes B`$ is the form in
which W. F. Stinespring (Proc. AMS 6 (1955) 211) stated the theorem, in a far more general setting
than the finite-dimensional one proved here (historical attribution; further reading: V. Paulsen,
*Completely Bounded Maps and Operator Algebras*, Cambridge 2002, Ch. 4). For a map between direct
sums, let $`\iota'_d`$ be the inclusion of the $`d`$-th summand of $`\bigoplus_d\mathbb C^{k'_d}`$ and
$`V_c=\sum_{d,l}e_l\otimes\iota'_dK^{dc}_l:\mathbb C^{k_c}\to\mathbb C^r\otimes\bigoplus_d\mathbb C^{k'_d}`$.
Then $`(e_l^*\otimes1)V_c=\sum_d\iota'_dK^{dc}_l`$, so the $`(d,d)`$ block of
$`\sum_c\mathrm{Tr}_E(V_cA_cV_c^*)=\sum_{c,l}\big(\sum_d\iota'_dK^{dc}_l\big)A_c\big(\sum_{d'}\iota'_{d'}K^{d'c}_l\big)^*`$
is $`\sum_{c,l}K^{dc}_lA_cK^{dc*}_l=\Phi(A)_d`$; and, as $`\iota'^*_d\iota'_{d'}=\delta_{dd'}1`$,
$`V_c^*V_c=\sum_{d,l}K^{dc*}_lK^{dc}_l`$, so $`\Phi`$ is TP iff every $`V_c`$ is an isometry
(Proposition 6.4). This is how `random_channel` and `random_inter_sector_channel` draw channels
(Section 6.5).

*In the library:* `channels.random_channel` and `channels.random_inter_sector_channel` — [tutorial, Section 6.5](../tutorial/06_quantum_channels.md#s-6-5).

<a id="r-6-6"></a>

**Lemma 6.6.** If $`X,Y\in M_{p\times n}`$ and $`XX^*=YY^*`$, then $`Y=XW`$ for a unitary $`W\in M_n`$.

*Proof.* $`\|X^*z\|^2=\langle z,XX^*z\rangle=\|Y^*z\|^2`$ for all $`z`$, so $`X^*z\mapsto Y^*z`$ is a
well-defined linear isometry from the range of $`X^*`$ onto the range of $`Y^*`$. As in the proof of
Theorem 1.8 it extends to a unitary $`U`$ with $`UX^*=Y^*`$; take $`W=U^*`$. $`\square`$

<a id="r-6-7"></a>

**Theorem 6.7 (unitary freedom).** Let $`K_1,\dots,K_n`$ and $`L_1,\dots,L_n`$ be in $`M_{k'\times k}`$
(pad the shorter family with zeros). Then $`\sum_lK_lAK_l^*=\sum_jL_jAL_j^*`$ for all $`A`$ iff
$`L_j=\sum_lu_{jl}K_l`$ for a unitary $`u\in M_n`$.

*Proof.* If $`L_j=\sum_lu_{jl}K_l`$, then
$`\sum_jL_jAL_j^*=\sum_{l,l'}(u^*u)_{l'l}K_lAK_{l'}^*=\sum_lK_lAK_l^*`$. Conversely, equal maps have
equal Choi matrices, $`\sum_lv_lv_l^*=\sum_jw_jw_j^*`$ with $`v_l=\sum_ie_i\otimes K_le_i`$ and
$`w_j=\sum_ie_i\otimes L_je_i`$. Lemma 6.6, applied to the matrices with columns $`v_l`$ and $`w_j`$, gives
$`w_j=\sum_lW_{lj}v_l`$; as $`K\mapsto\sum_ie_i\otimes Ke_i`$ is linear and injective,
$`L_j=\sum_lW_{lj}K_l`$, and $`u=W^T`$ is unitary. $`\square`$

*Computation:* Theorem 6.7 with $`n=3`$ — [tutorial, Section 6.2](../tutorial/06_quantum_channels.md#s-6-2).

> **Physical meaning** (physical interpretation, not a theorem of this course). Every channel is an
> isometric interaction of the system with an environment
> in a fixed state, followed by ignoring the environment. Measuring the environment in the basis
> $`e_l`$ gives outcome $`l`$ with probability $`\mathrm{Tr}\,K_l\rho K_l^*`$ and leaves the state
> $`K_l\rho K_l^*/\mathrm{Tr}(K_l\rho K_l^*)`$: the Kraus operators are the possible "quantum jumps".
> Another basis gives other jumps but the same channel (Theorem 6.7), and the smallest environment
> has dimension equal to the Kraus rank, as shown after Theorem 6.5.

<a id="s-6-3"></a>

## 6.3 Duals with respect to a trace

**Definition.** Let $`\mathrm{tr}_w`$, $`\mathrm{tr}_v`$ be faithful traces on $`\mathcal M`$,
$`\mathcal N`$. The **dual** of a linear map $`\Phi:\mathcal M\to\mathcal N`$ with respect to $`(w,v)`$ is
the linear map $`\Phi^\sharp:\mathcal N\to\mathcal M`$ with
```math
\mathrm{tr}_v\big(\Phi(x)\,y\big)=\mathrm{tr}_w\big(x\,\Phi^\sharp(y)\big)\qquad(x\in\mathcal M,\ y\in\mathcal N).
```
It exists and is unique: for fixed $`y`$, $`x\mapsto\mathrm{tr}_v(\Phi(x)y)`$ is a linear functional,
so by Theorem 3.2 it equals $`\mathrm{tr}_w(\rho x)=\mathrm{tr}_w(x\rho)`$ for exactly one
$`\rho=:\Phi^\sharp(y)`$, and uniqueness makes $`\Phi^\sharp`$ linear. For blunt traces we write $`\Phi^*`$.
Let $`z_w=\sum_cw_cP_c`$ and $`z_v=\sum_dv_dQ_d`$ (central and invertible).

<a id="r-6-8"></a>

**Proposition 6.8.** If $`\Phi`$ has Kraus operators $`K^{dc}_l`$ (Proposition 6.4), then
```math
\Phi^\sharp(y)_c=\sum_d\frac{v_d}{w_c}\sum_lK^{dc*}_l\,y_d\,K^{dc}_l,\qquad\text{i.e.}\qquad
\Phi^\sharp(y)=z_w^{-1}\,\Phi^*(z_v\,y).
```
In particular $`\Phi^\sharp`$ is CP, with Kraus operators $`(v_d/w_c)^{1/2}K^{dc*}_l`$.

*Proof.* By cyclicity of the trace in each block,
```math
\mathrm{tr}_v(\Phi(x)y)=\sum_dv_d\sum_{c,l}\mathrm{Tr}\big(K^{dc}_lx_cK^{dc*}_ly_d\big)
=\sum_cw_c\mathrm{Tr}\Big(x_c\sum_d\frac{v_d}{w_c}\sum_lK^{dc*}_ly_dK^{dc}_l\Big),
```
and uniqueness gives the formula. For $`v=w=1`$ it reads $`\Phi^*(y)_c=\sum_{d,l}K^{dc*}_ly_dK^{dc}_l`$;
the factors $`v_d`$ and $`1/w_c`$ are multiplication by $`z_v`$ and $`z_w^{-1}`$. $`\square`$

$`\Phi`$ is **TP with respect to $`(w,v)`$** if $`\mathrm{tr}_v\,\Phi(x)=\mathrm{tr}_w(x)`$ for all $`x`$.

<a id="r-6-9"></a>

**Corollary 6.9.** (a) $`\Phi`$ is TP with respect to $`(w,v)`$ iff $`\Phi^\sharp(1)=1`$, i.e.
$`\sum_{d,l}(v_d/w_c)K^{dc*}_lK^{dc}_l=1_{k_c}`$ for every $`c`$. (b) $`\Phi`$ is unital iff $`\Phi^\sharp`$ is
TP with respect to $`(v,w)`$. (c) The dual of $`\Phi^\sharp`$ with respect to $`(v,w)`$ is $`\Phi`$. (d) If
$`\Phi:\mathcal M\to\mathcal M`$ is sector preserving and $`v=w`$, then $`\Phi^\sharp=\Phi^*`$ for every $`w`$,
and $`\Phi`$ is TP for one trace iff it is TP for all.

*Proof.* (a) $`\mathrm{tr}_v\,\Phi(x)=\mathrm{tr}_w(x\Phi^\sharp(1))`$, and
$`\mathrm{tr}_w(x(\Phi^\sharp(1)-1))=0`$ for all $`x`$ iff $`\Phi^\sharp(1)=1`$ (take
$`x=(\Phi^\sharp(1)-1)^*`$; $`\mathrm{tr}_w`$ is faithful); then use Proposition 6.8. (b)
$`\mathrm{tr}_w\,\Phi^\sharp(y)=\mathrm{tr}_v(\Phi(1)y)`$ equals $`\mathrm{tr}_v(y)`$ for all
$`y`$ iff $`\Phi(1)=1`$. (c) $`\mathrm{tr}_w(\Phi^\sharp(y)x)=\mathrm{tr}_v(y\Phi(x))`$ by the trace
property; use uniqueness. (d) Only the pairs $`d=c`$ occur, where $`v_d/w_c=1`$. $`\square`$

**Densities in the Heisenberg picture.** Let $`\Psi:\mathcal N\to\mathcal M`$ be linear, $`\omega`$ a linear
functional on $`\mathcal M`$ with density $`\rho`$ for $`\mathrm{tr}_w`$ (Theorem 3.2), and
$`\Psi^\sharp:\mathcal M\to\mathcal N`$ the dual of $`\Psi`$ with respect to $`(v,w)`$. Then
$`\omega(\Psi(y))=\mathrm{tr}_w(\Psi(y)\rho)=\mathrm{tr}_v(y\Psi^\sharp(\rho))=\mathrm{tr}_v(\Psi^\sharp(\rho)y)`$,
so $`\omega\circ\Psi`$ has density $`\Psi^\sharp(\rho)`$ for $`\mathrm{tr}_v`$; if $`\Psi`$ is unital, $`\Psi^\sharp`$
is TP with respect to $`(w,v)`$ by Corollary 6.9(b) (with the roles of $`\mathcal M`$ and $`\mathcal N`$
exchanged).

For maps that change the sector, Corollary 6.9(d) fails. Let $`\mathcal M=\mathcal N=M_1\oplus M_2`$ (sectors $`0`$ and $`1`$) and let $`\Phi`$ be
the map called **particle loss** (see the box below),
$`\Phi(x)=(x_0+\mathrm{Tr}\,x_1)\oplus0`$, with Kraus operators $`K^{00}=1`$, $`K^{01}_l=\langle l|`$.
It is TP for the blunt traces. For $`\mathrm{tr}_w`$ with $`w=v=(1,\tfrac12)`$, $`\Phi^\sharp(1)=1\oplus2\cdot1_2`$
and $`\mathrm{tr}_v\,\Phi(x)=x_0+\mathrm{Tr}\,x_1\ne x_0+\tfrac12\mathrm{Tr}\,x_1=\mathrm{tr}_w(x)`$.

*In the library:* `Channel.adjoint()` is $`\Phi^*`$ and `InterSectorChannel.adjoint(trace=...)` is $`\Phi^\sharp`$; checks of Proposition 6.8, Corollary 6.9(a) and particle loss — [tutorial, Section 6.3](../tutorial/06_quantum_channels.md#s-6-3).

> **Physical meaning** (physical interpretation, not a theorem of this course). The trace-independent object is the Heisenberg map
> $`\Psi:\mathcal N\to\mathcal M`$ on observables: a state $`\omega`$ before the process becomes
> $`\omega\circ\Psi`$ after it, and $`\Psi`$ is CP and unital. Densities depend on the trace
> (Theorem 3.2), and the densities transform by the dual of $`\Psi`$ for the chosen traces (the
> paragraph after Corollary 6.9). If the algebras act with multiplicities, the physical trace has
> weights $`w_c=m_c`$ (Section 2.3). In the example, sector $`0`$ is read as "no particle" and sector $`1`$ as
> "one particle on two sites"; the map then removes the particle. It conserves probability for the
> blunt trace, which counts states, not for `Tr_norm`, which weighs sectors equally.

<a id="s-6-4"></a>

## 6.4 Decoherence

The library's standard channels act inside each block ($`p,\gamma\in[0,1]`$):

| channel | $`\Phi(\rho)_c`$ | Kraus operators in block $`c`$ |
|---|---|---|
| `dephasing_channel` | $`(1-p)\rho_c+p\mathrm{diag}\,\rho_c`$ | $`\sqrt{1-p}\,1`$, $`\sqrt p\,E_{ii}`$ |
| `depolarizing_channel` | $`(1-p)\rho_c+p\mathrm{Tr}(\rho_c)\,1/k_c`$ | $`\sqrt{1-p}\,1`$, $`\sqrt{p/k_c}\,E_{ij}`$ |
| `center_expectation` | $`\mathrm{Tr}(\rho_c)\,1/k_c`$ | $`\sqrt{1/k_c}\,E_{ij}`$ |
| `amplitude_damping_channel` | levels $`i>0`$ decay to $`0`$ | $`E_{00}+\sqrt{1-\gamma}\sum_{i>0}E_{ii}`$, $`\sqrt\gamma\,E_{0i}`$ |
| `random_mixed_unitary_channel` | $`\sum_jp_jU_j\rho U_j^*`$ | $`\sqrt{p_j}\,U_j`$ |

The Kraus forms are checked by direct computation, e.g.
$`\sum_{ij}\frac1kE_{ij}\rho E_{ji}=\frac1k\sum_{ij}\rho_{jj}E_{ii}`$. All five are TP, and all but
amplitude damping with $`\gamma>0`$ (when some $`k_c\ge2`$) are unital (Proposition 6.4: compute $`\sum_lK_l^*K_l`$ and
$`\sum_lK_lK_l^*`$; for amplitude damping in a block of size $`k`$ they are $`1`$ and
$`1+(k-1)\gamma E_{00}-\gamma\sum_{i>0}E_{ii}`$, and $`\sum_jp_j=1`$ for the mixed unitary channel).
Explicitly, amplitude damping maps $`\rho_{00}\mapsto\rho_{00}+\gamma\sum_{i>0}\rho_{ii}`$,
$`\rho_{ii}\mapsto(1-\gamma)\rho_{ii}`$ and $`\rho_{ij}\mapsto(1-\gamma)\rho_{ij}`$ for $`i\ne j`$ both
positive, and $`\rho_{0j}\mapsto\sqrt{1-\gamma}\,\rho_{0j}`$, $`\rho_{j0}\mapsto\sqrt{1-\gamma}\,\rho_{j0}`$ for $`j>0`$. The conditional expectation $`E`$ onto the centre
(Proposition 2.8) is depolarization with $`p=1`$. Unlike the measurement of the sector,
$`X\mapsto\sum_cP_cXP_c`$, which is the identity on $`\mathcal M`$ (Section 3.3), it keeps only $`p_c`$.

*In the library:* the five channels of the table are functions of `channels` — [tutorial, Section 6.4](../tutorial/06_quantum_channels.md#s-6-4).

<a id="r-6-10"></a>

**Theorem 6.10 (fixed points of unital channels).** Let $`\Phi(X)=\sum_lK_lXK_l^*`$ with $`K_l\in\mathcal M`$
be TP and unital. Then $`\Phi(X)=X`$ iff $`X`$ commutes with every $`K_l`$ and $`K_l^*`$. In particular the
fixed points form a von Neumann subalgebra of $`\mathcal M`$.

*Proof.* If $`X`$ commutes with all $`K_l`$, then $`\Phi(X)=X\sum_lK_lK_l^*=X`$. Conversely, $`\Phi(X^*)=\Phi(X)^*`$,
so the self-adjoint parts $`(X+X^*)/2`$ and $`(X-X^*)/2i`$ of a fixed point are fixed points; let
$`X=X^*`$. With $`Y_l=K_lX-XK_l`$,
```math
\sum_lY_lY_l^*=\sum_l\big(K_lX^2K_l^*-K_lXK_l^*X-XK_lXK_l^*+XK_lK_l^*X\big)
=\Phi(X^2)-\Phi(X)X-X\Phi(X)+X\Phi(1)X=\Phi(X^2)-X^2,
```
using $`Y_l^*=XK_l^*-K_l^*X`$, $`\Phi(X)=X`$ and $`\Phi(1)=1`$. Its trace is $`0`$ because $`\Phi`$ is TP. As
$`\mathrm{Tr}\,Y_lY_l^*=\|Y_l\|_2^2\ge0`$, every $`Y_l=0`$; taking adjoints, $`X`$ commutes with $`K_l^*`$ too. These $`X`$ form the commutant of a self-adjoint set,
intersected with $`\mathcal M`$: a $`*`$-algebra containing $`1`$, hence a von Neumann algebra
(Theorem 2.2). $`\square`$

So dephasing ($`p>0`$) keeps exactly the diagonal of every block, and depolarization ($`p>0`$) and $`E`$
keep exactly the centre. Without unitality the theorem fails: amplitude damping fixes $`E_{00}`$,
which does not commute with the jump $`E_{01}`$. For $`0<\gamma\le1`$ and $`k\ge2`$ its fixed points in a block
and the commutant of its Kraus operators both have dimension $`1`$, but they are $`\mathbb CE_{00}`$ and
$`\mathbb C1`$: by the explicit form above, a fixed point has $`\rho_{ij}=0`$ for $`i\ne j`$ and
$`\rho_{ii}=0`$ for $`i>0`$ (these entries are multiplied by $`\sqrt{1-\gamma}\ne1`$ or $`1-\gamma\ne1`$), and every
multiple of $`E_{00}`$ is fixed; and anything commuting with all $`E_{0i}`$, $`E_{i0}`$ ($`i>0`$) commutes
with their products $`E_{ij}=E_{i0}E_{0j}`$ and $`E_{00}=E_{01}E_{10}`$, hence with $`M_k`$, so it is a multiple
of $`1`$ (Chapter 2, Exercise 1).

*Computation:* fixed points and commutants for the five channels — [tutorial, Section 6.4](../tutorial/06_quantum_channels.md#s-6-4).

<a id="r-6-11"></a>

**Theorem 6.11 (unital channels increase entropy).** Let $`\Phi:\mathcal M\to\mathcal N`$ be positive,
TP and unital, and let $`\rho`$ and $`\Phi(\rho)`$ be states with eigenvalues $`\lambda_i`$ and $`\mu_j`$.
Then $`\mu=\Delta\lambda`$ for a doubly stochastic matrix $`\Delta`$, so $`\sum_jf(\mu_j)\ge\sum_if(\lambda_i)`$ for
every concave $`f`$ on $`[0,1]`$. In particular $`S(\Phi(\rho))\ge S(\rho)`$ and
$`\mathrm{Tr}\,\Phi(\rho)^2\le\mathrm{Tr}\,\rho^2`$.

*Proof.* Diagonalise each block: $`\rho=\sum_i\lambda_i|e_i\rangle\langle e_i|`$ and
$`\Phi(\rho)=\sum_j\mu_j|f_j\rangle\langle f_j|`$ with orthonormal bases of $`\bigoplus_c\mathbb C^{k_c}`$
and $`\bigoplus_d\mathbb C^{k'_d}`$ whose vectors lie in one summand each, so that
$`|e_i\rangle\langle e_i|\in\mathcal M`$. Then $`\mu_j=\langle f_j,\Phi(\rho)f_j\rangle=\sum_i\Delta_{ji}\lambda_i`$
with $`\Delta_{ji}=\langle f_j,\Phi(|e_i\rangle\langle e_i|)f_j\rangle\ge0`$,
$`\sum_j\Delta_{ji}=\mathrm{Tr}\,\Phi(|e_i\rangle\langle e_i|)=1`$ and
$`\sum_i\Delta_{ji}=\langle f_j,\Phi(1)f_j\rangle=1`$ (so $`\Delta`$ is square: sum all entries). By the
finite Jensen inequality $`f(\sum_i\Delta_{ji}\lambda_i)\ge\sum_i\Delta_{ji}f(\lambda_i)`$ (induction on
the number of terms from the definition of concavity, as in the solution of Exercise 2 in Chapter 3),
$`\sum_jf(\mu_j)\ge\sum_j\sum_i\Delta_{ji}f(\lambda_i)=\sum_if(\lambda_i)`$. Take $`f(t)=-t\log t`$
($`f(0)=0`$) and $`f(t)=-t^2`$, which are concave on $`[0,1]`$ by Lemma 3.4: $`t\log t`$ and $`t^2`$ are
continuous on $`[0,1]`$ with second derivatives $`1/t>0`$ and $`2>0`$ on $`(0,1)`$. $`\square`$

*Computation:* entropy and purity under the unital channels, and the entropy drop under amplitude damping — [tutorial, Section 6.4](../tutorial/06_quantum_channels.md#s-6-4).

The classical part of a map between sectors is its **transition matrix**
$`T_{dc}(\rho)=\mathrm{Tr}\,\Phi_{dc}(\rho_c)/\mathrm{Tr}\,\rho_c`$ ($`\rho_c\ne0`$); the library's
`transition_matrix()` is $`T=T(1)`$, for the maximally mixed state of every input sector.

<a id="r-6-12"></a>

**Proposition 6.12.** Let $`\Phi:\mathcal M\to\mathcal N`$ be positive and TP, $`\rho`$ a state with
$`p_c=\mathrm{Tr}\,\rho_c>0`$, and $`q_d=\mathrm{Tr}\,\Phi(\rho)_d`$. Then $`T(\rho)`$ is column
stochastic and $`q=T(\rho)p`$. With the conditional expectations $`E_{\mathcal M}`$, $`E_{\mathcal N}`$ onto
the centres (Proposition 2.8),
```math
E_{\mathcal N}\Phi E_{\mathcal M}(\rho)=\bigoplus_d\frac{(Tp)_d}{k'_d}\,1_{k'_d},\qquad T=T(1).
```

*Proof.* $`T_{dc}\ge0`$ as $`\Phi_{dc}`$ is positive, and
$`\sum_dT_{dc}=\mathrm{Tr}\,\Phi(\iota_c\rho_c)/\mathrm{Tr}\,\rho_c=1`$. Next
$`q_d=\sum_c\mathrm{Tr}\,\Phi_{dc}(\rho_c)=\sum_cT_{dc}(\rho)p_c`$. Finally
$`E_{\mathcal M}\rho=\bigoplus_cp_c1_{k_c}/k_c`$, so the trace of
$`\Phi(E_{\mathcal M}\rho)_d=\sum_cp_c\Phi_{dc}(1_{k_c}/k_c)`$ is $`\sum_cT_{dc}p_c`$, and $`E_{\mathcal N}`$
replaces block $`d`$ by its trace times $`1/k'_d`$. $`\square`$

*In the library:* the transition matrix is `InterSectorChannel.transition_matrix()`; a check of Proposition 6.12 — [tutorial, Section 6.4](../tutorial/06_quantum_channels.md#s-6-4).

**Qubit coherences.** On $`M_2`$ let $`D_p`$ be dephasing and $`A_\gamma`$ amplitude damping. By the table
and the explicit form above, $`D_p`$ multiplies $`\rho_{01}`$ and $`\rho_{10}`$ by $`1-p`$ and fixes the
diagonal, while $`A_\gamma`$ multiplies $`\rho_{01}`$ and $`\rho_{10}`$ by $`\sqrt{1-\gamma}`$ and maps the
diagonal to $`(\rho_{00}+\gamma\rho_{11},(1-\gamma)\rho_{11})`$, independently of the coherences. Hence
$`D_pA_\gamma=A_\gamma D_p`$: both multiply the coherences by $`(1-p)\sqrt{1-\gamma}`$ and act on the diagonal
as $`A_\gamma`$. In particular, for $`T_1,T_\varphi>0`$, $`t\ge0`$, $`1-p=e^{-t/T_\varphi}`$ and $`1-\gamma=e^{-t/T_1}`$,
the composition multiplies the coherences by $`e^{-t/T_2}`$ with
```math
\frac1{T_2}=\frac1{2T_1}+\frac1{T_\varphi},
```
and these families compose additively in $`t`$: $`D_{p(t)}D_{p(s)}=D_{p(t+s)}`$ and
$`A_{\gamma(t)}A_{\gamma(s)}=A_{\gamma(t+s)}`$, since the factors $`e^{-t/T}`$ multiply and the trace is preserved.

*Computation:* the coherence factor and $`1/T_2=1/(2T_1)+1/T_\varphi`$ — [tutorial, Section 6.4](../tutorial/06_quantum_channels.md#s-6-4).

> **Physical meaning** (physical interpretation, not a theorem of this course). For a qubit,
> dephasing models pure dephasing with time $`T_\varphi`$ and amplitude damping the decay of the
> excited level with time $`T_1`$; the computation above then gives the coherence time $`T_2`$ of a qubit
> subject to both. Theorem 6.10 is the algebraic
> form of decoherence-free information: what commutes with all jump operators survives (Chapter 10).
> Proposition 6.12 says that, seen through the centre, a map between sectors is a classical Markov
> matrix on the charges.

<a id="s-6-5"></a>

## 6.5 Random channels

A **Haar-random isometry** $`V:\mathbb C^k\to\mathbb C^N`$ is, as in the library (`random_isometry`),
the factor $`Q`$ of $`G=QR`$, where $`G`$ is an $`N\times k`$ Ginibre matrix (independent entries with
density $`\pi^{-1}e^{-|z|^2}`$), $`Q^*Q=1`$, and $`R`$ is upper triangular with positive diagonal (on the
CPU the library takes a LAPACK factorisation $`G=Q'R'`$ and multiplies the $`i`$-th column of $`Q'`$ by
$`R'_{ii}/|R'_{ii}|`$, which gives this factorisation, as in Theorem 5.3, whenever all $`R'_{ii}\ne0`$, i.e.
whenever $`G`$ has rank $`k`$, which holds almost surely by the proof of Lemma 6.13; for large batches on a
GPU it uses Gram–Schmidt, which gives it directly).
`random_channel(alg, r)` draws such $`V_c:\mathbb C^{k_c}\to\mathbb C^r\otimes\mathbb C^{k_c}`$
independently in every block and uses the Kraus operators of Theorem 6.5.

*In the library:* `random_isometry` and `channels.random_channel` — [tutorial, Section 6.5](../tutorial/06_quantum_channels.md#s-6-5).

<a id="r-6-13"></a>

**Lemma 6.13.** $`V`$ has the law of the first $`k`$ columns of a Haar unitary of $`U(N)`$. In
particular, for every unitary $`U\in M_N`$, $`UV`$ has the same law as $`V`$.

*Proof.* Let $`Z`$ be an $`N\times N`$ Ginibre matrix with first $`k`$ columns $`G`$, and $`E\in M_{N\times k}`$
the first $`k`$ columns of $`1_N`$. By Theorem 5.3, $`Z`$ is invertible almost surely and $`Z=Q_ZR_Z`$ with
$`Q_Z`$ Haar distributed. As $`R_Z`$ is upper triangular, $`G=ZE=(Q_ZE)(E^*R_ZE)`$ is a factorisation of
the required type, and it is the only one: if $`G=QR=Q'R'`$, then $`R^*R=G^*G=R'^*R'`$, so $`R'R^{-1}`$
is unitary, upper triangular with positive diagonal, hence $`1`$ (proof of Lemma 5.2). So $`V=Q_ZE`$,
and $`UV=(UQ_Z)E`$ has the law of $`V`$ because $`UQ_Z`$ has that of $`Q_Z`$. $`\square`$

<a id="r-6-14"></a>

**Lemma 6.14.** Let $`u`$ be a random unit vector in $`\mathbb C^N`$ whose law is invariant under all
unitaries. Then $`\mathbb E[u_a\bar u_b]=\delta_{ab}/N`$ and
```math
\mathbb E\big[u_a\bar u_b\,u_c\bar u_d\big]=\frac{\delta_{ab}\delta_{cd}+\delta_{ad}\delta_{cb}}{N(N+1)} .
```

*Proof.* Invariance under $`\mathrm{diag}(e^{i\theta_1},\dots,e^{i\theta_N})`$ multiplies the
moments by $`e^{i(\theta_a-\theta_b)}`$ and $`e^{i(\theta_a-\theta_b+\theta_c-\theta_d)}`$, so they vanish
unless $`a=b`$, resp. $`\{a,c\}=\{b,d\}`$ as multisets. By invariance under permutation matrices
$`\mathbb E|u_a|^2=1/N`$, and only $`\alpha=\mathbb E|u_1|^4`$ and $`\beta=\mathbb E|u_1|^2|u_2|^2`$ remain.
First, $`\sum_c|u_c|^2=1`$ gives $`\alpha+(N-1)\beta=\mathbb E|u_1|^2=1/N`$. Second, $`(u_1+u_2)/\sqrt2`$ is
the first coordinate of $`Uu`$ for a unitary $`U`$, so $`\mathbb E|u_1+u_2|^4=4\alpha`$. With
$`s=|u_1|^2+|u_2|^2`$ and $`z=u_1\bar u_2`$, $`|u_1+u_2|^4=s^2+2s(z+\bar z)+(z+\bar z)^2`$; the
invariances $`u_1\mapsto-u_1`$ and $`u_1\mapsto iu_1`$ give $`\mathbb E[s(z+\bar z)]=0`$ and
$`\mathbb Ez^2=0`$, so $`4\alpha=(2\alpha+2\beta)+2\beta`$, i.e. $`\alpha=2\beta`$. Hence
$`\beta=1/(N(N+1))`$ and $`\alpha=2\beta`$, which is the formula in all cases. (This used $`N\ge2`$; for
$`N=1`$, $`|u_1|=1`$ and both sides equal $`1`$.) $`\square`$

<a id="r-6-15"></a>

**Theorem 6.15 (output purity of a random channel).** Let $`\psi`$ be a pure state in a sector of size
$`k`$ and $`\Phi`$ drawn by `random_channel` with Kraus rank $`r`$. Then
```math
\mathbb E\mathrm{Tr}\,\Phi(\psi)^2=\frac{k+r}{kr+1},\qquad\mathbb E\,\Phi(\psi)=\frac1k1_k .
```

*Proof.* In the sector of $`\psi`$, $`\Phi(\psi)=\mathrm{Tr}_E|u\rangle\langle u|`$ with $`u=V\psi`$, a
unit vector whose law is unitarily invariant by Lemma 6.13. Index $`\mathbb C^r\otimes\mathbb C^k`$ by
pairs $`(e,x)`$. Then $`\Phi(\psi)_{xy}=\sum_eu_{ex}\bar u_{ey}`$, whose mean is $`r\delta_{xy}/(kr)`$, and
$`\mathrm{Tr}\,\Phi(\psi)^2=\sum_{x,y,e,f}u_{ex}\bar u_{ey}u_{fy}\bar u_{fx}`$. By Lemma 6.14 with
$`N=kr`$, the first term forces $`x=y`$ ($`kr^2`$ index tuples) and the second $`e=f`$ ($`k^2r`$ tuples), so
the mean is $`(kr^2+k^2r)/(kr(kr+1))`$. $`\square`$

*Computation:* Theorem 6.15 by Monte Carlo for $`r=1,2,4`$ — [tutorial, Section 6.5](../tutorial/06_quantum_channels.md#s-6-5).

For $`r=1`$ the channel is unitary and the purity is $`1`$; as $`r\to\infty`$ it tends to $`1/k`$, the
purity of the maximally mixed state. The proof used only the unitary invariance of the law of $`u`$,
so it also gives Lubkin's mean purity of the reduced state of a random pure state $`u=g/\|g\|`$,
$`g`$ Gaussian (whose law is unitarily invariant by Theorem A.20(b)), which Chapter 7, Theorem 7.11
proves with Gaussian integrals; the numbers $`\alpha,\beta`$ of Lemma 6.14 are those of Corollary 5.11
for the entries of a Haar unitary, and Corollary A.21 recomputes them for $`u=g/\|g\|`$.

> **Physical meaning** (physical interpretation, not a theorem of this course). A random channel
> models generic, unstructured noise from an environment of
> dimension $`r`$. Already $`r=4`$ brings the purity of a pure qutrit state to $`7/13\approx0.54`$ on
> average (Theorem 6.15), against $`1`$ without noise and the minimum $`1/3`$ (by Cauchy–Schwarz,
> $`\mathrm{Tr}\,\rho^2\ge(\mathrm{Tr}\,\rho)^2/k`$). In a superselected system the noise acts in each sector
> separately, and the relevant dimension is that of the sector, not of the whole space.

<a id="s-6-6"></a>

## 6.6 What to compare with

*In the library:* comparison with other software — [tutorial, Section 6.6](../tutorial/06_quantum_channels.md#s-6-6).

<a id="x-6"></a>

## Exercises

1. (Werner–Holevo) Show that $`\Phi(\rho)=\big(\mathrm{Tr}(\rho)1-\rho^T\big)/(k-1)`$ on $`M_k`$,
   $`k\ge2`$, is a unital channel with Kraus rank $`k(k-1)/2`$, although it is built from the
   transpose.
2. Show that for $`0<p\le1`$ the dephasing channel on $`M_k`$ has Kraus rank $`k`$, one less than the
   number of Kraus operators the library stores.
3. Let $`X`$ be the cyclic shift on $`\mathbb C^k`$ and $`Z=\mathrm{diag}(1,\omega,\dots,\omega^{k-1})`$,
   $`\omega=e^{2\pi i/k}`$. Show that $`\frac1{k^2}\sum_{a,b=0}^{k-1}X^aZ^b\rho Z^{-b}X^{-a}=\mathrm{Tr}(\rho)1/k`$:
   the conditional expectation onto the centre of $`M_k`$ is a mixed-unitary channel.
4. Show that $`(\Psi\circ\Phi)^\sharp=\Phi^\sharp\circ\Psi^\sharp`$ when the same trace is used on the
   middle algebra.

*Computation:* the four exercises numerically — [tutorial, Exercises](../tutorial/06_quantum_channels.md#x-6).

<details>
<summary>Solutions</summary>

1. With $`F`$ the swap, $`J(\theta)=\sum_{ij}E_{ij}\otimes E_{ji}=F`$ and $`J(\rho\mapsto\mathrm{Tr}(\rho)1)=1\otimes1`$,
   so $`J(\Phi)=(1-F)/(k-1)=2P_-/(k-1)`$, where $`P_-=(1-F)/2`$ is the projection onto the
   antisymmetric vectors ($`F=F^*`$, $`F^2=1`$), of rank $`\mathrm{Tr}\,P_-=(k^2-\mathrm{Tr}\,F)/2`$ with
   $`\mathrm{Tr}\,F=\sum_{ij}\langle e_i\otimes e_j,e_j\otimes e_i\rangle=k`$. So $`J\ge0`$ with rank $`k(k-1)/2`$ (Theorem 6.2,
   Corollary 6.3), and $`\mathrm{Tr}\,\Phi(\rho)=\mathrm{Tr}\,\rho`$, $`\Phi(1)=1`$ directly.

2. $`J=(1-p)|\Omega\rangle\langle\Omega|+p\sum_iE_{ii}\otimes E_{ii}`$ vanishes outside the span of the
   $`e_i\otimes e_i`$, and on this $`k`$-dimensional space it is the matrix $`(1-p)\mathbf 1\mathbf 1^T+p1`$
   ($`\mathbf 1`$ the vector of ones), with eigenvalues $`p`$ and $`(1-p)k+p`$, all positive. So
   $`\mathrm{rank}\,J=k`$, while the library stores $`1+k_{\max}`$ operators.

3. $`\mathrm{Tr}((X^aZ^b)^*X^{a'}Z^{b'})=k\,\delta_{aa'}\delta_{bb'}`$: $`Z^{-b}X^{a'-a}Z^{b'}`$ has zero
   diagonal unless $`a=a'`$, and $`\mathrm{Tr}\,Z^{b'-b}=\sum_j\omega^{j(b'-b)}=k\delta_{bb'}`$. So the
   $`B_{ab}=X^aZ^b/\sqrt k`$ are an orthonormal basis of $`M_k`$ for $`\mathrm{Tr}(A^*B)`$, and any
   two orthonormal bases are related by a unitary; by Theorem 6.7, $`\sum_mB_m\rho B_m^*`$ is the
   same for all orthonormal bases. For the basis $`E_{ab}`$ it is
   $`\sum_{ab}E_{ab}\rho E_{ba}=\mathrm{Tr}(\rho)1`$. Divide by $`k`$.

4. $`\mathrm{tr}(\Psi(\Phi(x))z)=\mathrm{tr}(\Phi(x)\Psi^\sharp(z))=\mathrm{tr}(x\,\Phi^\sharp(\Psi^\sharp(z)))`$
   for all $`x`$ and $`z`$; uniqueness of the dual gives the claim.

</details>

**Next:** [Chapter 7](07_tensor_products.md) treats tensor products and entanglement: the partial trace —
a channel between different algebras whose dual is the embedding $`A\mapsto A\otimes1`$ — Page's
formula, and the fusion of charges.
