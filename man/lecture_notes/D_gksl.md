**English** · [Español](es/D_gksl.md) · [Français](fr/D_gksl.md) · [Deutsch](de/D_gksl.md) · [中文](zh/D_gksl.md) · [日本語](ja/D_gksl.md) · [Русский](ru/D_gksl.md)

# Appendix D. Generators of quantum dynamical semigroups

Chapter 8 proved that a generator of Gorini–Kossakowski–Sudarshan–Lindblad (GKSL) form produces a
quantum dynamical semigroup (Theorem 8.6) and used the converse (Theorem 8.7) and its version for
$`\mathcal M=\bigoplus_cM_{k_c}`$. This appendix proves both in finite dimensions, together with the
freedom in the GKSL form and a criterion for the uniqueness of the stationary state. The tool is the
Choi matrix of Chapter 6: a generator is in general not completely positive, but its Choi matrix is positive on
the orthogonal complement of the maximally entangled vector. The same kind of argument — positivity
of a family of operators, then a decomposition — proves in the last section the theorem on operator
algebra error correction stated in Chapter 10 (Theorem 10.7), with an explicit recovery channel.

**You will learn** why the generator of a completely positive semigroup is conditionally completely
positive, how to read off $`H`$ and the jump operators from the Choi matrix of a generator, in what sense
they are unique, why every generator on an algebra with sectors has the form used in Section 8.3, when
the stationary state is unique, and which algebras of observables a recovery channel can protect.

Throughout, $`M_n=M_n(\mathbb C)`$, $`E_{ij}=|i\rangle\langle j|`$ are its matrix units, a
**superoperator** is a linear map on $`M_n`$ (or on $`\mathcal M`$), and
```math
\Omega=\sum_{i=1}^ne_i\otimes e_i\in\mathbb C^n\otimes\mathbb C^n,\qquad
J(\Psi)=\sum_{i,j=1}^nE_{ij}\otimes\Psi(E_{ij})=[\Psi(E_{ij})]_{i,j=1}^n
```
are the maximally entangled vector and the **Choi matrix** of a superoperator $`\Psi`$ (Theorem 6.2;
here $`\Psi`$ need not be positive). For vectors $`v,w`$, $`|v\rangle\langle w|`$ is the operator
$`x\mapsto\langle w,x\rangle v`$. A superoperator $`\Psi`$ is **Hermiticity preserving** if
$`\Psi(A^*)=\Psi(A)^*`$ for all $`A`$. Quantum dynamical semigroups and their generators are those of
Section 8.2 (Lemma 8.5): $`T_t=e^{t\mathcal L}`$ with channels $`T_t`$.

*Computation:* the numerical checks of this appendix — [tutorial, Appendix D](../tutorial/D_gksl.md).

<a id="s-D-1"></a>

## D.1 Choi matrices of the building blocks

<a id="r-D-1"></a>

**Lemma D.1.** Let $`\Psi`$ be a superoperator on $`M_n`$.

(a) $`\Psi\mapsto J(\Psi)`$ is a linear bijection from the superoperators on $`M_n`$ onto $`M_n(M_n)=M_{n^2}`$.
(b) Every $`v\in\mathbb C^n\otimes\mathbb C^n`$ is $`v=(1\otimes K)\Omega`$ for exactly one $`K\in M_n`$, and
$`\langle(1\otimes K)\Omega,(1\otimes L)\Omega\rangle=\mathrm{Tr}(K^*L)`$. In particular
$`\langle\Omega,(1\otimes K)\Omega\rangle=\mathrm{Tr}\,K`$.
(c) For $`K,L\in M_n`$ the superoperator $`A\mapsto KAL^*`$ has the Choi matrix
$`|(1\otimes K)\Omega\rangle\langle(1\otimes L)\Omega|`$. In particular $`J(\mathrm{id})=|\Omega\rangle\langle\Omega|`$.
(d) $`\Psi`$ is Hermiticity preserving iff $`J(\Psi)`$ is self-adjoint. Every CP map is Hermiticity preserving.

*Proof.* (a) The $`(i,j)`$ block of $`J(\Psi)`$ is $`\Psi(E_{ij})`$ and the $`E_{ij}`$ form a basis of $`M_n`$, so
$`J`$ is linear and injective, and every block matrix $`[X_{ij}]`$ is the Choi matrix of the superoperator
with $`E_{ij}\mapsto X_{ij}`$. (b) $`(1\otimes K)\Omega=\sum_ie_i\otimes Ke_i`$ has the $`i`$-th block $`Ke_i`$, the
$`i`$-th column of $`K`$, which determines $`K`$, and every $`v`$ arises in this way. Then
$`\langle\sum_ie_i\otimes Ke_i,\sum_je_j\otimes Le_j\rangle=\sum_i\langle Ke_i,Le_i\rangle=\mathrm{Tr}(K^*L)`$;
for the last claim replace $`(K,L)`$ by $`(1,K)`$. (c) The $`(i,j)`$ block of the right-hand side is
$`Ke_i(Le_j)^*=KE_{ij}L^*`$. (d) The $`(i,j)`$ block of $`J(\Psi)^*`$ is $`\Psi(E_{ji})^*`$. If $`\Psi`$ is
Hermiticity preserving, $`\Psi(E_{ji})^*=\Psi(E_{ji}^*)=\Psi(E_{ij})`$. Conversely, if
$`\Psi(E_{ji})^*=\Psi(E_{ij})`$ for all $`i,j`$ and $`A=\sum a_{ij}E_{ij}`$, then
$`\Psi(A)^*=\sum\bar a_{ij}\Psi(E_{ij})^*=\sum\bar a_{ij}\Psi(E_{ji})=\Psi(A^*)`$. A CP map has a Kraus form
(Theorem 6.2), and $`(\sum_lK_lAK_l^*)^*=\sum_lK_lA^*K_l^*`$. $`\square`$

<a id="r-D-2"></a>

**Lemma D.2 (decomposition).** Let $`q_1,\dots,q_C`$ be orthonormal vectors in
$`\mathbb C^n\otimes\mathbb C^n`$, written $`q_c=(1\otimes Q_c)\Omega`$, let $`W=\mathrm{span}\{q_c\}`$, and let
$`\Psi`$ be a Hermiticity preserving superoperator with $`\langle v,J(\Psi)v\rangle\ge0`$ for all $`v\perp W`$.
Put $`J=J(\Psi)`$ and $`\xi_c=Jq_c-\frac12\sum_d\langle q_d,Jq_c\rangle q_d`$, and let $`K_c\in M_n`$ be defined by
$`(1\otimes K_c)\Omega=\xi_c`$. Then there are $`L_1,\dots,L_r\in M_n`$, $`r=\mathrm{rank}(P_W^\perp JP_W^\perp)`$
with $`P_W^\perp`$ the projection onto $`W^\perp`$, such that
```math
\Psi(A)=\sum_{j=1}^rL_jAL_j^*+\sum_{c=1}^C\big(K_cAQ_c^*+Q_cAK_c^*\big)\qquad(A\in M_n),
```
$`\mathrm{Tr}(Q_c^*L_j)=0`$ for all $`c,j`$, and $`\mathrm{Tr}(L_j^*L_k)=0`$ for $`j\neq k`$.

*Proof.* $`J`$ is self-adjoint by Lemma D.1(d). Let $`P=\sum_c|q_c\rangle\langle q_c|`$ be the projection onto $`W`$
and $`P^\perp=P_W^\perp=1-P`$. Expanding $`J=(P+P^\perp)J(P+P^\perp)`$ and using $`PJP^\perp=PJ-PJP`$ and
$`P^\perp JP=JP-PJP`$,
```math
J=P^\perp JP^\perp+\big(PJ+JP-PJP\big).
```
The first term is self-adjoint, and $`\langle v,P^\perp JP^\perp v\rangle=\langle P^\perp v,JP^\perp v\rangle\ge0`$
because $`P^\perp v\perp W`$; so $`P^\perp JP^\perp\ge0`$. By Theorem 1.3,
$`P^\perp JP^\perp=\sum_{j=1}^r|l_j\rangle\langle l_j|`$, where the $`l_j`$ are the eigenvectors of the nonzero
eigenvalues $`\mu_j`$ scaled by $`\sqrt{\mu_j}`$: they are mutually orthogonal and lie in the range of
$`P^\perp JP^\perp`$, which is orthogonal to $`W`$. Write $`l_j=(1\otimes L_j)\Omega`$; by Lemma D.1(b),
$`\mathrm{Tr}(Q_c^*L_j)=\langle q_c,l_j\rangle=0`$ and $`\mathrm{Tr}(L_j^*L_k)=\langle l_j,l_k\rangle=0`$ for $`j\ne k`$,
and by Lemma D.1(c), $`P^\perp JP^\perp`$ is the Choi matrix of $`A\mapsto\sum_jL_jAL_j^*`$. For the second term,
$`JP=\sum_c|Jq_c\rangle\langle q_c|`$ and $`PJP=\sum_{c,d}\langle q_d,Jq_c\rangle|q_d\rangle\langle q_c|`$, so
```math
\sum_c|\xi_c\rangle\langle q_c|=JP-\tfrac12PJP,\qquad\sum_c|q_c\rangle\langle\xi_c|=PJ-\tfrac12PJP
```
(the second is the adjoint of the first; $`PJP`$ is self-adjoint), and their sum is $`PJ+JP-PJP`$. By
Lemma D.1(c), $`|\xi_c\rangle\langle q_c|`$ and $`|q_c\rangle\langle\xi_c|`$ are the Choi matrices of
$`A\mapsto K_cAQ_c^*`$ and $`A\mapsto Q_cAK_c^*`$. So both sides of the claimed identity have the Choi matrix $`J`$,
and Lemma D.1(a) gives the identity. $`\square`$

<a id="s-D-2"></a>

## D.2 Conditional complete positivity

A superoperator $`\Psi`$ on $`M_n`$ is **conditionally completely positive** (CCP) if
$`\langle v,J(\Psi)v\rangle\ge0`$ for every $`v\perp\Omega`$. By Lemma D.1(b) these $`v`$ are the vectors
$`(1\otimes K)\Omega`$ with $`\mathrm{Tr}\,K=0`$.

<a id="r-D-3"></a>

**Theorem D.3.** For a superoperator $`\mathcal L`$ on $`M_n`$ the following are equivalent:

1. $`e^{t\mathcal L}`$ is CP for every $`t\ge0`$;
2. $`\mathcal L`$ is Hermiticity preserving and conditionally completely positive;
3. there are a CP map $`\Phi`$ on $`M_n`$ and $`G\in M_n`$ with $`\mathcal L(A)=\Phi(A)+GA+AG^*`$ for all $`A`$.

In 3 one can take $`\Phi(A)=\sum_{j=1}^rL_jAL_j^*`$ with $`r\le n^2-1`$, $`\mathrm{Tr}\,L_j=0`$ and
$`\mathrm{Tr}(L_j^*L_k)=0`$ for $`j\neq k`$, and $`G`$ with $`\mathrm{Tr}\,G\in\mathbb R`$.

*Proof.* (1)$`\Rightarrow`$(2). By Lemma 8.1, $`(e^{t\mathcal L}-\mathrm{id})/t\to\mathcal L`$ as $`t\to0`$, and $`J`$ is
linear on a finite-dimensional space, hence continuous, so
$`J(e^{t\mathcal L})=J(\mathrm{id})+t\,J\big((e^{t\mathcal L}-\mathrm{id})/t\big)`$ with the last Choi matrix tending to
$`J(\mathcal L)`$. For $`v\perp\Omega`$, $`\langle v,J(\mathrm{id})v\rangle=|\langle\Omega,v\rangle|^2=0`$ (Lemma D.1(c)), so
for $`t>0`$, as $`J(e^{t\mathcal L})\ge0`$ (Theorem 6.2),
```math
0\le\frac1t\langle v,J(e^{t\mathcal L})v\rangle=\big\langle v,J\big((e^{t\mathcal L}-\mathrm{id})/t\big)v\big\rangle
\;\xrightarrow[t\to0]{}\;\langle v,J(\mathcal L)v\rangle .
```
Each $`e^{t\mathcal L}`$ is CP, hence Hermiticity preserving (Lemma D.1(d)), so
$`\mathcal L(A^*)=\lim_t(e^{t\mathcal L}(A)^*-A^*)/t=\mathcal L(A)^*`$.

(2)$`\Rightarrow`$(3). Apply Lemma D.2 with $`C=1`$ and $`q_1=\Omega/\sqrt n`$, i.e. $`Q_1=1/\sqrt n`$:
$`\mathcal L(A)=\sum_jL_jAL_j^*+GA+AG^*`$ with $`G=K_1/\sqrt n`$, and $`\Phi=\sum_jL_j\cdot L_j^*`$ is CP
(Theorem 6.2). Further $`\mathrm{Tr}\,L_j=\sqrt n\mathrm{Tr}(Q_1^*L_j)=0`$, $`r\le\mathrm{rank}\,P_W^\perp=n^2-1`$,
and by Lemma D.1(b)
$`\mathrm{Tr}\,G=\langle\Omega,\xi_1\rangle/\sqrt n=\langle q_1,\xi_1\rangle=\frac12\langle q_1,J(\mathcal L)q_1\rangle\in\mathbb R`$,
as $`J(\mathcal L)`$ is self-adjoint.

(3)$`\Rightarrow`$(1). Steps 1–4 of the proof of Theorem 8.6 use only that the generator is
$`\mathcal L_0+\Phi`$ with $`\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)`$ for some matrix
$`H_{\mathrm{eff}}`$ and $`\Phi`$ CP. With $`H_{\mathrm{eff}}=iG`$, $`\mathcal L_0(\rho)=G\rho+\rho G^*`$, and these
steps show that $`e^{t\mathcal L}=\sum_{m\ge0}S_m(t)`$ is CP. $`\square`$

*Computation:* Theorem D.3 on a GKSL generator on $`M_3`$, and the example $`\theta-\mathrm{id}`$ below — [tutorial, Section D.2](../tutorial/D_gksl.md#s-D-2). *In the library:* the matrix of a GKSL generator, `dynamics.lindblad_superoperator` — [tutorial, Section 8.2](../tutorial/08_dynamics.md#s-8-2).

Positivity alone would not do. Let $`n\ge2`$, let $`\theta`$ be the transpose on $`M_n`$ and $`\mathcal L=\theta-\mathrm{id}`$. As $`\theta^2=\mathrm{id}`$,
```math
e^{t\mathcal L}=e^{-t}\big(\cosh t\,\mathrm{id}+\sinh t\,\theta\big)=\tfrac{1+e^{-2t}}2\,\mathrm{id}+\tfrac{1-e^{-2t}}2\,\theta ,
```
a convex combination of positive trace preserving maps (Proposition 6.1). But $`J(\theta)=[E_{ji}]=\sum_{i,j}E_{ij}\otimes E_{ji}`$ is the swap $`e_k\otimes e_l\mapsto e_l\otimes e_k`$, so it maps
$`v=e_1\otimes e_2-e_2\otimes e_1`$, which is orthogonal to $`\Omega`$, to $`-v`$ (as in the proof of Proposition 6.1); with
$`\|v\|^2=2`$ and $`\langle v,J(\mathrm{id})v\rangle=|\langle\Omega,v\rangle|^2=0`$,
$`\langle v,J(\mathcal L)v\rangle=-2`$: $`\mathcal L`$ is not CCP, and indeed
$`\langle v,J(e^{t\mathcal L})v\rangle=-(1-e^{-2t})<0`$ for every $`t>0`$. A semigroup of positive trace preserving maps
need not consist of channels.

<a id="s-D-3"></a>

## D.3 The GKSL theorem

<a id="r-D-4"></a>

**Theorem D.4 (Gorini–Kossakowski–Sudarshan; Lindblad; Theorem 8.7).** A superoperator $`\mathcal L`$ on
$`M_n`$ is the generator of a quantum dynamical semigroup on $`M_n`$ iff there are $`H=H^*`$ and
$`L_1,\dots,L_r`$ in $`M_n`$ with
```math
\mathcal L(\rho)=-i[H,\rho]+\sum_{j=1}^r\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
```
One can choose $`r\le n^2-1`$, $`\mathrm{Tr}\,H=0`$, $`\mathrm{Tr}\,L_j=0`$, and $`\mathrm{Tr}(L_j^*L_k)=0`$ for $`j\neq k`$.

*Proof.* "If" is Theorem 8.6. Conversely, let $`(T_t)`$ be a quantum dynamical semigroup with generator
$`\mathcal L`$, so $`T_t=e^{t\mathcal L}`$ (Lemma 8.5). As every $`T_t`$ is CP, Theorem D.3 gives
$`\mathcal L(\rho)=\sum_jL_j\rho L_j^*+G\rho+\rho G^*`$ with $`r\le n^2-1`$ traceless, mutually orthogonal $`L_j`$ and
$`\mathrm{Tr}\,G\in\mathbb R`$. Put $`H=\frac i2(G-G^*)`$ and $`R=\frac12(G+G^*)`$, both self-adjoint; then
$`G=-iH+R`$, $`G^*=iH+R`$, and $`G\rho+\rho G^*=-i[H,\rho]+\{R,\rho\}`$. Moreover
$`\mathrm{Tr}\,H=\frac i2(\mathrm{Tr}\,G-\overline{\mathrm{Tr}\,G})=0`$. Since every $`T_t`$ is TP,
$`\mathrm{Tr}\,\mathcal L(\rho)=\frac d{dt}\mathrm{Tr}\,T_t(\rho)\big|_{t=0}=0`$ (Lemma 8.1), and by cyclicity
```math
0=\mathrm{Tr}\,\mathcal L(\rho)=\mathrm{Tr}\Big(\rho\,\Big(\sum_jL_j^*L_j+2R\Big)\Big)\qquad\text{for all }\rho .
```
With $`\rho=X=\sum_jL_j^*L_j+2R`$, which is self-adjoint, $`\mathrm{Tr}(X^*X)=0`$, so $`X=0`$ and
$`R=-\frac12\sum_jL_j^*L_j`$, which is the GKSL form. $`\square`$

*Computation:* the proof of Theorem D.4 run as an algorithm, which recovers $`H`$ and the jumps from $`\mathcal L`$ — [tutorial, Section D.3](../tutorial/D_gksl.md#s-D-3).

<a id="r-D-5"></a>

**Theorem D.5 (freedom in the GKSL form).** Let $`(H;L_1,\dots,L_r)`$ and $`(H';L'_1,\dots,L'_s)`$ be GKSL
data ($`H,H'`$ self-adjoint) with generators $`\mathcal L`$, $`\mathcal L'`$ as in Theorem D.4.

(a) $`\mathcal L'=\mathcal L`$ in each of the following cases:
(i) $`L'_j=L_j`$ and $`H'=H+c1`$, $`c\in\mathbb R`$;
(ii) $`H'=H`$ and, after padding the shorter family with zeros, $`L'_j=\sum_ku_{jk}L_k`$ with a unitary $`u`$;
(iii) $`L'_j=L_j+a_j1`$ with $`a_j\in\mathbb C`$, and $`H'=H-\frac i2\sum_j(\bar a_jL_j-a_jL_j^*)`$, which is self-adjoint.

(b) If $`\mathrm{Tr}\,H=\mathrm{Tr}\,H'=0`$, all $`L_j,L'_j`$ are traceless and $`\mathcal L=\mathcal L'`$, then $`H=H'`$ and,
after padding, $`L'_j=\sum_ku_{jk}L_k`$ with a unitary $`u`$.

(c) Hence two GKSL data have the same generator iff one is obtained from the other by finitely many
steps (i)–(iii). In particular the operator $`H`$ of the normalised form of Theorem D.4 and its
dissipative part $`\Phi(\rho)=\sum_jL_j\rho L_j^*`$ are determined by $`\mathcal L`$.

*Proof.* (a) (i) $`[c1,\rho]=0`$. (ii) By Theorem 6.7, $`\Phi'=\sum_jL'_j\cdot L_j'^*`$ equals
$`\Phi=\sum_jL_j\cdot L_j^*`$; the dual does not depend on the Kraus operators (Section 6.2), so
$`\sum_jL_j'^*L'_j=\Phi^*(1)=\sum_jL_j^*L_j`$. (iii) For one jump $`L`$ and $`a\in\mathbb C`$, expanding gives
```math
(L+a)\rho(L+a)^*-\tfrac12\{(L+a)^*(L+a),\rho\}-\Big(L\rho L^*-\tfrac12\{L^*L,\rho\}\Big)=\tfrac12[B,\rho],
\qquad B=\bar aL-aL^* ,
```
since the terms $`\bar aL\rho+a\rho L^*`$ from the first product combine with
$`-\frac12(aL^*\rho+a\rho L^*+\bar aL\rho+\bar a\rho L)`$ and the terms with $`|a|^2`$ cancel. As $`B^*=-B`$, the
matrix $`\frac i2B`$ is self-adjoint, so $`H'=H-\frac i2\sum_jB_j`$ is self-adjoint, and
$`-i[H'-H,\rho]=-\frac12\sum_j[B_j,\rho]`$ cancels the change of the dissipator.

(b) Let $`\Phi=\sum_jL_j\cdot L_j^*`$ and $`G=-iH-\frac12\Phi^*(1)`$, so that $`\mathcal L(\rho)=\Phi(\rho)+G\rho+\rho G^*`$.
By Lemma D.1(c),
$`J(\mathcal L)=\sum_j|l_j\rangle\langle l_j|+|g\rangle\langle\Omega|+|\Omega\rangle\langle g|`$ with $`l_j=(1\otimes L_j)\Omega`$ and
$`g=(1\otimes G)\Omega`$, and $`l_j\perp\Omega`$ because $`\mathrm{Tr}\,L_j=0`$ (Lemma D.1(b)). With $`P^\perp`$ the projection
onto $`\Omega^\perp`$, therefore $`P^\perp J(\mathcal L)P^\perp=\sum_j|l_j\rangle\langle l_j|=J(\Phi)`$, and likewise
$`P^\perp J(\mathcal L')P^\perp=J(\Phi')`$. So $`J(\Phi)=J(\Phi')`$, $`\Phi=\Phi'`$ (Lemma D.1(a)), and Theorem 6.7 gives the
unitary $`u`$. Then $`\Phi^*(1)=\Phi'^*(1)`$, hence $`[H-H',\rho]=0`$ for all $`\rho`$. A matrix $`X`$ that commutes with all
$`E_{ij}`$ is a multiple of $`1`$ (compare the entries of $`XE_{ij}`$ and $`E_{ij}X`$), so $`H-H'=\lambda1`$, and
$`\lambda=0`$ by taking traces.

(c) By (a), the steps (i)–(iii) do not change the generator, and the inverse of each step is a step of the
same kind ((iii) with $`-a_j`$, since $`\bar a_jL'_j-a_jL_j'^*=\bar a_jL_j-a_jL_j^*`$). Any data are brought to
normalised form by (iii) with $`a_j=-\mathrm{Tr}\,L_j/n`$ followed by (i) with $`c=-\mathrm{Tr}\,H''/n`$, $`H''`$ the
Hamiltonian after (iii); two normalised data with the same generator are related by (ii), by (b). $`\square`$

*Computation:* Theorem D.5: the normalised form by (iii) and (i), and the uniqueness in (b) — [tutorial, Section D.3](../tutorial/D_gksl.md#s-D-3).

<a id="r-D-6"></a>

**Theorem D.6 (uniqueness of the stationary state).** Let $`\mathcal L`$ be the GKSL generator of Theorem D.4
with data $`(H;L_1,\dots,L_r)`$. If $`\mathcal L`$ has a faithful stationary state $`\sigma`$ ($`\sigma>0`$) and the only
matrices commuting with $`H`$ and all $`L_j`$, $`L_j^*`$ are the multiples of $`1`$, then $`\ker\mathcal L=\mathbb C\sigma`$:
$`\sigma`$ is the only stationary state. Moreover $`\ker\mathcal L^*=\mathbb C1`$.

*Proof.* The dual for the trace pairing, $`\mathrm{Tr}(\mathcal L(\rho)A)=\mathrm{Tr}(\rho\,\mathcal L^*(A))`$, is
$`\mathcal L^*(A)=i[H,A]+\sum_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})`$: by cyclicity
$`\mathrm{Tr}(-i[H,\rho]A)=\mathrm{Tr}(\rho\,i[H,A])`$, $`\mathrm{Tr}(L\rho L^*A)=\mathrm{Tr}(\rho L^*AL)`$ and
$`\mathrm{Tr}(\{L^*L,\rho\}A)=\mathrm{Tr}(\rho\{L^*L,A\})`$ (this is Exercise 1 of Chapter 8).
A direct expansion, in which the Hamiltonian terms cancel because $`A\mapsto i[H,A]`$ is a derivation, gives
for every $`X`$
```math
\mathcal L^*(X^*X)-X^*\mathcal L^*(X)-\mathcal L^*(X^*)X=\sum_j[L_j,X]^*[L_j,X]:
```
for one jump, of the six anticommutator terms four cancel in pairs and the other two,
$`\frac12X^*L^*LX+\frac12X^*L^*LX`$, add up, so the left side is $`L^*X^*XL-X^*L^*XL-L^*X^*LX+X^*L^*LX`$, which is
$`(XL-LX)^*(XL-LX)`$. Let $`X\in\ker\mathcal L^*`$. Then $`\mathcal L^*(X^*)=\mathcal L^*(X)^*=0`$, as $`\mathcal L^*`$
is Hermiticity preserving, so $`\mathcal L^*(X^*X)=\sum_j[L_j,X]^*[L_j,X]=:Y\ge0`$, and
$`\mathrm{Tr}(\sigma Y)=\mathrm{Tr}(\mathcal L(\sigma)X^*X)=0`$. Since $`\sigma>0`$,
$`\mathrm{Tr}(\sigma Y)=\mathrm{Tr}(\sigma^{1/2}Y\sigma^{1/2})`$ with $`\sigma^{1/2}Y\sigma^{1/2}\ge0`$ vanishes only if
$`\sigma^{1/2}Y\sigma^{1/2}=0`$, i.e. $`Y=0`$; taking traces, $`\sum_j\mathrm{Tr}([L_j,X]^*[L_j,X])=0`$ with nonnegative
terms, so $`[L_j,X]=0`$ for all $`j`$. The same argument for $`X^*\in\ker\mathcal L^*`$
gives $`[L_j,X^*]=0`$, i.e. $`[L_j^*,X]=0`$. Then $`X`$ commutes with $`L_j^*L_j`$, the dissipative part of $`\mathcal L^*(X)`$
vanishes, and $`0=\mathcal L^*(X)=i[H,X]`$. By hypothesis $`X\in\mathbb C1`$; conversely $`\mathcal L^*(1)=0`$. So
$`\ker\mathcal L^*=\mathbb C1`$. In the bases $`(E_{ij})`$ and $`(E_{ji})`$ of $`M_n`$, which are dual for the pairing
$`(\rho,A)\mapsto\mathrm{Tr}(\rho A)`$, the matrix of $`\mathcal L^*`$ is the transpose of that of $`\mathcal L`$; so both have the
same rank, and $`\dim\ker\mathcal L=\dim\ker\mathcal L^*=1`$. A state in $`\mathbb C\sigma`$ with trace $`1`$ is $`\sigma`$. $`\square`$

This is the finite-dimensional case of a criterion of Frigerio [*Commun. Math. Phys.* 63, 269 (1978)],
used in Section 8.5. Both hypotheses are needed. In Exercise 4 of Chapter 8 the state $`\tfrac12\,1`$ is faithful
and stationary, but $`H`$ and $`\sigma_z`$ commute with all diagonal matrices, and the stationary states span the
two-dimensional space of diagonal matrices. Conversely, on $`M_3`$ let $`H=0`$, $`L_1=|e_1\rangle\langle e_3|`$,
$`L_2=|e_2\rangle\langle e_3|`$. The $`*`$-algebra generated by $`L_1,L_2`$ contains $`L_1L_2^*=E_{12}`$, $`L_1^*L_1=E_{33}`$ and hence
all matrix units, so only the multiples of $`1`$ commute with $`L_j,L_j^*`$; but every $`\rho=P\rho P`$ with $`P=E_{11}+E_{22}`$
is stationary, as $`L_j\rho=0=\rho L_j^*`$ and $`\{E_{33},\rho\}=0`$, so $`\dim\ker\mathcal L\ge4`$; by Theorem D.6 no stationary
state is faithful.

*Computation:* Theorem D.6 for the generator of Section D.2, and the second example — [tutorial, Section D.3](../tutorial/D_gksl.md#s-D-3).

<a id="s-D-4"></a>

## D.4 Algebras with sectors

Let $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}\subseteq M_D`$ be the block-diagonal matrices, $`D=\sum_ck_c`$, with
minimal central projections $`P_c`$, and let $`E(X)=\sum_cP_cXP_c`$ be the conditional expectation of $`M_D`$ onto
$`\mathcal M`$. As in the proof of Proposition 8.8, $`E(B)=B`$, $`E(XB)=E(X)B`$ and $`E(BX)=BE(X)`$ for $`B\in\mathcal M`$,
and $`E(X^*)=E(X)^*`$. For a CP map $`\Phi:\mathcal M\to\mathcal M`$, $`\Phi^*`$ is its blunt dual (Proposition 6.8).

<a id="r-D-7"></a>

**Theorem D.7.** Let $`\mathcal L`$ be the generator of a quantum dynamical semigroup $`(T_t)`$ on $`\mathcal M`$. There
are a CP map $`\Phi:\mathcal M\to\mathcal M`$ and $`H=H^*\in\mathcal M`$ with
```math
\mathcal L(\rho)=\Phi(\rho)-\tfrac12\{\Phi^*(1),\rho\}-i[H,\rho]\qquad(\rho\in\mathcal M).
```
By Proposition 6.4, $`\Phi(\rho)_d=\sum_{c,l}K^{dc}_l\rho_cK^{dc*}_l`$ with Kraus operators between the sectors;
conversely every such $`\mathcal L`$ generates a quantum dynamical semigroup on $`\mathcal M`$ (Proposition 8.8(1)).
So the generators of Section 8.3 (where $`\Phi=\sum_jg_jL_j\cdot L_j^*+\sum_mg_m\Phi_m`$) are exactly the
generators of quantum dynamical semigroups on $`\mathcal M`$.

*Proof.* *Step 1: a CP map on $`M_D`$.* By Proposition 6.4,
$`T_t(\rho)_d=\sum_{c,l}K^{dc}_l(t)\rho_cK^{dc}_l(t)^*`$. Let $`N^{dc}_l(t)\in M_D`$ have the block $`K^{dc}_l(t)`$ in
position $`(d,c)`$ and zeros elsewhere. For $`X\in M_D`$, $`N^{dc}_lXN^{dc*}_l`$ is the block
$`K^{dc}_lX_{cc}K^{dc*}_l`$ in position $`(d,d)`$, where $`X_{cc}`$ is the $`c`$-th diagonal block, i.e. the $`c`$-th block
of $`E(X)`$. Hence $`\hat T_t=T_t\circ E=\sum_{d,c,l}N^{dc}_l\cdot N^{dc*}_l`$ is CP on $`M_D`$ (Theorem 6.2), with
$`\hat T_0=E`$. By Lemma 8.5, $`(\hat T_t-E)/t=\big((T_t-\mathrm{id})/t\big)\circ E\to\mathcal L\circ E`$.

*Step 2: positivity off the sector vectors.* Let $`p_c=(1\otimes P_c)\Omega\in\mathbb C^D\otimes\mathbb C^D`$ and
$`W=\mathrm{span}\{p_c\}`$. By Lemma D.1(c), $`J(E)=\sum_c|p_c\rangle\langle p_c|`$, so for $`v\perp W`$ and $`t>0`$,
$`0\le\frac1t\langle v,J(\hat T_t)v\rangle=\langle v,J((\hat T_t-E)/t)v\rangle\to\langle v,J(\mathcal L\circ E)v\rangle`$.
Every $`T_t`$ is Hermiticity preserving (its Kraus form), hence so is $`\mathcal L`$ (as in Theorem D.3), and so are $`E`$
and $`\mathcal L\circ E`$.

*Step 3: decomposition.* By Lemma D.1(b), $`\langle p_c,p_d\rangle=\mathrm{Tr}(P_cP_d)=\delta_{cd}k_c`$, so
$`q_c=p_c/\sqrt{k_c}=(1\otimes Q_c)\Omega`$ with $`Q_c=P_c/\sqrt{k_c}`$ are orthonormal. Lemma D.2 for
$`\Psi=\mathcal L\circ E`$ gives a CP map $`\hat\Phi=\sum_j\hat L_j\cdot\hat L_j^*`$ on $`M_D`$ and $`K_c\in M_D`$ with
```math
\mathcal L(E(A))=\hat\Phi(A)+\sum_c\big(K_cAP_c+P_cAK_c^*\big)/\sqrt{k_c}\qquad(A\in M_D).
```
For $`A\in\mathcal M`$, $`E(A)=A`$ and $`AP_c=P_cA`$, so the sum is $`\tilde GA+A\tilde G^*`$ with
$`\tilde G=\sum_cK_cP_c/\sqrt{k_c}`$.

*Step 4: compression.* For $`A\in\mathcal M`$, $`\mathcal L(A)\in\mathcal M`$, so
$`\mathcal L(A)=E(\mathcal L(A))=E(\hat\Phi(A))+GA+AG^*`$ with $`G=E(\tilde G)\in\mathcal M`$. The map
$`\Phi=E\circ\hat\Phi`$ on $`\mathcal M`$ is
$`\Phi(A)_d=\sum_{c,j}(P_d\hat L_jP_c)A_c(P_d\hat L_jP_c)^*`$: it is CP by Proposition 6.4, with Kraus operators the
blocks of the $`\hat L_j`$.

*Step 5: trace preservation.* As in Theorem D.4, $`\mathrm{Tr}\,\mathcal L(A)=0`$ for $`A\in\mathcal M`$, and
$`\mathrm{Tr}\,\Phi(A)=\mathrm{Tr}(A\Phi^*(1))`$, so $`\mathrm{Tr}(AX)=0`$ for all $`A\in\mathcal M`$, where
$`X=\Phi^*(1)+G+G^*\in\mathcal M`$ is self-adjoint. With $`A=X`$, $`X=0`$ (the blunt trace is faithful, Theorem 3.1).
With $`H=\frac i2(G-G^*)\in\mathcal M`$, $`G=-iH-\frac12\Phi^*(1)`$ and $`G\rho+\rho G^*=-i[H,\rho]-\frac12\{\Phi^*(1),\rho\}`$. $`\square`$

The vectors $`p_c`$ replace $`\Omega=\sum_cp_c`$ because $`\hat T_0=E`$, not $`\mathrm{id}`$: $`J(\mathcal L\circ E)`$ need not be
positive on $`\Omega^\perp`$, and in the example below it is not.

*Computation:* Theorem D.7 on $`M_1\oplus M_2\subseteq M_3`$, following the five steps of the proof; there $`J(\mathcal L\circ E)`$ is not positive on $`\Omega^\perp`$ — [tutorial, Section D.4](../tutorial/D_gksl.md#s-D-4). *In the library:* generators on an algebra with sectors, `dynamics.lindblad_rhs` with `channels.random_inter_sector_channel` — [tutorial, Section 8.3](../tutorial/08_dynamics.md#s-8-3).

<a id="r-D-8"></a>

**Corollary D.8 (Heisenberg picture).** Let $`(S_t)_{t\ge0}`$ be a continuous family of unital CP maps on
$`\mathcal M`$ with $`S_0=\mathrm{id}`$ and $`S_{s+t}=S_sS_t`$. Its generator is $`A\mapsto\Psi(A)+K^*A+AK`$ with
$`\Psi`$ CP on $`\mathcal M`$, $`K\in\mathcal M`$ and $`\Psi(1)+K+K^*=0`$.

*Proof.* Let $`T_t=S_t^*`$ (blunt duals). They are CP (Proposition 6.8) and TP (Corollary 6.9(b)),
$`T_0=\mathrm{id}`$, and $`T_{s+t}=(S_tS_s)^*=S_s^*S_t^*=T_sT_t`$, because
$`\mathrm{Tr}(S_tS_s(x)y)=\mathrm{Tr}(xS_s^*S_t^*(y))`$ and $`S_tS_s=S_{s+t}`$. As $`\Psi\mapsto\Psi^*`$ is linear on a
finite-dimensional space, $`t\mapsto T_t`$ is continuous and $`(T_t-\mathrm{id})/t\to\mathcal L_S^*`$, $`\mathcal L_S`$ the
generator of $`(S_t)`$. Theorem D.7 gives $`\mathcal L_S^*(\rho)=\Phi(\rho)+G\rho+\rho G^*`$ with $`G+G^*=-\Phi^*(1)`$, and by
cyclicity, $`\mathrm{Tr}((G\rho+\rho G^*)A)=\mathrm{Tr}(\rho(G^*A+AG))`$, so by Corollary 6.9(c)
$`\mathcal L_S(A)=\Phi^*(A)+G^*A+AG`$. Take $`\Psi=\Phi^*`$ (CP by Proposition 6.8) and $`K=G`$. $`\square`$

This is the finite-dimensional case of the theorem of Christensen and Evans quoted after Theorem 8.7.

*Computation:* the generator of Corollary D.8 as the dual of $`\mathcal L`$ — [tutorial, Section D.4](../tutorial/D_gksl.md#s-D-4).

<a id="s-D-5"></a>

## D.5 Operator algebra error correction

We use the setting of Section 10.3. $`\Phi(\rho)=\sum_{i=1}^rE_i\rho E_i^*`$ is a channel on $`M_D`$, $`\mathcal C=P\mathbb C^D`$
a code, $`B(\mathcal C)=\{X\in M_D:X=PXP\}`$ the operators on $`\mathcal C`$, $`\mathfrak A\subseteq B(\mathcal C)`$ a $`*`$-algebra
with unit $`P`$ (a von Neumann algebra on $`\mathcal C`$ by Theorem 2.2), and $`\mathfrak A'`$ its commutant in
$`B(\mathcal C)`$. $`\mathfrak A`$ is **correctable** for $`\Phi`$ if there is a channel $`\mathcal R`$ on $`M_D`$ with
$`P(\mathcal R\circ\Phi)^*(A)P=A`$ for all $`A\in\mathfrak A`$. If $`\mathcal R`$ has the Kraus operators $`R_m`$, then
$`\mathcal R\circ\Phi`$ has the Kraus operators $`R_mE_i`$, and $`\sum_mR_m^*R_m=1`$.

<a id="r-D-9"></a>

**Theorem D.9 (Bény–Kempf–Kribs; Theorem 10.7).** $`\mathfrak A`$ is correctable for $`\Phi`$ iff
$`PE_i^*E_jP\in\mathfrak A'`$ for all $`i,j`$. In that case the channel $`\mathcal R`$ constructed in the proof corrects $`\mathfrak A`$.

*Proof.* *Necessity.* Let $`\mathcal R`$ correct $`\mathfrak A`$ and put $`W_{mi}=R_mE_iP`$, so that
$`\sum_{m,i}W_{mi}^*XW_{mi}=P(\mathcal R\circ\Phi)^*(X)P`$ for all $`X\in M_D`$. This equals $`X`$ for $`X\in\mathfrak A`$, and $`P`$ for
$`X=1`$, because $`(\mathcal R\circ\Phi)^*`$ is unital (Corollary 6.9(a)). Let $`p\in\mathfrak A`$ be a projection; then $`p\le P`$,
and $`\sum W_{mi}^*(1-p)W_{mi}=P-p`$. For $`x=px`$,
$`\sum_{m,i}\|(1-p)W_{mi}x\|^2=\langle x,(P-p)x\rangle=0`$, and for $`x=(P-p)x`$,
$`\sum_{m,i}\|pW_{mi}x\|^2=\langle x,px\rangle=0`$. Hence $`(1-p)W_{mi}p=0`$ and $`pW_{mi}(P-p)=0`$, and since
$`W_{mi}=W_{mi}P`$,
```math
W_{mi}p=pW_{mi}p=pW_{mi}P=pW_{mi}.
```
Every self-adjoint $`A\in\mathfrak A`$ is a linear combination of projections in $`\mathfrak A`$: by Theorem 1.3 on $`\mathcal C`$,
$`A=\sum_k\lambda_kp_k`$ with distinct $`\lambda_k`$ and orthogonal projections $`p_k`$, $`\sum_kp_k=P`$, and
$`p_k=\prod_{l\neq k}(A-\lambda_lP)/(\lambda_k-\lambda_l)\in\mathfrak A`$ (apply both sides to the eigenvectors in $`\mathcal C`$ and to $`\mathcal C^\perp`$; the empty product is $`P`$). Every element
of $`\mathfrak A`$ is $`A_1+iA_2`$ with $`A_1,A_2\in\mathfrak A`$ self-adjoint. So $`W_{mi}A=AW_{mi}`$ for all $`A\in\mathfrak A`$, and,
taking adjoints, $`W_{mi}^*A=AW_{mi}^*`$. As $`\sum_mW_{mi}^*W_{mj}=PE_i^*(\sum_mR_m^*R_m)E_jP=PE_i^*E_jP`$,
```math
A\,PE_i^*E_jP=\sum_mW_{mi}^*AW_{mj}=PE_i^*E_jP\,A\qquad(A\in\mathfrak A),
```
i.e. $`PE_i^*E_jP\in\mathfrak A'`$.

*Sufficiency.* By Theorem 2.6 for $`\mathfrak A`$ on $`\mathcal C`$ there is an orthonormal basis $`|c;a,b\rangle`$
($`c=1,\dots,C`$, $`a=1,\dots,k_c`$, $`b=1,\dots,m_c`$) of $`\mathcal C`$ — the images of $`e_a\otimes f_b`$ in the $`c`$-th summand under
the unitary $`V`$ — in which $`\mathfrak A`$ consists of the operators $`\sum_{c,a,a',b}(A_c)_{aa'}|c;a,b\rangle\langle c;a',b|`$,
$`A_c\in M_{k_c}`$, and $`\mathfrak A'`$ of the operators $`\sum_{c,a,b,b'}(Y_c)_{bb'}|c;a,b\rangle\langle c;a,b'|`$, $`Y_c\in M_{m_c}`$.
The hypothesis therefore gives matrices $`\Lambda^c_{ij}\in M_{m_c}`$ with
```math
\langle c;a,b|E_i^*E_j|d;a',b'\rangle=\delta_{cd}\,\delta_{aa'}\,(\Lambda^c_{ij})_{bb'} .
```
*Step 1.* For each $`c`$ and each pair $`\alpha=(i,b)`$ let $`F^c_\alpha:\mathbb C^{k_c}\to\mathbb C^D`$,
$`F^c_\alpha e_a=E_i|c;a,b\rangle`$, and $`\Lambda^c_{\alpha\beta}=(\Lambda^c_{ij})_{bb'}`$ for $`\beta=(j,b')`$. The display says
$`F^{c*}_\alpha F^d_\beta=\delta_{cd}\Lambda^c_{\alpha\beta}1_{k_c}`$: a Knill–Laflamme condition for the "errors" $`F^c_\alpha`$ on the
code $`\mathbb C^{k_c}`$, with in addition orthogonality between the sectors.

*Step 2.* As in the proof of Theorem 10.5, $`\Lambda^c`$ is self-adjoint, since
$`\Lambda^c_{\beta\alpha}1=F^{c*}_\beta F^c_\alpha=(F^{c*}_\alpha F^c_\beta)^*=\bar\Lambda^c_{\alpha\beta}1`$, and positive, since
$`\sum\bar z_\alpha z_\beta\Lambda^c_{\alpha\beta}1=(\sum z_\alpha F^c_\alpha)^*(\sum z_\beta F^c_\beta)\ge0`$. Write
$`\Lambda^c=U^c\mathrm{diag}(d^c_\kappa)U^{c*}`$ with $`U^c`$ unitary and $`d^c_\kappa\ge0`$, and
$`G^c_\kappa=\sum_\alpha u^c_{\alpha\kappa}F^c_\alpha`$. Then $`G^{c*}_\kappa G^d_\lambda=\delta_{cd}\delta_{\kappa\lambda}d^c_\kappa1`$, so
$`G^c_\kappa=0`$ if $`d^c_\kappa=0`$, and otherwise $`V^c_\kappa=G^c_\kappa/\sqrt{d^c_\kappa}`$ is an isometry with
$`V^{c*}_\kappa V^d_\lambda=\delta_{cd}\delta_{\kappa\lambda}1`$. Inverting, $`F^c_\alpha=\sum_\kappa\bar u^c_{\alpha\kappa}G^c_\kappa =\sum_{\kappa:d^c_\kappa>0}\bar u^c_{\alpha\kappa}\sqrt{d^c_\kappa}\,V^c_\kappa`$, hence
```math
V^{c'*}_\kappa F^d_\beta=\delta_{c'd}\,\bar u^d_{\beta\kappa}\sqrt{d^d_\kappa}\,1,\qquad
F^{c*}_\alpha V^{c'}_\kappa=\delta_{cc'}\,u^c_{\alpha\kappa}\sqrt{d^c_\kappa}\,1 .
```
*Step 3: the recovery.* The $`V^c_\kappa V^{c*}_\kappa`$ ($`d^c_\kappa>0`$) are mutually orthogonal projections; let $`Q`$ be their
sum. Let $`\iota_c:\mathbb C^{k_c}\to\mathbb C^D`$, $`\iota_ce_a=|c;a,1\rangle`$, an isometry, and
```math
\mathcal R(\sigma)=\sum_{c}\sum_{\kappa:d^c_\kappa>0}\iota_cV^{c*}_\kappa\,\sigma\,V^c_\kappa\iota_c^*+(1-Q)\sigma(1-Q).
```
It is a channel: $`\sum_{c,\kappa}V^c_\kappa\iota_c^*\iota_cV^{c*}_\kappa+(1-Q)^2=Q+1-Q=1`$. It measures the sector $`c`$ and
the syndrome $`\kappa`$, undoes the error on $`\mathbb C^{k_c}`$ and writes the result into the copy $`b=1`$.

*Step 4: correctness.* $`\mathcal R^*(X)=\sum_{c,\kappa}V^c_\kappa\iota_c^*X\iota_cV^{c*}_\kappa+(1-Q)X(1-Q)`$, and for
$`A\in\mathfrak A`$ with blocks $`A_c`$, $`\iota_c^*A\iota_c=A_c`$. Since every $`F^d_\beta`$ is a combination of the $`V^d_\kappa`$,
$`(1-Q)F^d_\beta=0`$. With Step 2, for $`\alpha=(i,b)`$, $`\beta=(j,b')`$,
```math
\langle c;a,b|E_i^*\mathcal R^*(A)E_j|d;a',b'\rangle=\sum_{c',\kappa}\langle e_a,F^{c*}_\alpha V^{c'}_\kappa A_{c'}V^{c'*}_\kappa F^d_\beta e_{a'}\rangle
=\delta_{cd}\sum_\kappa u^c_{\alpha\kappa}d^c_\kappa\bar u^c_{\beta\kappa}(A_c)_{aa'}=\delta_{cd}(\Lambda^c_{ij})_{bb'}(A_c)_{aa'} .
```
Taking $`j=i`$ and summing over $`i`$: by the display with $`\sum_iE_i^*E_i=1`$, $`\sum_i\Lambda^c_{ii}=1_{m_c}`$, so
```math
\langle c;a,b|(\mathcal R\circ\Phi)^*(A)|d;a',b'\rangle=\delta_{cd}\,\delta_{bb'}(A_c)_{aa'}=\langle c;a,b|A|d;a',b'\rangle .
```
Both $`P(\mathcal R\circ\Phi)^*(A)P`$ and $`A`$ lie in $`B(\mathcal C)`$ and have the same matrix in the basis $`|c;a,b\rangle`$, so they are
equal. $`\square`$

*Computation:* Theorem D.9 and the recovery of its proof on a hybrid code in $`\mathbb C^8`$ — [tutorial, Section D.5](../tutorial/D_gksl.md#s-D-5).

<a id="r-D-10"></a>

**Corollary D.10.** (a) For $`\mathfrak A=B(\mathcal C)`$, $`\mathfrak A'=\mathbb CP`$, and Theorem D.9 is the Knill–Laflamme
theorem 10.5 ($`C=1`$, $`m_1=1`$; Step 3 is the recovery of Theorem 10.5). (b) For $`\mathfrak A=\mathrm{span}\{P_1,\dots,P_C\}`$ with
orthogonal projections $`P_c\ne0`$, $`\sum_cP_c=P`$, it is the equivalence of 2 and 3 in Proposition 10.6. (c) (*Noiseless
subsystems.*) In the notation of Proposition 10.1(a), let $`\mathcal C=V(\mathbb C^{k_c}\otimes\mathbb C^{m_c})`$ be the $`c`$-th
summand for the interaction algebra $`\mathcal A`$ of $`\Phi`$, and $`\mathfrak A=V(1_{k_c}\otimes M_{m_c})V^*`$ on $`\mathcal C`$, an
algebra of the form $`M_{m}\otimes1_{k}`$ after exchanging the two factors. Then $`\mathfrak A`$ is corrected by
$`\mathcal R=\mathrm{id}`$, and $`PE_i^*E_jP=V(E_{i,c}^*E_{j,c}\otimes1_{m_c})V^*\in\mathfrak A'`$ ($`P=VV^*`$).

*Proof.* (a) $`\mathfrak A'=\mathbb CP`$ is the case $`C=1`$, $`m_1=1`$ of the description of $`\mathfrak A'`$ above, and
$`PE_i^*E_jP\in\mathbb CP`$ is the Knill–Laflamme condition. Correctability of $`B(\mathcal C)`$ is correction in the sense of
Section 10.2: if $`\mathcal R(\Phi(\rho))=\rho`$ for all states $`\rho`$ supported in $`\mathcal C`$, then
$`\mathrm{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(A)P)=\mathrm{Tr}(\rho A)`$ for all of them, and these $`\rho`$ span $`B(\mathcal C)`$; the
converse is shown in Section 10.3. With $`C=1`$, $`m_1=1`$ the Kraus operators $`\iota_1V^{1*}_\kappa`$ of Step 3 are those
of Theorem 10.5, because $`\iota_1\iota_1^*=P`$. (b) $`T\in B(\mathcal C)`$ commutes with all $`P_c`$ iff
$`P_cTP_d=0`$ for $`c\ne d`$, which is condition 2 of Proposition 10.6. (c) Every $`A\in\mathfrak A`$, extended by $`0`$ to the
other summands, lies in $`\mathcal A'=V(\bigoplus_c1_{k_c}\otimes M_{m_c})V^*`$ (Theorem 2.6), so $`\Phi^*(A)=A`$ by Proposition
10.1(a) and $`PAP=A`$. The formula for $`PE_i^*E_jP`$ follows from $`V^*E_iV=\bigoplus_cE_{i,c}\otimes1_{m_c}`$, and
$`E_{i,c}^*E_{j,c}\otimes1_{m_c}`$ commutes with $`1_{k_c}\otimes M_{m_c}`$. $`\square`$

By Theorem 2.6 every $`\mathfrak A`$ is $`*`$-isomorphic to some $`\bigoplus_cM_{k_c}`$. When $`C\ge2`$ and some $`k_c\ge2`$, a
correctable $`\mathfrak A`$ protects the sector label $`c`$ (classical information) together with the quantum information of
each $`M_{k_c}`$; such codes are called **hybrid** codes.

The results of this appendix are due to V. Gorini, A. Kossakowski and E. C. G. Sudarshan
[*J. Math. Phys.* 17, 821 (1976)] and G. Lindblad [*Commun. Math. Phys.* 48, 119 (1976)] (Theorems D.3–D.5),
E. Christensen and D. E. Evans [*J. London Math. Soc.* 20, 358 (1979)] (Corollary D.8 in general von Neumann
algebras), A. Frigerio (Theorem D.6), and C. Bény, A. Kempf and D. W. Kribs [*Phys. Rev. Lett.* 98, 100502
(2007); *Phys. Rev. A* 76, 042303 (2007)] (Theorem D.9); the proofs above are complete for the finite-dimensional
statements used in the course.
