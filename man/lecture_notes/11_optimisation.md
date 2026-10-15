**English** · [Español](es/11_optimisation.md) · [Français](fr/11_optimisation.md) · [Deutsch](de/11_optimisation.md) · [中文](zh/11_optimisation.md) · [日本語](ja/11_optimisation.md) · [Русский](ru/11_optimisation.md)

# 11. Optimisation and trace inequalities

A trace inequality is a statement about the extreme values of a functional of operators. Section 3.7 showed one
that is plausible and false (Proposition 3.9); a sampling run ([tutorial, Section 3.7](../tutorial/03_traces_and_states.md#s-3-7)) did not find the violation numerically.
This chapter explains why sampling fails in high dimension and replaces it by optimisation over exactly
parametrised constraint sets; where the extremes are known in closed form, we prove them, so that the optimiser is
checked against a theorem.

**You will learn** how unitaries and spectra with a prescribed contrast are parametrised, and why every
admissible operator is reached; von Neumann's trace inequality and the exact extremes it gives; a trace inequality
with its equality condition; why Haar samples are, in mean square, far from the extremes; and what a multistart
search proves. The numerical checks of this chapter, and the optimiser of the
library, are in [Tutorial 11](../tutorial/11_optimisation.md).

Throughout, $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$ with $`D=\sum_ck_c\ge2`$, $`\mathrm{Tr}`$ is the
blunt trace (Section 3.1) and $`\|A\|_1=\sum_c\|A_c\|_1`$, the trace norm of $`A`$ as a block-diagonal matrix in
$`M_D(\mathbb C)`$. The singular values of $`A\in M_k(\mathbb C)`$ are $`s_1(A)\ge\dots\ge s_k(A)`$; $`x^\downarrow`$ and
$`x^\uparrow`$ are the decreasing and increasing rearrangements of a real vector $`x`$. $`\Delta`$ is the Michelson
contrast (Section 3.6), and for $`\delta\in[0,1)`$
```math
\ell(\delta)=\frac{1-\delta}{1+\delta},\qquad 1-\ell(\delta)=\frac{2\delta}{1+\delta};
```
a positive operator with $`\lambda_{\max}=1`$ has contrast $`\ell(\lambda_{\min})`$, so it has contrast $`\delta`$ iff
$`\lambda_{\min}=\ell(\delta)`$, since $`\ell`$ is a strictly decreasing involution of $`[0,1]`$. The examples are the three functionals of the companion paper
(section *From sampling to optimisation*):

| | functional | variables |
|---|---|---|
| Exp. 1 | $`z=\lvert\mathrm{Tr}(XUY)\rvert-\mathrm{Tr}(XY)`$ | $`X,Y\ge0`$, $`U\in\mathcal M`$ unitary |
| Exp. 2 | $`z=\mathrm{Tr}\lvert XY\rvert-\mathrm{Tr}(\lvert X\rvert\,\lvert Y\rvert)`$ | $`X=UX_0`$, $`Y=VY_0`$, $`X_0,Y_0\ge0`$, $`U,V`$ unitary |
| Exp. 3 | $`z=\mathrm{Tr}(Y\lvert X\rvert Y)-\mathrm{Tr}\lvert YXY\rvert`$ | $`X=X^*`$, $`Y\ge0`$ |

with $`|X|`$, $`|Y|`$ normalised to $`\lambda_{\max}=1`$ and their contrasts as parameters. The **envelopes** are
$`\overline z(\delta)=\sup\{z:\Delta(|X|)=\delta\}`$ and $`\underline z(\delta)=\inf\{z:\Delta(|X|)=\delta\}`$, all other
variables being free.

<a id="s-11-1"></a>

## 11.1 Constraint sets and their parametrisations

**Definition.** For $`\delta\in[0,1)`$ and $`\ell=\ell(\delta)`$, the **contrast shell** is
```math
\mathcal P_\delta=\{X\in\mathcal M:X\ge0,\ \lambda_{\max}(X)=1,\ \Delta(X)=\delta\}
=\{X=X^*:\mathrm{spec}\,X\subseteq[\ell,1],\ 1,\ell\in\mathrm{spec}\,X\}.
```
The **slots** of $`\mathcal M`$ are the $`D`$ pairs $`(c,i)`$, $`i\le k_c`$; for $`\lambda\in\mathbb R^D`$,
$`\mathrm{diag}(\lambda)=\bigoplus_c\mathrm{diag}(\lambda_{c,1},\dots,\lambda_{c,k_c})`$.

**The parametrisation.** For a prescribed contrast $`\delta`$, a search draws, per start, two distinct slots
$`s_\pm=(c_\pm,i_\pm)`$ ($`s_+`$ uniformly among the $`D`$ slots, then $`s_-`$ uniformly among the $`D-1`$ others, so every
ordered pair of distinct slots has probability $`1/(D(D-1))`$) and a Haar-random unitary $`U_0\in\mathcal M`$
(Lemma 11.10), and maps $`\theta\in\mathbb R^D`$ and
$`G=\bigoplus_cG_c`$ to
```math
X(\theta,G)=U\mathrm{diag}(\lambda)U^*,\qquad U=U_0\,e^{G-G^*},\qquad
\lambda_{s_+}=1,\quad\lambda_{s_-}=\ell,\quad\lambda_s=\ell+(1-\ell)\tfrac{1+\sin\theta_s}2\ \text{otherwise}.
```
Unitaries alone are parametrised by $`G\mapsto U_0e^{G-G^*}`$. Unlike a sigmoid, $`(1+\sin\theta)/2`$ reaches $`0`$ and $`1`$
at finite $`\theta`$, which is convenient because the extremal spectra of this chapter (Corollary 11.7,
Proposition 11.14) sit at $`\ell`$ and $`1`$. A free contrast in $`[a,b]`$ is one more parameter,
$`\delta=a+(b-a)(1+\sin\eta)/2`$, for instance with $`[a,b]=[0,0.999]`$. A value
attained with this restriction is also attained without it, so it is still a lower bound for an unrestricted
supremum (an upper bound for an unrestricted infimum).

*In the library:* this parametrisation is `optimize.PositiveParam`, and $`G\mapsto U_0e^{G-G^*}`$ is `optimize.UnitaryParam` — [tutorial, Section 11.1](../tutorial/11_optimisation.md#s-11-1).

<a id="r-11-1"></a>

**Lemma 11.1.** For every unitary $`U_0\in\mathcal M`$, $`G\mapsto U_0e^{G-G^*}`$ maps onto the unitary group of
$`\mathcal M`$.

*Proof.* Let $`U`$ be unitary. By Theorem 1.3 (blockwise) and Corollary 1.4(2),
$`U_0^*U=\bigoplus_cW_c\mathrm{diag}(e^{i\varphi_{c,j}})_jW_c^*`$ with $`W_c`$ unitary, $`\varphi_{c,j}\in\mathbb R`$. Then
$`S=\bigoplus_cW_c\mathrm{diag}(i\varphi_{c,j})_jW_c^*`$ satisfies $`S^*=-S`$ and $`e^S=U_0^*U`$ (the series
$`\sum_nS^n/n!`$ is computed term by term in the eigenbasis), and $`G=S/2`$ gives $`G-G^*=S`$. $`\square`$

<a id="r-11-2"></a>

**Proposition 11.2 (the parametrisation is onto).** For fixed $`\delta`$, $`s_\pm`$ and $`U_0`$, the image of
$`(\theta,G)\mapsto X(\theta,G)`$ is $`\{X\in\mathcal P_\delta:1\in\mathrm{spec}\,X_{c_+},\ \ell\in\mathrm{spec}\,X_{c_-}\}`$,
and $`\mathcal P_\delta`$ is the union of these images over all pairs of distinct slots.

*Proof.* The image lies in this set: $`X`$ has eigenvalues $`\lambda_s\in[\ell,1]`$, with $`\lambda_{s_\pm}`$ in blocks
$`c_\pm`$. Conversely, diagonalise $`X_c=W_c\mathrm{diag}(\mu_c)W_c^*`$ and permute each $`\mu_c`$ (replace $`W_c`$ by
$`W_c\Pi_c`$, $`\Pi_c`$ a permutation matrix) so that $`\mu_{s_+}=1`$ and $`\mu_{s_-}=\ell`$. This is possible: if $`c_+=c_-`$
and $`\delta>0`$, then $`1\ne\ell`$ are two eigenvalues of the same block and occupy different positions; if $`\delta=0`$,
all eigenvalues are $`1`$. If $`\ell<1`$, put $`\theta_s=\arcsin\big(2(\mu_s-\ell)/(1-\ell)-1\big)`$ for the other slots,
so that $`\lambda=\mu`$ (if $`\ell=1`$, any $`\theta`$). Lemma 11.1 gives $`G`$ with $`U_0e^{G-G^*}=\bigoplus_cW_c`$.
*Union:* for $`X\in\mathcal P_\delta`$ choose blocks $`c_\pm`$ containing the eigenvalues $`1`$, $`\ell`$ and distinct slots
in them; if $`c_+=c_-`$ and $`\delta>0`$ the block has two different eigenvalues, so $`k_{c_+}\ge2`$; if $`\delta=0`$, $`X=1`$
and any two slots will do. $`\square`$

*Computation:* Proposition 11.2 constructively: the parameters of the proof reproduce a target operator — [tutorial, Section 11.1](../tutorial/11_optimisation.md#s-11-1).

Since every pair of slots is drawn with positive probability, every point of the shell lies in the image of a
start with positive probability, and (in exact arithmetic) the constraint holds at every iterate.
For self-adjoint $`X`$ (Exp. 3) one uses $`X=U\mathrm{diag}(\sigma\lambda)U^*`$ with signs $`\sigma\in\{\pm1\}^D`$ drawn per start
(per block, the number of minus signs is uniform in $`\{0,\dots,k_c\}`$ and their positions are uniformly random, so
every sign pattern has positive probability) and changed only by the sign flips of Section 11.5. Then
$`|X|=U\mathrm{diag}(\lambda)U^*`$ (Theorem 1.5(e), as $`\lambda>0`$), and the union over signs and slots is
$`\{X=X^*:|X|\in\mathcal P_\delta\}`$: for such $`X`$, diagonalise $`X_c=W_c\mathrm{diag}(\nu_c)W_c^*`$, so that
$`|X|_c=W_c\mathrm{diag}(|\nu_c|)W_c^*`$ with $`|\nu_{c,i}|\ge\ell>0`$, and apply the proof of Proposition 11.2 to
$`\mu=|\nu|`$ with $`\sigma_s=\nu_s/|\nu_s|`$. Only $`S=G-G^*`$ matters.

*In the library:* `optimize.SelfAdjointParam` — [tutorial, Section 11.1](../tutorial/11_optimisation.md#s-11-1).

The chart $`S\mapsto U_0e^S`$ is well behaved near $`S=0`$ only. Let $`\mathfrak u=\{S\in\mathcal M:S^*=-S\}`$ and
$`\mathcal M_h=\{H\in\mathcal M:H^*=H\}`$, real subspaces with $`\mathcal M=\mathfrak u\oplus\mathcal M_h`$
($`X=\frac{X-X^*}2+\frac{X+X^*}2`$). *Near $`0`$, $`S\mapsto U_0e^S`$ is a homeomorphism of an open ball $`N\ni0`$ of
$`\mathfrak u`$ onto a neighbourhood of $`U_0`$ in the unitary group.* Indeed, $`\Psi(S,H)=U_0e^Se^H`$ is a $`C^1`$ map
$`\mathfrak u\times\mathcal M_h\to\mathcal M`$ (the series $`\sum_{m,n\ge0}U_0S^mH^n/(m!\,n!)`$ and the series of the
partial derivatives of its terms with respect to the real coordinates converge uniformly on bounded sets), and
$`e^{tS}e^{tH}=1+t(S+H)+O(t^2)`$ shows that its differential at $`0`$ is $`(S,H)\mapsto U_0(S+H)`$, a bijection onto
$`\mathcal M`$. By the inverse function theorem there are open balls $`N\ni0`$ in $`\mathfrak u`$ and $`N_h\ni0`$ in
$`\mathcal M_h`$ such that $`\Psi`$ is a homeomorphism of $`N\times N_h`$ onto an open set $`O\subseteq\mathcal M`$. If
$`\Psi(S,H)`$ is unitary, then so is $`e^H=e^{-S}U_0^*\Psi(S,H)`$ (in an eigenbasis of the normal matrix $`S`$, whose
eigenvalues are imaginary, $`e^{\pm S}`$ are diagonal with entries of modulus $`1`$ and inverse to each other); but $`e^H`$
has the positive eigenvalues $`e^{h_i}`$, $`h_i`$ the eigenvalues of $`H`$, so $`e^{h_i}=1`$ (Corollary 1.4(2)) and $`H=0`$.
Hence $`O\cap U(\mathcal M)=U_0e^N`$, where $`U(\mathcal M)`$ is the unitary group of $`\mathcal M`$: a neighbourhood of
$`U_0`$ in $`U(\mathcal M)`$, and $`S\mapsto U_0e^S`$ is a continuous bijection $`N\to U_0e^N`$ whose inverse, the first
component of $`\Psi^{-1}`$, is continuous. Far from $`0`$ the chart is not injective: $`e^{2\pi iP}=1=e^0`$ for every
projection $`P\in\mathcal M`$. A search therefore **re-centres** the chart after every round ($`U_0\leftarrow U_0e^{G-G^*}`$,
$`G\leftarrow0`$). Two other charts have the same differential $`S\mapsto U_0S`$ at $`0`$. The Cayley chart
$`S\mapsto U_0C(-S/2)`$ (Exercise 1): $`(1+tS/2)(1-tS/2)^{-1}=(1+tS/2)(1+tS/2+O(t^2))=1+tS+O(t^2)`$. The polar retraction,
$`S\mapsto`$ the unitary factor of $`A=U_0(1+S)`$ (Theorem 1.8): $`A^*A=(1-S)(1+S)=1+S^*S\ge1`$ is invertible, so $`|A|`$ is
invertible and the unitary factor is $`A|A|^{-1}`$; for $`tS`$ in place of $`S`$, $`|A|^{-1}=(1+t^2S^*S)^{-1/2}=1+O(t^2)`$
(in an eigenbasis of $`S^*S`$), so the factor is $`U_0(1+tS)+O(t^2)`$.

> **Physical meaning** (physical interpretation, not a theorem of this course). The Michelson contrast is the
> visibility $`(I_{\max}-I_{\min})/(I_{\max}+I_{\min})`$ of interference fringes. For $`0\le X\le1`$ with
> $`\lambda_{\max}=1`$ — a filter, an attenuator, a detector efficiency — it fixes the dynamic range
> $`\lambda_{\min}/\lambda_{\max}`$, the inverse of a condition number.

<a id="s-11-2"></a>

## 11.2 Von Neumann's trace inequality

A matrix $`(c_{ij})`$ with $`c_{ij}\ge0`$ and all row and column sums $`1`$ is **doubly stochastic**; for a unitary $`G`$,
$`(|G_{ij}|^2)`$ and $`(|G_{ji}|^2)`$ are doubly stochastic, as rows and columns of $`G`$ are unit vectors.

<a id="r-11-3"></a>

**Lemma 11.3 (rearrangement).** Let $`p_1\ge\dots\ge p_k`$, $`q_1\ge\dots\ge q_k`$ be real and $`(c_{ij})`$ doubly
stochastic. Then $`\sum_ip_iq_{k+1-i}\le\sum_{i,j}p_iq_jc_{ij}\le\sum_ip_iq_i`$.

*Proof.* *Upper bound.* For $`m\le k`$ let $`r_j=\sum_{i\le m}c_{ij}\in[0,1]`$, so $`\sum_jr_j=m`$. Since $`q_j\ge q_m`$ for
$`j\le m`$ and $`q_j\le q_m`$ for $`j>m`$,
```math
\sum_{j\le m}q_j-\sum_jq_jr_j=\sum_{j\le m}q_j(1-r_j)-\sum_{j>m}q_jr_j\ge q_m\Big(m-\sum_jr_j\Big)=0 .
```
With $`x_i=\sum_jq_jc_{ij}`$ this says $`\sum_{i\le m}x_i\le\sum_{j\le m}q_j`$, with equality for $`m=k`$. By Abel
summation and $`p_m-p_{m+1}\ge0`$,
```math
\sum_ip_ix_i=p_k\sum_{i\le k}x_i+\sum_{m<k}(p_m-p_{m+1})\sum_{i\le m}x_i
\le p_k\sum_{j\le k}q_j+\sum_{m<k}(p_m-p_{m+1})\sum_{j\le m}q_j=\sum_ip_iq_i .
```
*Lower bound.* Apply the upper bound to $`p`$, the decreasing $`q'_j=-q_{k+1-j}`$ and the doubly stochastic
$`c'_{ij}=c_{i,k+1-j}`$: $`-\sum_{ij}p_iq_jc_{ij}=\sum_{ij}p_iq'_jc'_{ij}\le\sum_ip_iq'_i=-\sum_ip_iq_{k+1-i}`$. $`\square`$

Every $`A\in M_k(\mathbb C)`$ has a **singular value decomposition** $`A=W\Sigma V^*`$, $`W,V`$ unitary,
$`\Sigma=\mathrm{diag}(s_i(A))`$: write $`A=W'|A|`$ (Theorem 1.8), $`|A|=V\Sigma V^*`$ (Theorem 1.3), $`W=W'V`$.

<a id="r-11-4"></a>

**Theorem 11.4 (von Neumann's trace inequality).** For $`A,B\in M_k(\mathbb C)`$,
$`|\mathrm{Tr}(AB)|\le\sum_is_i(A)\,s_i(B)`$.

*Proof.* Let $`A=W_A\Sigma_AV_A^*`$, $`B=W_B\Sigma_BV_B^*`$, $`a_i=s_i(A)`$, $`b_j=s_j(B)`$, and $`G=V_A^*W_B`$, $`H=V_B^*W_A`$.
By cyclicity $`\mathrm{Tr}(AB)=\mathrm{Tr}(\Sigma_AG\Sigma_BH)=\sum_{i,j}a_iG_{ij}b_jH_{ji}`$, so, by
$`|g||h|\le\frac12(|g|^2+|h|^2)`$ and $`a_ib_j\ge0`$,
```math
|\mathrm{Tr}(AB)|\le\sum_{i,j}a_ib_j\,d_{ij},\qquad d_{ij}=\tfrac12\big(|G_{ij}|^2+|H_{ji}|^2\big).
```
$`(d_{ij})`$ is doubly stochastic, and Lemma 11.3 gives $`\sum_{ij}a_ib_jd_{ij}\le\sum_ia_ib_i`$. $`\square`$

<a id="r-11-5"></a>

**Corollary 11.5.** (a) For self-adjoint $`P,Q\in M_k(\mathbb C)`$ with eigenvalues $`p,q`$:
$`\sum_ip^\downarrow_iq^\uparrow_i\le\mathrm{Tr}(PQ)\le\sum_ip^\downarrow_iq^\downarrow_i`$.
(b) For $`P,Q\ge0`$: $`\mathrm{Tr}(PQ)\le\|PQ\|_1\le\sum_ip^\downarrow_iq^\downarrow_i`$.

*Proof.* (a) As in the proof of Theorem 3.5, $`\mathrm{Tr}(PQ)=\sum_{ij}p_iq_jc_{ij}`$ with the doubly stochastic
$`c_{ij}=|\langle u_i,v_j\rangle|^2`$ for eigenbases ordered decreasingly; apply Lemma 11.3 (its lower bound is
$`\sum_ip^\downarrow_iq^\downarrow_{k+1-i}=\sum_ip^\downarrow_iq^\uparrow_i`$). (b) The first inequality is
Theorem 1.9. By Theorem 1.9 and Theorem 11.4, $`\|PQ\|_1=\max_U|\mathrm{Tr}(P\,QU)|\le\sum_is_i(P)s_i(QU)`$, and
$`s_i(QU)=q^\downarrow_i`$ because $`|QU|=(U^*Q^2U)^{1/2}=U^*QU`$ (Theorem 1.6). $`\square`$

*Computation:* Theorem 11.4 and Corollary 11.5 on random matrices — [tutorial, Section 11.2](../tutorial/11_optimisation.md#s-11-2).

So for $`X,Y\ge0`$ in $`M_k(\mathbb C)`$ with spectra $`x,y`$, the violation $`\|YX\|_1-\mathrm{Tr}(XY)`$ of
Proposition 3.9 is at most $`\sum_ix^\downarrow_i(y^\downarrow_i-y^\uparrow_i)`$ (Corollary 11.5(b) for $`(Y,X)`$ and
the lower bound of (a)). For Exp. 2 this number is the exact answer:

<a id="r-11-6"></a>

**Proposition 11.6 (Exp. 2 for fixed spectra).** In Exp. 2, $`|X|=X_0`$, $`|Y|=Y_0`$ and
$`z=\|X_0VY_0\|_1-\mathrm{Tr}(X_0Y_0)`$. If the spectra $`x_c`$ of $`(X_0)_c`$ and $`y_c`$ of $`(Y_0)_c`$ are fixed and the
eigenbases and $`V`$ vary, then
```math
\sup z=-\inf z=\sum_c\sum_ix^\downarrow_{c,i}\big(y^\downarrow_{c,i}-y^\uparrow_{c,i}\big),
```
and both are attained.

*Proof.* $`X^*X=X_0U^*UX_0=X_0^2`$, so $`|X|=X_0`$ (Theorem 1.6); likewise $`|Y|=Y_0`$, and $`|UA|=|A|`$ gives
$`\mathrm{Tr}|XY|=\|X_0VY_0\|_1`$. All variables act blockwise and $`z`$ is a sum over blocks, so let $`C=1`$.
Put $`Y'=VY_0V^*\ge0`$ (spectrum $`y`$); since $`|AV|=V^*|A|V`$, $`\|X_0VY_0\|_1=\|X_0Y'V\|_1=\|X_0Y'\|_1`$. Corollary 11.5
for the pairs $`(X_0,Y')`$ and $`(X_0,Y_0)`$, with $`\|X_0Y'\|_1\ge\mathrm{Tr}(X_0Y')`$, gives
```math
\sum_ix^\downarrow_iy^\uparrow_i-\sum_ix^\downarrow_iy^\downarrow_i\le\|X_0Y'\|_1-\mathrm{Tr}(X_0Y_0)\le\sum_ix^\downarrow_iy^\downarrow_i-\sum_ix^\downarrow_iy^\uparrow_i .
```
Let $`J`$ be the reversal permutation, $`Je_i=e_{k+1-i}`$; then $`J\mathrm{diag}(w)=\mathrm{diag}(w_k,\dots,w_1)J`$, and
the trace norm of $`\mathrm{diag}(w)J`$ is $`\sum_i|w_i|`$. With $`X_0=\mathrm{diag}(x^\downarrow)`$ and $`V=J`$:
for $`Y_0=\mathrm{diag}(y^\uparrow)`$, $`X_0VY_0=\mathrm{diag}(x^\downarrow_iy^\downarrow_i)J`$ and
$`\mathrm{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\uparrow_i`$ (the upper bound); for $`Y_0=\mathrm{diag}(y^\downarrow)`$,
$`X_0VY_0=\mathrm{diag}(x^\downarrow_iy^\uparrow_i)J`$ and $`\mathrm{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\downarrow_i`$
(the lower bound). $`\square`$

<a id="r-11-7"></a>

**Corollary 11.7 (envelope of Exp. 2).** Let $`\ell_X,\ell_Y\in(0,1]`$ and some $`k_c\ge2`$. Over all $`X_0,Y_0\ge0`$
with $`\lambda_{\max}=1`$, $`\lambda_{\min}(X_0)=\ell_X`$, $`\lambda_{\min}(Y_0)=\ell_Y`$, and all unitaries,
```math
\sup z=-\inf z=\sum_c\Big\lfloor\frac{k_c}2\Big\rfloor(1-\ell_X)(1-\ell_Y).
```

*Proof.* In one block write $`x_i=x^\downarrow_{c,i}`$, $`y_i=y^\downarrow_{c,i}`$, $`k=k_c`$, so $`y^\uparrow_{c,i}=y_{k+1-i}`$.
Pairing $`i`$ with $`k+1-i`$, $`x_i(y_i-y_{k+1-i})+x_{k+1-i}(y_{k+1-i}-y_i)=(x_i-x_{k+1-i})(y_i-y_{k+1-i})`$, and for odd
$`k`$ the middle index contributes $`x_i(y_i-y_i)=0`$. Hence
```math
\sum_ix_i(y_i-y_{k+1-i})=\sum_{i\le\lfloor k/2\rfloor}(x_i-x_{k+1-i})(y_i-y_{k+1-i})\le\Big\lfloor\frac k2\Big\rfloor(1-\ell_X)(1-\ell_Y),
```
as $`x_i\in[\ell_X,1]`$, $`y_i\in[\ell_Y,1]`$ are decreasing, so both factors lie in $`[0,1-\ell_X]`$ and $`[0,1-\ell_Y]`$.
Equality holds if in every block the $`\lfloor k_c/2\rfloor`$ largest eigenvalues are $`1`$ and the $`\lfloor k_c/2\rfloor`$
smallest are $`\ell_X`$ (resp. $`\ell_Y`$), the middle one (odd $`k_c`$) being arbitrary in $`[\ell_X,1]`$ (resp.
$`[\ell_Y,1]`$); these spectra have global extremes $`1`$ and $`\ell_X`$ ($`\ell_Y`$) because some $`k_c\ge2`$. Proposition 11.6
concludes. $`\square`$

Here $`1-\ell_X=2\delta/(1+\delta)`$ for $`\Delta(|X|)=\delta`$; for odd $`k_c`$ the middle eigenvalue is irrelevant.

*Computation:* the extremal configurations of Proposition 11.6 and the envelope of Corollary 11.7 for $`k=3`$, $`C=2`$ — [tutorial, Section 11.2](../tutorial/11_optimisation.md#s-11-2); the envelope found by a multistart search — [tutorial, Section 11.5](../tutorial/11_optimisation.md#s-11-5).

<a id="s-11-3"></a>

## 11.3 A trace inequality with an equality condition

For $`X=X^*`$ let $`X_\pm=f_\pm(X)`$ with $`f_+(t)=\max(t,0)`$, $`f_-(t)=\max(-t,0)`$. Then $`X_\pm\ge0`$, $`X=X_+-X_-`$,
$`|X|=X_++X_-`$ and $`X_+X_-=0`$ (Theorem 1.5(c), and 1.5(e) for $`|X|=(X^2)^{1/2}=h(X)`$ with
$`h(t)=\sqrt{t^2}=|t|=f_+(t)+f_-(t)`$; Corollary 1.4(3) for positivity).

<a id="r-11-8"></a>

**Lemma 11.8.** For $`A,B\ge0`$ in $`M_n(\mathbb C)`$, $`\|A-B\|_1\le\mathrm{Tr}\,A+\mathrm{Tr}\,B`$, with equality iff
$`AB=0`$.

*Proof.* Let $`H=A-B`$ and $`Q=g(H)`$ with $`g=1`$ on $`(0,\infty)`$, $`g=0`$ elsewhere; $`Q`$ is a projection, as $`g=\bar g=g^2`$
(Theorem 1.5(c)). As $`|t|=t(2g(t)-1)`$, $`|H|=H(2Q-1)`$ (Theorem 1.5(c),(e)), and, using
$`\mathrm{Tr}(AR)=\mathrm{Tr}(RAR)`$ for the projections $`R=Q,1-Q`$ (cyclicity),
```math
\|H\|_1=\mathrm{Tr}\big(A(2Q-1)\big)-\mathrm{Tr}\big(B(2Q-1)\big)
=\mathrm{Tr}\,A+\mathrm{Tr}\,B-2\mathrm{Tr}\big((1-Q)A(1-Q)\big)-2\mathrm{Tr}(QBQ).
```
The subtracted traces are $`\ge0`$, and $`\mathrm{Tr}(RAR)=\|A^{1/2}R\|_2^2`$ vanishes iff $`AR=0`$; so equality holds
iff $`A(1-Q)=0`$ and $`BQ=0`$, and then $`AB=AQB=A(BQ)^*=0`$. Conversely, if $`AB=0`$ then $`BA=0`$,
$`(A-B)^2=A^2+B^2=(A+B)^2`$, and Theorem 1.6 gives $`|A-B|=A+B`$. $`\square`$

<a id="r-11-9"></a>

**Proposition 11.9 (Exp. 3).** For $`X=X^*`$ and $`Y\ge0`$ in $`\mathcal M`$,
$`\mathrm{Tr}|YXY|\le\mathrm{Tr}(Y|X|Y)`$, with equality iff $`YX_+Y^2X_-Y=0`$, in particular if $`XY=YX`$.

*Proof.* In $`M_D(\mathbb C)`$, $`A=YX_+Y`$ and $`B=YX_-Y`$ are positive (Corollary 1.7 with $`T=X_\pm^{1/2}Y`$), $`YXY=A-B`$
and $`\mathrm{Tr}(Y|X|Y)=\mathrm{Tr}\,A+\mathrm{Tr}\,B`$; apply Lemma 11.8. If $`XY=YX`$, then $`Y`$ commutes with
$`X_\pm`$ (Theorem 1.5(d)), and $`AB=Y^4X_+X_-=0`$. $`\square`$

So $`z\ge0`$ in Exp. 3, and $`\underline z=0`$ (take $`Y=1`$). The equality condition is weaker than commutation: for
$`X=\mathrm{diag}(1,-1,0)`$ and $`Y`$ with rows $`(1,0,a)`$, $`(0,1,0)`$, $`(a,0,1)`$, $`-1<a<1`$,
$`AB=(Ye_1)\langle Ye_1,Ye_2\rangle(Ye_2)^*=0`$, while $`XY\ne YX`$ for $`a\ne0`$.

*Computation:* Proposition 11.9 on random pairs, and the example of equality without commutation — [tutorial, Section 11.3](../tutorial/11_optimisation.md#s-11-3).

<a id="s-11-4"></a>

## 11.4 Sampling versus optimisation

How large is $`|\mathrm{Tr}(AU)|`$ typically, compared with its maximum $`\|A\|_1`$? We sample $`U`$ as
the factor $`Q`$ of the QR decomposition $`Z=QR`$, $`R_{ii}>0`$, of a **Ginibre matrix** $`Z`$ (independent entries with
density $`\pi^{-1}e^{-|z|^2}`$; Mezzadri's recipe), independently in every block.

<a id="r-11-10"></a>

**Lemma 11.10.** $`Q`$ is Haar distributed on $`U(k)`$; in particular, for every fixed unitary $`V\in M_k(\mathbb C)`$, $`VQ`$ has
the law of $`Q`$. The unitary $`\bigoplus_cQ_c`$ with independent blocks is Haar distributed on the unitary group of
$`\mathcal M`$.

*Proof.* This is Theorem 5.3 (the left invariance is the equivariance step of its proof; the law is identified
with the Haar measure by the uniqueness Lemma 5.1, whose measure-theoretic step is Corollary A.2) and
Corollary 5.4. $`\square`$

*In the library:* this sampler is `TypeIAlgebra.random_unitary_operator` — [tutorial, Section 11.4](../tutorial/11_optimisation.md#s-11-4).

<a id="r-11-11"></a>

**Proposition 11.11.** For Haar $`U\in U(k)`$ and $`A\in M_k(\mathbb C)`$: $`\mathbb E\mathrm{Tr}(AU)=0`$ and
$`\mathbb E|\mathrm{Tr}(AU)|^2=\|A\|_2^2/k`$.

*Proof.* $`\mathbb EU_{ij}=0`$ by Corollary 5.5. Proposition 5.7 (with $`C=1`$) for $`X=|e_j\rangle\langle e_m|`$ gives the
Weingarten formula of order one,
$`\mathbb E[U_{ij}\bar U_{lm}]=\mathbb E(UXU^*)_{il}=\delta_{il}\delta_{jm}/k`$; for $`k\ge2`$ it also follows from the
order-two formula, Theorem 5.10, by setting $`i_2=i'_2`$, $`j_2=j'_2`$ and summing over both, since
$`\sum_{i_2,j_2}|U_{i_2j_2}|^2=k`$. With $`\mathrm{Tr}(AU)=\sum_{ij}A_{ji}U_{ij}`$,
$`\mathbb E|\mathrm{Tr}(AU)|^2=\sum_{ijlm}A_{ji}\bar A_{ml}\mathbb E[U_{ij}\bar U_{lm}]=\frac1k\sum_{ij}|A_{ji}|^2`$. $`\square`$

<a id="r-11-12"></a>

**Corollary 11.12.** Let $`X,Y\ge0`$ in $`M_k(\mathbb C)`$ have spectra in $`[a,b]`$, $`a>0`$, and $`U`$ be Haar. Then
```math
\frac{\big(\mathbb E|\mathrm{Tr}(XUY)|^2\big)^{1/2}}{\max_V|\mathrm{Tr}(XVY)|}=\frac{\|YX\|_2}{\sqrt k\,\|YX\|_1}\le\min\Big(\frac1{\sqrt k},\frac{b^2}{a^2k}\Big),
\qquad
\mathbb P\big(|\mathrm{Tr}(XUY)|\ge\mathrm{Tr}(XY)\big)\le\frac{b^4}{a^4k^2}.
```

*Proof.* The equality is Proposition 11.11 for $`A=YX`$ (as $`\mathrm{Tr}(XUY)=\mathrm{Tr}(AU)`$) and
Proposition 3.9. As $`\|A\|_2^2=\sum_is_i^2\le(\sum_is_i)^2=\|A\|_1^2`$ (Section 1.4), the ratio is $`\le1/\sqrt k`$.
Since $`\|Av\|=\||A|v\|`$ (Section 1.4), diagonalising $`|A|`$ gives $`s_k(A)=\min_{\|v\|=1}\|Av\|`$ and
$`s_1(A)=\max_{\|v\|=1}\|Av\|`$; and $`\|Xv\|^2=\langle v,X^2v\rangle\in[a^2,b^2]`$ for unit $`v`$ (Corollary 1.4(3)
for $`X^2-a^2`$ and $`b^2-X^2`$), likewise for $`Y`$. Hence $`a^2\le s_k(A)\le s_1(A)\le b^2`$, $`\|A\|_2^2\le kb^4`$ and
$`\|A\|_1\ge ka^2`$. Finally $`\mathrm{Tr}(XY)=\mathrm{Tr}(X^{1/2}YX^{1/2})\ge a\mathrm{Tr}\,X\ge ka^2`$, and
Markov's inequality for $`|\mathrm{Tr}(XUY)|^2`$ bounds the probability by $`(kb^4/k)/(ka^2)^2`$. $`\square`$

So $`N`$ Haar samples find a violation of $`|\mathrm{Tr}(XUY)|\le\mathrm{Tr}(XY)`$ with probability at most
$`Nb^4/(a^4k^2)`$ (a union of $`N`$ events), although one exists for every non-commuting pair (Proposition 3.9). The bound
is informative only for $`k>\sqrt N\,b^2/a^2`$ (for $`b/a=6`$ and $`N=2000`$: $`k>1609`$).

*Computation:* Proposition 11.11 by Monte Carlo, and the ratio of Corollary 11.12 against the best of $`2000`$ samples for $`k\le32`$ — [tutorial, Section 11.4](../tutorial/11_optimisation.md#s-11-4).

<a id="r-11-13"></a>

**Lemma 11.13.** If $`k\ge2`$, then $`\min_U|\mathrm{Tr}(AU)|=0`$ for every $`A\in M_k(\mathbb C)`$.

*Proof.* Let $`A=W\Sigma V^*`$ and $`P`$ the cyclic shift, $`Pe_i=e_{(i\bmod k)+1}`$, whose diagonal vanishes for $`k\ge2`$. For
$`U=VPW^*`$, $`\mathrm{Tr}(AU)=\mathrm{Tr}(\Sigma P)=0`$. $`\square`$

So $`\inf_Uz=-\mathrm{Tr}(XY)`$ and $`\sup_Uz=\|YX\|_1-\mathrm{Tr}(XY)`$ (Proposition 3.9). By Proposition 11.11 the root mean square of
$`|\mathrm{Tr}(XUY)|`$ is $`\|YX\|_2/\sqrt k`$; with the bounds in the proof of Corollary 11.12 (spectra in $`[a,b]`$),
$`\mathbb E|\mathrm{Tr}(XUY)|^2\le b^4`$ and $`\mathrm{Tr}(XY)\ge ka^2`$, so Markov's inequality gives
$`\mathbb P\big(z\ge-(1-\varepsilon)\mathrm{Tr}(XY)\big)\le b^4/(\varepsilon^2a^4k^2)`$ for $`\varepsilon>0`$: for large $`k`$,
$`z`$ is close to its infimum with high probability.

*Computation:* both extremes over $`U`$ by gradient ascent and descent, against $`20\,000`$ Haar samples — [tutorial, Section 11.4](../tutorial/11_optimisation.md#s-11-4).

> **Physical meaning** (physical interpretation, not a theorem of this course). Corollary 11.12 is a weak, provable
> instance of what statistical mechanics calls **typicality**: in high dimension almost every state or evolution is
> expected to give nearly the same value of a macroscopic quantity. A sharp inequality is decided by fine-tuned
> configurations — aligned eigenbases, two-level spectra — that a random unitary (a random circuit, a random
> scattering matrix) rarely produces.

<a id="s-11-5"></a>

## 11.5 Multistart search for envelopes

A **multistart search** runs a gradient method, in rounds, from many independent starts, each with its own
slots, signs, $`U_0`$ and contrast; it re-centres the unitary charts after every round and returns, per start, the
best value of the objective at the iterates. The signs of the self-adjoint parametrisation are discrete; between
rounds they are changed by single sign flips at random slots, each kept if it improves the objective.

*In the library:* `optimize.extremize` — [tutorial, Section 11.5](../tutorial/11_optimisation.md#s-11-5).

In exact arithmetic every iterate is admissible (Proposition 11.2 and the remarks after it), so every returned value
is the value of $`z`$ at an admissible point; in floating point it is that value up to rounding. Hence a maximisation gives a **lower bound** for the supremum, a minimisation an upper bound for
the infimum, and nothing more — a local search can stop at a local extremum. The other direction needs a proof,
like Corollary 11.7 and Proposition 11.9, or a convex relaxation (Section 11.6). A computed envelope is checked
(i) against closed forms; (ii) against configurations built by hand, such as direct sums of smaller solutions, which
it must not fall below; (iii) against structured samples (two-level spectra, permutation bases), whose values are
further attained values.

<a id="r-11-14"></a>

**Proposition 11.14 (lower envelope of Exp. 1).** Let all $`k_c\ge2`$. Over $`X\in\mathcal P_\delta`$, $`Y\ge0`$ with
$`\lambda_{\max}(Y)=1`$ and unitaries $`U\in\mathcal M`$, $`\ \min z=-\big(D-2\delta/(1+\delta)\big)`$.

*Proof.* $`Y\le1`$ gives $`X^{1/2}YX^{1/2}\le X`$, so $`\mathrm{Tr}(XY)\le\mathrm{Tr}\,X\le(D-1)+\ell(\delta)`$, as one
eigenvalue of $`X`$ is $`\ell(\delta)`$ and the others are $`\le1`$; hence $`z\ge-\mathrm{Tr}(XY)\ge-(D-2\delta/(1+\delta))`$.
Equality holds for $`X`$ with eigenvalues $`\ell(\delta)`$ once and $`1`$ otherwise, $`Y=1`$, and $`U=\bigoplus_cU_c`$ with
$`\mathrm{Tr}(X_cU_c)=0`$ (Lemma 11.13). $`\square`$

Checks of type (ii) can use direct sums: $`z`$ of Exp. 3 is additive over direct sums, and the contrast of a direct sum is the largest
contrast of its summands (Exercise 3), so, for the envelopes of Exp. 3 with $`\Delta(Y)`$ free in a fixed interval $`[0,b]`$,
```math
\overline z_{M_2\oplus M_2}(\delta)\ \ge\ \overline z_{M_2}(\delta)+\sup_{\delta'\le\delta}\overline z_{M_2}(\delta');
```
in particular the sum of two values attained on $`M_2`$ at contrasts $`\delta`$ and $`\delta'\le\delta`$ is a value
attained on $`M_2\oplus M_2`$ at contrast $`\delta`$.

*Computation:* Proposition 11.14 on $`M_2\oplus M_2`$, Corollary 11.7 on $`M_3\oplus M_3`$ and Proposition 11.9 by multistart search, with a check of type (ii) — [tutorial, Section 11.5](../tutorial/11_optimisation.md#s-11-5).

For blocks of size $`2`$ the upper envelopes of Exp. 1 and Exp. 3 are known exactly. In $`M_2(\mathbb C)`$ the
eigenbases enter only through one number: for orthonormal eigenbases $`e_1,e_2`$ and $`f_1,f_2`$ the doubly stochastic
matrix $`(|\langle e_i,f_j\rangle|^2)`$ has rows $`(\gamma,1-\gamma)`$ and $`(1-\gamma,\gamma)`$, $`\gamma=|\langle e_1,f_1\rangle|^2\in[0,1]`$, and
every $`\gamma`$ occurs for real bases ($`f_1=\cos t\,e_1+\sin t\,e_2`$, $`f_2=-\sin t\,e_1+\cos t\,e_2`$, $`\gamma=\cos^2t`$).

<a id="r-11-15"></a>

**Proposition 11.15 (Exp. 1 and Exp. 3 on $`M_2`$).** Let $`Y\ge0`$ in $`M_2(\mathbb C)`$ have eigenvalues
$`y_1\ge y_2\ge0`$ with eigenvectors $`f_1,f_2`$, and $`y_1>0`$.

(a) Let $`X\ge0`$ have eigenvalues $`x_1\ge x_2\ge0`$, $`x_1>0`$, with eigenvectors $`e_1,e_2`$. For every unitary $`U`$,
```math
\lvert\mathrm{Tr}(XUY)\rvert-\mathrm{Tr}(XY)\ \le\ \frac{(x_1-x_2)^2(y_1-y_2)^2}{4(x_1+x_2)(y_1+y_2)} .
```
If $`x_1>x_2`$ and $`y_1>y_2`$, equality holds iff $`|\mathrm{Tr}(XUY)|=\|YX\|_1`$ and
$`\gamma=\gamma_1:=\frac14+\frac{x_1y_2+x_2y_1}{2(x_1+x_2)(y_1+y_2)}`$; otherwise the bound is $`0`$, and it is attained iff
$`|\mathrm{Tr}(XUY)|=\|YX\|_1`$.

(b) Let $`X=X^*`$, and let $`|X|`$ have eigenvalues $`x_1\ge x_2\ge0`$, $`x_1>0`$. Then
```math
\mathrm{Tr}(Y|X|Y)-\mathrm{Tr}|YXY|\ \le\ \frac{2x_1x_2}{x_1+x_2}\,(y_1-y_2)^2 .
```
If $`X`$ has eigenvalues of both signs and $`y_1>y_2`$, let $`e_1`$, $`e_2`$ be eigenvectors of $`X`$ for its eigenvalues
of modulus $`x_1`$ and $`x_2`$; then equality holds iff $`\gamma=\gamma_3:=\frac{x_1y_2+x_2y_1}{(x_1+x_2)(y_1+y_2)}`$.

For all eigenvalues as above both bounds are attained, by real $`X`$, $`Y`$ and $`U`$ (in (b) with $`X=\mathrm{diag}(x_1,-x_2)`$); in
particular they are the same in $`M_2(\mathbb R)`$.

*Proof.* Write $`z`$ for the left side, and put $`S=x_1y_2+x_2y_1`$, $`p=(x_1-x_2)(y_1-y_2)\ge0`$ and
$`q=(x_1+x_2)(y_1+y_2)>0`$, so that $`q-2S=p`$.

(a) By Theorem 1.9, $`|\mathrm{Tr}(XUY)|=|\mathrm{Tr}(YX\,U)|\le\|YX\|_1`$, with equality for $`U=W^*`$, $`YX=W|YX|`$.
As in the proof of Corollary 11.5, $`\mathrm{Tr}(PQ)=\sum_{i,j}a_ib_j|\langle e_i,f_j\rangle|^2`$ for $`P=\sum_ia_i|e_i\rangle\langle e_i|`$,
$`Q=\sum_jb_j|f_j\rangle\langle f_j|`$; hence
$`\mathrm{Tr}(XY)=(x_1y_2+x_2y_1)+\gamma(x_1-x_2)(y_1-y_2)=S+\gamma p`$ and
$`\mathrm{Tr}(X^2Y^2)=(x_1^2y_2^2+x_2^2y_1^2)+\gamma(x_1^2-x_2^2)(y_1^2-y_2^2)`$. The squares $`s_1^2,s_2^2`$ of the singular
values of $`YX`$ are the eigenvalues of $`(YX)^*YX=XY^2X`$, so $`s_1^2+s_2^2=\mathrm{Tr}(X^2Y^2)`$ and
$`s_1s_2=|\det(YX)|=x_1x_2y_1y_2`$. Therefore
```math
\|YX\|_1^2=(s_1+s_2)^2=\mathrm{Tr}(X^2Y^2)+2x_1x_2y_1y_2=S^2+\gamma pq .
```
For $`r\ge0`$, $`\sqrt r\le r/q+q/4`$, with equality iff $`\sqrt r=q/2`$ (expand $`(\sqrt r-q/2)^2\ge0`$). With $`r=S^2+\gamma pq`$,
```math
z\ \le\ \|YX\|_1-\mathrm{Tr}(XY)\ \le\ \frac{S^2+\gamma pq}q+\frac q4-S-\gamma p=\frac{(q/2-S)^2}q=\frac{p^2}{4q}.
```
Equality holds iff $`|\mathrm{Tr}(XUY)|=\|YX\|_1`$ and $`S^2+\gamma pq=q^2/4`$, i.e. $`\gamma pq=(q/2-S)(q/2+S)=\frac p2(q/2+S)`$.
If $`p>0`$ this says $`\gamma=(q/2+S)/(2q)=\gamma_1`$, which lies in $`[\frac14,\frac12]`$ as $`0\le S\le q/2`$; if $`p=0`$ it holds for
every $`\gamma`$.

(b) Replacing $`X`$ by $`-X`$ changes neither side. If $`X\ge0`$, then $`YXY\ge0`$, $`|YXY|=YXY`$ and $`z=0`$;
likewise for $`X\le0`$. Otherwise $`x_2>0`$ and, after replacing $`X`$ by $`-X`$ if necessary, $`Xe_1=x_1e_1`$,
$`Xe_2=-x_2e_2`$. Let $`u=\langle e_1,Y^2e_1\rangle=\gamma y_1^2+(1-\gamma)y_2^2`$ and $`v=\langle e_2,Y^2e_2\rangle=s-u`$, $`s=y_1^2+y_2^2`$.
Then $`\mathrm{Tr}(Y|X|Y)=\mathrm{Tr}(|X|Y^2)=x_1u+x_2v=:A`$. The self-adjoint $`YXY`$ has determinant
$`\det X(\det Y)^2=-x_1x_2y_1^2y_2^2\le0`$, so its eigenvalues are $`\mu_+\ge0\ge\mu_-`$ and
```math
\mathrm{Tr}|YXY|=\mu_+-\mu_-=\sqrt{(\mu_++\mu_-)^2-4\mu_+\mu_-}=\sqrt{L^2+K^2},\qquad L=\mathrm{Tr}(XY^2)=x_1u-x_2v,\quad K=2\sqrt{x_1x_2}\,y_1y_2 .
```
Let $`\alpha=(x_1-x_2)/(x_1+x_2)`$ and $`\beta=2\sqrt{x_1x_2}/(x_1+x_2)>0`$, so $`\alpha^2+\beta^2=1`$. Since
$`A=(x_1-x_2)u+x_2s`$ and $`L=(x_1+x_2)u-x_2s`$, we have $`A-\alpha L=x_2s(1+\alpha)=2x_1x_2s/(x_1+x_2)`$ and
$`\beta K=4x_1x_2y_1y_2/(x_1+x_2)`$, hence $`A-\alpha L-\beta K=\frac{2x_1x_2}{x_1+x_2}(y_1-y_2)^2=:m`$. The
Cauchy–Schwarz inequality in $`\mathbb R^2`$ gives $`\sqrt{L^2+K^2}\ge\alpha L+\beta K`$, so
```math
z=A-\sqrt{L^2+K^2}\ \le\ A-\alpha L-\beta K=m,
```
with equality iff $`(L,K)`$ is a non-negative multiple of $`(\alpha,\beta)`$, i.e. (as $`\beta>0`$, $`K\ge0`$) iff
$`\beta L=\alpha K`$, i.e. $`(x_1+x_2)u-x_2s=(x_1-x_2)y_1y_2`$. This determines $`u`$, and
$`u-y_2^2=(y_1-y_2)(x_1y_2+x_2y_1)/(x_1+x_2)`$; for $`y_1>y_2`$, $`u=\gamma y_1^2+(1-\gamma)y_2^2`$ gives $`\gamma=\gamma_3\in[0,\frac12]`$.

*Real maximisers.* Take $`e_i`$ the standard basis, $`f_1,f_2`$ the rotated real basis with $`\cos^2t=\gamma_1`$ (resp. $`\gamma_3`$),
and, in (a), $`U=VW^{\mathrm T}`$ from a real singular value decomposition $`YX=W\Sigma V^{\mathrm T}`$ (diagonalise
the real symmetric $`(YX)^{\mathrm T}YX`$), so that $`\mathrm{Tr}(UYX)=\mathrm{Tr}\,\Sigma=\|YX\|_1`$. $`\square`$

<a id="r-11-16"></a>

**Corollary 11.16 (envelopes of Exp. 1 and Exp. 3 for blocks of size $`2`$).** Let all $`k_c=2`$ (so $`D=2C`$),
$`\delta\in[0,1)`$, $`\delta'\in[0,1]`$, $`\ell=\ell(\delta)`$, $`\ell'=\ell(\delta')`$ (with $`\ell(1)=0`$), and let $`Y`$
range over the positive operators with $`\lambda_{\max}(Y)=1`$ and $`\Delta(Y)=\delta'`$, i.e. $`\lambda_{\min}(Y)=\ell'`$.

(a) Exp. 1, over $`X\in\mathcal P_\delta`$, such $`Y`$ and unitaries $`U`$:
$`\ \sup z=C\,\frac{\delta^2\delta'^2}{(1+\delta)(1+\delta')}`$ and $`\inf z=-(2C-1+\ell\ell')`$.

(b) Exp. 3, over $`X=X^*`$ with $`|X|\in\mathcal P_\delta`$ and such $`Y`$:
$`\ \sup z=(C-\delta)\Big(\frac{2\delta'}{1+\delta'}\Big)^2`$ and $`\inf z=0`$.

All four values are attained. Consequently, with the other contrast free in $`[0,d]`$, $`d\in(0,1)`$: the upper
envelope of Exp. 1 is $`C\delta^2d^2/((1+\delta)(1+d))`$ along $`\Delta(X)=\delta`$ and (by the symmetry of (a))
$`C\delta'^2d^2/((1+\delta')(1+d))`$ along $`\Delta(Y)=\delta'`$; that of Exp. 3 is $`(C-\delta)\big(2d/(1+d)\big)^2`$ along $`\Delta(|X|)=\delta`$ and
$`C\big(2\delta'/(1+\delta')\big)^2`$ along $`\Delta(Y)=\delta'`$, whatever $`d`$. As $`d\to1`$ the first two increase to
$`C\delta^2/(2(1+\delta))`$ and $`C-\delta`$, the values at $`\Delta(Y)=1`$ (singular $`Y`$).

*Proof.* Let $`x_{c,1}\ge x_{c,2}`$ be the eigenvalues of $`X_c`$ (in (b): of $`|X|_c`$) and $`y_{c,1}\ge y_{c,2}`$ those of
$`Y_c`$; they lie in $`[\ell,1]`$ and $`[\ell',1]`$, and $`\ell`$ is the smaller one in its block. Note
$`1-\ell=2\delta/(1+\delta)`$, $`1+\ell=2/(1+\delta)`$, so $`(1-\ell)^2/(1+\ell)=2\delta^2/(1+\delta)`$ and
$`2\ell/(1+\ell)=1-\delta`$; likewise for $`\ell'`$.

(a) *Supremum.* By the triangle inequality and Theorem 1.9 in each block,
$`|\mathrm{Tr}(XUY)|\le\sum_c|\mathrm{Tr}(X_cU_cY_c)|\le\sum_c\|Y_cX_c\|_1`$, so by Proposition 11.15(a) (a block with
$`Y_c=0`$ contributes $`-\mathrm{Tr}(X_cY_c)=0`$)
```math
z\le\sum_c\frac{x_{c,1}y_{c,1}}4\,g\Big(\frac{x_{c,2}}{x_{c,1}}\Big)\,g\Big(\frac{y_{c,2}}{y_{c,1}}\Big),\qquad g(r)=\frac{(1-r)^2}{1+r}.
```
$`g`$ decreases on $`[0,1]`$, $`x_{c,2}/x_{c,1}\ge x_{c,2}\ge\ell`$, $`y_{c,2}/y_{c,1}\ge\ell'`$ and $`x_{c,1},y_{c,1}\le1`$, so
every term is at most $`g(\ell)g(\ell')/4=\delta^2\delta'^2/((1+\delta)(1+\delta'))`$. Equality: in every block
$`X_c=\mathrm{diag}(1,\ell)`$, $`Y_c`$ with eigenvalues $`1,\ell'`$ at $`\gamma=\gamma_1`$, and $`U_c=W_c^*`$, $`Y_cX_c=W_c|Y_cX_c|`$;
then $`\mathrm{Tr}(X_cU_cY_c)=\|Y_cX_c\|_1\ge0`$ in every block, so the moduli add.
*Infimum.* $`z\ge-\mathrm{Tr}(XY)`$, with equality for $`U=\bigoplus_cU_c`$, $`\mathrm{Tr}(Y_cX_cU_c)=0`$ (Lemma 11.13). By
Corollary 11.5(a), $`\mathrm{Tr}(XY)\le\sum_c(x_{c,1}y_{c,1}+x_{c,2}y_{c,2})`$. If $`\ell`$ and $`\ell'`$ lie in the same
block, it contributes at most $`1+\ell\ell'`$ and every other block at most $`2`$; otherwise the two blocks contribute at
most $`1+\ell`$ and $`1+\ell'`$, and $`2C-2+\ell+\ell'=2C-1+\ell\ell'-(1-\ell)(1-\ell')`$. So
$`\mathrm{Tr}(XY)\le2C-1+\ell\ell'`$, with equality for $`X=\mathrm{diag}(1,\ell)\oplus1\oplus\dots\oplus1`$ and
$`Y=\mathrm{diag}(1,\ell')\oplus1\oplus\dots\oplus1`$.

(b) Traces, trace norms and $`|\cdot|`$ act blockwise (Exercise 3), so $`z=\sum_cz_c`$, and by Proposition 11.15(b)
$`z_c\le H_c(y_{c,1}-y_{c,2})^2\le H_c(1-\ell')^2`$ with $`H_c=2x_{c,1}x_{c,2}/(x_{c,1}+x_{c,2})\le x_{c,1}\le1`$. In a block
$`c_-`$ in which $`|X|`$ has the eigenvalue $`\ell`$, $`x_{c_-,2}=\ell`$ and $`H_{c_-}=2\ell/(1+\ell/x_{c_-,1})\le2\ell/(1+\ell)=1-\delta`$.
Hence $`z\le(C-1+1-\delta)(1-\ell')^2`$. Equality: $`X_1=\mathrm{diag}(1,-\ell)`$, $`X_c=\mathrm{diag}(1,-1)`$ for $`c\ge2`$,
and every $`Y_c`$ with eigenvalues $`1,\ell'`$ at $`\gamma=\gamma_3`$ of its block (for $`c\ge2`$, $`\gamma_3=\frac12`$). The infimum is
Proposition 11.9, with equality for commuting $`X`$, $`Y`$.

*Free contrast.* The supremum in (a) increases in $`\delta'`$ (strictly if $`\delta>0`$), and the one in (b) increases
strictly in $`\delta'`$ and decreases in $`\delta`$. Over the free contrast in $`[0,d]`$ the suprema are therefore
attained at $`d`$ (along $`\Delta(Y)`$ in Exp. 3: at $`\Delta(|X|)=0`$), and the limits as $`d\to1`$ are reached only at
$`\Delta(Y)=1`$ (in Exp. 1 if $`\delta>0`$). $`\square`$

For $`C=1`$ and $`d=0.999`$ these are the closed forms that the companion paper first observed numerically (in real
algebras). For every $`C`$ the upper envelopes are the values of direct sums of $`C`$ solutions on $`M_2`$
($`C-\delta=(1-\delta)+(C-1)\cdot1`$ in Exp. 3): on $`M_2\oplus M_2`$ the direct-sum inequality above is an equality. Section 11.7 extends these suprema to blocks of
any size (Theorem 11.20).

*Computation:* the maximisers of Proposition 11.15, the envelopes of Corollary 11.16 for $`C=1,2`$ and their limit
$`d\to1`$, and a multistart search against them — [tutorial, Section 11.5](../tutorial/11_optimisation.md#s-11-5).

<a id="s-11-6"></a>

## 11.6 What to compare with

*In the library:* comparison with other software, among them SDP solvers, whose convex relaxations bound a supremum from above — [tutorial, Section 11.6](../tutorial/11_optimisation.md#s-11-6).

<a id="s-11-7"></a>

## 11.7 Envelopes for blocks of any size

The suprema of Corollary 11.16 extend to blocks of every size. Let
```math
N=\sum_c\Big\lfloor\frac{k_c}2\Big\rfloor,
```
the largest number of mutually orthogonal two-dimensional subspaces that fit into the sectors. The upper envelopes are
attained by direct sums of the maximisers of Proposition 11.15 on $`N`$ blocks of size $`2`$ and of $`1\times1`$ blocks.
The proof reduces the spectra to two values, then the operators to blocks of size $`1`$ and $`2`$ (Lemma 11.17), where
Proposition 11.15 applies — except in Exp. 3 when the eigenvector of $`|X|`$ for $`\lambda_{\min}`$ is spread over several
blocks; Proposition 11.19 settles this case by a linearisation.

*Vertex argument.* A convex function $`f`$ on a face $`F=\{x\in[s,t]^n:x_i=t\ (i\in I),\ x_j=s\ (j\in J)\}`$ of a cube
attains its maximum over $`F`$ at a vertex of $`F`$, a point with all coordinates in $`\{s,t\}`$: every point of $`F`$ is
a convex combination of these vertices (write one free coordinate as a convex combination of $`s`$ and $`t`$, and
induct on the number of free coordinates), and the value of $`f`$ there is at most the largest value at the vertices.

<a id="r-11-17"></a>

**Lemma 11.17 (two projections).** For projections $`P,Q\in M_n(\mathbb C)`$, $`\mathbb C^n`$ is an orthogonal sum of
subspaces of dimension $`1`$ or $`2`$, each invariant under $`P`$ and $`Q`$, such that each two-dimensional one has an
orthonormal basis $`f,f'`$ with $`P=|f\rangle\langle f|`$ and $`Q=|q\rangle\langle q|`$ on it,
$`q=\cos\theta\,f+e^{i\psi}\sin\theta\,f'`$, $`\theta\in[0,\pi/2]`$. For real $`P,Q`$ everything can be chosen real.

*Proof.* $`PQP`$ is self-adjoint on $`\mathrm{ran}\,P`$; let $`v_1,\dots,v_r`$ be an orthonormal eigenbasis,
$`PQPv_i=\mu_iv_i`$ (Theorem 1.3), and $`W_i=\mathrm{span}\{v_i,Qv_i\}`$. Then $`PQv_i=PQPv_i=\mu_iv_i`$, so $`W_i`$ is
invariant under $`P`$ and (as $`Q^2=Q`$) under $`Q`$. For $`i\ne j`$,
$`\langle v_i,Qv_j\rangle=\langle Qv_i,Qv_j\rangle=\langle v_i,PQPv_j\rangle=0`$, so the $`W_i`$ are mutually orthogonal.
Their sum contains $`\mathrm{ran}\,P`$; its orthogonal complement $`W`$ is invariant under $`Q`$ (both projections are
self-adjoint) and contained in $`\ker P`$, and an orthonormal eigenbasis of $`Q`$ on $`W`$ gives one-dimensional summands.
A $`W_i`$ of dimension $`1`$ is spanned by $`v_i`$. If $`\dim W_i=2`$, put $`f=v_i`$ and let $`f'`$ be the normalised
component of $`Qv_i`$ orthogonal to $`f`$; then $`P=|f\rangle\langle f|`$ on $`W_i`$, and $`Q`$ maps $`W_i`$ onto
$`\mathrm{span}\{Qv_i\}`$, so $`Q=|q\rangle\langle q|`$ there for a unit vector $`q=\alpha f+\beta f'`$; multiplying $`q`$ by a
phase makes $`\alpha\ge0`$, and $`\alpha=\cos\theta`$, $`\beta=e^{i\psi}\sin\theta`$ with $`\theta\in[0,\pi/2]`$. In the real
case the eigenbasis, $`f'`$ and $`q`$ are real and $`e^{i\psi}=\pm1`$. $`\square`$

<a id="r-11-18"></a>

**Lemma 11.18 (Exp. 3: monotonicity and convexity).** Let $`X=X^*`$ be invertible, $`P`$ the projection onto its
positive part ($`X_+=PX`$), $`S=2P-1`$, $`Y\ge0`$, and $`z(X,Y)=\mathrm{Tr}(Y|X|Y)-\mathrm{Tr}|YXY|`$.

(a) If $`|X|\le|X'|`$, where $`|X'|`$ commutes with $`P`$ and $`X'=S|X'|`$, then $`z(X,Y)\le z(X',Y)`$.

(b) For fixed $`X`$, $`z`$ is a convex function of $`K=Y^2`$: for $`K_1,K_2\ge0`$ and $`t\in[0,1]`$,
$`z\big(X,(tK_1+(1-t)K_2)^{1/2}\big)\le t\,z(X,K_1^{1/2})+(1-t)\,z(X,K_2^{1/2})`$.

*Proof.* With $`A=YX_+Y`$ and $`B=YX_-Y`$, $`z=\mathrm{Tr}\,A+\mathrm{Tr}\,B-\|A-B\|_1`$ (proof of Proposition 11.9).

(a) $`X_+=P|X|=|X|^{1/2}P|X|^{1/2}`$ and $`X_-=(1-P)|X|`$, and $`P|X'|-P|X|=P(|X'|-|X|)P\ge0`$, likewise with $`1-P`$; so
$`A\le A'`$ and $`B\le B'`$ for the operators $`A',B'`$ of $`X'`$ (Corollary 1.7). For $`R\ge0`$,
$`\|A+R-B\|_1\le\|A-B\|_1+\mathrm{Tr}\,R`$ (triangle inequality, $`\|R\|_1=\mathrm{Tr}\,R`$), so
$`\mathrm{Tr}\,A+\mathrm{Tr}\,B-\|A-B\|_1`$ does not decrease when $`A`$ or $`B`$ increases.

(b) For self-adjoint $`H`$ and self-adjoint $`T`$ with $`-1\le T\le1`$, $`|\mathrm{Tr}(TH)|\le\|H\|_1`$: the operators
$`T\pm i(1-T^2)^{1/2}`$ are unitary (Theorem 1.5) with mean $`T`$, and Theorem 1.9 applies to each. Equality holds for
$`T=2Q-1`$, $`Q`$ the projection onto the positive part of $`H`$ (proof of Lemma 11.8). With $`T=1-2E`$,
```math
z=2\min_{0\le E\le1}\mathrm{Tr}\big(EA+(1-E)B\big).
```
Let first $`Y`$ be invertible. Then $`N=|X|^{1/2}Y`$ is invertible, $`A=N^*PN`$, $`B=N^*(1-P)N`$, and $`E\mapsto F=NEN^*`$
maps $`\{0\le E\le1\}`$ onto $`\{0\le F\le G\}`$, $`G=NN^*=|X|^{1/2}Y^2|X|^{1/2}`$ (Corollary 1.7, both ways). As
$`\mathrm{Tr}(EA)=\mathrm{Tr}(PF)`$ and $`\mathrm{Tr}((1-E)B)=\mathrm{Tr}((1-P)(G-F))`$,
$`z=2\min\{\mathrm{Tr}(PF+(1-P)(G-F)):0\le F\le G\}`$. If $`0\le F_i\le G_i`$, then $`tF_1+(1-t)F_2`$ is admissible for
$`tG_1+(1-t)G_2`$, and the objective is linear in $`(F,G)`$; so $`z`$ is convex in $`G`$, which is affine in $`Y^2`$. For
singular $`K_1`$ or $`K_2`$ apply this to $`K_i+\varepsilon`$ and let $`\varepsilon\to0`$ ($`z`$ is continuous in $`Y`$, and
$`K\mapsto K^{1/2}`$ is continuous). $`\square`$

<a id="r-11-19"></a>

**Proposition 11.19 (Exp. 3 on $`m`$ blocks).** Let $`X=X_+\oplus(-1)`$ on $`\mathbb C^m\oplus\mathbb C^m`$, where
$`0<X_+\le1`$ has eigenvalues $`x_1,\dots,x_m`$; let $`C=\mathrm{diag}(\cos\theta_j)`$, $`\Sigma=\mathrm{diag}(\sin\theta_j)`$, $`Q`$
the projection onto the range of the isometry $`v=\binom C\Sigma:\mathbb C^m\to\mathbb C^{2m}`$, $`b\in[0,1)`$ and
$`Y=b+(1-b)Q`$. Then
```math
\mathrm{Tr}(Y|X|Y)-\mathrm{Tr}|YXY|\ \le\ (1-b)^2\sum_{i=1}^m\frac{2x_i}{1+x_i},
```
with equality for $`X_+=\mathrm{diag}(x_1,\dots,x_m)`$ and $`\cos^2\theta_i=(1+bx_i)/((1+x_i)(1+b))`$.

*Proof.* Two facts from linear algebra: for square $`R,T`$, $`\det(\lambda-RT)=\det(\lambda-TR)`$ (for invertible $`R`$,
$`RT=R(TR)R^{-1}`$; in general by continuity in $`R`$), and for a block matrix with invertible lower right block $`E`$,
$`\det\begin{pmatrix}R&S\\T&E\end{pmatrix}=\det E\,\det(R-SE^{-1}T)`$ (multiply on the left by
$`\begin{pmatrix}1&-SE^{-1}\\0&1\end{pmatrix}`$).

Let $`\gamma=1-b^2`$, $`K=1+X_+`$ and $`A=C^2`$, so $`\Sigma^2=1-A`$; $`Q=vv^*`$ has the blocks $`C^2,\ C\Sigma,\ C\Sigma,\ \Sigma^2`$,
and $`Y^2=b^2+\gamma Q`$ (as $`Q^2=Q`$). For $`\lambda>0`$ the lower right block of $`\lambda-XY^2`$ is the invertible
diagonal matrix $`E=\lambda+b^2+\gamma\Sigma^2=\lambda+1-\gamma A`$, and since all of $`C,\Sigma,E`$ are diagonal and
$`(b^2+\gamma A)E-\gamma^2A\Sigma^2=\lambda(b^2+\gamma A)+b^2`$,
```math
\det(\lambda-XY^2)=\det E\,\det\Big(\lambda-X_+\big(\lambda(b^2+\gamma A)+b^2\big)E^{-1}\Big)
=\det\big((\lambda+1)(\lambda-b^2X_+)-\lambda\gamma KA\big)=p(\lambda),
```
where $`p(\lambda)=\det\big((\lambda+1)(\lambda-b^2X_+)-\lambda\gamma K^{1/2}AK^{1/2}\big)`$ (conjugate the matrix
$`K^{-1}(\lambda+1)(\lambda-b^2X_+)-\lambda\gamma A`$ by $`K^{1/2}`$). The self-adjoint matrix
```math
\mathcal L=\begin{pmatrix}\gamma K^{1/2}AK^{1/2}+b^2X_+-1&bX_+^{1/2}\\ bX_+^{1/2}&0\end{pmatrix}
```
has, by the same formula with lower right block $`\lambda`$,
$`\det(\lambda-\mathcal L)=\lambda^m\det(\lambda+1-\gamma K^{1/2}AK^{1/2}-b^2X_+-b^2X_+/\lambda)=p(\lambda)`$. Since
$`\det(\lambda-YXY)=\det(\lambda-XY^2)`$, the characteristic polynomials of $`YXY`$ and $`\mathcal L`$ agree for
$`\lambda>0`$, hence everywhere, and these self-adjoint matrices have the same eigenvalues: $`\mathrm{Tr}|YXY|=\|\mathcal L\|_1`$.
Let $`\Delta=(1-X_+)(1+X_+)^{-1}`$ and $`\Gamma=2X_+^{1/2}(1+X_+)^{-1}`$ (functions of $`X_+`$, so they commute, and
$`\Delta^2+\Gamma^2=1`$), and
```math
T=\begin{pmatrix}-\Delta&\Gamma\\ \Gamma&\Delta\end{pmatrix},\qquad T^2=1 .
```
$`T`$ is self-adjoint with eigenvalues $`\pm1`$, so $`\|\mathcal L\|_1\ge\mathrm{Tr}(T\mathcal L)`$ (proof of Lemma 11.18(b)).
Now $`\mathrm{Tr}(Y|X|Y)=\mathrm{Tr}(|X|Y^2)=b^2\mathrm{Tr}\,X_++\gamma\mathrm{Tr}(X_+A)+mb^2+\gamma\mathrm{Tr}(1-A)`$ and, as
$`K^{1/2}\Delta K^{1/2}=1-X_+`$,
$`\mathrm{Tr}(T\mathcal L)=-\gamma\mathrm{Tr}((1-X_+)A)-\mathrm{Tr}(\Delta(b^2X_+-1))+2b\,\mathrm{Tr}(\Gamma X_+^{1/2})`$. The terms
containing $`A`$ cancel, and in an eigenbasis of $`X_+`$
```math
\mathrm{Tr}(Y|X|Y)-\mathrm{Tr}|YXY|\le\sum_{i=1}^m\Big(b^2x_i+1+\frac{1-x_i}{1+x_i}(b^2x_i-1)-\frac{4bx_i}{1+x_i}\Big)
=(1-b)^2\sum_{i=1}^m\frac{2x_i}{1+x_i}.
```
For diagonal $`X_+`$, $`X`$ and $`Y`$ are direct sums of the $`2\times2`$ blocks on $`\mathrm{span}\{e_i,e_{m+i}\}`$, with
$`X=\mathrm{diag}(x_i,-1)`$ there, and Proposition 11.15(b) gives the value $`\frac{2x_i}{1+x_i}(1-b)^2`$ when
$`\gamma_3=\sin^2\theta_i`$, i.e. at the stated angle. $`\square`$

Since $`\mathcal L`$ is affine in $`A=PQP`$ ($`P`$ the projection onto the first summand), $`z`$ is a concave function of
$`A`$, and $`T`$ is the sign of $`\mathcal L`$ at the maximiser. For fixed angles, however, the maximum over the
eigenvector of $`X_+`$ for its smallest eigenvalue is in general not attained inside a single block (a numerical
observation); only the joint maximum is.

<a id="r-11-20"></a>

**Theorem 11.20 (upper envelopes of Exp. 1 and Exp. 3).** Let $`N\ge1`$, $`\delta\in[0,1)`$, $`\delta'\in[0,1]`$,
$`\ell=\ell(\delta)`$, $`\ell'=\ell(\delta')`$ (with $`\ell(1)=0`$), and let $`Y`$ range over the positive operators with
$`\lambda_{\max}(Y)=1`$ and $`\Delta(Y)=\delta'`$.

(a) Exp. 1, over $`X\in\mathcal P_\delta`$, such $`Y`$ and unitaries $`U\in\mathcal M`$:
$`\ \sup z=N\,\frac{\delta^2\delta'^2}{(1+\delta)(1+\delta')}`$.

(b) Exp. 3, over $`X=X^*`$ with $`|X|\in\mathcal P_\delta`$ and such $`Y`$:
$`\ \sup z=(N-\delta)\Big(\frac{2\delta'}{1+\delta'}\Big)^2`$ if every $`k_c`$ is even, and
$`\ \sup z=N\Big(\frac{2\delta'}{1+\delta'}\Big)^2`$ if some $`k_c`$ is odd.

Both suprema are attained, by real operators, and are the same in $`\bigoplus_cM_{k_c}(\mathbb R)`$.

*Proof.* Upper bounds in $`\mathcal M`$ hold a fortiori for real operators, and the maximisers below are real.

(a) By Theorem 1.9 in each block, $`|\mathrm{Tr}(XUY)|\le\|YX\|_1`$, with equality for $`U=W^*`$, $`YX=W|YX|`$, a unitary of
$`\mathcal M`$; so $`\sup_Uz=F(X,Y):=\|YX\|_1-\mathrm{Tr}(XY)`$. For a fixed orthonormal eigenbasis $`(e_s)`$ of $`X`$ (adapted to
the sectors), $`F`$ is a convex function of the eigenvalues $`x\in[\ell,1]^D`$ (a norm of a linear function of $`x`$ minus a
linear function), and $`X\in\mathcal P_\delta`$ means that $`x`$ lies on a face $`\{x_s=1,\ x_{s'}=\ell\}`$; by the vertex
argument $`F(X,Y)\le F(X',Y)`$ with $`X'=\ell+(1-\ell)P`$ for a projection $`P\in\mathcal M`$. Likewise
$`F(X',Y)\le F(X',Y')`$ with $`Y'=\ell'+(1-\ell')Q`$. Apply Lemma 11.17 to $`P_c,Q_c`$ in every sector: $`X'`$ and $`Y'`$ are
direct sums over subspaces of dimension $`1`$ and $`2`$, and so are $`\|Y'X'\|_1`$ and $`\mathrm{Tr}(X'Y')`$. A one-dimensional
summand contributes $`xy-xy=0`$. On a two-dimensional one, $`X'`$ has eigenvalues $`1,\ell`$ and $`Y'`$ eigenvalues $`1,\ell'`$,
so by Proposition 11.15(a) it contributes at most $`(1-\ell)^2(1-\ell')^2/(4(1+\ell)(1+\ell'))=\delta^2\delta'^2/((1+\delta)(1+\delta'))`$
(see the proof of Corollary 11.16). A sector of size $`k_c`$ contains at most $`\lfloor k_c/2\rfloor`$ of them. Equality: $`N`$
blocks with $`X=\mathrm{diag}(1,\ell)`$ and $`Y`$ with eigenvalues $`1,\ell'`$ at $`\gamma=\gamma_1`$, $`X=Y=1`$ on the remaining
one-dimensional blocks, and the unitary $`W^*`$.

(b) $`X`$ is invertible ($`\ell>0`$), and $`z(-X,Y)=z(X,Y)`$. Let $`P`$ be the projection onto the positive part of $`X`$
and $`S=2P-1`$. $`|X|`$ commutes with $`P`$ (both are functions of $`X`$), so it has a unit eigenvector $`e`$ for the
eigenvalue $`\ell`$ in one sector and in $`\mathrm{ran}\,P`$ or $`\ker P`$; replacing $`X`$ by $`-X`$ if necessary,
$`e\in\mathrm{ran}\,P`$. Then $`|X|\le|X'|:=1-(1-\ell)|e\rangle\langle e|`$, $`X'=S|X'|=S-(1-\ell)|e\rangle\langle e|`$ has
$`|X'|\in\mathcal P_\delta`$ ($`D\ge2`$), and $`z(X,Y)\le z(X',Y)`$ by Lemma 11.18(a). By Lemma 11.18(b) and the vertex argument
for the eigenvalues of $`Y^2`$ in $`[\ell'^2,1]`$, $`z(X',Y)\le z(X',Y')`$ with $`Y'=\ell'+(1-\ell')Q`$ for a projection
$`Q\in\mathcal M`$. Apply Lemma 11.17 to $`P_c,Q_c`$ in every sector.
In a sector without $`e`$, $`X'=S`$ is a direct sum over the subspaces; a one-dimensional summand contributes $`0`$ (there
$`X'`$ and $`Y'`$ commute, Proposition 11.9), and a two-dimensional one at most $`(1-\ell')^2`$ by Proposition 11.15(b)
with $`x_1=x_2=1`$: in all at most $`\lfloor k_c/2\rfloor(1-\ell')^2`$.
In the sector of $`e`$, Lemma 11.18(a) with $`|X'|\le1`$ gives the same bound, and so $`\lfloor k_c/2\rfloor(1-\ell')^2`$ if
$`k_c`$ is odd, and $`(k_c/2-1)(1-\ell')^2`$ if the sector contains fewer than $`k_c/2`$ two-dimensional subspaces.
Otherwise its $`m=k_c/2`$ two-dimensional subspaces fill it, with bases $`f_j,f_j'`$, $`f_j\in\mathrm{ran}\,P`$,
$`f_j'\in\ker P`$. The unitary that multiplies $`f_j'`$ by $`e^{-i\psi_j}`$ commutes with $`X'`$ ($`X'=-1`$ on $`\ker P`$) and
removes the phases of Lemma 11.17, and in the basis $`(f_1,\dots,f_m,f_1',\dots,f_m')`$ the sector is the setting of
Proposition 11.19 with $`X_+=1-(1-\ell)|e\rangle\langle e|`$, whose eigenvalues are $`\ell,1,\dots,1`$; as
$`2\ell/(1+\ell)=1-\delta`$, it contributes at most $`(1-\ell')^2(1-\delta+m-1)=(k_c/2-\delta)(1-\ell')^2`$.
Summing over the sectors, $`z\le(N-\delta)(1-\ell')^2`$ if every $`k_c`$ is even (note $`N-1\le N-\delta`$) and
$`z\le N(1-\ell')^2`$ otherwise; $`1-\ell'=2\delta'/(1+\delta')`$.
Equality: if every $`k_c`$ is even, $`X=\mathrm{diag}(\ell,-1)`$ on one block and $`\mathrm{diag}(1,-1)`$ on the other
$`N-1`$ blocks, and $`Y`$ with eigenvalues $`1,\ell'`$ at $`\gamma=\gamma_3`$ of each block (Proposition 11.15(b)); if some $`k_c`$
is odd, $`X=\mathrm{diag}(1,-1)`$ on $`N`$ blocks, the eigenvalue $`\ell`$ on a one-dimensional block of an odd sector and
$`1`$ on the other one-dimensional blocks, $`Y`$ as before on the blocks and $`1`$ on the one-dimensional blocks. $`\square`$

For $`k_c=2`$ this is the supremum part of Corollary 11.16 ($`N=C`$). For $`C`$ sectors of size $`k`$ and the other contrast
free in $`[0,d]`$, $`d\in(0,1)`$, the suprema increase in $`\delta'`$ and the upper envelopes are: in Exp. 1
$`C\lfloor k/2\rfloor\delta^2d^2/((1+\delta)(1+d))`$, along either contrast; in Exp. 3 along $`\Delta(|X|)=\delta`$,
$`(Ck/2-\delta)\big(2d/(1+d)\big)^2`$ for even $`k`$ and $`C\lfloor k/2\rfloor\big(2d/(1+d)\big)^2`$ for odd $`k`$; along
$`\Delta(Y)=\delta'`$, attained at $`|X|=1`$, $`C\lfloor k/2\rfloor\big(2\delta'/(1+\delta')\big)^2`$. For even $`k`$ these are the
values of the direct sums of Section 11.5: the direct-sum inequality is an equality, in both functionals. For odd $`k`$ the
direct sums with a one-dimensional summand are optimal in Exp. 1, but not in Exp. 3: the eigenvalue $`\ell`$ of $`|X|`$ sits on
a one-dimensional block, where it costs nothing, and the envelope does not depend on $`\delta`$.

*Computation:* Proposition 11.19 (spectra of $`YXY`$ and $`\mathcal L`$, the bound and its equality case), the maximisers of
Theorem 11.20 for blocks of several sizes, and multistart searches that reach the envelopes from structured starting
points — [tutorial, Section 11.7](../tutorial/11_optimisation.md#s-11-7).

<a id="x-11"></a>

## Exercises

1. (**Cayley transform.**) For $`S^*=-S`$ show that $`C(S)=(1-S)(1+S)^{-1}`$ is unitary, and that $`S\mapsto C(S)`$ is a
   bijection onto the unitaries without eigenvalue $`-1`$. So, unlike the exponential map (Lemma 11.1), a Cayley chart
   $`U_0C(S)`$ misses $`-U_0`$.
2. (**Real Haar measure.**) For $`O`$ Haar-distributed on $`O(k)`$ (QR of a real Ginibre matrix) and real $`A`$, show
   $`\mathbb E\mathrm{Tr}(AO)^2=\|A\|_2^2/k`$, and check it with a real algebra.
3. (**Direct sums.**) Show that $`z`$ of Exp. 3 is additive over direct sums and that
   $`\Delta(A\oplus A')=\max(\Delta(A),\Delta(A'))`$ for $`A,A'\ge0`$ with $`\lambda_{\max}=1`$. Deduce the direct-sum
   inequality of Section 11.5 (with $`\Delta(Y)`$ free in a fixed interval $`[0,b]`$).

<details>
<summary>Solutions</summary>

1. $`S`$ is normal and $`iS`$ is self-adjoint, so $`S`$ has imaginary eigenvalues (Corollary 1.4(2)) and $`1\pm S`$ are
   invertible; they commute, and $`(1+S)^*=1-S`$. Hence
   $`C^*C=(1-S)^{-1}(1+S)(1-S)(1+S)^{-1}=1`$, and $`1+C=2(1+S)^{-1}`$ is invertible, so $`-1\notin\mathrm{spec}\,C`$.
   Conversely, for a unitary $`U`$ with $`1+U`$ invertible, $`S=(1-U)(1+U)^{-1}`$ is skew-adjoint:
   $`S^*=(1+U^*)^{-1}(1-U^*)=(U+1)^{-1}U\,U^*(U-1)=(U+1)^{-1}(U-1)=-S`$, as $`U\pm1`$ commute. For commuting $`C`$ and $`S`$,
   $`C(1+S)=1-S`$ is equivalent to $`S(1+C)=1-C`$ (both say $`C+CS+S=1`$). Now $`C(S)`$ commutes with $`S`$, and
   $`S=(1-U)(1+U)^{-1}`$ commutes with $`U`$; so $`U=C(S)`$ iff $`S=(1-U)(1+U)^{-1}`$: the two maps are inverse to each other.

2. By Theorem B.17 (also Lemma C.2) the sign-corrected QR factor $`O`$ of a real Gaussian matrix, which the library
   samples for a real algebra, is Haar distributed on $`O(k)`$, so $`VO`$ has the law of $`O`$ for $`V\in O(k)`$.
   Left multiplication by $`\mathrm{diag}(\pm1)`$ gives $`\mathbb E[O_{ij}O_{lm}]=0`$ for $`i\ne l`$, permutations show
   that $`\mathbb E[O_{ij}O_{im}]`$ does not depend on $`i`$, and $`\sum_iO_{ij}O_{im}=\delta_{jm}`$; so
   $`\mathbb E[O_{ij}O_{lm}]=\delta_{il}\delta_{jm}/k`$ and $`\mathbb E\mathrm{Tr}(AO)^2=\sum_{ij}A_{ji}^2/k`$.

3. Traces and trace norms of block-diagonal operators are sums over blocks, and $`|X\oplus X'|=|X|\oplus|X'|`$;
   $`\lambda_{\min}(A\oplus A')=\min(\lambda_{\min}(A),\lambda_{\min}(A'))`$, and $`\Delta`$ decreases with $`\lambda_{\min}`$.
   So if $`(X,Y)`$ is admissible at $`\delta`$ and $`(X',Y')`$ at $`\delta'\le\delta`$, then $`(X\oplus X',Y\oplus Y')`$ is
   admissible at $`\delta`$, with $`z=z(X,Y)+z(X',Y')`$; take suprema over both pairs and over $`\delta'\le\delta`$.

</details>

*Computation:* Exercises 1 and 2 numerically — [tutorial, Exercises](../tutorial/11_optimisation.md#x-11).

**Next:** [Chapter 12](12_equivariant_learning.md) treats equivariant machine learning: layers that commute with a symmetry are
elements of a commutant (Chapter 2), built from the same batched, differentiable operators.
