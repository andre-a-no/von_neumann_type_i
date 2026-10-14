**English** · [Español](es/03_traces_and_states.md) · [Français](fr/03_traces_and_states.md) · [Deutsch](de/03_traces_and_states.md) · [中文](zh/03_traces_and_states.md) · [日本語](ja/03_traces_and_states.md) · [Русский](ru/03_traces_and_states.md)

# 3. Traces, states and entropy

On $`M_n(\mathbb C)`$ there is essentially one trace. On $`\bigoplus_cM_{k_c}`$ there is a
$`C`$-parameter family, and the choice matters: densities, entropies and the duals of maps
all depend on it. This chapter classifies the traces, represents states by densities,
proves the basic entropy inequalities, and ends with a plausible trace inequality that
is false.

**You will learn** the classification of the traces on $`\bigoplus_cM_{k_c}`$ and how the densities
of a state with respect to different traces are related, the classical part of a state, the
decomposition of the von Neumann entropy over sectors, Klein's inequality and its consequences,
and the Michelson contrast. The numerical checks of this chapter are in
[Tutorial 3](../tutorial/03_traces_and_states.md).

Throughout, $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$ with all $`k_c\ge1`$,
$`D=\sum_ck_c`$, $`P_c`$ the minimal central projections, and $`A=\bigoplus_cA_c`$.

<a id="s-3-1"></a>

## 3.1 Traces

A **trace** on $`\mathcal M`$ is a linear functional $`\tau`$ with $`\tau(AB)=\tau(BA)`$ for all
$`A,B`$. It is **positive** if $`\tau(A^*A)\ge0`$ for all $`A`$, and **faithful** if moreover
$`\tau(A^*A)=0`$ implies $`A=0`$.

<a id="r-3-1"></a>

**Theorem 3.1.** Every trace on $`\mathcal M`$ has the form
$`\mathrm{tr}_w(A)=\sum_cw_c\mathrm{Tr}\,A_c`$ with $`w\in\mathbb C^C`$. It is positive
iff all $`w_c\ge0`$ and faithful iff all $`w_c>0`$.

*Proof.* Let $`E^c_{ij}`$ be the matrix units of the $`c`$-th block. For a trace $`\tau`$,
$`\tau(E^c_{ij})=\tau(E^c_{i1}E^c_{1j})=\tau(E^c_{1j}E^c_{i1})=\delta_{ij}\tau(E^c_{11})`$.
Put $`w_c=\tau(E^c_{11})`$; by linearity $`\tau(A)=\sum_c\sum_i(A_c)_{ii}w_c`$. Conversely
every $`\mathrm{tr}_w`$ is a trace. Positivity: $`E^c_{11}=(E^c_{11})^*E^c_{11}`$, so
$`w_c=\tau\big((E^c_{11})^*E^c_{11}\big)\ge0`$ is necessary, and $`w_c>0`$ is necessary for faithfulness, as
$`E^c_{11}\neq0`$. Conversely $`\mathrm{tr}_w(A^*A)=\sum_cw_c\|A_c\|_2^2`$, which is $`\ge0`$ if all $`w_c\ge0`$, and
vanishes only for $`A=0`$ if all $`w_c>0`$. $`\square`$

The library implements three faithful traces:

| name | weights $`w_c`$ | value on $`1`$ | meaning |
|---|---|---|---|
| `Tr_blunt` | $`1`$ | $`D`$ | trace of $`\bigoplus_c\mathbb C^{k_c}`$ (multiplicities $`m_c=1`$) |
| `Tr_norm` | $`1/k_c`$ | $`C`$ | every block normalised |
| `tau_vN` | $`1/(Ck_c)`$ | $`1`$ | the tracial state with equal sector weights $`1/C`$ |

By Section 2.3, the trace of a space on which $`\mathcal M`$ acts with multiplicities $`m_c`$
is $`\mathrm{tr}_w`$ with $`w_c=m_c`$. A faithful trace is thus a choice of "how much each
sector counts".

*In the library:* the traces `Tr_blunt`, `Tr_norm`, `tau_vN`, their values on $`1`$ and the trace property of Theorem 3.1 — [tutorial, Section 3.1](../tutorial/03_traces_and_states.md#s-3-1).

<a id="s-3-2"></a>

## 3.2 States and densities

A **state** is a linear functional $`\omega`$ on $`\mathcal M`$ that is positive
($`\omega(A^*A)\ge0`$) and normalised ($`\omega(1)=1`$).

<a id="r-3-2"></a>

**Theorem 3.2.** Let $`\mathrm{tr}_w`$ be a faithful trace. For every linear functional
$`\omega`$ there is a unique $`\rho_w\in\mathcal M`$ with $`\omega(A)=\mathrm{tr}_w(\rho_wA)`$ for all
$`A`$. The functional $`\omega`$ is positive iff $`\rho_w\ge0`$, and a state iff moreover
$`\mathrm{tr}_w(\rho_w)=1`$. For two faithful traces, the densities are related blockwise by
```math
w_c\,(\rho_w)_c=w'_c\,(\rho_{w'})_c .
```

*Proof.* $`\langle X,Y\rangle_w=\mathrm{tr}_w(X^*Y)`$ is an inner product on $`\mathcal M`$
(Theorem 3.1). By the Riesz representation theorem in the finite-dimensional space
$`\mathcal M`$, there is a unique $`R`$ with $`\omega(A)=\langle R,A\rangle_w=\mathrm{tr}_w(R^*A)`$; set
$`\rho_w=R^*`$. If $`\rho_w\ge0`$, then
$`\omega(A^*A)=\sum_cw_c\mathrm{Tr}(A_c\,(\rho_w)_c\,A_c^*)\ge0`$, since
$`A_c(\rho_w)_cA_c^*\ge0`$ (Corollary 1.7 with $`T=(\rho_w)_c^{1/2}A_c^*`$), a positive matrix has trace
$`\ge0`$ (its eigenvalues are $`\ge0`$, Corollary 1.4(3)) and $`w_c>0`$. Conversely, if
$`\omega`$ is positive, take $`x\in\mathbb C^{k_c}`$ and $`A`$ the rank-one operator
$`|x\rangle\langle x|`$ in block $`c`$ ($`x\neq0`$); then $`A=A^*A/\|x\|^2`$ and
$`0\le\omega(A)=w_c\langle x,(\rho_w)_cx\rangle`$, so $`(\rho_w)_c\ge0`$. The normalisation is
$`\omega(1)=\mathrm{tr}_w(\rho_w)`$. Finally both $`\sum_cw_c\mathrm{Tr}((\rho_w)_cA_c)`$ and
$`\sum_cw'_c\mathrm{Tr}((\rho_{w'})_cA_c)`$ equal $`\omega(A)`$ for all $`A`$, so $`\bigoplus_cw_c(\rho_w)_c`$ and
$`\bigoplus_cw'_c(\rho_{w'})_c`$ are both densities of $`\omega`$ for the blunt trace ($`w_c=1`$), and by uniqueness
$`w_c(\rho_w)_c=w'_c(\rho_{w'})_c`$. $`\square`$

*In the library:* a state is a `DensityMatrix`, stored by its blunt density; the conversions of Theorem 3.2 are `DensityMatrix.density` and `DensityMatrix.from_density` — [tutorial, Section 3.2](../tutorial/03_traces_and_states.md#s-3-2).

> **Physical meaning** (physical interpretation, not a theorem of this course). $`\omega(A)`$ is the expectation value of the observable $`A`$.
> The density $`\rho`$ of the blunt trace is the usual density matrix of quantum
> mechanics on $`\bigoplus_c\mathbb C^{k_c}`$. Densities with respect to other traces appear
> when the same algebra acts with multiplicities (Chapter 2), for example on the
> symmetric and the non-symmetric part of a spin system.

<a id="s-3-3"></a>

## 3.3 The classical part of a state

The restriction of $`\omega`$ to the centre is a probability distribution on the sectors,
$`p_c=\omega(P_c)=\mathrm{Tr}\,\rho_c`$. If $`p_c>0`$, the **conditional state** in sector $`c`$
has density $`\sigma_c=\rho_c/p_c`$ on $`\mathbb C^{k_c}`$, and
```math
\rho=\bigoplus_cp_c\,\sigma_c
```
(if $`p_c=0`$ then $`\rho_c=0`$, since the eigenvalues of $`\rho_c\ge0`$ are $`\ge0`$ with sum
$`\mathrm{Tr}\,\rho_c=0`$, Corollary 1.4(3), hence all $`0`$, and $`\rho_c`$ is unitarily diagonalisable,
Theorem 1.3; the term is omitted).
So every state on $`\mathcal M`$ is a convex combination of states of the sectors. The map
$`\rho\mapsto\sum_cP_c\rho P_c`$ is the identity on $`\mathcal M`$, because every element of
$`\mathcal M`$ is block diagonal.

> **Physical meaning** (physical interpretation, not a theorem of this course). A state on
> $`\mathcal M`$ is a classical mixture of quantum states in the sectors. The map
> $`\rho\mapsto\sum_cP_c\rho P_c`$ is the measurement of "which sector", and it does not disturb
> any state of $`\mathcal M`$. This is the operational content of a superselection rule.

*Computation:* the sector probabilities $`p_c`$ (`DensityMatrix.sector_probabilities`) and the measurement $`\rho\mapsto\sum_cP_c\rho P_c`$ (`channels.measurement_channel`), which does not disturb $`\rho`$ — [tutorial, Section 3.3](../tutorial/03_traces_and_states.md#s-3-3).

<a id="s-3-4"></a>

## 3.4 Entropy

The **von Neumann entropy** of a state with (blunt) density $`\rho`$ is
$`S(\rho)=-\mathrm{Tr}\,\rho\log\rho=-\sum_i\lambda_i\log\lambda_i`$ over the eigenvalues of
$`\rho`$, with $`0\log0=0`$; that is, $`\rho\log\rho:=f(\rho)`$ for $`f(t)=t\log t`$, $`f(0)=0`$, and
$`\mathrm{Tr}\,f(\rho)=\sum_if(\lambda_i)`$ by the definition of $`f(\rho)=U\,\mathrm{diag}(f(\lambda_i))\,U^*`$
(Section 1.3).

<a id="r-3-3"></a>

**Theorem 3.3 (entropy of a mixture of sectors).** With $`p_c`$ and $`\sigma_c`$ as above
(sectors with $`p_c=0`$ omitted),
```math
S(\rho)=H(p)+\sum_cp_c\,S(\sigma_c),\qquad H(p)=-\sum_cp_c\log p_c .
```

*Proof.* The nonzero eigenvalues of $`\rho`$ (with multiplicity) are the numbers $`p_c\mu`$ for the nonzero
eigenvalues $`\mu`$ of the $`\sigma_c`$ (sectors with $`p_c>0`$); zero eigenvalues contribute $`0`$
to both sides. Hence
$`-\sum_c\sum_{\mu>0}p_c\mu\log(p_c\mu)=-\sum_cp_c\log p_c\sum_{\mu>0}\mu-\sum_cp_c\sum_{\mu>0}\mu\log\mu`$,
and $`\sum_{\mu>0}\mu=\mathrm{Tr}\,\sigma_c=1`$. $`\square`$

*Computation:* Theorem 3.3 for random states, with `DensityMatrix.von_neumann_entropy` — [tutorial, Section 3.4](../tutorial/03_traces_and_states.md#s-3-4).

> **Physical meaning** (physical interpretation, not a theorem of this course). The formula splits the uncertainty about the system into two parts:
> $`H(p)`$, the classical uncertainty about the charge (which sector the system is in), and
> $`\sum_cp_cS(\sigma_c)`$, the average quantum uncertainty within the sectors. Under a
> superselection rule the first part is the entropy of a classical random variable, the charge.

The entropy relative to another trace is $`S_w(\rho)=-\mathrm{tr}_w(\rho_w\log\rho_w)`$, with the
same convention $`\rho_w\log\rho_w:=f(\rho_w)`$. By
Theorem 3.2, $`(\rho_w)_c=\rho_c/w_c`$, so
```math
S_w(\rho)=-\sum_cw_c\mathrm{Tr}\Big(\frac{\rho_c}{w_c}\log\frac{\rho_c}{w_c}\Big)
=S(\rho)+\sum_cp_c\log w_c .
```
The entropy is a property of the state *and* the trace: the same state has different
entropies with respect to `Tr_blunt` and `tau_vN`.

*Computation:* the entropy relative to $`\tau_{\mathrm{vN}}`$ and the formula for $`S_w(\rho)`$ — [tutorial, Section 3.4](../tutorial/03_traces_and_states.md#s-3-4).

<a id="s-3-5"></a>

## 3.5 Relative entropy and Klein's inequality

For states with densities $`\rho,\sigma`$, $`\sigma>0`$, the **relative entropy** is
$`S(\rho\|\sigma)=\mathrm{Tr}\,\rho\log\rho-\mathrm{Tr}\,\rho\log\sigma`$, where
$`\mathrm{Tr}\,\rho\log\rho:=\sum_i\lambda_i\log\lambda_i`$ over the eigenvalues of $`\rho`$, with
$`0\log0=0`$ (and likewise in each block). It does not depend on the trace: by Theorem 3.2,
$`(\rho_w)_c=\rho_c/w_c`$ and $`(\sigma_w)_c=\sigma_c/w_c`$, so the eigenvalues of $`(\rho_w)_c`$ are those of
$`\rho_c`$ divided by $`w_c`$, $`\log(\sigma_c/w_c)=\log\sigma_c-(\log w_c)1`$ (if $`\sigma_c=U\,\mathrm{diag}(s_i)\,U^*`$, then
$`\sigma_c/w_c=U\,\mathrm{diag}(s_i/w_c)\,U^*`$, and Theorem 1.5(a)), and both $`\mathrm{tr}_w(\rho_w\log\rho_w)`$ and
$`\mathrm{tr}_w(\rho_w\log\sigma_w)`$ differ from $`\mathrm{Tr}\,\rho\log\rho`$ and
$`\mathrm{Tr}\,\rho\log\sigma`$ by the same amount $`-\sum_cp_c\log w_c`$ (the logarithm of $`\rho`$
is only needed on its support, by the convention $`0\log0=0`$). The difference is unchanged.

<a id="r-3-4"></a>

**Lemma 3.4 (tangent inequality).** Let $`f:I\to\mathbb R`$ be convex on an interval $`I`$ and
differentiable at $`b\in I`$ (one-sided if $`b`$ is an endpoint of $`I`$). Then $`f(a)\ge f(b)+(a-b)f'(b)`$ for
all $`a\in I`$, strictly for $`a\neq b`$ if $`f`$ is strictly convex, i.e. if
$`f((1-t)x+ty)<(1-t)f(x)+tf(y)`$ for $`x\neq y`$ in $`I`$ and $`t\in(0,1)`$. Moreover, if $`f`$ is continuous on
$`I`$ and twice differentiable with $`f''>0`$ on the interior of $`I`$, then $`f`$ is strictly convex.

*Proof.* For $`t\in(0,1]`$, convexity gives $`f(b+t(a-b))\le(1-t)f(b)+tf(a)`$, i.e.
$`\big(f(b+t(a-b))-f(b)\big)/t\le f(a)-f(b)`$; letting $`t\to0`$ gives the inequality. If $`f`$
is strictly convex and $`a\neq b`$, let $`m=(a+b)/2`$; then
$`f(a)-f(b)>2\big(f(m)-f(b)\big)\ge2(m-b)f'(b)=(a-b)f'(b)`$, the first step by strict
convexity at $`m`$ and the second by the inequality just proved, applied to $`m`$. For the last
claim let $`x<y`$ in $`I`$, $`t\in(0,1)`$ and $`z=(1-t)x+ty`$. By the mean value theorem there are
$`\xi_1\in(x,z)`$ and $`\xi_2\in(z,y)`$ with $`f(z)-f(x)=t(y-x)f'(\xi_1)`$ and
$`f(y)-f(z)=(1-t)(y-x)f'(\xi_2)`$, and $`f'(\xi_1)<f'(\xi_2)`$ by the mean value theorem for $`f'`$ on
$`[\xi_1,\xi_2]`$, as $`f''>0`$ there. Multiplying by $`t(1-t)(y-x)>0`$ gives
$`(1-t)\big(f(z)-f(x)\big)<t\big(f(y)-f(z)\big)`$, i.e. $`f(z)<(1-t)f(x)+tf(y)`$. $`\square`$

<a id="r-3-5"></a>

**Theorem 3.5 (Klein's inequality).** Let $`f:I\to\mathbb R`$ be convex on an interval $`I`$,
and $`X,Y`$ self-adjoint with spectra in $`I`$, with $`f`$ differentiable at every eigenvalue of $`Y`$
(one-sided at an endpoint of $`I`$).
Then $`\mathrm{Tr}[f(X)-f(Y)-(X-Y)f'(Y)]\ge0`$. If $`f`$ is strictly convex, equality holds
only for $`X=Y`$.

*Proof.* Let $`X=\sum_ia_i|u_i\rangle\langle u_i|`$ and $`Y=\sum_jb_j|v_j\rangle\langle v_j|`$ with
orthonormal eigenbases, and $`c_{ij}=|\langle u_i,v_j\rangle|^2`$, so $`\sum_ic_{ij}=\sum_jc_{ij}=1`$.
Then $`\mathrm{Tr}\,f(X)=\sum_{ij}c_{ij}f(a_i)`$, $`\mathrm{Tr}\,f(Y)=\sum_{ij}c_{ij}f(b_j)`$ and
$`\mathrm{Tr}\,Xf'(Y)=\sum_{ij}c_{ij}a_if'(b_j)`$, $`\mathrm{Tr}\,Yf'(Y)=\sum_{ij}c_{ij}b_jf'(b_j)`$.
Hence the left side equals $`\sum_{ij}c_{ij}[f(a_i)-f(b_j)-(a_i-b_j)f'(b_j)]`$, a sum of
nonnegative terms by Lemma 3.4. If $`f`$ is strictly convex and the sum is zero, then
$`a_i=b_j`$ whenever $`c_{ij}>0`$, so $`Xv_j=\sum_ia_i\langle u_i,v_j\rangle u_i=b_jv_j`$ for every $`j`$,
and $`X=Y`$. $`\square`$

<a id="r-3-6"></a>

**Corollary 3.6.** $`S(\rho\|\sigma)\ge0`$, with equality iff $`\rho=\sigma`$. In particular
$`S(\rho)\le\log D`$, with equality iff $`\rho=1/D`$.

*Proof.* Apply Theorem 3.5 to $`f(t)=t\log t`$, $`f(0)=0`$, which is strictly convex on
$`I=[0,\infty)`$ by Lemma 3.4 ($`f`$ is continuous at $`0`$, as $`t\log t\to0`$ for $`t\to0^+`$, and
$`f''(t)=1/t>0`$ on $`(0,\infty)`$) and
differentiable on $`(0,\infty)\supseteq\mathrm{spec}\,\sigma`$ with $`f'(t)=\log t+1`$, with
$`X=\rho`$, $`Y=\sigma`$: the left side is
$`\mathrm{Tr}\,\rho\log\rho-\mathrm{Tr}\,\sigma\log\sigma-\mathrm{Tr}(\rho-\sigma)(\log\sigma+1) =S(\rho\|\sigma)`$, since $`\mathrm{Tr}\,\rho=\mathrm{Tr}\,\sigma=1`$. For the second claim, take
$`\sigma=1/D`$: $`0\le S(\rho\|1/D)=-S(\rho)+\log D`$. $`\square`$

*Computation:* $`S(\rho)\le\log D`$ and the maximally mixed state — [tutorial, Section 3.5](../tutorial/03_traces_and_states.md#s-3-5).

<a id="r-3-7"></a>

**Theorem 3.7 (chain rule; monotonicity under the conditional expectation onto the centre).** Let $`\rho,\sigma`$ be states with
$`\sigma>0`$, $`p_c=\mathrm{Tr}\,\rho_c`$, $`q_c=\mathrm{Tr}\,\sigma_c>0`$, $`\hat\sigma_c=\sigma_c/q_c`$, and
for $`p_c>0`$, $`\hat\rho_c=\rho_c/p_c`$. Then, with sums over the $`c`$ with $`p_c>0`$ (the other
sectors contribute $`0`$ by $`0\log0=0`$),
```math
S(\rho\|\sigma)=\sum_cp_c\log\frac{p_c}{q_c}+\sum_cp_c\,S(\hat\rho_c\|\hat\sigma_c)
\ \ge\ \sum_cp_c\log\frac{p_c}{q_c}=S(E\rho\|E\sigma),
```
where $`E`$ is the conditional expectation onto the centre (Proposition 2.8).

*Proof.* For $`p_c>0`$, $`\mathrm{Tr}\,\rho_c\log\rho_c=p_c\log p_c+p_c\mathrm{Tr}\,\hat\rho_c\log\hat\rho_c`$
(eigenvalues $`p_c\mu`$, as in Theorem 3.3) and $`\log\sigma_c=(\log q_c)1+\log\hat\sigma_c`$ (diagonalise
$`\hat\sigma_c=U\,\mathrm{diag}(s_i)\,U^*`$, so $`\sigma_c=U\,\mathrm{diag}(q_cs_i)\,U^*`$, and use Theorem 1.5(a)), so
$`\mathrm{Tr}\,\rho_c\log\sigma_c=p_c\log q_c+p_c\mathrm{Tr}\,\hat\rho_c\log\hat\sigma_c`$ and
$`\mathrm{Tr}\,\rho_c\log\rho_c-\mathrm{Tr}\,\rho_c\log\sigma_c=p_c\log(p_c/q_c)+p_cS(\hat\rho_c\|\hat\sigma_c)`$.
If $`p_c=0`$, then $`\rho_c=0`$ and the sector contributes $`0`$.
Sum over $`c`$; the second sum is $`\ge0`$ by Corollary 3.6, applied in $`M_{k_c}`$ to the states
$`\hat\rho_c`$ and $`\hat\sigma_c>0`$. Finally
$`E\rho=\bigoplus_c(p_c/k_c)1_{k_c}`$ and $`E\sigma=\bigoplus_c(q_c/k_c)1_{k_c}>0`$ are states; $`E\rho`$ is the density of
the state $`\omega\circ E`$ (where $`\omega(A)=\mathrm{Tr}\,\rho A`$), since
$`\mathrm{Tr}(\rho E(A))=\sum_cp_c\mathrm{Tr}(A_c)/k_c=\mathrm{Tr}(E(\rho)A)`$. And
$`S(E\rho\|E\sigma)=\sum_ck_c\frac{p_c}{k_c}\big(\log\frac{p_c}{k_c}-\log\frac{q_c}{k_c}\big)=\sum_cp_c\log\frac{p_c}{q_c}`$, the eigenvalues of
$`E\rho`$ in block $`c`$ being $`p_c/k_c`$ ($`k_c`$ times) and sectors with $`p_c=0`$ contributing $`0`$. $`\square`$

*Computation:* Theorem 3.7 for random states, with `states.relative_entropy` and `channels.center_expectation` — [tutorial, Section 3.5](../tutorial/03_traces_and_states.md#s-3-5).

So replacing both states by their images under $`E`$ can only decrease the relative entropy; what
remains is the Kullback–Leibler divergence $`\sum_cp_c\log(p_c/q_c)`$ of the sector distributions.

> **Physical meaning** (physical interpretation, not a theorem of this course). Relative entropy
> measures how well two states can be distinguished. Forgetting the quantum information inside the
> sectors can only make two states harder to distinguish; what remains is the classical
> distinguishability of the charge distributions.

<a id="s-3-6"></a>

## 3.6 The Michelson contrast

For $`A\ge0`$, $`A\neq0`$, the **Michelson contrast** is
$`\Delta(A)=(\lambda_{\max}-\lambda_{\min})/(\lambda_{\max}+\lambda_{\min})\in[0,1]`$, with the
extreme eigenvalues taken over all sectors (and $`\Delta(0)=0`$). It measures how far $`A`$ is
from a multiple of the identity: $`\Delta(A)=0`$ iff $`A=\lambda1`$.

<a id="r-3-8"></a>

**Proposition 3.8.** If $`\Phi:\mathcal M\to\mathcal M`$ is linear, positive and unital
($`\Phi(1)=1`$), then $`\Delta(\Phi(A))\le\Delta(A)`$ for all $`A\ge0`$.

*Proof.* $`\lambda_{\min}1\le A\le\lambda_{\max}1`$ (Corollary 1.4(3) for $`A-\lambda_{\min}1`$ and
$`\lambda_{\max}1-A`$), and positivity, linearity and unitality give
$`\lambda_{\min}1\le\Phi(A)\le\lambda_{\max}1`$. Hence $`\Phi(A)`$ is self-adjoint (a positive matrix is,
Corollary 1.4(3)),
$`\Phi(A)\ge0`$ (as $`\lambda_{\min}\ge0`$), and by Corollary 1.4(3) $`\lambda_{\min}(\Phi(A))\ge\lambda_{\min}`$ and
$`\lambda_{\max}(\Phi(A))\le\lambda_{\max}`$. If $`\Phi(A)=0`$ the claim is trivial. Otherwise $`A\neq0`$, both
$`\lambda_{\max}`$ and $`\lambda_{\max}(\Phi(A))`$ are $`>0`$, and with $`r=\lambda_{\min}/\lambda_{\max}\in[0,1]`$,
$`\Delta=(1-r)/(1+r)`$ is decreasing in $`r`$, and $`r(\Phi(A))\ge r(A)`$. $`\square`$

*Computation:* Proposition 3.8 for a random mixed-unitary channel and for $`E`$, with `Operator.michelson_contrast` — [tutorial, Section 3.6](../tutorial/03_traces_and_states.md#s-3-6).

<a id="s-3-7"></a>

## 3.7 A plausible inequality that is false

For positive $`X,Y`$ in $`M_k`$ and a unitary $`U`$, is $`|\mathrm{Tr}(XUY)|\le\mathrm{Tr}(XY)`$? It
holds for $`U=1`$ and looks like a Cauchy–Schwarz-type statement. It is false whenever $`X`$
and $`Y`$ do not commute.

<a id="r-3-9"></a>

**Proposition 3.9.** Let $`X,Y\ge0`$ in $`M_k`$. Then $`\max_U|\mathrm{Tr}(XUY)|=\|YX\|_1\ge\mathrm{Tr}(XY)`$, the maximum taken over the unitary
$`U\in M_k`$, with equality iff $`XY=YX`$.

*Proof.* $`\mathrm{Tr}(XUY)=\mathrm{Tr}(YX\,U)`$, so the maximum is $`\|YX\|_1`$ by Theorem 1.9,
and $`\|YX\|_1\ge|\mathrm{Tr}(YX)|=\mathrm{Tr}(XY)`$, since
$`\mathrm{Tr}(XY)=\mathrm{Tr}(X^{1/2}YX^{1/2})\ge0`$, as
$`X^{1/2}YX^{1/2}=(Y^{1/2}X^{1/2})^*(Y^{1/2}X^{1/2})\ge0`$ (Corollary 1.7). If $`XY=YX`$, then
$`X^{1/2}`$ commutes with $`Y`$ (Theorem 1.5(d)), so $`YX=X^{1/2}YX^{1/2}\ge0`$ and
$`\|YX\|_1=\mathrm{Tr}(YX)`$ (for $`P\ge0`$, $`\|P\|_1=\mathrm{Tr}\,P`$, Section 1.4). Conversely
suppose $`\|YX\|_1=\mathrm{Tr}(YX)`$. If $`|YX|=0`$ then $`YX=0=(YX)^*=XY`$. Otherwise write
$`YX=W|YX|`$ (Theorem 1.8) and $`R=|YX|^{1/2}\neq0`$. By cyclicity,
$`\mathrm{Tr}(YX)=\mathrm{Tr}(RWR)=\langle R,WR\rangle_{HS}`$, while
$`\|R\|_2^2=\|WR\|_2^2=\mathrm{Tr}|YX|=\|YX\|_1`$. So $`\langle R,WR\rangle_{HS}=\|R\|_2\|WR\|_2`$: equality
holds in the Cauchy–Schwarz inequality, hence $`WR=\lambda R`$ for a scalar $`\lambda`$, and
$`\lambda\|R\|_2^2=\langle R,WR\rangle_{HS}=\|R\|_2^2`$ gives $`\lambda=1`$. So $`W|YX|^{1/2}=|YX|^{1/2}`$ and $`YX=W|YX|=|YX|\ge0`$. A positive
matrix is self-adjoint, so $`YX=(YX)^*=XY`$. $`\square`$

By Theorem 1.9, for $`YX=W|YX|`$ the unitary $`U=W^*`$ attains the maximum $`\|YX\|_1`$. A random
unitary does worse on average: for a Haar-distributed $`U`$ (Theorem 5.3), Proposition 11.11 with $`A=YX`$ gives
$`\big(\mathbb E|\mathrm{Tr}(XUY)|^2\big)^{1/2}=\|YX\|_2/\sqrt k\le\|YX\|_1/\sqrt k=\max_V|\mathrm{Tr}(XVY)|/\sqrt k`$
(as $`\sum_is_i^2\le(\sum_is_i)^2`$), and Corollary 11.12 shows that for $`X,Y`$ with spectra in
$`[a,b]`$, $`a>0`$, a single sample reaches $`\mathrm{Tr}(XY)`$ with probability at most $`b^4/(a^4k^2)`$
(informative only for $`k>b^2/a^2`$).

*Computation:* Proposition 3.9 for a non-commuting pair in $`M_6`$, the maximiser $`W^*`$ and a random search over Haar unitaries — [tutorial, Section 3.7](../tutorial/03_traces_and_states.md#s-3-7).

<a id="s-3-8"></a>

## 3.8 What to compare with

*In the library:* comparison with other software — [tutorial, Section 3.8](../tutorial/03_traces_and_states.md#s-3-8).

<a id="x-3"></a>

## Exercises

1. Show that the density of the tracial state $`\tau_{\mathrm{vN}}`$ (equal sector weights) with respect to
   $`\tau_{\mathrm{vN}}`$ is $`1`$, and compute its blunt density and its entropy $`S`$.
2. Let $`\mathrm{tr}_w`$ be faithful. Show that $`S_w(\rho)\le\log\mathrm{tr}_w(1)`$, with equality iff
   $`\rho_w=1/\mathrm{tr}_w(1)`$.
3. Show that the measurement channel $`\rho\mapsto\sum_cP_c\rho P_c`$ on $`M_D`$ (not on
   $`\mathcal M`$!) does disturb states, and compute $`S(\rho)`$ before and after for a
   superposition of two sectors.

<details>
<summary>Solutions</summary>

1. $`\tau_{\mathrm{vN}}(A)=\tau_{\mathrm{vN}}(1\cdot A)`$, so $`\rho_{\tau}=1`$ by uniqueness. By
   Theorem 3.2 the blunt density is $`\bigoplus_c1_{k_c}/(Ck_c)`$, with $`p_c=1/C`$ and
   $`\sigma_c=1/k_c`$, so $`S=\log C+\frac1C\sum_c\log k_c`$ by Theorem 3.3.

2. Sums run over the $`c`$ with $`p_c>0`$. Then
   $`S_w(\rho)=S(\rho)+\sum_cp_c\log w_c=-\sum_cp_c\log p_c+\sum_cp_cS(\sigma_c)+\sum_cp_c\log w_c    \le\sum_cp_c\log\frac{k_cw_c}{p_c}\le\log\sum_{p_c>0}k_cw_c\le\log\mathrm{tr}_w(1)`$, using
   $`S(\sigma_c)\le\log k_c`$ (Corollary 3.6 in $`M_{k_c}`$) and the finite Jensen inequality
   $`\sum_{c=1}^mp_c\log x_c\le\log\sum_{c=1}^mp_cx_c`$ ($`p_c>0`$, $`\sum_cp_c=1`$, $`x_c>0`$), with
   equality only if all $`x_c`$ are equal, applied to $`x_c=k_cw_c/p_c`$. Proof by induction on $`m`$:
   for $`m=1`$ there is nothing to prove. For $`m\ge2`$ put $`q=1-p_m\in(0,1)`$ and
   $`\bar x=\sum_{c<m}(p_c/q)x_c`$; by the induction hypothesis
   $`\sum_{c<m}(p_c/q)\log x_c\le\log\bar x`$, with equality only if $`x_1=\dots=x_{m-1}`$, and by
   strict concavity of $`\log`$ (Lemma 3.4 for $`-\log`$, whose second derivative is $`1/t^2>0`$)
   $`q\log\bar x+p_m\log x_m\le\log(q\bar x+p_mx_m)`$,
   with equality only if $`\bar x=x_m`$. Multiplying the first inequality by $`q`$ and adding,
   $`\sum_cp_c\log x_c\le\log\sum_cp_cx_c`$, and equality forces $`x_1=\dots=x_{m-1}=\bar x=x_m`$.
   Equality throughout
   forces every $`p_c>0`$ (as all $`k_cw_c>0`$), $`\sigma_c=1/k_c`$ and $`k_cw_c/p_c`$ constant, i.e.
   $`p_c=k_cw_c/\mathrm{tr}_w(1)`$, i.e. $`(\rho_w)_c=\rho_c/w_c=1_{k_c}/\mathrm{tr}_w(1)`$.
   Conversely, for $`\rho_w=1/\mathrm{tr}_w(1)`$,
   $`S_w(\rho)=-\mathrm{tr}_w(1)\cdot\frac1{\mathrm{tr}_w(1)}\log\frac1{\mathrm{tr}_w(1)}=\log\mathrm{tr}_w(1)`$.

3. In $`M_D`$ a pure superposition $`\psi=(\psi_1\oplus\psi_2)/\sqrt2`$ of unit vectors in two
   sectors has $`S=0`$; after the measurement the state is
   $`\frac12|\psi_1\rangle\langle\psi_1|\oplus\frac12|\psi_2\rangle\langle\psi_2|`$ with $`S=\log2`$. Restricted to
   $`\mathcal M`$, however, the vector state of $`\psi`$ coincides with the mixture:
   $`\langle\psi,A\psi\rangle=\tfrac12\langle\psi_1,A_1\psi_1\rangle+\tfrac12\langle\psi_2,A_2\psi_2\rangle`$
   for all $`A\in\mathcal M`$. No element of $`\mathcal M`$ detects the coherence (physically, this is
   the content of a superselection rule).

</details>

*Computation:* Exercise 1 numerically — [tutorial, Exercises](../tutorial/03_traces_and_states.md#x-3).

**Next:** [Chapter 4](04_symmetries.md) treats symmetries and representations: how a
group action on a physical system produces the sectors, and how the library builds
the corresponding algebras for spin systems and qubit registers.
