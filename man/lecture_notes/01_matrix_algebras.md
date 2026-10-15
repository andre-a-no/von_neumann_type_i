**English** · [Español](es/01_matrix_algebras.md) · [Français](fr/01_matrix_algebras.md) · [Deutsch](de/01_matrix_algebras.md) · [中文](zh/01_matrix_algebras.md) · [日本語](ja/01_matrix_algebras.md) · [Русский](ru/01_matrix_algebras.md)

# 1. Matrix algebras and functional calculus

This chapter fixes notation and proves the facts about $`M_n(\mathbb C)`$ that every
later chapter uses: the spectral theorem, the functional calculus and its derivative, the polar
decomposition and the norms defined by singular values. The library code and the numerical
checks of these results are in [Tutorial 1](../tutorial/01_matrix_algebras.md).

**Conventions.** Inner products $`\langle x,y\rangle=\sum_i\bar x_iy_i`$ are linear in the
second argument. $`A^*`$ is the conjugate transpose, so $`\langle x,Ay\rangle=\langle A^*x,y\rangle`$.
$`e_1,\dots,e_n`$ is the standard basis, $`1`$ the identity matrix, $`\mathrm{Tr}`$ the
usual trace.

<a id="s-1-1"></a>

## 1.1 The algebra $`M_n(\mathbb C)`$

$`M_n(\mathbb C)`$ is the algebra of complex $`n\times n`$ matrices with the involution
$`A\mapsto A^*`$ and the operator norm $`\|A\|=\sup_{\|x\|=1}\|Ax\|`$. An element $`A`$ is

| name | definition |
|---|---|
| normal | $`A^*A=AA^*`$ |
| self-adjoint | $`A=A^*`$ |
| positive, $`A\ge0`$ | $`\langle x,Ax\rangle\ge0`$ for all $`x`$ |
| unitary | $`A^*A=AA^*=1`$ |
| projection | $`A=A^*=A^2`$ |

We write $`A\le B`$ if $`B-A\ge0`$, and $`A>0`$ if $`A\ge0`$ and $`A`$ is invertible.

*In the library:* operators with a prescribed spectrum, self-adjoint or positive, `TypeIAlgebra.operator_from_eigenvalues` — [tutorial, Section 1.5](../tutorial/01_matrix_algebras.md#s-1-5); the operator norm, `operator_norm` — [tutorial, Section 1.6](../tutorial/01_matrix_algebras.md#s-1-6).

<a id="r-1-1"></a>

**Lemma 1.1 (C\*-identity).** $`\|A^*\|=\|A\|`$ and $`\|A^*A\|=\|A\|^2`$.

*Proof.* $`\|Ax\|=\sup_{\|y\|=1}|\langle y,Ax\rangle|=\sup_{\|y\|=1}|\langle A^*y,x\rangle|`$,
so $`\|A\|=\sup_{\|x\|=\|y\|=1}|\langle A^*y,x\rangle|=\|A^*\|`$. Next,
$`\|Ax\|^2=\langle x,A^*Ax\rangle\le\|A^*A\|`$ for $`\|x\|=1`$, hence $`\|A\|^2\le\|A^*A\|\le\|A^*\|\,\|A\|=\|A\|^2`$. $`\square`$

> **Physical meaning** (physical interpretation, not a theorem of this course). In quantum
> mechanics a system with $`n`$ distinguishable configurations is described by $`\mathbb C^n`$.
> Self-adjoint matrices are *observables* (energy, spin, position on a lattice); their
> eigenvalues are the possible results of a measurement. Projections are yes/no questions. Unitary
> matrices are symmetries and time evolutions, $`e^{-itH}`$ for a Hamiltonian $`H`$.
> Positive matrices of trace one are *states* (density matrices, Chapter 3).

<a id="s-1-2"></a>

## 1.2 The spectral theorem

<a id="r-1-2"></a>

**Lemma 1.2 (Schur).** For every $`A\in M_n(\mathbb C)`$ there is a unitary $`U`$ with
$`U^*AU`$ upper triangular.

*Proof.* Induction on $`n`$; $`n=1`$ is trivial. The characteristic polynomial of $`A`$
has a root $`\lambda`$ (fundamental theorem of algebra), so $`Av=\lambda v`$ for some
unit vector $`v`$. Complete $`v`$ to an orthonormal basis $`v,w_2,\dots,w_n`$ and let
$`V`$ be the unitary with these columns. Then
$`V^*AV=\begin{pmatrix}\lambda&*\\0&A'\end{pmatrix}`$ with $`A'\in M_{n-1}(\mathbb C)`$.
By induction $`U'^*A'U'`$ is upper triangular for a unitary $`U'`$, and
$`U=V\,(1\oplus U')`$ does the job. $`\square`$

<a id="r-1-3"></a>

**Theorem 1.3 (spectral theorem).** $`A`$ is normal if and only if
$`A=U\,\mathrm{diag}(\lambda_1,\dots,\lambda_n)\,U^*`$ for a unitary $`U`$. The
$`\lambda_i`$ are the eigenvalues of $`A`$ with multiplicity.

*Proof.* If $`A=U\Lambda U^*`$ with $`\Lambda`$ diagonal, then
$`A^*A=U\bar\Lambda\Lambda U^*=U\Lambda\bar\Lambda U^*=AA^*`$. Conversely, let $`A`$
be normal and $`T=U^*AU`$ upper triangular (Lemma 1.2); $`T`$ is normal. Compare the
$`(1,1)`$ entries: $`(T^*T)_{11}=|t_{11}|^2`$ and $`(TT^*)_{11}=\sum_j|t_{1j}|^2`$, so
$`t_{1j}=0`$ for $`j>1`$. Then $`T=t_{11}\oplus T'`$ with $`T'`$ normal and upper
triangular, and induction shows that $`T`$ is diagonal. The diagonal of $`U^*AU`$
consists of the roots of $`\det(A-\lambda)=\det(T-\lambda)`$. $`\square`$

*In the library:* Theorem 1.3 read backwards, $`A=U\Lambda U^*`$ with prescribed eigenvalues and a Haar-random $`U`$, is `TypeIAlgebra.operator_from_eigenvalues` — [tutorial, Section 1.5](../tutorial/01_matrix_algebras.md#s-1-5).

<a id="r-1-4"></a>

**Corollary 1.4.** Let $`A`$ be normal with eigenvalues $`\lambda_i`$. Then

1. $`\|A\|=\max_i|\lambda_i|`$;
2. $`A`$ is self-adjoint iff all $`\lambda_i\in\mathbb R`$; unitary iff all $`|\lambda_i|=1`$;
   a projection iff all $`\lambda_i\in\{0,1\}`$;
3. $`A\ge0`$ iff $`A`$ is self-adjoint and all $`\lambda_i\ge0`$. Moreover every positive matrix is
   self-adjoint, hence normal, so (3) characterises positivity among *all* matrices and the
   functional calculus below applies to every $`A\ge0`$.

*Proof.* With $`y=U^*x`$, $`\|Ax\|^2=\sum_i|\lambda_i|^2|y_i|^2`$ and $`\|y\|=\|x\|`$,
which gives (1), with equality for $`x=Ue_i`$, $`|\lambda_i|`$ maximal. In (2), $`A^*=U\bar\Lambda U^*`$,
$`A^*A=U|\Lambda|^2U^*`$ and $`A^2=U\Lambda^2U^*`$, and the conditions hold iff they hold
entrywise for $`\Lambda`$. For (3), first: every positive $`A`$ is self-adjoint, since
$`\langle x,(A-A^*)x\rangle=\langle x,Ax\rangle-\overline{\langle x,Ax\rangle}=0`$ for all $`x`$,
and a matrix $`D`$ with $`\langle x,Dx\rangle=0`$ for all $`x`$ vanishes (polarisation:
$`4\langle x,Dy\rangle=\sum_{k=0}^3 i^{-k}\langle x+i^ky,D(x+i^ky)\rangle`$). Then
$`\langle x,Ax\rangle=\sum_i\lambda_i|y_i|^2`$, which is $`\ge0`$ for all $`x`$ iff all $`\lambda_i\ge0`$. $`\square`$

<a id="s-1-3"></a>

## 1.3 The functional calculus

Let $`A=U\Lambda U^*`$ be normal and $`f`$ a complex function on $`\mathrm{spec}\,A`$. Define
```math
f(A)=U\,\mathrm{diag}\big(f(\lambda_1),\dots,f(\lambda_n)\big)\,U^*.
```

<a id="r-1-5"></a>

**Theorem 1.5.** (a) $`f(A)`$ does not depend on the choice of $`U`$. (b) For a
polynomial $`p`$, $`p(A)`$ is the usual polynomial in $`A`$. (c) $`f\mapsto f(A)`$ is
linear and multiplicative, $`\bar f(A)=f(A)^*`$, and $`\mathrm{spec}\,f(A)=f(\mathrm{spec}\,A)`$.
(d) If $`B`$ commutes with $`A`$, then $`B`$ commutes with $`f(A)`$. (e) If $`g`$ is defined on
$`f(\mathrm{spec}\,A)`$, then $`g(f(A))=(g\circ f)(A)`$.

*Proof.* (b) $`A^k=U\Lambda^kU^*`$, so $`p(A)=Up(\Lambda)U^*`$ by linearity. (a) The
spectrum is finite, so Lagrange interpolation gives a polynomial $`p`$ with $`p=f`$ on
$`\mathrm{spec}\,A`$: if $`\mu_1,\dots,\mu_r`$ are the distinct eigenvalues, take
$`p(t)=\sum_if(\mu_i)\prod_{j\neq i}(t-\mu_j)/(\mu_i-\mu_j)`$. Then $`U\,f(\Lambda)\,U^*=U\,p(\Lambda)\,U^*=p(A)`$ for *every*
unitary $`U`$ diagonalising $`A`$, and $`p(A)`$ is defined without $`U`$. (c) holds for the
diagonal matrices $`f(\Lambda)`$ and is preserved by conjugation with $`U`$; the spectrum
of $`Uf(\Lambda)U^*`$ is the diagonal of $`f(\Lambda)`$. (d) $`B`$ commutes with every
polynomial in $`A`$, in particular with $`p(A)=f(A)`$. (e) $`f(A)=Uf(\Lambda)U^*`$ is normal and
diagonalised by the same $`U`$, so by (a) $`g(f(A))=U\,g(f(\Lambda))\,U^*=(g\circ f)(A)`$. $`\square`$

<a id="r-1-6"></a>

**Theorem 1.6 (square root).** For $`A\ge0`$ there is exactly one $`B\ge0`$ with
$`B^2=A`$, namely $`B=\sqrt A`$ defined by the functional calculus.

*Proof.* $`\sqrt A\ge0`$ and $`(\sqrt A)^2=A`$ by Theorem 1.5(c) and Corollary 1.4(3).
Let $`B\ge0`$ with $`B^2=A`$, $`B=V\,\mathrm{diag}(\mu_j)\,V^*`$, $`\mu_j\ge0`$. Then
$`A=V\,\mathrm{diag}(\mu_j^2)\,V^*`$, so $`\mathrm{spec}\,A=\{\mu_j^2\}`$. Let $`p`$ be a
polynomial with $`p(t)=\sqrt t`$ on $`\mathrm{spec}\,A`$. Then
$`p(A)=\sqrt A`$ and $`p(A)=p(B^2)=V\,\mathrm{diag}(p(\mu_j^2))\,V^*=V\,\mathrm{diag}(\mu_j)\,V^*=B`$. $`\square`$

*Computation:* Theorems 1.5 and 1.6 checked with `apply_function`, `sqrt`, `log`, `expm` and `power` — [tutorial, Section 1.5](../tutorial/01_matrix_algebras.md#s-1-5).

<a id="r-1-7"></a>

**Corollary 1.7.** $`A\ge0`$ iff $`A=T^*T`$ for some $`T`$. (If $`A\ge0`$ take $`T=\sqrt A`$;
conversely $`\langle x,T^*Tx\rangle=\|Tx\|^2\ge0`$.)

**Derivatives.** Gradient methods differentiate $`f(A)`$ with respect to $`A`$. For self-adjoint $`A`$ the
derivative is an entrywise product with a matrix of difference quotients of $`f`$.

<a id="r-1-11"></a>

**Theorem 1.11 (Daleckii–Krein).** Let $`A=U\,\mathrm{diag}(\lambda_1,\dots,\lambda_n)\,U^*`$ be self-adjoint, $`I`$ an
open interval containing $`\mathrm{spec}\,A`$ and $`f:I\to\mathbb C`$ continuously differentiable. Let
```math
f^{[1]}(x,y)=\frac{f(x)-f(y)}{x-y}\quad(x\ne y),\qquad f^{[1]}(x,x)=f'(x)
```
be the **divided difference** of $`f`$ and $`L_{ij}=f^{[1]}(\lambda_i,\lambda_j)`$ the **Loewner matrix**. Then, for
every self-adjoint $`E`$,
```math
\frac{d}{dt}\Big|_{t=0}f(A+tE)=U\big(L\circ U^*EU\big)U^*,
```
where $`(X\circ Y)_{ij}=X_{ij}Y_{ij}`$ is the entrywise (Schur) product.

*Proof.* By the fundamental theorem of calculus, $`f^{[1]}(x,y)=\int_0^1f'(y+s(x-y))\,ds`$ for all $`x,y\in I`$ (for
$`x=y`$ both sides are $`f'(x)`$), so $`f^{[1]}`$ is continuous on $`I\times I`$. Reordering the columns of $`U`$ permutes
the rows and columns of $`L`$ and of $`U^*EU`$ alike and does not change the right-hand side, so we may assume
$`\lambda_1\le\dots\le\lambda_n`$. Let $`\delta>0`$ with $`[\lambda_1-\delta,\lambda_n+\delta]\subset I`$. For
$`0<|t|\,\|E\|<\delta`$ write $`A+tE=V\,\mathrm{diag}(\mu_1,\dots,\mu_n)\,V^*`$ with $`\mu_1\le\dots\le\mu_n`$; then
$`|\mu_j-\lambda_j|\le|t|\,\|E\|`$ (Theorem A.10), so all $`\mu_j\in I`$. Put $`W=U^*V`$ (unitary) and $`\tilde E=U^*EU`$.
As $`U^*f(A)=\mathrm{diag}(f(\lambda_i))\,U^*`$ and $`f(A+tE)V=V\,\mathrm{diag}(f(\mu_j))`$,
```math
\big(U^*(f(A+tE)-f(A))V\big)_{ij}=\big(f(\mu_j)-f(\lambda_i)\big)W_{ij}=f^{[1]}(\lambda_i,\mu_j)\,(\mu_j-\lambda_i)\,W_{ij},
```
and the same computation for $`f(x)=x`$ gives $`(\mu_j-\lambda_i)W_{ij}=t\,(\tilde EW)_{ij}`$. With
$`(L_t)_{ij}=f^{[1]}(\lambda_i,\mu_j)`$ this is the exact identity
```math
\frac1t\big(f(A+tE)-f(A)\big)=U\big((L_t\circ\tilde EW)\,W^*\big)U^*.
```
As $`t\to0`$, $`L_t\to L`$ entrywise, because $`f^{[1]}`$ is continuous and $`\mu_j\to\lambda_j`$. Let $`D_{ij}=W_{ij}`$ if
$`\lambda_i=\lambda_j`$, $`D_{ij}=0`$ otherwise, and $`R=W-D`$. If $`\lambda_i\ne\lambda_j`$ and
$`|t|\,\|E\|\le|\lambda_i-\lambda_j|/2`$, then $`|\mu_j-\lambda_i|\ge|\lambda_i-\lambda_j|/2`$ and
$`|(\tilde EW)_{ij}|\le\|\tilde EW\|=\|E\|`$, so $`|W_{ij}|\le2|t|\,\|E\|/|\lambda_i-\lambda_j|`$. Hence $`R\to0`$ and
$`DD^*=(W-R)(W-R)^*\to WW^*=1`$. Since $`L_{ij}=L_{ik}`$ whenever $`\lambda_k=\lambda_j`$, and $`D_{kj}=0`$ unless
$`\lambda_k=\lambda_j`$, we have $`(L\circ\tilde ED)_{ij}=\sum_kL_{ik}\tilde E_{ik}D_{kj}=((L\circ\tilde E)D)_{ij}`$. Therefore
```math
(L_t\circ\tilde EW)W^*-L\circ\tilde E=\big((L_t-L)\circ\tilde EW\big)W^*+(L\circ\tilde ER)W^*+(L\circ\tilde ED)R^*+(L\circ\tilde E)(DD^*-1),
```
and each term tends to $`0`$, as the entries of $`W`$, $`D`$ and $`\tilde E`$ are bounded by $`1`$, $`1`$ and $`\|E\|`$. $`\square`$

For real $`f`$, $`L`$ is real, and $`\mathrm{Tr}\big(G^*U(L\circ U^*EU)U^*\big)=\mathrm{Tr}(\Gamma^*E)`$ with
$`\Gamma=U(L\circ U^*GU)U^*`$, for every $`G`$ (write both sides as sums over $`i,j`$ in the basis $`U`$): the derivative
of $`A\mapsto\mathrm{Re}\,\mathrm{Tr}(G^*f(A))`$ in the direction $`E`$ is $`\mathrm{Re}\,\mathrm{Tr}(\Gamma^*E)`$.

*In the library:* `apply_function` is differentiated (by autograd) with Theorem 1.11; the Loewner matrix is
`functional.divided_differences`.

**Rounding.** For close $`\lambda_i,\lambda_j`$ the quotient $`(f(\lambda_i)-f(\lambda_j))/(\lambda_i-\lambda_j)`$ loses
accuracy by cancellation, and for $`\lambda_i=\lambda_j`$ it is not defined. The derivatives are then the better data.

<a id="r-1-12"></a>

**Proposition 1.12 (midpoint value).** Let $`a<b`$, $`h=b-a`$, and $`f`$ three times continuously differentiable
on $`[a,b]`$ (complex values allowed). Then
```math
f^{[1]}(a,b)-\frac{f'(a)+f'(b)}2=-\frac{h^2}2\int_0^1s(1-s)\,f'''(a+sh)\,ds,
\qquad
\Big|f^{[1]}(a,b)-\frac{f'(a)+f'(b)}2\Big|\le\frac{h^2}{12}\max_{[a,b]}|f'''|,
```
with equality for polynomials of degree $`3`$.

*Proof.* Let $`g(s)=f'(a+sh)`$, so that $`f^{[1]}(a,b)=\frac1h\int_a^bf'=\int_0^1g`$, and $`\varphi(s)=s(1-s)/2`$, so that
$`\varphi(0)=\varphi(1)=0`$, $`\varphi'(0)=\frac12=-\varphi'(1)`$ and $`\varphi''=-1`$. Integrating by parts twice,
$`\int_0^1\varphi g''=-\int_0^1\varphi'g'=-[\varphi'g]_0^1+\int_0^1\varphi''g=\frac{g(0)+g(1)}2-\int_0^1g`$. As
$`g''(s)=h^2f'''(a+sh)`$, this is the identity; $`\int_0^1\varphi=\frac1{12}`$ gives the bound, and for a cubic
$`f'''`$ is constant. $`\square`$

<a id="r-1-13"></a>

**Proposition 1.13 (quotient or midpoint).** Assume floating-point arithmetic with unit roundoff $`u\le1/20`$:
for floating-point numbers $`x,y`$ and $`\circ\in\{+,-,\times,/\}`$ the computed result is $`(x\circ y)(1+\theta)`$
with $`|\theta|\le u`$ (no overflow or underflow). Let $`a\ne b`$ be real floating-point numbers, $`h=|a-b|`$,
$`x=f^{[1]}(a,b)`$, and let $`\hat f_a,\hat f_b`$ be computed values with $`|\hat f_a-f(a)|\le\delta_a`$ and
$`|\hat f_b-f(b)|\le\delta_b`$.

(a) The computed quotient $`\hat q`$ of $`\hat f_a-\hat f_b`$ and $`a-b`$ satisfies
$`|\hat q-x|\le R:=(\delta_a+\delta_b)/h+5u|\hat q|`$.

(b) Let $`|\hat m-x|\le T`$ (for instance, $`\hat m`$ the computed midpoint value and $`T`$ the bound of
Proposition 1.12 plus its rounding error). Put $`y=\hat m`$ if $`|\hat q-\hat m|\le2R`$, and $`y=\hat q`$ otherwise.
Then $`|y-x|\le\min(T,3R)`$.

(c) If $`\delta_a+\delta_b\le\Delta`$ and $`T=Mh^2+T_0`$, then $`|y-x|\le(3\Delta)^{2/3}M^{1/3}+T_0+15u|\hat q|`$, for
every $`h`$. The two error terms $`3\Delta/h`$ and $`Mh^2`$ balance at $`h_*=(3\Delta/M)^{1/3}`$.

*Proof.* (a) Let $`\tilde q=(\hat f_a-\hat f_b)/(a-b)`$. The three operations give
$`\hat q=\tilde q\,(1+\theta_1)(1+\theta_3)/(1+\theta_2)=\tilde q\,(1+\theta)`$ with $`|\theta_i|\le u`$, hence
$`-3u\le\frac{(1-u)^2}{1+u}-1\le\theta\le\frac{(1+u)^2}{1-u}-1=\frac{3u+u^2}{1-u}\le4u`$. So
$`|\hat q-\tilde q|=|\theta|\,|\tilde q|\le\frac{4u}{1-4u}|\hat q|\le5u|\hat q|`$, and
$`|\tilde q-x|=|(\hat f_a-f(a))-(\hat f_b-f(b))|/h\le(\delta_a+\delta_b)/h`$. (b) If $`y=\hat m`$, then $`|y-x|\le T`$
and $`|y-x|\le|\hat m-\hat q|+|\hat q-x|\le3R`$. If $`y=\hat q`$, then $`|y-x|\le R`$, and
$`T\ge|\hat m-x|\ge|\hat m-\hat q|-|\hat q-x|>2R-R=R`$. (c) $`\min(T,3R)\le\min(Mh^2,3\Delta/h)+T_0+15u|\hat q|`$; for
$`h\ge h_*`$, $`3\Delta/h\le3\Delta/h_*`$, for $`h<h_*`$, $`Mh^2<Mh_*^2`$, and $`3\Delta/h_*=Mh_*^2=(3\Delta)^{2/3}M^{1/3}`$. $`\square`$

The proof holds verbatim for complex-valued $`f`$ (real and imaginary parts are rounded separately, and $`a-b`$ is
real; $`\theta`$ is then complex with the same bound). The rule (b) needs $`R`$, which is computed from the data, but
not $`T`$: without any knowledge of $`f'''`$, its error is at most three times the smaller of the two bounds. A fixed
threshold, midpoint value if and only if $`h\le h_*`$, would need $`M=\max|f'''|/12`$. The threshold
$`h\le\varepsilon^{1/3}\max(|a|,|b|)`$ ($`\varepsilon=2u`$) of earlier versions of the library is of this kind, with
$`M`$ of the order of $`\Delta/\max(|a|,|b|)^3`$: right for $`x^p`$ or $`\log x`$, whose derivatives scale with $`|x|`$,
wrong for functions that vary on a scale unrelated to $`|x|`$. For $`f(x)=e^{x-c}`$ at $`a\approx c=10^4`$ it took the
midpoint value up to $`h\approx0.06`$, with relative error up to $`h^2/12\approx3\cdot10^{-4}`$ (Proposition 1.12 is an
equality up to a factor $`e^{h}`$ there); the same happens for $`e^{-i\tau x}`$ and a Hamiltonian with a large
constant shift.

*In the library:* `functional.divided_differences` applies rule (b) to every pair of eigenvalues, with
$`\delta_a=\varepsilon\,(|\hat f_a|+|a\hat f'_a|)`$ and $`4\varepsilon|\hat q|`$ in place of $`5u|\hat q|`$. The choice of
$`\delta_a`$ is an assumption on the evaluation of $`f`$, not a theorem: relative error at most $`\varepsilon`$, at an
argument perturbed by at most $`\varepsilon|a|`$ (to first order; the second term covers functions such as
$`e^{-i\tau x}`$, which round $`\tau x`$). The same rule is applied to the complex eigenvalues of normal matrices,
which the proof does not cover. The midpoint value is also taken for ties and for pairs that the eigensolver does
not resolve, $`h\le n\varepsilon\max_l|\lambda_l|`$: there its error is at most $`(n\varepsilon\max_l|\lambda_l|)^2\max|f'''|/12`$
by Proposition 1.12, while the quotient would be meaningless for a discontinuous $`f`$ such as an indicator function.

*Computation:* Theorem 1.11 against a difference quotient at a degenerate spectrum, and Propositions 1.12 and 1.13
for $`e^{x-c}`$ with $`c=10^4`$ — [tutorial, Section 1.5](../tutorial/01_matrix_algebras.md#s-1-5).

<a id="s-1-4"></a>

## 1.4 Polar decomposition and singular values

For any $`A`$, the matrix $`A^*A`$ is positive (Corollary 1.7); its square root
$`|A|=(A^*A)^{1/2}`$ is the **absolute value** of $`A`$.

<a id="r-1-8"></a>

**Theorem 1.8 (polar decomposition).** Every $`A\in M_n(\mathbb C)`$ can be written
$`A=W|A|`$ with $`W`$ unitary.

*Proof.* For every $`x`$, $`\||A|x\|^2=\langle x,|A|^2x\rangle=\langle x,A^*Ax\rangle=\|Ax\|^2`$.
Hence $`W_0:|A|x\mapsto Ax`$ is a well-defined linear isometry from the range
$`R_1`$ of $`|A|`$ onto the range $`R_2`$ of $`A`$: if $`|A|x=|A|x'`$ then
$`\|A(x-x')\|=\||A|(x-x')\|=0`$. In particular $`\dim R_1=\dim R_2`$, so the orthogonal
complements also have equal dimension; let $`W_1:R_1^\perp\to R_2^\perp`$ be any
isometry between them (map an orthonormal basis to an orthonormal basis). Let
$`W=W_0\oplus W_1`$, i.e. $`W(u_1+u_2)=W_0u_1+W_1u_2`$ for $`u_1\in R_1`$, $`u_2\in R_1^\perp`$. As
$`W_0u_1\in R_2\perp R_2^\perp\ni W_1u_2`$, Pythagoras gives $`\|Wu\|^2=\|u_1\|^2+\|u_2\|^2=\|u\|^2`$
for all $`u`$. Hence $`\langle u,(W^*W-1)u\rangle=0`$ for all $`u`$, so $`W^*W=1`$ by polarisation (proof of
Corollary 1.4(3)); then $`W`$ is invertible with $`W^{-1}=W^*`$, so also $`WW^*=1`$: $`W`$ is unitary, and
$`W|A|x=W_0|A|x=Ax`$. $`\square`$

*Computation:* Theorem 1.8 on a Ginibre matrix, with $`|A|`$ computed by `abs` — [tutorial, Section 1.6](../tutorial/01_matrix_algebras.md#s-1-6).

The eigenvalues $`s_1\ge\dots\ge s_n\ge0`$ of $`|A|`$ are the **singular values** of
$`A`$. Three norms are used throughout:
```math
\|A\|=s_1,\qquad \|A\|_2=\Big(\sum_is_i^2\Big)^{1/2}=(\mathrm{Tr}\,A^*A)^{1/2},
\qquad\|A\|_1=\sum_is_i=\mathrm{Tr}|A|.
```
The first identity holds because $`\|Ax\|=\||A|x\|`$ for every $`x`$, so
$`\|A\|=\||A|\|=s_1`$ by Corollary 1.4(1); the second because
$`\mathrm{Tr}\,A^*A=\mathrm{Tr}|A|^2=\sum_is_i^2`$.

<a id="r-1-9"></a>

**Theorem 1.9 (duality of the trace norm).**
$`\|A\|_1=\max_{U\text{ unitary}}|\mathrm{Tr}(AU)|`$, and the maximum is attained at
$`U=W^*`$, where $`A=W|A|`$. In particular $`|\mathrm{Tr}\,A|\le\|A\|_1`$.

*Proof.* $`\mathrm{Tr}(AW^*)=\mathrm{Tr}(W|A|W^*)=\mathrm{Tr}|A|`$. For any unitary $`U`$,
by cyclicity and Cauchy–Schwarz for $`\langle X,Y\rangle_{HS}=\mathrm{Tr}(X^*Y)`$,
```math
|\mathrm{Tr}(AU)|=|\mathrm{Tr}(|A|^{1/2}\,UW\,|A|^{1/2})|
=|\langle|A|^{1/2},\,UW|A|^{1/2}\rangle_{HS}|
\le\||A|^{1/2}\|_2\,\|UW|A|^{1/2}\|_2=\mathrm{Tr}|A|,
```
since $`\|UW X\|_2=\|X\|_2`$ for unitary $`UW`$ and $`\||A|^{1/2}\|_2^2=\mathrm{Tr}|A|`$. $`\square`$

In particular the three quantities above are norms. $`\|\cdot\|`$ is an operator norm;
$`\|\cdot\|_2`$ is the norm of the inner product $`\langle X,Y\rangle_{HS}`$; and for
$`\|\cdot\|_1`$, Theorem 1.9 gives
$`\|A+B\|_1=\max_U|\mathrm{Tr}(AU)+\mathrm{Tr}(BU)|\le\|A\|_1+\|B\|_1`$ and
$`\|cA\|_1=\max_U|c\mathrm{Tr}(AU)|=|c|\,\|A\|_1`$, while
$`\|A\|_1=0`$ forces $`|A|=0`$ and $`A=W|A|=0`$. For $`P\ge0`$, $`|P|=(P^2)^{1/2}=P`$ by Theorem 1.6,
so $`\|P\|_1=\mathrm{Tr}\,P`$.

*Computation:* the three norms (`operator_norm`, `frobenius_norm`, `trace_norm`) and Theorem 1.9 — [tutorial, Section 1.6](../tutorial/01_matrix_algebras.md#s-1-6).

<a id="s-1-5"></a>

## 1.5 Operators in the library

*In the library:* the algebra $`M_n(\mathbb C)`$ as `TypeIAlgebra([n], [n])`, how its operators are stored, and the functional calculus `apply_function` — [tutorial, Section 1.5](../tutorial/01_matrix_algebras.md#s-1-5).

<a id="s-1-6"></a>

## 1.6 Norms and the trace-norm duality

*Computation:* the norms, the polar decomposition and the duality of Section 1.4, on a Ginibre matrix — [tutorial, Section 1.6](../tutorial/01_matrix_algebras.md#s-1-6).

<a id="s-1-7"></a>

## 1.7 Batches: many operators at once

As an example of a quantity that many random operators at once estimate, we compute an exact moment of
the **Gaussian unitary ensemble** (GUE).

Let $`G`$ have independent entries with independent real and imaginary parts
$`\sim N(0,\tfrac12)`$, so $`\mathbb E|G_{ij}|^2=1`$, and let $`H=(G+G^*)/\sqrt2`$.

<a id="r-1-10"></a>

**Proposition 1.10.** $`\mathbb E\,\mathrm{Tr}\,H^4=2n^3+n`$. Hence for
$`X=H/\sqrt n`$, $`\ \mathbb E\,\tfrac1n\mathrm{Tr}\,X^4=2+\tfrac1{n^2}`$.

*Proof.* *Isserlis' (Wick's) formula* (Theorem A.24, case $`m=4`$). The real and imaginary parts
of the entries of $`G`$ are independent centred Gaussian variables, so they are the coordinates of
one centred Gaussian vector, and the entries of $`H`$ and of $`\bar H`$ are complex linear combinations
of them. Theorem A.24 then gives, for any four such variables $`Z_1,\dots,Z_4`$ (repetitions allowed),
```math
\mathbb E[Z_1Z_2Z_3Z_4]=\mathbb E[Z_1Z_2]\,\mathbb E[Z_3Z_4]+\mathbb E[Z_1Z_3]\,\mathbb E[Z_2Z_4]+\mathbb E[Z_1Z_4]\,\mathbb E[Z_2Z_3],
```
the sum over the three pairings of $`\{1,2,3,4\}`$ of products of pair expectations (without
complex conjugation). The proof in Appendix A compares the power series of the moment generating
function with the Gaussian moments and does not differentiate under the expectation.

*The moment.* The pair expectations are
$`\mathbb E[H_{ij}H_{kl}]=\tfrac12\mathbb E[(G_{ij}+\bar G_{ji})(G_{kl}+\bar G_{lk})]=\delta_{il}\delta_{jk}`$,
because $`\mathbb E[G_{ij}G_{kl}]=0`$ and $`\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}`$ (Proposition A.16(a);
for $`(i,j)\neq(k,l)`$ the entries are independent and centred, and for $`G_{ij}=x+iy`$ one has
$`\mathbb E[G_{ij}^2]=\mathbb E[x^2-y^2]+2i\,\mathbb E[xy]=0`$ and $`\mathbb E|G_{ij}|^2=\tfrac12+\tfrac12=1`$).
In $`\mathrm{Tr}\,H^4=\sum_{ijkl}H_{ij}H_{jk}H_{kl}H_{li}`$ the pairings give:
$`\{12\}\{34\}`$: $`\delta_{ik}\delta_{ki}`$, summing to $`n^3`$ (free indices $`i=k`$, $`j`$, $`l`$);
$`\{14\}\{23\}`$: $`\delta_{jl}\delta_{jl}`$, summing to $`n^3`$;
$`\{13\}\{24\}`$: $`\delta_{il}\delta_{jk}\delta_{ji}\delta_{kl}`$, forcing $`i=j=k=l`$, summing to $`n`$.
Dividing by $`n\cdot n^2`$ gives the second statement. $`\square`$

The two *non-crossing* pairings $`\{12\}\{34\}`$ and $`\{14\}\{23\}`$ (the pairs can be drawn as
non-intersecting arcs above the points $`1,2,3,4`$) give the leading term $`2`$, for every $`n`$; the
crossing pairing $`\{13\}\{24\}`$ gives the $`1/n^2`$ correction.

*Computation:* Proposition 1.10 by Monte Carlo, on a batch of $`200\,000`$ samples — [tutorial, Section 1.7](../tutorial/01_matrix_algebras.md#s-1-7).

<a id="s-1-8"></a>

## 1.8 What to compare with

*In the library:* reference implementations to compare with — [tutorial, Section 1.8](../tutorial/01_matrix_algebras.md#s-1-8).

<a id="x-1"></a>

## Exercises

1. Show that $`(fg)(A)=f(A)g(A)`$ for a normal $`A`$ and functions $`f,g`$ on $`\mathrm{spec}\,A`$.
2. Show that $`\|A\|_1=\|A^*\|_1`$ and $`\|AB\|_1\le\|A\|\,\|B\|_1`$, using Theorem 1.9.
3. The **geometric mean** of $`A,B>0`$ is $`A\#B=A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}`$.
   Prove that $`A\#B`$ is the unique positive solution $`X`$ of $`XA^{-1}X=B`$, and deduce
   $`A\#B=B\#A`$.

<details>
<summary>Solutions</summary>

1. By Theorem 1.5(c).

2. $`|\mathrm{Tr}(A^*U)|=|\overline{\mathrm{Tr}(U^*A)}|=|\mathrm{Tr}(AU^*)|`$, and
   $`U\mapsto U^*`$ is a bijection of the unitary group, so the maxima in Theorem 1.9
   agree. For the second claim, let $`B=V|B|`$. For unitary $`U`$,
   $`|\mathrm{Tr}(ABU)|=|\mathrm{Tr}(|B|^{1/2}UAV|B|^{1/2})|    \le\||B|^{1/2}\|_2\,\|UAV|B|^{1/2}\|_2\le\|A\|\mathrm{Tr}|B|`$, using
   $`\|YX\|_2\le\|Y\|\,\|X\|_2`$ (apply $`\|Yx\|\le\|Y\|\,\|x\|`$ to each column of $`X`$) and
   $`\|UAV\|=\|A\|`$ for unitaries $`U,V`$ (they preserve the length of vectors).

3. Write $`T=A^{-1/2}XA^{-1/2}`$. Then $`XA^{-1}X=B`$ iff $`T^2=A^{-1/2}BA^{-1/2}`$, and
   $`X\ge0`$ iff $`T\ge0`$, since $`\langle x,Xx\rangle=\langle A^{1/2}x,TA^{1/2}x\rangle`$ and
   $`A^{1/2}`$ is invertible (here $`A^{\pm1/2}`$ are given by the functional calculus with $`t^{\pm1/2}`$ on
   $`\mathrm{spec}\,A\subset(0,\infty)`$, so $`A^{-1/2}A^{1/2}=1`$ by Theorem 1.5(c)). Also
   $`A^{-1/2}BA^{-1/2}\ge0`$, since $`\langle x,A^{-1/2}BA^{-1/2}x\rangle=\langle A^{-1/2}x,BA^{-1/2}x\rangle\ge0`$.
   By Theorem 1.6, $`T=(A^{-1/2}BA^{-1/2})^{1/2}`$ is the unique
   positive solution, i.e. $`X=A\#B`$. Every solution $`X`$ of $`XA^{-1}X=B`$ is invertible,
   because $`\det(X)^2\det(A)^{-1}=\det B\neq0`$. So $`XA^{-1}X=B`$ is equivalent to
   $`XB^{-1}X=A`$ (from $`A^{-1}=X^{-1}BX^{-1}`$, invert both sides), and uniqueness gives
   $`A\#B=B\#A`$.

</details>

*Computation:* Exercises 1 and 3 numerically — [tutorial, Exercises](../tutorial/01_matrix_algebras.md#x-1).

**Next:** [Chapter 2](02_von_neumann_algebras.md) replaces $`M_n`$ by a direct sum of
matrix algebras and proves that every finite-dimensional von Neumann algebra has
this form.
