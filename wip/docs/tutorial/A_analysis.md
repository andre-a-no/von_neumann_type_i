**English** · [Español](es/A_analysis.md) · [Français](fr/A_analysis.md) · [Deutsch](de/A_analysis.md) · [中文](zh/A_analysis.md) · [日本語](ja/A_analysis.md) · [Русский](ru/A_analysis.md)

# Appendix A. Facts from analysis and probability

The chapters use a handful of facts from analysis and probability that go beyond the prerequisites
of the course but are not theorems of their own: zero sets of polynomials are null sets, the ordered
eigenvalues of a Hermitian matrix depend continuously on it, non-degenerate linear images of Gaussian vectors
have densities, the direction of a Gaussian vector is uniform on the sphere and independent of its length,
Gaussian moments are sums over pairings, and the spacings of uniform points are normalised
exponential variables. This appendix proves them, so that the chapters can cite them by number.

**Prerequisites used.** Linear algebra; real analysis (limits, continuity, compactness and the
Bolzano–Weierstrass theorem in $\mathbb R^n$, power series in one variable, the Riemann integral of a
continuous function, the change of variables formula); measure and integration (Lebesgue measure,
product measures and Fubini–Tonelli, monotone and dominated convergence, push-forward measures,
densities, independence); elementary probability (the law $N(0,s^2)$ with density
$(2\pi s^2)^{-1/2}e^{-y^2/2s^2}$ and moment generating function $\mathbb Ee^{aY}=e^{a^2s^2/2}$). For an
invertible linear map $L$ of $\mathbb R^m$ and a Borel set $E$, $\lambda_m(LE)=|\det L|\,\lambda_m(E)$, where
$\lambda_m$ is Lebesgue measure; this is the linear case of the change of variables formula.

Throughout, "measurable" means Borel measurable, $|x|$ is the Euclidean norm of $x\in\mathbb R^n$ or
$\mathbb C^n$, $\|A\|=\sup_{|x|=1}|Ax|$ the operator norm and $\|A\|_2=(\operatorname{Tr}A^*A)^{1/2}$ the
Frobenius norm.

## A.1 Measures and integrals

**Lemma A.1 (Dynkin).** Let $\mu,\nu$ be finite measures on a measurable space $(\Omega,\mathcal F)$ with
$\mu(\Omega)=\nu(\Omega)$, and let $\mathcal P\subseteq\mathcal F$ be closed under finite intersections and
generate $\mathcal F$ as a $\sigma$-algebra. If $\mu(P)=\nu(P)$ for all $P\in\mathcal P$, then $\mu=\nu$.

*Proof.* Call $\mathcal D\subseteq\mathcal F$ a *Dynkin system* if $\Omega\in\mathcal D$, $B\setminus A\in\mathcal D$
whenever $A\subseteq B$ are in $\mathcal D$, and $\bigcup_kA_k\in\mathcal D$ for every increasing sequence
$A_k\in\mathcal D$. The family $\{A\in\mathcal F:\mu(A)=\nu(A)\}$ is a Dynkin system, since
$\mu(B\setminus A)=\mu(B)-\mu(A)$ for finite measures and measures are continuous along increasing
unions. Let $d$ be the intersection of all Dynkin systems containing $\mathcal P$, itself a Dynkin system.
The family $\mathcal D_1=\{A\in d:A\cap P\in d\text{ for all }P\in\mathcal P\}$ is a Dynkin system (for
$A\subseteq A'$, $(A'\setminus A)\cap P=(A'\cap P)\setminus(A\cap P)$) and contains $\mathcal P$, so
$\mathcal D_1=d$. Then $\mathcal D_2=\{A\in d:A\cap B\in d\text{ for all }B\in d\}$ is a Dynkin system that
contains $\mathcal P$ (because $\mathcal D_1=d$), so $\mathcal D_2=d$: $d$ is closed under finite intersections.
Hence $d$ is a $\sigma$-algebra: $\Omega\setminus A\in d$, $A\cup B=\Omega\setminus\big((\Omega\setminus A)\cap(\Omega\setminus B)\big)\in d$,
and a countable union is the increasing union of finite unions. So
$\mathcal F=\sigma(\mathcal P)\subseteq d\subseteq\{\mu=\nu\}$. $\square$

**Corollary A.2.** Two finite Borel measures $\mu,\nu$ on a metric space $M$ with $\int f\,d\mu=\int f\,d\nu$
for every bounded continuous $f:M\to\mathbb R$ are equal.

*Proof.* With $f=1$, $\mu(M)=\nu(M)$. For a nonempty closed set $K$ (for $K=\emptyset$ there is nothing to prove), the continuous functions
$f_k(x)=\max(0,1-k\operatorname{dist}(x,K))$ take values in $[0,1]$ and decrease to the indicator of $K$
(for $x\notin K$, $\operatorname{dist}(x,K)>0$). By dominated convergence $\mu(K)=\nu(K)$. The closed sets are
closed under finite intersections and generate the Borel sets; apply Lemma A.1. $\square$

Lemma 5.1 uses exactly this argument.

**Lemma A.3 (Lebesgue measure on a vector space).** Let $W$ be a real vector space of dimension $m$. For a
linear bijection $\varphi:\mathbb R^m\to W$ put $\lambda_\varphi(E)=\lambda_m(\varphi^{-1}E)$ for Borel $E\subseteq W$.

(a) For another linear bijection $\psi$, $\lambda_\psi=|\det(\psi^{-1}\varphi)|\,\lambda_\varphi$. So the null sets of
$W$, and whether a random vector in $W$ has a density, do not depend on the choice of $\varphi$.

(b) For every linear bijection $L:W\to W$, $\lambda_\varphi(LE)=|\det L|\,\lambda_\varphi(E)$. If $L$ is orthogonal
for some inner product on $W$, then $|\det L|=1$ and $L$ preserves $\lambda_\varphi$.

(c) If a random vector $X$ in $W$ has the density $f$ with respect to $\lambda_\varphi$, then $LX$ has the density
$y\mapsto f(L^{-1}y)/|\det L|$.

*Proof.* (a) $\lambda_\psi(E)=\lambda_m\big((\psi^{-1}\varphi)\varphi^{-1}E\big)=|\det(\psi^{-1}\varphi)|\,\lambda_m(\varphi^{-1}E)$.
(b) $\lambda_\varphi(LE)=\lambda_m\big((\varphi^{-1}L\varphi)\varphi^{-1}E\big)=|\det(\varphi^{-1}L\varphi)|\,\lambda_\varphi(E)$, and
$\det(\varphi^{-1}L\varphi)=\det L$ does not depend on the basis. In an orthonormal basis an orthogonal $L$ has a
matrix $O$ with $O^TO=1$, so $(\det L)^2=1$. (c) By (b), the push-forward of $\lambda_\varphi$ under $L$ is
$|\det L|^{-1}\lambda_\varphi$, so for Borel $E$,
$\mathbb P(LX\in E)=\int1_E(Lx)f(x)\,d\lambda_\varphi(x)=|\det L|^{-1}\int1_E(y)f(L^{-1}y)\,d\lambda_\varphi(y)$. $\square$

*Conventions.* "Lebesgue measure" on a real vector space means any $\lambda_\varphi$; by (a) the choice never
matters for null sets and densities. On $\mathbb C^n$ and on $M_{a\times b}(\mathbb C)$ we use the real and imaginary
parts of the entries as coordinates; they are orthonormal for the real inner product
$\operatorname{Re}\operatorname{Tr}(X^*Y)$, for which every map $X\mapsto VXW$ with $V,W$ unitary is orthogonal
(it preserves $\operatorname{Tr}(X^*Y)$). On the Hermitian matrices $\operatorname{Herm}_n$ we use the coordinates
$h_{ii}$ and $\operatorname{Re}h_{ij},\operatorname{Im}h_{ij}$ ($i<j$); they are not orthonormal for the real inner
product $\operatorname{Tr}(XY)$, but by (b) the maps $H\mapsto UHU^*$ ($U$ unitary), which are orthogonal for it,
still preserve Lebesgue measure. Likewise, on the real symmetric matrices $\operatorname{Sym}_n(\mathbb R)$ we use the
coordinates $h_{ii}$ and $h_{ij}$ ($i<j$), and the maps $H\mapsto OHO^T$ ($O\in O(n)$), orthogonal for $\operatorname{Tr}(XY)$,
preserve Lebesgue measure.

**Lemma A.4 (integrals of vector-valued functions).** Let $V$ be a finite-dimensional real or complex vector
space and $F:[a,b]\to V$ continuous. Define $\int_a^bF(t)\,dt$ coordinatewise in a basis (real and imaginary
parts separately). Then:

(a) for partitions $a=t_0<\dots<t_N=b$ with tags $\tau_i\in[t_{i-1},t_i]$, the Riemann sums
$\sum_iF(\tau_i)(t_i-t_{i-1})$ converge to $\int_a^bF$ as $\max_i(t_i-t_{i-1})\to0$; in particular the integral
does not depend on the basis;

(b) $L\int_a^bF=\int_a^bLF$ for every (real-)linear map $L:V\to V'$ into a finite-dimensional space $V'$;

(c) $\|\int_a^bF\|\le\int_a^b\|F(t)\|\,dt$ for every norm on $V$;

(d) if all $F(t)$ lie in a closed convex set $K$, then $\frac1{b-a}\int_a^bF\in K$; if $K$ is moreover a cone,
$\int_a^bF\in K$;

(e) $t\mapsto\int_a^tF$ is differentiable with derivative $F(t)$; $\int_a^bF=\int_a^cF+\int_c^bF$, and
$\int_a^bF(s+c)\,ds=\int_{a+c}^{b+c}F(s)\,ds$.

*Proof.* Each coordinate of $F$ is a continuous real function, so (a) and (e) hold coordinatewise; (b) follows
from (a), since $L$ maps Riemann sums to Riemann sums and is continuous. A norm is continuous
($|\|x\|-\|y\||\le\|x-y\|\le\sum_k|x_k-y_k|\,\|e_k\|$), so $t\mapsto\|F(t)\|$ is continuous, and (c) follows from
(a) and the triangle inequality for Riemann sums. For (d) take the partitions into $N$ equal parts: the
Riemann sums divided by $b-a$ are convex combinations $\frac1N\sum_iF(\tau_i)\in K$, and $K$ is closed. $\square$

This justifies Lemma 8.4(c) (Riemann sums of CP maps) and the Cesàro means
$\sigma_\tau=\frac1\tau\int_0^\tau T_s\rho\,ds$ of Proposition 8.12: $\sigma_\tau$ is a state by (d), and
$\|T_t\sigma_\tau-\sigma_\tau\|_1\le2t/\tau$ by (b), (c) and (e).

## A.2 Zero sets of analytic functions

For a multi-index $\alpha\in\mathbb N^n$ write $|\alpha|=\alpha_1+\dots+\alpha_n$ and $y^\alpha=y_1^{\alpha_1}\cdots y_n^{\alpha_n}$.
For $a\in\mathbb R^n$ and $0<\rho\le\infty$ let $Q_\rho(a)=\{x\in\mathbb R^n:\max_i|x_i-a_i|<\rho\}$ be the open
cube ($Q_\infty(a)=\mathbb R^n$).

**Definition.** A function $f:Q_\rho(a)\to\mathbb R$ is a **power series on $Q_\rho(a)$** if
$f(x)=\sum_\alpha c_\alpha(x-a)^\alpha$ with $\sum_\alpha|c_\alpha|r^{|\alpha|}<\infty$ for every $0<r<\rho$.

Polynomials are power series on $\mathbb R^n$ (finitely many $c_\alpha\neq0$); in several variables the course applies
Theorem A.6 only to polynomials. (A function that is only locally a power series on a connected open set would need, in
addition, the identity theorem in several variables, which the course does not use.) On the closed cube of radius
$r<\rho$ the terms are bounded by $|c_\alpha|r^{|\alpha|}$, so the series converges absolutely and uniformly
there (Weierstrass M-test): its sum does not depend on the order of the terms, and $f$ is continuous.

**Lemma A.5 (one variable).** Let $g(t)=\sum_kc_k(t-a)^k$ be a power series on $I=(a-\rho,a+\rho)$.

(a) For $b\in I$ and $|t-b|<\rho_b=\rho-|b-a|$, $g(t)=\sum_jd_j(t-b)^j$ with
$d_j=\sum_{k\ge j}\binom kjc_k(b-a)^{k-j}$, a power series on $(b-\rho_b,b+\rho_b)$.

(b) If some $d_j\neq0$, then $g$ has no zero in a punctured neighbourhood of $b$.

(c) If the zeros of $g$ accumulate at a point of $I$, then $g=0$ on $I$. If $g\not\equiv0$ on $I$, its zero set is
countable, hence a Lebesgue null set. If $g\equiv0$ on $I$, then all $c_k=0$.

*Proof.* (a) For $0\le s<\rho_b$,
$\sum_k\sum_{j\le k}|c_k|\binom kjs^j|b-a|^{k-j}=\sum_k|c_k|(s+|b-a|)^k<\infty$, so the double series obtained by
expanding $(t-a)^k=\big((t-b)+(b-a)\big)^k$ binomially converges absolutely for $|t-b|\le s$ and may be
summed in either order; summing over $k$ first gives the $d_j$, with $\sum_j|d_j|s^j<\infty$.
(b) Let $j_0$ be the least $j$ with $d_j\ne0$. Then $g(t)=(t-b)^{j_0}h(t)$ with
$h(t)=\sum_{j\ge j_0}d_j(t-b)^{j-j_0}$, which is continuous near $b$ (M-test) with $h(b)=d_{j_0}\neq0$; so
$h\neq0$ near $b$. (c) Let $A$ be the set of $b\in I$ near which $g$ vanishes identically; $A$ is open. If the
zeros accumulate at $b\in I$, then by (b) all $d_j$ vanish and $b\in A$ by (a). If $b\in I$ is a limit of
points of $A$, then the zeros accumulate at $b$, so $b\in A$: $A$ is closed in $I$. As $I$ is connected, $A$ is
empty or $I$. If $g\not\equiv0$, then $A=\emptyset$, so no point of $I$ is an accumulation point of zeros; every
compact interval $J\subset I$ then contains finitely many zeros (an infinite subset of $J$ has an
accumulation point in $J$ by Bolzano–Weierstrass), and $I$ is a countable union of such $J$. If $g\equiv0$,
(b) at $b=a$ (where $d_j=c_j$) shows that all $c_k=0$. $\square$

**Theorem A.6 (zero sets are null).** Let $f$ be a power series on a cube $Q\subseteq\mathbb R^n$. If $f$ vanishes
on a subset of $Q$ of positive Lebesgue measure, then $f=0$ on $Q$. Equivalently: if $f\not\equiv0$ on $Q$,
then $\{x\in Q:f(x)=0\}$ is a Lebesgue null set.

*Proof.* Induction on $n$; the case $n=1$ is Lemma A.5(c). Let $n\ge2$, write $x=(x',t)$ with
$x'\in\mathbb R^{n-1}$, $Q=Q'\times I$ and $a=(a',a_n)$, and group the terms by the power of $t-a_n$:
$$
f(x',t)=\sum_kh_k(x')(t-a_n)^k,\qquad h_k(x')=\sum_{\alpha'}c_{(\alpha',k)}(x'-a')^{\alpha'} .
$$
This is allowed: for $x'\in Q'$ and $\max_i|x'_i-a'_i|\le r<\rho$,
$\sum_k\sum_{\alpha'}|c_{(\alpha',k)}||x'-a'|^{\alpha'}r^k\le\sum_\alpha|c_\alpha|r^{|\alpha|}<\infty$ (with
$|y|^{\alpha'}=\prod_i|y_i|^{\alpha'_i}$). The same bound shows that each $h_k$ is a power series on $Q'$ and that, for
fixed $x'\in Q'$, $g_{x'}(t)=f(x',t)$ is a power series on $I$ with coefficients $h_k(x')$.

Let $Z=\{x\in Q:f(x)=0\}$, closed in $Q$ since $f$ is continuous, and assume $\lambda_n(Z)>0$. By Tonelli's
theorem, $\lambda_n(Z)=\int_{Q'}\lambda_1(Z_{x'})\,dx'$ with $Z_{x'}=\{t\in I:f(x',t)=0\}$, so the set
$F=\{x'\in Q':\lambda_1(Z_{x'})>0\}$ has positive $(n-1)$-dimensional measure. For $x'\in F$, Lemma A.5(c) gives
$g_{x'}\equiv0$ on $I$ and then $h_k(x')=0$ for all $k$. So each $h_k$ vanishes on $F$, and by the induction
hypothesis $h_k=0$ on $Q'$. Hence $f(x',t)=\sum_kh_k(x')(t-a_n)^k=0$ on $Q$. $\square$

**Corollary A.7.** (a) A real polynomial $p$ on $\mathbb R^n$ that is not identically zero vanishes only on a
Lebesgue null set. A polynomial with a nonzero coefficient is not identically zero.

(b) Let $W$ be a finite-dimensional real vector space and $p:W\to\mathbb R$ a polynomial in the coordinates of
some basis (then it is one in every basis), not identically zero. If a random vector $X$ in $W$ has a density,
then $\mathbb P(p(X)=0)=0$.

(c) A proper linear subspace $U\subsetneq W$ is a null set, and a random vector with a density lies outside
$U$ almost surely. In particular, if $E\neq0$ is a subspace of $\mathbb C^n$ and $v$ has a density on $\mathbb C^n$
(e.g. a Gaussian start vector), then $v$ has a nonzero component in $E$ almost surely.

(d) Let $X$ and $Y$ be independent random vectors in finite-dimensional real spaces $W_1$ and $W_2$, $Y$ with a
density, and $N\subseteq W_1\times W_2$ a Borel set whose sections $N_x=\{y:(x,y)\in N\}$ are null sets for all
$x$. Then $\mathbb P((X,Y)\in N)=0$.

*Proof.* (a) The first statement is Theorem A.6 with $Q=\mathbb R^n$. For the second, induction on $n$: write
$p(x',t)=\sum_kp_k(x')t^k$; some $p_k$ has a nonzero coefficient, so $p_k(x'_0)\neq0$ for some $x'_0$ by
induction (for $n=1$ the $p_k$ are constants and there is no $x'$); then $t\mapsto p(x'_0,t)$ is a nonzero polynomial in one variable, which has finitely many roots.
(b) A change of basis is linear, so $p$ remains a polynomial, and $\{p=0\}$ is null by (a) and Lemma A.3(a);
then $\mathbb P(p(X)=0)=\int_{\{p=0\}}f\,d\lambda=0$. (c) Choose a nonzero linear functional $\ell$ vanishing on $U$
and apply (b) to $p=\ell$. For the example, $v\mapsto P_Ev$ ($P_E$ the orthogonal projection) is a nonzero real-linear
map, and $\{P_Ev=0\}$ is its kernel, a proper subspace. (d) The law of $(X,Y)$ is the product of the laws, so by
Tonelli's theorem $\mathbb P((X,Y)\in N)=\int\mathbb P(Y\in N_x)\,d\mathbb P_X(x)=0$. $\square$

The course applies Corollary A.7 to $\det Z$ for a Ginibre matrix $Z$ ($|\det Z|^2$ is a polynomial in the real
and imaginary parts of the entries, equal to $1$ at $Z=1$; Theorem 5.3, Section 1.6), to the Gram determinant of
Proposition 4.7(a), to the rank of Gaussian matrices (Chapter 12, Exercise 4), and to start vectors of the
Lanczos method (Chapter 9). Part (d) is the step "as $Y$ is independent of $K$, Tonelli's theorem gives
probability $0$" in Proposition 4.7(b).

```python
import math
import numpy as np

rng = np.random.default_rng(0)


def mc(samples, exact, label, kind="exact"):
    """Monte Carlo mean with standard error; the reference value must lie within 5 standard errors."""
    samples = np.asarray(samples, dtype=float)
    mean, err = samples.mean(), samples.std() / math.sqrt(samples.size)
    print(f"{label}: {mean:.4f} +- {err:.4f}   ({kind} {exact:.4f})")
    assert abs(mean - exact) <= 5 * err + 1e-12, label


# Corollary A.7: det Z = 0 has probability 0 for a real Gaussian 3 x 3 matrix, so by continuity of the measure
# P(|det Z| < eps) -> 0 as eps -> 0. Numerical observation (not a theorem of this appendix): the fractions
# below shrink by about a factor 10 per factor 10 in eps.
d = np.abs(np.linalg.det(rng.standard_normal((200_000, 3, 3))))
fractions = [(d < eps).mean() for eps in (1e-1, 1e-2, 1e-3)]
print("P(|det Z| < eps) for eps = 0.1, 0.01, 0.001:", fractions)
assert fractions[1] < 0.2 * fractions[0] and fractions[2] < 0.2 * fractions[1]
```

## A.3 Eigenvalues of Hermitian matrices

For $A\in\operatorname{Herm}_n$ let $\lambda_1(A)\le\dots\le\lambda_n(A)$ be its eigenvalues with multiplicity, and
$v_1,\dots,v_n$ an orthonormal basis of eigenvectors (Theorem 1.3, Corollary 1.4). For a unit vector $x$,
$\langle x,Ax\rangle=\sum_i\lambda_i|\langle v_i,x\rangle|^2$ with $\sum_i|\langle v_i,x\rangle|^2=1$, so
$\lambda_1(A)\le\langle x,Ax\rangle\le\lambda_n(A)$; and $\|A\|=\max_i|\lambda_i(A)|$ (Corollary 1.4).

**Theorem A.8 (Courant–Fischer).** For $A\in\operatorname{Herm}_n$ and $1\le j\le n$,
$$
\lambda_j(A)=\min_{\dim U=j}\ \max_{x\in U,\,|x|=1}\langle x,Ax\rangle=\max_{\dim U=n-j+1}\ \min_{x\in U,\,|x|=1}\langle x,Ax\rangle ,
$$
over the subspaces $U\subseteq\mathbb C^n$ of the given dimension. (The inner maxima and minima exist because the
unit sphere of $U$ is compact; the outer ones are attained.)

*Proof.* Let $\dim U=j$ and $W=\operatorname{span}(v_j,\dots,v_n)$. As $\dim U+\dim W=n+1>n$, there is a unit vector
$x\in U\cap W$, and $\langle x,Ax\rangle=\sum_{i\ge j}\lambda_i|\langle v_i,x\rangle|^2\ge\lambda_j$. So the maximum over $U$
is $\ge\lambda_j$. For $U=\operatorname{span}(v_1,\dots,v_j)$, $\langle x,Ax\rangle=\sum_{i\le j}\lambda_i|\langle v_i,x\rangle|^2\le\lambda_j$,
with equality at $x=v_j$; so the minimum over $U$ is attained and equals $\lambda_j$. For the second formula
apply the first to $-A$, whose eigenvalues are $\lambda_j(-A)=-\lambda_{n+1-j}(A)$:
$-\lambda_j(A)=\lambda_{n+1-j}(-A)=\min_{\dim U=n+1-j}\max_x(-\langle x,Ax\rangle)=-\max_{\dim U=n+1-j}\min_x\langle x,Ax\rangle$. $\square$

**Corollary A.9 (compressions).** Let $Q\in M_{n\times m}(\mathbb C)$ with $Q^*Q=1_m$ and $A\in\operatorname{Herm}_n$. Then
$\lambda_j(Q^*AQ)=\min\max\langle x,Ax\rangle$ over the $j$-dimensional $U\subseteq\operatorname{ran}Q$ and unit $x\in U$, and
$$
\lambda_j(A)\le\lambda_j(Q^*AQ)\le\lambda_{j+n-m}(A)\qquad(1\le j\le m).
$$

*Proof.* $U'\mapsto QU'$ is a bijection from the $j$-dimensional subspaces of $\mathbb C^m$ onto those of
$\operatorname{ran}Q$, $|Qy|=|y|$ and $\langle y,Q^*AQy\rangle=\langle Qy,AQy\rangle$; so Theorem A.8 for $Q^*AQ$ gives the
formula. Minimising over the larger family of all $j$-dimensional $U\subseteq\mathbb C^n$ gives $\lambda_j(A)$, which is
therefore smaller. Likewise the max–min formula of $\lambda_j(Q^*AQ)$ runs over $(m-j+1)$-dimensional
$U\subseteq\operatorname{ran}Q$, and over all $U\subseteq\mathbb C^n$ of that dimension it gives the larger number
$\lambda_{n-(m-j+1)+1}(A)=\lambda_{j+n-m}(A)$. $\square$

**Theorem A.10 (Weyl).** For $A,E\in\operatorname{Herm}_n$ and all $j$,
$$
\lambda_j(A)+\lambda_1(E)\le\lambda_j(A+E)\le\lambda_j(A)+\lambda_n(E),
\qquad\text{hence}\qquad
|\lambda_j(A+E)-\lambda_j(A)|\le\|E\|\le\|E\|_2 .
$$

*Proof.* For unit $x$, $\langle x,(A+E)x\rangle\le\langle x,Ax\rangle+\lambda_n(E)$. Take the maximum over the unit vectors of a
$j$-dimensional $U$ and then the minimum over $U$ (Theorem A.8): $\lambda_j(A+E)\le\lambda_j(A)+\lambda_n(E)$. Applied to
$A+E$ and $-E$ this gives $\lambda_j(A)\le\lambda_j(A+E)+\lambda_n(-E)=\lambda_j(A+E)-\lambda_1(E)$. Hence
$|\lambda_j(A+E)-\lambda_j(A)|\le\max(|\lambda_1(E)|,|\lambda_n(E)|)=\|E\|$, and $\|E\|^2=\max_i\lambda_i(E)^2\le\sum_i\lambda_i(E)^2=\|E\|_2^2$. $\square$

**Corollary A.11.** Each $\lambda_j:\operatorname{Herm}_n\to\mathbb R$ is $1$-Lipschitz for $\|\cdot\|$, hence continuous and
Borel. For a random Hermitian matrix $H$, $(\lambda_1(H),\dots,\lambda_n(H))$ is a random vector, and so are its spacings and
spacing ratios (on the Borel set where their denominators are nonzero); the eigenvalue vectors of independent random Hermitian matrices are independent. The matrices with
$n$ distinct eigenvalues form an open set.

*Proof.* Theorem A.10 with $E=B-A$. Functions of independent random variables are independent. If
$\lambda_j(A)<\lambda_{j+1}(A)$ for all $j$, the strict inequalities persist for $\|B-A\|<\frac12\min_j(\lambda_{j+1}(A)-\lambda_j(A))$. $\square$

**Lemma A.12 (Lagrange formula).** Let $A\in\operatorname{Herm}_n$ have the distinct eigenvalues $\mu_1<\dots<\mu_r$ and the
spectral projections $P_1,\dots,P_r$ (orthogonal projections onto the eigenspaces, so $A=\sum_i\mu_iP_i$,
$P_iP_l=\delta_{il}P_i$, $\sum_iP_i=1$). Then $p(A)=\sum_ip(\mu_i)P_i$ for every polynomial $p$, and
$$
P_i=\prod_{l\neq i}\frac{A-\mu_l}{\mu_i-\mu_l} .
$$

*Proof.* $A^k=\sum_i\mu_i^kP_i$ by induction on $k$ (using $P_iP_l=\delta_{il}P_i$), and $p(A)=\sum_ip(\mu_i)P_i$ by linearity.
The polynomial $\ell_i(t)=\prod_{l\neq i}(t-\mu_l)/(\mu_i-\mu_l)$ has $\ell_i(\mu_l)=\delta_{il}$, so $\ell_i(A)=P_i$. $\square$

**Proposition A.13 (continuity of eigenvalues and spectral projections).** Let $S$ be a subset of a finite-dimensional
real vector space, $A:S\to\operatorname{Herm}_n$ continuous, and suppose that $A(x)$ has exactly $r$ distinct eigenvalues
$\mu_1(x)<\dots<\mu_r(x)$ for every $x\in S$. Then their multiplicities are locally constant on $S$, and the $\mu_i$ and
the spectral projections $P_i(x)$ of $A(x)$ are continuous on $S$. In particular, on the open set of Hermitian
matrices with simple spectrum, the eigenvalues and the rank-one spectral projections are continuous.

*Proof.* If $r=1$, $A(x)=\mu_1(x)1$ with $\mu_1=\operatorname{Tr}A/n$ and $P_1=1$. Let $r\ge2$, fix $x_0\in S$ and let
$d=\min_i(\mu_{i+1}(x_0)-\mu_i(x_0))>0$. Let $J_i$ be the set of indices $j$ with $\lambda_j(A(x_0))=\mu_i(x_0)$, a block of
consecutive integers. By continuity of $A$ there is a neighbourhood $N$ of $x_0$ in $S$ with $\|A(x)-A(x_0)\|<d/2$ for
$x\in N$; by Theorem A.10, $\lambda_j(A(x))$ lies in the interval $I_i=(\mu_i(x_0)-d/2,\mu_i(x_0)+d/2)$ for $j\in J_i$, and
these $r$ intervals are disjoint. So $A(x)$ has at least one distinct eigenvalue in each $I_i$, and, having exactly $r$,
exactly one: $\lambda_j(A(x))$ is the same for all $j\in J_i$. Hence on $N$ the multiplicities are the constants $|J_i|$ and
$\mu_i(x)=\lambda_{\max J_i}(A(x))$ is continuous (Corollary A.11). By Lemma A.12,
$P_i(x)=\prod_{l\neq i}(A(x)-\mu_l(x))/(\mu_i(x)-\mu_l(x))$, continuous on $N$ because the denominators stay away from $0$.
The last statement is the case $r=n$, $S$ the open set of Corollary A.11. $\square$

In Proposition 4.7, $S$ is the set where the Gram determinant $\delta$ is nonzero, on which $K$ has exactly $N$ distinct
eigenvalues; so $K\mapsto\Pi_i(K)$ is continuous there and the event of Proposition 4.7(b) is Borel. Theorem 5.15
uses Corollary A.11, and Theorem 9.14 uses Corollary A.9.

```python
def herm(n, batch):
    G = rng.standard_normal((batch, n, n)) + 1j * rng.standard_normal((batch, n, n))
    return (G + G.conj().transpose(0, 2, 1)) / 2


n = 6
A, E = herm(n, 2000), 0.1 * herm(n, 2000)
lam_A, lam_AE, lam_E = np.linalg.eigvalsh(A), np.linalg.eigvalsh(A + E), np.linalg.eigvalsh(E)   # ascending
normE = np.linalg.norm(E, ord=2, axis=(1, 2))                                    # operator norm
assert (np.abs(lam_AE - lam_A).max(axis=1) <= normE + 1e-12).all()               # Theorem A.10
assert (lam_AE >= lam_A + lam_E[:, :1] - 1e-12).all() and (lam_AE <= lam_A + lam_E[:, -1:] + 1e-12).all()

a0, j = A[0], 3                          # Theorem A.8: every 3-dim subspace has max <x, Ax> >= lambda_3
for _ in range(200):
    Q, _ = np.linalg.qr(rng.standard_normal((n, j)) + 1j * rng.standard_normal((n, j)))
    assert np.linalg.eigvalsh(Q.conj().T @ a0 @ Q)[-1] >= lam_A[0, j - 1] - 1e-12
m = 4                                    # Corollary A.9: interlacing for a compression to 4 dimensions
Q, _ = np.linalg.qr(rng.standard_normal((n, m)) + 1j * rng.standard_normal((n, m)))
theta = np.linalg.eigvalsh(Q.conj().T @ a0 @ Q)
assert (theta >= lam_A[0, :m] - 1e-12).all() and (theta <= lam_A[0, n - m:] + 1e-12).all()

U, _ = np.linalg.qr(rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)))
mu, mult = np.array([-1.0, 0.5, 2.0]), [2, 3, 1]          # Lemma A.12 with multiplicities 2, 3, 1
spec = np.repeat(mu, mult)
A1 = U @ np.diag(spec) @ U.conj().T
for i, m_i in enumerate(mult):
    P = np.eye(n, dtype=complex)
    for l in range(len(mu)):
        if l != i:
            P = P @ (A1 - mu[l] * np.eye(n)) / (mu[i] - mu[l])
    Vi = U[:, spec == mu[i]]
    assert np.abs(P - Vi @ Vi.conj().T).max() < 1e-10 and abs(np.trace(P).real - m_i) < 1e-10
```

## A.4 Gaussian vectors

A **standard Gaussian vector** in $\mathbb R^m$ is $g=(g_1,\dots,g_m)$ with independent $N(0,1)$ coordinates. Its law
is the product of the laws of the coordinates, so by Tonelli's theorem it has the density
$\varphi_m(x)=(2\pi)^{-m/2}e^{-|x|^2/2}$.

**Proposition A.14 (orthogonal invariance).** (a) $Og$ has the law of $g$ for every orthogonal $O\in O(m)$.

(b) Let $W$ be a real inner product space of dimension $m$. For an orthonormal basis $(w_i)$ the law of
$\sum_ig_iw_i$ does not depend on the basis; call such a vector a **standard Gaussian vector of $W$**. Its coordinates
in every orthonormal basis are independent $N(0,1)$, its law is invariant under the orthogonal maps of $W$, and if
$W=W_1\oplus W_2$ orthogonally, its components are independent standard Gaussian vectors of $W_1$ and $W_2$.

(c) For $v\in\mathbb R^m$, $\langle v,g\rangle$ has the law $N(0,|v|^2)$; in particular $\sum_ia_ig_i\sim N(0,\sum_ia_i^2)$.

*Proof.* (a) By Lemma A.3(c), $Og$ has the density $\varphi_m(O^{-1}y)/|\det O|=\varphi_m(y)$. (b) If $w'_i=\sum_jO_{ji}w_j$
is another orthonormal basis, $O$ is orthogonal and $\sum_ig_iw'_i=\sum_j(Og)_jw_j$, which has the law of
$\sum_jg_jw_j$ by (a). The coordinates of $X=\sum_ig_iw_i$ in another orthonormal basis $(w'_i)$ are $O^Tg$, again
standard by (a). An orthogonal map $L$ maps $(w_i)$ to an orthonormal basis, and $LX=\sum_ig_iLw_i$. For the last
claim take $(w_i)$ adapted to $W_1\oplus W_2$. (c) For $v=0$ both sides are the point mass at $0$. For $v\neq0$
complete $v/|v|$ to an orthonormal basis; by (b) the coordinate $\langle v/|v|,g\rangle$ is $N(0,1)$, so
$\langle v,g\rangle=|v|\,\langle v/|v|,g\rangle$ has the density $(2\pi|v|^2)^{-1/2}e^{-y^2/2|v|^2}$ (Lemma A.3(c) with $m=1$). $\square$

**Definition.** A **centred Gaussian vector** in a finite-dimensional real vector space $W$ is $X=Ag$ for a linear
map $A:\mathbb R^M\to W$ and a standard Gaussian vector $g$ of $\mathbb R^M$. In coordinates ($W=\mathbb R^m$) its
**covariance** is $C=\mathbb E[XX^T]=AA^T$. It is **non-degenerate** if $A$ is onto, equivalently $C$ is positive definite
($\langle x,Cx\rangle=|A^Tx|^2$).

Every random vector of the course that is built linearly from independent Gaussian entries — complex Gaussian
vectors, Ginibre matrices, GOE and GUE matrices, their images under twirls — is a centred Gaussian vector.

**Theorem A.15 (linear images and densities).** Let $X=Ag$ be a centred Gaussian vector in $W$.

(a) For every linear $L:W\to W'$, $LX$ is a centred Gaussian vector; it is non-degenerate if $X$ is and $L$ is onto.

(b) If $X$ is non-degenerate, it has a density with respect to Lebesgue measure on $W$; in coordinates
$W=\mathbb R^m$,
$$
f_X(x)=(2\pi)^{-m/2}(\det C)^{-1/2}\,e^{-\langle x,C^{-1}x\rangle/2}.
$$
If $X$ is degenerate, it lies in the proper subspace $\operatorname{ran}A$, a null set (Corollary A.7(c)), and has no density.

(c) The law of $X$ depends only on its covariance.

*Proof.* (a) $LX=(LA)g$, and $LA$ is onto if $A$ and $L$ are. (b) Let $K=\ker A$, so
$\mathbb R^M=K^\perp\oplus K$. By Proposition A.14(b), $g=g_1+g_2$ with $g_1$ a standard Gaussian vector of
$K^\perp$, and $X=Ag_1=Bg_1$ with $B=A|_{K^\perp}$, which is injective and, as $A$ is onto, bijective
$K^\perp\to W$; so $\dim K^\perp=m$. In an orthonormal basis of $K^\perp$ completed by one of $K$, the matrix of $A$ is
$(B\ 0)$, so $C=AA^T=BB^T$. The coordinates of $g_1$ have the density $\varphi_m$, and by Lemma A.3(c) $X$ has
the density $\varphi_m(B^{-1}x)/|\det B|$, with $|B^{-1}x|^2=\langle x,(BB^T)^{-1}x\rangle$ and $|\det B|=(\det C)^{1/2}$. If
$X$ is degenerate and had a density $f$, then $1=\mathbb P(X\in\operatorname{ran}A)=\int_{\operatorname{ran}A}f\,d\lambda=0$.
(c) $\ker C=\ker A^T$ (if $AA^Tx=0$ then $|A^Tx|^2=0$), so $\operatorname{ran}C=(\ker C)^\perp=(\ker A^T)^\perp=\operatorname{ran}A$
contains $X$. Let $d=\operatorname{rank}C$ (if $d=0$, then $X=0$), choose an orthonormal basis $p_1,\dots,p_d$ of
$\operatorname{ran}C$ (it depends on $C$ only) and let $P$ be the $d\times m$ matrix with rows $p_k^T$, so that $PP^T=1_d$ and
$P^TP$ is the orthogonal projection onto $\operatorname{ran}C$. Then $PX=(PA)g$ with $PA$ onto $\mathbb R^d$, a non-degenerate
centred Gaussian vector with covariance $PCP^T$; by (b) its law depends only on $C$, and so does the law of
$X=P^T(PX)$. $\square$

For instance, in Chapter 4 the twirl $P$ is a linear projection of $M_n(\mathbb C)$ onto $\mathcal A'$, so
$Z\mapsto P(Z)$ and $Z\mapsto P(Z+Z^*)$ are real-linear surjections onto $\mathcal A'$ and onto $\mathcal A'_h$ (for
$H\in\mathcal A'_h$, $P(\frac H2+\frac H2)=H$); for a Ginibre matrix $Z$ (non-degenerate by Proposition A.16), $Y=P(Z)$
and $K=P(Z+Z^*)$ have densities on these spaces by Theorem A.15. A Gaussian combination $\sum_k\xi_kB_k$ of a basis
$(B_k)$ of $\mathcal A'$ is the image of the standard Gaussian vector $\xi$ under a linear bijection.

**Proposition A.16 (complex Gaussian vectors, Ginibre matrices).** (a) Let $z\in\mathbb C^n$ have independent
coordinates whose real and imaginary parts are independent $N(0,\frac12)$ (a **standard complex Gaussian vector**).
Then $\sqrt2z$ is the standard Gaussian vector of $\mathbb C^n\cong\mathbb R^{2n}$ with the inner product
$\operatorname{Re}\langle u,w\rangle$; $z$ has the density $\pi^{-n}e^{-|z|^2}$; $Uz$ has the law of $z$ for every
$U\in U(n)$, so the coordinates $\langle w_i,z\rangle$ in any orthonormal basis are again independent standard complex
Gaussians; and $\mathbb E[z_i\bar z_j]=\delta_{ij}$, $\mathbb E[z_iz_j]=0$.

(b) An $a\times b$ Ginibre matrix $G$ (independent standard complex Gaussian entries) has the density
$\pi^{-ab}e^{-\|G\|_2^2}$ on $M_{a\times b}(\mathbb C)$, and $VGW$ has the law of $G$ for all $V\in U(a)$, $W\in U(b)$.

(c) An $a\times b$ real matrix $G$ with independent $N(0,1)$ entries has the density $(2\pi)^{-ab/2}e^{-\|G\|_2^2/2}$, and
$OGO'$ has the law of $G$ for all $O\in O(a)$, $O'\in O(b)$.

*Proof.* (a) The $2n$ numbers $\sqrt2\operatorname{Re}z_i,\sqrt2\operatorname{Im}z_i$ are independent $N(0,1)$, and they are the
coordinates of $\sqrt2z$ in the orthonormal basis $(e_i,ie_i)$. By Lemma A.3(c) with $L=2^{-1/2}$, $z$ has the density
$(2\pi)^{-n}e^{-|\sqrt2z|^2/2}2^{n}=\pi^{-n}e^{-|z|^2}$. A unitary $U$ is real-linear and preserves
$\operatorname{Re}\langle u,w\rangle$, hence orthogonal; apply Proposition A.14(b). The coordinates in the basis $(w_i)$ are
$W^*z$ with $W$ the unitary of columns $w_i$. Finally, with $z_j=x_j+iy_j$,
$\mathbb E[z_i\bar z_j]=\mathbb E[x_ix_j+y_iy_j]+i\,\mathbb E[y_ix_j-x_iy_j]=\delta_{ij}$ and $\mathbb E[z_iz_j]=\mathbb E[x_ix_j-y_iy_j]+i\,\mathbb E[x_iy_j+y_ix_j]=0$.
(b) $M_{a\times b}(\mathbb C)\cong\mathbb C^{ab}$ with $|G|=\|G\|_2$, and $X\mapsto VXW$ is unitary for
$\operatorname{Tr}(X^*Y)$ (conventions after Lemma A.3); apply (a). (c) The same with $\mathbb R^{ab}$ and Proposition A.14. $\square$

**Proposition A.17 (GOE and GUE).** Let $H=(G+G^T)/2$ with $G$ a real $n\times n$ matrix with independent $N(0,1)$ entries
(GOE), or $H=(G+G^*)/\sqrt2$ with $G$ a Ginibre matrix (GUE). Then the entries $h_{ij}$, $i\le j$, are independent,
$h_{ii}\sim N(0,1)$, and for $i<j$, $h_{ij}\sim N(0,\frac12)$ (GOE), resp. $\operatorname{Re}h_{ij},\operatorname{Im}h_{ij}$ are
independent $N(0,\frac12)$ (GUE). The density of $H$ in the coordinates fixed after Lemma A.3 is
$$
(2\pi)^{-n/2}\pi^{-n(n-1)/4}e^{-\operatorname{Tr}H^2/2}\ \ \text{(GOE)},\qquad
(2\pi)^{-n/2}\pi^{-n(n-1)/2}e^{-\operatorname{Tr}H^2/2}\ \ \text{(GUE)},
$$
and $OHO^T$ ($O\in O(n)$), resp. $UHU^*$ ($U\in U(n)$), has the law of $H$.

*Proof.* GOE: $G$ is the standard Gaussian vector of $M_n(\mathbb R)$ with $\operatorname{Tr}(X^TY)$, and the matrices
$e_{ii}$, $(e_{ij}\pm e_{ji})/\sqrt2$ ($i<j$) form an orthonormal basis. By Proposition A.14(b) the coordinates
$G_{ii}$, $(G_{ij}+G_{ji})/\sqrt2$, $(G_{ij}-G_{ji})/\sqrt2$ are independent $N(0,1)$. Now $h_{ii}=G_{ii}$ and
$h_{ij}=2^{-1/2}\cdot(G_{ij}+G_{ji})/\sqrt2$. GUE: $\sqrt2G$ is the standard Gaussian vector of $M_n(\mathbb C)$ with
$\operatorname{Re}\operatorname{Tr}(X^*Y)$ (Proposition A.16), with the orthonormal basis $e_{ii}$, $ie_{ii}$,
$(e_{ij}\pm e_{ji})/\sqrt2$, $i(e_{ij}\pm e_{ji})/\sqrt2$; the coordinates of $\sqrt2G$ are $\sqrt2\operatorname{Re}G_{ii}$,
$\sqrt2\operatorname{Im}G_{ii}$, $\operatorname{Re}(G_{ij}\pm G_{ji})$, $\operatorname{Im}(G_{ij}\pm G_{ji})$, independent $N(0,1)$. Now
$h_{ii}=\sqrt2\operatorname{Re}G_{ii}$, $\operatorname{Re}h_{ij}=2^{-1/2}\operatorname{Re}(G_{ij}+G_{ji})$ and
$\operatorname{Im}h_{ij}=2^{-1/2}\operatorname{Im}(G_{ij}-G_{ji})$. In both cases the density is the product of the
one-dimensional densities, $(2\pi)^{-1/2}e^{-h_{ii}^2/2}$ and $\pi^{-1/2}e^{-t^2}$ for each real coordinate $t$ of variance
$\frac12$, and $\operatorname{Tr}H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2$. The map $H\mapsto UHU^*$ is linear, invertible and
orthogonal for $\operatorname{Tr}(XY)$, so by Lemma A.3(b),(c) $UHU^*$ has the density $f(U^*HU)=f(H)$, as $f$ depends on
$\operatorname{Tr}H^2$ only; the same for $O$. $\square$

```python
Amat = rng.standard_normal((2, 4))                       # Theorem A.15: X = A g in R^2, A of rank 2
C = Amat @ Amat.T
X = rng.standard_normal((400_000, 4)) @ Amat.T
h = 1 / 400                                              # midpoint rule for the density on [0, 1]^2
grid = (np.arange(400) + 0.5) * h
xx = np.stack(np.meshgrid(grid, grid, indexing='ij'), -1).reshape(-1, 2)
dens = np.exp(-0.5 * np.einsum('ni,ij,nj->n', xx, np.linalg.inv(C), xx)) / (2 * math.pi * math.sqrt(np.linalg.det(C)))
mc(((X >= 0) & (X <= 1)).all(axis=1), dens.sum() * h * h, "P(X in [0,1]^2) against the integral of f_X",
   kind="midpoint rule")

nn, B = 4, 100_000                                       # Proposition A.17
G = rng.standard_normal((B, nn, nn))
goe = (G + G.transpose(0, 2, 1)) / 2
Gc = (rng.standard_normal((B, nn, nn)) + 1j * rng.standard_normal((B, nn, nn))) / math.sqrt(2)
gue = (Gc + Gc.conj().transpose(0, 2, 1)) / math.sqrt(2)
mc(np.einsum('bij,bji->b', goe, goe), nn * (nn + 1) / 2, "E Tr H^2, GOE")      # n + 2 * n(n-1)/2 * 1/2
mc(np.einsum('bij,bji->b', gue, gue).real, nn ** 2, "E Tr H^2, GUE")           # n + 2 * n(n-1)/2 * 1
O, _ = np.linalg.qr(rng.standard_normal((nn, nn)))
rot = O @ goe @ O.T                                                             # has the law of the GOE
mc(rot[:, 0, 0] ** 4, 3.0, "E (O H O^T)_11^4")
mc(rot[:, 0, 1] ** 2, 0.5, "E (O H O^T)_12^2")
```

## A.5 The direction of a Gaussian vector

Let $S^{n-1}=\{x\in\mathbb R^n:|x|=1\}$, $V_n=\lambda_n(\{|x|\le1\})$, and for a Borel set $B\subseteq S^{n-1}$ let
$C(B)=\{x\in\mathbb R^n:0<|x|\le1,\ x/|x|\in B\}$, a Borel set ($x\mapsto(|x|,x/|x|)$ is continuous on $\mathbb R^n\setminus0$).

**Theorem A.18 (polar coordinates).** For Borel sets $A\subseteq(0,\infty)$ and $B\subseteq S^{n-1}$,
$$
\lambda_n\big(\{x\neq0:|x|\in A,\ x/|x|\in B\}\big)=n\,\lambda_n(C(B))\int_Ar^{n-1}\,dr .
$$
Consequently, for Borel $h:(0,\infty)\to[0,\infty]$ and $\phi:S^{n-1}\to[0,\infty]$, with $\mu(B)=\lambda_n(C(B))$,
$$
\int_{\mathbb R^n}h(|x|)\,\phi(x/|x|)\,dx=n\int_0^\infty h(r)r^{n-1}\,dr\int_{S^{n-1}}\phi\,d\mu .
$$

*Proof.* Fix $B$ and let $\nu_1(A)$ and $\nu_2(A)$ be the two sides, Borel measures on $(0,\infty)$. The dilation
$x\mapsto rx$ has determinant $r^n$ and maps $C(B)$ onto $\{0<|x|\le r,\ x/|x|\in B\}$, since $rx/|rx|=x/|x|$; so
$\nu_1((0,r])=r^n\lambda_n(C(B))=\nu_2((0,r])$. For $R>0$ the restrictions of $\nu_1,\nu_2$ to $(0,R]$ are finite
(at most $\lambda_n(\{|x|\le R\})$) and agree on the intervals $(0,r]$, $r\le R$, which are closed under intersections,
contain $(0,R]$ and generate the Borel sets of $(0,R]$ (they give $(s,r]=(0,r]\setminus(0,s]$, and every relatively open
subset of $(0,R]$ is a countable union of such intervals). By Lemma A.1, $\nu_1=\nu_2$ on $(0,R]$; let $R\to\infty$. The
integral formula holds for $h=1_A$, $\phi=1_B$ (the identity just proved; $\{0\}$ is null), hence for nonnegative simple
$h,\phi$ by linearity, and for general $h,\phi$ by monotone convergence along simple $h_k\uparrow h$, $\phi_k\uparrow\phi$. $\square$

**Corollary A.19 (radial integrals).** For Borel $h\ge0$ on $(0,\infty)$,
$\int_{\mathbb R^n}h(|x|)\,dx=s_n\int_0^\infty h(r)r^{n-1}\,dr$ with $s_n=nV_n=(2\pi)^{n/2}/I_n$,
$I_n=\int_0^\infty r^{n-1}e^{-r^2/2}dr$. Here $I_1=\sqrt{\pi/2}$, $I_2=1$, $I_{n+2}=nI_n$, so $s_1=2$, $s_2=2\pi$, $s_3=4\pi$,
$s_4=2\pi^2$; in general $s_n=2\pi^{n/2}/\Gamma(n/2)$.

*Proof.* The formula is Theorem A.18 with $\phi=1$, $\mu(S^{n-1})=V_n$. For $h(r)=e^{-r^2/2}$ the left side is
$\big(\int_{\mathbb R}e^{-t^2/2}dt\big)^n=(2\pi)^{n/2}$ (Fubini), which gives $s_n$. $I_1=\frac12\sqrt{2\pi}$, $I_2=1$, and
integration by parts gives $I_{n+2}=\big[-r^ne^{-r^2/2}\big]_0^\infty+n\int_0^\infty r^{n-1}e^{-r^2/2}dr=nI_n$; the substitution
$r=\sqrt{2u}$ gives $I_n=2^{n/2-1}\Gamma(n/2)$. $\square$

With $n=2$ and $n=3$ these are the polar and spherical coordinates of Proposition 5.13.

**Definition.** The **uniform distribution** on $S^{n-1}$ is $\sigma_{n-1}(B)=\lambda_n(C(B))/V_n$. On the unit sphere of
$\mathbb C^N\cong\mathbb R^{2N}$ it is $\sigma_{2N-1}$.

**Theorem A.20 (polar decomposition of radial laws).** Let $X$ be a random vector in $\mathbb R^n$ with a density of the
form $f(|x|)$. Then $X\neq0$ almost surely, $|X|$ has the density $s_nr^{n-1}f(r)$ on $(0,\infty)$, $X/|X|$ has the law
$\sigma_{n-1}$, and $|X|$ and $X/|X|$ are independent. Moreover $\sigma_{n-1}(OB)=\sigma_{n-1}(B)$ for every $O\in O(n)$.
In particular:

(a) for a standard Gaussian vector $g$ of $\mathbb R^n$, $g/|g|\sim\sigma_{n-1}$ independently of $|g|$;

(b) for a standard complex Gaussian vector $z\in\mathbb C^N$, $\psi=z/|z|\sim\sigma_{2N-1}$, a law invariant under $U(N)$,
independent of $|z|$; and $|z|^2$ has the density $t^{N-1}e^{-t}/(N-1)!$ on $(0,\infty)$, so $\mathbb E|z|^4=N(N+1)$.

*Proof.* By Theorem A.18 with $h=f\cdot1_A$, $\phi=1_B$:
$\mathbb P(|X|\in A,\ X/|X|\in B)=n\,\lambda_n(C(B))\int_Af(r)r^{n-1}dr$. With $A=(0,\infty)$ and $B=S^{n-1}$ this is
$1=nV_n\int_0^\infty f(r)r^{n-1}dr$ (and $\mathbb P(X=0)=0$). Taking $A=(0,\infty)$ gives $\mathbb P(X/|X|\in B)=\sigma_{n-1}(B)$,
taking $B=S^{n-1}$ gives $\mathbb P(|X|\in A)=s_n\int_Af(r)r^{n-1}dr$, and the general case is the product of the two.
For $O\in O(n)$, $OC(B)=C(OB)$ and $\lambda_n(OE)=\lambda_n(E)$ (Lemma A.3(b)). (a) $\varphi_n(x)=(2\pi)^{-n/2}e^{-|x|^2/2}$ is radial.
(b) $z=g/\sqrt2$ with $g$ standard in $\mathbb R^{2N}$ (Proposition A.16), so $z/|z|=g/|g|$ and $|z|=|g|/\sqrt2$; unitaries are
orthogonal on $\mathbb R^{2N}$. The density of $|g|$ is $s_{2N}(2\pi)^{-N}r^{2N-1}e^{-r^2/2}$ with $s_{2N}=2\pi^N/(N-1)!$, and the
substitution $t=r^2/2$ gives the density of $|z|^2$; then $\mathbb E|z|^4=\int_0^\infty t^{N+1}e^{-t}dt/(N-1)!=(N+1)N$. $\square$

Lemma 7.9 is part (b). For $N=1$, $\sigma_1$ is normalised arc length, that is, the push-forward of
$\frac1{2\pi}\lambda_1$ on $[0,2\pi)$ under $\theta\mapsto e^{i\theta}$: for Borel $B\subseteq S^1$ let
$\Theta(B)=\{\theta\in(0,2\pi):e^{i\theta}\in B\}$. The map $(r,\theta)\mapsto re^{i\theta}$ is a $C^1$ bijection of the open set
$(0,1)\times(0,2\pi)$ onto the open unit disc minus the segment $[0,1)$, with Jacobian $r$, and $C(B)$ differs from the
image of $(0,1)\times\Theta(B)$ by a subset of that segment and of the unit circle, a null set. By the change of variables
formula $\lambda_2(C(B))=\int_{\Theta(B)}\int_0^1r\,dr\,d\theta=\frac12\lambda_1(\Theta(B))$, and with $V_2=\pi$,
$\sigma_1(B)=\lambda_1(\Theta(B))/2\pi$ (the point $\theta=0$ is a null set).

**Corollary A.21 (moments of a uniform unit vector).** If $u\sim\sigma_{2N-1}$ on the unit sphere of $\mathbb C^N$, then
$\mathbb E[u_a\bar u_b]=\delta_{ab}/N$, $\mathbb E|u_1|^4=\frac2{N(N+1)}$ and $\mathbb E|u_1|^2|u_2|^2=\frac1{N(N+1)}$ ($N\ge2$). If
$u\sim\sigma_{n-1}$ on $S^{n-1}\subset\mathbb R^n$, then $\mathbb Eu_1^2=\frac1n$, $\mathbb Eu_1^4=\frac3{n(n+2)}$ and
$\mathbb Eu_1^2u_2^2=\frac1{n(n+2)}$ ($n\ge2$).

*Proof.* Take $u=z/|z|$ as in Theorem A.20(b). Since $z_a\bar z_b=|z|^2u_a\bar u_b$ and $|z|$, $u$ are independent,
$\delta_{ab}=\mathbb E[z_a\bar z_b]=\mathbb E|z|^2\,\mathbb E[u_a\bar u_b]=N\,\mathbb E[u_a\bar u_b]$ ($\mathbb E|z|^2=N$ by
Proposition A.16(a)). Likewise $\mathbb E|z_1|^4=N(N+1)\mathbb E|u_1|^4$ and
$\mathbb E|z_1|^2|z_2|^2=N(N+1)\mathbb E|u_1|^2|u_2|^2$, where $\mathbb E|z_1|^2|z_2|^2=1$ by independence and
$\mathbb E|z_1|^4=\mathbb E(x^2+y^2)^2=\frac34+2\cdot\frac14+\frac34=2$ for independent $x,y\sim N(0,\frac12)$. (For
$Y\sim N(0,s^2)$, integration by parts gives $\mathbb EY^4=3s^2\,\mathbb EY^2=3s^4$, since
$y^3\cdot ye^{-y^2/2s^2}=-s^2y^3\frac d{dy}e^{-y^2/2s^2}$.) In the real case $\mathbb E|g|^2=n$,
$\mathbb E|g|^4=s_n(2\pi)^{-n/2}I_{n+4}=I_{n+4}/I_n=n(n+2)$ (Corollary A.19), $\mathbb Eg_1^4=3$, $\mathbb Eg_1^2g_2^2=1$, and
$|g|^2$, $|g|^4$ are independent of $u=g/|g|$ (Theorem A.20(a)). $\square$

These agree with Lemma 6.14, which derives them from unitary invariance alone. Conversely, invariance determines the law:

**Proposition A.22 (uniqueness of the invariant law).** Let $u$ be a random unit vector in $\mathbb C^N$ whose law is
invariant under all unitaries. Then $u\sim\sigma_{2N-1}$. For $N=1$: every rotation-invariant Borel probability measure
on the circle is normalised arc length.

*Proof.* Let $Q$ be Haar distributed on $U(N)$ (Theorem 5.3) and independent of $u$. For a fixed unit vector $v$ choose
$W\in U(N)$ with $We_1=v$. By Lemma 5.1, $QW$ is Haar distributed, so $Qv=(QW)e_1$ has the law of $Qe_1$. In
Theorem 5.3, $Q$ is the Gram–Schmidt factor of a Ginibre matrix (Lemma 5.2), so $Qe_1=z/|z|$ for its first column $z$, which has
the law $\sigma_{2N-1}$ by Theorem A.20(b). Conditioning on $u$ (Fubini, independence), $\mathbb P(Qu\in B)=\sigma_{2N-1}(B)$
for every Borel $B$. Conditioning on $Q$ instead, and using that $Uu$ has the law of $u$ for each fixed unitary $U$,
$\mathbb P(Qu\in B)=\mathbb P(u\in B)$. A probability measure on the circle is the law of the identity map, and
$U(1)$ acts by rotations. $\square$

```python
N, B = 3, 400_000                                   # Theorem A.20 (b), Corollary A.21, in C^3
z = (rng.standard_normal((B, N)) + 1j * rng.standard_normal((B, N))) / math.sqrt(2)
R2 = (np.abs(z) ** 2).sum(axis=1)
u = z / np.sqrt(R2)[:, None]
t = 2.0
p_gamma = 1 - math.exp(-t) * sum(t ** k / math.factorial(k) for k in range(N))   # P(|z|^2 <= t) for Gamma(N)
mc(R2 <= t, p_gamma, "P(|z|^2 <= 2)")
mc(R2 ** 2, N * (N + 1), "E|z|^4")
mc(np.abs(u[:, 0]) ** 4, 2 / (N * (N + 1)), "E|u_1|^4")
mc(np.abs(u[:, 0] * u[:, 1]) ** 2, 1 / (N * (N + 1)), "E|u_1|^2 |u_2|^2")
mc((R2 <= t) * np.abs(u[:, 0]) ** 4, p_gamma * 2 / (N * (N + 1)), "E[1{|z|^2 <= 2} |u_1|^4], independence")
mc((u[:, 0] * u[:, 1].conj()).real, 0.0, "E Re u_1 conj(u_2)")

x3 = rng.uniform(-1, 1, (B, 3))                     # Corollary A.19: V_3 = s_3 / 3 = 4 pi / 3
mc((x3 ** 2).sum(axis=1) <= 1, math.pi / 6, "P(uniform point of [-1,1]^3 in the unit ball)")
g3 = rng.standard_normal((B, 3))
u3 = g3 / np.linalg.norm(g3, axis=1, keepdims=True)
mc(u3[:, 0] ** 4, 3 / 15, "E u_1^4 on S^2")
phase = np.angle(z[:, 0])                           # sigma_1 is normalised arc length
mc((phase >= 0) & (phase < 1), 1 / (2 * math.pi), "P(arg z_1 in [0, 1))")
```

## A.6 Gaussian moments: the Isserlis–Wick formula

**Lemma A.23.** Let $Y$ be a centred Gaussian vector in $\mathbb R^n$ with covariance $C=(c_{ij})$. Then
$c_{ij}=\mathbb E[Y_iY_j]$, $\mathbb Ee^{\langle t,Y\rangle}=e^{\langle t,Ct\rangle/2}$ and $\mathbb Ee^{|\langle t,Y\rangle|}<\infty$ for all
$t\in\mathbb R^n$, and all moments $\mathbb E|Y_i|^k$ are finite.

*Proof.* Write $Y=Ag$. Then $\mathbb E[Y_iY_j]=\sum_{k,l}A_{ik}A_{jl}\mathbb E[g_kg_l]=(AA^T)_{ij}$. With $a=A^Tt$,
$\langle t,Y\rangle=\sum_ka_kg_k$, and by independence
$\mathbb Ee^{\langle t,Y\rangle}=\prod_k\mathbb Ee^{a_kg_k}=\prod_ke^{a_k^2/2}=e^{|A^Tt|^2/2}=e^{\langle t,Ct\rangle/2}$. Finally
$e^{|s|}\le e^s+e^{-s}$ and $|Y_i|^k\le k!\,e^{|Y_i|}$. $\square$

A **pairing** of $\{1,\dots,2p\}$ is a partition into $p$ two-element blocks; there are $(2p-1)!!=(2p)!/(2^pp!)$ of them
(the partner of $1$ can be chosen in $2p-1$ ways, and the remaining $2p-2$ elements are then paired; induction on $p$).

**Theorem A.24 (Isserlis–Wick).** Let $Y_1,\dots,Y_m$ be coordinates of a centred Gaussian vector (repetitions
allowed). If $m$ is odd, $\mathbb E[Y_1\cdots Y_m]=0$; if $m=2p$,
$$
\mathbb E[Y_1\cdots Y_{2p}]=\sum_\pi\prod_{\{a,b\}\in\pi}\mathbb E[Y_aY_b],
$$
the sum over all pairings $\pi$ of $\{1,\dots,2p\}$. The same formula holds for complex random variables
$Z_1,\dots,Z_m$ that are complex linear combinations of the coordinates of one centred Gaussian vector (for example the
entries of a Ginibre or GUE matrix and their complex conjugates), with $\mathbb E[Z_aZ_b]$ (no conjugation) in the pair factors.

*Proof.* The vector $(Y_1,\dots,Y_m)$ is the image of the given Gaussian vector under a linear map, hence centred
Gaussian (Theorem A.15(a)); let $C$ be its covariance. Fix $t\in\mathbb R^m$ and put $S=\langle t,Y\rangle$. For $s\in\mathbb R$ the
partial sums of $\sum_k(sS)^k/k!$ are bounded by $e^{|sS|}$, which is integrable by Lemma A.23, so by dominated
convergence and Lemma A.23,
$$
\sum_{k\ge0}\frac{s^k}{k!}\mathbb E[S^k]=\mathbb Ee^{sS}=e^{s^2q/2}=\sum_{p\ge0}\frac{s^{2p}q^p}{2^pp!},\qquad q=\langle t,Ct\rangle .
$$
Both sides are power series on $\mathbb R$ in the sense of Section A.2: on the right this is clear, and on the left
$\sum_k\frac{\rho^k}{k!}|\mathbb E[S^k]|\le\mathbb E\sum_k\frac{\rho^k|S|^k}{k!}=\mathbb Ee^{\rho|S|}<\infty$ for every $\rho>0$
(monotone convergence and Lemma A.23). So their coefficients agree (Lemma A.5(c) applied to the difference). Hence $\mathbb E\langle t,Y\rangle^m/m!$ equals $0$ for odd $m$ and $\langle t,Ct\rangle^p/(2^pp!)$ for $m=2p$, for every
$t\in\mathbb R^m$. Both sides are polynomials in $t$, so their coefficients agree (Corollary A.7(a) applied to the difference).
In $\mathbb E\langle t,Y\rangle^m=\sum_{i_1,\dots,i_m}t_{i_1}\cdots t_{i_m}\mathbb E[Y_{i_1}\cdots Y_{i_m}]$ the monomial $t_1t_2\cdots t_m$
comes from the $m!$ sequences $(i_1,\dots,i_m)$ that are permutations of $(1,\dots,m)$, so its coefficient in
$\mathbb E\langle t,Y\rangle^m/m!$ is $\mathbb E[Y_1\cdots Y_m]$. For $m$ odd it is therefore $0$. For $m=2p$,
$q^p=\sum c_{i_1j_1}\cdots c_{i_pj_p}\,t_{i_1}t_{j_1}\cdots t_{i_p}t_{j_p}$ over all sequences of $p$ ordered pairs; the monomial
$t_1\cdots t_{2p}$ arises exactly when the unordered pairs $\{i_k,j_k\}$ form a pairing $\pi$, and each pairing arises from
$p!\,2^p$ sequences (order of the pairs, order within each pair), each contributing $\prod_{\{a,b\}\in\pi}c_{ab}$ as $C$ is
symmetric. Dividing by $2^pp!$ gives the formula, since $c_{ab}=\mathbb E[Y_aY_b]$.
For the complex statement write $Z_a=\sum_i\alpha_{ai}Y_i$. Both sides are multilinear in $(Z_1,\dots,Z_m)$ — the left by
linearity of $\mathbb E$, the right because each $Z_a$ occurs in exactly one factor of each term — and they agree when each
$Z_a$ is one of the $Y_i$; expand. $\square$

For example $\mathbb E|z|^4=\mathbb E[z\bar zz\bar z]=1\cdot1+0+1\cdot1=2$ for a standard complex Gaussian $z$. Proposition 1.10
and Lemma 7.10 use the case $m=4$; Theorem A.24 proves it without differentiating under the expectation.

```python
def pairings(idx):
    """All pairings of the list idx, as lists of pairs."""
    if not idx:
        yield []
        return
    for i in range(1, len(idx)):
        for rest in pairings(idx[1:i] + idx[i + 1:]):
            yield [(idx[0], idx[i])] + rest


def wick(cov, idx):                                       # sum over pairings of the positions of idx
    return sum(math.prod(cov[idx[a], idx[b]] for a, b in p) for p in pairings(list(range(len(idx)))))


assert [len(list(pairings(list(range(2 * p))))) for p in range(1, 5)] == [1, 3, 15, 105]    # (2p-1)!!
Amat = rng.standard_normal((3, 3)) / math.sqrt(3)
cov = Amat @ Amat.T
Y = rng.standard_normal((1_000_000, 3)) @ Amat.T
for idx in ((0, 1, 2, 0), (0, 0, 0, 0), (0, 1, 1, 2, 2, 0)):
    mc(np.prod(Y[:, list(idx)], axis=1), wick(cov, idx), f"E Y{idx}")
mc(np.prod(Y[:, [0, 1, 2]], axis=1), 0.0, "E Y(0, 1, 2), odd")
M = np.array([[1, 1j], [1, -1j]])                  # (z, conj z) = M (x, y) with x, y independent N(0, 1/2)
cz = M @ M.T / 2                                     # E[Z_a Z_b] = [[0, 1], [1, 0]]
assert abs(wick(cz, (0, 1, 0, 1)) - 2) < 1e-12                           # E|z|^4 = 2
```

## A.7 Exponential variables and uniform spacings

A random variable $e$ is **exponential** ($e\sim\operatorname{Exp}(1)$) if it has the density $e^{-y}$ on $(0,\infty)$;
$e/\beta$ is exponential with rate $\beta$. Let $\Delta_n=\{x\in\mathbb R^n:0<x_1<\dots<x_n<1\}$.

**Theorem A.25.** Let $e_0,\dots,e_n$ ($n\ge1$) be independent $\operatorname{Exp}(1)$ variables, $T=e_0+\dots+e_n$ and
$X_k=(e_0+\dots+e_{k-1})/T$ for $k=1,\dots,n$. Then $(X_1,\dots,X_n)$ has the density $n!$ on $\Delta_n$, $T$ has the density
$t^ne^{-t}/n!$ on $(0,\infty)$, and $(X_1,\dots,X_n)$ and $T$ are independent. In particular
$(e_0,\dots,e_n)/T$ is independent of $T$.

*Proof.* The map $\Phi(x,t)=t\,(x_1,\ x_2-x_1,\ \dots,\ x_n-x_{n-1},\ 1-x_n)$ is a bijection from the open set
$\Delta_n\times(0,\infty)$ onto $(0,\infty)^{n+1}$, with inverse $y\mapsto\big(((y_0+\dots+y_{k-1})/\sum_iy_i)_{k=1}^n,\ \sum_iy_i\big)$;
both are $C^1$. In the Jacobian matrix of $\Phi$ the column of $\partial/\partial t$ has entry sum $1$ and each column
$\partial/\partial x_j$ (entries $t$ and $-t$) has entry sum $0$. Replacing the last row by the sum of all rows (which does
not change the determinant) gives the row $(0,\dots,0,1)$, and the remaining $n\times n$ block is lower triangular with
diagonal $t$; so $|\det D\Phi|=t^n$. The vector $e=(e_0,\dots,e_n)$ has the density $e^{-\sum_iy_i}$, and
$\sum_i\Phi(x,t)_i=t$. By the change of variables formula, for Borel $E\subseteq\Delta_n\times(0,\infty)$,
$$
\mathbb P\big((X,T)\in E\big)=\mathbb P\big(e\in\Phi(E)\big)=\int_{\Phi(E)}e^{-\sum_iy_i}dy=\int_Ee^{-t}t^n\,dx\,dt .
$$
The density $e^{-t}t^n=n!\cdot\frac{t^ne^{-t}}{n!}$ is a product, and $\int_0^\infty t^ne^{-t}dt=n!$, so
$\lambda_n(\Delta_n)=1/n!$; by Fubini $\mathbb P(X\in A,T\in A')=\mathbb P(X\in A)\mathbb P(T\in A')$. Finally $(e_0,\dots,e_n)/T$ is
$(X_1,X_2-X_1,\dots,1-X_n)$, a function of $X$. $\square$

**Corollary A.26 (uniform spacings).** Let $U_1,\dots,U_n$ be independent and uniform on $(0,1)$, with order statistics
$U_{(1)}<\dots<U_{(n)}$ and spacings $s_0=U_{(1)}$, $s_k=U_{(k+1)}-U_{(k)}$ ($1\le k\le n-1$), $s_n=1-U_{(n)}$. Then:

(a) ties have probability $0$, and $(U_{(1)},\dots,U_{(n)})$ has the density $n!$ on $\Delta_n$;

(b) $(s_0,\dots,s_n)$ has the law of $(e_0,\dots,e_n)/T$ with $e_i$, $T$ as in Theorem A.25;

(c) $\mathbb Es_i=\frac1{n+1}$, $\mathbb Es_i^2=\frac2{(n+1)(n+2)}$ and $\mathbb E[s_is_j]=\frac1{(n+1)(n+2)}$ for $i\neq j$;

(d) for $i\neq j$, $s_i/s_j$ has the law of $e_i/e_j$, with $\mathbb P(e_i/e_j\le t)=\frac t{1+t}$ ($t>0$); so for
$0\le i\le n-1$ the spacing ratio $\min(s_i,s_{i+1})/\max(s_i,s_{i+1})$ has the law of $\min(e_i,e_{i+1})/\max(e_i,e_{i+1})$
found in Proposition 5.14(a), $\mathbb P(r\le t)=\frac{2t}{1+t}$ on $[0,1]$, with mean $2\ln2-1$.

For $n$ independent uniform points on an interval $[a,b]$, the affine map $u\mapsto a+(b-a)u$ multiplies all spacings
by $b-a$: the spacings divided by $b-a$ satisfy (b) and (c), and the ratios in (d) are unchanged.

*Proof.* (a) The set where two coordinates coincide is a finite union of hyperplanes, a null set (Corollary A.7(c)). Off it
there is exactly one permutation $\pi$ with $U_{\pi(1)}<\dots<U_{\pi(n)}$, so for Borel $E\subseteq\Delta_n$ the event
$\{(U_{(1)},\dots,U_{(n)})\in E\}$ is, up to a null set, the disjoint union over $\pi$ of
$\{(U_{\pi(1)},\dots,U_{\pi(n)})\in E\}$. Each has probability $\lambda_n(E)$, because $(U_1,\dots,U_n)$ has the density
$1_{(0,1)^n}$ and permuting coordinates preserves both $\lambda_n$ and $(0,1)^n$. So the probability is $n!\,\lambda_n(E)$.
(b) By (a) and Theorem A.25 the vectors $(U_{(k)})_k$ and $(X_k)_k$ have the same law, and the spacings are the same linear
function of either. (c) Let $w=(e_0,\dots,e_n)/T$, independent of $T$ (Theorem A.25), with $\mathbb ET=n+1$ and
$\mathbb ET^2=\int_0^\infty t^{n+2}e^{-t}dt/n!=(n+1)(n+2)$. From $e_i=Tw_i$: $1=\mathbb Ee_i=\mathbb ET\,\mathbb Ew_i$,
$2=\mathbb Ee_i^2=\mathbb ET^2\,\mathbb Ew_i^2$ and $1=\mathbb E[e_ie_j]=\mathbb ET^2\,\mathbb E[w_iw_j]$. (d) $s_i/s_j=w_i/w_j=e_i/e_j$, and by
Tonelli $\mathbb P(e_i\le te_j)=\int_0^\infty e^{-y}(1-e^{-ty})\,dy=1-\frac1{1+t}$. The spacing ratio is a function of
$s_i/s_{i+1}$ (namely $\min(x,1/x)$ at $x=s_i/s_{i+1}$, defined almost surely since all $e_k>0$), so it has the law of the
same function of $e_i/e_{i+1}$, which is the ratio $r$ of Proposition 5.14(a) for the independent $\operatorname{Exp}(1)$
variables $e_i,e_{i+1}$. $\square$

In particular, consider the "Poisson levels" $\lambda_k=e_0+\dots+e_{k-1}$, $k=1,\dots,n+1$ (independent exponential
spacings $e_0,\dots,e_n$), and $n$ sorted uniform points with their $n+1$ spacings $s_0,\dots,s_n$. The vector of
consecutive spacing ratios $r_i=\min(d_i,d_{i+1})/\max(d_i,d_{i+1})$, $0\le i<n$, has the same joint law for
$d=(e_0,\dots,e_n)$ and for $d=(s_0,\dots,s_n)$, because the ratios do not change when all $d_i$ are divided by $T$, and
$(e_0,\dots,e_n)/T$ has the law of $(s_0,\dots,s_n)$ (part (b)); and by Theorem A.25 the first $n$ Poisson levels divided by the $(n+1)$-st have exactly the law of $n$ sorted
uniform points.

```python
n, B = 5, 200_000                                       # Corollary A.26
levels = np.sort(rng.uniform(0, 1, (B, n)), axis=1)
s = np.diff(np.concatenate([np.zeros((B, 1)), levels, np.ones((B, 1))], axis=1), axis=1)    # s_0, ..., s_n
e = rng.exponential(1.0, (B, n + 1))
w = e / e.sum(axis=1, keepdims=True)                    # normalised exponentials
for name, sp in (("uniform spacings", s), ("normalised exponentials", w)):
    mc(sp[:, 2], 1 / (n + 1), f"{name}: E s_2")
    mc(sp[:, 2] ** 2, 2 / ((n + 1) * (n + 2)), f"{name}: E s_2^2")
    mc(sp[:, 0] * sp[:, n], 1 / ((n + 1) * (n + 2)), f"{name}: E s_0 s_n")
    for tt in (0.5, 2.0):
        mc(sp[:, 1] / sp[:, 2] <= tt, tt / (1 + tt), f"{name}: P(s_1 / s_2 <= {tt})")
r = np.minimum(s[:, 1:-2], s[:, 2:-1]) / np.maximum(s[:, 1:-2], s[:, 2:-1])    # ratios of interior spacings
mc(r.mean(axis=1), 2 * math.log(2) - 1, "<r>, uniform levels")
T = e.sum(axis=1)
mc(w[:, 0] * T ** 2, n + 2, "E[w_0 T^2] = E w_0 E T^2 (independence)")
```

**Back to the** [contents](README.md).
