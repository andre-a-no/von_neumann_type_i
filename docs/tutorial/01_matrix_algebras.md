**English** · [Español](es/01_matrix_algebras.md) · [Français](fr/01_matrix_algebras.md) · [Deutsch](de/01_matrix_algebras.md) · [中文](zh/01_matrix_algebras.md) · [日本語](ja/01_matrix_algebras.md) · [Русский](ru/01_matrix_algebras.md)

# 1. Matrix algebras and functional calculus

This chapter fixes notation and proves the facts about $M_n(\mathbb C)$ that every
later chapter uses: the spectral theorem, the functional calculus, the polar
decomposition and the norms defined by singular values. Each result is followed by
the corresponding library code and a numerical check.

**You will learn** how to build operators with a prescribed spectrum, how to apply a
function to an operator, how batches of operators are stored, and how to test a
computation against an exact result.

**Conventions.** Inner products $\langle x,y\rangle=\sum_i\bar x_iy_i$ are linear in the
second argument. $A^*$ is the conjugate transpose, so $\langle x,Ay\rangle=\langle A^*x,y\rangle$.
$e_1,\dots,e_n$ is the standard basis, $1$ the identity matrix, $\operatorname{Tr}$ the
usual trace.

## 1.1 The algebra $M_n(\mathbb C)$

$M_n(\mathbb C)$ is the algebra of complex $n\times n$ matrices with the involution
$A\mapsto A^*$ and the operator norm $\|A\|=\sup_{\|x\|=1}\|Ax\|$. An element $A$ is

| name | definition |
|---|---|
| normal | $A^*A=AA^*$ |
| self-adjoint | $A=A^*$ |
| positive, $A\ge0$ | $\langle x,Ax\rangle\ge0$ for all $x$ |
| unitary | $A^*A=AA^*=1$ |
| projection | $A=A^*=A^2$ |

We write $A\le B$ if $B-A\ge0$, and $A>0$ if $A\ge0$ and $A$ is invertible.

**Lemma 1.1 (C\*-identity).** $\|A^*\|=\|A\|$ and $\|A^*A\|=\|A\|^2$.

*Proof.* $\|Ax\|=\sup_{\|y\|=1}|\langle y,Ax\rangle|=\sup_{\|y\|=1}|\langle A^*y,x\rangle|$,
so $\|A\|=\sup_{\|x\|=\|y\|=1}|\langle A^*y,x\rangle|=\|A^*\|$. Next,
$\|Ax\|^2=\langle x,A^*Ax\rangle\le\|A^*A\|$ for $\|x\|=1$, hence $\|A\|^2\le\|A^*A\|\le\|A^*\|\,\|A\|=\|A\|^2$. $\square$

> **Physical meaning.** In quantum mechanics a system with $n$ distinguishable
> configurations is described by $\mathbb C^n$. Self-adjoint matrices are
> *observables* (energy, spin, position on a lattice); their eigenvalues are the
> possible results of a measurement. Projections are yes/no questions. Unitary
> matrices are symmetries and time evolutions, $e^{-itH}$ for a Hamiltonian $H$.
> Positive matrices of trace one are *states* (density matrices, Chapter 3).

## 1.2 The spectral theorem

**Lemma 1.2 (Schur).** For every $A\in M_n(\mathbb C)$ there is a unitary $U$ with
$U^*AU$ upper triangular.

*Proof.* Induction on $n$; $n=1$ is trivial. The characteristic polynomial of $A$
has a root $\lambda$ (fundamental theorem of algebra), so $Av=\lambda v$ for some
unit vector $v$. Complete $v$ to an orthonormal basis $v,w_2,\dots,w_n$ and let
$V$ be the unitary with these columns. Then
$V^*AV=\begin{pmatrix}\lambda&*\\0&A'\end{pmatrix}$ with $A'\in M_{n-1}(\mathbb C)$.
By induction $U'^*A'U'$ is upper triangular for a unitary $U'$, and
$U=V\,(1\oplus U')$ does the job. $\square$

**Theorem 1.3 (spectral theorem).** $A$ is normal if and only if
$A=U\,\mathrm{diag}(\lambda_1,\dots,\lambda_n)\,U^*$ for a unitary $U$. The
$\lambda_i$ are the eigenvalues of $A$ with multiplicity.

*Proof.* If $A=U\Lambda U^*$ with $\Lambda$ diagonal, then
$A^*A=U\bar\Lambda\Lambda U^*=U\Lambda\bar\Lambda U^*=AA^*$. Conversely, let $A$
be normal and $T=U^*AU$ upper triangular (Lemma 1.2); $T$ is normal. Compare the
$(1,1)$ entries: $(T^*T)_{11}=|t_{11}|^2$ and $(TT^*)_{11}=\sum_j|t_{1j}|^2$, so
$t_{1j}=0$ for $j>1$. Then $T=t_{11}\oplus T'$ with $T'$ normal and upper
triangular, and induction shows that $T$ is diagonal. The diagonal of $U^*AU$
consists of the roots of $\det(A-\lambda)=\det(T-\lambda)$. $\square$

**Corollary 1.4.** Let $A$ be normal with eigenvalues $\lambda_i$. Then

1. $\|A\|=\max_i|\lambda_i|$;
2. $A$ is self-adjoint iff all $\lambda_i\in\mathbb R$; unitary iff all $|\lambda_i|=1$;
   a projection iff all $\lambda_i\in\{0,1\}$;
3. $A\ge0$ iff $A$ is self-adjoint and all $\lambda_i\ge0$.

*Proof.* With $y=U^*x$, $\|Ax\|^2=\sum_i|\lambda_i|^2|y_i|^2$ and $\|y\|=\|x\|$,
which gives (1), with equality for $x=Ue_i$, $|\lambda_i|$ maximal. In (2), $A^*=U\bar\Lambda U^*$,
$A^*A=U|\Lambda|^2U^*$ and $A^2=U\Lambda^2U^*$, and the conditions hold iff they hold
entrywise for $\Lambda$. For (3), first: every positive $A$ is self-adjoint, since
$\langle x,(A-A^*)x\rangle=\langle x,Ax\rangle-\overline{\langle x,Ax\rangle}=0$ for all $x$,
and a matrix $D$ with $\langle x,Dx\rangle=0$ for all $x$ vanishes (polarisation:
$4\langle x,Dy\rangle=\sum_{k=0}^3 i^{-k}\langle x+i^ky,D(x+i^ky)\rangle$). Then
$\langle x,Ax\rangle=\sum_i\lambda_i|y_i|^2$, which is $\ge0$ for all $x$ iff all $\lambda_i\ge0$. $\square$

## 1.3 The functional calculus

Let $A=U\Lambda U^*$ be normal and $f$ a complex function on $\operatorname{spec}A$. Define
$$
f(A)=U\,\mathrm{diag}\big(f(\lambda_1),\dots,f(\lambda_n)\big)\,U^*.
$$

**Theorem 1.5.** (a) $f(A)$ does not depend on the choice of $U$. (b) For a
polynomial $p$, $p(A)$ is the usual polynomial in $A$. (c) $f\mapsto f(A)$ is
linear and multiplicative, $\bar f(A)=f(A)^*$, and $\operatorname{spec}f(A)=f(\operatorname{spec}A)$.
(d) If $B$ commutes with $A$, then $B$ commutes with $f(A)$. (e) If $g$ is defined on
$f(\operatorname{spec}A)$, then $g(f(A))=(g\circ f)(A)$.

*Proof.* (b) $A^k=U\Lambda^kU^*$, so $p(A)=Up(\Lambda)U^*$ by linearity. (a) The
spectrum is finite, so Lagrange interpolation gives a polynomial $p$ with $p=f$ on
$\operatorname{spec}A$. Then $U\,f(\Lambda)\,U^*=U\,p(\Lambda)\,U^*=p(A)$ for *every*
unitary $U$ diagonalising $A$, and $p(A)$ is defined without $U$. (c) holds for the
diagonal matrices $f(\Lambda)$ and is preserved by conjugation with $U$; the spectrum
of $Uf(\Lambda)U^*$ is the diagonal of $f(\Lambda)$. (d) $B$ commutes with every
polynomial in $A$, in particular with $p(A)=f(A)$. (e) $f(A)=Uf(\Lambda)U^*$ is normal and
diagonalised by the same $U$, so by (a) $g(f(A))=U\,g(f(\Lambda))\,U^*=(g\circ f)(A)$. $\square$

**Theorem 1.6 (square root).** For $A\ge0$ there is exactly one $B\ge0$ with
$B^2=A$, namely $B=\sqrt A$ defined by the functional calculus.

*Proof.* $\sqrt A\ge0$ and $(\sqrt A)^2=A$ by Theorem 1.5(c) and Corollary 1.4(3).
Let $B\ge0$ with $B^2=A$, $B=V\,\mathrm{diag}(\mu_j)\,V^*$, $\mu_j\ge0$. Then
$A=V\,\mathrm{diag}(\mu_j^2)\,V^*$, so $\operatorname{spec}A=\{\mu_j^2\}$. Let $p$ be a
polynomial with $p(t)=\sqrt t$ on $\operatorname{spec}A$. Then
$p(A)=\sqrt A$ and $p(A)=p(B^2)=V\,\mathrm{diag}(p(\mu_j^2))\,V^*=V\,\mathrm{diag}(\mu_j)\,V^*=B$. $\square$

**Corollary 1.7.** $A\ge0$ iff $A=T^*T$ for some $T$. (If $A\ge0$ take $T=\sqrt A$;
conversely $\langle x,T^*Tx\rangle=\|Tx\|^2\ge0$.)

## 1.4 Polar decomposition and singular values

For any $A$, the matrix $A^*A$ is positive (Corollary 1.7); its square root
$|A|=(A^*A)^{1/2}$ is the **absolute value** of $A$.

**Theorem 1.8 (polar decomposition).** Every $A\in M_n(\mathbb C)$ can be written
$A=W|A|$ with $W$ unitary.

*Proof.* For every $x$, $\||A|x\|^2=\langle x,|A|^2x\rangle=\langle x,A^*Ax\rangle=\|Ax\|^2$.
Hence $W_0:|A|x\mapsto Ax$ is a well-defined linear isometry from the range
$R_1$ of $|A|$ onto the range $R_2$ of $A$: if $|A|x=|A|x'$ then
$\|A(x-x')\|=\||A|(x-x')\|=0$. In particular $\dim R_1=\dim R_2$, so the orthogonal
complements also have equal dimension; let $W_1:R_1^\perp\to R_2^\perp$ be any
isometry between them (map an orthonormal basis to an orthonormal basis). Then
$W=W_0\oplus W_1$ is unitary and $W|A|x=W_0|A|x=Ax$. $\square$

The eigenvalues $s_1\ge\dots\ge s_n\ge0$ of $|A|$ are the **singular values** of
$A$. Three norms are used throughout:
$$
\|A\|=s_1,\qquad \|A\|_2=\Big(\sum_is_i^2\Big)^{1/2}=(\operatorname{Tr}A^*A)^{1/2},
\qquad\|A\|_1=\sum_is_i=\operatorname{Tr}|A|.
$$
The first identity holds because $\|Ax\|=\||A|x\|$ for every $x$, so
$\|A\|=\||A|\|=s_1$ by Corollary 1.4(1); the second because
$\operatorname{Tr}A^*A=\operatorname{Tr}|A|^2=\sum_is_i^2$.

**Theorem 1.9 (duality of the trace norm).**
$\|A\|_1=\max_{U\text{ unitary}}|\operatorname{Tr}(AU)|$, and the maximum is attained at
$U=W^*$, where $A=W|A|$. In particular $|\operatorname{Tr}A|\le\|A\|_1$.

*Proof.* $\operatorname{Tr}(AW^*)=\operatorname{Tr}(W|A|W^*)=\operatorname{Tr}|A|$. For any unitary $U$,
by cyclicity and Cauchy–Schwarz for $\langle X,Y\rangle_{HS}=\operatorname{Tr}(X^*Y)$,
$$
|\operatorname{Tr}(AU)|=|\operatorname{Tr}(|A|^{1/2}\,UW\,|A|^{1/2})|
=|\langle|A|^{1/2},\,UW|A|^{1/2}\rangle_{HS}|
\le\||A|^{1/2}\|_2\,\|UW|A|^{1/2}\|_2=\operatorname{Tr}|A|,
$$
since $\|UW X\|_2=\|X\|_2$ for unitary $UW$ and $\||A|^{1/2}\|_2^2=\operatorname{Tr}|A|$. $\square$

In particular the three quantities above are norms. $\|\cdot\|$ is an operator norm;
$\|\cdot\|_2$ is the norm of the inner product $\langle X,Y\rangle_{HS}$; and for
$\|\cdot\|_1$, Theorem 1.9 gives
$\|A+B\|_1=\max_U|\operatorname{Tr}(AU)+\operatorname{Tr}(BU)|\le\|A\|_1+\|B\|_1$, while
$\|A\|_1=0$ forces $|A|=0$ and $A=W|A|=0$. For $P\ge0$, $|P|=(P^2)^{1/2}=P$ by Theorem 1.6,
so $\|P\|_1=\operatorname{Tr}P$.

## 1.5 Operators in the library

The library works with algebras $\bigoplus_cM_{k_c}(\mathbb C)$ (Chapter 2). A single
block is the case $C=1$: `TypeIAlgebra([n], [n])` is $M_n(\mathbb C)$. (The two lists
are the nominal and the actual block sizes; for now they are equal.)

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra

torch.manual_seed(0)
M3 = TypeIAlgebra([3], [3], precision='double')       # the algebra M_3(C), complex128

# a positive operator with eigenvalues 1, 2, 4 and Haar-random eigenvectors
A = M3.operator_from_eigenvalues(lambda d: torch.tensor([[1.0, 2.0, 4.0]]),
                                 force_positive=True)
print(A.matrix.shape)          # (batch, C, k_max, k_max) = (1, 1, 3, 3)
print(A.eigenvalues()[0])      # one tensor per block, eigenvalues in ascending order
```

Every operator is stored as one tensor of shape `(batch, C, k_max, k_max)`: a batch
of operators, each a list of `C` blocks. Here `batch = C = 1`. The construction is
exactly Theorem 1.3 read backwards: $A=U\Lambda U^*$ with $U$ Haar-distributed. The
flag `force_positive=True` is checked against the eigenvalues (Corollary 1.4(3)), so
the library knows that `A` is positive without testing it again.

The functional calculus for self-adjoint operators is `apply_function` (it diagonalises each
block with `eigh`). The function $f$ must act elementwise: it receives the eigenvalues of all
sectors of one size $k$ as one tensor of shape `(batch, n·k)`, $n$ the number of such sectors,
and returns a tensor of the same shape. The common functions have their own methods. We check
Theorems 1.5 and 1.6:

```python
def close(X, Y, tol=1e-10):
    """Maximal entry of |X - Y| for two operators (or tensors) is below tol."""
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

R = A.sqrt()
assert close(R @ R, A)                                  # Theorem 1.6
assert close(A.log().expm(), A)                         # exp(log A) = A, Theorem 1.5(e)
assert close(A.apply_function(lambda w: w ** 3), A @ A @ A)   # Theorem 1.5(b)
assert close(A.power(-1), A.inverse())                  # f(t) = 1/t
spec = A.apply_function(torch.exp).eigenvalues()[0]     # spec f(A) = f(spec A)
assert torch.allclose(spec, torch.exp(torch.tensor([[1.0, 2.0, 4.0]], dtype=spec.dtype)))
```

`R @ R` does not compute anything yet: products, sums and functions build a recipe,
which is evaluated when `.matrix` or a number such as a trace is requested.

## 1.6 Norms and the trace-norm duality

A general (not normal) matrix enters through `alg.operator(tensor)`. We check the
three norms against singular values computed directly by PyTorch, the polar
decomposition (Theorem 1.8) and the duality (Theorem 1.9):

```python
G = torch.randn(1, 1, 3, 3, dtype=torch.complex128)    # Ginibre matrix: no symmetry at all
X = M3.operator(G)
s = torch.linalg.svdvals(G[0, 0])                       # singular values, descending

assert abs(X.operator_norm().item() - s[0].item()) < 1e-12
assert abs(X.frobenius_norm().item() - s.norm().item()) < 1e-12
assert abs(X.trace_norm().item() - s.sum().item()) < 1e-12

absX = X.abs()
W = X @ absX.inverse()                                  # |X| is invertible (almost surely, for Ginibre)
assert close(W.adjoint() @ W, M3.identity())            # W is unitary
assert abs((X @ W.adjoint()).trace.item() - s.sum().item()) < 1e-10    # Tr(X W*) = ||X||_1
assert abs(X.trace.item()) <= X.trace_norm().item()     # |Tr X| <= ||X||_1
```

## 1.7 Batches: many operators at once

The batch axis is the reason the library exists: random matrix theory, Monte Carlo
and multistart optimisation need thousands of operators, and one batched call is much
faster than a Python loop. As an example we check an exact moment of the **Gaussian
unitary ensemble** (GUE).

Let $G$ have independent entries with independent real and imaginary parts
$\sim N(0,\tfrac12)$, so $\mathbb E|G_{ij}|^2=1$, and let $H=(G+G^*)/\sqrt2$.

**Proposition 1.10.** $\mathbb E\,\operatorname{Tr}H^4=2n^3+n$. Hence for
$X=H/\sqrt n$, $\ \mathbb E\,\tfrac1n\operatorname{Tr}X^4=2+\tfrac1{n^2}$.

*Proof.* *Isserlis' (Wick's) formula.* Let $Y_1,\dots,Y_4$ be jointly Gaussian centred
real random variables (repetitions allowed) with $c_{ij}=\mathbb E[Y_iY_j]$. Their joint
moment generating function is $\mathbb E\,e^{\sum_it_iY_i}=e^{q(t)/2}$ with
$q(t)=\sum_{ij}c_{ij}t_it_j$ (the Gaussian moment generating function, a standard fact of
probability). Differentiating under the expectation is allowed: the moment generating function
is finite for every $t$, and for $|t_i|\le1$ all partial derivatives of $e^{\sum_it_iY_i}$ up to
order $4$ are bounded by $(1+\sum_i|Y_i|)^4e^{\sum_i|Y_i|}$, which is integrable because
$e^{\sum_i|Y_i|}\le\sum_{s\in\{\pm1\}^4}e^{\sum_is_iY_i}$; so dominated convergence applies. Apply $\partial^4/\partial t_1\cdots\partial t_4$ at $t=0$: on the left this
gives $\mathbb E[Y_1Y_2Y_3Y_4]$; on the right only the term $q^2/8$ of the exponential
series contains $t_1t_2t_3t_4$, with coefficient $c_{12}c_{34}+c_{13}c_{24}+c_{14}c_{23}$.
Hence $\mathbb E[Y_1Y_2Y_3Y_4]$ is the sum over the three pairings of products of pair
expectations. By multilinearity the same holds for complex variables that are complex
linear combinations of jointly Gaussian real ones, such as the entries of $H$.

*The moment.* The pair expectations are
$\mathbb E[H_{ij}H_{kl}]=\tfrac12\mathbb E[(G_{ij}+\bar G_{ji})(G_{kl}+\bar G_{lk})]=\delta_{il}\delta_{jk}$,
because $\mathbb E[G_{ij}G_{kl}]=0$ and $\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}$.
In $\operatorname{Tr}H^4=\sum_{ijkl}H_{ij}H_{jk}H_{kl}H_{li}$ the pairings give:
$\{12\}\{34\}$: $\delta_{ik}\delta_{ki}$, summing to $n^3$ (free indices $i=k$, $j$, $l$);
$\{14\}\{23\}$: $\delta_{jl}\delta_{jl}$, summing to $n^3$;
$\{13\}\{24\}$: $\delta_{il}\delta_{jk}\delta_{ji}\delta_{kl}$, forcing $i=j=k=l$, summing to $n$.
Dividing by $n\cdot n^2$ gives the second statement. $\square$

The two *planar* pairings give the Catalan number $2$ of the semicircle law; the
crossing pairing gives the $1/n^2$ correction.

```python
n, B = 4, 200_000
G = torch.complex(torch.randn(B, 1, n, n), torch.randn(B, 1, n, n)).to(torch.complex128) / math.sqrt(2)
M4 = TypeIAlgebra([n], [n], precision='double')
X = M4.operator((G + G.conj().transpose(-2, -1)) / math.sqrt(2 * n), is_self_adjoint=True)

X2 = X @ X
m4 = (X2 @ X2).trace.real / n                           # shape (B,): one value per sample
mean, err = m4.mean().item(), m4.std().item() / math.sqrt(B)
print(f"E tr X^4 / n = {mean:.4f} +- {err:.4f}   (exact {2 + 1 / n**2:.4f})")
assert abs(mean - (2 + 1 / n**2)) < 5 * err
```

A Monte Carlo check is a statistical test: the estimate must agree with the exact
value to within a few standard errors. This is the pattern for every check against
known results in this tutorial: an exact statement, a numerical estimate, and an
explicit tolerance.

## 1.8 What to compare with

On single matrices, `numpy.linalg` and `scipy.linalg` (`sqrtm`, `expm`, `logm`,
`polar`) are the reference implementations, and the checks above can be repeated
with them. The library adds three things that matter later: batching (Section 1.7),
operators that carry a block structure (Chapter 2), and automatic differentiation
through every operation.

## Exercises

1. Check $(fg)(A)=f(A)g(A)$ for $f=\sin$, $g=\exp$ on a random self-adjoint $A$.
2. Show that $\|A\|_1=\|A^*\|_1$ and $\|AB\|_1\le\|A\|\,\|B\|_1$, using Theorem 1.9.
3. The **geometric mean** of $A,B>0$ is $A\#B=A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}$.
   Prove that $A\#B$ is the unique positive solution $X$ of $XA^{-1}X=B$, deduce
   $A\#B=B\#A$, and check both numerically.

<details>
<summary>Solutions</summary>

1. By Theorem 1.5(c).

```python
Hs = M3.operator_from_eigenvalues(lambda d: torch.randn(1, d), force_self_adjoint=True)
assert close(Hs.apply_function(lambda w: torch.sin(w) * torch.exp(w)),
             Hs.apply_function(torch.sin) @ Hs.apply_function(torch.exp))
```

2. $|\operatorname{Tr}(A^*U)|=|\overline{\operatorname{Tr}(U^*A)}|=|\operatorname{Tr}(AU^*)|$, and
   $U\mapsto U^*$ is a bijection of the unitary group, so the maxima in Theorem 1.9
   agree. For the second claim, let $B=V|B|$. For unitary $U$,
   $|\operatorname{Tr}(ABU)|=|\operatorname{Tr}(|B|^{1/2}UAV|B|^{1/2})|
   \le\||B|^{1/2}\|_2\,\|UAV|B|^{1/2}\|_2\le\|A\|\operatorname{Tr}|B|$, using
   $\|YX\|_2\le\|Y\|\,\|X\|_2$ (apply $\|Yx\|\le\|Y\|\,\|x\|$ to each column of $X$) and
   $\|UAV\|=\|A\|$ for unitaries $U,V$ (they preserve the length of vectors).

3. Write $T=A^{-1/2}XA^{-1/2}$. Then $XA^{-1}X=B$ iff $T^2=A^{-1/2}BA^{-1/2}$, and
   $X\ge0$ iff $T\ge0$, since $\langle x,Xx\rangle=\langle A^{1/2}x,TA^{1/2}x\rangle$ and
   $A^{1/2}$ is invertible. By Theorem 1.6, $T=(A^{-1/2}BA^{-1/2})^{1/2}$ is the unique
   positive solution, i.e. $X=A\#B$. Every solution $X$ of $XA^{-1}X=B$ is invertible,
   because $\det(X)^2\det(A)^{-1}=\det B\neq0$. So $XA^{-1}X=B$ is equivalent to
   $XB^{-1}X=A$ (from $A^{-1}=X^{-1}BX^{-1}$, invert both sides), and uniqueness gives
   $A\#B=B\#A$.

```python
def geometric_mean(A, B):
    a, ai = A.sqrt(), A.power(-0.5)
    return a @ (ai @ B @ ai).sqrt() @ a

P = M3.operator_from_eigenvalues(lambda d: 0.5 + torch.rand(1, d), force_positive=True)
Q = M3.operator_from_eigenvalues(lambda d: 0.5 + torch.rand(1, d), force_positive=True)
GM = geometric_mean(P, Q)
assert close(GM, geometric_mean(Q, P), 1e-9)
assert close(GM @ P.inverse() @ GM, Q, 1e-9)
```

</details>

**Next:** [Chapter 2](02_von_neumann_algebras.md) replaces $M_n$ by a direct sum of
matrix algebras and proves that every finite-dimensional von Neumann algebra has
this form.
