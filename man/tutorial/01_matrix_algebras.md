**English** · [Español](es/01_matrix_algebras.md) · [Français](fr/01_matrix_algebras.md) · [Deutsch](de/01_matrix_algebras.md) · [中文](zh/01_matrix_algebras.md) · [日本語](ja/01_matrix_algebras.md) · [Русский](ru/01_matrix_algebras.md)

# Tutorial 1. Matrix algebras and functional calculus

This tutorial accompanies [Chapter 1](../lecture_notes/01_matrix_algebras.md) of the lecture notes, which proves the spectral theorem,
the functional calculus, the polar decomposition and the duality of the trace norm in $`M_n(\mathbb C)`$.
Here each of these results is checked with the library.

**You will learn** how to build operators with a prescribed spectrum, how to apply a
function to an operator, how batches of operators are stored, and how to test a
computation against an exact result.

<a id="s-1-5"></a>

## 1.5 Operators in the library

*Lecture notes:* [Section 1.2](../lecture_notes/01_matrix_algebras.md#s-1-2) and [Section 1.3](../lecture_notes/01_matrix_algebras.md#s-1-3): [Theorem 1.3](../lecture_notes/01_matrix_algebras.md#r-1-3), [Corollary 1.4](../lecture_notes/01_matrix_algebras.md#r-1-4), [Theorem 1.5](../lecture_notes/01_matrix_algebras.md#r-1-5), [Theorem 1.6](../lecture_notes/01_matrix_algebras.md#r-1-6), [Theorem 1.11](../lecture_notes/01_matrix_algebras.md#r-1-11), [Proposition 1.12](../lecture_notes/01_matrix_algebras.md#r-1-12), [Proposition 1.13](../lecture_notes/01_matrix_algebras.md#r-1-13).

The library works with algebras $`\bigoplus_cM_{k_c}(\mathbb C)`$ (Chapter 2). A single
block is the case $`C=1`$: `TypeIAlgebra([n], [n])` is $`M_n(\mathbb C)`$. (The two lists
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
exactly Theorem 1.3 read backwards: $`A=U\Lambda U^*`$ with $`U`$ Haar-distributed (the library
draws $`U`$ by the QR construction of Theorem 5.3). The flag `force_positive=True` is checked
against the prescribed eigenvalues (the library raises an error if one of them is below
$`-10^{-12}`$); by Corollary 1.4(3) the result is then positive (exactly so when the prescribed
eigenvalues are $`\ge0`$, as here; up to this tolerance in general), and the library records this
instead of testing it again. For the same reason `A.eigenvalues()` above returns the prescribed
spectrum, sorted, without computing anything; it is not a check. The check of
$`\mathrm{spec}\,f(A)=f(\mathrm{spec}\,A)`$ below does compute: `apply_function` diagonalises the
matrix of $`A`$, and the eigenvalues of the resulting matrix $`e^A`$ are computed from scratch.

The functional calculus for self-adjoint operators is `apply_function` (it diagonalises each
block with `eigh`). The function $`f`$ must act elementwise: it receives the eigenvalues of all
sectors of one size $`k`$ as one tensor of shape `(batch, n·k)`, $`n`$ the number of such sectors,
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

**Derivatives.** By [Theorem 1.11](../lecture_notes/01_matrix_algebras.md#r-1-11) (Daleckii–Krein), the derivative of
$`f(A)`$ in a self-adjoint direction $`E`$ is $`U(L\circ U^*EU)U^*`$, with the Loewner matrix
$`L_{ij}=f^{[1]}(\lambda_i,\lambda_j)`$ of divided differences ($`f'(\lambda_i)`$ for equal eigenvalues).
`functional.divided_differences` computes $`L`$ from the eigenvalues and the values of $`f`$ and $`f'`$, and autograd
differentiates `apply_function` with it. We check the formula at a degenerate spectrum against a central difference
quotient of the dense exponential, whose error is $`O(s^2)`$ plus rounding of order $`\varepsilon/s`$:

```python
from torch_vn_algebra import functional

lam = torch.tensor([1.0, 1.0, 2.0], dtype=torch.float64)          # a degenerate spectrum
Q = torch.linalg.qr(torch.randn(3, 3, dtype=torch.complex128))[0]  # a unitary
A0 = Q @ torch.diag(lam).to(torch.complex128) @ Q.mH
E = torch.randn(3, 3, dtype=torch.complex128)
E = (E + E.mH) / 2                                                 # a self-adjoint direction
L = functional.divided_differences(lam, torch.exp(lam), torch.exp(lam))
D = Q @ (L * (Q.mH @ E @ Q)) @ Q.mH                                # Theorem 1.11 for f = exp
s = 1e-6
fd = (torch.linalg.matrix_exp(A0 + s * E) - torch.linalg.matrix_exp(A0 - s * E)) / (2 * s)
print(f"{(D - fd).abs().max().item():.1e}")
assert (D - fd).abs().max() < 1e-8
```

For close eigenvalues $`a,b`$ at distance $`h`$ the quotient $`\hat q`$ loses accuracy by cancellation, and the
midpoint value $`\hat m=(f'(a)+f'(b))/2`$ errs by at most $`T=h^2\max|f'''|/12`$ plus rounding
([Proposition 1.12](../lecture_notes/01_matrix_algebras.md#r-1-12)). `divided_differences` takes $`\hat m`$ if and only if
$`|\hat q-\hat m|\le2R`$, with the rounding bound $`R`$ of [Proposition 1.13](../lecture_notes/01_matrix_algebras.md#r-1-13)
and the evaluation error $`\varepsilon(|\hat f_a|+|a\hat f'_a|)`$; its error is then at most $`\min(T,3R)`$. We check
this for $`f(x)=e^{x-c}`$ near $`c=10^4`$, where $`f^{[1]}(a,b)=e^{a-c}(e^h-1)/h`$ is computed accurately with `expm1`.
The last distance, $`h=0.06`$, is below the threshold $`\varepsilon^{1/3}\max(|a|,|b|)`$ of earlier versions, which took
the midpoint value there:

```python
eps = torch.finfo(torch.float64).eps
c = 1e4
f = lambda x: torch.exp(x - c)
for h in (1e-12, 1e-9, 1e-6, 1e-4, 1e-2, 0.06):
    w = torch.tensor([c, c + h], dtype=torch.float64)
    h = (w[1] - w[0]).item()                                       # the distance actually stored
    fw = f(w)                                                      # f' = f
    exact = math.expm1(h) / h                                      # f[a, b], as e^(a - c) = 1
    q = ((fw[1] - fw[0]) / h).item()                               # quotient
    m = (fw.sum() / 2).item()                                      # midpoint value
    y = functional.divided_differences(w, fw, fw)[0, 1].item()     # the library's choice
    R = eps * (fw.abs() + (w * fw).abs()).sum().item() / h + 4 * eps * abs(q)
    T = h ** 2 * math.exp(h) / 12 + 2 * eps                        # Proposition 1.12, plus rounding
    print(f"h = {h:.0e}: quotient {abs(q - exact):.0e}, midpoint {abs(m - exact):.0e}, library {abs(y - exact):.0e}")
    assert abs(y - exact) <= min(T, 3 * R) + 2 * eps              # Proposition 1.13(b)
```

*Numerical observation:* the quotient errs by up to $`5\cdot10^{-10}`$ for $`h\le10^{-6}`$ and by $`10^{-15}`$ at
$`h=0.06`$; the midpoint value is exact to rounding for $`h\le10^{-9}`$ and errs by $`3\cdot10^{-4}`$ at $`h=0.06`$. The
library has the smaller of the two errors in every row but $`h=10^{-4}`$, where it takes the midpoint value ($`8\cdot10^{-10}`$) although
the quotient is better ($`10^{-12}`$): $`R`$ allows for an argument perturbed by $`\varepsilon|a|\approx2\cdot10^{-12}`$,
which is pessimistic for `torch.exp` but needed for functions such as $`e^{-i\tau x}`$, which round $`\tau x`$. The
guarantee $`\min(T,3R)`$ holds in every row.

<a id="s-1-6"></a>

## 1.6 Norms and the trace-norm duality

*Lecture notes:* [Section 1.4](../lecture_notes/01_matrix_algebras.md#s-1-4): [Theorem 1.8](../lecture_notes/01_matrix_algebras.md#r-1-8), the three norms, [Theorem 1.9](../lecture_notes/01_matrix_algebras.md#r-1-9).

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
W = X @ absX.inverse()                                  # |X| is invertible almost surely, see below
assert close(W.adjoint() @ W, M3.identity())            # W is unitary
assert abs((X @ W.adjoint()).trace.item() - s.sum().item()) < 1e-10    # Tr(X W*) = ||X||_1
assert abs(X.trace.item()) <= X.trace_norm().item()     # |Tr X| <= ||X||_1
```

`torch.randn` with a complex dtype draws independent entries with independent real and imaginary
parts $`\sim N(0,\tfrac12)`$, so `G` is a Ginibre matrix, which has a density (Proposition A.16(b)).
$`|X|`$ is invertible almost surely: $`(\det|X|)^2=\det(|X|^2)=\det(X^*X)=|\det X|^2`$, and $`|\det X|^2`$ is a polynomial
in the real and imaginary parts of the entries that equals $`1`$ at $`X=1`$, so it vanishes with
probability zero (Corollary A.7(b); this is also the first step of Theorem 5.3).

<a id="s-1-7"></a>

## 1.7 Batches: many operators at once

*Lecture notes:* [Section 1.7](../lecture_notes/01_matrix_algebras.md#s-1-7), [Proposition 1.10](../lecture_notes/01_matrix_algebras.md#r-1-10).

The batch axis is the reason the library exists: random matrix theory, Monte Carlo
and multistart optimisation need thousands of operators, and one batched call replaces a
Python loop over them. As an example we check an exact moment of the **Gaussian
unitary ensemble** (GUE).

By [lecture notes, Proposition 1.10](../lecture_notes/01_matrix_algebras.md#r-1-10), if $`G`$ has independent entries with independent real and imaginary
parts $`\sim N(0,\tfrac12)`$ and $`X=(G+G^*)/\sqrt{2n}`$, then $`\mathbb E\,\tfrac1n\mathrm{Tr}\,X^4=2+\tfrac1{n^2}`$.

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
value to within a few standard errors. The variable $`\tfrac1n\mathrm{Tr}\,X^4`$ is a polynomial in
Gaussian variables, so it has finite variance $`\sigma^2`$ (Lemma A.23: all Gaussian moments are finite;
products are handled by $`|y_1\cdots y_r|\le\sum_i|y_i|^r`$),
and by Chebyshev's inequality the mean of $`B`$ independent samples deviates from the exact value
by more than $`5\sigma/\sqrt B`$ with probability at most $`1/25`$. The code uses the sample standard
deviation in place of $`\sigma`$, and the fixed seed makes the check reproducible. This is the
pattern for every check against proved results in this tutorial: an exact statement, a numerical
estimate, and an explicit tolerance.

<a id="s-1-8"></a>

## 1.8 What to compare with

On single matrices, `numpy.linalg` and `scipy.linalg` (`sqrtm`, `expm`, `logm`,
`polar`) are the reference implementations, and the checks above can be repeated
with them. The library adds three things that matter later: batching (Section 1.7),
operators that carry a block structure (Chapter 2), and automatic differentiation
(PyTorch's autograd) through its operations.

<a id="x-1"></a>

## Exercises

[Exercises 1 and 3](../lecture_notes/01_matrix_algebras.md#x-1) of the lecture notes, numerically:
- Check $`(fg)(A)=f(A)g(A)`$ for $`f=\sin`$, $`g=\exp`$ on a random self-adjoint $`A`$ (Exercise 1).
- The **geometric mean** of $`A,B>0`$ is $`A\#B=A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}`$, the unique positive
  solution $`X`$ of $`XA^{-1}X=B`$ (Exercise 3). Check this equation and $`A\#B=B\#A`$ numerically.

<details>
<summary>Solutions</summary>

```python
Hs = M3.operator_from_eigenvalues(lambda d: torch.randn(1, d), force_self_adjoint=True)
assert close(Hs.apply_function(lambda w: torch.sin(w) * torch.exp(w)),
             Hs.apply_function(torch.sin) @ Hs.apply_function(torch.exp))
```

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

**Next:** [Tutorial 2](02_von_neumann_algebras.md): commutants, the centre and the algebra $`\bigoplus_cM_{k_c}(\mathbb C)`$ in the library.
