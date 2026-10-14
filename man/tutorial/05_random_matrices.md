**English** · [Español](es/05_random_matrices.md) · [Français](fr/05_random_matrices.md) · [Deutsch](de/05_random_matrices.md) · [中文](zh/05_random_matrices.md) · [日本語](ja/05_random_matrices.md) · [Русский](ru/05_random_matrices.md)

# Tutorial 5. Random matrices in an algebra

This tutorial accompanies [Chapter 5](../lecture_notes/05_random_matrices.md) of the lecture notes, which constructs the Haar measure
on the unitary group of $`\mathcal M=\bigoplus_cM_{k_c}`$ ([Corollary 5.4](../lecture_notes/05_random_matrices.md#r-5-4)), computes its first and second
moments ([Proposition 5.7](../lecture_notes/05_random_matrices.md#r-5-7), [Theorem 5.10](../lecture_notes/05_random_matrices.md#r-5-10)) and shows why level statistics must be
computed sector by sector ([Theorem 5.15](../lecture_notes/05_random_matrices.md#r-5-15)). We sample Haar unitaries and Gaussian ensembles with
the library and check these results by Monte Carlo: a mean must agree with the exact value within $`5`$
standard errors. Results without a theorem behind them are marked as numerical observations.

Throughout, $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$ with $`D=\sum_ck_c`$, $`P_c`$ are the
minimal central projections, and $`E`$ is the conditional expectation onto the centre of
Proposition 2.8. $`U(k)`$ is the group of unitary $`k\times k`$ matrices and
$`U(\mathcal M)=\{U\in\mathcal M:U^*U=1\}`$ the unitary group of $`\mathcal M`$. Random variables
live on a probability space with expectation $`\mathbb E`$; the expectation of a random matrix is
taken entrywise. A *standard complex Gaussian* variable has independent real and imaginary parts
$`\sim N(0,\tfrac12)`$, so that $`\mathbb E|z|^2=1`$ (this is what `torch.randn` draws for a complex
dtype, according to its documentation; the first code block checks the second moments); a *Ginibre matrix* has independent standard complex Gaussian entries.

<a id="s-5-1"></a>

## 5.1 Haar measure

*Lecture notes:* [Section 5.1](../lecture_notes/05_random_matrices.md#s-5-1); [Theorem 5.3](../lecture_notes/05_random_matrices.md#r-5-3) (the unitary factor of the QR decomposition of a Ginibre matrix, with the phases fixed by $`\Lambda=\mathrm{diag}(R'_{ii}/|R'_{ii}|)`$, is Haar distributed), [Corollary 5.4](../lecture_notes/05_random_matrices.md#r-5-4) (the blocks of a Haar unitary of $`U(\mathcal M)`$ are independent and Haar distributed) and [Corollary 5.5](../lecture_notes/05_random_matrices.md#r-5-5) ($`\mathbb EU_{ij}=0`$, $`\mathbb E|U_{ij}|^2=1/k`$). The Haar measure of $`SU(k)`$ is the law of $`e^{-i\theta/k}U`$, $`\theta=\arg\det U`$.

Without the correction $`\Lambda`$ of [Theorem 5.3](../lecture_notes/05_random_matrices.md#r-5-3) the law of $`Q'`$ depends on the sign convention
of the QR routine. For the routine called by `torch.linalg.qr` (LAPACK) the code below finds
$`\mathbb E\mathrm{Re}\,\mathrm{Tr}\,Q'`$ many standard errors away from $`0`$, whereas a Haar unitary has
$`\mathbb EU=0`$ (Corollary 5.5): so this $`Q'`$ is not Haar distributed. This is a statement about that
routine, checked by the code, not a theorem.

The library draws $`U(\mathcal M)`$-Haar unitaries with `random_unitary_operator`, one batched QR
decomposition per block size for all sectors and samples at once (on a GPU, Gram–Schmidt, which
gives the same $`Q`$ by the uniqueness in Lemma 5.2). The measure `'haar'` is Theorem 5.3 in every
block; for a complex algebra `'haar_su'` multiplies by $`e^{-i\theta/k}`$ as above; `'coe'` is the
ensemble $`S=W^TW`$ of Exercise 3 (the library's other measures are not used in this course). We check
Corollaries 5.4 and 5.5, and that the library's sampler is exactly the corrected QR of Theorem 5.3:

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, SpinChain, channels

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def mc(samples, exact, name):
    """Monte Carlo mean of real samples with its standard error; must agree within 5 standard errors."""
    mean, err = samples.mean().item(), samples.std().item() / math.sqrt(samples.numel())
    print(f"{name}: {mean:.4f} +- {err:.4f}   (theory {exact:.4f})")
    assert abs(mean - exact) <= 5 * err + 1e-12
    return mean

alg = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')        # M_1 (+) M_2 (+) M_3, D = 6
B = 100_000
U = alg.random_unitary_operator(B)                                   # Haar on U(1) x U(2) x U(3)
assert close(U.adjoint() @ U, alg.identity())                        # unitary elements of the algebra
u = [U.matrix[:, c, :k, :k] for c, k in enumerate(alg.k_factors)]   # the blocks U_c
for c, k in enumerate(alg.k_factors):                                # Corollary 5.5
    mc(u[c][:, 0, 0].real, 0.0, f"E Re (U_{c})_11")
    mc(u[c][:, 0, 0].abs() ** 2, 1 / k, f"E |(U_{c})_11|^2")
mc(u[1][:, 0, 0].abs() ** 2 * u[2][:, 0, 0].abs() ** 2, 1 / 6, "E |(U_1)_11|^2 |(U_2)_11|^2")   # independent blocks

z = torch.randn(B, dtype=dt)                                         # torch.randn draws standard complex Gaussians:
mc(z.real ** 2, 0.5, "E (Re z)^2")                                   # Re z, Im z ~ N(0, 1/2), uncorrelated
mc(z.imag ** 2, 0.5, "E (Im z)^2")
mc(z.real * z.imag, 0.0, "E Re z Im z")

SU = alg.random_unitary_operator(4, measure='haar_su').matrix
assert torch.allclose(torch.linalg.det(SU[:, 2, :3, :3]), torch.ones(4, dtype=dt))

torch.manual_seed(1)
Z = torch.randn(B, 2, 2, dtype=dt)                                   # Ginibre matrices
torch.manual_seed(1)
V = alg.random_unitary(2, batch_size=B)                              # the library, same random numbers
Q, R = torch.linalg.qr(Z)                                            # LAPACK: diag(R) real, either sign
d = R.diagonal(dim1=-2, dim2=-1)
assert close(V, Q * (d / d.abs()).unsqueeze(-2))                     # Theorem 5.3: V = Q Lambda
raw = Q.diagonal(dim1=-2, dim2=-1).sum(-1).real                      # Re Tr Q without the correction
print(f"E Re Tr Q without the phase correction: {raw.mean():.4f} +- {raw.std() / math.sqrt(B):.4f}")
assert abs(raw.mean()) > 20 * raw.std() / math.sqrt(B)               # far from 0: not Haar
```

*Numerical observation.* Without the correction the code finds
$`\mathbb E\mathrm{Re}\,\mathrm{Tr}\,Q'\approx-0.84`$ for $`k=2`$, whereas a Haar unitary has $`\mathbb EU=0`$.

<a id="s-5-2"></a>

## 5.2 Moments and the Weingarten calculus

*Lecture notes:* [Section 5.2](../lecture_notes/05_random_matrices.md#s-5-2); [Proposition 5.7](../lecture_notes/05_random_matrices.md#r-5-7) ($`\mathbb E[UXU^*]=E(X)`$), [Lemma 5.8](../lecture_notes/05_random_matrices.md#r-5-8) (the commutant of all $`V\otimes V`$ is $`\mathrm{span}\{1,F\}`$, $`F`$ the swap), [Theorem 5.10](../lecture_notes/05_random_matrices.md#r-5-10) (the Weingarten formula of order two) and [Corollary 5.11](../lecture_notes/05_random_matrices.md#r-5-11).

We check Proposition 5.7 against `channels.center_expectation`, Lemma 5.8 by linear algebra (the
commutant of three sampled $`V\otimes V`$, computed as a null space, contains $`F`$ and is $`2`$-dimensional, so by
Lemma 5.8 it is already $`\mathrm{span}\{1,F\}`$), and Theorem 5.10 and Corollary 5.11 by Monte Carlo with the
$`U(3)`$ block of the samples above:

```python
X = alg.from_blocks([torch.randn(k, k, dtype=dt) for k in alg.k_factors])      # a fixed element
T = (U @ X @ U.adjoint()).matrix                                    # B samples of U X U*
mean, err = T.mean(0), T.std(0) / math.sqrt(B)
EX = channels.center_expectation(alg)(X).matrix[0]
random = err > 1e-12                                                 # entries that are not constant
z = ((mean - EX).abs()[random] / err[random]).max().item()
print(f"twirl over U(M) vs E(X): largest deviation {z:.2f} standard errors")
assert ((mean - EX).abs() <= 5 * err + 1e-12).all()                 # Proposition 5.7

def commutant_dim(ops, tol=1e-9):
    """Dimension of {X : XA = AX for all A in ops}, as the null space of X -> (AX - XA)_A."""
    d = ops[0].shape[0]
    eye = torch.eye(d, dtype=dt)
    ops = [A.contiguous() for A in ops]                              # torch.kron needs contiguous input
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    return int((torch.linalg.svdvals(L) < tol).sum())

k = 3
F = torch.zeros(k * k, k * k, dtype=dt)
for i in range(k):
    for j in range(k):
        F[j * k + i, i * k + j] = 1                                  # F (e_i (x) e_j) = e_j (x) e_i
VV = [torch.kron(W, W) for W in alg.random_unitary(k, batch_size=3).contiguous()]
assert all(close(F @ W, W @ F) for W in VV)
assert commutant_dim(VV) == 2                                        # Lemma 5.8: span{1, F}

u3 = u[2]                                                            # Haar on U(3)
mc(u3[:, 0, 0].abs() ** 4, 2 / (k * (k + 1)), "E|U_11|^4")
mc(u3[:, 0, 0].abs() ** 2 * u3[:, 0, 1].abs() ** 2, 1 / (k * (k + 1)), "E|U_11|^2 |U_12|^2")
mc((u3[:, 0, 0] * u3[:, 1, 1] * (u3[:, 0, 1] * u3[:, 1, 0]).conj()).real,
   -1 / ((k - 1) * k * (k + 1)), "E U_11 U_22 conj(U_12 U_21)")

n_s = 20_000                                                         # Theorem 5.10 for a random X
X2 = torch.randn(k * k, k * k, dtype=dt)
trX, trFX = X2.diagonal().sum(), (F @ X2).diagonal().sum()
alpha, beta = (trX - trFX / k) / (k ** 2 - 1), (trFX - trX / k) / (k ** 2 - 1)
UU = torch.einsum('bij,bkl->bikjl', u3[:n_s], u3[:n_s]).reshape(n_s, k * k, k * k)   # U (x) U
T2 = UU @ X2 @ UU.mH
mean, err = T2.mean(0), T2.std(0) / math.sqrt(n_s)
assert ((mean - alpha * torch.eye(k * k, dtype=dt) - beta * F).abs() <= 5 * err).all()
```

<a id="s-5-3"></a>

## 5.3 Level statistics

*Lecture notes:* [Section 5.3](../lecture_notes/05_random_matrices.md#s-5-3): the GOE $`H=(G+G^T)/2`$, the GUE $`H=(G+G^*)/\sqrt2`$, and the spacing ratio $`r_i=\min(s_i,s_{i+1})/\max(s_i,s_{i+1})`$ of the spacings $`s_i=\lambda_{i+1}-\lambda_i`$; [Proposition 5.13](../lecture_notes/05_random_matrices.md#r-5-13) (the Wigner surmise for $`2\times2`$ matrices) and [Proposition 5.14](../lecture_notes/05_random_matrices.md#r-5-14) (the mean ratio is $`2\ln2-1`$ for independent levels, $`4-2\sqrt3`$ and $`\frac{2\sqrt3}\pi-\frac12`$ for the $`3\times3`$ GOE and GUE, where $`\mathbb Er=I_1/I_0`$ with $`I_j=\int_0^1x^jf_\beta(x)\,dx`$).

The function `SpinChain.level_spacing_ratio` (in `chains`) returns the mean ratio of each row of
levels. We check Propositions 5.13 and 5.14 and record the GOE value for $`n=64`$:

```python
ratio = SpinChain.level_spacing_ratio
n2 = 200_000
alg_r = TypeIAlgebra([2], [2], complex_valued=False, precision='double')     # real symmetric 2 x 2
alg_c = TypeIAlgebra([2], [2], precision='double')
G = torch.randn(n2, 1, 2, 2, dtype=torch.float64)
Gc = torch.randn(n2, 1, 2, 2, dtype=dt)
goe2 = alg_r.operator((G + G.mT) / 2, is_self_adjoint=True)
gue2 = alg_c.operator((Gc + Gc.mH) / math.sqrt(2), is_self_adjoint=True)
cdf = {1: lambda t: 1 - math.exp(-math.pi * t ** 2 / 4),
       2: lambda t: math.erf(2 * t / math.sqrt(math.pi)) - 4 * t / math.pi * math.exp(-4 * t ** 2 / math.pi)}
for beta, H, mean_s in ((1, goe2, math.sqrt(math.pi)), (2, gue2, 4 / math.sqrt(math.pi))):
    lam = H.eigenvalues()[0]                                         # (n2, 2), ascending
    S = (lam[:, 1] - lam[:, 0]) / mean_s
    mc(S, 1.0, f"beta={beta}: E S")
    for t in (0.25, 1.0, 2.0):
        mc((S <= t).double(), cdf[beta](t), f"beta={beta}: P(S <= {t})")

levels = torch.rand(4000, 50, dtype=torch.float64)                   # independent levels
mc(ratio(levels), 2 * math.log(2) - 1, "<r>, independent uniform levels")      # Proposition 5.14(a)

n3 = 100_000                                                         # Proposition 5.14(b): 3 x 3
G = torch.randn(n3, 1, 3, 3, dtype=torch.float64)
Gc = torch.randn(n3, 1, 3, 3, dtype=dt)
goe3 = TypeIAlgebra([3], [3], complex_valued=False, precision='double').operator((G + G.mT) / 2, is_self_adjoint=True)
gue3 = TypeIAlgebra([3], [3], precision='double').operator((Gc + Gc.mH) / math.sqrt(2), is_self_adjoint=True)
h = 1e-4
xs = (torch.arange(int(1 / h), dtype=torch.float64) + 0.5) * h      # midpoint rule on [0, 1]
for beta, H, exact in ((1, goe3, 4 - 2 * math.sqrt(3)), (2, gue3, 2 * math.sqrt(3) / math.pi - 0.5)):
    f = (xs + xs ** 2) ** beta / (1 + xs + xs ** 2) ** (1 + 1.5 * beta)
    assert abs((xs * f).sum() / f.sum() - exact) < 1e-7              # E r = I_1 / I_0
    r3 = ratio(H.eigenvalues()[0])                                   # one ratio per matrix
    mc(r3, exact, f"beta={beta}: <r>, 3 x 3")
    mc((r3 <= 0.5).double(), (f * (xs <= 0.5)).sum().item() / f.sum().item(), f"beta={beta}: P(r <= 1/2), 3 x 3")

n = 64                                                               # numerical observation, no theorem
alg64 = TypeIAlgebra([n], [n], complex_valued=False, precision='double')
G = torch.randn(800, 1, n, n, dtype=torch.float64)
lam = alg64.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()[0]
r_goe = ratio(lam[:, n // 4: 3 * n // 4])                            # 800 samples of the mean ratio
print(f"<r>, GOE(64), middle half of the spectrum: {r_goe.mean():.4f} +- {r_goe.std() / math.sqrt(800):.4f}")

def same_mean(a, b, name):
    """Two independent Monte Carlo samples of the same law: means agree within 5 combined standard errors."""
    ea, eb = a.std().item() / math.sqrt(a.numel()), b.std().item() / math.sqrt(b.numel())
    print(f"{name}: {a.mean():.4f} vs {b.mean():.4f} +- {math.hypot(ea, eb):.4f}")
    assert abs(a.mean().item() - b.mean().item()) <= 5 * math.hypot(ea, eb)
```

*Numerical observation.* The code finds $`\langle r\rangle\approx0.53`$ for the middle half of the spectrum
of the GOE with $`n=64`$, below the exact $`3\times3`$ value $`4-2\sqrt3\approx0.536`$ and far above the value
$`2\ln2-1\approx0.386`$ of independent levels. This number is an output of the computation, not a
theorem; below it serves as the reference for GOE blocks of the same size $`64`$, whose spacing ratios have
exactly the same law.

<a id="s-5-4"></a>

## 5.4 Why sectors must not be mixed

*Lecture notes:* [Section 5.4](../lecture_notes/05_random_matrices.md#s-5-4); by [Theorem 5.15](../lecture_notes/05_random_matrices.md#r-5-15) levels of different independent blocks do not repel: for two GOE blocks the smallest spacing $`\delta`$ of the merged spectrum has $`P(\delta<\varepsilon)`$ at least of order $`\varepsilon`$, whereas inside one GOE block it is at most of order $`\varepsilon^2`$.

With more sectors the effect is stronger. We do not prove a limit theorem for many sectors; the
following experiment is a numerical observation. In the library the blocks are separate from the start:
`eigenvalues()` returns one spectrum per sector, and mixing them is an explicit `torch.cat`. We put $`C`$
independent GOE blocks of size $`64`$ in one algebra and compare the mean ratio per sector with that of the
merged spectrum. Each sector is a GOE matrix of size $`64`$, so its ratios have exactly the law of `r_goe`
above, and the code tests this; the behaviour of the merged spectrum is only observed:

```python
k, merged_r = 64, {}
for C in (1, 2, 4, 8, 16):
    alg_C = TypeIAlgebra([k] * C, [k] * C, complex_valued=False, precision='double')
    G = torch.randn(800 // C, C, k, k, dtype=torch.float64)
    lam = alg_C.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()    # one spectrum per sector
    per_sector = torch.stack([ratio(l[:, k // 4: 3 * k // 4]) for l in lam])
    merged = torch.cat(lam, dim=-1).sort(dim=-1).values                       # the spectrum on C^D
    r_m = ratio(merged[:, k * C // 4: 3 * k * C // 4])
    merged_r[C] = r_m.mean().item()
    print(f"C = {C:2d}: <r> per sector {per_sector.mean():.4f}, merged {merged_r[C]:.4f}"
          f" +- {r_m.std() / math.sqrt(len(r_m)):.4f}")
    same_mean(per_sector.flatten(), r_goe, f"C = {C:2d}: sectors vs GOE(64)")   # every sector is GOE(64)
assert merged_r[1] > merged_r[2] > merged_r[4] > merged_r[16]       # numerical observation
assert abs(merged_r[16] - (2 * math.log(2) - 1)) < 0.01              # numerical observation
```

*Numerical observation.* Every sector has $`\langle r\rangle\approx0.53`$, the GOE(64) value of Section 5.3
(this agreement is a test of a fact: the sectors are GOE(64) matrices). The merged spectrum drops to
$`0.42`$ for two sectors and to about $`0.39`$, within $`0.01`$ of the value $`2\ln2-1`$ of independent levels
(Proposition 5.14(a)), for $`8`$ and $`16`$ sectors. That the merged spectrum approaches independent levels as
$`C`$ grows is what the code observes; it is not proved in this course.

The same happens in a physical model. The Heisenberg chain with random fields (Chapter 9) conserves
the magnetisation, so its Hamiltonian lies in $`\bigoplus_{N=0}^LM_{\binom LN}`$ (Proposition 9.1). Its blocks are *not*
independent — all of them depend on the same random fields — so Theorem 5.15 does not apply as stated;
the experiment below shows the same effect, as a numerical observation.

```python
chain = SpinChain(10, complex_valued=False)                         # sectors N = 0, ..., 10
lam = chain.random_field_heisenberg(W=0.5, batch_size=20).eigenvalues()
k5 = chain.sector_dim(5)                                            # 252
r_sector = ratio(lam[5][:, k5 // 4: 3 * k5 // 4]).mean().item()
spectrum = torch.cat(lam, dim=-1).sort(dim=-1).values               # all 1024 levels mixed
r_all = ratio(spectrum[:, 256:768]).mean().item()
print(f"Heisenberg chain, L = 10: <r> in sector N = 5: {r_sector:.4f}, all sectors mixed: {r_all:.4f}")
assert r_sector > 0.49 and r_all < 0.42                             # numerical observation
```

*Numerical observation.* Inside the sector $`N=5`$ the code finds $`\langle r\rangle\approx0.53`$, close to the
GOE(64) value above, and the spectrum with all sectors mixed gives $`0.39`$, close to $`2\ln2-1`$.

<a id="s-5-5"></a>

## 5.5 What to compare with

Haar unitaries: `scipy.stats.unitary_group` and `ortho_group`, QuTiP's Haar-random unitaries, and any
`numpy`/`torch` QR followed by the phase correction.
Weingarten calculus: the symbolic packages IntU (Mathematica; Z. Puchała, J. A. Miszczak), RTNI
(Mathematica and Python; M. Fukuda, R. König, I. Nechita, J. Phys. A 52 (2019), 425303) and the
Python package Haarpy compute moments of every order symbolically; Theorem 5.10 is their order-two
case and Theorem 5.12 (proved in Appendix B) the general formula they implement. Level statistics are
usually computed with a few lines of `numpy`, and exact-diagonalisation codes such as QuSpin supply spectra sector by sector. What the library
adds is the unitary group of an algebra as one object: independent Haar blocks of all sizes drawn in
batched calls, the twirl identified with the conditional expectation
`channels.center_expectation`, and spectra returned per sector, so that mixing sectors is never
the default.

<a id="x-5"></a>

## Exercises

The numerical parts of the [exercises](../lecture_notes/05_random_matrices.md#x-5) of the lecture notes:
1. For $`U`$ Haar on $`U(\mathcal M)`$, $`\mathcal M=M_1\oplus M_2\oplus M_3`$, check $`\mathbb E|\mathrm{Tr}\,U|^2=C=3`$ and
   $`\mathbb E|\mathrm{Tr}\,U|^4=17`$.
2. Compute the Weingarten matrix of order $`3`$ for $`k=3`$ ([Theorem 5.12](../lecture_notes/05_random_matrices.md#r-5-12)) numerically and check
   $`\mathbb E|U_{11}|^6=\frac1{10}`$.
3. The **circular orthogonal ensemble** (`measure='coe'`) is $`S=W^TW`$ with $`W`$ Haar on $`U(k)`$. Check that the
   library's COE is $`W^TW`$, that $`S`$ is unitary and symmetric, that $`\mathbb E\mathrm{Tr}\,S^m=0`$ for $`m=1,2,3`$
   and that the mean density of the eigenphases is uniform. Compute $`\langle r\rangle`$ for the eigenphases of
   the COE with $`k=32`$ and compare it with the GOE value of Section 5.3 (a numerical observation).

<details>
<summary>Solutions</summary>

*Exercise 1.*

```python
tr = U.trace                                                        # Tr U = sum_c Tr U_c, shape (B,)
mc(tr.abs() ** 2, 3.0, "E|Tr U|^2")
mc(tr.abs() ** 4, 17.0, "E|Tr U|^4")
```

*Exercise 2.*

```python
import itertools

def n_cycles(p):                                                    # number of cycles of a permutation
    seen, count = set(), 0
    for i in range(len(p)):
        count += i not in seen
        while i not in seen:
            seen.add(i)
            i = p[i]
    return count

perms = list(itertools.permutations(range(3)))                      # S_3; s.index(j) = s^{-1}(j)
Gram = torch.tensor([[3.0 ** n_cycles([s.index(p[i]) for i in range(3)]) for p in perms] for s in perms],
                    dtype=torch.float64)                            # k^{#cycles(s^{-1} p)}, k = 3
Wg3 = torch.linalg.inv(Gram)                                        # the Weingarten matrix, n = k = 3
assert abs(Wg3.sum().item() - 6 / (3 * 4 * 5)) < 1e-12
mc(u3[:, 0, 0].abs() ** 6, 0.1, "E|U_11|^6")
```

*Exercise 3.*

   *Numerical observation.* The code finds $`\langle r\rangle\approx0.53`$ for the COE with $`k=32`$, within
   $`0.015`$ of the GOE(64) value of Section 5.3. This agreement is an output of the computation; no
   theorem of this course relates the spacing statistics of the COE to those of the GOE (for the
   circular ensembles see F. J. Dyson, J. Math. Phys. 3 (1962), 140–156).

```python
alg_s = TypeIAlgebra([32], [32], precision='double')
torch.manual_seed(2)
S = alg_s.random_unitary_operator(400, measure='coe').matrix[:, 0]
torch.manual_seed(2)
W = alg_s.random_unitary(32, batch_size=400)                         # same random numbers
assert close(S, W.mT @ W)                                            # the library's COE is W^T W
assert close(S, S.mT) and close(S @ S.mH, torch.eye(32, dtype=dt))
for m in (1, 2, 3):                                                  # E Tr S^m = 0
    trm = torch.linalg.matrix_power(S, m).diagonal(dim1=-2, dim2=-1).sum(-1)
    mc(trm.real, 0.0, f"E Re Tr S^{m}")
    mc(trm.imag, 0.0, f"E Im Tr S^{m}")
phases = torch.linalg.eigvals(S).angle().sort(dim=-1).values         # in (-pi, pi]
mc(((phases > 0) & (phases <= 1)).sum(-1).double(), 32 / (2 * math.pi), "E #eigenphases in (0, 1]")
s = torch.diff(torch.cat([phases, phases[:, :1] + 2 * math.pi], dim=-1), dim=-1)   # around the circle
r = torch.minimum(s, s.roll(-1, -1)) / torch.maximum(s, s.roll(-1, -1))
r_coe = r.mean(-1)
print(f"<r>, COE(32) eigenphases: {r_coe.mean():.4f}; GOE(64), Section 5.3: {r_goe.mean():.4f}")
assert abs(r_coe.mean() - r_goe.mean()) < 0.015                      # numerical observation
```

</details>

**Next:** [Tutorial 6](06_quantum_channels.md): quantum channels in the library.
