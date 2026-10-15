**English** · [Español](es/11_optimisation.md) · [Français](fr/11_optimisation.md) · [Deutsch](de/11_optimisation.md) · [中文](zh/11_optimisation.md) · [日本語](ja/11_optimisation.md) · [Русский](ru/11_optimisation.md)

# Tutorial 11. Optimisation and trace inequalities

This tutorial accompanies [Chapter 11](../lecture_notes/11_optimisation.md) of the lecture notes, which proves the extremes of three
trace functionals and explains why sampling unitaries fails to find them. Here the parametrisations of
constraint sets (`optimize.PositiveParam`, `UnitaryParam`, `SelfAdjointParam`) and the multistart optimiser
(`optimize.extremize`) of the library are checked against those theorems. With $`\Delta`$ the Michelson contrast
and $`\ell(\delta)=(1-\delta)/(1+\delta)`$, the three functionals are:

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

*Lecture notes:* [Section 11.1](../lecture_notes/11_optimisation.md#s-11-1); [Lemma 11.1](../lecture_notes/11_optimisation.md#r-11-1), [Proposition 11.2](../lecture_notes/11_optimisation.md#r-11-2).

`optimize.PositiveParam` with prescribed contrast $`\delta`$ implements the parametrisation of
[lecture notes, Section 11.1](../lecture_notes/11_optimisation.md#s-11-1): per start it draws two distinct slots $`s_\pm`$ (`spectrum.pin_max`,
`spectrum.pin_min`; every ordered pair has probability $`1/(D(D-1))`$) and a Haar-random unitary $`U_0`$
(`unitary.base`, drawn by `random_unitary_operator`, [Lemma 11.10](../lecture_notes/11_optimisation.md#r-11-10)), and maps the parameters $`\theta`$
(`spectrum.theta`) and $`G`$ (`unitary.G`) to
```math
X(\theta,G)=U\mathrm{diag}(\lambda)U^*,\qquad U=U_0\,e^{G-G^*},\qquad
\lambda_{s_+}=1,\quad\lambda_{s_-}=\ell(\delta),\quad\lambda_s=\ell+(1-\ell)\tfrac{1+\sin\theta_s}2\ \text{otherwise}.
```
`UnitaryParam` is the map $`G\mapsto U_0e^{G-G^*}`$ alone, and `SelfAdjointParam` adds signs,
$`X=U\mathrm{diag}(\sigma\lambda)U^*`$. By [Proposition 11.2](../lecture_notes/11_optimisation.md#r-11-2) the image is the set of operators of
the contrast shell $`\mathcal P_\delta`$ with $`1`$ and $`\ell(\delta)`$ in the drawn slots. A free contrast in $`[a,b]`$ is one more
parameter; the library's default is $`[a,b]=[0,0.999]`$, and the code below keeps it.

The code checks the constraint, and Proposition 11.2 constructively: the parameters of the proof reproduce a target.

```python
import math
import scipy.linalg
import torch
from torch_vn_algebra import TypeIAlgebra, optimize as opt

torch.manual_seed(0)
alg = TypeIAlgebra([3, 3], [3, 3], precision='double')          # M_3 (+) M_3, D = 6
ell = lambda delta: (1 - delta) / (1 + delta)                     # lambda_min for contrast delta
d = torch.tensor([0.0, 0.25, 0.5, 0.9], dtype=torch.float64)
x = opt.PositiveParam(alg, 4, delta=d).operator()                 # four starts, one contrast each
assert torch.allclose(x.michelson_contrast, d, atol=1e-12)        # the constraint holds exactly
assert torch.allclose(x.lambda_max, torch.ones(4, dtype=torch.float64), atol=1e-12)

P = opt.PositiveParam(alg, 1, delta=torch.tensor([0.5], dtype=torch.float64))
lo = ell(0.5)
mu = lo + (1 - lo) * torch.rand(1, 2, 3, dtype=torch.float64)    # target spectrum in [lo, 1] ...
(cp, ip), (cm, im) = P.spectrum.pin_max[0].tolist(), P.spectrum.pin_min[0].tolist()
mu[0, cp, ip], mu[0, cm, im] = 1.0, lo                            # ... with 1 and lo in the drawn slots
W = alg.random_unitary_operator(1).matrix                         # target eigenvectors
target = W @ torch.diag_embed(mu.to(W.dtype)) @ W.mH
with torch.no_grad():                                             # the parameters of the proof
    P.spectrum.theta.copy_(torch.asin((2 * (mu - lo) / (1 - lo) - 1).clamp(-1, 1)))
    V = P.unitary.base.mH @ W                                     # U0^* W = exp(S), G = S / 2
    S = torch.stack([torch.from_numpy(scipy.linalg.logm(V[0, c].numpy())) for c in range(2)])
    P.unitary.G.copy_(S[None] / 2)
err = (P.operator().matrix - target).abs().max().item()
print(f"slots of 1 and lo: {(cp, ip)}, {(cm, im)};  reconstruction error {err:.1e}")
assert err < 1e-10
```

<a id="s-11-2"></a>

## 11.2 Von Neumann's trace inequality

*Lecture notes:* [Section 11.2](../lecture_notes/11_optimisation.md#s-11-2); [Theorem 11.4](../lecture_notes/11_optimisation.md#r-11-4), [Corollary 11.5](../lecture_notes/11_optimisation.md#r-11-5), [Proposition 11.6](../lecture_notes/11_optimisation.md#r-11-6), [Corollary 11.7](../lecture_notes/11_optimisation.md#r-11-7).

By [lecture notes, Theorem 11.4](../lecture_notes/11_optimisation.md#r-11-4) (von Neumann's trace inequality),
$`|\mathrm{Tr}(AB)|\le\sum_is_i(A)s_i(B)`$, and by [Corollary 11.5](../lecture_notes/11_optimisation.md#r-11-5),
$`\sum_ip^\downarrow_iq^\uparrow_i\le\mathrm{Tr}(PQ)\le\sum_ip^\downarrow_iq^\downarrow_i`$ for self-adjoint $`P,Q`$ with
eigenvalues $`p,q`$, and $`\mathrm{Tr}(PQ)\le\|PQ\|_1\le\sum_ip^\downarrow_iq^\downarrow_i`$ for $`P,Q\ge0`$. For Exp. 2 with
fixed spectra, [Proposition 11.6](../lecture_notes/11_optimisation.md#r-11-6) gives
$`\sup z=-\inf z=\sum_c\sum_ix^\downarrow_{c,i}(y^\downarrow_{c,i}-y^\uparrow_{c,i})`$, attained for
$`X_0=\mathrm{diag}(x^\downarrow)`$, $`V=J`$ the reversal permutation and $`Y_0=\mathrm{diag}(y^\uparrow)`$ (resp.
$`\mathrm{diag}(y^\downarrow)`$), and [Corollary 11.7](../lecture_notes/11_optimisation.md#r-11-7) gives the envelope
$`\sup z=-\inf z=\sum_c\lfloor k_c/2\rfloor(1-\ell_X)(1-\ell_Y)`$.

Here $`1-\ell_X=2\delta/(1+\delta)`$ for $`\Delta(|X|)=\delta`$; for odd $`k_c`$ the middle eigenvalue is irrelevant. The code
checks Theorem 11.4, Corollary 11.5 and the extremal configurations of Proposition 11.6 for $`k=3`$, $`C=2`$.

```python
M5 = TypeIAlgebra([5], [5], precision='double')
A = M5.operator(torch.randn(2000, 1, 5, 5, dtype=torch.complex128))         # Ginibre, no structure
B = M5.operator(torch.randn(2000, 1, 5, 5, dtype=torch.complex128))
sA, sB = torch.linalg.svdvals(A.matrix[:, 0]), torch.linalg.svdvals(B.matrix[:, 0])   # decreasing
assert ((A @ B).trace.abs() <= (sA * sB).sum(-1) + 1e-12).all()             # Theorem 11.4
Ph, Qh = (M5.operator_from_eigenvalues(lambda k: torch.randn(2000, k), batch_size=2000, force_self_adjoint=True)
          for _ in range(2))
p, q = Ph.eigenvalues()[0], Qh.eigenvalues()[0]                             # increasing order
t = (Ph @ Qh).trace.real
assert (t <= (p * q).sum(-1) + 1e-12).all() and (t >= (p * q.flip(-1)).sum(-1) - 1e-12).all()   # Cor. 11.5(a)
Pp, Qp = Ph @ Ph, Qh @ Qh                                                   # positive, spectra p^2, q^2
nuc = (Pp @ Qp).trace_norm()
assert ((Pp @ Qp).trace.real <= nuc + 1e-12).all()                          # Cor. 11.5(b)
assert (nuc <= ((p ** 2).sort(-1)[0] * (q ** 2).sort(-1)[0]).sum(-1) + 1e-10).all()

alg3 = TypeIAlgebra([3, 3], [3, 3], precision='double')
x_dn = torch.tensor([1.0, 0.6, ell(0.5)], dtype=torch.float64)              # Delta = 0.5, middle value 0.6
y_dn = torch.tensor([1.0, 0.7, ell(0.8)], dtype=torch.float64)              # Delta = 0.8, middle value 0.7
J = torch.eye(3, dtype=torch.float64).flip(0)                               # reversal permutation
blocks = lambda m: alg3.from_blocks([m, m])                                 # the same block in both sectors

def z_exp2(X0, V, Y0):
    return ((X0 @ V @ Y0).trace_norm() - (X0 @ Y0).trace.real).item()

X0 = blocks(torch.diag(x_dn))
z_max = z_exp2(X0, blocks(J), blocks(torch.diag(y_dn.flip(0))))
z_min = z_exp2(X0, blocks(J), blocks(torch.diag(y_dn)))
exact = 2 * (x_dn * (y_dn - y_dn.flip(0))).sum().item()                     # Proposition 11.6
envelope = 2 * (1 - ell(0.5)) * (1 - ell(0.8))                              # Corollary 11.7, floor(3/2) = 1
print(f"Exp. 2, k = 3, C = 2: z_max = {z_max:.6f}, z_min = {z_min:.6f}, formula {envelope:.6f}")
assert abs(z_max - exact) < 1e-12 and abs(z_min + exact) < 1e-12 and abs(exact - envelope) < 1e-12
```

<a id="s-11-3"></a>

## 11.3 A trace inequality with an equality condition

*Lecture notes:* [Section 11.3](../lecture_notes/11_optimisation.md#s-11-3); [Lemma 11.8](../lecture_notes/11_optimisation.md#r-11-8), [Proposition 11.9](../lecture_notes/11_optimisation.md#r-11-9).

By [lecture notes, Proposition 11.9](../lecture_notes/11_optimisation.md#r-11-9), $`z=\mathrm{Tr}(Y|X|Y)-\mathrm{Tr}|YXY|\ge0`$ in Exp. 3,
with equality iff $`YX_+Y^2X_-Y=0`$, in particular if $`XY=YX`$. For $`X=\mathrm{diag}(1,-1,0)`$ and $`Y`$ with rows
$`(1,0,a)`$, $`(0,1,0)`$, $`(a,0,1)`$, $`-1<a<1`$, equality holds although $`XY\ne YX`$ for $`a\ne0`$. The code checks both
(with $`a=0.5`$); `X.abs()` is $`|X|`$ and `trace_norm()` is $`\|\cdot\|_1`$.

```python
def z_exp3(X, Y):
    return (Y @ X.abs() @ Y).trace.real - (Y @ X @ Y).trace_norm()

Xs = alg3.operator_from_eigenvalues(lambda k: torch.randn(500, k), batch_size=500, force_self_adjoint=True)
Yp = alg3.operator_from_eigenvalues(lambda k: torch.rand(500, k), batch_size=500, force_positive=True)
assert (z_exp3(Xs, Yp) >= -1e-12).all()                                     # Proposition 11.9

M3 = TypeIAlgebra([3], [3], precision='double')
Xe = M3.operator(torch.diag(torch.tensor([1.0, -1.0, 0.0]))[None, None])
Ye = M3.operator(torch.tensor([[1.0, 0, 0.5], [0, 1, 0], [0.5, 0, 1]])[None, None])
print(f"equality without commutation: z = {z_exp3(Xe, Ye).item():.1e}, "
      f"||XY - YX|| = {(Xe @ Ye - Ye @ Xe).operator_norm().item():.2f}")
assert abs(z_exp3(Xe, Ye).item()) < 1e-12 and (Xe @ Ye - Ye @ Xe).operator_norm().item() > 0.1
```

<a id="s-11-4"></a>

## 11.4 Sampling versus optimisation

*Lecture notes:* [Section 11.4](../lecture_notes/11_optimisation.md#s-11-4); [Lemma 11.10](../lecture_notes/11_optimisation.md#r-11-10), [Proposition 11.11](../lecture_notes/11_optimisation.md#r-11-11), [Corollary 11.12](../lecture_notes/11_optimisation.md#r-11-12), [Lemma 11.13](../lecture_notes/11_optimisation.md#r-11-13).

How large is $`|\mathrm{Tr}(AU)|`$ typically, compared with its maximum $`\|A\|_1`$? The library samples $`U`$ as
the factor $`Q`$ of the QR decomposition $`Z=QR`$, $`R_{ii}>0`$, of a **Ginibre matrix** $`Z`$ (independent entries with
density $`\pi^{-1}e^{-|z|^2}`$; Mezzadri's recipe), independently in every block.

By [lecture notes, Proposition 11.11](../lecture_notes/11_optimisation.md#r-11-11), $`\mathbb E|\mathrm{Tr}(AU)|^2=\|A\|_2^2/k`$ for Haar
$`U\in U(k)`$, and by [Corollary 11.12](../lecture_notes/11_optimisation.md#r-11-12), for $`X,Y\ge0`$ with spectra in $`[a,b]`$, $`a>0`$, the root mean square of
$`|\mathrm{Tr}(XUY)|`$ is at most $`\min(1/\sqrt k,b^2/(a^2k))`$ times its maximum over $`U`$, and
$`\mathbb P\big(|\mathrm{Tr}(XUY)|\ge\mathrm{Tr}(XY)\big)\le b^4/(a^4k^2)`$.

So $`N`$ Haar samples find a violation of $`|\mathrm{Tr}(XUY)|\le\mathrm{Tr}(XY)`$ with probability at most
$`Nb^4/(a^4k^2)`$ (a union of $`N`$ events), although one exists for every non-commuting pair (Proposition 3.9). The bound
is informative only for $`k>\sqrt N\,b^2/a^2`$ (for $`b/a=6`$ and $`N=2000`$: $`k>1609`$); numerically, samples miss the
violation already for much smaller $`k`$ (below). The code checks Proposition 11.11 by Monte Carlo (the sample mean
lies within five standard errors of the exact value, for this seeded run) and prints, for spectra uniform in
$`[0.2,1.2]`$ ($`b/a=6`$), the ratio of Corollary 11.12, $`\sup z`$ and the best $`z`$ of $`2000`$ samples, both relative to
$`\mathrm{Tr}(XY)`$.

```python
def random_pair(alg_k):
    pos = lambda k: 0.2 + torch.rand(1, k)                                  # spectrum in [a, b] = [0.2, 1.2]
    return tuple(alg_k.operator_from_eigenvalues(pos, force_positive=True) for _ in range(2))

M6 = TypeIAlgebra([6], [6], precision='double')
Xf, Yf = random_pair(M6)
N = 20_000
t2 = (Yf @ Xf @ M6.random_unitary_operator(N)).trace.abs() ** 2
m2, err, exact2 = t2.mean().item(), t2.std().item() / math.sqrt(N), (Yf @ Xf).frobenius_norm().item() ** 2 / 6
print(f"E|Tr(AU)|^2 = {m2:.4f} +- {err:.4f}   (exact ||A||_2^2 / k = {exact2:.4f})")
assert abs(m2 - exact2) < 5 * err                                           # Proposition 11.11

print(" k   rms/max   sup z/TrXY   best of 2000 samples")
for k in (2, 4, 8, 16, 32):
    Mk = TypeIAlgebra([k], [k], precision='double')
    Xk, Yk = random_pair(Mk)
    trXY, nuc = (Xk @ Yk).trace.real.item(), (Yk @ Xk).trace_norm().item()
    ratio = (Yk @ Xk).frobenius_norm().item() / (math.sqrt(k) * nuc)
    best = (Xk @ Mk.random_unitary_operator(2000) @ Yk).trace.abs().max().item()
    print(f"{k:2d}   {ratio:.4f}    {(nuc - trXY) / trXY:.2e}     {(best - trXY) / trXY:+.3f}")
    assert ratio <= min(1 / math.sqrt(k), 36 / k) + 1e-12                   # Corollary 11.12
    assert nuc > trXY and best < nuc and (k < 4 or best < trXY)            # from k = 4 on: no violation found
```

*Numerical observation* (the printed table; asserted only for this seeded run): the ratio stays below the proved bound
$`\min(1/\sqrt k,b^2/(a^2k))`$ and roughly halves whenever $`k`$ doubles, i.e. it falls like $`1/k`$, and from $`k=4`$ on the
best of $`2000`$ samples does not even reach $`z=0`$. A gradient search, in contrast, comes close to both extremes in the
run below (again a numerical observation). For the pair `Xf, Yf` both extremes over $`U`$ are known:

By [lecture notes, Lemma 11.13](../lecture_notes/11_optimisation.md#r-11-13), $`\min_U|\mathrm{Tr}(AU)|=0`$ for every $`A\in M_k(\mathbb C)`$, $`k\ge2`$.

So $`\inf_Uz=-\mathrm{Tr}(XY)`$ and $`\sup_Uz=\|YX\|_1-\mathrm{Tr}(XY)`$ (Proposition 3.9). One call of `extremize`
approximates both, with $`12`$ starts maximising $`z`$ and $`12`$ maximising $`-z`$ (to relative accuracy $`10^{-3}`$ in this
run, asserted below). By Proposition 11.11 the root mean square of
$`|\mathrm{Tr}(XUY)|`$ is $`\|YX\|_2/\sqrt k`$; with the bounds in the proof of Corollary 11.12 (spectra in $`[a,b]`$),
$`\mathbb E|\mathrm{Tr}(XUY)|^2\le b^4`$ and $`\mathrm{Tr}(XY)\ge ka^2`$, so Markov's inequality gives
$`\mathbb P\big(z\ge-(1-\varepsilon)\mathrm{Tr}(XY)\big)\le b^4/(\varepsilon^2a^4k^2)`$ for $`\varepsilon>0`$: for large $`k`$,
$`z`$ is close to its infimum with high probability.
*Numerical observation* (asserted for this seeded run): sampling comes within $`1\,\%`$ of the infimum, but does not
even reach $`z=0`$, while the supremum is positive:

```python
trXY = (Xf @ Yf).trace.real
sup_exact, inf_exact = ((Yf @ Xf).trace_norm() - trXY).item(), -trXY.item()
U6 = opt.UnitaryParam(M6, 24)                                               # 24 starts
pm = torch.tensor([1.0, -1.0], dtype=torch.float64).repeat_interleave(12)   # 12 maximise z, 12 maximise -z
best = opt.extremize(lambda: pm * ((Xf @ U6.operator() @ Yf).trace.abs() - trXY), [U6],
                     maximize=True, rounds=6, steps=30).reshape(2, 12).max(1)[0]
best_up, best_down = best[0].item(), -best[1].item()
z_mc = (Xf @ M6.random_unitary_operator(N) @ Yf).trace.abs() - trXY
print(f"sup z: exact {sup_exact:.6e}, ascent {best_up:.6e}, best sample {z_mc.max().item():+.4f}")
print(f"inf z: exact {inf_exact:.6f}, descent {best_down:.6f}, worst sample {z_mc.min().item():+.6f}")
assert best_up <= sup_exact + 1e-12 and best_down >= inf_exact - 1e-12    # attained values: one-sided
assert sup_exact - best_up < 1e-3 * sup_exact and best_down - inf_exact < 1e-3 * abs(inf_exact)
assert z_mc.min().item() - inf_exact < 1e-2 * abs(inf_exact)               # sampling finds the typical side ...
assert z_mc.max().item() < 0                                               # ... but misses the violation
```

<a id="s-11-5"></a>

## 11.5 Multistart search for envelopes

*Lecture notes:* [Section 11.5](../lecture_notes/11_optimisation.md#s-11-5); [Proposition 11.14](../lecture_notes/11_optimisation.md#r-11-14), [Corollary 11.7](../lecture_notes/11_optimisation.md#r-11-7), [Proposition 11.9](../lecture_notes/11_optimisation.md#r-11-9), [Proposition 11.15](../lecture_notes/11_optimisation.md#r-11-15), [Corollary 11.16](../lecture_notes/11_optimisation.md#r-11-16).

`optimize.extremize(objective, params, maximize)` runs Adam in `rounds` rounds of `steps` steps (learning rate
`lr * decay**round`) on the loss $`\mp\sum_b\text{objective}_b`$, re-centres the unitary charts after every round,
and returns **per start** the largest (smallest, for `maximize=False`) value of the objective among those computed at
the gradient iterates and once at the final point. Every batch element is an independent start with its own slots,
signs, $`U_0`$ and contrast: Adam acts coordinatewise, so if the $`b`$-th entry of the objective depends only on the
$`b`$-th parameters (as in all objectives below), the starts do not interact. The signs of `SelfAdjointParam` are
discrete: between rounds, `discrete_step` tries eight single sign flips per start (each at a uniformly random slot),
one after the other, and keeps each flip that improves the objective.

In exact arithmetic every iterate is admissible (Proposition 11.2 and the remarks after it), so every returned value
is the value of $`z`$ at an admissible point; in floating point it is that value up to rounding (the asserts below
allow $`10^{-12}`$). Hence a maximisation gives a **lower bound** for the supremum, a minimisation an upper bound for
the infimum, and nothing more — a local search can stop at a local extremum. The other direction needs a proof,
like Corollary 11.7 and Proposition 11.9, or a convex relaxation (Section 11.6). A computed envelope is checked
(i) against closed forms; (ii) against configurations built by hand, such as direct sums of smaller solutions, which
it must not fall below; (iii) against structured samples (two-level spectra, permutation bases), whose values are
further attained values.

By [lecture notes, Proposition 11.14](../lecture_notes/11_optimisation.md#r-11-14), if all $`k_c\ge2`$, the lower envelope of Exp. 1 is
$`\underline z(\delta)=-\big(D-2\delta/(1+\delta)\big)`$.

The code checks Proposition 11.14 on $`M_2\oplus M_2`$, Corollary 11.7 on $`M_3\oplus M_3`$ (both contrasts prescribed,
$`k=3`$ odd) and Proposition 11.9; suprema and infima share one call, as above. The one-sided asserts (no returned
value beyond a proved extreme) are consequences of these theorems; the asserts that the search comes within the
stated tolerance of the extreme are numerical observations of this seeded run. The keyword `is_self_adjoint` selects
the algorithm for the trace norm (singular values, or eigenvalues for $`YXY`$) without a test at every step. The last
check is of type (ii): $`z`$ of Exp. 3 is additive over direct sums, and the contrast of a direct sum is the largest
contrast of its summands (Exercise 3), so, for the envelopes of Exp. 3 with $`\Delta(Y)`$ free in $`[0,0.999]`$ as in
the code,
```math
\overline z_{M_2\oplus M_2}(\delta)\ \ge\ \overline z_{M_2}(\delta)+\sup_{\delta'\le\delta}\overline z_{M_2}(\delta');
```
in particular the sum of two values attained on $`M_2`$ at contrasts $`\delta`$ and $`\delta'\le\delta`$ is a value
attained on $`M_2\oplus M_2`$ at contrast $`\delta`$.

```python
torch.manual_seed(1)
S = 6                                                                       # starts per grid point
grid = torch.tensor([0.0, 0.5, 0.9], dtype=torch.float64)
alg22 = TypeIAlgebra([2, 2], [2, 2], precision='double')                    # D = 4

# Exp. 1, lower envelope along Delta(X); Delta(Y) free in [0, 0.999]
X1, Y1, U1 = (opt.PositiveParam(alg22, 3 * S, delta=grid.repeat_interleave(S)),
              opt.PositiveParam(alg22, 3 * S), opt.UnitaryParam(alg22, 3 * S))
def z1():
    x, y = X1.operator(), Y1.operator()
    return (x @ U1.operator() @ y).trace.abs() - (x @ y).trace.real
inf1 = opt.extremize(z1, [X1, Y1, U1], maximize=False, rounds=4, steps=30).reshape(3, S).min(1)[0]
exact1 = -(4 - 2 * grid / (1 + grid))                                       # Proposition 11.14
print("Exp. 1  inf:", [f"{v:.4f}" for v in inf1.tolist()], " exact:", [f"{v:.4f}" for v in exact1.tolist()])
assert (inf1 >= exact1 - 1e-12).all() and (inf1 - exact1).max() < 4e-3

# Exp. 2 on M_3 (+) M_3, Delta(X0) = 0.5, Delta(Y0) = 0.8: sup and inf in one call
sign = torch.tensor([1.0, -1.0], dtype=torch.float64).repeat_interleave(S)
X2, Y2, V2 = (opt.PositiveParam(alg3, 2 * S, delta=torch.full((2 * S,), 0.5, dtype=torch.float64)),
              opt.PositiveParam(alg3, 2 * S, delta=torch.full((2 * S,), 0.8, dtype=torch.float64)),
              opt.UnitaryParam(alg3, 2 * S))
def z2():
    x0, y0 = X2.operator(), Y2.operator()
    M = alg3.operator((x0 @ V2.operator() @ y0).matrix, is_self_adjoint=False)
    return sign * (M.trace_norm() - (x0 @ y0).trace.real)
best2 = opt.extremize(z2, [X2, Y2, V2], maximize=True, rounds=4, steps=30).reshape(2, S).max(1)[0]
sup2, inf2 = best2[0].item(), -best2[1].item()
print(f"Exp. 2  sup {sup2:.6f}, inf {inf2:.6f}, exact +-{envelope:.6f}")
assert sup2 <= envelope + 1e-12 and inf2 >= -envelope - 1e-12               # Corollary 11.7: one-sided ...
assert envelope - sup2 < 1e-3 * envelope and inf2 + envelope < 1e-3 * envelope   # ... and reached

# Exp. 3: maximise sign * z at the contrasts `deltas` of |X|; Delta(Y) free in [0, 0.999]
def exp3(alg_, deltas, sign):
    X3, Y3 = opt.SelfAdjointParam(alg_, len(deltas), delta=deltas), opt.PositiveParam(alg_, len(deltas))
    def z3():
        y = Y3.operator()
        M = alg_.operator((y @ X3.operator() @ y).matrix, is_self_adjoint=True)
        return sign * ((y @ X3.abs_operator() @ y).trace.real - M.trace_norm())
    return opt.extremize(z3, [X3, Y3], maximize=True, rounds=4, steps=30).reshape(2, S).max(1)[0]

best3 = exp3(alg22, torch.full((2 * S,), 0.5, dtype=torch.float64), sign)
sup_M22, inf3 = best3[0].item(), -best3[1].item()
s_M2 = exp3(TypeIAlgebra([2], [2], precision='double'),
            torch.tensor([0.0, 0.5], dtype=torch.float64).repeat_interleave(S), 1.0)
direct_sum = (s_M2[1] + s_M2.max()).item()                                 # delta = 0.5, best delta' <= 0.5
closed = [(1 - dl) * (2 * 0.999 / 1.999) ** 2 for dl in (0.0, 0.5)]        # Corollary 11.16, d = 0.999
print(f"Exp. 3  inf {inf3:.1e};  sup for M_2 at delta = 0, 0.5: {s_M2[0]:.5f}, {s_M2[1]:.5f}"
      f"  (closed form {closed[0]:.5f}, {closed[1]:.5f})")
print(f"        sup for M_2 (+) M_2 at 0.5: search {sup_M22:.5f}, direct sum {direct_sum:.5f}")
assert inf3 >= -1e-12 and inf3 < 1e-6                                       # Proposition 11.9
assert sup_M22 >= direct_sum - 5e-3 * direct_sum                           # (ii), up to 0.5 %
assert (s_M2 <= torch.tensor(closed, dtype=torch.float64) + 1e-12).all()   # Corollary 11.16: one-sided
assert sup_M22 <= (2 - 0.5) * (2 * 0.999 / 1.999) ** 2 + 1e-12
```

*Numerical observations* (the printed values of this run). The search on $`M_2\oplus M_2`$ falls short of the direct
sum by about $`0.2\,\%`$: it has not fully converged, and the best lower bound for $`\overline z_{M_2\oplus M_2}(0.5)`$ is
the larger number — detecting this is the purpose of check (ii). The values for $`M_2`$ agree to $`10^{-4}`$ with
$`(1-\delta)\big(2d/(1+d)\big)^2=(1-\delta)(1-\ell(d))^2`$, where $`d=0.999`$ is the largest admitted contrast of $`Y`$. This
closed form, first observed numerically in the companion paper, is the upper envelope of
[Corollary 11.16](../lecture_notes/11_optimisation.md#r-11-16) (below); the code asserts only the proved direction (no
value above it). By the same corollary the direct sum is optimal: $`\overline z_{M_2\oplus M_2}(0.5)=1.5\,(2d/(1+d))^2`$.
For larger blocks we know no closed form for the upper envelopes of Exp. 1 and Exp. 3; computed curves are then lower
bounds (attained values), nothing more.

By [lecture notes, Proposition 11.15](../lecture_notes/11_optimisation.md#r-11-15), on $`M_2`$, for spectra
$`x_1\ge x_2`$ of $`X`$ (of $`|X|`$ in Exp. 3) and $`y_1\ge y_2`$ of $`Y`$, $`\sup z`$ is
$`(x_1-x_2)^2(y_1-y_2)^2/(4(x_1+x_2)(y_1+y_2))`$ in Exp. 1 and $`\frac{2x_1x_2}{x_1+x_2}(y_1-y_2)^2`$ in Exp. 3. With
$`e_1`$ an eigenvector of $`X`$ for its eigenvalue of modulus $`x_1`$ and $`f_1`$ one of $`Y`$ for $`y_1`$, the supremum is
attained when $`\gamma=|\langle e_1,f_1\rangle|^2`$ is $`\gamma_1=\frac14+\frac S{2q}`$ (Exp. 1, with
$`|\mathrm{Tr}(XUY)|=\|YX\|_1`$) or $`\gamma_3=\frac Sq`$ (Exp. 3, $`X=\mathrm{diag}(x_1,-x_2)`$), where $`S=x_1y_2+x_2y_1`$
and $`q=(x_1+x_2)(y_1+y_2)`$. By
[Corollary 11.16](../lecture_notes/11_optimisation.md#r-11-16), for $`C`$ blocks of size $`2`$ with $`\Delta(|X|)=\delta`$ and
$`\Delta(Y)=\delta'`$, Exp. 1 has $`\sup z=C\delta^2\delta'^2/((1+\delta)(1+\delta'))`$ and
$`\inf z=-(2C-1+\ell(\delta)\ell(\delta'))`$, and Exp. 3 has $`\sup z=(C-\delta)\big(2\delta'/(1+\delta')\big)^2`$. Both
suprema increase with $`\delta'`$: with $`\Delta(Y)`$ free in $`[0,d]`$ the upper envelopes are the values at $`\delta'=d`$,
and as $`d\to1`$ they tend to $`C\delta^2/(2(1+\delta))`$ and $`C-\delta`$, the values for a singular $`Y`$
($`\Delta(Y)=1`$).

The code builds the real maximisers of the proof for $`C=1,2`$ (`from_blocks`), checks the contrasts and the values,
checks the bounds of Proposition 11.15 on random complex pairs, attains the limits $`d\to1`$ with rank-one projections
in every block, and runs a multistart search for Exp. 1 on $`M_2`$ along $`\Delta(X)`$ with $`\Delta(Y)`$ free in
$`[0,0.999]`$, which must stay below the closed form (a theorem) and, in this run, comes within $`10^{-3}`$ of it.

```python
def rot(gamma):                                                             # real basis f with |<e_1, f_1>|^2 = gamma
    t = math.acos(math.sqrt(gamma))
    return torch.tensor([[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]], dtype=torch.float64)

def gammas(x1, x2, y1, y2):                                                 # gamma_1 (Exp. 1), gamma_3 (Exp. 3)
    S, q = x1 * y2 + x2 * y1, (x1 + x2) * (y1 + y2)
    return 0.25 + S / (2 * q), S / q

env1 = lambda C, dl, dp: C * dl ** 2 * dp ** 2 / ((1 + dl) * (1 + dp))     # Corollary 11.16(a)
env3 = lambda C, dl, dp: (C - dl) * (2 * dp / (1 + dp)) ** 2               # Corollary 11.16(b)

def maximisers(C, dl, dp):
    """Exp. 1: (X, U, Y) and Exp. 3: (X, Y) attaining Corollary 11.16 on C blocks of size 2."""
    a2 = TypeIAlgebra([2] * C, [2] * C, complex_valued=False, precision='double')
    lx, ly = ell(dl), ell(dp)
    R = rot(gammas(1.0, lx, 1.0, ly)[0])
    X0 = torch.diag(torch.tensor([1.0, lx], dtype=torch.float64))
    Y0 = R @ torch.diag(torch.tensor([1.0, ly], dtype=torch.float64)) @ R.T
    W, _, Vh = torch.linalg.svd(Y0 @ X0)                                    # Y X = W Sigma V^T, U = V W^T
    exp1 = [a2.from_blocks([m] * C) for m in (X0, Vh.T @ W.T, Y0)]
    xs = [lx] + [1.0] * (C - 1)                                             # |X|: {1, ell} in block 1, {1, 1} else
    Rs = [rot(gammas(1.0, x2, 1.0, ly)[1]) for x2 in xs]
    X3 = a2.from_blocks([torch.diag(torch.tensor([1.0, -x2], dtype=torch.float64)) for x2 in xs])
    Y3 = a2.from_blocks([Rc @ torch.diag(torch.tensor([1.0, ly], dtype=torch.float64)) @ Rc.T for Rc in Rs])
    return exp1, (X3, Y3)

print(" C  delta  delta'   Exp. 1 z / formula      Exp. 3 z / formula")
for C in (1, 2):
    for dl, dp in ((0.5, 0.8), (0.9, 0.3), (0.4, 1.0)):
        (X, U, Y), (X3, Y3) = maximisers(C, dl, dp)
        assert abs(X.michelson_contrast.item() - dl) < 1e-12 and abs(Y.michelson_contrast.item() - dp) < 1e-12
        assert abs(X3.abs().michelson_contrast.item() - dl) < 1e-12 and abs(Y3.michelson_contrast.item() - dp) < 1e-12
        z1v = ((X @ U @ Y).trace.abs() - (X @ Y).trace.real).item()
        z3v = ((Y3 @ X3.abs() @ Y3).trace.real - (Y3 @ X3 @ Y3).trace_norm()).item()
        print(f" {C}   {dl:.1f}    {dp:.1f}    {z1v:.6f} / {env1(C, dl, dp):.6f}    {z3v:.6f} / {env3(C, dl, dp):.6f}")
        assert abs(z1v - env1(C, dl, dp)) < 1e-12 and abs(z3v - env3(C, dl, dp)) < 1e-12
print("limits d -> 1 at delta = 0.4:", env1(1, 0.4, 1.0), 0.4 ** 2 / (2 * 1.4), env3(1, 0.4, 1.0), 1 - 0.4)

# Proposition 11.15 on random complex pairs in M_2: no pair exceeds the bounds
M2c = TypeIAlgebra([2], [2], precision='double')
x = torch.rand(4000, 2, dtype=torch.float64).sort(-1, descending=True)[0]
y = torch.rand(4000, 2, dtype=torch.float64).sort(-1, descending=True)[0]
V, Wc = M2c.random_unitary_operator(4000).matrix, M2c.random_unitary_operator(4000).matrix
Xr = M2c.operator(V @ torch.diag_embed(x.to(V.dtype)).unsqueeze(1) @ V.mH, is_positive=True)
Xi = M2c.operator(V @ torch.diag_embed((x * torch.tensor([1.0, -1.0])).to(V.dtype)).unsqueeze(1) @ V.mH,
                  is_self_adjoint=True)
Yr = M2c.operator(Wc @ torch.diag_embed(y.to(V.dtype)).unsqueeze(1) @ Wc.mH, is_positive=True)
(x1, x2), (y1, y2) = x.unbind(-1), y.unbind(-1)
sup1 = (Yr @ Xr).trace_norm() - (Xr @ Yr).trace.real                        # sup over U (Theorem 1.9)
z3r = (Yr @ Xi.abs() @ Yr).trace.real - M2c.operator((Yr @ Xi @ Yr).matrix, is_self_adjoint=True).trace_norm()
assert (sup1 <= (x1 - x2) ** 2 * (y1 - y2) ** 2 / (4 * (x1 + x2) * (y1 + y2)) + 1e-12).all()
assert (z3r <= 2 * x1 * x2 / (x1 + x2) * (y1 - y2) ** 2 + 1e-12).all()

# multistart search for Exp. 1 on M_2 along Delta(X), Delta(Y) free in [0, 0.999]
torch.manual_seed(2)
M2r = TypeIAlgebra([2], [2], complex_valued=False, precision='double')
gx = torch.tensor([0.3, 0.6, 0.9], dtype=torch.float64)
Xp, Yp1, Up = opt.PositiveParam(M2r, 3 * S, delta=gx.repeat_interleave(S)), opt.PositiveParam(M2r, 3 * S), \
    opt.UnitaryParam(M2r, 3 * S)
def z1m():
    xm, ym = Xp.operator(), Yp1.operator()
    return (xm @ Up.operator() @ ym).trace.abs() - (xm @ ym).trace.real
sup1m = opt.extremize(z1m, [Xp, Yp1, Up], maximize=True, rounds=4, steps=30).reshape(3, S).max(1)[0]
exact1m = env1(1, gx, 0.999)                                                # Corollary 11.16(c)
print("Exp. 1 sup on M_2:", [f"{v:.5f}" for v in sup1m.tolist()], " closed form:", [f"{v:.5f}" for v in exact1m.tolist()])
assert (sup1m <= exact1m + 1e-12).all() and ((exact1m - sup1m) / exact1m).max() < 1e-3
```

<a id="s-11-6"></a>

## 11.6 What to compare with

**Manopt**, **pymanopt** and **geoopt** provide matrix manifolds (Stiefel manifolds, orthogonal and unitary groups)
with retractions such as QR, polar or Cayley, and Riemannian solvers (conjugate gradients, trust regions, Riemannian
Adam), an alternative to Adam in a re-centred chart (we make no performance comparison here). **cvxpy** passes
convex problems in operator variables to SDP solvers, which return an optimum together with a dual certificate,
both up to the solver tolerance; a relaxation — a convex problem over a larger feasible set — bounds the supremum
from above, the direction a local search cannot give. **NCTSSOS**, **ncpol2sdpa** and **NCSOStools** bound eigenvalues and traces of noncommutative
polynomials in all dimensions at once (Section 0.3). The library adds exact parametrisations of the spectral and
contrast constraints of $`\bigoplus_cM_{k_c}`$ (the extremes may sit in any sector), unitaries that respect the
blocks, and one batched computation for many starts and constraint values, on a CPU or a GPU. It produces candidate
extremisers and lower bounds for suprema; proofs or relaxations must supply the rest.

<a id="x-11"></a>

## Exercises

[Exercises 1 and 2](../lecture_notes/11_optimisation.md#x-11) of the lecture notes, numerically: the Cayley transform
$`C(S)=(1-S)(1+S)^{-1}`$ of a skew-adjoint $`S`$ is unitary, has no eigenvalue $`-1`$, and $`S=(1-C)(1+C)^{-1}`$
(Exercise 1); for $`O`$ Haar-distributed on $`O(k)`$, sampled by `random_unitary_operator` of a real algebra
(`complex_valued=False`), $`\mathbb E\mathrm{Tr}(AO)^2=\|A\|_2^2/k`$ (Exercise 2).

<details>
<summary>Solutions</summary>

```python
Sk = torch.randn(4, 4, dtype=torch.complex128)
Sk = Sk - Sk.mH                                                   # skew-adjoint
I4 = torch.eye(4, dtype=torch.complex128)
Cay = (I4 - Sk) @ torch.linalg.inv(I4 + Sk)
assert torch.allclose(Cay.mH @ Cay, I4, atol=1e-12)
assert torch.allclose((I4 - Cay) @ torch.linalg.inv(I4 + Cay), Sk, atol=1e-10)
assert torch.linalg.eigvals(Cay).add(1).abs().min() > 1e-6        # -1 is never an eigenvalue
```

```python
Mr = TypeIAlgebra([5], [5], complex_valued=False, precision='double')
Ar = Mr.operator(torch.randn(1, 1, 5, 5, dtype=torch.float64))
tr2 = (Ar @ Mr.random_unitary_operator(20_000)).trace ** 2
assert abs(tr2.mean().item() - Ar.frobenius_norm().item() ** 2 / 5) < 5 * tr2.std().item() / math.sqrt(20_000)
```

</details>

**Next:** [Tutorial 12](12_equivariant_learning.md): equivariant layers as elements of a commutant, in the library.
