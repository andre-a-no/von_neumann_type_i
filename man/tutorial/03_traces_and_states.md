**English** · [Español](es/03_traces_and_states.md) · [Français](fr/03_traces_and_states.md) · [Deutsch](de/03_traces_and_states.md) · [中文](zh/03_traces_and_states.md) · [日本語](ja/03_traces_and_states.md) · [Русский](ru/03_traces_and_states.md)

# Tutorial 3. Traces, states and entropy

This tutorial accompanies [Chapter 3](../lecture_notes/03_traces_and_states.md) of the lecture notes, which classifies the traces
on $`\mathcal M=\bigoplus_cM_{k_c}`$ ([Theorem 3.1](../lecture_notes/03_traces_and_states.md#r-3-1)), represents states by densities
([Theorem 3.2](../lecture_notes/03_traces_and_states.md#r-3-2)) and proves the basic entropy inequalities. We use the three traces of the
library, convert densities between them, and check the theorems numerically on $`M_2\oplus M_3\oplus M_5`$.

<a id="s-3-1"></a>

## 3.1 Traces

*Lecture notes:* [Section 3.1](../lecture_notes/03_traces_and_states.md#s-3-1); by [Theorem 3.1](../lecture_notes/03_traces_and_states.md#r-3-1) every trace is $`\mathrm{tr}_w(A)=\sum_cw_c\mathrm{Tr}\,A_c`$.

The library implements three faithful traces:

| name | weights $`w_c`$ | value on $`1`$ | meaning |
|---|---|---|---|
| `Tr_blunt` | $`1`$ | $`D`$ | trace of $`\bigoplus_c\mathbb C^{k_c}`$ (multiplicities $`m_c=1`$) |
| `Tr_norm` | $`1/k_c`$ | $`C`$ | every block normalised |
| `tau_vN` | $`1/(Ck_c)`$ | $`1`$ | the tracial state with equal sector weights $`1/C`$ |

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, DensityMatrix, states, channels

torch.manual_seed(0)
alg = TypeIAlgebra([2, 3, 5], [2, 3, 5], precision='double')    # M_2 (+) M_3 (+) M_5, D = 10, C = 3
U = alg.random_unitary_operator(4)
A = alg.operator_from_eigenvalues(lambda d: torch.randn(4, d), batch_size=4, force_self_adjoint=True)
B = U @ A                                                       # a batch of four non-normal operators
one = alg.identity()

for f, value in ((lambda X: X.Tr_blunt(), 10), (lambda X: X.Tr_norm(), 3), (lambda X: X.tau_vN(), 1)):
    assert abs(f(one).item() - value) < 1e-12                   # normalisation
    assert (f(A @ B) - f(B @ A)).abs().max() < 1e-10            # trace property
assert (A.Tr_norm() - 3 * A.tau_vN()).abs().max() < 1e-12       # Tr_norm = C tau_vN
```

<a id="s-3-2"></a>

## 3.2 States and densities

*Lecture notes:* [Section 3.2](../lecture_notes/03_traces_and_states.md#s-3-2); [Theorem 3.2](../lecture_notes/03_traces_and_states.md#r-3-2) gives the density $`\rho_w`$ of a state for every faithful trace, with $`w_c(\rho_w)_c=w'_c(\rho_{w'})_c`$.

The library stores a state by its density $`\rho=\rho_{\mathrm{blunt}}`$ with respect to
`Tr_blunt` (class `DensityMatrix`); `density(trace)` and `DensityMatrix.from_density`
convert. The value of the state does not depend on the convention:

```python
rho = states.random_density_matrix(alg, batch_size=4)          # Hilbert-Schmidt random states
for trace, tr in (('blunt', lambda X: X.Tr_blunt()), ('norm', lambda X: X.Tr_norm()),
                  ('tau_vN', lambda X: X.tau_vN())):
    d = alg.operator(rho.density(trace))
    assert (tr(d @ A) - (rho @ A).Tr_blunt()).abs().max() < 1e-10       # omega(A) for all traces
    assert (tr(d) - 1).abs().max() < 1e-10                               # tr_w(rho_w) = 1
    back = DensityMatrix.from_density(alg, rho.density(trace), trace=trace)
    assert (back.matrix - rho.matrix).abs().max() < 1e-12
```

<a id="s-3-3"></a>

## 3.3 The classical part of a state

*Lecture notes:* [Section 3.3](../lecture_notes/03_traces_and_states.md#s-3-3): the sector probabilities $`p_c=\mathrm{Tr}\,\rho_c`$, and the map $`\rho\mapsto\sum_cP_c\rho P_c`$, which is the identity on $`\mathcal M`$.

```python
p = rho.sector_probabilities()                                   # shape (batch, C)
assert torch.allclose(p.sum(-1), torch.ones(4, dtype=p.dtype))
P = [alg.central([1.0 if c == d else 0.0 for d in range(alg.C)]) for c in range(alg.C)]
measure_sector = channels.measurement_channel(P)
assert (measure_sector(rho).matrix - rho.matrix).abs().max() < 1e-12   # no disturbance
```

<a id="s-3-4"></a>

## 3.4 Entropy

*Lecture notes:* [Section 3.4](../lecture_notes/03_traces_and_states.md#s-3-4); [Theorem 3.3](../lecture_notes/03_traces_and_states.md#r-3-3) gives $`S(\rho)=H(p)+\sum_cp_cS(\sigma_c)`$, and the entropy relative to $`\mathrm{tr}_w`$ is $`S_w(\rho)=S(\rho)+\sum_cp_c\log w_c`$.

```python
S = rho.von_neumann_entropy()
sigma = [rho.matrix[:, c, :k, :k] / p[:, c, None, None] for c, k in enumerate(alg.k_factors)]
S_sec = [-(w * torch.log(w)).sum(-1) for w in (torch.linalg.eigvalsh(s) for s in sigma)]
H = -(p * torch.log(p)).sum(-1)
assert torch.allclose(S, H + sum(p[:, c] * S_sec[c] for c in range(alg.C)))    # Theorem 3.3

w_tau = torch.tensor([1 / (alg.C * k) for k in alg.k_factors], dtype=torch.float64)
d_tau = alg.operator(rho.density('tau_vN'), is_positive=True)
S_tau = -(d_tau @ d_tau.log()).tau_vN().real
assert torch.allclose(S_tau, S + (p * torch.log(w_tau)).sum(-1))
```

<a id="s-3-5"></a>

## 3.5 Relative entropy and Klein's inequality

*Lecture notes:* [Section 3.5](../lecture_notes/03_traces_and_states.md#s-3-5); [Corollary 3.6](../lecture_notes/03_traces_and_states.md#r-3-6) ($`S(\rho)\le\log D`$, with equality for $`\rho=1/D`$) and [Theorem 3.7](../lecture_notes/03_traces_and_states.md#r-3-7) ($`S(\rho\|\sigma)\ge S(E\rho\|E\sigma)=\sum_cp_c\log(p_c/q_c)`$ for the conditional expectation $`E`$ onto the centre).

```python
sigma_state = states.random_density_matrix(alg, batch_size=4)
q = sigma_state.sector_probabilities()
D_rel = states.relative_entropy(rho, sigma_state)
kl = (p * torch.log(p / q)).sum(-1)
E = channels.center_expectation(alg)
assert (D_rel >= kl - 1e-12).all()                                       # Theorem 3.7
assert torch.allclose(states.relative_entropy(E(rho), E(sigma_state)), kl)
assert (rho.von_neumann_entropy() <= math.log(10) + 1e-12).all()         # Corollary 3.6
assert abs(states.maximally_mixed_state(alg).von_neumann_entropy().item() - math.log(10)) < 1e-12
```

<a id="s-3-6"></a>

## 3.6 The Michelson contrast

*Lecture notes:* [Section 3.6](../lecture_notes/03_traces_and_states.md#s-3-6); by [Proposition 3.8](../lecture_notes/03_traces_and_states.md#r-3-8) a positive unital map does not increase the Michelson contrast, which the library computes as `Operator.michelson_contrast`.

```python
A_pos = alg.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(4, d), batch_size=4, force_positive=True)
Phi = channels.random_mixed_unitary_channel(alg, 3, batch_size=4)       # unital
assert (Phi(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
assert (E(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
```

<a id="s-3-7"></a>

## 3.7 A plausible inequality that is false

*Lecture notes:* [Section 3.7](../lecture_notes/03_traces_and_states.md#s-3-7); by [Proposition 3.9](../lecture_notes/03_traces_and_states.md#r-3-9), for $`X,Y\ge0`$ in $`M_k`$, $`\max_U|\mathrm{Tr}(XUY)|=\|YX\|_1\ge\mathrm{Tr}(XY)`$, with equality iff $`XY=YX`$.

The library finds the maximiser with the polar decomposition: for $`YX=W|YX|`$ the unitary $`U=W^*`$
attains $`\|YX\|_1`$ ([lecture notes, Theorem 1.9](../lecture_notes/01_matrix_algebras.md#r-1-9)). Random search over
Haar-distributed unitaries (the library samples $`U`$ in this way) is worse on average
([lecture notes, Section 3.7](../lecture_notes/03_traces_and_states.md#s-3-7)); the probability bound $`b^4/(a^4k^2)`$ recalled there says nothing
in the run below, where $`b/a=6`$ and $`k=6`$. The code
computes the exact supremum for one non-commuting pair with $`k=6`$ and the best of $`20000`$ Haar
samples, both minus $`\mathrm{Tr}(XY)`$. By Proposition 3.9 no sample can exceed the supremum; how far
the best sample stays below it is a numerical observation of this run, not a theorem.

```python
Mk = TypeIAlgebra([6], [6], precision='double')
X = Mk.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
Y = Mk.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
assert (X @ Y - Y @ X).matrix.abs().max().item() > 1e-6             # X and Y do not commute
YX = Y @ X
W = YX @ YX.abs().inverse()                                          # polar part of YX
gap = (X @ W.adjoint() @ Y).trace.abs() - (X @ Y).trace.real
assert abs(gap.item() - (YX.trace_norm() - (X @ Y).trace.real).item()) < 1e-10   # Theorem 1.9
assert gap.item() > 0                                  # Proposition 3.9: the inequality fails
samples = Mk.random_unitary_operator(20_000)
best = (X @ samples @ Y).trace.abs().max() - (X @ Y).trace.real
print(f"sup - Tr(XY) = {gap.item():.4f}; numerical observation: best of 20000 Haar samples "
      f"- Tr(XY) = {best.item():.4f}")
assert best.item() <= gap.item() + 1e-10               # Proposition 3.9: no sample beats the sup
```

<a id="s-3-8"></a>

## 3.8 What to compare with

QuTiP (`entropy_vn`, `entropy_relative`) and other quantum-information packages (e.g.
`toqito`) compute entropies and relative entropies of density matrices on a
single Hilbert space, with the trace of that space; on $`\bigoplus_c\mathbb C^{k_c}`$ this is the
blunt trace here (up to the base of the logarithm, which differs between packages). What they do not provide is the family of traces of Theorem 3.1 and
the conversions of Theorem 3.2, which become necessary as soon as an algebra acts with
multiplicities or the sectors must be weighted differently.

<a id="x-3"></a>

## Exercises

[Exercise 1](../lecture_notes/03_traces_and_states.md#x-3) of the lecture notes: the density of the tracial state $`\tau_{\mathrm{vN}}`$
(`states.tracial_state`) with respect to $`\tau_{\mathrm{vN}}`$ is $`1`$, and its entropy is
$`S=\log C+\frac1C\sum_c\log k_c`$. Check both for $`M_2\oplus M_3\oplus M_5`$.

<details>
<summary>Solutions</summary>

```python
tau = states.tracial_state(alg)
assert (tau.density('tau_vN') - one.matrix).abs().max() < 1e-12
assert abs(tau.von_neumann_entropy().item() - (math.log(3) + sum(math.log(k) for k in (2, 3, 5)) / 3)) < 1e-12
```

</details>

**Next:** [Tutorial 4. Symmetries and representations](04_symmetries.md)
