**English** · [Español](es/03_traces_and_states.md) · [Français](fr/03_traces_and_states.md) · [Deutsch](de/03_traces_and_states.md) · [中文](zh/03_traces_and_states.md) · [日本語](ja/03_traces_and_states.md) · [Русский](ru/03_traces_and_states.md)

# 3. Traces, states and entropy

On $M_n(\mathbb C)$ there is essentially one trace. On $\bigoplus_cM_{k_c}$ there is a
$C$-parameter family, and the choice matters: densities, entropies and the duals of maps
all depend on it. This chapter classifies the traces, represents states by densities,
proves the basic entropy inequalities, and ends with a plausible trace inequality that
is false.

**You will learn** the three traces of the library and how to convert between them,
the classical part of a state, the decomposition of the von Neumann entropy over
sectors, Klein's inequality and its consequences, and the Michelson contrast.

Throughout, $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ with all $k_c\ge1$,
$D=\sum_ck_c$, $P_c$ the minimal central projections, and $A=\bigoplus_cA_c$.

## 3.1 Traces

A **trace** on $\mathcal M$ is a linear functional $\tau$ with $\tau(AB)=\tau(BA)$ for all
$A,B$. It is **positive** if $\tau(A^*A)\ge0$ for all $A$, and **faithful** if moreover
$\tau(A^*A)=0$ implies $A=0$.

**Theorem 3.1.** Every trace on $\mathcal M$ has the form
$\operatorname{tr}_w(A)=\sum_cw_c\operatorname{Tr}A_c$ with $w\in\mathbb C^C$. It is positive
iff all $w_c\ge0$ and faithful iff all $w_c>0$.

*Proof.* Let $E^c_{ij}$ be the matrix units of the $c$-th block. For a trace $\tau$,
$\tau(E^c_{ij})=\tau(E^c_{i1}E^c_{1j})=\tau(E^c_{1j}E^c_{i1})=\delta_{ij}\tau(E^c_{11})$.
Put $w_c=\tau(E^c_{11})$; by linearity $\tau(A)=\sum_c\sum_i(A_c)_{ii}w_c$. Conversely
every $\operatorname{tr}_w$ is a trace. Positivity: $E^c_{11}=(E^c_{11})^*E^c_{11}$, so
$w_c\ge0$ is necessary, and $\operatorname{tr}_w(A^*A)=\sum_cw_c\|A_c\|_2^2$ shows sufficiency
and the faithfulness criterion. $\square$

The library implements three faithful traces:

| name | weights $w_c$ | value on $1$ | meaning |
|---|---|---|---|
| `Tr_blunt` | $1$ | $D$ | trace of $\bigoplus_c\mathbb C^{k_c}$ (multiplicities $m_c=1$) |
| `Tr_norm` | $1/k_c$ | $C$ | every block normalised |
| `tau_vN` | $1/(Ck_c)$ | $1$ | the tracial state: every sector has weight $1/C$ |

By Chapter 2, the trace of a space on which $\mathcal M$ acts with multiplicities $m_c$
is $\operatorname{tr}_w$ with $w_c=m_c$. A faithful trace is thus a choice of "how much each
sector counts".

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

## 3.2 States and densities

A **state** is a linear functional $\omega$ on $\mathcal M$ that is positive
($\omega(A^*A)\ge0$) and normalised ($\omega(1)=1$).

**Theorem 3.2.** Let $\operatorname{tr}_w$ be a faithful trace. For every linear functional
$\omega$ there is a unique $\rho_w\in\mathcal M$ with $\omega(A)=\operatorname{tr}_w(\rho_wA)$ for all
$A$. The functional $\omega$ is positive iff $\rho_w\ge0$, and a state iff moreover
$\operatorname{tr}_w(\rho_w)=1$. For two faithful traces, the densities are related blockwise by
$$
w_c\,(\rho_w)_c=w'_c\,(\rho_{w'})_c .
$$

*Proof.* $\langle X,Y\rangle_w=\operatorname{tr}_w(X^*Y)$ is an inner product on $\mathcal M$
(Theorem 3.1). By the Riesz representation theorem in the finite-dimensional space
$\mathcal M$, there is a unique $R$ with $\omega(A)=\langle R,A\rangle_w=\operatorname{tr}_w(R^*A)$; set
$\rho_w=R^*$. If $\rho_w\ge0$, then
$\omega(A^*A)=\sum_cw_c\operatorname{Tr}(A_c\,(\rho_w)_c\,A_c^*)\ge0$, since
$A_c(\rho_w)_cA_c^*\ge0$ (Corollary 1.7 with $C=(\rho_w)_c^{1/2}A_c^*$). Conversely, if
$\omega$ is positive, take $x\in\mathbb C^{k_c}$ and $A$ the rank-one operator
$|x\rangle\langle x|$ in block $c$; then $A=A^*A/\|x\|^2$ and
$0\le\omega(A)=w_c\langle x,(\rho_w)_cx\rangle$, so $(\rho_w)_c\ge0$. The normalisation is
$\omega(1)=\operatorname{tr}_w(\rho_w)$. Finally both $\sum_cw_c\operatorname{Tr}((\rho_w)_cA_c)$ and
$\sum_cw'_c\operatorname{Tr}((\rho_{w'})_cA_c)$ equal $\omega(A)$ for all $A$, and by uniqueness
$w_c(\rho_w)_c=w'_c(\rho_{w'})_c$. $\square$

The library stores a state by its density $\rho=\rho_{\mathrm{blunt}}$ with respect to
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

> **Physical meaning.** $\omega(A)$ is the expectation value of the observable $A$.
> The density $\rho$ of the blunt trace is the usual density matrix of quantum
> mechanics on $\bigoplus_c\mathbb C^{k_c}$. Densities with respect to other traces appear
> when the same algebra acts with multiplicities (Chapter 2), for example on the
> symmetric and the non-symmetric part of a spin system.

## 3.3 The classical part of a state

The restriction of $\omega$ to the centre is a probability distribution on the sectors,
$p_c=\omega(P_c)=\operatorname{Tr}\rho_c$. If $p_c>0$, the **conditional state** in sector $c$
has density $\sigma_c=\rho_c/p_c$ on $\mathbb C^{k_c}$, and
$$
\rho=\bigoplus_cp_c\,\sigma_c .
$$
A state on $\mathcal M$ is a classical mixture of quantum states in the sectors. Measuring
"which sector" does not disturb it: the measurement channel
$\rho\mapsto\sum_cP_c\rho P_c$ is the identity on $\mathcal M$, because every element of
$\mathcal M$ is block diagonal. This is the operational content of a superselection rule.

```python
p = rho.sector_probabilities()                                   # shape (batch, C)
assert torch.allclose(p.sum(-1), torch.ones(4, dtype=p.dtype))
P = [alg.central([1.0 if c == d else 0.0 for d in range(alg.C)]) for c in range(alg.C)]
measure_sector = channels.measurement_channel(P)
assert (measure_sector(rho).matrix - rho.matrix).abs().max() < 1e-12   # no disturbance
```

## 3.4 Entropy

The **von Neumann entropy** of a state with (blunt) density $\rho$ is
$S(\rho)=-\operatorname{Tr}\rho\log\rho=-\sum_i\lambda_i\log\lambda_i$ over the eigenvalues of
$\rho$, with $0\log0=0$.

**Theorem 3.3 (entropy of a mixture of sectors).** With $p_c$ and $\sigma_c$ as above
(sectors with $p_c=0$ omitted),
$$
S(\rho)=H(p)+\sum_cp_c\,S(\sigma_c),\qquad H(p)=-\sum_cp_c\log p_c .
$$

*Proof.* The eigenvalues of $\rho$ are the numbers $p_c\mu$ for the eigenvalues $\mu$ of
the $\sigma_c$. Hence
$-\sum_c\sum_\mu p_c\mu\log(p_c\mu)=-\sum_cp_c\log p_c\sum_\mu\mu-\sum_cp_c\sum_\mu\mu\log\mu$,
and $\sum_\mu\mu=\operatorname{Tr}\sigma_c=1$. $\square$

The entropy relative to another trace is $S_w(\rho)=-\operatorname{tr}_w(\rho_w\log\rho_w)$. By
Theorem 3.2, $(\rho_w)_c=\rho_c/w_c$, so
$$
S_w(\rho)=-\sum_cw_c\operatorname{Tr}\Big(\frac{\rho_c}{w_c}\log\frac{\rho_c}{w_c}\Big)
=S(\rho)+\sum_cp_c\log w_c .
$$
The entropy is a property of the state *and* the trace: the same state has different
entropies with respect to `Tr_blunt` and `tau_vN`.

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

## 3.5 Relative entropy and Klein's inequality

For states with densities $\rho,\sigma$, $\sigma>0$, the **relative entropy** is
$S(\rho\|\sigma)=\operatorname{Tr}\rho(\log\rho-\log\sigma)$. It does not depend on the trace:
$\log(\rho_w)_c-\log(\sigma_w)_c=\log\rho_c-\log\sigma_c$ because the factor $1/w_c$ cancels,
and $\operatorname{tr}_w(\rho_wX)=\operatorname{Tr}(\rho X)$.

**Theorem 3.4 (Klein's inequality).** Let $f$ be convex and differentiable on an interval
$I$ and $X,Y$ self-adjoint with spectra in $I$. Then
$\operatorname{Tr}[f(X)-f(Y)-(X-Y)f'(Y)]\ge0$. If $f$ is strictly convex, equality holds only
for $X=Y$.

*Proof.* Let $X=\sum_ia_i|u_i\rangle\langle u_i|$ and $Y=\sum_jb_j|v_j\rangle\langle v_j|$ with
orthonormal eigenbases, and $c_{ij}=|\langle u_i,v_j\rangle|^2$, so $\sum_ic_{ij}=\sum_jc_{ij}=1$.
Then $\operatorname{Tr}f(X)=\sum_{ij}c_{ij}f(a_i)$, $\operatorname{Tr}f(Y)=\sum_{ij}c_{ij}f(b_j)$ and
$\operatorname{Tr}Xf'(Y)=\sum_{ij}c_{ij}a_if'(b_j)$, $\operatorname{Tr}Yf'(Y)=\sum_{ij}c_{ij}b_jf'(b_j)$.
Hence the left side equals $\sum_{ij}c_{ij}[f(a_i)-f(b_j)-(a_i-b_j)f'(b_j)]$, a sum of
nonnegative terms by convexity. If $f$ is strictly convex and the sum is zero, then
$a_i=b_j$ whenever $c_{ij}>0$, so $Xv_j=\sum_ia_i\langle u_i,v_j\rangle u_i=b_jv_j$ for every $j$,
and $X=Y$. $\square$

**Corollary 3.5.** $S(\rho\|\sigma)\ge0$, with equality iff $\rho=\sigma$. In particular
$S(\rho)\le\log D$, with equality iff $\rho=1/D$.

*Proof.* Apply Theorem 3.4 to $f(t)=t\log t$ (strictly convex on $[0,\infty)$,
$f'(t)=\log t+1$; at $t=0$ only $f$ itself enters, and the computation holds by continuity)
with $X=\rho$, $Y=\sigma$: the left side is
$\operatorname{Tr}\rho\log\rho-\operatorname{Tr}\sigma\log\sigma-\operatorname{Tr}(\rho-\sigma)(\log\sigma+1)
=S(\rho\|\sigma)$, since $\operatorname{Tr}\rho=\operatorname{Tr}\sigma=1$. For the second claim, take
$\sigma=1/D$: $0\le S(\rho\|1/D)=-S(\rho)+\log D$. $\square$

**Theorem 3.6 (chain rule; monotonicity under the centre).** Let $q_c=\operatorname{Tr}\sigma_c$
and $\hat\rho_c=\rho_c/p_c$, $\hat\sigma_c=\sigma_c/q_c$. Then
$$
S(\rho\|\sigma)=\sum_cp_c\log\frac{p_c}{q_c}+\sum_cp_c\,S(\hat\rho_c\|\hat\sigma_c)
\ \ge\ \sum_cp_c\log\frac{p_c}{q_c}=S(E\rho\|E\sigma),
$$
where $E$ is the conditional expectation onto the centre (Proposition 2.7).

*Proof.* Blockwise, $\log\rho_c=\log p_c+\log\hat\rho_c$ and similarly for $\sigma$, so
$\operatorname{Tr}\rho_c(\log\rho_c-\log\sigma_c)=p_c\log(p_c/q_c)+p_cS(\hat\rho_c\|\hat\sigma_c)$.
Sum over $c$; the second sum is $\ge0$ by Corollary 3.5 in each block. Finally
$E\rho=\bigoplus_c(p_c/k_c)1_{k_c}$, $E\sigma=\bigoplus_c(q_c/k_c)1_{k_c}$, and
$\operatorname{Tr}E\rho(\log E\rho-\log E\sigma)=\sum_ck_c\frac{p_c}{k_c}\log\frac{p_c}{q_c}$. $\square$

So forgetting the quantum information inside the sectors can only decrease the
distinguishability of two states; what remains is the classical Kullback–Leibler
divergence of the sector distributions.

```python
sigma_state = states.random_density_matrix(alg, batch_size=4)
q = sigma_state.sector_probabilities()
D_rel = states.relative_entropy(rho, sigma_state)
kl = (p * torch.log(p / q)).sum(-1)
E = channels.center_expectation(alg)
assert (D_rel >= kl - 1e-12).all()                                       # Theorem 3.6
assert torch.allclose(states.relative_entropy(E(rho), E(sigma_state)), kl)
assert (rho.von_neumann_entropy() <= math.log(10) + 1e-12).all()         # Corollary 3.5
assert abs(states.maximally_mixed_state(alg).von_neumann_entropy().item() - math.log(10)) < 1e-12
```

## 3.6 The Michelson contrast

For $A\ge0$, $A\neq0$, the **Michelson contrast** is
$\Delta(A)=(\lambda_{\max}-\lambda_{\min})/(\lambda_{\max}+\lambda_{\min})\in[0,1]$, with the
extreme eigenvalues taken over all sectors (and $\Delta(0)=0$). It measures how far $A$ is
from a multiple of the identity: $\Delta(A)=0$ iff $A=\lambda1$.

**Proposition 3.7.** If $\Phi:\mathcal M\to\mathcal M$ is linear, positive and unital
($\Phi(1)=1$), then $\Delta(\Phi(A))\le\Delta(A)$ for all $A\ge0$.

*Proof.* $\lambda_{\min}1\le A\le\lambda_{\max}1$, and positivity and unitality give
$\lambda_{\min}1\le\Phi(A)\le\lambda_{\max}1$. Hence
$\lambda_{\min}(\Phi(A))\ge\lambda_{\min}$ and $\lambda_{\max}(\Phi(A))\le\lambda_{\max}$. If
$\Phi(A)=0$ the claim is trivial. Otherwise, with $r=\lambda_{\min}/\lambda_{\max}$,
$\Delta=(1-r)/(1+r)$ is decreasing in $r$, and $r(\Phi(A))\ge r(A)$. $\square$

```python
A_pos = alg.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(4, d), batch_size=4, force_positive=True)
Phi = channels.random_mixed_unitary_channel(alg, 3, batch_size=4)       # unital
assert (Phi(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
assert (E(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
```

## 3.7 A plausible inequality that is false

For positive $X,Y$ in $M_k$ and a unitary $U$, is $|\operatorname{Tr}(XUY)|\le\operatorname{Tr}(XY)$? It
holds for $U=1$ and looks like a Cauchy–Schwarz-type statement. It is false whenever $X$
and $Y$ do not commute.

**Proposition 3.8.** $\max_U|\operatorname{Tr}(XUY)|=\|YX\|_1\ge\operatorname{Tr}(XY)$, with equality iff
$XY=YX$.

*Proof.* $\operatorname{Tr}(XUY)=\operatorname{Tr}(YX\,U)$, so the maximum is $\|YX\|_1$ by Theorem 1.9,
and $\|YX\|_1\ge|\operatorname{Tr}(YX)|=\operatorname{Tr}(XY)$. If $XY=YX$, then $YX\ge0$ (it equals
$X^{1/2}YX^{1/2}$, as $X^{1/2}$ commutes with $Y$ by Theorem 1.5(d)), so
$\|YX\|_1=\operatorname{Tr}(YX)$. Conversely suppose $\|YX\|_1=\operatorname{Tr}(YX)$ and write
$YX=W|YX|$. By Theorem 1.9, $\operatorname{Tr}(YX\cdot1)=\|YX\|_1$ means equality in the
Cauchy–Schwarz step for $U=1$:
$|YX|^{1/2}$ and $W|YX|^{1/2}$ are proportional with a nonnegative factor, which by
comparing norms is $1$, so $W|YX|^{1/2}=|YX|^{1/2}$ and $YX=W|YX|=|YX|\ge0$. A positive
matrix is self-adjoint, so $YX=(YX)^*=XY$. $\square$

The library finds the maximiser with the polar decomposition; Monte Carlo sampling of
$U$ almost never comes close to it for large $k$ — a central theme of the companion paper.

```python
Mk = TypeIAlgebra([6], [6], precision='double')
X = Mk.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
Y = Mk.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
YX = Y @ X
W = YX @ YX.abs().inverse()                                          # polar part of YX
gap = (X @ W.adjoint() @ Y).trace.abs() - (X @ Y).trace.real
assert abs(gap.item() - (YX.trace_norm() - (X @ Y).trace.real).item()) < 1e-10
assert gap.item() > 0                                                # the inequality fails
samples = Mk.random_unitary_operator(20_000)
best = (X @ samples @ Y).trace.abs().max() - (X @ Y).trace.real
print(f"exact sup {gap.item():.4f}, best of 20000 Haar samples {best.item():.4f}")
assert best.item() < gap.item()
```

## 3.8 What to compare with

QuTiP (`entropy_vn`, `entropy_relative`) and the quantum-information packages (e.g.
`qutip`, `toqito`) compute entropies and relative entropies of density matrices on a
single Hilbert space; on $\bigoplus_c\mathbb C^{k_c}$ they agree with the blunt-trace
conventions here. What they do not provide is the family of traces of Theorem 3.1 and
the conversions of Theorem 3.2, which become necessary as soon as an algebra acts with
multiplicities or the sectors must be weighted differently.

## Exercises

1. Show that the density of the tracial state $\tau_{\mathrm{vN}}$ with respect to
   $\tau_{\mathrm{vN}}$ is $1$, and compute its blunt density and its entropy $S$.
2. Show that $S_w(\rho)\le\log\operatorname{tr}_w(1)$, with equality iff
   $\rho_w=1/\operatorname{tr}_w(1)$.
3. Show that the measurement channel $\rho\mapsto\sum_cP_c\rho P_c$ on $M_D$ (not on
   $\mathcal M$!) does disturb states, and compute $S(\rho)$ before and after for a
   superposition of two sectors.

<details>
<summary>Solutions</summary>

1. $\tau_{\mathrm{vN}}(A)=\tau_{\mathrm{vN}}(1\cdot A)$, so $\rho_{\tau}=1$ by uniqueness. By
   Theorem 3.2 the blunt density is $\bigoplus_c1_{k_c}/(Ck_c)$, with $p_c=1/C$ and
   $\sigma_c=1/k_c$, so $S=\log C+\frac1C\sum_c\log k_c$ by Theorem 3.3.

```python
tau = states.tracial_state(alg)
assert (tau.density('tau_vN') - one.matrix).abs().max() < 1e-12
assert abs(tau.von_neumann_entropy().item() - (math.log(3) + sum(math.log(k) for k in (2, 3, 5)) / 3)) < 1e-12
```

2. $S_w(\rho)=S(\rho)+\sum_cp_c\log w_c=-\sum_cp_c\log p_c+\sum_cp_cS(\sigma_c)+\sum_cp_c\log w_c
   \le\sum_cp_c\log\frac{k_cw_c}{p_c}\le\log\sum_ck_cw_c=\log\operatorname{tr}_w(1)$, using
   $S(\sigma_c)\le\log k_c$ (Corollary 3.5) and Jensen's inequality for the concave
   $\log$. Equality requires $\sigma_c=1/k_c$ and $k_cw_c/p_c$ constant, i.e.
   $p_c=k_cw_c/\operatorname{tr}_w(1)$, which means $(\rho_w)_c=\rho_c/w_c=1_{k_c}/\operatorname{tr}_w(1)$.

3. In $M_D$ a pure superposition $\psi=(\psi_1\oplus\psi_2)/\sqrt2$ of unit vectors in two
   sectors has $S=0$; after the measurement the state is
   $\frac12|\psi_1\rangle\langle\psi_1|\oplus\frac12|\psi_2\rangle\langle\psi_2|$ with $S=\log2$. Such a
   superposition is not a state *of $\mathcal M$*: restricted to $\mathcal M$ it already is
   the mixture, which is why superselection rules forbid coherent superpositions of sectors.

</details>

**Next:** Chapter 4 (in preparation) treats symmetries and representations: how a
group action on a physical system produces the sectors, and how the library builds
the corresponding algebras for spin systems and qubit registers.
