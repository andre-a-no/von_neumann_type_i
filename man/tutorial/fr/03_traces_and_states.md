[English](../03_traces_and_states.md) · [Español](../es/03_traces_and_states.md) · **Français** · [Deutsch](../de/03_traces_and_states.md) · [中文](../zh/03_traces_and_states.md) · [日本語](../ja/03_traces_and_states.md) · [Русский](../ru/03_traces_and_states.md)

# Tutoriel 3. Traces, états et entropie

Ce tutoriel accompagne le [chapitre 3](../../lecture_notes/fr/03_traces_and_states.md) des notes de cours, qui classifie les traces
sur $`\mathcal M=\bigoplus_cM_{k_c}`$ ([théorème 3.1](../../lecture_notes/fr/03_traces_and_states.md#r-3-1)), représente les états par des densités
([théorème 3.2](../../lecture_notes/fr/03_traces_and_states.md#r-3-2)) et démontre les inégalités d'entropie de base. Nous utilisons les trois traces de la
bibliothèque, convertissons les densités de l'une à l'autre, et vérifions numériquement les théorèmes sur $`M_2\oplus M_3\oplus M_5`$.

<a id="s-3-1"></a>

## 3.1 Traces

*Notes de cours :* [section 3.1](../../lecture_notes/fr/03_traces_and_states.md#s-3-1) ; d'après le [théorème 3.1](../../lecture_notes/fr/03_traces_and_states.md#r-3-1), toute trace est de la forme $`\mathrm{tr}_w(A)=\sum_cw_c\mathrm{Tr}\,A_c`$.

La bibliothèque implémente trois traces fidèles :

| nom | poids $`w_c`$ | valeur en $`1`$ | signification |
|---|---|---|---|
| `Tr_blunt` | $`1`$ | $`D`$ | trace de $`\bigoplus_c\mathbb C^{k_c}`$ (multiplicités $`m_c=1`$) |
| `Tr_norm` | $`1/k_c`$ | $`C`$ | chaque bloc normalisé |
| `tau_vN` | $`1/(Ck_c)`$ | $`1`$ | l'état tracial de poids égaux $`1/C`$ sur les secteurs |

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

## 3.2 États et densités

*Notes de cours :* [section 3.2](../../lecture_notes/fr/03_traces_and_states.md#s-3-2) ; le [théorème 3.2](../../lecture_notes/fr/03_traces_and_states.md#r-3-2) donne la densité $`\rho_w`$ d'un état pour toute trace fidèle, avec $`w_c(\rho_w)_c=w'_c(\rho_{w'})_c`$.

La bibliothèque stocke un état par sa densité $`\rho=\rho_{\mathrm{blunt}}`$ relativement à
la trace brute `Tr_blunt` (classe `DensityMatrix`) ; `density(trace)` et
`DensityMatrix.from_density` effectuent les conversions. La valeur de l'état ne dépend pas
de la convention :

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

## 3.3 La partie classique d'un état

*Notes de cours :* [section 3.3](../../lecture_notes/fr/03_traces_and_states.md#s-3-3) : les probabilités des secteurs $`p_c=\mathrm{Tr}\,\rho_c`$, et l'application $`\rho\mapsto\sum_cP_c\rho P_c`$, qui est l'identité sur $`\mathcal M`$.

```python
p = rho.sector_probabilities()                                   # shape (batch, C)
assert torch.allclose(p.sum(-1), torch.ones(4, dtype=p.dtype))
P = [alg.central([1.0 if c == d else 0.0 for d in range(alg.C)]) for c in range(alg.C)]
measure_sector = channels.measurement_channel(P)
assert (measure_sector(rho).matrix - rho.matrix).abs().max() < 1e-12   # no disturbance
```

<a id="s-3-4"></a>

## 3.4 Entropie

*Notes de cours :* [section 3.4](../../lecture_notes/fr/03_traces_and_states.md#s-3-4) ; le [théorème 3.3](../../lecture_notes/fr/03_traces_and_states.md#r-3-3) donne $`S(\rho)=H(p)+\sum_cp_cS(\sigma_c)`$, et l'entropie relativement à $`\mathrm{tr}_w`$ est $`S_w(\rho)=S(\rho)+\sum_cp_c\log w_c`$.

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

## 3.5 Entropie relative et inégalité de Klein

*Notes de cours :* [section 3.5](../../lecture_notes/fr/03_traces_and_states.md#s-3-5) ; [corollaire 3.6](../../lecture_notes/fr/03_traces_and_states.md#r-3-6) ($`S(\rho)\le\log D`$, avec égalité pour $`\rho=1/D`$) et [théorème 3.7](../../lecture_notes/fr/03_traces_and_states.md#r-3-7) ($`S(\rho\|\sigma)\ge S(E\rho\|E\sigma)=\sum_cp_c\log(p_c/q_c)`$ pour l'espérance conditionnelle $`E`$ sur le centre).

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

## 3.6 Le contraste de Michelson

*Notes de cours :* [section 3.6](../../lecture_notes/fr/03_traces_and_states.md#s-3-6) ; d'après la [proposition 3.8](../../lecture_notes/fr/03_traces_and_states.md#r-3-8), une application positive et unifère n'augmente pas le contraste de Michelson, que la bibliothèque calcule par `Operator.michelson_contrast`.

```python
A_pos = alg.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(4, d), batch_size=4, force_positive=True)
Phi = channels.random_mixed_unitary_channel(alg, 3, batch_size=4)       # unital
assert (Phi(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
assert (E(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
```

<a id="s-3-7"></a>

## 3.7 Une inégalité plausible mais fausse

*Notes de cours :* [section 3.7](../../lecture_notes/fr/03_traces_and_states.md#s-3-7) ; d'après la [proposition 3.9](../../lecture_notes/fr/03_traces_and_states.md#r-3-9), pour $`X,Y\ge0`$ dans $`M_k`$, $`\max_U|\mathrm{Tr}(XUY)|=\|YX\|_1\ge\mathrm{Tr}(XY)`$, avec égalité si et seulement si $`XY=YX`$.

La bibliothèque trouve le point où le maximum est atteint grâce à la décomposition polaire : pour $`YX=W|YX|`$, l'unitaire $`U=W^*`$
atteint $`\|YX\|_1`$ ([notes de cours, théorème 1.9](../../lecture_notes/fr/01_matrix_algebras.md#r-1-9)). Une recherche aléatoire parmi des
unitaires distribués selon la loi de Haar (c'est ainsi que la bibliothèque tire $`U`$) fait moins bien en moyenne
([notes de cours, section 3.7](../../lecture_notes/fr/03_traces_and_states.md#s-3-7)) ; la borne de probabilité $`b^4/(a^4k^2)`$ qui y est rappelée ne dit rien
dans l'exécution ci-dessous, où $`b/a=6`$ et $`k=6`$. Le code calcule le supremum exact
pour un couple qui ne commute pas avec $`k=6`$ ainsi que le meilleur de $`20000`$ tirages de Haar,
tous deux diminués de $`\mathrm{Tr}(XY)`$. D'après la proposition 3.9, aucun tirage ne peut
dépasser le supremum ; l'écart entre le meilleur tirage et le supremum est une observation
numérique de cette exécution, pas un théorème.

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

## 3.8 Points de comparaison

QuTiP (`entropy_vn`, `entropy_relative`) et d'autres bibliothèques d'information quantique
(par exemple `toqito`) calculent les entropies et entropies relatives de matrices densité sur
un seul espace de Hilbert, avec la trace de cet espace ; sur $`\bigoplus_c\mathbb C^{k_c}`$, c'est
la trace brute d'ici (à la base du logarithme près, qui varie d'une bibliothèque à l'autre).
Ce qu'elles ne fournissent pas, c'est
la famille de traces du théorème 3.1 et les conversions du théorème 3.2, qui deviennent
nécessaires dès qu'une algèbre agit avec des multiplicités ou que les secteurs doivent être
pondérés différemment.

<a id="x-3"></a>

## Exercices

[Exercice 1](../../lecture_notes/fr/03_traces_and_states.md#x-3) des notes de cours : la densité de l'état tracial $`\tau_{\mathrm{vN}}`$
(`states.tracial_state`) relativement à $`\tau_{\mathrm{vN}}`$ est $`1`$, et son entropie vaut
$`S=\log C+\frac1C\sum_c\log k_c`$. Vérifier les deux pour $`M_2\oplus M_3\oplus M_5`$.

<details>
<summary>Solutions</summary>

```python
tau = states.tracial_state(alg)
assert (tau.density('tau_vN') - one.matrix).abs().max() < 1e-12
assert abs(tau.von_neumann_entropy().item() - (math.log(3) + sum(math.log(k) for k in (2, 3, 5)) / 3)) < 1e-12
```

</details>

**Suite :** [Tutoriel 4. Symétries et représentations](04_symmetries.md)
