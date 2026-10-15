[English](../03_traces_and_states.md) · **Español** · [Français](../fr/03_traces_and_states.md) · [Deutsch](../de/03_traces_and_states.md) · [中文](../zh/03_traces_and_states.md) · [日本語](../ja/03_traces_and_states.md) · [Русский](../ru/03_traces_and_states.md)

# Tutorial 3. Trazas, estados y entropía

Este tutorial acompaña al [capítulo 3](../../lecture_notes/es/03_traces_and_states.md) de los apuntes de clase, que clasifica las trazas
sobre $`\mathcal M=\bigoplus_cM_{k_c}`$ ([teorema 3.1](../../lecture_notes/es/03_traces_and_states.md#r-3-1)), representa los estados mediante densidades
([teorema 3.2](../../lecture_notes/es/03_traces_and_states.md#r-3-2)) y demuestra las desigualdades básicas de entropía. Usamos las tres trazas de la
biblioteca, convertimos densidades de una a otra y comprobamos numéricamente los teoremas sobre $`M_2\oplus M_3\oplus M_5`$.

<a id="s-3-1"></a>

## 3.1 Trazas

*Apuntes:* [sección 3.1](../../lecture_notes/es/03_traces_and_states.md#s-3-1); por el [teorema 3.1](../../lecture_notes/es/03_traces_and_states.md#r-3-1), toda traza es $`\mathrm{tr}_w(A)=\sum_cw_c\mathrm{Tr}\,A_c`$.

La biblioteca implementa tres trazas fieles:

| nombre | pesos $`w_c`$ | valor en $`1`$ | significado |
|---|---|---|---|
| `Tr_blunt` | $`1`$ | $`D`$ | traza de $`\bigoplus_c\mathbb C^{k_c}`$ (multiplicidades $`m_c=1`$) |
| `Tr_norm` | $`1/k_c`$ | $`C`$ | cada bloque normalizado |
| `tau_vN` | $`1/(Ck_c)`$ | $`1`$ | el estado tracial con pesos iguales $`1/C`$ en los sectores |

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

## 3.2 Estados y densidades

*Apuntes:* [sección 3.2](../../lecture_notes/es/03_traces_and_states.md#s-3-2); el [teorema 3.2](../../lecture_notes/es/03_traces_and_states.md#r-3-2) da la densidad $`\rho_w`$ de un estado para toda traza fiel, con $`w_c(\rho_w)_c=w'_c(\rho_{w'})_c`$.

La biblioteca almacena un estado mediante su densidad $`\rho=\rho_{\mathrm{blunt}}`$ respecto de
`Tr_blunt` (clase `DensityMatrix`); `density(trace)` y `DensityMatrix.from_density`
realizan la conversión. El valor del estado no depende del convenio:

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

## 3.3 La parte clásica de un estado

*Apuntes:* [sección 3.3](../../lecture_notes/es/03_traces_and_states.md#s-3-3): las probabilidades de los sectores $`p_c=\mathrm{Tr}\,\rho_c`$, y la aplicación $`\rho\mapsto\sum_cP_c\rho P_c`$, que es la identidad sobre $`\mathcal M`$.

```python
p = rho.sector_probabilities()                                   # shape (batch, C)
assert torch.allclose(p.sum(-1), torch.ones(4, dtype=p.dtype))
P = [alg.central([1.0 if c == d else 0.0 for d in range(alg.C)]) for c in range(alg.C)]
measure_sector = channels.measurement_channel(P)
assert (measure_sector(rho).matrix - rho.matrix).abs().max() < 1e-12   # no disturbance
```

<a id="s-3-4"></a>

## 3.4 Entropía

*Apuntes:* [sección 3.4](../../lecture_notes/es/03_traces_and_states.md#s-3-4); el [teorema 3.3](../../lecture_notes/es/03_traces_and_states.md#r-3-3) da $`S(\rho)=H(p)+\sum_cp_cS(\sigma_c)`$, y la entropía respecto de $`\mathrm{tr}_w`$ es $`S_w(\rho)=S(\rho)+\sum_cp_c\log w_c`$.

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

## 3.5 Entropía relativa y desigualdad de Klein

*Apuntes:* [sección 3.5](../../lecture_notes/es/03_traces_and_states.md#s-3-5); [corolario 3.6](../../lecture_notes/es/03_traces_and_states.md#r-3-6) ($`S(\rho)\le\log D`$, con igualdad para $`\rho=1/D`$) y [teorema 3.7](../../lecture_notes/es/03_traces_and_states.md#r-3-7) ($`S(\rho\|\sigma)\ge S(E\rho\|E\sigma)=\sum_cp_c\log(p_c/q_c)`$ para la esperanza condicional $`E`$ sobre el centro).

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

## 3.6 El contraste de Michelson

*Apuntes:* [sección 3.6](../../lecture_notes/es/03_traces_and_states.md#s-3-6); por la [proposición 3.8](../../lecture_notes/es/03_traces_and_states.md#r-3-8), una aplicación positiva y unital no aumenta el contraste de Michelson, que la biblioteca calcula como `Operator.michelson_contrast`.

```python
A_pos = alg.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(4, d), batch_size=4, force_positive=True)
Phi = channels.random_mixed_unitary_channel(alg, 3, batch_size=4)       # unital
assert (Phi(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
assert (E(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
```

<a id="s-3-7"></a>

## 3.7 Una desigualdad verosímil que es falsa

*Apuntes:* [sección 3.7](../../lecture_notes/es/03_traces_and_states.md#s-3-7); por la [proposición 3.9](../../lecture_notes/es/03_traces_and_states.md#r-3-9), para $`X,Y\ge0`$ en $`M_k`$, $`\max_U|\mathrm{Tr}(XUY)|=\|YX\|_1\ge\mathrm{Tr}(XY)`$, con igualdad si y solo si $`XY=YX`$.

La biblioteca encuentra el maximizador mediante la descomposición polar: para $`YX=W|YX|`$ la unitaria $`U=W^*`$
alcanza $`\|YX\|_1`$ ([apuntes de clase, teorema 1.9](../../lecture_notes/es/01_matrix_algebras.md#r-1-9)). La búsqueda aleatoria sobre
unitarias distribuidas según la medida de Haar (la biblioteca muestrea $`U`$ de este modo) es peor en promedio
([apuntes de clase, sección 3.7](../../lecture_notes/es/03_traces_and_states.md#s-3-7)); la cota de probabilidad $`b^4/(a^4k^2)`$ que allí se recuerda no dice nada
en la ejecución siguiente, donde $`b/a=6`$ y $`k=6`$. El código
calcula el supremo exacto para un par que no conmuta con $`k=6`$ y la mejor de $`20000`$ muestras de
Haar, ambos menos $`\mathrm{Tr}(XY)`$. Por la proposición 3.9, ninguna muestra puede superar el supremo; cuánto
queda la mejor muestra por debajo de él es una observación numérica de esta ejecución, no un teorema.

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

## 3.8 Con qué comparar

QuTiP (`entropy_vn`, `entropy_relative`) y otros paquetes de información cuántica (por
ejemplo, `toqito`) calculan entropías y entropías relativas de matrices densidad sobre un
único espacio de Hilbert, con la traza de ese espacio; sobre $`\bigoplus_c\mathbb C^{k_c}`$ esta es la
traza `blunt` de aquí (salvo por la base del logaritmo, que difiere de un paquete a otro). Lo que no proporcionan es la familia de
trazas del teorema 3.1 ni las conversiones del teorema 3.2, que se vuelven necesarias en
cuanto un álgebra actúa con multiplicidades o los sectores deben ponderarse de forma
distinta.

<a id="x-3"></a>

## Ejercicios

[Ejercicio 1](../../lecture_notes/es/03_traces_and_states.md#x-3) de los apuntes de clase: la densidad del estado tracial $`\tau_{\mathrm{vN}}`$
(`states.tracial_state`) respecto de $`\tau_{\mathrm{vN}}`$ es $`1`$, y su entropía es
$`S=\log C+\frac1C\sum_c\log k_c`$. Compruebe ambas cosas para $`M_2\oplus M_3\oplus M_5`$.

<details>
<summary>Soluciones</summary>

```python
tau = states.tracial_state(alg)
assert (tau.density('tau_vN') - one.matrix).abs().max() < 1e-12
assert abs(tau.von_neumann_entropy().item() - (math.log(3) + sum(math.log(k) for k in (2, 3, 5)) / 3)) < 1e-12
```

</details>

**Siguiente:** [Tutorial 4. Simetrías y representaciones](04_symmetries.md)
