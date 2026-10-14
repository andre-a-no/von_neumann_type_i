[English](../03_traces_and_states.md) · [Español](../es/03_traces_and_states.md) · [Français](../fr/03_traces_and_states.md) · **Deutsch** · [中文](../zh/03_traces_and_states.md) · [日本語](../ja/03_traces_and_states.md) · [Русский](../ru/03_traces_and_states.md)

# Tutorial 3. Spuren, Zustände und Entropie

Dieses Tutorial begleitet [Kapitel 3](../../lecture_notes/de/03_traces_and_states.md) der Vorlesungsnotizen, das die Spuren
auf $\mathcal M=\bigoplus_cM_{k_c}$ klassifiziert ([Satz 3.1](../../lecture_notes/de/03_traces_and_states.md#r-3-1)), Zustände durch Dichten darstellt
([Satz 3.2](../../lecture_notes/de/03_traces_and_states.md#r-3-2)) und die grundlegenden Entropieungleichungen beweist. Wir verwenden die drei Spuren der
Bibliothek, rechnen Dichten zwischen ihnen um und überprüfen die Sätze numerisch auf $M_2\oplus M_3\oplus M_5$.

<a id="s-3-1"></a>

## 3.1 Spuren

*Vorlesungsnotizen:* [Abschnitt 3.1](../../lecture_notes/de/03_traces_and_states.md#s-3-1); nach [Satz 3.1](../../lecture_notes/de/03_traces_and_states.md#r-3-1) ist jede Spur von der Form $\mathrm{tr}_w(A)=\sum_cw_c\mathrm{Tr}\,A_c$.

Die Bibliothek implementiert drei treue Spuren:

| Name | Gewichte $w_c$ | Wert auf $1$ | Bedeutung |
|---|---|---|---|
| `Tr_blunt` | $1$ | $D$ | Spur von $\bigoplus_c\mathbb C^{k_c}$ (Vielfachheiten $m_c=1$) |
| `Tr_norm` | $1/k_c$ | $C$ | jeder Block normiert |
| `tau_vN` | $1/(Ck_c)$ | $1$ | der Spurzustand mit gleichen Sektorgewichten $1/C$ |

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

## 3.2 Zustände und Dichten

*Vorlesungsnotizen:* [Abschnitt 3.2](../../lecture_notes/de/03_traces_and_states.md#s-3-2); [Satz 3.2](../../lecture_notes/de/03_traces_and_states.md#r-3-2) liefert die Dichte $\rho_w$ eines Zustands für jede treue Spur, mit $w_c(\rho_w)_c=w'_c(\rho_{w'})_c$.

Die Bibliothek speichert einen Zustand durch seine Dichte $\rho=\rho_{\mathrm{blunt}}$ bezüglich
`Tr_blunt` (Klasse `DensityMatrix`); `density(trace)` und `DensityMatrix.from_density`
rechnen um. Der Wert des Zustands hängt nicht von der Konvention ab:

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

## 3.3 Der klassische Teil eines Zustands

*Vorlesungsnotizen:* [Abschnitt 3.3](../../lecture_notes/de/03_traces_and_states.md#s-3-3): die Sektorwahrscheinlichkeiten $p_c=\mathrm{Tr}\,\rho_c$ und die Abbildung $\rho\mapsto\sum_cP_c\rho P_c$, die auf $\mathcal M$ die Identität ist.

```python
p = rho.sector_probabilities()                                   # shape (batch, C)
assert torch.allclose(p.sum(-1), torch.ones(4, dtype=p.dtype))
P = [alg.central([1.0 if c == d else 0.0 for d in range(alg.C)]) for c in range(alg.C)]
measure_sector = channels.measurement_channel(P)
assert (measure_sector(rho).matrix - rho.matrix).abs().max() < 1e-12   # no disturbance
```

<a id="s-3-4"></a>

## 3.4 Entropie

*Vorlesungsnotizen:* [Abschnitt 3.4](../../lecture_notes/de/03_traces_and_states.md#s-3-4); [Satz 3.3](../../lecture_notes/de/03_traces_and_states.md#r-3-3) liefert $S(\rho)=H(p)+\sum_cp_cS(\sigma_c)$, und die Entropie bezüglich $\mathrm{tr}_w$ ist $S_w(\rho)=S(\rho)+\sum_cp_c\log w_c$.

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

## 3.5 Relative Entropie und die Kleinsche Ungleichung

*Vorlesungsnotizen:* [Abschnitt 3.5](../../lecture_notes/de/03_traces_and_states.md#s-3-5); [Korollar 3.6](../../lecture_notes/de/03_traces_and_states.md#r-3-6) ($S(\rho)\le\log D$, mit Gleichheit für $\rho=1/D$) und [Satz 3.7](../../lecture_notes/de/03_traces_and_states.md#r-3-7) ($S(\rho\|\sigma)\ge S(E\rho\|E\sigma)=\sum_cp_c\log(p_c/q_c)$ für die bedingte Erwartung $E$ auf das Zentrum).

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

## 3.6 Der Michelson-Kontrast

*Vorlesungsnotizen:* [Abschnitt 3.6](../../lecture_notes/de/03_traces_and_states.md#s-3-6); nach [Proposition 3.8](../../lecture_notes/de/03_traces_and_states.md#r-3-8) vergrößert eine positive unitale Abbildung den Michelson-Kontrast nicht, den die Bibliothek als `Operator.michelson_contrast` berechnet.

```python
A_pos = alg.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(4, d), batch_size=4, force_positive=True)
Phi = channels.random_mixed_unitary_channel(alg, 3, batch_size=4)       # unital
assert (Phi(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
assert (E(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
```

<a id="s-3-7"></a>

## 3.7 Eine plausible Ungleichung, die falsch ist

*Vorlesungsnotizen:* [Abschnitt 3.7](../../lecture_notes/de/03_traces_and_states.md#s-3-7); nach [Proposition 3.9](../../lecture_notes/de/03_traces_and_states.md#r-3-9) gilt für $X,Y\ge0$ in $M_k$ $\max_U|\mathrm{Tr}(XUY)|=\|YX\|_1\ge\mathrm{Tr}(XY)$, mit Gleichheit genau dann, wenn $XY=YX$.

Die Bibliothek findet die Maximalstelle mithilfe der Polarzerlegung: Für $YX=W|YX|$ erreicht das unitäre $U=W^*$
den Wert $\|YX\|_1$ ([Vorlesungsnotizen, Satz 1.9](../../lecture_notes/de/01_matrix_algebras.md#r-1-9)). Zufällige Suche über
Haar-verteilte unitäre Matrizen (die Bibliothek zieht $U$ auf diese Weise) ist im Mittel schlechter
([Vorlesungsnotizen, Abschnitt 3.7](../../lecture_notes/de/03_traces_and_states.md#s-3-7)); die dort angegebene Wahrscheinlichkeitsschranke $b^4/(a^4k^2)$ sagt
im folgenden Lauf nichts, da $b/a=6$ und $k=6$. Der Code berechnet das exakte Supremum für ein nicht
kommutierendes Paar mit $k=6$ und das Beste aus $20000$ Haar-Stichproben, jeweils abzüglich
$\mathrm{Tr}(XY)$. Nach Proposition 3.9 kann keine Stichprobe das Supremum übertreffen; wie
weit die beste Stichprobe darunter bleibt, ist eine numerische Beobachtung dieses Laufs, kein
Satz.

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

## 3.8 Womit man vergleichen sollte

QuTiP (`entropy_vn`, `entropy_relative`) und andere Pakete für Quanteninformation (z. B.
`toqito`) berechnen Entropien und relative Entropien von Dichtematrizen auf einem
einzelnen Hilbertraum, mit der Spur dieses Raumes; auf $\bigoplus_c\mathbb C^{k_c}$ ist das die
hiesige ungewichtete Spur (bis auf die Basis des Logarithmus, die sich zwischen den Paketen
unterscheidet). Was sie nicht bereitstellen, ist die Familie
von Spuren aus Satz 3.1 und die Umrechnungen aus Satz 3.2, die notwendig werden, sobald
eine Algebra mit Vielfachheiten wirkt oder die Sektoren unterschiedlich gewichtet werden
müssen.

<a id="x-3"></a>

## Aufgaben

[Aufgabe 1](../../lecture_notes/de/03_traces_and_states.md#x-3) der Vorlesungsnotizen: Die Dichte des Spurzustands $\tau_{\mathrm{vN}}$
(`states.tracial_state`) bezüglich $\tau_{\mathrm{vN}}$ ist $1$, und seine Entropie ist
$S=\log C+\frac1C\sum_c\log k_c$. Man überprüfe beides für $M_2\oplus M_3\oplus M_5$.

<details>
<summary>Lösungen</summary>

```python
tau = states.tracial_state(alg)
assert (tau.density('tau_vN') - one.matrix).abs().max() < 1e-12
assert abs(tau.von_neumann_entropy().item() - (math.log(3) + sum(math.log(k) for k in (2, 3, 5)) / 3)) < 1e-12
```

</details>

**Weiter:** [Tutorial 4. Symmetrien und Darstellungen](04_symmetries.md)
