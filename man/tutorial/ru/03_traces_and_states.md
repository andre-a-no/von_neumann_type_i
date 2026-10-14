[English](../03_traces_and_states.md) · [Español](../es/03_traces_and_states.md) · [Français](../fr/03_traces_and_states.md) · [Deutsch](../de/03_traces_and_states.md) · [中文](../zh/03_traces_and_states.md) · [日本語](../ja/03_traces_and_states.md) · **Русский**

# Руководство 3. Следы, состояния и энтропия

Это руководство сопровождает [главу 3](../../lecture_notes/ru/03_traces_and_states.md) конспекта лекций, где классифицируются следы
на $\mathcal M=\bigoplus_cM_{k_c}$ ([теорема 3.1](../../lecture_notes/ru/03_traces_and_states.md#r-3-1)), состояния представляются плотностями
([теорема 3.2](../../lecture_notes/ru/03_traces_and_states.md#r-3-2)) и доказываются основные энтропийные неравенства. Мы используем три следа
библиотеки, переводим плотности из одного соглашения в другое и численно проверяем теоремы на $M_2\oplus M_3\oplus M_5$.

<a id="s-3-1"></a>

## 3.1 Следы

*Конспект:* [раздел 3.1](../../lecture_notes/ru/03_traces_and_states.md#s-3-1); по [теореме 3.1](../../lecture_notes/ru/03_traces_and_states.md#r-3-1) всякий след имеет вид $\mathrm{tr}_w(A)=\sum_cw_c\mathrm{Tr}\,A_c$.

В библиотеке реализованы три точных следа:

| название | веса $w_c$ | значение на $1$ | смысл |
|---|---|---|---|
| `Tr_blunt` | $1$ | $D$ | след $\bigoplus_c\mathbb C^{k_c}$ (кратности $m_c=1$) |
| `Tr_norm` | $1/k_c$ | $C$ | каждый блок нормирован |
| `tau_vN` | $1/(Ck_c)$ | $1$ | следовое состояние с равными весами секторов $1/C$ |

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

## 3.2 Состояния и плотности

*Конспект:* [раздел 3.2](../../lecture_notes/ru/03_traces_and_states.md#s-3-2); [теорема 3.2](../../lecture_notes/ru/03_traces_and_states.md#r-3-2) даёт плотность $\rho_w$ состояния для любого точного следа, причём $w_c(\rho_w)_c=w'_c(\rho_{w'})_c$.

Библиотека хранит состояние через его плотность $\rho=\rho_{\mathrm{blunt}}$ относительно
`Tr_blunt` (класс `DensityMatrix`); для перехода между соглашениями служат `density(trace)` и
`DensityMatrix.from_density`. Значение состояния от соглашения не зависит:

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

## 3.3 Классическая часть состояния

*Конспект:* [раздел 3.3](../../lecture_notes/ru/03_traces_and_states.md#s-3-3): вероятности секторов $p_c=\mathrm{Tr}\,\rho_c$ и отображение $\rho\mapsto\sum_cP_c\rho P_c$, тождественное на $\mathcal M$.

```python
p = rho.sector_probabilities()                                   # shape (batch, C)
assert torch.allclose(p.sum(-1), torch.ones(4, dtype=p.dtype))
P = [alg.central([1.0 if c == d else 0.0 for d in range(alg.C)]) for c in range(alg.C)]
measure_sector = channels.measurement_channel(P)
assert (measure_sector(rho).matrix - rho.matrix).abs().max() < 1e-12   # no disturbance
```

<a id="s-3-4"></a>

## 3.4 Энтропия

*Конспект:* [раздел 3.4](../../lecture_notes/ru/03_traces_and_states.md#s-3-4); [теорема 3.3](../../lecture_notes/ru/03_traces_and_states.md#r-3-3) даёт $S(\rho)=H(p)+\sum_cp_cS(\sigma_c)$, а энтропия относительно $\mathrm{tr}_w$ равна $S_w(\rho)=S(\rho)+\sum_cp_c\log w_c$.

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

## 3.5 Относительная энтропия и неравенство Клейна

*Конспект:* [раздел 3.5](../../lecture_notes/ru/03_traces_and_states.md#s-3-5); [следствие 3.6](../../lecture_notes/ru/03_traces_and_states.md#r-3-6) ($S(\rho)\le\log D$, с равенством при $\rho=1/D$) и [теорема 3.7](../../lecture_notes/ru/03_traces_and_states.md#r-3-7) ($S(\rho\|\sigma)\ge S(E\rho\|E\sigma)=\sum_cp_c\log(p_c/q_c)$ для условного ожидания $E$ на центр).

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

## 3.6 Контраст Майкельсона

*Конспект:* [раздел 3.6](../../lecture_notes/ru/03_traces_and_states.md#s-3-6); по [предложению 3.8](../../lecture_notes/ru/03_traces_and_states.md#r-3-8) положительное унитальное отображение не увеличивает контраст Майкельсона, который библиотека вычисляет как `Operator.michelson_contrast`.

```python
A_pos = alg.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(4, d), batch_size=4, force_positive=True)
Phi = channels.random_mixed_unitary_channel(alg, 3, batch_size=4)       # unital
assert (Phi(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
assert (E(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
```

<a id="s-3-7"></a>

## 3.7 Правдоподобное неравенство, которое неверно

*Конспект:* [раздел 3.7](../../lecture_notes/ru/03_traces_and_states.md#s-3-7); по [предложению 3.9](../../lecture_notes/ru/03_traces_and_states.md#r-3-9) для $X,Y\ge0$ в $M_k$ имеем $\max_U|\mathrm{Tr}(XUY)|=\|YX\|_1\ge\mathrm{Tr}(XY)$, причём равенство имеет место тогда и только тогда, когда $XY=YX$.

Библиотека находит точку максимума с помощью полярного разложения: для $YX=W|YX|$
унитарная матрица $U=W^*$ достигает $\|YX\|_1$ ([конспект лекций, теорема 1.9](../../lecture_notes/ru/01_matrix_algebras.md#r-1-9)). Случайный поиск по
унитарным матрицам, распределённым по Хаару (библиотека выбирает $U$ именно так), в среднем хуже
([конспект лекций, раздел 3.7](../../lecture_notes/ru/03_traces_and_states.md#s-3-7)); приведённая там оценка вероятности $b^4/(a^4k^2)$ ничего не даёт
в запуске ниже, где $b/a=6$ и $k=6$. Код
вычисляет точный супремум для одной некоммутирующей пары с $k=6$ и лучшую из $20000$ выборок по Хаару,
обе величины за вычетом $\mathrm{Tr}(XY)$. По предложению 3.9 ни одна выборка не может превзойти супремум; насколько
лучшая выборка остаётся ниже него — это *численное наблюдение* данного запуска, а не теорема.

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

## 3.8 С чем сравнивать

QuTiP (`entropy_vn`, `entropy_relative`) и другие пакеты квантовой теории информации (например,
`toqito`) вычисляют энтропии и относительные энтропии матриц плотности на
одном гильбертовом пространстве, со следом этого пространства; на $\bigoplus_c\mathbb C^{k_c}$ это
след `Tr_blunt` в наших обозначениях (с точностью до основания логарифма, которое в разных пакетах различно). Чего они не предоставляют, так это семейства следов из теоремы 3.1 и
переходов из теоремы 3.2, которые становятся необходимы, как только алгебра действует с
кратностями или секторы нужно взвешивать по-разному.

<a id="x-3"></a>

## Упражнения

[Упражнение 1](../../lecture_notes/ru/03_traces_and_states.md#x-3) конспекта лекций: плотность следового состояния $\tau_{\mathrm{vN}}$
(`states.tracial_state`) относительно $\tau_{\mathrm{vN}}$ равна $1$, а его энтропия равна
$S=\log C+\frac1C\sum_c\log k_c$. Проверьте оба утверждения для $M_2\oplus M_3\oplus M_5$.

<details>
<summary>Решения</summary>

```python
tau = states.tracial_state(alg)
assert (tau.density('tau_vN') - one.matrix).abs().max() < 1e-12
assert abs(tau.von_neumann_entropy().item() - (math.log(3) + sum(math.log(k) for k in (2, 3, 5)) / 3)) < 1e-12
```

</details>

**Далее:** [Руководство 4. Симметрии и представления](04_symmetries.md)
