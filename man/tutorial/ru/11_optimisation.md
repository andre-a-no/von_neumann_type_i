[English](../11_optimisation.md) · [Español](../es/11_optimisation.md) · [Français](../fr/11_optimisation.md) · [Deutsch](../de/11_optimisation.md) · [中文](../zh/11_optimisation.md) · [日本語](../ja/11_optimisation.md) · **Русский**

# Руководство 11. Оптимизация и следовые неравенства

Это руководство сопровождает [главу 11](../../lecture_notes/ru/11_optimisation.md) конспекта лекций, в которой доказываются экстремумы трёх следовых
функционалов и объясняется, почему выборка унитарных операторов их не находит. Здесь параметризации множеств ограничений
(`optimize.PositiveParam`, `UnitaryParam`, `SelfAdjointParam`) и мультистартовый оптимизатор (`optimize.extremize`)
библиотеки сверяются с этими теоремами. Пусть $`\Delta`$ — контраст Майкельсона и $`\ell(\delta)=(1-\delta)/(1+\delta)`$;
три функционала таковы:

| | функционал | переменные |
|---|---|---|
| Эксп. 1 | $`z=\lvert\mathrm{Tr}(XUY)\rvert-\mathrm{Tr}(XY)`$ | $`X,Y\ge0`$, $`U\in\mathcal M`$ унитарен |
| Эксп. 2 | $`z=\mathrm{Tr}\lvert XY\rvert-\mathrm{Tr}(\lvert X\rvert\,\lvert Y\rvert)`$ | $`X=UX_0`$, $`Y=VY_0`$, $`X_0,Y_0\ge0`$, $`U,V`$ унитарны |
| Эксп. 3 | $`z=\mathrm{Tr}(Y\lvert X\rvert Y)-\mathrm{Tr}\lvert YXY\rvert`$ | $`X=X^*`$, $`Y\ge0`$ |

где $`|X|`$, $`|Y|`$ нормированы условием $`\lambda_{\max}=1`$, а их контрасты служат параметрами. **Огибающие** — это
$`\overline z(\delta)=\sup\{z:\Delta(|X|)=\delta\}`$ и $`\underline z(\delta)=\inf\{z:\Delta(|X|)=\delta\}`$, где все остальные
переменные свободны.

<a id="s-11-1"></a>

## 11.1 Множества ограничений и их параметризации

*Конспект:* [раздел 11.1](../../lecture_notes/ru/11_optimisation.md#s-11-1); [лемма 11.1](../../lecture_notes/ru/11_optimisation.md#r-11-1), [предложение 11.2](../../lecture_notes/ru/11_optimisation.md#r-11-2).

`optimize.PositiveParam` с заданным контрастом $`\delta`$ реализует параметризацию из [раздела 11.1](../../lecture_notes/ru/11_optimisation.md#s-11-1)
конспекта лекций: для каждого старта он выбирает два различных слота $`s_\pm`$ (`spectrum.pin_max`, `spectrum.pin_min`; каждая
упорядоченная пара имеет вероятность $`1/(D(D-1))`$) и случайный по Хаару унитарный $`U_0`$ (`unitary.base`, выбирается функцией
`random_unitary_operator`, [лемма 11.10](../../lecture_notes/ru/11_optimisation.md#r-11-10)) и отображает параметры $`\theta`$ (`spectrum.theta`) и $`G`$ (`unitary.G`) в
```math
X(\theta,G)=U\mathrm{diag}(\lambda)U^*,\qquad U=U_0\,e^{G-G^*},\qquad
\lambda_{s_+}=1,\quad\lambda_{s_-}=\ell(\delta),\quad\lambda_s=\ell+(1-\ell)\tfrac{1+\sin\theta_s}2\ \text{otherwise}.
```
`UnitaryParam` — это одно лишь отображение $`G\mapsto U_0e^{G-G^*}`$, а `SelfAdjointParam` добавляет знаки:
$`X=U\mathrm{diag}(\sigma\lambda)U^*`$. По [предложению 11.2](../../lecture_notes/ru/11_optimisation.md#r-11-2) образ — это множество операторов оболочки контраста $`\mathcal P_\delta`$
со значениями $`1`$ и $`\ell(\delta)`$ в выбранных слотах. Свободный контраст в $`[a,b]`$ — это ещё один параметр; по умолчанию
в библиотеке $`[a,b]=[0,0.999]`$, и код ниже это значение сохраняет.

Код проверяет ограничение, а также конструктивно проверяет предложение 11.2: параметры из доказательства воспроизводят заданный оператор.

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

## 11.2 Неравенство фон Неймана для следов

*Конспект:* [раздел 11.2](../../lecture_notes/ru/11_optimisation.md#s-11-2); [теорема 11.4](../../lecture_notes/ru/11_optimisation.md#r-11-4), [следствие 11.5](../../lecture_notes/ru/11_optimisation.md#r-11-5), [предложение 11.6](../../lecture_notes/ru/11_optimisation.md#r-11-6), [следствие 11.7](../../lecture_notes/ru/11_optimisation.md#r-11-7).

Согласно [теореме 11.4](../../lecture_notes/ru/11_optimisation.md#r-11-4) конспекта лекций (неравенство фон Неймана для следов), $`|\mathrm{Tr}(AB)|\le\sum_is_i(A)s_i(B)`$, а по
[следствию 11.5](../../lecture_notes/ru/11_optimisation.md#r-11-5) $`\sum_ip^\downarrow_iq^\uparrow_i\le\mathrm{Tr}(PQ)\le\sum_ip^\downarrow_iq^\downarrow_i`$ для самосопряжённых $`P,Q`$ с собственными значениями $`p,q`$ и $`\mathrm{Tr}(PQ)\le\|PQ\|_1\le\sum_ip^\downarrow_iq^\downarrow_i`$ для
$`P,Q\ge0`$. Для эксп. 2 с фиксированными спектрами [предложение 11.6](../../lecture_notes/ru/11_optimisation.md#r-11-6) даёт $`\sup z=-\inf z=\sum_c\sum_ix^\downarrow_{c,i}(y^\downarrow_{c,i}-y^\uparrow_{c,i})`$, причём значение достигается
при $`X_0=\mathrm{diag}(x^\downarrow)`$, $`V=J`$ (обращающая перестановка) и $`Y_0=\mathrm{diag}(y^\uparrow)`$ (соответственно $`\mathrm{diag}(y^\downarrow)`$), а [следствие 11.7](../../lecture_notes/ru/11_optimisation.md#r-11-7) даёт
огибающую $`\sup z=-\inf z=\sum_c\lfloor k_c/2\rfloor(1-\ell_X)(1-\ell_Y)`$.

Здесь $`1-\ell_X=2\delta/(1+\delta)`$ при $`\Delta(|X|)=\delta`$; при нечётном $`k_c`$ среднее собственное значение роли не играет. Код
проверяет теорему 11.4, следствие 11.5 и экстремальные конфигурации предложения 11.6 для $`k=3`$, $`C=2`$.

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

## 11.3 Следовое неравенство с условием равенства

*Конспект:* [раздел 11.3](../../lecture_notes/ru/11_optimisation.md#s-11-3); [лемма 11.8](../../lecture_notes/ru/11_optimisation.md#r-11-8), [предложение 11.9](../../lecture_notes/ru/11_optimisation.md#r-11-9).

Согласно [предложению 11.9](../../lecture_notes/ru/11_optimisation.md#r-11-9) конспекта лекций, в эксп. 3 $`z=\mathrm{Tr}(Y|X|Y)-\mathrm{Tr}|YXY|\ge0`$, причём равенство имеет место тогда и только
тогда, когда $`YX_+Y^2X_-Y=0`$, в частности если $`XY=YX`$. Для $`X=\mathrm{diag}(1,-1,0)`$ и $`Y`$ со строками $`(1,0,a)`$, $`(0,1,0)`$, $`(a,0,1)`$, $`-1<a<1`$, равенство выполняется,
хотя $`XY\ne YX`$ при $`a\ne0`$. Код проверяет и то и другое (при $`a=0.5`$); `X.abs()` — это $`|X|`$, а `trace_norm()` — $`\|\cdot\|_1`$.

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

## 11.4 Выборка против оптимизации

*Конспект:* [раздел 11.4](../../lecture_notes/ru/11_optimisation.md#s-11-4); [лемма 11.10](../../lecture_notes/ru/11_optimisation.md#r-11-10), [предложение 11.11](../../lecture_notes/ru/11_optimisation.md#r-11-11), [следствие 11.12](../../lecture_notes/ru/11_optimisation.md#r-11-12), [лемма 11.13](../../lecture_notes/ru/11_optimisation.md#r-11-13).

Насколько велико обычно $`|\mathrm{Tr}(AU)|`$ по сравнению с его максимумом $`\|A\|_1`$? Библиотека выбирает $`U`$ как
множитель $`Q`$ в QR-разложении $`Z=QR`$, $`R_{ii}>0`$, **матрицы Жинибра** $`Z`$ (независимые элементы с
плотностью $`\pi^{-1}e^{-|z|^2}`$; рецепт Меццадри), независимо в каждом блоке.

Согласно [предложению 11.11](../../lecture_notes/ru/11_optimisation.md#r-11-11) конспекта лекций, $`\mathbb E|\mathrm{Tr}(AU)|^2=\|A\|_2^2/k`$ для распределённой по Хаару $`U\in U(k)`$, а по
[следствию 11.12](../../lecture_notes/ru/11_optimisation.md#r-11-12) для $`X,Y\ge0`$ со спектрами в $`[a,b]`$, $`a>0`$, среднеквадратичное значение $`|\mathrm{Tr}(XUY)|`$ не превосходит
его максимума по $`U`$, умноженного на $`\min(1/\sqrt k,b^2/(a^2k))`$, и $`\mathbb P\big(|\mathrm{Tr}(XUY)|\ge\mathrm{Tr}(XY)\big)\le b^4/(a^4k^2)`$.

Итак, $`N`$ выборок по Хаару находят нарушение неравенства $`|\mathrm{Tr}(XUY)|\le\mathrm{Tr}(XY)`$ с вероятностью не более
$`Nb^4/(a^4k^2)`$ (объединение $`N`$ событий), хотя нарушение существует для каждой некоммутирующей пары (предложение 3.9). Эта
оценка информативна лишь при $`k>\sqrt N\,b^2/a^2`$ (для $`b/a=6`$ и $`N=2000`$: $`k>1609`$); численно выборки пропускают
нарушение уже при гораздо меньших $`k`$ (см. ниже). Код проверяет предложение 11.11 методом Монте-Карло (выборочное среднее
отличается от точного значения меньше чем на пять стандартных ошибок — для этого запуска с фиксированным зерном) и для спектров, равномерно распределённых в $`[0.2,1.2]`$ ($`b/a=6`$), печатает
отношение из следствия 11.12, $`\sup z`$ и лучшее $`z`$ из $`2000`$ выборок, оба относительно $`\mathrm{Tr}(XY)`$.

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

*Численное наблюдение* (напечатанная таблица; утверждается только для этого запуска с фиксированным зерном): отношение остаётся
ниже доказанной оценки $`\min(1/\sqrt k,b^2/(a^2k))`$ и примерно уменьшается вдвое при каждом удвоении $`k`$, т. е. убывает как $`1/k`$,
а начиная с $`k=4`$ лучшая из $`2000`$ выборок не достигает даже $`z=0`$. Градиентный поиск, напротив, подходит близко к обоим
экстремумам в запуске ниже (это тоже численное наблюдение). Для пары `Xf, Yf` оба экстремума по $`U`$ известны:

Согласно [лемме 11.13](../../lecture_notes/ru/11_optimisation.md#r-11-13) конспекта лекций, $`\min_U|\mathrm{Tr}(AU)|=0`$ для всякой $`A\in M_k(\mathbb C)`$, $`k\ge2`$.

Поэтому $`\inf_Uz=-\mathrm{Tr}(XY)`$ и $`\sup_Uz=\|YX\|_1-\mathrm{Tr}(XY)`$ (предложение 3.9). Один вызов `extremize`
приближает оба, причём $`12`$ стартов максимизируют $`z`$, а $`12`$ — максимизируют $`-z`$ (с относительной точностью $`10^{-3}`$ в этом
запуске, что утверждается ниже). По предложению 11.11 среднеквадратичное значение
$`|\mathrm{Tr}(XUY)|`$ равно $`\|YX\|_2/\sqrt k`$; с оценками из доказательства следствия 11.12 (спектры в $`[a,b]`$)
$`\mathbb E|\mathrm{Tr}(XUY)|^2\le b^4`$ и $`\mathrm{Tr}(XY)\ge ka^2`$, так что неравенство Маркова даёт
$`\mathbb P\big(z\ge-(1-\varepsilon)\mathrm{Tr}(XY)\big)\le b^4/(\varepsilon^2a^4k^2)`$ для $`\varepsilon>0`$: при больших $`k`$
величина $`z`$ с высокой вероятностью близка к своему инфимуму.
*Численное наблюдение* (утверждается для этого запуска с фиксированным зерном): выборка подходит к инфимуму ближе чем на $`1\,\%`$,
но не достигает даже $`z=0`$, тогда как супремум положителен:

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

## 11.5 Мультистартовый поиск огибающих

*Конспект:* [раздел 11.5](../../lecture_notes/ru/11_optimisation.md#s-11-5); [предложение 11.14](../../lecture_notes/ru/11_optimisation.md#r-11-14), [следствие 11.7](../../lecture_notes/ru/11_optimisation.md#r-11-7), [предложение 11.9](../../lecture_notes/ru/11_optimisation.md#r-11-9), [предложение 11.15](../../lecture_notes/ru/11_optimisation.md#r-11-15), [следствие 11.16](../../lecture_notes/ru/11_optimisation.md#r-11-16).

`optimize.extremize(objective, params, maximize)` запускает Adam на `rounds` раундов по `steps` шагов (скорость обучения
`lr * decay**round`) с функцией потерь $`\mp\sum_b\text{objective}_b`$, после каждого раунда перецентрирует унитарные карты и
возвращает **для каждого старта** наибольшее (при `maximize=False` — наименьшее) из значений целевой функции, вычисленных в
градиентных итерациях и один раз в конечной точке. Каждый элемент пакета — независимый старт со своими слотами, знаками, $`U_0`$ и
контрастом: Adam действует покоординатно, поэтому если $`b`$-я компонента целевой функции зависит только от $`b`$-х параметров (как во
всех целевых функциях ниже), старты не взаимодействуют. Знаки
`SelfAdjointParam` дискретны: между раундами `discrete_step` для каждого старта пробует восемь одиночных перебросов знака (каждый
в равномерно случайном слоте), один за другим, и сохраняет каждый переброс, улучшающий целевую функцию.

В точной арифметике каждая итерация допустима (предложение 11.2 и замечания после него), поэтому каждое возвращённое значение
есть значение $`z`$ в допустимой точке; в арифметике с плавающей точкой — это значение с точностью до округления (утверждения ниже
допускают $`10^{-12}`$). Следовательно, максимизация даёт
**нижнюю оценку** супремума, минимизация — верхнюю оценку инфимума, и не более того: локальный поиск
может остановиться в локальном экстремуме. Оценка в другую сторону требует доказательства, как следствие 11.7 и предложение 11.9, или
выпуклой релаксации (раздел 11.6). Вычисленная огибающая проверяется (i) по замкнутым формулам; (ii) по
построенным вручную конфигурациям, например прямым суммам решений меньшего размера, ниже которых она опускаться не должна; (iii) по
структурированным выборкам (двухуровневые спектры, базисы-перестановки), значения которых — тоже достигаемые значения.

Согласно [предложению 11.14](../../lecture_notes/ru/11_optimisation.md#r-11-14) конспекта лекций, если все $`k_c\ge2`$, то нижняя огибающая эксп. 1 равна $`\underline z(\delta)=-\big(D-2\delta/(1+\delta)\big)`$.

Код проверяет предложение 11.14 на $`M_2\oplus M_2`$, следствие 11.7 на $`M_3\oplus M_3`$ (оба контраста заданы,
$`k=3`$ нечётно) и предложение 11.9; супремумы и инфимумы находятся одним вызовом, как выше. Односторонние утверждения (ни одно
возвращённое значение не выходит за доказанный экстремум) являются следствиями этих теорем; утверждения о том, что поиск подходит
к экстремуму с указанной точностью, — численные наблюдения этого запуска с фиксированным зерном. Ключевое слово `is_self_adjoint` выбирает
алгоритм для следовой нормы (сингулярные числа или собственные значения для $`YXY`$) без проверки на каждом шаге. Последняя
проверка относится к типу (ii): $`z`$ из эксп. 3 аддитивна относительно прямых сумм, а контраст прямой суммы равен наибольшему
из контрастов слагаемых (упражнение 3), поэтому для огибающих эксп. 3 с $`\Delta(Y)`$, свободным в $`[0,0.999]`$, как в коде,
```math
\overline z_{M_2\oplus M_2}(\delta)\ \ge\ \overline z_{M_2}(\delta)+\sup_{\delta'\le\delta}\overline z_{M_2}(\delta');
```
в частности, сумма двух значений, достигнутых на $`M_2`$ при контрастах $`\delta`$ и $`\delta'\le\delta`$, является значением,
достигаемым на $`M_2\oplus M_2`$ при контрасте $`\delta`$.

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

*Численные наблюдения* (напечатанные значения этого запуска). Поиск на $`M_2\oplus M_2`$ отстаёт от прямой суммы примерно на
$`0.2\,\%`$: он сошёлся не полностью, и
наилучшая нижняя оценка для $`\overline z_{M_2\oplus M_2}(0.5)`$ — большее из двух чисел; обнаруживать такое — назначение проверки
(ii). Значения для $`M_2`$ совпадают с точностью до $`10^{-4}`$ с $`(1-\delta)\big(2d/(1+d)\big)^2=(1-\delta)(1-\ell(d))^2`$, где
$`d=0.999`$ — наибольший допустимый контраст $`Y`$. Эта замкнутая формула, впервые наблюдавшаяся численно в
сопроводительной статье, является верхней огибающей из [следствия 11.16](../../lecture_notes/ru/11_optimisation.md#r-11-16)
(ниже); код утверждает только доказанное направление (ни одного значения выше неё). По тому же следствию прямая сумма
оптимальна: $`\overline z_{M_2\oplus M_2}(0.5)=1.5\,(2d/(1+d))^2`$. Для блоков большего размера замкнутая форма верхних
огибающих эксп. 1 и эксп. 3 нам неизвестна; вычисленные кривые тогда являются нижними оценками (достигнутыми значениями),
и не более того.

Согласно [предложению 11.15](../../lecture_notes/ru/11_optimisation.md#r-11-15) конспекта лекций, на $`M_2`$ для спектров
$`x_1\ge x_2`$ матрицы $`X`$ (матрицы $`|X|`$ в эксп. 3) и $`y_1\ge y_2`$ матрицы $`Y`$ супремум $`\sup z`$ равен
$`(x_1-x_2)^2(y_1-y_2)^2/(4(x_1+x_2)(y_1+y_2))`$ в эксп. 1 и $`\frac{2x_1x_2}{x_1+x_2}(y_1-y_2)^2`$ в эксп. 3. Если
$`e_1`$ — собственный вектор $`X`$ для собственного значения с модулем $`x_1`$, а $`f_1`$ — собственный вектор $`Y`$ для $`y_1`$,
то супремум достигается, когда $`\gamma=|\langle e_1,f_1\rangle|^2`$ равно $`\gamma_1=\frac14+\frac S{2q}`$ (эксп. 1, при
$`|\mathrm{Tr}(XUY)|=\|YX\|_1`$) или $`\gamma_3=\frac Sq`$ (эксп. 3, $`X=\mathrm{diag}(x_1,-x_2)`$), где $`S=x_1y_2+x_2y_1`$
и $`q=(x_1+x_2)(y_1+y_2)`$. По
[следствию 11.16](../../lecture_notes/ru/11_optimisation.md#r-11-16) для $`C`$ блоков размера $`2`$ при $`\Delta(|X|)=\delta`$ и
$`\Delta(Y)=\delta'`$ в эксп. 1 $`\sup z=C\delta^2\delta'^2/((1+\delta)(1+\delta'))`$ и
$`\inf z=-(2C-1+\ell(\delta)\ell(\delta'))`$, а в эксп. 3 $`\sup z=(C-\delta)\big(2\delta'/(1+\delta')\big)^2`$. Оба
супремума возрастают по $`\delta'`$: при $`\Delta(Y)`$, свободном в $`[0,d]`$, верхние огибающие равны значениям при $`\delta'=d`$,
а при $`d\to1`$ стремятся к $`C\delta^2/(2(1+\delta))`$ и $`C-\delta`$ — значениям для вырожденного $`Y`$
($`\Delta(Y)=1`$).

Код строит вещественные точки максимума из доказательства для $`C=1,2`$ (`from_blocks`), проверяет контрасты и значения,
проверяет оценки предложения 11.15 на случайных комплексных парах, достигает пределов $`d\to1`$ с проекторами ранга один
в каждом блоке и запускает мультистартовый поиск для эксп. 1 на $`M_2`$ вдоль $`\Delta(X)`$ при $`\Delta(Y)`$, свободном в
$`[0,0.999]`$; результат не должен превышать замкнутую формулу (теорему) и в этом запуске отличается от неё менее чем на
$`10^{-3}`$.

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

## 11.6 С чем сравнивать

**Manopt**, **pymanopt** и **geoopt** предоставляют матричные многообразия (многообразия Штифеля, ортогональные и унитарные группы)
с ретракциями — QR, полярной или Кэли — и римановы решатели (сопряжённые градиенты, доверительные области, риманов
Adam) — альтернативу Adam в перецентрируемой карте (сравнения производительности мы здесь не проводим). **cvxpy** передаёт
выпуклые задачи с операторными переменными SDP-решателям, которые возвращают оптимум вместе с двойственным сертификатом, оба —
с точностью до допуска решателя; релаксация — выпуклая задача на большем допустимом множестве — оценивает супремум сверху, т. е. в
том направлении, которое локальный поиск дать не может. **NCTSSOS**, **ncpol2sdpa** и **NCSOStools** оценивают собственные значения и следы некоммутативных
многочленов сразу во всех размерностях (раздел 0.3). Библиотека добавляет к этому точные параметризации спектральных ограничений
и ограничений на контраст для $`\bigoplus_cM_{k_c}`$ (экстремумы могут находиться в любом секторе), унитарные операторы, согласованные с
блоками, и одно пакетное вычисление для многих стартов и значений ограничений на CPU или GPU. Она даёт кандидатов в
экстремальные точки и нижние оценки супремумов; остальное должны дать доказательства или релаксации.

<a id="x-11"></a>

## Упражнения

[Упражнения 1 и 2](../../lecture_notes/ru/11_optimisation.md#x-11) конспекта лекций, численно: преобразование Кэли $`C(S)=(1-S)(1+S)^{-1}`$ косоэрмитовой $`S`$ унитарно,
не имеет собственного значения $`-1`$, и $`S=(1-C)(1+C)^{-1}`$ (упражнение 1); для $`O`$, распределённой по Хаару на $`O(k)`$ и выбранной
функцией `random_unitary_operator` вещественной алгебры (`complex_valued=False`), $`\mathbb E\mathrm{Tr}(AO)^2=\|A\|_2^2/k`$ (упражнение 2).

<details>
<summary>Решения</summary>

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

**Далее:** [Руководство 12](12_equivariant_learning.md): эквивариантные слои как элементы коммутанта в библиотеке.
