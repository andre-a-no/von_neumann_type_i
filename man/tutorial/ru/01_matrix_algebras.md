[English](../01_matrix_algebras.md) · [Español](../es/01_matrix_algebras.md) · [Français](../fr/01_matrix_algebras.md) · [Deutsch](../de/01_matrix_algebras.md) · [中文](../zh/01_matrix_algebras.md) · [日本語](../ja/01_matrix_algebras.md) · **Русский**

# Руководство 1. Матричные алгебры и функциональное исчисление

Это руководство сопровождает [главу 1](../../lecture_notes/ru/01_matrix_algebras.md) конспекта лекций, в которой доказываются спектральная
теорема, функциональное исчисление, полярное разложение и двойственность для следовой нормы в $M_n(\mathbb C)$.
Здесь каждый из этих результатов проверяется с помощью библиотеки.

**Вы узнаете**, как строить операторы с заданным спектром, как применять функцию к
оператору, как хранятся пакеты операторов и как проверять вычисление по точному
результату.

<a id="s-1-5"></a>

## 1.5 Операторы в библиотеке

*Конспект:* [раздел 1.2](../../lecture_notes/ru/01_matrix_algebras.md#s-1-2) и [раздел 1.3](../../lecture_notes/ru/01_matrix_algebras.md#s-1-3): [теорема 1.3](../../lecture_notes/ru/01_matrix_algebras.md#r-1-3), [следствие 1.4](../../lecture_notes/ru/01_matrix_algebras.md#r-1-4), [теорема 1.5](../../lecture_notes/ru/01_matrix_algebras.md#r-1-5), [теорема 1.6](../../lecture_notes/ru/01_matrix_algebras.md#r-1-6).

Библиотека работает с алгебрами $\bigoplus_cM_{k_c}(\mathbb C)$ (глава 2). Один
блок — это случай $C=1$: `TypeIAlgebra([n], [n])` есть $M_n(\mathbb C)$. (Два списка —
это номинальные и фактические размеры блоков; пока они совпадают.)

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

Каждый оператор хранится как один тензор формы `(batch, C, k_max, k_max)`: пакет
операторов, каждый из которых — список из `C` блоков. Здесь `batch = C = 1`. Построение —
это в точности теорема 1.3, прочитанная в обратную сторону: $A=U\Lambda U^*$ с $U$,
распределённой по мере Хаара (библиотека порождает $U$ QR-конструкцией из теоремы 5.3). Флаг
`force_positive=True` проверяется по заданным собственным значениям (библиотека выдаёт ошибку, если
одно из них меньше $-10^{-12}$); по следствию 1.4(3) результат тогда положителен (точно — когда заданные
собственные значения $\ge0$, как здесь; в общем случае — с точностью до этого допуска), и библиотека
запоминает это, а не проверяет повторно. По той же причине `A.eigenvalues()` выше возвращает заданный
спектр, отсортированный, ничего не вычисляя; это не проверка. Проверка
$\mathrm{spec}\,f(A)=f(\mathrm{spec}\,A)$ ниже действительно вычисляет: `apply_function` диагонализует
матрицу $A$, а собственные значения полученной матрицы $e^A$ вычисляются заново.

Функциональное исчисление для самосопряжённых операторов реализует `apply_function` (метод
диагонализует каждый блок с помощью `eigh`). Функция $f$ должна действовать поэлементно: она
получает собственные значения всех секторов одного размера $k$ одним тензором формы
`(batch, n·k)`, где $n$ — число таких секторов, и возвращает тензор той же формы. Для
распространённых функций есть отдельные методы. Проверим теоремы 1.5 и 1.6:

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

`R @ R` пока ничего не вычисляет: произведения, суммы и функции строят «рецепт»,
который вычисляется, когда запрашивается `.matrix` или число, например след.

<a id="s-1-6"></a>

## 1.6 Нормы и двойственность для следовой нормы

*Конспект:* [раздел 1.4](../../lecture_notes/ru/01_matrix_algebras.md#s-1-4): [теорема 1.8](../../lecture_notes/ru/01_matrix_algebras.md#r-1-8), три нормы, [теорема 1.9](../../lecture_notes/ru/01_matrix_algebras.md#r-1-9).

Произвольная (не обязательно нормальная) матрица вводится через `alg.operator(tensor)`.
Сверим три нормы с сингулярными числами, вычисленными непосредственно в PyTorch, а также
проверим полярное разложение (теорема 1.8) и двойственность (теорема 1.9):

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

`torch.randn` с комплексным dtype порождает независимые элементы с независимыми вещественными и мнимыми
частями $\sim N(0,\tfrac12)$, так что `G` — матрица Жинибра, и она имеет плотность (предложение A.16(b)).
$|X|$ обратима почти наверное: $(\det|X|)^2=\det(|X|^2)=\det(X^*X)=|\det X|^2$, а $|\det X|^2$ — многочлен
от вещественных и мнимых частей элементов, равный $1$ при $X=1$, поэтому он обращается в нуль с
вероятностью нуль (следствие A.7(b); это также первый шаг теоремы 5.3).

<a id="s-1-7"></a>

## 1.7 Пакеты: много операторов сразу

*Конспект:* [раздел 1.7](../../lecture_notes/ru/01_matrix_algebras.md#s-1-7), [предложение 1.10](../../lecture_notes/ru/01_matrix_algebras.md#r-1-10).

Ось пакета — причина, по которой существует библиотека: теории случайных матриц,
методу Монте-Карло и оптимизации с многими начальными точками нужны тысячи операторов, и
один пакетный вызов заменяет цикл на Python по ним. В качестве примера проверим точный
момент **гауссовского унитарного ансамбля** (GUE).

Согласно [конспекту лекций, предложение 1.10](../../lecture_notes/ru/01_matrix_algebras.md#r-1-10), если элементы $G$ независимы, с независимыми
вещественными и мнимыми частями $\sim N(0,\tfrac12)$, и $X=(G+G^*)/\sqrt{2n}$, то $\mathbb E\,\tfrac1n\mathrm{Tr}\,X^4=2+\tfrac1{n^2}$.

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

Проверка методом Монте-Карло — это статистический тест: оценка должна совпадать с точным
значением с точностью до нескольких стандартных ошибок. Величина $\tfrac1n\mathrm{Tr}\,X^4$ — многочлен от
гауссовых величин, поэтому у неё конечная дисперсия $\sigma^2$ (лемма A.23: все гауссовы моменты конечны;
произведения оцениваются через $|y_1\cdots y_r|\le\sum_i|y_i|^r$),
и по неравенству Чебышёва среднее $B$ независимых выборок отклоняется от точного значения
больше чем на $5\sigma/\sqrt B$ с вероятностью не более $1/25$. Код использует выборочное стандартное
отклонение вместо $\sigma$, а фиксированное начальное значение генератора делает проверку воспроизводимой.
Так устроена в этом курсе каждая проверка по доказанным результатам: точное утверждение, численная оценка и явно заданный
допуск.

<a id="s-1-8"></a>

## 1.8 С чем сравнивать

Для отдельных матриц эталонными реализациями служат `numpy.linalg` и `scipy.linalg` (`sqrtm`, `expm`, `logm`,
`polar`), и приведённые выше проверки можно повторить с их помощью. Библиотека добавляет
три вещи, важные в дальнейшем: пакетную обработку (раздел 1.7), операторы, несущие
блочную структуру (глава 2), и автоматическое дифференцирование (autograd из PyTorch) через её операции.

<a id="x-1"></a>

## Упражнения

[Упражнения 1 и 3](../../lecture_notes/ru/01_matrix_algebras.md#x-1) конспекта лекций, численно:
- Проверьте равенство $(fg)(A)=f(A)g(A)$ для $f=\sin$, $g=\exp$ на случайной самосопряжённой $A$ (упражнение 1).
- **Геометрическое среднее** матриц $A,B>0$, $A\#B=A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}$, — единственное
  положительное решение $X$ уравнения $XA^{-1}X=B$ (упражнение 3). Проверьте численно это уравнение и
  равенство $A\#B=B\#A$.

<details>
<summary>Решения</summary>

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

**Далее:** [Руководство 2](02_von_neumann_algebras.md): коммутанты, центр и алгебра $\bigoplus_cM_{k_c}(\mathbb C)$ в библиотеке.
