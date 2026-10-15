[English](../04_symmetries.md) · [Español](../es/04_symmetries.md) · [Français](../fr/04_symmetries.md) · [Deutsch](../de/04_symmetries.md) · [中文](../zh/04_symmetries.md) · [日本語](../ja/04_symmetries.md) · **Русский**

# Руководство 4. Симметрии и представления

Это руководство сопровождает [главу 4](../../lecture_notes/ru/04_symmetries.md) конспекта лекций, в которой развивается теория представлений,
стоящая за алгебрами симметрий, превращается в рандомизированный алгоритм, раскладывающий алгебру фон Неймана, порождённую
произвольным списком матриц ([предложение 4.9](../../lecture_notes/ru/04_symmetries.md#r-4-9)), применяется к $`U(1)`$ и $`SU(2)`$, а секторы алгебры,
порождённой вещественными матрицами, разделяются на вещественный, комплексный и кватернионный тип
([теорема 4.20](../../lecture_notes/ru/04_symmetries.md#r-4-20)), а также анализируется разложение зашумлённых образующих
([теорема 4.25](../../lecture_notes/ru/04_symmetries.md#r-4-25)). Мы вычисляем характеры, запускаем разложение для групп и для образующих, строим
спиновые цепочки по секторам с помощью `SpinChain` и `su2.SU2Chain`, вычисляем вещественную структуру разложения и сертифицируем разложение зашумлённых образующих.

<a id="s-4-2"></a>

## 4.2 Изотипическое разложение и коммутант

*Конспект:* [раздел 4.2](../../lecture_notes/ru/04_symmetries.md#s-4-2); [теорема 4.3](../../lecture_notes/ru/04_symmetries.md#r-4-3), [теорема 4.4](../../lecture_notes/ru/04_symmetries.md#r-4-4), [следствие 4.5](../../lecture_notes/ru/04_symmetries.md#r-4-5).

В [разделе 4.2](../../lecture_notes/ru/04_symmetries.md#s-4-2) конспекта лекций вручную вычислены характеры трёх неприводимых представлений
$`S_3`$ — тривиального, знакового и стандартного на $`(1,1,1)^\perp\subseteq\mathbb C^3`$ — и представления $`\pi`$, переставляющего три
сомножителя $`(\mathbb C^2)^{\otimes3}`$, $`\chi_\pi(\sigma)=2^{\#\text{cycles}(\sigma)}`$ (таблица ниже), а из
[следствия 4.5](../../lecture_notes/ru/04_symmetries.md#r-4-5) получены кратности $`m_{\rm triv}=4`$, $`m_{\rm sgn}=0`$, $`m_{\rm std}=2`$ и $`\dim\pi(S_3)'=20`$. Код
пересчитывает таблицу: он находит $`P_\sigma`$ как сужение $`\pi(\sigma)`$ на инвариантную линейную оболочку
$`|100\rangle,|010\rangle,|001\rangle`$ (на которой $`\pi(\sigma)`$ переставляет позиции единственной $`1`$, лемма B.1(a)), а
$`\#\text{cycles}(\sigma)=\dim\ker(P_\sigma-1)`$ (векторы, постоянные на циклах).

| | тождественная ($`1`$) | транспозиции ($`3`$) | $`3`$-циклы ($`2`$) |
|---|---|---|---|
| $`\chi_{\rm triv}`$ | $`1`$ | $`1`$ | $`1`$ |
| $`\mathrm{sgn}`$ | $`1`$ | $`-1`$ | $`1`$ |
| $`\chi_{\rm std}`$ | $`2`$ | $`0`$ | $`-1`$ |
| $`\chi_\pi`$ | $`8`$ | $`4`$ | $`2`$ |

```python
import math
import torch
import torch_vn_algebra.decompose as dc

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    return (torch.as_tensor(X) - torch.as_tensor(Y)).abs().max().item() < tol

G = dc.permutation_representation(3, 2)                  # (6, 8, 8): pi(sigma) permutes three qubits
chi_pi = G.diagonal(dim1=-2, dim2=-1).sum(-1).real       # characters Tr pi(sigma)
P3 = G[:, [4, 2, 1]][:, :, [4, 2, 1]].real              # on |100>, |010>, |001>: S_3 permuting C^3
cycles = 3 - torch.linalg.matrix_rank(P3 - torch.eye(3, dtype=torch.float64)).double()   # dim ker(P - 1)
assert close(chi_pi, 2.0 ** cycles)
irreps = {'trivial': torch.ones(6, dtype=torch.float64), 'sign': (-1.0) ** (3 - cycles),
          'standard': P3.diagonal(dim1=-2, dim2=-1).sum(-1) - 1}                 # fixed points - 1
chis = torch.stack(list(irreps.values()))
assert close(chis @ chis.T / 6, torch.eye(3, dtype=torch.float64))    # Corollary 4.5(a): all irreducible
mult = chis @ chi_pi / 6                                               # Corollary 4.5(b)
print(dict(zip(irreps, mult.tolist())), (chi_pi ** 2).mean().item())
assert close(mult, torch.tensor([4.0, 0.0, 2.0], dtype=torch.float64))
assert dc.commutant_dimension(G) == 20                                # Corollary 4.5(c): M_4 (+) M_2
```

Это разложение из раздела 2.4, полученное теперь из таблицы по шесть чисел на представление.

<a id="s-4-3"></a>

## 4.3 От образующих к алгебре: рандомизированный алгоритм

*Конспект:* [раздел 4.3](../../lecture_notes/ru/04_symmetries.md#s-4-3); [лемма 4.6](../../lecture_notes/ru/04_symmetries.md#r-4-6), [предложение 4.7](../../lecture_notes/ru/04_symmetries.md#r-4-7), [лемма 4.8](../../lecture_notes/ru/04_symmetries.md#r-4-8), [предложение 4.9](../../lecture_notes/ru/04_symmetries.md#r-4-9).

Рандомизированный алгоритм из конспекта лекций — это `decompose.wedderburn` с интерфейсами `from_group` (все элементы
конечной группы; выборка из коммутанта делается твирлингом `decompose.twirl`) и `from_generators` (произвольный список
матриц; коммутант — ядро, `decompose.commutant_basis`). Он расщепляет $`\mathbb C^n`$ на собственные подпространства
случайного $`K\in\mathcal A'_h`$, группирует и выравнивает их с помощью блоков случайного $`Y\in\mathcal A'`$ и возвращает
`Decomposition` с унитарной матрицей `V`, числами $`k_c`$ (`k_factors`) и $`m_c`$ (`multiplicities`). В арифметике с плавающей
точкой он считает равными собственные значения, отличающиеся меньше чем на допуск, и дальше расщепляет полученные
подпространства, сжимая на них новые случайные элементы коммутанта ([лемма 4.8(a)](../../lecture_notes/ru/04_symmetries.md#r-4-8)).

`Decomposition.verify()` сообщает относительные расстояния образующих (и новых элементов коммутанта) до
$`\mathcal B`$ (соответственно $`\mathcal B'`$) и сравнивает $`\sum_cm_c^2`$ с $`\dim\mathcal A'`$ (из характеров для группы,
следствие 4.5(c); из ядра для образующих). Запустим его для $`S_3`$ и трёх кубитов, а затем заглянем внутрь:
собственные подпространства одного случайного $`K`$ и сплетающие операторы одного случайного $`Y`$.

```python
D = dc.from_group(G)                                     # the twirl samples pi(S_3)'
print(D)
assert D.k_factors == (1, 2) and D.multiplicities == (4, 2)            # Section 2.4: k = (1, 2), m = (4, 2)
assert D.verify()['ok']
rho = D.to_algebra(G).matrix                             # (6, 2, 2, 2): the blocks rho_c(sigma), padded
assert close(rho[:, 0, 0, 0], torch.ones(6))                           # sector 1: trivial representation
assert close(rho[:, 1].diagonal(dim1=-2, dim2=-1).sum(-1), chis[2])    # sector 2: the standard one
schur = torch.einsum('gij,gkl->ijkl', rho[:, 1], rho[:, 1].conj()) / 6
assert close(schur.reshape(4, 4), torch.eye(4) / 2)                    # Theorem 4.4 with k = 2
assert close(D.from_algebra(D.to_algebra(G)), G)

Z = torch.randn(8, 8, dtype=dt)
K = dc.twirl(G, Z + Z.mH)                                # a random element of pi(S_3)'_h
w, Q = torch.linalg.eigh(K)
cuts = [0] + [i + 1 for i in range(7) if w[i + 1] - w[i] > 1e-8] + [8]
spaces = [Q[:, a:b] for a, b in zip(cuts[:-1], cuts[1:])]
assert sorted(W.shape[1] for W in spaces) == [1, 1, 1, 1, 2, 2]      # Lemma 4.6: m_c pieces of size k_c
assert all(close(G @ W, W @ (W.mH @ G @ W)) for W in spaces)           # invariant subspaces
Y = dc.twirl(G, torch.randn(8, 8, dtype=dt))
Wa, Wb = [W for W in spaces if W.shape[1] == 2]
s = torch.linalg.svdvals(Wa.mH @ Y @ Wb)                               # Lemma 4.8(b): a multiple of a unitary
assert s[1] > 1e-3 and abs(s[0] - s[1]) < 1e-10
assert all(close(W.mH @ Y @ Wa, 0) for W in spaces if W.shape[1] == 1)   # different summands: zero
```

Основная стоимость — коммутант: одна выборка твирлинга стоит $`|G|`$ произведений матриц $`n\times n`$. Для образующих
библиотека сначала диагонализует эрмитов элемент $`R=U\mathrm{diag}(\lambda)U^*`$ алгебры $`\mathcal A`$ (случайную
вещественную комбинацию эрмитовых частей образующих плюс кратное симметризованного произведения двух таких
комбинаций). Каждый $`X\in\mathcal A'`$ коммутирует с $`R`$, а значит, и с его спектральными проекторами (теорема 1.5(d)), так что
$`U^*XU`$ блочно-диагонален относительно собственных подпространств $`R`$ размерностей $`d_\lambda`$; это накладывается
заранее, и задача о ядре имеет $`N_R=\sum_\lambda d_\lambda^2`$ неизвестных вместо $`n^2`$. Объединение нескольких
собственных значений лишь увеличивает блоки, так что у сокращённой задачи то же множество решений $`\mathcal A'`$.
В разложении $`V^*RV=\bigoplus_cr_c\otimes1_{m_c}`$ с эрмитовыми $`r_c\in M_{k_c}`$; если $`\lambda`$ имеет кратность
$`\mu_c(\lambda)`$ в $`r_c`$, то $`d_\lambda=\sum_c\mu_c(\lambda)m_c`$, и, поскольку $`\mu^2\ge\mu`$ для целых $`\mu\ge0`$,
```math
N_R=\sum_\lambda\Big(\sum_c\mu_c(\lambda)m_c\Big)^2\ge\sum_{\lambda,c}\mu_c(\lambda)^2m_c^2\ge\sum_{\lambda,c}\mu_c(\lambda)m_c^2=\sum_ck_cm_c^2,
```
с равенством тогда и только тогда, когда собственные значения всех $`r_c`$ простые и попарно различны. Для полного спина
$`L=3,4,5`$ спинов код ниже находит равенство (*численное наблюдение*). Даже минимум велик: для восьми спинов
$`\frac12`$ ($`n=256`$) $`\sum_ck_cm_c^2=4900`$, а сокращённая система — матрица $`(n^2+N_R)\times N_R`$ — содержит $`3.5\cdot10^8`$
комплексных элементов. Поэтому для спиновых цепочек в разделе 4.6 используется модуль `su2`, которому ответ известен заранее.

```python
for n in (3, 4, 5):
    Dn = dc.from_generators(dc.total_spin(n))
    _, info = dc.commutant_basis(dc.total_spin(n), return_info=True)
    N_min = sum(k * m * m for k, m in zip(Dn.k_factors, Dn.multiplicities))
    assert info['unknowns'] >= N_min                                  # N_R >= sum_c k_c m_c^2
    print(n, info['unknowns'], N_min)                                # observed: equality
```

<a id="s-4-4"></a>

## 4.4 $`U(1)`$ и зарядовые секторы

*Конспект:* [раздел 4.4](../../lecture_notes/ru/04_symmetries.md#s-4-4); [предложение 4.10](../../lecture_notes/ru/04_symmetries.md#r-4-10).

Для $`n`$ кубитов зарядом служит оператор числа частиц $`\hat N=\sum_jn_j`$, $`n_j=|1\rangle\langle1|`$ на кубите $`j`$, а операторы,
сохраняющие его, образуют $`U(1)'\cong\bigoplus_{w=0}^nM_{\binom nw}(\mathbb C)`$ ([конспект лекций, раздел 4.4](../../lecture_notes/ru/04_symmetries.md#s-4-4)),
а это `chains.SpinChain(n).algebra` (бит $`1`$ = спин вверх, $`\hat N=S^z+\frac n2`$; в `decompose.total_spin` спин вверх — это,
напротив, *первый* базисный вектор $`\mathbb C^2`$, что переименовывает базис, но не меняет секторы). Всякий XXZ-гамильтониан
лежит в $`U(1)'`$. Код находит секторы по одному лишь $`\hat N`$ и диагонализует XXZ-кольцо из четырёх спинов по секторам.

```python
import functools
from torch_vn_algebra import SpinChain

L = 4
Sx, Sy, Sz = dc.total_spin(L)                              # (16, 16) each, site 0 most significant
N_op = Sz + L / 2 * torch.eye(2 ** L, dtype=dt)            # number of up spins
DU = dc.from_generators(N_op)                              # the algebra generated by U(1)
print(DU)
chain = SpinChain(L, boundary='periodic')                  # U(1)' = (+)_w M_binom(L, w)
assert set(DU.k_factors) == {1} and sorted(DU.multiplicities) == sorted(chain.algebra.k_factors)
assert sorted(chain.algebra.k_factors) == sorted(math.comb(L, w) for w in range(L + 1))

sx, sy, sz = dc.spin_matrices(0.5)
def two_site(a, b, i, j, n=L):                             # a on site i, b on site j of n sites
    ops = [a if k == i else b if k == j else torch.eye(2, dtype=dt) for k in range(n)]
    return functools.reduce(torch.kron, ops)

H = sum(two_site(sx, sx, i, j) + two_site(sy, sy, i, j) + 0.5 * two_site(sz, sz, i, j)
        for i, j in chain.bonds)                           # XXZ ring, Delta = 0.5
assert close(H @ N_op, N_op @ H)                           # [H, N] = 0 ...
assert close(DU.from_commutant(DU.to_commutant(H)), H)     # ... so H is block diagonal (Proposition 4.10)
Hc = chain.xxz(J=1.0, Delta=0.5)                           # the same H, built sector by sector
E_sec = torch.cat([torch.linalg.eigvalsh(Hc.matrix[0, N, :k, :k]) for N, k in enumerate(chain.algebra.k_factors)])
assert close(torch.sort(E_sec).values, torch.linalg.eigvalsh(H))
```

<a id="s-4-5"></a>

## 4.5 $`SU(2)`$: спин

*Конспект:* [раздел 4.5](../../lecture_notes/ru/04_symmetries.md#s-4-5); [теорема 4.11](../../lecture_notes/ru/04_symmetries.md#r-4-11), [предложение 4.13](../../lecture_notes/ru/04_symmetries.md#r-4-13), [теорема 4.14](../../lecture_notes/ru/04_symmetries.md#r-4-14). Спин $`s`$ — это `decompose.spin_matrices(s)`, а полный спин $`L`$ спинов $`\frac12`$ — это `decompose.total_spin(L)`.

При $`L=4`$ [конспект лекций, предложение 4.13](../../lecture_notes/ru/04_symmetries.md#r-4-13) даёт $`\{S_x,S_y,S_z\}''\cong M_5\oplus M_3\oplus M_1`$
с кратностями $`1,3,2`$, размерности $`35`$, и коммутант $`\{S_x,S_y,S_z\}'\cong M_1\oplus M_3\oplus M_2`$ размерности $`14`$.
Библиотека находит то же самое по одним лишь трём образующим (секторы упорядочены по $`(k_c,m_c)`$), а также
$`\mathbf S^2=S(S+1)`$ на спине $`S`$ ([теорема 4.11](../../lecture_notes/ru/04_symmetries.md#r-4-11)):

```python
S4 = dc.total_spin(4)
assert close(S4[0] @ S4[1] - S4[1] @ S4[0], 1j * S4[2])                # [S_x, S_y] = i S_z
D4 = dc.from_generators(S4)
print(D4)
assert D4.k_factors == (1, 3, 5) and D4.multiplicities == (2, 3, 1)     # M_1, M_3, M_5 with m = 2, 3, 1
assert D4.verify()['ok'] and D4.dims == (35, 14)
m_S = {S: math.comb(4, 2 - S) - (math.comb(4, 1 - S) if S < 2 else 0) for S in (0, 1, 2)}
assert [m_S[(k - 1) // 2] for k in D4.k_factors] == list(D4.multiplicities)    # Proposition 4.13
casimir = D4.to_algebra((S4 @ S4).sum(0))                                  # S^2 = S(S+1) on spin S
for c, k in enumerate(D4.k_factors):
    assert close(casimir.matrix[0, c, :k, :k], (k - 1) / 2 * (k + 1) / 2 * torch.eye(k))
```

Согласно [конспекту лекций, теорема 4.14](../../lecture_notes/ru/04_symmetries.md#r-4-14) (двойственность Шура–Вейля), перестановки узлов
порождают коммутант полного спина: $`\mathrm{span}\,\pi(S_L)=\{S_x,S_y,S_z\}'`$. Для $`L\le5`$ код сравнивает оба
пространства численно: в силу части (a) теоремы, $`\mathrm{span}\,\pi(S_L)\subseteq\{S_x,S_y,S_z\}'`$, достаточно, чтобы они
имели одинаковую размерность.

```python
for n in range(2, 6):
    comm = dc.commutant_basis(dc.total_spin(n)).reshape(-1, 4 ** n)          # a basis of {S_a}'
    perm = dc.permutation_representation(n, 2).reshape(-1, 4 ** n)          # all n! permutations
    dims = [torch.linalg.matrix_rank(A).item() for A in (comm, perm, torch.cat([comm, perm]))]
    print(n, dims)
    assert dims[0] == dims[1] == dims[2]                                     # span pi(S_n) = {S_a}'
Dp = dc.from_group(dc.permutation_representation(4, 2))                   # pi(S_4)'' = {S_a}'
print(Dp)
assert sorted(zip(Dp.multiplicities, Dp.k_factors)) == sorted(zip(D4.k_factors, D4.multiplicities))
```

Напечатанные размерности — $`2,5,14,42`$, числа Каталана $`\binom{2L}L/(L+1)`$; следствие B.10 доказывает
$`\dim\{S_a\}'=\sum_Sm_S^2=\binom{2L}L/(L+1)`$ для всех $`L`$. Последняя строка — двойственность
в действии: у $`\pi(S_4)''`$ блоки $`k=(1,2,3)`$ с кратностями $`(5,1,3)`$ — пары $`(k_c,m_c)`$ полного спина с
переставленными ролями. В общем случае (следствие B.9) пространство кратности спина $`S`$ несёт
неприводимое представление $`S_L`$ размерности $`m_S`$, и эти представления попарно неэквивалентны; здесь
размерности равны $`1,3,2`$ для $`S=2,1,0`$.

<a id="s-4-6"></a>

## 4.6 Спиновые системы

*Конспект:* [раздел 4.6](../../lecture_notes/ru/04_symmetries.md#s-4-6); [предложение 4.15](../../lecture_notes/ru/04_symmetries.md#r-4-15), [теорема 4.16](../../lecture_notes/ru/04_symmetries.md#r-4-16), [следствие 4.17](../../lecture_notes/ru/04_symmetries.md#r-4-17).

Согласно [конспекту лекций, следствие 4.17](../../lecture_notes/ru/04_symmetries.md#r-4-17), ортонормированные базисы векторов старшего веса
$`L`$ спинов $`\frac12`$ можно строить узел за узлом с коэффициентами Клебша–Гордана из [теоремы 4.16](../../lecture_notes/ru/04_symmetries.md#r-4-16).
Пометив каждый базисный вектор последовательностью промежуточных спинов $`S_1=\frac12,S_2,\dots,S_L=S`$, получаем базис
Яманучи, который строит `su2.SU2Chain` (его исходный код использует ровно эти коэффициенты; код ниже проверяет лишь выводы
следствия 4.17), внутри $`U(1)`$-сектора $`\hat N=\frac L2+S`$ класса `SpinChain` (чей порядок узлов — перенумерация тензорных
сомножителей, коммутирующая с полным спином по лемме B.1(b)). Код проверяет теорему 4.16 на плотных матрицах и
сравнивает с `su2.clebsch_gordan` библиотеки, а затем проверяет выводы следствия 4.17 для базисов восьми спинов. Кроме
того, мы проверяем [предложение 4.15(a)](../../lecture_notes/ru/04_symmetries.md#r-4-15) на кольце из четырёх спинов с помощью `D4` и сравниваем
размеры блоков для восьми спинов.

```python
from torch_vn_algebra import su2

def ket(j, m):                                            # |j, m> of Theorem 4.11 (zero if |m| > j)
    v = torch.zeros(int(2 * j + 1), dtype=dt)
    if abs(m) <= j:
        v[int(j - m)] = 1
    return v

up, dn = ket(0.5, 0.5), ket(0.5, -0.5)
for j in (0.5, 1.0, 1.5, 2.0):                            # Theorem 4.16
    Sj = dc.spin_matrices(j)
    sp_j = [math.sqrt(j * (j + 1) - (j - k) * (j - k + 1)) for k in range(1, int(2 * j) + 1)]
    assert close(Sj[0] + 1j * Sj[1], torch.diag(torch.tensor(sp_j, dtype=dt), 1))     # S_+ of Theorem 4.11
    assert close(Sj[2], torch.diag(torch.tensor([j - k for k in range(int(2 * j) + 1)], dtype=dt)))
    Sc = torch.stack([torch.kron(a, torch.eye(2, dtype=dt)) + torch.kron(torch.eye(len(a), dtype=dt), b)
                      for a, b in zip(Sj, dc.spin_matrices(0.5))])
    Sp, Sm = Sc[0] + 1j * Sc[1], Sc[0] - 1j * Sc[1]
    cols = []
    for J in (j + 0.5, j - 0.5):
        for M in [J - k for k in range(int(2 * J + 1))]:
            al, be = math.sqrt((j + M + 0.5) / (2 * j + 1)), math.sqrt((j - M + 0.5) / (2 * j + 1))
            p, q = (al, be) if J > j else (-be, al)
            x = p * torch.kron(ket(j, M - 0.5), up) + q * torch.kron(ket(j, M + 0.5), dn)
            assert close(Sc[2] @ x, M * x)
            if M == J:
                assert close(Sp @ x, 0)                               # highest weight
            else:
                assert close(Sm @ cols[-1], math.sqrt((J + M + 1) * (J - M)) * x)   # lowering, (c)
            for s2, e in ((0.5, up), (-0.5, dn)):                     # the library's Racah formula agrees
                cg = su2.clebsch_gordan(j, M - s2, 0.5, s2, J, M) if abs(M - s2) <= j else 0.0
                assert abs(cg - torch.vdot(torch.kron(ket(j, M - s2), e), x).real.item()) < 1e-12
            cols.append(x)
    U = torch.stack(cols, 1)
    assert U.shape == (2 * len(Sj[0]), 2 * len(Sj[0])) and close(U.mH @ U, torch.eye(U.shape[1]))   # (a)

def swap(i, j, n):
    p = list(range(n))
    p[i], p[j] = j, i
    return dc.permutation_representation(n, 2, [p])[0]

ring = [(i, (i + 1) % 4) for i in range(4)]
H4 = sum(swap(i, j, 4) / 2 - torch.eye(16, dtype=dt) / 4 for i, j in ring)        # Dirac's identity
assert close(H4, sum(two_site(a, a, i, j) for a in (sx, sy, sz) for i, j in ring))
h = D4.to_commutant(H4)                                   # blocks H_S of sizes m_S = 2, 3, 1
assert close(D4.from_commutant(h), H4)                    # H4 lies in {S_a}'
levels = torch.cat([torch.linalg.eigvalsh(h.matrix[0, c, :m, :m]).repeat_interleave(k)
                    for c, (k, m) in enumerate(zip(D4.k_factors, D4.multiplicities))])
assert close(torch.sort(levels).values, torch.linalg.eigvalsh(H4))   # level of H_S repeated 2S+1 times

u1, sym = SpinChain(8, boundary='periodic'), su2.SU2Chain(8, boundary='periodic')
for S, m in zip(sym.spins, sym.multiplicities):           # Corollary 4.17: orthonormal highest-weight bases
    B, N = sym.hw_basis(S), int(4 + S)                     # columns in the U(1) sector N = L/2 + S
    assert B.shape == (math.comb(8, N), m) and close(B.T @ B, torch.eye(m, dtype=B.dtype))
    assert N == 8 or close(sym.raise_(B, N), 0)            # S_+ B = 0 (S_z = S holds in sector N)
H_u1, H_su2 = u1.xxz(J=1.0, Delta=1.0), sym.heisenberg(J=1.0)
E_u1 = torch.sort(torch.cat([torch.linalg.eigvalsh(H_u1.matrix[0, N, :k, :k])
                             for N, k in enumerate(u1.algebra.k_factors)])).values
E_su2, S_lev, _ = sym.full_spectrum(H_su2, return_quantum_numbers=True)   # level of block S, 2S+1 times
E_su2, S_lev = E_su2[0], S_lev[0]
assert close(E_u1, E_su2)
assert S_lev[0] == 0 and (S_lev[1:4] == 1).all() and E_su2[4] - E_su2[3] > 1e-3     # singlet, then triplet
assert E_su2[1] - E_su2[0] > 1e-3 and close(E_su2[1:4], E_su2[1].expand(3))       # simple singlet (Theorem 9.9)
assert close(sym.restrict(H_u1).matrix, H_su2.matrix)      # Proposition 4.15(b)
flops = lambda ks: sum(k ** 3 for k in ks)
print(u1.algebra.k_factors, sym.multiplicities, flops(u1.algebra.k_factors), flops(sym.multiplicities))
print(f"ground state energy {E_su2[0].item():.6f}, first excited {E_su2[1].item():.6f} (x3)")
```

Для кольца из восьми спинов блоки $`U(1)`$ имеют размеры $`\binom8w=1,8,28,56,70,56,28,8,1`$ (раздел 4.4), блоки $`SU(2)`$ —
размеры $`m_S=14,28,20,7,1`$ для $`S=0,\dots,4`$ (предложение 4.13). Стоимость диагонализации, $`\sum k^3`$, падает с
$`739162`$ до $`33040`$, в $`22`$ раза; размерность алгебры, $`\sum k^2`$, — с
$`\sum_w\binom8w\binom8{8-w}=\binom{16}8=12870`$ (выбор $`8`$ из $`16`$ элементов, $`w`$ из них среди первых $`8`$) до
$`\sum_Sm_S^2=1430=\binom{16}8/9`$ (следствие B.10). По теореме 9.9 (кольцо чётной длины и $`J>0`$) основное
состояние — единственный синглет и $`E_0(0)<E_0(1)`$. *Численное наблюдение* (вычислено кодом выше): энергия основного
состояния равна $`-3.651093`$, а первый возбуждённый уровень — триплет при $`-3.128419`$.

<a id="s-4-7"></a>

## 4.7 Вещественные алгебры: вещественный, комплексный и кватернионный тип

*Конспект:* [раздел 4.7](../../lecture_notes/ru/04_symmetries.md#s-4-7). Для $`*`$-алгебры $`\mathcal A`$, замкнутой относительно комплексного
сопряжения, с разложением $`V`$, столбцовые блоки $`V_c`$ которого принадлежат секторам, [теорема 4.20](../../lecture_notes/ru/04_symmetries.md#r-4-20)
даёт перестановку $`\sigma`$ секторов с $`\bar P_c=P_{\sigma(c)}`$, унитарные $`J_c`$, $`M_c`$ с $`\overline{V_c}=V_{\sigma(c)}(J_c\otimes M_c)`$
и, при $`\sigma(c)=c`$, знак $`\varepsilon_c`$ с $`J_c\bar J_c=\varepsilon_c1`$; типы и вещественная алгебра $`\mathcal A_{\mathbb R}`$ описаны
в [теореме 4.21](../../lecture_notes/ru/04_symmetries.md#r-4-21), вещественная лемма Шура — это [теорема 4.22](../../lecture_notes/ru/04_symmetries.md#r-4-22), а вещественные
элементы коммутанта описывает [предложение 4.23](../../lecture_notes/ru/04_symmetries.md#r-4-23).

`Decomposition.real_structure()` вычисляет $`\sigma`$, $`J_c`$, $`M_c`$ и типы для разложения, возвращённого
`wedderburn` (`from_generators`, `from_group`). Она проверяет, что сопряжённые к хранимым образующим лежат в
$`V(\bigoplus_cM_{k_c}\otimes1)V^*`$ (относительная невязка `'conjugation'`), и строит унитарную $`P=V^*\bar V`$, чей блок столбцов
$`c`$ есть $`V^*\overline{V_c}`$: по теореме 4.20(a) он сосредоточен на строках сектора $`\sigma(c)`$, где равен
$`J_c\otimes M_c`$. Библиотека берёт в качестве $`\sigma(c)`$ сектор, несущий наибольшую часть блока, требует
$`(k_{\sigma(c)},m_{\sigma(c)})=(k_c,m_c)`$ и оценивает остаток (`'leak'`); она переставляет блок в
$`k_c^2\times m_c^2`$-матрицу с элементами $`T_{(a,b),(j,l)}=T_{(a,j),(b,l)}`$, которая для $`T=J\otimes M`$ есть матрица ранга один
$`\mathrm{vec}\,J\,(\mathrm{vec}\,M)^T`$, считывает $`J_c`$ и $`M_c`$ с её старших сингулярных векторов (с
$`\|J_c\|_F^2=k_c`$) и записывает расстояние блока от $`J_c\otimes M_c`$ (`'kron'`). Если одна из этих трёх невязок
превышает допуск, размеры различаются или партнёры не образуют инволюцию, она выбрасывает `ValueError` (алгебра
не замкнута относительно сопряжения). При $`\sigma(c)=c`$ она полагает $`\varepsilon_c`$ равным знаку $`\mathrm{Re}\,\mathrm{Tr}(J_c\bar J_c)/k_c`$ —
в точной арифметике это ровно $`\pm1`$, так что два кандидата отстоят на $`2`$, — и записывает как
`'sign'` максимум по этим секторам величин $`\|J_c\bar J_c-\varepsilon_c\|_F/\sqrt{k_c}`$ и
$`\|M_c\bar M_c-\varepsilon_c\|_F/\sqrt{m_c}`$, увеличенный по меньшей мере до $`1`$, если некоторое $`\varepsilon_c=-1`$ при нечётном $`k_c`$ или $`m_c`$
(что невозможно по теореме 4.20(b)). Результат имеет `types`, `partner` ($`\sigma`$),
`J`, `M`, `residuals` и `ok`, означающее, что все невязки не превосходят допуска (по умолчанию — порога
принятия разложения); тогда и $`|\mathrm{Tr}(J_c\bar J_c)/k_c-\varepsilon_c|\le\|J_c\bar J_c-\varepsilon_c\|_F/\sqrt{k_c}`$
(Коши–Буняковский–Шварц) не превосходит допуска.

`Decomposition.real_basis()` возвращает вещественную ортогональную $`n\times n`$-матрицу $`O`$ и список, описывающий её блоки столбцов (тип, секторы,
диапазон столбцов): для
сектора вещественного типа — столбцы $`O_c=V_c(W_J\otimes W_M)`$ из теоремы 4.21(a), так что $`O_c^TaO_c=\alpha_c(a)\otimes1_{m_c}`$
с вещественной $`\alpha_c(a)`$; для сектора кватернионного типа — столбцы $`V_cW`$ с $`J_c\otimes M_c=WW^T`$ и
$`(J_c\otimes M_c)\bar W=W`$ (лемма 4.19(c) с $`\varepsilon=1`$, так как $`(J_c\otimes M_c)\overline{(J_c\otimes M_c)}=\varepsilon_c^2=1`$),
вещественные по вычислению из теоремы 4.21(a); и для пары $`c,\sigma(c)`$ (с $`c<\sigma(c)`$) — столбцы $`\sqrt2\mathrm{Re}\,V_c`$, за которыми следуют
$`\sqrt2\mathrm{Im}\,V_c`$. Они ортонормированы: $`V_c^*\overline{V_c}=V_c^*V_{\sigma(c)}(J_c\otimes M_c)=0`$ и
$`V_c^T\overline{V_c}=\overline{V_c^*V_c}=1`$ вместе с их сопряжёнными дают $`(\mathrm{Re}\,V_c)^T\mathrm{Re}\,V_c=(\mathrm{Im}\,V_c)^T\mathrm{Im}\,V_c=\frac12`$ и
$`(\mathrm{Re}\,V_c)^T\mathrm{Im}\,V_c=0`$; столбцы разных блоков лежат в образах
ортогональных проекторов $`P_{\mathcal O}=\sum_{e\in\mathcal O}P_e`$ разных орбит $`\mathcal O`$ перестановки $`\sigma`$ (так как
$`\overline{V_c}=V_{\sigma(c)}(J_c\otimes M_c)`$), поэтому $`O`$ ортогональна. Беря вещественную и мнимую части $`aV_c=V_cx`$, $`x=\pi_c(a)\otimes1_{m_c}`$,
для вещественной $`a`$, получаем
```math
O^TaO=\begin{pmatrix}\mathrm{Re}\,x&\mathrm{Im}\,x\\-\mathrm{Im}\,x&\mathrm{Re}\,x\end{pmatrix}\quad\text{on these columns.}
```
Множители Такаги $`W`$ с $`WW^T=J`$ симметричной унитарной $`J=X+iY`$ вычисляются по ортогональной
диагонализации $`X+tY`$ для случайного $`t`$ ($`X`$ и $`Y`$ — коммутирующие вещественные симметричные матрицы, так как мнимая часть
$`J\bar J=JJ^*=1`$ равна $`YX-XY`$); библиотека принимает результат по невязкам $`\|WW^T-J\|_F/\|J\|_F`$ и
$`\|\mathrm{Im}(V_cW)\|_F/\sqrt{k_cm_c}`$ и выдаёт предупреждение (`DecompositionWarning`), если одна из них превышает
$`\max(\text{tolerance},10^3\,\epsilon_{\rm mach}\,n)`$.

В [разделе 4.7](../../lecture_notes/ru/04_symmetries.md#s-4-7) конспекта лекций разобраны три алгебры, порождённые вещественными матрицами:
(i) $`L_a\otimes1_2`$ на $`\mathbb R^3\otimes\mathbb R^2`$ с образующими $`(L_a)_{bc}=-\epsilon_{abc}`$ поворотов $`\mathbb R^3`$: один сектор
вещественного типа, $`k=3`$, $`m=2`$; (ii) $`\Omega_2\otimes1_2`$ на $`\mathbb R^2\otimes\mathbb R^2`$: пара комплексного типа, $`k=1`$,
$`m=2`$; (iii) левые умножения $`1_2\otimes L_q`$ на кватернионы на $`\mathbb H^2=\mathbb R^8`$: один сектор кватернионного типа,
$`k=2`$, $`m=4`$. Код вычисляет разложение и вещественную структуру каждой из них и проверяет теорему 4.20(a),(b) и блочный вид
$`O^TaO`$ для вещественного типа (теорема 4.21(a)) и для пары.

```python
def qmul(p, q):                                           # H = CD(C), Section E.6: (a, b)(c, d)
    (a1, b1), (a2, b2) = p, q
    return (a1 * a2 - b2.conjugate() * b1, b2 * a1 + b1 * a2.conjugate())

def q_left(q):                                            # left multiplication by q on H = R^4
    cols = []                                             # coordinates (Re a, Im a, Re b, Im b)
    for e in torch.eye(4, dtype=torch.float64).tolist():
        u, v = qmul(q, (complex(e[0], e[1]), complex(e[2], e[3])))
        cols.append([u.real, u.imag, v.real, v.imag])
    return torch.tensor(cols, dtype=torch.float64).T

L_so3 = torch.zeros(3, 3, 3, dtype=torch.float64)         # (L_a)_bc = -epsilon_abc
for i, (b1, b2) in enumerate([(1, 2), (2, 0), (0, 1)]):
    L_so3[i, b1, b2], L_so3[i, b2, b1] = -1.0, 1.0
R_rot = torch.tensor([[[0.0, -1.0], [1.0, 0.0]]])         # Omega_2
L_H = torch.stack([q_left(u) for u in ((1j, 0), (0, 1), (0, 1j))])            # L_I, L_J, L_K
assert close(L_H[0] @ L_H[1], L_H[2]) and close(q_left((2 - 1j, 3j)).T, q_left((2 + 1j, -3j)))   # IJ = K, L_qbar = L_q^T
eye2 = torch.eye(2, dtype=torch.float64)
real_examples = {'R': (torch.kron(L_so3, eye2[None]), (3,), (2,)),            # L_a (x) 1_2 on R^3 (x) R^2
                 'C': (torch.kron(R_rot, eye2[None]), (1, 1), (2, 2)),        # Omega_2 (x) 1_2
                 'H': (torch.kron(eye2[None], L_H), (2,), (4,))}              # L_q on H^2 = R^8
dec_real = {}
for t, (gens, ks, ms) in real_examples.items():
    D_t = dc.from_generators(gens)
    rs_t, (O_t, blocks_t) = D_t.real_structure(), D_t.real_basis()
    print(t, D_t, rs_t)
    assert D_t.verify()['ok'] and rs_t.ok and (D_t.k_factors, D_t.multiplicities) == (ks, ms)
    assert set(rs_t.types) == {t} and close(O_t.T @ O_t, torch.eye(D_t.n))     # real orthogonal basis
    for c, d in enumerate(rs_t.partner):                                      # Theorem 4.20
        o_c, o_d, w = D_t.offsets[c], D_t.offsets[d], ks[c] * ms[c]
        assert close(D_t.V[:, o_c:o_c + w].conj(), D_t.V[:, o_d:o_d + w] @ torch.kron(rs_t.J[c], rs_t.M[c]))
        if d == c:                                                            # J conj(J) = M conj(M) = +-1
            eps = 1.0 if t == 'R' else -1.0
            assert close(rs_t.J[c] @ rs_t.J[c].conj(), eps * torch.eye(ks[c]))
            assert close(rs_t.M[c] @ rs_t.M[c].conj(), eps * torch.eye(ms[c]))
    dec_real[t] = D_t, rs_t, O_t

D_t, rs_t, O_t = dec_real['R']
for g in real_examples['R'][0]:                           # Theorem 4.21(a): O^T a O = alpha (x) 1_2, alpha real
    G_t = O_t.T @ g @ O_t
    assert close(G_t, torch.kron(G_t[0::2, 0::2], eye2))
D_t, rs_t, O_t = dec_real['C']                            # the pair: [[Re x, Im x], [-Im x, Re x]]
x_t = torch.kron(D_t.to_algebra(real_examples['C'][0]).matrix[0, 0, :1, :1], torch.eye(2, dtype=dt))
assert close(O_t.T @ real_examples['C'][0][0] @ O_t,
             torch.cat([torch.cat([x_t.real, x_t.imag], 1), torch.cat([-x_t.imag, x_t.real], 1)]))
```

Код проверяет теорему 4.22 на неприводимых вещественных модулях $`\mathbb R^3`$, $`\mathbb R^2`$, $`\mathbb R^4`$ трёх
примеров (типы и кратности, найденные библиотекой, — вещественный с $`m=1`$, комплексная пара с $`m=1`$, кватернионный
с $`m=2`$ — делают их неприводимыми по теореме 4.22(c)), вычисляя вещественный
коммутант как ядро вещественной линейной системы $`Xs=sX`$, $`Xs^T=s^TX`$, а затем предложение 4.23 для случайного
вещественного симметричного элемента коммутанта в каждом примере.

```python
def real_commutant_dim(gens):                             # dim_R {X real : X s = s X, X s^T = s^T X}
    gens = torch.cat([gens, gens.mT])
    eye_n = torch.eye(gens.shape[-1], dtype=torch.float64)
    eqs = torch.cat([torch.kron(eye_n, s) - torch.kron(s.T.contiguous(), eye_n) for s in gens])
    return gens.shape[-1] ** 2 - torch.linalg.matrix_rank(eqs).item()

for gens, types, ms, dim_D in ((L_so3, ('R',), (1,), 1), (R_rot, ('C', 'C'), (1, 1), 2), (L_H, ('H',), (2,), 4)):
    D_irr = dc.from_generators(gens)                      # Theorem 4.22(c): R^3, R^2, R^4 are irreducible
    assert D_irr.real_structure().types == types and D_irr.multiplicities == ms
    assert real_commutant_dim(gens) == dim_D              # Theorem 4.22(b): D_U = R, C, H

gen_X = torch.Generator().manual_seed(1)
for t, (gens, ks, ms) in real_examples.items():           # Proposition 4.23
    D_t, rs_t, O_t = dec_real[t]
    B_t = dc.commutant_basis(gens)
    X_t = torch.einsum('r,rij->ij', torch.randn(len(B_t), dtype=dt, generator=gen_X), B_t).real
    X_t = X_t + X_t.T                                     # real symmetric, in A' (Lemma 4.18(b))
    assert close(X_t @ gens, gens @ X_t)
    x_t = D_t.to_commutant(X_t).matrix[0]
    ev = [torch.linalg.eigvalsh(x_t[c, :m, :m]) for c, m in enumerate(ms)]
    if t == 'R':                                          # (a): O^T X O = 1_3 (x) y, y real symmetric
        Y_t = O_t.T @ X_t @ O_t
        assert close(Y_t, torch.kron(torch.eye(3, dtype=torch.float64), Y_t[:2, :2]))
    if t == 'C':                                          # (c): the two sectors have the same spectrum
        assert close(ev[0], ev[1])
    if t == 'H':                                          # (b): Kramers pairs, each level 2 k_c = 4 times in X
        assert close(ev[0][0::2], ev[0][1::2])
        assert close(torch.linalg.eigvalsh(X_t), ev[0].repeat_interleave(2))
        print('quaternionic block x_c:', [round(v, 4) for v in ev[0].tolist()])
```

<a id="s-4-8"></a>

## 4.8 С чем сравнивать

Для конечных групп GAP и SageMath точно вычисляют таблицы характеров, и тогда следствие 4.5 даёт
кратности; RepLAB (MATLAB; D. Rosset, F. Montealegre-Mora, J.-D. Bancal) численно раскладывает
представления конечных и компактных групп, тоже с помощью случайных элементов коммутанта, а
SymbolicWedderburn.jl (Julia) вычисляет разложения Веддерберна для симметрийной редукции
полуопределённых программ. В физике QuSpin и другие программы точной диагонализации приводят к блочно-
диагональному виду $`U(1)`$ и решёточные симметрии, ITensor и TeNPy работают с абелевыми зарядами в тензорных
сетях, а TensorKit.jl обрабатывает $`SU(2)`$ с помощью деревьев слияния. Модуль `decompose` библиотеки — плотный
метод для небольших матриц ($`n`$ до нескольких сотен); он принимает любой список матриц (элементы группы
или образующие), сертифицирует свой результат (предложение 4.9) и возвращает алгебру и её коммутант как
объекты `TypeIAlgebra` с `to_algebra`, `from_algebra`, `to_commutant`, `from_commutant`, так что главы 3, 5–8
применимы сразу. Для цепочек спинов $`\frac12`$ `su2.SU2Chain` даёт блоки $`SU(2)`$ напрямую (следствие 4.17); он принимает
до $`L=18`$ спинов — предел, определяемый памятью для базисов.

<a id="s-4-9"></a>

## 4.9 Устойчивость численного разложения

*Конспект:* [раздел 4.9](../../lecture_notes/ru/04_symmetries.md#s-4-9); [лемма 4.24](../../lecture_notes/ru/04_symmetries.md#r-4-24), [теорема 4.25](../../lecture_notes/ru/04_symmetries.md#r-4-25), [предложение 4.26](../../lecture_notes/ru/04_symmetries.md#r-4-26), [предложение 4.27](../../lecture_notes/ru/04_symmetries.md#r-4-27), [теорема 4.28](../../lecture_notes/ru/04_symmetries.md#r-4-28).

Коммутант зашумлённых образующих определяется допуском: [теорема 4.25](../../lecture_notes/ru/04_symmetries.md#r-4-25) утверждает, что порог $`\theta`$ для сингулярных
чисел коммутаторного отображения $`X\mapsto([H_j,X])_j`$ даёт точную размерность $`d`$ коммутанта, если $`2\eta<\theta<\gamma-2\eta`$, где
$`\eta`$ — шум, а $`\gamma`$ — щель точного отображения, и что оценка $`2\eta`$ для возмущённых нулевых сингулярных чисел достигается. Первый
блок вычисляет сингулярные числа примера из теоремы 4.25(c), $`\tilde H_1=\mathrm{diag}(1,1,0)+\eta(E_{12}+E_{21})`$, и определяет коммутант
с помощью `commutant_basis` ниже окна, внутри него и выше (библиотека нормирует образующие на единичную норму Фробениуса, отсюда
перемасштабированный `tol`).

```python
import warnings

def commutator_svals(Hs):
    """Singular values, in increasing order, of X -> ([H_j, X])_j with Hilbert-Schmidt norms (row-major vec)."""
    n = Hs.shape[-1]
    eye = torch.eye(n, dtype=dt)
    L = torch.cat([torch.kron(H, eye) - torch.kron(eye, H.T.contiguous()) for H in Hs])
    return torch.linalg.svdvals(L).flip(0)

eta = 1e-3
H1 = torch.diag(torch.tensor([1.0, 1.0, 0.0], dtype=dt))
H1t = H1.clone()
H1t[0, 1] = H1t[1, 0] = eta                               # a symmetry broken at the level eta
s = commutator_svals(H1t[None])
expected = torch.tensor([0, 0, 0, 2 * eta, 2 * eta, 1 - eta, 1 - eta, 1 + eta, 1 + eta], dtype=torch.float64)
assert close(s, expected, 1e-12)                          # Theorem 4.25(a) and (c): the bound 2 eta is attained
assert close(commutator_svals(H1[None])[4:6], [0, 1])      # exact generator: d = 5, gap 1
scale = torch.linalg.matrix_norm(H1t).item()              # theta = tol for unit-norm generators
for theta, dim in ((0.5 * eta, 3), (20 * eta, 5), (2.0, 9)):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', dc.DecompositionWarning)
        assert dc.commutant_basis(H1t, tol=theta / scale).shape[0] == dim
print('theta < 2 eta: dim 3 (symmetry lost); in the window: 5; above the gap: 9')
```

Далее — три спиновых оператора трёх спинов $`\frac12`$ ($`k=(2,4)`$, $`m=(2,1)`$, $`d=5`$, щель $`\gamma=1/\sqrt3`$ после нормировки),
возмущённые эрмитовым шумом возрастающей величины $`\eta`$. Для каждого уровня код проверяет теорему 4.25(a) на сингулярных числах
$`\sigma_d`$ и $`\sigma_{d+1}`$, помещает `tol` в середину наблюдаемого окна (в логарифмическом масштабе), проверяет размерность и оценку
$`\|\tilde\Pi-\Pi\|\le4\eta/\gamma`$ теоремы 4.25(b) и выполняет разложение. С `tol` по умолчанию $`=\epsilon_{\rm mach}^{1/2}`$ шум
принимается всерьёз, как только он превышает допуск: зашумлённые образующие порождают $`M_8`$ — это правильный ответ для них, но не та
симметрия, которую ищут. (Отфильтрованное здесь `DecompositionWarning` сообщает, что некоторое сингулярное число отстоит от разреза менее
чем в $`\epsilon_{\rm mach}^{-1/8}`$ раз.)

```python
S = dc.total_spin(3)
unit = lambda S: S / torch.linalg.matrix_norm(S)[:, None, None]     # the Hermitian generators the library uses
sv = commutator_svals(unit(S))
d, gamma = 5, sv[5].item()
assert sv[4] < 1e-12 and abs(gamma - 1 / math.sqrt(3)) < 1e-12
B0 = dc.commutant_basis(S).reshape(d, -1)
Pi0 = B0.T @ B0.conj()                                    # orthogonal projection onto vec(A')
gen = torch.Generator().manual_seed(1)
print(' eta       sigma_d   sigma_d+1  tol       (k, m)            default tol')
for eta in (1e-8, 1e-6, 1e-4, 1e-3):
    E = torch.randn(S.shape, dtype=dt, generator=gen)
    E = E + E.mH
    St = S + eta * E / torch.linalg.matrix_norm(E, ord=2)[:, None, None]
    eta_n = torch.linalg.matrix_norm(unit(St) - unit(S), ord=2).norm().item()   # eta of Theorem 4.25
    st = commutator_svals(unit(St))
    assert st[d - 1] <= 2 * eta_n and st[d] >= gamma - 2 * eta_n                # Theorem 4.25(a)
    tol = math.sqrt(st[d - 1] * st[d]) / max(1.0, st[-1].item())
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', dc.DecompositionWarning)
        Bt = dc.commutant_basis(St, tol=tol).reshape(d, -1)                     # Theorem 4.25(b): dimension d
        Dt = dc.from_generators(St, tol=tol)
        Dd = dc.from_generators(St)
    assert torch.linalg.matrix_norm(Bt.T @ Bt.conj() - Pi0, ord=2) <= 4 * eta_n / gamma
    assert (Dt.k_factors, Dt.multiplicities) == ((2, 4), (2, 1))
    print(f'{eta_n:.1e}  {st[d - 1]:.1e}   {st[d]:.3f}      {tol:.1e}   {Dt.k_factors}, {Dt.multiplicities}'
          f'   {Dd.k_factors}, {Dd.multiplicities}')
```

`Decomposition.certificate()` вычисляет на результате [предложение 4.26](../../lecture_notes/ru/04_symmetries.md#r-4-26): обратную ошибку $`r`$ (расстояние нормированных
образующих до вычисленной алгебры $`\tilde{\mathcal B}`$), щель $`\gamma_{\mathcal B}`$ коммутаторного отображения их проекций $`H_j=E(\tilde H_j)`$ —
вычисленную по блокам и сверенную здесь с плотными сингулярными числами — и окно $`(2r,\gamma_{\mathcal B}-2r)`$ порогов, сертифицирующих
$`\dim\mathcal A'=\sum_cm_c^2`$. Результат `trusted`, если фактически использованный порог лежит в этом окне; иначе сертификат выдаёт
предупреждение. Допуск ниже шума даёт $`M_8`$, сертифицированную для этих образующих, но со щелью порядка шума: модель с большей
симметрией может находиться на расстоянии $`\gamma_{\mathcal B}/2`$. Допуск выше шума, но ниже $`2r`$ даёт правильную структуру, за которую
сертификат поручиться не может.

```python
c = dc.from_generators(S).certificate()
assert c.certified and c.trusted and c.residuals['backward'] < 1e-12 and abs(c.gap - gamma) < 1e-10
E = torch.randn(S.shape, dtype=dt, generator=gen)
St = S + 1e-4 * (E + E.mH)
with warnings.catch_warnings():
    warnings.simplefilter('ignore', dc.DecompositionWarning)
    Dt = dc.from_generators(St, tol=1e-2)
    D_low = dc.from_generators(St, tol=1e-12)                          # below the noise
    D_mid = dc.from_generators(St, tol=1e-3)                           # above the noise, below 2 r
ct = Dt.certificate()
print(ct, 'window', [f'{w:.1e}' for w in ct.window], 'threshold', f"{ct.rank['threshold']:.1e}")
assert ct.trusted and 1e-5 < ct.residuals['backward'] < 1e-2 and ct.window[0] < ct.rank['threshold'] < ct.window[1]
sb = commutator_svals(Dt.conditional_expectation(unit(St)))           # the projected generators H_j
assert sb[d - 1] < 1e-12 and abs(sb[d] - ct.gap) < 1e-10              # Proposition 4.26(a), (b)
c_low = D_low.certificate()
assert D_low.k_factors == (8,) and c_low.certified and c_low.gap < 1e-3
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter('always')
    c_mid = D_mid.certificate()
assert D_mid.k_factors == (2, 4) and not c_mid.trusted and caught
print('tol 1e-12:', D_low.k_factors, f'gap {c_low.gap:.1e};  tol 1e-3:', D_mid.k_factors, 'trusted:', c_mid.trusted)
```

Наконец, [предложение 4.27](../../lecture_notes/ru/04_symmetries.md#r-4-27): для стандартного гауссова вектора $`G`$ в $`\mathcal A'_h`$ блоки равны $`\kappa_c=k_c^{-1/2}H_c`$ с
независимыми матрицами GUE, и вероятность того, что собственные значения разных секторов сблизятся меньше чем на $`x`$, не более чем
порядка $`x`$ (здесь оценка равна $`\sqrt{2/\pi}\,x\,m_1m_2(\frac1{k_1m_1}+\frac1{k_2m_2})^{-1/2}`$), тогда как внутри блока $`2\times2`$
расстояние имеет порядок $`x^3`$. Проверка методом Монте-Карло с $`2\cdot10^4`$ выборками в ортонормированном базисе
`D.basis('commutant')` пространства $`\mathcal A'`$:

```python
D = dc.from_generators(S)
coef = torch.randn(20000, d, dtype=dt, generator=gen)
Z = torch.einsum('br,rij->bij', coef, D.basis('commutant'))         # standard complex Gaussian vectors of A'
x_c = D.to_commutant((Z + Z.mH) / math.sqrt(2)).matrix              # blocks of G: kappa_1 (2 x 2), kappa_2 (1 x 1)
assert abs((x_c[:, 0, :2, :2].abs() ** 2).sum((-2, -1)).mean().item() - 4 / 2) < 0.05   # E Tr kappa_1^2 = m_1^2 / k_1
ev = torch.linalg.eigvalsh(x_c[:, 0, :2, :2])
cross = (ev - x_c[:, 1, :1, 0].real).abs().min(1).values           # delta_x
(k1, k2), (m1, m2) = D.k_factors, D.multiplicities
for x in (0.01, 0.03, 0.1):
    bound = math.sqrt(2 / math.pi) * x * m1 * m2 / math.sqrt(1 / (k1 * m1) + 1 / (k2 * m2))
    p = (cross < x).double().mean().item()
    print(f'x = {x}:  P(delta_x < x) = {p:.4f} <= {bound:.4f}')
    assert p <= bound
p_in = (ev[:, 1] - ev[:, 0] < 0.2).double().mean().item()
assert p_in <= k1 ** 1.5 * 0.2 ** 3 / (6 * math.sqrt(math.pi))      # Proposition 4.27(c), m = 2
```

<a id="x-4"></a>

## Упражнения

[Упражнения](../../lecture_notes/ru/04_symmetries.md#x-4) конспекта лекций, проверенные численно (решения ниже идут в этом порядке):
1. $`\dim\pi(S_L)'=\binom{d^2+L-1}L`$ при $`d=2`$ и $`L=3,4`$ (значения $`20`$ и $`35`$);
2. $`1\otimes\frac12\cong\frac12\oplus\frac32`$ и $`1\otimes1\cong0\oplus1\oplus2`$ с помощью `from_generators`;
3. коммутант $`e^{i\theta\hat N}`$ на четырёх кубитах при $`\theta=1`$ (заряды $`0,\dots,4`$ дают различные $`e^{iq}`$, так как
   $`0<|q-q'|\le4<2\pi`$) и при $`\theta=\pi`$;
4. эффект Зеемана: спектр $`H+hS_z`$ для гейзенберговского кольца из восьми спинов с помощью `su2.SU2Chain` и в
   $`U(1)`$-секторах;
5. секторы регулярного представления $`Q_8`$ с помощью `real_structure`.

<details>
<summary>Решения</summary>

```python
for n, value in ((3, 20), (4, 35)):
    assert dc.commutant_dimension(dc.permutation_representation(n, 2)) == math.comb(n + 3, n) == value
```

```python
def coupled(a, b):                                        # S^(a) (x) 1 + 1 (x) S^(b)
    Sa, Sb = dc.spin_matrices(a), dc.spin_matrices(b)
    ea, eb = torch.eye(Sa.shape[-1], dtype=dt), torch.eye(Sb.shape[-1], dtype=dt)
    return torch.stack([torch.kron(x, eb) + torch.kron(ea, y) for x, y in zip(Sa, Sb)])
for a, b, ks in ((1.0, 0.5, (2, 4)), (1.0, 1.0, (1, 3, 5))):  # 1 (x) 1/2 = 1/2 (+) 3/2, 1 (x) 1 = 0 (+) 1 (+) 2
    D_cg = dc.from_generators(coupled(a, b))
    assert D_cg.k_factors == ks and set(D_cg.multiplicities) == {1}
```

```python
for theta, m in ((1.0, (1, 1, 4, 4, 6)), (math.pi, (8, 8))):   # theta = 1: e^{iq}, q = 0..4, distinct as 4 < 2 pi
    D_theta = dc.from_generators(torch.linalg.matrix_exp(1j * theta * N_op))
    assert D_theta.multiplicities == m and set(D_theta.k_factors) == {1}
```

```python
assert not close(H4 @ Sx - Sx @ H4 + 0.3 * (Sz @ Sx - Sx @ Sz), 0)          # [H4 + 0.3 S_z, S_x] != 0
E_field = sym.full_spectrum(sym.heisenberg(J=1.0, h=0.3))[0]                # E + h M on every multiplet
H_field = u1.xxz(J=1.0, Delta=1.0, h=0.3)
E_u1_field = torch.cat([torch.linalg.eigvalsh(H_field.matrix[0, N, :k, :k])
                        for N, k in enumerate(u1.algebra.k_factors)])
assert close(torch.sort(E_u1_field).values, E_field)
```

```python
Q8 = [(s * a, s * b) for a, b in ((1, 0), (1j, 0), (0, 1), (0, 1j)) for s in (1, -1)]
key = lambda q: tuple(round(v) for z in q for v in (complex(z).real, complex(z).imag))
pos = {key(q): i for i, q in enumerate(Q8)}
reg = torch.zeros(8, 8, 8, dtype=torch.float64)           # pi(g) e_h = e_{gh}
for i, g in enumerate(Q8):
    for j, h in enumerate(Q8):
        reg[i, pos[key(qmul(g, h))], j] = 1.0
D_Q8 = dc.from_group(reg)
rs_Q8 = D_Q8.real_structure()
assert rs_Q8.ok and sorted(zip(rs_Q8.types, D_Q8.k_factors, D_Q8.multiplicities)) == [('H', 2, 2)] + [('R', 1, 1)] * 4
assert D_Q8.dims[0] == 8                                  # dim_R A_R = dim_C A = |Q_8| (Theorem 4.21)
```

</details>

**Далее:** [Руководство 5](05_random_matrices.md): случайные матрицы в алгебре — унитарные матрицы Хаара, исчисление Вайнгартена и статистика уровней в библиотеке.
