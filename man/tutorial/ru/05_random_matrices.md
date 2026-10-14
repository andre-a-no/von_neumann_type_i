[English](../05_random_matrices.md) · [Español](../es/05_random_matrices.md) · [Français](../fr/05_random_matrices.md) · [Deutsch](../de/05_random_matrices.md) · [中文](../zh/05_random_matrices.md) · [日本語](../ja/05_random_matrices.md) · **Русский**

# Руководство 5. Случайные матрицы в алгебре

Это руководство сопровождает [главу 5](../../lecture_notes/ru/05_random_matrices.md) конспекта лекций, в которой строится мера Хаара на унитарной
группе $`\mathcal M=\bigoplus_cM_{k_c}`$ ([следствие 5.4](../../lecture_notes/ru/05_random_matrices.md#r-5-4)), вычисляются её первые и вторые моменты
([предложение 5.7](../../lecture_notes/ru/05_random_matrices.md#r-5-7), [теорема 5.10](../../lecture_notes/ru/05_random_matrices.md#r-5-10)) и объясняется, почему статистику уровней нужно вычислять по
секторам ([теорема 5.15](../../lecture_notes/ru/05_random_matrices.md#r-5-15)). Мы порождаем с помощью библиотеки унитарные матрицы Хаара и гауссовы ансамбли
и проверяем эти результаты методом Монте-Карло: среднее должно совпадать с точным значением с точностью до $`5`$
стандартных ошибок. Результаты, за которыми не стоит теорема, помечены как численные наблюдения.

Всюду далее $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$, $`D=\sum_ck_c`$, $`P_c`$ —
минимальные центральные проекторы, а $`E`$ — условное ожидание на центр из
предложения 2.8. $`U(k)`$ — группа унитарных матриц размера $`k\times k`$, а
$`U(\mathcal M)=\{U\in\mathcal M:U^*U=1\}`$ — унитарная группа $`\mathcal M`$. Случайные величины
заданы на вероятностном пространстве с математическим ожиданием $`\mathbb E`$; математическое ожидание случайной матрицы
берётся поэлементно. У *стандартной комплексной гауссовой* величины вещественная и мнимая части независимы
и $`\sim N(0,\tfrac12)`$, так что $`\mathbb E|z|^2=1`$ (согласно документации, именно такие величины порождает `torch.randn` для комплексного
dtype; первый блок кода проверяет вторые моменты); *матрица Жинибра* имеет независимые стандартные комплексные гауссовы элементы.

<a id="s-5-1"></a>

## 5.1 Мера Хаара

*Конспект:* [раздел 5.1](../../lecture_notes/ru/05_random_matrices.md#s-5-1); [теорема 5.3](../../lecture_notes/ru/05_random_matrices.md#r-5-3) (унитарный множитель QR-разложения матрицы Жинибра с фазами, фиксированными при помощи $`\Lambda=\mathrm{diag}(R'_{ii}/|R'_{ii}|)`$, распределён по Хаару), [следствие 5.4](../../lecture_notes/ru/05_random_matrices.md#r-5-4) (блоки распределённого по Хаару унитарного элемента $`U(\mathcal M)`$ независимы и распределены по Хаару) и [следствие 5.5](../../lecture_notes/ru/05_random_matrices.md#r-5-5) ($`\mathbb EU_{ij}=0`$, $`\mathbb E|U_{ij}|^2=1/k`$). Мера Хаара на $`SU(k)`$ — это распределение $`e^{-i\theta/k}U`$, $`\theta=\arg\det U`$.

Без поправки $`\Lambda`$ из [теоремы 5.3](../../lecture_notes/ru/05_random_matrices.md#r-5-3) распределение $`Q'`$ зависит от соглашения о знаках в процедуре QR.
Для процедуры, вызываемой `torch.linalg.qr` (LAPACK), код ниже находит, что
$`\mathbb E\mathrm{Re}\,\mathrm{Tr}\,Q'`$ отстоит от $`0`$ на много стандартных ошибок, тогда как у унитарной матрицы Хаара
$`\mathbb EU=0`$ (следствие 5.5): значит, эта $`Q'`$ не распределена по Хаару. Это утверждение об этой
процедуре, проверяемое кодом, а не теорема.

Библиотека порождает распределённые по Хаару унитарные элементы $`U(\mathcal M)`$ функцией `random_unitary_operator`: одно пакетное
QR-разложение на каждый размер блока сразу для всех секторов и выборок (на GPU — процесс Грама–Шмидта, который
даёт ту же $`Q`$ в силу единственности в лемме 5.2). Мера `'haar'` — это теорема 5.3 в каждом
блоке; для комплексной алгебры `'haar_su'` домножает на $`e^{-i\theta/k}`$, как выше; `'coe'` — это
ансамбль $`S=W^TW`$ из упражнения 3 (остальные меры библиотеки в этом курсе не используются). Проверим
следствия 5.4 и 5.5, а также то, что сэмплер библиотеки — это в точности исправленное QR-разложение теоремы 5.3:

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, SpinChain, channels

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def mc(samples, exact, name):
    """Monte Carlo mean of real samples with its standard error; must agree within 5 standard errors."""
    mean, err = samples.mean().item(), samples.std().item() / math.sqrt(samples.numel())
    print(f"{name}: {mean:.4f} +- {err:.4f}   (theory {exact:.4f})")
    assert abs(mean - exact) <= 5 * err + 1e-12
    return mean

alg = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')        # M_1 (+) M_2 (+) M_3, D = 6
B = 100_000
U = alg.random_unitary_operator(B)                                   # Haar on U(1) x U(2) x U(3)
assert close(U.adjoint() @ U, alg.identity())                        # unitary elements of the algebra
u = [U.matrix[:, c, :k, :k] for c, k in enumerate(alg.k_factors)]   # the blocks U_c
for c, k in enumerate(alg.k_factors):                                # Corollary 5.5
    mc(u[c][:, 0, 0].real, 0.0, f"E Re (U_{c})_11")
    mc(u[c][:, 0, 0].abs() ** 2, 1 / k, f"E |(U_{c})_11|^2")
mc(u[1][:, 0, 0].abs() ** 2 * u[2][:, 0, 0].abs() ** 2, 1 / 6, "E |(U_1)_11|^2 |(U_2)_11|^2")   # independent blocks

z = torch.randn(B, dtype=dt)                                         # torch.randn draws standard complex Gaussians:
mc(z.real ** 2, 0.5, "E (Re z)^2")                                   # Re z, Im z ~ N(0, 1/2), uncorrelated
mc(z.imag ** 2, 0.5, "E (Im z)^2")
mc(z.real * z.imag, 0.0, "E Re z Im z")

SU = alg.random_unitary_operator(4, measure='haar_su').matrix
assert torch.allclose(torch.linalg.det(SU[:, 2, :3, :3]), torch.ones(4, dtype=dt))

torch.manual_seed(1)
Z = torch.randn(B, 2, 2, dtype=dt)                                   # Ginibre matrices
torch.manual_seed(1)
V = alg.random_unitary(2, batch_size=B)                              # the library, same random numbers
Q, R = torch.linalg.qr(Z)                                            # LAPACK: diag(R) real, either sign
d = R.diagonal(dim1=-2, dim2=-1)
assert close(V, Q * (d / d.abs()).unsqueeze(-2))                     # Theorem 5.3: V = Q Lambda
raw = Q.diagonal(dim1=-2, dim2=-1).sum(-1).real                      # Re Tr Q without the correction
print(f"E Re Tr Q without the phase correction: {raw.mean():.4f} +- {raw.std() / math.sqrt(B):.4f}")
assert abs(raw.mean()) > 20 * raw.std() / math.sqrt(B)               # far from 0: not Haar
```

*Численное наблюдение.* Без поправки код находит
$`\mathbb E\mathrm{Re}\,\mathrm{Tr}\,Q'\approx-0.84`$ при $`k=2`$, тогда как у унитарной матрицы Хаара $`\mathbb EU=0`$.

<a id="s-5-2"></a>

## 5.2 Моменты и исчисление Вайнгартена

*Конспект:* [раздел 5.2](../../lecture_notes/ru/05_random_matrices.md#s-5-2); [предложение 5.7](../../lecture_notes/ru/05_random_matrices.md#r-5-7) ($`\mathbb E[UXU^*]=E(X)`$), [лемма 5.8](../../lecture_notes/ru/05_random_matrices.md#r-5-8) (коммутант всех $`V\otimes V`$ равен $`\mathrm{span}\{1,F\}`$, где $`F`$ — оператор перестановки), [теорема 5.10](../../lecture_notes/ru/05_random_matrices.md#r-5-10) (формула Вайнгартена второго порядка) и [следствие 5.11](../../lecture_notes/ru/05_random_matrices.md#r-5-11).

Проверим предложение 5.7 сравнением с `channels.center_expectation`, лемму 5.8 — линейной алгеброй (коммутант
трёх случайно выбранных $`V\otimes V`$, вычисленный как нуль-пространство, содержит $`F`$ и двумерен, так что по
лемме 5.8 он уже равен $`\mathrm{span}\{1,F\}`$), а теорему 5.10 и следствие 5.11 — методом Монте-Карло на
блоке $`U(3)`$ из полученных выше выборок:

```python
X = alg.from_blocks([torch.randn(k, k, dtype=dt) for k in alg.k_factors])      # a fixed element
T = (U @ X @ U.adjoint()).matrix                                    # B samples of U X U*
mean, err = T.mean(0), T.std(0) / math.sqrt(B)
EX = channels.center_expectation(alg)(X).matrix[0]
random = err > 1e-12                                                 # entries that are not constant
z = ((mean - EX).abs()[random] / err[random]).max().item()
print(f"twirl over U(M) vs E(X): largest deviation {z:.2f} standard errors")
assert ((mean - EX).abs() <= 5 * err + 1e-12).all()                 # Proposition 5.7

def commutant_dim(ops, tol=1e-9):
    """Dimension of {X : XA = AX for all A in ops}, as the null space of X -> (AX - XA)_A."""
    d = ops[0].shape[0]
    eye = torch.eye(d, dtype=dt)
    ops = [A.contiguous() for A in ops]                              # torch.kron needs contiguous input
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    return int((torch.linalg.svdvals(L) < tol).sum())

k = 3
F = torch.zeros(k * k, k * k, dtype=dt)
for i in range(k):
    for j in range(k):
        F[j * k + i, i * k + j] = 1                                  # F (e_i (x) e_j) = e_j (x) e_i
VV = [torch.kron(W, W) for W in alg.random_unitary(k, batch_size=3).contiguous()]
assert all(close(F @ W, W @ F) for W in VV)
assert commutant_dim(VV) == 2                                        # Lemma 5.8: span{1, F}

u3 = u[2]                                                            # Haar on U(3)
mc(u3[:, 0, 0].abs() ** 4, 2 / (k * (k + 1)), "E|U_11|^4")
mc(u3[:, 0, 0].abs() ** 2 * u3[:, 0, 1].abs() ** 2, 1 / (k * (k + 1)), "E|U_11|^2 |U_12|^2")
mc((u3[:, 0, 0] * u3[:, 1, 1] * (u3[:, 0, 1] * u3[:, 1, 0]).conj()).real,
   -1 / ((k - 1) * k * (k + 1)), "E U_11 U_22 conj(U_12 U_21)")

n_s = 20_000                                                         # Theorem 5.10 for a random X
X2 = torch.randn(k * k, k * k, dtype=dt)
trX, trFX = X2.diagonal().sum(), (F @ X2).diagonal().sum()
alpha, beta = (trX - trFX / k) / (k ** 2 - 1), (trFX - trX / k) / (k ** 2 - 1)
UU = torch.einsum('bij,bkl->bikjl', u3[:n_s], u3[:n_s]).reshape(n_s, k * k, k * k)   # U (x) U
T2 = UU @ X2 @ UU.mH
mean, err = T2.mean(0), T2.std(0) / math.sqrt(n_s)
assert ((mean - alpha * torch.eye(k * k, dtype=dt) - beta * F).abs() <= 5 * err).all()
```

<a id="s-5-3"></a>

## 5.3 Статистика уровней

*Конспект:* [раздел 5.3](../../lecture_notes/ru/05_random_matrices.md#s-5-3): GOE $`H=(G+G^T)/2`$, GUE $`H=(G+G^*)/\sqrt2`$ и отношение соседних интервалов $`r_i=\min(s_i,s_{i+1})/\max(s_i,s_{i+1})`$ для интервалов $`s_i=\lambda_{i+1}-\lambda_i`$; [предложение 5.13](../../lecture_notes/ru/05_random_matrices.md#r-5-13) (гипотеза Вигнера для матриц $`2\times2`$) и [предложение 5.14](../../lecture_notes/ru/05_random_matrices.md#r-5-14) (среднее отношение равно $`2\ln2-1`$ для независимых уровней, $`4-2\sqrt3`$ и $`\frac{2\sqrt3}\pi-\frac12`$ для GOE и GUE размера $`3\times3`$, где $`\mathbb Er=I_1/I_0`$ и $`I_j=\int_0^1x^jf_\beta(x)\,dx`$).

Функция `SpinChain.level_spacing_ratio` (в `chains`) возвращает среднее отношение для каждой строки
уровней. Проверим предложения 5.13 и 5.14 и зафиксируем значение для GOE при $`n=64`$:

```python
ratio = SpinChain.level_spacing_ratio
n2 = 200_000
alg_r = TypeIAlgebra([2], [2], complex_valued=False, precision='double')     # real symmetric 2 x 2
alg_c = TypeIAlgebra([2], [2], precision='double')
G = torch.randn(n2, 1, 2, 2, dtype=torch.float64)
Gc = torch.randn(n2, 1, 2, 2, dtype=dt)
goe2 = alg_r.operator((G + G.mT) / 2, is_self_adjoint=True)
gue2 = alg_c.operator((Gc + Gc.mH) / math.sqrt(2), is_self_adjoint=True)
cdf = {1: lambda t: 1 - math.exp(-math.pi * t ** 2 / 4),
       2: lambda t: math.erf(2 * t / math.sqrt(math.pi)) - 4 * t / math.pi * math.exp(-4 * t ** 2 / math.pi)}
for beta, H, mean_s in ((1, goe2, math.sqrt(math.pi)), (2, gue2, 4 / math.sqrt(math.pi))):
    lam = H.eigenvalues()[0]                                         # (n2, 2), ascending
    S = (lam[:, 1] - lam[:, 0]) / mean_s
    mc(S, 1.0, f"beta={beta}: E S")
    for t in (0.25, 1.0, 2.0):
        mc((S <= t).double(), cdf[beta](t), f"beta={beta}: P(S <= {t})")

levels = torch.rand(4000, 50, dtype=torch.float64)                   # independent levels
mc(ratio(levels), 2 * math.log(2) - 1, "<r>, independent uniform levels")      # Proposition 5.14(a)

n3 = 100_000                                                         # Proposition 5.14(b): 3 x 3
G = torch.randn(n3, 1, 3, 3, dtype=torch.float64)
Gc = torch.randn(n3, 1, 3, 3, dtype=dt)
goe3 = TypeIAlgebra([3], [3], complex_valued=False, precision='double').operator((G + G.mT) / 2, is_self_adjoint=True)
gue3 = TypeIAlgebra([3], [3], precision='double').operator((Gc + Gc.mH) / math.sqrt(2), is_self_adjoint=True)
h = 1e-4
xs = (torch.arange(int(1 / h), dtype=torch.float64) + 0.5) * h      # midpoint rule on [0, 1]
for beta, H, exact in ((1, goe3, 4 - 2 * math.sqrt(3)), (2, gue3, 2 * math.sqrt(3) / math.pi - 0.5)):
    f = (xs + xs ** 2) ** beta / (1 + xs + xs ** 2) ** (1 + 1.5 * beta)
    assert abs((xs * f).sum() / f.sum() - exact) < 1e-7              # E r = I_1 / I_0
    r3 = ratio(H.eigenvalues()[0])                                   # one ratio per matrix
    mc(r3, exact, f"beta={beta}: <r>, 3 x 3")
    mc((r3 <= 0.5).double(), (f * (xs <= 0.5)).sum().item() / f.sum().item(), f"beta={beta}: P(r <= 1/2), 3 x 3")

n = 64                                                               # numerical observation, no theorem
alg64 = TypeIAlgebra([n], [n], complex_valued=False, precision='double')
G = torch.randn(800, 1, n, n, dtype=torch.float64)
lam = alg64.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()[0]
r_goe = ratio(lam[:, n // 4: 3 * n // 4])                            # 800 samples of the mean ratio
print(f"<r>, GOE(64), middle half of the spectrum: {r_goe.mean():.4f} +- {r_goe.std() / math.sqrt(800):.4f}")

def same_mean(a, b, name):
    """Two independent Monte Carlo samples of the same law: means agree within 5 combined standard errors."""
    ea, eb = a.std().item() / math.sqrt(a.numel()), b.std().item() / math.sqrt(b.numel())
    print(f"{name}: {a.mean():.4f} vs {b.mean():.4f} +- {math.hypot(ea, eb):.4f}")
    assert abs(a.mean().item() - b.mean().item()) <= 5 * math.hypot(ea, eb)
```

*Численное наблюдение.* Для средней половины спектра GOE при $`n=64`$ код находит $`\langle r\rangle\approx0.53`$ —
меньше точного значения $`4-2\sqrt3\approx0.536`$ для $`3\times3`$ и намного больше значения
$`2\ln2-1\approx0.386`$ для независимых уровней. Это число — результат вычисления, а не
теорема; ниже оно служит эталоном для блоков GOE того же размера $`64`$, отношения интервалов которых имеют
в точности то же распределение.

<a id="s-5-4"></a>

## 5.4 Почему секторы нельзя смешивать

*Конспект:* [раздел 5.4](../../lecture_notes/ru/05_random_matrices.md#s-5-4); по [теореме 5.15](../../lecture_notes/ru/05_random_matrices.md#r-5-15) уровни разных независимых блоков не отталкиваются: для двух блоков GOE у наименьшего интервала $`\delta`$ объединённого спектра вероятность $`P(\delta<\varepsilon)`$ имеет порядок не меньше $`\varepsilon`$, тогда как внутри одного блока GOE — порядок не больше $`\varepsilon^2`$.

При большем числе секторов эффект сильнее. Предельную теорему для многих секторов мы не доказываем; следующий
эксперимент — численное наблюдение. В библиотеке блоки разделены с самого начала:
`eigenvalues()` возвращает по одному спектру на сектор, а смешать их можно только явным `torch.cat`. Поместим $`C`$
независимых блоков GOE размера $`64`$ в одну алгебру и сравним среднее отношение интервалов по секторам с его значением для
объединённого спектра. Каждый сектор — матрица GOE размера $`64`$, поэтому его отношения имеют в точности то же распределение, что и `r_goe`
выше, и код это проверяет; поведение объединённого спектра лишь наблюдается:

```python
k, merged_r = 64, {}
for C in (1, 2, 4, 8, 16):
    alg_C = TypeIAlgebra([k] * C, [k] * C, complex_valued=False, precision='double')
    G = torch.randn(800 // C, C, k, k, dtype=torch.float64)
    lam = alg_C.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()    # one spectrum per sector
    per_sector = torch.stack([ratio(l[:, k // 4: 3 * k // 4]) for l in lam])
    merged = torch.cat(lam, dim=-1).sort(dim=-1).values                       # the spectrum on C^D
    r_m = ratio(merged[:, k * C // 4: 3 * k * C // 4])
    merged_r[C] = r_m.mean().item()
    print(f"C = {C:2d}: <r> per sector {per_sector.mean():.4f}, merged {merged_r[C]:.4f}"
          f" +- {r_m.std() / math.sqrt(len(r_m)):.4f}")
    same_mean(per_sector.flatten(), r_goe, f"C = {C:2d}: sectors vs GOE(64)")   # every sector is GOE(64)
assert merged_r[1] > merged_r[2] > merged_r[4] > merged_r[16]       # numerical observation
assert abs(merged_r[16] - (2 * math.log(2) - 1)) < 0.01              # numerical observation
```

*Численное наблюдение.* В каждом секторе $`\langle r\rangle\approx0.53`$ — значение GOE(64) из раздела 5.3
(это совпадение — проверка факта: секторы являются матрицами GOE(64)). Для объединённого спектра значение падает до
$`0.42`$ при двух секторах и примерно до $`0.39`$, т. е. до значения $`2\ln2-1`$ для независимых уровней с точностью $`0.01`$
(предложение 5.14(a)), при $`8`$ и $`16`$ секторах. То, что с ростом $`C`$ объединённый спектр приближается к независимым уровням,
код лишь наблюдает; в этом курсе это не доказывается.

То же происходит и в физической модели. Цепочка Гейзенберга со случайными полями (глава 9) сохраняет
намагниченность, поэтому её гамильтониан лежит в $`\bigoplus_{N=0}^LM_{\binom LN}`$ (предложение 9.1). Её блоки *не*
независимы — все они зависят от одних и тех же случайных полей, — поэтому теорема 5.15 в данной формулировке неприменима;
эксперимент ниже показывает тот же эффект в качестве численного наблюдения.

```python
chain = SpinChain(10, complex_valued=False)                         # sectors N = 0, ..., 10
lam = chain.random_field_heisenberg(W=0.5, batch_size=20).eigenvalues()
k5 = chain.sector_dim(5)                                            # 252
r_sector = ratio(lam[5][:, k5 // 4: 3 * k5 // 4]).mean().item()
spectrum = torch.cat(lam, dim=-1).sort(dim=-1).values               # all 1024 levels mixed
r_all = ratio(spectrum[:, 256:768]).mean().item()
print(f"Heisenberg chain, L = 10: <r> in sector N = 5: {r_sector:.4f}, all sectors mixed: {r_all:.4f}")
assert r_sector > 0.49 and r_all < 0.42                             # numerical observation
```

*Численное наблюдение.* Внутри сектора $`N=5`$ код находит $`\langle r\rangle\approx0.53`$, что близко к значению
GOE(64) выше, а спектр со смешанными секторами даёт $`0.39`$, что близко к $`2\ln2-1`$.

<a id="s-5-5"></a>

## 5.5 С чем сравнивать

Унитарные матрицы Хаара: `scipy.stats.unitary_group` и `ortho_group`, случайные по Хаару унитарные матрицы QuTiP,
а также любое QR-разложение `numpy`/`torch` с последующей поправкой фаз.
Исчисление Вайнгартена: символьные пакеты IntU (Mathematica; Z. Puchała, J. A. Miszczak), RTNI
(Mathematica и Python; M. Fukuda, R. König, I. Nechita, J. Phys. A 52 (2019), 425303) и
пакет Python Haarpy символьно вычисляют моменты любого порядка; теорема 5.10 — их случай второго
порядка, а теорема 5.12 (доказанная в приложении B) — общая формула, которую они реализуют. Статистику уровней обычно вычисляют несколькими строками `numpy`,
а программы точной диагонализации, такие как QuSpin, выдают спектры по секторам. Библиотека
добавляет унитарную группу алгебры как единый объект: независимые блоки Хаара всех размеров, порождаемые
пакетными вызовами, твирлинг, отождествлённый с условным ожиданием
`channels.center_expectation`, и спектры, возвращаемые по секторам, так что смешивание секторов никогда не
происходит по умолчанию.

<a id="x-5"></a>

## Упражнения

Численные части [упражнений](../../lecture_notes/ru/05_random_matrices.md#x-5) конспекта лекций:
1. Для $`U`$, распределённой по Хаару на $`U(\mathcal M)`$, $`\mathcal M=M_1\oplus M_2\oplus M_3`$, проверьте, что $`\mathbb E|\mathrm{Tr}\,U|^2=C=3`$ и
   $`\mathbb E|\mathrm{Tr}\,U|^4=17`$.
2. Вычислите численно матрицу Вайнгартена порядка $`3`$ при $`k=3`$ ([теорема 5.12](../../lecture_notes/ru/05_random_matrices.md#r-5-12)) и проверьте, что
   $`\mathbb E|U_{11}|^6=\frac1{10}`$.
3. **Круговой ортогональный ансамбль** (`measure='coe'`) — это $`S=W^TW`$, где $`W`$ распределена по Хаару на $`U(k)`$. Проверьте,
   что COE библиотеки — это $`W^TW`$, что $`S`$ унитарна и симметрична, что $`\mathbb E\mathrm{Tr}\,S^m=0`$ при $`m=1,2,3`$
   и что средняя плотность собственных фаз равномерна. Вычислите $`\langle r\rangle`$ для собственных фаз COE при $`k=32`$
   и сравните со значением для GOE из раздела 5.3 (численное наблюдение).

<details>
<summary>Решения</summary>

*Упражнение 1.*

```python
tr = U.trace                                                        # Tr U = sum_c Tr U_c, shape (B,)
mc(tr.abs() ** 2, 3.0, "E|Tr U|^2")
mc(tr.abs() ** 4, 17.0, "E|Tr U|^4")
```

*Упражнение 2.*

```python
import itertools

def n_cycles(p):                                                    # number of cycles of a permutation
    seen, count = set(), 0
    for i in range(len(p)):
        count += i not in seen
        while i not in seen:
            seen.add(i)
            i = p[i]
    return count

perms = list(itertools.permutations(range(3)))                      # S_3; s.index(j) = s^{-1}(j)
Gram = torch.tensor([[3.0 ** n_cycles([s.index(p[i]) for i in range(3)]) for p in perms] for s in perms],
                    dtype=torch.float64)                            # k^{#cycles(s^{-1} p)}, k = 3
Wg3 = torch.linalg.inv(Gram)                                        # the Weingarten matrix, n = k = 3
assert abs(Wg3.sum().item() - 6 / (3 * 4 * 5)) < 1e-12
mc(u3[:, 0, 0].abs() ** 6, 0.1, "E|U_11|^6")
```

*Упражнение 3.*

   *Численное наблюдение.* Для COE при $`k=32`$ код находит $`\langle r\rangle\approx0.53`$, в пределах
   $`0.015`$ от значения GOE(64) из раздела 5.3. Это совпадение — результат вычисления; никакая
   теорема этого курса не связывает статистику интервалов COE со статистикой GOE (о
   круговых ансамблях см. F. J. Dyson, J. Math. Phys. 3 (1962), 140–156).

```python
alg_s = TypeIAlgebra([32], [32], precision='double')
torch.manual_seed(2)
S = alg_s.random_unitary_operator(400, measure='coe').matrix[:, 0]
torch.manual_seed(2)
W = alg_s.random_unitary(32, batch_size=400)                         # same random numbers
assert close(S, W.mT @ W)                                            # the library's COE is W^T W
assert close(S, S.mT) and close(S @ S.mH, torch.eye(32, dtype=dt))
for m in (1, 2, 3):                                                  # E Tr S^m = 0
    trm = torch.linalg.matrix_power(S, m).diagonal(dim1=-2, dim2=-1).sum(-1)
    mc(trm.real, 0.0, f"E Re Tr S^{m}")
    mc(trm.imag, 0.0, f"E Im Tr S^{m}")
phases = torch.linalg.eigvals(S).angle().sort(dim=-1).values         # in (-pi, pi]
mc(((phases > 0) & (phases <= 1)).sum(-1).double(), 32 / (2 * math.pi), "E #eigenphases in (0, 1]")
s = torch.diff(torch.cat([phases, phases[:, :1] + 2 * math.pi], dim=-1), dim=-1)   # around the circle
r = torch.minimum(s, s.roll(-1, -1)) / torch.maximum(s, s.roll(-1, -1))
r_coe = r.mean(-1)
print(f"<r>, COE(32) eigenphases: {r_coe.mean():.4f}; GOE(64), Section 5.3: {r_goe.mean():.4f}")
assert abs(r_coe.mean() - r_goe.mean()) < 0.015                      # numerical observation
```

</details>

**Далее:** [Руководство 6](06_quantum_channels.md): квантовые каналы в библиотеке.
