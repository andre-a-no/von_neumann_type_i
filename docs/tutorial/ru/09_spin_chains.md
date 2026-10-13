[English](../09_spin_chains.md) · [Español](../es/09_spin_chains.md) · [Français](../fr/09_spin_chains.md) · [Deutsch](../de/09_spin_chains.md) · [中文](../zh/09_spin_chains.md) · [日本語](../ja/09_spin_chains.md) · **Русский**

# 9. Спиновые цепочки

Цепочка из $L$ спинов $\frac12$ живёт в $(\mathbb C^2)^{\otimes L}$, и плотная линейная алгебра в размерности $2^L$
останавливается на $L\approx14$. Симметрии отодвигают этот предел: намагниченность разбивает пространство на
$L+1$ секторов (глава 4), трансляции ещё раз разбивают каждый сектор, а $SU(2)$ ещё раз уменьшает блоки. В этой главе
эти секторы строятся для XXZ-цепочки, доказывается, что они говорят об основных состояниях, библиотека
проверяется по свободным фермионам и анзацу Бете, и мы выходим за пределы плотных матриц (методы Крылова) и
чистых цепочек (беспорядок).

**Вы узнаете**, почему алгебра цепочки есть $\bigoplus_NM_{\binom LN}$ и почему башня этих алгебр имеет
треугольник Паскаля в качестве диаграммы Браттели и индекс $4\cos^2\frac\pi{2L+4}$, как подсчитываются секторы
импульса, что такое правило знаков Маршалла и каковы его следствия для импульса и спина основных состояний
(Либ–Маттис), как решается XX-цепочка преобразованием Йордана–Вигнера, что такое двухмагнонные уравнения Бете,
почему метод Ланцоша сходится монотонно и как измерять статистику уровней и запутанность свободных фермионов.

Всюду далее узлы — это $0,\dots,L-1$, базис $(\mathbb C^2)^{\otimes L}$ состоит из
конфигураций $|b\rangle=|b_0\dots b_{L-1}\rangle$, где $b_i=1$ для спина вверх (соглашение библиотеки: бит
$i$ метки — это узел $i$), $s^a_i$ ($a=x,y,z$) — спиновые матрицы $\frac12\sigma_a$ на узле $i$,
$s^\pm_i=s^x_i\pm is^y_i$, $n_i=s^z_i+\frac12$, $S^a=\sum_is^a_i$, а $N=S^z+\frac L2$ — число спинов,
направленных вверх. Цепочка либо **открытая** (связи $(i,i+1)$, $i<L-1$), либо **периодическая** (ещё и связь
$(L-1,0)$; тогда $L\ge3$, индексы по модулю $L$). $J>0$ (антиферромагнетик), $\hbar=1$, $\log$ — натуральный логарифм.

## 9.1 XXZ-цепочка и её $U(1)$-секторы

**XXZ-цепочка** с анизотропией $\Delta\in\mathbb R$ и полями $h_i\in\mathbb R$ — это
$$
H=J\sum_{(i,j)}\Big[\tfrac12\big(s^+_is^-_j+s^-_is^+_j\big)+\Delta\,s^z_is^z_j\Big]+\sum_ih_is^z_i ,
$$
где сумма берётся по связям; $\frac12(s^+_is^-_j+s^-_is^+_j)=s^x_is^x_j+s^y_is^y_j$ перемещает спин вверх вдоль связи.
При $\Delta=1$ это цепочка Гейзенберга $J\sum\mathbf s_i\cdot\mathbf s_j$, при $\Delta=0$ — XX-цепочка.

**Предложение 9.1.** $[H,S^z]=0$. Следовательно, $H\in U(1)'\cong\bigoplus_{N=0}^LM_{\binom LN}(\mathbb C)$, причём
блок $N$ действует на линейной оболочке конфигураций с $N$ спинами вверх.

*Доказательство.* Диагональные слагаемые коммутируют с диагональным $S^z$. При $i\neq j$ оператор $s^+_is^-_j$ коммутирует с
$s^z_l$, $l\ne i,j$, и $[s^z_i+s^z_j,s^+_is^-_j]=[s^z_i,s^+_i]s^-_j+s^+_i[s^z_j,s^-_j]=s^+_is^-_j-s^+_is^-_j=0$ в силу
$[s^z,s^\pm]=\pm s^\pm$ (раздел 4.5); для $s^-_is^+_j$ переходим к сопряжённым. Остальное — предложение 4.10(2). $\square$

`SpinChain(L).xxz` строит блоки непосредственно по конфигурациям каждого сектора (`chain.labels[N]`).
Сравним их с полной $2^L\times2^L$-матрицей, построенной из кронекеровых произведений:

```python
import functools, itertools, math, time
import numpy as np
import torch
from scipy.optimize import brentq
from torch_vn_algebra import SpinChain, SU2Chain, DensityMatrix, inclusions, krylov

torch.manual_seed(0)
f64, c128 = torch.float64, torch.complex128
close = lambda a, b, tol=1e-10: bool((torch.as_tensor(a) - torch.as_tensor(b)).abs().max() < tol)
def op(L, ops):                             # product of one-site operators {site: 2 x 2}; site 0 = lowest bit
    return functools.reduce(torch.kron, [ops.get(j, torch.eye(2, dtype=f64)) for j in reversed(range(L))])

sp = torch.tensor([[0., 0.], [1., 0.]], dtype=f64)                    # s^+ : |0> = down -> |1> = up
sm, sz = sp.T.contiguous(), torch.diag(torch.tensor([-0.5, 0.5], dtype=f64))
L = 6
chain = SpinChain(L, boundary='periodic')
h = 0.3 * torch.randn(L, dtype=f64)                                  # random fields
H_full = sum(0.5 * (op(L, {i: sp, j: sm}) + op(L, {i: sm, j: sp})) + 0.7 * op(L, {i: sz, j: sz})
             for i, j in chain.bonds) + sum(h[i] * op(L, {i: sz}) for i in range(L))   # J = 1, Delta = 0.7
Sz = sum(op(L, {i: sz}) for i in range(L))
assert close(H_full @ Sz, Sz @ H_full)                               # Proposition 9.1
H = chain.xxz(J=1.0, Delta=0.7, h=h)                                 # an element of (+)_N M_binom(6, N)
for N, k in enumerate(chain.algebra.k_factors):
    lab = chain.labels[N]                                            # configurations of sector N
    assert close(H.matrix[0, N, :k, :k], H_full[lab][:, lab])
```

**Башня цепочек.** Пусть $\mathcal A_L=\bigoplus_NM_{\binom LN}$ действует на $L$ узлах. Добавление узла,
$x\mapsto x\otimes1_2$, отображает $\mathcal A_L$ в $\mathcal A_{L+1}$, поскольку $x\otimes1$ коммутирует с
$S^z\otimes1+1\otimes s^z$. Это даёт башню $\mathcal A_0\subset\mathcal A_1\subset\mathcal A_2\subset\cdots$.

**Предложение 9.2.** Включение $\mathcal A_L\subset\mathcal A_{L+1}$ имеет $(L+2)\times(L+1)$-матрицу
включения $\Lambda_L$ с $\lambda_{N,N}=\lambda_{N,N-1}=1$ и всеми остальными элементами, равными $0$; размеры
блоков подчиняются правилу Паскаля $\binom{L+1}N=\binom LN+\binom L{N-1}$, диаграмма Браттели башни — это
треугольник Паскаля, а двудольный граф матрицы $\Lambda_L$ — путь с $2L+3$ вершинами.

*Доказательство.* Конфигурации с $N$ спинами вверх на $L+1$ узлах — это конфигурации с $N$ спинами вверх на первых
$L$ узлах и новым спином вниз, а также конфигурации с $N-1$ спинами вверх и новым спином вверх. На этих двух подпространствах
$x\otimes1$ действует как $x_N$ и $x_{N-1}$, поэтому блок $N$ алгебры $\mathcal A_{L+1}$ содержит блоки $N$ и $N-1$
алгебры $\mathcal A_L$ по одному разу; сложение размерностей даёт правило Паскаля. В двудольном графе сектор $N$
алгебры $\mathcal A_{L+1}$ соединён с секторами $N-1$ и $N$ алгебры $\mathcal A_L$, что даёт путь
$0'-0-1'-1-\cdots-L'-L-(L+1)'$ (штрихи — для $\mathcal A_{L+1}$). $\square$

Путь от вершины диаграммы Браттели до сектора $N$ на уровне $L$ — это конфигурация $b$; таких путей
$\binom LN$, столько же, каков размер блока. Для связного включения с матрицей включения $\Lambda$ в качестве
**индекса** берём $\|\Lambda\|^2$ — индекс марковского следа (V. Jones, *Invent. Math.* 72 (1983); F. Goodman,
P. de la Harpe, V. Jones, *Coxeter graphs and towers of algebras*, Springer 1989, Ch. 2–3); для
$M_k\otimes1\subset M_k\otimes M_m$ это $m^2$ из раздела 7.2. Его вычисляет `Inclusion.index`.

**Лемма 9.3.** Матрица смежности $A_n$ пути с вершинами $1,\dots,n$ имеет собственные значения
$2\cos\frac{\pi k}{n+1}$, $k=1,\dots,n$, с собственными векторами $v_k(j)=\sin\frac{\pi jk}{n+1}$; следовательно,
$\|A_n\|=2\cos\frac\pi{n+1}$.

*Доказательство.* Положим $v_k(0)=v_k(n+1)=0$. Тогда $(A_nv_k)(j)=v_k(j-1)+v_k(j+1)=2\cos\frac{\pi k}{n+1}\,v_k(j)$ в силу
$\sin(a-b)+\sin(a+b)=2\sin a\cos b$, и $v_k\ne0$, так как $v_k(1)=\sin\frac{\pi k}{n+1}\neq0$. Эти $n$ собственных значений
различны ($\cos$ строго убывает на $(0,\pi)$), поэтому других нет, а норма
симметричной матрицы $A_n$ — наибольший модуль собственного значения, достигаемый при $k=1$ и $k=n$. $\square$

**Следствие 9.4.** $[\mathcal A_{L+1}:\mathcal A_L]=\|\Lambda_L\|^2=4\cos^2\frac\pi{2L+4}$; эта величина возрастает к $4$.

*Доказательство.* Если упорядочить вершины вдоль пути из предложения 9.2, матрица смежности
$B=\begin{pmatrix}0&\Lambda_L\\\Lambda_L^T&0\end{pmatrix}$ становится $A_{2L+3}$ (подобие с помощью перестановки). По лемме 1.1
$\|B\|^2=\|B^2\|$ и $\|\Lambda_L\Lambda_L^T\|=\|\Lambda_L^T\Lambda_L\|=\|\Lambda_L\|^2$, а $B^2=\Lambda_L\Lambda_L^T\oplus\Lambda_L^T\Lambda_L$, поэтому
$\|B\|=\|\Lambda_L\|$; лемма 9.3 при $n=2L+3$ даёт значение. $\square$

Числа $4\cos^2\frac\pi n$, $n\ge3$, — это в точности те значения меньше $4$, которые может принимать индекс подфактора
(Jones, указ. соч.); здесь чётные $n$ появляются по элементарной причине — из-за леммы 9.3.

```python
tower = inclusions.chain_tower(8)                                    # A_0 c A_1 c ... c A_8
print(inclusions.bratteli_diagram(tower).levels[:6])                 # Pascal's triangle
for L_, inc in enumerate(tower):
    sizes = [math.comb(L_ + 1, N) for N in range(L_ + 2)]
    assert inc.M.k_factors == sizes and inc.is_connected
    assert list(inc.matrix @ np.array(inc.N.k_factors)) == sizes      # unitality = Pascal's rule
    n = 2 * L_ + 3
    A_path = np.eye(n, k=1) + np.eye(n, k=-1)
    assert abs(np.linalg.eigvalsh(A_path).max() - 2 * math.cos(math.pi / (n + 1))) < 1e-12   # Lemma 9.3
    assert abs(inc.index() - 4 * math.cos(math.pi / (2 * L_ + 4)) ** 2) < 1e-12           # Corollary 9.4
```

> **Физический смысл.** $S^z$ сохраняется всякий раз, когда взаимодействие инвариантно относительно вращений вокруг
> оси $z$; $N$ считает спины вверх, или частицы после отображения Йордана–Вигнера из раздела 9.4. В башне каждый
> зарядовый сектор более длинной цепочки питается двумя секторами более короткой.

## 9.2 Трансляции, импульс и правило знаков Маршалла

На периодической цепочке пусть $T$ — трансляция $T|b_0b_1\dots b_{L-1}\rangle=|b_{L-1}b_0\dots b_{L-2}\rangle$,
переводящая узел $i$ в узел $i+1$. $T$ — матрица перестановки, поэтому она унитарна, $T^L=1$ и
$Ts^a_iT^*=s^a_{i+1}$. Следовательно, $TS^zT^*=S^z$ и $THT^*=H$ для однородного поля $h_i=h$: $T$ переводит связи в связи.
Её собственные значения — $e^{ik}$ с **импульсами** $k=2\pi m/L$, $m=0,\dots,L-1$.

**Предложение 9.5 (секторы импульса).** Зафиксируем $N$. Для конфигурации $a$ с $N$ спинами вверх пусть $R_a$ —
размер её орбиты $\{T^ja\}$. (a) $R_a$ делит $L$; для $k=2\pi m/L$ с $mR_a\equiv0\pmod L$ вектор
$|a,k\rangle=\sum_{j=0}^{R_a-1}e^{-ikj}T^j|a\rangle$ ненулевой и $T|a,k\rangle=e^{ik}|a,k\rangle$, причём эти
$R_a$ векторов образуют ортогональный базис линейной оболочки орбиты. (b) Собственное значение $e^{2\pi im/L}$ оператора $T$ на
секторе $N$ имеет кратность $d_{N,m}=\#\{\text{орбит }:mR_a\equiv0\pmod L\}$, и всякий оператор, коммутирующий с
$S^z$ и $T$, — например, трансляционно-инвариантная XXZ-цепочка, — блочно-диагонален и лежит в $\bigoplus_N\bigoplus_mM_{d_{N,m}}$.

*Доказательство.* (a) Стабилизатор $a$ в циклической группе $\{T^j\}\cong\mathbb Z_L$ — подгруппа, и $R_a$ — её
индекс, делитель $L$; векторы $T^j|a\rangle$, $0\le j<R_a$, ортонормированы и $T^{R_a}|a\rangle=|a\rangle$. Если
$e^{ikR_a}=1$, то $T|a,k\rangle=e^{ik}\sum_{j=1}^{R_a}e^{-ikj}T^j|a\rangle=e^{ik}|a,k\rangle$, поскольку слагаемое с $j=R_a$ равно
слагаемому с $j=0$; коэффициент при $|a\rangle$ равен $1$. Допустимые $m\in\{0,\dots,L-1\}$ — кратные
$L/R_a$, их $R_a$ штук, с различными собственными значениями унитарного $T$, поэтому векторы ортогональны и порождают
$R_a$-мерную линейную оболочку орбиты. (b) Сектор $N$ — ортогональная сумма линейных оболочек его орбит. Оператор $X$,
коммутирующий с $S^z$ и $T$, коммутирует с их спектральными проекторами $P_N$, $Q_m$ (теорема 1.5(d)), поэтому
$X=\sum_{N,m}P_NQ_mXP_NQ_m$. $\square$

`xxz_momentum` строит блоки в нормированном базисе $|a,k\rangle$, по одному представителю $a$ на орбиту. Мы
независимо подсчитываем орбиты и сравниваем спектры с плотным сектором:

```python
L = 12
ring = SpinChain(L, boundary='periodic')
rot = lambda x: ((x << 1) | (x >> (L - 1))) & ((1 << L) - 1)        # T on labels: site i -> i + 1
for N in range(L + 1):
    orbits = {frozenset(itertools.accumulate(range(L - 1), lambda x, _: rot(x), initial=a))
              for a in ring.labels[N].tolist()}                       # {a, Ta, ..., T^{L-1} a}
    d = [sum((m * len(o)) % L == 0 for o in orbits) for m in range(L)]
    assert d == ring.momentum_dims(N)                                   # Proposition 9.5(b)
print(ring.momentum_dims(6))
Hk = ring.xxz_momentum(J=1.0, Delta=1.0, h=0.2, sector=6)               # 12 blocks instead of one of size 924
E_k = torch.sort(torch.cat([e[0] for e in Hk.eigenvalues()])).values
assert close(E_k, torch.linalg.eigvalsh(ring.xxz(1.0, 1.0, h=0.2, sector=6).matrix[0, 0]))
```

Сектор $N=6$ из $924$ состояний распадается на блоки размеров $80,75,78,76,78,75,80,75,78,76,78,75$.

**Лемма 9.6 (Перрон–Фробениус, симметричный случай).** Пусть $M$ — вещественная симметричная матрица с $M_{ij}\ge0$ при
$i\ne j$, граф которой ($i\sim j$ тогда и только тогда, когда $M_{ij}>0$) связен. Тогда её наибольшее собственное значение $\lambda$ простое,
и у него есть собственный вектор, все компоненты которого $>0$.

*Доказательство.* $\lambda=\max_{\|x\|=1}\langle x,Mx\rangle$, и максимум достигается в точности на единичных векторах
$\lambda$-собственного подпространства (разложите по ортонормированному собственному базису). Для такого $x$ вектор модулей $|x|$
удовлетворяет $\langle|x|,M|x|\rangle\ge\operatorname{Re}\langle x,Mx\rangle=\lambda$, поскольку диагональные слагаемые совпадают, а
$M_{ij}|x_i||x_j|\ge\operatorname{Re}M_{ij}\bar x_ix_j$ при $i\ne j$; значит, $|x|$ — собственный вектор. Если $|x|_i=0$, то
$0=(\lambda-M_{ii})|x|_i=\sum_{j\ne i}M_{ij}|x|_j$ влечёт $|x|_j=0$ для всех соседей $j$ вершины $i$, и в силу
связности $x=0$. Итак, ни у какого ненулевого собственного вектора нет нулевой компоненты; но собственное подпространство размерности $\ge2$
содержит ненулевой вектор с $x_1=0$. Следовательно, $\lambda$ простое, и $|x|>0$ — собственный вектор. $\square$

**Теорема 9.7 (правило знаков Маршалла).** Пусть цепочка открытая или периодическая с чётным $L$; пусть $B$ —
нечётные узлы, а $N_B(b)$ — число спинов вверх конфигурации $b$ на $B$; пусть $J>0$, а $\Delta,h_i\in\mathbb R$ произвольны.
В каждом секторе $N$ наименьшее собственное значение $H_N$ простое, и его собственный вектор можно выбрать так, что
$(-1)^{N_B(b)}\psi(b)>0$ для каждой конфигурации $b$ с $N$ спинами вверх.

*Доказательство.* Пусть $V|b\rangle=(-1)^{N_B(b)}|b\rangle$; $V$ диагонален, $V=V^*=V^{-1}$. Каждая связь соединяет чётный и
нечётный узлы (для кольца — потому что $L$ чётно), поэтому $s^+_is^-_j$ меняет $N_B$ на $\pm1$ и
$Vs^+_is^-_jV=-s^+_is^-_j$. Таким образом, $M=-VH_NV$ имеет внедиагональные элементы $J/2>0$ между конфигурациями, отличающимися
обменом двух спинов некоторой связи, и $0$ в остальных случаях. Её граф связен: обмены
соседей $(i,i+1)$, $i<L-1$, порождают все перестановки узлов, а значит, соединяют любые две конфигурации
с $N$ спинами вверх (обмен одинаковых спинов ничего не меняет). По лемме 9.6 наибольшее собственное значение $M$ — минус
наименьшее собственное значение $H_N$ — простое и имеет собственный вектор $\varphi>0$; возьмём $\psi=V\varphi$. $\square$

**Следствие 9.8 (импульс основного состояния).** На периодической цепочке с чётным $L$, $J>0$ и однородным полем
основное состояние сектора $N$ имеет импульс $k=0$ при чётном $N$ и $k=\pi$ при нечётном $N$.

*Доказательство.* $T$ коммутирует с $H_N$, а основное состояние $\psi=V\varphi$ простое, поэтому $T\psi=e^{ik}\psi$. Далее,
$(T\psi)(b)=\psi(T^{-1}b)$, и у $T^{-1}b$ на узле $i$ стоит спин $b_{i+1}$; поскольку $L$ чётно, $i+1$ чётно при нечётном
$i$, так что $N_B(T^{-1}b)=N-N_B(b)$ и $(VT\psi)(b)=(-1)^{N_B(b)}(-1)^{N-N_B(b)}\varphi(T^{-1}b)=(-1)^N\varphi(T^{-1}b)$. Кроме того,
$VT\psi=e^{ik}\varphi$. И $\varphi$, и $\varphi\circ T^{-1}$ положительны, поэтому $e^{ik}=(-1)^N$. $\square$

```python
for bc in ('open', 'periodic'):
    c = SpinChain(10, bc, complex_valued=False)
    for N in range(11):
        w, V = torch.linalg.eigh(c.xxz(1.0, 0.7, h=0.3 * torch.randn(10, dtype=f64), sector=N).matrix[0, 0])
        x = (-1.0) ** c.bits(N)[:, 1::2].sum(-1) * V[:, 0]                  # (-1)^{N_B(b)} psi(b)
        assert (x * x.sum().sign()).min() > 0 and (len(w) == 1 or w[1] - w[0] > 1e-8)   # Theorem 9.7
for N in range(1, L):
    lowest = [e[0, 0].item() for e in ring.xxz_momentum(1.0, 0.5, sector=N).eigenvalues()]
    m0 = int(np.argmin(lowest))
    assert m0 == (0 if N % 2 == 0 else L // 2)                              # Corollary 9.8: k = 0 or pi
    assert sorted(lowest)[1] - lowest[m0] > 1e-8
```

> **Физический смысл.** В антиферромагнетике решётка фиксирует знаки основного состояния: перемещение
> спина через связь меняет знак. Благодаря этому квантовый метод Монте-Карло свободен от проблемы знака для двудольных
> антиферромагнетиков; для фрустрированных цепочек, где некоторая связь соединяет узлы одной подрешётки, это не так.

## 9.3 $SU(2)$ при $\Delta=1$

При $\Delta=1$ и $h=0$ гамильтониан $H=J\sum\mathbf s_i\cdot\mathbf s_j$ коммутирует с $S^x,S^y,S^z$ (предложение 4.15), поэтому
$H\cong\bigoplus_S1_{2S+1}\otimes H_S$ с блоками $H_S\in M_{m_S}$, $m_S=\binom L{L/2-S}-\binom L{L/2-S-1}$
(предложение 4.13), и модуль `su2` строит $H_S$ на векторах старшего веса (предложение 4.15(b)).
Пусть $E_0(S)$ — наименьшее собственное значение $H_S$.

**Теорема 9.9 (упорядочение Либа–Маттиса для цепочек).** Пусть цепочка такая же, как в теореме 9.7, $\Delta=1$, $h_i=0$,
$J>0$. Для каждого полного спина $S$ системы из $L$ спинов основное состояние сектора $S^z=S$ имеет полный спин $S$. Следовательно,
$E_0(S)$ — простое собственное значение $H_S$, и $E_0(S_{\min})<E_0(S_{\min}+1)<\dots<E_0(\frac L2)$, где
$S_{\min}\in\{0,\frac12\}$: основное состояние — единственный мультиплет наименьшего спина.

*Доказательство.* Пусть $\psi$ — основное состояние сектора $S^z=S$, простое по теореме 9.7. $\mathbf S^2$ коммутирует с
$H$ и сохраняет сектор, поэтому $\mathbf S^2\psi$ лежит в одномерном основном подпространстве: $\psi$ имеет некоторый полный
спин $S'$. Пусть $\chi$ — произведение синглетов $\frac1{\sqrt2}(|{\uparrow\downarrow}\rangle-|{\downarrow\uparrow}\rangle)$ на
парах $(0,1),(2,3),\dots,(L-2S-2,L-2S-1)$, а остальные $2S$ узлов направлены вверх. $S^+$ уничтожает каждый сомножитель,
и $S^z\chi=S\chi$, поэтому $\mathbf S^2\chi=(S^-S^++(S^z)^2+S^z)\chi=S(S+1)\chi$ (раздел 4.5). На каждой паре спин вверх
находится на чётном узле (коэффициент $+$) или на нечётном ($-$), поэтому $(-1)^{N_B(b)}\chi(b)$ имеет один и тот же знак для
всех $b$ из его носителя; с $V,\varphi$, как в теореме 9.7, $\langle\psi,\chi\rangle=\langle\varphi,V\chi\rangle\neq0$. Собственные векторы
$\mathbf S^2$ с различными собственными значениями ортогональны, поэтому $S'=S$. По предложению 4.15(a) и следствию 4.12 каждый
мультиплет спина $S'\ge S$ пересекает сектор $S^z=S$ по одному вектору, так что спектр $H$ там — объединение
спектров $H_{S'}$, $S'\ge S$, с учётом кратностей. Его наименьшее собственное значение простое, а собственный вектор имеет
спин $S$, поэтому это $E_0(S)$, простое в $H_S$ и не являющееся собственным значением $H_{S'}$, $S'>S$: $E_0(S)<E_0(S')$. $\square$

Либ и Маттис доказали больше: для связной двудольной решётки со спинами подрешёток $S_A,S_B$ и константами связи,
антиферромагнитными между подрешётками и ферромагнитными внутри них, $E_0(S)$ возрастает при $S\ge|S_A-S_B|$,
а основное состояние имеет спин $|S_A-S_B|$ (E. Lieb, D. Mattis, *J. Math. Phys.* 3 (1962) 749; приводится без
доказательства). **Синглет-триплетная щель** $E_0(1)-E_0(0)$ положительна по теореме 9.9; на кольце
$L(E_0(1)-E_0(0))\to\pi^2/2$ с логарифмическими поправками (J. des Cloizeaux, J. Pearson, *Phys. Rev.* 128
(1962) 2131; I. Affleck, D. Gepner, H. Schulz, T. Ziman, *J. Phys. A* 22 (1989) 511; приводится без доказательства).
Сравним размеры блоков для двенадцати спинов, проверим теорему 9.9 по $U(1)$-секторам и вычислим щель:

```python
sym = SU2Chain(12, boundary='periodic', complex_valued=False)
H_su2 = sym.heisenberg(J=1.0)                                       # (+)_S M_{m_S}, S = 0, ..., 6
u1_sizes = [math.comb(12, N) for N in range(13)]
print(sym.multiplicities, max(u1_sizes), sum(k ** 3 for k in u1_sizes), sum(m ** 3 for m in sym.multiplicities))
E0_S = sym.lowest_energies(H_su2)[0]                                # E_0(S), S = 0, ..., 6
E0_N = [torch.linalg.eigvalsh(ring.xxz(1.0, 1.0, sector=6 + S).matrix[0, 0])[0] for S in range(7)]
assert close(E0_S, torch.stack(E0_N))                               # Theorem 9.9: sector S^z = S has spin S
assert (E0_S[1:] > E0_S[:-1] + 1e-8).all()                          # E_0(S) strictly increasing
gap = {}
for n in (6, 8, 10, 12):
    s = SU2Chain(n, 'periodic', complex_valued=False)
    e = [torch.linalg.eigvalsh(s.heisenberg(spin=S).matrix[0, 0])[0].item() for S in (0, 1)]
    gap[n] = n * (e[1] - e[0])
print({n: round(g, 4) for n, g in gap.items()}, round(math.pi ** 2 / 2, 4))
assert 0 < gap[6] < gap[8] < gap[10] < gap[12] < math.pi ** 2 / 2
```

Блоки $SU(2)$ имеют размеры $m_S=132,297,275,154,54,11,1$ против $924$ у наибольшего блока $U(1)$, и
стоимость диагонализации $\sum k^3$ падает с $2.0\cdot10^9$ до $5.3\cdot10^7$. Произведение $L$ на щель растёт медленно,
$4.108,4.181,4.232,4.270$ для $L=6,\dots,12$, приближаясь к $\pi^2/2\approx4.935$.

**Корреляции.** Спины узлов образуют тензорный оператор ранга один, и по **теореме Вигнера–Эккарта** их
матричные элементы между мультиплетами равны коэффициенту Клебша–Гордана, умноженному на **приведённый матричный элемент**
$\langle S',a\|\mathbf s_i\|S,b\rangle$, не зависящий от магнитных квантовых чисел (A. R. Edmonds, *Angular
Momentum in Quantum Mechanics*, Princeton 1957, Ch. 5; приводится без доказательства). `spin_correlation` использует её, чтобы
строить $\mathbf s_i\cdot\mathbf s_j$ на каждом $H_S$ только по векторам старшего веса; метод `'direct'`
вместо этого проецирует оператор связи. Сравним оба способа с основным состоянием, выписанным в $U(1)$-секторе:

```python
gs = sym.ground_state(H_su2)
S0, coef = gs.spin.item(), gs.coefficients[0, :sym.multiplicity(gs.spin.item())]
pairs = [(0, r) for r in range(12)]
corr_we = sym.expectation(sym.spin_correlation(pairs), S0, coef)
corr_dir = sym.expectation(sym.spin_correlation(pairs, method='direct'), S0, coef)
psi = sym.embed(S0, coef)                                            # the ground state in the sector N = 6
lab, z = sym.u1.labels[6], sym.u1.bits(6).to(f64) - 0.5
def s_dot_s(i, j):                                                   # <psi| s_i . s_j |psi>, from the definition
    if i == j:
        return 0.75
    flip = ((lab >> i) & 1) != ((lab >> j) & 1)
    partner = sym.u1.index_of[lab[flip] ^ ((1 << i) | (1 << j))]
    return ((psi ** 2 * z[:, i] * z[:, j]).sum() + 0.5 * (psi[flip] * psi[partner]).sum()).item()
assert S0 == 0 and close(corr_we, corr_dir)
assert close(corr_we, torch.tensor([s_dot_s(0, r) for r in range(12)], dtype=f64))
assert ((-1.0) ** torch.arange(1, 12) * corr_we[1:] > 0).all()       # signs (-1)^r
print([round(x, 4) for x in corr_we.tolist()])
```

Корреляции $\langle\mathbf s_0\cdot\mathbf s_r\rangle=0.75,-0.4489,0.1878,-0.1658,0.1209,-0.1224,0.1070,\dots$
чередуются по знаку, как того требуют теоремы 9.7 и 9.9: синглет инвариантен относительно вращений, которые переставляют
оси, и $\langle s^+_0s^+_r\rangle=0$ в секторе $S^z=0$, поэтому $\langle\mathbf s_0\cdot\mathbf s_r\rangle=3\langle s^x_0s^x_r\rangle=
\frac34\langle s^+_0s^-_r+s^-_0s^+_r\rangle=\frac34\sum_{b_0\ne b_r}\psi(b)\psi(b')$, где $b'$ получается из $b$ обменом
спинов узлов $0$ и $r$; так как $N_B(b')-N_B(b)$ равно $\pm1$ при нечётном $r$ и $0$ при чётном $r$, каждое слагаемое имеет знак
$(-1)^r$. Убывание медленное, как $(-1)^r\sqrt{\log r}/r$ при $L\to\infty$ (Affleck et al., указ. соч.; приводится
без доказательства).

## 9.4 Точные проверки: свободные фермионы и анзац Бете

**Йордан–Вигнер.** На открытой цепочке положим $c_j=\big(\prod_{l<j}(-1)^{n_l}\big)s^-_j$.

**Теорема 9.10 (XX-цепочка свободна).** (a) $\{c_i,c_j^*\}=\delta_{ij}$, $\{c_i,c_j\}=0$ и $c_j^*c_j=n_j$.
(b) Открытая XX-цепочка есть $H=\frac J2\sum_{j<L-1}(s^+_js^-_{j+1}+s^-_js^+_{j+1})=\sum_{i,j}h_{ij}c_i^*c_j$ с
$h=\frac J2A_L$. (c) Пусть $u_k$ — ортонормированные вещественные собственные векторы $h$ с собственными значениями
$\varepsilon_k=J\cos\frac{\pi k}{L+1}$ (лемма 9.3) и $b_k=\sum_ju_k(j)c_j$. Для $K=\{k_1<\dots<k_N\}$ векторы
$b_K^*|0\rangle=b_{k_1}^*\cdots b_{k_N}^*|0\rangle$ ($|0\rangle$ — все спины вниз) образуют ортонормированный собственный базис сектора
$N$ с энергиями $\sum_{k\in K}\varepsilon_k$.

*Доказательство.* (a) Операторы на разных узлах коммутируют, а на одном узле $(-1)^n$ антикоммутирует с $s^\pm$. При
$i<j$ операторы $c_i$ и $c_j^*$ (или $c_j$) встречаются на узле $i$ как $s^-_i$ и $(-1)^{n_i}$ — одна антикоммутация — и
в остальном коммутируют, поэтому они антикоммутируют. При $i=j$ квадраты струн равны $1$ и
$c_jc_j^*+c_j^*c_j=s^-_js^+_j+s^+_js^-_j=1$, $c_j^*c_j=s^+_js^-_j=n_j$. (b) $c_j^*c_{j+1}=s^+_j(-1)^{n_j}s^-_{j+1}=s^+_js^-_{j+1}$,
поскольку $s^+_j\neq0$ только на $n_j=0$; далее переходим к сопряжённым. (c) $\{b_k,b_l^*\}=\sum_ju_k(j)u_l(j)=\delta_{kl}$ и
$\{b_k,b_l\}=0$, а $H=\sum_k\varepsilon_kb_k^*b_k$ по спектральному разложению $h$. Из канонических антикоммутационных соотношений (КАС)
$[b_l^*b_l,b_k^*]=\delta_{lk}b_k^*$, поэтому $[H,b_k^*]=\varepsilon_kb_k^*$; вместе с $H|0\rangle=0$ это даёт
$Hb_K^*|0\rangle=(\sum_K\varepsilon_k)b_K^*|0\rangle$. Поскольку $b_k|0\rangle=0$ (каждый $c_j$ опускает спин), перенос операторов уничтожения
вправо с помощью КАС даёт $\langle b_{K'}^*0,b_K^*0\rangle=\delta_{KK'}$ по индукции по $N$. Каждый $c_j^*$ увеличивает
$N$ на единицу, поэтому $b_K^*|0\rangle$ лежит в секторе $|K|$: $\binom LN$ ортонормированных векторов в пространстве такой же размерности. $\square$

```python
L = 10
xx = SpinChain(L)
eps = [math.cos(math.pi * k / (L + 1)) for k in range(1, L + 1)]
for N in range(L + 1):
    E = torch.linalg.eigvalsh(xx.xxz(J=1.0, Delta=0.0, sector=N).matrix[0, 0])
    E_ff = sorted(sum(K) for K in itertools.combinations(eps, N))     # sums over N-element subsets
    assert close(E, torch.tensor(E_ff, dtype=f64))                    # all 1024 levels, Theorem 9.10
    assert abs(E[0].item() - sum(sorted(eps)[:N])) < 1e-12           # ground state: the N lowest modes
```

**Магноны.** На кольце с $h=0$ полностью поляризованное состояние $|0\rangle$ имеет энергию $E_F=J\Delta L/4$. Один спин вверх
в точке $x$ разрывает две связи: $H|x\rangle=(E_F-J\Delta)|x\rangle+\frac J2(|x-1\rangle+|x+1\rangle)$ ($L\ge3$).

**Предложение 9.11 (один магнон).** В секторе $N=1$ кольца вектор $|k\rangle=\sum_xe^{-ikx}|x\rangle$,
$k=2\pi m/L$, удовлетворяет $T|k\rangle=e^{ik}|k\rangle$ и $H|k\rangle=\big(E_F+J(\cos k-\Delta)\big)|k\rangle$.

*Доказательство.* $\sum_xe^{-ikx}(|x-1\rangle+|x+1\rangle)=\sum_y(e^{-ik(y+1)}+e^{-ik(y-1)})|y\rangle=2\cos k\,|k\rangle$ с использованием
$e^{ikL}=1$; та же перенумерация даёт $T|k\rangle=e^{ik}|k\rangle$. $\square$

Два магнона рассеиваются. Обозначим через $|x_1,x_2\rangle$, $x_1<x_2$, состояние со спинами вверх в $x_1,x_2$; пара соседей
разрывает две связи, разнесённая пара — четыре.

**Предложение 9.12 (два магнона; уравнения Бете).** Пусть $k_1,k_2\in\mathbb C$, $K=k_1+k_2$, и
$(A_{12},A_{21})\neq0$ удовлетворяют (i) $A_{12}(1+e^{iK}-2\Delta e^{ik_2})+A_{21}(1+e^{iK}-2\Delta e^{ik_1})=0$ и
(ii) $A_{21}=e^{ik_2L}A_{12}$, $A_{12}=e^{ik_1L}A_{21}$. Пусть
$\psi(x_1,x_2)=A_{12}e^{i(k_1x_1+k_2x_2)}+A_{21}e^{i(k_2x_1+k_1x_2)}$ на $\mathbb Z^2$ и
$\Psi=\sum_{0\le x_1<x_2<L}\psi(x_1,x_2)|x_1,x_2\rangle$. Если $\Psi\neq0$, то
$H\Psi=\big(E_F+J(\cos k_1+\cos k_2-2\Delta)\big)\Psi$. Если знаменатели не обращаются в нуль, (i) и (ii) означают
$$
e^{ik_1L}=-\frac{1+e^{iK}-2\Delta e^{ik_1}}{1+e^{iK}-2\Delta e^{ik_2}},\qquad
e^{ik_2L}=-\frac{1+e^{iK}-2\Delta e^{ik_2}}{1+e^{iK}-2\Delta e^{ik_1}} .
$$

*Доказательство.* Пусть $E=E_F+J(\cos k_1+\cos k_2-2\Delta)$. (1) Для всех $(x_1,x_2)\in\mathbb Z^2$ у каждой плоской волны
$\psi(x_1\pm1,x_2)$ и $\psi(x_1,x_2\pm1)$ в сумме дают $2(\cos k_1+\cos k_2)\psi(x_1,x_2)$. (2) Подстановкой
$\psi(x,x)+\psi(x+1,x+1)-2\Delta\psi(x,x+1)$ равно $e^{iKx}$, умноженному на левую часть (i), поэтому обращается в нуль при всех $x$.
(3) По (ii) $\psi(x_2,x_1+L)=A_{12}e^{ik_2L}e^{i(k_2x_1+k_1x_2)}+A_{21}e^{ik_1L}e^{i(k_1x_1+k_2x_2)}=\psi(x_1,x_2)$, и,
применяя это дважды, $\psi(x_1+L,x_2+L)=\psi(x_1,x_2)$. Представители конфигурации $\{y_1<y_2\}$ в
$D=\{x_1<x_2<x_1+L\}$ — это $(y_1+nL,y_2+nL)$ и $(y_2+nL,y_1+(n+1)L)$, поэтому $\psi$ — функция конфигураций.
Пусть $(x_1,x_2)\in D$. Если $2\le x_2-x_1\le L-2$, все четыре перескока остаются в $D$, а диагональный элемент равен $E_F-2J\Delta$, поэтому
$(H\Psi)(x_1,x_2)=(E_F-2J\Delta)\psi+\frac J2\cdot2(\cos k_1+\cos k_2)\psi=E\psi$ по (1). Если $x_2=x_1+1$, перескоки
$(x_1+1,x_2)$, $(x_1,x_2-1)$ заблокированы, а диагональный элемент равен $E_F-J\Delta$; по (1) и (2)
$(H\Psi)(x_1,x_2)=(E_F-J\Delta)\psi+\frac J2\big(2(\cos k_1+\cos k_2)\psi-\psi(x_1,x_1)-\psi(x_1+1,x_1+1)\big)=E\psi$.
Если $x_2=x_1+L-1$, используем представителя $(x_2,x_1+L)$, относящегося к предыдущему случаю ($L\ge3$). $\square$

Предложение 9.12 говорит лишь о том, что решения дают собственные векторы; оно не утверждает, что все собственные векторы
получаются таким образом. Запишем $k_{1,2}=K/2\pm q$ и $c=\cos\frac K2$. Тогда $1+e^{iK}-2\Delta e^{ik_{1,2}}=2e^{iK/2}(c-\Delta e^{\pm iq})$,
поэтому $A_{12}=c-\Delta e^{iq}$, $A_{21}=-(c-\Delta e^{-iq})$ решают (i). Перемножая два уравнения (ii), получаем
$e^{iKL}=1$, $K=2\pi n/L$, и тогда (ii) равносильно $F(q)=(-1)^ne^{iqL}(c-\Delta e^{-iq})+c-\Delta e^{iq}=0$ —
одному уравнению с одной неизвестной. Мы ищем его корни
$q\in(0,\pi)$, где $e^{-iqL/2}F(q)$ вещественно или чисто мнимо, поскольку $c-\Delta e^{\pm iq}$ комплексно сопряжены
(состояния рассеяния), и на прямых $q=i\kappa$, $\pi+i\kappa$, $\kappa>0$, где $F$ вещественна (связанные состояния);
для $L=9$, $\Delta=0.6$ строим $\Psi$ для каждого корня и сравниваем с сектором $N=2$ из $36$ состояний:

```python
L, Delta = 9, 0.6
ring9 = SpinChain(L, boundary='periodic')
H2 = ring9.xxz(J=1.0, Delta=Delta, sector=2).matrix[0, 0]
x1, x2 = np.triu_indices(L, 1)                                       # all pairs x1 < x2
rows = ring9.index_of[torch.tensor((1 << x1) | (1 << x2))]
def roots(g, grid):                                                  # sign changes of g on the grid, refined
    v = np.array([g(t) for t in grid])
    return [brentq(g, grid[i], grid[i + 1], xtol=1e-15) for i in np.nonzero(v[:-1] * v[1:] < 0)[0]]
E_bethe, vecs, n_bound = [], [], 0
for n in range(L):
    K, c = 2 * np.pi * n / L, math.cos(np.pi * n / L)
    F = lambda q: (-1) ** n * np.exp(1j * q * L) * (c - Delta * np.exp(-1j * q)) + c - Delta * np.exp(1j * q)
    G = lambda q: (lambda z: z.real + z.imag)(np.exp(-0.5j * q * L) * F(q))    # real, or imaginary
    bound = [s + 1j * t for s in (0, np.pi)                          # q = i kappa and pi + i kappa
             for t in roots(lambda t: F(s + 1j * t).real, np.linspace(1e-3, 20, 4000))]
    n_bound += len(bound)
    for q in roots(G, (np.arange(1, 2000) + 0.5 ** 0.5) * np.pi / 2000) + bound:
        psi = np.exp(0.5j * K * (x1 + x2)) * ((c - Delta * np.exp(1j * q)) * np.exp(1j * q * (x1 - x2))
                                              - (c - Delta * np.exp(-1j * q)) * np.exp(-1j * q * (x1 - x2)))
        v = torch.zeros(36, dtype=c128)
        v[rows] = torch.tensor(psi / np.linalg.norm(psi))
        E = 9 * Delta / 4 + 2 * (c * np.cos(q) - Delta)                 # E_F + J (cos k1 + cos k2 - 2 Delta)
        assert abs(E.imag) < 1e-12 and close(H2 @ v, E.real * v, 1e-9)   # Proposition 9.12
        E_bethe.append(E.real)
        vecs.append(v)
assert torch.linalg.matrix_rank(torch.stack(vecs)) == 36             # a complete eigenbasis
assert close(torch.sort(torch.tensor(E_bethe)).values, torch.linalg.eigvalsh(H2))
print(len(E_bethe), n_bound)
ring8, Phi = SpinChain(8, 'periodic'), torch.zeros(28, dtype=c128)     # even L: Phi = sum_x (-1)^x |x, x+1>
for x in range(8):
    Phi[ring8.index_of[(1 << x) | (1 << (x + 1) % 8)]] = (-1) ** x
assert close(ring8.xxz(1.0, Delta, sector=2).matrix[0, 0] @ Phi, (8 * Delta / 4 - Delta) * Phi)   # E_F - J Delta
```

При $L=9$ все $36$ собственных состояний получаются из уравнений Бете: $34$ состояния рассеяния и $2$ связанных состояния.
Полнота здесь проверена, а не доказана, и при чётном $L$ она нарушается. Если $\Delta\ne0$, то $\Phi=\sum_x(-1)^x|x,x+1\rangle$
удовлетворяет $T\Phi=-\Phi$ и $H\Phi=(E_F-J\Delta)\Phi$ (перескоки в $|x-1,x+1\rangle$ и $|x,x+2\rangle$ попарно сокращаются).
Вектор Бете удовлетворяет $T\Psi=e^{-iK}\Psi$, так как $\psi(x_1-1,x_2-1)=e^{-iK}\psi(x_1,x_2)$; если $e^{iK}=-1$, то
$\cos k_1+\cos k_2=2\cos\frac K2\cos\frac{k_1-k_2}2=0$ и его энергия равна $E_F-2J\Delta$. Поэтому $\Phi$ ортогонален всем
векторам Бете (последние строки кода проверяют $H\Phi$ для $L=8$). При $\Delta=0.6$ и $4\le L\le13$ описанный выше
поиск находит все состояния при нечётном $L$ и все состояния, кроме $\Phi$, при чётном $L$. Для $N>2$ магнонов анзац содержит $N!$ слагаемых и $N$ связанных уравнений Бете (H. Bethe, *Z. Phys.* 71 (1931)
205; M. Karbach, G. Müller, *Computers in Physics* 11 (1997) 36).

> **Физический смысл.** Магнон — это спиновая волна: один перевёрнутый спин, делокализованный по кольцу с
> импульсом $k$. Два магнона рассеиваются со сдвигом фазы (вещественное $q$) или образуют связанное состояние (комплексное $q$). Идея Бете
> состояла в том, что для $N$ магнонов рассеяние факторизуется в двухчастичные фазы; его уравнения дают, например,
> энергию основного состояния на узел $\frac14-\log2$ бесконечной цепочки Гейзенберга (L. Hulthén, *Ark. Mat. Astron.
> Fys.* 26A (1938) no. 11; приводится без доказательства).

## 9.5 Методы Ланцоша и Крылова

Сектор $N=L/2$ имеет размерность $12870$ при $L=16$ и $184756$ при $L=20$, но $H$ разрежен ($O(L)$ ненулевых
элементов в строке, `xxz_sparse`), а методам Крылова нужны только произведения $Hv$. Пусть $H$ — эрмитова матрица размера $n$ и
$v\ne0$. **Крыловское подпространство** — это $\mathcal K_m=\operatorname{span}\{v,Hv,\dots,H^{m-1}v\}$.

**Предложение 9.13 (Ланцош).** Пусть $\dim\mathcal K_m=m$. Положим $q_1=v/\|v\|$ и при $j<m$ пусть $q_{j+1}$ —
нормированная компонента $Hq_j$, ортогональная $\mathcal K_j$. Тогда $q_1,\dots,q_m$ — ортонормированный базис
$\mathcal K_m$; при $Q_m=(q_1\cdots q_m)$ матрица $T_m=Q_m^*HQ_m$ вещественная симметричная трёхдиагональная с диагональю
$\alpha_j=\langle q_j,Hq_j\rangle$ и внедиагональными элементами $\beta_j=\langle q_{j+1},Hq_j\rangle>0$; и $HQ_m=Q_mT_m+r_me_m^T$, где
$r_m\perp\mathcal K_m$ обращается в нуль тогда и только тогда, когда $\mathcal K_{m+1}=\mathcal K_m$ (иначе $r_m=\beta_mq_{m+1}$, $q_{m+1}=r_m/\|r_m\|$).

*Доказательство.* Пусть $\operatorname{span}\{q_1,\dots,q_j\}=\mathcal K_j$ (верно при $j=1$). Тогда $H\mathcal K_j\subseteq\mathcal K_j+\mathbb CHq_j$,
так как $Hq_i\in\mathcal K_{i+1}$; поэтому $\mathcal K_{j+1}=\mathcal K_j+\mathbb CH^jv\subseteq\mathcal K_j+\mathbb CHq_j\subseteq\mathcal K_{j+1}$ и
$\mathcal K_{j+1}=\mathcal K_j+\mathbb CHq_j$. При $j<m$ размерность растёт, поэтому компонента $Hq_j$, ортогональная
$\mathcal K_j$, ненулевая, и индукция продолжается; при $j=m$ эта компонента, $r_m$, обращается в нуль тогда и только тогда, когда
$\mathcal K_{m+1}=\mathcal K_m$. При $i>j+1$ имеем $q_i\perp\mathcal K_{j+1}\ni Hq_j$, поэтому $(T_m)_{ij}=0$, и $(T_m)_{ji}=0$, так как $T_m$
эрмитова. $\alpha_j$ вещественно, а $\beta_j=\langle q_{j+1},Hq_j\rangle$ — норма компоненты, $>0$. Таким образом,
$Hq_j=\beta_{j-1}q_{j-1}+\alpha_jq_j+\beta_jq_{j+1}$ при $j<m$ и $Hq_m=\beta_{m-1}q_{m-1}+\alpha_mq_m+r_m$. $\square$

Собственные значения $\theta^{(m)}_1\le\dots\le\theta^{(m)}_m$ матрицы $T_m$ называются **значениями Ритца**.

**Теорема 9.14 (монотонная сходимость).** Пусть $\lambda_1\le\dots\le\lambda_n$ — собственные значения $H$ и
$\dim\mathcal K_m=m$. Тогда $\theta^{(m)}_j\ge\lambda_j$, и $\theta^{(m)}_j\ge\theta^{(m+1)}_j$, если $\dim\mathcal K_{m+1}=m+1$
($j\le m$). Если у $v$ есть ненулевая компонента в собственном подпространстве $\lambda_1$ и $\mathcal K_{m+1}=\mathcal K_m$, то
$\theta^{(m)}_1=\lambda_1$.

*Доказательство.* По теореме Куранта–Фишера (стандартная линейная алгебра) $\theta^{(m)}_j=\min_{U}\max_{x\in U,\|x\|=1}\langle x,Hx\rangle$
по $j$-мерным $U\subseteq\mathcal K_m$, поскольку $Q_m$ — изометрия на $\mathcal K_m$ и
$\langle y,T_my\rangle=\langle Q_my,HQ_my\rangle$. Та же формула по $U\subseteq\mathcal K_{m+1}$ даёт $\theta^{(m+1)}_j$, а
по $U\subseteq\mathbb C^n$ — $\lambda_j$; минимизация по большему множеству подпространств даёт меньшее значение. Если
$\mathcal K_{m+1}=\mathcal K_m$, то $\mathcal K_m$ инвариантно относительно $H$ и содержит $p(H)v$ для каждого многочлена $p$,
в частности $P_1v\ne0$ ($P_1$ — спектральный проектор $\lambda_1$, многочлен от $H$ по интерполяции
Лагранжа, как в теореме 1.5); поэтому $\theta^{(m)}_1\le\langle x,Hx\rangle=\lambda_1$ при $x=P_1v/\|P_1v\|\in\mathcal K_m$. $\square$

Скорость сходимости $\theta^{(m)}_1\to\lambda_1$ зависит от щели $\lambda_2-\lambda_1$ относительно ширины спектра
(оценки Каниэля–Пейджа–Саада через многочлены Чебышёва; Y. Saad, *Numerical Methods for Large Eigenvalue
Problems*, SIAM 2011, Ch. 6; приводится без доказательства). В арифметике с плавающей точкой $q_j$ теряют ортогональность; `krylov.lanczos`
выполняет полную реортогонализацию. Для XX-цепочки точная энергия основного состояния известна (теорема 9.10):

```python
L, N = 12, 6
xx = SpinChain(L, complex_valued=False)
Hs = xx.xxz_sparse(J=1.0, Delta=0.0, sector=N)                      # sparse: diagonal + hopping
alpha, beta, Q = krylov.lanczos(Hs.matvec, torch.randn(1, xx.sector_dim(N), dtype=f64), 40)
T = torch.diag(alpha[0]) + torch.diag(beta[0], 1) + torch.diag(beta[0], -1)
assert close(Q[0] @ Q[0].T, torch.eye(40, dtype=f64)) and close(Q[0] @ Hs.to_dense()[0] @ Q[0].T, T)
ritz = torch.tensor([torch.linalg.eigvalsh(T[:m, :m])[0].item() for m in range(1, 41)])   # T_m = leading block
E_xx = sum(sorted(math.cos(math.pi * k / (L + 1)) for k in range(1, L + 1))[:N])
assert (ritz[1:] <= ritz[:-1] + 1e-12).all() and (ritz >= E_xx - 1e-12).all()             # Theorem 9.14
print([f"{(ritz[m - 1] - E_xx).item():.1e}" for m in (5, 10, 20, 40)])
```

Ошибка наименьшего значения Ритца монотонно убывает: $0.42$, $3.1\cdot10^{-3}$ и $3.1\cdot10^{-6}$ после
$5$, $10$ и $20$ шагов, а после $40$ шагов достигает уровня ошибок округления, при размерности $924$. `krylov.ground_state`
перезапускает метод Ланцоша, пока невязка $\|H\psi-E\psi\|$ не станет меньше $10^{-10}$. Проверим его по энергиям кольца
Гейзенберга, приведённым с девятью знаками в `scripts/chains.py` (эталонные значения сопутствующей статьи), — плотным методом
до $L=12$ и методом Ланцоша при $L=16$:

```python
RING = {4: -2.0, 6: -2.802775638, 8: -3.651093409, 10: -4.515446354, 12: -5.387390917,
        14: -6.263549534, 16: -7.142296361}
for n in (4, 6, 8, 10, 12):
    e0 = torch.linalg.eigvalsh(SpinChain(n, 'periodic').heisenberg(sector=n // 2).matrix[0, 0])[0].item()
    assert abs(e0 - RING[n]) < 5e-10
    print(n, round(e0 / n, 6))
t0 = time.time()
E16, psi16, res = krylov.ground_state(SpinChain(16, 'periodic', complex_valued=False).xxz_sparse(sector=8),
                                      generator=torch.Generator().manual_seed(0))
print(f"L = 16: E0 = {E16.item():.10f}, residual {res.item():.1e}, {time.time() - t0:.2f} s; "
      f"infinite chain: {0.25 - math.log(2):.6f}")
assert abs(E16.item() - RING[16]) < 5e-10 and res.item() < 1e-9
```

При $L=16$ (размерность $12870$) метод Ланцоша работает заметно меньше секунды на одном потоке CPU. Энергия на узел,
$-0.5,-0.467129,-0.456387,-0.451545,-0.448949$ для $L=4,\dots,12$, приближается снизу к значению $\frac14-\log2=-0.443147$
для бесконечной цепочки.

**Эволюция во времени.** Для $e^{-i\tau H}v$ проецируем на $\mathcal K_m$: $w_\tau=\|v\|\,Q_me^{-i\tau T_m}e_1$.

**Предложение 9.15.** Пусть $\dim\mathcal K_m=m$. Для каждого $\tau$ имеем $\|w_\tau\|=\|v\|$ и $\langle w_\tau,Hw_\tau\rangle=\langle v,Hv\rangle$. Если
$\mathcal K_m$ инвариантно относительно $H$, то $w_\tau=e^{-i\tau H}v$.

*Доказательство.* $Q_m$ — изометрия, а $e^{-i\tau T_m}$ унитарна ($T_m$ вещественная симметричная), что даёт утверждение о норме.
Поскольку $Q_m^*HQ_m=T_m$ коммутирует с $e^{-i\tau T_m}$,
$\langle w_\tau,Hw_\tau\rangle=\|v\|^2\langle e^{-i\tau T_m}e_1,T_me^{-i\tau T_m}e_1\rangle=\|v\|^2\langle e_1,T_me_1\rangle=\langle v,Hv\rangle$. Если
$\mathcal K_m$ инвариантно, то $HQ_m=Q_mT_m$ (предложение 9.13, $r_m=0$), поэтому $H^jQ_m=Q_mT_m^j$ и $e^{-i\tau H}Q_m=Q_me^{-i\tau T_m}$; применяем к
$e_1$, учитывая $Q_me_1=v/\|v\|$. $\square$

Итак, каждый шаг Крылова сохраняет норму и энергию, точен он или нет; для точности нужно, чтобы $m$ было велико по сравнению с
$\tau\|H\|$, и `krylov.evolve` разбивает интервал времени на части и оценивает ошибку апостериори
(M. Hochbruck, C. Lubich, *SIAM J. Numer. Anal.* 34 (1997) 1911; приводится без доказательства). Выполним квенч из неелевского
состояния $|{\uparrow\downarrow\uparrow\downarrow}\cdots\rangle$ на кольце Гейзенберга и проследим за шахматной намагниченностью
$m_s=\frac2L\sum_i(-1)^i\langle s^z_i\rangle$, сравнивая с плотным пропагатором при $L=12$:

```python
times, m_s = [0.0, 0.5, 1.0, 2.0, 4.0], {}
stag = lambda c, N, out: [((p[0].abs() ** 2) @ (c.bits(N).to(f64) - 0.5) @ torch.tensor(
    [(-1.0) ** i for i in range(c.L)], dtype=f64) * 2 / c.L).item() for p in out]
for n in (12, 16):
    c = SpinChain(n, 'periodic')
    Hn = c.xxz_sparse(J=1.0, Delta=1.0, sector=n // 2)
    psi0 = c.vector_in_sector('10' * (n // 2))[None]
    out = krylov.evolve(Hn, psi0, times)                               # (times, batch, dim)
    energy = torch.stack([(p.conj() * Hn.matvec(p)).sum(-1).real for p in out])
    assert close(out.norm(dim=-1), 1.0) and close(energy, energy[0])   # Proposition 9.15
    m_s[n] = stag(c, n // 2, out)
    print(n, [round(x, 4) for x in m_s[n]])
    if n == 12:
        w, U = torch.linalg.eigh(c.xxz(1.0, 1.0, sector=6).matrix[0, 0])
        for t, p in zip(times, out):
            assert close(p[0], U @ (torch.exp(-1j * t * w) * (U.mH @ psi0[0])))
print(f"{max(abs(a - b) for a, b in zip(m_s[12][:4], m_s[16][:4])):.1e}")   # L = 12 vs 16 up to t = 2
```

$m_s$ убывает от $1$ через $0.770$ и $0.279$ до $-0.183$ при $t=2$, для обоих размеров с точностью до $1.7\cdot10^{-5}$; при
$t=4$ значения различаются ($0.026$ против $0.035$): становится заметен конечный размер меньшего кольца.

## 9.6 Беспорядок: статистика уровней и запутанность

**Цепочка Гейзенберга в случайном поле** имеет независимые поля $h_i$, равномерно распределённые на $[-W,W]$. Она сохраняет $S^z$
(предложение 9.1) и в общем случае больше ничего: поля нарушают трансляции, отражения, переворот спинов и
$SU(2)$. Поэтому статистика уровней вычисляется отдельно в каждом секторе $N$ (теорема 5.15 и раздел 5.4: смешивание секторов имитирует
независимые уровни). При слабом беспорядке уровни отталкиваются, как в GOE, $\langle r\rangle\approx0.5307$ (раздел 5.3);
при сильном беспорядке собственные состояния локализуются (**многочастичная локализация**), а уровни становятся независимыми,
$\langle r\rangle=2\log2-1\approx0.3863$ (предложение 5.14). Вычислим $\langle r\rangle$ в средней половине сектора
$N=L/2$ с одной стандартной ошибкой по независимым реализациям (каждая даёт одно среднее отношение):

```python
gen = torch.Generator().manual_seed(1)
r_mean = {}
for n, R in ((8, 200), (10, 100)):
    c = SpinChain(n, 'periodic', complex_valued=False)
    k = c.sector_dim(n // 2)
    for W in (0.5, 1, 2, 4, 8, 16):
        E = torch.linalg.eigvalsh(c.random_field_heisenberg(W, R, sector=n // 2, generator=gen).matrix[:, 0])
        r = SpinChain.level_spacing_ratio(E[:, k // 4: 3 * k // 4])     # one value per realisation
        r_mean[n, W] = (r.mean().item(), (r.std() / math.sqrt(R)).item())
    print(n, '  '.join(f"W={W}: {m:.3f}+-{s:.3f}" for (n_, W), (m, s) in r_mean.items() if n_ == n))
assert abs(r_mean[10, 0.5][0] - 0.5307) < 0.015 and abs(r_mean[10, 16][0] - 0.3863) < 0.015
assert r_mean[10, 1][0] > r_mean[10, 2][0] > r_mean[10, 4][0] > r_mean[10, 8][0]
```

При $L=10$ среднее отношение меняется от $0.530$ ($W=0.5$) через $0.484$ ($W=2$) и $0.410$ ($W=4$) до $0.384$
($W=16$) со стандартными ошибками $0.003$–$0.004$. При $L=8$ значение при $W=0.5$ ($0.507$) меньше значения при
$W=1$ ($0.525$): при $W=0$ цепочка интегрируема и обладает всеми перечисленными выше симметриями, а для значения GOE нужен
беспорядок, большой по сравнению с расстоянием между уровнями. Кроссовер широк и смещается с ростом $L$; конечноразмерный
скейлинг помещает его около $W\approx3.7$ (D. Luitz, N. Laflorencie, F. Alet, *Phys. Rev. B* 91 (2015) 081103), а
вопрос о том, выживает ли локализованная фаза при $L\to\infty$, остаётся дискуссионным (J. Šuntajs, J. Bonča, T. Prosen, L. Vidmar,
*Phys. Rev. E* 102 (2020) 062144).

**Запутанность свободных фермионов.** Для собственных состояний открытой XX-цепочки энтропия запутанности
первых $n$ узлов определяется матрицей $n\times n$.

**Теорема 9.16 (формула Пешеля).** Пусть $\Phi=b_{k_1}^*\cdots b_{k_N}^*|0\rangle$, как в теореме 9.10,
$C_{ij}=\langle\Phi,c_i^*c_j\Phi\rangle$ и $C_A=(C_{ij})_{i,j<n}$. Тогда $0\le C_A\le1$, а энтропия запутанности
узлов $A=\{0,\dots,n-1\}$, где $\nu$ пробегает собственные значения $C_A$ с учётом кратности, равна
$$
S_A=\sum_\nu \eta(\nu),\qquad \eta(\nu)=-\nu\log\nu-(1-\nu)\log(1-\nu)\quad(0\log0=0).
$$

*Доказательство.* Для вещественного вектора $x$ положим $b(x)=\sum_jx_jc_j$; тогда $b(x)^*b(y)^*=-b(y)^*b(x)^*$ и
$\{b(x),b(y)^*\}=\langle x,y\rangle$. Пусть $U\subseteq\mathbb R^L$ — линейная оболочка векторов $u_{k_a}$, а $P_A$, $P_B=1-P_A$ — координатные проекторы на
$\{0,\dots,n-1\}$ и его дополнение. Диагонализация формы $\langle P_Ax,P_Ay\rangle$ на $U$ даёт ортонормированный базис
$\phi_1,\dots,\phi_N$ пространства $U$ с $\langle P_A\phi_a,P_A\phi_b\rangle=\nu_a\delta_{ab}$, а значит,
$\langle P_B\phi_a,P_B\phi_b\rangle=(1-\nu_a)\delta_{ab}$. Поскольку отображение $(x_1,\dots,x_N)\mapsto b(x_1)^*\cdots b(x_N)^*|0\rangle$
полилинейно и кососимметрично, $\Phi=\pm\prod_a\big(\sqrt{\nu_a}f_a^*+\sqrt{1-\nu_a}\,g_a^*\big)|0\rangle$ с
$f_a=b(P_A\phi_a/\sqrt{\nu_a})$, $g_a=b(P_B\phi_a/\sqrt{1-\nu_a})$ (слагаемые с нулевым коэффициентом опускаются).
Раскрывая скобки, получаем $\Phi=\sum_I\pm\prod_{a\in I}\sqrt{\nu_a}\prod_{a\notin I}\sqrt{1-\nu_a}\,F_I^*G_{I^c}^*|0\rangle$ по
подмножествам $I$, где $F_I^*$, $G_{I^c}^*$ — упорядоченные произведения. При $j\ge n$ струна оператора $c_j^*$ содержит всё
$A$: $c_j^*=(-1)^{N_A}\otimes\tilde c_j^*$ с операторами Йордана–Вигнера $\tilde c_j$ узлов $\ge n$, тогда как при
$i<n$ имеем $c_i^*=\hat c_i^*\otimes1$. Так как $(-1)^{N_A}$ равен $1$ на $|0_A\rangle$,
$F_I^*G_{I^c}^*|0\rangle=\xi_I\otimes\chi_{I^c}$ с $\xi_I=\hat F_I^*|0_A\rangle$, $\chi_{I^c}=\tilde G_{I^c}^*|0_B\rangle$, и эти семейства
ортонормированы в силу вычисления с КАС из теоремы 9.10. Это разложение Шмидта (лемма 7.7) с
весами $p_I=\prod_{a\in I}\nu_a\prod_{a\notin I}(1-\nu_a)$ — распределением независимых бернуллиевских величин, поэтому
$S_A=-\sum_Ip_I\log p_I=\sum_a\eta(\nu_a)$. Наконец, $C_{ij}=\sum_a\phi_a(i)\phi_a(j)$
(разложите $c_i^*c_j$ по $b(\phi_a)$), поэтому $C_A=\sum_a(P_A\phi_a)(P_A\phi_a)^T$ имеет собственные значения $\nu_a$ на
ортогональных векторах $P_A\phi_a$ и $0$ на остальном, а $\eta(0)=0$. $\square$

`entanglement_entropy` вычисляет редуцированное состояние сектор за сектором (оно коммутирует с намагниченностью
блока, предложение 7.8) и его энтропию. Основное состояние сектора $N$ — это $\Phi$ с $N$ нижними
модами; оно простое, поскольку $\varepsilon_k$ различны:

```python
L, N = 12, 6
xx = SpinChain(L)
w, V = torch.linalg.eigh(xx.xxz(J=1.0, Delta=0.0, sector=N).matrix[0, 0])
rho = DensityMatrix(xx.sector_algebra(N), matrix=torch.outer(V[:, 0], V[:, 0].conj())[None, None])
eps_k, u = np.linalg.eigh(np.eye(L, k=1) / 2 + np.eye(L, k=-1) / 2)  # h = A_L / 2, ascending
C = u[:, :N] @ u[:, :N].T                                            # correlation matrix of the N lowest modes
bern = lambda nu: -(nu * np.log(nu) + (1 - nu) * np.log(1 - nu))
S_peschel = [bern(np.clip(np.linalg.eigvalsh(C[:n, :n]), 1e-15, 1 - 1e-15)).sum() for n in range(1, L)]
S_lib = [xx.entanglement_entropy(rho, n, sector=N).item() for n in range(1, L)]
assert np.allclose(S_lib, S_peschel, atol=1e-9)                       # Theorem 9.16
print([round(s, 4) for s in S_lib[:6]])
```

Энтропии $0.6931,0.5199,0.7456,0.6283,0.7719,0.6560$ для $n=1,\dots,6$ осциллируют в зависимости от чётности $n$
и растут как $\frac16\log L$ при фиксированном $n/L$ — закон конформной теории поля для открытой цепочки (I. Peschel,
*J. Phys. A* 36 (2003) L205; P. Calabrese, J. Cardy, *J. Stat. Mech.* (2004) P06002; приводится без доказательства).

## 9.7 С чем сравнивать

QuSpin (Python) строит гамильтонианы спиновых цепочек в секторах намагниченности, импульса, чётности и
инверсии спинов и использует разреженные решатели задачи на собственные значения и матричные экспоненты из SciPy; специализированные программы точной диагонализации
(XDiag, SPINPACK, HΦ) оптимизированы для больших разреженных задач, и с решёточными симметриями и на больших машинах
такие программы достигли примерно 40–50 спинов; ITensor, TeNPy и TensorKit.jl (с деревьями слияния $SU(2)$) обрабатывают
длинные цепочки с помощью матричных произведений состояний. То, что добавляет библиотека, скромно и конкретно: цепочка как
`TypeIAlgebra` $\bigoplus_NM_{\binom LN}$, так что состояния, каналы, энтропии и динамика глав 3–8
применимы сектор за сектором; башня включений с её индексом; блоки $SU(2)$ из векторов старшего веса с
корреляциями по Вигнеру–Эккарту; пакеты реализаций беспорядка в одном вызове; и редуцированные состояния, никогда не покидающие
секторов. В ней нет секторов отражения и переворота спинов, секторы импульса есть только для однородных XXZ-цепочек, плотные
методы работают до $L\approx16$, методы Крылова — до $L\approx24$ (около $26$ на большом GPU), тензорных сетей нет.

## Упражнения

1. (Секторы импульса вручную.) Для $L=4$ и $N=2$ найдите орбиты шести конфигураций под действием $T$
   и размерности $d_{2,m}$, $m=0,\dots,3$.
2. (Правила сумм.) Пусть $\psi$ — основное состояние кольца Гейзенберга с чётным $L$. Покажите, что
   $E_0=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle$ и $\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle=0$.
3. (Слишком много симметрии.) Покажите, что в секторе $N=L/2$ вектор, все коэффициенты которого равны $1$, кратен
   $(S^-)^{L/2}|{\uparrow\cdots\uparrow}\rangle$ и является собственным вектором кольца Гейзенберга с энергией
   $JL/4$, а также что каждое крыловское подпространство $\mathcal K_m(H,v)$ лежит в любом $H$-инвариантном подпространстве, содержащем $v$.
   Что возвращает `krylov.ground_state`, если начать с этого вектора?

<details>
<summary>Решения</summary>

1. $0011\to1001\to1100\to0110\to0011$ — одна орбита с $R=4$, а $0101\to1010\to0101$ — одна с $R=2$. По
   предложению 9.5 первая даёт вклад в каждое $m$, вторая — в $m$ с $2m\equiv0\pmod4$, т. е. $m=0,2$:
   $d_{2,m}=2,1,2,1$.

```python
assert SpinChain(4, 'periodic').momentum_dims(2) == [2, 1, 2, 1]
```

2. Основное состояние — простой синглет (теорема 9.9), и $T$ коммутирует с $H$, поэтому $T\psi=e^{ik}\psi$. Так как
   $\mathbf s_i\cdot\mathbf s_{i+1}=T^i(\mathbf s_0\cdot\mathbf s_1)T^{*i}$, величина $\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle$ одна и та же для всех $L$ связей, и
   $E_0=J\sum_i\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle$. Так как $\mathbf S^2\psi=0$,
   $0=\langle\mathbf S^2\rangle=\sum_{i,j}\langle\mathbf s_i\cdot\mathbf s_j\rangle=L\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle$, снова в силу трансляционной
   инвариантности.

```python
assert abs(12 * corr_we[1].item() - gs.energy.item()) < 1e-10 and abs(corr_we.sum().item()) < 1e-10
```

3. $S^-$ опускает один из спинов, направленных вверх; применяя его $L/2$ раз к полностью поляризованному состоянию, мы приходим к каждой
   конфигурации с $L/2$ спинами вверх $(L/2)!$ способами, все с коэффициентом $1$. У $|{\uparrow\cdots\uparrow}\rangle$ все
   $L$ связей параллельны, энергия равна $JL/4$, а $S^-$ коммутирует с $H$, поэтому вектор — собственный с той же
   энергией. Если $W$ инвариантно и $v\in W$, то $H^jv\in W$ для всех $j$. Здесь $\mathcal K_m=\mathbb Cv$ при всех $m$, и
   `ground_state` возвращает $E=JL/4$ — верхний край спектра (каждое $\mathbf s_i\cdot\mathbf s_j\le\frac14$); случайный начальный
   вектор избегает этого с вероятностью единица.

```python
c10 = SpinChain(10, 'periodic', complex_valued=False)
E_sym, _, _ = krylov.ground_state(c10.xxz_sparse(1.0, 1.0, sector=5), v0=torch.ones(1, 252, dtype=f64))
assert abs(E_sym.item() - 2.5) < 1e-10 and E_sym.item() > RING[10] + 7     # far from E_0
```

</details>

**Далее:** в [главе 10](10_quantum_computing.md) квантовые вычисления рассматриваются с помощью секторов: подпространства,
свободные от декогеренции, и бесшумные подсистемы, коррекция ошибок для операторных алгебр и миксеры, сохраняющие
вес Хэмминга, — алгебра $\bigoplus_wM_{\binom nw}$ из этой главы в новой роли.
