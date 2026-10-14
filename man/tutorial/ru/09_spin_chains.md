[English](../09_spin_chains.md) · [Español](../es/09_spin_chains.md) · [Français](../fr/09_spin_chains.md) · [Deutsch](../de/09_spin_chains.md) · [中文](../zh/09_spin_chains.md) · [日本語](../ja/09_spin_chains.md) · **Русский**

# Руководство 9. Спиновые цепочки

Это руководство сопровождает [главу 9](../../lecture_notes/ru/09_spin_chains.md) конспекта лекций, в которой строятся секторы XXZ-цепочки —
намагниченности ([предложение 9.1](../../lecture_notes/ru/09_spin_chains.md#r-9-1)), импульса ([предложение 9.5](../../lecture_notes/ru/09_spin_chains.md#r-9-5)) и полного спина — и доказывается, что
они говорят об основных состояниях ([теорема 9.7](../../lecture_notes/ru/09_spin_chains.md#r-9-7), [теорема 9.9](../../lecture_notes/ru/09_spin_chains.md#r-9-9)). Мы строим цепочки с помощью `SpinChain` и
`SU2Chain`, проверяем библиотеку на свободных фермионах и анзаце Бете, используем методы Ланцоша и Крылова из
`krylov` для секторов, слишком больших для плотных матриц, и вычисляем статистику уровней и запутанность.

Всюду далее узлы — это $0,\dots,L-1$, базис $(\mathbb C^2)^{\otimes L}$ состоит из
конфигураций $|b\rangle=|b_0\dots b_{L-1}\rangle$, где $b_i=1$ для спина вверх (соглашение библиотеки: бит
$i$ метки — это узел $i$), $s^a_i$ ($a=x,y,z$) — спиновые матрицы $\frac12\sigma_a$ на узле $i$,
$s^\pm_i=s^x_i\pm is^y_i$, $n_i=s^z_i+\frac12$, $S^a=\sum_is^a_i$, а $N=S^z+\frac L2$ — число спинов,
направленных вверх. Цепочка либо **открытая** (связи $(i,i+1)$, $i<L-1$), либо **периодическая** (ещё и связь
$(L-1,0)$; тогда $L\ge3$, индексы по модулю $L$). $J>0$ (антиферромагнетик), $\hbar=1$, $\log$ — натуральный логарифм.

<a id="s-9-1"></a>

## 9.1 XXZ-цепочка и её $U(1)$-секторы

*Конспект:* [раздел 9.1](../../lecture_notes/ru/09_spin_chains.md#s-9-1); [предложение 9.1](../../lecture_notes/ru/09_spin_chains.md#r-9-1) ($H$ сохраняет $S^z$), [предложение 9.2](../../lecture_notes/ru/09_spin_chains.md#r-9-2), [лемма 9.3](../../lecture_notes/ru/09_spin_chains.md#r-9-3) и [следствие 9.4](../../lecture_notes/ru/09_spin_chains.md#r-9-4) (башня цепочек и её индекс).

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

Башня цепочек $\mathcal A_0\subset\mathcal A_1\subset\cdots$: `inclusions.chain_tower` строит включения,
`inclusions.bratteli_diagram` — их диаграмму Браттели, а `Inclusion.index` — индекс $\|\Lambda\|^2$. Код проверяет
правило Паскаля, [лемму 9.3](../../lecture_notes/ru/09_spin_chains.md#r-9-3) и $[\mathcal A_{L+1}:\mathcal A_L]=4\cos^2\frac\pi{2L+4}$ ([следствие 9.4](../../lecture_notes/ru/09_spin_chains.md#r-9-4)).

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

<a id="s-9-2"></a>

## 9.2 Трансляции, импульс и правило знаков Маршалла

*Конспект:* [раздел 9.2](../../lecture_notes/ru/09_spin_chains.md#s-9-2); [предложение 9.5](../../lecture_notes/ru/09_spin_chains.md#r-9-5) (секторы импульса), [теорема 9.7](../../lecture_notes/ru/09_spin_chains.md#r-9-7) (правило знаков Маршалла), [следствие 9.8](../../lecture_notes/ru/09_spin_chains.md#r-9-8) (импульс основного состояния).

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

По правилу знаков Маршалла ([теорема 9.7](../../lecture_notes/ru/09_spin_chains.md#r-9-7)) для открытых и периодических цепочек со случайными полями
$(-1)^{N_B(b)}\psi(b)$ имеет один знак в каждом секторе, а по [следствию 9.8](../../lecture_notes/ru/09_spin_chains.md#r-9-8) нижний уровень сектора $N$
кольца лежит в блоке импульса $k=0$ при чётном $N$ и $k=\pi$ при нечётном $N$:

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

<a id="s-9-3"></a>

## 9.3 $SU(2)$ при $\Delta=1$

*Конспект:* [раздел 9.3](../../lecture_notes/ru/09_spin_chains.md#s-9-3); [теорема 9.9](../../lecture_notes/ru/09_spin_chains.md#r-9-9) (Либ–Маттис). При $\Delta=1$ и $h=0$ имеем $H\cong\bigoplus_S1_{2S+1}\otimes H_S$, и `SU2Chain` строит блоки $H_S$ на векторах старшего веса.

Сравним размеры блоков для двенадцати спинов, проверим теорему 9.9 по $U(1)$-секторам, вычислим щель и проверим
теоремы 9.7 и 9.9 на двудольном графе, который не является цепочкой ($n_A=5$, $n_B=3$, $S_*=1$, случайные константы связи):

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
print({n: round(g, 4) for n, g in gap.items()})
assert 0 < gap[6] < gap[8] < gap[10] < gap[12]                      # > 0: Theorem 9.9; increase: numerical
n8, A_set = 8, (0, 1, 2, 3, 4)                                      # A = {0, ..., 4}, B = {5, 6, 7}
edges = [(0, 5), (1, 5), (1, 6), (2, 6), (3, 6), (4, 7), (0, 7), (0, 1), (2, 3), (5, 6)]   # connected
g8 = torch.Generator().manual_seed(2)
J8 = [(1 if (i in A_set) != (j in A_set) else -1) * (0.5 + torch.rand(1, generator=g8).item()) for i, j in edges]
dot = lambda i, j: 0.5 * (op(n8, {i: sp, j: sm}) + op(n8, {i: sm, j: sp})) + op(n8, {i: sz, j: sz})
H8 = sum(Jij * dot(i, j) for Jij, (i, j) in zip(J8, edges))          # J_ij > 0 between, < 0 within A, B
Sp8, Sz8 = sum(op(n8, {i: sp}) for i in range(n8)), sum(op(n8, {i: sz}) for i in range(n8))
S2 = Sp8.T @ Sp8 + Sz8 @ Sz8 + Sz8                                  # S^2 = S^- S^+ + (S^z)^2 + S^z
n_up = torch.tensor([bin(b).count('1') for b in range(2 ** n8)])
n_B = torch.tensor([bin(b >> 5).count('1') for b in range(2 ** n8)])  # up spins on B
for N in range(4, 9):                                               # sectors S^z = M = N - 4 >= 0
    idx = torch.nonzero(n_up == N).flatten()
    w, U = torch.linalg.eigh(H8[idx][:, idx])
    x, spin = (-1.0) ** n_B[idx] * U[:, 0], max(N - 4, 1)            # expected spin max(M, S_*)
    assert (x * x.sum().sign()).min() > 0 and (len(w) == 1 or w[1] - w[0] > 1e-8)          # Theorem 9.7
    assert abs(U[:, 0] @ S2[idx][:, idx] @ U[:, 0] - spin * (spin + 1)) < 1e-10           # Theorem 9.9(a)
w8 = torch.linalg.eigvalsh(H8)
assert w8[2] - w8[0] < 1e-10 and w8[3] - w8[0] > 1e-8                # ground space: one triplet, S_* = 1
```

Блоки $SU(2)$ имеют размеры $m_S=132,297,275,154,54,11,1$ против $924$ у наибольшего блока $U(1)$, и
стоимость диагонализации $\sum k^3$ падает с $2.0\cdot10^9$ до $5.3\cdot10^7$. *Численное наблюдение:* произведение $L$ на щель
растёт медленно, $4.108,4.181,4.232,4.270$ для $L=6,\dots,12$; о его поведении при $L\to\infty$ мы ничего не доказываем.

**Корреляции.** Спины узлов образуют векторные операторы, а их сферические компоненты — тензорные операторы ранга
один (лемма B.15(a),(b)). По **теореме Вигнера–Эккарта** (теорема B.16) их матричные элементы между двумя мультиплетами
равны коэффициенту Клебша–Гордана, умноженному на **приведённый матричный элемент**, не зависящий от магнитных
квантовых чисел, и $\mathbf s_i\cdot\mathbf s_j=\sum_q(-1)^qs^{(i)}_qs^{(j)}_{-q}$ (лемма B.15(a)). `spin_correlation` использует это, чтобы
строить $\mathbf s_i\cdot\mathbf s_j$ на каждом $H_S$ только по векторам старшего веса; метод `'direct'` вместо этого сжимает
оператор связи на векторы старшего веса, что корректно по предложению 4.15(b), так как $\mathbf s_i\cdot\mathbf s_j\in\{S_a\}'$.
Корректность реализации первого метода здесь не доказывается; код проверяет, что оба метода согласуются друг с другом и
с определением, вычисленным на основном состоянии, выписанном в $U(1)$-секторе:

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
чередуются по знаку, как того требуют теоремы 9.7 и 9.9 для любого чётного кольца ([конспект лекций, раздел 9.3](../../lecture_notes/ru/09_spin_chains.md#s-9-3)).

<a id="s-9-4"></a>

## 9.4 Точные проверки: свободные фермионы и анзац Бете

*Конспект:* [раздел 9.4](../../lecture_notes/ru/09_spin_chains.md#s-9-4); по [теореме 9.10](../../lecture_notes/ru/09_spin_chains.md#r-9-10) энергии сектора $N$ открытой XX-цепочки — это суммы $N$ различных $\varepsilon_k=J\cos\frac{\pi k}{L+1}$; [предложение 9.12](../../lecture_notes/ru/09_spin_chains.md#r-9-12) (два магнона, уравнения Бете).

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

По [предложению 9.12](../../lecture_notes/ru/09_spin_chains.md#r-9-12) каждый корень $q$ функции $F(q)=(-1)^ne^{iqL}(c-\Delta e^{-iq})+c-\Delta e^{iq}$, где
$K=2\pi n/L$ и $c=\cos\frac K2$, даёт собственный вектор $\Psi$ сектора $N=2$ кольца
([конспект лекций, раздел 9.4](../../lecture_notes/ru/09_spin_chains.md#s-9-4)). Код ищет корни $q\in(0,\pi)$ (состояния рассеяния) и на прямых
$q=i\kappa$, $\pi+i\kappa$, $\kappa>0$ (связанные состояния), строит $\Psi$ для каждого корня при $L=9$, $\Delta=0.6$ и сравнивает
с сектором $N=2$ из $36$ состояний:

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

*Численное наблюдение:* при $L=9$ и $\Delta=0.6$ поиск находит $36$ корней, $34$ состояния рассеяния и $2$ связанных
состояния, векторы которых — линейно независимые собственные векторы (с точностью до округления): все $36$ собственных
состояний получаются из уравнений Бете. Полнота для других $L$ не доказана, а при чётном $L$ и $\Delta\ne0$ она
нарушается ([конспект лекций, раздел 9.4](../../lecture_notes/ru/09_spin_chains.md#s-9-4)): вектор $\Phi=\sum_x(-1)^x|x,x+1\rangle$ ортогонален всем векторам Бете,
и последние строки кода проверяют $H\Phi=(E_F-J\Delta)\Phi$ при $L=8$.

<a id="s-9-5"></a>

## 9.5 Методы Ланцоша и Крылова

*Конспект:* [раздел 9.5](../../lecture_notes/ru/09_spin_chains.md#s-9-5); [предложение 9.13](../../lecture_notes/ru/09_spin_chains.md#r-9-13) (Ланцош: $T_m=Q_m^*HQ_m$ трёхдиагональна), [теорема 9.14](../../lecture_notes/ru/09_spin_chains.md#r-9-14) (значения Ритца, собственные значения $T_m$, удовлетворяют $\theta^{(m)}_j\ge\lambda_j$ и убывают с ростом $m$), [предложение 9.15](../../lecture_notes/ru/09_spin_chains.md#r-9-15) (эволюция во времени методом Крылова).

Сектор $N=L/2$ имеет размерность $12870$ при $L=16$ и $184756$ при $L=20$, но $H$ разрежен ($O(L)$ ненулевых
элементов в строке, `xxz_sparse`), а методам Крылова нужны только произведения $Hv$. Пусть $H$ — эрмитова матрица размера $n$ и
$v\ne0$. **Крыловское подпространство** — это $\mathcal K_m=\mathrm{span}\{v,Hv,\dots,H^{m-1}v\}$.

Теорема 9.14 ничего не говорит о том, как быстро $\theta^{(m)}_1$ приближается к $\lambda_1$; мы лишь наблюдаем это ниже.
В арифметике с плавающей точкой $q_j$ теряют ортогональность; `krylov.lanczos` выполняет полную реортогонализацию, и код
проверяет, что вычисленные $q_j$ ортонормированы и $Q_m^*HQ_m=T_m$ с точностью до округления. Для открытой XX-цепочки
точная энергия основного состояния сектора $N$ — сумма $N$ наименьших $\varepsilon_k$ (теорема 9.10(c)):

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

*Численное наблюдение:* ошибка наименьшего значения Ритца монотонно убывает (как и должно быть по теореме 9.14):
$0.42$, $3.1\cdot10^{-3}$ и $3.1\cdot10^{-6}$ после $5$, $10$ и $20$ шагов, а после $40$ шагов достигает уровня ошибок
округления, при размерности $924$. `krylov.ground_state` перезапускает метод Ланцоша (со случайного гауссовского вектора,
если не задан `v0`), пока невязка $\|H\psi-E\psi\|$ его нормированной пары Ритца $(E,\psi)$ не станет меньше
$10^{-10}\max(1,|E|)$ (не более $20$ перезапусков, иначе выдаётся предупреждение), и возвращает невязку. Малая невязка
гарантирует близкое собственное значение, но не то, что оно наименьшее ([конспект лекций, раздел 9.5](../../lecture_notes/ru/09_spin_chains.md#s-9-5)), поэтому мы сравниваем с плотной диагонализацией, которая
служит нашим единственным эталоном (вычисляемым здесь же, так что значения численные): сектор $N=L/2$ до $L=12$, а при
$L=16$ — блок импульса $k=0$ размерности $810$, который по следствию 9.8 содержит основное состояние сектора $N=8$. При
$L=4$ значение $-2J$ точное ([конспект лекций, раздел 9.5](../../lecture_notes/ru/09_spin_chains.md#s-9-5)).

```python
E_dense = {}                                                         # numerical references, dense
for n in (4, 6, 8, 10, 12):
    E_dense[n] = torch.linalg.eigvalsh(SpinChain(n, 'periodic').heisenberg(sector=n // 2).matrix[0, 0])[0].item()
    print(n, round(E_dense[n] / n, 6))
assert abs(E_dense[4] + 2) < 1e-12                                   # exact value for L = 4
ring16 = SpinChain(16, 'periodic')
k0 = ring16.momentum_dims(8)[0]                                      # dimension of the block k = 0
H16 = ring16.xxz_momentum(J=1.0, Delta=1.0, sector=8).matrix[0, 0, :k0, :k0]
E_dense[16] = torch.linalg.eigvalsh(H16)[0].item()                   # ground state of N = 8 (Corollary 9.8)
t0 = time.time()
E16, psi16, res = krylov.ground_state(SpinChain(16, 'periodic', complex_valued=False).xxz_sparse(sector=8),
                                      generator=torch.Generator().manual_seed(0))
print(f"L = 16: E0 = {E16.item():.10f}, residual {res.item():.1e}, {time.time() - t0:.2f} s; k = 0 block: {k0}")
assert abs(E16.item() - E_dense[16]) < 1e-9 and res.item() < 1e-9
```

При $L=16$ (размерность $12870$) метод Ланцоша работает заметно меньше секунды на одном потоке CPU. *Численное
наблюдение:* энергия на узел, $-0.5,-0.467129,-0.456387,-0.451545,-0.448949$ для $L=4,\dots,12$, возрастает с ростом $L$;
о её пределе мы ничего не доказываем.

**Эволюция во времени.** Для $e^{-i\tau H}v$ проецируем на $\mathcal K_m$: $w_\tau=\|v\|\,Q_me^{-i\tau T_m}e_1$.

По [предложению 9.15](../../lecture_notes/ru/09_spin_chains.md#r-9-15) $\|w_\tau\|=\|v\|$ и $\langle w_\tau,Hw_\tau\rangle=\langle v,Hv\rangle$ для любого $\tau$, и
$w_\tau=e^{-i\tau H}v$, если $\mathcal K_m$ $H$-инвариантно.

Итак, каждый шаг Крылова сохраняет норму и энергию, точен он или нет; о точности предложение 9.15 ничего не говорит,
если только $\mathcal K_m$ не инвариантно. `krylov.evolve` разбивает интервал времени на короткие шаги, каждый из которых
выполняется в крыловском подпространстве с $m=30$, и выдаёт предупреждение, когда эвристический индикатор ошибки велик;
никакой оценки ошибки для него мы не доказываем, а вместо этого сравниваем его с $e^{-itH}$, вычисленным по плотному
спектральному разложению (теорема 1.5). Выполним квенч из неелевского
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

*Численное наблюдение:* $m_s$ убывает от $1$ через $0.770$ и $0.279$ до $-0.183$ при $t=2$, одинаково для обоих размеров
с точностью до $1.7\cdot10^{-5}$; при $t=4$ значения различаются ($0.026$ против $0.035$).

<a id="s-9-6"></a>

## 9.6 Беспорядок: статистика уровней и запутанность

*Конспект:* [раздел 9.6](../../lecture_notes/ru/09_spin_chains.md#s-9-6); [теорема 9.16](../../lecture_notes/ru/09_spin_chains.md#r-9-16) (формула Пешеля).

**Цепочка Гейзенберга в случайном поле** имеет независимые поля $h_i$, равномерно распределённые на $[-W,W]$. Она сохраняет $S^z$
(предложение 9.1), поэтому статистика уровней вычисляется отдельно в каждом секторе $N$ (теорема 5.15 и раздел 5.4:
смешивание секторов имитирует независимые уровни). Если два поля различны, то $THT^*-H=\sum_i(h_{i-1}-h_i)s^z_i\ne0$, так что
трансляции больше не являются симметрией. Вычислим $\langle r\rangle$ в средней половине сектора $N=L/2$ с одной стандартной
ошибкой по независимым реализациям (каждая даёт одно среднее отношение) и сравним с двумя эталонами: $2\log2-1\approx0.3863$,
точным средним отношением для независимых экспоненциальных расстояний (предложение 5.14(a)), и средним отношением для
матриц GOE (в смысле определения из раздела 5.3; код выбирает $G+G^T$, удвоенную матрицу GOE, что не меняет ни одного
отношения) того же размера $252$, что и сектор при $L=10$, выбираемых здесь же. Эти сравнения — численные наблюдения, а
не теоремы: о спектре случайной цепочки мы ничего не доказываем.

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
X = torch.randn(100, k, k, dtype=f64, generator=gen)                  # 100 GOE matrices of size k = 252
r_goe = SpinChain.level_spacing_ratio(torch.linalg.eigvalsh(X + X.mT)[:, k // 4: 3 * k // 4])
print(f"GOE({k}): {r_goe.mean().item():.3f}+-{(r_goe.std() / 10).item():.3f}, Poisson: {2 * math.log(2) - 1:.4f}")
assert abs(r_mean[10, 0.5][0] - r_goe.mean().item()) < 0.015               # numerical observation
assert abs(r_mean[10, 16][0] - (2 * math.log(2) - 1)) < 0.015               # numerical observation
assert r_mean[10, 1][0] > r_mean[10, 2][0] > r_mean[10, 4][0] > r_mean[10, 8][0]
```

*Численное наблюдение:* при $L=10$ среднее отношение меняется от $0.530$ ($W=0.5$), близкого к выборочному значению GOE
$0.533$, через $0.484$ ($W=2$) и $0.410$ ($W=4$) до $0.384$ ($W=16$), близкого к $2\log2-1$, со стандартными ошибками
$0.003$–$0.004$. При $L=8$ значение при $W=0.5$ ($0.507$) меньше значения при $W=1$ ($0.525$). При $W=0$ кольцо
коммутирует с $T$ (раздел 9.2) и с полным спином (раздел 9.3), поэтому его сектор $N=L/2$ распадается на блоки; поля
связывают эти блоки только через $\sum_ih_is^z_i$, операторная норма которого не больше $\frac{WL}2$.

> **Физический смысл** (физическая интерпретация, а не теорема этого курса). Физики интерпретируют изменение
> $\langle r\rangle$ как кроссовер от уровней, отталкивающихся как уровни GOE (слабый беспорядок), к уровням, ведущим себя как
> независимые, и приписывают его **многочастичной локализации** собственных состояний (сильный беспорядок). Кроссовер
> смещается с ростом $L$; где он находится при $L\to\infty$ и существует ли локализованная фаза вообще — открытые
> вопросы физической литературы (D. Luitz, N. Laflorencie, F. Alet, *Phys. Rev. B* 91 (2015) 081103;
> J. Šuntajs, J. Bonča, T. Prosen, L. Vidmar, *Phys. Rev. E* 102 (2020) 062144). Провал при $L=8$, $W=0.5$ объясняют
> остатком дополнительной блочной структуры при $W=0$ (ср. раздел 5.4).

**Запутанность свободных фермионов.** Для собственных состояний открытой XX-цепочки энтропия запутанности
первых $n$ узлов определяется матрицей $n\times n$.

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

*Численное наблюдение:* энтропии $0.6931,0.5199,0.7456,0.6283,0.7719,0.6560$ для $n=1,\dots,6$ осциллируют в
зависимости от чётности $n$. Теорема 9.16 принадлежит I. Peschel (*J. Phys. A* 36 (2003) L205).

<a id="s-9-7"></a>

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

<a id="x-9"></a>

## Упражнения

[Упражнения 1–3](../../lecture_notes/ru/09_spin_chains.md#x-9) конспекта лекций: размерности $d_{2,m}=2,1,2,1$ секторов импульса при
$L=4$, $N=2$; правила сумм $E_0=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle$ и $\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle=0$ для
основного состояния кольца Гейзенберга; и `krylov.ground_state`, запущенный с вектора, все коэффициенты которого равны $1$,
в секторе $N=L/2$: он возвращает $E=JL/4$ — верхний край спектра. Проверьте всё это численно.

<details>
<summary>Решения</summary>

```python
assert SpinChain(4, 'periodic').momentum_dims(2) == [2, 1, 2, 1]
```

```python
assert abs(12 * corr_we[1].item() - gs.energy.item()) < 1e-10 and abs(corr_we.sum().item()) < 1e-10
```

```python
c10 = SpinChain(10, 'periodic', complex_valued=False)
E_sym, _, _ = krylov.ground_state(c10.xxz_sparse(1.0, 1.0, sector=5), v0=torch.ones(1, 252, dtype=f64))
assert abs(E_sym.item() - 2.5) < 1e-10 and E_sym.item() > E_dense[10] + 7  # far from E_0
```

</details>

**Далее:** [Руководство 10. Квантовые вычисления](10_quantum_computing.md).
