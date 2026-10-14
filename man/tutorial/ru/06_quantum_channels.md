[English](../06_quantum_channels.md) · [Español](../es/06_quantum_channels.md) · [Français](../fr/06_quantum_channels.md) · [Deutsch](../de/06_quantum_channels.md) · [中文](../zh/06_quantum_channels.md) · [日本語](../ja/06_quantum_channels.md) · **Русский**

# Руководство 6. Квантовые каналы

Это руководство сопровождает [главу 6](../../lecture_notes/ru/06_quantum_channels.md) конспекта лекций, в которой доказываются описания вполне
положительных отображений между прямыми суммами матричных алгебр — матрицы Чоя, операторы Крауса, дилатация
Стайнспринга — и изучаются сопряжённые отображения, неподвижные точки, матрицы переходов и случайные каналы.
Здесь мы видим, как библиотека хранит каналы, и численно проверяем эти результаты.

Всюду далее $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ и
$\mathcal N=\bigoplus_{d=1}^{C'}M_{k'_d}(\mathbb C)$ с минимальными центральными проекторами $P_c$ и
$Q_d$; $A=\bigoplus_cA_c$; $M_{k'\times k}$ — пространство комплексных матриц размера $k'\times k$;
$E_{ij}=|i\rangle\langle j|$ — матричные единицы; $\mathrm{Tr}$ — след `Tr_blunt`
$\sum_c\mathrm{Tr}\,A_c$, а $\mathrm{tr}_w=\sum_cw_c\mathrm{Tr}\,A_c$ — точный
след (теорема 3.1). Блочные матрицы действуют на $\mathbb C^n\otimes\mathbb C^k$, блочный индекс идёт первым.

<a id="s-6-1"></a>

## 6.1 Вполне положительные отображения и теорема Чоя

*Конспект:* [раздел 6.1](../../lecture_notes/ru/06_quantum_channels.md#s-6-1): [предложение 6.1](../../lecture_notes/ru/06_quantum_channels.md#r-6-1) (транспонирование), [теорема 6.2](../../lecture_notes/ru/06_quantum_channels.md#r-6-2) (Чой), [следствие 6.3](../../lecture_notes/ru/06_quantum_channels.md#r-6-3) (ранг Крауса) и [предложение 6.4](../../lecture_notes/ru/06_quantum_channels.md#r-6-4) (отображения между прямыми суммами).

В библиотеке два класса. `Channel` хранит **сохраняющие секторы** отображения ($C'=C$, $k'=k$,
$\Phi_{dc}=0$ при $d\ne c$), операторы Крауса которых $K_l=\bigoplus_cK^{cc}_l$ лежат в $\mathcal M$, в виде
тензора формы `(batch, r, C, k_max, k_max)`; `InterSectorChannel` хранит общее отображение из
[предложения 6.4](../../lecture_notes/ru/06_quantum_channels.md#r-6-4) в виде тензора формы `(batch, r, C', C, k'_max, k_max)`. `Channel.choi()` возвращает
$J(\Phi_{cc})$ для каждого сектора, а `Channel.superoperator()` — матрицу $S_c$ отображения $\Phi_{cc}$ на
векторах, развёрнутых по строкам, $\mathrm{vec}(X)_{(a,b)}=X_{ab}$. Обе содержат числа
$\Phi(E_{ij})_{ab}=J_{(i,a),(j,b)}=S_{(a,b),(i,j)}$, но в разном порядке. `Channel.from_superoperator`
выбрасывает `ValueError`, если матрица Чоя не эрмитова, $\|J-J^*\|>2\cdot10^{-10}s$ (в двойной точности), где
$s=\max(1,\|J\|)$ со спектральной нормой; в противном случае он заменяет $J$ его эрмитовой частью $(J+J^*)/2$,
диагонализует её, как в доказательстве [теоремы 6.2](../../lecture_notes/ru/06_quantum_channels.md#r-6-2), выбрасывает `ValueError`, если какое-то собственное значение меньше $-10^{-10}s$,
и полагает остальные отрицательные собственные значения равными $0$. Так как $J^*=J$ тогда и только тогда, когда $\Phi(X^*)=\Phi(X)^*$
для всех $X$ ([конспект лекций, раздел 6.1](../../lecture_notes/ru/06_quantum_channels.md#s-6-1)), с точностью до этих допусков на округление он отвергает в точности
отображения, не сохраняющие сопряжение, а среди остальных — в точности те, что не являются CP;
например, отвергаются $X\mapsto iX$ и транспонирование.

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, states, channels
from torch_vn_algebra.channels import Channel, InterSectorChannel

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def unit(i, j, k):
    X = torch.zeros(k, k, dtype=dt)
    X[i, j] = 1
    return X

def choi(f, k):
    """Choi matrix sum_ij E_ij (x) f(E_ij) of a linear map f on k x k matrices."""
    return sum(torch.kron(unit(i, j, k), f(unit(i, j, k)).contiguous()) for i in range(k) for j in range(k))

def to_superoperator(J, k):
    """S[(a,b),(i,j)] = J[(i,a),(j,b)], with a leading batch axis."""
    return J.reshape(k, k, k, k).permute(1, 3, 0, 2).reshape(1, k * k, k * k)

ev = torch.linalg.eigvalsh(choi(lambda X: X.T, 2))
print(ev)                                                 # the Choi matrix of the transpose is the swap
assert abs(ev[0] + 1) < 1e-12                             # Proposition 6.1

alg = TypeIAlgebra([2, 3], [2, 3], precision='double')    # M_2 (+) M_3
Phi = channels.random_channel(alg, kraus_rank=2)          # Kraus tensor of shape (1, 2, 2, 3, 3)

def block(Phi, c):
    """The block Phi_cc of a sector-preserving map, as a function on k_c x k_c matrices."""
    ks = Phi.algebra.k_factors
    def f(X):
        blocks = [torch.zeros(k, k, dtype=dt) for k in ks]
        blocks[c] = X
        return Phi(Phi.algebra.from_blocks(blocks)).matrix[0, c, :ks[c], :ks[c]]
    return f

J, S = Phi.choi(), Phi.superoperator()                    # one matrix per sector, (1, k_c^2, k_c^2)
rho = states.random_density_matrix(alg)
out = Phi(rho)
for c, k in enumerate(alg.k_factors):
    assert close(J[c][0], choi(block(Phi, c), k))         # the convention of Theorem 6.2
    w = torch.linalg.eigvalsh(J[c][0])
    assert w.min() > -1e-12 and int((w > 1e-10).sum()) == 2      # J >= 0, Kraus rank 2 (Corollary 6.3)
    assert close(S[c][0], to_superoperator(J[c][0], k)[0])
    assert close(S[c][0] @ rho.matrix[0, c, :k, :k].reshape(-1), out.matrix[0, c, :k, :k].reshape(-1))
assert Phi.is_trace_preserving()
assert close(Channel.from_superoperator(alg, S)(rho), out)       # back from the superoperator
try:
    Channel.from_superoperator(alg, [to_superoperator(choi(lambda X: X.T, k), k) for k in (2, 3)])
    raise AssertionError("the transpose was accepted")
except ValueError as e:
    print(e)                                              # superoperator is not completely positive
try:
    Channel.from_superoperator(alg, [1j * torch.eye(k * k, dtype=dt)[None] for k in (2, 3)])
    raise AssertionError("X -> iX was accepted")
except ValueError as e:
    print(e)                                              # superoperator does not preserve Hermiticity
```

<a id="s-6-2"></a>

## 6.2 Дилатация Стайнспринга и свобода в выборе операторов Крауса

*Конспект:* [раздел 6.2](../../lecture_notes/ru/06_quantum_channels.md#s-6-2): [теорема 6.5](../../lecture_notes/ru/06_quantum_channels.md#r-6-5) (Стайнспринг) и [теорема 6.7](../../lecture_notes/ru/06_quantum_channels.md#r-6-7) (унитарная свобода). Для блока $c=1$ канала из раздела 6.1 код проверяет, что $V$ — изометрия, что $\Phi(\rho)=\mathrm{Tr}_E(V\rho V^*)$ и $\Phi^*(B)=V^*(1\otimes B)V$, и что операторы Крауса $\sum_lu_{jl}K_l$ задают тот же канал.

В коде индекс окружения у $V$ — первый индекс тензора Крауса:

```python
K = Phi.kraus[0, :, 1, :3, :3]                            # the two Kraus operators of block c = 1
V = K.reshape(2 * 3, 3)                                   # V = sum_l e_l (x) K_l : C^3 -> C^2 (x) C^3
assert close(V.conj().T @ V, torch.eye(3, dtype=dt))      # TP <=> V is an isometry
r1 = rho.matrix[0, 1, :3, :3]
TrE = torch.einsum('lalb->ab', (V @ r1 @ V.conj().T).reshape(2, 3, 2, 3))
assert close(TrE, out.matrix[0, 1, :3, :3])               # Phi(rho) = Tr_E(V rho V^*)
B1 = torch.randn(3, 3, dtype=dt)
Bop = alg.from_blocks([torch.zeros(2, 2, dtype=dt), B1])
assert close(V.conj().T @ torch.kron(torch.eye(2, dtype=dt), B1) @ V,
             Phi.adjoint()(Bop).matrix[0, 1, :3, :3])    # Phi^*(B) = V^*(1 (x) B)V

u = alg.random_unitary(3)                                 # Theorem 6.7 with n = 3: pad with a zero operator
Kpad = torch.cat([Phi.kraus, torch.zeros_like(Phi.kraus[:, :1])], dim=1)
PhiL = Channel(alg, torch.einsum('jl,blcxy->bjcxy', u, Kpad))
assert close(PhiL(rho), out) and all(close(a, b) for a, b in zip(PhiL.choi(), J))
```

<a id="s-6-3"></a>

## 6.3 Сопряжённые относительно следа

*Конспект:* [раздел 6.3](../../lecture_notes/ru/06_quantum_channels.md#s-6-3): сопряжённое отображение $\Phi^\sharp$ относительно $(w,v)$, [предложение 6.8](../../lecture_notes/ru/06_quantum_channels.md#r-6-8) и [следствие 6.9](../../lecture_notes/ru/06_quantum_channels.md#r-6-9).

В библиотеке `Channel.adjoint()` — это $\Phi^*$, а `InterSectorChannel.adjoint(trace=...)` — это
$\Phi^\sharp$ для одного и того же именованного следа на обеих алгебрах (`'blunt'`, `'norm'`, `'tau_vN'`; веса
$1/(Ck_c)$ и $1/(C'k'_d)$ следа `tau_vN` учитывают ещё и числа секторов). Проверим [предложение 6.8](../../lecture_notes/ru/06_quantum_channels.md#r-6-8) для
случайного отображения $M_2\oplus M_3\to M_1\oplus M_2\oplus M_2$ и трёх следов, а также [следствие 6.9(a)](../../lecture_notes/ru/06_quantum_channels.md#r-6-9)
для потери частицы $\Phi(x)=(x_0+\mathrm{Tr}\,x_1)\oplus0$ на $M_1\oplus M_2$, которая является TP для следов blunt,
но не для `Tr_norm` ([конспект лекций, раздел 6.3](../../lecture_notes/ru/06_quantum_channels.md#s-6-3)):

```python
alg_out = TypeIAlgebra([1, 2, 2], [1, 2, 2], precision='double')
Psi = channels.random_inter_sector_channel(alg, alg_out, kraus_rank=2)     # Kraus tensor (1, 2, 3, 2, 2, 3)
x = alg.operator(torch.randn(1, 2, 3, 3, dtype=dt))
y = alg_out.operator(torch.randn(1, 3, 2, 2, dtype=dt))
weights = {'blunt': lambda k, C: 1.0, 'norm': lambda k, C: 1 / k, 'tau_vN': lambda k, C: 1 / (C * k)}
trace = {'blunt': lambda X: X.Tr_blunt(), 'norm': lambda X: X.Tr_norm(), 'tau_vN': lambda X: X.tau_vN()}

def z(a, t):                                              # the central element sum_c w_c P_c of trace t
    return a.central(torch.tensor([weights[t](k, a.C) for k in a.k_factors], dtype=torch.float64))

for t in trace:
    dual = Psi.adjoint(trace=t)
    assert (trace[t](Psi(x) @ y) - trace[t](x @ dual(y))).abs().max() < 1e-10   # the definition
    assert close(dual(y), z(alg, t).inverse() @ Psi.adjoint()(z(alg_out, t) @ y))   # Proposition 6.8
dev = [(Psi.adjoint(trace=t)(alg_out.identity()) - alg.identity()).matrix.abs().max().item() for t in trace]
print(dev)                         # TP only for the blunt traces: the other duals are not unital
assert Psi.is_trace_preserving() and dev[0] < 1e-12 and min(dev[1:]) > 0.1

alg01 = TypeIAlgebra([1, 2], [1, 2], precision='double')
loss = InterSectorChannel.from_blocks(alg01, alg01, {(0, 0): [torch.ones(1, 1)],
                                                     (0, 1): [torch.tensor([[1.0, 0.0]]), torch.tensor([[0.0, 1.0]])]})
assert loss.is_trace_preserving()
assert close(loss.adjoint(trace='norm')(alg01.identity()), alg01.central([1.0, 2.0]))
s = states.random_density_matrix(alg01)
assert abs(loss(s).Tr_norm() - (s @ alg01.central([1.0, 2.0])).Tr_norm()) < 1e-12   # Corollary 6.9(a)
print(s.Tr_norm().real.item(), loss(s).Tr_norm().real.item())      # Tr_norm is not conserved
```

<a id="s-6-4"></a>

## 6.4 Декогеренция

*Конспект:* [раздел 6.4](../../lecture_notes/ru/06_quantum_channels.md#s-6-4): [теорема 6.10](../../lecture_notes/ru/06_quantum_channels.md#r-6-10) (неподвижные точки унитальных каналов), [теорема 6.11](../../lecture_notes/ru/06_quantum_channels.md#r-6-11) (унитальные каналы увеличивают энтропию) и [предложение 6.12](../../lecture_notes/ru/06_quantum_channels.md#r-6-12) (матрица переходов).

Стандартные каналы библиотеки действуют внутри каждого блока ($p,\gamma\in[0,1]$):

| канал | $\Phi(\rho)_c$ | операторы Крауса в блоке $c$ |
|---|---|---|
| `dephasing_channel` | $(1-p)\rho_c+p\mathrm{diag}\,\rho_c$ | $\sqrt{1-p}\,1$, $\sqrt p\,E_{ii}$ |
| `depolarizing_channel` | $(1-p)\rho_c+p\mathrm{Tr}(\rho_c)\,1/k_c$ | $\sqrt{1-p}\,1$, $\sqrt{p/k_c}\,E_{ij}$ |
| `center_expectation` | $\mathrm{Tr}(\rho_c)\,1/k_c$ | $\sqrt{1/k_c}\,E_{ij}$ |
| `amplitude_damping_channel` | уровни $i>0$ распадаются в $0$ | $E_{00}+\sqrt{1-\gamma}\sum_{i>0}E_{ii}$, $\sqrt\gamma\,E_{0i}$ |
| `random_mixed_unitary_channel` | $\sum_jp_jU_j\rho U_j^*$ | $\sqrt{p_j}\,U_j$ |

Для каждого канала из таблицы на $M_2\oplus M_3$ код вычисляет размерность пространства его неподвижных точек
(ядра $S_c-1$) и размерность коммутанта его операторов Крауса: для унитальных каналов они совпадают
([теорема 6.10](../../lecture_notes/ru/06_quantum_channels.md#r-6-10)), а для амплитудного затухания размерности равны, но пространства различны. Затем он
проверяет [теорему 6.11](../../lecture_notes/ru/06_quantum_channels.md#r-6-11) на $100$ случайных состояниях и [предложение 6.12](../../lecture_notes/ru/06_quantum_channels.md#r-6-12) для отображения `Psi`
из раздела 6.3 и для потери частицы.

```python
def fixed_dims(Phi):
    """dim {X in M_{k_c} : Phi_cc(X) = X} for every sector: the null space of S_c - 1."""
    return [int((torch.linalg.svdvals(S_c[0] - torch.eye(k * k, dtype=dt)) < 1e-9).sum())
            for S_c, k in zip(Phi.superoperator(), Phi.algebra.k_factors)]

def commutant_dims(Phi):
    """dim of the commutant of {K_l, K_l^*} in every block (row-major vec, as in Chapter 2)."""
    out = []
    for c, k in enumerate(Phi.algebra.k_factors):
        I, Ks = torch.eye(k, dtype=dt), Phi.kraus[0, :, c, :k, :k]
        ops = [A.contiguous() for A in torch.cat([Ks, Ks.mH])]
        L = torch.cat([torch.kron(A, I) - torch.kron(I, A.T.contiguous()) for A in ops])
        out.append(k * k - int((torch.linalg.svdvals(L) > 1e-9).sum()))
    return out

zoo = {'dephasing': channels.dephasing_channel(alg, 0.3),
       'depolarizing': channels.depolarizing_channel(alg, 0.3),
       'centre': channels.center_expectation(alg),
       'mixed unitary': channels.random_mixed_unitary_channel(alg, 2),
       'amplitude damping': channels.amplitude_damping_channel(alg, 0.4)}
for name, ch in zoo.items():
    print(f"{name:18s} unital {ch.is_unital()!s:5s} fixed {fixed_dims(ch)} commutant {commutant_dims(ch)}")
    assert ch.is_trace_preserving()
    if ch.is_unital():
        assert fixed_dims(ch) == commutant_dims(ch)                         # Theorem 6.10
assert fixed_dims(zoo['dephasing']) == [2, 3] and fixed_dims(zoo['depolarizing']) == [1, 1]
E00 = alg.from_blocks([unit(0, 0, 2), unit(0, 0, 3)])
E01 = alg.from_blocks([unit(0, 1, 2), unit(0, 1, 3)])
assert close(zoo['amplitude damping'](E00), E00) and not close(E00 @ E01, E01 @ E00)   # equal dims, different spaces
assert fixed_dims(zoo['amplitude damping']) == commutant_dims(zoo['amplitude damping']) == [1, 1]

rhos = states.random_density_matrix(alg, batch_size=100)
for name in ('dephasing', 'depolarizing', 'centre', 'mixed unitary'):                 # Theorem 6.11
    assert (zoo[name](rhos).von_neumann_entropy() >= rhos.von_neumann_entropy() - 1e-12).all()
    assert (zoo[name](rhos).purity() <= rhos.purity() + 1e-12).all()
mm = states.maximally_mixed_state(alg)
S_after = channels.amplitude_damping_channel(alg, 1.0)(mm).von_neumann_entropy().item()
print(mm.von_neumann_entropy().item(), S_after)          # not unital: the entropy drops
assert abs(S_after + 0.4 * math.log(0.4) + 0.6 * math.log(0.6)) < 1e-12

T = Psi.transition_matrix()[0]                            # (C', C) = (3, 2)
p, q = rho.sector_probabilities()[0], Psi(rho).sector_probabilities()[0]
assert torch.allclose(T.sum(0), torch.ones(2, dtype=T.dtype))                       # Proposition 6.12
assert torch.allclose(q, Psi.transition_matrix(rho)[0] @ p)
classical = channels.center_expectation(alg_out) @ Psi @ channels.center_expectation(alg)
Tp = T @ p
assert close(classical(rho), alg_out.central([Tp[d] / k for d, k in enumerate(alg_out.k_factors)]))
assert torch.allclose(loss.transition_matrix()[0], torch.tensor([[1.0, 1.0], [0.0, 0.0]], dtype=T.dtype))
```

При полном амплитудном затухании энтропия максимально смешанного состояния $M_2\oplus M_3$ падает
с $\log5\approx1.609$ до $H(\tfrac25,\tfrac35)\approx0.673$: неунитальный канал может очищать состояние.

На $M_2$ дефазировка $D_p$ и амплитудное затухание $A_\gamma$ коммутируют, и их композиция умножает
когерентности на $(1-p)\sqrt{1-\gamma}$; при $1-p=e^{-t/T_\varphi}$ и $1-\gamma=e^{-t/T_1}$ это даёт
$1/T_2=1/(2T_1)+1/T_\varphi$ ([конспект лекций, раздел 6.4](../../lecture_notes/ru/06_quantum_channels.md#s-6-4), *Когерентности кубита*). Код проверяет эти
тождества и законы композиции $D_{p(t)}D_{p(s)}=D_{p(t+s)}$, $A_{\gamma(t)}A_{\gamma(s)}=A_{\gamma(t+s)}$:

```python
q2 = TypeIAlgebra([2], [2], precision='double')
p_, g_ = 0.3, 0.4
Dp, Ag = channels.dephasing_channel(q2, p_), channels.amplitude_damping_channel(q2, g_)
r2 = states.random_density_matrix(q2)
assert close(Dp(Ag(r2)), Ag(Dp(r2)))                                         # D_p A_gamma = A_gamma D_p
coh = Dp(Ag(r2)).matrix[0, 0, 0, 1] / r2.matrix[0, 0, 0, 1]
assert abs(coh - (1 - p_) * math.sqrt(1 - g_)) < 1e-12                       # coherence factor
T1, Tphi, t1, t2 = 2.0, 3.0, 0.7, 1.1
gam = lambda t: 1 - math.exp(-t / T1)
assert close(channels.amplitude_damping_channel(q2, gam(t1))(channels.amplitude_damping_channel(q2, gam(t2))(r2)),
             channels.amplitude_damping_channel(q2, gam(t1 + t2))(r2))       # A_gamma(t) A_gamma(s) = A_gamma(t+s)
pdeph = lambda t: 1 - math.exp(-t / Tphi)
assert close(channels.dephasing_channel(q2, pdeph(t1))(channels.dephasing_channel(q2, pdeph(t2))(r2)),
             channels.dephasing_channel(q2, pdeph(t1 + t2))(r2))           # D_p(t) D_p(s) = D_p(t+s)
both = channels.dephasing_channel(q2, pdeph(t1))(channels.amplitude_damping_channel(q2, gam(t1))(r2))
coh_t = (both.matrix[0, 0, 0, 1] / r2.matrix[0, 0, 0, 1]).item()
assert abs(coh_t - math.exp(-t1 * (1 / (2 * T1) + 1 / Tphi))) < 1e-12       # 1/T_2 = 1/(2T_1) + 1/T_phi
```

<a id="s-6-5"></a>

## 6.5 Случайные каналы

*Конспект:* [раздел 6.5](../../lecture_notes/ru/06_quantum_channels.md#s-6-5): [лемма 6.13](../../lecture_notes/ru/06_quantum_channels.md#r-6-13) и [теорема 6.15](../../lecture_notes/ru/06_quantum_channels.md#r-6-15) (выходная чистота случайного канала).

**Изометрия, случайная по Хаару**, $V:\mathbb C^k\to\mathbb C^N$ — это, как и в библиотеке (`random_isometry`),
множитель $Q$ в $G=QR$, где $G$ — матрица Жинибра размера $N\times k$ (независимые элементы с
плотностью $\pi^{-1}e^{-|z|^2}$), $Q^*Q=1$, а $R$ — верхнетреугольная с положительной диагональю (на
CPU библиотека берёт разложение LAPACK $G=Q'R'$ и умножает $i$-й столбец $Q'$ на
$R'_{ii}/|R'_{ii}|$, что, как в теореме 5.3, даёт это разложение, если все $R'_{ii}\ne0$, т. е.
если $G$ имеет ранг $k$, а это выполняется почти наверное по доказательству леммы 6.13; для больших пакетов на
GPU она использует процесс Грама–Шмидта, который даёт его непосредственно).
`random_channel(alg, r)` независимо в каждом блоке генерирует такие $V_c:\mathbb C^{k_c}\to\mathbb C^r\otimes\mathbb C^{k_c}$
и использует операторы Крауса из теоремы 6.5.

По [теореме 6.15 конспекта лекций](../../lecture_notes/ru/06_quantum_channels.md#r-6-15), для чистого состояния $\psi$ в секторе размера $k$ и канала $\Phi$,
порождённого `random_channel` с рангом Крауса $r$, $\mathbb E\mathrm{Tr}\,\Phi(\psi)^2=(k+r)/(kr+1)$ и
$\mathbb E\,\Phi(\psi)=1_k/k$. Оценим обе величины по $20\,000$ выборкам со стандартной ошибкой среднего:

```python
B = 20_000
for r in (1, 2, 4):
    Phi_r = channels.random_channel(alg, kraus_rank=r, batch_size=B)
    for c, k in enumerate(alg.k_factors):
        psi = torch.randn(k, 1, dtype=dt)
        psi = psi / psi.norm()
        out_r = Phi_r(alg.from_blocks([psi @ psi.conj().T if d == c else torch.zeros(kk, kk, dtype=dt)
                                       for d, kk in enumerate(alg.k_factors)]))     # a batch of B outputs
        purity = (out_r @ out_r).trace.real
        mean, err = purity.mean().item(), purity.std().item() / math.sqrt(B)
        exact = (k + r) / (k * r + 1)
        print(f"k = {k}, r = {r}: E Tr Phi(psi)^2 = {mean:.4f} +- {err:.4f}   (exact {exact:.4f})")
        assert abs(mean - exact) < 5 * err + 1e-12
        avg = out_r.matrix.mean(0)[c, :k, :k]
        assert (avg - torch.eye(k, dtype=dt) / k).abs().max() < 0.02          # E Phi(psi) = 1/k
```

<a id="s-6-6"></a>

## 6.6 С чем сравнивать

QuTiP (`to_choi`, `to_kraus`, `to_super`, `to_stinespring`), Qiskit (`qiskit.quantum_info`:
классы `Kraus`, `Choi`, `SuperOp`, `Stinespring` и `random_quantum_channel`; модели шума, такие как
`amplitude_damping_error` в Qiskit Aer), toqito (`toqito.channel_ops`, `toqito.channel_props`) и
QETLAB (`ChoiMatrix`, `KrausOperators`, `ApplyMap`, `RandomSuperoperator`) переходят между
представлениями из разделов 6.1–6.2, проверяют полную положительность и (большинство из них) генерируют случайные
каналы — для отображений $M_n\to M_m$ между полными матричными алгебрами.
Их соглашения различаются (развёртка по столбцам или по строкам, порядок сомножителей в матрице Чоя), поэтому
сравнивайте соглашения, прежде чем сравнивать числа. Сохраняющий секторы канал на $\mathcal M$ можно
рассматривать там как канал на $M_D$, $D=\sum_ck_c$, с супероператорами размера $D^2\times D^2$
вместо $k_c^2\times k_c^2$ на сектор. Библиотека добавляет саму блочную структуру: операторы Крауса,
хранимые по секторам и пакетами, отображения между разными алгебрами (`InterSectorChannel`),
сопряжённые относительно каждого из следов главы 3 и матрицу переходов между секторами.

<a id="x-6"></a>

## Упражнения

[Упражнения](../../lecture_notes/ru/06_quantum_channels.md#x-6) конспекта лекций, численно:

1. (Вернер–Холево) Проверьте с помощью `Channel.from_superoperator`, что $\Phi(\rho)=\big(\mathrm{Tr}(\rho)1-\rho^T\big)/(k-1)$
   на $M_3$ — унитальный канал с рангом Крауса $3$.
2. Проверьте, что канал дефазировки на $M_2\oplus M_3$ имеет в блоке $c$ ранг Крауса $k_c$, тогда как библиотека
   хранит $1+k_{\max}$ операторов Крауса.
3. Проверьте, что смешанно-унитарный канал с операторами Крауса $X^aZ^b/k$ — это `center_expectation` на $M_3$.
4. Проверьте $(\Psi\circ\Phi)^\sharp=\Phi^\sharp\circ\Psi^\sharp$ для трёх следов библиотеки.

<details>
<summary>Решения</summary>

*Упражнение 1.*

```python
k = 3
M3 = TypeIAlgebra([k], [k], precision='double')
J_WH = choi(lambda X: (X.trace() * torch.eye(k, dtype=dt) - X.T) / (k - 1), k)
w = torch.linalg.eigvalsh(J_WH)
assert w.min() > -1e-12 and int((w > 1e-10).sum()) == k * (k - 1) // 2
WH = Channel.from_superoperator(M3, [to_superoperator(J_WH, k)])
assert WH.is_trace_preserving() and WH.is_unital()
r3 = states.random_density_matrix(M3).matrix[0, 0]
assert close(WH(M3.operator(r3[None, None])).matrix[0, 0], (torch.eye(k, dtype=dt) - r3.T) / (k - 1))
```

*Упражнение 2.*

```python
ranks = [int((torch.linalg.eigvalsh(J_c[0]) > 1e-10).sum()) for J_c in zoo['dephasing'].choi()]
print(ranks, zoo['dephasing'].kraus_rank)                 # [2, 3] versus 4 stored operators
assert ranks == alg.k_factors and zoo['dephasing'].kraus_rank == 1 + max(alg.k_factors)
```

*Упражнение 3.*

```python
Xs = torch.roll(torch.eye(k, dtype=dt), 1, dims=0)
Zs = torch.diag(torch.exp(2j * math.pi * torch.arange(k, dtype=torch.float64) / k)).to(dt)
W = [torch.linalg.matrix_power(Xs, a) @ torch.linalg.matrix_power(Zs, b) / k for a in range(k) for b in range(k)]
twirl = channels.kraus_channel([M3.operator(Wab[None, None]) for Wab in W])
r3op = states.random_density_matrix(M3)
assert close(twirl(r3op), channels.center_expectation(M3)(r3op))
```

*Упражнение 4.*

```python
Psi2 = channels.random_inter_sector_channel(alg_out, alg, kraus_rank=2)
for t in trace:
    assert close((Psi2 @ Psi).adjoint(trace=t)(x), Psi.adjoint(trace=t)(Psi2.adjoint(trace=t)(x)))
```

</details>

**Далее:** [Руководство 7](07_tensor_products.md): тензорные произведения, частичные следы и запутанность в библиотеке.
