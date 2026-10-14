[English](../08_dynamics.md) · [Español](../es/08_dynamics.md) · [Français](../fr/08_dynamics.md) · [Deutsch](../de/08_dynamics.md) · [中文](../zh/08_dynamics.md) · [日本語](../ja/08_dynamics.md) · **Русский**

# Руководство 8. Динамика

Это руководство сопровождает [главу 8](../../lecture_notes/ru/08_dynamics.md) конспекта лекций, в которой доказывается, что гамильтониан
из $\mathcal M$ сохраняет веса секторов ([теорема 8.3](../../lecture_notes/ru/08_dynamics.md#r-8-3)), что генераторы формы GKSL порождают
квантовые динамические полугруппы ([теорема 8.6](../../lecture_notes/ru/08_dynamics.md#r-8-6)), и рассматриваются скачки между секторами,
постселекция и стационарные состояния. Мы используем точные решатели и решатели Рунге–Кутты из `dynamics` для уравнений
Шрёдингера, фон Неймана и Линдблада и сверяем результаты со SciPy и с генераторами GKSL на $M_D$.

<a id="s-8-1"></a>

## 8.1 Замкнутые системы

*Конспект:* [раздел 8.1](../../lecture_notes/ru/08_dynamics.md#s-8-1); [предложение 8.2](../../lecture_notes/ru/08_dynamics.md#r-8-2) (пропагатор $U(t)=e^{-itH}$), [теорема 8.3](../../lecture_notes/ru/08_dynamics.md#r-8-3) (сохранение заряда).

По предложению 8.2(1) эволюция распадается на $C$ независимых задач размеров $k_c$, по одной на блок.
Библиотека вычисляет $U(t)$ с помощью одного спектрального разложения на каждый размер блока (`dynamics.propagator`), а
`dynamics.von_neumann` эволюционирует плотность сразу ко многим моментам времени в собственном базисе, где
$\rho(t)_{ij}=e^{-it(\lambda_i-\lambda_j)}\rho_{ij}$. Проверяем по `scipy.linalg.expm`:

```python
import math
import numpy as np
import scipy.integrate
import scipy.linalg as sla
import torch
from torch_vn_algebra import TypeIAlgebra, DensityMatrix, InterSectorChannel, channels, dynamics, states

torch.manual_seed(0)
f64, c128 = torch.float64, torch.complex128
alg = TypeIAlgebra([2, 3, 4], [2, 3, 4], precision='double')     # M_2 (+) M_3 (+) M_4, D = 9
H = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
t = 0.7
U = dynamics.propagator(H, t)                                      # exp(-itH), block by block
for c, k in enumerate(alg.k_factors):                              # independent reference
    ref = sla.expm(-1j * t * H.matrix[0, c, :k, :k].numpy())
    assert np.abs(U.matrix[0, c, :k, :k].numpy() - ref).max() < 1e-12
assert (U.adjoint() @ U - alg.identity()).matrix.abs().max() < 1e-12

rho0 = states.random_density_matrix(alg)
times = [0.0, 0.5, 1.0, 2.0]
traj = dynamics.von_neumann(rho0, H, times)                        # exact, all times at once
for s, rho_s in zip(times, traj):
    U_s = dynamics.propagator(H, s)
    assert (rho_s.matrix - (U_s @ rho0 @ U_s.adjoint()).matrix).abs().max() < 1e-12
    assert (rho_s.sector_probabilities() - rho0.sector_probabilities()).abs().max() < 1e-12   # Theorem 8.3
A = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
A_t = U.adjoint() @ A @ U                                          # Heisenberg picture
assert abs(rho0.expectation(A_t) - dynamics.von_neumann(rho0, H, [t])[0].expectation(A)) < 1e-12
```

Единичный вектор $\psi=\bigoplus_c\psi_c$ задаёт на $\mathcal M$ состояние $\bigoplus_c|\psi_c\rangle\langle\psi_c|$, которое
эволюционирует под действием $U(t)$ ([конспект лекций, раздел 8.1](../../lecture_notes/ru/08_dynamics.md#s-8-1)), поэтому `dynamics.schrodinger` (векторы в блочной раскладке `(batch, C, k_max)`) и `von_neumann`
должны на нём совпадать. Для зависящего от времени $H(t)$
библиотека использует классический метод Рунге–Кутты с фиксированным шагом (`dynamics.rk4`, `schrodinger_rk4`): шаг умножает
$\psi$ на многочлен от значений $H(t)\in\mathcal M$, поэтому он никогда не переносит амплитуду между секторами,
но он не унитарен: при постоянном $H$ шаг размера $h$ умножает собственный вектор $H$ с собственным значением
$\lambda$ на $T(-ih\lambda)$, где $T(z)=\sum_{k=0}^4z^k/k!$, а $|T(ix)|^2=1-\frac{x^6}{72}+\frac{x^8}{576}<1$ при
$0<|x|<2\sqrt2$ (раскройте скобки). Поэтому RK4 не сохраняет веса секторов точно: при постоянном $H$ компонента сектора
с ненулевым коэффициентом при собственном векторе с $0<|h\lambda|<2\sqrt2$ теряет вес на каждом шаге.
Насколько они меняются за один запуск, — численное наблюдение (см. ниже).
(Для одного большого сектора `krylov.evolve` нужны только произведения матрицы на вектор; см. [Руководство 9, раздел 9.5](../../tutorial/ru/09_spin_chains.md#s-9-5).)

```python
rows = torch.tensor([[1, 1, 0, 0], [1, 1, 1, 0], [1, 1, 1, 1]], dtype=c128)   # entries of C^2 (+) C^3 (+) C^4
psi0 = torch.randn(1, alg.C, alg.k_max, dtype=c128) * rows
psi0 = psi0 / psi0.norm()
psi = dynamics.schrodinger(psi0, H, times)                         # shape (T, batch, C, k_max)
p_vec = (psi.abs() ** 2).sum(-1)                                   # ||P_c psi(t)||^2
assert (p_vec - p_vec[0]).abs().max() < 1e-12
mix = DensityMatrix(alg, matrix=psi0[..., :, None] * psi0[..., None, :].conj())   # the state of psi0 on M
outer = psi[..., :, None] * psi[..., None, :].conj()
assert all((o - r.matrix).abs().max() < 1e-12 for o, r in zip(outer, dynamics.von_neumann(mix, H, times)))

V = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
H_of_t = lambda s: H.matrix + math.cos(3 * s) * V.matrix          # driven, in M at every time
psi_rk = dynamics.schrodinger_rk4(psi0, H_of_t, times, substeps=50, algebra=alg)
drift = ((psi_rk.abs() ** 2).sum(-1) - p_vec[0]).abs().max().item()
print(f"RK4, time-dependent H: largest change of a sector weight {drift:.1e}")
assert drift < 1e-8 and (dynamics.schrodinger_rk4(psi0, H, times, substeps=50) - psi).abs().max() < 1e-7

```

*Численное наблюдение*: в этом запуске веса секторов меняются на $1.1\cdot10^{-9}$ при RK4 (напечатанное
значение) и меньше чем на допуск $10^{-12}$ проверок `assert` при точных методах решения, для которых
теорема 8.3 даёт точное сохранение.

<a id="s-8-2"></a>

## 8.2 Квантовые динамические полугруппы

*Конспект:* [раздел 8.2](../../lecture_notes/ru/08_dynamics.md#s-8-2); [теорема 8.6](../../lecture_notes/ru/08_dynamics.md#r-8-6) (генераторы GKSL и разложение Дайсона), [теорема 8.7](../../lecture_notes/ru/08_dynamics.md#r-8-7).

В библиотеке четыре процедуры для уравнения Линдблада: `lindblad_rhs` возвращает правую часть в виде функции для
`rk4`; `lindblad_evolve` интегрирует методом RK4; `lindblad_superoperator` возвращает для каждого сектора
матрицу $k_c^2\times k_c^2$ отображения $\mathcal L$ на векторах $\mathrm{vec}\,\rho_c$, развёрнутых по строкам; а
`lindblad_channel` вычисляет $e^{t\mathcal L}$ точно, сектор за сектором, в виде `Channel`. Проверим
теорему 8.6 для произвольных скачков на $M_3$: матрица Чоя $e^{t\mathcal L}$ положительна (теорема
6.2), отображение является TP, выполняется полугрупповой закон, RK4 согласуется с экспонентой, а скачковая часть
$e^{t\mathcal L}-S_0(t)=\sum_{n\ge1}S_n(t)$ сама является CP:

```python
M3 = TypeIAlgebra([3], [3], precision='double')
H3 = M3.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
jumps = [M3.operator(torch.randn(1, 1, 3, 3, dtype=c128) / 2) for _ in range(2)]     # arbitrary L_j
T = lambda s: dynamics.lindblad_channel(M3, H3, jumps, t=s)        # exp(sL), exact
S = lambda s: T(s).superoperator()[0][0]                           # its 9 x 9 matrix on row-major vec
choi = lambda X: X.reshape(3, 3, 3, 3).permute(2, 0, 3, 1).reshape(9, 9)   # J_{(i,a),(j,b)} = S_{(a,b),(i,j)}
assert (choi(S(0.8)) - T(0.8).choi()[0][0]).abs().max() < 1e-12
assert torch.linalg.eigvalsh(choi(S(0.8))).min() > -1e-12         # completely positive
assert T(0.8).trace_preservation_error().max() < 1e-12             # trace preserving
assert (S(0.3) @ S(0.5) - S(0.8)).abs().max() < 1e-12              # semigroup

rho3 = states.random_density_matrix(M3)
rk = dynamics.lindblad_evolve(rho3, H3, jumps, [0.0, 0.8], substeps=40)[-1]
assert (rk.matrix - T(0.8)(rho3).matrix).abs().max() < 1e-7        # RK4 against exp(tL)
assert abs(dynamics.lindblad_rhs(H3, jumps)(0.0, rho3.matrix)[0, 0].trace()) < 1e-12   # Tr L(rho) = 0

H_eff = H3.matrix - 0.5j * sum(L.adjoint().matrix @ L.matrix for L in jumps)
S0 = channels.Channel(M3, torch.linalg.matrix_exp(-0.8j * H_eff)[:, None]).superoperator()[0][0]
jump_part = torch.linalg.eigvalsh(choi(S(0.8) - S0))
print(f"smallest eigenvalue of the Choi matrix of the jump part: {jump_part.min().item():.1e}")
assert jump_part.min() > -1e-12                                    # sum_{n >= 1} S_n(t) is CP
```

<a id="s-8-3"></a>

## 8.3 Уравнения Линдблада в алгебре с секторами

*Конспект:* [раздел 8.3](../../lecture_notes/ru/08_dynamics.md#s-8-3): генератор $\mathcal L$ со скачками внутри секторов и между ними, который реализует `dynamics.lindblad_rhs`; [предложение 8.8](../../lecture_notes/ru/08_dynamics.md#r-8-8), [предложение 8.9](../../lecture_notes/ru/08_dynamics.md#r-8-9) (веса секторов).

**Пример: потеря частиц** ([конспект лекций, раздел 8.3](../../lecture_notes/ru/08_dynamics.md#s-8-3)). Две бозонные моды, не более двух
частиц, $\mathcal M=M_1\oplus M_2\oplus M_3$ (секторы $N=0,1,2$), димер Бозе–Хаббарда
$H=J(a_1^*a_2+a_2^*a_1)+\frac U2\sum_jn_j(n_j-1)$, и каждая мода теряет частицы со скоростью $\gamma$: это один
`InterSectorChannel` с операторами Крауса $a_1,a_2$, ограниченными на $N\to N-1$. По
[предложению 8.9(2)](../../lecture_notes/ru/08_dynamics.md#r-8-9) веса секторов подчиняются уравнению $p'=Qp$ с $Q=R-\mathrm{diag}(\sum_dR_{dc})_c$ при скорости $1$,
$R$ = `transition_matrix()`, и при старте из $N=2$ следуют биномиальному закону $p_2=q^2$, $p_1=2q(1-q)$, $p_0=(1-q)^2$,
$q=e^{-\gamma t}$, что бы ни делали перескоки и взаимодействие. Код проверяет
$\Phi^*(1)=0\oplus1\oplus2$, т. е. $a_1^*a_1+a_2^*a_2=N$ на используемых блоках, и сверяет $p(t)$ с обеими формулами.

```python
fock = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double', charges=[0, 1, 2])   # N = 0, 1, 2
r2 = math.sqrt(2)
a1 = {(0, 1): [[1.0, 0.0]], (1, 2): [[r2, 0.0, 0.0], [0.0, 1.0, 0.0]]}           # a_1 : N -> N - 1
a2 = {(0, 1): [[0.0, 1.0]], (1, 2): [[0.0, 1.0, 0.0], [0.0, 0.0, r2]]}           # a_2 : N -> N - 1
blocks = lambda a: {dc: torch.tensor(K, dtype=f64) for dc, K in a.items()}
loss = InterSectorChannel.from_blocks(fock, fock, {dc: [blocks(a1)[dc], blocks(a2)[dc]] for dc in a1})
J, U_int = 1.0, 0.8                                                # hopping and interaction
H_bh = fock.from_blocks([torch.zeros(1, 1, dtype=f64), torch.tensor([[0.0, J], [J, 0.0]], dtype=f64),
                         torch.tensor([[U_int, r2 * J, 0.0], [r2 * J, 0.0, r2 * J], [0.0, r2 * J, U_int]], dtype=f64)])

assert (loss.adjoint()(fock.identity()).matrix - fock.central([0.0, 1.0, 2.0]).matrix).abs().max() < 1e-12
R = loss.transition_matrix()[0]                                    # r_dc (columns: from sector c)
Q = R - torch.diag(R.sum(0))                                       # Proposition 8.9(2), for rate 1
print(Q)

gamma = 0.5
psi2 = torch.randn(3, dtype=c128)
psi2 = psi2 / psi2.norm()
m0 = torch.zeros(1, 3, 3, 3, dtype=c128)
m0[0, 2] = torch.outer(psi2, psi2.conj())                          # a pure state with N = 2
ts = [0.0, 0.5, 1.0, 2.0, 4.0]
traj = dynamics.lindblad_evolve(DensityMatrix(fock, matrix=m0), H_bh, [loss], ts, rates=[gamma], substeps=100)
for s, rho_s in zip(ts, traj):
    q = math.exp(-gamma * s)
    p = rho_s.sector_probabilities()[0]
    assert (p - torch.tensor([(1 - q) ** 2, 2 * q * (1 - q), q ** 2], dtype=f64)).abs().max() < 1e-8
    assert (p - torch.from_numpy(sla.expm(gamma * s * Q.numpy())[:, 2])).abs().max() < 1e-8
neg = min(w.min().item() for rho_s in traj for w in rho_s.eigenvalues())
print(f"smallest eigenvalue of rho(t) under RK4: {neg:.1e}")
assert neg > -1e-7                                                 # positive up to the RK4 error

sm = torch.tensor([[0.0, 1.0], [0.0, 0.0]], dtype=c128)           # |g><e|; H = 0, rate 3, one step h = 1
decay = lambda s, X: 3.0 * (sm @ X @ sm.mH - 0.5 * (sm.mH @ sm @ X + X @ sm.mH @ sm))
step = dynamics.rk4(decay, torch.diag(torch.tensor([0.0, 1.0], dtype=c128)), [0.0, 1.0], substeps=1)[-1]
assert (step - torch.diag(torch.tensor([-3 / 8, 11 / 8], dtype=c128))).abs().max() < 1e-14   # not positive
```

RK4 не обязан сохранять положительность: для $\rho'=\mathcal L\rho$ шаг размера $h$ — это
$\rho\mapsto T(h\mathcal L)\rho$ с многочленом $T$ из раздела 8.1 (подставьте четыре стадии); это отображение
сохраняет след (так как $\mathrm{Tr}\,\mathcal L(X)=0$ для всех $X$), но не обязано быть положительным. Например, на $M_2$ с
базисом $|g\rangle,|e\rangle$, $H=0$ и одним скачком $|g\rangle\langle e|$ со скоростью $\gamma$ генератор $\mathcal L$ переводит
$\mathrm{diag}(a,b)$ в $\gamma b\mathrm{diag}(1,-1)$, поэтому один шаг с $h\gamma=3$ переводит $\mathrm{diag}(0,1)$ в
$\mathrm{diag}(1-T(-3),T(-3))=\mathrm{diag}(-\frac38,\frac{11}8)$. *Численное наблюдение*: в этом запуске наименьшее собственное значение
$\rho(t)$ равно $\approx-2\cdot10^{-8}$ (напечатанное значение, при шагах не больше $0.02$). `lindblad_channel` точен,
но работает только со скачками внутри секторов. Наконец, проверим предложение 8.8(2) для случайного отображения между
секторами и предложение 8.8(3) для потери частиц, сравнивая с генераторами GKSL на $M_6$, построенными средствами numpy:

```python
Phi = channels.random_inter_sector_channel(fock, fock, kraus_rank=2)
rho = states.random_density_matrix(fock)
on_M = dynamics.lindblad_rhs(None, [Phi], algebra=fock)(0.0, rho.matrix)[0].numpy()
ks, K = fock.k_factors, Phi.kraus[0].numpy()                       # K[l, d, c] = K_l^{dc}, zero padded
embed = lambda Y: sla.block_diag(*[Y[c, :k, :k] for c, k in enumerate(ks)])    # (+)_c Y_c in M_6
Ls = [np.block([[K[l, d, c, :kd, :kc] for c, kc in enumerate(ks)] for d, kd in enumerate(ks)])
      for l in range(2)]                                           # L_l = sum_{d,c} J_l^{dc}
X = embed(rho.matrix[0].numpy())
full = sum(L @ X @ L.conj().T - 0.5 * (L.conj().T @ L @ X + X @ L.conj().T @ L) for L in Ls)   # GKSL on M_6
E_full = full * embed(np.ones((3, 3, 3)))                          # E(X) = sum_d P_d X P_d
assert np.abs(E_full - embed(on_M)).max() < 1e-12                  # Proposition 8.8(2)
print(f"coherences between sectors created on M_6: {np.abs(full - E_full).max():.2f}")

zero = lambda d, c: np.zeros((ks[d], ks[c]))
L_loss = [np.block([[blocks(a)[(d, c)].numpy() if (d, c) in a else zero(d, c) for c in range(3)] for d in range(3)])
          for a in (a1, a2)]                                       # a_1, a_2 on M_6: one block per block row and column
full_loss = sum(L @ X @ L.T - 0.5 * (L.T @ L @ X + X @ L.T @ L) for L in L_loss)
loss_on_M = dynamics.lindblad_rhs(None, [loss], algebra=fock)(0.0, rho.matrix)[0].numpy()
assert np.abs(full_loss - embed(loss_on_M)).max() < 1e-12         # Proposition 8.8(3): no compression needed
```

<a id="s-8-4"></a>

## 8.4 Постселекция и эволюция без скачков

*Конспект:* [раздел 8.4](../../lecture_notes/ru/08_dynamics.md#s-8-4); [предложение 8.10](../../lecture_notes/ru/08_dynamics.md#r-8-10) (вероятность отсутствия скачков $w_0$), [предложение 8.11](../../lecture_notes/ru/08_dynamics.md#r-8-11) (постселекция на «ни одна частица не потеряна»).

Функция библиотеки `states.lueders_update(rho, P)` возвращает $P\rho P/\mathrm{Tr}(P\rho P)$ и
$\mathrm{Tr}(\rho P)$. Пусть теперь две моды теряют частицы с разными скоростями $\gamma_1>\gamma_2$; тогда
$A=\gamma_1n_1+\gamma_2n_2$ не централен (в секторе $N=2$ это $\mathrm{diag}(2\gamma_1,\gamma_1+\gamma_2,2\gamma_2)$ в
базисе $|20\rangle,|11\rangle,|02\rangle$), предложение 8.9(2) неприменимо, и постселектированное состояние на $N=2$ —
это нормированная эволюция под действием $H_{\mathrm{eff}}$, антиэрмитова часть $-\frac i2A$ которого не кратна $1$ в этом секторе. Проверим предложение 8.11 для $n=0$ и, вычислив
интеграл для $S_1$ по формуле Симпсона, для $n=1$:

```python
g1, g2 = 0.6, 0.2                                                  # mode 1 leaks faster
loss1 = InterSectorChannel.from_blocks(fock, fock, {dc: [K] for dc, K in blocks(a1).items()})
loss2 = InterSectorChannel.from_blocks(fock, fock, {dc: [K] for dc, K in blocks(a2).items()})
m1 = torch.zeros(1, 3, 3, 3, dtype=c128)
m1[0, 2, 0, 0] = 1.0                                               # |20>: both particles in mode 1
ts = [0.0, 0.5, 1.0, 2.0]
traj = dynamics.lindblad_evolve(DensityMatrix(fock, matrix=m1), H_bh, [loss1, loss2], ts, rates=[g1, g2],
                                substeps=100)

n_1 = [np.diag(x) for x in ([0.0], [1.0, 0.0], [2.0, 1.0, 0.0])]   # n_1 in each sector; n_2 = N - n_1
H_eff = [H_bh.matrix[0, N, :N + 1, :N + 1].numpy() - 0.5j * (g1 * n + g2 * (N * np.eye(N + 1) - n))
         for N, n in enumerate(n_1)]

def S0(s, N, X):                                  # e^{-isH_eff} X e^{isH_eff^*} in sector N, s a time or an array
    E = torch.linalg.matrix_exp(torch.from_numpy(-1j * np.asarray(s)[..., None, None] * H_eff[N])).numpy()
    return E @ X @ E.conj().swapaxes(-1, -2)

rho_20 = m1[0, 2].numpy()
for s, rho_s in zip(ts, traj):
    w0 = np.trace(S0(s, 2, rho_20)).real                           # probability of no loss
    post, prob = states.lueders_update(rho_s, fock.central([0.0, 0.0, 1.0]))   # post-select N = 2
    assert abs(prob.item() - w0) < 1e-8                            # Proposition 8.11, n = 0
    assert np.abs(post.matrix[0, 2].numpy() - S0(s, 2, rho_20) / w0).max() < 1e-7
closed = dynamics.von_neumann(DensityMatrix(fock, matrix=m1), H_bh, [2.0])[0]   # no losses at all
dist = 0.5 * np.abs(np.linalg.eigvalsh(S0(2.0, 2, rho_20) / w0 - closed.matrix[0, 2].numpy())).sum()
print(f"t = 2: P(no loss) = {w0:.3f} (without hopping: {math.exp(-2 * g1 * 2):.3f}); "
      f"trace distance post-selected/closed system: {dist:.3f}")

K1, K2 = blocks(a1)[(1, 2)].numpy(), blocks(a2)[(1, 2)].numpy()   # S_1(2) rho = int_0^2 S_0(2-s) Phi S_0(s) rho ds
jump = lambda X: g1 * K1 @ X @ K1.T + g2 * K2 @ X @ K2.T
s_grid = np.linspace(0.0, 2.0, 201)
S1 = scipy.integrate.simpson(S0(2.0 - s_grid, 1, jump(S0(s_grid, 2, rho_20))), x=s_grid, axis=0)
assert np.abs(S1 - traj[-1].matrix[0, 1, :2, :2].numpy()).max() < 1e-8     # Proposition 8.11, n = 1
```

*Численное наблюдение*: при $t=2$ вес $p_2(2)=w_0(2)$ сектора $N=2$ («без потерь») равен $0.230$ (напечатан
выше). Без перескоков ($J=0$) он был бы ровно $e^{-2\gamma_1t}=e^{-2.4}=0.091$, поскольку тогда $|20\rangle$ — собственный вектор
$H_{\mathrm{eff}}$ с собственным значением $U-i\gamma_1$ и $w_0(t)=|e^{-it(U-i\gamma_1)}|^2$ (предложение 8.10); так что
здесь веса секторов действительно зависят от гамильтониана. Постселектированное состояние находится на следовом расстоянии
$0.135$ (напечатано выше) от состояния замкнутой системы: обусловливание на $N=2$ меняет состояние сверх унитарной эволюции.

<a id="s-8-5"></a>

## 8.5 Стационарные состояния и релаксация

*Конспект:* [раздел 8.5](../../lecture_notes/ru/08_dynamics.md#s-8-5); [следствие 8.13](../../lecture_notes/ru/08_dynamics.md#r-8-13), [предложение 8.14](../../lecture_notes/ru/08_dynamics.md#r-8-14) (спонтанное излучение).

По [предложению 8.14](../../lecture_notes/ru/08_dynamics.md#r-8-14) для спонтанного излучения на $M_2$ ($H=\omega|e\rangle\langle e|$, один скачок
$|g\rangle\langle e|$ со скоростью $\gamma$) имеем $e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}}$, спектр
$\{0,-\gamma,-\gamma/2\pm i\omega\}$ и спектральную щель $\gamma/2$. Проверим спектр
с помощью `lindblad_superoperator`, канал и следовое расстояние — с помощью `lindblad_channel`, а
[следствие 8.13](../../lecture_notes/ru/08_dynamics.md#r-8-13) — для случайных скачков внутри секторов $M_2\oplus M_3$:

```python
q2 = TypeIAlgebra([2], [2], precision='double')                    # basis |g>, |e>
omega, gamma = 1.3, 0.4
H_at = q2.from_blocks([torch.tensor([[0.0, 0.0], [0.0, omega]], dtype=f64)])
sigma_minus = q2.from_blocks([torch.tensor([[0.0, 1.0], [0.0, 0.0]], dtype=f64)])   # |g><e|
gen = dynamics.lindblad_superoperator(q2, H_at, [sigma_minus], rates=[gamma])[0][0]
spec = sorted(torch.linalg.eigvals(gen).tolist(), key=lambda z: (z.real, z.imag))
print([f"{z.real:+.2f}{z.imag:+.2f}i" for z in spec])
assert np.allclose(spec, [-gamma, -gamma / 2 - 1j * omega, -gamma / 2 + 1j * omega, 0.0], atol=1e-12)

psi_at = torch.tensor([0.6, 0.8], dtype=c128)                      # rho_ee(0) = 0.64, rho_ge(0) = 0.48
rho_at = DensityMatrix(q2, matrix=torch.outer(psi_at, psi_at.conj())[None, None])
ground = DensityMatrix(q2, matrix=torch.tensor([[[[1.0, 0.0], [0.0, 0.0]]]], dtype=c128))
for s in (0.5, 2.0, 10.0):
    ch = dynamics.lindblad_channel(q2, H_at, [sigma_minus], t=s, rates=[gamma])
    UA = channels.unitary_channel(dynamics.propagator(H_at, s)).compose(
        channels.amplitude_damping_channel(q2, 1 - math.exp(-gamma * s)))
    assert (ch.superoperator()[0] - UA.superoperator()[0]).abs().max() < 1e-12
    r = ch(rho_at).matrix[0, 0]
    assert abs(r[1, 1].real - 0.64 * math.exp(-gamma * s)) < 1e-12
    assert abs(r[0, 1] - 0.48 * np.exp((1j * omega - gamma / 2) * s)) < 1e-12
    exact = math.sqrt(0.64 ** 2 * math.exp(-2 * gamma * s) + 0.48 ** 2 * math.exp(-gamma * s))
    assert abs(states.trace_distance(ch(rho_at), ground).item() - exact) < 1e-12

two = TypeIAlgebra([2, 3], [2, 3], precision='double')
H_two = two.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
L_two = [two.operator(torch.randn(1, 2, 3, 3, dtype=c128)) for _ in range(2)]
kernel = [int((torch.linalg.eigvals(G_c[0]).abs() < 1e-10).sum())
          for G_c in dynamics.lindblad_superoperator(two, H_two, L_two)]
print("stationary states per sector:", kernel)
assert all(dim >= 1 for dim in kernel)                             # Corollary 8.13
```

*Численное наблюдение*: для этих случайных скачков код находит ровно одно стационарное состояние в каждом секторе.
По теореме D.6, применённой посекторно ([конспект лекций, раздел 8.5](../../lecture_notes/ru/08_dynamics.md#s-8-5)), так и должно быть, если в каждом секторе
есть точное стационарное состояние и единственные матрицы, коммутирующие с $H_c$ и со всеми $L_{j,c}$, $L_{j,c}^*$, — это
кратные $1$; код проверяет оба эти предположения (в арифметике с плавающей точкой, для использованных выше скоростей $g_j=1$),
а затем проверяет вывод.

```python
for c, k in enumerate(two.k_factors):                              # Theorem D.6 in every sector
    G_c = dynamics.lindblad_superoperator(two, H_two, L_two)[c][0].numpy()   # on row-major vec(rho_c)
    w_c, V_c = np.linalg.eig(G_c)
    sigma_c = V_c[:, np.argmin(abs(w_c))].reshape(k, k)
    sigma_c = sigma_c / np.trace(sigma_c)
    assert np.abs(G_c @ sigma_c.reshape(-1)).max() < 1e-12            # stationary
    assert np.abs(sigma_c - sigma_c.conj().T).max() < 1e-10 and np.linalg.eigvalsh(sigma_c).min() > 1e-3   # faithful
    ops = [H_two.matrix[0, c, :k, :k].numpy()] + [L.matrix[0, c, :k, :k].numpy() for L in L_two]
    ops += [A.conj().T for A in ops[1:]]
    comm = np.concatenate([np.kron(A, np.eye(k)) - np.kron(np.eye(k), A.T) for A in ops])   # X -> [A, X]
    assert int((np.linalg.svd(comm, compute_uv=False) < 1e-9).sum()) == 1   # commutant = C 1
assert kernel == [1, 1]                                            # the conclusion of Theorem D.6
```

<a id="s-8-6"></a>

## 8.6 С чем сравнивать

QuTiP (`sesolve`, `mesolve` для уравнения Линдблада с адаптивными решателями ОДУ, `mcsolve` для
квантовых траекторий, `steadystate`, `liouvillian`), dynamiqs (JAX; пакетные и дифференцируемые
решатели уравнений Шрёдингера и Линдблада на GPU), QuantumOptics.jl и SciPy (`scipy.linalg.expm`, а также
`scipy.sparse.linalg.expm_multiply` для $e^{t\mathcal L}\rho$ с разреженным векторизованным генератором) решают
те же уравнения на одном гильбертовом пространстве, и для задач выше они должны давать те же числа.
На $\mathcal M$ им передают блочно-диагональные матрицы $D\times D$; векторизованный генератор тогда
действует на пространстве размерности $D^2$ (разреженный, но с секторами, связанными в одном решении), если только пользователь
не разбивает задачу на секторы вручную. Библиотека встраивает блочную структуру: точные пропагаторы
и каналы Линдблада сектор за сектором (супероператоры $k_c^2\times k_c^2$), скачки между секторами в виде
`InterSectorChannel` со сжатием из предложения 8.8, пакеты гамильтонианов и операторов
скачков (для усреднения по беспорядку) и автоматическое дифференцирование (последние два — также в dynamiqs). В ней
нет адаптивного выбора шага (RK4 работает с фиксированным шагом и лишь приближённо сохраняет положительность), нет
решателя для траекторий, и `lindblad_channel` не принимает скачков между секторами; для них эталоном служат QuTiP или dynamiqs на
$M_D$.

<a id="x-8"></a>

## Упражнения

[Упражнения 2 и 4](../../lecture_notes/ru/08_dynamics.md#x-8) конспекта лекций: собственные значения
$-\frac i4(\gamma_1+\gamma_2)\pm\sqrt{J^2-\delta^2}$, $\delta=(\gamma_1-\gamma_2)/4$, оператора $H_{\mathrm{eff}}$ в секторе $N=1$ из
раздела 8.4, с единственной жордановой клеткой размера $2$ при $J=|\delta|$; и спектр $\{0,0,-2\gamma\pm i\omega\}$
чистой дефазировки на $M_2$ ($L=\sigma_z$). Проверьте то и другое численно.

<details>
<summary>Решения</summary>

```python
lam = np.linalg.eigvals(H_eff[1])
delta = (g1 - g2) / 4
root = np.sqrt(J ** 2 - delta ** 2 + 0j)
assert np.allclose(sorted(lam, key=np.real), sorted(-0.25j * (g1 + g2) + np.array([root, -root]), key=np.real))
H_ep = np.array([[-0.5j * g1, delta], [delta, -0.5j * g2]])         # J = delta: the exceptional point
Y = H_ep + 0.25j * (g1 + g2) * np.eye(2)                           # H_ep minus its only eigenvalue
assert np.abs(Y @ Y).max() < 1e-15 and np.linalg.matrix_rank(Y) == 1   # one Jordan block of size 2
```

```python
deph = dynamics.lindblad_superoperator(q2, H_at, [q2.from_blocks([torch.diag(torch.tensor([1.0, -1.0], dtype=f64))])],
                                       rates=[gamma])[0][0]
ev = torch.linalg.eigvals(deph)
assert int((ev.abs() < 1e-12).sum()) == 2
assert np.allclose(sorted(ev.tolist(), key=lambda z: (z.real, z.imag))[:2],
                   [-2 * gamma - 1j * omega, -2 * gamma + 1j * omega], atol=1e-12)
```

</details>

**Далее:** [Руководство 9. Спиновые цепочки](09_spin_chains.md).
