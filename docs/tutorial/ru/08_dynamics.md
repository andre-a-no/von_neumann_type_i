[English](../08_dynamics.md) · [Español](../es/08_dynamics.md) · [Français](../fr/08_dynamics.md) · [Deutsch](../de/08_dynamics.md) · [中文](../zh/08_dynamics.md) · [日本語](../ja/08_dynamics.md) · **Русский**

# 8. Динамика

Замкнутая квантовая система эволюционирует под действием однопараметрической группы унитарных операторов, открытая — система,
связанная с ненаблюдаемым окружением, — под действием полугруппы каналов. В алгебре с
секторами у обеих эволюций есть особенность, невидимая на одной матричной алгебре: гамильтониан из
$\mathcal M$ никогда не может изменить сектор, поэтому заряды сохраняются, тогда как открытая система может
перескакивать между секторами, например теряя частицу. В этой главе доказываются основные факты
об обоих видах эволюции, выписывается уравнение Линдблада со скачками между секторами и показывается, что
«ни одна частица не потеряна» — это измерение заряда, и постселектированное состояние эволюционирует под действием
неэрмитова гамильтониана.

**Вы узнаете**, почему веса секторов являются интегралами движения, почему генераторы
формы Горини–Коссаковского–Сударшана–Линдблада (GKSL) порождают каналы, как скачки между секторами входят
в уравнение Линдблада и когда веса секторов подчиняются классическому кинетическому уравнению, что означают
эволюция без скачков и постселекция и как находить стационарные состояния и скорости релаксации.

Всюду далее $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ действует на
$\mathcal H=\bigoplus_c\mathbb C^{k_c}$, $D=\sum_ck_c$, $P_c$ — минимальные центральные проекторы,
состояния задаются своими плотностями $\rho$ относительно `Tr_blunt` с весами секторов $p_c=\operatorname{Tr}\rho_c$
(глава 3), $[A,B]=AB-BA$, $\{A,B\}=AB+BA$ и $\hbar=1$. **Супероператор** — это линейное отображение
на $\mathcal M$ или на $M_D(\mathbb C)$; все нормы супероператоров — операторные нормы, а потому
субмультипликативны. Для линейного отображения $B$ конечномерного пространства $e^{tB}=\sum_{n\ge0}t^nB^n/n!$.

## 8.1 Замкнутые системы

**Гамильтониан** — это самосопряжённый элемент $H\in\mathcal M$. Вектор эволюционирует по **уравнению
Шрёдингера** $i\psi'(t)=H\psi(t)$, плотность — по **уравнению фон Неймана** $\rho'(t)=-i[H,\rho(t)]$.
Оба уравнения линейны, с постоянными коэффициентами, на конечномерном пространстве, и
следующая лемма — это вся теория дифференциальных уравнений, которая нам нужна.

**Лемма 8.1.** Пусть $B$ — линейное отображение конечномерного нормированного пространства $X$. Тогда
$\frac d{dt}e^{tB}=Be^{tB}=e^{tB}B$, $e^{(s+t)B}=e^{sB}e^{tB}$, $\|e^{tB}\|\le e^{|t|\,\|B\|}$, и для
любого $x_0\in X$ функция $x(t)=e^{tB}x_0$ — единственная дифференцируемая функция $x$ на промежутке
$I\ni0$ с $x'=Bx$ и $x(0)=x_0$.

*Доказательство.* Ряд сходится абсолютно, с указанной оценкой нормы, и как степенной ряд по $t$ его
можно дифференцировать почленно, что даёт первую формулу. Далее,
$\frac d{dt}(e^{tB}e^{-tB})=Be^{tB}e^{-tB}-e^{tB}Be^{-tB}=0$, поэтому $e^{tB}e^{-tB}=1$. Если $x'=Bx$, то
$\frac d{dt}(e^{-tB}x(t))=e^{-tB}(x'-Bx)=0$, так что $e^{-tB}x(t)=x_0$ и $x(t)=e^{tB}x_0$. Наконец,
$t\mapsto e^{(s+t)B}x_0$ и $t\mapsto e^{tB}e^{sB}x_0$ — решения $y'=By$ с одним и тем же начальным значением. $\square$

**Предложение 8.2.** Пусть $H=H^*\in\mathcal M$ и $U(t)=e^{-itH}$, определённый с помощью функционального
исчисления из раздела 1.3 с $f_t(\lambda)=e^{-it\lambda}$.

1. $U(t)=\bigoplus_ce^{-itH_c}$ — унитарный элемент $\mathcal M$, $U(s+t)=U(s)U(t)$ и $U'(t)=-iHU(t)$;
   $U(t)$ равен сумме экспоненциального ряда для $-itH$.
2. $\psi(t)=U(t)\psi_0$ — единственное решение $i\psi'=H\psi$ с $\psi(0)=\psi_0$, а
   $\rho(t)=U(t)\rho_0U(t)^*$ — единственное решение в $\mathcal M$ уравнения $\rho'=-i[H,\rho]$ с
   $\rho(0)=\rho_0$. Если $\rho_0$ — состояние, то и каждое $\rho(t)$ — состояние.
3. (*Картина Гейзенберга*) $\operatorname{Tr}(\rho(t)A)=\operatorname{Tr}(\rho_0A(t))$, где
   $A(t)=U(t)^*AU(t)$, и $A'=i[H,A]$.

*Доказательство.* (1) Диагонализуем каждый блок: $H_c=V_c\Lambda_cV_c^*$; тогда $V=\bigoplus_cV_c$ диагонализует
$H$, и $U(t)=Ve^{-it\Lambda}V^*$ блочно-диагонален (по теореме 1.5(a) результат не зависит от
диагонализующей унитарной матрицы). Поскольку $|f_t|^2=1$ и $f_sf_t=f_{s+t}$, теорема 1.5(c) даёт унитарность и
групповой закон. Поэлементное дифференцирование $e^{-it\Lambda}$ даёт $U'=V(-i\Lambda)e^{-it\Lambda}V^*=-iHU$, так что
$U(t)$ — это ряд $e^{-itH}$ по лемме 8.1 ($X=M_D$, $B$ — умножение слева на $-iH$).
(2) $(U\psi_0)'=-iHU\psi_0$ и $(U\rho_0U^*)'=-iHU\rho_0U^*+U\rho_0U^*(iH)$, поскольку $(U')^*=iU^*H$.
Единственность — это лемма 8.1 с $B=-iH$ на $\mathcal H$ и соответственно $B=-i[H,\cdot\,]$ на $\mathcal M$. Положительность:
$\langle x,U\rho_0U^*x\rangle=\langle U^*x,\rho_0U^*x\rangle\ge0$; след — по цикличности. (3) Цикличность
следа и $(U^*AU)'=iU^*HAU-iU^*AHU=i[H,A(t)]$. $\square$

**Теорема 8.3 (сохранение заряда).** Пусть $t\mapsto H(t)$ принимает значения в множестве самосопряжённых
элементов $\mathcal M$ на промежутке $I$.

1. Если $\psi:I\to\mathcal H$ дифференцируема и $i\psi'=H(t)\psi$, то каждая норма $\|P_c\psi(t)\|$ постоянна.
2. Если $\rho:I\to\mathcal M$ дифференцируема и $\rho'=-i[H(t),\rho]$, то $\operatorname{Tr}(Z\rho(t))$
   постоянен для каждого центрального $Z$; в частности, постоянен каждый вес сектора $p_c(t)=\operatorname{Tr}(P_c\rho(t))$.
3. Если $H$ постоянен и $A\in\mathcal M$ коммутирует с $H$, то $\operatorname{Tr}(A\rho(t))$ постоянен.

*Доказательство.* (1) Скалярное произведение антилинейно по первому аргументу, поэтому
$\frac d{dt}\langle\psi,P_c\psi\rangle=\langle-iH\psi,P_c\psi\rangle+\langle\psi,-iP_cH\psi\rangle
=i\langle\psi,HP_c\psi\rangle-i\langle\psi,P_cH\psi\rangle=0$, так как $P_c$ централен. (2)
$\frac d{dt}\operatorname{Tr}(Z\rho)=-i\operatorname{Tr}(Z[H,\rho])=-i\operatorname{Tr}([Z,H]\rho)=0$ в силу
цикличности. (3) То же вычисление с $A$. $\square$

> **Физический смысл.** Гамильтониан из $\mathcal M$ коммутирует с зарядами (числом частиц,
> намагниченностью, квантовым числом симметрии), поэтому вероятность каждого значения заряда — интеграл
> движения, даже для зависящих от времени гамильтонианов, описывающих внешнее воздействие. С вычислительной точки зрения эволюция распадается на
> $C$ независимых задач: одна диагонализация стоит порядка $\sum_ck_c^3$ операций вместо $D^3$.

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

На $\mathcal M$ вектор $\psi=\bigoplus_c\psi_c$ задаёт состояние смеси
$\bigoplus_c|\psi_c\rangle\langle\psi_c|$ (глава 3, упражнение 3), поэтому `dynamics.schrodinger` (векторы в
блочной раскладке `(batch, C, k_max)`) и `von_neumann` должны на нём совпадать. Для зависящего от времени $H(t)$
библиотека использует классический метод Рунге–Кутты с фиксированным шагом (`dynamics.rk4`, `schrodinger_rk4`): шаг умножает
$\psi$ на многочлен от значений $H(t)\in\mathcal M$, поэтому он никогда не переносит амплитуду между секторами,
но он не является точно унитарным, и веса секторов дрейфуют на уровне ошибки дискретизации.
(Для одного большого сектора `krylov.evolve` нужны только произведения матрицы на вектор; см. главу 9.)

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

Веса секторов меняются на $1.1\cdot10^{-9}$ при RK4 и не меняются сверх ошибок округления при
точных методах решения.

## 8.2 Квантовые динамические полугруппы

CP- и TP-отображения и каналы определены в разделе 6.1 ([глава 6](06_quantum_channels.md)). Помимо
теоремы 6.2 и предложения 6.4 (отображение $A\mapsto\sum_lK_lAK_l^*$, в том числе между секторами, является CP) нам нужны
три свойства замкнутости.

**Лемма 8.4.** (a) Суммы, неотрицательные кратные и композиции CP-отображений являются CP. Если
$\Phi$ является CP на $M_D$ и переводит $*$-подалгебру $\mathcal N\subseteq M_D$ (например, $\mathcal M$ или
$P_c\mathcal MP_c\cong M_{k_c}$ с единицей $P_c$) в себя, то $\Phi|_{\mathcal N}$ является CP.
(b) Если CP-отображения $\Phi_m$ сходятся к $\Phi$ (т. е. $\Phi_m(A)\to\Phi(A)$ для всех $A$), то $\Phi$ является CP.
(c) Если $s\mapsto\Psi(s)$ — непрерывное отображение $[a,b]$ в множество CP-отображений, то $\int_a^b\Psi(s)\,ds$ является CP.

*Доказательство.* (a) Отображение $\Phi\mapsto\Phi_n$ линейно и $(\Phi\Psi)_n=\Phi_n\Psi_n$; суммы положительных матриц
положительны. Для ограничения: $M_n(\mathcal N)\subseteq M_n(M_D)$, и положительность в обеих алгебрах — одно и то же
условие (после перестановки базиса элемент $M_n(\mathcal N)$ является прямой суммой
блоков $[(A_{ij})_c]$ и нулей). (b) Положительные матрицы образуют замкнутое множество ($\langle x,Xx\rangle\ge0$ при каждом
$x$ — замкнутое условие), и $(\Phi_m)_n(X)\to\Phi_n(X)$. (c) Суммы Римана являются CP по (a) и
сходятся к интегралу, поскольку $\Psi$ равномерно непрерывно; применяем (b). $\square$

**Определение.** **Квантовая динамическая полугруппа** на $\mathcal M$ — это семейство $(T_t)_{t\ge0}$
каналов $\mathcal M\to\mathcal M$ с $T_0=\mathrm{id}$, $T_{s+t}=T_sT_t$, для которого $t\mapsto T_t$ непрерывно.

**Лемма 8.5.** Всякая непрерывная полугруппа супероператоров ($T_0=\mathrm{id}$, $T_{s+t}=T_sT_t$) имеет вид $T_t=e^{t\mathcal L}$
ровно для одного супероператора $\mathcal L$, её **генератора**, и $\mathcal L=\lim_{t\to0}(T_t-\mathrm{id})/t$.

*Доказательство.* $V_h=\frac1h\int_0^hT_s\,ds\to\mathrm{id}$ при $h\to0$, поэтому $V_h$ обратим при некотором фиксированном
малом $h$. Тогда $T_tV_h=\frac1h\int_t^{t+h}T_s\,ds$ дифференцируемо по $t$ с производной
$\frac1h(T_{t+h}-T_t)=T_t\frac{T_h-\mathrm{id}}h$. Следовательно, $T_t$ дифференцируемо и $T_t'=T_t\mathcal L$
с $\mathcal L=\frac{T_h-\mathrm{id}}hV_h^{-1}$, и $T_t=e^{t\mathcal L}$ по лемме 8.1 (где $B$ — композиция
справа с $\mathcal L$). Единственность: $\mathcal L=T'(0)$. $\square$

**Теорема 8.6 (генераторы GKSL).** Пусть $H=H^*$ и $L_1,\dots,L_J$ лежат в $M_D$, и
$$
\mathcal L(\rho)=-i[H,\rho]+\sum_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
$$
Тогда $(e^{t\mathcal L})_{t\ge0}$ — квантовая динамическая полугруппа на $M_D$. Точнее, положим
$H_{\mathrm{eff}}=H-\frac i2\sum_jL_j^*L_j$, $\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)$
и $\Phi(\rho)=\sum_jL_j\rho L_j^*$, так что $\mathcal L=\mathcal L_0+\Phi$, и определим
$$
S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},\qquad
S_n(t)=\int_0^tS_0(t-s)\,\Phi\,S_{n-1}(s)\,ds\quad(n\ge1).
$$
Тогда каждое $S_n(t)$ является CP, и $e^{t\mathcal L}=\sum_{n\ge0}S_n(t)$ (**разложение Дайсона**), причём
$\|S_n(t)\|\le e^{t\|\mathcal L_0\|}(t\|\Phi\|)^n/n!$.

*Доказательство.* *Шаг 1.* $X(t)=e^{-itH_{\mathrm{eff}}}\rho e^{itH_{\mathrm{eff}}^*}$ удовлетворяет $X'=\mathcal L_0(X)$
и $X(0)=\rho$, поэтому $S_0(t)=e^{t\mathcal L_0}$ по лемме 8.1; это CP-отображение с единственным оператором Крауса
$e^{-itH_{\mathrm{eff}}}$. *Шаг 2.* Обозначим $a=\|\mathcal L_0\|$, $b=\|\Phi\|$. По индукции $S_n$
непрерывны и
$\|S_n(t)\|\le\int_0^te^{(t-s)a}\,b\,e^{sa}\frac{(bs)^{n-1}}{(n-1)!}ds=e^{ta}\frac{(bt)^n}{n!}$, так что
ряд сходится равномерно на каждом $[0,T]$. *Шаг 3.* Поскольку
$S_n(t)=S_0(t)\int_0^te^{-s\mathcal L_0}\Phi S_{n-1}(s)\,ds$, правило дифференцирования произведения и основная
теорема анализа дают $S_n'=\mathcal L_0S_n+\Phi S_{n-1}$ ($n\ge1$) и $S_0'=\mathcal L_0S_0$. В силу
оценок ряд из производных сходится равномерно на $[0,T]$, поэтому $Y(t)=\sum_nS_n(t)$
дифференцируемо, $Y'=\mathcal L_0Y+\Phi Y=\mathcal LY$, $Y(0)=\mathrm{id}$, и $Y(t)=e^{t\mathcal L}$ по
лемме 8.1. *Шаг 4.* $\Phi$ является CP (теорема 6.2), поэтому по индукции подынтегральное выражение
$s\mapsto S_0(t-s)\Phi S_{n-1}(s)$ — непрерывное семейство CP-отображений, и $S_n(t)$ является CP (лемма 8.4);
таковы же частичные суммы и их предел $e^{t\mathcal L}$. *Шаг 5.* В силу цикличности
$\operatorname{Tr}\mathcal L(\rho)=\sum_j\big(\operatorname{Tr}(L_j\rho L_j^*)-\operatorname{Tr}(L_j^*L_j\rho)\big)=0$,
следовательно, $\frac d{dt}\operatorname{Tr}e^{t\mathcal L}\rho=\operatorname{Tr}\mathcal L(e^{t\mathcal L}\rho)=0$: каждое
$e^{t\mathcal L}$ является TP. Полугрупповой закон и непрерывность — это лемма 8.1. $\square$

**Теорема 8.7 (Горини–Коссаковский–Сударшан; Линдблад; приводится без доказательства).** Обратно, генератор любой квантовой
динамической полугруппы на $M_D$ имеет вид из теоремы 8.6. [V. Gorini, A. Kossakowski,
E. C. G. Sudarshan, *J. Math. Phys.* 17, 821 (1976); G. Lindblad, *Commun. Math. Phys.* 48, 119
(1976).] Для непрерывной по норме полугруппы унитальных CP-отображений на произвольной алгебре фон Неймана
$\mathcal M$ Кристенсен и Эванс [*J. London Math. Soc.* 20, 358 (1979)] показали, что генератор
имеет вид $A\mapsto\Psi(A)+K^*A+AK$, где $\Psi$ является CP и $K\in\mathcal M$.

На $\mathcal M=\bigoplus_cM_{k_c}$ применим это к полугруппе $(T_t^*)$ в картине Гейзенберга: её
отображения являются CP и унитальны (предложение 6.8, следствие 6.9(b)), а её генератор равен $\mathcal L^*$.
Переходя обратно к сопряжённым относительно `Tr_blunt`, получаем $\mathcal L(\rho)=\Phi(\rho)+G\rho+\rho G^*$, где $\Phi=\Psi^*$ является CP на
$\mathcal M$ (операторы Крауса между секторами, предложение 6.4) и $G=K$. Сохранение следа, $\operatorname{Tr}\big((\Phi^*(1)+G+G^*)\rho\big)=0$ для всех $\rho\in\mathcal M$, влечёт
$G+G^*=-\Phi^*(1)$, т. е. $G=-iH-\frac12\Phi^*(1)$ с $H=H^*\in\mathcal M$: это в точности генераторы
следующего раздела.

> **Физический смысл.** Система, слабо взаимодействующая с большим окружением, корреляции в котором
> быстро затухают (приближение Борна–Маркова), теряет память о прошлом: её состояние в момент $s+t$
> получается из состояния в момент $s$ каналом, зависящим только от $t$. Теорема 8.7 утверждает, что
> такая динамика всегда есть гамильтонова часть плюс «скачки» $L_j$ (испускание фотона, потеря частицы,
> дефазировка) со скоростями, задаваемыми нормировкой $L_j$.

В библиотеке четыре процедуры для уравнения Линдблада: `lindblad_rhs` возвращает правую часть в виде функции для
`rk4`; `lindblad_evolve` интегрирует методом RK4; `lindblad_superoperator` возвращает для каждого сектора
матрицу $k_c^2\times k_c^2$ отображения $\mathcal L$ на векторах $\operatorname{vec}\rho_c$, развёрнутых по строкам; а
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

## 8.3 Уравнения Линдблада в алгебре с секторами

На $\mathcal M$ бывают скачки двух видов. Скачок **внутри секторов** — это оператор
$L_j\in\mathcal M$. Скачок **между секторами** — это CP-отображение $\Phi:\mathcal M\to\mathcal M$ с операторами
Крауса $K^{dc}_l:\mathbb C^{k_c}\to\mathbb C^{k_d}$ (предложение 6.4, объект `InterSectorChannel`):
$\Phi(\rho)_d=\sum_{c,l}K^{dc}_l\rho_cK^{dc*}_l$. По предложению 6.8 его сопряжённое относительно `Tr_blunt` удовлетворяет
$\Phi^*(1)_c=\sum_{d,l}K^{dc*}_lK^{dc}_l$ — блочно-диагональный элемент $\mathcal M$. Со скоростями
$g_j,g_m>0$ генератор `lindblad_rhs` равен
$$
\mathcal L(\rho)=-i[H,\rho]+\sum_jg_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big)
+\sum_mg_m\Big(\Phi_m(\rho)-\tfrac12\{\Phi_m^*(1),\rho\}\Big).
$$

**Предложение 8.8.** Пусть $J^{dc}_l\in M_D$ — матрица с блоком $K^{dc}_l$ в позиции $(d,c)$
и нулями в остальных местах.

1. $\mathcal L$ переводит $\mathcal M$ в $\mathcal M$ и является ограничением на $\mathcal M$ генератора GKSL
   $\hat{\mathcal L}$ на $M_D$ со скачками $\sqrt{g_j}L_j$ и $\sqrt{g_m}J^{dc}_{m,l}$ (при всех $m,l,d,c$). Следовательно,
   $(e^{t\mathcal L})_{t\ge0}$ — квантовая динамическая полугруппа на $\mathcal M$.
2. Пусть теперь $\tilde{\mathcal L}$ — генератор GKSL на $M_D$, в котором блоки с одним и тем же индексом
   Крауса просуммированы в один скачок, $L_{m,l}=\sum_{d,c}J^{dc}_{m,l}$. Тогда
   $\mathcal L(\rho)=E\big(\tilde{\mathcal L}(\rho)\big)$ при $\rho\in\mathcal M$, где $E(X)=\sum_dP_dXP_d$ —
   условное ожидание $M_D$ на $\mathcal M$; вообще говоря, $\tilde{\mathcal L}(\rho)\notin\mathcal M$.

*Доказательство.* Поскольку $\rho=\sum_cP_c\rho P_c$ и $J^{dc}P_{c'}=\delta_{cc'}J^{dc}$,
$J^{dc}\rho J^{dc*}$ — это блок $K^{dc}\rho_cK^{dc*}$ в позиции $(d,d)$, а $J^{dc*}J^{dc}$ — это
$K^{dc*}K^{dc}$ в позиции $(c,c)$. Суммирование по $d,l$ даёт $\Phi(\rho)$ и $\Phi^*(1)$, что доказывает
тождество из (1), и все слагаемые лежат в $\mathcal M$. Значит, $\hat{\mathcal L}$ и все его степени переводят $\mathcal M$
в $\mathcal M$, $e^{t\hat{\mathcal L}}|_{\mathcal M}=e^{t\mathcal L}$, и это канал на $\mathcal M$ по
теореме 8.6 и лемме 8.4(a). (2) $L\rho L^*=\sum_{d,d',c}J^{dc}\rho J^{d'c*}$ для $L=\sum_{d,c}J^{dc}$, и
$E$ сохраняет слагаемые с $d=d'$. Для $B\in\mathcal M$ имеем $P_dXBP_d=P_dXP_dB$, поэтому $E(XB)=E(X)B$ и
$E(BX)=BE(X)$; значит, $E\{L^*L,\rho\}=\{E(L^*L),\rho\}$, и $E(L^*L)=\sum_cP_cL^*LP_c=\sum_{c,d}J^{dc*}J^{dc}$,
поскольку $J^{dc*}J^{d'c}=0$ при $d\ne d'$. Наконец, $E(-i[H,\rho])=-i[H,\rho]$. $\square$

> **Физический смысл.** Часть (2) говорит, что уравнение Линдблада на $\mathcal M$ — это то, что остаётся от
> уравнения на $M_D$, если когерентности между секторами, создаваемые скачками, сразу стираются, —
> в этом и состоит содержание правила суперотбора. При потере частиц ничего не стирается: оператор
> уничтожения $a_j$ понижает число частиц ровно на единицу, поэтому $a_j\rho a_j^*$ и $a_j^*a_j=n_j$
> блочно-диагональны, и уравнение на $\mathcal M$ — обычное.

**Предложение 8.9 (веса секторов).** Пусть $\rho(t)=e^{t\mathcal L}\rho_0$ и $p(t)$ — его веса секторов.

1. Если скачков между секторами нет, то каждое $p_c(t)$ постоянно.
2. Если для каждого $m$ и каждой пары $d\ne c$ выполнено $\sum_lK^{dc*}_{m,l}K^{dc}_{m,l}=r_{m,dc}1_{k_c}$ с числами
   $r_{m,dc}$, то $p'=Qp$, где $Q_{dc}=\sum_mg_mr_{m,dc}$ при $d\ne c$ и $Q_{cc}=-\sum_{d\ne c}Q_{dc}$.
   Следовательно, $p(t)=e^{tQ}p(0)$ не зависит ни от $H$, ни от $L_j$, ни от состояний внутри секторов.

*Доказательство.* $p_d'=\operatorname{Tr}(P_d\mathcal L(\rho))$. Гамильтоново слагаемое обращается в нуль, как в теореме 8.3.
Для $L_j\in\mathcal M$ имеем $\operatorname{Tr}(P_dL_j\rho L_j^*)-\operatorname{Tr}(P_dL_j^*L_j\rho)
=\operatorname{Tr}(L_j^*P_dL_j\rho)-\operatorname{Tr}(L_j^*L_jP_d\rho)=0$, так как $P_d$ централен. Для скачка
$\Phi$ между секторами $\operatorname{Tr}(P_d\Phi(\rho))=\sum_c\operatorname{Tr}\big(\sum_lK^{dc*}_lK^{dc}_l\rho_c\big)$
и $\operatorname{Tr}(P_d\tfrac12\{\Phi^*(1),\rho\})=\sum_{d'}\operatorname{Tr}\big(\sum_lK^{d'd*}_lK^{d'd}_l\rho_d\big)$;
слагаемые с $c=d$ и $d'=d$ сокращаются. Это доказывает (1), а при условии (2) остающиеся слагаемые
равны $\sum_{c\ne d}r_{dc}p_c-\sum_{d'\ne d}r_{d'd}p_d$. Лемма 8.1 даёт $p(t)=e^{tQ}p(0)$. $\square$

При условии (2) для $d\ne c$ число $r_{m,dc}$ совпадает с элементом $T_{dc}(\rho)$ матрицы
переходов $\Phi_m$ (определённой перед предложением 6.12) для любого состояния $\rho$. Поэтому для одного скачка со скоростью
$1$ имеем $Q=R-\operatorname{diag}(\sum_dR_{dc})_c$, где $R$ = `transition_matrix()` (диагональ $R$ сокращается).

**Пример: потеря частиц.** Две бозонные моды, не более двух частиц; секторы —
$N=0,1,2$ с базисами $|00\rangle$; $|10\rangle,|01\rangle$; $|20\rangle,|11\rangle,|02\rangle$, так что
$\mathcal M=M_1\oplus M_2\oplus M_3$. Димер Бозе–Хаббарда
$H=J(a_1^*a_2+a_2^*a_1)+\frac U2\sum_jn_j(n_j-1)$ коммутирует с $N$, и каждая мода теряет частицы со
скоростью $\gamma$: это один скачок между секторами с операторами Крауса $a_1,a_2$, ограниченными на $N\to N-1$.
Поскольку $a_1^*a_1+a_2^*a_2=N$, применимо предложение 8.9(2) с $Q_{N-1,N}=\gamma N$, и при старте из $N=2$
$$
p_2(t)=q^2,\qquad p_1(t)=2q(1-q),\qquad p_0(t)=(1-q)^2,\qquad q=e^{-\gamma t},
$$
— биномиальный закон: каждая частица независимо выживает с вероятностью $e^{-\gamma t}$, что бы ни делали
перескоки и взаимодействие (проверьте $p'=Qp$ дифференцированием).

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
```

RK4 не является точно вполне положительным: собственные значения $\rho(t)$ могут быть отрицательными на уровне
ошибки дискретизации (здесь $\approx-2\cdot10^{-8}$ при шагах не больше $0.02$); `lindblad_channel` точен,
но работает только со скачками внутри секторов. Наконец, проверим предложение 8.8(2) для случайного отображения между
секторами, сравнивая с генератором GKSL на $M_6$, построенным средствами numpy:

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
```

## 8.4 Постселекция и эволюция без скачков

Соберём скачки генератора из раздела 8.3 в $\Phi_{\mathrm{tot}}(\rho)=\sum_jg_jL_j\rho L_j^*+\sum_mg_m\Phi_m(\rho)$,
и пусть $A=\Phi_{\mathrm{tot}}^*(1)=\sum_jg_jL_j^*L_j+\sum_mg_m\Phi_m^*(1)\ge0$ в $\mathcal M$ и
$H_{\mathrm{eff}}=H-\frac i2A$. По предложению 8.8(1) и теореме 8.6 (её $S_n$ переводят $\mathcal M$ в
$\mathcal M$, поскольку $H_{\mathrm{eff}}\in\mathcal M$)
$$
\rho(t)=\sum_{n\ge0}S_n(t)\rho,\qquad S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},
\qquad S_n(t)=\int_0^tS_0(t-s)\,\Phi_{\mathrm{tot}}\,S_{n-1}(s)\,ds .
$$

**Предложение 8.10.** Для состояния $\rho$ положим $w_n(t)=\operatorname{Tr}S_n(t)\rho$. Тогда $w_n(t)\ge0$,
$\sum_nw_n(t)=1$, и $w_0$ не возрастает, причём $w_0'(t)=-\operatorname{Tr}\big(A\,S_0(t)\rho\big)$. Для
чистого состояния $\rho=|\psi\rangle\langle\psi|$ имеем $w_0(t)=\|e^{-itH_{\mathrm{eff}}}\psi\|^2$.

*Доказательство.* $S_n(t)$ является CP, а значит, положительно, и $\sum_nw_n=\operatorname{Tr}\rho(t)=1$. При
$X=S_0(t)\rho\ge0$ имеем $w_0'=\operatorname{Tr}(-iH_{\mathrm{eff}}X+iXH_{\mathrm{eff}}^*)=-\operatorname{Tr}(AX)$,
и $\operatorname{Tr}(AX)=\operatorname{Tr}(X^{1/2}AX^{1/2})\ge0$. Последнее утверждение — это
$\operatorname{Tr}|\phi\rangle\langle\phi|=\|\phi\|^2$. $\square$

> **Физический смысл.** Пусть детектор регистрирует каждый скачок — каждый испущенный фотон, каждый потерянный атом.
> Тогда $w_n(t)$ — вероятность ровно $n$ срабатываний к моменту $t$, $S_n(t)\rho/w_n(t)$ — состояние
> при условии $n$ срабатываний, а $S_0(t)\rho/w_0(t)$ — состояние при условии, что **срабатываний не было**: нулевой результат тоже меняет
> состояние — через неэрмитов $H_{\mathrm{eff}}$. Эта интерпретация в терминах считающего процесса [E. B. Davies,
> *Quantum Theory of Open Systems*, 1976; M. D. Srinivas and E. B. Davies, *Optica Acta* 28, 981 (1981)]
> лежит в основе метода квантовых траекторий (метода волновых функций Монте-Карло) [J. Dalibard, Y. Castin,
> K. Mølmer, *Phys. Rev. Lett.* 68, 580 (1992)]; здесь мы используем только её математическую часть.

Для потерь подсчёт срабатываний — то же самое, что измерение заряда:

**Предложение 8.11.** Пусть секторы занумерованы как $N=0,\dots,C-1$, скачков внутри
секторов нет, каждый скачок понижает заряд на единицу ($K^{dc}_{m,l}=0$, если только не $d=c-1$), и
$\rho=P_{N_0}\rho P_{N_0}$. Тогда $S_n(t)\rho=P_{N_0-n}\,\rho(t)\,P_{N_0-n}$ при $0\le n\le N_0$ и
$S_n(t)\rho=0$ при $n>N_0$. В частности, обновление Людерса состояния $\rho(t)$ с помощью $P_{N_0}$ («ни одна частица не потеряна»)
есть состояние без скачков:
$$
\frac{P_{N_0}\rho(t)P_{N_0}}{p_{N_0}(t)}=\frac{e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*}}{w_0(t)},
\qquad p_{N_0}(t)=w_0(t).
$$

*Доказательство.* $S_0(s)$ переводит операторы, сосредоточенные в секторе $N$, в операторы, сосредоточенные в секторе $N$, поскольку
$H_{\mathrm{eff}}\in\mathcal M$; $\Phi_{\mathrm{tot}}$ переводит сектор $N\ge1$ в сектор $N-1$ и аннулирует сектор $0$.
По индукции $S_n(s)\rho$ сосредоточено в секторе $N_0-n$ при всех $s$ и равно $0$ при $n>N_0$. Сжатие
$\rho(t)=\sum_nS_n(t)\rho$ с помощью $P_{N_0-n}$ оставляет в точности $n$-е слагаемое. $\square$

Функция библиотеки `states.lueders_update(rho, P)` возвращает $P\rho P/\operatorname{Tr}(P\rho P)$ и
$\operatorname{Tr}(\rho P)$. Пусть теперь две моды теряют частицы с разными скоростями $\gamma_1>\gamma_2$; тогда
$A=\gamma_1n_1+\gamma_2n_2$ не централен, предложение 8.9(2) неприменимо, и постселекция на
$N=2$ даёт по-настоящему неэрмитову эволюцию. Проверим предложение 8.11 для $n=0$ и, вычислив
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

При $t=2$ вероятность отсутствия потерь равна $0.230$; без перескоков она была бы $e^{-2\gamma_1t}=0.091$, так что
здесь веса секторов действительно зависят от гамильтониана. Постселектированное состояние находится на следовом расстоянии
$0.135$ от состояния замкнутой системы: то, что частица не потеряна, — это информация.

## 8.5 Стационарные состояния и релаксация

Состояние $\sigma$ называется **стационарным**, если $e^{t\mathcal L}\sigma=\sigma$ при всех $t\ge0$, или, что равносильно,
$\mathcal L(\sigma)=0$ (продифференцируйте в $t=0$; обратно, воспользуйтесь рядом).

**Предложение 8.12.** У всякой квантовой динамической полугруппы на $\mathcal M$ есть стационарное состояние.

*Доказательство.* Зафиксируем состояние $\rho$ и положим $\sigma_\tau=\frac1\tau\int_0^\tau T_s\rho\,ds$. Это состояние (предел
выпуклых комбинаций состояний; состояния образуют замкнутое множество). По линейности
$T_t\sigma_\tau-\sigma_\tau=\frac1\tau\big(\int_\tau^{\tau+t}-\int_0^t\big)T_s\rho\,ds$, и следовая норма этой разности
не превосходит $2t/\tau$, поскольку $\|T_s\rho\|_1=1$. Состояния образуют ограниченное замкнутое множество в конечномерном
пространстве, поэтому некоторая последовательность $\sigma_{\tau_m}$, $\tau_m\to\infty$, сходится к состоянию $\sigma$
(теорема Больцано–Вейерштрасса), и $T_t\sigma-\sigma=\lim_m(T_t\sigma_{\tau_m}-\sigma_{\tau_m})=0$. $\square$

**Следствие 8.13.** Если скачков между секторами нет, то стационарные состояния порождают пространство
размерности не меньше $C$: в каждом секторе есть своё стационарное состояние.

*Доказательство.* $\mathcal L$ переводит $P_c\mathcal MP_c\cong M_{k_c}$ в себя (все $H$, $L_j$ блочно-диагональны),
поэтому ограничение $e^{t\mathcal L}$ — квантовая динамическая полугруппа на $M_{k_c}$ (лемма 8.4(a)) со
стационарным состоянием $\sigma_c$ по предложению 8.12. Носители $\sigma_c$ не пересекаются. $\square$

**Предложение 8.14 (спонтанное излучение).** Пусть на $M_2$ с базисом $|g\rangle,|e\rangle$
$H=\omega|e\rangle\langle e|$, один скачок $L=|g\rangle\langle e|$ и скорость $\gamma>0$. Тогда
$$
\rho_{ee}(t)=e^{-\gamma t}\rho_{ee}(0),\qquad\rho_{ge}(t)=e^{(i\omega-\gamma/2)t}\rho_{ge}(0),
\qquad e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}},
$$
где $\mathcal U_t(\rho)=e^{-itH}\rho e^{itH}$, а $\mathcal A_p$ — канал амплитудного затухания с
параметром $p$ (раздел 6.4). Спектр $\mathcal L$ равен $\{0,-\gamma,-\gamma/2\pm i\omega\}$,
$|g\rangle\langle g|$ — единственное стационарное состояние, и
$\tfrac12\|\rho(t)-|g\rangle\langle g|\|_1=\big(\rho_{ee}(0)^2e^{-2\gamma t}+|\rho_{ge}(0)|^2e^{-\gamma t}\big)^{1/2}$.

*Доказательство.* Здесь $L^*L=|e\rangle\langle e|$ и $L\rho L^*=\rho_{ee}|g\rangle\langle g|$, поэтому поэлементно
$\rho_{ee}'=-\gamma\rho_{ee}$, $\rho_{gg}'=\gamma\rho_{ee}$ и $\rho_{ge}'=(i\omega-\gamma/2)\rho_{ge}$, что
даёт решение. В базисе $(\rho_{gg},\rho_{ge},\rho_{eg},\rho_{ee})$ матрица $\mathcal L$ верхнетреугольна
с диагональю $(0,-\gamma/2+i\omega,-\gamma/2-i\omega,-\gamma)$. Из $\mathcal L(\sigma)=0$ следует
$\sigma_{ee}=\sigma_{ge}=0$. Канал амплитудного затухания умножает $\rho_{ee}$ на $1-p$, а $\rho_{ge}$ — на
$\sqrt{1-p}$ (его операторы Крауса — $|g\rangle\langle g|+\sqrt{1-p}\,|e\rangle\langle e|$ и
$\sqrt p\,|g\rangle\langle e|$), а $\mathcal U_t$ умножает $\rho_{ge}$ на $e^{i\omega t}$; оба канала являются TP, что
определяет $\rho_{gg}$. Наконец, $\rho(t)-|g\rangle\langle g|$ — самосопряжённая матрица со следом нуль и диагональю
$(-\rho_{ee}(t),\rho_{ee}(t))$, поэтому её собственные значения равны $\pm(\rho_{ee}(t)^2+|\rho_{ge}(t)|^2)^{1/2}$. $\square$

Когерентность затухает вдвое медленнее населённости — это соотношение $T_2=2T_1$ из раздела 6.4 в отсутствие
чистой дефазировки;
наименьшая ненулевая скорость, **спектральная щель** $\gamma/2$, управляет релаксацией. Проверим спектр
с помощью `lindblad_superoperator`, канал и следовое расстояние — с помощью `lindblad_channel`, а
следствие 8.13 — для случайных скачков внутри секторов $M_2\oplus M_3$:

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

Для этих случайных скачков код находит ровно одно стационарное состояние в каждом секторе. Критерии
единственности — например, существование точного стационарного состояния вместе с $\{H,L_j,L_j^*\}'=\mathbb C1$
[A. Frigerio, *Commun. Math. Phys.* 63, 269 (1978)] — здесь не доказываются; упражнение 4 показывает, что
единственность может нарушаться даже в факторе.

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

## Упражнения

1. (Картина Гейзенберга) Покажите, что $\operatorname{Tr}(\mathcal L(\rho)A)=\operatorname{Tr}(\rho\,\mathcal L^*(A))$ для
   генератора из раздела 8.3, где
   $\mathcal L^*(A)=i[H,A]+\sum_jg_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})+\sum_mg_m(\Phi_m^*(A)-\frac12\{\Phi_m^*(1),A\})$,
   что $\mathcal L^*(1)=0$ и что $\mathcal L^*(P_c)=0$, если скачков между секторами нет.
2. (Исключительная точка) В секторе $N=1$ из раздела 8.4
   $H_{\mathrm{eff}}=\begin{pmatrix}-i\gamma_1/2&J\\J&-i\gamma_2/2\end{pmatrix}$. Найдите его собственные значения и покажите, что
   при $J=|\gamma_1-\gamma_2|/4>0$ он не диагонализуем.
3. Пусть каждый из $N_0$ бозонов в $m$ модах теряется со скоростью $\gamma$ (скачки $a_1,\dots,a_m$, ограниченные на
   $N\to N-1$). Покажите, что $p_N(t)=\binom{N_0}Nq^N(1-q)^{N_0-N}$ с $q=e^{-\gamma t}$ для любого гамильтониана
   из $\mathcal M$.
4. (Чистая дефазировка) Пусть на $M_2$ $H=\omega|e\rangle\langle e|$ и $L=\sigma_z=\operatorname{diag}(1,-1)$ со скоростью
   $\gamma$. Покажите, что населённости постоянны, $\rho_{ge}(t)=e^{(i\omega-2\gamma)t}\rho_{ge}(0)$ и что
   стационарные состояния порождают двумерное пространство диагональных матриц, хотя $C=1$.

<details>
<summary>Решения</summary>

1. В силу цикличности $\operatorname{Tr}(L\rho L^*A)=\operatorname{Tr}(\rho L^*AL)$,
   $\operatorname{Tr}(\{L^*L,\rho\}A)=\operatorname{Tr}(\rho\{L^*L,A\})$ и
   $\operatorname{Tr}(-i[H,\rho]A)=\operatorname{Tr}(\rho\,i[H,A])$; для скачков между секторами
   $\operatorname{Tr}(\Phi(\rho)A)=\operatorname{Tr}(\rho\Phi^*(A))$ — это определение $\Phi^*$
   (раздел 6.3). Тогда $\mathcal L^*(1)=\sum_jg_j(L_j^*L_j-L_j^*L_j)+\sum_mg_m(\Phi_m^*(1)-\Phi_m^*(1))=0$,
   что ещё раз выражает сохранение следа. Без скачков между секторами $P_c$ коммутирует с $H$ и
   всеми $L_j$, поэтому $L_j^*P_cL_j=L_j^*L_jP_c=\frac12\{L_j^*L_j,P_c\}$ и $\mathcal L^*(P_c)=0$; тогда
   $p_c'=\operatorname{Tr}(\rho\mathcal L^*(P_c))=0$ — второе доказательство предложения 8.9(1).

2. $H_{\mathrm{eff}}=-\frac i4(\gamma_1+\gamma_2)1+\begin{pmatrix}-i\delta&J\\J&i\delta\end{pmatrix}$ с
   $\delta=(\gamma_1-\gamma_2)/4$; у второй матрицы след $0$ и определитель $\delta^2-J^2$, поэтому
   собственные значения равны $-\frac i4(\gamma_1+\gamma_2)\pm\sqrt{J^2-\delta^2}$: две скорости затухания при $J<|\delta|$,
   одна скорость и две частоты при $J>|\delta|$. При $J=|\delta|>0$ вторая матрица $Y$ удовлетворяет
   $Y^2=(J^2-\delta^2)1=0$, но $Y\ne0$, поэтому $H_{\mathrm{eff}}$ имеет единственное собственное значение
   $-\frac i4(\gamma_1+\gamma_2)$ и не кратен $1$: он не диагонализуем.

```python
lam = np.linalg.eigvals(H_eff[1])
delta = (g1 - g2) / 4
root = np.sqrt(J ** 2 - delta ** 2 + 0j)
assert np.allclose(sorted(lam, key=np.real), sorted(-0.25j * (g1 + g2) + np.array([root, -root]), key=np.real))
H_ep = np.array([[-0.5j * g1, delta], [delta, -0.5j * g2]])         # J = delta: the exceptional point
Y = H_ep + 0.25j * (g1 + g2) * np.eye(2)                           # H_ep minus its only eigenvalue
assert np.abs(Y @ Y).max() < 1e-15 and np.linalg.matrix_rank(Y) == 1   # one Jordan block of size 2
```

3. На секторе $N$ имеем $\sum_ja_j^*a_j=N1$, и единственная пара секторов — $(N-1,N)$, поэтому предложение
   8.9(2) даёт $p_N'=-\gamma Np_N+\gamma(N+1)p_{N+1}$. С учётом $q'=-\gamma q$ производная биномиального
   слагаемого равна $-\gamma Np_N+\gamma(N_0-N)\binom{N_0}Nq^{N+1}(1-q)^{N_0-N-1}$, а
   $(N+1)\binom{N_0}{N+1}=(N_0-N)\binom{N_0}N$ показывает, что второе слагаемое равно $\gamma(N+1)p_{N+1}$. При $t=0$
   это единичная масса в точке $N_0$; единственность — это лемма 8.1.

4. $L\rho L^*-\frac12\{L^*L,\rho\}=\sigma_z\rho\sigma_z-\rho$ умножает внедиагональные элементы на $-2$ и
   обнуляет диагональ, а $-i[H,\rho]_{ge}=i\omega\rho_{ge}$. Поэтому $\mathcal L(\sigma)=0$ тогда и только тогда, когда $\sigma_{ge}=0$ (при
   $\gamma>0$), и спектр равен $\{0,0,-2\gamma\pm i\omega\}$. Поскольку $H$ и $\sigma_z$ коммутируют, а
   $X\mapsto\sigma_zX\sigma_z$ — инволюция, $e^{t\mathcal L}(\rho)=\frac{1+e^{-2\gamma t}}2U_t\rho U_t^*
   +\frac{1-e^{-2\gamma t}}2U_t\sigma_z\rho\sigma_zU_t^*$ с $U_t=e^{-itH}$: это унитальный канал, неподвижные
   точки которого — матрицы, коммутирующие с его операторами Крауса (теорема 6.10), т. е. диагональные.

```python
deph = dynamics.lindblad_superoperator(q2, H_at, [q2.from_blocks([torch.diag(torch.tensor([1.0, -1.0], dtype=f64))])],
                                       rates=[gamma])[0][0]
ev = torch.linalg.eigvals(deph)
assert int((ev.abs() < 1e-12).sum()) == 2
assert np.allclose(sorted(ev.tolist(), key=lambda z: (z.real, z.imag))[:2],
                   [-2 * gamma - 1j * omega, -2 * gamma + 1j * omega], atol=1e-12)
```

</details>

**Далее:** в [главе 9](09_spin_chains.md) рассматриваются спиновые цепочки: секторы намагниченности и импульса,
эволюция во времени методами Ланцоша и Крылова для больших секторов и беспорядок.
