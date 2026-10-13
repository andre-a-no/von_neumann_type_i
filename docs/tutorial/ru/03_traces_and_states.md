[English](../03_traces_and_states.md) · [Español](../es/03_traces_and_states.md) · [Français](../fr/03_traces_and_states.md) · [Deutsch](../de/03_traces_and_states.md) · [中文](../zh/03_traces_and_states.md) · [日本語](../ja/03_traces_and_states.md) · **Русский**

# 3. Следы, состояния и энтропия

На $M_n(\mathbb C)$ след по существу один. На $\bigoplus_cM_{k_c}$ имеется
$C$-параметрическое семейство следов, и выбор существен: от него зависят плотности,
энтропии и сопряжённые отображения. В этой главе классифицируются следы, состояния
представляются плотностями, доказываются основные энтропийные неравенства, а в конце
разбирается правдоподобное следовое неравенство, которое оказывается неверным.

**Вы узнаете** три следа библиотеки и способы перехода между ними,
классическую часть состояния, разложение энтропии фон Неймана по
секторам, неравенство Клейна и его следствия, а также контраст Майкельсона.

Всюду далее $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$, где все $k_c\ge1$,
$D=\sum_ck_c$, $P_c$ — минимальные центральные проекторы и $A=\bigoplus_cA_c$.

## 3.1 Следы

**Следом** на $\mathcal M$ называется линейный функционал $\tau$, для которого $\tau(AB)=\tau(BA)$ при всех
$A,B$. След **положителен**, если $\tau(A^*A)\ge0$ для всех $A$, и **точен**, если, кроме того,
из $\tau(A^*A)=0$ следует $A=0$.

**Теорема 3.1.** Всякий след на $\mathcal M$ имеет вид
$\operatorname{tr}_w(A)=\sum_cw_c\operatorname{Tr}A_c$ с $w\in\mathbb C^C$. Он положителен
тогда и только тогда, когда все $w_c\ge0$, и точен тогда и только тогда, когда все $w_c>0$.

*Доказательство.* Пусть $E^c_{ij}$ — матричные единицы $c$-го блока. Для следа $\tau$
$\tau(E^c_{ij})=\tau(E^c_{i1}E^c_{1j})=\tau(E^c_{1j}E^c_{i1})=\delta_{ij}\tau(E^c_{11})$.
Положим $w_c=\tau(E^c_{11})$; по линейности $\tau(A)=\sum_c\sum_i(A_c)_{ii}w_c$. Обратно,
всякий $\operatorname{tr}_w$ является следом. Положительность: $E^c_{11}=(E^c_{11})^*E^c_{11}$, поэтому
условие $w_c\ge0$ необходимо, а равенство $\operatorname{tr}_w(A^*A)=\sum_cw_c\|A_c\|_2^2$ доказывает достаточность
и критерий точности. $\square$

В библиотеке реализованы три точных следа:

| название | веса $w_c$ | значение на $1$ | смысл |
|---|---|---|---|
| `Tr_blunt` | $1$ | $D$ | след $\bigoplus_c\mathbb C^{k_c}$ (кратности $m_c=1$) |
| `Tr_norm` | $1/k_c$ | $C$ | каждый блок нормирован |
| `tau_vN` | $1/(Ck_c)$ | $1$ | следовое состояние с равными весами секторов $1/C$ |

Согласно главе 2, след пространства, на котором $\mathcal M$ действует с кратностями $m_c$,
есть $\operatorname{tr}_w$ с $w_c=m_c$. Таким образом, выбор точного следа — это выбор того, «какой
вес имеет каждый сектор».

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

## 3.2 Состояния и плотности

**Состоянием** называется линейный функционал $\omega$ на $\mathcal M$, который положителен
($\omega(A^*A)\ge0$) и нормирован ($\omega(1)=1$).

**Теорема 3.2.** Пусть $\operatorname{tr}_w$ — точный след. Для всякого линейного функционала
$\omega$ существует единственный $\rho_w\in\mathcal M$, для которого $\omega(A)=\operatorname{tr}_w(\rho_wA)$ при всех
$A$. Функционал $\omega$ положителен тогда и только тогда, когда $\rho_w\ge0$, и является
состоянием тогда и только тогда, когда, кроме того, $\operatorname{tr}_w(\rho_w)=1$. Для двух точных следов
плотности связаны поблочно соотношением
$$
w_c\,(\rho_w)_c=w'_c\,(\rho_{w'})_c .
$$

*Доказательство.* $\langle X,Y\rangle_w=\operatorname{tr}_w(X^*Y)$ — скалярное произведение на $\mathcal M$
(теорема 3.1). По теореме Рисса о представлении в конечномерном пространстве
$\mathcal M$ существует единственный $R$ с $\omega(A)=\langle R,A\rangle_w=\operatorname{tr}_w(R^*A)$; положим
$\rho_w=R^*$. Если $\rho_w\ge0$, то
$\omega(A^*A)=\sum_cw_c\operatorname{Tr}(A_c\,(\rho_w)_c\,A_c^*)\ge0$, поскольку
$A_c(\rho_w)_cA_c^*\ge0$ (следствие 1.7 с $T=(\rho_w)_c^{1/2}A_c^*$). Обратно, если
$\omega$ положителен, возьмём $x\in\mathbb C^{k_c}$ и в качестве $A$ — оператор ранга один
$|x\rangle\langle x|$ в блоке $c$; тогда $A=A^*A/\|x\|^2$ и
$0\le\omega(A)=w_c\langle x,(\rho_w)_cx\rangle$, так что $(\rho_w)_c\ge0$. Нормировка означает
$\omega(1)=\operatorname{tr}_w(\rho_w)$. Наконец, и $\sum_cw_c\operatorname{Tr}((\rho_w)_cA_c)$, и
$\sum_cw'_c\operatorname{Tr}((\rho_{w'})_cA_c)$ равны $\omega(A)$ для всех $A$, и в силу единственности
$w_c(\rho_w)_c=w'_c(\rho_{w'})_c$. $\square$

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

> **Физический смысл.** $\omega(A)$ — среднее значение наблюдаемой $A$.
> Плотность $\rho$ относительно следа `Tr_blunt` — это обычная матрица плотности квантовой
> механики на $\bigoplus_c\mathbb C^{k_c}$. Плотности относительно других следов возникают,
> когда та же алгебра действует с кратностями (глава 2), например на
> симметричной и несимметричной частях спиновой системы.

## 3.3 Классическая часть состояния

Ограничение $\omega$ на центр — это распределение вероятностей на секторах,
$p_c=\omega(P_c)=\operatorname{Tr}\rho_c$. Если $p_c>0$, то **условное состояние** в секторе $c$
имеет плотность $\sigma_c=\rho_c/p_c$ на $\mathbb C^{k_c}$, и
$$
\rho=\bigoplus_cp_c\,\sigma_c
$$
(если $p_c=0$, то $\rho_c=0$, поскольку $\rho_c\ge0$ и $\operatorname{Tr}\rho_c=0$, и это слагаемое
опускается).
Состояние на $\mathcal M$ — это классическая смесь квантовых состояний в секторах. Измерение
«в каком секторе» его не возмущает: канал измерения
$\rho\mapsto\sum_cP_c\rho P_c$ тождествен на $\mathcal M$, поскольку каждый элемент
$\mathcal M$ блочно-диагонален. В этом и состоит операциональное содержание правила суперотбора.

```python
p = rho.sector_probabilities()                                   # shape (batch, C)
assert torch.allclose(p.sum(-1), torch.ones(4, dtype=p.dtype))
P = [alg.central([1.0 if c == d else 0.0 for d in range(alg.C)]) for c in range(alg.C)]
measure_sector = channels.measurement_channel(P)
assert (measure_sector(rho).matrix - rho.matrix).abs().max() < 1e-12   # no disturbance
```

## 3.4 Энтропия

**Энтропия фон Неймана** состояния с плотностью $\rho$ (относительно `Tr_blunt`) равна
$S(\rho)=-\operatorname{Tr}\rho\log\rho=-\sum_i\lambda_i\log\lambda_i$, где сумма берётся по собственным значениям
$\rho$ и $0\log0=0$.

**Теорема 3.3 (энтропия смеси секторов).** При $p_c$ и $\sigma_c$, как выше
(секторы с $p_c=0$ опущены),
$$
S(\rho)=H(p)+\sum_cp_c\,S(\sigma_c),\qquad H(p)=-\sum_cp_c\log p_c .
$$

*Доказательство.* Ненулевые собственные значения $\rho$ — это числа $p_c\mu$, где $\mu$ пробегает ненулевые
собственные значения $\sigma_c$ (для секторов с $p_c>0$); нулевые собственные значения дают вклад $0$
в обе части. Следовательно,
$-\sum_c\sum_{\mu>0}p_c\mu\log(p_c\mu)=-\sum_cp_c\log p_c\sum_{\mu>0}\mu-\sum_cp_c\sum_{\mu>0}\mu\log\mu$,
и $\sum_{\mu>0}\mu=\operatorname{Tr}\sigma_c=1$. $\square$

> **Физический смысл.** Формула делит неопределённость системы на две части: $H(p)$ —
> классическую неопределённость заряда (в каком секторе находится система) и
> $\sum_cp_cS(\sigma_c)$ — среднюю квантовую неопределённость внутри секторов. При правиле
> суперотбора первая часть есть энтропия классической случайной величины — заряда.

Энтропия относительно другого следа равна $S_w(\rho)=-\operatorname{tr}_w(\rho_w\log\rho_w)$. По
теореме 3.2 $(\rho_w)_c=\rho_c/w_c$, поэтому
$$
S_w(\rho)=-\sum_cw_c\operatorname{Tr}\Big(\frac{\rho_c}{w_c}\log\frac{\rho_c}{w_c}\Big)
=S(\rho)+\sum_cp_c\log w_c .
$$
Энтропия — характеристика состояния *и* следа: одно и то же состояние имеет разные
энтропии относительно `Tr_blunt` и `tau_vN`.

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

## 3.5 Относительная энтропия и неравенство Клейна

Для состояний с плотностями $\rho,\sigma$, $\sigma>0$, **относительная энтропия** равна
$S(\rho\|\sigma)=\operatorname{Tr}\rho\log\rho-\operatorname{Tr}\rho\log\sigma$, где
$\operatorname{Tr}\rho\log\rho:=\sum_i\lambda_i\log\lambda_i$ по собственным значениям $\rho$ с
$0\log0=0$ (и аналогично в каждом блоке). Она не зависит от следа: по теореме 3.2
$(\rho_w)_c=\rho_c/w_c$ и $(\sigma_w)_c=\sigma_c/w_c$, так что собственные значения $(\rho_w)_c$ — это
собственные значения $\rho_c$, делённые на $w_c$, и как $\operatorname{tr}_w(\rho_w\log\rho_w)$, так и
$\operatorname{tr}_w(\rho_w\log\sigma_w)$ отличаются от $\operatorname{Tr}\rho\log\rho$ и
$\operatorname{Tr}\rho\log\sigma$ на одну и ту же величину $-\sum_cp_c\log w_c$ (логарифм $\rho$ нужен
лишь на её носителе, по соглашению $0\log0=0$). Разность не меняется.

**Лемма 3.4 (неравенство для касательной).** Пусть $f:I\to\mathbb R$ выпукла на промежутке $I$ и
дифференцируема в точке $b\in I$. Тогда $f(a)\ge f(b)+(a-b)f'(b)$ для всех $a\in I$, причём
при $a\neq b$ неравенство строгое, если $f$ строго выпукла.

*Доказательство.* Для $t\in(0,1]$ выпуклость даёт $f(b+t(a-b))\le(1-t)f(b)+tf(a)$, т. е.
$\big(f(b+t(a-b))-f(b)\big)/t\le f(a)-f(b)$; переходя к пределу при $t\to0$, получаем неравенство. Если $f$
строго выпукла и $a\neq b$, положим $m=(a+b)/2$; тогда
$f(a)-f(b)>2\big(f(m)-f(b)\big)\ge2(m-b)f'(b)=(a-b)f'(b)$: первый шаг — по строгой
выпуклости в точке $m$, второй — по только что доказанному неравенству, применённому к $m$. $\square$

**Теорема 3.5 (неравенство Клейна).** Пусть $f:I\to\mathbb R$ выпукла на промежутке $I$,
$X,Y$ самосопряжены и их спектры лежат в $I$, причём $f$ дифференцируема в каждом собственном значении $Y$.
Тогда $\operatorname{Tr}[f(X)-f(Y)-(X-Y)f'(Y)]\ge0$. Если $f$ строго выпукла, равенство достигается
только при $X=Y$.

*Доказательство.* Пусть $X=\sum_ia_i|u_i\rangle\langle u_i|$ и $Y=\sum_jb_j|v_j\rangle\langle v_j|$ —
разложения по ортонормированным собственным базисам, и $c_{ij}=|\langle u_i,v_j\rangle|^2$, так что $\sum_ic_{ij}=\sum_jc_{ij}=1$.
Тогда $\operatorname{Tr}f(X)=\sum_{ij}c_{ij}f(a_i)$, $\operatorname{Tr}f(Y)=\sum_{ij}c_{ij}f(b_j)$ и
$\operatorname{Tr}Xf'(Y)=\sum_{ij}c_{ij}a_if'(b_j)$, $\operatorname{Tr}Yf'(Y)=\sum_{ij}c_{ij}b_jf'(b_j)$.
Следовательно, левая часть равна $\sum_{ij}c_{ij}[f(a_i)-f(b_j)-(a_i-b_j)f'(b_j)]$ — сумме
неотрицательных слагаемых в силу леммы 3.4. Если $f$ строго выпукла и сумма равна нулю, то
$a_i=b_j$ всякий раз, когда $c_{ij}>0$, поэтому $Xv_j=\sum_ia_i\langle u_i,v_j\rangle u_i=b_jv_j$ для каждого $j$,
и $X=Y$. $\square$

**Следствие 3.6.** $S(\rho\|\sigma)\ge0$, причём равенство имеет место тогда и только тогда, когда $\rho=\sigma$. В частности,
$S(\rho)\le\log D$, причём равенство имеет место тогда и только тогда, когда $\rho=1/D$.

*Доказательство.* Применим теорему 3.5 к функции $f(t)=t\log t$, $f(0)=0$, которая строго выпукла на
$I=[0,\infty)$ ($f''(t)=1/t>0$ на $(0,\infty)$, и $f$ непрерывна в $0$) и
дифференцируема на $(0,\infty)\supseteq\operatorname{spec}\sigma$ с $f'(t)=\log t+1$, при
$X=\rho$, $Y=\sigma$: левая часть равна
$\operatorname{Tr}\rho\log\rho-\operatorname{Tr}\sigma\log\sigma-\operatorname{Tr}(\rho-\sigma)(\log\sigma+1)
=S(\rho\|\sigma)$, поскольку $\operatorname{Tr}\rho=\operatorname{Tr}\sigma=1$. Для второго утверждения возьмём
$\sigma=1/D$: $0\le S(\rho\|1/D)=-S(\rho)+\log D$. $\square$

**Теорема 3.7 (цепное правило; монотонность при условном ожидании на центр).** Пусть $\rho,\sigma$ — состояния,
$\sigma>0$, $p_c=\operatorname{Tr}\rho_c$, $q_c=\operatorname{Tr}\sigma_c>0$, $\hat\sigma_c=\sigma_c/q_c$ и
$\hat\rho_c=\rho_c/p_c$ при $p_c>0$. Тогда, если суммы берутся по тем $c$, для которых $p_c>0$ (остальные
секторы дают вклад $0$ в силу $0\log0=0$),
$$
S(\rho\|\sigma)=\sum_cp_c\log\frac{p_c}{q_c}+\sum_cp_c\,S(\hat\rho_c\|\hat\sigma_c)
\ \ge\ \sum_cp_c\log\frac{p_c}{q_c}=S(E\rho\|E\sigma),
$$
где $E$ — условное ожидание на центр (предложение 2.8).

*Доказательство.* При $p_c>0$ имеем $\operatorname{Tr}\rho_c\log\rho_c=p_c\log p_c+p_c\operatorname{Tr}\hat\rho_c\log\hat\rho_c$
(собственные значения $p_c\mu$, как в теореме 3.3) и $\log\sigma_c=(\log q_c)1+\log\hat\sigma_c$, поэтому
$\operatorname{Tr}\rho_c\log\sigma_c=p_c\log q_c+p_c\operatorname{Tr}\hat\rho_c\log\hat\sigma_c$ и
$\operatorname{Tr}\rho_c\log\rho_c-\operatorname{Tr}\rho_c\log\sigma_c=p_c\log(p_c/q_c)+p_cS(\hat\rho_c\|\hat\sigma_c)$.
Если $p_c=0$, то $\rho_c=0$, и вклад сектора равен $0$.
Просуммируем по $c$; вторая сумма $\ge0$ по следствию 3.6, применённому в каждом блоке. Наконец,
$E\rho=\bigoplus_c(p_c/k_c)1_{k_c}$, $E\sigma=\bigoplus_c(q_c/k_c)1_{k_c}$ и
$\operatorname{Tr}E\rho(\log E\rho-\log E\sigma)=\sum_ck_c\frac{p_c}{k_c}\log\frac{p_c}{q_c}$. $\square$

Итак, если забыть квантовую информацию внутри секторов, различимость двух состояний может
лишь уменьшиться; остаётся классическое расхождение Кульбака–Лейблера
распределений по секторам.

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

## 3.6 Контраст Майкельсона

Для $A\ge0$, $A\neq0$, **контраст Майкельсона** равен
$\Delta(A)=(\lambda_{\max}-\lambda_{\min})/(\lambda_{\max}+\lambda_{\min})\in[0,1]$, где
крайние собственные значения берутся по всем секторам (и $\Delta(0)=0$). Он измеряет, насколько $A$
далёк от кратного единицы: $\Delta(A)=0$ тогда и только тогда, когда $A=\lambda1$.

**Предложение 3.8.** Если $\Phi:\mathcal M\to\mathcal M$ линейно, положительно и унитально
($\Phi(1)=1$), то $\Delta(\Phi(A))\le\Delta(A)$ для всех $A\ge0$.

*Доказательство.* $\lambda_{\min}1\le A\le\lambda_{\max}1$, и из положительности и унитальности следует
$\lambda_{\min}1\le\Phi(A)\le\lambda_{\max}1$. Значит,
$\lambda_{\min}(\Phi(A))\ge\lambda_{\min}$ и $\lambda_{\max}(\Phi(A))\le\lambda_{\max}$. Если
$\Phi(A)=0$, утверждение тривиально. В противном случае при $r=\lambda_{\min}/\lambda_{\max}$
величина $\Delta=(1-r)/(1+r)$ убывает по $r$, а $r(\Phi(A))\ge r(A)$. $\square$

```python
A_pos = alg.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(4, d), batch_size=4, force_positive=True)
Phi = channels.random_mixed_unitary_channel(alg, 3, batch_size=4)       # unital
assert (Phi(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
assert (E(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
```

## 3.7 Правдоподобное неравенство, которое неверно

Верно ли для положительных $X,Y$ из $M_k$ и унитарной $U$, что $|\operatorname{Tr}(XUY)|\le\operatorname{Tr}(XY)$? Это
так при $U=1$, и утверждение похоже на неравенство типа Коши–Буняковского. Однако оно неверно всякий раз,
когда $X$ и $Y$ не коммутируют.

**Предложение 3.9.** Пусть $X,Y\ge0$ в $M_k$. Тогда $\max_U|\operatorname{Tr}(XUY)|=\|YX\|_1\ge\operatorname{Tr}(XY)$, где максимум берётся по
унитарным $U\in M_k$, причём равенство имеет место тогда и только тогда, когда $XY=YX$.

*Доказательство.* $\operatorname{Tr}(XUY)=\operatorname{Tr}(YX\,U)$, поэтому максимум равен $\|YX\|_1$ по теореме 1.9,
и $\|YX\|_1\ge|\operatorname{Tr}(YX)|=\operatorname{Tr}(XY)$, поскольку
$\operatorname{Tr}(XY)=\operatorname{Tr}(X^{1/2}YX^{1/2})\ge0$, так как
$X^{1/2}YX^{1/2}=(Y^{1/2}X^{1/2})^*(Y^{1/2}X^{1/2})\ge0$ (следствие 1.7). Если $XY=YX$, то
$X^{1/2}$ коммутирует с $Y$ (теорема 1.5(d)), поэтому $YX=X^{1/2}YX^{1/2}\ge0$ и
$\|YX\|_1=\operatorname{Tr}(YX)$ (для $P\ge0$ имеем $\|P\|_1=\operatorname{Tr}P$, раздел 1.4). Обратно,
пусть $\|YX\|_1=\operatorname{Tr}(YX)$. Если $|YX|=0$, то $YX=0=(YX)^*=XY$. В противном случае запишем
$YX=W|YX|$. По теореме 1.9 равенство $\operatorname{Tr}(YX\cdot1)=\|YX\|_1$ означает равенство в
шаге с неравенством Коши–Буняковского при $U=1$:
$|YX|^{1/2}$ и $W|YX|^{1/2}$ пропорциональны с неотрицательным множителем, который при
сравнении норм оказывается равным $1$, так что $W|YX|^{1/2}=|YX|^{1/2}$ и $YX=W|YX|=|YX|\ge0$. Положительная
матрица самосопряжена, поэтому $YX=(YX)^*=XY$. $\square$

Библиотека находит точку максимума с помощью полярного разложения. В сопроводительной статье
численно показано, что выборка $U$ методом Монте-Карло до максимума не дотягивает, и тем
сильнее, чем больше $k$, — это одна из её центральных тем; код ниже демонстрирует это для одной пары.
Причина в том, что $U=1$ уже близко к оптимуму, а для случайного по Хаару $U$ в большой
размерности $|\operatorname{Tr}(XUY)|$, как правило, гораздо меньше $\operatorname{Tr}(XY)$ (случайная
унитарная матрица перемешивает фазы), поэтому все выборки остаются далеко ниже точного супремума.

```python
Mk = TypeIAlgebra([6], [6], precision='double')
X = Mk.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
Y = Mk.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
YX = Y @ X
W = YX @ YX.abs().inverse()                                          # polar part of YX
gap = (X @ W.adjoint() @ Y).trace.abs() - (X @ Y).trace.real
assert abs(gap.item() - (YX.trace_norm() - (X @ Y).trace.real).item()) < 1e-10
assert gap.item() > 0                                                # the inequality fails
samples = Mk.random_unitary_operator(20_000)
best = (X @ samples @ Y).trace.abs().max() - (X @ Y).trace.real
print(f"exact sup {gap.item():.4f}, best of 20000 Haar samples {best.item():.4f}")
assert best.item() < gap.item()
```

## 3.8 С чем сравнивать

QuTiP (`entropy_vn`, `entropy_relative`) и пакеты квантовой теории информации (например,
`qutip`, `toqito`) вычисляют энтропии и относительные энтропии матриц плотности на
одном гильбертовом пространстве; на $\bigoplus_c\mathbb C^{k_c}$ они согласуются с принятыми здесь
соглашениями для `Tr_blunt`. Чего они не предоставляют, так это семейства следов из теоремы 3.1 и
переходов из теоремы 3.2, которые становятся необходимы, как только алгебра действует с
кратностями или секторы нужно взвешивать по-разному.

## Упражнения

1. Покажите, что плотность следового состояния $\tau_{\mathrm{vN}}$ (равные веса секторов) относительно
   $\tau_{\mathrm{vN}}$ равна $1$, и вычислите его плотность относительно `Tr_blunt` и его энтропию $S$.
2. Пусть $\operatorname{tr}_w$ точен. Покажите, что $S_w(\rho)\le\log\operatorname{tr}_w(1)$, причём равенство имеет место тогда и только тогда, когда
   $\rho_w=1/\operatorname{tr}_w(1)$.
3. Покажите, что канал измерения $\rho\mapsto\sum_cP_c\rho P_c$ на $M_D$ (а не на
   $\mathcal M$!) действительно возмущает состояния, и вычислите $S(\rho)$ до и после измерения для
   суперпозиции двух секторов.

<details>
<summary>Решения</summary>

1. $\tau_{\mathrm{vN}}(A)=\tau_{\mathrm{vN}}(1\cdot A)$, поэтому $\rho_{\tau}=1$ в силу единственности. По
   теореме 3.2 плотность относительно `Tr_blunt` равна $\bigoplus_c1_{k_c}/(Ck_c)$, причём $p_c=1/C$ и
   $\sigma_c=1/k_c$, так что $S=\log C+\frac1C\sum_c\log k_c$ по теореме 3.3.

```python
tau = states.tracial_state(alg)
assert (tau.density('tau_vN') - one.matrix).abs().max() < 1e-12
assert abs(tau.von_neumann_entropy().item() - (math.log(3) + sum(math.log(k) for k in (2, 3, 5)) / 3)) < 1e-12
```

2. Суммы берутся по тем $c$, для которых $p_c>0$. Тогда
   $S_w(\rho)=S(\rho)+\sum_cp_c\log w_c=-\sum_cp_c\log p_c+\sum_cp_cS(\sigma_c)+\sum_cp_c\log w_c
   \le\sum_cp_c\log\frac{k_cw_c}{p_c}\le\log\sum_{p_c>0}k_cw_c\le\log\operatorname{tr}_w(1)$; здесь использованы
   $S(\sigma_c)\le\log k_c$ (следствие 3.6) и конечное неравенство Йенсена
   $\sum_cp_c\log x_c\le\log\sum_cp_cx_c$ для вогнутого $\log$ (оно следует из
   определения вогнутости индукцией по числу слагаемых) с $x_c=k_cw_c/p_c$; это неравенство
   строгое, если не все $x_c$ равны, так как $\log$ строго вогнут. Равенство во всей цепочке
   влечёт, что все $p_c>0$, $\sigma_c=1/k_c$ и $k_cw_c/p_c$ постоянно, т. е.
   $p_c=k_cw_c/\operatorname{tr}_w(1)$, т. е. $(\rho_w)_c=\rho_c/w_c=1_{k_c}/\operatorname{tr}_w(1)$.
   Обратно, при $\rho_w=1/\operatorname{tr}_w(1)$
   $S_w(\rho)=-\operatorname{tr}_w(1)\cdot\frac1{\operatorname{tr}_w(1)}\log\frac1{\operatorname{tr}_w(1)}=\log\operatorname{tr}_w(1)$.

3. В $M_D$ чистая суперпозиция $\psi=(\psi_1\oplus\psi_2)/\sqrt2$ единичных векторов из двух
   секторов имеет $S=0$; после измерения состояние равно
   $\frac12|\psi_1\rangle\langle\psi_1|\oplus\frac12|\psi_2\rangle\langle\psi_2|$, и $S=\log2$. Однако при ограничении на
   $\mathcal M$ векторное состояние $\psi$ совпадает со смесью:
   $\langle\psi,A\psi\rangle=\tfrac12\langle\psi_1,A_1\psi_1\rangle+\tfrac12\langle\psi_2,A_2\psi_2\rangle$
   для всех $A\in\mathcal M$. Никакая наблюдаемая из $\mathcal M$ не обнаруживает когерентности — в этом и
   состоит содержание правила суперотбора.

</details>

**Далее:** в [главе 4](04_symmetries.md) рассматриваются симметрии и представления: как
действие группы на физической системе порождает секторы и как библиотека строит
соответствующие алгебры для спиновых систем и регистров кубитов.
