[English](../11_optimisation.md) · [Español](../es/11_optimisation.md) · [Français](../fr/11_optimisation.md) · [Deutsch](../de/11_optimisation.md) · [中文](../zh/11_optimisation.md) · [日本語](../ja/11_optimisation.md) · **Русский**

# 11. Оптимизация и следовые неравенства

Следовое неравенство — это утверждение об экстремальных значениях функционала от операторов. В разделе 3.7 было
показано правдоподобное, но неверное неравенство (предложение 3.9), а также выборочный запуск, который численно не нашёл
нарушения. В этой главе
объясняется, почему выборка в большой размерности не срабатывает, и она заменяется оптимизацией по точно
параметризованным множествам ограничений; там, где экстремумы известны в замкнутой форме, мы их доказываем, так что
результат оптимизатора сверяется с теоремой.

**Вы узнаете**, как библиотека параметризует унитарные операторы и спектры с заданным контрастом и почему при этом
достигается каждый допустимый оператор; неравенство фон Неймана для следов и точные экстремумы, которые оно даёт; следовое
неравенство вместе с условием равенства; почему выборки по Хаару в среднеквадратичном далеки от экстремумов; и что доказывает
мультистартовый поиск.

Всюду далее $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$, где $D=\sum_ck_c\ge2$, $\operatorname{Tr}$ — след
`Tr_blunt` (раздел 3.1), а $\|A\|_1=\sum_c\|A_c\|_1$ — следовая норма $A$ как блочно-диагональной матрицы из
$M_D(\mathbb C)$. Сингулярные числа $A\in M_k(\mathbb C)$ обозначаются $s_1(A)\ge\dots\ge s_k(A)$; $x^\downarrow$ и
$x^\uparrow$ — упорядочения вещественного вектора $x$ по убыванию и по возрастанию. $\Delta$ — контраст
Майкельсона (раздел 3.6), и для $\delta\in[0,1)$
$$
\ell(\delta)=\frac{1-\delta}{1+\delta},\qquad 1-\ell(\delta)=\frac{2\delta}{1+\delta};
$$
у положительного оператора с $\lambda_{\max}=1$ контраст равен $\ell(\lambda_{\min})$, поэтому его контраст равен $\delta$ тогда и только тогда, когда
$\lambda_{\min}=\ell(\delta)$, поскольку $\ell$ — строго убывающая инволюция отрезка $[0,1]$. Примерами служат три функционала
сопроводительной статьи (раздел *From sampling to optimisation*):

| | функционал | переменные |
|---|---|---|
| Эксп. 1 | $z=\lvert\operatorname{Tr}(XUY)\rvert-\operatorname{Tr}(XY)$ | $X,Y\ge0$, $U\in\mathcal M$ унитарен |
| Эксп. 2 | $z=\operatorname{Tr}\lvert XY\rvert-\operatorname{Tr}(\lvert X\rvert\,\lvert Y\rvert)$ | $X=UX_0$, $Y=VY_0$, $X_0,Y_0\ge0$, $U,V$ унитарны |
| Эксп. 3 | $z=\operatorname{Tr}(Y\lvert X\rvert Y)-\operatorname{Tr}\lvert YXY\rvert$ | $X=X^*$, $Y\ge0$ |

где $|X|$, $|Y|$ нормированы условием $\lambda_{\max}=1$, а их контрасты служат параметрами. **Огибающие** — это
$\overline z(\delta)=\sup\{z:\Delta(|X|)=\delta\}$ и $\underline z(\delta)=\inf\{z:\Delta(|X|)=\delta\}$, где все остальные
переменные свободны.

## 11.1 Множества ограничений и их параметризации

**Определение.** Для $\delta\in[0,1)$ и $\ell=\ell(\delta)$ **оболочкой контраста** называется
$$
\mathcal P_\delta=\{X\in\mathcal M:X\ge0,\ \lambda_{\max}(X)=1,\ \Delta(X)=\delta\}
=\{X=X^*:\operatorname{spec}X\subseteq[\ell,1],\ 1,\ell\in\operatorname{spec}X\}.
$$
**Слоты** алгебры $\mathcal M$ — это $D$ пар $(c,i)$, $i\le k_c$; для $\lambda\in\mathbb R^D$ положим
$\operatorname{diag}(\lambda)=\bigoplus_c\operatorname{diag}(\lambda_{c,1},\dots,\lambda_{c,k_c})$.

`optimize.PositiveParam` с заданным контрастом $\delta$ выбирает для каждого старта два различных слота
$s_\pm=(c_\pm,i_\pm)$ ($s_+$ равномерно среди $D$ слотов, затем $s_-$ равномерно среди $D-1$ остальных, так что каждая
упорядоченная пара различных слотов имеет вероятность $1/(D(D-1))$) и случайный по Хаару унитарный $U_0\in\mathcal M$
(`random_unitary_operator`, лемма 11.10) и отображает $\theta\in\mathbb R^D$ и
$G=\bigoplus_cG_c$ в
$$
X(\theta,G)=U\operatorname{diag}(\lambda)U^*,\qquad U=U_0\,e^{G-G^*},\qquad
\lambda_{s_+}=1,\quad\lambda_{s_-}=\ell,\quad\lambda_s=\ell+(1-\ell)\tfrac{1+\sin\theta_s}2\ \text{otherwise}.
$$
`UnitaryParam` — это одно лишь отображение $G\mapsto U_0e^{G-G^*}$. В отличие от сигмоиды, $(1+\sin\theta)/2$ достигает $0$ и $1$
при конечных $\theta$, что удобно, поскольку экстремальные спектры этой главы (следствие 11.7, предложение 11.14) лежат в
$\ell$ и $1$. Свободный контраст в $[a,b]$ — это ещё один параметр, $\delta=a+(b-a)(1+\sin\eta)/2$; по умолчанию в библиотеке
$[a,b]=[0,0.999]$, и код ниже это значение сохраняет. Значение, достигнутое при этом ограничении, достигается и без него,
поэтому оно по-прежнему является нижней оценкой для супремума без ограничения (верхней оценкой для инфимума без ограничения).

**Лемма 11.1.** Для всякого унитарного $U_0\in\mathcal M$ отображение $G\mapsto U_0e^{G-G^*}$ сюръективно на унитарную группу
$\mathcal M$.

*Доказательство.* Пусть $U$ унитарен. По теореме 1.3 (поблочно) и следствию 1.4(2)
$U_0^*U=\bigoplus_cW_c\operatorname{diag}(e^{i\varphi_{c,j}})_jW_c^*$, где $W_c$ унитарны, $\varphi_{c,j}\in\mathbb R$. Тогда
$S=\bigoplus_cW_c\operatorname{diag}(i\varphi_{c,j})_jW_c^*$ удовлетворяет $S^*=-S$ и $e^S=U_0^*U$ (ряд
$\sum_nS^n/n!$ вычисляется почленно в собственном базисе), а $G=S/2$ даёт $G-G^*=S$. $\square$

**Предложение 11.2 (параметризация сюръективна).** При фиксированных $\delta$, $s_\pm$ и $U_0$ образ отображения
$(\theta,G)\mapsto X(\theta,G)$ равен $\{X\in\mathcal P_\delta:1\in\operatorname{spec}X_{c_+},\ \ell\in\operatorname{spec}X_{c_-}\}$,
и $\mathcal P_\delta$ есть объединение этих образов по всем парам различных слотов.

*Доказательство.* Образ лежит в этом множестве: собственные значения $X$ равны $\lambda_s\in[\ell,1]$, причём $\lambda_{s_\pm}$ лежат в блоках
$c_\pm$. Обратно, диагонализуем $X_c=W_c\operatorname{diag}(\mu_c)W_c^*$ и переставим компоненты каждого $\mu_c$ (заменив $W_c$ на
$W_c\Pi_c$, где $\Pi_c$ — матрица перестановки) так, чтобы $\mu_{s_+}=1$ и $\mu_{s_-}=\ell$. Это возможно: если $c_+=c_-$
и $\delta>0$, то $1\ne\ell$ — два собственных значения одного блока, и они занимают разные позиции; если $\delta=0$,
все собственные значения равны $1$. Если $\ell<1$, положим $\theta_s=\arcsin\big(2(\mu_s-\ell)/(1-\ell)-1\big)$ для остальных слотов,
так что $\lambda=\mu$ (если $\ell=1$, то $\theta$ любое). Лемма 11.1 даёт $G$ с $U_0e^{G-G^*}=\bigoplus_cW_c$.
*Объединение:* для $X\in\mathcal P_\delta$ выберем блоки $c_\pm$, содержащие собственные значения $1$, $\ell$, и различные слоты
в них; если $c_+=c_-$ и $\delta>0$, то у блока два разных собственных значения, так что $k_{c_+}\ge2$; если $\delta=0$, то $X=1$
и подойдут любые два слота. $\square$

Поскольку каждая пара слотов выбирается с положительной вероятностью, каждая точка оболочки с положительной вероятностью
лежит в образе некоторого старта, а ограничение (в точной арифметике) выполняется на каждой итерации.
`SelfAdjointParam` использует $X=U\operatorname{diag}(\sigma\lambda)U^*$ со знаками $\sigma\in\{\pm1\}^D$, выбираемыми для каждого старта
(в каждом блоке число минусов равномерно распределено на $\{0,\dots,k_c\}$, а их позиции равномерно случайны, так что каждая
расстановка знаков имеет положительную вероятность) и меняемыми только перебросами знаков из раздела 11.5. Тогда
$|X|=U\operatorname{diag}(\lambda)U^*$ (теорема 1.5(e), так как $\lambda>0$), и объединение по знакам и слотам равно
$\{X=X^*:|X|\in\mathcal P_\delta\}$: для такого $X$ диагонализуем $X_c=W_c\operatorname{diag}(\nu_c)W_c^*$, так что
$|X|_c=W_c\operatorname{diag}(|\nu_c|)W_c^*$ с $|\nu_{c,i}|\ge\ell>0$, и применим доказательство предложения 11.2 к
$\mu=|\nu|$ с $\sigma_s=\nu_s/|\nu_s|$. Значение имеет лишь $S=G-G^*$.

Карта $S\mapsto U_0e^S$ ведёт себя хорошо лишь вблизи $S=0$. Пусть $\mathfrak u=\{S\in\mathcal M:S^*=-S\}$ и
$\mathcal M_h=\{H\in\mathcal M:H^*=H\}$ — вещественные подпространства, причём $\mathcal M=\mathfrak u\oplus\mathcal M_h$
($X=\frac{X-X^*}2+\frac{X+X^*}2$). *Вблизи $0$ отображение $S\mapsto U_0e^S$ есть гомеоморфизм открытого шара $N\ni0$ из
$\mathfrak u$ на окрестность $U_0$ в унитарной группе.* Действительно, $\Psi(S,H)=U_0e^Se^H$ — отображение класса $C^1$
$\mathfrak u\times\mathcal M_h\to\mathcal M$ (ряд $\sum_{m,n\ge0}U_0S^mH^n/(m!\,n!)$ и ряды частных производных его членов
по вещественным координатам сходятся равномерно на ограниченных множествах), и
$e^{tS}e^{tH}=1+t(S+H)+O(t^2)$ показывает, что его дифференциал в $0$ равен $(S,H)\mapsto U_0(S+H)$ — биекции на
$\mathcal M$. По теореме об обратной функции существуют открытые шары $N\ni0$ в $\mathfrak u$ и $N_h\ni0$ в
$\mathcal M_h$, такие что $\Psi$ — гомеоморфизм $N\times N_h$ на открытое множество $O\subseteq\mathcal M$. Если
$\Psi(S,H)$ унитарен, то унитарен и $e^H=e^{-S}U_0^*\Psi(S,H)$ (в собственном базисе нормальной матрицы $S$, собственные
значения которой чисто мнимые, $e^{\pm S}$ диагональны с элементами модуля $1$ и взаимно обратны); но $e^H$ имеет
положительные собственные значения $e^{h_i}$, где $h_i$ — собственные значения $H$, так что $e^{h_i}=1$ (следствие 1.4(2)) и $H=0$.
Следовательно, $O\cap U(\mathcal M)=U_0e^N$, где $U(\mathcal M)$ — унитарная группа $\mathcal M$: это окрестность
$U_0$ в $U(\mathcal M)$, и $S\mapsto U_0e^S$ — непрерывная биекция $N\to U_0e^N$, обратное к которой отображение, первая
компонента $\Psi^{-1}$, непрерывно. Вдали от $0$ карта не инъективна: $e^{2\pi iP}=1=e^0$ для любого
проектора $P\in\mathcal M$. Поэтому `extremize` после каждого раунда **перецентрирует** карту ($U_0\leftarrow U_0e^{G-G^*}$,
$G\leftarrow0$). Ещё две карты имеют в $0$ тот же дифференциал $S\mapsto U_0S$. Карта Кэли
$S\mapsto U_0C(-S/2)$ (упражнение 1): $(1+tS/2)(1-tS/2)^{-1}=(1+tS/2)(1+tS/2+O(t^2))=1+tS+O(t^2)$. Полярная ретракция,
$S\mapsto$ унитарный множитель $A=U_0(1+S)$ (теорема 1.8): $A^*A=(1-S)(1+S)=1+S^*S\ge1$ обратим, поэтому $|A|$
обратим и унитарный множитель равен $A|A|^{-1}$; при замене $S$ на $tS$ имеем $|A|^{-1}=(1+t^2S^*S)^{-1/2}=1+O(t^2)$
(в собственном базисе $S^*S$), так что множитель равен $U_0(1+tS)+O(t^2)$.

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

> **Физический смысл** (физическая интерпретация, а не теорема этого курса). Контраст Майкельсона — это видность $(I_{\max}-I_{\min})/(I_{\max}+I_{\min})$
> интерференционных полос. Для $0\le X\le1$ с $\lambda_{\max}=1$ — фильтра, аттенюатора, эффективности детектора —
> он задаёт динамический диапазон $\lambda_{\min}/\lambda_{\max}$, величину, обратную числу обусловленности.

## 11.2 Неравенство фон Неймана для следов

Матрица $(c_{ij})$ с $c_{ij}\ge0$, у которой все суммы по строкам и по столбцам равны $1$, называется **дважды стохастической**; для унитарной $G$
матрицы $(|G_{ij}|^2)$ и $(|G_{ji}|^2)$ дважды стохастические, так как строки и столбцы $G$ — единичные векторы.

**Лемма 11.3 (перестановочное неравенство).** Пусть $p_1\ge\dots\ge p_k$, $q_1\ge\dots\ge q_k$ вещественны, а $(c_{ij})$ дважды
стохастическая. Тогда $\sum_ip_iq_{k+1-i}\le\sum_{i,j}p_iq_jc_{ij}\le\sum_ip_iq_i$.

*Доказательство.* *Верхняя оценка.* Для $m\le k$ положим $r_j=\sum_{i\le m}c_{ij}\in[0,1]$, так что $\sum_jr_j=m$. Поскольку $q_j\ge q_m$ при
$j\le m$ и $q_j\le q_m$ при $j>m$,
$$
\sum_{j\le m}q_j-\sum_jq_jr_j=\sum_{j\le m}q_j(1-r_j)-\sum_{j>m}q_jr_j\ge q_m\Big(m-\sum_jr_j\Big)=0 .
$$
При $x_i=\sum_jq_jc_{ij}$ это означает $\sum_{i\le m}x_i\le\sum_{j\le m}q_j$, причём при $m=k$ имеет место равенство. Преобразованием
Абеля и в силу $p_m-p_{m+1}\ge0$
$$
\sum_ip_ix_i=p_k\sum_{i\le k}x_i+\sum_{m<k}(p_m-p_{m+1})\sum_{i\le m}x_i
\le p_k\sum_{j\le k}q_j+\sum_{m<k}(p_m-p_{m+1})\sum_{j\le m}q_j=\sum_ip_iq_i .
$$
*Нижняя оценка.* Применим верхнюю оценку к $p$, убывающему $q'_j=-q_{k+1-j}$ и дважды стохастической
$c'_{ij}=c_{i,k+1-j}$: $-\sum_{ij}p_iq_jc_{ij}=\sum_{ij}p_iq'_jc'_{ij}\le\sum_ip_iq'_i=-\sum_ip_iq_{k+1-i}$. $\square$

Всякая $A\in M_k(\mathbb C)$ имеет **сингулярное разложение** $A=W\Sigma V^*$, где $W,V$ унитарны,
$\Sigma=\operatorname{diag}(s_i(A))$: запишем $A=W'|A|$ (теорема 1.8), $|A|=V\Sigma V^*$ (теорема 1.3), $W=W'V$.

**Теорема 11.4 (неравенство фон Неймана для следов).** Для $A,B\in M_k(\mathbb C)$
$|\operatorname{Tr}(AB)|\le\sum_is_i(A)\,s_i(B)$.

*Доказательство.* Пусть $A=W_A\Sigma_AV_A^*$, $B=W_B\Sigma_BV_B^*$, $a_i=s_i(A)$, $b_j=s_j(B)$, $G=V_A^*W_B$, $H=V_B^*W_A$.
По цикличности $\operatorname{Tr}(AB)=\operatorname{Tr}(\Sigma_AG\Sigma_BH)=\sum_{i,j}a_iG_{ij}b_jH_{ji}$, поэтому в силу
$|g||h|\le\frac12(|g|^2+|h|^2)$ и $a_ib_j\ge0$
$$
|\operatorname{Tr}(AB)|\le\sum_{i,j}a_ib_j\,d_{ij},\qquad d_{ij}=\tfrac12\big(|G_{ij}|^2+|H_{ji}|^2\big).
$$
Матрица $(d_{ij})$ дважды стохастическая, и лемма 11.3 даёт $\sum_{ij}a_ib_jd_{ij}\le\sum_ia_ib_i$. $\square$

**Следствие 11.5.** (a) Для самосопряжённых $P,Q\in M_k(\mathbb C)$ с собственными значениями $p,q$:
$\sum_ip^\downarrow_iq^\uparrow_i\le\operatorname{Tr}(PQ)\le\sum_ip^\downarrow_iq^\downarrow_i$.
(b) Для $P,Q\ge0$: $\operatorname{Tr}(PQ)\le\|PQ\|_1\le\sum_ip^\downarrow_iq^\downarrow_i$.

*Доказательство.* (a) Как в доказательстве теоремы 3.5, $\operatorname{Tr}(PQ)=\sum_{ij}p_iq_jc_{ij}$ с дважды стохастической
матрицей $c_{ij}=|\langle u_i,v_j\rangle|^2$ для собственных базисов, упорядоченных по убыванию; применим лемму 11.3 (её нижняя оценка равна
$\sum_ip^\downarrow_iq^\downarrow_{k+1-i}=\sum_ip^\downarrow_iq^\uparrow_i$). (b) Первое неравенство — это
теорема 1.9. По теоремам 1.9 и 11.4 $\|PQ\|_1=\max_U|\operatorname{Tr}(P\,QU)|\le\sum_is_i(P)s_i(QU)$, и
$s_i(QU)=q^\downarrow_i$, так как $|QU|=(U^*Q^2U)^{1/2}=U^*QU$ (теорема 1.6). $\square$

Итак, для $X,Y\ge0$ из $M_k(\mathbb C)$ со спектрами $x,y$ нарушение $\|YX\|_1-\operatorname{Tr}(XY)$ из
предложения 3.9 не превосходит $\sum_ix^\downarrow_i(y^\downarrow_i-y^\uparrow_i)$ (следствие 11.5(b) для $(Y,X)$ и
нижняя оценка из (a)). Для эксп. 2 это число и есть точный ответ:

**Предложение 11.6 (эксп. 2 при фиксированных спектрах).** В эксп. 2 $|X|=X_0$, $|Y|=Y_0$ и
$z=\|X_0VY_0\|_1-\operatorname{Tr}(X_0Y_0)$. Если спектры $x_c$ матрицы $(X_0)_c$ и $y_c$ матрицы $(Y_0)_c$ фиксированы, а
собственные базисы и $V$ меняются, то
$$
\sup z=-\inf z=\sum_c\sum_ix^\downarrow_{c,i}\big(y^\downarrow_{c,i}-y^\uparrow_{c,i}\big),
$$
и оба значения достигаются.

*Доказательство.* $X^*X=X_0U^*UX_0=X_0^2$, поэтому $|X|=X_0$ (теорема 1.6); аналогично $|Y|=Y_0$, а $|UA|=|A|$ даёт
$\operatorname{Tr}|XY|=\|X_0VY_0\|_1$. Все переменные действуют поблочно, и $z$ — сумма по блокам, поэтому можно считать $C=1$.
Положим $Y'=VY_0V^*\ge0$ (спектр $y$); поскольку $|AV|=V^*|A|V$, $\|X_0VY_0\|_1=\|X_0Y'V\|_1=\|X_0Y'\|_1$. Следствие 11.5
для пар $(X_0,Y')$ и $(X_0,Y_0)$ вместе с $\|X_0Y'\|_1\ge\operatorname{Tr}(X_0Y')$ даёт
$$
\sum_ix^\downarrow_iy^\uparrow_i-\sum_ix^\downarrow_iy^\downarrow_i\le\|X_0Y'\|_1-\operatorname{Tr}(X_0Y_0)\le\sum_ix^\downarrow_iy^\downarrow_i-\sum_ix^\downarrow_iy^\uparrow_i .
$$
Пусть $J$ — обращающая перестановка, $Je_i=e_{k+1-i}$; тогда $J\operatorname{diag}(w)=\operatorname{diag}(w_k,\dots,w_1)J$, и
следовая норма $\operatorname{diag}(w)J$ равна $\sum_i|w_i|$. При $X_0=\operatorname{diag}(x^\downarrow)$ и $V=J$:
для $Y_0=\operatorname{diag}(y^\uparrow)$ имеем $X_0VY_0=\operatorname{diag}(x^\downarrow_iy^\downarrow_i)J$ и
$\operatorname{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\uparrow_i$ (верхняя оценка); для $Y_0=\operatorname{diag}(y^\downarrow)$ —
$X_0VY_0=\operatorname{diag}(x^\downarrow_iy^\uparrow_i)J$ и $\operatorname{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\downarrow_i$
(нижняя оценка). $\square$

**Следствие 11.7 (огибающая эксп. 2).** Пусть $\ell_X,\ell_Y\in(0,1]$ и некоторое $k_c\ge2$. По всем $X_0,Y_0\ge0$
с $\lambda_{\max}=1$, $\lambda_{\min}(X_0)=\ell_X$, $\lambda_{\min}(Y_0)=\ell_Y$ и всем унитарным операторам
$$
\sup z=-\inf z=\sum_c\Big\lfloor\frac{k_c}2\Big\rfloor(1-\ell_X)(1-\ell_Y).
$$

*Доказательство.* В одном блоке запишем $x_i=x^\downarrow_{c,i}$, $y_i=y^\downarrow_{c,i}$, $k=k_c$, так что $y^\uparrow_{c,i}=y_{k+1-i}$.
Объединяя в пары $i$ и $k+1-i$, получаем $x_i(y_i-y_{k+1-i})+x_{k+1-i}(y_{k+1-i}-y_i)=(x_i-x_{k+1-i})(y_i-y_{k+1-i})$, а при нечётном
$k$ средний индекс даёт вклад $x_i(y_i-y_i)=0$. Следовательно,
$$
\sum_ix_i(y_i-y_{k+1-i})=\sum_{i\le\lfloor k/2\rfloor}(x_i-x_{k+1-i})(y_i-y_{k+1-i})\le\Big\lfloor\frac k2\Big\rfloor(1-\ell_X)(1-\ell_Y),
$$
так как $x_i\in[\ell_X,1]$, $y_i\in[\ell_Y,1]$ убывают, так что оба множителя лежат в $[0,1-\ell_X]$ и $[0,1-\ell_Y]$
соответственно. Равенство достигается, если в каждом блоке $\lfloor k_c/2\rfloor$ наибольших
собственных значений равны $1$, а $\lfloor k_c/2\rfloor$ наименьших равны $\ell_X$ (соответственно $\ell_Y$), причём среднее
(при нечётном $k_c$) произвольно в $[\ell_X,1]$ (соответственно $[\ell_Y,1]$); глобальные экстремумы этих спектров равны $1$ и $\ell_X$ ($\ell_Y$), поскольку некоторое $k_c\ge2$. Предложение 11.6
завершает доказательство. $\square$

Здесь $1-\ell_X=2\delta/(1+\delta)$ при $\Delta(|X|)=\delta$; при нечётном $k_c$ среднее собственное значение роли не играет. Код
проверяет теорему 11.4, следствие 11.5 и экстремальные конфигурации предложения 11.6 для $k=3$, $C=2$.

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

## 11.3 Следовое неравенство с условием равенства

Для $X=X^*$ положим $X_\pm=f_\pm(X)$, где $f_+(t)=\max(t,0)$, $f_-(t)=\max(-t,0)$. Тогда $X_\pm\ge0$, $X=X_+-X_-$,
$|X|=X_++X_-$ и $X_+X_-=0$ (теорема 1.5(c), а также 1.5(e) для $|X|=(X^2)^{1/2}=h(X)$ с
$h(t)=\sqrt{t^2}=|t|=f_+(t)+f_-(t)$; следствие 1.4(3) для положительности).

**Лемма 11.8.** Для $A,B\ge0$ из $M_n(\mathbb C)$ $\|A-B\|_1\le\operatorname{Tr}A+\operatorname{Tr}B$, причём равенство имеет место тогда и только тогда, когда
$AB=0$.

*Доказательство.* Пусть $H=A-B$ и $Q=g(H)$, где $g=1$ на $(0,\infty)$ и $g=0$ в остальных точках; $Q$ — проектор, так как $g=\bar g=g^2$
(теорема 1.5(c)). Так как $|t|=t(2g(t)-1)$, имеем $|H|=H(2Q-1)$ (теорема 1.5(c),(e)), и, используя
$\operatorname{Tr}(AR)=\operatorname{Tr}(RAR)$ для проекторов $R=Q,1-Q$ (цикличность),
$$
\|H\|_1=\operatorname{Tr}\big(A(2Q-1)\big)-\operatorname{Tr}\big(B(2Q-1)\big)
=\operatorname{Tr}A+\operatorname{Tr}B-2\operatorname{Tr}\big((1-Q)A(1-Q)\big)-2\operatorname{Tr}(QBQ).
$$
Вычитаемые следы $\ge0$, а $\operatorname{Tr}(RAR)=\|A^{1/2}R\|_2^2$ обращается в нуль тогда и только тогда, когда $AR=0$; поэтому равенство имеет место
тогда и только тогда, когда $A(1-Q)=0$ и $BQ=0$, и тогда $AB=AQB=A(BQ)^*=0$. Обратно, если $AB=0$, то $BA=0$,
$(A-B)^2=A^2+B^2=(A+B)^2$, и теорема 1.6 даёт $|A-B|=A+B$. $\square$

**Предложение 11.9 (эксп. 3).** Для $X=X^*$ и $Y\ge0$ из $\mathcal M$
$\operatorname{Tr}|YXY|\le\operatorname{Tr}(Y|X|Y)$, причём равенство имеет место тогда и только тогда, когда $YX_+Y^2X_-Y=0$, в частности если $XY=YX$.

*Доказательство.* В $M_D(\mathbb C)$ операторы $A=YX_+Y$ и $B=YX_-Y$ положительны (следствие 1.7 с $T=X_\pm^{1/2}Y$), $YXY=A-B$
и $\operatorname{Tr}(Y|X|Y)=\operatorname{Tr}A+\operatorname{Tr}B$; применим лемму 11.8. Если $XY=YX$, то $Y$ коммутирует с
$X_\pm$ (теорема 1.5(d)), и $AB=Y^4X_+X_-=0$. $\square$

Итак, в эксп. 3 $z\ge0$ и $\underline z=0$ (возьмём $Y=1$). Условие равенства слабее коммутирования: для
$X=\operatorname{diag}(1,-1,0)$ и $Y$ со строками $(1,0,a)$, $(0,1,0)$, $(a,0,1)$, $|a|<1$,
$AB=(Ye_1)\langle Ye_1,Ye_2\rangle(Ye_2)^*=0$, тогда как $XY\ne YX$ при $a\ne0$.

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

## 11.4 Выборка против оптимизации

Насколько велико обычно $|\operatorname{Tr}(AU)|$ по сравнению с его максимумом $\|A\|_1$? Библиотека выбирает $U$ как
множитель $Q$ в QR-разложении $Z=QR$, $R_{ii}>0$, **матрицы Жинибра** $Z$ (независимые элементы с
плотностью $\pi^{-1}e^{-|z|^2}$; рецепт Меццадри), независимо в каждом блоке.

**Лемма 11.10.** $Q$ распределена по Хаару на $U(k)$; в частности, для всякой фиксированной унитарной $V\in M_k(\mathbb C)$ матрица $VQ$ имеет
то же распределение, что и $Q$. Унитарный оператор $\bigoplus_cQ_c$ с независимыми блоками распределён по Хаару на унитарной группе
$\mathcal M$.

*Доказательство.* Это теорема 5.3 (левая инвариантность — это шаг эквивариантности в её доказательстве; распределение
отождествляется с мерой Хаара леммой единственности 5.1, теоретико-мерный шаг которой — следствие A.2) и следствие 5.4. $\square$

**Предложение 11.11.** Для распределённой по Хаару $U\in U(k)$ и $A\in M_k(\mathbb C)$: $\mathbb E\operatorname{Tr}(AU)=0$ и
$\mathbb E|\operatorname{Tr}(AU)|^2=\|A\|_2^2/k$.

*Доказательство.* $\mathbb EU_{ij}=0$ по следствию 5.5. Предложение 5.7 (при $C=1$) для $X=|e_j\rangle\langle e_m|$ даёт
формулу Вайнгартена первого порядка
$\mathbb E[U_{ij}\bar U_{lm}]=\mathbb E(UXU^*)_{il}=\delta_{il}\delta_{jm}/k$; при $k\ge2$ она следует также из
формулы второго порядка, теоремы 5.10, если положить $i_2=i'_2$, $j_2=j'_2$ и просуммировать по обоим индексам, поскольку
$\sum_{i_2,j_2}|U_{i_2j_2}|^2=k$. Так как $\operatorname{Tr}(AU)=\sum_{ij}A_{ji}U_{ij}$,
$\mathbb E|\operatorname{Tr}(AU)|^2=\sum_{ijlm}A_{ji}\bar A_{ml}\mathbb E[U_{ij}\bar U_{lm}]=\frac1k\sum_{ij}|A_{ji}|^2$. $\square$

**Следствие 11.12.** Пусть спектры $X,Y\ge0$ из $M_k(\mathbb C)$ лежат в $[a,b]$, $a>0$, а $U$ распределена по Хаару. Тогда
$$
\frac{\big(\mathbb E|\operatorname{Tr}(XUY)|^2\big)^{1/2}}{\max_V|\operatorname{Tr}(XVY)|}=\frac{\|YX\|_2}{\sqrt k\,\|YX\|_1}\le\min\Big(\frac1{\sqrt k},\frac{b^2}{a^2k}\Big),
\qquad
\mathbb P\big(|\operatorname{Tr}(XUY)|\ge\operatorname{Tr}(XY)\big)\le\frac{b^4}{a^4k^2}.
$$

*Доказательство.* Равенство — это предложение 11.11 для $A=YX$ (так как $\operatorname{Tr}(XUY)=\operatorname{Tr}(AU)$) и
предложение 3.9. Так как $\|A\|_2^2=\sum_is_i^2\le(\sum_is_i)^2=\|A\|_1^2$ (раздел 1.4), отношение $\le1/\sqrt k$.
Поскольку $\|Av\|=\||A|v\|$ (раздел 1.4), диагонализация $|A|$ даёт $s_k(A)=\min_{\|v\|=1}\|Av\|$ и
$s_1(A)=\max_{\|v\|=1}\|Av\|$; а $\|Xv\|^2=\langle v,X^2v\rangle\in[a^2,b^2]$ для единичного $v$ (следствие 1.4(3)
для $X^2-a^2$ и $b^2-X^2$), и аналогично для $Y$. Следовательно, $a^2\le s_k(A)\le s_1(A)\le b^2$, $\|A\|_2^2\le kb^4$ и
$\|A\|_1\ge ka^2$. Наконец, $\operatorname{Tr}(XY)=\operatorname{Tr}(X^{1/2}YX^{1/2})\ge a\operatorname{Tr}X\ge ka^2$, и
неравенство Маркова для $|\operatorname{Tr}(XUY)|^2$ ограничивает вероятность величиной $(kb^4/k)/(ka^2)^2$. $\square$

Итак, $N$ выборок по Хаару находят нарушение неравенства $|\operatorname{Tr}(XUY)|\le\operatorname{Tr}(XY)$ с вероятностью не более
$Nb^4/(a^4k^2)$ (объединение $N$ событий), хотя нарушение существует для каждой некоммутирующей пары (предложение 3.9). Эта
оценка информативна лишь при $k>\sqrt N\,b^2/a^2$ (для $b/a=6$ и $N=2000$: $k>1609$); численно выборки пропускают
нарушение уже при гораздо меньших $k$ (см. ниже). Код проверяет предложение 11.11 методом Монте-Карло (выборочное среднее
отличается от точного значения меньше чем на пять стандартных ошибок — для этого запуска с фиксированным зерном) и для спектров, равномерно распределённых в $[0.2,1.2]$ ($b/a=6$), печатает
отношение из следствия 11.12, $\sup z$ и лучшее $z$ из $2000$ выборок, оба относительно $\operatorname{Tr}(XY)$.

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
ниже доказанной оценки $\min(1/\sqrt k,b^2/(a^2k))$ и примерно уменьшается вдвое при каждом удвоении $k$, т. е. убывает как $1/k$,
а начиная с $k=4$ лучшая из $2000$ выборок не достигает даже $z=0$. Градиентный поиск, напротив, подходит близко к обоим
экстремумам в запуске ниже (это тоже численное наблюдение). Для пары `Xf, Yf` оба экстремума по $U$ известны:

**Лемма 11.13.** Если $k\ge2$, то $\min_U|\operatorname{Tr}(AU)|=0$ для всякой $A\in M_k(\mathbb C)$.

*Доказательство.* Пусть $A=W\Sigma V^*$, а $P$ — циклический сдвиг, $Pe_i=e_{(i\bmod k)+1}$, диагональ которого при $k\ge2$ нулевая. Для
$U=VPW^*$ имеем $\operatorname{Tr}(AU)=\operatorname{Tr}(\Sigma P)=0$. $\square$

Поэтому $\inf_Uz=-\operatorname{Tr}(XY)$ и $\sup_Uz=\|YX\|_1-\operatorname{Tr}(XY)$ (предложение 3.9). Один вызов `extremize`
приближает оба, причём $12$ стартов максимизируют $z$, а $12$ — максимизируют $-z$ (с относительной точностью $10^{-3}$ в этом
запуске, что утверждается ниже). По предложению 11.11 среднеквадратичное значение
$|\operatorname{Tr}(XUY)|$ равно $\|YX\|_2/\sqrt k$; с оценками из доказательства следствия 11.12 (спектры в $[a,b]$)
$\mathbb E|\operatorname{Tr}(XUY)|^2\le b^4$ и $\operatorname{Tr}(XY)\ge ka^2$, так что неравенство Маркова даёт
$\mathbb P\big(z\ge-(1-\varepsilon)\operatorname{Tr}(XY)\big)\le b^4/(\varepsilon^2a^4k^2)$ для $\varepsilon>0$: при больших $k$
величина $z$ с высокой вероятностью близка к своему инфимуму.
*Численное наблюдение* (утверждается для этого запуска с фиксированным зерном): выборка подходит к инфимуму ближе чем на $1\,\%$,
но не достигает даже $z=0$, тогда как супремум положителен:

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

> **Физический смысл** (физическая интерпретация, а не теорема этого курса). Следствие 11.12 — слабый, но доказуемый
> пример того, что статистическая механика называет **типичностью**: в большой размерности почти каждое состояние или эволюция,
> как ожидается, даёт почти одно и то же значение макроскопической величины.
> Точное неравенство определяется тонко настроенными конфигурациями — согласованными собственными базисами, двухуровневыми спектрами, —
> которые случайный унитарный оператор (случайная схема, случайная матрица рассеяния) порождает редко.

## 11.5 Мультистартовый поиск огибающих

`optimize.extremize(objective, params, maximize)` запускает Adam на `rounds` раундов по `steps` шагов (скорость обучения
`lr * decay**round`) с функцией потерь $\mp\sum_b\text{objective}_b$, после каждого раунда перецентрирует унитарные карты и
возвращает **для каждого старта** наибольшее (при `maximize=False` — наименьшее) из значений целевой функции, вычисленных в
градиентных итерациях и один раз в конечной точке. Каждый элемент пакета — независимый старт со своими слотами, знаками, $U_0$ и
контрастом: Adam действует покоординатно, поэтому если $b$-я компонента целевой функции зависит только от $b$-х параметров (как во
всех целевых функциях ниже), старты не взаимодействуют. Знаки
`SelfAdjointParam` дискретны: между раундами `discrete_step` для каждого старта пробует восемь одиночных перебросов знака (каждый
в равномерно случайном слоте), один за другим, и сохраняет каждый переброс, улучшающий целевую функцию.

В точной арифметике каждая итерация допустима (предложение 11.2 и замечания после него), поэтому каждое возвращённое значение
есть значение $z$ в допустимой точке; в арифметике с плавающей точкой — это значение с точностью до округления (утверждения ниже
допускают $10^{-12}$). Следовательно, максимизация даёт
**нижнюю оценку** супремума, минимизация — верхнюю оценку инфимума, и не более того: локальный поиск
может остановиться в локальном экстремуме. Оценка в другую сторону требует доказательства, как следствие 11.7 и предложение 11.9, или
выпуклой релаксации (раздел 11.6). Вычисленная огибающая проверяется (i) по замкнутым формулам; (ii) по
построенным вручную конфигурациям, например прямым суммам решений меньшего размера, ниже которых она опускаться не должна; (iii) по
структурированным выборкам (двухуровневые спектры, базисы-перестановки), значения которых — тоже достигаемые значения.

**Предложение 11.14 (нижняя огибающая эксп. 1).** Пусть все $k_c\ge2$. По $X\in\mathcal P_\delta$, $Y\ge0$ с
$\lambda_{\max}(Y)=1$ и унитарным $U\in\mathcal M$ $\ \min z=-\big(D-2\delta/(1+\delta)\big)$.

*Доказательство.* Из $Y\le1$ следует $X^{1/2}YX^{1/2}\le X$, поэтому $\operatorname{Tr}(XY)\le\operatorname{Tr}X\le(D-1)+\ell(\delta)$, так как одно
собственное значение $X$ равно $\ell(\delta)$, а остальные $\le1$; следовательно, $z\ge-\operatorname{Tr}(XY)\ge-(D-2\delta/(1+\delta))$.
Равенство достигается для $X$ с собственными значениями $\ell(\delta)$ (однократно) и $1$ (остальные), $Y=1$ и $U=\bigoplus_cU_c$ с
$\operatorname{Tr}(X_cU_c)=0$ (лемма 11.13). $\square$

Код проверяет предложение 11.14 на $M_2\oplus M_2$, следствие 11.7 на $M_3\oplus M_3$ (оба контраста заданы,
$k=3$ нечётно) и предложение 11.9; супремумы и инфимумы находятся одним вызовом, как выше. Односторонние утверждения (ни одно
возвращённое значение не выходит за доказанный экстремум) являются следствиями этих теорем; утверждения о том, что поиск подходит
к экстремуму с указанной точностью, — численные наблюдения этого запуска с фиксированным зерном. Ключевое слово `is_self_adjoint` выбирает
алгоритм для следовой нормы (сингулярные числа или собственные значения для $YXY$) без проверки на каждом шаге. Последняя
проверка относится к типу (ii): $z$ из эксп. 3 аддитивна относительно прямых сумм, а контраст прямой суммы равен наибольшему
из контрастов слагаемых (упражнение 3), поэтому для огибающих эксп. 3 с $\Delta(Y)$, свободным в $[0,0.999]$, как в коде,
$$
\overline z_{M_2\oplus M_2}(\delta)\ \ge\ \overline z_{M_2}(\delta)+\sup_{\delta'\le\delta}\overline z_{M_2}(\delta');
$$
в частности, сумма двух значений, достигнутых на $M_2$ при контрастах $\delta$ и $\delta'\le\delta$, является значением,
достигаемым на $M_2\oplus M_2$ при контрасте $\delta$.

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
conj = [(1 - dl) * (2 * 0.999 / 1.999) ** 2 for dl in (0.0, 0.5)]          # unproven closed form, d = 0.999
print(f"Exp. 3  inf {inf3:.1e};  sup for M_2 at delta = 0, 0.5: {s_M2[0]:.5f}, {s_M2[1]:.5f}"
      f"  (closed form {conj[0]:.5f}, {conj[1]:.5f})")
print(f"        sup for M_2 (+) M_2 at 0.5: search {sup_M22:.5f}, direct sum {direct_sum:.5f}")
assert inf3 >= -1e-12 and inf3 < 1e-6                                       # Proposition 11.9
assert sup_M22 >= direct_sum - 5e-3 * direct_sum                           # (ii), up to 0.5 %
```

*Численные наблюдения* (напечатанные значения этого запуска). Поиск на $M_2\oplus M_2$ отстаёт от прямой суммы примерно на
$0.2\,\%$: он сошёлся не полностью, и
наилучшая нижняя оценка для $\overline z_{M_2\oplus M_2}(0.5)$ — большее из двух чисел; обнаруживать такое — назначение проверки
(ii). Значения для $M_2$ совпадают с точностью до $10^{-4}$ с $(1-\delta)\big(2d/(1+d)\big)^2=(1-\delta)(1-\ell(d))^2$, где
$d=0.999$ — наибольший допустимый контраст $Y$. Эта замкнутая формула **не** является теоремой этого курса: сопроводительная
статья наблюдала её численно (в вещественных алгебрах) и не доказала, а код её печатает, ничего о ней не утверждая. Для блоков
большего размера замкнутая форма верхних огибающих эксп. 1 и эксп. 3 нам неизвестна; вычисленные кривые тогда являются нижними
оценками (достигнутыми значениями), и не более того.

## 11.6 С чем сравнивать

**Manopt**, **pymanopt** и **geoopt** предоставляют матричные многообразия (многообразия Штифеля, ортогональные и унитарные группы)
с ретракциями — QR, полярной или Кэли — и римановы решатели (сопряжённые градиенты, доверительные области, риманов
Adam) — альтернативу Adam в перецентрируемой карте (сравнения производительности мы здесь не проводим). **cvxpy** передаёт
выпуклые задачи с операторными переменными SDP-решателям, которые возвращают оптимум вместе с двойственным сертификатом, оба —
с точностью до допуска решателя; релаксация — выпуклая задача на большем допустимом множестве — оценивает супремум сверху, т. е. в
том направлении, которое локальный поиск дать не может. **NCTSSOS**, **ncpol2sdpa** и **NCSOStools** оценивают собственные значения и следы некоммутативных
многочленов сразу во всех размерностях (раздел 0.3). Библиотека добавляет к этому точные параметризации спектральных ограничений
и ограничений на контраст для $\bigoplus_cM_{k_c}$ (экстремумы могут находиться в любом секторе), унитарные операторы, согласованные с
блоками, и одно пакетное вычисление для многих стартов и значений ограничений на CPU или GPU. Она даёт кандидатов в
экстремальные точки и нижние оценки супремумов; остальное должны дать доказательства или релаксации.

## Упражнения

1. (**Преобразование Кэли.**) Для $S^*=-S$ покажите, что $C(S)=(1-S)(1+S)^{-1}$ унитарна и что $S\mapsto C(S)$ —
   биекция на множество унитарных матриц без собственного значения $-1$. Значит, в отличие от экспоненциального отображения (лемма 11.1), карта Кэли
   $U_0C(S)$ не достигает $-U_0$.
2. (**Вещественная мера Хаара.**) Для $O$, распределённой по Хаару на $O(k)$ (QR вещественной матрицы Жинибра), и вещественной $A$ покажите, что
   $\mathbb E\operatorname{Tr}(AO)^2=\|A\|_2^2/k$, и проверьте это на вещественной алгебре.
3. (**Прямые суммы.**) Покажите, что $z$ из эксп. 3 аддитивна относительно прямых сумм и что
   $\Delta(A\oplus A')=\max(\Delta(A),\Delta(A'))$ для $A,A'\ge0$ с $\lambda_{\max}=1$. Выведите неравенство для прямых сумм
   из раздела 11.5 (при $\Delta(Y)$, свободном в фиксированном отрезке $[0,b]$).

<details>
<summary>Решения</summary>

1. $S$ нормальна, а $iS$ самосопряжена, поэтому собственные значения $S$ чисто мнимые (следствие 1.4(2)) и $1\pm S$
   обратимы; они коммутируют, и $(1+S)^*=1-S$. Следовательно,
   $C^*C=(1-S)^{-1}(1+S)(1-S)(1+S)^{-1}=1$, а $1+C=2(1+S)^{-1}$ обратима, так что $-1\notin\operatorname{spec}C$.
   Обратно, для унитарной $U$ с обратимой $1+U$ матрица $S=(1-U)(1+U)^{-1}$ антисамосопряжена:
   $S^*=(1+U^*)^{-1}(1-U^*)=(U+1)^{-1}U\,U^*(U-1)=(U+1)^{-1}(U-1)=-S$, так как $U\pm1$ коммутируют. Для коммутирующих $C$ и $S$
   равенство $C(1+S)=1-S$ равносильно $S(1+C)=1-C$ (оба означают $C+CS+S=1$). Далее, $C(S)$ коммутирует с $S$, а
   $S=(1-U)(1+U)^{-1}$ коммутирует с $U$; поэтому $U=C(S)$ тогда и только тогда, когда $S=(1-U)(1+U)^{-1}$: эти два отображения
   взаимно обратны.

```python
Sk = torch.randn(4, 4, dtype=torch.complex128)
Sk = Sk - Sk.mH                                                   # skew-adjoint
I4 = torch.eye(4, dtype=torch.complex128)
Cay = (I4 - Sk) @ torch.linalg.inv(I4 + Sk)
assert torch.allclose(Cay.mH @ Cay, I4, atol=1e-12)
assert torch.allclose((I4 - Cay) @ torch.linalg.inv(I4 + Cay), Sk, atol=1e-10)
assert torch.linalg.eigvals(Cay).add(1).abs().min() > 1e-6        # -1 is never an eigenvalue
```

2. По теореме B.17 (см. также лемму C.2) QR-множитель $O$ с исправленными знаками вещественной гауссовской матрицы, который
   библиотека выбирает для вещественной алгебры, распределён по Хаару на $O(k)$, поэтому $VO$ имеет то же распределение, что и $O$,
   для $V\in O(k)$.
   Умножение слева на $\operatorname{diag}(\pm1)$ даёт $\mathbb E[O_{ij}O_{lm}]=0$ при $i\ne l$, перестановки показывают,
   что $\mathbb E[O_{ij}O_{im}]$ не зависит от $i$, и $\sum_iO_{ij}O_{im}=\delta_{jm}$; поэтому
   $\mathbb E[O_{ij}O_{lm}]=\delta_{il}\delta_{jm}/k$ и $\mathbb E\operatorname{Tr}(AO)^2=\sum_{ij}A_{ji}^2/k$.

```python
Mr = TypeIAlgebra([5], [5], complex_valued=False, precision='double')
Ar = Mr.operator(torch.randn(1, 1, 5, 5, dtype=torch.float64))
tr2 = (Ar @ Mr.random_unitary_operator(20_000)).trace ** 2
assert abs(tr2.mean().item() - Ar.frobenius_norm().item() ** 2 / 5) < 5 * tr2.std().item() / math.sqrt(20_000)
```

3. Следы и следовые нормы блочно-диагональных операторов — это суммы по блокам, и $|X\oplus X'|=|X|\oplus|X'|$;
   $\lambda_{\min}(A\oplus A')=\min(\lambda_{\min}(A),\lambda_{\min}(A'))$, а $\Delta$ убывает с ростом $\lambda_{\min}$.
   Поэтому если пара $(X,Y)$ допустима при $\delta$, а $(X',Y')$ — при $\delta'\le\delta$, то $(X\oplus X',Y\oplus Y')$
   допустима при $\delta$, и $z=z(X,Y)+z(X',Y')$; осталось взять супремумы по обеим парам и по $\delta'\le\delta$.

</details>

**Далее:** в [главе 12](12_equivariant_learning.md) рассматривается эквивариантное машинное обучение: слои, коммутирующие с симметрией, являются
элементами коммутанта (глава 2) и строятся из тех же пакетных дифференцируемых операторов.
