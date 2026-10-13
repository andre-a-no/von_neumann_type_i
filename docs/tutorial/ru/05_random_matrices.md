[English](../05_random_matrices.md) · [Español](../es/05_random_matrices.md) · [Français](../fr/05_random_matrices.md) · [Deutsch](../de/05_random_matrices.md) · [中文](../zh/05_random_matrices.md) · [日本語](../ja/05_random_matrices.md) · **Русский**

# 5. Случайные матрицы в алгебре

Теория случайных матриц спрашивает, как выглядит «типичный» оператор. В алгебре
$\bigoplus_cM_{k_c}$ ответ должен учитывать секторы: естественные случайные унитарные операторы —
это унитарные элементы самой алгебры, а не большой матричной алгебры $M_D$, а у естественных случайных
гамильтонианов блоки независимы. В этой главе строится мера Хаара на унитарной
группе алгебры, точно вычисляются её первые и вторые моменты, а затем показывается, что
идёт не так в спектральной статистике, если секторы игнорировать.

**Вы узнаете**, почему QR-разложение гауссовой матрицы даёт унитарную матрицу, распределённую по Хаару (и почему
нужно фиксировать фазы), что усреднение по унитарной группе алгебры — это
условное ожидание на центр, формулу Вайнгартена второго порядка, гипотезу Вигнера
и отношение соседних интервалов, а также почему статистика уровней, вычисленная по спектру со смешанными
секторами, выглядит как статистика независимых уровней.

Всюду далее $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$, $D=\sum_ck_c$, $P_c$ —
минимальные центральные проекторы, а $E$ — условное ожидание на центр из
предложения 2.8. $U(k)$ — группа унитарных матриц размера $k\times k$, а
$U(\mathcal M)=\{U\in\mathcal M:U^*U=1\}$ — унитарная группа $\mathcal M$. Случайные величины
заданы на вероятностном пространстве с математическим ожиданием $\mathbb E$; математическое ожидание случайной матрицы
берётся поэлементно. У *стандартной комплексной гауссовой* величины вещественная и мнимая части независимы
и $\sim N(0,\tfrac12)$, так что $\mathbb E|z|^2=1$ (именно такие величины порождает `torch.randn` для комплексного
dtype); *матрица Жинибра* имеет независимые стандартные комплексные гауссовы элементы.

## 5.1 Мера Хаара

Пусть $G$ — компактная группа матриц (замкнутая подгруппа некоторой $U(n)$). Борелевская вероятностная
мера $\mu$ на $G$ называется **левоинвариантной**, если $\mu(gB)=\mu(B)$ для всех $g\in G$ и борелевских множеств
$B$, или, что равносильно, $\int f(gh)\,d\mu(h)=\int f\,d\mu$ для любой ограниченной борелевской функции $f$;
**правоинвариантность** определяется через $Bg$ и $f(hg)$.

**Лемма 5.1 (единственность).** Если $\mu$ — левоинвариантная, а $\nu$ — правоинвариантная борелевская
вероятностная мера на $G$, то $\mu=\nu$. Следовательно, на $G$ существует не более одной левоинвариантной борелевской
вероятностной меры, и если она существует, то она также правоинвариантна и инвариантна относительно $g\mapsto g^{-1}$.

*Доказательство.* Пусть $f:G\to\mathbb R$ непрерывна (а значит, ограничена). По теореме Фубини
$\int\!\big(\int f(gh)\,d\mu(h)\big)d\nu(g)=\int\!\big(\int f(gh)\,d\nu(g)\big)d\mu(h)$. Внутренний
интеграл слева равен $\int f\,d\mu$ по левоинвариантности, внутренний интеграл справа равен $\int f\,d\nu$
по правоинвариантности; поэтому $\int f\,d\mu=\int f\,d\nu$. Две борелевские вероятностные меры на метрическом
пространстве с одинаковыми интегралами от непрерывных функций совпадают: непрерывные функции
$\max(0,1-n\operatorname{dist}(\cdot,K))$ монотонно убывают к индикатору замкнутого множества $K$, так что меры
совпадают на замкнутых множествах (теорема о мажорируемой сходимости), а замкнутые множества порождают борелевские множества и замкнуты относительно
пересечений (теорема Дынкина о $\pi$–$\lambda$-системах). Если $\mu_1,\mu_2$ левоинвариантны, то
$\nu(B)=\mu_2(B^{-1})$ правоинвариантна, поскольку $(Bg)^{-1}=g^{-1}B^{-1}$; поэтому $\mu_1=\nu=\mu_2$,
и $\mu_1$ правоинвариантна и инвариантна относительно обращения. $\square$

**Определение.** **Мерой Хаара** группы $G$ называется её левоинвариантная борелевская вероятностная мера.
Случайный элемент с этим распределением называется **распределённым по Хаару**.

Мера Хаара есть у всякой компактной группы (см. G. B. Folland, *A Course in Abstract Harmonic
Analysis*, 2nd ed., CRC Press, 2016, Section 2.2; приводится без доказательства). Для $U(k)$
эта общая теорема нам не нужна: теорема 5.3 строит меру явно.

**Лемма 5.2 (QR-разложение).** Всякую обратимую $Z\in M_k(\mathbb C)$ можно записать в виде $Z=QR$,
где $Q$ унитарна, а $R$ верхнетреугольная с положительной диагональю. $Q$ и $R$ единственны и
непрерывно зависят от $Z$.

*Доказательство.* Применим процесс Грама–Шмидта к столбцам $z_1,\dots,z_k$ матрицы $Z$, которые линейно независимы:
$q_j=w_j/\|w_j\|$, где $w_j=z_j-\sum_{i<j}\langle q_i,z_j\rangle q_i$; здесь $w_j\neq0$, так как
$z_j\notin\operatorname{span}\{z_1,\dots,z_{j-1}\}=\operatorname{span}\{q_1,\dots,q_{j-1}\}$. Тогда
$z_j=\sum_{i<j}\langle q_i,z_j\rangle q_i+\|w_j\|q_j$, т. е. $Z=QR$ с $R_{jj}=\|w_j\|>0$, и
каждый шаг непрерывен на обратимых матрицах. Если $Q_1R_1=Q_2R_2$, то
$T=Q_2^*Q_1=R_2R_1^{-1}$ унитарна и верхнетреугольна с положительной диагональю. Обратная к ней
$T^*$ нижнетреугольна и одновременно верхнетреугольна (обратная к верхнетреугольной матрице
верхнетреугольна), поэтому $T$ диагональна; диагональная унитарная матрица с положительной диагональю равна $1$. $\square$

**Теорема 5.3 (унитарные матрицы Хаара из гауссовых матриц).** Пусть $Z$ — матрица Жинибра
размера $k\times k$. Тогда $Z$ обратима с вероятностью единица, и множитель $Q$ из леммы 5.2
распределён по Хаару на $U(k)$. Более того, если $Z=Q'R'$ — любое QR-разложение с унитарной $Q'$ и
верхнетреугольной $R'$ (например, возвращаемое LAPACK, у которого диагональ вещественна, но может иметь любой знак), то
$Q=Q'\Lambda$ с $\Lambda=\operatorname{diag}(R'_{ii}/|R'_{ii}|)$ распределена по Хаару.

*Доказательство.* *Обратимость.* $\det Z$ — ненулевой многочлен от $2k^2$ вещественных координат $Z$,
а множество нулей ненулевого вещественного многочлена имеет меру Лебега нуль (индукция по числу
переменных: запишем $p=\sum_jp_j(x')x_n^j$, где некоторый $p_j\neq0$; для $x'$ вне множества нулей
$p_j$, имеющего меру нуль, корней $x_n$ конечное число; завершаем по теореме Фубини). У $Z$ есть плотность, поэтому
$P(\det Z=0)=0$.

*Инвариантность гауссова распределения.* Плотность $Z$ относительно меры Лебега на
$M_k(\mathbb C)\cong\mathbb R^{2k^2}$ равна $\pi^{-k^2}e^{-\|Z\|_2^2}$. Для $V\in U(k)$ отображение
$Z\mapsto VZ$ вещественно-линейно и сохраняет $\|Z\|_2$ (раздел 1.4), поэтому оно является ортогональным отображением
$\mathbb R^{2k^2}$ с якобианом $1$ и сохраняет плотность: $VZ$ имеет то же распределение, что и $Z$.

*Эквивариантность.* Если $Z=QR$, то $VZ=(VQ)R$ — разложение леммы 5.2 для $VZ$, поэтому
$Q(VZ)=VQ(Z)$. Следовательно, для любого борелевского $B\subseteq U(k)$
$P(VQ(Z)\in B)=P(Q(VZ)\in B)=P(Q(Z)\in B)$: распределение $Q(Z)$ левоинвариантно, и по
лемме 5.1 это мера Хаара (заодно доказано, что мера Хаара на $U(k)$ существует).
Наконец, $Q'R'=(Q'\Lambda)(\Lambda^*R')$, а $\Lambda^*R'$ верхнетреугольна с диагональю
$|R'_{ii}|>0$, поэтому $Q'\Lambda=Q(Z)$ в силу единственности. $\square$

То же доказательство для вещественных гауссовых матриц даёт меру Хаара ортогональной группы $O(k)$.
Без поправки $\Lambda$ результат *не* распределён по Хаару: соглашение LAPACK о
знаках $R'_{ii}$ вносит смещение, которое демонстрирует код ниже
(F. Mezzadri, *How to generate random matrices from the classical compact groups*, Notices AMS 54
(2007), 592–604).

**Следствие 5.4 (мера Хаара алгебры).** $U(\mathcal M)=\{\bigoplus_cU_c:U_c\in U(k_c)\}\cong\prod_cU(k_c)$,
и $U=\bigoplus_cU_c$ распределена по Хаару на $U(\mathcal M)$ тогда и только тогда, когда блоки $U_c$ независимы
и каждый $U_c$ распределён по Хаару на $U(k_c)$.

*Доказательство.* $U^*U=\bigoplus_cU_c^*U_c$, поэтому $U$ унитарна тогда и только тогда, когда унитарен каждый блок. Пусть $U_c$
независимы и распределены по Хаару. Для $V=\bigoplus_cV_c\in U(\mathcal M)$ блоки $V_cU_c$ матрицы
$VU$ независимы (как функции независимых величин), и $V_cU_c$ распределён так же, как $U_c$, поэтому
$VU$ распределена так же, как $U$. Значит, распределение $U$ левоинвариантно, и по лемме 5.1 это мера
Хаара компактной группы $U(\mathcal M)$. Обратно, в силу единственности всякая распределённая по Хаару $U$ имеет
это распределение. $\square$

Группа $U(\mathcal M)$ гораздо меньше $U(D)$: при $C\ge2$ унитарная матрица Хаара из $U(D)$ лежит
в $\mathcal M$ с вероятностью нуль (по теореме 5.3 и построению Грама–Шмидта в
лемме 5.2 её первый столбец равен $z/\|z\|$ для гауссова вектора $z$, так что все её элементы почти наверное
ненулевые). Мера Хаара на $SU(k)$ получается в том же духе: если $U$ распределена по Хаару на
$U(k)$ и $\theta=\arg\det U$, то $V=e^{-i\theta/k}U\in SU(k)$, и для $W\in SU(k)$ равенство
$\det(WU)=\det U$ даёт $V(WU)=WV(U)$; поскольку $WU$ распределена так же, как $U$, распределение $V$
левоинвариантно на $SU(k)$, а значит, по лемме 5.1 это мера Хаара на $SU(k)$.

**Следствие 5.5 (первые моменты элементов).** Если $U$ распределена по Хаару на $U(k)$, то
$\mathbb EU_{ij}=0$ и $\mathbb E|U_{ij}|^2=1/k$.

*Доказательство.* $-1\in U(k)$, поэтому $-U$ распределена так же, как $U$, и $\mathbb EU=-\mathbb EU$. Для матрицы
перестановки $P$ матрица $PU$ распределена так же, как $U$, поэтому $\mathbb E|U_{ij}|^2$ не зависит от $i$; так как
$\sum_i|U_{ij}|^2=1$, она равна $1/k$. $\square$

Библиотека порождает распределённые по Хаару унитарные элементы $U(\mathcal M)$ функцией `random_unitary_operator`: одно пакетное
QR-разложение на каждый размер блока сразу для всех секторов и выборок (на GPU — процесс Грама–Шмидта, который
даёт ту же $Q$ в силу единственности в лемме 5.2). Мера `'haar'` — это теорема 5.3 в каждом
блоке; для комплексной алгебры `'haar_su'` домножает на $e^{-i\theta/k}$, как выше; `'coe'` и
`'cse'` — круговые ансамбли Дайсона (первый — предмет упражнения 3). Проверим следствия 5.4 и 5.5,
а также то, что сэмплер библиотеки — это в точности исправленное QR-разложение теоремы 5.3:

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

Без поправки $\mathbb E\operatorname{Re}\operatorname{Tr}Q\approx-0.84$ при $k=2$, тогда как у
унитарной матрицы Хаара $\mathbb EU=0$.

## 5.2 Моменты и исчисление Вайнгартена

**Лемма 5.6.** Если $T\in M_k(\mathbb C)$ коммутирует со всеми диагональными унитарными матрицами и всеми матрицами
перестановок, то $T\in\mathbb C1$. В частности, коммутант $U(k)$ равен $\mathbb C1$.

*Доказательство.* Для $\Delta=\operatorname{diag}(e^{i\theta_1},\dots,e^{i\theta_k})$ с попарно различными
$e^{i\theta_j}$ равенство $\Delta T=T\Delta$ даёт $(e^{i\theta_i}-e^{i\theta_j})T_{ij}=0$, поэтому $T$
диагональна. Для матрицы перестановки $P_\sigma e_j=e_{\sigma(j)}$ равенство $P_\sigma TP_\sigma^{-1}=T$ даёт
$T_{\sigma(j)\sigma(j)}=T_{jj}$, поэтому все диагональные элементы равны. $\square$

**Предложение 5.7 (твирлинг по Хаару — это условное ожидание на центр).** Пусть $U$
распределена по Хаару на $U(\mathcal M)$. Для любого $X\in\mathcal M$
$$
\mathbb E\,[UXU^*]=E(X)=\bigoplus_c\frac{\operatorname{Tr}X_c}{k_c}\,1_{k_c}.
$$

*Доказательство.* $(UXU^*)_c=U_cX_cU_c^*$, и $U_c$ распределена по Хаару на $U(k_c)$ (следствие 5.4). Положим
$T_c=\mathbb E[U_cX_cU_c^*]$ (элементы ограничены). Для $V\in U(k_c)$ матрица $VU_c$ распределена так же, как
$U_c$, поэтому $VT_cV^*=\mathbb E[(VU_c)X_c(VU_c)^*]=T_c$. По лемме 5.6 $T_c=\lambda1$, и
$k_c\lambda=\operatorname{Tr}T_c=\mathbb E\operatorname{Tr}(U_cX_cU_c^*)=\operatorname{Tr}X_c$. $\square$

Нормировка та же, что в предложении 2.8 и в `channels.center_expectation`. Если
$X\in\mathcal M$ коммутирует с $U(\mathcal M)$, то $UXU^*=X$, поэтому $X=E(X)$ централен; усреднение по
Хаару — непрерывный аналог усреднения по группе из упражнения 3 главы 2. Если же $U$
распределена по Хаару на $U(D)$, то то же доказательство (с $C=1$) даёт $\mathbb E[UXU^*]=\frac{\operatorname{Tr}X}{D}1$:
веса секторов теряются, и $UXU^*$ даже не лежит в $\mathcal M$.

> **Физический смысл.** Случайный по Хаару унитарный элемент алгебры — это случайная эволюция, соблюдающая
> правило суперотбора. Усреднение по нему перемешивает всё внутри каждого сектора и сохраняет
> ровно классическую часть: вес каждого заряда. Это «твирлинг» (twirl) квантовой
> теории информации, и это условное ожидание на центр из главы 2.

Для вторых моментов нужен коммутант $\{V\otimes V\}$ на $\mathbb C^k\otimes\mathbb C^k$.
Пусть $F$ — **оператор перестановки** (swap), $F(x\otimes y)=y\otimes x$; он самосопряжён и унитарен, и
$\operatorname{Tr}1=k^2$, $\operatorname{Tr}F=\sum_{i,j}\langle e_i\otimes e_j,e_j\otimes e_i\rangle=k$, $F^2=1$.

**Лемма 5.8.** Пусть $k\ge2$. Матрица $X$ на $\mathbb C^k\otimes\mathbb C^k$ коммутирует с $V\otimes V$
для всех $V\in U(k)$ тогда и только тогда, когда $X=a1+bF$ с $a,b\in\mathbb C$; при этом $1,F$ линейно независимы.

*Доказательство.* $F(V\otimes V)(x\otimes y)=Vy\otimes Vx=(V\otimes V)F(x\otimes y)$, поэтому $1$ и $F$ коммутируют
со всеми $V\otimes V$. Обратно, пусть $X$ коммутирует с ними; запишем
$X_{ij,pq}=\langle e_i\otimes e_j,X(e_p\otimes e_q)\rangle$.
*Диагональные унитарные матрицы:* $\Delta\otimes\Delta$ умножает $e_i\otimes e_j$ на $e^{i(\theta_i+\theta_j)}$,
поэтому $(e^{i(\theta_i+\theta_j)}-e^{i(\theta_p+\theta_q)})X_{ij,pq}=0$. Если мультимножества $\{i,j\}$ и
$\{p,q\}$ различны, то некоторый индекс $m$ встречается в первом на $n\in\{\pm1,\pm2\}$ раз чаще; при
$\theta_m=\pi/2$ и всех остальных $\theta=0$ скобка равна $e^{i(\theta_p+\theta_q)}(e^{in\pi/2}-1)\neq0$,
поэтому $X_{ij,pq}=0$. Следовательно, $X(e_i\otimes e_i)=a_ie_i\otimes e_i$, а при $i\neq j$
$X(e_i\otimes e_j)=b_{ij}\,e_i\otimes e_j+c_{ij}\,e_j\otimes e_i$.
*Перестановки:* применяя $X(P_\sigma\otimes P_\sigma)=(P_\sigma\otimes P_\sigma)X$ к $e_i\otimes e_j$, получаем
$a_{\sigma(i)}=a_i$, $b_{\sigma(i)\sigma(j)}=b_{ij}$, $c_{\sigma(i)\sigma(j)}=c_{ij}$. Перестановки действуют
транзитивно на упорядоченных парах различных индексов, поэтому $a_i=a'$, $b_{ij}=b'$, $c_{ij}=c'$; тем самым
$X=b'1+c'F$ на линейной оболочке $e_i\otimes e_j$ с $i\ne j$.
*Поворот:* пусть $R$ — унитарная матрица с $Re_1=(e_1+e_2)/\sqrt2$, $Re_2=(-e_1+e_2)/\sqrt2$, $Re_m=e_m$ при $m\ge3$.
Тогда $(R\otimes R)(e_1\otimes e_1)=\frac12(e_1\otimes e_1+e_1\otimes e_2+e_2\otimes e_1+e_2\otimes e_2)$, и
$X(R\otimes R)(e_1\otimes e_1)=\frac12\big(a'e_1\otimes e_1+(b'+c')(e_1\otimes e_2+e_2\otimes e_1)+a'e_2\otimes e_2\big)$
должно совпадать с $(R\otimes R)X(e_1\otimes e_1)=a'(R\otimes R)(e_1\otimes e_1)$. Сравнение коэффициентов при
$e_1\otimes e_2$ даёт $a'=b'+c'$, а это и есть значение $b'1+c'F$ на $e_i\otimes e_i$. Итак, $X=b'1+c'F$.
Наконец, $1$ и $F$ переводят $e_1\otimes e_2$ в линейно независимые векторы $e_1\otimes e_2$ и $e_2\otimes e_1$. $\square$

**Следствие 5.9 (двойственность Шура–Вейля второго порядка).** Пусть $k\ge2$ и $P_\pm=\frac12(1\pm F)$.
Это проекторы на симметрическое и антисимметрическое подпространства
$\mathrm{Sym}=\{x:Fx=x\}$ и $\mathrm{Alt}=\{x:Fx=-x\}$ размерностей $k(k+1)/2$ и $k(k-1)/2$.
Представление $V\mapsto V\otimes V$ группы $U(k)$ неприводимо на $\mathrm{Sym}$ и на $\mathrm{Alt}$,
и эти два представления неэквивалентны.

*Доказательство.* Из $F=F^*$ и $F^2=1$ следует, что $P_\pm$ — проекторы с $P_++P_-=1$ и образами $\mathrm{Sym}$,
$\mathrm{Alt}$; их размерности равны $\operatorname{Tr}P_\pm=(k^2\pm k)/2$. Линейная оболочка $\mathcal A$
множества $\{V\otimes V\}$ — $*$-алгебра ($1\otimes1$, произведения $VW\otimes VW$, сопряжённые $V^*\otimes V^*$) с
коммутантом $\mathcal A'=\operatorname{span}\{1,F\}=\operatorname{span}\{P_+,P_-\}$ по лемме 5.8. Пусть
$K\subseteq\mathrm{Sym}$ инвариантно относительно всех $V\otimes V$. По лемме 2.1 проектор на него лежит в $\mathcal A'$, поэтому он
равен $xP_++yP_-$ с $x,y\in\{0,1\}$ (это проектор); он обращается в нуль на $\mathrm{Alt}$, поэтому $y=0$, и
$K=0$ или $K=\mathrm{Sym}$. То же рассуждение применимо к $\mathrm{Alt}$. Если
$T:\mathrm{Sym}\to\mathrm{Alt}$ — сплетающий оператор для этих двух представлений, то $TP_+$ (как
оператор на $\mathbb C^k\otimes\mathbb C^k$) коммутирует с $\mathcal A$, поэтому $TP_+=xP_++yP_-$; он
обращается в нуль на $\mathrm{Alt}$, поэтому $y=0$, а его образ лежит в $\mathrm{Alt}\cap\mathrm{Sym}=0$, поэтому $T=0$. $\square$

**Теорема 5.10 (формула Вайнгартена второго порядка).** Пусть $k\ge2$, $U$ распределена по Хаару на $U(k)$,
и $X$ — матрица на $\mathbb C^k\otimes\mathbb C^k$. Тогда
$$
\mathbb E\,[(U\otimes U)X(U\otimes U)^*]=\alpha\,1+\beta F,\qquad
\alpha=\frac{\operatorname{Tr}X-\operatorname{Tr}(FX)/k}{k^2-1},\quad
\beta=\frac{\operatorname{Tr}(FX)-\operatorname{Tr}X/k}{k^2-1}.
$$
Равносильно, если положить $\mathrm{Wg}(e)=\frac1{k^2-1}$, $\mathrm{Wg}(\tau)=-\frac1{k(k^2-1)}$ на
группе $S_2=\{e,\tau\}$, то
$$
\mathbb E\big[U_{i_1j_1}U_{i_2j_2}\overline{U_{i'_1j'_1}U_{i'_2j'_2}}\big]
=\sum_{\sigma,\pi\in S_2}\delta_{i_1i'_{\sigma(1)}}\delta_{i_2i'_{\sigma(2)}}\,
\delta_{j_1j'_{\pi(1)}}\delta_{j_2j'_{\pi(2)}}\,\mathrm{Wg}(\sigma\pi^{-1}).
$$

*Доказательство.* Пусть $T=\mathbb E[(U\otimes U)X(U\otimes U)^*]$. Как в предложении 5.7, $T$ коммутирует со
всеми $V\otimes V$, поэтому $T=\alpha1+\beta F$ по лемме 5.8. Беря следы, получаем
$\operatorname{Tr}T=\operatorname{Tr}X$, а поскольку $F$ коммутирует с $U\otimes U$,
$\operatorname{Tr}(FT)=\mathbb E\operatorname{Tr}((U\otimes U)^*F(U\otimes U)X)=\operatorname{Tr}(FX)$.
Это даёт $k^2\alpha+k\beta=\operatorname{Tr}X$ и $k\alpha+k^2\beta=\operatorname{Tr}(FX)$ — линейную
систему с определителем $k^2(k^2-1)\neq0$, решение которой указано выше. Для
элементов возьмём $X=|e_{j_1}\otimes e_{j_2}\rangle\langle e_{j'_1}\otimes e_{j'_2}|$: тогда
$\langle e_{i_1}\otimes e_{i_2},T\,e_{i'_1}\otimes e_{i'_2}\rangle$ — левая часть, $\operatorname{Tr}X=\delta_{j_1j'_1}\delta_{j_2j'_2}$
(слагаемое $\pi=e$), $\operatorname{Tr}(FX)=\delta_{j_1j'_2}\delta_{j_2j'_1}$ (слагаемое $\pi=\tau$), а
$1$ и $F$ дают $\delta_{i_1i'_1}\delta_{i_2i'_2}$ (слагаемое $\sigma=e$) и
$\delta_{i_1i'_2}\delta_{i_2i'_1}$ ($\sigma=\tau$). Собирая слагаемые, получаем
$\alpha=\sum_\pi\mathrm{Wg}(e\pi^{-1})[\pi]$ и $\beta=\sum_\pi\mathrm{Wg}(\tau\pi^{-1})[\pi]$,
где $[\pi]$ обозначает произведение $j$-дельт, отвечающее $\pi$. $\square$

Матрица $\big(\mathrm{Wg}(\sigma^{-1}\pi)\big)_{\sigma,\pi}$ — обратная к матрице Грама
$\begin{pmatrix}k^2&k\\k&k^2\end{pmatrix}$ элементов $1$ и $F$ относительно
$\langle A,B\rangle=\operatorname{Tr}(A^*B)$; именно эта форма допускает обобщение.

**Следствие 5.11.** При $k\ge2$:
$$
\mathbb E|U_{11}|^4=\frac2{k(k+1)},\qquad
\mathbb E|U_{11}|^2|U_{12}|^2=\frac1{k(k+1)},\qquad
\mathbb E\big[U_{11}U_{22}\overline{U_{12}U_{21}}\big]=-\frac1{(k-1)k(k+1)}.
$$

*Доказательство.* В теореме 5.10: (i) все индексы равны $1$: каждое произведение дельт равно $1$, поэтому
момент равен $2(\mathrm{Wg}(e)+\mathrm{Wg}(\tau))=2\frac{k-1}{k(k^2-1)}$. (ii) $i=(1,1)$, $j=(1,2)$,
$i'=(1,1)$, $j'=(1,2)$: остаётся только $\pi=e$, что даёт $\mathrm{Wg}(e)+\mathrm{Wg}(\tau)$.
(iii) $i=(1,2)$, $j=(1,2)$, $i'=(1,2)$, $j'=(2,1)$: остаются только $\sigma=e$, $\pi=\tau$,
что даёт $\mathrm{Wg}(\tau)$. $\square$

Последний момент отрицателен: элементы унитарной матрицы Хаара не независимы, и
корреляции имеют относительный порядок $1/k$. Для алгебры элементы разных блоков
независимы (следствие 5.4), поэтому смешанные моменты по разным секторам распадаются в произведения моментов отдельных блоков.

Для моментов любого порядка годится то же рассуждение, как только известен коммутант
$\{V^{\otimes n}\}$; он натянут на перестановки $n$ тензорных сомножителей
(двойственность Шура–Вейля, см. R. Goodman, N. R. Wallach, *Symmetry, Representations, and Invariants*,
Springer, 2009, Chapter 9). Мы приводим результат без доказательства.

**Теорема 5.12 (формула Вайнгартена; без доказательства).** Пусть $n\ge1$ и $U$ распределена по Хаару на $U(k)$. Тогда
$$
\mathbb E\big[U_{i_1j_1}\cdots U_{i_nj_n}\overline{U_{i'_1j'_1}\cdots U_{i'_nj'_n}}\big]
=\sum_{\sigma,\pi\in S_n}\prod_{m=1}^n\delta_{i_mi'_{\sigma(m)}}\delta_{j_mj'_{\pi(m)}}\;\mathrm{Wg}(\sigma\pi^{-1},k),
$$
где $\mathrm{Wg}(\cdot,k)$ — функция на $S_n$, зависящая только от цикловой структуры. При $n\le k$
матрица $\big(\mathrm{Wg}(\sigma^{-1}\pi,k)\big)_{\sigma,\pi}$ обратна к матрице Грама
$\big(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)}\big)_{\sigma,\pi}$; при $n>k$ эта матрица вырождена,
и берётся её псевдообратная. Моменты с разным числом $U$ и $\overline U$ равны нулю. (D. Weingarten,
J. Math. Phys. 19 (1978), 999–1001; B. Collins, Int. Math. Res. Not. 2003, 953–982; B. Collins,
P. Śniady, Comm. Math. Phys. 264 (2006), 773–795.)

Проверим предложение 5.7 сравнением с `channels.center_expectation`, коммутант из леммы 5.8 — точной
линейной алгеброй, а теорему 5.10 и следствие 5.11 — методом Монте-Карло на блоке $U(3)$
из полученных выше выборок:

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

## 5.3 Статистика уровней

Пусть $H$ — случайная самосопряжённая матрица размера $n\times n$ с собственными значениями
$\lambda_1\le\dots\le\lambda_n$. Её **интервалы** (расстояния между соседними уровнями) — это $s_i=\lambda_{i+1}-\lambda_i$. Стандартных ансамблей
два: **GOE**, $H=(G+G^T)/2$, где $G$ — вещественная матрица с независимыми элементами $N(0,1)$, и
**GUE**, $H=(G+G^*)/\sqrt2$, где $G$ — матрица Жинибра (раздел 1.7). В обоих случаях диагональные элементы
распределены по $N(0,1)$, наддиагональные элементы — по $N(0,\tfrac12)$ в GOE, а в GUE их вещественная и мнимая части независимы
и распределены по $N(0,\tfrac12)$; наконец,
плотность $H$ пропорциональна $e^{-\operatorname{Tr}H^2/2}$.

**Предложение 5.13 (гипотеза Вигнера).** Пусть $n=2$, $s=\lambda_2-\lambda_1$ и $S=s/\mathbb Es$.
Для GOE ($\beta=1$) и GUE ($\beta=2$) величина $S$ имеет плотности
$$
p_1(S)=\frac\pi2\,S\,e^{-\pi S^2/4},\qquad p_2(S)=\frac{32}{\pi^2}\,S^2\,e^{-4S^2/\pi},
$$
с функциями распределения $P(S\le t)=1-e^{-\pi t^2/4}$ и
$P(S\le t)=\operatorname{erf}(2t/\sqrt\pi)-\frac{4t}\pi e^{-4t^2/\pi}$. В частности,
$P(S\le\varepsilon)$ имеет порядок $\varepsilon^{\beta+1}$ при $\varepsilon\to0$.

*Доказательство.* Для $H=\begin{pmatrix}a&b\\\bar b&c\end{pmatrix}$ собственные значения равны
$\frac{a+c}2\pm\frac12\sqrt{(a-c)^2+4|b|^2}$, поэтому $s=\sqrt{(a-c)^2+4|b|^2}$. *GOE:* $x=a-c$ и
$y=2b$ — независимые величины $N(0,2)$, и в полярных координатах
$P(s\le t)=\int_{x^2+y^2\le t^2}\frac{e^{-(x^2+y^2)/4}}{4\pi}\,dx\,dy=\int_0^t\frac r2e^{-r^2/4}dr=1-e^{-t^2/4}$.
Тогда $\mathbb Es=\int_0^\infty P(s>t)\,dt=\int_0^\infty e^{-t^2/4}dt=\sqrt\pi$ и
$P(S\le t)=P(s\le\sqrt\pi t)=1-e^{-\pi t^2/4}$; остаётся продифференцировать. *GUE:* $x=a-c$, $y=2\operatorname{Re}b$,
$y'=2\operatorname{Im}b$ — независимые $N(0,2)$, и в сферических координатах плотность
$s=\sqrt{x^2+y^2+y'^2}$ равна $4\pi r^2(4\pi)^{-3/2}e^{-r^2/4}=\frac{r^2e^{-r^2/4}}{2\sqrt\pi}$. Её среднее равно
$\frac1{2\sqrt\pi}\int_0^\infty r^3e^{-r^2/4}dr=\frac{8}{2\sqrt\pi}=\frac4{\sqrt\pi}$. Подстановка
$r=4S/\sqrt\pi$ даёт $p_2$, а интегрирование $p_2$ по частям
($\int x^2e^{-ax^2}dx=-\frac{x}{2a}e^{-ax^2}+\frac1{2a}\int e^{-ax^2}dx$ с $a=4/\pi$) даёт
функцию распределения. Поведение при малых $\varepsilon$ следует из $p_\beta(S)\sim cS^\beta$. $\square$

Показатель считает размерности: вещественная симметрическая матрица $2\times2$ вырождена тогда и только тогда, когда выполнены два вещественных
условия ($a=c$, $b=0$), а комплексная эрмитова — когда выполнены три. При любом $n$ совместная
плотность собственных значений пропорциональна $\prod_{i<j}|\lambda_i-\lambda_j|^\beta e^{-\sum_i\lambda_i^2/2}$, и
при $n\to\infty$ распределение интервалов в объёме спектра сходится к пределу (Годена–Меты), который близок
к гипотезе Вигнера, но не совпадает с ней (M. L. Mehta, *Random Matrices*, 3rd ed., Elsevier, 2004; приводится без
доказательства).

> **Физический смысл.** Вигнер предположил, что энергетические уровни сложных квантовых систем (тяжёлых
> ядер) статистически устроены как собственные значения случайного гамильтониана с нужной симметрией:
> вещественного симметрического при наличии симметрии относительно обращения времени ($\beta=1$) и комплексного эрмитова в противном случае
> ($\beta=2$). Тогда уровни *отталкиваются*: два уровня одной симметрии почти никогда не сближаются
> (фон Нейман и Вигнер, 1929). Отталкивание — типичный признак квантового хаоса, независимые
> уровни (без отталкивания) — признак интегрируемых или локализованных систем (гипотезы Богигаса–Джаннони–Шмита и
> Берри–Табора).

Интервалы зависят от локальной плотности уровней, которая меняется вдоль спектра; чтобы сравнивать их
с $p_\beta$, нужна **развёртка** (unfolding) — замена переменных $\lambda\mapsto\bar N(\lambda)$ с помощью
гладкого приближения $\bar N$ считающей функции, после которой средний интервал равен $1$. Выбор
$\bar N$ — источник смещения. Отношение соседних интервалов позволяет его избежать.

**Определение.** **Отношение соседних интервалов** — это $r_i=\min(s_i,s_{i+1})/\max(s_i,s_{i+1})\in[0,1]$
(V. Oganesyan, D. A. Huse, Phys. Rev. B 75 (2007), 155111). Оно не меняется при
$\lambda\mapsto a\lambda+b$, $a>0$, поэтому не требует развёртки там, где плотность примерно постоянна
на протяжении двух интервалов.

**Предложение 5.14 (независимые уровни).** Пусть $s,s'$ независимы и имеют одно и то же экспоненциальное
распределение (как два соседних интервала пуассоновского процесса). Тогда величина $r=\min(s,s')/\max(s,s')$ имеет
на $[0,1]$ функцию распределения $P(r\le t)=\frac{2t}{1+t}$ и плотность $\frac2{(1+t)^2}$, и
$$
\mathbb Er=2\ln2-1\approx0.3863 .
$$

*Доказательство.* Масштабированием сведём к интенсивности $1$. Для $x=s/s'$ имеем $P(x\le t)=\mathbb E\,P(s\le ts'\mid s')=\mathbb E(1-e^{-ts'})=1-\frac1{1+t}=\frac t{1+t}$.
При $t<1$ события $\{x\le t\}$ и $\{1/x\le t\}$ не пересекаются, $x$ и $1/x$ одинаково распределены
($s,s'$ перестановочны), а $\{r\le t\}$ — их объединение, поэтому $P(r\le t)=\frac{2t}{1+t}$. Тогда
$\mathbb Er=\int_0^1\frac{2t}{(1+t)^2}dt=2\Big[\ln(1+t)+\frac1{1+t}\Big]_0^1=2\ln2-1$. $\square$

То же самое точно верно для $n$ независимых равномерно распределённых уровней на отрезке, поскольку их интервалы
распределены как $(e_1,\dots,e_{n+1})/\sum_je_j$ с независимыми экспоненциальными $e_j$, так что каждое
отношение $s_i/s_{i+1}$ распределено как $e_i/e_{i+1}$ (W. Feller, *An Introduction to Probability Theory and Its
Applications*, Vol. II, Wiley, 1971, Chapter I; приводится без доказательства). Для GOE и GUE в
пределе $n\to\infty$ $\langle r\rangle=0.5307$ и $0.5996$ (численные значения; аналоги гипотезы Вигнера для матриц $3\times3$
дают $4-2\sqrt3\approx0.5359$ и $2\sqrt3/\pi-\frac12\approx0.6027$; Y. Y. Atas,
E. Bogomolny, O. Giraud, G. Roux, Phys. Rev. Lett. 110 (2013), 084101; приводится без доказательства).

Функция `SpinChain.level_spacing_ratio` (в `chains`) возвращает среднее отношение для каждой строки
уровней. Проверим предложения 5.13 и 5.14 и значение для GOE:

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
mc(ratio(levels), 2 * math.log(2) - 1, "<r>, independent uniform levels")

n = 64
alg64 = TypeIAlgebra([n], [n], complex_valued=False, precision='double')
G = torch.randn(800, 1, n, n, dtype=torch.float64)
lam = alg64.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()[0]
mc(ratio(lam[:, n // 4: 3 * n // 4]), 0.5307, "<r>, GOE(64), middle half of the spectrum")
```

## 5.4 Почему секторы нельзя смешивать

Пусть $H=\bigoplus_cH_c\in\mathcal M$ — случайный самосопряжённый элемент с **независимыми блоками**,
например блоками GOE. Как оператор на $\mathbb C^D$ он имеет спектр, равный объединению спектров
блоков. Внутри блока уровни отталкиваются; покажем, что уровни разных блоков — нет.

Предположим, что у каждого блока есть ограниченная непрерывная **плотность уровней** $\rho_c$:
$\mathbb E\sum_i\varphi(\lambda_i(H_c))=\int\varphi\,\rho_c$ для ограниченных борелевских $\varphi$ (так что
$\int\rho_c=k_c$). Для GOE и GUE это следует из приведённой выше совместной плотности.

**Теорема 5.15 (отсутствие отталкивания между секторами).** Для $c\neq d$ и любой ограниченной борелевской
$\psi:\mathbb R^2\to\mathbb R$
$$
\mathbb E\sum_{i,j}\psi\big(\lambda_i(H_c),\lambda_j(H_d)\big)=\iint\psi(x,y)\,\rho_c(x)\rho_d(y)\,dx\,dy .
$$
Следовательно, если $N_\varepsilon$ — число пар уровней из разных блоков на расстоянии
$<\varepsilon$, то $\mathbb EN_\varepsilon/(2\varepsilon)\to\sum_{c<d}\int\rho_c\rho_d$ при
$\varepsilon\to0$, а наименьший интервал $\delta$ объединённого спектра удовлетворяет неравенству
$$
P(\delta<\varepsilon)\ \ge\ \frac{\mathbb EN_\varepsilon}{\sum_{c<d}k_ck_d}
=\frac{2\varepsilon\sum_{c<d}\int\rho_c\rho_d}{\sum_{c<d}k_ck_d}+o(\varepsilon).
$$

*Доказательство.* Упорядоченные собственные значения — непрерывные функции матрицы, поэтому
$\Lambda_c=(\lambda_i(H_c))_i$ и $\Lambda_d$ — независимые случайные векторы. По теореме Фубини (с условием на
$\Lambda_d$) $\mathbb E\sum_{i,j}\psi(\Lambda_{c,i},\Lambda_{d,j})=\mathbb E\sum_j g(\Lambda_{d,j})$, где
$g(y)=\mathbb E\sum_i\psi(\Lambda_{c,i},y)=\int\psi(x,y)\rho_c(x)\,dx$ — ограниченная борелевская функция; остаётся применить
определение $\rho_d$ к $g$. При $\psi(x,y)=1\{|x-y|<\varepsilon\}$
$$
\frac{\mathbb EN_\varepsilon}{2\varepsilon}=\sum_{c<d}\int\rho_c(x)\Big(\frac1{2\varepsilon}\int_{x-\varepsilon}^{x+\varepsilon}\rho_d(y)\,dy\Big)dx ;
$$
выражение в скобках стремится к $\rho_d(x)$ в силу непрерывности и ограничено величиной $\sup\rho_d$, а $\rho_c$ интегрируема,
поэтому теорема о мажорируемой сходимости даёт предел. Наконец, $N_\varepsilon\le\sum_{c<d}k_ck_d$, поэтому
$\mathbb EN_\varepsilon\le\big(\sum_{c<d}k_ck_d\big)P(N_\varepsilon\ge1)$, а из $N_\varepsilon\ge1$ следует $\delta<\varepsilon$. $\square$

Итак, для двух блоков с перекрывающимися спектрами $P(\delta<\varepsilon)$ имеет порядок не меньше $\varepsilon$,
тогда как внутри блока GOE размера $2\times2$ она равна $1-e^{-\varepsilon^2/4}\le\varepsilon^2/4$
(предложение 5.13, ненормированный интервал), а для блока GOE размера $n\times n$ имеет порядок
$\varepsilon^2$ (из совместной плотности; здесь не доказывается). Тест на отталкивание уровней не проходит для
объединённого спектра, хотя каждый блок — идеальный GOE.

*Много секторов (вывод, а не доказательство).* Для стационарных спектров вероятность $\mathcal G(s)$
не найти ни одного уровня в интервале длины $s$ определяет плотность интервалов: $p(s)=\mathcal G''(s)$ при
единичном среднем интервале (Mehta, *op. cit.*). При наложении независимых спектров эти вероятности
пустых интервалов перемножаются. Если каждый из $C$ секторов даёт долю $1/C$ уровней и имеет
вероятность пустого интервала $\mathcal G_1$ в собственных развёрнутых единицах, то $\mathcal G(s)=\mathcal G_1(s/C)^C$, где
$\mathcal G_1(x)=1-x+\int_0^x(x-t)p_1(t)\,dt$. Дифференцируя дважды в $0$, получаем
$$
p(0)=1-\frac1C+\frac{p_1(0)}C=1-\frac1C ,
$$
аналог теоремы 5.15 для развёрнутых спектров, и $\mathcal G_1(s/C)^C\to e^{-s}$ — вероятность пустого интервала
для пуассоновского процесса: при $C\to\infty$ объединённый спектр выглядит как набор независимых уровней
(N. Rosenzweig, C. E. Porter, Phys. Rev. 120 (1960), 1698–1714; M. V. Berry, M. Robnik, J. Phys. A 17
(1984), 2413–2421; строгие варианты пуассоновского предела — это предельные теоремы для суперпозиций
разреженных точечных процессов, см. B. Grigelionis, Theory Probab. Appl. 8 (1963)).

В библиотеке блоки разделены с самого начала: `eigenvalues()` возвращает по одному спектру на
сектор, а смешать их можно только явным `torch.cat`. Поместим $C$ независимых блоков GOE размера $64$
в одну алгебру и сравним среднее отношение интервалов по секторам с его значением для объединённого спектра:

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
    assert abs(per_sector.mean() - 0.5307) < 0.01                    # every sector is GOE
assert merged_r[1] > merged_r[2] > merged_r[4] > merged_r[16]
assert abs(merged_r[16] - (2 * math.log(2) - 1)) < 0.01              # close to independent levels
```

В каждом секторе $\langle r\rangle\approx0.53$; для объединённого спектра значение падает до $0.42$ при двух секторах и
до $0.39$ — значения для независимых уровней — при $8$ и $16$.

То же происходит и в физической модели. Цепочка Гейзенберга со случайными полями (глава 9) сохраняет
намагниченность, поэтому её гамильтониан лежит в $\bigoplus_{N=0}^LM_{\binom LN}$; при слабом беспорядке
уровни внутри каждого сектора имеют статистику GOE (численно; A. Pal, D. A. Huse, Phys. Rev. B 82 (2010),
174411). Её блоки *не* независимы — все они зависят от одних и тех же случайных полей, — поэтому
теорема 5.15 в данной формулировке неприменима; эксперимент показывает тот же эффект.

```python
chain = SpinChain(10, complex_valued=False)                         # sectors N = 0, ..., 10
lam = chain.random_field_heisenberg(W=0.5, batch_size=20).eigenvalues()
k5 = chain.sector_dim(5)                                            # 252
r_sector = ratio(lam[5][:, k5 // 4: 3 * k5 // 4]).mean().item()
spectrum = torch.cat(lam, dim=-1).sort(dim=-1).values               # all 1024 levels mixed
r_all = ratio(spectrum[:, 256:768]).mean().item()
print(f"Heisenberg chain, L = 10: <r> in sector N = 5: {r_sector:.4f}, all sectors mixed: {r_all:.4f}")
assert r_sector > 0.49 and r_all < 0.42
```

Внутри сектора $N=5$ отношение равно $0.52$ (значение для конечного размера, немного меньше $0.5307$);
спектр со смешанными секторами даёт $0.39$.

> **Физический смысл.** Уровни с разными сохраняющимися квантовыми числами не взаимодействуют и поэтому
> свободно пересекаются; отталкиваются только уровни одной симметрии. Поэтому спектр, смешивающий секторы симметрии,
> показывает пуассоноподобную статистику, даже если динамика в каждом секторе полностью хаотична,
> и неразрешённую симметрию можно ошибочно принять за интегрируемость или локализацию. Профессиональное
> правило — разрешить все симметрии, прежде чем вычислять статистику уровней, — и есть утверждение, что
> вычислять следует в алгебре $\bigoplus_cM_{k_c}$, а не в $M_D$.

## 5.5 С чем сравнивать

Унитарные матрицы Хаара: `scipy.stats.unitary_group` (исправленное QR-разложение теоремы 5.3) и `ortho_group`,
случайные по Хаару унитарные матрицы QuTiP, а также любое QR-разложение `numpy`/`torch` с последующей поправкой фаз.
Исчисление Вайнгартена: символьные пакеты IntU (Mathematica; Z. Puchała, J. A. Miszczak), RTNI
(Mathematica и Python; M. Fukuda, R. König, I. Nechita, J. Phys. A 52 (2019), 425303) и
пакет Python Haarpy символьно вычисляют моменты любого порядка; формулы, доказанные в этой
главе, — их случай второго порядка. Статистику уровней обычно вычисляют несколькими строками `numpy`,
а программы точной диагонализации, такие как QuSpin, выдают спектры по секторам. Библиотека
добавляет унитарную группу алгебры как единый объект: независимые блоки Хаара всех размеров, порождаемые
пакетными вызовами, твирлинг, отождествлённый с условным ожиданием
`channels.center_expectation`, и спектры, возвращаемые по секторам, так что смешивание секторов никогда не
происходит по умолчанию.

## Упражнения

1. Пусть $U$ распределена по Хаару на $U(\mathcal M)$ и $\operatorname{Tr}U=\sum_c\operatorname{Tr}U_c$. Покажите, что
   $\mathbb E|\operatorname{Tr}U|^2=C$ и $\mathbb E|\operatorname{Tr}U|^4=\sum_cm(k_c)+2C(C-1)$, где
   $m(1)=1$ и $m(k)=2$ при $k\ge2$; сравните с $1$ и $2$ для унитарной матрицы Хаара из $U(D)$, $D\ge2$. Проверьте
   для $\mathcal M=M_1\oplus M_2\oplus M_3$.
2. С помощью теоремы 5.12 покажите, что $\mathbb E|U_{11}|^{2n}=\frac{n!}{k(k+1)\cdots(k+n-1)}$ при $n\le k$.
   Вычислите численно матрицу Вайнгартена порядка $3$ при $k=3$ и проверьте, что $\mathbb E|U_{11}|^6=\frac1{10}$.
3. **Круговой ортогональный ансамбль** (`measure='coe'`) — это $S=W^TW$, где $W$ распределена по Хаару на $U(k)$. Покажите,
   что $S$ унитарна и симметрична, что её распределение инвариантно относительно $S\mapsto V^TSV$ для любой
   $V\in U(k)$ и что средняя плотность её собственных фаз равномерна. Проверьте численно, что для её
   собственных фаз $\langle r\rangle\approx0.53$, как для GOE.

<details>
<summary>Решения</summary>

1. Обозначим $t_c=\operatorname{Tr}U_c$; величины $t_c$ независимы (следствие 5.4). Умножение $U_c$ на
   $e^{i\theta}$ не меняет его распределения и умножает $t_c$ на $e^{i\theta}$; поэтому момент
   $\mathbb E[t_at_b\bar t_c\bar t_d]$ равен нулю, если только каждый сектор не встречается в $(a,b)$ столько же раз, сколько в
   $(c,d)$, т. е. $\{a,b\}=\{c,d\}$ как мультимножества; аналогично $\mathbb E[t_a\bar t_b]=0$ при $a\neq b$. По
   предложению 5.7 $\mathbb E[(U_c)_{ii}\overline{(U_c)_{jj}}]=\langle e_i,\mathbb E[U_ce_ie_j^*U_c^*]e_j\rangle=\delta_{ij}/k_c$,
   поэтому $\mathbb E|t_c|^2=1$ и $\mathbb E|\operatorname{Tr}U|^2=C$. В $\mathbb E|\sum_at_a|^4$ слагаемые с
   $a=b$ дают $\sum_a\mathbb E|t_a|^4$, слагаемые с $a\neq b$ дают
   $2\sum_{a\neq b}\mathbb E|t_a|^2\,\mathbb E|t_b|^2=2C(C-1)$ (два способа выбрать $(c,d)$). При $k=1$ $|t|=1$.
   При $k\ge2$ $\mathbb E|t|^4=\sum_{i,j,i',j'}\mathbb E[U_{ii}U_{jj}\overline{U_{i'i'}U_{j'j'}}]$; в
   теореме 5.10 два слагаемых с $\sigma=\pi$ одним и тем же способом требуют $\{i',j'\}=\{i,j\}$ и после суммирования дают
   по $k^2\mathrm{Wg}(e)$, а два слагаемых с $\sigma\neq\pi$ требуют $i=j=i'=j'$ и дают
   по $k\,\mathrm{Wg}(\tau)$: всего $2\big(\frac{k^2}{k^2-1}-\frac1{k^2-1}\big)=2$. Случай $C=1$,
   $k=D$ даёт $1$ и $2$ для $U(D)$. Для $M_1\oplus M_2\oplus M_3$:
   $\mathbb E|\operatorname{Tr}U|^4=1+2+2+12=17$ против $2$ для унитарной матрицы Хаара из $U(6)$.

```python
tr = U.trace                                                        # Tr U = sum_c Tr U_c, shape (B,)
mc(tr.abs() ** 2, 3.0, "E|Tr U|^2")
mc(tr.abs() ** 4, 17.0, "E|Tr U|^4")
```

2. Когда все индексы равны $1$, каждая дельта равна $1$, поэтому
   $\mathbb E|U_{11}|^{2n}=\sum_{\sigma,\pi}\mathrm{Wg}(\sigma\pi^{-1})=\mathbf 1^TW\mathbf 1$, где
   $W=\Gamma^{-1}=(\mathrm{Wg}(\sigma^{-1}\pi))$ и $\Gamma=(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)})$. Здесь мы воспользовались тем, что
   $\mathrm{Wg}$ зависит только от цикловой структуры, и тем, что $\sigma\pi^{-1}=\sigma(\pi^{-1}\sigma)\sigma^{-1}$
   сопряжена с $\pi^{-1}\sigma=(\sigma^{-1}\pi)^{-1}$, которая имеет ту же цикловую структуру, что и $\sigma^{-1}\pi$. Каждая
   строчная сумма $\Gamma$ равна $c_n(k)=\sum_{\rho\in S_n}k^{\#\mathrm{cycles}(\rho)}$, поэтому $\Gamma\mathbf 1=c_n(k)\mathbf 1$,
   $W\mathbf 1=\mathbf 1/c_n(k)$ и $\mathbf 1^TW\mathbf 1=n!/c_n(k)$. Наконец,
   $c_n(x)=x(x+1)\cdots(x+n-1)$ по индукции: перестановка $\{1,\dots,n\}$ получается из перестановки
   $\{1,\dots,n-1\}$ либо добавлением $n$ как неподвижной точки (на один цикл больше, множитель $x$), либо вставкой $n$
   после одного из $n-1$ элементов в его цикле ($n-1$ способов, число циклов то же); это биекция, поэтому
   $c_n(x)=(x+n-1)c_{n-1}(x)$.

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

3. $S^*S=W^*\bar WW^TW=1$, поскольку $W^T$ унитарна, и $S^T=W^TW=S$. Для $V\in U(k)$ матрица $WV$ распределена
   так же, как $W$ (правоинвариантность, лемма 5.1), и $(WV)^T(WV)=V^TSV$. При $V=e^{i\theta/2}1$ матрица $e^{i\theta}S$
   распределена так же, как $S$, поэтому распределение множества собственных фаз инвариантно относительно поворотов, а значит, инвариантна и
   средняя плотность, которая поэтому равномерна. В частности, интервалы не требуют развёртки. Локальная
   статистика COE совпадает со статистикой GOE при больших $k$ (F. J. Dyson, J. Math. Phys. 3
   (1962), 140–156; Mehta, *op. cit.*; приводится без доказательства).

```python
alg_s = TypeIAlgebra([32], [32], precision='double')
S = alg_s.random_unitary_operator(400, measure='coe').matrix[:, 0]
assert close(S, S.mT) and close(S @ S.mH, torch.eye(32, dtype=dt))
phases = torch.linalg.eigvals(S).angle().sort(dim=-1).values
s = torch.diff(torch.cat([phases, phases[:, :1] + 2 * math.pi], dim=-1), dim=-1)   # around the circle
r = torch.minimum(s, s.roll(-1, -1)) / torch.maximum(s, s.roll(-1, -1))
mc(r.mean(-1), 0.5307, "<r>, COE(32) eigenphases")
```

</details>

**Далее:** в [главе 6](06_quantum_channels.md) рассматриваются квантовые каналы: представления Крауса и Стайнспринга,
сопряжённые относительно следа и отображения между секторами. Условное
ожидание $E$, встретившееся здесь как среднее по Хаару, — первый пример.
