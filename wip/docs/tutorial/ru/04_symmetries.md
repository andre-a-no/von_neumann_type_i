[English](../04_symmetries.md) · [Español](../es/04_symmetries.md) · [Français](../fr/04_symmetries.md) · [Deutsch](../de/04_symmetries.md) · [中文](../zh/04_symmetries.md) · [日本語](../ja/04_symmetries.md) · **Русский**

# 4. Симметрии и представления

В главе 2 было показано, что группа симметрий порождает алгебру фон Неймана — линейную оболочку
операторов симметрии — и её коммутант, симметричные наблюдаемые (предложение 2.7), и один пример был
разложен вручную. Эта глава даёт стоящую за этим теорию представлений, превращает её в алгоритм,
который библиотека применяет к любому списку матриц, и применяет его к двум симметриям спиновых систем:
$U(1)$ (число частиц или намагниченность) и $SU(2)$ (повороты спинов).

**Вы узнаете** полную приводимость и лемму Шура, то, как разложение на неприводимые определяет
коммутант и порождённую алгебру, как характеры подсчитывают кратности, почему случайный элемент
коммутанта находит разложение численно, что такое зарядовые секторы, неприводимые спиновые
представления и разложение $L$ спинов $\frac12$, а также почему у $SU(2)$-инвариантных гамильтонианов
блоки гораздо меньше, чем у $U(1)$-инвариантных.

Всюду далее $S'$ — коммутант множества $S\subseteq M_n(\mathbb C)$ (глава 2). Если $S^*=S$, то $S''$ —
**алгебра фон Неймана, порождённая** $S$: это $*$-алгебра, содержащая $S$, и всякая $*$-алгебра
$\mathcal B\supseteq S$ удовлетворяет $\mathcal B=\mathcal B''\supseteq S''$ (теорема 2.2).
$\mathcal A'_h=\{K\in\mathcal A':K=K^*\}$ — вещественное векторное пространство. На
$\mathbb C^k\otimes\mathbb C^m$ используется кронекеровский порядок из главы 7 ($e_a\otimes f_j$ — базисный
вектор с номером $am+j$, нумерация с $0$), так что $1_k\otimes x$ — блочно-диагональная матрица с $k$
копиями $x$. $|G|$ — порядок конечной группы, и $\hbar=1$.

## 4.1 Представления и лемма Шура

**Определение.** Пусть $X$ — множество с инволюцией $x\mapsto x^*$ ($x^{**}=x$).
**Представлением** $X$ на $\mathbb C^k$ называется отображение $\rho:X\to M_k(\mathbb C)$ с
$\rho(x^*)=\rho(x)^*$. Подпространство $K$ **инвариантно**, если $\rho(x)K\subseteq K$ для всех $x$;
$\rho$ **неприводимо**, если $k\ge1$ и $0,\mathbb C^k$ — единственные инвариантные подпространства.
**Сплетающий оператор** из $\rho$ в представление $\rho'$ на $\mathbb C^{k'}$ — это $k'\times k$-матрица
$T$ с $T\rho(x)=\rho'(x)T$ для всех $x$; такие операторы образуют векторное пространство
$\operatorname{Hom}(\rho,\rho')$, и $\operatorname{Hom}(\rho,\rho)=\rho(X)'$. Представления $\rho$ и $\rho'$
**эквивалентны**, $\rho\cong\rho'$, если существует обратимый сплетающий оператор.

Встречаются три случая: группа $G$ с $g^*=g^{-1}$, где мы требуем, чтобы $\rho$ было гомоморфизмом
$G\to U(k)$ (**унитарное представление**; тогда $\rho(g^{-1})=\rho(g)^{-1}=\rho(g)^*$); $*$-алгебра со
своим сопряжением, где $\rho$ — унитальный $*$-гомоморфизм; и метки образующих с $x^*=x$, например
$\{x,y,z\}$, представляемые наборами эрмитовых матриц, такими как спиновые матрицы (раздел 4.5). Если $K$
инвариантно и столбцы $Q\in M_{k\times d}$ образуют ортонормированный базис $K$, то $\rho(x)Q=QQ^*\rho(x)Q$,
так что $\rho_K(x)=Q^*\rho(x)Q$ — снова представление (гомоморфизм, если таково $\rho$), **сужение**
$\rho$ на $K$.

**Лемма 4.1 (полная приводимость).** Пусть $\rho$ — представление $X$ на $\mathbb C^n$.
(a) Если $K$ инвариантно, то инвариантно и $K^\perp$. (b) $\mathbb C^n$ есть ортогональная прямая сумма
инвариантных подпространств, на каждом из которых $\rho$ неприводимо.

*Доказательство.* (a) Для $y\in K^\perp$ и $z\in K$ имеем $\langle z,\rho(x)y\rangle=\langle\rho(x^*)z,y\rangle=0$,
поскольку $\rho(x^*)z\in K$. (b) Индукция по $n$. Ненулевое инвариантное подпространство $K_1$ минимальной
размерности неприводимо; $K_1^\perp$ инвариантно по (a), а инвариантные подпространства сужения на
$K_1^\perp$ являются инвариантными подпространствами $\rho$, так что к нему применимо предположение
индукции. $\square$

Для конечной группы унитарность ничего не стоит: если $\rho(g)$ лишь обратимы, то
$P=\sum_g\rho(g)^*\rho(g)>0$ удовлетворяет $\rho(h)^*P\rho(h)=P$, так что $P^{1/2}\rho(h)P^{-1/2}$ унитарен
(трюк Вейля).

**Лемма 4.2 (лемма Шура).** Пусть $\rho,\rho'$ — неприводимые представления $X$ на
$\mathbb C^k,\mathbb C^{k'}$ и $T\in\operatorname{Hom}(\rho,\rho')$.
(a) $T=0$ или $T$ обратим. (b) $\operatorname{Hom}(\rho,\rho)=\mathbb C1$.
(c) $T^*\in\operatorname{Hom}(\rho',\rho)$, и если $T\neq0$, то $T=\lambda U$ с $\lambda>0$ и унитарным $U$.
Следовательно, $\dim\operatorname{Hom}(\rho,\rho')$ равна $1$, если $\rho\cong\rho'$, и $0$ иначе, а
эквивалентные неприводимые представления унитарно эквивалентны.

*Доказательство.* (a) и (b) были доказаны для групп после предложения 2.7, и доказательство использует
лишь то, что $\ker T$ инвариантно относительно $\rho$ ($T\rho(x)v=\rho'(x)Tv=0$, если $Tv=0$), а образ $T$ —
относительно $\rho'$: поэтому $\ker T\in\{0,\mathbb C^k\}$ и $\operatorname{ran}T\in\{0,\mathbb C^{k'}\}$; а для
(b) $T-\lambda1$, где $\lambda$ — собственное значение $T$, — необратимый сплетающий оператор, значит, $0$.
(c) Сопряжение равенства $T\rho(x^*)=\rho'(x^*)T$ даёт $\rho(x)T^*=T^*\rho'(x)$. Тогда
$T^*T\in\operatorname{Hom}(\rho,\rho)$, поэтому $T^*T=\lambda^21$ с $\lambda>0$, если $T\neq0$ ($T^*T\ge0$
ненулевой), и $U=T/\lambda$ квадратна по (a) и удовлетворяет $U^*U=1$. Наконец, если
$T_1,T_2\in\operatorname{Hom}(\rho,\rho')$ ненулевые, то $T_1^{-1}T_2\in\operatorname{Hom}(\rho,\rho)=\mathbb C1$. $\square$

Если все $\rho(x)$ попарно коммутируют (абелева группа), то неприводимое $\rho$ одномерно: каждое
$\rho(x)$ лежит в $\operatorname{Hom}(\rho,\rho)=\mathbb C1$, так что всякое подпространство инвариантно.

## 4.2 Изотипическое разложение и коммутант

**Теорема 4.3 (изотипическое разложение).** Пусть $\rho$ — представление $X$ на $\mathbb C^n$.
Существуют неприводимые попарно неэквивалентные представления $\rho_1,\dots,\rho_C$ множества $X$ на
$\mathbb C^{k_c}$, целые числа $m_c\ge1$ с $\sum_ck_cm_c=n$ и унитарный оператор
$V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n$, такие что
$V^*\rho(x)V=\bigoplus_c\rho_c(x)\otimes1_{m_c}$ для всех $x$. Для всякого такого разложения
$$
V^*\rho(X)'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),\qquad
V^*\rho(X)''V=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},
$$
и $m_c=\dim\operatorname{Hom}(\rho_c,\rho)$; неприводимое представление $\sigma$ множества $X$ имеет
$\operatorname{Hom}(\sigma,\rho)\neq0$ тогда и только тогда, когда $\sigma\cong\rho_c$ ровно для одного $c$.

*Доказательство.* *Существование.* Лемма 4.1 даёт попарно ортогональные инвариантные подпространства
$W_1,\dots,W_M$ с суммой $\mathbb C^n$ и неприводимыми сужениями $\sigma_i(x)=Q_i^*\rho(x)Q_i$. Разобьём $W_i$
на классы эквивалентных сужений и выберем в каждом классе $c$ представителя $\rho_c=\sigma_{i_0}$. Для
каждого другого $i$ из этого класса лемма 4.2(c) даёт унитарный $U_i$ с $\sigma_i(x)=U_i\rho_c(x)U_i^*$;
замена $Q_i$ на $Q_iU_i$ даёт $Q_i^*\rho(x)Q_i=\rho_c(x)$. Занумеруем элементы класса $c$ числами
$j=0,\dots,m_c-1$, обозначим их базисы $Q_{c,j}$ и положим $V(e_a\otimes f_j)=Q_{c,j}e_a$ на $c$-м слагаемом.
Так как $W_i$ инвариантны и ортогональны, $Q_{c,j}^*\rho(x)Q_{d,l}=\delta_{cd}\delta_{jl}\rho_c(x)$, а это и есть
требуемый вид.

*Коммутант.* $\rho(x)Q_{d,l}=Q_{d,l}\rho_d(x)$ и $Q_{c,j}^*\rho(x)=\rho_c(x)Q_{c,j}^*$ (инвариантность и
сопряжённое к ней равенство, применённое к $x^*$). Поскольку $V$ унитарен, $Y\in M_n$ коммутирует со всеми
$\rho(x)$ тогда и только тогда, когда каждый блок $Y_{cj,dl}=Q_{c,j}^*YQ_{d,l}$ удовлетворяет
$Y_{cj,dl}\rho_d(x)=\rho_c(x)Y_{cj,dl}$, т. е. лежит в $\operatorname{Hom}(\rho_d,\rho_c)$, которое равно $0$ при
$c\neq d$ и $\mathbb C1$ при $c=d$ (лемма 4.2). Итак, $Y\in\rho(X)'$ тогда и только тогда, когда
$Y_{cj,dl}=\delta_{cd}(y_c)_{jl}1_{k_c}$, т. е. $V^*YV=\bigoplus_c1_{k_c}\otimes y_c$.

*Бикоммутант.* Если $Z$ коммутирует с $\rho(X)'$, то он коммутирует с проекторами на слагаемые
(возьмём $y_c=1$ и $y_d=0$ при $d\neq c$), так что $V^*ZV=\bigoplus_cZ_c$, и $Z_c$ коммутирует с
$1_{k_c}\otimes M_{m_c}$; по лемме 2.5, шаг 5 (с переставленными тензорными сомножителями),
$Z_c\in M_{k_c}\otimes1_{m_c}$. Обратно, такие операторы коммутируют с $\rho(X)'$.

*Кратности.* Для $T\in\operatorname{Hom}(\sigma,\rho)$ с неприводимым $\sigma$ то же рассуждение показывает, что
$Q_{c,j}^*T\in\operatorname{Hom}(\sigma,\rho_c)$, а это пространство есть $\mathbb Ct_{c}$ с обратимым $t_c$, если
$\sigma\cong\rho_c$, и $0$ иначе; и $T=\sum_{c,j}Q_{c,j}Q_{c,j}^*T$. Так как $\rho_c$ попарно неэквивалентны,
$\sigma$ эквивалентно не более чем одному из них, скажем $\rho_c$ (иначе $T=0$), и тогда $T=\sum_ja_jQ_{c,j}t_c$;
обратно, каждое $Q_{c,j}t_c$ лежит в $\operatorname{Hom}(\sigma,\rho)$, и эти $m_c$ отображений имеют ортогональные
образы, так что $\dim\operatorname{Hom}(\sigma,\rho)=m_c$. $\square$

Для группы $\rho(G)''$ — линейная оболочка $\rho(G)$ (раздел 2.4), и теорема 4.3 — это предложение 2.7,
прочитанное в обратную сторону: там структурная теорема порождала неприводимые представления, здесь они
порождают структуру; в обоих случаях $\dim\rho(X)'=\sum_cm_c^2$ и $\dim\rho(X)''=\sum_ck_c^2$. Образ
$P_c=\sum_jQ_{c,j}Q_{c,j}^*$, минимального центрального проектора алгебры $\rho(X)''$, — это **изотипическая
компонента** представления $\rho_c$.

> **Физический смысл.** Если $\rho(X)$ — симметрия, то всякий симметричный гамильтониан $H\in\rho(X)'$
> имеет вид $V(\bigoplus_c1_{k_c}\otimes h_c)V^*$: каждое собственное значение блока $h_c$ размера
> $m_c\times m_c$ встречается $k_c$ раз. Уровни с $k_c$-кратным вырождением — это **мультиплеты**
> симметрии, а блоки, которые нужно диагонализовать, имеют размеры $m_c$, равные кратностям, а не размерам
> мультиплетов.

Для конечной группы кратности можно вычислить, ничего не раскладывая.

**Теорема 4.4 (ортогональность Шура).** Пусть $\rho,\rho'$ — неприводимые унитарные представления
конечной группы $G$ на $\mathbb C^k,\mathbb C^{k'}$. Тогда
$$
\frac1{|G|}\sum_{g\in G}\rho(g)_{ij}\,\overline{\rho'(g)_{i'j'}}=
\begin{cases}\delta_{ii'}\delta_{jj'}/k&\text{if }\rho=\rho',\\0&\text{if }\rho\not\cong\rho'.\end{cases}
$$

*Доказательство.* Для $k\times k'$-матрицы $Y$ положим $T=\frac1{|G|}\sum_g\rho(g)Y\rho'(g)^*$. Для $h\in G$
замена $g\to hg$ и равенство $\rho'(hg)^*=\rho'(g)^*\rho'(h)^{-1}$ дают $\rho(h)T=T\rho'(h)$, так что
$T\in\operatorname{Hom}(\rho',\rho)$. Если $\rho\not\cong\rho'$, то $T=0$; если $\rho=\rho'$, то $T=\lambda1$ с
$k\lambda=\operatorname{Tr}T=\operatorname{Tr}Y$ (лемма 4.2). Для $Y=|j\rangle\langle j'|$ элемент $T_{ii'}$ есть
левая часть, а $\operatorname{Tr}Y=\delta_{jj'}$. $\square$

**Характер** представления — это $\chi_\rho(g)=\operatorname{Tr}\rho(g)$; эквивалентные представления имеют
равные характеры, поскольку $\operatorname{Tr}(T\rho(g)T^{-1})=\operatorname{Tr}\rho(g)$. Положим
$\langle\varphi,\psi\rangle_G=\frac1{|G|}\sum_g\overline{\varphi(g)}\psi(g)$.

**Следствие 4.5 (характеры).** Пусть $G$ конечна.
(a) Для неприводимых унитарных $\rho,\rho'$: $\langle\chi_\rho,\chi_{\rho'}\rangle_G$ равно $1$, если $\rho\cong\rho'$, и $0$ иначе.
(b) Для унитарного представления $\pi$, разложенного как в теореме 4.3, и неприводимого $\rho$
кратность $\rho$ в $\pi$ равна
$$
m_\rho=\frac1{|G|}\sum_{g\in G}\overline{\chi_\rho(g)}\,\chi_\pi(g)
$$
($m_\rho=m_c$, если $\rho\cong\rho_c$, и $m_\rho=0$, если $\rho$ не эквивалентно ни одному $\rho_c$).
(c) $\dim\pi(G)'=\sum_cm_c^2=\frac1{|G|}\sum_{g\in G}|\chi_\pi(g)|^2$. В частности, $\pi$ неприводимо тогда и
только тогда, когда это число равно $1$.

*Доказательство.* (a) Заменим $\rho'$ эквивалентным представлением (с тем же характером). Просуммируем
теорему 4.4 по $i=j$, $i'=j'$ и перейдём к комплексному сопряжению: получится $\sum_{i,i'}\delta_{ii'}/k=1$,
если $\rho=\rho'$, и $0$, если $\rho\not\cong\rho'$. (b) След унитарно инвариантен и
$\operatorname{Tr}(a\otimes1_m)=m\operatorname{Tr}a$, поэтому $\chi_\pi=\sum_cm_c\chi_{\rho_c}$; применим (a). (c) По тому же
разложению и (a) $\langle\chi_\pi,\chi_\pi\rangle_G=\sum_cm_c^2$, а это $\dim\pi(G)'$ по теореме 4.3. $\square$

Второе доказательство (c) показывает, как библиотека делает выборки из коммутантов. **Твирлинг**
$P(Y)=\frac1{|G|}\sum_g\pi(g)Y\pi(g)^*$ идемпотентен с образом $\pi(G)'$ (глава 2, упражнение 3) и
самосопряжён относительно $\langle Y,Z\rangle_{HS}=\operatorname{Tr}(Y^*Z)$ (замена $g\to g^{-1}$), а значит, является
ортогональным проектором на $\pi(G)'$, и $\dim\pi(G)'$ — его след. Отображение $Y\mapsto AYB$ имеет след
$\operatorname{Tr}A\operatorname{Tr}B$ (его диагональные элементы равны $A_{ii}B_{jj}$), что даёт $\frac1{|G|}\sum_g|\chi_\pi(g)|^2$.

Повторим пример из раздела 2.4 с помощью характеров. След матрицы перестановки равен числу базисных
векторов, которые она оставляет на месте, а $\pi(\sigma)$ оставляет на месте $|b_1b_2b_3\rangle$ тогда и только
тогда, когда $b$ постоянно на циклах $\sigma$, поэтому $\chi_\pi(\sigma)=2^{\#\text{cycles}(\sigma)}$. У группы $S_3$ есть
тривиальное, знаковое и стандартное представления (перестановочное представление $P_\sigma$ на $\mathbb C^3$ за
вычетом тривиального) с характерами $1$, $\operatorname{sgn}\sigma=(-1)^{3-\#\text{cycles}}$ и $\operatorname{Tr}P_\sigma-1$;
они вещественны. Код находит $P_\sigma$ как сужение $\pi(\sigma)$ на линейную оболочку
$|100\rangle,|010\rangle,|001\rangle$, а $\#\text{cycles}(\sigma)=\dim\ker(P_\sigma-1)$ (векторы, постоянные на циклах).

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

Тривиальное представление встречается $4$ раза, знаковое — ни разу, стандартное — дважды; поскольку
$4\cdot1+2\cdot2=8$, больше ничего не встречается, и $\dim\pi(S_3)'=\frac16(64+3\cdot16+2\cdot4)=20$. Это
разложение из раздела 2.4, полученное теперь из таблицы по шесть чисел на представление.

## 4.3 От образующих к алгебре: рандомизированный алгоритм

Теорема 4.3 описывает разложение, но не вычисляет его: её доказательству нужны неприводимые
инвариантные подпространства. Мурота, Канно, Кодзима и Кодзима (*Japan J. Indust. Appl. Math.* 27 (2010)
125–160; см. также Maehara and Murota, *ibid.* 263–293) заметили, что их находят собственные подпространства
*случайного элемента коммутанта*. Пусть $\mathcal A$ — алгебра фон Неймана, порождённая данными матрицами
$S_1,\dots,S_J$ и их сопряжёнными, с разложением $V^*\mathcal AV=\bigoplus_cM_{k_c}\otimes1_{m_c}$,
$V^*\mathcal A'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$ (теорема 2.6); $V$ используется только в анализе.

**Лемма 4.6.** Пусть $K\in\mathcal A'_h$ и $V^*KV=\bigoplus_c1_{k_c}\otimes\kappa_c$. Каждое собственное
подпространство $K$ $\mathcal A$-инвариантно. Если все $N=\sum_cm_c$ собственных значений всех $\kappa_c$ вместе
попарно различны, то собственные подпространства $K$ — это $N$ подпространств $V(\mathbb C^{k_c}\otimes e)$
(в $c$-м слагаемом, $e$ — собственный вектор $\kappa_c$), и $\mathcal A$ действует на каждом из них
неприводимо, как $M_{k_c}$.

*Доказательство.* Если $Kv=\lambda v$ и $a\in\mathcal A$, то $Kav=aKv=\lambda av$. Запишем
$\kappa_c=\sum_j\mu_{cj}|e_{cj}\rangle\langle e_{cj}|$ с ортонормированными $e_{cj}$ (теорема 1.3). Тогда
$V^*KV=\sum_{c,j}\mu_{cj}R_{cj}$, где $R_{cj}$ — проектор на $\mathbb C^{k_c}\otimes e_{cj}$; проекторы
$R_{cj}$ ортогональны с суммой $1$, а $\mu_{cj}$ различны, так что это спектральные проекторы. На
$\mathbb C^{k_c}\otimes e_{cj}$ элемент $V^*aV=\bigoplus_da_d\otimes1$ действует как $a_c$, а $a_c$ пробегают
всю $M_{k_c}$, у которой нет инвариантных подпространств, кроме $0$ и $\mathbb C^{k_c}$. $\square$

**Предложение 4.7 (общее положение).** Пусть $K$ и $Y$ — независимые случайные элементы $\mathcal A'_h$ и
$\mathcal A'$, законы распределения которых имеют плотности относительно меры Лебега на этих вещественных
векторных пространствах. С вероятностью единица: (a) выполнено условие леммы 4.6; (b) $\Pi_iY\Pi_j\neq0$
для любых двух спектральных проекторов $\Pi_i,\Pi_j$ оператора $K$, образы которых лежат в одном и том же
слагаемом.

*Доказательство.* Мы используем следующее: *ненулевой вещественный многочлен $p\not\equiv0$ на $\mathbb R^d$
обращается в нуль лишь на множестве меры Лебега нуль.* При $d=1$ у него конечное число корней. При $d>1$
запишем $p=\sum_jp_j(x')t^j$ с $x'\in\mathbb R^{d-1}$ и некоторым $p_j\not\equiv0$; вне множества меры нуль
$\{p_j=0\}$ (индукция) сечение $\{t:p(x',t)=0\}$ конечно, и теорема Тонелли даёт меру $0$.

(a) Эрмитова матрица $\kappa=\bigoplus_c\kappa_c\in M_N$ вещественно-линейно зависит от $K$. Определитель Грама
$\delta(K)=\det[\operatorname{Tr}\kappa^{i+j}]_{i,j=0}^{N-1}$ системы $1,\kappa,\dots,\kappa^{N-1}$ относительно
$\langle\cdot,\cdot\rangle_{HS}$ — вещественный многочлен от координат $K$, отличный от нуля тогда и только тогда,
когда эти степени линейно независимы. Так как $p(\kappa)=U\operatorname{diag}(p(\lambda_l))U^*$ для
$\kappa=U\operatorname{diag}(\lambda)U^*$, ненулевой многочлен степени $<N$ аннулирует $\kappa$ тогда и только тогда,
когда он обращается в нуль на $\operatorname{spec}\kappa$, что возможно тогда и только тогда, когда у $\kappa$
меньше $N$ различных собственных значений (возьмём $\prod_\mu(t-\mu)$ по ним). Итак, $\delta(K)\neq0$ — это
условие леммы 4.6, и $\delta\not\equiv0$: он отличен от нуля в точке
$V(\bigoplus_c1\otimes\operatorname{diag}(\text{различные вещественные}))V^*\in\mathcal A'_h$. Следовательно, $\{\delta=0\}$ — множество меры нуль.

(b) На открытом множестве, где выполнено (a), спектральные проекторы — многочлены от $K$, коэффициенты
которых — непрерывные функции собственных значений (доказательство теоремы 1.5(a)). Зафиксируем такое $K$.
По лемме 4.6 $V^*\Pi_iV$ и $V^*\Pi_jV$ проектируют на $\mathbb C^{k_c}\otimes e_i$ и $\mathbb C^{k_c}\otimes e_j$, и
для $V^*YV=\bigoplus_d1\otimes y_d$ имеем $\Pi_iY\Pi_j=0$ тогда и только тогда, когда $\langle e_i,y_ce_j\rangle=0$.
Отображение $Y\mapsto\langle e_i,y_ce_j\rangle$ вещественно-линейно $\mathcal A'\to\mathbb C$ и ненулевое (его
значение на $V(1\otimes|e_i\rangle\langle e_j|)V^*$ равно $1$), так что его ядро — множество меры нуль. Так как
$Y$ не зависит от $K$, теорема Тонелли даёт вероятность $0$. $\square$

**Лемма 4.8 (выравнивание копий).** Пусть $W_i,W_j$ — инвариантные подпространства, на которых $\mathcal A$
действует неприводимо, с ортонормированными базисами $Q_i,Q_j$ и сужениями $\sigma_i(a)=Q_i^*aQ_i$, и пусть
$Y\in\mathcal A'$. Тогда $T=Q_i^*YQ_j\in\operatorname{Hom}(\sigma_j,\sigma_i)$. Если $T\neq0$, то $T=\lambda U$ с
$\lambda>0$ и унитарным $U$, и $Q_j'=Q_jU^*$ удовлетворяет $Q_j'^*aQ_j'=\sigma_i(a)$ для всех $a\in\mathcal A$.
Если $W_i,W_j$ лежат в разных слагаемых ($W_i\subseteq P_c\mathbb C^n$, $W_j\subseteq P_d\mathbb C^n$, $c\neq d$),
то $T=0$.

*Доказательство.* Как в доказательстве теоремы 4.3, $aQ_j=Q_j\sigma_j(a)$ и $Q_i^*a=\sigma_i(a)Q_i^*$, поэтому
$T\sigma_j(a)=Q_i^*YaQ_j=Q_i^*aYQ_j=\sigma_i(a)T$. Лемма 4.2 даёт $T=\lambda U$ и
$U\sigma_j(a)U^*=\sigma_i(a)$, а это $Q_j'^*aQ_j'$. Наконец, $P_c\in\mathcal A\cap\mathcal A'$, $Q_i=P_cQ_i$ и
$P_cQ_j=0$, так что $T=Q_i^*P_cYQ_j=Q_i^*YP_cQ_j=0$. $\square$

Теперь алгоритм (`decompose.wedderburn` с интерфейсами `from_group` и `from_generators`) таков:

1. *Выборка из коммутанта.* Для всех элементов конечной группы $Y=P(Z)$ и $K=P(Z+Z^*)$ с твирлингом
   $P$ и гауссовскими матрицами $Z$: это образы гауссовских законов при линейных сюръекциях на $\mathcal A'$ и
   $\mathcal A'_h$, а значит, они имеют плотности. Для образующих $Y,K$ — гауссовские линейные комбинации
   базиса $\mathcal A'=\{X:[S_j,X]=[S_j^*,X]=0\}$, ядра линейного отображения на $M_n$.
2. *Расщепление.* Собственные подпространства $K$ неприводимы (лемма 4.6, предложение 4.7(a)); их
   размерности — это $k_c$.
3. *Группировка и выравнивание.* Блоки $Q_i^*YQ_j$ ненулевые ровно внутри одного слагаемого
   (предложение 4.7(b), лемма 4.8); классы дают секторы, их размеры — $m_c$, а полярные множители
   выравнивают копии.
4. *Сборка* $V$ в кронекеровском порядке, так что $V^*aV=\bigoplus_ca_c\otimes1_{m_c}$.

В арифметике с плавающей точкой близкие собственные значения разных копий держатся вместе и
расщепляются сжатием нового случайного элемента коммутанта (docstring модуля `decompose`). В любом
случае результат в конце сертифицируется:

**Предложение 4.9 (сертификат).** Пусть $\mathcal A$ порождена $S_1,\dots,S_J$ и их сопряжёнными, и
$\mathcal B=V(\bigoplus_cM_{k_c}\otimes1_{m_c})V^*$ с унитарным $V$. Если все $S_j\in\mathcal B$ и
$\dim\mathcal A'=\sum_cm_c^2$, то $\mathcal A=\mathcal B$.

*Доказательство.* $\mathcal B$ — $*$-алгебра, содержащая все $S_j,S_j^*$, поэтому $\mathcal A\subseteq\mathcal B$ и
$\mathcal B'\subseteq\mathcal A'$. Теорема 4.3, применённая к представлению $\mathcal B$ на $\mathbb C^n$ (блочные
отображения $b\mapsto b_c$ неприводимы и неэквивалентны, поскольку $c$-й центральный проектор переводится
одним из них в $1$, а остальными — в $0$), даёт $\dim\mathcal B'=\sum_cm_c^2=\dim\mathcal A'$. Значит,
$\mathcal B'=\mathcal A'$ и $\mathcal B=\mathcal B''=\mathcal A''=\mathcal A$. $\square$

`Decomposition.verify()` сообщает относительные расстояния образующих (и новых элементов коммутанта) до
$\mathcal B$ (соответственно $\mathcal B'$) и сравнивает $\sum_cm_c^2$ с $\dim\mathcal A'$ (из характеров для группы,
следствие 4.5(c); из ядра для образующих). Запустим его для $S_3$ и трёх кубитов, а затем заглянем внутрь:
собственные подпространства одного случайного $K$ и сплетающие операторы одного случайного $Y$.

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
s = torch.linalg.svdvals(Wa.mH @ Y @ Wb)                               # Lemma 4.8: a multiple of a unitary
assert s[1] > 1e-3 and abs(s[0] - s[1]) < 1e-10
assert all(close(W.mH @ Y @ Wa, 0) for W in spaces if W.shape[1] == 1)   # different summands: zero
```

Основная стоимость — коммутант: $O(|G|n^3)$ на одну выборку твирлинга; для образующих библиотека с
помощью случайного элемента $\mathcal A$ сводит задачу о ядре к $N=\sum_ck_cm_c^2$ неизвестным, но всё равно
требует памяти $O(n^2N)$ — около $10$ ГБ для полного спина восьми спинов $\frac12$ ($n=256$, $N=4900$).
Поэтому для спиновых цепочек в разделе 4.6 используется модуль `su2`, которому ответ известен заранее.

## 4.4 $U(1)$ и зарядовые секторы

**Предложение 4.10.** Пусть $Q$ — эрмитова матрица на $\mathbb C^n$ с целыми собственными значениями $q$
(**зарядами**), спектральными проекторами $P_q$ и $d_q=\operatorname{rank}P_q$, и $U_\theta=e^{i\theta Q}$; мы пишем
$U(1)'=\{U_\theta:\theta\in\mathbb R\}'$ и $U(1)''$ для его коммутанта.

1. $\theta\mapsto U_\theta$ — унитарное представление $U(1)=\mathbb R/2\pi\mathbb Z$, и всякое непрерывное
   унитарное представление $U(1)$ на $\mathbb C^n$ имеет такой вид.
2. Для $X\in M_n$ следующие условия эквивалентны: (i) $XU_\theta=U_\theta X$ для всех $\theta$; (ii) $[X,Q]=0$;
   (iii) $X=\sum_qP_qXP_q$. Следовательно, $U(1)'\cong\bigoplus_qM_{d_q}(\mathbb C)$ (блочно-диагональные
   операторы) и $U(1)''=\operatorname{span}\{P_q\}$: все $k_c=1$, а кратности равны $d_q$.
3. (*Сохраняющийся заряд.*) Для эрмитова $H$ равенство $e^{itH}Qe^{-itH}=Q$ выполнено при всех $t$ тогда и
   только тогда, когда $[H,Q]=0$, т. е. когда $H\in U(1)'$. Тогда каждый вес сектора
   $\operatorname{Tr}(P_q\rho(t))$ для $\rho(t)=e^{-itH}\rho e^{itH}$ постоянен.

*Доказательство.* (1) $U_\theta=\sum_qe^{iq\theta}P_q$ унитарен, $U_\theta U_{\theta'}=U_{\theta+\theta'}$ (теорема 1.5(c)) и
$U_{2\pi}=1$. Обратно, для непрерывного представления $(U_\theta)_{\theta\ge0}$ — непрерывная полугруппа
матриц, поэтому $U_\theta=e^{\theta B}$ при $\theta\ge0$ по лемме 8.5 (доказана в главе 8; доказательство
использует лишь конечномерность), а при $\theta<0$ имеем $U_\theta=U_{-\theta}^{-1}=e^{\theta B}$. Дифференцирование
$U_\theta^*U_\theta=1$ в $0$ даёт $B^*+B=0$, так что $B=iQ$ с эрмитовым $Q$, а $U_{2\pi}=1$ влечёт
$e^{2\pi iq}=1$, $q\in\mathbb Z$, для каждого собственного значения. (2) (i)$\Rightarrow$(ii): продифференцируем
в $\theta=0$. (ii)$\Rightarrow$(iii): $X$ коммутирует с $P_q=f(Q)$ (теорема 1.5(d)), поэтому
$P_qXP_{q'}=XP_qP_{q'}=0$ при $q\neq q'$. (iii)$\Rightarrow$(i): $XU_\theta=\sum_qe^{iq\theta}P_qXP_q=U_\theta X$.
Блочно-диагональные операторы образуют $\bigoplus_qM_{d_q}\otimes1_1$, коммутант которой есть
$\bigoplus_q1_{d_q}\otimes M_1=\operatorname{span}\{P_q\}$ (теорема 2.6).
(3) $\frac d{dt}e^{itH}Qe^{-itH}=ie^{itH}[H,Q]e^{-itH}$ (предложение 8.2) обращается в нуль при всех $t$ тогда
и только тогда, когда $[H,Q]=0$. Тогда $e^{itH}$ коммутирует с $P_q$ (теорема 1.5(d), дважды), и
$\operatorname{Tr}(P_q\rho(t))=\operatorname{Tr}(e^{itH}P_qe^{-itH}\rho)=\operatorname{Tr}(P_q\rho)$. $\square$

Таким образом, неприводимые представления $U(1)$ — это характеры $\theta\mapsto e^{iq\theta}$, $q\in\mathbb Z$
(одномерные, как у всякой абелевой группы), а кратность заряда $q$ равна $d_q$. При тензорном умножении
заряды складываются: $e^{i\theta(Q_1\otimes1+1\otimes Q_2)}=e^{i\theta Q_1}\otimes e^{i\theta Q_2}$ действует на образе
$P_{q_1}\otimes P_{q_2}$ с зарядом $q_1+q_2$. Это правило слияния $(q_1,q_2)\mapsto q_1+q_2$ для тензорных
произведений со слиянием из главы 7.

Для $n$ кубитов зарядом служит **оператор числа частиц** $\hat N=\sum_jn_j$, $n_j=|1\rangle\langle1|$ на кубите
$j$: $\hat N|b\rangle=|b|\,|b\rangle$ с весом Хэмминга $|b|$. Его секторы имеют размерность $\binom nw$, так что
операторы, сохраняющие вес Хэмминга, образуют
$$
U(1)'\cong\bigoplus_{w=0}^nM_{\binom nw}(\mathbb C),
$$
а это `chains.SpinChain(n).algebra` — произведение со слиянием $n$ копий $M_1\oplus M_1$ с зарядами $0,1$
(бит $1$ = спин вверх, $\hat N=S^z+\frac n2$; в `decompose.total_spin` спин вверх — это, напротив, *первый*
базисный вектор $\mathbb C^2$, что переименовывает базис, но не меняет секторы). Всякий XXZ-гамильтониан
$\sum_{\langle ij\rangle}J(s^x_is^x_j+s^y_is^y_j+\Delta s^z_is^z_j)$ сохраняет $\hat N$, потому что член перескока
$s^x_is^x_j+s^y_is^y_j=\frac12(s^+_is^-_j+s^-_is^+_j)$ переносит спин вверх с одного узла на другой.

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

> **Физический смысл.** Симметрия $U(1)$ — это сохранение числа частиц, заряда или намагниченности.
> Предложение 4.10(3) — теорема Нётер в этой ситуации: заряд сохраняется тогда и только тогда, когда
> гамильтониан симметричен, и тогда динамика никогда не покидает сектор (теорема 8.3). Если вдобавок
> ничто не может измерить когерентности между секторами, то секторы являются суперотборными секторами (глава 3).

## 4.5 $SU(2)$: спин

**Определение.** **Спиновое представление** на $\mathbb C^d$ — это тройка эрмитовых матриц
$S_x,S_y,S_z$ с
$$
[S_x,S_y]=iS_z,\qquad[S_y,S_z]=iS_x,\qquad[S_z,S_x]=iS_y,
$$
т. е. представление меток $\{x,y,z\}$ (раздел 4.1) с этими соотношениями. Основной пример —
$s_a=\frac12\sigma_a$ на $\mathbb C^2$ (матрицы Паули). Если $(S_a)$ и $(T_a)$ — спиновые представления, то
таково и $S_a\otimes1+1\otimes T_a$ (перекрёстные члены коммутируют); поэтому **полный спин**
$S_a=\sum_is_a^{(i)}$ системы $L$ спинов $\frac12$ (`decompose.total_spin(L)`) — спиновое представление. Повороты
$e^{-i\theta\,\mathbf n\cdot\mathbf S}$ коммутируют с $X$ при всех $\theta,\mathbf n$ тогда и только тогда, когда $X$
коммутирует с тремя $S_a$ (дифференцирование в $\theta=0$; обратно — теорема 1.5(d)), так что коммутант
«действия $SU(2)$» есть $\{S_x,S_y,S_z\}'$, и по теореме 2.2 повороты и $S_a$ порождают одну и ту же
алгебру. Для $L$ спинов $\frac12$ повороты — это $V^{\otimes L}$, $V\in SU(2)$: всякий такой $V$ имеет вид
$e^{-i\theta\,\mathbf n\cdot\mathbf s}$, и $e^{-i\theta\,\mathbf n\cdot\mathbf S}=(e^{-i\theta\,\mathbf n\cdot\mathbf s})^{\otimes L}$.

Положим $S_\pm=S_x\pm iS_y$, так что $S_-=S_+^*$, и $\mathbf S^2=S_x^2+S_y^2+S_z^2$. Из соотношений следует
$$
[S_z,S_\pm]=\pm S_\pm,\qquad[S_+,S_-]=2S_z,\qquad\mathbf S^2=S_-S_++S_z^2+S_z,
$$
и $\mathbf S^2$ коммутирует с каждым $S_a$ (например,
$[\mathbf S^2,S_z]=S_x[S_x,S_z]+[S_x,S_z]S_x+S_y[S_y,S_z]+[S_y,S_z]S_y=-i\{S_x,S_y\}+i\{S_y,S_x\}=0$).

**Теорема 4.11 (неприводимые спиновые представления).** Для всякого $s\in\{0,\frac12,1,\frac32,\dots\}$
матрицы `decompose.spin_matrices(s)` — на $\mathbb C^{2s+1}$ с базисом $|m\rangle$, $m=s,s-1,\dots,-s$,
$$
S_z|m\rangle=m|m\rangle,\qquad S_+|m\rangle=\sqrt{s(s+1)-m(m+1)}\;|m+1\rangle,\qquad S_-=S_+^*,
$$
образуют неприводимое спиновое представление (**спин $s$**) с $\mathbf S^2=s(s+1)1$. Всякое неприводимое
спиновое представление унитарно эквивалентно ровно одному из них.

*Доказательство.* *Спин $s$ — спиновое представление.* $[S_z,S_+]=S_+$ очевидно.
$S_+S_-|m\rangle=(s(s+1)-(m-1)m)|m\rangle$ и $S_-S_+|m\rangle=(s(s+1)-m(m+1))|m\rangle$, поэтому $[S_+,S_-]=2S_z$;
эти два соотношения вместе с сопряжёнными к ним эквивалентны трём определяющим.
$\mathbf S^2|m\rangle=(s(s+1)-m(m+1)+m^2+m)|m\rangle=s(s+1)|m\rangle$. *Неприводимость:* инвариантное $W\neq0$
инвариантно относительно $S_z$, собственные значения которого простые, поэтому содержит некоторый
$|m\rangle$; коэффициенты $S_\pm$ обращаются в нуль лишь при $m=\pm s$, так что применение $S_\pm$ даёт все $|m'\rangle$.

*Классификация.* Пусть $(S_a)$ неприводимо на $\mathbb C^d$, $\lambda$ — наибольшее собственное значение $S_z$ и
$v_0$ — единичный собственный вектор. $S_zS_+v_0=S_+(S_z+1)v_0=(\lambda+1)S_+v_0$, поэтому $S_+v_0=0$ в силу
максимальности. Положим $v_j=S_-^jv_0$; тогда $S_zv_j=(\lambda-j)v_j$, и по индукции, используя
$S_+S_-=S_-S_++2S_z$,
$$
S_+v_j=j(2\lambda-j+1)\,v_{j-1}
$$
($S_+v_j=S_-S_+v_{j-1}+2S_zv_{j-1}=[(j-1)(2\lambda-j+2)+2(\lambda-j+1)]v_{j-1}$). Ненулевые $v_j$ — собственные векторы
$S_z$ с различными собственными значениями, поэтому существует последний ненулевой, $v_r$. Линейная оболочка
$v_0,\dots,v_r$ инвариантна относительно $S_z,S_\pm$, а значит, совпадает с $\mathbb C^d$, и $d=r+1$. Из
$0=S_+v_{r+1}=(r+1)(2\lambda-r)v_r$ получаем $\lambda=r/2=:s$. Наконец,
$\|v_j\|^2=\langle v_{j-1},S_+v_j\rangle=j(2s-j+1)\|v_{j-1}\|^2$, так что единичные векторы
$|s-j\rangle=v_j/\|v_j\|$ удовлетворяют $S_+|m\rangle=\sqrt{(s-m)(s+m+1)}\,|m+1\rangle$ с $m=s-j$, а это формула
выше. Разным $s$ отвечают разные размерности. $\square$

**Следствие 4.12.** Всякое спиновое представление на $\mathbb C^n$ унитарно эквивалентно
$\bigoplus_s(\text{spin }s)\otimes1_{m_s}$, и при $d_M=\dim\ker(S_z-M)$
$$
m_s=d_s-d_{s+1}.
$$

*Доказательство.* Теорема 4.3 и теорема 4.11. Каждая копия спина $s'$ содержит ровно один собственный
вектор $S_z$ с собственным значением $M$, если $|M|\le s'$ и $s'-M\in\mathbb Z$, и ни одного иначе; поэтому
$d_M=\sum_{s'\ge|M|,\,s'-M\in\mathbb Z}m_{s'}$, и разность двух соседних значений равна $m_s$. $\square$

**Предложение 4.13 ($L$ спинов $\frac12$).** Полный спин на $(\mathbb C^2)^{\otimes L}$ содержит спин $S$ тогда
и только тогда, когда $S\in\{\frac L2,\frac L2-1,\dots\}$, $S\ge0$, с кратностью
$$
m_S=\binom L{\frac L2-S}-\binom L{\frac L2-S-1}\qquad\Big(\tbinom L{-1}=0\Big),
$$
так что алгебра, порождённая $S_x,S_y,S_z$, есть $\bigoplus_SM_{2S+1}\otimes1_{m_S}$, а её коммутант —
$\bigoplus_S1_{2S+1}\otimes M_{m_S}$.

*Доказательство.* $S_z|b\rangle=(|b|-\frac L2)|b\rangle$, поэтому $d_M=\binom L{L/2+M}=\binom L{L/2-M}$, если
$\frac L2+M\in\{0,\dots,L\}$, и $0$ иначе. Применим следствие 4.12 и теорему 4.3. Значения положительны, поскольку
$\binom L{j+1}/\binom Lj=\frac{L-j}{j+1}>1$ при $j<\frac{L-1}2$. $\square$

При $L=2$: $\frac12\otimes\frac12=1\oplus0$ (**Клебш–Гордан**): триплет $S_-^j|{\uparrow\uparrow}\rangle$, $j=0,1,2$,
симметричный относительно перестановки $F$ двух спинов, и синглет
$\frac1{\sqrt2}(|{\uparrow\downarrow}\rangle-|{\downarrow\uparrow}\rangle)$, антисимметричный и аннулируемый $S_+$ и $S_z$.
При $L=4$: $m_2=1$, $m_1=4-1=3$, $m_0=6-4=2$, так что
$$
\{S_x,S_y,S_z\}''\cong M_5\oplus M_3\oplus M_1\ \text{ with multiplicities }1,3,2,\qquad
\{S_x,S_y,S_z\}'\cong M_1\oplus M_3\oplus M_2,
$$
размерностей $35$ и $14$ — это числа, найденные в главе 10. Библиотека находит то же самое по одним лишь
трём образующим (секторы упорядочены по $(k_c,m_c)$):

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

Какие операторы коммутируют со всеми поворотами? Перестановки узлов: сопряжение посредством $\pi(\sigma)$
переводит $s^{(i)}_a$ в $s^{(\sigma(i))}_a$, так что оставляет $\sum_is^{(i)}_a$ на месте. Обратное утверждение —
классическая теорема.

**Теорема 4.14 (двойственность Шура–Вейля).** Пусть $\pi$ — представление $S_L$, переставляющее сомножители
$(\mathbb C^2)^{\otimes L}$. (a) $\operatorname{span}\pi(S_L)\subseteq\{S_x,S_y,S_z\}'$. (b) Имеет место равенство;
эквивалентно (теорема 2.2), алгебра, порождённая полным спином, есть $\pi(S_L)'$.

Часть (a) только что доказана. Часть (b) приводится без доказательства для общего $L$; она верна для
$(\mathbb C^d)^{\otimes L}$ и диагонального действия $U(d)$ (W. Fulton and J. Harris, *Representation
Theory*, Springer 1991, Lecture 6; R. Goodman and N. Wallach, *Symmetry, Representations, and Invariants*,
Springer 2009, Section 9.1). Для $L\le5$ мы её проверяем: в силу (a) достаточно, чтобы оба пространства
имели одинаковую размерность. При $L=2$ это лемма 5.8 с $k=2$, поскольку
$\{S_a\}'=\{V\otimes V:V\in SU(2)\}'=\{V\otimes V:V\in U(2)\}'$ (фаза действует как скаляр).

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

Напечатанные размерности — $2,5,14,42$, числа Каталана $\binom{2L}L/(L+1)$. Последняя строка — двойственность
в действии: у $\pi(S_4)''$ блоки $k=(1,2,3)$ с кратностями $(5,1,3)$ — пары $(k_c,m_c)$ полного спина с
переставленными ролями. По теореме 4.14(b) и предложению 2.7 пространство кратности спина $S$ несёт
неприводимое представление $S_L$, здесь размерностей $1,3,2$ для $S=2,1,0$.

> **Физический смысл.** Мультиплет спина $S$ состоит из $2S+1$ состояний равной энергии для любого
> вращательно-инвариантного гамильтониана, различаемых значением $S^z=M$. Два спина $\frac12$ образуют
> триплет и синглет; четыре спина — один квинтет, три триплета и два синглета. Пространства кратности
> двух синглетов и трёх триплетов — это свободный от декогеренции кубит и бесшумный кутрит из главы 10.

## 4.6 Спиновые системы

**Модель Гейзенберга** на $L$ спинах $\frac12$ — это $H=\sum_{i<j}J_{ij}\,\mathbf s^{(i)}\cdot\mathbf s^{(j)}$ с
вещественными константами связи $J_{ij}$.

**Предложение 4.15.** (a) (*Тождество Дирака*) $\mathbf s^{(i)}\cdot\mathbf s^{(j)}=\frac12F_{ij}-\frac14$, где $F_{ij}$
переставляет спины $i$ и $j$. Следовательно, всякий гейзенберговский гамильтониан лежит в
$\operatorname{span}\pi(S_L)\subseteq\{S_a\}'$, и $V^*HV=\bigoplus_S1_{2S+1}\otimes H_S$ с блоками $H_S\in M_{m_S}$ в
разложении из предложения 4.13: каждое собственное значение $H_S$ — это $(2S+1)$-кратно вырожденный уровень $H$.
(b) Пространство **векторов старшего веса** $\ker S_+\cap\ker(S_z-S)$ имеет размерность $m_S$, инвариантно
относительно $H$, и $H$ действует на нём как $H_S$ (с точностью до унитарной замены базиса).

*Доказательство.* (a) На двух спинах $\mathbf s^{(1)}\cdot\mathbf s^{(2)}=\frac12(\mathbf S^2-\mathbf s^{(1)2}-\mathbf s^{(2)2})=\frac12\mathbf S^2-\frac34$
(теорема 4.11 с $s=\frac12$). На триплете $\mathbf S^2=2$ и $F=1$ (он натянут на
$S_-^j|{\uparrow\uparrow}\rangle$, а $F$ коммутирует с $S_-$ и оставляет $|{\uparrow\uparrow}\rangle$ на месте); на синглете
$\mathbf S^2=0$ и $F=-1$. Обе части совпадают на этих четырёх векторах. Для узлов $i,j$ из $L$ берём тензорное
произведение с единицей. Остальное — теорема 4.14(a) и теорема 4.3. (b) В разложении из следствия 4.12
вектор с $S_z=S$ имеет вид $\sum_{S'\ge S}|S',S\rangle\otimes w_{S'}$, и $S_+|S',S\rangle\neq0$ при $S'>S$ по теореме
4.11; поэтому $\ker S_+\cap\ker(S_z-S)=|S,S\rangle\otimes\mathbb C^{m_S}$, и на нём $1\otimes H_S$ действует как $H_S$. $\square$

Часть (b) — это то, как модуль `su2` (спиновые цепочки, глава 9) строит $H_S$ без численного разложения:
ортонормированный базис векторов старшего веса строится сложением моментов по Клебшу–Гордану, по одному
узлу за раз, внутри $U(1)$-сектора $\hat N=\frac L2+S$ класса `SpinChain`. Проверим (a) на кольце из четырёх
спинов с помощью `D4`, а затем сравним размеры блоков для восьми спинов.

```python
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

from torch_vn_algebra import su2
u1, sym = SpinChain(8, boundary='periodic'), su2.SU2Chain(8, boundary='periodic')
H_u1, H_su2 = u1.xxz(J=1.0, Delta=1.0), sym.heisenberg(J=1.0)
E_u1 = torch.sort(torch.cat([torch.linalg.eigvalsh(H_u1.matrix[0, N, :k, :k])
                             for N, k in enumerate(u1.algebra.k_factors)])).values
E_su2, S_lev, _ = sym.full_spectrum(H_su2, return_quantum_numbers=True)   # level of block S, 2S+1 times
E_su2, S_lev = E_su2[0], S_lev[0]
assert close(E_u1, E_su2)
assert S_lev[0] == 0 and (S_lev[1:4] == 1).all() and E_su2[4] - E_su2[3] > 1e-3     # singlet, then triplet
assert close(sym.restrict(H_u1).matrix, H_su2.matrix)      # Proposition 4.15(b)
flops = lambda ks: sum(k ** 3 for k in ks)
print(u1.algebra.k_factors, sym.multiplicities, flops(u1.algebra.k_factors), flops(sym.multiplicities))
print(f"ground state energy {E_su2[0].item():.6f}, first excited {E_su2[1].item():.6f} (x3)")
```

Для кольца из восьми спинов блоки $U(1)$ имеют размеры $1,8,28,56,70,56,28,8,1$, блоки $SU(2)$ — размеры
$m_S=14,28,20,7,1$ для $S=0,\dots,4$. Стоимость диагонализации, $\sum k^3$, падает с $739162$ до $33040$, в
$22$ раза; размерность алгебры, $\sum k^2$, — с $\binom{16}8=12870$ до $1430$. Основное состояние — синглет с
энергией $-3.651093$, первый возбуждённый уровень — триплет при $-3.128419$.

> **Физический смысл.** Модель Гейзенберга описывает обменное взаимодействие локализованных
> электронных спинов в магнитных диэлектриках; $J_{ij}>0$ благоприятствует антипараллельным соседям. Её спектр
> состоит из мультиплетов, а магнитное поле $h\,S^z$ расщепляет каждый мультиплет на $2S+1$ эквидистантных
> уровней (упражнение 4). Для констант связи $J_{ij}>0$ между двумя подрешётками двудольной решётки,
> $J_{ij}\le0$ внутри них и подрешёток одинакового размера основное состояние — синглет (теорема Либа–Маттиса,
> E. Lieb and D. Mattis, *J. Math. Phys.* 3 (1962) 749; приводится без доказательства). Кольцо из восьми спинов
> выше — пример этого.

## 4.7 С чем сравнивать

Для конечных групп GAP и SageMath точно вычисляют таблицы характеров, и тогда следствие 4.5 даёт
кратности; RepLAB (MATLAB; D. Rosset, F. Montealegre-Mora, J.-D. Bancal) численно раскладывает
представления конечных и компактных групп, тоже с помощью случайных элементов коммутанта, а
SymbolicWedderburn.jl (Julia) вычисляет разложения Веддерберна для симметрийной редукции
полуопределённых программ. В физике QuSpin и другие программы точной диагонализации приводят к блочно-
диагональному виду $U(1)$ и решёточные симметрии, ITensor и TeNPy работают с абелевыми зарядами в тензорных
сетях, а TensorKit.jl обрабатывает $SU(2)$ с помощью деревьев слияния. Модуль `decompose` библиотеки — плотный
метод для небольших матриц ($n$ до нескольких сотен); он принимает любой список матриц (элементы группы
или образующие), сертифицирует свой результат (предложение 4.9) и возвращает алгебру и её коммутант как
объекты `TypeIAlgebra` с `to_algebra`, `from_algebra`, `to_commutant`, `from_commutant`, так что главы 3, 5–8
применимы сразу. Для цепочек спинов $\frac12$ `su2.SU2Chain` даёт блоки $SU(2)$ напрямую, для всех значений
спина при $L$ до $16$.

## Упражнения

1. Покажите, что для представления $S_L$ на $(\mathbb C^d)^{\otimes L}$
   $\chi_\pi(\sigma)=d^{\#\text{cycles}(\sigma)}$ и $\dim\pi(S_L)'=\binom{d^2+L-1}L$. Проверьте $20$ (раздел 2.4) и
   $35$ ($=\dim\{S_a\}''$ для четырёх спинов).
2. (*Клебш–Гордан.*) Покажите, что $j_1\otimes j_2\cong\bigoplus_jj$, $j=|j_1-j_2|,|j_1-j_2|+1,\dots,j_1+j_2$, для
   спинового представления $S^{(j_1)}_a\otimes1+1\otimes S^{(j_2)}_a$ — неабелево правило слияния $SU(2)$, упомянутое
   в главе 7, — и проверьте $1\otimes\frac12$ и $1\otimes1$ с помощью `from_generators`.
3. Пусть $\theta/2\pi$ иррационально. Покажите, что $\{U_\theta\}'=U(1)'$ для одного-единственного
   $U_\theta=e^{i\theta Q}$, но что $\{e^{i\pi\hat N}\}'\cong M_{2^{n-1}}\oplus M_{2^{n-1}}$ для $n$ кубитов.
4. (*Эффект Зеемана.*) Пусть $H$ — гейзенберговский гамильтониан на $L\ge1$ спинах и $H_h=H+hS_z$. Покажите,
   что $H_h$ коммутирует с $S_z$ и $\mathbf S^2$, но при $h\neq0$ не коммутирует с $S_x$, и что его собственные
   значения равны $E+hM$, $M=-S,\dots,S$, для собственных значений $E$ блока $H_S$. Проверьте это для кольца
   из восьми спинов.

<details>
<summary>Решения</summary>

1. Как и при $d=2$: $\pi(\sigma)$ — матрица перестановки, оставляющая $|b_1\dots b_L\rangle$ на месте тогда и
   только тогда, когда $b$ постоянно на циклах, а таких $b$ ровно $d^{\#\text{cycles}}$. По следствию 4.5(c)
   $\dim\pi(S_L)'=\frac1{L!}\sum_\sigma(d^2)^{\#\text{cycles}(\sigma)}=\frac{x(x+1)\cdots(x+L-1)}{L!}$ с $x=d^2$
   (глава 5, упражнение 2), а это $\binom{d^2+L-1}L$: $\binom63=20$ и $\binom74=35$. То, что $35$ равно также
   $\dim\{S_a\}''$, — это теорема 4.14 для $L=4$.

```python
for n, value in ((3, 20), (4, 35)):
    assert dc.commutant_dimension(dc.permutation_representation(n, 2)) == math.comb(n + 3, n) == value
```

2. Пусть $j_1\ge j_2$. $S_z$ имеет собственное значение $m_1+m_2$ на $|m_1\rangle\otimes|m_2\rangle$, поэтому при
   $M\in j_1+j_2+\mathbb Z$ число $d_M$ равно количеству $m_2\in\{-j_2,\dots,j_2\}$ с $|M-m_2|\le j_1$: все $2j_2+1$,
   если $|M|\le j_1-j_2$, и $j_1+j_2-|M|+1$, если $j_1-j_2\le|M|\le j_1+j_2$ (при $M\ge0$ условие имеет вид
   $m_2\ge M-j_1$), и $0$ за этими пределами. По следствию 4.12 $m_j=d_j-d_{j+1}$ равно $1$ при
   $j_1-j_2\le j\le j_1+j_2$ и $(2j_2+1)-(2j_2+1)=0$ при $j<j_1-j_2$.

```python
def coupled(a, b):                                        # S^(a) (x) 1 + 1 (x) S^(b)
    Sa, Sb = dc.spin_matrices(a), dc.spin_matrices(b)
    ea, eb = torch.eye(Sa.shape[-1], dtype=dt), torch.eye(Sb.shape[-1], dtype=dt)
    return torch.stack([torch.kron(x, eb) + torch.kron(ea, y) for x, y in zip(Sa, Sb)])
for a, b, ks in ((1.0, 0.5, (2, 4)), (1.0, 1.0, (1, 3, 5))):  # 1 (x) 1/2 = 1/2 (+) 3/2, 1 (x) 1 = 0 (+) 1 (+) 2
    D_cg = dc.from_generators(coupled(a, b))
    assert D_cg.k_factors == ks and set(D_cg.multiplicities) == {1}
```

3. $U_\theta=\sum_qe^{iq\theta}P_q$, и $e^{iq\theta}=e^{iq'\theta}$ тогда и только тогда, когда $(q-q')\theta\in2\pi\mathbb Z$,
   что при иррациональном $\theta/2\pi$ влечёт $q=q'$. Поэтому собственные подпространства $U_\theta$ — это
   образы $P_q$, а оператор коммутирует с нормальной матрицей $U_\theta$ тогда и только тогда, когда он сохраняет
   её собственные подпространства (теорема 1.5(d) и интерполяция Лагранжа, как в разделе 2.2), т. е. тогда и
   только тогда, когда $X=\sum_qP_qXP_q$; по предложению 4.10 это $U(1)'$.
   При $\theta=\pi$ собственное значение $(-1)^w$ видит лишь чётность веса Хэмминга: два собственных
   подпространства размерности $2^{n-1}$ (при $n\ge1$; чётные и нечётные подмножества), и коммутант есть
   $M_{2^{n-1}}\oplus M_{2^{n-1}}$.

```python
for theta, m in ((1.0, (1, 1, 4, 4, 6)), (math.pi, (8, 8))):
    D_theta = dc.from_generators(torch.linalg.matrix_exp(1j * theta * N_op))
    assert D_theta.multiplicities == m and set(D_theta.k_factors) == {1}
```

4. $S_z$ и $\mathbf S^2$ коммутируют с $H$ (предложение 4.15(a)) и с $S_z$, а значит, и с $H_h$; и
   $[H_h,S_x]=h[S_z,S_x]=ihS_y\neq0$. В разложении из предложения 4.13
   $V^*S_zV=\bigoplus_SS^{(S)}_z\otimes1_{m_S}$ и $V^*HV=\bigoplus_S1_{2S+1}\otimes H_S$ (предложение 4.15), поэтому на
   слагаемом спина $S$ имеем $V^*H_hV=hS^{(S)}_z\otimes1+1\otimes H_S$ — сумму коммутирующих операторов,
   диагонализуемую произведениями собственных векторов, с собственными значениями $hM+E$.

```python
assert not close(H4 @ Sx - Sx @ H4 + 0.3 * (Sz @ Sx - Sx @ Sz), 0)          # [H4 + 0.3 S_z, S_x] != 0
E_field = sym.full_spectrum(sym.heisenberg(J=1.0, h=0.3))[0]                # E + h M on every multiplet
H_field = u1.xxz(J=1.0, Delta=1.0, h=0.3)
E_u1_field = torch.cat([torch.linalg.eigvalsh(H_field.matrix[0, N, :k, :k])
                        for N, k in enumerate(u1.algebra.k_factors)])
assert close(torch.sort(E_u1_field).values, E_field)
```

</details>

**Далее:** в [главе 5](05_random_matrices.md) случайные матрицы выбираются в алгебре: унитарные матрицы
Хаара группы $\prod_cU(k_c)$, исчисление Вайнгартена — случай второго порядка которого опирается на
коммутант $\{V\otimes V\}$, простейший пример двойственности Шура–Вейля, — и статистика уровней по
секторам.
