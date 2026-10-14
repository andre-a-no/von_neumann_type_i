[English](../12_equivariant_learning.md) · [Español](../es/12_equivariant_learning.md) · [Français](../fr/12_equivariant_learning.md) · [Deutsch](../de/12_equivariant_learning.md) · [中文](../zh/12_equivariant_learning.md) · [日本語](../ja/12_equivariant_learning.md) · **Русский**

# 12. Эквивариантное машинное обучение

Слой нейронной сети называется **эквивариантным**, если преобразовать вход симметрией и затем применить
слой — то же самое, что применить слой и затем преобразовать выход: перенумерация вершин графа точно так же
перенумеровывает его признаки. Для линейных слоёв это коммутационное соотношение, поэтому эквивариантные
слои образуют коммутант, который по главе 2 с точностью до замены базиса есть алгебра $\bigoplus_cM_{m_c}(\mathbb C)$.
В этой главе подсчитываются их параметры, они строятся и обучаются с помощью `torch_vn_algebra.decompose`, а
структура алгебры используется для унитарных, положительных и спектральных слоёв, эквивариантных по построению.

**Вы узнаете**, почему эквивариантные линейные отображения между двумя представлениями — это блочные матрицы
с размерами блоков, заданными кратностями; как подсчитывать их по орбитам или по характерам (DeepSets и числа
Белла для графовых сетей); почему усреднение по группе есть ортогональный проектор на эквивариантные
отображения; как функциональное исчисление и параметризации унитарных операторов действуют на коммутанте
поблочно; и как подсчитываются инвариантные векторы, инвариантные квадратичные формы и инвариантные многочлены
любой степени (формула Молина).

Всюду далее $G$ — конечная группа, а **представление** $G$ на $\mathbb C^n$ — гомоморфизм
$\pi:G\to M_n(\mathbb C)$ с унитарными значениями (это не ограничение: обратимые $\pi(g)$ унитарны относительно скалярного произведения
$\sum_g\langle\pi(g)x,\pi(g)y\rangle$, раздел 4.1). Его **характер** — это $\chi_\pi(g)=\operatorname{Tr}\pi(g)$. Неприводимые и
эквивалентные представления и **кратность** $m_\rho(\pi)$ неприводимого $\rho$ в $\pi$ ($0$, если $\rho$
не встречается) понимаются как в разделе 2.4, где доказано, что $m_\rho(\pi)$ не зависит от разложения, и
доказана лемма Шура; подробно и то и другое разбирается в [главе 4](04_symmetries.md) (разделы 4.1–4.2). На прямоугольных матрицах
$\langle S,T\rangle=\operatorname{Tr}(S^*T)$ — скалярное произведение **Гильберта–Шмидта**. В коде пакет векторов
— это тензор формы `(batch, n)`, на который матрица $T$ действует как `x @ T.T`.

## 12.1 Эквивариантные линейные отображения

**Определение.** Пусть $\pi_1,\pi_2$ — представления $G$ на $\mathbb C^{n_1}$, $\mathbb C^{n_2}$. Матрица
$T\in M_{n_2\times n_1}(\mathbb C)$ называется **эквивариантной** (**сплетающим оператором**), если $T\pi_1(g)=\pi_2(g)T$ для всех
$g\in G$. Такие матрицы образуют векторное пространство $\operatorname{Hom}_G(\pi_1,\pi_2)$; при $\pi_1=\pi_2=\pi$ это
коммутант, $\operatorname{Hom}_G(\pi,\pi)=\pi(G)'$.

Чаще всего в машинном обучении встречаются представления, переставляющие координаты. Пусть $G$ действует на
конечном множестве $\Omega$; **перестановочное представление** на $\mathbb C^\Omega$ — это $\pi(g)e_x=e_{gx}$. Примеры:
$S_n$, переставляющая $n$ точек множества ($\Omega=[n]=\{1,\dots,n\}$), и $S_n$, перенумеровывающая вершины
графа и действующая на $\Omega=[n]^2$, т. е. на матрице смежности по правилу $X\mapsto P_\sigma XP_\sigma^T$.

**Лемма 12.1 (орбиты).** Пусть $\pi_i$ — перестановочные представления $G$ на $\mathbb C^{\Omega_i}$. Тогда
$T\in\operatorname{Hom}_G(\pi_1,\pi_2)$ в том и только в том случае, когда $T_{gy,gx}=T_{yx}$ для всех $g$, $y\in\Omega_2$, $x\in\Omega_1$.
Индикаторные матрицы орбит $G$ на $\Omega_2\times\Omega_1$ образуют базис, и
$\dim\operatorname{Hom}_G(\pi_1,\pi_2)$ равна числу этих орбит.

*Доказательство.* $\pi_2(g)T\pi_1(g)^*e_x=\pi_2(g)Te_{g^{-1}x}=\sum_yT_{y,g^{-1}x}e_{gy}$, поэтому элемент
$\pi_2(g)T\pi_1(g)^*$ с индексами $(gy,x)$ равен $T_{y,g^{-1}x}$. Следовательно, $T$ эквивариантна, т. е. $\pi_2(g)T\pi_1(g)^*=T$
для всех $g$, тогда и только тогда, когда $T_{gy,x}=T_{y,g^{-1}x}$ для всех $g,y,x$, а это (после замены $x$ на $gx$) и
есть сформулированное условие: $T$ постоянна на орбитах. Индикаторные матрицы орбит имеют непересекающиеся
носители, поэтому линейно независимы, а всякая функция, постоянная на орбитах, есть их линейная комбинация. $\square$

**Следствие 12.2.** (a) (*DeepSets.*) Для $S_n$ на $\mathbb C^n$, $n\ge2$, коммутант равен
$\operatorname{span}\{1,J\}$, где $J$ — матрица из одних единиц, а $\mathbb C^n=\mathbb C\mathbf 1\oplus\mathbf 1^\perp$ —
сумма двух неэквивалентных неприводимых представлений: тривиального и **стандартного** размерности
$n-1$. (b) (*Графовые сети.*) Размерность пространства эквивариантных отображений из $\mathbb C^{[n]^p}$ в $\mathbb C^{[n]^q}$ (тензоров
порядков $p$ и $q$ на $n$ вершинах) равна числу разбиений множества из $p+q$ элементов не более чем на $n$
блоков; при $n\ge p+q$ это **число Белла** $B_{p+q}$, число всех разбиений множества из $p+q$ элементов
($B_2=2$, $B_3=5$, $B_4=15$).

*Доказательство.* (b) Два набора $x,y\in[n]^{p+q}$ лежат в одной орбите тогда и только тогда, когда у них одинаковый **шаблон равенств** —
разбиение множества позиций $\{1,\dots,p+q\}$ по правилу $i\sim j\Leftrightarrow x_i=x_j$. Шаблон, очевидно,
инвариантен. Обратно, если у $x,y$ одинаковый шаблон, то $x_i\mapsto y_i$ — корректно определённая инъекция
$\{x_i\}$ в $[n]$, которая продолжается до перестановки $\sigma$ с $\sigma x=y$. Всякое разбиение не более чем на
$n$ блоков является шаблоном (присвоим блокам различные значения). Остаётся применить лемму 12.1. При $n\ge p+q$ ограничение на число
блоков пусто. Подсчёт по числу блоков: $B_2=1+1$; $B_3=1+3+1$ (разбиение $3$-элементного множества на два блока
определяется своим одноэлементным блоком); $B_4=1+7+6+1$ (два блока: $4$ разбиения с размерами $3+1$ и $3$ с размерами
$2+2$; три блока: $\binom42=6$ способов выбрать пару). (a) При $p=q=1$ на $[n]^2$
две орбиты — диагональ и её дополнение — с индикаторами $1$ и $J-1$. Поэтому
$\dim\pi(S_n)'=2=\sum_cm_c^2$ (теорема 2.6, предложение 2.7), откуда следует, что секторов два и в обоих $m_c=1$, т. е. имеются два
неэквивалентных неприводимых слагаемых, и $\pi(S_n)'$ порождается их проекторами $P_1,P_2$. Проектор на
инвариантную прямую $\mathbb C\mathbf 1$ лежит в $\pi(S_n)'$ (лемма 2.1), поэтому он равен $P_1$ или $P_2$: эта прямая,
несущая тривиальное представление, — один сектор, а $\mathbf 1^\perp$ — другой. $\square$

Слой DeepSets $x\mapsto\lambda x+\gamma\,(\mathbf 1^Tx)\mathbf 1$ (Zaheer et al. 2017) — это (a); базисы
Maron et al. (2019) для графовых сетей — это (b). Код подсчитывает орбиты перебором, сравнивает результат с
численно вычисленным коммутантом и раскладывает оба представления $S_4$:

```python
import itertools
import torch
import torch_vn_algebra.decompose as dc

torch.manual_seed(0)
dt = torch.complex128
n = 4
perms = list(itertools.permutations(range(n)))                   # S_4, |G| = 24

def perm_matrix(p):
    """pi(sigma) e_x = e_{sigma(x)} on C^n."""
    P = torch.zeros(n, n, dtype=dt)
    P[list(p), list(range(n))] = 1
    return P

G_node = torch.stack([perm_matrix(p) for p in perms])            # S_4 on C^4: vertex features
G_edge = torch.stack([torch.kron(P, P) for P in G_node])         # S_4 on C^16: X -> P X P^T, row-major

def n_orbits(tuples):
    """Number of orbits of S_4 on tuples of vertices: the canonical form is the smallest image."""
    return len({min(tuple(p[i] for i in t) for p in perms) for t in tuples})

bell = {r: n_orbits(itertools.product(range(n), repeat=r)) for r in (2, 3, 4)}
print(bell)
assert bell == {2: 2, 3: 5, 4: 15}                               # Corollary 12.2(b)
for G, r in ((G_node, 2), (G_edge, 4)):                          # Lemma 12.1: dim pi(G)' = number of orbits
    assert dc.commutant_basis(G, group=True).shape[0] == dc.commutant_dimension(G) == bell[r]
for X in dc.commutant_basis(G_node, group=True):                 # Corollary 12.2(a): X = a 1 + b (J - 1)
    off = X[~torch.eye(n, dtype=torch.bool)]
    assert (X.diagonal() - X[0, 0]).abs().max() < 1e-12 and (off - off[0]).abs().max() < 1e-12

D_node, D_edge = dc.from_group(G_node), dc.from_group(G_edge)
print(D_node)
print(D_edge)
assert (D_node.k_factors, D_node.multiplicities) == ((1, 3), (1, 1))
assert (D_edge.k_factors, D_edge.multiplicities) == ((1, 2, 3, 3), (2, 1, 1, 3))
```

Итак, в разложении, вычисленном `dc.from_group` (его результат проверяют утверждения assert), $\mathbb C^{16}$ содержит
представление размерности $1$ дважды — это тривиальное представление, кратность которого по предложению 12.7(a) и
лемме 12.1 равна числу $2$ орбит на $[4]^2$ (код раздела 12.4 определяет этот сектор), — представление размерности $2$
один раз и два представления размерности $3$ — один раз и три раза: $16=1\cdot2+2\cdot1+3\cdot1+3\cdot3$, а
$15=2^2+1^2+1^2+3^2$ согласуется с подсчётом орбит по лемме 12.1 и с $\dim\pi(G)'=\sum_cm_c^2$ (теорема 2.6). Для двух
представлений всё сводится к одному коммутанту:

**Лемма 12.3 (прямая сумма).** Пусть $\pi=\pi_1\oplus\pi_2$ на $\mathbb C^{n_1}\oplus\mathbb C^{n_2}$ с
изометрическими вложениями $Q_i$ и проекторами $P_i=Q_iQ_i^*$. Тогда $P_1,P_2\in\pi(G)'$,
$\operatorname{Hom}_G(\pi_1,\pi_2)=Q_2^*\,\pi(G)'\,Q_1$, а $T\mapsto Q_2TQ_1^*$ — линейная биекция
$\operatorname{Hom}_G(\pi_1,\pi_2)$ на $P_2\,\pi(G)'P_1$.

*Доказательство.* $\pi(g)Q_i=Q_i\pi_i(g)$ и $Q_i^*\pi(g)=\pi_i(g)Q_i^*$, поскольку $\pi(g)$ блочно-диагональна; следовательно,
$P_i\pi(g)=Q_i\pi_i(g)Q_i^*=\pi(g)P_i$. Если $X\in\pi(G)'$, то
$Q_2^*XQ_1\pi_1(g)=Q_2^*X\pi(g)Q_1=Q_2^*\pi(g)XQ_1=\pi_2(g)Q_2^*XQ_1$. Если $T$ эквивариантна, то
$\pi(g)Q_2TQ_1^*=Q_2\pi_2(g)TQ_1^*=Q_2T\pi_1(g)Q_1^*=Q_2TQ_1^*\pi(g)$, так что $Q_2TQ_1^*\in\pi(G)'$, эта матрица равна
$P_2(Q_2TQ_1^*)P_1$, и $Q_2^*(Q_2TQ_1^*)Q_1=T$. Обратно, всякая $Y=P_2XP_1$ с $X\in\pi(G)'$ лежит в
$\pi(G)'$, и $Y=Q_2(Q_2^*YQ_1)Q_1^*$. $\square$

**Теорема 12.4 (эквивариантные отображения — блочные матрицы).** Пусть $\rho_1,\dots,\rho_C$ — попарно неэквивалентные
неприводимые представления, встречающиеся в $\pi_1\oplus\pi_2$, $k_c=\dim\rho_c$, $a_c=m_{\rho_c}(\pi_1)$,
$b_c=m_{\rho_c}(\pi_2)$. Существуют унитарные операторы $V_1:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{a_c}\to\mathbb C^{n_1}$
и $V_2:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{b_c}\to\mathbb C^{n_2}$ с
$V_i^*\pi_i(g)V_i=\bigoplus_c\rho_c(g)\otimes1$, такие что
$$
\operatorname{Hom}_G(\pi_1,\pi_2)=\Big\{V_2\Big(\bigoplus_{c=1}^C1_{k_c}\otimes y_c\Big)V_1^*:\ y_c\in M_{b_c\times a_c}(\mathbb C)\Big\},
$$
и отображение $(y_c)_c\mapsto T$ инъективно. Следовательно, $\operatorname{Hom}_G(\pi_1,\pi_2)\cong\bigoplus_cM_{b_c\times a_c}(\mathbb C)$,
$$
\dim\operatorname{Hom}_G(\pi_1,\pi_2)=\sum_\rho m_\rho(\pi_1)\,m_\rho(\pi_2),
$$
где сумма берётся по всем неприводимым представлениям с точностью до эквивалентности; при $\pi_1=\pi_2$ это
$\pi(G)'\cong\bigoplus_cM_{m_c}$ (предложение 2.7), и $\operatorname{Hom}_G(\pi_1,\pi_2)=0$ тогда и только тогда, когда никакое
неприводимое представление не встречается в обоих.

*Доказательство.* Применим теорему 2.6 и предложение 2.7 к $\pi=\pi_1\oplus\pi_2$: существует унитарный $V$ с
$V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$ и $V^*\pi(G)'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$. В силу единственности
кратностей (раздел 2.4) всякое неприводимое представление, встречающееся в $\pi$, эквивалентно ровно одному
$\rho_c$, так что $\rho_c$ — представители классов из формулировки. Любые другие представители имеют вид
$\rho_c'=w_c\rho_cw_c^*$ с унитарными $w_c$ (лемма 4.2(c)), и замена $V$ на $V\bigoplus_c(w_c^*\otimes1)$ превращает
первую нормальную форму в $\bigoplus_c\rho'_c(g)\otimes1_{m_c}$ и не меняет вторую, так как $w_c^*\otimes1$
коммутирует с $1\otimes M_{m_c}$; поэтому можно использовать $\rho_c$ из предложения 2.7. Так как
$X\mapsto V^*XV$ — $*$-изоморфизм и $P_1\in\pi(G)'$ (лемма 12.3), имеем $V^*P_1V=\bigoplus_c1\otimes p_c$ с
проекторами $p_c$. Выберем унитарные $u_c$ с $u_c^*p_cu_c=1_{a_c}\oplus0_{b_c}$ (теорема 1.3), где пока
$a_c=\operatorname{rank}p_c$, $b_c=m_c-a_c$, и заменим $V$ на $V\bigoplus_c(1\otimes u_c)$; это не меняет ни одной из
нормальных форм, так как $1\otimes u_c$ коммутирует с $\rho_c(g)\otimes1$. Пусть $\iota_c:\mathbb C^{a_c}\to\mathbb C^{m_c}$ и
$\kappa_c:\mathbb C^{b_c}\to\mathbb C^{m_c}$ — вложения первых $a_c$ и последних $b_c$ координат, и
$$
V_1=Q_1^*V\bigoplus_c(1\otimes\iota_c),\qquad V_2=Q_2^*V\bigoplus_c(1\otimes\kappa_c).
$$
$W_1=V\bigoplus_c(1\otimes\iota_c)$ — изометрия с образом $\operatorname{range}P_1=\operatorname{range}Q_1$, поэтому
$V_1^*V_1=W_1^*P_1W_1=1$, и $V_1$ унитарен, так как $\sum_ck_ca_c=\operatorname{rank}P_1=n_1$. Кроме того,
$\pi_1(g)V_1=Q_1^*\pi(g)W_1=Q_1^*V\big(\bigoplus_c\rho_c(g)\otimes1_{m_c}\big)\bigoplus_c(1\otimes\iota_c)
=V_1\bigoplus_c\rho_c(g)\otimes1_{a_c}$. Значит, $\pi_1\cong\bigoplus_c\rho_c\otimes1_{a_c}$ с неэквивалентными
неприводимыми $\rho_c$, и $a_c=m_{\rho_c}(\pi_1)$ в силу единственности кратностей (раздел 2.4); аналогично для
$V_2$ и $b_c$. Наконец, для $X=V\big(\bigoplus_c1\otimes x_c\big)V^*\in\pi(G)'$
$P_2XP_1=V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*$, где $y_c=\kappa_c^*x_c\iota_c$ — левый нижний
угол $x_c$ размера $b_c\times a_c$, произвольная матрица. По лемме 12.3
$\operatorname{Hom}_G(\pi_1,\pi_2)=Q_2^*P_2\,\pi(G)'P_1Q_1$ состоит из матриц
$Q_2^*V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*Q_1=V_2\big(\bigoplus_c1\otimes y_c\big)V_1^*$, и
отображение $(y_c)\mapsto T$ инъективно, поскольку $V_1,V_2$ унитарны. Для неприводимого $\rho$, не встречающегося в
$\pi_1\oplus\pi_2$, $m_\rho(\pi_1)m_\rho(\pi_2)=0$. $\square$

Для следствия 12.2(b) при $n=4$ в вычисленном разложении (оно печатается и проверяется ниже): признаки вершин $\pi_1$
имеют $a=(1,0,0,1)$ в четырёх секторах $\pi_1\oplus\pi_2$ (размерностей $1,2,3,3$; первый — тривиальный, поскольку по
следствию 12.2(a) единственная одномерная составляющая $\pi_1$ тривиальна), признаки рёбер $\pi_2$ имеют
$b=(2,1,1,3)$, так что $\dim\operatorname{Hom}_G=1\cdot2+1\cdot3=5=B_3$, как и доказывает следствие 12.2(b).

> **Физический смысл** (физическая интерпретация, а не теорема этого курса). Теорема 12.4 — линейно-алгебраическое ядро **правил отбора** и скалярный случай
> теоремы Вигнера–Эккарта: оператор, коммутирующий с симметрией, может связывать лишь состояния одного и того же
> типа симметрии, а внутри типа действует как $1_{k_c}\otimes y_c$, т. е. через «приведённый матричный элемент» $y_c$,
> не чувствительный к положению внутри мультиплета. (Операторам, которые сами преобразуются по нетривиальному
> представлению, например векторным, нужны ещё коэффициенты Клебша–Гордана.) Симметричный гамильтониан —
> это случай $\pi_1=\pi_2$; его вырождения $k_c$ равны размерностям неприводимых представлений.

## 12.2 Эквивариантный слой средствами библиотеки

Доказательство теоремы 12.4 является алгоритмом, коль скоро известно разложение из теоремы 2.6; `dc.from_group`
вычисляет его рандомизированным методом из раздела 4.3. При $\pi_1=\pi_2$ это просто теорема 2.6: веса — это блоки
$w_c\in M_{m_c}$, т. е. оператор абстрактного коммутанта `D.commutant_algebra()` $=\bigoplus_cM_{m_c}$, а
слой — это $x\mapsto V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*x$, вычисляемый `D.from_commutant`. Модуль
ниже — код главы, а не часть библиотеки.

```python
class EquivariantLinear(torch.nn.Module):
    """x -> T x, T = V ((+)_c 1_{k_c} (x) w_c) V^* in pi(G)'; one weight block w_c of size m_c x m_c per sector."""
    def __init__(self, D):
        super().__init__()
        self.D, self.CA = D, D.commutant_algebra()
        self.w = torch.nn.ParameterList(torch.nn.Parameter(torch.randn(m, m, dtype=dt) / m ** 0.5)
                                        for m in D.multiplicities)

    def matrix(self):
        return self.D.from_commutant(self.CA.from_blocks(list(self.w)))[0]

    def forward(self, x):                                        # x: (batch, n)
        return x @ self.matrix().T

def equivariance_error(layer, G_in, G_out):
    """max_g |layer(pi_in(g) x) - pi_out(g) layer(x)| on a random batch."""
    x = torch.randn(8, G_in.shape[-1], dtype=dt)
    with torch.no_grad():
        return max((layer(x @ g.T) - layer(x) @ h.T).abs().max().item() for g, h in zip(G_in, G_out))

edge_layer = EquivariantLinear(D_edge)
n_params = sum(p.numel() for p in edge_layer.parameters())
err = equivariance_error(edge_layer, G_edge, G_edge)
print(f"edge -> edge: {n_params} parameters, equivariance error {err:.1e}")
assert n_params == D_edge.commutant_dim == bell[4] and err < 1e-12
```

При $\pi_1\ne\pi_2$ доказательство воспроизводится строка за строкой: разложить $\pi_1\oplus\pi_2$, вычислить блоки $p_c$
проектора $P_1$ с помощью `to_commutant`, диагонализовать их и оставить левые нижние углы. След $p_c$ равен
$a_c=m_c(\pi_1)$.

```python
class EquivariantMap(torch.nn.Module):
    """x -> T x with T in Hom_G(pi_1, pi_2) (Theorem 12.4): one weight block y_c of size b_c x a_c per sector."""
    def __init__(self, G1, G2):
        super().__init__()
        self.n1 = n1 = G1.shape[-1]
        self.D = dc.from_group(torch.stack([torch.block_diag(g1, g2) for g1, g2 in zip(G1, G2)]))
        P1 = torch.block_diag(torch.eye(n1, dtype=dt), torch.zeros(G2.shape[-1], G2.shape[-1], dtype=dt))
        p = self.D.to_commutant(P1).matrix[0]                    # blocks p_c of P_1, padded
        self.u, self.a = [], []
        for c, m in enumerate(self.D.multiplicities):
            lam, u = torch.linalg.eigh(p[c, :m, :m])             # eigenvalues 0, ..., 0, 1, ..., 1
            self.a.append(round(lam.sum().item()))               # a_c = rank p_c = m_c(pi_1)
            self.u.append(u.flip(-1))                            # u_c^* p_c u_c = 1_{a_c} (+) 0_{b_c}
        self.y = torch.nn.ParameterList(torch.nn.Parameter(torch.randn(m - a, a, dtype=dt))
                                        for m, a in zip(self.D.multiplicities, self.a))

    def matrix(self):
        blocks = []
        for u, a, y in zip(self.u, self.a, self.y):
            x = torch.zeros(u.shape, dtype=dt)
            x[a:, :a] = y                                        # kappa_c y_c iota_c^*
            blocks.append(u @ x @ u.mH)
        T = self.D.from_commutant(self.D.commutant_algebra().from_blocks(blocks))[0]
        return T[self.n1:, :self.n1]                             # Q_2^* (...) Q_1

    def forward(self, x):
        return x @ self.matrix().T

node_to_edge = EquivariantMap(G_node, G_edge)
n_params = sum(p.numel() for p in node_to_edge.parameters())
err = equivariance_error(node_to_edge, G_node, G_edge)
print(f"node -> edge: sectors k = {node_to_edge.D.k_factors}, a = {node_to_edge.a}, "
      f"{n_params} parameters, equivariance error {err:.1e}")
assert n_params == bell[3] and err < 1e-12
corners = node_to_edge.D.basis('commutant')[:, n:, :n].reshape(-1, 16 * n)   # Lemma 12.3: Q_2^* pi(G)' Q_1
assert torch.linalg.matrix_rank(corners).item() == bell[3]
```

Второй способ сделать слой эквивариантным — усреднить его по группе. **Твирлинг** — это
$$
\mathcal P(W)=\frac1{|G|}\sum_{g\in G}\pi_2(g)\,W\,\pi_1(g)^*,\qquad W\in M_{n_2\times n_1}(\mathbb C).
$$

**Теорема 12.5 (твирлинг — ортогональный проектор).** $\mathcal P$ — ортогональный относительно скалярного произведения Гильберта–Шмидта
проектор $M_{n_2\times n_1}(\mathbb C)$ на $\operatorname{Hom}_G(\pi_1,\pi_2)$, и
$$
\dim\operatorname{Hom}_G(\pi_1,\pi_2)=\frac1{|G|}\sum_{g\in G}\overline{\chi_{\pi_1}(g)}\,\chi_{\pi_2}(g).
$$

*Доказательство.* Для $h\in G$ $\pi_2(h)\mathcal P(W)\pi_1(h)^*=\frac1{|G|}\sum_g\pi_2(hg)W\pi_1(hg)^*=\mathcal P(W)$, так как
$g\mapsto hg$ — биекция $G$; поэтому $\mathcal P(W)$ эквивариантна. Если $W$ эквивариантна, то каждое слагаемое равно
$W$, так что $\mathcal P(W)=W$. Следовательно, $\mathcal P^2=\mathcal P$, и образ равен $\operatorname{Hom}_G(\pi_1,\pi_2)$. Отображение
самосопряжено: $\langle S,\mathcal P(W)\rangle=\frac1{|G|}\sum_g\operatorname{Tr}\big((\pi_2(g)^*S\pi_1(g))^*W\big)$
$=\langle\mathcal P(S),W\rangle$, где мы сделали замену $g\mapsto g^{-1}$ и воспользовались $\pi_i(g^{-1})=\pi_i(g)^*$.
Самосопряжённый идемпотент есть ортогональный проектор на свой образ. Его след равен его рангу. Линейное отображение
$W\mapsto AWB$ на $M_{n_2\times n_1}$ имеет след $\operatorname{Tr}A\operatorname{Tr}B$ (коэффициент при
матричной единице $E_{ij}$ в $AE_{ij}B$ равен $A_{ii}B_{jj}$); при $A=\pi_2(g)$, $B=\pi_1(g)^*$ это даёт формулу. $\square$

Вместе с теоремой 12.4 получаем $\sum_\rho m_\rho(\pi_1)m_\rho(\pi_2)=|G|^{-1}\sum_g\overline{\chi_{\pi_1}(g)}\chi_{\pi_2}(g)$;
при $\pi_1=\pi_2$ это формула для $\dim\pi(G)'$ через характеры из раздела 4.2, которую
вычисляет `dc.commutant_dimension`. Параметризация $W\mapsto\mathcal P(W)$ достигает всякого
эквивариантного отображения, $\mathcal P(W)$ — ближайшее к $W$ эквивариантное отображение, а так как $\mathcal P$ самосопряжён,
градиент по $W$ есть $\mathcal P$ от градиента по $T=\mathcal P(W)$. Но
$n_1n_2-\dim\operatorname{Hom}_G$ параметров тратятся впустую на ядро ($256-15$ для признаков рёбер), а одно
вычисление стоит $O(|G|n_1n_2(n_1+n_2))$ для плотных $\pi(g)$, причём $|G|=n!$ для $S_n$. Библиотека вычисляет тот же
проектор по разложению, `D.conditional_expectation(W, 'commutant')`, за четыре произведения матриц $n\times n$.

```python
W = torch.randn(16, 16, dtype=dt)
PW = dc.twirl(G_edge, W)
S = torch.randn(16, 16, dtype=dt)
hs = lambda A, B: (A.conj() * B).sum()                           # Hilbert-Schmidt inner product
T_layer = edge_layer.matrix().detach()
assert (PW - D_edge.conditional_expectation(W, 'commutant')[0]).abs().max() < 1e-12
assert (dc.twirl(G_edge, T_layer) - T_layer).abs().max() < 1e-12           # P = id on pi(G)'
assert abs(hs(S, PW) - hs(dc.twirl(G_edge, S), W)) < 1e-10                  # P is self-adjoint
assert abs(hs(T_layer, W - PW)) < 1e-10                                     # W - P(W) is orthogonal to pi(G)'
chi_node, chi_edge = (G.diagonal(dim1=-2, dim2=-1).sum(-1) for G in (G_node, G_edge))
assert abs((chi_node.conj() * chi_edge).mean() - bell[3]) < 1e-12           # Theorem 12.5: 5
```

Слои выше линейны. Сеть чередует их с нелинейностями, и те тоже должны быть эквивариантными.
Для перестановочного представления функция, применяемая к каждой координате, эквивариантна, поскольку
$(\pi(g)x)_{gy}=x_y$; для других представлений это не так, и используют нелинейности, действующие на нормы или на
коэффициенты в базисе, который $G$ переставляет (упражнение 3).

> **Физический смысл.** Твирлинг — это условное ожидание на коммутант (упражнение 3 главы 2):
> усреднение взаимодействия по всем операциям симметрии оставляет в точности его симметричную часть.
> «Симметризация» гамильтониана или взаимодействия суммированием по группе симметрий — тот же самый проектор.

## 12.3 Слои с ограничениями и спектральные слои

Веса $w=(w_c)_c$ слоя `EquivariantLinear` живут в $*$-алгебре, поэтому к ним применима любая конструкция глав 1 и 11.
Положим $\Phi(w)=V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*$ для $w\in\bigoplus_cM_{m_c}$.

**Предложение 12.6.** (a) $\Phi$ — $*$-изоморфизм $\bigoplus_cM_{m_c}$ на $\pi(G)'$, и $\Phi(w)$
самосопряжён, унитарен или положителен тогда и только тогда, когда таков каждый $w_c$. (b) Для нормального $w$ и функции $f$ на
$\bigcup_c\operatorname{spec}w_c$ $f(\Phi(w))=\Phi\big((f(w_c))_c\big)$; в частности, $f(X)\in\pi(G)'$ для всякого
нормального $X\in\pi(G)'$. (c) Эквивариантные унитарные операторы — это в точности $\Phi(u_0e^{s-s^*})$, где $s$ пробегает
$\bigoplus_cM_{m_c}$, при любом фиксированном унитарном $u_0\in\bigoplus_cM_{m_c}$; в частности, это в точности экспоненты $e^S$
антисамосопряжённых $S\in\pi(G)'$.

*Доказательство.* (a) $\Phi$ линейно, унитально, мультипликативно и сохраняет $*$, поскольку $(1\otimes x)(1\otimes y)=1\otimes xy$
и $(1\otimes x)^*=1\otimes x^*$; оно инъективно и по теореме 2.6 сюръективно на $\pi(G)'$. Следовательно, $\Phi(w)^*=\Phi(w)$ тогда и только тогда, когда
$w^*=w$, и $\Phi(w)^*\Phi(w)=1=\Phi(w)\Phi(w)^*$ тогда и только тогда, когда $w^*w=1=ww^*$. Если все $w_c=t_c^*t_c$, то
$\Phi(w)=\Phi(t)^*\Phi(t)\ge0$ (следствие 1.7). Обратно, если $\Phi(w)\ge0$, то для $v\in\mathbb C^{m_c}$ и
вектора $\xi=V(e_1\otimes v)$ в секторе $c$ имеем $0\le\langle\xi,\Phi(w)\xi\rangle=\langle v,w_cv\rangle$.
(b) Пусть $w_c=u_c\Lambda_cu_c^*$ (теорема 1.3). Тогда $\Phi(w)=U\big(\bigoplus_c1\otimes\Lambda_c\big)U^*$ с
унитарным $U=V\bigoplus_c(1\otimes u_c)$ — диагонализация, в которой собственные значения $w_c$ повторены $k_c$ раз; по
определению и теореме 1.5(a) $f(\Phi(w))=U\big(\bigoplus_c1\otimes f(\Lambda_c)\big)U^*=\Phi\big((f(w_c))_c\big)$.
(c) По (a) эквивариантные унитарные операторы — это $\Phi(u)$ с унитарным $u$ из $\bigoplus_cM_{m_c}$, а такие $u$ — это в точности
$u_0e^{s-s^*}$ по лемме 11.1. При $u_0=1$ из (b) получаем $\Phi(e^{s-s^*})=e^{\Phi(s-s^*)}$, а $\Phi(s-s^*)$ пробегает
антисамосопряжённые элементы $\pi(G)'$. $\square$

Итак, ортогональный (унитарный) эквивариантный слой — это `optimize.UnitaryParam` на `D.commutant_algebra()` с последующим
$\Phi$; по (c) эта параметризация достигает всякого эквивариантного унитарного оператора, и каждая итерация градиентного метода
точно унитарна и точно эквивариантна. Аналогично `PositiveParam` даёт положительные эквивариантные слои, нормированные как
в главе 11 (спектр в $[\ell,1]$ с максимумом $1$; по (a) общий положительный слой имеет вид $\Phi(t^*t)$), а
спектральный слой $f(L)$ самосопряжённого эквивариантного генератора $L$, например диффузия $e^{-tL}$, требует
спектральных разложений стоимостью $\sum_cm_c^3$ (здесь $37$) вместо $n^3=4096$, с градиентами по
формуле Далецкого–Крейна (`apply_function`). Код подгоняет эквивариантный унитарный оператор под $32$ пары вход–выход и сверяет (b)
и её производную по $t$, $\frac d{dt}e^{-tL}=-Le^{-tL}$, с плотными матрицами.

```python
from torch_vn_algebra import optimize as opt

CA = D_edge.commutant_algebra()                                  # M_2 (+) M_1 (+) M_1 (+) M_3
U_target = D_edge.from_commutant(CA.random_unitary_operator(1))[0]
X = torch.randn(32, 16, dtype=dt)
Y = X @ U_target.T

U = opt.UnitaryParam(CA, 1)                                      # u_0 Haar-random, s = 0
adam = torch.optim.Adam(U.parameters(), lr=0.05)
for step in range(150):
    U_full = D_edge.from_commutant(U.operator())[0]
    loss = (X @ U_full.T - Y).abs().pow(2).mean()
    if step == 0:
        loss0 = loss.item()
    adam.zero_grad()
    loss.backward()
    adam.step()
U_full = D_edge.from_commutant(U.operator())[0].detach()
unit_err = (U_full.mH @ U_full - torch.eye(16, dtype=dt)).abs().max().item()
equi_err = max((g @ U_full - U_full @ g).abs().max().item() for g in G_edge)
print(f"loss {loss0:.3f} -> {loss.item():.1e};  |U*U - 1| = {unit_err:.1e}, |[U, pi(g)]| = {equi_err:.1e}")
assert loss.item() < 1e-3 * loss0 and unit_err < 1e-12 and equi_err < 1e-12

h = torch.randn(1, CA.C, CA.k_max, CA.k_max, dtype=dt)
h = CA.operator(h + h.mH, is_self_adjoint=True)                  # a self-adjoint generator h = (h_c)
L = D_edge.from_commutant(h)[0]
t = torch.tensor(0.7, dtype=torch.float64, requires_grad=True)
E = D_edge.from_commutant(h.apply_function(lambda w: torch.exp(-t * w)))[0]   # e^{-tL}, in the small algebra
Z = torch.randn(16, 16, dtype=dt)
(Z.conj() * E).sum().real.backward()                             # d/dt Re <Z, e^{-tL}>
dE = -L @ E.detach()
assert (E - torch.linalg.matrix_exp(-0.7 * L)).abs().max() < 1e-12          # Proposition 12.6(b)
assert abs(t.grad.item() - (Z.conj() * dE).sum().real.item()) < 1e-10
```

> **Физический смысл.** Симметричный гамильтониан $H\in\pi(G)'$ порождает симметричную эволюцию во времени
> $e^{-itH}\in\pi(G)'$ (глава 8), которая никогда не смешивает секторы и действует только на пространствах кратности.
> Эквивариантный унитарный слой — это такая эволюция, но не заданная, а выученная.

## 12.4 Инварианты

Вектор $x$ называется **инвариантным**, если $\pi(g)x=x$ для всех $g$; функция $q$ на $\mathbb C^n$ инвариантна, если
$q(\pi(g)x)=q(x)$. Инвариантные признаки — это то, что может выдавать инвариантная сеть (скажем, классификатор графов).

**Предложение 12.7.** (a) Инвариантные векторы образуют образ ортогонального проектора
$P_0=|G|^{-1}\sum_g\pi(g)$; их размерность равна $m_{\mathrm{triv}}(\pi)=|G|^{-1}\sum_g\chi_\pi(g)$, и если тривиальное
представление встречается — как $\rho_{c_0}$, — то $P_0$ есть минимальный центральный проектор $P_{c_0}$ алгебры $\pi(G)''$.
(b) Эрмитова форма $q_A(x)=\langle x,Ax\rangle$, $A=A^*$, инвариантна тогда и только тогда, когда $A\in\pi(G)'$; инвариантные эрмитовы
формы образуют вещественное векторное пространство размерности $\sum_cm_c^2$.
(c) Однородные инвариантные многочлены степени $2$ от координат $x_1,\dots,x_n$ образуют пространство размерности
$$
\frac1{|G|}\sum_{g\in G}\frac{\chi_\pi(g)^2+\chi_\pi(g^2)}2 .
$$
Для перестановочного представления это число орбит $G$ на мультимножествах $\{x,y\}$ из двух элементов $\Omega$.

*Доказательство.* (a) Инвариантные векторы — это $\operatorname{Hom}_G(\mathrm{triv},\pi)$, где $\mathrm{triv}$ — тривиальное
представление на $\mathbb C=M_{1\times1}$; для него твирлинг равен $\mathcal P(x)=P_0x$, так что теорема 12.5 даёт
проектор и размерность, а теорема 12.4 даёт $m_{\mathrm{triv}}(\pi)$. В базисе теоремы 2.6 запишем
инвариантный вектор как $\xi=\sum_c\xi_c$ с $\xi_c=\sum_jv_{c,j}\otimes e_j\in\mathbb C^{k_c}\otimes\mathbb C^{m_c}$.
Так как $\pi(g)$ сохраняет секторы, $\rho_c(g)v_{c,j}=v_{c,j}$ для всех $g$. Ненулевой $v_{c,j}$ порождает инвариантную
прямую, несущую тривиальное представление, поэтому тогда $\rho_c$, будучи неприводимым, тривиально, и $c=c_0$. Следовательно,
инвариантные векторы — это векторы сектора $c_0$, где $\pi(g)$ действует как $1$: образ $P_{c_0}$, который
поэтому и есть $P_0$.
(b) $q_A(\pi(g)x)=\langle x,\pi(g)^*A\pi(g)x\rangle$, а матрица $D$ с $\langle x,Dx\rangle=0$ для всех $x$ равна нулю
(поляризация, доказательство следствия 1.4). Поэтому $q_A$ инвариантна тогда и только тогда, когда $\pi(g)^*A\pi(g)=A$ для всех $g$, т. е. $A\in\pi(G)'$.
По предложению 12.6(a) самосопряжённые элементы $\pi(G)'$ — это $\Phi(w)$ с эрмитовыми $w_c$; они образуют вещественное пространство
размерности $\sum_cm_c^2$.
(c) Квадратичный многочлен $p(x)=x^TBx$ определяется своей симметричной матрицей $B$. Положим $(g\cdot p)(x)=p(\pi(g)^{-1}x)$;
тогда $g\cdot(h\cdot p)=(gh)\cdot p$, и у $g\cdot p$ матрица $\sigma(g)(B)=A^TBA$, $A=\pi(g)^{-1}$.
Среднее $M=|G|^{-1}\sum_g\sigma(g)$ удовлетворяет $\sigma(h)M=M$ и $MB=B$ для инвариантных $B$, поэтому (как в доказательстве
теоремы 12.5) оно идемпотентно с образом, состоящим из инвариантных многочленов, а размерность их пространства равна
$\operatorname{Tr}M=|G|^{-1}\sum_g\operatorname{Tr}\sigma(g)$. Матрица $A$ унитарна, следовательно, $A=R\Lambda R^{-1}$ с
$\Lambda=\operatorname{diag}(\lambda)$ (теорема 1.3). Симметричные матрицы
$B_{ij}=(R^{-1})^T(e_ie_j^T+e_je_i^T)R^{-1}$, $i\le j$, образуют базис, и
$A^TB_{ij}A=(R^{-1})^T\Lambda(e_ie_j^T+e_je_i^T)\Lambda R^{-1}=\lambda_i\lambda_jB_{ij}$. Поэтому
$\operatorname{Tr}\sigma(g)=\sum_{i\le j}\lambda_i\lambda_j=\frac12\big((\operatorname{Tr}A)^2+\operatorname{Tr}A^2\big)
=\frac12\big(\chi_\pi(g^{-1})^2+\chi_\pi(g^{-2})\big)$, и суммирование по $g^{-1}$ вместо $g$ даёт формулу.
Для перестановочного представления $G$ переставляет мономы $x_ax_b$, поэтому многочлен инвариантен тогда и только тогда, когда его
коэффициенты постоянны на орбитах мономов. $\square$

Для графов инвариантные признаки рёбер порождаются индикаторами диагонали и её дополнения
($B_2=2$), а инвариантные линейные считывания — это суммы диагональных и внедиагональных элементов. Код
находит тривиальный сектор по блокам $\rho_c(g)$ (`to_algebra`) и проверяет (a), (b) и оба подсчёта из (c).

```python
rho = D_edge.to_algebra(G_edge).matrix                           # blocks rho_c(g), shape (24, C, 3, 3)
c0 = [c for c, k in enumerate(D_edge.k_factors) if k == 1 and (rho[:, c, 0, 0] - 1).abs().max() < 1e-12]
P0 = G_edge.mean(0)
assert c0 == [0] and (P0 - D_edge.central_projections()[0]).abs().max() < 1e-12      # (a)
assert abs(chi_edge.mean() - D_edge.multiplicities[0]) < 1e-12
diag = torch.eye(n, dtype=dt).reshape(16)
assert (P0 @ diag - diag).abs().max() < 1e-12 and (P0 @ (1 - diag) - (1 - diag)).abs().max() < 1e-12

w = [torch.randn(m, m, dtype=dt) for m in D_edge.multiplicities]
A = D_edge.from_commutant(CA.from_blocks([b + b.mH for b in w]))[0]                   # (b)
x = torch.randn(16, dtype=dt)
q = lambda v: (v.conj() @ A @ v).real
assert (A - A.mH).abs().max() < 1e-12
assert max(abs(q(g @ x) - q(x)).item() for g in G_edge) < 1e-12

def degree2(G):                                                  # (c): characters
    chi, chi_sq = (H.diagonal(dim1=-2, dim2=-1).sum(-1) for H in (G, G @ G))
    return ((chi ** 2 + chi_sq) / 2).mean().real.item()

pairs = list(itertools.product(range(n), repeat=2))
multisets = [(a, b) for i, a in enumerate(pairs) for b in pairs[i:]]
orbits2 = len({min(tuple(sorted(tuple(p[v] for v in e) for e in (a, b))) for p in perms) for a, b in multisets})
print(f"invariant quadratic polynomials: vertices {degree2(G_node):.4f}, edges {degree2(G_edge):.4f} "
      f"(orbits of S_4 on monomials: {orbits2})")
assert abs(degree2(G_node) - 2) < 1e-12 and abs(degree2(G_edge) - orbits2) < 1e-12 and orbits2 == 11
```

На $\mathbb C^4$ это $\sum_ix_i^2$ и $\sum_{i<j}x_ix_j$. Все степени подсчитываются формулой Молина
$\sum_d\dim\mathbb C[x]^G_dt^d=|G|^{-1}\sum_g\det(1-t\pi(g))^{-1}$, коэффициент которой при $t^2$ есть (c) (приводится без
доказательства; B. Sturmfels, *Algorithms in Invariant Theory*, Springer 2008, теорема 2.2.1).

> **Физический смысл.** В теории фазовых переходов Ландау свободная энергия — инвариантный многочлен от
> вещественного параметра порядка $x$, на котором $G$ действует вещественными ортогональными матрицами. Её квадратичная часть — это $x^TBx$ с $B$ из
> коммутанта (так как $\pi(g)^T=\pi(g)^{-1}$), поэтому по предложению 12.6(a) она блочно-диагональна по секторам:
> симметричная фаза теряет устойчивость в том секторе, где собственное значение $B$ первым меняет знак, и (в общем случае)
> неприводимое представление этого сектора определяет симметрию упорядоченной фазы.

## 12.5 С чем сравнивать

Эквивариантные слои для конкретных групп предоставляются зрелыми библиотеками. **escnn** (Cesa, Lang и Weiler 2022,
преемник e2cnn) строит управляемые (steerable) свёртки для подгрупп $O(2)$ и $O(3)$ из неприводимых
представлений и гармонических базисов ядер; **e3nn** (Geiger и Smidt 2022) строит $E(3)$-эквивариантные сети
из сферических гармоник и тензорных произведений Клебша–Гордана. **EMLP** (Finzi, Welling и Wilson 2021) работает с любой
матричной группой: он налагает $\pi_2(g)T=T\pi_1(g)$ для дискретных образующих и $d\pi_2(X)T=Td\pi_1(X)$ для образующих
алгебры Ли и вычисляет базис решений как численное ядро — это линейная задача
`dc.commutant_basis`. Групповые свёртки (Cohen и Welling 2016) — случай регулярного представления,
а Kondor и Trivedi (2018) показывают, что для отображений признаков на однородных пространствах компактной группы эквивариантные
линейные слои являются обобщёнными свёртками. Для перестановок DeepSets (Zaheer et al. 2017) и базисы
Maron et al. (2019) — это следствие 12.2.

Библиотека не предоставляет свёрток на решётках, управляемых ядер на непрерывных пространствах, разложений над
$\mathbb R$ (вещественный, комплексный и кватернионный типы; упражнение 2 рассматривает вещественные матрицы) и не рассчитана на большие $n$: $V$ — плотная
матрица $n\times n$, и разложение стоит $O(|G|n^3)$ на одну случайную выборку для группы, заданной всеми своими элементами
(образующие алгебры Ли обрабатываются через `dc.from_generators`; раздел 4.3). Что она добавляет, так это
**структуру** пространства решений: базис коммутанта — это лишь список матриц, тогда как разложение
отождествляет коммутант с $*$-алгеброй $\bigoplus_cM_{m_c}$. Произведения, сопряжения, унитарность, положительность и
функциональное исчисление тогда действуют поблочно на маленьких матрицах, с градиентами (предложение 12.6), а
`to_commutant`/`from_commutant` переводят из одной картины в другую.

## Упражнения

1. (**DeepSets с каналами.**) $S_n$ ($n\ge2$) действует на $\mathbb C^n\otimes\mathbb C^d$ посредством $P_\sigma\otimes1_d$. Покажите,
   что эквивариантные отображения из $\mathbb C^n\otimes\mathbb C^{d_1}$ в $\mathbb C^n\otimes\mathbb C^{d_2}$ — это в точности
   $1\otimes A+J\otimes B$ с $A,B\in M_{d_2\times d_1}$, и проверьте это с помощью `EquivariantMap`.
2. (**Вещественные слои.**) Пусть матрицы $\pi_1,\pi_2$ вещественны (например, у перестановочных представлений). Покажите, что вещественные
   эквивариантные матрицы образуют вещественное векторное пространство размерности $\dim_{\mathbb C}\operatorname{Hom}_G(\pi_1,\pi_2)$ и
   что $T\mapsto\operatorname{Re}T$ отображает $\operatorname{Hom}_G(\pi_1,\pi_2)$ на него.
3. (**Нелинейности.**) (a) Покажите, что применение функции $f:\mathbb C\to\mathbb C$ к каждой координате коммутирует
   с любым перестановочным представлением. (b) Покажите, что при $f(t)=\max(t,0)$ оно не коммутирует с поворотом
   $\mathbb R^2$ на $\pi/2$. (c) Покажите, что $x\mapsto\varphi(\|x\|)x$ коммутирует с любым унитарным
   представлением.
4. (**Сколько выборок определяют слой?**) Пусть $T=\Phi(w)\in\pi(G)'$ и $x_1,\dots,x_N\in\mathbb C^n$. Запишем
   $V^*x_i$ в секторе $c$ как $\sum_{a,j}(X_{i,c})_{aj}e_a\otimes e_j$ с $X_{i,c}\in M_{k_c\times m_c}$, и пусть
   $S_c\in M_{Nk_c\times m_c}$ — матрица, составленная из $X_{1,c},\dots,X_{N,c}$, поставленных друг под другом. Покажите, что эквивариантные $T$ с
   $Tx_1=\dots=Tx_N=0$ образуют пространство размерности $\sum_cm_c(m_c-\operatorname{rank}S_c)$; для гауссовских $x_i$ это
   почти наверное $\sum_cm_c\max(0,m_c-Nk_c)$. Для представления на рёбрах две выборки определяют эквивариантный
   слой, тогда как произвольной матрице $16\times16$ нужно шестнадцать.

<details>
<summary>Решения</summary>

1. $P_\sigma J=J=JP_\sigma$, поэтому эти отображения эквивариантны. Они линейно независимы: если
   $1\otimes A+J\otimes B=0$, то внедиагональные блоки дают $B=0$, а затем $A=0$. По следствию 12.2(a) и
   теореме 2.6 $P_\sigma\otimes1_d\cong(\mathrm{triv}\otimes1_d)\oplus(\mathrm{std}\otimes1_d)$, так что оба неприводимых
   представления имеют кратность $d$, и теорема 12.4 даёт $\dim\operatorname{Hom}_G=2d_1d_2$ — размерность
   пространства отображений вида $1\otimes A+J\otimes B$.

```python
G_in = torch.stack([torch.kron(P, torch.eye(2, dtype=dt)) for P in G_node])     # d_1 = 2
G_out = torch.stack([torch.kron(P, torch.eye(3, dtype=dt)) for P in G_node])    # d_2 = 3
deepsets = EquivariantMap(G_in, G_out)
assert sum(p.numel() for p in deepsets.parameters()) == 2 * 2 * 3
assert equivariance_error(deepsets, G_in, G_out) < 1e-12
T = deepsets.matrix().detach().reshape(n, 3, n, 2)                # blocks T[i, :, j, :] = delta_ij A + B
A_, B_ = T[0, :, 0] - T[0, :, 1], T[0, :, 1]
assert all((T[i, :, j] - (A_ if i == j else 0) - B_).abs().max() < 1e-12 for i in range(n) for j in range(n))
```

2. Если $T$ эквивариантна, то эквивариантна и $\bar T$, так как $\overline{\pi_i(g)}=\pi_i(g)$. Следовательно,
   $\operatorname{Re}T=(T+\bar T)/2$ и $\operatorname{Im}T=(T-\bar T)/2i$ вещественны и эквивариантны, поэтому отображение
   $T\mapsto\operatorname{Re}T$ сюръективно на вещественные эквивариантные матрицы (на них оно тождественно), и
   $\operatorname{Hom}_G(\pi_1,\pi_2)=H_{\mathbb R}\oplus iH_{\mathbb R}$, где $H_{\mathbb R}$ — вещественные эквивариантные
   матрицы. Тогда вещественный базис $H_{\mathbb R}$ является комплексным базисом $\operatorname{Hom}_G$: он порождает всё пространство, а
   комплексное соотношение $\sum_j(\alpha_j+i\beta_j)T_j=0$ с вещественными $T_j$ распадается на вещественную и мнимую части. В
   лемме 12.1 индикаторы орбит образуют именно такой вещественный базис. (Блоки $y_c$ из теоремы 12.4 комплексны даже для вещественных
   данных: $V$ в общем случае комплексен, а параметризация через $\operatorname{Re}T$ использует вдвое больше вещественных параметров,
   чем нужно.)

3. (a) $(\pi(g)f(x))_{gy}=f(x)_y=f(x_y)=f((\pi(g)x)_{gy})$. (b) Для $R(x_1,x_2)=(-x_2,x_1)$ и $x=(1,-1)$
   $f(Rx)=f(1,1)=(1,1)$, но $Rf(x)=R(1,0)=(0,1)$. (c) $\|\pi(g)x\|=\|x\|$, поэтому
   $\varphi(\|\pi(g)x\|)\pi(g)x=\pi(g)\varphi(\|x\|)x$.

```python
relu = lambda v: v.clamp(min=0)
R = torch.tensor([[0.0, -1.0], [1.0, 0.0]], dtype=torch.float64)
v = torch.tensor([1.0, -1.0], dtype=torch.float64)
assert (relu(R @ v) - R @ relu(v)).abs().max() == 1                             # (b)
xr = torch.randn(16, dtype=torch.float64)
assert all(torch.equal(relu(g.real @ xr), g.real @ relu(xr)) for g in G_edge)   # (a): exact
gate = lambda v: torch.tanh(v.norm()) * v
assert max((gate(g @ x) - g @ gate(x)).abs().max().item() for g in G_edge) < 1e-12   # (c)
```

4. $(1\otimes w_c)(e_a\otimes e_j)=\sum_le_a\otimes(w_c)_{lj}e_l$, поэтому у $V^*Tx_i$ в секторе $c$ матрица коэффициентов
   $X_{i,c}w_c^T$. Следовательно, $Tx_i=0$ для всех $i$ тогда и только тогда, когда $S_cw_c^T=0$ для всех $c$, т. е. когда каждый столбец $w_c^T$ лежит в
   $\ker S_c$ — пространстве размерности $m_c-\operatorname{rank}S_c$: $w_c$ с этим свойством образуют пространство
   размерности $m_c(m_c-\operatorname{rank}S_c)$, а $\Phi$ инъективно. Для гауссовских $x_i$ вектор $V^*x_i$ гауссовский
   ($V$ унитарен), поэтому у $S_c$ независимые гауссовские элементы, и почти наверное её ранг равен $\min(Nk_c,m_c)$: для
   ведущей квадратной подматрицы этого размера $|\det|^2$ — вещественный многочлен от вещественных и мнимых частей
   элементов, не равный тождественно нулю (он равен $1$, когда подматрица единичная), а такой многочлен обращается в нуль
   лишь на множестве лебеговой меры нуль (доказано в разделе 4.3). Для представления на рёбрах $(k,m)=(1,2),(2,1),(3,1),(3,3)$, и
   вклад даёт только сектор $(1,2)$: $2\cdot(2-1)=2$ при $N=1$ и $0$ при $N=2$.

```python
basis = D_edge.basis('commutant')                                # (15, 16, 16), Hilbert-Schmidt orthonormal
for N, lost in ((1, 2), (2, 0)):
    Xs = torch.randn(16, N, dtype=dt)
    M = torch.stack([(B @ Xs).reshape(-1) for B in basis], dim=1)               # w -> (T x_1, ..., T x_N)
    assert torch.linalg.matrix_rank(M).item() == 15 - lost
```

</details>

**Что дальше.** В последней, [13-й главе](13_jordan_algebras.md), обсуждается, почему квантовая механика использует *-алгебры, а не йордановы алгебры. В [главе 4](04_symmetries.md) развивается
используемая здесь теория представлений: лемма Шура, характеры и рандомизированный алгоритм, лежащий в основе `dc.from_group`
(разделы 4.1–4.3). [Глава 9](09_spin_chains.md) применяет те же секторы к спиновым цепочкам, а мера Хаара на
$\prod_cU(m_c)$ из [главы 5](05_random_matrices.md) даёт естественную случайную инициализацию эквивариантного
унитарного слоя. По математике см. J.-P. Serre, *Linear Representations of Finite Groups* (Springer 1977) и
W. Fulton and J. Harris, *Representation Theory* (Springer 1991); по эквивариантному обучению — M. Bronstein, J. Bruna,
T. Cohen and P. Veličković, *Geometric Deep Learning* (2021) и M. Weiler et al., *Equivariant and Coordinate
Independent Convolutional Networks* (2023); по алгоритму разложения — ссылки в
`torch_vn_algebra/decompose.py`.
