[English](../07_tensor_products.md) · [Español](../es/07_tensor_products.md) · [Français](../fr/07_tensor_products.md) · [Deutsch](../de/07_tensor_products.md) · [中文](../zh/07_tensor_products.md) · [日本語](../ja/07_tensor_products.md) · **Русский**

# 7. Тензорные произведения и запутанность

В этой главе изучается тензорное произведение двух алгебр с секторами (физически — алгебра двух систем
вместе; см. блоки *Физический смысл*). При наличии секторов произведений два: обычное, в котором каждая
часть сохраняет свой собственный заряд, и произведение *со слиянием*, в котором сохраняется лишь полный заряд.
Мы доказываем строение обоих, определяем частичный след как сопряжённое к вложению отображение и изучаем запутанность.

**Вы узнаете**, как перемножаются секторы и в каком порядке их располагает библиотека, почему частичный
след является каналом, а после нормировки — условным ожиданием конечного индекса, каково редуцированное
состояние для состояний, сохраняющих заряд, как запутанность распадается на числовую и конфигурационную
части, а также среднюю чистоту и среднюю энтропию случайных состояний.

Всюду далее $\mathcal M_1=\bigoplus_{c=1}^{C_1}M_{k_c}(\mathbb C)$ действует на
$H_1=\bigoplus_c\mathbb C^{k_c}$ с минимальными центральными проекторами $P_c$, а
$\mathcal M_2=\bigoplus_{d=1}^{C_2}M_{l_d}(\mathbb C)$ действует на $H_2=\bigoplus_d\mathbb C^{l_d}$
с минимальными центральными проекторами $Q_d$. $\operatorname{Tr}$ — след пространства, на котором
действует оператор; на $\mathcal M_i$ это след `Tr_blunt` из главы 3. Для $X$ на
$\mathbb C^k\otimes\mathbb C^l$ пишем $X_{(a,b),(a',b')}=\langle e_a\otimes e_b,X\,e_{a'}\otimes e_{b'}\rangle$
и отождествляем $\mathbb C^k\otimes\mathbb C^l$ с $\mathbb C^{kl}$ по правилу $e_a\otimes e_b\mapsto e_{(a-1)l+b}$
(соглашение Кронекера), так что $(A\otimes B)_{(a,b),(a',b')}=A_{aa'}B_{bb'}$.

## 7.1 Тензорное произведение двух алгебр

Для $*$-алгебр $\mathcal M_1$ на $H_1$ и $\mathcal M_2$ на $H_2$ **тензорным произведением**
$\mathcal M_1\otimes\mathcal M_2$ называется линейная оболочка операторов $A\otimes B$ на
$H_1\otimes H_2$, $A\in\mathcal M_1$, $B\in\mathcal M_2$. Поскольку
$(A\otimes B)(A'\otimes B')=AA'\otimes BB'$, $(A\otimes B)^*=A^*\otimes B^*$ и $1\otimes1=1$, это
$*$-алгебра, а значит, алгебра фон Неймана (теорема 2.2).

**Теорема 7.1.** (a) $H_1\otimes H_2=\bigoplus_{c,d}\mathbb C^{k_c}\otimes\mathbb C^{l_d}$ — ортогональная
сумма образов проекторов $P_c\otimes Q_d$, и
$$
\mathcal M_1\otimes\mathcal M_2=\bigoplus_{c,d}M_{k_c}\otimes M_{l_d}\cong\bigoplus_{c,d}M_{k_cl_d}(\mathbb C):
$$
она состоит из всех операторов $X=\sum_{c,d}(P_c\otimes Q_d)X(P_c\otimes Q_d)$ с произвольными блоками
на $\mathbb C^{k_c}\otimes\mathbb C^{l_d}$. (b) Её минимальные центральные проекторы — это $C_1C_2$
проекторов $P_c\otimes Q_d$, и $\dim(\mathcal M_1\otimes\mathcal M_2)=\dim\mathcal M_1\dim\mathcal M_2$.
(c) Для следов $\operatorname{tr}_w$ на $\mathcal M_1$ и $\operatorname{tr}_v$ на $\mathcal M_2$ (теорема 3.1)
$\operatorname{tr}_w(A)\operatorname{tr}_v(B)=\operatorname{tr}_{w\otimes v}(A\otimes B)$, где
$(w\otimes v)_{cd}=w_cv_d$. В частности, `Tr_blunt`, `Tr_norm` и `tau_vN` произведения равны
произведениям одноимённых следов сомножителей.

*Доказательство.* (a) Проекторы $P_c\otimes Q_d$ попарно ортогональны, их сумма равна $1\otimes1$, и
образ $P_c\otimes Q_d$ есть $\mathbb C^{k_c}\otimes\mathbb C^{l_d}$. Так как $P_c$, $Q_d$ центральны,
$(P_c\otimes Q_d)(A\otimes B)=AP_c\otimes BQ_d=(A\otimes B)(P_c\otimes Q_d)$, поэтому каждый элемент
$\mathcal M_1\otimes\mathcal M_2$ блочно-диагонален. Обратно, для матричных единиц
$E^c_{aa'}\in\mathcal M_1$ и $F^d_{bb'}\in\mathcal M_2$ блоков $c$ и $d$ имеем
$E^c_{aa'}\otimes F^d_{bb'}=|e_a\otimes e_b\rangle\langle e_{a'}\otimes e_{b'}|$ внутри
$\mathbb C^{k_c}\otimes\mathbb C^{l_d}$; это матричные единицы $M_{k_cl_d}(\mathbb C)$, так что всякий
оператор с носителем в блоке $(c,d)$ лежит в $\mathcal M_1\otimes\mathcal M_2$. (b) Центр
прямой суммы полных матричных алгебр порождается проекторами на слагаемые (глава 2,
упражнение 1), и $\sum_{c,d}k_c^2l_d^2=\sum_ck_c^2\sum_dl_d^2$. (c)
$\operatorname{Tr}(X\otimes Y)=\sum_{a,b}X_{aa}Y_{bb}=\operatorname{Tr}X\operatorname{Tr}Y$, поэтому
$\operatorname{tr}_{w\otimes v}(A\otimes B)=\sum_{c,d}w_cv_d\operatorname{Tr}A_c\operatorname{Tr}B_d=\operatorname{tr}_w(A)\operatorname{tr}_v(B)$.
Веса $1$, $1/k_c$, $1/(C_1k_c)$ трёх следов (раздел 3.1) перемножаются в $1$,
$1/(k_cl_d)$, $1/(C_1C_2k_cl_d)$ — веса тех же следов на произведении. $\square$

`tensor_product(alg1, alg2)` строит эту алгебру. Сектор $(c,d)$ имеет номер $s=cC_2+d$ (в коде
нумерация с $0$) и заряд $(q_c,r_d)$, а внутри него $e_a\otimes e_b$ имеет номер $a\,l_d+b$.
`kron(A, B, alg12)` вычисляет $A\otimes B$. Порядок секторов отличается от кронекеровского порядка в
$H_1\otimes H_2$ перестановкой базиса, которую явно строит приведённая ниже функция `kron_rows`; с её помощью
мы сравниваем результаты с кронекеровскими произведениями полных матриц.

```python
import math
import numpy as np
import torch
from torch_vn_algebra import (TypeIAlgebra, DensityMatrix, states, composite,
                              tensor_product, fused_tensor_product, kron, partial_trace)

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def full(X):
    """An operator of (+)_c M_{k_c} as a block-diagonal matrix on (+)_c C^{k_c}, shape (batch, D, D)."""
    return torch.stack([torch.block_diag(*[x[c, :k, :k] for c, k in enumerate(X.algebra.k_factors)])
                        for x in X.matrix])

def kron_rows(alg12):
    """For every sector of a product algebra, the positions of its basis vectors in H1 (x) H2 (Kronecker order)."""
    k1, k2 = alg12.tensor_factors[0].k_factors, alg12.tensor_factors[1].k_factors
    return [torch.tensor([(sum(k1[:c]) + a) * sum(k2) + sum(k2[:d]) + b
                          for c, d, off in pairs for a in range(k1[c]) for b in range(k2[d])])
            for pairs in alg12.product_blocks]              # pair blocks (c, d, offset), by offset

def full12(X):
    """An operator of a (plain or fused) product algebra as a matrix on H1 (x) H2 in Kronecker order."""
    f1, f2 = X.algebra.tensor_factors
    N = sum(f1.k_factors) * sum(f2.k_factors)
    out = torch.zeros(X.matrix.shape[0], N, N, dtype=X.matrix.dtype)
    for s, r in enumerate(kron_rows(X.algebra)):
        out[:, r[:, None], r[None, :]] = X.matrix[:, s, :len(r), :len(r)]
    return out

def random_operator(alg):
    return alg.from_blocks([torch.randn(k, k, dtype=dt) for k in alg.k_factors])

alg1 = TypeIAlgebra([1, 2], [1, 2], precision='double')            # C (+) M_2
alg2 = TypeIAlgebra([2, 3], [2, 3], precision='double')            # M_2 (+) M_3
alg12 = tensor_product(alg1, alg2)
print(alg12.k_factors, alg12.charges, sum(k * k for k in alg12.k_factors))     # 4 sectors, dimension 65
assert sum(k * k for k in alg12.k_factors) == 5 * 13                           # Theorem 7.1(b)

A, A2, B, B2 = random_operator(alg1), random_operator(alg1), random_operator(alg2), random_operator(alg2)
AB = kron(A, B, alg12)
assert close(full12(AB)[0], torch.kron(full(A)[0], full(B)[0]))               # Theorem 7.1(a)
assert close(AB @ kron(A2, B2, alg12), kron(A @ A2, B @ B2, alg12))
for tr in ('Tr_blunt', 'Tr_norm', 'tau_vN'):                                   # Theorem 7.1(c)
    assert (getattr(AB, tr)() - getattr(A, tr)() * getattr(B, tr)()).abs().max() < 1e-10
```

> **Физический смысл** (физическая интерпретация, а не теорема этого курса). Если каждая из двух систем подчиняется своему правилу суперотбора — две коробки,
> которые не могут обмениваться частицами, — то у составной системы есть сектор для каждой *пары* зарядов,
> и никакая наблюдаемая не связывает разные пары. Составной сектор входит с весом, равным произведению весов
> (теорема 7.1(c)).

## 7.2 Частичный след

**Определение.** Пусть $\rho\in\mathcal M_1\otimes\mathcal M_2$. **Частичным следом**
$\operatorname{Tr}_2\rho\in\mathcal M_1$ называется элемент, для которого
$$
\operatorname{Tr}\big((\operatorname{Tr}_2\rho)\,A\big)=\operatorname{Tr}\big(\rho\,(A\otimes1)\big)\qquad\text{for all }A\in\mathcal M_1 .
$$
Он существует и единствен по теореме 3.2, применённой к линейному функционалу
$A\mapsto\operatorname{Tr}(\rho(A\otimes1))$ и точному следу $\operatorname{Tr}$ алгебры $\mathcal M_1$:
частичный след сопряжён (относительно следов `Tr_blunt`) к вложению $A\mapsto A\otimes1$.
$\operatorname{Tr}_1$ определяется так же; для состояния $\rho$ элементы $\operatorname{Tr}_2\rho$ и
$\operatorname{Tr}_1\rho$ называются **редуцированными состояниями**. Для полных алгебр $\mathcal M_i=B(H_i)$ это
**частичный след гильбертова пространства** $\operatorname{Tr}_{H_2}:B(H_1\otimes H_2)\to B(H_1)$.

**Предложение 7.2.** (a) Для $\rho=\bigoplus_{c,d}\rho_{cd}$
$$
(\operatorname{Tr}_2\rho)_c=\sum_d\operatorname{tr}_{l_d}\rho_{cd},\qquad
(\operatorname{tr}_lX)_{aa'}=\sum_{b=1}^lX_{(a,b),(a',b)}\quad(X\text{ on }\mathbb C^k\otimes\mathbb C^l),
$$
и $\operatorname{Tr}_2\rho=\operatorname{Tr}_{H_2}\rho$, где $(\operatorname{Tr}_{H_2}X)_{ii'}=\sum_jX_{(i,j),(i',j)}$.
(b) $\operatorname{Tr}_2(A\otimes B)=(\operatorname{Tr}B)\,A$.
(c) $\operatorname{Tr}_2\big((A\otimes1)X(A'\otimes1)\big)=A\,(\operatorname{Tr}_2X)\,A'$ для $A,A'\in\mathcal M_1$.
(d) $\operatorname{Tr}(\operatorname{Tr}_2X)=\operatorname{Tr}X$, и $\operatorname{Tr}_2$ переводит состояния в состояния.

*Доказательство.* (a) При $A=\bigoplus_cA_c$ имеем $A\otimes1=\bigoplus_{c,d}A_c\otimes1_{l_d}$ и
$$
\operatorname{Tr}\big(\rho(A\otimes1)\big)=\sum_{c,d}\sum_{a,a',b}(\rho_{cd})_{(a,b),(a',b)}(A_c)_{a'a}
=\sum_c\operatorname{Tr}\Big(\Big(\sum_d\operatorname{tr}_{l_d}\rho_{cd}\Big)A_c\Big),
$$
так что элемент справа удовлетворяет определению и равен $\operatorname{Tr}_2\rho$. То же
вычисление для $B(H_i)$ даёт формулу для $\operatorname{Tr}_{H_2}$; для блочно-диагональной $\rho$
элементы $\rho_{(i,j),(i',j)}$ с $i,i'$ из разных блоков равны нулю, так что $\operatorname{Tr}_{H_2}\rho$ имеет
блоки из (a). (b) $\operatorname{tr}_l(X\otimes Y)=(\operatorname{Tr}Y)X$. (c) Для $B\in\mathcal M_1$
$\operatorname{Tr}\big(\operatorname{Tr}_2((A\otimes1)X(A'\otimes1))B\big)=\operatorname{Tr}\big(X(A'BA\otimes1)\big)
=\operatorname{Tr}\big(A(\operatorname{Tr}_2X)A'B\big)$; остаётся воспользоваться единственностью. (d) Возьмём $A=1$. Если $\rho$ —
состояние, то $A\mapsto\operatorname{Tr}(\rho(A\otimes1))$ положительно, так как $A^*A\otimes1=(A\otimes1)^*(A\otimes1)$, и
равно $1$ при $A=1$, поэтому его плотность $\operatorname{Tr}_2\rho$ — состояние (теорема 3.2). $\square$

Каналам посвящена глава 6; здесь нам нужно лишь следующее определение: линейное отображение
$\Phi$ **вполне положительно**, если $\mathrm{id}_n\otimes\Phi$ положительно на
$M_n(\mathbb C)\otimes(\cdot)$ при каждом $n$, и является **каналом**, если оно к тому же сохраняет след.

**Теорема 7.3.** Пусть $j$ пробегает стандартный базис $H_2$ и $K_j=1\otimes\langle j|:H_1\otimes H_2\to H_1$,
$K_j(x\otimes y)=\langle j,y\rangle x$. Тогда
$$
\operatorname{Tr}_{H_2}X=\sum_jK_jXK_j^*,\qquad\sum_jK_j^*K_j=1,\qquad\sum_jK_j^*AK_j=A\otimes1 .
$$
Следовательно, $\operatorname{Tr}_2$ — канал из $\mathcal M_1\otimes\mathcal M_2$ в $\mathcal M_1$, и
сопряжённое к нему относительно следов `Tr_blunt` отображение есть $A\mapsto A\otimes1$.

*Доказательство.* $K_j^*x=x\otimes e_j$, поэтому $\langle e_a,K_jXK_j^*e_{a'}\rangle=X_{(a,j),(a',j)}$; суммирование по $j$
даёт предложение 7.2(a). Далее, $\sum_jK_j^*K_j(x\otimes y)=\sum_j\langle j,y\rangle x\otimes e_j=x\otimes y$
и $\sum_jK_j^*AK_j=\sum_jA\otimes|j\rangle\langle j|=A\otimes1$. Для $Y\ge0$
$(\mathrm{id}_n\otimes\operatorname{Tr}_{H_2})(Y)=\sum_j(1_n\otimes K_j)Y(1_n\otimes K_j)^*\ge0$, поскольку
$\langle x,TYT^*x\rangle=\langle T^*x,YT^*x\rangle\ge0$; это сохраняется и на подалгебре
$\mathcal M_1\otimes\mathcal M_2$. Наконец, $\operatorname{Tr}\sum_jK_jXK_j^*A=\operatorname{Tr}X\sum_jK_j^*AK_j$,
что при $A=1$ означает сохранение следа, а в общем случае — двойственность. $\square$

`partial_trace(X, keep=1)` переводит каждый блок $\rho_{cd}$ к форме $(k_c,l_d,k_c,l_d)$ и берёт след;
`composite.partial_trace_channel(alg12, keep=1)` — канал из теоремы 7.3 с операторами Крауса
$1_{k_c}\otimes\langle j|$ из сектора $(c,d)$ в сектор $c$. Проверим оба по определению,
независимым вычислением через `numpy.einsum` и по двойственности:

```python
rho = states.random_density_matrix(alg12, batch_size=3)
rho1, rho2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
assert isinstance(rho1, DensityMatrix)                                      # states give states
R = full12(rho).numpy().reshape(3, 3, 5, 3, 5)                              # indices (batch, i, j, i', j')
assert np.abs(full(rho1).numpy() - np.einsum('xijkj->xik', R)).max() < 1e-12
assert np.abs(full(rho2).numpy() - np.einsum('xijil->xjl', R)).max() < 1e-12
one2 = alg2.identity()
assert (states.expectation(rho, kron(A, one2, alg12)) - states.expectation(rho1, A)).abs().max() < 1e-12
assert close(partial_trace(AB, keep=1), A * B.trace.item())                 # Proposition 7.2(b)

Ptr = composite.partial_trace_channel(alg12, keep=1)
print(Ptr)
assert Ptr.is_trace_preserving()
assert close(Ptr(rho), rho1)
assert close(Ptr.adjoint('blunt')(A), kron(A, one2, alg12))                 # Theorem 7.3: the dual
```

> **Физический смысл** (физическая интерпретация, а не теорема этого курса). Редуцированное состояние $\operatorname{Tr}_2\rho$ — это всё, что узнаёт
> наблюдатель, способный измерять лишь первую систему, — по определению оно воспроизводит каждое среднее значение
> наблюдаемой $A\otimes1$. Составные состояния с одинаковыми редуцированными состояниями различаются корреляциями.

**Частичный след как условное ожидание.** Пусть $\mathcal M_2=M_m(\mathbb C)$, так что
секторы $\mathcal M_1\otimes M_m$ совпадают с секторами $\mathcal M_1$ (теорема 7.1); пусть
$\iota(A)=A\otimes1_m$ и $\mathcal E=\frac1m\operatorname{Tr}_2:\mathcal M_1\otimes M_m\to\mathcal M_1$.
Пусть $e_{ab}$ — матричные единицы $M_m$.

**Предложение 7.4.** (a) $\iota\circ\mathcal E$ — условное ожидание из $\mathcal M_1\otimes M_m$
на $\mathcal M_1\otimes1$ (раздел 2.6), и для всякого точного следа $\operatorname{tr}_w$ алгебры
$\mathcal M_1\otimes M_m$ оно единственное, сохраняющее $\operatorname{tr}_w$. (b) $m^2$ элементов
$u_{ab}=\sqrt m\,(1\otimes e_{ab})$ удовлетворяют равенству $X=\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*X))$ при всех $X$,
и $\sum_{a,b}u_{ab}u_{ab}^*=m^2\,1$. (c) Для всякого семейства $(v_k)$ с
$X=\sum_kv_k\,\iota(\mathcal E(v_k^*X))$ при всех $X$ выполнено $\sum_kv_kv_k^*=m^2\,1$.

*Доказательство.* (a) $\mathcal E$ положительно (теорема 7.3), $\iota$ тоже
($T^*T\otimes1=(T\otimes1)^*(T\otimes1)$); $\mathcal E(A\otimes1)=A$ по предложению 7.2(b), а
свойство бимодульности — это предложение 7.2(c). Будем обозначать через $\operatorname{tr}_w$ и след с теми же
весами $w_c$ на $\mathcal M_1$; тогда $\operatorname{tr}_w(Y\otimes1)=m\operatorname{tr}_w(Y)$. Определение
$\operatorname{Tr}_2$, применённое к $AP_c$, умноженное на $w_c$ и просуммированное по $c$, даёт
$\operatorname{tr}_w(X(A\otimes1))=\operatorname{tr}_w((\operatorname{Tr}_2X)A)$, а при $A=1$ получаем
$\operatorname{tr}_w(\iota\mathcal E(X))=\operatorname{tr}_w(X)$. Если $F=\iota\circ F_0$ — другое условное
ожидание на $\mathcal M_1\otimes1$, сохраняющее $\operatorname{tr}_w$, то, поскольку $F(X(A\otimes1))=F(X)(A\otimes1)$,
$\operatorname{tr}_w(X(A\otimes1))=\operatorname{tr}_w(F(X)(A\otimes1))=m\operatorname{tr}_w(F_0(X)A)$, поэтому
$\operatorname{tr}_w((\operatorname{Tr}_2X-mF_0(X))A)=0$ для всех $A\in\mathcal M_1$, и $F_0=\mathcal E$ в силу точности
(возьмём $A=(\operatorname{Tr}_2X-mF_0(X))^*$). (b) Запишем $X=\sum_{a,b}X_{ab}\otimes e_{ab}$ с
$X_{ab}\in\mathcal M_1$ (теорема 7.1(a)). По предложению 7.2(b)
$\mathcal E(u_{ab}^*X)=\frac{\sqrt m}m\sum_{b'}\operatorname{Tr}_2(X_{ab'}\otimes e_{bb'})=X_{ab}/\sqrt m$, поэтому
$u_{ab}\,\iota(\mathcal E(u_{ab}^*X))=X_{ab}\otimes e_{ab}$; кроме того,
$\sum_{a,b}u_{ab}u_{ab}^*=m\sum_{a,b}1\otimes e_{aa}=m^21$. (c) $\operatorname{Tr}_2(Y^*)=(\operatorname{Tr}_2Y)^*$ по
предложению 7.2(a), так что сопряжение разложения для $Y=u_{ab}$ даёт
$u_{ab}^*=\sum_k\iota(\mathcal E(u_{ab}^*v_k))\,v_k^*$, и по (b)
$\sum_{a,b}u_{ab}u_{ab}^*=\sum_k\big(\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*v_k))\big)v_k^*=\sum_kv_kv_k^*$. $\square$

Семейство, как в (c), называется *квазибазисом* для $\mathcal E$, а $\sum_kv_kv_k^*$, которое по (c) от
него не зависит и по (b) существует, — **индексом** $\operatorname{Ind}\mathcal E=m^2$ отображения $\mathcal E$ (это
определение; название следует Y. Watatani, Mem. Amer. Math. Soc. 83, no. 424, 1990). **Матрица включения**
включения $\bigoplus_cM_{k_c}$ в $\bigoplus_sM_{n_s}$ имеет элементы $\Lambda_{sc}$ — числа копий блока $c$ подалгебры
в блоке $s$ большей алгебры (глава 9 использует её для башен алгебр). Для фактора $M_k\otimes1\subset M_k\otimes M_m$:
если упорядочить базис $e_a\otimes e_b$ пространства $\mathbb C^{km}$ сначала по $b$, то $A\otimes1_m=\operatorname{diag}(A,\dots,A)$
с $m$ копиями, так что $\Lambda=(m)$ и $m^2=\|\Lambda\|^2=\dim(M_k\otimes M_m)/\dim M_k$. (Дополнительная
литература, здесь не используемая: V. F. R. Jones, Invent. Math. 72, 1, 1983, об индексе субфакторов, и
F. M. Goodman, P. de la Harpe and V. F. R. Jones, *Coxeter graphs and towers of algebras*, Springer,
1989, о матрицах включения.) Для общей $\mathcal M_1$ включение есть прямая сумма
$C_1$ включений $M_{k_c}\otimes1\subset M_{k_c}\otimes M_m$, $\Lambda=m1_{C_1}$. Модуль
`torch_vn_algebra.inclusions` реализует включения, их условные ожидания и индексы;
`from_tensor_product(alg1, m)` строит данное включение:

```python
from torch_vn_algebra import inclusions
m = 3
inc = inclusions.from_tensor_product(alg1, m)                # alg1 (x) 1 c alg1 (x) M_3
X = random_operator(inc.M)
for trace in ('Tr_blunt', 'Tr_norm', 'tau_vN'):
    E = inc.conditional_expectation(trace.replace('Tr_', ''))
    assert close(E(X), partial_trace(X, keep=1) * (1 / m))   # Proposition 7.4(a): the same E for every trace
    assert abs(getattr(inc.embed(E(X)), trace)() - getattr(X, trace)()).max() < 1e-10
    assert close(E.index(), inc.M.identity() * m**2)         # Proposition 7.4(c)

def unit(a, b):
    e = torch.zeros(m, m, dtype=dt)
    e[a, b] = 1
    return math.sqrt(m) * kron(alg1.identity(), inc.M.tensor_factors[1].from_blocks([e]), inc.M)

u = [unit(a, b) for a in range(m) for b in range(m)]         # the quasi-basis of Proposition 7.4(b)
assert close(sum((v @ inc.embed(E(v.adjoint() @ X))).matrix for v in u), X)
assert close(sum((v @ v.adjoint()).matrix for v in u), inc.M.identity() * m**2)
print(inc.matrix.tolist(), inc.component_indices())           # Lambda = 3 * 1_2, index 9 per component
```

## 7.3 Слияние зарядов

Пусть теперь сектор $c$ алгебры $\mathcal M_1$ несёт заряд $q_c$, сектор $d$ алгебры $\mathcal M_2$ —
заряд $r_d$, а **правило слияния** $f$ задаёт полный заряд $f(q_c,r_d)$ — сложение для
чисел частиц, сложение по модулю $2$ для фермионной чётности или любое другое отображение. Положим
$$
H_q=\bigoplus_{f(q_c,r_d)=q}\mathbb C^{k_c}\otimes\mathbb C^{l_d},\qquad
\Pi_q=\sum_{f(q_c,r_d)=q}P_c\otimes Q_d,\qquad n_q=\sum_{f(q_c,r_d)=q}k_cl_d .
$$
**Произведение со слиянием** $\mathcal M_1\otimes_f\mathcal M_2=\{\Pi_q\}'$ состоит из операторов на
$H_1\otimes H_2$, сохраняющих полный заряд.

**Предложение 7.5.** (a) $\mathcal M_1\otimes_f\mathcal M_2=\bigoplus_qB(H_q)\cong\bigoplus_qM_{n_q}(\mathbb C)$
с минимальными центральными проекторами $\Pi_q$. (b) $\mathcal M_1\otimes\mathcal M_2\subseteq\mathcal M_1\otimes_f\mathcal M_2$,
причём равенство имеет место тогда и только тогда, когда каждое $H_q$ содержит лишь одну пару $(c,d)$. (c) Если заряды —
вещественные числа, а $f$ — сложение, то $\mathcal M_1\otimes_f\mathcal M_2=\{N\}'$, где
$N=N_1\otimes1+1\otimes N_2$, $N_1=\sum_cq_cP_c$, $N_2=\sum_dr_dQ_d$.

*Доказательство.* (a) Каждая пара лежит ровно в одном $H_q$, поэтому $\Pi_q$ — попарно ортогональные
проекторы с суммой $1$, и $X$ коммутирует со всеми ними тогда и только тогда, когда $X=\sum_q\Pi_qX\Pi_q$
с произвольными блоками в $B(H_q)$; центр — как в теореме 7.1(b). (b) Каждый элемент
$\mathcal M_1\otimes\mathcal M_2$ коммутирует с каждым $P_c\otimes Q_d$, а значит, и с их суммами $\Pi_q$.
Если каждое $H_q$ — это одно $\mathbb C^{k_c}\otimes\mathbb C^{l_d}$, то обе алгебры равны
$\bigoplus_{c,d}B(\mathbb C^{k_c}\otimes\mathbb C^{l_d})$; если $H_q$ содержит две пары, то $B(H_q)$ содержит
операторы, переводящие блок одной пары в блок другой, а они не блочно-диагональны. (c)
$N=\sum_{c,d}(q_c+r_d)P_c\otimes Q_d=\sum_qq\,\Pi_q$, так что ненулевые $\Pi_q$ — спектральные
проекторы $N$. Они являются многочленами от $N$ (доказательство теоремы 1.5(a)), а $N$ — их линейная комбинация,
поэтому $X$ коммутирует с $N$ тогда и только тогда, когда он коммутирует с каждым $\Pi_q$. $\square$

`fused_tensor_product(alg1, alg2, fuse)` (по умолчанию `fuse` — сложение) упорядочивает секторы по
полному заряду, а пары $(c,d)$ внутри сектора — лексикографически, каждую с базисом $e_a\otimes e_b$.
Пример: бозоны с жёсткой сердцевиной (hard-core bosons). На $L$ узлах, $(\mathbb C^2)^{\otimes L}$ с числом заполнения
$N=\sum_in^{(i)}$, $n=\operatorname{diag}(0,1)$, собственное подпространство $N=q$ натянуто на $\binom Lq$ векторов
произведения базиса с $q$ занятыми узлами, поэтому операторы $\{N\}'$, сохраняющие число частиц, образуют
$\bigoplus_qM_{\binom Lq}$ (предложение 4.10(2)): для двух узлов $M_1\oplus M_2\oplus M_1$ (заряды $0,1,2$), для
трёх узлов $M_1\oplus M_3\oplus M_3\oplus M_1$. По предложению 7.5(c) произведение со слиянием есть $\{N\}'$ для
полного числа заполнения пяти узлов с точностью до порядка базиса, и действительно
$n_q=\sum_c\binom2c\binom3{q-c}=\binom5q$ (тождество Вандермонда: $q$-элементное подмножество пяти узлов
состоит из $c$ первых двух и $q-c$ последних трёх).

```python
bos1 = TypeIAlgebra([1, 2, 1], [1, 2, 1], precision='double', charges=[0, 1, 2])
bos2 = TypeIAlgebra([1, 3, 3, 1], [1, 3, 3, 1], precision='double', charges=[0, 1, 2, 3])
F = fused_tensor_product(bos1, bos2)                         # U(1): total particle number
P = tensor_product(bos1, bos2)                               # particle numbers conserved separately
print(F.charges, F.k_factors, sum(k * k for k in F.k_factors), sum(k * k for k in P.k_factors))
assert F.k_factors == [math.comb(5, q) for q in range(6)]
assert sum(k * k for k in F.k_factors) == 252 and sum(k * k for k in P.k_factors) == 120

N1, N2 = full(bos1.central([0.0, 1.0, 2.0]))[0], full(bos2.central([0.0, 1.0, 2.0, 3.0]))[0]
N = torch.kron(N1, torch.eye(8, dtype=dt)) + torch.kron(torch.eye(4, dtype=dt), N2)
X = full12(F.operator(torch.randn(1, F.C, F.k_max, F.k_max, dtype=dt)))[0]
assert close(X @ N, N @ X)                                   # Proposition 7.5(c)
A, B = random_operator(bos1), random_operator(bos2)
assert close(full12(kron(A, B, F)), full12(kron(A, B, P)))   # Proposition 7.5(b)
```

Алгебра со слиянием имеет размерность $\sum_q\binom5q^2=\binom{10}5=252$ (снова тождество Вандермонда), обычное
произведение — $(1+4+1)(1+9+9+1)=120$ (теорема 7.1(b)). По теореме 7.1(a) и предложению 7.5(b) обычное произведение
состоит из тех элементов произведения со слиянием, которые блочно-диагональны относительно пар $(c,d)$; остальные
$132$ измерения натянуты на операторы, переводящие блок пары $(c,d)$ в блок другой пары $(c',d')$ того же $H_q$, $c'\ne c$
(а значит, $d'\ne d$): они меняют заряд каждой части, сохраняя полный, т. е. здесь перемещают частицы между двумя частями.

Состояние на $\mathcal M_1\otimes_f\mathcal M_2$ может иметь когерентности между блоками разных пар, поэтому
предложение 7.2 неприменимо. Пусть $E_1(X)=\sum_cP_cXP_c$ на $B(H_1)$ (условное ожидание из раздела
2.6); оно отображает на $\mathcal M_1$, и $\operatorname{Tr}(E_1(X)A)=\operatorname{Tr}(XA)$ для
$A\in\mathcal M_1$, поскольку $A=\sum_cP_cAP_c$.

**Теорема 7.6.** Пусть $\rho$ — плотность состояния $\omega$ на $\mathcal M_1\otimes_f\mathcal M_2$.
Плотность его ограничения $A\mapsto\omega(A\otimes1)$ на $\mathcal M_1$ равна
$$
\rho_1=E_1\big(\operatorname{Tr}_{H_2}\rho\big).
$$
Блок $(c,c')$ матрицы $\operatorname{Tr}_{H_2}\rho$ равен $\sum_d\operatorname{tr}_{l_d}\rho_{(c,d),(c',d)}$, где сумма берётся по
тем $d$, для которых $f(q_c,r_d)=f(q_{c'},r_d)$, а $\rho_{(c,d),(c',d)}$ — блок $\rho$ из пары $(c',d)$ в пару $(c,d)$.
Следовательно, $E_1$ действует тривиально на все $\operatorname{Tr}_{H_2}\rho$ тогда и только тогда, когда
не существует $c\neq c'$ и $d$ с $f(q_c,r_d)=f(q_{c'},r_d)$; например, если $f$ — сложение и
все $q_c$ попарно различны.

*Доказательство.* $A\otimes1\in\mathcal M_1\otimes_f\mathcal M_2$ (предложение 7.5(b)) и
$\omega(A\otimes1)=\operatorname{Tr}(\rho(A\otimes1))=\operatorname{Tr}\big((\operatorname{Tr}_{H_2}\rho)A\big)
=\operatorname{Tr}\big(E_1(\operatorname{Tr}_{H_2}\rho)A\big)$. Так как $E_1(\operatorname{Tr}_{H_2}\rho)\in\mathcal M_1$, это и
есть плотность в силу единственности в теореме 3.2. Формула для блоков — это предложение 7.2(a) для
$\operatorname{Tr}_{H_2}$ с $i$ в блоке $c$, $i'$ в блоке $c'$ и $j$ в блоке $d$; блок
$\rho$ между $(c',d)$ и $(c,d)$ равен нулю, если только обе пары не лежат в одном $H_q$. Если $c\neq c'$
и $d$ с $f(q_c,r_d)=f(q_{c'},r_d)=q$ существуют, возьмём единичные векторы $x\in\mathbb C^{k_c}$,
$x'\in\mathbb C^{k_{c'}}$, $y\in\mathbb C^{l_d}$ и $\psi=(x+x')\otimes y/\sqrt2\in H_q$. Тогда
$|\psi\rangle\langle\psi|\in B(H_q)$, а $\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|=\frac12|x+x'\rangle\langle x+x'|$
имеет ненулевой блок $\frac12|x\rangle\langle x'|$, который $E_1$ удаляет. $\square$

По блочной формуле теоремы 7.6 $(\rho_1)_c=\sum_d\operatorname{tr}_{l_d}\rho_{(c,d),(c,d)}$. Обе реализации в библиотеке
построены так, чтобы вычислять это отображение при любом правиле слияния: `partial_trace` берёт след
только диагональных блоков пар $\rho_{(c,d),(c,d)}$, а `partial_trace_channel` использует отдельные индексы
Крауса для пар одного сектора с одинаковым $c$, так что сопряжённое к нему отображение снова $A\mapsto A\otimes1$. Код
проверяет обе, для двух правил слияния ниже и обоих выборов сохраняемого сомножителя, сравнивая с
$E_1(\operatorname{Tr}_{H_2}\rho)$, вычисленным независимо с помощью `einsum`, а сопряжённое — с $A\otimes1$. При
сложении $E_1$ в нашем примере ничего не делает; при правиле чётности $f(q,r)=q+r\bmod2$ заряды $0$
и $2$ алгебры `bos1` сливаются одинаково с каждым $r$, и $E_1$ удаляет настоящие когерентности:

```python
G = fused_tensor_product(bos1, bos2, fuse=lambda q, r: (q + r) % 2)       # only the parity is conserved
print(G.k_factors)

def block_mask(alg):
    """Ones in the diagonal blocks of (+)_c C^{k_c}: X * block_mask(alg) = sum_c P_c X P_c."""
    return full(alg.operator(torch.ones(1, alg.C, alg.k_max, alg.k_max, dtype=dt)))

for alg, coherent in ((F, False), (G, True)):
    rho = states.random_density_matrix(alg, batch_size=2)
    R = full12(rho).reshape(2, 4, 8, 4, 8)
    for keep, X_H, alg_k in ((1, torch.einsum('xijkj->xik', R), bos1), (2, torch.einsum('xijil->xjl', R), bos2)):
        E_X = X_H * block_mask(alg_k)                                     # E(Tr_H rho), Theorem 7.6
        assert close(full(partial_trace(rho, keep=keep)), E_X)
        assert close(full(composite.partial_trace_channel(alg, keep=keep)(rho)), E_X)
        lost = (X_H - E_X).abs().max().item()                             # coherences removed by E
        print(f"keep={keep}: largest coherence removed: {lost:.3e}")
        assert (lost > 1e-3) if coherent else (lost == 0)
    assert close(composite.partial_trace_channel(alg, keep=1).adjoint('blunt')(A), kron(A, bos2.identity(), alg))
```

> **Физический смысл** (физическая интерпретация, а не теорема этого курса). В БКШ-описании сверхпроводника сохраняется лишь чётность числа электронов,
> поэтому глобальное состояние может быть суперпозицией чисел частиц $0$ и $2$ в некоторой
> области. Наблюдатель, чьи наблюдаемые сохраняют число частиц в этой области, такой когерентности
> не видит: для него редуцированное состояние — это блочно-диагональная часть (pinching) из теоремы 7.6 (математическое
> содержание этого блока), смесь чисел частиц.

## 7.4 Запутанность чистых состояний

**Лемма 7.7 (разложение Шмидта).** Всякий единичный вектор $\psi\in\mathbb C^a\otimes\mathbb C^b$ можно
записать в виде $\psi=\sum_{i=1}^rs_i\,u_i\otimes v_i$, где $s_i>0$, $\sum_is_i^2=1$, а
$u_i\in\mathbb C^a$, $v_i\in\mathbb C^b$ ортонормированы. Тогда $\operatorname{Tr}_2|\psi\rangle\langle\psi|=\sum_is_i^2|u_i\rangle\langle u_i|$
и $\operatorname{Tr}_1|\psi\rangle\langle\psi|=\sum_is_i^2|v_i\rangle\langle v_i|$: оба редуцированных состояния
имеют одни и те же ненулевые собственные значения $s_i^2$.

*Доказательство.* Пусть $\Psi$ — $a\times b$-матрица $\Psi_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,\psi\rangle$. По
спектральной теореме $\Psi\Psi^*=\sum_{i=1}^a\lambda_i|u_i\rangle\langle u_i|$ с ортонормированным базисом
$u_i$ и $\lambda_i=\|\Psi^*u_i\|^2\ge0$, занумерованными так, что $\lambda_i>0$ ровно при $i\le r$. Положим
$s_i=\sqrt{\lambda_i}$ и $w_i=\Psi^*u_i/s_i$ при $i\le r$; тогда
$\langle w_i,w_j\rangle=\langle u_i,\Psi\Psi^*u_j\rangle/(s_is_j)=\delta_{ij}$, а $\Psi^*u_i=0$ при $i>r$.
Следовательно, $\Psi=\sum_i|u_i\rangle\langle u_i|\Psi=\sum_{i\le r}|u_i\rangle\langle\Psi^*u_i|=\sum_{i\le r}s_i|u_i\rangle\langle w_i|$,
т. е. $\Psi_{\alpha\beta}=\sum_is_i(u_i)_\alpha\overline{(w_i)_\beta}$: это и есть утверждение с $v_i=\bar w_i$, причём
$\sum_is_i^2=\operatorname{Tr}\Psi\Psi^*=\|\psi\|^2=1$. Наконец,
$|\psi\rangle\langle\psi|=\sum_{i,j}s_is_j|u_i\rangle\langle u_j|\otimes|v_i\rangle\langle v_j|$, и
предложение 7.2(b) даёт редуцированные состояния. $\square$

**Энтропия запутанности** вектора $\psi$ равна
$E(\psi)=S(\operatorname{Tr}_2|\psi\rangle\langle\psi|)=S(\operatorname{Tr}_1|\psi\rangle\langle\psi|)=-\sum_is_i^2\log s_i^2$.
Она обращается в нуль тогда и только тогда, когда $r=1$ (каждый член $-s_i^2\log s_i^2$ неотрицателен и равен $0$ лишь при
$s_i^2=1$), т. е. когда $\psi=u\otimes v$ (тогда $\Psi=uv^T$ имеет ранг $1$, а $r=\operatorname{rank}\Psi\Psi^*=\operatorname{rank}\Psi$), и
$E(\psi)\le\log\min(a,b)$, причём равенство имеет место тогда и только тогда, когда меньшее из редуцированных состояний
максимально смешанное (следствие 3.6, применённое к обоим редуцированным состояниям).

При наличии зарядов мы берём единичный вектор $\psi$ с определённым полным зарядом, т. е. лежащий в одном секторе $H_q$
произведения со слиянием (суперпозиция разных $H_q$ даёт на $\mathcal M_1\otimes_f\mathcal M_2$ то же состояние, что и
смесь, глава 3, упражнение 3); запишем $\psi=\sum_{(c,d)}\psi_{cd}$ с
$\psi_{cd}\in\mathbb C^{k_c}\otimes\mathbb C^{l_d}$.

**Предложение 7.8 (запутанность с разрешением по симметрии).** Пусть $f$ — сложение, все $q_c$
попарно различны и все $r_d$ попарно различны, так что для каждого $c$ не более одного $d=d(c)$ удовлетворяет
$q_c+r_d=q$. Положим $\psi_c=\psi_{c\,d(c)}$ (или $0$), $p_c=\|\psi_c\|^2$ и при $p_c>0$
$\hat\psi_c=\psi_c/\sqrt{p_c}$. Тогда $\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|$ лежит в $\mathcal M_1$ (он
коммутирует с $N_1$) и равен $\bigoplus_cp_c\operatorname{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|$;
редуцированное состояние на $\mathcal M_2$ имеет вероятность сектора $p_c$ в $d(c)$; и
$$
E(\psi)=H(p)+\sum_cp_c\,E(\hat\psi_c),\qquad H(p)=-\sum_cp_c\log p_c .
$$

*Доказательство.* По формуле для блоков из теоремы 7.6 блок $(c,c')$ матрицы $\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|$
есть сумма $\operatorname{tr}_{l_d}|\psi_{cd}\rangle\langle\psi_{c'd}|$ по тем $d$, для которых
$q_c+r_d=q=q_{c'}+r_d$, а это влечёт $c=c'$ и $d=d(c)$. Значит, редуцированное состояние равно
$\bigoplus_c\operatorname{Tr}_2|\psi_c\rangle\langle\psi_c|$ со следами блоков $p_c$; поменяв ролями
сомножители, получаем утверждение для $\mathcal M_2$. По лемме 7.7 $E(\psi)$ — энтропия этого
состояния, и теорема 3.3 с $\sigma_c=\operatorname{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|$, $S(\sigma_c)=E(\hat\psi_c)$
даёт формулу. $\square$

$H(p)$ — это **числовая запутанность** (запутанность по числу частиц), а $\sum_cp_cE(\hat\psi_c)$ —
**конфигурационная запутанность**. Для случайного состояния $2$ частиц на пяти узлах вычислим $E(\psi)$ по
частичным следам с обеих сторон и по сингулярным числам $4\times8$-матрицы $\Psi$:

```python
s = F.charges.index(2)                                      # the sector of two particles, n_q = 10
n = F.k_factors[s]
psi = torch.randn(n, dtype=dt)
psi = psi / psi.norm()
mat = torch.zeros(1, F.C, F.k_max, F.k_max, dtype=dt)
mat[0, s, :n, :n] = torch.outer(psi, psi.conj())
rho = DensityMatrix(F, matrix=mat)
r1, r2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)

psi_H = torch.zeros(4 * 8, dtype=dt)
psi_H[kron_rows(F)[s]] = psi                                # psi as a vector of H1 (x) H2
lam = torch.linalg.svdvals(psi_H.reshape(4, 8)) ** 2        # squared Schmidt coefficients
lam = lam[lam > 1e-14]
E_psi = -(lam * lam.log()).sum().item()
assert abs(r1.von_neumann_entropy().item() - E_psi) < 1e-10
assert abs(r2.von_neumann_entropy().item() - E_psi) < 1e-10  # Lemma 7.7

def entropy(sigma):
    w = torch.linalg.eigvalsh(sigma).clamp(min=1e-300)
    return -(w * w.log()).sum().item()

p = r1.sector_probabilities()[0]                            # c = 0, 1, 2 particles in the first part
assert torch.allclose(r2.sector_probabilities()[0, [2, 1, 0]], p)   # d(c) = 2 - c
H_num = -(p * p.log()).sum().item()
S_conf = sum(p[c].item() * entropy(r1.matrix[0, c, :k, :k] / p[c]) for c, k in enumerate(bos1.k_factors))
print(f"E = {E_psi:.4f} = H(p) {H_num:.4f} + configurational {S_conf:.4f}")
assert abs(E_psi - H_num - S_conf) < 1e-10                  # Proposition 7.8
```

> **Физический смысл** (физическая интерпретация, а не теорема этого курса). Стороны, связанные правилом суперотбора по числу частиц, могут использовать
> как ресурс лишь конфигурационную часть: это *запутанность частиц* H. M. Wiseman и J. A.
> Vaccaro (Phys. Rev. Lett. 91, 097902, 2003). $H(p)$ измеряют в холодных атомах как *числовую
> энтропию* (A. Lukin et al., Science 364, 256, 2019).

## 7.5 Случайные состояния и формула Пейджа

Пусть координаты $g\in\mathbb C^N$ независимы, а их вещественные и мнимые части независимы и распределены
по закону $N(0,\tfrac12)$ (стандартный комплексный гауссов вектор, предложение A.16(a)), так что $\mathbb E|g_i|^2=1$;
$g\neq0$ почти наверное (теорема A.20). **Случайное по Хаару чистое состояние** — это $\psi=g/\|g\|$. Его
распределение унитарно инвариантно: $Ug$ распределён так же, как $g$, для любой унитарной $U$ (предложение A.16(a)), и
$Ug/\|Ug\|=U\psi$. По предложению A.22 это единственное унитарно инвариантное распределение на единичной сфере, и,
как показано в его доказательстве, это распределение первого столбца распределённой по Хаару унитарной матрицы (теорема 5.3).

**Лемма 7.9.** $\|g\|$ и $\psi=g/\|g\|$ независимы.

*Доказательство.* (Это часть теоремы A.20(b).) Пусть $h\ge0$ и $\phi\ge0$ — борелевские функции на
$(0,\infty)$ и на единичной сфере $S^{2N-1}$ пространства $\mathbb C^N\cong\mathbb R^{2N}$. Поскольку $g$ имеет плотность
$\pi^{-N}e^{-|x|^2}$ (предложение A.16(a)), полярные координаты (теорема A.18 с $n=2N$ и определённой там
мерой $\mu$) дают
$$
\mathbb E\big[h(\|g\|)\phi(\psi)\big]=\pi^{-N}\!\int_{\mathbb R^{2N}}\!h(|x|)e^{-|x|^2}\phi\Big(\frac x{|x|}\Big)dx
=2N\pi^{-N}\int_0^\infty h(t)e^{-t^2}t^{2N-1}dt\int_{S^{2N-1}}\phi\,d\mu=:\alpha(h)\beta(\phi).
$$
Беря $\phi=1$, затем $h=1$, затем обе функции равными $1$, получаем $\mathbb Eh(\|g\|)=\alpha(h)\beta(1)$,
$\mathbb E\phi(\psi)=\alpha(1)\beta(\phi)$ и $\alpha(1)\beta(1)=1$, откуда
$\mathbb E[h(\|g\|)\phi(\psi)]=\mathbb Eh(\|g\|)\,\mathbb E\phi(\psi)$. $\square$

**Лемма 7.10.** $\mathbb E\|g\|^4=N(N+1)$. Если $G$ — $a\times b$-матрица с независимыми элементами, распределёнными
так же, как $g_i$ (матрица *Жинибра*), то $\mathbb E\operatorname{Tr}(GG^*)^2=ab(a+b)$.

*Доказательство.* При $g_i=x+iy$ с независимыми $x,y\sim N(0,\frac12)$ имеем $\mathbb Ex^4=3(\mathbb Ex^2)^2=\frac34$
(теорема A.24: три разбиения на пары), поэтому
$\mathbb E|g_i|^4=\mathbb Ex^4+2\mathbb Ex^2\mathbb Ey^2+\mathbb Ey^4=\tfrac34+\tfrac12+\tfrac34=2$ и
$\mathbb E\|g\|^4=\sum_i\mathbb E|g_i|^4+\sum_{i\ne j}\mathbb E|g_i|^2\mathbb E|g_j|^2=2N+N(N-1)$. Далее,
$\operatorname{Tr}(GG^*)^2=\sum_{i,i',j,j'}G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}$. Элементы и сопряжённые к ним
являются комплексными линейными комбинациями координат одного центрированного гауссова вектора (вещественных и
мнимых частей всех элементов), причём $\mathbb E[G_{ij}G_{kl}]=\mathbb E[\bar G_{ij}\bar G_{kl}]=0$ и
$\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}$ (предложение A.16(a)). По формуле Иссерлиса–Вика
(теорема A.24, комплексная форма, с тремя разбиениями четырёх сомножителей на пары)
$\mathbb E[G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}]=\delta_{ii'}\delta_{i'i}+0+\delta_{jj'}\delta_{jj'}=\delta_{ii'}+\delta_{jj'}$.
Суммирование по $i,i'\le a$ и $j,j'\le b$ даёт $ab^2+a^2b$. $\square$

**Теорема 7.11 (E. Lubkin, J. Math. Phys. 19, 1028, 1978).** Для случайного по Хаару
$\psi\in\mathbb C^a\otimes\mathbb C^b$ и $\rho_A=\operatorname{Tr}_2|\psi\rangle\langle\psi|$
$$
\mathbb E\operatorname{Tr}\rho_A^2=\frac{a+b}{ab+1}.
$$

*Доказательство.* Пусть $G_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,g\rangle$ — матрица Жинибра; матрица вектора $\psi$ есть
$\Psi=G/\|g\|$, и по предложению 7.2(a) $(\rho_A)_{\alpha\alpha'}=\sum_\beta\Psi_{\alpha\beta}\bar\Psi_{\alpha'\beta}$, т. е.
$\rho_A=\Psi\Psi^*=GG^*/\|g\|^2$. Следовательно,
$\operatorname{Tr}(GG^*)^2=\|g\|^4\operatorname{Tr}\rho_A^2$, где $\operatorname{Tr}\rho_A^2\in[0,1]$ — функция от
$\psi$. По лемме 7.9 $\mathbb E\operatorname{Tr}(GG^*)^2=\mathbb E\|g\|^4\,\mathbb E\operatorname{Tr}\rho_A^2$, и по лемме 7.10
$\mathbb E\operatorname{Tr}\rho_A^2=ab(a+b)/(ab(ab+1))$. $\square$

В главе 6 то же среднее получается из четвёртых моментов унитарно инвариантного единичного вектора
(лемма 6.14) как выходная чистота случайного канала (теорема 6.15 и замечание после неё).

**Следствие 7.12.** $\log\frac{ab+1}{a+b}\le\mathbb E\,S(\rho_A)\le\log\min(a,b)$.

*Доказательство.* По следствию 3.6 $S(\rho_A)\le\log a$, а по лемме 7.7
$S(\rho_A)=S(\operatorname{Tr}_1|\psi\rangle\langle\psi|)\le\log b$. Конечное неравенство Йенсена для вогнутого
$\log$ (решение упражнения 2 главы 3) с весами $\lambda_i>0$ — ненулевыми собственными значениями
$\rho_A$ — даёт $\sum_i\lambda_i\log\lambda_i\le\log\sum_i\lambda_i^2$, т. е. $S(\rho_A)\ge-\log X$ с
$X=\operatorname{Tr}\rho_A^2\in(0,1]$. Пусть $m=\mathbb EX=\frac{a+b}{ab+1}$ (теорема 7.11). Лемма 3.4 для выпуклой
функции $-\log$ на $(0,\infty)$ в точке $m$ даёт $-\log X\ge-\log m-(X-m)/m$; обе части ограничены
($X\ge1/a$ по неравенству Коши–Шварца), и взятие математических ожиданий даёт $\mathbb E[-\log X]\ge-\log m$. $\square$

При фиксированном $a$ и $b\to\infty$ границы сходятся: средняя энтропия малой части случайного чистого
состояния стремится к своему максимуму $\log a$. Точное среднее даётся формулой Пейджа, доказанной в приложении C.

**Теорема 7.13 (формула Пейджа).** При $a\le b$
$$
\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}.
$$
Кроме того, $\log a-\frac{a-1}{ab+1}-\frac{a-1}{2b}\le\mathbb E\,S(\rho_A)\le\log a-\frac{a-1}{2b}$, так что
$\big|\mathbb E\,S(\rho_A)-\log a+\frac a{2b}\big|\le\frac1{2b}$; для двух кубитов $\mathbb E\,S(\rho_A)=\frac13$.

*Доказательство.* Это теорема C.17 вместе со следствием C.18. Доказательство в приложении C исходит из
$\rho_A=GG^*/\operatorname{Tr}GG^*$ (доказательство теоремы 7.11), вычисляет совместную плотность собственных значений
матрицы Уишарта $GG^*$ (теорема C.13) и матрицы $\rho_A$ (теорема C.14) и интегрирует $x\log x$
по плотности уровней Лагерра (леммы C.15 и C.16). $\square$

Формула была высказана как гипотеза D. N. Page (Phys. Rev. Lett. 71, 1291, 1993) и доказана S. K.
Foong и S. Kanno (Phys. Rev. Lett. 72, 1148, 1994), J. Sánchez-Ruiz (Phys. Rev. E 52, 5653, 1995)
и S. Sen (Phys. Rev. Lett. 77, 1, 1996).

При $a>b$ роли $a$ и $b$ меняются местами. Действительно, унитарное отображение $e_\alpha\otimes e_\beta\mapsto e_\beta\otimes e_\alpha$
пространства $\mathbb C^a\otimes\mathbb C^b$ на $\mathbb C^b\otimes\mathbb C^a$ переставляет координаты $g$, поэтому оно переводит
случайное по Хаару чистое состояние $\psi$ в случайное по Хаару чистое состояние $\psi'$ пространства $\mathbb C^b\otimes\mathbb C^a$;
его матрица есть $\Psi^T$, так что $\operatorname{Tr}_2|\psi'\rangle\langle\psi'|=\Psi^T\bar\Psi=\operatorname{Tr}_1|\psi\rangle\langle\psi|$
(предложение 7.2(a)), а энтропия этого состояния равна $S(\rho_A)$ по лемме 7.7. Следовательно, при $a'=\min(a,b)$ и
$b'=\max(a,b)$ теорема 7.13 утверждает, что средняя энтропия запутанности меньше своего максимума $\log a'$ не менее чем на
$\frac{a'-1}{2b'}$ и не более чем на $\frac{a'-1}{2b'}+\frac{a'-1}{a'b'+1}$, а это меньше одного ната: это $0$ при $b'=1$, а
иначе $\frac{a'-1}{2b'}\le\frac12$ и $\frac{a'-1}{a'b'+1}<\frac1{b'}\le\frac12$. Таким образом, при фиксированном $ab$, как
функция от $a$, среднее отличается от $\log\min(a,b)=\min(\log a,\log(ab)-\log a)$ меньше чем на один нат.

Распределение $\rho_A$ — это **индуцированная мера** (название принадлежит K. Życzkowski и H.-J. Sommers,
J. Phys. A 34, 7111, 2001; теорема C.14 даёт плотность его собственных значений).
`states.random_density_matrix(alg, rank=r)` составляет $G_cG_c^*$ из независимых блоков Жинибра $G_c$
размера $k_c\times r$ и делит на $\sum_c\operatorname{Tr}G_cG_c^*$. Для одного блока это распределение $\rho_A$
выше с $b=r$ (доказательство теоремы 7.11). Для нескольких блоков это распределение плотности ограничения на
$\mathcal M_1$ векторного состояния случайного по Хаару чистого состояния $\psi=g/\|g\|$ пространства
$H_1\otimes\mathbb C^r$. Действительно, пусть $G_{i\beta}=\langle e_i\otimes e_\beta,g\rangle$ — матрица Жинибра размера
$D\times r$ ($D=\sum_ck_c$), а $G_c$ — её строки в блоке $c$: независимые матрицы Жинибра размера $k_c\times r$.
Возьмём $\mathcal M_2=M_r$ (один сектор) и постоянное правило слияния; тогда полный заряд единствен,
$\Pi=1$ и $\mathcal M_1\otimes_f\mathcal M_2=B(H_1\otimes\mathbb C^r)$. По теореме 7.6 ограничение имеет
плотность $E_1(\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|)$, а $\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|=GG^*/\|g\|^2$
(предложение 7.2(a), как в доказательстве теоремы 7.11), и её диагональные блоки равны $G_cG_c^*/\|g\|^2$,
где $\|g\|^2=\sum_c\operatorname{Tr}G_cG_c^*$.

**Предложение 7.14.** Для $\rho$, выбранной с помощью `random_density_matrix(alg, rank=r)` на $\bigoplus_cM_{k_c}$,
$D=\sum_ck_c$,
$$
\mathbb E\operatorname{Tr}\rho^2=\frac{\sum_ck_c(k_c+r)}{D(Dr+1)} .
$$

*Доказательство.* $\operatorname{Tr}\rho^2=\sum_c\operatorname{Tr}(G_cG_c^*)^2/\|g\|^4$, где $g\in\mathbb C^{Dr}$ — набор
всех элементов; это функция от $g/\|g\|$ со значениями в $[0,1]$. Как в теореме 7.11, лемма 7.9 даёт
$\sum_c\mathbb E\operatorname{Tr}(G_cG_c^*)^2=\mathbb E\|g\|^4\,\mathbb E\operatorname{Tr}\rho^2$, а лемма 7.10 даёт
$\mathbb E\operatorname{Tr}(G_cG_c^*)^2=k_cr(k_c+r)$ и $\mathbb E\|g\|^4=Dr(Dr+1)$, так что
$\mathbb E\operatorname{Tr}\rho^2=\sum_ck_cr(k_c+r)/(Dr(Dr+1))$. $\square$

Проверим теоремы 7.11 и 7.13 на $10^4$ чистых состояниях $\mathbb C^3\otimes\mathbb C^6$, редуцированных с помощью
`partial_trace`, и на индуцированной мере, а предложение 7.14 — для $M_1\oplus M_2\oplus M_3$, $r=2$:

```python
a, b, n_samples = 3, 6, 10_000
g = torch.randn(n_samples, a * b, dtype=dt)                 # complex Gaussian: E|g_i|^2 = 1
psi = g / g.norm(dim=-1, keepdim=True)
AB_alg = tensor_product(TypeIAlgebra([a], [a], precision='double'), TypeIAlgebra([b], [b], precision='double'))
pure = DensityMatrix(AB_alg, matrix=(psi[:, :, None] * psi.conj()[:, None, :])[:, None],
                     validate=False)                        # pure by construction
rho_A = partial_trace(pure, keep=1)
purity_exact = (a + b) / (a * b + 1)
page = sum(1 / k for k in range(b + 1, a * b + 1)) - (a - 1) / (2 * b)

def mean_se(x):
    return x.mean().item(), x.std().item() / math.sqrt(x.numel())

for name, rho_s in (('reduced pure states', rho_A),
                    ('induced measure', states.random_density_matrix(rho_A.algebra, 40_000, rank=b))):
    (m_p, e_p), (m_S, e_S) = mean_se(rho_s.purity()), mean_se(rho_s.von_neumann_entropy())
    print(f"{name}: E Tr rho^2 = {m_p:.4f} +- {e_p:.4f} (exact {purity_exact:.4f}), "
          f"E S = {m_S:.4f} +- {e_S:.4f} (Page {page:.4f})")
    assert abs(m_p - purity_exact) < 5 * e_p and abs(m_S - page) < 5 * e_S     # Theorems 7.11, 7.13
    assert math.log((a * b + 1) / (a + b)) < m_S < math.log(a)                # Corollary 7.12
assert math.log(a) - (a - 1) / (a * b + 1) - (a - 1) / (2 * b) <= page <= math.log(a) - (a - 1) / (2 * b)
assert abs(sum(1 / k for k in (3, 4)) - 1 / 4 - 1 / 3) < 1e-15                # two qubits: E S = 1/3

alg3 = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')
m_p, e_p = mean_se(states.random_density_matrix(alg3, 40_000, rank=2).purity())
print(f"E Tr rho^2 = {m_p:.4f} +- {e_p:.4f} (exact {26 / 78:.4f})")
assert abs(m_p - 26 / 78) < 5 * e_p                                          # Proposition 7.14
```

> **Физический смысл** (физическая интерпретация, а не теорема этого курса). Средняя
> энтропия запутанности случайного чистого состояния как функция размера части при фиксированной полной
> размерности (кривая, описанная после теоремы 7.13) называется *кривой Пейджа*. В физике её используют,
> чтобы описать, как информация покидает чёрную дыру с излучением Хокинга, и почему части замкнутой системы
> выглядят тепловыми; это интерпретации, а не результаты этого курса.

## 7.6 С чем сравнивать

На одном гильбертовом пространстве QuTiP (`Qobj.ptrace`, `entropy_vn`, `rand_ket`, `rand_dm` с распределением
Жинибра и заданным рангом) и `numpy` (reshape и `einsum`, как выше) вычисляют частичные следы,
энтропии и индуцированные состояния для обычных произведений (проверки с `einsum` выше сравнивают такое
вычисление с библиотекой). Библиотеки
тензорных сетей с сохраняющимися зарядами — ITensor (блочно-разреженные тензоры с квантовыми числами), TeNPy (заряженные
ноги, спектры запутанности по зарядам) и, для неабелевых симметрий, TensorKit.jl (деревья слияния) —
сливают секторы тем же способом для больших систем и дают запутанность с разрешением по симметрии.
Библиотека добавляет алгебраический уровень: произведения являются алгебрами фон Неймана с секторами и следами,
правило слияния — любое отображение меток, частичный след — канал с сопряжённым для каждого следа
и правильным ограничением для произведений со слиянием (теоремы 7.3, 7.6), и всё это работает пакетно. В ней
нет тензорных сетей, а `fused_tensor_product` принимает правило с одним исходом для каждой пары зарядов, как
для абелевых симметрий; неабелево правило, как для $SU(2)$, где пара спинов $j_1,j_2$ складывается в
несколько полных спинов $|j_1-j_2|,\dots,j_1+j_2$ (глава 4, упражнение 2), здесь произведением со слиянием не является (библиотека рассматривает $SU(2)$ для цепочек
спинов $\frac12$ отдельно, глава 4).

## Упражнения

1. Вычислите сопряжённое к $\operatorname{Tr}_2$ отображение относительно $\tau_{\mathrm{vN}}$ как на
   $\mathcal M_1\otimes\mathcal M_2$, так и на $\mathcal M_1$, и проверьте его с помощью
   `partial_trace_channel(alg12).adjoint('tau_vN')`.
2. (Субаддитивность.) Пусть $\rho$ — состояние на $\mathcal M_1\otimes\mathcal M_2$ с обратимыми редуцированными
   состояниями $\rho_1,\rho_2$. Покажите, что $S(\rho)\le S(\rho_1)+S(\rho_2)$, причём равенство имеет место тогда и только тогда, когда
   $\rho=\rho_1\otimes\rho_2$.
3. Пусть `site` — это $\mathbb C\oplus\mathbb C$ с зарядами $0,1$ (пустой или занятый узел), слитое с
   самим собой по сложению. Вычислите $E(\psi)$ и две его части для $\psi=(|01\rangle+|10\rangle)/\sqrt2$ — одной
   частицы, разделённой между двумя узлами, — и опишите ограничение векторного состояния
   $(|00\rangle+|11\rangle)/\sqrt2$ на алгебру со слиянием.

<details>
<summary>Решения</summary>

1. Веса равны $1/(C_1k_c)$ на $\mathcal M_1$ и $1/(C_1C_2k_cl_d)$ на произведении. Для
   $A\in\mathcal M_1$ по предложению 7.2(a)
   $\tau_{\mathrm{vN}}((\operatorname{Tr}_2X)A)=\sum_c\frac1{C_1k_c}\sum_d\operatorname{Tr}(X_{cd}(A_c\otimes1))
   =\sum_{c,d}\frac1{C_1C_2k_cl_d}\operatorname{Tr}\big(X_{cd}\,C_2l_d(A_c\otimes1)\big)$, так что сопряжённое отображение есть
   $A\mapsto\bigoplus_{c,d}C_2l_d\,A_c\otimes1_{l_d}=A\otimes Z$ с центральным элементом
   $Z=\bigoplus_dC_2l_d1_{l_d}$ алгебры $\mathcal M_2$. Оно унитально лишь при $Z=1$, т. е. $\mathcal M_2=\mathbb C$,
   и не унитально в примере ($Z=4\cdot1_2\oplus6\cdot1_3$). Для `Tr_norm` то же вычисление даёт
   $A\otimes\bigoplus_dl_d1_{l_d}$, что унитально лишь при всех $l_d=1$. Следы `Tr_blunt` всегда делают
   $\operatorname{Tr}_2$ сопряжённым к $A\mapsto A\otimes1$; в общем случае, по тому же вычислению, это верно для
   $\operatorname{tr}_w$ на $\mathcal M_1$ и весов $w_{cd}=w_c$ на произведении.

```python
dual_tau = composite.partial_trace_channel(alg12, keep=1).adjoint('tau_vN')
A = random_operator(alg1)
Z = alg2.central([alg2.C * l for l in alg2.k_factors])
assert close(dual_tau(A), kron(A, Z, alg12))
```

2. $\rho_1\otimes\rho_2$ — состояние (теорема 7.1(c)), обратимое, с собственными векторами $u_i\otimes v_j$ и
   собственными значениями $\lambda_i\mu_j$, поэтому $\log(\rho_1\otimes\rho_2)=\log\rho_1\otimes1+1\otimes\log\rho_2$
   (теорема 1.5(a)). По определению частичного следа
   $\operatorname{Tr}\rho(\log\rho_1\otimes1)=\operatorname{Tr}\rho_1\log\rho_1$, и аналогично для второго
   сомножителя. Следовательно, $S(\rho\|\rho_1\otimes\rho_2)=S(\rho_1)+S(\rho_2)-S(\rho)$, а это $\ge0$, причём
   равенство имеет место тогда и только тогда, когда $\rho=\rho_1\otimes\rho_2$ (следствие 3.6). Эта разность — взаимная информация.

```python
rho = states.random_density_matrix(alg12, batch_size=4)
rho1, rho2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
mutual = rho1.von_neumann_entropy() + rho2.von_neumann_entropy() - rho.von_neumann_entropy()
assert (mutual > 0).all()
assert torch.allclose(mutual, states.relative_entropy(rho, kron(rho1, rho2, alg12)))
```

3. Секторы — это $N=0,1,2$ размеров $1,2,1$. $\psi$ лежит в $N=1$, причём $\psi_0=|01\rangle/\sqrt2$,
   $\psi_1=|10\rangle/\sqrt2$, так что $p=(\frac12,\frac12)$, и, поскольку блоки одномерны,
   $E(\hat\psi_c)=0$: по предложению 7.8 $E(\psi)=H(p)=\log2$ целиком является числовой запутанностью, а
   конфигурационная часть равна $0$. Вектор $(|00\rangle+|11\rangle)/\sqrt2$ — суперпозиция $N=0$ и
   $N=2$; на блочно-диагональной алгебре со слиянием его векторное состояние — смесь
   $\frac12|00\rangle\langle00|+\frac12|11\rangle\langle11|$ состояний-произведений (глава 3, упражнение 3).

```python
site = TypeIAlgebra([1, 1], [1, 1], precision='double', charges=[0, 1])
S2 = fused_tensor_product(site, site)
assert S2.k_factors == [1, 2, 1]
mat = torch.zeros(1, 3, 2, 2, dtype=dt)
mat[0, 1] = 0.5                                       # |psi><psi| with psi = (|01> + |10>)/sqrt 2
r1 = partial_trace(DensityMatrix(S2, matrix=mat), keep=1)
assert torch.allclose(r1.sector_probabilities()[0], torch.tensor([0.5, 0.5], dtype=torch.float64))
assert abs(r1.von_neumann_entropy().item() - math.log(2)) < 1e-12
```

</details>

**Далее:** в [главе 8](08_dynamics.md) рассматривается динамика: уравнения Шрёдингера и Линдблада в
алгебре с секторами, переходы между секторами и постселекция.
