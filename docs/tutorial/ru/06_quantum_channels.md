[English](../06_quantum_channels.md) · [Español](../es/06_quantum_channels.md) · [Français](../fr/06_quantum_channels.md) · [Deutsch](../de/06_quantum_channels.md) · [中文](../zh/06_quantum_channels.md) · [日本語](../ja/06_quantum_channels.md) · **Русский**

# 6. Квантовые каналы

Квантовый канал — это вполне положительное линейное отображение, сохраняющее след. Примеры: сопряжение
унитарным элементом, условное ожидание на центр, отображение $\rho\mapsto\operatorname{Tr}_E(V\rho V^*)$, заданное
изометрией $V$, и отображения, переносящие вес из одного сектора в другой; их физическое
прочтение (эволюция, измерение, результат которого забыт, окружение, потеря
частицы) обсуждается только в блоках *Физический смысл*. В этой главе
доказываются стандартные описания таких отображений — матрицы Чоя, операторы Крауса, дилатация
Стайнспринга — для отображений между конечномерными алгебрами фон Неймана, а затем изучается то, что
специфично для прямых сумм: отображения между секторами, сопряжённые относительно следов из главы 3
и классическая часть канала. В главах 2 и 3 каналы библиотеки уже использовались; здесь
излагается их теория.

**Вы узнаете**, почему положительности недостаточно, как переходить между матрицами Чоя, операторами
Крауса, супероператорами и изометриями, почему «сохранение следа» зависит от следа, как только отображение
меняет сектор, какие операторы унитальный канал оставляет неподвижными и что случайный канал в среднем делает с
чистым состоянием.

Всюду далее $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ и
$\mathcal N=\bigoplus_{d=1}^{C'}M_{k'_d}(\mathbb C)$ с минимальными центральными проекторами $P_c$ и
$Q_d$; $A=\bigoplus_cA_c$; $M_{k'\times k}$ — пространство комплексных матриц размера $k'\times k$;
$E_{ij}=|i\rangle\langle j|$ — матричные единицы; $\operatorname{Tr}$ — след `Tr_blunt`
$\sum_c\operatorname{Tr}A_c$, а $\operatorname{tr}_w=\sum_cw_c\operatorname{Tr}A_c$ — точный
след (теорема 3.1). Блочные матрицы действуют на $\mathbb C^n\otimes\mathbb C^k$, блочный индекс идёт первым.

## 6.1 Вполне положительные отображения и теорема Чоя

Линейное отображение $\Phi:\mathcal M\to\mathcal N$ **положительно**, если из $A\ge0$ следует $\Phi(A)\ge0$.
Матрицы $[A_{ij}]$ размера $n\times n$ с элементами из $\mathcal M$ снова образуют прямую сумму матричных
алгебр, $M_n(\mathcal M)\cong\bigoplus_cM_n(M_{k_c})$ посредством $[A_{ij}]\mapsto([(A_{ij})_c]_{ij})_c$,
так что $[A_{ij}]\ge0$ означает, что каждая блочная матрица $[(A_{ij})_c]\in M_{nk_c}$ положительна. $\Phi$
называется **$n$-положительным**, если $\Phi_n:[A_{ij}]\mapsto[\Phi(A_{ij})]$ положительно, и **вполне положительным**
(CP), если оно $n$-положительно при всех $n$. Оно **сохраняет след** (TP), если
$\operatorname{Tr}\Phi(A)=\operatorname{Tr}A$, **унитально**, если $\Phi(1)=1$, и называется **каналом**, если оно
CP и TP. Мы постоянно используем то, что из $X\ge0$ следует $TXT^*\ge0$ для любой прямоугольной $T$, поскольку
$\langle y,TXT^*y\rangle=\langle T^*y,XT^*y\rangle$.

**Предложение 6.1.** При $k\ge2$ транспонирование $\theta(A)=A^T$ на $M_k$ положительно, но не
$2$-положительно.

*Доказательство.* Для $A\ge0$ имеем $\langle x,A^Tx\rangle=\sum_{ij}\bar x_iA_{ji}x_j=\langle\bar x,A\bar x\rangle\ge0$.
Для $\Omega=e_1\otimes e_1+e_2\otimes e_2\in\mathbb C^2\otimes\mathbb C^k$
$|\Omega\rangle\langle\Omega|=[E_{ij}]_{i,j=1}^2\ge0$, а $\theta_2$ переводит этот оператор в
$F=[E_{ji}]_{i,j=1}^2$, который переводит $e_1\otimes e_2-e_2\otimes e_1$ в противоположный вектор. $\square$

**Теорема 6.2 (Чой).** Для линейного отображения $\Phi:M_k\to M_{k'}$ следующие условия равносильны:
(a) $\Phi$ является CP; (b) $\Phi$ является $k$-положительным; (c) **матрица Чоя**
$$
J(\Phi)=\sum_{i,j=1}^kE_{ij}\otimes\Phi(E_{ij})=[\Phi(E_{ij})]_{i,j=1}^k\in M_k(M_{k'})
$$
положительна; (d) существуют $K_1,\dots,K_r\in M_{k'\times k}$ (**операторы Крауса**), для которых
$\Phi(A)=\sum_lK_lAK_l^*$, причём можно взять $r=\operatorname{rank}J(\Phi)\le kk'$.

*Доказательство.* (a)$\Rightarrow$(b) тривиально. (b)$\Rightarrow$(c): для
$\Omega=\sum_ie_i\otimes e_i\in\mathbb C^k\otimes\mathbb C^k$
$|\Omega\rangle\langle\Omega|=[E_{ij}]\ge0$ и $J(\Phi)=\Phi_k([E_{ij}])$. (c)$\Rightarrow$(d): по
теореме 1.3 $J(\Phi)=\sum_{l=1}^rv_lv_l^*$ с $r=\operatorname{rank}J(\Phi)$ (собственные векторы
ненулевых собственных значений, умноженные на квадратные корни из них). Всякий $v\in\mathbb C^k\otimes\mathbb C^{k'}$
записывается в виде $v=\sum_ie_i\otimes Ke_i$ ровно с одной $K\in M_{k'\times k}$ ($Ke_i$ — $i$-й блок
$v$), и $(i,j)$-блок $vv^*$ равен $Ke_ie_j^*K^*=KE_{ij}K^*$. Записав так $v_l$ с
$K_l$ и сравнив блоки, получаем $\Phi(E_{ij})=\sum_lK_lE_{ij}K_l^*$, откуда (d) следует по линейности.
(d)$\Rightarrow$(a): $\Phi_n([A_{ij}])=\sum_l(1_n\otimes K_l)[A_{ij}](1_n\otimes K_l)^*\ge0$. $\square$

**Следствие 6.3.** Если $\Phi(A)=\sum_{l=1}^sL_lAL_l^*$, то $s\ge\operatorname{rank}J(\Phi)$. Поэтому
минимальное число операторов Крауса, **ранг Крауса**, равно $\operatorname{rank}J(\Phi)$.

*Доказательство.* Если прочитать последнее вычисление в обратном порядке, $J(\Phi)=\sum_{l=1}^sw_lw_l^*$ с
$w_l=\sum_ie_i\otimes L_le_i$ — сумма $s$ матриц ранга $\le1$. $\square$

Отображение между прямыми суммами — это матрица отображений между блоками. Пусть $\iota_c:M_{k_c}\to\mathcal M$
помещает матрицу в блок $c$, $\pi_d:\mathcal N\to M_{k'_d}$ берёт блок $d$, и
$\Phi_{dc}=\pi_d\circ\Phi\circ\iota_c:M_{k_c}\to M_{k'_d}$. Поскольку $A=\sum_c\iota_c(A_c)$,
$\Phi(A)_d=\sum_c\Phi_{dc}(A_c)$.

**Предложение 6.4.** Отображение $\Phi:\mathcal M\to\mathcal N$ является $n$-положительным тогда и только тогда, когда таково каждое $\Phi_{dc}$. Следовательно,
$\Phi$ является CP тогда и только тогда, когда существуют $K^{dc}_l\in M_{k'_d\times k_c}$, $l=1,\dots,r$, такие, что
$$
\Phi(A)_d=\sum_c\sum_{l=1}^rK^{dc}_l\,A_c\,K^{dc*}_l ,
$$
и достаточно $r=\max_{d,c}\operatorname{rank}J(\Phi_{dc})$. $\Phi$ является TP тогда и только тогда, когда
$\sum_{d,l}K^{dc*}_lK^{dc}_l=1_{k_c}$ для каждого $c$.

*Доказательство.* Если все $\Phi_{dc}$ $n$-положительны и $[A_{ij}]\ge0$, то каждая $X^c=[(A_{ij})_c]\ge0$, и
блок $d$ матрицы $\Phi_n([A_{ij}])$ равен $\sum_c(\Phi_{dc})_n(X^c)\ge0$. Обратно, для $[X_{ij}]\ge0$ в
$M_n(M_{k_c})$ имеем $[\iota_c(X_{ij})]\ge0$ в $M_n(\mathcal M)$ (её блок $c$ равен $[X_{ij}]$, остальные
нулевые), а блок $d$ её образа равен $[\Phi_{dc}(X_{ij})]$. Форма Крауса следует из теоремы 6.2
для каждой пары $(d,c)$, дополненной нулевыми операторами до общего $r$. Наконец, в силу цикличности
$\operatorname{Tr}\Phi(A)=\sum_c\operatorname{Tr}\big(A_c\sum_{d,l}K^{dc*}_lK^{dc}_l\big)$, и это равно
$\sum_c\operatorname{Tr}A_c$ для всех $A$ тогда и только тогда, когда каждая сумма равна $1_{k_c}$ (возьмите матричные единицы:
$\operatorname{Tr}(E_{ij}X)=X_{ji}$). $\square$

В библиотеке два класса. `Channel` хранит **сохраняющие секторы** отображения ($C'=C$, $k'=k$,
$\Phi_{dc}=0$ при $d\ne c$), операторы Крауса которых $K_l=\bigoplus_cK^{cc}_l$ лежат в $\mathcal M$, в виде
тензора формы `(batch, r, C, k_max, k_max)`; `InterSectorChannel` хранит общее отображение из
предложения 6.4 в виде тензора формы `(batch, r, C', C, k'_max, k_max)`. `Channel.choi()` возвращает
$J(\Phi_{cc})$ для каждого сектора, а `Channel.superoperator()` — матрицу $S_c$ отображения $\Phi_{cc}$ на
векторах, развёрнутых по строкам, $\operatorname{vec}(X)_{(a,b)}=X_{ab}$. Обе содержат числа
$\Phi(E_{ij})_{ab}=J_{(i,a),(j,b)}=S_{(a,b),(i,j)}$, но в разном порядке. `Channel.from_superoperator`
выбрасывает `ValueError`, если матрица Чоя не эрмитова, $\|J-J^*\|>2\cdot10^{-10}s$ (в двойной точности), где
$s=\max(1,\|J\|)$ со спектральной нормой; в противном случае он заменяет $J$ его эрмитовой частью $(J+J^*)/2$,
диагонализует её, как в доказательстве теоремы 6.2, выбрасывает `ValueError`, если какое-то собственное значение меньше $-10^{-10}s$,
и обнуляет остальные отрицательные собственные значения. Если $\Phi(X^*)=\Phi(X)^*$ для всех $X$, то $J^*=J$ (блок
$(i,j)$ матрицы $J^*$ равен $\Phi(E_{ji})^*=\Phi(E_{ij})$); обратно, из $J^*=J$ следует $\Phi(E_{ji})^*=\Phi(E_{ij})$,
а значит, по линейности $\Phi(X^*)=\Phi(X)^*$. Итак, с точностью до этих допусков на округление он отвергает в точности
отображения, не сохраняющие сопряжение, а среди остальных, по теореме 6.2, в точности те, что не являются CP;
например, $X\mapsto iX$ отвергается.

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

> **Физический смысл** (физическая интерпретация, а не теорема этого курса). Система никогда не бывает одна: она может быть скоррелирована с нетронутой
> вспомогательной системой с $n$ уровнями. Физический процесс действует на пару как $\Phi_n=\mathrm{id}_n\otimes\Phi$ и должен
> переводить состояния в состояния, поэтому $\Phi$ должно быть вполне положительным. Транспонирование положительно, но
> переводит запутанное состояние $|\Omega\rangle\langle\Omega|/2$ в оператор с собственным значением $-1/2$.
> Матрица Чоя — это образ максимально запутанного состояния: одно состояние решает вопрос о полной положительности.

## 6.2 Дилатация Стайнспринга и свобода в выборе операторов Крауса

Для $\Phi(A)=\sum_{l=1}^rK_lAK_l^*$ из $M_k$ в $M_{k'}$ положим $\Phi^*(B)=\sum_lK_l^*BK_l$; в силу
цикличности $\operatorname{Tr}(\Phi(A)B)=\operatorname{Tr}(A\Phi^*(B))$, так что $\Phi^*$ — сопряжённое к
$\Phi$ относительно следа (раздел 6.3) и не зависит от выбора $K_l$ (если
$\operatorname{Tr}(AX)=\operatorname{Tr}(AX')$ для всех $A$, то $X=X'$: возьмите $A=E_{ji}$). На
$\mathbb C^r\otimes\mathbb C^{k'}$ (**окружение** $\mathbb C^r$ идёт первым) положим
$\operatorname{Tr}_EX=\sum_l(e_l^*\otimes1)X(e_l\otimes1)$.

**Теорема 6.5 (Стайнспринг, конечномерный случай).** Пусть
$V=\sum_le_l\otimes K_l:\mathbb C^k\to\mathbb C^r\otimes\mathbb C^{k'}$. Тогда
$$
\Phi(A)=\operatorname{Tr}_E(VAV^*),\qquad\Phi^*(B)=V^*(1_r\otimes B)V,
$$
и $\Phi$ является TP тогда и только тогда, когда $V$ — изометрия. Обратно, всякое $V$ задаёт CP-отображение первой формулой.

*Доказательство.* $(e_l^*\otimes1)V=K_l$, поэтому $\operatorname{Tr}_E(VAV^*)=\sum_lK_lAK_l^*$; обратно,
$\operatorname{Tr}_E(VAV^*)$ — форма Крауса с $K_l=(e_l^*\otimes1)V$. Далее,
$V^*(1\otimes B)V=\sum_{l,l'}\langle e_l,e_{l'}\rangle K_l^*BK_{l'}=\Phi^*(B)$ и
$V^*V=\sum_lK_l^*K_l$, что равно $1$ тогда и только тогда, когда $\Phi$ является TP (предложение 6.4). $\square$

Поэтому наименьшее $r$, при котором $\Phi(A)=\operatorname{Tr}_E(VAV^*)$ для некоторого
$V:\mathbb C^k\to\mathbb C^r\otimes\mathbb C^{k'}$, — это ранг Крауса $\operatorname{rank}J(\Phi)$: всякое
такое $V$ даёт форму Крауса с $r$ операторами, так что $r\ge\operatorname{rank}J(\Phi)$ (следствие 6.3), а
$r=\operatorname{rank}J(\Phi)$ достигается по теореме 6.2(d).

Формула $\Phi^*(B)=V^*\pi(B)V$ с представлением $\pi(B)=1_r\otimes B$ — это та форма, в которой
W. F. Stinespring (Proc. AMS 6 (1955) 211) сформулировал теорему, причём в гораздо более общей ситуации,
чем доказанный здесь конечномерный случай (историческая атрибуция; для дальнейшего чтения: V. Paulsen,
*Completely Bounded Maps and Operator Algebras*, Cambridge 2002, Ch. 4). Для отображения между прямыми
суммами пусть $\iota'_d$ — вложение $d$-го слагаемого суммы $\bigoplus_d\mathbb C^{k'_d}$ и
$V_c=\sum_{d,l}e_l\otimes\iota'_dK^{dc}_l:\mathbb C^{k_c}\to\mathbb C^r\otimes\bigoplus_d\mathbb C^{k'_d}$.
Тогда $(e_l^*\otimes1)V_c=\sum_d\iota'_dK^{dc}_l$, поэтому блок $(d,d)$ оператора
$\sum_c\operatorname{Tr}_E(V_cA_cV_c^*)=\sum_{c,l}\big(\sum_d\iota'_dK^{dc}_l\big)A_c\big(\sum_{d'}\iota'_{d'}K^{d'c}_l\big)^*$
равен $\sum_{c,l}K^{dc}_lA_cK^{dc*}_l=\Phi(A)_d$; а так как $\iota'^*_d\iota'_{d'}=\delta_{dd'}1$, то
$V_c^*V_c=\sum_{d,l}K^{dc*}_lK^{dc}_l$, и $\Phi$ является TP тогда и только тогда, когда каждое $V_c$ — изометрия
(предложение 6.4). Именно так `random_channel` и `random_inter_sector_channel` генерируют каналы
(раздел 6.5).

**Лемма 6.6.** Если $X,Y\in M_{p\times n}$ и $XX^*=YY^*$, то $Y=XW$ для некоторой унитарной $W\in M_n$.

*Доказательство.* $\|X^*z\|^2=\langle z,XX^*z\rangle=\|Y^*z\|^2$ для всех $z$, поэтому $X^*z\mapsto Y^*z$ —
корректно определённая линейная изометрия образа $X^*$ на образ $Y^*$. Как в доказательстве
теоремы 1.8, она продолжается до унитарной $U$ с $UX^*=Y^*$; возьмём $W=U^*$. $\square$

**Теорема 6.7 (унитарная свобода).** Пусть $K_1,\dots,K_n$ и $L_1,\dots,L_n$ лежат в $M_{k'\times k}$
(более короткое семейство дополняется нулями). Тогда $\sum_lK_lAK_l^*=\sum_jL_jAL_j^*$ для всех $A$ тогда и только тогда, когда
$L_j=\sum_lu_{jl}K_l$ для некоторой унитарной $u\in M_n$.

*Доказательство.* Если $L_j=\sum_lu_{jl}K_l$, то
$\sum_jL_jAL_j^*=\sum_{l,l'}(u^*u)_{l'l}K_lAK_{l'}^*=\sum_lK_lAK_l^*$. Обратно, у равных отображений
равны матрицы Чоя: $\sum_lv_lv_l^*=\sum_jw_jw_j^*$ с $v_l=\sum_ie_i\otimes K_le_i$ и
$w_j=\sum_ie_i\otimes L_je_i$. Лемма 6.6, применённая к матрицам со столбцами $v_l$ и $w_j$, даёт
$w_j=\sum_lW_{lj}v_l$; так как $K\mapsto\sum_ie_i\otimes Ke_i$ линейно и инъективно,
$L_j=\sum_lW_{lj}K_l$, и $u=W^T$ унитарна. $\square$

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

> **Физический смысл** (физическая интерпретация, а не теорема этого курса). Всякий канал — это
> изометрическое взаимодействие системы с окружением,
> находящимся в фиксированном состоянии, после которого окружение игнорируется. Измерение окружения в базисе
> $e_l$ даёт исход $l$ с вероятностью $\operatorname{Tr}K_l\rho K_l^*$ и оставляет состояние
> $K_l\rho K_l^*/\operatorname{Tr}(K_l\rho K_l^*)$: операторы Крауса — это возможные «квантовые скачки».
> Другой базис даёт другие скачки, но тот же канал (теорема 6.7), а наименьшее окружение
> имеет размерность, равную рангу Крауса, как показано после теоремы 6.5.

## 6.3 Сопряжённые относительно следа

**Определение.** Пусть $\operatorname{tr}_w$, $\operatorname{tr}_v$ — точные следы на $\mathcal M$,
$\mathcal N$. **Сопряжённым** к линейному отображению $\Phi:\mathcal M\to\mathcal N$ относительно $(w,v)$ называется
линейное отображение $\Phi^\sharp:\mathcal N\to\mathcal M$, для которого
$$
\operatorname{tr}_v\big(\Phi(x)\,y\big)=\operatorname{tr}_w\big(x\,\Phi^\sharp(y)\big)\qquad(x\in\mathcal M,\ y\in\mathcal N).
$$
Оно существует и единственно: при фиксированном $y$ отображение $x\mapsto\operatorname{tr}_v(\Phi(x)y)$ — линейный функционал,
поэтому по теореме 3.2 он равен $\operatorname{tr}_w(\rho x)=\operatorname{tr}_w(x\rho)$ ровно для одного
$\rho=:\Phi^\sharp(y)$, а из единственности следует линейность $\Phi^\sharp$. Для следов `Tr_blunt` пишем $\Phi^*$.
Пусть $z_w=\sum_cw_cP_c$ и $z_v=\sum_dv_dQ_d$ (центральные и обратимые).

**Предложение 6.8.** Если у $\Phi$ операторы Крауса $K^{dc}_l$ (предложение 6.4), то
$$
\Phi^\sharp(y)_c=\sum_d\frac{v_d}{w_c}\sum_lK^{dc*}_l\,y_d\,K^{dc}_l,\qquad\text{i.e.}\qquad
\Phi^\sharp(y)=z_w^{-1}\,\Phi^*(z_v\,y).
$$
В частности, $\Phi^\sharp$ является CP с операторами Крауса $(v_d/w_c)^{1/2}K^{dc*}_l$.

*Доказательство.* В силу цикличности следа в каждом блоке
$$
\operatorname{tr}_v(\Phi(x)y)=\sum_dv_d\sum_{c,l}\operatorname{Tr}\big(K^{dc}_lx_cK^{dc*}_ly_d\big)
=\sum_cw_c\operatorname{Tr}\Big(x_c\sum_d\frac{v_d}{w_c}\sum_lK^{dc*}_ly_dK^{dc}_l\Big),
$$
и единственность даёт формулу. При $v=w=1$ она принимает вид $\Phi^*(y)_c=\sum_{d,l}K^{dc*}_ly_dK^{dc}_l$;
множители $v_d$ и $1/w_c$ — это умножение на $z_v$ и $z_w^{-1}$. $\square$

$\Phi$ **сохраняет след относительно $(w,v)$**, если $\operatorname{tr}_v\Phi(x)=\operatorname{tr}_w(x)$ для всех $x$.

**Следствие 6.9.** (a) $\Phi$ является TP относительно $(w,v)$ тогда и только тогда, когда $\Phi^\sharp(1)=1$, т. е.
$\sum_{d,l}(v_d/w_c)K^{dc*}_lK^{dc}_l=1_{k_c}$ для каждого $c$. (b) $\Phi$ унитально тогда и только тогда, когда $\Phi^\sharp$
является TP относительно $(v,w)$. (c) Сопряжённое к $\Phi^\sharp$ относительно $(v,w)$ есть $\Phi$. (d) Если
$\Phi:\mathcal M\to\mathcal M$ сохраняет секторы и $v=w$, то $\Phi^\sharp=\Phi^*$ при любом $w$,
и $\Phi$ является TP для одного следа тогда и только тогда, когда оно TP для всех.

*Доказательство.* (a) $\operatorname{tr}_v\Phi(x)=\operatorname{tr}_w(x\Phi^\sharp(1))$, и
$\operatorname{tr}_w(x(\Phi^\sharp(1)-1))=0$ для всех $x$ тогда и только тогда, когда $\Phi^\sharp(1)=1$ (возьмите
$x=(\Phi^\sharp(1)-1)^*$; $\operatorname{tr}_w$ точен); затем примените предложение 6.8. (b)
$\operatorname{tr}_w\Phi^\sharp(y)=\operatorname{tr}_v(\Phi(1)y)$ равно $\operatorname{tr}_v(y)$ для всех
$y$ тогда и только тогда, когда $\Phi(1)=1$. (c) $\operatorname{tr}_w(\Phi^\sharp(y)x)=\operatorname{tr}_v(y\Phi(x))$ по свойству
следа; воспользуйтесь единственностью. (d) Встречаются только пары $d=c$, для которых $v_d/w_c=1$. $\square$

**Плотности в картине Гейзенберга.** Пусть $\Psi:\mathcal N\to\mathcal M$ линейно, $\omega$ — линейный
функционал на $\mathcal M$ с плотностью $\rho$ относительно $\operatorname{tr}_w$ (теорема 3.2), а
$\Psi^\sharp:\mathcal M\to\mathcal N$ — сопряжённое к $\Psi$ относительно $(v,w)$. Тогда
$\omega(\Psi(y))=\operatorname{tr}_w(\Psi(y)\rho)=\operatorname{tr}_v(y\Psi^\sharp(\rho))=\operatorname{tr}_v(\Psi^\sharp(\rho)y)$,
так что $\omega\circ\Psi$ имеет плотность $\Psi^\sharp(\rho)$ относительно $\operatorname{tr}_v$; если $\Psi$ унитально,
то $\Psi^\sharp$ является TP относительно $(w,v)$ по следствию 6.9(b) (с переставленными ролями $\mathcal M$ и
$\mathcal N$).

Для отображений, меняющих сектор, следствие 6.9(d) неверно. Пусть $\mathcal M=\mathcal N=M_1\oplus M_2$ (секторы $0$ и $1$), и пусть $\Phi$ —
отображение, называемое **потерей частицы** (см. блок ниже),
$\Phi(x)=(x_0+\operatorname{Tr}x_1)\oplus0$, с операторами Крауса $K^{00}=1$, $K^{01}_l=\langle l|$.
Оно является TP для следов `Tr_blunt`. Для `Tr_norm` $w=v=(1,\tfrac12)$, $\Phi^\sharp(1)=1\oplus2\cdot1_2$
и $\operatorname{tr}_v\Phi(x)=x_0+\operatorname{Tr}x_1\ne x_0+\tfrac12\operatorname{Tr}x_1=\operatorname{tr}_w(x)$.
В библиотеке `Channel.adjoint()` — это $\Phi^*$, а `InterSectorChannel.adjoint(trace=...)` — это
$\Phi^\sharp$ для одного и того же именованного следа на обеих алгебрах (`'blunt'`, `'norm'`, `'tau_vN'`; веса
$1/(Ck_c)$ и $1/(C'k'_d)$ следа `tau_vN` учитывают ещё и числа секторов). Проверим предложение 6.8 для
случайного отображения $M_2\oplus M_3\to M_1\oplus M_2\oplus M_2$ и трёх следов:

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

> **Физический смысл** (физическая интерпретация, а не теорема этого курса). Объект, не зависящий от следа, — это гейзенберговское отображение
> $\Psi:\mathcal N\to\mathcal M$ на наблюдаемых: состояние $\omega$ до процесса становится
> $\omega\circ\Psi$ после него, и $\Psi$ является CP и унитальным. Плотности зависят от следа
> (теорема 3.2), и плотности преобразуются сопряжённым к $\Psi$ относительно выбранных следов (абзац после
> следствия 6.9). Если алгебры действуют с кратностями, физический след имеет веса $w_c=m_c$ (раздел 2.3).
> В примере сектор $0$ читается как «нет частиц», а сектор $1$ — как «одна частица на двух узлах»; тогда
> отображение удаляет частицу. Оно сохраняет вероятность для следа `Tr_blunt`, который считает состояния,
> но не для `Tr_norm`, который взвешивает секторы одинаково.

## 6.4 Декогеренция

Стандартные каналы библиотеки действуют внутри каждого блока ($p,\gamma\in[0,1]$):

| канал | $\Phi(\rho)_c$ | операторы Крауса в блоке $c$ |
|---|---|---|
| `dephasing_channel` | $(1-p)\rho_c+p\operatorname{diag}\rho_c$ | $\sqrt{1-p}\,1$, $\sqrt p\,E_{ii}$ |
| `depolarizing_channel` | $(1-p)\rho_c+p\operatorname{Tr}(\rho_c)\,1/k_c$ | $\sqrt{1-p}\,1$, $\sqrt{p/k_c}\,E_{ij}$ |
| `center_expectation` | $\operatorname{Tr}(\rho_c)\,1/k_c$ | $\sqrt{1/k_c}\,E_{ij}$ |
| `amplitude_damping_channel` | уровни $i>0$ распадаются в $0$ | $E_{00}+\sqrt{1-\gamma}\sum_{i>0}E_{ii}$, $\sqrt\gamma\,E_{0i}$ |
| `random_mixed_unitary_channel` | $\sum_jp_jU_j\rho U_j^*$ | $\sqrt{p_j}\,U_j$ |

Формы Крауса проверяются прямым вычислением, например
$\sum_{ij}\frac1kE_{ij}\rho E_{ji}=\frac1k\sum_{ij}\rho_{jj}E_{ii}$. Все пять каналов являются TP, и все, кроме
амплитудного затухания с $\gamma>0$, унитальны (предложение 6.4: вычислите $\sum_lK_l^*K_l$ и
$\sum_lK_lK_l^*$; для амплитудного затухания в блоке размера $k$ это $1$ и
$1+(k-1)\gamma E_{00}-\gamma\sum_{i>0}E_{ii}$, а для смешанно-унитарного канала $\sum_jp_j=1$).
Явно, амплитудное затухание переводит $\rho_{00}\mapsto\rho_{00}+\gamma\sum_{i>0}\rho_{ii}$,
$\rho_{ii}\mapsto(1-\gamma)\rho_{ii}$ и $\rho_{ij}\mapsto(1-\gamma)\rho_{ij}$ при положительных
$i\ne j$, а также $\rho_{0j}\mapsto\sqrt{1-\gamma}\,\rho_{0j}$, $\rho_{j0}\mapsto\sqrt{1-\gamma}\,\rho_{j0}$ при $j>0$. Условное ожидание $E$ на центр
(предложение 2.8) — это деполяризация с $p=1$. В отличие от измерения сектора,
$X\mapsto\sum_cP_cXP_c$, которое тождественно на $\mathcal M$ (раздел 3.3), оно сохраняет лишь $p_c$.

**Теорема 6.10 (неподвижные точки унитальных каналов).** Пусть $\Phi(X)=\sum_lK_lXK_l^*$ с $K_l\in\mathcal M$
является TP и унитальным. Тогда $\Phi(X)=X$ тогда и только тогда, когда $X$ коммутирует со всеми $K_l$ и $K_l^*$. В частности,
неподвижные точки образуют подалгебру фон Неймана в $\mathcal M$.

*Доказательство.* Если $X$ коммутирует со всеми $K_l$, то $\Phi(X)=X\sum_lK_lK_l^*=X$. Обратно, $\Phi(X^*)=\Phi(X)^*$,
поэтому самосопряжённые части $(X+X^*)/2$ и $(X-X^*)/2i$ неподвижной точки — тоже неподвижные точки; пусть
$X=X^*$. При $Y_l=K_lX-XK_l$
$$
\sum_lY_lY_l^*=\sum_l\big(K_lX^2K_l^*-K_lXK_l^*X-XK_lXK_l^*+XK_lK_l^*X\big)
=\Phi(X^2)-\Phi(X)X-X\Phi(X)+X\Phi(1)X=\Phi(X^2)-X^2,
$$
где использованы $Y_l^*=XK_l^*-K_l^*X$, $\Phi(X)=X$ и $\Phi(1)=1$. След этого выражения равен $0$, поскольку $\Phi$ является TP. Так как
$\operatorname{Tr}Y_lY_l^*=\|Y_l\|_2^2\ge0$, все $Y_l=0$;
переходя к сопряжённым, получаем, что $X$ коммутирует и с $K_l^*$. Такие $X$ образуют коммутант самосопряжённого множества,
пересечённый с $\mathcal M$: это $*$-алгебра, содержащая $1$, а значит, алгебра фон Неймана
(теорема 2.2). $\square$

Итак, дефазировка ($p>0$) сохраняет в точности диагональ каждого блока, а деполяризация ($p>0$) и $E$
сохраняют в точности центр. Без унитальности теорема неверна: амплитудное затухание оставляет на месте $E_{00}$,
который не коммутирует со скачком $E_{01}$. При $0<\gamma\le1$ и $k\ge2$ его неподвижные точки в блоке
и коммутант его операторов Крауса имеют размерность $1$, но это $\mathbb CE_{00}$ и
$\mathbb C1$ соответственно: по явной формуле выше у неподвижной точки $\rho_{ij}=0$ при $i\ne j$ и
$\rho_{ii}=0$ при $i>0$ (эти элементы умножаются на $\sqrt{1-\gamma}\ne1$ или $1-\gamma\ne1$), и всякое
кратное $E_{00}$ неподвижно; а всё, что коммутирует со всеми $E_{0i}$, $E_{i0}$ ($i>0$), коммутирует
и с их произведениями $E_{ij}=E_{i0}E_{0j}$ и $E_{00}=E_{01}E_{10}$, а значит, с $M_k$, и поэтому кратно
$1$ (глава 2, упражнение 1).

**Теорема 6.11 (унитальные каналы увеличивают энтропию).** Пусть $\Phi:\mathcal M\to\mathcal N$ положительно,
является TP и унитально, а $\rho$ и $\Phi(\rho)$ — состояния с собственными значениями $\lambda_i$ и $\mu_j$.
Тогда $\mu=\Delta\lambda$ для некоторой дважды стохастической матрицы $\Delta$, поэтому $\sum_jf(\mu_j)\ge\sum_if(\lambda_i)$ для
любой вогнутой $f$ на $[0,1]$. В частности, $S(\Phi(\rho))\ge S(\rho)$ и
$\operatorname{Tr}\Phi(\rho)^2\le\operatorname{Tr}\rho^2$.

*Доказательство.* Диагонализуем каждый блок: $\rho=\sum_i\lambda_i|e_i\rangle\langle e_i|$ и
$\Phi(\rho)=\sum_j\mu_j|f_j\rangle\langle f_j|$ с ортонормированными базисами $\bigoplus_c\mathbb C^{k_c}$
и $\bigoplus_d\mathbb C^{k'_d}$, каждый вектор которых лежит в одном слагаемом, так что
$|e_i\rangle\langle e_i|\in\mathcal M$. Тогда $\mu_j=\langle f_j,\Phi(\rho)f_j\rangle=\sum_i\Delta_{ji}\lambda_i$
с $\Delta_{ji}=\langle f_j,\Phi(|e_i\rangle\langle e_i|)f_j\rangle\ge0$,
$\sum_j\Delta_{ji}=\operatorname{Tr}\Phi(|e_i\rangle\langle e_i|)=1$ и
$\sum_i\Delta_{ji}=\langle f_j,\Phi(1)f_j\rangle=1$ (поэтому $\Delta$ квадратная: сложите все элементы). По
конечному неравенству Йенсена $f(\sum_i\Delta_{ji}\lambda_i)\ge\sum_i\Delta_{ji}f(\lambda_i)$ (индукция по
числу слагаемых из определения вогнутости, как в решении упражнения 2 главы 3)
$\sum_jf(\mu_j)\ge\sum_j\sum_i\Delta_{ji}f(\lambda_i)=\sum_if(\lambda_i)$. Возьмите $f(t)=-t\log t$
($f(0)=0$) и $f(t)=-t^2$, которые вогнуты на $[0,1]$ по лемме 3.4: $t\log t$ и $t^2$
непрерывны на $[0,1]$ и имеют на $(0,1)$ вторые производные $1/t>0$ и $2>0$. $\square$

Классическая часть отображения между секторами — это его **матрица переходов**
$T_{dc}(\rho)=\operatorname{Tr}\Phi_{dc}(\rho_c)/\operatorname{Tr}\rho_c$ ($\rho_c\ne0$); метод библиотеки
`transition_matrix()` возвращает $T=T(1)$, т. е. матрицу для максимально смешанного состояния каждого входного сектора.

**Предложение 6.12.** Пусть $\Phi:\mathcal M\to\mathcal N$ положительно и является TP, $\rho$ — состояние с
$p_c=\operatorname{Tr}\rho_c>0$ и $q_d=\operatorname{Tr}\Phi(\rho)_d$. Тогда $T(\rho)$ стохастична по
столбцам и $q=T(\rho)p$. С условными ожиданиями $E_{\mathcal M}$, $E_{\mathcal N}$ на
центры (предложение 2.8)
$$
E_{\mathcal N}\Phi E_{\mathcal M}(\rho)=\bigoplus_d\frac{(Tp)_d}{k'_d}\,1_{k'_d},\qquad T=T(1).
$$

*Доказательство.* $T_{dc}\ge0$, так как $\Phi_{dc}$ положительно, и
$\sum_dT_{dc}=\operatorname{Tr}\Phi(\iota_c\rho_c)/\operatorname{Tr}\rho_c=1$. Далее,
$q_d=\sum_c\operatorname{Tr}\Phi_{dc}(\rho_c)=\sum_cT_{dc}(\rho)p_c$. Наконец,
$E_{\mathcal M}\rho=\bigoplus_cp_c1_{k_c}/k_c$, поэтому след
$\Phi(E_{\mathcal M}\rho)_d=\sum_cp_c\Phi_{dc}(1_{k_c}/k_c)$ равен $\sum_cT_{dc}p_c$, а $E_{\mathcal N}$
заменяет блок $d$ его следом, умноженным на $1/k'_d$. $\square$

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

**Когерентности кубита.** На $M_2$ пусть $D_p$ — дефазировка, а $A_\gamma$ — амплитудное затухание. По таблице
и явной формуле выше $D_p$ умножает $\rho_{01}$ и $\rho_{10}$ на $1-p$ и оставляет на месте
диагональ, а $A_\gamma$ умножает $\rho_{01}$ и $\rho_{10}$ на $\sqrt{1-\gamma}$ и переводит
диагональ в $(\rho_{00}+\gamma\rho_{11},(1-\gamma)\rho_{11})$ независимо от когерентностей. Следовательно,
$D_pA_\gamma=A_\gamma D_p$: оба умножают когерентности на $(1-p)\sqrt{1-\gamma}$ и действуют на диагонали
как $A_\gamma$. В частности, при $T_1,T_\varphi>0$, $t\ge0$, $1-p=e^{-t/T_\varphi}$ и $1-\gamma=e^{-t/T_1}$
композиция умножает когерентности на $e^{-t/T_2}$, где
$$
\frac1{T_2}=\frac1{2T_1}+\frac1{T_\varphi},
$$
и эти семейства композируются аддитивно по $t$: $D_{p(t)}D_{p(s)}=D_{p(t+s)}$ и
$A_{\gamma(t)}A_{\gamma(s)}=A_{\gamma(t+s)}$, поскольку множители $e^{-t/T}$ перемножаются, а след сохраняется.

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

> **Физический смысл** (физическая интерпретация, а не теорема этого курса). Для кубита
> дефазировка моделирует чистую дефазировку со временем $T_\varphi$, а амплитудное затухание — распад
> возбуждённого уровня со временем $T_1$; тогда вычисление выше даёт время когерентности $T_2$ кубита,
> подверженного обоим процессам. Теорема 6.10 — алгебраическая
> форма информации, свободной от декогеренции: выживает то, что коммутирует со всеми операторами скачков (глава 10).
> Предложение 6.12 утверждает, что, если смотреть через центр, отображение между секторами — это классическая марковская
> матрица на зарядах.

## 6.5 Случайные каналы

**Изометрия, случайная по Хаару**, $V:\mathbb C^k\to\mathbb C^N$ — это, как и в библиотеке (`random_isometry`),
множитель $Q$ в $G=QR$, где $G$ — матрица Жинибра размера $N\times k$ (независимые элементы с
плотностью $\pi^{-1}e^{-|z|^2}$), $Q^*Q=1$, а $R$ — верхнетреугольная с положительной диагональю (на
CPU библиотека берёт разложение LAPACK $G=Q'R'$ и умножает $i$-й столбец $Q'$ на
$R'_{ii}/|R'_{ii}|$, что, как в теореме 5.3, даёт это разложение, если все $R'_{ii}\ne0$, т. е.
если $G$ имеет ранг $k$, а это выполняется почти наверное по доказательству леммы 6.13; для больших пакетов на
GPU она использует процесс Грама–Шмидта, который даёт его непосредственно).
`random_channel(alg, r)` независимо в каждом блоке генерирует такие $V_c:\mathbb C^{k_c}\to\mathbb C^r\otimes\mathbb C^{k_c}$
и использует операторы Крауса из теоремы 6.5.

**Лемма 6.13.** $V$ распределена как первые $k$ столбцов унитарной матрицы, случайной по Хаару в $U(N)$. В
частности, для любой унитарной $U\in M_N$ матрица $UV$ имеет то же распределение, что и $V$.

*Доказательство.* Пусть $Z$ — матрица Жинибра размера $N\times N$, первые $k$ столбцов которой образуют $G$, и $E\in M_{N\times k}$ —
первые $k$ столбцов $1_N$. По теореме 5.3 $Z$ почти наверное обратима и $Z=Q_ZR_Z$, где
$Q_Z$ распределена по Хаару. Так как $R_Z$ верхнетреугольна, $G=ZE=(Q_ZE)(E^*R_ZE)$ — разложение
требуемого вида, и оно единственно: если $G=QR=Q'R'$, то $R^*R=G^*G=R'^*R'$, поэтому $R'R^{-1}$
унитарна, верхнетреугольна и имеет положительную диагональ, а значит, равна $1$ (доказательство леммы 5.2). Итак, $V=Q_ZE$,
и $UV=(UQ_Z)E$ распределена как $V$, поскольку $UQ_Z$ распределена как $Q_Z$. $\square$

**Лемма 6.14.** Пусть $u$ — случайный единичный вектор в $\mathbb C^N$, распределение которого инвариантно относительно всех
унитарных преобразований. Тогда $\mathbb E[u_a\bar u_b]=\delta_{ab}/N$ и
$$
\mathbb E\big[u_a\bar u_b\,u_c\bar u_d\big]=\frac{\delta_{ab}\delta_{cd}+\delta_{ad}\delta_{cb}}{N(N+1)} .
$$

*Доказательство.* Инвариантность относительно $\operatorname{diag}(e^{i\theta_1},\dots,e^{i\theta_N})$ умножает
моменты на $e^{i(\theta_a-\theta_b)}$ и $e^{i(\theta_a-\theta_b+\theta_c-\theta_d)}$, так что они обращаются в нуль,
если только не $a=b$, соответственно не $\{a,c\}=\{b,d\}$ как мультимножества. В силу инвариантности относительно матриц перестановок
$\mathbb E|u_a|^2=1/N$, и остаются лишь $\alpha=\mathbb E|u_1|^4$ и $\beta=\mathbb E|u_1|^2|u_2|^2$.
Во-первых, $\sum_c|u_c|^2=1$ даёт $\alpha+(N-1)\beta=\mathbb E|u_1|^2=1/N$. Во-вторых, $(u_1+u_2)/\sqrt2$ —
первая координата $Uu$ для некоторой унитарной $U$, поэтому $\mathbb E|u_1+u_2|^4=4\alpha$. При
$s=|u_1|^2+|u_2|^2$ и $z=u_1\bar u_2$ имеем $|u_1+u_2|^4=s^2+2s(z+\bar z)+(z+\bar z)^2$;
инвариантность относительно $u_1\mapsto-u_1$ и $u_1\mapsto iu_1$ даёт $\mathbb E[s(z+\bar z)]=0$ и
$\mathbb Ez^2=0$, так что $4\alpha=(2\alpha+2\beta)+2\beta$, т. е. $\alpha=2\beta$. Следовательно,
$\beta=1/(N(N+1))$ и $\alpha=2\beta$, что и даёт формулу во всех случаях. (Здесь использовалось $N\ge2$; при
$N=1$ $|u_1|=1$ и обе части равны $1$.) $\square$

**Теорема 6.15 (выходная чистота случайного канала).** Пусть $\psi$ — чистое состояние в секторе размера
$k$, а $\Phi$ сгенерирован функцией `random_channel` с рангом Крауса $r$. Тогда
$$
\mathbb E\operatorname{Tr}\Phi(\psi)^2=\frac{k+r}{kr+1},\qquad\mathbb E\,\Phi(\psi)=\frac1k1_k .
$$

*Доказательство.* В секторе $\psi$ имеем $\Phi(\psi)=\operatorname{Tr}_E|u\rangle\langle u|$ с $u=V\psi$ —
единичным вектором, распределение которого унитарно инвариантно по лемме 6.13. Занумеруем базис $\mathbb C^r\otimes\mathbb C^k$
парами $(e,x)$. Тогда $\Phi(\psi)_{xy}=\sum_eu_{ex}\bar u_{ey}$, среднее этой величины равно $r\delta_{xy}/(kr)$, и
$\operatorname{Tr}\Phi(\psi)^2=\sum_{x,y,e,f}u_{ex}\bar u_{ey}u_{fy}\bar u_{fx}$. По лемме 6.14 с
$N=kr$ первое слагаемое требует $x=y$ ($kr^2$ наборов индексов), а второе — $e=f$ ($k^2r$ наборов), так что
среднее равно $(kr^2+k^2r)/(kr(kr+1))$. $\square$

При $r=1$ канал унитарен и чистота равна $1$; при $r\to\infty$ она стремится к $1/k$ —
чистоте максимально смешанного состояния. Доказательство использовало лишь унитарную инвариантность распределения $u$,
поэтому оно даёт и среднюю чистоту Лубкина для редуцированного состояния случайного чистого состояния $u=g/\|g\|$,
$g$ гауссов (его распределение унитарно инвариантно по теореме A.20(b)), которую глава 7 (теорема 7.11)
доказывает с помощью гауссовых интегралов; числа $\alpha,\beta$ из леммы 6.14 — это числа из следствия 5.11
для элементов унитарной матрицы, случайной по Хаару, а следствие A.21 вычисляет их заново для $u=g/\|g\|$.

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

> **Физический смысл** (физическая интерпретация, а не теорема этого курса). Случайный канал
> моделирует типичный неструктурированный шум от окружения
> размерности $r$. Уже $r=4$ снижает чистоту чистого состояния кутрита в среднем до $7/13\approx0.54$
> (теорема 6.15), тогда как без шума она равна $1$, а минимум равен $1/3$ (по неравенству Коши–Буняковского
> $\operatorname{Tr}\rho^2\ge(\operatorname{Tr}\rho)^2/k$). В системе с суперотбором шум действует в каждом секторе
> отдельно, и существенна размерность сектора, а не всего пространства.

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

## Упражнения

1. (Вернер–Холево) Покажите, что $\Phi(\rho)=\big(\operatorname{Tr}(\rho)1-\rho^T\big)/(k-1)$ на $M_k$,
   $k\ge2$, — унитальный канал с рангом Крауса $k(k-1)/2$, хотя он построен из
   транспонирования. Проверьте это с помощью `Channel.from_superoperator`.
2. Покажите, что при $0<p\le1$ канал дефазировки на $M_k$ имеет ранг Крауса $k$, на единицу меньше
   числа операторов Крауса, которые хранит библиотека.
3. Пусть $X$ — циклический сдвиг на $\mathbb C^k$ и $Z=\operatorname{diag}(1,\omega,\dots,\omega^{k-1})$,
   $\omega=e^{2\pi i/k}$. Покажите, что $\frac1{k^2}\sum_{a,b=0}^{k-1}X^aZ^b\rho Z^{-b}X^{-a}=\operatorname{Tr}(\rho)1/k$:
   условное ожидание на центр $M_k$ является смешанно-унитарным каналом.
4. Покажите, что $(\Psi\circ\Phi)^\sharp=\Phi^\sharp\circ\Psi^\sharp$, если на промежуточной алгебре используется
   один и тот же след, и проверьте это для трёх следов библиотеки.

<details>
<summary>Решения</summary>

1. Если $F$ — оператор перестановки, то $J(\theta)=\sum_{ij}E_{ij}\otimes E_{ji}=F$ и $J(\rho\mapsto\operatorname{Tr}(\rho)1)=1\otimes1$,
   так что $J(\Phi)=(1-F)/(k-1)=2P_-/(k-1)$, где $P_-=(1-F)/2$ — проектор на
   антисимметричные векторы ($F=F^*$, $F^2=1$) ранга $\operatorname{Tr}P_-=(k^2-\operatorname{Tr}F)/2$, где
   $\operatorname{Tr}F=\sum_{ij}\langle e_i\otimes e_j,e_j\otimes e_i\rangle=k$. Значит, $J\ge0$ и имеет ранг $k(k-1)/2$ (теорема 6.2,
   следствие 6.3), а $\operatorname{Tr}\Phi(\rho)=\operatorname{Tr}\rho$ и $\Phi(1)=1$ проверяются непосредственно.

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

2. $J=(1-p)|\Omega\rangle\langle\Omega|+p\sum_iE_{ii}\otimes E_{ii}$ обращается в нуль вне линейной оболочки
   векторов $e_i\otimes e_i$, а на этом $k$-мерном пространстве является матрицей $(1-p)\mathbf 1\mathbf 1^T+p1$
   ($\mathbf 1$ — вектор из единиц) с собственными значениями $p$ и $(1-p)k+p$, которые все положительны. Поэтому
   $\operatorname{rank}J=k$, тогда как библиотека хранит $1+k_{\max}$ операторов.

```python
ranks = [int((torch.linalg.eigvalsh(J_c[0]) > 1e-10).sum()) for J_c in zoo['dephasing'].choi()]
print(ranks, zoo['dephasing'].kraus_rank)                 # [2, 3] versus 4 stored operators
assert ranks == alg.k_factors and zoo['dephasing'].kraus_rank == 1 + max(alg.k_factors)
```

3. $\operatorname{Tr}((X^aZ^b)^*X^{a'}Z^{b'})=k\,\delta_{aa'}\delta_{bb'}$: у $Z^{-b}X^{a'-a}Z^{b'}$ нулевая
   диагональ, если $a\ne a'$, и $\operatorname{Tr}Z^{b'-b}=\sum_j\omega^{j(b'-b)}=k\delta_{bb'}$. Значит,
   $B_{ab}=X^aZ^b/\sqrt k$ образуют ортонормированный базис $M_k$ относительно $\operatorname{Tr}(A^*B)$, а любые
   два ортонормированных базиса связаны унитарным преобразованием; по теореме 6.7 сумма $\sum_mB_m\rho B_m^*$
   одна и та же для всех ортонормированных базисов. Для базиса $E_{ab}$ она равна
   $\sum_{ab}E_{ab}\rho E_{ba}=\operatorname{Tr}(\rho)1$. Остаётся разделить на $k$.

```python
Xs = torch.roll(torch.eye(k, dtype=dt), 1, dims=0)
Zs = torch.diag(torch.exp(2j * math.pi * torch.arange(k, dtype=torch.float64) / k)).to(dt)
W = [torch.linalg.matrix_power(Xs, a) @ torch.linalg.matrix_power(Zs, b) / k for a in range(k) for b in range(k)]
twirl = channels.kraus_channel([M3.operator(Wab[None, None]) for Wab in W])
r3op = states.random_density_matrix(M3)
assert close(twirl(r3op), channels.center_expectation(M3)(r3op))
```

4. $\operatorname{tr}(\Psi(\Phi(x))z)=\operatorname{tr}(\Phi(x)\Psi^\sharp(z))=\operatorname{tr}(x\,\Phi^\sharp(\Psi^\sharp(z)))$
   для всех $x$ и $z$; единственность сопряжённого даёт утверждение.

```python
Psi2 = channels.random_inter_sector_channel(alg_out, alg, kraus_rank=2)
for t in trace:
    assert close((Psi2 @ Psi).adjoint(trace=t)(x), Psi.adjoint(trace=t)(Psi2.adjoint(trace=t)(x)))
```

</details>

**Далее:** в [главе 7](07_tensor_products.md) рассматриваются тензорные произведения и запутанность: частичный след —
канал между разными алгебрами, сопряжённым к которому является вложение $A\mapsto A\otimes1$, — формула
Пейджа и слияние зарядов.
