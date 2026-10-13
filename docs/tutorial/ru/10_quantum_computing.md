[English](../10_quantum_computing.md) · [Español](../es/10_quantum_computing.md) · [Français](../fr/10_quantum_computing.md) · [Deutsch](../de/10_quantum_computing.md) · [中文](../zh/10_quantum_computing.md) · [日本語](../ja/10_quantum_computing.md) · **Русский**

# 10. Квантовые вычисления

Шум на квантовом компьютере — это канал, а канал порождает алгебру: $*$-алгебру своих операторов
Крауса. Её структура (глава 2) показывает, какая информация переживает шум нетронутой; условие на
произведения операторов Крауса, в котором участвует другая алгебра, показывает, какую информацию
может вернуть восстановление. Последняя часть посвящена алгоритмам, сохраняющим заряд: их состояния
живут в алгебре зарядовых секторов, а преобладающий шум переводит их из сектора в сектор.

**Вы узнаете**, почему у унитального канала и сопряжённого к нему одни и те же неподвижные точки
(коммутант шума), теорему Нилла–Лафламма с полным доказательством, как защитить только метку
сектора и как XY-миксер квантового анзаца чередующихся операторов, амплитудное затухание и
постселекция выглядят в алгебре $\bigoplus_wM_{\binom nw}$ секторов с фиксированным весом Хэмминга.

Всюду далее **канал** на $M_D(\mathbb C)$ — это отображение $\Phi(\rho)=\sum_{i=1}^rE_i\rho E_i^*$ с
**операторами Крауса**, удовлетворяющими $\sum_iE_i^*E_i=1$ (сохранение следа); он **унитален**, если,
кроме того, $\sum_iE_iE_i^*=1$. Его **сопряжённое** $\Phi^*(T)=\sum_iE_i^*TE_i$ удовлетворяет
$\operatorname{Tr}(\Phi(\rho)T)=\operatorname{Tr}(\rho\,\Phi^*(T))$ (раздел 6.2). (По теореме 6.2 всякое вполне
положительное отображение имеет такой вид; здесь каналы *задаются* операторами Крауса.)
**Алгебра взаимодействия** $\mathcal A$ канала $\Phi$ — это $*$-алгебра, порождённая $1$ и $E_i$, т. е. линейная
оболочка всех произведений $E_i$ и $E_i^*$; матрица коммутирует с $\mathcal A$ тогда и только тогда, когда
она коммутирует со всеми $E_i$ и $E_i^*$. Состояние $\rho$ **сосредоточено** в подпространстве $K$ с
проектором $P_K$, если $\rho=P_K\rho P_K$. Для $n$ кубитов $D=2^n$; $X,Y,Z$ — матрицы Паули,
$Z=\operatorname{diag}(1,-1)$, $X_j$ — это $X$ на кубите $j$ (кубиты нумеруются $1,\dots,n$; в разделе 10.4 мы
следуем библиотеке, узлы $0,\dots,n-1$), $|b\rangle$, $b\in\{0,1\}^n$, — вычислительный базис, $|b|$ —
**вес Хэмминга**, а $P_w$ — проектор на $H_w=\operatorname{span}\{|b\rangle:|b|=w\}$.

## 10.1 Алгебры шума, бесшумные подсистемы и неподвижные точки

По теореме 2.6 существует унитарная $V$, для которой $V^*\mathcal AV=\bigoplus_cM_{k_c}\otimes1_{m_c}$ и
$V^*\mathcal A'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$: каждый оператор Крауса действует только на первый
тензорный сомножитель каждого слагаемого.

**Предложение 10.1 (бесшумные подсистемы; подпространства, свободные от декогеренции).** (a) Запишем
$V^*E_iV=\bigoplus_cE_{i,c}\otimes1_{m_c}$ и пусть $\rho=V(\sigma\otimes\tau)V^*$ помещено в слагаемое
$c$, где $\sigma$ и $\tau$ — состояния на $\mathbb C^{k_c}$ и $\mathbb C^{m_c}$. Тогда
$\Phi(\rho)=V(\Phi_c(\sigma)\otimes\tau)V^*$ с каналом $\Phi_c(\sigma)=\sum_iE_{i,c}\sigma E_{i,c}^*$ на
$M_{k_c}$: состояние $\tau$ не затрагивается. Кроме того, $\Phi^*(T)=T$ для всякого $T\in\mathcal A'$.
(b) Если $K$ — подпространство, для которого $E_iP_K=e_iP_K$, $e_i\in\mathbb C$, при всех $i$, то
$\Phi(\rho)=\rho$ для всякого состояния $\rho$, сосредоточенного в $K$.

*Доказательство.* (a) $\Phi(\rho)=V\sum_i(E_{i,c}\sigma E_{i,c}^*\otimes\tau)V^*$, и $\Phi_c$ сохраняет след,
поскольку слагаемое $c$ в $V^*(\sum_iE_i^*E_i)V=1$ равно $(\sum_iE_{i,c}^*E_{i,c})\otimes1_{m_c}$. Для
$T\in\mathcal A'$ имеем $\Phi^*(T)=T\sum_iE_i^*E_i=T$. (b) $P_KE_i^*=(E_iP_K)^*=\bar e_iP_K$, поэтому
$E_i\rho E_i^*=E_iP_K\rho P_KE_i^*=|e_i|^2\rho$, и $\sum_i|e_i|^2P_K=P_K\sum_iE_i^*E_iP_K=P_K$. $\square$

Слагаемое с $k_c=1$ — это случай (b): **подпространство, свободное от декогеренции** (Занарди–Расетти 1997,
Лидар–Чуан–Уэйли 1998). При $k_c>1$ сомножитель $\mathbb C^{m_c}$ — **бесшумная подсистема**
(Нилл–Лафламм–Виола 2000): шум действует, но только на другой сомножитель. Предложение 10.1 верно
для любого канала; для унитальных каналов сохраняющиеся наблюдаемые — это *в точности* наблюдаемые из $\mathcal A'$:

**Теорема 10.2 (неподвижные точки унитальных каналов).** Пусть канал $\Phi$ унитален. Тогда
$$
\operatorname{Fix}(\Phi)=\{T\in M_D(\mathbb C):\Phi(T)=T\}=\mathcal A'=\operatorname{Fix}(\Phi^*).
$$
В частности, $\operatorname{Fix}(\Phi)$ — алгебра фон Неймана, $*$-изоморфная $\bigoplus_cM_{m_c}$.

*Доказательство.* Первое равенство — это теорема 6.10 для алгебры $M_D(\mathbb C)$ (один блок): $T$ —
неподвижная точка тогда и только тогда, когда $T$ коммутирует со всеми $E_i$ и $E_i^*$, т. е. когда $T\in\mathcal A'$.
Сопряжённое отображение $\Phi^*$ — снова унитальный канал: его операторы Крауса $E_i^*$ удовлетворяют
$\sum_iE_iE_i^*=1$ (сохранение следа) и $\sum_iE_i^*E_i=1$ (унитальность) и порождают ту же $*$-алгебру
$\mathcal A$. Поэтому теорема 6.10, применённая к $\Phi^*$, даёт $\operatorname{Fix}(\Phi^*)=\mathcal A'$. Последнее
утверждение — это теорема 2.6. $\square$

Для любого канала $\mathcal A'\subseteq\operatorname{Fix}(\Phi^*)$ (предложение 10.1(a)), но без
унитальности теорема неверна: для амплитудного затухания $\operatorname{Fix}(\Phi)=\mathbb C|0\rangle\langle0|$, тогда как
$\mathcal A'=\operatorname{Fix}(\Phi^*)=\mathbb C1$ (раздел 6.4 и упражнение 1). (Крибс 2003,
Ариас–Геондя–Гаддер 2002; в бесконечномерном случае нужны дополнительные предположения.)

**Следствие 10.3 (коллективная дефазировка).** Пусть $U=e^{-i\theta Z_{\mathrm{tot}}}$,
$Z_{\mathrm{tot}}=\sum_jZ_j$, причём $e^{2i\theta d}\neq1$ при $d=1,\dots,n$, и
$\Phi(\rho)=\tfrac12(U\rho U^*+U^*\rho U)$. Тогда $\mathcal A=\operatorname{span}\{P_0,\dots,P_n\}$, каждое $H_w$ —
подпространство, свободное от декогеренции, и $\operatorname{Fix}(\Phi)=\bigoplus_wB(H_w)\cong\bigoplus_wM_{\binom nw}$
размерности $\sum_w\binom nw^2$.

*Доказательство.* $Z_{\mathrm{tot}}|b\rangle=(n-2|b|)|b\rangle$, поэтому $U=\sum_we^{-i\theta(n-2w)}P_w$ с попарно
различными собственными значениями (их отношения равны $e^{2i\theta(w-w')}$, $0<|w-w'|\le n$). Многочлены
Лагранжа из доказательства теоремы 1.5(a) дают $P_w=p_w(U)\in\mathcal A$; обратно, $U,U^*$ лежат в алгебре
$\operatorname{span}\{P_w\}$. Значит, $\mathcal A=\operatorname{span}\{P_w\}$ ($k_c=1$, $m_c=\binom nw$), и к каждому
$H_w$ применимо предложение 10.1(b). Если $T$ коммутирует со всеми $P_w$, то $P_vTP_{v'}=TP_vP_{v'}=0$ при
$v\neq v'$, т. е. $T$ блочно-диагональна, и обратно; остальное даёт теорема 10.2. $\square$

Коммутанты — это ядра отображений $T\mapsto(AT-TA)_A$ (глава 2), неподвижные точки — ядро
$\sum_iE_i\otimes\bar E_i-1$ (матрицы $\Phi-\mathrm{id}$ на матрицах, векторизованных по строкам). При $n=4$
три пространства из следствия имеют размерности $5$, $70=1+16+36+16+1$ и $70$:

```python
import itertools
import math
import torch

torch.manual_seed(0)
dt = torch.complex128
I2, X = torch.eye(2, dtype=dt), torch.tensor([[0, 1], [1, 0]], dtype=dt)
Y, Z = torch.tensor([[0, -1j], [1j, 0]], dtype=dt), torch.tensor([[1, 0], [0, -1]], dtype=dt)

def on(op, j, n):                             # op on qubit j of n (j = 0 is the leftmost factor)
    out = torch.ones(1, 1, dtype=dt)
    for i in range(n):
        out = torch.kron(out, op if i == j else I2)
    return out

def null_space(L, tol=1e-9):                  # orthonormal columns spanning {x : L x = 0}
    w, V = torch.linalg.eigh(L.mH @ L)
    return V[:, w < tol]

def commutant_basis(ops):
    d, ops = ops[0].shape[0], [A.resolve_conj().contiguous() for A in ops]
    eye = torch.eye(d, dtype=dt)
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    return [v.reshape(d, d) for v in null_space(L).T]

def fixed_points(kraus):                      # vec(E X E^*) = (E (x) conj E) vec(X)
    d = kraus[0].shape[0]
    S = sum(torch.kron(E.resolve_conj(), E.conj().resolve_conj()) for E in kraus)
    return [v.reshape(d, d) for v in null_space(S - torch.eye(d * d, dtype=dt)).T]

n, theta = 4, 0.3
U = torch.linalg.matrix_exp(-1j * theta * sum(on(Z, j, n) for j in range(n)))
kraus = [U / math.sqrt(2), U.mH / math.sqrt(2)]                  # unital: U U^* = U^* U = 1
comm = commutant_basis(kraus + [E.mH for E in kraus])             # A'
A_int, fix = commutant_basis(comm), fixed_points(kraus)           # A = A'' (Theorem 2.2), Fix(Phi)
print(len(A_int), len(comm), len(fix))
assert (len(A_int), len(comm), len(fix)) == (5, 70, 70)
span = torch.linalg.svdvals(torch.stack([F.reshape(-1) for F in comm + fix]))
assert int((span > 1e-9).sum()) == 70                             # Fix(Phi) = A' (Theorem 10.2)
```

$\bigoplus_wM_{\binom nw}$ — это алгебра `SpinChain` (глава 9), секторы которой нумеруются числом спинов,
направленных вверх, т. е. весом Хэмминга. Там операторы Крауса центральны,
$U=\sum_we^{-i\theta(n-2w)}P_w$, и канал тождествен: библиотека хранит бесшумную часть.

```python
from torch_vn_algebra import SpinChain, TypeIAlgebra, DensityMatrix, channels, states, dynamics

M = SpinChain(n).algebra                                          # (+)_w M_binom(4, w)
print(M.k_factors)
assert sum(k * k for k in M.k_factors) == len(fix)
phases = [torch.exp(torch.tensor(-1j * theta * (n - 2 * w), dtype=dt)) for w in range(n + 1)]
Phi_M = channels.kraus_channel([M.central([z / math.sqrt(2) for z in phases]),
                                M.central([z.conj() / math.sqrt(2) for z in phases])])
rho = states.random_density_matrix(M, batch_size=3)
assert (Phi_M(rho).matrix - rho.matrix).abs().max() < 1e-12      # no state of M is disturbed
```

> **Физический смысл.** Флуктуирующее поле, одинаковое на всех кубитах, сдвигает фазу $|b\rangle$
> пропорционально $n-2|b|$. Состояния с одним и тем же весом Хэмминга приобретают общую фазу, поэтому $|01\rangle,|10\rangle$
> кодируют кубит, невосприимчивый к этому шуму (ионы в ловушках: Kielpinski et al., *Science* 291, 2001), тогда как
> когерентности между разными весами затухают: шум навязывает правило суперотбора (раздел 3.3).

**Коллективная декогеренция.** Если поле флуктуирует ещё и по направлению, операторы Крауса — это
повороты $e^{-i\theta\,\mathbf n\cdot\mathbf S}$, $S_\alpha=\tfrac12\sum_j\sigma_{\alpha,j}$, $\sigma\in\{X,Y,Z\}$,
лежащие в $*$-алгебре $\mathcal B$, порождённой $1,S_x,S_y,S_z$ (теорема 1.5(a)).

**Предложение 10.4 (синглетное подпространство).** Пусть $K=\ker S_x\cap\ker S_y\cap\ker S_z$. Всякий канал,
операторы Крауса которого лежат в $\mathcal B$, оставляет неподвижным всякое состояние, сосредоточенное в $K$.

*Доказательство.* $\mathcal B$ натянута на слова из $S_x,S_y,S_z$ (и $1$); слова длины $\ge1$ аннулируют
$K$, поэтому каждый элемент $\mathcal B$ действует на $K$ как скаляр. Остаётся применить предложение 10.1(b). $\square$

При $n=4$ код находит $\dim K=2$, $\dim\mathcal B'=14$, $\dim\mathcal B=35$, как и должно быть для $\mathcal B\cong M_5\oplus M_3\oplus M_1$
с кратностями $1,3,2$ (полный спин $2,1,0$), что выводится в главе 4 (готовится): четыре
кубита несут свободный от декогеренции кубит ($k=1$, $m=2$) и бесшумный кутрит ($k=3$, $m=3$).

```python
S = [sum(on(P, j, n) for j in range(n)) / 2 for P in (X, Y, Z)]
comm_S = commutant_basis(S)
print(len(comm_S), len(commutant_basis(comm_S)))
assert len(comm_S) == 14 and len(commutant_basis(comm_S)) == 35
K = null_space(torch.cat(S))                                      # columns: a basis of K
assert K.shape[1] == 2
axes = torch.nn.functional.normalize(torch.randn(3, 3, dtype=torch.float64), dim=1)
rotations = [torch.linalg.matrix_exp(-1j * t * sum(c * Sa for c, Sa in zip(ax, S))) / math.sqrt(3)
             for t, ax in zip(torch.rand(3, dtype=torch.float64) * math.pi, axes)]
apply = lambda kr, r: sum(E @ r @ E.mH for E in kr)
for G in (K @ torch.randn(2, 2, dtype=dt), torch.randn(16, 16, dtype=dt)):  # supported in K / generic
    r = G @ G.mH / (G @ G.mH).trace()
    assert ((apply(rotations, r) - r).abs().max() < 1e-12) == (G.shape[1] == 2)
```

## 10.2 Коррекция ошибок: теорема Нилла–Лафламма

**Код** — это подпространство $\mathcal C\subseteq\mathbb C^D$ с проектором $P$. Он **исправляет** $\Phi$, если
существует канал $\mathcal R$ (**восстановление**), для которого $\mathcal R(\Phi(\rho))=\rho$ при всяком состоянии $\rho$,
сосредоточенном в $\mathcal C$.

**Теорема 10.5 (Нилл–Лафламм).** $\mathcal C$ исправляет $\Phi$ тогда и только тогда, когда существуют числа $\lambda_{ij}$,
для которых
$$
PE_i^*E_jP=\lambda_{ij}P\qquad\text{for all }i,j.
$$

*Доказательство.* *Необходимость.* Пусть $R_m$ — операторы Крауса канала $\mathcal R$; тогда у $\mathcal R\circ\Phi$
операторы Крауса $R_mE_i$. Для единичного вектора $\psi\in\mathcal C$ имеем
$\sum_{m,i}R_mE_i|\psi\rangle\langle\psi|E_i^*R_m^*=|\psi\rangle\langle\psi|$, и при $\phi\perp\psi$ среднее в
состоянии $\phi$ даёт $\sum_{m,i}|\langle\phi,R_mE_i\psi\rangle|^2=0$, так что $R_mE_i\psi\in\mathbb C\psi$.
Следовательно, $R_mE_iP=c_{mi}P$: линейное отображение $T$ с $Tx\in\mathbb Cx$ для всех $x$ кратно $1$ (если
$Tx=ax$, $Ty=by$, $T(x+y)=c(x+y)$ для линейно независимых $x,y$, то $a=c=b$). Так как $\sum_mR_m^*R_m=1$,
$PE_i^*E_jP=\sum_m(R_mE_iP)^*(R_mE_jP)=\sum_m\bar c_{mi}c_{mj}P$.

*Достаточность.* Матрица $\Lambda=(\lambda_{ij})$ эрмитова, так как $\lambda_{ji}P=(PE_i^*E_jP)^*=\bar\lambda_{ij}P$, и
положительна, так как $\sum_{ij}\bar c_ic_j\lambda_{ij}P=PF^*FP\ge0$ для $F=\sum_jc_jE_j$. Запишем $\Lambda=UDU^*$,
где $U$ унитарна, $D=\operatorname{diag}(d_k)$, $d_k\ge0$ (теорема 1.3), и положим $F_k=\sum_iu_{ik}E_i$. Тогда
$\sum_kF_k\rho F_k^*=\sum_{ij}(UU^*)_{ij}E_i\rho E_j^*=\Phi(\rho)$ и $PF_k^*F_lP=(U^*\Lambda U)_{kl}P=d_k\delta_{kl}P$.
Если $d_k=0$, то $\|F_kPx\|^2=0$ для всех $x$, поэтому $F_kP=0$. Если $d_k>0$, положим $V_k=F_kP/\sqrt{d_k}$; тогда
$V_k^*V_l=\delta_{kl}P$, поэтому $V_kP=V_k$, операторы $Q_k=V_kV_k^*$ — попарно ортогональные проекторы и
$Q=\sum_kQ_k$ — проектор. Восстановление
$$
\mathcal R(\sigma)=\sum_{k:d_k>0}V_k^*\sigma V_k+(1-Q)\sigma(1-Q)
$$
является каналом, поскольку $\sum_kV_kV_k^*+(1-Q)^2=1$. Для $\rho=P\rho P$ имеем $\Phi(\rho)=\sum_{d_k>0}d_kV_k\rho V_k^*$;
так как $V_l^*V_k=\delta_{kl}P$ и $(1-Q)V_k=(1-Q)Q_kV_k=0$, получаем $\mathcal R(\Phi(\rho))=\sum_kd_k\rho$, причём
$\sum_kd_k=\operatorname{Tr}\Lambda=1$, поскольку $\sum_i\lambda_{ii}P=P\sum_iE_i^*E_iP=P$. $\square$

(Нилл–Лафламм 1997; Нильсен–Чанг, теорема 10.1.) $\mathcal R$ измеряет **синдром** — то, какое из
ортогональных подпространств $Q_k\mathbb C^D$ содержит состояние, — и поворачивает состояние обратно в код.

**Пример: код против инверсий бита.** $\mathcal C=\operatorname{span}\{|000\rangle,|111\rangle\}$ и «не более одной
инверсии бита»: $E_0=\sqrt{p_0}\,1$, $E_j=\sqrt{p_j}\,X_j$. При $i\neq j$ оператор $X_iX_j$ (с $X_0=1$) инвертирует один или два
бита и переводит $\mathcal C$ в $\mathcal C^\perp$; поэтому $\Lambda=\operatorname{diag}(p_0,\dots,p_3)$, $V_k=X_kP$,
$Q=1$ и $\mathcal R(\sigma)=\sum_{k=0}^3PX_k\sigma X_kP$. Если каждый кубит инвертируется независимо с
вероятностью $p$, условие нарушается, и $\mathcal R$ переводит $|0_L\rangle=|000\rangle$ в $|1_L\rangle$, когда
происходит $\ge2$ инверсий, и оставляет его на месте в противном случае: верность (fidelity) состояния $|0_L\rangle$ после
восстановления равна вероятности не более чем одной инверсии, $(1-p)^3+3p(1-p)^2=1-3p^2+2p^3=0.972$ при $p=0.1$.
В библиотеке $M_8$ — алгебра с одним блоком, а `R @ Phi` — композиция $\mathcal R\circ\Phi$:

```python
M8 = TypeIAlgebra([8], [8], precision='double')
op8 = lambda T: M8.operator(T.reshape(1, 1, 8, 8))
state8 = lambda v: DensityMatrix(M8, matrix=torch.outer(v, v.conj()).reshape(1, 1, 8, 8))
ket = lambda bits: torch.eye(8, dtype=dt)[int(bits, 2)]
def kl_matrix(kraus, P, tol=1e-12):
    """(lambda_ij) with P E_i^* E_j P = lambda_ij P, or None if the Knill-Laflamme condition fails."""
    lam = torch.tensor([[(P @ a.mH @ b @ P).trace() / P.trace() for b in kraus] for a in kraus])
    ok = all((P @ a.mH @ b @ P - lam[i, j] * P).abs().max() < tol
             for i, a in enumerate(kraus) for j, b in enumerate(kraus))
    return lam if ok else None
def independent(op, prob):                    # op on each qubit independently with probability prob
    return [math.sqrt(prob ** sum(s) * (1 - prob) ** (3 - sum(s)))
            * torch.linalg.multi_dot([on(op, j, 3) if s[j] else torch.eye(8, dtype=dt) for j in range(3)])
            for s in itertools.product([0, 1], repeat=3)]

P = torch.outer(ket('000'), ket('000')) + torch.outer(ket('111'), ket('111'))
flips = [torch.eye(8, dtype=dt)] + [on(X, j, 3) for j in range(3)]
E_one = [math.sqrt(pk) * Xk for pk, Xk in zip((0.7, 0.1, 0.1, 0.1), flips)]
assert torch.allclose(kl_matrix(E_one, P), torch.diag(torch.tensor([0.7, 0.1, 0.1, 0.1], dtype=dt)))
Phi_one = channels.kraus_channel([op8(E) for E in E_one])
R = channels.kraus_channel([op8(P @ Xk) for Xk in flips])        # syndrome measurement + correction
assert R.is_trace_preserving()
a, b = torch.nn.functional.normalize(torch.randn(2, dtype=dt), dim=0)
rho_L = state8(a * ket('000') + b * ket('111'))
assert abs((R @ Phi_one)(rho_L).fidelity(rho_L).item() - 1) < 1e-12       # Theorem 10.5
p = 0.1                                                           # independent flips: KL fails
E_ind = independent(X, p)
assert kl_matrix(E_ind, P) is None
zero_L = state8(ket('000'))
F = (R @ channels.kraus_channel([op8(E) for E in E_ind]))(zero_L).fidelity(zero_L).item()
print(f"fidelity after recovery {F:.4f}, exact {1 - 3 * p**2 + 2 * p**3:.4f}")
assert abs(F - (1 - 3 * p**2 + 2 * p**3)) < 1e-12
```

## 10.3 Коррекция ошибок для операторных алгебр

Картина Гейзенберга позволяет защищать меньше, чем весь код. Пусть $\mathfrak A$ — алгебра фон
Неймана на $\mathcal C=P\mathbb C^D$ ($*$-алгебра операторов на $\mathcal C$, содержащая $P$), состоящая из наблюдаемых,
которые нужно защитить. $\mathfrak A$ **корректируема** для $\Phi$, если существует канал
$\mathcal R$, для которого
$$
P\,(\mathcal R\circ\Phi)^*(A)\,P=A\qquad\text{for all }A\in\mathfrak A,
$$
т. е. $\operatorname{Tr}(\mathcal R(\Phi(\rho))A)=\operatorname{Tr}(\rho A)$ для всех состояний $\rho$, сосредоточенных в $\mathcal C$,
и всех $A\in\mathfrak A$. При $\mathfrak A=B(\mathcal C)$ это коррекция в смысле раздела 10.2: если
$\sigma=\mathcal R(\Phi(\rho))$ удовлетворяет $P\sigma P=\rho$, то $\operatorname{Tr}(1-P)\sigma(1-P)=0$, поэтому $\sigma^{1/2}(1-P)=0$ и
$\sigma=\rho$. Противоположный крайний случай — коммутативная алгебра: защищается лишь классическая метка.

**Предложение 10.6 (коррекция метки сектора).** Пусть $P=\sum_{c=1}^CP_c$ с попарно ортогональными
проекторами $P_c\neq0$ и $\mathfrak A=\operatorname{span}\{P_1,\dots,P_C\}$. Следующие условия равносильны:
1. существуют $N_1,\dots,N_C\ge0$, $\sum_cN_c=1$, такие, что $\operatorname{Tr}(N_c\Phi(\rho))=1$ для каждого $c$ и
   каждого состояния $\rho$, сосредоточенного в $P_c\mathbb C^D$;
2. $P_cE_i^*E_jP_d=0$ для всех $i,j$ и $c\neq d$, т. е. $PE_i^*E_jP\in\mathfrak A'$ — коммутанту
   $\mathfrak A$ в $B(\mathcal C)$;
3. $\mathfrak A$ корректируема для $\Phi$.

*Доказательство.* $T\in B(\mathcal C)$ коммутирует со всеми $P_c$ тогда и только тогда, когда $P_cTP_d=0$ при $c\neq d$ (как в следствии 10.3).

(2)$\Rightarrow$(1). Подпространства $S_c=\operatorname{span}\{E_iP_cx\}$ попарно ортогональны, так как
$\langle E_iP_cx,E_jP_dy\rangle=\langle x,P_cE_i^*E_jP_dy\rangle=0$. Пусть $N_c$ — проектор на $S_c$ при
$c<C$ и $N_C=1-\sum_{c<C}N_c$ — тоже проектор. Для $\rho$, сосредоточенного в $P_c\mathbb C^D$, образ
$\Phi(\rho)=\sum_iE_iP_c\rho P_cE_i^*$ лежит в $S_c$, а $N_c$ тождествен на $S_c$ (при $c=C$ —
потому что $S_C\perp S_{c'}$, $c'<C$), поэтому $\operatorname{Tr}(N_c\Phi(\rho))=\operatorname{Tr}\Phi(\rho)=1$.
(1)$\Rightarrow$(2). Для единичного вектора $x\in P_c\mathbb C^D$ и $d\neq c$ числа
$\operatorname{Tr}(N_{d'}\Phi(|x\rangle\langle x|))\ge0$ в сумме дают $1$, и число с $d'=c$ равно $1$, поэтому
$0=\operatorname{Tr}(N_d\Phi(|x\rangle\langle x|))=\sum_i\|N_d^{1/2}E_ix\|^2$; значит, $N_dE_iP_c=0$ и, после суммирования по
$d$, $N_cE_iP_c=E_iP_c$. Отсюда $P_cE_i^*E_jP_d=P_cE_i^*N_cE_jP_d=0$ при $c\neq d$.
(1),(2)$\Rightarrow$(3). Пусть $\omega_c=P_c/\operatorname{Tr}P_c$ и $\mathcal R(\sigma)=\sum_c\operatorname{Tr}(N_c\sigma)\,\omega_c$
(«измерить, затем приготовить») — канал с операторами Крауса $\sqrt\mu\,|u\rangle\langle f|N_c^{1/2}$, где $(\mu,u)$
пробегает собственные пары $\omega_c$, а $f$ — векторы ортонормированного базиса. Так как
$\operatorname{Tr}(\mathcal R(\sigma)P_d)=\operatorname{Tr}(N_d\sigma)$, имеем $\mathcal R^*(P_d)=N_d$. Используя $N_dE_iP_c=\delta_{cd}E_iP_c$
(из предыдущего шага), а затем 2, получаем
$$
P(\mathcal R\circ\Phi)^*(P_d)P=P\,\Phi^*(N_d)\,P=\sum_i\sum_{c,c'}P_{c'}E_i^*N_dE_iP_c=\sum_i\sum_{c'}P_{c'}E_i^*E_iP_d=P_d\Big(\sum_iE_i^*E_i\Big)P_d=P_d .
$$
(3)$\Rightarrow$(1). Пусть $N_c=\mathcal R^*(P_c)$ при $c<C$ и $N_C=\mathcal R^*(1-\sum_{c<C}P_c)$; они
положительны и в сумме дают $\mathcal R^*(1)=1$. Для $\rho$, сосредоточенного в $P_c\mathbb C^D$, и $d<C$ имеем
$\operatorname{Tr}(N_d\Phi(\rho))=\operatorname{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(P_d)P)=\operatorname{Tr}(\rho P_d)=\delta_{cd}$, а значит, и
$\operatorname{Tr}(N_C\Phi(\rho))=\delta_{cC}$. $\square$

**Теорема 10.7 (Бени–Кемпф–Крибс 2007; приводится без доказательства).** Алгебра фон Неймана $\mathfrak A$
на $\mathcal C$ корректируема для $\Phi$ тогда и только тогда, когда $PE_i^*E_jP\in\mathfrak A'$ для всех $i,j$, где
$\mathfrak A'$ — коммутант $\mathfrak A$ в $B(\mathcal C)$.

При $\mathfrak A=B(\mathcal C)$ имеем $\mathfrak A'=\mathbb CP$, и это теорема 10.5; предложение 10.6 — это
случай коммутативной $\mathfrak A$. При
$\mathfrak A=M_k\otimes1_m$ на $\mathcal C=\mathbb C^k\otimes\mathbb C^m$ получаются операторные (подсистемные)
коды Крибса–Лафламма–Пулена (*Phys. Rev. Lett.* 94, 180501, 2005), которые содержат бесшумные
подсистемы раздела 10.1 с $\mathcal R=\mathrm{id}$. Доказательство см. в работах C. Bény, A. Kempf, D. W.
Kribs, *Phys. Rev. Lett.* 98, 100502 (2007) и *Phys. Rev. A* 76, 042303 (2007).

**Пример.** Код повторения при независимых инверсиях фазы (вероятность $q$ на кубит), за которыми
следует не более одной инверсии бита. Операторы Крауса кратны $X_kZ^s$, $s\in\{0,1\}^3$, и каждое
$E_i^*E_j$ инвертирует не более двух битов, поэтому условие 2 выполнено для $P_0=|000\rangle\langle000|$,
$P_1=|111\rangle\langle111|$: классический бит выживает и считывается голосованием по большинству,
$N_0=\sum_{|b|\le1}|b\rangle\langle b|$, $N_1=1-N_0$. Кубит не выживает: $PZ_1P$ не кратно $P$,
и после восстановления $\mathcal R$ из раздела 10.2 состояние $|{+_L}\rangle=(|000\rangle+|111\rangle)/\sqrt2$ претерпевает логическую инверсию фазы
$Z_1Z_2Z_3$ ровно тогда, когда произошло нечётное число инверсий фазы; его верность равна
$(1+(1-2q)^3)/2=0.8645$ при $q=0.05$.

```python
q = 0.05
Phi2 = Phi_one @ channels.kraus_channel([op8(E) for E in independent(Z, q)])   # dephasing, then a flip
E2 = [Phi2.kraus[0, i, 0] for i in range(Phi2.kraus_rank)]        # its 32 Kraus operators
P0, P1 = torch.outer(ket('000'), ket('000')), torch.outer(ket('111'), ket('111'))
assert kl_matrix(E2, P) is None                                   # the qubit is not protected
assert max((P0 @ e.mH @ f @ P1).abs().max().item() for e in E2 for f in E2) == 0      # condition 2
N0 = torch.diag(torch.tensor([float(bin(i).count('1') <= 1) for i in range(8)], dtype=dt))
for N, Pd in ((N0, P0), (torch.eye(8, dtype=dt) - N0, P1)):       # P Phi^*(N_d) P = P_d
    assert (P @ sum(e.mH @ N @ e for e in E2) @ P - Pd).abs().max() < 1e-12
votes = states.born_probabilities(Phi2(rho_L), [op8(N0), op8(torch.eye(8, dtype=dt) - N0)])[0]
assert torch.allclose(votes, torch.stack([abs(a) ** 2, abs(b) ** 2]))      # the classical bit survives
plus_L = state8((ket('000') + ket('111')) / math.sqrt(2))
F_plus = (R @ Phi2)(plus_L).fidelity(plus_L).item()
print(f"fidelity of |+_L> {F_plus:.4f}, exact {(1 + (1 - 2 * q)**3) / 2:.4f}")
assert abs(F_plus - (1 + (1 - 2 * q) ** 3) / 2) < 1e-12
```

> **Физический смысл.** Защитить коммутативную алгебру — значит защитить классическую информацию (бит,
> метку сектора суперотбора), а это гораздо дешевле, чем защищать кубит: шум, коммутирующий с
> проекторами секторов, безвреден. «Гибридные коды» соответствуют $\bigoplus_cM_{k_c}$.

## 10.4 Миксеры с ограничениями

Квантовый анзац чередующихся операторов (Hadfield et al., *Algorithms* 12, 34, 2019) приготавливает
$|\beta,\gamma\rangle=e^{-i\beta_pH_M}e^{-i\gamma_pH_C}\cdots e^{-i\beta_1H_M}e^{-i\gamma_1H_C}|b_0\rangle$, где $H_C$ — диагональный
гамильтониан стоимости, а $H_M$ — **миксер**. Для ограничения $|b|=w$ (бисекция графа, портфели из $w$
активов) берут $|b_0|=w$ и **XY-миксер** на кольце, $H_M=\tfrac14\sum_{(i,j)}(X_iX_j+Y_iY_j)$.

**Предложение 10.8.** $H_M$, $H_C$, $e^{-i\beta H_M}$ и $e^{-i\gamma H_C}$ лежат в
$\mathcal M=\{P_0,\dots,P_n\}'=\bigoplus_wB(H_w)\cong\bigoplus_wM_{\binom nw}$, поэтому $|\beta,\gamma\rangle\in H_w$:
анзац никогда не покидает допустимое подпространство.

*Доказательство.* На кубитах $i,j$ имеем $X\otimes X+Y\otimes Y=2(|01\rangle\langle10|+|10\rangle\langle01|)$ (используйте $Y|0\rangle=i|1\rangle$,
$Y|1\rangle=-i|0\rangle$): $X_iX_j+Y_iY_j$ меняет местами биты $i,j$, если они различны, и аннулирует $|b\rangle$
в противном случае, поэтому переводит каждое $H_w$ в себя. Диагональный $H_C$ коммутирует с диагональными $P_w$,
как и экспоненты (теорема 1.5(d)). $\square$

В библиотеке $\mathcal M$ — это `SpinChain(n).algebra` (бит $1$ = спин вверх, так что $S^z_i=-\tfrac12Z_i$), а $H_M$ —
её XXZ-гамильтониан $J\sum_{(i,j)}\big(\tfrac12(S^+_iS^-_j+S^-_iS^+_j)+\Delta S^z_iS^z_j\big)$ с $J=1$,
$\Delta=0$, поскольку $S^+_iS^-_j+S^-_iS^+_j=\tfrac12(X_iX_j+Y_iY_j)$; мы сверяем его с плотной
матрицей $64\times64$. В качестве $H_C$ берём бисекцию графа из двух треугольников $\{0,1,2\}$, $\{3,4,5\}$, соединённых
рёбрами $(2,3)$, $(5,0)$: $H_C=\sum_{(i,j)\in E}\tfrac12(1-Z_iZ_j)$, $Z_iZ_j=4S^z_iS^z_j$, считает разрезанные
рёбра. Лучшая бисекция разрезает $2$ ребра; мы начинаем с $|101010\rangle$, которое разрезает $6$. Наконец, весь
анзац пересчитывается с плотными матрицами $64\times64$ независимо от библиотеки.

```python
n, w = 6, 3
chain = SpinChain(n, 'periodic')
M = chain.algebra                                                 # (+)_w M_binom(6, w), padded to 20
H_M = chain.xxz(J=1.0, Delta=0.0)
H_dense = sum(on(X, i, n) @ on(X, j, n) + on(Y, i, n) @ on(Y, j, n) for i, j in chain.bonds) / 4
spec = torch.sort(torch.cat([e[0] for e in H_M.eigenvalues()]))[0]
assert torch.allclose(spec, torch.linalg.eigvalsh(H_dense))       # same spectrum, sector by sector
edges = chain.bonds + [(0, 2), (3, 5)]                            # two triangles joined by (2,3), (5,0)
H_C = M.zero()
for i, j in edges:
    H_C = H_C + (M.identity() - chain.sz(i) @ chain.sz(j) * 4.0) * 0.5
layers = [dynamics.propagator(H_M, beta) @ dynamics.propagator(H_C, gamma)
          for gamma, beta in ((0.6, 1.0), (1.4, 1.4))]
rho0 = ideal = chain.basis_state('101010')
for U_l in layers:
    ideal = channels.unitary_channel(U_l)(ideal)
assert abs(ideal.sector_probabilities()[0, w].item() - 1) < 1e-12           # Proposition 10.8

HC_dense = sum(torch.eye(64, dtype=dt) - on(Z, i, n) @ on(Z, j, n) for i, j in edges) / 2
psi = torch.eye(64, dtype=dt)[int('101010', 2)]                   # the same circuit on 6 full qubits
for gamma, beta in ((0.6, 1.0), (1.4, 1.4)):
    psi = torch.linalg.matrix_exp(-1j * beta * H_dense) @ torch.linalg.matrix_exp(-1j * gamma * HC_dense) @ psi
weight = torch.tensor([bin(b).count('1') for b in range(64)])
assert abs(psi[weight == w].norm().item() - 1) < 1e-12           # Proposition 10.8, dense
assert abs((psi.conj() @ HC_dense @ psi).real.item() - states.expectation(ideal, H_C).real.item()) < 1e-10
```

**Амплитудное затухание** кубита $i$ ($|1\rangle\to|0\rangle$ с вероятностью $g$) имеет операторы Крауса
$K_0^{(i)}=|0\rangle\langle0|_i+\sqrt{1-g}\,|1\rangle\langle1|_i$ и $K_1^{(i)}=\sqrt g\,|0\rangle\langle1|_i$; последний понижает
вес: это отображение между секторами, `chain.site_amplitude_damping(i, g)` (объект `InterSectorChannel`).

**Предложение 10.9 (утечка и постселекция).** Пусть $\rho_0$ сосредоточено в $H_w$, $U_l$ — унитарные
элементы $\mathcal M$, $\mathcal D$ — затухание всех кубитов и $\rho_l=\mathcal D(U_l\rho_{l-1}U_l^*)$.
Тогда каждое $\rho_l$ блочно-диагонально с весами $\le w$ и
$$
P_w\rho_pP_w=(1-g)^{wp}\;U\rho_0U^*,\qquad U=U_p\cdots U_1 :
$$
вес $w$ обнаруживается с вероятностью $(1-g)^{wp}$, а состояние после постселекции совпадает с бесшумным.

*Доказательство.* Операторы Крауса канала $\mathcal D$ — это $K_s=\prod_iK^{(i)}_{s_i}$, $s\in\{0,1\}^n$; $K_s$ переводит
$H_v$ в $H_{v-|s|}$, т. е. $P_vK_s=K_sP_{v+|s|}$ ($P_v=0$ при $v>n$). Если $\sigma$ блочно-диагональна с
весами $\le w$, то $P_vK_s\sigma K_s^*P_{v'}=K_sP_{v+|s|}\sigma P_{v'+|s|}K_s^*$ обращается в нуль при $v\neq v'$ и
при $v>w-|s|$: $\mathcal D(\sigma)$ блочно-диагонально с весами $\le w$, и его блок веса $w$ возникает
только от $s=0$. Так как $K_0|b\rangle=(1-g)^{|b|/2}|b\rangle$, имеем $P_w\mathcal D(\sigma)P_w=(1-g)^wP_w\sigma P_w$. Операторы $U_l$
коммутируют с $P_v$; остаётся провести индукцию по $l$. $\square$

Канал `InterSectorChannel` сохраняет только диагональные блоки своего выхода; здесь других блоков нет,
так что ничего не теряется. При $g=0.05$ и $p=2$ вероятность принятия равна $0.95^6=0.7351$:

```python
g = 0.05
damping = [chain.site_amplitude_damping(i, g) for i in range(n)]
noisy = rho0
for U_l in layers:
    noisy = channels.unitary_channel(U_l)(noisy)
    for D in damping:                                             # one qubit after the other
        noisy = D(noisy)
print(noisy.sector_probabilities()[0])                            # weight leaks to w = 2, 1, 0
P_w = M.central([1.0 if v == w else 0.0 for v in range(n + 1)])
post, accept = noisy.condition_on(P_w)                            # measure the weight, keep w = 3
cut = lambda r: states.expectation(r, H_C).real.item()
print(f"P(weight 3) {accept.item():.4f}, fidelity {noisy.fidelity(ideal).item():.4f}, "
      f"post-selected {post.fidelity(ideal).item():.4f}; "
      f"cut: ideal {cut(ideal):.4f}, noisy {cut(noisy):.4f}, post-selected {cut(post):.4f}")
assert abs(accept.item() - (1 - g) ** (w * len(layers))) < 1e-12              # Proposition 10.9
assert abs(post.fidelity(ideal).item() - 1) < 1e-10 and abs(cut(post) - cut(ideal)) < 1e-10
```

Недопустимая часть снижает ожидаемый разрез с $4.3789$ до $4.3123$ — мнимое улучшение, которое
ввело бы в заблуждение оптимизатор углов; после постселекции значение точное.

**Непрерывное время.** То же верно для уравнения Линдблада (теорема 8.6, [глава 8](08_dynamics.md))
$$
\dot\rho=-i[H,\rho]+\gamma_1\sum_i\Big(S^-_i\rho S^+_i-\tfrac12\{S^+_iS^-_i,\rho\}\Big)+\gamma_2\sum_k\Big(L_k\rho L_k^*-\tfrac12\{L_k^*L_k,\rho\}\Big)
$$
с $H\in\mathcal M$ и скачками $L_k\in\mathcal M$ внутри секторов, например чистой дефазировкой $L_k=S^z_k$.
Здесь $\sum_iS^+_iS^-_i=\sum_i|1\rangle\langle1|_i=N$ — **оператор числа**, $NP_v=vP_v$: он централен.

**Предложение 10.10.** Если $\rho(0)$ сосредоточено в $H_w$, то $P_w\rho(t)P_w=e^{-\gamma_1wt}\sigma(t)$,
где $\sigma$ — решение того же уравнения с $\gamma_1=0$ и $\sigma(0)=\rho(0)$. Следовательно, вероятность принятия равна
$e^{-\gamma_1wt}$, а состояние после постселекции $\sigma(t)$ не зависит от $\gamma_1$.

*Доказательство.* Пусть $B_{vv'}=P_v\rho P_{v'}$. Так как $H,L_k,N$ коммутируют с $P_v$ и $P_vS^-_i=S^-_iP_{v+1}$,
имеем $\dot B_{vv'}=-i[H,B_{vv'}]+\gamma_1\sum_iS^-_iB_{v+1,v'+1}S^+_i-\tfrac{\gamma_1}2(v+v')B_{vv'}+\gamma_2\mathcal L_2(B_{vv'})$,
где $\mathcal L_2$ — дефазирующая часть. Блоки с $\max(v,v')>w$ удовлетворяют замкнутой линейной системе с
нулевыми начальными значениями, поэтому они равны нулю (единственность в лемме 8.1). Тогда $\dot B_{ww}=-i[H,B_{ww}]-\gamma_1wB_{ww}+\gamma_2\mathcal L_2(B_{ww})$, так что
$e^{\gamma_1wt}B_{ww}$ решает уравнение для $\sigma$ с тем же начальным значением и потому совпадает с $\sigma$
(лемма 8.1). Это уравнение снова имеет форму GKSL, поэтому $\operatorname{Tr}\sigma(t)=1$ (теорема 8.6). $\square$

Сравните с предложением 8.11, которое допускает нецентральную сумму $\sum_jL_j^*L_j$, но не допускает скачков внутри секторов: здесь
член потерь $\gamma_1N$ централен, поэтому множитель отсутствия скачков — скаляр $e^{-\gamma_1wt}$, и дефазировка
внутри сектора допустима.

При $\gamma_2=0$ имеем $\sigma(t)=e^{-iHt}\rho(0)e^{iHt}$: постселекция идеальна. В сопутствующей статье
(`examples/paper/ex9_qaoa_xy_mixer.py`) состояние $|111000\rangle$ эволюционирует под действием $H_M$ до $T=2$; мы воспроизводим
её строку $\gamma_1=\gamma_2=0.1$: $P(w=3)=e^{-0.6}=0.5488$, верность $0.4805$, после постселекции $0.8755$ (как при $\gamma_1=0$).

```python
T = 2.0
jumps = [chain.lowering(i) for i in range(n)] + [chain.sz(i) for i in range(n)]
rho0 = chain.basis_state('111000')                                # the initial state of the paper
ideal_T = dynamics.von_neumann(rho0, H_M, [T])[0]
post_T = {}
for g1, g2 in ((0.1, 0.0), (0.0, 0.1), (0.1, 0.1)):          # the last one is printed
    rho_T = dynamics.lindblad_evolve(rho0, H_M, jumps, [0.0, T], rates=[g1] * n + [g2] * n, substeps=80)[-1]
    post_T[g1, g2], accept_T = rho_T.condition_on(P_w)
    assert abs(accept_T.item() - math.exp(-w * g1 * T)) < 1e-8                 # Proposition 10.10
print(f"P(weight 3) {accept_T.item():.4f}, fidelity {rho_T.fidelity(ideal_T).item():.4f}, "
      f"post-selected {post_T[0.1, 0.1].fidelity(ideal_T).item():.4f}")
assert abs(post_T[0.1, 0.0].fidelity(ideal_T).item() - 1) < 1e-6               # gamma_2 = 0: perfect
assert post_T[0.1, 0.1].trace_distance(post_T[0.0, 0.1]).item() < 1e-6       # independent of gamma_1
```

> **Физический смысл.** Энергетическая релаксация ($T_1$), одна из главных ошибок сверхпроводящих кубитов, —
> это скачок в нижний сектор. Измерение веса Хэмминга — центральный проектор, безвредный внутри
> сектора (раздел 3.3), — обнаруживает её, а постселекция устраняет её ценой доли принятых
> запусков. Дефазировка коммутирует с проекторами секторов и не обнаруживается.

## 10.5 С чем сравнивать

Qiskit (`qiskit.quantum_info`, модели шума Aer, `QAOAAnsatz` с пользовательским миксером), Cirq, PennyLane
(`qml.qaoa.xy_mixer`, `default.mixed`) и QuTiP/`qutip-qip` моделируют всё это на $n$ полных кубитах
и являются подходящими инструментами для схем; stim сэмплирует клиффордовы схемы с шумом паулиевского типа
на тысячах кубитов (амплитудное затухание не является паулиевским каналом). Нам не известна функция в этих
пакетах, которая вычисляла бы алгебру взаимодействия канала, её коммутант или всю алгебру неподвижных точек
(`steadystate` в QuTiP возвращает одно стационарное состояние линдбладиана); о численных блочных разложениях
матричных $*$-алгебр см. Murota et al., *Japan J. Indust. Appl. Math.* 27 (2010). Библиотека добавляет хранение в
$\bigoplus_wM_{\binom nw}$, отображения между секторами и постселекцию как центральный проектор, всё в пакетном режиме.
Одному сектору нужно $\binom nw^2$ элементов вместо $4^n$ ($1.7\cdot10^8$ против $4.3\cdot10^9$ при $n=16$,
$w=8$), но вся алгебра дополняется до $(n+1)\binom n{\lfloor n/2\rfloor}^2$ элементов ($2.8\cdot10^9$), что даёт
лишь скромную экономию. Библиотека не является симулятором схем и не содержит формализма стабилизаторов.

## Упражнения

1. (Теореме 10.2 нужна унитальность; ср. раздел 6.4.) Для амплитудного затухания одного кубита,
   $E_0=|0\rangle\langle0|+\sqrt{1-g}\,|1\rangle\langle1|$, $E_1=\sqrt g\,|0\rangle\langle1|$, $0<g<1$, покажите, что
   $\mathcal A'=\mathbb C1$, $\operatorname{Fix}(\Phi)=\mathbb C|0\rangle\langle0|$ и $\operatorname{Fix}(\Phi^*)=\mathbb C1$.
2. Докажите, что $\sum_w\binom nw^2=\binom{2n}n$. Каковы $\mathcal A$ и
   $\operatorname{Fix}(\Phi)$ для канала из следствия 10.3 при $\theta=\pi/2$ (что нарушает его условие)?
   Проверьте $\dim\operatorname{Fix}(\Phi)$ при $n=4$.
3. (Линейность коррекции ошибок.) Пусть $PE_a^*E_bP=\lambda_{ab}P$ для конечного числа операторов $E_a$ (не
   обязательно операторов Крауса канала), и пусть $\mathcal R$ построено по $\Lambda$, как в доказательстве
   теоремы 10.5. Покажите, что $\mathcal R$ исправляет всякий канал, операторы Крауса которого лежат в
   $\operatorname{span}\{E_a\}$. Выведите отсюда, что код против инверсий бита исправляет когерентную ошибку $e^{-i\phi X_1}$.

<details>
<summary>Решения</summary>

1. $E_1$, $E_1^*$ и их произведения дают все матричные единицы, поэтому $\mathcal A=M_2$ и $\mathcal A'=\mathbb C1$.
   Для $T=(x_{kl})$ имеем $\Phi(T)=\begin{pmatrix}x_{00}+gx_{11}&\sqrt{1-g}\,x_{01}\\\sqrt{1-g}\,x_{10}&(1-g)x_{11}\end{pmatrix}=T$
   тогда и только тогда, когда $x_{11}=x_{01}=x_{10}=0$, и
   $\Phi^*(T)=\begin{pmatrix}x_{00}&\sqrt{1-g}\,x_{01}\\\sqrt{1-g}\,x_{10}&(1-g)x_{11}+gx_{00}\end{pmatrix}$,
   что равно $T$ тогда и только тогда, когда $T=x_{00}1$. ($\Phi$ не унитален, а $\Phi^*$ унитален, но не сохраняет след.)

```python
E_ad = [torch.tensor([[1, 0], [0, math.sqrt(0.7)]], dtype=dt), torch.tensor([[0, math.sqrt(0.3)], [0, 0]], dtype=dt)]
fix_ad = fixed_points(E_ad)
assert len(fix_ad) == 1 and len(commutant_basis(E_ad + [E.mH for E in E_ad])) == 1
assert torch.allclose(fix_ad[0] / fix_ad[0][0, 0], torch.tensor([[1, 0], [0, 0]], dtype=dt))
assert len(fixed_points([E.mH for E in E_ad])) == 1                     # Fix(Phi^*) = C 1
```

2. Пересчитаем $n$-элементные подмножества $\{1,\dots,2n\}$ по числу $w$ элементов в первой половине:
   $\binom{2n}n=\sum_w\binom nw\binom n{n-w}=\sum_w\binom nw^2$. При $\theta=\pi/2$ оператор
   $U=(-i)^n\sum_w(-1)^wP_w$ различает лишь чётность $w$, поэтому $\mathcal A=\operatorname{span}\{P_{\mathrm{even}},P_{\mathrm{odd}}\}$,
   где $P_{\mathrm{even}}=\sum_{w\text{ чётно}}P_w$, и по теореме 10.2 $\operatorname{Fix}(\Phi)$ состоит из
   операторов, блочно-диагональных относительно чётности: размерность $2\cdot(2^{n-1})^2=128$ при $n=4$.

```python
U2 = torch.linalg.matrix_exp(-1j * (math.pi / 2) * sum(on(Z, j, 4) for j in range(4)))
assert len(fixed_points([U2 / math.sqrt(2), U2.mH / math.sqrt(2)])) == 128
assert sum(math.comb(4, v) ** 2 for v in range(5)) == math.comb(8, 4) == 70
```

3. При $F_k=\sum_au_{ak}E_a$ имеем $\operatorname{span}\{E_a\}=\operatorname{span}\{F_k\}$, поэтому каждый оператор Крауса канала имеет вид
   $G_m=\sum_kh_{mk}F_k$, и $G_mP=\sum_{d_k>0}h_{mk}\sqrt{d_k}V_k$. Так как $V_l^*V_k=\delta_{lk}P$ и $(1-Q)V_k=0$,
   получаем $\mathcal R(G_m\rho G_m^*)=\sum_l|h_{ml}|^2d_l\,\rho$ для $\rho=P\rho P$, так что $\mathcal R(\Phi(\rho))=c\rho$; взяв следы,
   находим $c=1$. Для $E_a\in\{1,X_1,X_2,X_3\}$ имеем $\Lambda=1$, $V_k=X_kP$, и $\mathcal R$ — восстановление из раздела
   10.2; $e^{-i\phi X_1}=\cos\phi\,1-i\sin\phi\,X_1$ лежит в линейной оболочке.

```python
for phi in (0.2, 1.0, 2.5):
    G1 = torch.linalg.matrix_exp(-1j * phi * on(X, 0, 3))
    out = (R @ channels.kraus_channel([op8(G1)]))(rho_L)
    assert (out.matrix - rho_L.matrix).abs().max() < 1e-12       # corrected exactly
```

</details>

**Далее:** в [главе 11](11_optimisation.md) рассматриваются оптимизация по унитарным группам алгебры и неравенства для следов.
