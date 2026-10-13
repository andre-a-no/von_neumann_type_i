[English](../D_gksl.md) · [Español](../es/D_gksl.md) · [Français](../fr/D_gksl.md) · [Deutsch](../de/D_gksl.md) · [中文](../zh/D_gksl.md) · [日本語](../ja/D_gksl.md) · **Русский**

# Приложение D. Генераторы квантовых динамических полугрупп

В главе 8 доказано, что генератор в форме Горини–Коссаковского–Сударшана–Линдблада (GKSL) порождает
квантовую динамическую полугруппу (теорема 8.6), и использовано обратное утверждение (теорема 8.7) и его вариант для
$\mathcal M=\bigoplus_cM_{k_c}$. В этом приложении оба утверждения доказываются в конечной размерности, вместе с
описанием свободы в выборе формы GKSL и критерием единственности стационарного состояния. Инструментом служит
матрица Чоя из главы 6: генератор, вообще говоря, не является вполне положительным, но его матрица Чоя положительна на
ортогональном дополнении максимально запутанного вектора. Рассуждение того же рода — положительность
семейства операторов, затем разложение — доказывает в последнем разделе теорему о коррекции ошибок для операторных
алгебр, сформулированную в главе 10 (теорема 10.7), с явным каналом восстановления.

**Вы узнаете**, почему генератор полугруппы вполне положительных отображений условно вполне
положителен, как прочитать $H$ и операторы скачков по матрице Чоя генератора, в каком смысле
они единственны, почему всякий генератор на алгебре с секторами имеет вид, использованный в разделе 8.3, когда
стационарное состояние единственно и какие алгебры наблюдаемых может защитить канал восстановления.

Всюду далее $M_n=M_n(\mathbb C)$, $E_{ij}=|i\rangle\langle j|$ — его матричные единицы,
**супероператор** — это линейное отображение на $M_n$ (или на $\mathcal M$), а
$$
\Omega=\sum_{i=1}^ne_i\otimes e_i\in\mathbb C^n\otimes\mathbb C^n,\qquad
J(\Psi)=\sum_{i,j=1}^nE_{ij}\otimes\Psi(E_{ij})=[\Psi(E_{ij})]_{i,j=1}^n
$$
— максимально запутанный вектор и **матрица Чоя** супероператора $\Psi$ (теорема 6.2;
здесь $\Psi$ не обязан быть положительным). Для векторов $v,w$ через $|v\rangle\langle w|$ обозначается оператор
$x\mapsto\langle w,x\rangle v$. Супероператор $\Psi$ **сохраняет эрмитовость**, если
$\Psi(A^*)=\Psi(A)^*$ для всех $A$. Квантовые динамические полугруппы и их генераторы — те же, что в
разделе 8.2 (лемма 8.5): $T_t=e^{t\mathcal L}$, где $T_t$ — каналы.

## D.1 Матрицы Чоя составных элементов

**Лемма D.1.** Пусть $\Psi$ — супероператор на $M_n$.

(a) $\Psi\mapsto J(\Psi)$ — линейная биекция множества супероператоров на $M_n$ на $M_n(M_n)=M_{n^2}$.
(b) Всякий $v\in\mathbb C^n\otimes\mathbb C^n$ представляется в виде $v=(1\otimes K)\Omega$ ровно с одним $K\in M_n$, и
$\langle(1\otimes K)\Omega,(1\otimes L)\Omega\rangle=\operatorname{Tr}(K^*L)$. В частности,
$\langle\Omega,(1\otimes K)\Omega\rangle=\operatorname{Tr}K$.
(c) Для $K,L\in M_n$ супероператор $A\mapsto KAL^*$ имеет матрицу Чоя
$|(1\otimes K)\Omega\rangle\langle(1\otimes L)\Omega|$. В частности, $J(\mathrm{id})=|\Omega\rangle\langle\Omega|$.
(d) $\Psi$ сохраняет эрмитовость тогда и только тогда, когда $J(\Psi)$ самосопряжена. Всякое CP-отображение сохраняет эрмитовость.

*Доказательство.* (a) Блок $(i,j)$ матрицы $J(\Psi)$ равен $\Psi(E_{ij})$, а $E_{ij}$ образуют базис $M_n$, поэтому
$J$ линейно и инъективно, и всякая блочная матрица $[X_{ij}]$ — матрица Чоя супероператора
$E_{ij}\mapsto X_{ij}$. (b) У $(1\otimes K)\Omega=\sum_ie_i\otimes Ke_i$ $i$-й блок равен $Ke_i$, т. е.
$i$-му столбцу $K$; эти столбцы определяют $K$, и всякий $v$ получается таким образом. Далее,
$\langle\sum_ie_i\otimes Ke_i,\sum_je_j\otimes Le_j\rangle=\sum_i\langle Ke_i,Le_i\rangle=\operatorname{Tr}(K^*L)$;
для последнего утверждения заменим $(K,L)$ на $(1,K)$. (c) Блок $(i,j)$ правой части равен
$Ke_i(Le_j)^*=KE_{ij}L^*$. (d) Блок $(i,j)$ матрицы $J(\Psi)^*$ равен $\Psi(E_{ji})^*$. Если $\Psi$
сохраняет эрмитовость, то $\Psi(E_{ji})^*=\Psi(E_{ji}^*)=\Psi(E_{ij})$. Обратно, если
$\Psi(E_{ji})^*=\Psi(E_{ij})$ для всех $i,j$ и $A=\sum a_{ij}E_{ij}$, то
$\Psi(A)^*=\sum\bar a_{ij}\Psi(E_{ij})^*=\sum\bar a_{ij}\Psi(E_{ji})=\Psi(A^*)$. CP-отображение имеет форму Крауса
(теорема 6.2), и $(\sum_lK_lAK_l^*)^*=\sum_lK_lA^*K_l^*$. $\square$

**Лемма D.2 (разложение).** Пусть $q_1,\dots,q_C$ — ортонормированные векторы в
$\mathbb C^n\otimes\mathbb C^n$, записанные как $q_c=(1\otimes Q_c)\Omega$, пусть $W=\operatorname{span}\{q_c\}$, и пусть
$\Psi$ — сохраняющий эрмитовость супероператор, для которого $\langle v,J(\Psi)v\rangle\ge0$ при всех $v\perp W$.
Положим $J=J(\Psi)$ и $\xi_c=Jq_c-\frac12\sum_d\langle q_d,Jq_c\rangle q_d$, и определим $K_c\in M_n$ условием
$(1\otimes K_c)\Omega=\xi_c$. Тогда существуют $L_1,\dots,L_r\in M_n$, где $r=\operatorname{rank}(P_W^\perp JP_W^\perp)$,
а $P_W^\perp$ — проектор на $W^\perp$, такие, что
$$
\Psi(A)=\sum_{j=1}^rL_jAL_j^*+\sum_{c=1}^C\big(K_cAQ_c^*+Q_cAK_c^*\big)\qquad(A\in M_n),
$$
$\operatorname{Tr}(Q_c^*L_j)=0$ для всех $c,j$ и $\operatorname{Tr}(L_j^*L_k)=0$ при $j\neq k$.

*Доказательство.* $J$ самосопряжена по лемме D.1(d). Пусть $P=\sum_c|q_c\rangle\langle q_c|$ — проектор на $W$,
и $P^\perp=P_W^\perp=1-P$. Раскрывая $J=(P+P^\perp)J(P+P^\perp)$ и используя $PJP^\perp=PJ-PJP$ и
$P^\perp JP=JP-PJP$, получаем
$$
J=P^\perp JP^\perp+\big(PJ+JP-PJP\big).
$$
Первое слагаемое самосопряжено, и $\langle v,P^\perp JP^\perp v\rangle=\langle P^\perp v,JP^\perp v\rangle\ge0$,
поскольку $P^\perp v\perp W$; значит, $P^\perp JP^\perp\ge0$. По теореме 1.3
$P^\perp JP^\perp=\sum_{j=1}^r|l_j\rangle\langle l_j|$, где $l_j$ — собственные векторы, отвечающие ненулевым
собственным значениям $\mu_j$, умноженные на $\sqrt{\mu_j}$: они попарно ортогональны и лежат в образе
$P^\perp JP^\perp$, который ортогонален $W$. Запишем $l_j=(1\otimes L_j)\Omega$; по лемме D.1(b)
$\operatorname{Tr}(Q_c^*L_j)=\langle q_c,l_j\rangle=0$ и $\operatorname{Tr}(L_j^*L_k)=\langle l_j,l_k\rangle=0$ при $j\ne k$,
а по лемме D.1(c) $P^\perp JP^\perp$ — матрица Чоя отображения $A\mapsto\sum_jL_jAL_j^*$. Для второго слагаемого:
$JP=\sum_c|Jq_c\rangle\langle q_c|$ и $PJP=\sum_{c,d}\langle q_d,Jq_c\rangle|q_d\rangle\langle q_c|$, поэтому
$$
\sum_c|\xi_c\rangle\langle q_c|=JP-\tfrac12PJP,\qquad\sum_c|q_c\rangle\langle\xi_c|=PJ-\tfrac12PJP
$$
(второе равенство сопряжено первому; $PJP$ самосопряжена), а их сумма равна $PJ+JP-PJP$. По
лемме D.1(c) $|\xi_c\rangle\langle q_c|$ и $|q_c\rangle\langle\xi_c|$ — матрицы Чоя отображений
$A\mapsto K_cAQ_c^*$ и $A\mapsto Q_cAK_c^*$. Итак, обе части доказываемого тождества имеют матрицу Чоя $J$,
и лемма D.1(a) даёт тождество. $\square$

## D.2 Условная вполне положительность

Супероператор $\Psi$ на $M_n$ называется **условно вполне положительным** (CCP), если
$\langle v,J(\Psi)v\rangle\ge0$ для всякого $v\perp\Omega$. По лемме D.1(b) эти $v$ — в точности векторы
$(1\otimes K)\Omega$ с $\operatorname{Tr}K=0$.

**Теорема D.3.** Для супероператора $\mathcal L$ на $M_n$ следующие условия равносильны:

1. $e^{t\mathcal L}$ является CP при всяком $t\ge0$;
2. $\mathcal L$ сохраняет эрмитовость и условно вполне положителен;
3. существуют CP-отображение $\Phi$ на $M_n$ и $G\in M_n$, для которых $\mathcal L(A)=\Phi(A)+GA+AG^*$ при всех $A$.

В п. 3 можно взять $\Phi(A)=\sum_{j=1}^rL_jAL_j^*$ с $r\le n^2-1$, $\operatorname{Tr}L_j=0$ и
$\operatorname{Tr}(L_j^*L_k)=0$ при $j\neq k$, а $G$ — с $\operatorname{Tr}G\in\mathbb R$.

*Доказательство.* (1)$\Rightarrow$(2). По лемме 8.1 $(e^{t\mathcal L}-\mathrm{id})/t\to\mathcal L$ при $t\to0$, а $J$
линейно на конечномерном пространстве и, следовательно, непрерывно, поэтому
$J(e^{t\mathcal L})=J(\mathrm{id})+t\,J\big((e^{t\mathcal L}-\mathrm{id})/t\big)$, где последняя матрица Чоя стремится к
$J(\mathcal L)$. Для $v\perp\Omega$ имеем $\langle v,J(\mathrm{id})v\rangle=|\langle\Omega,v\rangle|^2=0$ (лемма D.1(c)), поэтому
при $t>0$, поскольку $J(e^{t\mathcal L})\ge0$ (теорема 6.2),
$$
0\le\frac1t\langle v,J(e^{t\mathcal L})v\rangle=\big\langle v,J\big((e^{t\mathcal L}-\mathrm{id})/t\big)v\big\rangle
\;\xrightarrow[t\to0]{}\;\langle v,J(\mathcal L)v\rangle .
$$
Каждое $e^{t\mathcal L}$ является CP, а значит, сохраняет эрмитовость (лемма D.1(d)), поэтому
$\mathcal L(A^*)=\lim_t(e^{t\mathcal L}(A)^*-A^*)/t=\mathcal L(A)^*$.

(2)$\Rightarrow$(3). Применим лемму D.2 с $C=1$ и $q_1=\Omega/\sqrt n$, т. е. $Q_1=1/\sqrt n$:
$\mathcal L(A)=\sum_jL_jAL_j^*+GA+AG^*$ с $G=K_1/\sqrt n$, и $\Phi=\sum_jL_j\cdot L_j^*$ является CP
(теорема 6.2). Далее, $\operatorname{Tr}L_j=\sqrt n\operatorname{Tr}(Q_1^*L_j)=0$, $r\le\operatorname{rank}P_W^\perp=n^2-1$,
и по лемме D.1(b)
$\operatorname{Tr}G=\langle\Omega,\xi_1\rangle/\sqrt n=\langle q_1,\xi_1\rangle=\frac12\langle q_1,J(\mathcal L)q_1\rangle\in\mathbb R$,
так как $J(\mathcal L)$ самосопряжена.

(3)$\Rightarrow$(1). Шаги 1–4 доказательства теоремы 8.6 используют лишь то, что генератор равен
$\mathcal L_0+\Phi$, где $\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)$ для некоторой матрицы
$H_{\mathrm{eff}}$, а $\Phi$ является CP. При $H_{\mathrm{eff}}=iG$ имеем $\mathcal L_0(\rho)=G\rho+\rho G^*$, и
эти шаги показывают, что $e^{t\mathcal L}=\sum_{m\ge0}S_m(t)$ является CP. $\square$

Одной положительности недостаточно. Пусть $n\ge2$, $\theta$ — транспонирование на $M_n$ и $\mathcal L=\theta-\mathrm{id}$. Так как $\theta^2=\mathrm{id}$,
$$
e^{t\mathcal L}=e^{-t}\big(\cosh t\,\mathrm{id}+\sinh t\,\theta\big)=\tfrac{1+e^{-2t}}2\,\mathrm{id}+\tfrac{1-e^{-2t}}2\,\theta ,
$$
— выпуклая комбинация положительных отображений, сохраняющих след (предложение 6.1). Но $J(\theta)=[E_{ji}]=\sum_{i,j}E_{ij}\otimes E_{ji}$ —
это перестановка $e_k\otimes e_l\mapsto e_l\otimes e_k$, поэтому она переводит ортогональный $\Omega$ вектор
$v=e_1\otimes e_2-e_2\otimes e_1$ в $-v$ (как в доказательстве предложения 6.1); так как
$\|v\|^2=2$ и $\langle v,J(\mathrm{id})v\rangle=|\langle\Omega,v\rangle|^2=0$, получаем
$\langle v,J(\mathcal L)v\rangle=-2$: $\mathcal L$ не является CCP, и действительно
$\langle v,J(e^{t\mathcal L})v\rangle=-(1-e^{-2t})<0$ при всяком $t>0$. Полугруппа положительных отображений, сохраняющих след,
не обязана состоять из каналов.

Проверим теорему D.3 на генераторе GKSL с произвольными скачками на $M_3$, построенном функцией `lindblad_superoperator`
(матрица $S$ отображения $\mathcal L$ на развёрнутом по строкам $\operatorname{vec}\rho$, так что $J_{(i,a),(j,b)}=S_{(a,b),(i,j)}$, как в
разделе 6.1), и на $\theta-\mathrm{id}$:

```python
import math
import numpy as np
import scipy.linalg as sla
import torch
from torch_vn_algebra import TypeIAlgebra, channels, dynamics

torch.manual_seed(0)
rng = np.random.default_rng(0)
f64, c128 = torch.float64, torch.complex128

def choi_of(S, n):
    """J_{(i,a),(j,b)} = S_{(a,b),(i,j)}: the Choi matrix of the map with row-major matrix S (Section 6.1)."""
    return S.reshape(n, n, n, n).transpose(2, 0, 3, 1).reshape(n * n, n * n)

def superop(f, n):
    """Row-major matrix of a linear map f on n x n matrices."""
    E = np.eye(n * n).reshape(n * n, n, n)
    return np.stack([f(X).reshape(-1) for X in E], axis=1)

vec_of = lambda K: K.T.reshape(-1)                        # (1 (x) K) Omega = sum_i e_i (x) K e_i
op_of = lambda v, n: v.reshape(n, n).T                    # its inverse (Lemma D.1(b))
cplx = lambda *s: rng.standard_normal(s) + 1j * rng.standard_normal(s)

n = 3
M3 = TypeIAlgebra([n], [n], precision='double')
H = cplx(n, n)
H = (H + H.conj().T) / 2
Ls = [cplx(n, n) / 2 for _ in range(2)]                   # arbitrary jumps, not traceless
to_op = lambda X: M3.operator(torch.from_numpy(X).reshape(1, 1, n, n))
S = dynamics.lindblad_superoperator(M3, to_op(H), [to_op(L) for L in Ls])[0][0].numpy()
J = choi_of(S, n)
q = vec_of(np.eye(n)) / math.sqrt(n)                      # the unit vector Omega / sqrt(n)
Pperp = np.eye(n * n) - np.outer(q, q)
assert np.abs(J - J.conj().T).max() < 1e-12               # Hermiticity preserving (Lemma D.1(d))
print(f"min eigenvalue of J(L): {np.linalg.eigvalsh(J).min():.2f}; "
      f"on the complement of Omega: {np.linalg.eigvalsh(Pperp @ J @ Pperp).min():.1e}")
assert np.linalg.eigvalsh(J).min() < -0.1                 # L itself is not CP ...
assert np.linalg.eigvalsh(Pperp @ J @ Pperp).min() > -1e-12        # ... but CCP (Theorem D.3)
for t in (0.01, 0.1, 1.0):
    assert np.linalg.eigvalsh(choi_of(sla.expm(t * S), n)).min() > -1e-12

L_theta = superop(lambda X: X.T, n) - np.eye(n * n)     # L = transpose - id: Hermiticity preserving, not CCP
anti = vec_of(np.eye(n)[:, [1]] @ np.eye(n)[[0], :] - np.eye(n)[:, [0]] @ np.eye(n)[[1], :])   # e1(x)e2 - e2(x)e1
assert abs(anti @ q) < 1e-15 and abs(anti @ choi_of(L_theta, n) @ anti + 2) < 1e-12
X = cplx(n, n)
for t in (0.01, 0.1, 1.0):
    T_t = sla.expm(t * L_theta)
    b = (1 - math.exp(-2 * t)) / 2
    assert abs(anti @ choi_of(T_t, n) @ anti + (1 - math.exp(-2 * t))) < 1e-12        # <v, J(e^{tL}) v>
    assert np.abs(T_t - ((1 - b) * np.eye(n * n) + b * superop(lambda Y: Y.T, n))).max() < 1e-12
    assert abs(np.trace((T_t @ X.reshape(-1)).reshape(n, n)) - np.trace(X)) < 1e-12   # trace preserving
    assert np.linalg.eigvalsh(choi_of(T_t, n)).min() < -0.9 * b                     # but not CP
```

## D.3 Теорема GKSL

**Теорема D.4 (Горини–Коссаковский–Сударшан; Линдблад; теорема 8.7).** Супероператор $\mathcal L$ на
$M_n$ является генератором квантовой динамической полугруппы на $M_n$ тогда и только тогда, когда существуют $H=H^*$ и
$L_1,\dots,L_r$ в $M_n$, для которых
$$
\mathcal L(\rho)=-i[H,\rho]+\sum_{j=1}^r\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
$$
Можно выбрать $r\le n^2-1$, $\operatorname{Tr}H=0$, $\operatorname{Tr}L_j=0$ и $\operatorname{Tr}(L_j^*L_k)=0$ при $j\neq k$.

*Доказательство.* Часть «если» — это теорема 8.6. Обратно, пусть $(T_t)$ — квантовая динамическая полугруппа с генератором
$\mathcal L$, так что $T_t=e^{t\mathcal L}$ (лемма 8.5). Поскольку каждое $T_t$ является CP, теорема D.3 даёт
$\mathcal L(\rho)=\sum_jL_j\rho L_j^*+G\rho+\rho G^*$ с $r\le n^2-1$ бесследовыми попарно ортогональными $L_j$ и
$\operatorname{Tr}G\in\mathbb R$. Положим $H=\frac i2(G-G^*)$ и $R=\frac12(G+G^*)$; обе матрицы самосопряжены, и
$G=-iH+R$, $G^*=iH+R$, $G\rho+\rho G^*=-i[H,\rho]+\{R,\rho\}$. Кроме того,
$\operatorname{Tr}H=\frac i2(\operatorname{Tr}G-\overline{\operatorname{Tr}G})=0$. Поскольку каждое $T_t$ является TP,
$\operatorname{Tr}\mathcal L(\rho)=\frac d{dt}\operatorname{Tr}T_t(\rho)\big|_{t=0}=0$ (лемма 8.1), и по цикличности
$$
0=\operatorname{Tr}\mathcal L(\rho)=\operatorname{Tr}\Big(\rho\,\Big(\sum_jL_j^*L_j+2R\Big)\Big)\qquad\text{for all }\rho .
$$
При $\rho=X=\sum_jL_j^*L_j+2R$, которая самосопряжена, получаем $\operatorname{Tr}(X^*X)=0$, поэтому $X=0$ и
$R=-\frac12\sum_jL_j^*L_j$, что и есть форма GKSL. $\square$

**Теорема D.5 (свобода в форме GKSL).** Пусть $(H;L_1,\dots,L_r)$ и $(H';L'_1,\dots,L'_s)$ — данные
GKSL ($H,H'$ самосопряжены) с генераторами $\mathcal L$, $\mathcal L'$, как в теореме D.4.

(a) $\mathcal L'=\mathcal L$ в каждом из следующих случаев:
(i) $L'_j=L_j$ и $H'=H+c1$, $c\in\mathbb R$;
(ii) $H'=H$ и, после дополнения более короткого семейства нулями, $L'_j=\sum_ku_{jk}L_k$ с унитарной $u$;
(iii) $L'_j=L_j+a_j1$ с $a_j\in\mathbb C$ и $H'=H-\frac i2\sum_j(\bar a_jL_j-a_jL_j^*)$; эта матрица самосопряжена.

(b) Если $\operatorname{Tr}H=\operatorname{Tr}H'=0$, все $L_j,L'_j$ бесследовые и $\mathcal L=\mathcal L'$, то $H=H'$ и,
после дополнения нулями, $L'_j=\sum_ku_{jk}L_k$ с унитарной $u$.

(c) Следовательно, два набора данных GKSL имеют один и тот же генератор тогда и только тогда, когда один получается из
другого конечным числом шагов (i)–(iii). В частности, оператор $H$ нормированной формы из теоремы D.4 и её
диссипативная часть $\Phi(\rho)=\sum_jL_j\rho L_j^*$ определяются генератором $\mathcal L$.

*Доказательство.* (a) (i) $[c1,\rho]=0$. (ii) По теореме 6.7 $\Phi'=\sum_jL'_j\cdot L_j'^*$ совпадает с
$\Phi=\sum_jL_j\cdot L_j^*$; сопряжённое отображение не зависит от выбора операторов Крауса (раздел 6.2), поэтому
$\sum_jL_j'^*L'_j=\Phi^*(1)=\sum_jL_j^*L_j$. (iii) Для одного скачка $L$ и $a\in\mathbb C$ раскрытие скобок даёт
$$
(L+a)\rho(L+a)^*-\tfrac12\{(L+a)^*(L+a),\rho\}-\Big(L\rho L^*-\tfrac12\{L^*L,\rho\}\Big)=\tfrac12[B,\rho],
\qquad B=\bar aL-aL^* ,
$$
поскольку слагаемые $\bar aL\rho+a\rho L^*$ из первого произведения складываются с
$-\frac12(aL^*\rho+a\rho L^*+\bar aL\rho+\bar a\rho L)$, а слагаемые с $|a|^2$ сокращаются. Так как $B^*=-B$,
матрица $\frac i2B$ самосопряжена, поэтому $H'=H-\frac i2\sum_jB_j$ самосопряжена, и
$-i[H'-H,\rho]=-\frac12\sum_j[B_j,\rho]$ компенсирует изменение диссипатора.

(b) Пусть $\Phi=\sum_jL_j\cdot L_j^*$ и $G=-iH-\frac12\Phi^*(1)$, так что $\mathcal L(\rho)=\Phi(\rho)+G\rho+\rho G^*$.
По лемме D.1(c)
$J(\mathcal L)=\sum_j|l_j\rangle\langle l_j|+|g\rangle\langle\Omega|+|\Omega\rangle\langle g|$, где $l_j=(1\otimes L_j)\Omega$ и
$g=(1\otimes G)\Omega$, причём $l_j\perp\Omega$, поскольку $\operatorname{Tr}L_j=0$ (лемма D.1(b)). Обозначая через $P^\perp$ проектор
на $\Omega^\perp$, получаем $P^\perp J(\mathcal L)P^\perp=\sum_j|l_j\rangle\langle l_j|=J(\Phi)$ и аналогично
$P^\perp J(\mathcal L')P^\perp=J(\Phi')$. Значит, $J(\Phi)=J(\Phi')$, $\Phi=\Phi'$ (лемма D.1(a)), и теорема 6.7 даёт
унитарную $u$. Тогда $\Phi^*(1)=\Phi'^*(1)$, откуда $[H-H',\rho]=0$ для всех $\rho$. Матрица $X$, коммутирующая со всеми
$E_{ij}$, кратна $1$ (сравните элементы $XE_{ij}$ и $E_{ij}X$), поэтому $H-H'=\lambda1$, и
$\lambda=0$, как видно из взятия следа.

(c) По (a) шаги (i)–(iii) не меняют генератор, и шаг, обратный к каждому шагу, — шаг того же
вида ((iii) с $-a_j$, так как $\bar a_jL'_j-a_jL_j'^*=\bar a_jL_j-a_jL_j^*$). Любые данные приводятся к
нормированной форме шагом (iii) с $a_j=-\operatorname{Tr}L_j/n$ и затем шагом (i) с $c=-\operatorname{Tr}H''/n$, где $H''$ —
гамильтониан после (iii); два нормированных набора данных с одним и тем же генератором связаны шагом (ii) в силу (b). $\square$

Доказательство теоремы D.4 — это алгоритм: $\Phi$ получается из сжатия $J(\mathcal L)$ на $\Omega^\perp$, затем
$G$ — из $\xi_1$, затем $H$. Применим его к генератору выше, восстановим генератор с помощью библиотеки и
сравним с нормированной формой исходных данных, вычисленной по теореме D.5(iii) и (i):

```python
def gksl_from_generator(S, n, tol=1e-10):
    """(H, [L_j]) with Tr H = Tr L_j = 0 from the matrix of a GKSL generator: the proof of Theorem D.4."""
    J = choi_of(S, n)
    q = vec_of(np.eye(n)) / math.sqrt(n)
    Pperp = np.eye(n * n) - np.outer(q, q)
    w, U = np.linalg.eigh(Pperp @ J @ Pperp)              # = J(Phi) >= 0
    jumps = [math.sqrt(x) * op_of(U[:, k], n) for k, x in enumerate(w) if x > tol]
    xi = J @ q - 0.5 * (q.conj() @ J @ q) * q
    G = op_of(xi, n) / math.sqrt(n)
    return 0.5j * (G - G.conj().T), jumps

H_rec, L_rec = gksl_from_generator(S, n)
S_rec = dynamics.lindblad_superoperator(M3, to_op(H_rec), [to_op(L) for L in L_rec])[0][0].numpy()
assert np.abs(S_rec - S).max() < 1e-12                    # the same generator
assert abs(np.trace(H_rec)) < 1e-12 and all(abs(np.trace(L)) < 1e-12 for L in L_rec)

a = [-np.trace(L) / n for L in Ls]                        # Theorem D.5(iii): make the jumps traceless
L_norm = [L + aj * np.eye(n) for L, aj in zip(Ls, a)]
H_norm = H - 0.5j * sum(np.conj(aj) * L - aj * L.conj().T for L, aj in zip(Ls, a))
H_norm = H_norm - np.trace(H_norm) / n * np.eye(n)        # Theorem D.5(i)
assert np.abs(H_norm - H_norm.conj().T).max() < 1e-12
assert np.abs(H_rec - H_norm).max() < 1e-12               # H is unique (Theorem D.5(b))
jump_choi = lambda Ks: sum(np.outer(vec_of(K), vec_of(K).conj()) for K in Ks)
assert np.abs(jump_choi(L_rec) - jump_choi(L_norm)).max() < 1e-12      # so is sum_j L_j . L_j^*
Wr, Wn = np.stack([vec_of(K) for K in L_rec], 1), np.stack([vec_of(K) for K in L_norm], 1)
u = np.linalg.lstsq(Wr, Wn, rcond=None)[0]                # L_norm_j = sum_k u_kj L_rec_k
assert np.abs(Wr @ u - Wn).max() < 1e-12 and np.abs(u.conj().T @ u - np.eye(2)).max() < 1e-12   # unitary
print("eigenvalues of the normalised H:", np.round(np.linalg.eigvalsh(H_rec), 4), " jumps:", len(L_rec))
```

**Теорема D.6 (единственность стационарного состояния).** Пусть $\mathcal L$ — генератор GKSL из теоремы D.4
с данными $(H;L_1,\dots,L_r)$. Если у $\mathcal L$ есть точное стационарное состояние $\sigma$ ($\sigma>0$) и единственные
матрицы, коммутирующие с $H$ и всеми $L_j$, $L_j^*$, — это кратные $1$, то $\ker\mathcal L=\mathbb C\sigma$:
$\sigma$ — единственное стационарное состояние. Кроме того, $\ker\mathcal L^*=\mathbb C1$.

*Доказательство.* Сопряжённое относительно следового спаривания, $\operatorname{Tr}(\mathcal L(\rho)A)=\operatorname{Tr}(\rho\,\mathcal L^*(A))$, равно
$\mathcal L^*(A)=i[H,A]+\sum_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})$: по цикличности
$\operatorname{Tr}(-i[H,\rho]A)=\operatorname{Tr}(\rho\,i[H,A])$, $\operatorname{Tr}(L\rho L^*A)=\operatorname{Tr}(\rho L^*AL)$ и
$\operatorname{Tr}(\{L^*L,\rho\}A)=\operatorname{Tr}(\rho\{L^*L,A\})$ (это упражнение 1 главы 8).
Прямое вычисление, в котором гамильтоновы слагаемые сокращаются, поскольку $A\mapsto i[H,A]$ — дифференцирование, даёт
для всякого $X$
$$
\mathcal L^*(X^*X)-X^*\mathcal L^*(X)-\mathcal L^*(X^*)X=\sum_j[L_j,X]^*[L_j,X]:
$$
для одного скачка из шести слагаемых с антикоммутаторами четыре попарно сокращаются, а оставшиеся два,
$\frac12X^*L^*LX+\frac12X^*L^*LX$, складываются, так что левая часть равна $L^*X^*XL-X^*L^*XL-L^*X^*LX+X^*L^*LX$, т. е.
$(XL-LX)^*(XL-LX)$. Пусть $X\in\ker\mathcal L^*$. Тогда $\mathcal L^*(X^*)=\mathcal L^*(X)^*=0$, так как $\mathcal L^*$
сохраняет эрмитовость, поэтому $\mathcal L^*(X^*X)=\sum_j[L_j,X]^*[L_j,X]=:Y\ge0$ и
$\operatorname{Tr}(\sigma Y)=\operatorname{Tr}(\mathcal L(\sigma)X^*X)=0$. Поскольку $\sigma>0$, величина
$\operatorname{Tr}(\sigma Y)=\operatorname{Tr}(\sigma^{1/2}Y\sigma^{1/2})$ при $\sigma^{1/2}Y\sigma^{1/2}\ge0$ обращается в нуль, лишь если
$\sigma^{1/2}Y\sigma^{1/2}=0$, т. е. $Y=0$; взяв след, получаем $\sum_j\operatorname{Tr}([L_j,X]^*[L_j,X])=0$ с неотрицательными
слагаемыми, поэтому $[L_j,X]=0$ для всех $j$. То же рассуждение для $X^*\in\ker\mathcal L^*$
даёт $[L_j,X^*]=0$, т. е. $[L_j^*,X]=0$. Тогда $X$ коммутирует с $L_j^*L_j$, диссипативная часть $\mathcal L^*(X)$
обращается в нуль, и $0=\mathcal L^*(X)=i[H,X]$. По предположению $X\in\mathbb C1$; обратно, $\mathcal L^*(1)=0$. Итак,
$\ker\mathcal L^*=\mathbb C1$. В базисах $(E_{ij})$ и $(E_{ji})$ пространства $M_n$, двойственных относительно спаривания
$(\rho,A)\mapsto\operatorname{Tr}(\rho A)$, матрица $\mathcal L^*$ — транспонированная матрица $\mathcal L$; поэтому у них
одинаковый ранг, и $\dim\ker\mathcal L=\dim\ker\mathcal L^*=1$. Состояние из $\mathbb C\sigma$ со следом $1$ — это $\sigma$. $\square$

Это конечномерный случай критерия Фриджерио [*Commun. Math. Phys.* 63, 269 (1978)],
используемого в разделе 8.5. Оба предположения необходимы. В упражнении 4 главы 8 состояние $\tfrac12\,1$ точное
и стационарное, но $H$ и $\sigma_z$ коммутируют со всеми диагональными матрицами, и стационарные состояния порождают
двумерное пространство диагональных матриц. Обратно, пусть на $M_3$ $H=0$, $L_1=|e_1\rangle\langle e_3|$,
$L_2=|e_2\rangle\langle e_3|$. $*$-алгебра, порождённая $L_1,L_2$, содержит $L_1L_2^*=E_{12}$ и $L_1^*L_1=E_{33}$, а значит, и
все матричные единицы, поэтому с $L_j,L_j^*$ коммутируют лишь кратные $1$; но всякое $\rho=P\rho P$ с $P=E_{11}+E_{22}$
стационарно, так как $L_j\rho=0=\rho L_j^*$ и $\{E_{33},\rho\}=0$, поэтому $\dim\ker\mathcal L\ge4$; по теореме D.6 ни одно стационарное
состояние не является точным. Для случайного генератора выше код проверяет предположения и вывод
(матрица $\mathcal L^*$ на развёрнутых по строкам векторах $A^T$ равна $S^T$: $\operatorname{vec}(\mathcal L^*(A)^T)=S^T\operatorname{vec}(A^T)$), а также второй пример:

```python
def kernel_dim(S, tol=1e-9):
    return int((np.linalg.svd(S, compute_uv=False) < tol).sum())

def commutant_dim(ops, n):                    # dim {X : XA = AX for all A in ops}
    return kernel_dim(np.concatenate([np.kron(A, np.eye(n)) - np.kron(np.eye(n), A.T) for A in ops]))

w, V = np.linalg.eig(S)
sigma = V[:, np.argmin(abs(w))].reshape(n, n)             # row-major: vec(X)_{(a,b)} = X_ab
sigma = sigma / np.trace(sigma)
assert np.abs(S @ sigma.reshape(-1)).max() < 1e-12        # stationary
assert np.abs(sigma - sigma.conj().T).max() < 1e-10 and np.linalg.eigvalsh(sigma).min() > 1e-3   # faithful
assert commutant_dim([H] + Ls + [L.conj().T for L in Ls], n) == 1
assert kernel_dim(S) == 1 and kernel_dim(S.T) == 1        # Theorem D.6
assert np.abs(S.T @ np.eye(n).reshape(-1)).max() < 1e-12  # ker L^* = C 1

e3 = np.eye(3)
L_abs = [np.outer(e3[0], e3[2]), np.outer(e3[1], e3[2])]       # L_1 = |e_1><e_3|, L_2 = |e_2><e_3|
S_abs = superop(lambda R: sum(L @ R @ L.T - 0.5 * (L.T @ L @ R + R @ L.T @ L) for L in L_abs), 3)
assert commutant_dim(L_abs + [L.T for L in L_abs], 3) == 1
print("second example: dim ker L =", kernel_dim(S_abs))
assert kernel_dim(S_abs) == 4                             # no faithful stationary state, no uniqueness
```

## D.4 Алгебры с секторами

Пусть $\mathcal M=\bigoplus_{c=1}^CM_{k_c}\subseteq M_D$ — блочно-диагональные матрицы, $D=\sum_ck_c$, с
минимальными центральными проекторами $P_c$, и пусть $E(X)=\sum_cP_cXP_c$ — условное ожидание $M_D$ на
$\mathcal M$. Как в доказательстве предложения 8.8, $E(B)=B$, $E(XB)=E(X)B$ и $E(BX)=BE(X)$ для $B\in\mathcal M$,
и $E(X^*)=E(X)^*$. Для CP-отображения $\Phi:\mathcal M\to\mathcal M$ через $\Phi^*$ обозначается сопряжённое к нему относительно `Tr_blunt` (предложение 6.8).

**Теорема D.7.** Пусть $\mathcal L$ — генератор квантовой динамической полугруппы $(T_t)$ на $\mathcal M$. Существуют
CP-отображение $\Phi:\mathcal M\to\mathcal M$ и $H=H^*\in\mathcal M$, для которых
$$
\mathcal L(\rho)=\Phi(\rho)-\tfrac12\{\Phi^*(1),\rho\}-i[H,\rho]\qquad(\rho\in\mathcal M).
$$
По предложению 6.4 $\Phi(\rho)_d=\sum_{c,l}K^{dc}_l\rho_cK^{dc*}_l$ с операторами Крауса между секторами;
обратно, всякий такой $\mathcal L$ порождает квантовую динамическую полугруппу на $\mathcal M$ (предложение 8.8(1)).
Итак, генераторы из раздела 8.3 (где $\Phi=\sum_jg_jL_j\cdot L_j^*+\sum_mg_m\Phi_m$) — это в точности
генераторы квантовых динамических полугрупп на $\mathcal M$.

*Доказательство.* *Шаг 1: CP-отображение на $M_D$.* По предложению 6.4
$T_t(\rho)_d=\sum_{c,l}K^{dc}_l(t)\rho_cK^{dc}_l(t)^*$. Пусть $N^{dc}_l(t)\in M_D$ — матрица с блоком $K^{dc}_l(t)$ в
позиции $(d,c)$ и нулями в остальных позициях. Для $X\in M_D$ матрица $N^{dc}_lXN^{dc*}_l$ — это блок
$K^{dc}_lX_{cc}K^{dc*}_l$ в позиции $(d,d)$, где $X_{cc}$ — $c$-й диагональный блок, т. е. $c$-й блок
$E(X)$. Следовательно, $\hat T_t=T_t\circ E=\sum_{d,c,l}N^{dc}_l\cdot N^{dc*}_l$ является CP на $M_D$ (теорема 6.2), причём
$\hat T_0=E$. По лемме 8.5 $(\hat T_t-E)/t=\big((T_t-\mathrm{id})/t\big)\circ E\to\mathcal L\circ E$.

*Шаг 2: положительность вне секторных векторов.* Пусть $p_c=(1\otimes P_c)\Omega\in\mathbb C^D\otimes\mathbb C^D$ и
$W=\operatorname{span}\{p_c\}$. По лемме D.1(c) $J(E)=\sum_c|p_c\rangle\langle p_c|$, поэтому для $v\perp W$ и $t>0$
$0\le\frac1t\langle v,J(\hat T_t)v\rangle=\langle v,J((\hat T_t-E)/t)v\rangle\to\langle v,J(\mathcal L\circ E)v\rangle$.
Каждое $T_t$ сохраняет эрмитовость (в силу своей формы Крауса), значит, и $\mathcal L$ тоже (как в теореме D.3), как и $E$
и $\mathcal L\circ E$.

*Шаг 3: разложение.* По лемме D.1(b) $\langle p_c,p_d\rangle=\operatorname{Tr}(P_cP_d)=\delta_{cd}k_c$, поэтому
векторы $q_c=p_c/\sqrt{k_c}=(1\otimes Q_c)\Omega$ с $Q_c=P_c/\sqrt{k_c}$ ортонормированы. Лемма D.2 для
$\Psi=\mathcal L\circ E$ даёт CP-отображение $\hat\Phi=\sum_j\hat L_j\cdot\hat L_j^*$ на $M_D$ и $K_c\in M_D$, для которых
$$
\mathcal L(E(A))=\hat\Phi(A)+\sum_c\big(K_cAP_c+P_cAK_c^*\big)/\sqrt{k_c}\qquad(A\in M_D).
$$
Для $A\in\mathcal M$ имеем $E(A)=A$ и $AP_c=P_cA$, поэтому сумма равна $\tilde GA+A\tilde G^*$ с
$\tilde G=\sum_cK_cP_c/\sqrt{k_c}$.

*Шаг 4: сжатие.* Для $A\in\mathcal M$ имеем $\mathcal L(A)\in\mathcal M$, поэтому
$\mathcal L(A)=E(\mathcal L(A))=E(\hat\Phi(A))+GA+AG^*$ с $G=E(\tilde G)\in\mathcal M$. Отображение
$\Phi=E\circ\hat\Phi$ на $\mathcal M$ имеет вид
$\Phi(A)_d=\sum_{c,j}(P_d\hat L_jP_c)A_c(P_d\hat L_jP_c)^*$: оно является CP по предложению 6.4, а его операторы Крауса —
блоки матриц $\hat L_j$.

*Шаг 5: сохранение следа.* Как в теореме D.4, $\operatorname{Tr}\mathcal L(A)=0$ для $A\in\mathcal M$, и
$\operatorname{Tr}\Phi(A)=\operatorname{Tr}(A\Phi^*(1))$, поэтому $\operatorname{Tr}(AX)=0$ для всех $A\in\mathcal M$, где
$X=\Phi^*(1)+G+G^*\in\mathcal M$ самосопряжена. При $A=X$ получаем $X=0$ (след `Tr_blunt` точен, теорема 3.1).
При $H=\frac i2(G-G^*)\in\mathcal M$ имеем $G=-iH-\frac12\Phi^*(1)$ и $G\rho+\rho G^*=-i[H,\rho]-\frac12\{\Phi^*(1),\rho\}$. $\square$

Векторы $p_c$ заменяют $\Omega=\sum_cp_c$, потому что $\hat T_0=E$, а не $\mathrm{id}$: матрица $J(\mathcal L\circ E)$ не обязана быть
положительной на $\Omega^\perp$, и в примере ниже она таковой не является.

**Следствие D.8 (картина Гейзенберга).** Пусть $(S_t)_{t\ge0}$ — непрерывное семейство унитальных CP-отображений на
$\mathcal M$ с $S_0=\mathrm{id}$ и $S_{s+t}=S_sS_t$. Его генератор имеет вид $A\mapsto\Psi(A)+K^*A+AK$, где
$\Psi$ является CP на $\mathcal M$, $K\in\mathcal M$ и $\Psi(1)+K+K^*=0$.

*Доказательство.* Пусть $T_t=S_t^*$ (сопряжённые относительно `Tr_blunt`). Они являются CP (предложение 6.8) и TP (следствие 6.9(b)),
$T_0=\mathrm{id}$ и $T_{s+t}=(S_tS_s)^*=S_s^*S_t^*=T_sT_t$, поскольку
$\operatorname{Tr}(S_tS_s(x)y)=\operatorname{Tr}(xS_s^*S_t^*(y))$ и $S_tS_s=S_{s+t}$. Так как $\Psi\mapsto\Psi^*$ линейно на
конечномерном пространстве, $t\mapsto T_t$ непрерывно и $(T_t-\mathrm{id})/t\to\mathcal L_S^*$, где $\mathcal L_S$ —
генератор $(S_t)$. Теорема D.7 даёт $\mathcal L_S^*(\rho)=\Phi(\rho)+G\rho+\rho G^*$ с $G+G^*=-\Phi^*(1)$, а по
цикличности $\operatorname{Tr}((G\rho+\rho G^*)A)=\operatorname{Tr}(\rho(G^*A+AG))$, поэтому по следствию 6.9(c)
$\mathcal L_S(A)=\Phi^*(A)+G^*A+AG$. Возьмём $\Psi=\Phi^*$ (CP по предложению 6.8) и $K=G$. $\square$

Это конечномерный случай теоремы Кристенсена и Эванса, процитированной после теоремы 8.7. Проверим
теорему D.7 на $M_1\oplus M_2\subseteq M_3$ для генератора `lindblad_rhs` со скачком внутри
секторов и случайным скачком между ними, следуя пяти шагам доказательства:

```python
alg = TypeIAlgebra([1, 2], [1, 2], precision='double')    # M = M_1 (+) M_2 in M_3
D, ks = 3, alg.k_factors
P = [np.diag([1.0, 0, 0]), np.diag([0, 1.0, 1.0])]
E = lambda X: sum(Pc @ X @ Pc for Pc in P)               # the conditional expectation onto M
def blocks(X):                                            # dense block-diagonal 3 x 3 -> library layout
    out = torch.zeros(1, 2, 2, 2, dtype=c128)
    out[0, 0, :1, :1] = torch.from_numpy(X[:1, :1])
    out[0, 1] = torch.from_numpy(X[1:, 1:])
    return out
dense = lambda Y: sla.block_diag(Y[0, 0, :1, :1].numpy(), Y[0, 1].numpy())
Hs = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
inside = alg.operator(torch.randn(1, 2, 2, 2, dtype=c128))       # a jump inside the sectors
between = channels.random_inter_sector_channel(alg, alg, kraus_rank=2)
rhs = dynamics.lindblad_rhs(Hs, [inside, between], rates=[0.7, 1.3], algebra=alg)
L_M = lambda X: dense(rhs(0.0, blocks(X)))               # the generator on M (X block diagonal)
J = choi_of(superop(lambda X: L_M(E(X)), D), D)          # J(L o E) on M_3
p = [vec_of(Pc) for Pc in P]
Q = np.stack([pc / np.linalg.norm(pc) for pc in p], 1)    # orthonormal q_c = p_c / sqrt(k_c)
Pw = np.eye(D * D) - Q @ Q.T                              # projection onto the complement of W
assert np.linalg.eigvalsh(Pw @ J @ Pw).min() > -1e-12     # Step 2
w, U = np.linalg.eigh(Pw @ J @ Pw)                        # Step 3 (Lemma D.2)
jumps = [math.sqrt(x) * op_of(U[:, k], D) for k, x in enumerate(w) if x > 1e-10]
xis = [J @ Q[:, c] - 0.5 * sum((Q[:, d] @ J @ Q[:, c]) * Q[:, d] for d in range(2)) for c in range(2)]
G = E(sum(op_of(xi, D) @ Pc / math.sqrt(k) for xi, Pc, k in zip(xis, P, ks)))   # Step 4
Phi = lambda X: E(sum(K @ X @ K.conj().T for K in jumps))
Phi_dual_one = E(sum(K.conj().T @ K for K in jumps))      # Phi^*(1)
H_M = 0.5j * (G - G.conj().T)
assert np.abs(G + G.conj().T + Phi_dual_one).max() < 1e-12      # Step 5
for _ in range(3):
    X = E(cplx(D, D))
    gksl = Phi(X) - 0.5 * (Phi_dual_one @ X + X @ Phi_dual_one) - 1j * (H_M @ X - X @ H_M)
    assert np.abs(L_M(X) - gksl).max() < 1e-12            # Theorem D.7
    dual = lambda A: E(sum(K.conj().T @ A @ K for K in jumps)) + G.conj().T @ A + A @ G   # Corollary D.8
    assert abs(np.trace(L_M(X) @ (Y := E(cplx(D, D)))) - np.trace(X @ dual(Y))) < 1e-12
Pomega = np.eye(D * D) - np.outer(sum(p), sum(p)) / D     # projection onto the complement of Omega
print(f"numerical observation: min eigenvalue of J(L o E) on the complement of Omega: "
      f"{np.linalg.eigvalsh(Pomega @ J @ Pomega).min():.2f}")
assert np.linalg.eigvalsh(Pomega @ J @ Pomega).min() < -0.1
```

## D.5 Коррекция ошибок для операторных алгебр

Используем постановку раздела 10.3. $\Phi(\rho)=\sum_{i=1}^rE_i\rho E_i^*$ — канал на $M_D$, $\mathcal C=P\mathbb C^D$ —
код, $B(\mathcal C)=\{X\in M_D:X=PXP\}$ — операторы на $\mathcal C$, $\mathfrak A\subseteq B(\mathcal C)$ — $*$-алгебра
с единицей $P$ (по теореме 2.2 — алгебра фон Неймана на $\mathcal C$), а $\mathfrak A'$ — её коммутант в
$B(\mathcal C)$. $\mathfrak A$ **корректируема** для $\Phi$, если существует канал $\mathcal R$ на $M_D$, для которого
$P(\mathcal R\circ\Phi)^*(A)P=A$ при всех $A\in\mathfrak A$. Если у $\mathcal R$ операторы Крауса $R_m$, то
у $\mathcal R\circ\Phi$ операторы Крауса $R_mE_i$, и $\sum_mR_m^*R_m=1$.

**Теорема D.9 (Бени–Кемпф–Крибс; теорема 10.7).** $\mathfrak A$ корректируема для $\Phi$ тогда и только тогда, когда
$PE_i^*E_jP\in\mathfrak A'$ для всех $i,j$. В этом случае канал $\mathcal R$, построенный в доказательстве, корректирует $\mathfrak A$.

*Доказательство.* *Необходимость.* Пусть $\mathcal R$ корректирует $\mathfrak A$; положим $W_{mi}=R_mE_iP$, так что
$\sum_{m,i}W_{mi}^*XW_{mi}=P(\mathcal R\circ\Phi)^*(X)P$ для всех $X\in M_D$. Это равно $X$ при $X\in\mathfrak A$ и $P$ при
$X=1$, поскольку $(\mathcal R\circ\Phi)^*$ унитально (следствие 6.9(b)). Пусть $p\in\mathfrak A$ — проектор; тогда $p\le P$
и $\sum W_{mi}^*(1-p)W_{mi}=P-p$. Для $x=px$
$\sum_{m,i}\|(1-p)W_{mi}x\|^2=\langle x,(P-p)x\rangle=0$, а для $x=(P-p)x$
$\sum_{m,i}\|pW_{mi}x\|^2=\langle x,px\rangle=0$. Следовательно, $(1-p)W_{mi}p=0$ и $pW_{mi}(P-p)=0$, и, поскольку
$W_{mi}=W_{mi}P$,
$$
W_{mi}p=pW_{mi}p=pW_{mi}P=pW_{mi}.
$$
Всякий самосопряжённый $A\in\mathfrak A$ — линейная комбинация проекторов из $\mathfrak A$: по теореме 1.3 на $\mathcal C$
$A=\sum_k\lambda_kp_k$ с различными $\lambda_k$ и ортогональными проекторами $p_k$, $\sum_kp_k=P$, и
$p_k=\prod_{l\neq k}(A-\lambda_lP)/(\lambda_k-\lambda_l)\in\mathfrak A$ (примените обе части к собственным векторам в $\mathcal C$ и к $\mathcal C^\perp$; пустое произведение равно $P$). Всякий элемент
$\mathfrak A$ имеет вид $A_1+iA_2$ с самосопряжёнными $A_1,A_2\in\mathfrak A$. Поэтому $W_{mi}A=AW_{mi}$ для всех $A\in\mathfrak A$, и,
переходя к сопряжённым, $W_{mi}^*A=AW_{mi}^*$. Так как $\sum_mW_{mi}^*W_{mj}=PE_i^*(\sum_mR_m^*R_m)E_jP=PE_i^*E_jP$,
$$
A\,PE_i^*E_jP=\sum_mW_{mi}^*AW_{mj}=PE_i^*E_jP\,A\qquad(A\in\mathfrak A),
$$
т. е. $PE_i^*E_jP\in\mathfrak A'$.

*Достаточность.* По теореме 2.6 для $\mathfrak A$ на $\mathcal C$ существует ортонормированный базис $|c;a,b\rangle$
($c=1,\dots,C$, $a=1,\dots,k_c$, $b=1,\dots,m_c$) пространства $\mathcal C$ — образы $e_a\otimes f_b$ из $c$-го слагаемого при
унитарном $V$, — в котором $\mathfrak A$ состоит из операторов $\sum_{c,a,a',b}(A_c)_{aa'}|c;a,b\rangle\langle c;a',b|$,
$A_c\in M_{k_c}$, а $\mathfrak A'$ — из операторов $\sum_{c,a,b,b'}(Y_c)_{bb'}|c;a,b\rangle\langle c;a,b'|$, $Y_c\in M_{m_c}$.
Поэтому предположение даёт матрицы $\Lambda^c_{ij}\in M_{m_c}$, для которых
$$
\langle c;a,b|E_i^*E_j|d;a',b'\rangle=\delta_{cd}\,\delta_{aa'}\,(\Lambda^c_{ij})_{bb'} .
$$
*Шаг 1.* Для каждого $c$ и каждой пары $\alpha=(i,b)$ пусть $F^c_\alpha:\mathbb C^{k_c}\to\mathbb C^D$,
$F^c_\alpha e_a=E_i|c;a,b\rangle$, и $\Lambda^c_{\alpha\beta}=(\Lambda^c_{ij})_{bb'}$ для $\beta=(j,b')$. Выписанная формула означает, что
$F^{c*}_\alpha F^d_\beta=\delta_{cd}\Lambda^c_{\alpha\beta}1_{k_c}$: это условие Нилла–Лафламма для «ошибок» $F^c_\alpha$ на
коде $\mathbb C^{k_c}$, дополненное ортогональностью между секторами.

*Шаг 2.* Как в доказательстве теоремы 10.5, $\Lambda^c$ самосопряжена, так как
$\Lambda^c_{\beta\alpha}1=F^{c*}_\beta F^c_\alpha=(F^{c*}_\alpha F^c_\beta)^*=\bar\Lambda^c_{\alpha\beta}1$, и положительна, так как
$\sum\bar z_\alpha z_\beta\Lambda^c_{\alpha\beta}1=(\sum z_\alpha F^c_\alpha)^*(\sum z_\beta F^c_\beta)\ge0$. Запишем
$\Lambda^c=U^c\operatorname{diag}(d^c_\kappa)U^{c*}$ с унитарной $U^c$ и $d^c_\kappa\ge0$, и положим
$G^c_\kappa=\sum_\alpha u^c_{\alpha\kappa}F^c_\alpha$. Тогда $G^{c*}_\kappa G^d_\lambda=\delta_{cd}\delta_{\kappa\lambda}d^c_\kappa1$, поэтому
$G^c_\kappa=0$ при $d^c_\kappa=0$, а иначе $V^c_\kappa=G^c_\kappa/\sqrt{d^c_\kappa}$ — изометрия с
$V^{c*}_\kappa V^d_\lambda=\delta_{cd}\delta_{\kappa\lambda}1$. Обращая, получаем $F^c_\alpha=\sum_\kappa\bar u^c_{\alpha\kappa}G^c_\kappa
=\sum_{\kappa:d^c_\kappa>0}\bar u^c_{\alpha\kappa}\sqrt{d^c_\kappa}\,V^c_\kappa$, откуда
$$
V^{c'*}_\kappa F^d_\beta=\delta_{c'd}\,\bar u^d_{\beta\kappa}\sqrt{d^d_\kappa}\,1,\qquad
F^{c*}_\alpha V^{c'}_\kappa=\delta_{cc'}\,u^c_{\alpha\kappa}\sqrt{d^c_\kappa}\,1 .
$$
*Шаг 3: восстановление.* Матрицы $V^c_\kappa V^{c*}_\kappa$ ($d^c_\kappa>0$) — попарно ортогональные проекторы; пусть $Q$ — их
сумма. Пусть $\iota_c:\mathbb C^{k_c}\to\mathbb C^D$, $\iota_ce_a=|c;a,1\rangle$, — изометрия, и
$$
\mathcal R(\sigma)=\sum_{c}\sum_{\kappa:d^c_\kappa>0}\iota_cV^{c*}_\kappa\,\sigma\,V^c_\kappa\iota_c^*+(1-Q)\sigma(1-Q).
$$
Это канал: $\sum_{c,\kappa}V^c_\kappa\iota_c^*\iota_cV^{c*}_\kappa+(1-Q)^2=Q+1-Q=1$. Он измеряет сектор $c$ и
синдром $\kappa$, устраняет ошибку на $\mathbb C^{k_c}$ и записывает результат в копию $b=1$.

*Шаг 4: корректность.* $\mathcal R^*(X)=\sum_{c,\kappa}V^c_\kappa\iota_c^*X\iota_cV^{c*}_\kappa+(1-Q)X(1-Q)$, и для
$A\in\mathfrak A$ с блоками $A_c$ имеем $\iota_c^*A\iota_c=A_c$. Поскольку каждый $F^d_\beta$ — линейная комбинация $V^d_\kappa$,
$(1-Q)F^d_\beta=0$. С учётом шага 2 для $\alpha=(i,b)$, $\beta=(j,b')$
$$
\langle c;a,b|E_i^*\mathcal R^*(A)E_j|d;a',b'\rangle=\sum_{c',\kappa}\langle e_a,F^{c*}_\alpha V^{c'}_\kappa A_{c'}V^{c'*}_\kappa F^d_\beta e_{a'}\rangle
=\delta_{cd}\sum_\kappa u^c_{\alpha\kappa}d^c_\kappa\bar u^c_{\beta\kappa}(A_c)_{aa'}=\delta_{cd}(\Lambda^c_{ij})_{bb'}(A_c)_{aa'} .
$$
Положим $j=i$ и просуммируем по $i$: по выписанной выше формуле, так как $\sum_iE_i^*E_i=1$, имеем $\sum_i\Lambda^c_{ii}=1_{m_c}$, поэтому
$$
\langle c;a,b|(\mathcal R\circ\Phi)^*(A)|d;a',b'\rangle=\delta_{cd}\,\delta_{bb'}(A_c)_{aa'}=\langle c;a,b|A|d;a',b'\rangle .
$$
Оба оператора $P(\mathcal R\circ\Phi)^*(A)P$ и $A$ лежат в $B(\mathcal C)$ и имеют одинаковую матрицу в базисе $|c;a,b\rangle$, поэтому они
равны. $\square$

**Следствие D.10.** (a) Для $\mathfrak A=B(\mathcal C)$ имеем $\mathfrak A'=\mathbb CP$, и теорема D.9 — это теорема
Нилла–Лафламма 10.5 ($C=1$, $m_1=1$; шаг 3 — восстановление из теоремы 10.5). (b) Для $\mathfrak A=\operatorname{span}\{P_1,\dots,P_C\}$ с
попарно ортогональными проекторами $P_c\ne0$, $\sum_cP_c=P$, это равносильность условий 2 и 3 в предложении 10.6. (c) (*Бесшумные
подсистемы.*) В обозначениях предложения 10.1(a) пусть $\mathcal C=V(\mathbb C^{k_c}\otimes\mathbb C^{m_c})$ — $c$-е
слагаемое для алгебры взаимодействия $\mathcal A$ канала $\Phi$, и $\mathfrak A=V(1_{k_c}\otimes M_{m_c})V^*$ на $\mathcal C$ —
алгебра вида $M_{m}\otimes1_{k}$ после перестановки двух сомножителей. Тогда $\mathfrak A$ корректируется
каналом $\mathcal R=\mathrm{id}$, и $PE_i^*E_jP=V(E_{i,c}^*E_{j,c}\otimes1_{m_c})V^*\in\mathfrak A'$ ($P=VV^*$).

*Доказательство.* (a) $\mathfrak A'=\mathbb CP$ — это случай $C=1$, $m_1=1$ приведённого выше описания $\mathfrak A'$, а
$PE_i^*E_jP\in\mathbb CP$ — условие Нилла–Лафламма. Корректируемость $B(\mathcal C)$ — это коррекция в смысле
раздела 10.2: если $\mathcal R(\Phi(\rho))=\rho$ для всех состояний $\rho$, сосредоточенных в $\mathcal C$, то
$\operatorname{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(A)P)=\operatorname{Tr}(\rho A)$ для всех них, а эти $\rho$ порождают $B(\mathcal C)$;
обратное доказано в разделе 10.3. При $C=1$, $m_1=1$ операторы Крауса $\iota_1V^{1*}_\kappa$ из шага 3 совпадают с операторами
из теоремы 10.5, поскольку $\iota_1\iota_1^*=P$. (b) $T\in B(\mathcal C)$ коммутирует со всеми $P_c$ тогда и только тогда, когда
$P_cTP_d=0$ при $c\ne d$, а это условие 2 предложения 10.6. (c) Всякий $A\in\mathfrak A$, продолженный нулём на
остальные слагаемые, лежит в $\mathcal A'=V(\bigoplus_c1_{k_c}\otimes M_{m_c})V^*$ (теорема 2.6), поэтому $\Phi^*(A)=A$ по предложению
10.1(a), и $PAP=A$. Формула для $PE_i^*E_jP$ следует из $V^*E_iV=\bigoplus_cE_{i,c}\otimes1_{m_c}$, а
$E_{i,c}^*E_{j,c}\otimes1_{m_c}$ коммутирует с $1_{k_c}\otimes M_{m_c}$. $\square$

По теореме 2.6 всякая $\mathfrak A$ $*$-изоморфна некоторой $\bigoplus_cM_{k_c}$. Когда $C\ge2$ и некоторое $k_c\ge2$,
корректируемая $\mathfrak A$ защищает метку сектора $c$ (классическую информацию) вместе с квантовой информацией
каждой $M_{k_c}$; такие коды называются **гибридными**. В качестве примера возьмём $\mathfrak A=(M_2\otimes1_2)\oplus M_1$ на случайном
пятимерном коде в $\mathbb C^8$ и шум, который действует на калибровочные сомножители $\mathbb C^{m_c}$ посредством
$\sum_iB_{c,i}^*B_{c,i}=1_{m_c}$, переводит секторы в ортогональные подпространства, перемешивает $\mathcal C^\perp$ и затем смешивается
случайной унитарной матрицей (теорема 6.7). Условие Нилла–Лафламма нарушено, условие теоремы D.9 выполнено, и
восстановление из доказательства восстанавливает каждый $A\in\mathfrak A$, но не генерический оператор на $\mathcal C$:

```python
D, dims = 8, [(2, 2), (1, 1)]                             # A = (M_2 (x) 1_2) (+) (M_1 (x) 1_1), dim C = 5
U8 = np.linalg.qr(cplx(D, D))[0]
Vc = np.linalg.qr(cplx(D, 5))[0]                          # orthonormal basis |c; a, b> of the code
W = [Vc[:, :4], Vc[:, 4:]]                                # sector c: columns indexed by (a, b), row-major
Pcode = sum(Wc @ Wc.conj().T for Wc in W)
Y = [U8[:, :4], U8[:, 4:6]]                               # orthogonal images C^{k_c} (x) C^{s_c}, s_c = 2
r = 3
B = [np.linalg.qr(cplx(2 * r, m))[0].reshape(r, 2, m) for _, m in dims]   # sum_i B_i^* B_i = 1_{m_c}
E_code = [sum(Yc @ np.kron(np.eye(k), Bc[i]) @ Wc.conj().T for Yc, Bc, Wc, (k, m) in zip(Y, B, W, dims))
          for i in range(r)]
kraus = E_code + [np.linalg.qr(cplx(D, D))[0] @ (np.eye(D) - Pcode)]
mix = np.linalg.qr(cplx(r + 1, r + 1))[0]                 # Theorem 6.7: an equivalent family
kraus = [sum(mix[i, l] * K for l, K in enumerate(kraus)) for i in range(r + 1)]
assert np.abs(sum(K.conj().T @ K for K in kraus) - np.eye(D)).max() < 1e-12      # a channel

units = []                                                 # matrix units of A
for Wc, (k, m) in zip(W, dims):
    for a1 in range(k):
        for a2 in range(k):
            units.append(Wc @ np.kron(np.eye(k)[:, [a1]] @ np.eye(k)[[a2], :], np.eye(m)) @ Wc.conj().T)
comm = lambda X: max(np.abs(X @ A - A @ X).max() for A in units)
assert max(comm(Pcode @ Ki.conj().T @ Kj @ Pcode) for Ki in kraus for Kj in kraus) < 1e-12   # Theorem D.9
KL_dev = max(np.abs(Pcode @ Ki.conj().T @ Kj @ Pcode - np.trace(Pcode @ Ki.conj().T @ Kj) / 5 * Pcode).max()
             for Ki in kraus for Kj in kraus)
print(f"Knill-Laflamme condition violated by {KL_dev:.2f}")
assert KL_dev > 0.05

R_kraus, Qsum = [], np.zeros((D, D), complex)             # Steps 1-3 of the proof
for Wc, (k, m) in zip(W, dims):
    F = [Ki @ Wc[:, b::m] for Ki in kraus for b in range(m)]   # F_(i,b) e_a = E_i |c; a, b>
    Lam = np.array([[(Fa.conj().T @ Fb)[0, 0] for Fb in F] for Fa in F])
    assert max(np.abs(Fa.conj().T @ Fb - Lam[x, y] * np.eye(k)).max()
               for x, Fa in enumerate(F) for y, Fb in enumerate(F)) < 1e-12      # Step 1
    d, u = np.linalg.eigh(Lam)
    iota = Wc[:, 0::m]                                    # e_a -> |c; a, 1>
    for kk in range(len(d)):
        if d[kk] > 1e-12:
            Vk = sum(u[x, kk] * Fx for x, Fx in enumerate(F)) / math.sqrt(d[kk])
            R_kraus.append(iota @ Vk.conj().T)
            Qsum += Vk @ Vk.conj().T
R_kraus.append(np.eye(D) - Qsum)
assert np.abs(sum(R.conj().T @ R for R in R_kraus) - np.eye(D)).max() < 1e-10    # R is a channel
dual = lambda X: sum(Ki.conj().T @ R.conj().T @ X @ R @ Ki for R in R_kraus for Ki in kraus)
for _ in range(3):
    A = sum(complex(z) * Ue for z, Ue in zip(cplx(len(units)), units))
    assert np.abs(Pcode @ dual(A) @ Pcode - A).max() < 1e-10                       # A is corrected
B_gen = Pcode @ cplx(D, D) @ Pcode
print(f"a generic operator on the code: error {np.abs(Pcode @ dual(B_gen) @ Pcode - B_gen).max():.2f}")
assert np.abs(Pcode @ dual(B_gen) @ Pcode - B_gen).max() > 0.1
assert max(comm(R @ Ki @ Pcode) for R in R_kraus for Ki in kraus) < 1e-10   # proof of necessity: W_mi in A'
```

Результаты этого приложения принадлежат V. Gorini, A. Kossakowski и E. C. G. Sudarshan
[*J. Math. Phys.* 17, 821 (1976)] и G. Lindblad [*Commun. Math. Phys.* 48, 119 (1976)] (теоремы D.3–D.5),
E. Christensen и D. E. Evans [*J. London Math. Soc.* 20, 358 (1979)] (следствие D.8 в произвольных алгебрах
фон Неймана), A. Frigerio (теорема D.6) и C. Bény, A. Kempf и D. W. Kribs [*Phys. Rev. Lett.* 98, 100502
(2007); *Phys. Rev. A* 76, 042303 (2007)] (теорема D.9); приведённые выше доказательства полны для конечномерных
утверждений, используемых в курсе.
