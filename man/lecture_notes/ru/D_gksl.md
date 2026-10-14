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

*Вычисление:* численные проверки этого приложения — [руководство, приложение D](../../tutorial/ru/D_gksl.md).

<a id="s-D-1"></a>

## D.1 Матрицы Чоя составных элементов

<a id="r-D-1"></a>

**Лемма D.1.** Пусть $\Psi$ — супероператор на $M_n$.

(a) $\Psi\mapsto J(\Psi)$ — линейная биекция множества супероператоров на $M_n$ на $M_n(M_n)=M_{n^2}$.
(b) Всякий $v\in\mathbb C^n\otimes\mathbb C^n$ представляется в виде $v=(1\otimes K)\Omega$ ровно с одним $K\in M_n$, и
$\langle(1\otimes K)\Omega,(1\otimes L)\Omega\rangle=\mathrm{Tr}(K^*L)$. В частности,
$\langle\Omega,(1\otimes K)\Omega\rangle=\mathrm{Tr}\,K$.
(c) Для $K,L\in M_n$ супероператор $A\mapsto KAL^*$ имеет матрицу Чоя
$|(1\otimes K)\Omega\rangle\langle(1\otimes L)\Omega|$. В частности, $J(\mathrm{id})=|\Omega\rangle\langle\Omega|$.
(d) $\Psi$ сохраняет эрмитовость тогда и только тогда, когда $J(\Psi)$ самосопряжена. Всякое CP-отображение сохраняет эрмитовость.

*Доказательство.* (a) Блок $(i,j)$ матрицы $J(\Psi)$ равен $\Psi(E_{ij})$, а $E_{ij}$ образуют базис $M_n$, поэтому
$J$ линейно и инъективно, и всякая блочная матрица $[X_{ij}]$ — матрица Чоя супероператора
$E_{ij}\mapsto X_{ij}$. (b) У $(1\otimes K)\Omega=\sum_ie_i\otimes Ke_i$ $i$-й блок равен $Ke_i$, т. е.
$i$-му столбцу $K$; эти столбцы определяют $K$, и всякий $v$ получается таким образом. Далее,
$\langle\sum_ie_i\otimes Ke_i,\sum_je_j\otimes Le_j\rangle=\sum_i\langle Ke_i,Le_i\rangle=\mathrm{Tr}(K^*L)$;
для последнего утверждения заменим $(K,L)$ на $(1,K)$. (c) Блок $(i,j)$ правой части равен
$Ke_i(Le_j)^*=KE_{ij}L^*$. (d) Блок $(i,j)$ матрицы $J(\Psi)^*$ равен $\Psi(E_{ji})^*$. Если $\Psi$
сохраняет эрмитовость, то $\Psi(E_{ji})^*=\Psi(E_{ji}^*)=\Psi(E_{ij})$. Обратно, если
$\Psi(E_{ji})^*=\Psi(E_{ij})$ для всех $i,j$ и $A=\sum a_{ij}E_{ij}$, то
$\Psi(A)^*=\sum\bar a_{ij}\Psi(E_{ij})^*=\sum\bar a_{ij}\Psi(E_{ji})=\Psi(A^*)$. CP-отображение имеет форму Крауса
(теорема 6.2), и $(\sum_lK_lAK_l^*)^*=\sum_lK_lA^*K_l^*$. $\square$

<a id="r-D-2"></a>

**Лемма D.2 (разложение).** Пусть $q_1,\dots,q_C$ — ортонормированные векторы в
$\mathbb C^n\otimes\mathbb C^n$, записанные как $q_c=(1\otimes Q_c)\Omega$, пусть $W=\mathrm{span}\{q_c\}$, и пусть
$\Psi$ — сохраняющий эрмитовость супероператор, для которого $\langle v,J(\Psi)v\rangle\ge0$ при всех $v\perp W$.
Положим $J=J(\Psi)$ и $\xi_c=Jq_c-\frac12\sum_d\langle q_d,Jq_c\rangle q_d$, и определим $K_c\in M_n$ условием
$(1\otimes K_c)\Omega=\xi_c$. Тогда существуют $L_1,\dots,L_r\in M_n$, где $r=\mathrm{rank}(P_W^\perp JP_W^\perp)$,
а $P_W^\perp$ — проектор на $W^\perp$, такие, что
$$
\Psi(A)=\sum_{j=1}^rL_jAL_j^*+\sum_{c=1}^C\big(K_cAQ_c^*+Q_cAK_c^*\big)\qquad(A\in M_n),
$$
$\mathrm{Tr}(Q_c^*L_j)=0$ для всех $c,j$ и $\mathrm{Tr}(L_j^*L_k)=0$ при $j\neq k$.

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
$\mathrm{Tr}(Q_c^*L_j)=\langle q_c,l_j\rangle=0$ и $\mathrm{Tr}(L_j^*L_k)=\langle l_j,l_k\rangle=0$ при $j\ne k$,
а по лемме D.1(c) $P^\perp JP^\perp$ — матрица Чоя отображения $A\mapsto\sum_jL_jAL_j^*$. Для второго слагаемого:
$JP=\sum_c|Jq_c\rangle\langle q_c|$ и $PJP=\sum_{c,d}\langle q_d,Jq_c\rangle|q_d\rangle\langle q_c|$, поэтому
$$
\sum_c|\xi_c\rangle\langle q_c|=JP-\tfrac12PJP,\qquad\sum_c|q_c\rangle\langle\xi_c|=PJ-\tfrac12PJP
$$
(второе равенство сопряжено первому; $PJP$ самосопряжена), а их сумма равна $PJ+JP-PJP$. По
лемме D.1(c) $|\xi_c\rangle\langle q_c|$ и $|q_c\rangle\langle\xi_c|$ — матрицы Чоя отображений
$A\mapsto K_cAQ_c^*$ и $A\mapsto Q_cAK_c^*$. Итак, обе части доказываемого тождества имеют матрицу Чоя $J$,
и лемма D.1(a) даёт тождество. $\square$

<a id="s-D-2"></a>

## D.2 Условная вполне положительность

Супероператор $\Psi$ на $M_n$ называется **условно вполне положительным** (CCP), если
$\langle v,J(\Psi)v\rangle\ge0$ для всякого $v\perp\Omega$. По лемме D.1(b) эти $v$ — в точности векторы
$(1\otimes K)\Omega$ с $\mathrm{Tr}\,K=0$.

<a id="r-D-3"></a>

**Теорема D.3.** Для супероператора $\mathcal L$ на $M_n$ следующие условия равносильны:

1. $e^{t\mathcal L}$ является CP при всяком $t\ge0$;
2. $\mathcal L$ сохраняет эрмитовость и условно вполне положителен;
3. существуют CP-отображение $\Phi$ на $M_n$ и $G\in M_n$, для которых $\mathcal L(A)=\Phi(A)+GA+AG^*$ при всех $A$.

В п. 3 можно взять $\Phi(A)=\sum_{j=1}^rL_jAL_j^*$ с $r\le n^2-1$, $\mathrm{Tr}\,L_j=0$ и
$\mathrm{Tr}(L_j^*L_k)=0$ при $j\neq k$, а $G$ — с $\mathrm{Tr}\,G\in\mathbb R$.

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
(теорема 6.2). Далее, $\mathrm{Tr}\,L_j=\sqrt n\mathrm{Tr}(Q_1^*L_j)=0$, $r\le\mathrm{rank}\,P_W^\perp=n^2-1$,
и по лемме D.1(b)
$\mathrm{Tr}\,G=\langle\Omega,\xi_1\rangle/\sqrt n=\langle q_1,\xi_1\rangle=\frac12\langle q_1,J(\mathcal L)q_1\rangle\in\mathbb R$,
так как $J(\mathcal L)$ самосопряжена.

(3)$\Rightarrow$(1). Шаги 1–4 доказательства теоремы 8.6 используют лишь то, что генератор равен
$\mathcal L_0+\Phi$, где $\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)$ для некоторой матрицы
$H_{\mathrm{eff}}$, а $\Phi$ является CP. При $H_{\mathrm{eff}}=iG$ имеем $\mathcal L_0(\rho)=G\rho+\rho G^*$, и
эти шаги показывают, что $e^{t\mathcal L}=\sum_{m\ge0}S_m(t)$ является CP. $\square$

*Вычисление:* теорема D.3 на генераторе GKSL на $M_3$ и пример $\theta-\mathrm{id}$ ниже — [руководство, раздел D.2](../../tutorial/ru/D_gksl.md#s-D-2). *В библиотеке:* матрица генератора GKSL, `dynamics.lindblad_superoperator` — [руководство, раздел 8.2](../../tutorial/ru/08_dynamics.md#s-8-2).

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

<a id="s-D-3"></a>

## D.3 Теорема GKSL

<a id="r-D-4"></a>

**Теорема D.4 (Горини–Коссаковский–Сударшан; Линдблад; теорема 8.7).** Супероператор $\mathcal L$ на
$M_n$ является генератором квантовой динамической полугруппы на $M_n$ тогда и только тогда, когда существуют $H=H^*$ и
$L_1,\dots,L_r$ в $M_n$, для которых
$$
\mathcal L(\rho)=-i[H,\rho]+\sum_{j=1}^r\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
$$
Можно выбрать $r\le n^2-1$, $\mathrm{Tr}\,H=0$, $\mathrm{Tr}\,L_j=0$ и $\mathrm{Tr}(L_j^*L_k)=0$ при $j\neq k$.

*Доказательство.* Часть «если» — это теорема 8.6. Обратно, пусть $(T_t)$ — квантовая динамическая полугруппа с генератором
$\mathcal L$, так что $T_t=e^{t\mathcal L}$ (лемма 8.5). Поскольку каждое $T_t$ является CP, теорема D.3 даёт
$\mathcal L(\rho)=\sum_jL_j\rho L_j^*+G\rho+\rho G^*$ с $r\le n^2-1$ бесследовыми попарно ортогональными $L_j$ и
$\mathrm{Tr}\,G\in\mathbb R$. Положим $H=\frac i2(G-G^*)$ и $R=\frac12(G+G^*)$; обе матрицы самосопряжены, и
$G=-iH+R$, $G^*=iH+R$, $G\rho+\rho G^*=-i[H,\rho]+\{R,\rho\}$. Кроме того,
$\mathrm{Tr}\,H=\frac i2(\mathrm{Tr}\,G-\overline{\mathrm{Tr}\,G})=0$. Поскольку каждое $T_t$ является TP,
$\mathrm{Tr}\,\mathcal L(\rho)=\frac d{dt}\mathrm{Tr}\,T_t(\rho)\big|_{t=0}=0$ (лемма 8.1), и по цикличности
$$
0=\mathrm{Tr}\,\mathcal L(\rho)=\mathrm{Tr}\Big(\rho\,\Big(\sum_jL_j^*L_j+2R\Big)\Big)\qquad\text{for all }\rho .
$$
При $\rho=X=\sum_jL_j^*L_j+2R$, которая самосопряжена, получаем $\mathrm{Tr}(X^*X)=0$, поэтому $X=0$ и
$R=-\frac12\sum_jL_j^*L_j$, что и есть форма GKSL. $\square$

*Вычисление:* доказательство теоремы D.4, выполненное как алгоритм, восстанавливающий $H$ и операторы скачков по $\mathcal L$ — [руководство, раздел D.3](../../tutorial/ru/D_gksl.md#s-D-3).

<a id="r-D-5"></a>

**Теорема D.5 (свобода в форме GKSL).** Пусть $(H;L_1,\dots,L_r)$ и $(H';L'_1,\dots,L'_s)$ — данные
GKSL ($H,H'$ самосопряжены) с генераторами $\mathcal L$, $\mathcal L'$, как в теореме D.4.

(a) $\mathcal L'=\mathcal L$ в каждом из следующих случаев:
(i) $L'_j=L_j$ и $H'=H+c1$, $c\in\mathbb R$;
(ii) $H'=H$ и, после дополнения более короткого семейства нулями, $L'_j=\sum_ku_{jk}L_k$ с унитарной $u$;
(iii) $L'_j=L_j+a_j1$ с $a_j\in\mathbb C$ и $H'=H-\frac i2\sum_j(\bar a_jL_j-a_jL_j^*)$; эта матрица самосопряжена.

(b) Если $\mathrm{Tr}\,H=\mathrm{Tr}\,H'=0$, все $L_j,L'_j$ бесследовые и $\mathcal L=\mathcal L'$, то $H=H'$ и,
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
$g=(1\otimes G)\Omega$, причём $l_j\perp\Omega$, поскольку $\mathrm{Tr}\,L_j=0$ (лемма D.1(b)). Обозначая через $P^\perp$ проектор
на $\Omega^\perp$, получаем $P^\perp J(\mathcal L)P^\perp=\sum_j|l_j\rangle\langle l_j|=J(\Phi)$ и аналогично
$P^\perp J(\mathcal L')P^\perp=J(\Phi')$. Значит, $J(\Phi)=J(\Phi')$, $\Phi=\Phi'$ (лемма D.1(a)), и теорема 6.7 даёт
унитарную $u$. Тогда $\Phi^*(1)=\Phi'^*(1)$, откуда $[H-H',\rho]=0$ для всех $\rho$. Матрица $X$, коммутирующая со всеми
$E_{ij}$, кратна $1$ (сравните элементы $XE_{ij}$ и $E_{ij}X$), поэтому $H-H'=\lambda1$, и
$\lambda=0$, как видно из взятия следа.

(c) По (a) шаги (i)–(iii) не меняют генератор, и шаг, обратный к каждому шагу, — шаг того же
вида ((iii) с $-a_j$, так как $\bar a_jL'_j-a_jL_j'^*=\bar a_jL_j-a_jL_j^*$). Любые данные приводятся к
нормированной форме шагом (iii) с $a_j=-\mathrm{Tr}\,L_j/n$ и затем шагом (i) с $c=-\mathrm{Tr}\,H''/n$, где $H''$ —
гамильтониан после (iii); два нормированных набора данных с одним и тем же генератором связаны шагом (ii) в силу (b). $\square$

*Вычисление:* теорема D.5: нормированная форма через (iii) и (i) и единственность в (b) — [руководство, раздел D.3](../../tutorial/ru/D_gksl.md#s-D-3).

<a id="r-D-6"></a>

**Теорема D.6 (единственность стационарного состояния).** Пусть $\mathcal L$ — генератор GKSL из теоремы D.4
с данными $(H;L_1,\dots,L_r)$. Если у $\mathcal L$ есть точное стационарное состояние $\sigma$ ($\sigma>0$) и единственные
матрицы, коммутирующие с $H$ и всеми $L_j$, $L_j^*$, — это кратные $1$, то $\ker\mathcal L=\mathbb C\sigma$:
$\sigma$ — единственное стационарное состояние. Кроме того, $\ker\mathcal L^*=\mathbb C1$.

*Доказательство.* Сопряжённое относительно следового спаривания, $\mathrm{Tr}(\mathcal L(\rho)A)=\mathrm{Tr}(\rho\,\mathcal L^*(A))$, равно
$\mathcal L^*(A)=i[H,A]+\sum_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})$: по цикличности
$\mathrm{Tr}(-i[H,\rho]A)=\mathrm{Tr}(\rho\,i[H,A])$, $\mathrm{Tr}(L\rho L^*A)=\mathrm{Tr}(\rho L^*AL)$ и
$\mathrm{Tr}(\{L^*L,\rho\}A)=\mathrm{Tr}(\rho\{L^*L,A\})$ (это упражнение 1 главы 8).
Прямое вычисление, в котором гамильтоновы слагаемые сокращаются, поскольку $A\mapsto i[H,A]$ — дифференцирование, даёт
для всякого $X$
$$
\mathcal L^*(X^*X)-X^*\mathcal L^*(X)-\mathcal L^*(X^*)X=\sum_j[L_j,X]^*[L_j,X]:
$$
для одного скачка из шести слагаемых с антикоммутаторами четыре попарно сокращаются, а оставшиеся два,
$\frac12X^*L^*LX+\frac12X^*L^*LX$, складываются, так что левая часть равна $L^*X^*XL-X^*L^*XL-L^*X^*LX+X^*L^*LX$, т. е.
$(XL-LX)^*(XL-LX)$. Пусть $X\in\ker\mathcal L^*$. Тогда $\mathcal L^*(X^*)=\mathcal L^*(X)^*=0$, так как $\mathcal L^*$
сохраняет эрмитовость, поэтому $\mathcal L^*(X^*X)=\sum_j[L_j,X]^*[L_j,X]=:Y\ge0$ и
$\mathrm{Tr}(\sigma Y)=\mathrm{Tr}(\mathcal L(\sigma)X^*X)=0$. Поскольку $\sigma>0$, величина
$\mathrm{Tr}(\sigma Y)=\mathrm{Tr}(\sigma^{1/2}Y\sigma^{1/2})$ при $\sigma^{1/2}Y\sigma^{1/2}\ge0$ обращается в нуль, лишь если
$\sigma^{1/2}Y\sigma^{1/2}=0$, т. е. $Y=0$; взяв след, получаем $\sum_j\mathrm{Tr}([L_j,X]^*[L_j,X])=0$ с неотрицательными
слагаемыми, поэтому $[L_j,X]=0$ для всех $j$. То же рассуждение для $X^*\in\ker\mathcal L^*$
даёт $[L_j,X^*]=0$, т. е. $[L_j^*,X]=0$. Тогда $X$ коммутирует с $L_j^*L_j$, диссипативная часть $\mathcal L^*(X)$
обращается в нуль, и $0=\mathcal L^*(X)=i[H,X]$. По предположению $X\in\mathbb C1$; обратно, $\mathcal L^*(1)=0$. Итак,
$\ker\mathcal L^*=\mathbb C1$. В базисах $(E_{ij})$ и $(E_{ji})$ пространства $M_n$, двойственных относительно спаривания
$(\rho,A)\mapsto\mathrm{Tr}(\rho A)$, матрица $\mathcal L^*$ — транспонированная матрица $\mathcal L$; поэтому у них
одинаковый ранг, и $\dim\ker\mathcal L=\dim\ker\mathcal L^*=1$. Состояние из $\mathbb C\sigma$ со следом $1$ — это $\sigma$. $\square$

Это конечномерный случай критерия Фриджерио [*Commun. Math. Phys.* 63, 269 (1978)],
используемого в разделе 8.5. Оба предположения необходимы. В упражнении 4 главы 8 состояние $\tfrac12\,1$ точное
и стационарное, но $H$ и $\sigma_z$ коммутируют со всеми диагональными матрицами, и стационарные состояния порождают
двумерное пространство диагональных матриц. Обратно, пусть на $M_3$ $H=0$, $L_1=|e_1\rangle\langle e_3|$,
$L_2=|e_2\rangle\langle e_3|$. $*$-алгебра, порождённая $L_1,L_2$, содержит $L_1L_2^*=E_{12}$ и $L_1^*L_1=E_{33}$, а значит, и
все матричные единицы, поэтому с $L_j,L_j^*$ коммутируют лишь кратные $1$; но всякое $\rho=P\rho P$ с $P=E_{11}+E_{22}$
стационарно, так как $L_j\rho=0=\rho L_j^*$ и $\{E_{33},\rho\}=0$, поэтому $\dim\ker\mathcal L\ge4$; по теореме D.6 ни одно стационарное
состояние не является точным.

*Вычисление:* теорема D.6 для генератора из раздела D.2 и второй пример — [руководство, раздел D.3](../../tutorial/ru/D_gksl.md#s-D-3).

<a id="s-D-4"></a>

## D.4 Алгебры с секторами

Пусть $\mathcal M=\bigoplus_{c=1}^CM_{k_c}\subseteq M_D$ — блочно-диагональные матрицы, $D=\sum_ck_c$, с
минимальными центральными проекторами $P_c$, и пусть $E(X)=\sum_cP_cXP_c$ — условное ожидание $M_D$ на
$\mathcal M$. Как в доказательстве предложения 8.8, $E(B)=B$, $E(XB)=E(X)B$ и $E(BX)=BE(X)$ для $B\in\mathcal M$,
и $E(X^*)=E(X)^*$. Для CP-отображения $\Phi:\mathcal M\to\mathcal M$ через $\Phi^*$ обозначается сопряжённое к нему относительно `Tr_blunt` (предложение 6.8).

<a id="r-D-7"></a>

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
$W=\mathrm{span}\{p_c\}$. По лемме D.1(c) $J(E)=\sum_c|p_c\rangle\langle p_c|$, поэтому для $v\perp W$ и $t>0$
$0\le\frac1t\langle v,J(\hat T_t)v\rangle=\langle v,J((\hat T_t-E)/t)v\rangle\to\langle v,J(\mathcal L\circ E)v\rangle$.
Каждое $T_t$ сохраняет эрмитовость (в силу своей формы Крауса), значит, и $\mathcal L$ тоже (как в теореме D.3), как и $E$
и $\mathcal L\circ E$.

*Шаг 3: разложение.* По лемме D.1(b) $\langle p_c,p_d\rangle=\mathrm{Tr}(P_cP_d)=\delta_{cd}k_c$, поэтому
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

*Шаг 5: сохранение следа.* Как в теореме D.4, $\mathrm{Tr}\,\mathcal L(A)=0$ для $A\in\mathcal M$, и
$\mathrm{Tr}\,\Phi(A)=\mathrm{Tr}(A\Phi^*(1))$, поэтому $\mathrm{Tr}(AX)=0$ для всех $A\in\mathcal M$, где
$X=\Phi^*(1)+G+G^*\in\mathcal M$ самосопряжена. При $A=X$ получаем $X=0$ (след `Tr_blunt` точен, теорема 3.1).
При $H=\frac i2(G-G^*)\in\mathcal M$ имеем $G=-iH-\frac12\Phi^*(1)$ и $G\rho+\rho G^*=-i[H,\rho]-\frac12\{\Phi^*(1),\rho\}$. $\square$

Векторы $p_c$ заменяют $\Omega=\sum_cp_c$, потому что $\hat T_0=E$, а не $\mathrm{id}$: матрица $J(\mathcal L\circ E)$ не обязана быть
положительной на $\Omega^\perp$, и в примере ниже она таковой не является.

*Вычисление:* теорема D.7 на $M_1\oplus M_2\subseteq M_3$ по пяти шагам доказательства; там $J(\mathcal L\circ E)$ не положительна на $\Omega^\perp$ — [руководство, раздел D.4](../../tutorial/ru/D_gksl.md#s-D-4). *В библиотеке:* генераторы на алгебре с секторами, `dynamics.lindblad_rhs` с `channels.random_inter_sector_channel` — [руководство, раздел 8.3](../../tutorial/ru/08_dynamics.md#s-8-3).

<a id="r-D-8"></a>

**Следствие D.8 (картина Гейзенберга).** Пусть $(S_t)_{t\ge0}$ — непрерывное семейство унитальных CP-отображений на
$\mathcal M$ с $S_0=\mathrm{id}$ и $S_{s+t}=S_sS_t$. Его генератор имеет вид $A\mapsto\Psi(A)+K^*A+AK$, где
$\Psi$ является CP на $\mathcal M$, $K\in\mathcal M$ и $\Psi(1)+K+K^*=0$.

*Доказательство.* Пусть $T_t=S_t^*$ (сопряжённые относительно `Tr_blunt`). Они являются CP (предложение 6.8) и TP (следствие 6.9(b)),
$T_0=\mathrm{id}$ и $T_{s+t}=(S_tS_s)^*=S_s^*S_t^*=T_sT_t$, поскольку
$\mathrm{Tr}(S_tS_s(x)y)=\mathrm{Tr}(xS_s^*S_t^*(y))$ и $S_tS_s=S_{s+t}$. Так как $\Psi\mapsto\Psi^*$ линейно на
конечномерном пространстве, $t\mapsto T_t$ непрерывно и $(T_t-\mathrm{id})/t\to\mathcal L_S^*$, где $\mathcal L_S$ —
генератор $(S_t)$. Теорема D.7 даёт $\mathcal L_S^*(\rho)=\Phi(\rho)+G\rho+\rho G^*$ с $G+G^*=-\Phi^*(1)$, а по
цикличности $\mathrm{Tr}((G\rho+\rho G^*)A)=\mathrm{Tr}(\rho(G^*A+AG))$, поэтому по следствию 6.9(c)
$\mathcal L_S(A)=\Phi^*(A)+G^*A+AG$. Возьмём $\Psi=\Phi^*$ (CP по предложению 6.8) и $K=G$. $\square$

Это конечномерный случай теоремы Кристенсена и Эванса, процитированной после теоремы 8.7.

*Вычисление:* генератор из следствия D.8 как двойственный к $\mathcal L$ — [руководство, раздел D.4](../../tutorial/ru/D_gksl.md#s-D-4).

<a id="s-D-5"></a>

## D.5 Коррекция ошибок для операторных алгебр

Используем постановку раздела 10.3. $\Phi(\rho)=\sum_{i=1}^rE_i\rho E_i^*$ — канал на $M_D$, $\mathcal C=P\mathbb C^D$ —
код, $B(\mathcal C)=\{X\in M_D:X=PXP\}$ — операторы на $\mathcal C$, $\mathfrak A\subseteq B(\mathcal C)$ — $*$-алгебра
с единицей $P$ (по теореме 2.2 — алгебра фон Неймана на $\mathcal C$), а $\mathfrak A'$ — её коммутант в
$B(\mathcal C)$. $\mathfrak A$ **корректируема** для $\Phi$, если существует канал $\mathcal R$ на $M_D$, для которого
$P(\mathcal R\circ\Phi)^*(A)P=A$ при всех $A\in\mathfrak A$. Если у $\mathcal R$ операторы Крауса $R_m$, то
у $\mathcal R\circ\Phi$ операторы Крауса $R_mE_i$, и $\sum_mR_m^*R_m=1$.

<a id="r-D-9"></a>

**Теорема D.9 (Бени–Кемпф–Крибс; теорема 10.7).** $\mathfrak A$ корректируема для $\Phi$ тогда и только тогда, когда
$PE_i^*E_jP\in\mathfrak A'$ для всех $i,j$. В этом случае канал $\mathcal R$, построенный в доказательстве, корректирует $\mathfrak A$.

*Доказательство.* *Необходимость.* Пусть $\mathcal R$ корректирует $\mathfrak A$; положим $W_{mi}=R_mE_iP$, так что
$\sum_{m,i}W_{mi}^*XW_{mi}=P(\mathcal R\circ\Phi)^*(X)P$ для всех $X\in M_D$. Это равно $X$ при $X\in\mathfrak A$ и $P$ при
$X=1$, поскольку $(\mathcal R\circ\Phi)^*$ унитально (следствие 6.9(a)). Пусть $p\in\mathfrak A$ — проектор; тогда $p\le P$
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
$\Lambda^c=U^c\mathrm{diag}(d^c_\kappa)U^{c*}$ с унитарной $U^c$ и $d^c_\kappa\ge0$, и положим
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

*Вычисление:* теорема D.9 и восстановление из её доказательства на гибридном коде в $\mathbb C^8$ — [руководство, раздел D.5](../../tutorial/ru/D_gksl.md#s-D-5).

<a id="r-D-10"></a>

**Следствие D.10.** (a) Для $\mathfrak A=B(\mathcal C)$ имеем $\mathfrak A'=\mathbb CP$, и теорема D.9 — это теорема
Нилла–Лафламма 10.5 ($C=1$, $m_1=1$; шаг 3 — восстановление из теоремы 10.5). (b) Для $\mathfrak A=\mathrm{span}\{P_1,\dots,P_C\}$ с
попарно ортогональными проекторами $P_c\ne0$, $\sum_cP_c=P$, это равносильность условий 2 и 3 в предложении 10.6. (c) (*Бесшумные
подсистемы.*) В обозначениях предложения 10.1(a) пусть $\mathcal C=V(\mathbb C^{k_c}\otimes\mathbb C^{m_c})$ — $c$-е
слагаемое для алгебры взаимодействия $\mathcal A$ канала $\Phi$, и $\mathfrak A=V(1_{k_c}\otimes M_{m_c})V^*$ на $\mathcal C$ —
алгебра вида $M_{m}\otimes1_{k}$ после перестановки двух сомножителей. Тогда $\mathfrak A$ корректируется
каналом $\mathcal R=\mathrm{id}$, и $PE_i^*E_jP=V(E_{i,c}^*E_{j,c}\otimes1_{m_c})V^*\in\mathfrak A'$ ($P=VV^*$).

*Доказательство.* (a) $\mathfrak A'=\mathbb CP$ — это случай $C=1$, $m_1=1$ приведённого выше описания $\mathfrak A'$, а
$PE_i^*E_jP\in\mathbb CP$ — условие Нилла–Лафламма. Корректируемость $B(\mathcal C)$ — это коррекция в смысле
раздела 10.2: если $\mathcal R(\Phi(\rho))=\rho$ для всех состояний $\rho$, сосредоточенных в $\mathcal C$, то
$\mathrm{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(A)P)=\mathrm{Tr}(\rho A)$ для всех них, а эти $\rho$ порождают $B(\mathcal C)$;
обратное доказано в разделе 10.3. При $C=1$, $m_1=1$ операторы Крауса $\iota_1V^{1*}_\kappa$ из шага 3 совпадают с операторами
из теоремы 10.5, поскольку $\iota_1\iota_1^*=P$. (b) $T\in B(\mathcal C)$ коммутирует со всеми $P_c$ тогда и только тогда, когда
$P_cTP_d=0$ при $c\ne d$, а это условие 2 предложения 10.6. (c) Всякий $A\in\mathfrak A$, продолженный нулём на
остальные слагаемые, лежит в $\mathcal A'=V(\bigoplus_c1_{k_c}\otimes M_{m_c})V^*$ (теорема 2.6), поэтому $\Phi^*(A)=A$ по предложению
10.1(a), и $PAP=A$. Формула для $PE_i^*E_jP$ следует из $V^*E_iV=\bigoplus_cE_{i,c}\otimes1_{m_c}$, а
$E_{i,c}^*E_{j,c}\otimes1_{m_c}$ коммутирует с $1_{k_c}\otimes M_{m_c}$. $\square$

По теореме 2.6 всякая $\mathfrak A$ $*$-изоморфна некоторой $\bigoplus_cM_{k_c}$. Когда $C\ge2$ и некоторое $k_c\ge2$,
корректируемая $\mathfrak A$ защищает метку сектора $c$ (классическую информацию) вместе с квантовой информацией
каждой $M_{k_c}$; такие коды называются **гибридными**.

Результаты этого приложения принадлежат V. Gorini, A. Kossakowski и E. C. G. Sudarshan
[*J. Math. Phys.* 17, 821 (1976)] и G. Lindblad [*Commun. Math. Phys.* 48, 119 (1976)] (теоремы D.3–D.5),
E. Christensen и D. E. Evans [*J. London Math. Soc.* 20, 358 (1979)] (следствие D.8 в произвольных алгебрах
фон Неймана), A. Frigerio (теорема D.6) и C. Bény, A. Kempf и D. W. Kribs [*Phys. Rev. Lett.* 98, 100502
(2007); *Phys. Rev. A* 76, 042303 (2007)] (теорема D.9); приведённые выше доказательства полны для конечномерных
утверждений, используемых в курсе.
