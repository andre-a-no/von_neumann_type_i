[English](../B_schur_weyl.md) · [Español](../es/B_schur_weyl.md) · [Français](../fr/B_schur_weyl.md) · [Deutsch](../de/B_schur_weyl.md) · [中文](../zh/B_schur_weyl.md) · [日本語](../ja/B_schur_weyl.md) · **Русский**

# Приложение B. Двойственность Шура–Вейля и исчисление Вайнгартена

Унитарная группа $U(d)$ действует на $(\mathbb C^d)^{\otimes n}$ посредством $V\mapsto V^{\otimes n}$, а симметрическая
группа $S_n$ — перестановками тензорных сомножителей. Каждое из двух семейств операторов является коммутантом
другого (**двойственность Шура–Вейля**). В этом приложении двойственность доказывается полностью, и из неё выводится
формула Вайнгартена для моментов унитарной матрицы Хаара любого порядка. Попутно доказываются факты, которые главы
4, 5, 9 и 11 берут отсюда: теорема 4.14(b) для всех $L$, неприводимость представлений $S_L$ на пространствах
кратности полного спина, числа Каталана из раздела 4.5, теорема 5.12, теорема Вигнера–Эккарта для $SU(2)$,
используемая в разделе 9.3, и мера Хаара ортогональной группы, используемая в главе 11.

Путь элементарен: тождество поляризации отождествляет коммутант перестановок с линейной оболочкой всех
$A^{\otimes n}$, $A\in M_d(\mathbb C)$; рассуждение с преобразованием Кэли показывает, что унитарные $V^{\otimes n}$ порождают
ту же линейную оболочку; теорема о бикоммутанте (теорема 2.2) превращает это в двойственность. Формула Вайнгартена
затем следует из того, что усреднение по Хаару есть ортогональная проекция на коммутант.

Всюду $d,n\ge1$, $\mathcal H_n=(\mathbb C^d)^{\otimes n}$ с ортонормированным базисом
$e_{\mathbf i}=e_{i_1}\otimes\dots\otimes e_{i_n}$, где мультииндекс $\mathbf i=(i_1,\dots,i_n)\in[d]^n$ понимается как отображение
$[n]\to[d]$, $m\mapsto i_m$ ($[d]=\{1,\dots,d\}$). На $\operatorname{End}\mathcal H_n=M_{d^n}(\mathbb C)$ задано
скалярное произведение Гильберта–Шмидта $\langle X,Y\rangle=\operatorname{Tr}(X^*Y)$. Для $X\in M_d$ и $1\le m\le n$
$X^{(m)}=1\otimes\dots\otimes X\otimes\dots\otimes1$ ($X$ на позиции $m$) и $d\Pi(X)=\sum_{m=1}^nX^{(m)}$. Обозначим
$$
\mathcal T_n=\operatorname{span}\{A^{\otimes n}:A\in M_d(\mathbb C)\},\qquad
\mathcal U_n=\operatorname{span}\{V^{\otimes n}:V\in U(d)\}.
$$
Произведения кронекеровых произведений вычисляются посомножительно, $(A_1\otimes\dots\otimes A_n)(B_1\otimes\dots\otimes B_n)=A_1B_1\otimes\dots\otimes A_nB_n$,
и $(A_1\otimes\dots\otimes A_n)^*=A_1^*\otimes\dots\otimes A_n^*$. Перестановки компонуются как отображения,
$(\sigma\tau)(m)=\sigma(\tau(m))$, а $\#\text{cycles}(\sigma)$ — число циклов $\sigma$, включая неподвижные точки.

## B.1 Перестановки тензорных сомножителей

Для $\sigma\in S_n$ положим $\pi(\sigma)e_{\mathbf i}=e_{\mathbf i\circ\sigma^{-1}}$, т. е. базисный вектор из позиции $m$ переходит в позицию
$\sigma(m)$. Это соглашение `decompose.permutation_representation` и глав 2 и 4.

**Лемма B.1.** (a) $\pi$ — унитарное представление $S_n$, и
$\pi(\sigma)(x_1\otimes\dots\otimes x_n)=x_{\sigma^{-1}(1)}\otimes\dots\otimes x_{\sigma^{-1}(n)}$ для всех $x_m\in\mathbb C^d$.
(b) $\pi(\sigma)(A_1\otimes\dots\otimes A_n)\pi(\sigma)^*=A_{\sigma^{-1}(1)}\otimes\dots\otimes A_{\sigma^{-1}(n)}$ для всех
$A_m\in M_d$. В частности, $\pi(\sigma)$ коммутирует с каждым $A^{\otimes n}$ и каждым $d\Pi(X)$, и
$\pi(\sigma)X^{(m)}\pi(\sigma)^*=X^{(\sigma(m))}$.
(c) $\operatorname{Tr}\pi(\sigma)=d^{\#\text{cycles}(\sigma)}$ и $\langle\pi(\sigma),\pi(\tau)\rangle=d^{\#\text{cycles}(\sigma^{-1}\tau)}$.
Число циклов инвариантно относительно сопряжения и обращения, и две перестановки сопряжены тогда и только тогда,
когда у них одинаковая цикловая структура.
(d) $n!$ операторов $\pi(\sigma)$ линейно независимы тогда и только тогда, когда $n\le d$.

*Доказательство.* (a) $\mathbf i\mapsto\mathbf i\circ\sigma^{-1}$ — биекция $[d]^n$, поэтому $\pi(\sigma)$ переставляет векторы
ортонормированного базиса и унитарен. $\pi(\sigma)\pi(\tau)e_{\mathbf i}=\pi(\sigma)e_{\mathbf i\circ\tau^{-1}}=e_{\mathbf i\circ\tau^{-1}\circ\sigma^{-1}}=\pi(\sigma\tau)e_{\mathbf i}$.
Разлагая $x_m=\sum_ax_m(a)e_a$, получаем $x_1\otimes\dots\otimes x_n=\sum_{\mathbf i}\prod_mx_m(i_m)\,e_{\mathbf i}$; применение $\pi(\sigma)$ и
подстановка $\mathbf j=\mathbf i\circ\sigma^{-1}$ дают $\sum_{\mathbf j}\prod_mx_m(j_{\sigma(m)})\,e_{\mathbf j}=\sum_{\mathbf j}\prod_{p}x_{\sigma^{-1}(p)}(j_p)\,e_{\mathbf j}$,
а это и есть заявленный вектор-произведение.
(b) По (a), на векторе-произведении: $\pi(\sigma)^*=\pi(\sigma^{-1})$ помещает $x_{\sigma(m)}$ в позицию $m$, оператор применяет
там $A_m$, а $\pi(\sigma)$ переносит позицию $m$ в позицию $\sigma(m)$; поэтому в позиции $p$ в итоге оказывается
$A_{\sigma^{-1}(p)}x_p$. Векторы-произведения порождают $\mathcal H_n$. Для $X^{(m)}$ сомножитель $X$ переходит из позиции $m$
в позицию $\sigma(m)$, так что $d\Pi(X)=\sum_mX^{(m)}$ остаётся на месте.
(c) $\pi(\sigma)$ — матрица перестановки, и $e_{\mathbf i}$ неподвижен тогда и только тогда, когда $i_{\sigma^{-1}(p)}=i_p$ для всех
$p$, т. е. когда $\mathbf i$ постоянен на циклах $\sigma$: таких $\mathbf i$ ровно $d^{\#\text{cycles}(\sigma)}$. Далее,
$\langle\pi(\sigma),\pi(\tau)\rangle=\operatorname{Tr}(\pi(\sigma^{-1})\pi(\tau))=\operatorname{Tr}\pi(\sigma^{-1}\tau)$. Циклы
$\kappa\rho\kappa^{-1}$ — образы циклов $\rho$ под действием $\kappa$, а $\rho^{-1}$ проходит те же циклы в обратном
направлении. Обратно, если у $\rho,\rho'$ одинаковая цикловая структура, запишем их циклы друг под другом,
упорядочив по длине, и пусть $\kappa$ переводит каждый элемент первой строки в элемент под ним; тогда
$\kappa\rho\kappa^{-1}=\rho'$.
(d) Если $n\le d$, положим $\mathbf u=(1,2,\dots,n)$. Векторы $\pi(\sigma)e_{\mathbf u}=e_{\mathbf u\circ\sigma^{-1}}$ — различные базисные
векторы для различных $\sigma$, поэтому применение $\sum_\sigma c_\sigma\pi(\sigma)=0$ к $e_{\mathbf u}$ даёт все $c_\sigma=0$. Если $n>d$,
положим $\mathcal A=\sum_\sigma\operatorname{sgn}(\sigma)\pi(\sigma)$. У каждого $\mathbf i\in[d]^n$ есть $i_a=i_b$ для некоторых $a\neq b$; для
транспозиции $t=(a\,b)$ имеем $\pi(t)e_{\mathbf i}=e_{\mathbf i}$ и $\mathcal A\pi(t)=\sum_\sigma\operatorname{sgn}(\sigma)\pi(\sigma t)=-\mathcal A$.
Следовательно, $\mathcal Ae_{\mathbf i}=\mathcal A\pi(t)e_{\mathbf i}=-\mathcal Ae_{\mathbf i}=0$ для всех $\mathbf i$: $\mathcal A=0$ — нетривиальное
линейное соотношение. $\square$

## B.2 Поляризация: коммутант перестановок

**Лемма B.2 (поляризация).** Для $A_1,\dots,A_n\in M_d(\mathbb C)$ при $|\varepsilon|=\sum_m\varepsilon_m$
$$
\sum_{\sigma\in S_n}A_{\sigma(1)}\otimes\dots\otimes A_{\sigma(n)}=
\sum_{\varepsilon\in\{0,1\}^n}(-1)^{n-|\varepsilon|}\Big(\sum_{m=1}^n\varepsilon_mA_m\Big)^{\otimes n}.
$$

*Доказательство.* По полилинейности $\big(\sum_m\varepsilon_mA_m\big)^{\otimes n}=\sum_f\varepsilon_{f(1)}\cdots\varepsilon_{f(n)}\,A_{f(1)}\otimes\dots\otimes A_{f(n)}$,
где сумма берётся по всем отображениям $f:[n]\to[n]$. Поскольку $\varepsilon_m^r=\varepsilon_m$ при $r\ge1$, коэффициент при
слагаемом $f$ в правой части равен
$$
\sum_{\varepsilon}\prod_{m=1}^n(-1)^{1-\varepsilon_m}\prod_{m\in f([n])}\varepsilon_m
=\prod_{m\notin f([n])}\Big(\sum_{\varepsilon_m=0,1}(-1)^{1-\varepsilon_m}\Big)\prod_{m\in f([n])}\Big(\sum_{\varepsilon_m=0,1}(-1)^{1-\varepsilon_m}\varepsilon_m\Big).
$$
Множители первого рода равны $-1+1=0$, второго — $1$. Поэтому коэффициент равен $1$, если $f$ сюръективно, т. е.
является перестановкой, и $0$ иначе. $\square$

**Симметризатор** — это $\mathcal S(X)=\frac1{n!}\sum_\sigma\pi(\sigma)X\pi(\sigma)^*$ на $\operatorname{End}\mathcal H_n$.

**Предложение B.3.** $\mathcal S$ — ортогональная проекция $\operatorname{End}\mathcal H_n$ на $\pi(S_n)'$, и
$$
\pi(S_n)'=\mathcal T_n,\qquad\dim\pi(S_n)'=\binom{d^2+n-1}n .
$$

*Доказательство.* $\mathcal S$ — твирлинг конечной группы $S_n$ в представлении $\pi$, а значит, ортогональная
проекция на $\pi(S_n)'$ (теорема 12.5 с $\pi_1=\pi_2=\pi$; см. также раздел 4.2). По лемме B.1(b)
$\mathcal T_n\subseteq\pi(S_n)'$. Обратно, $\operatorname{End}\mathcal H_n$ порождается произведениями
$A_1\otimes\dots\otimes A_n$ матричных единиц, и по лемме B.1(b) и подстановке $\sigma\to\sigma^{-1}$
$\mathcal S(A_1\otimes\dots\otimes A_n)=\frac1{n!}\sum_\sigma A_{\sigma(1)}\otimes\dots\otimes A_{\sigma(n)}$, что лежит в $\mathcal T_n$ по
лемме B.2. Итак, $\pi(S_n)'=\mathcal S(\operatorname{End}\mathcal H_n)\subseteq\mathcal T_n$.

*Размерность.* Пусть $P=[d]\times[d]$ и $E_\alpha=E_{\alpha_1}\otimes\dots\otimes E_{\alpha_n}$ для $\alpha\in P^n$ с матричными единицами
$E_{(a,b)}=|e_a\rangle\langle e_b|$: это базис $\operatorname{End}\mathcal H_n$. По лемме B.1(b)
$\mathcal S(E_\alpha)=\frac1{n!}\sum_\sigma E_{\alpha\circ\sigma^{-1}}$, что равно $|O_\alpha|^{-1}\sum_{\beta\in O_\alpha}E_\beta$ для орбиты
$O_\alpha=\{\alpha\circ\sigma\}$ (каждый элемент орбиты встречается при одном и том же числе $n!/|O_\alpha|$ перестановок $\sigma$).
Суммы по разным орбитам содержат непересекающиеся множества базисных элементов, поэтому образ $\mathcal S$ имеет
базис, занумерованный орбитами. Два набора лежат в одной орбите тогда и только тогда, когда каждый элемент $P$
встречается в них одинаково часто: орбиты — это мультимножества размера $n$ из $N=d^2$ элементов $P$. Если
упорядочить $P$ как $1,\dots,N$, такое мультимножество есть последовательность $1\le c_1\le\dots\le c_n\le N$, а
$c_m\mapsto c_m+m-1$ — биекция на строго возрастающие последовательности в $[N+n-1]$, т. е. на $n$-элементные
подмножества; их $\binom{N+n-1}n$. $\square$

## B.3 Двойственность Шура–Вейля

**Лемма B.4 (многочлены, обращающиеся в нуль на унитарной группе).** Пусть $p:M_d(\mathbb C)\to\mathbb C$ — многочлен от
$d^2$ элементов матрицы, однородный степени $n\ge1$. Если $p(V)=0$ для всех $V\in U(d)$, то $p=0$.

*Доказательство.* *Шаг 1 (преобразование Кэли).* Для $A\in M_d$ положим $q(A)=p\big((1+iA)\operatorname{adj}(1-iA)\big)$, где
$\operatorname{adj}$ — присоединённая (союзная) матрица; $q$ — многочлен от элементов $A$. Пусть $H$ эрмитова. Её
собственные значения вещественны (следствие 1.4), поэтому $1-iH$ обратима, $\operatorname{adj}(1-iH)=\det(1-iH)(1-iH)^{-1}$
(правило Крамера), и в силу однородности $q(H)=\det(1-iH)^n\,p(C)$ с $C=(1+iH)(1-iH)^{-1}$. Матрицы $1\pm iH$ и обратные к
ним являются функциями от $H$ и коммутируют (теорема 1.5), а $(1+iH)^*=1-iH$, поэтому
$C^*C=(1+iH)^{-1}(1-iH)(1+iH)(1-iH)^{-1}=1$. Следовательно, $q(H)=0$.

*Шаг 2 (от эрмитовых матриц ко всем).* Выберем базис $B_1,\dots,B_{d^2}$ вещественного векторного пространства эрмитовых
матриц. Он является также базисом $M_d(\mathbb C)$ над $\mathbb C$, поскольку всякая $A$ равна
$\frac{A+A^*}2+i\frac{A-A^*}{2i}$, а размерности совпадают. $r(t)=q(\sum_kt_kB_k)$ — многочлен от $t\in\mathbb C^{d^2}$,
обращающийся в нуль на $\mathbb R^{d^2}$. Все коэффициенты такого многочлена равны $0$: от одной переменной у него
бесконечно много корней; для $N$ переменных запишем $r=\sum_jr_j(t')t_N^j$ — при каждом вещественном $t'$ многочлен
от $t_N$ обращается в нуль на $\mathbb R$, поэтому каждый $r_j$ обращается в нуль на $\mathbb R^{N-1}$, и применима индукция.
Значит, $q(A)=r(t)=0$ для всякой $A=\sum_kt_kB_k\in M_d$.

*Шаг 3 (обращение преобразования).* Пусть $B\in M_d$ и $1+B$ обратима, $A=-i(1+B)^{-1}(B-1)$. Тогда
$1-iA=1-(1+B)^{-1}(B-1)=2(1+B)^{-1}$ и $1+iA=(1+B)^{-1}2B$; так как $B$ коммутирует с $1+B$,
$(1+iA)(1-iA)^{-1}=(1+B)^{-1}B(1+B)=B$. Поэтому $0=q(A)=\det(1-iA)^np(B)$ с $\det(1-iA)\neq0$: $p(B)=0$.

*Шаг 4.* Для любой $B\neq0$ и $0<t<1/\|B\|$ число $-1$ не является собственным значением $tB$ (если $tBv=-v$, то
$\|v\|=t\|Bv\|\le t\|B\|\|v\|<\|v\|$, кроме случая $v=0$), поэтому $t^np(B)=p(tB)=0$. Наконец, $p(0)=0$ в силу
однородности. $\square$

**Лемма B.5 (тождества Ньютона).** Пусть $x_1,\dots,x_n$ — попарно коммутирующие элементы алгебры, $x_S=\prod_{m\in S}x_m$,
$e_k=\sum_{|S|=k}x_S$ ($e_0=1$) и $p_j=\sum_mx_m^j$. Тогда $k\,e_k=\sum_{j=1}^k(-1)^{j-1}e_{k-j}\,p_j$ при $1\le k\le n$.

*Доказательство.* Для $1\le j\le k$ положим $a_j=\sum x_Sx_m^j$, где сумма берётся по парам $(S,m)$ с $|S|=k-j$ и $m\notin S$,
и $a_{k+1}=0$. В $e_{k-j}p_j=\sum_{|S|=k-j}\sum_mx_Sx_m^j$ слагаемые с $m\notin S$ дают $a_j$, а слагаемые с $m\in S$,
записанные через $S'=S\setminus\{m\}$, равны $x_{S'}x_m^{j+1}$ и дают $a_{j+1}$. Поэтому знакопеременная сумма
телескопически сворачивается в $a_1$, а $a_1=\sum_{|T|=k}|T|\,x_T=k\,e_k$ (каждое $k$-элементное $T$ возникает при $k$
выборах $m\in T$). $\square$

**Теорема B.6 (двойственность Шура–Вейля).** Пусть $d,n\ge1$.
(a) $\pi(S_n)'=\mathcal T_n=\mathcal U_n=\operatorname{span}\{V^{\otimes n}:V\in SU(d)\}$. Это пространство — $*$-алгебра, а значит, это
алгебра, порождённая $V^{\otimes n}$, $V\in U(d)$.
(b) Для всякого множества $\mathcal X\subseteq M_d$, линейная оболочка которого есть $M_d$, алгебра с единицей, порождённая
$\{d\Pi(X):X\in\mathcal X\}$, есть $\mathcal T_n$.
(c) $\{V^{\otimes n}:V\in U(d)\}'=\{V^{\otimes n}:V\in SU(d)\}'=\{d\Pi(X):X\in M_d\}'=\operatorname{span}\pi(S_n)$.

*Доказательство.* (a) $\pi(S_n)'=\mathcal T_n$ — это предложение B.3, и $\mathcal U_n\subseteq\mathcal T_n$. Предположим, что
$\mathcal U_n\neq\mathcal T_n$. Тогда некоторое $A_0^{\otimes n}\notin\mathcal U_n$; пусть $Y\neq0$ — компонента $A_0^{\otimes n}$,
ортогональная $\mathcal U_n$. Функция $p(A)=\langle Y,A^{\otimes n}\rangle$ — однородный многочлен степени $n$ от элементов $A$
(каждый элемент $A^{\otimes n}$ есть произведение $n$ элементов $A$), она обращается в нуль на $U(d)$, а
$p(A_0)=\|Y\|^2\neq0$, что противоречит лемме B.4.
Для $V\in U(d)$ выберем $\theta$ с $\det V=e^{i\theta}$; тогда $W=e^{-i\theta/d}V\in SU(d)$ и $V^{\otimes n}=e^{in\theta/d}W^{\otimes n}$,
так что две линейные оболочки совпадают. $\mathcal U_n$ содержит $1=1^{\otimes n}$, $V^{\otimes n}W^{\otimes n}=(VW)^{\otimes n}$ и
$(V^{\otimes n})^*=(V^*)^{\otimes n}$, поэтому является $*$-алгеброй, а всякая алгебра, содержащая $V^{\otimes n}$, содержит их
линейную оболочку.

(b) Пусть $\mathcal A$ — эта алгебра. $d\Pi$ линейно, поэтому $\mathcal A$ содержит $d\Pi(X)$ для всех $X\in M_d$. По лемме B.1(b)
каждое $d\Pi(X)$ лежит в алгебре с единицей $\pi(S_n)'=\mathcal T_n$, так что $\mathcal A\subseteq\mathcal T_n$. Обратно, зафиксируем
$A\in M_d$ и применим лемму B.5 к попарно коммутирующим $x_m=A^{(m)}$: тогда $e_n=A^{\otimes n}$ и $p_j=d\Pi(A^j)\in\mathcal A$, и
по индукции по $k$ каждое $e_k=\frac1k\sum_j(-1)^{j-1}e_{k-j}p_j$ лежит в $\mathcal A$. Поэтому $A^{\otimes n}\in\mathcal A$ и
$\mathcal T_n\subseteq\mathcal A$.

(c) Коммутант множества совпадает с коммутантом порождённой им алгебры. По (a) и (b) все три множества
порождают $\mathcal T_n=\pi(S_n)'$, коммутант которой есть $\pi(S_n)''$. Линейная оболочка $\pi(S_n)$ — $*$-алгебра
($\pi(\sigma)\pi(\tau)=\pi(\sigma\tau)$, $\pi(\sigma)^*=\pi(\sigma^{-1})$, $\pi(e)=1$), поэтому $\pi(S_n)''=\operatorname{span}\pi(S_n)$ по
теореме 2.2. $\square$

**Следствие B.7 (теорема 4.14(b)).** Для полного спина $S_a=\sum_is_a^{(i)}$ системы $L$ спинов $\frac12$ алгебра, порождённая
$S_x,S_y,S_z$, есть $\{S_a\}''=\pi(S_L)'$, и $\{S_a\}'=\operatorname{span}\pi(S_L)$.

*Доказательство.* $M_2=\operatorname{span}\{1,s_x,s_y,s_z\}$ (матрицы Паули и $1$ образуют базис), $d\Pi(1)=L\cdot1$ и $d\Pi(s_a)=S_a$.
По теореме B.6(b) алгебра с единицей $\mathcal A$, порождённая $S_x,S_y,S_z$, есть $\mathcal T_L=\pi(S_L)'$. Поскольку $S_a$
эрмитовы, $\mathcal A$ — $*$-алгебра (сопряжённое к произведению образующих есть произведение в обратном порядке), поэтому
$\mathcal A=\mathcal A''\supseteq\{S_a\}''$ по теореме 2.2, и $\{S_a\}''\supseteq\mathcal A$, поскольку $\{S_a\}''$ — алгебра с единицей,
содержащая $S_a$. Наконец, $\{S_a\}'=\{S_a\}'''=\pi(S_L)''=\operatorname{span}\pi(S_L)$ (раздел 2.1 и теорема B.6(c)). $\square$

По теореме B.6(a) при $d=2$ та же алгебра является линейной оболочкой поворотов $V^{\otimes L}$, $V\in SU(2)$; для этого не
нужна экспоненциальная параметризация $SU(2)$, упомянутая в разделе 4.5.

Код сравнивает для малых $d$ и $n$ четыре размерности: $\operatorname{span}\pi(S_n)$, коммутанта двух случайных
$V^{\otimes n}$ (он содержит $\{V^{\otimes n}\}'=\operatorname{span}\pi(S_n)$ и совпадает с ним во всех напечатанных случаях) и
линейных оболочек множества случайных $A^{\otimes n}$ и $V^{\otimes n}$ — с $\binom{d^2+n-1}n$. Он также проверяет, что
`permutation_representation` библиотеки следует соглашению леммы B.1(a), лемму B.2 и рекурсию Ньютона из доказательства
теоремы B.6(b).

```python
import functools
import itertools
import math
import torch
import torch_vn_algebra.decompose as dc

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    return (torch.as_tensor(X) - torch.as_tensor(Y)).abs().max().item() < tol

def haar(d, B):                                        # Theorem 5.3: QR of Ginibre matrices, phases fixed
    Q, R = torch.linalg.qr(torch.randn(B, d, d, dtype=dt))
    r = R.diagonal(dim1=-2, dim2=-1)
    return (Q * (r / r.abs()).unsqueeze(-2)).contiguous()

def power(A, n):                                       # A^{(x) n} for a batch A of shape (B, d, d)
    out = A
    for _ in range(n - 1):
        out = torch.einsum('bij,bkl->bikjl', out, A).reshape(A.shape[0], out.shape[1] * A.shape[1], -1)
    return out

def rank(X, tol=1e-9):                                 # dimension of the span of a batch of matrices
    return torch.linalg.matrix_rank(X.reshape(X.shape[0], -1), atol=tol).item()

def commutant_dim(ops, tol=1e-9):                      # dim {X : XA = AX for all A in ops}
    eye = torch.eye(ops.shape[-1], dtype=dt)
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    return int((torch.linalg.svdvals(L) < tol).sum())

for d, n in ((2, 2), (2, 3), (2, 4), (3, 2), (3, 3)):
    P = dc.permutation_representation(n, d)            # the n! operators pi(sigma)
    N_sym = math.comb(d * d + n - 1, n)                # Proposition B.3
    span_A = rank(power(torch.randn(N_sym + 5, d, d, dtype=dt), n))
    span_V = rank(power(haar(d, N_sym + 5), n))
    comm_V = commutant_dim(power(haar(d, 2), n))       # commutant of two random V^{(x) n}
    print(d, n, rank(P), comm_V, span_A, span_V, N_sym)
    assert span_A == span_V == N_sym == commutant_dim(P)                     # Theorem B.6(a)
    assert comm_V == rank(P)              # numerical observation: two V suffice here (B.6(c) gives >=)
    assert rank(P) == (math.factorial(n) if n <= d else math.comb(2 * n, n) // (n + 1))   # B.1(d), B.10

A = torch.randn(3, 2, 2, dtype=dt)                     # Lemma B.2 for n = 3
lhs = sum(functools.reduce(torch.kron, [A[s[0]], A[s[1]], A[s[2]]]) for s in itertools.permutations(range(3)))
rhs = sum((-1) ** (3 - sum(eps)) * functools.reduce(torch.kron, [sum(e * a for e, a in zip(eps, A))] * 3)
          for eps in itertools.product((0, 1), repeat=3))
assert close(lhs, rhs)

xs = torch.randn(3, 2, dtype=dt)                       # Lemma B.1(a): pi(s) x_1 (x) x_2 (x) x_3 = (x) x_{s^-1(p)}
for s in itertools.permutations(range(3)):
    inv = [s.index(m) for m in range(3)]
    Ps = dc.permutation_representation(3, 2, [list(s)])[0]
    assert close(Ps @ functools.reduce(torch.kron, list(xs)), functools.reduce(torch.kron, [xs[inv[p]] for p in range(3)]))

def d_pi(X, n):                                        # sum_m 1 (x) .. (x) X (x) .. (x) 1
    eye = torch.eye(X.shape[0], dtype=dt)
    return sum(functools.reduce(torch.kron, [X if l == m else eye for l in range(n)]) for m in range(n))

n, A = 4, torch.randn(2, 2, dtype=dt)
e = [torch.eye(2 ** n, dtype=dt)]
for k in range(1, n + 1):                              # Lemma B.5 with x_m = A^(m), p_j = dPi(A^j)
    e.append(sum((-1) ** (j - 1) * e[k - j] @ d_pi(torch.linalg.matrix_power(A, j), n)
                 for j in range(1, k + 1)) / k)
assert close(e[n], functools.reduce(torch.kron, [A] * n), 1e-8)
```

## B.4 Разложение $(\mathbb C^d)^{\otimes n}$

**Теорема B.8.** Существуют $C\ge1$, целые $k_c,m_c\ge1$ с $\sum_ck_cm_c=d^n$ и унитарный оператор
$W:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathcal H_n$ такие, что
$$
W^*\pi(\sigma)W=\bigoplus_{c=1}^C\rho_c(\sigma)\otimes1_{m_c},\qquad W^*V^{\otimes n}W=\bigoplus_{c=1}^C1_{k_c}\otimes\tau_c(V)
\qquad(\sigma\in S_n,\ V\in U(d)),
$$
где $\rho_c$ — попарно неэквивалентные неприводимые унитарные представления $S_n$ на $\mathbb C^{k_c}$, а
$\tau_c$ — попарно неэквивалентные неприводимые унитарные представления $U(d)$ на $\mathbb C^{m_c}$. Кроме того,
$W^*\operatorname{span}\pi(S_n)W=\bigoplus_cM_{k_c}\otimes1_{m_c}$, $W^*\mathcal U_nW=\bigoplus_c1_{k_c}\otimes M_{m_c}$,
$$
\sum_ck_c^2=\dim\operatorname{span}\pi(S_n),\qquad\sum_cm_c^2=\binom{d^2+n-1}n,
$$
и $\sum_ck_c^2=n!$ при $n\le d$.

*Доказательство.* Предложение 2.7 для конечной группы $S_n$ и теорема 2.6 для $\mathcal M=\operatorname{span}\pi(S_n)=\pi(S_n)''$ дают
$W$, представления $\rho_c$ с указанными свойствами и $W^*\pi(S_n)'W=\bigoplus_c1_{k_c}\otimes M_{m_c}$. По теореме B.6(a)
$\pi(S_n)'=\mathcal U_n$. Поэтому $W^*V^{\otimes n}W=\bigoplus_c1\otimes\tau_c(V)$ для некоторых $\tau_c(V)\in M_{m_c}$, и $V\mapsto\tau_c(V)$ —
унитарное представление, поскольку таково $V\mapsto V^{\otimes n}$. Так как $V^{\otimes n}$ порождают $\mathcal U_n$, линейное отображение
$\sum_la_lV_l^{\otimes n}\mapsto(\sum_la_l\tau_c(V_l))_c$ сюръективно на $\bigoplus_cM_{m_c}$; неприводимость и попарная
неэквивалентность $\tau_c$ следуют дословно так же, как в доказательстве предложения 2.7. Размерности — это теорема 2.6,
предложение B.3 и лемма B.1(d). $\square$

При $d=2$ секторы отождествляются по полному спину:

**Следствие B.9 ($L$ спинов $\frac12$).** Пусть $V$ — унитарная матрица с $V^*S_aV=\bigoplus_SS^{(S)}_a\otimes1_{m_S}$ для полного
спина $L$ спинов $\frac12$ (следствие 4.12, предложение 4.13). Тогда
$V^*\pi(\sigma)V=\bigoplus_S1_{2S+1}\otimes\rho_S(\sigma)$, где $\rho_S$ — попарно неэквивалентные неприводимые унитарные
представления $S_L$ размерности
$$
m_S=\binom Lj-\binom L{j-1},\qquad j=\tfrac L2-S\in\{0,1,\dots,\lfloor L/2\rfloor\}.
$$
В частности, $\dim\operatorname{span}\pi(S_L)=\sum_Sm_S^2$ и $\dim\pi(S_L)'=\sum_S(2S+1)^2=\binom{L+3}3$.

*Доказательство.* По теореме 4.3 $V^*\{S_a\}'V=\bigoplus_S1_{2S+1}\otimes M_{m_S}$, а $\{S_a\}'=\operatorname{span}\pi(S_L)$ по следствию B.7. Поэтому
$V^*\pi(\sigma)V=\bigoplus_S1\otimes\rho_S(\sigma)$ с унитарным представлением $\rho_S$, и отображение
$\sum_\sigma a_\sigma\pi(\sigma)\mapsto(\sum_\sigma a_\sigma\rho_S(\sigma))_S$ сюръективно на $\bigoplus_SM_{m_S}$; неприводимость и
неэквивалентность следуют, как в доказательстве предложения 2.7. Значения $m_S$ — это предложение 4.13. По теореме 4.3
$\dim\{S_a\}'=\sum_Sm_S^2$ и $\dim\{S_a\}''=\sum_S(2S+1)^2$; по следствию B.7 $\{S_a\}'=\operatorname{span}\pi(S_L)$ и
$\{S_a\}''=\pi(S_L)'$, размерность которого равна $\binom{4+L-1}L=\binom{L+3}3$ по предложению B.3 при $d=2$. $\square$

**Следствие B.10 (числа Каталана).** При $L\ge1$ $\displaystyle\sum_Sm_S^2=\binom{2L}L\frac1{L+1}$; это
$\dim\{S_a\}'=\dim\operatorname{span}\pi(S_L)$ на $(\mathbb C^2)^{\otimes L}$.

*Доказательство.* Положим $N=2^L$ и пусть $\operatorname{vec}:M_N\to\mathbb C^N\otimes\mathbb C^N$ — линейная биекция
$|x\rangle\langle y|\mapsto x\otimes\bar y$. Тогда $\operatorname{vec}(AX-XA)=(A\otimes1-1\otimes A^T)\operatorname{vec}X$, поскольку
$A|x\rangle\langle y|B=|Ax\rangle\langle B^*y|\mapsto Ax\otimes B^T\bar y$. Поэтому $\{S_a\}'$ отображается на общее ядро
$S_a\otimes1-1\otimes S_a^T$, $a=x,y,z$, где $S_a^T=\sum_i(s_a^T)^{(i)}$. Из $s_x^T=s_x$, $s_y^T=-s_y$, $s_z^T=s_z$ и
$\sigma_ys_x\sigma_y=-s_x$, $\sigma_ys_y\sigma_y=s_y$, $\sigma_ys_z\sigma_y=-s_z$ получаем $\sigma_y(-s_a^T)\sigma_y=s_a$, так что унитарная
матрица $Y=1\otimes\sigma_y^{\otimes L}$ даёт $Y(S_a\otimes1-1\otimes S_a^T)Y^*=S_a\otimes1+1\otimes S_a$ — полный спин $2L$ спинов
$\frac12$ на $(\mathbb C^2)^{\otimes2L}=\mathbb C^N\otimes\mathbb C^N$. В разложении следствия 4.12 для $2L$ спинов вектор,
аннулируемый всеми тремя компонентами, аннулируется и $\mathbf S^2$, который равен $S(S+1)\ne0$ на слагаемых спина
$S\ge\frac12$; обратно, компоненты обращаются в нуль на слагаемом спина $0$. Поэтому общее ядро — это слагаемое
размерности $m_0(2L)=\binom{2L}L-\binom{2L}{L-1}=\binom{2L}L\big(1-\frac L{L+1}\big)$ (предложение 4.13). Первое равенство
размерностей — это теорема 4.3, второе — следствие B.7. $\square$

Код раскладывает $\pi(S_3)$ на $(\mathbb C^3)^{\otimes3}$ и проверяет размерности теоремы B.8, затем раскладывает
$\pi(S_L)$ на $(\mathbb C^2)^{\otimes L}$ при $L=5,6$ и сравнивает пары $(k_c,m_c)$ с $(m_S,2S+1)$ из следствия B.9.

```python
D33 = dc.from_group(dc.permutation_representation(3, 3))          # S_3 on (C^3)^{(x)3}
print(D33)
assert D33.verify()['ok']
assert sum(k * k for k in D33.k_factors) == 6 == math.factorial(3)            # Theorem B.8, n <= d
assert sum(m * m for m in D33.multiplicities) == math.comb(11, 3)              # 165

for L in (5, 6):
    DL = dc.from_group(dc.permutation_representation(L, 2))
    pairs = sorted((math.comb(L, j) - (math.comb(L, j - 1) if j else 0), L - 2 * j + 1)
                   for j in range(L // 2 + 1))                                 # (m_S, 2S + 1), S = L/2 - j
    print(L, DL.k_factors, DL.multiplicities)
    assert sorted(zip(DL.k_factors, DL.multiplicities)) == pairs               # Corollary B.9
    assert sum(k * k for k in DL.k_factors) == math.comb(2 * L, L) // (L + 1)  # Corollary B.10
    assert sum(m * m for m in DL.multiplicities) == math.comb(L + 3, 3)
```

## B.5 Исчисление Вайнгартена

**Лемма B.11 (твирлинг по Хаару — это ортогональная проекция на коммутант).** Пусть $U$ распределена по Хаару на
$U(d)$ и $\Phi(X)=\mathbb E[U^{\otimes n}X(U^{\otimes n})^*]$ для $X\in\operatorname{End}\mathcal H_n$. Тогда $\Phi$ — ортогональная
проекция $\operatorname{End}\mathcal H_n$ на $\{V^{\otimes n}:V\in U(d)\}'=\operatorname{span}\pi(S_n)$.

*Доказательство.* Элементы $U^{\otimes n}X(U^{\otimes n})^*$ — многочлены от элементов $U$ и $\bar U$, которые по модулю не
превосходят $1$, поэтому $\Phi$ — корректно определённое линейное отображение. Для $V\in U(d)$ матрица $VU$ имеет то же
распределение, что и $U$, поэтому
$V^{\otimes n}\Phi(X)(V^{\otimes n})^*=\mathbb E[(VU)^{\otimes n}X((VU)^{\otimes n})^*]=\Phi(X)$: образ лежит в коммутанте. Если $Y$ лежит в
коммутанте, то $U^{\otimes n}Y(U^{\otimes n})^*=Y$ и $\Phi(Y)=Y$. Итак, $\Phi$ идемпотентно, и его образ — коммутант. По
лемме 5.1 $U^*=U^{-1}$ имеет то же распределение, что и $U$, следовательно,
$\langle Y,\Phi(X)\rangle=\mathbb E\operatorname{Tr}\big(((U^*)^{\otimes n}YU^{\otimes n})^*X\big)=\langle\Phi(Y),X\rangle$: $\Phi$ самосопряжено.
Самосопряжённый идемпотент — это ортогональная проекция на свой образ
($\langle X-\Phi(X),\Phi(Z)\rangle=\langle\Phi(X-\Phi(X)),Z\rangle=0$). Коммутант есть $\operatorname{span}\pi(S_n)$ по теореме B.6(c). $\square$

Чтобы записать эту проекцию через порождающее семейство $\pi(\sigma)$, которое линейно зависимо при $n>d$, используем
**псевдообратную** к эрмитовой матрице $G$: $G^+=f(G)$ с $f(\lambda)=1/\lambda$ при $\lambda\neq0$ и $f(0)=0$
(функциональное исчисление, теорема 1.5). Она эрмитова, $G^+=G^{-1}$, если $G$ обратима, $GG^+=\mathbf 1_{\neq0}(G)$ —
проекция на $(\ker G)^\perp$, $G^+v=v/\lambda$, если $Gv=\lambda v$ с $\lambda\neq0$, и
$(R^*GR)^+=R^*G^+R$ для унитарной $R$ (если $G=Q\Lambda Q^*$, то $R^*GR=(R^*Q)\Lambda(R^*Q)^*$).

**Лемма B.12.** Пусть $v_1,\dots,v_N$ — векторы конечномерного пространства со скалярным произведением $\mathcal K$,
$G_{ab}=\langle v_a,v_b\rangle$ — их матрица Грама и $W=G^+$. Тогда ортогональная проекция на $\operatorname{span}\{v_a\}$ есть
$x\mapsto\sum_{a,b}v_a\,W_{ab}\,\langle v_b,x\rangle$.

*Доказательство.* Пусть $T:\mathbb C^N\to\mathcal K$, $Tc=\sum_ac_av_a$. Тогда $T^*x=(\langle v_b,x\rangle)_b$, $G=T^*T$, а
рассматриваемое отображение есть $Q=TG^+T^*$. Так как $\langle c,Gc\rangle=\|Tc\|^2$, имеем $\ker G=\ker T$, поэтому
$TG^+G=T(1-\Pi)=T$ для проекции $\Pi$ на $\ker T$. Следовательно, $QT=T$: $Q$ тождественно на $\operatorname{ran}T$. $Q$
обращается в нуль на $(\operatorname{ran}T)^\perp=\ker T^*$, и $Q^*=Q$. Поэтому $Q$ — ортогональная проекция на
$\operatorname{ran}T$. $\square$

**Теорема B.13 (формула Вайнгартена любого порядка).** Пусть $n,d\ge1$, $U$ распределена по Хаару на $U(d)$, $G$ —
матрица размера $n!\times n!$ с $G_{\sigma\tau}=d^{\#\text{cycles}(\sigma^{-1}\tau)}$ ($\sigma,\tau\in S_n$) и $W=G^+$.
(a) $W_{\sigma\tau}=\mathrm{Wg}(\sigma^{-1}\tau)$ для функции $\mathrm{Wg}(\rho)=\mathrm{Wg}(\rho,d)=W_{e\rho}$ на $S_n$, которая зависит
только от цикловой структуры $\rho$; в частности, $\mathrm{Wg}(\rho^{-1})=\mathrm{Wg}(\rho)$.
(b) $G$ обратима тогда и только тогда, когда $n\le d$, и тогда $W=G^{-1}$.
(c) Для $X\in\operatorname{End}\mathcal H_n$
$\mathbb E[U^{\otimes n}X(U^{\otimes n})^*]=\sum_{\sigma,\tau}\mathrm{Wg}(\sigma^{-1}\tau)\operatorname{Tr}(\pi(\tau)^*X)\,\pi(\sigma)$.
(d) Для всех $\mathbf i,\mathbf j,\mathbf i',\mathbf j'\in[d]^n$
$$
\mathbb E\big[U_{i_1j_1}\cdots U_{i_nj_n}\overline{U_{i'_1j'_1}\cdots U_{i'_nj'_n}}\big]
=\sum_{\sigma,\tau\in S_n}\prod_{m=1}^n\delta_{i_mi'_{\sigma(m)}}\delta_{j_mj'_{\tau(m)}}\;\mathrm{Wg}(\sigma\tau^{-1}),
$$
т. е. сумма берётся по парам с $\mathbf i=\mathbf i'\circ\sigma$ и $\mathbf j=\mathbf j'\circ\tau$.
(e) Моменты с разным числом $p\neq q$ множителей $U$ и $\bar U$ равны нулю.

Это теорема 5.12 (с $k=d$), теперь доказанная.

*Доказательство.* (a) Запишем $g(\rho)=d^{\#\text{cycles}(\rho)}$, так что $G_{\sigma\tau}=g(\sigma^{-1}\tau)$, и пусть $\Lambda_\kappa,R_\kappa$ —
матрицы перестановок на $\mathbb C^{S_n}$ с $\Lambda_\kappa e_\sigma=e_{\kappa\sigma}$, $R_\kappa e_\sigma=e_{\sigma\kappa}$. Тогда
$(\Lambda_\kappa^*G\Lambda_\kappa)_{\sigma\tau}=g(\sigma^{-1}\kappa^{-1}\kappa\tau)=G_{\sigma\tau}$ и
$(R_\kappa^*GR_\kappa)_{\sigma\tau}=g(\kappa^{-1}\sigma^{-1}\tau\kappa)=G_{\sigma\tau}$ (лемма B.1(c)). По свойствам
псевдообратной $W$ обладает теми же двумя инвариантностями: $W_{\kappa\sigma,\kappa\tau}=W_{\sigma\tau}$ и
$W_{\sigma\kappa,\tau\kappa}=W_{\sigma\tau}$. Первая при $\kappa=\sigma^{-1}$ даёт $W_{\sigma\tau}=\mathrm{Wg}(\sigma^{-1}\tau)$; вторая
даёт $\mathrm{Wg}(\kappa^{-1}\rho\kappa)=\mathrm{Wg}(\rho)$, так что $\mathrm{Wg}$ постоянна на классах сопряжённости, а это цикловые
структуры (лемма B.1(c)); у $\rho^{-1}$ та же цикловая структура, что и у $\rho$.
(b) $G$ — матрица Грама операторов $\pi(\sigma)$ (лемма B.1(c)), и с $T$ из доказательства леммы B.12 матрица $G=T^*T$
обратима тогда и только тогда, когда $T$ инъективно, т. е. когда $\pi(\sigma)$ линейно независимы: по лемме B.1(d) это
равносильно $n\le d$.
(c) Леммы B.11 и B.12 с $\mathcal K=\operatorname{End}\mathcal H_n$ и $v_\sigma=\pi(\sigma)$.
(d) Элементы кронекерова произведения равны $\langle e_{\mathbf i},U^{\otimes n}e_{\mathbf j}\rangle=\prod_mU_{i_mj_m}$, поэтому момент
равен $\langle e_{\mathbf i},\Phi(X)e_{\mathbf i'}\rangle$ для $X=|e_{\mathbf j}\rangle\langle e_{\mathbf j'}|$. Здесь
$\operatorname{Tr}(\pi(\tau)^*X)=\langle\pi(\tau)e_{\mathbf j'},e_{\mathbf j}\rangle=\delta_{\mathbf j,\mathbf j'\circ\tau^{-1}}$ и
$\langle e_{\mathbf i},\pi(\sigma)e_{\mathbf i'}\rangle=\delta_{\mathbf i,\mathbf i'\circ\sigma^{-1}}$. По (c) момент равен
$\sum_{\sigma,\tau}\delta_{\mathbf i,\mathbf i'\circ\sigma^{-1}}\delta_{\mathbf j,\mathbf j'\circ\tau^{-1}}\mathrm{Wg}(\sigma^{-1}\tau)$; подстановка
$\sigma\to\sigma^{-1}$, $\tau\to\tau^{-1}$ даёт формулу.
(e) Для вещественного $\theta$ имеем $e^{i\theta}1\in U(d)$, поэтому $e^{i\theta}U$ имеет то же распределение, что и $U$, и
момент равен самому себе, умноженному на $e^{i(p-q)\theta}$. Возьмём $\theta$ с $e^{i(p-q)\theta}\neq1$. $\square$

При $n=2\le d$ имеем $G=\begin{pmatrix}d^2&d\\d&d^2\end{pmatrix}$, и $G^{-1}$ даёт $\mathrm{Wg}(e)=\frac1{d^2-1}$,
$\mathrm{Wg}(\tau)=-\frac1{d(d^2-1)}$: это теорема 5.10. При $d=1$ все $G_{\sigma\tau}=1$, $W=\frac1{(n!)^2}G$, и момент
фазы $U=U_{11}$ равен $\sum_{\sigma,\tau}\frac1{(n!)^2}=1=|U_{11}|^{2n}$, как и должно быть.

**Следствие B.14.** Для всех $n,d\ge1$ $\mathbb E|U_{11}|^{2n}=\dfrac{n!}{d(d+1)\cdots(d+n-1)}$.

*Доказательство.* Если все индексы в теореме B.13(d) равны $1$, момент равен $\sum_{\sigma,\tau}\mathrm{Wg}(\sigma\tau^{-1})=\sum_{\sigma,\tau}W_{\sigma^{-1},\tau^{-1}}=\mathbf 1^TW\mathbf 1$ (часть (a)). Каждая
строчная сумма $G$ равна $c_n(d)=\sum_{\rho\in S_n}d^{\#\text{cycles}(\rho)}\ge d^n>0$, поэтому $G\mathbf 1=c_n(d)\mathbf 1$,
$W\mathbf 1=\mathbf 1/c_n(d)$, и момент равен $n!/c_n(d)$. Наконец, $c_n(x)=x(x+1)\cdots(x+n-1)$ (глава 5, упражнение 2,
решение). $\square$

Это распространяет упражнение 2 главы 5 на случай $n>d$, когда $G$ вырождена. Код вычисляет $\mathrm{Wg}$ через
псевдообратную, проверяет (a), (b) и теорему 5.10 и сравнивает весь тензор моментов порядка $n=3$ при $d=2$ — случай
с вырожденной $G$ — и твирлинг случайной $X$ со средними по методу Монте-Карло по унитарным матрицам Хаара.

```python
def n_cycles(p):                                        # number of cycles of a permutation in one-line form
    seen, count = set(), 0
    for i in range(len(p)):
        count += i not in seen
        while i not in seen:
            seen.add(i)
            i = p[i]
    return count

def compose(s, t):                                     # (s t)(m) = s(t(m))
    return tuple(s[m] for m in t)

def inverse(s):
    return tuple(sorted(range(len(s)), key=lambda m: s[m]))

def weingarten(n, d):
    perms = list(itertools.permutations(range(n)))
    G = torch.tensor([[float(d) ** n_cycles(compose(inverse(s), t)) for t in perms] for s in perms],
                     dtype=torch.float64)
    return perms, G, torch.linalg.pinv(G, hermitian=True)

for n, d in ((2, 2), (3, 2), (3, 3), (4, 2), (4, 3)):
    perms, G, W = weingarten(n, d)
    pos = {p: a for a, p in enumerate(perms)}
    Wg = {p: W[0, pos[p]].item() for p in perms}        # perms[0] is the identity
    assert all(abs(W[a, b] - Wg[compose(inverse(s), t)]) < 1e-12
               for a, s in enumerate(perms) for b, t in enumerate(perms))       # Theorem B.13(a)
    assert all(abs(Wg[p] - Wg[compose(compose(inverse(r), p), r)]) < 1e-12 for p in perms for r in perms)
    assert (torch.linalg.matrix_rank(G).item() == math.factorial(n)) == (n <= d)  # Theorem B.13(b)
for d in (2, 3, 5):                                     # Theorem 5.10
    _, _, W = weingarten(2, d)
    assert close(W[0], torch.tensor([1 / (d * d - 1), -1 / (d * (d * d - 1))], dtype=torch.float64), 1e-12)

n, d, B = 3, 2, 100_000
perms, G, W = weingarten(n, d)
P = dc.permutation_representation(n, d, [list(p) for p in perms])        # pi(sigma), same order as W
Phi = lambda X: torch.einsum('st,tij,ij,skl->kl', W.to(dt), P.conj(), X, P)   # Theorem B.13(c)
U3 = power(haar(d, B), n)                                                 # B samples of U^{(x)3}
moments = torch.einsum('bx,by->xy', U3.reshape(B, -1), U3.reshape(B, -1).conj()) / B
second = torch.einsum('bx,by->xy', U3.reshape(B, -1).abs() ** 2, U3.reshape(B, -1).abs() ** 2) / B
err = (second - moments.abs() ** 2).clamp(min=0).sqrt() / math.sqrt(B)   # standard errors
N = d ** n
exact = torch.stack([Phi(torch.outer(torch.eye(N, dtype=dt)[j], torch.eye(N, dtype=dt)[jp]))
                     for j in range(N) for jp in range(N)])               # [(j, j'), i, i']
exact = exact.reshape(N, N, N, N).permute(2, 0, 3, 1).reshape(N * N, N * N)   # [(i, j), (i', j')]
z = ((moments - exact).abs() / (err + 1e-15)).max().item()
print(f"{N ** 4} moments of order 3 for d = 2: largest deviation {z:.2f} standard errors")
assert ((moments - exact).abs() <= 5 * err + 1e-12).all()                 # Theorem B.13(d), singular G
assert abs(exact[0, 0].item() - 1 / 4) < 1e-12                            # Corollary B.14: 3!/(2*3*4)

X = torch.randn(N, N, dtype=dt)
T = Phi(X)
coef = torch.linalg.lstsq(P.reshape(len(perms), -1).T, T.reshape(-1, 1)).solution
assert close(Phi(T), T) and close(torch.einsum('s,sij->ij', coef[:, 0], P), T)   # a projection into span pi(S_n)
Vn = power(haar(d, 1), n)[0]
assert close(Vn @ T, T @ Vn)                                             # = {V^{(x)n}}'
assert close(Phi(P[3]), P[3])
Y = U3 @ X @ U3.mH                                                        # samples of U X U*
assert ((Y.mean(0) - T).abs() <= 5 * Y.std(0) / math.sqrt(B) + 1e-12).all()   # Lemma B.11
```

## B.6 Тензорные операторы и теорема Вигнера–Эккарта для $SU(2)$

В разделе 9.3 спиновые корреляции строятся с помощью теоремы Вигнера–Эккарта. Докажем её для спиновых
представлений из раздела 4.5. Спин $k$ — это представление $S^{(k)}_a$ из теоремы 4.11 на $\mathbb C^{2k+1}$ с базисом
$|k,q\rangle$, $q=k,\dots,-k$; напомним, что $S^{(k)}_\pm|k,q\rangle=\sqrt{k(k+1)-q(q\pm1)}\,|k,q\pm1\rangle$ (для $S_-=S_+^*$ это
следует из формулы для $S_+$).

**Определение.** Пусть $(S_a)$ — спиновое представление на $\mathbb C^N$ и $k\in\{0,\frac12,1,\dots\}$. **Тензорный оператор
ранга $k$** — это семейство $(T_q)_{q=-k,\dots,k}$ матриц $N\times N$ с
$$
[S_z,T_q]=q\,T_q,\qquad[S_\pm,T_q]=\sqrt{k(k+1)-q(q\pm1)}\;T_{q\pm1}
$$
(коэффициент обращается в нуль при $q=\pm k$, где $T_{\pm(k+1)}$ не определено). **Векторный оператор** — это тройка
$(T_x,T_y,T_z)$ с $[S_a,T_b]=i\sum_c\varepsilon_{abc}T_c$, т. е. $[S_x,T_y]=iT_z$, $[S_x,T_z]=-iT_y$ и циклически, и
$[S_a,T_a]=0$.

**Лемма B.15.** (a) Если $(T_a)$ — векторный оператор, то его **сферические компоненты** $T_{\pm1}=\mp\frac1{\sqrt2}(T_x\pm iT_y)$,
$T_0=T_z$ образуют тензорный оператор ранга $1$, и для двух векторных операторов
$\sum_aT_aT'_a=\sum_q(-1)^qT_qT'_{-q}$.
(b) Для $L$ спинов $\frac12$ и полного спина $S_a$ каждый спин узла $(s^{(i)}_a)_a$ и сам $(S_a)_a$ являются векторными операторами.
(c) $(T_q)$ — тензорный оператор ранга $k$ тогда и только тогда, когда отображение $\Psi:\mathbb C^N\otimes\mathbb C^{2k+1}\to\mathbb C^N$,
$\Psi(v\otimes|k,q\rangle)=T_qv$, удовлетворяет $\Psi\,(S_a\otimes1+1\otimes S^{(k)}_a)=S_a\Psi$ для $a=x,y,z$.

*Доказательство.* (a) Из соотношений $[S_z,T_x\pm iT_y]=\pm(T_x\pm iT_y)$, $[S_\pm,T_z]=\mp(T_x\pm iT_y)$,
$[S_\pm,T_x\mp iT_y]=\pm2T_z$ и $[S_\pm,T_x\pm iT_y]=0$ (например, $[S_x+iS_y,T_x-iT_y]=-i[S_x,T_y]+i[S_y,T_x]=T_z+T_z$). В
сферических компонентах: $[S_z,T_q]=qT_q$, $[S_\pm,T_0]=\sqrt2\,T_{\pm1}$, $[S_\pm,T_{\mp1}]=\sqrt2\,T_0$, $[S_\pm,T_{\pm1}]=0$,
а это определяющие соотношения для $k=1$ ($\sqrt{2-q(q\pm1)}$ равно $\sqrt2$ при $q=0$ и $q=\mp1$ и $0$ при $q=\pm1$).
Для произведения $T_{+1}T'_{-1}+T_{-1}T'_{+1}=-\frac12\big((T_x+iT_y)(T'_x-iT'_y)+(T_x-iT_y)(T'_x+iT'_y)\big)=-(T_xT'_x+T_yT'_y)$.
(b) $[S_a,s^{(i)}_b]=[s^{(i)}_a,s^{(i)}_b]$ (остальные узлы коммутируют), а матрицы спина $\frac12$ удовлетворяют спиновым
соотношениям; для $S_a$ просуммируем по $i$. (c) $\Psi(S_a\otimes1+1\otimes S^{(k)}_a)(v\otimes|k,q\rangle)=T_qS_av+\sum_{q'}(S^{(k)}_a)_{q'q}T_{q'}v$,
поэтому сплетающее соотношение для $S_a$ означает $[S_a,T_q]=\sum_{q'}(S^{(k)}_a)_{q'q}T_{q'}$ для всех $q$. Для $S_z,S_\pm$ это
определяющие соотношения, а соотношения для $S_x,S_y,S_z$ и для $S_z,S_\pm$ являются линейными комбинациями друг друга. $\square$

**Теорема B.16 (Вигнер–Эккарт).** Пусть $(S_a)$ — спиновое представление на $\mathbb C^N$, $V$ — унитарная матрица с
$V^*S_aV=\bigoplus_SS^{(S)}_a\otimes1_{m_S}$ (следствие 4.12), $|S,M;b\rangle=V(|S,M\rangle\otimes f_b)$ — соответствующий
ортонормированный базис ($f_b$ — стандартный базис $\mathbb C^{m_S}$) и $(T_q)$ — тензорный оператор ранга $k$. Пусть $S,S'$
встречаются ($m_S,m_{S'}\ge1$).
(a) Если $S'\notin\{|S-k|,|S-k|+1,\dots,S+k\}$, то $\langle S',M';a|T_q|S,M;b\rangle=0$ для всех $M',M,q,a,b$.
(b) В противном случае сплетающие операторы из спина $S\otimes$ спин $k$ в спин $S'$ образуют одномерное пространство; для
каждого ненулевого такого $C=C^{S'}_{S,k}:\mathbb C^{2S+1}\otimes\mathbb C^{2k+1}\to\mathbb C^{2S'+1}$ существует единственная
матрица $R=\langle S'\|T\|S\rangle$ размера $m_{S'}\times m_S$, **приведённый матричный элемент**, такая, что
$$
\langle S',M';a|\,T_q\,|S,M;b\rangle=\langle S',M'|\,C\,(|S,M\rangle\otimes|k,q\rangle)\rangle\;R_{ab}\qquad\text{for all }M',M,q,a,b .
$$
(c) Элемент равен нулю, если не выполнено $M'=M+q$. В случае (b) существует $q$ с
$\langle S',S'|C(|S,S'-q\rangle\otimes|k,q\rangle)\rangle\neq0$, поэтому $R$, а вместе с ней и весь блок $T$ между мультиплетами
$S$ и $S'$, определяется матричными элементами $\langle S',S';a|T_q|S,S'-q;b\rangle$ между векторами старшего веса $S'$ и
одним вектором каждой копии $S$.

*Доказательство.* Зафиксируем $a\le m_{S'}$, $b\le m_S$ и пусть $J:\mathbb C^{2S+1}\to\mathbb C^N$, $Jx=V(x\otimes f_b)$ (в слагаемом $S$),
а $K:\mathbb C^N\to\mathbb C^{2S'+1}$ — сопряжённое к $x\mapsto V(x\otimes f_a)$ (в слагаемом $S'$). Обозначим через
$\alpha\in\{x,y,z\}$ индекс компоненты. Из вида $V^*S_\alpha V$ следует $S_\alpha J=JS^{(S)}_\alpha$, а переходя к сопряжённым в том же
соотношении для $S'$ (спиновые матрицы эрмитовы), получаем $KS_\alpha=S^{(S')}_\alpha K$. Пусть $t=K\Psi(J\otimes1)$ с $\Psi$ из
леммы B.15(c); её матричные элементы равны $\langle S',M'|t(|S,M\rangle\otimes|k,q\rangle)\rangle=\langle S',M';a|T_q|S,M;b\rangle$, и
$$
t\,(S^{(S)}_\alpha\otimes1+1\otimes S^{(k)}_\alpha)=K\Psi(S_\alpha\otimes1+1\otimes S^{(k)}_\alpha)(J\otimes1)=KS_\alpha\Psi(J\otimes1)=S^{(S')}_\alpha\,t,
$$
так что $t$ — сплетающий оператор из $\rho=$ спин $S\otimes$ спин $k$ (спиновое представление
$S^{(S)}_\alpha\otimes1+1\otimes S^{(k)}_\alpha$, раздел 4.5) в $\sigma=$ спин $S'$. Переход к сопряжённым — линейная биекция между
$\operatorname{Hom}(\rho,\sigma)$ и $\operatorname{Hom}(\sigma,\rho)$ (сопряжённое к $t\rho(\alpha)=\sigma(\alpha)t$ есть
$\rho(\alpha)t^*=t^*\sigma(\alpha)$, так как все участвующие матрицы эрмитовы). По следствию 4.12 и упражнению 2 главы 4
(доказанному в его решении только с помощью следствия 4.12) $\rho$ унитарно эквивалентно $\bigoplus_j(\text{spin }j)$, по одной
копии каждого $j\in\{|S-k|,|S-k|+1,\dots,S+k\}$. Эти спины неприводимы и попарно неэквивалентны (теорема 4.11), так что это
разложение в смысле теоремы 4.3, и её последнее утверждение даёт $\dim\operatorname{Hom}(\sigma,\rho)=1$, если $S'$ — одно из
этих $j$, и $\operatorname{Hom}(\sigma,\rho)=0$ иначе (тогда спин $S'$ не эквивалентен ни одному из них по теореме 4.11).
Это даёт (a), а в случае (b) $t=R_{ab}C$ для единственного числа $R_{ab}$.
(c) Сплетающее соотношение для $S_z$ даёт $M'\langle S',M'|t(|S,M\rangle\otimes|k,q\rangle)\rangle=(M+q)\langle S',M'|t(|S,M\rangle\otimes|k,q\rangle)\rangle$.
Образ $C\neq0$ — ненулевое инвариантное подпространство неприводимого спина $S'$, а значит, всё пространство, поэтому
$\langle S',S'|C(|S,M\rangle\otimes|k,q\rangle)\rangle\neq0$ для некоторого базисного вектора, и тогда $M=S'-q$ по тому же соотношению
для сплетающего оператора $C$. Разделим на это число. $\square$

По лемме B.15(a) блок между мультиплетами $S$ и $S'$ скалярного произведения вроде
$\mathbf s^{(i)}\cdot\mathbf s^{(j)}=\sum_q(-1)^q s^{(i)}_qs^{(j)}_{-q}$ равен
$\sum_{S''}\sum_q(-1)^q(\text{блок }S''\to S'\text{ из }s^{(i)}_q)(\text{блок }S\to S''\text{ из }s^{(j)}_{-q})$, и каждый множитель задаётся
теоремой B.16 через приведённый матричный элемент и отображение Клебша–Гордана. Это принцип, на котором основана
`su2.spin_correlation`; корректность её реализации глава 9 проверяет численно, сравнивая с методом
`'direct'`. Код проверяет лемму B.15 для спинов узлов шести спинов $\frac12$ и теорему B.16 в бескоординатной
форме: в разложении, вычисленном библиотекой (базисы секторов которого не обязаны быть $|S,M\rangle$), блок
$(T_q)_q$ между двумя мультиплетами, как матрица со строками $(M',M,q)$ и столбцами $(a,b)$, имеет ранг не больше единицы
и равен нулю вне диапазона из (a).

```python
L = 6
S6 = dc.total_spin(L)
Sp, Sm = S6[0] + 1j * S6[1], S6[0] - 1j * S6[1]
D6 = dc.from_generators(S6)
print(D6)
off = [0]
for k, m in zip(D6.k_factors, D6.multiplicities):
    off.append(off[-1] + k * m)
sx, sy, sz = dc.spin_matrices(0.5)
site = lambda x, i: functools.reduce(torch.kron, [x if l == i else torch.eye(2, dtype=dt) for l in range(L)])
spherical = lambda X, Y, Z: {1: -(X + 1j * Y) / math.sqrt(2), 0: Z, -1: (X - 1j * Y) / math.sqrt(2)}
for i in (0, 3):
    T = spherical(site(sx, i), site(sy, i), site(sz, i))
    zero = torch.zeros_like(T[0])
    for q in (-1, 0, 1):                                                   # Lemma B.15(a), (b)
        assert close(S6[2] @ T[q] - T[q] @ S6[2], q * T[q])
        assert close(Sp @ T[q] - T[q] @ Sp, math.sqrt(2 - q * (q + 1)) * T.get(q + 1, zero))
        assert close(Sm @ T[q] - T[q] @ Sm, math.sqrt(2 - q * (q - 1)) * T.get(q - 1, zero))
    for c, (kc, mc) in enumerate(zip(D6.k_factors, D6.multiplicities)):
        for e, (ke, me) in enumerate(zip(D6.k_factors, D6.multiplicities)):
            Vc, Ve = D6.V[:, off[c]:off[c + 1]], D6.V[:, off[e]:off[e + 1]]
            block = torch.stack([(Ve.mH @ T[q] @ Vc).reshape(ke, me, kc, mc) for q in (-1, 0, 1)])
            sv = torch.linalg.svdvals(block.permute(0, 1, 3, 2, 4).reshape(3 * ke * kc, me * mc))
            S, S2 = (kc - 1) / 2, (ke - 1) / 2
            if abs(S - S2) <= 1 and S + S2 >= 1:
                assert sv.numel() == 1 or sv[1] < 1e-10                        # Theorem B.16(b): rank one
            else:
                assert sv[0] < 1e-10                                           # Theorem B.16(a)
X, Y, Z = (site(s, 0) for s in (sx, sy, sz))
X2, Y2, Z2 = (site(s, 1) for s in (sx, sy, sz))
T, T2 = spherical(X, Y, Z), spherical(X2, Y2, Z2)
assert close(X @ X2 + Y @ Y2 + Z @ Z2, sum((-1) ** q * T[q] @ T2[-q] for q in (-1, 0, 1)))   # Lemma B.15(a)
```

## B.7 Меры Хаара, используемые в курсе

Курс использует меры Хаара следующих групп: $U(k)$ (теорема 5.3), $U(1)$ (случай $k=1$; глава 5,
упражнение 3), $U(\mathcal M)\cong\prod_cU(k_c)$ (следствие 5.4), $SU(k)$ (раздел 5.1, после следствия 5.4), ортогональной
группы $O(k)$ (раздел 5.1, глава 11, упражнение 2, и вещественные алгебры) и диагональных групп $T$ из леммы C.2.
Единственность для всех них — это лемма 5.1, поскольку каждая из них — замкнутая подгруппа некоторой $U(n)$;
существование в каждом случае устанавливается явной конструкцией (для $T$ — в лемме C.2), и конструкция для $O(k)$
приведена ниже. Никакая общая теорема существования меры Хаара компактной группы нигде в курсе не используется.

**Теорема B.17 (мера Хаара ортогональной группы).** Пусть $Z$ — вещественная матрица $k\times k$ с независимыми
элементами $N(0,1)$. Тогда $Z$ обратима с вероятностью единица, $Z=OR$ с ортогональной $O$ и вещественной верхнетреугольной
$R$ с положительной диагональю, и $O$ распределена по Хаару на $O(k)$. Если $Z=O'R'$ — любое вещественное QR-разложение
($O'$ ортогональна, $R'$ вещественная верхнетреугольная), то $O'\operatorname{diag}(\operatorname{sgn}R'_{ii})=O$. Для
распределённой по Хаару $O$ и $X\in M_k(\mathbb C)$ $\mathbb E[O]=0$ и $\mathbb E[OXO^T]=\frac{\operatorname{Tr}X}k1$; в частности,
$\mathbb E[O_{ij}O_{lm}]=\delta_{il}\delta_{jm}/k$.

*Доказательство.* $\det Z$ — многочлен от $k^2$ вещественных элементов, не равный тождественно нулю ($\det1=1$), а $Z$ имеет
плотность, поэтому $\det Z\neq0$ почти наверное (следствие A.7(b), как в доказательстве теоремы 5.3). Процесс Грама–Шмидта,
применённый к вещественным столбцам, даёт вещественные $q_j$ и вещественные коэффициенты, т. е. вещественное разложение
$Z=OR$ с $R_{jj}>0$; в силу единственности в лемме 5.2 это разложение леммы 5.2, так что $O$ непрерывно зависит от $Z$.
Плотность $(2\pi)^{-k^2/2}e^{-\|Z\|_2^2/2}$ матрицы $Z$ инвариантна относительно $Z\mapsto VZ$ для $V\in O(k)$ — линейного
отображения $\mathbb R^{k^2}$, сохраняющего $\|Z\|_2^2=\operatorname{Tr}(Z^TZ)$, а значит, ортогонального, с якобианом $1$;
поэтому $VZ$ имеет то же распределение, что и $Z$ (предложение A.14(a)). Так как $VZ=(VO)R$, единственность даёт
$O(VZ)=VO(Z)$, так что распределение $O$ левоинвариантно на компактной группе $O(k)$, и применима лемма 5.1. Поправка
знаков — это конец доказательства теоремы 5.3 с $\operatorname{sgn}R'_{ii}=R'_{ii}/|R'_{ii}|$ ($R'_{ii}\neq0$, так как $Z$
обратима). Поскольку $-1\in O(k)$, $-O$ имеет то же распределение, что и $O$, так что $\mathbb E[O]=-\mathbb E[O]=0$ (элементы
по модулю не превосходят $1$). Наконец, $T=\mathbb E[OXO^T]$ удовлетворяет $VTV^T=T$ для $V\in O(k)$, поскольку $VO$ имеет то же
распределение, что и $O$. Взяв в качестве $V$ диагональную матрицу с $-1$ на позиции $i$ и $1$ на остальных, получаем
$T_{ij}=-T_{ij}$ при $j\neq i$; с матрицами перестановок — что все $T_{ii}$ равны (как в лемме 5.6). Поэтому $T=\lambda1$ и
$k\lambda=\operatorname{Tr}T=\operatorname{Tr}X$. Для $X=|e_j\rangle\langle e_m|$ имеем $(OXO^T)_{il}=O_{ij}O_{lm}$. $\square$

Код проверяет методом Монте-Карло вторые моменты, $\mathbb E\operatorname{Tr}O=0$ и $\mathbb E(\operatorname{Tr}O)^2=\sum_{i,l}\mathbb E[O_{ii}O_{ll}]=1$.
Он также фиксирует *численное наблюдение* о процедуре QR в `torch.linalg.qr` (утверждение об этой процедуре, а не
теорема): без поправки знаков у её множителя $Q'$ величина $\mathbb E\operatorname{Tr}Q'$ отстоит от $0$ на много
стандартных ошибок, так что он не распределён по Хаару.

```python
k, B = 3, 100_000
Zr = torch.randn(B, k, k, dtype=torch.float64)
Q, R = torch.linalg.qr(Zr)                              # LAPACK: diag(R) of either sign
O = Q * torch.sign(R.diagonal(dim1=-2, dim2=-1)).unsqueeze(-2)     # Theorem B.17
assert close(O.mT @ O, torch.eye(k, dtype=torch.float64).expand(B, k, k))
assert (R.diagonal(dim1=-2, dim2=-1) * torch.sign(R.diagonal(dim1=-2, dim2=-1)) > 0).all()
second = torch.einsum('bij,blm->ijlm', O, O) / B                   # E O_ij O_lm
expected = torch.einsum('il,jm->ijlm', torch.eye(k, dtype=torch.float64), torch.eye(k, dtype=torch.float64)) / k
err = (torch.einsum('bij,blm->ijlm', O ** 2, O ** 2) / B - second ** 2).clamp(min=0).sqrt() / math.sqrt(B)
assert ((second - expected).abs() <= 5 * err).all()
tr, tr_raw = O.diagonal(dim1=-2, dim2=-1).sum(-1), Q.diagonal(dim1=-2, dim2=-1).sum(-1)
print(f"E Tr O = {tr.mean():.4f}, E (Tr O)^2 = {(tr ** 2).mean():.4f};  without the sign correction "
      f"E Tr Q = {tr_raw.mean():.4f}")
assert abs(tr.mean()) < 5 * tr.std() / math.sqrt(B)                # E Tr O = 0
assert abs((tr ** 2).mean() - 1) < 5 * (tr ** 2).std() / math.sqrt(B)   # E (Tr O)^2 = sum_i 1/k = 1
assert abs(tr_raw.mean()) > 20 * tr_raw.std() / math.sqrt(B)       # LAPACK's signs: not Haar
```

**Назад к:** [главе 4](04_symmetries.md) (теорема 4.14), [главе 5](05_random_matrices.md) (теорема 5.12),
[главе 9](09_spin_chains.md) (корреляции).
