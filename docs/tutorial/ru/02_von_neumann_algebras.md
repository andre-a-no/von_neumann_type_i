[English](../02_von_neumann_algebras.md) · [Español](../es/02_von_neumann_algebras.md) · [Français](../fr/02_von_neumann_algebras.md) · [Deutsch](../de/02_von_neumann_algebras.md) · [中文](../zh/02_von_neumann_algebras.md) · [日本語](../ja/02_von_neumann_algebras.md) · **Русский**

# 2. Конечномерные алгебры фон Неймана

Алгебра фон Неймана — это алгебра операторов, замкнутая относительно взятия «двойного
коммутанта». В конечномерном случае это чисто алгебраическое условие, и основная теорема
этой главы утверждает, что всякая такая алгебра с точностью до унитарной замены базиса
есть прямая сумма полных матричных алгебр, каждая из которых повторена с некоторой
кратностью. Именно на этой структуре построена библиотека.

**Вы узнаете** теорему о бикоммутанте, структурную теорему с полным доказательством, то,
как группа симметрий порождает такую алгебру, что такое центр и как всё это представлено в
библиотеке.

Всюду далее *$*$-алгебра на $\mathbb C^n$* означает комплексное линейное подпространство
$\mathcal M\subseteq M_n(\mathbb C)$, замкнутое относительно умножения и сопряжения и
содержащее единицу $1$.

## 2.1 Коммутанты и теорема о бикоммутанте

Для множества $S\subseteq M_n(\mathbb C)$ его **коммутант** — это
$S'=\{X\in M_n(\mathbb C):XA=AX\text{ для всех }A\in S\}$ (все $X$, коммутирующие со всеми $A\in S$). Он всегда является алгеброй с единицей;
если $S$ замкнуто относительно сопряжения ($S^*=S$), то $S'$ — $*$-алгебра, поскольку из
$XA=AX$ для всех $A\in S$ следует $X^*A^*=A^*X^*$ для всех $A\in S$. Каждый элемент $S$
коммутирует с каждым элементом $S'$, поэтому $S\subseteq S''$; применяя это к $S'$, получаем
$S'\subseteq S'''$, а переходя к коммутантам во включении $S\subseteq S''$ (что обращает включения),
получаем $S'''\subseteq S'$. Следовательно, $S'=S'''$.

Когда алгебра действует на подпространстве $E$ пространства $\mathbb C^n$, мы отождествляем $E$ с
$\mathbb C^{\dim E}$ посредством ортонормированного базиса, так что к нему применимы все
результаты, сформулированные для $\mathbb C^n$.

**Определение.** **Алгеброй фон Неймана** на $\mathbb C^n$ называется $*$-алгебра
$\mathcal M$, для которой $\mathcal M=\mathcal M''$. Её **центр** — это
$Z(\mathcal M)=\mathcal M\cap\mathcal M'$. **Фактор** — это алгебра фон Неймана с
$Z(\mathcal M)=\mathbb C1$.

В бесконечномерном случае $\mathcal M=\mathcal M''$ — топологическое условие (замкнутость в
слабой операторной топологии). В конечномерном случае оно выполняется автоматически:

**Лемма 2.1.** Пусть $\mathcal M$ — $*$-алгебра на $\mathbb C^n$, а $K\subseteq\mathbb C^n$ —
подпространство с $\mathcal MK\subseteq K$. Тогда ортогональный проектор $P$ на $K$ лежит
в $\mathcal M'$.

*Доказательство.* Для $A\in\mathcal M$ имеем $AK\subseteq K$ и $A^*K\subseteq K$ (так как $A^*\in\mathcal M$).
Второе включение даёт $AK^\perp\subseteq K^\perp$: для $y\in K^\perp$, $x\in K$
$\langle x,Ay\rangle=\langle A^*x,y\rangle=0$. Значит, $A$ переводит $K$ в $K$ и $K^\perp$
в $K^\perp$, что и означает $AP=PAP=PA$. $\square$

**Теорема 2.2 (теорема о бикоммутанте).** Всякая $*$-алгебра $\mathcal M$ на
$\mathbb C^n$ удовлетворяет равенству $\mathcal M=\mathcal M''$.

*Доказательство.* Имеем $\mathcal M\subseteq\mathcal M''$; пусть $T\in\mathcal M''$.

*Шаг 1: один вектор.* Для $x\in\mathbb C^n$ подпространство $K=\mathcal Mx$ инвариантно
относительно $\mathcal M$ и содержит $x=1x$. По лемме 2.1 его проектор $P$ лежит в $\mathcal M'$,
поэтому $TP=PT$ и $Tx=TPx=PTx\in K$: найдётся $A\in\mathcal M$ с $Tx=Ax$.

*Шаг 2: все векторы сразу.* На $\mathbb C^n\otimes\mathbb C^n=(\mathbb C^n)^{\oplus n}$
рассмотрим $*$-алгебру $\tilde{\mathcal M}=\{A\oplus\dots\oplus A:A\in\mathcal M\}$. Блочная
матрица $[X_{ij}]_{i,j=1}^n$ коммутирует со всеми $A\oplus\dots\oplus A$ тогда и только тогда, когда $X_{ij}A=AX_{ij}$
для всех $i,j$, поэтому $\tilde{\mathcal M}'=\{[X_{ij}]:X_{ij}\in\mathcal M'\}$. Оператор
$\tilde T=T\oplus\dots\oplus T$ коммутирует с каждой такой $[X_{ij}]$, поскольку
$TX_{ij}=X_{ij}T$. Следовательно, $\tilde T\in\tilde{\mathcal M}''$. Применим шаг 1 к $\tilde{\mathcal M}$,
$\tilde T$ и вектору $\xi=(e_1,\dots,e_n)$: найдётся $A\in\mathcal M$ с
$(Te_1,\dots,Te_n)=(Ae_1,\dots,Ae_n)$, т. е. $T=A\in\mathcal M$. $\square$

Итак, в конечномерном случае «алгебра фон Неймана» и «$*$-алгебра матриц с единицей» —
одно и то же. Тем не менее теорема о бикоммутанте задаёт правильную точку зрения:
физически алгебра наблюдаемых очень часто *задаётся* как коммутант — как множество
операторов, коммутирующих с симметрией (раздел 2.4).

## 2.2 Центр

**Лемма 2.3.** Пусть $\mathcal Z$ — коммутативная $*$-алгебра на $\mathbb C^n$. Существуют
ненулевые попарно ортогональные проекторы $P_1,\dots,P_C\in\mathcal Z$ с
$\sum_cP_c=1$, такие, что $\mathcal Z=\operatorname{span}\{P_1,\dots,P_C\}$. Это в точности
минимальные ненулевые проекторы $\mathcal Z$, и они единственны с точностью до порядка.

*Доказательство.* Всякий $z\in\mathcal Z$ нормален ($z^*\in\mathcal Z$ коммутирует с $z$). Пусть
$z_1,\dots,z_r$ порождают $\mathcal Z$ как линейное пространство. Разложим $\mathbb C^n$ в сумму
общих собственных подпространств. Собственные подпространства $z_1$ попарно ортогональны и в
сумме дают $\mathbb C^n$ (теорема 1.3); они инвариантны относительно $z_2,\dots,z_r$ и их
сопряжённых (которые коммутируют с $z_1$), поэтому сужение каждого
$z_j$ на собственное подпространство $z_1$ снова нормально, и каждое собственное
подпространство $z_1$ можно разложить на собственные подпространства $z_2$, и так далее.
В результате получаем ортогональное разложение $\mathbb C^n=E_1\oplus\dots\oplus E_C$ в сумму
ненулевых подпространств, на каждом из которых всякий
$z\in\mathcal Z$ действует как скаляр $\chi_c(z)$; при этом подпространства с одной и той же
функцией $\chi_c:\mathcal Z\to\mathbb C$ мы объединили, так что $\chi_c\neq\chi_d$ при $c\neq d$. Пусть
$P_c$ — проектор на $E_c$. Тогда $z=\sum_c\chi_c(z)P_c$ для любого $z\in\mathcal Z$.

Осталось показать, что $P_c\in\mathcal Z$. Для $d\neq c$ выберем $z_{cd}\in\mathcal Z$ с
$\chi_c(z_{cd})\neq\chi_d(z_{cd})$. Произведение
$\prod_{d\neq c}\big(z_{cd}-\chi_d(z_{cd})1\big)/\big(\chi_c(z_{cd})-\chi_d(z_{cd})\big)$
лежит в $\mathcal Z$, действует как $1$ на $E_c$ и как $0$ на каждом $E_d$, $d\neq c$; оно равно
$P_c$ (при $C=1$ имеем $P_1=1$: пустое произведение равно $1$). Минимальность и единственность: у проектора $\sum_c a_cP_c$ все $a_c\in\{0,1\}$, поэтому
ненулевые проекторы $\mathcal Z$ — это суммы $P_c$ по непустым множествам индексов. $\square$

Для $\mathcal Z=Z(\mathcal M)$ проекторы $P_c$ называются **минимальными центральными
проекторами** алгебры $\mathcal M$. Всякий $A\in\mathcal M$ раскладывается в сумму
$A=\sum_cAP_c$ с $AP_c=P_cAP_c$, поэтому
$\mathcal M=\bigoplus_c\mathcal MP_c$ (ортогональная прямая сумма алгебр, действующих на
$E_c=P_c\mathbb C^n$).

**Лемма 2.4.** Каждая $\mathcal M_c=\mathcal MP_c$, рассматриваемая как алгебра на $E_c$ с единицей $P_c$,
является фактором.

*Доказательство.* Пусть $y\in\mathcal M_c$ коммутирует с $\mathcal M_c$. Для $A\in\mathcal M$
$Ay=AP_cy=(AP_c)y=y(AP_c)=yP_cA=yA$, поскольку $P_c$ централен и $y=yP_c=P_cy$.
Значит, $y\in Z(\mathcal M)$, $y=\sum_d\chi_d(y)P_d$ и $y=yP_c=\chi_c(y)P_c$. $\square$

## 2.3 Структурная теорема

**Лемма 2.5 (факторы).** Пусть $\mathcal N$ — фактор на пространстве $E$. Существуют целые
$k,m\ge1$ и унитарный оператор $V:\mathbb C^k\otimes\mathbb C^m\to E$, для которых
$V^*\mathcal NV=M_k(\mathbb C)\otimes1_m$. Более того, $V^*\mathcal N'V=1_k\otimes M_m(\mathbb C)$.

*Доказательство.* *Шаг 1: минимальный проектор.* $\mathcal N$ содержит ненулевые проекторы (например, $1$).
Пусть $e\in\mathcal N$ — ненулевой проектор минимального ранга. Утверждается, что
$e\mathcal Ne=\mathbb Ce$. Пусть $x=x^*\in e\mathcal Ne$. Каждый спектральный проектор $q$ оператора $x$,
отвечающий ненулевому собственному значению $\lambda$, равен $q=p(x)$, где $p$ — многочлен
Лагранжа с $p(\lambda)=1$ и $p(\mu)=0$ для $\mu\in(\operatorname{spec}x\cup\{0\})\setminus\{\lambda\}$
(доказательство теоремы 1.5(a)). Поскольку $p(0)=0$ и $x=ex=xe$, имеем $q\in\mathcal N$ и $q=eqe$, так что
образ $q$ содержится в образе $e$. В силу минимальности $q=e$ или $q=0$; следовательно, у $x$
не более одного ненулевого собственного значения, и соответствующее собственное
подпространство совпадает с образом $e$, так что
$x=\lambda e$. Произвольный элемент $e\mathcal Ne$ имеет вид $x_1+ix_2$, где
$x_{1,2}=x_{1,2}^*\in e\mathcal Ne$.

*Шаг 2: $e$ порождает $E$.* Пусть $F=\operatorname{span}\mathcal Ne E$. Это подпространство инвариантно
относительно $\mathcal N$, поэтому его проектор $p_F$ лежит в $\mathcal N'$ (лемма 2.1). Оно инвариантно
и относительно $*$-алгебры $\mathcal N'$: для $b\in\mathcal N'$ имеем $bAe\xi=Aeb\xi\in F$; поэтому по лемме 2.1,
применённой к $\mathcal N'$, $p_F\in\mathcal N''=\mathcal N$ (теорема 2.2). Следовательно,
$p_F\in Z(\mathcal N)=\mathbb C1$, причём $p_F\neq0$ (так как $eE\subseteq F$), откуда $F=E$.

*Шаг 3: матричные единицы.* Назовём $v\in\mathcal N$ с $v^*v=e$ *$e$-изометрией*; тогда $vv^*$ —
проектор в $\mathcal N$. (Сначала $v=ve$: для любого $x$
$\|v(1-e)x\|^2=\langle(1-e)x,v^*v(1-e)x\rangle=\langle(1-e)x,e(1-e)x\rangle=0$. Отсюда
$vv^*v=ve=v$ и $(vv^*)^2=vev^*=vv^*$. Для $v_i$, введённых ниже, это даёт $e_iv_i=v_i$ и
$v_i^*e_i=v_i^*$.) Выберем
$e$-изометрии $v_1=e,v_2,\dots,v_k$ с попарно ортогональными образами $e_i=v_iv_i^*$, где $k$
максимально (образы — ортогональные подпространства одной и той же размерности $\operatorname{rank}e$,
поэтому $k\le\dim E/\operatorname{rank}e$). Утверждается, что $\sum_ie_i=1$. Иначе $g=1-\sum_ie_i\neq0$ — проектор в
$\mathcal N$. По шагу 2 $gAe\neq0$ для некоторого $A\in\mathcal N$ (если $gAe=0$ для всех $A$, то
$gAe\xi=0$ для всех $A,\xi$, и $g$ обращается в нуль на $F=E$). Положим $y=gAe$. Тогда $y^*y\in e\mathcal Ne$ положителен и отличен от нуля,
поэтому $y^*y=\lambda e$ с $\lambda>0$, и $v=y/\sqrt\lambda$ — $e$-изометрия с образом
в $gE$, ортогональным всем $e_i$, что противоречит максимальности.

При $i\neq j$ имеем $v_i^*v_j=v_i^*e_ie_jv_j=0$. Положим $E_{ij}=v_iv_j^*$. Тогда
$E_{ij}E_{kl}=v_iv_j^*v_kv_l^*=\delta_{jk}v_iev_l^*=\delta_{jk}E_{il}$,
$E_{ij}^*=E_{ji}$ и $\sum_iE_{ii}=1$. Для $A\in\mathcal N$
$A=\sum_{ij}e_iAe_j=\sum_{ij}v_i(v_i^*Av_j)v_j^*$ и $v_i^*Av_j\in e\mathcal Ne=\mathbb Ce$,
так что $A=\sum_{ij}a_{ij}E_{ij}$ с $a_{ij}\in\mathbb C$. Тем самым $\mathcal N=\operatorname{span}\{E_{ij}\}$.

*Шаг 4: унитарный оператор.* Пусть $m=\operatorname{rank}e$, а $f_1,\dots,f_m$ — ортонормированный базис
$eE$. Векторы $v_if_r$ ортонормированы:
$\langle v_if_r,v_jf_s\rangle=\langle f_r,v_i^*v_jf_s\rangle=\delta_{ij}\langle f_r,ef_s\rangle=\delta_{ij}\delta_{rs}$.
Они порождают $E$, поскольку $E=\bigoplus_ie_iE$ и $e_iE=v_ieE$. Положим $V(|i\rangle\otimes f_r)=v_if_r$, где $|1\rangle,\dots,|k\rangle$ —
стандартный базис $\mathbb C^k$, а $\mathbb C^m$ отождествлено с $eE$ посредством $f_r$.
Тогда $E_{ij}v_lf_r=\delta_{jl}v_if_r$, т. е. $V^*E_{ij}V=|i\rangle\langle j|\otimes1_m$, и
$V^*\mathcal NV=M_k(\mathbb C)\otimes1_m$.

*Шаг 5: коммутант.* Запишем $X\in M_k\otimes M_m$ как блочную матрицу
$X=\sum_{ij}|i\rangle\langle j|\otimes X_{ij}$. При фиксированных $a,b$
$(|a\rangle\langle b|\otimes1)X=\sum_j|a\rangle\langle j|\otimes X_{bj}$ и
$X(|a\rangle\langle b|\otimes1)=\sum_i|i\rangle\langle b|\otimes X_{ia}$. Если они равны, то
сравнение блоков в позиции $(a,j)$ при $j\neq b$ даёт $X_{bj}=0$, а в позиции
$(a,b)$ — $X_{bb}=X_{aa}$. Поскольку $a,b$ произвольны, $X=1_k\otimes Y$ с $Y=X_{11}$;
обратно, всякий $1_k\otimes Y$ коммутирует с $M_k\otimes1_m$. $\square$

**Теорема 2.6 (строение конечномерных алгебр фон Неймана).** Пусть
$\mathcal M$ — алгебра фон Неймана на $\mathbb C^n$. Существуют $C\ge1$, целые
$k_c,m_c\ge1$ с $\sum_ck_cm_c=n$ и унитарный оператор
$V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n$, такие, что
$$
V^*\mathcal MV=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},\qquad
V^*\mathcal M'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),
$$
и $V^*Z(\mathcal M)V=\operatorname{span}\{V^*P_1V,\dots,V^*P_CV\}$, где $P_1,\dots,P_C$ —
минимальные центральные проекторы алгебры $\mathcal M$ (лемма 2.3), а $V^*P_cV$ — проектор на
$c$-е слагаемое. Пары $(k_c,m_c)$ единственны с точностью до порядка: $P_c$ единственны (лемма 2.3),
$k_c^2=\dim\mathcal MP_c$ и $m_c=\operatorname{rank}P_c/k_c$. В частности, $\mathcal M\cong\bigoplus_cM_{k_c}(\mathbb C)$ как $*$-алгебры, и
$\dim\mathcal M=\sum_ck_c^2$, $\dim\mathcal M'=\sum_cm_c^2$.

*Доказательство.* По лемме 2.3 и замечанию после неё $\mathcal M=\bigoplus_c\mathcal M_c$, где
$\mathcal M_c$ действует на $E_c$, и каждая $\mathcal M_c$ — фактор (лемма 2.4). Лемма 2.5
даёт унитарные $V_c$ с $V_c^*\mathcal M_cV_c=M_{k_c}\otimes1_{m_c}$; положим
$V=\bigoplus_cV_c$. Для коммутанта: $X\in\mathcal M'$ коммутирует с каждым $P_c$, поэтому
$X=\sum_cP_cXP_c$, и $P_cXP_c$ коммутирует с $\mathcal M_c$; обратно, всякая такая сумма
коммутирует с $\mathcal M$. Лемма 2.5 описывает $\mathcal M_c'$ на $E_c$. $\square$

Числа $k_c$ описывают саму алгебру. **Кратности** $m_c$ описывают, как она
расположена в $M_n(\mathbb C)$: алгебраической структуре они не видны, но определяют
сужение следа пространства $\mathbb C^n$ на $\mathcal M$:
$$
\operatorname{Tr}_{\mathbb C^n}(A)=\sum_cm_c\operatorname{Tr}(A_c),\qquad A=V\Big(\bigoplus_cA_c\otimes1_{m_c}\Big)V^*.
$$
Поэтому библиотека работает с семейством следов (глава 3).

> **Физический смысл.** $P_c$ — это **суперотборные секторы**, или **секторы
> симметрии**: никакая наблюдаемая из $\mathcal M$ не связывает разные секторы, поэтому
> измерение «в каком секторе» совместимо с любым другим измерением. Центр — это
> *классическая* часть системы. Кратность $m_c$ — это вырождение: каждое
> собственное значение блока $A_c$ встречается в $\mathbb C^n$ с кратностью, в $m_c$ раз большей, чем в $A_c$.

## 2.4 Пример: группа симметрий и её коммутант

Пусть конечная группа $G$ действует на $\mathbb C^n$ посредством унитарного представления $\pi$.
Операторы, коммутирующие с симметрией, образуют $\pi(G)'$ — алгебру фон Неймана; её
коммутант $\pi(G)''$ есть линейная оболочка $\pi(g)$ (теорема 2.2, применённая к этой оболочке,
которая является $*$-алгеброй: она содержит $\pi(e)=1$, произведения $\pi(g)\pi(h)=\pi(gh)$ и сопряжённые
$\pi(g)^*=\pi(g^{-1})$).

**Предложение 2.7.** Применим теорему 2.6 к $\mathcal M=\pi(G)''$. Тогда
$V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$, где каждое $\rho_c$ — неприводимое
унитарное представление $G$ на $\mathbb C^{k_c}$, причём $\rho_c,\rho_d$ неэквивалентны при
$c\neq d$. Следовательно,
$$
\pi\cong\bigoplus_c\rho_c\otimes1_{m_c},\qquad
\pi(G)''\cong\bigoplus_cM_{\dim\rho_c},\qquad
\pi(G)'\cong\bigoplus_cM_{m_c}:
$$
размеры блоков симметричных операторов $\pi(G)'$ равны кратностям неприводимых
представлений, а роли $k_c$ и $m_c$ у алгебры и её коммутанта меняются местами.

*Доказательство.* $\pi(g)\in\mathcal M$, поэтому $V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$ для некоторых
$\rho_c(g)\in M_{k_c}$, и $g\mapsto\rho_c(g)$ — унитарное представление, поскольку таково $\pi$.
Так как $\pi(g)$ порождают $\mathcal M$ как линейное пространство, $\rho_c(g)$ порождают $M_{k_c}$; более того,
отображение $\sum_ga_g\pi(g)\mapsto\big(\sum_ga_g\rho_c(g)\big)_c$ сюръективно на $\bigoplus_cM_{k_c}$.
*Неприводимость:* подпространство $\mathbb C^{k_c}$, инвариантное относительно всех $\rho_c(g)$,
инвариантно и относительно их линейной оболочки $M_{k_c}$, а значит, равно $0$ или $\mathbb C^{k_c}$.
*Неэквивалентность:* если
$\rho_d(g)=T\rho_c(g)T^{-1}$ для всех $g$ и некоторого обратимого $T$, $c\neq d$, то каждый элемент
$(x_{c'})_{c'}$ образа удовлетворяет $x_d=Tx_cT^{-1}$, что исключает элемент с
$x_c=1$, $x_d=0$, — противоречие с сюръективностью. Последние два изоморфизма — это
теорема 2.6. $\square$

Разложение представления на неприводимые единственно с точностью до эквивалентности, так что
$m_c$ — *те самые* кратности в любом разложении. *Доказательство.* Пусть $\rho,\rho'$ неприводимы
и $T$ — сплетающий оператор, $T\rho(g)=\rho'(g)T$ для всех $g$. Ядро и образ $T$ инвариантны,
поэтому $T=0$ или $T$ обратим (лемма Шура); если $\rho'=\rho$, то же верно для $T-\lambda1$, где
$\lambda$ — собственное значение $T$, так что $T=\lambda1$. Следовательно, сплетающие операторы из
$\rho$ в $\rho'$ образуют пространство размерности $1$, если $\rho\cong\rho'$, и $0$ иначе. Если
$\pi\cong\bigoplus_j\sigma_j$ с неприводимыми $\sigma_j$ (повторения допускаются), то сплетающие
операторы из неприводимого $\rho$ в $\pi$ образуют пространство размерности
$\#\{j:\sigma_j\cong\rho\}$ — числа, зависящего только от $\pi$ и $\rho$. Для разложения из
предложения 2.7, где $\rho_c$ попарно неэквивалентны, оно равно $m_c$ при $\rho\cong\rho_c$ и $0$,
если $\rho$ не эквивалентно ни одному $\rho_c$. $\square$

Возьмём $G=S_3$, переставляющую три сомножителя $(\mathbb C^2)^{\otimes3}$, т. е. три
кубита или три спина $\tfrac12$. Разложение таково: симметрическое подпространство (полный спин
$\tfrac32$, размерность $4$) несёт тривиальное представление с кратностью $4$; его дополнение
(размерность $4$) несёт двумерное стандартное представление с кратностью $2$; знаковое
представление не встречается. Поэтому $\pi(S_3)''\cong M_1\oplus M_2$ (размерность $5$) и
$\pi(S_3)'\cong M_4\oplus M_2$ (размерность $20$).

*Доказательство.* Запишем $\pi(S_3)''\cong\bigoplus_cM_{k_c}\otimes1_{m_c}$ (теорема 2.6). Тогда
$\sum_ck_cm_c=8$, $\sum_ck_c^2=\dim\pi(S_3)''=5$ и $\sum_cm_c^2=\dim\pi(S_3)'=20$; обе размерности
код ниже вычисляет точными методами линейной алгебры (коммутант множества $S$ — ядро линейного
отображения $X\mapsto(AX-XA)_{A\in S}$). Число $20$ можно подсчитать и вручную: $X$ коммутирует со
всеми $\pi(g)$ тогда и только тогда, когда его элементы $X_{xy}$, $x,y\in\{0,1\}^3$, постоянны на
орбитах действия $S_3$ на парах $(x,y)$, т. е. на тройках элементов $\{0,1\}^2$, а таких орбит
$\binom63=20$ (мультимножества размера $3$ из $4$ элементов). Из $\sum_ck_c^2=5$ следует $k=(1,2)$
или $k=(1,1,1,1,1)$. В первом случае $m_1+2m_2=8$ и $m_1^2+m_2^2=20$ дают $m=(4,2)$. Во втором
пять натуральных чисел с суммой $8$ и суммой квадратов $20$ обязаны быть $m=(4,1,1,1,1)$; но тогда
$\pi(S_3)''$ была бы коммутативной, тогда как транспозиции $(12)$ и $(23)$ не коммутируют:
$\pi((12))\pi((23))|001\rangle=|100\rangle$, но $\pi((23))\pi((12))|001\rangle=|010\rangle$. Итак,
$k=(1,2)$ и $m=(4,2)$. Одномерное представление $\rho_1$ — тривиальное или знаковое. Оно не
знаковое: каждый базисный вектор $|b_1b_2b_3\rangle$ имеет два равных бита, поэтому неподвижен
относительно некоторой транспозиции $t$, и $A=\sum_g\operatorname{sgn}(g)\pi(g)$, для которого
$A\pi(t)=-A$, даёт $A|b_1b_2b_3\rangle=A\pi(t)|b_1b_2b_3\rangle=-A|b_1b_2b_3\rangle$; значит,
$A=0$, и знаковое представление не встречается. Итак, $\rho_1$ тривиально, его $4$ копии
порождают симметрическое подпространство, а $\rho_2$ — неприводимое двумерное представление
(предложение 2.7), т. е. стандартное. $\square$

```python
import itertools
import torch

torch.manual_seed(0)
dt = torch.complex128

def permutation_operator(perm, n_qubits=3):
    """Unitary on (C^2)^{(x)3} that permutes the tensor factors."""
    dim = 2 ** n_qubits
    P = torch.zeros(dim, dim, dtype=dt)
    for bits in itertools.product([0, 1], repeat=n_qubits):
        src = int(''.join(map(str, bits)), 2)
        dst = int(''.join(str(bits[perm[i]]) for i in range(n_qubits)), 2)
        P[dst, src] = 1
    return P

group = [permutation_operator(p) for p in itertools.permutations(range(3))]

def commutant_basis(ops, tol=1e-10):
    """Orthonormal basis (Hilbert-Schmidt) of {X : X A = A X for all A in ops}."""
    d = ops[0].shape[0]
    eye = torch.eye(d, dtype=dt)
    # vec(AX - XA) = (A (x) 1 - 1 (x) A^T) vec(X) for row-major vec
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    _, s, Vh = torch.linalg.svd(L)
    null = Vh[(s > tol).sum():].conj().resolve_conj()   # rows span the null space
    return [v.reshape(d, d) for v in null]

def span_dim(ops, tol=1e-10):
    return int((torch.linalg.svdvals(torch.stack([A.reshape(-1) for A in ops])) > tol).sum())

comm = commutant_basis(group)
print(len(comm), span_dim(group), len(commutant_basis(comm)))
assert len(comm) == 20                                   # pi(S3)'  = M_4 (+) M_2
assert span_dim(group) == 5                              # pi(S3)'' = M_1 (+) M_2
assert len(commutant_basis(comm)) == 5                   # double commutant theorem
```

Кратности видны в спектре. Усреднение по группе
$X\mapsto\frac1{|G|}\sum_g\pi(g)X\pi(g)^*$ переводит любой оператор в коммутант
(упражнение 3). Поэтому случайный самосопряжённый элемент $\pi(S_3)'$ (почти наверное) имеет $4$ простых
собственных значения на симметрическом подпространстве и $2$ двукратно вырожденных на его
дополнении, поскольку блок $\pi(S_3)'$ там равен $1_2\otimes M_2$; код проверяет вырожденные пары.

```python
def group_average(X):
    return sum(g @ X @ g.conj().T for g in group) / len(group)

H0 = torch.randn(8, 8, dtype=dt)
Hsym = group_average(H0 + H0.conj().T)                  # a random symmetric Hamiltonian
assert all(torch.allclose(g @ Hsym, Hsym @ g) for g in group)

P_sym = sum(group) / len(group)                         # projection onto the symmetric subspace
w, Q = torch.linalg.eigh(P_sym)
Q_rest = Q[:, w.real < 0.5]                             # orthonormal basis of the complement (dim 4)
e_rest = torch.linalg.eigvalsh(Q_rest.conj().T @ Hsym @ Q_rest)
print(e_rest)
assert Q_rest.shape[1] == 4
assert torch.allclose(e_rest[0], e_rest[1]) and torch.allclose(e_rest[2], e_rest[3])
```

## 2.5 Алгебра в библиотеке

Библиотека хранит абстрактную алгебру $\bigoplus_cM_{k_c}(\mathbb C)$ из теоремы 2.6 на
$\bigoplus_c\mathbb C^{k_c}$, без кратностей: оператор — это список его
блоков $A_c$, дополненных до общего размера. Коммутант из примера — это
`TypeIAlgebra([4, 2], [4, 2])`.

```python
from torch_vn_algebra import TypeIAlgebra, channels

alg = TypeIAlgebra([4, 2], [4, 2], precision='double')    # M_4 (+) M_2
print(alg.k_factors, alg.C)

A = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
print(A.matrix.shape)              # (1, 2, 4, 4): block 1 occupies the top-left 2 x 2 corner
assert A.matrix[0, 1, 2:, :].abs().max() == 0 and A.matrix[0, 1, :, 2:].abs().max() == 0
```

Центр состоит из элементов $\sum_cz_cP_c$; их строит `alg.central`, и они
коммутируют со всем:

```python
def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

P0, P1 = alg.central([1.0, 0.0]), alg.central([0.0, 1.0])   # minimal central projections
Z = alg.central([2.0, -1.0])
assert P0.is_projection and P1.is_projection
assert close(P0 + P1, alg.identity())
assert close(Z @ A, A @ Z)
assert close(A @ P0 + A @ P1, A)                         # A = sum_c A P_c
```

## 2.6 Условные ожидания

**Условным ожиданием** из $\mathcal M$ на подалгебру фон Неймана
$\mathcal N\subseteq\mathcal M$ называется линейное отображение $E:\mathcal M\to\mathcal N$, которое положительно
($A\ge0\Rightarrow E(A)\ge0$), удовлетворяет $E(B)=B$ для $B\in\mathcal N$ и является
гомоморфизмом $\mathcal N$-бимодулей: $E(B_1AB_2)=B_1E(A)B_2$ для $B_1,B_2\in\mathcal N$.

**Предложение 2.8.** Для $\mathcal M=\bigoplus_cM_{k_c}$ и $\mathcal N=Z(\mathcal M)$ отображение
$$
E(A)=\sum_c\frac{\operatorname{Tr}A_c}{k_c}\,P_c
$$
является условным ожиданием. Оно сохраняет всякий след вида
$\sum_cw_c\operatorname{Tr}A_c$.

*Доказательство.* Линейность очевидна; из $A\ge0$ следует $\operatorname{Tr}A_c\ge0$, поэтому $E(A)\ge0$. Для
$B=\sum_cb_cP_c$ имеем $E(B)=\sum_cb_ck_cP_c/k_c=B$. Для $B_i=\sum_cb^{(i)}_cP_c$
$\operatorname{Tr}(B_1AB_2)_c=b^{(1)}_cb^{(2)}_c\operatorname{Tr}A_c$, что даёт свойство
бимодульности. Наконец, $\operatorname{Tr}(E(A)_c)=\operatorname{Tr}A_c$ для каждого $c$. $\square$

В библиотеке это отображение — канал `channels.center_expectation(alg)`. Применяя его
к операторам, получаем классическую (усреднённую по секторам) часть:

```python
E = channels.center_expectation(alg)
EA = E(A)
assert close(E(EA), EA)                                  # idempotent
assert close(E(Z @ A @ P1), Z @ EA @ P1)                 # bimodule property
A2 = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
assert close(EA @ A2, A2 @ EA)                           # E(A) commutes with a fresh element
for c in range(alg.C):                                   # Tr E(A)_c = Tr A_c
    assert abs(EA.matrix[0, c].diagonal().sum() - A.matrix[0, c].diagonal().sum()) < 1e-10
```

Усреднение по группе из раздела 2.4 — это условное ожидание на коммутант
$\pi(G)'$ внутри $M_n(\mathbb C)$, и оно сохраняет $\operatorname{Tr}$ (упражнение 3).
(В общем случае условные ожидания не единственны: на $\mathbb C1$ каждое состояние
$\varphi$ задаёт своё, $A\mapsto\varphi(A)1$.) Оба примера воплощают одну и ту же идею:
проекцию на операторы, которые «не видят» некоторой структуры.

> **Физический смысл.** На полной матричной алгебре $M_D(\mathbb C)$, $D=\sum_ck_c$, отображение
> $X\mapsto\sum_cP_cXP_c$ — условное ожидание на $\mathcal M$: оно стирает когерентности между
> суперотборными секторами, что и происходит, когда заряд (сектор) измеряют, не записывая
> результат. Условное ожидание $E$ на центр идёт дальше: оно сохраняет лишь классическую
> информацию — вес $p_c$ каждого сектора — и заменяет состояние внутри каждого сектора
> максимально смешанным. Состояние $\omega$ и состояние $\omega\circ E$ совпадают на всякой
> наблюдаемой, которая спрашивает только «в каком секторе?».

## 2.7 С чем сравнивать

Примеры раздела 2.4 можно проверить в любом пакете линейной алгебры;
программы для теории представлений (GAP, SageMath) дают разложение
$\pi\cong\bigoplus\rho_\lambda\otimes1_{m_\lambda}$ в символьном виде. Физические программы
точной диагонализации (QuSpin, модули симметрий QuTiP) используют ту же блочную структуру,
чтобы диагонализовать гамильтонианы посекторно. Библиотека берёт эту структуру за
отправную точку: когда блоки известны, все дальнейшие вычисления — функции,
следы, состояния, каналы, динамика — остаются внутри алгебры.

## Упражнения

1. Покажите, что центр $M_k(\mathbb C)$ равен $\mathbb C1$ (воспользуйтесь матричными единицами
   $|i\rangle\langle j|$), и что $Z(\bigoplus_cM_{k_c})$ имеет размерность $C$.
2. Пусть $S$ — циклический сдвиг на $\mathbb C^3$, $Se_i=e_{i+1\bmod3}$. Найдите
   $\{S\}'$ и $\{S\}''$ и проверьте размерности численно.
3. Для унитарного представления $\pi$ конечной группы покажите, что
   $E(X)=\frac1{|G|}\sum_g\pi(g)X\pi(g)^*$ — условное ожидание
   $M_n(\mathbb C)$ на $\pi(G)'$, сохраняющее $\operatorname{Tr}$.

<details>
<summary>Решения</summary>

1. Если $X$ коммутирует со всеми $|i\rangle\langle j|$, то
   $X|i\rangle\langle j|=|i\rangle\langle j|X$; применяя обе части к $e_j$ и сравнивая,
   получаем $Xe_i=X_{jj}e_i$ для всех $i,j$, так что $X=X_{11}1$. В прямой сумме элемент
   $\bigoplus X_c$ централен тогда и только тогда, когда централен каждый $X_c$, поэтому центр равен
   $\{\bigoplus_cz_c1_{k_c}\}$ и имеет размерность $C$.

2. $S$ унитарен с собственными значениями $1,\omega,\omega^2$, $\omega=e^{2\pi i/3}$, каждое из которых простое.
   Оператор, коммутирующий с $S$, сохраняет его собственные подпространства, поэтому в
   собственном базисе он диагонален: $\{S\}'\cong\mathbb C^3$ — это циркулянтные матрицы. Тогда и
   $\{S\}''\cong\mathbb C^3$, поскольку в собственном базисе $S$ матрица, коммутирующая со всеми
   диагональными матрицами, диагональна (сравните элементы $XD$ и $DX$ для $D$ с попарно различными
   диагональными элементами); она порождается $1,S,S^2$ (определитель Вандермонда: собственные значения различны).

```python
S = torch.roll(torch.eye(3, dtype=dt), 1, dims=0)
assert len(commutant_basis([S])) == 3
assert len(commutant_basis(commutant_basis([S]))) == 3
assert span_dim([torch.eye(3, dtype=dt), S, S @ S]) == 3
```

3. $E$ линейно и положительно (это сумма отображений $X\mapsto UXU^*$). Для $h\in G$
   $\pi(h)E(X)\pi(h)^*=\frac1{|G|}\sum_g\pi(hg)X\pi(hg)^*=E(X)$, поэтому $E(X)\in\pi(G)'$. Если
   $X\in\pi(G)'$, то каждое слагаемое равно $X$, так что $E(X)=X$. Для $B_1,B_2\in\pi(G)'$
   $\pi(g)B_1XB_2\pi(g)^*=B_1\pi(g)X\pi(g)^*B_2$, что даёт свойство бимодульности.
   Наконец, $\operatorname{Tr}(\pi(g)X\pi(g)^*)=\operatorname{Tr}X$ для каждого $g$, поэтому
   $\operatorname{Tr}E(X)=\operatorname{Tr}X$.

</details>

**Далее:** в [главе 3](03_traces_and_states.md) классифицируются следы на
$\bigoplus_cM_{k_c}$, вводятся состояния и плотности и доказываются энтропийные и
следовые неравенства, используемые в дальнейшем.
