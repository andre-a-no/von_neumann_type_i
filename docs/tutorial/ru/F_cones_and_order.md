[English](../F_cones_and_order.md) · [Español](../es/F_cones_and_order.md) · [Français](../fr/F_cones_and_order.md) · [Deutsch](../de/F_cones_and_order.md) · [中文](../zh/F_cones_and_order.md) · [日本語](../ja/F_cones_and_order.md) · **Русский**

# Приложение F. Симметрические конусы, порядковые изоморфизмы и йордановы отображения

В главе 13 порядок на наблюдаемых был связан с их йордановым произведением тремя теоремами: теоремой
Кёхера–Винберга (теорема 13.8), теоремой Кадисона (теорема 13.10) и описанием йордановых
автоморфизмов $M_n(\mathbb C)$ (предложение 13.11 доказывает его при $n=2$). В этом приложении первая доказывается в полной
общности, а две другие — для алгебр $\bigoplus_cM_{k_c}(\mathbb C)$ курса при всех размерах блоков. Попутно
доказываются факты о евклидовых йордановых алгебрах, нужные в разделе 13.6: спектральная теорема,
самодвойственность конуса квадратов и $P(y)Q\subseteq Q$.

**Вы узнаете**, почему всякая евклидова йорданова алгебра степенно ассоциативна и обладает спектральной теоремой; почему её конус
квадратов самодвойствен и однороден; как йорданово произведение восстанавливается по одному лишь самодвойственному однородному конусу —
через третью производную логарифма его характеристической функции; почему унитальная линейная биекция, сохраняющая
положительность в обе стороны, сохраняет квадраты; и почему всякое такое отображение $\bigoplus_cM_{k_c}$ переставляет
одинаковые блоки и действует на каждом блоке как $UxU^*$ или как $Ux^TU^*$.

Помимо линейной алгебры, анализа и результатов глав 1, 6 и 13, на которые мы ссылаемся по номерам, в разделе F.3 используется интегрирование
по Лебегу (лемма Фату, мажорированная сходимость и дифференцирование под знаком интеграла, теорема Фубини–Тонелли,
формула линейной замены переменных); больше ничего не предполагается. Спектральная теорема из приложения E с йордановыми реперами уточняет
теорему F.3 ниже; ничто здесь от неё не зависит.
**Евклидово пространство** — это конечномерное вещественное векторное пространство $V\ne\{0\}$ со скалярным произведением
$\langle\cdot,\cdot\rangle$; $\|x\|$ — норма (одинарные черты $|x|$ для элемента $V$ будут обозначать абсолютную
величину из теоремы F.3, но никогда не норму), $\|A\|$ — операторная норма линейного отображения $A$ пространства $V$, $A^T$ — его сопряжённое
($\langle Ax,y\rangle=\langle x,A^Ty\rangle$), а $A\ge0$ ($A>0$) означает, что $A=A^T$ положительно полуопределён
(определён). $e^A=\exp A=\sum_{k\ge0}A^k/k!$.

## F.1 Евклидовы йордановы алгебры

Как и в разделе 13.1, **йорданова алгебра** — это вещественное векторное пространство с билинейным коммутативным произведением $\circ$, таким, что
$x\circ(x^2\circ y)=x^2\circ(x\circ y)$, т. е. $[L(x),L(x^2)]=0$, где $L(x)y=x\circ y$ и $x^2=x\circ x$; кроме того,
$P(x)=2L(x)^2-L(x^2)$.

**Определение.** **Евклидова йорданова алгебра** — это конечномерная йорданова алгебра $V$ с единицей $e$
($e\circ x=x$) и скалярным произведением, относительно которого каждый оператор $L(x)$ симметричен: $\langle x\circ y,z\rangle=\langle y,x\circ z\rangle$.
Её **конус квадратов** — это $Q=\{x^2:x\in V\}$.

Примеры: $M_{sa}$ (самосопряжённая часть $\bigoplus_cM_{k_c}(\mathbb C)$) и $\operatorname{Herm}_n(\mathbb R)$
с $\langle x,y\rangle=\operatorname{Tr}(xy)$, которое вещественно и удовлетворяет $\operatorname{Tr}(x^2)=\sum_{ij}|x_{ij}|^2$ (лемма 13.1;
$\operatorname{Tr}((x\circ y)z)=\frac12\operatorname{Tr}(xyz+yxz)=\operatorname{Re}\operatorname{Tr}(xyz)$, поскольку
$\operatorname{Tr}(yxz)=\overline{\operatorname{Tr}((yxz)^*)}=\overline{\operatorname{Tr}(zxy)}=\overline{\operatorname{Tr}(xyz)}$,
инвариантно относительно циклических перестановок и по той же причине относительно $x\leftrightarrow z$, а значит, симметрично по
$x,y,z$); спин-факторы $V_m$ (йордановы алгебры по лемме 13.3) со скалярным произведением
$\langle x,y\rangle=\operatorname{tr}(x\circ y)=2(ab+u\cdot v)$ для $x=(a,u)$, $y=(b,v)$, поскольку
$\operatorname{tr}((x\circ y)\circ z)=2(abc+c\,u\cdot v+a\,v\cdot w+b\,u\cdot w)$ для $z=(c,w)$ симметрично по
$x,y,z$; и $\mathbb R^N$ с покоординатным произведением (упражнение 1).

**Лемма F.1.** В любой йордановой алгебре для всех $a,b,d,x,y$:

1. $[L(a),L(b\circ d)]+[L(b),L(d\circ a)]+[L(d),L(a\circ b)]=0$;
2. $L(x^2\circ y)=L(x^2)L(y)-2L(x)L(y)L(x)+2L(x\circ y)L(x)$.

*Доказательство.* (1) Для $s,t,r\in\mathbb R$ и $x=sa+tb+rd$ равенство $[L(x),L(x^2)]=0$ — полиномиальное тождество по $(s,t,r)$, поэтому все
коэффициенты равны нулю; коэффициент при $str$ равен удвоенной левой части (1). (2) Положим $a=b=x$ в (1):
$2[L(x),L(x\circ d)]+[L(d),L(x^2)]=0$. Применённое к $y$, это означает
$2x\circ((x\circ d)\circ y)-2(x\circ d)\circ(x\circ y)+d\circ(x^2\circ y)-x^2\circ(d\circ y)=0$, т. е.
$\big(2L(x)L(y)L(x)-2L(x\circ y)L(x)+L(x^2\circ y)-L(x^2)L(y)\big)d=0$ для всех $d$. $\square$

**Предложение F.2 (степенная ассоциативность).** Пусть $V$ — йорданова алгебра с единицей $e$, $x\in V$, $x^0=e$ и
$x^{k+1}=x\circ x^k$. Тогда $x^p\circ x^q=x^{p+q}$ для всех $p,q\ge0$, и каждый оператор $L(x^p)$ — многочлен от
коммутирующих операторов $L(x)$ и $L(x^2)$. Следовательно, $p(x)\circ q(x)=(pq)(x)$ для вещественных многочленов $p,q$, и операторы
$L(p(x))$ коммутируют друг с другом.

*Доказательство.* $L(x)$ и $L(x^2)$ коммутируют в силу йорданова тождества. Докажем индукцией по $N\ge2$: (A$_N$)
$x^i\circ x^j=x^{i+j}$ при $i+j\le N$ и (B$_N$) $L(x^i)$ — многочлен от $L(x),L(x^2)$ при $i\le N$. Оба утверждения
очевидны при $N=2$. Предположим их для $N$, и пусть $i+j=N+1$ с $i,j\ge1$ (при $i=0$ или $j=0$ доказывать нечего).
Поскольку $i\le N$, $L(x^i)$ коммутирует с $L(x)$, так что по (A$_N$)
$$
x^i\circ x^j=L(x^i)L(x)x^{j-1}=L(x)L(x^i)x^{j-1}=x\circ x^N=x^{N+1}.
$$
Это (A$_{N+1}$). Лемма F.1(2) с $y=x^{N-1}$ и $x^2\circ x^{N-1}=x^{N+1}$ даёт
$L(x^{N+1})=L(x^2)L(x^{N-1})-2L(x)L(x^{N-1})L(x)+2L(x^N)L(x)$, что и есть (B$_{N+1}$). Последнее предложение следует из
билинейности. $\square$

**Теорема F.3 (спектральная теорема).** Пусть $V$ — евклидова йорданова алгебра и $x\in V$. Существуют единственные
$k\ge1$, вещественные числа $\lambda_1>\dots>\lambda_k$ и ненулевые $c_1,\dots,c_k\in V$, для которых
$$
c_i\circ c_j=\delta_{ij}c_i,\qquad \textstyle\sum_ic_i=e,\qquad x=\sum_i\lambda_ic_i .
$$
Элементы $c_i$ — многочлены от $x$, и $p(x)=\sum_ip(\lambda_i)c_i$ для любого вещественного многочлена $p$.

Мы называем $(c_i)$ **полной системой ортогональных идемпотентов**, а $\lambda_i$ — **собственными значениями** $x$, и полагаем
$f(x)=\sum_if(\lambda_i)c_i$ для $f:\{\lambda_1,\dots,\lambda_k\}\to\mathbb R$; например, $|x|$, $x_\pm=\max(\pm x,0)$,
$\sqrt x$ (если все $\lambda_i\ge0$), $\exp x$ и $\log x$ (если все $\lambda_i>0$).

*Доказательство.* Существование. $\mathbb R[x]=\{p(x)\}$ — подпространство, которое содержит $e$ и по предложению F.2 удовлетворяет
$p(x)\circ q(x)=(pq)(x)$. Поэтому $M=L(x)|_{\mathbb R[x]}$ переводит $\mathbb R[x]$ в себя, симметричен и удовлетворяет
$p(M)q(x)=p(x)\circ q(x)$. Симметричный оператор диагонализуем с вещественными собственными значениями; пусть $\lambda_1>\dots>\lambda_k$
— собственные значения $M$, а $m(t)=\prod_i(t-\lambda_i)$ — его минимальный многочлен. Тогда $m(x)=m(M)e=0$. Для интерполяционных
многочленов Лагранжа $\ell_i(t)=\prod_{j\ne i}\frac{t-\lambda_j}{\lambda_i-\lambda_j}$ многочлены
$\ell_i\ell_j-\delta_{ij}\ell_i$, $\sum_i\ell_i-1$ и $\sum_i\lambda_i\ell_i-t$ обращаются в нуль в каждой точке $\lambda_l$, поэтому они
кратны $m$; подстановка $x$ показывает, что $c_i=\ell_i(x)$ обладают тремя указанными свойствами. Если $c_i=0$, то
$\ell_i(M)q(x)=\ell_i(x)\circ q(x)=0$ для всех $q$, так что $\ell_i(M)=0$, хотя $\deg\ell_i<\deg m$; следовательно, $c_i\ne0$.
Наконец, $p(x)=\sum_ip(\lambda_i)c_i$: это верно для $p=1$ и $p=t$, и $x\circ\sum_i\lambda_i^rc_i=\sum_i\lambda_i^{r+1}c_i$.

Единственность. Пусть $x=\sum_j\mu_jd_j$ с различными $\mu_j$ и полной системой ненулевых ортогональных идемпотентов $d_j$.
Как и выше, $p(x)=\sum_jp(\mu_j)d_j$. Элементы $d_j$ линейно независимы (умножьте $\sum_ja_jd_j=0$ на $d_i$), поэтому
$p(x)=0$ тогда и только тогда, когда $p(\mu_j)=0$ для всех $j$. Значит, из $m(x)=0$ следует $\{\mu_j\}\subseteq\{\lambda_i\}$, а
$\tilde m(t)=\prod_j(t-\mu_j)$ удовлетворяет $\tilde m(M)q(x)=\tilde m(x)\circ q(x)=0$, так что $m$ делит $\tilde m$ и $\{\lambda_i\}\subseteq\{\mu_j\}$.
Тогда $d_j=\ell_j(x)=c_j$. $\square$

**Лемма F.4 (разложение Пирса).** Пусть $V$ — евклидова йорданова алгебра и $c$ — идемпотент ($c^2=c$).

1. $L(c)(2L(c)-I)(L(c)-I)=0$. Следовательно, $V$ — ортогональная сумма собственных подпространств $V(c,\nu)=\{z:c\circ z=\nu z\}$,
   $\nu\in\{0,\frac12,1\}$, и $L(c)\ge0$.
2. $V(c,1)\circ V(c,0)=\{0\}$; $V(c,\nu)\circ V(c,\nu)\subseteq V(c,\nu)$ и $V(c,\nu)\circ V(c,\frac12)\subseteq V(c,\frac12)$
   для $\nu\in\{0,1\}$; $V(c,\frac12)\circ V(c,\frac12)\subseteq V(c,0)+V(c,1)$.
3. $s=2c-e$ удовлетворяет $s^2=e$, и $P(s)$ — йорданов автоморфизм $V$: он равен $+1$ на $V(c,0)+V(c,1)$, $-1$ на
   $V(c,\frac12)$, и $P(s)(z\circ w)=P(s)z\circ P(s)w$.
4. Пусть $x=\sum_{i=1}^k\lambda_ic_i$, как в теореме F.3. Операторы $L(c_i)$ коммутируют, и $V$ — ортогональная сумма
   подпространств $V_{ii}=\{z:c_i\circ z=z\}$ и $V_{ij}=\{z:c_i\circ z=c_j\circ z=\frac12z\}$ ($i<j$). Для
   $a=\sum_la_lc_l$ ($a_l\in\mathbb R$) оператор $L(a)$ действует на $V_{ij}$ ($i\le j$) как скаляр $\frac12(a_i+a_j)$, а $P(a)$ —
   как скаляр $a_ia_j$.

*Доказательство.* (1) Лемма F.1(2) с $x=y=c$ даёт $L(c)=L(c)^2-2L(c)^3+2L(c)^2$, т. е. $2L(c)^3-3L(c)^2+L(c)=0$, а это и есть
указанное разложение на множители. Поэтому собственные значения симметричного оператора $L(c)$ лежат в $\{0,\frac12,1\}$.

(2) Лемма F.1(1) с $a=b=c$, $d=z\in V(c,\nu)$ даёт $2\nu[L(c),L(z)]+[L(z),L(c)]=0$, так что $L(z)$ коммутирует с $L(c)$
при $\nu\in\{0,1\}$ и тогда сохраняет каждое $V(c,\mu)$. Это даёт второе и третье включения, а также
$V(c,1)\circ V(c,0)\subseteq V(c,0)\cap V(c,1)=\{0\}$. Для $z,w\in V(c,\frac12)$ применим лемму F.1(1) с $(a,b,d)=(z,w,c)$
к $c$: так как $L(w\circ c)=\frac12L(w)$ и $L(c\circ z)=\frac12L(z)$, первые два слагаемых сокращаются, и
$[L(c),L(u)]c=L(c)^2u-L(c)u=0$ для $u=z\circ w$. Если записать $u=u_1+u_{1/2}+u_0$ согласно (1), это означает $-\frac14u_{1/2}=0$.

(3) $s^2=4c-4c+e=e$, оператор $L(s)=2L(c)-I$ равен $1,0,-1$ на $V(c,1),V(c,\frac12),V(c,0)$, и $P(s)=2L(s)^2-I$. Запишем
$z=z'+z''$, $w=w'+w''$ с $z',w'\in V(c,0)+V(c,1)$ и $z'',w''\in V(c,\frac12)$. По (2) $z'\circ w'+z''\circ w''$ лежит в
$V(c,0)+V(c,1)$, а $z'\circ w''+z''\circ w'$ — в $V(c,\frac12)$, поэтому
$P(s)(z\circ w)=(z'-z'')\circ(w'-w'')=P(s)z\circ P(s)w$.

(4) Элементы $c_l$ — многочлены от $x$, поэтому операторы $L(c_l)$ коммутируют (предложение F.2); будучи симметричными, они обладают общим
ортонормированным собственным базисом. На общем собственном векторе $L(c_l)$ имеет собственное значение $\varepsilon_l\in\{0,\frac12,1\}$ по (1), и
$\sum_l\varepsilon_l=1$, поскольку $\sum_lL(c_l)=L(e)=I$. Поэтому либо одно $\varepsilon_i=1$, а остальные равны нулю, либо
$\varepsilon_i=\varepsilon_j=\frac12$ для одной пары $i<j$, а остальные равны нулю. Итак, $V$ — ортогональная сумма
общих собственных подпространств $W_{ii}$ (где $L(c_i)=1$, $L(c_l)=0$ при $l\ne i$) и $W_{ij}$ (где $L(c_i)=L(c_j)=\frac12$,
$L(c_l)=0$ в остальных случаях), и $W_{ij}=V_{ij}$: вектор с $c_i\circ z=z$ (соответственно $c_i\circ z=c_j\circ z=\frac12z$) не имеет
компонент в других общих собственных подпространствах, на которых $L(c_i)$ (соответственно пара $(L(c_i),L(c_j))$) принимает другое
значение. Поскольку
$a^2=\sum_la_l^2c_l$, на $V_{ij}$ ($i<j$) имеем $L(a)=\frac12(a_i+a_j)$, $L(a^2)=\frac12(a_i^2+a_j^2)$ и
$P(a)=\frac12(a_i+a_j)^2-\frac12(a_i^2+a_j^2)=a_ia_j$; на $V_{ii}$ имеем $L(a)=a_i$ и $P(a)=2a_i^2-a_i^2$. $\square$

## F.2 Конус квадратов симметричен

Пусть $\Omega$ — непустой открытый **выпуклый конус** в евклидовом пространстве $V$ ($tx\in\Omega$ для $x\in\Omega$, $t>0$). Его
(открытый) **двойственный конус** — это $\Omega^*=\{y:\langle x,y\rangle>0\text{ для всех }x\in\bar\Omega\setminus\{0\}\}$; $\Omega$
**самодвойствен**, если $\Omega^*=\Omega$. Положим $G(\Omega)=\{g\in GL(V):g\Omega=\Omega\}$; $\Omega$ **однороден**, если
$G(\Omega)$ действует на $\Omega$ транзитивно. **Симметрический конус** — это самодвойственный однородный конус. Если $\Omega_1\subseteq\Omega_2$,
то $\Omega_2^*\subseteq\Omega_1^*$.

**Теорема F.5 (самодвойственность).** Пусть $V$ — евклидова йорданова алгебра и $y=\sum_i\lambda_ic_i$ (теорема F.3).

1. Следующие условия равносильны: (i) $y\in Q$; (ii) все $\lambda_i\ge0$; (iii) $L(y)\ge0$; (iv) $\langle y,x\rangle\ge0$
   для всех $x\in Q$. В частности, $Q$ — замкнутый выпуклый конус.
2. Следующие условия равносильны: (i) $y\in\Omega:=\operatorname{int}Q$; (ii) все $\lambda_i>0$; (iii) $L(y)>0$; (iv)
   $\langle y,x\rangle>0$ для всех $x\in Q\setminus\{0\}$. Кроме того, $e\in\Omega$, $\bar\Omega=Q$ и $\Omega^*=\Omega$:
   $\Omega$ — самодвойственный открытый выпуклый конус.
3. Если $x,y\in Q$ и $\langle x,y\rangle=0$, то $x\circ y=0$.
4. Из $x^2=0$ следует $x=0$, и $Q\cap(-Q)=\{0\}$.

*Доказательство.* Мы используем $\langle y,x^2\rangle=\langle y\circ x,x\rangle=\langle L(y)x,x\rangle$ и
$\langle x^2,e\rangle=\langle x,x\circ e\rangle=\|x\|^2$; последнее вместе с (1)(ii) даёт (4), поскольку
собственные значения $-y$ — это $-\lambda_i$ в силу единственности.

(1) (i)$\Rightarrow$(iii): если $y=x^2$ с $x=\sum_i\mu_id_i$, то $y=\sum_i\mu_i^2d_i$ и
$L(y)=\sum_i\mu_i^2L(d_i)\ge0$ по лемме F.4(1). (iii)$\Leftrightarrow$(iv) в силу первого тождества.
(iii)$\Rightarrow$(ii): $0\le\langle L(y)c_i,c_i\rangle=\lambda_i\|c_i\|^2$. (ii)$\Rightarrow$(i):
$y=(\sum_i\sqrt{\lambda_i}c_i)^2$. По (iv) $Q$ — пересечение замкнутых полупространств, проходящих через $0$.

(2) (ii)$\Rightarrow$(iii): $L(y)=\sum_i\lambda_iL(c_i)\ge(\min_i\lambda_i)\sum_iL(c_i)=(\min_i\lambda_i)I$.
(iii)$\Rightarrow$(iv): ненулевой $x\in Q$ равен $z^2$ с $z\ne0$, и $\langle y,x\rangle=\langle L(y)z,z\rangle>0$.
(iv)$\Rightarrow$(ii): $\langle y,c_i\rangle=\lambda_i\|c_i\|^2$ и $c_i=c_i^2\in Q\setminus\{0\}$. (iii)$\Rightarrow$(i): если
$L(y)\ge mI$ с $m>0$, то $L(y')\ge(m-\|L\|\,\|y-y'\|)I$, где $\|L\|$ — норма линейного отображения
$b\mapsto L(b)$, так что некоторая окрестность $y$ состоит из элементов с
$L>0$, которые лежат в $Q$ по (1). (i)$\Rightarrow$(iii): $L(y)\ge0$ по (1); если $\langle L(y)z,z\rangle=0$ при $z\ne0$, то
$x=z^2\ne0$ и $\langle y,x\rangle=0$; при малых $\varepsilon>0$ имеем $y-\varepsilon x\in Q$, но
$\langle y-\varepsilon x,x\rangle=-\varepsilon\|x\|^2<0$, что противоречит (1)(iv). Далее, $L(e)=I>0$, и для $y\in Q$
элементы $y+\varepsilon e\in\Omega$ стремятся к $y$, так что $\bar\Omega=Q$; тогда $\Omega^*=\{y:\langle y,x\rangle>0\ \forall x\in Q\setminus\{0\}\}=\Omega$
по (iv). Множество $\Omega=\{L>0\}$ выпукло.

(3) Пусть $x=\sum_i\lambda_ic_i$ с $\lambda_i\ge0$. Тогда $0=\sum_i\lambda_i\langle c_i,y\rangle$, причём все слагаемые $\ge0$ по (1).
Если $\lambda_i>0$, то $0=\langle y,c_i\rangle=\langle L(y)c_i,c_i\rangle$ при $L(y)\ge0$, так что $L(y)c_i=0$. Следовательно,
$x\circ y=\sum_i\lambda_i\,c_i\circ y=0$. $\square$

**Теорема F.6 (однородность и квадратичное представление).** Пусть $V$ — евклидова йорданова алгебра, $Q$ и
$\Omega$ — как в теореме F.5.

1. Для каждого $b\in V$: $e^{L(b)}Q=Q$, $e^{L(b)}\Omega=\Omega$ и $e^{L(b)}e=\exp b$.
2. Каждый $y\in\Omega$ равен $\exp(\log y)$. Следовательно, $\Omega$ однороден и является симметрическим конусом.
3. Для $y\in\Omega$ имеем $P(y)=e^{2L(\log y)}$; поэтому $P(y)Q=Q$, $P(y)\Omega=\Omega$ и $P(y)e=y^2$.
4. $P(y)Q\subseteq Q$ для всех $y\in V$.

*Доказательство.* (1) Пусть $0\le s<1/\|L(b)\|$ (любое $s\ge0$, если $L(b)=0$). Симметричный оператор $I-sL(b)$ имеет собственные значения
$\ge1-s\|L(b)\|>0$. Пусть $x\in Q$ и $z=(I-sL(b))^{-1}x=\sum_j\nu_jd_j$ (теорема F.3). Используя $z\circ d_j=\nu_jd_j$,
$\langle d_l,d_j\rangle=\langle d_l,d_j\circ d_j\rangle=\langle d_l\circ d_j,d_j\rangle=\delta_{lj}\|d_j\|^2$ и
$\langle b,d_j\rangle=\langle b,d_j^2\rangle=\langle L(b)d_j,d_j\rangle\le\|L(b)\|\,\|d_j\|^2$, получаем
$$
0\le\langle x,d_j\rangle=\langle z-s\,b\circ z,d_j\rangle=\nu_j\|d_j\|^2-s\langle b,z\circ d_j\rangle
=\nu_j\big(\|d_j\|^2-s\langle b,d_j\rangle\big),
$$
где левое неравенство — это теорема F.5(1)(iv). Скобка $\ge(1-s\|L(b)\|)\|d_j\|^2>0$, поэтому все $\nu_j\ge0$ и
$z\in Q$. Следовательно, для $t\ge0$ и $n>t\|L(b)\|$ имеем $(I-\frac tnL(b))^{-n}Q\subseteq Q$. При $n\to\infty$ этот оператор стремится
к $e^{tL(b)}$ (в ортонормированном собственном базисе $L(b)$: $(1-t\beta/n)^{-n}\to e^{t\beta}$), а $Q$ замкнут, так что
$e^{tL(b)}Q\subseteq Q$. Заменяя $b$ на $-b$, получаем
$e^{-L(b)}Q\subseteq Q$, так что $e^{L(b)}Q=Q$, и гомеоморфизм $e^{L(b)}$ отображает $\operatorname{int}Q$ на себя.
Наконец, $L(b)^re=b^r$ (йордановы степени), поэтому $e^{L(b)}e=\sum_rb^r/r!=\sum_ie^{\beta_i}d_i=\exp b$ для $b=\sum_i\beta_id_i$
по теореме F.3.

(2) Если $y=\sum_i\lambda_ic_i$ с $\lambda_i>0$, то $\log y=\sum_i(\log\lambda_i)c_i$ и $\exp\log y=y$. Поэтому
$y=e^{L(\log y)}e$, где $e^{L(\log y)}\in G(\Omega)$ по (1); самодвойственность — это теорема F.5(2).

(3) По лемме F.4(4) для $c_i$ элемента $y$ оба оператора $P(y)$ и $e^{2L(\log y)}$ действуют на $V_{ij}$ как
$\lambda_i\lambda_j=e^{\log\lambda_i+\log\lambda_j}$. Остаётся воспользоваться (1) и $P(y)e=2y\circ y-y^2=y^2$.

(4) Пусть $y=\sum_i\lambda_ic_i$, $c=\sum_{\lambda_i\ge0}c_i$ (идемпотент), $s=2c-e=\sum_i\varepsilon_ic_i$ с
$\varepsilon_i=\pm1$ и $|y|=\sum_i|\lambda_i|c_i$. По лемме F.4(4) операторы $P(s)$, $P(|y|)$ и $P(y)$ действуют на $V_{ij}$ как
$\varepsilon_i\varepsilon_j$, $|\lambda_i\lambda_j|$ и $\lambda_i\lambda_j$, поэтому $P(y)=P(s)P(|y|)$. Для $\eta>0$
имеем $|y|+\eta e\in\Omega$, так что $P(|y|+\eta e)Q\subseteq Q$ по (3); $P(\cdot)$ — многочлен, а $Q$ замкнут, поэтому
$P(|y|)Q\subseteq Q$ ($\eta\to0$). $P(s)$ — йорданов автоморфизм (лемма F.4(3)), так что $P(s)x^2=(P(s)x)^2\in Q$. $\square$

В разделе 13.6 неравенство $\operatorname{tr}|P(y)x|\le\operatorname{tr}P(y)|x|$ в $V_m$ было выведено из трёх свойств. Они выполняются в
любой евклидовой йордановой алгебре для линейной формы $\tau(x)=\langle x,e\rangle$, которая совпадает с $\operatorname{Tr}$ на
$M_{sa}$ и с $\operatorname{tr}$ на $V_m$ (при указанных выше скалярных произведениях):

**Следствие F.7 (предложение 13.14 в любой евклидовой йордановой алгебре).** Для $x,y$ из евклидовой йордановой алгебры
$\tau(|P(y)x|)\le\tau(P(y)|x|)$, причём равенство имеет место, если $x\in Q$ или $-x\in Q$.

*Доказательство.* $x_\pm\in Q$, $x=x_+-x_-$, $|x|=x_++x_-$ (теорема F.3). $A=P(y)x_+$ и $B=P(y)x_-$ лежат в $Q$ (теорема F.6(4)),
и $z=P(y)x=A-B$. Пусть $z=\sum_i\mu_ic_i$, $c=\sum_{\mu_i>0}c_i$; тогда $|z|=z\circ(2c-e)$ и
$\tau(a\circ b)=\langle a,b\rangle$. Поскольку $\langle Q,Q\rangle\ge0$ (теорема F.5(1)) и $c,e-c\in Q$,
$$
\tau(|z|)=\langle A,c\rangle-\langle A,e-c\rangle-\langle B,c\rangle+\langle B,e-c\rangle\le\tau(A)+\tau(B)=\tau(P(y)|x|).
$$
Если $x\in Q$, то $x_-=0$ и $z=A\in Q$, так что $|z|=z=P(y)|x|$; аналогично при $-x\in Q$. $\square$

Код проверяет теорему F.5 для $\operatorname{Herm}_4(\mathbb C)$ и для конуса Лоренца $Q\subseteq V_4$ (отрицательное
спаривание — это свидетель $c_i$ из (1), $\langle y,c_i\rangle=\lambda_i\|c_i\|^2<0$), теорему F.6 в $V_4$, где в координатах
$(a,u)$ операторы $L(b)$ и $P(y)$ — симметричные матрицы $5\times5$, и следствие F.7 в
$\operatorname{Herm}_5(\mathbb R)$:

```python
import math
import torch
from torch_vn_algebra import SpinFactor

torch.manual_seed(0)
dt = torch.complex128
f64 = torch.float64
close = lambda A, B, tol=1e-10: (A - B).abs().max().item() < tol

G = torch.randn(2, 500, 4, 4, dtype=dt)
X, Y = G[0] @ G[0].mH, G[1] @ G[1].mH                             # squares in Herm_4(C)
assert ((X @ Y).diagonal(dim1=-2, dim2=-1).sum(-1).real > -1e-10).all()   # <Q, Q> >= 0
Z = torch.randn(500, 4, 4, dtype=dt)
Z = Z + Z.mH
w, U = torch.linalg.eigh(Z)
c1 = U[:, :, :1] @ U[:, :, :1].mH                                  # eigenprojection of the smallest eigenvalue
pair = (Z @ c1).diagonal(dim1=-2, dim2=-1).sum(-1).real            # <Z, c1> = lambda_min ||c1||^2
assert close(pair, w[:, 0]) and (pair < 0).float().mean() > 0.9   # witnesses of Z not in Q

V = SpinFactor(4)                                                  # <x, y> = tr(x o y) = V.inner
x, y = V.random(2000, cone=True), V.random(2000, cone=True)
assert (V.inner(x, y) > -1e-12).all()
z = V.random(2000)
lam, idem = V.spectral_decomposition(z)                            # (lambda-, lambda+), (c-, c+); ||c-||^2 = 1
assert close(V.inner(z, idem[:, 0]), lam[:, 0])
assert ((lam[:, 0] < 0) == ~V.in_cone(z, tol=0.0)).all()

I5 = torch.eye(5, dtype=f64)
Lmat = lambda b: V.product(b.expand(5, 5), I5).T                   # columns b o e_k
Pmat = lambda b: V.quadratic_representation(b.expand(5, 5), I5).T
for _ in range(20):
    b, yi, yg = V.random(), V.random(cone=True), V.random()
    E = torch.linalg.matrix_exp(Lmat(b))
    assert close(E @ V.unit(), V.exp(b))                           # Theorem F.6(1)
    assert V.in_cone((E @ x.T).T, tol=1e-9).all()
    assert close(Pmat(yi), torch.linalg.matrix_exp(2 * Lmat(V.log(yi))))      # Theorem F.6(3)
    lam_g, idem_g = V.spectral_decomposition(yg)
    s = (torch.where(lam_g >= 0, 1.0, -1.0)[:, None] * idem_g).sum(0)
    assert close(Pmat(yg), Pmat(s) @ Pmat(V.apply_function(yg, torch.abs)))  # P(y) = P(s) P(|y|)
    assert V.in_cone((Pmat(yg) @ x.T).T, tol=1e-9).all()                     # Theorem F.6(4)

xs = torch.randn(300, 5, 5, dtype=f64)
xs, ys = xs + xs.mT, torch.randn(300, 5, 5, dtype=f64)
ys = ys + ys.mT                                                    # Corollary F.7 in Herm_5(R), tau = Tr
absm = lambda A: (lambda w, U: U @ torch.diag_embed(w.abs()) @ U.mT)(*torch.linalg.eigh(A))
lhs = absm(ys @ xs @ ys).diagonal(dim1=-2, dim2=-1).sum(-1)
rhs = (ys @ absm(xs) @ ys).diagonal(dim1=-2, dim2=-1).sum(-1)
assert (rhs - lhs > -1e-9 * rhs).all()
```

## F.3 Теорема Кёхера–Винберга

В разделе F.2 доказано, что внутренность конуса квадратов евклидовой йордановой алгебры — симметрический конус. Обратное утверждение
восстанавливает йорданово произведение по конусу. Вплоть до теоремы F.12 $\Omega$ — самодвойственный открытый выпуклый конус в
евклидовом пространстве $V$ размерности $N$, а $dy$ — мера Лебега в ортонормированных координатах (у ортогональных замен
координат $|\det|=1$).

**Лемма F.8.** Пусть $\Omega$ самодвойствен.

1. $0\notin\Omega$; $\langle x,y\rangle\ge0$ для $x,y\in\bar\Omega$; для каждого $x\in\Omega$ существует $\delta>0$ с
   $\langle x,y\rangle\ge\delta\|y\|$ для всех $y\in\bar\Omega$.
2. Если $x\in\bar\Omega\setminus\Omega$, то существует $u\in\bar\Omega\setminus\{0\}$ с $\langle x,u\rangle=0$. Если $y\in\Omega$,
   $u\in\bar\Omega$ и $t\ge0$, то $y+tu\in\Omega$.
3. $\Omega$ — внутренность $\bar\Omega$.
4. $G(\Omega)$ — замкнутая подгруппа $GL(V)$, и из $g\in G(\Omega)$ следует $g^T\in G(\Omega)$.

*Доказательство.* (1) Если $0\in\Omega$, то $\Omega$ содержит шар с центром в $0$, а значит, $\Omega=V$; но $V^*=\emptyset$. Для
$x,y\in\Omega$ имеем $\langle x,y\rangle>0$, так как $y\in\Omega^*$; по непрерывности $\langle x,y\rangle\ge0$ на $\bar\Omega$. Если
замкнутый шар $\bar B(x,\delta)$ лежит в $\Omega$ и $y\in\bar\Omega\setminus\{0\}$, то
$\langle x-\delta y/\|y\|,y\rangle\ge0$. (2) Из $x\notin\Omega=\Omega^*$ получаем $u\in\bar\Omega\setminus\{0\}$ с
$\langle u,x\rangle\le0$, а $\langle u,x\rangle\ge0$ по (1). Для $x'\in\bar\Omega\setminus\{0\}$
$\langle x',y+tu\rangle\ge\langle x',y\rangle>0$, так что $y+tu\in\Omega^*=\Omega$. (3) Если $\bar B(y,\varepsilon)\subseteq\bar\Omega$,
то для $x\in\bar\Omega\setminus\{0\}$ имеем $\langle x,y-\varepsilon x/\|x\|\rangle\ge0$, так что $\langle x,y\rangle\ge\varepsilon\|x\|>0$ и
$y\in\Omega^*=\Omega$. (4) Для $g\in G(\Omega)$ отображение $g$ переводит $\bar\Omega\setminus\{0\}$ на себя, поэтому для $y\in\Omega$ и
$x\in\bar\Omega\setminus\{0\}$ имеем $\langle x,g^Ty\rangle=\langle gx,y\rangle>0$: $g^T\Omega\subseteq\Omega$. То же для $g^{-1}$
даёт $g^T\Omega=\Omega$. Если $g_k\in G(\Omega)$ сходятся к $g\in GL(V)$, то множество $g\Omega\subseteq\bar\Omega$ открыто, так что
$g\Omega\subseteq\Omega$ по (3); аналогично $g^{-1}=\lim g_k^{-1}$ переводит $\Omega$ в $\Omega$, так что $g\Omega=\Omega$. $\square$

**Характеристическая функция** конуса $\Omega$ — это $\varphi(x)=\int_{\Omega^*}e^{-\langle x,y\rangle}dy=\int_\Omega e^{-\langle x,y\rangle}dy$.
Положим $f=\log\varphi$, пусть $H(x)$ — симметричный оператор с $\langle H(x)u,v\rangle=D^2f(x)[u,v]$, и
$x^*=-\nabla f(x)$, т. е. $Df(x)[u]=-\langle x^*,u\rangle$.

**Лемма F.9 (характеристическая функция).** Пусть $\Omega$ самодвойствен.

1. $\varphi$ конечна, положительна и принадлежит классу $C^\infty$ на $\Omega$, и её можно дифференцировать под знаком интеграла.
2. $\varphi(tx)=t^{-N}\varphi(x)$ при $t>0$ и $\varphi(gx)=|\det g|^{-1}\varphi(x)$ при $g\in G(\Omega)$.
3. $H(x)>0$ и $x^*\in\Omega$ для $x\in\Omega$.
4. $H(x)x=x^*$ и $D^3f(x)[x,u,v]=-2\langle H(x)u,v\rangle$; для $g\in G(\Omega)$ имеем $(gx)^*=g^{-T}x^*$ и
   $H(gx)=g^{-T}H(x)g^{-1}$.
5. $\varphi(x)\ge c\|x\|^{-N}$ с некоторой константой $c>0$, и $\varphi(x_k)\to\infty$, если $x_k\in\Omega$ сходятся к точке из
   $\bar\Omega\setminus\Omega$.

*Доказательство.* (1) Пусть $x\in\Omega$ и $\delta$ — как в лемме F.8(1). Для $\|x'-x\|<\delta/2$ и $y\in\Omega$
$\langle x',y\rangle\ge\frac\delta2\|y\|$, поэтому производные $k$-го порядка подынтегральной функции по $x'$,
$(-1)^k\prod_{i}\langle u_i,y\rangle e^{-\langle x',y\rangle}$, ограничены функцией $\prod_i\|u_i\|\,\|y\|^ke^{-\delta\|y\|/2}$, которая
интегрируема ($\|y\|\ge N^{-1/2}\sum_i|y_i|$ и теорема Фубини). Теорема о мажорированной сходимости позволяет дифференцировать под знаком
интеграла. (2) Сделаем замену $y=y'/t$, переводящую $\Omega$ на $\Omega$, и замену $y'=g^Ty$, переводящую $\Omega$ на $\Omega$ по лемме
F.8(4), причём $\langle gx,y\rangle=\langle x,y'\rangle$ и $dy=|\det g|^{-1}dy'$.
(3) Пусть $\mu_x$ — вероятностная мера $e^{-\langle x,y\rangle}dy/\varphi(x)$ на $\Omega$. По (1)
$Df(x)[u]=-\mathbb E\langle u,y\rangle$, а $D^2f(x)[u,u]=\mathbb E\langle u,y\rangle^2-(\mathbb E\langle u,y\rangle)^2$ —
дисперсия $\langle u,y\rangle$. Она обращается в нуль, только если $\langle u,y\rangle$ равно константе $\gamma$ $\mu_x$-почти
всюду, т. е. почти всюду по Лебегу на $\Omega$ (плотность $\mu_x$ положительна); множество
$\{y\in\Omega:\langle u,y\rangle\ne\gamma\}$ открыто и имеет меру нуль, а значит, пусто, так что линейная функция $\langle u,\cdot\rangle$
постоянна на шаре, откуда $u=0$. Кроме того, $x^*=\mathbb E\,y$, поэтому $\langle z,x^*\rangle=\mathbb E\langle z,y\rangle>0$
для $z\in\bar\Omega\setminus\{0\}$: $x^*\in\Omega^*=\Omega$.
(4) По (2) $f(tx)=f(x)-N\log t$, так что $Df(x)[x]=-N$ на $\Omega$. Дифференцируя по $x$ в направлении $u$, получаем
$D^2f(x)[x,u]+Df(x)[u]=0$, т. е. $\langle H(x)x,u\rangle=\langle x^*,u\rangle$; ещё раз в направлении $v$:
$D^3f(x)[x,u,v]+2D^2f(x)[u,v]=0$. Дифференцируя $f(gx)=f(x)-\log|\det g|$ один и два раза, получаем
$\langle(gx)^*,gu\rangle=\langle x^*,u\rangle$ и $\langle H(gx)gu,gv\rangle=\langle H(x)u,v\rangle$.
(5) $\langle x,y\rangle\le\|x\|\|y\|\le1$ при $\|y\|\le1/\|x\|$, поэтому $\varphi(x)\ge e^{-1}\operatorname{vol}(\Omega\cap B(0,1/\|x\|))=c\|x\|^{-N}$
с $c=e^{-1}\operatorname{vol}(\Omega\cap B(0,1))>0$. Пусть $x_k\to x_0\in\bar\Omega\setminus\Omega$. По лемме Фату
$\liminf\varphi(x_k)\ge\int_\Omega e^{-\langle x_0,y\rangle}dy$. Возьмём $u$, как в лемме F.8(2), шар $B(y_1,r)\subseteq\Omega$
и непересекающиеся шары $B_l=B(y_1+2lru/\|u\|,r)$, $l\ge0$; они лежат в $\Omega$ по лемме F.8(2), и на $B_l$
$\langle x_0,y\rangle=\langle x_0,y-2lru/\|u\|\rangle\le\|x_0\|(\|y_1\|+r)$. Поэтому интеграл бесконечен. $\square$

**Лемма F.10 (однопараметрические группы).** Пусть $G$ — замкнутая подгруппа $GL(V)$ и
$\mathfrak g=\{Z:e^{tZ}\in G\text{ для всех }t\in\mathbb R\}$.

1. Если $c:[0,s_0)\to G$ удовлетворяет $c(s)=I+sZ+o(s)$ при $s\to0^+$, то $Z\in\mathfrak g$.
2. $\mathfrak g$ — алгебра Ли: линейное подпространство, замкнутое относительно $[Z,W]=ZW-WZ$.

*Доказательство.* (1) Зафиксируем $t>0$ и положим $A_n=c(t/n)$, $B_n=e^{tZ/n}$. Тогда $\|A_n-B_n\|=o(1/n)$ и
$\|A_n\|,\|B_n\|\le1+K/n$ при некотором $K$ и больших $n$, так что
$$
\|A_n^n-e^{tZ}\|=\Big\|\sum_{r=0}^{n-1}A_n^r(A_n-B_n)B_n^{n-1-r}\Big\|\le n\,e^K\|A_n-B_n\|\to0 .
$$
$A_n^n\in G$, а $G$ замкнута в $GL(V)$, поэтому $e^{tZ}\in G$; кроме того, $e^{-tZ}=(e^{tZ})^{-1}\in G$. (2) Для кратных утверждение очевидно.
Применим (1) к $c(s)=e^{sZ}e^{sW}=I+s(Z+W)+O(s^2)$ и к $c(s)=e^{aZ}e^{aW}e^{-aZ}e^{-aW}$ с $a=\sqrt s$: если разложить каждый
множитель до второго порядка, члены порядка $a$ сокращаются, а члены порядка $a^2$ в сумме дают $ZW-WZ$, так что
$c(s)=I+s[Z,W]+O(s^{3/2})$. $\square$

**Предложение F.11 (базовая точка).** Пусть $\Omega$ — симметрический конус, $G=G(\Omega)$, $\mathfrak g$ — как в лемме F.10.

1. $x^{**}=x$ и $H(x^*)=H(x)^{-1}$ для $x\in\Omega$.
2. Существует $e\in\Omega$ с $e^*=e$, и для всякого такого $e$ выполнено $H(e)=I$.
3. Для $x\in\Omega$ и любого $g\in G$ с $ge=x$: $H(x)^{-1}=gg^T\in G$.
4. Если $Z\in\mathfrak g$ и $Ze=0$, то $Z^T=-Z$.

*Доказательство.* (1) Отображение $x\mapsto x^*$ гладкое из $\Omega$ в $\Omega$ (лемма F.9(3)), с производной $-H(x)$. Поскольку
$g^{-T}\in G$ (лемма F.8(4)), лемма F.9(2,4) даёт
$f(gx)+f((gx)^*)=f(x)-\log|\det g|+f(g^{-T}x^*)=f(x)+f(x^*)$, так что $F(x)=f(x)+f(x^*)$ постоянна на орбите
$G x=\Omega$. Следовательно, $0=DF(x)[u]=-\langle x^*,u\rangle+\langle x^{**},H(x)u\rangle$, т. е. $H(x)x^{**}=x^*=H(x)x$
(лемма F.9(4)), и $x^{**}=x$, так как $H(x)$ обратим. Дифференцируя $x^{**}=x$, получаем $H(x^*)H(x)=I$.
(2) $F_2(x)=f(x)+\frac12\|x\|^2$ удовлетворяет $F_2(x)\ge\log c-N\log\|x\|+\frac12\|x\|^2$ (лемма F.9(5)), что стремится к $\infty$ при
$\|x\|\to\infty$. Поэтому $S=\{x\in\Omega:F_2(x)\le F_2(x_1)\}$ ($x_1\in\Omega$) ограничено; оно замкнуто в $V$, потому что предельная
точка в $\bar\Omega\setminus\Omega$ сделала бы $f$ неограниченной на $S$ (лемма F.9(5)). Непрерывная функция $F_2$ достигает минимума
на компакте $S$ в некоторой точке $e$ открытого множества $\Omega$, где $0=\nabla F_2(e)=-e^*+e$. Тогда $H(e)=H(e^*)=H(e)^{-1}$
по (1), так что $H(e)^2=I$, и $H(e)=I$, поскольку $H(e)>0$.
(3) $H(x)=H(ge)=g^{-T}g^{-1}$ по лемме F.9(4) и (2), а $gg^T\in G$ по лемме F.8(4).
(4) $e^{tZ}\in G$ оставляет $e$ на месте, поэтому $I=H(e)=e^{-tZ^T}H(e)e^{-tZ}=e^{-tZ^T}e^{-tZ}$; продифференцируем при $t=0$. $\square$

**Теорема F.12 (Кёхер–Винберг).** Пусть $\Omega$ — непустой открытый выпуклый конус в евклидовом пространстве $V$. Тогда $\Omega$ —
симметрический конус тогда и только тогда, когда на $V$ существует произведение $\circ$, такое, что $(V,\circ,\langle\cdot,\cdot\rangle)$ —
евклидова йорданова алгебра, а $\Omega$ — внутренность её конуса квадратов. Для симметрического конуса в качестве единицы можно взять любой
$e\in\Omega$ с $e^*=e$, а в качестве произведения —
$$
\langle u\circ v,w\rangle=-\tfrac12\,D^3\log\varphi(e)[u,v,w].
$$

*Доказательство.* Часть «если» — это теоремы F.5(2) и F.6(2). Обратно, пусть $\Omega$ симметричен, $e$ — как в предложении F.11(2);
определим $L(u)$ равенством $\langle L(u)v,w\rangle=-\frac12D^3f(e)[u,v,w]$ и положим $u\circ v=L(u)v$.

*Шаг 1: коммутативность, симметричность, единица.* $D^3f(e)$ симметрична по трём аргументам, поэтому $u\circ v=v\circ u$ и
$L(u)^T=L(u)$. По лемме F.9(4) и $H(e)=I$ имеем $\langle e\circ v,w\rangle=\langle H(e)v,w\rangle=\langle v,w\rangle$: $e$ — единица.

*Шаг 2: $L(u)\in\mathfrak g$.* $R(x)=H(x)^{-1}$ лежит в $G$ (предложение F.11(3)), гладко и $R(e)=I$. Его
производная равна $DR(e)[u]=-H(e)^{-1}DH(e)[u]H(e)^{-1}=-DH(e)[u]=2L(u)$, поскольку $\langle DH(e)[u]v,w\rangle=D^3f(e)[u,v,w]$.
Лемма F.10(1) для $c(s)=R(e+su)$ ($s\ge0$ мало, так что $e+su\in\Omega$) и замкнутой подгруппы $G$ (лемма F.8(4))
даёт $2L(u)\in\mathfrak g$, а значит, $L(u)\in\mathfrak g$.

*Шаг 3: единственность.* Если $Y\in\mathfrak g$ симметричен, то $Y=L(Ye)$: симметричный $Y-L(Ye)\in\mathfrak g$ аннулирует $e$
($L(w)e=w\circ e=w$), поэтому он и антисимметричен по предложению F.11(4), а значит, равен $0$.

*Шаг 4: внутренние дифференцирования.* Пусть $D=[L(x),L(y)]$. Он лежит в $\mathfrak g$ (лемма F.10(2)), $D^T=-D$ и
$De=x\circ y-y\circ x=0$. Для $z\in V$ оператор $[D,L(z)]\in\mathfrak g$ симметричен и $[D,L(z)]e=Dz-L(z)De=Dz$, поэтому
$[D,L(z)]=L(Dz)$ по шагу 3; применяя это к $w$, получаем $D(z\circ w)=Dz\circ w+z\circ Dw$.

*Шаг 5: йорданово тождество.* Для $D=[L(y),L(x)]$ имеем $Dx=y\circ x^2-x\circ(x\circ y)$ и
$D(x^2)=y\circ(x\circ x^2)-x\circ(y\circ x^2)$. Шаг 4 даёт $D(x^2)=2x\circ Dx$, т. е.
$y\circ(x\circ x^2)=3x\circ(x^2\circ y)-2x\circ(x\circ(x\circ y))$ для всех $y$:
$$
L(x\circ x^2)=3L(x)L(x^2)-2L(x)^3 .
$$
Все $L$ симметричны, поэтому транспонирование даёт $L(x\circ x^2)=3L(x^2)L(x)-2L(x)^3$, и $[L(x),L(x^2)]=0$. Итак, $V$ —
евклидова йорданова алгебра; пусть $\Omega_J$ — внутренность её конуса квадратов.

*Шаг 6: $\Omega=\Omega_J$.* Для $b\in V$ имеем $e^{L(b)}\in G$ по шагу 2, поэтому $\exp b=e^{L(b)}e\in\Omega$ (теорема F.6(1)). По
теореме F.6(2) $\Omega_J\subseteq\Omega$. Оба конуса самодвойственны относительно $\langle\cdot,\cdot\rangle$ (теорема F.5(2)), поэтому
$\Omega=\Omega^*\subseteq\Omega_J^*=\Omega_J$. $\square$

Доказательство следует структуре книги J. Faraut, A. Korányi, *Analysis on Symmetric Cones* (Oxford 1994), Chapters I–III;
там, где авторы используют структуру группы Ли на $G(\Omega)$, лемма F.10 и предложение F.11(3) непосредственно дают нужные
однопараметрические группы. Теорему доказали M. Koecher (1957) и E. B. Vinberg (1960).

**Пример F.13 (конус полуопределённого программирования).** Пусть $V=\operatorname{Herm}_n(\mathbb R)$, $\langle x,y\rangle=\operatorname{Tr}(xy)$,
$N=n(n+1)/2$ и $\Omega$ — множество положительно определённых матриц, симметрический конус по теоремам F.5 и F.6. Тогда
$\varphi(x)=\varphi(1)\det(x)^{-(n+1)/2}$, $e=\kappa1$ с $\kappa=\sqrt{(n+1)/2}$, и произведение из теоремы F.12 равно
$u\circ_\kappa v=\frac1{2\kappa}(uv+vu)$, так что $x\mapsto x/\kappa$ — изоморфизм на обычную йорданову алгебру.

*Доказательство.* Спектральное разложение $x=\sum_i\lambda_ip_i$ симметричной матрицы по её собственным проекторам — это
разложение из теоремы F.3 (единственность; $p_i\circ p_j=p_ip_j$, поскольку $p_i$ коммутируют), поэтому по теореме F.5(2)
внутренность конуса квадратов — действительно множество положительно определённых матриц. Для $a\in\Omega$ отображение $P(a)x=axa$
(предложение 13.2(1)) лежит в $G(\Omega)$ (теорема F.6(3)). В ортонормированном собственном базисе $a$
матрицы $E_{ii}$ и $(E_{ij}+E_{ji})/\sqrt2$ ($i<j$) образуют ортонормированный базис из собственных векторов $P(a)$ с
собственными значениями $a_ia_j$, поэтому $\det P(a)=\prod_{i\le j}a_ia_j=\det(a)^{n+1}$ (каждое $a_i$ встречается $n+1$ раз). Лемма F.9(2) даёт
$\varphi(a^2)=\varphi(P(a)1)=\det(a)^{-(n+1)}\varphi(1)$, и каждый $x\in\Omega$ равен $a^2$ с $a=\sqrt x$. С помощью
$D\log\det(x)[u]=\operatorname{Tr}(x^{-1}u)$ (из $\det(x+su)=\det x\,\det(1+sx^{-1}u)=\det x\,(1+s\operatorname{Tr}(x^{-1}u)+O(s^2))$)
и $D(x^{-1})[u]=-x^{-1}ux^{-1}$ получаем
$x^*=\kappa^2x^{-1}$; поэтому $e^*=e$ означает $e^2=\kappa^2 1$ при положительно определённой $e$, т. е. $e=\kappa1$; $D^2f(x)[v,w]=\kappa^2\operatorname{Tr}(x^{-1}vx^{-1}w)$ и
$D^3f(e)[u,v,w]=-\kappa^2\kappa^{-3}\operatorname{Tr}(uvw+vuw)=-\frac2\kappa\operatorname{Tr}((u\circ v)w)$. $\square$

Код вычисляет $D^3f$ в точке $e$ автоматическим дифференцированием $f(x)=-\frac{n+1}2\log\det x$ в ортонормированных
координатах $\operatorname{Herm}_3(\mathbb R)$ и проверяет утверждения предложения F.11 и теоремы F.12:

```python
from torch.autograd.functional import jacobian, hessian

n = 3
pairs = [(i, j) for i in range(n) for j in range(i, n)]
Nd = len(pairs)                                                    # N = 6
Bs = torch.zeros(Nd, n, n, dtype=f64)                              # orthonormal basis for Tr(xy)
for k, (i, j) in enumerate(pairs):
    Bs[k, i, j] = Bs[k, j, i] = 1.0 if i == j else 1 / math.sqrt(2)
mat = lambda v: torch.einsum('k,kij->ij', v, Bs)
vec = lambda A: torch.einsum('kij,ij->k', Bs, A)
f = lambda v: -(n + 1) / 2 * torch.logdet(mat(v))                  # log phi up to a constant (Example F.13)
kappa = math.sqrt((n + 1) / 2)
e = vec(kappa * torch.eye(n, dtype=f64))
assert close(-jacobian(f, e), e)                                   # e* = e
assert close(hessian(f, e), torch.eye(Nd, dtype=f64))             # H(e) = I (Proposition F.11(2))
T = jacobian(lambda v: hessian(f, v, create_graph=True), e)        # T[a, b, c] = D^3 f(e)[B_a, B_b, B_c]
prod = lambda u, v: torch.einsum('a,b,abc->c', u, v, -0.5 * T)     # Theorem F.12
u, v = torch.randn(2, Nd, dtype=f64)
assert close(prod(e, v), v)                                        # unit
assert close(prod(u, v), vec(mat(u) @ mat(v) + mat(v) @ mat(u)) / (2 * kappa))
u2 = prod(u, u)
assert close(prod(u, prod(u2, v)), prod(u2, prod(u, v)))           # Jordan identity
```

## F.4 Порядковые изоморфизмы: теорема Кадисона в конечной размерности

Пусть $\mathcal M=\bigoplus_cM_{k_c}(\mathbb C)$ и $\mathcal N=\bigoplus_dM_{k'_d}(\mathbb C)$ — алгебры блочно-диагональных матриц;
$x\ge0$ означает, что $x=x^*$ с неотрицательными собственными значениями (равносильно, каждый блок положительно полуопределён), и
$[0,1]=\{a\in\mathcal M_{sa}:0\le a\le1\}$. **Проектор** — это $p=p^*=p^2$. Линейное отображение $\Phi$ **положительно**, если
из $x\ge0$ следует $\Phi(x)\ge0$; **йорданов $*$-изоморфизм** — это линейная биекция с $\Phi(x^*)=\Phi(x)^*$ и
$\Phi(x\circ y)=\Phi(x)\circ\Phi(y)$ (раздел 13.5). Точка $a$ выпуклого множества $K$ называется **крайней**, если она не является
серединой $\frac12(b+b')$ двух различных точек $b,b'\in K$ (равносильно, не имеет вида $tb+(1-t)b'$ с $b\ne b'$
и $0<t<1$: такая $a$ — середина отрезка между $a\pm\varepsilon(b-b')\in K$ при малом $\varepsilon>0$).

**Лемма F.14.** В $\mathcal M$:

1. Крайние точки выпуклого множества $[0,1]$ — это в точности проекторы.
2. Для проекторов $p,q$: $pq=0$ тогда и только тогда, когда $p+q\le1$, и тогда и только тогда, когда $p\circ q=0$.
3. Для проекторов $p,q$: $p\le q$ тогда и только тогда, когда $pq=p$, и тогда и только тогда, когда $p\circ q=p$. Минимальные ненулевые проекторы — это
   проекторы ранга один (в каком-либо блоке).
4. Проектор $p$ централен тогда и только тогда, когда $pxp+(1-p)x(1-p)=x$ для всех $x\in\mathcal M$. Центральные проекторы — это суммы
   блочных единиц $P_c$, а минимальные ненулевые среди них — это $P_c$.

*Доказательство.* (1) Пусть $p$ — проектор и $p=\frac12(a+b)$ с $a,b\in[0,1]$. Для $\xi\in\ker p$
$0=\langle\xi,a\xi\rangle+\langle\xi,b\xi\rangle$, где оба слагаемых $\ge0$, так что $\|a^{1/2}\xi\|^2=0$ и $a\xi=0$; для $\xi$ из
образа $p$ то же рассуждение с $1-a,1-b$ даёт $a\xi=\xi$. Поэтому $a=p$. Обратно, если $a\in[0,1]$ не проектор,
то у него есть собственное значение $\lambda\in(0,1)$; со спектральным проектором $q_\lambda\in\mathcal M$ (многочленом от $a$) и
$\varepsilon=\min(\lambda,1-\lambda)$ имеем $a\pm\varepsilon q_\lambda\in[0,1]$, и $a$ — их середина.
(2) Если $pq=0$, то $qp=(pq)^*=0$, $p+q$ — проектор и $1-(p+q)=(1-p-q)^2\ge0$; кроме того, $p\circ q=0$. Если
$p+q\le1$, то для любого $\xi$ имеем $\|qp\xi\|^2=\langle p\xi,qp\xi\rangle\le\langle p\xi,(1-p)p\xi\rangle=0$, так что $qp=0=pq$. Если
$p\circ q=0$, то $pq=-qp$, поэтому $pq=p(pq)=-pqp=(qp)p=qp$ и $pq=-pq=0$.
(3) Из $pq=p$ следуют $qp=p$, $p\circ q=p$ и $q-p=(q-p)^2\ge0$. Если $p\circ q=p$, то умножение $pq+qp=2p$ на $p$ слева
и справа даёт $pq+pqp=2p=pqp+qp$, так что $pq=qp=p$. Если $p\le q$, то $p+(1-q)\le1$, и (2) даёт
$p(1-q)=0$. Проектор ранга $\ge2$ строго мажорирует проектор $\xi\xi^*$ на прямую в его образе, лежащую внутри
одного блока ($\xi$ — единичный вектор в образе одного из блоков $p$). Проектор $p$ ранга один лежит в одном блоке, и если
$0\ne q\le p$, то образ $q=pq$ содержится в образе $p$, так что $q=p$.
(4) $x=pxp+px(1-p)+(1-p)xp+(1-p)x(1-p)$, и $px(1-p)+(1-p)xp=0$ тогда и только тогда, когда оба слагаемых равны нулю (умножьте на $p$ и $1-p$),
т. е. тогда и только тогда, когда $px=pxp=xp$. Центральный элемент коммутирует со всеми матричными единицами каждого блока, поэтому на блоке
$c$ он равен скаляру $z_c$, причём для проектора $z_c\in\{0,1\}$. $\square$

**Теорема F.15 (Кадисон, конечная размерность).** Пусть $\Phi:\mathcal M\to\mathcal N$ — линейная биекция с
$\Phi(1)=1$, причём $\Phi$ и $\Phi^{-1}$ положительны. Тогда $\Phi$ — йорданов $*$-изоморфизм. Обратно, всякий
йорданов $*$-изоморфизм унитален, и он сам и обратный к нему положительны.

*Доказательство.* Каждый $x\in\mathcal M_{sa}$ равен $x_+-x_-$ с $x_\pm\ge0$, поэтому $\Phi(x)$ самосопряжён; записывая $x=h+ik$ с
самосопряжёнными $h,k$, получаем $\Phi(x^*)=\Phi(h)-i\Phi(k)=\Phi(x)^*$. То же верно для $\Phi^{-1}$, так что $\Phi$ — вещественно-линейная
биекция $\mathcal M_{sa}\to\mathcal N_{sa}$ с $x\ge0\iff\Phi(x)\ge0$. Поскольку $\Phi(1)=1$, $\Phi$ отображает $[0,1]$ на $[0,1]$.
Линейная биекция между выпуклыми множествами переводит крайние точки на крайние точки (если $\Phi(a)=\frac12(b+b')$, то
$a=\frac12(\Phi^{-1}b+\Phi^{-1}b')$), поэтому по лемме F.14(1) $\Phi$ переводит проекторы на проекторы, а по лемме
F.14(2) $pq=0\iff p+q\le1\iff\Phi(p)+\Phi(q)\le1\iff\Phi(p)\Phi(q)=0$. Для $x=\sum_i\lambda_ip_i\in\mathcal M_{sa}$ (спектральное
разложение, $p_ip_j=0$ при $i\ne j$) $\Phi(x)^2=\sum_i\lambda_i^2\Phi(p_i)=\Phi(x^2)$. Поляризация
$x\circ y=\frac12((x+y)^2-x^2-y^2)$ даёт $\Phi(x\circ y)=\Phi(x)\circ\Phi(y)$ на $\mathcal M_{sa}$, а комплексная
билинейность — на $\mathcal M=\mathcal M_{sa}+i\mathcal M_{sa}$.

Обратно, $\Phi(1)\circ\Phi(y)=\Phi(y)$ для всех $y$, и $\Phi$ сюръективно, поэтому $\Phi(1)$ — единица для $\circ$ и
$\Phi(1)=\Phi(1)\circ1=1$. Каждый $x\ge0$ равен $h^2$ с самосопряжённым $h=\sqrt x$ (теорема 1.6), и
$\Phi(h^2)=\Phi(h)^2\ge0$, поскольку $\Phi(h)$ самосопряжён; $\Phi^{-1}$ — тоже йорданов $*$-изоморфизм. $\square$

Кадисон доказал теорему для всех унитальных C\*-алгебр (R. V. Kadison, *A generalized Schwarz inequality and algebraic
invariants for operator algebras*, Ann. of Math. 56 (1952) 494–503); этот вариант в курсе не используется.

```python
def blocks_rand(ks, psd=False):
    out = []
    for k in ks:
        A = torch.randn(k, k, dtype=dt)
        out.append(A @ A.mH if psd else A + A.mH)
    return out

ks = (2, 2, 3)                                                     # M = M_2 (+) M_2 (+) M_3
Us = [torch.linalg.qr(torch.randn(k, k, dtype=dt))[0] for k in ks]
def Phi(x):                                                        # swap the M_2 blocks, transpose two blocks
    x1, x2, x3 = x
    return [Us[0] @ x2.mT @ Us[0].mH, Us[1] @ x1 @ Us[1].mH, Us[2] @ x3.mT @ Us[2].mH]
def Phi_inv(y):
    y1, y2, y3 = y
    return [Us[1].mH @ y2 @ Us[1], (Us[0].mH @ y1 @ Us[0]).mT, (Us[2].mH @ y3 @ Us[2]).mT]
is_psd = lambda x: all(torch.linalg.eigvalsh(b)[0] > -1e-10 for b in x)
one = [torch.eye(k, dtype=dt) for k in ks]
assert all(close(a, b) for a, b in zip(Phi(one), one))            # unital
for _ in range(50):
    p = blocks_rand(ks, psd=True)
    assert is_psd(Phi(p)) and is_psd(Phi_inv(p))                   # positive in both directions
x, y = blocks_rand(ks), blocks_rand(ks)
assert all(close(a, b) for a, b in zip(Phi_inv(Phi(x)), x))
jordan = lambda a, b: [0.5 * (s @ t + t @ s) for s, t in zip(a, b)]
assert all(close(a, b) for a, b in zip(Phi(jordan(x, y)), jordan(Phi(x), Phi(y))))   # Theorem F.15
assert not close(Phi([s @ t for s, t in zip(x, y)])[0], (Phi(x)[0] @ Phi(y)[0]))      # not multiplicative
```

## F.5 Йордановы $*$-изоморфизмы $\bigoplus_cM_{k_c}$

**Лемма F.16.** Пусть $\Phi:\mathcal M\to\mathcal N$ — йорданов $*$-изоморфизм. Тогда $\Phi(1)=1$,
$\Phi(xyx)=\Phi(x)\Phi(y)\Phi(x)$, и $\Phi$ биективно отображает проекторы $\mathcal M$ на проекторы $\mathcal N$,
сохраняя и отражая $pq=0$, $p\le q$, минимальность, центральность и минимальность среди центральных проекторов.

*Доказательство.* $\Phi(1)=1$, как в теореме F.15. В любой ассоциативной алгебре $xyx=2x\circ(x\circ y)-(x\circ x)\circ y$, что
даёт второе утверждение. Равенство $p=p^*=p\circ p$ сохраняется, как и характеризации из леммы F.14(2,3,4), которые
используют лишь $\circ$, тройное произведение $xyx$, $1$ и кванторы по всей алгебре; $\Phi^{-1}$ — тоже йорданов
$*$-изоморфизм. $\square$

**Теорема F.17 (йордановы $*$-автоморфизмы $M_n$).** Всякий йорданов $*$-изоморфизм $\Phi:M_n(\mathbb C)\to M_n(\mathbb C)$
имеет вид $x\mapsto UxU^*$ или $x\mapsto Ux^TU^*$ с унитарной $U$. При $n\ge2$ отображения первого вида мультипликативны, а второго —
антимультипликативны, так что никакое отображение не имеет обоих видов.

*Доказательство.* При $n=1$ имеем $\Phi=\mathrm{id}$. Пусть $n\ge2$ и $E_{ij}$ — матричные единицы. По лемме F.16 $\Phi(E_{ii})$ —
попарно ортогональные проекторы ранга один, $\Phi(E_{ii})=\xi_i\xi_i^*$ с ортонормированными $\xi_i$; с унитарной $W$,
$We_i=\xi_i$, йорданов $*$-изоморфизм $\Phi_1(x)=W^*\Phi(x)W$ удовлетворяет $\Phi_1(E_{ii})=E_{ii}$.

Пусть $i\ne j$. Пространство $S_{ij}=\operatorname{span}\{E_{ij},E_{ji}\}$ совпадает с $\{x:pxp=x,\ E_{ii}xE_{ii}=E_{jj}xE_{jj}=0\}$ при
$p=E_{ii}+E_{jj}$, поэтому $\Phi_1(S_{ij})\subseteq S_{ij}$ по лемме F.16. Запишем $\Phi_1(E_{ij})=\alpha E_{ij}+\beta E_{ji}$. Из
$E_{ij}\circ E_{ij}=0$: $0=\Phi_1(E_{ij})^2=\alpha\beta(E_{ii}+E_{jj})$, так что ровно одно из чисел $\alpha,\beta$ отлично от нуля
($\Phi_1$ инъективно). Положим $\varepsilon(i,j)=+$, если $\Phi_1(E_{ij})\in\mathbb CE_{ij}$, и $-$, если $\Phi_1(E_{ij})\in\mathbb CE_{ji}$;
$\varepsilon(j,i)=\varepsilon(i,j)$, поскольку $\Phi_1(E_{ji})=\Phi_1(E_{ij})^*$. Для различных $i,j,k$
$E_{ij}\circ E_{jk}=\frac12E_{ik}\ne0$; если $\varepsilon(i,j)=+$ и $\varepsilon(j,k)=-$, то
$\Phi_1(E_{ij})\circ\Phi_1(E_{jk})\in\mathbb C\,E_{ij}\circ E_{kj}=\{0\}$ — противоречие; аналогично для $(-,+)$. Итак,
$\varepsilon(i,j)=\varepsilon(j,k)$ всякий раз, когда $i,j,k$ различны, и $\varepsilon$ постоянно: две пары либо имеют общий
индекс, либо связаны через третью пару. Если $\varepsilon=-$, заменим $\Phi_1$ на $\Phi_1\circ T$ ($T$ — транспонирование,
йорданов $*$-автоморфизм по предложению 13.9); у него $\varepsilon=+$, и оно по-прежнему оставляет $E_{ii}$ на месте.

Теперь $\Phi_1(E_{ij})=\alpha_{ij}E_{ij}$ для всех $i,j$, причём $\alpha_{ii}=1$ и $\alpha_{ji}=\bar\alpha_{ij}$. Из
$E_{ij}\circ E_{ji}=\frac12(E_{ii}+E_{jj})$: $|\alpha_{ij}|^2=1$; из $E_{ij}\circ E_{jk}=\frac12E_{ik}$ ($i,j,k$ различны):
$\alpha_{ij}\alpha_{jk}=\alpha_{ik}$. При $u_i=\alpha_{i1}$ отсюда $\alpha_{ij}=\alpha_{i1}\alpha_{1j}=u_i\bar u_j$ (также при
$i=j$, $i=1$ или $j=1$). Поэтому $\Phi_1(x)=DxD^*$ с $D=\operatorname{diag}(u_i)$, и $\Phi(x)=WDxD^*W^*$ или
$\Phi(x)=WDx^TD^*W^*$. Если бы оба вида совпадали, $T$ было бы мультипликативным, но $(E_{12}E_{21})^T=E_{11}\ne E_{22}=E_{12}^TE_{21}^T$.
$\square$

Тем самым доказан результат N. Jacobson, C. E. Rickart (Trans. Amer. Math. Soc. 69 (1950) 479–502) и I. N.
Herstein (Trans. Amer. Math. Soc. 81 (1956) 331–341), цитированный в разделе 13.5, — в конечномерном случае и для всех
$n$; предложение 13.11 — это случай $n=2$. Йорданов автоморфизм $\operatorname{Herm}_n(\mathbb C)$ продолжается до
йорданова $*$-автоморфизма $M_n$ по формуле $\Phi(x+iy)=\Phi(x)+i\Phi(y)$, так что он тоже имеет один из двух видов.

**Теорема F.18 (прямые суммы).** Пусть $\Phi:\bigoplus_{c=1}^CM_{k_c}\to\bigoplus_{d=1}^{C'}M_{k'_d}$ — йорданов
$*$-изоморфизм. Тогда $C'=C$, и существуют биекция $\sigma$ с $k'_{\sigma(c)}=k_c$, унитарные $U_c\in M_{k_c}$ и
множество блоков $S$, такие, что
$$
\Phi(x)_{\sigma(c)}=U_cx_cU_c^*\ (c\notin S),\qquad \Phi(x)_{\sigma(c)}=U_cx_c^TU_c^*\ (c\in S).
$$
Обратно, всякое такое отображение — йорданов $*$-изоморфизм, и $S$ единственно с точностью до блоков с $k_c=1$.

*Доказательство.* По лемме F.16 и лемме F.14(4) $\Phi(P_c)=P'_{\sigma(c)}$ для некоторой биекции $\sigma$ между минимальными центральными
проекторами, так что $C=C'$. Если $x=P_cxP_c$, то $\Phi(x)=P'_{\sigma(c)}\Phi(x)P'_{\sigma(c)}$, так что $\Phi$ переводит блок $c$ в блок
$\sigma(c)$, а $\Phi^{-1}$ — блок $\sigma(c)$ в блок $c$; следовательно, $k_c^2=k'^2_{\sigma(c)}$, и ограничение есть
йорданов $*$-изоморфизм $M_{k_c}\to M_{k_c}$, по теореме F.17 одного из двух видов. Поскольку $x=\sum_cP_cxP_c$, это
определяет $\Phi$. Обратное верно поблочно (предложение 13.9), а единственность $S$ следует из теоремы F.17. $\square$

**Следствие F.19 (унитальные порядковые изоморфизмы).** Унитальные линейные биекции $\Phi:\mathcal M\to\mathcal N$, для которых
$\Phi$ и $\Phi^{-1}$ положительны, — это в точности отображения из теоремы F.18. Такое $\Phi$ является $2$-положительным тогда и только тогда, когда оно
вполне положительно, тогда и только тогда, когда $S$ не содержит блоков с $k_c\ge2$, и тогда и только тогда, когда $\Phi$ — $*$-изоморфизм.

*Доказательство.* Первое предложение — это теоремы F.15 и F.18. Если в $S$ нет блоков с $k_c\ge2$, то каждое блочное отображение имеет вид
$x\mapsto U_cxU_c^*$ (при $k_c=1$ имеем $x^T=x$), так что $\Phi$ вполне положительно (предложение 6.4) и мультипликативно. Если
$c\in S$ с $k_c\ge2$ и $\Phi$ является $2$-положительным, то таково и блочное отображение $\Phi_{\sigma(c)c}(a)=U_ca^TU_c^*$
(предложение 6.4), а также $a\mapsto U_c^*\Phi_{\sigma(c)c}(a)U_c=a^T$, поскольку из $[A_{ij}]\ge0$ следует
$[U^*A_{ij}U]=(1\otimes U)^*[A_{ij}](1\otimes U)\ge0$; это противоречит предложению 6.1. Такое $\Phi$ и не мультипликативно —
по теореме F.17. $\square$

Итак, унитальный порядковый изоморфизм определяется, с точностью до унитарных матриц, тем, какие одинаковые блоки он переставляет и какие блоки
транспонирует, и по следствию F.19 именно транспонированные блоки мешают полной положительности. Код строит йорданов
автоморфизм $x\mapsto Ux^TU^*$ алгебры $M_3$ и восстанавливает $U$ (с точностью до фазы), следуя шагам доказательства теоремы F.17:

```python
k = 3
E = lambda i, j: torch.zeros(k, k, dtype=dt).index_put_((torch.tensor(i), torch.tensor(j)), torch.tensor(1.0, dtype=dt))
U0 = torch.linalg.qr(torch.randn(k, k, dtype=dt))[0]
Psi = lambda A: U0 @ A.mT @ U0.mH                                  # the unknown Jordan *-automorphism

W = torch.stack([torch.linalg.eigh(Psi(E(i, i)))[1][:, -1] for i in range(k)], dim=1)   # Psi(E_ii) = xi_i xi_i^*
Psi1 = lambda A: W.mH @ Psi(A) @ W
assert all(close(Psi1(E(i, i)), E(i, i)) for i in range(k))
transposed = Psi1(E(0, 1))[1, 0].abs().item() > 0.5               # epsilon = - : Psi1(E_01) is a multiple of E_10
Psi2 = (lambda A: Psi1(A.mT)) if transposed else Psi1
alpha = torch.stack([Psi2(E(i, 0))[i, 0] for i in range(k)])       # u_i = alpha_{i1}
for i in range(k):
    for j in range(k):
        assert close(Psi2(E(i, j)), alpha[i] * alpha[j].conj() * E(i, j))
Ufound = W @ torch.diag(alpha)
X = torch.randn(k, k, dtype=dt)
assert transposed and close(Psi(X), Ufound @ X.mT @ Ufound.mH)     # Theorem F.17, second form
phase = (U0.mH @ Ufound)[0, 0]
assert close(U0.mH @ Ufound, phase * torch.eye(k, dtype=dt)) and abs(abs(phase) - 1) < 1e-12
print("U recovered up to the phase", complex(phase))
```

## Упражнения

1. Покажите, что $\mathbb R^N$ с покоординатным произведением и стандартным скалярным произведением — евклидова йорданова алгебра,
   конус квадратов которой есть $\mathbb R^N_{\ge0}$, и выведите из следствия F.19, что унитальные порядковые изоморфизмы
   $\mathbb C^N=\bigoplus_{c=1}^NM_1$ — это перестановки координат.
2. (**Положительности недостаточно.**) Для $0<t<1$ положим $\Psi_t(x)=(1-t)x+\frac tn\operatorname{Tr}(x)1$ на $M_n$, $n\ge2$. Покажите, что
   $\Psi_t$ — унитальная положительная биекция, что $\Psi_t^{-1}$ не положительно и что $\Psi_t$ не переводит
   проекторы в проекторы.

<details>
<summary>Решения</summary>

1. Произведение коммутативно и ассоциативно, а значит, йорданово, с единицей $(1,\dots,1)$, и $L(x)$ — диагональная матрица
   $\operatorname{diag}(x)$, которая симметрична. Квадраты — это векторы с неотрицательными координатами. По
   следствию F.19 при всех $k_c=1$ имеем $\Phi(x)_{\sigma(c)}=U_cx_cU_c^*=x_c$ (унитарные элементы $M_1$ — фазы), так что
   $\Phi$ переставляет координаты; обратно, перестановки являются унитальными порядковыми изоморфизмами.

2. $\Psi_t(1)=(1-t)1+t1=1$; для $x\ge0$ оба слагаемых $\ge0$. Поскольку $\operatorname{Tr}\Psi_t(x)=\operatorname{Tr}x$,
   имеем $\Psi_t^{-1}(y)=\frac1{1-t}\big(y-\frac tn\operatorname{Tr}(y)1\big)$. Для проектора $p$ ранга один $\Psi_t^{-1}(p)$
   имеет собственное значение $-\frac{t}{n(1-t)}<0$ на $\ker p$, а $\Psi_t(p)$ — собственное значение $\frac tn\notin\{0,1\}$. Это
   согласуется с теоремой F.15: именно положительность обратного отображения заставляет переводить проекторы в проекторы.

```python
nn, t = 3, 0.25
p = E(0, 0)
Psi_t = lambda A: (1 - t) * A + t / nn * torch.trace(A) * torch.eye(nn, dtype=dt)
Psi_t_inv = lambda A: (A - t / nn * torch.trace(A) * torch.eye(nn, dtype=dt)) / (1 - t)
assert close(Psi_t_inv(Psi_t(X)), X)
assert abs(torch.linalg.eigvalsh(Psi_t_inv(p))[0].item() + t / (nn * (1 - t))) < 1e-12
assert not close(Psi_t(p) @ Psi_t(p), Psi_t(p))                   # not a projection
```

</details>

**Назад к:** [главе 13](13_jordan_algebras.md), разделы 13.4–13.6: теорема F.12 доказывает теорему 13.8, теорема F.15
доказывает теорему 13.10 для алгебр $\bigoplus_cM_{k_c}(\mathbb C)$, теорема F.17 распространяет предложение 13.11 на все
$n$, а следствие F.7 распространяет предложение 13.14 на любую евклидову йорданову алгебру.
