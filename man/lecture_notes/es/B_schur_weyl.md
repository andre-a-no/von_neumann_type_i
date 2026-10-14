[English](../B_schur_weyl.md) · **Español** · [Français](../fr/B_schur_weyl.md) · [Deutsch](../de/B_schur_weyl.md) · [中文](../zh/B_schur_weyl.md) · [日本語](../ja/B_schur_weyl.md) · [Русский](../ru/B_schur_weyl.md)

# Apéndice B. Dualidad de Schur–Weyl y cálculo de Weingarten

El grupo unitario $U(d)$ actúa sobre $(\mathbb C^d)^{\otimes n}$ mediante $V\mapsto V^{\otimes n}$, y el grupo simétrico
$S_n$ permutando los factores tensoriales. Cada una de las dos familias de operadores es el conmutante de la otra
(**dualidad de Schur–Weyl**). Este apéndice demuestra la dualidad por completo y deduce de ella la fórmula de
Weingarten para los momentos de una unitaria de Haar de cualquier orden. Por el camino demuestra los hechos que los
capítulos 4, 5, 9 y 11 toman de aquí: el teorema 4.14(b) para todo $L$, la irreducibilidad de las representaciones de
$S_L$ sobre los espacios de multiplicidad del espín total, los números de Catalan de la sección 4.5, el teorema 5.12,
el teorema de Wigner–Eckart para $SU(2)$ usado en la sección 9.3 y la medida de Haar del grupo ortogonal usada en el
capítulo 11.

El camino es elemental: una identidad de polarización identifica el conmutante de las permutaciones con el subespacio
generado por todos los $A^{\otimes n}$, $A\in M_d(\mathbb C)$; un argumento con la transformada de Cayley muestra que las
unitarias $V^{\otimes n}$ generan el mismo espacio; el teorema del doble conmutante (teorema 2.2) convierte esto en la
dualidad. La fórmula de Weingarten se sigue entonces porque el promedio de Haar es la proyección ortogonal sobre el
conmutante.

Las comprobaciones numéricas de este apéndice están en [Tutorial, apéndice B](../../tutorial/es/B_schur_weyl.md).

En todo el apéndice, $d,n\ge1$, $\mathcal H_n=(\mathbb C^d)^{\otimes n}$ con la base ortonormal
$e_{\mathbf i}=e_{i_1}\otimes\dots\otimes e_{i_n}$, donde un multiíndice $\mathbf i=(i_1,\dots,i_n)\in[d]^n$ se lee como una
aplicación $[n]\to[d]$, $m\mapsto i_m$ ($[d]=\{1,\dots,d\}$). $\mathrm{End}\,\mathcal H_n=M_{d^n}(\mathbb C)$ lleva el
producto escalar de Hilbert–Schmidt $\langle X,Y\rangle=\mathrm{Tr}(X^*Y)$. Para $X\in M_d$ y $1\le m\le n$,
$X^{(m)}=1\otimes\dots\otimes X\otimes\dots\otimes1$ ($X$ en la posición $m$) y $d\Pi(X)=\sum_{m=1}^nX^{(m)}$. Escribimos
$$
\mathcal T_n=\mathrm{span}\{A^{\otimes n}:A\in M_d(\mathbb C)\},\qquad
\mathcal U_n=\mathrm{span}\{V^{\otimes n}:V\in U(d)\}.
$$
Los productos de productos de Kronecker se calculan factor a factor, $(A_1\otimes\dots\otimes A_n)(B_1\otimes\dots\otimes B_n)=A_1B_1\otimes\dots\otimes A_nB_n$,
y $(A_1\otimes\dots\otimes A_n)^*=A_1^*\otimes\dots\otimes A_n^*$. Las permutaciones se componen como aplicaciones,
$(\sigma\tau)(m)=\sigma(\tau(m))$, y $\#\text{cycles}(\sigma)$ cuenta los ciclos de $\sigma$, incluidos los puntos fijos.

<a id="s-B-1"></a>

## B.1 Permutaciones de los factores tensoriales

Para $\sigma\in S_n$ sea $\pi(\sigma)e_{\mathbf i}=e_{\mathbf i\circ\sigma^{-1}}$, es decir, el vector de la base en la casilla
$m$ pasa a la casilla $\sigma(m)$. Es el convenio de `decompose.permutation_representation` y de los capítulos 2 y 4.

*En la biblioteca:* los operadores $\pi(\sigma)$ los construye `decompose.permutation_representation` — [tutorial, sección B.3](../../tutorial/es/B_schur_weyl.md#s-B-3).

<a id="r-B-1"></a>

**Lema B.1.** (a) $\pi$ es una representación unitaria de $S_n$, y
$\pi(\sigma)(x_1\otimes\dots\otimes x_n)=x_{\sigma^{-1}(1)}\otimes\dots\otimes x_{\sigma^{-1}(n)}$ para todos los $x_m\in\mathbb C^d$.
(b) $\pi(\sigma)(A_1\otimes\dots\otimes A_n)\pi(\sigma)^*=A_{\sigma^{-1}(1)}\otimes\dots\otimes A_{\sigma^{-1}(n)}$ para todos los
$A_m\in M_d$. En particular, $\pi(\sigma)$ conmuta con todo $A^{\otimes n}$ y todo $d\Pi(X)$, y
$\pi(\sigma)X^{(m)}\pi(\sigma)^*=X^{(\sigma(m))}$.
(c) $\mathrm{Tr}\,\pi(\sigma)=d^{\#\text{cycles}(\sigma)}$ y $\langle\pi(\sigma),\pi(\tau)\rangle=d^{\#\text{cycles}(\sigma^{-1}\tau)}$.
El número de ciclos es invariante por conjugación e inversión, y dos permutaciones son conjugadas si y solo si
tienen el mismo tipo de ciclos.
(d) Los $n!$ operadores $\pi(\sigma)$ son linealmente independientes si y solo si $n\le d$.

*Demostración.* (a) $\mathbf i\mapsto\mathbf i\circ\sigma^{-1}$ es una biyección de $[d]^n$, así que $\pi(\sigma)$ permuta una base
ortonormal y es unitario. $\pi(\sigma)\pi(\tau)e_{\mathbf i}=\pi(\sigma)e_{\mathbf i\circ\tau^{-1}}=e_{\mathbf i\circ\tau^{-1}\circ\sigma^{-1}}=\pi(\sigma\tau)e_{\mathbf i}$.
Desarrollando $x_m=\sum_ax_m(a)e_a$, $x_1\otimes\dots\otimes x_n=\sum_{\mathbf i}\prod_mx_m(i_m)\,e_{\mathbf i}$; aplicando $\pi(\sigma)$ y
sustituyendo $\mathbf j=\mathbf i\circ\sigma^{-1}$ se obtiene $\sum_{\mathbf j}\prod_mx_m(j_{\sigma(m)})\,e_{\mathbf j}=\sum_{\mathbf j}\prod_{p}x_{\sigma^{-1}(p)}(j_p)\,e_{\mathbf j}$,
que es el vector producto afirmado.
(b) Por (a), sobre un vector producto: $\pi(\sigma)^*=\pi(\sigma^{-1})$ pone $x_{\sigma(m)}$ en la casilla $m$, el operador
aplica allí $A_m$, y $\pi(\sigma)$ lleva la casilla $m$ a la casilla $\sigma(m)$; así la casilla $p$ acaba conteniendo
$A_{\sigma^{-1}(p)}x_p$. Los vectores producto generan $\mathcal H_n$. Para $X^{(m)}$, el factor $X$ pasa de la casilla $m$ a
la casilla $\sigma(m)$, de modo que $d\Pi(X)$ queda fijo.
(c) $\pi(\sigma)$ es una matriz de permutación, y $e_{\mathbf i}$ queda fijo si y solo si $i_{\sigma^{-1}(p)}=i_p$ para todo
$p$, es decir, si y solo si $\mathbf i$ es constante sobre los ciclos de $\sigma$: hay $d^{\#\text{cycles}(\sigma)}$ tales
$\mathbf i$. Entonces $\langle\pi(\sigma),\pi(\tau)\rangle=\mathrm{Tr}(\pi(\sigma^{-1})\pi(\tau))=\mathrm{Tr}\,\pi(\sigma^{-1}\tau)$.
Los ciclos de $\kappa\rho\kappa^{-1}$ son las imágenes por $\kappa$ de los ciclos de $\rho$, y $\rho^{-1}$ recorre los mismos
ciclos en sentido inverso. Recíprocamente, si $\rho,\rho'$ tienen el mismo tipo de ciclos, escríbanse sus ciclos uno
debajo del otro, ordenados por longitud, y sea $\kappa$ la aplicación que lleva cada entrada de la primera fila a la
entrada que tiene debajo; entonces $\kappa\rho\kappa^{-1}=\rho'$.
(d) Si $n\le d$, sea $\mathbf u=(1,2,\dots,n)$. Los vectores $\pi(\sigma)e_{\mathbf u}=e_{\mathbf u\circ\sigma^{-1}}$ son vectores de la
base distintos para $\sigma$ distintas, así que $\sum_\sigma c_\sigma\pi(\sigma)=0$ aplicado a $e_{\mathbf u}$ da todos los
$c_\sigma=0$. Si $n>d$, sea $\mathcal A=\sum_\sigma\mathrm{sgn}(\sigma)\pi(\sigma)$. Todo $\mathbf i\in[d]^n$ tiene $i_a=i_b$ para
algún $a\neq b$; para la trasposición $t=(a\,b)$, $\pi(t)e_{\mathbf i}=e_{\mathbf i}$ y
$\mathcal A\pi(t)=\sum_\sigma\mathrm{sgn}(\sigma)\pi(\sigma t)=-\mathcal A$. Por tanto
$\mathcal Ae_{\mathbf i}=\mathcal A\pi(t)e_{\mathbf i}=-\mathcal Ae_{\mathbf i}=0$ para todo $\mathbf i$: $\mathcal A=0$ es una
relación lineal no trivial. $\square$

*Cálculo:* el lema B.1(a) y (d) — [tutorial, sección B.3](../../tutorial/es/B_schur_weyl.md#s-B-3).

<a id="s-B-2"></a>

## B.2 Polarización: el conmutante de las permutaciones

<a id="r-B-2"></a>

**Lema B.2 (polarización).** Para $A_1,\dots,A_n\in M_d(\mathbb C)$, con $|\varepsilon|=\sum_m\varepsilon_m$,
$$
\sum_{\sigma\in S_n}A_{\sigma(1)}\otimes\dots\otimes A_{\sigma(n)}=
\sum_{\varepsilon\in\{0,1\}^n}(-1)^{n-|\varepsilon|}\Big(\sum_{m=1}^n\varepsilon_mA_m\Big)^{\otimes n}.
$$

*Demostración.* Por multilinealidad, $\big(\sum_m\varepsilon_mA_m\big)^{\otimes n}=\sum_f\varepsilon_{f(1)}\cdots\varepsilon_{f(n)}\,A_{f(1)}\otimes\dots\otimes A_{f(n)}$,
donde la suma recorre todas las aplicaciones $f:[n]\to[n]$. Como $\varepsilon_m^r=\varepsilon_m$ para $r\ge1$, el
coeficiente del término $f$ en el miembro derecho es
$$
\sum_{\varepsilon}\prod_{m=1}^n(-1)^{1-\varepsilon_m}\prod_{m\in f([n])}\varepsilon_m
=\prod_{m\notin f([n])}\Big(\sum_{\varepsilon_m=0,1}(-1)^{1-\varepsilon_m}\Big)\prod_{m\in f([n])}\Big(\sum_{\varepsilon_m=0,1}(-1)^{1-\varepsilon_m}\varepsilon_m\Big).
$$
El primer tipo de factor vale $-1+1=0$, el segundo vale $1$. Así, el coeficiente es $1$ si $f$ es sobreyectiva, es
decir, una permutación, y $0$ en otro caso. $\square$

*Cálculo:* el lema B.2 para $n=3$ — [tutorial, sección B.3](../../tutorial/es/B_schur_weyl.md#s-B-3).

El **simetrizador** es $\mathcal S(X)=\frac1{n!}\sum_\sigma\pi(\sigma)X\pi(\sigma)^*$ sobre $\mathrm{End}\,\mathcal H_n$.

<a id="r-B-3"></a>

**Proposición B.3.** $\mathcal S$ es la proyección ortogonal de $\mathrm{End}\,\mathcal H_n$ sobre $\pi(S_n)'$, y
$$
\pi(S_n)'=\mathcal T_n,\qquad\dim\pi(S_n)'=\binom{d^2+n-1}n .
$$

*Demostración.* $\mathcal S$ es el twirl del grupo finito $S_n$ en la representación $\pi$, luego la proyección
ortogonal sobre $\pi(S_n)'$ (teorema 12.5 con $\pi_1=\pi_2=\pi$; véase también la sección 4.2). Por el lema B.1(b),
$\mathcal T_n\subseteq\pi(S_n)'$. Recíprocamente, $\mathrm{End}\,\mathcal H_n$ está generado por los productos
$A_1\otimes\dots\otimes A_n$ de unidades matriciales, y por el lema B.1(b) y la sustitución $\sigma\to\sigma^{-1}$,
$\mathcal S(A_1\otimes\dots\otimes A_n)=\frac1{n!}\sum_\sigma A_{\sigma(1)}\otimes\dots\otimes A_{\sigma(n)}$, que está en $\mathcal T_n$
por el lema B.2. Así, $\pi(S_n)'=\mathcal S(\mathrm{End}\,\mathcal H_n)\subseteq\mathcal T_n$.

*Dimensión.* Sean $P=[d]\times[d]$ y $E_\alpha=E_{\alpha_1}\otimes\dots\otimes E_{\alpha_n}$ para $\alpha\in P^n$, con las unidades
matriciales $E_{(a,b)}=|e_a\rangle\langle e_b|$: una base de $\mathrm{End}\,\mathcal H_n$. Por el lema B.1(b),
$\mathcal S(E_\alpha)=\frac1{n!}\sum_\sigma E_{\alpha\circ\sigma^{-1}}$, que es $|O_\alpha|^{-1}\sum_{\beta\in O_\alpha}E_\beta$ para la órbita
$O_\alpha=\{\alpha\circ\sigma\}$ (cada elemento de la órbita aparece para el mismo número $n!/|O_\alpha|$ de $\sigma$). Las
sumas sobre órbitas distintas involucran conjuntos disjuntos de elementos de la base, así que la imagen de
$\mathcal S$ tiene una base indexada por las órbitas. Dos tuplas están en la misma órbita si y solo si cada elemento de
$P$ aparece el mismo número de veces en ambas: las órbitas son los multiconjuntos de tamaño $n$ formados con los
$N=d^2$ elementos de $P$. Ordenando $P$ como $1,\dots,N$, tal multiconjunto es una sucesión
$1\le c_1\le\dots\le c_n\le N$, y $c_m\mapsto c_m+m-1$ es una biyección sobre las sucesiones estrictamente crecientes
en $[N+n-1]$, es decir, los subconjuntos de $n$ elementos; hay $\binom{N+n-1}n$ de ellos. $\square$

*Cálculo:* $\dim\pi(S_n)'=\binom{d^2+n-1}n$ para $d$ y $n$ pequeños — [tutorial, sección B.3](../../tutorial/es/B_schur_weyl.md#s-B-3).

<a id="s-B-3"></a>

## B.3 Dualidad de Schur–Weyl

<a id="r-B-4"></a>

**Lema B.4 (polinomios que se anulan sobre el grupo unitario).** Sea $p:M_d(\mathbb C)\to\mathbb C$ un polinomio en las
$d^2$ entradas, homogéneo de grado $n\ge1$. Si $p(V)=0$ para todo $V\in U(d)$, entonces $p=0$.

*Demostración.* *Paso 1 (transformada de Cayley).* Para $A\in M_d$ sea $q(A)=p\big((1+iA)\mathrm{adj}(1-iA)\big)$,
donde $\mathrm{adj}$ es la matriz adjunta clásica (adjugada); $q$ es un polinomio en las entradas de $A$. Sea $H$
hermítica. Sus valores propios son reales (corolario 1.4), así que $1-iH$ es invertible,
$\mathrm{adj}(1-iH)=\det(1-iH)(1-iH)^{-1}$ (regla de Cramer), y por homogeneidad $q(H)=\det(1-iH)^n\,p(C)$ con
$C=(1+iH)(1-iH)^{-1}$. Las matrices $1\pm iH$ y sus inversas son funciones de $H$ y conmutan (teorema 1.5), y
$(1+iH)^*=1-iH$, de modo que $C^*C=(1+iH)^{-1}(1-iH)(1+iH)(1-iH)^{-1}=1$. Por tanto $q(H)=0$.

*Paso 2 (de las hermíticas a todas las matrices).* Elíjase una base $B_1,\dots,B_{d^2}$ del espacio vectorial real de
las matrices hermíticas. Es también una base de $M_d(\mathbb C)$ sobre $\mathbb C$, ya que toda $A$ es
$\frac{A+A^*}2+i\frac{A-A^*}{2i}$ y las dimensiones coinciden. $r(t)=q(\sum_kt_kB_k)$ es un polinomio en
$t\in\mathbb C^{d^2}$ que se anula sobre $\mathbb R^{d^2}$. Un polinomio así tiene todos sus coeficientes nulos: en una
variable tiene infinitas raíces; para $N$ variables se escribe $r=\sum_jr_j(t')t_N^j$ — para cada $t'$ real el
polinomio en $t_N$ se anula sobre $\mathbb R$, así que cada $r_j$ se anula sobre $\mathbb R^{N-1}$, y se aplica
inducción. Por tanto $q(A)=r(t)=0$ para toda $A=\sum_kt_kB_k\in M_d$.

*Paso 3 (inversión de la transformada).* Sea $B\in M_d$ con $1+B$ invertible, y $A=-i(1+B)^{-1}(B-1)$. Entonces
$1-iA=1-(1+B)^{-1}(B-1)=2(1+B)^{-1}$ y $1+iA=(1+B)^{-1}2B$; como $B$ conmuta con $1+B$,
$(1+iA)(1-iA)^{-1}=(1+B)^{-1}B(1+B)=B$. Así, $0=q(A)=\det(1-iA)^np(B)$ con $\det(1-iA)\neq0$: $p(B)=0$.

*Paso 4.* Para cualquier $B\neq0$ y $0<t<1/\|B\|$, $-1$ no es valor propio de $tB$ (si $tBv=-v$, entonces
$\|v\|=t\|Bv\|\le t\|B\|\|v\|<\|v\|$ salvo que $v=0$), así que $t^np(B)=p(tB)=0$. Finalmente, $p(0)=0$ por
homogeneidad. $\square$

<a id="r-B-5"></a>

**Lema B.5 (identidades de Newton).** Sean $x_1,\dots,x_n$ elementos de un álgebra que conmutan dos a dos,
$x_S=\prod_{m\in S}x_m$, $e_k=\sum_{|S|=k}x_S$ ($e_0=1$) y $p_j=\sum_mx_m^j$. Entonces
$k\,e_k=\sum_{j=1}^k(-1)^{j-1}e_{k-j}\,p_j$ para $1\le k\le n$.

*Demostración.* Para $1\le j\le k$ sea $a_j=\sum x_Sx_m^j$, sumado sobre los pares $(S,m)$ con $|S|=k-j$ y $m\notin S$,
y $a_{k+1}=0$. En $e_{k-j}p_j=\sum_{|S|=k-j}\sum_mx_Sx_m^j$, los términos con $m\notin S$ dan $a_j$, y los términos con
$m\in S$, escritos con $S'=S\setminus\{m\}$, son $x_{S'}x_m^{j+1}$ y dan $a_{j+1}$. Así la suma alternada es telescópica
y vale $a_1$, y $a_1=\sum_{|T|=k}|T|\,x_T=k\,e_k$ (cada conjunto $T$ de $k$ elementos surge de las $k$ elecciones de
$m\in T$). $\square$

*Cálculo:* el lema B.5 para $x_m=A^{(m)}$ y $n=4$, que recupera $A^{\otimes4}$ a partir de los $d\Pi(A^j)$ — [tutorial, sección B.3](../../tutorial/es/B_schur_weyl.md#s-B-3).

<a id="r-B-6"></a>

**Teorema B.6 (dualidad de Schur–Weyl).** Sean $d,n\ge1$.
(a) $\pi(S_n)'=\mathcal T_n=\mathcal U_n=\mathrm{span}\{V^{\otimes n}:V\in SU(d)\}$. Este espacio es una $*$-álgebra, así
que es el álgebra generada por los $V^{\otimes n}$, $V\in U(d)$.
(b) Para todo conjunto $\mathcal X\subseteq M_d$ cuyo subespacio generado es $M_d$, el álgebra unital generada por
$\{d\Pi(X):X\in\mathcal X\}$ es $\mathcal T_n$.
(c) $\{V^{\otimes n}:V\in U(d)\}'=\{V^{\otimes n}:V\in SU(d)\}'=\{d\Pi(X):X\in M_d\}'=\mathrm{span}\,\pi(S_n)$.

*Demostración.* (a) $\pi(S_n)'=\mathcal T_n$ es la proposición B.3, y $\mathcal U_n\subseteq\mathcal T_n$. Supongamos
$\mathcal U_n\neq\mathcal T_n$. Entonces algún $A_0^{\otimes n}\notin\mathcal U_n$; sea $Y\neq0$ la componente de
$A_0^{\otimes n}$ ortogonal a $\mathcal U_n$. La función $p(A)=\langle Y,A^{\otimes n}\rangle$ es un polinomio homogéneo de
grado $n$ en las entradas de $A$ (cada entrada de $A^{\otimes n}$ es un producto de $n$ entradas de $A$), se anula
sobre $U(d)$, y $p(A_0)=\|Y\|^2\neq0$, en contradicción con el lema B.4.
Para $V\in U(d)$ elíjase $\theta$ con $\det V=e^{i\theta}$; entonces $W=e^{-i\theta/d}V\in SU(d)$ y
$V^{\otimes n}=e^{in\theta/d}W^{\otimes n}$, así que los dos subespacios generados coinciden. $\mathcal U_n$ contiene
$1=1^{\otimes n}$, $V^{\otimes n}W^{\otimes n}=(VW)^{\otimes n}$ y $(V^{\otimes n})^*=(V^*)^{\otimes n}$, así que es una
$*$-álgebra, y toda álgebra que contenga los $V^{\otimes n}$ contiene el subespacio que generan.

(b) Sea $\mathcal A$ dicha álgebra. $d\Pi$ es lineal, así que $\mathcal A$ contiene $d\Pi(X)$ para todo $X\in M_d$. Por el
lema B.1(b), todo $d\Pi(X)$ está en el álgebra unital $\pi(S_n)'=\mathcal T_n$, así que $\mathcal A\subseteq\mathcal T_n$.
Recíprocamente, fijemos $A\in M_d$ y apliquemos el lema B.5 a los $x_m=A^{(m)}$, que conmutan dos a dos: entonces
$e_n=A^{\otimes n}$ y $p_j=d\Pi(A^j)\in\mathcal A$, y por inducción sobre $k$ todo $e_k=\frac1k\sum_j(-1)^{j-1}e_{k-j}p_j$
está en $\mathcal A$. Así, $A^{\otimes n}\in\mathcal A$ y $\mathcal T_n\subseteq\mathcal A$.

(c) El conmutante de un conjunto es igual al conmutante del álgebra que genera. Por (a) y (b), los tres conjuntos
generan $\mathcal T_n=\pi(S_n)'$, cuyo conmutante es $\pi(S_n)''$. El subespacio generado por $\pi(S_n)$ es una
$*$-álgebra ($\pi(\sigma)\pi(\tau)=\pi(\sigma\tau)$, $\pi(\sigma)^*=\pi(\sigma^{-1})$, $\pi(e)=1$), así que
$\pi(S_n)''=\mathrm{span}\,\pi(S_n)$ por el teorema 2.2. $\square$

*Cálculo:* las dimensiones del teorema B.6(a) y (c) para $d$ y $n$ pequeños — [tutorial, sección B.3](../../tutorial/es/B_schur_weyl.md#s-B-3).

<a id="r-B-7"></a>

**Corolario B.7 (teorema 4.14(b)).** Para el espín total $S_a=\sum_is_a^{(i)}$ de $L$ espines $\frac12$, el álgebra
generada por $S_x,S_y,S_z$ es $\{S_a\}''=\pi(S_L)'$, y $\{S_a\}'=\mathrm{span}\,\pi(S_L)$.

*Demostración.* $M_2=\mathrm{span}\{1,s_x,s_y,s_z\}$ (las matrices de Pauli y $1$ forman una base),
$d\Pi(1)=L\cdot1$ y $d\Pi(s_a)=S_a$. Por el teorema B.6(b), el álgebra unital $\mathcal A$ generada por $S_x,S_y,S_z$ es
$\mathcal T_L=\pi(S_L)'$. Como los $S_a$ son hermíticos, $\mathcal A$ es una $*$-álgebra (el adjunto de un producto de
generadores es el producto en orden inverso), así que $\mathcal A=\mathcal A''\supseteq\{S_a\}''$ por el teorema 2.2, y
$\{S_a\}''\supseteq\mathcal A$ porque $\{S_a\}''$ es un álgebra unital que contiene los $S_a$. Finalmente,
$\{S_a\}'=\{S_a\}'''=\pi(S_L)''=\mathrm{span}\,\pi(S_L)$ (sección 2.1 y teorema B.6(c)). $\square$

Por el teorema B.6(a) para $d=2$, la misma álgebra está generada linealmente por las rotaciones $V^{\otimes L}$,
$V\in SU(2)$; esto no necesita la parametrización exponencial de $SU(2)$ mencionada en la sección 4.5.

<a id="s-B-4"></a>

## B.4 La descomposición de $(\mathbb C^d)^{\otimes n}$

<a id="r-B-8"></a>

**Teorema B.8.** Existen $C\ge1$, enteros $k_c,m_c\ge1$ con $\sum_ck_cm_c=d^n$ y una unitaria
$W:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathcal H_n$ tales que
$$
W^*\pi(\sigma)W=\bigoplus_{c=1}^C\rho_c(\sigma)\otimes1_{m_c},\qquad W^*V^{\otimes n}W=\bigoplus_{c=1}^C1_{k_c}\otimes\tau_c(V)
\qquad(\sigma\in S_n,\ V\in U(d)),
$$
donde las $\rho_c$ son representaciones unitarias irreducibles de $S_n$ sobre $\mathbb C^{k_c}$, no equivalentes dos a
dos, y las $\tau_c$ representaciones unitarias irreducibles de $U(d)$ sobre $\mathbb C^{m_c}$, no equivalentes dos a dos.
Además, $W^*\mathrm{span}\,\pi(S_n)W=\bigoplus_cM_{k_c}\otimes1_{m_c}$, $W^*\mathcal U_nW=\bigoplus_c1_{k_c}\otimes M_{m_c}$,
$$
\sum_ck_c^2=\dim\mathrm{span}\,\pi(S_n),\qquad\sum_cm_c^2=\binom{d^2+n-1}n,
$$
y $\sum_ck_c^2=n!$ si $n\le d$.

*Demostración.* La proposición 2.7 para el grupo finito $S_n$ y el teorema 2.6 para
$\mathcal M=\mathrm{span}\,\pi(S_n)=\pi(S_n)''$ dan $W$, las $\rho_c$ con las propiedades enunciadas y
$W^*\pi(S_n)'W=\bigoplus_c1_{k_c}\otimes M_{m_c}$. Por el teorema B.6(a), $\pi(S_n)'=\mathcal U_n$. Así,
$W^*V^{\otimes n}W=\bigoplus_c1\otimes\tau_c(V)$ para ciertos $\tau_c(V)\in M_{m_c}$, y $V\mapsto\tau_c(V)$ es una
representación unitaria porque $V\mapsto V^{\otimes n}$ lo es. Como los $V^{\otimes n}$ generan $\mathcal U_n$, la
aplicación lineal $\sum_la_lV_l^{\otimes n}\mapsto(\sum_la_l\tau_c(V_l))_c$ es sobreyectiva sobre $\bigoplus_cM_{m_c}$; la
irreducibilidad y la no equivalencia de las $\tau_c$ se siguen palabra por palabra como en la demostración de la
proposición 2.7. Las dimensiones son el teorema 2.6, la proposición B.3 y el lema B.1(d). $\square$

*Cálculo:* el teorema B.8 para $S_3$ sobre $(\mathbb C^3)^{\otimes3}$, descompuesto por `decompose.from_group` — [tutorial, sección B.4](../../tutorial/es/B_schur_weyl.md#s-B-4).

Para $d=2$ los sectores se identifican por el espín total:

<a id="r-B-9"></a>

**Corolario B.9 ($L$ espines $\frac12$).** Sea $V$ una unitaria con $V^*S_aV=\bigoplus_SS^{(S)}_a\otimes1_{m_S}$ para el
espín total de $L$ espines $\frac12$ (corolario 4.12, proposición 4.13). Entonces
$V^*\pi(\sigma)V=\bigoplus_S1_{2S+1}\otimes\rho_S(\sigma)$, donde las $\rho_S$ son representaciones unitarias irreducibles
de $S_L$, no equivalentes dos a dos, de dimensión
$$
m_S=\binom Lj-\binom L{j-1},\qquad j=\tfrac L2-S\in\{0,1,\dots,\lfloor L/2\rfloor\}.
$$
En particular, $\dim\mathrm{span}\,\pi(S_L)=\sum_Sm_S^2$ y $\dim\pi(S_L)'=\sum_S(2S+1)^2=\binom{L+3}3$.

*Demostración.* Por el teorema 4.3, $V^*\{S_a\}'V=\bigoplus_S1_{2S+1}\otimes M_{m_S}$, y $\{S_a\}'=\mathrm{span}\,\pi(S_L)$
por el corolario B.7. Así, $V^*\pi(\sigma)V=\bigoplus_S1\otimes\rho_S(\sigma)$ con una representación unitaria $\rho_S$, y la
aplicación $\sum_\sigma a_\sigma\pi(\sigma)\mapsto(\sum_\sigma a_\sigma\rho_S(\sigma))_S$ es sobreyectiva sobre $\bigoplus_SM_{m_S}$;
la irreducibilidad y la no equivalencia se siguen como en la demostración de la proposición 2.7. Los valores de
$m_S$ son la proposición 4.13. Por el teorema 4.3, $\dim\{S_a\}'=\sum_Sm_S^2$ y $\dim\{S_a\}''=\sum_S(2S+1)^2$; por el
corolario B.7, $\{S_a\}'=\mathrm{span}\,\pi(S_L)$ y $\{S_a\}''=\pi(S_L)'$, cuya dimensión es
$\binom{4+L-1}L=\binom{L+3}3$ por la proposición B.3 con $d=2$. $\square$

*Cálculo:* los pares $(m_S,2S+1)$ para $L=5,6$ — [tutorial, sección B.4](../../tutorial/es/B_schur_weyl.md#s-B-4).

<a id="r-B-10"></a>

**Corolario B.10 (números de Catalan).** Para $L\ge1$, $\displaystyle\sum_Sm_S^2=\binom{2L}L\frac1{L+1}$; esto es
$\dim\{S_a\}'=\dim\mathrm{span}\,\pi(S_L)$ sobre $(\mathbb C^2)^{\otimes L}$.

*Demostración.* Escribamos $N=2^L$ y sea $\mathrm{vec}:M_N\to\mathbb C^N\otimes\mathbb C^N$ la biyección lineal
$|x\rangle\langle y|\mapsto x\otimes\bar y$. Entonces $\mathrm{vec}(AX-XA)=(A\otimes1-1\otimes A^T)\mathrm{vec}\,X$, ya que
$A|x\rangle\langle y|B=|Ax\rangle\langle B^*y|\mapsto Ax\otimes B^T\bar y$. Así, $\{S_a\}'$ se aplica sobre el núcleo común de
$S_a\otimes1-1\otimes S_a^T$, $a=x,y,z$, donde $S_a^T=\sum_i(s_a^T)^{(i)}$. De $s_x^T=s_x$, $s_y^T=-s_y$, $s_z^T=s_z$ y
$\sigma_ys_x\sigma_y=-s_x$, $\sigma_ys_y\sigma_y=s_y$, $\sigma_ys_z\sigma_y=-s_z$ obtenemos $\sigma_y(-s_a^T)\sigma_y=s_a$, así que la
unitaria $Y=1\otimes\sigma_y^{\otimes L}$ da $Y(S_a\otimes1-1\otimes S_a^T)Y^*=S_a\otimes1+1\otimes S_a$, el espín total de $2L$
espines $\frac12$ sobre $(\mathbb C^2)^{\otimes2L}=\mathbb C^N\otimes\mathbb C^N$. En la descomposición del corolario 4.12 para
$2L$ espines, un vector anulado por las tres componentes es anulado por $\mathbf S^2$, que vale $S(S+1)\ne0$ sobre los
sumandos de espín $S\ge\frac12$; recíprocamente, las componentes se anulan sobre el sumando de espín $0$. Así, el
núcleo común es ese sumando, de dimensión $m_0(2L)=\binom{2L}L-\binom{2L}{L-1}=\binom{2L}L\big(1-\frac L{L+1}\big)$
(proposición 4.13). La primera igualdad de dimensiones es el teorema 4.3, la segunda el corolario B.7. $\square$

*Cálculo:* el corolario B.10 para $L=5,6$ — [tutorial, sección B.4](../../tutorial/es/B_schur_weyl.md#s-B-4).

<a id="s-B-5"></a>

## B.5 El cálculo de Weingarten

<a id="r-B-11"></a>

**Lema B.11 (el twirl de Haar es la proyección ortogonal sobre el conmutante).** Sea $U$ distribuida según la medida
de Haar en $U(d)$ y $\Phi(X)=\mathbb E[U^{\otimes n}X(U^{\otimes n})^*]$ para $X\in\mathrm{End}\,\mathcal H_n$. Entonces
$\Phi$ es la proyección ortogonal de $\mathrm{End}\,\mathcal H_n$ sobre
$\{V^{\otimes n}:V\in U(d)\}'=\mathrm{span}\,\pi(S_n)$.

*Demostración.* Las entradas de $U^{\otimes n}X(U^{\otimes n})^*$ son polinomios en las entradas de $U$ y $\bar U$, que
están acotadas por $1$, así que $\Phi$ es una aplicación lineal bien definida. Para $V\in U(d)$, $VU$ tiene la ley de
$U$, así que $V^{\otimes n}\Phi(X)(V^{\otimes n})^*=\mathbb E[(VU)^{\otimes n}X((VU)^{\otimes n})^*]=\Phi(X)$: la imagen está
contenida en el conmutante. Si $Y$ está en el conmutante, $U^{\otimes n}Y(U^{\otimes n})^*=Y$ y $\Phi(Y)=Y$. Así, $\Phi$ es
idempotente con imagen el conmutante. Por el lema 5.1, $U^*=U^{-1}$ tiene la ley de $U$, luego
$\langle Y,\Phi(X)\rangle=\mathbb E\mathrm{Tr}\big(((U^*)^{\otimes n}YU^{\otimes n})^*X\big)=\langle\Phi(Y),X\rangle$: $\Phi$ es
autoadjunta. Un idempotente autoadjunto es la proyección ortogonal sobre su imagen
($\langle X-\Phi(X),\Phi(Z)\rangle=\langle\Phi(X-\Phi(X)),Z\rangle=0$). El conmutante es $\mathrm{span}\,\pi(S_n)$ por el
teorema B.6(c). $\square$

*Cálculo:* el lema B.11 por Monte Carlo sobre unitarias de Haar — [tutorial, sección B.5](../../tutorial/es/B_schur_weyl.md#s-B-5).

Para escribir esta proyección en términos de la familia generadora $\pi(\sigma)$, que es linealmente dependiente para
$n>d$, usamos la **pseudoinversa** de una matriz hermítica $G$: $G^+=f(G)$ con $f(\lambda)=1/\lambda$ para
$\lambda\neq0$ y $f(0)=0$ (cálculo funcional, teorema 1.5). Es hermítica, $G^+=G^{-1}$ si $G$ es invertible,
$GG^+=\mathbf 1_{\neq0}(G)$ es la proyección sobre $(\ker G)^\perp$, $G^+v=v/\lambda$ si $Gv=\lambda v$ con $\lambda\neq0$, y
$(R^*GR)^+=R^*G^+R$ para $R$ unitaria (si $G=Q\Lambda Q^*$, entonces $R^*GR=(R^*Q)\Lambda(R^*Q)^*$).

<a id="r-B-12"></a>

**Lema B.12.** Sean $v_1,\dots,v_N$ vectores de un espacio con producto escalar de dimensión finita $\mathcal K$,
$G_{ab}=\langle v_a,v_b\rangle$ su matriz de Gram y $W=G^+$. Entonces la proyección ortogonal sobre
$\mathrm{span}\{v_a\}$ es $x\mapsto\sum_{a,b}v_a\,W_{ab}\,\langle v_b,x\rangle$.

*Demostración.* Sea $T:\mathbb C^N\to\mathcal K$, $Tc=\sum_ac_av_a$. Entonces $T^*x=(\langle v_b,x\rangle)_b$, $G=T^*T$, y la
aplicación en cuestión es $Q=TG^+T^*$. Como $\langle c,Gc\rangle=\|Tc\|^2$, $\ker G=\ker T$, así que $TG^+G=T(1-\Pi)=T$
para la proyección $\Pi$ sobre $\ker T$. Por tanto $QT=T$: $Q$ es la identidad sobre $\mathrm{ran}\,T$. $Q$ se anula
sobre $(\mathrm{ran}\,T)^\perp=\ker T^*$, y $Q^*=Q$. Así, $Q$ es la proyección ortogonal sobre $\mathrm{ran}\,T$. $\square$

<a id="r-B-13"></a>

**Teorema B.13 (fórmula de Weingarten de cualquier orden).** Sean $n,d\ge1$, $U$ distribuida según la medida de Haar
en $U(d)$, $G$ la matriz $n!\times n!$ $G_{\sigma\tau}=d^{\#\text{cycles}(\sigma^{-1}\tau)}$ ($\sigma,\tau\in S_n$) y $W=G^+$.
(a) $W_{\sigma\tau}=\mathrm{Wg}(\sigma^{-1}\tau)$ para la función $\mathrm{Wg}(\rho)=\mathrm{Wg}(\rho,d)=W_{e\rho}$ sobre $S_n$, que
solo depende del tipo de ciclos de $\rho$; en particular, $\mathrm{Wg}(\rho^{-1})=\mathrm{Wg}(\rho)$.
(b) $G$ es invertible si y solo si $n\le d$, y entonces $W=G^{-1}$.
(c) Para $X\in\mathrm{End}\,\mathcal H_n$,
$\mathbb E[U^{\otimes n}X(U^{\otimes n})^*]=\sum_{\sigma,\tau}\mathrm{Wg}(\sigma^{-1}\tau)\mathrm{Tr}(\pi(\tau)^*X)\,\pi(\sigma)$.
(d) Para todos los $\mathbf i,\mathbf j,\mathbf i',\mathbf j'\in[d]^n$,
$$
\mathbb E\big[U_{i_1j_1}\cdots U_{i_nj_n}\overline{U_{i'_1j'_1}\cdots U_{i'_nj'_n}}\big]
=\sum_{\sigma,\tau\in S_n}\prod_{m=1}^n\delta_{i_mi'_{\sigma(m)}}\delta_{j_mj'_{\tau(m)}}\;\mathrm{Wg}(\sigma\tau^{-1}),
$$
es decir, la suma recorre los pares con $\mathbf i=\mathbf i'\circ\sigma$ y $\mathbf j=\mathbf j'\circ\tau$.
(e) Los momentos con números distintos $p\neq q$ de factores $U$ y $\bar U$ se anulan.

Este es el teorema 5.12 (con $k=d$), ahora demostrado.

*Demostración.* (a) Escribamos $g(\rho)=d^{\#\text{cycles}(\rho)}$, de modo que $G_{\sigma\tau}=g(\sigma^{-1}\tau)$, y sean
$\Lambda_\kappa,R_\kappa$ las matrices de permutación sobre $\mathbb C^{S_n}$ con $\Lambda_\kappa e_\sigma=e_{\kappa\sigma}$,
$R_\kappa e_\sigma=e_{\sigma\kappa}$. Entonces $(\Lambda_\kappa^*G\Lambda_\kappa)_{\sigma\tau}=g(\sigma^{-1}\kappa^{-1}\kappa\tau)=G_{\sigma\tau}$
y $(R_\kappa^*GR_\kappa)_{\sigma\tau}=g(\kappa^{-1}\sigma^{-1}\tau\kappa)=G_{\sigma\tau}$ (lema B.1(c)). Por las propiedades de la
pseudoinversa, $W$ tiene las mismas dos invariancias: $W_{\kappa\sigma,\kappa\tau}=W_{\sigma\tau}$ y
$W_{\sigma\kappa,\tau\kappa}=W_{\sigma\tau}$. La primera con $\kappa=\sigma^{-1}$ da $W_{\sigma\tau}=\mathrm{Wg}(\sigma^{-1}\tau)$; la
segunda da $\mathrm{Wg}(\kappa^{-1}\rho\kappa)=\mathrm{Wg}(\rho)$, así que $\mathrm{Wg}$ es constante sobre las clases de
conjugación, que son los tipos de ciclos (lema B.1(c)); $\rho^{-1}$ tiene el tipo de ciclos de $\rho$.
(b) $G$ es la matriz de Gram de los $\pi(\sigma)$ (lema B.1(c)), y con $T$ como en la demostración del lema B.12,
$G=T^*T$ es invertible si y solo si $T$ es inyectiva, es decir, si y solo si los $\pi(\sigma)$ son linealmente
independientes: si y solo si $n\le d$, por el lema B.1(d).
(c) Lema B.11 y lema B.12 con $\mathcal K=\mathrm{End}\,\mathcal H_n$ y $v_\sigma=\pi(\sigma)$.
(d) Las entradas de un producto de Kronecker son $\langle e_{\mathbf i},U^{\otimes n}e_{\mathbf j}\rangle=\prod_mU_{i_mj_m}$, así que
el momento es $\langle e_{\mathbf i},\Phi(X)e_{\mathbf i'}\rangle$ para $X=|e_{\mathbf j}\rangle\langle e_{\mathbf j'}|$. Aquí
$\mathrm{Tr}(\pi(\tau)^*X)=\langle\pi(\tau)e_{\mathbf j'},e_{\mathbf j}\rangle=\delta_{\mathbf j,\mathbf j'\circ\tau^{-1}}$ y
$\langle e_{\mathbf i},\pi(\sigma)e_{\mathbf i'}\rangle=\delta_{\mathbf i,\mathbf i'\circ\sigma^{-1}}$. Por (c), el momento es
$\sum_{\sigma,\tau}\delta_{\mathbf i,\mathbf i'\circ\sigma^{-1}}\delta_{\mathbf j,\mathbf j'\circ\tau^{-1}}\mathrm{Wg}(\sigma^{-1}\tau)$; la
sustitución $\sigma\to\sigma^{-1}$, $\tau\to\tau^{-1}$ da la fórmula.
(e) Para $\theta$ real, $e^{i\theta}1\in U(d)$, así que $e^{i\theta}U$ tiene la ley de $U$, y el momento es igual a
$e^{i(p-q)\theta}$ por sí mismo. Tómese $\theta$ con $e^{i(p-q)\theta}\neq1$. $\square$

Para $n=2\le d$, $G=\begin{pmatrix}d^2&d\\d&d^2\end{pmatrix}$ y $G^{-1}$ da $\mathrm{Wg}(e)=\frac1{d^2-1}$,
$\mathrm{Wg}(\tau)=-\frac1{d(d^2-1)}$: el teorema 5.10. Para $d=1$, todo $G_{\sigma\tau}=1$, $W=\frac1{(n!)^2}G$, y el
momento de la fase $U=U_{11}$ es $\sum_{\sigma,\tau}\frac1{(n!)^2}=1=|U_{11}|^{2n}$, como debe ser.

*Cálculo:* el teorema B.13, con $\mathrm{Wg}$ calculada a partir de la pseudoinversa, incluidos los momentos de orden $3$ para $d=2$, donde $G$ es singular — [tutorial, sección B.5](../../tutorial/es/B_schur_weyl.md#s-B-5).

<a id="r-B-14"></a>

**Corolario B.14.** Para todos los $n,d\ge1$, $\mathbb E|U_{11}|^{2n}=\dfrac{n!}{d(d+1)\cdots(d+n-1)}$.

*Demostración.* Con todos los índices iguales a $1$ en el teorema B.13(d), el momento es $\sum_{\sigma,\tau}\mathrm{Wg}(\sigma\tau^{-1})=\sum_{\sigma,\tau}W_{\sigma^{-1},\tau^{-1}}=\mathbf 1^TW\mathbf 1$
(parte (a)). Toda suma de fila de $G$ es $c_n(d)=\sum_{\rho\in S_n}d^{\#\text{cycles}(\rho)}\ge d^n>0$, así que
$G\mathbf 1=c_n(d)\mathbf 1$, $W\mathbf 1=\mathbf 1/c_n(d)$ y el momento es $n!/c_n(d)$. Finalmente,
$c_n(x)=x(x+1)\cdots(x+n-1)$ (capítulo 5, ejercicio 2, solución). $\square$

Esto extiende el ejercicio 2 del capítulo 5 a $n>d$, donde $G$ es singular.

*Cálculo:* el corolario B.14 para $n=3$, $d=2$ — [tutorial, sección B.5](../../tutorial/es/B_schur_weyl.md#s-B-5).

<a id="s-B-6"></a>

## B.6 Operadores tensoriales y el teorema de Wigner–Eckart para $SU(2)$

La sección 9.3 construye correlaciones de espín a partir del teorema de Wigner–Eckart. Lo demostramos para las
representaciones de espín de la sección 4.5. El espín $k$ es la representación $S^{(k)}_a$ del teorema 4.11 sobre
$\mathbb C^{2k+1}$ con base $|k,q\rangle$, $q=k,\dots,-k$; recuérdese que
$S^{(k)}_\pm|k,q\rangle=\sqrt{k(k+1)-q(q\pm1)}\,|k,q\pm1\rangle$ (para $S_-=S_+^*$ esto se sigue de la fórmula para $S_+$).

**Definición.** Sean $(S_a)$ una representación de espín sobre $\mathbb C^N$ y $k\in\{0,\frac12,1,\dots\}$. Un
**operador tensorial de rango $k$** es una familia $(T_q)_{q=-k,\dots,k}$ de matrices $N\times N$ con
$$
[S_z,T_q]=q\,T_q,\qquad[S_\pm,T_q]=\sqrt{k(k+1)-q(q\pm1)}\;T_{q\pm1}
$$
(el coeficiente se anula para $q=\pm k$, donde $T_{\pm(k+1)}$ no está definido). Un **operador vectorial** es una
terna $(T_x,T_y,T_z)$ con $[S_a,T_b]=i\sum_c\varepsilon_{abc}T_c$, es decir, $[S_x,T_y]=iT_z$, $[S_x,T_z]=-iT_y$ y
cíclicamente, y $[S_a,T_a]=0$.

<a id="r-B-15"></a>

**Lema B.15.** (a) Si $(T_a)$ es un operador vectorial, sus **componentes esféricas**
$T_{\pm1}=\mp\frac1{\sqrt2}(T_x\pm iT_y)$, $T_0=T_z$ forman un operador tensorial de rango $1$, y para dos operadores
vectoriales $\sum_aT_aT'_a=\sum_q(-1)^qT_qT'_{-q}$.
(b) Para $L$ espines $\frac12$ y el espín total $S_a$, cada espín de sitio $(s^{(i)}_a)_a$ y el propio $(S_a)_a$ son
operadores vectoriales.
(c) $(T_q)$ es un operador tensorial de rango $k$ si y solo si la aplicación
$\Psi:\mathbb C^N\otimes\mathbb C^{2k+1}\to\mathbb C^N$, $\Psi(v\otimes|k,q\rangle)=T_qv$, satisface
$\Psi\,(S_a\otimes1+1\otimes S^{(k)}_a)=S_a\Psi$ para $a=x,y,z$.

*Demostración.* (a) De las relaciones, $[S_z,T_x\pm iT_y]=\pm(T_x\pm iT_y)$, $[S_\pm,T_z]=\mp(T_x\pm iT_y)$,
$[S_\pm,T_x\mp iT_y]=\pm2T_z$ y $[S_\pm,T_x\pm iT_y]=0$ (p. ej. $[S_x+iS_y,T_x-iT_y]=-i[S_x,T_y]+i[S_y,T_x]=T_z+T_z$). En
componentes esféricas: $[S_z,T_q]=qT_q$, $[S_\pm,T_0]=\sqrt2\,T_{\pm1}$, $[S_\pm,T_{\mp1}]=\sqrt2\,T_0$, $[S_\pm,T_{\pm1}]=0$,
que son las relaciones de definición para $k=1$ ($\sqrt{2-q(q\pm1)}$ vale $\sqrt2$ para $q=0$ y $q=\mp1$, y $0$ para
$q=\pm1$). Para el producto,
$T_{+1}T'_{-1}+T_{-1}T'_{+1}=-\frac12\big((T_x+iT_y)(T'_x-iT'_y)+(T_x-iT_y)(T'_x+iT'_y)\big)=-(T_xT'_x+T_yT'_y)$.
(b) $[S_a,s^{(i)}_b]=[s^{(i)}_a,s^{(i)}_b]$ (los demás sitios conmutan), y las matrices de espín $\frac12$ satisfacen las
relaciones de espín; para $S_a$ se suma sobre $i$. (c)
$\Psi(S_a\otimes1+1\otimes S^{(k)}_a)(v\otimes|k,q\rangle)=T_qS_av+\sum_{q'}(S^{(k)}_a)_{q'q}T_{q'}v$, así que la relación de
entrelazamiento para $S_a$ dice $[S_a,T_q]=\sum_{q'}(S^{(k)}_a)_{q'q}T_{q'}$ para todo $q$. Para $S_z,S_\pm$ estas son las
relaciones de definición, y las relaciones para $S_x,S_y,S_z$ y para $S_z,S_\pm$ son combinaciones lineales unas de
otras. $\square$

*Cálculo:* el lema B.15(a) y (b) para los espines de sitio de seis espines $\frac12$ — [tutorial, sección B.6](../../tutorial/es/B_schur_weyl.md#s-B-6).

<a id="r-B-16"></a>

**Teorema B.16 (Wigner–Eckart).** Sean $(S_a)$ una representación de espín sobre $\mathbb C^N$, $V$ una unitaria con
$V^*S_aV=\bigoplus_SS^{(S)}_a\otimes1_{m_S}$ (corolario 4.12), $|S,M;b\rangle=V(|S,M\rangle\otimes f_b)$ la base ortonormal
correspondiente ($f_b$ la base canónica de $\mathbb C^{m_S}$), y $(T_q)$ un operador tensorial de rango $k$. Supongamos
que $S,S'$ aparecen ($m_S,m_{S'}\ge1$).
(a) Si $S'\notin\{|S-k|,|S-k|+1,\dots,S+k\}$, entonces $\langle S',M';a|T_q|S,M;b\rangle=0$ para todos los $M',M,q,a,b$.
(b) En otro caso, los operadores de entrelazamiento de espín $S\otimes$ espín $k$ en espín $S'$ forman un espacio de
dimensión uno; para todo tal $C=C^{S'}_{S,k}:\mathbb C^{2S+1}\otimes\mathbb C^{2k+1}\to\mathbb C^{2S'+1}$ no nulo existe una única
matriz $m_{S'}\times m_S$ $R=\langle S'\|T\|S\rangle$, el **elemento de matriz reducido**, con
$$
\langle S',M';a|\,T_q\,|S,M;b\rangle=\langle S',M'|\,C\,(|S,M\rangle\otimes|k,q\rangle)\rangle\;R_{ab}\qquad\text{for all }M',M,q,a,b .
$$
(c) El elemento se anula salvo que $M'=M+q$. En el caso (b) existe un $q$ con
$\langle S',S'|C(|S,S'-q\rangle\otimes|k,q\rangle)\rangle\neq0$, de modo que $R$, y con él todo el bloque de $T$ entre los
multipletes $S$ y $S'$, queda determinado por los elementos de matriz $\langle S',S';a|T_q|S,S'-q;b\rangle$ entre los
vectores de peso máximo de $S'$ y un vector de cada copia de $S$.

*Demostración.* Fijemos $a\le m_{S'}$, $b\le m_S$, y sean $J:\mathbb C^{2S+1}\to\mathbb C^N$, $Jx=V(x\otimes f_b)$ (en el
sumando $S$), y $K:\mathbb C^N\to\mathbb C^{2S'+1}$ el adjunto de $x\mapsto V(x\otimes f_a)$ (en el sumando $S'$). Escribimos
$\alpha\in\{x,y,z\}$ para el índice de componente. Por la forma de $V^*S_\alpha V$, $S_\alpha J=JS^{(S)}_\alpha$ y, tomando
adjuntos en la misma relación para $S'$ (las matrices de espín son hermíticas), $KS_\alpha=S^{(S')}_\alpha K$. Sea
$t=K\Psi(J\otimes1)$ con $\Psi$ del lema B.15(c); sus elementos de matriz son
$\langle S',M'|t(|S,M\rangle\otimes|k,q\rangle)\rangle=\langle S',M';a|T_q|S,M;b\rangle$, y
$$
t\,(S^{(S)}_\alpha\otimes1+1\otimes S^{(k)}_\alpha)=K\Psi(S_\alpha\otimes1+1\otimes S^{(k)}_\alpha)(J\otimes1)=KS_\alpha\Psi(J\otimes1)=S^{(S')}_\alpha\,t,
$$
así que $t$ es un operador de entrelazamiento de $\rho=$ espín $S\otimes$ espín $k$ (la representación de espín
$S^{(S)}_\alpha\otimes1+1\otimes S^{(k)}_\alpha$, sección 4.5) en $\sigma=$ espín $S'$. Tomar adjuntos es una biyección lineal
entre $\mathrm{Hom}(\rho,\sigma)$ y $\mathrm{Hom}(\sigma,\rho)$ (el adjunto de $t\rho(\alpha)=\sigma(\alpha)t$ es
$\rho(\alpha)t^*=t^*\sigma(\alpha)$, pues todas las matrices involucradas son hermíticas). Por el corolario 4.12 y el
ejercicio 2 del capítulo 4 (demostrado en su solución a partir del corolario 4.12 únicamente), $\rho$ es
unitariamente equivalente a $\bigoplus_j(\text{espín }j)$, una copia de cada $j\in\{|S-k|,|S-k|+1,\dots,S+k\}$. Estos espines
son irreducibles y no equivalentes dos a dos (teorema 4.11), así que se trata de una descomposición como en el
teorema 4.3, y su última afirmación da $\dim\mathrm{Hom}(\sigma,\rho)=1$ si $S'$ es uno de estos $j$, y
$\mathrm{Hom}(\sigma,\rho)=0$ en otro caso (el espín $S'$ no es entonces equivalente a ninguno de ellos, por el
teorema 4.11). Esto da (a), y en el caso (b) $t=R_{ab}C$ para un único número $R_{ab}$.
(c) La relación de entrelazamiento para $S_z$ da $M'\langle S',M'|t(|S,M\rangle\otimes|k,q\rangle)\rangle=(M+q)\langle S',M'|t(|S,M\rangle\otimes|k,q\rangle)\rangle$.
La imagen de $C\neq0$ es un subespacio invariante no nulo del espín irreducible $S'$, luego es todo el espacio, así
que $\langle S',S'|C(|S,M\rangle\otimes|k,q\rangle)\rangle\neq0$ para algún vector de la base, y entonces $M=S'-q$ por la
misma relación para el operador de entrelazamiento $C$. Divídase por este número. $\square$

*Cálculo:* el teorema B.16(a) y (b) en una forma independiente de la base, en la descomposición calculada por `decompose.from_generators` — [tutorial, sección B.6](../../tutorial/es/B_schur_weyl.md#s-B-6).

Por el lema B.15(a), el bloque entre los multipletes $S$ y $S'$ de un producto escalar como
$\mathbf s^{(i)}\cdot\mathbf s^{(j)}=\sum_q(-1)^q s^{(i)}_qs^{(j)}_{-q}$ es
$\sum_{S''}\sum_q(-1)^q(\text{bloque }S''\to S'\text{ de }s^{(i)}_q)(\text{bloque }S\to S''\text{ de }s^{(j)}_{-q})$, y cada factor viene
dado por el teorema B.16 mediante un elemento de matriz reducido y una aplicación de Clebsch–Gordan.

*En la biblioteca:* este es el principio en que se basa `su2.spin_correlation`, que el capítulo 9 comprueba frente al método `'direct'` — [tutorial, sección 9.3](../../tutorial/es/09_spin_chains.md#s-9-3).

<a id="s-B-7"></a>

## B.7 Medidas de Haar usadas en el curso

El curso usa las medidas de Haar de los grupos siguientes: $U(k)$ (teorema 5.3), $U(1)$ (el caso $k=1$; capítulo 5,
ejercicio 3), $U(\mathcal M)\cong\prod_cU(k_c)$ (corolario 5.4), $SU(k)$ (sección 5.1, tras el corolario 5.4), el grupo
ortogonal $O(k)$ (sección 5.1, capítulo 11, ejercicio 2, y álgebras reales) y los grupos diagonales $T$ del lema C.2.
La unicidad para todos ellos es el lema 5.1, ya que cada uno es un subgrupo cerrado de algún $U(n)$; la existencia es
una construcción explícita en cada caso (para $T$, en el lema C.2), y la construcción para $O(k)$ es la siguiente. En
ninguna parte del curso se usa un teorema general de existencia de la medida de Haar de un grupo compacto.

<a id="r-B-17"></a>

**Teorema B.17 (medida de Haar del grupo ortogonal).** Sea $Z$ una matriz real $k\times k$ con entradas $N(0,1)$
independientes. Entonces $Z$ es invertible con probabilidad uno, $Z=OR$ con $O$ ortogonal y $R$ real triangular
superior con diagonal positiva, y $O$ está distribuida según la medida de Haar en $O(k)$. Si $Z=O'R'$ es cualquier
descomposición QR real ($O'$ ortogonal, $R'$ real triangular superior), entonces
$O'\mathrm{diag}(\mathrm{sgn}\,R'_{ii})=O$. Para $O$ distribuida según la medida de Haar y
$X\in M_k(\mathbb C)$, $\mathbb E[O]=0$ y $\mathbb E[OXO^T]=\frac{\mathrm{Tr}\,X}k1$; en particular,
$\mathbb E[O_{ij}O_{lm}]=\delta_{il}\delta_{jm}/k$.

*Demostración.* $\det Z$ es un polinomio en las $k^2$ entradas reales, no idénticamente nulo ($\det1=1$), y $Z$ tiene
densidad, así que $\det Z\neq0$ casi seguramente (corolario A.7(b), como en la demostración del teorema 5.3).
Gram–Schmidt aplicado a columnas reales produce $q_j$ reales y coeficientes reales, es decir, una descomposición real
$Z=OR$ con $R_{jj}>0$; por la unicidad del lema 5.2 es la descomposición del lema 5.2, así que $O$ depende
continuamente de $Z$. La densidad $(2\pi)^{-k^2/2}e^{-\|Z\|_2^2/2}$ de $Z$ es invariante por $Z\mapsto VZ$ para
$V\in O(k)$, una aplicación lineal de $\mathbb R^{k^2}$ que conserva $\|Z\|_2^2=\mathrm{Tr}(Z^TZ)$, luego ortogonal con
jacobiano $1$; así, $VZ$ tiene la ley de $Z$ (proposición A.14(a)). Como $VZ=(VO)R$, la unicidad da $O(VZ)=VO(Z)$, de
modo que la ley de $O$ es invariante por la izquierda sobre el grupo compacto $O(k)$, y se aplica el lema 5.1. La
corrección de signos es el final de la demostración del teorema 5.3, con $\mathrm{sgn}\,R'_{ii}=R'_{ii}/|R'_{ii}|$
($R'_{ii}\neq0$ porque $Z$ es invertible). Como $-1\in O(k)$, $-O$ tiene la ley de $O$, así que
$\mathbb E[O]=-\mathbb E[O]=0$ (las entradas están acotadas por $1$). Finalmente, $T=\mathbb E[OXO^T]$ satisface
$VTV^T=T$ para $V\in O(k)$, ya que $VO$ tiene la ley de $O$. Con $V$ la matriz diagonal con $-1$ en la posición $i$ y
$1$ en las demás, $T_{ij}=-T_{ij}$ para $j\neq i$; con matrices de permutación, todos los $T_{ii}$ son iguales (como en
el lema 5.6). Así, $T=\lambda1$, y $k\lambda=\mathrm{Tr}\,T=\mathrm{Tr}\,X$. Para $X=|e_j\rangle\langle e_m|$,
$(OXO^T)_{il}=O_{ij}O_{lm}$. $\square$

*Cálculo:* los segundos momentos del teorema B.17 por Monte Carlo, y por qué hace falta la corrección de signos — [tutorial, sección B.7](../../tutorial/es/B_schur_weyl.md#s-B-7).

**Volver a:** [capítulo 4](04_symmetries.md) (teorema 4.14), [capítulo 5](05_random_matrices.md) (teorema 5.12),
[capítulo 9](09_spin_chains.md) (correlaciones).
