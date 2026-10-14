[English](../F_cones_and_order.md) · **Español** · [Français](../fr/F_cones_and_order.md) · [Deutsch](../de/F_cones_and_order.md) · [中文](../zh/F_cones_and_order.md) · [日本語](../ja/F_cones_and_order.md) · [Русский](../ru/F_cones_and_order.md)

# Apéndice F. Conos simétricos, isomorfismos de orden y aplicaciones de Jordan

El capítulo 13 relacionó el orden de los observables con su producto de Jordan mediante tres teoremas: el teorema de
Koecher–Vinberg (teorema 13.8), el teorema de Kadison (teorema 13.10) y la descripción de los automorfismos de Jordan
de $M_n(\mathbb C)$ (la proposición 13.11 la demuestra para $n=2$). Este apéndice demuestra el primero con total
generalidad y los otros dos para las álgebras $\bigoplus_cM_{k_c}(\mathbb C)$ del curso, para todos los tamaños de
bloque. Por el camino demuestra los hechos sobre álgebras de Jordan euclídeas que necesita la sección 13.6: un teorema
espectral, la autodualidad del cono de cuadrados y $P(y)Q\subseteq Q$.

**Aprenderá** por qué toda álgebra de Jordan euclídea es asociativa en potencias y tiene un teorema espectral; por qué
su cono de cuadrados es autodual y homogéneo; cómo se recupera un producto de Jordan a partir de un cono autodual y
homogéneo solamente, mediante la tercera derivada del logaritmo de su función característica; por qué una biyección
lineal unital que conserva la positividad en ambos sentidos conserva los cuadrados; y por qué toda aplicación de este
tipo de $\bigoplus_cM_{k_c}$ permuta bloques iguales y actúa sobre cada bloque mediante $UxU^*$ o mediante $Ux^TU^*$.

Además de álgebra lineal, cálculo y resultados de los capítulos 1, 6 y 13 citados por su número, la sección F.3 usa la
integración de Lebesgue (el lema de Fatou, la convergencia dominada y la derivación bajo el signo integral,
Fubini–Tonelli, la fórmula del cambio lineal de variables); no se supone nada más. El teorema espectral del apéndice E,
con marcos de Jordan, refina el teorema F.3 de más abajo; nada de lo que sigue depende de él. Un
**espacio euclídeo** es un espacio vectorial real de dimensión finita $V\ne\{0\}$ con un producto escalar
$\langle\cdot,\cdot\rangle$; $\|x\|$ es la norma (las barras simples $|x|$ de un elemento de $V$ denotarán el valor
absoluto del teorema F.3, nunca la norma), $\|A\|$ la norma de operador de una aplicación lineal $A$ de $V$, $A^T$ su
adjunta ($\langle Ax,y\rangle=\langle x,A^Ty\rangle$), y $A\ge0$ ($A>0$) significa que $A=A^T$ es semidefinida positiva
(definida positiva). $e^A=\exp A=\sum_{k\ge0}A^k/k!$.

## F.1 Álgebras de Jordan euclídeas

Como en la sección 13.1, un **álgebra de Jordan** es un espacio vectorial real con un producto bilineal y conmutativo
$\circ$ tal que $x\circ(x^2\circ y)=x^2\circ(x\circ y)$, es decir, $[L(x),L(x^2)]=0$, donde $L(x)y=x\circ y$ y
$x^2=x\circ x$; y $P(x)=2L(x)^2-L(x^2)$.

**Definición.** Un **álgebra de Jordan euclídea** es un álgebra de Jordan $V$ de dimensión finita con una unidad $e$
($e\circ x=x$) y un producto escalar para el cual todo $L(x)$ es simétrico: $\langle x\circ y,z\rangle=\langle y,x\circ z\rangle$.
Su **cono de cuadrados** es $Q=\{x^2:x\in V\}$.

Ejemplos: $M_{sa}$ (la parte autoadjunta de $\bigoplus_cM_{k_c}(\mathbb C)$) y $\operatorname{Herm}_n(\mathbb R)$
con $\langle x,y\rangle=\operatorname{Tr}(xy)$, que es real y cumple $\operatorname{Tr}(x^2)=\sum_{ij}|x_{ij}|^2$ (lema 13.1;
$\operatorname{Tr}((x\circ y)z)=\frac12\operatorname{Tr}(xyz+yxz)=\operatorname{Re}\operatorname{Tr}(xyz)$, ya que
$\operatorname{Tr}(yxz)=\overline{\operatorname{Tr}((yxz)^*)}=\overline{\operatorname{Tr}(zxy)}=\overline{\operatorname{Tr}(xyz)}$,
es invariante bajo permutaciones cíclicas y, por la misma razón, bajo $x\leftrightarrow z$, luego simétrico en
$x,y,z$); los factores de espín $V_m$ (álgebras de Jordan por el lema 13.3) con el producto escalar
$\langle x,y\rangle=\operatorname{tr}(x\circ y)=2(ab+u\cdot v)$ para $x=(a,u)$, $y=(b,v)$, ya que
$\operatorname{tr}((x\circ y)\circ z)=2(abc+c\,u\cdot v+a\,v\cdot w+b\,u\cdot w)$ para $z=(c,w)$ es simétrico en
$x,y,z$; y $\mathbb R^N$ con el producto coordenada a coordenada (ejercicio 1).

**Lema F.1.** En toda álgebra de Jordan, para todos $a,b,d,x,y$:

1. $[L(a),L(b\circ d)]+[L(b),L(d\circ a)]+[L(d),L(a\circ b)]=0$;
2. $L(x^2\circ y)=L(x^2)L(y)-2L(x)L(y)L(x)+2L(x\circ y)L(x)$.

*Demostración.* (1) Para $s,t,r\in\mathbb R$ y $x=sa+tb+rd$, $[L(x),L(x^2)]=0$ es una identidad polinómica en $(s,t,r)$,
así que cada coeficiente se anula; el coeficiente de $str$ es el doble del lado izquierdo de (1). (2) Pongamos $a=b=x$
en (1): $2[L(x),L(x\circ d)]+[L(d),L(x^2)]=0$. Aplicado a $y$ esto se lee
$2x\circ((x\circ d)\circ y)-2(x\circ d)\circ(x\circ y)+d\circ(x^2\circ y)-x^2\circ(d\circ y)=0$, es decir,
$\big(2L(x)L(y)L(x)-2L(x\circ y)L(x)+L(x^2\circ y)-L(x^2)L(y)\big)d=0$ para todo $d$. $\square$

**Proposición F.2 (asociatividad de potencias).** Sean $V$ un álgebra de Jordan con unidad $e$, $x\in V$, $x^0=e$ y
$x^{k+1}=x\circ x^k$. Entonces $x^p\circ x^q=x^{p+q}$ para todos $p,q\ge0$, y todo $L(x^p)$ es un polinomio en los
operadores conmutantes $L(x)$ y $L(x^2)$. Por tanto, $p(x)\circ q(x)=(pq)(x)$ para polinomios reales $p,q$, y los
operadores $L(p(x))$ conmutan entre sí.

*Demostración.* $L(x)$ y $L(x^2)$ conmutan por la identidad de Jordan. Demostramos por inducción sobre $N\ge2$: (A$_N$)
$x^i\circ x^j=x^{i+j}$ para $i+j\le N$, y (B$_N$) $L(x^i)$ es un polinomio en $L(x),L(x^2)$ para $i\le N$. Ambas son
claras para $N=2$. Supongámoslas para $N$ y sea $i+j=N+1$ con $i,j\ge1$ (para $i=0$ o $j=0$ no hay nada que demostrar).
Como $i\le N$, $L(x^i)$ conmuta con $L(x)$, así que por (A$_N$)
$$
x^i\circ x^j=L(x^i)L(x)x^{j-1}=L(x)L(x^i)x^{j-1}=x\circ x^N=x^{N+1}.
$$
Esto es (A$_{N+1}$). El lema F.1(2) con $y=x^{N-1}$ y $x^2\circ x^{N-1}=x^{N+1}$ da
$L(x^{N+1})=L(x^2)L(x^{N-1})-2L(x)L(x^{N-1})L(x)+2L(x^N)L(x)$, que es (B$_{N+1}$). La última frase se sigue por
bilinealidad. $\square$

**Teorema F.3 (teorema espectral).** Sean $V$ un álgebra de Jordan euclídea y $x\in V$. Existen únicos
$k\ge1$, números reales $\lambda_1>\dots>\lambda_k$ y $c_1,\dots,c_k\in V$ no nulos con
$$
c_i\circ c_j=\delta_{ij}c_i,\qquad \textstyle\sum_ic_i=e,\qquad x=\sum_i\lambda_ic_i .
$$
Los $c_i$ son polinomios en $x$, y $p(x)=\sum_ip(\lambda_i)c_i$ para todo polinomio real $p$.

Llamamos a $(c_i)$ un **sistema completo de idempotentes ortogonales** y a los $\lambda_i$ los **valores propios** de
$x$, y ponemos $f(x)=\sum_if(\lambda_i)c_i$ para $f:\{\lambda_1,\dots,\lambda_k\}\to\mathbb R$; por ejemplo $|x|$,
$x_\pm=\max(\pm x,0)$, $\sqrt x$ (si todos los $\lambda_i\ge0$), $\exp x$ y $\log x$ (si todos los $\lambda_i>0$).

*Demostración.* Existencia. $\mathbb R[x]=\{p(x)\}$ es un subespacio que contiene a $e$ y que, por la proposición F.2,
cumple $p(x)\circ q(x)=(pq)(x)$. Así, $M=L(x)|_{\mathbb R[x]}$ lleva $\mathbb R[x]$ en sí mismo, es simétrico y cumple
$p(M)q(x)=p(x)\circ q(x)$. Un operador simétrico es diagonalizable con valores propios reales; sean $\lambda_1>\dots>\lambda_k$
los de $M$ y $m(t)=\prod_i(t-\lambda_i)$ su polinomio mínimo. Entonces $m(x)=m(M)e=0$. Con los polinomios de Lagrange
$\ell_i(t)=\prod_{j\ne i}\frac{t-\lambda_j}{\lambda_i-\lambda_j}$, los polinomios
$\ell_i\ell_j-\delta_{ij}\ell_i$, $\sum_i\ell_i-1$ y $\sum_i\lambda_i\ell_i-t$ se anulan en todo $\lambda_l$, así que son
múltiplos de $m$; evaluando en $x$ se ve que los $c_i=\ell_i(x)$ tienen las tres propiedades. Si $c_i=0$, entonces
$\ell_i(M)q(x)=\ell_i(x)\circ q(x)=0$ para todo $q$, luego $\ell_i(M)=0$ aunque $\deg\ell_i<\deg m$; por tanto $c_i\ne0$.
Finalmente $p(x)=\sum_ip(\lambda_i)c_i$: es cierto para $p=1$ y $p=t$, y $x\circ\sum_i\lambda_i^rc_i=\sum_i\lambda_i^{r+1}c_i$.

Unicidad. Sea $x=\sum_j\mu_jd_j$ con $\mu_j$ distintos y un sistema completo de idempotentes ortogonales no nulos $d_j$.
Como antes, $p(x)=\sum_jp(\mu_j)d_j$. Los $d_j$ son linealmente independientes (multiplíquese $\sum_ja_jd_j=0$ por $d_i$),
así que $p(x)=0$ si y solo si $p(\mu_j)=0$ para todo $j$. Así, $m(x)=0$ da $\{\mu_j\}\subseteq\{\lambda_i\}$, y
$\tilde m(t)=\prod_j(t-\mu_j)$ cumple $\tilde m(M)q(x)=\tilde m(x)\circ q(x)=0$, luego $m$ divide a $\tilde m$ y $\{\lambda_i\}\subseteq\{\mu_j\}$.
Entonces $d_j=\ell_j(x)=c_j$. $\square$

**Lema F.4 (descomposición de Peirce).** Sean $V$ un álgebra de Jordan euclídea y $c$ un idempotente ($c^2=c$).

1. $L(c)(2L(c)-I)(L(c)-I)=0$. Por tanto, $V$ es la suma ortogonal de los subespacios propios $V(c,\nu)=\{z:c\circ z=\nu z\}$,
   $\nu\in\{0,\frac12,1\}$, y $L(c)\ge0$.
2. $V(c,1)\circ V(c,0)=\{0\}$; $V(c,\nu)\circ V(c,\nu)\subseteq V(c,\nu)$ y $V(c,\nu)\circ V(c,\frac12)\subseteq V(c,\frac12)$
   para $\nu\in\{0,1\}$; $V(c,\frac12)\circ V(c,\frac12)\subseteq V(c,0)+V(c,1)$.
3. $s=2c-e$ cumple $s^2=e$, y $P(s)$ es un automorfismo de Jordan de $V$: vale $+1$ sobre $V(c,0)+V(c,1)$ y $-1$ sobre
   $V(c,\frac12)$, y $P(s)(z\circ w)=P(s)z\circ P(s)w$.
4. Sea $x=\sum_{i=1}^k\lambda_ic_i$ como en el teorema F.3. Los $L(c_i)$ conmutan, y $V$ es la suma ortogonal de los
   subespacios $V_{ii}=\{z:c_i\circ z=z\}$ y $V_{ij}=\{z:c_i\circ z=c_j\circ z=\frac12z\}$ ($i<j$). Para
   $a=\sum_la_lc_l$ ($a_l\in\mathbb R$), $L(a)$ actúa sobre $V_{ij}$ ($i\le j$) como el escalar $\frac12(a_i+a_j)$, y $P(a)$
   como el escalar $a_ia_j$.

*Demostración.* (1) El lema F.1(2) con $x=y=c$ da $L(c)=L(c)^2-2L(c)^3+2L(c)^2$, es decir, $2L(c)^3-3L(c)^2+L(c)=0$, que
es la factorización enunciada. Así, los valores propios del operador simétrico $L(c)$ están en $\{0,\frac12,1\}$.

(2) El lema F.1(1) con $a=b=c$, $d=z\in V(c,\nu)$ da $2\nu[L(c),L(z)]+[L(z),L(c)]=0$, luego $L(z)$ conmuta con $L(c)$
si $\nu\in\{0,1\}$ y entonces conserva todo $V(c,\mu)$. Esto da la segunda y la tercera inclusión, y
$V(c,1)\circ V(c,0)\subseteq V(c,0)\cap V(c,1)=\{0\}$. Para $z,w\in V(c,\frac12)$, apliquemos el lema F.1(1) con $(a,b,d)=(z,w,c)$
a $c$: como $L(w\circ c)=\frac12L(w)$ y $L(c\circ z)=\frac12L(z)$, los dos primeros términos se cancelan, y
$[L(c),L(u)]c=L(c)^2u-L(c)u=0$ para $u=z\circ w$. Escribiendo $u=u_1+u_{1/2}+u_0$ según (1), esto es $-\frac14u_{1/2}=0$.

(3) $s^2=4c-4c+e=e$, $L(s)=2L(c)-I$ vale $1,0,-1$ sobre $V(c,1),V(c,\frac12),V(c,0)$, y $P(s)=2L(s)^2-I$. Escribamos
$z=z'+z''$, $w=w'+w''$ con $z',w'\in V(c,0)+V(c,1)$ y $z'',w''\in V(c,\frac12)$. Por (2), $z'\circ w'+z''\circ w''$ está en
$V(c,0)+V(c,1)$ y $z'\circ w''+z''\circ w'$ en $V(c,\frac12)$, así que
$P(s)(z\circ w)=(z'-z'')\circ(w'-w'')=P(s)z\circ P(s)w$.

(4) Los $c_l$ son polinomios en $x$, así que los $L(c_l)$ conmutan (proposición F.2); al ser simétricos, tienen una base
ortonormal común de vectores propios. Sobre un vector propio común, $L(c_l)$ tiene un valor propio $\varepsilon_l\in\{0,\frac12,1\}$
por (1), y $\sum_l\varepsilon_l=1$ porque $\sum_lL(c_l)=L(e)=I$. Así, o bien un $\varepsilon_i=1$ y los demás se anulan, o
bien $\varepsilon_i=\varepsilon_j=\frac12$ para un par $i<j$ y los demás se anulan. Por tanto, $V$ es la suma ortogonal de
los subespacios propios comunes $W_{ii}$ (donde $L(c_i)=1$, $L(c_l)=0$ para $l\ne i$) y $W_{ij}$ (donde $L(c_i)=L(c_j)=\frac12$,
$L(c_l)=0$ en los demás casos), y $W_{ij}=V_{ij}$: un vector con $c_i\circ z=z$ (resp. $c_i\circ z=c_j\circ z=\frac12z$) no
tiene componente en los otros subespacios propios comunes, sobre los cuales $L(c_i)$ (resp. el par $(L(c_i),L(c_j))$)
toma un valor distinto. Como
$a^2=\sum_la_l^2c_l$, sobre $V_{ij}$ ($i<j$) se tiene $L(a)=\frac12(a_i+a_j)$, $L(a^2)=\frac12(a_i^2+a_j^2)$ y
$P(a)=\frac12(a_i+a_j)^2-\frac12(a_i^2+a_j^2)=a_ia_j$; sobre $V_{ii}$, $L(a)=a_i$ y $P(a)=2a_i^2-a_i^2$. $\square$

## F.2 El cono de cuadrados es simétrico

Sea $\Omega$ un **cono convexo** abierto no vacío en un espacio euclídeo $V$ ($tx\in\Omega$ para $x\in\Omega$, $t>0$). Su
**cono dual** (abierto) es $\Omega^*=\{y:\langle x,y\rangle>0\text{ para todo }x\in\bar\Omega\setminus\{0\}\}$; $\Omega$ es
**autodual** si $\Omega^*=\Omega$. $G(\Omega)=\{g\in GL(V):g\Omega=\Omega\}$, y $\Omega$ es **homogéneo** si
$G(\Omega)$ actúa transitivamente sobre $\Omega$. Un **cono simétrico** es uno autodual y homogéneo. Si $\Omega_1\subseteq\Omega_2$,
entonces $\Omega_2^*\subseteq\Omega_1^*$.

**Teorema F.5 (autodualidad).** Sean $V$ un álgebra de Jordan euclídea e $y=\sum_i\lambda_ic_i$ (teorema F.3).

1. Son equivalentes: (i) $y\in Q$; (ii) todos los $\lambda_i\ge0$; (iii) $L(y)\ge0$; (iv) $\langle y,x\rangle\ge0$
   para todo $x\in Q$. En particular, $Q$ es un cono convexo cerrado.
2. Son equivalentes: (i) $y\in\Omega:=\operatorname{int}Q$; (ii) todos los $\lambda_i>0$; (iii) $L(y)>0$; (iv)
   $\langle y,x\rangle>0$ para todo $x\in Q\setminus\{0\}$. Además, $e\in\Omega$, $\bar\Omega=Q$ y $\Omega^*=\Omega$:
   $\Omega$ es un cono convexo abierto autodual.
3. Si $x,y\in Q$ y $\langle x,y\rangle=0$, entonces $x\circ y=0$.
4. $x^2=0$ implica $x=0$, y $Q\cap(-Q)=\{0\}$.

*Demostración.* Usamos $\langle y,x^2\rangle=\langle y\circ x,x\rangle=\langle L(y)x,x\rangle$ y
$\langle x^2,e\rangle=\langle x,x\circ e\rangle=\|x\|^2$; la segunda da (4) junto con (1)(ii), ya que los valores propios
de $-y$ son los $-\lambda_i$ por la unicidad.

(1) (i)$\Rightarrow$(iii): si $y=x^2$ con $x=\sum_i\mu_id_i$, entonces $y=\sum_i\mu_i^2d_i$ y
$L(y)=\sum_i\mu_i^2L(d_i)\ge0$ por el lema F.4(1). (iii)$\Leftrightarrow$(iv) por la primera identidad.
(iii)$\Rightarrow$(ii): $0\le\langle L(y)c_i,c_i\rangle=\lambda_i\|c_i\|^2$. (ii)$\Rightarrow$(i):
$y=(\sum_i\sqrt{\lambda_i}c_i)^2$. Por (iv), $Q$ es una intersección de semiespacios cerrados que pasan por $0$.

(2) (ii)$\Rightarrow$(iii): $L(y)=\sum_i\lambda_iL(c_i)\ge(\min_i\lambda_i)\sum_iL(c_i)=(\min_i\lambda_i)I$.
(iii)$\Rightarrow$(iv): un $x\in Q$ no nulo es $z^2$ con $z\ne0$, y $\langle y,x\rangle=\langle L(y)z,z\rangle>0$.
(iv)$\Rightarrow$(ii): $\langle y,c_i\rangle=\lambda_i\|c_i\|^2$ y $c_i=c_i^2\in Q\setminus\{0\}$. (iii)$\Rightarrow$(i): si
$L(y)\ge mI$ con $m>0$, entonces $L(y')\ge(m-\|L\|\,\|y-y'\|)I$, donde $\|L\|$ es la norma de la aplicación lineal
$b\mapsto L(b)$, de modo que un entorno de $y$ consta de elementos con
$L>0$, que están en $Q$ por (1). (i)$\Rightarrow$(iii): $L(y)\ge0$ por (1); si $\langle L(y)z,z\rangle=0$ con $z\ne0$, entonces
$x=z^2\ne0$ y $\langle y,x\rangle=0$; para $\varepsilon>0$ pequeño, $y-\varepsilon x\in Q$ pero
$\langle y-\varepsilon x,x\rangle=-\varepsilon\|x\|^2<0$, lo que contradice (1)(iv). A continuación, $L(e)=I>0$, y para $y\in Q$,
$y+\varepsilon e\in\Omega$ tiende a $y$, así que $\bar\Omega=Q$; entonces $\Omega^*=\{y:\langle y,x\rangle>0\ \forall x\in Q\setminus\{0\}\}=\Omega$
por (iv). $\Omega=\{L>0\}$ es convexo.

(3) Sea $x=\sum_i\lambda_ic_i$ con $\lambda_i\ge0$. Entonces $0=\sum_i\lambda_i\langle c_i,y\rangle$ con todos los términos $\ge0$
por (1). Si $\lambda_i>0$, entonces $0=\langle y,c_i\rangle=\langle L(y)c_i,c_i\rangle$ con $L(y)\ge0$, luego $L(y)c_i=0$. Por
tanto, $x\circ y=\sum_i\lambda_i\,c_i\circ y=0$. $\square$

**Teorema F.6 (homogeneidad y representación cuadrática).** Sean $V$ un álgebra de Jordan euclídea, y $Q$ y
$\Omega$ como en el teorema F.5.

1. Para todo $b\in V$: $e^{L(b)}Q=Q$, $e^{L(b)}\Omega=\Omega$ y $e^{L(b)}e=\exp b$.
2. Todo $y\in\Omega$ es igual a $\exp(\log y)$. Por tanto, $\Omega$ es homogéneo, y es un cono simétrico.
3. Para $y\in\Omega$, $P(y)=e^{2L(\log y)}$; así, $P(y)Q=Q$, $P(y)\Omega=\Omega$ y $P(y)e=y^2$.
4. $P(y)Q\subseteq Q$ para todo $y\in V$.

*Demostración.* (1) Sea $0\le s<1/\|L(b)\|$ (cualquier $s\ge0$ si $L(b)=0$). El operador simétrico $I-sL(b)$ tiene
valores propios $\ge1-s\|L(b)\|>0$. Sean $x\in Q$ y $z=(I-sL(b))^{-1}x=\sum_j\nu_jd_j$ (teorema F.3). Usando $z\circ d_j=\nu_jd_j$,
$\langle d_l,d_j\rangle=\langle d_l,d_j\circ d_j\rangle=\langle d_l\circ d_j,d_j\rangle=\delta_{lj}\|d_j\|^2$ y
$\langle b,d_j\rangle=\langle b,d_j^2\rangle=\langle L(b)d_j,d_j\rangle\le\|L(b)\|\,\|d_j\|^2$,
$$
0\le\langle x,d_j\rangle=\langle z-s\,b\circ z,d_j\rangle=\nu_j\|d_j\|^2-s\langle b,z\circ d_j\rangle
=\nu_j\big(\|d_j\|^2-s\langle b,d_j\rangle\big),
$$
donde la desigualdad de la izquierda es el teorema F.5(1)(iv). El corchete es $\ge(1-s\|L(b)\|)\|d_j\|^2>0$, así que todos
los $\nu_j\ge0$ y $z\in Q$. Por tanto, para $t\ge0$ y $n>t\|L(b)\|$, $(I-\frac tnL(b))^{-n}Q\subseteq Q$. Cuando $n\to\infty$
este operador tiende a $e^{tL(b)}$ (en una base ortonormal de vectores propios de $L(b)$: $(1-t\beta/n)^{-n}\to e^{t\beta}$),
y $Q$ es cerrado, luego $e^{tL(b)}Q\subseteq Q$. Con $-b$ en lugar de $b$,
$e^{-L(b)}Q\subseteq Q$, así que $e^{L(b)}Q=Q$, y el homeomorfismo $e^{L(b)}$ lleva $\operatorname{int}Q$ sobre sí mismo.
Finalmente, $L(b)^re=b^r$ (potencias de Jordan), así que $e^{L(b)}e=\sum_rb^r/r!=\sum_ie^{\beta_i}d_i=\exp b$ para
$b=\sum_i\beta_id_i$, por el teorema F.3.

(2) Si $y=\sum_i\lambda_ic_i$ con $\lambda_i>0$, entonces $\log y=\sum_i(\log\lambda_i)c_i$ y $\exp\log y=y$. Así,
$y=e^{L(\log y)}e$ con $e^{L(\log y)}\in G(\Omega)$ por (1); la autodualidad es el teorema F.5(2).

(3) Por el lema F.4(4) para los $c_i$ de $y$, tanto $P(y)$ como $e^{2L(\log y)}$ actúan sobre $V_{ij}$ como
$\lambda_i\lambda_j=e^{\log\lambda_i+\log\lambda_j}$. Úsese (1), y $P(y)e=2y\circ y-y^2=y^2$.

(4) Sean $y=\sum_i\lambda_ic_i$, $c=\sum_{\lambda_i\ge0}c_i$ (un idempotente), $s=2c-e=\sum_i\varepsilon_ic_i$ con
$\varepsilon_i=\pm1$, y $|y|=\sum_i|\lambda_i|c_i$. Por el lema F.4(4), $P(s)$, $P(|y|)$ y $P(y)$ actúan sobre $V_{ij}$ como
$\varepsilon_i\varepsilon_j$, $|\lambda_i\lambda_j|$ y $\lambda_i\lambda_j$, luego $P(y)=P(s)P(|y|)$. Para $\eta>0$,
$|y|+\eta e\in\Omega$, así que $P(|y|+\eta e)Q\subseteq Q$ por (3); $P(\cdot)$ es un polinomio y $Q$ es cerrado, luego
$P(|y|)Q\subseteq Q$ ($\eta\to0$). $P(s)$ es un automorfismo de Jordan (lema F.4(3)), así que $P(s)x^2=(P(s)x)^2\in Q$. $\square$

La sección 13.6 demostró $\operatorname{tr}|P(y)x|\le\operatorname{tr}P(y)|x|$ en $V_m$ a partir de tres propiedades. Estas
valen en toda álgebra de Jordan euclídea para la forma lineal $\tau(x)=\langle x,e\rangle$, que es $\operatorname{Tr}$ sobre
$M_{sa}$ y $\operatorname{tr}$ sobre $V_m$ (con los productos escalares de arriba):

**Corolario F.7 (proposición 13.14 en toda álgebra de Jordan euclídea).** Para $x,y$ en un álgebra de Jordan euclídea,
$\tau(|P(y)x|)\le\tau(P(y)|x|)$, con igualdad si $x\in Q$ o $-x\in Q$.

*Demostración.* $x_\pm\in Q$, $x=x_+-x_-$, $|x|=x_++x_-$ (teorema F.3). $A=P(y)x_+$ y $B=P(y)x_-$ están en $Q$
(teorema F.6(4)), y $z=P(y)x=A-B$. Sean $z=\sum_i\mu_ic_i$ y $c=\sum_{\mu_i>0}c_i$; entonces $|z|=z\circ(2c-e)$, y
$\tau(a\circ b)=\langle a,b\rangle$. Como $\langle Q,Q\rangle\ge0$ (teorema F.5(1)) y $c,e-c\in Q$,
$$
\tau(|z|)=\langle A,c\rangle-\langle A,e-c\rangle-\langle B,c\rangle+\langle B,e-c\rangle\le\tau(A)+\tau(B)=\tau(P(y)|x|).
$$
Si $x\in Q$, entonces $x_-=0$ y $z=A\in Q$, luego $|z|=z=P(y)|x|$; análogamente si $-x\in Q$. $\square$

El código comprueba el teorema F.5 para $\operatorname{Herm}_4(\mathbb C)$ y para el cono de Lorentz $Q\subseteq V_4$ (el
emparejamiento negativo es el testigo $c_i$ de (1), $\langle y,c_i\rangle=\lambda_i\|c_i\|^2<0$), el teorema F.6 en $V_4$, donde en las coordenadas
$(a,u)$ los operadores $L(b)$ y $P(y)$ son matrices simétricas $5\times5$, y el corolario F.7 en
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

## F.3 El teorema de Koecher–Vinberg

La sección F.2 demostró que el interior del cono de cuadrados de un álgebra de Jordan euclídea es simétrico. El
recíproco reconstruye el producto de Jordan a partir del cono. Hasta el teorema F.12, $\Omega$ es un cono convexo abierto
autodual en un espacio euclídeo $V$ de dimensión $N$, y $dy$ es la medida de Lebesgue en coordenadas ortonormales (los
cambios ortogonales de coordenadas tienen $|\det|=1$).

**Lema F.8.** Sea $\Omega$ autodual.

1. $0\notin\Omega$; $\langle x,y\rangle\ge0$ para $x,y\in\bar\Omega$; para cada $x\in\Omega$ existe $\delta>0$ con
   $\langle x,y\rangle\ge\delta\|y\|$ para todo $y\in\bar\Omega$.
2. Si $x\in\bar\Omega\setminus\Omega$, existe $u\in\bar\Omega\setminus\{0\}$ con $\langle x,u\rangle=0$. Si $y\in\Omega$,
   $u\in\bar\Omega$ y $t\ge0$, entonces $y+tu\in\Omega$.
3. $\Omega$ es el interior de $\bar\Omega$.
4. $G(\Omega)$ es un subgrupo cerrado de $GL(V)$, y $g\in G(\Omega)$ implica $g^T\in G(\Omega)$.

*Demostración.* (1) Si $0\in\Omega$, entonces $\Omega$ contiene una bola alrededor de $0$, luego $\Omega=V$; pero
$V^*=\emptyset$. Para $x,y\in\Omega$, $\langle x,y\rangle>0$ porque $y\in\Omega^*$; por continuidad, $\langle x,y\rangle\ge0$
sobre $\bar\Omega$. Si la bola cerrada $\bar B(x,\delta)$ está contenida en $\Omega$ e $y\in\bar\Omega\setminus\{0\}$, entonces
$\langle x-\delta y/\|y\|,y\rangle\ge0$. (2) $x\notin\Omega=\Omega^*$ da $u\in\bar\Omega\setminus\{0\}$ con
$\langle u,x\rangle\le0$, y $\langle u,x\rangle\ge0$ por (1). Para $x'\in\bar\Omega\setminus\{0\}$,
$\langle x',y+tu\rangle\ge\langle x',y\rangle>0$, así que $y+tu\in\Omega^*=\Omega$. (3) Si $\bar B(y,\varepsilon)\subseteq\bar\Omega$,
entonces para $x\in\bar\Omega\setminus\{0\}$, $\langle x,y-\varepsilon x/\|x\|\rangle\ge0$, luego $\langle x,y\rangle\ge\varepsilon\|x\|>0$ e
$y\in\Omega^*=\Omega$. (4) Para $g\in G(\Omega)$, $g$ lleva $\bar\Omega\setminus\{0\}$ sobre sí mismo, así que para $y\in\Omega$ y
$x\in\bar\Omega\setminus\{0\}$, $\langle x,g^Ty\rangle=\langle gx,y\rangle>0$: $g^T\Omega\subseteq\Omega$. Lo mismo para $g^{-1}$
da $g^T\Omega=\Omega$. Si $g_k\in G(\Omega)$ convergen a $g\in GL(V)$, entonces $g\Omega\subseteq\bar\Omega$ es abierto, luego
$g\Omega\subseteq\Omega$ por (3); del mismo modo $g^{-1}=\lim g_k^{-1}$ lleva $\Omega$ en $\Omega$, así que $g\Omega=\Omega$. $\square$

La **función característica** de $\Omega$ es $\varphi(x)=\int_{\Omega^*}e^{-\langle x,y\rangle}dy=\int_\Omega e^{-\langle x,y\rangle}dy$.
Pongamos $f=\log\varphi$, sea $H(x)$ el operador simétrico con $\langle H(x)u,v\rangle=D^2f(x)[u,v]$, y
$x^*=-\nabla f(x)$, es decir, $Df(x)[u]=-\langle x^*,u\rangle$.

**Lema F.9 (función característica).** Sea $\Omega$ autodual.

1. $\varphi$ es finita, positiva y $C^\infty$ sobre $\Omega$, y puede derivarse bajo el signo integral.
2. $\varphi(tx)=t^{-N}\varphi(x)$ para $t>0$, y $\varphi(gx)=|\det g|^{-1}\varphi(x)$ para $g\in G(\Omega)$.
3. $H(x)>0$ y $x^*\in\Omega$ para $x\in\Omega$.
4. $H(x)x=x^*$ y $D^3f(x)[x,u,v]=-2\langle H(x)u,v\rangle$; para $g\in G(\Omega)$, $(gx)^*=g^{-T}x^*$ y
   $H(gx)=g^{-T}H(x)g^{-1}$.
5. $\varphi(x)\ge c\|x\|^{-N}$ para una constante $c>0$, y $\varphi(x_k)\to\infty$ si los $x_k\in\Omega$ convergen a un punto de
   $\bar\Omega\setminus\Omega$.

*Demostración.* (1) Sean $x\in\Omega$ y $\delta$ como en el lema F.8(1). Para $\|x'-x\|<\delta/2$ e $y\in\Omega$,
$\langle x',y\rangle\ge\frac\delta2\|y\|$, así que las derivadas $k$-ésimas del integrando en $x'$,
$(-1)^k\prod_{i}\langle u_i,y\rangle e^{-\langle x',y\rangle}$, están acotadas por $\prod_i\|u_i\|\,\|y\|^ke^{-\delta\|y\|/2}$, que es
integrable ($\|y\|\ge N^{-1/2}\sum_i|y_i|$ y Fubini). La convergencia dominada permite derivar bajo el signo integral.
(2) Sustitúyase $y=y'/t$, que lleva $\Omega$ sobre $\Omega$; e $y'=g^Ty$, que lleva $\Omega$ sobre $\Omega$ por el lema
F.8(4), con $\langle gx,y\rangle=\langle x,y'\rangle$ y $dy=|\det g|^{-1}dy'$.
(3) Sea $\mu_x$ la medida de probabilidad $e^{-\langle x,y\rangle}dy/\varphi(x)$ sobre $\Omega$. Por (1),
$Df(x)[u]=-\mathbb E\langle u,y\rangle$ y $D^2f(x)[u,u]=\mathbb E\langle u,y\rangle^2-(\mathbb E\langle u,y\rangle)^2$, la
varianza de $\langle u,y\rangle$. Se anula solo si $\langle u,y\rangle$ es igual a una constante $\gamma$ en $\mu_x$-casi
todo punto, es decir, en casi todo punto de $\Omega$ respecto de Lebesgue (la densidad de $\mu_x$ es positiva); el conjunto
$\{y\in\Omega:\langle u,y\rangle\ne\gamma\}$ es abierto y nulo, luego vacío, así que la función lineal $\langle u,\cdot\rangle$
es constante sobre una bola, lo que obliga a $u=0$. Y $x^*=\mathbb E\,y$, así que $\langle z,x^*\rangle=\mathbb E\langle z,y\rangle>0$
para $z\in\bar\Omega\setminus\{0\}$: $x^*\in\Omega^*=\Omega$.
(4) Por (2), $f(tx)=f(x)-N\log t$, luego $Df(x)[x]=-N$ sobre $\Omega$. Derivando en $x$ en la dirección $u$:
$D^2f(x)[x,u]+Df(x)[u]=0$, es decir, $\langle H(x)x,u\rangle=\langle x^*,u\rangle$; una vez más en la dirección $v$:
$D^3f(x)[x,u,v]+2D^2f(x)[u,v]=0$. Derivar $f(gx)=f(x)-\log|\det g|$ una y dos veces da
$\langle(gx)^*,gu\rangle=\langle x^*,u\rangle$ y $\langle H(gx)gu,gv\rangle=\langle H(x)u,v\rangle$.
(5) $\langle x,y\rangle\le\|x\|\|y\|\le1$ para $\|y\|\le1/\|x\|$, así que $\varphi(x)\ge e^{-1}\operatorname{vol}(\Omega\cap B(0,1/\|x\|))=c\|x\|^{-N}$
con $c=e^{-1}\operatorname{vol}(\Omega\cap B(0,1))>0$. Sea $x_k\to x_0\in\bar\Omega\setminus\Omega$. Por el lema de Fatou,
$\liminf\varphi(x_k)\ge\int_\Omega e^{-\langle x_0,y\rangle}dy$. Tomemos $u$ como en el lema F.8(2), una bola $B(y_1,r)\subseteq\Omega$
y las bolas disjuntas $B_l=B(y_1+2lru/\|u\|,r)$, $l\ge0$; están contenidas en $\Omega$ por el lema F.8(2), y sobre $B_l$,
$\langle x_0,y\rangle=\langle x_0,y-2lru/\|u\|\rangle\le\|x_0\|(\|y_1\|+r)$. Así, la integral es infinita. $\square$

**Lema F.10 (grupos uniparamétricos).** Sean $G$ un subgrupo cerrado de $GL(V)$ y
$\mathfrak g=\{Z:e^{tZ}\in G\text{ para todo }t\in\mathbb R\}$.

1. Si $c:[0,s_0)\to G$ cumple $c(s)=I+sZ+o(s)$ cuando $s\to0^+$, entonces $Z\in\mathfrak g$.
2. $\mathfrak g$ es un álgebra de Lie: un subespacio lineal cerrado bajo $[Z,W]=ZW-WZ$.

*Demostración.* (1) Fijemos $t>0$ y pongamos $A_n=c(t/n)$, $B_n=e^{tZ/n}$. Entonces $\|A_n-B_n\|=o(1/n)$ y
$\|A_n\|,\|B_n\|\le1+K/n$ para algún $K$ y $n$ grande, así que
$$
\|A_n^n-e^{tZ}\|=\Big\|\sum_{r=0}^{n-1}A_n^r(A_n-B_n)B_n^{n-1-r}\Big\|\le n\,e^K\|A_n-B_n\|\to0 .
$$
$A_n^n\in G$ y $G$ es cerrado en $GL(V)$, luego $e^{tZ}\in G$; también $e^{-tZ}=(e^{tZ})^{-1}\in G$. (2) Los múltiplos son
claros. Apliquemos (1) a $c(s)=e^{sZ}e^{sW}=I+s(Z+W)+O(s^2)$, y a $c(s)=e^{aZ}e^{aW}e^{-aZ}e^{-aW}$ con $a=\sqrt s$:
desarrollando cada factor hasta segundo orden, los términos de orden $a$ se cancelan y los de orden $a^2$ suman $ZW-WZ$,
así que $c(s)=I+s[Z,W]+O(s^{3/2})$. $\square$

**Proposición F.11 (el punto base).** Sean $\Omega$ un cono simétrico, $G=G(\Omega)$ y $\mathfrak g$ como en el lema F.10.

1. $x^{**}=x$ y $H(x^*)=H(x)^{-1}$ para $x\in\Omega$.
2. Existe $e\in\Omega$ con $e^*=e$, y todo $e$ así cumple $H(e)=I$.
3. Para $x\in\Omega$ y cualquier $g\in G$ con $ge=x$: $H(x)^{-1}=gg^T\in G$.
4. Si $Z\in\mathfrak g$ y $Ze=0$, entonces $Z^T=-Z$.

*Demostración.* (1) La aplicación $x\mapsto x^*$ es suave de $\Omega$ en $\Omega$ (lema F.9(3)), con derivada $-H(x)$. Como
$g^{-T}\in G$ (lema F.8(4)), el lema F.9(2,4) da
$f(gx)+f((gx)^*)=f(x)-\log|\det g|+f(g^{-T}x^*)=f(x)+f(x^*)$, así que $F(x)=f(x)+f(x^*)$ es constante sobre la órbita
$G x=\Omega$. Por tanto, $0=DF(x)[u]=-\langle x^*,u\rangle+\langle x^{**},H(x)u\rangle$, es decir, $H(x)x^{**}=x^*=H(x)x$
(lema F.9(4)), y $x^{**}=x$ porque $H(x)$ es invertible. Derivar $x^{**}=x$ da $H(x^*)H(x)=I$.
(2) $F_2(x)=f(x)+\frac12\|x\|^2$ cumple $F_2(x)\ge\log c-N\log\|x\|+\frac12\|x\|^2$ (lema F.9(5)), que tiende a $\infty$
cuando $\|x\|\to\infty$. Así, $S=\{x\in\Omega:F_2(x)\le F_2(x_1)\}$ ($x_1\in\Omega$) es acotado; es cerrado en $V$, porque un
punto límite en $\bar\Omega\setminus\Omega$ haría que $f$ no estuviera acotada sobre $S$ (lema F.9(5)). La función continua
$F_2$ alcanza su mínimo sobre el compacto $S$ en algún $e$ del abierto $\Omega$, donde $0=\nabla F_2(e)=-e^*+e$. Entonces
$H(e)=H(e^*)=H(e)^{-1}$ por (1), así que $H(e)^2=I$, y $H(e)=I$ porque $H(e)>0$.
(3) $H(x)=H(ge)=g^{-T}g^{-1}$ por el lema F.9(4) y (2), y $gg^T\in G$ por el lema F.8(4).
(4) $e^{tZ}\in G$ fija $e$, así que $I=H(e)=e^{-tZ^T}H(e)e^{-tZ}=e^{-tZ^T}e^{-tZ}$; derívese en $t=0$. $\square$

**Teorema F.12 (Koecher–Vinberg).** Sea $\Omega$ un cono convexo abierto no vacío en un espacio euclídeo $V$. Entonces
$\Omega$ es un cono simétrico si y solo si existe un producto $\circ$ sobre $V$ tal que $(V,\circ,\langle\cdot,\cdot\rangle)$
es un álgebra de Jordan euclídea y $\Omega$ es el interior de su cono de cuadrados. Para un cono simétrico puede tomarse
como unidad cualquier $e\in\Omega$ con $e^*=e$, y el producto
$$
\langle u\circ v,w\rangle=-\tfrac12\,D^3\log\varphi(e)[u,v,w].
$$

*Demostración.* La implicación «si» la dan los teoremas F.5(2) y F.6(2). Recíprocamente, sean $\Omega$ simétrico y $e$ como en la
proposición F.11(2), y definamos $L(u)$ por $\langle L(u)v,w\rangle=-\frac12D^3f(e)[u,v,w]$ y $u\circ v=L(u)v$.

*Paso 1: conmutatividad, simetría, unidad.* $D^3f(e)$ es simétrica en sus tres argumentos, así que $u\circ v=v\circ u$ y
$L(u)^T=L(u)$. Por el lema F.9(4) y $H(e)=I$, $\langle e\circ v,w\rangle=\langle H(e)v,w\rangle=\langle v,w\rangle$: $e$ es una unidad.

*Paso 2: $L(u)\in\mathfrak g$.* $R(x)=H(x)^{-1}$ está en $G$ (proposición F.11(3)), es suave y $R(e)=I$. Su
derivada es $DR(e)[u]=-H(e)^{-1}DH(e)[u]H(e)^{-1}=-DH(e)[u]=2L(u)$, ya que $\langle DH(e)[u]v,w\rangle=D^3f(e)[u,v,w]$.
El lema F.10(1) para $c(s)=R(e+su)$ ($s\ge0$ pequeño, de modo que $e+su\in\Omega$), con el subgrupo cerrado $G$ (lema F.8(4)),
da $2L(u)\in\mathfrak g$, luego $L(u)\in\mathfrak g$.

*Paso 3: unicidad.* Si $Y\in\mathfrak g$ es simétrico, entonces $Y=L(Ye)$: el simétrico $Y-L(Ye)\in\mathfrak g$ anula a $e$
($L(w)e=w\circ e=w$), así que también es antisimétrico por la proposición F.11(4), luego $0$.

*Paso 4: derivaciones interiores.* Sea $D=[L(x),L(y)]$. Está en $\mathfrak g$ (lema F.10(2)), $D^T=-D$ y
$De=x\circ y-y\circ x=0$. Para $z\in V$, $[D,L(z)]\in\mathfrak g$ es simétrico y $[D,L(z)]e=Dz-L(z)De=Dz$, así que
$[D,L(z)]=L(Dz)$ por el paso 3; aplicado a $w$: $D(z\circ w)=Dz\circ w+z\circ Dw$.

*Paso 5: la identidad de Jordan.* Para $D=[L(y),L(x)]$, $Dx=y\circ x^2-x\circ(x\circ y)$ y
$D(x^2)=y\circ(x\circ x^2)-x\circ(y\circ x^2)$. El paso 4 da $D(x^2)=2x\circ Dx$, es decir,
$y\circ(x\circ x^2)=3x\circ(x^2\circ y)-2x\circ(x\circ(x\circ y))$ para todo $y$:
$$
L(x\circ x^2)=3L(x)L(x^2)-2L(x)^3 .
$$
Todos los $L$ son simétricos, así que trasponiendo se obtiene $L(x\circ x^2)=3L(x^2)L(x)-2L(x)^3$, y $[L(x),L(x^2)]=0$.
Así, $V$ es un álgebra de Jordan euclídea; sea $\Omega_J$ el interior de su cono de cuadrados.

*Paso 6: $\Omega=\Omega_J$.* Para $b\in V$, $e^{L(b)}\in G$ por el paso 2, así que $\exp b=e^{L(b)}e\in\Omega$ (teorema F.6(1)).
Por el teorema F.6(2), $\Omega_J\subseteq\Omega$. Ambos conos son autoduales para $\langle\cdot,\cdot\rangle$ (teorema F.5(2)),
luego $\Omega=\Omega^*\subseteq\Omega_J^*=\Omega_J$. $\square$

La demostración sigue la estructura de J. Faraut, A. Korányi, *Analysis on Symmetric Cones* (Oxford 1994), capítulos
I–III; donde ellos usan la estructura de grupo de Lie de $G(\Omega)$, el lema F.10 y la proposición F.11(3) producen
directamente los grupos uniparamétricos necesarios. M. Koecher (1957) y E. B. Vinberg (1960) demostraron el teorema.

**Ejemplo F.13 (el cono de la programación semidefinida).** Sean $V=\operatorname{Herm}_n(\mathbb R)$, $\langle x,y\rangle=\operatorname{Tr}(xy)$,
$N=n(n+1)/2$ y $\Omega$ las matrices definidas positivas, simétrico por los teoremas F.5 y F.6. Entonces
$\varphi(x)=\varphi(1)\det(x)^{-(n+1)/2}$, $e=\kappa1$ con $\kappa=\sqrt{(n+1)/2}$, y el producto del teorema F.12 es
$u\circ_\kappa v=\frac1{2\kappa}(uv+vu)$, así que $x\mapsto x/\kappa$ es un isomorfismo sobre el álgebra de Jordan usual.

*Demostración.* La descomposición espectral $x=\sum_i\lambda_ip_i$ de una matriz simétrica en sus proyecciones propias es
la descomposición del teorema F.3 (unicidad; $p_i\circ p_j=p_ip_j$ porque los $p_i$ conmutan), así que, por el teorema
F.5(2), el interior del cono de cuadrados es en efecto el conjunto de las matrices definidas positivas. Para $a\in\Omega$,
$P(a)x=axa$ (proposición 13.2(1)) está en $G(\Omega)$ (teorema F.6(3)). En una base ortonormal de vectores propios de $a$,
las matrices $E_{ii}$ y $(E_{ij}+E_{ji})/\sqrt2$ ($i<j$) forman una base ortonormal de vectores propios de $P(a)$ con
valores propios $a_ia_j$, así que $\det P(a)=\prod_{i\le j}a_ia_j=\det(a)^{n+1}$ (cada $a_i$ aparece $n+1$ veces). El lema
F.9(2) da $\varphi(a^2)=\varphi(P(a)1)=\det(a)^{-(n+1)}\varphi(1)$, y todo $x\in\Omega$ es $a^2$ con $a=\sqrt x$. Con
$D\log\det(x)[u]=\operatorname{Tr}(x^{-1}u)$ (a partir de $\det(x+su)=\det x\,\det(1+sx^{-1}u)=\det x\,(1+s\operatorname{Tr}(x^{-1}u)+O(s^2))$)
y $D(x^{-1})[u]=-x^{-1}ux^{-1}$:
$x^*=\kappa^2x^{-1}$; así, $e^*=e$ significa $e^2=\kappa^2 1$ con $e$ definida positiva, es decir, $e=\kappa1$; $D^2f(x)[v,w]=\kappa^2\operatorname{Tr}(x^{-1}vx^{-1}w)$, y
$D^3f(e)[u,v,w]=-\kappa^2\kappa^{-3}\operatorname{Tr}(uvw+vuw)=-\frac2\kappa\operatorname{Tr}((u\circ v)w)$. $\square$

El código calcula $D^3f$ en $e$ por diferenciación automática de $f(x)=-\frac{n+1}2\log\det x$ en coordenadas
ortonormales de $\operatorname{Herm}_3(\mathbb R)$ y comprueba las afirmaciones de la proposición F.11 y del teorema F.12:

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

## F.4 Isomorfismos de orden: el teorema de Kadison en dimensión finita

Sean $\mathcal M=\bigoplus_cM_{k_c}(\mathbb C)$ y $\mathcal N=\bigoplus_dM_{k'_d}(\mathbb C)$, como matrices diagonales por
bloques; $x\ge0$ significa $x=x^*$ con valores propios no negativos (de forma equivalente, todo bloque es semidefinido
positivo), y $[0,1]=\{a\in\mathcal M_{sa}:0\le a\le1\}$. Una **proyección** es $p=p^*=p^2$. Una aplicación lineal $\Phi$ es
**positiva** si $x\ge0$ implica $\Phi(x)\ge0$; un **$*$-isomorfismo de Jordan** es una biyección lineal con
$\Phi(x^*)=\Phi(x)^*$ y $\Phi(x\circ y)=\Phi(x)\circ\Phi(y)$ (sección 13.5). Un punto $a$ de un conjunto convexo $K$ es
**extremo** si no es el punto medio $\frac12(b+b')$ de dos puntos distintos $b,b'\in K$ (de forma equivalente, no es de la
forma $tb+(1-t)b'$ con $b\ne b'$ y $0<t<1$: un tal $a$ es el punto medio de $a\pm\varepsilon(b-b')\in K$ para
$\varepsilon>0$ pequeño).

**Lema F.14.** En $\mathcal M$:

1. Los puntos extremos del conjunto convexo $[0,1]$ son exactamente las proyecciones.
2. Para proyecciones $p,q$: $pq=0$ si y solo si $p+q\le1$ si y solo si $p\circ q=0$.
3. Para proyecciones $p,q$: $p\le q$ si y solo si $pq=p$ si y solo si $p\circ q=p$. Las proyecciones no nulas minimales
   son las proyecciones de rango uno (en algún bloque).
4. Una proyección $p$ es central si y solo si $pxp+(1-p)x(1-p)=x$ para todo $x\in\mathcal M$. Las proyecciones centrales
   son las sumas de las unidades de bloque $P_c$, y las no nulas minimales son las $P_c$.

*Demostración.* (1) Sean $p$ una proyección y $p=\frac12(a+b)$ con $a,b\in[0,1]$. Para $\xi\in\ker p$,
$0=\langle\xi,a\xi\rangle+\langle\xi,b\xi\rangle$ con ambos términos $\ge0$, así que $\|a^{1/2}\xi\|^2=0$ y $a\xi=0$; para
$\xi$ en la imagen de $p$, lo mismo con $1-a,1-b$ da $a\xi=\xi$. Así, $a=p$. Recíprocamente, si $a\in[0,1]$ no es una
proyección, tiene un valor propio $\lambda\in(0,1)$; con la proyección espectral $q_\lambda\in\mathcal M$ (un polinomio en
$a$) y $\varepsilon=\min(\lambda,1-\lambda)$, $a\pm\varepsilon q_\lambda\in[0,1]$ y $a$ es su punto medio.
(2) Si $pq=0$, entonces $qp=(pq)^*=0$, $p+q$ es una proyección y $1-(p+q)=(1-p-q)^2\ge0$; además $p\circ q=0$. Si
$p+q\le1$, entonces para todo $\xi$, $\|qp\xi\|^2=\langle p\xi,qp\xi\rangle\le\langle p\xi,(1-p)p\xi\rangle=0$, luego $qp=0=pq$. Si
$p\circ q=0$, entonces $pq=-qp$, así que $pq=p(pq)=-pqp=(qp)p=qp$, y $pq=-pq=0$.
(3) $pq=p$ implica $qp=p$, $p\circ q=p$ y $q-p=(q-p)^2\ge0$. Si $p\circ q=p$, multiplicar $pq+qp=2p$ por $p$ por la
izquierda y por la derecha da $pq+pqp=2p=pqp+qp$, luego $pq=qp=p$. Si $p\le q$, entonces $p+(1-q)\le1$ y (2) da
$p(1-q)=0$. Una proyección de rango $\ge2$ domina estrictamente la proyección $\xi\xi^*$ sobre una recta de su imagen
dentro de un bloque ($\xi$ un vector unitario de la imagen de un bloque de $p$). Una proyección de rango uno $p$ está en
un solo bloque, y si $0\ne q\le p$, entonces $q=pq$ tiene imagen contenida en la de $p$, luego $q=p$.
(4) $x=pxp+px(1-p)+(1-p)xp+(1-p)x(1-p)$, y $px(1-p)+(1-p)xp=0$ si y solo si ambos términos se anulan (multiplíquese por
$p$ y $1-p$), si y solo si $px=pxp=xp$. Un elemento central conmuta con todas las unidades matriciales de cada bloque, así
que es un escalar $z_c$ sobre el bloque $c$, y $z_c\in\{0,1\}$ para una proyección. $\square$

**Teorema F.15 (Kadison, dimensión finita).** Sea $\Phi:\mathcal M\to\mathcal N$ una biyección lineal con
$\Phi(1)=1$ tal que $\Phi$ y $\Phi^{-1}$ son positivas. Entonces $\Phi$ es un $*$-isomorfismo de Jordan. Recíprocamente,
todo $*$-isomorfismo de Jordan es unital, y tanto él como su inverso son positivos.

*Demostración.* Todo $x\in\mathcal M_{sa}$ es $x_+-x_-$ con $x_\pm\ge0$, así que $\Phi(x)$ es autoadjunto; escribiendo
$x=h+ik$ con $h,k$ autoadjuntos, $\Phi(x^*)=\Phi(h)-i\Phi(k)=\Phi(x)^*$. Lo mismo vale para $\Phi^{-1}$, así que $\Phi$ es
una biyección real-lineal $\mathcal M_{sa}\to\mathcal N_{sa}$ con $x\ge0\iff\Phi(x)\ge0$. Como $\Phi(1)=1$, $\Phi$ lleva
$[0,1]$ sobre $[0,1]$. Una biyección lineal entre conjuntos convexos lleva puntos extremos sobre puntos extremos (si
$\Phi(a)=\frac12(b+b')$, entonces $a=\frac12(\Phi^{-1}b+\Phi^{-1}b')$), así que por el lema F.14(1) $\Phi$ lleva proyecciones
sobre proyecciones y, por el lema F.14(2), $pq=0\iff p+q\le1\iff\Phi(p)+\Phi(q)\le1\iff\Phi(p)\Phi(q)=0$. Para
$x=\sum_i\lambda_ip_i\in\mathcal M_{sa}$ (descomposición espectral, $p_ip_j=0$ para $i\ne j$),
$\Phi(x)^2=\sum_i\lambda_i^2\Phi(p_i)=\Phi(x^2)$. La polarización, $x\circ y=\frac12((x+y)^2-x^2-y^2)$, da
$\Phi(x\circ y)=\Phi(x)\circ\Phi(y)$ sobre $\mathcal M_{sa}$, y la bilinealidad compleja la da sobre
$\mathcal M=\mathcal M_{sa}+i\mathcal M_{sa}$.

Recíprocamente, $\Phi(1)\circ\Phi(y)=\Phi(y)$ para todo $y$ y $\Phi$ es sobreyectiva, así que $\Phi(1)$ es una unidad para
$\circ$, y $\Phi(1)=\Phi(1)\circ1=1$. Todo $x\ge0$ es $h^2$ con $h=\sqrt x$ autoadjunto (teorema 1.6), y
$\Phi(h^2)=\Phi(h)^2\ge0$ porque $\Phi(h)$ es autoadjunto; $\Phi^{-1}$ también es un $*$-isomorfismo de Jordan. $\square$

Kadison demostró el teorema para todas las C\*-álgebras unitales (R. V. Kadison, *A generalized Schwarz inequality and
algebraic invariants for operator algebras*, Ann. of Math. 56 (1952) 494–503); esa versión no se usa en el curso.

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

## F.5 $*$-isomorfismos de Jordan de $\bigoplus_cM_{k_c}$

**Lema F.16.** Sea $\Phi:\mathcal M\to\mathcal N$ un $*$-isomorfismo de Jordan. Entonces $\Phi(1)=1$,
$\Phi(xyx)=\Phi(x)\Phi(y)\Phi(x)$, y $\Phi$ lleva biyectivamente las proyecciones de $\mathcal M$ sobre las de $\mathcal N$,
conservando y reflejando $pq=0$, $p\le q$, la minimalidad, la centralidad y la minimalidad entre las proyecciones
centrales.

*Demostración.* $\Phi(1)=1$ como en el teorema F.15. En toda álgebra asociativa $xyx=2x\circ(x\circ y)-(x\circ x)\circ y$,
lo que da la segunda afirmación. $p=p^*=p\circ p$ se conserva, y también las caracterizaciones del lema F.14(2,3,4), que
usan solo $\circ$, el producto triple $xyx$, $1$ y cuantificadores sobre toda el álgebra; $\Phi^{-1}$ también es un
$*$-isomorfismo de Jordan. $\square$

**Teorema F.17 ($*$-automorfismos de Jordan de $M_n$).** Todo $*$-isomorfismo de Jordan $\Phi:M_n(\mathbb C)\to M_n(\mathbb C)$
es $x\mapsto UxU^*$ o $x\mapsto Ux^TU^*$ para una unitaria $U$. Para $n\ge2$ las primeras aplicaciones son
multiplicativas y las segundas antimultiplicativas, así que ninguna aplicación es de ambas formas.

*Demostración.* Para $n=1$, $\Phi=\mathrm{id}$. Sean $n\ge2$ y $E_{ij}$ las unidades matriciales. Por el lema F.16, los
$\Phi(E_{ii})$ son proyecciones de rango uno ortogonales dos a dos, $\Phi(E_{ii})=\xi_i\xi_i^*$ con $\xi_i$ ortonormales; con
la unitaria $W$, $We_i=\xi_i$, el $*$-isomorfismo de Jordan $\Phi_1(x)=W^*\Phi(x)W$ cumple $\Phi_1(E_{ii})=E_{ii}$.

Sea $i\ne j$. El espacio $S_{ij}=\operatorname{span}\{E_{ij},E_{ji}\}$ es $\{x:pxp=x,\ E_{ii}xE_{ii}=E_{jj}xE_{jj}=0\}$ con
$p=E_{ii}+E_{jj}$, así que $\Phi_1(S_{ij})\subseteq S_{ij}$ por el lema F.16. Escribamos $\Phi_1(E_{ij})=\alpha E_{ij}+\beta E_{ji}$. De
$E_{ij}\circ E_{ij}=0$: $0=\Phi_1(E_{ij})^2=\alpha\beta(E_{ii}+E_{jj})$, así que exactamente uno de $\alpha,\beta$ es no nulo
($\Phi_1$ es inyectiva). Sea $\varepsilon(i,j)=+$ si $\Phi_1(E_{ij})\in\mathbb CE_{ij}$ y $-$ si $\Phi_1(E_{ij})\in\mathbb CE_{ji}$;
$\varepsilon(j,i)=\varepsilon(i,j)$ porque $\Phi_1(E_{ji})=\Phi_1(E_{ij})^*$. Para $i,j,k$ distintos,
$E_{ij}\circ E_{jk}=\frac12E_{ik}\ne0$; si $\varepsilon(i,j)=+$ y $\varepsilon(j,k)=-$, entonces
$\Phi_1(E_{ij})\circ\Phi_1(E_{jk})\in\mathbb C\,E_{ij}\circ E_{kj}=\{0\}$, una contradicción; análogamente para $(-,+)$. Así,
$\varepsilon(i,j)=\varepsilon(j,k)$ siempre que $i,j,k$ sean distintos, y $\varepsilon$ es constante: dos pares o bien
comparten un índice o bien están enlazados a través de un tercer par. Si $\varepsilon=-$, sustitúyase $\Phi_1$ por
$\Phi_1\circ T$ ($T$ la traspuesta, un $*$-automorfismo de Jordan por la proposición 13.9), que tiene $\varepsilon=+$ y
sigue fijando los $E_{ii}$.

Ahora $\Phi_1(E_{ij})=\alpha_{ij}E_{ij}$ para todos $i,j$, con $\alpha_{ii}=1$ y $\alpha_{ji}=\bar\alpha_{ij}$. De
$E_{ij}\circ E_{ji}=\frac12(E_{ii}+E_{jj})$: $|\alpha_{ij}|^2=1$; de $E_{ij}\circ E_{jk}=\frac12E_{ik}$ ($i,j,k$ distintos):
$\alpha_{ij}\alpha_{jk}=\alpha_{ik}$. Con $u_i=\alpha_{i1}$ esto da $\alpha_{ij}=\alpha_{i1}\alpha_{1j}=u_i\bar u_j$ (también para
$i=j$, $i=1$ o $j=1$). Así, $\Phi_1(x)=DxD^*$ con $D=\operatorname{diag}(u_i)$, y $\Phi(x)=WDxD^*W^*$ o
$\Phi(x)=WDx^TD^*W^*$. Si ambas formas coincidieran, $T$ sería multiplicativa, pero $(E_{12}E_{21})^T=E_{11}\ne E_{22}=E_{12}^TE_{21}^T$.
$\square$

Esto demuestra el resultado de N. Jacobson, C. E. Rickart (Trans. Amer. Math. Soc. 69 (1950) 479–502) e I. N.
Herstein (Trans. Amer. Math. Soc. 81 (1956) 331–341) citado en la sección 13.5, en el caso de dimensión finita, para
todo $n$; la proposición 13.11 es el caso $n=2$. Un automorfismo de Jordan de $\operatorname{Herm}_n(\mathbb C)$ se
extiende a un $*$-automorfismo de Jordan de $M_n$ mediante $\Phi(x+iy)=\Phi(x)+i\Phi(y)$, así que también es de una de
las dos formas.

**Teorema F.18 (sumas directas).** Sea $\Phi:\bigoplus_{c=1}^CM_{k_c}\to\bigoplus_{d=1}^{C'}M_{k'_d}$ un
$*$-isomorfismo de Jordan. Entonces $C'=C$, y existen una biyección $\sigma$ con $k'_{\sigma(c)}=k_c$, unitarias
$U_c\in M_{k_c}$ y un conjunto $S$ de bloques tales que
$$
\Phi(x)_{\sigma(c)}=U_cx_cU_c^*\ (c\notin S),\qquad \Phi(x)_{\sigma(c)}=U_cx_c^TU_c^*\ (c\in S).
$$
Recíprocamente, toda aplicación de este tipo es un $*$-isomorfismo de Jordan, y $S$ es único salvo bloques con $k_c=1$.

*Demostración.* Por el lema F.16 y el lema F.14(4), $\Phi(P_c)=P'_{\sigma(c)}$ para una biyección $\sigma$ entre las
proyecciones centrales minimales, así que $C=C'$. Si $x=P_cxP_c$, entonces $\Phi(x)=P'_{\sigma(c)}\Phi(x)P'_{\sigma(c)}$, así
que $\Phi$ lleva el bloque $c$ en el bloque $\sigma(c)$, y $\Phi^{-1}$ lleva el bloque $\sigma(c)$ en el bloque $c$; por
tanto, $k_c^2=k'^2_{\sigma(c)}$, y la restricción es un $*$-isomorfismo de Jordan $M_{k_c}\to M_{k_c}$, de una de las dos
formas por el teorema F.17. Como $x=\sum_cP_cxP_c$, esto determina $\Phi$. El recíproco vale bloque a bloque
(proposición 13.9), y la unicidad de $S$ se sigue del teorema F.17. $\square$

**Corolario F.19 (isomorfismos de orden unitales).** Las biyecciones lineales unitales $\Phi:\mathcal M\to\mathcal N$ tales
que $\Phi$ y $\Phi^{-1}$ son positivas son exactamente las aplicaciones del teorema F.18. Una tal $\Phi$ es $2$-positiva
si y solo si es completamente positiva, si y solo si $S$ no contiene ningún bloque con $k_c\ge2$, si y solo si $\Phi$ es
un $*$-isomorfismo.

*Demostración.* La primera frase son los teoremas F.15 y F.18. Si $S$ no tiene ningún bloque con $k_c\ge2$, toda
aplicación de bloque es $x\mapsto U_cxU_c^*$ (para $k_c=1$, $x^T=x$), así que $\Phi$ es completamente positiva
(proposición 6.4) y multiplicativa. Si $c\in S$ con $k_c\ge2$ y $\Phi$ es $2$-positiva, entonces también lo es la
aplicación de bloque $\Phi_{\sigma(c)c}(a)=U_ca^TU_c^*$ (proposición 6.4), y también $a\mapsto U_c^*\Phi_{\sigma(c)c}(a)U_c=a^T$,
ya que $[A_{ij}]\ge0$ implica $[U^*A_{ij}U]=(1\otimes U)^*[A_{ij}](1\otimes U)\ge0$; esto contradice la proposición 6.1.
Una tal $\Phi$ tampoco es multiplicativa, por el teorema F.17. $\square$

Así, un isomorfismo de orden unital está determinado, salvo las unitarias, por qué bloques iguales intercambia y qué
bloques traspone, y por el corolario F.19 los bloques traspuestos son exactamente lo que impide la positividad completa.
El código construye el automorfismo de Jordan $x\mapsto Ux^TU^*$ de $M_3$ y recupera $U$ (salvo una fase) siguiendo los
pasos de la demostración del teorema F.17:

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

## Ejercicios

1. Demuestre que $\mathbb R^N$ con el producto coordenada a coordenada y el producto escalar estándar es un álgebra de
   Jordan euclídea cuyo cono de cuadrados es $\mathbb R^N_{\ge0}$, y deduzca del corolario F.19 que los isomorfismos de
   orden unitales de $\mathbb C^N=\bigoplus_{c=1}^NM_1$ son las permutaciones de coordenadas.
2. (**Ser positiva no basta.**) Para $0<t<1$ sea $\Psi_t(x)=(1-t)x+\frac tn\operatorname{Tr}(x)1$ sobre $M_n$, $n\ge2$.
   Demuestre que $\Psi_t$ es una biyección positiva unital, que $\Psi_t^{-1}$ no es positiva y que $\Psi_t$ no lleva
   proyecciones en proyecciones.

<details>
<summary>Soluciones</summary>

1. El producto es conmutativo y asociativo, luego de Jordan, con unidad $(1,\dots,1)$, y $L(x)$ es la matriz diagonal
   $\operatorname{diag}(x)$, que es simétrica. Los cuadrados son los vectores con coordenadas no negativas. Por el
   corolario F.19 con todos los $k_c=1$, $\Phi(x)_{\sigma(c)}=U_cx_cU_c^*=x_c$ (las unitarias de $M_1$ son fases), así
   que $\Phi$ permuta las coordenadas; recíprocamente, las permutaciones son isomorfismos de orden unitales.

2. $\Psi_t(1)=(1-t)1+t1=1$; para $x\ge0$ ambos términos son $\ge0$. Como $\operatorname{Tr}\Psi_t(x)=\operatorname{Tr}x$,
   $\Psi_t^{-1}(y)=\frac1{1-t}\big(y-\frac tn\operatorname{Tr}(y)1\big)$. Para una proyección de rango uno $p$, $\Psi_t^{-1}(p)$
   tiene el valor propio $-\frac{t}{n(1-t)}<0$ sobre $\ker p$, y $\Psi_t(p)$ tiene el valor propio $\frac tn\notin\{0,1\}$.
   Esto es coherente con el teorema F.15: la positividad de la inversa es lo que obliga a llevar proyecciones en
   proyecciones.

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

**Volver a:** [capítulo 13](13_jordan_algebras.md), secciones 13.4–13.6: el teorema F.12 demuestra el teorema 13.8, el
teorema F.15 demuestra el teorema 13.10 para las álgebras $\bigoplus_cM_{k_c}(\mathbb C)$, el teorema F.17 extiende la
proposición 13.11 a todo $n$, y el corolario F.7 extiende la proposición 13.14 a toda álgebra de Jordan euclídea.
