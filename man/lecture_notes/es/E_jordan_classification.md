[English](../E_jordan_classification.md) · **Español** · [Français](../fr/E_jordan_classification.md) · [Deutsch](../de/E_jordan_classification.md) · [中文](../zh/E_jordan_classification.md) · [日本語](../ja/E_jordan_classification.md) · [Русский](../ru/E_jordan_classification.md)

# Apéndice E. Clasificación de las álgebras de Jordan euclídeas

El capítulo 13 utilizó la clasificación de Jordan–von Neumann–Wigner de las álgebras de Jordan formalmente reales de
dimensión finita. Este apéndice la demuestra por completo, a partir únicamente de la identidad de Jordan y del álgebra
lineal en dimensión finita. El camino es el clásico (J. Faraut, A. Korányi, *Analysis on Symmetric Cones*, capítulos
II–V; K. McCrimmon, *A Taste of Jordan Algebras*): la asociatividad de potencias y un lema espectral proporcionan
idempotentes y una forma traza asociativa y positiva; la descomposición de Peirce respecto de un marco de Jordan divide
el álgebra en piezas $J_{ij}$ cuyos productos obedecen reglas estrictas; en un álgebra simple de rango $r\ge3$ estas
reglas convierten el espacio de Peirce $J_{13}$ en un álgebra de composición, que por el teorema de Hurwitz es
$\mathbb R,\mathbb C,\mathbb H$ o $\mathbb O$, y el álgebra entera es entonces el $\mathrm{Herm}_r$ de ella. La
sección E.6 demuestra además el teorema de Frobenius sobre las álgebras de división reales asociativas, y la sección
E.8 demuestra que el álgebra de Albert $\mathrm{Herm}_3(\mathbb O)$ es un álgebra de Jordan que no se sumerge en
ninguna álgebra asociativa.

**Aprenderá** por qué la asociatividad de potencias, la positividad de la forma traza y la existencia de una unidad se
siguen solo de la realidad formal; las reglas de Peirce y cómo obligan a las coordenadas de un álgebra simple a
multiplicarse como números; la demostración del teorema de Hurwitz mediante la duplicación de Cayley–Dickson; y por qué
los octoniones dan una nueva álgebra simple solo en rango $3$.

*Cálculo:* las comprobaciones numéricas de este apéndice — [Tutorial, apéndice E](../../tutorial/es/E_jordan_classification.md).

**Convenciones.** Todas las álgebras son espacios vectoriales reales de dimensión finita con un producto bilineal. Un
**álgebra de Jordan** es como en la sección 13.1: $x\circ y=y\circ x$ y $(x^2\circ y)\circ x=x^2\circ(y\circ x)$,
$x^2=x\circ x$; no se supone la existencia de unidad. $L(x)y=x\circ y$, $[A,B]=AB-BA$, y $\mathrm{Tr}$ es la traza
de una aplicación lineal. Una **subálgebra** es un subespacio cerrado bajo $\circ$, un **ideal** un subespacio $I$ con
$J\circ I\subseteq I$, y $J$ es **simple** si $J\ne0$ y sus únicos ideales son $0$ y $J$. Un **idempotente** es un $c$
con $c^2=c$; los idempotentes $c,c'$ son **ortogonales** si $c\circ c'=0$. $J$ es **formalmente real** si
$\sum_ix_i^2=0$ implica que todos los $x_i=0$, y **euclídea** si tiene una unidad $1$ y un producto escalar con
$\langle x\circ y,z\rangle=\langle y,x\circ z\rangle$. Las identidades polinómicas se **linealizan** así: si
$p(t)=\sum_kp_kt^k$ tiene coeficientes vectoriales y $p(t)=0$ para todo $t$ real, entonces todos los $p_k=0$ (se evalúa
en $\deg p+1$ puntos distintos y se invierte la matriz de Vandermonde).

<a id="s-E-1"></a>

## E.1 Identidades de Jordan linealizadas y asociatividad de potencias

<a id="r-E-1"></a>

**Lema E.1.** En toda álgebra de Jordan, para todos los $x,y,z,w$:

1. $[L(x),L(x^2)]=0$ y $[L(x^2),L(z)]+2[L(x\circ z),L(x)]=0$;
2. $L(x^2\circ y)=2L(x\circ y)L(x)+L(x^2)L(y)-2L(x)L(y)L(x)$;
3. $L((x\circ w)\circ y)=L(w\circ y)L(x)+L(x\circ y)L(w)+L(x\circ w)L(y)-L(w)L(y)L(x)-L(x)L(y)L(w)$;
4. $L((x\circ w)\circ y)-L(x\circ(w\circ y))=[L(w),[L(x),L(y)]]$; en particular
   $\mathrm{Tr}\,L((x\circ w)\circ y)=\mathrm{Tr}\,L(x\circ(w\circ y))$.

*Demostración.* (1) La identidad de Jordan dice que $L(x)L(x^2)y=L(x^2)L(x)y$. Sustituyendo $x$ por $x+tz$, el
coeficiente de $t$ es $[L(z),L(x^2)]+2[L(x),L(x\circ z)]=0$. (2) Sustituyendo $x$ por $x+tz$ en la propia identidad de
Jordan, el coeficiente de $t$ es
$$
2((x\circ z)\circ y)\circ x+(x^2\circ y)\circ z=2(x\circ z)\circ(y\circ x)+x^2\circ(y\circ z).
$$
Leídos como aplicaciones lineales de $z$, los cuatro términos son $2L(x)L(y)L(x)z$, $L(x^2\circ y)z$, $2L(x\circ y)L(x)z$
y $L(x^2)L(y)z$. (3) Se sustituye $x$ por $x+sw$ en (2) y se toma el coeficiente de $s$, dividido por $2$. (4) Escríbase
(3) una vez tal cual y otra con $x$ e $y$ intercambiados. Los tres primeros términos de los miembros derechos coinciden, y
$(y\circ w)\circ x=x\circ(w\circ y)$, así que la diferencia es
$-L(w)L(y)L(x)-L(x)L(y)L(w)+L(w)L(x)L(y)+L(y)L(x)L(w)=[L(w),[L(x),L(y)]]$, un conmutador, cuya traza es $0$.
$\square$

*Cálculo:* (4) en $\mathrm{Herm}_3(\mathbb C)$ — [tutorial, sección E.3](../../tutorial/es/E_jordan_classification.md#s-E-3).

<a id="r-E-2"></a>

**Lema E.2 (asociatividad de potencias).** Sean $x^1=x$ y $x^{k+1}=x\circ x^k$. Entonces $x^p\circ x^q=x^{p+q}$ para
todos los $p,q\ge1$, y los operadores $L(x^k)$ conmutan dos a dos. Por tanto, para polinomios reales $f,g$ sin término
constante y $f(x):=\sum_kf_kx^k$, $(fg)(x)=f(x)\circ g(x)$, y $\mathbb R[x]=\{f(x)\}$ es una subálgebra asociativa y
conmutativa.

*Demostración.* $L(x)$ y $L(x^2)$ conmutan por el lema E.1(1); sea $\mathcal A$ el álgebra conmutativa de operadores
que generan. Demostramos por inducción sobre $n\ge2$: $x^p\circ x^q=x^{p+q}$ siempre que $p+q\le n$, y
$L(x^k)\in\mathcal A$ para $k\le n$. Para $n=2$ es la definición. Supongámoslo para $n$, y sea $p+q=n+1$. Si $q=1$,
$x^p\circ x=x^{n+1}$ por conmutatividad. Si $q\ge2$, entonces
$x^p\circ x^q=L(x^p)L(x)x^{q-1}=L(x)L(x^p)x^{q-1}=L(x)x^{n}=x^{n+1}$. A continuación, el lema E.1(2) con $y=x^{n-1}$ y
$x^2\circ x^{n-1}=x^{n+1}$ da $L(x^{n+1})=2L(x^n)L(x)+L(x^2)L(x^{n-1})-2L(x)^2L(x^{n-1})\in\mathcal A$. El resto se
sigue por bilinealidad. $\square$

*En la biblioteca:* el producto de Jordan y las potencias de Jordan de operadores, `jordan.jordan_product`, `jordan.jordan_power` — [tutorial, sección 13.1](../../tutorial/es/13_jordan_algebras.md#s-13-1).

<a id="s-E-2"></a>

## E.2 Idempotentes y la descomposición de Peirce

<a id="r-E-3"></a>

**Lema E.3 (lema espectral).** Sean $J$ formalmente real y $x\in J$. Existen $k\ge0$, reales no nulos distintos
$\lambda_1,\dots,\lambda_k$ e idempotentes no nulos y ortogonales dos a dos $c_1,\dots,c_k\in\mathbb R[x]$ con
$$
x=\sum_{i=1}^k\lambda_ic_i,\qquad x^n=\sum_{i=1}^k\lambda_i^nc_i\quad(n\ge1),\qquad\mathbb R[x]=\mathrm{span}\{c_1,\dots,c_k\}.
$$
En particular, toda álgebra de Jordan formalmente real no nula contiene un idempotente no nulo.

*Demostración.* La realidad formal con un solo término dice que $z^2=0$ implica $z=0$. Sea
$K=\{q\in\mathbb R[t]:(tq)(x)=0\}$. Por el lema E.2, $K$ es un ideal de $\mathbb R[t]$: si $q\in K$ y $s=s_0+s_1$ con
$s_1(0)=0$, entonces $(tqs)(x)=s_0(tq)(x)+(tq)(x)\circ s_1(x)=0$. Como $x,x^2,\dots,x^{\dim J+1}$ son linealmente
dependientes, $K\ne0$, así que $K=\mu\,\mathbb R[t]$ para un $\mu$ mónico. (i) $\mu(0)\neq0$: si $\mu=th$, entonces
$w=(th)(x)$ cumple $w^2=(t\cdot\mu h)(x)=0$, pues $\mu h\in K$; así $w=0$, $h\in K$ y $\mu\mid h$, lo que es imposible
por el grado. (ii) $\mu$ no tiene factores repetidos: si $\mu=g^2h$ con $\deg g\ge1$, entonces $w=(tgh)(x)$ cumple
$w^2=(t\cdot\mu\cdot th)(x)=0$, luego $gh\in K$, imposible por el grado. (iii) $\mu$ no tiene raíces no reales: si
$\mu=((t-\alpha)^2+\beta^2)h$ con $\beta\ne0$, entonces $a=(t(t-\alpha)h)(x)$ y $b=(\beta th)(x)$ cumplen
$a^2+b^2=(t\cdot\mu\cdot th)(x)=0$, luego $b=0$ y $h\in K$, imposible. Por tanto $\mu=\prod_{i=1}^k(t-\lambda_i)$ con
$\lambda_i$ reales, no nulos y distintos. Pongamos
$$
p_i(t)=\frac t{\lambda_i}\prod_{j\ne i}\frac{t-\lambda_j}{\lambda_i-\lambda_j},\qquad c_i=p_i(x).
$$
El polinomio $p_ip_l-\delta_{il}p_i$ es $t$ por un polinomio que se anula en cada $\lambda_j$, es decir, $t$ por un
elemento de $K$; así que $c_i\circ c_l=\delta_{il}c_i$. Del mismo modo $t-\sum_i\lambda_ip_i=t\big(1-\sum_i\lambda_ip_i/t\big)$
y el corchete se anula en cada $\lambda_j$, luego $x=\sum_i\lambda_ic_i$. Si $c_i=0$, entonces $p_i/t\in K$, de grado
$k-1<\deg\mu$; así que $c_i\ne0$. La fórmula para $x^n$ se sigue de la ortogonalidad, y muestra que
$\mathbb R[x]\subseteq\mathrm{span}\{c_i\}$.
$\square$

*Cálculo:* los idempotentes $p_i(x)$, obtenidos solo a partir de potencias de Jordan, en $\mathrm{Herm}_3(\mathbb C)$ — [tutorial, sección E.3](../../tutorial/es/E_jordan_classification.md#s-E-3).

<a id="r-E-4"></a>

**Proposición E.4 (descomposición de Peirce respecto de un idempotente).** Sea $c$ un idempotente de un álgebra de
Jordan $J$.

1. $L(c)(2L(c)-1)(L(c)-1)=0$. Por tanto $J=J_1(c)\oplus J_{1/2}(c)\oplus J_0(c)$ con
   $J_\alpha(c)=\{x:c\circ x=\alpha x\}$, y $c\in J_1(c)$.
2. $[L(c),L(z)]=0$ para $z\in J_0(c)\cup J_1(c)$, y $[L(c),L(y^2)]=0$ para $y\in J_{1/2}(c)$.
3. $J_1\circ J_1\subseteq J_1$, $J_0\circ J_0\subseteq J_0$, $J_1\circ J_0=0$, $(J_1+J_0)\circ J_{1/2}\subseteq J_{1/2}$
   y $J_{1/2}\circ J_{1/2}\subseteq J_1+J_0$ (todo respecto de $c$).
4. Si $c'$ es un idempotente ortogonal a $c$, entonces $c+c'$ es un idempotente, $[L(c),L(c')]=0$ y
   $J_1(c')\subseteq J_0(c)$.

*Demostración.* (1) El lema E.1(2) con $x=y=c$ da $L(c)=3L(c)^2-2L(c)^3$. Una aplicación lineal anulada por un
polinomio de raíces simples es diagonalizable, con valores propios entre esas raíces. (2) El lema E.1(1) con $x=c$ da
$[L(c),L(z)]=2[L(c),L(c\circ z)]$, es decir, $(1-2\alpha)[L(c),L(z)]=0$ para $z\in J_\alpha(c)$. El lema E.1(1) con
$x=y$ y $z=c$ da $[L(y^2),L(c)]+2[L(y\circ c),L(y)]=[L(y^2),L(c)]+[L(y),L(y)]=[L(y^2),L(c)]=0$. (3) Por (2), $L(z)$
lleva cada espacio propio de $L(c)$ en sí mismo si $z\in J_0\cup J_1$. Esto da la primera, la segunda y la cuarta
regla, y $J_1\circ J_0\subseteq J_1\cap J_0=0$. Para $y\in J_{1/2}$, $L(y^2)$ conmuta con $L(c)$, así que
$c\circ y^2=y^2\circ c$ es un vector propio de $L(c)$ de valor propio $1$, como lo es $c$; escribiendo
$y^2=a_1+a_{1/2}+a_0$, $c\circ y^2=a_1+\frac12a_{1/2}$ está en $J_1$, luego $a_{1/2}=0$. Polarizando ($y\to y+y'$) se
obtiene la última regla. (4) $(c+c')^2=c+2c\circ c'+c'=c+c'$, y $c'\in J_0(c)$, así que $L(c)$ y $L(c')$ conmutan por
(2) y son simultáneamente diagonalizables; en particular $J_1(c')$ está generado por vectores propios comunes. Sobre un
vector propio común con valores propios $\alpha$ para $L(c)$ y $1$ para $L(c')$, $L(c+c')$ tiene el valor propio
$\alpha+1\in\{0,\frac12,1\}$, luego $\alpha=0$. $\square$

*Cálculo:* (1) y los valores propios de $L(c)$ en $\mathrm{Herm}_3(\mathbb C)$ — [tutorial, sección E.3](../../tutorial/es/E_jordan_classification.md#s-E-3).

<a id="s-E-3"></a>

## E.3 La forma traza; euclídea = formalmente real

<a id="r-E-5"></a>

**Teorema E.5.** Sean $J$ formalmente real y $\tau(x,y)=\mathrm{Tr}\,L(x\circ y)$.

1. $\tau$ es simétrica, bilineal y **asociativa**: $\tau(x\circ y,z)=\tau(x,y\circ z)$.
2. $\tau$ es definida positiva; para todo idempotente no nulo $c$, $\tau(c,c)=\mathrm{Tr}\,L(c)\ge1$.
3. $J$ tiene un elemento unidad.
4. $J$ es euclídea, con el producto escalar $\tau$. Recíprocamente, toda álgebra de Jordan euclídea es formalmente
   real.

*Demostración.* (1) Lema E.1(4). (2) Por la proposición E.4(1), $\mathrm{Tr}\,L(c)=\dim J_1(c)+\frac12\dim J_{1/2}(c)\ge1$,
pues $0\ne c\in J_1(c)$. Para $x=\sum_i\lambda_ic_i$ como en el lema E.3, $x^2=\sum_i\lambda_i^2c_i$ y
$\tau(x,x)=\mathrm{Tr}\,L(x^2)=\sum_i\lambda_i^2\mathrm{Tr}\,L(c_i)$, que es $>0$ salvo si $k=0$, es decir,
$x=0$. (3) Sean $J\neq0$ y $e$ un idempotente no nulo (lema E.3) con $\dim J_1(e)$ máxima. Si $J_0(e)\ne0$, es una
subálgebra formalmente real (proposición E.4(3)) y contiene un idempotente no nulo $f$; $f$ es ortogonal a $e$, así
que $e+f$ es un idempotente, $J_1(e)\subseteq J_1(e+f)$ (pues $f\circ J_1(e)\subseteq J_0\circ J_1=0$) y
$f\in J_1(e+f)\setminus J_1(e)$, en contradicción con la maximalidad. Así $J_0(e)=0$. Para $y\in J_{1/2}(e)$ se tiene
entonces $y^2\in J_1(e)$, y la asociatividad da
$\frac12\tau(y,y)=\tau(e\circ y,y)=\tau(e,y^2)=\mathrm{Tr}\,L(e\circ y^2)=\mathrm{Tr}\,L(y^2)=\tau(y,y)$, luego
$y=0$. Por tanto $J=J_1(e)$ y $e$ es una unidad. (4) Por (1)–(3), $\tau(x\circ y,z)=\tau(y,x\circ z)$. Recíprocamente,
si $J$ es euclídea y $\sum_ix_i^2=0$, entonces $0=\sum_i\langle x_i\circ x_i,1\rangle=\sum_i\langle x_i,x_i\rangle$.
$\square$

*Cálculo:* (1) y (2) en $\mathrm{Herm}_3(\mathbb C)$ — [tutorial, sección E.3](../../tutorial/es/E_jordan_classification.md#s-E-3); la positividad de $\tau$ para cada tipo — [tutorial, sección E.9](../../tutorial/es/E_jordan_classification.md#s-E-9).

Así, «formalmente real» y «euclídea» describen la misma clase de álgebras de dimensión finita. En adelante $J$ es
formalmente real, $1$ es su unidad y $\tau$ su forma traza.

<a id="s-E-4"></a>

## E.4 Marcos de Jordan y el teorema espectral

Un idempotente no nulo es **primitivo** si no es suma de dos idempotentes ortogonales no nulos. Un **marco de Jordan**
es una familia $(c_1,\dots,c_r)$ de idempotentes primitivos ortogonales dos a dos con $\sum_ic_i=1$.

<a id="r-E-6"></a>

**Lema E.6.** (1) Si $c$ es primitivo, entonces $J_1(c)=\mathbb Rc$. (2) Los idempotentes no nulos ortogonales dos a
dos son linealmente independientes; todo idempotente no nulo es suma de idempotentes primitivos ortogonales dos a dos;
y $J$ tiene un marco de Jordan.

*Demostración.* (1) $J_1(c)$ es una subálgebra formalmente real con unidad $c$ (proposición E.4(3)). Para
$x\in J_1(c)$ escribamos $x=\sum_i\lambda_ic_i$ según el lema E.3, con $c_i\in\mathbb R[x]\subseteq J_1(c)$. Entonces
$c-c_i$ es un idempotente ortogonal a $c_i$ ($(c-c_i)^2=c-2c_i+c_i$ y $c_i\circ(c-c_i)=0$), así que $c-c_i=0$ por la
primitividad. Los $c_i$ distintos son ortogonales, así que hay a lo sumo uno, y $x\in\mathbb Rc$. (2) Si
$\sum_ia_if_i=0$, el producto con $f_j$ da $a_jf_j=0$. Así, un idempotente $e$ es suma de un número máximo posible de
idempotentes no nulos ortogonales dos a dos. Si uno de ellos, $f$, fuera $g+g'$ con $g,g'$ idempotentes no nulos
ortogonales, entonces $g=f\circ g\in J_1(f)\subseteq J_0(f')$ para cualquier otro miembro $f'$ (proposición E.4(4)), y
lo mismo $g'$, así que $f$ podría sustituirse por $g,g'$. Tómese $e=1$. $\square$

<a id="r-E-7"></a>

**Teorema E.7 (teorema espectral).** Para todo $x\in J$ existen un marco de Jordan $(c_1,\dots,c_r)$ y números reales
$\lambda_1,\dots,\lambda_r$ (no necesariamente distintos) con $x=\sum_i\lambda_ic_i$.

*Demostración.* Escribamos $x=\sum_{i=1}^k\lambda_ie_i$ según el lema E.3 y añadamos $e_0=1-\sum_ie_i$, $\lambda_0=0$;
$e_0$ es un idempotente ortogonal a todos los $e_i$. Descompóngase cada $e_i$ no nulo en idempotentes primitivos
ortogonales dos a dos (lema E.6(2)). Una pieza $f$ de $e_i$ está en $J_1(e_i)\subseteq J_0(e_l)$ para $l\ne i$, y una
pieza de $e_l$ está en $J_1(e_l)$; como $J_0(e_l)\circ J_1(e_l)=0$, todas las piezas son ortogonales dos a dos. Su suma
es $\sum_ie_i=1$. $\square$

Que todos los marcos de Jordan de $J$ tienen el mismo número de elementos se demuestra en el corolario E.28.

<a id="r-E-8"></a>

**Teorema E.8 (descomposición de Peirce respecto de un marco).** Sean $(c_1,\dots,c_r)$ un marco de Jordan,
$t_i=\tau(c_i,c_i)$ y, para $i\ne j$, $J_{ij}=J_{ji}=J_{1/2}(c_i)\cap J_{1/2}(c_j)$.

1. Los $L(c_k)$ conmutan, $J=\bigoplus_i\mathbb Rc_i\oplus\bigoplus_{i<j}J_{ij}$, y
   $c_k\circ x=\frac12(\delta_{ki}+\delta_{kj})x$ para $x\in J_{ij}$.
2. Para $i,j,k$ distintos: $J_{ij}\circ J_{jk}\subseteq J_{ik}$; y $J_{ij}\circ J_{kl}=0$ si $\{i,j\}\cap\{k,l\}=\emptyset$.
3. Existe una forma cuadrática definida positiva $N$ sobre $J_{ij}$ con $x^2=N(x)(c_i+c_j)$; con su forma polar
   $\langle\cdot,\cdot\rangle$, $x\circ y=\langle x,y\rangle(c_i+c_j)$ y $\tau(x,y)=(t_i+t_j)\langle x,y\rangle$ para
   $x,y\in J_{ij}$. Si $J_{ij}\ne0$, entonces $t_i=t_j$.

*Demostración.* (1) Los $L(c_k)$ conmutan por la proposición E.4(4), son diagonalizables y suman $L(1)=\mathrm{id}$.
Un vector propio común tiene valores propios $\alpha_k\in\{0,\frac12,1\}$ con $\sum_k\alpha_k=1$: o bien un
$\alpha_i=1$ y los demás se anulan, y entonces el vector está en $J_1(c_i)=\mathbb Rc_i$ (lema E.6(1)); o bien
$\alpha_i=\alpha_j=\frac12$ para un par y los demás se anulan. Así $J_{ij}$ es exactamente el espacio propio común del
segundo tipo, y $J$ es la suma directa de los espacios propios comunes. (2) Sean $x\in J_{ij}$ e $y\in J_{jk}$.
Respecto de $c_i$, $x\in J_{1/2}$ e $y\in J_0$, así que $x\circ y\in J_{1/2}(c_i)$; respecto de $c_k$, igualmente
$x\circ y\in J_{1/2}(c_k)$. Por tanto $x\circ y\in J_{ik}$. Si en cambio $y\in J_{kl}$, $x\circ y$ está en
$J_{1/2}(c_p)$ para $p=i,j,k,l$; no existe un vector propio común con cuatro valores propios $\frac12$, así que
$x\circ y=0$. (3) Para $x,y\in J_{ij}$, $x\circ y\in J_0(c_p)+J_1(c_p)$ para $p=i,j$ y $x\circ y\in J_0(c_k)$ para
$k\notin\{i,j\}$ (proposición E.4(3)); por (1), $x\circ y\in\mathbb Rc_i+\mathbb Rc_j$. Escribamos
$x^2=\alpha c_i+\beta c_j$. La identidad de Jordan con $y=c_i$ da $(x^2\circ c_i)\circ x=\frac12\alpha x$ y
$x^2\circ(c_i\circ x)=\frac12x^2\circ x=\frac14(\alpha+\beta)x$, así que $\alpha=\beta$ si $x\ne0$; pongamos
$N(x)=\alpha$. Como $\tau(c_i,1)=\mathrm{Tr}\,L(c_i)=t_i$, $\tau(x,x)=\tau(x^2,1)=N(x)(t_i+t_j)$; así que $N$ es una
forma cuadrática, definida positiva por el teorema E.5(2), y la polarización da las fórmulas para $x\circ y$ y $\tau$.
Por último, $N(x)t_i=\tau(x^2,c_i)=\tau(x,x\circ c_i)=\frac12\tau(x,x)$, y lo mismo con $j$; así que $t_i=t_j$ si existe
algún $x\ne0$. $\square$

*Cálculo:* los espacios de Peirce como espacios propios comunes de los $L(c_i)$, para cada tipo — [tutorial, sección E.9](../../tutorial/es/E_jordan_classification.md#s-E-9).

<a id="s-E-5"></a>

## E.5 Álgebras simples y rango dos

<a id="r-E-9"></a>

**Teorema E.9.** $J$ es una suma directa $J^{(1)}\oplus\dots\oplus J^{(s)}$ de ideales simples con
$J^{(p)}\circ J^{(q)}=0$ para $p\ne q$. Cada $J^{(p)}$ es un álgebra de Jordan formalmente real.

*Demostración.* Si $I$ es un ideal, también lo es su complemento ortogonal $I^\perp$ para $\tau$:
$\tau(J\circ I^\perp,I)=\tau(I^\perp,J\circ I)=0$. Por tanto $J=I\oplus I^\perp$ e $I\circ I^\perp\subseteq I\cap I^\perp=0$,
así que $J$ es la suma directa de las álgebras $I$ e $I^\perp$, y un ideal de un sumando es un ideal de $J$. Sea $I$ un
ideal no nulo de dimensión mínima; entonces $I$ es simple, y se concluye por inducción sobre la dimensión, aplicada al
álgebra formalmente real $I^\perp$. $\square$

<a id="r-E-10"></a>

**Lema E.10.** Sean $(c_i)$ un marco de Jordan, $i,j,k,l$ distintos, $u\in J_{ij}$, $w\in J_{jk}$ y $z\in J_{kl}$.
Entonces (1) $u\circ(u\circ w)=\frac14N(u)w$, y (2) $(u\circ w)\circ z=u\circ(w\circ z)$.

*Demostración.* (1) Aplíquese el lema E.1(1) con $x=u$ a $c_j$. Con $L(u^2)=N(u)L(c_i+c_j)$, $c_j\circ w=\frac12w$,
$c_j\circ u=\frac12u$ y $c_j\circ(u\circ w)=0$ (pues $u\circ w\in J_{ik}$):
$[L(u^2),L(w)]c_j=\frac14N(u)w-\frac12N(u)w$ y $2[L(u\circ w),L(u)]c_j=(u\circ w)\circ u$; su suma se anula.
(2) Aplíquese el lema E.1(3), con $(x,w,y)$ sustituido por $(u,z,w)$, a $c_j$. Aquí $u\circ z=0$ (teorema E.8(2)),
$c_j\circ u=\frac12u$ y $c_j\circ z=0$, así que $0=\frac12(z\circ w)\circ u-\frac12z\circ(w\circ u)$. $\square$

<a id="r-E-11"></a>

**Proposición E.11.** Sea $(c_1,\dots,c_r)$ un marco de Jordan de $J$, y únanse $i\neq j$ por una arista si
$J_{ij}\ne0$.

1. $J$ es simple si y solo si este grafo es conexo.
2. Si $J$ es simple y $r\ge2$, todos los $J_{ij}$, $i\ne j$, tienen la misma dimensión $d\ge1$, todos los $t_i$ valen
   $t=1+\frac{(r-1)d}2$, y $\dim J=r+d\,\frac{r(r-1)}2=rt$.

*Demostración.* (1) Si el grafo no es conexo, sean $S$ una componente conexa e $I$ el subespacio generado por los
$c_i$ y los $J_{ij}$ con $i,j\in S$, e $I'$ lo mismo para el complemento. Entonces $J=I\oplus I'$ (todos los $J_{ij}$
con $i\in S$, $j\notin S$ se anulan), $I\circ I\subseteq I$ e $I\circ I'=0$ por el teorema E.8; así que $I$ es un
ideal, $\ne0,J$. Recíprocamente, sean el grafo conexo e $I\ne0$ un ideal. $I$ es invariante por los $L(c_k)$, que
conmutan y son diagonalizables, así que es la suma de sus intersecciones con los espacios propios comunes. Si
$0\ne x\in I\cap J_{ij}$, entonces $c_i\circ x^2=N(x)c_i\in I$; así que algún $c_i$ está en $I$. Si $c_i\in I$ y
$0\ne x\in J_{ij}$, entonces $x=2c_i\circ x\in I$ y $c_j=c_j\circ x^2/N(x)\in I$. Por conexión, todos los
$c_j\in I$, luego $1\in I$ e $I=J$. (2) Sean $u\in J_{ij}$ con $N(u)=1$ y $k\notin\{i,j\}$. Por el teorema E.8(2),
$2L(u)$ lleva $J_{jk}$ a $J_{ik}$ y $J_{ik}$ a $J_{jk}$, y por el lema E.10(1) (también con $i,j$ intercambiados) su
cuadrado es la identidad sobre ambos. Así $\dim J_{ik}=\dim J_{jk}$ siempre que $J_{ij}\ne0$. Por tanto, las aristas
$ij$, $jk$ fuerzan una arista $ik$; un grafo conexo es entonces completo, y cambiando un índice cada vez se ve que todas
las $\dim J_{ij}$ son iguales. Por último, $L(c_i)$ tiene el valor propio $1$ sobre $\mathbb Rc_i$, $\frac12$ sobre los
$r-1$ espacios $J_{ij}$ y $0$ en el resto, así que $t_i=\mathrm{Tr}\,L(c_i)=1+\frac{(r-1)d}2$. $\square$

*Cálculo:* la fórmula de la dimensión de (2) para cada tipo — [tutorial, sección E.9](../../tutorial/es/E_jordan_classification.md#s-E-9).

En el resto de la sección, $J$ es simple, $(c_1,\dots,c_r)$ es un marco de Jordan con $r\ge2$, y para $x\in J_{ij}$,
$y\in J_{jk}$ ($i,j,k$ distintos) escribimos
$$
m(x,y)=2\,x\circ y\in J_{ik}.
$$
Por el teorema E.8(3) y la proposición E.11(2), $\tau(x,y)=2t\langle x,y\rangle$ sobre cada $J_{ij}$.

<a id="r-E-12"></a>

**Lema E.12.** Para $i,j,k$ distintos y $x\in J_{ij}$, $y\in J_{jk}$, $z\in J_{ik}$:

1. $N(m(x,y))=N(x)N(y)$ y $m(x,m(x,y))=N(x)y$;
2. $\langle m(x,y),z\rangle=\langle y,m(x,z)\rangle=\langle x,m(y,z)\rangle$;
3. $m(m(x,y),v)=m(x,m(y,v))$ para $l\notin\{i,j,k\}$ y $v\in J_{kl}$.

*Demostración.* (2) Los tres números son $\tau(x\circ y,z)/t$, $\tau(y,x\circ z)/t$ y $\tau(x,y\circ z)/t$, iguales por
asociatividad y conmutatividad. (1) La segunda fórmula es el lema E.10(1). Entonces
$2tN(m(x,y))=4\tau(x\circ y,x\circ y)=4\tau(y,x\circ(x\circ y))=N(x)\tau(y,y)=2tN(x)N(y)$. (3) Lema E.10(2). $\square$

<a id="r-E-13"></a>

**Teorema E.13 (rango uno y dos).** Sea $J$ simple con un marco de Jordan de $r$ elementos. Si $r=1$, entonces
$J=\mathbb Rc_1\cong\mathbb R$. Si $r=2$, entonces $J\cong V_{d+1}$ con $d=\dim J_{12}\ge1$ (factores de espín como en
la sección 13.2).

*Demostración.* $r=1$: teorema E.8(1). $r=2$: $J=\mathbb R1\oplus W$ con $W=\mathbb Rf\oplus J_{12}$, $f=c_1-c_2$. Por
el teorema E.8, $f\circ f=1$, $f\circ x=\frac12x-\frac12x=0$ y $x\circ y=\langle x,y\rangle1$ para $x,y\in J_{12}$. Así
$w\circ w'=B(w,w')1$ para $w,w'\in W$, donde $B(\alpha f+x,\beta f+y)=\alpha\beta+\langle x,y\rangle$ es un producto
escalar. Si $U:W\to\mathbb R^{d+1}$ es una isometría, $a1+w\mapsto(a,Uw)$ lleva
$(a1+w)\circ(b1+w')=(ab+B(w,w'))1+aw'+bw$ al producto de $V_{d+1}$. $\square$

*En la biblioteca:* el factor de espín $V_m$ es `SpinFactor` — [tutorial, sección 13.2](../../tutorial/es/13_jordan_algebras.md#s-13-2).

<a id="s-E-6"></a>

## E.6 Álgebras de composición: los teoremas de Hurwitz y de Frobenius

Un **álgebra de composición** es un álgebra $A$ con unidad $1$ y una forma cuadrática definida positiva $N$ tal que
$N(xy)=N(x)N(y)$. Escribimos $\langle\cdot,\cdot\rangle$ para la forma polar ($N(x)=\langle x,x\rangle$),
$\mathrm{Re}\,x=\langle x,1\rangle$ y $\bar x=2\mathrm{Re}(x)1-x$. Como $N(1)=1$ (lema E.14(1), cuya
demostración no usa la conjugación), esta es la reflexión respecto de la recta $\mathbb R1$, así que $x\mapsto\bar x$ es
una involución isométrica y $\mathrm{Re}\,\bar x=\mathrm{Re}\,x$. Un isomorfismo de álgebras de composición es
un isomorfismo de álgebras isométrico. El **asociador** es $(x,y,z)=(xy)z-x(yz)$.

<a id="r-E-14"></a>

**Lema E.14.** En un álgebra de composición, para todos los $x,y,z,w$:

1. $N(1)=1$, $\langle xy,xz\rangle=N(x)\langle y,z\rangle$, $\langle yx,zx\rangle=\langle y,z\rangle N(x)$ y
   $\langle xy,wz\rangle+\langle wy,xz\rangle=2\langle x,w\rangle\langle y,z\rangle$;
2. $\langle xy,z\rangle=\langle y,\bar xz\rangle=\langle x,z\bar y\rangle$;
3. $\bar x(xy)=N(x)y=(yx)\bar x$; por tanto $\bar xx=x\bar x=N(x)1$, $x^2-2\mathrm{Re}(x)x+N(x)1=0$,
   $xy+yx=2\mathrm{Re}(x)y+2\mathrm{Re}(y)x-2\langle x,y\rangle1$,
   $\bar x(wy)+\bar w(xy)=2\langle x,w\rangle y$ y $(yx)\bar w+(yw)\bar x=2\langle x,w\rangle y$;
4. $\overline{xy}=\bar y\bar x$, $\mathrm{Re}(xy)=\mathrm{Re}(yx)$ y
   $\mathrm{Re}((xy)z)=\mathrm{Re}(x(yz))$;
5. $x(xy)=x^2y$ y $(yx)x=yx^2$; el asociador es alternado (cambia de signo al intercambiar dos argumentos).

*Demostración.* (1) $N(1)=N(1)^2>0$. Polarícese $N(xy)=N(x)N(y)$ en $y$, respectivamente en $x$, y después
$\langle xy,xz\rangle=N(x)\langle y,z\rangle$ en $x$. (2) La última identidad de (1) con $w=1$ se lee
$\langle xy,z\rangle+\langle y,xz\rangle=2\mathrm{Re}(x)\langle y,z\rangle$, es decir,
$\langle xy,z\rangle=\langle y,\bar xz\rangle$; la polarización de $\langle yx,zx\rangle=\langle y,z\rangle N(x)$ en
$x$, con una de las dos variables igual a $1$, da $\langle yx,z\rangle=\langle y,z\bar x\rangle$. (3) Por (2) y (1),
$\langle\bar x(xy),z\rangle=\langle xy,xz\rangle=N(x)\langle y,z\rangle$ para todo $z$, y análogamente por la derecha.
Póngase $y=1$ y sustitúyase $\bar x=2\mathrm{Re}(x)1-x$; después, polarícese. (4) Por (2) y la isometría de la
conjugación,
$\langle\bar y\bar x,z\rangle=\langle\bar x,yz\rangle=\langle\bar x\bar z,y\rangle=\langle\bar z,xy\rangle=\langle z,\overline{xy}\rangle$.
A continuación, $\mathrm{Re}(xy)=\langle y,\bar x\rangle=\langle\bar y,x\rangle=\mathrm{Re}(yx)$, y
$\mathrm{Re}((xy)z)=\langle xy,\bar z\rangle=\langle y,\bar x\bar z\rangle=\langle yz,\bar x\rangle=\mathrm{Re}(x(yz))$.
(5) Por (3), $x(xy)=(2\mathrm{Re}(x)-\bar x)(xy)=(2\mathrm{Re}(x)x-N(x))y=x^2y$, y análogamente por la
derecha. Así $(x,x,y)=(y,x,x)=0$; polarizando, el asociador cambia de signo al intercambiar los dos primeros o los dos
últimos argumentos, y estas transposiciones generan todas las permutaciones. $\square$

*Cálculo:* $N(xy)=N(x)N(y)$, (2), (4) y (5) sobre octoniones aleatorios — [tutorial, sección E.6](../../tutorial/es/E_jordan_classification.md#s-E-6).

**La construcción de Cayley–Dickson.** Para un álgebra $B$ con unidad, involución $b\mapsto\bar b$ y forma cuadrática
$N$, sea $\mathrm{CD}(B)=B\times B$ con
$$
(a,b)(c,d)=(ac-\bar db,\;da+b\bar c),\qquad N(a,b)=N(a)+N(b),\qquad\overline{(a,b)}=(\bar a,-b).
$$
Partiendo de $\mathbb R$ ($\bar x=x$, $N(x)=x^2$) definimos $\mathbb C=\mathrm{CD}(\mathbb R)$,
$\mathbb H=\mathrm{CD}(\mathbb C)$ y $\mathbb O=\mathrm{CD}(\mathbb H)$, de dimensiones $1,2,4,8$, donde
cada duplicación usa la involución y la forma recién definidas (por el lema E.15(1) más abajo, son la conjugación y la
norma del álgebra de composición). En $\mathbb C$, $(a,b)(c,d)=(ac-bd,ad+bc)$: son los números complejos $a+bi$. $B$ es
la subálgebra $B\times0$.

*Cálculo:* el producto de Cayley–Dickson sobre $\mathbb R^{2^k}$ (la función `cd`) — [tutorial, sección E.6](../../tutorial/es/E_jordan_classification.md#s-E-6).

<a id="r-E-15"></a>

**Lema E.15.** (1) Si $B$ es un álgebra de composición asociativa, entonces $\mathrm{CD}(B)$ es un álgebra de
composición con unidad $(1,0)$ y conjugación $\overline{(a,b)}=(\bar a,-b)$. (2) $\mathbb C$ es conmutativa y
asociativa. $M(a,b)=\begin{pmatrix}a&-b\\\bar b&\bar a\end{pmatrix}$ es un homomorfismo de álgebras inyectivo
$\mathbb H\to M_2(\mathbb C)$ con $M(\bar q)=M(q)^*$; así que $\mathbb H$ es asociativa; no es conmutativa. $\mathbb O$
no es asociativa. (3) $\mathbb R,\mathbb C,\mathbb H,\mathbb O$ son álgebras de composición, y
$\mathbb R,\mathbb C,\mathbb H$ son álgebras de división: todo $x\ne0$ tiene el inverso bilateral $\bar x/N(x)$, así
que $xy=0$ implica $x=0$ o $y=0$.

*Demostración.* (1) $N((a,b)(c,d))=N(a)N(c)+N(d)N(b)+N(d)N(a)+N(b)N(c)-2\langle ac,\bar db\rangle+2\langle da,b\bar c\rangle$,
y por el lema E.14(2) $\langle ac,\bar db\rangle=\langle d(ac),b\rangle$ y $\langle da,b\bar c\rangle=\langle(da)c,b\rangle$,
que coinciden cuando $B$ es asociativa. $(1,0)$ es una unidad, y $2\langle(a,b),(1,0)\rangle(1,0)-(a,b)=(\bar a,-b)$.
(2) Desarrollando y usando que $\mathbb C$ es conmutativa, $M(a,b)M(c,d)=M(ac-\bar db,da+b\bar c)$; $M$ es lineal e
inyectiva, y $M(\bar a,-b)=M(a,b)^*$. En $\mathbb H$, $(i,0)(0,1)=(0,i)$ pero $(0,1)(i,0)=(0,-i)$. En
$\mathbb O=\mathbb H\times\mathbb H$ sean $p=((i,0),0)$, $q=((0,1),0)$, $s=(0,1)$; entonces $(pq)s=(0,(0,i))$ y
$p(qs)=(0,(0,-i))$. (3) Por (1), (2) e inducción. Si $x\ne0$, $x(\bar x/N(x))=(\bar x/N(x))x=1$ por el lema E.14(3);
en un álgebra asociativa, $xy=0$ da entonces $y=(\bar x/N(x))(xy)=0$. $\square$

*Cálculo:* la no asociatividad de $\mathbb O$ y el homomorfismo $M$ — [tutorial, sección E.6](../../tutorial/es/E_jordan_classification.md#s-E-6).

<a id="r-E-16"></a>

**Lema E.16 (duplicación).** Sean $A$ un álgebra de composición, $B\subseteq A$ una subálgebra con $1\in B$, e $i\in A$
con $i\perp B$ y $N(i)=1$. Para $a,b,c,d\in B$: $Bi\perp B$, $N(a+bi)=N(a)+N(b)$ y
$$
(a+bi)(c+di)=(ac-\bar db)+(da+b\bar c)\,i .
$$
Así $B\oplus Bi$ es una subálgebra de dimensión $2\dim B$, y $a+bi\mapsto(h(a),h(b))$ es un isomorfismo
$B\oplus Bi\cong\mathrm{CD}(D)$ siempre que $h:B\to D$ sea un isomorfismo de álgebras de composición. Además, $B$
es asociativa; si $B\oplus Bi$ es asociativa, $B$ es conmutativa; si $B\oplus Bi$ es conmutativa, $B=\mathbb R1$.

*Demostración.* $B$ contiene $\bar b=2\mathrm{Re}(b)1-b$ junto con $b$. Por el lema E.14,
$\langle bi,a\rangle=\langle i,\bar ba\rangle=0$, $N(bi)=N(b)$, $\bar i=-i$, y $ci+ic=2\mathrm{Re}(c)i$, es decir,
$ci=i\bar c$, para $c\in B$. Las identidades del lema E.14(3) dan, sucesivamente: con $(x,w,y)=(\bar a,i,b)$,
$a(ib)=i(\bar ab)$; como $ib=\bar bi$ e $i(\bar ab)=\overline{\bar ab}\,i=(\bar ba)i$, sustituir $b$ por $\bar b$ da
$a(bi)=(ba)i$. Con $(y,x,w)=(a,i,\bar c)$, $(ai)c=(a\bar c)i$. Con $x=i$, $i(iy)=-y$. Con $(x,w,y)=(ai,i,\bar d)$, donde
$\overline{ai}=-ai$ (pues $ai\perp1$) y $\langle ai,i\rangle=\mathrm{Re}\,a$ (lema E.14(1)),
$$
-(ai)(i\bar d)-i((ai)\bar d)=2\mathrm{Re}(a)\bar d .
$$
Aquí $i\bar d=di$ e $i((ai)\bar d)=i((ad)i)=i(i\,\overline{ad})=-\bar d\bar a$, así que
$(ai)(di)=\bar d\bar a-2\mathrm{Re}(a)\bar d=-\bar da$. Desarrollando $(a+bi)(c+di)$ con estas cuatro reglas se
obtiene la fórmula, y $N$ y los productos se corresponden con los de $\mathrm{CD}(D)$ mediante $h$. Como $N$ es
multiplicativa sobre $B\oplus Bi\subseteq A$, el cálculo de la demostración del lema E.15(1) da
$\langle d(ac)-(da)c,b\rangle=0$ para todo $b\in B$, así que $B$ es asociativa. Si $B\oplus Bi$ es asociativa,
$(ab)i=a(bi)=(ba)i$, luego $ab=ba$. Si es conmutativa, $ai=ia=\bar ai$, luego $a=\bar a\in\mathbb R1$. $\square$

<a id="r-E-17"></a>

**Teorema E.17 (Hurwitz).** Toda álgebra de composición es isomorfa a exactamente una de
$\mathbb R,\mathbb C,\mathbb H,\mathbb O$.

*Demostración.* Sea $B_0=\mathbb R1\cong\mathbb R$. Mientras $B_k\ne A$, elíjase $i_{k+1}\perp B_k$ con $N(i_{k+1})=1$ y
póngase $B_{k+1}=B_k\oplus B_ki_{k+1}$; por el lema E.16, $B_k\cong\mathrm{CD}^k(\mathbb R)$ tiene dimensión $2^k$.
Si existiera $B_4$, el lema E.16 haría $B_3$ asociativa, luego $B_2$ conmutativa, luego $B_1=\mathbb R1$, que tiene
dimensión $1\ne2$. Así $A=B_k$ para algún $k\le3$. Las cuatro álgebras tienen dimensiones distintas. $\square$

<a id="r-E-18"></a>

**Lema E.18.** Sea $A$ un álgebra de composición. (1) Dos elementos cualesquiera de $A$ están en una subálgebra
asociativa. (2) (**Identidad de Moufang**) $(ux)(yu)=u((xy)u)=(u(xy))u$ para todos los $u,x,y$.

*Demostración.* (1) Sean $p,q\in A$; podemos suponer $p\notin\mathbb R1$ (si no, se intercambian, y si ambos son reales
se usa $\mathbb R1$). $B=\mathrm{span}\{1,p\}$ es una subálgebra ($p^2=2\mathrm{Re}(p)p-N(p)1$), un álgebra
de composición de dimensión $2$, así que $B\cong\mathbb C$ (teorema E.17). Si $q\in B$ hemos terminado; si no, sea $i$
la componente normalizada de $q$ ortogonal a $B$. Por el lema E.16, $B\oplus Bi\cong\mathrm{CD}(\mathbb C)=\mathbb H$
es asociativa y contiene a $p,q$. (2) El asociador es alternado (lema E.14(5)), así que $(u,xy,u)=0$, que es la segunda
igualdad. Sea $X=(u,x,y)$. $\mathrm{Re}\,X=0$ por el lema E.14(4), y $\langle u,X\rangle=\langle u\bar y,ux\rangle-N(u)\langle1,xy\rangle
=N(u)(\langle\bar y,x\rangle-\langle\bar x,y\rangle)=0$; así que la identidad cuadrática polarizada del lema E.14(3) da
$uX+Xu=2\mathrm{Re}(u)X$, es decir, $\bar uX=Xu$. Desarrollando las definiciones se ve que en toda álgebra
$$
(wx,y,z)-(w,xy,z)+(w,x,yz)=w(x,y,z)+(w,x,y)z .
$$
Con $(w,x,y,z)=(u,u,x,y)$, y $(u^2,x,y)=2\mathrm{Re}(u)X$: $(u,ux,y)=2\mathrm{Re}(u)X-uX=\bar uX=Xu$. Con
$(w,x,y,z)=(u,x,y,u)$, $(ux,y,u)=(u,ux,y)=Xu$ y $(x,y,u)=X$ (permutaciones cíclicas): $(u,x,yu)=uX+Xu-Xu=uX$. Por tanto
$(ux)(yu)-u((xy)u)=(u,x,yu)-u(x,y,u)=uX-uX=0$. $\square$

*Cálculo:* la identidad de Moufang sobre octoniones aleatorios — [tutorial, sección E.6](../../tutorial/es/E_jordan_classification.md#s-E-6).

<a id="r-E-19"></a>

**Teorema E.19 (Frobenius).** Toda álgebra real asociativa $A$ de dimensión finita, con unidad $1\ne0$ y sin divisores
de cero, es isomorfa a $\mathbb R$, $\mathbb C$ o $\mathbb H$.

*Demostración.* (i) Para $x\in A$ las potencias $1,x,x^2,\dots$ son linealmente dependientes; sea $\mu$ el polinomio real
mónico de grado mínimo con $\mu(x)=0$. Si $\mu=fg$ con $\deg f,\deg g\ge1$, entonces $f(x)g(x)=0$ obliga a $f(x)=0$ o
$g(x)=0$; así que $\mu$ es irreducible y, por tanto (teorema fundamental del álgebra), $\mu=t-\alpha$ o
$\mu=(t-\alpha)^2+\beta^2$, $\beta>0$. Sea $V=\{v\in A:v^2\in(-\infty,0]\cdot1\}$. Entonces $A=\mathbb R1+V$ (en el
segundo caso $x-\alpha\in V$), y $\mathbb R1\cap V=0$. (ii) $V$ es un subespacio: es cerrado bajo escalares, y para
$u,v\in V$ demostramos que $u+v\in V$. Si $u,v$ son linealmente dependientes, $u+v$ es múltiplo de uno de ellos; así
que supongamos que son independientes. Entonces $1,u,v$ son independientes: si $1=\alpha u+\gamma v$ con, por ejemplo,
$\alpha\ne0$, entonces $u^2=(1-2\gamma v+\gamma^2v^2)/\alpha^2\in\mathbb R1$ obliga a $\gamma=0$ y
$u\in\mathbb R1\cap V=0$. Por (i), $(u\pm v)^2=p_\pm(u\pm v)+q_\pm1$ con $p_\pm,q_\pm$ reales; sumando,
$2u^2+2v^2=(p_++p_-)u+(p_+-p_-)v+(q_++q_-)1$, luego $p_\pm=0$. Si $q_+>0$, entonces
$(u+v-\sqrt{q_+})(u+v+\sqrt{q_+})=0$ y $u+v\in\mathbb R1$, lo cual es imposible; así que $u+v\in V$. (iii) Para
$u,v\in V$, $uv+vu=(u+v)^2-u^2-v^2=-2\langle u,v\rangle1$ define una forma bilineal simétrica con
$\langle u,u\rangle=-u^2\ge0$, y $\langle u,u\rangle=0$ solo si $u^2=0$, es decir, $u=0$: un producto escalar. Los
$u,v\in V$ ortogonales anticonmutan. (iv) Si $V=0$, $A=\mathbb R$. Si $\dim V=1$, $V=\mathbb Ri$ con $i^2=-1$, y
$A\cong\mathbb C$. Si $\dim V\ge2$, tómense $i,j\in V$ ortonormales y $k=ij$. Entonces $k^2=ijij=-i^2j^2=-1$, así que
$k\in V$, $ik+ki=-j+j=0$ y $jk+kj=i-i=0$: $i,j,k$ son ortonormales, con $ij=-ji=k$, $jk=-kj=i$, $ki=-ik=j$. Si hubiera
un vector unitario $l\in V$ ortogonal a $i,j,k$, entonces $lk=(li)j=-i(lj)=(ij)l=kl$ y $lk=-kl$, luego $kl=0$, lo cual
es imposible. Por tanto $A=\mathrm{span}\{1,i,j,k\}$. En $\mathbb H$, los elementos $I=(i,0)$, $J=(0,1)$,
$K=(0,i)$ tienen la misma tabla de multiplicación (el lema E.15(2) da $IJ=K=-JI$, e $I^2=J^2=-1$, $K=IJ$ determinan el
resto por asociatividad), así que $A\cong\mathbb H$. $\square$

<a id="s-E-7"></a>

## E.7 Coordinatización

Para $D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}$, $d=\dim D$ y $n\ge1$, sea
$\mathrm{Herm}_n(D)=\{X\in M_n(D):X_{ji}=\overline{X_{ij}}\}$ con $X\circ Y=\frac12(XY+YX)$, donde
$(XY)_{ik}=\sum_jX_{ij}Y_{jk}$ (para definir esto no hace falta la asociatividad de $D$). Las entradas diagonales cumplen
$\bar x=x$, luego son reales. Sean $E_{ij}$ las unidades matriciales y, para $i\ne j$ y $a\in D$,
$a[ij]=aE_{ij}+\bar aE_{ji}$; así $a[ji]=\bar a[ij]$.

<a id="r-E-20"></a>

**Lema E.20.** $\mathrm{Herm}_n(D)$ es cerrado bajo $\circ$, tiene la unidad $\sum_iE_{ii}$ y dimensión
$n+d\frac{n(n-1)}2$, y, para $a,b\in D$, con $i,j,k$ arbitrarios en las dos primeras fórmulas e $i,j,k,l$ distintos en las demás ($a[ij]$ siempre con $i\ne j$):
$$
E_{ii}\circ E_{jj}=\delta_{ij}E_{ii},\quad E_{kk}\circ a[ij]=\tfrac12(\delta_{ki}+\delta_{kj})a[ij],\quad
a[ij]\circ b[ij]=\langle a,b\rangle(E_{ii}+E_{jj}),\quad a[ij]\circ b[jk]=\tfrac12(ab)[ik],\quad a[ij]\circ b[kl]=0 .
$$
Además, $(X|Y)=\mathrm{Re}\,\mathrm{tr}(XY)$ es un producto escalar con $(X\circ Y|Z)=(Y|X\circ Z)$.

*Demostración.* Multiplicación de matrices: $a[ij]\,b[ij]=a\bar bE_{ii}+\bar abE_{jj}$, y
$a\bar b+b\bar a=\bar ab+\bar ba=2\langle a,b\rangle$ (polarícese $x\bar x=\bar xx=N(x)$); $a[ij]\,b[jk]=abE_{ik}$ y
$b[jk]\,a[ij]=\bar b\bar aE_{ki}=\overline{ab}E_{ki}$; los demás productos son inmediatos. Todos los productos de
elementos de la base están en $\mathrm{Herm}_n(D)$, luego también todos los productos. A continuación,
$(X|X)=\sum_{i,j}\mathrm{Re}(X_{ij}\overline{X_{ij}})=\sum_{i,j}N(X_{ij})>0$ para $X\ne0$. Sea
$T(X,Y,Z)=\mathrm{Re}\,\mathrm{tr}((XY)Z)=\sum_{i,j,k}\mathrm{Re}((X_{ij}Y_{jk})Z_{ki})$. Por el lema
E.14(4), $\mathrm{Re}\,\mathrm{tr}(AB)=\mathrm{Re}\,\mathrm{tr}(BA)$, $T(X,Y,Z)=T(Z,X,Y)$ y, tomando
conjugados ($\mathrm{Re}\,\bar p=\mathrm{Re}\,p$, $\overline{(pq)r}=\bar r(\bar q\bar p)$) y usando que $X,Y,Z$
son hermíticas, $T(X,Y,Z)=\sum\mathrm{Re}(Z_{ik}(Y_{kj}X_{ji}))=T(Z,Y,X)$. Así $T$ es simétrica, y
$(X\circ Y|Z)=\frac12(T(X,Y,Z)+T(Y,X,Z))=T(X,Y,Z)$ es simétrica en $X,Y,Z$. $\square$

*Cálculo:* $\mathrm{Herm}_n(D)$ mediante constantes de estructura (la función `herm`) — [tutorial, sección E.3](../../tutorial/es/E_jordan_classification.md#s-E-3); el producto escalar $(X|Y)$ en $\mathrm{Herm}_3(\mathbb O)$ — [tutorial, sección E.8](../../tutorial/es/E_jordan_classification.md#s-E-8).

<a id="r-E-21"></a>

**Lema E.21 (orientaciones).** Sean $P_{12},P_{13},P_{23}$ espacios con producto escalar, $P_{ji}=P_{ij}$, y para cada
terna ordenada $(i,j,k)$ de elementos distintos de $\{1,2,3\}$ sea $m:P_{ij}\times P_{jk}\to P_{ik}$ bilineal, de modo
que $m(x,y)=m(y,x)$ (el miembro derecho es la aplicación para $(k,j,i)$) y
$\langle m(x,y),z\rangle=\langle y,m(x,z)\rangle$ para $x\in P_{ij}$, $y\in P_{jk}$, $z\in P_{ik}$. Sean $D$ un álgebra
de composición e $\iota_{ij}:D\to P_{ij}$ ($i<j$) biyecciones lineales isométricas con
$m(\iota_{12}a,\iota_{23}b)=\iota_{13}(ab)$; póngase $\iota_{ji}(a)=\iota_{ij}(\bar a)$. Entonces
$m(\iota_{ij}a,\iota_{jk}b)=\iota_{ik}(ab)$ para todas las ternas ordenadas y todos los $a,b\in D$.

Las hipótesis se cumplen (i) para tres índices de un marco de Jordan de un $J$ simple, con los espacios de Peirce y
$m=2\circ$ (conmutatividad y lema E.12(2)); (ii) en $\mathrm{Herm}_3(D)$, con $P_{ij}=\{a[ij]\}$,
$\langle a[ij],b[ij]\rangle=\langle a,b\rangle$ (bien definido, pues la conjugación es isométrica) y $m=2\circ$, ya que
$m(a[ij],b[jk])=(ab)[ik]$ y, para $c[ik]$, $\langle ab,c\rangle=\langle b,\bar ac\rangle$ (lema E.14(2)).

*Demostración.* Si la afirmación vale para $(i,j,k)$, vale para $(k,j,i)$:
$m(\iota_{kj}a,\iota_{ji}b)=m(\iota_{ij}\bar b,\iota_{jk}\bar a)=\iota_{ik}(\bar b\bar a)=\iota_{ik}(\overline{ab})=\iota_{ki}(ab)$.
Así que basta tratar $(2,1,3)$ y $(1,3,2)$, dado $(1,2,3)$. Para todo $c\in D$, usando la adjunción, $(1,2,3)$ y el lema
E.14(2),
$$
\langle m(\iota_{21}a,\iota_{13}b),\iota_{23}c\rangle=\langle\iota_{13}b,m(\iota_{12}\bar a,\iota_{23}c)\rangle=\langle b,\bar ac\rangle=\langle ab,c\rangle=\langle\iota_{23}(ab),\iota_{23}c\rangle,
$$
$$
\langle m(\iota_{13}a,\iota_{32}b),\iota_{12}c\rangle=\langle\iota_{13}a,m(\iota_{12}c,\iota_{23}\bar b)\rangle=\langle a,c\bar b\rangle=\langle ab,c\rangle=\langle\iota_{12}(ab),\iota_{12}c\rangle,
$$
y $\iota_{23},\iota_{12}$ son sobreyectivas. (En la segunda línea la adjunción se usa para $y=\iota_{32}b$ en el papel
de $x$, tras $m(x,y)=m(y,x)$.) $\square$

<a id="r-E-22"></a>

**Teorema E.22 (coordinatización).** Sea $J$ simple con un marco de Jordan $(c_1,\dots,c_r)$, $r\ge3$. Entonces existen
un álgebra de composición $D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}$ de dimensión $d$ y un isomorfismo de álgebras
$\Phi:\mathrm{Herm}_r(D)\to J$ con $\Phi(E_{ii})=c_i$. Si $r\ge4$, $D$ es asociativa, es decir, $D\ne\mathbb O$.

*Demostración.* *Paso 1 (un triángulo).* Elíjanse $e_{12}\in J_{12}$, $e_{23}\in J_{23}$ con $N=1$ (proposición E.11).
Sobre $D'=J_{13}$ defínase $p\cdot q=m(m(p,e_{23}),m(e_{12},q))$, usando los productos $J_{13}\times J_{32}\to J_{12}$,
$J_{21}\times J_{13}\to J_{23}$ y $J_{12}\times J_{23}\to J_{13}$. Por el lema E.12(1), $1'=m(e_{12},e_{23})$ es una
unidad: $p\cdot1'=m(m(p,e_{23}),e_{23})=p$ y $1'\cdot q=m(e_{12},m(e_{12},q))=q$; y $N(p\cdot q)=N(p)N(q)$. Así
$(D',N)$ es un álgebra de composición, y por el teorema E.17 existe un isomorfismo $h:D\to D'$ con
$D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}$. Pongamos $\iota_{13}=h$, $\iota_{12}=m(h(\cdot),e_{23})$ e
$\iota_{23}=m(e_{12},h(\cdot))$; son biyecciones isométricas (lema E.12(1)) y
$m(\iota_{12}a,\iota_{23}b)=h(a)\cdot h(b)=\iota_{13}(ab)$. Por el lema E.21, para $k=3$:

> $(\mathrm T_k)$ existen biyecciones isométricas $\iota_{ij}:D\to J_{ij}$, $i<j\le k$, tales que, con $\iota_{ji}(a)=\iota_{ij}(\bar a)$,
> $m(\iota_{ij}a,\iota_{jl}b)=\iota_{il}(ab)$ para todos los $i,j,l\le k$ distintos y todos los $a,b\in D$.

*Paso 2 (añadir un índice).* Supongamos $(\mathrm T_{k-1})$ con $4\le k\le r$. Elíjase $v\in J_{1k}$ con $N(v)=1$ y
póngase $\iota_{ik}(a)=m(\iota_{i1}a,v)$ para $2\le i<k$ e $\iota_{1k}(a)=m(\iota_{12}a,\iota_{2k}1)$, biyecciones
isométricas por el lema E.12(1). Usando el lema E.12(3) (cada vez con cuatro índices distintos) y $(\mathrm T_{k-1})$:
para $i,j\in\{2,\dots,k-1\}$ distintos,
$$
m(\iota_{ij}a,\iota_{jk}b)=m(\iota_{ij}a,m(\iota_{j1}b,v))=m(m(\iota_{ij}a,\iota_{j1}b),v)=m(\iota_{i1}(ab),v)=\iota_{ik}(ab).
$$
En particular $\iota_{2k}(b)=m(\iota_{23}b,\iota_{3k}1)$ y, por tanto,
$$
m(\iota_{12}a,\iota_{2k}b)=m(m(\iota_{12}a,\iota_{23}b),\iota_{3k}1)=m(\iota_{13}(ab),\iota_{3k}1)
=m(m(\iota_{12}(ab),\iota_{23}1),\iota_{3k}1)=m(\iota_{12}(ab),\iota_{2k}1)=\iota_{1k}(ab).
$$
Para $3\le j<k$, $\iota_{1j}a=m(\iota_{12}a,\iota_{2j}1)$, así que
$m(\iota_{1j}a,\iota_{jk}b)=m(\iota_{12}a,m(\iota_{2j}1,\iota_{jk}b))=m(\iota_{12}a,\iota_{2k}b)=\iota_{1k}(ab)$. Así la
relación vale para todas las ternas $i<j<k$, y el lema E.21, aplicado a cada triángulo $\{i,j,k\}$, da $(\mathrm T_k)$.

*Paso 3 (asociatividad).* Si $r\ge4$, por $(\mathrm T_4)$ y el lema E.12(3),
$\iota_{14}((ab)c)=m(m(\iota_{12}a,\iota_{23}b),\iota_{34}c)=m(\iota_{12}a,m(\iota_{23}b,\iota_{34}c))=\iota_{14}(a(bc))$.

*Paso 4 (el isomorfismo).* Por $(\mathrm T_r)$, defínase la aplicación lineal
$\Phi(\sum_i\alpha_iE_{ii}+\sum_{i<j}a_{ij}[ij])=\sum_i\alpha_ic_i+\sum_{i<j}\iota_{ij}(a_{ij})$; es biyectiva por el
teorema E.8(1), y $\Phi(a[ij])=\iota_{ij}(a)$ para todos los $i\ne j$. Por el teorema E.8 y $(\mathrm T_r)$,
$c_i\circ c_j=\delta_{ij}c_i$, $c_k\circ\iota_{ij}(a)=\frac12(\delta_{ki}+\delta_{kj})\iota_{ij}(a)$,
$\iota_{ij}(a)\circ\iota_{ij}(b)=\langle a,b\rangle(c_i+c_j)$, $\iota_{ij}(a)\circ\iota_{jk}(b)=\frac12\iota_{ik}(ab)$ e
$\iota_{ij}(a)\circ\iota_{kl}(b)=0$: la tabla del lema E.20. Por bilinealidad, $\Phi(X\circ Y)=\Phi(X)\circ\Phi(Y)$.
$\square$

<a id="s-E-8"></a>

## E.8 El álgebra de Albert

<a id="r-E-23"></a>

**Lema E.23 (automorfismos de $\mathrm{Herm}_3(D)$).** Sean $D$ un álgebra de composición y $f_{12},f_{23},f_{13}$
biyecciones lineales isométricas de $D$ con $f_{13}(ab)=f_{12}(a)f_{23}(b)$ para todos los $a,b$. Entonces
$$
\Phi\Big(\sum_i\alpha_iE_{ii}+\sum_{i<j}a_{ij}[ij]\Big)=\sum_i\alpha_iE_{ii}+\sum_{i<j}f_{ij}(a_{ij})[ij]
$$
es un automorfismo de $(\mathrm{Herm}_3(D),\circ)$. Para todo $u\in D$ con $N(u)=1$,
$(f_{12},f_{23},f_{13})=(L_u,R_u,B_u)$ con $L_ua=ua$, $R_ua=au$, $B_ua=uau$ es una terna de este tipo.

*Demostración.* Aplíquese el lema E.21 en el ejemplo (ii) a $\iota_{ij}(a)=f_{ij}(a)[ij]$: entonces
$\Phi(a[ij])=\iota_{ij}(a)$ para todos los $i\ne j$, y $2\iota_{ij}(a)\circ\iota_{jk}(b)=\iota_{ik}(ab)$ para todos los
$i,j,k$ distintos. Como $\Phi$ fija los $E_{ii}$ y los $f_{ij}$ son isométricos, $\Phi$ conserva cada producto de la
tabla del lema E.20 y, por tanto, $\circ$. En cuanto a la terna, $L_u,R_u,B_u$ son isométricas ($N(ua)=N(a)$) y
biyectivas, y $u(ab)u=(ua)(bu)$ es la identidad de Moufang (lema E.18(2)). $\square$

*Cálculo:* el automorfismo $(L_u,R_u,B_u)$ para un octonión aleatorio $u$ de norma uno — [tutorial, sección E.8](../../tutorial/es/E_jordan_classification.md#s-E-8).

<a id="r-E-24"></a>

**Teorema E.24.** (1) Para $D\in\{\mathbb R,\mathbb C,\mathbb H\}$ y $n\ge1$, $\mathrm{Herm}_n(D)$ es un álgebra
de Jordan euclídea, y $\mathrm{Herm}_n(\mathbb H)$ es isomorfa a una subálgebra de Jordan de
$\mathrm{Herm}_{2n}(\mathbb C)$. (2) El **álgebra de Albert** $\mathrm{Herm}_3(\mathbb O)$ (dimensión $27$) y
$\mathrm{Herm}_2(\mathbb O)$ son álgebras de Jordan euclídeas. (3) En cada una de ellas, $(E_{11},\dots,E_{nn})$
es un marco de Jordan con espacios de Peirce $J_{ij}=\{a[ij]:a\in D\}$ de dimensión $d=\dim D$; el álgebra es simple;
$\mathrm{Herm}_1(D)\cong\mathbb R$ y $\mathrm{Herm}_2(D)\cong V_{d+1}$ para
$D=\mathbb R,\mathbb C,\mathbb H,\mathbb O$.

*Demostración.* (1) $M_n(D)$ es un álgebra real asociativa (lema E.15(2)), así que, por el lema E.20 y el lema 13.1,
$\mathrm{Herm}_n(D)$ es un álgebra de Jordan, euclídea por el lema E.20. Sustituir cada entrada
$X_{ij}\in\mathbb H$ por el bloque $2\times2$ $M(X_{ij})$ es un homomorfismo de álgebras inyectivo
$M_n(\mathbb H)\to M_{2n}(\mathbb C)$ (lema E.15(2)) que lleva matrices hermíticas a matrices hermíticas, pues
$M(\bar q)=M(q)^*$; así que conserva $\circ$.

(2) Sea $J(x,y)=(x^2\circ y)\circ x-x^2\circ(y\circ x)$ sobre $\mathrm{Herm}_3(\mathbb O)$; un automorfismo $\Phi$
cumple $J(\Phi x,\Phi y)=\Phi J(x,y)$. Sea $x=\sum_i\alpha_iE_{ii}+p[12]+q[23]+s[13]$. Componemos tres automorfismos de
la forma del lema E.23 de modo que las entradas en $(1,3)$ y $(2,3)$ se vuelvan reales. (a) Si $s\notin\mathbb R$,
escríbase $s=|s|(\cos\theta+v\sin\theta)$ con $v\perp1$, $N(v)=1$, de modo que $v^2=-1$ (lema E.14(3)); para
$u=\cos\frac\theta2-v\sin\frac\theta2$, calculando en la subálgebra conmutativa
$\mathrm{span}\{1,v\}\cong\mathbb C$, $B_us=u^2s=|s|$. (b) Ahora $s$ es real. Elíjase $u\perp1$ con $N(u)=1$ y
$u\perp q$ (posible porque $\dim\mathbb O=8$); entonces $B_us=su^2=-s$ es real y
$\mathrm{Re}(qu)=\langle u,\bar q\rangle=0$, pues $\bar q\in\mathrm{span}\{1,q\}$. (c) Ahora $q\perp1$; si
$q\ne0$, tómese $u=q/|q|$: $qu=q^2/|q|=-|q|$ y $B_us=-s$ son reales. Así existe un automorfismo $\Phi$ tal que $\Phi x$
tiene entradas reales salvo, quizá, $p'$ en $(1,2)$. Como $J(\Phi x,\cdot)$ es lineal, basta demostrar
$J(\Phi x,y)=0$ para $y=E_{ii}$ e $y=c[ij]$, $c\in\mathbb O$. Por el lema E.18(1), $p'$ y $c$ están en una subálgebra
asociativa $Q$; entonces $\Phi x,y\in\mathrm{Herm}_3(Q)$, una subálgebra de $\mathrm{Herm}_3(\mathbb O)$
contenida en el álgebra asociativa $M_3(Q)$, y $J(\Phi x,y)=0$ por el lema 13.1. Por tanto
$J(x,\Phi^{-1}y)=\Phi^{-1}J(\Phi x,y)=0$ para todo $y$: la identidad de Jordan se cumple, y el lema E.20 hace euclídea
el álgebra. Las matrices con tercera fila y tercera columna nulas forman una subálgebra (lema E.20) con unidad
$E_{11}+E_{22}$, que es $\mathrm{Herm}_2(\mathbb O)$.

(3) Por el lema E.20, los $E_{ii}$ son idempotentes ortogonales dos a dos de suma $1$; $L(E_{kk})$ actúa sobre $E_{jj}$
por $\delta_{kj}$ y sobre $a[ij]$ por $\frac12(\delta_{ki}+\delta_{kj})$, así que $J_1(E_{ii})=\mathbb RE_{ii}$ y $J_{1/2}(E_{ii})\cap
J_{1/2}(E_{jj})=\{a[ij]\}$. $E_{ii}$ es primitivo: si $E_{ii}=g+g'$ con idempotentes ortogonales no nulos, entonces
$g=E_{ii}\circ g\in\mathbb RE_{ii}$, luego $g=E_{ii}$ y $g'=0$. El grafo de la proposición E.11 es completo, así que el
álgebra es simple. Para $n=1,2$ aplíquese el teorema E.13. $\square$

*Cálculo:* la identidad de Jordan del álgebra de Albert, y los pasos (a)–(c) de la demostración de (2) sobre un elemento aleatorio — [tutorial, sección E.8](../../tutorial/es/E_jordan_classification.md#s-E-8).

<a id="r-E-25"></a>

**Corolario E.25.** Para $n\ge4$, $\mathrm{Herm}_n(\mathbb O)$ no es un álgebra de Jordan.

*Demostración.* Si lo fuera, se aplicaría el lema E.10(2): su demostración usa solo la identidad de Jordan y los
productos del lema E.20 (con $c_k=E_{kk}$). Entonces
$\frac14((ab)c)[14]=(a[12]\circ b[23])\circ c[34]=a[12]\circ(b[23]\circ c[34])=\frac14(a(bc))[14]$ para todos los
$a,b,c\in\mathbb O$, en contradicción con el lema E.15(2). $\square$

*Cálculo:* el fallo de la identidad de Jordan en $\mathrm{Herm}_4(\mathbb O)$ — [tutorial, sección E.8](../../tutorial/es/E_jordan_classification.md#s-E-8).

<a id="r-E-26"></a>

**Teorema E.26 (Albert).** Sea $D$ un álgebra de composición. Si existen un álgebra real asociativa $A$ y una aplicación
lineal inyectiva $\varphi:\mathrm{Herm}_3(D)\to A$ con $\varphi(x\circ y)=\frac12(\varphi(x)\varphi(y)+\varphi(y)\varphi(x))$,
entonces $D$ es asociativa. En particular, el álgebra de Albert es **excepcional**: no es isomorfa a ninguna subálgebra
de Jordan de un álgebra asociativa, como $M_N(\mathbb C)$.

*Demostración.* Pongamos $e_i=\varphi(E_{ii})$ y $e=e_1+e_2+e_3=\varphi(1)$. Entonces $e_i^2=e_i$, y $e_ie_j+e_je_i=0$
para $i\ne j$; multiplicando por $e_i$ por la izquierda, respectivamente por la derecha, se obtiene
$e_ie_j+e_ie_je_i=0=e_ie_je_i+e_je_i$, luego $e_ie_j=e_je_i=0$. Para $x=\varphi(X)$, $ex+xe=2x$ da $exe=ex=xe=x$, así
que $x=\sum_{p,q}x_{pq}$ con $x_{pq}=e_pxe_q$. Para $X=a[ij]$, las relaciones $e_i\circ x=e_j\circ x=\frac12x$ y
$e_k\circ x=0$ ($k\ne i,j$) dicen $\frac12([p=i]+[q=i])x_{pq}=\frac12x_{pq}$, etc., lo que deja solo $x_{ij}$ y
$x_{ji}$. Defínase, para pares ordenados, $\alpha_{ij}(a)=e_i\varphi(a[ij])e_j$. Para $i,j,k$ distintos,
$\varphi(a[ij])\varphi(b[jk])=\alpha_{ij}(a)\alpha_{jk}(b)$ es el único término en $e_iAe_k$ de
$\varphi(a[ij])\varphi(b[jk])+\varphi(b[jk])\varphi(a[ij])=\varphi((ab)[ik])$, así que
$$
\alpha_{ij}(a)\,\alpha_{jk}(b)=\alpha_{ik}(ab).
$$
De $\varphi(1[13])^2=\varphi(E_{11}+E_{33})$ obtenemos $\alpha_{13}(1)\alpha_{31}(1)=e_1$, $\alpha_{31}(1)\alpha_{13}(1)=e_3$;
de $\varphi(a[13])^2=N(a)(e_1+e_3)$, $\alpha_{13}(a)\alpha_{31}(\bar a)=N(a)e_1$, así que $\alpha_{13}$ es inyectiva.
Ahora
$$
\alpha_{13}(ay)=\alpha_{12}(a)\alpha_{23}(y)=\alpha_{13}(a)\alpha_{32}(1)\alpha_{21}(1)\alpha_{13}(y)=\alpha_{13}(a)\alpha_{31}(1)\alpha_{13}(y),
$$
y por tanto, con $\pi=\alpha_{31}(1)$ y la asociatividad en $A$,
$\alpha_{13}((ab)c)=\alpha_{13}(a)\pi\alpha_{13}(b)\pi\alpha_{13}(c)=\alpha_{13}(a(bc))$. Así $(ab)c=a(bc)$. Por el
lema E.15(2), $\mathbb O$ no es asociativa. $\square$

El argumento original de Albert (A. A. Albert, Ann. of Math. 35 (1934) 65–73) es de este tipo.

<a id="s-E-9"></a>

## E.9 La clasificación

<a id="r-E-27"></a>

**Teorema E.27 (Jordan–von Neumann–Wigner).**

1. Toda álgebra de Jordan formalmente real (equivalentemente, euclídea) de dimensión finita es suma directa de ideales
   simples, y toda álgebra simple de este tipo es isomorfa a $\mathrm{Herm}_n(\mathbb R)$, $\mathrm{Herm}_n(\mathbb C)$ o
   $\mathrm{Herm}_n(\mathbb H)$ para algún $n\ge1$, a un factor de espín $V_m$ con $m\ge2$, o al álgebra de Albert
   $\mathrm{Herm}_3(\mathbb O)$.
2. Las álgebras $\mathbb R$; $V_m$ ($m\ge2$); $\mathrm{Herm}_n(D)$ ($D=\mathbb R,\mathbb C,\mathbb H$, $n\ge3$);
   $\mathrm{Herm}_3(\mathbb O)$ son simples, formalmente reales y no isomorfas dos a dos, así que toda álgebra de
   Jordan simple formalmente real es isomorfa a exactamente una de ellas. Además, $\mathrm{Herm}_1(D)\cong\mathbb R$ y
   $\mathrm{Herm}_2(D)\cong V_{\dim D+1}$ para $D=\mathbb R,\mathbb C,\mathbb H,\mathbb O$.
3. En un álgebra de Jordan simple formalmente real, todos los marcos de Jordan tienen el mismo número $r$ de elementos y
   todos sus espacios de Peirce $J_{ij}$ la misma dimensión $d$: $r=1$ para $\mathbb R$, $(r,d)=(2,m-1)$ para $V_m$ y
   $(n,\dim D)$ para $\mathrm{Herm}_n(D)$, $n\ge2$.

*Demostración.* (1) Por el teorema E.9, $J$ es suma directa de ideales simples formalmente reales. Cada uno de ellos
tiene un marco de Jordan de algún tamaño $r$ (lema E.6(2)); es $\cong\mathbb R=\mathrm{Herm}_1(\mathbb R)$ si
$r=1$, $\cong V_{d+1}$ si $r=2$ (teorema E.13), y $\cong\mathrm{Herm}_r(D)$ con $D\ne\mathbb O$ si $r\ge4$ y con
$D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}$ si $r=3$ (teorema E.22).

(2) Para $\mathrm{Herm}_n(D)$ véase el teorema E.24. $V_m$ es un álgebra de Jordan por el lema 13.3, y euclídea
para $\langle(a,u),(b,v)\rangle=ab+u\cdot v$, pues $\langle x\circ y,z\rangle=abc+c\,u\cdot v+a\,v\cdot w+b\,u\cdot w$ es
simétrica en $x=(a,u)$, $y=(b,v)$, $z=(c,w)$. Para un vector unitario $w$, $c_{\pm w}=\frac12(1,\pm w)$ es un marco de
Jordan (teorema 13.4(1); $c_{\pm w}$ es primitivo porque una suma $c_{w'}+c_{-w'}$ de dos idempotentes ortogonales no
nulos es $1$) con $J_{12}=\{(0,v):v\perp w\}$, de dimensión $m-1\ge1$ (las dimensiones suman $m+1$); así que $V_m$ es
simple por la proposición E.11. Para separar las álgebras, sea $s(J)$ el mayor número de idempotentes no nulos
ortogonales dos a dos (finito por el lema E.6(2)); $s$ y $\dim$ son invariantes por isomorfismo. Claramente
$s(\mathbb R)=1$. En $V_m$ los idempotentes no nulos son $1$ y los $c_w$, con $1\circ c_w=c_w$ y $c_w\circ c_{w'}=0$
solo para $w'=-w$ (teorema 13.4(1)); así que $s(V_m)=2$. Para $D\in\{\mathbb R,\mathbb C,\mathbb H\}$,
$s(\mathrm{Herm}_n(D))=n$: los $E_{ii}$ dan $s\ge n$. Recíprocamente, los idempotentes ortogonales $p_k\ne0$
cumplen $p_kp_l=0$ para $k\ne l$ (como en la demostración del teorema E.26), así que los $D$-submódulos por la derecha
$V_k=p_kD^n$ de $D^n$ son no nulos e independientes: si $\sum_kv_k=0$ con $v_k\in V_k$, entonces
$v_l=p_l\sum_kv_k=0$. Un $V_k$ no nulo contiene algún $v\ne0$ y, con él, $vD$, de dimensión real $d$ ($D$ no tiene
divisores de cero); así que $sd\le nd$. Para $\mathrm{Herm}_3(\mathbb O)$, los $E_{ii}$ dan $s\ge3$; y $s$
idempotentes no nulos ortogonales dos a dos, junto con $1$ menos su suma, se refinan como en la demostración del
teorema E.7 hasta un marco de Jordan de $r'\ge s$ elementos. Por la proposición E.11(2) y el teorema E.22,
$27=r'+d'\frac{r'(r'-1)}2$ con $d'\ge1$, y $d'\in\{1,2,4\}$ si $r'\ge4$. Para $r'=4,5,6$ el miembro derecho está en
$\{10,16,28\}$, $\{15,25,45\}$, $\{21,36,66\}$, y para $r'\ge7$ supera $27$; así que $s\le r'=3$. Los pares
$(s,\dim)$ son $(1,1)$, $(2,m+1)$, $(n,\frac{n(n+1)}2)$, $(n,n^2)$, $(n,n(2n-1))$ para $n\ge3$, y $(3,27)$, todos
distintos, pues $\frac{n(n+1)}2<n^2<n(2n-1)$ y $6,9,15\ne27$. Los isomorfismos en rango $1$ y $2$ son el teorema
E.24(3).

(3) Si $F$ es un marco de Jordan de $r_F$ elementos de un $J$ simple, la demostración de (1), aplicada con $F$, da un
isomorfismo de $J$ con un modelo de rango $r_F$ de (2), y el $s$ de ese modelo es $r_F$. Así $r_F=s(J)$ para todo
marco, y entonces $d=(\dim J-r)/\binom r2$ para $r\ge2$ (proposición E.11(2)). $\square$

*Cálculo:* los marcos de Jordan, las dimensiones de Peirce $(r,d)$ de (3) y $\dim J=r+d\,\frac{r(r-1)}2$ para cada tipo — [tutorial, sección E.9](../../tutorial/es/E_jordan_classification.md#s-E-9).

<a id="r-E-28"></a>

**Corolario E.28 (rango y traza).** En toda álgebra de Jordan formalmente real, todos los marcos de Jordan tienen el
mismo número de elementos, el **rango**. Si $J$ es simple de rango $r$, entonces
$\mathrm{tr}\,x=\frac r{\dim J}\mathrm{Tr}\,L(x)$ cumple $\mathrm{tr}\,c=1$ para todo idempotente primitivo
$c$, por tanto $\mathrm{tr}(\sum_i\lambda_ic_i)=\sum_i\lambda_i$ para toda descomposición del teorema E.7, y
$\mathrm{tr}(x\circ y)$ es un producto escalar asociativo.

*Demostración.* Sean $J=\bigoplus_pJ^{(p)}$ como en el teorema E.9 y $c$ primitivo. Sus componentes
$c^{(p)}\in J^{(p)}$ son idempotentes ortogonales dos a dos ($c^2=\sum_p(c^{(p)})^2$), así que exactamente una es no
nula. Por tanto, los miembros de un marco de Jordan de $J$ que están en $J^{(p)}$ forman un marco de Jordan de
$J^{(p)}$ (suman la componente $1^{(p)}$ de $1$), y el rango de $J$ es la suma de los rangos de los $J^{(p)}$ (teorema
E.27(3)). En un $J$ simple, un $c$ primitivo pertenece a un marco de Jordan (refínese $1-c$ como en la demostración del
teorema E.7), así que $\mathrm{Tr}\,L(c)=t=\dim J/r$ por la proposición E.11(2) (para $r=1$, $J=\mathbb Rc$ y
$\mathrm{Tr}\,L(c)=1$). Por último, $\mathrm{tr}(x\circ y)=\frac r{\dim J}\tau(x,y)$ (teorema E.5). $\square$

<a id="s-E-10"></a>

## E.10 Qué demuestra este apéndice para el capítulo 13

| enunciado del capítulo 13 | demostrado aquí |
|---|---|
| formalmente real $\Leftrightarrow$ euclídea (sección 13.4) | teorema E.5 |
| teorema 13.7: descomposición en ideales simples, la lista, «exactamente una» | teoremas E.9, E.27 |
| $\mathrm{Herm}_2(\mathbb K)\cong V_{\dim\mathbb K+1}$ para $\mathbb K=\mathbb R,\mathbb C,\mathbb H,\mathbb O$ | teorema E.24(3) |
| $\mathrm{Herm}_r(\mathbb H)$ se sumerge en $\mathrm{Herm}_{2r}(\mathbb C)$ | teorema E.24(1) |
| el álgebra de Albert es un álgebra de Jordan y es excepcional | teoremas E.24(2), E.26 |
| $\mathrm{Herm}_r(\mathbb O)$ es un álgebra de Jordan solo para $r\le3$ | teorema E.24(2), corolario E.25 |
| teorema espectral con marcos de Jordan, rango, $\mathrm{tr}\,z=\sum_i\mu_i$, asociatividad de $\mathrm{tr}(x\circ y)$ (sección 13.6) | teorema E.7, corolario E.28 |

**Notas.** La clasificación se debe a P. Jordan, J. von Neumann y E. Wigner, *On an algebraic generalization of the
quantum mechanical formalism*, Ann. of Math. 35 (1934) 29–64, y la excepcionalidad de $\mathrm{Herm}_3(\mathbb O)$
a A. A. Albert (mismo volumen, 65–73). La demostración anterior sigue la estructura de Faraut–Korányi, *Analysis on
Symmetric Cones* (Oxford 1994), capítulos II–V; la demostración del teorema de Hurwitz por duplicación (Nachr. Ges.
Wiss. Göttingen 1898, 309–316) es como en J. H. Conway, D. A. Smith, *On Quaternions and Octonions* (2003) y
T. A. Springer, F. D. Veldkamp, *Octonions, Jordan Algebras and Exceptional Groups* (2000), y la demostración del
teorema de Frobenius (1878) es la de R. S. Palais, Amer. Math. Monthly 75 (1968) 366–368. Lecturas adicionales (no se
usan en ninguna parte del curso): C. M. Glennie, Pacific J. Math. 16 (1966) 47–59, para un enfoque distinto de la
excepcionalidad mediante identidades polinómicas, y *A Taste of Jordan Algebras* (2004) de McCrimmon para la historia y
la teoría general.
