[English](../A_analysis.md) · **Español** · [Français](../fr/A_analysis.md) · [Deutsch](../de/A_analysis.md) · [中文](../zh/A_analysis.md) · [日本語](../ja/A_analysis.md) · [Русский](../ru/A_analysis.md)

# Apéndice A. Hechos de análisis y de probabilidad

Los capítulos utilizan unos cuantos hechos de análisis y de probabilidad que van más allá de los requisitos
previos del curso, pero que no son teoremas propios: los conjuntos de ceros de los polinomios tienen medida nula,
los valores propios ordenados de una matriz hermítica dependen continuamente de ella, las imágenes lineales no
degeneradas de vectores gaussianos tienen densidad, la dirección de un vector gaussiano es uniforme sobre la esfera e
independiente de su longitud, los momentos gaussianos son sumas sobre emparejamientos, y los espaciamientos de puntos
uniformes son variables exponenciales normalizadas. Este apéndice los demuestra, para que los capítulos puedan citarlos
por su número.

**Requisitos previos utilizados.** Álgebra lineal; análisis real (límites, continuidad, compacidad y el teorema de
Bolzano–Weierstrass en $\mathbb R^n$, series de potencias en una variable, la integral de Riemann de una función
continua, la fórmula del cambio de variables); medida e integración (medida de Lebesgue, medidas producto y
Fubini–Tonelli, convergencia monótona y dominada, medidas imagen, densidades, independencia); probabilidad elemental
(la ley $N(0,s^2)$, con densidad $(2\pi s^2)^{-1/2}e^{-y^2/2s^2}$ y función generatriz de momentos
$\mathbb Ee^{aY}=e^{a^2s^2/2}$). Para una aplicación lineal invertible $L$ de $\mathbb R^m$ y un boreliano $E$,
$\lambda_m(LE)=|\det L|\,\lambda_m(E)$, donde $\lambda_m$ es la medida de Lebesgue; este es el caso lineal de la fórmula
del cambio de variables.

En todo lo que sigue, «medible» significa medible Borel, $|x|$ es la norma euclídea de $x\in\mathbb R^n$ o
$\mathbb C^n$, $\|A\|=\sup_{|x|=1}|Ax|$ la norma de operador y $\|A\|_2=(\operatorname{Tr}A^*A)^{1/2}$ la norma de
Frobenius.

## A.1 Medidas e integrales

**Lema A.1 (Dynkin).** Sean $\mu,\nu$ medidas finitas sobre un espacio medible $(\Omega,\mathcal F)$ con
$\mu(\Omega)=\nu(\Omega)$, y sea $\mathcal P\subseteq\mathcal F$ estable por intersecciones finitas y generador de
$\mathcal F$ como $\sigma$-álgebra. Si $\mu(P)=\nu(P)$ para todo $P\in\mathcal P$, entonces $\mu=\nu$.

*Demostración.* Llamemos *sistema de Dynkin* a una familia $\mathcal D\subseteq\mathcal F$ tal que $\Omega\in\mathcal D$,
$B\setminus A\in\mathcal D$ siempre que $A\subseteq B$ estén en $\mathcal D$, y $\bigcup_kA_k\in\mathcal D$ para toda sucesión
creciente $A_k\in\mathcal D$. La familia $\{A\in\mathcal F:\mu(A)=\nu(A)\}$ es un sistema de Dynkin, ya que
$\mu(B\setminus A)=\mu(B)-\mu(A)$ para medidas finitas y las medidas son continuas a lo largo de uniones crecientes. Sea
$d$ la intersección de todos los sistemas de Dynkin que contienen a $\mathcal P$, que es a su vez un sistema de Dynkin.
La familia $\mathcal D_1=\{A\in d:A\cap P\in d\text{ para todo }P\in\mathcal P\}$ es un sistema de Dynkin (para
$A\subseteq A'$, $(A'\setminus A)\cap P=(A'\cap P)\setminus(A\cap P)$) y contiene a $\mathcal P$, luego
$\mathcal D_1=d$. Entonces $\mathcal D_2=\{A\in d:A\cap B\in d\text{ para todo }B\in d\}$ es un sistema de Dynkin que
contiene a $\mathcal P$ (porque $\mathcal D_1=d$), luego $\mathcal D_2=d$: $d$ es estable por intersecciones finitas.
Por tanto $d$ es una $\sigma$-álgebra: $\Omega\setminus A\in d$, $A\cup B=\Omega\setminus\big((\Omega\setminus A)\cap(\Omega\setminus B)\big)\in d$,
y una unión numerable es la unión creciente de uniones finitas. Así,
$\mathcal F=\sigma(\mathcal P)\subseteq d\subseteq\{\mu=\nu\}$. $\square$

**Corolario A.2.** Dos medidas de Borel finitas $\mu,\nu$ sobre un espacio métrico $M$ tales que $\int f\,d\mu=\int f\,d\nu$
para toda $f:M\to\mathbb R$ continua y acotada son iguales.

*Demostración.* Con $f=1$, $\mu(M)=\nu(M)$. Para un cerrado no vacío $K$ (para $K=\emptyset$ no hay nada que demostrar), las funciones continuas
$f_k(x)=\max(0,1-k\operatorname{dist}(x,K))$ toman valores en $[0,1]$ y decrecen hacia la indicatriz de $K$
(para $x\notin K$, $\operatorname{dist}(x,K)>0$). Por convergencia dominada, $\mu(K)=\nu(K)$. Los cerrados son
estables por intersecciones finitas y generan los borelianos; se aplica el lema A.1. $\square$

El lema 5.1 utiliza exactamente este argumento.

**Lema A.3 (medida de Lebesgue sobre un espacio vectorial).** Sea $W$ un espacio vectorial real de dimensión $m$. Para
una biyección lineal $\varphi:\mathbb R^m\to W$ se pone $\lambda_\varphi(E)=\lambda_m(\varphi^{-1}E)$ para todo boreliano $E\subseteq W$.

(a) Para otra biyección lineal $\psi$, $\lambda_\psi=|\det(\psi^{-1}\varphi)|\,\lambda_\varphi$. Así, los conjuntos nulos de
$W$, y el hecho de que un vector aleatorio en $W$ tenga densidad, no dependen de la elección de $\varphi$.

(b) Para toda biyección lineal $L:W\to W$, $\lambda_\varphi(LE)=|\det L|\,\lambda_\varphi(E)$. Si $L$ es ortogonal
para algún producto escalar sobre $W$, entonces $|\det L|=1$ y $L$ conserva $\lambda_\varphi$.

(c) Si un vector aleatorio $X$ en $W$ tiene la densidad $f$ respecto de $\lambda_\varphi$, entonces $LX$ tiene la densidad
$y\mapsto f(L^{-1}y)/|\det L|$.

*Demostración.* (a) $\lambda_\psi(E)=\lambda_m\big((\psi^{-1}\varphi)\varphi^{-1}E\big)=|\det(\psi^{-1}\varphi)|\,\lambda_m(\varphi^{-1}E)$.
(b) $\lambda_\varphi(LE)=\lambda_m\big((\varphi^{-1}L\varphi)\varphi^{-1}E\big)=|\det(\varphi^{-1}L\varphi)|\,\lambda_\varphi(E)$, y
$\det(\varphi^{-1}L\varphi)=\det L$ no depende de la base. En una base ortonormal, una $L$ ortogonal tiene una
matriz $O$ con $O^TO=1$, luego $(\det L)^2=1$. (c) Por (b), la medida imagen de $\lambda_\varphi$ por $L$ es
$|\det L|^{-1}\lambda_\varphi$, así que, para todo boreliano $E$,
$\mathbb P(LX\in E)=\int1_E(Lx)f(x)\,d\lambda_\varphi(x)=|\det L|^{-1}\int1_E(y)f(L^{-1}y)\,d\lambda_\varphi(y)$. $\square$

*Convenciones.* «Medida de Lebesgue» sobre un espacio vectorial real significa cualquier $\lambda_\varphi$; por (a), la elección
nunca importa para los conjuntos nulos ni para las densidades. Sobre $\mathbb C^n$ y sobre $M_{a\times b}(\mathbb C)$ usamos como
coordenadas las partes real e imaginaria de las entradas; son ortonormales para el producto escalar real
$\operatorname{Re}\operatorname{Tr}(X^*Y)$, para el cual toda aplicación $X\mapsto VXW$ con $V,W$ unitarias es ortogonal
(conserva $\operatorname{Tr}(X^*Y)$). Sobre las matrices hermíticas $\operatorname{Herm}_n$ usamos las coordenadas
$h_{ii}$ y $\operatorname{Re}h_{ij},\operatorname{Im}h_{ij}$ ($i<j$); no son ortonormales para el producto escalar real
$\operatorname{Tr}(XY)$, pero por (b) las aplicaciones $H\mapsto UHU^*$ ($U$ unitaria), que son ortogonales para él,
conservan igualmente la medida de Lebesgue. Del mismo modo, sobre las matrices simétricas reales $\operatorname{Sym}_n(\mathbb R)$ usamos las
coordenadas $h_{ii}$ y $h_{ij}$ ($i<j$), y las aplicaciones $H\mapsto OHO^T$ ($O\in O(n)$), ortogonales para $\operatorname{Tr}(XY)$,
conservan la medida de Lebesgue.

**Lema A.4 (integrales de funciones con valores vectoriales).** Sea $V$ un espacio vectorial real o complejo de dimensión
finita y $F:[a,b]\to V$ continua. Se define $\int_a^bF(t)\,dt$ coordenada a coordenada en una base (partes real e
imaginaria por separado). Entonces:

(a) para particiones $a=t_0<\dots<t_N=b$ con puntos intermedios $\tau_i\in[t_{i-1},t_i]$, las sumas de Riemann
$\sum_iF(\tau_i)(t_i-t_{i-1})$ convergen hacia $\int_a^bF$ cuando $\max_i(t_i-t_{i-1})\to0$; en particular, la integral
no depende de la base;

(b) $L\int_a^bF=\int_a^bLF$ para toda aplicación ($\mathbb R$-)lineal $L:V\to V'$ en un espacio $V'$ de dimensión finita;

(c) $\|\int_a^bF\|\le\int_a^b\|F(t)\|\,dt$ para toda norma sobre $V$;

(d) si todos los $F(t)$ están en un convexo cerrado $K$, entonces $\frac1{b-a}\int_a^bF\in K$; si además $K$ es un cono,
$\int_a^bF\in K$;

(e) $t\mapsto\int_a^tF$ es derivable con derivada $F(t)$; $\int_a^bF=\int_a^cF+\int_c^bF$, y
$\int_a^bF(s+c)\,ds=\int_{a+c}^{b+c}F(s)\,ds$.

*Demostración.* Cada coordenada de $F$ es una función real continua, así que (a) y (e) valen coordenada a coordenada;
(b) se sigue de (a), pues $L$ transforma sumas de Riemann en sumas de Riemann y es continua. Una norma es continua
($|\|x\|-\|y\||\le\|x-y\|\le\sum_k|x_k-y_k|\,\|e_k\|$), luego $t\mapsto\|F(t)\|$ es continua, y (c) se sigue de
(a) y de la desigualdad triangular para las sumas de Riemann. Para (d) se toman las particiones en $N$ partes iguales:
las sumas de Riemann divididas por $b-a$ son combinaciones convexas $\frac1N\sum_iF(\tau_i)\in K$, y $K$ es cerrado. $\square$

Esto justifica el lema 8.4(c) (sumas de Riemann de aplicaciones CP) y las medias de Cesàro
$\sigma_\tau=\frac1\tau\int_0^\tau T_s\rho\,ds$ de la proposición 8.12: $\sigma_\tau$ es un estado por (d), y
$\|T_t\sigma_\tau-\sigma_\tau\|_1\le2t/\tau$ por (b), (c) y (e).

## A.2 Conjuntos de ceros de funciones analíticas

Para un multiíndice $\alpha\in\mathbb N^n$ escribimos $|\alpha|=\alpha_1+\dots+\alpha_n$ e $y^\alpha=y_1^{\alpha_1}\cdots y_n^{\alpha_n}$.
Para $a\in\mathbb R^n$ y $0<\rho\le\infty$, sea $Q_\rho(a)=\{x\in\mathbb R^n:\max_i|x_i-a_i|<\rho\}$ el cubo
abierto ($Q_\infty(a)=\mathbb R^n$).

**Definición.** Una función $f:Q_\rho(a)\to\mathbb R$ es una **serie de potencias sobre $Q_\rho(a)$** si
$f(x)=\sum_\alpha c_\alpha(x-a)^\alpha$ con $\sum_\alpha|c_\alpha|r^{|\alpha|}<\infty$ para todo $0<r<\rho$.

Los polinomios son series de potencias sobre $\mathbb R^n$ (con un número finito de $c_\alpha\neq0$); en varias variables, el curso
aplica el teorema A.6 solo a polinomios. (Una función que solo es localmente una serie de potencias sobre un abierto conexo
requeriría además el principio de identidad en varias variables, que el curso no utiliza.) Sobre el cubo cerrado de radio
$r<\rho$ los términos están acotados por $|c_\alpha|r^{|\alpha|}$, así que la serie converge absoluta y uniformemente
allí (criterio M de Weierstrass): su suma no depende del orden de los términos, y $f$ es continua.

**Lema A.5 (una variable).** Sea $g(t)=\sum_kc_k(t-a)^k$ una serie de potencias sobre $I=(a-\rho,a+\rho)$.

(a) Para $b\in I$ y $|t-b|<\rho_b=\rho-|b-a|$, $g(t)=\sum_jd_j(t-b)^j$ con
$d_j=\sum_{k\ge j}\binom kjc_k(b-a)^{k-j}$, una serie de potencias sobre $(b-\rho_b,b+\rho_b)$.

(b) Si algún $d_j\neq0$, entonces $g$ no tiene ceros en un entorno perforado de $b$.

(c) Si los ceros de $g$ se acumulan en un punto de $I$, entonces $g=0$ sobre $I$. Si $g\not\equiv0$ sobre $I$, su conjunto
de ceros es numerable y, por tanto, de medida de Lebesgue nula. Si $g\equiv0$ sobre $I$, entonces todos los $c_k=0$.

*Demostración.* (a) Para $0\le s<\rho_b$,
$\sum_k\sum_{j\le k}|c_k|\binom kjs^j|b-a|^{k-j}=\sum_k|c_k|(s+|b-a|)^k<\infty$, así que la serie doble que se obtiene
desarrollando $(t-a)^k=\big((t-b)+(b-a)\big)^k$ por el binomio converge absolutamente para $|t-b|\le s$ y puede
sumarse en cualquier orden; sumando primero sobre $k$ se obtienen los $d_j$, con $\sum_j|d_j|s^j<\infty$.
(b) Sea $j_0$ el menor $j$ con $d_j\ne0$. Entonces $g(t)=(t-b)^{j_0}h(t)$ con
$h(t)=\sum_{j\ge j_0}d_j(t-b)^{j-j_0}$, que es continua cerca de $b$ (criterio M) con $h(b)=d_{j_0}\neq0$; luego
$h\neq0$ cerca de $b$. (c) Sea $A$ el conjunto de los $b\in I$ cerca de los cuales $g$ se anula idénticamente; $A$ es abierto.
Si los ceros se acumulan en $b\in I$, entonces por (b) todos los $d_j$ se anulan y $b\in A$ por (a). Si $b\in I$ es límite de
puntos de $A$, entonces los ceros se acumulan en $b$, luego $b\in A$: $A$ es cerrado en $I$. Como $I$ es conexo, $A$ es
vacío o igual a $I$. Si $g\not\equiv0$, entonces $A=\emptyset$, así que ningún punto de $I$ es punto de acumulación de ceros;
todo intervalo compacto $J\subset I$ contiene entonces un número finito de ceros (un subconjunto infinito de $J$ tiene un
punto de acumulación en $J$ por Bolzano–Weierstrass), e $I$ es unión numerable de tales $J$. Si $g\equiv0$,
(b) en $b=a$ (donde $d_j=c_j$) muestra que todos los $c_k=0$. $\square$

**Teorema A.6 (los conjuntos de ceros son nulos).** Sea $f$ una serie de potencias sobre un cubo $Q\subseteq\mathbb R^n$. Si $f$
se anula sobre un subconjunto de $Q$ de medida de Lebesgue positiva, entonces $f=0$ sobre $Q$. De forma equivalente: si
$f\not\equiv0$ sobre $Q$, entonces $\{x\in Q:f(x)=0\}$ tiene medida de Lebesgue nula.

*Demostración.* Inducción sobre $n$; el caso $n=1$ es el lema A.5(c). Sea $n\ge2$; escribamos $x=(x',t)$ con
$x'\in\mathbb R^{n-1}$, $Q=Q'\times I$ y $a=(a',a_n)$, y agrupemos los términos según la potencia de $t-a_n$:
$$
f(x',t)=\sum_kh_k(x')(t-a_n)^k,\qquad h_k(x')=\sum_{\alpha'}c_{(\alpha',k)}(x'-a')^{\alpha'} .
$$
Esto es lícito: para $x'\in Q'$ y $\max_i|x'_i-a'_i|\le r<\rho$,
$\sum_k\sum_{\alpha'}|c_{(\alpha',k)}||x'-a'|^{\alpha'}r^k\le\sum_\alpha|c_\alpha|r^{|\alpha|}<\infty$ (con
$|y|^{\alpha'}=\prod_i|y_i|^{\alpha'_i}$). La misma cota muestra que cada $h_k$ es una serie de potencias sobre $Q'$ y que, para
$x'\in Q'$ fijo, $g_{x'}(t)=f(x',t)$ es una serie de potencias sobre $I$ con coeficientes $h_k(x')$.

Sea $Z=\{x\in Q:f(x)=0\}$, cerrado en $Q$ porque $f$ es continua, y supongamos $\lambda_n(Z)>0$. Por el teorema de
Tonelli, $\lambda_n(Z)=\int_{Q'}\lambda_1(Z_{x'})\,dx'$ con $Z_{x'}=\{t\in I:f(x',t)=0\}$, así que el conjunto
$F=\{x'\in Q':\lambda_1(Z_{x'})>0\}$ tiene medida $(n-1)$-dimensional positiva. Para $x'\in F$, el lema A.5(c) da
$g_{x'}\equiv0$ sobre $I$ y luego $h_k(x')=0$ para todo $k$. Así, cada $h_k$ se anula sobre $F$ y, por la hipótesis de
inducción, $h_k=0$ sobre $Q'$. Por tanto $f(x',t)=\sum_kh_k(x')(t-a_n)^k=0$ sobre $Q$. $\square$

**Corolario A.7.** (a) Un polinomio real $p$ sobre $\mathbb R^n$ que no es idénticamente nulo solo se anula sobre un conjunto
de medida de Lebesgue nula. Un polinomio con algún coeficiente no nulo no es idénticamente nulo.

(b) Sea $W$ un espacio vectorial real de dimensión finita y $p:W\to\mathbb R$ un polinomio en las coordenadas de
alguna base (entonces lo es en todas las bases), no idénticamente nulo. Si un vector aleatorio $X$ en $W$ tiene densidad,
entonces $\mathbb P(p(X)=0)=0$.

(c) Un subespacio lineal propio $U\subsetneq W$ es un conjunto nulo, y un vector aleatorio con densidad está fuera de
$U$ casi seguramente. En particular, si $E\neq0$ es un subespacio de $\mathbb C^n$ y $v$ tiene densidad sobre $\mathbb C^n$
(por ejemplo, un vector inicial gaussiano), entonces $v$ tiene una componente no nula en $E$ casi seguramente.

(d) Sean $X$ e $Y$ vectores aleatorios independientes en espacios reales de dimensión finita $W_1$ y $W_2$, $Y$ con
densidad, y $N\subseteq W_1\times W_2$ un boreliano cuyas secciones $N_x=\{y:(x,y)\in N\}$ son conjuntos nulos para todo
$x$. Entonces $\mathbb P((X,Y)\in N)=0$.

*Demostración.* (a) La primera afirmación es el teorema A.6 con $Q=\mathbb R^n$. Para la segunda, inducción sobre $n$: se escribe
$p(x',t)=\sum_kp_k(x')t^k$; algún $p_k$ tiene un coeficiente no nulo, así que $p_k(x'_0)\neq0$ para algún $x'_0$ por
inducción (para $n=1$ los $p_k$ son constantes y no hay $x'$); entonces $t\mapsto p(x'_0,t)$ es un polinomio no nulo en una variable, que tiene un número finito de raíces.
(b) Un cambio de base es lineal, así que $p$ sigue siendo un polinomio, y $\{p=0\}$ es nulo por (a) y el lema A.3(a);
entonces $\mathbb P(p(X)=0)=\int_{\{p=0\}}f\,d\lambda=0$. (c) Se elige una forma lineal no nula $\ell$ que se anule sobre $U$
y se aplica (b) a $p=\ell$. Para el ejemplo, $v\mapsto P_Ev$ ($P_E$ la proyección ortogonal) es una aplicación
$\mathbb R$-lineal no nula, y $\{P_Ev=0\}$ es su núcleo, un subespacio propio. (d) La ley de $(X,Y)$ es el producto de las leyes, así que por
el teorema de Tonelli $\mathbb P((X,Y)\in N)=\int\mathbb P(Y\in N_x)\,d\mathbb P_X(x)=0$. $\square$

El curso aplica el corolario A.7 a $\det Z$ para una matriz de Ginibre $Z$ ($|\det Z|^2$ es un polinomio en las partes real
e imaginaria de las entradas, igual a $1$ en $Z=1$; teorema 5.3, sección 1.6), al determinante de Gram de la
proposición 4.7(a), al rango de matrices gaussianas (capítulo 12, ejercicio 4) y a los vectores iniciales del
método de Lanczos (capítulo 9). La parte (d) es el paso «como $Y$ es independiente de $K$, el teorema de Tonelli da
probabilidad $0$» de la proposición 4.7(b).

```python
import math
import numpy as np

rng = np.random.default_rng(0)


def mc(samples, exact, label, kind="exact"):
    """Monte Carlo mean with standard error; the reference value must lie within 5 standard errors."""
    samples = np.asarray(samples, dtype=float)
    mean, err = samples.mean(), samples.std() / math.sqrt(samples.size)
    print(f"{label}: {mean:.4f} +- {err:.4f}   ({kind} {exact:.4f})")
    assert abs(mean - exact) <= 5 * err + 1e-12, label


# Corollary A.7: det Z = 0 has probability 0 for a real Gaussian 3 x 3 matrix, so by continuity of the measure
# P(|det Z| < eps) -> 0 as eps -> 0. Numerical observation (not a theorem of this appendix): the fractions
# below shrink by about a factor 10 per factor 10 in eps.
d = np.abs(np.linalg.det(rng.standard_normal((200_000, 3, 3))))
fractions = [(d < eps).mean() for eps in (1e-1, 1e-2, 1e-3)]
print("P(|det Z| < eps) for eps = 0.1, 0.01, 0.001:", fractions)
assert fractions[1] < 0.2 * fractions[0] and fractions[2] < 0.2 * fractions[1]
```

## A.3 Valores propios de matrices hermíticas

Para $A\in\operatorname{Herm}_n$ sean $\lambda_1(A)\le\dots\le\lambda_n(A)$ sus valores propios contados con multiplicidad, y
$v_1,\dots,v_n$ una base ortonormal de vectores propios (teorema 1.3, corolario 1.4). Para un vector unitario $x$,
$\langle x,Ax\rangle=\sum_i\lambda_i|\langle v_i,x\rangle|^2$ con $\sum_i|\langle v_i,x\rangle|^2=1$, luego
$\lambda_1(A)\le\langle x,Ax\rangle\le\lambda_n(A)$; y $\|A\|=\max_i|\lambda_i(A)|$ (corolario 1.4).

**Teorema A.8 (Courant–Fischer).** Para $A\in\operatorname{Herm}_n$ y $1\le j\le n$,
$$
\lambda_j(A)=\min_{\dim U=j}\ \max_{x\in U,\,|x|=1}\langle x,Ax\rangle=\max_{\dim U=n-j+1}\ \min_{x\in U,\,|x|=1}\langle x,Ax\rangle ,
$$
sobre los subespacios $U\subseteq\mathbb C^n$ de la dimensión indicada. (Los máximos y mínimos interiores existen porque la
esfera unidad de $U$ es compacta; los exteriores se alcanzan.)

*Demostración.* Sean $\dim U=j$ y $W=\operatorname{span}(v_j,\dots,v_n)$. Como $\dim U+\dim W=n+1>n$, existe un vector unitario
$x\in U\cap W$, y $\langle x,Ax\rangle=\sum_{i\ge j}\lambda_i|\langle v_i,x\rangle|^2\ge\lambda_j$. Así, el máximo sobre $U$
es $\ge\lambda_j$. Para $U=\operatorname{span}(v_1,\dots,v_j)$, $\langle x,Ax\rangle=\sum_{i\le j}\lambda_i|\langle v_i,x\rangle|^2\le\lambda_j$,
con igualdad en $x=v_j$; luego el mínimo sobre $U$ se alcanza y vale $\lambda_j$. Para la segunda fórmula
se aplica la primera a $-A$, cuyos valores propios son $\lambda_j(-A)=-\lambda_{n+1-j}(A)$:
$-\lambda_j(A)=\lambda_{n+1-j}(-A)=\min_{\dim U=n+1-j}\max_x(-\langle x,Ax\rangle)=-\max_{\dim U=n+1-j}\min_x\langle x,Ax\rangle$. $\square$

**Corolario A.9 (compresiones).** Sean $Q\in M_{n\times m}(\mathbb C)$ con $Q^*Q=1_m$ y $A\in\operatorname{Herm}_n$. Entonces
$\lambda_j(Q^*AQ)=\min\max\langle x,Ax\rangle$ sobre los $U\subseteq\operatorname{ran}Q$ de dimensión $j$ y los $x\in U$ unitarios, y
$$
\lambda_j(A)\le\lambda_j(Q^*AQ)\le\lambda_{j+n-m}(A)\qquad(1\le j\le m).
$$

*Demostración.* $U'\mapsto QU'$ es una biyección de los subespacios de dimensión $j$ de $\mathbb C^m$ sobre los de
$\operatorname{ran}Q$, $|Qy|=|y|$ y $\langle y,Q^*AQy\rangle=\langle Qy,AQy\rangle$; así, el teorema A.8 para $Q^*AQ$ da la
fórmula. Minimizando sobre la familia más amplia de todos los $U\subseteq\mathbb C^n$ de dimensión $j$ se obtiene $\lambda_j(A)$, que es
por tanto menor. Del mismo modo, la fórmula max–min de $\lambda_j(Q^*AQ)$ recorre los $U\subseteq\operatorname{ran}Q$ de
dimensión $m-j+1$, y sobre todos los $U\subseteq\mathbb C^n$ de esa dimensión da el número mayor
$\lambda_{n-(m-j+1)+1}(A)=\lambda_{j+n-m}(A)$. $\square$

**Teorema A.10 (Weyl).** Para $A,E\in\operatorname{Herm}_n$ y todo $j$,
$$
\lambda_j(A)+\lambda_1(E)\le\lambda_j(A+E)\le\lambda_j(A)+\lambda_n(E),
\qquad\text{hence}\qquad
|\lambda_j(A+E)-\lambda_j(A)|\le\|E\|\le\|E\|_2 .
$$

*Demostración.* Para $x$ unitario, $\langle x,(A+E)x\rangle\le\langle x,Ax\rangle+\lambda_n(E)$. Tomando el máximo sobre los vectores
unitarios de un $U$ de dimensión $j$ y luego el mínimo sobre $U$ (teorema A.8): $\lambda_j(A+E)\le\lambda_j(A)+\lambda_n(E)$. Aplicado a
$A+E$ y $-E$, esto da $\lambda_j(A)\le\lambda_j(A+E)+\lambda_n(-E)=\lambda_j(A+E)-\lambda_1(E)$. Por tanto
$|\lambda_j(A+E)-\lambda_j(A)|\le\max(|\lambda_1(E)|,|\lambda_n(E)|)=\|E\|$, y $\|E\|^2=\max_i\lambda_i(E)^2\le\sum_i\lambda_i(E)^2=\|E\|_2^2$. $\square$

**Corolario A.11.** Cada $\lambda_j:\operatorname{Herm}_n\to\mathbb R$ es $1$-lipschitziana para $\|\cdot\|$, luego continua y
boreliana. Para una matriz hermítica aleatoria $H$, $(\lambda_1(H),\dots,\lambda_n(H))$ es un vector aleatorio, y también lo son
sus espaciamientos y sus cocientes de espaciamientos (sobre el boreliano donde sus denominadores no se anulan); los vectores de valores propios de matrices hermíticas aleatorias independientes son independientes. Las matrices con
$n$ valores propios distintos forman un abierto.

*Demostración.* Teorema A.10 con $E=B-A$. Las funciones de variables aleatorias independientes son independientes. Si
$\lambda_j(A)<\lambda_{j+1}(A)$ para todo $j$, las desigualdades estrictas persisten para $\|B-A\|<\frac12\min_j(\lambda_{j+1}(A)-\lambda_j(A))$. $\square$

**Lema A.12 (fórmula de Lagrange).** Sea $A\in\operatorname{Herm}_n$ con valores propios distintos $\mu_1<\dots<\mu_r$ y
proyecciones espectrales $P_1,\dots,P_r$ (proyecciones ortogonales sobre los subespacios propios, de modo que $A=\sum_i\mu_iP_i$,
$P_iP_l=\delta_{il}P_i$, $\sum_iP_i=1$). Entonces $p(A)=\sum_ip(\mu_i)P_i$ para todo polinomio $p$, y
$$
P_i=\prod_{l\neq i}\frac{A-\mu_l}{\mu_i-\mu_l} .
$$

*Demostración.* $A^k=\sum_i\mu_i^kP_i$ por inducción sobre $k$ (usando $P_iP_l=\delta_{il}P_i$), y $p(A)=\sum_ip(\mu_i)P_i$ por linealidad.
El polinomio $\ell_i(t)=\prod_{l\neq i}(t-\mu_l)/(\mu_i-\mu_l)$ cumple $\ell_i(\mu_l)=\delta_{il}$, luego $\ell_i(A)=P_i$. $\square$

**Proposición A.13 (continuidad de los valores propios y de las proyecciones espectrales).** Sea $S$ un subconjunto de un espacio
vectorial real de dimensión finita, $A:S\to\operatorname{Herm}_n$ continua, y supongamos que $A(x)$ tiene exactamente $r$ valores propios
distintos $\mu_1(x)<\dots<\mu_r(x)$ para todo $x\in S$. Entonces sus multiplicidades son localmente constantes sobre $S$, y las $\mu_i$ y
las proyecciones espectrales $P_i(x)$ de $A(x)$ son continuas sobre $S$. En particular, sobre el abierto de las matrices
hermíticas de espectro simple, los valores propios y las proyecciones espectrales de rango uno son continuos.

*Demostración.* Si $r=1$, $A(x)=\mu_1(x)1$ con $\mu_1=\operatorname{Tr}A/n$ y $P_1=1$. Sea $r\ge2$; fijemos $x_0\in S$ y sea
$d=\min_i(\mu_{i+1}(x_0)-\mu_i(x_0))>0$. Sea $J_i$ el conjunto de los índices $j$ con $\lambda_j(A(x_0))=\mu_i(x_0)$, un bloque de
enteros consecutivos. Por la continuidad de $A$ existe un entorno $N$ de $x_0$ en $S$ con $\|A(x)-A(x_0)\|<d/2$ para
$x\in N$; por el teorema A.10, $\lambda_j(A(x))$ está en el intervalo $I_i=(\mu_i(x_0)-d/2,\mu_i(x_0)+d/2)$ para $j\in J_i$, y
estos $r$ intervalos son disjuntos. Así, $A(x)$ tiene al menos un valor propio distinto en cada $I_i$ y, como tiene exactamente $r$,
exactamente uno: $\lambda_j(A(x))$ es el mismo para todos los $j\in J_i$. Por tanto, sobre $N$ las multiplicidades son las constantes $|J_i|$ y
$\mu_i(x)=\lambda_{\max J_i}(A(x))$ es continua (corolario A.11). Por el lema A.12,
$P_i(x)=\prod_{l\neq i}(A(x)-\mu_l(x))/(\mu_i(x)-\mu_l(x))$, continua sobre $N$ porque los denominadores se mantienen alejados de $0$.
La última afirmación es el caso $r=n$, con $S$ el abierto del corolario A.11. $\square$

En la proposición 4.7, $S$ es el conjunto donde el determinante de Gram $\delta$ no se anula, sobre el cual $K$ tiene exactamente $N$
valores propios distintos; así, $K\mapsto\Pi_i(K)$ es continua allí y el suceso de la proposición 4.7(b) es boreliano. El teorema 5.15
utiliza el corolario A.11, y el teorema 9.14 utiliza el corolario A.9.

```python
def herm(n, batch):
    G = rng.standard_normal((batch, n, n)) + 1j * rng.standard_normal((batch, n, n))
    return (G + G.conj().transpose(0, 2, 1)) / 2


n = 6
A, E = herm(n, 2000), 0.1 * herm(n, 2000)
lam_A, lam_AE, lam_E = np.linalg.eigvalsh(A), np.linalg.eigvalsh(A + E), np.linalg.eigvalsh(E)   # ascending
normE = np.linalg.norm(E, ord=2, axis=(1, 2))                                    # operator norm
assert (np.abs(lam_AE - lam_A).max(axis=1) <= normE + 1e-12).all()               # Theorem A.10
assert (lam_AE >= lam_A + lam_E[:, :1] - 1e-12).all() and (lam_AE <= lam_A + lam_E[:, -1:] + 1e-12).all()

a0, j = A[0], 3                          # Theorem A.8: every 3-dim subspace has max <x, Ax> >= lambda_3
for _ in range(200):
    Q, _ = np.linalg.qr(rng.standard_normal((n, j)) + 1j * rng.standard_normal((n, j)))
    assert np.linalg.eigvalsh(Q.conj().T @ a0 @ Q)[-1] >= lam_A[0, j - 1] - 1e-12
m = 4                                    # Corollary A.9: interlacing for a compression to 4 dimensions
Q, _ = np.linalg.qr(rng.standard_normal((n, m)) + 1j * rng.standard_normal((n, m)))
theta = np.linalg.eigvalsh(Q.conj().T @ a0 @ Q)
assert (theta >= lam_A[0, :m] - 1e-12).all() and (theta <= lam_A[0, n - m:] + 1e-12).all()

U, _ = np.linalg.qr(rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)))
mu, mult = np.array([-1.0, 0.5, 2.0]), [2, 3, 1]          # Lemma A.12 with multiplicities 2, 3, 1
spec = np.repeat(mu, mult)
A1 = U @ np.diag(spec) @ U.conj().T
for i, m_i in enumerate(mult):
    P = np.eye(n, dtype=complex)
    for l in range(len(mu)):
        if l != i:
            P = P @ (A1 - mu[l] * np.eye(n)) / (mu[i] - mu[l])
    Vi = U[:, spec == mu[i]]
    assert np.abs(P - Vi @ Vi.conj().T).max() < 1e-10 and abs(np.trace(P).real - m_i) < 1e-10
```

## A.4 Vectores gaussianos

Un **vector gaussiano estándar** en $\mathbb R^m$ es $g=(g_1,\dots,g_m)$ con coordenadas $N(0,1)$ independientes. Su ley
es el producto de las leyes de las coordenadas, así que, por el teorema de Tonelli, tiene la densidad
$\varphi_m(x)=(2\pi)^{-m/2}e^{-|x|^2/2}$.

**Proposición A.14 (invariancia ortogonal).** (a) $Og$ tiene la ley de $g$ para toda $O\in O(m)$ ortogonal.

(b) Sea $W$ un espacio real con producto escalar de dimensión $m$. Para una base ortonormal $(w_i)$, la ley de
$\sum_ig_iw_i$ no depende de la base; un tal vector se llama **vector gaussiano estándar de $W$**. Sus coordenadas
en toda base ortonormal son $N(0,1)$ independientes, su ley es invariante por las aplicaciones ortogonales de $W$, y si
$W=W_1\oplus W_2$ ortogonalmente, sus componentes son vectores gaussianos estándar independientes de $W_1$ y de $W_2$.

(c) Para $v\in\mathbb R^m$, $\langle v,g\rangle$ tiene la ley $N(0,|v|^2)$; en particular $\sum_ia_ig_i\sim N(0,\sum_ia_i^2)$.

*Demostración.* (a) Por el lema A.3(c), $Og$ tiene la densidad $\varphi_m(O^{-1}y)/|\det O|=\varphi_m(y)$. (b) Si $w'_i=\sum_jO_{ji}w_j$
es otra base ortonormal, $O$ es ortogonal y $\sum_ig_iw'_i=\sum_j(Og)_jw_j$, que tiene la ley de
$\sum_jg_jw_j$ por (a). Las coordenadas de $X=\sum_ig_iw_i$ en otra base ortonormal $(w'_i)$ son $O^Tg$, de nuevo
estándar por (a). Una aplicación ortogonal $L$ transforma $(w_i)$ en una base ortonormal, y $LX=\sum_ig_iLw_i$. Para la última
afirmación se toma $(w_i)$ adaptada a $W_1\oplus W_2$. (c) Para $v=0$ ambos lados son la masa puntual en $0$. Para $v\neq0$
se completa $v/|v|$ hasta una base ortonormal; por (b) la coordenada $\langle v/|v|,g\rangle$ es $N(0,1)$, luego
$\langle v,g\rangle=|v|\,\langle v/|v|,g\rangle$ tiene la densidad $(2\pi|v|^2)^{-1/2}e^{-y^2/2|v|^2}$ (lema A.3(c) con $m=1$). $\square$

**Definición.** Un **vector gaussiano centrado** en un espacio vectorial real de dimensión finita $W$ es $X=Ag$ para una aplicación
lineal $A:\mathbb R^M\to W$ y un vector gaussiano estándar $g$ de $\mathbb R^M$. En coordenadas ($W=\mathbb R^m$), su
**covarianza** es $C=\mathbb E[XX^T]=AA^T$. Es **no degenerado** si $A$ es sobreyectiva, o, de forma equivalente, si $C$ es definida
positiva ($\langle x,Cx\rangle=|A^Tx|^2$).

Todo vector aleatorio del curso construido linealmente a partir de entradas gaussianas independientes —vectores gaussianos
complejos, matrices de Ginibre, matrices del GOE y del GUE, sus imágenes por twirls— es un vector gaussiano centrado.

**Teorema A.15 (imágenes lineales y densidades).** Sea $X=Ag$ un vector gaussiano centrado en $W$.

(a) Para toda aplicación lineal $L:W\to W'$, $LX$ es un vector gaussiano centrado; es no degenerado si $X$ lo es y $L$ es sobreyectiva.

(b) Si $X$ es no degenerado, tiene densidad respecto de la medida de Lebesgue sobre $W$; en coordenadas
$W=\mathbb R^m$,
$$
f_X(x)=(2\pi)^{-m/2}(\det C)^{-1/2}\,e^{-\langle x,C^{-1}x\rangle/2}.
$$
Si $X$ es degenerado, está en el subespacio propio $\operatorname{ran}A$, un conjunto nulo (corolario A.7(c)), y no tiene densidad.

(c) La ley de $X$ solo depende de su covarianza.

*Demostración.* (a) $LX=(LA)g$, y $LA$ es sobreyectiva si $A$ y $L$ lo son. (b) Sea $K=\ker A$, de modo que
$\mathbb R^M=K^\perp\oplus K$. Por la proposición A.14(b), $g=g_1+g_2$ con $g_1$ un vector gaussiano estándar de
$K^\perp$, y $X=Ag_1=Bg_1$ con $B=A|_{K^\perp}$, que es inyectiva y, como $A$ es sobreyectiva, biyectiva
$K^\perp\to W$; luego $\dim K^\perp=m$. En una base ortonormal de $K^\perp$ completada con una de $K$, la matriz de $A$ es
$(B\ 0)$, luego $C=AA^T=BB^T$. Las coordenadas de $g_1$ tienen la densidad $\varphi_m$ y, por el lema A.3(c), $X$ tiene
la densidad $\varphi_m(B^{-1}x)/|\det B|$, con $|B^{-1}x|^2=\langle x,(BB^T)^{-1}x\rangle$ y $|\det B|=(\det C)^{1/2}$. Si
$X$ fuera degenerado y tuviera una densidad $f$, entonces $1=\mathbb P(X\in\operatorname{ran}A)=\int_{\operatorname{ran}A}f\,d\lambda=0$.
(c) $\ker C=\ker A^T$ (si $AA^Tx=0$, entonces $|A^Tx|^2=0$), así que $\operatorname{ran}C=(\ker C)^\perp=(\ker A^T)^\perp=\operatorname{ran}A$
contiene a $X$. Sean $d=\operatorname{rank}C$ (si $d=0$, entonces $X=0$), $p_1,\dots,p_d$ una base ortonormal de
$\operatorname{ran}C$ (que solo depende de $C$) y $P$ la matriz $d\times m$ de filas $p_k^T$, de modo que $PP^T=1_d$ y
$P^TP$ es la proyección ortogonal sobre $\operatorname{ran}C$. Entonces $PX=(PA)g$ con $PA$ sobreyectiva sobre $\mathbb R^d$, un vector
gaussiano centrado no degenerado de covarianza $PCP^T$; por (b) su ley solo depende de $C$, y lo mismo ocurre con la ley de
$X=P^T(PX)$. $\square$

Por ejemplo, en el capítulo 4 el twirl $P$ es una proyección lineal de $M_n(\mathbb C)$ sobre $\mathcal A'$, así que
$Z\mapsto P(Z)$ y $Z\mapsto P(Z+Z^*)$ son sobreyecciones $\mathbb R$-lineales sobre $\mathcal A'$ y sobre $\mathcal A'_h$ (para
$H\in\mathcal A'_h$, $P(\frac H2+\frac H2)=H$); para una matriz de Ginibre $Z$ (no degenerada por la proposición A.16), $Y=P(Z)$
y $K=P(Z+Z^*)$ tienen densidades sobre estos espacios por el teorema A.15. Una combinación gaussiana $\sum_k\xi_kB_k$ de una base
$(B_k)$ de $\mathcal A'$ es la imagen del vector gaussiano estándar $\xi$ por una biyección lineal.

**Proposición A.16 (vectores gaussianos complejos, matrices de Ginibre).** (a) Sea $z\in\mathbb C^n$ con coordenadas
independientes cuyas partes real e imaginaria son $N(0,\frac12)$ independientes (un **vector gaussiano complejo estándar**).
Entonces $\sqrt2z$ es el vector gaussiano estándar de $\mathbb C^n\cong\mathbb R^{2n}$ con el producto escalar
$\operatorname{Re}\langle u,w\rangle$; $z$ tiene la densidad $\pi^{-n}e^{-|z|^2}$; $Uz$ tiene la ley de $z$ para toda
$U\in U(n)$, así que las coordenadas $\langle w_i,z\rangle$ en cualquier base ortonormal son de nuevo gaussianas complejas
estándar independientes; y $\mathbb E[z_i\bar z_j]=\delta_{ij}$, $\mathbb E[z_iz_j]=0$.

(b) Una matriz de Ginibre $G$ de tamaño $a\times b$ (entradas gaussianas complejas estándar independientes) tiene la densidad
$\pi^{-ab}e^{-\|G\|_2^2}$ sobre $M_{a\times b}(\mathbb C)$, y $VGW$ tiene la ley de $G$ para todas $V\in U(a)$, $W\in U(b)$.

(c) Una matriz real $G$ de tamaño $a\times b$ con entradas $N(0,1)$ independientes tiene la densidad $(2\pi)^{-ab/2}e^{-\|G\|_2^2/2}$, y
$OGO'$ tiene la ley de $G$ para todas $O\in O(a)$, $O'\in O(b)$.

*Demostración.* (a) Los $2n$ números $\sqrt2\operatorname{Re}z_i,\sqrt2\operatorname{Im}z_i$ son $N(0,1)$ independientes, y son las
coordenadas de $\sqrt2z$ en la base ortonormal $(e_i,ie_i)$. Por el lema A.3(c) con $L=2^{-1/2}$, $z$ tiene la densidad
$(2\pi)^{-n}e^{-|\sqrt2z|^2/2}2^{n}=\pi^{-n}e^{-|z|^2}$. Una unitaria $U$ es $\mathbb R$-lineal y conserva
$\operatorname{Re}\langle u,w\rangle$, luego es ortogonal; se aplica la proposición A.14(b). Las coordenadas en la base $(w_i)$ son
$W^*z$, con $W$ la unitaria de columnas $w_i$. Finalmente, con $z_j=x_j+iy_j$,
$\mathbb E[z_i\bar z_j]=\mathbb E[x_ix_j+y_iy_j]+i\,\mathbb E[y_ix_j-x_iy_j]=\delta_{ij}$ y $\mathbb E[z_iz_j]=\mathbb E[x_ix_j-y_iy_j]+i\,\mathbb E[x_iy_j+y_ix_j]=0$.
(b) $M_{a\times b}(\mathbb C)\cong\mathbb C^{ab}$ con $|G|=\|G\|_2$, y $X\mapsto VXW$ es unitaria para
$\operatorname{Tr}(X^*Y)$ (convenciones tras el lema A.3); se aplica (a). (c) Igual, con $\mathbb R^{ab}$ y la proposición A.14. $\square$

**Proposición A.17 (GOE y GUE).** Sea $H=(G+G^T)/2$ con $G$ una matriz real $n\times n$ de entradas $N(0,1)$ independientes
(GOE), o $H=(G+G^*)/\sqrt2$ con $G$ una matriz de Ginibre (GUE). Entonces las entradas $h_{ij}$, $i\le j$, son independientes,
$h_{ii}\sim N(0,1)$ y, para $i<j$, $h_{ij}\sim N(0,\frac12)$ (GOE), o bien $\operatorname{Re}h_{ij},\operatorname{Im}h_{ij}$ son
$N(0,\frac12)$ independientes (GUE). La densidad de $H$ en las coordenadas fijadas tras el lema A.3 es
$$
(2\pi)^{-n/2}\pi^{-n(n-1)/4}e^{-\operatorname{Tr}H^2/2}\ \ \text{(GOE)},\qquad
(2\pi)^{-n/2}\pi^{-n(n-1)/2}e^{-\operatorname{Tr}H^2/2}\ \ \text{(GUE)},
$$
y $OHO^T$ ($O\in O(n)$), respectivamente $UHU^*$ ($U\in U(n)$), tiene la ley de $H$.

*Demostración.* GOE: $G$ es el vector gaussiano estándar de $M_n(\mathbb R)$ con $\operatorname{Tr}(X^TY)$, y las matrices
$e_{ii}$, $(e_{ij}\pm e_{ji})/\sqrt2$ ($i<j$) forman una base ortonormal. Por la proposición A.14(b), las coordenadas
$G_{ii}$, $(G_{ij}+G_{ji})/\sqrt2$, $(G_{ij}-G_{ji})/\sqrt2$ son $N(0,1)$ independientes. Ahora bien, $h_{ii}=G_{ii}$ y
$h_{ij}=2^{-1/2}\cdot(G_{ij}+G_{ji})/\sqrt2$. GUE: $\sqrt2G$ es el vector gaussiano estándar de $M_n(\mathbb C)$ con
$\operatorname{Re}\operatorname{Tr}(X^*Y)$ (proposición A.16), con la base ortonormal $e_{ii}$, $ie_{ii}$,
$(e_{ij}\pm e_{ji})/\sqrt2$, $i(e_{ij}\pm e_{ji})/\sqrt2$; las coordenadas de $\sqrt2G$ son $\sqrt2\operatorname{Re}G_{ii}$,
$\sqrt2\operatorname{Im}G_{ii}$, $\operatorname{Re}(G_{ij}\pm G_{ji})$, $\operatorname{Im}(G_{ij}\pm G_{ji})$, $N(0,1)$ independientes. Ahora bien,
$h_{ii}=\sqrt2\operatorname{Re}G_{ii}$, $\operatorname{Re}h_{ij}=2^{-1/2}\operatorname{Re}(G_{ij}+G_{ji})$ y
$\operatorname{Im}h_{ij}=2^{-1/2}\operatorname{Im}(G_{ij}-G_{ji})$. En ambos casos la densidad es el producto de las
densidades unidimensionales, $(2\pi)^{-1/2}e^{-h_{ii}^2/2}$ y $\pi^{-1/2}e^{-t^2}$ para cada coordenada real $t$ de varianza
$\frac12$, y $\operatorname{Tr}H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2$. La aplicación $H\mapsto UHU^*$ es lineal, invertible y
ortogonal para $\operatorname{Tr}(XY)$, así que, por el lema A.3(b),(c), $UHU^*$ tiene la densidad $f(U^*HU)=f(H)$, ya que $f$ solo depende de
$\operatorname{Tr}H^2$; lo mismo para $O$. $\square$

```python
Amat = rng.standard_normal((2, 4))                       # Theorem A.15: X = A g in R^2, A of rank 2
C = Amat @ Amat.T
X = rng.standard_normal((400_000, 4)) @ Amat.T
h = 1 / 400                                              # midpoint rule for the density on [0, 1]^2
grid = (np.arange(400) + 0.5) * h
xx = np.stack(np.meshgrid(grid, grid, indexing='ij'), -1).reshape(-1, 2)
dens = np.exp(-0.5 * np.einsum('ni,ij,nj->n', xx, np.linalg.inv(C), xx)) / (2 * math.pi * math.sqrt(np.linalg.det(C)))
mc(((X >= 0) & (X <= 1)).all(axis=1), dens.sum() * h * h, "P(X in [0,1]^2) against the integral of f_X",
   kind="midpoint rule")

nn, B = 4, 100_000                                       # Proposition A.17
G = rng.standard_normal((B, nn, nn))
goe = (G + G.transpose(0, 2, 1)) / 2
Gc = (rng.standard_normal((B, nn, nn)) + 1j * rng.standard_normal((B, nn, nn))) / math.sqrt(2)
gue = (Gc + Gc.conj().transpose(0, 2, 1)) / math.sqrt(2)
mc(np.einsum('bij,bji->b', goe, goe), nn * (nn + 1) / 2, "E Tr H^2, GOE")      # n + 2 * n(n-1)/2 * 1/2
mc(np.einsum('bij,bji->b', gue, gue).real, nn ** 2, "E Tr H^2, GUE")           # n + 2 * n(n-1)/2 * 1
O, _ = np.linalg.qr(rng.standard_normal((nn, nn)))
rot = O @ goe @ O.T                                                             # has the law of the GOE
mc(rot[:, 0, 0] ** 4, 3.0, "E (O H O^T)_11^4")
mc(rot[:, 0, 1] ** 2, 0.5, "E (O H O^T)_12^2")
```

## A.5 La dirección de un vector gaussiano

Sean $S^{n-1}=\{x\in\mathbb R^n:|x|=1\}$, $V_n=\lambda_n(\{|x|\le1\})$ y, para un boreliano $B\subseteq S^{n-1}$,
$C(B)=\{x\in\mathbb R^n:0<|x|\le1,\ x/|x|\in B\}$, que es boreliano ($x\mapsto(|x|,x/|x|)$ es continua sobre $\mathbb R^n\setminus0$).

**Teorema A.18 (coordenadas polares).** Para borelianos $A\subseteq(0,\infty)$ y $B\subseteq S^{n-1}$,
$$
\lambda_n\big(\{x\neq0:|x|\in A,\ x/|x|\in B\}\big)=n\,\lambda_n(C(B))\int_Ar^{n-1}\,dr .
$$
En consecuencia, para $h:(0,\infty)\to[0,\infty]$ y $\phi:S^{n-1}\to[0,\infty]$ borelianas, con $\mu(B)=\lambda_n(C(B))$,
$$
\int_{\mathbb R^n}h(|x|)\,\phi(x/|x|)\,dx=n\int_0^\infty h(r)r^{n-1}\,dr\int_{S^{n-1}}\phi\,d\mu .
$$

*Demostración.* Fijemos $B$ y sean $\nu_1(A)$ y $\nu_2(A)$ los dos lados, medidas de Borel sobre $(0,\infty)$. La dilatación
$x\mapsto rx$ tiene determinante $r^n$ y transforma $C(B)$ en $\{0<|x|\le r,\ x/|x|\in B\}$, pues $rx/|rx|=x/|x|$; así,
$\nu_1((0,r])=r^n\lambda_n(C(B))=\nu_2((0,r])$. Para $R>0$ las restricciones de $\nu_1,\nu_2$ a $(0,R]$ son finitas
(a lo sumo $\lambda_n(\{|x|\le R\})$) y coinciden sobre los intervalos $(0,r]$, $r\le R$, que son estables por intersección,
contienen $(0,R]$ y generan los borelianos de $(0,R]$ (dan $(s,r]=(0,r]\setminus(0,s]$, y todo abierto relativo
de $(0,R]$ es unión numerable de tales intervalos). Por el lema A.1, $\nu_1=\nu_2$ sobre $(0,R]$; se hace tender $R\to\infty$. La
fórmula integral vale para $h=1_A$, $\phi=1_B$ (la identidad que se acaba de demostrar; $\{0\}$ es nulo), luego para $h,\phi$ simples
no negativas por linealidad, y para $h,\phi$ generales por convergencia monótona a lo largo de funciones simples $h_k\uparrow h$, $\phi_k\uparrow\phi$. $\square$

**Corolario A.19 (integrales radiales).** Para $h\ge0$ boreliana sobre $(0,\infty)$,
$\int_{\mathbb R^n}h(|x|)\,dx=s_n\int_0^\infty h(r)r^{n-1}\,dr$ con $s_n=nV_n=(2\pi)^{n/2}/I_n$,
$I_n=\int_0^\infty r^{n-1}e^{-r^2/2}dr$. Aquí $I_1=\sqrt{\pi/2}$, $I_2=1$, $I_{n+2}=nI_n$, luego $s_1=2$, $s_2=2\pi$, $s_3=4\pi$,
$s_4=2\pi^2$; en general $s_n=2\pi^{n/2}/\Gamma(n/2)$.

*Demostración.* La fórmula es el teorema A.18 con $\phi=1$, $\mu(S^{n-1})=V_n$. Para $h(r)=e^{-r^2/2}$ el lado izquierdo es
$\big(\int_{\mathbb R}e^{-t^2/2}dt\big)^n=(2\pi)^{n/2}$ (Fubini), lo que da $s_n$. $I_1=\frac12\sqrt{2\pi}$, $I_2=1$, y
la integración por partes da $I_{n+2}=\big[-r^ne^{-r^2/2}\big]_0^\infty+n\int_0^\infty r^{n-1}e^{-r^2/2}dr=nI_n$; el cambio
$r=\sqrt{2u}$ da $I_n=2^{n/2-1}\Gamma(n/2)$. $\square$

Con $n=2$ y $n=3$ estas son las coordenadas polares y esféricas de la proposición 5.13.

**Definición.** La **distribución uniforme** sobre $S^{n-1}$ es $\sigma_{n-1}(B)=\lambda_n(C(B))/V_n$. Sobre la esfera unidad de
$\mathbb C^N\cong\mathbb R^{2N}$ es $\sigma_{2N-1}$.

**Teorema A.20 (descomposición polar de las leyes radiales).** Sea $X$ un vector aleatorio en $\mathbb R^n$ con una densidad de la
forma $f(|x|)$. Entonces $X\neq0$ casi seguramente, $|X|$ tiene la densidad $s_nr^{n-1}f(r)$ sobre $(0,\infty)$, $X/|X|$ tiene la ley
$\sigma_{n-1}$, y $|X|$ y $X/|X|$ son independientes. Además, $\sigma_{n-1}(OB)=\sigma_{n-1}(B)$ para toda $O\in O(n)$.
En particular:

(a) para un vector gaussiano estándar $g$ de $\mathbb R^n$, $g/|g|\sim\sigma_{n-1}$, independientemente de $|g|$;

(b) para un vector gaussiano complejo estándar $z\in\mathbb C^N$, $\psi=z/|z|\sim\sigma_{2N-1}$, una ley invariante por $U(N)$,
independiente de $|z|$; y $|z|^2$ tiene la densidad $t^{N-1}e^{-t}/(N-1)!$ sobre $(0,\infty)$, luego $\mathbb E|z|^4=N(N+1)$.

*Demostración.* Por el teorema A.18 con $h=f\cdot1_A$, $\phi=1_B$:
$\mathbb P(|X|\in A,\ X/|X|\in B)=n\,\lambda_n(C(B))\int_Af(r)r^{n-1}dr$. Con $A=(0,\infty)$ y $B=S^{n-1}$ esto es
$1=nV_n\int_0^\infty f(r)r^{n-1}dr$ (y $\mathbb P(X=0)=0$). Tomando $A=(0,\infty)$ se obtiene $\mathbb P(X/|X|\in B)=\sigma_{n-1}(B)$,
tomando $B=S^{n-1}$ se obtiene $\mathbb P(|X|\in A)=s_n\int_Af(r)r^{n-1}dr$, y el caso general es el producto de ambos.
Para $O\in O(n)$, $OC(B)=C(OB)$ y $\lambda_n(OE)=\lambda_n(E)$ (lema A.3(b)). (a) $\varphi_n(x)=(2\pi)^{-n/2}e^{-|x|^2/2}$ es radial.
(b) $z=g/\sqrt2$ con $g$ estándar en $\mathbb R^{2N}$ (proposición A.16), luego $z/|z|=g/|g|$ y $|z|=|g|/\sqrt2$; las unitarias son
ortogonales sobre $\mathbb R^{2N}$. La densidad de $|g|$ es $s_{2N}(2\pi)^{-N}r^{2N-1}e^{-r^2/2}$ con $s_{2N}=2\pi^N/(N-1)!$, y el
cambio $t=r^2/2$ da la densidad de $|z|^2$; entonces $\mathbb E|z|^4=\int_0^\infty t^{N+1}e^{-t}dt/(N-1)!=(N+1)N$. $\square$

El lema 7.9 es la parte (b). Para $N=1$, $\sigma_1$ es la longitud de arco normalizada, es decir, la medida imagen de
$\frac1{2\pi}\lambda_1$ sobre $[0,2\pi)$ por $\theta\mapsto e^{i\theta}$: para un boreliano $B\subseteq S^1$ sea
$\Theta(B)=\{\theta\in(0,2\pi):e^{i\theta}\in B\}$. La aplicación $(r,\theta)\mapsto re^{i\theta}$ es una biyección $C^1$ del abierto
$(0,1)\times(0,2\pi)$ sobre el disco unidad abierto menos el segmento $[0,1)$, con jacobiano $r$, y $C(B)$ difiere de la
imagen de $(0,1)\times\Theta(B)$ en un subconjunto de ese segmento y de la circunferencia unidad, un conjunto nulo. Por la fórmula
del cambio de variables, $\lambda_2(C(B))=\int_{\Theta(B)}\int_0^1r\,dr\,d\theta=\frac12\lambda_1(\Theta(B))$ y, con $V_2=\pi$,
$\sigma_1(B)=\lambda_1(\Theta(B))/2\pi$ (el punto $\theta=0$ es un conjunto nulo).

**Corolario A.21 (momentos de un vector unitario uniforme).** Si $u\sim\sigma_{2N-1}$ sobre la esfera unidad de $\mathbb C^N$, entonces
$\mathbb E[u_a\bar u_b]=\delta_{ab}/N$, $\mathbb E|u_1|^4=\frac2{N(N+1)}$ y $\mathbb E|u_1|^2|u_2|^2=\frac1{N(N+1)}$ ($N\ge2$). Si
$u\sim\sigma_{n-1}$ sobre $S^{n-1}\subset\mathbb R^n$, entonces $\mathbb Eu_1^2=\frac1n$, $\mathbb Eu_1^4=\frac3{n(n+2)}$ y
$\mathbb Eu_1^2u_2^2=\frac1{n(n+2)}$ ($n\ge2$).

*Demostración.* Tomemos $u=z/|z|$ como en el teorema A.20(b). Como $z_a\bar z_b=|z|^2u_a\bar u_b$ y $|z|$, $u$ son independientes,
$\delta_{ab}=\mathbb E[z_a\bar z_b]=\mathbb E|z|^2\,\mathbb E[u_a\bar u_b]=N\,\mathbb E[u_a\bar u_b]$ ($\mathbb E|z|^2=N$ por la
proposición A.16(a)). Del mismo modo, $\mathbb E|z_1|^4=N(N+1)\mathbb E|u_1|^4$ y
$\mathbb E|z_1|^2|z_2|^2=N(N+1)\mathbb E|u_1|^2|u_2|^2$, donde $\mathbb E|z_1|^2|z_2|^2=1$ por independencia y
$\mathbb E|z_1|^4=\mathbb E(x^2+y^2)^2=\frac34+2\cdot\frac14+\frac34=2$ para $x,y\sim N(0,\frac12)$ independientes. (Para
$Y\sim N(0,s^2)$, la integración por partes da $\mathbb EY^4=3s^2\,\mathbb EY^2=3s^4$, ya que
$y^3\cdot ye^{-y^2/2s^2}=-s^2y^3\frac d{dy}e^{-y^2/2s^2}$.) En el caso real, $\mathbb E|g|^2=n$,
$\mathbb E|g|^4=s_n(2\pi)^{-n/2}I_{n+4}=I_{n+4}/I_n=n(n+2)$ (corolario A.19), $\mathbb Eg_1^4=3$, $\mathbb Eg_1^2g_2^2=1$, y
$|g|^2$, $|g|^4$ son independientes de $u=g/|g|$ (teorema A.20(a)). $\square$

Estos valores coinciden con el lema 6.14, que los deduce solo de la invariancia unitaria. Recíprocamente, la invariancia determina la ley:

**Proposición A.22 (unicidad de la ley invariante).** Sea $u$ un vector unitario aleatorio en $\mathbb C^N$ cuya ley es
invariante por todas las unitarias. Entonces $u\sim\sigma_{2N-1}$. Para $N=1$: toda medida de probabilidad boreliana sobre la
circunferencia invariante por rotaciones es la longitud de arco normalizada.

*Demostración.* Sea $Q$ una unitaria de Haar en $U(N)$ (teorema 5.3), independiente de $u$. Para un vector unitario fijo $v$ se elige
$W\in U(N)$ con $We_1=v$. Por el lema 5.1, $QW$ es de Haar, luego $Qv=(QW)e_1$ tiene la ley de $Qe_1$. En el
teorema 5.3, $Q$ es el factor de Gram–Schmidt de una matriz de Ginibre (lema 5.2), así que $Qe_1=z/|z|$ para su primera columna $z$, que tiene
la ley $\sigma_{2N-1}$ por el teorema A.20(b). Condicionando respecto de $u$ (Fubini, independencia), $\mathbb P(Qu\in B)=\sigma_{2N-1}(B)$
para todo boreliano $B$. Condicionando en cambio respecto de $Q$, y usando que $Uu$ tiene la ley de $u$ para cada unitaria fija $U$,
$\mathbb P(Qu\in B)=\mathbb P(u\in B)$. Una medida de probabilidad sobre la circunferencia es la ley de la aplicación identidad, y
$U(1)$ actúa por rotaciones. $\square$

```python
N, B = 3, 400_000                                   # Theorem A.20 (b), Corollary A.21, in C^3
z = (rng.standard_normal((B, N)) + 1j * rng.standard_normal((B, N))) / math.sqrt(2)
R2 = (np.abs(z) ** 2).sum(axis=1)
u = z / np.sqrt(R2)[:, None]
t = 2.0
p_gamma = 1 - math.exp(-t) * sum(t ** k / math.factorial(k) for k in range(N))   # P(|z|^2 <= t) for Gamma(N)
mc(R2 <= t, p_gamma, "P(|z|^2 <= 2)")
mc(R2 ** 2, N * (N + 1), "E|z|^4")
mc(np.abs(u[:, 0]) ** 4, 2 / (N * (N + 1)), "E|u_1|^4")
mc(np.abs(u[:, 0] * u[:, 1]) ** 2, 1 / (N * (N + 1)), "E|u_1|^2 |u_2|^2")
mc((R2 <= t) * np.abs(u[:, 0]) ** 4, p_gamma * 2 / (N * (N + 1)), "E[1{|z|^2 <= 2} |u_1|^4], independence")
mc((u[:, 0] * u[:, 1].conj()).real, 0.0, "E Re u_1 conj(u_2)")

x3 = rng.uniform(-1, 1, (B, 3))                     # Corollary A.19: V_3 = s_3 / 3 = 4 pi / 3
mc((x3 ** 2).sum(axis=1) <= 1, math.pi / 6, "P(uniform point of [-1,1]^3 in the unit ball)")
g3 = rng.standard_normal((B, 3))
u3 = g3 / np.linalg.norm(g3, axis=1, keepdims=True)
mc(u3[:, 0] ** 4, 3 / 15, "E u_1^4 on S^2")
phase = np.angle(z[:, 0])                           # sigma_1 is normalised arc length
mc((phase >= 0) & (phase < 1), 1 / (2 * math.pi), "P(arg z_1 in [0, 1))")
```

## A.6 Momentos gaussianos: la fórmula de Isserlis–Wick

**Lema A.23.** Sea $Y$ un vector gaussiano centrado en $\mathbb R^n$ con covarianza $C=(c_{ij})$. Entonces
$c_{ij}=\mathbb E[Y_iY_j]$, $\mathbb Ee^{\langle t,Y\rangle}=e^{\langle t,Ct\rangle/2}$ y $\mathbb Ee^{|\langle t,Y\rangle|}<\infty$ para todo
$t\in\mathbb R^n$, y todos los momentos $\mathbb E|Y_i|^k$ son finitos.

*Demostración.* Escribamos $Y=Ag$. Entonces $\mathbb E[Y_iY_j]=\sum_{k,l}A_{ik}A_{jl}\mathbb E[g_kg_l]=(AA^T)_{ij}$. Con $a=A^Tt$,
$\langle t,Y\rangle=\sum_ka_kg_k$ y, por independencia,
$\mathbb Ee^{\langle t,Y\rangle}=\prod_k\mathbb Ee^{a_kg_k}=\prod_ke^{a_k^2/2}=e^{|A^Tt|^2/2}=e^{\langle t,Ct\rangle/2}$. Finalmente,
$e^{|s|}\le e^s+e^{-s}$ y $|Y_i|^k\le k!\,e^{|Y_i|}$. $\square$

Un **emparejamiento** de $\{1,\dots,2p\}$ es una partición en $p$ bloques de dos elementos; hay $(2p-1)!!=(2p)!/(2^pp!)$ de ellos
(la pareja de $1$ puede elegirse de $2p-1$ maneras, y después se emparejan los $2p-2$ elementos restantes; inducción sobre $p$).

**Teorema A.24 (Isserlis–Wick).** Sean $Y_1,\dots,Y_m$ coordenadas de un vector gaussiano centrado (se admiten
repeticiones). Si $m$ es impar, $\mathbb E[Y_1\cdots Y_m]=0$; si $m=2p$,
$$
\mathbb E[Y_1\cdots Y_{2p}]=\sum_\pi\prod_{\{a,b\}\in\pi}\mathbb E[Y_aY_b],
$$
donde la suma recorre todos los emparejamientos $\pi$ de $\{1,\dots,2p\}$. La misma fórmula vale para variables aleatorias complejas
$Z_1,\dots,Z_m$ que son combinaciones lineales complejas de las coordenadas de un mismo vector gaussiano centrado (por ejemplo, las
entradas de una matriz de Ginibre o del GUE y sus conjugadas complejas), con $\mathbb E[Z_aZ_b]$ (sin conjugación) en los factores de cada pareja.

*Demostración.* El vector $(Y_1,\dots,Y_m)$ es la imagen del vector gaussiano dado por una aplicación lineal, luego es gaussiano
centrado (teorema A.15(a)); sea $C$ su covarianza. Fijemos $t\in\mathbb R^m$ y pongamos $S=\langle t,Y\rangle$. Para $s\in\mathbb R$, las
sumas parciales de $\sum_k(sS)^k/k!$ están acotadas por $e^{|sS|}$, que es integrable por el lema A.23, así que, por convergencia
dominada y el lema A.23,
$$
\sum_{k\ge0}\frac{s^k}{k!}\mathbb E[S^k]=\mathbb Ee^{sS}=e^{s^2q/2}=\sum_{p\ge0}\frac{s^{2p}q^p}{2^pp!},\qquad q=\langle t,Ct\rangle .
$$
Ambos lados son series de potencias sobre $\mathbb R$ en el sentido de la sección A.2: a la derecha esto es claro, y a la izquierda
$\sum_k\frac{\rho^k}{k!}|\mathbb E[S^k]|\le\mathbb E\sum_k\frac{\rho^k|S|^k}{k!}=\mathbb Ee^{\rho|S|}<\infty$ para todo $\rho>0$
(convergencia monótona y lema A.23). Así, sus coeficientes coinciden (lema A.5(c) aplicado a la diferencia). Por tanto, $\mathbb E\langle t,Y\rangle^m/m!$ vale $0$ para $m$ impar y $\langle t,Ct\rangle^p/(2^pp!)$ para $m=2p$, para todo
$t\in\mathbb R^m$. Ambos lados son polinomios en $t$, así que sus coeficientes coinciden (corolario A.7(a) aplicado a la diferencia).
En $\mathbb E\langle t,Y\rangle^m=\sum_{i_1,\dots,i_m}t_{i_1}\cdots t_{i_m}\mathbb E[Y_{i_1}\cdots Y_{i_m}]$, el monomio $t_1t_2\cdots t_m$
procede de las $m!$ sucesiones $(i_1,\dots,i_m)$ que son permutaciones de $(1,\dots,m)$, así que su coeficiente en
$\mathbb E\langle t,Y\rangle^m/m!$ es $\mathbb E[Y_1\cdots Y_m]$. Para $m$ impar es, por tanto, $0$. Para $m=2p$,
$q^p=\sum c_{i_1j_1}\cdots c_{i_pj_p}\,t_{i_1}t_{j_1}\cdots t_{i_p}t_{j_p}$ sobre todas las sucesiones de $p$ pares ordenados; el monomio
$t_1\cdots t_{2p}$ aparece exactamente cuando los pares no ordenados $\{i_k,j_k\}$ forman un emparejamiento $\pi$, y cada emparejamiento procede de
$p!\,2^p$ sucesiones (orden de los pares, orden dentro de cada par), cada una de las cuales contribuye con $\prod_{\{a,b\}\in\pi}c_{ab}$, pues $C$ es
simétrica. Dividiendo por $2^pp!$ se obtiene la fórmula, ya que $c_{ab}=\mathbb E[Y_aY_b]$.
Para el enunciado complejo, escribamos $Z_a=\sum_i\alpha_{ai}Y_i$. Ambos lados son multilineales en $(Z_1,\dots,Z_m)$ —el izquierdo por
la linealidad de $\mathbb E$, el derecho porque cada $Z_a$ aparece en exactamente un factor de cada término— y coinciden cuando cada
$Z_a$ es uno de los $Y_i$; basta desarrollar. $\square$

Por ejemplo, $\mathbb E|z|^4=\mathbb E[z\bar zz\bar z]=1\cdot1+0+1\cdot1=2$ para una gaussiana compleja estándar $z$. La proposición 1.10
y el lema 7.10 utilizan el caso $m=4$; el teorema A.24 lo demuestra sin derivar bajo la esperanza.

```python
def pairings(idx):
    """All pairings of the list idx, as lists of pairs."""
    if not idx:
        yield []
        return
    for i in range(1, len(idx)):
        for rest in pairings(idx[1:i] + idx[i + 1:]):
            yield [(idx[0], idx[i])] + rest


def wick(cov, idx):                                       # sum over pairings of the positions of idx
    return sum(math.prod(cov[idx[a], idx[b]] for a, b in p) for p in pairings(list(range(len(idx)))))


assert [len(list(pairings(list(range(2 * p))))) for p in range(1, 5)] == [1, 3, 15, 105]    # (2p-1)!!
Amat = rng.standard_normal((3, 3)) / math.sqrt(3)
cov = Amat @ Amat.T
Y = rng.standard_normal((1_000_000, 3)) @ Amat.T
for idx in ((0, 1, 2, 0), (0, 0, 0, 0), (0, 1, 1, 2, 2, 0)):
    mc(np.prod(Y[:, list(idx)], axis=1), wick(cov, idx), f"E Y{idx}")
mc(np.prod(Y[:, [0, 1, 2]], axis=1), 0.0, "E Y(0, 1, 2), odd")
M = np.array([[1, 1j], [1, -1j]])                  # (z, conj z) = M (x, y) with x, y independent N(0, 1/2)
cz = M @ M.T / 2                                     # E[Z_a Z_b] = [[0, 1], [1, 0]]
assert abs(wick(cz, (0, 1, 0, 1)) - 2) < 1e-12                           # E|z|^4 = 2
```

## A.7 Variables exponenciales y espaciamientos uniformes

Una variable aleatoria $e$ es **exponencial** ($e\sim\operatorname{Exp}(1)$) si tiene la densidad $e^{-y}$ sobre $(0,\infty)$;
$e/\beta$ es exponencial de tasa $\beta$. Sea $\Delta_n=\{x\in\mathbb R^n:0<x_1<\dots<x_n<1\}$.

**Teorema A.25.** Sean $e_0,\dots,e_n$ ($n\ge1$) variables $\operatorname{Exp}(1)$ independientes, $T=e_0+\dots+e_n$ y
$X_k=(e_0+\dots+e_{k-1})/T$ para $k=1,\dots,n$. Entonces $(X_1,\dots,X_n)$ tiene la densidad $n!$ sobre $\Delta_n$, $T$ tiene la densidad
$t^ne^{-t}/n!$ sobre $(0,\infty)$, y $(X_1,\dots,X_n)$ y $T$ son independientes. En particular,
$(e_0,\dots,e_n)/T$ es independiente de $T$.

*Demostración.* La aplicación $\Phi(x,t)=t\,(x_1,\ x_2-x_1,\ \dots,\ x_n-x_{n-1},\ 1-x_n)$ es una biyección del abierto
$\Delta_n\times(0,\infty)$ sobre $(0,\infty)^{n+1}$, con inversa $y\mapsto\big(((y_0+\dots+y_{k-1})/\sum_iy_i)_{k=1}^n,\ \sum_iy_i\big)$;
ambas son $C^1$. En la matriz jacobiana de $\Phi$, la columna de $\partial/\partial t$ tiene suma de entradas $1$ y cada columna
$\partial/\partial x_j$ (entradas $t$ y $-t$) tiene suma de entradas $0$. Sustituyendo la última fila por la suma de todas las filas (lo que
no cambia el determinante) se obtiene la fila $(0,\dots,0,1)$, y el bloque $n\times n$ restante es triangular inferior con
diagonal $t$; así, $|\det D\Phi|=t^n$. El vector $e=(e_0,\dots,e_n)$ tiene la densidad $e^{-\sum_iy_i}$, y
$\sum_i\Phi(x,t)_i=t$. Por la fórmula del cambio de variables, para todo boreliano $E\subseteq\Delta_n\times(0,\infty)$,
$$
\mathbb P\big((X,T)\in E\big)=\mathbb P\big(e\in\Phi(E)\big)=\int_{\Phi(E)}e^{-\sum_iy_i}dy=\int_Ee^{-t}t^n\,dx\,dt .
$$
La densidad $e^{-t}t^n=n!\cdot\frac{t^ne^{-t}}{n!}$ es un producto, y $\int_0^\infty t^ne^{-t}dt=n!$, luego
$\lambda_n(\Delta_n)=1/n!$; por Fubini, $\mathbb P(X\in A,T\in A')=\mathbb P(X\in A)\mathbb P(T\in A')$. Finalmente, $(e_0,\dots,e_n)/T$ es
$(X_1,X_2-X_1,\dots,1-X_n)$, una función de $X$. $\square$

**Corolario A.26 (espaciamientos uniformes).** Sean $U_1,\dots,U_n$ independientes y uniformes sobre $(0,1)$, con estadísticos de orden
$U_{(1)}<\dots<U_{(n)}$ y espaciamientos $s_0=U_{(1)}$, $s_k=U_{(k+1)}-U_{(k)}$ ($1\le k\le n-1$), $s_n=1-U_{(n)}$. Entonces:

(a) los empates tienen probabilidad $0$, y $(U_{(1)},\dots,U_{(n)})$ tiene la densidad $n!$ sobre $\Delta_n$;

(b) $(s_0,\dots,s_n)$ tiene la ley de $(e_0,\dots,e_n)/T$, con $e_i$, $T$ como en el teorema A.25;

(c) $\mathbb Es_i=\frac1{n+1}$, $\mathbb Es_i^2=\frac2{(n+1)(n+2)}$ y $\mathbb E[s_is_j]=\frac1{(n+1)(n+2)}$ para $i\neq j$;

(d) para $i\neq j$, $s_i/s_j$ tiene la ley de $e_i/e_j$, con $\mathbb P(e_i/e_j\le t)=\frac t{1+t}$ ($t>0$); así, para
$0\le i\le n-1$, el cociente de espaciamientos $\min(s_i,s_{i+1})/\max(s_i,s_{i+1})$ tiene la ley de $\min(e_i,e_{i+1})/\max(e_i,e_{i+1})$
obtenida en la proposición 5.14(a), $\mathbb P(r\le t)=\frac{2t}{1+t}$ sobre $[0,1]$, con media $2\ln2-1$.

Para $n$ puntos uniformes independientes en un intervalo $[a,b]$, la aplicación afín $u\mapsto a+(b-a)u$ multiplica todos los espaciamientos
por $b-a$: los espaciamientos divididos por $b-a$ cumplen (b) y (c), y los cocientes de (d) no cambian.

*Demostración.* (a) El conjunto donde dos coordenadas coinciden es una unión finita de hiperplanos, un conjunto nulo (corolario A.7(c)). Fuera
de él existe exactamente una permutación $\pi$ con $U_{\pi(1)}<\dots<U_{\pi(n)}$, así que, para todo boreliano $E\subseteq\Delta_n$, el suceso
$\{(U_{(1)},\dots,U_{(n)})\in E\}$ es, salvo un conjunto nulo, la unión disjunta sobre $\pi$ de los sucesos
$\{(U_{\pi(1)},\dots,U_{\pi(n)})\in E\}$. Cada uno tiene probabilidad $\lambda_n(E)$, porque $(U_1,\dots,U_n)$ tiene la densidad
$1_{(0,1)^n}$ y permutar las coordenadas conserva tanto $\lambda_n$ como $(0,1)^n$. Así, la probabilidad es $n!\,\lambda_n(E)$.
(b) Por (a) y el teorema A.25, los vectores $(U_{(k)})_k$ y $(X_k)_k$ tienen la misma ley, y los espaciamientos son la misma función
lineal de cualquiera de ellos. (c) Sea $w=(e_0,\dots,e_n)/T$, independiente de $T$ (teorema A.25), con $\mathbb ET=n+1$ y
$\mathbb ET^2=\int_0^\infty t^{n+2}e^{-t}dt/n!=(n+1)(n+2)$. De $e_i=Tw_i$: $1=\mathbb Ee_i=\mathbb ET\,\mathbb Ew_i$,
$2=\mathbb Ee_i^2=\mathbb ET^2\,\mathbb Ew_i^2$ y $1=\mathbb E[e_ie_j]=\mathbb ET^2\,\mathbb E[w_iw_j]$. (d) $s_i/s_j=w_i/w_j=e_i/e_j$ y, por
Tonelli, $\mathbb P(e_i\le te_j)=\int_0^\infty e^{-y}(1-e^{-ty})\,dy=1-\frac1{1+t}$. El cociente de espaciamientos es una función de
$s_i/s_{i+1}$ (a saber, $\min(x,1/x)$ en $x=s_i/s_{i+1}$, definida casi seguramente pues todos los $e_k>0$), así que tiene la ley de la
misma función de $e_i/e_{i+1}$, que es el cociente $r$ de la proposición 5.14(a) para las variables $\operatorname{Exp}(1)$ independientes
$e_i,e_{i+1}$. $\square$

En particular, consideremos los «niveles de Poisson» $\lambda_k=e_0+\dots+e_{k-1}$, $k=1,\dots,n+1$ (espaciamientos exponenciales
independientes $e_0,\dots,e_n$), y $n$ puntos uniformes ordenados con sus $n+1$ espaciamientos $s_0,\dots,s_n$. El vector de los
cocientes de espaciamientos consecutivos $r_i=\min(d_i,d_{i+1})/\max(d_i,d_{i+1})$, $0\le i<n$, tiene la misma ley conjunta para
$d=(e_0,\dots,e_n)$ y para $d=(s_0,\dots,s_n)$, porque los cocientes no cambian cuando todos los $d_i$ se dividen por $T$, y
$(e_0,\dots,e_n)/T$ tiene la ley de $(s_0,\dots,s_n)$ (parte (b)); y, por el teorema A.25, los $n$ primeros niveles de Poisson divididos por el $(n+1)$-ésimo tienen exactamente la ley de $n$ puntos
uniformes ordenados.

```python
n, B = 5, 200_000                                       # Corollary A.26
levels = np.sort(rng.uniform(0, 1, (B, n)), axis=1)
s = np.diff(np.concatenate([np.zeros((B, 1)), levels, np.ones((B, 1))], axis=1), axis=1)    # s_0, ..., s_n
e = rng.exponential(1.0, (B, n + 1))
w = e / e.sum(axis=1, keepdims=True)                    # normalised exponentials
for name, sp in (("uniform spacings", s), ("normalised exponentials", w)):
    mc(sp[:, 2], 1 / (n + 1), f"{name}: E s_2")
    mc(sp[:, 2] ** 2, 2 / ((n + 1) * (n + 2)), f"{name}: E s_2^2")
    mc(sp[:, 0] * sp[:, n], 1 / ((n + 1) * (n + 2)), f"{name}: E s_0 s_n")
    for tt in (0.5, 2.0):
        mc(sp[:, 1] / sp[:, 2] <= tt, tt / (1 + tt), f"{name}: P(s_1 / s_2 <= {tt})")
r = np.minimum(s[:, 1:-2], s[:, 2:-1]) / np.maximum(s[:, 1:-2], s[:, 2:-1])    # ratios of interior spacings
mc(r.mean(axis=1), 2 * math.log(2) - 1, "<r>, uniform levels")
T = e.sum(axis=1)
mc(w[:, 0] * T ** 2, n + 2, "E[w_0 T^2] = E w_0 E T^2 (independence)")
```

**Volver al** [índice](README.md).
