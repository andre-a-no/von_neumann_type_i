[English](../08_dynamics.md) · **Español** · [Français](../fr/08_dynamics.md) · [Deutsch](../de/08_dynamics.md) · [中文](../zh/08_dynamics.md) · [日本語](../ja/08_dynamics.md) · [Русский](../ru/08_dynamics.md)

# 8. Dinámica

En mecánica cuántica, un sistema cerrado se modeliza mediante un grupo uniparamétrico de unitarias; uno
abierto —un sistema acoplado a un entorno que no se observa—, mediante un semigrupo de canales. En un álgebra con
sectores ambos tienen un rasgo invisible en una sola álgebra de matrices: un hamiltoniano de
$`\mathcal M`$ nunca puede cambiar el sector, de modo que las cargas se conservan, mientras que un
sistema abierto puede saltar entre sectores, por ejemplo al perder una partícula. Este capítulo
demuestra los hechos básicos sobre ambos tipos de evolución, escribe la ecuación de Lindblad con
saltos entre sectores y muestra que "no se perdió ninguna partícula" es una medición de la carga
cuyo estado postseleccionado evoluciona según un hamiltoniano no hermítico.

**Aprenderá** por qué los pesos de los sectores son constantes del movimiento, por qué los
generadores en forma de Gorini–Kossakowski–Sudarshan–Lindblad (GKSL) producen canales, cómo entran
los saltos entre sectores en la ecuación de Lindblad y cuándo los pesos de los sectores obedecen una
ecuación de tasas clásica, qué significan la evolución sin saltos y la postselección, y cómo hallar
los estados estacionarios y las tasas de relajación. Las comprobaciones numéricas de este capítulo
están en el [Tutorial 8](../../tutorial/es/08_dynamics.md).

En todo el capítulo, $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$ actúa sobre
$`\mathcal H=\bigoplus_c\mathbb C^{k_c}`$, $`D=\sum_ck_c`$, $`P_c`$ son las proyecciones centrales
minimales, los estados se dan mediante sus densidades `blunt` $`\rho`$ con pesos de los sectores
$`p_c=\mathrm{Tr}\,\rho_c`$ (capítulo 3), $`[A,B]=AB-BA`$, $`\{A,B\}=AB+BA`$ y $`\hbar=1`$. Un
**superoperador** es una aplicación lineal sobre $`\mathcal M`$ o sobre $`M_D(\mathbb C)`$; todas las
normas de superoperadores son normas de operador y, por tanto, submultiplicativas. Para una
aplicación lineal $`B`$ sobre un espacio de dimensión finita, $`e^{tB}=\sum_{n\ge0}t^nB^n/n!`$.

<a id="s-8-1"></a>

## 8.1 Sistemas cerrados

Un **hamiltoniano** es un $`H\in\mathcal M`$ autoadjunto. Un vector evoluciona según la **ecuación de
Schrödinger** $`i\psi'(t)=H\psi(t)`$, y una densidad según la **ecuación de von Neumann**
$`\rho'(t)=-i[H,\rho(t)]`$. Ambas son ecuaciones lineales con coeficientes constantes sobre un espacio
de dimensión finita, y el lema siguiente es toda la teoría de ecuaciones diferenciales que
necesitamos.

<a id="r-8-1"></a>

**Lema 8.1.** Sea $`B`$ una aplicación lineal sobre un espacio normado de dimensión finita $`X`$.
Entonces $`\frac d{dt}e^{tB}=Be^{tB}=e^{tB}B`$, $`e^{(s+t)B}=e^{sB}e^{tB}`$,
$`\|e^{tB}\|\le e^{|t|\,\|B\|}`$ y, para todo $`x_0\in X`$, la aplicación $`x(t)=e^{tB}x_0`$ es la única
$`x`$ derivable en un intervalo $`I\ni0`$ con $`x'=Bx`$ y $`x(0)=x_0`$.

*Demostración.* La serie converge absolutamente, con la cota de la norma, y como serie de potencias
en $`t`$ puede derivarse término a término, lo que da la primera fórmula. Entonces
$`\frac d{dt}(e^{tB}e^{-tB})=Be^{tB}e^{-tB}-e^{tB}Be^{-tB}=0`$, luego $`e^{tB}e^{-tB}=1`$. Si $`x'=Bx`$,
entonces $`\frac d{dt}(e^{-tB}x(t))=e^{-tB}(x'-Bx)=0`$, luego $`e^{-tB}x(t)=x_0`$ y $`x(t)=e^{tB}x_0`$.
Por último, $`t\mapsto e^{(s+t)B}x_0`$ y $`t\mapsto e^{tB}e^{sB}x_0`$ resuelven $`y'=By`$ con el mismo
valor inicial. $`\square`$

<a id="r-8-2"></a>

**Proposición 8.2.** Sea $`H=H^*\in\mathcal M`$ y $`U(t)=e^{-itH}`$, definida mediante el cálculo
funcional de la sección 1.3 con $`f_t(\lambda)=e^{-it\lambda}`$.

1. $`U(t)=\bigoplus_ce^{-itH_c}`$ es una unitaria de $`\mathcal M`$, $`U(s+t)=U(s)U(t)`$ y
   $`U'(t)=-iHU(t)`$; $`U(t)`$ coincide con la serie exponencial de $`-itH`$.
2. $`\psi(t)=U(t)\psi_0`$ es la única solución de $`i\psi'=H\psi`$ con $`\psi(0)=\psi_0`$, y
   $`\rho(t)=U(t)\rho_0U(t)^*`$ es la única solución en $`\mathcal M`$ de $`\rho'=-i[H,\rho]`$ con
   $`\rho(0)=\rho_0`$. Si $`\rho_0`$ es un estado, también lo es cada $`\rho(t)`$.
3. (*Imagen de Heisenberg*) $`\mathrm{Tr}(\rho(t)A)=\mathrm{Tr}(\rho_0A(t))`$ con
   $`A(t)=U(t)^*AU(t)`$, y $`A'=i[H,A]`$.

*Demostración.* (1) Diagonalicemos cada bloque, $`H_c=V_c\Lambda_cV_c^*`$; entonces $`V=\bigoplus_cV_c`$
diagonaliza $`H`$, y $`U(t)=Ve^{-it\Lambda}V^*`$ es diagonal por bloques (por el teorema 1.5(a), el
resultado no depende de la unitaria que diagonaliza). Como $`|f_t|^2=1`$ y $`f_sf_t=f_{s+t}`$, el
teorema 1.5(c) da la unitariedad y la ley de grupo. Derivando $`e^{-it\Lambda}`$ entrada a entrada se
obtiene $`U'=V(-i\Lambda)e^{-it\Lambda}V^*=-iHU`$, luego $`U(t)`$ es la serie $`e^{-itH}`$ por el lema 8.1
($`X=M_D`$, $`B`$ la multiplicación por la izquierda por $`-iH`$). (2) $`(U\psi_0)'=-iHU\psi_0`$, y
$`(U\rho_0U^*)'=-iHU\rho_0U^*+U\rho_0U^*(iH)`$, porque $`(U')^*=iU^*H`$. La unicidad es el lema 8.1 con
$`B=-iH`$ sobre $`\mathcal H`$, o bien con $`B=-i[H,\cdot\,]`$ sobre $`\mathcal M`$. Positividad:
$`\langle x,U\rho_0U^*x\rangle=\langle U^*x,\rho_0U^*x\rangle\ge0`$; la traza, por ciclicidad. (3)
Ciclicidad de la traza, y $`(U^*AU)'=iU^*HAU-iU^*AHU=i[H,A(t)]`$. $`\square`$

*En la biblioteca:* $`U(t)=e^{-itH}`$, bloque a bloque, es `dynamics.propagator`; las soluciones de la proposición 8.2(2) son `dynamics.schrodinger` y `dynamics.von_neumann` — [tutorial, sección 8.1](../../tutorial/es/08_dynamics.md#s-8-1).

<a id="r-8-3"></a>

**Teorema 8.3 (conservación de la carga).** Sea $`t\mapsto H(t)`$, definida en un intervalo $`I`$, con
valores en los elementos autoadjuntos de $`\mathcal M`$.

1. Si $`\psi:I\to\mathcal H`$ es derivable con $`i\psi'=H(t)\psi`$, entonces cada $`\|P_c\psi(t)\|`$ es
   constante.
2. Si $`\rho:I\to\mathcal M`$ es derivable con $`\rho'=-i[H(t),\rho]`$, entonces
   $`\mathrm{Tr}(Z\rho(t))`$ es constante para todo $`Z`$ central; en particular, lo es cada peso
   de sector $`p_c(t)=\mathrm{Tr}(P_c\rho(t))`$.
3. Si $`H`$ es constante y $`A\in\mathcal M`$ conmuta con $`H`$, entonces $`\mathrm{Tr}(A\rho(t))`$ es
   constante.

*Demostración.* (1) El producto escalar es antilineal en el primer argumento, así que
$`\frac d{dt}\langle\psi,P_c\psi\rangle=\langle-iH\psi,P_c\psi\rangle+\langle\psi,-iP_cH\psi\rangle =i\langle\psi,HP_c\psi\rangle-i\langle\psi,P_cH\psi\rangle=0`$,
pues $`P_c`$ es central. (2)
$`\frac d{dt}\mathrm{Tr}(Z\rho)=-i\mathrm{Tr}(Z[H,\rho])=-i\mathrm{Tr}([Z,H]\rho)=0`$
por ciclicidad. (3) El mismo cálculo con $`A`$. $`\square`$

> **Significado físico** (interpretación física, no un teorema de este curso). Un hamiltoniano de
> $`\mathcal M`$ conmuta con las cargas (número de partículas, magnetización, un número cuántico de
> simetría), de modo que la probabilidad de cada carga es una constante del movimiento, incluso para
> hamiltonianos forzados y dependientes del tiempo.

Sobre $`\mathcal M`$, un vector unitario $`\psi=\bigoplus_c\psi_c`$ define el estado de la mezcla
$`\bigoplus_c|\psi_c\rangle\langle\psi_c|`$: todo $`A\in\mathcal M`$ es diagonal por bloques, así que
$`\langle\psi,A\psi\rangle=\sum_c\langle\psi_c,A_c\psi_c\rangle=\mathrm{Tr}\big(A\bigoplus_c|\psi_c\rangle\langle\psi_c|\big)`$
(como en el capítulo 3, ejercicio 3). Si $`\psi(t)`$ resuelve $`i\psi'=H\psi`$, entonces
$`\bigoplus_c|\psi_c(t)\rangle\langle\psi_c(t)|=U(t)\big(\bigoplus_c|\psi_c(0)\rangle\langle\psi_c(0)|\big)U(t)^*`$, ya que $`U(t)`$ es
diagonal por bloques.

*Cálculo:* el teorema 8.3 y la coincidencia de las evoluciones de $`\psi`$ y de $`\bigoplus_c|\psi_c\rangle\langle\psi_c|`$, con los solucionadores exactos y con el solucionador de Runge–Kutta `dynamics.schrodinger_rk4` para un hamiltoniano forzado — [tutorial, sección 8.1](../../tutorial/es/08_dynamics.md#s-8-1).

<a id="s-8-2"></a>

## 8.2 Semigrupos dinámicos cuánticos

Las aplicaciones CP y TP y los canales se definen en la sección 6.1 ([capítulo
6](06_quantum_channels.md)). Además del teorema 6.2 y de la proposición 6.4 (una aplicación
$`A\mapsto\sum_lK_lAK_l^*`$, también entre sectores, es CP), necesitamos tres propiedades de clausura.

<a id="r-8-4"></a>

**Lema 8.4.** (a) Las sumas, los múltiplos no negativos y las composiciones de aplicaciones CP son
CP. Sea $`\Phi`$ CP sobre $`\mathcal A=M_D`$ o sobre $`\mathcal A=\mathcal M`$, y sea
$`\mathcal N\subseteq\mathcal A`$ igual a $`\mathcal M`$ (si $`\mathcal A=M_D`$) o a una de las esquinas
$`P_c\mathcal MP_c\cong M_{k_c}`$ (con unidad $`P_c`$). Si $`\Phi`$ lleva $`\mathcal N`$ en sí misma, entonces
$`\Phi|_{\mathcal N}`$ es CP. (b) Si unas aplicaciones CP $`\Phi_m`$ convergen a $`\Phi`$ (es decir,
$`\Phi_m(A)\to\Phi(A)`$ para todo $`A`$), entonces $`\Phi`$ es CP. (c) Si $`s\mapsto\Psi(s)`$ es una
aplicación continua de $`[a,b]`$ en las aplicaciones CP, entonces $`\int_a^b\Psi(s)\,ds`$ es CP.

*Demostración.* (a) $`\Phi\mapsto\Phi_n`$ es lineal y $`(\Phi\Psi)_n=\Phi_n\Psi_n`$; las sumas de
matrices positivas son positivas. Para la restricción: $`M_n(\mathcal N)\subseteq M_n(\mathcal A)`$, y
la positividad es la misma condición en ambos. En efecto, por la sección 6.1, un elemento $`[A_{ij}]`$ de
$`M_n(\mathcal M)`$ es positivo si y solo si lo es cada matriz por bloques $`[(A_{ij})_c]`$; para
$`\mathcal A=M_D`$ esto es la positividad de $`[A_{ij}]\in M_{nD}`$, porque tras una permutación de la base
$`[A_{ij}]`$ es la suma directa de estos bloques (todos los $`A_{ij}`$ son diagonales por bloques); y un
elemento de $`M_n(P_c\mathcal MP_c)`$ tiene el único bloque no nulo $`[(A_{ij})_c]`$, cuya positividad es
la positividad en $`M_n(M_{k_c})`$. (b) Las matrices positivas forman un conjunto cerrado
($`\langle x,Xx\rangle\ge0`$ para cada $`x`$ es una condición cerrada), y $`(\Phi_m)_n(X)\to\Phi_n(X)`$. (c)
Las sumas de Riemann son CP por (a) y convergen a la integral (lema A.4(a)); se aplica (b). $`\square`$

**Definición.** Un **semigrupo dinámico cuántico** sobre $`\mathcal M`$ es una familia $`(T_t)_{t\ge0}`$
de canales $`\mathcal M\to\mathcal M`$ con $`T_0=\mathrm{id}`$, $`T_{s+t}=T_sT_t`$ y $`t\mapsto T_t`$
continua.

<a id="r-8-5"></a>

**Lema 8.5.** Todo semigrupo continuo $`(T_t)_{t\ge0}`$ de superoperadores ($`T_0=\mathrm{id}`$,
$`T_{s+t}=T_sT_t`$) es de la forma $`T_t=e^{t\mathcal L}`$ para exactamente un superoperador $`\mathcal L`$, su
**generador**, y $`\mathcal L=\lim_{t\to0}(T_t-\mathrm{id})/t`$. Lo mismo vale, con la misma
demostración, para los semigrupos continuos de aplicaciones lineales sobre cualquier espacio de
dimensión finita (por ejemplo, de matrices).

*Demostración.* (Las integrales de funciones continuas con valores superoperadores son las del lema
A.4.) $`\|V_h-\mathrm{id}\|\le\max_{0\le s\le h}\|T_s-\mathrm{id}\|\to0`$ para
$`V_h=\frac1h\int_0^hT_s\,ds`$ (lema A.4(c)), luego $`\det V_h\neq0`$ para algún $`h`$ pequeño fijo, por la
continuidad del determinante. Ahora bien,
$`T_tV_h=\frac1h\int_0^hT_{t+s}\,ds=\frac1h\int_t^{t+h}T_s\,ds`$ (lema A.4(b),(e)) es derivable en $`t`$ con
derivada $`\frac1h(T_{t+h}-T_t)=T_t\frac{T_h-\mathrm{id}}h`$. Por tanto $`T_t`$ es derivable y
$`T_t'=T_t\mathcal L`$ con $`\mathcal L=\frac{T_h-\mathrm{id}}hV_h^{-1}`$, y $`T_t=e^{t\mathcal L}`$ por
el lema 8.1 (con $`B`$ la composición por la derecha con $`\mathcal L`$). Unicidad: $`\mathcal L=T'(0)`$.
$`\square`$

<a id="r-8-6"></a>

**Teorema 8.6 (generadores GKSL).** Sean $`H=H^*`$ y $`L_1,\dots,L_J`$ en $`M_D`$, y
```math
\mathcal L(\rho)=-i[H,\rho]+\sum_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
```
Entonces $`(e^{t\mathcal L})_{t\ge0}`$ es un semigrupo dinámico cuántico sobre $`M_D`$. Más
precisamente, pongamos $`H_{\mathrm{eff}}=H-\frac i2\sum_jL_j^*L_j`$,
$`\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)`$ y
$`\Phi(\rho)=\sum_jL_j\rho L_j^*`$, de modo que $`\mathcal L=\mathcal L_0+\Phi`$, y definamos
```math
S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},\qquad
S_n(t)=\int_0^tS_0(t-s)\,\Phi\,S_{n-1}(s)\,ds\quad(n\ge1).
```
Entonces cada $`S_n(t)`$ es CP, y $`e^{t\mathcal L}=\sum_{n\ge0}S_n(t)`$ (el **desarrollo de Dyson**),
con $`\|S_n(t)\|\le e^{t\|\mathcal L_0\|}(t\|\Phi\|)^n/n!`$.

*Demostración.* *Paso 1.* $`X(t)=e^{-itH_{\mathrm{eff}}}\rho e^{itH_{\mathrm{eff}}^*}`$ cumple
$`X'=\mathcal L_0(X)`$ y $`X(0)=\rho`$, luego $`S_0(t)=e^{t\mathcal L_0}`$ por el lema 8.1; es CP con el
único operador de Kraus $`e^{-itH_{\mathrm{eff}}}`$. *Paso 2.* Escribamos $`a=\|\mathcal L_0\|`$,
$`b=\|\Phi\|`$. Por inducción, las $`S_n`$ son continuas y
$`\|S_n(t)\|\le\int_0^te^{(t-s)a}\,b\,e^{sa}\frac{(bs)^{n-1}}{(n-1)!}ds=e^{ta}\frac{(bt)^n}{n!}`$, de
modo que la serie converge uniformemente en todo $`[0,T]`$. *Paso 3.* Como
$`S_n(t)=S_0(t)\int_0^te^{-s\mathcal L_0}\Phi S_{n-1}(s)\,ds`$, la regla del producto y el teorema
fundamental del cálculo (lema A.4(e)) dan $`S_n'=\mathcal L_0S_n+\Phi S_{n-1}`$ ($`n\ge1`$) y $`S_0'=\mathcal L_0S_0`$.
Por las cotas, la serie de las derivadas converge uniformemente en $`[0,T]`$, luego
$`Y(t)=\sum_nS_n(t)`$ es derivable con $`Y'=\mathcal L_0Y+\Phi Y=\mathcal LY`$, $`Y(0)=\mathrm{id}`$, e
$`Y(t)=e^{t\mathcal L}`$ por el lema 8.1. *Paso 4.* $`\Phi`$ es CP (teorema 6.2), así que, por
inducción, el integrando $`s\mapsto S_0(t-s)\Phi S_{n-1}(s)`$ es una familia continua de aplicaciones
CP y $`S_n(t)`$ es CP (lema 8.4); también lo son las sumas parciales y su límite $`e^{t\mathcal L}`$.
*Paso 5.* Por ciclicidad,
$`\mathrm{Tr}\,\mathcal L(\rho)=\sum_j\big(\mathrm{Tr}(L_j\rho L_j^*)-\mathrm{Tr}(L_j^*L_j\rho)\big)=0`$,
luego
$`\frac d{dt}\mathrm{Tr}\,e^{t\mathcal L}\rho=\mathrm{Tr}\,\mathcal L(e^{t\mathcal L}\rho)=0`$:
cada $`e^{t\mathcal L}`$ es TP. La ley de semigrupo y la continuidad son el lema 8.1. $`\square`$

*Cálculo:* el teorema 8.6 para saltos arbitrarios sobre $`M_3`$, incluida la positividad completa de la parte de saltos $`\sum_{n\ge1}S_n(t)`$ — [tutorial, sección 8.2](../../tutorial/es/08_dynamics.md#s-8-2). *En la biblioteca:* $`e^{t\mathcal L}`$ es `dynamics.lindblad_channel`, y su matriz $`\mathcal L`$ es `dynamics.lindblad_superoperator` — [tutorial, sección 8.2](../../tutorial/es/08_dynamics.md#s-8-2).

<a id="r-8-7"></a>

**Teorema 8.7 (Gorini–Kossakowski–Sudarshan; Lindblad).** (1) Recíprocamente, el generador de todo
semigrupo dinámico cuántico sobre $`M_D`$ tiene la forma del teorema 8.6, con $`J\le D^2-1`$ saltos. Se
puede elegir $`\mathrm{Tr}\,H=0`$, $`\mathrm{Tr}\,L_j=0`$ y $`\mathrm{Tr}(L_j^*L_k)=0`$ para
$`j\neq k`$; con estas normalizaciones, $`H`$ y $`\Phi=\sum_jL_j\cdot L_j^*`$ están determinados por
$`\mathcal L`$, y los $`L_j`$ lo están salvo $`L'_j=\sum_ku_{jk}L_k`$ con una unitaria $`u`$ (tras completar con
ceros la familia más corta).

(2) El generador de todo semigrupo dinámico cuántico sobre $`\mathcal M=\bigoplus_cM_{k_c}`$ es
```math
\mathcal L(\rho)=\Phi(\rho)-\tfrac12\{\Phi^*(1),\rho\}-i[H,\rho]
```
con una aplicación CP $`\Phi:\mathcal M\to\mathcal M`$ (operadores de Kraus entre los sectores,
proposición 6.4) y $`H=H^*\in\mathcal M`$: la forma de la sección 8.3, con un único salto $`\Phi_1=\Phi`$
entre sectores y tasa $`1`$.

*Demostración.* (1) es el teorema D.4 junto con las afirmaciones de unicidad del teorema D.5(b),(c), y
(2) es el teorema D.7, ambos demostrados en el [apéndice D](D_gksl.md). $`\square`$

El paso clave (teorema D.3) es que la matriz de Choi del generador de un semigrupo de aplicaciones CP
es positiva sobre el complemento ortogonal del vector máximamente entrelazado; la propiedad de
semigrupo por sí sola, con aplicaciones $`T_t`$ positivas, no bastaría (el ejemplo que sigue al teorema
D.3). Para (2), las aplicaciones $`T_t\circ E`$ con la esperanza condicional $`E(X)=\sum_cP_cXP_c`$ son CP
sobre $`M_D`$; la matriz de Choi de la derivada $`\mathcal L\circ E`$ en $`t=0`$ es positiva sobre el
complemento ortogonal de los vectores $`(1\otimes P_c)\Omega`$, lo que da una descomposición de
$`\mathcal L\circ E`$ (lema D.2), que $`E`$ comprime de vuelta a $`\mathcal M`$. En la imagen de Heisenberg
(corolario D.8), el generador de un semigrupo continuo de aplicaciones CP unitales sobre $`\mathcal M`$ es
$`A\mapsto\Psi(A)+K^*A+AK`$ con $`\Psi`$ CP y $`K\in\mathcal M`$. Históricamente, (1) se debe a V. Gorini,
A. Kossakowski, E. C. G. Sudarshan [*J. Math. Phys.* 17, 821 (1976)] y G. Lindblad [*Commun. Math.
Phys.* 48, 119 (1976)], y la forma en la imagen de Heisenberg, para semigrupos continuos en norma sobre
álgebras de von Neumann generales, a E. Christensen y D. E. Evans [*J. London Math. Soc.* 20, 358
(1979)].

> **Significado físico** (interpretación física, no un teorema de este curso). Se espera que un
> sistema que interactúa débilmente con un entorno grande cuyas correlaciones decaen rápidamente (la
> aproximación de Born–Markov) pierda la memoria de su pasado: su estado en el instante $`s+t`$ se
> obtiene entonces del estado en el instante $`s`$ mediante un canal que solo depende de $`t`$. *Si* es
> así, el teorema 8.7 dice que la dinámica es una parte hamiltoniana más unos "saltos" $`L_j`$ (emisión
> de fotones, pérdida de partículas, desfase) con tasas fijadas por la normalización de los $`L_j`$.

<a id="s-8-3"></a>

## 8.3 Ecuaciones de Lindblad en un álgebra con sectores

Sobre $`\mathcal M`$ hay dos tipos de saltos. Un salto **dentro de los sectores** es un operador
$`L_j\in\mathcal M`$. Un salto **entre sectores** es una aplicación CP $`\Phi:\mathcal M\to\mathcal M`$
con operadores de Kraus $`K^{dc}_l:\mathbb C^{k_c}\to\mathbb C^{k_d}`$ (proposición 6.4, un
`InterSectorChannel`): $`\Phi(\rho)_d=\sum_{c,l}K^{dc}_l\rho_cK^{dc*}_l`$. Por la proposición 6.8, su
dual `blunt` cumple $`\Phi^*(1)_c=\sum_{d,l}K^{dc*}_lK^{dc}_l`$, un elemento diagonal por bloques de
$`\mathcal M`$. Con tasas $`g_j,g_m>0`$, el generador de `lindblad_rhs` es
```math
\mathcal L(\rho)=-i[H,\rho]+\sum_jg_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big)
+\sum_mg_m\Big(\Phi_m(\rho)-\tfrac12\{\Phi_m^*(1),\rho\}\Big).
```
Por la proposición 8.8(1) siguiente, genera un semigrupo dinámico cuántico sobre $`\mathcal M`$;
recíprocamente, por el teorema 8.7(2) (teorema D.7), todo semigrupo dinámico cuántico sobre $`\mathcal M`$
tiene un generador de esta forma.

*En la biblioteca:* un salto entre sectores es un `InterSectorChannel`, y $`\mathcal L(\rho)`$ es `dynamics.lindblad_rhs` — [tutorial, sección 8.3](../../tutorial/es/08_dynamics.md#s-8-3).

<a id="r-8-8"></a>

**Proposición 8.8.** Sea $`J^{dc}_l\in M_D`$ la matriz con el bloque $`K^{dc}_l`$ en la posición $`(d,c)`$
y ceros en el resto.

1. $`\mathcal L`$ lleva $`\mathcal M`$ en $`\mathcal M`$ y es la restricción a $`\mathcal M`$ del generador
   GKSL $`\hat{\mathcal L}`$ sobre $`M_D`$ con saltos $`\sqrt{g_j}L_j`$ y $`\sqrt{g_m}J^{dc}_{m,l}`$ (todos
   los $`m,l,d,c`$). Por tanto, $`(e^{t\mathcal L})_{t\ge0}`$ es un semigrupo dinámico cuántico sobre
   $`\mathcal M`$.
2. Sea en cambio $`\tilde{\mathcal L}`$ el generador GKSL sobre $`M_D`$ en el que los bloques de un
   mismo índice de Kraus se suman en un único salto, $`L_{m,l}=\sum_{d,c}J^{dc}_{m,l}`$ (con el mismo
   hamiltoniano, los mismos saltos $`\sqrt{g_j}L_j`$ y los saltos $`\sqrt{g_m}L_{m,l}`$). Entonces
   $`\mathcal L(\rho)=E\big(\tilde{\mathcal L}(\rho)\big)`$ para $`\rho\in\mathcal M`$, donde
   $`E(X)=\sum_dP_dXP_d`$ es la esperanza condicional de $`M_D`$ sobre $`\mathcal M`$; en general,
   $`\tilde{\mathcal L}(\rho)\notin\mathcal M`$.
3. Si cada $`L_{m,l}`$ tiene a lo sumo un bloque no nulo en cada fila de bloques y en cada columna de
   bloques (por ejemplo, si solo lleva el sector $`c`$ al sector $`c-1`$, como los operadores de
   aniquilación del ejemplo siguiente), entonces $`\tilde{\mathcal L}(\rho)=\mathcal L(\rho)`$ para
   $`\rho\in\mathcal M`$.

*Demostración.* Como $`\rho=\sum_cP_c\rho P_c`$ y $`J^{dc}P_{c'}=\delta_{cc'}J^{dc}`$,
$`J^{dc}\rho J^{dc*}`$ es el bloque $`K^{dc}\rho_cK^{dc*}`$ en la posición $`(d,d)`$, y $`J^{dc*}J^{dc}`$ es
$`K^{dc*}K^{dc}`$ en la posición $`(c,c)`$. Sumando sobre $`d,l`$ se obtienen $`\Phi(\rho)`$ y $`\Phi^*(1)`$,
lo que prueba la identidad de (1), y todos los términos están en $`\mathcal M`$. Así,
$`\hat{\mathcal L}`$ y todas sus potencias llevan $`\mathcal M`$ en $`\mathcal M`$,
$`e^{t\hat{\mathcal L}}|_{\mathcal M}=e^{t\mathcal L}`$, y esto es un canal de $`\mathcal M`$ por el
teorema 8.6 y el lema 8.4(a). (2) $`L\rho L^*=\sum_{d,d',c}J^{dc}\rho J^{d'c*}`$ para
$`L=\sum_{d,c}J^{dc}`$, y $`E`$ conserva los términos con $`d=d'`$. Para $`B\in\mathcal M`$,
$`P_dXBP_d=P_dXP_dB`$, luego $`E(XB)=E(X)B`$ y $`E(BX)=BE(X)`$; así $`E\{L^*L,\rho\}=\{E(L^*L),\rho\}`$, y
$`E(L^*L)=\sum_cP_cL^*LP_c=\sum_{c,d}J^{dc*}J^{dc}`$ porque $`J^{dc*}J^{d'c}=0`$ para $`d\ne d'`$. Por
último, $`E(-i[H,\rho])=-i[H,\rho]`$, y los términos de los saltos $`\sqrt{g_j}L_j\in\mathcal M`$ están
en $`\mathcal M`$. (3) Sea $`L=\sum_{(d,c)\in S}J^{dc}`$, donde el conjunto $`S`$ de pares contiene a lo sumo
un par con un $`d`$ dado y a lo sumo uno con un $`c`$ dado. Para $`\rho\in\mathcal M`$,
$`J^{dc}\rho J^{d'c'*}=J^{dc}P_c\rho P_{c'}J^{d'c'*}`$ se anula salvo si $`c=c'`$, y entonces $`d=d'`$; y
$`J^{dc*}J^{d'c'}`$ se anula salvo si $`d=d'`$, y entonces $`c=c'`$. Así, $`L\rho L^*`$ y $`L^*L`$ son diagonales
por bloques, todos los términos de $`\tilde{\mathcal L}(\rho)`$ están en $`\mathcal M`$, y
$`\tilde{\mathcal L}(\rho)=E(\tilde{\mathcal L}(\rho))=\mathcal L(\rho)`$ por (2). $`\square`$

> **Significado físico** (interpretación física, no un teorema de este curso). La parte (2) dice que
> la ecuación de Lindblad sobre $`\mathcal M`$ es lo que queda de una ecuación sobre $`M_D`$ cuando las
> coherencias entre sectores que crean los saltos se borran de inmediato: el contenido de una regla de
> superselección. En la pérdida de partículas no se borra nada (parte (3)): la ecuación de Lindblad
> usual con los operadores de aniquilación como saltos ya es una ecuación sobre $`\mathcal M`$.

*Cálculo:* la proposición 8.8(2) para una aplicación aleatoria entre los sectores y la proposición 8.8(3) para la pérdida de partículas, frente a generadores GKSL sobre $`M_6`$ — [tutorial, sección 8.3](../../tutorial/es/08_dynamics.md#s-8-3).

<a id="r-8-9"></a>

**Proposición 8.9 (pesos de los sectores).** Sea $`\rho(t)=e^{t\mathcal L}\rho_0`$ y $`p(t)`$ sus pesos
de los sectores.

1. Si no hay saltos entre sectores ($`K^{dc}_{m,l}=0`$ para todos $`m,l`$ y $`d\ne c`$), cada $`p_c(t)`$ es
   constante.
2. Si para todo $`m`$ y todo par $`d\ne c`$ se tiene $`\sum_lK^{dc*}_{m,l}K^{dc}_{m,l}=r_{m,dc}1_{k_c}`$
   con números $`r_{m,dc}`$, entonces $`p'=Qp`$ con $`Q_{dc}=\sum_mg_mr_{m,dc}`$ para $`d\ne c`$ y
   $`Q_{cc}=-\sum_{d\ne c}Q_{dc}`$. Por tanto, $`p(t)=e^{tQ}p(0)`$ no depende de $`H`$, ni de los $`L_j`$,
   ni de los estados dentro de los sectores.

*Demostración.* $`p_d'=\mathrm{Tr}(P_d\mathcal L(\rho))`$ por el lema 8.1. El término hamiltoniano se anula como
en el teorema 8.3. Para $`L_j\in\mathcal M`$,
$`\mathrm{Tr}(P_dL_j\rho L_j^*)-\mathrm{Tr}(P_dL_j^*L_j\rho) =\mathrm{Tr}(L_j^*P_dL_j\rho)-\mathrm{Tr}(L_j^*L_jP_d\rho)=0`$,
pues $`P_d`$ es central. Para un salto $`\Phi`$ entre sectores,
$`\mathrm{Tr}(P_d\Phi(\rho))=\sum_c\mathrm{Tr}\big(\sum_lK^{dc*}_lK^{dc}_l\rho_c\big)`$ y
$`\mathrm{Tr}(P_d\tfrac12\{\Phi^*(1),\rho\})=\sum_{d'}\mathrm{Tr}\big(\sum_lK^{d'd*}_lK^{d'd}_l\rho_d\big)`$;
los términos $`c=d`$ y $`d'=d`$ se cancelan. Esto prueba (1), y bajo la hipótesis de (2) los términos
restantes son $`\sum_{c\ne d}r_{dc}p_c-\sum_{d'\ne d}r_{d'd}p_d`$. El lema 8.1 da $`p(t)=e^{tQ}p(0)`$.
$`\square`$

Bajo la hipótesis de (2), para $`d\ne c`$, $`r_{m,dc}`$ es la entrada $`T_{dc}(\rho)`$ de la matriz de
transición de $`\Phi_m`$ (definida antes de la proposición 6.12) para todo estado $`\rho`$ con $`p_c>0`$, ya
que
$`\mathrm{Tr}\,\Phi_{m,dc}(\rho_c)=\mathrm{Tr}\big(\sum_lK^{dc*}_{m,l}K^{dc}_{m,l}\rho_c\big)=r_{m,dc}p_c`$.
Así, para un único salto con tasa $`1`$, $`Q=R-\mathrm{diag}(\sum_dR_{dc})_c`$ con $`R`$ = `transition_matrix()`
(la diagonal de $`R`$ se cancela).

*En la biblioteca:* $`R`$ es `InterSectorChannel.transition_matrix` — [tutorial, sección 8.3](../../tutorial/es/08_dynamics.md#s-8-3).

**Ejemplo: pérdida de partículas.** Dos modos bosónicos con a lo sumo dos partículas; los sectores
son $`N=0,1,2`$ con bases $`|00\rangle`$; $`|10\rangle,|01\rangle`$; $`|20\rangle,|11\rangle,|02\rangle`$,
de modo que $`\mathcal M=M_1\oplus M_2\oplus M_3`$. El dímero de Bose–Hubbard
$`H=J(a_1^*a_2+a_2^*a_1)+\frac U2\sum_jn_j(n_j-1)`$ conmuta con $`N`$, y cada modo pierde partículas con
tasa $`\gamma`$: un salto entre sectores con operadores de Kraus $`a_1,a_2`$, restringidos a $`N\to N-1`$.
Como $`a_1^*a_1+a_2^*a_2=N`$, se aplica la proposición 8.9(2) con $`Q_{N-1,N}=\gamma N`$ y, partiendo de
$`N=2`$,
```math
p_2(t)=q^2,\qquad p_1(t)=2q(1-q),\qquad p_0(t)=(1-q)^2,\qquad q=e^{-\gamma t},
```
la ley binomial de parámetros $`2`$ y $`q`$ (la ley del número de supervivientes si cada partícula
sobrevive independientemente con probabilidad $`q`$), hagan lo que hagan el efecto túnel entre modos
(hopping) y la interacción: ambos miembros resuelven $`p'=Qp`$ (derívese) con $`p(0)=(0,0,1)`$, así que
coinciden por el lema 8.1.

*Cálculo:* el ejemplo, con $`\Phi^*(1)=0\oplus1\oplus2`$ y $`p(t)`$ comprobado frente a ambas fórmulas — [tutorial, sección 8.3](../../tutorial/es/08_dynamics.md#s-8-3).

<a id="s-8-4"></a>

## 8.4 Postselección y evolución sin saltos

Reunamos los saltos del generador de la sección 8.3 en
$`\Phi_{\mathrm{tot}}(\rho)=\sum_jg_jL_j\rho L_j^*+\sum_mg_m\Phi_m(\rho)`$, con
$`A=\Phi_{\mathrm{tot}}^*(1)=\sum_jg_jL_j^*L_j+\sum_mg_m\Phi_m^*(1)\ge0`$ en $`\mathcal M`$, y
$`H_{\mathrm{eff}}=H-\frac i2A`$. Apliquemos el teorema 8.6 al generador $`\hat{\mathcal L}`$ sobre $`M_D`$
de la proposición 8.8(1): por la demostración de la proposición 8.8, su parte de saltos lleva
$`\mathcal M`$ en $`\mathcal M`$, donde coincide con $`\Phi_{\mathrm{tot}}`$, y su suma $`\sum\hat L^*\hat L`$
sobre los saltos es $`A`$, así que su $`H_{\mathrm{eff}}`$ es $`H-\frac i2A\in\mathcal M`$. Por tanto, sus
$`S_n(t)`$ llevan $`\mathcal M`$ en $`\mathcal M`$ (por inducción sobre $`n`$), y para $`\rho\in\mathcal M`$
```math
\rho(t)=\sum_{n\ge0}S_n(t)\rho,\qquad S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},
\qquad S_n(t)=\int_0^tS_0(t-s)\,\Phi_{\mathrm{tot}}\,S_{n-1}(s)\,ds .
```

<a id="r-8-10"></a>

**Proposición 8.10.** Para un estado $`\rho`$, sea $`w_n(t)=\mathrm{Tr}\,S_n(t)\rho`$. Entonces
$`w_n(t)\ge0`$, $`\sum_nw_n(t)=1`$, y $`w_0`$ es no creciente con
$`w_0'(t)=-\mathrm{Tr}\big(A\,S_0(t)\rho\big)`$. Para un estado puro
$`\rho=|\psi\rangle\langle\psi|`$, $`w_0(t)=\|e^{-itH_{\mathrm{eff}}}\psi\|^2`$.

*Demostración.* $`S_n(t)`$ es CP, luego positiva, y $`\sum_nw_n=\mathrm{Tr}\,\rho(t)=1`$. Con
$`X=S_0(t)\rho\ge0`$,
$`w_0'=\mathrm{Tr}(-iH_{\mathrm{eff}}X+iXH_{\mathrm{eff}}^*)=-\mathrm{Tr}(AX)`$, y
$`\mathrm{Tr}(AX)=\mathrm{Tr}(X^{1/2}AX^{1/2})\ge0`$. La última afirmación es
$`\mathrm{Tr}|\phi\rangle\langle\phi|=\|\phi\|^2`$. $`\square`$

> **Significado físico** (interpretación física, no un teorema de este curso). Supongamos que un detector registra cada salto: cada fotón emitido, cada
> átomo perdido. Entonces $`w_n(t)`$ es la probabilidad de exactamente $`n`$ clics hasta el instante
> $`t`$, $`S_n(t)\rho/w_n(t)`$ el estado condicionado a $`n`$ clics, y $`S_0(t)\rho/w_0(t)`$ el estado
> condicionado a **ningún clic**: un resultado nulo también cambia el estado, a través del
> $`H_{\mathrm{eff}}`$ no hermítico. Esta interpretación como proceso de conteo [E. B. Davies,
> *Quantum Theory of Open Systems*, 1976; M. D. Srinivas y E. B. Davies, *Optica Acta* 28, 981
> (1981)] es la base del método de las trayectorias cuánticas (función de onda de Monte Carlo) [J.
> Dalibard, Y. Castin, K. Mølmer, *Phys. Rev. Lett.* 68, 580 (1992)]; aquí solo usamos su parte
> matemática.

En el caso de las pérdidas, los términos de Dyson son las componentes de sector de $`\rho(t)`$ (en el
lenguaje del recuadro: contar los clics es lo mismo que medir la carga). La **actualización de
Lüders** de un estado $`\rho`$ por una proyección $`P`$ con $`\mathrm{Tr}(P\rho)>0`$ es el estado
$`P\rho P/\mathrm{Tr}(P\rho P)`$ (nótese que $`\mathrm{Tr}(P\rho P)=\mathrm{Tr}(P\rho)`$).

*En la biblioteca:* la actualización de Lüders es `states.lueders_update` — [tutorial, sección 8.4](../../tutorial/es/08_dynamics.md#s-8-4).

<a id="r-8-11"></a>

**Proposición 8.11.** Etiquetemos los sectores por $`N=0,\dots,C-1`$, supongamos que no hay saltos
dentro de los sectores, que todo salto reduce la carga en uno ($`K^{dc}_{m,l}=0`$ salvo si $`d=c-1`$) y
que $`\rho=P_{N_0}\rho P_{N_0}`$. Entonces $`S_n(t)\rho=P_{N_0-n}\,\rho(t)\,P_{N_0-n}`$ para
$`0\le n\le N_0`$ y $`S_n(t)\rho=0`$ para $`n>N_0`$. En particular, $`p_{N_0}(t)=w_0(t)>0`$, y la
actualización de Lüders de
$`\rho(t)`$ por $`P_{N_0}`$ ("no se perdió ninguna partícula") es el estado sin saltos:
```math
\frac{P_{N_0}\rho(t)P_{N_0}}{p_{N_0}(t)}=\frac{e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*}}{w_0(t)},
\qquad p_{N_0}(t)=w_0(t).
```

*Demostración.* $`S_0(s)`$ lleva operadores con soporte en el sector $`N`$ a operadores con soporte en
el sector $`N`$, ya que $`H_{\mathrm{eff}}\in\mathcal M`$; $`\Phi_{\mathrm{tot}}`$ lleva el sector $`N\ge1`$
al sector $`N-1`$ y anula el sector $`0`$. Por inducción, $`S_n(s)\rho`$ tiene soporte en el sector
$`N_0-n`$ para todo $`s`$, y es $`0`$ para $`n>N_0`$. Comprimir $`\rho(t)=\sum_nS_n(t)\rho`$ con $`P_{N_0-n}`$
conserva exactamente el $`n`$-ésimo término. Por último, $`w_0(t)>0`$, porque $`e^{-itH_{\mathrm{eff}}}`$
es invertible (lema 8.1), así que $`S_0(t)\rho`$ es una matriz positiva no nula. $`\square`$

*Cálculo:* la proposición 8.11 para $`n=0`$ y $`n=1`$, para dos modos que pierden partículas con tasas distintas — [tutorial, sección 8.4](../../tutorial/es/08_dynamics.md#s-8-4).

<a id="s-8-5"></a>

## 8.5 Estados estacionarios y relajación

Un estado $`\sigma`$ es **estacionario** si $`e^{t\mathcal L}\sigma=\sigma`$ para todo $`t\ge0`$, o,
equivalentemente, si $`\mathcal L(\sigma)=0`$ (derívese en $`t=0`$; recíprocamente, úsese la serie).

<a id="r-8-12"></a>

**Proposición 8.12.** Todo semigrupo dinámico cuántico sobre $`\mathcal M`$ tiene un estado
estacionario.

*Demostración.* Fijemos un estado $`\rho`$ y sea $`\sigma_\tau=\frac1\tau\int_0^\tau T_s\rho\,ds`$. Es
un estado por el lema A.4(d) (los estados forman un conjunto convexo cerrado). Por el lema
A.4(b),(e),
$`T_t\sigma_\tau-\sigma_\tau=\frac1\tau\big(\int_\tau^{\tau+t}-\int_0^t\big)T_s\rho\,ds`$, cuya norma
traza es a lo sumo $`2t/\tau`$ por el lema A.4(c), ya que $`\|T_s\rho\|_1=1`$. Los estados forman un conjunto cerrado y
acotado en un espacio de dimensión finita, así que alguna sucesión $`\sigma_{\tau_m}`$,
$`\tau_m\to\infty`$, converge a un estado $`\sigma`$ (Bolzano–Weierstrass), y
$`T_t\sigma-\sigma=\lim_m(T_t\sigma_{\tau_m}-\sigma_{\tau_m})=0`$. $`\square`$

<a id="r-8-13"></a>

**Corolario 8.13.** Si no hay saltos entre sectores (como en la proposición 8.9(1)), los estados
estacionarios generan un espacio de dimensión al menos $`C`$: cada sector tiene su propio estado
estacionario.

*Demostración.* $`\mathcal L`$ lleva $`P_c\mathcal MP_c\cong M_{k_c}`$ en sí mismo ($`H`$ y los $`L_j`$ están
en $`\mathcal M`$, y $`\Phi_m(X)=\sum_lK^{cc}_{m,l}X_cK^{cc*}_{m,l}`$ y $`\Phi_m^*(1)`$ son diagonales por
bloques), luego también lo hace cada $`e^{t\mathcal L}`$. Su restricción es CP (lema 8.4(a)) y TP, así
que es un semigrupo dinámico cuántico sobre $`M_{k_c}`$, con un estado estacionario $`\sigma_c`$ por la
proposición 8.12; como elemento de $`\mathcal M`$, $`\sigma_c`$ es un estado estacionario de
$`e^{t\mathcal L}`$. Los $`\sigma_c`$ están en bloques distintos, luego son linealmente independientes.
$`\square`$

El teorema D.6 (el caso de dimensión finita de un criterio de A. Frigerio) da una condición
suficiente de unicidad: si un generador GKSL sobre $`M_k`$ con datos $`(H;L_1,\dots,L_r)`$ tiene un
estado estacionario fiel $`\sigma>0`$ y las únicas matrices que conmutan con $`H`$ y con todos los $`L_j`$,
$`L_j^*`$ son los múltiplos de $`1`$, entonces $`\sigma`$ es el único estado estacionario. Sin saltos entre
sectores, la restricción de $`\mathcal L`$ a $`P_c\mathcal MP_c\cong M_{k_c}`$ es el generador GKSL con
datos $`(H_c;\sqrt{g_1}L_{1,c},\sqrt{g_2}L_{2,c},\dots)`$, así que el teorema se aplica sector a sector.
El ejercicio 4 muestra que la unicidad puede fallar incluso en un factor: allí falla la segunda
hipótesis.

*Cálculo:* el corolario 8.13 y las hipótesis del teorema D.6 para saltos aleatorios dentro de los sectores de $`M_2\oplus M_3`$ — [tutorial, sección 8.5](../../tutorial/es/08_dynamics.md#s-8-5).

<a id="r-8-14"></a>

**Proposición 8.14 (emisión espontánea).** Sobre $`M_2`$ con base $`|g\rangle,|e\rangle`$, sean
$`H=\omega|e\rangle\langle e|`$, un salto $`L=|g\rangle\langle e|`$ y una tasa $`\gamma>0`$. Entonces
```math
\rho_{ee}(t)=e^{-\gamma t}\rho_{ee}(0),\qquad\rho_{ge}(t)=e^{(i\omega-\gamma/2)t}\rho_{ge}(0),
\qquad e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}},
```
donde $`\mathcal U_t(\rho)=e^{-itH}\rho e^{itH}`$ y $`\mathcal A_p`$ es el canal de amortiguamiento de
amplitud de parámetro $`p`$ (sección 6.4). El espectro de $`\mathcal L`$ es
$`\{0,-\gamma,-\gamma/2\pm i\omega\}`$, $`|g\rangle\langle g|`$ es el único estado estacionario, y
$`\tfrac12\|\rho(t)-|g\rangle\langle g|\|_1=\big(\rho_{ee}(0)^2e^{-2\gamma t}+|\rho_{ge}(0)|^2e^{-\gamma t}\big)^{1/2}`$.

*Demostración.* Aquí $`L^*L=|e\rangle\langle e|`$ y $`L\rho L^*=\rho_{ee}|g\rangle\langle g|`$, luego,
entrada a entrada, $`\rho_{ee}'=-\gamma\rho_{ee}`$, $`\rho_{gg}'=\gamma\rho_{ee}`$ y
$`\rho_{ge}'=(i\omega-\gamma/2)\rho_{ge}`$, lo que da la solución. En la base
$`(\rho_{gg},\rho_{ge},\rho_{eg},\rho_{ee})`$, la matriz de $`\mathcal L`$ es triangular superior con
diagonal $`(0,-\gamma/2+i\omega,-\gamma/2-i\omega,-\gamma)`$. $`\mathcal L(\sigma)=0`$ obliga a
$`\sigma_{ee}=\sigma_{ge}=0`$. El canal de amortiguamiento de amplitud multiplica $`\rho_{ee}`$ por
$`1-p`$ y $`\rho_{ge}`$ por $`\sqrt{1-p}`$ (sus operadores de Kraus son
$`|g\rangle\langle g|+\sqrt{1-p}\,|e\rangle\langle e|`$ y $`\sqrt p\,|g\rangle\langle e|`$), y
$`\mathcal U_t`$ multiplica $`\rho_{ge}`$ por $`e^{i\omega t}`$; ambos son TP, lo que fija $`\rho_{gg}`$.
Por último, $`\rho(t)-|g\rangle\langle g|`$ es autoadjunta, de traza nula y con diagonal
$`(-\rho_{ee}(t),\rho_{ee}(t))`$, así que sus valores propios son
$`\pm(\rho_{ee}(t)^2+|\rho_{ge}(t)|^2)^{1/2}`$. $`\square`$

La coherencia decae a la mitad de la tasa de la población: con $`T_1=1/\gamma`$ y $`T_2=2/\gamma`$, esta es
la relación $`\frac1{T_2}=\frac1{2T_1}+\frac1{T_\varphi}`$ de la sección 6.4 sin el término de desfase. La
tasa de decaimiento no nula más lenta, la **brecha espectral**
$`\min\{-\mathrm{Re}\,\mu:\mu\ne0\text{ un valor propio de }\mathcal L\}=\gamma/2`$, es aquí la tasa
de relajación: por la última fórmula, la distancia a $`|g\rangle\langle g|`$ es a lo sumo
$`e^{-\gamma t/2}`$ (pues $`|\rho_{ge}|^2\le\rho_{gg}\rho_{ee}`$ por positividad,
$`\rho_{ee}^2+|\rho_{ge}|^2\le\rho_{ee}\le1`$), y $`e^{\gamma t/2}`$ veces ella tiende a $`|\rho_{ge}(0)|`$
cuando $`t\to\infty`$.

*Cálculo:* el espectro de $`\mathcal L`$, la factorización $`e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}}`$ y la distancia de traza de la proposición 8.14, con `dynamics.lindblad_superoperator` y `dynamics.lindblad_channel` — [tutorial, sección 8.5](../../tutorial/es/08_dynamics.md#s-8-5).

<a id="s-8-6"></a>

## 8.6 Con qué comparar

*Cálculo:* la comparación con otro software (QuTiP, dynamiqs, QuantumOptics.jl, SciPy) — [tutorial, sección 8.6](../../tutorial/es/08_dynamics.md#s-8-6).

<a id="x-8"></a>

## Ejercicios

1. (Imagen de Heisenberg) Demuestre que
   $`\mathrm{Tr}(\mathcal L(\rho)A)=\mathrm{Tr}(\rho\,\mathcal L^*(A))`$ para el generador
   de la sección 8.3, con
   $`\mathcal L^*(A)=i[H,A]+\sum_jg_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})+\sum_mg_m(\Phi_m^*(A)-\frac12\{\Phi_m^*(1),A\})`$,
   que $`\mathcal L^*(1)=0`$, y que $`\mathcal L^*(P_c)=0`$ si no hay saltos entre sectores.
2. (Punto excepcional) En el sector $`N=1`$ de la sección 8.4,
   $`H_{\mathrm{eff}}=\begin{pmatrix}-i\gamma_1/2&J\\J&-i\gamma_2/2\end{pmatrix}`$. Halle sus valores
   propios y demuestre que para $`J=|\gamma_1-\gamma_2|/4>0`$ no es diagonalizable.
3. Supongamos que $`N_0`$ bosones en $`m`$ modos se pierden cada uno con tasa $`\gamma`$ (saltos
   $`a_1,\dots,a_m`$ restringidos a $`N\to N-1`$). Demuestre que $`p_N(t)=\binom{N_0}Nq^N(1-q)^{N_0-N}`$
   con $`q=e^{-\gamma t}`$, para todo hamiltoniano de $`\mathcal M`$.
4. (Desfase puro) Sobre $`M_2`$, sean $`H=\omega|e\rangle\langle e|`$ y
   $`L=\sigma_z=\mathrm{diag}(1,-1)`$ con tasa $`\gamma`$. Demuestre que las poblaciones son
   constantes, que $`\rho_{ge}(t)=e^{(i\omega-2\gamma)t}\rho_{ge}(0)`$, y que los estados
   estacionarios generan el espacio bidimensional de las matrices diagonales aunque $`C=1`$.

<details>
<summary>Soluciones</summary>

1. Por ciclicidad, $`\mathrm{Tr}(L\rho L^*A)=\mathrm{Tr}(\rho L^*AL)`$,
   $`\mathrm{Tr}(\{L^*L,\rho\}A)=\mathrm{Tr}(\rho\{L^*L,A\})`$ y
   $`\mathrm{Tr}(-i[H,\rho]A)=\mathrm{Tr}(\rho\,i[H,A])`$; para los saltos entre sectores,
   $`\mathrm{Tr}(\Phi(\rho)A)=\mathrm{Tr}(\rho\Phi^*(A))`$ es la definición de $`\Phi^*`$
   (sección 6.3). Entonces
   $`\mathcal L^*(1)=\sum_jg_j(L_j^*L_j-L_j^*L_j)+\sum_mg_m(\Phi_m^*(1)-\Phi_m^*(1))=0`$, que es una
   vez más la preservación de la traza. Sin saltos entre sectores, $`P_c`$ conmuta con $`H`$ y con todos
   los $`L_j`$, luego $`L_j^*P_cL_j=L_j^*L_jP_c=\frac12\{L_j^*L_j,P_c\}`$, cada $`\Phi_m`$ tiene solo los bloques $`K^{cc}_{m,l}`$, luego
   $`\Phi_m^*(P_c)=\Phi_m^*(1)P_c=\frac12\{\Phi_m^*(1),P_c\}`$, y $`\mathcal L^*(P_c)=0`$;
   entonces $`p_c'=\mathrm{Tr}(\rho\mathcal L^*(P_c))=0`$, una segunda demostración de la
   proposición 8.9(1).

2. $`H_{\mathrm{eff}}=-\frac i4(\gamma_1+\gamma_2)1+\begin{pmatrix}-i\delta&J\\J&i\delta\end{pmatrix}`$
   con $`\delta=(\gamma_1-\gamma_2)/4`$; la segunda matriz tiene traza $`0`$ y determinante
   $`\delta^2-J^2`$, así que los valores propios son
   $`-\frac i4(\gamma_1+\gamma_2)\pm\sqrt{J^2-\delta^2}`$: dos tasas de decaimiento para $`J<|\delta|`$,
   una tasa y dos frecuencias para $`J>|\delta|`$. Para $`J=|\delta|>0`$, la segunda matriz $`Y`$ cumple
   $`Y^2=(J^2-\delta^2)1=0`$ pero $`Y\ne0`$, de modo que $`H_{\mathrm{eff}}`$ tiene el único valor propio
   $`-\frac i4(\gamma_1+\gamma_2)`$ y no es múltiplo de $`1`$: no es diagonalizable.

3. En el sector $`N`$, $`\sum_ja_j^*a_j=N1`$, y el único par de sectores es $`(N-1,N)`$, así que la
   proposición 8.9(2) da $`p_N'=-\gamma Np_N+\gamma(N+1)p_{N+1}`$. Con $`q'=-\gamma q`$, la derivada del
   término binomial es $`-\gamma Np_N+\gamma(N_0-N)\binom{N_0}Nq^{N+1}(1-q)^{N_0-N-1}`$, y
   $`(N+1)\binom{N_0}{N+1}=(N_0-N)\binom{N_0}N`$ muestra que el segundo término es
   $`\gamma(N+1)p_{N+1}`$. En $`t=0`$ es la masa puntual en $`N_0`$; la unicidad es el lema 8.1.

4. $`L\rho L^*-\frac12\{L^*L,\rho\}=\sigma_z\rho\sigma_z-\rho`$ multiplica las entradas fuera de la
   diagonal por $`-2`$ y anula la diagonal, y $`-i[H,\rho]_{ge}=i\omega\rho_{ge}`$. Así,
   $`\mathcal L(\sigma)=0`$ si y solo si $`\sigma_{ge}=0`$ (para $`\gamma>0`$), y el espectro es
   $`\{0,0,-2\gamma\pm i\omega\}`$. Con $`U_t=e^{-itH}`$,
   $`e^{t\mathcal L}(\rho)=\frac{1+e^{-2\gamma t}}2U_t\rho U_t^*+\frac{1-e^{-2\gamma t}}2U_t\sigma_z\rho\sigma_zU_t^*`$:
   como $`\sigma_z`$ conmuta con $`U_t`$ y $`\sigma_z^2=1`$, derivando se ve que el segundo miembro $`R(t)`$
   cumple $`R'=-i[H,R]+\gamma(\sigma_zR\sigma_z-R)=\mathcal L(R)`$ y $`R(0)=\rho`$, así que se aplica el
   lema 8.1. Para $`t>0`$ es un canal unital cuyos puntos fijos son las matrices que conmutan con sus
   operadores de Kraus (teorema 6.10), es decir, las diagonales.

</details>

*Cálculo:* los ejercicios 2 y 4, numéricamente — [tutorial, Ejercicios](../../tutorial/es/08_dynamics.md#x-8).

**Siguiente:** el [capítulo 9](09_spin_chains.md) trata las cadenas de espines: sectores de
magnetización y de momento, evolución temporal de Lanczos y de Krylov para sectores grandes, y
desorden.
