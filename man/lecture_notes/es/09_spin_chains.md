[English](../09_spin_chains.md) · **Español** · [Français](../fr/09_spin_chains.md) · [Deutsch](../de/09_spin_chains.md) · [中文](../zh/09_spin_chains.md) · [日本語](../ja/09_spin_chains.md) · [Русский](../ru/09_spin_chains.md)

# 9. Cadenas de espines

Una cadena de $`L`$ espines $`\frac12`$ vive en $`(\mathbb C^2)^{\otimes L}`$, y el álgebra lineal densa en dimensión
$`2^L`$ se vuelve costosa más allá de $`L\approx14`$ (una matriz densa tiene $`4^L`$ entradas). Las simetrías desplazan ese límite: la magnetización divide el espacio en
$`L+1`$ sectores (capítulo 4), las traslaciones dividen de nuevo cada sector y $`SU(2)`$ reduce los bloques una vez
más. Este capítulo construye estos sectores para la cadena XXZ, demuestra lo que dicen sobre los estados
fundamentales, contrasta la biblioteca con los fermiones libres y el ansatz de Bethe, y va más allá de las
matrices densas (métodos de Krylov) y de las cadenas limpias (desorden).

**Aprenderá** por qué el álgebra de una cadena es $`\bigoplus_NM_{\binom LN}`$ y por qué la torre de estas
álgebras tiene como diagrama de Bratteli el triángulo de Pascal e índice $`4\cos^2\frac\pi{2L+4}`$, cómo se
cuentan los sectores de momento, la regla de signos de Marshall y sus consecuencias para el momento y el espín
de los estados fundamentales (Lieb–Mattis), la solución de Jordan–Wigner de la cadena XX, las ecuaciones de
Bethe de dos magnones, por qué Lanczos converge monótonamente, cuán preciso es un paso de Krylov en el tiempo, y cómo medir la estadística de niveles y el
entrelazamiento de fermiones libres. Las comprobaciones numéricas de este capítulo están en el
[Tutorial 9](../../tutorial/es/09_spin_chains.md).

En todo el capítulo, los sitios son $`0,\dots,L-1`$, la base de $`(\mathbb C^2)^{\otimes L}`$ está formada por las
configuraciones $`|b\rangle=|b_0\dots b_{L-1}\rangle`$ con $`b_i=1`$ para el espín hacia arriba (la convención de la
biblioteca: el bit $`i`$ de una etiqueta es el sitio $`i`$), $`s^a_i`$ ($`a=x,y,z`$) son las matrices de espín
$`\frac12\sigma_a`$ en el sitio $`i`$, $`s^\pm_i=s^x_i\pm is^y_i`$, $`n_i=s^z_i+\frac12`$, $`S^a=\sum_is^a_i`$, y
$`N=S^z+\frac L2`$ es el número de espines hacia arriba. Una cadena es **abierta** (enlaces $`(i,i+1)`$, $`i<L-1`$) o
**periódica** (también el enlace $`(L-1,0)`$; entonces $`L\ge3`$ y los índices se toman mód $`L`$). $`J>0`$
(antiferromagnética), $`\hbar=1`$, y $`\log`$ es el logaritmo natural.

<a id="s-9-1"></a>

## 9.1 La cadena XXZ y sus sectores $`U(1)`$

La **cadena XXZ** con anisotropía $`\Delta\in\mathbb R`$ y campos $`h_i\in\mathbb R`$ es
```math
H=J\sum_{(i,j)}\Big[\tfrac12\big(s^+_is^-_j+s^-_is^+_j\big)+\Delta\,s^z_is^z_j\Big]+\sum_ih_is^z_i ,
```
con la suma sobre los enlaces; $`\frac12(s^+_is^-_j+s^-_is^+_j)=s^x_is^x_j+s^y_is^y_j`$ desplaza un espín hacia
arriba a lo largo del enlace. Para $`\Delta=1`$ es la cadena de Heisenberg $`J\sum\mathbf s_i\cdot\mathbf s_j`$; para
$`\Delta=0`$, la cadena XX.

<a id="r-9-1"></a>

**Proposición 9.1.** $`[H,S^z]=0`$. Por tanto, $`H\in U(1)'\cong\bigoplus_{N=0}^LM_{\binom LN}(\mathbb C)`$, donde el
bloque $`N`$ actúa sobre el subespacio generado por las configuraciones con $`N`$ espines hacia arriba.

*Demostración.* Los términos diagonales conmutan con la matriz diagonal $`S^z`$. Para $`i\neq j`$, $`s^+_is^-_j`$ conmuta
con $`s^z_l`$, $`l\ne i,j`$, y $`[s^z_i+s^z_j,s^+_is^-_j]=[s^z_i,s^+_i]s^-_j+s^+_i[s^z_j,s^-_j]=s^+_is^-_j-s^+_is^-_j=0`$ por
$`[s^z,s^\pm]=\pm s^\pm`$ (sección 4.5); para $`s^-_is^+_j`$ se toman adjuntos. Por tanto, $`H`$ conmuta
con $`N=S^z+\frac L2`$, que tiene los valores propios enteros $`0,\dots,L`$, con espacios propios generados por las
configuraciones con $`N`$ espines hacia arriba, de dimensión $`\binom LN`$; el resto es la proposición 4.10(2) con
$`Q=N`$. $`\square`$

*Cálculo:* la proposición 9.1, y los bloques $`H_N`$ que `SpinChain.xxz` construye a partir de las configuraciones de cada sector (`chain.labels[N]`), frente a la matriz completa $`2^L\times2^L`$ — [tutorial, sección 9.1](../../tutorial/es/09_spin_chains.md#s-9-1).

**La torre de cadenas.** Sea $`\mathcal A_L=\bigoplus_NM_{\binom LN}`$, que actúa sobre $`L`$ sitios. Añadir un sitio,
$`x\mapsto x\otimes1_2`$, aplica $`\mathcal A_L`$ en $`\mathcal A_{L+1}`$, porque $`x\otimes1`$ conmuta con
$`S^z\otimes1+1\otimes s^z`$. Esto da una torre $`\mathcal A_0\subset\mathcal A_1\subset\mathcal A_2\subset\cdots`$.

<a id="r-9-2"></a>

**Proposición 9.2.** La inclusión $`\mathcal A_L\subset\mathcal A_{L+1}`$ tiene la matriz de inclusión
$`\Lambda_L`$ de tamaño $`(L+2)\times(L+1)`$ con $`\lambda_{N,N}=\lambda_{N,N-1}=1`$ y todas las demás entradas $`0`$; los
tamaños de los bloques obedecen la regla de Pascal $`\binom{L+1}N=\binom LN+\binom L{N-1}`$, el diagrama de Bratteli
de la torre es el triángulo de Pascal y el grafo bipartito de $`\Lambda_L`$ es un camino con $`2L+3`$ vértices.

*Demostración.* Las configuraciones con $`N`$ espines hacia arriba en $`L+1`$ sitios son las que tienen $`N`$ espines
hacia arriba en los primeros $`L`$ sitios y el nuevo espín hacia abajo, y las que tienen $`N-1`$ y el nuevo espín
hacia arriba. Sobre estos dos subespacios, $`x\otimes1`$ actúa como $`x_N`$ y $`x_{N-1}`$, de modo que el bloque $`N`$ de
$`\mathcal A_{L+1}`$ contiene los bloques $`N`$ y $`N-1`$ de $`\mathcal A_L`$ una vez cada uno; las dimensiones suman según
la regla de Pascal. En el grafo bipartito, el sector $`N`$ de $`\mathcal A_{L+1}`$ está unido a $`N-1`$ y $`N`$ de
$`\mathcal A_L`$, lo que da el camino $`0'-0-1'-1-\cdots-L'-L-(L+1)'`$ (con primas para $`\mathcal A_{L+1}`$). $`\square`$

Un camino desde la cima del diagrama de Bratteli hasta el sector $`N`$ del nivel $`L`$ es una configuración $`b`$; hay
$`\binom LN`$ de ellos, el tamaño del bloque. Para una inclusión conexa con matriz de inclusión $`\Lambda`$
**definimos** el **índice** como $`\|\Lambda\|^2`$ (norma de operador); para $`M_k\otimes1\subset M_k\otimes M_m`$ la matriz
de inclusión es la matriz $`1\times1`$ $`(m)`$, así que el índice es el $`m^2`$ de la sección 7.2. `Inclusion.index` lo
calcula. (El nombre se toma de la teoría del índice de V. Jones, *Invent. Math.* 72 (1983) 1; en este curso, el
índice de una inclusión de álgebras de dimensión finita es esta definición y nada más.)

*En la biblioteca:* el índice es `Inclusion.index` — [tutorial, sección 9.1](../../tutorial/es/09_spin_chains.md#s-9-1).

<a id="r-9-3"></a>

**Lema 9.3.** La matriz de adyacencia $`A_n`$ del camino con vértices $`1,\dots,n`$ tiene los valores propios
$`2\cos\frac{\pi k}{n+1}`$, $`k=1,\dots,n`$, con vectores propios $`v_k(j)=\sin\frac{\pi jk}{n+1}`$; por tanto,
$`\|A_n\|=2\cos\frac\pi{n+1}`$.

*Demostración.* Pongamos $`v_k(0)=v_k(n+1)=0`$. Entonces $`(A_nv_k)(j)=v_k(j-1)+v_k(j+1)=2\cos\frac{\pi k}{n+1}\,v_k(j)`$
por $`\sin(a-b)+\sin(a+b)=2\sin a\cos b`$, y $`v_k\ne0`$ porque $`v_k(1)=\sin\frac{\pi k}{n+1}\neq0`$. Los $`n`$ valores
propios son distintos ($`\cos`$ es estrictamente decreciente en $`(0,\pi)`$), luego son todos, y la norma de la matriz
simétrica $`A_n`$ es el mayor módulo, alcanzado en $`k=1`$ y $`k=n`$. $`\square`$

<a id="r-9-4"></a>

**Corolario 9.4.** $`[\mathcal A_{L+1}:\mathcal A_L]=\|\Lambda_L\|^2=4\cos^2\frac\pi{2L+4}`$, que crece hacia $`4`$.

*Demostración.* Ordenando los vértices a lo largo del camino de la proposición 9.2, la matriz de adyacencia
$`B=\begin{pmatrix}0&\Lambda_L\\\Lambda_L^T&0\end{pmatrix}`$ se convierte en $`A_{2L+3}`$ (una semejanza por permutación). Por el
lema 1.1, $`\|B\|^2=\|B^2\|`$ y $`\|\Lambda_L\Lambda_L^T\|=\|\Lambda_L^T\Lambda_L\|=\|\Lambda_L\|^2`$, y $`B^2=\Lambda_L\Lambda_L^T\oplus\Lambda_L^T\Lambda_L`$, así que
$`\|B\|=\|\Lambda_L\|`$; el lema 9.3 con $`n=2L+3`$ da el valor. $`\square`$

*Cálculo:* la proposición 9.2, el lema 9.3 y el corolario 9.4 para $`L\le8`$, con `inclusions.chain_tower` e `inclusions.bratteli_diagram` — [tutorial, sección 9.1](../../tutorial/es/09_spin_chains.md#s-9-1).

> **Significado físico** (interpretación física, no un teorema de este curso). $`S^z`$ se conserva siempre que la interacción sea invariante bajo rotaciones alrededor
> del eje $`z`$; $`N`$ cuenta los espines hacia arriba, o las partículas tras la transformación de Jordan–Wigner de
> la sección 9.4. En la torre, cada sector de carga de la cadena más larga se alimenta de dos sectores de la
> más corta.

<a id="s-9-2"></a>

## 9.2 Traslaciones, momento y la regla de signos de Marshall

En una cadena periódica, sea $`T`$ la traslación $`T|b_0b_1\dots b_{L-1}\rangle=|b_{L-1}b_0\dots b_{L-2}\rangle`$, que lleva
el sitio $`i`$ al sitio $`i+1`$. $`T`$ es una matriz de permutación, luego unitaria, $`T^L=1`$ y $`Ts^a_iT^*=s^a_{i+1}`$. Por
tanto, $`TS^zT^*=S^z`$ y $`THT^*=H`$ para un campo uniforme $`h_i=h`$: $`T`$ lleva enlaces a enlaces. Por tanto, sus
valores propios son de la forma $`e^{ik}`$ con **momentos** $`k=2\pi m/L`$, $`m=0,\dots,L-1`$.

<a id="r-9-5"></a>

**Proposición 9.5 (sectores de momento).** Fijemos $`N`$. Para una configuración $`a`$ con $`N`$ espines hacia arriba,
sea $`R_a`$ el tamaño de su órbita $`\{T^ja\}`$. (a) $`R_a`$ divide a $`L`$; para $`k=2\pi m/L`$ con $`mR_a\equiv0\pmod L`$, el
vector $`|a,k\rangle=\sum_{j=0}^{R_a-1}e^{-ikj}T^j|a\rangle`$ es no nulo y $`T|a,k\rangle=e^{ik}|a,k\rangle`$, y estos $`R_a`$
vectores forman una base ortogonal del subespacio generado por la órbita. (b) El valor propio $`e^{2\pi im/L}`$ de
$`T`$ sobre el sector $`N`$ tiene multiplicidad $`d_{N,m}=\#\{\text{órbitas }:mR_a\equiv0\pmod L\}`$, y todo operador que
conmute con $`S^z`$ y $`T`$ —por ejemplo, la cadena XXZ invariante por traslaciones— es diagonal por bloques, en
$`\bigoplus_N\bigoplus_mM_{d_{N,m}}`$.

*Demostración.* (a) El estabilizador de $`a`$ en el grupo cíclico $`\{T^j\}\cong\mathbb Z_L`$ es un subgrupo, y $`R_a`$ es
su índice, un divisor de $`L`$; los vectores $`T^j|a\rangle`$, $`0\le j<R_a`$, son ortonormales y $`T^{R_a}|a\rangle=|a\rangle`$. Si
$`e^{ikR_a}=1`$, entonces $`T|a,k\rangle=e^{ik}\sum_{j=1}^{R_a}e^{-ikj}T^j|a\rangle=e^{ik}|a,k\rangle`$, ya que el término $`j=R_a`$
es igual al término $`j=0`$; el coeficiente de $`|a\rangle`$ es $`1`$. Los $`m\in\{0,\dots,L-1\}`$ admisibles son los
múltiplos de $`L/R_a`$, $`R_a`$ en total, con valores propios distintos del unitario $`T`$, de modo que los vectores son
ortogonales y generan el subespacio de dimensión $`R_a`$ generado por la órbita. (b) El sector $`N`$ es la suma
ortogonal de los subespacios generados por sus órbitas, cada uno invariante por $`T`$, y por (a) el valor propio
$`e^{2\pi im/L}`$ aparece una vez en el subespacio generado por la órbita de $`a`$ si $`mR_a\equiv0\pmod L`$, y ninguna
en otro caso. Un $`X`$ que conmuta con $`S^z`$ y $`T`$ conmuta con sus
proyecciones espectrales $`P_N`$, $`Q_m`$ (teorema 1.5(d)), así que $`X=\sum_{N,m}P_NQ_mXP_NQ_m`$. $`\square`$

*Cálculo:* la proposición 9.5(b) para $`L=12`$, y los bloques de momento que construye `SpinChain.xxz_momentum` — [tutorial, sección 9.2](../../tutorial/es/09_spin_chains.md#s-9-2).

<a id="r-9-6"></a>

**Lema 9.6 (Perron–Frobenius, caso simétrico).** Sea $`M`$ una matriz real simétrica con $`M_{ij}\ge0`$ para $`i\ne j`$
cuyo grafo ($`i\sim j`$ si y solo si $`M_{ij}>0`$) es conexo. Entonces su mayor valor propio $`\lambda`$ es simple, con un
vector propio cuyas entradas son todas $`>0`$.

*Demostración.* $`\lambda=\max_{\|x\|=1}\langle x,Mx\rangle`$, y el máximo se alcanza exactamente en los vectores unitarios
del espacio propio de $`\lambda`$ (desarrollando en una base ortonormal de vectores propios). Para tal $`x`$, el vector
$`|x|`$ de los módulos cumple $`\langle|x|,M|x|\rangle\ge\mathrm{Re}\langle x,Mx\rangle=\lambda`$, porque los términos
diagonales coinciden y $`M_{ij}|x_i||x_j|\ge\mathrm{Re}\,M_{ij}\bar x_ix_j`$ para $`i\ne j`$; así que $`|x|`$ es un vector
propio. Si $`|x|_i=0`$, entonces $`0=(\lambda-M_{ii})|x|_i=\sum_{j\ne i}M_{ij}|x|_j`$ obliga a $`|x|_j=0`$ para todos los
vecinos $`j`$ de $`i`$, y por conexión $`x=0`$. Así, ningún vector propio no nulo tiene una entrada nula; pero un espacio
propio de dimensión $`\ge2`$ contiene un vector no nulo con $`x_1=0`$. Por tanto, $`\lambda`$ es simple y $`|x|>0`$ es un
vector propio. $`\square`$

Demostramos la regla de signos para cualquier disposición bipartita de los acoplamientos; la cadena es un caso
particular.

<a id="r-9-7"></a>

**Teorema 9.7 (regla de signos de Marshall).** Sea $`A\sqcup B`$ una partición de los sitios $`\{0,\dots,L-1\}`$ en dos
conjuntos, y sea
```math
H=\sum_{i<j}\Big[\tfrac{J_{ij}}2\big(s^+_is^-_j+s^-_is^+_j\big)+K_{ij}\,s^z_is^z_j\Big]+\sum_ih_is^z_i
```
con $`K_{ij},h_i`$ reales y $`J_{ij}`$ reales tales que $`J_{ij}\ge0`$ si $`i,j`$ están en conjuntos distintos, $`J_{ij}\le0`$ si
están en el mismo conjunto, y el grafo $`G`$ sobre los sitios con una arista $`i\sim j`$ si y solo si $`J_{ij}\ne0`$ es
conexo. Sea $`N_B(b)`$ el número de espines hacia arriba de $`b`$ en $`B`$. En cada sector $`N`$, el menor valor propio de
$`H_N`$ es simple, y su vector propio puede elegirse con $`(-1)^{N_B(b)}\psi(b)>0`$ para toda configuración $`b`$ con $`N`$
espines hacia arriba.

En particular, esto vale para la cadena XXZ con $`J>0`$ y $`\Delta,h_i`$ arbitrarios si la cadena es abierta, o
periódica con $`L`$ par, con $`A`$ los sitios pares y $`B`$ los impares: entonces $`J_{ij}=J`$, $`K_{ij}=J\Delta`$ sobre los
enlaces y $`0`$ en otro caso, todo enlace une un sitio par con uno impar (en el anillo porque $`L`$ es par), y $`G`$
contiene el camino $`0-1-\cdots-(L-1)`$.

*Demostración.* Sea $`V|b\rangle=(-1)^{N_B(b)}|b\rangle`$; $`V`$ es diagonal y $`V=V^*=V^{-1}`$. El operador $`s^+_is^-_j`$ lleva
un espín hacia arriba de $`j`$ a $`i`$; cambia $`N_B`$ en $`\pm1`$ si $`i,j`$ están en conjuntos distintos y lo deja igual
en otro caso, así que $`Vs^+_is^-_jV=\mp s^+_is^-_j`$ según el caso, y lo mismo para $`s^-_is^+_j`$. Dos configuraciones
distintas $`b,b'`$ están unidas por un tal término si y solo si $`b'`$ se obtiene de $`b`$ intercambiando los espines
opuestos de dos sitios $`i,j`$, y entonces $`i,j`$ quedan determinados por $`b,b'`$ (los sitios en que difieren). Por
tanto, $`M=-VH_NV`$ tiene las entradas fuera de la diagonal $`J_{ij}/2`$ (conjuntos distintos) o $`-J_{ij}/2`$ (mismo
conjunto), es decir, $`|J_{ij}|/2\ge0`$ en ambos casos, positivas si y solo si $`i\sim j`$ en $`G`$, y $`0`$ entre
configuraciones que difieren de otro modo. El grafo de $`M`$ es conexo. En efecto, escribamos $`(i\,j)`$ para el
intercambio de los espines de los sitios $`i,j`$; para $`a,b,c`$ distintos, $`(a\,c)=(a\,b)(b\,c)(a\,b)`$, así que, por
inducción a lo largo de un camino en $`G`$ de $`i`$ a $`j`$, todo $`(i\,j)`$ es una composición de intercambios a lo largo
de aristas de $`G`$, y cada uno de esos intercambios o bien fija una configuración (espines iguales) o bien recorre
una arista del grafo de $`M`$. Si $`b\ne b'`$ tienen $`N`$ espines hacia arriba, hay sitios $`i`$, hacia arriba en $`b`$ y
hacia abajo en $`b'`$, y $`j`$, hacia abajo en $`b`$ y hacia arriba en $`b'`$; $`(i\,j)b`$ difiere de $`b'`$ en dos sitios menos
que $`b`$, así que un número finito de tales pasos lleva de $`b`$ a $`b'`$. Por el lema 9.6, el mayor valor propio de
$`M`$ —el opuesto del menor de $`H_N`$— es simple, con un vector propio $`\varphi>0`$; tómese $`\psi=V\varphi`$. $`\square`$

<a id="r-9-8"></a>

**Corolario 9.8 (momento del estado fundamental).** En una cadena periódica con $`L`$ par, $`J>0`$ y un campo
uniforme, el estado fundamental del sector $`N`$ tiene momento $`k=0`$ si $`N`$ es par y $`k=\pi`$ si $`N`$ es impar.

*Demostración.* $`T`$ conmuta con $`H_N`$ y el estado fundamental $`\psi=V\varphi`$ es simple, así que $`T\psi=e^{ik}\psi`$.
Ahora bien, $`(T\psi)(b)=\psi(T^{-1}b)`$, y $`T^{-1}b`$ tiene en el sitio $`i`$ el espín $`b_{i+1}`$; como $`L`$ es par, $`i+1`$
es par para $`i`$ impar, luego $`N_B(T^{-1}b)=N-N_B(b)`$ y
$`(VT\psi)(b)=(-1)^{N_B(b)}(-1)^{N-N_B(b)}\varphi(T^{-1}b)=(-1)^N\varphi(T^{-1}b)`$. Además, $`VT\psi=e^{ik}\varphi`$. Tanto
$`\varphi`$ como $`\varphi\circ T^{-1}`$ son positivos, así que $`e^{ik}=(-1)^N`$. $`\square`$

*Cálculo:* el teorema 9.7 para cadenas abiertas y periódicas con campos aleatorios, y el corolario 9.8 para $`L=12`$ — [tutorial, sección 9.2](../../tutorial/es/09_spin_chains.md#s-9-2).

> **Significado físico** (interpretación física, no un teorema de este curso). En un antiferromagneto, la red
> fija los signos del estado fundamental: mover un espín a través de un enlace cambia el signo. Los físicos lo
> usan para ejecutar Monte Carlo cuántico sin problema del signo en antiferromagnetos bipartitos; en las cadenas
> frustradas, en las que un acoplamiento antiferromagnético une dos sitios del mismo conjunto, la hipótesis del
> teorema 9.7 no se cumple.

<a id="s-9-3"></a>

## 9.3 $`SU(2)`$ en $`\Delta=1`$

Para $`\Delta=1`$ y $`h=0`$, $`H=J\sum\mathbf s_i\cdot\mathbf s_j`$ conmuta con $`S^x,S^y,S^z`$ (proposición 4.15), así que
$`H\cong\bigoplus_S1_{2S+1}\otimes H_S`$ con bloques $`H_S\in M_{m_S}`$, $`m_S=\binom L{L/2-S}-\binom L{L/2-S-1}`$
(proposición 4.13), y el módulo `su2` construye $`H_S`$ sobre los vectores de peso máximo (proposición 4.15(b)).
Lo mismo vale para todo $`H=\sum_{i<j}J_{ij}\,\mathbf s_i\cdot\mathbf s_j`$ con $`J_{ij}`$ reales: por la proposición 4.15(a),
cada $`\mathbf s_i\cdot\mathbf s_j`$ está en $`\{S_a\}'`$. Sea $`E_0(S)`$ el menor valor propio de $`H_S`$, para los espines $`S`$ que
aparecen.

*En la biblioteca:* los bloques $`H_S`$ los construye `SU2Chain` — [tutorial, sección 9.3](../../tutorial/es/09_spin_chains.md#s-9-3).

<a id="r-9-9"></a>

**Teorema 9.9 (Lieb–Mattis).** Sea $`H=\sum_{i<j}J_{ij}\,\mathbf s_i\cdot\mathbf s_j`$ con $`A,B,J_{ij}`$ como en el teorema 9.7
(así que $`K_{ij}=J_{ij}`$, $`h_i=0`$), $`n_A=|A|\ge n_B=|B|`$ y $`S_*=\frac{n_A-n_B}2`$.
(a) Para $`S_*\le S\le\frac L2`$ (con $`S-\frac L2`$ entero), el estado fundamental del sector $`S^z=S`$ tiene espín total
$`S`$; para $`0\le M\le S_*`$ (con $`M-\frac L2`$ entero), el estado fundamental del sector $`S^z=M`$ tiene espín total $`S_*`$.
(b) Por tanto, $`E_0(S)`$ es un valor propio simple de $`H_S`$ para $`S\ge S_*`$, $`E_0(S_*)<E_0(S_*+1)<\dots<E_0(\frac L2)`$,
y $`E_0(S)>E_0(S_*)`$ para $`S<S_*`$: el espacio fundamental de $`H`$ es un único multiplete, de espín $`S_*`$.
Para la cadena de Heisenberg con $`J>0`$, abierta o periódica con $`L`$ par, $`A`$ y $`B`$ son los sitios pares e impares,
$`S_*=0`$ para $`L`$ par y $`S_*=\frac12`$ para $`L`$ impar (cadena abierta), y $`E_0(S_*)<E_0(S_*+1)<\dots<E_0(\frac L2)`$.

*Demostración.* (a) Sean $`M\ge0`$ con $`M-\frac L2\in\mathbb Z`$, $`s=\max(M,S_*)`$, y $`\psi`$ el estado fundamental del
sector $`S^z=M`$, simple por el teorema 9.7: $`\psi=V\varphi`$ con $`\varphi>0`$. $`\mathbf S^2`$ conmuta con $`H`$ y con $`S^z`$, así
que preserva el espacio fundamental unidimensional del sector: $`\psi`$ tiene un espín total $`S'`$. Construimos un
vector $`\chi`$ del sector con $`\mathbf S^2\chi=s(s+1)\chi`$ y $`\langle\psi,\chi\rangle\ne0`$; como los vectores propios de
$`\mathbf S^2`$ para valores propios distintos son ortogonales, esto da $`S'=s`$. Pongamos $`p=\frac L2-s`$; entonces
$`0\le p\le n_B`$, porque $`s\ge S_*`$ significa $`p\le\frac{n_A+n_B}2-\frac{n_A-n_B}2`$. Elijamos $`a_1,\dots,a_p\in A`$ distintos y
$`b_1,\dots,b_p\in B`$ distintos, sea $`R`$ el conjunto de los $`2s`$ sitios restantes, y sea
```math
\chi_0=\bigotimes_{l=1}^p\tfrac1{\sqrt2}\big(|{\uparrow}_{a_l}{\downarrow}_{b_l}\rangle-|{\downarrow}_{a_l}{\uparrow}_{b_l}\rangle\big)\otimes|{\uparrow\cdots\uparrow}\rangle_R ,
\qquad \chi=(S^-)^{s-M}\chi_0 .
```
Cada singlete es anulado por $`s^\pm_{a_l}+s^\pm_{b_l}`$ y los espines hacia arriba de $`R`$ por $`s^+_i`$, así que
$`S^+\chi_0=0`$; con $`S^z\chi_0=s\chi_0`$ esto da $`\mathbf S^2\chi_0=(S^-S^++(S^z)^2+S^z)\chi_0=s(s+1)\chi_0`$ (sección 4.5), y
$`\mathbf S^2\chi=s(s+1)\chi`$ porque $`\mathbf S^2`$ conmuta con $`S^-`$. Escribiendo $`S^-=\sum_l(s^-_{a_l}+s^-_{b_l})+S^-_R`$, los
sumandos conmutan y los primeros anulan los singletes, así que $`\chi`$ es el producto de los singletes con
$`(S^-_R)^{s-M}|{\uparrow\cdots\uparrow}\rangle_R=(s-M)!\sum_D|D\rangle`$, la suma sobre las configuraciones de $`R`$ con
exactamente $`s-M`$ espines hacia abajo (existen, pues $`s-M\le2s=|R|`$); así, $`\chi\ne0`$, $`S^z\chi=M\chi`$, y los
coeficientes de $`\chi`$ son reales. Su signo: un singlete contribuye $`+`$ cuando su espín hacia arriba está en $`A`$ y
$`-`$ cuando está en $`B`$, y la parte de $`R`$ es positiva. Si $`M=s`$, los sitios de $`R`$ están todos hacia arriba; si
$`M<s`$, entonces $`s=S_*`$, $`p=n_B`$ y $`R\subseteq A`$. En ambos casos, $`N_B(b)`$ es el número de singletes con el espín
hacia arriba en $`B`$ más una constante, $`|R\cap B|`$ (si $`M=s`$) o $`0`$, así que $`(-1)^{N_B(b)}\chi(b)`$ tiene el mismo
signo para todos los $`b`$ con $`\chi(b)\ne0`$. Por tanto,
$`\langle\psi,\chi\rangle=\langle\varphi,V\chi\rangle=\pm\sum_b\varphi(b)|\chi(b)|\ne0`$.

(b) Por la proposición 4.15(a), el corolario 4.12 y el teorema 4.11, un multiplete de espín $`S'`$ corta al sector
$`S^z=M`$ en un vector si $`S'\ge|M|`$ y en ninguno en otro caso, así que el espectro de $`H`$ sobre el sector $`S^z=M\ge0`$
es la unión de los espectros de los $`H_{S'}`$, $`S'\ge M`$, con multiplicidades. Por (a), su menor valor propio es
simple y pertenece al espín $`s=\max(M,S_*)`$: así, $`E_0(s)`$ es simple en $`H_s`$ y $`E_0(s)<E_0(S')`$ para todo $`S'\ge M`$,
$`S'\ne s`$. Para $`M=S\ge S_*`$ esto da la simplicidad y $`E_0(S)<E_0(S')`$ para $`S'>S`$; para $`M=S<S_*`$ da
$`E_0(S_*)<E_0(S)`$. Así, $`E_0(S_*)`$ es el menor valor propio de $`H`$, alcanzado solo en $`H_{S_*}`$, donde es simple.
Para la cadena, los sitios pares superan en número a los impares en $`0`$ ($`L`$ par) o en $`1`$ ($`L`$ impar). $`\square`$

El teorema 9.7 se remonta a W. Marshall (*Proc. Roy. Soc. A* 232 (1955) 48) y el teorema 9.9 a E. Lieb y
D. Mattis (*J. Math. Phys.* 3 (1962) 749). La **brecha singlete–triplete** $`E_0(1)-E_0(0)`$ del anillo de
Heisenberg par es positiva por el teorema 9.9.

*Cálculo:* el teorema 9.9 frente a los sectores $`U(1)`$ para $`L=12`$, la brecha singlete–triplete, y los teoremas 9.7 y 9.9 en un grafo bipartito que no es una cadena — [tutorial, sección 9.3](../../tutorial/es/09_spin_chains.md#s-9-3).

**Correlaciones.** En todo anillo de Heisenberg par, las correlaciones del estado fundamental
$`\langle\mathbf s_0\cdot\mathbf s_r\rangle`$ alternan de signo; los teoremas 9.7 y 9.9 lo imponen. Sea $`\psi`$ el estado fundamental
real normalizado, un singlete (teorema 9.9) y, por tanto, el estado fundamental del sector $`S^z=0`$. De
$`\mathbf S^2\psi=0`$ se sigue $`\|S^\pm\psi\|^2=\langle\psi,S^\mp S^\pm\psi\rangle=\langle\psi,(\mathbf S^2-(S^z)^2\mp S^z)\psi\rangle=0`$, así que
$`S^\pm\psi=0`$ y $`\langle\psi,[S^+,[S^-,X]]\psi\rangle=0`$ para todo $`X`$ (cada uno de los cuatro términos tiene $`S^\pm`$ junto a
$`\psi`$). Para $`X=s^z_0s^z_r`$, $`r\ne0`$, las relaciones $`[s^\pm,s^z]=\mp s^\pm`$, $`[s^+,s^-]=2s^z`$ dan
$`[S^-,X]=s^-_0s^z_r+s^z_0s^-_r`$ y $`[S^+,[S^-,X]]=4s^z_0s^z_r-s^+_0s^-_r-s^-_0s^+_r`$. Por tanto,
$`\langle s^z_0s^z_r\rangle=\frac14\langle s^+_0s^-_r+s^-_0s^+_r\rangle`$ y
```math
\langle\mathbf s_0\cdot\mathbf s_r\rangle=\langle s^z_0s^z_r\rangle+\tfrac12\langle s^+_0s^-_r+s^-_0s^+_r\rangle
=\tfrac34\sum_{b_0\ne b_r}\psi(b)\psi(b') ,
```
donde $`b'`$ se obtiene de $`b`$ intercambiando los espines de los sitios $`0`$ y $`r`$. Por el teorema 9.7,
$`\psi(b)\psi(b')`$ tiene el signo de $`(-1)^{N_B(b)+N_B(b')}`$, y $`N_B(b')-N_B(b)`$ vale $`\pm1`$ para $`r`$ impar y $`0`$ para
$`r`$ par; la suma no es vacía, así que $`\langle\mathbf s_0\cdot\mathbf s_r\rangle`$ es no nulo y tiene el signo $`(-1)^r`$ para
$`0<r<L`$.

*Cálculo:* las correlaciones $`\langle\mathbf s_0\cdot\mathbf s_r\rangle`$ para $`L=12`$, con `SU2Chain.spin_correlation`, que usa el teorema de Wigner–Eckart (teorema B.16) — [tutorial, sección 9.3](../../tutorial/es/09_spin_chains.md#s-9-3).

<a id="s-9-4"></a>

## 9.4 Comprobaciones exactas: fermiones libres y el ansatz de Bethe

**Jordan–Wigner.** En una cadena abierta, pongamos $`c_j=\big(\prod_{l<j}(-1)^{n_l}\big)s^-_j`$.

<a id="r-9-10"></a>

**Teorema 9.10 (la cadena XX es libre).** (a) $`\{c_i,c_j^*\}=\delta_{ij}`$, $`\{c_i,c_j\}=0`$ y $`c_j^*c_j=n_j`$.
(b) La cadena XX abierta es $`H=\frac J2\sum_{j<L-1}(s^+_js^-_{j+1}+s^-_js^+_{j+1})=\sum_{i,j}h_{ij}c_i^*c_j`$ con
$`h=\frac J2A_L`$. (c) Sean $`u_k`$ vectores propios reales ortonormales de $`h`$ con valores propios
$`\varepsilon_k=J\cos\frac{\pi k}{L+1}`$ (lema 9.3) y $`b_k=\sum_ju_k(j)c_j`$. Para $`K=\{k_1<\dots<k_N\}`$, los vectores
$`b_K^*|0\rangle=b_{k_1}^*\cdots b_{k_N}^*|0\rangle`$ ($`|0\rangle`$ con todos los espines hacia abajo) forman una base
ortonormal de vectores propios del sector $`N`$, con energías $`\sum_{k\in K}\varepsilon_k`$.

*Demostración.* (a) Los operadores en sitios distintos conmutan, y en un mismo sitio $`(-1)^n`$ anticonmuta con
$`s^\pm`$. Para $`i<j`$, $`c_i`$ y $`c_j^*`$ (o $`c_j`$) se encuentran en el sitio $`i`$ como $`s^-_i`$ y $`(-1)^{n_i}`$ —una
anticonmutación— y por lo demás conmutan, así que anticonmutan. Para $`i=j`$, las cadenas de signos tienen cuadrado
$`1`$ y $`c_jc_j^*+c_j^*c_j=s^-_js^+_j+s^+_js^-_j=1`$, $`c_j^*c_j=s^+_js^-_j=n_j`$, $`c_j^2=(s^-_j)^2=0`$; el caso $`i>j`$ se
sigue intercambiando los papeles de $`i,j`$ (y tomando adjuntos). (b)
$`c_j^*c_{j+1}=s^+_j(-1)^{n_j}s^-_{j+1}=s^+_js^-_{j+1}`$, porque $`s^+_j\neq0`$ solo sobre $`n_j=0`$; se toman adjuntos.
(c) $`\{b_k,b_l^*\}=\sum_ju_k(j)u_l(j)=\delta_{kl}`$ y $`\{b_k,b_l\}=0`$, y $`H=\sum_k\varepsilon_kb_k^*b_k`$ por la
descomposición espectral de $`h`$. De las CAR, $`[b_l^*b_l,b_k^*]=\delta_{lk}b_k^*`$, así que $`[H,b_k^*]=\varepsilon_kb_k^*`$;
con $`H|0\rangle=0`$ esto da $`Hb_K^*|0\rangle=(\sum_K\varepsilon_k)b_K^*|0\rangle`$. Como $`b_k|0\rangle=0`$ (cada $`c_j`$ baja un
espín), llevar los aniquiladores a la derecha con las CAR da $`\langle b_{K'}^*0,b_K^*0\rangle=\delta_{KK'}`$ por inducción
sobre $`N`$. Cada $`c_j^*`$ aumenta $`N`$ en uno, así que $`b_K^*|0\rangle`$ está en el sector $`|K|`$: $`\binom LN`$ vectores
ortonormales en un espacio de esa dimensión. $`\square`$

*Cálculo:* los $`1024`$ niveles de la cadena XX abierta con $`L=10`$, frente al teorema 9.10 — [tutorial, sección 9.4](../../tutorial/es/09_spin_chains.md#s-9-4).

**Magnones.** En el anillo con $`h=0`$, el estado totalmente polarizado $`|0\rangle`$ tiene energía $`E_F=J\Delta L/4`$. Un
espín hacia arriba en $`x`$ rompe dos enlaces: $`H|x\rangle=(E_F-J\Delta)|x\rangle+\frac J2(|x-1\rangle+|x+1\rangle)`$ ($`L\ge3`$).

<a id="r-9-11"></a>

**Proposición 9.11 (un magnón).** En el sector $`N=1`$ del anillo, $`|k\rangle=\sum_xe^{-ikx}|x\rangle`$, $`k=2\pi m/L`$,
cumple $`T|k\rangle=e^{ik}|k\rangle`$ y $`H|k\rangle=\big(E_F+J(\cos k-\Delta)\big)|k\rangle`$.

*Demostración.* $`\sum_xe^{-ikx}(|x-1\rangle+|x+1\rangle)=\sum_y(e^{-ik(y+1)}+e^{-ik(y-1)})|y\rangle=2\cos k\,|k\rangle`$, usando
$`e^{ikL}=1`$; la misma reindexación da $`T|k\rangle=e^{ik}|k\rangle`$. $`\square`$

Dos magnones se dispersan. Escribamos $`|x_1,x_2\rangle`$, $`x_1<x_2`$, para los espines hacia arriba en $`x_1,x_2`$; un
par de vecinos rompe dos enlaces, un par separado cuatro.

<a id="r-9-12"></a>

**Proposición 9.12 (dos magnones; ecuaciones de Bethe).** Sean $`k_1,k_2\in\mathbb C`$, $`K=k_1+k_2`$, y
$`(A_{12},A_{21})\neq0`$ con (i) $`A_{12}(1+e^{iK}-2\Delta e^{ik_2})+A_{21}(1+e^{iK}-2\Delta e^{ik_1})=0`$ y
(ii) $`A_{21}=e^{ik_2L}A_{12}`$, $`A_{12}=e^{ik_1L}A_{21}`$. Sean
$`\psi(x_1,x_2)=A_{12}e^{i(k_1x_1+k_2x_2)}+A_{21}e^{i(k_2x_1+k_1x_2)}`$ sobre $`\mathbb Z^2`$ y
$`\Psi=\sum_{0\le x_1<x_2<L}\psi(x_1,x_2)|x_1,x_2\rangle`$. Si $`\Psi\neq0`$, entonces
$`H\Psi=\big(E_F+J(\cos k_1+\cos k_2-2\Delta)\big)\Psi`$. Si los denominadores no se anulan, (i) y (ii) dicen
```math
e^{ik_1L}=-\frac{1+e^{iK}-2\Delta e^{ik_1}}{1+e^{iK}-2\Delta e^{ik_2}},\qquad
e^{ik_2L}=-\frac{1+e^{iK}-2\Delta e^{ik_2}}{1+e^{iK}-2\Delta e^{ik_1}} .
```

*Demostración.* Sea $`E=E_F+J(\cos k_1+\cos k_2-2\Delta)`$. (1) Para todo $`(x_1,x_2)\in\mathbb Z^2`$, cada onda plana da
valores $`\psi(x_1\pm1,x_2)`$ y $`\psi(x_1,x_2\pm1)`$ que suman $`2(\cos k_1+\cos k_2)\psi(x_1,x_2)`$. (2) Sustituyendo,
$`\psi(x,x)+\psi(x+1,x+1)-2\Delta\psi(x,x+1)`$ es $`e^{iKx}`$ por el miembro izquierdo de (i), así que se anula para todo
$`x`$. (3) Por (ii), $`\psi(x_2,x_1+L)=A_{12}e^{ik_2L}e^{i(k_2x_1+k_1x_2)}+A_{21}e^{ik_1L}e^{i(k_1x_1+k_2x_2)}=\psi(x_1,x_2)`$, y
aplicando esto dos veces, $`\psi(x_1+L,x_2+L)=\psi(x_1,x_2)`$. Los representantes de una configuración $`\{y_1<y_2\}`$
en $`D=\{x_1<x_2<x_1+L\}`$ son $`(y_1+nL,y_2+nL)`$ y $`(y_2+nL,y_1+(n+1)L)`$, así que $`\psi`$ es una función de las
configuraciones. Sea $`(x_1,x_2)\in D`$. Si $`2\le x_2-x_1\le L-2`$, los cuatro saltos permanecen en $`D`$ y el término
diagonal es $`E_F-2J\Delta`$, así que $`(H\Psi)(x_1,x_2)=(E_F-2J\Delta)\psi+\frac J2\cdot2(\cos k_1+\cos k_2)\psi=E\psi`$ por (1).
Si $`x_2=x_1+1`$, los saltos $`(x_1+1,x_2)`$, $`(x_1,x_2-1)`$ están bloqueados y el término diagonal es $`E_F-J\Delta`$; por
(1) y (2),
$`(H\Psi)(x_1,x_2)=(E_F-J\Delta)\psi+\frac J2\big(2(\cos k_1+\cos k_2)\psi-\psi(x_1,x_1)-\psi(x_1+1,x_1+1)\big)=E\psi`$.
Si $`x_2=x_1+L-1`$, se usa el representante $`(x_2,x_1+L)`$, que es del tipo anterior ($`L\ge3`$). $`\square`$

La proposición 9.12 solo dice que las soluciones dan vectores propios; no dice que todos los vectores propios
surjan así. Escribamos $`k_{1,2}=K/2\pm q`$ y $`c=\cos\frac K2`$. Entonces
$`1+e^{iK}-2\Delta e^{ik_{1,2}}=2e^{iK/2}(c-\Delta e^{\pm iq})`$, así que $`A_{12}=c-\Delta e^{iq}`$,
$`A_{21}=-(c-\Delta e^{-iq})`$ resuelven (i). Por (ii) y $`(A_{12},A_{21})\ne0`$, ambos coeficientes son no nulos, así
que multiplicar las dos ecuaciones (ii) da $`e^{iKL}=1`$, $`K=2\pi n/L`$, y
entonces (ii) equivale a $`F(q)=(-1)^ne^{iqL}(c-\Delta e^{-iq})+c-\Delta e^{iq}=0`$, una ecuación en una variable.
Buscamos sus raíces $`q\in(0,\pi)`$, donde $`e^{-iqL/2}F(q)`$ es real o imaginario porque $`c-\Delta e^{\pm iq}`$ son
complejos conjugados (estados de dispersión), y sobre las rectas $`q=i\kappa`$, $`\pi+i\kappa`$, $`\kappa>0`$, donde $`F`$ es
real (estados ligados).

Para $`L`$ par y $`\Delta\ne0`$, los vectores de Bethe no son completos. Si $`\Delta\ne0`$,
$`\Phi=\sum_x(-1)^x|x,x+1\rangle`$ cumple $`T\Phi=-\Phi`$ y $`H\Phi=(E_F-J\Delta)\Phi`$ (los saltos a $`|x-1,x+1\rangle`$ y
$`|x,x+2\rangle`$ se cancelan por pares). Un vector de Bethe cumple $`T\Psi=e^{-iK}\Psi`$, ya que
$`\psi(x_1-1,x_2-1)=e^{-iK}\psi(x_1,x_2)`$. Si $`e^{-iK}\ne-1`$, $`\Psi\perp\Phi`$ porque son vectores propios del unitario $`T`$
para valores propios distintos; si $`e^{iK}=-1`$, entonces $`\cos k_1+\cos k_2=2\cos\frac K2\cos\frac{k_1-k_2}2=0`$, la
energía de $`\Psi`$ es $`E_F-2J\Delta\ne E_F-J\Delta`$, y $`\Psi\perp\Phi`$ como vectores propios de $`H`$ para valores propios
distintos. Así, $`\Phi`$ es ortogonal a todos los vectores de Bethe de la proposición 9.12, que por tanto nunca
generan el sector $`N=2`$ de un anillo par. El cálculo
de dos magnones es el caso $`N=2`$ del ansatz de H. Bethe (*Z. Phys.* 71 (1931) 205); no usamos ningún otro caso.

*Cálculo:* las raíces de $`F`$ para $`L=9`$, $`\Delta=0.6`$, que dan los $`36`$ estados del sector $`N=2`$, y el vector $`\Phi`$ para $`L=8`$ — [tutorial, sección 9.4](../../tutorial/es/09_spin_chains.md#s-9-4).

> **Significado físico** (interpretación física, no un teorema de este curso). Un magnón es una onda de espín:
> un espín invertido deslocalizado sobre el anillo con momento $`k`$. Dos magnones se dispersan con un desfase (con
> $`q`$ real) o forman un estado ligado (con $`q`$ complejo).

<a id="s-9-5"></a>

## 9.5 Lanczos y métodos de Krylov

El sector $`N=L/2`$ tiene dimensión $`12870`$ para $`L=16`$ y $`184756`$ para $`L=20`$, pero $`H`$ es dispersa ($`O(L)`$ entradas
no nulas por fila, `xxz_sparse`), y los métodos de Krylov solo necesitan productos $`Hv`$. Sea $`H`$ hermítica de
tamaño $`n`$ y $`v\ne0`$. El **espacio de Krylov** es $`\mathcal K_m=\mathrm{span}\{v,Hv,\dots,H^{m-1}v\}`$.

<a id="r-9-13"></a>

**Proposición 9.13 (Lanczos).** Supongamos $`\dim\mathcal K_m=m`$. Sea $`q_1=v/\|v\|`$ y, para $`j<m`$, sea $`q_{j+1}`$ la
componente normalizada de $`Hq_j`$ ortogonal a $`\mathcal K_j`$. Entonces $`q_1,\dots,q_m`$ es una base ortonormal de
$`\mathcal K_m`$; con $`Q_m=(q_1\cdots q_m)`$, $`T_m=Q_m^*HQ_m`$ es real, simétrica y tridiagonal, con diagonal
$`\alpha_j=\langle q_j,Hq_j\rangle`$ y subdiagonal $`\beta_j=\langle q_{j+1},Hq_j\rangle>0`$; y $`HQ_m=Q_mT_m+r_me_m^T`$, donde
$`r_m\perp\mathcal K_m`$ se anula si y solo si $`\mathcal K_{m+1}=\mathcal K_m`$ (en otro caso $`r_m=\beta_mq_{m+1}`$,
$`q_{m+1}=r_m/\|r_m\|`$).

*Demostración.* Supongamos $`\mathrm{span}\{q_1,\dots,q_j\}=\mathcal K_j`$ (cierto para $`j=1`$). Entonces
$`H\mathcal K_j\subseteq\mathcal K_j+\mathbb CHq_j`$, pues $`Hq_i\in\mathcal K_{i+1}`$; así que
$`\mathcal K_{j+1}=\mathcal K_j+\mathbb CH^jv\subseteq\mathcal K_j+\mathbb CHq_j\subseteq\mathcal K_{j+1}`$, y
$`\mathcal K_{j+1}=\mathcal K_j+\mathbb CHq_j`$. Para $`j<m`$ la dimensión crece, de modo que la componente de $`Hq_j`$
ortogonal a $`\mathcal K_j`$ es no nula y la inducción continúa; para $`j=m`$ esa componente, $`r_m`$, se anula si y solo
si $`\mathcal K_{m+1}=\mathcal K_m`$. Para $`i>j+1`$, $`q_i\perp\mathcal K_{j+1}\ni Hq_j`$, así que $`(T_m)_{ij}=0`$, y
$`(T_m)_{ji}=0`$ porque $`T_m`$ es hermítica. $`\alpha_j`$ es real, y $`\beta_j=\langle q_{j+1},Hq_j\rangle`$ es la norma de la
componente, $`>0`$. Así, $`Hq_j=\beta_{j-1}q_{j-1}+\alpha_jq_j+\beta_jq_{j+1}`$ para $`j<m`$ y
$`Hq_m=\beta_{m-1}q_{m-1}+\alpha_mq_m+r_m`$. $`\square`$

Los valores propios $`\theta^{(m)}_1\le\dots\le\theta^{(m)}_m`$ de $`T_m`$ son los **valores de Ritz**.

<a id="r-9-14"></a>

**Teorema 9.14 (convergencia monótona).** Sean $`\lambda_1\le\dots\le\lambda_n`$ los valores propios de $`H`$ y
$`\dim\mathcal K_m=m`$. Entonces $`\theta^{(m)}_j\ge\lambda_j`$, y $`\theta^{(m)}_j\ge\theta^{(m+1)}_j`$ si
$`\dim\mathcal K_{m+1}=m+1`$ ($`j\le m`$). Si $`v`$ tiene una componente no nula en el espacio propio de $`\lambda_1`$ y
$`\mathcal K_{m+1}=\mathcal K_m`$, entonces $`\theta^{(m)}_1=\lambda_1`$.

*Demostración.* Por el corolario A.9 (Courant–Fischer para compresiones, teorema A.8) aplicado a $`Q=Q_m`$,
$`\theta^{(m)}_j=\min_{U}\max_{x\in U,\|x\|=1}\langle x,Hx\rangle`$ sobre los $`U\subseteq\mathcal K_m=\mathrm{ran}\,Q_m`$ de
dimensión $`j`$. La misma fórmula sobre $`U\subseteq\mathcal K_{m+1}`$ da $`\theta^{(m+1)}_j`$ y sobre $`U\subseteq\mathbb C^n`$ da
$`\lambda_j`$ (teorema A.8); minimizar sobre más subespacios da un valor menor. Si $`\mathcal K_{m+1}=\mathcal K_m`$, entonces $`\mathcal K_m`$ es $`H`$-invariante y contiene
$`p(H)v`$ para todo polinomio $`p`$, en particular $`P_1v\ne0`$ ($`P_1`$ la proyección espectral de $`\lambda_1`$, un
polinomio en $`H`$ por interpolación de Lagrange, como en el teorema 1.5); así que
$`\theta^{(m)}_1\le\langle x,Hx\rangle=\lambda_1`$ en $`x=P_1v/\|P_1v\|\in\mathcal K_m`$. $`\square`$

Un residuo pequeño certifica un valor propio cercano: para $`H`$ hermítica, un número $`E`$ y $`\|\psi\|=1`$,
$`\|H\psi-E\psi\|^2=\sum_i|\lambda_i-E|^2|\langle u_i,\psi\rangle|^2\ge\min_i|\lambda_i-E|^2`$ ($`u_i`$ una base ortonormal de vectores
propios). No certifica que ese valor propio sea el menor. Para $`L=4`$, la energía del estado fundamental del anillo
de Heisenberg es exacta: con $`\mathbf S_A=\mathbf s_0+\mathbf s_2`$, $`\mathbf S_B=\mathbf s_1+\mathbf s_3`$, el
anillo es $`H=J\,\mathbf S_A\cdot\mathbf S_B=\frac J2(\mathbf S^2-\mathbf S_A^2-\mathbf S_B^2)\ge\frac J2(0-2-2)=-2J`$, ya que
$`\mathbf S^2\ge0`$ y $`\mathbf S_A^2,\mathbf S_B^2`$ tienen los valores propios $`0,2`$ (proposición 4.13 con $`L=2`$); la igualdad se
da en un vector de espín $`0`$ del producto de los dos tripletes, que existe por el capítulo 4, ejercicio 2
($`1\otimes1\cong2\oplus1\oplus0`$).

*Cálculo:* el teorema 9.14 para la cadena XX abierta, con `krylov.lanczos` sobre el sector disperso de `SpinChain.xxz_sparse`, y estados fundamentales del anillo de Heisenberg hasta $`L=16`$ con `krylov.ground_state` y su residuo — [tutorial, sección 9.5](../../tutorial/es/09_spin_chains.md#s-9-5).

**Evolución temporal.** Para $`e^{-i\tau H}v`$ se proyecta sobre $`\mathcal K_m`$: $`w_\tau=\|v\|\,Q_me^{-i\tau T_m}e_1`$.

<a id="r-9-15"></a>

**Proposición 9.15.** Sea $`\dim\mathcal K_m=m`$. Para todo $`\tau`$, $`\|w_\tau\|=\|v\|`$ y
$`\langle w_\tau,Hw_\tau\rangle=\langle v,Hv\rangle`$. Si $`\mathcal K_m`$ es $`H`$-invariante, $`w_\tau=e^{-i\tau H}v`$.

*Demostración.* $`Q_m`$ es una isometría y $`e^{-i\tau T_m}`$ es unitaria ($`T_m`$ es real simétrica), lo que da la norma.
Como $`Q_m^*HQ_m=T_m`$ conmuta con $`e^{-i\tau T_m}`$,
$`\langle w_\tau,Hw_\tau\rangle=\|v\|^2\langle e^{-i\tau T_m}e_1,T_me^{-i\tau T_m}e_1\rangle=\|v\|^2\langle e_1,T_me_1\rangle=\langle v,Hv\rangle`$. Si
$`\mathcal K_m`$ es invariante, $`HQ_m=Q_mT_m`$ (proposición 9.13, $`r_m=0`$), así que $`H^jQ_m=Q_mT_m^j`$ y, sumando la
serie exponencial (proposición 8.2(1)), $`e^{-i\tau H}Q_m=Q_me^{-i\tau T_m}`$; se aplica a $`e_1`$, con $`Q_me_1=v/\|v\|`$. $`\square`$

*Cálculo:* la proposición 9.15 para `krylov.evolve`, que avanza por pasos en espacios de Krylov, en un cambio súbito (quench) desde el estado de Néel, frente al propagador denso — [tutorial, sección 9.5](../../tutorial/es/09_spin_chains.md#s-9-5).

**Cotas de error.** Salvo que $`\mathcal K_m`$ sea invariante, la proposición 9.15 no dice nada sobre la precisión de
$`w_\tau`$. El error depende de lo bien que los polinomios de grado $`<m`$ aproximan $`e^{-i\tau x}`$ en el espectro, y los
polinomios de Chebyshev proporcionan la aproximación. Los **polinomios de Chebyshev** $`\mathcal T_0=1`$, $`\mathcal T_1(y)=y`$,
$`\mathcal T_{k+1}(y)=2y\mathcal T_k(y)-\mathcal T_{k-1}(y)`$ tienen grado $`k`$ y cumplen $`\mathcal T_k(\cos\vartheta)=\cos k\vartheta`$ (por inducción, a
partir de $`\cos(k+1)\vartheta+\cos(k-1)\vartheta=2\cos\vartheta\cos k\vartheta`$), luego $`|\mathcal T_k|\le1`$ en $`[-1,1]`$. Escribimos
```math
S_m(x)=\sum_{j\ge m}\frac{x^j}{j!}\le\frac{x^m}{m!}\Big(1-\frac x{m+1}\Big)^{-1}\qquad(0\le x<m+1),
```
donde la cota se sigue de $`m!/(m+i)!\le(m+1)^{-i}`$ y la serie geométrica.

<a id="r-9-17"></a>

**Lema 9.17 (desarrollo de Chebyshev de $`e^{-i\rho y}`$).** Para $`\rho`$ real e $`y\in[-1,1]`$,
```math
e^{-i\rho y}=a_0+2\sum_{k\ge1}a_k\mathcal T_k(y),\qquad a_k=(-i)^k\sum_{l\ge0}\frac{(-1)^l(\rho/2)^{k+2l}}{l!\,(k+l)!},\qquad|a_k|\le\frac{(|\rho|/2)^k}{k!}.
```
Por tanto, para todo $`m\ge1`$, el polinomio $`p=a_0+2\sum_{k=1}^{m-1}a_k\mathcal T_k`$ de grado $`\le m-1`$ cumple
$`|e^{-i\rho y}-p(y)|\le2S_m(|\rho|/2)`$ en $`[-1,1]`$.

(Las sumas que definen $`a_k`$ son las funciones de Bessel, $`a_k=(-i)^kJ_k(\rho)`$; la demostración no lo necesita.)

*Demostración.* *Desarrollo.* Para $`y=\cos\vartheta`$, $`\cos^j\vartheta=2^{-j}(e^{i\vartheta}+e^{-i\vartheta})^j=2^{-j}\sum_{l=0}^j\binom jle^{i(j-2l)\vartheta}`$, luego
```math
e^{-i\rho\cos\vartheta}=\sum_{j\ge0}\frac{(-i\rho)^j}{j!}\,2^{-j}\sum_{l=0}^j\binom jl\,e^{i(j-2l)\vartheta}.
```
La serie doble converge absolutamente (los módulos de sus términos suman $`e^{|\rho|}`$), así que puede sumarse en
cualquier orden, en particular según $`k=j-2l`$. Para $`k\ge0`$ los términos con $`j-2l=k`$ tienen $`j=k+2l`$, y como
$`\binom{k+2l}l/(k+2l)!=1/(l!\,(k+l)!)`$ suman $`a_ke^{ik\vartheta}`$; para $`k<0`$, $`\binom jl=\binom j{j-l}`$ muestra que el
coeficiente de $`e^{ik\vartheta}`$ es $`a_{|k|}`$. Por tanto
$`e^{-i\rho\cos\vartheta}=a_0+\sum_{k\ge1}a_k(e^{ik\vartheta}+e^{-ik\vartheta})=a_0+2\sum_{k\ge1}a_k\mathcal T_k(\cos\vartheta)`$.

*Cota.* Fijemos $`k\ge0`$. Sea $`\vartheta`$ uniforme en $`[0,2\pi]`$; sea $`\sigma=1`$ si $`k=0`$, y para $`k\ge1`$ sea $`\sigma`$
independiente de $`\vartheta`$ con densidad $`k(1-t)^{k-1}`$ en $`[0,1]`$. Pongamos $`Z=\sqrt\sigma\cos\vartheta\in[-1,1]`$.
Entonces $`\mathbb E\cos^{2l}\vartheta=4^{-l}\binom{2l}l`$ (en el desarrollo binomial anterior solo el término constante tiene
media no nula), y $`\mathbb E\sigma^l=k\int_0^1t^l(1-t)^{k-1}dt=l!\,k!/(k+l)!`$ (integrando por partes $`l`$ veces,
$`\int_0^1t^l(1-t)^{k-1}dt=\frac lk\int_0^1t^{l-1}(1-t)^kdt`$; para $`k=0`$ ambos lados valen $`1`$). Por independencia,
$`\mathbb EZ^{2l}=(2l)!\,k!/(4^l\,l!\,(k+l)!)`$. La serie $`\cos(\rho Z)=\sum_l(-1)^l(\rho Z)^{2l}/(2l)!`$ converge uniformemente,
pues $`|Z|\le1`$, así que su media puede tomarse término a término:
```math
\frac{(\rho/2)^k}{k!}\,\mathbb E\cos(\rho Z)=\sum_{l\ge0}\frac{(-1)^l\rho^{2l}}{(2l)!}\cdot\frac{(\rho/2)^k\,(2l)!}{4^l\,l!\,(k+l)!}=\sum_{l\ge0}\frac{(-1)^l(\rho/2)^{k+2l}}{l!\,(k+l)!}=i^ka_k.
```
Como $`|\cos|\le1`$, $`|a_k|\le(|\rho|/2)^k/k!`$.

*Truncación.* $`|e^{-i\rho y}-p(y)|=2\,|\sum_{k\ge m}a_k\mathcal T_k(y)|\le2\sum_{k\ge m}(|\rho|/2)^k/k!=2S_m(|\rho|/2)`$. $`\square`$

<a id="r-9-18"></a>

**Teorema 9.18 (error del propagador de Krylov).** Sean $`\dim\mathcal K_m=m\ge2`$ y $`\|v\|=1`$, $`Q_m,T_m,r_m`$ como en la
proposición 9.13, $`\beta_m=\|r_m\|`$, $`\theta_1\le\dots\le\theta_m`$ los valores de Ritz, $`c_T=(\theta_1+\theta_m)/2`$ y
$`w_T=(\theta_m-\theta_1)/2`$. Sean $`w_\tau=Q_me^{-i\tau T_m}e_1`$ y $`\varphi(t)=\langle e_m,e^{-itT_m}e_1\rangle`$, y sean
$`c\in\mathbb R`$, $`w\ge0`$ con $`\mathrm{spec}\,H\subseteq[c-w,c+w]`$. Para $`\tau`$ real:

(a) (*invariancia por desplazamiento*) para $`s`$ real, los vectores de Lanczos de $`H+s`$ y $`v`$ son los de $`H`$ y $`v`$,
$`T_m`$ pasa a ser $`T_m+s`$, $`r_m`$ no cambia, y $`e^{-i\tau H}v`$ y $`w_\tau`$ quedan ambos multiplicados por $`e^{-i\tau s}`$; el
error $`\|e^{-i\tau H}v-w_\tau\|`$ y las cotas (b)–(d) no cambian;

(b) (*a posteriori*) $`e^{-i\tau H}v-w_\tau=-i\int_0^\tau\varphi(t)\,e^{-i(\tau-t)H}r_m\,dt`$, luego
$`\|e^{-i\tau H}v-w_\tau\|\le\beta_m\,\big|\int_0^\tau|\varphi(t)|\,dt\big|`$;

(c) (*cota calculable*) $`w_T>0`$, $`|\varphi(t)|\le2S_{m-1}(|t|w_T/2)`$, y
```math
\|e^{-i\tau H}v-w_\tau\|\le\frac{4\beta_m}{w_T}\,S_m\Big(\frac{|\tau|\,w_T}2\Big);
```

(d) (*a priori*) $`\beta_m\le w`$ y $`w_T\le w`$, y con $`\rho=|\tau|\,w`$
```math
\|e^{-i\tau H}v-w_\tau\|\le\frac{4\beta_m}w\,S_m\Big(\frac\rho2\Big)\le4S_m\Big(\frac\rho2\Big)\le\frac4{m!}\Big(\frac\rho2\Big)^m\Big(1-\frac\rho{2m+2}\Big)^{-1}\qquad(\rho<2m+2).
```

*Demostración.* (a) $`(H+s)^jv`$ es $`H^jv`$ más una combinación de los $`H^iv`$, $`i<j`$, y recíprocamente, así que
$`\mathcal K_j(H+s,v)=\mathcal K_j(H,v)`$ para todo $`j`$. La componente de $`(H+s)q_j=Hq_j+sq_j`$ ortogonal a $`\mathcal K_j\ni q_j`$ es
la de $`Hq_j`$, de modo que la construcción de la proposición 9.13 da los mismos $`q_j`$; entonces $`Q_m^*(H+s)Q_m=T_m+s`$ y
$`(H+s)Q_m-Q_m(T_m+s)=HQ_m-Q_mT_m`$. Por el teorema 1.5(c), $`e^{-i\tau(X+s)}=e^{-i\tau s}e^{-i\tau X}`$ para $`X=H`$ y
$`X=T_m`$. Finalmente, $`\beta_m`$, $`|\varphi|`$ y $`w_T`$ no cambian, y $`\mathrm{spec}(H+s)\subseteq[c+s-w,c+s+w]`$.

(b) Por la proposición 8.2(1) para $`T_m`$, $`\frac d{dt}w_t=-iQ_mT_me^{-itT_m}e_1`$, y $`Q_mT_m=HQ_m-r_me_m^T`$ (proposición
9.13), luego $`\frac d{dt}w_t=-iHw_t+i\varphi(t)r_m`$. Como $`\frac d{dt}e^{itH}=iHe^{itH}`$ (proposición 8.2(1)), la regla del
producto da $`\frac d{dt}(e^{itH}w_t)=i\varphi(t)e^{itH}r_m`$. Integrando de $`0`$ a $`\tau`$ (lema A.4(e)), con $`w_0=Q_me_1=v`$,
$`e^{i\tau H}w_\tau-v=i\int_0^\tau\varphi(t)e^{itH}r_m\,dt`$. Se aplica la unitaria $`e^{-i\tau H}`$ (lema A.4(b)) y se toman
normas (lema A.4(c)).

(c) $`T_m`$ no es múltiplo de $`1`$ ($`\beta_1>0`$), así que sus valores propios no son todos iguales (teorema 1.3) y
$`w_T>0`$. Como $`T_m`$ es tridiagonal, $`T_m^je_1`$ se anula por debajo de su entrada $`j+1`$, luego $`\langle e_m,p(T_m)e_1\rangle=0`$
para todo polinomio $`p`$ de grado $`\le m-2`$. Por el lema 9.17, con $`m-1`$ en lugar de $`m`$ y $`\rho=tw_T`$, hay un
polinomio $`\tilde p`$ de grado $`\le m-2`$ con $`|e^{-itw_Ty}-\tilde p(y)|\le2S_{m-1}(|t|w_T/2)`$ en $`[-1,1]`$; entonces
$`p(x)=e^{-itc_T}\tilde p((x-c_T)/w_T)`$ tiene grado $`\le m-2`$ y $`|e^{-itx}-p(x)|\le2S_{m-1}(|t|w_T/2)`$ en cada valor de
Ritz. Por el teorema 1.5 y el corolario 1.4(1),
$`|\varphi(t)|=|\langle e_m,(e^{-itT_m}-p(T_m))e_1\rangle|\le\max_i|e^{-it\theta_i}-p(\theta_i)|`$. La serie $`S_{m-1}`$ tiene términos
no negativos y converge uniformemente en intervalos acotados, así que puede integrarse término a término:
$`\int_0^{|\tau|}2S_{m-1}(tw_T/2)\,dt=\frac4{w_T}S_m(|\tau|w_T/2)`$; con (b) esta es la cota.

(d) $`Q_m(T_m-c)e_m=Q_mQ_m^*(H-c)q_m`$ es la proyección ortogonal de $`(H-c)q_m`$ sobre $`\mathcal K_m`$, y
$`r_m=(H-c)q_m-Q_m(T_m-c)e_m`$ es la componente ortogonal, así que $`\beta_m\le\|(H-c)q_m\|\le\|H-c\|\le w`$ (corolario 1.4(1)).
Para un vector propio unitario $`y`$ de $`T_m`$, el valor de Ritz $`\langle y,T_my\rangle=\langle Q_my,HQ_my\rangle`$ está en
$`[c-w,c+w]`$, luego $`w_T\le w`$. Ahora $`\frac{4\beta_m}{w_T}S_m(x)=2\beta_m|\tau|\,S_m(x)/x`$ con $`x=|\tau|w_T/2\le\rho/2`$, y
$`S_m(x)/x=\sum_{j\ge m}x^{j-1}/j!`$ es no decreciente; esto da la primera desigualdad, $`\beta_m\le w`$ la segunda y la
cota de $`S_m`$ la tercera. $`\square`$

La cota (c) usa solo los datos de la ejecución de Lanczos, $`\beta_m`$ y los valores de Ritz. Los valores de Ritz están
dentro del espectro (teorema 9.14) y pueden subestimar mucho su semianchura $`w`$, para un vector inicial con
componentes pequeñas sobre los vectores propios extremos; pero ninguna estimación de $`w`$ entra en (c), porque la
representación (b) solo involucra a $`T_m`$. Por (d) y $`m!\ge(m/e)^m`$, el error es a lo sumo del orden de
$`(e\rho/2m)^m`$: es pequeño en cuanto $`\rho=|\tau|w`$ está bastante por debajo de $`2m/e\approx0.74\,m`$, así que los
pasos pueden ser largos. La cantidad $`\beta_m|\varphi(\tau)|`$, que versiones anteriores de la biblioteca sumaban como
estimación del error, es por la demostración de (b) la norma del defecto $`\frac d{d\tau}w_\tau+iHw_\tau`$: el integrando de
(b) en el extremo, una tasa y no un error. Para $`\tau`$ pequeño, $`\varphi(t)\approx(-it)^{m-1}\beta_1\cdots\beta_{m-1}/(m-1)!`$, así
que supera la cota (b) por un factor de aproximadamente $`m/|\tau|`$; para $`\tau`$ grande, $`\varphi`$ oscila y puede
anularse en el extremo sin que el error se anule.

<a id="r-9-19"></a>

**Corolario 9.19 (regla de paso).** Sean $`\tau=\tau_1+\dots+\tau_N`$ con $`\tau_k>0`$, $`\psi_0=v\ne0`$,
$`\psi_k=\|\psi_{k-1}\|\,Q_me^{-i\tau_kT_m}e_1`$ la aproximación de Krylov de $`e^{-i\tau_kH}\psi_{k-1}`$ en
$`\mathcal K_m(H,\psi_{k-1})`$, y $`\varepsilon_k`$ la cota (c) de este paso multiplicada por $`\|\psi_{k-1}\|`$ ($`\varepsilon_k=0`$ si
este espacio de Krylov tiene dimensión $`<m`$). Entonces $`\|e^{-i(\tau_1+\dots+\tau_k)H}v-\psi_k\|\le\varepsilon_1+\dots+\varepsilon_k`$ para
todo $`k`$. En particular, si $`\varepsilon_k\le\mathrm{tol}\cdot\|v\|\,\tau_k/\tau`$ en cada paso, el error tras $`k`$ pasos es a lo
sumo $`\mathrm{tol}\cdot\|v\|\,(\tau_1+\dots+\tau_k)/\tau`$. Para $`m\ge2`$ la cota de un paso de longitud $`h`$, dividida por
$`h`$, es no decreciente en $`h`$ y tiende a $`0`$ cuando $`h\to0`$, así que el paso admisible más largo es positivo.

*Demostración.* Sea $`d_k=e^{-i(\tau_1+\dots+\tau_k)H}v-\psi_k`$. Entonces $`d_k=e^{-i\tau_kH}d_{k-1}+(e^{-i\tau_kH}\psi_{k-1}-\psi_k)`$,
y $`e^{-i\tau_kH}`$ es unitaria, así que $`\|d_k\|\le\|d_{k-1}\|+\varepsilon_k`$ por el teorema 9.18(c) aplicado a
$`\psi_{k-1}/\|\psi_{k-1}\|`$; si el espacio de Krylov tiene dimensión $`j<m`$, coincide con $`\mathcal K_j`$ y es invariante, y
el paso es exacto por la proposición 9.15. Todas las $`\|\psi_k\|=\|v\|`$ por la proposición 9.15. La cota de un paso de
longitud $`h`$ dividida por $`h`$ es $`2\beta_m\|\psi_{k-1}\|S_m(x)/x`$ con $`x=hw_T/2`$, no decreciente en $`h`$, y
$`S_m(x)/x\le\frac{x^{m-1}}{m!}(1-\frac x{m+1})^{-1}\to0`$ para $`m\ge2`$. $`\square`$

*En la biblioteca:* `krylov.expm_multiply_step` devuelve la cota (c) (`krylov.krylov_error_bound`) con cada paso.
`krylov.evolve` (por defecto `dt=None`) construye el espacio de Krylov en el estado actual y da el paso más largo $`h`$,
hasta el siguiente instante de salida, cuya cota (c) es a lo sumo $`\mathrm{tol}\cdot\|\psi_0\|\,h/(t_{\mathrm{end}}-t_0)`$; por el
corolario 9.19 cada estado devuelto está entonces a distancia a lo sumo $`\mathrm{tol}\cdot\|\psi_0\|`$ de
$`e^{-i(t-t_0)H}\psi_0`$. Esto vale en aritmética exacta: los errores de redondeo no están incluidos (crecen con el número
de pasos, y con $`\varepsilon\|H\|`$ por unidad de tiempo para un gran desplazamiento constante de $`H`$, que arrastra cada
producto matriz–vector). Por el teorema 9.18(a) los pasos no dependen de un desplazamiento constante. Con un `dt`
fijo, la suma de las cotas se compara con $`\mathrm{tol}\cdot\|\psi_0\|`$, con un aviso si es mayor.

*Cálculo:* el lema 9.17, el teorema 9.18 y el corolario 9.19 frente al propagador denso, para espectros con cúmulos,
un gran desplazamiento y una gran norma, y los pasos de `krylov.evolve` — [tutorial, sección 9.5](../../tutorial/es/09_spin_chains.md#s-9-5).

<a id="s-9-6"></a>

## 9.6 Desorden: estadística de niveles y entrelazamiento

La **cadena de Heisenberg con campo aleatorio** tiene campos independientes $`h_i`$ uniformes en $`[-W,W]`$. Conserva
$`S^z`$ (proposición 9.1), así que la estadística de niveles se calcula por sector $`N`$ (teorema 5.15 y sección 5.4:
mezclar sectores simula niveles independientes). Si dos campos difieren, $`THT^*-H=\sum_i(h_{i-1}-h_i)s^z_i\ne0`$, así
que las traslaciones dejan de ser una simetría.

*Cálculo:* el cociente medio $`\langle r\rangle`$ de la cadena con campo aleatorio, comparado con matrices del GOE y con niveles independientes (observaciones numéricas) — [tutorial, sección 9.6](../../tutorial/es/09_spin_chains.md#s-9-6).

**Entrelazamiento de fermiones libres.** Para los estados propios de la cadena XX abierta, la entropía de
entrelazamiento de los primeros $`n`$ sitios está determinada por una matriz $`n\times n`$.

<a id="r-9-16"></a>

**Teorema 9.16 (fórmula de Peschel).** Sea $`\Phi=b_{k_1}^*\cdots b_{k_N}^*|0\rangle`$ como en el teorema 9.10,
$`C_{ij}=\langle\Phi,c_i^*c_j\Phi\rangle`$ y $`C_A=(C_{ij})_{i,j<n}`$. Entonces $`0\le C_A\le1`$, y la entropía de
entrelazamiento de los sitios $`A=\{0,\dots,n-1\}`$ es, con $`\nu`$ recorriendo los valores propios de $`C_A`$ contados
con multiplicidad,
```math
S_A=\sum_\nu \eta(\nu),\qquad \eta(\nu)=-\nu\log\nu-(1-\nu)\log(1-\nu)\quad(0\log0=0).
```

*Demostración.* Para un vector real $`x`$, pongamos $`b(x)=\sum_jx_jc_j`$; entonces $`b(x)^*b(y)^*=-b(y)^*b(x)^*`$ y
$`\{b(x),b(y)^*\}=\langle x,y\rangle`$. Sean $`U\subseteq\mathbb R^L`$ el subespacio generado por los $`u_{k_a}`$, y $`P_A`$,
$`P_B=1-P_A`$ las proyecciones de coordenadas sobre $`\{0,\dots,n-1\}`$ y su complemento. Diagonalizar la forma
$`\langle P_Ax,P_Ay\rangle`$ sobre $`U`$ da una base ortonormal $`\phi_1,\dots,\phi_N`$ de $`U`$ con
$`\langle P_A\phi_a,P_A\phi_b\rangle=\nu_a\delta_{ab}`$, donde $`\nu_a=\|P_A\phi_a\|^2\in[0,1]`$, y, por tanto, $`\langle P_B\phi_a,P_B\phi_b\rangle=(1-\nu_a)\delta_{ab}`$. Como
$`(x_1,\dots,x_N)\mapsto b(x_1)^*\cdots b(x_N)^*|0\rangle`$ es multilineal y alternada,
$`\Phi=\pm\prod_a\big(\sqrt{\nu_a}f_a^*+\sqrt{1-\nu_a}\,g_a^*\big)|0\rangle`$ con $`f_a=b(P_A\phi_a/\sqrt{\nu_a})`$,
$`g_a=b(P_B\phi_a/\sqrt{1-\nu_a})`$ (omitiendo los términos con coeficiente nulo). Desarrollando,
$`\Phi=\sum_I\pm\prod_{a\in I}\sqrt{\nu_a}\prod_{a\notin I}\sqrt{1-\nu_a}\,F_I^*G_{I^c}^*|0\rangle`$ sobre los subconjuntos
$`I`$, con $`F_I^*`$, $`G_{I^c}^*`$ los productos ordenados. Para $`j\ge n`$, la cadena de signos de $`c_j^*`$ contiene todo
$`A`$: $`c_j^*=(-1)^{N_A}\otimes\tilde c_j^*`$ con operadores de Jordan–Wigner $`\tilde c_j`$ de los sitios $`\ge n`$,
mientras que para $`i<n`$, $`c_i^*=\hat c_i^*\otimes1`$. Como $`(-1)^{N_A}`$ vale $`1`$ sobre $`|0_A\rangle`$,
$`F_I^*G_{I^c}^*|0\rangle=\xi_I\otimes\chi_{I^c}`$ con $`\xi_I=\hat F_I^*|0_A\rangle`$, $`\chi_{I^c}=\tilde G_{I^c}^*|0_B\rangle`$, que son
familias ortonormales por el cálculo con las CAR del teorema 9.10. Esto es una descomposición de Schmidt
(lema 7.7) con pesos $`p_I=\prod_{a\in I}\nu_a\prod_{a\notin I}(1-\nu_a)`$ (los signos se incorporan a $`\xi_I`$), la
ley de variables de Bernoulli independientes, así que, como
$`\log p_I=\sum_{a\in I}\log\nu_a+\sum_{a\notin I}\log(1-\nu_a)`$ y $`\sum_{I\ni a}p_I=\nu_a`$, $`S_A=-\sum_Ip_I\log p_I=\sum_a\eta(\nu_a)`$. Por último, $`C_{ij}=\sum_a\phi_a(i)\phi_a(j)`$
(desarrollando $`c_i^*c_j`$ en los $`b(\phi_a)`$), así que $`C_A=\sum_a(P_A\phi_a)(P_A\phi_a)^T`$ tiene los valores propios
$`\nu_a`$ sobre los vectores ortogonales $`P_A\phi_a`$ y $`0`$ en el resto, y $`\eta(0)=0`$. $`\square`$

El teorema 9.16 se debe a I. Peschel (*J. Phys. A* 36 (2003) L205).

*Cálculo:* el teorema 9.16 para los estados fundamentales de la cadena XX abierta con $`L=12`$, frente a `SpinChain.entanglement_entropy` — [tutorial, sección 9.6](../../tutorial/es/09_spin_chains.md#s-9-6).

<a id="s-9-7"></a>

## 9.7 Con qué comparar

*Cálculo:* la comparación con otro software (QuSpin, códigos de diagonalización exacta, redes tensoriales) — [tutorial, sección 9.7](../../tutorial/es/09_spin_chains.md#s-9-7).

<a id="x-9"></a>

## Ejercicios

1. (Sectores de momento a mano) Para $`L=4`$ y $`N=2`$, halle las órbitas de las seis configuraciones bajo $`T`$ y las
   dimensiones $`d_{2,m}`$, $`m=0,\dots,3`$.
2. (Reglas de suma) Sea $`\psi`$ el estado fundamental del anillo de Heisenberg con $`L`$ par. Demuestre que
   $`E_0=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle`$ y $`\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle=0`$.
3. (Demasiada simetría) Demuestre que en el sector $`N=L/2`$ el vector con todos los coeficientes iguales a $`1`$ es
   un múltiplo de $`(S^-)^{L/2}|{\uparrow\cdots\uparrow}\rangle`$ y un vector propio del anillo de Heisenberg con
   energía $`JL/4`$, y que todo espacio de Krylov $`\mathcal K_m(H,v)`$ está contenido en cada subespacio
   $`H`$-invariante que contenga a $`v`$. ¿Qué devuelve `krylov.ground_state` cuando se inicia desde este vector?

<details>
<summary>Soluciones</summary>

1. $`0011\to1001\to1100\to0110\to0011`$ es una órbita con $`R=4`$, y $`0101\to1010\to0101`$ una con $`R=2`$. Por la
   proposición 9.5, la primera contribuye a todo $`m`$, la segunda a los $`m`$ con $`2m\equiv0\pmod4`$, es decir,
   $`m=0,2`$: $`d_{2,m}=2,1,2,1`$.

2. El estado fundamental es un singlete simple (teorema 9.9) y $`T`$ conmuta con $`H`$, así que $`T\psi=e^{ik}\psi`$.
   Como $`\mathbf s_i\cdot\mathbf s_{i+1}=T^i(\mathbf s_0\cdot\mathbf s_1)T^{*i}`$, $`\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle`$ es el
   mismo para los $`L`$ enlaces, y $`E_0=J\sum_i\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle`$.
   Como $`\mathbf S^2\psi=0`$,
   $`0=\langle\mathbf S^2\rangle=\sum_{i,j}\langle\mathbf s_i\cdot\mathbf s_j\rangle=L\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle`$, de nuevo
   por invariancia por traslaciones.

3. $`S^-`$ baja uno de los espines hacia arriba; aplicándolo $`L/2`$ veces al estado totalmente polarizado se llega a
   cada configuración con $`L/2`$ espines hacia arriba en $`(L/2)!`$ órdenes, todos con coeficiente $`1`$.
   $`|{\uparrow\cdots\uparrow}\rangle`$ tiene los $`L`$ enlaces paralelos, energía $`JL/4`$, y $`S^-`$ conmuta con $`H`$, así que
   el vector es un vector propio con la misma energía. Si $`W`$ es invariante y $`v\in W`$, entonces $`H^jv\in W`$ para
   todo $`j`$. Aquí $`\mathcal K_m=\mathbb Cv`$ para todo $`m`$, y `ground_state` devuelve $`E=JL/4`$, el extremo superior
   del espectro (cada $`\mathbf s_i\cdot\mathbf s_j=\frac12F_{ij}-\frac14\le\frac14`$ por la proposición 4.15(a), ya que $`F_{ij}`$ es
   una matriz de permutación simétrica con valores propios $`\pm1`$). Un vector inicial aleatorio evita esta trampa:
   un $`v`$ gaussiano tiene densidad, así que casi seguramente tiene una componente no nula $`P_1v`$ en el espacio
   propio de $`\lambda_1`$ (corolario A.7(c); para un $`v`$ real se aplica en $`W=\mathbb R^n`$ al subespacio estricto
   $`\{v:P_1v=0\}`$, estricto porque la matriz real $`P_1\ne0`$ tiene una columna no nula), y entonces, por el
   teorema 9.14, un espacio de Krylov $`H`$-invariante $`\mathcal K_m(H,v)`$ tiene como menor valor de Ritz $`\lambda_1`$,
   mientras que aquí el espacio invariante $`\mathcal K_1=\mathbb Cv`$ tiene el valor de Ritz $`JL/4`$.

</details>

*Cálculo:* los ejercicios 1–3, numéricamente — [tutorial, Ejercicios](../../tutorial/es/09_spin_chains.md#x-9).

**Siguiente:** el [capítulo 10](10_quantum_computing.md) trata la computación cuántica con sectores: subespacios
libres de decoherencia y subsistemas sin ruido, corrección de errores en álgebras de operadores, y mezcladores
que respetan un peso de Hamming conservado: el álgebra $`\bigoplus_wM_{\binom nw}`$ de este capítulo en un nuevo
papel.
