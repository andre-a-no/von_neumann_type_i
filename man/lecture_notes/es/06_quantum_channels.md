[English](../06_quantum_channels.md) · **Español** · [Français](../fr/06_quantum_channels.md) · [Deutsch](../de/06_quantum_channels.md) · [中文](../zh/06_quantum_channels.md) · [日本語](../ja/06_quantum_channels.md) · [Русский](../ru/06_quantum_channels.md)

# 6. Canales cuánticos

Un canal cuántico es una aplicación lineal completamente positiva que preserva la traza. La
conjugación por una unitaria, la esperanza condicional sobre el centro, la aplicación
$`\rho\mapsto\mathrm{Tr}_E(V\rho V^*)`$ dada por una isometría $`V`$ y las aplicaciones que trasladan
peso de un sector a otro son ejemplos; su lectura física (una evolución, una medición cuyo resultado
se olvida, un entorno, la pérdida de una partícula) se discute solo en los recuadros *Significado
físico*. Este capítulo demuestra las descripciones estándar de tales aplicaciones —matrices de Choi,
operadores de Kraus, la dilatación de Stinespring— para aplicaciones entre álgebras de von Neumann
de dimensión finita, y estudia después lo que tienen de especial las sumas directas: las
aplicaciones entre sectores, los duales respecto de las trazas del capítulo 3 y la parte clásica
de un canal. Los capítulos 2 y 3 ya usaron los canales de la biblioteca; aquí está su teoría.

**Aprenderá** por qué la positividad no basta, cómo pasar de matrices de Choi a operadores de
Kraus, superoperadores e isometrías y viceversa, por qué "preservar la traza" depende de la traza
en cuanto una aplicación cambia de sector, qué operadores deja fijos un canal unital y qué le
hace en promedio un canal aleatorio a un estado puro.

Las comprobaciones numéricas de este capítulo están en el [Tutorial 6](../../tutorial/es/06_quantum_channels.md).

En todo el capítulo, $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$ y
$`\mathcal N=\bigoplus_{d=1}^{C'}M_{k'_d}(\mathbb C)`$ con proyecciones centrales minimales $`P_c`$ y
$`Q_d`$; $`A=\bigoplus_cA_c`$; $`M_{k'\times k}`$ es el espacio de las matrices complejas $`k'\times k`$;
$`E_{ij}=|i\rangle\langle j|`$ son las unidades matriciales; $`\mathrm{Tr}`$ es la traza `blunt`
$`\sum_c\mathrm{Tr}\,A_c`$, y $`\mathrm{tr}_w=\sum_cw_c\mathrm{Tr}\,A_c`$ una traza fiel
(teorema 3.1). Las matrices por bloques actúan sobre $`\mathbb C^n\otimes\mathbb C^k`$, con el índice
de bloque en primer lugar.

<a id="s-6-1"></a>

## 6.1 Aplicaciones completamente positivas y teorema de Choi

Una aplicación lineal $`\Phi:\mathcal M\to\mathcal N`$ es **positiva** si $`A\ge0`$ implica $`\Phi(A)\ge0`$.
Las matrices $`n\times n`$ $`[A_{ij}]`$ con entradas en $`\mathcal M`$ forman de nuevo una suma directa de
álgebras de matrices, $`M_n(\mathcal M)\cong\bigoplus_cM_n(M_{k_c})`$ mediante
$`[A_{ij}]\mapsto([(A_{ij})_c]_{ij})_c`$, de modo que $`[A_{ij}]\ge0`$ significa que toda matriz por
bloques $`[(A_{ij})_c]\in M_{nk_c}`$ es positiva. $`\Phi`$ es **$`n`$-positiva** si
$`\Phi_n:[A_{ij}]\mapsto[\Phi(A_{ij})]`$ es positiva, y **completamente positiva** (CP) si es
$`n`$-positiva para todo $`n`$. **Preserva la traza** (TP) si
$`\mathrm{Tr}\,\Phi(A)=\mathrm{Tr}\,A`$, es **unital** si $`\Phi(1)=1`$, y es un **canal** si es
CP y TP. Usaremos repetidamente que $`X\ge0`$ implica $`TXT^*\ge0`$ para toda matriz rectangular $`T`$, ya
que $`\langle y,TXT^*y\rangle=\langle T^*y,XT^*y\rangle`$.

<a id="r-6-1"></a>

**Proposición 6.1.** Para $`k\ge2`$, la trasposición $`\theta(A)=A^T`$ sobre $`M_k`$ es positiva pero no
$`2`$-positiva.

*Demostración.* Para $`A\ge0`$, $`\langle x,A^Tx\rangle=\sum_{ij}\bar x_iA_{ji}x_j=\langle\bar x,A\bar x\rangle\ge0`$.
Con $`\Omega=e_1\otimes e_1+e_2\otimes e_2\in\mathbb C^2\otimes\mathbb C^k`$,
$`|\Omega\rangle\langle\Omega|=[E_{ij}]_{i,j=1}^2\ge0`$, y $`\theta_2`$ la lleva a
$`F=[E_{ji}]_{i,j=1}^2`$, que envía $`e_1\otimes e_2-e_2\otimes e_1`$ a su opuesto. $`\square`$

*Cálculo:* la matriz de Choi de la trasposición tiene el valor propio $`-1`$ — [tutorial, sección 6.1](../../tutorial/es/06_quantum_channels.md#s-6-1).

<a id="r-6-2"></a>

**Teorema 6.2 (Choi).** Para una aplicación lineal $`\Phi:M_k\to M_{k'}`$ son equivalentes:
(a) $`\Phi`$ es CP; (b) $`\Phi`$ es $`k`$-positiva; (c) la **matriz de Choi**
```math
J(\Phi)=\sum_{i,j=1}^kE_{ij}\otimes\Phi(E_{ij})=[\Phi(E_{ij})]_{i,j=1}^k\in M_k(M_{k'})
```
es positiva; (d) existen $`K_1,\dots,K_r\in M_{k'\times k}`$ (**operadores de Kraus**) con
$`\Phi(A)=\sum_lK_lAK_l^*`$, y se puede tomar $`r=\mathrm{rank}\,J(\Phi)\le kk'`$.

*Demostración.* (a)$`\Rightarrow`$(b) es trivial. (b)$`\Rightarrow`$(c): para
$`\Omega=\sum_ie_i\otimes e_i\in\mathbb C^k\otimes\mathbb C^k`$,
$`|\Omega\rangle\langle\Omega|=[E_{ij}]\ge0`$ y $`J(\Phi)=\Phi_k([E_{ij}])`$. (c)$`\Rightarrow`$(d): por el
teorema 1.3, $`J(\Phi)=\sum_{l=1}^rv_lv_l^*`$ con $`r=\mathrm{rank}\,J(\Phi)`$ (los vectores propios de
los valores propios no nulos, multiplicados por sus raíces cuadradas). Todo $`v\in\mathbb C^k\otimes\mathbb C^{k'}`$
se escribe $`v=\sum_ie_i\otimes Ke_i`$ para exactamente una $`K\in M_{k'\times k}`$ ($`Ke_i`$ es el
$`i`$-ésimo bloque de $`v`$), y el bloque $`(i,j)`$ de $`vv^*`$ es $`Ke_ie_j^*K^*=KE_{ij}K^*`$. Escribiendo así
cada $`v_l`$ con $`K_l`$ y comparando bloques se obtiene $`\Phi(E_{ij})=\sum_lK_lE_{ij}K_l^*`$, y de ahí (d)
por linealidad.
(d)$`\Rightarrow`$(a): $`\Phi_n([A_{ij}])=\sum_l(1_n\otimes K_l)[A_{ij}](1_n\otimes K_l)^*\ge0`$. $`\square`$

<a id="r-6-3"></a>

**Corolario 6.3.** Si $`\Phi(A)=\sum_{l=1}^sL_lAL_l^*`$, entonces $`s\ge\mathrm{rank}\,J(\Phi)`$. Así, el
número mínimo de operadores de Kraus, el **rango de Kraus**, es $`\mathrm{rank}\,J(\Phi)`$.

*Demostración.* Leyendo el último cálculo al revés, $`J(\Phi)=\sum_{l=1}^sw_lw_l^*`$ con
$`w_l=\sum_ie_i\otimes L_le_i`$, una suma de $`s`$ matrices de rango $`\le1`$. $`\square`$

*Cálculo:* la matriz de Choi y el rango de Kraus de un canal aleatorio — [tutorial, sección 6.1](../../tutorial/es/06_quantum_channels.md#s-6-1).

Una aplicación entre sumas directas es una matriz de aplicaciones entre bloques. Sean
$`\iota_c:M_{k_c}\to\mathcal M`$ la aplicación que coloca una matriz en el bloque $`c`$,
$`\pi_d:\mathcal N\to M_{k'_d}`$ la que toma el bloque $`d`$, y
$`\Phi_{dc}=\pi_d\circ\Phi\circ\iota_c:M_{k_c}\to M_{k'_d}`$. Como $`A=\sum_c\iota_c(A_c)`$,
$`\Phi(A)_d=\sum_c\Phi_{dc}(A_c)`$.

<a id="r-6-4"></a>

**Proposición 6.4.** $`\Phi:\mathcal M\to\mathcal N`$ es $`n`$-positiva si y solo si lo es cada
$`\Phi_{dc}`$. Por tanto, $`\Phi`$ es CP si y solo si existen $`K^{dc}_l\in M_{k'_d\times k_c}`$,
$`l=1,\dots,r`$, con
```math
\Phi(A)_d=\sum_c\sum_{l=1}^rK^{dc}_l\,A_c\,K^{dc*}_l ,
```
y basta con $`r=\max_{d,c}\mathrm{rank}\,J(\Phi_{dc})`$. $`\Phi`$ es TP si y solo si
$`\sum_{d,l}K^{dc*}_lK^{dc}_l=1_{k_c}`$ para todo $`c`$.

*Demostración.* Si todas las $`\Phi_{dc}`$ son $`n`$-positivas y $`[A_{ij}]\ge0`$, entonces cada
$`X^c=[(A_{ij})_c]\ge0`$, y el bloque $`d`$ de $`\Phi_n([A_{ij}])`$ es $`\sum_c(\Phi_{dc})_n(X^c)\ge0`$.
Recíprocamente, para $`[X_{ij}]\ge0`$ en $`M_n(M_{k_c})`$, se tiene $`[\iota_c(X_{ij})]\ge0`$ en
$`M_n(\mathcal M)`$ (su bloque $`c`$ es $`[X_{ij}]`$ y los demás se anulan), y el bloque $`d`$ de su imagen es
$`[\Phi_{dc}(X_{ij})]`$. La forma de Kraus se deduce del teorema 6.2 para cada par $`(d,c)`$,
completando con operadores nulos hasta un $`r`$ común. Por último, por ciclicidad,
$`\mathrm{Tr}\,\Phi(A)=\sum_c\mathrm{Tr}\big(A_c\sum_{d,l}K^{dc*}_lK^{dc}_l\big)`$, que es igual a
$`\sum_c\mathrm{Tr}\,A_c`$ para todo $`A`$ si y solo si cada suma es $`1_{k_c}`$ (tómense unidades
matriciales: $`\mathrm{Tr}(E_{ij}X)=X_{ji}`$). $`\square`$

Una aplicación lineal $`\Phi:M_k\to M_{k'}`$ cumple $`\Phi(X^*)=\Phi(X)^*`$ para todo $`X`$ si y solo si
$`J(\Phi)^*=J(\Phi)`$: el bloque $`(i,j)`$ de $`J^*`$ es $`\Phi(E_{ji})^*`$, que es igual a $`\Phi(E_{ij})`$ para
todos $`i,j`$ si y solo si $`\Phi`$ preserva los adjuntos, por linealidad.

*En la biblioteca:* las aplicaciones que preservan los sectores son un `Channel`, y las aplicaciones generales de la proposición 6.4, un `InterSectorChannel`; `Channel.choi()`, `Channel.superoperator()` y `Channel.from_superoperator` pasan de una representación a otra — [tutorial, sección 6.1](../../tutorial/es/06_quantum_channels.md#s-6-1).

> **Significado físico** (interpretación física, no un teorema de este curso). Un sistema nunca está solo: puede estar correlacionado con una ancilla
> intacta de $`n`$ niveles. Un proceso físico actúa como $`\Phi_n=\mathrm{id}_n\otimes\Phi`$ sobre el par
> y debe llevar estados a estados, de modo que $`\Phi`$ debe ser completamente positiva. La
> trasposición es positiva, pero lleva el estado entrelazado $`|\Omega\rangle\langle\Omega|/2`$ a un
> operador con valor propio $`-1/2`$. La matriz de Choi es la imagen del estado máximamente
> entrelazado: un solo estado decide la positividad completa.

<a id="s-6-2"></a>

## 6.2 Dilatación de Stinespring y libertad en los operadores de Kraus

Para $`\Phi(A)=\sum_{l=1}^rK_lAK_l^*`$ de $`M_k`$ en $`M_{k'}`$, sea $`\Phi^*(B)=\sum_lK_l^*BK_l`$; por
ciclicidad $`\mathrm{Tr}(\Phi(A)B)=\mathrm{Tr}(A\Phi^*(B))`$, de modo que $`\Phi^*`$ es el dual
de $`\Phi`$ respecto de la traza (sección 6.3) y no depende de la elección de las $`K_l`$ (si
$`\mathrm{Tr}(AX)=\mathrm{Tr}(AX')`$ para todo $`A`$, entonces $`X=X'`$: tómese $`A=E_{ji}`$). Sobre
$`\mathbb C^r\otimes\mathbb C^{k'}`$ (con el **entorno** $`\mathbb C^r`$ en primer lugar) sea
$`\mathrm{Tr}_E\,X=\sum_l(e_l^*\otimes1)X(e_l\otimes1)`$.

<a id="r-6-5"></a>

**Teorema 6.5 (Stinespring, dimensión finita).** Sea
$`V=\sum_le_l\otimes K_l:\mathbb C^k\to\mathbb C^r\otimes\mathbb C^{k'}`$. Entonces
```math
\Phi(A)=\mathrm{Tr}_E(VAV^*),\qquad\Phi^*(B)=V^*(1_r\otimes B)V,
```
y $`\Phi`$ es TP si y solo si $`V`$ es una isometría. Recíprocamente, todo $`V`$ define mediante la
primera fórmula una aplicación CP.

*Demostración.* $`(e_l^*\otimes1)V=K_l`$, luego $`\mathrm{Tr}_E(VAV^*)=\sum_lK_lAK_l^*`$;
recíprocamente, $`\mathrm{Tr}_E(VAV^*)`$ es una forma de Kraus con $`K_l=(e_l^*\otimes1)V`$. Además,
$`V^*(1\otimes B)V=\sum_{l,l'}\langle e_l,e_{l'}\rangle K_l^*BK_{l'}=\Phi^*(B)`$, y
$`V^*V=\sum_lK_l^*K_l`$, que es $`1`$ si y solo si $`\Phi`$ es TP (proposición 6.4). $`\square`$

Así, el menor $`r`$ para el que $`\Phi(A)=\mathrm{Tr}_E(VAV^*)`$ con algún
$`V:\mathbb C^k\to\mathbb C^r\otimes\mathbb C^{k'}`$ es el rango de Kraus $`\mathrm{rank}\,J(\Phi)`$: todo
tal $`V`$ da una forma de Kraus con $`r`$ operadores, luego $`r\ge\mathrm{rank}\,J(\Phi)`$ (corolario 6.3),
y $`r=\mathrm{rank}\,J(\Phi)`$ se alcanza por el teorema 6.2(d).

*Cálculo:* el teorema 6.5 para un canal aleatorio — [tutorial, sección 6.2](../../tutorial/es/06_quantum_channels.md#s-6-2).

La fórmula $`\Phi^*(B)=V^*\pi(B)V`$ con la representación $`\pi(B)=1_r\otimes B`$ es la forma en que
W. F. Stinespring (Proc. AMS 6 (1955) 211) enunció el teorema, en un marco mucho más general que el
de dimensión finita demostrado aquí (atribución histórica; lectura adicional: V. Paulsen,
*Completely Bounded Maps and Operator Algebras*, Cambridge 2002, cap. 4). Para una aplicación entre
sumas directas, sean $`\iota'_d`$ la inclusión del $`d`$-ésimo sumando de $`\bigoplus_d\mathbb C^{k'_d}`$ y
$`V_c=\sum_{d,l}e_l\otimes\iota'_dK^{dc}_l:\mathbb C^{k_c}\to\mathbb C^r\otimes\bigoplus_d\mathbb C^{k'_d}`$.
Entonces $`(e_l^*\otimes1)V_c=\sum_d\iota'_dK^{dc}_l`$, de modo que el bloque $`(d,d)`$ de
$`\sum_c\mathrm{Tr}_E(V_cA_cV_c^*)=\sum_{c,l}\big(\sum_d\iota'_dK^{dc}_l\big)A_c\big(\sum_{d'}\iota'_{d'}K^{d'c}_l\big)^*`$
es $`\sum_{c,l}K^{dc}_lA_cK^{dc*}_l=\Phi(A)_d`$; y, como $`\iota'^*_d\iota'_{d'}=\delta_{dd'}1`$,
$`V_c^*V_c=\sum_{d,l}K^{dc*}_lK^{dc}_l`$, así que $`\Phi`$ es TP si y solo si cada $`V_c`$ es una isometría
(proposición 6.4). Así es como `random_channel` y `random_inter_sector_channel` generan canales
(sección 6.5).

*En la biblioteca:* `channels.random_channel` y `channels.random_inter_sector_channel` — [tutorial, sección 6.5](../../tutorial/es/06_quantum_channels.md#s-6-5).

<a id="r-6-6"></a>

**Lema 6.6.** Si $`X,Y\in M_{p\times n}`$ y $`XX^*=YY^*`$, entonces $`Y=XW`$ para una unitaria $`W\in M_n`$.

*Demostración.* $`\|X^*z\|^2=\langle z,XX^*z\rangle=\|Y^*z\|^2`$ para todo $`z`$, luego $`X^*z\mapsto Y^*z`$
es una isometría lineal bien definida de la imagen de $`X^*`$ sobre la imagen de $`Y^*`$. Como en la
demostración del teorema 1.8, se extiende a una unitaria $`U`$ con $`UX^*=Y^*`$; tómese $`W=U^*`$. $`\square`$

<a id="r-6-7"></a>

**Teorema 6.7 (libertad unitaria).** Sean $`K_1,\dots,K_n`$ y $`L_1,\dots,L_n`$ en $`M_{k'\times k}`$
(se completa con ceros la familia más corta). Entonces $`\sum_lK_lAK_l^*=\sum_jL_jAL_j^*`$ para todo
$`A`$ si y solo si $`L_j=\sum_lu_{jl}K_l`$ para una unitaria $`u\in M_n`$.

*Demostración.* Si $`L_j=\sum_lu_{jl}K_l`$, entonces
$`\sum_jL_jAL_j^*=\sum_{l,l'}(u^*u)_{l'l}K_lAK_{l'}^*=\sum_lK_lAK_l^*`$. Recíprocamente, aplicaciones
iguales tienen matrices de Choi iguales, $`\sum_lv_lv_l^*=\sum_jw_jw_j^*`$ con
$`v_l=\sum_ie_i\otimes K_le_i`$ y $`w_j=\sum_ie_i\otimes L_je_i`$. El lema 6.6, aplicado a las matrices
de columnas $`v_l`$ y $`w_j`$, da $`w_j=\sum_lW_{lj}v_l`$; como $`K\mapsto\sum_ie_i\otimes Ke_i`$ es lineal e
inyectiva, $`L_j=\sum_lW_{lj}K_l`$, y $`u=W^T`$ es unitaria. $`\square`$

*Cálculo:* el teorema 6.7 con $`n=3`$ — [tutorial, sección 6.2](../../tutorial/es/06_quantum_channels.md#s-6-2).

> **Significado físico** (interpretación física, no un teorema de este curso). Todo canal es una interacción isométrica del sistema con un entorno en
> un estado fijo, seguida de ignorar el entorno. Medir el entorno en la base $`e_l`$ da el resultado
> $`l`$ con probabilidad $`\mathrm{Tr}\,K_l\rho K_l^*`$ y deja el estado
> $`K_l\rho K_l^*/\mathrm{Tr}(K_l\rho K_l^*)`$: los operadores de Kraus son los posibles "saltos
> cuánticos". Otra base da otros saltos pero el mismo canal (teorema 6.7), y el entorno más pequeño
> tiene dimensión igual al rango de Kraus, como se muestra tras el teorema 6.5.

<a id="s-6-3"></a>

## 6.3 Duales respecto de una traza

**Definición.** Sean $`\mathrm{tr}_w`$, $`\mathrm{tr}_v`$ trazas fieles sobre $`\mathcal M`$,
$`\mathcal N`$. El **dual** de una aplicación lineal $`\Phi:\mathcal M\to\mathcal N`$ respecto de $`(w,v)`$
es la aplicación lineal $`\Phi^\sharp:\mathcal N\to\mathcal M`$ con
```math
\mathrm{tr}_v\big(\Phi(x)\,y\big)=\mathrm{tr}_w\big(x\,\Phi^\sharp(y)\big)\qquad(x\in\mathcal M,\ y\in\mathcal N).
```
Existe y es único: para $`y`$ fijo, $`x\mapsto\mathrm{tr}_v(\Phi(x)y)`$ es una forma lineal, luego
por el teorema 3.2 es igual a $`\mathrm{tr}_w(\rho x)=\mathrm{tr}_w(x\rho)`$ para exactamente
un $`\rho=:\Phi^\sharp(y)`$, y la unicidad hace que $`\Phi^\sharp`$ sea lineal. Para las trazas `blunt`
escribimos $`\Phi^*`$. Sean $`z_w=\sum_cw_cP_c`$ y $`z_v=\sum_dv_dQ_d`$ (centrales e invertibles).

<a id="r-6-8"></a>

**Proposición 6.8.** Si $`\Phi`$ tiene operadores de Kraus $`K^{dc}_l`$ (proposición 6.4), entonces
```math
\Phi^\sharp(y)_c=\sum_d\frac{v_d}{w_c}\sum_lK^{dc*}_l\,y_d\,K^{dc}_l,\qquad\text{i.e.}\qquad
\Phi^\sharp(y)=z_w^{-1}\,\Phi^*(z_v\,y).
```
En particular, $`\Phi^\sharp`$ es CP, con operadores de Kraus $`(v_d/w_c)^{1/2}K^{dc*}_l`$.

*Demostración.* Por la ciclicidad de la traza en cada bloque,
```math
\mathrm{tr}_v(\Phi(x)y)=\sum_dv_d\sum_{c,l}\mathrm{Tr}\big(K^{dc}_lx_cK^{dc*}_ly_d\big)
=\sum_cw_c\mathrm{Tr}\Big(x_c\sum_d\frac{v_d}{w_c}\sum_lK^{dc*}_ly_dK^{dc}_l\Big),
```
y la unicidad da la fórmula. Para $`v=w=1`$ se lee $`\Phi^*(y)_c=\sum_{d,l}K^{dc*}_ly_dK^{dc}_l`$; los
factores $`v_d`$ y $`1/w_c`$ son la multiplicación por $`z_v`$ y por $`z_w^{-1}`$. $`\square`$

$`\Phi`$ **preserva la traza respecto de $`(w,v)`$** si $`\mathrm{tr}_v\,\Phi(x)=\mathrm{tr}_w(x)`$
para todo $`x`$.

<a id="r-6-9"></a>

**Corolario 6.9.** (a) $`\Phi`$ es TP respecto de $`(w,v)`$ si y solo si $`\Phi^\sharp(1)=1`$, es decir,
$`\sum_{d,l}(v_d/w_c)K^{dc*}_lK^{dc}_l=1_{k_c}`$ para todo $`c`$. (b) $`\Phi`$ es unital si y solo si
$`\Phi^\sharp`$ es TP respecto de $`(v,w)`$. (c) El dual de $`\Phi^\sharp`$ respecto de $`(v,w)`$ es $`\Phi`$.
(d) Si $`\Phi:\mathcal M\to\mathcal M`$ preserva los sectores y $`v=w`$, entonces $`\Phi^\sharp=\Phi^*`$
para todo $`w`$, y $`\Phi`$ es TP para una traza si y solo si lo es para todas.

*Demostración.* (a) $`\mathrm{tr}_v\,\Phi(x)=\mathrm{tr}_w(x\Phi^\sharp(1))`$, y
$`\mathrm{tr}_w(x(\Phi^\sharp(1)-1))=0`$ para todo $`x`$ si y solo si $`\Phi^\sharp(1)=1`$ (tómese
$`x=(\Phi^\sharp(1)-1)^*`$; $`\mathrm{tr}_w`$ es fiel); después se usa la proposición 6.8. (b)
$`\mathrm{tr}_w\,\Phi^\sharp(y)=\mathrm{tr}_v(\Phi(1)y)`$ es igual a $`\mathrm{tr}_v(y)`$
para todo $`y`$ si y solo si $`\Phi(1)=1`$. (c) $`\mathrm{tr}_w(\Phi^\sharp(y)x)=\mathrm{tr}_v(y\Phi(x))`$
por la propiedad de traza; úsese la unicidad. (d) Solo aparecen los pares $`d=c`$, para los que
$`v_d/w_c=1`$. $`\square`$

**Densidades en la imagen de Heisenberg.** Sean $`\Psi:\mathcal N\to\mathcal M`$ lineal, $`\omega`$ una
forma lineal sobre $`\mathcal M`$ con densidad $`\rho`$ para $`\mathrm{tr}_w`$ (teorema 3.2), y
$`\Psi^\sharp:\mathcal M\to\mathcal N`$ el dual de $`\Psi`$ respecto de $`(v,w)`$. Entonces
$`\omega(\Psi(y))=\mathrm{tr}_w(\Psi(y)\rho)=\mathrm{tr}_v(y\Psi^\sharp(\rho))=\mathrm{tr}_v(\Psi^\sharp(\rho)y)`$,
de modo que $`\omega\circ\Psi`$ tiene densidad $`\Psi^\sharp(\rho)`$ para $`\mathrm{tr}_v`$; si $`\Psi`$ es
unital, $`\Psi^\sharp`$ es TP respecto de $`(w,v)`$ por el corolario 6.9(b) (con los papeles de
$`\mathcal M`$ y $`\mathcal N`$ intercambiados).

Para las aplicaciones que cambian de sector, el corolario 6.9(d) falla. Sean
$`\mathcal M=\mathcal N=M_1\oplus M_2`$ (sectores $`0`$ y $`1`$) y $`\Phi`$ la aplicación llamada **pérdida de
partículas** (véase el recuadro de abajo), $`\Phi(x)=(x_0+\mathrm{Tr}\,x_1)\oplus0`$, con operadores de Kraus $`K^{00}=1`$,
$`K^{01}_l=\langle l|`$. Es TP para las trazas blunt. Para $`\mathrm{tr}_w`$ con $`w=v=(1,\tfrac12)`$,
$`\Phi^\sharp(1)=1\oplus2\cdot1_2`$ y
$`\mathrm{tr}_v\,\Phi(x)=x_0+\mathrm{Tr}\,x_1\ne x_0+\tfrac12\mathrm{Tr}\,x_1=\mathrm{tr}_w(x)`$.

*En la biblioteca:* `Channel.adjoint()` es $`\Phi^*`$ e `InterSectorChannel.adjoint(trace=...)` es $`\Phi^\sharp`$; comprobaciones de la proposición 6.8, del corolario 6.9(a) y de la pérdida de partículas — [tutorial, sección 6.3](../../tutorial/es/06_quantum_channels.md#s-6-3).

> **Significado físico** (interpretación física, no un teorema de este curso). El objeto independiente de la traza es la aplicación de Heisenberg
> $`\Psi:\mathcal N\to\mathcal M`$ sobre los observables: un estado $`\omega`$ antes del proceso se
> convierte en $`\omega\circ\Psi`$ después de él, y $`\Psi`$ es CP y unital. Las densidades dependen de
> la traza (teorema 3.2), y las densidades se transforman mediante el dual de $`\Psi`$ para las trazas
> elegidas (el párrafo que sigue al corolario 6.9). Si las álgebras actúan con multiplicidades, la
> traza física tiene pesos $`w_c=m_c`$ (sección 2.3). En el ejemplo, el sector $`0`$ se lee como "ninguna
> partícula" y el sector $`1`$ como "una partícula en dos sitios"; la aplicación elimina entonces la
> partícula. Conserva la probabilidad para la traza `blunt`, que cuenta estados, y no para
> `Tr_norm`, que pondera los sectores por igual.

<a id="s-6-4"></a>

## 6.4 Decoherencia

Los canales estándar de la biblioteca actúan dentro de cada bloque ($`p,\gamma\in[0,1]`$):

| canal | $`\Phi(\rho)_c`$ | operadores de Kraus en el bloque $`c`$ |
|---|---|---|
| `dephasing_channel` | $`(1-p)\rho_c+p\mathrm{diag}\,\rho_c`$ | $`\sqrt{1-p}\,1`$, $`\sqrt p\,E_{ii}`$ |
| `depolarizing_channel` | $`(1-p)\rho_c+p\mathrm{Tr}(\rho_c)\,1/k_c`$ | $`\sqrt{1-p}\,1`$, $`\sqrt{p/k_c}\,E_{ij}`$ |
| `center_expectation` | $`\mathrm{Tr}(\rho_c)\,1/k_c`$ | $`\sqrt{1/k_c}\,E_{ij}`$ |
| `amplitude_damping_channel` | los niveles $`i>0`$ decaen a $`0`$ | $`E_{00}+\sqrt{1-\gamma}\sum_{i>0}E_{ii}`$, $`\sqrt\gamma\,E_{0i}`$ |
| `random_mixed_unitary_channel` | $`\sum_jp_jU_j\rho U_j^*`$ | $`\sqrt{p_j}\,U_j`$ |

Las formas de Kraus se comprueban por cálculo directo, p. ej.
$`\sum_{ij}\frac1kE_{ij}\rho E_{ji}=\frac1k\sum_{ij}\rho_{jj}E_{ii}`$. Los cinco son TP, y todos salvo
el amortiguamiento de amplitud con $`\gamma>0`$ (si algún $`k_c\ge2`$) son unitales (proposición 6.4: calcúlense
$`\sum_lK_l^*K_l`$ y $`\sum_lK_lK_l^*`$; para el amortiguamiento de amplitud en un bloque de tamaño $`k`$
valen $`1`$ y $`1+(k-1)\gamma E_{00}-\gamma\sum_{i>0}E_{ii}`$, y $`\sum_jp_j=1`$ para el canal unitario
mixto). Explícitamente, el amortiguamiento de amplitud lleva
$`\rho_{00}\mapsto\rho_{00}+\gamma\sum_{i>0}\rho_{ii}`$, $`\rho_{ii}\mapsto(1-\gamma)\rho_{ii}`$ y
$`\rho_{ij}\mapsto(1-\gamma)\rho_{ij}`$ para $`i\ne j`$ ambos positivos, y
$`\rho_{0j}\mapsto\sqrt{1-\gamma}\,\rho_{0j}`$, $`\rho_{j0}\mapsto\sqrt{1-\gamma}\,\rho_{j0}`$ para $`j>0`$. La
esperanza condicional $`E`$ sobre el centro (proposición 2.8) es la despolarización con $`p=1`$. A diferencia de la medición del sector,
$`X\mapsto\sum_cP_cXP_c`$, que es la identidad sobre $`\mathcal M`$ (sección 3.3), solo conserva los $`p_c`$.

*En la biblioteca:* los cinco canales de la tabla son funciones de `channels` — [tutorial, sección 6.4](../../tutorial/es/06_quantum_channels.md#s-6-4).

<a id="r-6-10"></a>

**Teorema 6.10 (puntos fijos de los canales unitales).** Sea $`\Phi(X)=\sum_lK_lXK_l^*`$ con
$`K_l\in\mathcal M`$ TP y unital. Entonces $`\Phi(X)=X`$ si y solo si $`X`$ conmuta con todas las $`K_l`$ y
$`K_l^*`$. En particular, los puntos fijos forman una subálgebra de von Neumann de $`\mathcal M`$.

*Demostración.* Si $`X`$ conmuta con todas las $`K_l`$, entonces $`\Phi(X)=X\sum_lK_lK_l^*=X`$.
Recíprocamente, $`\Phi(X^*)=\Phi(X)^*`$, de modo que las partes autoadjuntas $`(X+X^*)/2`$ y
$`(X-X^*)/2i`$ de un punto fijo son puntos fijos; sea $`X=X^*`$. Con $`Y_l=K_lX-XK_l`$,
```math
\sum_lY_lY_l^*=\sum_l\big(K_lX^2K_l^*-K_lXK_l^*X-XK_lXK_l^*+XK_lK_l^*X\big)
=\Phi(X^2)-\Phi(X)X-X\Phi(X)+X\Phi(1)X=\Phi(X^2)-X^2,
```
usando $`Y_l^*=XK_l^*-K_l^*X`$, $`\Phi(X)=X`$ y $`\Phi(1)=1`$. Su traza es $`0`$ porque $`\Phi`$ es TP. Como
$`\mathrm{Tr}\,Y_lY_l^*=\|Y_l\|_2^2\ge0`$, toda $`Y_l=0`$;
tomando adjuntos, $`X`$ conmuta también con $`K_l^*`$. Estos $`X`$ forman el conmutante de un conjunto
autoadjunto, intersecado con $`\mathcal M`$: una $`*`$-álgebra que contiene a $`1`$, y por tanto un álgebra
de von Neumann (teorema 2.2). $`\square`$

Así, el desfase ($`p>0`$) conserva exactamente la diagonal de cada bloque, y la despolarización
($`p>0`$) y $`E`$ conservan exactamente el centro. Sin unitalidad el teorema falla: el amortiguamiento
de amplitud fija $`E_{00}`$, que no conmuta con el salto $`E_{01}`$. Para $`0<\gamma\le1`$ y $`k\ge2`$, sus puntos fijos
en un bloque y el conmutante de sus operadores de Kraus tienen ambos dimensión $`1`$, pero son
$`\mathbb CE_{00}`$ y $`\mathbb C1`$: por la forma explícita de arriba, un punto fijo tiene $`\rho_{ij}=0`$
para $`i\ne j`$ y $`\rho_{ii}=0`$ para $`i>0`$ (estas entradas se multiplican por $`\sqrt{1-\gamma}\ne1`$ o
$`1-\gamma\ne1`$), y todo múltiplo de $`E_{00}`$ es fijo; y lo que conmuta con todas las $`E_{0i}`$, $`E_{i0}`$
($`i>0`$) conmuta con sus productos $`E_{ij}=E_{i0}E_{0j}`$ y $`E_{00}=E_{01}E_{10}`$, y por tanto con $`M_k`$,
así que es un múltiplo de $`1`$ (capítulo 2, ejercicio 1).

*Cálculo:* puntos fijos y conmutantes de los cinco canales — [tutorial, sección 6.4](../../tutorial/es/06_quantum_channels.md#s-6-4).

<a id="r-6-11"></a>

**Teorema 6.11 (los canales unitales aumentan la entropía).** Sea $`\Phi:\mathcal M\to\mathcal N`$
positiva, TP y unital, y sean $`\rho`$ y $`\Phi(\rho)`$ estados con valores propios $`\lambda_i`$ y $`\mu_j`$.
Entonces $`\mu=\Delta\lambda`$ para una matriz doblemente estocástica $`\Delta`$, de modo que
$`\sum_jf(\mu_j)\ge\sum_if(\lambda_i)`$ para toda $`f`$ cóncava sobre $`[0,1]`$. En particular,
$`S(\Phi(\rho))\ge S(\rho)`$ y $`\mathrm{Tr}\,\Phi(\rho)^2\le\mathrm{Tr}\,\rho^2`$.

*Demostración.* Diagonalicemos cada bloque: $`\rho=\sum_i\lambda_i|e_i\rangle\langle e_i|`$ y
$`\Phi(\rho)=\sum_j\mu_j|f_j\rangle\langle f_j|`$ con bases ortonormales de $`\bigoplus_c\mathbb C^{k_c}`$
y de $`\bigoplus_d\mathbb C^{k'_d}`$ cuyos vectores están cada uno en un solo sumando, de modo que
$`|e_i\rangle\langle e_i|\in\mathcal M`$. Entonces $`\mu_j=\langle f_j,\Phi(\rho)f_j\rangle=\sum_i\Delta_{ji}\lambda_i`$
con $`\Delta_{ji}=\langle f_j,\Phi(|e_i\rangle\langle e_i|)f_j\rangle\ge0`$,
$`\sum_j\Delta_{ji}=\mathrm{Tr}\,\Phi(|e_i\rangle\langle e_i|)=1`$ y
$`\sum_i\Delta_{ji}=\langle f_j,\Phi(1)f_j\rangle=1`$ (así que $`\Delta`$ es cuadrada: súmense todas las
entradas). Por la desigualdad de Jensen finita $`f(\sum_i\Delta_{ji}\lambda_i)\ge\sum_i\Delta_{ji}f(\lambda_i)`$
(por inducción sobre el número de términos a partir de la definición de concavidad, como en la
solución del ejercicio 2 del capítulo 3),
$`\sum_jf(\mu_j)\ge\sum_j\sum_i\Delta_{ji}f(\lambda_i)=\sum_if(\lambda_i)`$. Tómense $`f(t)=-t\log t`$
($`f(0)=0`$) y $`f(t)=-t^2`$, que son cóncavas sobre $`[0,1]`$ por el lema 3.4: $`t\log t`$ y $`t^2`$ son
continuas sobre $`[0,1]`$ con segundas derivadas $`1/t>0`$ y $`2>0`$ sobre $`(0,1)`$. $`\square`$

*Cálculo:* la entropía y la pureza bajo los canales unitales, y la caída de la entropía bajo el amortiguamiento de amplitud — [tutorial, sección 6.4](../../tutorial/es/06_quantum_channels.md#s-6-4).

La parte clásica de una aplicación entre sectores es su **matriz de transición**
$`T_{dc}(\rho)=\mathrm{Tr}\,\Phi_{dc}(\rho_c)/\mathrm{Tr}\,\rho_c`$ ($`\rho_c\ne0`$); el método
`transition_matrix()` de la biblioteca es $`T=T(1)`$, para el estado máximamente mixto de cada sector
de entrada.

<a id="r-6-12"></a>

**Proposición 6.12.** Sean $`\Phi:\mathcal M\to\mathcal N`$ positiva y TP, $`\rho`$ un estado con
$`p_c=\mathrm{Tr}\,\rho_c>0`$, y $`q_d=\mathrm{Tr}\,\Phi(\rho)_d`$. Entonces $`T(\rho)`$ es estocástica
por columnas y $`q=T(\rho)p`$. Con las esperanzas condicionales $`E_{\mathcal M}`$, $`E_{\mathcal N}`$ sobre
los centros (proposición 2.8),
```math
E_{\mathcal N}\Phi E_{\mathcal M}(\rho)=\bigoplus_d\frac{(Tp)_d}{k'_d}\,1_{k'_d},\qquad T=T(1).
```

*Demostración.* $`T_{dc}\ge0`$ porque $`\Phi_{dc}`$ es positiva, y
$`\sum_dT_{dc}=\mathrm{Tr}\,\Phi(\iota_c\rho_c)/\mathrm{Tr}\,\rho_c=1`$. Además,
$`q_d=\sum_c\mathrm{Tr}\,\Phi_{dc}(\rho_c)=\sum_cT_{dc}(\rho)p_c`$. Por último,
$`E_{\mathcal M}\rho=\bigoplus_cp_c1_{k_c}/k_c`$, de modo que la traza de
$`\Phi(E_{\mathcal M}\rho)_d=\sum_cp_c\Phi_{dc}(1_{k_c}/k_c)`$ es $`\sum_cT_{dc}p_c`$, y $`E_{\mathcal N}`$
sustituye el bloque $`d`$ por su traza multiplicada por $`1/k'_d`$. $`\square`$

*En la biblioteca:* la matriz de transición es `InterSectorChannel.transition_matrix()`; una comprobación de la proposición 6.12 — [tutorial, sección 6.4](../../tutorial/es/06_quantum_channels.md#s-6-4).

**Coherencias de un cúbit.** Sobre $`M_2`$, sean $`D_p`$ el desfase y $`A_\gamma`$ el amortiguamiento de
amplitud. Por la tabla y la forma explícita de arriba, $`D_p`$ multiplica $`\rho_{01}`$ y $`\rho_{10}`$ por
$`1-p`$ y fija la diagonal, mientras que $`A_\gamma`$ multiplica $`\rho_{01}`$ y $`\rho_{10}`$ por
$`\sqrt{1-\gamma}`$ y lleva la diagonal a $`(\rho_{00}+\gamma\rho_{11},(1-\gamma)\rho_{11})`$,
independientemente de las coherencias. Por tanto, $`D_pA_\gamma=A_\gamma D_p`$: ambas composiciones
multiplican las coherencias por $`(1-p)\sqrt{1-\gamma}`$ y actúan sobre la diagonal como $`A_\gamma`$. En
particular, para $`T_1,T_\varphi>0`$, $`t\ge0`$, $`1-p=e^{-t/T_\varphi}`$ y $`1-\gamma=e^{-t/T_1}`$, la composición
multiplica las coherencias por $`e^{-t/T_2}`$ con
```math
\frac1{T_2}=\frac1{2T_1}+\frac1{T_\varphi},
```
y estas familias se componen aditivamente en $`t`$: $`D_{p(t)}D_{p(s)}=D_{p(t+s)}`$ y
$`A_{\gamma(t)}A_{\gamma(s)}=A_{\gamma(t+s)}`$, ya que los factores $`e^{-t/T}`$ se multiplican y la traza se
preserva.

*Cálculo:* el factor de coherencia y $`1/T_2=1/(2T_1)+1/T_\varphi`$ — [tutorial, sección 6.4](../../tutorial/es/06_quantum_channels.md#s-6-4).

> **Significado físico** (interpretación física, no un teorema de este curso). Para un cúbit, el
> desfase modela el desfase puro con tiempo $`T_\varphi`$ y el amortiguamiento de amplitud el
> decaimiento del nivel excitado con tiempo $`T_1`$; el cálculo anterior da entonces el tiempo de
> coherencia $`T_2`$ de un cúbit sometido a ambos. El teorema 6.10 es la forma
> algebraica de la información libre de decoherencia: lo que conmuta con todos los operadores de
> salto sobrevive (capítulo 10). La proposición 6.12 dice que, visto a través del centro, una
> aplicación entre sectores es una matriz de Markov clásica sobre las cargas.

<a id="s-6-5"></a>

## 6.5 Canales aleatorios

Una **isometría Haar-aleatoria** $`V:\mathbb C^k\to\mathbb C^N`$ es, como en la biblioteca
(`random_isometry`), el factor $`Q`$ de $`G=QR`$, donde $`G`$ es una matriz de Ginibre $`N\times k`$
(entradas independientes con densidad $`\pi^{-1}e^{-|z|^2}`$), $`Q^*Q=1`$, y $`R`$ es triangular superior
con diagonal positiva (en la CPU la biblioteca toma una factorización de LAPACK $`G=Q'R'`$ y multiplica la
$`i`$-ésima columna de $`Q'`$ por $`R'_{ii}/|R'_{ii}|`$, lo que da esta factorización, como en el
teorema 5.3, siempre que todos los $`R'_{ii}\ne0`$, es decir, siempre que $`G`$ tenga rango $`k`$, lo que
ocurre casi seguramente por la demostración del lema 6.13; para lotes grandes en una GPU usa
Gram–Schmidt, que la da directamente). `random_channel(alg, r)` genera tales
$`V_c:\mathbb C^{k_c}\to\mathbb C^r\otimes\mathbb C^{k_c}`$ independientemente en cada bloque y usa los
operadores de Kraus del teorema 6.5.

*En la biblioteca:* `random_isometry` y `channels.random_channel` — [tutorial, sección 6.5](../../tutorial/es/06_quantum_channels.md#s-6-5).

<a id="r-6-13"></a>

**Lema 6.13.** $`V`$ tiene la ley de las primeras $`k`$ columnas de una unitaria de Haar de $`U(N)`$. En
particular, para toda unitaria $`U\in M_N`$, $`UV`$ tiene la misma ley que $`V`$.

*Demostración.* Sean $`Z`$ una matriz de Ginibre $`N\times N`$ cuyas primeras $`k`$ columnas son $`G`$, y
$`E\in M_{N\times k}`$ las primeras $`k`$ columnas de $`1_N`$. Por el teorema 5.3, $`Z`$ es invertible casi
seguramente y $`Z=Q_ZR_Z`$ con $`Q_Z`$ distribuida según Haar. Como $`R_Z`$ es triangular superior,
$`G=ZE=(Q_ZE)(E^*R_ZE)`$ es una factorización del tipo requerido, y es la única: si $`G=QR=Q'R'`$,
entonces $`R^*R=G^*G=R'^*R'`$, luego $`R'R^{-1}`$ es unitaria y triangular superior con diagonal
positiva, y por tanto igual a $`1`$ (demostración del lema 5.2). Así, $`V=Q_ZE`$, y $`UV=(UQ_Z)E`$ tiene
la ley de $`V`$ porque $`UQ_Z`$ tiene la de $`Q_Z`$. $`\square`$

<a id="r-6-14"></a>

**Lema 6.14.** Sea $`u`$ un vector unitario aleatorio de $`\mathbb C^N`$ cuya ley es invariante bajo
todas las unitarias. Entonces $`\mathbb E[u_a\bar u_b]=\delta_{ab}/N`$ y
```math
\mathbb E\big[u_a\bar u_b\,u_c\bar u_d\big]=\frac{\delta_{ab}\delta_{cd}+\delta_{ad}\delta_{cb}}{N(N+1)} .
```

*Demostración.* La invariancia bajo $`\mathrm{diag}(e^{i\theta_1},\dots,e^{i\theta_N})`$
multiplica los momentos por $`e^{i(\theta_a-\theta_b)}`$ y $`e^{i(\theta_a-\theta_b+\theta_c-\theta_d)}`$,
de modo que se anulan salvo si $`a=b`$, resp. $`\{a,c\}=\{b,d\}`$ como multiconjuntos. Por la invariancia
bajo matrices de permutación, $`\mathbb E|u_a|^2=1/N`$, y solo quedan $`\alpha=\mathbb E|u_1|^4`$ y
$`\beta=\mathbb E|u_1|^2|u_2|^2`$. En primer lugar, $`\sum_c|u_c|^2=1`$ da
$`\alpha+(N-1)\beta=\mathbb E|u_1|^2=1/N`$. En segundo lugar, $`(u_1+u_2)/\sqrt2`$ es la primera
coordenada de $`Uu`$ para una unitaria $`U`$, luego $`\mathbb E|u_1+u_2|^4=4\alpha`$. Con
$`s=|u_1|^2+|u_2|^2`$ y $`z=u_1\bar u_2`$, $`|u_1+u_2|^4=s^2+2s(z+\bar z)+(z+\bar z)^2`$; las
invariancias $`u_1\mapsto-u_1`$ y $`u_1\mapsto iu_1`$ dan $`\mathbb E[s(z+\bar z)]=0`$ y
$`\mathbb Ez^2=0`$, luego $`4\alpha=(2\alpha+2\beta)+2\beta`$, es decir, $`\alpha=2\beta`$. Por tanto,
$`\beta=1/(N(N+1))`$ y $`\alpha=2\beta`$, que es la fórmula en todos los casos. (Esto ha usado $`N\ge2`$;
para $`N=1`$, $`|u_1|=1`$ y ambos lados valen $`1`$.) $`\square`$

<a id="r-6-15"></a>

**Teorema 6.15 (pureza de salida de un canal aleatorio).** Sean $`\psi`$ un estado puro en un sector
de tamaño $`k`$ y $`\Phi`$ generado por `random_channel` con rango de Kraus $`r`$. Entonces
```math
\mathbb E\mathrm{Tr}\,\Phi(\psi)^2=\frac{k+r}{kr+1},\qquad\mathbb E\,\Phi(\psi)=\frac1k1_k .
```

*Demostración.* En el sector de $`\psi`$, $`\Phi(\psi)=\mathrm{Tr}_E|u\rangle\langle u|`$ con
$`u=V\psi`$, un vector unitario cuya ley es unitariamente invariante por el lema 6.13. Indexemos
$`\mathbb C^r\otimes\mathbb C^k`$ por pares $`(e,x)`$. Entonces $`\Phi(\psi)_{xy}=\sum_eu_{ex}\bar u_{ey}`$,
cuya media es $`r\delta_{xy}/(kr)`$, y
$`\mathrm{Tr}\,\Phi(\psi)^2=\sum_{x,y,e,f}u_{ex}\bar u_{ey}u_{fy}\bar u_{fx}`$. Por el lema 6.14 con
$`N=kr`$, el primer término impone $`x=y`$ ($`kr^2`$ tuplas de índices) y el segundo $`e=f`$ ($`k^2r`$
tuplas), de modo que la media es $`(kr^2+k^2r)/(kr(kr+1))`$. $`\square`$

*Cálculo:* el teorema 6.15 mediante Monte Carlo para $`r=1,2,4`$ — [tutorial, sección 6.5](../../tutorial/es/06_quantum_channels.md#s-6-5).

Para $`r=1`$ el canal es unitario y la pureza es $`1`$; cuando $`r\to\infty`$ tiende a $`1/k`$, la pureza
del estado máximamente mixto. La demostración solo ha usado la invariancia unitaria de la ley de
$`u`$, así que también da la pureza media de Lubkin del estado reducido de un estado puro aleatorio
$`u=g/\|g\|`$, con $`g`$ gaussiano (cuya ley es unitariamente invariante por el teorema A.20(b)), que
el capítulo 7, teorema 7.11, demuestra con integrales gaussianas; los números $`\alpha,\beta`$ del
lema 6.14 son los del corolario 5.11 para las entradas de una unitaria de Haar, y el corolario A.21
los recalcula para $`u=g/\|g\|`$.

> **Significado físico** (interpretación física, no un teorema de este curso). Un canal aleatorio modela un ruido genérico, sin estructura, procedente de
> un entorno de dimensión $`r`$. Ya con $`r=4`$ la pureza de un estado puro de un cútrit baja en promedio
> a $`7/13\approx0.54`$ (teorema 6.15), frente a $`1`$ sin ruido y al mínimo $`1/3`$ (por
> Cauchy–Schwarz, $`\mathrm{Tr}\,\rho^2\ge(\mathrm{Tr}\,\rho)^2/k`$). En un sistema con superselección el
> ruido actúa en cada sector por separado, y la dimensión relevante es la del sector, no la de todo el
> espacio.

<a id="s-6-6"></a>

## 6.6 Con qué comparar

*En la biblioteca:* comparación con otro software — [tutorial, sección 6.6](../../tutorial/es/06_quantum_channels.md#s-6-6).

<a id="x-6"></a>

## Ejercicios

1. (Werner–Holevo) Demuestre que $`\Phi(\rho)=\big(\mathrm{Tr}(\rho)1-\rho^T\big)/(k-1)`$ sobre
   $`M_k`$, $`k\ge2`$, es un canal unital con rango de Kraus $`k(k-1)/2`$, aunque esté construido a partir
   de la trasposición.
2. Demuestre que para $`0<p\le1`$ el canal de desfase sobre $`M_k`$ tiene rango de Kraus $`k`$, uno menos
   que el número de operadores de Kraus que almacena la biblioteca.
3. Sean $`X`$ el desplazamiento cíclico sobre $`\mathbb C^k`$ y $`Z=\mathrm{diag}(1,\omega,\dots,\omega^{k-1})`$,
   $`\omega=e^{2\pi i/k}`$. Demuestre que $`\frac1{k^2}\sum_{a,b=0}^{k-1}X^aZ^b\rho Z^{-b}X^{-a}=\mathrm{Tr}(\rho)1/k`$:
   la esperanza condicional sobre el centro de $`M_k`$ es un canal unitario mixto.
4. Demuestre que $`(\Psi\circ\Phi)^\sharp=\Phi^\sharp\circ\Psi^\sharp`$ cuando se usa la misma traza en
   el álgebra intermedia.

*Cálculo:* los cuatro ejercicios, numéricamente — [tutorial, Ejercicios](../../tutorial/es/06_quantum_channels.md#x-6).

<details>
<summary>Soluciones</summary>

1. Con $`F`$ el operador de intercambio, $`J(\theta)=\sum_{ij}E_{ij}\otimes E_{ji}=F`$ y $`J(\rho\mapsto\mathrm{Tr}(\rho)1)=1\otimes1`$,
   luego $`J(\Phi)=(1-F)/(k-1)=2P_-/(k-1)`$, donde $`P_-=(1-F)/2`$ es la proyección sobre los vectores
   antisimétricos ($`F=F^*`$, $`F^2=1`$), de rango $`\mathrm{Tr}\,P_-=(k^2-\mathrm{Tr}\,F)/2`$ con
   $`\mathrm{Tr}\,F=\sum_{ij}\langle e_i\otimes e_j,e_j\otimes e_i\rangle=k`$. Así, $`J\ge0`$ con rango $`k(k-1)/2`$ (teorema 6.2,
   corolario 6.3), y $`\mathrm{Tr}\,\Phi(\rho)=\mathrm{Tr}\,\rho`$, $`\Phi(1)=1`$ directamente.

2. $`J=(1-p)|\Omega\rangle\langle\Omega|+p\sum_iE_{ii}\otimes E_{ii}`$ se anula fuera del subespacio
   generado por los $`e_i\otimes e_i`$, y sobre este espacio de dimensión $`k`$ es la matriz
   $`(1-p)\mathbf 1\mathbf 1^T+p1`$ ($`\mathbf 1`$ el vector de unos), con valores propios $`p`$ y
   $`(1-p)k+p`$, todos positivos. Así, $`\mathrm{rank}\,J=k`$, mientras que la biblioteca almacena
   $`1+k_{\max}`$ operadores.

3. $`\mathrm{Tr}((X^aZ^b)^*X^{a'}Z^{b'})=k\,\delta_{aa'}\delta_{bb'}`$: $`Z^{-b}X^{a'-a}Z^{b'}`$ tiene
   diagonal nula salvo si $`a=a'`$, y $`\mathrm{Tr}\,Z^{b'-b}=\sum_j\omega^{j(b'-b)}=k\delta_{bb'}`$. Así,
   las $`B_{ab}=X^aZ^b/\sqrt k`$ forman una base ortonormal de $`M_k`$ para $`\mathrm{Tr}(A^*B)`$, y dos
   bases ortonormales cualesquiera están relacionadas por una unitaria; por el teorema 6.7,
   $`\sum_mB_m\rho B_m^*`$ es el mismo para todas las bases ortonormales. Para la base $`E_{ab}`$ vale
   $`\sum_{ab}E_{ab}\rho E_{ba}=\mathrm{Tr}(\rho)1`$. Divídase entre $`k`$.

4. $`\mathrm{tr}(\Psi(\Phi(x))z)=\mathrm{tr}(\Phi(x)\Psi^\sharp(z))=\mathrm{tr}(x\,\Phi^\sharp(\Psi^\sharp(z)))`$
   para todos los $`x`$ y $`z`$; la unicidad del dual da la afirmación.

</details>

**Siguiente:** el [capítulo 7](07_tensor_products.md) trata los productos tensoriales y el entrelazamiento:
la traza parcial —un canal entre álgebras distintas cuyo dual es la inclusión $`A\mapsto A\otimes1`$—,
la fórmula de Page y la fusión de cargas.
