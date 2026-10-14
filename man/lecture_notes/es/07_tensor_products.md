[English](../07_tensor_products.md) · **Español** · [Français](../fr/07_tensor_products.md) · [Deutsch](../de/07_tensor_products.md) · [中文](../zh/07_tensor_products.md) · [日本語](../ja/07_tensor_products.md) · [Русский](../ru/07_tensor_products.md)

# 7. Productos tensoriales y entrelazamiento

Este capítulo estudia el producto tensorial de dos álgebras con sectores (físicamente, el álgebra
de dos sistemas considerados conjuntamente; véanse los recuadros *Significado físico*). Con sectores
hay dos productos: el simple, en el que cada parte conserva su propia carga, y el *fusionado*, en el
que solo se conserva la carga total. Demostramos la estructura de ambos, definimos la traza parcial
como el dual de una inclusión y estudiamos el entrelazamiento.

**Aprenderá** cómo se multiplican los sectores, por qué la traza
parcial es un canal y, normalizada, una esperanza condicional de índice finito, cuál es el estado
reducido de un estado que conserva la carga, cómo el entrelazamiento se divide en una parte de
número y una parte configuracional, y la pureza y la entropía medias de los estados aleatorios. Las
comprobaciones numéricas de este capítulo están en el [Tutorial 7](../../tutorial/es/07_tensor_products.md).

En todo el capítulo, $`\mathcal M_1=\bigoplus_{c=1}^{C_1}M_{k_c}(\mathbb C)`$ actúa sobre
$`H_1=\bigoplus_c\mathbb C^{k_c}`$ con proyecciones centrales minimales $`P_c`$, y
$`\mathcal M_2=\bigoplus_{d=1}^{C_2}M_{l_d}(\mathbb C)`$ actúa sobre $`H_2=\bigoplus_d\mathbb C^{l_d}`$
con proyecciones centrales minimales $`Q_d`$. $`\mathrm{Tr}`$ es la traza del espacio sobre el que
actúa un operador; sobre $`\mathcal M_i`$ es la traza `blunt` del capítulo 3. Para $`X`$ sobre
$`\mathbb C^k\otimes\mathbb C^l`$ escribimos $`X_{(a,b),(a',b')}=\langle e_a\otimes e_b,X\,e_{a'}\otimes e_{b'}\rangle`$,
e identificamos $`\mathbb C^k\otimes\mathbb C^l`$ con $`\mathbb C^{kl}`$ mediante $`e_a\otimes e_b\mapsto e_{(a-1)l+b}`$
(el convenio de Kronecker), de modo que $`(A\otimes B)_{(a,b),(a',b')}=A_{aa'}B_{bb'}`$.

<a id="s-7-1"></a>

## 7.1 El producto tensorial de dos álgebras

Para $`*`$-álgebras $`\mathcal M_1`$ sobre $`H_1`$ y $`\mathcal M_2`$ sobre $`H_2`$, el **producto tensorial**
$`\mathcal M_1\otimes\mathcal M_2`$ es el subespacio generado por los operadores $`A\otimes B`$ sobre
$`H_1\otimes H_2`$, con $`A\in\mathcal M_1`$, $`B\in\mathcal M_2`$. Como
$`(A\otimes B)(A'\otimes B')=AA'\otimes BB'`$, $`(A\otimes B)^*=A^*\otimes B^*`$ y $`1\otimes1=1`$, es
una $`*`$-álgebra y, por tanto, un álgebra de von Neumann (teorema 2.2).

<a id="r-7-1"></a>

**Teorema 7.1.** (a) $`H_1\otimes H_2=\bigoplus_{c,d}\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$, suma
ortogonal de las imágenes de las proyecciones $`P_c\otimes Q_d`$, y
```math
\mathcal M_1\otimes\mathcal M_2=\bigoplus_{c,d}M_{k_c}\otimes M_{l_d}\cong\bigoplus_{c,d}M_{k_cl_d}(\mathbb C):
```
está formada por todos los operadores $`X=\sum_{c,d}(P_c\otimes Q_d)X(P_c\otimes Q_d)`$, con bloques
arbitrarios sobre los $`\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$. (b) Sus proyecciones centrales
minimales son las $`C_1C_2`$ proyecciones $`P_c\otimes Q_d`$, y
$`\dim(\mathcal M_1\otimes\mathcal M_2)=\dim\mathcal M_1\dim\mathcal M_2`$.
(c) Para trazas $`\mathrm{tr}_w`$ sobre $`\mathcal M_1`$ y $`\mathrm{tr}_v`$ sobre $`\mathcal M_2`$ (teorema 3.1),
$`\mathrm{tr}_w(A)\mathrm{tr}_v(B)=\mathrm{tr}_{w\otimes v}(A\otimes B)`$ con
$`(w\otimes v)_{cd}=w_cv_d`$. En particular, `Tr_blunt`, `Tr_norm` y `tau_vN` del producto son los
productos de las mismas trazas de los factores.

*Demostración.* (a) Las $`P_c\otimes Q_d`$ son proyecciones mutuamente ortogonales de suma
$`1\otimes1`$, y la imagen de $`P_c\otimes Q_d`$ es $`\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$. Como $`P_c`$ y
$`Q_d`$ son centrales, $`(P_c\otimes Q_d)(A\otimes B)=AP_c\otimes BQ_d=(A\otimes B)(P_c\otimes Q_d)`$, de
modo que todo elemento de $`\mathcal M_1\otimes\mathcal M_2`$ es diagonal por bloques. Recíprocamente,
con las unidades matriciales $`E^c_{aa'}\in\mathcal M_1`$ y $`F^d_{bb'}\in\mathcal M_2`$ de los bloques
$`c`$ y $`d`$, se tiene $`E^c_{aa'}\otimes F^d_{bb'}=|e_a\otimes e_b\rangle\langle e_{a'}\otimes e_{b'}|`$ dentro de
$`\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$; estas son las unidades matriciales de $`M_{k_cl_d}(\mathbb C)`$,
así que todo operador con soporte en el bloque $`(c,d)`$ pertenece a $`\mathcal M_1\otimes\mathcal M_2`$.
(b) El centro de una suma directa de álgebras de matrices completas está generado por las
proyecciones sobre los sumandos (capítulo 2, ejercicio 1), y
$`\sum_{c,d}k_c^2l_d^2=\sum_ck_c^2\sum_dl_d^2`$. (c)
$`\mathrm{Tr}(X\otimes Y)=\sum_{a,b}X_{aa}Y_{bb}=\mathrm{Tr}\,X\mathrm{Tr}\,Y`$, luego
$`\mathrm{tr}_{w\otimes v}(A\otimes B)=\sum_{c,d}w_cv_d\mathrm{Tr}\,A_c\mathrm{Tr}\,B_d=\mathrm{tr}_w(A)\mathrm{tr}_v(B)`$.
Los pesos $`1`$, $`1/k_c`$, $`1/(C_1k_c)`$ de las tres trazas (sección 3.1) se multiplican y dan $`1`$,
$`1/(k_cl_d)`$, $`1/(C_1C_2k_cl_d)`$, los pesos de las mismas trazas sobre el producto. $`\square`$

*En la biblioteca:* el producto es `tensor_product` y $`A\otimes B`$ es `kron`; una comprobación del teorema 7.1 sobre $`(\mathbb C\oplus M_2)\otimes(M_2\oplus M_3)`$ — [tutorial, sección 7.1](../../tutorial/es/07_tensor_products.md#s-7-1).

> **Significado físico** (interpretación física, no un teorema de este curso). Si cada uno de dos
> sistemas obedece su propia regla de superselección
> —dos cajas que no pueden intercambiar partículas—, el sistema conjunto tiene un sector por cada
> *par* de cargas, y ningún observable conecta pares distintos. Un sector conjunto cuenta con el
> producto de los pesos (teorema 7.1(c)).

<a id="s-7-2"></a>

## 7.2 La traza parcial

**Definición.** Sea $`\rho\in\mathcal M_1\otimes\mathcal M_2`$. La **traza parcial**
$`\mathrm{Tr}_2\,\rho\in\mathcal M_1`$ es el elemento tal que
```math
\mathrm{Tr}\big((\mathrm{Tr}_2\,\rho)\,A\big)=\mathrm{Tr}\big(\rho\,(A\otimes1)\big)\qquad\text{for all }A\in\mathcal M_1 .
```
Existe y es única por el teorema 3.2, aplicado a la forma lineal
$`A\mapsto\mathrm{Tr}(\rho(A\otimes1))`$ y a la traza fiel $`\mathrm{Tr}`$ de $`\mathcal M_1`$: la
traza parcial es el dual, respecto de las trazas `blunt`, de la inclusión $`A\mapsto A\otimes1`$.
$`\mathrm{Tr}_1`$ se define del mismo modo; para un estado $`\rho`$, $`\mathrm{Tr}_2\,\rho`$ y
$`\mathrm{Tr}_1\,\rho`$ son los **estados reducidos**. Para las álgebras completas $`\mathcal M_i=B(H_i)`$ se
trata de la **traza parcial del espacio de Hilbert** $`\mathrm{Tr}_{H_2}:B(H_1\otimes H_2)\to B(H_1)`$.

<a id="r-7-2"></a>

**Proposición 7.2.** (a) Para $`\rho=\bigoplus_{c,d}\rho_{cd}`$,
```math
(\mathrm{Tr}_2\,\rho)_c=\sum_d\mathrm{tr}_{l_d}\,\rho_{cd},\qquad
(\mathrm{tr}_l\,X)_{aa'}=\sum_{b=1}^lX_{(a,b),(a',b)}\quad(X\text{ on }\mathbb C^k\otimes\mathbb C^l),
```
y $`\mathrm{Tr}_2\,\rho=\mathrm{Tr}_{H_2}\,\rho`$, donde $`(\mathrm{Tr}_{H_2}\,X)_{ii'}=\sum_jX_{(i,j),(i',j)}`$.
(b) $`\mathrm{Tr}_2(A\otimes B)=(\mathrm{Tr}\,B)\,A`$.
(c) $`\mathrm{Tr}_2\big((A\otimes1)X(A'\otimes1)\big)=A\,(\mathrm{Tr}_2\,X)\,A'`$ para $`A,A'\in\mathcal M_1`$.
(d) $`\mathrm{Tr}(\mathrm{Tr}_2\,X)=\mathrm{Tr}\,X`$, y $`\mathrm{Tr}_2`$ transforma estados en estados.

*Demostración.* (a) Con $`A=\bigoplus_cA_c`$, se tiene $`A\otimes1=\bigoplus_{c,d}A_c\otimes1_{l_d}`$ y
```math
\mathrm{Tr}\big(\rho(A\otimes1)\big)=\sum_{c,d}\sum_{a,a',b}(\rho_{cd})_{(a,b),(a',b)}(A_c)_{a'a}
=\sum_c\mathrm{Tr}\Big(\Big(\sum_d\mathrm{tr}_{l_d}\,\rho_{cd}\Big)A_c\Big),
```
así que el elemento de la derecha satisface la definición y es $`\mathrm{Tr}_2\,\rho`$. El mismo
cálculo para $`B(H_i)`$ da la fórmula de $`\mathrm{Tr}_{H_2}`$; para $`\rho`$ diagonal por bloques, las
entradas $`\rho_{(i,j),(i',j)}`$ con $`i,i'`$ en bloques distintos se anulan, de modo que
$`\mathrm{Tr}_{H_2}\,\rho`$ tiene los bloques de (a). (b) $`\mathrm{tr}_l(X\otimes Y)=(\mathrm{Tr}\,Y)X`$.
(c) Para $`B\in\mathcal M_1`$,
$`\mathrm{Tr}\big(\mathrm{Tr}_2((A\otimes1)X(A'\otimes1))B\big)=\mathrm{Tr}\big(X(A'BA\otimes1)\big) =\mathrm{Tr}\big(A(\mathrm{Tr}_2\,X)A'B\big)`$; se concluye por unicidad. (d) Tomemos $`A=1`$. Si $`\rho`$ es
un estado, $`A\mapsto\mathrm{Tr}(\rho(A\otimes1))`$ es positiva, ya que $`A^*A\otimes1=(A\otimes1)^*(A\otimes1)`$,
y vale $`1`$ en $`A=1`$, así que su densidad $`\mathrm{Tr}_2\,\rho`$ es un estado (teorema 3.2). $`\square`$

*En la biblioteca:* $`\mathrm{Tr}_2`$ es `partial_trace`; una comprobación de la proposición 7.2(a), (b) frente a `numpy.einsum` — [tutorial, sección 7.2](../../tutorial/es/07_tensor_products.md#s-7-2).

Los canales son el tema del capítulo 6; aquí solo necesitamos esta definición: una aplicación
lineal $`\Phi`$ es **completamente positiva** si $`\mathrm{id}_n\otimes\Phi`$ es positiva sobre
$`M_n(\mathbb C)\otimes(\cdot)`$ para todo $`n`$, y es un **canal** si además conserva la traza.

<a id="r-7-3"></a>

**Teorema 7.3.** Sea $`j`$ un índice que recorre la base canónica de $`H_2`$ y sea
$`K_j=1\otimes\langle j|:H_1\otimes H_2\to H_1`$, $`K_j(x\otimes y)=\langle j,y\rangle x`$. Entonces
```math
\mathrm{Tr}_{H_2}\,X=\sum_jK_jXK_j^*,\qquad\sum_jK_j^*K_j=1,\qquad\sum_jK_j^*AK_j=A\otimes1 .
```
Por tanto, $`\mathrm{Tr}_2`$ es un canal de $`\mathcal M_1\otimes\mathcal M_2`$ en $`\mathcal M_1`$, y su
dual respecto de las trazas `blunt` es $`A\mapsto A\otimes1`$.

*Demostración.* $`K_j^*x=x\otimes e_j`$, luego $`\langle e_a,K_jXK_j^*e_{a'}\rangle=X_{(a,j),(a',j)}`$; sumando
sobre $`j`$ se obtiene la proposición 7.2(a). Además,
$`\sum_jK_j^*K_j(x\otimes y)=\sum_j\langle j,y\rangle x\otimes e_j=x\otimes y`$
y $`\sum_jK_j^*AK_j=\sum_jA\otimes|j\rangle\langle j|=A\otimes1`$. Para $`Y\ge0`$,
$`(\mathrm{id}_n\otimes\mathrm{Tr}_{H_2})(Y)=\sum_j(1_n\otimes K_j)Y(1_n\otimes K_j)^*\ge0`$, porque
$`\langle x,TYT^*x\rangle=\langle T^*x,YT^*x\rangle\ge0`$; esto se mantiene en la subálgebra
$`\mathcal M_1\otimes\mathcal M_2`$. Por último,
$`\mathrm{Tr}\,\sum_jK_jXK_j^*A=\mathrm{Tr}\,X\sum_jK_j^*AK_j`$, lo que para $`A=1`$ es la conservación
de la traza y, en general, la dualidad. $`\square`$

*En la biblioteca:* el canal del teorema 7.3 es `composite.partial_trace_channel`, con su dual — [tutorial, sección 7.2](../../tutorial/es/07_tensor_products.md#s-7-2).

> **Significado físico** (interpretación física, no un teorema de este curso). El estado reducido
> $`\mathrm{Tr}_2\,\rho`$ es todo lo que averigua un observador que solo puede medir el primer sistema; la definición dice que reproduce todos los
> valores esperados de los observables $`A\otimes1`$. Los estados conjuntos con estados reducidos
> iguales difieren en sus correlaciones.

**La traza parcial como esperanza condicional.** Sea $`\mathcal M_2=M_m(\mathbb C)`$, de modo que los
sectores de $`\mathcal M_1\otimes M_m`$ son los de $`\mathcal M_1`$ (teorema 7.1), sean
$`\iota(A)=A\otimes1_m`$ y $`\mathcal E=\frac1m\mathrm{Tr}_2:\mathcal M_1\otimes M_m\to\mathcal M_1`$,
y sean $`e_{ab}`$ las unidades matriciales de $`M_m`$.

<a id="r-7-4"></a>

**Proposición 7.4.** (a) $`\iota\circ\mathcal E`$ es una esperanza condicional de
$`\mathcal M_1\otimes M_m`$ sobre $`\mathcal M_1\otimes1`$ (sección 2.6) y, para toda traza fiel
$`\mathrm{tr}_w`$ de $`\mathcal M_1\otimes M_m`$, es la única que conserva $`\mathrm{tr}_w`$. (b) Los
$`m^2`$ elementos $`u_{ab}=\sqrt m\,(1\otimes e_{ab})`$ satisfacen
$`X=\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*X))`$ para todo $`X`$, y
$`\sum_{a,b}u_{ab}u_{ab}^*=m^2\,1`$. (c) Toda familia $`(v_k)`$ con
$`X=\sum_kv_k\,\iota(\mathcal E(v_k^*X))`$ para todo $`X`$ cumple $`\sum_kv_kv_k^*=m^2\,1`$.

*Demostración.* (a) $`\mathcal E`$ es positiva (teorema 7.3) y también lo es $`\iota`$
($`T^*T\otimes1=(T\otimes1)^*(T\otimes1)`$); $`\mathcal E(A\otimes1)=A`$ por la proposición 7.2(b), y la
propiedad de bimódulo es la proposición 7.2(c). Escribamos también $`\mathrm{tr}_w`$ para la traza
con los mismos pesos $`w_c`$ sobre $`\mathcal M_1`$; entonces $`\mathrm{tr}_w(Y\otimes1)=m\mathrm{tr}_w(Y)`$. La
definición de $`\mathrm{Tr}_2`$ aplicada a $`AP_c`$, multiplicada por $`w_c`$ y sumada sobre $`c`$, da
$`\mathrm{tr}_w(X(A\otimes1))=\mathrm{tr}_w((\mathrm{Tr}_2\,X)A)`$, y $`A=1`$ da
$`\mathrm{tr}_w(\iota\mathcal E(X))=\mathrm{tr}_w(X)`$. Si $`F=\iota\circ F_0`$ es otra esperanza
condicional sobre $`\mathcal M_1\otimes1`$ que conserva $`\mathrm{tr}_w`$, entonces, como
$`F(X(A\otimes1))=F(X)(A\otimes1)`$,
$`\mathrm{tr}_w(X(A\otimes1))=\mathrm{tr}_w(F(X)(A\otimes1))=m\mathrm{tr}_w(F_0(X)A)`$, luego
$`\mathrm{tr}_w((\mathrm{Tr}_2\,X-mF_0(X))A)=0`$ para todo $`A\in\mathcal M_1`$, y $`F_0=\mathcal E`$ por
fidelidad (tómese $`A=(\mathrm{Tr}_2\,X-mF_0(X))^*`$). (b) Escribamos $`X=\sum_{a,b}X_{ab}\otimes e_{ab}`$ con
$`X_{ab}\in\mathcal M_1`$ (teorema 7.1(a)). Por la proposición 7.2(b),
$`\mathcal E(u_{ab}^*X)=\frac{\sqrt m}m\sum_{b'}\mathrm{Tr}_2(X_{ab'}\otimes e_{bb'})=X_{ab}/\sqrt m`$, luego
$`u_{ab}\,\iota(\mathcal E(u_{ab}^*X))=X_{ab}\otimes e_{ab}`$; y
$`\sum_{a,b}u_{ab}u_{ab}^*=m\sum_{a,b}1\otimes e_{aa}=m^21`$. (c) $`\mathrm{Tr}_2(Y^*)=(\mathrm{Tr}_2\,Y)^*`$ por
la proposición 7.2(a), así que el adjunto del desarrollo de $`Y=u_{ab}`$ es
$`u_{ab}^*=\sum_k\iota(\mathcal E(u_{ab}^*v_k))\,v_k^*`$ y, por (b),
$`\sum_{a,b}u_{ab}u_{ab}^*=\sum_k\big(\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*v_k))\big)v_k^*=\sum_kv_kv_k^*`$. $`\square`$

Una familia como la de (c) es una *cuasibase* de $`\mathcal E`$, y $`\sum_kv_kv_k^*`$, que por (c) no
depende de ella y por (b) existe, es el **índice** $`\mathrm{Ind}\,\mathcal E=m^2`$ de $`\mathcal E`$ (una
definición; el nombre sigue a Y. Watatani, Mem. Amer. Math. Soc. 83, no. 424, 1990). La
**matriz de inclusión** de una inclusión de $`\bigoplus_cM_{k_c}`$ en $`\bigoplus_sM_{n_s}`$ tiene las entradas
$`\Lambda_{sc}`$, el número de copias del bloque $`c`$ de la subálgebra en el bloque $`s`$ de la mayor
(el capítulo 9 la usa para las torres de álgebras). Para un factor, $`M_k\otimes1\subset M_k\otimes M_m`$:
ordenando la base $`e_a\otimes e_b`$ de $`\mathbb C^{km}`$ primero por $`b`$, $`A\otimes1_m=\mathrm{diag}(A,\dots,A)`$ con
$`m`$ copias, así que $`\Lambda=(m)`$, y $`m^2=\|\Lambda\|^2=\dim(M_k\otimes M_m)/\dim M_k`$. (Lectura adicional,
no usada aquí: V. F. R. Jones, Invent. Math. 72, 1, 1983, sobre el índice de los subfactores, y
F. M. Goodman, P. de la Harpe y V. F. R. Jones, *Coxeter graphs and towers of algebras*, Springer,
1989, sobre las matrices de inclusión.) Para $`\mathcal M_1`$ general, la inclusión es la suma directa de
las $`C_1`$ inclusiones $`M_{k_c}\otimes1\subset M_{k_c}\otimes M_m`$, $`\Lambda=m1_{C_1}`$.

*En la biblioteca:* `inclusions.from_tensor_product`, con la esperanza condicional, la cuasibase y el índice de la proposición 7.4 — [tutorial, sección 7.2](../../tutorial/es/07_tensor_products.md#s-7-2).

<a id="s-7-3"></a>

## 7.3 Fusión de cargas

Supongamos ahora que el sector $`c`$ de $`\mathcal M_1`$ lleva una carga $`q_c`$, el sector $`d`$ de
$`\mathcal M_2`$ una carga $`r_d`$, y que una **regla de fusión** $`f`$ da la carga total $`f(q_c,r_d)`$: la
suma para números de partículas, la suma módulo $`2`$ para la paridad fermiónica, o cualquier otra
aplicación. Pongamos
```math
H_q=\bigoplus_{f(q_c,r_d)=q}\mathbb C^{k_c}\otimes\mathbb C^{l_d},\qquad
\Pi_q=\sum_{f(q_c,r_d)=q}P_c\otimes Q_d,\qquad n_q=\sum_{f(q_c,r_d)=q}k_cl_d .
```
El **producto fusionado** $`\mathcal M_1\otimes_f\mathcal M_2=\{\Pi_q\}'`$ está formado por los operadores
sobre $`H_1\otimes H_2`$ que conservan la carga total.

<a id="r-7-5"></a>

**Proposición 7.5.** (a) $`\mathcal M_1\otimes_f\mathcal M_2=\bigoplus_qB(H_q)\cong\bigoplus_qM_{n_q}(\mathbb C)`$,
con proyecciones centrales minimales $`\Pi_q`$. (b) $`\mathcal M_1\otimes\mathcal M_2\subseteq\mathcal M_1\otimes_f\mathcal M_2`$,
con igualdad si y solo si cada $`H_q`$ contiene un único par $`(c,d)`$. (c) Si las cargas son números
reales y $`f`$ es la suma, entonces $`\mathcal M_1\otimes_f\mathcal M_2=\{N\}'`$ con
$`N=N_1\otimes1+1\otimes N_2`$, $`N_1=\sum_cq_cP_c`$, $`N_2=\sum_dr_dQ_d`$.

*Demostración.* (a) Cada par está exactamente en un $`H_q`$, así que las $`\Pi_q`$ son proyecciones
mutuamente ortogonales de suma $`1`$, y $`X`$ conmuta con todas ellas si y solo si
$`X=\sum_q\Pi_qX\Pi_q`$, con bloques arbitrarios en $`B(H_q)`$; el centro es como en el teorema 7.1(b).
(b) Todo elemento de $`\mathcal M_1\otimes\mathcal M_2`$ conmuta con cada $`P_c\otimes Q_d`$ y, por tanto,
con sus sumas $`\Pi_q`$. Si cada $`H_q`$ es un único $`\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$, ambas álgebras
son $`\bigoplus_{c,d}B(\mathbb C^{k_c}\otimes\mathbb C^{l_d})`$; si $`H_q`$ contiene dos pares, $`B(H_q)`$
contiene operadores que llevan el bloque de un par al del otro, y que no son diagonales por bloques.
(c) $`N=\sum_{c,d}(q_c+r_d)P_c\otimes Q_d=\sum_qq\,\Pi_q`$, de modo que las $`\Pi_q`$ no nulas son las
proyecciones espectrales de $`N`$. Son polinomios en $`N`$ (demostración del teorema 1.5(a)) y $`N`$ es
una combinación de ellas, así que $`X`$ conmuta con $`N`$ si y solo si conmuta con cada $`\Pi_q`$. $`\square`$

**Ejemplo: bosones de núcleo duro.** Sobre $`L`$ sitios, $`(\mathbb C^2)^{\otimes L}`$ con el
número de ocupación $`N=\sum_in^{(i)}`$, $`n=\mathrm{diag}(0,1)`$, el subespacio propio $`N=q`$ está generado
por los $`\binom Lq`$ vectores de la base producto con $`q`$ sitios ocupados, así que los operadores que
conservan el número $`\{N\}'`$ forman $`\bigoplus_qM_{\binom Lq}`$ (proposición 4.10(2)): para dos sitios
$`M_1\oplus M_2\oplus M_1`$ (cargas $`0,1,2`$), para tres sitios $`M_1\oplus M_3\oplus M_3\oplus M_1`$. Por la
proposición 7.5(c), el producto fusionado es $`\{N\}'`$ para el número de ocupación total de cinco
sitios, salvo el orden de la base, y en efecto $`n_q=\sum_c\binom2c\binom3{q-c}=\binom5q`$ (identidad de
Vandermonde: un subconjunto de $`q`$ elementos de cinco sitios consta de $`c`$ de los dos primeros y
$`q-c`$ de los tres últimos).

*En la biblioteca:* el producto fusionado es `fused_tensor_product`; los bosones de núcleo duro y la proposición 7.5 — [tutorial, sección 7.3](../../tutorial/es/07_tensor_products.md#s-7-3).

El álgebra fusionada tiene dimensión $`\sum_q\binom5q^2=\binom{10}5=252`$ (de nuevo Vandermonde), el
producto simple $`(1+4+1)(1+9+9+1)=120`$ (teorema 7.1(b)). Por el teorema 7.1(a) y la proposición 7.5(b),
el producto simple consta de los elementos del fusionado que son diagonales por bloques para los
pares $`(c,d)`$; las $`132`$ dimensiones restantes están generadas por los operadores que llevan un bloque
de par $`(c,d)`$ a otro bloque de par $`(c',d')`$ del mismo $`H_q`$, $`c'\ne c`$ (luego $`d'\ne d`$): cambian la
carga de cada parte conservando la total, es decir, aquí mueven partículas entre las dos partes.

Un estado de $`\mathcal M_1\otimes_f\mathcal M_2`$ puede tener coherencias entre bloques de pares
distintos, así que la proposición 7.2 no se aplica. Sea $`E_1(X)=\sum_cP_cXP_c`$ sobre $`B(H_1)`$ (la
esperanza condicional de la sección 2.6); su imagen es $`\mathcal M_1`$, y
$`\mathrm{Tr}(E_1(X)A)=\mathrm{Tr}(XA)`$ para $`A\in\mathcal M_1`$, porque $`A=\sum_cP_cAP_c`$.

<a id="r-7-6"></a>

**Teorema 7.6.** Sea $`\rho`$ la densidad de un estado $`\omega`$ sobre $`\mathcal M_1\otimes_f\mathcal M_2`$.
La densidad de su restricción $`A\mapsto\omega(A\otimes1)`$ a $`\mathcal M_1`$ es
```math
\rho_1=E_1\big(\mathrm{Tr}_{H_2}\,\rho\big).
```
El bloque $`(c,c')`$ de $`\mathrm{Tr}_{H_2}\,\rho`$ es $`\sum_d\mathrm{tr}_{l_d}\,\rho_{(c,d),(c',d)}`$, sumado sobre
los $`d`$ con $`f(q_c,r_d)=f(q_{c'},r_d)`$, donde $`\rho_{(c,d),(c',d)}`$ es el bloque de $`\rho`$ del par
$`(c',d)`$ al par $`(c,d)`$. Por tanto, $`E_1`$ actúa trivialmente sobre todos los $`\mathrm{Tr}_{H_2}\,\rho`$
si y solo si no existen $`c\neq c'`$ y $`d`$ con $`f(q_c,r_d)=f(q_{c'},r_d)`$; por ejemplo, si $`f`$ es la
suma y los $`q_c`$ son distintos dos a dos.

*Demostración.* $`A\otimes1\in\mathcal M_1\otimes_f\mathcal M_2`$ (proposición 7.5(b)) y
$`\omega(A\otimes1)=\mathrm{Tr}(\rho(A\otimes1))=\mathrm{Tr}\big((\mathrm{Tr}_{H_2}\,\rho)A\big) =\mathrm{Tr}\big(E_1(\mathrm{Tr}_{H_2}\,\rho)A\big)`$. Como $`E_1(\mathrm{Tr}_{H_2}\,\rho)\in\mathcal M_1`$, es
la densidad por la unicidad del teorema 3.2. La fórmula por bloques es la proposición 7.2(a) para
$`\mathrm{Tr}_{H_2}`$, con $`i`$ en el bloque $`c`$, $`i'`$ en el bloque $`c'`$ y $`j`$ en el bloque $`d`$; el bloque
de $`\rho`$ entre $`(c',d)`$ y $`(c,d)`$ se anula salvo que ambos pares estén en el mismo $`H_q`$. Si
existen $`c\neq c'`$ y $`d`$ con $`f(q_c,r_d)=f(q_{c'},r_d)=q`$, tomemos vectores unitarios
$`x\in\mathbb C^{k_c}`$, $`x'\in\mathbb C^{k_{c'}}`$, $`y\in\mathbb C^{l_d}`$ y
$`\psi=(x+x')\otimes y/\sqrt2\in H_q`$. Entonces $`|\psi\rangle\langle\psi|\in B(H_q)`$, y
$`\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|=\frac12|x+x'\rangle\langle x+x'|`$ tiene el bloque no nulo
$`\frac12|x\rangle\langle x'|`$, que $`E_1`$ elimina. $`\square`$

Por la fórmula por bloques del teorema 7.6, $`(\rho_1)_c=\sum_d\mathrm{tr}_{l_d}\,\rho_{(c,d),(c,d)}`$.

*En la biblioteca:* `partial_trace` y `composite.partial_trace_channel` calculan $`\rho_1`$ para cualquier regla de fusión; una comprobación del teorema 7.6 para la suma y para la regla de paridad — [tutorial, sección 7.3](../../tutorial/es/07_tensor_products.md#s-7-3).

> **Significado físico** (interpretación física, no un teorema de este curso). En la descripción
> BCS de un superconductor solo se conserva la paridad
> del número de electrones, de modo que el estado global puede superponer los números de
> partículas $`0`$ y $`2`$ de una región. Un observador cuyos observables conservan el número de
> partículas de la región no puede ver esta coherencia: para él, el estado reducido es el estado
> diagonal por bloques del teorema 7.6 (el contenido matemático de este recuadro), una mezcla de
> números de partículas.

<a id="s-7-4"></a>

## 7.4 Entrelazamiento de estados puros

<a id="r-7-7"></a>

**Lema 7.7 (descomposición de Schmidt).** Todo vector unitario $`\psi\in\mathbb C^a\otimes\mathbb C^b`$
puede escribirse $`\psi=\sum_{i=1}^rs_i\,u_i\otimes v_i`$ con $`s_i>0`$, $`\sum_is_i^2=1`$ y
$`u_i\in\mathbb C^a`$, $`v_i\in\mathbb C^b`$ ortonormales. Entonces
$`\mathrm{Tr}_2|\psi\rangle\langle\psi|=\sum_is_i^2|u_i\rangle\langle u_i|`$
y $`\mathrm{Tr}_1|\psi\rangle\langle\psi|=\sum_is_i^2|v_i\rangle\langle v_i|`$: los dos estados reducidos
tienen los mismos valores propios no nulos $`s_i^2`$.

*Demostración.* Sea $`\Psi`$ la matriz $`a\times b`$ dada por $`\Psi_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,\psi\rangle`$.
Por el teorema espectral, $`\Psi\Psi^*=\sum_{i=1}^a\lambda_i|u_i\rangle\langle u_i|`$ con una base
ortonormal $`u_i`$ y $`\lambda_i=\|\Psi^*u_i\|^2\ge0`$, ordenados de modo que $`\lambda_i>0`$ exactamente
para $`i\le r`$. Pongamos $`s_i=\sqrt{\lambda_i}`$ y $`w_i=\Psi^*u_i/s_i`$ para $`i\le r`$; entonces
$`\langle w_i,w_j\rangle=\langle u_i,\Psi\Psi^*u_j\rangle/(s_is_j)=\delta_{ij}`$, y $`\Psi^*u_i=0`$ para $`i>r`$.
Por tanto $`\Psi=\sum_i|u_i\rangle\langle u_i|\Psi=\sum_{i\le r}|u_i\rangle\langle\Psi^*u_i|=\sum_{i\le r}s_i|u_i\rangle\langle w_i|`$,
es decir, $`\Psi_{\alpha\beta}=\sum_is_i(u_i)_\alpha\overline{(w_i)_\beta}`$: es la afirmación con
$`v_i=\bar w_i`$, y $`\sum_is_i^2=\mathrm{Tr}\,\Psi\Psi^*=\|\psi\|^2=1`$. Por último,
$`|\psi\rangle\langle\psi|=\sum_{i,j}s_is_j|u_i\rangle\langle u_j|\otimes|v_i\rangle\langle v_j|`$, y la
proposición 7.2(b) da los estados reducidos. $`\square`$

La **entropía de entrelazamiento** de $`\psi`$ es
$`E(\psi)=S(\mathrm{Tr}_2|\psi\rangle\langle\psi|)=S(\mathrm{Tr}_1|\psi\rangle\langle\psi|)=-\sum_is_i^2\log s_i^2`$.
Se anula si y solo si $`r=1`$ (cada término $`-s_i^2\log s_i^2`$ es $`\ge0`$, y vale $`0`$ solo para
$`s_i^2=1`$), es decir, si y solo si $`\psi=u\otimes v`$ (entonces $`\Psi=uv^T`$ tiene rango $`1`$, y
$`r=\mathrm{rank}\,\Psi\Psi^*=\mathrm{rank}\,\Psi`$), y $`E(\psi)\le\log\min(a,b)`$, con igualdad si y solo si
el estado reducido del factor más pequeño es máximamente mixto (corolario 3.6, aplicado a ambos
estados reducidos).

Con cargas tomamos un vector unitario $`\psi`$ con una carga total definida, es decir, en un sector
$`H_q`$ del producto fusionado (una superposición de distintos $`H_q`$ da sobre
$`\mathcal M_1\otimes_f\mathcal M_2`$ el mismo estado que una mezcla, capítulo 3, ejercicio 3); escribamos
$`\psi=\sum_{(c,d)}\psi_{cd}`$ con $`\psi_{cd}\in\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$.

<a id="r-7-8"></a>

**Proposición 7.8 (entrelazamiento resuelto por simetría).** Sea $`f`$ la suma, y supongamos que los
$`q_c`$ son distintos dos a dos y que los $`r_d`$ son distintos dos a dos, de modo que para cada $`c`$ a
lo sumo un $`d=d(c)`$ cumple $`q_c+r_d=q`$. Pongamos $`\psi_c=\psi_{c\,d(c)}`$ (o $`0`$),
$`p_c=\|\psi_c\|^2`$ y, para $`p_c>0`$, $`\hat\psi_c=\psi_c/\sqrt{p_c}`$. Entonces
$`\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|`$ pertenece a $`\mathcal M_1`$ (conmuta con $`N_1`$) y es igual a
$`\bigoplus_cp_c\mathrm{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|`$; el estado reducido sobre $`\mathcal M_2`$
tiene la probabilidad de sector $`p_c`$ en $`d(c)`$; y
```math
E(\psi)=H(p)+\sum_cp_c\,E(\hat\psi_c),\qquad H(p)=-\sum_cp_c\log p_c .
```

*Demostración.* Por la fórmula por bloques del teorema 7.6, el bloque $`(c,c')`$ de
$`\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|`$ es la suma de $`\mathrm{tr}_{l_d}|\psi_{cd}\rangle\langle\psi_{c'd}|`$
sobre los $`d`$ con $`q_c+r_d=q=q_{c'}+r_d`$, lo que obliga a que $`c=c'`$ y $`d=d(c)`$. Así, el estado
reducido es $`\bigoplus_c\mathrm{Tr}_2|\psi_c\rangle\langle\psi_c|`$, con trazas de bloque $`p_c`$;
intercambiando los papeles de los factores se obtiene el enunciado para $`\mathcal M_2`$. Por el
lema 7.7, $`E(\psi)`$ es la entropía de este estado, y el teorema 3.3 con
$`\sigma_c=\mathrm{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|`$, $`S(\sigma_c)=E(\hat\psi_c)`$, da la
fórmula. $`\square`$

$`H(p)`$ es el **entrelazamiento de número** y $`\sum_cp_cE(\hat\psi_c)`$ el **entrelazamiento
configuracional**.

*Cálculo:* el lema 7.7 y la proposición 7.8 para un estado aleatorio de dos partículas en cinco sitios — [tutorial, sección 7.4](../../tutorial/es/07_tensor_products.md#s-7-4).

> **Significado físico** (interpretación física, no un teorema de este curso). Las partes sujetas
> a la regla de superselección del número de partículas solo pueden usar como recurso la parte configuracional: el *entrelazamiento de partículas* de
> H. M. Wiseman y J. A. Vaccaro (Phys. Rev. Lett. 91, 097902, 2003). $`H(p)`$ se mide en átomos fríos
> como la *entropía de número* (A. Lukin et al., Science 364, 256, 2019).

<a id="s-7-5"></a>

## 7.5 Estados aleatorios y la fórmula de Page

Sea $`g\in\mathbb C^N`$ con coordenadas independientes cuyas partes real e imaginaria son
independientes con ley $`N(0,\tfrac12)`$ (un vector gaussiano complejo estándar, proposición A.16(a)), de
modo que $`\mathbb E|g_i|^2=1`$; $`g\neq0`$ casi seguramente (teorema A.20). Un **estado puro aleatorio
de Haar** es $`\psi=g/\|g\|`$. Su ley es invariante unitaria: $`Ug`$ tiene la ley de $`g`$ para toda unitaria
$`U`$ (proposición A.16(a)), y $`Ug/\|Ug\|=U\psi`$. Por la proposición A.22 es la única ley invariante
unitaria sobre la esfera unidad y, como se ve en la demostración de esa proposición, es la ley de la
primera columna de una unitaria de Haar (teorema 5.3).

<a id="r-7-9"></a>

**Lema 7.9.** $`\|g\|`$ y $`\psi=g/\|g\|`$ son independientes.

*Demostración.* (Es parte del teorema A.20(b).) Sean $`h\ge0`$ y $`\phi\ge0`$ funciones de Borel sobre
$`(0,\infty)`$ y sobre la esfera unidad $`S^{2N-1}`$ de $`\mathbb C^N\cong\mathbb R^{2N}`$, respectivamente. Como $`g`$
tiene la densidad $`\pi^{-N}e^{-|x|^2}`$ (proposición A.16(a)), las coordenadas polares (teorema A.18 con
$`n=2N`$ y la medida $`\mu`$ allí definida) dan
```math
\mathbb E\big[h(\|g\|)\phi(\psi)\big]=\pi^{-N}\!\int_{\mathbb R^{2N}}\!h(|x|)e^{-|x|^2}\phi\Big(\frac x{|x|}\Big)dx
=2N\pi^{-N}\int_0^\infty h(t)e^{-t^2}t^{2N-1}dt\int_{S^{2N-1}}\phi\,d\mu=:\alpha(h)\beta(\phi).
```
Tomando $`\phi=1`$, luego $`h=1`$ y luego ambas, se obtiene $`\mathbb Eh(\|g\|)=\alpha(h)\beta(1)`$,
$`\mathbb E\phi(\psi)=\alpha(1)\beta(\phi)`$ y $`\alpha(1)\beta(1)=1`$, luego
$`\mathbb E[h(\|g\|)\phi(\psi)]=\mathbb Eh(\|g\|)\,\mathbb E\phi(\psi)`$. $`\square`$

<a id="r-7-10"></a>

**Lema 7.10.** $`\mathbb E\|g\|^4=N(N+1)`$. Si $`G`$ es una matriz $`a\times b`$ con entradas independientes
distribuidas como los $`g_i`$ (una matriz *de Ginibre*), entonces $`\mathbb E\mathrm{Tr}(GG^*)^2=ab(a+b)`$.

*Demostración.* Con $`g_i=x+iy`$ y $`x,y`$ independientes con ley $`N(0,\frac12)`$, $`\mathbb Ex^4=3(\mathbb Ex^2)^2=\frac34`$
(teorema A.24: tres emparejamientos), luego
$`\mathbb E|g_i|^4=\mathbb Ex^4+2\mathbb Ex^2\mathbb Ey^2+\mathbb Ey^4=\tfrac34+\tfrac12+\tfrac34=2`$ y
$`\mathbb E\|g\|^4=\sum_i\mathbb E|g_i|^4+\sum_{i\ne j}\mathbb E|g_i|^2\mathbb E|g_j|^2=2N+N(N-1)`$. Además,
$`\mathrm{Tr}(GG^*)^2=\sum_{i,i',j,j'}G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}`$. Las entradas y sus
conjugadas son combinaciones lineales complejas de las coordenadas de un mismo vector gaussiano
centrado (las partes real e imaginaria de todas las entradas), con
$`\mathbb E[G_{ij}G_{kl}]=\mathbb E[\bar G_{ij}\bar G_{kl}]=0`$ y $`\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}`$
(proposición A.16(a)). Por la fórmula de Isserlis–Wick (teorema A.24, forma compleja, con los tres
emparejamientos de los cuatro factores),
$`\mathbb E[G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}]=\delta_{ii'}\delta_{i'i}+0+\delta_{jj'}\delta_{jj'}=\delta_{ii'}+\delta_{jj'}`$.
Sumando sobre $`i,i'\le a`$ y $`j,j'\le b`$ se obtiene $`ab^2+a^2b`$. $`\square`$

<a id="r-7-11"></a>

**Teorema 7.11 (E. Lubkin, J. Math. Phys. 19, 1028, 1978).** Para $`\psi\in\mathbb C^a\otimes\mathbb C^b`$
aleatorio de Haar y $`\rho_A=\mathrm{Tr}_2|\psi\rangle\langle\psi|`$,
```math
\mathbb E\mathrm{Tr}\,\rho_A^2=\frac{a+b}{ab+1}.
```

*Demostración.* Sea $`G_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,g\rangle`$, una matriz de Ginibre; la
matriz de $`\psi`$ es $`\Psi=G/\|g\|`$ y, por la proposición 7.2(a),
$`(\rho_A)_{\alpha\alpha'}=\sum_\beta\Psi_{\alpha\beta}\bar\Psi_{\alpha'\beta}`$, es decir,
$`\rho_A=\Psi\Psi^*=GG^*/\|g\|^2`$. Por tanto
$`\mathrm{Tr}(GG^*)^2=\|g\|^4\mathrm{Tr}\,\rho_A^2`$, donde $`\mathrm{Tr}\,\rho_A^2\in[0,1]`$ es función de
$`\psi`$. Por el lema 7.9, $`\mathbb E\mathrm{Tr}(GG^*)^2=\mathbb E\|g\|^4\,\mathbb E\mathrm{Tr}\,\rho_A^2`$, y por el
lema 7.10, $`\mathbb E\mathrm{Tr}\,\rho_A^2=ab(a+b)/(ab(ab+1))`$. $`\square`$

El capítulo 6 obtiene la misma media a partir de los cuartos momentos de un vector unitario
de ley invariante unitaria (lema 6.14), como la pureza de salida de un canal aleatorio (teorema 6.15 y la
observación que le sigue).

<a id="r-7-12"></a>

**Corolario 7.12.** $`\log\frac{ab+1}{a+b}\le\mathbb E\,S(\rho_A)\le\log\min(a,b)`$.

*Demostración.* Por el corolario 3.6, $`S(\rho_A)\le\log a`$, y por el lema 7.7,
$`S(\rho_A)=S(\mathrm{Tr}_1|\psi\rangle\langle\psi|)\le\log b`$. La desigualdad de Jensen finita para la
función cóncava $`\log`$ (solución del ejercicio 2 del capítulo 3) con los pesos $`\lambda_i>0`$, los
valores propios no nulos de $`\rho_A`$, da $`\sum_i\lambda_i\log\lambda_i\le\log\sum_i\lambda_i^2`$, es decir,
$`S(\rho_A)\ge-\log X`$ con $`X=\mathrm{Tr}\,\rho_A^2\in(0,1]`$. Sea $`m=\mathbb EX=\frac{a+b}{ab+1}`$ (teorema 7.11).
El lema 3.4 para la función convexa $`-\log`$ sobre $`(0,\infty)`$ en el punto $`m`$ da
$`-\log X\ge-\log m-(X-m)/m`$; ambos miembros están acotados ($`X\ge1/a`$ por Cauchy–Schwarz), y tomando
esperanzas se obtiene $`\mathbb E[-\log X]\ge-\log m`$. $`\square`$

Para $`a`$ fijo y $`b\to\infty`$ las cotas se juntan: la entropía media de una parte pequeña de un
estado puro aleatorio tiende a su máximo $`\log a`$. La media exacta es la fórmula de Page,
demostrada en el apéndice C.

<a id="r-7-13"></a>

**Teorema 7.13 (fórmula de Page).** Para $`a\le b`$,
```math
\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}.
```
Además, $`\log a-\frac{a-1}{ab+1}-\frac{a-1}{2b}\le\mathbb E\,S(\rho_A)\le\log a-\frac{a-1}{2b}`$, de modo que
$`\big|\mathbb E\,S(\rho_A)-\log a+\frac a{2b}\big|\le\frac1{2b}`$; para dos cúbits, $`\mathbb E\,S(\rho_A)=\frac13`$.

*Demostración.* Es el teorema C.17 junto con el corolario C.18. La demostración del apéndice C parte
de $`\rho_A=GG^*/\mathrm{Tr}\,GG^*`$ (demostración del teorema 7.11), calcula la densidad conjunta de los
valores propios de la matriz de Wishart $`GG^*`$ (teorema C.13) y de $`\rho_A`$ (teorema C.14), e integra
$`x\log x`$ frente a la densidad de niveles de Laguerre (lemas C.15 y C.16). $`\square`$

*Cálculo:* los teoremas 7.11 y 7.13 y el corolario 7.12 en $`10^4`$ estados puros aleatorios de $`\mathbb C^3\otimes\mathbb C^6`$ — [tutorial, sección 7.5](../../tutorial/es/07_tensor_products.md#s-7-5).

La fórmula fue conjeturada por D. N. Page (Phys. Rev. Lett. 71, 1291, 1993) y demostrada por S. K.
Foong y S. Kanno (Phys. Rev. Lett. 72, 1148, 1994), J. Sánchez-Ruiz (Phys. Rev. E 52, 5653, 1995)
y S. Sen (Phys. Rev. Lett. 77, 1, 1996).

Para $`a>b`$ se intercambian los papeles de $`a`$ y $`b`$. En efecto, la unitaria
$`e_\alpha\otimes e_\beta\mapsto e_\beta\otimes e_\alpha`$ de $`\mathbb C^a\otimes\mathbb C^b`$ sobre $`\mathbb C^b\otimes\mathbb C^a`$
permuta las coordenadas de $`g`$, así que transforma un estado puro aleatorio de Haar $`\psi`$ en un
estado puro aleatorio de Haar $`\psi'`$ de $`\mathbb C^b\otimes\mathbb C^a`$; su matriz es $`\Psi^T`$, luego
$`\mathrm{Tr}_2|\psi'\rangle\langle\psi'|=\Psi^T\bar\Psi=\mathrm{Tr}_1|\psi\rangle\langle\psi|`$ (proposición 7.2(a)),
cuya entropía es $`S(\rho_A)`$ por el lema 7.7. Por tanto, con $`a'=\min(a,b)`$ y $`b'=\max(a,b)`$, el
teorema 7.13 dice que la entropía de entrelazamiento media se queda por debajo de su máximo
$`\log a'`$ en al menos $`\frac{a'-1}{2b'}`$ y a lo sumo $`\frac{a'-1}{2b'}+\frac{a'-1}{a'b'+1}`$, lo que es menos de
un nat: es $`0`$ para $`b'=1`$ y, en otro caso, $`\frac{a'-1}{2b'}\le\frac12`$ y $`\frac{a'-1}{a'b'+1}<\frac1{b'}\le\frac12`$.
Así, para $`ab`$ fijo, como función de $`a`$, la media está a menos de un nat de
$`\log\min(a,b)=\min(\log a,\log(ab)-\log a)`$.

La ley de $`\rho_A`$ es la **medida inducida** (el nombre es el de K. Życzkowski y H.-J. Sommers,
J. Phys. A 34, 7111, 2001; el teorema C.14 da la densidad de sus valores propios).
`states.random_density_matrix(alg, rank=r)` forma $`G_cG_c^*`$ con bloques de Ginibre independientes
$`G_c`$ de tamaño $`k_c\times r`$ y divide por $`\sum_c\mathrm{Tr}\,G_cG_c^*`$. Para un solo bloque, es la ley
de $`\rho_A`$ anterior con $`b=r`$ (demostración del teorema 7.11). Para varios bloques, es la ley de la
densidad de la restricción a $`\mathcal M_1`$ del estado vectorial de un estado puro aleatorio de Haar
$`\psi=g/\|g\|`$ de $`H_1\otimes\mathbb C^r`$. En efecto, sea $`G_{i\beta}=\langle e_i\otimes e_\beta,g\rangle`$, una
matriz de Ginibre $`D\times r`$ ($`D=\sum_ck_c`$), y sean $`G_c`$ sus filas del bloque $`c`$: matrices de
Ginibre $`k_c\times r`$ independientes. Tomemos $`\mathcal M_2=M_r`$ (un solo sector) y una regla de fusión
constante; entonces hay una única carga total, $`\Pi=1`$, y
$`\mathcal M_1\otimes_f\mathcal M_2=B(H_1\otimes\mathbb C^r)`$. Por el teorema 7.6, la restricción tiene la
densidad $`E_1(\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|)`$, y
$`\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|=GG^*/\|g\|^2`$ (proposición 7.2(a), como en la demostración del
teorema 7.11), cuyos bloques diagonales son $`G_cG_c^*/\|g\|^2`$, con $`\|g\|^2=\sum_c\mathrm{Tr}\,G_cG_c^*`$.

*En la biblioteca:* `states.random_density_matrix(alg, rank=r)` — [tutorial, sección 7.5](../../tutorial/es/07_tensor_products.md#s-7-5).

<a id="r-7-14"></a>

**Proposición 7.14.** Para $`\rho`$ generada por `random_density_matrix(alg, rank=r)` sobre
$`\bigoplus_cM_{k_c}`$, con $`D=\sum_ck_c`$,
```math
\mathbb E\mathrm{Tr}\,\rho^2=\frac{\sum_ck_c(k_c+r)}{D(Dr+1)} .
```

*Demostración.* $`\mathrm{Tr}\,\rho^2=\sum_c\mathrm{Tr}(G_cG_c^*)^2/\|g\|^4`$, con $`g\in\mathbb C^{Dr}`$ el conjunto
de todas las entradas, es una función de $`g/\|g\|`$ con valores en $`[0,1]`$. Como en el teorema 7.11,
el lema 7.9 da $`\sum_c\mathbb E\mathrm{Tr}(G_cG_c^*)^2=\mathbb E\|g\|^4\,\mathbb E\mathrm{Tr}\,\rho^2`$, y el lema 7.10 da
$`\mathbb E\mathrm{Tr}(G_cG_c^*)^2=k_cr(k_c+r)`$ y $`\mathbb E\|g\|^4=Dr(Dr+1)`$, así que
$`\mathbb E\mathrm{Tr}\,\rho^2=\sum_ck_cr(k_c+r)/(Dr(Dr+1))`$. $`\square`$

*Cálculo:* la proposición 7.14 para $`M_1\oplus M_2\oplus M_3`$, $`r=2`$ — [tutorial, sección 7.5](../../tutorial/es/07_tensor_products.md#s-7-5).

> **Significado físico** (interpretación física, no un teorema de este curso). La entropía de
> entrelazamiento media de un estado puro aleatorio, representada frente al tamaño de la parte con
> dimensión total fija (la curva descrita tras el teorema 7.13), se llama *curva de Page*. En física
> se usa para describir cómo sale la información de un agujero negro en la radiación de Hawking y
> por qué las partes de un sistema cerrado parecen térmicas; estas son interpretaciones, no
> resultados de este curso.

<a id="s-7-6"></a>

## 7.6 Con qué comparar

*En la biblioteca:* comparación con otro software — [tutorial, sección 7.6](../../tutorial/es/07_tensor_products.md#s-7-6).

<a id="x-7"></a>

## Ejercicios

1. Calcule el dual de $`\mathrm{Tr}_2`$ respecto de $`\tau_{\mathrm{vN}}`$ tanto sobre
   $`\mathcal M_1\otimes\mathcal M_2`$ como sobre $`\mathcal M_1`$.
2. (Subaditividad.) Sea $`\rho`$ un estado sobre $`\mathcal M_1\otimes\mathcal M_2`$ con estados reducidos
   invertibles $`\rho_1,\rho_2`$. Demuestre que $`S(\rho)\le S(\rho_1)+S(\rho_2)`$, con igualdad si y
   solo si $`\rho=\rho_1\otimes\rho_2`$.
3. Sea `site` el álgebra $`\mathbb C\oplus\mathbb C`$ con cargas $`0,1`$ (un sitio vacío u ocupado),
   fusionada consigo misma mediante la suma. Calcule $`E(\psi)`$ y sus dos partes para
   $`\psi=(|01\rangle+|10\rangle)/\sqrt2`$, una partícula compartida por dos sitios, y describa la
   restricción del estado vectorial de $`(|00\rangle+|11\rangle)/\sqrt2`$ al álgebra fusionada.

*Cálculo:* los tres ejercicios, numéricamente — [tutorial, Ejercicios](../../tutorial/es/07_tensor_products.md#x-7).

<details>
<summary>Soluciones</summary>

1. Los pesos son $`1/(C_1k_c)`$ sobre $`\mathcal M_1`$ y $`1/(C_1C_2k_cl_d)`$ sobre el producto. Para
   $`A\in\mathcal M_1`$, por la proposición 7.2(a),
   $`\tau_{\mathrm{vN}}((\mathrm{Tr}_2\,X)A)=\sum_c\frac1{C_1k_c}\sum_d\mathrm{Tr}(X_{cd}(A_c\otimes1))    =\sum_{c,d}\frac1{C_1C_2k_cl_d}\mathrm{Tr}\big(X_{cd}\,C_2l_d(A_c\otimes1)\big)`$, así que el dual es
   $`A\mapsto\bigoplus_{c,d}C_2l_d\,A_c\otimes1_{l_d}=A\otimes Z`$ con el elemento central
   $`Z=\bigoplus_dC_2l_d1_{l_d}`$ de $`\mathcal M_2`$. Es unital solo si $`Z=1`$, es decir, $`\mathcal M_2=\mathbb C`$,
   y no lo es en el ejemplo ($`Z=4\cdot1_2\oplus6\cdot1_3`$). Para `Tr_norm` el mismo cálculo da
   $`A\otimes\bigoplus_dl_d1_{l_d}`$, unital solo si todos los $`l_d=1`$. Las trazas `blunt` siempre hacen de
   $`\mathrm{Tr}_2`$ el dual de $`A\mapsto A\otimes1`$; en general lo hacen $`\mathrm{tr}_w`$ sobre
   $`\mathcal M_1`$ y los pesos $`w_{cd}=w_c`$ sobre el producto, por el mismo cálculo.

2. $`\rho_1\otimes\rho_2`$ es un estado (teorema 7.1(c)), invertible, con vectores propios
   $`u_i\otimes v_j`$ y valores propios $`\lambda_i\mu_j`$, luego
   $`\log(\rho_1\otimes\rho_2)=\log\rho_1\otimes1+1\otimes\log\rho_2`$ (teorema 1.5(a)). Por la
   definición de la traza parcial,
   $`\mathrm{Tr}\,\rho(\log\rho_1\otimes1)=\mathrm{Tr}\,\rho_1\log\rho_1`$, y análogamente para el segundo
   factor. Por tanto $`S(\rho\|\rho_1\otimes\rho_2)=S(\rho_1)+S(\rho_2)-S(\rho)`$, que es $`\ge0`$, con
   igualdad si y solo si $`\rho=\rho_1\otimes\rho_2`$ (corolario 3.6). Esta diferencia es la
   información mutua.

3. Los sectores son $`N=0,1,2`$, de tamaños $`1,2,1`$. $`\psi`$ está en $`N=1`$, con
   $`\psi_0=|01\rangle/\sqrt2`$, $`\psi_1=|10\rangle/\sqrt2`$, luego $`p=(\frac12,\frac12)`$ y, como los
   bloques son unidimensionales, $`E(\hat\psi_c)=0`$: por la proposición 7.8, $`E(\psi)=H(p)=\log2`$
   es todo entrelazamiento de número, y la parte configuracional es $`0`$. El vector
   $`(|00\rangle+|11\rangle)/\sqrt2`$ superpone $`N=0`$ y $`N=2`$; sobre el álgebra fusionada, diagonal
   por bloques, su estado vectorial es la mezcla
   $`\frac12|00\rangle\langle00|+\frac12|11\rangle\langle11|`$ de estados producto (capítulo 3,
   ejercicio 3).

</details>

**Siguiente:** el [capítulo 8](08_dynamics.md) trata la dinámica: las ecuaciones de Schrödinger y de
Lindblad en un álgebra con sectores, los saltos entre sectores y la postselección.
