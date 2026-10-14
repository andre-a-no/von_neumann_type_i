[English](../01_matrix_algebras.md) · **Español** · [Français](../fr/01_matrix_algebras.md) · [Deutsch](../de/01_matrix_algebras.md) · [中文](../zh/01_matrix_algebras.md) · [日本語](../ja/01_matrix_algebras.md) · [Русский](../ru/01_matrix_algebras.md)

# 1. Álgebras de matrices y cálculo funcional

Este capítulo fija la notación y demuestra los hechos sobre $`M_n(\mathbb C)`$ que se usan en
todos los capítulos posteriores: el teorema espectral, el cálculo funcional, la
descomposición polar y las normas definidas por los valores singulares. El código de la
biblioteca y las comprobaciones numéricas de estos resultados están en el
[Tutorial 1](../../tutorial/es/01_matrix_algebras.md).

**Convenios.** Los productos escalares $`\langle x,y\rangle=\sum_i\bar x_iy_i`$ son lineales en
el segundo argumento. $`A^*`$ es la traspuesta conjugada, de modo que $`\langle x,Ay\rangle=\langle A^*x,y\rangle`$.
$`e_1,\dots,e_n`$ es la base canónica, $`1`$ la matriz identidad y $`\mathrm{Tr}`$ la traza
usual.

<a id="s-1-1"></a>

## 1.1 El álgebra $`M_n(\mathbb C)`$

$`M_n(\mathbb C)`$ es el álgebra de las matrices complejas $`n\times n`$ con la involución
$`A\mapsto A^*`$ y la norma de operador $`\|A\|=\sup_{\|x\|=1}\|Ax\|`$. Un elemento $`A`$ es

| nombre | definición |
|---|---|
| normal | $`A^*A=AA^*`$ |
| autoadjunto | $`A=A^*`$ |
| positivo, $`A\ge0`$ | $`\langle x,Ax\rangle\ge0`$ para todo $`x`$ |
| unitario | $`A^*A=AA^*=1`$ |
| proyección | $`A=A^*=A^2`$ |

Escribimos $`A\le B`$ si $`B-A\ge0`$, y $`A>0`$ si $`A\ge0`$ y $`A`$ es invertible.

*En la biblioteca:* operadores con un espectro prescrito, autoadjuntos o positivos, `TypeIAlgebra.operator_from_eigenvalues` — [tutorial, sección 1.5](../../tutorial/es/01_matrix_algebras.md#s-1-5); la norma de operador, `operator_norm` — [tutorial, sección 1.6](../../tutorial/es/01_matrix_algebras.md#s-1-6).

<a id="r-1-1"></a>

**Lema 1.1 (identidad C\*).** $`\|A^*\|=\|A\|`$ y $`\|A^*A\|=\|A\|^2`$.

*Demostración.* $`\|Ax\|=\sup_{\|y\|=1}|\langle y,Ax\rangle|=\sup_{\|y\|=1}|\langle A^*y,x\rangle|`$,
luego $`\|A\|=\sup_{\|x\|=\|y\|=1}|\langle A^*y,x\rangle|=\|A^*\|`$. Además,
$`\|Ax\|^2=\langle x,A^*Ax\rangle\le\|A^*A\|`$ para $`\|x\|=1`$, de donde $`\|A\|^2\le\|A^*A\|\le\|A^*\|\,\|A\|=\|A\|^2`$. $`\square`$

> **Significado físico** (interpretación física, no un teorema de este curso). En mecánica
> cuántica, un sistema con $`n`$ configuraciones distinguibles se describe mediante $`\mathbb C^n`$.
> Las matrices autoadjuntas son los *observables* (energía, espín, posición en una red); sus
> valores propios son los resultados posibles de una medición. Las proyecciones son preguntas de respuesta sí/no.
> Las matrices unitarias son las simetrías y las evoluciones temporales, $`e^{-itH}`$ para un
> hamiltoniano $`H`$. Las matrices positivas de traza uno son los *estados* (matrices
> densidad, capítulo 3).

<a id="s-1-2"></a>

## 1.2 El teorema espectral

<a id="r-1-2"></a>

**Lema 1.2 (Schur).** Para toda $`A\in M_n(\mathbb C)`$ existe una unitaria $`U`$ tal que
$`U^*AU`$ es triangular superior.

*Demostración.* Inducción sobre $`n`$; el caso $`n=1`$ es trivial. El polinomio característico
de $`A`$ tiene una raíz $`\lambda`$ (teorema fundamental del álgebra), luego $`Av=\lambda v`$ para
algún vector unitario $`v`$. Completamos $`v`$ hasta una base ortonormal $`v,w_2,\dots,w_n`$ y
sea $`V`$ la unitaria cuyas columnas son estos vectores. Entonces
$`V^*AV=\begin{pmatrix}\lambda&*\\0&A'\end{pmatrix}`$ con $`A'\in M_{n-1}(\mathbb C)`$.
Por inducción, $`U'^*A'U'`$ es triangular superior para cierta unitaria $`U'`$, y
$`U=V\,(1\oplus U')`$ cumple lo pedido. $`\square`$

<a id="r-1-3"></a>

**Teorema 1.3 (teorema espectral).** $`A`$ es normal si y solo si
$`A=U\,\mathrm{diag}(\lambda_1,\dots,\lambda_n)\,U^*`$ para alguna unitaria $`U`$. Los
$`\lambda_i`$ son los valores propios de $`A`$ contados con multiplicidad.

*Demostración.* Si $`A=U\Lambda U^*`$ con $`\Lambda`$ diagonal, entonces
$`A^*A=U\bar\Lambda\Lambda U^*=U\Lambda\bar\Lambda U^*=AA^*`$. Recíprocamente, sea $`A`$
normal y $`T=U^*AU`$ triangular superior (lema 1.2); $`T`$ es normal. Comparando las
entradas $`(1,1)`$: $`(T^*T)_{11}=|t_{11}|^2`$ y $`(TT^*)_{11}=\sum_j|t_{1j}|^2`$, luego
$`t_{1j}=0`$ para $`j>1`$. Así $`T=t_{11}\oplus T'`$ con $`T'`$ normal y triangular superior, y
por inducción $`T`$ es diagonal. La diagonal de $`U^*AU`$ está formada por las raíces de
$`\det(A-\lambda)=\det(T-\lambda)`$. $`\square`$

*En la biblioteca:* el teorema 1.3 leído al revés, $`A=U\Lambda U^*`$ con valores propios prescritos y $`U`$ aleatoria según la medida de Haar, es `TypeIAlgebra.operator_from_eigenvalues` — [tutorial, sección 1.5](../../tutorial/es/01_matrix_algebras.md#s-1-5).

<a id="r-1-4"></a>

**Corolario 1.4.** Sea $`A`$ normal con valores propios $`\lambda_i`$. Entonces

1. $`\|A\|=\max_i|\lambda_i|`$;
2. $`A`$ es autoadjunta si y solo si todos los $`\lambda_i\in\mathbb R`$; unitaria si y solo si
   todos los $`|\lambda_i|=1`$; una proyección si y solo si todos los $`\lambda_i\in\{0,1\}`$;
3. $`A\ge0`$ si y solo si $`A`$ es autoadjunta y todos los $`\lambda_i\ge0`$. Además, toda matriz positiva
   es autoadjunta, luego normal, de modo que (3) caracteriza la positividad entre *todas* las matrices
   y el cálculo funcional de más abajo se aplica a toda $`A\ge0`$.

*Demostración.* Con $`y=U^*x`$, se tiene $`\|Ax\|^2=\sum_i|\lambda_i|^2|y_i|^2`$ y $`\|y\|=\|x\|`$,
lo que da (1), con igualdad para $`x=Ue_i`$ con $`|\lambda_i|`$ máximo. En (2), $`A^*=U\bar\Lambda U^*`$,
$`A^*A=U|\Lambda|^2U^*`$ y $`A^2=U\Lambda^2U^*`$, y las condiciones se cumplen si y solo si se
cumplen entrada a entrada para $`\Lambda`$. Para (3), en primer lugar: toda $`A`$ positiva es
autoadjunta, ya que
$`\langle x,(A-A^*)x\rangle=\langle x,Ax\rangle-\overline{\langle x,Ax\rangle}=0`$ para todo $`x`$,
y una matriz $`D`$ con $`\langle x,Dx\rangle=0`$ para todo $`x`$ es nula (polarización:
$`4\langle x,Dy\rangle=\sum_{k=0}^3 i^{-k}\langle x+i^ky,D(x+i^ky)\rangle`$). Entonces
$`\langle x,Ax\rangle=\sum_i\lambda_i|y_i|^2`$, que es $`\ge0`$ para todo $`x`$ si y solo si todos los $`\lambda_i\ge0`$. $`\square`$

<a id="s-1-3"></a>

## 1.3 El cálculo funcional

Sea $`A=U\Lambda U^*`$ normal y $`f`$ una función compleja definida en $`\mathrm{spec}\,A`$. Se define
```math
f(A)=U\,\mathrm{diag}\big(f(\lambda_1),\dots,f(\lambda_n)\big)\,U^*.
```

<a id="r-1-5"></a>

**Teorema 1.5.** (a) $`f(A)`$ no depende de la elección de $`U`$. (b) Para un polinomio $`p`$,
$`p(A)`$ es el polinomio usual en $`A`$. (c) $`f\mapsto f(A)`$ es lineal y multiplicativa,
$`\bar f(A)=f(A)^*`$, y $`\mathrm{spec}\,f(A)=f(\mathrm{spec}\,A)`$.
(d) Si $`B`$ conmuta con $`A`$, entonces $`B`$ conmuta con $`f(A)`$. (e) Si $`g`$ está definida en
$`f(\mathrm{spec}\,A)`$, entonces $`g(f(A))=(g\circ f)(A)`$.

*Demostración.* (b) $`A^k=U\Lambda^kU^*`$, luego $`p(A)=Up(\Lambda)U^*`$ por linealidad. (a) El
espectro es finito, así que la interpolación de Lagrange proporciona un polinomio $`p`$ con
$`p=f`$ en $`\mathrm{spec}\,A`$: si $`\mu_1,\dots,\mu_r`$ son los valores propios distintos, tómese
$`p(t)=\sum_if(\mu_i)\prod_{j\neq i}(t-\mu_j)/(\mu_i-\mu_j)`$. Entonces $`U\,f(\Lambda)\,U^*=U\,p(\Lambda)\,U^*=p(A)`$ para
*toda* unitaria $`U`$ que diagonalice $`A`$, y $`p(A)`$ está definido sin recurrir a $`U`$. (c) se
cumple para las matrices diagonales $`f(\Lambda)`$ y se conserva al conjugar por $`U`$; el
espectro de $`Uf(\Lambda)U^*`$ es la diagonal de $`f(\Lambda)`$. (d) $`B`$ conmuta con todo
polinomio en $`A`$, en particular con $`p(A)=f(A)`$. (e) $`f(A)=Uf(\Lambda)U^*`$ es normal y la
diagonaliza la misma $`U`$, luego por (a) $`g(f(A))=U\,g(f(\Lambda))\,U^*=(g\circ f)(A)`$. $`\square`$

<a id="r-1-6"></a>

**Teorema 1.6 (raíz cuadrada).** Para $`A\ge0`$ existe exactamente una $`B\ge0`$ con
$`B^2=A`$, a saber, $`B=\sqrt A`$ definida mediante el cálculo funcional.

*Demostración.* $`\sqrt A\ge0`$ y $`(\sqrt A)^2=A`$ por el teorema 1.5(c) y el corolario 1.4(3).
Sea $`B\ge0`$ con $`B^2=A`$, $`B=V\,\mathrm{diag}(\mu_j)\,V^*`$, $`\mu_j\ge0`$. Entonces
$`A=V\,\mathrm{diag}(\mu_j^2)\,V^*`$, luego $`\mathrm{spec}\,A=\{\mu_j^2\}`$. Sea $`p`$ un
polinomio con $`p(t)=\sqrt t`$ en $`\mathrm{spec}\,A`$. Entonces
$`p(A)=\sqrt A`$ y $`p(A)=p(B^2)=V\,\mathrm{diag}(p(\mu_j^2))\,V^*=V\,\mathrm{diag}(\mu_j)\,V^*=B`$. $`\square`$

*Cálculo:* los teoremas 1.5 y 1.6, comprobados con `apply_function`, `sqrt`, `log`, `expm` y `power` — [tutorial, sección 1.5](../../tutorial/es/01_matrix_algebras.md#s-1-5).

<a id="r-1-7"></a>

**Corolario 1.7.** $`A\ge0`$ si y solo si $`A=T^*T`$ para alguna $`T`$. (Si $`A\ge0`$, tómese
$`T=\sqrt A`$; recíprocamente, $`\langle x,T^*Tx\rangle=\|Tx\|^2\ge0`$.)

<a id="s-1-4"></a>

## 1.4 Descomposición polar y valores singulares

Para cualquier $`A`$, la matriz $`A^*A`$ es positiva (corolario 1.7); su raíz cuadrada
$`|A|=(A^*A)^{1/2}`$ es el **valor absoluto** de $`A`$.

<a id="r-1-8"></a>

**Teorema 1.8 (descomposición polar).** Toda $`A\in M_n(\mathbb C)`$ puede escribirse como
$`A=W|A|`$ con $`W`$ unitaria.

*Demostración.* Para todo $`x`$, $`\||A|x\|^2=\langle x,|A|^2x\rangle=\langle x,A^*Ax\rangle=\|Ax\|^2`$.
Por tanto $`W_0:|A|x\mapsto Ax`$ es una isometría lineal bien definida de la imagen
$`R_1`$ de $`|A|`$ sobre la imagen $`R_2`$ de $`A`$: si $`|A|x=|A|x'`$, entonces
$`\|A(x-x')\|=\||A|(x-x')\|=0`$. En particular $`\dim R_1=\dim R_2`$, de modo que los
complementos ortogonales también tienen la misma dimensión; sea $`W_1:R_1^\perp\to R_2^\perp`$
una isometría cualquiera entre ellos (que envíe una base ortonormal en una base
ortonormal). Sea $`W=W_0\oplus W_1`$, es decir, $`W(u_1+u_2)=W_0u_1+W_1u_2`$ para $`u_1\in R_1`$,
$`u_2\in R_1^\perp`$. Como $`W_0u_1\in R_2\perp R_2^\perp\ni W_1u_2`$, el teorema de Pitágoras da
$`\|Wu\|^2=\|u_1\|^2+\|u_2\|^2=\|u\|^2`$ para todo $`u`$. Por tanto $`\langle u,(W^*W-1)u\rangle=0`$ para
todo $`u`$, luego $`W^*W=1`$ por polarización (demostración del corolario 1.4(3)); entonces $`W`$ es
invertible con $`W^{-1}=W^*`$, así que también $`WW^*=1`$: $`W`$ es unitaria, y
$`W|A|x=W_0|A|x=Ax`$. $`\square`$

*Cálculo:* el teorema 1.8 sobre una matriz de Ginibre, con $`|A|`$ calculado por `abs` — [tutorial, sección 1.6](../../tutorial/es/01_matrix_algebras.md#s-1-6).

Los valores propios $`s_1\ge\dots\ge s_n\ge0`$ de $`|A|`$ son los **valores singulares** de
$`A`$. A lo largo del curso se usan tres normas:
```math
\|A\|=s_1,\qquad \|A\|_2=\Big(\sum_is_i^2\Big)^{1/2}=(\mathrm{Tr}\,A^*A)^{1/2},
\qquad\|A\|_1=\sum_is_i=\mathrm{Tr}|A|.
```
La primera identidad se cumple porque $`\|Ax\|=\||A|x\|`$ para todo $`x`$, luego
$`\|A\|=\||A|\|=s_1`$ por el corolario 1.4(1); la segunda, porque
$`\mathrm{Tr}\,A^*A=\mathrm{Tr}|A|^2=\sum_is_i^2`$.

<a id="r-1-9"></a>

**Teorema 1.9 (dualidad de la norma traza).**
$`\|A\|_1=\max_{U\text{ unitaria}}|\mathrm{Tr}(AU)|`$, y el máximo se alcanza en
$`U=W^*`$, donde $`A=W|A|`$. En particular $`|\mathrm{Tr}\,A|\le\|A\|_1`$.

*Demostración.* $`\mathrm{Tr}(AW^*)=\mathrm{Tr}(W|A|W^*)=\mathrm{Tr}|A|`$. Para cualquier unitaria $`U`$,
por la ciclicidad de la traza y la desigualdad de Cauchy–Schwarz para $`\langle X,Y\rangle_{HS}=\mathrm{Tr}(X^*Y)`$,
```math
|\mathrm{Tr}(AU)|=|\mathrm{Tr}(|A|^{1/2}\,UW\,|A|^{1/2})|
=|\langle|A|^{1/2},\,UW|A|^{1/2}\rangle_{HS}|
\le\||A|^{1/2}\|_2\,\|UW|A|^{1/2}\|_2=\mathrm{Tr}|A|,
```
ya que $`\|UW X\|_2=\|X\|_2`$ por ser $`UW`$ unitaria, y $`\||A|^{1/2}\|_2^2=\mathrm{Tr}|A|`$. $`\square`$

En particular, las tres cantidades anteriores son normas. $`\|\cdot\|`$ es una norma de
operador; $`\|\cdot\|_2`$ es la norma asociada al producto escalar $`\langle X,Y\rangle_{HS}`$; y
para $`\|\cdot\|_1`$, el teorema 1.9 da
$`\|A+B\|_1=\max_U|\mathrm{Tr}(AU)+\mathrm{Tr}(BU)|\le\|A\|_1+\|B\|_1`$ y
$`\|cA\|_1=\max_U|c\mathrm{Tr}(AU)|=|c|\,\|A\|_1`$, mientras que
$`\|A\|_1=0`$ obliga a que $`|A|=0`$ y $`A=W|A|=0`$. Para $`P\ge0`$, $`|P|=(P^2)^{1/2}=P`$ por el
teorema 1.6, luego $`\|P\|_1=\mathrm{Tr}\,P`$.

*Cálculo:* las tres normas (`operator_norm`, `frobenius_norm`, `trace_norm`) y el teorema 1.9 — [tutorial, sección 1.6](../../tutorial/es/01_matrix_algebras.md#s-1-6).

<a id="s-1-5"></a>

## 1.5 Operadores en la biblioteca

*En la biblioteca:* el álgebra $`M_n(\mathbb C)`$ como `TypeIAlgebra([n], [n])`, cómo se almacenan sus operadores y el cálculo funcional `apply_function` — [tutorial, sección 1.5](../../tutorial/es/01_matrix_algebras.md#s-1-5).

<a id="s-1-6"></a>

## 1.6 Normas y la dualidad de la norma traza

*Cálculo:* las normas, la descomposición polar y la dualidad de la sección 1.4, sobre una matriz de Ginibre — [tutorial, sección 1.6](../../tutorial/es/01_matrix_algebras.md#s-1-6).

<a id="s-1-7"></a>

## 1.7 Lotes: muchos operadores a la vez

Como ejemplo de una magnitud que se estima con muchos operadores aleatorios a la vez, calculamos un momento
exacto del **conjunto gaussiano unitario** (GUE).

Sea $`G`$ una matriz de entradas independientes, con partes real e imaginaria independientes
$`\sim N(0,\tfrac12)`$, de modo que $`\mathbb E|G_{ij}|^2=1`$, y sea $`H=(G+G^*)/\sqrt2`$.

<a id="r-1-10"></a>

**Proposición 1.10.** $`\mathbb E\,\mathrm{Tr}\,H^4=2n^3+n`$. Por tanto, para
$`X=H/\sqrt n`$, $`\ \mathbb E\,\tfrac1n\mathrm{Tr}\,X^4=2+\tfrac1{n^2}`$.

*Demostración.* *Fórmula de Isserlis (de Wick)* (teorema A.24, caso $`m=4`$). Las partes real e
imaginaria de las entradas de $`G`$ son variables gaussianas centradas independientes, así que son
las coordenadas de un único vector gaussiano centrado, y las entradas de $`H`$ y de $`\bar H`$ son
combinaciones lineales complejas de ellas. El teorema A.24 da entonces, para cuatro variables
cualesquiera de este tipo $`Z_1,\dots,Z_4`$ (se admiten repeticiones),
```math
\mathbb E[Z_1Z_2Z_3Z_4]=\mathbb E[Z_1Z_2]\,\mathbb E[Z_3Z_4]+\mathbb E[Z_1Z_3]\,\mathbb E[Z_2Z_4]+\mathbb E[Z_1Z_4]\,\mathbb E[Z_2Z_3],
```
la suma, sobre los tres emparejamientos de $`\{1,2,3,4\}`$, de los productos de las esperanzas por
pares (sin conjugación compleja). La demostración del apéndice A compara la serie de potencias de
la función generatriz de momentos con los momentos gaussianos y no deriva bajo la esperanza.

*El momento.* Las esperanzas por pares son
$`\mathbb E[H_{ij}H_{kl}]=\tfrac12\mathbb E[(G_{ij}+\bar G_{ji})(G_{kl}+\bar G_{lk})]=\delta_{il}\delta_{jk}`$,
porque $`\mathbb E[G_{ij}G_{kl}]=0`$ y $`\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}`$ (proposición A.16(a);
para $`(i,j)\neq(k,l)`$ las entradas son independientes y centradas, y para $`G_{ij}=x+iy`$ se tiene
$`\mathbb E[G_{ij}^2]=\mathbb E[x^2-y^2]+2i\,\mathbb E[xy]=0`$ y $`\mathbb E|G_{ij}|^2=\tfrac12+\tfrac12=1`$).
En $`\mathrm{Tr}\,H^4=\sum_{ijkl}H_{ij}H_{jk}H_{kl}H_{li}`$ los emparejamientos dan:
$`\{12\}\{34\}`$: $`\delta_{ik}\delta_{ki}`$, cuya suma es $`n^3`$ (índices libres $`i=k`$, $`j`$, $`l`$);
$`\{14\}\{23\}`$: $`\delta_{jl}\delta_{jl}`$, cuya suma es $`n^3`$;
$`\{13\}\{24\}`$: $`\delta_{il}\delta_{jk}\delta_{ji}\delta_{kl}`$, que obliga a $`i=j=k=l`$ y cuya suma es $`n`$.
Dividiendo entre $`n\cdot n^2`$ se obtiene la segunda afirmación. $`\square`$

Los dos emparejamientos *no cruzados* $`\{12\}\{34\}`$ y $`\{14\}\{23\}`$ (los pares pueden dibujarse
como arcos que no se cortan por encima de los puntos $`1,2,3,4`$) dan el término principal $`2`$, para
todo $`n`$; el emparejamiento cruzado $`\{13\}\{24\}`$ da la corrección $`1/n^2`$.

*Cálculo:* la proposición 1.10 por Monte Carlo, sobre un lote de $`200\,000`$ muestras — [tutorial, sección 1.7](../../tutorial/es/01_matrix_algebras.md#s-1-7).

<a id="s-1-8"></a>

## 1.8 Con qué comparar

*En la biblioteca:* implementaciones de referencia con las que comparar — [tutorial, sección 1.8](../../tutorial/es/01_matrix_algebras.md#s-1-8).

<a id="x-1"></a>

## Ejercicios

1. Demuestre que $`(fg)(A)=f(A)g(A)`$ para una $`A`$ normal y funciones $`f,g`$ sobre $`\mathrm{spec}\,A`$.
2. Demuestre que $`\|A\|_1=\|A^*\|_1`$ y que $`\|AB\|_1\le\|A\|\,\|B\|_1`$, usando el teorema 1.9.
3. La **media geométrica** de $`A,B>0`$ es $`A\#B=A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}`$.
   Demuestre que $`A\#B`$ es la única solución positiva $`X`$ de $`XA^{-1}X=B`$ y deduzca que
   $`A\#B=B\#A`$.

<details>
<summary>Soluciones</summary>

1. Por el teorema 1.5(c).

2. $`|\mathrm{Tr}(A^*U)|=|\overline{\mathrm{Tr}(U^*A)}|=|\mathrm{Tr}(AU^*)|`$, y
   $`U\mapsto U^*`$ es una biyección del grupo unitario, luego los máximos del teorema 1.9
   coinciden. Para la segunda afirmación, sea $`B=V|B|`$. Para $`U`$ unitaria,
   $`|\mathrm{Tr}(ABU)|=|\mathrm{Tr}(|B|^{1/2}UAV|B|^{1/2})|    \le\||B|^{1/2}\|_2\,\|UAV|B|^{1/2}\|_2\le\|A\|\mathrm{Tr}|B|`$, usando
   $`\|YX\|_2\le\|Y\|\,\|X\|_2`$ (aplíquese $`\|Yx\|\le\|Y\|\,\|x\|`$ a cada columna de $`X`$) y
   $`\|UAV\|=\|A\|`$ para unitarias $`U,V`$ (conservan la longitud de los vectores).

3. Escribamos $`T=A^{-1/2}XA^{-1/2}`$. Entonces $`XA^{-1}X=B`$ si y solo si $`T^2=A^{-1/2}BA^{-1/2}`$, y
   $`X\ge0`$ si y solo si $`T\ge0`$, ya que $`\langle x,Xx\rangle=\langle A^{1/2}x,TA^{1/2}x\rangle`$ y
   $`A^{1/2}`$ es invertible (aquí $`A^{\pm1/2}`$ vienen dadas por el cálculo funcional con $`t^{\pm1/2}`$
   en $`\mathrm{spec}\,A\subset(0,\infty)`$, luego $`A^{-1/2}A^{1/2}=1`$ por el teorema 1.5(c)).
   Además, $`A^{-1/2}BA^{-1/2}\ge0`$, ya que
   $`\langle x,A^{-1/2}BA^{-1/2}x\rangle=\langle A^{-1/2}x,BA^{-1/2}x\rangle\ge0`$.
   Por el teorema 1.6, $`T=(A^{-1/2}BA^{-1/2})^{1/2}`$ es la única
   solución positiva, es decir, $`X=A\#B`$. Toda solución $`X`$ de $`XA^{-1}X=B`$ es invertible,
   porque $`\det(X)^2\det(A)^{-1}=\det B\neq0`$. Así, $`XA^{-1}X=B`$ equivale a
   $`XB^{-1}X=A`$ (de $`A^{-1}=X^{-1}BX^{-1}`$, invirtiendo ambos miembros), y la unicidad da
   $`A\#B=B\#A`$.

</details>

*Cálculo:* los ejercicios 1 y 3, numéricamente — [tutorial, Ejercicios](../../tutorial/es/01_matrix_algebras.md#x-1).

**Siguiente:** el [capítulo 2](02_von_neumann_algebras.md) sustituye $`M_n`$ por una suma
directa de álgebras de matrices y demuestra que toda álgebra de von Neumann de dimensión
finita es de esta forma.
