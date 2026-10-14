[English](../00_overview.md) · **Español** · [Français](../fr/00_overview.md) · [Deutsch](../de/00_overview.md) · [中文](../zh/00_overview.md) · [日本語](../ja/00_overview.md) · [Русский](../ru/00_overview.md)

# 0. De qué trata este curso

Las álgebras de von Neumann se estudian sobre todo con lápiz y papel. En dimensión finita, cada
una de ellas es $*$-isomorfa a una suma directa de álgebras de matrices completas (teorema 2.6),
de modo que en principio todo se puede calcular; pero cuando un álgebra así se trata
numéricamente, a menudo se aplana en una única matriz grande, y su estructura (los sectores, el
centro, las múltiples trazas) desaparece del cálculo. Este curso y la biblioteca
`torch_vn_algebra` adoptan el enfoque contrario: **el álgebra misma es el tipo de datos**, y los
cálculos en ella se realizan de forma eficiente.

## 0.1 La observación

Un álgebra de von Neumann $\mathcal M$ sobre $\mathbb C^n$ es, salvo un cambio de base unitario,
$$
\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c}
$$
(teorema 2.6). Como álgebra abstracta es $\bigoplus_cM_{k_c}(\mathbb C)$: un elemento es una lista
de $C$ bloques cuadrados. Esta es exactamente la forma de un tensor con un *eje de sectores*: un
elemento de $\mathcal M$ es un arreglo de forma $(C,k,k)$ (con los bloques rellenados con ceros
hasta $k=\max_ck_c$; así lo almacena la biblioteca, sección 2.5), y $B$ elementos (muestras de
Monte Carlo, puntos de partida de una optimización, pasos de tiempo, realizaciones del desorden)
forman un arreglo de forma $(B,C,k,k)$. Las bibliotecas modernas de arreglos (PyTorch en CPU y
GPU) ofrecen rutinas por lotes (productos de matrices, descomposiciones espectrales,
descomposiciones QR) precisamente para estos «lotes de matrices pequeñas», que también usan las
redes neuronales.

Bajo esta identificación, las operaciones del álgebra son operaciones por lotes sobre arreglos.
La correspondencia es exacta (cada una de las seis primeras filas de la tabla es una identidad,
demostrada más abajo a partir de resultados del curso); solo las rutinas numéricas que ejecutan
las operaciones trabajan en coma flotante:

| en el álgebra $\mathcal M$ | sobre el arreglo de forma $(B,C,k,k)$ |
|---|---|
| producto $AB$, suma, adjunto | producto matricial por lotes, suma, traspuesta conjugada sobre $(B,C)$ |
| cálculo funcional $f(A)$ (capítulo 1) | descomposición espectral por lotes de los bloques |
| una traza $\sum_cw_c\operatorname{Tr}A_c$ (capítulo 3) | suma de la diagonal de cada bloque, suma ponderada sobre el eje de sectores |
| el centro, las proyecciones centrales (capítulo 2) | un número por sector |
| esperanza condicional sobre el centro | traza de cada bloque dividida entre $k_c$, propagada de nuevo al bloque |
| medida de Haar en el grupo unitario $\prod_cU(k_c)$ | descomposición QR por lotes de bloques gaussianos, con la corrección de fase del teorema 5.3 |
| estados, canales, dinámica (capítulos 3, 6, 8) | construidos a partir de las operaciones anteriores |

*Demostraciones de las seis primeras filas.* Por el teorema 2.6 existe un unitario $V$ con
$\mathcal M=V\big(\bigoplus_cM_{k_c}\otimes1_{m_c}\big)V^*$, y
$\varphi\big(V(\bigoplus_cA_c\otimes1_{m_c})V^*\big)=\bigoplus_cA_c$ es un $*$-isomorfismo de $\mathcal M$ sobre
$\bigoplus_cM_{k_c}(\mathbb C)$ (bien definido, pues $A_c\otimes1_{m_c}$ determina $A_c$; multiplicativo, pues
$(A_c\otimes1)(B_c\otimes1)=A_cB_c\otimes1$). La biblioteca almacena $\varphi(A)$ (sección 2.5).
(1) En $\bigoplus_cM_{k_c}$ las sumas, los productos y los adjuntos se toman bloque a bloque.
(2) Sean $A\in\mathcal M$ normal, $\varphi(A)=\bigoplus_cA_c$ y $f$ una función sobre $\operatorname{spec}A$. Cada $A_c$ es
normal, ya que $\varphi(A)^*\varphi(A)=\varphi(A^*A)=\varphi(AA^*)=\varphi(A)\varphi(A)^*$ se cumple bloque a bloque. Si
$A_cx=\lambda x$ con $x\neq0$, entonces $V(0\oplus\dots\oplus x\otimes e_1\oplus\dots\oplus0)$ ($x\otimes e_1$ en el
$c$-ésimo sumando) es un vector propio de $A$ para $\lambda$, de modo que $\operatorname{spec}A_c\subseteq\operatorname{spec}A$. Sea $p$ un
polinomio con $p=f$ sobre $\operatorname{spec}A$ (interpolación de Lagrange, como en la demostración del teorema 1.5). Para todo
$X$ normal con $\operatorname{spec}X\subseteq\operatorname{spec}A$ y $X=U\Lambda U^*$, se tiene $f(X)=Uf(\Lambda)U^*=Up(\Lambda)U^*=p(X)$ por
el teorema 1.5(b); así $f(A)=p(A)$, $f(A_c)=p(A_c)$, y como $\varphi$ es un homomorfismo unital,
$$
\varphi(f(A))=p(\varphi(A))=\bigoplus_cp(A_c)=\bigoplus_cf(A_c),\qquad f(A_c)=U_cf(\Lambda_c)U_c^*
$$
para cualquier diagonalización $A_c=U_c\Lambda_cU_c^*$ (teorema 1.3): una descomposición espectral por bloque.
(3) Las trazas de $\bigoplus_cM_{k_c}$ son los funcionales $\sum_cw_c\operatorname{Tr}A_c$ (teorema 3.1).
(4) $Z(\mathcal M)$ está generado por las proyecciones centrales minimales $P_c$, y $\varphi(P_c)$ es $1_{k_c}$ en el bloque
$c$ y $0$ en los demás (teorema 2.6), de modo que $\varphi(Z(\mathcal M))=\{\bigoplus_cz_c1_{k_c}:z\in\mathbb C^C\}$.
(5) Proposición 2.8. (6) Corolario 5.4 (bloques independientes, cada uno con la medida de Haar en $U(k_c)$) y teorema 5.3
(cada bloque a partir de la descomposición QR de una matriz gaussiana). La última fila es solo una indicación: los estados
se representan mediante densidades (teorema 3.2), los canales entre sectores mediante operadores de Kraus entre bloques
(proposición 6.4), y la dinámica cerrada mediante $e^{-itH}=\bigoplus_ce^{-itH_c}$ (proposición 8.2).

Dos consecuencias hacen que esto sea algo más que una cuestión de contabilidad.

- **Coste.** Una descomposición espectral en $\mathcal M$ se calcula bloque a bloque (fila 2 de la
  tabla). Los solucionadores densos de valores propios estándar (LAPACK, al que PyTorch llama en
  una CPU) usan un número de operaciones en coma flotante que crece como el cubo del tamaño de la
  matriz; esto es una descripción de esos algoritmos, no un resultado de este curso. En el modelo
  en el que un problema de valores propios $k\times k$ cuesta una constante por $k^3$, una
  descomposición espectral en $\mathcal M$ cuesta $\sum_ck_c^3$ frente a $D^3$ para la inmersión en
  $M_D(\mathbb C)$, $D=\sum_ck_c$; con $C$ bloques iguales de tamaño $k$, el cociente es
  $(Ck)^3/(Ck^3)=C^2$. Los productos se ejecutan sobre el tensor rellenado; el producto escolar de
  dos matrices $k\times k$ usa $k^3$ multiplicaciones ($k^2$ entradas, cada una una suma de $k$
  productos), de modo que un producto en $\mathcal M$ usa $C\,k_{\max}^3$, que coincide con
  $\sum_ck_c^3$ cuando todos los bloques tienen el mismo tamaño. Un elemento se almacena como
  $C\,k_{\max}^2$ números frente a $D^2$; con bloques iguales, el cociente es $(Ck)^2/(Ck^2)=C$. Lo
  que hace el código: los productos, las sumas y las trazas son operaciones tensoriales únicas
  sobre todo el arreglo $(B,C,k_{\max},k_{\max})$; las descomposiciones espectrales hermíticas, los
  valores singulares y las descomposiciones QR se aplican a la pila de los $B\cdot n_k$ bloques de
  cada tamaño distinto $k$ ($n_k$ es el número de sectores de ese tamaño) en una sola llamada por
  lotes (una única llamada cuando todos los sectores tienen el mismo tamaño). En una CPU con
  $t\ge2$ hilos de PyTorch, la biblioteca divide una llamada así en $t$ trozos que se ejecutan en
  hilos paralelos, siempre que la pila tenga al menos $2t$ bloques y (número de bloques)$\,\cdot k^3\ge2^{18}$
  (módulo `cost`, función `_run`); los resultados son las mismas llamadas a LAPACK sobre las mismas
  matrices. El objetivo es mantener ocupada una GPU, o todos los núcleos de una CPU.
- **Corrección por construcción.** Todo arreglo de forma $(C,k_{\max},k_{\max})$ con relleno de
  ceros es un elemento de $\mathcal M$, así que todo lo que se calcula a partir de los bloques
  permanece en $\mathcal M$. Los unitarios aleatorios de la biblioteca se generan mediante la
  construcción del teorema 5.3 de forma independiente en cada bloque (`random_unitary_operator`
  llama al muestreador por bloque `random_unitary` una vez por tamaño de bloque, y la sección 5.1
  comprueba `random_unitary` frente a esta construcción), lo que da la medida de Haar de
  $\prod_cU(k_c)$ por el corolario 5.4 (para números aleatorios gaussianos ideales); un unitario de
  Haar de $U(D)$ pertenece a este subgrupo con probabilidad cero cuando $C\ge2$ (sección 5.1, tras
  el corolario 5.4). Un paso de gradiente modifica los bloques y por tanto no puede salirse del
  álgebra. Y las trazas son los funcionales $\sum_cw_c\operatorname{Tr}A_c$ del teorema 3.1 con
  pesos elegidos explícitamente (y las densidades y las esperanzas condicionales se toman respecto
  de ellas), no la traza de algún $\mathbb C^n$ sobre el que $\mathcal M$ resulte actuar, que
  fijaría $w_c=m_c$ (sección 2.3).

## 0.2 Un primer cálculo

Comparamos el mismo cálculo realizado en $\mathcal M=\bigoplus_{c=1}^{64}M_{16}(\mathbb C)$ y en su
inmersión en $M_{1024}(\mathbb C)$: la raíz cuadrada de un elemento positivo aleatorio, y su traza.
`TypeIAlgebra(n_factors, k_factors)` recibe dos listas de tamaños de bloque: los tamaños nominales
$n_c$ y los tamaños $k_c\le n_c$ de los bloques que realmente se usan (sección 1.5). Toda álgebra
que el curso construye con `TypeIAlgebra` tiene dos listas iguales.

```python
import time
import torch
from torch_vn_algebra import TypeIAlgebra

torch.manual_seed(0)
C, k = 64, 16
alg = TypeIAlgebra([k] * C, [k] * C, precision='double')       # M = (+)_{c=1}^{64} M_16(C)

A = alg.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(1, d), force_positive=True)
blocks = A.matrix                                                # build A before timing

t0 = time.perf_counter()
root_blocks = A.sqrt().matrix                                    # (1, 64, 16, 16): 64 small eigenproblems
t_blocks = time.perf_counter() - t0

big = torch.block_diag(*blocks[0])                               # the same element in M_1024(C)
t0 = time.perf_counter()
w, V = torch.linalg.eigh(big)                                    # one large eigenproblem
root_big = (V * w.sqrt()) @ V.conj().T
t_big = time.perf_counter() - t0

assert torch.allclose(torch.block_diag(*root_blocks[0]), root_big, atol=1e-10)   # row 2: f(A) = (+)_c f(A_c)
assert abs(A.sqrt().Tr_blunt().item() - root_big.trace().item()) < 1e-9
assert big.numel() == C * blocks.numel()                          # storage: a factor C = 64
print(f"in the algebra: {t_blocks * 1e3:.1f} ms, as one 1024 x 1024 matrix: {t_big * 1e3:.1f} ms")
```

Los resultados coinciden salvo errores de redondeo (fila 2 de la tabla), y el cociente de
almacenamiento es $C$, como se contó arriba. Los 64 bloques van al solucionador de valores propios
en una sola llamada por lotes sobre una pila de forma $(64,16,16)$; con $2\le t\le32$ hilos de
PyTorch la biblioteca la divide en $t$ trozos ($64\cdot16^3=2^{18}$ alcanza exactamente el umbral
del módulo `cost`), con un solo hilo no lo hace. *Observación numérica* (los tiempos dependen del
hardware; el código los imprime): en dos CPU de 4 núcleos, la versión en el álgebra tardó 3–10 ms
y la versión con la matriz grande 400–760 ms, un factor entre varias decenas y unos pocos cientos
(en la segunda máquina tanto con 4 hilos como con uno). El factor medido queda muy por debajo del
cociente $C^2=4096$ del modelo de coste cúbico: un recuento de operaciones es un modelo del tiempo
de ejecución, no el tiempo de ejecución, y en estas máquinas los problemas de valores propios
pequeños se ejecutaron a un ritmo de operaciones por segundo mucho menor que el grande. En el
modelo de coste, el cociente $C^2$ crece con el número de sectores. Con un lote de $B$ elementos,
la versión en el álgebra sigue siendo una sola llamada por lotes; la versión con la matriz grande
tiene que resolver $B$ problemas de valores propios de tamaño $1024$.

## 0.3 Qué existe y qué es nuevo

Hemos buscado software que calcule numéricamente *dentro* de álgebras de von Neumann de dimensión
finita. Lo que encontramos se divide en cuatro grupos, todos útiles y todos con fines distintos:

- **Software matricial general** (NumPy/SciPy, QuTiP y otras herramientas cuánticas) trabaja con
  una matriz cada vez. La estructura de bloques, si la hay, debe mantenerse a mano.
- **Bibliotecas de redes tensoriales simétricas** almacenan los tensores como arreglos dispersos
  por bloques sobre los sectores de carga de una simetría (abeliana en YASTN [1], TeNPy e ITensor;
  TensorKit también trata las no abelianas); YASTN tiene un backend de PyTorch para GPU y
  diferenciación automática. Técnicamente son los parientes más cercanos. Su objeto es una red
  tensorial que hay que contraer, no un álgebra de operadores: hasta donde hemos visto en su
  documentación, no proporcionan la familia de trazas, estados y densidades, el centro y las
  esperanzas condicionales, las medidas de Haar en $\prod_cU(k_c)$ ni los canales entre sectores.
- **La diagonalización numérica por bloques** de $*$-álgebras de matrices [2, 3, 4] encuentra la
  descomposición $\mathcal M\cong\bigoplus_cM_{k_c}\otimes1_{m_c}$ de un álgebra dada, normalmente
  para reducir el tamaño de programas semidefinidos con simetría. Esto es complementario: esos
  algoritmos encuentran los bloques; esta biblioteca calcula en ellos. Su paso de entrada (módulo
  `decompose`, sección 4.3) usa la idea de [2, 3] de que los subespacios propios de un elemento
  aleatorio del conmutante encuentran los bloques; la sección 4.3 demuestra que en aritmética
  exacta esto tiene éxito con probabilidad uno (lema 4.6, proposición 4.7, lema 4.8), y la
  proposición 4.9 es el criterio que `Decomposition.verify()` comprueba numéricamente. En álgebra
  computacional existen descomposiciones simbólicas exactas de álgebras de grupo (por ejemplo, el
  paquete Wedderga de GAP).
- **La optimización polinómica no conmutativa** (NCSOStools, ncpol2sdpa) acota valores propios y
  trazas de polinomios en variables no conmutativas, para todos los tamaños de matriz a la vez,
  mediante relajaciones semidefinidas [5]. Responde a otra pregunta (una cota válida para todos
  los tamaños de matriz) con otros medios (relajación convexa).

Este panorama es una descripción del software tal como lo encontramos, no un teorema. Hasta donde
sabemos, no existe ninguna biblioteca general en la que un álgebra de von Neumann de dimensión
finita sea el tipo de datos básico para el cálculo numérico por lotes y diferenciable. La novedad
no es, por tanto, matemática nueva, sino **cálculo eficiente**: la identificación de la §0.1,
implementada de modo que se preserve la estructura del álgebra y las operaciones se ejecuten como
núcleos por lotes (una llamada por cada tamaño de bloque distinto en las descomposiciones, dividida
en un trozo por hilo en una CPU de varios núcleos por encima del umbral de tamaño de la §0.1). Si
conoce software relacionado, nos alegrará saberlo.

## 0.4 Cómo leer este curso

| capítulos | contenido |
|---|---|
| 1–3 | las matemáticas: álgebras de matrices, el teorema de estructura, trazas, estados y entropía, con demostraciones completas |
| 4–8 | de dónde surgen estas álgebras y qué se calcula en ellas: simetrías, matrices aleatorias, canales, productos tensoriales, dinámica |
| 9–12 | aplicaciones: cadenas de espines, computación cuántica, optimización, aprendizaje automático equivariante |
| 13 | fundamentos: álgebras de Jordan y por qué la mecánica cuántica usa *-álgebras |
| apéndices | demostraciones completas de los resultados auxiliares que los capítulos citan por su número (por ejemplo, el teorema A.24) |

Cada uno de los capítulos 1–13 enuncia resultados con sus demostraciones, muestra el código de la
biblioteca, lo contrasta numéricamente con resultados demostrados en el curso e indica las
herramientas estándar con las que compararlo. Los requisitos previos (álgebra lineal, análisis
real, teoría de la medida, probabilidad elemental) se enumeran en el
[README](README.md#requisitos-previos); todo lo demás se demuestra en el curso, y los resultados de
otros capítulos y de los apéndices se citan por su número. Los números que el curso no demuestra
(por ejemplo, los tiempos de ejecución) aparecen solo como *observaciones numéricas* del propio
código del capítulo, y los recuadros titulados *Significado físico* son interpretación física, no
teoremas de este curso.

## Referencias

1. M. M. Rams, G. Wójtowicz, A. Sinha, J. Hasik, *YASTN: Yet another symmetric tensor networks; A
   Python library for Abelian symmetric tensor network calculations*, SciPost Phys. Codebases 52
   (2025).
2. K. Murota, Y. Kanno, M. Kojima, S. Kojima, *A numerical algorithm for block-diagonal
   decomposition of matrix $*$-algebras with application to semidefinite programming*, Japan J.
   Indust. Appl. Math. 27 (2010) 125–160 (preprint: Technical Report METR 2007-52, University of
   Tokyo, 2007).
3. T. Maehara, K. Murota, *A numerical algorithm for block-diagonal decomposition of matrix
   $*$-algebras with general irreducible components*, Japan J. Indust. Appl. Math. 27 (2010)
   263–293.
4. E. de Klerk, C. Dobre, D. V. Pasechnik, *Numerical block diagonalization of matrix
   $*$-algebras with application to semidefinite programming*, Math. Program. 129 (2011) 91–111.
5. I. Klep, J. Povh, *Semidefinite programming and sums of hermitian squares of noncommutative
   polynomials*, J. Pure Appl. Algebra 214 (2010) 740–749.

**Siguiente:** [capítulo 1](01_matrix_algebras.md).
