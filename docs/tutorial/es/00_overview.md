[English](../00_overview.md) · **Español** · [Français](../fr/00_overview.md) · [Deutsch](../de/00_overview.md) · [中文](../zh/00_overview.md) · [日本語](../ja/00_overview.md) · [Русский](../ru/00_overview.md)

# 0. De qué trata este curso

Las álgebras de von Neumann se estudian casi exclusivamente con lápiz y papel. En dimensión
finita, cada una de ellas es una suma directa de álgebras de matrices completas (capítulo 2), de
modo que en principio todo se puede calcular; pero en la práctica el álgebra suele aplanarse en
una única matriz grande, y su estructura (los sectores, el centro, las múltiples trazas)
desaparece del cálculo. Este curso y la biblioteca `torch_vn_algebra` adoptan el enfoque
contrario: **el álgebra misma es el tipo de datos**, y los cálculos en ella se realizan de forma
eficiente.

## 0.1 La observación

Un álgebra de von Neumann de dimensión finita es, salvo un cambio de base unitario,
$$
\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c}
$$
(teorema 2.6). Como álgebra abstracta es $\bigoplus_cM_{k_c}(\mathbb C)$: un elemento es una lista
de $C$ bloques cuadrados. Esta es exactamente la forma de un tensor con un *eje de sectores*: un
elemento de $\mathcal M$ es un arreglo de forma $(C,k,k)$ (con los bloques rellenados hasta
$k=\max_ck_c$), y $B$ elementos (muestras de Monte Carlo, puntos de partida de una optimización,
pasos de tiempo, realizaciones del desorden) forman un arreglo de forma $(B,C,k,k)$. Las
bibliotecas modernas de arreglos (PyTorch en CPU y GPU) están diseñadas precisamente para estos
«lotes de matrices pequeñas», porque las redes neuronales los necesitan.

Bajo esta identificación, las operaciones del álgebra son operaciones por lotes sobre arreglos,
de manera exacta y sin aproximación alguna:

| en el álgebra $\mathcal M$ | sobre el arreglo de forma $(B,C,k,k)$ |
|---|---|
| producto $AB$, suma, adjunto | producto matricial por lotes, suma, traspuesta conjugada sobre $(B,C)$ |
| cálculo funcional $f(A)$ (capítulo 1) | descomposición espectral por lotes de los bloques |
| una traza $\sum_cw_c\operatorname{Tr}A_c$ (capítulo 3) | suma de la diagonal de cada bloque, suma ponderada sobre el eje de sectores |
| el centro, las proyecciones centrales (capítulo 2) | un número por sector |
| esperanza condicional sobre el centro | traza de cada bloque dividida entre $k_c$, propagada de nuevo al bloque |
| medida de Haar en el grupo unitario $\prod_cU(k_c)$ | descomposición QR por lotes de bloques gaussianos (capítulo 5) |
| estados, canales, dinámica (capítulos 3, 6, 8) | las mismas operaciones, compuestas |

Dos consecuencias hacen que esto sea algo más que una cuestión de contabilidad.

- **Coste.** Una descomposición espectral en $\mathcal M$ cuesta $\sum_ck_c^3$ operaciones (los
  bloques se agrupan por tamaño), frente a $D^3$ para la inmersión en $M_D(\mathbb C)$,
  $D=\sum_ck_c$; con $C$ bloques iguales, esto supone un factor $C^2$. Los productos y las sumas
  se ejecutan sobre el tensor rellenado y cuestan $C\,k_{\max}^3$, que coincide con $\sum_ck_c^3$
  cuando todos los bloques tienen el mismo tamaño. Los productos, las sumas y las trazas procesan
  los $B\cdot C$ bloques en una sola llamada; las descomposiciones (valores propios, valores
  singulares) se ejecutan como una llamada por lotes por cada tamaño de bloque distinto (una sola
  llamada cuando todos los sectores tienen el mismo tamaño), y en una CPU de varios núcleos la
  biblioteca divide cada una de esas llamadas en un trozo por hilo. Esto mantiene ocupada una
  GPU, o todos los núcleos de una CPU.
- **Corrección por construcción.** Todo lo que se calcula permanece en $\mathcal M$: los unitarios
  aleatorios están distribuidos según la medida de Haar en $\prod_cU(k_c)$ (un conjunto de medida
  nula dentro de $U(D)$, que nunca se alcanza muestreando $U(D)$), los pasos de gradiente no se
  salen del álgebra, y la traza, las densidades y las esperanzas condicionales se refieren al
  álgebra y no a una inmersión accidental.

## 0.2 Un primer cálculo

Comparamos el mismo cálculo realizado en $\mathcal M=\bigoplus_{c=1}^{64}M_{16}(\mathbb C)$ y en su
inmersión en $M_{1024}(\mathbb C)$: la raíz cuadrada de un elemento positivo aleatorio, y su traza.
`TypeIAlgebra(n_factors, k_factors)` recibe dos listas de tamaños de bloque: los tamaños nominales
$n_c$ y los tamaños $k_c\le n_c$ de los bloques que realmente se usan; en este curso ambas listas
son siempre iguales (sección 1.5).

```python
import time
import torch
from torch_vn_algebra import TypeIAlgebra

torch.manual_seed(0)
C, k = 64, 16
alg = TypeIAlgebra([k] * C, [k] * C, precision='double')       # M = (+)_{c=1}^{64} M_16(C)

A = alg.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(1, d), force_positive=True)

t0 = time.perf_counter()
root_blocks = A.sqrt().matrix                                    # (1, 64, 16, 16): 64 small eigenproblems
t_blocks = time.perf_counter() - t0

big = torch.block_diag(*A.matrix[0])                             # the same element in M_1024(C)
t0 = time.perf_counter()
w, V = torch.linalg.eigh(big)                                    # one large eigenproblem
root_big = (V * w.sqrt()) @ V.conj().T
t_big = time.perf_counter() - t0

assert torch.allclose(torch.block_diag(*root_blocks[0]), root_big, atol=1e-10)   # same result
assert abs(A.sqrt().Tr_blunt().item() - root_big.trace().item()) < 1e-9
print(f"in the algebra: {t_blocks * 1e3:.1f} ms, as one 1024 x 1024 matrix: {t_big * 1e3:.1f} ms")
```

Los resultados coinciden salvo errores de redondeo. La versión en el álgebra es decenas de veces
más rápida en una CPU (unos 10 ms frente a 400 ms en una máquina de 4 núcleos; las cifras exactas
dependen del hardware, y el código las imprime). Los 64 bloques van al solucionador de valores
propios en una sola llamada por lotes, que en una CPU de varios núcleos la biblioteca divide en un
trozo por hilo. Aun así, la aceleración es mucho menor de lo que sugiere el recuento de
operaciones $C^2=4096$, porque una sola llamada a un solucionador de valores propios grande
aprovecha mejor el hardware que muchas pequeñas. La diferencia crece con el número de sectores,
y el ahorro de memoria es un factor $C$. Con un lote de $B$ elementos, la versión en el álgebra
sigue siendo una sola llamada por lotes; la versión con la matriz grande se convierte en $B$
problemas de valores propios grandes.

## 0.3 Qué existe y qué es nuevo

Hemos buscado software que calcule numéricamente *dentro* de álgebras de von Neumann de dimensión
finita. Lo que encontramos se divide en cuatro grupos, todos útiles y todos con fines distintos:

- **Software matricial general** (NumPy/SciPy, QuTiP y otras herramientas cuánticas) trabaja con
  una matriz cada vez. La estructura de bloques, si la hay, debe mantenerse a mano.
- **Bibliotecas de redes tensoriales simétricas** almacenan los tensores como arreglos dispersos
  por bloques sobre los sectores de carga de una simetría abeliana, algunas con un backend de
  PyTorch para GPU y diferenciación automática (por ejemplo, YASTN [1]; también TeNPy, ITensor,
  TensorKit). Técnicamente son los parientes más cercanos. Su objeto es una red tensorial que hay
  que contraer, no un álgebra de operadores: no proporcionan la familia de trazas, estados y
  densidades, el centro y las esperanzas condicionales, las medidas de Haar en $\prod_cU(k_c)$ ni
  los canales entre sectores.
- **La diagonalización numérica por bloques** de $*$-álgebras de matrices [2, 3, 4] encuentra la
  descomposición $\mathcal M\cong\bigoplus_cM_{k_c}\otimes1_{m_c}$ de un álgebra dada, normalmente
  para reducir el tamaño de programas semidefinidos con simetría. Esto es complementario: esos
  algoritmos encuentran los bloques; esta biblioteca calcula en ellos e implementa uno de ellos como
  paso de entrada (`decompose`, capítulo 4). En álgebra computacional
  existen descomposiciones simbólicas exactas de álgebras de grupo (por ejemplo, el paquete
  Wedderga de GAP).
- **La optimización polinómica no conmutativa** (NCSOStools, ncpol2sdpa) acota valores propios y
  trazas de polinomios en variables no conmutativas, para todos los tamaños de matriz a la vez,
  mediante relajaciones semidefinidas [5]. Responde a otra pregunta (una cota válida para todas
  las álgebras) con otros medios (relajación convexa).

Hasta donde sabemos, no existe ninguna biblioteca general en la que un álgebra de von Neumann de
dimensión finita sea el tipo de datos básico para el cálculo numérico por lotes y diferenciable.
La novedad no es, por tanto, matemática nueva, sino **cálculo eficiente**: la identificación de la
§0.1, implementada de modo que se preserve la estructura del álgebra y cada operación se ejecute
como núcleos por lotes (uno por cada tamaño de bloque distinto en las descomposiciones, dividido en un trozo por hilo
en una CPU de varios núcleos). Si conoce software relacionado, nos alegrará saberlo.

## 0.4 Cómo leer este curso

| capítulos | contenido |
|---|---|
| 1–3 | las matemáticas: álgebras de matrices, el teorema de estructura, trazas, estados y entropía, con demostraciones completas |
| 4–8 | de dónde surgen estas álgebras y qué se calcula en ellas: simetrías, matrices aleatorias, canales, productos tensoriales, dinámica |
| 9–12 | aplicaciones: cadenas de espines, computación cuántica, optimización, aprendizaje automático equivariante |

Cada capítulo enuncia resultados con sus demostraciones, muestra el código de la biblioteca, lo
contrasta con un resultado exacto e indica las herramientas estándar con las que compararlo.

## Referencias

1. M. M. Rams, G. Wójtowicz, A. Sinha, J. Hasik, *YASTN: Yet another symmetric tensor networks; A
   Python library for Abelian symmetric tensor network calculations*, SciPost Phys. Codebases 52
   (2025).
2. K. Murota, Y. Kanno, M. Kojima, S. Kojima, *A numerical algorithm for block-diagonal
   decomposition of matrix $*$-algebras*, Technical Report METR 2007-52, University of Tokyo (2007).
3. T. Maehara, K. Murota, *A numerical algorithm for block-diagonal decomposition of matrix
   $*$-algebras, Part II: general algorithm*, Optimization Online (2009).
4. E. de Klerk, C. Dobre, D. V. Pasechnik, *Numerical block diagonalization of matrix
   $*$-algebras with application to semidefinite programming*, Math. Program. 129 (2011) 91–111.
5. I. Klep, J. Povh, *Semidefinite programming and sums of hermitian squares of noncommutative
   polynomials*, J. Pure Appl. Algebra 214 (2010) 740–749.

**Siguiente:** [capítulo 1](01_matrix_algebras.md).
