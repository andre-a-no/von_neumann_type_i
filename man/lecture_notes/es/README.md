[English](../README.md) · **Español** · [Français](../fr/README.md) · [Deutsch](../de/README.md) · [中文](../zh/README.md) · [日本語](../ja/README.md) · [Русский](../ru/README.md)

# Álgebras de von Neumann de dimensión finita: apuntes de clase

Este curso expone la teoría de las álgebras de von Neumann de dimensión finita —las álgebras
$`\bigoplus_cM_{k_c}(\mathbb C)`$— junto con sus aplicaciones en física cuántica, información
cuántica, teoría de matrices aleatorias y aprendizaje automático.

Está escrito a la vez para dos tipos de lectores:

- **Los matemáticos** encontrarán un texto autocontenido: todo teorema, fórmula, constante y
  desigualdad que se enuncia se demuestra en el curso —en el propio capítulo o en uno de los
  apéndices A–F, citado por su número— a partir de los requisitos previos enumerados más
  abajo. Los números que solo se calculan (por ejemplo, la estadística de niveles para $`n`$
  grande) se señalan como *observación numérica* y nunca se usan como premisas ni como
  objetivos de las pruebas. Las referencias son atribución histórica o lectura complementaria,
  nunca un sustituto de una demostración.
- **Los físicos y los ingenieros** encontrarán el significado físico de cada noción en
  recuadros titulados *Significado físico*. Estos recuadros son interpretación física, no teoremas del curso:
  nada de lo que se demuestra en otras partes depende de ellos. El código operativo que pueden
  adaptar está en el [tutorial](../../tutorial/es/README.md).

Los capítulos 1–13 tienen la misma estructura: definiciones y resultados con demostraciones completas,
ejemplos resueltos y ejercicios con soluciones. Después de cada resultado que la biblioteca `torch_vn_algebra`
calcula o comprueba, una línea *Cálculo:* o *En la biblioteca:* enlaza con el lugar del
[tutorial](../../tutorial/es/README.md) donde se hace.

## Contenido

| | capítulo | resultados principales |
|---|---|---|
| 0 | [De qué trata este curso](00_overview.md) | el álgebra como tipo de datos: estructura de bloques = tensores por lotes; coste; qué existe y qué es nuevo |
| 1 | [Álgebras de matrices y cálculo funcional](01_matrix_algebras.md) | teorema espectral, cálculo funcional, descomposición polar, dualidad de la norma traza, un momento del GUE |
| 2 | [Álgebras de von Neumann de dimensión finita](02_von_neumann_algebras.md) | teorema del doble conmutante, teorema de estructura, centro, grupos de simetría y conmutantes, esperanzas condicionales |
| 3 | [Trazas, estados y entropía](03_traces_and_states.md) | clasificación de las trazas, densidades, entropía de los sectores, desigualdad de Klein, monotonía bajo la esperanza condicional sobre el centro, contraste de Michelson |
| 4 | [Simetrías y representaciones](04_symmetries.md) | lema de Schur, sectores de carga, $`U(1)`$ y $`SU(2)`$, sistemas de espines, tipo real, complejo y cuaterniónico |
| 5 | [Matrices aleatorias en un álgebra](05_random_matrices.md) | medida de Haar en $`\prod_cU(k_c)`$, cálculo de Weingarten, estadística de niveles y por qué no deben mezclarse los sectores |
| 6 | [Canales cuánticos](06_quantum_channels.md) | Kraus y Stinespring, duales respecto de una traza, aplicaciones entre sectores, decoherencia |
| 7 | [Productos tensoriales y entrelazamiento](07_tensor_products.md) | traza parcial, fórmula de Page, fusión de cargas |
| 8 | [Dinámica](08_dynamics.md) | ecuaciones de Schrödinger y de Lindblad, saltos entre sectores, postselección |
| 9 | [Cadenas de espines](09_spin_chains.md) | sectores de magnetización y de momento, Lanczos, comprobaciones con el ansatz de Bethe, desorden |
| 10 | [Computación cuántica](10_quantum_computing.md) | subespacios libres de decoherencia, corrección de errores en álgebras de operadores, mezcladores con restricciones |
| 11 | [Optimización y desigualdades de trazas](11_optimisation.md) | parametrizaciones con restricciones, muestreo frente a optimización |
| 12 | [Aprendizaje automático equivariante](12_equivariant_learning.md) | capas equivariantes como conmutantes |
| 13 | [Álgebras de Jordan: por qué la mecánica cuántica usa *-álgebras](13_jordan_algebras.md) | producto de Jordan y representación cuadrática, factores de espín y el cono de Lorentz, la clasificación de Jordan–von Neumann–Wigner, la traspuesta y la traspuesta parcial, una desigualdad de trazas en forma de Jordan |
| A | [Hechos de análisis y de probabilidad](A_analysis.md) | series de potencias y conjuntos nulos, vectores gaussianos, medidas esféricas, desigualdades de valores propios, Isserlis–Wick |
| B | [Dualidad de Schur–Weyl y cálculo de Weingarten](B_schur_weyl.md) | dualidad de Schur–Weyl, funciones de Weingarten, teorema de Wigner–Eckart, medida de Haar en $`O(k)`$ |
| C | [Densidades de valores propios y la fórmula de Page](C_eigenvalue_densities.md) | fórmula de integración de Weyl, GOE/GUE, densidades de Wishart e inducidas, fórmula de Page |
| D | [Generadores de semigrupos dinámicos cuánticos](D_gksl.md) | el teorema GKSL, generadores en sumas directas, teorema de Frigerio, Bény–Kempf–Kribs |
| E | [Clasificación de las álgebras de Jordan euclídeas](E_jordan_classification.md) | descomposición de Peirce, teoremas de Hurwitz y de Frobenius, coordinatización, el álgebra de Albert, el teorema de Jordan–von Neumann–Wigner |
| F | [Conos simétricos, isomorfismos de orden y aplicaciones de Jordan](F_cones_and_order.md) | conos de cuadrados, el teorema de Koecher–Vinberg, el teorema de Kadison, automorfismos de Jordan de $`\bigoplus_cM_{k_c}`$ |

<a id="prerequisites"></a>

## Requisitos previos

El curso presupone únicamente lo siguiente y demuestra todo lo demás.

- Álgebra lineal sobre $`\mathbb R`$ y $`\mathbb C`$ (de dimensión finita) y el teorema
  fundamental del álgebra.
- Análisis real y cálculo en varias variables estándar: límites, continuidad, compacidad en
  $`\mathbb R^n`$, diferenciación, el teorema del valor medio, la fórmula de Taylor, los
  teoremas de la función inversa y de la función implícita.
- Teoría estándar de la medida y de la integración: medida de Lebesgue, medidas producto y
  Fubini–Tonelli, convergencia dominada y monótona, cambio de variables en $`\mathbb R^n`$,
  medidas imagen, densidades, independencia.
- Probabilidad elemental: esperanza, varianza, independencia, la distribución gaussiana y su
  función generatriz de momentos, las desigualdades de Markov y de Chebyshov, la ley de los
  grandes números.

No se presupone ningún conocimiento de física ni de análisis funcional.

## Lecturas complementarias

M. Takesaki, *Theory of Operator Algebras I* (Springer, 2002); V. Jones,
*von Neumann algebras*, apuntes de clase; D. Petz, *Quantum Information Theory and Quantum
Statistics* (Springer, 2008); J. Watrous, *The Theory of Quantum Information* (Cambridge,
2018); M. L. Mehta, *Random Matrices* (Elsevier, 2004).

La biblioteca `torch_vn_algebra` y su uso en el contexto de estos apuntes son el tema del
[tutorial](../../tutorial/es/README.md). El artículo que acompaña a la biblioteca es
[`man/article/main.pdf`](https://github.com/andre-a-no/von_neumann_type_i/blob/main/man/article/main.pdf).
