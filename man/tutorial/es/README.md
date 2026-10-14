[English](../README.md) · **Español** · [Français](../fr/README.md) · [Deutsch](../de/README.md) · [中文](../zh/README.md) · [日本語](../ja/README.md) · [Русский](../ru/README.md)

# Tutorial: torch_vn_algebra en el contexto de los apuntes de clase

Este tutorial muestra cómo calcular con las álgebras de von Neumann de dimensión finita $\bigoplus_cM_{k_c}(\mathbb C)$
en la biblioteca `torch_vn_algebra`, en el contexto de los [apuntes de clase](../../lecture_notes/es/README.md). Cada tutorial acompaña
al capítulo de los apuntes de clase con el mismo número: muestra el código de la biblioteca para los resultados
demostrados allí, comprueba el código frente a un resultado exacto (un teorema, una fórmula cerrada o un cálculo
independiente) y menciona, en una sección *Con qué comparar*, las herramientas habituales para la misma tarea. Cada
sección empieza con una línea *Apuntes:* que enlaza con la sección y los resultados que usa. Los números que solo se
calculan (por ejemplo, la estadística de niveles para $n$ grande) se señalan como *observación numérica* y nunca se
usan como premisas ni como objetivos de las pruebas. Además de los [requisitos previos](../../lecture_notes/es/README.md#prerequisites)
de los apuntes de clase, el tutorial presupone Python.

Todo el código del tutorial se ejecuta en la batería de pruebas (`tests/test_man.py`), de
modo que cada ejemplo funciona con la versión actual de la biblioteca.

## Cómo ejecutar el código

```bash
pip install -e .                # from the root of the repository
python -m pytest tests/test_man.py -q
```

Todos los bloques de código de un capítulo pueden pegarse, en orden, en una misma sesión de
Python o en un cuaderno de Jupyter. Los bloques de un mismo capítulo comparten sus variables.

## Contenido

| | capítulo | resultados principales |
|---|---|---|
| 0 | [De qué trata este curso](00_overview.md) | el álgebra como tipo de datos: estructura de bloques = tensores por lotes; coste; qué existe y qué es nuevo |
| 1 | [Álgebras de matrices y cálculo funcional](01_matrix_algebras.md) | teorema espectral, cálculo funcional, descomposición polar, dualidad de la norma traza, un momento del GUE |
| 2 | [Álgebras de von Neumann de dimensión finita](02_von_neumann_algebras.md) | teorema del doble conmutante, teorema de estructura, centro, grupos de simetría y conmutantes, esperanzas condicionales |
| 3 | [Trazas, estados y entropía](03_traces_and_states.md) | clasificación de las trazas, densidades, entropía de los sectores, desigualdad de Klein, monotonía bajo la esperanza condicional sobre el centro, contraste de Michelson |
| 4 | [Simetrías y representaciones](04_symmetries.md) | lema de Schur, sectores de carga, $U(1)$ y $SU(2)$, sistemas de espines, tipo real, complejo y cuaterniónico |
| 5 | [Matrices aleatorias en un álgebra](05_random_matrices.md) | medida de Haar en $\prod_cU(k_c)$, cálculo de Weingarten, estadística de niveles y por qué no deben mezclarse los sectores |
| 6 | [Canales cuánticos](06_quantum_channels.md) | Kraus y Stinespring, duales respecto de una traza, aplicaciones entre sectores, decoherencia |
| 7 | [Productos tensoriales y entrelazamiento](07_tensor_products.md) | traza parcial, fórmula de Page, fusión de cargas |
| 8 | [Dinámica](08_dynamics.md) | ecuaciones de Schrödinger y de Lindblad, saltos entre sectores, postselección |
| 9 | [Cadenas de espines](09_spin_chains.md) | sectores de magnetización y de momento, Lanczos, comprobaciones con el ansatz de Bethe, desorden |
| 10 | [Computación cuántica](10_quantum_computing.md) | subespacios libres de decoherencia, corrección de errores en álgebras de operadores, mezcladores con restricciones |
| 11 | [Optimización y desigualdades de trazas](11_optimisation.md) | parametrizaciones con restricciones, muestreo frente a optimización |
| 12 | [Aprendizaje automático equivariante](12_equivariant_learning.md) | capas equivariantes como conmutantes |
| 13 | [Álgebras de Jordan: por qué la mecánica cuántica usa *-álgebras](13_jordan_algebras.md) | producto de Jordan y representación cuadrática, factores de espín y el cono de Lorentz, la clasificación de Jordan–von Neumann–Wigner, la traspuesta y la traspuesta parcial, una desigualdad de trazas en forma de Jordan |
| A | [Hechos de análisis y de probabilidad](A_analysis.md) | series de potencias y conjuntos nulos, vectores gaussianos, medidas esféricas, desigualdades de valores propios, Isserlis–Wick |
| B | [Dualidad de Schur–Weyl y cálculo de Weingarten](B_schur_weyl.md) | dualidad de Schur–Weyl, funciones de Weingarten, teorema de Wigner–Eckart, medida de Haar en $O(k)$ |
| C | [Densidades de valores propios y la fórmula de Page](C_eigenvalue_densities.md) | fórmula de integración de Weyl, GOE/GUE, densidades de Wishart e inducidas, fórmula de Page |
| D | [Generadores de semigrupos dinámicos cuánticos](D_gksl.md) | el teorema GKSL, generadores en sumas directas, teorema de Frigerio, Bény–Kempf–Kribs |
| E | [Clasificación de las álgebras de Jordan euclídeas](E_jordan_classification.md) | descomposición de Peirce, teoremas de Hurwitz y de Frobenius, coordinatización, el álgebra de Albert, el teorema de Jordan–von Neumann–Wigner |
| F | [Conos simétricos, isomorfismos de orden y aplicaciones de Jordan](F_cones_and_order.md) | conos de cuadrados, el teorema de Koecher–Vinberg, el teorema de Kadison, automorfismos de Jordan de $\bigoplus_cM_{k_c}$ |

Las matemáticas —definiciones, teoremas y demostraciones completas— están en los [apuntes de clase](../../lecture_notes/es/README.md),
junto con los requisitos previos, las lecturas complementarias y el artículo que acompaña a la biblioteca.
