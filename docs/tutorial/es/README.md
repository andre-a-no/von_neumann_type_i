[English](../README.md) · **Español** · [Français](../fr/README.md) · [Deutsch](../de/README.md) · [中文](../zh/README.md) · [日本語](../ja/README.md) · [Русский](../ru/README.md)

# Álgebras de von Neumann de dimensión finita: un curso con `torch_vn_algebra`

Este curso expone la teoría de las álgebras de von Neumann de dimensión finita —las álgebras
$\bigoplus_cM_{k_c}(\mathbb C)$— junto con sus aplicaciones en física cuántica, información
cuántica, teoría de matrices aleatorias y aprendizaje automático, y muestra en cada paso cómo
calcular con ellas en `torch_vn_algebra`.

Está escrito a la vez para dos tipos de lectores:

- **Los matemáticos** encontrarán demostraciones completas de todos los resultados que se
  utilizan. No se cita nada sin demostración, salvo el álgebra lineal estándar, el teorema
  fundamental del álgebra y hechos elementales de análisis y de probabilidad (límites,
  convexidad, la función generatriz de momentos gaussiana); cuando un resultado se enuncia
  antes de demostrarse, o solo se aplica, se indica explícitamente.
- **Los físicos y los ingenieros** encontrarán el significado físico de cada noción en
  recuadros titulados *Significado físico*, y código operativo que pueden adaptar.

Los capítulos 1–3 tienen la misma estructura: teoría con demostraciones; el código de la
biblioteca para cada resultado; una comprobación del código frente a un resultado exacto (un
teorema, una fórmula cerrada o un cálculo independiente); una sección *Con qué comparar* que
menciona las herramientas habituales para la misma tarea; y ejercicios con soluciones.

Todo el código del curso se ejecuta en la batería de pruebas (`tests/test_tutorial.py`), de
modo que cada ejemplo funciona con la versión actual de la biblioteca.

## Contenido

| | capítulo | resultados principales |
|---|---|---|
| 0 | [De qué trata este curso](00_overview.md) | el álgebra como tipo de datos: estructura de bloques = tensores por lotes; coste; qué existe y qué es nuevo |
| 1 | [Álgebras de matrices y cálculo funcional](01_matrix_algebras.md) | teorema espectral, cálculo funcional, descomposición polar, dualidad de la norma traza, un momento del GUE |
| 2 | [Álgebras de von Neumann de dimensión finita](02_von_neumann_algebras.md) | teorema del doble conmutante, teorema de estructura, centro, grupos de simetría y conmutantes, esperanzas condicionales |
| 3 | [Trazas, estados y entropía](03_traces_and_states.md) | clasificación de las trazas, densidades, entropía de los sectores, desigualdad de Klein, monotonía bajo la esperanza condicional sobre el centro, contraste de Michelson |
| 4 | Simetrías y representaciones *(en preparación)* | lema de Schur, sectores de carga, $U(1)$ y $SU(2)$, sistemas de espines |
| 5 | [Matrices aleatorias en un álgebra](05_random_matrices.md) | medida de Haar en $\prod_cU(k_c)$, cálculo de Weingarten, estadística de niveles y por qué no deben mezclarse los sectores |
| 6 | [Canales cuánticos](06_quantum_channels.md) | Kraus y Stinespring, duales respecto de una traza, aplicaciones entre sectores, decoherencia |
| 7 | [Productos tensoriales y entrelazamiento](07_tensor_products.md) | traza parcial, fórmula de Page, fusión de cargas |
| 8 | [Dinámica](08_dynamics.md) | ecuaciones de Schrödinger y de Lindblad, saltos entre sectores, postselección |
| 9 | Cadenas de espines *(en preparación)* | sectores de magnetización y de momento, Lanczos, comprobaciones con el ansatz de Bethe, desorden |
| 10 | Computación cuántica *(en preparación)* | subespacios libres de decoherencia, corrección de errores en álgebras de operadores, mezcladores con restricciones |
| 11 | Optimización y desigualdades de trazas *(en preparación)* | parametrizaciones con restricciones, muestreo frente a optimización |
| 12 | Aprendizaje automático equivariante *(en preparación)* | capas equivariantes como conmutantes |

## Requisitos previos

Álgebra lineal sobre $\mathbb C$ (valores propios, productos escalares, matrices unitarias),
probabilidad básica y Python. No se presupone ningún conocimiento de física ni de análisis
funcional.

## Cómo ejecutar el código

```bash
pip install -e .                # from the root of the repository
python -m pytest tests/test_tutorial.py -q
```

Todos los bloques de código de un capítulo pueden pegarse, en orden, en una misma sesión de
Python o en un cuaderno de Jupyter. Los bloques de un mismo capítulo comparten sus variables.

## Lecturas complementarias

M. Takesaki, *Theory of Operator Algebras I* (Springer, 2002); V. Jones,
*von Neumann algebras*, apuntes de clase; D. Petz, *Quantum Information Theory and Quantum
Statistics* (Springer, 2008); J. Watrous, *The Theory of Quantum Information* (Cambridge,
2018); M. L. Mehta, *Random Matrices* (Elsevier, 2004). Para el trasfondo físico en lenguaje
matemático, véase [Física para matemáticos](../../physics_for_mathematicians.es.md). El artículo
que acompaña a la biblioteca es [`paper/main.pdf`](https://github.com/andre-a-no/von_neumann_type_i/blob/main/paper/main.pdf).
