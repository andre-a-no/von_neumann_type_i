[English](../00_overview.md) · **Español** · [Français](../fr/00_overview.md) · [Deutsch](../de/00_overview.md) · [中文](../zh/00_overview.md) · [日本語](../ja/00_overview.md) · [Русский](../ru/00_overview.md)

# 0. De qué trata este curso

Las álgebras de von Neumann se estudian sobre todo con lápiz y papel. En dimensión finita, cada
una de ellas es $`*`$-isomorfa a una suma directa de álgebras de matrices completas (teorema 2.6),
de modo que en principio todo se puede calcular; pero cuando un álgebra así se trata
numéricamente, a menudo se aplana en una única matriz grande, y su estructura (los sectores, el
centro, las múltiples trazas) desaparece del cálculo. Este curso y la biblioteca
`torch_vn_algebra` adoptan el enfoque contrario: **el álgebra misma es el tipo de datos**, y los
cálculos en ella se realizan de forma eficiente.

<a id="s-0-1"></a>

## 0.1 La observación

Un álgebra de von Neumann $`\mathcal M`$ sobre $`\mathbb C^n`$ es, salvo un cambio de base unitario,
```math
\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c}
```
(teorema 2.6). Como álgebra abstracta es $`\bigoplus_cM_{k_c}(\mathbb C)`$: un elemento es una lista
de $`C`$ bloques cuadrados. Esta es exactamente la forma de un tensor con un *eje de sectores*: un
elemento de $`\mathcal M`$ es un arreglo de forma $`(C,k,k)`$ (con los bloques rellenados con ceros
hasta $`k=\max_ck_c`$; así lo almacena la biblioteca, sección 2.5), y $`B`$ elementos (muestras de
Monte Carlo, puntos de partida de una optimización, pasos de tiempo, realizaciones del desorden)
forman un arreglo de forma $`(B,C,k,k)`$. Las bibliotecas modernas de arreglos (PyTorch en CPU y
GPU) ofrecen rutinas por lotes (productos de matrices, descomposiciones espectrales,
descomposiciones QR) precisamente para estos «lotes de matrices pequeñas», que también usan las
redes neuronales.

*En la biblioteca:* un elemento de $`\mathcal M`$ es un operador de un `TypeIAlgebra`, almacenado como un único tensor de forma $`(B,C,k_{\max},k_{\max})`$ — [tutorial, sección 1.5](../../tutorial/es/01_matrix_algebras.md#s-1-5) y [tutorial, sección 2.5](../../tutorial/es/02_von_neumann_algebras.md#s-2-5).

Bajo esta identificación, las operaciones del álgebra son operaciones por lotes sobre arreglos.
La correspondencia es exacta (cada una de las seis primeras filas de la tabla es una identidad,
demostrada más abajo a partir de resultados del curso); solo las rutinas numéricas que ejecutan
las operaciones trabajan en coma flotante:

| en el álgebra $`\mathcal M`$ | sobre el arreglo de forma $`(B,C,k,k)`$ |
|---|---|
| producto $`AB`$, suma, adjunto | producto matricial por lotes, suma, traspuesta conjugada sobre $`(B,C)`$ |
| cálculo funcional $`f(A)`$ (capítulo 1) | descomposición espectral por lotes de los bloques |
| una traza $`\sum_cw_c\mathrm{Tr}\,A_c`$ (capítulo 3) | suma de la diagonal de cada bloque, suma ponderada sobre el eje de sectores |
| el centro, las proyecciones centrales (capítulo 2) | un número por sector |
| esperanza condicional sobre el centro | traza de cada bloque dividida entre $`k_c`$, propagada de nuevo al bloque |
| medida de Haar en el grupo unitario $`\prod_cU(k_c)`$ | descomposición QR por lotes de bloques gaussianos, con la corrección de fase del teorema 5.3 |
| estados, canales, dinámica (capítulos 3, 6, 8) | construidos a partir de las operaciones anteriores |

*Demostraciones de las seis primeras filas.* Por el teorema 2.6 existe un unitario $`V`$ con
$`\mathcal M=V\big(\bigoplus_cM_{k_c}\otimes1_{m_c}\big)V^*`$, y
$`\varphi\big(V(\bigoplus_cA_c\otimes1_{m_c})V^*\big)=\bigoplus_cA_c`$ es un $`*`$-isomorfismo de $`\mathcal M`$ sobre
$`\bigoplus_cM_{k_c}(\mathbb C)`$ (bien definido, pues $`A_c\otimes1_{m_c}`$ determina $`A_c`$; multiplicativo, pues
$`(A_c\otimes1)(B_c\otimes1)=A_cB_c\otimes1`$). La biblioteca almacena $`\varphi(A)`$ (sección 2.5).
(1) En $`\bigoplus_cM_{k_c}`$ las sumas, los productos y los adjuntos se toman bloque a bloque.
(2) Sean $`A\in\mathcal M`$ normal, $`\varphi(A)=\bigoplus_cA_c`$ y $`f`$ una función sobre $`\mathrm{spec}\,A`$. Cada $`A_c`$ es
normal, ya que $`\varphi(A)^*\varphi(A)=\varphi(A^*A)=\varphi(AA^*)=\varphi(A)\varphi(A)^*`$ se cumple bloque a bloque. Si
$`A_cx=\lambda x`$ con $`x\neq0`$, entonces $`V(0\oplus\dots\oplus x\otimes e_1\oplus\dots\oplus0)`$ ($`x\otimes e_1`$ en el
$`c`$-ésimo sumando) es un vector propio de $`A`$ para $`\lambda`$, de modo que $`\mathrm{spec}\,A_c\subseteq\mathrm{spec}\,A`$. Sea $`p`$ un
polinomio con $`p=f`$ sobre $`\mathrm{spec}\,A`$ (interpolación de Lagrange, como en la demostración del teorema 1.5). Para todo
$`X`$ normal con $`\mathrm{spec}\,X\subseteq\mathrm{spec}\,A`$ y $`X=U\Lambda U^*`$, se tiene $`f(X)=Uf(\Lambda)U^*=Up(\Lambda)U^*=p(X)`$ por
el teorema 1.5(b); así $`f(A)=p(A)`$, $`f(A_c)=p(A_c)`$, y como $`\varphi`$ es un homomorfismo unital,
```math
\varphi(f(A))=p(\varphi(A))=\bigoplus_cp(A_c)=\bigoplus_cf(A_c),\qquad f(A_c)=U_cf(\Lambda_c)U_c^*
```
para cualquier diagonalización $`A_c=U_c\Lambda_cU_c^*`$ (teorema 1.3): una descomposición espectral por bloque.
(3) Las trazas de $`\bigoplus_cM_{k_c}`$ son los funcionales $`\sum_cw_c\mathrm{Tr}\,A_c`$ (teorema 3.1).
(4) $`Z(\mathcal M)`$ está generado por las proyecciones centrales minimales $`P_c`$, y $`\varphi(P_c)`$ es $`1_{k_c}`$ en el bloque
$`c`$ y $`0`$ en los demás (teorema 2.6), de modo que $`\varphi(Z(\mathcal M))=\{\bigoplus_cz_c1_{k_c}:z\in\mathbb C^C\}`$.
(5) Proposición 2.8. (6) Corolario 5.4 (bloques independientes, cada uno con la medida de Haar en $`U(k_c)`$) y teorema 5.3
(cada bloque a partir de la descomposición QR de una matriz gaussiana). La última fila es solo una indicación: los estados
se representan mediante densidades (teorema 3.2), los canales entre sectores mediante operadores de Kraus entre bloques
(proposición 6.4), y la dinámica cerrada mediante $`e^{-itH}=\bigoplus_ce^{-itH_c}`$ (proposición 8.2).

*En la biblioteca:* los productos y funciones como `sqrt` y `apply_function` — [tutorial, sección 1.5](../../tutorial/es/01_matrix_algebras.md#s-1-5); las trazas, p. ej. `Tr_blunt` — [tutorial, sección 3.1](../../tutorial/es/03_traces_and_states.md#s-3-1); el centro, `TypeIAlgebra.central` — [tutorial, sección 2.5](../../tutorial/es/02_von_neumann_algebras.md#s-2-5); la esperanza condicional sobre el centro, `channels.center_expectation` — [tutorial, sección 2.6](../../tutorial/es/02_von_neumann_algebras.md#s-2-6); los unitarios de Haar, `random_unitary_operator` — [tutorial, sección 5.1](../../tutorial/es/05_random_matrices.md#s-5-1).

*Cálculo:* la fila 2 de la tabla sobre $`\bigoplus_{c=1}^{64}M_{16}(\mathbb C)`$ — [tutorial, sección 0.2](../../tutorial/es/00_overview.md#s-0-2); lo que cuesta calcular bloque a bloque — [tutorial, sección 0.1](../../tutorial/es/00_overview.md#s-0-1).

<a id="s-0-2"></a>

## 0.2 Un primer cálculo

*Cálculo:* la raíz cuadrada y la traza de un elemento positivo aleatorio de $`\bigoplus_{c=1}^{64}M_{16}(\mathbb C)`$, calculadas bloque a bloque y en $`M_{1024}(\mathbb C)`$ — [tutorial, sección 0.2](../../tutorial/es/00_overview.md#s-0-2).

<a id="s-0-3"></a>

## 0.3 Qué existe y qué es nuevo

*En la biblioteca:* el software relacionado y lo que la biblioteca le añade — [tutorial, sección 0.3](../../tutorial/es/00_overview.md#s-0-3).

<a id="s-0-4"></a>

## 0.4 Cómo leer este curso

| capítulos | contenido |
|---|---|
| 1–3 | las matemáticas: álgebras de matrices, el teorema de estructura, trazas, estados y entropía, con demostraciones completas |
| 4–8 | de dónde surgen estas álgebras y qué se calcula en ellas: simetrías, matrices aleatorias, canales, productos tensoriales, dinámica |
| 9–12 | aplicaciones: cadenas de espines, computación cuántica, optimización, aprendizaje automático equivariante |
| 13 | fundamentos: álgebras de Jordan y por qué la mecánica cuántica usa *-álgebras |
| apéndices | demostraciones completas de los resultados auxiliares que los capítulos citan por su número (por ejemplo, el teorema A.24) |

Cada uno de los capítulos 1–13 enuncia resultados con demostraciones completas. Los requisitos previos
(álgebra lineal, análisis real, teoría de la medida, probabilidad elemental) se enumeran en el
[README](../../lecture_notes/es/README.md#prerequisites); todo lo demás se demuestra en los apuntes de clase, y los resultados de
otros capítulos y de los apéndices se citan por su número. Los recuadros titulados *Significado físico* son
interpretación física, no teoremas de los apuntes de clase.

Cada capítulo tiene un gemelo en el [tutorial](../../tutorial/es/00_overview.md), que muestra el código de la biblioteca, lo contrasta numéricamente con los resultados demostrados aquí e indica las herramientas estándar con las que compararlo. Las líneas *Cálculo:* y *En la biblioteca:* que siguen a un resultado enlazan con el lugar del tutorial donde se calcula o se usa.

**Siguiente:** [capítulo 1](01_matrix_algebras.md).
