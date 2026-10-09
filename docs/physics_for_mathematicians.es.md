[English](physics_for_mathematicians.md) · **Español** · [Français](physics_for_mathematicians.fr.md) · [Deutsch](physics_for_mathematicians.de.md) · [中文](physics_for_mathematicians.zh.md) · [日本語](physics_for_mathematicians.ja.md) · [Русский](physics_for_mathematicians.ru.md)

# Física para matemáticos: un diccionario para `torch_vn_algebra`

*Traducción de la versión en inglés; en caso de discrepancia, prevalece la versión en inglés.*

Esta guía está dirigida a lectores que vienen de las álgebras de operadores y quieren saber qué significan
los modelos físicos que aparecen en la biblioteca y en el artículo. Todo se formula en términos del álgebra
M = ⊕_c M_{k_c}(C) que actúa sobre H = ⊕_c C^{k_c}. Para cada noción física indicamos qué es en M y
dónde se encuentra en la biblioteca.

---

## 0. La idea principal

Una simetría física es una representación unitaria π de un grupo G sobre H. Todo lo que es compatible con
la simetría (el hamiltoniano, la evolución temporal, los observables sujetos a una regla de superselección)
está en el **conmutante** π(G)'. Por el teorema del doble conmutante y el lema de Schur,

    H ≅ ⊕_λ V_λ ⊗ C^{m_λ},     π(G)' ≅ ⊕_λ M_{m_λ}(C),

de modo que π(G)' es un álgebra de von Neumann de tipo I de dimensión finita. Sus sectores son las
representaciones irreducibles λ, los tamaños de los bloques son las multiplicidades m_λ y el centro está
formado por funciones de las «cargas». Así pues, un físico que modela un sistema con una simetría calcula
en un álgebra de la forma M = ⊕_c M_{k_c}(C), normalmente sin llamarla así. Para una única magnitud conservada Q (un
operador autoadjunto), los sectores son los autovalores de Q y los bloques son sus autoespacios.

Esto responde también a la pregunta «¿por qué no usar una única matriz grande?». La inclusión M ⊂ M_D
(D = Σk_c) es exacta, pero se pierde
- **el coste**: Σk_c³ en lugar de D³; para C bloques iguales esto supone un factor C² en tiempo;
- **la medida**: el grupo unitario de M es Π U(k_c), un conjunto de medida nula en U(D). Un unitario
  aleatorio según la medida de Haar en U(D) no está en M, y un paso de gradiente en M_D se sale de M;
- **las trazas**: en M_D la traza es única salvo un factor; en M las trazas fieles forman la familia
  de C parámetros Σ w_c Tr_c. La preservación de la traza y las aplicaciones duales dependen de los pesos;
- **el centro**: las proyecciones de sector son la parte clásica del sistema (sección 2). Las
  probabilidades condicionadas y la posselección son operaciones con el centro, y con una sola matriz
  hay que reconstruirlas a partir de los índices de bloque.

El caso C = 1 también está soportado: es álgebra lineal por lotes en M_n.

---

## 1. Estados, observables, medidas

| física | matemáticas | en la biblioteca |
|---|---|---|
| observable (energía, espín, número de partículas) | A ∈ M autoadjunto | `Operator` |
| resultados posibles de una medida | espectro de A | `eigenvalues`, `eigh` |
| estado (mixto) | funcional positivo normal ω, ω(1) = 1; densidad ρ ≥ 0, Tr ρ = 1 | `DensityMatrix` |
| estado puro, «función de onda» | vector ψ en un único sector C^{k_c}, ρ = \|ψ⟩⟨ψ\| (una proyección minimal de M; un ψ repartido entre varios sectores da sobre M el estado mixto Σ_c P_c\|ψ⟩⟨ψ\|P_c) | `vector_in_sector`, `basis_state` |
| valor esperado | ω(A) = Tr(ρA) | `expectation` |
| probabilidad del resultado λ | ω(P_λ), con P_λ la proyección espectral | `apply_function` |
| estado tras una medida (regla de Lüders) | ρ ↦ PρP / Tr(Pρ) | `lueders_update`, `condition_on` |
| entropía de von Neumann | S(ρ) = −Tr ρ log ρ | `entropy` |
| temperatura, estado térmico | ρ = e^{−βH}/Tr e^{−βH} (estado de Gibbs = estado KMS) | `gibbs_state` |
| temperatura infinita | ρ = 1/D, D = Σk_c (Tr_blunt normalizada) | `maximally_mixed_state` |
| el mismo peso para cada sector | densidad de τ_vN, ⊕_c 1_c/(C k_c) | `tracial_state` |

**Significado físico de las tres trazas.**
- Tr_blunt (pesos 1) es la traza física sobre H. Normalizada, es el estado a temperatura infinita:
  los D estados de la base son equiprobables y el sector c tiene probabilidad k_c/D.
- Tr_norm (pesos 1/k_c): el estado uniforme (microcanónico) dentro de cada sector por separado.
- τ_vN (pesos 1/(C k_c)): la mezcla uniforme de los estados microcanónicos de los sectores. Cada
  *sector* (no cada estado de la base) tiene probabilidad 1/C.

La restricción de cualquier estado al centro es una distribución de probabilidad clásica sobre los sectores
(`sector_probabilities`). Esto es lo que los físicos llaman «la probabilidad de la carga q».

---

## 2. Superselección, decoherencia, lo clásico y lo cuántico

- **Regla de superselección** (Wick–Wightman–Wigner): las superposiciones de estados con distinta carga
  (carga eléctrica, paridad fermiónica) no son observables. Matemáticamente, los observables admisibles
  forman el conmutante de la carga, es decir, el álgebra M en lugar de todo B(H).
- **El centro de M** está formado por magnitudes que pueden medirse sin perturbar el sistema y que conmutan
  con todo: las variables «clásicas». Si todos los k_c = 1, el sistema es puramente clásico (un álgebra
  conmutativa, es decir, un espacio de probabilidad sobre C puntos). Si C = 1, es puramente cuántico. El
  caso general es híbrido.
- **Decoherencia / einselección** (Zurek): la interacción con un entorno suprime los bloques no diagonales
  ρ_{cd}, c ≠ d. Esto es la esperanza condicional E: B(H) → M, ρ ↦ Σ P_c ρ P_c, y el álgebra efectiva
  de observables se reduce a M. En la biblioteca los operadores ya están en M, de modo que los cálculos
  parten del resultado de esta esperanza condicional.

---

## 3. Dinámica de sistemas cerrados

| física | matemáticas | en la biblioteca |
|---|---|---|
| hamiltoniano | H ∈ M autoadjunto | — |
| ecuación de Schrödinger i dψ/dt = Hψ | ψ(t) = e^{−itH}ψ(0) | `dynamics.schrodinger`, `propagator` |
| ecuación de von Neumann dρ/dt = −i[H, ρ] | ρ(t) = U ρ U*, U = e^{−itH} | `von_neumann` |
| imagen de Heisenberg | A(t) = U* A U (un automorfismo de M) | `solve_operator_ode` |
| magnitud conservada | Q con [H, Q] = 0 | sectores |

Si H ∈ M, entonces e^{−itH} ∈ M y la dinámica nunca sale de los sectores. Cada sector evoluciona de forma
independiente y la exponencial se calcula bloque a bloque.

---

## 4. Sistemas abiertos y canales

- Un **canal cuántico** es una aplicación completamente positiva que preserva la traza (imagen de
  Schrödinger). Su adjunta (imagen de Heisenberg) es completamente positiva y unital. Forma de Kraus:
  Φ(ρ) = Σ_i K_i ρ K_i*. En la biblioteca: `Channel`.
- Un **canal entre sectores** (`InterSectorChannel`) describe un proceso que cambia la carga: pérdida de
  partículas, desintegración, un electrón que sale por efecto túnel. Sus bloques de Kraus son
  K^{dc}: C^{k_c} → C^{m_d}. La preservación de la traza y la aplicación dual dependen de los pesos de la
  traza en la entrada y en la salida; de ahí los factores w_out/w_in.
- La **ecuación de Lindblad** (GKSL) genera un semigrupo de canales e^{tL}:
  dρ/dt = −i[H, ρ] + Σ_j γ_j (L_j ρ L_j* − ½{L_j* L_j, ρ}).
  Los L_j son los «operadores de salto»; cada uno describe una forma en que el entorno actúa sobre el sistema.
  - **T1** (relajación de energía, amortiguamiento de amplitud): L = σ⁻ = |0⟩⟨1|, un qubit decae del
    estado excitado al fundamental. El número de excitaciones disminuye en una unidad, así que se trata de
    un salto **entre** sectores. Ejemplo: `site_amplitude_damping`, `lowering`.
  - **Desfase puro** (*pure dephasing*, tiempo T_φ): L = σ^z. Es diagonal y no cambia el sector; solo destruye las fases.
    Se trata de un salto **dentro** de un sector.
  - T1 y T2, con 1/T2 = 1/(2T1) + 1/T_φ, son los tiempos de coherencia estándar que se indican para cualquier ordenador cuántico.

---

## 5. Cadenas de espines

- **Espín ½ en un sitio**: H_site = C², base |↑⟩, |↓⟩ (o |1⟩, |0⟩). Matrices de Pauli σ^x, σ^y, σ^z;
  S^a = σ^a/2; S^± = S^x ± iS^y son los operadores de subida y de bajada.
- **Una cadena de L espines**: H = (C²)^{⊗L}, dim = 2^L. Este crecimiento exponencial es la razón por la
  que los físicos necesitan las simetrías.
- **Hamiltoniano XXZ**: H = Σ_⟨ij⟩ [J/2 (S_i⁺S_j⁻ + S_i⁻S_j⁺) + JΔ S_i^z S_j^z] + Σ_i h_i S_i^z.
  El primer término intercambia ↑↓ vecinos (un «salto de partícula»), el segundo es su interacción y el
  tercero un campo externo. Δ = 1 es el modelo de Heisenberg (un imán); Δ = 0 equivale a fermiones libres
  (mediante la transformación de Jordan–Wigner), lo que proporciona comprobaciones exactas.
- **Carga conservada**: S^z_tot, o N = número de ↑ (el peso de Hamming de una cadena de bits). Todos los
  términos de H la conservan, de modo que H ∈ ⊕_{N=0}^{L} M_{C(L,N)}: los bloques tienen tamaños binomiales.
- **Fusión de sectores** (`fused_tensor_product`): bajo el producto tensorial las cargas se suman. Es la
  regla de fusión de U(1): (⊕_q M_{n_q}) ⊗ (⊕_r M_{m_r}) ⊂ ⊕_s M_{Σ_{q+r=s} n_q m_r}. `SpinChain` construye
  su base de esta manera, sitio a sitio.
- **Momento** (`momentum_algebra`): un anillo tiene además la simetría de traslación T (el grupo Z_L). Los
  autoespacios conjuntos (N, k), con e^{ik} autovalor de T, dan una descomposición más fina, con bloques
  unas L veces más pequeños.
- El **estado fundamental** es el autovector de H con el menor autovalor, es decir, el comportamiento a
  temperatura cero. Para bloques grandes se obtiene con el método de Lanczos (`krylov`).
- Qué se comprueba en el artículo y por qué son pruebas estándar:
  - energías del estado fundamental del anillo de Heisenberg (valores exactos conocidos, ansatz de Bethe);
  - la **regla de signos de Marshall**: para el modelo de Heisenberg antiferromagnético (J > 0) en una red
    bipartita (un anillo de longitud par), el estado más bajo de cada sector N tiene signos
    (−1)^{número de ↑ en una subred}. Es un teorema, lo que la convierte en una prueba limpia;
  - **anidamiento SU(2)**: el modelo de Heisenberg tiene simetría SU(2) completa, de modo que para N < L/2
    el espectro del sector N está contenido en el del sector N+1 (multipletes);
  - **disolución de una pared de dominio**: |↑…↑↓…↓⟩ con Δ = 0 se ensancha con un perfil conocido;
  - **fórmula de Peschel**: la entropía de entrelazamiento de fermiones libres a partir de la matriz de
    correlación.
- **Desorden y caos** (`random_field_heisenberg`, `level_spacing_ratio`): campos aleatorios
  h_i ∈ [−W, W]. Para un desorden moderado (0 < W ≲ 2; con W = 0 la cadena sin desorden es integrable
  por ansatz de Bethe) el sistema es «caótico» y los niveles se repelen como en el GOE, ⟨r⟩ ≈ 0.53. Para W grande el sistema se localiza (localización de muchos cuerpos, MBL) y los niveles son
  independientes como en un proceso de Poisson, ⟨r⟩ ≈ 0.386. La estadística de niveles solo tiene sentido
  **dentro de un sector**: mezclar sectores produce artificialmente estadística de Poisson. Este es un buen
  argumento físico para trabajar en el álgebra con sus sectores.

---

## 6. Entrelazamiento y subsistemas

- Un sistema compuesto A+B: H_A ⊗ H_B, álgebra M_A ⊗ M_B.
- El **estado reducido** ρ_A = Tr_B ρ es la restricción del funcional ω a la subálgebra M_A ⊗ 1.
  En términos de densidades es la traza parcial, la aplicación dual de la inclusión A ↦ A ⊗ 1 (la
  esperanza condicional sobre M_A ⊗ 1 es X ↦ (Tr_B X) ⊗ 1/d_B). En la biblioteca:
  `partial_trace`, `reduced_state`.
- La **entropía de entrelazamiento** S(ρ_A) de un ρ puro mide la correlación cuántica entre A y B.
- **Fórmula de Page**: la entropía media exacta de un subsistema de un estado puro aleatorio. Sirve para
  comprobar que el muestreador sigue la distribución de Haar.
- **Entrelazamiento resuelto por simetría**: si el estado es simétrico, [ρ, Q_A + Q_B] = 0 (p. ej., un autoestado
  de la carga total), ρ_A es diagonal por bloques y su entropía se descompone en una parte «clásica» (la entropía de la distribución sobre el centro) y una
  parte «cuántica» (la entropía media de los bloques). Es un uso directo del centro.

---

## 7. Computación cuántica

- **Qubit** = C², n qubits = (C²)^{⊗n}; una **puerta** es un unitario; un **circuito** es un producto de
  puertas; una **medida** en la base computacional usa las proyecciones espectrales de σ^z.
- **QAOA** (*quantum approximate / alternating operator ansatz*): se alternan e^{−iγ H_C} (H_C es
  diagonal y codifica la función objetivo) y e^{−iβ H_M} (el «mezclador»). Para problemas con la
  restricción «exactamente k unos» (elegir k elementos de entre n) se usa el **mezclador XY**
  Σ(X_iX_j + Y_iY_j) = 2Σ(σ⁺σ⁻ + h.c.). Este conserva el peso de Hamming, de modo que todo el circuito
  está en ⊕_N M_{C(n,N)} y las soluciones factibles son exactamente el sector N = k.
- **Ruido**: T1 saca el estado del sector factible (se pierde un uno); el desfase no.
- **Posselección**: se mide el peso de Hamming y se descartan las ejecuciones incorrectas. Es la regla de
  Lüders con una proyección **central** P_k: ρ ↦ P_kρP_k / Tr(P_kρ), con probabilidad de éxito ω(P_k).
  Es lo que hace `ex9_qaoa_xy_mixer.py`. Con ruido T1 puro y tasas iguales en todos los qubits, la
  posselección devuelve exactamente el estado ideal (F = 1): todo salto T1 sale del sector y el
  amortiguamiento de la evolución sin saltos es el mismo para todos los estados del sector (con tasas
  distintas, F queda ligeramente por debajo de 1).
- Subespacios libres de decoherencia y corrección de errores en álgebras de operadores: la información
  protegida se codifica en un bloque (o en un factor tensorial de un bloque) de un álgebra que el ruido no
  afecta. También este es el lenguaje de las álgebras de tipo I (Knill–Laflamme–Viola, Bény–Kempf–Kribs).

---

## 8. Matrices aleatorias, reales y complejas

Esto está directamente relacionado con el parámetro `complex_valued`.
- Una **simetría antiunitaria** T (la inversión temporal, o una combinación como una reflexión compuesta
  con la conjugación compleja) determina el cuerpo. Si T² = 1, H es **real** en una base adecuada: GOE, conjuntos ortogonales, COE. Basta con aritmética
  real, que reduce la memoria a la mitad y es más rápida (véase la tabla de tiempos del artículo).
- Sin ninguna simetría de este tipo (flujos magnéticos orbitales, amplitudes de salto complejas), H es
  complejo: GUE, CUE. Los números complejos son imprescindibles. Un campo de Zeeman h_i S^z_i por sí solo
  mantiene H real.
- Si T² = −1 (espín semientero con acoplamiento espín–órbita), el caso es cuaterniónico: GSE, CSE.
- El modelo XXZ con campos reales es real. Los sectores de momento con k ∉ {0, π} son complejos en la
  base de ondas planas, por lo que requieren aritmética compleja; con simetría de reflexión, su estadística
  de niveles sigue siendo GOE. La dinámica e^{−itH} es siempre compleja, incluso para H real. Esta es la tabla de «qué basta y cuándo»
  del artículo.
- La doble precisión es necesaria cuando importan diferencias pequeñas: espectros degenerados, identidades
  comprobadas hasta 1e−12, evoluciones temporales largas, Lanczos con muchos pasos. La precisión simple
  (y TF32 en los tensor cores) basta para estadísticas sobre muchas muestras.

---

## 9. Aprendizaje automático

Una capa lineal equivariante es una aplicación W: V → V que conmuta con una representación de un grupo
(traslaciones en las convoluciones, permutaciones en las redes de grafos, rotaciones en las moléculas). Por
el lema de Schur, el espacio de tales W es ⊕_λ M_{m_λ} con multiplicidades m_λ (sobre C; para capas reales
los bloques son M_{m_λ}(D_λ) con D_λ = R, C o H, la misma tricotomía que en la sección 8), de modo que una
capa equivariante se parametriza mediante un elemento de un álgebra de la forma ⊕_c M_{k_c}(C). La biblioteca proporciona
inicialización aleatoria con la medida adecuada, restricciones espectrales y gradientes para este tipo de
parametrizaciones.

---

## 10. Alcance: qué es y qué no es la biblioteca

- La biblioteca **no** produce física nueva. Es una herramienta validada para cálculos en álgebras con
  sectores. Los modelos físicos del artículo son sistemas conocidos con respuestas conocidas, que se usan
  para comprobar la herramienta.
- Para un único sistema grande (un solo estado fundamental con L = 30, redes tensoriales, simulación de
  circuitos grandes) los paquetes especializados (QuSpin, ITensor/TeNPy, Qiskit) son la mejor opción.
- Su nicho son muchos operadores de tamaño moderado con estructura de sectores explícita: conjuntos de
  desorden, estadísticas sobre estados y canales aleatorios, operadores extremales bajo restricciones
  espectrales, canales entre sectores con trazas coherentes, todo ello con diferenciación automática y GPU.
- La pregunta matemática original (desigualdades de trazas, contraste de Michelson) es un ejemplo de
  problema extremal en el que el muestreo de Monte Carlo falla y la optimización con restricciones funciona.

---

## Lecturas complementarias

- M. A. Nielsen, I. L. Chuang, *Quantum Computation and Quantum Information* — cap. 2 (postulados),
  cap. 8 (canales, T1/T2).
- M. M. Wilde, *Quantum Information Theory* — canales y estados, en un enfoque cercano a las álgebras de
  operadores.
- A. W. Sandvik, "Computational studies of quantum spin systems", AIP Conf. Proc. 1297 (2010) —
  sectores y diagonalización exacta en la práctica.
- H.-P. Breuer, F. Petruccione, *The Theory of Open Quantum Systems* — la ecuación de Lindblad.
- F. Haake, *Quantum Signatures of Chaos* — matrices aleatorias, inversión temporal, GOE/GUE/GSE.
- S. D. Bartlett, T. Rudolph, R. W. Spekkens, Rev. Mod. Phys. 79, 555 (2007) — superselección y sistemas
  de referencia, accesible para matemáticos.
