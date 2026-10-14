[English](../12_equivariant_learning.md) · **Español** · [Français](../fr/12_equivariant_learning.md) · [Deutsch](../de/12_equivariant_learning.md) · [中文](../zh/12_equivariant_learning.md) · [日本語](../ja/12_equivariant_learning.md) · [Русский](../ru/12_equivariant_learning.md)

# 12. Aprendizaje automático equivariante

Una capa de una red neuronal es **equivariante** si transformar la entrada mediante una simetría y aplicar después
la capa da el mismo resultado que aplicar la capa y transformar después la salida: reetiquetar los vértices de un
grafo reetiqueta sus características de la misma manera. Para capas lineales esto es una relación de conmutación,
de modo que las capas equivariantes forman un conmutante, que por el capítulo 2 es, salvo un cambio de base, el
álgebra $\bigoplus_cM_{m_c}(\mathbb C)$. Este capítulo cuenta sus parámetros, las construye y las entrena con
`torch_vn_algebra.decompose`, y usa la estructura de álgebra para obtener capas unitarias, positivas y espectrales
que son equivariantes por construcción.

**Aprenderá** por qué las aplicaciones lineales equivariantes entre dos representaciones son matrices por bloques
cuyos tamaños vienen dados por las multiplicidades; cómo contarlas mediante órbitas o mediante caracteres (DeepSets y
los números de Bell de las redes de grafos); por qué promediar sobre el grupo es la proyección ortogonal sobre las
aplicaciones equivariantes; cómo el cálculo funcional y las parametrizaciones unitarias actúan bloque a bloque sobre
el conmutante; y cómo se cuentan los vectores invariantes, las formas cuadráticas invariantes y los polinomios
invariantes de cualquier grado (fórmula de Molien).

Las capas se construyen y se entrenan con la biblioteca en el [Tutorial 12](../../tutorial/es/12_equivariant_learning.md), que también comprueba numéricamente los resultados de este capítulo.

En todo el capítulo, $G$ es un grupo finito, y una **representación** de $G$ sobre $\mathbb C^n$ es un homomorfismo
$\pi:G\to M_n(\mathbb C)$ con valores unitarios (no es una restricción: las $\pi(g)$ invertibles son unitarias para el
producto escalar $\sum_g\langle\pi(g)x,\pi(g)y\rangle$, sección 4.1). Su **carácter** es $\chi_\pi(g)=\mathrm{Tr}\,\pi(g)$.
Las representaciones irreducibles y equivalentes y la **multiplicidad** $m_\rho(\pi)$ de una $\rho$ irreducible en
$\pi$ ($0$ si $\rho$ no aparece) son como en la sección 2.4, que demuestra que $m_\rho(\pi)$ no depende de la
descomposición y demuestra el lema de Schur; el [capítulo 4](04_symmetries.md) trata ambas cosas con detalle
(secciones 4.1–4.2). Sobre matrices rectangulares, $\langle S,T\rangle=\mathrm{Tr}(S^*T)$ es el producto escalar
de **Hilbert–Schmidt**.

<a id="s-12-1"></a>

## 12.1 Aplicaciones lineales equivariantes

**Definición.** Sean $\pi_1,\pi_2$ representaciones de $G$ sobre $\mathbb C^{n_1}$, $\mathbb C^{n_2}$. Una matriz
$T\in M_{n_2\times n_1}(\mathbb C)$ es **equivariante** (un **operador de entrelazamiento**) si $T\pi_1(g)=\pi_2(g)T$
para todo $g\in G$. Estas matrices forman un espacio vectorial $\mathrm{Hom}_G(\pi_1,\pi_2)$; para
$\pi_1=\pi_2=\pi$ es el conmutante, $\mathrm{Hom}_G(\pi,\pi)=\pi(G)'$.

Las representaciones más comunes en aprendizaje automático permutan coordenadas. Supongamos que $G$ actúa sobre un
conjunto finito $\Omega$; la **representación de permutación** sobre $\mathbb C^\Omega$ es $\pi(g)e_x=e_{gx}$.
Ejemplos: $S_n$ permutando los $n$ puntos de un conjunto ($\Omega=[n]=\{1,\dots,n\}$), y $S_n$ reetiquetando los
vértices de un grafo, actuando sobre $\Omega=[n]^2$, es decir, sobre la matriz de adyacencia mediante
$X\mapsto P_\sigma XP_\sigma^T$.

<a id="r-12-1"></a>

**Lema 12.1 (órbitas).** Sean $\pi_i$ las representaciones de permutación de $G$ sobre $\mathbb C^{\Omega_i}$.
Entonces $T\in\mathrm{Hom}_G(\pi_1,\pi_2)$ si y solo si $T_{gy,gx}=T_{yx}$ para todos $g$, $y\in\Omega_2$,
$x\in\Omega_1$. Las matrices indicadoras de las órbitas de $G$ sobre $\Omega_2\times\Omega_1$ forman una base, y
$\dim\mathrm{Hom}_G(\pi_1,\pi_2)$ es el número de esas órbitas.

*Demostración.* $\pi_2(g)T\pi_1(g)^*e_x=\pi_2(g)Te_{g^{-1}x}=\sum_yT_{y,g^{-1}x}e_{gy}$, luego la entrada $(gy,x)$ de
$\pi_2(g)T\pi_1(g)^*$ es $T_{y,g^{-1}x}$. Por tanto, $T$ es equivariante, es decir, $\pi_2(g)T\pi_1(g)^*=T$ para todo
$g$, si y solo si $T_{gy,x}=T_{y,g^{-1}x}$ para todos $g,y,x$, lo cual (sustituyendo $x$ por $gx$) es la condición
del enunciado: $T$ es constante en las órbitas. Las matrices indicadoras de las órbitas tienen soportes disjuntos,
luego son linealmente independientes, y toda función constante en las órbitas es combinación de ellas. $\square$

*Cálculo:* el lema 12.1 para $S_4$ sobre las características de vértices y de aristas, contando órbitas — [tutorial, sección 12.1](../../tutorial/es/12_equivariant_learning.md#s-12-1).

<a id="r-12-2"></a>

**Corolario 12.2.** (a) (*DeepSets.*) Para $S_n$ sobre $\mathbb C^n$, $n\ge2$, el conmutante es
$\mathrm{span}\{1,J\}$, con $J$ la matriz cuyas entradas son todas iguales a uno, y
$\mathbb C^n=\mathbb C\mathbf 1\oplus\mathbf 1^\perp$ es suma de dos representaciones irreducibles no equivalentes, la
trivial y la **estándar**, de dimensión $n-1$. (b) (*Redes de grafos.*) Las aplicaciones equivariantes de
$\mathbb C^{[n]^p}$ en $\mathbb C^{[n]^q}$ (tensores de orden $p$ y $q$ sobre $n$ vértices) se cuentan mediante las
particiones de un conjunto de $p+q$ elementos en a lo sumo $n$ bloques; para $n\ge p+q$ su número es el
**número de Bell** $B_{p+q}$, el número de todas las particiones de un conjunto de $p+q$ elementos ($B_2=2$,
$B_3=5$, $B_4=15$).

*Demostración.* (b) Dos tuplas $x,y\in[n]^{p+q}$ están en la misma órbita si y solo si tienen el mismo **patrón de
igualdades**, la partición de las posiciones $\{1,\dots,p+q\}$ dada por $i\sim j\Leftrightarrow x_i=x_j$. El patrón
es claramente invariante. Recíprocamente, si $x,y$ tienen el mismo patrón, $x_i\mapsto y_i$ es una inyección bien
definida de $\{x_i\}$ en $[n]$, que se extiende a una permutación $\sigma$ con $\sigma x=y$. Toda partición con a lo
sumo $n$ bloques es un patrón (asígnense valores distintos a los bloques). El lema 12.1 permite concluir. Para
$n\ge p+q$ la cota sobre el número de bloques es vacía. Contando según el número de bloques: $B_2=1+1$;
$B_3=1+3+1$ (una partición de un conjunto de $3$ elementos en dos bloques está determinada por su bloque unitario);
$B_4=1+7+6+1$ (dos bloques: $4$ de tamaños $3+1$ y $3$ de tamaños $2+2$; tres bloques: las $\binom42=6$ elecciones
del par). (a) Para $p=q=1$ hay dos órbitas sobre $[n]^2$, la diagonal y su complemento, con indicadoras $1$ y $J-1$. Así,
$\dim\pi(S_n)'=2=\sum_cm_c^2$ (teorema 2.6, proposición 2.7), lo que obliga a que haya dos sectores con $m_c=1$, es
decir, dos sumandos irreducibles no equivalentes, y $\pi(S_n)'$ está generado por sus proyecciones $P_1,P_2$. La
proyección sobre la recta invariante $\mathbb C\mathbf 1$ está en $\pi(S_n)'$ (lema 2.1), luego es $P_1$ o $P_2$: la
recta, que porta la representación trivial, es un sector, y $\mathbf 1^\perp$ es el otro. $\square$

La capa DeepSets $x\mapsto\lambda x+\gamma\,(\mathbf 1^Tx)\mathbf 1$ (Zaheer et al. 2017) es (a); las bases de
Maron et al. (2019) para redes de grafos son (b).

*Cálculo:* el corolario 12.2 para $n=4$ ($B_2=2$, $B_3=5$, $B_4=15$), y la descomposición de las características de vértices y de aristas de $S_4$ mediante `dc.from_group` (características de aristas: sectores de dimensiones $k=(1,2,3,3)$ con multiplicidades $m=(2,1,1,3)$) — [tutorial, sección 12.1](../../tutorial/es/12_equivariant_learning.md#s-12-1).

Para dos representaciones, todo se reduce a un solo conmutante:

<a id="r-12-3"></a>

**Lema 12.3 (suma directa).** Sea $\pi=\pi_1\oplus\pi_2$ sobre $\mathbb C^{n_1}\oplus\mathbb C^{n_2}$, con las
inclusiones isométricas $Q_i$ y las proyecciones $P_i=Q_iQ_i^*$. Entonces $P_1,P_2\in\pi(G)'$,
$\mathrm{Hom}_G(\pi_1,\pi_2)=Q_2^*\,\pi(G)'\,Q_1$, y $T\mapsto Q_2TQ_1^*$ es una biyección lineal de
$\mathrm{Hom}_G(\pi_1,\pi_2)$ sobre $P_2\,\pi(G)'P_1$.

*Demostración.* $\pi(g)Q_i=Q_i\pi_i(g)$ y $Q_i^*\pi(g)=\pi_i(g)Q_i^*$, porque $\pi(g)$ es diagonal por bloques; por
tanto, $P_i\pi(g)=Q_i\pi_i(g)Q_i^*=\pi(g)P_i$. Si $X\in\pi(G)'$, entonces
$Q_2^*XQ_1\pi_1(g)=Q_2^*X\pi(g)Q_1=Q_2^*\pi(g)XQ_1=\pi_2(g)Q_2^*XQ_1$. Si $T$ es equivariante, entonces
$\pi(g)Q_2TQ_1^*=Q_2\pi_2(g)TQ_1^*=Q_2T\pi_1(g)Q_1^*=Q_2TQ_1^*\pi(g)$, luego $Q_2TQ_1^*\in\pi(G)'$, es igual a
$P_2(Q_2TQ_1^*)P_1$, y $Q_2^*(Q_2TQ_1^*)Q_1=T$. Recíprocamente, todo $Y=P_2XP_1$ con $X\in\pi(G)'$ está en
$\pi(G)'$ y $Y=Q_2(Q_2^*YQ_1)Q_1^*$. $\square$

*Cálculo:* las esquinas $Q_2^*\,\pi(G)'\,Q_1$ para las características de vértices y de aristas — [tutorial, sección 12.2](../../tutorial/es/12_equivariant_learning.md#s-12-2).

<a id="r-12-4"></a>

**Teorema 12.4 (las aplicaciones equivariantes son matrices por bloques).** Sean $\rho_1,\dots,\rho_C$ las
representaciones irreducibles, dos a dos no equivalentes, que aparecen en $\pi_1\oplus\pi_2$, $k_c=\dim\rho_c$,
$a_c=m_{\rho_c}(\pi_1)$, $b_c=m_{\rho_c}(\pi_2)$. Existen unitarias
$V_1:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{a_c}\to\mathbb C^{n_1}$ y
$V_2:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{b_c}\to\mathbb C^{n_2}$ con
$V_i^*\pi_i(g)V_i=\bigoplus_c\rho_c(g)\otimes1$, tales que
$$
\mathrm{Hom}_G(\pi_1,\pi_2)=\Big\{V_2\Big(\bigoplus_{c=1}^C1_{k_c}\otimes y_c\Big)V_1^*:\ y_c\in M_{b_c\times a_c}(\mathbb C)\Big\},
$$
y $(y_c)_c\mapsto T$ es inyectiva. Por tanto, $\mathrm{Hom}_G(\pi_1,\pi_2)\cong\bigoplus_cM_{b_c\times a_c}(\mathbb C)$,
$$
\dim\mathrm{Hom}_G(\pi_1,\pi_2)=\sum_\rho m_\rho(\pi_1)\,m_\rho(\pi_2),
$$
donde la suma recorre todas las representaciones irreducibles salvo equivalencia; para $\pi_1=\pi_2$ esto es
$\pi(G)'\cong\bigoplus_cM_{m_c}$ (proposición 2.7), y $\mathrm{Hom}_G(\pi_1,\pi_2)=0$ si y solo si ninguna
representación irreducible aparece en ambas.

*Demostración.* Apliquemos el teorema 2.6 y la proposición 2.7 a $\pi=\pi_1\oplus\pi_2$: existe una unitaria $V$ con
$V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$ y $V^*\pi(G)'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$. Por la unicidad
de las multiplicidades (sección 2.4), toda representación irreducible que aparece en $\pi$ es equivalente a
exactamente una $\rho_c$, de modo que las $\rho_c$ son representantes de las clases del enunciado. Cualesquiera
otros representantes son $\rho_c'=w_c\rho_cw_c^*$ con unitarias $w_c$ (lema 4.2(c)), y sustituir $V$ por
$V\bigoplus_c(w_c^*\otimes1)$ convierte la primera forma normal en $\bigoplus_c\rho'_c(g)\otimes1_{m_c}$ y deja la
segunda sin cambios, pues $w_c^*\otimes1$ conmuta con $1\otimes M_{m_c}$; así, podemos usar las $\rho_c$ de la
proposición 2.7. Como
$X\mapsto V^*XV$ es un $*$-isomorfismo y $P_1\in\pi(G)'$ (lema 12.3), $V^*P_1V=\bigoplus_c1\otimes p_c$ con
proyecciones $p_c$. Elijamos unitarias $u_c$ con $u_c^*p_cu_c=1_{a_c}\oplus0_{b_c}$ (teorema 1.3), donde ahora
$a_c=\mathrm{rank}\,p_c$, $b_c=m_c-a_c$, y sustituyamos $V$ por $V\bigoplus_c(1\otimes u_c)$; esto no cambia
ninguna de las dos formas normales, pues $1\otimes u_c$ conmuta con $\rho_c(g)\otimes1$. Sean
$\iota_c:\mathbb C^{a_c}\to\mathbb C^{m_c}$ y $\kappa_c:\mathbb C^{b_c}\to\mathbb C^{m_c}$ las inclusiones de las
primeras $a_c$ y de las últimas $b_c$ coordenadas, y
$$
V_1=Q_1^*V\bigoplus_c(1\otimes\iota_c),\qquad V_2=Q_2^*V\bigoplus_c(1\otimes\kappa_c).
$$
$W_1=V\bigoplus_c(1\otimes\iota_c)$ es una isometría con imagen $\mathrm{range}\,P_1=\mathrm{range}\,Q_1$, luego
$V_1^*V_1=W_1^*P_1W_1=1$, y $V_1$ es unitaria porque $\sum_ck_ca_c=\mathrm{rank}\,P_1=n_1$. Además,
$\pi_1(g)V_1=Q_1^*\pi(g)W_1=Q_1^*V\big(\bigoplus_c\rho_c(g)\otimes1_{m_c}\big)\bigoplus_c(1\otimes\iota_c)
=V_1\bigoplus_c\rho_c(g)\otimes1_{a_c}$. Así, $\pi_1\cong\bigoplus_c\rho_c\otimes1_{a_c}$ con $\rho_c$ irreducibles
no equivalentes, y $a_c=m_{\rho_c}(\pi_1)$ por la unicidad de las multiplicidades (sección 2.4); análogamente para
$V_2$ y $b_c$. Por último, para $X=V\big(\bigoplus_c1\otimes x_c\big)V^*\in\pi(G)'$,
$P_2XP_1=V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*$, donde $y_c=\kappa_c^*x_c\iota_c$ es la esquina
inferior izquierda de tamaño $b_c\times a_c$ de $x_c$, una matriz arbitraria. Por el lema 12.3,
$\mathrm{Hom}_G(\pi_1,\pi_2)=Q_2^*P_2\,\pi(G)'P_1Q_1$ consta de las matrices
$Q_2^*V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*Q_1=V_2\big(\bigoplus_c1\otimes y_c\big)V_1^*$, y
$(y_c)\mapsto T$ es inyectiva porque $V_1,V_2$ son unitarias. Una $\rho$ irreducible que no aparece en
$\pi_1\oplus\pi_2$ tiene $m_\rho(\pi_1)m_\rho(\pi_2)=0$. $\square$

*En la biblioteca:* la demostración del teorema 12.4 como algoritmo, las capas `EquivariantLinear` y `EquivariantMap` construidas sobre `dc.from_group`, `D.commutant_algebra()` y `D.from_commutant`; de las características de vértices a las de aristas de $S_4$, $\dim\mathrm{Hom}_G=1\cdot2+1\cdot3=5=B_3$ — [tutorial, sección 12.2](../../tutorial/es/12_equivariant_learning.md#s-12-2).

> **Significado físico** (interpretación física, no un teorema de este curso). El teorema 12.4 es el núcleo de
> álgebra lineal de las **reglas de selección** y el caso escalar del teorema de Wigner–Eckart: un operador que conmuta con la simetría solo puede conectar estados del mismo
> tipo de simetría, y dentro de un tipo actúa como $1_{k_c}\otimes y_c$, es decir, mediante un «elemento de matriz
> reducido» $y_c$ que es ciego a la posición dentro del multiplete. (Los operadores que se transforman ellos mismos
> según una representación no trivial, como los vectoriales, necesitan además los coeficientes de Clebsch–Gordan.)
> Un hamiltoniano simétrico es el caso $\pi_1=\pi_2$; sus degeneraciones $k_c$ son las dimensiones de las
> representaciones irreducibles.

<a id="s-12-2"></a>

## 12.2 Una capa equivariante con la biblioteca

*En la biblioteca:* capas equivariantes parametrizadas por los bloques del teorema 12.4 — [tutorial, sección 12.2](../../tutorial/es/12_equivariant_learning.md#s-12-2).

La segunda manera de hacer equivariante una capa es promediarla sobre el grupo. El **twirl** es
$$
\mathcal P(W)=\frac1{|G|}\sum_{g\in G}\pi_2(g)\,W\,\pi_1(g)^*,\qquad W\in M_{n_2\times n_1}(\mathbb C).
$$

<a id="r-12-5"></a>

**Teorema 12.5 (el twirl es la proyección ortogonal).** $\mathcal P$ es la proyección ortogonal, para el producto
escalar de Hilbert–Schmidt, de $M_{n_2\times n_1}(\mathbb C)$ sobre $\mathrm{Hom}_G(\pi_1,\pi_2)$, y
$$
\dim\mathrm{Hom}_G(\pi_1,\pi_2)=\frac1{|G|}\sum_{g\in G}\overline{\chi_{\pi_1}(g)}\,\chi_{\pi_2}(g).
$$

*Demostración.* Para $h\in G$, $\pi_2(h)\mathcal P(W)\pi_1(h)^*=\frac1{|G|}\sum_g\pi_2(hg)W\pi_1(hg)^*=\mathcal P(W)$,
pues $g\mapsto hg$ es una biyección de $G$; así, $\mathcal P(W)$ es equivariante. Si $W$ es equivariante, cada
término es igual a $W$, luego $\mathcal P(W)=W$. Por tanto, $\mathcal P^2=\mathcal P$ con imagen
$\mathrm{Hom}_G(\pi_1,\pi_2)$. Es autoadjunta:
$\langle S,\mathcal P(W)\rangle=\frac1{|G|}\sum_g\mathrm{Tr}\big((\pi_2(g)^*S\pi_1(g))^*W\big)$
$=\langle\mathcal P(S),W\rangle$, donde hemos sustituido $g\mapsto g^{-1}$ y usado $\pi_i(g^{-1})=\pi_i(g)^*$. Un
idempotente autoadjunto es la proyección ortogonal sobre su imagen. Su traza es su rango. La aplicación lineal
$W\mapsto AWB$ sobre $M_{n_2\times n_1}$ tiene traza $\mathrm{Tr}\,A\mathrm{Tr}\,B$ (el coeficiente de la
unidad matricial $E_{ij}$ en $AE_{ij}B$ es $A_{ii}B_{jj}$); con $A=\pi_2(g)$, $B=\pi_1(g)^*$ esto da la
fórmula. $\square$

Junto con el teorema 12.4,
$\sum_\rho m_\rho(\pi_1)m_\rho(\pi_2)=|G|^{-1}\sum_g\overline{\chi_{\pi_1}(g)}\chi_{\pi_2}(g)$; para $\pi_1=\pi_2$ esta
es la fórmula de caracteres para $\dim\pi(G)'$ de la sección 4.2. La
parametrización $W\mapsto\mathcal P(W)$ alcanza toda aplicación equivariante, $\mathcal P(W)$ es la aplicación
equivariante más próxima a $W$, y, como $\mathcal P$ es autoadjunta, el gradiente respecto de $W$ es $\mathcal P$
aplicada al gradiente respecto de $T=\mathcal P(W)$. Pero $n_1n_2-\dim\mathrm{Hom}_G$ de sus parámetros se
desperdician en el núcleo.

*Cálculo:* el teorema 12.5 sobre las características de aristas y la fórmula de dimensión de las características de vértices a las de aristas — [tutorial, sección 12.2](../../tutorial/es/12_equivariant_learning.md#s-12-2). *En la biblioteca:* el twirl `dc.twirl`, la misma proyección calculada a partir de la descomposición, `D.conditional_expectation(W, 'commutant')`, y la fórmula de caracteres `dc.commutant_dimension` — [tutorial, sección 12.2](../../tutorial/es/12_equivariant_learning.md#s-12-2).

Las capas anteriores son lineales. Una red las alterna con no linealidades, y estas también han de ser
equivariantes. Para una representación de permutación, una función aplicada a cada coordenada es equivariante, ya
que $(\pi(g)x)_{gy}=x_y$; para otras representaciones no lo es, y se usan no linealidades que actúan sobre normas o
sobre los coeficientes en una base que $G$ permuta (ejercicio 3).

> **Significado físico** (interpretación física, no un teorema de este curso). El twirl es la esperanza condicional sobre el conmutante (capítulo 2, ejercicio 3):
> promediar un acoplamiento sobre todas las operaciones de simetría conserva exactamente su parte simétrica.
> «Simetrizar» un hamiltoniano o una interacción sumando sobre el grupo de simetría es la misma proyección.

<a id="s-12-3"></a>

## 12.3 Capas restringidas y espectrales

Los pesos $w=(w_c)_c$ de una capa equivariante $x\mapsto Tx$, $T\in\pi(G)'$, viven en la $*$-álgebra
$\bigoplus_cM_{m_c}$ (teorema 12.4), de modo que toda construcción de los capítulos 1 y 11 se les aplica. Escribamos
$\Phi(w)=V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*$ para $w\in\bigoplus_cM_{m_c}$.

<a id="r-12-6"></a>

**Proposición 12.6.** (a) $\Phi$ es un $*$-isomorfismo de $\bigoplus_cM_{m_c}$ sobre $\pi(G)'$, y $\Phi(w)$ es
autoadjunto, unitario o positivo si y solo si lo es cada $w_c$. (b) Para $w$ normal y una función $f$ sobre
$\bigcup_c\mathrm{spec}\,w_c$, $f(\Phi(w))=\Phi\big((f(w_c))_c\big)$; en particular, $f(X)\in\pi(G)'$ para todo
$X\in\pi(G)'$ normal. (c) Las unitarias equivariantes son las $\Phi(u_0e^{s-s^*})$, con $s$ recorriendo
$\bigoplus_cM_{m_c}$, para cualquier unitaria fija $u_0\in\bigoplus_cM_{m_c}$; en particular, son las exponenciales
$e^S$ de los $S\in\pi(G)'$ antiautoadjuntos.

*Demostración.* (a) $\Phi$ es lineal, unital, multiplicativa y conserva la involución porque
$(1\otimes x)(1\otimes y)=1\otimes xy$ y $(1\otimes x)^*=1\otimes x^*$; es inyectiva, y sobreyectiva sobre $\pi(G)'$
por el teorema 2.6. Por tanto, $\Phi(w)^*=\Phi(w)$ si y solo si $w^*=w$, y $\Phi(w)^*\Phi(w)=1=\Phi(w)\Phi(w)^*$ si y
solo si $w^*w=1=ww^*$. Si todas las $w_c=t_c^*t_c$, entonces $\Phi(w)=\Phi(t)^*\Phi(t)\ge0$ (corolario 1.7).
Recíprocamente, si $\Phi(w)\ge0$, entonces para $v\in\mathbb C^{m_c}$ y el vector $\xi=V(e_1\otimes v)$ del sector
$c$, $0\le\langle\xi,\Phi(w)\xi\rangle=\langle v,w_cv\rangle$.
(b) Sea $w_c=u_c\Lambda_cu_c^*$ (teorema 1.3). Entonces $\Phi(w)=U\big(\bigoplus_c1\otimes\Lambda_c\big)U^*$ con la
unitaria $U=V\bigoplus_c(1\otimes u_c)$, una diagonalización con los valores propios de $w_c$ repetidos $k_c$ veces;
por definición y por el teorema 1.5(a), $f(\Phi(w))=U\big(\bigoplus_c1\otimes f(\Lambda_c)\big)U^*=\Phi\big((f(w_c))_c\big)$.
(c) Por (a), las unitarias equivariantes son las $\Phi(u)$ con $u$ unitaria en $\bigoplus_cM_{m_c}$, y estas son las
$u_0e^{s-s^*}$ por el lema 11.1. Con $u_0=1$, (b) da $\Phi(e^{s-s^*})=e^{\Phi(s-s^*)}$, y $\Phi(s-s^*)$ recorre los
elementos antiautoadjuntos de $\pi(G)'$. $\square$

*En la biblioteca:* capas equivariantes unitarias y positivas con `optimize.UnitaryParam` y `PositiveParam` sobre `D.commutant_algebra()`, capas espectrales con `apply_function`; el entrenamiento de una capa unitaria equivariante, y (b) para $e^{-tL}$ con su derivada respecto de $t$ — [tutorial, sección 12.3](../../tutorial/es/12_equivariant_learning.md#s-12-3).

> **Significado físico** (interpretación física, no un teorema de este curso). Un hamiltoniano simétrico
> $H\in\pi(G)'$ genera una evolución temporal simétrica $e^{-itH}\in\pi(G)'$ (capítulo 8), que nunca mezcla
> sectores y actúa solo sobre los espacios de multiplicidad. Una capa unitaria equivariante es una evolución de este
> tipo, aprendida en lugar de dada.

<a id="s-12-4"></a>

## 12.4 Invariantes

Un vector $x$ es **invariante** si $\pi(g)x=x$ para todo $g$; una función $q$ sobre $\mathbb C^n$ es invariante si
$q(\pi(g)x)=q(x)$. Las características invariantes son lo que puede devolver una red invariante (por ejemplo, un
clasificador de grafos).

<a id="r-12-7"></a>

**Proposición 12.7.** (a) Los vectores invariantes son la imagen de la proyección ortogonal
$P_0=|G|^{-1}\sum_g\pi(g)$; su dimensión es $m_{\mathrm{triv}}(\pi)=|G|^{-1}\sum_g\chi_\pi(g)$, y si la representación
trivial aparece, como $\rho_{c_0}$, entonces $P_0$ es la proyección central minimal $P_{c_0}$ de $\pi(G)''$.
(b) Una forma hermítica $q_A(x)=\langle x,Ax\rangle$, $A=A^*$, es invariante si y solo si $A\in\pi(G)'$; las formas
hermíticas invariantes forman un espacio vectorial real de dimensión $\sum_cm_c^2$.
(c) Los polinomios invariantes homogéneos de grado $2$ en las coordenadas $x_1,\dots,x_n$ forman un espacio de
dimensión
$$
\frac1{|G|}\sum_{g\in G}\frac{\chi_\pi(g)^2+\chi_\pi(g^2)}2 .
$$
Para una representación de permutación, es el número de órbitas de $G$ sobre los multiconjuntos $\{x,y\}$ de dos
elementos de $\Omega$.

*Demostración.* (a) Los vectores invariantes son $\mathrm{Hom}_G(\mathrm{triv},\pi)$, con $\mathrm{triv}$ la
representación trivial sobre $\mathbb C=M_{1\times1}$; para ella el twirl es $\mathcal P(x)=P_0x$, de modo que el
teorema 12.5 da la proyección y la dimensión, y el teorema 12.4 da $m_{\mathrm{triv}}(\pi)$. En la base del
teorema 2.6, escribamos un vector invariante como $\xi=\sum_c\xi_c$ con
$\xi_c=\sum_jv_{c,j}\otimes e_j\in\mathbb C^{k_c}\otimes\mathbb C^{m_c}$. Como $\pi(g)$ conserva los sectores,
$\rho_c(g)v_{c,j}=v_{c,j}$ para todo $g$. Un $v_{c,j}$ no nulo genera una recta invariante que porta la
representación trivial, de modo que entonces $\rho_c$, por ser irreducible, es trivial y $c=c_0$. Por tanto, los
vectores invariantes son los vectores del sector $c_0$, donde $\pi(g)$ actúa como $1$: la imagen de $P_{c_0}$, que
es, en consecuencia, $P_0$.
(b) $q_A(\pi(g)x)=\langle x,\pi(g)^*A\pi(g)x\rangle$, y una matriz $D$ con $\langle x,Dx\rangle=0$ para todo $x$ es
nula (polarización, demostración del corolario 1.4). Así, $q_A$ es invariante si y solo si $\pi(g)^*A\pi(g)=A$ para
todo $g$, es decir, $A\in\pi(G)'$. Por la proposición 12.6(a), los elementos autoadjuntos de $\pi(G)'$ son los
$\Phi(w)$ con $w_c$ hermíticas, un espacio real de dimensión $\sum_cm_c^2$.
(c) Un polinomio cuadrático $p(x)=x^TBx$ está determinado por su matriz simétrica $B$. Pongamos
$(g\cdot p)(x)=p(\pi(g)^{-1}x)$; entonces $g\cdot(h\cdot p)=(gh)\cdot p$, y $g\cdot p$ tiene la matriz
$\sigma(g)(B)=A^TBA$, $A=\pi(g)^{-1}$. El promedio $M=|G|^{-1}\sum_g\sigma(g)$ cumple $\sigma(h)M=M$ y $MB=B$ para $B$
invariante, de modo que (como en la demostración del teorema 12.5) es idempotente con imagen los polinomios
invariantes, y la dimensión de estos es $\mathrm{Tr}\,M=|G|^{-1}\sum_g\mathrm{Tr}\,\sigma(g)$ (un idempotente
es $1\oplus0$ en una base adaptada a imagen $\oplus$ núcleo, luego su traza es su rango). $A$ es unitaria,
luego $A=R\Lambda R^{-1}$ con $\Lambda=\mathrm{diag}(\lambda)$ (teorema 1.3). Las matrices simétricas
$B_{ij}=(R^{-1})^T(e_ie_j^T+e_je_i^T)R^{-1}$, $i\le j$, forman una base, y
$A^TB_{ij}A=(R^{-1})^T\Lambda(e_ie_j^T+e_je_i^T)\Lambda R^{-1}=\lambda_i\lambda_jB_{ij}$. Así,
$\mathrm{Tr}\,\sigma(g)=\sum_{i\le j}\lambda_i\lambda_j=\frac12\big((\mathrm{Tr}\,A)^2+\mathrm{Tr}\,A^2\big)
=\frac12\big(\chi_\pi(g^{-1})^2+\chi_\pi(g^{-2})\big)$, y sumar sobre $g^{-1}$ en lugar de $g$ da la fórmula. Para
una representación de permutación, $G$ permuta los monomios $x_ax_b$, de modo que un polinomio es invariante si y
solo si sus coeficientes son constantes en las órbitas de monomios. $\square$

*Cálculo:* (a), (b) y los dos recuentos de (c) para las características de aristas de $S_4$ — [tutorial, sección 12.4](../../tutorial/es/12_equivariant_learning.md#s-12-4).

Para grafos, las características de aristas invariantes están generadas por las indicadoras de la diagonal y de su
complemento ($B_2=2$), y las lecturas lineales invariantes son la suma de las entradas diagonales y la de las no
diagonales. Sobre las características de vértices $\mathbb C^4$, (c) cuenta las órbitas de $S_4$ sobre los
multiconjuntos de dos vértices: hay dos.

Sobre $\mathbb C^4$ son $\sum_ix_i^2$ y $\sum_{i<j}x_ix_j$. La demostración de (c) funciona en todos los grados:

<a id="r-12-8"></a>

**Proposición 12.8 (fórmula de Molien).** Sea $\mathbb C[x]_d$ el espacio de los polinomios homogéneos de grado $d$
en $x_1,\dots,x_n$ y $\mathbb C[x]^G_d$ el subespacio de los invariantes, $p(\pi(g)x)=p(x)$ para todo $g$. Para
$|t|<1$,
$$
\sum_{d\ge0}\dim\mathbb C[x]^G_d\,t^d=\frac1{|G|}\sum_{g\in G}\frac1{\det\big(1-t\pi(g)\big)} .
$$
Para una representación de permutación, $\dim\mathbb C[x]^G_d$ es el número de órbitas de $G$ sobre los monomios de
grado $d$ (los multiconjuntos de $d$ elementos de $\Omega$). El coeficiente de $t^2$ es la proposición 12.7(c).

*Demostración.* Para una matriz invertible $M$, sea $\phi_M(p)(x)=p(Mx)$ (sustitución). $\phi_M$ es lineal, aplica
$\mathbb C[x]_d$ en sí mismo, y $\phi_N\circ\phi_M=\phi_{MN}$, luego $\phi_M$ es biyectiva con inversa
$\phi_{M^{-1}}$. Por tanto, $\sigma(g)=\phi_{\pi(g)^{-1}}$ cumple $\sigma(g)\sigma(h)=\sigma(gh)$, y $p$ es invariante si
y solo si $\sigma(g)p=p$ para todo $g$. Como en la demostración de la proposición 12.7(c),
$M_d=|G|^{-1}\sum_g\sigma(g)$ sobre $\mathbb C[x]_d$ es idempotente con imagen $\mathbb C[x]^G_d$, luego
$\dim\mathbb C[x]^G_d=\mathrm{Tr}\,M_d=|G|^{-1}\sum_g\mathrm{Tr}\,\sigma_d(g)$, con $\sigma_d(g)$ la restricción de
$\sigma(g)$ a $\mathbb C[x]_d$. Fijemos $g$ y escribamos la unitaria $A=\pi(g)^{-1}$ como $A=R\Lambda R^{-1}$,
$\Lambda=\mathrm{diag}(\lambda_1,\dots,\lambda_n)$, $|\lambda_i|=1$ (teorema 1.3). Con las formas lineales
$\ell_i(x)=(R^{-1}x)_i$, los productos $\ell^\alpha=\prod_i\ell_i^{\alpha_i}=\phi_{R^{-1}}(x^\alpha)$, $|\alpha|=d$, forman
una base de $\mathbb C[x]_d$ (la imagen de la base de monomios por la biyección $\phi_{R^{-1}}$), y como
$R^{-1}A=\Lambda R^{-1}$,
$$
(\sigma(g)\ell^\alpha)(x)=\prod_i\big((R^{-1}Ax)_i\big)^{\alpha_i}=\prod_i\big(\lambda_i\ell_i(x)\big)^{\alpha_i}=\lambda^\alpha\ell^\alpha(x).
$$
Así, $\mathrm{Tr}\,\sigma_d(g)=\sum_{|\alpha|=d}\lambda^\alpha$. Para $|t|<1$, el producto de las $n$ series
absolutamente convergentes $(1-t\lambda_i)^{-1}=\sum_{a\ge0}t^a\lambda_i^a$ puede desarrollarse y reagruparse por
grados:
$$
\frac1{\det(1-tA)}=\prod_{i=1}^n\frac1{1-t\lambda_i}=\sum_\alpha t^{|\alpha|}\lambda^\alpha=\sum_{d\ge0}t^d\mathrm{Tr}\,\sigma_d(g).
$$
Promediar sobre $g$ y sustituir $g$ por $g^{-1}$ (de modo que $A=\pi(g)$) da la fórmula; ambos miembros convergen,
pues $\dim\mathbb C[x]^G_d\le\dim\mathbb C[x]_d=\binom{n+d-1}d$. Para una representación de permutación,
$(\pi(g)^{-1}x)_y=x_{gy}$, de modo que $\sigma(g)$ aplica el monomio $\prod_jx_{y_j}$ en el monomio $\prod_jx_{gy_j}$:
$G$ permuta los monomios de grado $d$, un polinomio es invariante si y solo si sus coeficientes son constantes en las
órbitas, y las sumas sobre las órbitas forman una base de $\mathbb C[x]^G_d$. $\square$

(Es el teorema de Molien de 1897; para profundizar, véase B. Sturmfels, *Algorithms in Invariant Theory*, Springer
2008, sección 2.2.) Para $S_4$ sobre $\mathbb C^4$, los recuentos $1,1,2,3,5$ para $d\le4$ son las particiones de $d$
en a lo sumo cuatro partes (las órbitas de multiconjuntos están determinadas por sus multiplicidades, ordenadas).

*Cálculo:* la serie de Molien hasta $t^4$ frente a las órbitas de monomios, para las características de vértices y de aristas — [tutorial, sección 12.4](../../tutorial/es/12_equivariant_learning.md#s-12-4).

> **Significado físico** (interpretación física, no un teorema de este curso). En la teoría de Landau de las
> transiciones de fase, la energía libre es un polinomio invariante en un parámetro de orden real $x$, sobre el que
> $G$ actúa mediante matrices ortogonales reales. Su parte cuadrática es $x^TBx$ con $B$ en el conmutante (pues
> $\pi(g)^T=\pi(g)^{-1}$), de modo que, por la proposición 12.6(a), es diagonal por bloques sobre los sectores. La
> teoría de Landau postula entonces que la fase simétrica se vuelve inestable en el sector en el que un valor propio
> de $B$ cambia de signo en primer lugar, y que, genéricamente, la representación irreducible de ese sector determina
> la simetría de la fase ordenada.

<a id="s-12-5"></a>

## 12.5 Con qué comparar

*En la biblioteca:* comparación con otro software (escnn, e3nn, EMLP) y lo que aporta la biblioteca — [tutorial, sección 12.5](../../tutorial/es/12_equivariant_learning.md#s-12-5).

<a id="x-12"></a>

## Ejercicios

1. (**DeepSets con canales.**) $S_n$ ($n\ge2$) actúa sobre $\mathbb C^n\otimes\mathbb C^d$ mediante
   $P_\sigma\otimes1_d$. Demuestre que las aplicaciones equivariantes de $\mathbb C^n\otimes\mathbb C^{d_1}$ en
   $\mathbb C^n\otimes\mathbb C^{d_2}$ son las $1\otimes A+J\otimes B$ con $A,B\in M_{d_2\times d_1}$.
2. (**Capas reales.**) Supongamos que $\pi_1,\pi_2$ tienen matrices reales (por ejemplo, representaciones de
   permutación). Demuestre que las matrices equivariantes reales forman un espacio vectorial real de dimensión
   $\dim_{\mathbb C}\mathrm{Hom}_G(\pi_1,\pi_2)$, y que $T\mapsto\mathrm{Re}\,T$ aplica
   $\mathrm{Hom}_G(\pi_1,\pi_2)$ sobre dicho espacio.
3. (**No linealidades.**) (a) Demuestre que aplicar una función $f:\mathbb C\to\mathbb C$ a cada coordenada conmuta
   con toda representación de permutación. (b) Demuestre que no conmuta con la rotación de ángulo $\pi/2$ de
   $\mathbb R^2$ para $f(t)=\max(t,0)$. (c) Demuestre que $x\mapsto\varphi(\|x\|)x$ conmuta con toda representación
   unitaria.
4. (**¿Cuántas muestras determinan una capa?**) Sean $T=\Phi(w)\in\pi(G)'$ y $x_1,\dots,x_N\in\mathbb C^n$. Escriba
   $V^*x_i$ en el sector $c$ como $\sum_{a,j}(X_{i,c})_{aj}e_a\otimes e_j$ con $X_{i,c}\in M_{k_c\times m_c}$, y sea
   $S_c\in M_{Nk_c\times m_c}$ la matriz formada apilando $X_{1,c},\dots,X_{N,c}$. Demuestre que las $T$
   equivariantes con $Tx_1=\dots=Tx_N=0$ forman un espacio de dimensión $\sum_cm_c(m_c-\mathrm{rank}\,S_c)$;
   para vectores $x_i$ gaussianos complejos estándar independientes esto es $\sum_cm_c\max(0,m_c-Nk_c)$ casi
   seguramente. Deduzca que una $T$ equivariante está determinada por $Tx_1,\dots,Tx_N$ (entre las aplicaciones
   equivariantes) si y solo si esta dimensión es $0$, y que para la representación de aristas dos muestras
   gaussianas determinan una capa equivariante casi seguramente y una sola muestra nunca lo hace, mientras que una
   matriz $16\times16$ arbitraria nunca está determinada por menos de dieciséis muestras.

*Cálculo:* los ejercicios 1, 3 y 4, numéricamente — [tutorial, Ejercicios](../../tutorial/es/12_equivariant_learning.md#x-12).

<details>
<summary>Soluciones</summary>

1. $P_\sigma J=J=JP_\sigma$, luego estas aplicaciones son equivariantes. Son linealmente independientes: si
   $1\otimes A+J\otimes B=0$, los bloques fuera de la diagonal dan $B=0$, y entonces $A=0$. Por el corolario 12.2(a)
   y el teorema 2.6, $P_\sigma\otimes1_d\cong(\mathrm{triv}\otimes1_d)\oplus(\mathrm{std}\otimes1_d)$, de modo que
   ambas representaciones irreducibles tienen multiplicidad $d$, y el teorema 12.4 da
   $\dim\mathrm{Hom}_G=2d_1d_2$, la dimensión del espacio de las $1\otimes A+J\otimes B$.

2. Si $T$ es equivariante, también lo es $\bar T$, ya que $\overline{\pi_i(g)}=\pi_i(g)$. Por tanto,
   $\mathrm{Re}\,T=(T+\bar T)/2$ y $\mathrm{Im}\,T=(T-\bar T)/2i$ son reales y equivariantes, de modo que la
   aplicación $T\mapsto\mathrm{Re}\,T$ es sobreyectiva sobre las matrices equivariantes reales (es la identidad
   sobre ellas), y $\mathrm{Hom}_G(\pi_1,\pi_2)=H_{\mathbb R}\oplus iH_{\mathbb R}$, donde $H_{\mathbb R}$ denota
   las matrices equivariantes reales. Una base real de $H_{\mathbb R}$ es entonces una base compleja de
   $\mathrm{Hom}_G$: lo genera, y una relación compleja $\sum_j(\alpha_j+i\beta_j)T_j=0$ con $T_j$ reales se
   separa en parte real y parte imaginaria. En el lema 12.1, las indicadoras de las órbitas son una base real de este
   tipo. (Los bloques $y_c$ del teorema 12.4 son complejos incluso para datos reales, ya que la $V$ calculada por la
   biblioteca no tiene por qué ser real; parametrizar $\mathrm{Re}\,T$ mediante $y_c$ complejos usa el doble de
   parámetros reales de los necesarios.)

3. (a) $(\pi(g)f(x))_{gy}=f(x)_y=f(x_y)=f((\pi(g)x)_{gy})$. (b) Para $R(x_1,x_2)=(-x_2,x_1)$ y $x=(1,-1)$,
   $f(Rx)=f(1,1)=(1,1)$, pero $Rf(x)=R(1,0)=(0,1)$. (c) $\|\pi(g)x\|=\|x\|$, luego
   $\varphi(\|\pi(g)x\|)\pi(g)x=\pi(g)\varphi(\|x\|)x$.

4. $(1\otimes w_c)(e_a\otimes e_j)=\sum_le_a\otimes(w_c)_{lj}e_l$, luego $V^*Tx_i$ tiene, en el sector $c$, la matriz
   de coeficientes $X_{i,c}w_c^T$. Por tanto, $Tx_i=0$ para todo $i$ si y solo si $S_cw_c^T=0$ para todo $c$, es
   decir, si y solo si cada columna de $w_c^T$ está en $\ker S_c$, un espacio de dimensión
   $m_c-\mathrm{rank}\,S_c$: las $w_c$ con esta propiedad forman un espacio de dimensión
   $m_c(m_c-\mathrm{rank}\,S_c)$, y $\Phi$ es inyectiva. Para $x_i$ gaussianos, las coordenadas de $V^*x_i$ son
   las de $x_i$ en la base ortonormal de las columnas de $V$, de nuevo gaussianas complejas estándar independientes
   (proposición A.16(a)), e independientes en $i$; así, $S_c$ es una matriz de Ginibre $Nk_c\times m_c$, que tiene
   densidad (proposición A.16(b)). Su rango es a lo sumo $\min(Nk_c,m_c)$, con igualdad salvo que la submatriz
   cuadrada superior izquierda de ese tamaño sea singular; $|\det|^2$ de esa submatriz es un polinomio real en las
   partes real e imaginaria de las entradas, no idénticamente nulo (vale $1$ cuando la submatriz es la identidad), de
   modo que por el corolario A.7(b) el rango es $\min(Nk_c,m_c)$ casi seguramente. Para la representación de
   aristas, $(k,m)=(1,2),(2,1),(3,1),(3,3)$ (la descomposición calculada), solo contribuye el sector $(1,2)$:
   $2\cdot(2-1)=2$ para $N=1$, y $0$ para $N=2$. Una matriz $16\times16$ arbitraria es el caso del grupo trivial
   (un sector, $k=1$, $m=16$): la dimensión es $16(16-\mathrm{rank}\,S)\ge16(16-N)>0$ para cualquier elección de
   $N<16$ muestras, pues $\mathrm{rank}\,S\le N$, y $0$ casi seguramente para $N=16$ muestras gaussianas. Como
   $T\mapsto(Tx_1,\dots,Tx_N)$ es lineal, dos aplicaciones equivariantes coinciden en las muestras si y solo si su
   diferencia está en el espacio que se acaba de calcular, lo que da el criterio. Para $N=1$ el espacio tiene
   dimensión al menos $2$ para *todo* $x_1$, porque $\mathrm{rank}\,S_c\le Nk_c=1<m_c=2$ en el sector $(1,2)$.

</details>

**Cómo seguir.** El último capítulo, el [capítulo 13](13_jordan_algebras.md), explica por qué la mecánica cuántica usa *-álgebras y no álgebras de Jordan. El [capítulo 4](04_symmetries.md) desarrolla la teoría de
representaciones usada aquí: el lema de Schur, los caracteres y el algoritmo aleatorizado en que se basa
`dc.from_group` (secciones 4.1–4.3). El [capítulo 9](09_spin_chains.md) aplica los mismos sectores a las cadenas de
espines, y la medida de Haar sobre $\prod_cU(m_c)$ del [capítulo 5](05_random_matrices.md) proporciona la
inicialización aleatoria natural de una capa unitaria equivariante. Para las matemáticas, véanse J.-P. Serre,
*Linear Representations of Finite Groups* (Springer 1977) y W. Fulton y J. Harris, *Representation Theory*
(Springer 1991); para el aprendizaje equivariante, M. Bronstein, J. Bruna, T. Cohen y P. Veličković, *Geometric Deep
Learning* (2021), y M. Weiler et al., *Equivariant and Coordinate Independent Convolutional Networks* (2023); para el
algoritmo de descomposición, las referencias de `torch_vn_algebra/decompose.py`.
