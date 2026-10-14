[English](../04_symmetries.md) · **Español** · [Français](../fr/04_symmetries.md) · [Deutsch](../de/04_symmetries.md) · [中文](../zh/04_symmetries.md) · [日本語](../ja/04_symmetries.md) · [Русский](../ru/04_symmetries.md)

# 4. Simetrías y representaciones

El capítulo 2 mostró que un grupo de simetría produce un álgebra de von Neumann, el subespacio
generado por los operadores de simetría, y su conmutante, los observables simétricos
(proposición 2.7), y descompuso un ejemplo a mano. Este capítulo aporta la teoría de
representaciones que hay detrás, la convierte en un algoritmo que la biblioteca ejecuta sobre
cualquier lista de matrices, y la aplica a las dos simetrías de los sistemas de espines: $U(1)$ (un
número de partículas o una magnetización) y $SU(2)$ (las rotaciones de los espines).

**Aprenderá** la reducibilidad completa y el lema de Schur, cómo la descomposición en irreducibles
determina el conmutante y el álgebra generada, cómo los caracteres cuentan multiplicidades, por qué un
elemento aleatorio del conmutante encuentra numéricamente la descomposición, qué son los sectores de
carga, las representaciones irreducibles de espín y cómo se descomponen $L$ espines $\frac12$, por qué
los hamiltonianos invariantes por $SU(2)$ tienen bloques mucho más pequeños que los invariantes por
$U(1)$, y cómo la conjugación compleja clasifica los sectores de un álgebra generada por matrices reales en
tipo real, complejo y cuaterniónico. Las comprobaciones numéricas de este capítulo están en el
[Tutorial 4](../../tutorial/es/04_symmetries.md).

En todo el capítulo, $S'$ es el conmutante de $S\subseteq M_n(\mathbb C)$ (capítulo 2). Si $S^*=S$, $S''$ es
el **álgebra de von Neumann generada por** $S$: una $*$-álgebra que contiene $S$, y toda $*$-álgebra
$\mathcal B\supseteq S$ satisface $\mathcal B=\mathcal B''\supseteq S''$ (teorema 2.2). $\mathcal A'_h=\{K\in\mathcal A':K=K^*\}$ es un
espacio vectorial real. $\mathbb C^k\otimes\mathbb C^m$ lleva el orden de Kronecker del capítulo 7 ($e_a\otimes f_j$
es el vector de la base $am+j$, contando desde $0$), de modo que $1_k\otimes x$ es diagonal por bloques con
$k$ copias de $x$. $|G|$ es el orden de un grupo finito, y $\hbar=1$.

<a id="s-4-1"></a>

## 4.1 Representaciones y lema de Schur

**Definición.** Sea $X$ un conjunto con una involución $x\mapsto x^*$ ($x^{**}=x$). Una
**representación** de $X$ sobre $\mathbb C^k$ es una aplicación $\rho:X\to M_k(\mathbb C)$ con
$\rho(x^*)=\rho(x)^*$. Un subespacio $K$ es **invariante** si $\rho(x)K\subseteq K$ para todo $x$; $\rho$ es
**irreducible** si $k\ge1$ y $0,\mathbb C^k$ son los únicos subespacios invariantes. Un
**operador de entrelazamiento** de $\rho$ en una representación $\rho'$ sobre $\mathbb C^{k'}$ es una matriz
$T$ de tamaño $k'\times k$ con $T\rho(x)=\rho'(x)T$ para todo $x$; forman un espacio vectorial
$\mathrm{Hom}(\rho,\rho')$, y $\mathrm{Hom}(\rho,\rho)=\rho(X)'$. $\rho$ y $\rho'$ son **equivalentes**,
$\rho\cong\rho'$, si existe un operador de entrelazamiento invertible.

Se presentan tres casos: un grupo $G$ con $g^*=g^{-1}$, donde exigimos que $\rho$ sea un homomorfismo
$G\to U(k)$ (una **representación unitaria**; entonces $\rho(g^{-1})=\rho(g)^{-1}=\rho(g)^*$); una $*$-álgebra con
su adjunción, donde $\rho$ es un $*$-homomorfismo unital; y etiquetas de generadores con $x^*=x$, por ejemplo
$\{x,y,z\}$, representadas por tuplas de matrices hermíticas, como las matrices de espín (sección 4.5). Si $K$
es invariante y las columnas de $Q\in M_{k\times d}$ forman una base ortonormal de $K$, entonces
$\rho(x)Q=QQ^*\rho(x)Q$, de modo que $\rho_K(x)=Q^*\rho(x)Q$ es de nuevo una representación (un homomorfismo si
$\rho$ lo es), la **restricción** de $\rho$ a $K$.

<a id="r-4-1"></a>

**Lema 4.1 (reducibilidad completa).** Sea $\rho$ una representación de $X$ sobre $\mathbb C^n$.
(a) Si $K$ es invariante, también lo es $K^\perp$. (b) $\mathbb C^n$ es suma directa ortogonal de subespacios
invariantes sobre cada uno de los cuales $\rho$ es irreducible.

*Demostración.* (a) Para $y\in K^\perp$ y $z\in K$, $\langle z,\rho(x)y\rangle=\langle\rho(x^*)z,y\rangle=0$, pues
$\rho(x^*)z\in K$. (b) Inducción sobre $n$. Un subespacio invariante no nulo $K_1$ de dimensión mínima es
irreducible; $K_1^\perp$ es invariante por (a), y los subespacios invariantes de la restricción a
$K_1^\perp$ son subespacios invariantes de $\rho$, de modo que se le aplica la hipótesis de inducción. $\square$

Para un grupo finito, la unitariedad no cuesta nada: si los $\rho(g)$ son meramente invertibles,
$P=\sum_g\rho(g)^*\rho(g)$ satisface $P\ge\rho(e)^*\rho(e)=1$, luego $P>0$, y $\rho(h)^*P\rho(h)=\sum_g\rho(gh)^*\rho(gh)=P$
(sustitúyase $g\to gh^{-1}$); por tanto,
$(P^{1/2}\rho(h)P^{-1/2})^*(P^{1/2}\rho(h)P^{-1/2})=P^{-1/2}\rho(h)^*P\rho(h)P^{-1/2}=1$, y
$h\mapsto P^{1/2}\rho(h)P^{-1/2}$ es una representación unitaria equivalente a $\rho$ (truco de Weyl).

<a id="r-4-2"></a>

**Lema 4.2 (lema de Schur).** Sean $\rho,\rho'$ representaciones irreducibles de $X$ sobre
$\mathbb C^k,\mathbb C^{k'}$ y $T\in\mathrm{Hom}(\rho,\rho')$.
(a) $T=0$ o $T$ es invertible. (b) $\mathrm{Hom}(\rho,\rho)=\mathbb C1$.
(c) $T^*\in\mathrm{Hom}(\rho',\rho)$, y si $T\neq0$, entonces $T=\lambda U$ con $\lambda>0$ y $U$ unitaria.
Por tanto, $\dim\mathrm{Hom}(\rho,\rho')$ es $1$ si $\rho\cong\rho'$ y $0$ en otro caso, y las representaciones
irreducibles equivalentes son unitariamente equivalentes.

*Demostración.* (a) y (b) se demostraron para grupos tras la proposición 2.7, y la demostración solo usa
que $\ker T$ es invariante por $\rho$ ($T\rho(x)v=\rho'(x)Tv=0$ si $Tv=0$) y la imagen de $T$ por
$\rho'$: así, $\ker T\in\{0,\mathbb C^k\}$ y $\mathrm{ran}\,T\in\{0,\mathbb C^{k'}\}$; y para (b), $T-\lambda1$,
con $\lambda$ un valor propio de $T$, es un operador de entrelazamiento no invertible, luego $0$. (c) El adjunto
de $T\rho(x^*)=\rho'(x^*)T$ es $\rho(x)T^*=T^*\rho'(x)$. Entonces $T^*T\in\mathrm{Hom}(\rho,\rho)$, luego
$T^*T=\lambda^21$ con $\lambda>0$ si $T\neq0$ ($T^*T\ge0$ es no nulo), y $U=T/\lambda$ es cuadrada por
(a) con $U^*U=1$. Finalmente, si $T_1,T_2\in\mathrm{Hom}(\rho,\rho')$ son no nulos, entonces
$T_1^{-1}T_2\in\mathrm{Hom}(\rho,\rho)=\mathbb C1$. $\square$

Si todos los $\rho(x)$ conmutan entre sí (un grupo abeliano), una $\rho$ irreducible es unidimensional:
todo $\rho(x)$ está en $\mathrm{Hom}(\rho,\rho)=\mathbb C1$, de modo que todo subespacio es invariante.

<a id="s-4-2"></a>

## 4.2 La descomposición isotípica y el conmutante

<a id="r-4-3"></a>

**Teorema 4.3 (descomposición isotípica).** Sea $\rho$ una representación de $X$ sobre $\mathbb C^n$.
Existen representaciones irreducibles, no equivalentes dos a dos, $\rho_1,\dots,\rho_C$ de $X$ sobre
$\mathbb C^{k_c}$, enteros $m_c\ge1$ con $\sum_ck_cm_c=n$, y una unitaria
$V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n$ con
$V^*\rho(x)V=\bigoplus_c\rho_c(x)\otimes1_{m_c}$ para todo $x$. Para toda descomposición de este tipo
$$
V^*\rho(X)'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),\qquad
V^*\rho(X)''V=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},
$$
y $m_c=\dim\mathrm{Hom}(\rho_c,\rho)$; una representación irreducible $\sigma$ de $X$ cumple
$\mathrm{Hom}(\sigma,\rho)\neq0$ si y solo si $\sigma\cong\rho_c$ para exactamente un $c$.

*Demostración.* *Existencia.* El lema 4.1 da subespacios invariantes ortogonales dos a dos $W_1,\dots,W_M$ de
suma $\mathbb C^n$ y con restricciones irreducibles $\sigma_i(x)=Q_i^*\rho(x)Q_i$. Agrupamos los $W_i$ en clases
de restricciones equivalentes y elegimos en cada clase $c$ una referencia $\rho_c=\sigma_{i_0}$. Para cada otro
$i$ de la clase, el lema 4.2(c) da una unitaria $U_i$ con $\sigma_i(x)=U_i\rho_c(x)U_i^*$; sustituir $Q_i$
por $Q_iU_i$ hace que $Q_i^*\rho(x)Q_i=\rho_c(x)$. Numeramos los miembros de la clase $c$ con $j=0,\dots,m_c-1$,
escribimos $Q_{c,j}$ para sus bases y ponemos $V(e_a\otimes f_j)=Q_{c,j}e_a$ en el $c$-ésimo sumando. Como los
$W_i$ son invariantes y ortogonales, $Q_{c,j}^*\rho(x)Q_{d,l}=\delta_{cd}\delta_{jl}\rho_c(x)$, que es la forma
afirmada.

*Conmutante.* Sea ahora $V$ cualquier descomposición como en el enunciado, y sea $Q_{c,j}$ la isometría
$n\times k_c$ $e_a\mapsto V(e_a\otimes f_j)$ (en el $c$-ésimo sumando); para la $V$ recién construida son las $Q_{c,j}$
de arriba. La identidad $V^*\rho(x)V=\bigoplus_c\rho_c(x)\otimes1_{m_c}$ dice que $\rho(x)Q_{d,l}=Q_{d,l}\rho_d(x)$, y su
adjunta aplicada a $x^*$ da $Q_{c,j}^*\rho(x)=\rho_c(x)Q_{c,j}^*$. Como $V$ es unitaria, $Y\in M_n$ conmuta con todos los $\rho(x)$ si y solo si
cada bloque $Y_{cj,dl}=Q_{c,j}^*YQ_{d,l}$ satisface $Y_{cj,dl}\rho_d(x)=\rho_c(x)Y_{cj,dl}$, es decir, pertenece a
$\mathrm{Hom}(\rho_d,\rho_c)$, que es $0$ para $c\neq d$ y $\mathbb C1$ para $c=d$ (lema 4.2). Así,
$Y\in\rho(X)'$ si y solo si $Y_{cj,dl}=\delta_{cd}(y_c)_{jl}1_{k_c}$, es decir, $V^*YV=\bigoplus_c1_{k_c}\otimes y_c$.

*Doble conmutante.* Si $Z$ conmuta con $\rho(X)'$, conmuta con las proyecciones sobre los sumandos
(tómese $y_c=1$ e $y_d=0$ para $d\neq c$), luego $V^*ZV=\bigoplus_cZ_c$, y $Z_c$ conmuta con
$1_{k_c}\otimes M_{m_c}$; por el lema 2.5, paso 5 (con los dos factores tensoriales intercambiados),
$Z_c\in M_{k_c}\otimes1_{m_c}$. Recíprocamente, estos operadores conmutan con $\rho(X)'$.

*Multiplicidades.* Para $T\in\mathrm{Hom}(\sigma,\rho)$ con $\sigma$ irreducible, el mismo argumento muestra
que $Q_{c,j}^*T\in\mathrm{Hom}(\sigma,\rho_c)$, que es $\mathbb Ct_{c}$ con $t_c$ invertible si $\sigma\cong\rho_c$
y $0$ en otro caso; y $T=\sum_{c,j}Q_{c,j}Q_{c,j}^*T$. Como las $\rho_c$ son no equivalentes dos a dos, $\sigma$
es equivalente a lo sumo a una de ellas, digamos $\rho_c$ (si no, $T=0$), y entonces $T=\sum_ja_jQ_{c,j}t_c$;
recíprocamente, cada $Q_{c,j}t_c$ pertenece a $\mathrm{Hom}(\sigma,\rho)$, y estas $m_c$ aplicaciones tienen
imágenes ortogonales, luego $\dim\mathrm{Hom}(\sigma,\rho)=m_c$. $\square$

*En la biblioteca:* una descomposición como en el teorema 4.3 la calculan `decompose.from_group` y `decompose.from_generators` — [tutorial, sección 4.3](../../tutorial/es/04_symmetries.md#s-4-3).

Para un grupo, $\rho(G)''$ es el subespacio generado por $\rho(G)$ (sección 2.4), y el teorema 4.3 es la
proposición 2.7 leída en sentido inverso: allí el teorema de estructura producía las representaciones
irreducibles, aquí son ellas las que producen la estructura; en ambos casos, $\dim\rho(X)'=\sum_cm_c^2$ y
$\dim\rho(X)''=\sum_ck_c^2$. La imagen de $P_c=\sum_jQ_{c,j}Q_{c,j}^*$, una proyección central minimal de
$\rho(X)''$, es la **componente isotípica** de $\rho_c$.
Por el teorema 4.3, todo $H\in\rho(X)'$ hermítico es $V(\bigoplus_c1_{k_c}\otimes h_c)V^*$ con
$h_c\in M_{m_c}$ hermíticas, y $1_{k_c}\otimes h_c$ es diagonal por bloques con $k_c$ copias de $h_c$: el espectro de $H$
es la unión de los espectros de los bloques $h_c$ de tamaño $m_c\times m_c$, y cada valor propio de $h_c$ se repite
$k_c$ veces.

> **Significado físico** (interpretación física, no un teorema de este curso). Si $\rho(X)$ es una simetría
> y $H\in\rho(X)'$ el hamiltoniano, los niveles $k_c$ veces degenerados que se acaban de describir son los
> **multipletes** de la simetría, y los bloques que hay que diagonalizar tienen los tamaños $m_c$ de las
> multiplicidades, no los tamaños de los multipletes.

Para un grupo finito, las multiplicidades pueden calcularse sin descomponer nada.

<a id="r-4-4"></a>

**Teorema 4.4 (ortogonalidad de Schur).** Sean $\rho,\rho'$ representaciones unitarias irreducibles de un
grupo finito $G$ sobre $\mathbb C^k,\mathbb C^{k'}$. Entonces
$$
\frac1{|G|}\sum_{g\in G}\rho(g)_{ij}\,\overline{\rho'(g)_{i'j'}}=
\begin{cases}\delta_{ii'}\delta_{jj'}/k&\text{if }\rho=\rho',\\0&\text{if }\rho\not\cong\rho'.\end{cases}
$$

*Demostración.* Para una matriz $Y$ de tamaño $k\times k'$, sea $T=\frac1{|G|}\sum_g\rho(g)Y\rho'(g)^*$. Para
$h\in G$, la sustitución $g\to hg$ y la igualdad $\rho'(hg)^*=\rho'(g)^*\rho'(h)^{-1}$ dan $\rho(h)T=T\rho'(h)$, luego
$T\in\mathrm{Hom}(\rho',\rho)$. Si $\rho\not\cong\rho'$, $T=0$; si $\rho=\rho'$, $T=\lambda1$ con
$k\lambda=\mathrm{Tr}\,T=\mathrm{Tr}\,Y$ (lema 4.2). Para $Y=|j\rangle\langle j'|$, la entrada $T_{ii'}$
es el miembro izquierdo, y $\mathrm{Tr}\,Y=\delta_{jj'}$. $\square$

El **carácter** de una representación es $\chi_\rho(g)=\mathrm{Tr}\,\rho(g)$; las representaciones
equivalentes tienen caracteres iguales, pues $\mathrm{Tr}(T\rho(g)T^{-1})=\mathrm{Tr}\,\rho(g)$. Escribimos
$\langle\varphi,\psi\rangle_G=\frac1{|G|}\sum_g\overline{\varphi(g)}\psi(g)$.

<a id="r-4-5"></a>

**Corolario 4.5 (caracteres).** Sea $G$ finito.
(a) Para $\rho,\rho'$ unitarias irreducibles: $\langle\chi_\rho,\chi_{\rho'}\rangle_G$ es $1$ si $\rho\cong\rho'$ y $0$ en otro caso.
(b) Para una representación unitaria $\pi$ descompuesta como en el teorema 4.3 y una $\rho$ irreducible, la
multiplicidad de $\rho$ en $\pi$ es
$$
m_\rho=\frac1{|G|}\sum_{g\in G}\overline{\chi_\rho(g)}\,\chi_\pi(g)
$$
($m_\rho=m_c$ si $\rho\cong\rho_c$, y $m_\rho=0$ si $\rho$ no es equivalente a ninguna $\rho_c$).
(c) $\dim\pi(G)'=\sum_cm_c^2=\frac1{|G|}\sum_{g\in G}|\chi_\pi(g)|^2$. En particular, $\pi$ es irreducible si y solo
si este número es $1$.

*Demostración.* (a) Sustituimos $\rho'$ por una representación equivalente (mismo carácter). Sumamos el
teorema 4.4 sobre $i=j$, $i'=j'$ y conjugamos: el resultado es $\sum_{i,i'}\delta_{ii'}/k=1$ si $\rho=\rho'$, y $0$
si $\rho\not\cong\rho'$. (b) La traza es invariante unitaria y $\mathrm{Tr}(a\otimes1_m)=m\mathrm{Tr}\,a$, luego
$\chi_\pi=\sum_cm_c\chi_{\rho_c}$; se aplica (a). (c) Por el mismo desarrollo y (a),
$\langle\chi_\pi,\chi_\pi\rangle_G=\sum_cm_c^2$, que es $\dim\pi(G)'$ por el teorema 4.3. $\square$

*Cálculo:* el corolario 4.5 para $S_3$ permutando tres cúbits — [tutorial, sección 4.2](../../tutorial/es/04_symmetries.md#s-4-2); el teorema 4.4 para la representación estándar de $S_3$ — [tutorial, sección 4.3](../../tutorial/es/04_symmetries.md#s-4-3).

Una segunda demostración de (c) es la forma en que la biblioteca muestrea conmutantes. El **twirl**
$P(Y)=\frac1{|G|}\sum_g\pi(g)Y\pi(g)^*$ es idempotente con imagen $\pi(G)'$ (capítulo 2, ejercicio 3) y
autoadjunto para $\langle Y,Z\rangle_{HS}=\mathrm{Tr}(Y^*Z)$ (sustitúyase $g\to g^{-1}$), luego es la proyección
ortogonal sobre $\pi(G)'$, y $\dim\pi(G)'$ es su traza. La aplicación $Y\mapsto AYB$ tiene traza
$\mathrm{Tr}\,A\mathrm{Tr}\,B$ (sus entradas diagonales son $A_{ii}B_{jj}$), lo que da $\frac1{|G|}\sum_g|\chi_\pi(g)|^2$.

*En la biblioteca:* el twirl es `decompose.twirl` — [tutorial, sección 4.3](../../tutorial/es/04_symmetries.md#s-4-3); $\dim\pi(G)'$ es `decompose.commutant_dimension` — [tutorial, sección 4.2](../../tutorial/es/04_symmetries.md#s-4-2).

Rehacemos con caracteres el ejemplo de la sección 2.4. Una matriz de permutación tiene como traza el número
de vectores de la base que deja fijos, y $\pi(\sigma)$ fija $|b_1b_2b_3\rangle$ si y solo si $b$ es constante en
los ciclos de $\sigma$ (lema B.1(c)), luego $\chi_\pi(\sigma)=2^{\#\text{cycles}(\sigma)}$. Tenemos a mano tres
representaciones unitarias de $S_3$. Sea $P_\sigma e_i=e_{\sigma(i)}$ la representación de permutación sobre $\mathbb C^3$
(un homomorfismo en las matrices de permutación). La representación **trivial** es $\sigma\mapsto1$. El **signo** es
$\mathrm{sgn}\,\sigma=\det P_\sigma$, un homomorfismo en $\{\pm1\}\subseteq U(1)$ porque el determinante es
multiplicativo y el determinante de una matriz de permutación es $\pm1$. La recta $\mathbb C(1,1,1)$ es invariante,
luego también lo es $E=(1,1,1)^\perp$ (lema 4.1(a)); la representación **estándar** es la restricción de $P$ a $E$.
En una base ortonormal adaptada a $\mathbb C(1,1,1)\oplus E$, $P_\sigma=1\oplus\mathrm{std}(\sigma)$, luego
$\chi_{\rm std}(\sigma)=\mathrm{Tr}\,P_\sigma-1$, el número de puntos fijos menos uno. Los seis elementos de $S_3$ son
la identidad, tres transposiciones y dos $3$-ciclos, con $3,2,1$ ciclos y $3,1,0$ puntos fijos; una transposición
intercambia dos columnas de $1_3$, luego su signo es $-1$, y un $3$-ciclo es producto de dos transposiciones
($(1\,2\,3)=(1\,3)(1\,2)$), luego su signo es $+1$. Esto da la tabla

| | identidad ($1$) | transposiciones ($3$) | $3$-ciclos ($2$) |
|---|---|---|---|
| $\chi_{\rm triv}$ | $1$ | $1$ | $1$ |
| $\mathrm{sgn}$ | $1$ | $-1$ | $1$ |
| $\chi_{\rm std}$ | $2$ | $0$ | $-1$ |
| $\chi_\pi$ | $8$ | $4$ | $2$ |

(así, $\mathrm{sgn}\,\sigma=(-1)^{3-\#\text{cycles}(\sigma)}$; todos los caracteres son reales). Con
$\langle\varphi,\psi\rangle_G$ sumado clase por clase:
$\langle\chi_{\rm triv},\chi_{\rm triv}\rangle_G=\frac{1+3+2}6=1$, $\langle\mathrm{sgn},\mathrm{sgn}\,\rangle_G=\frac{1+3+2}6=1$ y
$\langle\chi_{\rm std},\chi_{\rm std}\rangle_G=\frac{4+0+2}6=1$, luego las tres representaciones son irreducibles por el
corolario 4.5(c), y son no equivalentes dos a dos porque sus caracteres difieren. El corolario 4.5(b) da las
multiplicidades en $\pi$:
$$
m_{\rm triv}=\tfrac{8+3\cdot4+2\cdot2}6=4,\qquad m_{\rm sgn}=\tfrac{8-3\cdot4+2\cdot2}6=0,\qquad
m_{\rm std}=\tfrac{2\cdot8+0-2\cdot2}6=2 .
$$
Así, la representación trivial y la estándar son equivalentes a constituyentes $\rho_c$ de $\pi$ (teorema 4.3) con
$m_c=4$ y $m_c=2$; suman $4\cdot1+2\cdot2=8=\sum_ck_cm_c$ dimensiones, luego no hay ningún otro constituyente:
$k=(1,2)$, $m=(4,2)$, y $\dim\pi(S_3)'=4^2+2^2=\frac16(64+3\cdot16+2\cdot4)=20$ (corolario 4.5(c)).

Es la descomposición de la sección 2.4, obtenida ahora a partir de una tabla de seis números por representación.

*Cálculo:* esta tabla, recalculada a partir de $\pi$ — [tutorial, sección 4.2](../../tutorial/es/04_symmetries.md#s-4-2).

<a id="s-4-3"></a>

## 4.3 De los generadores al álgebra: un algoritmo aleatorizado

El teorema 4.3 describe la descomposición, pero no la calcula: su demostración necesita los subespacios
invariantes irreducibles. Murota, Kanno, Kojima y Kojima (*Japan J. Indust. Appl. Math.* 27 (2010) 125–160;
véase también Maehara y Murota, *ibid.* 263–293) observaron que los subespacios propios de un *elemento
aleatorio del conmutante* los encuentran. Sea $\mathcal A$ el álgebra de von Neumann generada por unas matrices
dadas $S_1,\dots,S_J$ y sus adjuntas, con una descomposición $V^*\mathcal AV=\bigoplus_cM_{k_c}\otimes1_{m_c}$,
$V^*\mathcal A'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$ (teorema 2.6); solo el análisis usa $V$.

<a id="r-4-6"></a>

**Lema 4.6.** Sea $K\in\mathcal A'_h$ con $V^*KV=\bigoplus_c1_{k_c}\otimes\kappa_c$. Todo subespacio propio de
$K$ es $\mathcal A$-invariante. Si los $N=\sum_cm_c$ valores propios de todos los $\kappa_c$ juntos son distintos
dos a dos, los subespacios propios de $K$ son los $N$ subespacios $V(\mathbb C^{k_c}\otimes e)$ (en el $c$-ésimo
sumando, con $e$ un vector propio de $\kappa_c$), y $\mathcal A$ actúa irreduciblemente sobre cada uno de ellos,
como $M_{k_c}$.

*Demostración.* Si $Kv=\lambda v$ y $a\in\mathcal A$, entonces $Kav=aKv=\lambda av$. Escribimos
$\kappa_c=\sum_j\mu_{cj}|e_{cj}\rangle\langle e_{cj}|$ con $e_{cj}$ ortonormales (teorema 1.3). Entonces
$V^*KV=\sum_{c,j}\mu_{cj}R_{cj}$, donde $R_{cj}$ es la proyección sobre $\mathbb C^{k_c}\otimes e_{cj}$; las
$R_{cj}$ son ortogonales con suma $1$ y los $\mu_{cj}$ distintos, luego estas son las proyecciones espectrales.
Sobre $\mathbb C^{k_c}\otimes e_{cj}$ el elemento $V^*aV=\bigoplus_da_d\otimes1$ actúa como $a_c$, y los $a_c$
recorren $M_{k_c}$, que no tiene más subespacios invariantes que $0$ y $\mathbb C^{k_c}$. $\square$

<a id="r-4-7"></a>

**Proposición 4.7 (genericidad).** Sean $K$ e $Y$ elementos aleatorios independientes de $\mathcal A'_h$ y de
$\mathcal A'$ cuyas leyes tienen densidades respecto de la medida de Lebesgue en estos espacios vectoriales
reales (lema A.3). Con probabilidad uno: (a) se cumple la hipótesis del lema 4.6; (b) $\Pi_iY\Pi_j\neq0$ para
cualesquiera dos proyecciones espectrales distintas $\Pi_i,\Pi_j$ de $K$ cuyas imágenes estén en el mismo sumando.

*Demostración.* (a) La matriz hermítica $\kappa=\bigoplus_c\kappa_c\in M_N$ depende de forma real-lineal de $K$. El
determinante de Gram $\delta(K)=\det[\mathrm{Tr}\,\kappa^{i+j}]_{i,j=0}^{N-1}$ de $1,\kappa,\dots,\kappa^{N-1}$ para
$\langle\cdot,\cdot\rangle_{HS}$ es un polinomio real en las coordenadas de $K$, no nulo si y solo si estas potencias
son linealmente independientes. Como $p(\kappa)=U\mathrm{diag}(p(\lambda_l))U^*$ para
$\kappa=U\mathrm{diag}(\lambda)U^*$, un polinomio no nulo de grado $<N$ anula a $\kappa$ si y solo si se anula
en $\mathrm{spec}\,\kappa$, lo cual es posible si y solo si $\kappa$ tiene menos de $N$ valores propios distintos
(tómese $\prod_\mu(t-\mu)$ sobre ellos). Así, $\delta(K)\neq0$ es la hipótesis del lema 4.6, y $\delta\not\equiv0$:
es no nulo en $V(\bigoplus_c1\otimes\mathrm{diag}(\text{reales distintos}))V^*\in\mathcal A'_h$. Por tanto,
$\mathbb P(\delta(K)=0)=0$ por el corolario A.7(b).

(b) El conjunto $O=\{\delta\neq0\}\subseteq\mathcal A'_h$ es abierto ($\delta$ es continua). Para $K\in O$, $K$ tiene
exactamente $N$ valores propios distintos $\mu_1(K)<\dots<\mu_N(K)$ (los de $\kappa$, con multiplicidades $k_c$), así que
sus proyecciones espectrales $\Pi_1(K),\dots,\Pi_N(K)$, numeradas por valor propio creciente, son continuas en $O$
(proposición A.13). Fijemos $l\neq l'$ y sea $\mathcal N_{ll'}$ el conjunto de los $(K,Y)\in O\times\mathcal A'$ con
$P_c\Pi_l(K)=\Pi_l(K)$ y $P_c\Pi_{l'}(K)=\Pi_{l'}(K)$ para algún $c$ ($P_c$ las proyecciones centrales) y
$\Pi_l(K)Y\Pi_{l'}(K)=0$; está definido por ecuaciones continuas sobre el conjunto de Borel $O\times\mathcal A'$, luego es
de Borel. Fijemos $K\in O$ para el que las dos imágenes están en el sumando $c$. Por el lema 4.6, $V^*\Pi_lV$ y
$V^*\Pi_{l'}V$ proyectan sobre $\mathbb C^{k_c}\otimes e$ y $\mathbb C^{k_c}\otimes e'$ con vectores propios ortonormales
$e,e'$ de $\kappa_c$, y para $V^*YV=\bigoplus_d1\otimes y_d$, $\Pi_lY\Pi_{l'}=0$ si y solo si
$\langle e,y_ce'\rangle=0$. La aplicación $Y\mapsto\langle e,y_ce'\rangle$ es real-lineal $\mathcal A'\to\mathbb C$ y no
nula (su valor en $V(0\oplus\dots\oplus1\otimes|e\rangle\langle e'|\oplus\dots\oplus0)V^*$ es $1$), luego su núcleo es un
subespacio propio, un conjunto nulo (corolario A.7(c)). Así, toda sección $\{Y:(K,Y)\in\mathcal N_{ll'}\}$ es nula (es
vacía si las imágenes están en sumandos distintos) y, como $Y$ es independiente de $K$,
$\mathbb P((K,Y)\in\mathcal N_{ll'})=0$ por el corolario A.7(d). El suceso de (b) solo falla en
$\{\delta(K)=0\}\cup\bigcup_{l\neq l'}\{(K,Y)\in\mathcal N_{ll'}\}$, una unión finita de sucesos de probabilidad $0$.
$\square$

<a id="r-4-8"></a>

**Lema 4.8 (compresiones; alineación de las copias).** (a) Sea $W$ un subespacio $\mathcal A$-invariante con base
ortonormal $Q\in M_{n\times d}$ y $\sigma(a)=Q^*aQ$ la restricción. Entonces $\sigma(\mathcal A)$ es una $*$-álgebra,
$\{Q^*YQ:Y\in\mathcal A'\}=\sigma(\mathcal A)'$ y $\{Q^*KQ:K\in\mathcal A'_h\}=\sigma(\mathcal A)'_h$. $\mathcal A$ actúa
irreduciblemente sobre $W$ si y solo si $Q^*YQ\in\mathbb C1$ para todo $Y\in\mathcal A'$. Si no es así, y $X$ es un
elemento aleatorio de $\mathcal A'$ o de $\mathcal A'_h$ con densidad, entonces $Q^*XQ\notin\mathbb C1$ con probabilidad uno.

(b) Sean $W_i,W_j$ subespacios invariantes sobre los que $\mathcal A$ actúa
irreduciblemente, con bases ortonormales $Q_i,Q_j$ y restricciones $\sigma_i(a)=Q_i^*aQ_i$, y sea $Y\in\mathcal A'$.
Entonces $T=Q_i^*YQ_j\in\mathrm{Hom}(\sigma_j,\sigma_i)$. Si $T\neq0$, entonces $T=\lambda U$ con $\lambda>0$, $U$
unitaria, y $Q_j'=Q_jU^*$ satisface $Q_j'^*aQ_j'=\sigma_i(a)$ para todo $a\in\mathcal A$. Si $W_i,W_j$ están en
sumandos distintos ($W_i\subseteq P_c\mathbb C^n$, $W_j\subseteq P_d\mathbb C^n$, $c\neq d$), entonces $T=0$.

*Demostración.* (a) Como en la demostración del teorema 4.3 ($W$ es invariante por $a$ y por $a^*$), $aQ=Q\sigma(a)$ y
$Q^*a=\sigma(a)Q^*$ para $a\in\mathcal A$. Por tanto, $\sigma(ab)=Q^*aQ\sigma(b)=\sigma(a)\sigma(b)$, $\sigma(a^*)=\sigma(a)^*$ y
$\sigma(1)=Q^*Q=1_d$, de modo que $\sigma$ es un $*$-homomorfismo unital y $\sigma(\mathcal A)$ es una $*$-álgebra; y para
$Y\in\mathcal A'$, $Q^*YQ\sigma(a)=Q^*YaQ=Q^*aYQ=\sigma(a)Q^*YQ$. Recíprocamente, para $T\in\sigma(\mathcal A)'$ sea
$Y=QTQ^*$; entonces $aY=Q\sigma(a)TQ^*=QT\sigma(a)Q^*=Ya$, luego $Y\in\mathcal A'$, con $Q^*YQ=T$, e $Y$ es hermítica si
$T$ lo es. Si $\sigma$ es irreducible, $\sigma(\mathcal A)'=\mathrm{Hom}(\sigma,\sigma)=\mathbb C1$ (lema 4.2(b)). Si
$E\subseteq\mathbb C^d$ es un subespacio $\sigma$-invariante distinto de $0,\mathbb C^d$, también lo es $E^\perp$ (lema
4.1(a)), y la proyección ortogonal sobre $E$ es un elemento no escalar de $\sigma(\mathcal A)'_h$. En ese caso,
$\mathbb C1\cap\sigma(\mathcal A)'$ y $\mathbb R1\cap\sigma(\mathcal A)'_h$ son subespacios reales propios, sus preimágenes
por las aplicaciones real-lineales sobreyectivas $X\mapsto Q^*XQ$ son subespacios propios de $\mathcal A'$ y de
$\mathcal A'_h$, y estos son conjuntos nulos (corolario A.7(c)).

(b) Como en (a), $aQ_j=Q_j\sigma_j(a)$ y $Q_i^*a=\sigma_i(a)Q_i^*$, luego
$T\sigma_j(a)=Q_i^*YaQ_j=Q_i^*aYQ_j=\sigma_i(a)T$. El lema 4.2 da $T=\lambda U$ y
$U\sigma_j(a)U^*=\sigma_i(a)$, que es $Q_j'^*aQ_j'$. Finalmente, $P_c\in\mathcal A\cap\mathcal A'$, $Q_i=P_cQ_i$ y
$P_cQ_j=0$, luego $T=Q_i^*P_cYQ_j=Q_i^*YP_cQ_j=0$. $\square$

*Cálculo:* el lema 4.6 y el lema 4.8(b) para un elemento aleatorio de $\pi(S_3)'$ — [tutorial, sección 4.3](../../tutorial/es/04_symmetries.md#s-4-3).

El algoritmo (`decompose.wedderburn`, con las interfaces `from_group` y `from_generators`) es ahora:

1. *Muestrear el conmutante.* Para todos los elementos de un grupo finito, $Y=P(Z)$ y $K=P(Z'+Z'^*)$ con el twirl
   $P$ y matrices de Ginibre independientes $Z,Z'$. $P$ es una proyección lineal de $M_n$ sobre $\mathcal A'$ con
   $P(Z)^*=P(Z^*)$, así que $Z\mapsto P(Z)$ y $Z'\mapsto P(Z'+Z'^*)$ son sobreyecciones real-lineales sobre $\mathcal A'$ y
   sobre $\mathcal A'_h$ (para $H\in\mathcal A'_h$, $P(\frac H2+\frac H2)=H$); una matriz de Ginibre es un vector gaussiano
   centrado no degenerado de $M_n$ (proposición A.16), luego $Y$ y $K$ son vectores gaussianos centrados no degenerados
   de $\mathcal A'$ y de $\mathcal A'_h$ y tienen densidades (teorema A.15). Para generadores,
   $\mathcal A'=\{X:[S_j,X]=[S_j^*,X]=0\}$ es el núcleo de una aplicación lineal de $M_n$ (la biblioteca usa las partes
   hermíticas $\frac12(S_j+S_j^*)$, $\frac1{2i}(S_j-S_j^*)$, que generan el mismo espacio complejo que $S_j,S_j^*$ y
   tienen, por tanto, el mismo conmutante); con una base $(B_r)_{r\le D}$ de él, $Y=\sum_rc_rB_r$ y $K=Y'+Y'^*$, con $Y'$
   una copia independiente de $Y$ y $c_r$ gaussianas complejas estándar independientes: imágenes de un vector
   gaussiano no degenerado de $\mathbb C^D$ por la biyección real-lineal $c\mapsto\sum_rc_rB_r$ sobre $\mathcal A'$ y por la
   sobreyección real-lineal $c\mapsto X+X^*$ ($X=\sum_rc_rB_r$) sobre $\mathcal A'_h$, y por tanto de nuevo con densidades.
   Cada elemento aleatorio es una extracción nueva, independiente de las anteriores; la biblioteca los reescala a norma
   unidad, lo que no cambia ningún subespacio propio ni la anulación de ningún bloque.
2. *Separar.* Los subespacios propios de $K$ son irreducibles (lema 4.6, proposición 4.7(a)); sus dimensiones
   son los $k_c$.
3. *Agrupar y alinear.* Los bloques $Q_i^*YQ_j$ entre subespacios propios de igual dimensión son no nulos exactamente
   dentro de un sumando (proposición 4.7(b), lema 4.8(b)); las clases de la relación «unidos por un bloque no nulo» dan
   los sectores, sus tamaños los $m_c$, y los factores polares alinean las copias con una copia de referencia por clase
   (lema 4.8(b); una copia alineada con una copia alineada está alineada con la referencia).
4. *Ensamblar* $V$ en el orden de Kronecker, de modo que $V^*aV=\bigoplus_ca_c\otimes1_{m_c}$.

Así, en aritmética exacta, con los subespacios propios de $K$ separados exactamente, el algoritmo devuelve con
probabilidad uno una descomposición como en el teorema 4.3: por el lema 4.6, el sumando $c$ contiene $m_c$ subespacios
propios de dimensión $k_c$, y tras la alineación la construcción es la parte de existencia de la demostración del
teorema 4.3. En coma flotante, la biblioteca trata como iguales los valores propios más cercanos que una tolerancia,
de modo que un «subespacio propio» calculado $W$ puede ser suma de varios subespacios propios de $K$; sigue siendo
$\mathcal A$-invariante (lema 4.6), y la biblioteca comprime en él un nuevo $X\in\mathcal A'_h$ aleatorio. Por el lema
4.8(a), $Q^*XQ$ es un múltiplo de $1$ siempre que $\mathcal A$ actúe irreduciblemente sobre $W$, y en otro caso no lo es
con probabilidad uno; entonces sus subespacios propios, invariantes por $\sigma(\mathcal A)$ (lema 4.6 aplicado a
$\sigma(\mathcal A)$), dividen $W$ en subespacios invariantes más pequeños, y la biblioteca repite hasta que toda
compresión es escalar. Pase lo que pase en coma flotante, el resultado se certifica al final:

<a id="r-4-9"></a>

**Proposición 4.9 (certificado).** Sea $\mathcal A$ generada por $S_1,\dots,S_J$ y sus adjuntas, y
$\mathcal B=V(\bigoplus_cM_{k_c}\otimes1_{m_c})V^*$ con $V$ unitaria. Si todas las $S_j\in\mathcal B$ y
$\dim\mathcal A'=\sum_cm_c^2$, entonces $\mathcal A=\mathcal B$.

*Demostración.* $\mathcal B$ es una $*$-álgebra que contiene todas las $S_j,S_j^*$, luego $\mathcal A\subseteq\mathcal B$ y
$\mathcal B'\subseteq\mathcal A'$. El teorema 4.3, aplicado a la representación de $\mathcal B$ sobre $\mathbb C^n$ (las
aplicaciones de bloque $b\mapsto b_c$ son irreducibles, y no equivalentes, pues la $c$-ésima proyección central se
aplica en $1$ por una de ellas y en $0$ por las demás), da $\dim\mathcal B'=\sum_cm_c^2=\dim\mathcal A'$. Así,
$\mathcal B'=\mathcal A'$ y $\mathcal B=\mathcal B''=\mathcal A''=\mathcal A$. $\square$

*En la biblioteca:* el algoritmo, `decompose.wedderburn` (interfaces `from_group`, `from_generators`), el certificado de la proposición 4.9, `Decomposition.verify()`, y su coste — [tutorial, sección 4.3](../../tutorial/es/04_symmetries.md#s-4-3).

<a id="s-4-4"></a>

## 4.4 $U(1)$ y sectores de carga

<a id="r-4-10"></a>

**Proposición 4.10.** Sea $Q$ una matriz hermítica sobre $\mathbb C^n$ con valores propios enteros $q$ (las
**cargas**), proyecciones espectrales $P_q$ y $d_q=\mathrm{rank}\,P_q$, y sea $U_\theta=e^{i\theta Q}$; escribimos
$U(1)'=\{U_\theta:\theta\in\mathbb R\}'$ y $U(1)''$ para su conmutante.

1. $\theta\mapsto U_\theta$ es una representación unitaria de $U(1)=\mathbb R/2\pi\mathbb Z$, y toda representación
   unitaria continua de $U(1)$ sobre $\mathbb C^n$ es de esta forma.
2. Para $X\in M_n$ son equivalentes: (i) $XU_\theta=U_\theta X$ para todo $\theta$; (ii) $[X,Q]=0$;
   (iii) $X=\sum_qP_qXP_q$. Por tanto, $U(1)'\cong\bigoplus_qM_{d_q}(\mathbb C)$ (los operadores diagonales por
   bloques) y $U(1)''=\mathrm{span}\{P_q\}$: todos los $k_c=1$, y las multiplicidades son los $d_q$.
3. (*Carga conservada.*) Para un $H$ hermítico, $e^{itH}Qe^{-itH}=Q$ para todo $t$ si y solo si $[H,Q]=0$, es
   decir, si y solo si $H\in U(1)'$. Entonces todo peso de sector $\mathrm{Tr}(P_q\rho(t))$ de
   $\rho(t)=e^{-itH}\rho e^{itH}$ es constante.

*Demostración.* (1) $U_\theta=\sum_qe^{iq\theta}P_q$ es unitaria, $U_\theta U_{\theta'}=U_{\theta+\theta'}$ (teorema 1.5(c)) y
$U_{2\pi}=1$. Recíprocamente, para una representación continua, $(U_\theta)_{\theta\ge0}$ es un semigrupo continuo
de matrices, luego $U_\theta=e^{\theta B}$ para $\theta\ge0$ por el lema 8.5 (demostrado en el capítulo 8; la
demostración solo usa la dimensión finita), y para $\theta<0$, $U_\theta=U_{-\theta}^{-1}=e^{\theta B}$. Derivando
$U_\theta^*U_\theta=1$ en $0$ (lema 8.1) se obtiene $B^*+B=0$, luego $B=iQ$ con $Q$ hermítica. La serie exponencial
$e^{i\theta Q}$ coincide con el cálculo funcional (proposición 8.2(1) con $H=-Q$), luego $U_{2\pi}=\sum_qe^{2\pi iq}P_q=1$
obliga a $e^{2\pi iq}=1$, $q\in\mathbb Z$, para todo valor propio $q$ de $Q$. (2) (i)$\Rightarrow$(ii): se deriva en $\theta=0$.
(ii)$\Rightarrow$(iii): $X$ conmuta con $P_q=f(Q)$ (teorema 1.5(d)), luego $P_qXP_{q'}=XP_qP_{q'}=0$ para
$q\neq q'$. (iii)$\Rightarrow$(i): $XU_\theta=\sum_qe^{iq\theta}P_qXP_q=U_\theta X$. Los operadores diagonales por bloques
forman $\bigoplus_qM_{d_q}\otimes1_1$ (en una base ortonormal adaptada a las imágenes de los $P_q$), cuyo conmutante es
$\bigoplus_q1_{d_q}\otimes M_1=\mathrm{span}\{P_q\}$: un operador que conmuta con ellos conmuta con cada $P_q$,
luego es diagonal por bloques, y cada bloque conmuta con $M_{d_q}$, luego es un escalar (lema 2.5, paso 5, con $m=1$).
(3) $\frac d{dt}e^{itH}Qe^{-itH}=ie^{itH}[H,Q]e^{-itH}$ (proposición 8.2) se anula para todo $t$ si y
solo si $[H,Q]=0$. Entonces $e^{itH}$ conmuta con $P_q$ (teorema 1.5(d), dos veces), y
$\mathrm{Tr}(P_q\rho(t))=\mathrm{Tr}(e^{itH}P_qe^{-itH}\rho)=\mathrm{Tr}(P_q\rho)$. $\square$

Las representaciones irreducibles de $U(1)$ son, pues, los caracteres $\theta\mapsto e^{iq\theta}$, $q\in\mathbb Z$
(unidimensionales, como para todo grupo abeliano), y la multiplicidad de la carga $q$ es $d_q$. Las cargas se suman
en los productos tensoriales: si $u_k$ y $w_l$ son bases propias ortonormales de $Q_1$ y $Q_2$ con valores propios
$\alpha_k$ y $\beta_l$, la base ortonormal $u_k\otimes w_l$ consta de vectores propios de $Q_1\otimes1+1\otimes Q_2$ (valor
propio $\alpha_k+\beta_l$) y de $e^{i\theta Q_1}\otimes e^{i\theta Q_2}$ (valor propio $e^{i\theta\alpha_k}e^{i\theta\beta_l}$), luego,
por la definición del cálculo funcional,
$$
e^{i\theta(Q_1\otimes1+1\otimes Q_2)}=e^{i\theta Q_1}\otimes e^{i\theta Q_2},
$$
y actúa sobre la imagen de $P_{q_1}\otimes P_{q_2}$ con carga $q_1+q_2$. Esta es la regla de fusión $(q_1,q_2)\mapsto q_1+q_2$ de los
productos tensoriales fusionados del capítulo 7.

Para $n$ cúbits, la carga es el **operador número** $\hat N=\sum_jn_j$, con $n_j=|1\rangle\langle1|$ sobre el
cúbit $j$: $\hat N|b\rangle=|b|\,|b\rangle$ con el peso de Hamming $|b|$. Sus sectores tienen dimensión
$\binom nw$, de modo que los operadores que conservan el peso de Hamming forman
$$
U(1)'\cong\bigoplus_{w=0}^nM_{\binom nw}(\mathbb C),
$$
que es el producto fusionado de $n$ copias de $M_1\oplus M_1$ con cargas $0,1$
(bit $1$ = espín hacia arriba, $\hat N=S^z+\frac n2$). Todo hamiltoniano XXZ
$\sum_{\langle ij\rangle}J(s^x_is^x_j+s^y_is^y_j+\Delta s^z_is^z_j)$ conserva $\hat N$: $s^z_is^z_j$ es diagonal en la base
$|b\rangle$, y el término de salto $s^x_is^x_j+s^y_is^y_j=\frac12(s^+_is^-_j+s^-_is^+_j)$ (desarróllese
$s^\pm=s^x\pm is^y$; los términos cruzados se cancelan) lleva $|b\rangle$ a $0$ o a la configuración en la que un espín
hacia arriba ha pasado de uno de los sitios $i,j$ al otro, luego lleva cada subespacio propio de $\hat N$ en sí mismo.

*Cálculo:* la proposición 4.10 para cuatro espines: los sectores de carga hallados a partir únicamente de $\hat N$, y un anillo XXZ diagonalizado sector por sector — [tutorial, sección 4.4](../../tutorial/es/04_symmetries.md#s-4-4). *En la biblioteca:* $U(1)'$ para $n$ cúbits es `SpinChain(n).algebra`, y el hamiltoniano XXZ es `SpinChain.xxz` — [tutorial, sección 4.4](../../tutorial/es/04_symmetries.md#s-4-4).

> **Significado físico** (interpretación física, no un teorema de este curso). La simetría $U(1)$ es la conservación de un número de partículas, de una carga o de
> una magnetización. La proposición 4.10(3) es el teorema de Noether en este contexto: la carga se conserva si y
> solo si el hamiltoniano es simétrico, y entonces la dinámica nunca abandona un sector (teorema 8.3). Si además
> nada puede medir coherencias entre sectores, los sectores son sectores de superselección (capítulo 3).

<a id="s-4-5"></a>

## 4.5 $SU(2)$: el espín

**Definición.** Una **representación de espín** sobre $\mathbb C^d$ es una terna de matrices hermíticas
$S_x,S_y,S_z$ con
$$
[S_x,S_y]=iS_z,\qquad[S_y,S_z]=iS_x,\qquad[S_z,S_x]=iS_y,
$$
es decir, una representación de las etiquetas $\{x,y,z\}$ (sección 4.1) con estas relaciones. El ejemplo básico es
$s_a=\frac12\sigma_a$ sobre $\mathbb C^2$ (matrices de Pauli). Si $(S_a)$ y $(T_a)$ son representaciones de espín,
también lo es $S_a\otimes1+1\otimes T_a$ (los términos cruzados conmutan); por tanto, el **espín total**
$S_a=\sum_is_a^{(i)}$ de $L$ espines $\frac12$ (`decompose.total_spin(L)`) es una de ellas. Las rotaciones
$e^{-i\theta\,\mathbf n\cdot\mathbf S}$ ($\theta\in\mathbb R$, $\mathbf n\in\mathbb R^3$ un vector unitario, cálculo funcional)
conmutan con $X$ para todos los $\theta,\mathbf n$ si y solo si $X$ conmuta con las tres $S_a$ (derívese en $\theta=0$,
proposición 8.2(1); recíprocamente, teorema 1.5(d)), de modo que el conmutante de «la acción de $SU(2)$» es
$\{S_x,S_y,S_z\}'$, y por el teorema 2.2 las rotaciones y las $S_a$ generan la misma álgebra (ambos conjuntos son
cerrados por adjunción). Para $L$ espines $\frac12$, las rotaciones son las $V^{\otimes L}$, $V\in SU(2)$. En efecto,
$\mathbf n\cdot\mathbf S=\sum_ih^{(i)}$ con $h=\mathbf n\cdot\mathbf s$ y, como para las cargas en la sección 4.4
(inducción sobre $L$), $e^{-i\theta\,\mathbf n\cdot\mathbf S}=(e^{-i\theta h})^{\otimes L}$, con $e^{-i\theta h}\in SU(2)$ (sus
valores propios son $e^{\mp i\theta/2}$, pues $h$ tiene traza nula y $h^2=\frac14$). Recíprocamente, una $V\in SU(2)$ es
normal con valores propios de módulo $1$ y producto $1$, luego $V=W\mathrm{diag}(e^{i\varphi},e^{-i\varphi})W^*$ con
$W$ unitaria, y $V=e^{iA}$ para la matriz hermítica de traza nula $A=W\mathrm{diag}(\varphi,-\varphi)W^*$. Una matriz
hermítica $2\times2$ de traza nula es una combinación real $A=2\,\mathbf c\cdot\mathbf s$ de las matrices de Pauli, luego
$V=e^{-i\theta\,\mathbf n\cdot\mathbf s}$ con $\theta=2|\mathbf c|$, $\mathbf n=-\mathbf c/|\mathbf c|$ (y $\theta=0$ si
$\mathbf c=0$). El teorema B.6(a),(c) y el corolario B.7 (con $d=2$) muestran directamente que las $V^{\otimes L}$,
$V\in SU(2)$, generan linealmente el álgebra generada por las $S_a$.

Pongamos $S_\pm=S_x\pm iS_y$, de modo que $S_-=S_+^*$, y $\mathbf S^2=S_x^2+S_y^2+S_z^2$. Las relaciones dan
$$
[S_z,S_\pm]=\pm S_\pm,\qquad[S_+,S_-]=2S_z,\qquad\mathbf S^2=S_-S_++S_z^2+S_z,
$$
y $\mathbf S^2$ conmuta con todas las $S_a$ (por ejemplo,
$[\mathbf S^2,S_z]=S_x[S_x,S_z]+[S_x,S_z]S_x+S_y[S_y,S_z]+[S_y,S_z]S_y=-i\{S_x,S_y\}+i\{S_y,S_x\}=0$).

<a id="r-4-11"></a>

**Teorema 4.11 (representaciones de espín irreducibles).** Para todo $s\in\{0,\frac12,1,\frac32,\dots\}$, las
matrices de `decompose.spin_matrices(s)` —sobre $\mathbb C^{2s+1}$ con base $|m\rangle$, $m=s,s-1,\dots,-s$,
$$
S_z|m\rangle=m|m\rangle,\qquad S_+|m\rangle=\sqrt{s(s+1)-m(m+1)}\;|m+1\rangle,\qquad S_-=S_+^*,
$$
forman una representación de espín irreducible (**espín $s$**) con $\mathbf S^2=s(s+1)1$. Toda representación de
espín irreducible es unitariamente equivalente a exactamente una de ellas.

*Demostración.* *El espín $s$ es una representación de espín.* $[S_z,S_+]=S_+$ es claro.
$S_+S_-|m\rangle=(s(s+1)-(m-1)m)|m\rangle$ y $S_-S_+|m\rangle=(s(s+1)-m(m+1))|m\rangle$, luego $[S_+,S_-]=2S_z$; estas
dos relaciones y sus adjuntas equivalen a las tres relaciones de la definición.
$\mathbf S^2|m\rangle=(s(s+1)-m(m+1)+m^2+m)|m\rangle=s(s+1)|m\rangle$. *Irreducible:* un $W\neq0$ invariante es
invariante por $S_z$, cuyos valores propios son simples, luego contiene algún $|m\rangle$; los coeficientes de
$S_\pm$ solo se anulan en $m=\pm s$, así que aplicando $S_\pm$ se obtienen todos los $|m'\rangle$.

*Clasificación.* Sea $(S_a)$ irreducible sobre $\mathbb C^d$, $\lambda$ el mayor valor propio de $S_z$ y $v_0$ un
vector propio unitario. $S_zS_+v_0=S_+(S_z+1)v_0=(\lambda+1)S_+v_0$, luego $S_+v_0=0$ por maximalidad. Sea
$v_j=S_-^jv_0$; entonces $S_zv_j=(\lambda-j)v_j$ y, por inducción, usando $S_+S_-=S_-S_++2S_z$,
$$
S_+v_j=j(2\lambda-j+1)\,v_{j-1}
$$
($S_+v_j=S_-S_+v_{j-1}+2S_zv_{j-1}=[(j-1)(2\lambda-j+2)+2(\lambda-j+1)]v_{j-1}$). Los $v_j$ no nulos son vectores
propios de $S_z$ para valores propios distintos, luego hay un último no nulo, $v_r$. El subespacio generado por
$v_0,\dots,v_r$ es invariante por $S_z,S_\pm$, luego es igual a $\mathbb C^d$, y $d=r+1$. De
$0=S_+v_{r+1}=(r+1)(2\lambda-r)v_r$ se sigue $\lambda=r/2=:s$. Finalmente,
$\|v_j\|^2=\langle v_{j-1},S_+v_j\rangle=j(2s-j+1)\|v_{j-1}\|^2$, de modo que los vectores unitarios
$|s-j\rangle=v_j/\|v_j\|$ satisfacen $S_+|m\rangle=\sqrt{(s-m)(s+m+1)}\,|m+1\rangle$ con $m=s-j$, que es la fórmula de
arriba. Valores distintos de $s$ tienen dimensiones distintas. $\square$

*En la biblioteca:* el espín $s$ es `decompose.spin_matrices(s)`, comprobado frente al teorema 4.11 en [tutorial, sección 4.6](../../tutorial/es/04_symmetries.md#s-4-6); el espín total de $L$ espines $\frac12$ es `decompose.total_spin(L)` — [tutorial, sección 4.5](../../tutorial/es/04_symmetries.md#s-4-5).

<a id="r-4-12"></a>

**Corolario 4.12.** Toda representación de espín sobre $\mathbb C^n$ es unitariamente equivalente a
$\bigoplus_s(\text{espín }s)\otimes1_{m_s}$ y, con $d_M=\dim\ker(S_z-M)$,
$$
m_s=d_s-d_{s+1}.
$$

*Demostración.* Teorema 4.3 y teorema 4.11. Cada copia del espín $s'$ contiene exactamente un vector propio de
$S_z$ para el valor propio $M$ si $|M|\le s'$ y $s'-M\in\mathbb Z$, y ninguno en otro caso; así,
$d_M=\sum_{s'\ge|M|,\,s'-M\in\mathbb Z}m_{s'}$, y la diferencia de dos valores consecutivos es $m_s$. $\square$

<a id="r-4-13"></a>

**Proposición 4.13 ($L$ espines $\frac12$).** El espín total sobre $(\mathbb C^2)^{\otimes L}$ contiene el espín $S$
si y solo si $S\in\{\frac L2,\frac L2-1,\dots\}$, $S\ge0$, con multiplicidad
$$
m_S=\binom L{\frac L2-S}-\binom L{\frac L2-S-1}\qquad\Big(\tbinom L{-1}=0\Big),
$$
de modo que el álgebra generada por $S_x,S_y,S_z$ es $\bigoplus_SM_{2S+1}\otimes1_{m_S}$ y su conmutante
$\bigoplus_S1_{2S+1}\otimes M_{m_S}$.

*Demostración.* $S_z=\sum_is_z^{(i)}$ es diagonal en la base producto, con valor propio $u-\frac L2$ sobre un producto
de $u$ espines hacia arriba y $L-u$ espines hacia abajo, luego $d_M=\binom L{L/2+M}=\binom L{L/2-M}$ si
$\frac L2+M\in\{0,\dots,L\}$ y $0$ en otro caso. Se aplican el corolario 4.12 y el teorema 4.3. Los valores son
positivos porque $\binom L{j+1}/\binom Lj=\frac{L-j}{j+1}>1$ para $j<\frac{L-1}2$. $\square$

Para $L=2$: $\frac12\otimes\frac12=1\oplus0$ (**Clebsch–Gordan**): el triplete $S_-^j|{\uparrow\uparrow}\rangle$,
$j=0,1,2$, simétrico bajo el intercambio $F$ de los dos espines, y el singlete
$\frac1{\sqrt2}(|{\uparrow\downarrow}\rangle-|{\downarrow\uparrow}\rangle)$, antisimétrico y aniquilado por $S_+$ y $S_z$.
Para $L=4$: $m_2=1$, $m_1=4-1=3$, $m_0=6-4=2$, luego
$$
\{S_x,S_y,S_z\}''\cong M_5\oplus M_3\oplus M_1\ \text{ with multiplicities }1,3,2,\qquad
\{S_x,S_y,S_z\}'\cong M_1\oplus M_3\oplus M_2,
$$
de dimensiones $35$ y $14$, los números hallados en el capítulo 10.

*Cálculo:* la descomposición para $L=4$ a partir únicamente de los tres generadores, con la proposición 4.13 y $\mathbf S^2=S(S+1)$ sobre el espín $S$ — [tutorial, sección 4.5](../../tutorial/es/04_symmetries.md#s-4-5).

¿Qué operadores conmutan con todas las rotaciones? Las permutaciones de los sitios lo hacen: la conjugación por
$\pi(\sigma)$ lleva $s^{(i)}_a$ a $s^{(\sigma(i))}_a$ (lema B.1(b)), luego deja fijo $\sum_is^{(i)}_a$. Recíprocamente,
ningún otro operador lo hace:

<a id="r-4-14"></a>

**Teorema 4.14 (dualidad de Schur–Weyl).** Sea $\pi$ la representación de $S_L$ que permuta los factores de
$(\mathbb C^2)^{\otimes L}$. (a) $\mathrm{span}\,\pi(S_L)\subseteq\{S_x,S_y,S_z\}'$. (b) Se da la igualdad;
equivalentemente (teorema 2.2), el álgebra generada por el espín total es $\pi(S_L)'$.

*Demostración.* La parte (a) se acaba de demostrar. La parte (b) es el corolario B.7, el caso $d=2$ del teorema B.6,
que demuestra la dualidad para $(\mathbb C^d)^{\otimes L}$ y la acción diagonal de $U(d)$ para todos los $d$ y $L$.
$\square$

(Para profundizar: W. Fulton y J. Harris, *Representation Theory*, Springer 1991, lección 6; R. Goodman y
N. Wallach, *Symmetry, Representations, and Invariants*, Springer 2009, sección 9.1.)

*Cálculo:* el teorema 4.14 para $L\le5$, las dimensiones $2,5,14,42$ de $\{S_a\}'$ (corolario B.10), y $\pi(S_4)''$, cuyos pares $(k_c,m_c)$ son los del espín total con los papeles intercambiados (corolario B.9) — [tutorial, sección 4.5](../../tutorial/es/04_symmetries.md#s-4-5).

> **Significado físico** (interpretación física, no un teorema de este curso). Un multiplete de espín $S$ consta de $2S+1$ estados de igual energía para todo
> hamiltoniano invariante por rotaciones, que se distinguen por $S^z=M$. Dos espines $\frac12$ forman un triplete y
> un singlete; cuatro espines forman un quintuplete, tres tripletes y dos singletes. Los espacios de multiplicidad
> de los dos singletes y de los tres tripletes son el cúbit libre de decoherencia y el cútrit sin ruido del
> capítulo 10.

<a id="s-4-6"></a>

## 4.6 Sistemas de espines

El **modelo de Heisenberg** sobre $L$ espines $\frac12$ es $H=\sum_{i<j}J_{ij}\,\mathbf s^{(i)}\cdot\mathbf s^{(j)}$ con
acoplamientos reales $J_{ij}$.

<a id="r-4-15"></a>

**Proposición 4.15.** (a) (*Identidad de Dirac*) $\mathbf s^{(i)}\cdot\mathbf s^{(j)}=\frac12F_{ij}-\frac14$, donde $F_{ij}$
intercambia los espines $i$ y $j$. Por tanto, todo hamiltoniano de Heisenberg pertenece a
$\mathrm{span}\,\pi(S_L)\subseteq\{S_a\}'$, y $V^*HV=\bigoplus_S1_{2S+1}\otimes H_S$ con bloques $H_S\in M_{m_S}$ en la
descomposición de la proposición 4.13: cada valor propio de $H_S$ es un nivel $(2S+1)$ veces degenerado de $H$.
(b) El espacio de los **vectores de peso máximo** $\ker S_+\cap\ker(S_z-S)$ tiene dimensión $m_S$, es invariante
por $H$, y $H$ actúa sobre él como $H_S$ (salvo un cambio de base unitario).

*Demostración.* (a) Sobre dos espines,
$\mathbf s^{(1)}\cdot\mathbf s^{(2)}=\frac12(\mathbf S^2-\mathbf s^{(1)2}-\mathbf s^{(2)2})=\frac12\mathbf S^2-\frac34$
(teorema 4.11 con $s=\frac12$). Los vectores del triplete $S_-^j|{\uparrow\uparrow}\rangle$ son $|{\uparrow\uparrow}\rangle$,
$|{\downarrow\uparrow}\rangle+|{\uparrow\downarrow}\rangle$, $2|{\downarrow\downarrow}\rangle$, y junto con el singlete forman una base
ortogonal de $\mathbb C^2\otimes\mathbb C^2$. Sobre el triplete, $\mathbf S^2=2$ y $F=1$: $\mathbf S^2=S_-S_++S_z^2+S_z$ vale
$0+1+1=2$ sobre $|{\uparrow\uparrow}\rangle$, y $\mathbf S^2$ y $F$ conmutan con $S_-$, mientras que $F$ deja fijo
$|{\uparrow\uparrow}\rangle$. Sobre el singlete, $\mathbf S^2=0$ (está aniquilado por $S_+$ y $S_z$) y $F=-1$. Ambos miembros
coinciden sobre estos cuatro vectores. Para sitios $i,j$ de $L$, se toma el producto tensorial con la identidad. El
resto es el teorema 4.14(a) y el teorema 4.3. (b) Demostramos más: para cualquier representación de espín
descompuesta como en el corolario 4.12 y cualquier $\mu\in\frac12\mathbb Z$, $\ker S_+\cap\ker(S_z-\mu)$ es (la imagen
por la unitaria del corolario 4.12 de) $|\mu,\mu\rangle\otimes\mathbb C^{m_\mu}$ si $\mu\ge0$, y $0$ si $\mu<0$. En efecto,
en esa base un vector con $S_z=\mu$ es $\sum_{S'}|S',\mu\rangle\otimes w_{S'}$ sobre los $S'\ge|\mu|$ con
$S'-\mu\in\mathbb Z$, y por el teorema 4.11 $S_+$ lo lleva a
$\sum_{S'}\sqrt{(S'-\mu)(S'+\mu+1)}\,|S',\mu+1\rangle\otimes w_{S'}$, donde $S'+\mu+1\ge1$; luego se anula si y solo si
$w_{S'}=0$ para todo $S'\neq\mu$. Para el $H$ de Heisenberg, $1\otimes H_S$ actúa sobre $|S,S\rangle\otimes\mathbb C^{m_S}$
como $H_S$. $\square$

*Cálculo:* la proposición 4.15 sobre un anillo de cuatro espines y sobre un anillo de ocho espines — [tutorial, sección 4.6](../../tutorial/es/04_symmetries.md#s-4-6).

La parte (b) es la forma en que el módulo `su2` (cadenas de espines, capítulo 9) construye $H_S$ sin
descomposición numérica: construye una base ortonormal de los vectores de peso máximo, añadiendo un sitio cada vez,
con los coeficientes de Clebsch–Gordan de $j\otimes\frac12$. Los deducimos descendiendo desde el peso máximo. Por el
teorema 4.11 y $S_-=S_+^*$, el espín $j$ satisface
$$
S_\pm|j,m\rangle=\sqrt{(j\mp m)(j\pm m+1)}\;|j,m\pm1\rangle ,
$$
y ponemos $|j,m\rangle=0$ para $|m|>j$ (coherente con la fórmula, cuyo coeficiente se anula en $m=\pm j$).

<a id="r-4-16"></a>

**Teorema 4.16 (coeficientes de Clebsch–Gordan de $j\otimes\frac12$).** Sea $j\in\{\frac12,1,\frac32,\dots\}$, sea
$S_a=S^{(j)}_a\otimes1+1\otimes s_a$ sobre $\mathbb C^{2j+1}\otimes\mathbb C^2$ (el espín $j$ del teorema 4.11, con base
$|j,m\rangle$, acoplado a un espín $\frac12$ con $|{\uparrow}\rangle=|\frac12,\frac12\rangle$,
$|{\downarrow}\rangle=|\frac12,-\frac12\rangle$), y pongamos $\alpha_M=\sqrt{(j+M+\frac12)/(2j+1)}$,
$\beta_M=\sqrt{(j-M+\frac12)/(2j+1)}$. Para $J=j\pm\frac12$ y $M=J,J-1,\dots,-J$ sea
$$
\begin{aligned}
|j+\tfrac12,M\rangle\!\rangle&=\alpha_M\,|j,M-\tfrac12\rangle\otimes|{\uparrow}\rangle+\beta_M\,|j,M+\tfrac12\rangle\otimes|{\downarrow}\rangle,\\
|j-\tfrac12,M\rangle\!\rangle&=-\beta_M\,|j,M-\tfrac12\rangle\otimes|{\uparrow}\rangle+\alpha_M\,|j,M+\tfrac12\rangle\otimes|{\downarrow}\rangle .
\end{aligned}
$$
(a) Estos $2(2j+1)$ vectores forman una base ortonormal de $\mathbb C^{2j+1}\otimes\mathbb C^2$.
(b) $|j+\frac12,j+\frac12\rangle\!\rangle=|j,j\rangle\otimes|{\uparrow}\rangle$, y $|j-\frac12,j-\frac12\rangle\!\rangle$ es, salvo un factor de
módulo uno, el único vector unitario $x$ con $S_zx=(j-\frac12)x$ y $S_+x=0$; el factor se fija mediante un coeficiente
positivo de $|j,j\rangle\otimes|{\downarrow}\rangle$ (la convención de Condon–Shortley).
(c) $S_z|J,M\rangle\!\rangle=M|J,M\rangle\!\rangle$ y $S_-|J,M\rangle\!\rangle=\sqrt{(J+M)(J-M+1)}\,|J,M-1\rangle\!\rangle$: los estados
$|J,M\rangle\!\rangle$ se obtienen de $|J,J\rangle\!\rangle$ por descenso y normalización.
(d) Por tanto, $S_a|J,M\rangle\!\rangle=\sum_{M'}(S^{(J)}_a)_{M'M}|J,M'\rangle\!\rangle$ para $a=x,y,z$: la unitaria
$|J,M\rangle\mapsto|J,M\rangle\!\rangle$ lleva espín $(j+\frac12)\oplus$ espín $(j-\frac12)$ sobre $S_a$, luego
$j\otimes\frac12\cong(j+\frac12)\oplus(j-\frac12)$, y los coeficientes de Clebsch–Gordan
$\langle j,m;\frac12,\sigma|J,M\rangle$ (las coordenadas de $|J,M\rangle\!\rangle$ en la base producto) son $\pm\alpha_M,\beta_M$ como
se muestra.

*Demostración.* (a) Si un vector de la base que aparece en la fórmula es $0$, también lo es su coeficiente:
$|j,M-\frac12\rangle=0$ solo para $M=-j-\frac12$, donde $\alpha_M=0$, y $|j,M+\frac12\rangle=0$ solo para $M=j+\frac12$,
donde $\beta_M=0$. Así, cada vector tiene norma $\alpha_M^2+\beta_M^2=1$. Vectores con $M$ distinto involucran vectores
distintos de la base producto, y para $M$ igual el producto escalar es $-\alpha_M\beta_M+\beta_M\alpha_M=0$. Hay
$(2j+2)+2j=2(2j+1)$ de ellos.
(b) La primera afirmación es $\alpha_{j+1/2}=1$, $\beta_{j+1/2}=0$. El subespacio propio $S_z=j-\frac12$ está generado por
$|j,j-1\rangle\otimes|{\uparrow}\rangle$ y $|j,j\rangle\otimes|{\downarrow}\rangle$ y, con $s_+|{\uparrow}\rangle=0$,
$s_+|{\downarrow}\rangle=|{\uparrow}\rangle$,
$$
S_+\big(p\,|j,j-1\rangle\otimes|{\uparrow}\rangle+q\,|j,j\rangle\otimes|{\downarrow}\rangle\big)=(p\sqrt{2j}+q)\,|j,j\rangle\otimes|{\uparrow}\rangle ,
$$
que se anula si y solo si $q=-\sqrt{2j}\,p$; con $|p|^2+|q|^2=1$ y $q>0$ (luego $p<0$) esto es
$p=-\beta_{j-1/2}=-1/\sqrt{2j+1}$, $q=\alpha_{j-1/2}=\sqrt{2j/(2j+1)}$. (c) $S_z$: ambos términos tienen
$S_z=(M\mp\frac12)\pm\frac12=M$. Para $S_-$, sea
$x_M=p\,|j,M-\frac12\rangle\otimes|{\uparrow}\rangle+q\,|j,M+\frac12\rangle\otimes|{\downarrow}\rangle$; como
$s_-|{\uparrow}\rangle=|{\downarrow}\rangle$ y $s_-|{\downarrow}\rangle=0$,
$$
S_-x_M=p\sqrt{(j+M-\tfrac12)(j-M+\tfrac32)}\;|j,M-\tfrac32\rangle\otimes|{\uparrow}\rangle
+\Big(p+q\sqrt{(j+M+\tfrac12)(j-M+\tfrac12)}\Big)|j,M-\tfrac12\rangle\otimes|{\downarrow}\rangle .
$$
Usamos $\beta_M\sqrt{(j+M+\frac12)(j-M+\frac12)}=(j-M+\frac12)\alpha_M$ y
$\alpha_M\sqrt{(j+M+\frac12)(j-M+\frac12)}=(j+M+\frac12)\beta_M$.
Para $J=j+\frac12$, $(p,q)=(\alpha_M,\beta_M)$ y $c=\sqrt{(J+M)(J-M+1)}=\sqrt{(j+M+\frac12)(j-M+\frac32)}$: el primer
coeficiente es $\sqrt{(j+M+\frac12)(j+M-\frac12)(j-M+\frac32)/(2j+1)}=c\,\alpha_{M-1}$, el segundo
$\alpha_M(j-M+\frac32)=c\,\beta_{M-1}$.
Para $J=j-\frac12$, $(p,q)=(-\beta_M,\alpha_M)$ y $c=\sqrt{(j+M-\frac12)(j-M+\frac12)}$: el primer coeficiente es
$-\sqrt{(j-M+\frac12)(j+M-\frac12)(j-M+\frac32)/(2j+1)}=-c\,\beta_{M-1}$, el segundo $\beta_M(j+M-\frac12)=c\,\alpha_{M-1}$.
En ambos casos $S_-|J,M\rangle\!\rangle=c\,|J,M-1\rangle\!\rangle$ para $M>-J$. Para $M=-J$ ambos coeficientes se anulan
($p=\alpha_M=0$ si $J=j+\frac12$; $j+M-\frac12=0$ si $J=j-\frac12$), luego $S_-|J,-J\rangle\!\rangle=0$. (d) Por (a) y (c), en
la base ortonormal $|J,M\rangle\!\rangle$ las matrices de $S_z$ y $S_-$ son las de espín $(j+\frac12)\oplus$ espín
$(j-\frac12)$ (teorema 4.11 y la fórmula de arriba); lo mismo vale para la matriz de $S_+=S_-^*$, y por tanto para las
de $S_x=\frac12(S_++S_-)$ y $S_y=\frac1{2i}(S_+-S_-)$. $\square$

*Cálculo:* el teorema 4.16 para $j\le2$, frente a matrices densas y frente a `su2.clebsch_gordan` — [tutorial, sección 4.6](../../tutorial/es/04_symmetries.md#s-4-6).

El ejercicio 2 extiende la descomposición (pero no los coeficientes) a $j_1\otimes j_2$. Para la biblioteca solo
importan los dos vectores de peso máximo del teorema 4.16, y funcionan en cualquier representación de espín,
irreducible o no:

<a id="r-4-17"></a>

**Corolario 4.17 (vectores de peso máximo sitio a sitio).** Sea $(T_a)$ una representación de espín sobre $E$, sea
$\mathcal H_\mu(E)=\ker T_+\cap\ker(T_z-\mu)$, con una base ortonormal $B_\mu$ ($B_\mu=\emptyset$ si $\mathcal H_\mu(E)=0$),
y sea $S_a=T_a\otimes1+1\otimes s_a$ sobre $E\otimes\mathbb C^2$. Para todo $S\in\{0,\frac12,1,\dots\}$, los vectores
$$
v\otimes|{\uparrow}\rangle\quad(v\in B_{S-1/2}),\qquad
\sqrt{\tfrac{2S+1}{2S+2}}\;v\otimes|{\downarrow}\rangle-\tfrac1{\sqrt{(2S+1)(2S+2)}}\;T_-v\otimes|{\uparrow}\rangle\quad(v\in B_{S+1/2})
$$
forman una base ortonormal de $\ker S_+\cap\ker(S_z-S)$ (la primera familia es vacía para $S=0$). Partiendo de un
espín ($E=\mathbb C^2$, $B_{1/2}=\{|{\uparrow}\rangle\}$, todos los demás $B_\mu$ vacíos) y añadiendo $L-1$ espines, esto da
bases ortonormales de los espacios de vectores de peso máximo de $L$ espines $\frac12$, de dimensiones $m_S$
(proposición 4.15(b)).

*Demostración.* Escribimos $a=\sqrt{(2S+1)/(2S+2)}$, $b=-1/\sqrt{(2S+1)(2S+2)}$,
$u(v)=a\,v\otimes|{\downarrow}\rangle+b\,T_-v\otimes|{\uparrow}\rangle$. Para $v\in\mathcal H_{S+1/2}(E)$,
$T_+T_-v=(T_-T_++2T_z)v=(2S+1)v$, y $a+(2S+1)b=0$. *Peso máximo:* $S_z$ vale $S$ sobre ambas familias
($T_z=S\mp\frac12$ sobre $v$ y $S-\frac12$ sobre $T_-v$); $S_+(v\otimes|{\uparrow}\rangle)=T_+v\otimes|{\uparrow}\rangle=0$, y
$S_+u(v)=a\,T_+v\otimes|{\downarrow}\rangle+a\,v\otimes|{\uparrow}\rangle+b\,T_+T_-v\otimes|{\uparrow}\rangle=(a+(2S+1)b)\,v\otimes|{\uparrow}\rangle=0$.
*Ortonormalidad:* $\langle v\otimes|{\uparrow}\rangle,u(v')\rangle=b\langle v,T_-v'\rangle=b\langle T_+v,v'\rangle=0$ para
$v\in B_{S-1/2}$, y
$\langle u(v),u(v')\rangle=a^2\langle v,v'\rangle+b^2\langle v,T_+T_-v'\rangle=(a^2+(2S+1)b^2)\langle v,v'\rangle=\langle v,v'\rangle$.
*Generación:* sea $S_zw=Sw$ y $S_+w=0$, y escribamos $w=x\otimes|{\uparrow}\rangle+y\otimes|{\downarrow}\rangle$. Entonces
$T_zx=(S-\frac12)x$, $T_zy=(S+\frac12)y$ y $0=S_+w=(T_+x+y)\otimes|{\uparrow}\rangle+T_+y\otimes|{\downarrow}\rangle$, luego $T_+y=0$,
es decir, $y\in\mathcal H_{S+1/2}(E)$, y $T_+x=-y$. Entonces $x'=x+\frac1{2S+1}T_-y$ satisface $T_+x'=-y+y=0$, luego
$x'\in\mathcal H_{S-1/2}(E)$, que es $0$ para $S=0$ (demostración de la proposición 4.15(b)). Como $b/a=-1/(2S+1)$,
$w=x'\otimes|{\uparrow}\rangle+a^{-1}u(y)$, con $u$ extendida linealmente, está en el subespacio generado por las dos
familias. Para $L$ espines se aplica esto con $E=(\mathbb C^2)^{\otimes(n-1)}$ y su espín total $T_a$, $n=2,\dots,L$;
sobre un espín, $\mathcal H_{1/2}=\mathbb C|{\uparrow}\rangle$ y todos los demás $\mathcal H_\mu=0$ (el único espín es
$\frac12$). $\square$

La segunda familia consta de los vectores $|S,S\rangle\!\rangle$ del teorema 4.16 con $j=S+\frac12$, escritos dentro de la
copia del espín $j$ generada por los $T_-^kv$ (demostración del teorema 4.11, $|j,j-1\rangle=T_-v/\sqrt{2j}$), y la primera
familia es $|j+\frac12,j+\frac12\rangle\!\rangle$ con $j=S-\frac12$. Etiquetar cada vector de la base por su sucesión de espines
intermedios $S_1=\frac12,S_2,\dots,S_L=S$, $|S_n-S_{n-1}|=\frac12$, da la base de Yamanouchi.

*En la biblioteca:* `su2.SU2Chain` construye la base de Yamanouchi con estos coeficientes, dentro del sector $U(1)$ $\hat N=\frac L2+S$ de `SpinChain`; el corolario 4.17 se comprueba para ocho espines, con los bloques $U(1)$ y $SU(2)$ del anillo de Heisenberg y su coste — [tutorial, sección 4.6](../../tutorial/es/04_symmetries.md#s-4-6).

> **Significado físico** (interpretación física, no un teorema de este curso). El modelo de Heisenberg describe la
> interacción de intercambio de espines electrónicos localizados en aislantes magnéticos; $J_{ij}>0$ favorece vecinos
> antiparalelos. Su espectro consta de multipletes, y un campo magnético $h\,S^z$ desdobla cada multiplete en $2S+1$
> niveles equidistantes (ejercicio 4). El teorema 9.9 (acoplamientos bipartitos; para cadenas y anillos
> antiferromagnéticos de longitud par, un estado fundamental singlete) se remonta a E. Lieb y D. Mattis
> (*J. Math. Phys.* 3 (1962) 749).

<a id="s-4-7"></a>

## 4.7 Álgebras reales: tipo real, complejo y cuaterniónico

Las simetrías vienen dadas a menudo por matrices reales: las rotaciones de $\mathbb R^3$, las permutaciones,
hamiltonianos que son reales en alguna base. El álgebra compleja que generan se descompone por el teorema 2.6, pero
las matrices *reales* que contiene forman un álgebra real, cuyas piezas simples son de tres clases, $M_k(\mathbb R)$,
$M_k(\mathbb C)$ y $M_k(\mathbb H)$. Las deducimos de la descomposición compleja y de un dato más: lo que la
conjugación compleja hace con los sectores.

Para $a\in M_n(\mathbb C)$, sea $\bar a$ la conjugada compleja entrada a entrada (análogamente para vectores, y
$\bar S=\{\bar s:s\in S\}$ para conjuntos). La conjugación es antilineal y biyectiva, $\overline{ab}=\bar a\bar b$,
$\overline{a^*}=\bar a^*$ y $\bar{\bar a}=a$; $\bar a=a$ si y solo si $a$ es real, y entonces $a^*=a^T$. Para un
conjunto $\mathcal B$ de matrices (o de vectores), sea $\mathcal B_{\mathbb R}=\mathcal B\cap M_n(\mathbb R)$ (o
$\mathcal B\cap\mathbb R^n$). El producto escalar es $\langle u,v\rangle=\sum_i\bar u_iv_i$, de modo que
$\langle\bar u,\bar v\rangle=\overline{\langle u,v\rangle}=\langle v,u\rangle$.

<a id="r-4-18"></a>

**Lema 4.18 (partes reales).** (a) Sea $\mathcal B$ un subespacio complejo de $M_n(\mathbb C)$ o de $\mathbb C^n$ con
$\bar{\mathcal B}=\mathcal B$. Entonces $\mathcal B=\mathcal B_{\mathbb R}\oplus i\mathcal B_{\mathbb R}$ y
$\dim_{\mathbb R}\mathcal B_{\mathbb R}=\dim_{\mathbb C}\mathcal B$.
(b) Sea $\mathcal A$ una $*$-álgebra sobre $\mathbb C^n$ con $\bar{\mathcal A}=\mathcal A$. Entonces
$\mathcal A_{\mathbb R}$ es un álgebra real con unidad, cerrada para la transposición;
$\overline{\mathcal A'}=\mathcal A'$; y $(\mathcal A')_{\mathbb R}$ es el conjunto de las $X$ reales con $Xa=aX$ para
todo $a\in\mathcal A_{\mathbb R}$.
(c) Si $S\subseteq M_n(\mathbb C)$ satisface $S^*=S$ y $\bar S=S$ —por ejemplo,
$S=\{S_1,\dots,S_J,S_1^T,\dots,S_J^T\}$ con $S_j$ reales, o un grupo de matrices ortogonales reales—, entonces
$\overline{S''}=S''$.

*Demostración.* (a) Para $b\in\mathcal B$, $\mathrm{Re}\,b=\frac12(b+\bar b)$ y
$\mathrm{Im}\,b=\frac1{2i}(b-\bar b)$ son elementos reales de $\mathcal B$ con
$b=\mathrm{Re}\,b+i\mathrm{Im}\,b$, y un elemento de $\mathcal B_{\mathbb R}\cap i\mathcal B_{\mathbb R}$ es
real y puramente imaginario, luego $0$. Por tanto, una base $b_1,\dots,b_r$ del espacio real $\mathcal B_{\mathbb R}$
genera $\mathcal B$ sobre $\mathbb C$, y es linealmente independiente sobre $\mathbb C$: si
$\sum_j(\alpha_j+i\beta_j)b_j=0$ con $\alpha_j,\beta_j$ reales, las partes real e imaginaria dan
$\sum_j\alpha_jb_j=\sum_j\beta_jb_j=0$. (b) $\mathcal A_{\mathbb R}$ contiene a $1$ y es cerrada para productos, y
para $a\mapsto a^T=a^*$. $Y$ conmuta con todos los $t\in T$ si y solo si $\bar Y$ conmuta con todos los $\bar t$, así
que $\overline{T'}=\bar T'$ para todo conjunto $T$, y $\overline{\mathcal A'}=\bar{\mathcal A}'=\mathcal A'$. Una $X$
real que conmuta con $\mathcal A_{\mathbb R}$ conmuta con $\mathcal A_{\mathbb R}+i\mathcal A_{\mathbb R}=\mathcal A$
(por (a)). (c) Aplíquese $\overline{T'}=\bar T'$ dos veces: $\overline{S''}=\bar S''=S''$. $\square$

**Notación fija.** En el resto de esta sección, $\mathcal A$ es una $*$-álgebra sobre $\mathbb C^n$ con
$\bar{\mathcal A}=\mathcal A$, y $V$ es una descomposición como en el teorema 2.6. $V_c$ es el bloque de columnas
$n\times k_cm_c$ de $V$ que corresponde al $c$-ésimo sumando, de modo que $V_c^*V_d=\delta_{cd}1$ y $V_cV_c^*=P_c$,
las proyecciones centrales minimales, y
$$
aV_c=V_c\big(\pi_c(a)\otimes1_{m_c}\big)\quad(a\in\mathcal A),\qquad XV_c=V_c\big(1_{k_c}\otimes x_c\big)\quad(X\in\mathcal A'),
$$
con $\pi_c(a)\in M_{k_c}$ y $x_c\in M_{m_c}$. Cada $\pi_c$ es un $*$-homomorfismo unital de $\mathcal A$ sobre
$M_{k_c}$: para $x\in M_{k_c}$, el elemento $V_c(x\otimes1)V_c^*$ de $\mathcal A$ tiene $\pi_c=x$ y $\pi_e=0$ para
$e\neq c$. Por último, $\Omega_2=\begin{pmatrix}0&-1\\1&0\end{pmatrix}$ y $\Omega_{2r}=1_r\otimes\Omega_2$ (diagonal
por bloques con $r$ copias de $\Omega_2$).

Una aplicación $\theta:E\to E$ sobre un subespacio $E\subseteq\mathbb C^k$ es **antiunitaria** si
$\theta(\alpha u+\beta v)=\bar\alpha\theta u+\bar\beta\theta v$ y
$\langle\theta u,\theta v\rangle=\langle v,u\rangle$. Para una unitaria $J\in U(k)$, $\theta_Jv=J\bar v$ es
antiunitaria sobre $\mathbb C^k$ ($\langle J\bar u,J\bar v\rangle=\langle\bar u,\bar v\rangle=\langle v,u\rangle$),
con $\theta_J^2v=J\bar Jv$.

<a id="r-4-19"></a>

**Lema 4.19 (aplicaciones antiunitarias).** Sea $\theta$ antiunitaria sobre $E\neq0$ con
$\theta^2=\varepsilon\,\mathrm{id}$, $\varepsilon\in\{1,-1\}$.
(a) Si $\varepsilon=1$, $E$ tiene una base ortonormal de vectores $w$ con $\theta w=w$.
(b) Si $\varepsilon=-1$, entonces $\langle u,\theta u\rangle=0$ para todo $u\in E$; para todo vector unitario $u$,
$\mathrm{span}\{u,\theta u\}$ es un plano $\theta$-invariante; y $E$ tiene una base ortonormal
$w_1,\theta w_1,\dots,w_r,\theta w_r$. En particular, $\dim E=2r$ es par.
(c) Sea $J\in U(k)$ con $J\bar J=\varepsilon1$. Si $\varepsilon=1$, existe una unitaria $W$ con $J\bar W=W$, es decir,
$J=WW^T$. Si $\varepsilon=-1$, $k$ es par y existe una unitaria $W$ con $J\bar W=W\Omega_k$, es decir,
$J=W\Omega_kW^T$.

*Demostración.* (a) Sea $F=\{v\in E:\theta v=v\}$, un subespacio real. Todo $v\in E$ es $v_1+iv_2$ con
$v_1=\frac12(v+\theta v)$ y $v_2=\frac1{2i}(v-\theta v)$, y ambos están en $F$
($\theta v_2=\frac1{-2i}(\theta v-v)=v_2$). Para $u,v\in F$,
$\langle u,v\rangle=\langle\theta u,\theta v\rangle=\langle v,u\rangle$ es real, así que $\langle\cdot,\cdot\rangle$
es un producto escalar real sobre $F$; sea $w_1,\dots,w_s$ una base ortonormal de $F$ (Gram–Schmidt). Los $w_j$ son
ortonormales en $E$, luego linealmente independientes sobre $\mathbb C$, y generan $E$ sobre $\mathbb C$ porque
$E=F+iF$. (b) Aplicando $\langle\theta x,\theta y\rangle=\langle y,x\rangle$ a $x=\theta u$, $y=u$:
$\langle u,\theta u\rangle=\langle\theta^2u,\theta u\rangle=-\langle u,\theta u\rangle$. Así, para un vector unitario
$u$, $\theta u$ es un vector unitario ortogonal a $u$, y $E_1=\mathrm{span}\{u,\theta u\}$ es invariante:
$\theta(\alpha u+\beta\theta u)=\bar\alpha\theta u-\bar\beta u$. El complemento ortogonal de $E_1$ en $E$ también es
invariante: para $x\perp E_1$ e $y\in E_1$, $y=-\theta(\theta y)$ y
$\langle\theta x,y\rangle=-\langle\theta x,\theta(\theta y)\rangle=-\langle\theta y,x\rangle=0$, pues
$\theta y\in E_1$. La inducción sobre $\dim E$ da la base. (c) $\theta_J$ es antiunitaria sobre $\mathbb C^k$ con
$\theta_J^2=J\bar J=\varepsilon$. Sea $W$ la unitaria cuyas columnas son la base de (a), respectivamente la base
$w_1,\theta_Jw_1,w_2,\theta_Jw_2,\dots$ de (b). Las columnas de $J\bar W$ son las imágenes por $\theta_J$ de las
columnas de $W$. En el caso (a) esto dice $J\bar W=W$. En el caso (b) las imágenes son $\theta_Jw_j$ y
$\theta_J^2w_j=-w_j$, que son las columnas $2j-1$ y $2j$ de $W\Omega_k$ (pues $\Omega_2e_1=e_2$, $\Omega_2e_2=-e_1$).
Como $\bar W^{-1}=\overline{W^*}=W^T$, $J=WW^T$, respectivamente $J=W\Omega_kW^T$. $\square$

<a id="r-4-20"></a>

**Teorema 4.20 (conjugación de los sectores).** (a) Existen una permutación $\sigma$ de los sectores con
$\sigma\circ\sigma=\mathrm{id}$ y unitarias $J_c\in U(k_c)$, $M_c\in U(m_c)$ tales que $\bar P_c=P_{\sigma(c)}$,
$k_{\sigma(c)}=k_c$, $m_{\sigma(c)}=m_c$ y, para todo $c$ y todo $a\in\mathcal A$,
$$
\overline{V_c}=V_{\sigma(c)}\,(J_c\otimes M_c),\qquad\pi_{\sigma(c)}(\bar a)=J_c\,\overline{\pi_c(a)}\,J_c^*.
$$
$J_c\otimes M_c$ es única, y $J_c,M_c$ son únicas salvo $(J_c,M_c)\mapsto(\lambda J_c,\bar\lambda M_c)$ con
$|\lambda|=1$.
(b) Si $\sigma(c)=c$, existe $\varepsilon_c\in\{1,-1\}$ con $J_c\bar J_c=\varepsilon_c1_{k_c}$ y
$M_c\bar M_c=\varepsilon_c1_{m_c}$. Si $\varepsilon_c=-1$, entonces $k_c$ y $m_c$ son pares.
(c) $\sigma$ y los $\varepsilon_c$ no dependen de la elección de $V$.
(d) Para $y\in M_{m_c}$,
$\overline{V_c(1\otimes y)V_c^*}=V_{\sigma(c)}\big(1\otimes M_c\bar yM_c^*\big)V_{\sigma(c)}^*$. Por tanto, si
$X\in\mathcal A'$ tiene los bloques $x_c$, entonces $\bar X\in\mathcal A'$ tiene los bloques
$(\bar X)_{\sigma(c)}=M_c\bar x_cM_c^*$: sobre el conmutante, $M_c$ desempeña el papel que $J_c$ desempeña sobre
$\mathcal A$.

*Demostración.* (a) La conjugación aplica $\mathcal A$ sobre $\mathcal A$ y conserva sumas, productos y adjuntos, así
que aplica el centro $Z(\mathcal A)$ sobre sí mismo, y el conjunto de las proyecciones no nulas de $Z(\mathcal A)$
biyectivamente sobre sí mismo, conservando el orden $p\le q$ ($pq=p$). Por tanto, permuta las minimales, que son las
$P_c$ (lema 2.3): $\bar P_c=P_{\sigma(c)}$, y $\sigma\circ\sigma=\mathrm{id}$ porque $\bar{\bar P}_c=P_c$. Sea
$d=\sigma(c)$. La imagen de $\overline{V_c}$ es la conjugada de la imagen de $P_c$, que es la imagen de
$\bar P_c=P_d$; así que $\overline{V_c}=P_d\overline{V_c}=V_dT$ con $T=V_d^*\overline{V_c}$, y
$T^*T=\overline{V_c}^*P_d\overline{V_c}=\overline{V_c^*V_c}=1$. Como
$k_cm_c=\mathrm{rank}\,P_c=\mathrm{rank}\,\bar P_c=k_dm_d$, la isometría $T$ es cuadrada, luego unitaria.
Conjugando $aV_c=V_c(\pi_c(a)\otimes1)$ y usando $\bar aV_d=V_d(\pi_d(\bar a)\otimes1)$ se obtiene
$V_d(\pi_d(\bar a)\otimes1)T=V_dT(\overline{\pi_c(a)}\otimes1)$; multiplicando por $V_d^*$,
$$
\big(\pi_d(\bar a)\otimes1_{m_d}\big)\,T=T\,\big(\overline{\pi_c(a)}\otimes1_{m_c}\big)\qquad(a\in\mathcal A).
$$
Escríbase $T=\sum_{j,l}T_{jl}\otimes|j\rangle\langle l|$ con bloques $T_{jl}$ de tamaño $k_d\times k_c$. La igualdad
anterior dice $\pi_d(\bar a)T_{jl}=T_{jl}\overline{\pi_c(a)}$: cada $T_{jl}$ es un operador de entrelazamiento de
$\rho(a)=\overline{\pi_c(a)}$ en $\rho'(a)=\pi_d(\bar a)$. Estas son representaciones de $\mathcal A$ con su adjunción
en el sentido de la sección 4.1 ($\rho(a^*)=\rho(a)^*$, $\rho'(a^*)=\rho'(a)^*$), y son irreducibles, pues sus
imágenes son todo $M_{k_c}$ y todo $M_{k_d}$. Como $T\neq0$, algún $T_{jl}\neq0$, así que $k_c=k_d$ y
$\mathrm{Hom}(\rho,\rho')=\mathbb CJ_c$ para una unitaria $J_c$ (lema 4.2). Por tanto, $T_{jl}=\mu_{jl}J_c$, es
decir, $T=J_c\otimes M_c$ con $M_c=(\mu_{jl})$; $1=T^*T=1\otimes M_c^*M_c$ y $m_c=k_cm_c/k_c=m_d$, así que $M_c$ es
unitaria, y la igualdad se convierte en la fórmula para $\pi_d(\bar a)$. $T$ está determinada por $V$; si también
$T=J'\otimes M'$ con $J'$ unitaria, los bloques $T_{jl}=M'_{jl}J'$ muestran que $J'\in\mathrm{Hom}(\rho,\rho')$,
luego $J'=\lambda J_c$ con $|\lambda|=1$, y entonces $M'=\bar\lambda M_c$.

(b) Sean $\sigma(c)=c$ y $J=J_c$, $M=M_c$. Conjugando $\overline{V_c}=V_c(J\otimes M)$ se obtiene
$V_c=\overline{V_c}(\bar J\otimes\bar M)=V_c(J\bar J\otimes M\bar M)$, así que $J\bar J\otimes M\bar M=1$. Aplicando
dos veces la fórmula de (a), con $\overline{J^*}=J^T$,
$$
\pi_c(a)=\pi_c(\bar{\bar a})=J\,\overline{J\,\overline{\pi_c(a)}\,J^*}\,J^*=(J\bar J)\,\pi_c(a)\,(J\bar J)^*,
$$
de modo que la unitaria $J\bar J$ conmuta con $\pi_c(\mathcal A)=M_{k_c}$ y es un escalar $\lambda1$ (lema 2.5, paso
5, con $m=1$), $|\lambda|=1$. Entonces $\bar JJ=\overline{J\bar J}=\bar\lambda1$ y
$\lambda J=(J\bar J)J=J(\bar JJ)=\bar\lambda J$, así que $\lambda=\bar\lambda\in\{1,-1\}$; llamémoslo $\varepsilon_c$.
No cambia bajo $J\mapsto\lambda J$, y $\varepsilon_c1\otimes M\bar M=1$ da $M\bar M=\varepsilon_c1$. Si
$\varepsilon_c=-1$, entonces $(-1)^{k_c}=\det(J\bar J)=\det J\cdot\overline{\det J}=|\det J|^2>0$, así que $k_c$ es
par, y análogamente $m_c$.

(c) Las $P_c$ no dependen de $V$ (teorema 2.6), luego tampoco $\sigma$. Sea $V'$ otra descomposición, numerada de modo
que $V'_cV_c'^*=P_c$ (teorema 2.6), con homomorfismos $\pi'_c$ y datos $J'_c,M'_c$, y sea $\sigma(c)=c$. Entonces
$V'_c=P_cV'_c=V_cY$ con $Y=V_c^*V'_c$, una isometría cuadrada, luego unitaria, y
$(\pi_c(a)\otimes1)Y=Y(\pi'_c(a)\otimes1)$ para $a\in\mathcal A$. El argumento de (a), con las representaciones
irreducibles $\pi'_c,\pi_c$ en lugar de $\rho,\rho'$, da $Y=u\otimes w$ con unitarias $u,w$. Por tanto,
$\overline{V'_c}=\overline{V_c}(\bar u\otimes\bar w)=V_c(J_c\bar u\otimes M_c\bar w)=V'_c(u^*J_c\bar u\otimes w^*M_c\bar w)$,
así que $J'_c=\lambda u^*J_c\bar u$ con $|\lambda|=1$ por la unicidad de (a), y, como $\bar uu^T=\overline{uu^*}=1$,
$J'_c\bar J'_c=u^*J_c\bar uu^T\bar J_cu=u^*J_c\bar J_cu=\varepsilon_c1$.

(d) Por (a),
$\overline{V_c(1\otimes y)V_c^*}=V_{\sigma(c)}(J_c\otimes M_c)(1\otimes\bar y)(J_c\otimes M_c)^*V_{\sigma(c)}^*$, y
$(J_c\otimes M_c)(1\otimes\bar y)(J_c\otimes M_c)^*=J_cJ_c^*\otimes M_c\bar yM_c^*=1\otimes M_c\bar yM_c^*$. Súmese
sobre $c$ para $X=\sum_cV_c(1\otimes x_c)V_c^*$. $\square$

Un sector $c$ es de **tipo real** si $\sigma(c)=c$ y $\varepsilon_c=1$, de **tipo cuaterniónico** si $\sigma(c)=c$ y
$\varepsilon_c=-1$, y de **tipo complejo** si $\sigma(c)\neq c$; por el teorema 4.20(c), esta es una propiedad de
$\mathcal A$ y de $P_c$. Un sector con $\sigma(c)=c$ y $k_c$ impar o $m_c$ impar es de tipo real, por el teorema
4.20(b).

*En la biblioteca:* $\sigma$, $J_c$, $M_c$, $\varepsilon_c$ y los tipos de los sectores los calcula `Decomposition.real_structure()` — [tutorial, sección 4.7](../../tutorial/es/04_symmetries.md#s-4-7).

<a id="r-4-21"></a>

**Teorema 4.21 (el álgebra real).** Sea $\Gamma$ un conjunto de sectores que contiene todo sector $c$ con
$\sigma(c)=c$ y un sector de cada par $\{c,\sigma(c)\}$ con $\sigma(c)\neq c$. Para $\sigma(c)=c$, sea
$R_c=\{x\in M_{k_c}(\mathbb C):xJ_c=J_c\bar x\}$. Entonces $a\mapsto(\pi_c(a))_{c\in\Gamma}$ es un isomorfismo de
álgebras reales, que lleva la transpuesta a la adjunta,
$$
\mathcal A_{\mathbb R}\;\cong\bigoplus_{c\in\Gamma,\ \sigma(c)=c}R_c\;\oplus\bigoplus_{c\in\Gamma,\ \sigma(c)\neq c}M_{k_c}(\mathbb C),
$$
donde los últimos sumandos se consideran álgebras reales (**tipo complejo**). Además:
(a) (*tipo real*) Si $\varepsilon_c=1$, escríbase $J_c=W_JW_J^T$ y $M_c=W_MW_M^T$ con $J_c\bar W_J=W_J$,
$M_c\bar W_M=W_M$ (lema 4.19(c)). Entonces $x\mapsto W_J^*xW_J$ aplica $R_c$ isomorfamente sobre $M_{k_c}(\mathbb R)$.
Las columnas de $O_c=V_c(W_J\otimes W_M)$ son reales y ortonormales, generan la imagen de $P_c$ sobre $\mathbb C$, y
$aO_c=O_c(\alpha_c(a)\otimes1_{m_c})$ con la matriz real $\alpha_c(a)=W_J^*\pi_c(a)W_J$ para
$a\in\mathcal A_{\mathbb R}$.
(b) (*tipo cuaterniónico*) Si $\varepsilon_c=-1$, escríbase $J_c=W_J\Omega_{k_c}W_J^T$ con
$J_c\bar W_J=W_J\Omega_{k_c}$ (lema 4.19(c)). Entonces $x\mapsto\mathbb M^{-1}(W_J^*xW_J)$ es un isomorfismo
$R_c\cong M_{k_c/2}(\mathbb H)$ de álgebras reales que lleva la adjunta a la adjunta cuaterniónica
$(q^*)_{ij}=\bar q_{ji}$; aquí $\mathbb M:M_r(\mathbb H)\to M_{2r}(\mathbb C)$ sustituye cada entrada $q$ por el
bloque $2\times2$ $M(q)$ del lema E.15(2), y su imagen es $\{y\in M_{2r}(\mathbb C):y\Omega_{2r}=\Omega_{2r}\bar y\}$.
En particular, $\dim_{\mathbb R}\mathcal A_{\mathbb R}=\sum_ck_c^2=\dim_{\mathbb C}\mathcal A$.

*Demostración.* Para $a$ real, $\bar a=a$, y el teorema 4.20(a) se lee
$\pi_{\sigma(c)}(a)=J_c\overline{\pi_c(a)}J_c^*$. Para $\sigma(c)=c$ esto dice que $\pi_c(a)\in R_c$. Así,
$\Phi(a)=(\pi_c(a))_{c\in\Gamma}$ aplica $\mathcal A_{\mathbb R}$ en el lado derecho; es un homomorfismo de álgebras
reales, y $\Phi(a^T)=\Phi(a^*)=\Phi(a)^*$. Si $\Phi(a)=0$, la misma fórmula da $\pi_{\sigma(c)}(a)=0$ para
$c\in\Gamma$, así que todos los $\pi_e(a)=0$ y $a=0$. *Sobreyectividad:* sean $c\in\Gamma$, $x\in R_c$ si
$\sigma(c)=c$ y $x\in M_{k_c}$ en otro caso, y $a=V_c(x\otimes1)V_c^*\in\mathcal A$. Por el teorema 4.20(a) para el
sector $\sigma(e)$, $\pi_e(\bar a)=J_{\sigma(e)}\overline{\pi_{\sigma(e)}(a)}J_{\sigma(e)}^*$, que es $0$ salvo si
$\sigma(e)=c$. Si $\sigma(c)=c$, $\pi_c(\bar a)=J_c\bar xJ_c^*=x$ (definición de $R_c$), así que la matriz real
$\mathrm{Re}\,a=\frac12(a+\bar a)\in\mathcal A_{\mathbb R}$ tiene $\Phi(\mathrm{Re}\,a)=x$ en la posición $c$
y $0$ en las demás. Si $\sigma(c)\neq c$, $\pi_c(\bar a)=0$, y $2\mathrm{Re}\,a$ tiene $\pi_c=x$ y $\pi_e=0$ para
$e\notin\{c,\sigma(c)\}$; como $\sigma(c)\notin\Gamma$, $\Phi(2\mathrm{Re}\,a)$ es $x$ en la posición $c$ y $0$ en
las demás. Las sumas de estos elementos dan todo el lado derecho.

(a) Multiplíquese $xW_JW_J^T=W_JW_J^T\bar x$ por $W_J^*$ a la izquierda y por $\bar W_J$ a la derecha; como
$W_J^T\bar W_J=\overline{W_J^*W_J}=1$, se convierte en $W_J^*xW_J=W_J^T\bar x\bar W_J=\overline{W_J^*xW_J}$, y los
pasos se pueden invertir. Así, $R_c=W_JM_{k_c}(\mathbb R)W_J^*$, y $x\mapsto W_J^*xW_J$ es un isomorfismo de álgebras
reales compatible con la adjunta. Además,
$\overline{O_c}=\overline{V_c}(\bar W_J\otimes\bar W_M)=V_c(J_c\bar W_J\otimes M_c\bar W_M)=V_c(W_J\otimes W_M)=O_c$,
$O_c^*O_c=1$, y $aO_c=V_c(\pi_c(a)W_J\otimes W_M)=O_c(W_J^*\pi_c(a)W_J\otimes1)$, con $W_J^*\pi_c(a)W_J$ real para
$a\in\mathcal A_{\mathbb R}$.

(b) El mismo cálculo con $J_c=W_J\Omega W_J^T$ ($\Omega=\Omega_{k_c}$) muestra que $xJ_c=J_c\bar x$ si y solo si
$y=W_J^*xW_J$ satisface $y\Omega=\Omega\bar y$. Escríbase $y$ en bloques $2\times2$ $y_{ij}$; como $\Omega$ es
diagonal por bloques, $y\Omega=\Omega\bar y$ si y solo si $y_{ij}\Omega_2=\Omega_2\bar y_{ij}$ para todos los $i,j$, y
para $y_{ij}=\begin{pmatrix}\alpha&\beta\\\gamma&\delta\end{pmatrix}$ esto se lee
$\begin{pmatrix}\beta&-\alpha\\\delta&-\gamma\end{pmatrix}=\begin{pmatrix}-\bar\gamma&-\bar\delta\\\bar\alpha&\bar\beta\end{pmatrix}$,
es decir, $\delta=\bar\alpha$ y $\beta=-\bar\gamma$:
$y_{ij}=\begin{pmatrix}\alpha&-\bar\gamma\\\gamma&\bar\alpha\end{pmatrix}=M(\alpha,\bar\gamma)$ con
$M(a,b)=\begin{pmatrix}a&-b\\\bar b&\bar a\end{pmatrix}$. Así, la imagen de $\mathbb M$ es la indicada. Por el lema
E.15(2), $M$ es un homomorfismo real-lineal inyectivo $\mathbb H\to M_2(\mathbb C)$ con $M(\bar q)=M(q)^*$; como un
producto de matrices solo usa sumas y productos de entradas, $\mathbb M$ es un homomorfismo inyectivo de álgebras
reales con $\mathbb M(q^*)=\mathbb M(q)^*$ (como en la demostración del teorema E.24(1)). Por tanto,
$x\mapsto\mathbb M^{-1}(W_J^*xW_J)$ es el isomorfismo anunciado. *Dimensiones:* $\dim_{\mathbb R}R_c=k_c^2$ en ambos
casos ($\dim_{\mathbb R}M_r(\mathbb H)=4r^2$ con $2r=k_c$), y un par contribuye
$\dim_{\mathbb R}M_{k_c}(\mathbb C)=2k_c^2=k_c^2+k_{\sigma(c)}^2$. $\square$

Lo mismo vale para el conmutante, con los factores intercambiados: $(\mathcal A')_{\mathbb R}$ es la suma directa de
$M_{m_c}(\mathbb R)$, $M_{m_c/2}(\mathbb H)$ y $M_{m_c}(\mathbb C)$ sobre los mismos sectores, con el mismo tipo. En
efecto, $\mathcal A'$ es una $*$-álgebra con $\overline{\mathcal A'}=\mathcal A'$ (lema 4.18(b)), y la demostración
anterior solo usó tres hechos, que valen para $\mathcal A'$ con $X\mapsto x_c$ en lugar de $\pi_c$ y $M_c$ en lugar de
$J_c$: $X\mapsto x_c$ es un $*$-homomorfismo unital sobre $M_{m_c}$ ($V_c(1\otimes y)V_c^*\in\mathcal A'$ tiene el
bloque $y$ en el sector $c$ y $0$ en los demás), $(\bar X)_{\sigma(c)}=M_c\bar x_cM_c^*$ (teorema 4.20(d)), y
$M_c\bar M_c=\varepsilon_c1$ para $\sigma(c)=c$ (teorema 4.20(b)).

*En la biblioteca:* una base ortogonal real adaptada al teorema 4.21, `Decomposition.real_basis()` — [tutorial, sección 4.7](../../tutorial/es/04_symmetries.md#s-4-7).

Tres ejemplos. (i) Los generadores $(L_a)_{bc}=-\epsilon_{abc}$ de las
rotaciones de $\mathbb R^3$, que actúan sobre $\mathbb R^3\otimes\mathbb R^2$ como $L_a\otimes1_2$: hay un sector, $k=3$, $m=2$ ($iL_x,iL_y,iL_z$ es una representación de espín sobre $\mathbb C^3$ en la que
$iL_z$ tiene los valores propios simples $1,0,-1$, luego es el espín $1$ por el corolario 4.12); es de tipo real porque $k$ es impar, y
$\mathcal A_{\mathbb R}\cong M_3(\mathbb R)$. (ii) El generador de las rotaciones $\Omega_2$ de $\mathbb R^2$, que
actúa sobre $\mathbb R^2\otimes\mathbb R^2$ como $\Omega_2\otimes1_2$: como $\Omega_2^2=-1$ y $\Omega_2^*=-\Omega_2$,
el álgebra es $\mathrm{span}\{1,\Omega_2\otimes1\}$; los vectores propios $(1,\mp i)/\sqrt2$ de $\Omega_2$ para
$\pm i$, que son conjugados entre sí, dan dos sectores con $k=1$, $m=2$, intercambiados por la conjugación: un par de
tipo complejo, y $\mathcal A_{\mathbb R}=\{\alpha+\beta\,\Omega_2\otimes1\}\cong\mathbb C$. (iii) La multiplicación
por la izquierda $L_q$ por un cuaternión $q=(a,b)\in\mathbb H$ (sección E.6) sobre $\mathbb H=\mathbb R^4$, con
coordenadas $(\mathrm{Re}\,a,\mathrm{Im}\,a,\mathrm{Re}\,b,\mathrm{Im}\,b)$. $q\mapsto L_q$ es un
homomorfismo inyectivo de álgebras reales ($\mathbb H$ es asociativa, lema E.15(2)), y $L_{\bar q}=L_q^T$, porque
$\langle qx,y\rangle=\langle x,\bar qy\rangle$ (lema E.14(2)) para la norma $N(a,b)=|a|^2+|b|^2$, que es la norma
euclídea de $\mathbb R^4$. Así, el álgebra real generada por $L_I,L_J,L_K$ ($I=(i,0)$, $J=(0,1)$, $K=(0,i)$) y sus
transpuestas es $\{L_q:q\in\mathbb H\}\cong\mathbb H$; el álgebra compleja $\mathcal A$ es el subespacio complejo
generado por los productos de los generadores, $\{L_q\}+i\{L_q\}$, así que
$\mathcal A_{\mathbb R}=\{L_q\}\cong\mathbb H$ (lema 4.18(a)). Esto ya fija la descomposición:
$\sum_ck_c^2=\dim_{\mathbb C}\mathcal A=\dim_{\mathbb R}\mathcal A_{\mathbb R}=4$ (teorema 4.21), y no todos los $k_c$
valen $1$, pues entonces $\mathcal A$, y con ella $\mathcal A_{\mathbb R}\cong\mathbb H$, sería conmutativa; así que
hay un solo sector, con $k=2$ y $\sigma(c)=c$, y no es de tipo real, pues entonces el teorema 4.21(a) daría
$\mathcal A_{\mathbb R}\cong M_2(\mathbb R)$, que tiene divisores de cero, mientras que $\mathbb H$ no tiene ninguno
(lema E.15(3)): es de tipo cuaterniónico, con $m=4/2=2$ sobre $\mathbb R^4$. (Los cuaterniones de norma $1$ actúan
aquí como $SU(2)$: $M(q)$ es unitaria con determinante $|a|^2+|b|^2=1$, y toda
$U=\begin{pmatrix}\alpha&\beta\\\gamma&\delta\end{pmatrix}\in SU(2)$ es de esta forma, pues
$U^*=U^{-1}=\begin{pmatrix}\delta&-\beta\\-\gamma&\alpha\end{pmatrix}$ da $\delta=\bar\alpha$, $\beta=-\bar\gamma$,
así que $U=M(\alpha,\bar\gamma)$.) Sobre $\mathbb H^2=\mathbb R^8$, con $q$ actuando como $1_2\otimes L_q$, el mismo
argumento da un sector con $k=2$, $m=4$, de tipo cuaterniónico,
$\mathcal A_{\mathbb R}\cong M_1(\mathbb H)=\mathbb H$.

*Cálculo:* los tres ejemplos, con `Decomposition.real_structure()` y `Decomposition.real_basis()` — [tutorial, sección 4.7](../../tutorial/es/04_symmetries.md#s-4-7).

En la teoría de representaciones de un grupo por matrices reales, lo siguiente es la forma real del lema de Schur;
muestra que los tres tipos pueden leerse en la propia representación real.

<a id="r-4-22"></a>

**Teorema 4.22 (lema de Schur real).** Sea $U\subseteq\mathbb R^n$ un subespacio no nulo con $aU\subseteq U$ para todo
$a\in\mathcal A_{\mathbb R}$ que sea **irreducible**: $0$ y $U$ son sus únicos subespacios con esta propiedad. Sea
$D_U$ el conjunto de las aplicaciones real-lineales $T:U\to U$ con $Tau=aTu$ para todos los
$a\in\mathcal A_{\mathbb R}$, $u\in U$.
(a) $D_U$ es un álgebra real asociativa con unidad y sin divisores de cero, luego isomorfa a $\mathbb R$, $\mathbb C$
o $\mathbb H$ (teorema E.19).
(b) Para exactamente una órbita $\mathcal O=\{c,\sigma(c)\}$ de $\sigma$, $U$ está contenido en la imagen de
$P_{\mathcal O}=\sum_{e\in\mathcal O}P_e$, y $\dim_{\mathbb R}D_U=1,2,4$ y $\dim_{\mathbb R}U=k_c,2k_c,2k_c$ según que
$c$ sea de tipo real, complejo o cuaterniónico. Así, $D_U\cong\mathbb R,\mathbb C,\mathbb H$, respectivamente.
(c) Para toda órbita $\mathcal O$ existe un tal $U$ en la imagen de $P_{\mathcal O}$. El propio $\mathbb R^n$ es
irreducible si y solo si $\sigma$ tiene una sola órbita $\{c,\sigma(c)\}$ y $m_c=1$ (tipos real y complejo),
respectivamente $m_c=2$ (tipo cuaterniónico).

*Demostración.* (a) Para $T\in D_U$, $\ker T$ y $TU$ son subespacios invariantes, así que $T=0$ o $T$ es biyectiva, y
entonces $T^{-1}\in D_U$. $D_U$ es cerrada para la composición, que es asociativa, contiene $\mathrm{id}_U\neq0$, y
$ST=0$ con $S\neq0$ obliga a $T=S^{-1}ST=0$.

(b) *Subespacios invariantes y proyecciones reales.* Si $U'\subseteq\mathbb R^n$ es invariante, también lo es su
complemento ortogonal en $\mathbb R^n$ (para $y\perp U'$, $u\in U'$ y $a\in\mathcal A_{\mathbb R}$,
$u\cdot ay=a^Tu\cdot y=0$, pues $a^T\in\mathcal A_{\mathbb R}$), así que la proyección ortogonal $e$ sobre $U'$, un
idempotente real simétrico, conmuta con $\mathcal A_{\mathbb R}$; por el lema 4.18(b),
$e\in(\mathcal A')_{\mathbb R}$. Recíprocamente, $e\mathbb R^n$ es invariante para toda proyección
$e\in(\mathcal A')_{\mathbb R}$. Así, $U'\mapsto e$ es una biyección, que conserva la inclusión, de los subespacios
invariantes sobre el conjunto $\mathcal P$ de las proyecciones de $(\mathcal A')_{\mathbb R}$, y $U$ es irreducible si
y solo si su proyección $e$ es un elemento no nulo minimal de $\mathcal P$. Por el lema 4.18(a) aplicado a la imagen
$e\mathbb C^n$ (cuyos vectores reales son $e\mathbb R^n$), $\dim_{\mathbb R}U=\mathrm{rank}\,e$.

*$D_U\cong e(\mathcal A')_{\mathbb R}e$.* Para $T\in D_U$, sea $\hat T$ la matriz real $n\times n$ $x\mapsto T(ex)$.
Para $a\in\mathcal A_{\mathbb R}$, $\hat Tax=T(eax)=T(aex)=aT(ex)=a\hat Tx$, así que
$\hat T\in(\mathcal A')_{\mathbb R}$ (lema 4.18(b)) y $\hat T=e\hat Te$. Recíprocamente, para
$X\in(\mathcal A')_{\mathbb R}$, la restricción de $eXe$ a $U$ está en $D_U$. Estas aplicaciones son inversas entre sí
y multiplicativas. Además, $e\mathcal A'e$ es cerrado para la conjugación, y sus elementos reales son los
$e(\mathcal A')_{\mathbb R}e$ (si $eXe$ es real, $eXe=e\mathrm{Re}(X)e$ con
$\mathrm{Re}\,X\in(\mathcal A')_{\mathbb R}$), así que $\dim_{\mathbb R}D_U=\dim_{\mathbb C}e\mathcal A'e$ por el
lema 4.18(a).

*Proyecciones en la descomposición.* Una proyección $e\in\mathcal A'$ es $e=\sum_cV_c(1\otimes f_c)V_c^*$ con
proyecciones $f_c\in M_{m_c}$ (teorema 2.6). Entonces $\mathrm{rank}\,e=\sum_ck_c\mathrm{rank}\,f_c$,
$e\mathcal A'e=\sum_cV_c(1\otimes f_cM_{m_c}f_c)V_c^*$ tiene dimensión compleja $\sum_c(\mathrm{rank}\,f_c)^2$ y,
por el teorema 4.20(d), $\bar e=e$ si y solo si $f_{\sigma(c)}=M_c\bar f_cM_c^*$ para todo $c$.

*Elementos minimales.* Sea $e$ minimal en $\mathcal P$. Para cada órbita $\mathcal O$ de $\sigma$, $P_{\mathcal O}$ es
una proyección central real ($\bar P_{\mathcal O}=P_{\mathcal O}$ por el teorema 4.20(a)) en
$Z(\mathcal A)\subseteq\mathcal A'$, que conmuta con $e$, así que $eP_{\mathcal O}\in\mathcal P$ y
$eP_{\mathcal O}\le e$; como $e=\sum_{\mathcal O}eP_{\mathcal O}$, exactamente una órbita $\mathcal O$ tiene
$eP_{\mathcal O}=e$, es decir, $U=e\mathbb R^n$ está en la imagen de $P_{\mathcal O}$, y $f_c=0$ para
$c\notin\mathcal O$. Si $\mathcal O=\{c,d\}$ con $d\neq c$, entonces $f_d=M_c\bar f_cM_c^*$ tiene el rango de $f_c$;
para una proyección de rango uno $g\le f_c$, $g'=Y+\bar Y$ con $Y=V_c(1\otimes g)V_c^*$ es una suma de dos
proyecciones ortogonales de $\mathcal A'$ (en los sectores $c$ y $d$), real, y $g'\le e$ (pues $Y\le e$ y
$\bar Y\le\bar e=e$); así que $g'=e$, $\mathrm{rank}\,f_c=\mathrm{rank}\,f_d=1$,
$\dim_{\mathbb C}e\mathcal A'e=2$ y $\mathrm{rank}\,e=2k_c$. Si $\mathcal O=\{c\}$, sean $f=f_c$, $M=M_c$ y
$\theta v=M\bar v$, antiunitaria sobre $\mathbb C^{m_c}$ con $\theta^2=M\bar M=\varepsilon_c$ (teorema 4.20(b)). De
$fM=M\bar f$ se sigue $\theta(fv)=M\bar f\bar v=fM\bar v=f\theta v$, así que la imagen $F$ de $f$ es
$\theta$-invariante. Si $F_1\subseteq F$ es $\theta$-invariante con proyección $g$, entonces $M\bar gM^*$ es la
proyección sobre $M\bar F_1=\theta F_1=F_1$ ($\theta$ es inyectiva), así que $g=M\bar gM^*$ y
$V_c(1\otimes g)V_c^*\in\mathcal P$ está por debajo de $e$. Por minimalidad, $F$ no tiene más subespacios
$\theta$-invariantes que $0$ y $F$. Si $\varepsilon_c=1$, un vector unitario fijo por $\theta$ en $F$ (lema 4.19(a))
genera uno, así que $\mathrm{rank}\,f=1$, $\dim_{\mathbb C}e\mathcal A'e=1$ y $\mathrm{rank}\,e=k_c$. Si
$\varepsilon_c=-1$, el plano $\mathrm{span}\{u,\theta u\}$ de un vector unitario $u\in F$ es uno (lema 4.19(b)),
así que $\mathrm{rank}\,f=2$, $\dim_{\mathbb C}e\mathcal A'e=4$ y $\mathrm{rank}\,e=2k_c$. Como $\mathbb R$,
$\mathbb C$, $\mathbb H$ tienen dimensiones $1,2,4$, (a) da $D_U\cong\mathbb R,\mathbb C,\mathbb H$, respectivamente.

(c) $P_{\mathcal O}\in\mathcal P$ es no nula, y un elemento de $\mathcal P\setminus\{0\}$ por debajo de ella de rango
mínimo es minimal. $\mathbb R^n$ es irreducible si y solo si $e=1$ es minimal; sus bloques son $f_c=1_{m_c}$, así que,
por la demostración de (b), esto exige una sola órbita y $m_c=1,1,2$ en los tres casos. Recíprocamente, bajo estas
condiciones, un $e'\in\mathcal P$ no nulo tiene $f'_c\neq0$ para los sectores de la órbita
($f'_{\sigma(c)}=M_c\bar f'_cM_c^*$), así que $f'_c=1$ si $m_c=1$; en el caso cuaterniónico, la imagen de $f'_c$ es
$\theta$-invariante, luego de dimensión par (lema 4.19(b)), así que $f'_c=1_2$. Por tanto, $e'=1$. $\square$

Para un grupo $G$ que actúa por matrices ortogonales reales, $\mathcal A=\mathrm{span}\,\pi(G)$ es cerrada para la
conjugación (lema 4.18(c)), $\mathcal A_{\mathbb R}$ es el subespacio real generado por $\pi(G)$ (lema 4.18(a):
$\mathcal A_{\mathbb R}$ consta de las partes reales de los elementos de $\mathcal A$), y $D_U$ es el conmutante de
una representación real irreducible de $G$: es $\mathbb R$, $\mathbb C$ o $\mathbb H$. Por la demostración de (b), la
complexificación $U+iU=e\mathbb C^n=V_c(\mathbb C^{k_c}\otimes F)$ (más su conjugada en el caso de un par) lleva
$\pi_c$, $\pi_c\oplus\pi_{\sigma(c)}$ o $\pi_c\oplus\pi_c$ en los tres casos: es irreducible, la suma de dos
representaciones irreducibles no equivalentes, o la suma de dos equivalentes.

*Cálculo:* el teorema 4.22 sobre los módulos reales irreducibles $\mathbb R^3$, $\mathbb R^2$, $\mathbb R^4$ de los tres ejemplos — [tutorial, sección 4.7](../../tutorial/es/04_symmetries.md#s-4-7).

<a id="r-4-23"></a>

**Proposición 4.23 (elementos reales del conmutante).** Sea $X\in\mathcal A'$ real, con bloques $x_c$. Entonces
$x_{\sigma(c)}=M_c\bar x_cM_c^*$ para todo $c$. Si además $X=X^T$:
(a) para $c$ de tipo real, $y_c=W_M^*x_cW_M$ ($W_M$ como en el teorema 4.21(a)) es una matriz real simétrica, y
$XO_c=O_c(1_{k_c}\otimes y_c)$ con las columnas reales ortonormales $O_c$ del teorema 4.21(a);
(b) para $c$ de tipo cuaterniónico, todo valor propio de $x_c$ tiene multiplicidad par, de modo que todo valor propio
de $X$ sobre la imagen de $P_c$ tiene una multiplicidad divisible por $2k_c$ (**degeneración de Kramers**);
(c) para $c$ de tipo complejo, $x_c$ y $x_{\sigma(c)}$ tienen los mismos valores propios con las mismas
multiplicidades.
Sin $X=X^T$, $y_c$ en (a) sigue siendo real.

*Demostración.* $X=\bar X$, así que el teorema 4.20(d) da $x_{\sigma(c)}=M_c\bar x_cM_c^*$. Si $X=X^T$, entonces
$X^*=X$ y toda $x_c$ es hermítica. (a) Con $M_c=W_MW_M^T$, $M_c^*=\bar W_MW_M^*$ y $W_M^T\bar W_M=1$:
$\bar y_c=W_M^T\bar x_c\bar W_M=W_M^TM_c^*x_cM_c\bar W_M=W_M^*x_cW_M=y_c$, e $y_c$ es hermítica, luego real simétrica.
Además, $XO_c=XV_c(W_J\otimes W_M)=V_c(W_J\otimes x_cW_M)=O_c(1\otimes y_c)$. (b) Sea $\theta v=M_c\bar v$,
antiunitaria con $\theta^2=-1$ (teorema 4.20(b)). De $x_cM_c=M_c\bar x_c$ se sigue
$x_c\theta v=M_c\bar x_c\bar v=\theta(x_cv)$, así que, para un valor propio $\lambda\in\mathbb R$ de $x_c$, el
subespacio propio $E_\lambda$ satisface $\theta E_\lambda\subseteq E_\lambda$
($x_c\theta v=\theta(\lambda v)=\lambda\theta v$), y $\dim E_\lambda$ es par por el lema 4.19(b). Sobre la imagen de
$P_c$, $X$ actúa como $V_c(1_{k_c}\otimes x_c)V_c^*$, cuyos valores propios son los de $x_c$, cada uno repetido $k_c$
veces. (c) $\bar x_c$ tiene los valores propios $\bar\lambda=\lambda$ de $x_c$, con los vectores propios conjugados, y
$x_{\sigma(c)}$ es unitariamente equivalente a $\bar x_c$. $\square$

En el caso cuaterniónico, $y_c=W_M^*x_cW_M$, con $M_c=W_M\Omega_{m_c}W_M^T$ como en el lema 4.19(c), satisface
$y_c\Omega=\Omega\bar y_c$ (el cálculo del teorema 4.21(b) con $M_c$ en lugar de $J_c$), así que $y_c=\mathbb M(h)$
para una matriz cuaterniónica $h$; si $X=X^T$, entonces $h^*=h$, pues $\mathbb M(h^*)=y_c^*=y_c$ y $\mathbb M$ es
inyectiva.

*Cálculo:* la proposición 4.23 para un elemento simétrico real aleatorio del conmutante en cada uno de los tres ejemplos — [tutorial, sección 4.7](../../tutorial/es/04_symmetries.md#s-4-7).

> **Significado físico** (interpretación física, no un teorema de este curso). Si un hamiltoniano conmuta con una
> simetría y es real en una base en la que la simetría es real —el efecto de una simetría de inversión temporal, que
> allí actúa como conjugación compleja—, la proposición 4.23 dice qué aspecto tienen sus bloques: reales simétricos en
> los sectores de tipo real, hermíticos complejos con un compañero del mismo espectro en los pares de tipo complejo, y
> hermíticos cuaterniónicos con niveles doblemente degenerados en los sectores de tipo cuaterniónico, como ocurre con
> la inversión temporal de un número impar de espines $\frac12$, cuyo cuadrado es $-1$ (teorema de Kramers). La *vía
> triple* de Dyson (*threefold way*; F. J. Dyson, *J. Math. Phys.* 3 (1962) 1199) modela los bloques de un
> hamiltoniano genérico de cada clase mediante matrices aleatorias de la clase correspondiente —el GOE y el GUE de la
> sección 5.3, y el conjunto simpléctico gaussiano (GSE), que este curso no define— y predice en consecuencia la
> estadística de niveles de cada bloque.

<a id="s-4-8"></a>

## 4.8 Con qué comparar

*En la biblioteca:* la comparación con otro software — [tutorial, sección 4.8](../../tutorial/es/04_symmetries.md#s-4-8).

<a id="x-4"></a>

## Ejercicios

1. Demuestre que, para la representación de $S_L$ sobre $(\mathbb C^d)^{\otimes L}$,
   $\chi_\pi(\sigma)=d^{\#\text{cycles}(\sigma)}$ y $\dim\pi(S_L)'=\binom{d^2+L-1}L$. Compruebe $20$ (sección 2.4) y
   $35$ ($=\dim\{S_a\}''$ para cuatro espines).
2. (*Clebsch–Gordan.*) Demuestre que $j_1\otimes j_2\cong\bigoplus_jj$, $j=|j_1-j_2|,|j_1-j_2|+1,\dots,j_1+j_2$, para
   la representación de espín $S^{(j_1)}_a\otimes1+1\otimes S^{(j_2)}_a$ —la regla de fusión no abeliana de $SU(2)$
   que menciona el capítulo 7.
3. Sea $\theta/2\pi$ irracional. Demuestre que $\{U_\theta\}'=U(1)'$ para una sola $U_\theta=e^{i\theta Q}$, pero que
   $\{e^{i\pi\hat N}\}'\cong M_{2^{n-1}}\oplus M_{2^{n-1}}$ para $n$ cúbits.
4. (*Desdoblamiento Zeeman.*) Sea $H$ un hamiltoniano de Heisenberg sobre $L\ge1$ espines y $H_h=H+hS_z$. Demuestre
   que $H_h$ conmuta con $S_z$ y con $\mathbf S^2$ pero, para $h\neq0$, no con $S_x$, y que sus valores propios son
   $E+hM$, $M=-S,\dots,S$, para los valores propios $E$ de $H_S$.
5. (*El grupo de los cuaterniones.*) Sea $Q_8=\{\pm1,\pm I,\pm J,\pm K\}\subseteq\mathbb H$ con $I=(i,0)$, $J=(0,1)$,
   $K=(0,i)$ (sección E.6), y sea $\pi$ su representación regular sobre $\mathbb R^8$, $\pi(g)e_h=e_{gh}$. Demuestre
   que $\mathcal A=\mathrm{span}\,\pi(Q_8)$ tiene cuatro sectores de tipo real con $k=m=1$ y uno de tipo
   cuaterniónico con $k=m=2$, de modo que $\mathcal A_{\mathbb R}\cong\mathbb R^4\oplus\mathbb H$.

<details>
<summary>Soluciones</summary>

1. Por el lema B.1(c), $\pi(\sigma)$ es una matriz de permutación que deja fijo $|b_1\dots b_L\rangle$ si y solo si $b$
   es constante en los ciclos de $\sigma$, y hay $d^{\#\text{cycles}(\sigma)}$ tales $b$. Por el corolario 4.5(c),
   $\dim\pi(S_L)'=\frac1{L!}\sum_\sigma(d^2)^{\#\text{cycles}(\sigma)}=\frac{x(x+1)\cdots(x+L-1)}{L!}$ con $x=d^2$
   (la identidad $\sum_\sigma x^{\#\text{cycles}(\sigma)}=x(x+1)\cdots(x+L-1)$ se demuestra en la solución del capítulo 5,
   ejercicio 2), que es $\binom{d^2+L-1}L$ (la proposición B.3 demuestra la misma dimensión contando multiconjuntos):
   $\binom63=20$ y $\binom74=35$. Que $35$ sea también $\dim\{S_a\}''$ es el teorema 4.14(b) para $L=4$:
   $\{S_a\}''=\pi(S_4)'$.

2. $S_a=S^{(j_1)}_a\otimes1+1\otimes S^{(j_2)}_a$ es una representación de espín (sección 4.5), y la unitaria
   $x\otimes y\mapsto y\otimes x$ la entrelaza con la de $(j_2,j_1)$; así que sea $j_1\ge j_2$. $S_z$ tiene el valor propio
   $m_1+m_2$ sobre $|m_1\rangle\otimes|m_2\rangle$, luego sus valores propios están en $j_1+j_2+\mathbb Z$, y para tales $M$,
   $d_M$ es el número de $m_2\in\{-j_2,\dots,j_2\}$ con $|M-m_2|\le j_1$ (entonces $m_1=M-m_2\in j_1+\mathbb Z$
   automáticamente); $d_M=0$ para $M\notin j_1+j_2+\mathbb Z$. Para $M\ge0$ la condición se lee $M-j_1\le m_2\le M+j_1$, y
   $M+j_1\ge j_2$, así que $m_2$ va de $\max(-j_2,M-j_1)$ a $j_2$: $d_M=2j_2+1$ para $0\le M\le j_1-j_2$,
   $d_M=j_1+j_2-M+1$ para $j_1-j_2\le M\le j_1+j_2$, y $d_M=0$ para $M>j_1+j_2$. Por el corolario 4.12, la multiplicidad
   del espín $j\ge0$ es $m_j=d_j-d_{j+1}$. Vale $0$ si $j\notin j_1+j_2+\mathbb Z$ o $j>j_1+j_2$; $m_{j_1+j_2}=1-0=1$; para
   $j_1-j_2\le j<j_1+j_2$ (con $j\in j_1+j_2+\mathbb Z$, luego $j+1\le j_1+j_2$),
   $m_j=(j_1+j_2-j+1)-(j_1+j_2-j)=1$; y para $0\le j<j_1-j_2$ (luego $j+1\le j_1-j_2$), $m_j=(2j_2+1)-(2j_2+1)=0$.
   Por tanto, $j_1\otimes j_2\cong\bigoplus_{j=j_1-j_2}^{j_1+j_2}j$, cada espín una vez. Para $j_2=\frac12$ esto es el
   teorema 4.16(d), que además da los coeficientes.

3. $U_\theta=\sum_qe^{iq\theta}P_q$, y $e^{iq\theta}=e^{iq'\theta}$ si y solo si $(q-q')\theta\in2\pi\mathbb Z$, lo que para
   $\theta/2\pi$ irracional obliga a $q=q'$. Así, los subespacios propios de $U_\theta$ son las imágenes de los $P_q$, y
   un operador conmuta con la matriz normal $U_\theta$ si y solo si conserva sus subespacios propios (teorema 1.5(d) e
   interpolación de Lagrange como en la sección 2.2), es decir, si y solo si $X=\sum_qP_qXP_q$; esto es $U(1)'$ por la
   proposición 4.10. El argumento solo usa que los $e^{iq\theta}$ para las cargas $q$ que aparecen son distintos.
   Para $\theta=\pi$, el valor propio $(-1)^w$ solo ve la paridad del peso de Hamming: dos
   subespacios propios de dimensión $2^{n-1}$ (para $n\ge1$, las configuraciones de peso par y de peso impar; invertir
   el primer bit es una biyección entre ellas), y el conmutante es $M_{2^{n-1}}\oplus M_{2^{n-1}}$.

4. $S_z$ y $\mathbf S^2$ conmutan con $H$ (proposición 4.15(a)) y con $S_z$, luego con $H_h$; y
   $[H_h,S_x]=h[S_z,S_x]=ihS_y\neq0$. En la descomposición de la proposición 4.13,
   $V^*S_zV=\bigoplus_SS^{(S)}_z\otimes1_{m_S}$ y $V^*HV=\bigoplus_S1_{2S+1}\otimes H_S$ (proposición 4.15), luego sobre
   el sumando de espín $S$, $V^*H_hV=hS^{(S)}_z\otimes1+1\otimes H_S$, una suma de operadores que conmutan,
   diagonalizada por los productos de vectores propios, con valores propios $hM+E$.

5. Por la demostración del teorema E.19, $I^2=J^2=K^2=-1$ e $IJ=-JI=K$, $JK=-KJ=I$, $KI=-IK=J$, así que $Q_8$ es un
   grupo, y $h\mapsto gh$ es una biyección de $Q_8$: los $\pi(g)$ son matrices de permutación, una representación
   unitaria por matrices reales, y $\mathcal A=\pi(Q_8)''$ es cerrada para la conjugación (lema 4.18(c)).
   $\chi_\pi(g)$ es el número de $h$ con $gh=h$, es decir, $8$ para $g=1$ y $0$ en otro caso, así que, por el
   corolario 4.5(b), una $\rho$ unitaria irreducible tiene multiplicidad
   $m_\rho=\frac18\overline{\chi_\rho(1)}\cdot8=\dim\rho$. *Unidimensionales:* salvo signos, $I,J,K$ se multiplican
   como los elementos distintos de la identidad de $\mathbb Z_2\times\mathbb Z_2$, así que, para $s,t\in\{\pm1\}$,
   $\chi(\pm1)=1$, $\chi(\pm I)=s$, $\chi(\pm J)=t$, $\chi(\pm K)=st$ es un homomorfismo; son cuatro caracteres
   distintos, luego cuatro representaciones irreducibles no equivalentes, cada una con $m=1$. *Bidimensional:*
   $\rho(g)=M(g)$ (lema E.15(2)) es una representación unitaria ($M(g)^*M(g)=M(\bar gg)=1$) con $\chi_\rho(\pm1)=\pm2$
   y $\chi_\rho=0$ en los demás elementos ($\mathrm{Tr}\,M(a,b)=2\mathrm{Re}\,a$), así que
   $\langle\chi_\rho,\chi_\rho\rangle_{Q_8}=\frac{4+4}8=1$ y $\rho$ es irreducible (corolario 4.5(c)), con $m=2$. Para
   una descomposición $V$ como en el teorema 2.6, las $\pi_c\circ\pi$ son irreducibles y no equivalentes dos a dos
   (proposición 2.7), así que $V$ es una descomposición como en el teorema 4.3, y cada una de las cinco
   representaciones anteriores, al tener $m_\rho\ge1$, es equivalente a exactamente una $\pi_c\circ\pi$, y
   representaciones distintas a $c$ distintos (corolario 4.5(b)). Como $4\cdot1\cdot1+2\cdot2=8=\sum_ck_cm_c$, no hay
   más sectores. Por el teorema 4.20(c) podemos usar una descomposición con $\pi_c(\pi(g))=\rho_c(g)$ para estas
   representaciones (el lema 4.2(c) da unitarias $u$ con $\rho_c=u\,\pi_c(\pi(\cdot))\,u^*$; sustitúyase $V_c$ por
   $V_c(u^*\otimes1)$). Como $\pi(g)$ es real, el teorema 4.20(a) da
   $\pi_{\sigma(c)}(\pi(g))=J_c\overline{\rho_c(g)}J_c^*$. Para un carácter $\chi$ esto es $\chi(g)$, así que
   $\pi_{\sigma(c)}\circ\pi\cong\pi_c\circ\pi$, $\sigma(c)=c$ (proposición 2.7: sectores distintos llevan
   representaciones no equivalentes), y $J_c\bar J_c=|J_c|^2=1$: tipo real. Para $\rho=M$, $\sigma(c)=c$ (el único
   sector con $k=2$), y $\overline{M(g)}=\Omega_2^*M(g)\Omega_2$ (por la demostración del teorema 4.21(b),
   $M(q)\Omega_2=\Omega_2\overline{M(q)}$), así que $M(g)=(J_c\Omega_2^*)M(g)(J_c\Omega_2^*)^*$. Las matrices
   $M(1)=1$, $M(I)=\mathrm{diag}(i,-i)$, $M(J)=\Omega_2$, $M(K)=\begin{pmatrix}0&-i\\-i&0\end{pmatrix}$ son
   linealmente independientes, así que generan $M_2(\mathbb C)$; por tanto, $J_c\Omega_2^*=\lambda1$,
   $J_c=\lambda\Omega_2$ y $J_c\bar J_c=|\lambda|^2\Omega_2^2=-1$: tipo cuaterniónico. Por el teorema 4.21,
   $\mathcal A_{\mathbb R}\cong\mathbb R^4\oplus M_1(\mathbb H)$, de dimensión $8=|Q_8|$.

</details>

*Cálculo:* los ejercicios 1–5, numéricamente — [tutorial, Ejercicios](../../tutorial/es/04_symmetries.md#x-4).

**Siguiente:** el [capítulo 5](05_random_matrices.md) extrae matrices aleatorias en un álgebra: unitarias de Haar de
$\prod_cU(k_c)$, el cálculo de Weingarten —cuyo caso de orden dos se apoya en el conmutante de
$\{V\otimes V\}$, el caso más pequeño de la dualidad de Schur–Weyl— y la estadística de niveles sector por sector.
