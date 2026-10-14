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
carga, las representaciones irreducibles de espín y cómo se descomponen $L$ espines $\frac12$, y por qué
los hamiltonianos invariantes por $SU(2)$ tienen bloques mucho más pequeños que los invariantes por
$U(1)$.

En todo el capítulo, $S'$ es el conmutante de $S\subseteq M_n(\mathbb C)$ (capítulo 2). Si $S^*=S$, $S''$ es
el **álgebra de von Neumann generada por** $S$: una $*$-álgebra que contiene $S$, y toda $*$-álgebra
$\mathcal B\supseteq S$ satisface $\mathcal B=\mathcal B''\supseteq S''$ (teorema 2.2). $\mathcal A'_h=\{K\in\mathcal A':K=K^*\}$ es un
espacio vectorial real. $\mathbb C^k\otimes\mathbb C^m$ lleva el orden de Kronecker del capítulo 7 ($e_a\otimes f_j$
es el vector de la base $am+j$, contando desde $0$), de modo que $1_k\otimes x$ es diagonal por bloques con
$k$ copias de $x$. $|G|$ es el orden de un grupo finito, y $\hbar=1$.

## 4.1 Representaciones y lema de Schur

**Definición.** Sea $X$ un conjunto con una involución $x\mapsto x^*$ ($x^{**}=x$). Una
**representación** de $X$ sobre $\mathbb C^k$ es una aplicación $\rho:X\to M_k(\mathbb C)$ con
$\rho(x^*)=\rho(x)^*$. Un subespacio $K$ es **invariante** si $\rho(x)K\subseteq K$ para todo $x$; $\rho$ es
**irreducible** si $k\ge1$ y $0,\mathbb C^k$ son los únicos subespacios invariantes. Un
**operador de entrelazamiento** de $\rho$ en una representación $\rho'$ sobre $\mathbb C^{k'}$ es una matriz
$T$ de tamaño $k'\times k$ con $T\rho(x)=\rho'(x)T$ para todo $x$; forman un espacio vectorial
$\operatorname{Hom}(\rho,\rho')$, y $\operatorname{Hom}(\rho,\rho)=\rho(X)'$. $\rho$ y $\rho'$ son **equivalentes**,
$\rho\cong\rho'$, si existe un operador de entrelazamiento invertible.

Se presentan tres casos: un grupo $G$ con $g^*=g^{-1}$, donde exigimos que $\rho$ sea un homomorfismo
$G\to U(k)$ (una **representación unitaria**; entonces $\rho(g^{-1})=\rho(g)^{-1}=\rho(g)^*$); una $*$-álgebra con
su adjunción, donde $\rho$ es un $*$-homomorfismo unital; y etiquetas de generadores con $x^*=x$, por ejemplo
$\{x,y,z\}$, representadas por tuplas de matrices hermíticas, como las matrices de espín (sección 4.5). Si $K$
es invariante y las columnas de $Q\in M_{k\times d}$ forman una base ortonormal de $K$, entonces
$\rho(x)Q=QQ^*\rho(x)Q$, de modo que $\rho_K(x)=Q^*\rho(x)Q$ es de nuevo una representación (un homomorfismo si
$\rho$ lo es), la **restricción** de $\rho$ a $K$.

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

**Lema 4.2 (lema de Schur).** Sean $\rho,\rho'$ representaciones irreducibles de $X$ sobre
$\mathbb C^k,\mathbb C^{k'}$ y $T\in\operatorname{Hom}(\rho,\rho')$.
(a) $T=0$ o $T$ es invertible. (b) $\operatorname{Hom}(\rho,\rho)=\mathbb C1$.
(c) $T^*\in\operatorname{Hom}(\rho',\rho)$, y si $T\neq0$, entonces $T=\lambda U$ con $\lambda>0$ y $U$ unitaria.
Por tanto, $\dim\operatorname{Hom}(\rho,\rho')$ es $1$ si $\rho\cong\rho'$ y $0$ en otro caso, y las representaciones
irreducibles equivalentes son unitariamente equivalentes.

*Demostración.* (a) y (b) se demostraron para grupos tras la proposición 2.7, y la demostración solo usa
que $\ker T$ es invariante por $\rho$ ($T\rho(x)v=\rho'(x)Tv=0$ si $Tv=0$) y la imagen de $T$ por
$\rho'$: así, $\ker T\in\{0,\mathbb C^k\}$ y $\operatorname{ran}T\in\{0,\mathbb C^{k'}\}$; y para (b), $T-\lambda1$,
con $\lambda$ un valor propio de $T$, es un operador de entrelazamiento no invertible, luego $0$. (c) El adjunto
de $T\rho(x^*)=\rho'(x^*)T$ es $\rho(x)T^*=T^*\rho'(x)$. Entonces $T^*T\in\operatorname{Hom}(\rho,\rho)$, luego
$T^*T=\lambda^21$ con $\lambda>0$ si $T\neq0$ ($T^*T\ge0$ es no nulo), y $U=T/\lambda$ es cuadrada por
(a) con $U^*U=1$. Finalmente, si $T_1,T_2\in\operatorname{Hom}(\rho,\rho')$ son no nulos, entonces
$T_1^{-1}T_2\in\operatorname{Hom}(\rho,\rho)=\mathbb C1$. $\square$

Si todos los $\rho(x)$ conmutan entre sí (un grupo abeliano), una $\rho$ irreducible es unidimensional:
todo $\rho(x)$ está en $\operatorname{Hom}(\rho,\rho)=\mathbb C1$, de modo que todo subespacio es invariante.

## 4.2 La descomposición isotípica y el conmutante

**Teorema 4.3 (descomposición isotípica).** Sea $\rho$ una representación de $X$ sobre $\mathbb C^n$.
Existen representaciones irreducibles, no equivalentes dos a dos, $\rho_1,\dots,\rho_C$ de $X$ sobre
$\mathbb C^{k_c}$, enteros $m_c\ge1$ con $\sum_ck_cm_c=n$, y una unitaria
$V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n$ con
$V^*\rho(x)V=\bigoplus_c\rho_c(x)\otimes1_{m_c}$ para todo $x$. Para toda descomposición de este tipo
$$
V^*\rho(X)'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),\qquad
V^*\rho(X)''V=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},
$$
y $m_c=\dim\operatorname{Hom}(\rho_c,\rho)$; una representación irreducible $\sigma$ de $X$ cumple
$\operatorname{Hom}(\sigma,\rho)\neq0$ si y solo si $\sigma\cong\rho_c$ para exactamente un $c$.

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
$\operatorname{Hom}(\rho_d,\rho_c)$, que es $0$ para $c\neq d$ y $\mathbb C1$ para $c=d$ (lema 4.2). Así,
$Y\in\rho(X)'$ si y solo si $Y_{cj,dl}=\delta_{cd}(y_c)_{jl}1_{k_c}$, es decir, $V^*YV=\bigoplus_c1_{k_c}\otimes y_c$.

*Doble conmutante.* Si $Z$ conmuta con $\rho(X)'$, conmuta con las proyecciones sobre los sumandos
(tómese $y_c=1$ e $y_d=0$ para $d\neq c$), luego $V^*ZV=\bigoplus_cZ_c$, y $Z_c$ conmuta con
$1_{k_c}\otimes M_{m_c}$; por el lema 2.5, paso 5 (con los dos factores tensoriales intercambiados),
$Z_c\in M_{k_c}\otimes1_{m_c}$. Recíprocamente, estos operadores conmutan con $\rho(X)'$.

*Multiplicidades.* Para $T\in\operatorname{Hom}(\sigma,\rho)$ con $\sigma$ irreducible, el mismo argumento muestra
que $Q_{c,j}^*T\in\operatorname{Hom}(\sigma,\rho_c)$, que es $\mathbb Ct_{c}$ con $t_c$ invertible si $\sigma\cong\rho_c$
y $0$ en otro caso; y $T=\sum_{c,j}Q_{c,j}Q_{c,j}^*T$. Como las $\rho_c$ son no equivalentes dos a dos, $\sigma$
es equivalente a lo sumo a una de ellas, digamos $\rho_c$ (si no, $T=0$), y entonces $T=\sum_ja_jQ_{c,j}t_c$;
recíprocamente, cada $Q_{c,j}t_c$ pertenece a $\operatorname{Hom}(\sigma,\rho)$, y estas $m_c$ aplicaciones tienen
imágenes ortogonales, luego $\dim\operatorname{Hom}(\sigma,\rho)=m_c$. $\square$

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

**Teorema 4.4 (ortogonalidad de Schur).** Sean $\rho,\rho'$ representaciones unitarias irreducibles de un
grupo finito $G$ sobre $\mathbb C^k,\mathbb C^{k'}$. Entonces
$$
\frac1{|G|}\sum_{g\in G}\rho(g)_{ij}\,\overline{\rho'(g)_{i'j'}}=
\begin{cases}\delta_{ii'}\delta_{jj'}/k&\text{if }\rho=\rho',\\0&\text{if }\rho\not\cong\rho'.\end{cases}
$$

*Demostración.* Para una matriz $Y$ de tamaño $k\times k'$, sea $T=\frac1{|G|}\sum_g\rho(g)Y\rho'(g)^*$. Para
$h\in G$, la sustitución $g\to hg$ y la igualdad $\rho'(hg)^*=\rho'(g)^*\rho'(h)^{-1}$ dan $\rho(h)T=T\rho'(h)$, luego
$T\in\operatorname{Hom}(\rho',\rho)$. Si $\rho\not\cong\rho'$, $T=0$; si $\rho=\rho'$, $T=\lambda1$ con
$k\lambda=\operatorname{Tr}T=\operatorname{Tr}Y$ (lema 4.2). Para $Y=|j\rangle\langle j'|$, la entrada $T_{ii'}$
es el miembro izquierdo, y $\operatorname{Tr}Y=\delta_{jj'}$. $\square$

El **carácter** de una representación es $\chi_\rho(g)=\operatorname{Tr}\rho(g)$; las representaciones
equivalentes tienen caracteres iguales, pues $\operatorname{Tr}(T\rho(g)T^{-1})=\operatorname{Tr}\rho(g)$. Escribimos
$\langle\varphi,\psi\rangle_G=\frac1{|G|}\sum_g\overline{\varphi(g)}\psi(g)$.

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
si $\rho\not\cong\rho'$. (b) La traza es invariante unitaria y $\operatorname{Tr}(a\otimes1_m)=m\operatorname{Tr}a$, luego
$\chi_\pi=\sum_cm_c\chi_{\rho_c}$; se aplica (a). (c) Por el mismo desarrollo y (a),
$\langle\chi_\pi,\chi_\pi\rangle_G=\sum_cm_c^2$, que es $\dim\pi(G)'$ por el teorema 4.3. $\square$

Una segunda demostración de (c) es la forma en que la biblioteca muestrea conmutantes. El **twirl**
$P(Y)=\frac1{|G|}\sum_g\pi(g)Y\pi(g)^*$ es idempotente con imagen $\pi(G)'$ (capítulo 2, ejercicio 3) y
autoadjunto para $\langle Y,Z\rangle_{HS}=\operatorname{Tr}(Y^*Z)$ (sustitúyase $g\to g^{-1}$), luego es la proyección
ortogonal sobre $\pi(G)'$, y $\dim\pi(G)'$ es su traza. La aplicación $Y\mapsto AYB$ tiene traza
$\operatorname{Tr}A\operatorname{Tr}B$ (sus entradas diagonales son $A_{ii}B_{jj}$), lo que da $\frac1{|G|}\sum_g|\chi_\pi(g)|^2$.

Rehacemos con caracteres el ejemplo de la sección 2.4. Una matriz de permutación tiene como traza el número
de vectores de la base que deja fijos, y $\pi(\sigma)$ fija $|b_1b_2b_3\rangle$ si y solo si $b$ es constante en
los ciclos de $\sigma$ (lema B.1(c)), luego $\chi_\pi(\sigma)=2^{\#\text{cycles}(\sigma)}$. Tenemos a mano tres
representaciones unitarias de $S_3$. Sea $P_\sigma e_i=e_{\sigma(i)}$ la representación de permutación sobre $\mathbb C^3$
(un homomorfismo en las matrices de permutación). La representación **trivial** es $\sigma\mapsto1$. El **signo** es
$\operatorname{sgn}\sigma=\det P_\sigma$, un homomorfismo en $\{\pm1\}\subseteq U(1)$ porque el determinante es
multiplicativo y el determinante de una matriz de permutación es $\pm1$. La recta $\mathbb C(1,1,1)$ es invariante,
luego también lo es $E=(1,1,1)^\perp$ (lema 4.1(a)); la representación **estándar** es la restricción de $P$ a $E$.
En una base ortonormal adaptada a $\mathbb C(1,1,1)\oplus E$, $P_\sigma=1\oplus\operatorname{std}(\sigma)$, luego
$\chi_{\rm std}(\sigma)=\operatorname{Tr}P_\sigma-1$, el número de puntos fijos menos uno. Los seis elementos de $S_3$ son
la identidad, tres transposiciones y dos $3$-ciclos, con $3,2,1$ ciclos y $3,1,0$ puntos fijos; una transposición
intercambia dos columnas de $1_3$, luego su signo es $-1$, y un $3$-ciclo es producto de dos transposiciones
($(1\,2\,3)=(1\,3)(1\,2)$), luego su signo es $+1$. Esto da la tabla

| | identidad ($1$) | transposiciones ($3$) | $3$-ciclos ($2$) |
|---|---|---|---|
| $\chi_{\rm triv}$ | $1$ | $1$ | $1$ |
| $\operatorname{sgn}$ | $1$ | $-1$ | $1$ |
| $\chi_{\rm std}$ | $2$ | $0$ | $-1$ |
| $\chi_\pi$ | $8$ | $4$ | $2$ |

(así, $\operatorname{sgn}\sigma=(-1)^{3-\#\text{cycles}(\sigma)}$; todos los caracteres son reales). Con
$\langle\varphi,\psi\rangle_G$ sumado clase por clase:
$\langle\chi_{\rm triv},\chi_{\rm triv}\rangle_G=\frac{1+3+2}6=1$, $\langle\operatorname{sgn},\operatorname{sgn}\rangle_G=\frac{1+3+2}6=1$ y
$\langle\chi_{\rm std},\chi_{\rm std}\rangle_G=\frac{4+0+2}6=1$, luego las tres representaciones son irreducibles por el
corolario 4.5(c), y son no equivalentes dos a dos porque sus caracteres difieren. El corolario 4.5(b) da las
multiplicidades en $\pi$:
$$
$$
Así, la representación trivial y la estándar son equivalentes a constituyentes $\rho_c$ de $\pi$ (teorema 4.3) con
$m_c=4$ y $m_c=2$; suman $4\cdot1+2\cdot2=8=\sum_ck_cm_c$ dimensiones, luego no hay ningún otro constituyente:
$k=(1,2)$, $m=(4,2)$, y $\dim\pi(S_3)'=4^2+2^2=\frac16(64+3\cdot16+2\cdot4)=20$ (corolario 4.5(c)). El código
recalcula la tabla: obtiene $P_\sigma$ como la restricción de $\pi(\sigma)$ al subespacio invariante generado por
$|100\rangle,|010\rangle,|001\rangle$ (sobre el que $\pi(\sigma)$ permuta las posiciones del único $1$, lema B.1(a)), y
$\#\text{cycles}(\sigma)=\dim\ker(P_\sigma-1)$ (los vectores constantes en los ciclos).

```python
import math
import torch
import torch_vn_algebra.decompose as dc

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    return (torch.as_tensor(X) - torch.as_tensor(Y)).abs().max().item() < tol

G = dc.permutation_representation(3, 2)                  # (6, 8, 8): pi(sigma) permutes three qubits
chi_pi = G.diagonal(dim1=-2, dim2=-1).sum(-1).real       # characters Tr pi(sigma)
P3 = G[:, [4, 2, 1]][:, :, [4, 2, 1]].real              # on |100>, |010>, |001>: S_3 permuting C^3
cycles = 3 - torch.linalg.matrix_rank(P3 - torch.eye(3, dtype=torch.float64)).double()   # dim ker(P - 1)
assert close(chi_pi, 2.0 ** cycles)
irreps = {'trivial': torch.ones(6, dtype=torch.float64), 'sign': (-1.0) ** (3 - cycles),
          'standard': P3.diagonal(dim1=-2, dim2=-1).sum(-1) - 1}                 # fixed points - 1
chis = torch.stack(list(irreps.values()))
assert close(chis @ chis.T / 6, torch.eye(3, dtype=torch.float64))    # Corollary 4.5(a): all irreducible
mult = chis @ chi_pi / 6                                               # Corollary 4.5(b)
print(dict(zip(irreps, mult.tolist())), (chi_pi ** 2).mean().item())
assert close(mult, torch.tensor([4.0, 0.0, 2.0], dtype=torch.float64))
assert dc.commutant_dimension(G) == 20                                # Corollary 4.5(c): M_4 (+) M_2
```

Es la descomposición de la sección 2.4, obtenida ahora a partir de una tabla de seis números por representación.

## 4.3 De los generadores al álgebra: un algoritmo aleatorizado

El teorema 4.3 describe la descomposición, pero no la calcula: su demostración necesita los subespacios
invariantes irreducibles. Murota, Kanno, Kojima y Kojima (*Japan J. Indust. Appl. Math.* 27 (2010) 125–160;
véase también Maehara y Murota, *ibid.* 263–293) observaron que los subespacios propios de un *elemento
aleatorio del conmutante* los encuentran. Sea $\mathcal A$ el álgebra de von Neumann generada por unas matrices
dadas $S_1,\dots,S_J$ y sus adjuntas, con una descomposición $V^*\mathcal AV=\bigoplus_cM_{k_c}\otimes1_{m_c}$,
$V^*\mathcal A'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$ (teorema 2.6); solo el análisis usa $V$.

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

**Proposición 4.7 (genericidad).** Sean $K$ e $Y$ elementos aleatorios independientes de $\mathcal A'_h$ y de
$\mathcal A'$ cuyas leyes tienen densidades respecto de la medida de Lebesgue en estos espacios vectoriales
reales. Con probabilidad uno: (a) se cumple la hipótesis del lema 4.6; (b) $\Pi_iY\Pi_j\neq0$ para cualesquiera
dos proyecciones espectrales $\Pi_i,\Pi_j$ de $K$ cuyas imágenes estén en el mismo sumando.

*Demostración.* Usamos lo siguiente: *un polinomio real $p\not\equiv0$ sobre $\mathbb R^d$ solo se anula en un
conjunto de medida de Lebesgue nula.* Para $d=1$ tiene un número finito de raíces. Para $d>1$ escribimos
$p=\sum_jp_j(x')t^j$ con $x'\in\mathbb R^{d-1}$ y algún $p_j\not\equiv0$; fuera del conjunto nulo $\{p_j=0\}$
(inducción) la sección $\{t:p(x',t)=0\}$ es finita, y el teorema de Tonelli da medida $0$.

(a) La matriz hermítica $\kappa=\bigoplus_c\kappa_c\in M_N$ depende de forma real-lineal de $K$. El determinante
de Gram $\delta(K)=\det[\operatorname{Tr}\kappa^{i+j}]_{i,j=0}^{N-1}$ de $1,\kappa,\dots,\kappa^{N-1}$ para
$\langle\cdot,\cdot\rangle_{HS}$ es un polinomio real en las coordenadas de $K$, no nulo si y solo si estas potencias
son linealmente independientes. Como $p(\kappa)=U\operatorname{diag}(p(\lambda_l))U^*$ para
$\kappa=U\operatorname{diag}(\lambda)U^*$, un polinomio no nulo de grado $<N$ anula a $\kappa$ si y solo si se anula
en $\operatorname{spec}\kappa$, lo cual es posible si y solo si $\kappa$ tiene menos de $N$ valores propios distintos
(tómese $\prod_\mu(t-\mu)$ sobre ellos). Así, $\delta(K)\neq0$ es la hipótesis del lema 4.6, y $\delta\not\equiv0$:
es no nulo en $V(\bigoplus_c1\otimes\operatorname{diag}(\text{reales distintos}))V^*\in\mathcal A'_h$. Por tanto,
$\{\delta=0\}$ es un conjunto nulo.

(b) En el conjunto abierto donde se cumple (a), las proyecciones espectrales son polinomios en $K$ cuyos
coeficientes son funciones continuas de los valores propios (demostración del teorema 1.5(a)). Fijemos un tal
$K$. Por el lema 4.6, $V^*\Pi_iV$ y $V^*\Pi_jV$ proyectan sobre $\mathbb C^{k_c}\otimes e_i$ y
$\mathbb C^{k_c}\otimes e_j$, y para $V^*YV=\bigoplus_d1\otimes y_d$, $\Pi_iY\Pi_j=0$ si y solo si
$\langle e_i,y_ce_j\rangle=0$. La aplicación $Y\mapsto\langle e_i,y_ce_j\rangle$ es real-lineal
$\mathcal A'\to\mathbb C$ y no nula (su valor en $V(1\otimes|e_i\rangle\langle e_j|)V^*$ es $1$), luego su núcleo es
un conjunto nulo. Como $Y$ es independiente de $K$, el teorema de Tonelli da probabilidad $0$. $\square$

**Lema 4.8 (alineación de las copias).** Sean $W_i,W_j$ subespacios invariantes sobre los que $\mathcal A$ actúa
irreduciblemente, con bases ortonormales $Q_i,Q_j$ y restricciones $\sigma_i(a)=Q_i^*aQ_i$, y sea $Y\in\mathcal A'$.
Entonces $T=Q_i^*YQ_j\in\operatorname{Hom}(\sigma_j,\sigma_i)$. Si $T\neq0$, entonces $T=\lambda U$ con $\lambda>0$, $U$
unitaria, y $Q_j'=Q_jU^*$ satisface $Q_j'^*aQ_j'=\sigma_i(a)$ para todo $a\in\mathcal A$. Si $W_i,W_j$ están en
sumandos distintos ($W_i\subseteq P_c\mathbb C^n$, $W_j\subseteq P_d\mathbb C^n$, $c\neq d$), entonces $T=0$.

*Demostración.* Como en la demostración del teorema 4.3, $aQ_j=Q_j\sigma_j(a)$ y $Q_i^*a=\sigma_i(a)Q_i^*$, luego
$T\sigma_j(a)=Q_i^*YaQ_j=Q_i^*aYQ_j=\sigma_i(a)T$. El lema 4.2 da $T=\lambda U$ y
$U\sigma_j(a)U^*=\sigma_i(a)$, que es $Q_j'^*aQ_j'$. Finalmente, $P_c\in\mathcal A\cap\mathcal A'$, $Q_i=P_cQ_i$ y
$P_cQ_j=0$, luego $T=Q_i^*P_cYQ_j=Q_i^*YP_cQ_j=0$. $\square$

El algoritmo (`decompose.wedderburn`, con las interfaces `from_group` y `from_generators`) es ahora:

1. *Muestrear el conmutante.* Para todos los elementos de un grupo finito, $Y=P(Z)$ y $K=P(Z+Z^*)$ con el twirl
   $P$ y matrices gaussianas $Z$: imágenes de leyes gaussianas por aplicaciones lineales sobreyectivas sobre
   $\mathcal A'$ y $\mathcal A'_h$, y por tanto con densidades. Para generadores, $Y,K$ son combinaciones
   gaussianas de una base de $\mathcal A'=\{X:[S_j,X]=[S_j^*,X]=0\}$, el núcleo de una aplicación lineal de $M_n$.
2. *Separar.* Los subespacios propios de $K$ son irreducibles (lema 4.6, proposición 4.7(a)); sus dimensiones
   son los $k_c$.
3. *Agrupar y alinear.* Los bloques $Q_i^*YQ_j$ son no nulos exactamente dentro de un sumando (proposición 4.7(b),
   lema 4.8); las clases dan los sectores, sus tamaños los $m_c$, y los factores polares alinean las copias.
4. *Ensamblar* $V$ en el orden de Kronecker, de modo que $V^*aV=\bigoplus_ca_c\otimes1_{m_c}$.

En coma flotante, los valores propios cercanos de copias distintas se mantienen juntos y se separan mediante la
compresión de un nuevo elemento aleatorio del conmutante (docstring del módulo `decompose`). En cualquier caso,
el resultado se certifica al final:

**Proposición 4.9 (certificado).** Sea $\mathcal A$ generada por $S_1,\dots,S_J$ y sus adjuntas, y
$\mathcal B=V(\bigoplus_cM_{k_c}\otimes1_{m_c})V^*$ con $V$ unitaria. Si todas las $S_j\in\mathcal B$ y
$\dim\mathcal A'=\sum_cm_c^2$, entonces $\mathcal A=\mathcal B$.

*Demostración.* $\mathcal B$ es una $*$-álgebra que contiene todas las $S_j,S_j^*$, luego $\mathcal A\subseteq\mathcal B$ y
$\mathcal B'\subseteq\mathcal A'$. El teorema 4.3, aplicado a la representación de $\mathcal B$ sobre $\mathbb C^n$ (las
aplicaciones de bloque $b\mapsto b_c$ son irreducibles, y no equivalentes, pues la $c$-ésima proyección central se
aplica en $1$ por una de ellas y en $0$ por las demás), da $\dim\mathcal B'=\sum_cm_c^2=\dim\mathcal A'$. Así,
$\mathcal B'=\mathcal A'$ y $\mathcal B=\mathcal B''=\mathcal A''=\mathcal A$. $\square$

`Decomposition.verify()` informa de las distancias relativas de los generadores (y de nuevos elementos del
conmutante) a $\mathcal B$ (resp. $\mathcal B'$), y compara $\sum_cm_c^2$ con $\dim\mathcal A'$ (a partir de los
caracteres para un grupo, corolario 4.5(c); a partir del núcleo para generadores). Lo ejecutamos sobre $S_3$ y tres
cúbits, y luego miramos dentro: los subespacios propios de un $K$ aleatorio y los operadores de entrelazamiento de
un $Y$ aleatorio.

```python
D = dc.from_group(G)                                     # the twirl samples pi(S_3)'
print(D)
assert D.k_factors == (1, 2) and D.multiplicities == (4, 2)            # Section 2.4: k = (1, 2), m = (4, 2)
assert D.verify()['ok']
rho = D.to_algebra(G).matrix                             # (6, 2, 2, 2): the blocks rho_c(sigma), padded
assert close(rho[:, 0, 0, 0], torch.ones(6))                           # sector 1: trivial representation
assert close(rho[:, 1].diagonal(dim1=-2, dim2=-1).sum(-1), chis[2])    # sector 2: the standard one
schur = torch.einsum('gij,gkl->ijkl', rho[:, 1], rho[:, 1].conj()) / 6
assert close(schur.reshape(4, 4), torch.eye(4) / 2)                    # Theorem 4.4 with k = 2
assert close(D.from_algebra(D.to_algebra(G)), G)

Z = torch.randn(8, 8, dtype=dt)
K = dc.twirl(G, Z + Z.mH)                                # a random element of pi(S_3)'_h
w, Q = torch.linalg.eigh(K)
cuts = [0] + [i + 1 for i in range(7) if w[i + 1] - w[i] > 1e-8] + [8]
spaces = [Q[:, a:b] for a, b in zip(cuts[:-1], cuts[1:])]
assert sorted(W.shape[1] for W in spaces) == [1, 1, 1, 1, 2, 2]      # Lemma 4.6: m_c pieces of size k_c
assert all(close(G @ W, W @ (W.mH @ G @ W)) for W in spaces)           # invariant subspaces
Y = dc.twirl(G, torch.randn(8, 8, dtype=dt))
Wa, Wb = [W for W in spaces if W.shape[1] == 2]
s = torch.linalg.svdvals(Wa.mH @ Y @ Wb)                               # Lemma 4.8: a multiple of a unitary
assert s[1] > 1e-3 and abs(s[0] - s[1]) < 1e-10
assert all(close(W.mH @ Y @ Wa, 0) for W in spaces if W.shape[1] == 1)   # different summands: zero
```

El conmutante domina el coste: $O(|G|n^3)$ por cada muestra del twirl; para generadores, la biblioteca reduce el
problema del núcleo, con un elemento aleatorio de $\mathcal A$, a $N=\sum_ck_cm_c^2$ incógnitas, pero sigue
necesitando memoria $O(n^2N)$, unos $10$ GB para el espín total de ocho espines $\frac12$ ($n=256$, $N=4900$). Por
eso, para cadenas de espines, la sección 4.6 usa el módulo `su2`, que conoce la respuesta de antemano.

## 4.4 $U(1)$ y sectores de carga

**Proposición 4.10.** Sea $Q$ una matriz hermítica sobre $\mathbb C^n$ con valores propios enteros $q$ (las
**cargas**), proyecciones espectrales $P_q$ y $d_q=\operatorname{rank}P_q$, y sea $U_\theta=e^{i\theta Q}$; escribimos
$U(1)'=\{U_\theta:\theta\in\mathbb R\}'$ y $U(1)''$ para su conmutante.

1. $\theta\mapsto U_\theta$ es una representación unitaria de $U(1)=\mathbb R/2\pi\mathbb Z$, y toda representación
   unitaria continua de $U(1)$ sobre $\mathbb C^n$ es de esta forma.
2. Para $X\in M_n$ son equivalentes: (i) $XU_\theta=U_\theta X$ para todo $\theta$; (ii) $[X,Q]=0$;
   (iii) $X=\sum_qP_qXP_q$. Por tanto, $U(1)'\cong\bigoplus_qM_{d_q}(\mathbb C)$ (los operadores diagonales por
   bloques) y $U(1)''=\operatorname{span}\{P_q\}$: todos los $k_c=1$, y las multiplicidades son los $d_q$.
3. (*Carga conservada.*) Para un $H$ hermítico, $e^{itH}Qe^{-itH}=Q$ para todo $t$ si y solo si $[H,Q]=0$, es
   decir, si y solo si $H\in U(1)'$. Entonces todo peso de sector $\operatorname{Tr}(P_q\rho(t))$ de
   $\rho(t)=e^{-itH}\rho e^{itH}$ es constante.

*Demostración.* (1) $U_\theta=\sum_qe^{iq\theta}P_q$ es unitaria, $U_\theta U_{\theta'}=U_{\theta+\theta'}$ (teorema 1.5(c)) y
$U_{2\pi}=1$. Recíprocamente, para una representación continua, $(U_\theta)_{\theta\ge0}$ es un semigrupo continuo
de matrices, luego $U_\theta=e^{\theta B}$ para $\theta\ge0$ por el lema 8.5 (demostrado en el capítulo 8; la
demostración solo usa la dimensión finita), y para $\theta<0$, $U_\theta=U_{-\theta}^{-1}=e^{\theta B}$. Derivando
$U_\theta^*U_\theta=1$ en $0$ se obtiene $B^*+B=0$, luego $B=iQ$ con $Q$ hermítica, y $U_{2\pi}=1$ obliga a
$e^{2\pi iq}=1$, $q\in\mathbb Z$, para todo valor propio. (2) (i)$\Rightarrow$(ii): se deriva en $\theta=0$.
(ii)$\Rightarrow$(iii): $X$ conmuta con $P_q=f(Q)$ (teorema 1.5(d)), luego $P_qXP_{q'}=XP_qP_{q'}=0$ para
$q\neq q'$. (iii)$\Rightarrow$(i): $XU_\theta=\sum_qe^{iq\theta}P_qXP_q=U_\theta X$. Los operadores diagonales por bloques
forman $\bigoplus_qM_{d_q}\otimes1_1$, cuyo conmutante es $\bigoplus_q1_{d_q}\otimes M_1=\operatorname{span}\{P_q\}$
(teorema 2.6). (3) $\frac d{dt}e^{itH}Qe^{-itH}=ie^{itH}[H,Q]e^{-itH}$ (proposición 8.2) se anula para todo $t$ si y
solo si $[H,Q]=0$. Entonces $e^{itH}$ conmuta con $P_q$ (teorema 1.5(d), dos veces), y
$\operatorname{Tr}(P_q\rho(t))=\operatorname{Tr}(e^{itH}P_qe^{-itH}\rho)=\operatorname{Tr}(P_q\rho)$. $\square$

Las representaciones irreducibles de $U(1)$ son, pues, los caracteres $\theta\mapsto e^{iq\theta}$, $q\in\mathbb Z$
(unidimensionales, como para todo grupo abeliano), y la multiplicidad de la carga $q$ es $d_q$. Las cargas se suman
en los productos tensoriales: $e^{i\theta(Q_1\otimes1+1\otimes Q_2)}=e^{i\theta Q_1}\otimes e^{i\theta Q_2}$ actúa sobre la
imagen de $P_{q_1}\otimes P_{q_2}$ con carga $q_1+q_2$. Esta es la regla de fusión $(q_1,q_2)\mapsto q_1+q_2$ de los
productos tensoriales fusionados del capítulo 7.

Para $n$ cúbits, la carga es el **operador número** $\hat N=\sum_jn_j$, con $n_j=|1\rangle\langle1|$ sobre el
cúbit $j$: $\hat N|b\rangle=|b|\,|b\rangle$ con el peso de Hamming $|b|$. Sus sectores tienen dimensión
$\binom nw$, de modo que los operadores que conservan el peso de Hamming forman
$$
U(1)'\cong\bigoplus_{w=0}^nM_{\binom nw}(\mathbb C),
$$
que es `chains.SpinChain(n).algebra`, el producto fusionado de $n$ copias de $M_1\oplus M_1$ con cargas $0,1$
(bit $1$ = espín hacia arriba, $\hat N=S^z+\frac n2$; en `decompose.total_spin`, en cambio, el espín hacia arriba es
el *primer* vector de la base de $\mathbb C^2$, lo que reetiqueta la base pero no los sectores). Todo hamiltoniano XXZ
$\sum_{\langle ij\rangle}J(s^x_is^x_j+s^y_is^y_j+\Delta s^z_is^z_j)$ conserva $\hat N$, porque el término de salto
$s^x_is^x_j+s^y_is^y_j=\frac12(s^+_is^-_j+s^-_is^+_j)$ traslada un espín hacia arriba de un sitio a otro.

```python
import functools
from torch_vn_algebra import SpinChain

L = 4
Sx, Sy, Sz = dc.total_spin(L)                              # (16, 16) each, site 0 most significant
N_op = Sz + L / 2 * torch.eye(2 ** L, dtype=dt)            # number of up spins
DU = dc.from_generators(N_op)                              # the algebra generated by U(1)
print(DU)
chain = SpinChain(L, boundary='periodic')                  # U(1)' = (+)_w M_binom(L, w)
assert set(DU.k_factors) == {1} and sorted(DU.multiplicities) == sorted(chain.algebra.k_factors)
assert sorted(chain.algebra.k_factors) == sorted(math.comb(L, w) for w in range(L + 1))

sx, sy, sz = dc.spin_matrices(0.5)
def two_site(a, b, i, j, n=L):                             # a on site i, b on site j of n sites
    ops = [a if k == i else b if k == j else torch.eye(2, dtype=dt) for k in range(n)]
    return functools.reduce(torch.kron, ops)

H = sum(two_site(sx, sx, i, j) + two_site(sy, sy, i, j) + 0.5 * two_site(sz, sz, i, j)
        for i, j in chain.bonds)                           # XXZ ring, Delta = 0.5
assert close(H @ N_op, N_op @ H)                           # [H, N] = 0 ...
assert close(DU.from_commutant(DU.to_commutant(H)), H)     # ... so H is block diagonal (Proposition 4.10)
Hc = chain.xxz(J=1.0, Delta=0.5)                           # the same H, built sector by sector
E_sec = torch.cat([torch.linalg.eigvalsh(Hc.matrix[0, N, :k, :k]) for N, k in enumerate(chain.algebra.k_factors)])
assert close(torch.sort(E_sec).values, torch.linalg.eigvalsh(H))
```

> **Significado físico.** La simetría $U(1)$ es la conservación de un número de partículas, de una carga o de
> una magnetización. La proposición 4.10(3) es el teorema de Noether en este contexto: la carga se conserva si y
> solo si el hamiltoniano es simétrico, y entonces la dinámica nunca abandona un sector (teorema 8.3). Si además
> nada puede medir coherencias entre sectores, los sectores son sectores de superselección (capítulo 3).

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
$e^{-i\theta\,\mathbf n\cdot\mathbf S}$ conmutan con $X$ para todos los $\theta,\mathbf n$ si y solo si $X$ conmuta con
las tres $S_a$ (derívese en $\theta=0$; recíprocamente, teorema 1.5(d)), de modo que el conmutante de «la acción de
$SU(2)$» es $\{S_x,S_y,S_z\}'$, y por el teorema 2.2 las rotaciones y las $S_a$ generan la misma álgebra. Para $L$
espines $\frac12$, las rotaciones son las $V^{\otimes L}$, $V\in SU(2)$: toda $V$ de este tipo es alguna
$e^{-i\theta\,\mathbf n\cdot\mathbf s}$, y $e^{-i\theta\,\mathbf n\cdot\mathbf S}=(e^{-i\theta\,\mathbf n\cdot\mathbf s})^{\otimes L}$.

Pongamos $S_\pm=S_x\pm iS_y$, de modo que $S_-=S_+^*$, y $\mathbf S^2=S_x^2+S_y^2+S_z^2$. Las relaciones dan
$$
[S_z,S_\pm]=\pm S_\pm,\qquad[S_+,S_-]=2S_z,\qquad\mathbf S^2=S_-S_++S_z^2+S_z,
$$
y $\mathbf S^2$ conmuta con todas las $S_a$ (por ejemplo,
$[\mathbf S^2,S_z]=S_x[S_x,S_z]+[S_x,S_z]S_x+S_y[S_y,S_z]+[S_y,S_z]S_y=-i\{S_x,S_y\}+i\{S_y,S_x\}=0$).

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

**Corolario 4.12.** Toda representación de espín sobre $\mathbb C^n$ es unitariamente equivalente a
$\bigoplus_s(\text{espín }s)\otimes1_{m_s}$ y, con $d_M=\dim\ker(S_z-M)$,
$$
m_s=d_s-d_{s+1}.
$$

*Demostración.* Teorema 4.3 y teorema 4.11. Cada copia del espín $s'$ contiene exactamente un vector propio de
$S_z$ para el valor propio $M$ si $|M|\le s'$ y $s'-M\in\mathbb Z$, y ninguno en otro caso; así,
$d_M=\sum_{s'\ge|M|,\,s'-M\in\mathbb Z}m_{s'}$, y la diferencia de dos valores consecutivos es $m_s$. $\square$

**Proposición 4.13 ($L$ espines $\frac12$).** El espín total sobre $(\mathbb C^2)^{\otimes L}$ contiene el espín $S$
si y solo si $S\in\{\frac L2,\frac L2-1,\dots\}$, $S\ge0$, con multiplicidad
$$
m_S=\binom L{\frac L2-S}-\binom L{\frac L2-S-1}\qquad\Big(\tbinom L{-1}=0\Big),
$$
de modo que el álgebra generada por $S_x,S_y,S_z$ es $\bigoplus_SM_{2S+1}\otimes1_{m_S}$ y su conmutante
$\bigoplus_S1_{2S+1}\otimes M_{m_S}$.

*Demostración.* $S_z|b\rangle=(|b|-\frac L2)|b\rangle$, luego $d_M=\binom L{L/2+M}=\binom L{L/2-M}$ si
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
de dimensiones $35$ y $14$, los números hallados en el capítulo 10. La biblioteca encuentra lo mismo a partir
únicamente de los tres generadores (sectores ordenados por $(k_c,m_c)$):

```python
S4 = dc.total_spin(4)
assert close(S4[0] @ S4[1] - S4[1] @ S4[0], 1j * S4[2])                # [S_x, S_y] = i S_z
D4 = dc.from_generators(S4)
print(D4)
assert D4.k_factors == (1, 3, 5) and D4.multiplicities == (2, 3, 1)     # M_1, M_3, M_5 with m = 2, 3, 1
assert D4.verify()['ok'] and D4.dims == (35, 14)
m_S = {S: math.comb(4, 2 - S) - (math.comb(4, 1 - S) if S < 2 else 0) for S in (0, 1, 2)}
assert [m_S[(k - 1) // 2] for k in D4.k_factors] == list(D4.multiplicities)    # Proposition 4.13
casimir = D4.to_algebra((S4 @ S4).sum(0))                                  # S^2 = S(S+1) on spin S
for c, k in enumerate(D4.k_factors):
    assert close(casimir.matrix[0, c, :k, :k], (k - 1) / 2 * (k + 1) / 2 * torch.eye(k))
```

¿Qué operadores conmutan con todas las rotaciones? Las permutaciones de los sitios lo hacen: la conjugación por
$\pi(\sigma)$ lleva $s^{(i)}_a$ a $s^{(\sigma(i))}_a$, luego deja fijo $\sum_is^{(i)}_a$. El recíproco es un teorema
clásico.

**Teorema 4.14 (dualidad de Schur–Weyl).** Sea $\pi$ la representación de $S_L$ que permuta los factores de
$(\mathbb C^2)^{\otimes L}$. (a) $\operatorname{span}\pi(S_L)\subseteq\{S_x,S_y,S_z\}'$. (b) Se da la igualdad;
equivalentemente (teorema 2.2), el álgebra generada por el espín total es $\pi(S_L)'$.

La parte (a) se acaba de demostrar. La parte (b) se enuncia sin demostración para $L$ general; es válida para
$(\mathbb C^d)^{\otimes L}$ y la acción diagonal de $U(d)$ (W. Fulton y J. Harris, *Representation
Theory*, Springer 1991, lección 6; R. Goodman y N. Wallach, *Symmetry, Representations, and Invariants*,
Springer 2009, sección 9.1). Para $L\le5$ la comprobamos: por (a), basta con que los dos espacios tengan la misma
dimensión. Para $L=2$ es el lema 5.8 con $k=2$, pues $\{S_a\}'=\{V\otimes V:V\in SU(2)\}'=\{V\otimes V:V\in U(2)\}'$
(una fase actúa como un escalar).

```python
for n in range(2, 6):
    comm = dc.commutant_basis(dc.total_spin(n)).reshape(-1, 4 ** n)          # a basis of {S_a}'
    perm = dc.permutation_representation(n, 2).reshape(-1, 4 ** n)          # all n! permutations
    dims = [torch.linalg.matrix_rank(A).item() for A in (comm, perm, torch.cat([comm, perm]))]
    print(n, dims)
    assert dims[0] == dims[1] == dims[2]                                     # span pi(S_n) = {S_a}'
Dp = dc.from_group(dc.permutation_representation(4, 2))                   # pi(S_4)'' = {S_a}'
print(Dp)
assert sorted(zip(Dp.multiplicities, Dp.k_factors)) == sorted(zip(D4.k_factors, D4.multiplicities))
```

Las dimensiones impresas son $2,5,14,42$, los números de Catalan $\binom{2L}L/(L+1)$. La última línea es la
dualidad en acción: $\pi(S_4)''$ tiene bloques $k=(1,2,3)$ con multiplicidades $(5,1,3)$, los pares $(k_c,m_c)$ del
espín total con los papeles intercambiados. Por el teorema 4.14(b) y la proposición 2.7, el espacio de
multiplicidad del espín $S$ porta una representación irreducible de $S_L$, aquí de dimensiones $1,3,2$ para
$S=2,1,0$.

> **Significado físico.** Un multiplete de espín $S$ consta de $2S+1$ estados de igual energía para todo
> hamiltoniano invariante por rotaciones, que se distinguen por $S^z=M$. Dos espines $\frac12$ forman un triplete y
> un singlete; cuatro espines forman un quintuplete, tres tripletes y dos singletes. Los espacios de multiplicidad
> de los dos singletes y de los tres tripletes son el cúbit libre de decoherencia y el cútrit sin ruido del
> capítulo 10.

## 4.6 Sistemas de espines

El **modelo de Heisenberg** sobre $L$ espines $\frac12$ es $H=\sum_{i<j}J_{ij}\,\mathbf s^{(i)}\cdot\mathbf s^{(j)}$ con
acoplamientos reales $J_{ij}$.

**Proposición 4.15.** (a) (*Identidad de Dirac*) $\mathbf s^{(i)}\cdot\mathbf s^{(j)}=\frac12F_{ij}-\frac14$, donde $F_{ij}$
intercambia los espines $i$ y $j$. Por tanto, todo hamiltoniano de Heisenberg pertenece a
$\operatorname{span}\pi(S_L)\subseteq\{S_a\}'$, y $V^*HV=\bigoplus_S1_{2S+1}\otimes H_S$ con bloques $H_S\in M_{m_S}$ en la
descomposición de la proposición 4.13: cada valor propio de $H_S$ es un nivel $(2S+1)$ veces degenerado de $H$.
(b) El espacio de los **vectores de peso máximo** $\ker S_+\cap\ker(S_z-S)$ tiene dimensión $m_S$, es invariante
por $H$, y $H$ actúa sobre él como $H_S$ (salvo un cambio de base unitario).

*Demostración.* (a) Sobre dos espines,
$\mathbf s^{(1)}\cdot\mathbf s^{(2)}=\frac12(\mathbf S^2-\mathbf s^{(1)2}-\mathbf s^{(2)2})=\frac12\mathbf S^2-\frac34$
(teorema 4.11 con $s=\frac12$). Sobre el triplete, $\mathbf S^2=2$ y $F=1$ (está generado por
$S_-^j|{\uparrow\uparrow}\rangle$, y $F$ conmuta con $S_-$ y deja fijo $|{\uparrow\uparrow}\rangle$); sobre el singlete,
$\mathbf S^2=0$ y $F=-1$. Ambos miembros coinciden sobre estos cuatro vectores. Para sitios $i,j$ de $L$, se toma el
producto tensorial con la identidad. El resto es el teorema 4.14(a) y el teorema 4.3. (b) En la descomposición del
corolario 4.12, un vector con $S_z=S$ es $\sum_{S'\ge S}|S',S\rangle\otimes w_{S'}$, y $S_+|S',S\rangle\neq0$ para
$S'>S$ por el teorema 4.11; así, $\ker S_+\cap\ker(S_z-S)=|S,S\rangle\otimes\mathbb C^{m_S}$, sobre el que
$1\otimes H_S$ actúa como $H_S$. $\square$

La parte (b) es la forma en que el módulo `su2` (cadenas de espines, capítulo 9) construye $H_S$ sin
descomposición numérica: se construye una base ortonormal de los vectores de peso máximo mediante el acoplamiento
de Clebsch–Gordan, un sitio cada vez, dentro del sector $U(1)$ $\hat N=\frac L2+S$ de `SpinChain`. Comprobamos (a)
sobre un anillo de cuatro espines con `D4` y después comparamos los tamaños de los bloques para ocho espines.

```python
def swap(i, j, n):
    p = list(range(n))
    p[i], p[j] = j, i
    return dc.permutation_representation(n, 2, [p])[0]

ring = [(i, (i + 1) % 4) for i in range(4)]
H4 = sum(swap(i, j, 4) / 2 - torch.eye(16, dtype=dt) / 4 for i, j in ring)        # Dirac's identity
assert close(H4, sum(two_site(a, a, i, j) for a in (sx, sy, sz) for i, j in ring))
h = D4.to_commutant(H4)                                   # blocks H_S of sizes m_S = 2, 3, 1
assert close(D4.from_commutant(h), H4)                    # H4 lies in {S_a}'
levels = torch.cat([torch.linalg.eigvalsh(h.matrix[0, c, :m, :m]).repeat_interleave(k)
                    for c, (k, m) in enumerate(zip(D4.k_factors, D4.multiplicities))])
assert close(torch.sort(levels).values, torch.linalg.eigvalsh(H4))   # level of H_S repeated 2S+1 times

from torch_vn_algebra import su2
u1, sym = SpinChain(8, boundary='periodic'), su2.SU2Chain(8, boundary='periodic')
H_u1, H_su2 = u1.xxz(J=1.0, Delta=1.0), sym.heisenberg(J=1.0)
E_u1 = torch.sort(torch.cat([torch.linalg.eigvalsh(H_u1.matrix[0, N, :k, :k])
                             for N, k in enumerate(u1.algebra.k_factors)])).values
E_su2, S_lev, _ = sym.full_spectrum(H_su2, return_quantum_numbers=True)   # level of block S, 2S+1 times
E_su2, S_lev = E_su2[0], S_lev[0]
assert close(E_u1, E_su2)
assert S_lev[0] == 0 and (S_lev[1:4] == 1).all() and E_su2[4] - E_su2[3] > 1e-3     # singlet, then triplet
assert close(sym.restrict(H_u1).matrix, H_su2.matrix)      # Proposition 4.15(b)
flops = lambda ks: sum(k ** 3 for k in ks)
print(u1.algebra.k_factors, sym.multiplicities, flops(u1.algebra.k_factors), flops(sym.multiplicities))
print(f"ground state energy {E_su2[0].item():.6f}, first excited {E_su2[1].item():.6f} (x3)")
```

Para el anillo de ocho espines, los bloques de $U(1)$ tienen tamaños $1,8,28,56,70,56,28,8,1$, y los bloques de
$SU(2)$ tamaños $m_S=14,28,20,7,1$ para $S=0,\dots,4$. El coste de diagonalizar, $\sum k^3$, baja de $739162$ a
$33040$, en un factor $22$; la dimensión del álgebra, $\sum k^2$, de $\binom{16}8=12870$ a $1430$. El estado
fundamental es un singlete con energía $-3.651093$, y el primer nivel excitado un triplete en $-3.128419$.

> **Significado físico.** El modelo de Heisenberg describe la interacción de intercambio de espines electrónicos
> localizados en aislantes magnéticos; $J_{ij}>0$ favorece vecinos antiparalelos. Su espectro consta de
> multipletes, y un campo magnético $h\,S^z$ desdobla cada multiplete en $2S+1$ niveles equidistantes (ejercicio 4).
> Para acoplamientos $J_{ij}>0$ entre las dos subredes de una red bipartita, $J_{ij}\le0$ dentro de ellas, y
> subredes del mismo tamaño, el estado fundamental es un singlete (teorema de Lieb–Mattis, E. Lieb y D. Mattis,
> *J. Math. Phys.* 3 (1962) 749; se enuncia sin demostración). El anillo de ocho espines de arriba es un ejemplo.

## 4.7 Con qué comparar

Para grupos finitos, GAP y SageMath calculan tablas de caracteres de forma exacta, y el corolario 4.5 da entonces
las multiplicidades; RepLAB (MATLAB; D. Rosset, F. Montealegre-Mora, J.-D. Bancal) descompone numéricamente
representaciones de grupos finitos y compactos, también con elementos aleatorios del conmutante, y
SymbolicWedderburn.jl (Julia) calcula descomposiciones de Wedderburn para la reducción por simetría de programas
semidefinidos. En física, QuSpin y otros códigos de diagonalización exacta diagonalizan por bloques las simetrías
$U(1)$ y de red, ITensor y TeNPy manejan cargas abelianas en redes tensoriales, y TensorKit.jl trata $SU(2)$ con
árboles de fusión. El módulo `decompose` de la biblioteca es un método denso para matrices pequeñas ($n$ de hasta
unos pocos cientos); acepta cualquier lista de matrices (elementos de un grupo o generadores), certifica su
resultado (proposición 4.9) y devuelve el álgebra y su conmutante como objetos `TypeIAlgebra` con `to_algebra`,
`from_algebra`, `to_commutant`, `from_commutant`, de modo que los capítulos 3 y 5–8 se aplican de inmediato. Para
cadenas de espín $\frac12$, `su2.SU2Chain` da directamente los bloques de $SU(2)$, para todos los espines totales,
hasta $L=16$.

## Ejercicios

1. Demuestre que, para la representación de $S_L$ sobre $(\mathbb C^d)^{\otimes L}$,
   $\chi_\pi(\sigma)=d^{\#\text{cycles}(\sigma)}$ y $\dim\pi(S_L)'=\binom{d^2+L-1}L$. Compruebe $20$ (sección 2.4) y
   $35$ ($=\dim\{S_a\}''$ para cuatro espines).
2. (*Clebsch–Gordan.*) Demuestre que $j_1\otimes j_2\cong\bigoplus_jj$, $j=|j_1-j_2|,|j_1-j_2|+1,\dots,j_1+j_2$, para
   la representación de espín $S^{(j_1)}_a\otimes1+1\otimes S^{(j_2)}_a$ —la regla de fusión no abeliana de $SU(2)$
   que menciona el capítulo 7— y compruebe $1\otimes\frac12$ y $1\otimes1$ con `from_generators`.
3. Sea $\theta/2\pi$ irracional. Demuestre que $\{U_\theta\}'=U(1)'$ para una sola $U_\theta=e^{i\theta Q}$, pero que
   $\{e^{i\pi\hat N}\}'\cong M_{2^{n-1}}\oplus M_{2^{n-1}}$ para $n$ cúbits.
4. (*Desdoblamiento Zeeman.*) Sea $H$ un hamiltoniano de Heisenberg sobre $L\ge1$ espines y $H_h=H+hS_z$. Demuestre
   que $H_h$ conmuta con $S_z$ y con $\mathbf S^2$ pero, para $h\neq0$, no con $S_x$, y que sus valores propios son
   $E+hM$, $M=-S,\dots,S$, para los valores propios $E$ de $H_S$. Compruébelo para el anillo de ocho espines.

<details>
<summary>Soluciones</summary>

1. Como para $d=2$: $\pi(\sigma)$ es una matriz de permutación que deja fijo $|b_1\dots b_L\rangle$ si y solo si $b$
   es constante en los ciclos, y hay $d^{\#\text{cycles}}$ tales $b$. Por el corolario 4.5(c),
   $\dim\pi(S_L)'=\frac1{L!}\sum_\sigma(d^2)^{\#\text{cycles}(\sigma)}=\frac{x(x+1)\cdots(x+L-1)}{L!}$ con $x=d^2$
   (capítulo 5, ejercicio 2), que es $\binom{d^2+L-1}L$: $\binom63=20$ y $\binom74=35$. Que $35$ sea también
   $\dim\{S_a\}''$ es el teorema 4.14 para $L=4$.

```python
for n, value in ((3, 20), (4, 35)):
    assert dc.commutant_dimension(dc.permutation_representation(n, 2)) == math.comb(n + 3, n) == value
```

2. Sea $j_1\ge j_2$. $S_z$ tiene el valor propio $m_1+m_2$ sobre $|m_1\rangle\otimes|m_2\rangle$, luego para
   $M\in j_1+j_2+\mathbb Z$, $d_M$ es el número de $m_2\in\{-j_2,\dots,j_2\}$ con $|M-m_2|\le j_1$: los $2j_2+1$ si
   $|M|\le j_1-j_2$, y $j_1+j_2-|M|+1$ si $j_1-j_2\le|M|\le j_1+j_2$ (para $M\ge0$ la condición es
   $m_2\ge M-j_1$), y $0$ más allá. Por el corolario 4.12, $m_j=d_j-d_{j+1}$ vale $1$ para $j_1-j_2\le j\le j_1+j_2$ y
   $(2j_2+1)-(2j_2+1)=0$ para $j<j_1-j_2$.

```python
def coupled(a, b):                                        # S^(a) (x) 1 + 1 (x) S^(b)
    Sa, Sb = dc.spin_matrices(a), dc.spin_matrices(b)
    ea, eb = torch.eye(Sa.shape[-1], dtype=dt), torch.eye(Sb.shape[-1], dtype=dt)
    return torch.stack([torch.kron(x, eb) + torch.kron(ea, y) for x, y in zip(Sa, Sb)])
for a, b, ks in ((1.0, 0.5, (2, 4)), (1.0, 1.0, (1, 3, 5))):  # 1 (x) 1/2 = 1/2 (+) 3/2, 1 (x) 1 = 0 (+) 1 (+) 2
    D_cg = dc.from_generators(coupled(a, b))
    assert D_cg.k_factors == ks and set(D_cg.multiplicities) == {1}
```

3. $U_\theta=\sum_qe^{iq\theta}P_q$, y $e^{iq\theta}=e^{iq'\theta}$ si y solo si $(q-q')\theta\in2\pi\mathbb Z$, lo que para
   $\theta/2\pi$ irracional obliga a $q=q'$. Así, los subespacios propios de $U_\theta$ son las imágenes de los $P_q$, y
   un operador conmuta con la matriz normal $U_\theta$ si y solo si conserva sus subespacios propios (teorema 1.5(d) e
   interpolación de Lagrange como en la sección 2.2), es decir, si y solo si $X=\sum_qP_qXP_q$; esto es $U(1)'$ por la
   proposición 4.10. Para $\theta=\pi$, el valor propio $(-1)^w$ solo ve la paridad del peso de Hamming: dos
   subespacios propios de dimensión $2^{n-1}$ (para $n\ge1$, los subconjuntos pares e impares), y el conmutante es
   $M_{2^{n-1}}\oplus M_{2^{n-1}}$.

```python
for theta, m in ((1.0, (1, 1, 4, 4, 6)), (math.pi, (8, 8))):
    D_theta = dc.from_generators(torch.linalg.matrix_exp(1j * theta * N_op))
    assert D_theta.multiplicities == m and set(D_theta.k_factors) == {1}
```

4. $S_z$ y $\mathbf S^2$ conmutan con $H$ (proposición 4.15(a)) y con $S_z$, luego con $H_h$; y
   $[H_h,S_x]=h[S_z,S_x]=ihS_y\neq0$. En la descomposición de la proposición 4.13,
   $V^*S_zV=\bigoplus_SS^{(S)}_z\otimes1_{m_S}$ y $V^*HV=\bigoplus_S1_{2S+1}\otimes H_S$ (proposición 4.15), luego sobre
   el sumando de espín $S$, $V^*H_hV=hS^{(S)}_z\otimes1+1\otimes H_S$, una suma de operadores que conmutan,
   diagonalizada por los productos de vectores propios, con valores propios $hM+E$.

```python
assert not close(H4 @ Sx - Sx @ H4 + 0.3 * (Sz @ Sx - Sx @ Sz), 0)          # [H4 + 0.3 S_z, S_x] != 0
E_field = sym.full_spectrum(sym.heisenberg(J=1.0, h=0.3))[0]                # E + h M on every multiplet
H_field = u1.xxz(J=1.0, Delta=1.0, h=0.3)
E_u1_field = torch.cat([torch.linalg.eigvalsh(H_field.matrix[0, N, :k, :k])
                        for N, k in enumerate(u1.algebra.k_factors)])
assert close(torch.sort(E_u1_field).values, E_field)
```

</details>

**Siguiente:** el [capítulo 5](05_random_matrices.md) extrae matrices aleatorias en un álgebra: unitarias de Haar de
$\prod_cU(k_c)$, el cálculo de Weingarten —cuyo caso de orden dos se apoya en el conmutante de
$\{V\otimes V\}$, el caso más pequeño de la dualidad de Schur–Weyl— y la estadística de niveles sector por sector.
