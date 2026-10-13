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
el conmutante; y cómo se cuentan los vectores invariantes y las formas cuadráticas invariantes.

En todo el capítulo, $G$ es un grupo finito, y una **representación** de $G$ sobre $\mathbb C^n$ es un homomorfismo
$\pi:G\to M_n(\mathbb C)$ con valores unitarios (no es una restricción: las $\pi(g)$ invertibles son unitarias para el
producto escalar $\sum_g\langle\pi(g)x,\pi(g)y\rangle$, sección 4.1). Su **carácter** es $\chi_\pi(g)=\operatorname{Tr}\pi(g)$.
Las representaciones irreducibles y equivalentes y la **multiplicidad** $m_\rho(\pi)$ de una $\rho$ irreducible en
$\pi$ ($0$ si $\rho$ no aparece) son como en la sección 2.4, que demuestra que $m_\rho(\pi)$ no depende de la
descomposición y demuestra el lema de Schur; el [capítulo 4](04_symmetries.md) trata ambas cosas con detalle
(secciones 4.1–4.2). Sobre matrices rectangulares, $\langle S,T\rangle=\operatorname{Tr}(S^*T)$ es el producto escalar
de **Hilbert–Schmidt**. En el código, un lote de vectores es un tensor de forma `(batch, n)`, sobre el que una
matriz $T$ actúa como `x @ T.T`.

## 12.1 Aplicaciones lineales equivariantes

**Definición.** Sean $\pi_1,\pi_2$ representaciones de $G$ sobre $\mathbb C^{n_1}$, $\mathbb C^{n_2}$. Una matriz
$T\in M_{n_2\times n_1}(\mathbb C)$ es **equivariante** (un **operador de entrelazamiento**) si $T\pi_1(g)=\pi_2(g)T$
para todo $g\in G$. Estas matrices forman un espacio vectorial $\operatorname{Hom}_G(\pi_1,\pi_2)$; para
$\pi_1=\pi_2=\pi$ es el conmutante, $\operatorname{Hom}_G(\pi,\pi)=\pi(G)'$.

Las representaciones más comunes en aprendizaje automático permutan coordenadas. Supongamos que $G$ actúa sobre un
conjunto finito $\Omega$; la **representación de permutación** sobre $\mathbb C^\Omega$ es $\pi(g)e_x=e_{gx}$.
Ejemplos: $S_n$ permutando los $n$ puntos de un conjunto ($\Omega=[n]=\{1,\dots,n\}$), y $S_n$ reetiquetando los
vértices de un grafo, actuando sobre $\Omega=[n]^2$, es decir, sobre la matriz de adyacencia mediante
$X\mapsto P_\sigma XP_\sigma^T$.

**Lema 12.1 (órbitas).** Sean $\pi_i$ las representaciones de permutación de $G$ sobre $\mathbb C^{\Omega_i}$.
Entonces $T\in\operatorname{Hom}_G(\pi_1,\pi_2)$ si y solo si $T_{gy,gx}=T_{yx}$ para todos $g$, $y\in\Omega_2$,
$x\in\Omega_1$. Las matrices indicadoras de las órbitas de $G$ sobre $\Omega_2\times\Omega_1$ forman una base, y
$\dim\operatorname{Hom}_G(\pi_1,\pi_2)$ es el número de esas órbitas.

*Demostración.* $\pi_2(g)T\pi_1(g)^*e_x=\pi_2(g)Te_{g^{-1}x}=\sum_yT_{y,g^{-1}x}e_{gy}$, luego la entrada $(gy,x)$ de
$\pi_2(g)T\pi_1(g)^*$ es $T_{y,g^{-1}x}$. Por tanto, $T$ es equivariante, es decir, $\pi_2(g)T\pi_1(g)^*=T$ para todo
$g$, si y solo si $T_{gy,x}=T_{y,g^{-1}x}$ para todos $g,y,x$, lo cual (sustituyendo $x$ por $gx$) es la condición
del enunciado: $T$ es constante en las órbitas. Las matrices indicadoras de las órbitas tienen soportes disjuntos,
luego son linealmente independientes, y toda función constante en las órbitas es combinación de ellas. $\square$

**Corolario 12.2.** (a) (*DeepSets.*) Para $S_n$ sobre $\mathbb C^n$, $n\ge2$, el conmutante es
$\operatorname{span}\{1,J\}$, con $J$ la matriz cuyas entradas son todas iguales a uno, y
$\mathbb C^n=\mathbb C\mathbf 1\oplus\mathbf 1^\perp$ es suma de dos representaciones irreducibles no equivalentes, la
trivial y la **estándar**, de dimensión $n-1$. (b) (*Redes de grafos.*) Las aplicaciones equivariantes de
$\mathbb C^{[n]^p}$ en $\mathbb C^{[n]^q}$ (tensores de orden $p$ y $q$ sobre $n$ vértices) se cuentan mediante las
particiones de un conjunto de $p+q$ elementos en a lo sumo $n$ bloques; para $n\ge p+q$ su número es el
**número de Bell** $B_{p+q}$ ($B_2=2$, $B_3=5$, $B_4=15$).

*Demostración.* (b) Dos tuplas $x,y\in[n]^{p+q}$ están en la misma órbita si y solo si tienen el mismo **patrón de
igualdades**, la partición de las posiciones $\{1,\dots,p+q\}$ dada por $i\sim j\Leftrightarrow x_i=x_j$. El patrón
es claramente invariante. Recíprocamente, si $x,y$ tienen el mismo patrón, $x_i\mapsto y_i$ es una inyección bien
definida de $\{x_i\}$ en $[n]$, que se extiende a una permutación $\sigma$ con $\sigma x=y$. Toda partición con a lo
sumo $n$ bloques es un patrón (asígnense valores distintos a los bloques). El lema 12.1 permite concluir. (a) Para
$p=q=1$ hay dos órbitas sobre $[n]^2$, la diagonal y su complemento, con indicadoras $1$ y $J-1$. Así,
$\dim\pi(S_n)'=2=\sum_cm_c^2$ (teorema 2.6, proposición 2.7), lo que obliga a que haya dos sectores con $m_c=1$, es
decir, dos sumandos irreducibles no equivalentes, y $\pi(S_n)'$ está generado por sus proyecciones $P_1,P_2$. La
proyección sobre la recta invariante $\mathbb C\mathbf 1$ está en $\pi(S_n)'$ (lema 2.1), luego es $P_1$ o $P_2$: la
recta, que porta la representación trivial, es un sector, y $\mathbf 1^\perp$ es el otro. $\square$

La capa DeepSets $x\mapsto\lambda x+\gamma\,(\mathbf 1^Tx)\mathbf 1$ (Zaheer et al. 2017) es (a); las bases de
Maron et al. (2019) para redes de grafos son (b). El código cuenta las órbitas por fuerza bruta, las compara con el
conmutante calculado numéricamente y descompone ambas representaciones de $S_4$:

```python
import itertools
import torch
import torch_vn_algebra.decompose as dc

torch.manual_seed(0)
dt = torch.complex128
n = 4
perms = list(itertools.permutations(range(n)))                   # S_4, |G| = 24

def perm_matrix(p):
    """pi(sigma) e_x = e_{sigma(x)} on C^n."""
    P = torch.zeros(n, n, dtype=dt)
    P[list(p), list(range(n))] = 1
    return P

G_node = torch.stack([perm_matrix(p) for p in perms])            # S_4 on C^4: vertex features
G_edge = torch.stack([torch.kron(P, P) for P in G_node])         # S_4 on C^16: X -> P X P^T, row-major

def n_orbits(tuples):
    """Number of orbits of S_4 on tuples of vertices: the canonical form is the smallest image."""
    return len({min(tuple(p[i] for i in t) for p in perms) for t in tuples})

bell = {r: n_orbits(itertools.product(range(n), repeat=r)) for r in (2, 3, 4)}
print(bell)
assert bell == {2: 2, 3: 5, 4: 15}                               # Corollary 12.2(b)
for G, r in ((G_node, 2), (G_edge, 4)):                          # Lemma 12.1: dim pi(G)' = number of orbits
    assert dc.commutant_basis(G, group=True).shape[0] == dc.commutant_dimension(G) == bell[r]
for X in dc.commutant_basis(G_node, group=True):                 # Corollary 12.2(a): X = a 1 + b (J - 1)
    off = X[~torch.eye(n, dtype=torch.bool)]
    assert (X.diagonal() - X[0, 0]).abs().max() < 1e-12 and (off - off[0]).abs().max() < 1e-12

D_node, D_edge = dc.from_group(G_node), dc.from_group(G_edge)
print(D_node)
print(D_edge)
assert (D_node.k_factors, D_node.multiplicities) == ((1, 3), (1, 1))
assert (D_edge.k_factors, D_edge.multiplicities) == ((1, 2, 3, 3), (2, 1, 1, 3))
```

Así, $\mathbb C^{16}$ contiene la representación trivial dos veces, una de dimensión $2$ una vez, y dos de
dimensión $3$, una vez y tres veces: $16=1\cdot2+2\cdot1+3\cdot1+3\cdot3$ y $15=2^2+1^2+1^2+3^2$. Para dos
representaciones, todo se reduce a un solo conmutante:

**Lema 12.3 (suma directa).** Sea $\pi=\pi_1\oplus\pi_2$ sobre $\mathbb C^{n_1}\oplus\mathbb C^{n_2}$, con las
inclusiones isométricas $Q_i$ y las proyecciones $P_i=Q_iQ_i^*$. Entonces $P_1,P_2\in\pi(G)'$,
$\operatorname{Hom}_G(\pi_1,\pi_2)=Q_2^*\,\pi(G)'\,Q_1$, y $T\mapsto Q_2TQ_1^*$ es una biyección lineal de
$\operatorname{Hom}_G(\pi_1,\pi_2)$ sobre $P_2\,\pi(G)'P_1$.

*Demostración.* $\pi(g)Q_i=Q_i\pi_i(g)$ y $Q_i^*\pi(g)=\pi_i(g)Q_i^*$, porque $\pi(g)$ es diagonal por bloques; por
tanto, $P_i\pi(g)=Q_i\pi_i(g)Q_i^*=\pi(g)P_i$. Si $X\in\pi(G)'$, entonces
$Q_2^*XQ_1\pi_1(g)=Q_2^*X\pi(g)Q_1=Q_2^*\pi(g)XQ_1=\pi_2(g)Q_2^*XQ_1$. Si $T$ es equivariante, entonces
$\pi(g)Q_2TQ_1^*=Q_2\pi_2(g)TQ_1^*=Q_2T\pi_1(g)Q_1^*=Q_2TQ_1^*\pi(g)$, luego $Q_2TQ_1^*\in\pi(G)'$, es igual a
$P_2(Q_2TQ_1^*)P_1$, y $Q_2^*(Q_2TQ_1^*)Q_1=T$. Recíprocamente, todo $Y=P_2XP_1$ con $X\in\pi(G)'$ está en
$\pi(G)'$ y $Y=Q_2(Q_2^*YQ_1)Q_1^*$. $\square$

**Teorema 12.4 (las aplicaciones equivariantes son matrices por bloques).** Sean $\rho_1,\dots,\rho_C$ las
representaciones irreducibles, dos a dos no equivalentes, que aparecen en $\pi_1\oplus\pi_2$, $k_c=\dim\rho_c$,
$a_c=m_{\rho_c}(\pi_1)$, $b_c=m_{\rho_c}(\pi_2)$. Existen unitarias
$V_1:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{a_c}\to\mathbb C^{n_1}$ y
$V_2:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{b_c}\to\mathbb C^{n_2}$ con
$V_i^*\pi_i(g)V_i=\bigoplus_c\rho_c(g)\otimes1$, tales que
$$
\operatorname{Hom}_G(\pi_1,\pi_2)=\Big\{V_2\Big(\bigoplus_{c=1}^C1_{k_c}\otimes y_c\Big)V_1^*:\ y_c\in M_{b_c\times a_c}(\mathbb C)\Big\},
$$
y $(y_c)_c\mapsto T$ es inyectiva. Por tanto, $\operatorname{Hom}_G(\pi_1,\pi_2)\cong\bigoplus_cM_{b_c\times a_c}(\mathbb C)$,
$$
\dim\operatorname{Hom}_G(\pi_1,\pi_2)=\sum_\rho m_\rho(\pi_1)\,m_\rho(\pi_2),
$$
donde la suma recorre todas las representaciones irreducibles salvo equivalencia; para $\pi_1=\pi_2$ esto es
$\pi(G)'\cong\bigoplus_cM_{m_c}$ (proposición 2.7), y $\operatorname{Hom}_G(\pi_1,\pi_2)=0$ si y solo si ninguna
representación irreducible aparece en ambas.

*Demostración.* Apliquemos el teorema 2.6 y la proposición 2.7 a $\pi=\pi_1\oplus\pi_2$: existe una unitaria $V$ con
$V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$ y $V^*\pi(G)'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$. Como
$X\mapsto V^*XV$ es un $*$-isomorfismo y $P_1\in\pi(G)'$ (lema 12.3), $V^*P_1V=\bigoplus_c1\otimes p_c$ con
proyecciones $p_c$. Elijamos unitarias $u_c$ con $u_c^*p_cu_c=1_{a_c}\oplus0_{b_c}$ (teorema 1.3), donde ahora
$a_c=\operatorname{rank}p_c$, $b_c=m_c-a_c$, y sustituyamos $V$ por $V\bigoplus_c(1\otimes u_c)$; esto no cambia
ninguna de las dos formas normales, pues $1\otimes u_c$ conmuta con $\rho_c(g)\otimes1$. Sean
$\iota_c:\mathbb C^{a_c}\to\mathbb C^{m_c}$ y $\kappa_c:\mathbb C^{b_c}\to\mathbb C^{m_c}$ las inclusiones de las
primeras $a_c$ y de las últimas $b_c$ coordenadas, y
$$
V_1=Q_1^*V\bigoplus_c(1\otimes\iota_c),\qquad V_2=Q_2^*V\bigoplus_c(1\otimes\kappa_c).
$$
$W_1=V\bigoplus_c(1\otimes\iota_c)$ es una isometría con imagen $\operatorname{range}P_1=\operatorname{range}Q_1$, luego
$V_1^*V_1=W_1^*P_1W_1=1$, y $V_1$ es unitaria porque $\sum_ck_ca_c=\operatorname{rank}P_1=n_1$. Además,
$\pi_1(g)V_1=Q_1^*\pi(g)W_1=Q_1^*V\big(\bigoplus_c\rho_c(g)\otimes1_{m_c}\big)\bigoplus_c(1\otimes\iota_c)
=V_1\bigoplus_c\rho_c(g)\otimes1_{a_c}$. Así, $\pi_1\cong\bigoplus_c\rho_c\otimes1_{a_c}$ con $\rho_c$ irreducibles
no equivalentes, y $a_c=m_{\rho_c}(\pi_1)$ por la unicidad de las multiplicidades (sección 2.4); análogamente para
$V_2$ y $b_c$. Por último, para $X=V\big(\bigoplus_c1\otimes x_c\big)V^*\in\pi(G)'$,
$P_2XP_1=V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*$, donde $y_c=\kappa_c^*x_c\iota_c$ es la esquina
inferior izquierda de tamaño $b_c\times a_c$ de $x_c$, una matriz arbitraria. Por el lema 12.3,
$\operatorname{Hom}_G(\pi_1,\pi_2)=Q_2^*P_2\,\pi(G)'P_1Q_1$ consta de las matrices
$Q_2^*V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*Q_1=V_2\big(\bigoplus_c1\otimes y_c\big)V_1^*$, y
$(y_c)\mapsto T$ es inyectiva porque $V_1,V_2$ son unitarias. Una $\rho$ irreducible que no aparece en
$\pi_1\oplus\pi_2$ tiene $m_\rho(\pi_1)m_\rho(\pi_2)=0$. $\square$

Para el corolario 12.2(b) con $n=4$: las características de vértices $\pi_1$ tienen $a=(1,0,0,1)$ en los cuatro
sectores de $\pi_1\oplus\pi_2$ (la trivial, la de dimensión $2$ y las dos de dimensión $3$), y las características de
aristas $\pi_2$ tienen $b=(2,1,1,3)$, luego $\dim\operatorname{Hom}_G=1\cdot2+1\cdot3=5=B_3$.

> **Significado físico.** El teorema 12.4 es el núcleo de álgebra lineal de las **reglas de selección** y el caso
> escalar del teorema de Wigner–Eckart: un operador que conmuta con la simetría solo puede conectar estados del mismo
> tipo de simetría, y dentro de un tipo actúa como $1_{k_c}\otimes y_c$, es decir, mediante un «elemento de matriz
> reducido» $y_c$ que es ciego a la posición dentro del multiplete. (Los operadores que se transforman ellos mismos
> según una representación no trivial, como los vectoriales, necesitan además los coeficientes de Clebsch–Gordan.)
> Un hamiltoniano simétrico es el caso $\pi_1=\pi_2$; sus degeneraciones $k_c$ son las dimensiones de las
> representaciones irreducibles.

## 12.2 Una capa equivariante con la biblioteca

La demostración del teorema 12.4 es un algoritmo, una vez conocida la descomposición del teorema 2.6; `dc.from_group`
la calcula mediante el método aleatorizado de la sección 4.3. Para $\pi_1=\pi_2$ es simplemente el teorema 2.6: los
pesos son bloques $w_c\in M_{m_c}$, un operador del conmutante abstracto `D.commutant_algebra()` $=\bigoplus_cM_{m_c}$,
y la capa es $x\mapsto V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*x$, calculada por `D.from_commutant`. El módulo
siguiente es código del capítulo, no forma parte de la biblioteca.

```python
class EquivariantLinear(torch.nn.Module):
    """x -> T x, T = V ((+)_c 1_{k_c} (x) w_c) V^* in pi(G)'; one weight block w_c of size m_c x m_c per sector."""
    def __init__(self, D):
        super().__init__()
        self.D, self.CA = D, D.commutant_algebra()
        self.w = torch.nn.ParameterList(torch.nn.Parameter(torch.randn(m, m, dtype=dt) / m ** 0.5)
                                        for m in D.multiplicities)

    def matrix(self):
        return self.D.from_commutant(self.CA.from_blocks(list(self.w)))[0]

    def forward(self, x):                                        # x: (batch, n)
        return x @ self.matrix().T

def equivariance_error(layer, G_in, G_out):
    """max_g |layer(pi_in(g) x) - pi_out(g) layer(x)| on a random batch."""
    x = torch.randn(8, G_in.shape[-1], dtype=dt)
    with torch.no_grad():
        return max((layer(x @ g.T) - layer(x) @ h.T).abs().max().item() for g, h in zip(G_in, G_out))

edge_layer = EquivariantLinear(D_edge)
n_params = sum(p.numel() for p in edge_layer.parameters())
err = equivariance_error(edge_layer, G_edge, G_edge)
print(f"edge -> edge: {n_params} parameters, equivariance error {err:.1e}")
assert n_params == D_edge.commutant_dim == bell[4] and err < 1e-12
```

Para $\pi_1\ne\pi_2$ se sigue la demostración línea a línea: se descompone $\pi_1\oplus\pi_2$, se calculan los
bloques $p_c$ de $P_1$ con `to_commutant`, se diagonalizan y se conservan las esquinas inferiores izquierdas. La
traza de $p_c$ es $a_c=m_c(\pi_1)$.

```python
class EquivariantMap(torch.nn.Module):
    """x -> T x with T in Hom_G(pi_1, pi_2) (Theorem 12.4): one weight block y_c of size b_c x a_c per sector."""
    def __init__(self, G1, G2):
        super().__init__()
        self.n1 = n1 = G1.shape[-1]
        self.D = dc.from_group(torch.stack([torch.block_diag(g1, g2) for g1, g2 in zip(G1, G2)]))
        P1 = torch.block_diag(torch.eye(n1, dtype=dt), torch.zeros(G2.shape[-1], G2.shape[-1], dtype=dt))
        p = self.D.to_commutant(P1).matrix[0]                    # blocks p_c of P_1, padded
        self.u, self.a = [], []
        for c, m in enumerate(self.D.multiplicities):
            lam, u = torch.linalg.eigh(p[c, :m, :m])             # eigenvalues 0, ..., 0, 1, ..., 1
            self.a.append(round(lam.sum().item()))               # a_c = rank p_c = m_c(pi_1)
            self.u.append(u.flip(-1))                            # u_c^* p_c u_c = 1_{a_c} (+) 0_{b_c}
        self.y = torch.nn.ParameterList(torch.nn.Parameter(torch.randn(m - a, a, dtype=dt))
                                        for m, a in zip(self.D.multiplicities, self.a))

    def matrix(self):
        blocks = []
        for u, a, y in zip(self.u, self.a, self.y):
            x = torch.zeros(u.shape, dtype=dt)
            x[a:, :a] = y                                        # kappa_c y_c iota_c^*
            blocks.append(u @ x @ u.mH)
        T = self.D.from_commutant(self.D.commutant_algebra().from_blocks(blocks))[0]
        return T[self.n1:, :self.n1]                             # Q_2^* (...) Q_1

    def forward(self, x):
        return x @ self.matrix().T

node_to_edge = EquivariantMap(G_node, G_edge)
n_params = sum(p.numel() for p in node_to_edge.parameters())
err = equivariance_error(node_to_edge, G_node, G_edge)
print(f"node -> edge: sectors k = {node_to_edge.D.k_factors}, a = {node_to_edge.a}, "
      f"{n_params} parameters, equivariance error {err:.1e}")
assert n_params == bell[3] and err < 1e-12
corners = node_to_edge.D.basis('commutant')[:, n:, :n].reshape(-1, 16 * n)   # Lemma 12.3: Q_2^* pi(G)' Q_1
assert torch.linalg.matrix_rank(corners).item() == bell[3]
```

La segunda manera de hacer equivariante una capa es promediarla sobre el grupo. El **twirl** es
$$
\mathcal P(W)=\frac1{|G|}\sum_{g\in G}\pi_2(g)\,W\,\pi_1(g)^*,\qquad W\in M_{n_2\times n_1}(\mathbb C).
$$

**Teorema 12.5 (el twirl es la proyección ortogonal).** $\mathcal P$ es la proyección ortogonal, para el producto
escalar de Hilbert–Schmidt, de $M_{n_2\times n_1}(\mathbb C)$ sobre $\operatorname{Hom}_G(\pi_1,\pi_2)$, y
$$
\dim\operatorname{Hom}_G(\pi_1,\pi_2)=\frac1{|G|}\sum_{g\in G}\overline{\chi_{\pi_1}(g)}\,\chi_{\pi_2}(g).
$$

*Demostración.* Para $h\in G$, $\pi_2(h)\mathcal P(W)\pi_1(h)^*=\frac1{|G|}\sum_g\pi_2(hg)W\pi_1(hg)^*=\mathcal P(W)$,
pues $g\mapsto hg$ es una biyección de $G$; así, $\mathcal P(W)$ es equivariante. Si $W$ es equivariante, cada
término es igual a $W$, luego $\mathcal P(W)=W$. Por tanto, $\mathcal P^2=\mathcal P$ con imagen
$\operatorname{Hom}_G(\pi_1,\pi_2)$. Es autoadjunta:
$\langle S,\mathcal P(W)\rangle=\frac1{|G|}\sum_g\operatorname{Tr}\big((\pi_2(g)^*S\pi_1(g))^*W\big)$
$=\langle\mathcal P(S),W\rangle$, donde hemos sustituido $g\mapsto g^{-1}$ y usado $\pi_i(g^{-1})=\pi_i(g)^*$. Un
idempotente autoadjunto es la proyección ortogonal sobre su imagen. Su traza es su rango. La aplicación lineal
$W\mapsto AWB$ sobre $M_{n_2\times n_1}$ tiene traza $\operatorname{Tr}A\operatorname{Tr}B$ (el coeficiente de la
unidad matricial $E_{ij}$ en $AE_{ij}B$ es $A_{ii}B_{jj}$); con $A=\pi_2(g)$, $B=\pi_1(g)^*$ esto da la
fórmula. $\square$

Junto con el teorema 12.4,
$\sum_\rho m_\rho(\pi_1)m_\rho(\pi_2)=|G|^{-1}\sum_g\overline{\chi_{\pi_1}(g)}\chi_{\pi_2}(g)$; para $\pi_1=\pi_2$ esta
es la fórmula de caracteres para $\dim\pi(G)'$ de la sección 4.2, que evalúa `dc.commutant_dimension`. La
parametrización $W\mapsto\mathcal P(W)$ alcanza toda aplicación equivariante, $\mathcal P(W)$ es la aplicación
equivariante más próxima a $W$, y, como $\mathcal P$ es autoadjunta, el gradiente respecto de $W$ es $\mathcal P$
aplicada al gradiente respecto de $T=\mathcal P(W)$. Pero se desperdician $n_1n_2-\dim\operatorname{Hom}_G$
parámetros en el núcleo ($256-15$ para las características de aristas), y una evaluación cuesta
$O(|G|n_1n_2(n_1+n_2))$ para $\pi(g)$ densas, con $|G|=n!$ para $S_n$. La biblioteca calcula la misma proyección a
partir de la descomposición, `D.conditional_expectation(W, 'commutant')`, con cuatro productos $n\times n$.

```python
W = torch.randn(16, 16, dtype=dt)
PW = dc.twirl(G_edge, W)
S = torch.randn(16, 16, dtype=dt)
hs = lambda A, B: (A.conj() * B).sum()                           # Hilbert-Schmidt inner product
T_layer = edge_layer.matrix().detach()
assert (PW - D_edge.conditional_expectation(W, 'commutant')[0]).abs().max() < 1e-12
assert (dc.twirl(G_edge, T_layer) - T_layer).abs().max() < 1e-12           # P = id on pi(G)'
assert abs(hs(S, PW) - hs(dc.twirl(G_edge, S), W)) < 1e-10                  # P is self-adjoint
assert abs(hs(T_layer, W - PW)) < 1e-10                                     # W - P(W) is orthogonal to pi(G)'
chi_node, chi_edge = (G.diagonal(dim1=-2, dim2=-1).sum(-1) for G in (G_node, G_edge))
assert abs((chi_node.conj() * chi_edge).mean() - bell[3]) < 1e-12           # Theorem 12.5: 5
```

Las capas anteriores son lineales. Una red las alterna con no linealidades, y estas también han de ser
equivariantes. Para una representación de permutación, una función aplicada a cada coordenada es equivariante, ya
que $(\pi(g)x)_{gy}=x_y$; para otras representaciones no lo es, y se usan no linealidades que actúan sobre normas o
sobre los coeficientes en una base que $G$ permuta (ejercicio 3).

> **Significado físico.** El twirl es la esperanza condicional sobre el conmutante (capítulo 2, ejercicio 3):
> promediar un acoplamiento sobre todas las operaciones de simetría conserva exactamente su parte simétrica.
> «Simetrizar» un hamiltoniano o una interacción sumando sobre el grupo de simetría es la misma proyección.

## 12.3 Capas restringidas y espectrales

Los pesos $w=(w_c)_c$ de `EquivariantLinear` viven en una $*$-álgebra, de modo que toda construcción de los
capítulos 1 y 11 se les aplica. Escribamos $\Phi(w)=V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*$ para
$w\in\bigoplus_cM_{m_c}$.

**Proposición 12.6.** (a) $\Phi$ es un $*$-isomorfismo de $\bigoplus_cM_{m_c}$ sobre $\pi(G)'$, y $\Phi(w)$ es
autoadjunto, unitario o positivo si y solo si lo es cada $w_c$. (b) Para $w$ normal y una función $f$ sobre
$\bigcup_c\operatorname{spec}w_c$, $f(\Phi(w))=\Phi\big((f(w_c))_c\big)$; en particular, $f(X)\in\pi(G)'$ para todo
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

Así, una capa equivariante ortogonal (unitaria) es `optimize.UnitaryParam` sobre `D.commutant_algebra()` seguida de
$\Phi$; por (c), esta parametrización alcanza toda unitaria equivariante, y cada iterado de un método de gradiente es
exactamente unitario y exactamente equivariante. Del mismo modo, `PositiveParam` da capas equivariantes positivas,
normalizadas como en el capítulo 11 (espectro en $[\ell,1]$ con máximo $1$; por (a), una capa positiva general es
$\Phi(t^*t)$), y una capa espectral $f(L)$ de un generador equivariante autoadjunto $L$, como la difusión $e^{-tL}$,
requiere descomposiciones espectrales de coste $\sum_cm_c^3$ (aquí $37$) en lugar de $n^3=4096$, con gradientes
mediante la fórmula de Daleckii–Krein (`apply_function`). El código ajusta una unitaria equivariante a $32$ pares
entrada–salida, y contrasta (b) y su gradiente respecto de $t$, $\frac d{dt}e^{-tL}=-Le^{-tL}$, con matrices densas.

```python
from torch_vn_algebra import optimize as opt

CA = D_edge.commutant_algebra()                                  # M_2 (+) M_1 (+) M_1 (+) M_3
U_target = D_edge.from_commutant(CA.random_unitary_operator(1))[0]
X = torch.randn(32, 16, dtype=dt)
Y = X @ U_target.T

U = opt.UnitaryParam(CA, 1)                                      # u_0 Haar-random, s = 0
adam = torch.optim.Adam(U.parameters(), lr=0.05)
for step in range(150):
    U_full = D_edge.from_commutant(U.operator())[0]
    loss = (X @ U_full.T - Y).abs().pow(2).mean()
    if step == 0:
        loss0 = loss.item()
    adam.zero_grad()
    loss.backward()
    adam.step()
U_full = D_edge.from_commutant(U.operator())[0].detach()
unit_err = (U_full.mH @ U_full - torch.eye(16, dtype=dt)).abs().max().item()
equi_err = max((g @ U_full - U_full @ g).abs().max().item() for g in G_edge)
print(f"loss {loss0:.3f} -> {loss.item():.1e};  |U*U - 1| = {unit_err:.1e}, |[U, pi(g)]| = {equi_err:.1e}")
assert loss.item() < 1e-3 * loss0 and unit_err < 1e-12 and equi_err < 1e-12

h = torch.randn(1, CA.C, CA.k_max, CA.k_max, dtype=dt)
h = CA.operator(h + h.mH, is_self_adjoint=True)                  # a self-adjoint generator h = (h_c)
L = D_edge.from_commutant(h)[0]
t = torch.tensor(0.7, dtype=torch.float64, requires_grad=True)
E = D_edge.from_commutant(h.apply_function(lambda w: torch.exp(-t * w)))[0]   # e^{-tL}, in the small algebra
Z = torch.randn(16, 16, dtype=dt)
(Z.conj() * E).sum().real.backward()                             # d/dt Re <Z, e^{-tL}>
dE = -L @ E.detach()
assert (E - torch.linalg.matrix_exp(-0.7 * L)).abs().max() < 1e-12          # Proposition 12.6(b)
assert abs(t.grad.item() - (Z.conj() * dE).sum().real.item()) < 1e-10
```

> **Significado físico.** Un hamiltoniano simétrico $H\in\pi(G)'$ genera una evolución temporal simétrica
> $e^{-itH}\in\pi(G)'$ (capítulo 8), que nunca mezcla sectores y actúa solo sobre los espacios de multiplicidad. Una
> capa unitaria equivariante es una evolución de este tipo, aprendida en lugar de dada.

## 12.4 Invariantes

Un vector $x$ es **invariante** si $\pi(g)x=x$ para todo $g$; una función $q$ sobre $\mathbb C^n$ es invariante si
$q(\pi(g)x)=q(x)$. Las características invariantes son lo que puede devolver una red invariante (por ejemplo, un
clasificador de grafos).

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

*Demostración.* (a) Los vectores invariantes son $\operatorname{Hom}_G(\mathrm{triv},\pi)$, con $\mathrm{triv}$ la
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
invariantes, y la dimensión de estos es $\operatorname{Tr}M=|G|^{-1}\sum_g\operatorname{Tr}\sigma(g)$. $A$ es unitaria,
luego $A=R\Lambda R^{-1}$ con $\Lambda=\operatorname{diag}(\lambda)$ (teorema 1.3). Las matrices simétricas
$B_{ij}=(R^{-1})^T(e_ie_j^T+e_je_i^T)R^{-1}$, $i\le j$, forman una base, y
$A^TB_{ij}A=(R^{-1})^T\Lambda(e_ie_j^T+e_je_i^T)\Lambda R^{-1}=\lambda_i\lambda_jB_{ij}$. Así,
$\operatorname{Tr}\sigma(g)=\sum_{i\le j}\lambda_i\lambda_j=\frac12\big((\operatorname{Tr}A)^2+\operatorname{Tr}A^2\big)
=\frac12\big(\chi_\pi(g^{-1})^2+\chi_\pi(g^{-2})\big)$, y sumar sobre $g^{-1}$ en lugar de $g$ da la fórmula. Para
una representación de permutación, $G$ permuta los monomios $x_ax_b$, de modo que un polinomio es invariante si y
solo si sus coeficientes son constantes en las órbitas de monomios. $\square$

Para grafos, las características de aristas invariantes están generadas por las indicadoras de la diagonal y de su
complemento ($B_2=2$), y las lecturas lineales invariantes son la suma de las entradas diagonales y la de las no
diagonales. El código encuentra el sector trivial a partir de los bloques $\rho_c(g)$ (`to_algebra`) y comprueba
(a), (b) y los dos recuentos de (c).

```python
rho = D_edge.to_algebra(G_edge).matrix                           # blocks rho_c(g), shape (24, C, 3, 3)
c0 = [c for c, k in enumerate(D_edge.k_factors) if k == 1 and (rho[:, c, 0, 0] - 1).abs().max() < 1e-12]
P0 = G_edge.mean(0)
assert c0 == [0] and (P0 - D_edge.central_projections()[0]).abs().max() < 1e-12      # (a)
assert abs(chi_edge.mean() - D_edge.multiplicities[0]) < 1e-12
diag = torch.eye(n, dtype=dt).reshape(16)
assert (P0 @ diag - diag).abs().max() < 1e-12 and (P0 @ (1 - diag) - (1 - diag)).abs().max() < 1e-12

w = [torch.randn(m, m, dtype=dt) for m in D_edge.multiplicities]
A = D_edge.from_commutant(CA.from_blocks([b + b.mH for b in w]))[0]                   # (b)
x = torch.randn(16, dtype=dt)
q = lambda v: (v.conj() @ A @ v).real
assert (A - A.mH).abs().max() < 1e-12
assert max(abs(q(g @ x) - q(x)).item() for g in G_edge) < 1e-12

def degree2(G):                                                  # (c): characters
    chi, chi_sq = (H.diagonal(dim1=-2, dim2=-1).sum(-1) for H in (G, G @ G))
    return ((chi ** 2 + chi_sq) / 2).mean().real.item()

pairs = list(itertools.product(range(n), repeat=2))
multisets = [(a, b) for i, a in enumerate(pairs) for b in pairs[i:]]
orbits2 = len({min(tuple(sorted(tuple(p[v] for v in e) for e in (a, b))) for p in perms) for a, b in multisets})
print(f"invariant quadratic polynomials: vertices {degree2(G_node):.4f}, edges {degree2(G_edge):.4f} "
      f"(orbits of S_4 on monomials: {orbits2})")
assert abs(degree2(G_node) - 2) < 1e-12 and abs(degree2(G_edge) - orbits2) < 1e-12 and orbits2 == 11
```

Sobre $\mathbb C^4$ son $\sum_ix_i^2$ y $\sum_{i<j}x_ix_j$. Todos los grados se cuentan mediante la fórmula de Molien
$\sum_d\dim\mathbb C[x]^G_dt^d=|G|^{-1}\sum_g\det(1-t\pi(g))^{-1}$, cuyo coeficiente de $t^2$ es (c) (la enunciamos
sin demostración; B. Sturmfels, *Algorithms in Invariant Theory*, Springer 2008, teorema 2.2.1).

> **Significado físico.** En la teoría de Landau de las transiciones de fase, la energía libre es un polinomio
> invariante en un parámetro de orden real $x$, sobre el que $G$ actúa mediante matrices ortogonales reales. Su parte
> cuadrática es $x^TBx$ con $B$ en el conmutante (pues $\pi(g)^T=\pi(g)^{-1}$), de modo que, por la
> proposición 12.6(a), es diagonal por bloques sobre los sectores: la fase simétrica se vuelve inestable en el sector
> en el que un valor propio de $B$ cambia de signo en primer lugar, y (genéricamente) la representación irreducible
> de ese sector determina la simetría de la fase ordenada.

## 12.5 Con qué comparar

Hay bibliotecas maduras que proporcionan capas equivariantes para grupos concretos. **escnn** (Cesa, Lang y Weiler
2022, sucesora de e2cnn) construye convoluciones orientables (*steerable*) para subgrupos de $O(2)$ y $O(3)$ a partir
de representaciones irreducibles y de bases armónicas de núcleos; **e3nn** (Geiger y Smidt 2022) construye redes
$E(3)$-equivariantes a partir de armónicos esféricos y de productos tensoriales de Clebsch–Gordan. **EMLP** (Finzi,
Welling y Wilson 2021) trata cualquier grupo de matrices: impone $\pi_2(g)T=T\pi_1(g)$ para generadores discretos y
$d\pi_2(X)T=Td\pi_1(X)$ para generadores del álgebra de Lie, y calcula una base de las soluciones como un núcleo
numérico, que es el problema lineal de `dc.commutant_basis`. Las convoluciones de grupo (Cohen y Welling 2016) son el
caso de la representación regular, y Kondor y Trivedi (2018) muestran que, para mapas de características sobre
espacios homogéneos de un grupo compacto, las capas lineales equivariantes son convoluciones generalizadas. Para
permutaciones, DeepSets (Zaheer et al. 2017) y las bases de Maron et al. (2019) son el corolario 12.2.

La biblioteca no proporciona convoluciones sobre rejillas, núcleos orientables sobre espacios continuos,
descomposiciones sobre $\mathbb R$ (tipos real, complejo y cuaterniónico; el ejercicio 2 trata las matrices reales)
ni $n$ grande: $V$ es una matriz densa $n\times n$ y la descomposición cuesta $O(|G|n^3)$ por muestra aleatoria para
un grupo dado por todos sus elementos (los generadores de álgebras de Lie pasan por `dc.from_generators`;
sección 4.3). Lo que aporta es la **estructura** del espacio de soluciones: una base del conmutante es una lista de
matrices, mientras que la descomposición identifica el conmutante con la $*$-álgebra $\bigoplus_cM_{m_c}$. Los
productos, los adjuntos, la unitariedad, la positividad y el cálculo funcional actúan entonces bloque a bloque sobre
matrices pequeñas, con gradientes (proposición 12.6), y `to_commutant`/`from_commutant` traducen entre las dos
imágenes.

## Ejercicios

1. (**DeepSets con canales.**) $S_n$ ($n\ge2$) actúa sobre $\mathbb C^n\otimes\mathbb C^d$ mediante
   $P_\sigma\otimes1_d$. Demuestre que las aplicaciones equivariantes de $\mathbb C^n\otimes\mathbb C^{d_1}$ en
   $\mathbb C^n\otimes\mathbb C^{d_2}$ son las $1\otimes A+J\otimes B$ con $A,B\in M_{d_2\times d_1}$, y compruébelo
   con `EquivariantMap`.
2. (**Capas reales.**) Supongamos que $\pi_1,\pi_2$ tienen matrices reales (por ejemplo, representaciones de
   permutación). Demuestre que las matrices equivariantes reales forman un espacio vectorial real de dimensión
   $\dim_{\mathbb C}\operatorname{Hom}_G(\pi_1,\pi_2)$, y que $T\mapsto\operatorname{Re}T$ aplica
   $\operatorname{Hom}_G(\pi_1,\pi_2)$ sobre dicho espacio.
3. (**No linealidades.**) (a) Demuestre que aplicar una función $f:\mathbb C\to\mathbb C$ a cada coordenada conmuta
   con toda representación de permutación. (b) Demuestre que no conmuta con la rotación de ángulo $\pi/2$ de
   $\mathbb R^2$ para $f(t)=\max(t,0)$. (c) Demuestre que $x\mapsto\varphi(\|x\|)x$ conmuta con toda representación
   unitaria.
4. (**¿Cuántas muestras determinan una capa?**) Sean $T=\Phi(w)\in\pi(G)'$ y $x_1,\dots,x_N\in\mathbb C^n$. Escriba
   $V^*x_i$ en el sector $c$ como $\sum_{a,j}(X_{i,c})_{aj}e_a\otimes e_j$ con $X_{i,c}\in M_{k_c\times m_c}$, y sea
   $S_c\in M_{Nk_c\times m_c}$ la matriz formada apilando $X_{1,c},\dots,X_{N,c}$. Demuestre que las $T$
   equivariantes con $Tx_1=\dots=Tx_N=0$ forman un espacio de dimensión $\sum_cm_c(m_c-\operatorname{rank}S_c)$;
   para $x_i$ gaussianos esto es $\sum_cm_c\max(0,m_c-Nk_c)$ casi seguramente. Para la representación de aristas,
   dos muestras determinan una capa equivariante, mientras que una matriz $16\times16$ arbitraria necesita dieciséis.

<details>
<summary>Soluciones</summary>

1. $P_\sigma J=J=JP_\sigma$, luego estas aplicaciones son equivariantes. Son linealmente independientes: si
   $1\otimes A+J\otimes B=0$, los bloques fuera de la diagonal dan $B=0$, y entonces $A=0$. Por el corolario 12.2(a)
   y el teorema 2.6, $P_\sigma\otimes1_d\cong(\mathrm{triv}\otimes1_d)\oplus(\mathrm{std}\otimes1_d)$, de modo que
   ambas representaciones irreducibles tienen multiplicidad $d$, y el teorema 12.4 da
   $\dim\operatorname{Hom}_G=2d_1d_2$, la dimensión del espacio de las $1\otimes A+J\otimes B$.

```python
G_in = torch.stack([torch.kron(P, torch.eye(2, dtype=dt)) for P in G_node])     # d_1 = 2
G_out = torch.stack([torch.kron(P, torch.eye(3, dtype=dt)) for P in G_node])    # d_2 = 3
deepsets = EquivariantMap(G_in, G_out)
assert sum(p.numel() for p in deepsets.parameters()) == 2 * 2 * 3
assert equivariance_error(deepsets, G_in, G_out) < 1e-12
T = deepsets.matrix().detach().reshape(n, 3, n, 2)                # blocks T[i, :, j, :] = delta_ij A + B
A_, B_ = T[0, :, 0] - T[0, :, 1], T[0, :, 1]
assert all((T[i, :, j] - (A_ if i == j else 0) - B_).abs().max() < 1e-12 for i in range(n) for j in range(n))
```

2. Si $T$ es equivariante, también lo es $\bar T$, ya que $\overline{\pi_i(g)}=\pi_i(g)$. Por tanto,
   $\operatorname{Re}T=(T+\bar T)/2$ y $\operatorname{Im}T=(T-\bar T)/2i$ son reales y equivariantes, de modo que la
   aplicación $T\mapsto\operatorname{Re}T$ es sobreyectiva sobre las matrices equivariantes reales (es la identidad
   sobre ellas), y $\operatorname{Hom}_G(\pi_1,\pi_2)=H_{\mathbb R}\oplus iH_{\mathbb R}$, donde $H_{\mathbb R}$ denota
   las matrices equivariantes reales. Una base real de $H_{\mathbb R}$ es entonces una base compleja de
   $\operatorname{Hom}_G$: lo genera, y una relación compleja $\sum_j(\alpha_j+i\beta_j)T_j=0$ con $T_j$ reales se
   separa en parte real y parte imaginaria. En el lema 12.1, las indicadoras de las órbitas son una base real de este
   tipo. (Los bloques $y_c$ del teorema 12.4 son complejos incluso para datos reales: $V$ es compleja en general, y la
   parametrización mediante $\operatorname{Re}T$ usa el doble de parámetros reales de los necesarios.)

3. (a) $(\pi(g)f(x))_{gy}=f(x)_y=f(x_y)=f((\pi(g)x)_{gy})$. (b) Para $R(x_1,x_2)=(-x_2,x_1)$ y $x=(1,-1)$,
   $f(Rx)=f(1,1)=(1,1)$, pero $Rf(x)=R(1,0)=(0,1)$. (c) $\|\pi(g)x\|=\|x\|$, luego
   $\varphi(\|\pi(g)x\|)\pi(g)x=\pi(g)\varphi(\|x\|)x$.

```python
relu = lambda v: v.clamp(min=0)
R = torch.tensor([[0.0, -1.0], [1.0, 0.0]], dtype=torch.float64)
v = torch.tensor([1.0, -1.0], dtype=torch.float64)
assert (relu(R @ v) - R @ relu(v)).abs().max() == 1                             # (b)
xr = torch.randn(16, dtype=torch.float64)
assert all(torch.equal(relu(g.real @ xr), g.real @ relu(xr)) for g in G_edge)   # (a): exact
gate = lambda v: torch.tanh(v.norm()) * v
assert max((gate(g @ x) - g @ gate(x)).abs().max().item() for g in G_edge) < 1e-12   # (c)
```

4. $(1\otimes w_c)(e_a\otimes e_j)=\sum_le_a\otimes(w_c)_{lj}e_l$, luego $V^*Tx_i$ tiene, en el sector $c$, la matriz
   de coeficientes $X_{i,c}w_c^T$. Por tanto, $Tx_i=0$ para todo $i$ si y solo si $S_cw_c^T=0$ para todo $c$, es
   decir, si y solo si cada columna de $w_c^T$ está en $\ker S_c$, un espacio de dimensión
   $m_c-\operatorname{rank}S_c$: las $w_c$ con esta propiedad forman un espacio de dimensión
   $m_c(m_c-\operatorname{rank}S_c)$, y $\Phi$ es inyectiva. Para $x_i$ gaussianos, $V^*x_i$ es gaussiano ($V$ es
   unitaria), de modo que $S_c$ tiene entradas gaussianas independientes y rango $\min(Nk_c,m_c)$ casi seguramente:
   para la submatriz cuadrada superior izquierda de ese tamaño, $|\det|^2$ es un polinomio real en las partes real e
   imaginaria de las entradas, no idénticamente nulo (vale $1$ cuando la submatriz es la identidad), y un polinomio
   así solo se anula en un conjunto de medida de Lebesgue nula (demostrado en la sección 4.3). Para la
   representación de aristas, $(k,m)=(1,2),(2,1),(3,1),(3,3)$, solo contribuye el sector $(1,2)$: $2\cdot(2-1)=2$
   para $N=1$, y $0$ para $N=2$.

```python
basis = D_edge.basis('commutant')                                # (15, 16, 16), Hilbert-Schmidt orthonormal
for N, lost in ((1, 2), (2, 0)):
    Xs = torch.randn(16, N, dtype=dt)
    M = torch.stack([(B @ Xs).reshape(-1) for B in basis], dim=1)               # w -> (T x_1, ..., T x_N)
    assert torch.linalg.matrix_rank(M).item() == 15 - lost
```

</details>

**Cómo seguir.** Este es el último capítulo del curso. El [capítulo 4](04_symmetries.md) desarrolla la teoría de
representaciones usada aquí: el lema de Schur, los caracteres y el algoritmo aleatorizado en que se basa
`dc.from_group` (secciones 4.1–4.3). El capítulo 9 (en preparación) aplica los mismos sectores a las cadenas de
espines, y la medida de Haar sobre $\prod_cU(m_c)$ del [capítulo 5](05_random_matrices.md) proporciona la
inicialización aleatoria natural de una capa unitaria equivariante. Para las matemáticas, véanse J.-P. Serre,
*Linear Representations of Finite Groups* (Springer 1977) y W. Fulton y J. Harris, *Representation Theory*
(Springer 1991); para el aprendizaje equivariante, M. Bronstein, J. Bruna, T. Cohen y P. Veličković, *Geometric Deep
Learning* (2021), y M. Weiler et al., *Equivariant and Coordinate Independent Convolutional Networks* (2023); para el
algoritmo de descomposición, las referencias de `torch_vn_algebra/decompose.py`.
