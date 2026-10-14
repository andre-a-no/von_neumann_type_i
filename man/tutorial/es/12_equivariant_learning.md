[English](../12_equivariant_learning.md) · **Español** · [Français](../fr/12_equivariant_learning.md) · [Deutsch](../de/12_equivariant_learning.md) · [中文](../zh/12_equivariant_learning.md) · [日本語](../ja/12_equivariant_learning.md) · [Русский](../ru/12_equivariant_learning.md)

# Tutorial 12. Aprendizaje automático equivariante

Este tutorial acompaña al [capítulo 12](../../lecture_notes/es/12_equivariant_learning.md) de los apuntes de clase, que demuestra que las aplicaciones lineales
equivariantes entre dos representaciones de un grupo finito son matrices por bloques ([teorema 12.4](../../lecture_notes/es/12_equivariant_learning.md#r-12-4)), que
promediar sobre el grupo es la proyección ortogonal sobre ellas ([teorema 12.5](../../lecture_notes/es/12_equivariant_learning.md#r-12-5)), y cuenta invariantes.
Aquí construimos y entrenamos capas equivariantes con `torch_vn_algebra.decompose`, para $S_4$ actuando sobre las
características de vértices y de aristas de grafos, y comprobamos numéricamente estos resultados.

En el código, un lote de vectores es un tensor de forma `(batch, n)`, sobre el que una matriz $T$ actúa como `x @ T.T`.

<a id="s-12-1"></a>

## 12.1 Aplicaciones lineales equivariantes

*Apuntes:* [sección 12.1](../../lecture_notes/es/12_equivariant_learning.md#s-12-1): [lema 12.1](../../lecture_notes/es/12_equivariant_learning.md#r-12-1) (órbitas) y [corolario 12.2](../../lecture_notes/es/12_equivariant_learning.md#r-12-2) (DeepSets, números de Bell).

La capa DeepSets (Zaheer et al. 2017) es [apuntes de clase, corolario 12.2](../../lecture_notes/es/12_equivariant_learning.md#r-12-2)(a); las bases de Maron et al. (2019) para redes de grafos son (b).
El código cuenta las órbitas por fuerza bruta, las compara con el conmutante calculado numéricamente ([lema 12.1](../../lecture_notes/es/12_equivariant_learning.md#r-12-1))
y descompone ambas representaciones de $S_4$:

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

Así, en la descomposición calculada por `dc.from_group` (los asserts comprueban su resultado), $\mathbb C^{16}$
contiene una representación de dimensión $1$ dos veces —la trivial, cuya multiplicidad es el número $2$ de órbitas
sobre $[4]^2$ por la proposición 12.7(a) y el lema 12.1 (el código de la sección 12.4 identifica este sector)—, una
de dimensión $2$ una vez, y dos de dimensión $3$, una vez y tres veces:
$16=1\cdot2+2\cdot1+3\cdot1+3\cdot3$, y $15=2^2+1^2+1^2+3^2$ concuerda con el recuento de órbitas del lema 12.1 y
con $\dim\pi(G)'=\sum_cm_c^2$ (teorema 2.6). Para dos representaciones, todo se reduce a un solo conmutante:

Según los [apuntes de clase, lema 12.3](../../lecture_notes/es/12_equivariant_learning.md#r-12-3), $\mathrm{Hom}_G(\pi_1,\pi_2)=Q_2^*\,\pi(G)'\,Q_1$ para $\pi=\pi_1\oplus\pi_2$, y por el
[teorema 12.4](../../lecture_notes/es/12_equivariant_learning.md#r-12-4) es, en bases adecuadas, $\bigoplus_cM_{b_c\times a_c}(\mathbb C)$, de dimensión $\sum_\rho m_\rho(\pi_1)\,m_\rho(\pi_2)$.

<a id="s-12-2"></a>

## 12.2 Una capa equivariante con la biblioteca

*Apuntes:* [sección 12.2](../../lecture_notes/es/12_equivariant_learning.md#s-12-2): [teorema 12.5](../../lecture_notes/es/12_equivariant_learning.md#r-12-5) (el twirl); las capas siguen las demostraciones del [lema 12.3](../../lecture_notes/es/12_equivariant_learning.md#r-12-3) y del [teorema 12.4](../../lecture_notes/es/12_equivariant_learning.md#r-12-4).

La demostración del teorema 12.4 es un algoritmo, una vez conocida la descomposición del teorema 2.6; `dc.from_group`
la calcula mediante el método aleatorizado de la sección 4.3 (correcto con probabilidad uno en aritmética exacta;
aquí los asserts comprueban el resultado calculado). Para $\pi_1=\pi_2$ es simplemente el teorema 2.6: los
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

Para el corolario 12.2(b) con $n=4$, en la descomposición calculada (impresa y comprobada más abajo): las
características de vértices $\pi_1$ tienen $a=(1,0,0,1)$ en los cuatro sectores de $\pi_1\oplus\pi_2$ (de
dimensiones $1,2,3,3$; el primero es el trivial, pues por el corolario 12.2(a) el único constituyente de dimensión
uno de $\pi_1$ es el trivial), y las características de aristas $\pi_2$ tienen $b=(2,1,1,3)$, luego
$\dim\mathrm{Hom}_G=1\cdot2+1\cdot3=5=B_3$, como demuestra el corolario 12.2(b).

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
b = tuple(m - a for m, a in zip(node_to_edge.D.multiplicities, node_to_edge.a))
assert node_to_edge.D.k_factors == (1, 2, 3, 3) and node_to_edge.a == [1, 0, 0, 1] and b == (2, 1, 1, 3)
corners = node_to_edge.D.basis('commutant')[:, n:, :n].reshape(-1, 16 * n)   # Lemma 12.3: Q_2^* pi(G)' Q_1
assert torch.linalg.matrix_rank(corners).item() == bell[3]
```

La segunda manera de hacer equivariante una capa es promediarla sobre el grupo. El **twirl** es
$$
\mathcal P(W)=\frac1{|G|}\sum_{g\in G}\pi_2(g)\,W\,\pi_1(g)^*,\qquad W\in M_{n_2\times n_1}(\mathbb C).
$$

Según los [apuntes de clase, teorema 12.5](../../lecture_notes/es/12_equivariant_learning.md#r-12-5), $\mathcal P$ es la proyección ortogonal, para el producto escalar de Hilbert–Schmidt, sobre $\mathrm{Hom}_G(\pi_1,\pi_2)$, y
$\dim\mathrm{Hom}_G(\pi_1,\pi_2)=|G|^{-1}\sum_g\overline{\chi_{\pi_1}(g)}\chi_{\pi_2}(g)$; para $\pi_1=\pi_2$, `dc.commutant_dimension` evalúa esta
fórmula de caracteres. La parametrización $W\mapsto\mathcal P(W)$ desperdicia $n_1n_2-\dim\mathrm{Hom}_G$ parámetros en el núcleo ($256-15$ para las características de aristas), y una
evaluación cuesta $O(|G|n_1n_2(n_1+n_2))$ para $\pi(g)$ densas, con $|G|=n!$ para $S_n$. La biblioteca calcula la misma
proyección a partir de la descomposición, `D.conditional_expectation(W, 'commutant')`, con cuatro productos $n\times n$.
El código comprueba que `dc.twirl` y `D.conditional_expectation` coinciden, que $\mathcal P$ es la identidad sobre $\pi(G)'$
y autoadjunta, y la fórmula de dimensión de las características de vértices a las de aristas.

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

<a id="s-12-3"></a>

## 12.3 Capas restringidas y espectrales

*Apuntes:* [sección 12.3](../../lecture_notes/es/12_equivariant_learning.md#s-12-3): [proposición 12.6](../../lecture_notes/es/12_equivariant_learning.md#r-12-6).

Los pesos $w=(w_c)_c$ de `EquivariantLinear` viven en una $*$-álgebra, de modo que toda construcción de los
capítulos 1 y 11 se les aplica. Escribamos $\Phi(w)=V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*$ para
$w\in\bigoplus_cM_{m_c}$.

Según los [apuntes de clase, proposición 12.6](../../lecture_notes/es/12_equivariant_learning.md#r-12-6), (a) $\Phi$ es un $*$-isomorfismo de $\bigoplus_cM_{m_c}$ sobre $\pi(G)'$, y $\Phi(w)$ es
autoadjunto, unitario o positivo si y solo si lo es cada $w_c$; (b) $f(\Phi(w))=\Phi\big((f(w_c))_c\big)$ para $w$ normal; (c) las
unitarias equivariantes son las $\Phi(u_0e^{s-s^*})$.

Así, una capa equivariante ortogonal (unitaria) es `optimize.UnitaryParam` sobre `D.commutant_algebra()` seguida de
$\Phi$; por (c), esta parametrización alcanza toda unitaria equivariante, y en aritmética exacta cada iterado de un
método de gradiente es unitario y equivariante (el código comprueba ambas cosas hasta $10^{-12}$ tras el
entrenamiento). Del mismo modo, `PositiveParam` da capas equivariantes positivas,
normalizadas como en el capítulo 11 (espectro en $[\ell,1]$ con máximo $1$; por (a), una capa positiva general es
$\Phi(t^*t)$), y una capa espectral $f(L)$ de un generador equivariante autoadjunto $L$, como la difusión $e^{-tL}$,
requiere descomposiciones espectrales de los bloques $w_c$, de coste del orden de $\sum_cm_c^3$ (aquí $37$) en lugar
de $n^3=4096$; la biblioteca diferencia a través de `apply_function` (no demostramos nada sobre su paso hacia atrás;
el código comprueba una derivada más abajo). El código ajusta una unitaria equivariante a $32$ pares
entrada–salida con Adam —una observación numérica de esta ejecución, no un teorema: la pérdida disminuye en más de
un factor $10^3$, lo que comprueba el assert, y no se afirma ningún resultado de convergencia— y contrasta (b) y el
gradiente respecto de $t$ con matrices densas y con la fórmula $\frac d{dt}e^{-tL}=-Le^{-tL}$ (para
$L=U\mathrm{diag}(\lambda)U^*$, $e^{-tL}=U\mathrm{diag}(e^{-t\lambda_i})U^*$ por el teorema 1.5, y $U$ no
depende de $t$).

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

<a id="s-12-4"></a>

## 12.4 Invariantes

*Apuntes:* [sección 12.4](../../lecture_notes/es/12_equivariant_learning.md#s-12-4): [proposición 12.7](../../lecture_notes/es/12_equivariant_learning.md#r-12-7) (vectores invariantes y formas cuadráticas) y [proposición 12.8](../../lecture_notes/es/12_equivariant_learning.md#r-12-8) (fórmula de Molien).

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

Sobre $\mathbb C^4$ son $\sum_ix_i^2$ y $\sum_{i<j}x_ix_j$. La demostración de (c) funciona en todos los grados:

Según los [apuntes de clase, proposición 12.8](../../lecture_notes/es/12_equivariant_learning.md#r-12-8) (fórmula de Molien),
$\sum_{d\ge0}\dim\mathbb C[x]^G_d\,t^d=|G|^{-1}\sum_g\det\big(1-t\pi(g)\big)^{-1}$, y para una representación de permutación
$\dim\mathbb C[x]^G_d$ es el número de órbitas de $G$ sobre los monomios de grado $d$.

(Es el teorema de Molien de 1897; para profundizar, véase B. Sturmfels, *Algorithms in Invariant Theory*, Springer
2008, sección 2.2.) El código desarrolla el miembro derecho hasta $t^4$, multiplicando por una serie geométrica por
cada valor propio, y lo compara con las órbitas de monomios; para $S_4$ sobre $\mathbb C^4$, los recuentos
$1,1,2,3,5$ son las particiones de $d$ en a lo sumo cuatro partes (las órbitas de multiconjuntos están determinadas
por sus multiplicidades, ordenadas).

```python
def molien(G, d_max):
    """Coefficients of |G|^{-1} sum_g det(1 - t pi(g))^{-1} up to t^d_max: products of geometric series."""
    lam = torch.linalg.eigvals(G)                                # (|G|, n), |lambda| = 1
    c = torch.zeros(G.shape[0], d_max + 1, dtype=dt)
    c[:, 0] = 1
    for i in range(lam.shape[1]):                                # multiply by 1 / (1 - t lambda_i)
        for d in range(1, d_max + 1):
            c[:, d] = c[:, d] + lam[:, i] * c[:, d - 1]
    return c.mean(0)

def monomial_orbits(G, d):
    """Orbits of G on the monomials of degree d, i.e. on the multisets of d points of Omega."""
    act = [g.real.argmax(0).tolist() for g in G]                 # pi(g) e_x = e_{act[x]}
    mons = itertools.combinations_with_replacement(range(G.shape[-1]), d)
    return len({min(tuple(sorted(a[x] for x in m)) for a in act) for m in mons})

for G in (G_node, G_edge):
    series, counts = molien(G, 4), [monomial_orbits(G, d) for d in range(5)]
    print("Molien series:", [f"{v.real:.4f}" for v in series], " orbits of monomials:", counts)
    assert (series - torch.tensor(counts, dtype=dt)).abs().max() < 1e-9          # Proposition 12.8
assert [monomial_orbits(G_node, d) for d in range(5)] == [1, 1, 2, 3, 5] and monomial_orbits(G_edge, 2) == orbits2
```

<a id="s-12-5"></a>

## 12.5 Con qué comparar

Hay bibliotecas maduras que proporcionan capas equivariantes para grupos concretos. **escnn** (Cesa, Lang y Weiler
2022, sucesora de e2cnn) construye convoluciones orientables (*steerable*) para subgrupos de $O(2)$ y $O(3)$ a partir
de representaciones irreducibles y de bases armónicas de núcleos; **e3nn** (Geiger y Smidt 2022) construye redes
$E(3)$-equivariantes a partir de armónicos esféricos y de productos tensoriales de Clebsch–Gordan. **EMLP** (Finzi,
Welling y Wilson 2021) trata cualquier grupo de matrices: impone $\pi_2(g)T=T\pi_1(g)$ para generadores discretos y
$d\pi_2(X)T=Td\pi_1(X)$ para generadores del álgebra de Lie, y calcula una base de las soluciones como un núcleo
numérico, que es el problema lineal de `dc.commutant_basis`. Las convoluciones de grupo (Cohen y Welling 2016) son el
caso $\Omega_1=\Omega_2=G$ del lema 12.1, con $G$ actuando sobre sí mismo por multiplicación a la izquierda: $(y,x)$ e
$(y',x')$ están en la misma órbita si y solo si $x^{-1}y=x'^{-1}y'$ (en tal caso, $g=x'x^{-1}$ lleva un par al otro),
de modo que las aplicaciones equivariantes son $(Tv)(y)=\sum_xf(x^{-1}y)v(x)$, las convoluciones con funciones $f$
sobre $G$. Kondor y Trivedi (2018) estudian la cuestión análoga para espacios homogéneos de grupos compactos, que
este curso no trata. Para permutaciones, DeepSets (Zaheer et al. 2017) y las bases de Maron et al. (2019) son el corolario 12.2.

La biblioteca no proporciona convoluciones sobre rejillas, núcleos orientables sobre espacios continuos,
descomposiciones sobre $\mathbb R$ (tipos real, complejo y cuaterniónico; el ejercicio 2 trata las matrices reales)
ni $n$ grande: $V$ es una matriz densa $n\times n$ y la descomposición cuesta $O(|G|n^3)$ por muestra aleatoria para
un grupo dado por todos sus elementos (los generadores de álgebras de Lie pasan por `dc.from_generators`;
sección 4.3). Lo que aporta es la **estructura** del espacio de soluciones: una base del conmutante es una lista de
matrices, mientras que la descomposición identifica el conmutante con la $*$-álgebra $\bigoplus_cM_{m_c}$. Los
productos, los adjuntos, la unitariedad, la positividad y el cálculo funcional actúan entonces bloque a bloque sobre
matrices pequeñas (proposición 12.6), la biblioteca diferencia a través de ellos, y
`to_commutant`/`from_commutant` traducen entre las dos imágenes.

<a id="x-12"></a>

## Ejercicios

[Ejercicios](../../lecture_notes/es/12_equivariant_learning.md#x-12) 1, 3 y 4 de los apuntes de clase, numéricamente:

- *Ejercicio 1* (DeepSets con canales): compruebe con `EquivariantMap` que las aplicaciones equivariantes de
  $\mathbb C^n\otimes\mathbb C^{d_1}$ en $\mathbb C^n\otimes\mathbb C^{d_2}$ son las $1\otimes A+J\otimes B$.
- *Ejercicio 3* (no linealidades): ReLU sobre las coordenadas conmuta con la representación de permutación, pero no
  con la rotación de ángulo $\pi/2$ de $\mathbb R^2$; $x\mapsto\varphi(\|x\|)x$ conmuta con toda representación unitaria.
- *Ejercicio 4*: para la representación de aristas, una muestra gaussiana deja sin determinar un espacio de capas
  equivariantes de dimensión $2$; dos muestras determinan la capa.

<details>
<summary>Soluciones</summary>

*Ejercicio 1.*

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

*Ejercicio 3.*

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

*Ejercicio 4.*

```python
basis = D_edge.basis('commutant')                                # (15, 16, 16), Hilbert-Schmidt orthonormal
for N, lost in ((1, 2), (2, 0)):
    Xs = torch.randn(16, N, dtype=dt)
    M = torch.stack([(B @ Xs).reshape(-1) for B in basis], dim=1)               # w -> (T x_1, ..., T x_N)
    assert torch.linalg.matrix_rank(M).item() == 15 - lost
```

</details>

**Siguiente:** [Tutorial 13. Álgebras de Jordan: por qué la mecánica cuántica usa *-álgebras](13_jordan_algebras.md).
