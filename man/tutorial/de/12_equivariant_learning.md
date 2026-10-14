[English](../12_equivariant_learning.md) · [Español](../es/12_equivariant_learning.md) · [Français](../fr/12_equivariant_learning.md) · **Deutsch** · [中文](../zh/12_equivariant_learning.md) · [日本語](../ja/12_equivariant_learning.md) · [Русский](../ru/12_equivariant_learning.md)

# Tutorial 12. Äquivariantes maschinelles Lernen

Dieses Tutorial begleitet [Kapitel 12](../../lecture_notes/de/12_equivariant_learning.md) der Vorlesungsnotizen, das beweist, dass die äquivarianten linearen
Abbildungen zwischen zwei Darstellungen einer endlichen Gruppe Blockmatrizen sind ([Satz 12.4](../../lecture_notes/de/12_equivariant_learning.md#r-12-4)), dass
die Mittelung über die Gruppe die Orthogonalprojektion auf sie ist ([Satz 12.5](../../lecture_notes/de/12_equivariant_learning.md#r-12-5)), und das Invarianten abzählt.
Hier konstruieren und trainieren wir äquivariante Schichten mit `torch_vn_algebra.decompose`, für $S_4$, das auf den Knoten- und
Kantenmerkmalen von Graphen wirkt, und überprüfen diese Ergebnisse numerisch.

Im Code ist ein Batch von Vektoren ein Tensor der Form `(batch, n)`, auf den eine Matrix $T$ als `x @ T.T` wirkt.

<a id="s-12-1"></a>

## 12.1 Äquivariante lineare Abbildungen

*Vorlesungsnotizen:* [Abschnitt 12.1](../../lecture_notes/de/12_equivariant_learning.md#s-12-1): [Lemma 12.1](../../lecture_notes/de/12_equivariant_learning.md#r-12-1) (Bahnen) und [Korollar 12.2](../../lecture_notes/de/12_equivariant_learning.md#r-12-2) (DeepSets, Bell-Zahlen).

Die DeepSets-Schicht (Zaheer et al. 2017) ist [Vorlesungsnotizen, Korollar 12.2](../../lecture_notes/de/12_equivariant_learning.md#r-12-2)(a), die Basen von Maron et al. (2019) für Graphennetze sind (b).
Der Code zählt die Bahnen durch vollständige Aufzählung, vergleicht mit der numerisch berechneten Kommutante ([Lemma 12.1](../../lecture_notes/de/12_equivariant_learning.md#r-12-1))
und zerlegt beide Darstellungen von $S_4$:

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

In der von `dc.from_group` berechneten Zerlegung (die asserts prüfen ihre Ausgabe) enthält $\mathbb C^{16}$ also
eine Darstellung der Dimension $1$ zweimal – die triviale, deren Vielfachheit nach Proposition 12.7(a) und Lemma 12.1
die Anzahl $2$ der Bahnen auf $[4]^2$ ist (der Code in Abschnitt 12.4 identifiziert diesen Sektor) –, eine der
Dimension $2$ einmal und zwei der Dimension $3$, die eine einmal und die andere dreimal:
$16=1\cdot2+2\cdot1+3\cdot1+3\cdot3$, und $15=2^2+1^2+1^2+3^2$ stimmt mit der Bahnenzählung aus Lemma 12.1 und mit
$\dim\pi(G)'=\sum_cm_c^2$ (Satz 2.6) überein. Für zwei Darstellungen reduziert sich alles auf eine einzige
Kommutante:

Nach [Vorlesungsnotizen, Lemma 12.3](../../lecture_notes/de/12_equivariant_learning.md#r-12-3) ist $\mathrm{Hom}_G(\pi_1,\pi_2)=Q_2^*\,\pi(G)'\,Q_1$ für $\pi=\pi_1\oplus\pi_2$, und nach
[Satz 12.4](../../lecture_notes/de/12_equivariant_learning.md#r-12-4) ist es in geeigneten Basen $\bigoplus_cM_{b_c\times a_c}(\mathbb C)$, der Dimension $\sum_\rho m_\rho(\pi_1)\,m_\rho(\pi_2)$.

<a id="s-12-2"></a>

## 12.2 Eine äquivariante Schicht mit der Bibliothek

*Vorlesungsnotizen:* [Abschnitt 12.2](../../lecture_notes/de/12_equivariant_learning.md#s-12-2): [Satz 12.5](../../lecture_notes/de/12_equivariant_learning.md#r-12-5) (der Twirl); die Schichten folgen den Beweisen von [Lemma 12.3](../../lecture_notes/de/12_equivariant_learning.md#r-12-3) und [Satz 12.4](../../lecture_notes/de/12_equivariant_learning.md#r-12-4).

Der Beweis von Satz 12.4 ist ein Algorithmus, sobald die Zerlegung aus Satz 2.6 bekannt ist; `dc.from_group`
berechnet sie mit dem randomisierten Verfahren aus Abschnitt 4.3 (in exakter Arithmetik mit Wahrscheinlichkeit eins
korrekt; hier prüfen die asserts die berechnete Ausgabe). Für $\pi_1=\pi_2$ ist es einfach Satz 2.6: Die
Gewichte sind Blöcke $w_c\in M_{m_c}$, ein Operator der abstrakten Kommutante `D.commutant_algebra()`
$=\bigoplus_cM_{m_c}$, und die Schicht ist $x\mapsto V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*x$, berechnet von
`D.from_commutant`. Das folgende Modul ist Code des Kapitels, nicht Teil der Bibliothek.

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

Für $\pi_1\ne\pi_2$ folgt man dem Beweis Zeile für Zeile: $\pi_1\oplus\pi_2$ zerlegen, die Blöcke $p_c$ von $P_1$ mit
`to_commutant` berechnen, sie diagonalisieren und die linken unteren Ecken behalten. Die Spur von $p_c$ ist
$a_c=m_c(\pi_1)$.

Für Korollar 12.2(b) mit $n=4$, in der berechneten Zerlegung (unten ausgegeben und geprüft): Die Knotenmerkmale
$\pi_1$ haben $a=(1,0,0,1)$ in den vier Sektoren von $\pi_1\oplus\pi_2$ (Dimensionen $1,2,3,3$; der erste ist der
triviale, da nach Korollar 12.2(a) der einzige eindimensionale Bestandteil von $\pi_1$ trivial ist), die
Kantenmerkmale $\pi_2$ haben $b=(2,1,1,3)$, also $\dim\mathrm{Hom}_G=1\cdot2+1\cdot3=5=B_3$, wie
Korollar 12.2(b) beweist.

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

Die zweite Art, eine Schicht äquivariant zu machen, ist die Mittelung über die Gruppe. Der **Twirl** ist
$$
\mathcal P(W)=\frac1{|G|}\sum_{g\in G}\pi_2(g)\,W\,\pi_1(g)^*,\qquad W\in M_{n_2\times n_1}(\mathbb C).
$$

Nach [Vorlesungsnotizen, Satz 12.5](../../lecture_notes/de/12_equivariant_learning.md#r-12-5) ist $\mathcal P$ die bezüglich des Hilbert–Schmidt-Skalarprodukts orthogonale Projektion auf $\mathrm{Hom}_G(\pi_1,\pi_2)$, und
$\dim\mathrm{Hom}_G(\pi_1,\pi_2)=|G|^{-1}\sum_g\overline{\chi_{\pi_1}(g)}\chi_{\pi_2}(g)$; für $\pi_1=\pi_2$ wertet `dc.commutant_dimension` diese
Charakterformel aus. Die Parametrisierung $W\mapsto\mathcal P(W)$ verschwendet $n_1n_2-\dim\mathrm{Hom}_G$ Parameter an den Kern ($256-15$ für Kantenmerkmale), und eine
Auswertung kostet $O(|G|n_1n_2(n_1+n_2))$ für dichte $\pi(g)$, mit $|G|=n!$ für $S_n$. Die Bibliothek berechnet dieselbe
Projektion aus der Zerlegung, `D.conditional_expectation(W, 'commutant')`, mit vier $n\times n$-Produkten.
Der Code prüft, dass `dc.twirl` und `D.conditional_expectation` übereinstimmen, dass $\mathcal P$ auf $\pi(G)'$ die Identität
und selbstadjungiert ist, sowie die Dimensionsformel von Knoten- zu Kantenmerkmalen.

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

## 12.3 Eingeschränkte und spektrale Schichten

*Vorlesungsnotizen:* [Abschnitt 12.3](../../lecture_notes/de/12_equivariant_learning.md#s-12-3): [Proposition 12.6](../../lecture_notes/de/12_equivariant_learning.md#r-12-6).

Die Gewichte $w=(w_c)_c$ von `EquivariantLinear` leben in einer $*$-Algebra, also lässt sich jede Konstruktion der
Kapitel 1 und 11 auf sie anwenden. Man schreibe $\Phi(w)=V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*$ für
$w\in\bigoplus_cM_{m_c}$.

Nach [Vorlesungsnotizen, Proposition 12.6](../../lecture_notes/de/12_equivariant_learning.md#r-12-6) gilt: (a) $\Phi$ ist ein $*$-Isomorphismus von $\bigoplus_cM_{m_c}$ auf $\pi(G)'$, und $\Phi(w)$ ist
genau dann selbstadjungiert, unitär bzw. positiv, wenn jedes $w_c$ es ist; (b) $f(\Phi(w))=\Phi\big((f(w_c))_c\big)$ für normales $w$; (c) die
äquivarianten unitären Matrizen sind die $\Phi(u_0e^{s-s^*})$.

Eine orthogonale (unitäre) äquivariante Schicht ist also `optimize.UnitaryParam` auf `D.commutant_algebra()`, gefolgt
von $\Phi$; nach (c) erreicht diese Parametrisierung jede äquivariante unitäre Matrix, und in exakter Arithmetik ist
jede Iterierte eines Gradientenverfahrens unitär und äquivariant (der Code prüft beides nach dem Training bis auf
$10^{-12}$). Ebenso liefert `PositiveParam` positive äquivariante
Schichten, normiert wie in Kapitel 11 (Spektrum in $[\ell,1]$ mit Maximum $1$; nach (a) ist eine allgemeine positive
Schicht $\Phi(t^*t)$), und eine spektrale Schicht $f(L)$ eines selbstadjungierten äquivarianten Generators $L$, etwa
die Diffusion $e^{-tL}$, benötigt Eigenzerlegungen der Blöcke $w_c$ mit einem Aufwand der Größenordnung
$\sum_cm_c^3$ (hier $37$) statt $n^3=4096$; die Bibliothek differenziert durch `apply_function` hindurch (über ihren
Rückwärtsdurchlauf beweisen wir nichts; der Code prüft unten eine Ableitung). Der Code passt mit Adam eine
äquivariante unitäre Matrix an $32$ Eingabe-Ausgabe-Paare an – eine numerische Beobachtung dieses Laufs, kein Satz:
Der Verlust sinkt um mehr als einen Faktor $10^3$, was das assert prüft, und es wird keine Konvergenzaussage
behauptet – und prüft (b) sowie die Ableitung nach $t$ an dichten Matrizen und an der Formel
$\frac d{dt}e^{-tL}=-Le^{-tL}$ (für $L=U\mathrm{diag}(\lambda)U^*$ ist
$e^{-tL}=U\mathrm{diag}(e^{-t\lambda_i})U^*$ nach Satz 1.5, und $U$ hängt nicht von $t$ ab).

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

## 12.4 Invarianten

*Vorlesungsnotizen:* [Abschnitt 12.4](../../lecture_notes/de/12_equivariant_learning.md#s-12-4): [Proposition 12.7](../../lecture_notes/de/12_equivariant_learning.md#r-12-7) (invariante Vektoren und quadratische Formen) und [Proposition 12.8](../../lecture_notes/de/12_equivariant_learning.md#r-12-8) (Molien-Formel).

Für Graphen werden die invarianten Kantenmerkmale von den Indikatoren der Diagonale und ihres Komplements aufgespannt
($B_2=2$), und die invarianten linearen Ausleseabbildungen (*readouts*) sind die Summe der Diagonaleinträge und die
der Nichtdiagonaleinträge. Der Code findet den trivialen Sektor aus den Blöcken $\rho_c(g)$ (`to_algebra`) und prüft
(a), (b) und die beiden Abzählungen aus (c).

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

Auf $\mathbb C^4$ sind es $\sum_ix_i^2$ und $\sum_{i<j}x_ix_j$. Der Beweis von (c) funktioniert in jedem Grad:

Nach [Vorlesungsnotizen, Proposition 12.8](../../lecture_notes/de/12_equivariant_learning.md#r-12-8) (Molien-Formel) gilt
$\sum_{d\ge0}\dim\mathbb C[x]^G_d\,t^d=|G|^{-1}\sum_g\det\big(1-t\pi(g)\big)^{-1}$, und für eine Permutationsdarstellung ist
$\dim\mathbb C[x]^G_d$ die Anzahl der Bahnen von $G$ auf den Monomen vom Grad $d$.

(Das ist der Satz von Molien aus dem Jahr 1897; zur Vertiefung siehe B. Sturmfels, *Algorithms in Invariant
Theory*, Springer 2008, Abschnitt 2.2.) Der Code entwickelt die rechte Seite bis $t^4$, indem er mit einer
geometrischen Reihe pro Eigenwert multipliziert, und vergleicht mit den Bahnen der Monome; für $S_4$ auf
$\mathbb C^4$ sind die Anzahlen $1,1,2,3,5$ die Partitionen von $d$ in höchstens vier Teile (die Bahnen von Multimengen
sind durch ihre sortierten Vielfachheiten festgelegt).

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

## 12.5 Womit man vergleichen sollte

Äquivariante Schichten für spezielle Gruppen stellen ausgereifte Bibliotheken bereit. **escnn** (Cesa, Lang und
Weiler 2022, Nachfolger von e2cnn) baut steuerbare (*steerable*) Faltungen für Untergruppen von $O(2)$ und $O(3)$ aus
irreduziblen Darstellungen und harmonischen Basen von Faltungskernen; **e3nn** (Geiger und Smidt 2022) baut
$E(3)$-äquivariante Netze aus Kugelflächenfunktionen und Clebsch–Gordan-Tensorprodukten. **EMLP** (Finzi, Welling und
Wilson 2021) behandelt beliebige Matrixgruppen: Es fordert $\pi_2(g)T=T\pi_1(g)$ für diskrete Erzeuger und
$d\pi_2(X)T=Td\pi_1(X)$ für Erzeuger der Lie-Algebra und berechnet eine Basis der Lösungen als numerischen Nullraum –
das lineare Problem von `dc.commutant_basis`. Gruppenfaltungen (Cohen und Welling 2016) sind der Fall
$\Omega_1=\Omega_2=G$ von Lemma 12.1, wobei $G$ durch Linksmultiplikation auf sich selbst wirkt: $(y,x)$ und $(y',x')$
liegen genau dann in derselben Bahn, wenn $x^{-1}y=x'^{-1}y'$ (dann bildet $g=x'x^{-1}$ das eine Paar auf das andere
ab), also sind die äquivarianten Abbildungen $(Tv)(y)=\sum_xf(x^{-1}y)v(x)$, die Faltungen mit Funktionen $f$ auf $G$.
Kondor und Trivedi (2018) untersuchen die analoge Frage für homogene Räume kompakter Gruppen, die dieser Kurs nicht
behandelt. Für Permutationen sind DeepSets (Zaheer et al.
2017) und die Basen von Maron et al. (2019) Korollar 12.2.

Die Bibliothek bietet keine Faltungen auf Gittern, keine steuerbaren Kerne auf kontinuierlichen Räumen, keine
Zerlegungen über $\mathbb R$ (reeller, komplexer und quaternionischer Typ; Aufgabe 2 behandelt reelle Matrizen) und
kein großes $n$: $V$ ist eine dichte $n\times n$-Matrix, und die Zerlegung kostet $O(|G|n^3)$ pro Zufallsstichprobe
für eine durch all ihre Elemente gegebene Gruppe (Erzeuger von Lie-Algebren laufen über `dc.from_generators`;
Abschnitt 4.3). Was sie hinzufügt, ist die **Struktur** des Lösungsraums: Eine Basis der Kommutante ist eine Liste von
Matrizen, während die Zerlegung die Kommutante mit der $*$-Algebra $\bigoplus_cM_{m_c}$ identifiziert. Produkte,
Adjungierte, Unitarität, Positivität und der Funktionalkalkül wirken dann blockweise auf kleinen Matrizen
(Proposition 12.6), die Bibliothek differenziert durch sie hindurch, und `to_commutant`/`from_commutant` übersetzen
zwischen den beiden Bildern.

<a id="x-12"></a>

## Aufgaben

[Aufgaben](../../lecture_notes/de/12_equivariant_learning.md#x-12) 1, 3 und 4 der Vorlesungsnotizen, numerisch:

- *Aufgabe 1* (DeepSets mit Kanälen): Man prüfe mit `EquivariantMap`, dass die äquivarianten Abbildungen von
  $\mathbb C^n\otimes\mathbb C^{d_1}$ nach $\mathbb C^n\otimes\mathbb C^{d_2}$ die $1\otimes A+J\otimes B$ sind.
- *Aufgabe 3* (Nichtlinearitäten): Koordinatenweises ReLU vertauscht mit der Permutationsdarstellung, aber nicht
  mit der Drehung von $\mathbb R^2$ um $\pi/2$; $x\mapsto\varphi(\|x\|)x$ vertauscht mit jeder unitären Darstellung.
- *Aufgabe 4*: Für die Kantendarstellung lässt eine Gaußsche Stichprobe einen $2$-dimensionalen Raum äquivarianter
  Schichten unbestimmt, zwei Stichproben legen die Schicht fest.

<details>
<summary>Lösungen</summary>

*Aufgabe 1.*

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

*Aufgabe 3.*

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

*Aufgabe 4.*

```python
basis = D_edge.basis('commutant')                                # (15, 16, 16), Hilbert-Schmidt orthonormal
for N, lost in ((1, 2), (2, 0)):
    Xs = torch.randn(16, N, dtype=dt)
    M = torch.stack([(B @ Xs).reshape(-1) for B in basis], dim=1)               # w -> (T x_1, ..., T x_N)
    assert torch.linalg.matrix_rank(M).item() == 15 - lost
```

</details>

**Weiter:** [Tutorial 13. Jordan-Algebren: warum die Quantenmechanik *-Algebren verwendet](13_jordan_algebras.md).
