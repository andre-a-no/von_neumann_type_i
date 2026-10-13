[English](../12_equivariant_learning.md) · [Español](../es/12_equivariant_learning.md) · [Français](../fr/12_equivariant_learning.md) · **Deutsch** · [中文](../zh/12_equivariant_learning.md) · [日本語](../ja/12_equivariant_learning.md) · [Русский](../ru/12_equivariant_learning.md)

# 12. Äquivariantes maschinelles Lernen

Eine Schicht eines neuronalen Netzes heißt **äquivariant**, wenn man dasselbe Ergebnis erhält, ob man die Eingabe mit
einer Symmetrie transformiert und dann die Schicht anwendet oder die Schicht anwendet und dann die Ausgabe
transformiert: Ein Umnummerieren der Knoten eines Graphen nummeriert seine Merkmale auf dieselbe Weise um. Für
lineare Schichten ist das eine Vertauschungsrelation, also bilden die äquivarianten Schichten eine Kommutante, die
nach Kapitel 2 bis auf einen Basiswechsel die Algebra $\bigoplus_cM_{m_c}(\mathbb C)$ ist. Dieses Kapitel zählt ihre
Parameter, baut und trainiert sie mit `torch_vn_algebra.decompose` und nutzt die Algebrenstruktur für unitäre,
positive und spektrale Schichten, die nach Konstruktion äquivariant sind.

**Sie lernen**, warum die äquivarianten linearen Abbildungen zwischen zwei Darstellungen Blockmatrizen sind, deren
Größen durch die Vielfachheiten gegeben sind; wie man sie über Bahnen oder über Charaktere abzählt (DeepSets und die
Bell-Zahlen von Graphennetzen); warum die Mittelung über die Gruppe die Orthogonalprojektion auf die äquivarianten
Abbildungen ist; wie Funktionalkalkül und unitäre Parametrisierungen blockweise auf der Kommutante wirken; und wie
man invariante Vektoren und invariante quadratische Formen abzählt.

Im Folgenden ist $G$ eine endliche Gruppe, und eine **Darstellung** von $G$ auf $\mathbb C^n$ ist ein Homomorphismus
$\pi:G\to M_n(\mathbb C)$ mit unitären Werten (das ist keine Einschränkung: Invertierbare $\pi(g)$ sind unitär
bezüglich des Skalarprodukts $\sum_g\langle\pi(g)x,\pi(g)y\rangle$, Abschnitt 4.1). Ihr **Charakter** ist
$\chi_\pi(g)=\operatorname{Tr}\pi(g)$. Irreduzible und äquivalente Darstellungen sowie die **Vielfachheit**
$m_\rho(\pi)$ einer irreduziblen $\rho$ in $\pi$ ($0$, falls $\rho$ nicht auftritt) sind wie in Abschnitt 2.4 erklärt,
wo bewiesen wird, dass $m_\rho(\pi)$ nicht von der Zerlegung abhängt, und wo das Lemma von Schur bewiesen wird;
[Kapitel 4](04_symmetries.md) behandelt beides ausführlich (Abschnitte 4.1–4.2). Auf Rechteckmatrizen ist
$\langle S,T\rangle=\operatorname{Tr}(S^*T)$ das **Hilbert–Schmidt**-Skalarprodukt. Im Code ist ein Batch von Vektoren
ein Tensor der Form `(batch, n)`, auf den eine Matrix $T$ als `x @ T.T` wirkt.

## 12.1 Äquivariante lineare Abbildungen

**Definition.** Seien $\pi_1,\pi_2$ Darstellungen von $G$ auf $\mathbb C^{n_1}$, $\mathbb C^{n_2}$. Eine Matrix
$T\in M_{n_2\times n_1}(\mathbb C)$ heißt **äquivariant** (ein **Verflechtungsoperator**), wenn $T\pi_1(g)=\pi_2(g)T$
für alle $g\in G$ gilt. Diese Matrizen bilden einen Vektorraum $\operatorname{Hom}_G(\pi_1,\pi_2)$; für
$\pi_1=\pi_2=\pi$ ist er die Kommutante, $\operatorname{Hom}_G(\pi,\pi)=\pi(G)'$.

Die im maschinellen Lernen häufigsten Darstellungen permutieren Koordinaten. Die Gruppe $G$ wirke auf einer endlichen
Menge $\Omega$; die **Permutationsdarstellung** auf $\mathbb C^\Omega$ ist $\pi(g)e_x=e_{gx}$. Beispiele: $S_n$, die die
$n$ Punkte einer Menge permutiert ($\Omega=[n]=\{1,\dots,n\}$), und $S_n$, die die Knoten eines Graphen umnummeriert
und dabei auf $\Omega=[n]^2$ wirkt, d. h. auf die Adjazenzmatrix durch $X\mapsto P_\sigma XP_\sigma^T$.

**Lemma 12.1 (Bahnen).** Seien $\pi_i$ die Permutationsdarstellungen von $G$ auf $\mathbb C^{\Omega_i}$. Dann gilt
$T\in\operatorname{Hom}_G(\pi_1,\pi_2)$ genau dann, wenn $T_{gy,gx}=T_{yx}$ für alle $g$, $y\in\Omega_2$,
$x\in\Omega_1$. Die Indikatormatrizen der Bahnen von $G$ auf $\Omega_2\times\Omega_1$ bilden eine Basis, und
$\dim\operatorname{Hom}_G(\pi_1,\pi_2)$ ist die Anzahl dieser Bahnen.

*Beweis.* $\pi_2(g)T\pi_1(g)^*e_x=\pi_2(g)Te_{g^{-1}x}=\sum_yT_{y,g^{-1}x}e_{gy}$, also ist der $(gy,x)$-Eintrag von
$\pi_2(g)T\pi_1(g)^*$ gleich $T_{y,g^{-1}x}$. Daher ist $T$ genau dann äquivariant, d. h. $\pi_2(g)T\pi_1(g)^*=T$ für
alle $g$, wenn $T_{gy,x}=T_{y,g^{-1}x}$ für alle $g,y,x$ gilt, und das ist (mit $gx$ statt $x$) die angegebene
Bedingung: $T$ ist auf den Bahnen konstant. Die Indikatormatrizen der Bahnen haben disjunkte Träger, sind also linear
unabhängig, und jede auf den Bahnen konstante Funktion ist eine Linearkombination von ihnen. $\square$

**Korollar 12.2.** (a) (*DeepSets.*) Für $S_n$ auf $\mathbb C^n$, $n\ge2$, ist die Kommutante
$\operatorname{span}\{1,J\}$, wobei $J$ die Matrix ist, deren Einträge alle gleich eins sind, und
$\mathbb C^n=\mathbb C\mathbf 1\oplus\mathbf 1^\perp$ ist eine Summe zweier inäquivalenter irreduzibler
Darstellungen, der trivialen und der **Standarddarstellung** der Dimension $n-1$. (b) (*Graphennetze.*) Die
äquivarianten Abbildungen von $\mathbb C^{[n]^p}$ nach $\mathbb C^{[n]^q}$ (Tensoren der Stufe $p$ und $q$ auf $n$
Knoten) werden durch die Partitionen einer $(p+q)$-elementigen Menge in höchstens $n$ Blöcke abgezählt; für
$n\ge p+q$ ist ihre Anzahl die **Bell-Zahl** $B_{p+q}$ ($B_2=2$, $B_3=5$, $B_4=15$).

*Beweis.* (b) Zwei Tupel $x,y\in[n]^{p+q}$ liegen genau dann in derselben Bahn, wenn sie dasselbe
**Gleichheitsmuster** haben, d. h. dieselbe Partition der Positionen $\{1,\dots,p+q\}$ durch
$i\sim j\Leftrightarrow x_i=x_j$. Das Muster ist offensichtlich invariant. Haben umgekehrt $x,y$ dasselbe Muster, so
ist $x_i\mapsto y_i$ eine wohldefinierte Injektion von $\{x_i\}$ in $[n]$, die sich zu einer Permutation $\sigma$ mit
$\sigma x=y$ fortsetzt. Jede Partition mit höchstens $n$ Blöcken ist ein Muster (man gebe den Blöcken verschiedene
Werte). Lemma 12.1 schließt den Beweis ab. (a) Für $p=q=1$ gibt es zwei Bahnen auf $[n]^2$, die Diagonale und ihr
Komplement, mit den Indikatoren $1$ und $J-1$. Also ist $\dim\pi(S_n)'=2=\sum_cm_c^2$ (Satz 2.6, Proposition 2.7),
was zwei Sektoren mit $m_c=1$ erzwingt, d. h. zwei inäquivalente irreduzible Summanden, und $\pi(S_n)'$ wird von
deren Projektionen $P_1,P_2$ aufgespannt. Die Projektion auf die invariante Gerade $\mathbb C\mathbf 1$ liegt in
$\pi(S_n)'$ (Lemma 2.1), ist also $P_1$ oder $P_2$: Die Gerade, die die triviale Darstellung trägt, ist ein Sektor,
und $\mathbf 1^\perp$ ist der andere. $\square$

Die DeepSets-Schicht $x\mapsto\lambda x+\gamma\,(\mathbf 1^Tx)\mathbf 1$ (Zaheer et al. 2017) ist (a); die Basen von
Maron et al. (2019) für Graphennetze sind (b). Der Code zählt die Bahnen durch vollständige Aufzählung, vergleicht mit
der numerisch berechneten Kommutante und zerlegt beide Darstellungen von $S_4$:

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

Also enthält $\mathbb C^{16}$ die triviale Darstellung zweimal, eine der Dimension $2$ einmal und zwei der Dimension
$3$, die eine einmal und die andere dreimal: $16=1\cdot2+2\cdot1+3\cdot1+3\cdot3$ und $15=2^2+1^2+1^2+3^2$. Für zwei
Darstellungen reduziert sich alles auf eine einzige Kommutante:

**Lemma 12.3 (direkte Summe).** Sei $\pi=\pi_1\oplus\pi_2$ auf $\mathbb C^{n_1}\oplus\mathbb C^{n_2}$, mit den
isometrischen Inklusionen $Q_i$ und den Projektionen $P_i=Q_iQ_i^*$. Dann gilt $P_1,P_2\in\pi(G)'$,
$\operatorname{Hom}_G(\pi_1,\pi_2)=Q_2^*\,\pi(G)'\,Q_1$, und $T\mapsto Q_2TQ_1^*$ ist eine lineare Bijektion von
$\operatorname{Hom}_G(\pi_1,\pi_2)$ auf $P_2\,\pi(G)'P_1$.

*Beweis.* $\pi(g)Q_i=Q_i\pi_i(g)$ und $Q_i^*\pi(g)=\pi_i(g)Q_i^*$, weil $\pi(g)$ blockdiagonal ist; daher gilt
$P_i\pi(g)=Q_i\pi_i(g)Q_i^*=\pi(g)P_i$. Ist $X\in\pi(G)'$, so gilt
$Q_2^*XQ_1\pi_1(g)=Q_2^*X\pi(g)Q_1=Q_2^*\pi(g)XQ_1=\pi_2(g)Q_2^*XQ_1$. Ist $T$ äquivariant, so gilt
$\pi(g)Q_2TQ_1^*=Q_2\pi_2(g)TQ_1^*=Q_2T\pi_1(g)Q_1^*=Q_2TQ_1^*\pi(g)$, also $Q_2TQ_1^*\in\pi(G)'$; diese Matrix ist
gleich $P_2(Q_2TQ_1^*)P_1$, und $Q_2^*(Q_2TQ_1^*)Q_1=T$. Umgekehrt liegt jedes $Y=P_2XP_1$ mit $X\in\pi(G)'$ in
$\pi(G)'$, und $Y=Q_2(Q_2^*YQ_1)Q_1^*$. $\square$

**Satz 12.4 (äquivariante Abbildungen sind Blockmatrizen).** Seien $\rho_1,\dots,\rho_C$ die paarweise inäquivalenten
irreduziblen Darstellungen, die in $\pi_1\oplus\pi_2$ auftreten, $k_c=\dim\rho_c$, $a_c=m_{\rho_c}(\pi_1)$,
$b_c=m_{\rho_c}(\pi_2)$. Es gibt unitäre Abbildungen
$V_1:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{a_c}\to\mathbb C^{n_1}$ und
$V_2:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{b_c}\to\mathbb C^{n_2}$ mit
$V_i^*\pi_i(g)V_i=\bigoplus_c\rho_c(g)\otimes1$, sodass
$$
\operatorname{Hom}_G(\pi_1,\pi_2)=\Big\{V_2\Big(\bigoplus_{c=1}^C1_{k_c}\otimes y_c\Big)V_1^*:\ y_c\in M_{b_c\times a_c}(\mathbb C)\Big\},
$$
und $(y_c)_c\mapsto T$ ist injektiv. Folglich ist $\operatorname{Hom}_G(\pi_1,\pi_2)\cong\bigoplus_cM_{b_c\times a_c}(\mathbb C)$,
$$
\dim\operatorname{Hom}_G(\pi_1,\pi_2)=\sum_\rho m_\rho(\pi_1)\,m_\rho(\pi_2),
$$
wobei die Summe über alle irreduziblen Darstellungen bis auf Äquivalenz läuft; für $\pi_1=\pi_2$ ist das
$\pi(G)'\cong\bigoplus_cM_{m_c}$ (Proposition 2.7), und $\operatorname{Hom}_G(\pi_1,\pi_2)=0$ gilt genau dann, wenn
keine irreduzible Darstellung in beiden auftritt.

*Beweis.* Man wende Satz 2.6 und Proposition 2.7 auf $\pi=\pi_1\oplus\pi_2$ an: Es gibt eine unitäre Matrix $V$ mit
$V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$ und $V^*\pi(G)'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$. Da
$X\mapsto V^*XV$ ein $*$-Isomorphismus ist und $P_1\in\pi(G)'$ (Lemma 12.3), gilt $V^*P_1V=\bigoplus_c1\otimes p_c$ mit
Projektionen $p_c$. Man wähle unitäre $u_c$ mit $u_c^*p_cu_c=1_{a_c}\oplus0_{b_c}$ (Satz 1.3), wobei nun
$a_c=\operatorname{rank}p_c$, $b_c=m_c-a_c$, und ersetze $V$ durch $V\bigoplus_c(1\otimes u_c)$; das ändert keine der
beiden Normalformen, da $1\otimes u_c$ mit $\rho_c(g)\otimes1$ vertauscht. Seien
$\iota_c:\mathbb C^{a_c}\to\mathbb C^{m_c}$ und $\kappa_c:\mathbb C^{b_c}\to\mathbb C^{m_c}$ die Inklusionen der ersten
$a_c$ bzw. der letzten $b_c$ Koordinaten, und
$$
V_1=Q_1^*V\bigoplus_c(1\otimes\iota_c),\qquad V_2=Q_2^*V\bigoplus_c(1\otimes\kappa_c).
$$
$W_1=V\bigoplus_c(1\otimes\iota_c)$ ist eine Isometrie mit Bild $\operatorname{range}P_1=\operatorname{range}Q_1$, also
$V_1^*V_1=W_1^*P_1W_1=1$, und $V_1$ ist unitär, weil $\sum_ck_ca_c=\operatorname{rank}P_1=n_1$. Außerdem gilt
$\pi_1(g)V_1=Q_1^*\pi(g)W_1=Q_1^*V\big(\bigoplus_c\rho_c(g)\otimes1_{m_c}\big)\bigoplus_c(1\otimes\iota_c)
=V_1\bigoplus_c\rho_c(g)\otimes1_{a_c}$. Also ist $\pi_1\cong\bigoplus_c\rho_c\otimes1_{a_c}$ mit inäquivalenten
irreduziblen $\rho_c$, und $a_c=m_{\rho_c}(\pi_1)$ wegen der Eindeutigkeit der Vielfachheiten (Abschnitt 2.4); ebenso
für $V_2$ und $b_c$. Schließlich ist für $X=V\big(\bigoplus_c1\otimes x_c\big)V^*\in\pi(G)'$
$P_2XP_1=V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*$, wobei $y_c=\kappa_c^*x_c\iota_c$ die linke untere
$b_c\times a_c$-Ecke von $x_c$ ist, eine beliebige Matrix. Nach Lemma 12.3 besteht
$\operatorname{Hom}_G(\pi_1,\pi_2)=Q_2^*P_2\,\pi(G)'P_1Q_1$ aus den Matrizen
$Q_2^*V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*Q_1=V_2\big(\bigoplus_c1\otimes y_c\big)V_1^*$, und
$(y_c)\mapsto T$ ist injektiv, da $V_1,V_2$ unitär sind. Eine irreduzible $\rho$, die in $\pi_1\oplus\pi_2$ nicht
auftritt, hat $m_\rho(\pi_1)m_\rho(\pi_2)=0$. $\square$

Für Korollar 12.2(b) mit $n=4$: Die Knotenmerkmale $\pi_1$ haben $a=(1,0,0,1)$ in den vier Sektoren von
$\pi_1\oplus\pi_2$ (trivial, Dimension $2$, die beiden der Dimension $3$), die Kantenmerkmale $\pi_2$ haben
$b=(2,1,1,3)$, also $\dim\operatorname{Hom}_G=1\cdot2+1\cdot3=5=B_3$.

> **Physikalische Bedeutung.** Satz 12.4 ist der linear-algebraische Kern der **Auswahlregeln** und der skalare Fall
> des Wigner–Eckart-Theorems: Ein Operator, der mit der Symmetrie vertauscht, kann nur Zustände desselben
> Symmetrietyps verbinden, und innerhalb eines Typs wirkt er als $1_{k_c}\otimes y_c$, d. h. durch ein „reduziertes
> Matrixelement“ $y_c$, das blind für die Position innerhalb des Multipletts ist. (Operatoren, die sich selbst nach
> einer nichttrivialen Darstellung transformieren, etwa Vektoroperatoren, benötigen zusätzlich die
> Clebsch–Gordan-Koeffizienten.) Ein symmetrischer Hamiltonoperator ist der Fall $\pi_1=\pi_2$; seine Entartungen
> $k_c$ sind die Dimensionen der irreduziblen Darstellungen.

## 12.2 Eine äquivariante Schicht mit der Bibliothek

Der Beweis von Satz 12.4 ist ein Algorithmus, sobald die Zerlegung aus Satz 2.6 bekannt ist; `dc.from_group`
berechnet sie mit dem randomisierten Verfahren aus Abschnitt 4.3. Für $\pi_1=\pi_2$ ist es einfach Satz 2.6: Die
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

Die zweite Art, eine Schicht äquivariant zu machen, ist die Mittelung über die Gruppe. Der **Twirl** ist
$$
\mathcal P(W)=\frac1{|G|}\sum_{g\in G}\pi_2(g)\,W\,\pi_1(g)^*,\qquad W\in M_{n_2\times n_1}(\mathbb C).
$$

**Satz 12.5 (der Twirl ist die Orthogonalprojektion).** $\mathcal P$ ist die bezüglich des Hilbert–Schmidt-Skalarprodukts
orthogonale Projektion von $M_{n_2\times n_1}(\mathbb C)$ auf $\operatorname{Hom}_G(\pi_1,\pi_2)$, und
$$
\dim\operatorname{Hom}_G(\pi_1,\pi_2)=\frac1{|G|}\sum_{g\in G}\overline{\chi_{\pi_1}(g)}\,\chi_{\pi_2}(g).
$$

*Beweis.* Für $h\in G$ ist $\pi_2(h)\mathcal P(W)\pi_1(h)^*=\frac1{|G|}\sum_g\pi_2(hg)W\pi_1(hg)^*=\mathcal P(W)$, da
$g\mapsto hg$ eine Bijektion von $G$ ist; also ist $\mathcal P(W)$ äquivariant. Ist $W$ äquivariant, so ist jeder
Summand gleich $W$, also $\mathcal P(W)=W$. Daher gilt $\mathcal P^2=\mathcal P$ mit Bild
$\operatorname{Hom}_G(\pi_1,\pi_2)$. Die Abbildung ist selbstadjungiert:
$\langle S,\mathcal P(W)\rangle=\frac1{|G|}\sum_g\operatorname{Tr}\big((\pi_2(g)^*S\pi_1(g))^*W\big)$
$=\langle\mathcal P(S),W\rangle$, wobei wir $g\mapsto g^{-1}$ substituiert und $\pi_i(g^{-1})=\pi_i(g)^*$ benutzt haben.
Ein selbstadjungiertes Idempotent ist die Orthogonalprojektion auf sein Bild. Seine Spur ist sein Rang. Die lineare
Abbildung $W\mapsto AWB$ auf $M_{n_2\times n_1}$ hat die Spur $\operatorname{Tr}A\operatorname{Tr}B$ (der Koeffizient
der Matrixeinheit $E_{ij}$ in $AE_{ij}B$ ist $A_{ii}B_{jj}$); mit $A=\pi_2(g)$, $B=\pi_1(g)^*$ ergibt das die
Formel. $\square$

Zusammen mit Satz 12.4 folgt
$\sum_\rho m_\rho(\pi_1)m_\rho(\pi_2)=|G|^{-1}\sum_g\overline{\chi_{\pi_1}(g)}\chi_{\pi_2}(g)$; für $\pi_1=\pi_2$ ist das
die Charakterformel für $\dim\pi(G)'$ aus Abschnitt 4.2, die `dc.commutant_dimension` auswertet. Die Parametrisierung
$W\mapsto\mathcal P(W)$ erreicht jede äquivariante Abbildung, $\mathcal P(W)$ ist die zu $W$ nächstgelegene
äquivariante Abbildung, und da $\mathcal P$ selbstadjungiert ist, ist der Gradient nach $W$ das Bild unter
$\mathcal P$ des Gradienten nach $T=\mathcal P(W)$. Aber $n_1n_2-\dim\operatorname{Hom}_G$ Parameter werden an den
Kern verschwendet ($256-15$ für Kantenmerkmale), und eine Auswertung kostet $O(|G|n_1n_2(n_1+n_2))$ für dichte
$\pi(g)$, mit $|G|=n!$ für $S_n$. Die Bibliothek berechnet dieselbe Projektion aus der Zerlegung,
`D.conditional_expectation(W, 'commutant')`, mit vier $n\times n$-Produkten.

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

Die bisherigen Schichten sind linear. Ein Netz wechselt sie mit Nichtlinearitäten ab, und auch diese müssen
äquivariant sein. Für eine Permutationsdarstellung ist eine koordinatenweise angewandte Funktion äquivariant, da
$(\pi(g)x)_{gy}=x_y$; für andere Darstellungen ist sie es nicht, und man verwendet Nichtlinearitäten, die auf Normen
wirken oder auf den Koeffizienten bezüglich einer Basis, die $G$ permutiert (Aufgabe 3).

> **Physikalische Bedeutung.** Der Twirl ist die bedingte Erwartung auf die Kommutante (Kapitel 2, Aufgabe 3): Mittelt
> man eine Kopplung über alle Symmetrieoperationen, so bleibt genau ihr symmetrischer Anteil übrig. Einen
> Hamiltonoperator oder eine Wechselwirkung durch Summation über die Symmetriegruppe zu „symmetrisieren“ ist dieselbe
> Projektion.

## 12.3 Eingeschränkte und spektrale Schichten

Die Gewichte $w=(w_c)_c$ von `EquivariantLinear` leben in einer $*$-Algebra, also lässt sich jede Konstruktion der
Kapitel 1 und 11 auf sie anwenden. Man schreibe $\Phi(w)=V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*$ für
$w\in\bigoplus_cM_{m_c}$.

**Proposition 12.6.** (a) $\Phi$ ist ein $*$-Isomorphismus von $\bigoplus_cM_{m_c}$ auf $\pi(G)'$, und $\Phi(w)$ ist
genau dann selbstadjungiert, unitär bzw. positiv, wenn jedes $w_c$ es ist. (b) Für normales $w$ und eine Funktion $f$
auf $\bigcup_c\operatorname{spec}w_c$ gilt $f(\Phi(w))=\Phi\big((f(w_c))_c\big)$; insbesondere ist $f(X)\in\pi(G)'$ für
jedes normale $X\in\pi(G)'$. (c) Die äquivarianten unitären Matrizen sind die $\Phi(u_0e^{s-s^*})$, wobei $s$ ganz
$\bigoplus_cM_{m_c}$ durchläuft, für jedes feste unitäre $u_0\in\bigoplus_cM_{m_c}$; insbesondere sind es die
Exponentiale $e^S$ der schiefadjungierten $S\in\pi(G)'$.

*Beweis.* (a) $\Phi$ ist linear, unital, multiplikativ und mit der Involution verträglich, weil
$(1\otimes x)(1\otimes y)=1\otimes xy$ und $(1\otimes x)^*=1\otimes x^*$; sie ist injektiv und nach Satz 2.6 surjektiv
auf $\pi(G)'$. Daher gilt $\Phi(w)^*=\Phi(w)$ genau dann, wenn $w^*=w$, und $\Phi(w)^*\Phi(w)=1=\Phi(w)\Phi(w)^*$ genau
dann, wenn $w^*w=1=ww^*$. Sind alle $w_c=t_c^*t_c$, so ist $\Phi(w)=\Phi(t)^*\Phi(t)\ge0$ (Korollar 1.7). Ist
umgekehrt $\Phi(w)\ge0$, so gilt für $v\in\mathbb C^{m_c}$ und den Vektor $\xi=V(e_1\otimes v)$ im Sektor $c$
$0\le\langle\xi,\Phi(w)\xi\rangle=\langle v,w_cv\rangle$.
(b) Sei $w_c=u_c\Lambda_cu_c^*$ (Satz 1.3). Dann ist $\Phi(w)=U\big(\bigoplus_c1\otimes\Lambda_c\big)U^*$ mit der
unitären Matrix $U=V\bigoplus_c(1\otimes u_c)$, eine Diagonalisierung, in der die Eigenwerte von $w_c$ je $k_c$-mal
wiederholt auftreten; nach Definition und Satz 1.5(a) gilt
$f(\Phi(w))=U\big(\bigoplus_c1\otimes f(\Lambda_c)\big)U^*=\Phi\big((f(w_c))_c\big)$.
(c) Nach (a) sind die äquivarianten unitären Matrizen die $\Phi(u)$ mit unitärem $u$ in $\bigoplus_cM_{m_c}$, und
diese sind nach Lemma 11.1 die $u_0e^{s-s^*}$. Mit $u_0=1$ liefert (b) $\Phi(e^{s-s^*})=e^{\Phi(s-s^*)}$, und
$\Phi(s-s^*)$ durchläuft die schiefadjungierten Elemente von $\pi(G)'$. $\square$

Eine orthogonale (unitäre) äquivariante Schicht ist also `optimize.UnitaryParam` auf `D.commutant_algebra()`, gefolgt
von $\Phi$; nach (c) erreicht diese Parametrisierung jede äquivariante unitäre Matrix, und jede Iterierte eines
Gradientenverfahrens ist exakt unitär und exakt äquivariant. Ebenso liefert `PositiveParam` positive äquivariante
Schichten, normiert wie in Kapitel 11 (Spektrum in $[\ell,1]$ mit Maximum $1$; nach (a) ist eine allgemeine positive
Schicht $\Phi(t^*t)$), und eine spektrale Schicht $f(L)$ eines selbstadjungierten äquivarianten Generators $L$, etwa
die Diffusion $e^{-tL}$, benötigt Eigenzerlegungen mit Aufwand $\sum_cm_c^3$ (hier $37$) statt $n^3=4096$, mit
Gradienten nach der Daleckii–Krein-Formel (`apply_function`). Der Code passt eine äquivariante unitäre Matrix an $32$
Eingabe-Ausgabe-Paare an und prüft (b) sowie die Ableitung nach $t$, $\frac d{dt}e^{-tL}=-Le^{-tL}$, an dichten
Matrizen.

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

> **Physikalische Bedeutung.** Ein symmetrischer Hamiltonoperator $H\in\pi(G)'$ erzeugt eine symmetrische
> Zeitentwicklung $e^{-itH}\in\pi(G)'$ (Kapitel 8), die nie Sektoren mischt und nur auf den Vielfachheitsräumen wirkt.
> Eine äquivariante unitäre Schicht ist eine solche Entwicklung, gelernt statt vorgegeben.

## 12.4 Invarianten

Ein Vektor $x$ heißt **invariant**, wenn $\pi(g)x=x$ für alle $g$ gilt; eine Funktion $q$ auf $\mathbb C^n$ heißt
invariant, wenn $q(\pi(g)x)=q(x)$. Invariante Merkmale sind das, was ein invariantes Netz (etwa ein
Graphenklassifikator) ausgeben kann.

**Proposition 12.7.** (a) Die invarianten Vektoren bilden das Bild der Orthogonalprojektion
$P_0=|G|^{-1}\sum_g\pi(g)$; ihre Dimension ist $m_{\mathrm{triv}}(\pi)=|G|^{-1}\sum_g\chi_\pi(g)$, und wenn die triviale
Darstellung auftritt, als $\rho_{c_0}$, dann ist $P_0$ die minimale zentrale Projektion $P_{c_0}$ von $\pi(G)''$.
(b) Eine hermitesche Form $q_A(x)=\langle x,Ax\rangle$, $A=A^*$, ist genau dann invariant, wenn $A\in\pi(G)'$; die
invarianten hermiteschen Formen bilden einen reellen Vektorraum der Dimension $\sum_cm_c^2$.
(c) Die homogenen invarianten Polynome vom Grad $2$ in den Koordinaten $x_1,\dots,x_n$ bilden einen Raum der Dimension
$$
\frac1{|G|}\sum_{g\in G}\frac{\chi_\pi(g)^2+\chi_\pi(g^2)}2 .
$$
Für eine Permutationsdarstellung ist das die Anzahl der Bahnen von $G$ auf den Multimengen $\{x,y\}$ aus zwei
Elementen von $\Omega$.

*Beweis.* (a) Die invarianten Vektoren sind $\operatorname{Hom}_G(\mathrm{triv},\pi)$, mit der trivialen Darstellung
$\mathrm{triv}$ auf $\mathbb C=M_{1\times1}$; für sie ist der Twirl $\mathcal P(x)=P_0x$, also liefert Satz 12.5 die
Projektion und die Dimension, und Satz 12.4 liefert $m_{\mathrm{triv}}(\pi)$. In der Basis aus Satz 2.6 schreibe man
einen invarianten Vektor als $\xi=\sum_c\xi_c$ mit
$\xi_c=\sum_jv_{c,j}\otimes e_j\in\mathbb C^{k_c}\otimes\mathbb C^{m_c}$. Da $\pi(g)$ die Sektoren erhält, gilt
$\rho_c(g)v_{c,j}=v_{c,j}$ für alle $g$. Ein von null verschiedenes $v_{c,j}$ spannt eine invariante Gerade auf, die
die triviale Darstellung trägt; dann ist $\rho_c$ als irreduzible Darstellung trivial und $c=c_0$. Daher sind die
invarianten Vektoren die Vektoren des Sektors $c_0$, auf dem $\pi(g)$ als $1$ wirkt: das Bild von $P_{c_0}$, das
folglich $P_0$ ist.
(b) $q_A(\pi(g)x)=\langle x,\pi(g)^*A\pi(g)x\rangle$, und eine Matrix $D$ mit $\langle x,Dx\rangle=0$ für alle $x$
verschwindet (Polarisation, Beweis von Korollar 1.4). Also ist $q_A$ genau dann invariant, wenn $\pi(g)^*A\pi(g)=A$
für alle $g$, d. h. $A\in\pi(G)'$. Nach Proposition 12.6(a) sind die selbstadjungierten Elemente von $\pi(G)'$ die
$\Phi(w)$ mit hermiteschen $w_c$, ein reeller Raum der Dimension $\sum_cm_c^2$.
(c) Ein quadratisches Polynom $p(x)=x^TBx$ ist durch seine symmetrische Matrix $B$ bestimmt. Man setze
$(g\cdot p)(x)=p(\pi(g)^{-1}x)$; dann ist $g\cdot(h\cdot p)=(gh)\cdot p$, und $g\cdot p$ hat die Matrix
$\sigma(g)(B)=A^TBA$, $A=\pi(g)^{-1}$. Das Mittel $M=|G|^{-1}\sum_g\sigma(g)$ erfüllt $\sigma(h)M=M$ und $MB=B$ für
invariantes $B$, ist also (wie im Beweis von Satz 12.5) idempotent mit den invarianten Polynomen als Bild, und deren
Dimension ist $\operatorname{Tr}M=|G|^{-1}\sum_g\operatorname{Tr}\sigma(g)$. $A$ ist unitär, also $A=R\Lambda R^{-1}$
mit $\Lambda=\operatorname{diag}(\lambda)$ (Satz 1.3). Die symmetrischen Matrizen
$B_{ij}=(R^{-1})^T(e_ie_j^T+e_je_i^T)R^{-1}$, $i\le j$, bilden eine Basis, und
$A^TB_{ij}A=(R^{-1})^T\Lambda(e_ie_j^T+e_je_i^T)\Lambda R^{-1}=\lambda_i\lambda_jB_{ij}$. Also ist
$\operatorname{Tr}\sigma(g)=\sum_{i\le j}\lambda_i\lambda_j=\frac12\big((\operatorname{Tr}A)^2+\operatorname{Tr}A^2\big)
=\frac12\big(\chi_\pi(g^{-1})^2+\chi_\pi(g^{-2})\big)$, und Summation über $g^{-1}$ statt $g$ ergibt die Formel. Für
eine Permutationsdarstellung permutiert $G$ die Monome $x_ax_b$, also ist ein Polynom genau dann invariant, wenn seine
Koeffizienten auf den Bahnen der Monome konstant sind. $\square$

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

Auf $\mathbb C^4$ sind es $\sum_ix_i^2$ und $\sum_{i<j}x_ix_j$. Alle Grade zählt die Molien-Formel
$\sum_d\dim\mathbb C[x]^G_dt^d=|G|^{-1}\sum_g\det(1-t\pi(g))^{-1}$, deren Koeffizient von $t^2$ gerade (c) ist (ohne
Beweis angegeben; B. Sturmfels, *Algorithms in Invariant Theory*, Springer 2008, Theorem 2.2.1).

> **Physikalische Bedeutung.** In der Landau-Theorie der Phasenübergänge ist die freie Energie ein invariantes
> Polynom in einem reellen Ordnungsparameter $x$, auf den $G$ durch reelle orthogonale Matrizen wirkt. Ihr
> quadratischer Anteil ist $x^TBx$ mit $B$ in der Kommutante (da $\pi(g)^T=\pi(g)^{-1}$), also nach
> Proposition 12.6(a) blockdiagonal bezüglich der Sektoren: Die symmetrische Phase wird in dem Sektor instabil, in dem
> zuerst ein Eigenwert von $B$ das Vorzeichen wechselt, und (generisch) bestimmt die irreduzible Darstellung dieses
> Sektors die Symmetrie der geordneten Phase.

## 12.5 Womit man vergleichen sollte

Äquivariante Schichten für spezielle Gruppen stellen ausgereifte Bibliotheken bereit. **escnn** (Cesa, Lang und
Weiler 2022, Nachfolger von e2cnn) baut steuerbare (*steerable*) Faltungen für Untergruppen von $O(2)$ und $O(3)$ aus
irreduziblen Darstellungen und harmonischen Basen von Faltungskernen; **e3nn** (Geiger und Smidt 2022) baut
$E(3)$-äquivariante Netze aus Kugelflächenfunktionen und Clebsch–Gordan-Tensorprodukten. **EMLP** (Finzi, Welling und
Wilson 2021) behandelt beliebige Matrixgruppen: Es fordert $\pi_2(g)T=T\pi_1(g)$ für diskrete Erzeuger und
$d\pi_2(X)T=Td\pi_1(X)$ für Erzeuger der Lie-Algebra und berechnet eine Basis der Lösungen als numerischen Nullraum –
das lineare Problem von `dc.commutant_basis`. Gruppenfaltungen (Cohen und Welling 2016) sind der Fall der regulären
Darstellung, und Kondor und Trivedi (2018) zeigen, dass für Merkmalsabbildungen auf homogenen Räumen einer kompakten
Gruppe äquivariante lineare Schichten verallgemeinerte Faltungen sind. Für Permutationen sind DeepSets (Zaheer et al.
2017) und die Basen von Maron et al. (2019) Korollar 12.2.

Die Bibliothek bietet keine Faltungen auf Gittern, keine steuerbaren Kerne auf kontinuierlichen Räumen, keine
Zerlegungen über $\mathbb R$ (reeller, komplexer und quaternionischer Typ; Aufgabe 2 behandelt reelle Matrizen) und
kein großes $n$: $V$ ist eine dichte $n\times n$-Matrix, und die Zerlegung kostet $O(|G|n^3)$ pro Zufallsstichprobe
für eine durch all ihre Elemente gegebene Gruppe (Erzeuger von Lie-Algebren laufen über `dc.from_generators`;
Abschnitt 4.3). Was sie hinzufügt, ist die **Struktur** des Lösungsraums: Eine Basis der Kommutante ist eine Liste von
Matrizen, während die Zerlegung die Kommutante mit der $*$-Algebra $\bigoplus_cM_{m_c}$ identifiziert. Produkte,
Adjungierte, Unitarität, Positivität und der Funktionalkalkül wirken dann blockweise auf kleinen Matrizen, mit
Gradienten (Proposition 12.6), und `to_commutant`/`from_commutant` übersetzen zwischen den beiden Bildern.

## Aufgaben

1. (**DeepSets mit Kanälen.**) $S_n$ ($n\ge2$) wirkt auf $\mathbb C^n\otimes\mathbb C^d$ durch
   $P_\sigma\otimes1_d$. Man zeige, dass die äquivarianten Abbildungen von $\mathbb C^n\otimes\mathbb C^{d_1}$ nach
   $\mathbb C^n\otimes\mathbb C^{d_2}$ genau die $1\otimes A+J\otimes B$ mit $A,B\in M_{d_2\times d_1}$ sind, und
   prüfe dies mit `EquivariantMap`.
2. (**Reelle Schichten.**) Die Darstellungen $\pi_1,\pi_2$ mögen reelle Matrizen haben (z. B.
   Permutationsdarstellungen). Man zeige, dass die reellen äquivarianten Matrizen einen reellen Vektorraum der
   Dimension $\dim_{\mathbb C}\operatorname{Hom}_G(\pi_1,\pi_2)$ bilden und dass $T\mapsto\operatorname{Re}T$ den Raum
   $\operatorname{Hom}_G(\pi_1,\pi_2)$ auf ihn abbildet.
3. (**Nichtlinearitäten.**) (a) Man zeige, dass die koordinatenweise Anwendung einer Funktion
   $f:\mathbb C\to\mathbb C$ mit jeder Permutationsdarstellung vertauscht. (b) Man zeige, dass sie für
   $f(t)=\max(t,0)$ nicht mit der Drehung von $\mathbb R^2$ um $\pi/2$ vertauscht. (c) Man zeige, dass
   $x\mapsto\varphi(\|x\|)x$ mit jeder unitären Darstellung vertauscht.
4. (**Wie viele Stichproben legen eine Schicht fest?**) Seien $T=\Phi(w)\in\pi(G)'$ und
   $x_1,\dots,x_N\in\mathbb C^n$. Man schreibe $V^*x_i$ im Sektor $c$ als $\sum_{a,j}(X_{i,c})_{aj}e_a\otimes e_j$ mit
   $X_{i,c}\in M_{k_c\times m_c}$, und sei $S_c\in M_{Nk_c\times m_c}$ die Matrix, die durch Übereinanderstapeln von
   $X_{1,c},\dots,X_{N,c}$ entsteht. Man zeige, dass die äquivarianten $T$ mit $Tx_1=\dots=Tx_N=0$ einen Raum der
   Dimension $\sum_cm_c(m_c-\operatorname{rank}S_c)$ bilden; für Gaußsche $x_i$ ist dies fast sicher
   $\sum_cm_c\max(0,m_c-Nk_c)$. Für die Kantendarstellung legen zwei Stichproben eine äquivariante Schicht fest,
   während eine beliebige $16\times16$-Matrix sechzehn benötigt.

<details>
<summary>Lösungen</summary>

1. $P_\sigma J=J=JP_\sigma$, also sind diese Abbildungen äquivariant. Sie sind linear unabhängig: Aus
   $1\otimes A+J\otimes B=0$ folgt mit den Nichtdiagonalblöcken $B=0$ und dann $A=0$. Nach Korollar 12.2(a) und
   Satz 2.6 ist $P_\sigma\otimes1_d\cong(\mathrm{triv}\otimes1_d)\oplus(\mathrm{std}\otimes1_d)$, also haben beide
   irreduziblen Darstellungen die Vielfachheit $d$, und Satz 12.4 liefert $\dim\operatorname{Hom}_G=2d_1d_2$, die
   Dimension des Raums der $1\otimes A+J\otimes B$.

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

2. Ist $T$ äquivariant, so auch $\bar T$, da $\overline{\pi_i(g)}=\pi_i(g)$. Daher sind
   $\operatorname{Re}T=(T+\bar T)/2$ und $\operatorname{Im}T=(T-\bar T)/2i$ reell und äquivariant, also ist die
   Abbildung $T\mapsto\operatorname{Re}T$ surjektiv auf die reellen äquivarianten Matrizen (auf ihnen ist sie die
   Identität), und $\operatorname{Hom}_G(\pi_1,\pi_2)=H_{\mathbb R}\oplus iH_{\mathbb R}$, wobei $H_{\mathbb R}$ die
   reellen äquivarianten Matrizen bezeichnet. Eine reelle Basis von $H_{\mathbb R}$ ist dann eine komplexe Basis von
   $\operatorname{Hom}_G$: Sie spannt den Raum auf, und eine komplexe Relation $\sum_j(\alpha_j+i\beta_j)T_j=0$ mit
   reellen $T_j$ zerfällt in Real- und Imaginärteil. In Lemma 12.1 bilden die Bahnindikatoren eine solche reelle
   Basis. (Die Blöcke $y_c$ aus Satz 12.4 sind selbst für reelle Daten komplex: $V$ ist im Allgemeinen komplex, und
   die Parametrisierung über $\operatorname{Re}T$ verwendet doppelt so viele reelle Parameter wie nötig.)

3. (a) $(\pi(g)f(x))_{gy}=f(x)_y=f(x_y)=f((\pi(g)x)_{gy})$. (b) Für $R(x_1,x_2)=(-x_2,x_1)$ und $x=(1,-1)$ ist
   $f(Rx)=f(1,1)=(1,1)$, aber $Rf(x)=R(1,0)=(0,1)$. (c) $\|\pi(g)x\|=\|x\|$, also
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

4. $(1\otimes w_c)(e_a\otimes e_j)=\sum_le_a\otimes(w_c)_{lj}e_l$, also hat $V^*Tx_i$ im Sektor $c$ die
   Koeffizientenmatrix $X_{i,c}w_c^T$. Daher gilt $Tx_i=0$ für alle $i$ genau dann, wenn $S_cw_c^T=0$ für alle $c$,
   d. h. genau dann, wenn jede Spalte von $w_c^T$ in $\ker S_c$ liegt, einem Raum der Dimension
   $m_c-\operatorname{rank}S_c$: Die $w_c$ mit dieser Eigenschaft bilden einen Raum der Dimension
   $m_c(m_c-\operatorname{rank}S_c)$, und $\Phi$ ist injektiv. Für Gaußsche $x_i$ ist $V^*x_i$ Gaußsch ($V$ ist
   unitär), also hat $S_c$ unabhängige Gaußsche Einträge und fast sicher den Rang $\min(Nk_c,m_c)$: Für die linke
   obere quadratische Teilmatrix dieser Größe ist $|\det|^2$ ein reelles Polynom in den Real- und Imaginärteilen der
   Einträge, das nicht identisch verschwindet (es ist $1$, wenn die Teilmatrix die Einheitsmatrix ist), und ein solches
   Polynom verschwindet nur auf einer Lebesgue-Nullmenge (bewiesen in Abschnitt 4.3). Für die Kantendarstellung,
   $(k,m)=(1,2),(2,1),(3,1),(3,3)$, trägt nur der Sektor $(1,2)$ bei: $2\cdot(2-1)=2$ für $N=1$ und $0$ für $N=2$.

```python
basis = D_edge.basis('commutant')                                # (15, 16, 16), Hilbert-Schmidt orthonormal
for N, lost in ((1, 2), (2, 0)):
    Xs = torch.randn(16, N, dtype=dt)
    M = torch.stack([(B @ Xs).reshape(-1) for B in basis], dim=1)               # w -> (T x_1, ..., T x_N)
    assert torch.linalg.matrix_rank(M).item() == 15 - lost
```

</details>

**Wie es weitergeht.** Dies ist das letzte Kapitel des Kurses. [Kapitel 4](04_symmetries.md) entwickelt die hier
verwendete Darstellungstheorie: das Lemma von Schur, Charaktere und den randomisierten Algorithmus hinter
`dc.from_group` (Abschnitte 4.1–4.3). [Kapitel 9](09_spin_chains.md) wendet dieselben Sektoren auf Spinketten an, und das
Haar-Maß auf $\prod_cU(m_c)$ aus [Kapitel 5](05_random_matrices.md) liefert die natürliche zufällige Initialisierung
einer äquivarianten unitären Schicht. Zur Mathematik siehe J.-P. Serre, *Linear Representations of Finite Groups*
(Springer 1977) und W. Fulton und J. Harris, *Representation Theory* (Springer 1991); zum äquivarianten Lernen
M. Bronstein, J. Bruna, T. Cohen und P. Veličković, *Geometric Deep Learning* (2021), und M. Weiler et al.,
*Equivariant and Coordinate Independent Convolutional Networks* (2023); zum Zerlegungsalgorithmus die Literaturangaben
in `torch_vn_algebra/decompose.py`.
