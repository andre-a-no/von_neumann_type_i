[English](../10_quantum_computing.md) · [Español](../es/10_quantum_computing.md) · [Français](../fr/10_quantum_computing.md) · **Deutsch** · [中文](../zh/10_quantum_computing.md) · [日本語](../ja/10_quantum_computing.md) · [Русский](../ru/10_quantum_computing.md)

# Tutorial 10. Quantencomputing

Dieses Tutorial begleitet [Kapitel 10](../../lecture_notes/de/10_quantum_computing.md) der Vorlesungsnotizen, das Rauschen auf einem Quantencomputer
durch die $*$-Algebra seiner Kraus-Operatoren beschreibt. Wir berechnen Fixpunktalgebren und rauschfreie Teilsysteme
([Satz 10.2](../../lecture_notes/de/10_quantum_computing.md#r-10-2)), prüfen den Satz von Knill–Laflamme ([Satz 10.5](../../lecture_notes/de/10_quantum_computing.md#r-10-5)) und die Korrektur der Bezeichnung eines
Sektors ([Proposition 10.6](../../lecture_notes/de/10_quantum_computing.md#r-10-6)) am Wiederholungscode und simulieren einen Ansatz mit XY-Mixer, Amplitudendämpfung
und Postselektion in der Algebra $\bigoplus_wM_{\binom nw}$ der Hamming-Gewichts-Sektoren
([Propositionen 10.8](../../lecture_notes/de/10_quantum_computing.md#r-10-8)–[10.10](../../lecture_notes/de/10_quantum_computing.md#r-10-10)).

<a id="s-10-1"></a>

## 10.1 Rauschalgebren, rauschfreie Teilsysteme und Fixpunkte

*Vorlesungsnotizen:* [Abschnitt 10.1](../../lecture_notes/de/10_quantum_computing.md#s-10-1); [Satz 10.2](../../lecture_notes/de/10_quantum_computing.md#r-10-2), [Korollar 10.3](../../lecture_notes/de/10_quantum_computing.md#r-10-3), [Proposition 10.4](../../lecture_notes/de/10_quantum_computing.md#r-10-4).

Nach [Vorlesungsnotizen, Korollar 10.3](../../lecture_notes/de/10_quantum_computing.md#r-10-3) hat die kollektive Dephasierung $\Phi(\rho)=\tfrac12(U\rho U^*+U^*\rho U)$,
$U=e^{-i\theta Z_{\mathrm{tot}}}$, auf $n$ Qubits die Wechselwirkungsalgebra $\mathcal A=\mathrm{span}\{P_0,\dots,P_n\}$
($P_w$ die Projektion auf den Spann $H_w$ der Basisvektoren vom Hamming-Gewicht $w$), und
$\mathrm{Fix}(\Phi)=\mathcal A'=\bigoplus_wB(H_w)$ ([Satz 10.2](../../lecture_notes/de/10_quantum_computing.md#r-10-2)).

Kommutanten sind Nullräume von $T\mapsto(AT-TA)_A$ (Kapitel 2), Fixpunkte der von
$\sum_iE_i\otimes\bar E_i-1$ (der Matrix von $\Phi-\mathrm{id}$ auf zeilenweise vektorisierten
Matrizen). Für $n=4$ haben die drei Räume des Korollars die Dimensionen $5$, $70=1+16+36+16+1$ und
$70$:

```python
import itertools
import math
import torch

torch.manual_seed(0)
dt = torch.complex128
I2, X = torch.eye(2, dtype=dt), torch.tensor([[0, 1], [1, 0]], dtype=dt)
Y, Z = torch.tensor([[0, -1j], [1j, 0]], dtype=dt), torch.tensor([[1, 0], [0, -1]], dtype=dt)

def on(op, j, n):                             # op on qubit j of n (j = 0 is the leftmost factor)
    out = torch.ones(1, 1, dtype=dt)
    for i in range(n):
        out = torch.kron(out, op if i == j else I2)
    return out

def null_space(L, tol=1e-9):                  # orthonormal columns spanning {x : L x = 0}
    w, V = torch.linalg.eigh(L.mH @ L)
    return V[:, w < tol]

def commutant_basis(ops):
    d, ops = ops[0].shape[0], [A.resolve_conj().contiguous() for A in ops]
    eye = torch.eye(d, dtype=dt)
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    return [v.reshape(d, d) for v in null_space(L).T]

def fixed_points(kraus):                      # vec(E X E^*) = (E (x) conj E) vec(X)
    d = kraus[0].shape[0]
    S = sum(torch.kron(E.resolve_conj(), E.conj().resolve_conj()) for E in kraus)
    return [v.reshape(d, d) for v in null_space(S - torch.eye(d * d, dtype=dt)).T]

n, theta = 4, 0.3
U = torch.linalg.matrix_exp(-1j * theta * sum(on(Z, j, n) for j in range(n)))
kraus = [U / math.sqrt(2), U.mH / math.sqrt(2)]                  # unital: U U^* = U^* U = 1
comm = commutant_basis(kraus + [E.mH for E in kraus])             # A'
A_int, fix = commutant_basis(comm), fixed_points(kraus)           # A = A'' (Theorem 2.2), Fix(Phi)
print(len(A_int), len(comm), len(fix))
assert (len(A_int), len(comm), len(fix)) == (5, 70, 70)
span = torch.linalg.svdvals(torch.stack([F.reshape(-1) for F in comm + fix]))
assert int((span > 1e-9).sum()) == 70                             # Fix(Phi) = A' (Theorem 10.2)
```

$\bigoplus_wM_{\binom nw}$ ist die Algebra einer `SpinChain` (Kapitel 9), deren Sektoren die Anzahlen
der Spins nach oben sind, d. h. die Hamming-Gewichte. Dort sind die Kraus-Operatoren zentral,
$U=\sum_we^{-i\theta(n-2w)}P_w$, und der Kanal ist die Identität: Die Bibliothek speichert den
rauschfreien Teil.

```python
from torch_vn_algebra import SpinChain, TypeIAlgebra, DensityMatrix, channels, states, dynamics

M = SpinChain(n).algebra                                          # (+)_w M_binom(4, w)
print(M.k_factors)
assert sum(k * k for k in M.k_factors) == len(fix)
phases = [torch.exp(torch.tensor(-1j * theta * (n - 2 * w), dtype=dt)) for w in range(n + 1)]
Phi_M = channels.kraus_channel([M.central([z / math.sqrt(2) for z in phases]),
                                M.central([z.conj() / math.sqrt(2) for z in phases])])
rho = states.random_density_matrix(M, batch_size=3)
assert (Phi_M(rho).matrix - rho.matrix).abs().max() < 1e-12      # no state of M is disturbed
```

Kollektive Drehungen $e^{-i\theta\,\mathbf n\cdot\mathbf S}$, $S_\alpha=\tfrac12\sum_j\sigma_{\alpha,j}$, liegen in der von
$1,S_x,S_y,S_z$ erzeugten $*$-Algebra $\mathcal B$. Für $n=4$ zeigt [Abschnitt 10.1](../../lecture_notes/de/10_quantum_computing.md#s-10-1) der Vorlesungsnotizen $\dim\mathcal B=35$,
$\dim\mathcal B'=14$ und $\dim K=2$ für $K=\ker S_x\cap\ker S_y\cap\ker S_z$; der Code findet dieselben Zahlen und
prüft [Proposition 10.4](../../lecture_notes/de/10_quantum_computing.md#r-10-4): Ein Kanal mit zufälligen kollektiven Drehungen als Kraus-Operatoren lässt einen
Zustand mit Träger in $K$ fest und bewegt einen generischen.

```python
S = [sum(on(P, j, n) for j in range(n)) / 2 for P in (X, Y, Z)]
comm_S = commutant_basis(S)
print(len(comm_S), len(commutant_basis(comm_S)))
assert len(comm_S) == 14 and len(commutant_basis(comm_S)) == 35
K = null_space(torch.cat(S))                                      # columns: a basis of K
assert K.shape[1] == 2
axes = torch.nn.functional.normalize(torch.randn(3, 3, dtype=torch.float64), dim=1)
rotations = [torch.linalg.matrix_exp(-1j * t * sum(c * Sa for c, Sa in zip(ax, S))) / math.sqrt(3)
             for t, ax in zip(torch.rand(3, dtype=torch.float64) * math.pi, axes)]
apply = lambda kr, r: sum(E @ r @ E.mH for E in kr)
for G in (K @ torch.randn(2, 2, dtype=dt), torch.randn(16, 16, dtype=dt)):  # supported in K / generic
    r = G @ G.mH / (G @ G.mH).trace()
    assert ((apply(rotations, r) - r).abs().max() < 1e-12) == (G.shape[1] == 2)
```

<a id="s-10-2"></a>

## 10.2 Fehlerkorrektur: der Satz von Knill–Laflamme

*Vorlesungsnotizen:* [Abschnitt 10.2](../../lecture_notes/de/10_quantum_computing.md#s-10-2); [Satz 10.5](../../lecture_notes/de/10_quantum_computing.md#r-10-5).

Das Bitflip-Beispiel aus [Abschnitt 10.2](../../lecture_notes/de/10_quantum_computing.md#s-10-2) der Vorlesungsnotizen: Für $\mathcal C=\mathrm{span}\{|000\rangle,|111\rangle\}$ und
„höchstens ein Bitflip“ gilt die Knill–Laflamme-Bedingung aus [Satz 10.5](../../lecture_notes/de/10_quantum_computing.md#r-10-5) mit
$\Lambda=\mathrm{diag}(p_0,\dots,p_3)$, und die Wiederherstellung $\mathcal R(\sigma)=\sum_{k=0}^3PX_k\sigma X_kP$ stellt jeden
Zustand des Codes wieder her. Kippt jedes Qubit unabhängig mit Wahrscheinlichkeit $p$, so ist die Bedingung verletzt,
und die Fidelität von $|0_L\rangle$ nach der Wiederherstellung ist $1-3p^2+2p^3=0.972$ für $p=0.1$. Die Hilfsfunktion
`kl_matrix` liefert $\Lambda$ oder `None`, wenn die Bedingung verletzt ist. In der Bibliothek ist $M_8$ eine Algebra mit
einem Block, und `R @ Phi` ist die Komposition $\mathcal R\circ\Phi$:

```python
M8 = TypeIAlgebra([8], [8], precision='double')
op8 = lambda T: M8.operator(T.reshape(1, 1, 8, 8))
state8 = lambda v: DensityMatrix(M8, matrix=torch.outer(v, v.conj()).reshape(1, 1, 8, 8))
ket = lambda bits: torch.eye(8, dtype=dt)[int(bits, 2)]
def kl_matrix(kraus, P, tol=1e-12):
    """(lambda_ij) with P E_i^* E_j P = lambda_ij P, or None if the Knill-Laflamme condition fails."""
    lam = torch.tensor([[(P @ a.mH @ b @ P).trace() / P.trace() for b in kraus] for a in kraus])
    ok = all((P @ a.mH @ b @ P - lam[i, j] * P).abs().max() < tol
             for i, a in enumerate(kraus) for j, b in enumerate(kraus))
    return lam if ok else None
def independent(op, prob):                    # op on each qubit independently with probability prob
    return [math.sqrt(prob ** sum(s) * (1 - prob) ** (3 - sum(s)))
            * torch.linalg.multi_dot([on(op, j, 3) if s[j] else torch.eye(8, dtype=dt) for j in range(3)])
            for s in itertools.product([0, 1], repeat=3)]

P = torch.outer(ket('000'), ket('000')) + torch.outer(ket('111'), ket('111'))
flips = [torch.eye(8, dtype=dt)] + [on(X, j, 3) for j in range(3)]
E_one = [math.sqrt(pk) * Xk for pk, Xk in zip((0.7, 0.1, 0.1, 0.1), flips)]
assert torch.allclose(kl_matrix(E_one, P), torch.diag(torch.tensor([0.7, 0.1, 0.1, 0.1], dtype=dt)))
Phi_one = channels.kraus_channel([op8(E) for E in E_one])
R = channels.kraus_channel([op8(P @ Xk) for Xk in flips])        # syndrome measurement + correction
assert R.is_trace_preserving()
a, b = torch.nn.functional.normalize(torch.randn(2, dtype=dt), dim=0)
rho_L = state8(a * ket('000') + b * ket('111'))
assert abs((R @ Phi_one)(rho_L).fidelity(rho_L).item() - 1) < 1e-12       # Theorem 10.5
p = 0.1                                                           # independent flips: KL fails
E_ind = independent(X, p)
assert kl_matrix(E_ind, P) is None
zero_L = state8(ket('000'))
F = (R @ channels.kraus_channel([op8(E) for E in E_ind]))(zero_L).fidelity(zero_L).item()
print(f"fidelity after recovery {F:.4f}, exact {1 - 3 * p**2 + 2 * p**3:.4f}")
assert abs(F - (1 - 3 * p**2 + 2 * p**3)) < 1e-12
```

<a id="s-10-3"></a>

## 10.3 Fehlerkorrektur mit Operatoralgebren

*Vorlesungsnotizen:* [Abschnitt 10.3](../../lecture_notes/de/10_quantum_computing.md#s-10-3); [Proposition 10.6](../../lecture_notes/de/10_quantum_computing.md#r-10-6).

Das Beispiel aus [Abschnitt 10.3](../../lecture_notes/de/10_quantum_computing.md#s-10-3) der Vorlesungsnotizen: der Wiederholungscode unter unabhängigen Phasenflips
(Wahrscheinlichkeit $q$ pro Qubit), gefolgt von höchstens einem Bitflip. Das Qubit ist nicht geschützt (die
Knill–Laflamme-Bedingung ist verletzt), wohl aber das klassische Bit: Bedingung 2 aus [Proposition 10.6](../../lecture_notes/de/10_quantum_computing.md#r-10-6) gilt für
$P_0=|000\rangle\langle000|$, $P_1=|111\rangle\langle111|$, und der Mehrheitsentscheid $N_0=\sum_{|b|\le1}|b\rangle\langle b|$,
$N_1=1-N_0$ liest das Bit aus. Nach der Wiederherstellung aus Abschnitt 10.2 ist die Fidelität von $|{+_L}\rangle$ gleich
$(1+(1-2q)^3)/2=0.8645$ für $q=0.05$.

```python
q = 0.05
Phi2 = Phi_one @ channels.kraus_channel([op8(E) for E in independent(Z, q)])   # dephasing, then a flip
E2 = [Phi2.kraus[0, i, 0] for i in range(Phi2.kraus_rank)]        # its 32 Kraus operators
P0, P1 = torch.outer(ket('000'), ket('000')), torch.outer(ket('111'), ket('111'))
assert kl_matrix(E2, P) is None                                   # the qubit is not protected
assert max((P0 @ e.mH @ f @ P1).abs().max().item() for e in E2 for f in E2) == 0      # condition 2
N0 = torch.diag(torch.tensor([float(bin(i).count('1') <= 1) for i in range(8)], dtype=dt))
for N, Pd in ((N0, P0), (torch.eye(8, dtype=dt) - N0, P1)):       # P Phi^*(N_d) P = P_d
    assert (P @ sum(e.mH @ N @ e for e in E2) @ P - Pd).abs().max() < 1e-12
votes = states.born_probabilities(Phi2(rho_L), [op8(N0), op8(torch.eye(8, dtype=dt) - N0)])[0]
assert torch.allclose(votes, torch.stack([abs(a) ** 2, abs(b) ** 2]))      # the classical bit survives
plus_L = state8((ket('000') + ket('111')) / math.sqrt(2))
F_plus = (R @ Phi2)(plus_L).fidelity(plus_L).item()
print(f"fidelity of |+_L> {F_plus:.4f}, exact {(1 + (1 - 2 * q)**3) / 2:.4f}")
assert abs(F_plus - (1 + (1 - 2 * q) ** 3) / 2) < 1e-12
```

<a id="s-10-4"></a>

## 10.4 Mixer mit Nebenbedingungen

*Vorlesungsnotizen:* [Abschnitt 10.4](../../lecture_notes/de/10_quantum_computing.md#s-10-4); [Proposition 10.8](../../lecture_notes/de/10_quantum_computing.md#r-10-8), [Proposition 10.9](../../lecture_notes/de/10_quantum_computing.md#r-10-9), [Proposition 10.10](../../lecture_notes/de/10_quantum_computing.md#r-10-10).

Der Quantum Alternating Operator Ansatz wechselt $e^{-i\gamma H_C}$, $H_C$ ein diagonaler Kosten-Hamiltonoperator,
und $e^{-i\beta H_M}$ ab, $H_M=\tfrac14\sum_{(i,j)}(X_iX_j+Y_iY_j)$ der XY-Mixer auf einem Ring. Nach
[Vorlesungsnotizen, Proposition 10.8](../../lecture_notes/de/10_quantum_computing.md#r-10-8) liegen sie alle in $\mathcal M=\{P_0,\dots,P_n\}'\cong\bigoplus_wM_{\binom nw}$,
also verlässt der Ansatz den Hamming-Gewichts-Sektor seines Anfangszustands nie.

In der Bibliothek ist $\mathcal M$ die Algebra `SpinChain(n).algebra` (Bit $1$ = Spin nach oben,
also $S^z_i=-\tfrac12Z_i$), und $H_M$ ist ihr XXZ-Hamiltonoperator
$J\sum_{(i,j)}\big(\tfrac12(S^+_iS^-_j+S^-_iS^+_j)+\Delta S^z_iS^z_j\big)$ mit $J=1$, $\Delta=0$, da
$S^+_iS^-_j+S^-_iS^+_j=\tfrac12(X_iX_j+Y_iY_j)$; wir prüfen das an der dichten $64\times64$-Matrix. Für
$H_C$ nehmen wir die Graphbisektion der Dreiecke $\{0,1,2\}$, $\{3,4,5\}$, verbunden durch die Kanten
$(2,3)$, $(5,0)$: $H_C=\sum_{(i,j)\in E}\tfrac12(1-Z_iZ_j)$, $Z_iZ_j=4S^z_iS^z_j$, zählt die geschnittenen
Kanten. Die beste Bisektion schneidet $2$ Kanten (die beiden Dreiecke; jede andere Menge von drei Ecken
teilt beide Dreiecke, und ein geteiltes Dreieck hat zwei geschnittene Kanten); wir starten bei $|101010\rangle$, das $6$ schneidet.
Schließlich wird der ganze Ansatz mit dichten $64\times64$-Matrizen unabhängig von der Bibliothek
nachgerechnet.

```python
n, w = 6, 3
chain = SpinChain(n, 'periodic')
M = chain.algebra                                                 # (+)_w M_binom(6, w), padded to 20
H_M = chain.xxz(J=1.0, Delta=0.0)
H_dense = sum(on(X, i, n) @ on(X, j, n) + on(Y, i, n) @ on(Y, j, n) for i, j in chain.bonds) / 4
spec = torch.sort(torch.cat([e[0] for e in H_M.eigenvalues()]))[0]
assert torch.allclose(spec, torch.linalg.eigvalsh(H_dense))       # same spectrum, sector by sector
edges = chain.bonds + [(0, 2), (3, 5)]                            # two triangles joined by (2,3), (5,0)
H_C = M.zero()
for i, j in edges:
    H_C = H_C + (M.identity() - chain.sz(i) @ chain.sz(j) * 4.0) * 0.5
layers = [dynamics.propagator(H_M, beta) @ dynamics.propagator(H_C, gamma)
          for gamma, beta in ((0.6, 1.0), (1.4, 1.4))]
rho0 = ideal = chain.basis_state('101010')
for U_l in layers:
    ideal = channels.unitary_channel(U_l)(ideal)
assert abs(ideal.sector_probabilities()[0, w].item() - 1) < 1e-12           # Proposition 10.8

HC_dense = sum(torch.eye(64, dtype=dt) - on(Z, i, n) @ on(Z, j, n) for i, j in edges) / 2
psi = torch.eye(64, dtype=dt)[int('101010', 2)]                   # the same circuit on 6 full qubits
for gamma, beta in ((0.6, 1.0), (1.4, 1.4)):
    psi = torch.linalg.matrix_exp(-1j * beta * H_dense) @ torch.linalg.matrix_exp(-1j * gamma * HC_dense) @ psi
weight = torch.tensor([bin(b).count('1') for b in range(64)])
assert abs(psi[weight == w].norm().item() - 1) < 1e-12           # Proposition 10.8, dense
assert abs((psi.conj() @ HC_dense @ psi).real.item() - states.expectation(ideal, H_C).real.item()) < 1e-10
```

Die **Amplitudendämpfung** des Qubits $i$ ($|1\rangle\to|0\rangle$ mit Wahrscheinlichkeit $g$) hat die
Kraus-Operatoren $K_0^{(i)}=|0\rangle\langle0|_i+\sqrt{1-g}\,|1\rangle\langle1|_i$ und
$K_1^{(i)}=\sqrt g\,|0\rangle\langle1|_i$, der das Gewicht senkt: eine Abbildung zwischen Sektoren,
`chain.site_amplitude_damping(i, g)` (ein `InterSectorChannel`).

Nach [Vorlesungsnotizen, Proposition 10.9](../../lecture_notes/de/10_quantum_computing.md#r-10-9) wird, wenn auf jede Schicht des Schaltkreises die Dämpfung
aller Qubits folgt, das Gewicht $w$ des Anfangszustands nach $p$ Schichten noch mit Wahrscheinlichkeit $(1-g)^{wp}$
gefunden, und der postselektierte Zustand ist der rauschfreie.

Ein `InterSectorChannel` behält nur die Diagonalblöcke seiner Ausgabe; hier gibt es keine anderen
(der Beweis von Proposition 10.9 gilt wörtlich für die Dämpfung eines einzelnen Qubits, mit $s\in\{0,1\}$),
also geht nichts verloren. Mit $g=0.05$ und $p=2$ ist die Akzeptanz $0.95^6=0.7351$:

```python
g = 0.05
damping = [chain.site_amplitude_damping(i, g) for i in range(n)]
noisy = rho0
for U_l in layers:
    noisy = channels.unitary_channel(U_l)(noisy)
    for D in damping:                                             # one qubit after the other
        noisy = D(noisy)
print(noisy.sector_probabilities()[0])                            # weight leaks to w = 2, 1, 0
P_w = M.central([1.0 if v == w else 0.0 for v in range(n + 1)])
post, accept = noisy.condition_on(P_w)                            # measure the weight, keep w = 3
cut = lambda r: states.expectation(r, H_C).real.item()
print(f"P(weight 3) {accept.item():.4f}, fidelity {noisy.fidelity(ideal).item():.4f}, "
      f"post-selected {post.fidelity(ideal).item():.4f}; "
      f"cut: ideal {cut(ideal):.4f}, noisy {cut(noisy):.4f}, post-selected {cut(post):.4f}")
assert abs(accept.item() - (1 - g) ** (w * len(layers))) < 1e-12              # Proposition 10.9
assert abs(post.fidelity(ideal).item() - 1) < 1e-10 and abs(cut(post) - cut(ideal)) < 1e-10
```

Numerische Beobachtung (die ausgegebenen Werte): Der unzulässige Anteil senkt den erwarteten Schnitt
von $4.3789$ auf $4.3123$, eine scheinbare Verbesserung, die einen Optimierer der Winkel in die Irre
führen würde; nach der Postselektion ist der Wert exakt (Proposition 10.9).

Nach [Vorlesungsnotizen, Proposition 10.10](../../lecture_notes/de/10_quantum_computing.md#r-10-10) behält unter der Lindblad-Gleichung mit Amplitudendämpfung
(Rate $\gamma_1$) und Dephasierung $L_k=S^z_k$ (Rate $\gamma_2$) ein in $H_w$ startender Zustand das Gewicht $w$ mit
Wahrscheinlichkeit $e^{-\gamma_1wt}$, und der postselektierte Zustand hängt nicht von $\gamma_1$ ab; für $\gamma_2=0$
ist die Postselektion perfekt. Der
Begleitartikel (`examples/paper/ex9_qaoa_xy_mixer.py`) lässt $|111000\rangle$ unter $H_M$ bis $T=2$
evolvieren; wir reproduzieren seine Zeile $\gamma_1=\gamma_2=0.1$: $P(w=3)=e^{-0.6}=0.5488$ (Proposition
10.10) sowie, als numerische Beobachtung des folgenden Codes, die Fidelität $0.4805$ und die postselektierte
Fidelität $0.8755$; nach Proposition 10.10 ist letztere dieselbe wie für $\gamma_1=0$.

```python
T = 2.0
jumps = [chain.lowering(i) for i in range(n)] + [chain.sz(i) for i in range(n)]
rho0 = chain.basis_state('111000')                                # the initial state of the paper
ideal_T = dynamics.von_neumann(rho0, H_M, [T])[0]
post_T = {}
for g1, g2 in ((0.1, 0.0), (0.0, 0.1), (0.1, 0.1)):          # the last one is printed
    rho_T = dynamics.lindblad_evolve(rho0, H_M, jumps, [0.0, T], rates=[g1] * n + [g2] * n, substeps=80)[-1]
    post_T[g1, g2], accept_T = rho_T.condition_on(P_w)
    assert abs(accept_T.item() - math.exp(-w * g1 * T)) < 1e-8                 # Proposition 10.10
print(f"P(weight 3) {accept_T.item():.4f}, fidelity {rho_T.fidelity(ideal_T).item():.4f}, "
      f"post-selected {post_T[0.1, 0.1].fidelity(ideal_T).item():.4f}")
assert abs(post_T[0.1, 0.0].fidelity(ideal_T).item() - 1) < 1e-6               # gamma_2 = 0: perfect
assert post_T[0.1, 0.1].trace_distance(post_T[0.0, 0.1]).item() < 1e-6       # independent of gamma_1
```

<a id="s-10-5"></a>

## 10.5 Womit man vergleichen sollte

Qiskit (`qiskit.quantum_info`, Rauschmodelle von Aer, `QAOAAnsatz` mit einem eigenen Mixer), Cirq,
PennyLane (`qml.qaoa.xy_mixer`, `default.mixed`) und QuTiP/`qutip-qip` simulieren all dies auf $n$
vollen Qubits und sind die richtigen Werkzeuge für Schaltkreise; stim sampelt Clifford-Schaltkreise
mit Pauli-artigem Rauschen auf Tausenden von Qubits (Amplitudendämpfung ist kein Pauli-Kanal: Pauli-Kanäle sind unital, die Amplitudendämpfung nicht). Uns
ist in diesen Paketen keine Funktion bekannt, die die Wechselwirkungsalgebra eines Kanals, ihren
Kommutanten oder seine ganze Fixpunktalgebra berechnet (QuTiPs `steadystate` liefert einen
stationären Zustand eines Lindblad-Operators); zu numerischen Blockzerlegungen von
Matrix-$*$-Algebren siehe Murota et al., *Japan J. Indust. Appl. Math.* 27 (2010). Die Bibliothek
ergänzt die Speicherung in $\bigoplus_wM_{\binom nw}$, Abbildungen zwischen Sektoren und Postselektion
als zentrale Projektion, gestapelt. Ein Sektor braucht $\binom nw^2$ Einträge statt $4^n$
($1.7\cdot10^8$ gegenüber $4.3\cdot10^9$ für $n=16$, $w=8$), aber die ganze Algebra wird auf
$(n+1)\binom n{\lfloor n/2\rfloor}^2$ Einträge ($2.8\cdot10^9$) aufgefüllt, eine bescheidene Ersparnis. Sie
ist kein Schaltkreissimulator und hat keinen Stabilisatorformalismus.

<a id="x-10"></a>

## Aufgaben

[Aufgaben 1–3](../../lecture_notes/de/10_quantum_computing.md#x-10) der Vorlesungsnotizen, numerisch: Für die Amplitudendämpfung eines Qubits sind
$\mathcal A'$, $\mathrm{Fix}(\Phi)$ und $\mathrm{Fix}(\Phi^*)$ eindimensional (Aufgabe 1); für den Kanal aus
[Korollar 10.3](../../lecture_notes/de/10_quantum_computing.md#r-10-3) mit $\theta=\pi/2$ und $n=4$ ist $\dim\mathrm{Fix}(\Phi)=128$ (Aufgabe 2); die Wiederherstellung
$\mathcal R$ des Bitflip-Codes korrigiert den kohärenten Fehler $e^{-i\phi X_1}$ (Aufgabe 3).

<details>
<summary>Lösungen</summary>

```python
E_ad = [torch.tensor([[1, 0], [0, math.sqrt(0.7)]], dtype=dt), torch.tensor([[0, math.sqrt(0.3)], [0, 0]], dtype=dt)]
fix_ad = fixed_points(E_ad)
assert len(fix_ad) == 1 and len(commutant_basis(E_ad + [E.mH for E in E_ad])) == 1
assert torch.allclose(fix_ad[0] / fix_ad[0][0, 0], torch.tensor([[1, 0], [0, 0]], dtype=dt))
assert len(fixed_points([E.mH for E in E_ad])) == 1                     # Fix(Phi^*) = C 1
```

```python
U2 = torch.linalg.matrix_exp(-1j * (math.pi / 2) * sum(on(Z, j, 4) for j in range(4)))
assert len(fixed_points([U2 / math.sqrt(2), U2.mH / math.sqrt(2)])) == 128
assert sum(math.comb(4, v) ** 2 for v in range(5)) == math.comb(8, 4) == 70
```

```python
for phi in (0.2, 1.0, 2.5):
    G1 = torch.linalg.matrix_exp(-1j * phi * on(X, 0, 3))
    out = (R @ channels.kraus_channel([op8(G1)]))(rho_L)
    assert (out.matrix - rho_L.matrix).abs().max() < 1e-12       # corrected exactly
```

</details>

**Weiter:** [Tutorial 11](11_optimisation.md): Optimierung über die unitären Gruppen einer Algebra und Spurungleichungen.
