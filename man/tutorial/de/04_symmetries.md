[English](../04_symmetries.md) · [Español](../es/04_symmetries.md) · [Français](../fr/04_symmetries.md) · **Deutsch** · [中文](../zh/04_symmetries.md) · [日本語](../ja/04_symmetries.md) · [Русский](../ru/04_symmetries.md)

# Tutorial 4. Symmetrien und Darstellungen

Dieses Tutorial begleitet [Kapitel 4](../../lecture_notes/de/04_symmetries.md) der Vorlesungsnotizen, das die Darstellungstheorie hinter den
Symmetriealgebren entwickelt, daraus einen randomisierten Algorithmus gewinnt, der die von einer beliebigen Liste von
Matrizen erzeugte von-Neumann-Algebra zerlegt ([Proposition 4.9](../../lecture_notes/de/04_symmetries.md#r-4-9)), ihn auf $U(1)$ und $SU(2)$ anwendet und
die Sektoren einer von reellen Matrizen erzeugten Algebra in reellen, komplexen und quaternionischen Typ einteilt
([Satz 4.20](../../lecture_notes/de/04_symmetries.md#r-4-20)). Wir berechnen Charaktere, wenden die Zerlegung auf Gruppen und auf Erzeuger an, bauen
Spinketten sektorweise mit `SpinChain` und `su2.SU2Chain` auf und berechnen die reelle Struktur einer Zerlegung.

<a id="s-4-2"></a>

## 4.2 Die isotypische Zerlegung und die Kommutante

*Vorlesungsnotizen:* [Abschnitt 4.2](../../lecture_notes/de/04_symmetries.md#s-4-2); [Satz 4.3](../../lecture_notes/de/04_symmetries.md#r-4-3), [Satz 4.4](../../lecture_notes/de/04_symmetries.md#r-4-4), [Korollar 4.5](../../lecture_notes/de/04_symmetries.md#r-4-5).

[Abschnitt 4.2](../../lecture_notes/de/04_symmetries.md#s-4-2) der Vorlesungsnotizen berechnet von Hand die Charaktere dreier irreduzibler
Darstellungen von $S_3$ – der trivialen, des Signums und der Standarddarstellung auf $(1,1,1)^\perp\subseteq\mathbb C^3$ – und
der Darstellung $\pi$, die die drei Faktoren von $(\mathbb C^2)^{\otimes3}$ permutiert, $\chi_\pi(\sigma)=2^{\#\text{cycles}(\sigma)}$
(Tabelle unten), und gewinnt aus [Korollar 4.5](../../lecture_notes/de/04_symmetries.md#r-4-5) die Vielfachheiten $m_{\rm triv}=4$, $m_{\rm sgn}=0$,
$m_{\rm std}=2$ und $\dim\pi(S_3)'=20$. Der Code berechnet die Tabelle neu: Er gewinnt $P_\sigma$ als Einschränkung von
$\pi(\sigma)$ auf den invarianten Spann von $|100\rangle,|010\rangle,|001\rangle$ (auf dem $\pi(\sigma)$ die Positionen der
einzigen $1$ permutiert, Lemma B.1(a)), und $\#\text{cycles}(\sigma)=\dim\ker(P_\sigma-1)$ (die auf den Zyklen konstanten Vektoren).

| | Identität ($1$) | Transpositionen ($3$) | $3$-Zyklen ($2$) |
|---|---|---|---|
| $\chi_{\rm triv}$ | $1$ | $1$ | $1$ |
| $\mathrm{sgn}$ | $1$ | $-1$ | $1$ |
| $\chi_{\rm std}$ | $2$ | $0$ | $-1$ |
| $\chi_\pi$ | $8$ | $4$ | $2$ |

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

Das ist die Zerlegung aus Abschnitt 2.4, nun gewonnen aus einer Tabelle von sechs Zahlen pro Darstellung.

<a id="s-4-3"></a>

## 4.3 Von den Erzeugern zur Algebra: ein randomisierter Algorithmus

*Vorlesungsnotizen:* [Abschnitt 4.3](../../lecture_notes/de/04_symmetries.md#s-4-3); [Lemma 4.6](../../lecture_notes/de/04_symmetries.md#r-4-6), [Proposition 4.7](../../lecture_notes/de/04_symmetries.md#r-4-7), [Lemma 4.8](../../lecture_notes/de/04_symmetries.md#r-4-8), [Proposition 4.9](../../lecture_notes/de/04_symmetries.md#r-4-9).

Der randomisierte Algorithmus der Vorlesungsnotizen ist `decompose.wedderburn`, mit den Schnittstellen `from_group` (alle
Elemente einer endlichen Gruppe; die Kommutante wird mit dem Twirl `decompose.twirl` abgetastet) und `from_generators`
(eine beliebige Liste von Matrizen; die Kommutante ist ein Kern, `decompose.commutant_basis`). Er zerlegt $\mathbb C^n$
in die Eigenräume eines zufälligen $K\in\mathcal A'_h$, gruppiert sie und richtet sie mit den Blöcken eines zufälligen
$Y\in\mathcal A'$ aus und gibt eine `Decomposition` mit dem unitären `V`, den $k_c$ (`k_factors`) und den $m_c$
(`multiplicities`) zurück. In Gleitkommaarithmetik behandelt er Eigenwerte, die näher als eine Toleranz beieinanderliegen,
als gleich und zerlegt die entstehenden Unterräume weiter, indem er neue zufällige Elemente der Kommutante auf sie
komprimiert ([Lemma 4.8(a)](../../lecture_notes/de/04_symmetries.md#r-4-8)).

`Decomposition.verify()` meldet die relativen Abstände der Erzeuger (und neuer Elemente der Kommutante) von
$\mathcal B$ (bzw. $\mathcal B'$) und vergleicht $\sum_cm_c^2$ mit $\dim\mathcal A'$ (für eine Gruppe aus den
Charakteren, Korollar 4.5(c); für Erzeuger aus dem Kern). Wir wenden es auf $S_3$ und drei Qubits an und schauen
dann hinein: die Eigenräume eines zufälligen $K$ und die Verflechtungsoperatoren eines zufälligen $Y$.

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
s = torch.linalg.svdvals(Wa.mH @ Y @ Wb)                               # Lemma 4.8(b): a multiple of a unitary
assert s[1] > 1e-3 and abs(s[0] - s[1]) < 1e-10
assert all(close(W.mH @ Y @ Wa, 0) for W in spaces if W.shape[1] == 1)   # different summands: zero
```

Die Kommutante bestimmt die Kosten: Eine Twirl-Stichprobe kostet $|G|$ Produkte von $n\times n$-Matrizen. Für Erzeuger
diagonalisiert die Bibliothek zunächst ein hermitesches Element $R=U\mathrm{diag}(\lambda)U^*$ von $\mathcal A$
(eine zufällige reelle Linearkombination der hermiteschen Teile der Erzeuger plus ein Vielfaches eines
symmetrisierten Produkts zweier solcher Kombinationen). Jedes $X\in\mathcal A'$ vertauscht mit $R$, also mit dessen
Spektralprojektionen (Satz 1.5(d)), sodass $U^*XU$ blockdiagonal bezüglich der Eigenräume von $R$ der Dimensionen
$d_\lambda$ ist; dies wird vorab erzwungen, und das Kernproblem hat $N_R=\sum_\lambda d_\lambda^2$ statt $n^2$ Unbekannte.
Das Zusammenfassen mehrerer Eigenwerte vergrößert nur die Blöcke, sodass das reduzierte Problem dieselbe
Lösungsmenge $\mathcal A'$ hat. In der Zerlegung ist $V^*RV=\bigoplus_cr_c\otimes1_{m_c}$ mit hermiteschen
$r_c\in M_{k_c}$; hat $\lambda$ in $r_c$ die Vielfachheit $\mu_c(\lambda)$, so ist $d_\lambda=\sum_c\mu_c(\lambda)m_c$, und da
$\mu^2\ge\mu$ für ganze Zahlen $\mu\ge0$ gilt,
$$
N_R=\sum_\lambda\Big(\sum_c\mu_c(\lambda)m_c\Big)^2\ge\sum_{\lambda,c}\mu_c(\lambda)^2m_c^2\ge\sum_{\lambda,c}\mu_c(\lambda)m_c^2=\sum_ck_cm_c^2,
$$
mit Gleichheit genau dann, wenn die Eigenwerte aller $r_c$ einfach und paarweise verschieden sind. Für den
Gesamtspin von $L=3,4,5$ Spins findet der folgende Code Gleichheit (numerische Beobachtung). Selbst das Minimum ist
groß: Für acht Spins $\frac12$ ($n=256$) ist $\sum_ck_cm_c^2=4900$, und das reduzierte System, eine
$(n^2+N_R)\times N_R$-Matrix, hat $3.5\cdot10^8$ komplexe Einträge. Für Spinketten verwendet Abschnitt 4.6 daher das
Modul `su2`, das die Antwort im Voraus kennt.

```python
for n in (3, 4, 5):
    Dn = dc.from_generators(dc.total_spin(n))
    _, info = dc.commutant_basis(dc.total_spin(n), return_info=True)
    N_min = sum(k * m * m for k, m in zip(Dn.k_factors, Dn.multiplicities))
    assert info['unknowns'] >= N_min                                  # N_R >= sum_c k_c m_c^2
    print(n, info['unknowns'], N_min)                                # observed: equality
```

<a id="s-4-4"></a>

## 4.4 $U(1)$ und Ladungssektoren

*Vorlesungsnotizen:* [Abschnitt 4.4](../../lecture_notes/de/04_symmetries.md#s-4-4); [Proposition 4.10](../../lecture_notes/de/04_symmetries.md#r-4-10).

Für $n$ Qubits ist die Ladung der Teilchenzahloperator $\hat N=\sum_jn_j$, $n_j=|1\rangle\langle1|$ auf Qubit $j$, und die
Operatoren, die sie erhalten, bilden $U(1)'\cong\bigoplus_{w=0}^nM_{\binom nw}(\mathbb C)$
([Vorlesungsnotizen, Abschnitt 4.4](../../lecture_notes/de/04_symmetries.md#s-4-4)); das ist `chains.SpinChain(n).algebra` (Bit $1$ = Spin aufwärts,
$\hat N=S^z+\frac n2$; in `decompose.total_spin` ist Spin aufwärts dagegen der *erste* Basisvektor von $\mathbb C^2$, was die
Basis umbenennt, nicht aber die Sektoren). Jeder XXZ-Hamiltonoperator liegt in $U(1)'$. Der Code findet die Sektoren
allein aus $\hat N$ und diagonalisiert einen XXZ-Ring aus vier Spins sektorweise.

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

<a id="s-4-5"></a>

## 4.5 $SU(2)$: Spin

*Vorlesungsnotizen:* [Abschnitt 4.5](../../lecture_notes/de/04_symmetries.md#s-4-5); [Satz 4.11](../../lecture_notes/de/04_symmetries.md#r-4-11), [Proposition 4.13](../../lecture_notes/de/04_symmetries.md#r-4-13), [Satz 4.14](../../lecture_notes/de/04_symmetries.md#r-4-14). Spin $s$ ist `decompose.spin_matrices(s)`, und der Gesamtspin von $L$ Spins $\frac12$ ist `decompose.total_spin(L)`.

Für $L=4$ liefert [Vorlesungsnotizen, Proposition 4.13](../../lecture_notes/de/04_symmetries.md#r-4-13) $\{S_x,S_y,S_z\}''\cong M_5\oplus M_3\oplus M_1$
mit den Vielfachheiten $1,3,2$, der Dimension $35$, und die Kommutante $\{S_x,S_y,S_z\}'\cong M_1\oplus M_3\oplus M_2$, der
Dimension $14$. Die Bibliothek findet dasselbe allein aus den drei Erzeugern (Sektoren sortiert nach $(k_c,m_c)$), und
$\mathbf S^2=S(S+1)$ auf Spin $S$ ([Satz 4.11](../../lecture_notes/de/04_symmetries.md#r-4-11)):

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

Nach [Vorlesungsnotizen, Satz 4.14](../../lecture_notes/de/04_symmetries.md#r-4-14) (Schur-Weyl-Dualität) spannen die Permutationen der Plätze die
Kommutante des Gesamtspins auf: $\mathrm{span}\,\pi(S_L)=\{S_x,S_y,S_z\}'$. Für $L\le5$ vergleicht der Code die beiden
Räume numerisch: Wegen Teil (a) des Satzes, $\mathrm{span}\,\pi(S_L)\subseteq\{S_x,S_y,S_z\}'$, genügt es, dass sie
dieselbe Dimension haben.

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

Die ausgegebenen Dimensionen sind $2,5,14,42$, die Catalan-Zahlen $\binom{2L}L/(L+1)$; Korollar B.10 beweist
$\dim\{S_a\}'=\sum_Sm_S^2=\binom{2L}L/(L+1)$ für jedes $L$. Die letzte Zeile zeigt die
Dualität in Aktion: $\pi(S_4)''$ hat Blöcke $k=(1,2,3)$ mit Vielfachheiten $(5,1,3)$ – die Paare $(k_c,m_c)$ des
Gesamtspins mit vertauschten Rollen. Allgemein (Korollar B.9) trägt der Vielfachheitsraum von Spin $S$ eine
irreduzible Darstellung von $S_L$ der Dimension $m_S$, und diese sind paarweise inäquivalent; hier sind die
Dimensionen $1,3,2$ für $S=2,1,0$.

<a id="s-4-6"></a>

## 4.6 Spinsysteme

*Vorlesungsnotizen:* [Abschnitt 4.6](../../lecture_notes/de/04_symmetries.md#s-4-6); [Proposition 4.15](../../lecture_notes/de/04_symmetries.md#r-4-15), [Satz 4.16](../../lecture_notes/de/04_symmetries.md#r-4-16), [Korollar 4.17](../../lecture_notes/de/04_symmetries.md#r-4-17).

Nach [Vorlesungsnotizen, Korollar 4.17](../../lecture_notes/de/04_symmetries.md#r-4-17) lassen sich Orthonormalbasen der Höchstgewichtsvektoren von
$L$ Spins $\frac12$ Platz für Platz mit den Clebsch-Gordan-Koeffizienten aus [Satz 4.16](../../lecture_notes/de/04_symmetries.md#r-4-16) aufbauen.
Bezeichnet man jeden Basisvektor durch seine Folge von Zwischenspins $S_1=\frac12,S_2,\dots,S_L=S$, so erhält man die
Yamanouchi-Basis, die `su2.SU2Chain` aufbaut (ihr Quellcode verwendet genau diese Koeffizienten; der folgende Code prüft
nur die Aussagen von Korollar 4.17), innerhalb des $U(1)$-Sektors $\hat N=\frac L2+S$ von `SpinChain` (dessen Anordnung
der Plätze eine Umbenennung der Tensorfaktoren ist, die nach Lemma B.1(b) mit dem Gesamtspin vertauscht). Der Code prüft
Satz 4.16 an dichten Matrizen und gegen `su2.clebsch_gordan` der Bibliothek, dann die Aussagen von Korollar 4.17 für die
Basen von acht Spins. Außerdem prüfen wir [Proposition 4.15(a)](../../lecture_notes/de/04_symmetries.md#r-4-15) an einem Ring aus vier Spins mit
`D4` und vergleichen die Blockgrößen für acht Spins.

```python
from torch_vn_algebra import su2

def ket(j, m):                                            # |j, m> of Theorem 4.11 (zero if |m| > j)
    v = torch.zeros(int(2 * j + 1), dtype=dt)
    if abs(m) <= j:
        v[int(j - m)] = 1
    return v

up, dn = ket(0.5, 0.5), ket(0.5, -0.5)
for j in (0.5, 1.0, 1.5, 2.0):                            # Theorem 4.16
    Sj = dc.spin_matrices(j)
    sp_j = [math.sqrt(j * (j + 1) - (j - k) * (j - k + 1)) for k in range(1, int(2 * j) + 1)]
    assert close(Sj[0] + 1j * Sj[1], torch.diag(torch.tensor(sp_j, dtype=dt), 1))     # S_+ of Theorem 4.11
    assert close(Sj[2], torch.diag(torch.tensor([j - k for k in range(int(2 * j) + 1)], dtype=dt)))
    Sc = torch.stack([torch.kron(a, torch.eye(2, dtype=dt)) + torch.kron(torch.eye(len(a), dtype=dt), b)
                      for a, b in zip(Sj, dc.spin_matrices(0.5))])
    Sp, Sm = Sc[0] + 1j * Sc[1], Sc[0] - 1j * Sc[1]
    cols = []
    for J in (j + 0.5, j - 0.5):
        for M in [J - k for k in range(int(2 * J + 1))]:
            al, be = math.sqrt((j + M + 0.5) / (2 * j + 1)), math.sqrt((j - M + 0.5) / (2 * j + 1))
            p, q = (al, be) if J > j else (-be, al)
            x = p * torch.kron(ket(j, M - 0.5), up) + q * torch.kron(ket(j, M + 0.5), dn)
            assert close(Sc[2] @ x, M * x)
            if M == J:
                assert close(Sp @ x, 0)                               # highest weight
            else:
                assert close(Sm @ cols[-1], math.sqrt((J + M + 1) * (J - M)) * x)   # lowering, (c)
            for s2, e in ((0.5, up), (-0.5, dn)):                     # the library's Racah formula agrees
                cg = su2.clebsch_gordan(j, M - s2, 0.5, s2, J, M) if abs(M - s2) <= j else 0.0
                assert abs(cg - torch.vdot(torch.kron(ket(j, M - s2), e), x).real.item()) < 1e-12
            cols.append(x)
    U = torch.stack(cols, 1)
    assert U.shape == (2 * len(Sj[0]), 2 * len(Sj[0])) and close(U.mH @ U, torch.eye(U.shape[1]))   # (a)

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

u1, sym = SpinChain(8, boundary='periodic'), su2.SU2Chain(8, boundary='periodic')
for S, m in zip(sym.spins, sym.multiplicities):           # Corollary 4.17: orthonormal highest-weight bases
    B, N = sym.hw_basis(S), int(4 + S)                     # columns in the U(1) sector N = L/2 + S
    assert B.shape == (math.comb(8, N), m) and close(B.T @ B, torch.eye(m, dtype=B.dtype))
    assert N == 8 or close(sym.raise_(B, N), 0)            # S_+ B = 0 (S_z = S holds in sector N)
H_u1, H_su2 = u1.xxz(J=1.0, Delta=1.0), sym.heisenberg(J=1.0)
E_u1 = torch.sort(torch.cat([torch.linalg.eigvalsh(H_u1.matrix[0, N, :k, :k])
                             for N, k in enumerate(u1.algebra.k_factors)])).values
E_su2, S_lev, _ = sym.full_spectrum(H_su2, return_quantum_numbers=True)   # level of block S, 2S+1 times
E_su2, S_lev = E_su2[0], S_lev[0]
assert close(E_u1, E_su2)
assert S_lev[0] == 0 and (S_lev[1:4] == 1).all() and E_su2[4] - E_su2[3] > 1e-3     # singlet, then triplet
assert E_su2[1] - E_su2[0] > 1e-3 and close(E_su2[1:4], E_su2[1].expand(3))       # simple singlet (Theorem 9.9)
assert close(sym.restrict(H_u1).matrix, H_su2.matrix)      # Proposition 4.15(b)
flops = lambda ks: sum(k ** 3 for k in ks)
print(u1.algebra.k_factors, sym.multiplicities, flops(u1.algebra.k_factors), flops(sym.multiplicities))
print(f"ground state energy {E_su2[0].item():.6f}, first excited {E_su2[1].item():.6f} (x3)")
```

Für den Ring aus acht Spins haben die $U(1)$-Blöcke die Größen $\binom8w=1,8,28,56,70,56,28,8,1$ (Abschnitt 4.4), die
$SU(2)$-Blöcke die Größen $m_S=14,28,20,7,1$ für $S=0,\dots,4$ (Proposition 4.13). Der Aufwand der Diagonalisierung,
$\sum k^3$, sinkt von $739162$ auf $33040$, um einen Faktor $22$; die Dimension der Algebra, $\sum k^2$, von
$\sum_w\binom8w\binom8{8-w}=\binom{16}8=12870$ (man wähle $8$ von $16$ Elementen, $w$ davon unter den ersten $8$) auf
$\sum_Sm_S^2=1430=\binom{16}8/9$ (Korollar B.10). Nach Satz 9.9 (der Ring hat gerade Länge und $J>0$) ist der
Grundzustand ein einzelnes Singulett, und $E_0(0)<E_0(1)$. *Numerische Beobachtung* (berechnet vom obigen Code): Die
Grundzustandsenergie ist $-3.651093$, und das erste angeregte Niveau ist ein Triplett bei $-3.128419$.

<a id="s-4-7"></a>

## 4.7 Reelle Algebren: reeller, komplexer und quaternionischer Typ

*Vorlesungsnotizen:* [Abschnitt 4.7](../../lecture_notes/de/04_symmetries.md#s-4-7). Für eine unter komplexer Konjugation abgeschlossene
$*$-Algebra $\mathcal A$ mit einer Zerlegung $V$, deren Spaltenblöcke $V_c$ zu den Sektoren gehören, liefert
[Satz 4.20](../../lecture_notes/de/04_symmetries.md#r-4-20) die Permutation $\sigma$ der Sektoren mit $\bar P_c=P_{\sigma(c)}$, unitäre $J_c$, $M_c$ mit
$\overline{V_c}=V_{\sigma(c)}(J_c\otimes M_c)$ und, für $\sigma(c)=c$, das Vorzeichen $\varepsilon_c$ mit $J_c\bar J_c=\varepsilon_c1$;
die Typen und die reelle Algebra $\mathcal A_{\mathbb R}$ stehen in [Satz 4.21](../../lecture_notes/de/04_symmetries.md#r-4-21), das reelle Schursche
Lemma ist [Satz 4.22](../../lecture_notes/de/04_symmetries.md#r-4-22), und die reellen Elemente der Kommutante beschreibt
[Proposition 4.23](../../lecture_notes/de/04_symmetries.md#r-4-23).

`Decomposition.real_structure()` berechnet $\sigma$, $J_c$, $M_c$ und die Typen für eine von `wedderburn`
(`from_generators`, `from_group`) gelieferte Zerlegung. Sie prüft, dass die Konjugierten der gespeicherten Erzeuger in
$V(\bigoplus_cM_{k_c}\otimes1)V^*$ liegen (relatives Residuum `'conjugation'`), und bildet die unitäre Matrix
$P=V^*\bar V$, deren Spaltenblock $c$ gleich $V^*\overline{V_c}$ ist: Nach Satz 4.20(a) ist er auf den Zeilen des
Sektors $\sigma(c)$ getragen und dort gleich $J_c\otimes M_c$. Die Bibliothek nimmt als $\sigma(c)$ den Sektor, der
den größten Teil des Blocks trägt, verlangt $(k_{\sigma(c)},m_{\sigma(c)})=(k_c,m_c)$ und schätzt den Rest ab
(`'leak'`); sie ordnet den Block zur $k_c^2\times m_c^2$-Matrix mit den Einträgen $T_{(a,b),(j,l)}=T_{(a,j),(b,l)}$
um, die für $T=J\otimes M$ die Rang-eins-Matrix $\mathrm{vec}\,J\,(\mathrm{vec}\,M)^T$ ist, liest $J_c$ und
$M_c$ an ihren führenden Singulärvektoren ab (mit $\|J_c\|_F^2=k_c$) und speichert den Abstand des Blocks von
$J_c\otimes M_c$ (`'kron'`). Überschreitet eines dieser drei Residuen die Toleranz, passen die Größen nicht zusammen
oder bilden die Partner keine Involution, so löst sie `ValueError` aus (die Algebra ist nicht unter Konjugation
abgeschlossen). Für $\sigma(c)=c$ setzt sie $\varepsilon_c$ gleich dem Vorzeichen von
$\mathrm{Re}\,\mathrm{Tr}(J_c\bar J_c)/k_c$ – in exakter Arithmetik genau $\pm1$, sodass die beiden
Kandidaten den Abstand $2$ haben – und speichert als `'sign'` das Maximum über diese Sektoren von
$\|J_c\bar J_c-\varepsilon_c\|_F/\sqrt{k_c}$ und $\|M_c\bar M_c-\varepsilon_c\|_F/\sqrt{m_c}$, auf mindestens $1$
angehoben, wenn ein $\varepsilon_c=-1$ mit ungeradem $k_c$ oder $m_c$ auftritt (nach Satz 4.20(b) unmöglich). Das
Ergebnis hat `types`, `partner` ($\sigma$), `J`, `M`, `residuals` und `ok`; letzteres besagt, dass alle Residuen
höchstens gleich der Toleranz sind (standardmäßig die Annahmeschwelle der Zerlegung); dann ist auch
$|\mathrm{Tr}(J_c\bar J_c)/k_c-\varepsilon_c|\le\|J_c\bar J_c-\varepsilon_c\|_F/\sqrt{k_c}$ (Cauchy–Schwarz)
höchstens gleich der Toleranz.

`Decomposition.real_basis()` liefert eine reelle orthogonale $n\times n$-Matrix $O$ und eine Liste, die ihre
Spaltenblöcke beschreibt (Typ, Sektoren, Spaltenbereich): für einen Sektor vom reellen Typ die Spalten
$O_c=V_c(W_J\otimes W_M)$ aus Satz 4.21(a), sodass $O_c^TaO_c=\alpha_c(a)\otimes1_{m_c}$ mit reellem $\alpha_c(a)$;
für einen Sektor vom quaternionischen Typ die Spalten $V_cW$ mit $J_c\otimes M_c=WW^T$ und $(J_c\otimes M_c)\bar W=W$
(Lemma 4.19(c) mit $\varepsilon=1$, da $(J_c\otimes M_c)\overline{(J_c\otimes M_c)}=\varepsilon_c^2=1$), die nach der
Rechnung aus Satz 4.21(a) reell sind; und für ein Paar $c,\sigma(c)$ (mit $c<\sigma(c)$) die Spalten
$\sqrt2\mathrm{Re}\,V_c$, gefolgt von $\sqrt2\mathrm{Im}\,V_c$. Diese sind orthonormal:
$V_c^*\overline{V_c}=V_c^*V_{\sigma(c)}(J_c\otimes M_c)=0$ und $V_c^T\overline{V_c}=\overline{V_c^*V_c}=1$ ergeben
zusammen mit ihren Konjugierten
$(\mathrm{Re}\,V_c)^T\mathrm{Re}\,V_c=(\mathrm{Im}\,V_c)^T\mathrm{Im}\,V_c=\frac12$ und
$(\mathrm{Re}\,V_c)^T\mathrm{Im}\,V_c=0$; Spalten verschiedener Blöcke liegen in den Bildern der
Orthogonalprojektionen $P_{\mathcal O}=\sum_{e\in\mathcal O}P_e$ verschiedener Bahnen $\mathcal O$ von $\sigma$ (da
$\overline{V_c}=V_{\sigma(c)}(J_c\otimes M_c)$), also ist $O$ orthogonal. Real- und Imaginärteil von $aV_c=V_cx$,
$x=\pi_c(a)\otimes1_{m_c}$, ergeben für reelles $a$
$$
O^TaO=\begin{pmatrix}\mathrm{Re}\,x&\mathrm{Im}\,x\\-\mathrm{Im}\,x&\mathrm{Re}\,x\end{pmatrix}\quad\text{on these columns.}
$$
Die Takagi-Faktoren $W$ mit $WW^T=J$ einer symmetrischen unitären Matrix $J=X+iY$ werden aus einer orthogonalen
Diagonalisierung von $X+tY$ für ein zufälliges $t$ berechnet ($X$ und $Y$ sind vertauschende reelle symmetrische
Matrizen, da der Imaginärteil von $J\bar J=JJ^*=1$ gleich $YX-XY$ ist); die Bibliothek akzeptiert das Ergebnis anhand
seiner Residuen $\|WW^T-J\|_F/\|J\|_F$ und $\|\mathrm{Im}(V_cW)\|_F/\sqrt{k_cm_c}$ und warnt
(`DecompositionWarning`), wenn eines davon $\max(\text{Toleranz},10^3\,\epsilon_{\rm mach}\,n)$ überschreitet.

[Abschnitt 4.7](../../lecture_notes/de/04_symmetries.md#s-4-7) der Vorlesungsnotizen behandelt drei von reellen Matrizen erzeugte Algebren: (i)
$L_a\otimes1_2$ auf $\mathbb R^3\otimes\mathbb R^2$, mit den Erzeugern $(L_a)_{bc}=-\epsilon_{abc}$ der Drehungen des
$\mathbb R^3$: ein Sektor vom reellen Typ, $k=3$, $m=2$; (ii) $\Omega_2\otimes1_2$ auf $\mathbb R^2\otimes\mathbb R^2$: ein
Paar vom komplexen Typ, $k=1$, $m=2$; (iii) die Linksmultiplikationen $1_2\otimes L_q$ mit den Quaternionen auf
$\mathbb H^2=\mathbb R^8$: ein Sektor vom quaternionischen Typ, $k=2$, $m=4$. Der Code berechnet für jede die Zerlegung
und die reelle Struktur und prüft Satz 4.20(a),(b) sowie die Blockform von $O^TaO$ für den reellen Typ (Satz 4.21(a))
und für das Paar.

```python
def qmul(p, q):                                           # H = CD(C), Section E.6: (a, b)(c, d)
    (a1, b1), (a2, b2) = p, q
    return (a1 * a2 - b2.conjugate() * b1, b2 * a1 + b1 * a2.conjugate())

def q_left(q):                                            # left multiplication by q on H = R^4
    cols = []                                             # coordinates (Re a, Im a, Re b, Im b)
    for e in torch.eye(4, dtype=torch.float64).tolist():
        u, v = qmul(q, (complex(e[0], e[1]), complex(e[2], e[3])))
        cols.append([u.real, u.imag, v.real, v.imag])
    return torch.tensor(cols, dtype=torch.float64).T

L_so3 = torch.zeros(3, 3, 3, dtype=torch.float64)         # (L_a)_bc = -epsilon_abc
for i, (b1, b2) in enumerate([(1, 2), (2, 0), (0, 1)]):
    L_so3[i, b1, b2], L_so3[i, b2, b1] = -1.0, 1.0
R_rot = torch.tensor([[[0.0, -1.0], [1.0, 0.0]]])         # Omega_2
L_H = torch.stack([q_left(u) for u in ((1j, 0), (0, 1), (0, 1j))])            # L_I, L_J, L_K
assert close(L_H[0] @ L_H[1], L_H[2]) and close(q_left((2 - 1j, 3j)).T, q_left((2 + 1j, -3j)))   # IJ = K, L_qbar = L_q^T
eye2 = torch.eye(2, dtype=torch.float64)
real_examples = {'R': (torch.kron(L_so3, eye2[None]), (3,), (2,)),            # L_a (x) 1_2 on R^3 (x) R^2
                 'C': (torch.kron(R_rot, eye2[None]), (1, 1), (2, 2)),        # Omega_2 (x) 1_2
                 'H': (torch.kron(eye2[None], L_H), (2,), (4,))}              # L_q on H^2 = R^8
dec_real = {}
for t, (gens, ks, ms) in real_examples.items():
    D_t = dc.from_generators(gens)
    rs_t, (O_t, blocks_t) = D_t.real_structure(), D_t.real_basis()
    print(t, D_t, rs_t)
    assert D_t.verify()['ok'] and rs_t.ok and (D_t.k_factors, D_t.multiplicities) == (ks, ms)
    assert set(rs_t.types) == {t} and close(O_t.T @ O_t, torch.eye(D_t.n))     # real orthogonal basis
    for c, d in enumerate(rs_t.partner):                                      # Theorem 4.20
        o_c, o_d, w = D_t.offsets[c], D_t.offsets[d], ks[c] * ms[c]
        assert close(D_t.V[:, o_c:o_c + w].conj(), D_t.V[:, o_d:o_d + w] @ torch.kron(rs_t.J[c], rs_t.M[c]))
        if d == c:                                                            # J conj(J) = M conj(M) = +-1
            eps = 1.0 if t == 'R' else -1.0
            assert close(rs_t.J[c] @ rs_t.J[c].conj(), eps * torch.eye(ks[c]))
            assert close(rs_t.M[c] @ rs_t.M[c].conj(), eps * torch.eye(ms[c]))
    dec_real[t] = D_t, rs_t, O_t

D_t, rs_t, O_t = dec_real['R']
for g in real_examples['R'][0]:                           # Theorem 4.21(a): O^T a O = alpha (x) 1_2, alpha real
    G_t = O_t.T @ g @ O_t
    assert close(G_t, torch.kron(G_t[0::2, 0::2], eye2))
D_t, rs_t, O_t = dec_real['C']                            # the pair: [[Re x, Im x], [-Im x, Re x]]
x_t = torch.kron(D_t.to_algebra(real_examples['C'][0]).matrix[0, 0, :1, :1], torch.eye(2, dtype=dt))
assert close(O_t.T @ real_examples['C'][0][0] @ O_t,
             torch.cat([torch.cat([x_t.real, x_t.imag], 1), torch.cat([-x_t.imag, x_t.real], 1)]))
```

Der Code prüft Satz 4.22 an den irreduziblen reellen Moduln $\mathbb R^3$, $\mathbb R^2$, $\mathbb R^4$ der drei
Beispiele (die von der Bibliothek gefundenen Typen und Vielfachheiten – reell mit $m=1$, ein komplexes Paar mit $m=1$,
quaternionisch mit $m=2$ – machen sie nach Satz 4.22(c) irreduzibel), wobei er die reelle Kommutante als Nullraum der
reellen linearen Gleichungen $Xs=sX$, $Xs^T=s^TX$ berechnet, und danach Proposition 4.23 für ein zufälliges reelles
symmetrisches Element der Kommutante in jedem Beispiel.

```python
def real_commutant_dim(gens):                             # dim_R {X real : X s = s X, X s^T = s^T X}
    gens = torch.cat([gens, gens.mT])
    eye_n = torch.eye(gens.shape[-1], dtype=torch.float64)
    eqs = torch.cat([torch.kron(eye_n, s) - torch.kron(s.T.contiguous(), eye_n) for s in gens])
    return gens.shape[-1] ** 2 - torch.linalg.matrix_rank(eqs).item()

for gens, types, ms, dim_D in ((L_so3, ('R',), (1,), 1), (R_rot, ('C', 'C'), (1, 1), 2), (L_H, ('H',), (2,), 4)):
    D_irr = dc.from_generators(gens)                      # Theorem 4.22(c): R^3, R^2, R^4 are irreducible
    assert D_irr.real_structure().types == types and D_irr.multiplicities == ms
    assert real_commutant_dim(gens) == dim_D              # Theorem 4.22(b): D_U = R, C, H

gen_X = torch.Generator().manual_seed(1)
for t, (gens, ks, ms) in real_examples.items():           # Proposition 4.23
    D_t, rs_t, O_t = dec_real[t]
    B_t = dc.commutant_basis(gens)
    X_t = torch.einsum('r,rij->ij', torch.randn(len(B_t), dtype=dt, generator=gen_X), B_t).real
    X_t = X_t + X_t.T                                     # real symmetric, in A' (Lemma 4.18(b))
    assert close(X_t @ gens, gens @ X_t)
    x_t = D_t.to_commutant(X_t).matrix[0]
    ev = [torch.linalg.eigvalsh(x_t[c, :m, :m]) for c, m in enumerate(ms)]
    if t == 'R':                                          # (a): O^T X O = 1_3 (x) y, y real symmetric
        Y_t = O_t.T @ X_t @ O_t
        assert close(Y_t, torch.kron(torch.eye(3, dtype=torch.float64), Y_t[:2, :2]))
    if t == 'C':                                          # (c): the two sectors have the same spectrum
        assert close(ev[0], ev[1])
    if t == 'H':                                          # (b): Kramers pairs, each level 2 k_c = 4 times in X
        assert close(ev[0][0::2], ev[0][1::2])
        assert close(torch.linalg.eigvalsh(X_t), ev[0].repeat_interleave(2))
        print('quaternionic block x_c:', [round(v, 4) for v in ev[0].tolist()])
```

<a id="s-4-8"></a>

## 4.8 Womit man vergleichen sollte

Für endliche Gruppen berechnen GAP und SageMath Charaktertafeln exakt, und Korollar 4.5 liefert dann die
Vielfachheiten; RepLAB (MATLAB; D. Rosset, F. Montealegre-Mora, J.-D. Bancal) zerlegt Darstellungen endlicher und
kompakter Gruppen numerisch, ebenfalls mit zufälligen Elementen der Kommutante, und SymbolicWedderburn.jl (Julia)
berechnet Wedderburn-Zerlegungen zur Symmetriereduktion semidefiniter Programme. In der Physik bringen QuSpin und
andere Codes für exakte Diagonalisierung $U(1)$- und Gittersymmetrien in Blockdiagonalform, ITensor und TeNPy
führen abelsche Ladungen in Tensornetzwerken mit, und TensorKit.jl behandelt $SU(2)$ mit Fusionsbäumen. Das Modul
`decompose` der Bibliothek ist ein dichtes Verfahren für kleine Matrizen ($n$ bis zu einigen Hundert); es akzeptiert
jede Liste von Matrizen (Gruppenelemente oder Erzeuger), zertifiziert sein Ergebnis (Proposition 4.9) und gibt die
Algebra und ihre Kommutante als `TypeIAlgebra`-Objekte mit `to_algebra`, `from_algebra`, `to_commutant`,
`from_commutant` zurück, sodass die Kapitel 3 und 5–8 unmittelbar anwendbar sind. Für Spin-$\frac12$-Ketten liefert
`su2.SU2Chain` die $SU(2)$-Blöcke direkt (Korollar 4.17); es akzeptiert bis zu $L=18$ Spins, eine Grenze, die durch
den Speicherbedarf der Basen gesetzt ist.

<a id="x-4"></a>

## Aufgaben

Die [Aufgaben](../../lecture_notes/de/04_symmetries.md#x-4) der Vorlesungsnotizen, numerisch geprüft (die Lösungen unten folgen dieser Reihenfolge):
1. $\dim\pi(S_L)'=\binom{d^2+L-1}L$ für $d=2$ und $L=3,4$ (die Werte $20$ und $35$);
2. $1\otimes\frac12\cong\frac12\oplus\frac32$ und $1\otimes1\cong0\oplus1\oplus2$, mit `from_generators`;
3. die Kommutante von $e^{i\theta\hat N}$ auf vier Qubits für $\theta=1$ (die Ladungen $0,\dots,4$ ergeben verschiedene
   $e^{iq}$, da $0<|q-q'|\le4<2\pi$) und für $\theta=\pi$;
4. die Zeeman-Aufspaltung: das Spektrum von $H+hS_z$ für den Heisenberg-Ring aus acht Spins, mit `su2.SU2Chain` und in
   den $U(1)$-Sektoren;
5. die Sektoren der regulären Darstellung von $Q_8$, mit `real_structure`.

<details>
<summary>Lösungen</summary>

```python
for n, value in ((3, 20), (4, 35)):
    assert dc.commutant_dimension(dc.permutation_representation(n, 2)) == math.comb(n + 3, n) == value
```

```python
def coupled(a, b):                                        # S^(a) (x) 1 + 1 (x) S^(b)
    Sa, Sb = dc.spin_matrices(a), dc.spin_matrices(b)
    ea, eb = torch.eye(Sa.shape[-1], dtype=dt), torch.eye(Sb.shape[-1], dtype=dt)
    return torch.stack([torch.kron(x, eb) + torch.kron(ea, y) for x, y in zip(Sa, Sb)])
for a, b, ks in ((1.0, 0.5, (2, 4)), (1.0, 1.0, (1, 3, 5))):  # 1 (x) 1/2 = 1/2 (+) 3/2, 1 (x) 1 = 0 (+) 1 (+) 2
    D_cg = dc.from_generators(coupled(a, b))
    assert D_cg.k_factors == ks and set(D_cg.multiplicities) == {1}
```

```python
for theta, m in ((1.0, (1, 1, 4, 4, 6)), (math.pi, (8, 8))):   # theta = 1: e^{iq}, q = 0..4, distinct as 4 < 2 pi
    D_theta = dc.from_generators(torch.linalg.matrix_exp(1j * theta * N_op))
    assert D_theta.multiplicities == m and set(D_theta.k_factors) == {1}
```

```python
assert not close(H4 @ Sx - Sx @ H4 + 0.3 * (Sz @ Sx - Sx @ Sz), 0)          # [H4 + 0.3 S_z, S_x] != 0
E_field = sym.full_spectrum(sym.heisenberg(J=1.0, h=0.3))[0]                # E + h M on every multiplet
H_field = u1.xxz(J=1.0, Delta=1.0, h=0.3)
E_u1_field = torch.cat([torch.linalg.eigvalsh(H_field.matrix[0, N, :k, :k])
                        for N, k in enumerate(u1.algebra.k_factors)])
assert close(torch.sort(E_u1_field).values, E_field)
```

```python
Q8 = [(s * a, s * b) for a, b in ((1, 0), (1j, 0), (0, 1), (0, 1j)) for s in (1, -1)]
key = lambda q: tuple(round(v) for z in q for v in (complex(z).real, complex(z).imag))
pos = {key(q): i for i, q in enumerate(Q8)}
reg = torch.zeros(8, 8, 8, dtype=torch.float64)           # pi(g) e_h = e_{gh}
for i, g in enumerate(Q8):
    for j, h in enumerate(Q8):
        reg[i, pos[key(qmul(g, h))], j] = 1.0
D_Q8 = dc.from_group(reg)
rs_Q8 = D_Q8.real_structure()
assert rs_Q8.ok and sorted(zip(rs_Q8.types, D_Q8.k_factors, D_Q8.multiplicities)) == [('H', 2, 2)] + [('R', 1, 1)] * 4
assert D_Q8.dims[0] == 8                                  # dim_R A_R = dim_C A = |Q_8| (Theorem 4.21)
```

</details>

**Weiter:** [Tutorial 5](05_random_matrices.md): Zufallsmatrizen in einer Algebra – Haar-verteilte unitäre Matrizen, der Weingarten-Kalkül und Niveaustatistik in der Bibliothek.
