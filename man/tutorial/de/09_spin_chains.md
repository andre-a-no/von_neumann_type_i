[English](../09_spin_chains.md) · [Español](../es/09_spin_chains.md) · [Français](../fr/09_spin_chains.md) · **Deutsch** · [中文](../zh/09_spin_chains.md) · [日本語](../ja/09_spin_chains.md) · [Русский](../ru/09_spin_chains.md)

# Tutorial 9. Spinketten

Dieses Tutorial begleitet [Kapitel 9](../../lecture_notes/de/09_spin_chains.md) der Vorlesungsnotizen, das die Sektoren der XXZ-Kette
aufbaut – Magnetisierung ([Proposition 9.1](../../lecture_notes/de/09_spin_chains.md#r-9-1)), Impuls ([Proposition 9.5](../../lecture_notes/de/09_spin_chains.md#r-9-5)) und Gesamtspin – und beweist, was
sie über Grundzustände aussagen ([Satz 9.7](../../lecture_notes/de/09_spin_chains.md#r-9-7), [Satz 9.9](../../lecture_notes/de/09_spin_chains.md#r-9-9)). Wir bauen die Ketten mit `SpinChain` und
`SU2Chain`, prüfen die Bibliothek gegen freie Fermionen und den Bethe-Ansatz, verwenden die Lanczos- und
Krylov-Verfahren aus `krylov` für Sektoren, die für dichte Matrizen zu groß sind, und berechnen Niveaustatistik und
Verschränkung.

Im ganzen Kapitel sind die Plätze $`0,\dots,L-1`$, die Basis von $`(\mathbb C^2)^{\otimes L}`$ besteht aus den
Konfigurationen $`|b\rangle=|b_0\dots b_{L-1}\rangle`$ mit $`b_i=1`$ für Spin nach oben (Konvention der Bibliothek: Bit
$`i`$ eines Labels ist Platz $`i`$), $`s^a_i`$ ($`a=x,y,z`$) sind die Spinmatrizen $`\frac12\sigma_a`$ auf Platz $`i`$,
$`s^\pm_i=s^x_i\pm is^y_i`$, $`n_i=s^z_i+\frac12`$, $`S^a=\sum_is^a_i`$, und $`N=S^z+\frac L2`$ ist die Anzahl der Spins nach
oben. Eine Kette ist **offen** (Bindungen $`(i,i+1)`$, $`i<L-1`$) oder **periodisch** (zusätzlich die Bindung
$`(L-1,0)`$; dann ist $`L\ge3`$, und Indizes werden modulo $`L`$ genommen). $`J>0`$ (antiferromagnetisch), $`\hbar=1`$, und
$`\log`$ ist der natürliche Logarithmus.

<a id="s-9-1"></a>

## 9.1 Die XXZ-Kette und ihre $`U(1)`$-Sektoren

*Vorlesungsnotizen:* [Abschnitt 9.1](../../lecture_notes/de/09_spin_chains.md#s-9-1); [Proposition 9.1](../../lecture_notes/de/09_spin_chains.md#r-9-1) ($`H`$ erhält $`S^z`$), [Proposition 9.2](../../lecture_notes/de/09_spin_chains.md#r-9-2), [Lemma 9.3](../../lecture_notes/de/09_spin_chains.md#r-9-3) und [Korollar 9.4](../../lecture_notes/de/09_spin_chains.md#r-9-4) (der Turm der Ketten und sein Index).

`SpinChain(L).xxz` baut die Blöcke direkt aus den Konfigurationen jedes Sektors (`chain.labels[N]`). Wir
vergleichen sie mit der vollen $`2^L\times2^L`$-Matrix, die aus Kronecker-Produkten gebildet wird:

```python
import functools, itertools, math, time
import numpy as np
import torch
from scipy.optimize import brentq
from torch_vn_algebra import SpinChain, SU2Chain, DensityMatrix, inclusions, krylov

torch.manual_seed(0)
f64, c128 = torch.float64, torch.complex128
close = lambda a, b, tol=1e-10: bool((torch.as_tensor(a) - torch.as_tensor(b)).abs().max() < tol)
def op(L, ops):                             # product of one-site operators {site: 2 x 2}; site 0 = lowest bit
    return functools.reduce(torch.kron, [ops.get(j, torch.eye(2, dtype=f64)) for j in reversed(range(L))])

sp = torch.tensor([[0., 0.], [1., 0.]], dtype=f64)                    # s^+ : |0> = down -> |1> = up
sm, sz = sp.T.contiguous(), torch.diag(torch.tensor([-0.5, 0.5], dtype=f64))
L = 6
chain = SpinChain(L, boundary='periodic')
h = 0.3 * torch.randn(L, dtype=f64)                                  # random fields
H_full = sum(0.5 * (op(L, {i: sp, j: sm}) + op(L, {i: sm, j: sp})) + 0.7 * op(L, {i: sz, j: sz})
             for i, j in chain.bonds) + sum(h[i] * op(L, {i: sz}) for i in range(L))   # J = 1, Delta = 0.7
Sz = sum(op(L, {i: sz}) for i in range(L))
assert close(H_full @ Sz, Sz @ H_full)                               # Proposition 9.1
H = chain.xxz(J=1.0, Delta=0.7, h=h)                                 # an element of (+)_N M_binom(6, N)
for N, k in enumerate(chain.algebra.k_factors):
    lab = chain.labels[N]                                            # configurations of sector N
    assert close(H.matrix[0, N, :k, :k], H_full[lab][:, lab])
```

Der Turm $`\mathcal A_0\subset\mathcal A_1\subset\cdots`$ der Ketten: `inclusions.chain_tower` baut die Inklusionen,
`inclusions.bratteli_diagram` ihr Bratteli-Diagramm und `Inclusion.index` den Index $`\|\Lambda\|^2`$. Der Code prüft
die Pascalsche Regel, [Lemma 9.3](../../lecture_notes/de/09_spin_chains.md#r-9-3) und $`[\mathcal A_{L+1}:\mathcal A_L]=4\cos^2\frac\pi{2L+4}`$ ([Korollar 9.4](../../lecture_notes/de/09_spin_chains.md#r-9-4)).

```python
tower = inclusions.chain_tower(8)                                    # A_0 c A_1 c ... c A_8
print(inclusions.bratteli_diagram(tower).levels[:6])                 # Pascal's triangle
for L_, inc in enumerate(tower):
    sizes = [math.comb(L_ + 1, N) for N in range(L_ + 2)]
    assert inc.M.k_factors == sizes and inc.is_connected
    assert list(inc.matrix @ np.array(inc.N.k_factors)) == sizes      # unitality = Pascal's rule
    n = 2 * L_ + 3
    A_path = np.eye(n, k=1) + np.eye(n, k=-1)
    assert abs(np.linalg.eigvalsh(A_path).max() - 2 * math.cos(math.pi / (n + 1))) < 1e-12   # Lemma 9.3
    assert abs(inc.index() - 4 * math.cos(math.pi / (2 * L_ + 4)) ** 2) < 1e-12           # Corollary 9.4
```

<a id="s-9-2"></a>

## 9.2 Translationen, Impuls und die Marshallsche Vorzeichenregel

*Vorlesungsnotizen:* [Abschnitt 9.2](../../lecture_notes/de/09_spin_chains.md#s-9-2); [Proposition 9.5](../../lecture_notes/de/09_spin_chains.md#r-9-5) (Impulssektoren), [Satz 9.7](../../lecture_notes/de/09_spin_chains.md#r-9-7) (Marshallsche Vorzeichenregel), [Korollar 9.8](../../lecture_notes/de/09_spin_chains.md#r-9-8) (Impuls des Grundzustands).

`xxz_momentum` baut die Blöcke in der normierten Basis $`|a,k\rangle`$, mit einem Repräsentanten $`a`$ pro Bahn. Wir
zählen die Bahnen unabhängig ab und vergleichen die Spektren mit dem dichten Sektor:

```python
L = 12
ring = SpinChain(L, boundary='periodic')
rot = lambda x: ((x << 1) | (x >> (L - 1))) & ((1 << L) - 1)        # T on labels: site i -> i + 1
for N in range(L + 1):
    orbits = {frozenset(itertools.accumulate(range(L - 1), lambda x, _: rot(x), initial=a))
              for a in ring.labels[N].tolist()}                       # {a, Ta, ..., T^{L-1} a}
    d = [sum((m * len(o)) % L == 0 for o in orbits) for m in range(L)]
    assert d == ring.momentum_dims(N)                                   # Proposition 9.5(b)
print(ring.momentum_dims(6))
Hk = ring.xxz_momentum(J=1.0, Delta=1.0, h=0.2, sector=6)               # 12 blocks instead of one of size 924
E_k = torch.sort(torch.cat([e[0] for e in Hk.eigenvalues()])).values
assert close(E_k, torch.linalg.eigvalsh(ring.xxz(1.0, 1.0, h=0.2, sector=6).matrix[0, 0]))
```

Der Sektor $`N=6`$ mit $`924`$ Zuständen zerfällt in Blöcke der Größen $`80,75,78,76,78,75,80,75,78,76,78,75`$.

Nach der Marshallschen Vorzeichenregel ([Satz 9.7](../../lecture_notes/de/09_spin_chains.md#r-9-7)) hat $`(-1)^{N_B(b)}\psi(b)`$ für offene und periodische Ketten mit
Zufallsfeldern in jedem Sektor ein einziges Vorzeichen, und nach [Korollar 9.8](../../lecture_notes/de/09_spin_chains.md#r-9-8) liegt das tiefste Niveau des Sektors $`N`$
des Rings im Impulsblock $`k=0`$ für gerades $`N`$ und $`k=\pi`$ für ungerades $`N`$:

```python
for bc in ('open', 'periodic'):
    c = SpinChain(10, bc, complex_valued=False)
    for N in range(11):
        w, V = torch.linalg.eigh(c.xxz(1.0, 0.7, h=0.3 * torch.randn(10, dtype=f64), sector=N).matrix[0, 0])
        x = (-1.0) ** c.bits(N)[:, 1::2].sum(-1) * V[:, 0]                  # (-1)^{N_B(b)} psi(b)
        assert (x * x.sum().sign()).min() > 0 and (len(w) == 1 or w[1] - w[0] > 1e-8)   # Theorem 9.7
for N in range(1, L):
    lowest = [e[0, 0].item() for e in ring.xxz_momentum(1.0, 0.5, sector=N).eigenvalues()]
    m0 = int(np.argmin(lowest))
    assert m0 == (0 if N % 2 == 0 else L // 2)                              # Corollary 9.8: k = 0 or pi
    assert sorted(lowest)[1] - lowest[m0] > 1e-8
```

<a id="s-9-3"></a>

## 9.3 $`SU(2)`$ bei $`\Delta=1`$

*Vorlesungsnotizen:* [Abschnitt 9.3](../../lecture_notes/de/09_spin_chains.md#s-9-3); [Satz 9.9](../../lecture_notes/de/09_spin_chains.md#r-9-9) (Lieb-Mattis). Für $`\Delta=1`$ und $`h=0`$ ist $`H\cong\bigoplus_S1_{2S+1}\otimes H_S`$, und `SU2Chain` baut die Blöcke $`H_S`$ auf den Vektoren höchsten Gewichts.

Wir vergleichen die Blockgrößen für zwölf Spins, prüfen Satz 9.9 an den $`U(1)`$-Sektoren, berechnen die Lücke und
prüfen die Sätze 9.7 und 9.9 auf einem bipartiten Graphen, der keine Kette ist ($`n_A=5`$, $`n_B=3`$, $`S_*=1`$, zufällige
Kopplungen):

```python
sym = SU2Chain(12, boundary='periodic', complex_valued=False)
H_su2 = sym.heisenberg(J=1.0)                                       # (+)_S M_{m_S}, S = 0, ..., 6
u1_sizes = [math.comb(12, N) for N in range(13)]
print(sym.multiplicities, max(u1_sizes), sum(k ** 3 for k in u1_sizes), sum(m ** 3 for m in sym.multiplicities))
E0_S = sym.lowest_energies(H_su2)[0]                                # E_0(S), S = 0, ..., 6
E0_N = [torch.linalg.eigvalsh(ring.xxz(1.0, 1.0, sector=6 + S).matrix[0, 0])[0] for S in range(7)]
assert close(E0_S, torch.stack(E0_N))                               # Theorem 9.9: sector S^z = S has spin S
assert (E0_S[1:] > E0_S[:-1] + 1e-8).all()                          # E_0(S) strictly increasing
gap = {}
for n in (6, 8, 10, 12):
    s = SU2Chain(n, 'periodic', complex_valued=False)
    e = [torch.linalg.eigvalsh(s.heisenberg(spin=S).matrix[0, 0])[0].item() for S in (0, 1)]
    gap[n] = n * (e[1] - e[0])
print({n: round(g, 4) for n, g in gap.items()})
assert 0 < gap[6] < gap[8] < gap[10] < gap[12]                      # > 0: Theorem 9.9; increase: numerical
n8, A_set = 8, (0, 1, 2, 3, 4)                                      # A = {0, ..., 4}, B = {5, 6, 7}
edges = [(0, 5), (1, 5), (1, 6), (2, 6), (3, 6), (4, 7), (0, 7), (0, 1), (2, 3), (5, 6)]   # connected
g8 = torch.Generator().manual_seed(2)
J8 = [(1 if (i in A_set) != (j in A_set) else -1) * (0.5 + torch.rand(1, generator=g8).item()) for i, j in edges]
dot = lambda i, j: 0.5 * (op(n8, {i: sp, j: sm}) + op(n8, {i: sm, j: sp})) + op(n8, {i: sz, j: sz})
H8 = sum(Jij * dot(i, j) for Jij, (i, j) in zip(J8, edges))          # J_ij > 0 between, < 0 within A, B
Sp8, Sz8 = sum(op(n8, {i: sp}) for i in range(n8)), sum(op(n8, {i: sz}) for i in range(n8))
S2 = Sp8.T @ Sp8 + Sz8 @ Sz8 + Sz8                                  # S^2 = S^- S^+ + (S^z)^2 + S^z
n_up = torch.tensor([bin(b).count('1') for b in range(2 ** n8)])
n_B = torch.tensor([bin(b >> 5).count('1') for b in range(2 ** n8)])  # up spins on B
for N in range(4, 9):                                               # sectors S^z = M = N - 4 >= 0
    idx = torch.nonzero(n_up == N).flatten()
    w, U = torch.linalg.eigh(H8[idx][:, idx])
    x, spin = (-1.0) ** n_B[idx] * U[:, 0], max(N - 4, 1)            # expected spin max(M, S_*)
    assert (x * x.sum().sign()).min() > 0 and (len(w) == 1 or w[1] - w[0] > 1e-8)          # Theorem 9.7
    assert abs(U[:, 0] @ S2[idx][:, idx] @ U[:, 0] - spin * (spin + 1)) < 1e-10           # Theorem 9.9(a)
w8 = torch.linalg.eigvalsh(H8)
assert w8[2] - w8[0] < 1e-10 and w8[3] - w8[0] > 1e-8                # ground space: one triplet, S_* = 1
```

Die $`SU(2)`$-Blöcke haben die Größen $`m_S=132,297,275,154,54,11,1`$ gegenüber $`924`$ für den größten $`U(1)`$-Block, und
der Aufwand $`\sum k^3`$ der Diagonalisierung sinkt von $`2.0\cdot10^9`$ auf $`5.3\cdot10^7`$. *Numerische Beobachtung:*
$`L`$ mal die Lücke wächst langsam, $`4.108,4.181,4.232,4.270`$ für $`L=6,\dots,12`$; über ihr Verhalten für $`L\to\infty`$
beweisen wir nichts.

**Korrelationen.** Die Spins der Plätze bilden Vektoroperatoren und ihre sphärischen Komponenten Tensoroperatoren
vom Rang eins (Lemma B.15(a),(b)). Nach dem **Wigner-Eckart-Theorem** (Satz B.16) sind ihre Matrixelemente zwischen
zwei Multipletts ein Clebsch-Gordan-Koeffizient mal ein **reduziertes Matrixelement**, das nicht von den
magnetischen Quantenzahlen abhängt, und $`\mathbf s_i\cdot\mathbf s_j=\sum_q(-1)^qs^{(i)}_qs^{(j)}_{-q}`$ (Lemma B.15(a)).
`spin_correlation` nutzt dies, um $`\mathbf s_i\cdot\mathbf s_j`$ auf jedem $`H_S`$ allein aus Vektoren höchsten Gewichts
zu bauen; die Methode `'direct'` komprimiert stattdessen den Bindungsoperator auf die Vektoren höchsten Gewichts,
was nach Proposition 4.15(b) korrekt ist, da $`\mathbf s_i\cdot\mathbf s_j\in\{S_a\}'`$. Dass die Implementierung der ersten
Methode korrekt ist, wird hier nicht bewiesen; der Code prüft, dass beide untereinander und mit der Definition
übereinstimmen, ausgewertet auf dem im $`U(1)`$-Sektor ausgeschriebenen Grundzustand:

```python
gs = sym.ground_state(H_su2)
S0, coef = gs.spin.item(), gs.coefficients[0, :sym.multiplicity(gs.spin.item())]
pairs = [(0, r) for r in range(12)]
corr_we = sym.expectation(sym.spin_correlation(pairs), S0, coef)
corr_dir = sym.expectation(sym.spin_correlation(pairs, method='direct'), S0, coef)
psi = sym.embed(S0, coef)                                            # the ground state in the sector N = 6
lab, z = sym.u1.labels[6], sym.u1.bits(6).to(f64) - 0.5
def s_dot_s(i, j):                                                   # <psi| s_i . s_j |psi>, from the definition
    if i == j:
        return 0.75
    flip = ((lab >> i) & 1) != ((lab >> j) & 1)
    partner = sym.u1.index_of[lab[flip] ^ ((1 << i) | (1 << j))]
    return ((psi ** 2 * z[:, i] * z[:, j]).sum() + 0.5 * (psi[flip] * psi[partner]).sum()).item()
assert S0 == 0 and close(corr_we, corr_dir)
assert close(corr_we, torch.tensor([s_dot_s(0, r) for r in range(12)], dtype=f64))
assert ((-1.0) ** torch.arange(1, 12) * corr_we[1:] > 0).all()       # signs (-1)^r
print([round(x, 4) for x in corr_we.tolist()])
```

Die Korrelationen $`\langle\mathbf s_0\cdot\mathbf s_r\rangle=0.75,-0.4489,0.1878,-0.1658,0.1209,-0.1224,0.1070,\dots`$
wechseln das Vorzeichen, wie es die Sätze 9.7 und 9.9 für jeden geraden Ring erzwingen ([Vorlesungsnotizen, Abschnitt 9.3](../../lecture_notes/de/09_spin_chains.md#s-9-3)).

<a id="s-9-4"></a>

## 9.4 Exakte Prüfungen: freie Fermionen und der Bethe-Ansatz

*Vorlesungsnotizen:* [Abschnitt 9.4](../../lecture_notes/de/09_spin_chains.md#s-9-4); nach [Satz 9.10](../../lecture_notes/de/09_spin_chains.md#r-9-10) sind die Energien des Sektors $`N`$ der offenen XX-Kette die Summen von $`N`$ verschiedenen $`\varepsilon_k=J\cos\frac{\pi k}{L+1}`$; [Proposition 9.12](../../lecture_notes/de/09_spin_chains.md#r-9-12) (zwei Magnonen, Bethe-Gleichungen).

```python
L = 10
xx = SpinChain(L)
eps = [math.cos(math.pi * k / (L + 1)) for k in range(1, L + 1)]
for N in range(L + 1):
    E = torch.linalg.eigvalsh(xx.xxz(J=1.0, Delta=0.0, sector=N).matrix[0, 0])
    E_ff = sorted(sum(K) for K in itertools.combinations(eps, N))     # sums over N-element subsets
    assert close(E, torch.tensor(E_ff, dtype=f64))                    # all 1024 levels, Theorem 9.10
    assert abs(E[0].item() - sum(sorted(eps)[:N])) < 1e-12           # ground state: the N lowest modes
```

Nach [Proposition 9.12](../../lecture_notes/de/09_spin_chains.md#r-9-12) liefert jede Nullstelle $`q`$ von $`F(q)=(-1)^ne^{iqL}(c-\Delta e^{-iq})+c-\Delta e^{iq}`$, mit
$`K=2\pi n/L`$ und $`c=\cos\frac K2`$, einen Eigenvektor $`\Psi`$ des Sektors $`N=2`$ des Rings
([Vorlesungsnotizen, Abschnitt 9.4](../../lecture_notes/de/09_spin_chains.md#s-9-4)). Der Code sucht Nullstellen $`q\in(0,\pi)`$ (Streuzustände) und auf den Geraden
$`q=i\kappa`$, $`\pi+i\kappa`$, $`\kappa>0`$ (gebundene Zustände), baut $`\Psi`$ zu jeder Nullstelle für $`L=9`$, $`\Delta=0.6`$ und vergleicht
mit dem Sektor $`N=2`$ mit $`36`$ Zuständen:

```python
L, Delta = 9, 0.6
ring9 = SpinChain(L, boundary='periodic')
H2 = ring9.xxz(J=1.0, Delta=Delta, sector=2).matrix[0, 0]
x1, x2 = np.triu_indices(L, 1)                                       # all pairs x1 < x2
rows = ring9.index_of[torch.tensor((1 << x1) | (1 << x2))]
def roots(g, grid):                                                  # sign changes of g on the grid, refined
    v = np.array([g(t) for t in grid])
    return [brentq(g, grid[i], grid[i + 1], xtol=1e-15) for i in np.nonzero(v[:-1] * v[1:] < 0)[0]]
E_bethe, vecs, n_bound = [], [], 0
for n in range(L):
    K, c = 2 * np.pi * n / L, math.cos(np.pi * n / L)
    F = lambda q: (-1) ** n * np.exp(1j * q * L) * (c - Delta * np.exp(-1j * q)) + c - Delta * np.exp(1j * q)
    G = lambda q: (lambda z: z.real + z.imag)(np.exp(-0.5j * q * L) * F(q))    # real, or imaginary
    bound = [s + 1j * t for s in (0, np.pi)                          # q = i kappa and pi + i kappa
             for t in roots(lambda t: F(s + 1j * t).real, np.linspace(1e-3, 20, 4000))]
    n_bound += len(bound)
    for q in roots(G, (np.arange(1, 2000) + 0.5 ** 0.5) * np.pi / 2000) + bound:
        psi = np.exp(0.5j * K * (x1 + x2)) * ((c - Delta * np.exp(1j * q)) * np.exp(1j * q * (x1 - x2))
                                              - (c - Delta * np.exp(-1j * q)) * np.exp(-1j * q * (x1 - x2)))
        v = torch.zeros(36, dtype=c128)
        v[rows] = torch.tensor(psi / np.linalg.norm(psi))
        E = 9 * Delta / 4 + 2 * (c * np.cos(q) - Delta)                 # E_F + J (cos k1 + cos k2 - 2 Delta)
        assert abs(E.imag) < 1e-12 and close(H2 @ v, E.real * v, 1e-9)   # Proposition 9.12
        E_bethe.append(E.real)
        vecs.append(v)
assert torch.linalg.matrix_rank(torch.stack(vecs)) == 36             # a complete eigenbasis
assert close(torch.sort(torch.tensor(E_bethe)).values, torch.linalg.eigvalsh(H2))
print(len(E_bethe), n_bound)
ring8, Phi = SpinChain(8, 'periodic'), torch.zeros(28, dtype=c128)     # even L: Phi = sum_x (-1)^x |x, x+1>
for x in range(8):
    Phi[ring8.index_of[(1 << x) | (1 << (x + 1) % 8)]] = (-1) ** x
assert close(ring8.xxz(1.0, Delta, sector=2).matrix[0, 0] @ Phi, (8 * Delta / 4 - Delta) * Phi)   # E_F - J Delta
```

*Numerische Beobachtung:* Für $`L=9`$ und $`\Delta=0.6`$ findet die Suche $`36`$ Nullstellen, $`34`$ Streuzustände und $`2`$
gebundene Zustände, deren Vektoren linear unabhängige Eigenvektoren sind (bis auf Rundung): Alle $`36`$ Eigenzustände
stammen aus den Bethe-Gleichungen. Die Vollständigkeit ist für andere $`L`$ nicht bewiesen, und für gerades $`L`$ und
$`\Delta\ne0`$ ist sie verletzt ([Vorlesungsnotizen, Abschnitt 9.4](../../lecture_notes/de/09_spin_chains.md#s-9-4)): Der Vektor $`\Phi=\sum_x(-1)^x|x,x+1\rangle`$ ist zu allen
Bethe-Vektoren orthogonal, und die letzten Zeilen des Codes prüfen $`H\Phi=(E_F-J\Delta)\Phi`$ für $`L=8`$.

<a id="s-9-5"></a>

## 9.5 Lanczos- und Krylov-Verfahren

*Vorlesungsnotizen:* [Abschnitt 9.5](../../lecture_notes/de/09_spin_chains.md#s-9-5); [Proposition 9.13](../../lecture_notes/de/09_spin_chains.md#r-9-13) (Lanczos: $`T_m=Q_m^*HQ_m`$ ist tridiagonal), [Satz 9.14](../../lecture_notes/de/09_spin_chains.md#r-9-14) (die Ritz-Werte, die Eigenwerte von $`T_m`$, erfüllen $`\theta^{(m)}_j\ge\lambda_j`$ und fallen mit $`m`$), [Proposition 9.15](../../lecture_notes/de/09_spin_chains.md#r-9-15) (Krylov-Zeitentwicklung).

Der Sektor $`N=L/2`$ hat die Dimension $`12870`$ für $`L=16`$ und $`184756`$ für $`L=20`$, aber $`H`$ ist dünn besetzt ($`O(L)`$
von null verschiedene Einträge pro Zeile, `xxz_sparse`), und Krylov-Verfahren brauchen nur Produkte $`Hv`$. Sei $`H`$
hermitesch von der Größe $`n`$ und $`v\ne0`$. Der **Krylov-Raum** ist
$`\mathcal K_m=\mathrm{span}\{v,Hv,\dots,H^{m-1}v\}`$.

Satz 9.14 sagt nichts darüber, wie schnell sich $`\theta^{(m)}_1`$ an $`\lambda_1`$ annähert; wir beobachten es unten nur.
In Gleitkommaarithmetik verlieren die $`q_j`$ ihre Orthogonalität; `krylov.lanczos` reorthogonalisiert vollständig,
und der Code prüft, dass die berechneten $`q_j`$ orthonormal sind und $`Q_m^*HQ_m=T_m`$ bis auf Rundung gilt. Für die
offene XX-Kette ist die exakte Grundzustandsenergie des Sektors $`N`$ die Summe der $`N`$ kleinsten $`\varepsilon_k`$
(Satz 9.10(c)):

```python
L, N = 12, 6
xx = SpinChain(L, complex_valued=False)
Hs = xx.xxz_sparse(J=1.0, Delta=0.0, sector=N)                      # sparse: diagonal + hopping
alpha, beta, Q = krylov.lanczos(Hs.matvec, torch.randn(1, xx.sector_dim(N), dtype=f64), 40)
T = torch.diag(alpha[0]) + torch.diag(beta[0], 1) + torch.diag(beta[0], -1)
assert close(Q[0] @ Q[0].T, torch.eye(40, dtype=f64)) and close(Q[0] @ Hs.to_dense()[0] @ Q[0].T, T)
ritz = torch.tensor([torch.linalg.eigvalsh(T[:m, :m])[0].item() for m in range(1, 41)])   # T_m = leading block
E_xx = sum(sorted(math.cos(math.pi * k / (L + 1)) for k in range(1, L + 1))[:N])
assert (ritz[1:] <= ritz[:-1] + 1e-12).all() and (ritz >= E_xx - 1e-12).all()             # Theorem 9.14
print([f"{(ritz[m - 1] - E_xx).item():.1e}" for m in (5, 10, 20, 40)])
```

*Numerische Beobachtung:* Der Fehler des kleinsten Ritz-Werts fällt monoton (wie er nach Satz 9.14 muss), $`0.42`$,
$`3.1\cdot10^{-3}`$ und $`3.1\cdot10^{-6}`$ nach $`5`$, $`10`$ und $`20`$ Schritten, und erreicht nach $`40`$ Schritten das
Rundungsniveau, in Dimension $`924`$. `krylov.ground_state` startet Lanczos neu (von einem gaußschen Zufallsvektor,
sofern kein `v0` angegeben ist), bis das Residuum $`\|H\psi-E\psi\|`$ seines normierten Ritz-Paares $`(E,\psi)`$ unter
$`10^{-10}\max(1,|E|)`$ liegt (höchstens $`20`$ Neustarts, sonst mit einer Warnung), und gibt das Residuum zurück. Ein
kleines Residuum garantiert einen nahen Eigenwert, aber nicht, dass er der kleinste ist
([Vorlesungsnotizen, Abschnitt 9.5](../../lecture_notes/de/09_spin_chains.md#s-9-5)); daher vergleichen wir mit dichter
Diagonalisierung, unserer einzigen Referenz (hier berechnet, die Werte sind also numerisch): der Sektor $`N=L/2`$ bis
$`L=12`$ und für $`L=16`$ der Impulsblock $`k=0`$ der Dimension $`810`$, der nach Korollar 9.8 den Grundzustand des Sektors
$`N=8`$ enthält. Für $`L=4`$ ist der Wert $`-2J`$ exakt ([Vorlesungsnotizen, Abschnitt 9.5](../../lecture_notes/de/09_spin_chains.md#s-9-5)).

```python
E_dense = {}                                                         # numerical references, dense
for n in (4, 6, 8, 10, 12):
    E_dense[n] = torch.linalg.eigvalsh(SpinChain(n, 'periodic').heisenberg(sector=n // 2).matrix[0, 0])[0].item()
    print(n, round(E_dense[n] / n, 6))
assert abs(E_dense[4] + 2) < 1e-12                                   # exact value for L = 4
ring16 = SpinChain(16, 'periodic')
k0 = ring16.momentum_dims(8)[0]                                      # dimension of the block k = 0
H16 = ring16.xxz_momentum(J=1.0, Delta=1.0, sector=8).matrix[0, 0, :k0, :k0]
E_dense[16] = torch.linalg.eigvalsh(H16)[0].item()                   # ground state of N = 8 (Corollary 9.8)
t0 = time.time()
E16, psi16, res = krylov.ground_state(SpinChain(16, 'periodic', complex_valued=False).xxz_sparse(sector=8),
                                      generator=torch.Generator().manual_seed(0))
print(f"L = 16: E0 = {E16.item():.10f}, residual {res.item():.1e}, {time.time() - t0:.2f} s; k = 0 block: {k0}")
assert abs(E16.item() - E_dense[16]) < 1e-9 and res.item() < 1e-9
```

Für $`L=16`$ (Dimension $`12870`$) braucht Lanczos auf einem CPU-Thread deutlich weniger als eine Sekunde. *Numerische
Beobachtung:* Die Energie pro Platz, $`-0.5,-0.467129,-0.456387,-0.451545,-0.448949`$ für $`L=4,\dots,12`$, wächst mit $`L`$;
über ihren Grenzwert beweisen wir nichts.

**Zeitentwicklung.** Für $`e^{-i\tau H}v`$ projiziert man auf $`\mathcal K_m`$: $`w_\tau=\|v\|\,Q_me^{-i\tau T_m}e_1`$.

Nach [Proposition 9.15](../../lecture_notes/de/09_spin_chains.md#r-9-15) gilt $`\|w_\tau\|=\|v\|`$ und $`\langle w_\tau,Hw_\tau\rangle=\langle v,Hv\rangle`$ für jedes $`\tau`$, und
$`w_\tau=e^{-i\tau H}v`$, falls $`\mathcal K_m`$ $`H`$-invariant ist.

Jeder Krylov-Schritt erhält also Norm und Energie, ob er genau ist oder nicht; Proposition 9.15 sagt nichts
über die Genauigkeit, außer wenn $`\mathcal K_m`$ invariant ist. `krylov.evolve` unterteilt das Zeitintervall in kurze
Schritte, jeden in einem Krylov-Raum mit $`m=30`$, und warnt, wenn ein heuristischer Fehlerindikator groß ist; wir
beweisen dafür keine Fehlerschranke und vergleichen stattdessen mit $`e^{-itH}`$, berechnet aus einer dichten
Eigenzerlegung (Satz 1.5). Wir führen einen
Quench aus dem Néel-Zustand $`|{\uparrow\downarrow\uparrow\downarrow}\cdots\rangle`$ auf dem Heisenberg-Ring durch und
verfolgen die alternierende Magnetisierung $`m_s=\frac2L\sum_i(-1)^i\langle s^z_i\rangle`$, im Vergleich mit dem dichten
Propagator für $`L=12`$:

```python
times, m_s = [0.0, 0.5, 1.0, 2.0, 4.0], {}
stag = lambda c, N, out: [((p[0].abs() ** 2) @ (c.bits(N).to(f64) - 0.5) @ torch.tensor(
    [(-1.0) ** i for i in range(c.L)], dtype=f64) * 2 / c.L).item() for p in out]
for n in (12, 16):
    c = SpinChain(n, 'periodic')
    Hn = c.xxz_sparse(J=1.0, Delta=1.0, sector=n // 2)
    psi0 = c.vector_in_sector('10' * (n // 2))[None]
    out = krylov.evolve(Hn, psi0, times)                               # (times, batch, dim)
    energy = torch.stack([(p.conj() * Hn.matvec(p)).sum(-1).real for p in out])
    assert close(out.norm(dim=-1), 1.0) and close(energy, energy[0])   # Proposition 9.15
    m_s[n] = stag(c, n // 2, out)
    print(n, [round(x, 4) for x in m_s[n]])
    if n == 12:
        w, U = torch.linalg.eigh(c.xxz(1.0, 1.0, sector=6).matrix[0, 0])
        for t, p in zip(times, out):
            assert close(p[0], U @ (torch.exp(-1j * t * w) * (U.mH @ psi0[0])))
print(f"{max(abs(a - b) for a, b in zip(m_s[12][:4], m_s[16][:4])):.1e}")   # L = 12 vs 16 up to t = 2
```

*Numerische Beobachtung:* $`m_s`$ fällt von $`1`$ über $`0.770`$ und $`0.279`$ auf $`-0.183`$ bei $`t=2`$, für beide Größen bis
auf $`1.7\cdot10^{-5}`$ übereinstimmend; bei $`t=4`$ unterscheiden sie sich ($`0.026`$ gegenüber $`0.035`$).

<a id="s-9-6"></a>

## 9.6 Unordnung: Niveaustatistik und Verschränkung

*Vorlesungsnotizen:* [Abschnitt 9.6](../../lecture_notes/de/09_spin_chains.md#s-9-6); [Satz 9.16](../../lecture_notes/de/09_spin_chains.md#r-9-16) (Peschels Formel).

Die **Heisenberg-Kette im Zufallsfeld** hat unabhängige, auf $`[-W,W]`$ gleichverteilte Felder $`h_i`$. Sie erhält
$`S^z`$ (Proposition 9.1), daher wird die Niveaustatistik pro Sektor $`N`$ berechnet (Satz 5.15 und Abschnitt 5.4: Das
Mischen von Sektoren täuscht unabhängige Niveaus vor). Unterscheiden sich zwei Felder, so ist
$`THT^*-H=\sum_i(h_{i-1}-h_i)s^z_i\ne0`$, Translationen sind also keine Symmetrie mehr. Wir berechnen $`\langle r\rangle`$ in
der mittleren Hälfte des Sektors $`N=L/2`$, mit einem Standardfehler über die unabhängigen Realisierungen (jede
liefert ein mittleres Verhältnis), und vergleichen mit zwei Referenzen: $`2\log2-1\approx0.3863`$, dem exakten
mittleren Verhältnis unabhängiger exponentialverteilter Abstände (Proposition 5.14(a)), und dem mittleren
Verhältnis von GOE-Matrizen (wie in Abschnitt 5.3 definiert; der Code zieht $`G+G^T`$, das Doppelte einer
GOE-Matrix, was kein Verhältnis ändert) derselben Größe $`252`$ wie der Sektor für $`L=10`$, hier gezogen. Die Vergleiche
sind numerische Beobachtungen, keine Sätze: Über das Spektrum der zufälligen Kette beweisen wir nichts.

```python
gen = torch.Generator().manual_seed(1)
r_mean = {}
for n, R in ((8, 200), (10, 100)):
    c = SpinChain(n, 'periodic', complex_valued=False)
    k = c.sector_dim(n // 2)
    for W in (0.5, 1, 2, 4, 8, 16):
        E = torch.linalg.eigvalsh(c.random_field_heisenberg(W, R, sector=n // 2, generator=gen).matrix[:, 0])
        r = SpinChain.level_spacing_ratio(E[:, k // 4: 3 * k // 4])     # one value per realisation
        r_mean[n, W] = (r.mean().item(), (r.std() / math.sqrt(R)).item())
    print(n, '  '.join(f"W={W}: {m:.3f}+-{s:.3f}" for (n_, W), (m, s) in r_mean.items() if n_ == n))
X = torch.randn(100, k, k, dtype=f64, generator=gen)                  # 100 GOE matrices of size k = 252
r_goe = SpinChain.level_spacing_ratio(torch.linalg.eigvalsh(X + X.mT)[:, k // 4: 3 * k // 4])
print(f"GOE({k}): {r_goe.mean().item():.3f}+-{(r_goe.std() / 10).item():.3f}, Poisson: {2 * math.log(2) - 1:.4f}")
assert abs(r_mean[10, 0.5][0] - r_goe.mean().item()) < 0.015               # numerical observation
assert abs(r_mean[10, 16][0] - (2 * math.log(2) - 1)) < 0.015               # numerical observation
assert r_mean[10, 1][0] > r_mean[10, 2][0] > r_mean[10, 4][0] > r_mean[10, 8][0]
```

*Numerische Beobachtung:* Für $`L=10`$ geht das mittlere Verhältnis von $`0.530`$ ($`W=0.5`$), nahe am gezogenen GOE-Wert
$`0.533`$, über $`0.484`$ ($`W=2`$) und $`0.410`$ ($`W=4`$) bis $`0.384`$ ($`W=16`$), nahe an $`2\log2-1`$, mit Standardfehlern von
$`0.003`$–$`0.004`$. Für $`L=8`$ liegt der Wert bei $`W=0.5`$ ($`0.507`$) unter dem bei $`W=1`$ ($`0.525`$). Bei $`W=0`$ kommutiert
der Ring mit $`T`$ (Abschnitt 9.2) und mit dem Gesamtspin (Abschnitt 9.3), also zerfällt sein Sektor $`N=L/2`$ in
Blöcke; die Felder koppeln diese Blöcke nur über $`\sum_ih_is^z_i`$, mit Operatornorm höchstens $`\frac{WL}2`$.

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Physiker deuten die Änderung von
> $`\langle r\rangle`$ als Übergang von Niveaus, die sich wie die des GOE abstoßen (schwache Unordnung), zu Niveaus, die
> sich wie unabhängige verhalten, zurückgeführt auf eine **Vielteilchenlokalisierung** der Eigenzustände (starke
> Unordnung). Der Übergangsbereich driftet mit $`L`$; wo er für $`L\to\infty`$ liegt und ob es eine lokalisierte Phase
> überhaupt gibt, sind offene Fragen in der physikalischen Literatur (D. Luitz, N. Laflorencie, F. Alet,
> *Phys. Rev. B* 91 (2015) 081103; J. Šuntajs, J. Bonča, T. Prosen, L. Vidmar, *Phys. Rev. E* 102 (2020) 062144). Die
> Delle bei $`L=8`$, $`W=0.5`$ wird als Überrest der zusätzlichen Blockstruktur bei $`W=0`$ gedeutet (vergleiche
> Abschnitt 5.4).

**Verschränkung freier Fermionen.** Für die Eigenzustände der offenen XX-Kette ist die Verschränkungsentropie der
ersten $`n`$ Plätze durch eine $`n\times n`$-Matrix bestimmt.

`entanglement_entropy` berechnet den reduzierten Zustand Sektor für Sektor (er kommutiert mit der Magnetisierung
des Blocks, Proposition 7.8) und seine Entropie. Der Grundzustand des Sektors $`N`$ ist $`\Phi`$ mit den $`N`$
niedrigsten Moden, einfach, weil die $`\varepsilon_k`$ verschieden sind:

```python
L, N = 12, 6
xx = SpinChain(L)
w, V = torch.linalg.eigh(xx.xxz(J=1.0, Delta=0.0, sector=N).matrix[0, 0])
rho = DensityMatrix(xx.sector_algebra(N), matrix=torch.outer(V[:, 0], V[:, 0].conj())[None, None])
eps_k, u = np.linalg.eigh(np.eye(L, k=1) / 2 + np.eye(L, k=-1) / 2)  # h = A_L / 2, ascending
C = u[:, :N] @ u[:, :N].T                                            # correlation matrix of the N lowest modes
bern = lambda nu: -(nu * np.log(nu) + (1 - nu) * np.log(1 - nu))
S_peschel = [bern(np.clip(np.linalg.eigvalsh(C[:n, :n]), 1e-15, 1 - 1e-15)).sum() for n in range(1, L)]
S_lib = [xx.entanglement_entropy(rho, n, sector=N).item() for n in range(1, L)]
assert np.allclose(S_lib, S_peschel, atol=1e-9)                       # Theorem 9.16
print([round(s, 4) for s in S_lib[:6]])
```

*Numerische Beobachtung:* Die Entropien $`0.6931,0.5199,0.7456,0.6283,0.7719,0.6560`$ für $`n=1,\dots,6`$ oszillieren
mit der Parität von $`n`$. Satz 9.16 stammt von I. Peschel (*J. Phys. A* 36 (2003) L205).

<a id="s-9-7"></a>

## 9.7 Womit man vergleichen sollte

QuSpin (Python) baut Spinketten-Hamiltonoperatoren in Sektoren von Magnetisierung, Impuls, Parität und Spinumkehr
und verwendet die dünn besetzten Eigenwertlöser und Matrixexponentiale von SciPy; spezialisierte Codes für exakte
Diagonalisierung (XDiag, SPINPACK, HΦ) sind für große dünn besetzte Probleme optimiert, und mit Gittersymmetrien und
großen Rechnern haben solche Codes etwa 40–50 Spins erreicht; ITensor, TeNPy und TensorKit.jl (mit
$`SU(2)`$-Fusionsbäumen) behandeln lange Ketten mit Matrixproduktzuständen. Was die Bibliothek hinzufügt, ist
bescheiden und spezifisch: die Kette als `TypeIAlgebra` $`\bigoplus_NM_{\binom LN}`$, sodass die Zustände, Kanäle,
Entropien und Dynamiken der Kapitel 3–8 Sektor für Sektor gelten; den Inklusionsturm mit seinem Index;
$`SU(2)`$-Blöcke aus Vektoren höchsten Gewichts mit Wigner-Eckart-Korrelationen; Stapel von
Unordnungsrealisierungen in einem einzigen Aufruf; und reduzierte Zustände, die die Sektoren nie verlassen. Sie hat
keine Spiegelungs- oder Spinumkehrsektoren, Impulssektoren nur für homogene XXZ-Ketten, dichte Verfahren bis
$`L\approx16`$, Krylov-Verfahren bis $`L\approx24`$ (etwa $`26`$ auf einer großen GPU) und keine Tensornetzwerke.

<a id="x-9"></a>

## Aufgaben

[Aufgaben 1–3](../../lecture_notes/de/09_spin_chains.md#x-9) der Vorlesungsnotizen: die Dimensionen $`d_{2,m}=2,1,2,1`$ der Impulssektoren für
$`L=4`$, $`N=2`$; die Summenregeln $`E_0=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle`$ und $`\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle=0`$ für den
Grundzustand des Heisenberg-Rings; und `krylov.ground_state`, gestartet vom Vektor mit allen Koeffizienten gleich $`1`$
im Sektor $`N=L/2`$, das $`E=JL/4`$ liefert, das obere Ende des Spektrums. Man überprüfe dies numerisch.

<details>
<summary>Lösungen</summary>

```python
assert SpinChain(4, 'periodic').momentum_dims(2) == [2, 1, 2, 1]
```

```python
assert abs(12 * corr_we[1].item() - gs.energy.item()) < 1e-10 and abs(corr_we.sum().item()) < 1e-10
```

```python
c10 = SpinChain(10, 'periodic', complex_valued=False)
E_sym, _, _ = krylov.ground_state(c10.xxz_sparse(1.0, 1.0, sector=5), v0=torch.ones(1, 252, dtype=f64))
assert abs(E_sym.item() - 2.5) < 1e-10 and E_sym.item() > E_dense[10] + 7  # far from E_0
```

</details>

**Weiter:** [Tutorial 10. Quantencomputing](10_quantum_computing.md).
