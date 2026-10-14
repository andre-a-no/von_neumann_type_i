[English](../09_spin_chains.md) · [Español](../es/09_spin_chains.md) · **Français** · [Deutsch](../de/09_spin_chains.md) · [中文](../zh/09_spin_chains.md) · [日本語](../ja/09_spin_chains.md) · [Русский](../ru/09_spin_chains.md)

# Tutoriel 9. Chaînes de spins

Ce tutoriel accompagne le [chapitre 9](../../lecture_notes/fr/09_spin_chains.md) des notes de cours, qui construit les secteurs de la chaîne
XXZ — aimantation ([proposition 9.1](../../lecture_notes/fr/09_spin_chains.md#r-9-1)), impulsion ([proposition 9.5](../../lecture_notes/fr/09_spin_chains.md#r-9-5)) et spin total — et démontre ce
qu'ils disent des états fondamentaux ([théorème 9.7](../../lecture_notes/fr/09_spin_chains.md#r-9-7), [théorème 9.9](../../lecture_notes/fr/09_spin_chains.md#r-9-9)). Nous construisons les chaînes avec
`SpinChain` et `SU2Chain`, confrontons la bibliothèque aux fermions libres et à l'ansatz de Bethe, utilisons les
méthodes de Lanczos et de Krylov de `krylov` pour les secteurs trop grands pour des matrices denses, et calculons
des statistiques de niveaux et l'intrication.

Dans tout ce chapitre, les sites sont $0,\dots,L-1$, la base de $(\mathbb C^2)^{\otimes L}$ est formée des
configurations $|b\rangle=|b_0\dots b_{L-1}\rangle$ avec $b_i=1$ pour un spin vers le haut (convention de la
bibliothèque : le bit $i$ d'une étiquette est le site $i$), $s^a_i$ ($a=x,y,z$) sont les matrices de spin
$\frac12\sigma_a$ sur le site $i$, $s^\pm_i=s^x_i\pm is^y_i$, $n_i=s^z_i+\frac12$, $S^a=\sum_is^a_i$, et
$N=S^z+\frac L2$ est le nombre de spins vers le haut. Une chaîne est **ouverte** (liaisons $(i,i+1)$, $i<L-1$) ou
**périodique** (avec aussi la liaison $(L-1,0)$ ; alors $L\ge3$ et les indices sont pris modulo $L$). $J>0$
(antiferromagnétique), $\hbar=1$, et $\log$ est le logarithme népérien.

<a id="s-9-1"></a>

## 9.1 La chaîne XXZ et ses secteurs $U(1)$

*Notes de cours :* [section 9.1](../../lecture_notes/fr/09_spin_chains.md#s-9-1) ; [proposition 9.1](../../lecture_notes/fr/09_spin_chains.md#r-9-1) ($H$ conserve $S^z$), [proposition 9.2](../../lecture_notes/fr/09_spin_chains.md#r-9-2), [lemme 9.3](../../lecture_notes/fr/09_spin_chains.md#r-9-3) et [corollaire 9.4](../../lecture_notes/fr/09_spin_chains.md#r-9-4) (la tour des chaînes et son indice).

`SpinChain(L).xxz` construit les blocs directement à partir des configurations de chaque secteur
(`chain.labels[N]`). Nous les comparons à la matrice complète $2^L\times2^L$ construite par produits de
Kronecker :

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

La tour $\mathcal A_0\subset\mathcal A_1\subset\cdots$ des chaînes : `inclusions.chain_tower` construit les inclusions,
`inclusions.bratteli_diagram` leur diagramme de Bratteli, et `Inclusion.index` l'indice $\|\Lambda\|^2$. Le code vérifie
la règle de Pascal, le [lemme 9.3](../../lecture_notes/fr/09_spin_chains.md#r-9-3) et $[\mathcal A_{L+1}:\mathcal A_L]=4\cos^2\frac\pi{2L+4}$ ([corollaire 9.4](../../lecture_notes/fr/09_spin_chains.md#r-9-4)).

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

## 9.2 Translations, impulsion et règle des signes de Marshall

*Notes de cours :* [section 9.2](../../lecture_notes/fr/09_spin_chains.md#s-9-2) ; [proposition 9.5](../../lecture_notes/fr/09_spin_chains.md#r-9-5) (secteurs d'impulsion), [théorème 9.7](../../lecture_notes/fr/09_spin_chains.md#r-9-7) (règle des signes de Marshall), [corollaire 9.8](../../lecture_notes/fr/09_spin_chains.md#r-9-8) (impulsion de l'état fondamental).

`xxz_momentum` construit les blocs dans la base normalisée $|a,k\rangle$, avec un représentant $a$ par orbite. Nous
dénombrons les orbites indépendamment et comparons les spectres avec le secteur dense :

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

Le secteur $N=6$, de $924$ états, se découpe en blocs de tailles $80,75,78,76,78,75,80,75,78,76,78,75$.

D'après la règle des signes de Marshall ([théorème 9.7](../../lecture_notes/fr/09_spin_chains.md#r-9-7)), pour des chaînes ouvertes et périodiques en champs aléatoires,
$(-1)^{N_B(b)}\psi(b)$ a un seul signe dans chaque secteur, et d'après le [corollaire 9.8](../../lecture_notes/fr/09_spin_chains.md#r-9-8) le niveau le plus bas du secteur $N$ de
l'anneau se trouve dans le bloc d'impulsion $k=0$ pour $N$ pair et $k=\pi$ pour $N$ impair :

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

## 9.3 $SU(2)$ en $\Delta=1$

*Notes de cours :* [section 9.3](../../lecture_notes/fr/09_spin_chains.md#s-9-3) ; [théorème 9.9](../../lecture_notes/fr/09_spin_chains.md#r-9-9) (Lieb–Mattis). Pour $\Delta=1$ et $h=0$, $H\cong\bigoplus_S1_{2S+1}\otimes H_S$, et `SU2Chain` construit les blocs $H_S$ sur les vecteurs de plus haut poids.

Nous comparons les tailles des blocs pour douze spins, vérifions le théorème 9.9 à l'aide des secteurs $U(1)$,
calculons le gap et vérifions les théorèmes 9.7 et 9.9 sur un graphe biparti qui n'est pas une chaîne ($n_A=5$,
$n_B=3$, $S_*=1$, couplages aléatoires) :

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

Les blocs de $SU(2)$ ont pour tailles $m_S=132,297,275,154,54,11,1$, contre $924$ pour le plus grand bloc de
$U(1)$, et le coût $\sum k^3$ de la diagonalisation tombe de $2.0\cdot10^9$ à $5.3\cdot10^7$. *Observation numérique :*
$L$ fois le gap croît lentement, $4.108,4.181,4.232,4.270$ pour $L=6,\dots,12$ ; nous ne démontrons rien sur son
comportement quand $L\to\infty$.

**Corrélations.** Les spins des sites forment des opérateurs vectoriels, et leurs composantes sphériques des
opérateurs tensoriels de rang un (lemme B.15(a),(b)). D'après le **théorème de Wigner–Eckart** (théorème B.16), leurs
éléments de matrice entre deux multiplets sont un coefficient de Clebsch–Gordan multiplié par un **élément de
matrice réduit** qui ne dépend pas des nombres quantiques magnétiques, et
$\mathbf s_i\cdot\mathbf s_j=\sum_q(-1)^qs^{(i)}_qs^{(j)}_{-q}$ (lemme B.15(a)). `spin_correlation` s'en sert pour construire
$\mathbf s_i\cdot\mathbf s_j$ sur chaque $H_S$ à partir des seuls vecteurs de plus haut poids ; la méthode `'direct'`
comprime à la place l'opérateur de liaison sur les vecteurs de plus haut poids, ce qui est correct d'après la
proposition 4.15(b) puisque $\mathbf s_i\cdot\mathbf s_j\in\{S_a\}'$. La correction de l'implémentation de la première
méthode n'est pas démontrée ici ; le code vérifie que les deux concordent entre elles et avec la définition,
évaluées sur l'état fondamental écrit dans le secteur $U(1)$ :

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

Les corrélations $\langle\mathbf s_0\cdot\mathbf s_r\rangle=0.75,-0.4489,0.1878,-0.1658,0.1209,-0.1224,0.1070,\dots$
alternent en signe, comme l'imposent les théorèmes 9.7 et 9.9 pour tout anneau pair ([notes de cours, section 9.3](../../lecture_notes/fr/09_spin_chains.md#s-9-3)).

<a id="s-9-4"></a>

## 9.4 Vérifications exactes : fermions libres et ansatz de Bethe

*Notes de cours :* [section 9.4](../../lecture_notes/fr/09_spin_chains.md#s-9-4) ; d'après le [théorème 9.10](../../lecture_notes/fr/09_spin_chains.md#r-9-10), les énergies du secteur $N$ de la chaîne XX ouverte sont les sommes de $N$ valeurs distinctes $\varepsilon_k=J\cos\frac{\pi k}{L+1}$ ; [proposition 9.12](../../lecture_notes/fr/09_spin_chains.md#r-9-12) (deux magnons, équations de Bethe).

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

D'après la [proposition 9.12](../../lecture_notes/fr/09_spin_chains.md#r-9-12), chaque racine $q$ de $F(q)=(-1)^ne^{iqL}(c-\Delta e^{-iq})+c-\Delta e^{iq}$, avec
$K=2\pi n/L$ et $c=\cos\frac K2$, donne un vecteur propre $\Psi$ du secteur $N=2$ de l'anneau
([notes de cours, section 9.4](../../lecture_notes/fr/09_spin_chains.md#s-9-4)). Le code cherche des racines $q\in(0,\pi)$ (états de diffusion) et sur les droites
$q=i\kappa$, $\pi+i\kappa$, $\kappa>0$ (états liés), construit $\Psi$ pour chaque racine pour $L=9$, $\Delta=0.6$, et compare
avec le secteur $N=2$ de $36$ états :

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

*Observation numérique :* pour $L=9$ et $\Delta=0.6$, la recherche trouve $36$ racines, $34$ états de diffusion et
$2$ états liés, dont les vecteurs sont des vecteurs propres linéairement indépendants (aux erreurs d'arrondi près) :
les $36$ états propres proviennent tous des équations de Bethe. La complétude n'est pas démontrée pour d'autres
$L$, et pour $L$ pair et $\Delta\ne0$ elle est en défaut ([notes de cours, section 9.4](../../lecture_notes/fr/09_spin_chains.md#s-9-4)) : le vecteur
$\Phi=\sum_x(-1)^x|x,x+1\rangle$ est orthogonal à tous les vecteurs de Bethe, et les dernières lignes du code vérifient
$H\Phi=(E_F-J\Delta)\Phi$ pour $L=8$.

<a id="s-9-5"></a>

## 9.5 Lanczos et méthodes de Krylov

*Notes de cours :* [section 9.5](../../lecture_notes/fr/09_spin_chains.md#s-9-5) ; [proposition 9.13](../../lecture_notes/fr/09_spin_chains.md#r-9-13) (Lanczos : $T_m=Q_m^*HQ_m$ est tridiagonale), [théorème 9.14](../../lecture_notes/fr/09_spin_chains.md#r-9-14) (les valeurs de Ritz, valeurs propres de $T_m$, vérifient $\theta^{(m)}_j\ge\lambda_j$ et décroissent avec $m$), [proposition 9.15](../../lecture_notes/fr/09_spin_chains.md#r-9-15) (évolution temporelle de Krylov).

Le secteur $N=L/2$ est de dimension $12870$ pour $L=16$ et $184756$ pour $L=20$, mais $H$ est creuse ($O(L)$
coefficients non nuls par ligne, `xxz_sparse`), et les méthodes de Krylov n'ont besoin que de produits $Hv$. Soient
$H$ hermitienne de taille $n$ et $v\ne0$. L'**espace de Krylov** est
$\mathcal K_m=\mathrm{span}\{v,Hv,\dots,H^{m-1}v\}$.

Le théorème 9.14 ne dit rien de la vitesse à laquelle $\theta^{(m)}_1$ s'approche de $\lambda_1$ ; nous l'observons
seulement ci-dessous. En virgule flottante, les $q_j$ perdent leur orthogonalité ; `krylov.lanczos` réorthogonalise
complètement, et le code vérifie que les $q_j$ calculés sont orthonormés et que $Q_m^*HQ_m=T_m$ aux erreurs
d'arrondi près. Pour la chaîne XX ouverte, l'énergie exacte de l'état fondamental du secteur $N$ est la somme des
$N$ plus petits $\varepsilon_k$ (théorème 9.10(c)) :

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

*Observation numérique :* l'erreur de la plus petite valeur de Ritz décroît de façon monotone (comme il se doit,
d'après le théorème 9.14), $0.42$, $3.1\cdot10^{-3}$ et $3.1\cdot10^{-6}$ après $5$, $10$ et $20$ pas, et atteint le
niveau des erreurs d'arrondi après $40$, en dimension $924$. `krylov.ground_state` redémarre Lanczos (à partir d'un
vecteur aléatoire gaussien, sauf si `v0` est donné) jusqu'à ce que le résidu $\|H\psi-E\psi\|$ de sa paire de Ritz
normalisée $(E,\psi)$ soit inférieur à $10^{-10}\max(1,|E|)$ (au plus $20$ redémarrages, avec un avertissement
sinon), et renvoie le résidu. Un petit résidu certifie une valeur propre proche, mais pas qu'elle soit la plus petite
([notes de cours, section 9.5](../../lecture_notes/fr/09_spin_chains.md#s-9-5)) ; nous comparons
donc avec la diagonalisation dense, qui est notre seule référence (calculée ici, donc les valeurs sont
numériques) : le secteur $N=L/2$ jusqu'à $L=12$ et, pour $L=16$, le bloc d'impulsion $k=0$ de dimension $810$, qui
contient l'état fondamental du secteur $N=8$ d'après le corollaire 9.8. Pour $L=4$, la valeur $-2J$ est exacte
([notes de cours, section 9.5](../../lecture_notes/fr/09_spin_chains.md#s-9-5)).

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

Pour $L=16$ (dimension $12870$), Lanczos prend nettement moins d'une seconde sur un seul thread CPU. *Observation
numérique :* l'énergie par site, $-0.5,-0.467129,-0.456387,-0.451545,-0.448949$ pour $L=4,\dots,12$, croît avec
$L$ ; nous ne démontrons rien sur sa limite.

**Évolution temporelle.** Pour $e^{-i\tau H}v$, on projette sur $\mathcal K_m$ : $w_\tau=\|v\|\,Q_me^{-i\tau T_m}e_1$.

D'après la [proposition 9.15](../../lecture_notes/fr/09_spin_chains.md#r-9-15), $\|w_\tau\|=\|v\|$ et $\langle w_\tau,Hw_\tau\rangle=\langle v,Hv\rangle$ pour tout $\tau$, et
$w_\tau=e^{-i\tau H}v$ si $\mathcal K_m$ est $H$-invariant.

Ainsi chaque pas de Krylov conserve la norme et l'énergie, qu'il soit précis ou non ; la proposition 9.15 ne dit
rien de la précision, sauf si $\mathcal K_m$ est invariant. `krylov.evolve` subdivise l'intervalle de temps en pas
courts, chacun effectué dans un espace de Krylov avec $m=30$, et émet un avertissement lorsqu'un indicateur
d'erreur heuristique est grand ; nous ne démontrons aucune borne d'erreur pour lui et le comparons plutôt à
$e^{-itH}$ calculée à partir d'une décomposition spectrale dense (théorème 1.5). Nous effectuons une
trempe (*quench*) à partir de l'état de Néel $|{\uparrow\downarrow\uparrow\downarrow}\cdots\rangle$ sur l'anneau de
Heisenberg et suivons l'aimantation alternée $m_s=\frac2L\sum_i(-1)^i\langle s^z_i\rangle$, en comparant avec le
propagateur dense pour $L=12$ :

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

*Observation numérique :* $m_s$ décroît de $1$ à $-0.183$ en $t=2$, en passant par $0.770$ et $0.279$, de la même
façon pour les deux tailles à $1.7\cdot10^{-5}$ près ; en $t=4$, elles diffèrent ($0.026$ contre $0.035$).

<a id="s-9-6"></a>

## 9.6 Désordre : statistiques de niveaux et intrication

*Notes de cours :* [section 9.6](../../lecture_notes/fr/09_spin_chains.md#s-9-6) ; [théorème 9.16](../../lecture_notes/fr/09_spin_chains.md#r-9-16) (formule de Peschel).

La **chaîne de Heisenberg en champ aléatoire** a des champs indépendants $h_i$ uniformes sur $[-W,W]$. Elle conserve
$S^z$ (proposition 9.1), donc les statistiques de niveaux se calculent par secteur $N$ (théorème 5.15 et
section 5.4 : mélanger les secteurs simule des niveaux indépendants). Si deux champs diffèrent,
$THT^*-H=\sum_i(h_{i-1}-h_i)s^z_i\ne0$, donc les translations ne sont plus une symétrie. Nous calculons $\langle r\rangle$
dans la moitié centrale du secteur $N=L/2$, avec une erreur standard sur les réalisations indépendantes (chacune
donne un rapport moyen), et comparons à deux références : $2\log2-1\approx0.3863$, le rapport moyen exact
d'espacements exponentiels indépendants (proposition 5.14(a)), et le rapport moyen de matrices du GOE (au sens de
la section 5.3 ; le code tire $G+G^T$, le double d'une matrice du GOE, ce qui ne change aucun rapport) de même
taille $252$ que le secteur pour $L=10$, tirées ici. Ces comparaisons sont des observations numériques, pas des
théorèmes : nous ne démontrons rien sur le spectre de la chaîne aléatoire.

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

*Observation numérique :* pour $L=10$, le rapport moyen passe de $0.530$ ($W=0.5$), proche de la valeur GOE tirée
$0.533$, à $0.384$ ($W=16$), proche de $2\log2-1$, en passant par $0.484$ ($W=2$) et $0.410$ ($W=4$), avec des
erreurs standard de $0.003$–$0.004$. Pour $L=8$, la valeur en $W=0.5$ ($0.507$) est inférieure à celle en $W=1$
($0.525$). En $W=0$, l'anneau commute avec $T$ (section 9.2) et avec le spin total (section 9.3), donc son secteur
$N=L/2$ se découpe en blocs ; les champs ne couplent ces blocs que par $\sum_ih_is^z_i$, de norme d'opérateur au plus
$\frac{WL}2$.

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Les physiciens lisent la
> variation de $\langle r\rangle$ comme un crossover entre des niveaux qui se repoussent comme ceux du GOE (désordre
> faible) et des niveaux qui se comportent comme des niveaux indépendants, attribué à la **localisation à N corps**
> des états propres (désordre fort). Le crossover dérive avec $L$ ; sa position quand $L\to\infty$, et l'existence
> même d'une phase localisée, sont des questions ouvertes dans la littérature physique (D. Luitz, N. Laflorencie,
> F. Alet, *Phys. Rev. B* 91 (2015) 081103 ; J. Šuntajs, J. Bonča, T. Prosen, L. Vidmar, *Phys. Rev. E* 102 (2020)
> 062144). Le creux en $L=8$, $W=0.5$ est lu comme un vestige de la structure en blocs supplémentaire en $W=0$
> (comparer avec la section 5.4).

**Intrication des fermions libres.** Pour les états propres de la chaîne XX ouverte, l'entropie d'intrication des
$n$ premiers sites est déterminée par une matrice $n\times n$.

`entanglement_entropy` calcule l'état réduit secteur par secteur (il commute avec l'aimantation du bloc,
proposition 7.8) et son entropie. L'état fondamental du secteur $N$ est $\Phi$ avec les $N$ modes les plus bas,
simple parce que les $\varepsilon_k$ sont distincts :

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

*Observation numérique :* les entropies $0.6931,0.5199,0.7456,0.6283,0.7719,0.6560$ pour $n=1,\dots,6$ oscillent
avec la parité de $n$. Le théorème 9.16 est dû à I. Peschel (*J. Phys. A* 36 (2003) L205).

<a id="s-9-7"></a>

## 9.7 Points de comparaison

QuSpin (Python) construit des hamiltoniens de chaînes de spins dans des secteurs d'aimantation, d'impulsion, de
parité et d'inversion de spin, et utilise les solveurs creux de valeurs propres et les exponentielles de matrices
de SciPy ; des codes dédiés de diagonalisation exacte (XDiag, SPINPACK, HΦ) sont optimisés pour les grands
problèmes creux, et, avec les symétries du réseau et de grosses machines, de tels codes ont atteint environ
40–50 spins ; ITensor, TeNPy et TensorKit.jl (avec des arbres de fusion $SU(2)$) traitent les longues chaînes par
états produits de matrices. Ce qu'apporte la bibliothèque est modeste et spécifique : la chaîne en tant que
`TypeIAlgebra` $\bigoplus_NM_{\binom LN}$, de sorte que les états, canaux, entropies et dynamiques des chapitres 3–8
s'appliquent secteur par secteur ; la tour d'inclusions avec son indice ; des blocs $SU(2)$ issus des vecteurs de
plus haut poids, avec des corrélations par Wigner–Eckart ; des lots de réalisations du désordre en un seul appel ;
et des états réduits qui ne quittent jamais les secteurs. Elle n'a pas de secteurs de réflexion ni d'inversion de
spin, n'a de secteurs d'impulsion que pour les chaînes XXZ uniformes, des méthodes denses jusqu'à $L\approx16$, des
méthodes de Krylov jusqu'à $L\approx24$ (environ $26$ sur un gros GPU), et pas de réseaux de tenseurs.

<a id="x-9"></a>

## Exercices

[Exercices 1–3](../../lecture_notes/fr/09_spin_chains.md#x-9) des notes de cours : les dimensions $d_{2,m}=2,1,2,1$ des secteurs d'impulsion pour
$L=4$, $N=2$ ; les règles de somme $E_0=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle$ et $\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle=0$ pour
l'état fondamental de l'anneau de Heisenberg ; et `krylov.ground_state` lancé à partir du vecteur dont tous les coefficients
valent $1$ dans le secteur $N=L/2$, qui renvoie $E=JL/4$, le sommet du spectre. Les vérifier numériquement.

<details>
<summary>Solutions</summary>

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

**Suite :** [Tutoriel 10. Calcul quantique](10_quantum_computing.md).
