[English](../09_spin_chains.md) · **Español** · [Français](../fr/09_spin_chains.md) · [Deutsch](../de/09_spin_chains.md) · [中文](../zh/09_spin_chains.md) · [日本語](../ja/09_spin_chains.md) · [Русский](../ru/09_spin_chains.md)

# Tutorial 9. Cadenas de espines

Este tutorial acompaña al [capítulo 9](../../lecture_notes/es/09_spin_chains.md) de los apuntes de clase, que construye los sectores de la
cadena XXZ —magnetización ([proposición 9.1](../../lecture_notes/es/09_spin_chains.md#r-9-1)), momento ([proposición 9.5](../../lecture_notes/es/09_spin_chains.md#r-9-5)) y espín total— y demuestra lo
que dicen sobre los estados fundamentales ([teorema 9.7](../../lecture_notes/es/09_spin_chains.md#r-9-7), [teorema 9.9](../../lecture_notes/es/09_spin_chains.md#r-9-9)). Construimos las cadenas con
`SpinChain` y `SU2Chain`, contrastamos la biblioteca con los fermiones libres y el ansatz de Bethe, usamos los
métodos de Lanczos y de Krylov de `krylov` para sectores demasiado grandes para matrices densas, y calculamos la
estadística de niveles y el entrelazamiento.

En todo el capítulo, los sitios son $`0,\dots,L-1`$, la base de $`(\mathbb C^2)^{\otimes L}`$ está formada por las
configuraciones $`|b\rangle=|b_0\dots b_{L-1}\rangle`$ con $`b_i=1`$ para el espín hacia arriba (la convención de la
biblioteca: el bit $`i`$ de una etiqueta es el sitio $`i`$), $`s^a_i`$ ($`a=x,y,z`$) son las matrices de espín
$`\frac12\sigma_a`$ en el sitio $`i`$, $`s^\pm_i=s^x_i\pm is^y_i`$, $`n_i=s^z_i+\frac12`$, $`S^a=\sum_is^a_i`$, y
$`N=S^z+\frac L2`$ es el número de espines hacia arriba. Una cadena es **abierta** (enlaces $`(i,i+1)`$, $`i<L-1`$) o
**periódica** (también el enlace $`(L-1,0)`$; entonces $`L\ge3`$ y los índices se toman mód $`L`$). $`J>0`$
(antiferromagnética), $`\hbar=1`$, y $`\log`$ es el logaritmo natural.

<a id="s-9-1"></a>

## 9.1 La cadena XXZ y sus sectores $`U(1)`$

*Apuntes:* [sección 9.1](../../lecture_notes/es/09_spin_chains.md#s-9-1); [proposición 9.1](../../lecture_notes/es/09_spin_chains.md#r-9-1) ($`H`$ conserva $`S^z`$), [proposición 9.2](../../lecture_notes/es/09_spin_chains.md#r-9-2), [lema 9.3](../../lecture_notes/es/09_spin_chains.md#r-9-3) y [corolario 9.4](../../lecture_notes/es/09_spin_chains.md#r-9-4) (la torre de cadenas y su índice).

`SpinChain(L).xxz` construye los bloques directamente a partir de las configuraciones de cada sector
(`chain.labels[N]`). Los comparamos con la matriz completa $`2^L\times2^L`$ construida mediante productos de
Kronecker:

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

La torre $`\mathcal A_0\subset\mathcal A_1\subset\cdots`$ de las cadenas: `inclusions.chain_tower` construye las inclusiones,
`inclusions.bratteli_diagram` su diagrama de Bratteli, e `Inclusion.index` el índice $`\|\Lambda\|^2`$. El código comprueba
la regla de Pascal, el [lema 9.3](../../lecture_notes/es/09_spin_chains.md#r-9-3) y $`[\mathcal A_{L+1}:\mathcal A_L]=4\cos^2\frac\pi{2L+4}`$ ([corolario 9.4](../../lecture_notes/es/09_spin_chains.md#r-9-4)).

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

## 9.2 Traslaciones, momento y la regla de signos de Marshall

*Apuntes:* [sección 9.2](../../lecture_notes/es/09_spin_chains.md#s-9-2); [proposición 9.5](../../lecture_notes/es/09_spin_chains.md#r-9-5) (sectores de momento), [teorema 9.7](../../lecture_notes/es/09_spin_chains.md#r-9-7) (regla de signos de Marshall), [corolario 9.8](../../lecture_notes/es/09_spin_chains.md#r-9-8) (momento del estado fundamental).

`xxz_momentum` construye los bloques en la base normalizada $`|a,k\rangle`$, con un representante $`a`$ por órbita.
Contamos las órbitas de forma independiente y comparamos los espectros con el sector denso:

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

El sector $`N=6`$, de $`924`$ estados, se divide en bloques de tamaños $`80,75,78,76,78,75,80,75,78,76,78,75`$.

Por la regla de signos de Marshall ([teorema 9.7](../../lecture_notes/es/09_spin_chains.md#r-9-7)), para cadenas abiertas y periódicas con campos aleatorios,
$`(-1)^{N_B(b)}\psi(b)`$ tiene un solo signo en cada sector, y por el [corolario 9.8](../../lecture_notes/es/09_spin_chains.md#r-9-8) el nivel más bajo del sector $`N`$ del
anillo está en el bloque de momento $`k=0`$ para $`N`$ par y $`k=\pi`$ para $`N`$ impar:

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

## 9.3 $`SU(2)`$ en $`\Delta=1`$

*Apuntes:* [sección 9.3](../../lecture_notes/es/09_spin_chains.md#s-9-3); [teorema 9.9](../../lecture_notes/es/09_spin_chains.md#r-9-9) (Lieb–Mattis). Para $`\Delta=1`$ y $`h=0`$, $`H\cong\bigoplus_S1_{2S+1}\otimes H_S`$, y `SU2Chain` construye los bloques $`H_S`$ sobre los vectores de peso máximo.

Comparamos los tamaños de los bloques para doce espines, comprobamos el teorema 9.9 frente a los sectores $`U(1)`$,
calculamos la brecha y comprobamos los teoremas 9.7 y 9.9 en un grafo bipartito que no es una cadena ($`n_A=5`$,
$`n_B=3`$, $`S_*=1`$, acoplamientos aleatorios):

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

Los bloques de $`SU(2)`$ tienen tamaños $`m_S=132,297,275,154,54,11,1`$ frente a $`924`$ del mayor bloque de $`U(1)`$, y
el coste $`\sum k^3`$ de la diagonalización baja de $`2.0\cdot10^9`$ a $`5.3\cdot10^7`$. *Observación numérica:* $`L`$ veces
la brecha crece lentamente, $`4.108,4.181,4.232,4.270`$ para $`L=6,\dots,12`$; no demostramos nada sobre su
comportamiento cuando $`L\to\infty`$.

**Correlaciones.** Los espines de los sitios forman operadores vectoriales, y sus componentes esféricas,
operadores tensoriales de rango uno (lema B.15(a),(b)). Por el **teorema de Wigner–Eckart** (teorema B.16), sus
elementos de matriz entre dos multipletes son un coeficiente de Clebsch–Gordan por un **elemento de matriz
reducido** que no depende de los números cuánticos magnéticos, y
$`\mathbf s_i\cdot\mathbf s_j=\sum_q(-1)^qs^{(i)}_qs^{(j)}_{-q}`$ (lema B.15(a)). `spin_correlation` lo usa para construir
$`\mathbf s_i\cdot\mathbf s_j`$ sobre cada $`H_S`$ a partir únicamente de los vectores de peso máximo; el método `'direct'`
comprime en cambio el operador de enlace a los vectores de peso máximo, lo cual es correcto por la
proposición 4.15(b), ya que $`\mathbf s_i\cdot\mathbf s_j\in\{S_a\}'`$. Que la implementación del primer método es correcta no
se demuestra aquí; el código comprueba que ambos coinciden entre sí y con la definición, evaluados en el estado
fundamental escrito en el sector $`U(1)`$:

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

Las correlaciones $`\langle\mathbf s_0\cdot\mathbf s_r\rangle=0.75,-0.4489,0.1878,-0.1658,0.1209,-0.1224,0.1070,\dots`$
alternan de signo, como imponen los teoremas 9.7 y 9.9 para todo anillo par ([apuntes de clase, sección 9.3](../../lecture_notes/es/09_spin_chains.md#s-9-3)).

<a id="s-9-4"></a>

## 9.4 Comprobaciones exactas: fermiones libres y el ansatz de Bethe

*Apuntes:* [sección 9.4](../../lecture_notes/es/09_spin_chains.md#s-9-4); por el [teorema 9.10](../../lecture_notes/es/09_spin_chains.md#r-9-10), las energías del sector $`N`$ de la cadena XX abierta son las sumas de $`N`$ valores distintos $`\varepsilon_k=J\cos\frac{\pi k}{L+1}`$; [proposición 9.12](../../lecture_notes/es/09_spin_chains.md#r-9-12) (dos magnones, ecuaciones de Bethe).

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

Por la [proposición 9.12](../../lecture_notes/es/09_spin_chains.md#r-9-12), cada raíz $`q`$ de $`F(q)=(-1)^ne^{iqL}(c-\Delta e^{-iq})+c-\Delta e^{iq}`$, con
$`K=2\pi n/L`$ y $`c=\cos\frac K2`$, da un vector propio $`\Psi`$ del sector $`N=2`$ del anillo
([apuntes de clase, sección 9.4](../../lecture_notes/es/09_spin_chains.md#s-9-4)). El código busca raíces $`q\in(0,\pi)`$ (estados de dispersión) y sobre las rectas
$`q=i\kappa`$, $`\pi+i\kappa`$, $`\kappa>0`$ (estados ligados), construye $`\Psi`$ para cada raíz para $`L=9`$, $`\Delta=0.6`$, y compara
con el sector $`N=2`$ de $`36`$ estados:

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

*Observación numérica:* para $`L=9`$ y $`\Delta=0.6`$, la búsqueda encuentra $`36`$ raíces, $`34`$ estados de dispersión y
$`2`$ estados ligados, cuyos vectores son vectores propios linealmente independientes (salvo redondeo): los $`36`$
estados propios provienen de las ecuaciones de Bethe. La completitud no está demostrada para otros $`L`$, y para $`L`$
par y $`\Delta\ne0`$ falla ([apuntes de clase, sección 9.4](../../lecture_notes/es/09_spin_chains.md#s-9-4)): el vector $`\Phi=\sum_x(-1)^x|x,x+1\rangle`$ es ortogonal a todos
los vectores de Bethe, y las últimas líneas del código comprueban $`H\Phi=(E_F-J\Delta)\Phi`$ para $`L=8`$.

<a id="s-9-5"></a>

## 9.5 Lanczos y métodos de Krylov

*Apuntes:* [sección 9.5](../../lecture_notes/es/09_spin_chains.md#s-9-5); [proposición 9.13](../../lecture_notes/es/09_spin_chains.md#r-9-13) (Lanczos: $`T_m=Q_m^*HQ_m`$ es tridiagonal), [teorema 9.14](../../lecture_notes/es/09_spin_chains.md#r-9-14) (los valores de Ritz, los valores propios de $`T_m`$, cumplen $`\theta^{(m)}_j\ge\lambda_j`$ y decrecen con $`m`$), [proposición 9.15](../../lecture_notes/es/09_spin_chains.md#r-9-15) (evolución temporal de Krylov), [lema 9.17](../../lecture_notes/es/09_spin_chains.md#r-9-17), [teorema 9.18](../../lecture_notes/es/09_spin_chains.md#r-9-18) y [corolario 9.19](../../lecture_notes/es/09_spin_chains.md#r-9-19) (cotas de error y la regla de paso).

El sector $`N=L/2`$ tiene dimensión $`12870`$ para $`L=16`$ y $`184756`$ para $`L=20`$, pero $`H`$ es dispersa ($`O(L)`$ entradas
no nulas por fila, `xxz_sparse`), y los métodos de Krylov solo necesitan productos $`Hv`$. Sea $`H`$ hermítica de
tamaño $`n`$ y $`v\ne0`$. El **espacio de Krylov** es $`\mathcal K_m=\mathrm{span}\{v,Hv,\dots,H^{m-1}v\}`$.

El teorema 9.14 no dice nada sobre la rapidez con que $`\theta^{(m)}_1`$ se acerca a $`\lambda_1`$; solo la observamos
más abajo. En coma flotante, los $`q_j`$ pierden la ortogonalidad; `krylov.lanczos` reortogonaliza por completo, y
el código comprueba que los $`q_j`$ calculados son ortonormales y que $`Q_m^*HQ_m=T_m`$ salvo redondeo. Para la cadena
XX abierta, la energía exacta del estado fundamental del sector $`N`$ es la suma de los $`N`$ menores $`\varepsilon_k`$
(teorema 9.10(c)):

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

*Observación numérica:* el error del menor valor de Ritz cae monótonamente (como debe, por el teorema 9.14),
$`0.42`$, $`3.1\cdot10^{-3}`$ y $`3.1\cdot10^{-6}`$ tras $`5`$, $`10`$ y $`20`$ pasos, y alcanza el nivel de redondeo tras $`40`$, en
dimensión $`924`$. `krylov.ground_state` reinicia Lanczos (desde un vector aleatorio gaussiano, salvo que se dé
`v0`) hasta que el residuo $`\|H\psi-E\psi\|`$ de su par de Ritz normalizado $`(E,\psi)`$ queda por debajo de
$`10^{-10}\max(1,|E|)`$ (como mucho $`20`$ reinicios, con una advertencia en caso contrario), y devuelve el residuo.
Un residuo pequeño certifica un valor propio cercano, pero no que sea el menor ([apuntes de clase, sección 9.5](../../lecture_notes/es/09_spin_chains.md#s-9-5)), así que
comparamos con la diagonalización densa, que es
nuestra única referencia (calculada aquí, así que los valores son numéricos): el sector $`N=L/2`$ hasta $`L=12`$ y,
para $`L=16`$, el bloque de momento $`k=0`$ de dimensión $`810`$, que contiene el estado fundamental del sector $`N=8`$ por
el corolario 9.8. Para $`L=4`$ el valor $`-2J`$ es exacto ([apuntes de clase, sección 9.5](../../lecture_notes/es/09_spin_chains.md#s-9-5)).

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

Para $`L=16`$ (dimensión $`12870`$), Lanczos tarda bastante menos de un segundo en un hilo de CPU. *Observación
numérica:* la energía por sitio, $`-0.5,-0.467129,-0.456387,-0.451545,-0.448949`$ para $`L=4,\dots,12`$, crece con $`L`$;
no demostramos nada sobre su límite.

**Evolución temporal.** Para $`e^{-i\tau H}v`$ se proyecta sobre $`\mathcal K_m`$: $`w_\tau=\|v\|\,Q_me^{-i\tau T_m}e_1`$.

Por la [proposición 9.15](../../lecture_notes/es/09_spin_chains.md#r-9-15), $`\|w_\tau\|=\|v\|`$ y $`\langle w_\tau,Hw_\tau\rangle=\langle v,Hv\rangle`$ para todo $`\tau`$, y
$`w_\tau=e^{-i\tau H}v`$ si $`\mathcal K_m`$ es $`H`$-invariante.

Así, cada paso de Krylov conserva la norma y la energía, sea preciso o no. Su precisión es el objeto del
[teorema 9.18](../../lecture_notes/es/09_spin_chains.md#r-9-18): un paso de longitud $`\tau`$ yerra a lo sumo en
$`(4\beta_m/w_T)S_m(|\tau|w_T/2)`$, calculado solo a partir de la ejecución de Lanczos. `krylov.evolve` da cada paso, en
un espacio de Krylov con $`m=30`$, tan largo como esta cota permite para un error de a lo sumo $`\mathrm{tol}=10^{-8}`$
(el valor por defecto) en todo el intervalo de tiempo ([corolario 9.19](../../lecture_notes/es/09_spin_chains.md#r-9-19)).
Lo comparamos con $`e^{-itH}`$ calculada a partir de una descomposición espectral densa (teorema 1.5), y comprobamos
las propias cotas más abajo. Realizamos un
*quench* desde el estado de Néel $`|{\uparrow\downarrow\uparrow\downarrow}\cdots\rangle`$ en el anillo de Heisenberg y
seguimos la magnetización alternada $`m_s=\frac2L\sum_i(-1)^i\langle s^z_i\rangle`$, comparando con el propagador denso
para $`L=12`$:

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

*Observación numérica:* $`m_s`$ decae desde $`1`$, pasando por $`0.770`$ y $`0.279`$, hasta $`-0.183`$ en $`t=2`$, igual para
ambos tamaños salvo $`1.7\cdot10^{-5}`$; en $`t=4`$ difieren ($`0.026`$ frente a $`0.035`$).

**Cotas de error.** Por el [lema 9.17](../../lecture_notes/es/09_spin_chains.md#r-9-17), los coeficientes de Chebyshev de
$`e^{-i\rho y}`$ son $`a_k=(-i)^kJ_k(\rho)`$ con $`|a_k|\le(|\rho|/2)^k/k!`$. El [teorema 9.18](../../lecture_notes/es/09_spin_chains.md#r-9-18)
acota el error de un paso de Krylov por $`(4\beta_m/w_T)S_m(|\tau|w_T/2)`$ (c), el valor que devuelve
`krylov.expm_multiply_step`, y por $`4S_m(\rho/2)`$ con $`\rho=|\tau|w`$ (d), donde $`S_m(x)=\sum_{j\ge m}x^j/j!`$ y $`w`$ es la
semianchura del espectro. Comprobamos la cota de los coeficientes frente a las funciones de Bessel de SciPy, y ambas
cotas de error para un paso desde el estado de Néel del anillo $`L=12`$ con $`m=20`$, frente al propagador denso del
bloque anterior (`w`, `U`). El [corolario 9.19](../../lecture_notes/es/09_spin_chains.md#r-9-19) suma las cotas de pasos
sucesivos; y la regla de paso de `krylov.evolve` mantiene el error en el instante $`t`$ por debajo de
$`\mathrm{tol}\cdot t/t_{\mathrm{end}}`$, con los mismos pasos para $`H`$ y $`H+10^4`$ (teorema 9.18(a); el operador desplazado solo
cambia la diagonal de la forma dispersa):

```python
from scipy.special import jv, gammainc
S = lambda m, x: math.exp(x) * gammainc(m, x)                        # S_m(x) = sum_{j >= m} x^j / j!
for rho in (1.0, 10.0, 40.0):                                        # Lemma 9.17: |J_k(rho)| <= (rho/2)^k / k!
    assert all(abs(jv(k, rho)) <= (rho / 2) ** k / math.factorial(k) * (1 + 1e-12) for k in range(80))

c12 = SpinChain(12, 'periodic')
H12 = c12.xxz_sparse(J=1.0, Delta=1.0, sector=6)                    # w, U: its dense eigendecomposition (above)
half = ((w[-1] - w[0]) / 2).item()                                   # half-width of the spectrum
v = c12.vector_in_sector('10' * 6)[None]                             # the Neel state
exact = lambda t: U @ (torch.exp(-1j * t * w) * (U.mH @ v[0]))
m = 20
for tau in (0.5, 1.0, 2.0, 4.0, 8.0):
    out, bound = krylov.expm_multiply_step(H12, v, tau, m)            # bound: Theorem 9.18(c)
    err = (out[0] - exact(tau)).norm().item()
    apriori = 4 * S(m, tau * half / 2)                                # Theorem 9.18(d)
    print(f"tau w = {tau * half:4.1f}: error {err:.1e}, bound (c) {bound.item():.1e}, bound (d) {apriori:.1e}")
    assert err <= bound.item() + 1e-13 and bound.item() <= apriori

psi, total = v, 0.0                                                  # Corollary 9.19, five steps of length 1
for k in range(1, 6):
    psi, b = krylov.expm_multiply_step(H12, psi, 1.0, m)
    total += b.item()
    err = (psi[0] - exact(float(k))).norm().item()
    assert err <= total + 1e-12
print(f"after 5 steps: error {err:.1e}, sum of the bounds {total:.1e}")

times = [0.0, 1.0, 3.0, 10.0]
shifted = krylov.SparseSectorHamiltonian(H12.diag + 1e4, H12.offdiag, H12.dtype)    # H + 10^4
for s, Hn in ((0.0, H12), (1e4, shifted)):
    out = krylov.evolve(Hn, v, times, tol=1e-8)                      # the step rule of Corollary 9.19
    errs = [(p[0] - math.e ** (-1j * t * s) * exact(t)).norm().item() for t, p in zip(times, out)]
    print(f"shift {s:.0e}:", [f"{e:.1e}" for e in errs])
    assert all(e <= 1e-8 * t / times[-1] + 1e-10 for e, t in zip(errs, times))
```

*Observación numérica:* el error de un paso crece desde el nivel de redondeo ($`3\cdot10^{-15}`$ en $`\tau w=2.1`$) hasta
$`1.3\cdot10^{-8}`$ en $`\tau w=8.4`$ y $`0.5`$ en $`\tau w=34`$; la cota (c) lo supera en un factor de unos $`100`$ en el rango
útil, y (d) en otro factor $`4`$. Cinco pasos de longitud $`1`$ yerran en $`10^{-13}`$, por debajo de la suma
$`6\cdot10^{-12}`$ de sus cotas. `krylov.evolve` cumple aquí su tolerancia con amplio margen (errores por debajo de
$`10^{-13}`$); con el desplazamiento $`10^4`$ los errores son $`10^{-11}`$, los errores de redondeo de productos
matriz–vector de norma $`10^4`$, que ninguna cota del teorema 9.18 cubre.

<a id="s-9-6"></a>

## 9.6 Desorden: estadística de niveles y entrelazamiento

*Apuntes:* [sección 9.6](../../lecture_notes/es/09_spin_chains.md#s-9-6); [teorema 9.16](../../lecture_notes/es/09_spin_chains.md#r-9-16) (fórmula de Peschel).

La **cadena de Heisenberg con campo aleatorio** tiene campos independientes $`h_i`$ uniformes en $`[-W,W]`$. Conserva
$`S^z`$ (proposición 9.1), así que la estadística de niveles se calcula por sector $`N`$ (teorema 5.15 y sección 5.4:
mezclar sectores simula niveles independientes). Si dos campos difieren, $`THT^*-H=\sum_i(h_{i-1}-h_i)s^z_i\ne0`$, así
que las traslaciones dejan de ser una simetría. Calculamos $`\langle r\rangle`$ en la mitad central del sector $`N=L/2`$,
con un error estándar sobre las realizaciones independientes (cada una da un cociente medio), y comparamos con
dos referencias: $`2\log2-1\approx0.3863`$, el cociente medio exacto de espaciamientos exponenciales independientes
(proposición 5.14(a)), y el cociente medio de matrices del GOE (definidas como en la sección 5.3; el código
muestrea $`G+G^T`$, el doble de una matriz del GOE, lo que no cambia ningún cociente) del mismo tamaño $`252`$ que el
sector para $`L=10`$, muestreadas aquí. Las comparaciones son observaciones numéricas, no teoremas: no demostramos
nada sobre el espectro de la cadena aleatoria.

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

*Observación numérica:* para $`L=10`$, el cociente medio va de $`0.530`$ ($`W=0.5`$), cercano al valor muestreado del
GOE, $`0.533`$, pasando por $`0.484`$ ($`W=2`$) y $`0.410`$ ($`W=4`$), hasta $`0.384`$ ($`W=16`$), cercano a $`2\log2-1`$, con errores
estándar de $`0.003`$–$`0.004`$. Para $`L=8`$, el valor en $`W=0.5`$ ($`0.507`$) está por debajo del de $`W=1`$ ($`0.525`$). En
$`W=0`$ el anillo conmuta con $`T`$ (sección 9.2) y con el espín total (sección 9.3), así que su sector $`N=L/2`$ se
divide en bloques; los campos acoplan estos bloques solo a través de $`\sum_ih_is^z_i`$, de norma de operador a lo
sumo $`\frac{WL}2`$.

> **Significado físico** (interpretación física, no un teorema de este curso). Los físicos interpretan el cambio
> de $`\langle r\rangle`$ como una transición desde niveles que se repelen como los del GOE (desorden débil) hacia
> niveles que se comportan como independientes, atribuida a la **localización de muchos cuerpos** de los estados
> propios (desorden fuerte). La transición se desplaza con $`L`$; dónde se sitúa cuando $`L\to\infty`$, y si existe
> siquiera una fase localizada, son cuestiones abiertas en la literatura física (D. Luitz, N. Laflorencie,
> F. Alet, *Phys. Rev. B* 91 (2015) 081103; J. Šuntajs, J. Bonča, T. Prosen, L. Vidmar, *Phys. Rev. E* 102 (2020)
> 062144). La caída en $`L=8`$, $`W=0.5`$ se interpreta como un resto de la estructura de bloques adicional en $`W=0`$
> (compárese con la sección 5.4).

**Entrelazamiento de fermiones libres.** Para los estados propios de la cadena XX abierta, la entropía de
entrelazamiento de los primeros $`n`$ sitios está determinada por una matriz $`n\times n`$.

`entanglement_entropy` calcula el estado reducido sector a sector (conmuta con la magnetización del bloque,
proposición 7.8) y su entropía. El estado fundamental del sector $`N`$ es $`\Phi`$ con los $`N`$ modos más bajos, simple
porque los $`\varepsilon_k`$ son distintos:

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

*Observación numérica:* las entropías $`0.6931,0.5199,0.7456,0.6283,0.7719,0.6560`$ para $`n=1,\dots,6`$ oscilan con la
paridad de $`n`$. El teorema 9.16 se debe a I. Peschel (*J. Phys. A* 36 (2003) L205).

<a id="s-9-7"></a>

## 9.7 Con qué comparar

QuSpin (Python) construye hamiltonianos de cadenas de espines en sectores de magnetización, momento, paridad e
inversión de espín, y usa los resolvedores dispersos de valores propios y las exponenciales de matrices de SciPy;
los códigos dedicados de diagonalización exacta (XDiag, SPINPACK, HΦ) están optimizados para grandes problemas
dispersos, y con simetrías de red y máquinas grandes tales códigos han alcanzado unos 40–50 espines; ITensor,
TeNPy y TensorKit.jl (con árboles de fusión de $`SU(2)`$) tratan cadenas largas con estados de producto de
matrices. Lo que aporta la biblioteca es modesto y específico: la cadena como `TypeIAlgebra`
$`\bigoplus_NM_{\binom LN}`$, de modo que los estados, canales, entropías y dinámicas de los capítulos 3–8 se aplican
sector a sector; la torre de inclusiones con su índice; bloques de $`SU(2)`$ a partir de vectores de peso máximo
con correlaciones de Wigner–Eckart; lotes de realizaciones del desorden en una sola llamada; y estados reducidos
que nunca salen de los sectores. No tiene sectores de reflexión ni de inversión de espín, tiene sectores de
momento solo para cadenas XXZ uniformes, métodos densos hasta $`L\approx16`$, métodos de Krylov hasta $`L\approx24`$
(unos $`26`$ en una GPU grande), y no tiene redes tensoriales.

<a id="x-9"></a>

## Ejercicios

[Ejercicios 1–3](../../lecture_notes/es/09_spin_chains.md#x-9) de los apuntes de clase: las dimensiones $`d_{2,m}=2,1,2,1`$ de los sectores de momento para
$`L=4`$, $`N=2`$; las reglas de suma $`E_0=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle`$ y $`\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle=0`$ para el
estado fundamental del anillo de Heisenberg; y `krylov.ground_state` iniciado desde el vector con todos los coeficientes
iguales a $`1`$ en el sector $`N=L/2`$, que devuelve $`E=JL/4`$, el extremo superior del espectro. Compruébelos numéricamente.

<details>
<summary>Soluciones</summary>

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

**Siguiente:** [Tutorial 10. Computación cuántica](10_quantum_computing.md).
