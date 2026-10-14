**English** · [Español](es/09_spin_chains.md) · [Français](fr/09_spin_chains.md) · [Deutsch](de/09_spin_chains.md) · [中文](zh/09_spin_chains.md) · [日本語](ja/09_spin_chains.md) · [Русский](ru/09_spin_chains.md)

# Tutorial 9. Spin chains

This tutorial accompanies [Chapter 9](../lecture_notes/09_spin_chains.md) of the lecture notes, which builds the sectors of the XXZ
chain — magnetisation ([Proposition 9.1](../lecture_notes/09_spin_chains.md#r-9-1)), momentum ([Proposition 9.5](../lecture_notes/09_spin_chains.md#r-9-5)) and total spin — and proves what
they say about ground states ([Theorem 9.7](../lecture_notes/09_spin_chains.md#r-9-7), [Theorem 9.9](../lecture_notes/09_spin_chains.md#r-9-9)). We build the chains with `SpinChain` and
`SU2Chain`, check the library against free fermions and the Bethe ansatz, use the Lanczos and Krylov methods of
`krylov` for sectors too large for dense matrices, and compute level statistics and entanglement.

Throughout, the sites are $0,\dots,L-1$, the basis of $(\mathbb C^2)^{\otimes L}$ consists of the
configurations $|b\rangle=|b_0\dots b_{L-1}\rangle$ with $b_i=1$ for spin up (the library's convention: bit
$i$ of a label is site $i$), $s^a_i$ ($a=x,y,z$) are the spin matrices $\frac12\sigma_a$ on site $i$,
$s^\pm_i=s^x_i\pm is^y_i$, $n_i=s^z_i+\frac12$, $S^a=\sum_is^a_i$, and $N=S^z+\frac L2$ is the number of
up spins. A chain is **open** (bonds $(i,i+1)$, $i<L-1$) or **periodic** (also the bond $(L-1,0)$; then
$L\ge3$, indices mod $L$). $J>0$ (antiferromagnetic), $\hbar=1$, $\log$ is natural.

<a id="s-9-1"></a>

## 9.1 The XXZ chain and its $U(1)$ sectors

*Lecture notes:* [Section 9.1](../lecture_notes/09_spin_chains.md#s-9-1); [Proposition 9.1](../lecture_notes/09_spin_chains.md#r-9-1) ($H$ conserves $S^z$), [Proposition 9.2](../lecture_notes/09_spin_chains.md#r-9-2), [Lemma 9.3](../lecture_notes/09_spin_chains.md#r-9-3) and [Corollary 9.4](../lecture_notes/09_spin_chains.md#r-9-4) (the tower of chains and its index).

`SpinChain(L).xxz` builds the blocks directly from the configurations of each sector (`chain.labels[N]`).
We compare them with the full $2^L\times2^L$ matrix built from Kronecker products:

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

The tower $\mathcal A_0\subset\mathcal A_1\subset\cdots$ of the chains: `inclusions.chain_tower` builds the inclusions,
`inclusions.bratteli_diagram` their Bratteli diagram, and `Inclusion.index` the index $\|\Lambda\|^2$. The code checks
Pascal's rule, [Lemma 9.3](../lecture_notes/09_spin_chains.md#r-9-3) and $[\mathcal A_{L+1}:\mathcal A_L]=4\cos^2\frac\pi{2L+4}$ ([Corollary 9.4](../lecture_notes/09_spin_chains.md#r-9-4)).

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

## 9.2 Translations, momentum and the Marshall sign rule

*Lecture notes:* [Section 9.2](../lecture_notes/09_spin_chains.md#s-9-2); [Proposition 9.5](../lecture_notes/09_spin_chains.md#r-9-5) (momentum sectors), [Theorem 9.7](../lecture_notes/09_spin_chains.md#r-9-7) (Marshall sign rule), [Corollary 9.8](../lecture_notes/09_spin_chains.md#r-9-8) (momentum of the ground state).

`xxz_momentum` builds the blocks in the normalised basis $|a,k\rangle$, one representative $a$ per orbit. We
count orbits independently and compare spectra with the dense sector:

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

The sector $N=6$ of $924$ states splits into blocks of sizes $80,75,78,76,78,75,80,75,78,76,78,75$.

By the Marshall sign rule ([Theorem 9.7](../lecture_notes/09_spin_chains.md#r-9-7)), for open and periodic chains with random fields
$(-1)^{N_B(b)}\psi(b)$ has one sign in every sector, and by [Corollary 9.8](../lecture_notes/09_spin_chains.md#r-9-8) the lowest level of the sector $N$ of
the ring lies in the momentum block $k=0$ for even $N$ and $k=\pi$ for odd $N$:

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

## 9.3 $SU(2)$ at $\Delta=1$

*Lecture notes:* [Section 9.3](../lecture_notes/09_spin_chains.md#s-9-3); [Theorem 9.9](../lecture_notes/09_spin_chains.md#r-9-9) (Lieb–Mattis). For $\Delta=1$ and $h=0$, $H\cong\bigoplus_S1_{2S+1}\otimes H_S$, and `SU2Chain` builds the blocks $H_S$ on the highest-weight vectors.

We compare block sizes for twelve spins, check Theorem 9.9 against the $U(1)$ sectors, compute the gap, and
check Theorems 9.7 and 9.9 on a bipartite graph that is not a chain ($n_A=5$, $n_B=3$, $S_*=1$, random couplings):

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

The $SU(2)$ blocks have sizes $m_S=132,297,275,154,54,11,1$ against $924$ for the largest $U(1)$ block, and the
cost $\sum k^3$ of diagonalising drops from $2.0\cdot10^9$ to $5.3\cdot10^7$. *Numerical observation:* $L$ times
the gap grows slowly, $4.108,4.181,4.232,4.270$ for $L=6,\dots,12$; we prove nothing about its behaviour as
$L\to\infty$.

**Correlations.** The site spins form vector operators, and their spherical components tensor operators of rank
one (Lemma B.15(a),(b)). By the **Wigner–Eckart theorem** (Theorem B.16) their matrix elements between two
multiplets are a Clebsch–Gordan coefficient times a **reduced matrix element** that does not depend on the
magnetic quantum numbers, and $\mathbf s_i\cdot\mathbf s_j=\sum_q(-1)^qs^{(i)}_qs^{(j)}_{-q}$ (Lemma B.15(a)).
`spin_correlation` uses this to build $\mathbf s_i\cdot\mathbf s_j$ on every $H_S$ from highest-weight vectors alone; the
`'direct'` method instead compresses the bond operator to the highest-weight vectors, which is correct by
Proposition 4.15(b) since $\mathbf s_i\cdot\mathbf s_j\in\{S_a\}'$. That the implementation of the first method is correct
is not proved here; the code checks that both agree with each other and with the definition, evaluated on
the ground state written out in the $U(1)$ sector:

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

The correlations $\langle\mathbf s_0\cdot\mathbf s_r\rangle=0.75,-0.4489,0.1878,-0.1658,0.1209,-0.1224,0.1070,\dots$
alternate in sign, as Theorems 9.7 and 9.9 force for every even ring ([lecture notes, Section 9.3](../lecture_notes/09_spin_chains.md#s-9-3)).

<a id="s-9-4"></a>

## 9.4 Exact checks: free fermions and the Bethe ansatz

*Lecture notes:* [Section 9.4](../lecture_notes/09_spin_chains.md#s-9-4); by [Theorem 9.10](../lecture_notes/09_spin_chains.md#r-9-10) the energies of the sector $N$ of the open XX chain are the sums of $N$ distinct $\varepsilon_k=J\cos\frac{\pi k}{L+1}$; [Proposition 9.12](../lecture_notes/09_spin_chains.md#r-9-12) (two magnons, Bethe equations).

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

By [Proposition 9.12](../lecture_notes/09_spin_chains.md#r-9-12), every root $q$ of $F(q)=(-1)^ne^{iqL}(c-\Delta e^{-iq})+c-\Delta e^{iq}$, with
$K=2\pi n/L$ and $c=\cos\frac K2$, gives an eigenvector $\Psi$ of the sector $N=2$ of the ring
([lecture notes, Section 9.4](../lecture_notes/09_spin_chains.md#s-9-4)). The code looks for roots $q\in(0,\pi)$ (scattering states) and on the lines
$q=i\kappa$, $\pi+i\kappa$, $\kappa>0$ (bound states), builds $\Psi$ for each root for $L=9$, $\Delta=0.6$, and compares
with the sector $N=2$ of $36$ states:

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

*Numerical observation:* for $L=9$ and $\Delta=0.6$ the search finds $36$ roots, $34$ scattering states and $2$
bound states, whose vectors are linearly independent eigenvectors (up to rounding): all $36$ eigenstates come
from the Bethe equations. Completeness is not proved for other $L$, and for even $L$ and $\Delta\ne0$ it fails
([lecture notes, Section 9.4](../lecture_notes/09_spin_chains.md#s-9-4)): the vector $\Phi=\sum_x(-1)^x|x,x+1\rangle$ is orthogonal to all Bethe vectors,
and the last lines of the code check $H\Phi=(E_F-J\Delta)\Phi$ for $L=8$.

<a id="s-9-5"></a>

## 9.5 Lanczos and Krylov methods

*Lecture notes:* [Section 9.5](../lecture_notes/09_spin_chains.md#s-9-5); [Proposition 9.13](../lecture_notes/09_spin_chains.md#r-9-13) (Lanczos: $T_m=Q_m^*HQ_m$ is tridiagonal), [Theorem 9.14](../lecture_notes/09_spin_chains.md#r-9-14) (the Ritz values, the eigenvalues of $T_m$, satisfy $\theta^{(m)}_j\ge\lambda_j$ and decrease with $m$), [Proposition 9.15](../lecture_notes/09_spin_chains.md#r-9-15) (Krylov time evolution).

The sector $N=L/2$ has dimension $12870$ for $L=16$ and $184756$ for $L=20$, but $H$ is sparse ($O(L)$ nonzero
entries per row, `xxz_sparse`), and Krylov methods need only products $Hv$. Let $H$ be Hermitian of size $n$ and
$v\ne0$. The **Krylov space** is $\mathcal K_m=\mathrm{span}\{v,Hv,\dots,H^{m-1}v\}$.

Theorem 9.14 says nothing about how fast $\theta^{(m)}_1$ approaches $\lambda_1$; we only observe it below. In
floating point the $q_j$ lose orthogonality; `krylov.lanczos` re-orthogonalises fully, and the code checks that
the computed $q_j$ are orthonormal and $Q_m^*HQ_m=T_m$ to rounding. For the open XX chain the exact ground-state
energy of sector $N$ is the sum of the $N$ lowest $\varepsilon_k$ (Theorem 9.10(c)):

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

*Numerical observation:* the error of the lowest Ritz value falls monotonically (as it must, by Theorem 9.14),
$0.42$, $3.1\cdot10^{-3}$ and $3.1\cdot10^{-6}$ after $5$, $10$ and $20$ steps, and reaches rounding level after $40$,
in dimension $924$. `krylov.ground_state` restarts Lanczos (from a Gaussian random vector unless `v0` is given)
until the residual $\|H\psi-E\psi\|$ of its normalised Ritz pair $(E,\psi)$ is below $10^{-10}\max(1,|E|)$ (at
most $20$ restarts, with a warning otherwise), and returns the residual. A small residual certifies an eigenvalue
nearby, but not that it is the lowest one ([lecture notes, Section 9.5](../lecture_notes/09_spin_chains.md#s-9-5)), so we compare with dense
diagonalisation, which is our only reference (computed here, so the values are numerical): the sector $N=L/2$ up
to $L=12$, and for $L=16$ the momentum block $k=0$ of dimension $810$, which contains the ground state of the
sector $N=8$ by Corollary 9.8. For $L=4$ the value $-2J$ is exact ([lecture notes, Section 9.5](../lecture_notes/09_spin_chains.md#s-9-5)).

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

For $L=16$ (dimension $12870$) Lanczos takes well under a second on one CPU thread. *Numerical observation:*
the energy per site, $-0.5,-0.467129,-0.456387,-0.451545,-0.448949$ for $L=4,\dots,12$, increases with $L$; we
prove nothing about its limit.

**Time evolution.** For $e^{-i\tau H}v$ one projects onto $\mathcal K_m$: $w_\tau=\|v\|\,Q_me^{-i\tau T_m}e_1$.

By [Proposition 9.15](../lecture_notes/09_spin_chains.md#r-9-15), $\|w_\tau\|=\|v\|$ and $\langle w_\tau,Hw_\tau\rangle=\langle v,Hv\rangle$ for every $\tau$, and
$w_\tau=e^{-i\tau H}v$ if $\mathcal K_m$ is $H$-invariant.

So every Krylov step conserves norm and energy, accurate or not; Proposition 9.15 says nothing about
accuracy unless $\mathcal K_m$ is invariant. `krylov.evolve` subdivides the time interval into short steps, each
done in a Krylov space with $m=30$, and warns when a heuristic error indicator is large; we prove no error
bound for it and instead compare it with $e^{-itH}$ computed from a dense eigendecomposition (Theorem 1.5).
We quench the Néel
state $|{\uparrow\downarrow\uparrow\downarrow}\cdots\rangle$ on the Heisenberg ring and follow the staggered magnetisation
$m_s=\frac2L\sum_i(-1)^i\langle s^z_i\rangle$, comparing with the dense propagator for $L=12$:

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

*Numerical observation:* $m_s$ decays from $1$ through $0.770$ and $0.279$ to $-0.183$ at $t=2$, the same for both
sizes to $1.7\cdot10^{-5}$; at $t=4$ they differ ($0.026$ against $0.035$).

<a id="s-9-6"></a>

## 9.6 Disorder: level statistics and entanglement

*Lecture notes:* [Section 9.6](../lecture_notes/09_spin_chains.md#s-9-6); [Theorem 9.16](../lecture_notes/09_spin_chains.md#r-9-16) (Peschel's formula).

The **random-field Heisenberg chain** has independent fields $h_i$ uniform on $[-W,W]$. It conserves $S^z$
(Proposition 9.1), so level statistics are computed per sector $N$ (Theorem 5.15 and Section 5.4: mixing sectors
fakes independent levels). If two fields differ, $THT^*-H=\sum_i(h_{i-1}-h_i)s^z_i\ne0$, so translations are no longer
a symmetry. We compute $\langle r\rangle$ in the middle half of the sector $N=L/2$, with one standard error over the
independent realisations (each gives one mean ratio), and compare with two references: $2\log2-1\approx0.3863$,
the exact mean ratio of independent exponential spacings (Proposition 5.14(a)), and the mean ratio of GOE matrices
(as defined in Section 5.3; the code samples $G+G^T$, twice a GOE matrix, which leaves every ratio unchanged)
of the same size $252$ as the sector for $L=10$, sampled here. The comparisons are numerical
observations, not theorems: we prove nothing about the spectrum of the random chain.

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

*Numerical observation:* for $L=10$ the mean ratio goes from $0.530$ ($W=0.5$), close to the sampled GOE value
$0.533$, through $0.484$ ($W=2$) and $0.410$ ($W=4$) to $0.384$ ($W=16$), close to $2\log2-1$, with standard
errors $0.003$–$0.004$. For $L=8$ the value at $W=0.5$ ($0.507$) is below the one at $W=1$ ($0.525$). At $W=0$ the
ring commutes with $T$ (Section 9.2) and with the total spin (Section 9.3), so its sector $N=L/2$ splits into
blocks; the fields couple these blocks only through $\sum_ih_is^z_i$, of operator norm at most $\frac{WL}2$.

> **Physical meaning** (physical interpretation, not a theorem of this course). Physicists read the change of
> $\langle r\rangle$ as a crossover from levels that repel like those of the GOE (weak disorder) to levels that behave
> like independent ones, attributed to **many-body localisation** of the eigenstates (strong disorder). The
> crossover drifts with $L$; where it lies as $L\to\infty$, and whether a localised phase exists at all, are open
> questions in the physics literature (D. Luitz, N. Laflorencie, F. Alet, *Phys. Rev. B* 91 (2015) 081103;
> J. Šuntajs, J. Bonča, T. Prosen, L. Vidmar, *Phys. Rev. E* 102 (2020) 062144). The dip at $L=8$, $W=0.5$ is read
> as a remnant of the extra block structure at $W=0$ (compare Section 5.4).

**Entanglement of free fermions.** For the eigenstates of the open XX chain the entanglement entropy of
the first $n$ sites is determined by an $n\times n$ matrix.

`entanglement_entropy` computes the reduced state sector by sector (it commutes with the magnetisation of
the block, Proposition 7.8) and its entropy. The ground state of sector $N$ is $\Phi$ with the $N$ lowest
modes, simple because the $\varepsilon_k$ are distinct:

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

*Numerical observation:* the entropies $0.6931,0.5199,0.7456,0.6283,0.7719,0.6560$ for $n=1,\dots,6$ oscillate
with the parity of $n$. Theorem 9.16 is due to I. Peschel (*J. Phys. A* 36 (2003) L205).

<a id="s-9-7"></a>

## 9.7 What to compare with

QuSpin (Python) builds spin-chain Hamiltonians in sectors of magnetisation, momentum, parity and spin
inversion and uses SciPy's sparse eigensolvers and matrix exponentials; dedicated exact-diagonalisation codes
(XDiag, SPINPACK, HΦ) are optimised for large sparse problems, and with lattice symmetries and large machines
such codes have reached about 40–50 spins; ITensor, TeNPy and TensorKit.jl (with $SU(2)$ fusion trees) treat
long chains with matrix-product states. What the library adds is modest and specific: the chain as a
`TypeIAlgebra` $\bigoplus_NM_{\binom LN}$, so that the states, channels, entropies and dynamics of Chapters 3–8
apply sector by sector; the inclusion tower with its index; $SU(2)$ blocks from highest-weight vectors with
Wigner–Eckart correlations; batches of disorder realisations in one call; and reduced states that never leave
the sectors. It has no reflection or spin-flip sectors, momentum sectors only for uniform XXZ chains, dense
methods to $L\approx16$, Krylov methods to $L\approx24$ (about $26$ on a large GPU), and no tensor networks.

<a id="x-9"></a>

## Exercises

[Exercises 1–3](../lecture_notes/09_spin_chains.md#x-9) of the lecture notes: the dimensions $d_{2,m}=2,1,2,1$ of the momentum sectors for
$L=4$, $N=2$; the sum rules $E_0=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle$ and $\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle=0$ for the
ground state of the Heisenberg ring; and `krylov.ground_state` started from the vector with all coefficients $1$
in the sector $N=L/2$, which returns $E=JL/4$, the top of the spectrum. Check them numerically.

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

**Next:** [Tutorial 10. Quantum computing](10_quantum_computing.md).
