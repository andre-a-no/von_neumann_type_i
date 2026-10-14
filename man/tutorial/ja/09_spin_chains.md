[English](../09_spin_chains.md) · [Español](../es/09_spin_chains.md) · [Français](../fr/09_spin_chains.md) · [Deutsch](../de/09_spin_chains.md) · [中文](../zh/09_spin_chains.md) · **日本語** · [Русский](../ru/09_spin_chains.md)

# チュートリアル 9. スピン鎖

このチュートリアルは講義ノートの[第9章](../../lecture_notes/ja/09_spin_chains.md)に対応する。同章では XXZ 鎖のセクター——磁化
（[命題 9.1](../../lecture_notes/ja/09_spin_chains.md#r-9-1)）、運動量（[命題 9.5](../../lecture_notes/ja/09_spin_chains.md#r-9-5)）、全スピン——を構成し、それらが基底状態について何を語るかを証明する
（[定理 9.7](../../lecture_notes/ja/09_spin_chains.md#r-9-7)、[定理 9.9](../../lecture_notes/ja/09_spin_chains.md#r-9-9)）。ここでは `SpinChain` と `SU2Chain` で鎖を構成し、自由フェルミオンとベーテ仮説で
ライブラリを検証し、密な行列には大きすぎるセクターに `krylov` のランチョス法とクリロフ法を用い、準位統計と
エンタングルメントを計算する。

本章を通じて、サイトは $`0,\dots,L-1`$ とし、$`(\mathbb C^2)^{\otimes L}`$ の基底は配置
$`|b\rangle=|b_0\dots b_{L-1}\rangle`$ からなり、上向きスピンに対して $`b_i=1`$ とする（ライブラリの規約：ラベルの
第 $`i`$ ビットがサイト $`i`$ である）。$`s^a_i`$（$`a=x,y,z`$）はサイト $`i`$ 上のスピン行列 $`\frac12\sigma_a`$、
$`s^\pm_i=s^x_i\pm is^y_i`$、$`n_i=s^z_i+\frac12`$、$`S^a=\sum_is^a_i`$ とし、$`N=S^z+\frac L2`$ は上向きスピンの
個数である。鎖は **開いた** 鎖（ボンド $`(i,i+1)`$、$`i<L-1`$）か **周期的** な鎖（ボンド $`(L-1,0)`$ も含む。このとき
$`L\ge3`$ で、添字は $`L`$ を法とする）である。$`J>0`$（反強磁性）、$`\hbar=1`$ とし、$`\log`$ は自然対数である。

<a id="s-9-1"></a>

## 9.1 XXZ 鎖とその $`U(1)`$ セクター

*講義ノート：* [9.1節](../../lecture_notes/ja/09_spin_chains.md#s-9-1)；[命題 9.1](../../lecture_notes/ja/09_spin_chains.md#r-9-1)（$`H`$ は $`S^z`$ を保存する）、[命題 9.2](../../lecture_notes/ja/09_spin_chains.md#r-9-2)、[補題 9.3](../../lecture_notes/ja/09_spin_chains.md#r-9-3)、[系 9.4](../../lecture_notes/ja/09_spin_chains.md#r-9-4)（鎖の塔とその指数）。

`SpinChain(L).xxz` は各セクターの配置（`chain.labels[N]`）から直接ブロックを構成する。
これを、クロネッカー積から構成した完全な $`2^L\times2^L`$ 行列と比較する。

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

鎖の塔 $`\mathcal A_0\subset\mathcal A_1\subset\cdots`$：`inclusions.chain_tower` が包含を、
`inclusions.bratteli_diagram` がそのブラッテリ図を、`Inclusion.index` が指数 $`\|\Lambda\|^2`$ を与える。コードは
パスカルの法則、[補題 9.3](../../lecture_notes/ja/09_spin_chains.md#r-9-3)、および $`[\mathcal A_{L+1}:\mathcal A_L]=4\cos^2\frac\pi{2L+4}`$（[系 9.4](../../lecture_notes/ja/09_spin_chains.md#r-9-4)）を確かめる。

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

## 9.2 並進、運動量、マーシャルの符号則

*講義ノート：* [9.2節](../../lecture_notes/ja/09_spin_chains.md#s-9-2)；[命題 9.5](../../lecture_notes/ja/09_spin_chains.md#r-9-5)（運動量セクター）、[定理 9.7](../../lecture_notes/ja/09_spin_chains.md#r-9-7)（マーシャルの符号則）、[系 9.8](../../lecture_notes/ja/09_spin_chains.md#r-9-8)（基底状態の運動量）。

`xxz_momentum` は、各軌道から代表元 $`a`$ を一つずつとった正規化された基底 $`|a,k\rangle`$ でブロックを構成する。
軌道を独立に数え、スペクトルを密なセクターと比較する。

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

$`924`$ 個の状態からなるセクター $`N=6`$ は、大きさ $`80,75,78,76,78,75,80,75,78,76,78,75`$ のブロックに分かれる。

マーシャルの符号則（[定理 9.7](../../lecture_notes/ja/09_spin_chains.md#r-9-7)）により、ランダム磁場中の開いた鎖と周期的な鎖について
$`(-1)^{N_B(b)}\psi(b)`$ は各セクターで一つの符号をもち、[系 9.8](../../lecture_notes/ja/09_spin_chains.md#r-9-8) により環のセクター $`N`$ の最低準位は
$`N`$ が偶数なら運動量ブロック $`k=0`$ に、奇数なら $`k=\pi`$ にある：

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

## 9.3 $`\Delta=1`$ における $`SU(2)`$

*講義ノート：* [9.3節](../../lecture_notes/ja/09_spin_chains.md#s-9-3)；[定理 9.9](../../lecture_notes/ja/09_spin_chains.md#r-9-9)（リーブ–マティス）。$`\Delta=1`$ かつ $`h=0`$ のとき $`H\cong\bigoplus_S1_{2S+1}\otimes H_S`$ であり、`SU2Chain` は最高ウェイトベクトルの上でブロック $`H_S`$ を構成する。

十二個のスピンについてブロックの大きさを比較し、定理 9.9 を $`U(1)`$ セクターで確かめ、ギャップを計算し、さらに
鎖ではない二部グラフ（$`n_A=5`$、$`n_B=3`$、$`S_*=1`$、ランダムな結合）で定理 9.7 と 9.9 を確かめる。

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

$`SU(2)`$ ブロックの大きさは $`m_S=132,297,275,154,54,11,1`$ で、最大の $`U(1)`$ ブロックの $`924`$ と比べて小さく、
対角化のコスト $`\sum k^3`$ は $`2.0\cdot10^9`$ から $`5.3\cdot10^7`$ に減る。*数値的観察*：$`L`$ とギャップの積はゆっくりと増加し、
$`L=6,\dots,12`$ に対して $`4.108,4.181,4.232,4.270`$ である。$`L\to\infty`$ でのそのふるまいについては何も証明しない。

**相関.** 各サイトのスピンはベクトル作用素をなし、その球面成分は階数一のテンソル作用素をなす
（補題 B.15(a),(b)）。**ウィグナー–エッカートの定理**（定理 B.16）により、二つの多重項の間のそれらの行列要素は
クレブシュ–ゴルダン係数と、磁気量子数によらない **換算行列要素** との積であり、
$`\mathbf s_i\cdot\mathbf s_j=\sum_q(-1)^qs^{(i)}_qs^{(j)}_{-q}`$ である（補題 B.15(a)）。`spin_correlation` はこれを用いて、
最高ウェイトベクトルだけから各 $`H_S`$ 上の $`\mathbf s_i\cdot\mathbf s_j`$ を構成する。`'direct'` 法は代わりにボンド演算子を
最高ウェイトベクトルへ圧縮する。$`\mathbf s_i\cdot\mathbf s_j\in\{S_a\}'`$ であるから、これは命題 4.15(b) により正しい。
第一の方法の実装が正しいことはここでは証明しない。コードは、両者が互いに一致し、かつ $`U(1)`$ セクターで
書き下した基底状態上で定義どおりに評価した値とも一致することを確かめる。

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

相関 $`\langle\mathbf s_0\cdot\mathbf s_r\rangle=0.75,-0.4489,0.1878,-0.1658,0.1209,-0.1224,0.1070,\dots`$ は
符号が交互に変わる。これは定理 9.7 と 9.9 がすべての偶数の環に対して強制することである（[講義ノート、9.3節](../../lecture_notes/ja/09_spin_chains.md#s-9-3)）。

<a id="s-9-4"></a>

## 9.4 厳密な検証：自由フェルミオンとベーテ仮説

*講義ノート：* [9.4節](../../lecture_notes/ja/09_spin_chains.md#s-9-4)；[定理 9.10](../../lecture_notes/ja/09_spin_chains.md#r-9-10) により、開いた XX 鎖のセクター $`N`$ のエネルギーは相異なる $`N`$ 個の $`\varepsilon_k=J\cos\frac{\pi k}{L+1}`$ の和である；[命題 9.12](../../lecture_notes/ja/09_spin_chains.md#r-9-12)（二マグノン、ベーテ方程式）。

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

[命題 9.12](../../lecture_notes/ja/09_spin_chains.md#r-9-12) により、$`K=2\pi n/L`$、$`c=\cos\frac K2`$ として
$`F(q)=(-1)^ne^{iqL}(c-\Delta e^{-iq})+c-\Delta e^{iq}`$ の各根 $`q`$ は環のセクター $`N=2`$ の固有ベクトル $`\Psi`$ を与える
（[講義ノート、9.4節](../../lecture_notes/ja/09_spin_chains.md#s-9-4)）。コードは根を $`q\in(0,\pi)`$ で（散乱状態）、また直線
$`q=i\kappa`$、$`\pi+i\kappa`$（$`\kappa>0`$）の上で（束縛状態）探し、$`L=9`$、$`\Delta=0.6`$ に対して各根について $`\Psi`$ を構成し、
$`36`$ 個の状態からなるセクター $`N=2`$ と比較する：

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

*数値的観察*：$`L=9`$、$`\Delta=0.6`$ に対して探索は $`36`$ 個の根、すなわち散乱状態 $`34`$ 個と束縛状態 $`2`$ 個を見つけ、
それらのベクトルは（丸め誤差の範囲で）一次独立な固有ベクトルである：$`36`$ 個の固有状態すべてがベーテ方程式から
得られる。他の $`L`$ について完全性は証明されておらず、偶数の $`L`$ と $`\Delta\ne0`$ では完全性は成り立たない
（[講義ノート、9.4節](../../lecture_notes/ja/09_spin_chains.md#s-9-4)）：ベクトル $`\Phi=\sum_x(-1)^x|x,x+1\rangle`$ はすべてのベーテ・ベクトルと直交し、
コードの最後の数行は $`L=8`$ に対して $`H\Phi=(E_F-J\Delta)\Phi`$ を確かめる。

<a id="s-9-5"></a>

## 9.5 ランチョス法とクリロフ法

*講義ノート：* [9.5節](../../lecture_notes/ja/09_spin_chains.md#s-9-5)；[命題 9.13](../../lecture_notes/ja/09_spin_chains.md#r-9-13)（ランチョス法：$`T_m=Q_m^*HQ_m`$ は三重対角）、[定理 9.14](../../lecture_notes/ja/09_spin_chains.md#r-9-14)（リッツ値、すなわち $`T_m`$ の固有値は $`\theta^{(m)}_j\ge\lambda_j`$ を満たし、$`m`$ とともに減少する）、[命題 9.15](../../lecture_notes/ja/09_spin_chains.md#r-9-15)（クリロフ法による時間発展）。

セクター $`N=L/2`$ の次元は $`L=16`$ で $`12870`$、$`L=20`$ で $`184756`$ であるが、$`H`$ は疎であり（各行に $`O(L)`$ 個の非零
成分、`xxz_sparse`）、クリロフ法は積 $`Hv`$ だけを必要とする。$`H`$ をサイズ $`n`$ のエルミート行列とし、
$`v\ne0`$ とする。**クリロフ空間** は $`\mathcal K_m=\mathrm{span}\{v,Hv,\dots,H^{m-1}v\}`$ である。

定理 9.14 は $`\theta^{(m)}_1`$ がどれだけ速く $`\lambda_1`$ に近づくかについては何もいわない。それは下で観察するだけである。
浮動小数点演算では $`q_j`$ は直交性を失う。`krylov.lanczos` は完全な再直交化を行い、コードは計算された $`q_j`$ が
丸め誤差の範囲で正規直交であり $`Q_m^*HQ_m=T_m`$ であることを確かめる。開いた XX 鎖については、セクター $`N`$ の
厳密な基底状態エネルギーは $`N`$ 個の最小の $`\varepsilon_k`$ の和である（定理 9.10(c)）。

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

*数値的観察*：最小のリッツ値の誤差は（定理 9.14 により当然ながら）単調に減少し、$`5`$、$`10`$、$`20`$ ステップの後で
それぞれ $`0.42`$、$`3.1\cdot10^{-3}`$、$`3.1\cdot10^{-6}`$ となり、$`40`$ ステップの後には丸め誤差の水準に達する。次元は $`924`$ である。
`krylov.ground_state` は、正規化されたリッツ対 $`(E,\psi)`$ の残差 $`\|H\psi-E\psi\|`$ が $`10^{-10}\max(1,|E|)`$ を下回るまで
（`v0` が与えられなければガウス型のランダムベクトルから）ランチョス法を再始動し（再始動は高々 $`20`$ 回で、そうならなければ
警告を出す）、残差を返す。小さな残差は近くに固有値があることを保証するが、それが最小のものであることは保証しない
（[講義ノート、9.5節](../../lecture_notes/ja/09_spin_chains.md#s-9-5)）ので、密な対角化と比較する。これが唯一の参照である（ここで
計算するので、値は数値的なものである）：$`L=12`$ まではセクター $`N=L/2`$、$`L=16`$ では次元 $`810`$ の運動量ブロック $`k=0`$ で
あり、後者は系 9.8 によりセクター $`N=8`$ の基底状態を含む。$`L=4`$ では値 $`-2J`$ は厳密である（[講義ノート、9.5節](../../lecture_notes/ja/09_spin_chains.md#s-9-5)）。

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

$`L=16`$（次元 $`12870`$）では、ランチョス法は一つの CPU スレッドで一秒を十分に下回る時間で終わる。*数値的観察*：
サイトあたりのエネルギー $`-0.5,-0.467129,-0.456387,-0.451545,-0.448949`$（$`L=4,\dots,12`$）は $`L`$ とともに増加する。
その極限については何も証明しない。

**時間発展.** $`e^{-i\tau H}v`$ に対しては $`\mathcal K_m`$ に射影する：$`w_\tau=\|v\|\,Q_me^{-i\tau T_m}e_1`$。

[命題 9.15](../../lecture_notes/ja/09_spin_chains.md#r-9-15) により、任意の $`\tau`$ に対して $`\|w_\tau\|=\|v\|`$ かつ $`\langle w_\tau,Hw_\tau\rangle=\langle v,Hv\rangle`$ であり、
$`\mathcal K_m`$ が $`H`$ 不変ならば $`w_\tau=e^{-i\tau H}v`$ である。

したがって、クリロフ法の各ステップは、正確であろうとなかろうとノルムとエネルギーを保存する。$`\mathcal K_m`$ が不変で
ない限り、命題 9.15 は精度については何もいわない。`krylov.evolve` は時間区間を短いステップに分割し、各ステップを
$`m=30`$ のクリロフ空間で行い、発見的な誤差指標が大きいときには警告を出す。それに対する誤差限界は何も証明せず、代わりに
密な固有値分解（定理 1.5）から計算した $`e^{-itH}`$ と比較する。ハイゼンベルク環の上で
ネール状態 $`|{\uparrow\downarrow\uparrow\downarrow}\cdots\rangle`$ からクエンチを行い、交替磁化
$`m_s=\frac2L\sum_i(-1)^i\langle s^z_i\rangle`$ を追跡し、$`L=12`$ については密な伝播子と比較する。

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

*数値的観察*：$`m_s`$ は $`1`$ から $`0.770`$、$`0.279`$ を経て $`t=2`$ で $`-0.183`$ まで減衰し、二つのサイズで
$`1.7\cdot10^{-5}`$ の範囲で一致する。$`t=4`$ では両者は異なる（$`0.026`$ 対 $`0.035`$）。

<a id="s-9-6"></a>

## 9.6 乱れ：準位統計とエンタングルメント

*講義ノート：* [9.6節](../../lecture_notes/ja/09_spin_chains.md#s-9-6)；[定理 9.16](../../lecture_notes/ja/09_spin_chains.md#r-9-16)（ペシェルの公式）。

**ランダム磁場ハイゼンベルク鎖** は、$`[-W,W]`$ 上の一様分布に従う独立な磁場 $`h_i`$ をもつ。これは $`S^z`$ を保存する
（命題 9.1）ので、準位統計はセクター $`N`$ ごとに計算する（定理 5.15 と 5.4節：セクターを混ぜると独立な準位を偽造して
しまう）。二つの磁場が異なれば $`THT^*-H=\sum_i(h_{i-1}-h_i)s^z_i\ne0`$ であるから、並進はもはや対称性ではない。
セクター $`N=L/2`$ の中央の半分で $`\langle r\rangle`$ を計算し、独立な実現にわたる標準誤差を一つ添え（各実現が平均比を
一つ与える）、二つの参照値と比較する：独立な指数分布の間隔の厳密な平均比 $`2\log2-1\approx0.3863`$（命題 5.14(a)）と、
$`L=10`$ のセクターと同じ大きさ $`252`$ の GOE 行列（5.3節の定義による。コードは GOE 行列の二倍である $`G+G^T`$ を
サンプルするが、これはどの比も変えない）の平均比で、後者はここでサンプルする。これらの比較は数値的観察であって
定理ではない：ランダムな鎖のスペクトルについては何も証明しない。

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

*数値的観察*：$`L=10`$ では平均比は $`0.530`$（$`W=0.5`$、サンプルした GOE の値 $`0.533`$ に近い）から $`0.484`$（$`W=2`$）、
$`0.410`$（$`W=4`$）を経て $`0.384`$（$`W=16`$、$`2\log2-1`$ に近い）まで変化し、標準誤差は $`0.003`$–$`0.004`$ である。$`L=8`$ では
$`W=0.5`$ での値（$`0.507`$）が $`W=1`$ での値（$`0.525`$）より小さい。$`W=0`$ では環は $`T`$（9.2節）および全スピン（9.3節）と
可換であるから、そのセクター $`N=L/2`$ はブロックに分かれる。磁場はこれらのブロックを、作用素ノルムが高々 $`\frac{WL}2`$ の
$`\sum_ih_is^z_i`$ を通じてのみ結合する。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. 物理学者は $`\langle r\rangle`$ の変化を、GOE の準位のように
> 反発する準位（弱い乱れ）から独立な準位のようにふるまう準位へのクロスオーバーと読み、後者を固有状態の
> **多体局在**（強い乱れ）に帰する。クロスオーバーは $`L`$ とともにずれていく。$`L\to\infty`$ でそれがどこにあるか、
> そもそも局在相が存在するかどうかは、物理学の文献における未解決問題である（D. Luitz, N. Laflorencie, F. Alet,
> *Phys. Rev. B* 91 (2015) 081103; J. Šuntajs, J. Bonča, T. Prosen, L. Vidmar, *Phys. Rev. E* 102 (2020) 062144）。
> $`L=8`$、$`W=0.5`$ でのくぼみは、$`W=0`$ における余分なブロック構造の名残と読まれる（5.4節と比較せよ）。

**自由フェルミオンのエンタングルメント.** 開いた XX 鎖の固有状態について、最初の $`n`$ 個のサイトの
エンタングルメント・エントロピーは $`n\times n`$ 行列によって決まる。

`entanglement_entropy` は縮約状態をセクターごとに計算し（それはブロックの磁化と可換である、
命題 7.8）、そのエントロピーを求める。セクター $`N`$ の基底状態は $`N`$ 個の最低のモードからなる $`\Phi`$ であり、
$`\varepsilon_k`$ が相異なるので単純である。

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

*数値的観察*：エントロピー $`0.6931,0.5199,0.7456,0.6283,0.7719,0.6560`$（$`n=1,\dots,6`$）は $`n`$ の偶奇とともに振動する。
定理 9.16 は I. Peschel（*J. Phys. A* 36 (2003) L205）による。

<a id="s-9-7"></a>

## 9.7 比較の対象

QuSpin（Python）は磁化、運動量、パリティ、スピン反転のセクターでスピン鎖のハミルトニアンを構成し、
SciPy の疎な固有値ソルバーと行列指数関数を用いる。専用の厳密対角化コード
（XDiag、SPINPACK、HΦ）は大規模な疎問題向けに最適化されており、格子の対称性と大型計算機を用いて
そのようなコードは約 40–50 スピンに達している。ITensor、TeNPy、TensorKit.jl（$`SU(2)`$ の融合木をもつ）は
行列積状態によって長い鎖を扱う。このライブラリが付け加えるものは控えめで具体的である：鎖を
`TypeIAlgebra` $`\bigoplus_NM_{\binom LN}`$ として扱い、第3–8章の状態、チャネル、エントロピー、ダイナミクスを
セクターごとに適用できること。指数をもつ包含の塔。最高ウェイトベクトルから得られる $`SU(2)`$ ブロックと
ウィグナー–エッカートによる相関。一度の呼び出しで乱れの実現のバッチを扱えること。そしてセクターから決して
出ない縮約状態。反転やスピン反転のセクターはなく、運動量セクターは一様な XXZ 鎖に限られ、密な
方法は $`L\approx16`$ まで、クリロフ法は $`L\approx24`$ まで（大型 GPU では約 $`26`$ まで）であり、テンソルネットワークはない。

<a id="x-9"></a>

## 演習問題

講義ノートの[演習問題 1–3](../../lecture_notes/ja/09_spin_chains.md#x-9)：$`L=4`$、$`N=2`$ に対する運動量セクターの次元 $`d_{2,m}=2,1,2,1`$；
ハイゼンベルク環の基底状態に対する和則 $`E_0=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle`$ と $`\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle=0`$；
そしてセクター $`N=L/2`$ ですべての係数が $`1`$ のベクトルから始めた `krylov.ground_state`。これはスペクトルの最上端
$`E=JL/4`$ を返す。これらを数値的に確かめよ。

<details>
<summary>解答</summary>

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

**次へ：** [チュートリアル 10. 量子計算](10_quantum_computing.md)。
