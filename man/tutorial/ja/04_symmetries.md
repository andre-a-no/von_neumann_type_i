[English](../04_symmetries.md) · [Español](../es/04_symmetries.md) · [Français](../fr/04_symmetries.md) · [Deutsch](../de/04_symmetries.md) · [中文](../zh/04_symmetries.md) · **日本語** · [Русский](../ru/04_symmetries.md)

# チュートリアル 4. 対称性と表現

このチュートリアルは講義ノートの[第4章](../../lecture_notes/ja/04_symmetries.md)に対応する。同章では対称性の環の背後にある表現論を展開し、
それを任意の行列のリストが生成するフォン・ノイマン環を分解するランダム化アルゴリズムに変え（[命題 4.9](../../lecture_notes/ja/04_symmetries.md#r-4-9)）、
$U(1)$ と $SU(2)$ に適用し、実行列によって生成される環のセクターを実型・複素型・四元数型に分類する
（[定理 4.20](../../lecture_notes/ja/04_symmetries.md#r-4-20)）。ここでは指標を計算し、群と生成元に対して分解を実行し、`SpinChain` と `su2.SU2Chain` で
スピン鎖をセクターごとに構成し、分解の実構造を計算する。

<a id="s-4-2"></a>

## 4.2 等型分解と可換子環

*講義ノート：* [4.2節](../../lecture_notes/ja/04_symmetries.md#s-4-2)；[定理 4.3](../../lecture_notes/ja/04_symmetries.md#r-4-3)、[定理 4.4](../../lecture_notes/ja/04_symmetries.md#r-4-4)、[系 4.5](../../lecture_notes/ja/04_symmetries.md#r-4-5)。

講義ノートの[4.2節](../../lecture_notes/ja/04_symmetries.md#s-4-2)では、$S_3$ の三つの既約表現——自明表現、符号、および
$(1,1,1)^\perp\subseteq\mathbb C^3$ 上の標準表現——と、$(\mathbb C^2)^{\otimes3}$ の三つの因子を置換する表現 $\pi$ の指標
$\chi_\pi(\sigma)=2^{\#\text{cycles}(\sigma)}$ を手で計算し（下の表）、[系 4.5](../../lecture_notes/ja/04_symmetries.md#r-4-5) から重複度 $m_{\rm triv}=4$、$m_{\rm sgn}=0$、
$m_{\rm std}=2$ と $\dim\pi(S_3)'=20$ を得る。コードはこの表を計算し直す：$P_\sigma$ を
$|100\rangle,|010\rangle,|001\rangle$ の張る不変部分空間への $\pi(\sigma)$ の制限として求め（その上で $\pi(\sigma)$ はただ一つの $1$ の位置を置換する。補題 B.1(a)）、
$\#\text{cycles}(\sigma)=\dim\ker(P_\sigma-1)$（各巡回上で一定なベクトル）を用いる。

| | 単位元（$1$） | 互換（$3$） | $3$-巡回（$2$） |
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

これは2.4節の分解であり、今度は表現ごとに六つの数からなる表から得られた。

<a id="s-4-3"></a>

## 4.3 生成元から環へ：ランダム化アルゴリズム

*講義ノート：* [4.3節](../../lecture_notes/ja/04_symmetries.md#s-4-3)；[補題 4.6](../../lecture_notes/ja/04_symmetries.md#r-4-6)、[命題 4.7](../../lecture_notes/ja/04_symmetries.md#r-4-7)、[補題 4.8](../../lecture_notes/ja/04_symmetries.md#r-4-8)、[命題 4.9](../../lecture_notes/ja/04_symmetries.md#r-4-9)。

講義ノートのランダム化アルゴリズムは `decompose.wedderburn` であり、フロントエンドは `from_group`（有限群の全元。可換子環は
ツワール `decompose.twirl` でサンプリングする）と `from_generators`（任意の行列のリスト。可換子環は核 `decompose.commutant_basis`）である。
これは $\mathbb C^n$ をランダムな $K\in\mathcal A'_h$ の固有空間に分け、ランダムな $Y\in\mathcal A'$ のブロックを用いてそれらをまとめて整列させ、
ユニタリ `V`、$k_c$（`k_factors`）、$m_c$（`multiplicities`）をもつ `Decomposition` を返す。浮動小数点演算では、許容誤差より近い固有値を
等しいものとして扱い、得られた部分空間に可換子環の新しいランダムな元を圧縮することでそれらをさらに分割する
（[補題 4.8(a)](../../lecture_notes/ja/04_symmetries.md#r-4-8)）。

`Decomposition.verify()` は、生成元（および可換子環の新しい元）の $\mathcal B$（あるいは $\mathcal B'$）からの相対距離を報告し、$\sum_cm_c^2$ を $\dim\mathcal A'$（群に対しては指標から、系 4.5(c)。生成元に対しては核から）と比較する。$S_3$ と三つの量子ビットに対してこれを実行し、それから中身を見てみる：一つのランダムな $K$ の固有空間と、一つのランダムな $Y$ の絡作用素である。

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

計算量を支配するのは可換子環である：ツワールのサンプル一つに $n\times n$ 行列の積が $|G|$ 回かかる。生成元に対しては、
ライブラリはまず $\mathcal A$ のエルミートな元 $R=U\mathrm{diag}(\lambda)U^*$（生成元のエルミート部分のランダムな実線形結合に、
そのような二つの結合の対称化された積の定数倍を加えたもの）を対角化する。任意の $X\in\mathcal A'$ は $R$ と可換であり、したがってそのスペクトル射影とも可換である（定理 1.5(d)）から、
$U^*XU$ は $R$ の固有空間（次元 $d_\lambda$）に関してブロック対角である。これをあらかじめ課すことにより、核の問題の未知数は $n^2$ 個ではなく
$N_R=\sum_\lambda d_\lambda^2$ 個になる。いくつかの固有値をまとめてもブロックが大きくなるだけなので、縮小された問題は同じ解集合 $\mathcal A'$ をもつ。
分解において $V^*RV=\bigoplus_cr_c\otimes1_{m_c}$（エルミートな $r_c\in M_{k_c}$）である。$\lambda$ の $r_c$ における重複度を
$\mu_c(\lambda)$ とすると $d_\lambda=\sum_c\mu_c(\lambda)m_c$ であり、整数 $\mu\ge0$ に対して $\mu^2\ge\mu$ であるから
$$
N_R=\sum_\lambda\Big(\sum_c\mu_c(\lambda)m_c\Big)^2\ge\sum_{\lambda,c}\mu_c(\lambda)^2m_c^2\ge\sum_{\lambda,c}\mu_c(\lambda)m_c^2=\sum_ck_cm_c^2,
$$
であり、等号はすべての $r_c$ の固有値が単純かつ互いに相異なるときかつそのときに限り成り立つ。$L=3,4,5$ 個のスピンの全スピンについて、下のコードは等号を見いだす（*数値的観察*）。
最小値でさえ大きい：八つのスピン $\frac12$（$n=256$）では $\sum_ck_cm_c^2=4900$ であり、縮小された連立方程式は $(n^2+N_R)\times N_R$ 行列で、
$3.5\cdot10^8$ 個の複素成分をもつ。そのため4.6節ではスピン鎖に対して、答えをあらかじめ知っているモジュール `su2` を用いる。

```python
for n in (3, 4, 5):
    Dn = dc.from_generators(dc.total_spin(n))
    _, info = dc.commutant_basis(dc.total_spin(n), return_info=True)
    N_min = sum(k * m * m for k, m in zip(Dn.k_factors, Dn.multiplicities))
    assert info['unknowns'] >= N_min                                  # N_R >= sum_c k_c m_c^2
    print(n, info['unknowns'], N_min)                                # observed: equality
```

<a id="s-4-4"></a>

## 4.4 $U(1)$ と電荷セクター

*講義ノート：* [4.4節](../../lecture_notes/ja/04_symmetries.md#s-4-4)；[命題 4.10](../../lecture_notes/ja/04_symmetries.md#r-4-10)。

$n$ 個の量子ビットに対しては、電荷は個数作用素 $\hat N=\sum_jn_j$（$n_j=|1\rangle\langle1|$ は量子ビット $j$ 上に作用する）であり、
それを保存する作用素全体は $U(1)'\cong\bigoplus_{w=0}^nM_{\binom nw}(\mathbb C)$ をなす（[講義ノート、4.4節](../../lecture_notes/ja/04_symmetries.md#s-4-4)）。
これは `chains.SpinChain(n).algebra` である（ビット $1$ = 上向きスピン、$\hat N=S^z+\frac n2$。一方 `decompose.total_spin` では上向きスピンは
$\mathbb C^2$ の *第一の* 基底ベクトルであり、これは基底の番号を付け替えるがセクターは変えない）。任意の XXZ ハミルトニアンは $U(1)'$ に属する。
コードは $\hat N$ だけからセクターを求め、四つのスピンの XXZ 環をセクターごとに対角化する。

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

## 4.5 $SU(2)$：スピン

*講義ノート：* [4.5節](../../lecture_notes/ja/04_symmetries.md#s-4-5)；[定理 4.11](../../lecture_notes/ja/04_symmetries.md#r-4-11)、[命題 4.13](../../lecture_notes/ja/04_symmetries.md#r-4-13)、[定理 4.14](../../lecture_notes/ja/04_symmetries.md#r-4-14)。スピン $s$ は `decompose.spin_matrices(s)`、$L$ 個のスピン $\frac12$ の全スピンは `decompose.total_spin(L)` である。

$L=4$ に対して、[講義ノート、命題 4.13](../../lecture_notes/ja/04_symmetries.md#r-4-13) は重複度 $1,3,2$ をもつ $\{S_x,S_y,S_z\}''\cong M_5\oplus M_3\oplus M_1$
（次元 $35$）と、可換子環 $\{S_x,S_y,S_z\}'\cong M_1\oplus M_3\oplus M_2$（次元 $14$）を与える。ライブラリは三つの生成元だけから同じものを
見つけ（セクターは $(k_c,m_c)$ の順に並べられる）、スピン $S$ 上で $\mathbf S^2=S(S+1)$ となる（[定理 4.11](../../lecture_notes/ja/04_symmetries.md#r-4-11)）：

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

[講義ノート、定理 4.14](../../lecture_notes/ja/04_symmetries.md#r-4-14)（シューア・ワイル双対性）により、サイトの置換は全スピンの可換子環を張る：
$\mathrm{span}\,\pi(S_L)=\{S_x,S_y,S_z\}'$。$L\le5$ について、コードは二つの空間を数値的に比較する：定理の (a) により
$\mathrm{span}\,\pi(S_L)\subseteq\{S_x,S_y,S_z\}'$ であるから、それらの次元が等しければ十分である。

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

出力される次元は $2,5,14,42$ であり、カタラン数 $\binom{2L}L/(L+1)$ である。系 B.10 はすべての $L$ について
$\dim\{S_a\}'=\sum_Sm_S^2=\binom{2L}L/(L+1)$ を証明している。最後の行は双対性が働いている様子である：
$\pi(S_4)''$ はブロック $k=(1,2,3)$ と重複度 $(5,1,3)$ をもつ——全スピンの組 $(k_c,m_c)$ の役割を入れ替えたものである。一般に（系 B.9）スピン $S$ の重複度空間は次元 $m_S$ の $S_L$ の既約表現を担い、これらは互いに非同値である。ここでは $S=2,1,0$ に対してその次元は $1,3,2$ である。

<a id="s-4-6"></a>

## 4.6 スピン系

*講義ノート：* [4.6節](../../lecture_notes/ja/04_symmetries.md#s-4-6)；[命題 4.15](../../lecture_notes/ja/04_symmetries.md#r-4-15)、[定理 4.16](../../lecture_notes/ja/04_symmetries.md#r-4-16)、[系 4.17](../../lecture_notes/ja/04_symmetries.md#r-4-17)。

[講義ノート、系 4.17](../../lecture_notes/ja/04_symmetries.md#r-4-17) により、$L$ 個のスピン $\frac12$ の最高ウェイトベクトルの正規直交基底は、
[定理 4.16](../../lecture_notes/ja/04_symmetries.md#r-4-16) のクレブシュ・ゴルダン係数を用いて一サイトずつ構成できる。各基底ベクトルに中間スピンの列
$S_1=\frac12,S_2,\dots,S_L=S$ でラベルをつけると、`su2.SU2Chain` が構成するヤマノウチ基底が得られる（そのソースはまさにこれらの係数を用いる。
下のコードは系 4.17 の結論だけを確かめる）。それは `SpinChain` の $U(1)$ セクター $\hat N=\frac L2+S$ の内部にある（`SpinChain` のサイトの順序は
テンソル因子の付け替えであり、補題 B.1(b) により全スピンと可換である）。コードは定理 4.16 を密な行列とライブラリの `su2.clebsch_gordan` に
対して確かめ、それから八つのスピンの基底について系 4.17 の結論を確かめる。また[命題 4.15(a)](../../lecture_notes/ja/04_symmetries.md#r-4-15) を `D4` を用いて
四つのスピンの環で確かめ、八つのスピンについてブロックの大きさを比較する。

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

八つのスピンの環では、$U(1)$ のブロックの大きさは $\binom8w=1,8,28,56,70,56,28,8,1$（4.4節）、$SU(2)$ のブロックの大きさは
$S=0,\dots,4$ に対して $m_S=14,28,20,7,1$ である（命題 4.13）。対角化のコスト $\sum k^3$ は $739162$ から $33040$ へ、
$22$ 分の1に減る。環の次元 $\sum k^2$ は
$\sum_w\binom8w\binom8{8-w}=\binom{16}8=12870$（$16$ 個の元から $8$ 個を選び、そのうち $w$ 個を最初の $8$ 個から選ぶ）から
$\sum_Sm_S^2=1430=\binom{16}8/9$ に減る（系 B.10）。定理 9.9 により（環の長さは偶数で $J>0$）、基底状態はただ一つの一重項であり、
$E_0(0)<E_0(1)$ である。*数値的観察*（上のコードによる計算）：基底エネルギーは $-3.651093$ であり、第一励起準位は $-3.128419$ の三重項である。

<a id="s-4-7"></a>

## 4.7 実代数：実型、複素型、四元数型

*講義ノート：* [4.7節](../../lecture_notes/ja/04_symmetries.md#s-4-7)。複素共役で閉じた $*$-環 $\mathcal A$ と、列ブロック $V_c$ がセクターに属する分解 $V$ に対して、
[定理 4.20](../../lecture_notes/ja/04_symmetries.md#r-4-20) は $\bar P_c=P_{\sigma(c)}$ となるセクターの置換 $\sigma$、$\overline{V_c}=V_{\sigma(c)}(J_c\otimes M_c)$ となる
ユニタリ $J_c$、$M_c$、および $\sigma(c)=c$ のとき $J_c\bar J_c=\varepsilon_c1$ となる符号 $\varepsilon_c$ を与える；型と実代数 $\mathcal A_{\mathbb R}$ は
[定理 4.21](../../lecture_notes/ja/04_symmetries.md#r-4-21)、実シューアの補題は[定理 4.22](../../lecture_notes/ja/04_symmetries.md#r-4-22) であり、可換子環の実元は
[命題 4.23](../../lecture_notes/ja/04_symmetries.md#r-4-23) で記述される。

`Decomposition.real_structure()` は、`wedderburn`（`from_generators`、`from_group`）が返す分解に対して $\sigma$、$J_c$、$M_c$ と型を計算する。
保存されている生成元の共役が
$V(\bigoplus_cM_{k_c}\otimes1)V^*$ に属することを確かめ（相対残差 `'conjugation'`）、ユニタリ $P=V^*\bar V$ を作る。その列ブロック
$c$ は $V^*\overline{V_c}$ である：定理 4.20(a) により、それはセクター $\sigma(c)$ の行の上にのみ台をもち、そこで
$J_c\otimes M_c$ に等しい。ライブラリはブロックの最大部分を担うセクターを $\sigma(c)$ とし、
$(k_{\sigma(c)},m_{\sigma(c)})=(k_c,m_c)$ を要求し、残りを評価する（`'leak'`）；ブロックを、成分
$T_{(a,b),(j,l)}=T_{(a,j),(b,l)}$ をもつ $k_c^2\times m_c^2$ 行列に並べ替える。$T=J\otimes M$ に対してこれは階数一の行列
$\mathrm{vec}\,J\,(\mathrm{vec}\,M)^T$ であり、その主特異ベクトルから $J_c$ と $M_c$ を読み取り（$\|J_c\|_F^2=k_c$ とする）、ブロックと $J_c\otimes M_c$ の距離を記録する（`'kron'`）。これら三つの残差のいずれかが
許容誤差を超えるか、形が異なるか、相手が対合をなさないならば、`ValueError` を送出する（環が共役で閉じていない）。$\sigma(c)=c$ に対しては $\varepsilon_c$ を $\mathrm{Re}\,\mathrm{Tr}(J_c\bar J_c)/k_c$ の符号とし——
厳密な算術ではちょうど $\pm1$ であるから、二つの候補は $2$ 離れている——、これらのセクターにわたる
$\|J_c\bar J_c-\varepsilon_c\|_F/\sqrt{k_c}$ と
$\|M_c\bar M_c-\varepsilon_c\|_F/\sqrt{m_c}$ の最大値を `'sign'` として記録する。ただし、ある $\varepsilon_c=-1$ で $k_c$ または $m_c$ が奇数ならば
（定理 4.20(b) により不可能）、それを少なくとも $1$ に引き上げる。結果は `types`、`partner`（$\sigma$）、
`J`、`M`、`residuals`、`ok` をもち、`ok` はすべての残差が許容誤差（既定では分解の受理
閾値）以下であることを表す；そのとき $|\mathrm{Tr}(J_c\bar J_c)/k_c-\varepsilon_c|\le\|J_c\bar J_c-\varepsilon_c\|_F/\sqrt{k_c}$
（コーシー–シュワルツ）も許容誤差以下である。

`Decomposition.real_basis()` は実直交 $n\times n$ 行列 $O$ と、その列ブロックを記述するリスト（型、セクター、
列の範囲）を返す：
実型のセクターに対しては定理 4.21(a) の列 $O_c=V_c(W_J\otimes W_M)$ であり、実な $\alpha_c(a)$ により $O_c^TaO_c=\alpha_c(a)\otimes1_{m_c}$
となる；四元数型のセクターに対しては、$J_c\otimes M_c=WW^T$ かつ
$(J_c\otimes M_c)\bar W=W$ を満たす $W$ による列 $V_cW$ であり（$(J_c\otimes M_c)\overline{(J_c\otimes M_c)}=\varepsilon_c^2=1$ であるから $\varepsilon=1$ の補題 4.19(c)）、
これらは定理 4.21(a) の計算により実である；そして対 $c,\sigma(c)$（$c<\sigma(c)$）に対しては、列 $\sqrt2\mathrm{Re}\,V_c$ に続いて
$\sqrt2\mathrm{Im}\,V_c$ である。これらは正規直交である：$V_c^*\overline{V_c}=V_c^*V_{\sigma(c)}(J_c\otimes M_c)=0$ と
$V_c^T\overline{V_c}=\overline{V_c^*V_c}=1$ およびそれらの共役から $(\mathrm{Re}\,V_c)^T\mathrm{Re}\,V_c=(\mathrm{Im}\,V_c)^T\mathrm{Im}\,V_c=\frac12$ と
$(\mathrm{Re}\,V_c)^T\mathrm{Im}\,V_c=0$ が得られる；異なるブロックの列は、$\sigma$ の異なる軌道 $\mathcal O$ の
直交射影 $P_{\mathcal O}=\sum_{e\in\mathcal O}P_e$ の値域に属する（$\overline{V_c}=V_{\sigma(c)}(J_c\otimes M_c)$ であるから）ので、$O$ は直交行列である。実な $a$ に対し、$x=\pi_c(a)\otimes1_{m_c}$ として $aV_c=V_cx$ の実部と虚部をとると
$$
O^TaO=\begin{pmatrix}\mathrm{Re}\,x&\mathrm{Im}\,x\\-\mathrm{Im}\,x&\mathrm{Re}\,x\end{pmatrix}\quad\text{on these columns.}
$$
となる。対称なユニタリ $J=X+iY$ の $WW^T=J$ を満たす高木因子 $W$ は、ランダムな $t$ に対する $X+tY$ の直交
対角化から計算される（$J\bar J=JJ^*=1$ の虚部は $YX-XY$ であるから、$X$ と $Y$ は可換な実対称行列である）；ライブラリは残差 $\|WW^T-J\|_F/\|J\|_F$ と
$\|\mathrm{Im}(V_cW)\|_F/\sqrt{k_cm_c}$ によって結果を受理し、そのいずれかが
$\max(\text{tolerance},10^3\,\epsilon_{\rm mach}\,n)$ を超えると警告する（`DecompositionWarning`）。

講義ノートの[4.7節](../../lecture_notes/ja/04_symmetries.md#s-4-7)では、実行列によって生成される三つの環を調べる：(i) $\mathbb R^3\otimes\mathbb R^2$ 上の
$L_a\otimes1_2$（$(L_a)_{bc}=-\epsilon_{abc}$ は $\mathbb R^3$ の回転の生成元）：実型のセクターが一つ、$k=3$、$m=2$；(ii)
$\mathbb R^2\otimes\mathbb R^2$ 上の $\Omega_2\otimes1_2$：複素型の対、$k=1$、$m=2$；(iii) $\mathbb H^2=\mathbb R^8$ 上の四元数による左乗法
$1_2\otimes L_q$：四元数型のセクターが一つ、$k=2$、$m=4$。コードはそれぞれの分解と実構造を計算し、定理 4.20(a),(b) と、実型（定理 4.21(a)）
および対に対する $O^TaO$ のブロック形を確かめる。

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

コードは三つの例の既約な実加群 $\mathbb R^3$、$\mathbb R^2$、$\mathbb R^4$ の上で定理 4.22 を確かめる（ライブラリが見つける型と重複度、
すなわち $m=1$ の実型、$m=1$ の複素型の対、$m=2$ の四元数型により、定理 4.22(c) からこれらは既約である）。その際、実
可換子環を実線形方程式 $Xs=sX$、$Xs^T=s^TX$ の零空間として計算する。次に各例で可換子環のランダムな実対称元に対して命題 4.23 を確かめる。

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

## 4.8 比較の対象

有限群については、GAP と SageMath が指標表を厳密に計算し、系 4.5 により重複度が得られる。RepLAB（MATLAB；D. Rosset, F. Montealegre-Mora, J.-D. Bancal）は有限群およびコンパクト群の表現を数値的に分解し、これも可換子環のランダムな元を用いる。
SymbolicWedderburn.jl（Julia）は半正定値計画の対称性による簡約のためにウェダーバーン分解を計算する。物理では、QuSpin や他の厳密対角化コードが
$U(1)$ と格子の対称性によってブロック対角化を行い、ITensor と TeNPy はテンソルネットワークにおいてアーベル的な電荷を扱い、
TensorKit.jl は融合木によって $SU(2)$ を扱う。ライブラリの `decompose` は小さな行列（$n$ は数百程度まで）のための密な手法である。任意の行列のリスト（群の元または生成元）を受け付け、結果を証明書によって確認し（命題 4.9）、環とその可換子環を `to_algebra`、`from_algebra`、`to_commutant`、`from_commutant` をもつ `TypeIAlgebra`
オブジェクトとして返すので、第3章、第5–8章の内容がただちに適用できる。スピン $\frac12$ の鎖については、`su2.SU2Chain` が $SU(2)$ のブロックを直接与える（系 4.17）。受け付けるスピンの数は
$L=18$ までであり、この上限は基底のメモリによって決まる。

<a id="x-4"></a>

## 演習問題

講義ノートの[演習問題](../../lecture_notes/ja/04_symmetries.md#x-4)を数値的に確かめる（下の解答はこの順序である）：
1. $d=2$、$L=3,4$ に対する $\dim\pi(S_L)'=\binom{d^2+L-1}L$（値 $20$ と $35$）；
2. `from_generators` による $1\otimes\frac12\cong\frac12\oplus\frac32$ と $1\otimes1\cong0\oplus1\oplus2$；
3. 四つの量子ビット上の $e^{i\theta\hat N}$ の可換子環、$\theta=1$ の場合（$0<|q-q'|\le4<2\pi$ であるから、電荷 $0,\dots,4$ は相異なる
   $e^{iq}$ を与える）と $\theta=\pi$ の場合；
4. ゼーマン分裂：八つのスピンのハイゼンベルク環に対する $H+hS_z$ のスペクトルを、`su2.SU2Chain` と $U(1)$ セクターとで計算する；
5. `real_structure` による $Q_8$ の正則表現のセクター。

<details>
<summary>解答</summary>

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

**次へ：** [チュートリアル 5](05_random_matrices.md)：環の中のランダム行列——ライブラリにおけるハール・ユニタリ、ワインガルテン計算、準位統計。
