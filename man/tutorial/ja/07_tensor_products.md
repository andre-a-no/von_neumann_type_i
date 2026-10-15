[English](../07_tensor_products.md) · [Español](../es/07_tensor_products.md) · [Français](../fr/07_tensor_products.md) · [Deutsch](../de/07_tensor_products.md) · [中文](../zh/07_tensor_products.md) · **日本語** · [Русский](../ru/07_tensor_products.md)

# チュートリアル 7. テンソル積とエンタングルメント

このチュートリアルは講義ノートの[第7章](../../lecture_notes/ja/07_tensor_products.md)に対応する。同章ではセクターをもつ二つの環の通常のテンソル積と融合テンソル積の
構造を証明し、部分トレースを埋め込みの双対として定義し、エンタングルメントとランダムな状態を調べる。ここではこれらの積を
ライブラリで構成し、部分トレースとエンタングルメント・エントロピーを計算し、定理を数値的に確かめる。

本章を通じて、$`\mathcal M_1=\bigoplus_{c=1}^{C_1}M_{k_c}(\mathbb C)`$ は
$`H_1=\bigoplus_c\mathbb C^{k_c}`$ に極小中心射影 $`P_c`$ をもって作用し、
$`\mathcal M_2=\bigoplus_{d=1}^{C_2}M_{l_d}(\mathbb C)`$ は $`H_2=\bigoplus_d\mathbb C^{l_d}`$ に
極小中心射影 $`Q_d`$ をもって作用する。$`\mathrm{Tr}`$ は作用素が作用する空間のトレースであり、$`\mathcal M_i`$ 上では第3章の blunt トレースである。
$`\mathbb C^k\otimes\mathbb C^l`$ 上の $`X`$ に対し $`X_{(a,b),(a',b')}=\langle e_a\otimes e_b,X\,e_{a'}\otimes e_{b'}\rangle`$ と書き、
$`e_a\otimes e_b\mapsto e_{(a-1)l+b}`$ によって $`\mathbb C^k\otimes\mathbb C^l`$ を $`\mathbb C^{kl}`$ と同一視する
（クロネッカー規約）。したがって $`(A\otimes B)_{(a,b),(a',b')}=A_{aa'}B_{bb'}`$ である。

<a id="s-7-1"></a>

## 7.1 二つの環のテンソル積

*講義ノート：* [7.1節](../../lecture_notes/ja/07_tensor_products.md#s-7-1)：[定理 7.1](../../lecture_notes/ja/07_tensor_products.md#r-7-1)（$`\mathcal M_1\otimes\mathcal M_2`$ の構造）。

`tensor_product(alg1, alg2)` はこの環を構成する。セクター $`(c,d)`$ の番号は $`s=cC_2+d`$（コードでは番号は
$`0`$ から始まる）、電荷は $`(q_c,r_d)`$ であり、その内部で $`e_a\otimes e_b`$ の番号は $`a\,l_d+b`$ である。
`kron(A, B, alg12)` は $`A\otimes B`$ を計算する。セクターの順序は $`H_1\otimes H_2`$ のクロネッカー順序と基底の置換だけ異なる。この置換は下の `kron_rows` が明示的に構成し、それを用いて全行列のクロネッカー積と比較する。

```python
import math
import numpy as np
import torch
from torch_vn_algebra import (TypeIAlgebra, DensityMatrix, states, composite,
                              tensor_product, fused_tensor_product, kron, partial_trace)

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def full(X):
    """An operator of (+)_c M_{k_c} as a block-diagonal matrix on (+)_c C^{k_c}, shape (batch, D, D)."""
    return torch.stack([torch.block_diag(*[x[c, :k, :k] for c, k in enumerate(X.algebra.k_factors)])
                        for x in X.matrix])

def kron_rows(alg12):
    """For every sector of a product algebra, the positions of its basis vectors in H1 (x) H2 (Kronecker order)."""
    k1, k2 = alg12.tensor_factors[0].k_factors, alg12.tensor_factors[1].k_factors
    return [torch.tensor([(sum(k1[:c]) + a) * sum(k2) + sum(k2[:d]) + b
                          for c, d, off in pairs for a in range(k1[c]) for b in range(k2[d])])
            for pairs in alg12.product_blocks]              # pair blocks (c, d, offset), by offset

def full12(X):
    """An operator of a (plain or fused) product algebra as a matrix on H1 (x) H2 in Kronecker order."""
    f1, f2 = X.algebra.tensor_factors
    N = sum(f1.k_factors) * sum(f2.k_factors)
    out = torch.zeros(X.matrix.shape[0], N, N, dtype=X.matrix.dtype)
    for s, r in enumerate(kron_rows(X.algebra)):
        out[:, r[:, None], r[None, :]] = X.matrix[:, s, :len(r), :len(r)]
    return out

def random_operator(alg):
    return alg.from_blocks([torch.randn(k, k, dtype=dt) for k in alg.k_factors])

alg1 = TypeIAlgebra([1, 2], [1, 2], precision='double')            # C (+) M_2
alg2 = TypeIAlgebra([2, 3], [2, 3], precision='double')            # M_2 (+) M_3
alg12 = tensor_product(alg1, alg2)
print(alg12.k_factors, alg12.charges, sum(k * k for k in alg12.k_factors))     # 4 sectors, dimension 65
assert sum(k * k for k in alg12.k_factors) == 5 * 13                           # Theorem 7.1(b)

A, A2, B, B2 = random_operator(alg1), random_operator(alg1), random_operator(alg2), random_operator(alg2)
AB = kron(A, B, alg12)
assert close(full12(AB)[0], torch.kron(full(A)[0], full(B)[0]))               # Theorem 7.1(a)
assert close(AB @ kron(A2, B2, alg12), kron(A @ A2, B @ B2, alg12))
for tr in ('Tr_blunt', 'Tr_norm', 'tau_vN'):                                   # Theorem 7.1(c)
    assert (getattr(AB, tr)() - getattr(A, tr)() * getattr(B, tr)()).abs().max() < 1e-10
```

<a id="s-7-2"></a>

## 7.2 部分トレース

*講義ノート：* [7.2節](../../lecture_notes/ja/07_tensor_products.md#s-7-2)：$`\mathrm{Tr}_2`$ の定義、[命題 7.2](../../lecture_notes/ja/07_tensor_products.md#r-7-2)、[定理 7.3](../../lecture_notes/ja/07_tensor_products.md#r-7-3)（部分トレースはチャネルである）、[命題 7.4](../../lecture_notes/ja/07_tensor_products.md#r-7-4)（指数 $`m^2`$ の条件付き期待値）。

`partial_trace(X, keep=1)` は各ブロック $`\rho_{cd}`$ を $`(k_c,l_d,k_c,l_d)`$ の形に並べ替えてトレースをとる。
`composite.partial_trace_channel(alg12, keep=1)` は定理 7.3 のチャネルであり、そのクラウス作用素は
セクター $`(c,d)`$ からセクター $`c`$ への $`1_{k_c}\otimes\langle j|`$ である。両者を定義、独立な `numpy.einsum`、
および双対性と照合する。

```python
rho = states.random_density_matrix(alg12, batch_size=3)
rho1, rho2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
assert isinstance(rho1, DensityMatrix)                                      # states give states
R = full12(rho).numpy().reshape(3, 3, 5, 3, 5)                              # indices (batch, i, j, i', j')
assert np.abs(full(rho1).numpy() - np.einsum('xijkj->xik', R)).max() < 1e-12
assert np.abs(full(rho2).numpy() - np.einsum('xijil->xjl', R)).max() < 1e-12
one2 = alg2.identity()
assert (states.expectation(rho, kron(A, one2, alg12)) - states.expectation(rho1, A)).abs().max() < 1e-12
assert close(partial_trace(AB, keep=1), A * B.trace.item())                 # Proposition 7.2(b)

Ptr = composite.partial_trace_channel(alg12, keep=1)
print(Ptr)
assert Ptr.is_trace_preserving()
assert close(Ptr(rho), rho1)
assert close(Ptr.adjoint('blunt')(A), kron(A, one2, alg12))                 # Theorem 7.3: the dual
```

モジュール `torch_vn_algebra.inclusions` は包含とその条件付き期待値および指数を実装している。
`from_tensor_product(alg1, m)` は[命題 7.4](../../lecture_notes/ja/07_tensor_products.md#r-7-4) の包含 $`\mathcal M_1\otimes1\subset\mathcal M_1\otimes M_m`$ を、
包含行列 $`\Lambda=m1_{C_1}`$、指数 $`m^2`$ とともに構成する（[講義ノート、7.2節](../../lecture_notes/ja/07_tensor_products.md#s-7-2)）：

```python
from torch_vn_algebra import inclusions
m = 3
inc = inclusions.from_tensor_product(alg1, m)                # alg1 (x) 1 c alg1 (x) M_3
X = random_operator(inc.M)
for trace in ('Tr_blunt', 'Tr_norm', 'tau_vN'):
    E = inc.conditional_expectation(trace.replace('Tr_', ''))
    assert close(E(X), partial_trace(X, keep=1) * (1 / m))   # Proposition 7.4(a): the same E for every trace
    assert abs(getattr(inc.embed(E(X)), trace)() - getattr(X, trace)()).max() < 1e-10
    assert close(E.index(), inc.M.identity() * m**2)         # Proposition 7.4(c)

def unit(a, b):
    e = torch.zeros(m, m, dtype=dt)
    e[a, b] = 1
    return math.sqrt(m) * kron(alg1.identity(), inc.M.tensor_factors[1].from_blocks([e]), inc.M)

u = [unit(a, b) for a in range(m) for b in range(m)]         # the quasi-basis of Proposition 7.4(b)
assert close(sum((v @ inc.embed(E(v.adjoint() @ X))).matrix for v in u), X)
assert close(sum((v @ v.adjoint()).matrix for v in u), inc.M.identity() * m**2)
print(inc.matrix.tolist(), inc.component_indices())           # Lambda = 3 * 1_2, index 9 per component
```

<a id="s-7-3"></a>

## 7.3 電荷の融合

*講義ノート：* [7.3節](../../lecture_notes/ja/07_tensor_products.md#s-7-3)：[命題 7.5](../../lecture_notes/ja/07_tensor_products.md#r-7-5)（融合積）と[定理 7.6](../../lecture_notes/ja/07_tensor_products.md#r-7-6)（融合積の縮約状態）。

`fused_tensor_product(alg1, alg2, fuse)`（既定の `fuse` は加法）は、セクターを融合後の電荷で並べ、
セクター内の組 $`(c,d)`$ を辞書式に並べ、それぞれに基底 $`e_a\otimes e_b`$ を用いる。
例：ハードコア・ボソン（[講義ノート、7.3節](../../lecture_notes/ja/07_tensor_products.md#s-7-3)）。粒子数を保存する作用素は、2サイトでは $`M_1\oplus M_2\oplus M_1`$
（電荷 $`0,1,2`$）、3サイトでは $`M_1\oplus M_3\oplus M_3\oplus M_1`$ をなし、[命題 7.5(c)](../../lecture_notes/ja/07_tensor_products.md#r-7-5) により、それらの融合積は
5サイトの全占有数に対する $`\{N\}'`$ であって、$`n_q=\binom5q`$ である。

```python
bos1 = TypeIAlgebra([1, 2, 1], [1, 2, 1], precision='double', charges=[0, 1, 2])
bos2 = TypeIAlgebra([1, 3, 3, 1], [1, 3, 3, 1], precision='double', charges=[0, 1, 2, 3])
F = fused_tensor_product(bos1, bos2)                         # U(1): total particle number
P = tensor_product(bos1, bos2)                               # particle numbers conserved separately
print(F.charges, F.k_factors, sum(k * k for k in F.k_factors), sum(k * k for k in P.k_factors))
assert F.k_factors == [math.comb(5, q) for q in range(6)]
assert sum(k * k for k in F.k_factors) == 252 and sum(k * k for k in P.k_factors) == 120

N1, N2 = full(bos1.central([0.0, 1.0, 2.0]))[0], full(bos2.central([0.0, 1.0, 2.0, 3.0]))[0]
N = torch.kron(N1, torch.eye(8, dtype=dt)) + torch.kron(torch.eye(4, dtype=dt), N2)
X = full12(F.operator(torch.randn(1, F.C, F.k_max, F.k_max, dtype=dt)))[0]
assert close(X @ N, N @ X)                                   # Proposition 7.5(c)
A, B = random_operator(bos1), random_operator(bos2)
assert close(full12(kron(A, B, F)), full12(kron(A, B, P)))   # Proposition 7.5(b)
```

融合環の次元は $`\sum_q\binom5q^2=\binom{10}5=252`$（再びヴァンデルモンドの恒等式）、通常の積の次元は
$`(1+4+1)(1+9+9+1)=120`$ である（定理 7.1(b)）。定理 7.1(a) と命題 7.5(b) により、通常の積は融合積の元のうち組 $`(c,d)`$
についてブロック対角なものからなる。残りの $`132`$ 次元は、組のブロック $`(c,d)`$ を同じ $`H_q`$ の別の組のブロック
$`(c',d')`$, $`c'\ne c`$（したがって $`d'\ne d`$）へ写す作用素で張られる：これらは全電荷を保存しつつ各部分の電荷を変える、
すなわちここでは二つの部分の間で粒子を移動させる。

定理 7.6 のブロック公式により $`(\rho_1)_c=\sum_d\mathrm{tr}_{l_d}\,\rho_{(c,d),(c,d)}`$ である。ライブラリの
二つの実装はいずれも、任意の融合則に対してこの写像を計算するように作られている：`partial_trace` は
対角の組のブロック $`\rho_{(c,d),(c,d)}`$ だけのトレースをとり、`partial_trace_channel` は一つのセクター内で
$`c`$ が同じ組に別々のクラウス添字を用い、その双対がやはり $`A\mapsto A\otimes1`$ となるようにしている。コードは
両者を、下の二つの融合則と残す因子の両方の選び方について、`einsum` で独立に計算した
$`E_1(\mathrm{Tr}_{H_2}\,\rho)`$ と照合し、双対を $`A\otimes1`$ と照合する。
我々の例では、加法の場合 $`E_1`$ は何もしない。パリティ則 $`f(q,r)=q+r\bmod2`$ の場合には、`bos1` の電荷 $`0`$
と $`2`$ はどの $`r`$ とも同じように融合し、$`E_1`$ は本物のコヒーレンスを取り除く。

```python
G = fused_tensor_product(bos1, bos2, fuse=lambda q, r: (q + r) % 2)       # only the parity is conserved
print(G.k_factors)

def block_mask(alg):
    """Ones in the diagonal blocks of (+)_c C^{k_c}: X * block_mask(alg) = sum_c P_c X P_c."""
    return full(alg.operator(torch.ones(1, alg.C, alg.k_max, alg.k_max, dtype=dt)))

for alg, coherent in ((F, False), (G, True)):
    rho = states.random_density_matrix(alg, batch_size=2)
    R = full12(rho).reshape(2, 4, 8, 4, 8)
    for keep, X_H, alg_k in ((1, torch.einsum('xijkj->xik', R), bos1), (2, torch.einsum('xijil->xjl', R), bos2)):
        E_X = X_H * block_mask(alg_k)                                     # E(Tr_H rho), Theorem 7.6
        assert close(full(partial_trace(rho, keep=keep)), E_X)
        assert close(full(composite.partial_trace_channel(alg, keep=keep)(rho)), E_X)
        lost = (X_H - E_X).abs().max().item()                             # coherences removed by E
        print(f"keep={keep}: largest coherence removed: {lost:.3e}")
        assert (lost > 1e-3) if coherent else (lost == 0)
    assert close(composite.partial_trace_channel(alg, keep=1).adjoint('blunt')(A), kron(A, bos2.identity(), alg))
```

<a id="s-7-4"></a>

## 7.4 純粋状態のエンタングルメント

*講義ノート：* [7.4節](../../lecture_notes/ja/07_tensor_products.md#s-7-4)：[補題 7.7](../../lecture_notes/ja/07_tensor_products.md#r-7-7)（シュミット分解）と[命題 7.8](../../lecture_notes/ja/07_tensor_products.md#r-7-8)（対称性分解エンタングルメント）：$`E(\psi)=H(p)+\sum_cp_cE(\hat\psi_c)`$。ここで $`p_c`$ は第一の部分のセクター $`c`$ の確率である。

$`H(p)`$ を **粒子数エンタングルメント**、$`\sum_cp_cE(\hat\psi_c)`$ を **配置エンタングルメント** という。
5サイト上の $`2`$ 粒子のランダムな状態について、両側の部分トレースと $`4\times8`$ 行列 $`\Psi`$ の特異値から
$`E(\psi)`$ を計算する。

```python
s = F.charges.index(2)                                      # the sector of two particles, n_q = 10
n = F.k_factors[s]
psi = torch.randn(n, dtype=dt)
psi = psi / psi.norm()
mat = torch.zeros(1, F.C, F.k_max, F.k_max, dtype=dt)
mat[0, s, :n, :n] = torch.outer(psi, psi.conj())
rho = DensityMatrix(F, matrix=mat)
r1, r2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)

psi_H = torch.zeros(4 * 8, dtype=dt)
psi_H[kron_rows(F)[s]] = psi                                # psi as a vector of H1 (x) H2
lam = torch.linalg.svdvals(psi_H.reshape(4, 8)) ** 2        # squared Schmidt coefficients
lam = lam[lam > 1e-14]
E_psi = -(lam * lam.log()).sum().item()
assert abs(r1.von_neumann_entropy().item() - E_psi) < 1e-10
assert abs(r2.von_neumann_entropy().item() - E_psi) < 1e-10  # Lemma 7.7

def entropy(sigma):
    w = torch.linalg.eigvalsh(sigma).clamp(min=1e-300)
    return -(w * w.log()).sum().item()

p = r1.sector_probabilities()[0]                            # c = 0, 1, 2 particles in the first part
assert torch.allclose(r2.sector_probabilities()[0, [2, 1, 0]], p)   # d(c) = 2 - c
H_num = -(p * p.log()).sum().item()
S_conf = sum(p[c].item() * entropy(r1.matrix[0, c, :k, :k] / p[c]) for c, k in enumerate(bos1.k_factors))
print(f"E = {E_psi:.4f} = H(p) {H_num:.4f} + configurational {S_conf:.4f}")
assert abs(E_psi - H_num - S_conf) < 1e-10                  # Proposition 7.8
```

<a id="s-7-5"></a>

## 7.5 ランダムな状態とペイジの公式

*講義ノート：* [7.5節](../../lecture_notes/ja/07_tensor_products.md#s-7-5)：[定理 7.11](../../lecture_notes/ja/07_tensor_products.md#r-7-11)（Lubkin）、[系 7.12](../../lecture_notes/ja/07_tensor_products.md#r-7-12)、[定理 7.13](../../lecture_notes/ja/07_tensor_products.md#r-7-13)（ペイジの公式）、[命題 7.14](../../lecture_notes/ja/07_tensor_products.md#r-7-14)（誘導測度）。

`partial_trace` で縮約した $`\mathbb C^3\otimes\mathbb C^6`$ の $`10^4`$ 個の純粋状態と誘導測度について定理 7.11 と 7.13 を、
$`M_1\oplus M_2\oplus M_3`$, $`r=2`$ について命題 7.14 を確かめる。

```python
a, b, n_samples = 3, 6, 10_000
g = torch.randn(n_samples, a * b, dtype=dt)                 # complex Gaussian: E|g_i|^2 = 1
psi = g / g.norm(dim=-1, keepdim=True)
AB_alg = tensor_product(TypeIAlgebra([a], [a], precision='double'), TypeIAlgebra([b], [b], precision='double'))
pure = DensityMatrix(AB_alg, matrix=(psi[:, :, None] * psi.conj()[:, None, :])[:, None],
                     validate=False)                        # pure by construction
rho_A = partial_trace(pure, keep=1)
purity_exact = (a + b) / (a * b + 1)
page = sum(1 / k for k in range(b + 1, a * b + 1)) - (a - 1) / (2 * b)

def mean_se(x):
    return x.mean().item(), x.std().item() / math.sqrt(x.numel())

for name, rho_s in (('reduced pure states', rho_A),
                    ('induced measure', states.random_density_matrix(rho_A.algebra, 40_000, rank=b))):
    (m_p, e_p), (m_S, e_S) = mean_se(rho_s.purity()), mean_se(rho_s.von_neumann_entropy())
    print(f"{name}: E Tr rho^2 = {m_p:.4f} +- {e_p:.4f} (exact {purity_exact:.4f}), "
          f"E S = {m_S:.4f} +- {e_S:.4f} (Page {page:.4f})")
    assert abs(m_p - purity_exact) < 5 * e_p and abs(m_S - page) < 5 * e_S     # Theorems 7.11, 7.13
    assert math.log((a * b + 1) / (a + b)) < m_S < math.log(a)                # Corollary 7.12
assert math.log(a) - (a - 1) / (a * b + 1) - (a - 1) / (2 * b) <= page <= math.log(a) - (a - 1) / (2 * b)
assert abs(sum(1 / k for k in (3, 4)) - 1 / 4 - 1 / 3) < 1e-15                # two qubits: E S = 1/3

alg3 = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')
m_p, e_p = mean_se(states.random_density_matrix(alg3, 40_000, rank=2).purity())
print(f"E Tr rho^2 = {m_p:.4f} +- {e_p:.4f} (exact {26 / 78:.4f})")
assert abs(m_p - 26 / 78) < 5 * e_p                                          # Proposition 7.14
```

<a id="s-7-6"></a>

## 7.6 比較の対象

単一のヒルベルト空間上では、QuTiP（`Qobj.ptrace`, `entropy_vn`, `rand_ket`, ジニブル分布と階数を指定する `rand_dm`）と
`numpy`（上のように reshape と `einsum`）が通常の積について部分トレース、エントロピー、誘導状態を計算する
（上の `einsum` による検証は、そのような計算をライブラリと比較している）。保存電荷をもつテンソルネットワークのライブラリ――ITensor（量子数によるブロック疎テンソル）、
TeNPy（電荷付きの脚、電荷ごとのエンタングルメント・スペクトル）、そして非可換対称性に対しては TensorKit.jl（融合木）――は、
大きな系に対して同じ仕方でセクターを融合し、対称性分解エンタングルメントを与える。ライブラリが付け加えるのは代数的な層である：
積はセクターとトレースをもつフォン・ノイマン環であり、融合則はラベルの任意の写像であり、部分トレースは
各トレースに関する双対をもつチャネルで、融合積に対しては正しい制限を与え（定理 7.3, 7.6）、そしてそのすべてがバッチ化されている。
テンソルネットワークは備えておらず、`fused_tensor_product` が受け取るのは、可換な対称性の場合のように電荷の組ごとに結果が一つの規則である。
$`SU(2)`$ の場合のような非可換な規則――スピンの組 $`j_1,j_2`$ が複数の全スピン $`|j_1-j_2|,\dots,j_1+j_2`$ に合成される
（第4章、演習問題 2）――は、ここでの融合積ではない（ライブラリはスピン $`\frac12`$ 鎖に対する $`SU(2)`$ を別に扱う、第4章）。

<a id="x-7"></a>

## 演習問題

講義ノートの[演習問題](../../lecture_notes/ja/07_tensor_products.md#x-7)の数値的検証：

1. $`\tau_{\mathrm{vN}}`$ に関する $`\mathrm{Tr}_2`$ の双対
   `partial_trace_channel(alg12).adjoint('tau_vN')` が、$`Z=\bigoplus_dC_2l_d1_{l_d}`$ として $`A\mapsto A\otimes Z`$ であることを確かめよ。
2. （劣加法性.）`alg12` のランダムな状態について、相互情報量
   $`S(\rho_1)+S(\rho_2)-S(\rho)`$ が正であり、$`S(\rho\|\rho_1\otimes\rho_2)`$ に等しいことを確かめよ。
3. 二つのサイトで共有された一つの粒子 $`\psi=(|01\rangle+|10\rangle)/\sqrt2`$ について、
   $`E(\psi)=H(p)=\log2`$ を確かめよ。

<details>
<summary>解答</summary>

*演習問題 1.*

```python
dual_tau = composite.partial_trace_channel(alg12, keep=1).adjoint('tau_vN')
A = random_operator(alg1)
Z = alg2.central([alg2.C * l for l in alg2.k_factors])
assert close(dual_tau(A), kron(A, Z, alg12))
```

*演習問題 2.*

```python
rho = states.random_density_matrix(alg12, batch_size=4)
rho1, rho2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
mutual = rho1.von_neumann_entropy() + rho2.von_neumann_entropy() - rho.von_neumann_entropy()
assert (mutual > 0).all()
assert torch.allclose(mutual, states.relative_entropy(rho, kron(rho1, rho2, alg12)))
```

*演習問題 3.*

```python
site = TypeIAlgebra([1, 1], [1, 1], precision='double', charges=[0, 1])
S2 = fused_tensor_product(site, site)
assert S2.k_factors == [1, 2, 1]
mat = torch.zeros(1, 3, 2, 2, dtype=dt)
mat[0, 1] = 0.5                                       # |psi><psi| with psi = (|01> + |10>)/sqrt 2
r1 = partial_trace(DensityMatrix(S2, matrix=mat), keep=1)
assert torch.allclose(r1.sector_probabilities()[0], torch.tensor([0.5, 0.5], dtype=torch.float64))
assert abs(r1.von_neumann_entropy().item() - math.log(2)) < 1e-12
```

</details>

**次へ：** [チュートリアル 8](08_dynamics.md)：ライブラリにおけるシュレーディンガー方程式とリンドブラッド方程式。
