[English](../12_equivariant_learning.md) · [Español](../es/12_equivariant_learning.md) · [Français](../fr/12_equivariant_learning.md) · [Deutsch](../de/12_equivariant_learning.md) · [中文](../zh/12_equivariant_learning.md) · **日本語** · [Русский](../ru/12_equivariant_learning.md)

# チュートリアル 12. 同変機械学習

このチュートリアルは講義ノートの[第12章](../../lecture_notes/ja/12_equivariant_learning.md)に対応する。同章は、有限群の二つの表現の間の同変線形写像がブロック行列であること（[定理 12.4](../../lecture_notes/ja/12_equivariant_learning.md#r-12-4)）、
群上の平均がそれらの上への直交射影であること（[定理 12.5](../../lecture_notes/ja/12_equivariant_learning.md#r-12-5)）を証明し、不変量を数える。
ここでは `torch_vn_algebra.decompose` を用いて、グラフの頂点特徴量と辺特徴量に作用する $S_4$ に対し同変層を構成・学習し、
これらの結果を数値的に確かめる。

コードでは、ベクトルのバッチは形状 `(batch, n)` のテンソルであり、行列 $T$ はそれに `x @ T.T` として作用する。

<a id="s-12-1"></a>

## 12.1 同変線形写像

*講義ノート：* [12.1節](../../lecture_notes/ja/12_equivariant_learning.md#s-12-1)：[補題 12.1](../../lecture_notes/ja/12_equivariant_learning.md#r-12-1)（軌道）と[系 12.2](../../lecture_notes/ja/12_equivariant_learning.md#r-12-2)（DeepSets、ベル数）。

DeepSets 層（Zaheer et al. 2017）は[講義ノート、系 12.2](../../lecture_notes/ja/12_equivariant_learning.md#r-12-2)(a) であり、グラフネットワークのための Maron et al. (2019) の基底は (b) である。
次のコードは軌道を総当たりで数え、数値的に計算した可換子環と比較し（[補題 12.1](../../lecture_notes/ja/12_equivariant_learning.md#r-12-1)）、
$S_4$ の二つの表現を分解する。

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

したがって、`dc.from_group` が計算した分解において（assert がその出力を確かめる）、$\mathbb C^{16}$ は次元 $1$ の表現を二回——これは自明表現であり、その重複度は命題 12.7(a) と補題 12.1 により $[4]^2$ 上の軌道の個数 $2$ である（12.4節のコードがこのセクターを特定する）——、次元 $2$ の表現を一回、次元 $3$ の二つの表現をそれぞれ一回と三回含む：$16=1\cdot2+2\cdot1+3\cdot1+3\cdot3$ であり、$15=2^2+1^2+1^2+3^2$ は補題 12.1 による軌道の数え上げおよび $\dim\pi(G)'=\sum_cm_c^2$（定理 2.6）と一致する。二つの表現の場合も、すべては一つの可換子環に帰着する。

[講義ノート、補題 12.3](../../lecture_notes/ja/12_equivariant_learning.md#r-12-3) により、$\pi=\pi_1\oplus\pi_2$ に対して $\mathrm{Hom}_G(\pi_1,\pi_2)=Q_2^*\,\pi(G)'\,Q_1$ であり、
[定理 12.4](../../lecture_notes/ja/12_equivariant_learning.md#r-12-4) により、これは適当な基底のもとで $\bigoplus_cM_{b_c\times a_c}(\mathbb C)$ であって、その次元は $\sum_\rho m_\rho(\pi_1)\,m_\rho(\pi_2)$ である。

<a id="s-12-2"></a>

## 12.2 ライブラリによる同変層

*講義ノート：* [12.2節](../../lecture_notes/ja/12_equivariant_learning.md#s-12-2)：[定理 12.5](../../lecture_notes/ja/12_equivariant_learning.md#r-12-5)（ツワール）。層は[補題 12.3](../../lecture_notes/ja/12_equivariant_learning.md#r-12-3)と[定理 12.4](../../lecture_notes/ja/12_equivariant_learning.md#r-12-4)の証明に従う。

定理 2.6 の分解がわかれば、定理 12.4 の証明はそのままアルゴリズムである。`dc.from_group` は 4.3節の乱択法によってその分解を計算する（厳密な算術では確率 1 で正しい。ここでは assert が計算された出力を確かめる）。$\pi_1=\pi_2$ のときは定理 2.6 そのものである：重みはブロック
$w_c\in M_{m_c}$、すなわち抽象的な可換子環 `D.commutant_algebra()` $=\bigoplus_cM_{m_c}$ の作用素であり、層は $x\mapsto V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*x$ で、`D.from_commutant` によって計算される。以下のモジュールは本章のコードであり、ライブラリの一部ではない。

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

$\pi_1\ne\pi_2$ の場合は証明を一行ずつたどる：$\pi_1\oplus\pi_2$ を分解し、`to_commutant` で $P_1$ のブロック $p_c$ を計算し、それらを対角化して左下の隅を残す。$p_c$ のトレースは
$a_c=m_c(\pi_1)$ である。

$n=4$ の系 12.2(b) の場合、計算された分解において（下で出力し assert で確かめる）：$\pi_1\oplus\pi_2$ の四つのセクター（次元 $1,2,3,3$。最初のものは自明表現である。系 12.2(a) により $\pi_1$ の一次元の構成要素は自明表現のみだからである）において、頂点特徴量 $\pi_1$ は $a=(1,0,0,1)$、辺特徴量 $\pi_2$ は
$b=(2,1,1,3)$ をもつから、系 12.2(b) が証明するとおり $\dim\mathrm{Hom}_G=1\cdot2+1\cdot3=5=B_3$ である。

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
b = tuple(m - a for m, a in zip(node_to_edge.D.multiplicities, node_to_edge.a))
assert node_to_edge.D.k_factors == (1, 2, 3, 3) and node_to_edge.a == [1, 0, 0, 1] and b == (2, 1, 1, 3)
corners = node_to_edge.D.basis('commutant')[:, n:, :n].reshape(-1, 16 * n)   # Lemma 12.3: Q_2^* pi(G)' Q_1
assert torch.linalg.matrix_rank(corners).item() == bell[3]
```

層を同変にする第二の方法は、群上で平均することである。**ツワール** は次で定義される。
$$
\mathcal P(W)=\frac1{|G|}\sum_{g\in G}\pi_2(g)\,W\,\pi_1(g)^*,\qquad W\in M_{n_2\times n_1}(\mathbb C).
$$

[講義ノート、定理 12.5](../../lecture_notes/ja/12_equivariant_learning.md#r-12-5) により、$\mathcal P$ は $\mathrm{Hom}_G(\pi_1,\pi_2)$ へのヒルベルト・シュミット内積に関する直交射影であり、
$\dim\mathrm{Hom}_G(\pi_1,\pi_2)=|G|^{-1}\sum_g\overline{\chi_{\pi_1}(g)}\chi_{\pi_2}(g)$ である。$\pi_1=\pi_2$ のとき、`dc.commutant_dimension` はこの
指標公式を評価する。パラメータ化 $W\mapsto\mathcal P(W)$ は $n_1n_2-\dim\mathrm{Hom}_G$ 個のパラメータを核に浪費し（辺特徴量では $256-15$ 個）、密な $\pi(g)$ に対して一回の
評価に $O(|G|n_1n_2(n_1+n_2))$ の計算量がかかり、$S_n$ では $|G|=n!$ である。ライブラリは同じ射影を分解から
`D.conditional_expectation(W, 'commutant')` として、四回の $n\times n$ 行列積で計算する。
コードは `dc.twirl` と `D.conditional_expectation` が一致すること、$\mathcal P$ が $\pi(G)'$ 上で恒等写像であり
自己共役であること、および頂点特徴量から辺特徴量への次元公式を確かめる。

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

<a id="s-12-3"></a>

## 12.3 制約付きの層とスペクトル層

*講義ノート：* [12.3節](../../lecture_notes/ja/12_equivariant_learning.md#s-12-3)：[命題 12.6](../../lecture_notes/ja/12_equivariant_learning.md#r-12-6)。

`EquivariantLinear` の重み $w=(w_c)_c$ は $*$-環の中にあるので、第1章と第11章のすべての構成がそれらに適用できる。$w\in\bigoplus_cM_{m_c}$ に対して $\Phi(w)=V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*$ と書く。

[講義ノート、命題 12.6](../../lecture_notes/ja/12_equivariant_learning.md#r-12-6) により、(a) $\Phi$ は $\bigoplus_cM_{m_c}$ から $\pi(G)'$ の上への $*$-同型であり、$\Phi(w)$ が
自己共役、ユニタリ、あるいは正であるための必要十分条件は、すべての $w_c$ がそうであることである。(b) 正規な $w$ に対し $f(\Phi(w))=\Phi\big((f(w_c))_c\big)$ である。(c)
同変なユニタリは $\Phi(u_0e^{s-s^*})$ にほかならない。

したがって直交（ユニタリ）な同変層は、`D.commutant_algebra()` 上の `optimize.UnitaryParam` に $\Phi$ を続けたものである。(c) により、このパラメータ化はすべての同変ユニタリに到達し、厳密な算術では勾配法の各反復はユニタリかつ同変である（コードは学習後に両方を $10^{-12}$ の精度で確かめる）。同様に `PositiveParam` は正値な同変層を与え、それは第11章と同様に正規化される（スペクトルは $[\ell,1]$ にあり最大値は $1$。(a) により、一般の正値層は $\Phi(t^*t)$ である）。また、拡散 $e^{-tL}$ のような、自己共役な同変生成子 $L$ のスペクトル層 $f(L)$ には、ブロック $w_c$ の固有値分解が必要であり、その計算量は $n^3=4096$ ではなく $\sum_cm_c^3$ 程度（ここでは $37$）である。ライブラリは `apply_function` を通して微分する（その逆伝播については何も証明しない。下のコードは一つの導関数を確かめる）。次のコードは同変なユニタリを Adam で $32$ 組の入出力の組にフィットさせ——これはこの実行の*数値的観察*であり、定理ではない：損失は $10^3$ 倍より大きく減少し、それを assert が確かめるが、収束に関する主張は一切しない——、(b) と $t$ に関する勾配を、密行列および公式 $\frac d{dt}e^{-tL}=-Le^{-tL}$ と照合する（$L=U\mathrm{diag}(\lambda)U^*$ に対し、定理 1.5 により $e^{-tL}=U\mathrm{diag}(e^{-t\lambda_i})U^*$ であり、$U$ は $t$ によらない）。

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

<a id="s-12-4"></a>

## 12.4 不変量

*講義ノート：* [12.4節](../../lecture_notes/ja/12_equivariant_learning.md#s-12-4)：[命題 12.7](../../lecture_notes/ja/12_equivariant_learning.md#r-12-7)（不変ベクトルと二次形式）と[命題 12.8](../../lecture_notes/ja/12_equivariant_learning.md#r-12-8)（モリーンの公式）。

グラフに対しては、不変な辺特徴量は対角線とその補集合の指示ベクトルで張られ（$B_2=2$）、不変な線形読み出しは対角成分の和と非対角成分の和である。次のコードはブロック $\rho_c(g)$（`to_algebra`）から自明なセクターを見つけ、(a)、(b) および (c) の二通りの数え上げを確かめる。

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

$\mathbb C^4$ 上では、それらは $\sum_ix_i^2$ と $\sum_{i<j}x_ix_j$ である。(c) の証明はあらゆる次数で通用する。

[講義ノート、命題 12.8](../../lecture_notes/ja/12_equivariant_learning.md#r-12-8)（モリーンの公式）により、
$\sum_{d\ge0}\dim\mathbb C[x]^G_d\,t^d=|G|^{-1}\sum_g\det\big(1-t\pi(g)\big)^{-1}$ であり、置換表現に対しては
$\dim\mathbb C[x]^G_d$ は $d$ 次の単項式上の $G$ の軌道の個数である。

（これは 1897 年のモリーンの定理である。さらに学ぶには B. Sturmfels, *Algorithms in Invariant Theory*, Springer 2008, 2.2節を参照。）次のコードは右辺を固有値ごとに一つの等比級数を掛けることで $t^4$ まで展開し、単項式の軌道と比較する。$\mathbb C^4$ 上の $S_4$ では、個数 $1,1,2,3,5$ は $d$ を高々四つの部分に分ける分割の個数である（多重集合の軌道は、その重複度を並べ替えたものによって定まる）。

```python
def molien(G, d_max):
    """Coefficients of |G|^{-1} sum_g det(1 - t pi(g))^{-1} up to t^d_max: products of geometric series."""
    lam = torch.linalg.eigvals(G)                                # (|G|, n), |lambda| = 1
    c = torch.zeros(G.shape[0], d_max + 1, dtype=dt)
    c[:, 0] = 1
    for i in range(lam.shape[1]):                                # multiply by 1 / (1 - t lambda_i)
        for d in range(1, d_max + 1):
            c[:, d] = c[:, d] + lam[:, i] * c[:, d - 1]
    return c.mean(0)

def monomial_orbits(G, d):
    """Orbits of G on the monomials of degree d, i.e. on the multisets of d points of Omega."""
    act = [g.real.argmax(0).tolist() for g in G]                 # pi(g) e_x = e_{act[x]}
    mons = itertools.combinations_with_replacement(range(G.shape[-1]), d)
    return len({min(tuple(sorted(a[x] for x in m)) for a in act) for m in mons})

for G in (G_node, G_edge):
    series, counts = molien(G, 4), [monomial_orbits(G, d) for d in range(5)]
    print("Molien series:", [f"{v.real:.4f}" for v in series], " orbits of monomials:", counts)
    assert (series - torch.tensor(counts, dtype=dt)).abs().max() < 1e-9          # Proposition 12.8
assert [monomial_orbits(G_node, d) for d in range(5)] == [1, 1, 2, 3, 5] and monomial_orbits(G_edge, 2) == orbits2
```

<a id="s-12-5"></a>

## 12.5 比較の対象

特定の群に対する同変層は成熟したライブラリによって提供されている。**escnn**（Cesa, Lang and Weiler 2022、e2cnn の後継）は、既約表現と核の調和基底から $O(2)$ と $O(3)$ の部分群に対するステアラブル畳み込みを構成し、**e3nn**（Geiger and Smidt 2022）は球面調和関数とクレブシュ–ゴルダンのテンソル積から $E(3)$ 同変なネットワークを構成する。**EMLP**（Finzi, Welling and Wilson 2021）は任意の行列群を扱う：離散的な生成元に対して $\pi_2(g)T=T\pi_1(g)$ を、リー環の生成元に対して $d\pi_2(X)T=Td\pi_1(X)$ を課し、解の基底を数値的な零空間として計算する——これは `dc.commutant_basis` の線形問題である。群畳み込み（Cohen and Welling 2016）は、$G$ が左からの積によってそれ自身に作用する場合の、補題 12.1 で $\Omega_1=\Omega_2=G$ としたものである：$(y,x)$ と $(y',x')$ が同じ軌道にあるための必要十分条件は $x^{-1}y=x'^{-1}y'$ であり（このとき $g=x'x^{-1}$ が一方の組を他方に写す）、したがって同変写像は $(Tv)(y)=\sum_xf(x^{-1}y)v(x)$、すなわち $G$ 上の関数 $f$ との畳み込みである。Kondor and Trivedi (2018) はコンパクト群の等質空間に対する同様の問題を研究しているが、本講義ではそれを扱わない。置換については、DeepSets（Zaheer et al. 2017）と Maron et al. (2019) の基底が系 12.2 である。

本ライブラリは、格子上の畳み込み、連続空間上のステアラブル核、$\mathbb R$ 上の分解（実型、複素型、四元数型。演習問題 2 は実行列を扱う）、あるいは大きな $n$ を提供しない：$V$ は密な $n\times n$ 行列であり、すべての元で与えられた群に対して分解の計算量はランダムサンプル一つあたり $O(|G|n^3)$ である（リー環の生成元は `dc.from_generators` を通る。4.3節）。本ライブラリが付け加えるのは解空間の **構造** である：可換子環の基底は行列のリストにすぎないが、分解は可換子環を $*$-環 $\bigoplus_cM_{m_c}$ と同一視する。こうして積、随伴、ユニタリ性、正値性、関数計算が小さな行列の上でブロックごとに作用し（命題 12.6）、ライブラリはそれらを通して微分し、`to_commutant`/`from_commutant` が二つの描像の間を翻訳する。

<a id="x-12"></a>

## 演習問題

講義ノートの[演習問題](../../lecture_notes/ja/12_equivariant_learning.md#x-12) 1、3、4 の数値的検証：

- *演習問題 1*（チャネル付きの DeepSets）：$\mathbb C^n\otimes\mathbb C^{d_1}$ から $\mathbb C^n\otimes\mathbb C^{d_2}$ への同変写像が
  $1\otimes A+J\otimes B$ にほかならないことを `EquivariantMap` で確かめよ。
- *演習問題 3*（非線形性）：座標ごとの ReLU は置換表現と可換であるが、$\mathbb R^2$ の $\pi/2$ 回転とは
  可換でない。$x\mapsto\varphi(\|x\|)x$ はすべてのユニタリ表現と可換である。
- *演習問題 4*：辺の表現では、一つのガウス・サンプルでは同変層の $2$ 次元の空間が決定されずに残り、二つのサンプルで
  層が決定される。

<details>
<summary>解答</summary>

*演習問題 1.*

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

*演習問題 3.*

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

*演習問題 4.*

```python
basis = D_edge.basis('commutant')                                # (15, 16, 16), Hilbert-Schmidt orthonormal
for N, lost in ((1, 2), (2, 0)):
    Xs = torch.randn(16, N, dtype=dt)
    M = torch.stack([(B @ Xs).reshape(-1) for B in basis], dim=1)               # w -> (T x_1, ..., T x_N)
    assert torch.linalg.matrix_rank(M).item() == 15 - lost
```

</details>

**次へ：** [チュートリアル 13. ジョルダン代数：量子力学はなぜ *-環を用いるのか](13_jordan_algebras.md)。
