[English](../12_equivariant_learning.md) · [Español](../es/12_equivariant_learning.md) · [Français](../fr/12_equivariant_learning.md) · [Deutsch](../de/12_equivariant_learning.md) · [中文](../zh/12_equivariant_learning.md) · **日本語** · [Русский](../ru/12_equivariant_learning.md)

# 12. 同変機械学習

ニューラルネットワークの層が **同変** であるとは、入力を対称変換で変換してから層を適用しても、層を適用してから出力を変換しても同じ結果になることをいう：グラフの頂点の番号を付け替えると、その特徴量の番号も同じように付け替わる。線形層についてはこれは可換性の関係であるから、同変な層は可換子環をなし、第2章により、それは基底の取り替えを除いて環 $\bigoplus_cM_{m_c}(\mathbb C)$ である。本章では、それらのパラメータを数え、`torch_vn_algebra.decompose` を用いてそれらを構成・学習し、さらに環の構造を用いて、構成からして同変であるユニタリ層、正値層、スペクトル層を作る。

**本章で学ぶこと**：二つの表現の間の同変線形写像が、重複度で与えられる大きさのブロック行列であるのはなぜか。それらを軌道あるいは指標によってどう数えるか（DeepSets とグラフネットワークのベル数）。群上の平均化が同変写像への直交射影であるのはなぜか。関数計算とユニタリのパラメータ化が可換子環にブロックごとにどう作用するか。そして不変ベクトルと不変二次形式をどう数えるか。

本章を通じて、$G$ は有限群であり、$\mathbb C^n$ 上の $G$ の **表現** とはユニタリ値をとる準同型
$\pi:G\to M_n(\mathbb C)$ のことである（これは制約ではない：可逆な $\pi(g)$ は内積
$\sum_g\langle\pi(g)x,\pi(g)y\rangle$ に関してユニタリである。4.1節）。その **指標** は $\chi_\pi(g)=\operatorname{Tr}\pi(g)$ である。既約表現、同値な表現、および既約表現 $\rho$ の $\pi$ における **重複度** $m_\rho(\pi)$（$\rho$ が現れなければ $0$）は 2.4節のとおりであり、同節では $m_\rho(\pi)$ が分解によらないこととシューアの補題が証明されている。[第4章](04_symmetries.md) ではその両方を詳しく扱う（4.1–4.2節）。長方形行列の上で $\langle S,T\rangle=\operatorname{Tr}(S^*T)$ は **ヒルベルト・シュミット** 内積である。コードでは、ベクトルのバッチは形状 `(batch, n)` のテンソルであり、行列 $T$ はそれに `x @ T.T` として作用する。

## 12.1 同変線形写像

**定義.** $\pi_1,\pi_2$ を $\mathbb C^{n_1}$、$\mathbb C^{n_2}$ 上の $G$ の表現とする。行列
$T\in M_{n_2\times n_1}(\mathbb C)$ がすべての $g\in G$ に対して $T\pi_1(g)=\pi_2(g)T$ を満たすとき、$T$ は **同変** である（**絡作用素** である）という。それらはベクトル空間 $\operatorname{Hom}_G(\pi_1,\pi_2)$ をなし、$\pi_1=\pi_2=\pi$ のときこれは可換子環 $\operatorname{Hom}_G(\pi,\pi)=\pi(G)'$ である。

機械学習で最もよく現れる表現は座標を置換するものである。$G$ が有限集合 $\Omega$ に作用するとき、$\mathbb C^\Omega$ 上の **置換表現** は $\pi(g)e_x=e_{gx}$ である。例：集合の $n$ 個の点を置換する $S_n$（$\Omega=[n]=\{1,\dots,n\}$）、およびグラフの頂点の番号を付け替える $S_n$。後者は $\Omega=[n]^2$ に作用し、すなわち隣接行列に $X\mapsto P_\sigma XP_\sigma^T$ によって作用する。

**補題 12.1（軌道）.** $\pi_i$ を $\mathbb C^{\Omega_i}$ 上の $G$ の置換表現とする。このとき
$T\in\operatorname{Hom}_G(\pi_1,\pi_2)$ であるための必要十分条件は、すべての $g$、$y\in\Omega_2$、$x\in\Omega_1$ に対して $T_{gy,gx}=T_{yx}$ となることである。$\Omega_2\times\Omega_1$ 上の $G$ の軌道の指示行列は基底をなし、$\dim\operatorname{Hom}_G(\pi_1,\pi_2)$ はこれらの軌道の個数である。

*証明.* $\pi_2(g)T\pi_1(g)^*e_x=\pi_2(g)Te_{g^{-1}x}=\sum_yT_{y,g^{-1}x}e_{gy}$ であるから、$\pi_2(g)T\pi_1(g)^*$ の $(gy,x)$ 成分は $T_{y,g^{-1}x}$ である。よって $T$ が同変であること、すなわちすべての $g$ に対し $\pi_2(g)T\pi_1(g)^*=T$ であることは、すべての $g,y,x$ に対し $T_{gy,x}=T_{y,g^{-1}x}$ であることと同値であり、これは（$x$ を $gx$ で置き換えれば）述べた条件、すなわち $T$ が軌道上で一定であることにほかならない。軌道の指示行列は互いに交わらない台をもつので線形独立であり、軌道上で一定な関数はすべてそれらの線形結合である。$\square$

**系 12.2.** (a)（*DeepSets.*）$\mathbb C^n$ 上の $S_n$（$n\ge2$）に対し、可換子環は $\operatorname{span}\{1,J\}$（$J$ はすべての成分が $1$ の行列）であり、$\mathbb C^n=\mathbb C\mathbf 1\oplus\mathbf 1^\perp$ は互いに非同値な二つの既約表現、すなわち自明表現と次元 $n-1$ の **標準** 表現の和である。(b)（*グラフネットワーク.*）$\mathbb C^{[n]^p}$ から $\mathbb C^{[n]^q}$ への（$n$ 頂点上の $p$ 階および $q$ 階のテンソルの間の）同変写像の空間の次元は、$p+q$ 個の元からなる集合を高々 $n$ 個のブロックに分ける分割の個数に等しい。$n\ge p+q$ のとき、その個数は **ベル数** $B_{p+q}$ である（$B_2=2$、$B_3=5$、$B_4=15$）。

*証明.* (b) 二つの組 $x,y\in[n]^{p+q}$ が同じ軌道にあるための必要十分条件は、それらが同じ **等号パターン**、すなわち $i\sim j\Leftrightarrow x_i=x_j$ による位置の集合 $\{1,\dots,p+q\}$ の分割をもつことである。パターンが不変であることは明らかである。逆に $x,y$ が同じパターンをもてば、$x_i\mapsto y_i$ は $\{x_i\}$ から $[n]$ への矛盾なく定義された単射であり、これは $\sigma x=y$ を満たす置換 $\sigma$ に拡張される。高々 $n$ 個のブロックをもつ分割はすべてパターンである（各ブロックに相異なる値を与えればよい）。あとは補題 12.1 により結論が従う。(a) $p=q=1$ のとき、$[n]^2$ 上の軌道は対角線とその補集合の二つであり、その指示行列は $1$ と $J-1$ である。よって $\dim\pi(S_n)'=2=\sum_cm_c^2$（定理 2.6、命題 2.7）であり、これは $m_c=1$ の二つのセクター、すなわち互いに非同値な二つの既約な直和成分を強制し、$\pi(S_n)'$ はそれらの射影 $P_1,P_2$ で張られる。不変な直線 $\mathbb C\mathbf 1$ への射影は $\pi(S_n)'$ に属する（補題 2.1）から、それは $P_1$ または $P_2$ である：自明表現を担うこの直線が一つのセクターであり、$\mathbf 1^\perp$ がもう一つのセクターである。$\square$

DeepSets 層 $x\mapsto\lambda x+\gamma\,(\mathbf 1^Tx)\mathbf 1$（Zaheer et al. 2017）は (a) であり、グラフネットワークのための Maron et al. (2019) の基底は (b) である。次のコードは軌道を総当たりで数え、数値的に計算した可換子環と比較し、$S_4$ の二つの表現を分解する。

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

したがって $\mathbb C^{16}$ は自明表現を二回、次元 $2$ の表現を一回、次元 $3$ の二つの表現をそれぞれ一回と三回含む：$16=1\cdot2+2\cdot1+3\cdot1+3\cdot3$、$15=2^2+1^2+1^2+3^2$。二つの表現の場合も、すべては一つの可換子環に帰着する。

**補題 12.3（直和）.** $\pi=\pi_1\oplus\pi_2$ を $\mathbb C^{n_1}\oplus\mathbb C^{n_2}$ 上の表現とし、等長な埋め込みを $Q_i$、射影を $P_i=Q_iQ_i^*$ とする。このとき $P_1,P_2\in\pi(G)'$、
$\operatorname{Hom}_G(\pi_1,\pi_2)=Q_2^*\,\pi(G)'\,Q_1$ であり、$T\mapsto Q_2TQ_1^*$ は
$\operatorname{Hom}_G(\pi_1,\pi_2)$ から $P_2\,\pi(G)'P_1$ の上への線形全単射である。

*証明.* $\pi(g)$ はブロック対角であるから $\pi(g)Q_i=Q_i\pi_i(g)$ かつ $Q_i^*\pi(g)=\pi_i(g)Q_i^*$ であり、したがって
$P_i\pi(g)=Q_i\pi_i(g)Q_i^*=\pi(g)P_i$ である。$X\in\pi(G)'$ ならば
$Q_2^*XQ_1\pi_1(g)=Q_2^*X\pi(g)Q_1=Q_2^*\pi(g)XQ_1=\pi_2(g)Q_2^*XQ_1$ である。$T$ が同変ならば
$\pi(g)Q_2TQ_1^*=Q_2\pi_2(g)TQ_1^*=Q_2T\pi_1(g)Q_1^*=Q_2TQ_1^*\pi(g)$ であるから $Q_2TQ_1^*\in\pi(G)'$ であり、これは
$P_2(Q_2TQ_1^*)P_1$ に等しく、$Q_2^*(Q_2TQ_1^*)Q_1=T$ である。逆に、$X\in\pi(G)'$ による任意の $Y=P_2XP_1$ は
$\pi(G)'$ に属し、$Y=Q_2(Q_2^*YQ_1)Q_1^*$ である。$\square$

**定理 12.4（同変写像はブロック行列である）.** $\rho_1,\dots,\rho_C$ を $\pi_1\oplus\pi_2$ に現れる互いに非同値な既約表現とし、$k_c=\dim\rho_c$、$a_c=m_{\rho_c}(\pi_1)$、$b_c=m_{\rho_c}(\pi_2)$ とする。$V_i^*\pi_i(g)V_i=\bigoplus_c\rho_c(g)\otimes1$ を満たすユニタリ
$V_1:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{a_c}\to\mathbb C^{n_1}$ と
$V_2:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{b_c}\to\mathbb C^{n_2}$ が存在して
$$
\operatorname{Hom}_G(\pi_1,\pi_2)=\Big\{V_2\Big(\bigoplus_{c=1}^C1_{k_c}\otimes y_c\Big)V_1^*:\ y_c\in M_{b_c\times a_c}(\mathbb C)\Big\},
$$
となり、$(y_c)_c\mapsto T$ は単射である。したがって $\operatorname{Hom}_G(\pi_1,\pi_2)\cong\bigoplus_cM_{b_c\times a_c}(\mathbb C)$ であり、
$$
\dim\operatorname{Hom}_G(\pi_1,\pi_2)=\sum_\rho m_\rho(\pi_1)\,m_\rho(\pi_2),
$$
である。ここで和は同値を除いたすべての既約表現にわたる。$\pi_1=\pi_2$ のときこれは
$\pi(G)'\cong\bigoplus_cM_{m_c}$（命題 2.7）であり、$\operatorname{Hom}_G(\pi_1,\pi_2)=0$ であるための必要十分条件は、両方に現れる既約表現が存在しないことである。

*証明.* 定理 2.6 と命題 2.7 を $\pi=\pi_1\oplus\pi_2$ に適用する：
$V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$ かつ $V^*\pi(G)'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$ を満たすユニタリ $V$ が存在する。
$X\mapsto V^*XV$ は $*$-同型であり、$P_1\in\pi(G)'$（補題 12.3）であるから、射影 $p_c$ により $V^*P_1V=\bigoplus_c1\otimes p_c$ と書ける。$u_c^*p_cu_c=1_{a_c}\oplus0_{b_c}$ となるユニタリ $u_c$ を選び（定理 1.3）、ここではさしあたり
$a_c=\operatorname{rank}p_c$、$b_c=m_c-a_c$ とおき、$V$ を $V\bigoplus_c(1\otimes u_c)$ で置き換える。$1\otimes u_c$ は $\rho_c(g)\otimes1$ と可換であるから、これはどちらの標準形も変えない。$\iota_c:\mathbb C^{a_c}\to\mathbb C^{m_c}$ と
$\kappa_c:\mathbb C^{b_c}\to\mathbb C^{m_c}$ をそれぞれ最初の $a_c$ 個の座標と最後の $b_c$ 個の座標の埋め込みとし、
$$
V_1=Q_1^*V\bigoplus_c(1\otimes\iota_c),\qquad V_2=Q_2^*V\bigoplus_c(1\otimes\kappa_c).
$$
とおく。$W_1=V\bigoplus_c(1\otimes\iota_c)$ は値域 $\operatorname{range}P_1=\operatorname{range}Q_1$ をもつ等長作用素であるから
$V_1^*V_1=W_1^*P_1W_1=1$ であり、$\sum_ck_ca_c=\operatorname{rank}P_1=n_1$ であるから $V_1$ はユニタリである。さらに
$\pi_1(g)V_1=Q_1^*\pi(g)W_1=Q_1^*V\big(\bigoplus_c\rho_c(g)\otimes1_{m_c}\big)\bigoplus_c(1\otimes\iota_c)
=V_1\bigoplus_c\rho_c(g)\otimes1_{a_c}$ である。よって互いに非同値な既約表現 $\rho_c$ により $\pi_1\cong\bigoplus_c\rho_c\otimes1_{a_c}$ であり、重複度の一意性（2.4節）により $a_c=m_{\rho_c}(\pi_1)$ である。$V_2$ と $b_c$ についても同様である。最後に、$X=V\big(\bigoplus_c1\otimes x_c\big)V^*\in\pi(G)'$ に対して
$P_2XP_1=V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*$ である。ここで $y_c=\kappa_c^*x_c\iota_c$ は $x_c$ の左下の
$b_c\times a_c$ 隅であり、任意の行列である。補題 12.3 により、
$\operatorname{Hom}_G(\pi_1,\pi_2)=Q_2^*P_2\,\pi(G)'P_1Q_1$ は行列
$Q_2^*V\big(\bigoplus_c1\otimes\kappa_cy_c\iota_c^*\big)V^*Q_1=V_2\big(\bigoplus_c1\otimes y_c\big)V_1^*$ からなり、$V_1,V_2$ はユニタリであるから $(y_c)\mapsto T$ は単射である。$\pi_1\oplus\pi_2$ に現れない既約表現 $\rho$ については $m_\rho(\pi_1)m_\rho(\pi_2)=0$ である。$\square$

$n=4$ の系 12.2(b) の場合：$\pi_1\oplus\pi_2$ の四つのセクター（自明表現、次元 $2$ の表現、次元 $3$ の二つの表現）において、頂点特徴量 $\pi_1$ は $a=(1,0,0,1)$、辺特徴量 $\pi_2$ は
$b=(2,1,1,3)$ をもつから、$\dim\operatorname{Hom}_G=1\cdot2+1\cdot3=5=B_3$ である。

> **物理的意味.** 定理 12.4 は **選択則** の線形代数的な核心であり、ウィグナー・エッカートの定理のスカラーの場合である：対称性と可換な作用素は同じ対称性の型の状態しか結びつけることができず、一つの型の中では $1_{k_c}\otimes y_c$ として、すなわち多重項の中での位置を区別しない「換算行列要素」$y_c$ によって作用する。（ベクトルのように、それ自身が非自明な表現に従って変換する作用素には、さらにクレブシュ–ゴルダン係数が必要である。）対称なハミルトニアンは $\pi_1=\pi_2$ の場合であり、その縮退度 $k_c$ は既約表現の次元である。

## 12.2 ライブラリによる同変層

定理 2.6 の分解がわかれば、定理 12.4 の証明はそのままアルゴリズムである。`dc.from_group` は 4.3節の乱択法によってその分解を計算する。$\pi_1=\pi_2$ のときは定理 2.6 そのものである：重みはブロック
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
corners = node_to_edge.D.basis('commutant')[:, n:, :n].reshape(-1, 16 * n)   # Lemma 12.3: Q_2^* pi(G)' Q_1
assert torch.linalg.matrix_rank(corners).item() == bell[3]
```

層を同変にする第二の方法は、群上で平均することである。**ツワール** は次で定義される。
$$
\mathcal P(W)=\frac1{|G|}\sum_{g\in G}\pi_2(g)\,W\,\pi_1(g)^*,\qquad W\in M_{n_2\times n_1}(\mathbb C).
$$

**定理 12.5（ツワールは直交射影である）.** $\mathcal P$ は $M_{n_2\times n_1}(\mathbb C)$ から $\operatorname{Hom}_G(\pi_1,\pi_2)$ へのヒルベルト・シュミット内積に関する直交射影であり、
$$
\dim\operatorname{Hom}_G(\pi_1,\pi_2)=\frac1{|G|}\sum_{g\in G}\overline{\chi_{\pi_1}(g)}\,\chi_{\pi_2}(g).
$$
である。

*証明.* $h\in G$ に対し $\pi_2(h)\mathcal P(W)\pi_1(h)^*=\frac1{|G|}\sum_g\pi_2(hg)W\pi_1(hg)^*=\mathcal P(W)$ である。
$g\mapsto hg$ は $G$ の全単射だからである。よって $\mathcal P(W)$ は同変である。$W$ が同変ならばすべての項が $W$ に等しいので $\mathcal P(W)=W$ である。したがって $\mathcal P^2=\mathcal P$ であり、その値域は $\operatorname{Hom}_G(\pi_1,\pi_2)$ である。$\mathcal P$ は自己共役である：$\langle S,\mathcal P(W)\rangle=\frac1{|G|}\sum_g\operatorname{Tr}\big((\pi_2(g)^*S\pi_1(g))^*W\big)$
$=\langle\mathcal P(S),W\rangle$。ここで $g\mapsto g^{-1}$ と置換し、$\pi_i(g^{-1})=\pi_i(g)^*$ を用いた。自己共役な冪等作用素はその値域への直交射影である。そのトレースは階数に等しい。$M_{n_2\times n_1}$ 上の線形写像
$W\mapsto AWB$ のトレースは $\operatorname{Tr}A\operatorname{Tr}B$ である（$AE_{ij}B$ における行列単位 $E_{ij}$ の係数は $A_{ii}B_{jj}$ である）。$A=\pi_2(g)$、$B=\pi_1(g)^*$ とすれば公式が得られる。$\square$

定理 12.4 と合わせると $\sum_\rho m_\rho(\pi_1)m_\rho(\pi_2)=|G|^{-1}\sum_g\overline{\chi_{\pi_1}(g)}\chi_{\pi_2}(g)$ である。
$\pi_1=\pi_2$ のとき、これは 4.2節の $\dim\pi(G)'$ に対する指標公式であり、`dc.commutant_dimension` はこれを評価する。パラメータ化 $W\mapsto\mathcal P(W)$ はすべての同変写像に到達し、$\mathcal P(W)$ は $W$ に最も近い同変写像であり、$\mathcal P$ は自己共役であるから、$W$ に関する勾配は $T=\mathcal P(W)$ に関する勾配の $\mathcal P$ による像である。しかし
$n_1n_2-\dim\operatorname{Hom}_G$ 個のパラメータが核に浪費され（辺特徴量では $256-15$ 個）、密な $\pi(g)$ に対して一回の評価に $O(|G|n_1n_2(n_1+n_2))$ の計算量がかかり、$S_n$ では $|G|=n!$ である。ライブラリは同じ射影を分解から `D.conditional_expectation(W, 'commutant')` として、四回の $n\times n$ 行列積で計算する。

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

以上の層は線形である。ネットワークはそれらを非線形性と交互に重ねるので、非線形性もまた同変でなければならない。置換表現に対しては、各座標に関数を適用する操作は同変である。$(\pi(g)x)_{gy}=x_y$ だからである。他の表現ではそうではなく、ノルムに作用する非線形性や、$G$ が置換する基底の係数に作用する非線形性が用いられる（演習問題 3）。

> **物理的意味.** ツワールは可換子環への条件付き期待値である（第2章の演習問題 3）：結合をすべての対称操作について平均すると、ちょうどその対称な部分が残る。ハミルトニアンや相互作用を対称群にわたる和によって「対称化」することは、同じ射影である。

## 12.3 制約付きの層とスペクトル層

`EquivariantLinear` の重み $w=(w_c)_c$ は $*$-環の中にあるので、第1章と第11章のすべての構成がそれらに適用できる。$w\in\bigoplus_cM_{m_c}$ に対して $\Phi(w)=V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*$ と書く。

**命題 12.6.** (a) $\Phi$ は $\bigoplus_cM_{m_c}$ から $\pi(G)'$ の上への $*$-同型であり、$\Phi(w)$ が自己共役、ユニタリ、あるいは正であるための必要十分条件は、すべての $w_c$ がそうであることである。(b) 正規な $w$ と
$\bigcup_c\operatorname{spec}w_c$ 上の関数 $f$ に対し $f(\Phi(w))=\Phi\big((f(w_c))_c\big)$ である。特に、すべての正規な $X\in\pi(G)'$ に対し $f(X)\in\pi(G)'$ である。(c) 任意に固定したユニタリ $u_0\in\bigoplus_cM_{m_c}$ に対し、同変なユニタリは、$s$ が
$\bigoplus_cM_{m_c}$ 全体を動くときの $\Phi(u_0e^{s-s^*})$ にほかならない。特に、それらは歪自己共役な $S\in\pi(G)'$ の指数関数 $e^S$ にほかならない。

*証明.* (a) $(1\otimes x)(1\otimes y)=1\otimes xy$ かつ $(1\otimes x)^*=1\otimes x^*$ であるから、$\Phi$ は線形、単位的、乗法的で $*$ を保つ。また単射であり、定理 2.6 により $\pi(G)'$ の上への写像である。よって $\Phi(w)^*=\Phi(w)$ と
$w^*=w$ は同値であり、$\Phi(w)^*\Phi(w)=1=\Phi(w)\Phi(w)^*$ と $w^*w=1=ww^*$ は同値である。すべての $w_c=t_c^*t_c$ ならば
$\Phi(w)=\Phi(t)^*\Phi(t)\ge0$ である（系 1.7）。逆に $\Phi(w)\ge0$ ならば、$v\in\mathbb C^{m_c}$ とセクター $c$ のベクトル $\xi=V(e_1\otimes v)$ に対して $0\le\langle\xi,\Phi(w)\xi\rangle=\langle v,w_cv\rangle$ である。
(b) $w_c=u_c\Lambda_cu_c^*$ とする（定理 1.3）。このときユニタリ $U=V\bigoplus_c(1\otimes u_c)$ により $\Phi(w)=U\big(\bigoplus_c1\otimes\Lambda_c\big)U^*$ であり、これは $w_c$ の固有値がそれぞれ $k_c$ 回繰り返される対角化である。定義と定理 1.5(a) により $f(\Phi(w))=U\big(\bigoplus_c1\otimes f(\Lambda_c)\big)U^*=\Phi\big((f(w_c))_c\big)$ である。
(c) (a) により、同変なユニタリは $\bigoplus_cM_{m_c}$ のユニタリ $u$ による $\Phi(u)$ であり、補題 11.1 によりそれらの $u$ は
$u_0e^{s-s^*}$ にほかならない。$u_0=1$ とすると (b) により $\Phi(e^{s-s^*})=e^{\Phi(s-s^*)}$ であり、$\Phi(s-s^*)$ は
$\pi(G)'$ の歪自己共役な元全体を動く。$\square$

したがって直交（ユニタリ）な同変層は、`D.commutant_algebra()` 上の `optimize.UnitaryParam` に $\Phi$ を続けたものである。(c) により、このパラメータ化はすべての同変ユニタリに到達し、勾配法の各反復は厳密にユニタリかつ厳密に同変である。同様に `PositiveParam` は正値な同変層を与え、それは第11章と同様に正規化される（スペクトルは $[\ell,1]$ にあり最大値は $1$。(a) により、一般の正値層は $\Phi(t^*t)$ である）。また、拡散 $e^{-tL}$ のような、自己共役な同変生成子 $L$ のスペクトル層 $f(L)$ に必要な固有値分解の計算量は、$n^3=4096$ ではなく $\sum_cm_c^3$（ここでは $37$）であり、勾配はダレツキー・クレインの公式（`apply_function`）で与えられる。次のコードは同変なユニタリを $32$ 組の入出力の組にフィットさせ、(b) と、その $t$ に関する導関数 $\frac d{dt}e^{-tL}=-Le^{-tL}$ を密行列と照合する。

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

> **物理的意味.** 対称なハミルトニアン $H\in\pi(G)'$ は対称な時間発展
> $e^{-itH}\in\pi(G)'$ を生成する（第8章）。それはセクターを決して混ぜず、重複度空間にのみ作用する。同変なユニタリ層はそのような時間発展であり、与えられるのではなく学習されたものである。

## 12.4 不変量

ベクトル $x$ がすべての $g$ に対して $\pi(g)x=x$ を満たすとき **不変** であるといい、$\mathbb C^n$ 上の関数 $q$ が
$q(\pi(g)x)=q(x)$ を満たすとき不変であるという。不変な特徴量は、不変なネットワーク（たとえばグラフ分類器）が出力できるものである。

**命題 12.7.** (a) 不変ベクトルは直交射影
$P_0=|G|^{-1}\sum_g\pi(g)$ の値域であり、その次元は $m_{\mathrm{triv}}(\pi)=|G|^{-1}\sum_g\chi_\pi(g)$ である。また自明表現が $\rho_{c_0}$ として現れるならば、$P_0$ は $\pi(G)''$ の極小中心射影 $P_{c_0}$ である。
(b) エルミート形式 $q_A(x)=\langle x,Ax\rangle$（$A=A^*$）が不変であるための必要十分条件は $A\in\pi(G)'$ であり、不変なエルミート形式は次元 $\sum_cm_c^2$ の実ベクトル空間をなす。
(c) 座標 $x_1,\dots,x_n$ に関する $2$ 次の斉次不変多項式は、次元
$$
\frac1{|G|}\sum_{g\in G}\frac{\chi_\pi(g)^2+\chi_\pi(g^2)}2 .
$$
の空間をなす。置換表現に対しては、これは $\Omega$ の二つの元からなる多重集合 $\{x,y\}$ 上の $G$ の軌道の個数である。

*証明.* (a) 不変ベクトルとは $\operatorname{Hom}_G(\mathrm{triv},\pi)$ のことである。ここで $\mathrm{triv}$ は
$\mathbb C=M_{1\times1}$ 上の自明表現である。これに対してツワールは $\mathcal P(x)=P_0x$ であるから、定理 12.5 から射影と次元が、定理 12.4 から $m_{\mathrm{triv}}(\pi)$ が得られる。定理 2.6 の基底で、不変ベクトルを $\xi=\sum_c\xi_c$、
$\xi_c=\sum_jv_{c,j}\otimes e_j\in\mathbb C^{k_c}\otimes\mathbb C^{m_c}$ と書く。$\pi(g)$ はセクターを保つので、すべての $g$ に対し $\rho_c(g)v_{c,j}=v_{c,j}$ である。零でない $v_{c,j}$ は自明表現を担う不変な直線を張るから、そのとき $\rho_c$ は既約であることから自明であり、$c=c_0$ である。よって不変ベクトルはセクター $c_0$ のベクトル、すなわち $\pi(g)$ が $1$ として作用するところのベクトルであり、$P_{c_0}$ の値域をなす。したがって $P_{c_0}$ は $P_0$ である。
(b) $q_A(\pi(g)x)=\langle x,\pi(g)^*A\pi(g)x\rangle$ であり、すべての $x$ に対し $\langle x,Dx\rangle=0$ を満たす行列 $D$ は零である（偏極、系 1.4 の証明）。よって $q_A$ が不変であることは、すべての $g$ に対し $\pi(g)^*A\pi(g)=A$、すなわち $A\in\pi(G)'$ であることと同値である。命題 12.6(a) により、$\pi(G)'$ の自己共役元はエルミートな $w_c$ による $\Phi(w)$ であり、それらは次元 $\sum_cm_c^2$ の実空間をなす。
(c) 二次多項式 $p(x)=x^TBx$ はその対称行列 $B$ で定まる。$(g\cdot p)(x)=p(\pi(g)^{-1}x)$ とおくと
$g\cdot(h\cdot p)=(gh)\cdot p$ であり、$g\cdot p$ の行列は $\sigma(g)(B)=A^TBA$（$A=\pi(g)^{-1}$）である。平均
$M=|G|^{-1}\sum_g\sigma(g)$ は $\sigma(h)M=M$ を満たし、不変な $B$ に対し $MB=B$ を満たすから、（定理 12.5 の証明と同様に）値域が不変多項式全体である冪等写像であり、不変多項式の空間の次元は
$\operatorname{Tr}M=|G|^{-1}\sum_g\operatorname{Tr}\sigma(g)$ である。$A$ はユニタリであるから、
$\Lambda=\operatorname{diag}(\lambda)$ により $A=R\Lambda R^{-1}$ と書ける（定理 1.3）。対称行列
$B_{ij}=(R^{-1})^T(e_ie_j^T+e_je_i^T)R^{-1}$（$i\le j$）は基底をなし、
$A^TB_{ij}A=(R^{-1})^T\Lambda(e_ie_j^T+e_je_i^T)\Lambda R^{-1}=\lambda_i\lambda_jB_{ij}$ である。よって
$\operatorname{Tr}\sigma(g)=\sum_{i\le j}\lambda_i\lambda_j=\frac12\big((\operatorname{Tr}A)^2+\operatorname{Tr}A^2\big)
=\frac12\big(\chi_\pi(g^{-1})^2+\chi_\pi(g^{-2})\big)$ であり、$g$ の代わりに $g^{-1}$ について和をとれば公式が得られる。置換表現に対しては $G$ は単項式 $x_ax_b$ を置換するので、多項式が不変であるための必要十分条件は、その係数が単項式の軌道上で一定であることである。$\square$

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

$\mathbb C^4$ 上では、それらは $\sum_ix_i^2$ と $\sum_{i<j}x_ix_j$ である。すべての次数はモリーン（Molien）の公式
$\sum_d\dim\mathbb C[x]^G_dt^d=|G|^{-1}\sum_g\det(1-t\pi(g))^{-1}$ によって数えられ、その $t^2$ の係数が (c) である（証明なしに述べる。B. Sturmfels, *Algorithms in Invariant Theory*, Springer 2008, 定理 2.2.1）。

> **物理的意味.** ランダウの相転移の理論では、自由エネルギーは実の秩序変数 $x$ の不変多項式であり、$G$ は $x$ に実直交行列として作用する。その二次の部分は $x^TBx$ であり、$B$ は（$\pi(g)^T=\pi(g)^{-1}$ であるから）可換子環に属するので、命題 12.6(a) によりセクターについてブロック対角である：対称相は、$B$ の固有値が最初に符号を変えるセクターで不安定になり、（一般的には）そのセクターの既約表現が秩序相の対称性を決定する。

## 12.5 比較の対象

特定の群に対する同変層は成熟したライブラリによって提供されている。**escnn**（Cesa, Lang and Weiler 2022、e2cnn の後継）は、既約表現と核の調和基底から $O(2)$ と $O(3)$ の部分群に対するステアラブル畳み込みを構成し、**e3nn**（Geiger and Smidt 2022）は球面調和関数とクレブシュ–ゴルダンのテンソル積から $E(3)$ 同変なネットワークを構成する。**EMLP**（Finzi, Welling and Wilson 2021）は任意の行列群を扱う：離散的な生成元に対して $\pi_2(g)T=T\pi_1(g)$ を、リー環の生成元に対して $d\pi_2(X)T=Td\pi_1(X)$ を課し、解の基底を数値的な零空間として計算する——これは `dc.commutant_basis` の線形問題である。群畳み込み（Cohen and Welling 2016）は正則表現の場合であり、Kondor and Trivedi (2018) は、コンパクト群の等質空間上の特徴マップに対して、同変線形層が一般化された畳み込みであることを示している。置換については、DeepSets（Zaheer et al. 2017）と Maron et al. (2019) の基底が系 12.2 である。

本ライブラリは、格子上の畳み込み、連続空間上のステアラブル核、$\mathbb R$ 上の分解（実型、複素型、四元数型。演習問題 2 は実行列を扱う）、あるいは大きな $n$ を提供しない：$V$ は密な $n\times n$ 行列であり、すべての元で与えられた群に対して分解の計算量はランダムサンプル一つあたり $O(|G|n^3)$ である（リー環の生成元は `dc.from_generators` を通る。4.3節）。本ライブラリが付け加えるのは解空間の **構造** である：可換子環の基底は行列のリストにすぎないが、分解は可換子環を $*$-環 $\bigoplus_cM_{m_c}$ と同一視する。こうして積、随伴、ユニタリ性、正値性、関数計算が小さな行列の上でブロックごとに勾配付きで作用し（命題 12.6）、`to_commutant`/`from_commutant` が二つの描像の間を翻訳する。

## 演習問題

1. （**チャネル付きの DeepSets.**）$S_n$（$n\ge2$）は $\mathbb C^n\otimes\mathbb C^d$ に $P_\sigma\otimes1_d$ によって作用する。$\mathbb C^n\otimes\mathbb C^{d_1}$ から $\mathbb C^n\otimes\mathbb C^{d_2}$ への同変写像は $A,B\in M_{d_2\times d_1}$ による $1\otimes A+J\otimes B$ にほかならないことを示し、それを `EquivariantMap` で確かめよ。
2. （**実の層.**）$\pi_1,\pi_2$ が実行列をもつとする（たとえば置換表現）。実の同変行列は次元 $\dim_{\mathbb C}\operatorname{Hom}_G(\pi_1,\pi_2)$ の実ベクトル空間をなし、$T\mapsto\operatorname{Re}T$ は $\operatorname{Hom}_G(\pi_1,\pi_2)$ をその上へ写すことを示せ。
3. （**非線形性.**）(a) 関数 $f:\mathbb C\to\mathbb C$ を各座標に適用する操作は、すべての置換表現と可換であることを示せ。(b) $f(t)=\max(t,0)$ に対して、それが $\mathbb R^2$ の $\pi/2$ 回転と可換でないことを示せ。(c) $x\mapsto\varphi(\|x\|)x$ はすべてのユニタリ表現と可換であることを示せ。
4. （**層を決定するのに何個のサンプルが必要か？**）$T=\Phi(w)\in\pi(G)'$、$x_1,\dots,x_N\in\mathbb C^n$ とする。セクター $c$ における $V^*x_i$ を $X_{i,c}\in M_{k_c\times m_c}$ により $\sum_{a,j}(X_{i,c})_{aj}e_a\otimes e_j$ と書き、$S_c\in M_{Nk_c\times m_c}$ を行列 $X_{1,c},\dots,X_{N,c}$ を縦に積み重ねたものとする。$Tx_1=\dots=Tx_N=0$ を満たす同変な $T$ は次元 $\sum_cm_c(m_c-\operatorname{rank}S_c)$ の空間をなすことを示せ。ガウス分布に従う $x_i$ に対しては、これはほとんど確実に $\sum_cm_c\max(0,m_c-Nk_c)$ である。辺の表現では二つのサンプルで同変層が決定されるが、任意の $16\times16$ 行列には十六個のサンプルが必要である。

<details>
<summary>解答</summary>

1. $P_\sigma J=J=JP_\sigma$ であるから、これらの写像は同変である。それらは線形独立である：
   $1\otimes A+J\otimes B=0$ ならば、非対角ブロックから $B=0$、次いで $A=0$ が得られる。系 12.2(a) と定理 2.6 により $P_\sigma\otimes1_d\cong(\mathrm{triv}\otimes1_d)\oplus(\mathrm{std}\otimes1_d)$ であるから、二つの既約表現はともに重複度 $d$ をもち、定理 12.4 により $\dim\operatorname{Hom}_G=2d_1d_2$ となる。これは $1\otimes A+J\otimes B$ 全体のなす空間の次元である。

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

2. $T$ が同変ならば、$\overline{\pi_i(g)}=\pi_i(g)$ であるから $\bar T$ も同変である。よって
   $\operatorname{Re}T=(T+\bar T)/2$ と $\operatorname{Im}T=(T-\bar T)/2i$ は実かつ同変であり、したがって写像 $T\mapsto\operatorname{Re}T$ は実の同変行列全体の上への写像であり（その上では恒等写像である）、$H_{\mathbb R}$ を実の同変行列全体として $\operatorname{Hom}_G(\pi_1,\pi_2)=H_{\mathbb R}\oplus iH_{\mathbb R}$ である。すると $H_{\mathbb R}$ の実基底は $\operatorname{Hom}_G$ の複素基底である：それは全体を張り、実の $T_j$ に関する複素の線形関係 $\sum_j(\alpha_j+i\beta_j)T_j=0$ は実部と虚部に分かれる。補題 12.1 では軌道の指示行列がそのような実基底である。（定理 12.4 のブロック $y_c$ は実のデータに対しても複素である：$V$ は一般に複素であり、$\operatorname{Re}T$ によるパラメータ化は必要な数の二倍の実パラメータを用いる。）

3. (a) $(\pi(g)f(x))_{gy}=f(x)_y=f(x_y)=f((\pi(g)x)_{gy})$。(b) $R(x_1,x_2)=(-x_2,x_1)$ と $x=(1,-1)$ に対し、
   $f(Rx)=f(1,1)=(1,1)$ であるが $Rf(x)=R(1,0)=(0,1)$ である。(c) $\|\pi(g)x\|=\|x\|$ であるから
   $\varphi(\|\pi(g)x\|)\pi(g)x=\pi(g)\varphi(\|x\|)x$ である。

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

4. $(1\otimes w_c)(e_a\otimes e_j)=\sum_le_a\otimes(w_c)_{lj}e_l$ であるから、$V^*Tx_i$ はセクター $c$ において係数行列
   $X_{i,c}w_c^T$ をもつ。よってすべての $i$ に対して $Tx_i=0$ であることは、すべての $c$ に対して $S_cw_c^T=0$ であること、すなわち $w_c^T$ のすべての列が次元 $m_c-\operatorname{rank}S_c$ の空間 $\ker S_c$ に属することと同値である：この性質をもつ $w_c$ は次元 $m_c(m_c-\operatorname{rank}S_c)$ の空間をなし、$\Phi$ は単射である。ガウス分布に従う $x_i$ に対しては $V^*x_i$ もガウス分布に従い（$V$ はユニタリ）、したがって $S_c$ は独立なガウス成分をもち、ほとんど確実に階数 $\min(Nk_c,m_c)$ をもつ：その大きさの左上の正方小行列について、$|\det|^2$ は成分の実部と虚部の実多項式であって恒等的に零ではなく（小行列が単位行列のとき $1$ である）、そのような多項式はルベーグ零集合上でしか消えない（4.3節で証明されている）。辺の表現では $(k,m)=(1,2),(2,1),(3,1),(3,3)$ であり、寄与するのはセクター $(1,2)$ のみである：$N=1$ では $2\cdot(2-1)=2$、$N=2$ では $0$ である。

```python
basis = D_edge.basis('commutant')                                # (15, 16, 16), Hilbert-Schmidt orthonormal
for N, lost in ((1, 2), (2, 0)):
    Xs = torch.randn(16, N, dtype=dt)
    M = torch.stack([(B @ Xs).reshape(-1) for B in basis], dim=1)               # w -> (T x_1, ..., T x_N)
    assert torch.linalg.matrix_rank(M).item() == 15 - lost
```

</details>

**この先へ.** 本章が本講義の最後の章である。[第4章](04_symmetries.md) は、ここで用いた表現論、すなわちシューアの補題、指標、そして `dc.from_group` の背後にある乱択アルゴリズムを展開する（4.1–4.3節）。第9章（準備中）は同じセクターをスピン鎖に適用し、[第5章](05_random_matrices.md) の $\prod_cU(m_c)$ 上のハール測度は同変ユニタリ層の自然なランダム初期化を与える。数学については J.-P. Serre, *Linear Representations of Finite Groups* (Springer 1977) と W. Fulton and J. Harris, *Representation Theory* (Springer 1991) を、同変学習については M. Bronstein, J. Bruna, T. Cohen and P. Veličković, *Geometric Deep Learning* (2021) と M. Weiler et al., *Equivariant and Coordinate Independent Convolutional Networks* (2023) を、分解アルゴリズムについては `torch_vn_algebra/decompose.py` の参考文献を参照せよ。
