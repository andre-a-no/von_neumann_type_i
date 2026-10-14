[English](../03_traces_and_states.md) · [Español](../es/03_traces_and_states.md) · [Français](../fr/03_traces_and_states.md) · [Deutsch](../de/03_traces_and_states.md) · [中文](../zh/03_traces_and_states.md) · **日本語** · [Русский](../ru/03_traces_and_states.md)

# チュートリアル 3. トレース、状態、エントロピー

このチュートリアルは講義ノートの[第3章](../../lecture_notes/ja/03_traces_and_states.md)に対応する。同章は
$\mathcal M=\bigoplus_cM_{k_c}$ 上のトレースを分類し（[定理 3.1](../../lecture_notes/ja/03_traces_and_states.md#r-3-1)）、状態を密度で表し
（[定理 3.2](../../lecture_notes/ja/03_traces_and_states.md#r-3-2)）、基本的なエントロピー不等式を証明する。ここではライブラリの三つのトレースを用い、
それらの間で密度を変換し、$M_2\oplus M_3\oplus M_5$ 上で定理を数値的に確かめる。

<a id="s-3-1"></a>

## 3.1 トレース

*講義ノート：* [3.1節](../../lecture_notes/ja/03_traces_and_states.md#s-3-1)；[定理 3.1](../../lecture_notes/ja/03_traces_and_states.md#r-3-1) により、任意のトレースは $\mathrm{tr}_w(A)=\sum_cw_c\mathrm{Tr}\,A_c$ である。

ライブラリは三つの忠実なトレースを実装している。

| 名称 | 重み $w_c$ | $1$ での値 | 意味 |
|---|---|---|---|
| `Tr_blunt` | $1$ | $D$ | $\bigoplus_c\mathbb C^{k_c}$ のトレース（重複度 $m_c=1$） |
| `Tr_norm` | $1/k_c$ | $C$ | 各ブロックを正規化したもの |
| `tau_vN` | $1/(Ck_c)$ | $1$ | セクターの重みが等しく $1/C$ であるトレース状態 |

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, DensityMatrix, states, channels

torch.manual_seed(0)
alg = TypeIAlgebra([2, 3, 5], [2, 3, 5], precision='double')    # M_2 (+) M_3 (+) M_5, D = 10, C = 3
U = alg.random_unitary_operator(4)
A = alg.operator_from_eigenvalues(lambda d: torch.randn(4, d), batch_size=4, force_self_adjoint=True)
B = U @ A                                                       # a batch of four non-normal operators
one = alg.identity()

for f, value in ((lambda X: X.Tr_blunt(), 10), (lambda X: X.Tr_norm(), 3), (lambda X: X.tau_vN(), 1)):
    assert abs(f(one).item() - value) < 1e-12                   # normalisation
    assert (f(A @ B) - f(B @ A)).abs().max() < 1e-10            # trace property
assert (A.Tr_norm() - 3 * A.tau_vN()).abs().max() < 1e-12       # Tr_norm = C tau_vN
```

<a id="s-3-2"></a>

## 3.2 状態と密度

*講義ノート：* [3.2節](../../lecture_notes/ja/03_traces_and_states.md#s-3-2)；[定理 3.2](../../lecture_notes/ja/03_traces_and_states.md#r-3-2) は任意の忠実なトレースに対する状態の密度 $\rho_w$ を与え、$w_c(\rho_w)_c=w'_c(\rho_{w'})_c$ が成り立つ。

ライブラリは状態を
`Tr_blunt` に関する密度 $\rho=\rho_{\mathrm{blunt}}$ によって格納する（クラス `DensityMatrix`）。`density(trace)` と `DensityMatrix.from_density`
が変換を行う。状態の値は規約に依存しない。

```python
rho = states.random_density_matrix(alg, batch_size=4)          # Hilbert-Schmidt random states
for trace, tr in (('blunt', lambda X: X.Tr_blunt()), ('norm', lambda X: X.Tr_norm()),
                  ('tau_vN', lambda X: X.tau_vN())):
    d = alg.operator(rho.density(trace))
    assert (tr(d @ A) - (rho @ A).Tr_blunt()).abs().max() < 1e-10       # omega(A) for all traces
    assert (tr(d) - 1).abs().max() < 1e-10                               # tr_w(rho_w) = 1
    back = DensityMatrix.from_density(alg, rho.density(trace), trace=trace)
    assert (back.matrix - rho.matrix).abs().max() < 1e-12
```

<a id="s-3-3"></a>

## 3.3 状態の古典的部分

*講義ノート：* [3.3節](../../lecture_notes/ja/03_traces_and_states.md#s-3-3)：セクター確率 $p_c=\mathrm{Tr}\,\rho_c$、および $\mathcal M$ 上で恒等写像である写像 $\rho\mapsto\sum_cP_c\rho P_c$。

```python
p = rho.sector_probabilities()                                   # shape (batch, C)
assert torch.allclose(p.sum(-1), torch.ones(4, dtype=p.dtype))
P = [alg.central([1.0 if c == d else 0.0 for d in range(alg.C)]) for c in range(alg.C)]
measure_sector = channels.measurement_channel(P)
assert (measure_sector(rho).matrix - rho.matrix).abs().max() < 1e-12   # no disturbance
```

<a id="s-3-4"></a>

## 3.4 エントロピー

*講義ノート：* [3.4節](../../lecture_notes/ja/03_traces_and_states.md#s-3-4)；[定理 3.3](../../lecture_notes/ja/03_traces_and_states.md#r-3-3) は $S(\rho)=H(p)+\sum_cp_cS(\sigma_c)$ を与え、$\mathrm{tr}_w$ に関するエントロピーは $S_w(\rho)=S(\rho)+\sum_cp_c\log w_c$ である。

```python
S = rho.von_neumann_entropy()
sigma = [rho.matrix[:, c, :k, :k] / p[:, c, None, None] for c, k in enumerate(alg.k_factors)]
S_sec = [-(w * torch.log(w)).sum(-1) for w in (torch.linalg.eigvalsh(s) for s in sigma)]
H = -(p * torch.log(p)).sum(-1)
assert torch.allclose(S, H + sum(p[:, c] * S_sec[c] for c in range(alg.C)))    # Theorem 3.3

w_tau = torch.tensor([1 / (alg.C * k) for k in alg.k_factors], dtype=torch.float64)
d_tau = alg.operator(rho.density('tau_vN'), is_positive=True)
S_tau = -(d_tau @ d_tau.log()).tau_vN().real
assert torch.allclose(S_tau, S + (p * torch.log(w_tau)).sum(-1))
```

<a id="s-3-5"></a>

## 3.5 相対エントロピーとクラインの不等式

*講義ノート：* [3.5節](../../lecture_notes/ja/03_traces_and_states.md#s-3-5)；[系 3.6](../../lecture_notes/ja/03_traces_and_states.md#r-3-6)（$S(\rho)\le\log D$、等号は $\rho=1/D$ のとき）と[定理 3.7](../../lecture_notes/ja/03_traces_and_states.md#r-3-7)（中心の上への条件付き期待値 $E$ に対し $S(\rho\|\sigma)\ge S(E\rho\|E\sigma)=\sum_cp_c\log(p_c/q_c)$）。

```python
sigma_state = states.random_density_matrix(alg, batch_size=4)
q = sigma_state.sector_probabilities()
D_rel = states.relative_entropy(rho, sigma_state)
kl = (p * torch.log(p / q)).sum(-1)
E = channels.center_expectation(alg)
assert (D_rel >= kl - 1e-12).all()                                       # Theorem 3.7
assert torch.allclose(states.relative_entropy(E(rho), E(sigma_state)), kl)
assert (rho.von_neumann_entropy() <= math.log(10) + 1e-12).all()         # Corollary 3.6
assert abs(states.maximally_mixed_state(alg).von_neumann_entropy().item() - math.log(10)) < 1e-12
```

<a id="s-3-6"></a>

## 3.6 マイケルソン・コントラスト

*講義ノート：* [3.6節](../../lecture_notes/ja/03_traces_and_states.md#s-3-6)；[命題 3.8](../../lecture_notes/ja/03_traces_and_states.md#r-3-8) により、正かつ単位的な写像はマイケルソン・コントラストを増やさない。ライブラリはこれを `Operator.michelson_contrast` として計算する。

```python
A_pos = alg.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(4, d), batch_size=4, force_positive=True)
Phi = channels.random_mixed_unitary_channel(alg, 3, batch_size=4)       # unital
assert (Phi(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
assert (E(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
```

<a id="s-3-7"></a>

## 3.7 もっともらしいが偽である不等式

*講義ノート：* [3.7節](../../lecture_notes/ja/03_traces_and_states.md#s-3-7)；[命題 3.9](../../lecture_notes/ja/03_traces_and_states.md#r-3-9) により、$M_k$ の $X,Y\ge0$ に対し $\max_U|\mathrm{Tr}(XUY)|=\|YX\|_1\ge\mathrm{Tr}(XY)$ であり、等号は $XY=YX$ のときかつそのときに限り成り立つ。

ライブラリは極分解によって最大値を与える $U$ を求める：$YX=W|YX|$ に対しユニタリ $U=W^*$ が
$\|YX\|_1$ を達成する（[講義ノート、定理 1.9](../../lecture_notes/ja/01_matrix_algebras.md#r-1-9)）。ハール分布に従うユニタリ上のランダム探索
（ライブラリはこの方法で $U$ をサンプリングする）は平均的にはより悪い
（[講義ノート、3.7節](../../lecture_notes/ja/03_traces_and_states.md#s-3-7)）。そこで述べた確率の上界 $b^4/(a^4k^2)$ は、$b/a=6$ かつ $k=6$ である
以下の実行では何も言っていない。
コードは $k=6$ の可換でない一つの組について、厳密な上限と $20000$ 個のハール・サンプルの最良値を、ともに
$\mathrm{Tr}(XY)$ を引いて計算する。命題 3.9 により、どのサンプルも上限を超えることはできない。最良のサンプルが
上限よりどれだけ下にとどまるかは、この実行の*数値的観察*であり、定理ではない。

```python
Mk = TypeIAlgebra([6], [6], precision='double')
X = Mk.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
Y = Mk.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
assert (X @ Y - Y @ X).matrix.abs().max().item() > 1e-6             # X and Y do not commute
YX = Y @ X
W = YX @ YX.abs().inverse()                                          # polar part of YX
gap = (X @ W.adjoint() @ Y).trace.abs() - (X @ Y).trace.real
assert abs(gap.item() - (YX.trace_norm() - (X @ Y).trace.real).item()) < 1e-10   # Theorem 1.9
assert gap.item() > 0                                  # Proposition 3.9: the inequality fails
samples = Mk.random_unitary_operator(20_000)
best = (X @ samples @ Y).trace.abs().max() - (X @ Y).trace.real
print(f"sup - Tr(XY) = {gap.item():.4f}; numerical observation: best of 20000 Haar samples "
      f"- Tr(XY) = {best.item():.4f}")
assert best.item() <= gap.item() + 1e-10               # Proposition 3.9: no sample beats the sup
```

<a id="s-3-8"></a>

## 3.8 比較の対象

QuTiP（`entropy_vn`, `entropy_relative`）やその他の量子情報のパッケージ（例えば
`toqito`）は、単一のヒルベルト空間上の密度行列のエントロピーと相対エントロピーを、その空間のトレースを用いて計算する。$\bigoplus_c\mathbb C^{k_c}$ 上ではこれは本章の blunt トレースである（対数の底はパッケージによって異なりうるが、その違いを除いて）。それらが提供しないのは定理 3.1 のトレースの族と定理 3.2 の変換であり、これらは環が重複度をもって作用するとき、あるいはセクターに異なる重みをつけなければならないときに直ちに必要となる。

<a id="x-3"></a>

## 演習問題

講義ノートの[演習問題 1](../../lecture_notes/ja/03_traces_and_states.md#x-3)：トレース状態 $\tau_{\mathrm{vN}}$
（`states.tracial_state`）の $\tau_{\mathrm{vN}}$ に関する密度は $1$ であり、そのエントロピーは
$S=\log C+\frac1C\sum_c\log k_c$ である。$M_2\oplus M_3\oplus M_5$ について両方を確かめよ。

<details>
<summary>解答</summary>

```python
tau = states.tracial_state(alg)
assert (tau.density('tau_vN') - one.matrix).abs().max() < 1e-12
assert abs(tau.von_neumann_entropy().item() - (math.log(3) + sum(math.log(k) for k in (2, 3, 5)) / 3)) < 1e-12
```

</details>

**次へ：** [チュートリアル 4. 対称性と表現](04_symmetries.md)
