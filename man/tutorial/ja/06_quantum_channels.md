[English](../06_quantum_channels.md) · [Español](../es/06_quantum_channels.md) · [Français](../fr/06_quantum_channels.md) · [Deutsch](../de/06_quantum_channels.md) · [中文](../zh/06_quantum_channels.md) · **日本語** · [Русский](../ru/06_quantum_channels.md)

# チュートリアル 6. 量子チャネル

このチュートリアルは講義ノートの[第6章](../../lecture_notes/ja/06_quantum_channels.md)に対応する。同章では行列環の直和の間の完全正値写像の諸記述
— チョイ行列、クラウス作用素、スタインスプリング拡張 — を証明し、双対、不動点、遷移行列、ランダムなチャネルを調べる。
ここではライブラリがチャネルをどのように格納するかを見て、これらの結果を数値的に確かめる。

本章を通じて、$`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$ および
$`\mathcal N=\bigoplus_{d=1}^{C'}M_{k'_d}(\mathbb C)`$ とし、それぞれの極小中心射影を $`P_c`$ と
$`Q_d`$ とする。$`A=\bigoplus_cA_c`$、$`M_{k'\times k}`$ は複素 $`k'\times k`$ 行列の空間、
$`E_{ij}=|i\rangle\langle j|`$ は行列単位、$`\mathrm{Tr}`$ は blunt トレース
$`\sum_c\mathrm{Tr}\,A_c`$、$`\mathrm{tr}_w=\sum_cw_c\mathrm{Tr}\,A_c`$ は忠実な
トレース（定理 3.1）である。ブロック行列は $`\mathbb C^n\otimes\mathbb C^k`$ に作用し、ブロックの添字を先に書く。

<a id="s-6-1"></a>

## 6.1 完全正値写像とチョイの定理

*講義ノート：* [6.1節](../../lecture_notes/ja/06_quantum_channels.md#s-6-1)：[命題 6.1](../../lecture_notes/ja/06_quantum_channels.md#r-6-1)（転置）、[定理 6.2](../../lecture_notes/ja/06_quantum_channels.md#r-6-2)（チョイ）、[系 6.3](../../lecture_notes/ja/06_quantum_channels.md#r-6-3)（クラウス階数）、[命題 6.4](../../lecture_notes/ja/06_quantum_channels.md#r-6-4)（直和の間の写像）。

ライブラリには二つのクラスがある。`Channel` は **セクター保存** 写像（$`C'=C`$、$`k'=k`$、
$`d\ne c`$ に対し $`\Phi_{dc}=0`$）を格納する。そのクラウス作用素 $`K_l=\bigoplus_cK^{cc}_l`$ は $`\mathcal M`$ に属し、
形状 `(batch, r, C, k_max, k_max)` のテンソルとして格納される。`InterSectorChannel` は
[命題 6.4](../../lecture_notes/ja/06_quantum_channels.md#r-6-4) の一般の写像を形状 `(batch, r, C', C, k'_max, k_max)` で格納する。`Channel.choi()` は各セクターについて
$`J(\Phi_{cc})`$ を返し、`Channel.superoperator()` は行優先のベクトル
$`\mathrm{vec}(X)_{(a,b)}=X_{ab}`$ に作用する $`\Phi_{cc}`$ の行列 $`S_c`$ を返す。どちらも同じ数
$`\Phi(E_{ij})_{ab}=J_{(i,a),(j,b)}=S_{(a,b),(i,j)}`$ を異なる順序で含んでいる。`Channel.from_superoperator`
は、チョイ行列がエルミートでないとき、すなわち $`\|J-J^*\|>2\cdot10^{-10}s`$（倍精度の場合）のとき `ValueError` を送出する。
ここで $`s=\max(1,\|J\|)`$ であり、ノルムはスペクトルノルムである。そうでなければ $`J`$ をそのエルミート部分 $`(J+J^*)/2`$ で置き換え、
[定理 6.2](../../lecture_notes/ja/06_quantum_channels.md#r-6-2) の証明のようにこれを対角化し、ある固有値が $`-10^{-10}s`$ より小さければ `ValueError` を送出し、
残りの負の固有値を $`0`$ にする。$`J^*=J`$ はすべての $`X`$ に対し $`\Phi(X^*)=\Phi(X)^*`$ となることと同値であるから
（[講義ノート、6.1節](../../lecture_notes/ja/06_quantum_channels.md#s-6-1)）、これらの丸め誤差の許容範囲を除いて、このメソッドが拒否するのは
随伴を保たない写像と、それ以外の写像のうち CP でないものにちょうど一致する。
例えば $`X\mapsto iX`$ と転置は拒否される。

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, states, channels
from torch_vn_algebra.channels import Channel, InterSectorChannel

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def unit(i, j, k):
    X = torch.zeros(k, k, dtype=dt)
    X[i, j] = 1
    return X

def choi(f, k):
    """Choi matrix sum_ij E_ij (x) f(E_ij) of a linear map f on k x k matrices."""
    return sum(torch.kron(unit(i, j, k), f(unit(i, j, k)).contiguous()) for i in range(k) for j in range(k))

def to_superoperator(J, k):
    """S[(a,b),(i,j)] = J[(i,a),(j,b)], with a leading batch axis."""
    return J.reshape(k, k, k, k).permute(1, 3, 0, 2).reshape(1, k * k, k * k)

ev = torch.linalg.eigvalsh(choi(lambda X: X.T, 2))
print(ev)                                                 # the Choi matrix of the transpose is the swap
assert abs(ev[0] + 1) < 1e-12                             # Proposition 6.1

alg = TypeIAlgebra([2, 3], [2, 3], precision='double')    # M_2 (+) M_3
Phi = channels.random_channel(alg, kraus_rank=2)          # Kraus tensor of shape (1, 2, 2, 3, 3)

def block(Phi, c):
    """The block Phi_cc of a sector-preserving map, as a function on k_c x k_c matrices."""
    ks = Phi.algebra.k_factors
    def f(X):
        blocks = [torch.zeros(k, k, dtype=dt) for k in ks]
        blocks[c] = X
        return Phi(Phi.algebra.from_blocks(blocks)).matrix[0, c, :ks[c], :ks[c]]
    return f

J, S = Phi.choi(), Phi.superoperator()                    # one matrix per sector, (1, k_c^2, k_c^2)
rho = states.random_density_matrix(alg)
out = Phi(rho)
for c, k in enumerate(alg.k_factors):
    assert close(J[c][0], choi(block(Phi, c), k))         # the convention of Theorem 6.2
    w = torch.linalg.eigvalsh(J[c][0])
    assert w.min() > -1e-12 and int((w > 1e-10).sum()) == 2      # J >= 0, Kraus rank 2 (Corollary 6.3)
    assert close(S[c][0], to_superoperator(J[c][0], k)[0])
    assert close(S[c][0] @ rho.matrix[0, c, :k, :k].reshape(-1), out.matrix[0, c, :k, :k].reshape(-1))
assert Phi.is_trace_preserving()
assert close(Channel.from_superoperator(alg, S)(rho), out)       # back from the superoperator
try:
    Channel.from_superoperator(alg, [to_superoperator(choi(lambda X: X.T, k), k) for k in (2, 3)])
    raise AssertionError("the transpose was accepted")
except ValueError as e:
    print(e)                                              # superoperator is not completely positive
try:
    Channel.from_superoperator(alg, [1j * torch.eye(k * k, dtype=dt)[None] for k in (2, 3)])
    raise AssertionError("X -> iX was accepted")
except ValueError as e:
    print(e)                                              # superoperator does not preserve Hermiticity
```

<a id="s-6-2"></a>

## 6.2 スタインスプリング拡張とクラウス作用素の自由度

*講義ノート：* [6.2節](../../lecture_notes/ja/06_quantum_channels.md#s-6-2)：[定理 6.5](../../lecture_notes/ja/06_quantum_channels.md#r-6-5)（スタインスプリング）と[定理 6.7](../../lecture_notes/ja/06_quantum_channels.md#r-6-7)（ユニタリ自由度）。6.1節のチャネルのブロック $`c=1`$ について、コードは $`V`$ が等長作用素であること、$`\Phi(\rho)=\mathrm{Tr}_E(V\rho V^*)`$ と $`\Phi^*(B)=V^*(1\otimes B)V`$ が成り立つこと、そしてクラウス作用素 $`\sum_lu_{jl}K_l`$ が同じチャネルを与えることを確かめる。

コードでは、$`V`$ の環境の添字はクラウス・テンソルの最初の添字である：

```python
K = Phi.kraus[0, :, 1, :3, :3]                            # the two Kraus operators of block c = 1
V = K.reshape(2 * 3, 3)                                   # V = sum_l e_l (x) K_l : C^3 -> C^2 (x) C^3
assert close(V.conj().T @ V, torch.eye(3, dtype=dt))      # TP <=> V is an isometry
r1 = rho.matrix[0, 1, :3, :3]
TrE = torch.einsum('lalb->ab', (V @ r1 @ V.conj().T).reshape(2, 3, 2, 3))
assert close(TrE, out.matrix[0, 1, :3, :3])               # Phi(rho) = Tr_E(V rho V^*)
B1 = torch.randn(3, 3, dtype=dt)
Bop = alg.from_blocks([torch.zeros(2, 2, dtype=dt), B1])
assert close(V.conj().T @ torch.kron(torch.eye(2, dtype=dt), B1) @ V,
             Phi.adjoint()(Bop).matrix[0, 1, :3, :3])    # Phi^*(B) = V^*(1 (x) B)V

u = alg.random_unitary(3)                                 # Theorem 6.7 with n = 3: pad with a zero operator
Kpad = torch.cat([Phi.kraus, torch.zeros_like(Phi.kraus[:, :1])], dim=1)
PhiL = Channel(alg, torch.einsum('jl,blcxy->bjcxy', u, Kpad))
assert close(PhiL(rho), out) and all(close(a, b) for a, b in zip(PhiL.choi(), J))
```

<a id="s-6-3"></a>

## 6.3 トレースに関する双対

*講義ノート：* [6.3節](../../lecture_notes/ja/06_quantum_channels.md#s-6-3)：$`(w,v)`$ に関する双対 $`\Phi^\sharp`$、[命題 6.8](../../lecture_notes/ja/06_quantum_channels.md#r-6-8)、[系 6.9](../../lecture_notes/ja/06_quantum_channels.md#r-6-9)。

ライブラリでは、`Channel.adjoint()` は $`\Phi^*`$ であり、`InterSectorChannel.adjoint(trace=...)` は
両方の環で同じ名前のトレース（`'blunt'`、`'norm'`、`'tau_vN'`）を用いたときの $`\Phi^\sharp`$ である。`tau_vN` の重み
$`1/(Ck_c)`$ と $`1/(C'k'_d)`$ はセクターの個数にもよる。ランダムな写像
$`M_2\oplus M_3\to M_1\oplus M_2\oplus M_2`$ と三つのトレースについて[命題 6.8](../../lecture_notes/ja/06_quantum_channels.md#r-6-8) を確かめ、$`M_1\oplus M_2`$ 上の粒子損失
$`\Phi(x)=(x_0+\mathrm{Tr}\,x_1)\oplus0`$ について[系 6.9(a)](../../lecture_notes/ja/06_quantum_channels.md#r-6-9) を確かめる。粒子損失は blunt トレースに関しては TP であるが、
`Tr_norm` に関してはそうでない（[講義ノート、6.3節](../../lecture_notes/ja/06_quantum_channels.md#s-6-3)）：

```python
alg_out = TypeIAlgebra([1, 2, 2], [1, 2, 2], precision='double')
Psi = channels.random_inter_sector_channel(alg, alg_out, kraus_rank=2)     # Kraus tensor (1, 2, 3, 2, 2, 3)
x = alg.operator(torch.randn(1, 2, 3, 3, dtype=dt))
y = alg_out.operator(torch.randn(1, 3, 2, 2, dtype=dt))
weights = {'blunt': lambda k, C: 1.0, 'norm': lambda k, C: 1 / k, 'tau_vN': lambda k, C: 1 / (C * k)}
trace = {'blunt': lambda X: X.Tr_blunt(), 'norm': lambda X: X.Tr_norm(), 'tau_vN': lambda X: X.tau_vN()}

def z(a, t):                                              # the central element sum_c w_c P_c of trace t
    return a.central(torch.tensor([weights[t](k, a.C) for k in a.k_factors], dtype=torch.float64))

for t in trace:
    dual = Psi.adjoint(trace=t)
    assert (trace[t](Psi(x) @ y) - trace[t](x @ dual(y))).abs().max() < 1e-10   # the definition
    assert close(dual(y), z(alg, t).inverse() @ Psi.adjoint()(z(alg_out, t) @ y))   # Proposition 6.8
dev = [(Psi.adjoint(trace=t)(alg_out.identity()) - alg.identity()).matrix.abs().max().item() for t in trace]
print(dev)                         # TP only for the blunt traces: the other duals are not unital
assert Psi.is_trace_preserving() and dev[0] < 1e-12 and min(dev[1:]) > 0.1

alg01 = TypeIAlgebra([1, 2], [1, 2], precision='double')
loss = InterSectorChannel.from_blocks(alg01, alg01, {(0, 0): [torch.ones(1, 1)],
                                                     (0, 1): [torch.tensor([[1.0, 0.0]]), torch.tensor([[0.0, 1.0]])]})
assert loss.is_trace_preserving()
assert close(loss.adjoint(trace='norm')(alg01.identity()), alg01.central([1.0, 2.0]))
s = states.random_density_matrix(alg01)
assert abs(loss(s).Tr_norm() - (s @ alg01.central([1.0, 2.0])).Tr_norm()) < 1e-12   # Corollary 6.9(a)
print(s.Tr_norm().real.item(), loss(s).Tr_norm().real.item())      # Tr_norm is not conserved
```

<a id="s-6-4"></a>

## 6.4 デコヒーレンス

*講義ノート：* [6.4節](../../lecture_notes/ja/06_quantum_channels.md#s-6-4)：[定理 6.10](../../lecture_notes/ja/06_quantum_channels.md#r-6-10)（単位的チャネルの不動点）、[定理 6.11](../../lecture_notes/ja/06_quantum_channels.md#r-6-11)（単位的チャネルはエントロピーを増加させる）、[命題 6.12](../../lecture_notes/ja/06_quantum_channels.md#r-6-12)（遷移行列）。

ライブラリの標準的なチャネルは各ブロックの内部で作用する（$`p,\gamma\in[0,1]`$）：

| チャネル | $`\Phi(\rho)_c`$ | ブロック $`c`$ のクラウス作用素 |
|---|---|---|
| `dephasing_channel` | $`(1-p)\rho_c+p\mathrm{diag}\,\rho_c`$ | $`\sqrt{1-p}\,1`$, $`\sqrt p\,E_{ii}`$ |
| `depolarizing_channel` | $`(1-p)\rho_c+p\mathrm{Tr}(\rho_c)\,1/k_c`$ | $`\sqrt{1-p}\,1`$, $`\sqrt{p/k_c}\,E_{ij}`$ |
| `center_expectation` | $`\mathrm{Tr}(\rho_c)\,1/k_c`$ | $`\sqrt{1/k_c}\,E_{ij}`$ |
| `amplitude_damping_channel` | 準位 $`i>0`$ は $`0`$ へ減衰する | $`E_{00}+\sqrt{1-\gamma}\sum_{i>0}E_{ii}`$, $`\sqrt\gamma\,E_{0i}`$ |
| `random_mixed_unitary_channel` | $`\sum_jp_jU_j\rho U_j^*`$ | $`\sqrt{p_j}\,U_j`$ |

$`M_2\oplus M_3`$ 上の表の各チャネルについて、コードはその不動点（$`S_c-1`$ の零空間）の次元と、クラウス作用素の
可換子環の次元を計算する：単位的チャネルでは両者は一致し（[定理 6.10](../../lecture_notes/ja/06_quantum_channels.md#r-6-10)）、振幅減衰では次元は等しいが異なる空間である。
次に $`100`$ 個のランダムな状態で[定理 6.11](../../lecture_notes/ja/06_quantum_channels.md#r-6-11) を確かめ、6.3節の写像 `Psi` と粒子損失について[命題 6.12](../../lecture_notes/ja/06_quantum_channels.md#r-6-12) を確かめる。

```python
def fixed_dims(Phi):
    """dim {X in M_{k_c} : Phi_cc(X) = X} for every sector: the null space of S_c - 1."""
    return [int((torch.linalg.svdvals(S_c[0] - torch.eye(k * k, dtype=dt)) < 1e-9).sum())
            for S_c, k in zip(Phi.superoperator(), Phi.algebra.k_factors)]

def commutant_dims(Phi):
    """dim of the commutant of {K_l, K_l^*} in every block (row-major vec, as in Chapter 2)."""
    out = []
    for c, k in enumerate(Phi.algebra.k_factors):
        I, Ks = torch.eye(k, dtype=dt), Phi.kraus[0, :, c, :k, :k]
        ops = [A.contiguous() for A in torch.cat([Ks, Ks.mH])]
        L = torch.cat([torch.kron(A, I) - torch.kron(I, A.T.contiguous()) for A in ops])
        out.append(k * k - int((torch.linalg.svdvals(L) > 1e-9).sum()))
    return out

zoo = {'dephasing': channels.dephasing_channel(alg, 0.3),
       'depolarizing': channels.depolarizing_channel(alg, 0.3),
       'centre': channels.center_expectation(alg),
       'mixed unitary': channels.random_mixed_unitary_channel(alg, 2),
       'amplitude damping': channels.amplitude_damping_channel(alg, 0.4)}
for name, ch in zoo.items():
    print(f"{name:18s} unital {ch.is_unital()!s:5s} fixed {fixed_dims(ch)} commutant {commutant_dims(ch)}")
    assert ch.is_trace_preserving()
    if ch.is_unital():
        assert fixed_dims(ch) == commutant_dims(ch)                         # Theorem 6.10
assert fixed_dims(zoo['dephasing']) == [2, 3] and fixed_dims(zoo['depolarizing']) == [1, 1]
E00 = alg.from_blocks([unit(0, 0, 2), unit(0, 0, 3)])
E01 = alg.from_blocks([unit(0, 1, 2), unit(0, 1, 3)])
assert close(zoo['amplitude damping'](E00), E00) and not close(E00 @ E01, E01 @ E00)   # equal dims, different spaces
assert fixed_dims(zoo['amplitude damping']) == commutant_dims(zoo['amplitude damping']) == [1, 1]

rhos = states.random_density_matrix(alg, batch_size=100)
for name in ('dephasing', 'depolarizing', 'centre', 'mixed unitary'):                 # Theorem 6.11
    assert (zoo[name](rhos).von_neumann_entropy() >= rhos.von_neumann_entropy() - 1e-12).all()
    assert (zoo[name](rhos).purity() <= rhos.purity() + 1e-12).all()
mm = states.maximally_mixed_state(alg)
S_after = channels.amplitude_damping_channel(alg, 1.0)(mm).von_neumann_entropy().item()
print(mm.von_neumann_entropy().item(), S_after)          # not unital: the entropy drops
assert abs(S_after + 0.4 * math.log(0.4) + 0.6 * math.log(0.6)) < 1e-12

T = Psi.transition_matrix()[0]                            # (C', C) = (3, 2)
p, q = rho.sector_probabilities()[0], Psi(rho).sector_probabilities()[0]
assert torch.allclose(T.sum(0), torch.ones(2, dtype=T.dtype))                       # Proposition 6.12
assert torch.allclose(q, Psi.transition_matrix(rho)[0] @ p)
classical = channels.center_expectation(alg_out) @ Psi @ channels.center_expectation(alg)
Tp = T @ p
assert close(classical(rho), alg_out.central([Tp[d] / k for d, k in enumerate(alg_out.k_factors)]))
assert torch.allclose(loss.transition_matrix()[0], torch.tensor([[1.0, 1.0], [0.0, 0.0]], dtype=T.dtype))
```

完全な振幅減衰のもとで、$`M_2\oplus M_3`$ の最大混合状態のエントロピーは
$`\log5\approx1.609`$ から $`H(\tfrac25,\tfrac35)\approx0.673`$ へ下がる：単位的でないチャネルは状態を純粋化しうる。

$`M_2`$ 上で位相緩和 $`D_p`$ と振幅減衰 $`A_\gamma`$ は可換であり、その合成はコヒーレンスに $`(1-p)\sqrt{1-\gamma}`$ を掛ける。
$`1-p=e^{-t/T_\varphi}`$、$`1-\gamma=e^{-t/T_1}`$ とすれば $`1/T_2=1/(2T_1)+1/T_\varphi`$ が得られる（[講義ノート、6.4節](../../lecture_notes/ja/06_quantum_channels.md#s-6-4)、
*量子ビットのコヒーレンス*）。コードはこれらの恒等式と合成則 $`D_{p(t)}D_{p(s)}=D_{p(t+s)}`$、$`A_{\gamma(t)}A_{\gamma(s)}=A_{\gamma(t+s)}`$ を確かめる：

```python
q2 = TypeIAlgebra([2], [2], precision='double')
p_, g_ = 0.3, 0.4
Dp, Ag = channels.dephasing_channel(q2, p_), channels.amplitude_damping_channel(q2, g_)
r2 = states.random_density_matrix(q2)
assert close(Dp(Ag(r2)), Ag(Dp(r2)))                                         # D_p A_gamma = A_gamma D_p
coh = Dp(Ag(r2)).matrix[0, 0, 0, 1] / r2.matrix[0, 0, 0, 1]
assert abs(coh - (1 - p_) * math.sqrt(1 - g_)) < 1e-12                       # coherence factor
T1, Tphi, t1, t2 = 2.0, 3.0, 0.7, 1.1
gam = lambda t: 1 - math.exp(-t / T1)
assert close(channels.amplitude_damping_channel(q2, gam(t1))(channels.amplitude_damping_channel(q2, gam(t2))(r2)),
             channels.amplitude_damping_channel(q2, gam(t1 + t2))(r2))       # A_gamma(t) A_gamma(s) = A_gamma(t+s)
pdeph = lambda t: 1 - math.exp(-t / Tphi)
assert close(channels.dephasing_channel(q2, pdeph(t1))(channels.dephasing_channel(q2, pdeph(t2))(r2)),
             channels.dephasing_channel(q2, pdeph(t1 + t2))(r2))           # D_p(t) D_p(s) = D_p(t+s)
both = channels.dephasing_channel(q2, pdeph(t1))(channels.amplitude_damping_channel(q2, gam(t1))(r2))
coh_t = (both.matrix[0, 0, 0, 1] / r2.matrix[0, 0, 0, 1]).item()
assert abs(coh_t - math.exp(-t1 * (1 / (2 * T1) + 1 / Tphi))) < 1e-12       # 1/T_2 = 1/(2T_1) + 1/T_phi
```

<a id="s-6-5"></a>

## 6.5 ランダムなチャネル

*講義ノート：* [6.5節](../../lecture_notes/ja/06_quantum_channels.md#s-6-5)：[補題 6.13](../../lecture_notes/ja/06_quantum_channels.md#r-6-13) と[定理 6.15](../../lecture_notes/ja/06_quantum_channels.md#r-6-15)（ランダムなチャネルの出力純度）。

**ハール・ランダムな等長作用素** $`V:\mathbb C^k\to\mathbb C^N`$ とは、ライブラリ（`random_isometry`）と同様に、
$`G=QR`$ の因子 $`Q`$ のことである。ここで $`G`$ は $`N\times k`$ のジニブル行列（成分は独立で、密度
$`\pi^{-1}e^{-|z|^2}`$ をもつ）、$`Q^*Q=1`$、$`R`$ は対角成分が正の上三角行列である（CPU 上では、ライブラリは LAPACK の分解 $`G=Q'R'`$ をとり、$`Q'`$ の第 $`i`$ 列に
$`R'_{ii}/|R'_{ii}|`$ を掛ける。定理 5.3 と同様に、これはすべての $`R'_{ii}\ne0`$ のとき、すなわち $`G`$ の階数が $`k`$ のとき
この分解を与え、補題 6.13 の証明によりこれはほとんど確実に成り立つ。GPU 上の大きなバッチに対してはグラム・シュミットの
直交化を用い、これは直接この分解を与える）。
`random_channel(alg, r)` はこのような $`V_c:\mathbb C^{k_c}\to\mathbb C^r\otimes\mathbb C^{k_c}`$ を
各ブロックで独立に生成し、定理 6.5 のクラウス作用素を用いる。

講義ノートの[定理 6.15](../../lecture_notes/ja/06_quantum_channels.md#r-6-15)により、サイズ $`k`$ のセクターの純粋状態 $`\psi`$ と、`random_channel` でクラウス階数 $`r`$ で
生成した $`\Phi`$ に対し、$`\mathbb E\mathrm{Tr}\,\Phi(\psi)^2=(k+r)/(kr+1)`$ かつ $`\mathbb E\,\Phi(\psi)=1_k/k`$ である。
両者を $`20\,000`$ 個のサンプルから、平均の標準誤差とともに推定する：

```python
B = 20_000
for r in (1, 2, 4):
    Phi_r = channels.random_channel(alg, kraus_rank=r, batch_size=B)
    for c, k in enumerate(alg.k_factors):
        psi = torch.randn(k, 1, dtype=dt)
        psi = psi / psi.norm()
        out_r = Phi_r(alg.from_blocks([psi @ psi.conj().T if d == c else torch.zeros(kk, kk, dtype=dt)
                                       for d, kk in enumerate(alg.k_factors)]))     # a batch of B outputs
        purity = (out_r @ out_r).trace.real
        mean, err = purity.mean().item(), purity.std().item() / math.sqrt(B)
        exact = (k + r) / (k * r + 1)
        print(f"k = {k}, r = {r}: E Tr Phi(psi)^2 = {mean:.4f} +- {err:.4f}   (exact {exact:.4f})")
        assert abs(mean - exact) < 5 * err + 1e-12
        avg = out_r.matrix.mean(0)[c, :k, :k]
        assert (avg - torch.eye(k, dtype=dt) / k).abs().max() < 0.02          # E Phi(psi) = 1/k
```

<a id="s-6-6"></a>

## 6.6 比較の対象

QuTiP（`to_choi`, `to_kraus`, `to_super`, `to_stinespring`）、Qiskit（`qiskit.quantum_info`：
クラス `Kraus`, `Choi`, `SuperOp`, `Stinespring` と `random_quantum_channel`。Qiskit Aer の
`amplitude_damping_error` などの雑音モデル）、toqito（`toqito.channel_ops`, `toqito.channel_props`）および
QETLAB（`ChoiMatrix`, `KrausOperators`, `ApplyMap`, `RandomSuperoperator`）は、6.1–6.2節の表現の間の変換、
完全正値性の判定、そして（その多くは）ランダムなチャネルの生成を、
全行列環の間の写像 $`M_n\to M_m`$ について行う。
それらの規約は互いに異なる（列方向か行方向の積み重ね、チョイ行列の因子の順序）ので、
数値を比べる前に規約を比べること。$`\mathcal M`$ 上のセクター保存チャネルは、それらのパッケージでは
$`M_D`$（$`D=\sum_ck_c`$）上のチャネルとして扱え、その超演算子のサイズはセクターごとの $`k_c^2\times k_c^2`$
ではなく $`D^2\times D^2`$ となる。ライブラリが付け加えるのはブロック構造そのものである：セクターごとに、
バッチとして格納されたクラウス作用素、異なる環の間の写像（`InterSectorChannel`）、
第3章の各トレースに関する双対、そしてセクター間の遷移行列である。

<a id="x-6"></a>

## 演習問題

講義ノートの[演習問題](../../lecture_notes/ja/06_quantum_channels.md#x-6)の数値的検証：

1. （ワーナー–ホレボ）`Channel.from_superoperator` を用いて、$`M_3`$ 上の $`\Phi(\rho)=\big(\mathrm{Tr}(\rho)1-\rho^T\big)/(k-1)`$
   がクラウス階数 $`3`$ の単位的チャネルであることを確かめよ。
2. $`M_2\oplus M_3`$ 上の位相緩和チャネルのブロック $`c`$ におけるクラウス階数が $`k_c`$ である一方、ライブラリは
   $`1+k_{\max}`$ 個のクラウス作用素を格納することを確かめよ。
3. クラウス作用素 $`X^aZ^b/k`$ をもつ混合ユニタリ・チャネルが $`M_3`$ 上の `center_expectation` であることを確かめよ。
4. ライブラリの三つのトレースについて $`(\Psi\circ\Phi)^\sharp=\Phi^\sharp\circ\Psi^\sharp`$ を確かめよ。

<details>
<summary>解答</summary>

*演習問題 1.*

```python
k = 3
M3 = TypeIAlgebra([k], [k], precision='double')
J_WH = choi(lambda X: (X.trace() * torch.eye(k, dtype=dt) - X.T) / (k - 1), k)
w = torch.linalg.eigvalsh(J_WH)
assert w.min() > -1e-12 and int((w > 1e-10).sum()) == k * (k - 1) // 2
WH = Channel.from_superoperator(M3, [to_superoperator(J_WH, k)])
assert WH.is_trace_preserving() and WH.is_unital()
r3 = states.random_density_matrix(M3).matrix[0, 0]
assert close(WH(M3.operator(r3[None, None])).matrix[0, 0], (torch.eye(k, dtype=dt) - r3.T) / (k - 1))
```

*演習問題 2.*

```python
ranks = [int((torch.linalg.eigvalsh(J_c[0]) > 1e-10).sum()) for J_c in zoo['dephasing'].choi()]
print(ranks, zoo['dephasing'].kraus_rank)                 # [2, 3] versus 4 stored operators
assert ranks == alg.k_factors and zoo['dephasing'].kraus_rank == 1 + max(alg.k_factors)
```

*演習問題 3.*

```python
Xs = torch.roll(torch.eye(k, dtype=dt), 1, dims=0)
Zs = torch.diag(torch.exp(2j * math.pi * torch.arange(k, dtype=torch.float64) / k)).to(dt)
W = [torch.linalg.matrix_power(Xs, a) @ torch.linalg.matrix_power(Zs, b) / k for a in range(k) for b in range(k)]
twirl = channels.kraus_channel([M3.operator(Wab[None, None]) for Wab in W])
r3op = states.random_density_matrix(M3)
assert close(twirl(r3op), channels.center_expectation(M3)(r3op))
```

*演習問題 4.*

```python
Psi2 = channels.random_inter_sector_channel(alg_out, alg, kraus_rank=2)
for t in trace:
    assert close((Psi2 @ Psi).adjoint(trace=t)(x), Psi.adjoint(trace=t)(Psi2.adjoint(trace=t)(x)))
```

</details>

**次へ：** [チュートリアル 7](07_tensor_products.md)：ライブラリにおけるテンソル積、部分トレース、エンタングルメント。
