[English](../10_quantum_computing.md) · [Español](../es/10_quantum_computing.md) · [Français](../fr/10_quantum_computing.md) · [Deutsch](../de/10_quantum_computing.md) · [中文](../zh/10_quantum_computing.md) · **日本語** · [Русский](../ru/10_quantum_computing.md)

# チュートリアル 10. 量子計算

このチュートリアルは講義ノートの[第10章](../../lecture_notes/ja/10_quantum_computing.md)に対応する。同章は、量子計算機上のノイズを、そのクラウス作用素の
$`*`$-環によって記述する。ここでは不動点環とノイズレス部分系を計算し（[定理 10.2](../../lecture_notes/ja/10_quantum_computing.md#r-10-2)）、繰り返し符号でクニル–ラフラムの定理
（[定理 10.5](../../lecture_notes/ja/10_quantum_computing.md#r-10-5)）とセクターのラベルの訂正（[命題 10.6](../../lecture_notes/ja/10_quantum_computing.md#r-10-6)）を検証し、ハミング重みのセクターの環
$`\bigoplus_wM_{\binom nw}`$ の中で、振幅減衰と事後選択を伴う XY ミキサーのアンザッツをシミュレートする
（[命題 10.8](../../lecture_notes/ja/10_quantum_computing.md#r-10-8)–[10.10](../../lecture_notes/ja/10_quantum_computing.md#r-10-10)）。

<a id="s-10-1"></a>

## 10.1 ノイズ環、ノイズレス部分系、不動点

*講義ノート：* [10.1節](../../lecture_notes/ja/10_quantum_computing.md#s-10-1)；[定理 10.2](../../lecture_notes/ja/10_quantum_computing.md#r-10-2)、[系 10.3](../../lecture_notes/ja/10_quantum_computing.md#r-10-3)、[命題 10.4](../../lecture_notes/ja/10_quantum_computing.md#r-10-4)。

[講義ノート、系 10.3](../../lecture_notes/ja/10_quantum_computing.md#r-10-3) により、$`n`$ 量子ビット上の集団的位相緩和 $`\Phi(\rho)=\tfrac12(U\rho U^*+U^*\rho U)`$、
$`U=e^{-i\theta Z_{\mathrm{tot}}}`$ の相互作用環は $`\mathcal A=\mathrm{span}\{P_0,\dots,P_n\}`$ であり（$`P_w`$ はハミング重み $`w`$ の
基底ベクトルが張る部分空間 $`H_w`$ への射影）、$`\mathrm{Fix}(\Phi)=\mathcal A'=\bigoplus_wB(H_w)`$ である（[定理 10.2](../../lecture_notes/ja/10_quantum_computing.md#r-10-2)）。

可換子環は $`T\mapsto(AT-TA)_A`$ の零空間であり（第2章）、不動点は
$`\sum_iE_i\otimes\bar E_i-1`$（行優先でベクトル化した行列に対する $`\Phi-\mathrm{id}`$ の行列）の零空間である。$`n=4`$ に対し、
系の三つの空間の次元はそれぞれ $`5`$、$`70=1+16+36+16+1`$、$`70`$ である：

```python
import itertools
import math
import torch

torch.manual_seed(0)
dt = torch.complex128
I2, X = torch.eye(2, dtype=dt), torch.tensor([[0, 1], [1, 0]], dtype=dt)
Y, Z = torch.tensor([[0, -1j], [1j, 0]], dtype=dt), torch.tensor([[1, 0], [0, -1]], dtype=dt)

def on(op, j, n):                             # op on qubit j of n (j = 0 is the leftmost factor)
    out = torch.ones(1, 1, dtype=dt)
    for i in range(n):
        out = torch.kron(out, op if i == j else I2)
    return out

def null_space(L, tol=1e-9):                  # orthonormal columns spanning {x : L x = 0}
    w, V = torch.linalg.eigh(L.mH @ L)
    return V[:, w < tol]

def commutant_basis(ops):
    d, ops = ops[0].shape[0], [A.resolve_conj().contiguous() for A in ops]
    eye = torch.eye(d, dtype=dt)
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    return [v.reshape(d, d) for v in null_space(L).T]

def fixed_points(kraus):                      # vec(E X E^*) = (E (x) conj E) vec(X)
    d = kraus[0].shape[0]
    S = sum(torch.kron(E.resolve_conj(), E.conj().resolve_conj()) for E in kraus)
    return [v.reshape(d, d) for v in null_space(S - torch.eye(d * d, dtype=dt)).T]

n, theta = 4, 0.3
U = torch.linalg.matrix_exp(-1j * theta * sum(on(Z, j, n) for j in range(n)))
kraus = [U / math.sqrt(2), U.mH / math.sqrt(2)]                  # unital: U U^* = U^* U = 1
comm = commutant_basis(kraus + [E.mH for E in kraus])             # A'
A_int, fix = commutant_basis(comm), fixed_points(kraus)           # A = A'' (Theorem 2.2), Fix(Phi)
print(len(A_int), len(comm), len(fix))
assert (len(A_int), len(comm), len(fix)) == (5, 70, 70)
span = torch.linalg.svdvals(torch.stack([F.reshape(-1) for F in comm + fix]))
assert int((span > 1e-9).sum()) == 70                             # Fix(Phi) = A' (Theorem 10.2)
```

$`\bigoplus_wM_{\binom nw}`$ は `SpinChain`（第9章）の環であり、そのセクターは上向きスピンの
個数、すなわちハミング重みである。そこではクラウス作用素は中心に属し、
$`U=\sum_we^{-i\theta(n-2w)}P_w`$ となり、チャネルは恒等写像である：ライブラリはノイズのない部分を格納している。

```python
from torch_vn_algebra import SpinChain, TypeIAlgebra, DensityMatrix, channels, states, dynamics

M = SpinChain(n).algebra                                          # (+)_w M_binom(4, w)
print(M.k_factors)
assert sum(k * k for k in M.k_factors) == len(fix)
phases = [torch.exp(torch.tensor(-1j * theta * (n - 2 * w), dtype=dt)) for w in range(n + 1)]
Phi_M = channels.kraus_channel([M.central([z / math.sqrt(2) for z in phases]),
                                M.central([z.conj() / math.sqrt(2) for z in phases])])
rho = states.random_density_matrix(M, batch_size=3)
assert (Phi_M(rho).matrix - rho.matrix).abs().max() < 1e-12      # no state of M is disturbed
```

集団回転 $`e^{-i\theta\,\mathbf n\cdot\mathbf S}`$（$`S_\alpha=\tfrac12\sum_j\sigma_{\alpha,j}`$）は、$`1,S_x,S_y,S_z`$ が生成する $`*`$-環
$`\mathcal B`$ に属する。$`n=4`$ に対し、講義ノートの[10.1節](../../lecture_notes/ja/10_quantum_computing.md#s-10-1)は $`\dim\mathcal B=35`$、$`\dim\mathcal B'=14`$、および
$`K=\ker S_x\cap\ker S_y\cap\ker S_z`$ について $`\dim K=2`$ を示す。コードは同じ値を見いだし、[命題 10.4](../../lecture_notes/ja/10_quantum_computing.md#r-10-4) を検証する：
ランダムな集団回転をクラウス作用素とするチャネルは、$`K`$ に台をもつ状態を固定し、一般の状態は動かす。

```python
S = [sum(on(P, j, n) for j in range(n)) / 2 for P in (X, Y, Z)]
comm_S = commutant_basis(S)
print(len(comm_S), len(commutant_basis(comm_S)))
assert len(comm_S) == 14 and len(commutant_basis(comm_S)) == 35
K = null_space(torch.cat(S))                                      # columns: a basis of K
assert K.shape[1] == 2
axes = torch.nn.functional.normalize(torch.randn(3, 3, dtype=torch.float64), dim=1)
rotations = [torch.linalg.matrix_exp(-1j * t * sum(c * Sa for c, Sa in zip(ax, S))) / math.sqrt(3)
             for t, ax in zip(torch.rand(3, dtype=torch.float64) * math.pi, axes)]
apply = lambda kr, r: sum(E @ r @ E.mH for E in kr)
for G in (K @ torch.randn(2, 2, dtype=dt), torch.randn(16, 16, dtype=dt)):  # supported in K / generic
    r = G @ G.mH / (G @ G.mH).trace()
    assert ((apply(rotations, r) - r).abs().max() < 1e-12) == (G.shape[1] == 2)
```

<a id="s-10-2"></a>

## 10.2 誤り訂正：クニル–ラフラムの定理

*講義ノート：* [10.2節](../../lecture_notes/ja/10_quantum_computing.md#s-10-2)；[定理 10.5](../../lecture_notes/ja/10_quantum_computing.md#r-10-5)。

講義ノートの[10.2節](../../lecture_notes/ja/10_quantum_computing.md#s-10-2)のビット反転の例：$`\mathcal C=\mathrm{span}\{|000\rangle,|111\rangle\}`$ と
「高々一つのビット反転」に対し、[定理 10.5](../../lecture_notes/ja/10_quantum_computing.md#r-10-5) のクニル–ラフラム条件が $`\Lambda=\mathrm{diag}(p_0,\dots,p_3)`$ で成り立ち、
回復 $`\mathcal R(\sigma)=\sum_{k=0}^3PX_k\sigma X_kP`$ は符号のすべての状態を復元する。各量子ビットが確率 $`p`$ で独立に
反転する場合には条件は成り立たず、回復後の $`|0_L\rangle`$ の忠実度は $`p=0.1`$ に対して $`1-3p^2+2p^3=0.972`$ である。補助関数
`kl_matrix` は $`\Lambda`$ を返し、条件が成り立たなければ `None` を返す。ライブラリでは $`M_8`$ はブロックが一つの環であり、
`R @ Phi` は合成 $`\mathcal R\circ\Phi`$ である：

```python
M8 = TypeIAlgebra([8], [8], precision='double')
op8 = lambda T: M8.operator(T.reshape(1, 1, 8, 8))
state8 = lambda v: DensityMatrix(M8, matrix=torch.outer(v, v.conj()).reshape(1, 1, 8, 8))
ket = lambda bits: torch.eye(8, dtype=dt)[int(bits, 2)]
def kl_matrix(kraus, P, tol=1e-12):
    """(lambda_ij) with P E_i^* E_j P = lambda_ij P, or None if the Knill-Laflamme condition fails."""
    lam = torch.tensor([[(P @ a.mH @ b @ P).trace() / P.trace() for b in kraus] for a in kraus])
    ok = all((P @ a.mH @ b @ P - lam[i, j] * P).abs().max() < tol
             for i, a in enumerate(kraus) for j, b in enumerate(kraus))
    return lam if ok else None
def independent(op, prob):                    # op on each qubit independently with probability prob
    return [math.sqrt(prob ** sum(s) * (1 - prob) ** (3 - sum(s)))
            * torch.linalg.multi_dot([on(op, j, 3) if s[j] else torch.eye(8, dtype=dt) for j in range(3)])
            for s in itertools.product([0, 1], repeat=3)]

P = torch.outer(ket('000'), ket('000')) + torch.outer(ket('111'), ket('111'))
flips = [torch.eye(8, dtype=dt)] + [on(X, j, 3) for j in range(3)]
E_one = [math.sqrt(pk) * Xk for pk, Xk in zip((0.7, 0.1, 0.1, 0.1), flips)]
assert torch.allclose(kl_matrix(E_one, P), torch.diag(torch.tensor([0.7, 0.1, 0.1, 0.1], dtype=dt)))
Phi_one = channels.kraus_channel([op8(E) for E in E_one])
R = channels.kraus_channel([op8(P @ Xk) for Xk in flips])        # syndrome measurement + correction
assert R.is_trace_preserving()
a, b = torch.nn.functional.normalize(torch.randn(2, dtype=dt), dim=0)
rho_L = state8(a * ket('000') + b * ket('111'))
assert abs((R @ Phi_one)(rho_L).fidelity(rho_L).item() - 1) < 1e-12       # Theorem 10.5
p = 0.1                                                           # independent flips: KL fails
E_ind = independent(X, p)
assert kl_matrix(E_ind, P) is None
zero_L = state8(ket('000'))
F = (R @ channels.kraus_channel([op8(E) for E in E_ind]))(zero_L).fidelity(zero_L).item()
print(f"fidelity after recovery {F:.4f}, exact {1 - 3 * p**2 + 2 * p**3:.4f}")
assert abs(F - (1 - 3 * p**2 + 2 * p**3)) < 1e-12
```

<a id="s-10-3"></a>

## 10.3 作用素環による誤り訂正

*講義ノート：* [10.3節](../../lecture_notes/ja/10_quantum_computing.md#s-10-3)；[命題 10.6](../../lecture_notes/ja/10_quantum_computing.md#r-10-6)。

講義ノートの[10.3節](../../lecture_notes/ja/10_quantum_computing.md#s-10-3)の例：独立な位相反転（量子ビットごとの確率 $`q`$）を受け、その後に高々一つのビット反転を受ける
繰り返し符号。量子ビットは保護されない（クニル–ラフラム条件は成り立たない）が、古典ビットは保護される：
[命題 10.6](../../lecture_notes/ja/10_quantum_computing.md#r-10-6) の条件 2 が $`P_0=|000\rangle\langle000|`$、$`P_1=|111\rangle\langle111|`$ について成り立ち、多数決
$`N_0=\sum_{|b|\le1}|b\rangle\langle b|`$、$`N_1=1-N_0`$ がビットを読み出す。10.2節の回復の後、$`|{+_L}\rangle`$ の忠実度は
$`q=0.05`$ に対して $`(1+(1-2q)^3)/2=0.8645`$ である。

```python
q = 0.05
Phi2 = Phi_one @ channels.kraus_channel([op8(E) for E in independent(Z, q)])   # dephasing, then a flip
E2 = [Phi2.kraus[0, i, 0] for i in range(Phi2.kraus_rank)]        # its 32 Kraus operators
P0, P1 = torch.outer(ket('000'), ket('000')), torch.outer(ket('111'), ket('111'))
assert kl_matrix(E2, P) is None                                   # the qubit is not protected
assert max((P0 @ e.mH @ f @ P1).abs().max().item() for e in E2 for f in E2) == 0      # condition 2
N0 = torch.diag(torch.tensor([float(bin(i).count('1') <= 1) for i in range(8)], dtype=dt))
for N, Pd in ((N0, P0), (torch.eye(8, dtype=dt) - N0, P1)):       # P Phi^*(N_d) P = P_d
    assert (P @ sum(e.mH @ N @ e for e in E2) @ P - Pd).abs().max() < 1e-12
votes = states.born_probabilities(Phi2(rho_L), [op8(N0), op8(torch.eye(8, dtype=dt) - N0)])[0]
assert torch.allclose(votes, torch.stack([abs(a) ** 2, abs(b) ** 2]))      # the classical bit survives
plus_L = state8((ket('000') + ket('111')) / math.sqrt(2))
F_plus = (R @ Phi2)(plus_L).fidelity(plus_L).item()
print(f"fidelity of |+_L> {F_plus:.4f}, exact {(1 + (1 - 2 * q)**3) / 2:.4f}")
assert abs(F_plus - (1 + (1 - 2 * q) ** 3) / 2) < 1e-12
```

<a id="s-10-4"></a>

## 10.4 制約付きミキサー

*講義ノート：* [10.4節](../../lecture_notes/ja/10_quantum_computing.md#s-10-4)；[命題 10.8](../../lecture_notes/ja/10_quantum_computing.md#r-10-8)、[命題 10.9](../../lecture_notes/ja/10_quantum_computing.md#r-10-9)、[命題 10.10](../../lecture_notes/ja/10_quantum_computing.md#r-10-10)。

量子交互作用素アンザッツは $`e^{-i\gamma H_C}`$（$`H_C`$ は対角なコスト・ハミルトニアン）と $`e^{-i\beta H_M}`$
（$`H_M=\tfrac14\sum_{(i,j)}(X_iX_j+Y_iY_j)`$ はリング上の XY ミキサー）を交互に作用させる。[講義ノート、命題 10.8](../../lecture_notes/ja/10_quantum_computing.md#r-10-8) により、
これらはすべて $`\mathcal M=\{P_0,\dots,P_n\}'\cong\bigoplus_wM_{\binom nw}`$ に属するので、アンザッツが初期状態のハミング重みの
セクターを離れることはない。

ライブラリでは $`\mathcal M`$ は `SpinChain(n).algebra` である（ビット $`1`$ = 上向きスピンなので $`S^z_i=-\tfrac12Z_i`$）。$`H_M`$
はその XXZ ハミルトニアン $`J\sum_{(i,j)}\big(\tfrac12(S^+_iS^-_j+S^-_iS^+_j)+\Delta S^z_iS^z_j\big)`$ で $`J=1`$、
$`\Delta=0`$ としたものである。実際 $`S^+_iS^-_j+S^-_iS^+_j=\tfrac12(X_iX_j+Y_iY_j)`$ である。これを稠密な
$`64\times64`$ 行列と照合する。$`H_C`$ としては、辺 $`(2,3)`$、$`(5,0)`$ で結ばれた三角形 $`\{0,1,2\}`$、$`\{3,4,5\}`$
上のグラフ二分割をとる：$`H_C=\sum_{(i,j)\in E}\tfrac12(1-Z_iZ_j)`$（$`Z_iZ_j=4S^z_iS^z_j`$）はカットされる
辺の数を数える。最良の二分割は $`2`$ 本の辺をカットする（二つの三角形に分けるもの。三つの頂点からなる他のどの集合も両方の三角形を
分割し、分割された三角形は二本のカットされる辺をもつ）。我々は $`6`$ 本をカットする $`|101010\rangle`$ から出発する。最後に、
アンザッツ全体を稠密な $`64\times64`$ 行列で、ライブラリとは独立に計算し直す。

```python
n, w = 6, 3
chain = SpinChain(n, 'periodic')
M = chain.algebra                                                 # (+)_w M_binom(6, w), padded to 20
H_M = chain.xxz(J=1.0, Delta=0.0)
H_dense = sum(on(X, i, n) @ on(X, j, n) + on(Y, i, n) @ on(Y, j, n) for i, j in chain.bonds) / 4
spec = torch.sort(torch.cat([e[0] for e in H_M.eigenvalues()]))[0]
assert torch.allclose(spec, torch.linalg.eigvalsh(H_dense))       # same spectrum, sector by sector
edges = chain.bonds + [(0, 2), (3, 5)]                            # two triangles joined by (2,3), (5,0)
H_C = M.zero()
for i, j in edges:
    H_C = H_C + (M.identity() - chain.sz(i) @ chain.sz(j) * 4.0) * 0.5
layers = [dynamics.propagator(H_M, beta) @ dynamics.propagator(H_C, gamma)
          for gamma, beta in ((0.6, 1.0), (1.4, 1.4))]
rho0 = ideal = chain.basis_state('101010')
for U_l in layers:
    ideal = channels.unitary_channel(U_l)(ideal)
assert abs(ideal.sector_probabilities()[0, w].item() - 1) < 1e-12           # Proposition 10.8

HC_dense = sum(torch.eye(64, dtype=dt) - on(Z, i, n) @ on(Z, j, n) for i, j in edges) / 2
psi = torch.eye(64, dtype=dt)[int('101010', 2)]                   # the same circuit on 6 full qubits
for gamma, beta in ((0.6, 1.0), (1.4, 1.4)):
    psi = torch.linalg.matrix_exp(-1j * beta * H_dense) @ torch.linalg.matrix_exp(-1j * gamma * HC_dense) @ psi
weight = torch.tensor([bin(b).count('1') for b in range(64)])
assert abs(psi[weight == w].norm().item() - 1) < 1e-12           # Proposition 10.8, dense
assert abs((psi.conj() @ HC_dense @ psi).real.item() - states.expectation(ideal, H_C).real.item()) < 1e-10
```

量子ビット $`i`$ の **振幅減衰**（確率 $`g`$ で $`|1\rangle\to|0\rangle`$）はクラウス作用素
$`K_0^{(i)}=|0\rangle\langle0|_i+\sqrt{1-g}\,|1\rangle\langle1|_i`$ と $`K_1^{(i)}=\sqrt g\,|0\rangle\langle1|_i`$ をもち、後者は
重みを下げる：これはセクター間の写像 `chain.site_amplitude_damping(i, g)`（`InterSectorChannel`）である。

[講義ノート、命題 10.9](../../lecture_notes/ja/10_quantum_computing.md#r-10-9) により、回路の各層の後にすべての量子ビットの減衰が続く場合、$`p`$ 層の後でも
初期状態の重み $`w`$ は確率 $`(1-g)^{wp}`$ で見いだされ、事後選択された状態はノイズのない状態である。

`InterSectorChannel` は出力の対角ブロックだけを保持する。ここではそれ以外のブロックはない（命題 10.9 の証明は、$`s\in\{0,1\}`$ として一つの量子ビットの減衰に
そのまま適用できる）ので、何も失われない。$`g=0.05`$、$`p=2`$ のとき受理確率は $`0.95^6=0.7351`$ である：

```python
g = 0.05
damping = [chain.site_amplitude_damping(i, g) for i in range(n)]
noisy = rho0
for U_l in layers:
    noisy = channels.unitary_channel(U_l)(noisy)
    for D in damping:                                             # one qubit after the other
        noisy = D(noisy)
print(noisy.sector_probabilities()[0])                            # weight leaks to w = 2, 1, 0
P_w = M.central([1.0 if v == w else 0.0 for v in range(n + 1)])
post, accept = noisy.condition_on(P_w)                            # measure the weight, keep w = 3
cut = lambda r: states.expectation(r, H_C).real.item()
print(f"P(weight 3) {accept.item():.4f}, fidelity {noisy.fidelity(ideal).item():.4f}, "
      f"post-selected {post.fidelity(ideal).item():.4f}; "
      f"cut: ideal {cut(ideal):.4f}, noisy {cut(noisy):.4f}, post-selected {cut(post):.4f}")
assert abs(accept.item() - (1 - g) ** (w * len(layers))) < 1e-12              # Proposition 10.9
assert abs(post.fidelity(ideal).item() - 1) < 1e-10 and abs(cut(post) - cut(ideal)) < 1e-10
```

*数値的観察*（出力される値）：実行不可能な部分は期待カット数を $`4.3789`$ から $`4.3123`$ に下げる。これは見かけ上の改善であり、
角度の最適化器を誤らせるであろう。事後選択の後では値は厳密である（命題 10.9）。

[講義ノート、命題 10.10](../../lecture_notes/ja/10_quantum_computing.md#r-10-10) により、振幅減衰（レート $`\gamma_1`$）と位相緩和 $`L_k=S^z_k`$（レート $`\gamma_2`$）を
もつリンドブラッド方程式のもとで、$`H_w`$ から出発した状態は確率 $`e^{-\gamma_1wt}`$ で重み $`w`$ を保ち、事後選択された状態は
$`\gamma_1`$ に依存しない。$`\gamma_2=0`$ ならば事後選択は完全である。解説論文
（`examples/paper/ex9_qaoa_xy_mixer.py`）は $`|111000\rangle`$ を $`H_M`$ のもとで $`T=2`$ まで時間発展させる。我々はその
$`\gamma_1=\gamma_2=0.1`$ の行を再現する：$`P(w=3)=e^{-0.6}=0.5488`$（命題 10.10）であり、下のコードの *数値的観察* として、忠実度は $`0.4805`$、
事後選択後の忠実度は $`0.8755`$ である。命題 10.10 により、後者は $`\gamma_1=0`$ の場合と同じである。

```python
T = 2.0
jumps = [chain.lowering(i) for i in range(n)] + [chain.sz(i) for i in range(n)]
rho0 = chain.basis_state('111000')                                # the initial state of the paper
ideal_T = dynamics.von_neumann(rho0, H_M, [T])[0]
post_T = {}
for g1, g2 in ((0.1, 0.0), (0.0, 0.1), (0.1, 0.1)):          # the last one is printed
    rho_T = dynamics.lindblad_evolve(rho0, H_M, jumps, [0.0, T], rates=[g1] * n + [g2] * n, substeps=80)[-1]
    post_T[g1, g2], accept_T = rho_T.condition_on(P_w)
    assert abs(accept_T.item() - math.exp(-w * g1 * T)) < 1e-8                 # Proposition 10.10
print(f"P(weight 3) {accept_T.item():.4f}, fidelity {rho_T.fidelity(ideal_T).item():.4f}, "
      f"post-selected {post_T[0.1, 0.1].fidelity(ideal_T).item():.4f}")
assert abs(post_T[0.1, 0.0].fidelity(ideal_T).item() - 1) < 1e-6               # gamma_2 = 0: perfect
assert post_T[0.1, 0.1].trace_distance(post_T[0.0, 0.1]).item() < 1e-6       # independent of gamma_1
```

<a id="s-10-5"></a>

## 10.5 比較の対象

Qiskit（`qiskit.quantum_info`、Aer のノイズモデル、カスタム・ミキサー付きの `QAOAAnsatz`）、Cirq、PennyLane
（`qml.qaoa.xy_mixer`、`default.mixed`）、および QuTiP/`qutip-qip` は、これらすべてを $`n`$ 個の完全な量子ビット上でシミュレートし、
回路のための適切な道具である。stim はパウリ型ノイズをもつクリフォード回路を数千
量子ビットでサンプリングする（振幅減衰はパウリ・チャネルではない：パウリ・チャネルは単位的であるが、振幅減衰はそうではない）。我々の知る限り、これらのパッケージには、
チャネルの相互作用環、その可換子環、あるいはその不動点環全体を計算する関数はない（QuTiP の
`steadystate` はリンドブラッド生成子の定常状態を一つ返す）。行列 $`*`$-環の数値的なブロック分解については
Murota et al., *Japan J. Indust. Appl. Math.* 27 (2010) を参照。ライブラリが付け加えるのは、
$`\bigoplus_wM_{\binom nw}`$ での格納、セクター間の写像、中心射影としての事後選択であり、いずれもバッチ化されている。
一つのセクターは $`4^n`$ 個ではなく $`\binom nw^2`$ 個の成分を必要とする（$`n=16`$、
$`w=8`$ で $`1.7\cdot10^8`$ 対 $`4.3\cdot10^9`$）が、環全体は $`(n+1)\binom n{\lfloor n/2\rfloor}^2`$ 成分（$`2.8\cdot10^9`$）にパディングされるので、
節約は控えめである。ライブラリは回路シミュレータではなく、スタビライザー形式ももたない。

<a id="x-10"></a>

## 演習問題

講義ノートの[演習問題 1–3](../../lecture_notes/ja/10_quantum_computing.md#x-10)の数値的検証：一つの量子ビット上の振幅減衰について、$`\mathcal A'`$、
$`\mathrm{Fix}(\Phi)`$、$`\mathrm{Fix}(\Phi^*)`$ は一次元である（演習 1）。[系 10.3](../../lecture_notes/ja/10_quantum_computing.md#r-10-3) のチャネルで $`\theta=\pi/2`$、
$`n=4`$ としたものについて $`\dim\mathrm{Fix}(\Phi)=128`$ である（演習 2）。ビット反転符号の回復 $`\mathcal R`$ はコヒーレントな誤り
$`e^{-i\phi X_1}`$ を訂正する（演習 3）。

<details>
<summary>解答</summary>

```python
E_ad = [torch.tensor([[1, 0], [0, math.sqrt(0.7)]], dtype=dt), torch.tensor([[0, math.sqrt(0.3)], [0, 0]], dtype=dt)]
fix_ad = fixed_points(E_ad)
assert len(fix_ad) == 1 and len(commutant_basis(E_ad + [E.mH for E in E_ad])) == 1
assert torch.allclose(fix_ad[0] / fix_ad[0][0, 0], torch.tensor([[1, 0], [0, 0]], dtype=dt))
assert len(fixed_points([E.mH for E in E_ad])) == 1                     # Fix(Phi^*) = C 1
```

```python
U2 = torch.linalg.matrix_exp(-1j * (math.pi / 2) * sum(on(Z, j, 4) for j in range(4)))
assert len(fixed_points([U2 / math.sqrt(2), U2.mH / math.sqrt(2)])) == 128
assert sum(math.comb(4, v) ** 2 for v in range(5)) == math.comb(8, 4) == 70
```

```python
for phi in (0.2, 1.0, 2.5):
    G1 = torch.linalg.matrix_exp(-1j * phi * on(X, 0, 3))
    out = (R @ channels.kraus_channel([op8(G1)]))(rho_L)
    assert (out.matrix - rho_L.matrix).abs().max() < 1e-12       # corrected exactly
```

</details>

**次へ：** [チュートリアル 11](11_optimisation.md)：環のユニタリ群上の最適化と、トレース不等式。
