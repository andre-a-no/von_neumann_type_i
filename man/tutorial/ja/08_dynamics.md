[English](../08_dynamics.md) · [Español](../es/08_dynamics.md) · [Français](../fr/08_dynamics.md) · [Deutsch](../de/08_dynamics.md) · [中文](../zh/08_dynamics.md) · **日本語** · [Русский](../ru/08_dynamics.md)

# チュートリアル 8. ダイナミクス

このチュートリアルは講義ノートの[第8章](../../lecture_notes/ja/08_dynamics.md)に対応する。同章では、$`\mathcal M`$ のハミルトニアンが
セクターの重みを保存すること（[定理 8.3](../../lecture_notes/ja/08_dynamics.md#r-8-3)）、GKSL 形の生成子が量子動力学的半群を生み出すこと
（[定理 8.6](../../lecture_notes/ja/08_dynamics.md#r-8-6)）を証明し、セクター間のジャンプ、事後選択、定常状態を扱う。ここでは `dynamics` の
厳密なソルバーとルンゲ–クッタ・ソルバーを用いてシュレーディンガー方程式、フォン・ノイマン方程式、リンドブラッド方程式を解き、
結果を SciPy および $`M_D`$ 上の GKSL 生成子と照合する。

<a id="s-8-1"></a>

## 8.1 閉じた系

*講義ノート：* [8.1節](../../lecture_notes/ja/08_dynamics.md#s-8-1)；[命題 8.2](../../lecture_notes/ja/08_dynamics.md#r-8-2)（時間発展演算子 $`U(t)=e^{-itH}`$）、[定理 8.3](../../lecture_notes/ja/08_dynamics.md#r-8-3)（電荷の保存）。

命題 8.2(1) により、時間発展はブロックごとに一つずつ、サイズ $`k_c`$ の $`C`$ 個の独立な問題に分かれる。
ライブラリはブロックサイズごとに一回の固有値分解で $`U(t)`$ を計算し（`dynamics.propagator`）、
`dynamics.von_neumann` は固有基底において密度を多数の時刻へ一度に発展させる。固有基底では
$`\rho(t)_{ij}=e^{-it(\lambda_i-\lambda_j)}\rho_{ij}`$ である。`scipy.linalg.expm` と照合する：

```python
import math
import numpy as np
import scipy.integrate
import scipy.linalg as sla
import torch
from torch_vn_algebra import TypeIAlgebra, DensityMatrix, InterSectorChannel, channels, dynamics, states

torch.manual_seed(0)
f64, c128 = torch.float64, torch.complex128
alg = TypeIAlgebra([2, 3, 4], [2, 3, 4], precision='double')     # M_2 (+) M_3 (+) M_4, D = 9
H = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
t = 0.7
U = dynamics.propagator(H, t)                                      # exp(-itH), block by block
for c, k in enumerate(alg.k_factors):                              # independent reference
    ref = sla.expm(-1j * t * H.matrix[0, c, :k, :k].numpy())
    assert np.abs(U.matrix[0, c, :k, :k].numpy() - ref).max() < 1e-12
assert (U.adjoint() @ U - alg.identity()).matrix.abs().max() < 1e-12

rho0 = states.random_density_matrix(alg)
times = [0.0, 0.5, 1.0, 2.0]
traj = dynamics.von_neumann(rho0, H, times)                        # exact, all times at once
for s, rho_s in zip(times, traj):
    U_s = dynamics.propagator(H, s)
    assert (rho_s.matrix - (U_s @ rho0 @ U_s.adjoint()).matrix).abs().max() < 1e-12
    assert (rho_s.sector_probabilities() - rho0.sector_probabilities()).abs().max() < 1e-12   # Theorem 8.3
A = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
A_t = U.adjoint() @ A @ U                                          # Heisenberg picture
assert abs(rho0.expectation(A_t) - dynamics.von_neumann(rho0, H, [t])[0].expectation(A)) < 1e-12
```

単位ベクトル $`\psi=\bigoplus_c\psi_c`$ は $`\mathcal M`$ 上の状態 $`\bigoplus_c|\psi_c\rangle\langle\psi_c|`$ を定め、これは
$`U(t)`$ によって時間発展する（[講義ノート、8.1節](../../lecture_notes/ja/08_dynamics.md#s-8-1)）。したがって `dynamics.schrodinger`（ブロック配置 `(batch, C, k_max)` のベクトル）と `von_neumann`
はその上で一致しなければならない。時間依存の $`H(t)`$ に対しては、
ライブラリは固定刻みの古典的ルンゲ–クッタ法（`dynamics.rk4`、`schrodinger_rk4`）を用いる：一つのステップは
$`\psi`$ に $`H(t)\in\mathcal M`$ のいくつかの値の多項式を掛けるので、セクター間で振幅を移すことは決してないが、
ユニタリではない：定数の $`H`$ に対し、刻み幅 $`h`$ の一ステップは固有値 $`\lambda`$ に属する $`H`$ の固有ベクトルに
$`T(-ih\lambda)`$ を掛ける。ここで $`T(z)=\sum_{k=0}^4z^k/k!`$ であり、$`0<|x|<2\sqrt2`$ に対し
$`|T(ix)|^2=1-\frac{x^6}{72}+\frac{x^8}{576}<1`$ である（展開せよ）。したがって RK4 はセクターの重みを厳密には
保存しない：定数の $`H`$ に対し、$`0<|h\lambda|<2\sqrt2`$ を満たす固有ベクトル上に零でない係数をもつセクター成分は、
各ステップで重みを失う。一回の実行でそれらがどれだけ変化するかは数値的観察である（下記）。
（一つの大きなセクターに対しては、`krylov.evolve` は行列–ベクトル積しか必要としない。[チュートリアル 9、9.5節](../../tutorial/ja/09_spin_chains.md#s-9-5)を参照。）

```python
rows = torch.tensor([[1, 1, 0, 0], [1, 1, 1, 0], [1, 1, 1, 1]], dtype=c128)   # entries of C^2 (+) C^3 (+) C^4
psi0 = torch.randn(1, alg.C, alg.k_max, dtype=c128) * rows
psi0 = psi0 / psi0.norm()
psi = dynamics.schrodinger(psi0, H, times)                         # shape (T, batch, C, k_max)
p_vec = (psi.abs() ** 2).sum(-1)                                   # ||P_c psi(t)||^2
assert (p_vec - p_vec[0]).abs().max() < 1e-12
mix = DensityMatrix(alg, matrix=psi0[..., :, None] * psi0[..., None, :].conj())   # the state of psi0 on M
outer = psi[..., :, None] * psi[..., None, :].conj()
assert all((o - r.matrix).abs().max() < 1e-12 for o, r in zip(outer, dynamics.von_neumann(mix, H, times)))

V = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
H_of_t = lambda s: H.matrix + math.cos(3 * s) * V.matrix          # driven, in M at every time
psi_rk = dynamics.schrodinger_rk4(psi0, H_of_t, times, substeps=50, algebra=alg)
drift = ((psi_rk.abs() ** 2).sum(-1) - p_vec[0]).abs().max().item()
print(f"RK4, time-dependent H: largest change of a sector weight {drift:.1e}")
assert drift < 1e-8 and (dynamics.schrodinger_rk4(psi0, H, times, substeps=50) - psi).abs().max() < 1e-7

```

*数値的観察*：この実行では、セクターの重みは RK4 のもとでは $`1.1\cdot10^{-9}`$（出力される値）だけ変化し、
厳密なソルバーのもとでは assert の許容誤差 $`10^{-12}`$ 未満しか変化しない。厳密なソルバーについては、
定理 8.3 が厳密な保存を与える。

<a id="s-8-2"></a>

## 8.2 量子動力学的半群

*講義ノート：* [8.2節](../../lecture_notes/ja/08_dynamics.md#s-8-2)；[定理 8.6](../../lecture_notes/ja/08_dynamics.md#r-8-6)（GKSL 生成子とダイソン展開）、[定理 8.7](../../lecture_notes/ja/08_dynamics.md#r-8-7)。

ライブラリには四つのリンドブラッド・ルーチンがある：`lindblad_rhs` は右辺を
`rk4` 用の関数として返す。`lindblad_evolve` は RK4 で積分する。`lindblad_superoperator` は各セクターについて、
行優先の $`\mathrm{vec}\,\rho_c`$ に作用する $`\mathcal L`$ の $`k_c^2\times k_c^2`$ 行列を返す。そして
`lindblad_channel` は $`e^{t\mathcal L}`$ をセクターごとに厳密に計算し、`Channel` として返す。
$`M_3`$ 上の任意のジャンプについて定理 8.6 を確かめる：$`e^{t\mathcal L}`$ のチョイ行列は正であり（定理
6.2）、写像は TP であり、半群法則が成り立ち、RK4 は指数関数と一致し、ジャンプ部分
$`e^{t\mathcal L}-S_0(t)=\sum_{n\ge1}S_n(t)`$ それ自体が CP である：

```python
M3 = TypeIAlgebra([3], [3], precision='double')
H3 = M3.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
jumps = [M3.operator(torch.randn(1, 1, 3, 3, dtype=c128) / 2) for _ in range(2)]     # arbitrary L_j
T = lambda s: dynamics.lindblad_channel(M3, H3, jumps, t=s)        # exp(sL), exact
S = lambda s: T(s).superoperator()[0][0]                           # its 9 x 9 matrix on row-major vec
choi = lambda X: X.reshape(3, 3, 3, 3).permute(2, 0, 3, 1).reshape(9, 9)   # J_{(i,a),(j,b)} = S_{(a,b),(i,j)}
assert (choi(S(0.8)) - T(0.8).choi()[0][0]).abs().max() < 1e-12
assert torch.linalg.eigvalsh(choi(S(0.8))).min() > -1e-12         # completely positive
assert T(0.8).trace_preservation_error().max() < 1e-12             # trace preserving
assert (S(0.3) @ S(0.5) - S(0.8)).abs().max() < 1e-12              # semigroup

rho3 = states.random_density_matrix(M3)
rk = dynamics.lindblad_evolve(rho3, H3, jumps, [0.0, 0.8], substeps=40)[-1]
assert (rk.matrix - T(0.8)(rho3).matrix).abs().max() < 1e-7        # RK4 against exp(tL)
assert abs(dynamics.lindblad_rhs(H3, jumps)(0.0, rho3.matrix)[0, 0].trace()) < 1e-12   # Tr L(rho) = 0

H_eff = H3.matrix - 0.5j * sum(L.adjoint().matrix @ L.matrix for L in jumps)
S0 = channels.Channel(M3, torch.linalg.matrix_exp(-0.8j * H_eff)[:, None]).superoperator()[0][0]
jump_part = torch.linalg.eigvalsh(choi(S(0.8) - S0))
print(f"smallest eigenvalue of the Choi matrix of the jump part: {jump_part.min().item():.1e}")
assert jump_part.min() > -1e-12                                    # sum_{n >= 1} S_n(t) is CP
```

<a id="s-8-3"></a>

## 8.3 セクターをもつ環におけるリンドブラッド方程式

*講義ノート：* [8.3節](../../lecture_notes/ja/08_dynamics.md#s-8-3)：`dynamics.lindblad_rhs` が実装する、セクター内およびセクター間のジャンプをもつ生成子 $`\mathcal L`$；[命題 8.8](../../lecture_notes/ja/08_dynamics.md#r-8-8)、[命題 8.9](../../lecture_notes/ja/08_dynamics.md#r-8-9)（セクターの重み）。

**例：粒子損失**（[講義ノート、8.3節](../../lecture_notes/ja/08_dynamics.md#s-8-3)）. 粒子数が高々二の二つのボース・モード、
$`\mathcal M=M_1\oplus M_2\oplus M_3`$（セクター $`N=0,1,2`$）、ボース–ハバード二量体
$`H=J(a_1^*a_2+a_2^*a_1)+\frac U2\sum_jn_j(n_j-1)`$ を考え、各モードはレート $`\gamma`$ で粒子を失うとする：これはクラウス作用素
$`a_1,a_2`$ を $`N\to N-1`$ に制限した一つの `InterSectorChannel` である。
[命題 8.9(2)](../../lecture_notes/ja/08_dynamics.md#r-8-9) により、セクターの重みは $`p'=Qp`$ に従う。ここでレート $`1`$ に対し $`Q=R-\mathrm{diag}(\sum_dR_{dc})_c`$、
$`R`$ = `transition_matrix()` である。$`N=2`$ から出発すると、重みはホッピングと相互作用が何をしようとも二項分布 $`p_2=q^2`$、$`p_1=2q(1-q)`$、$`p_0=(1-q)^2`$、
$`q=e^{-\gamma t}`$ に従う。コードは
$`\Phi^*(1)=0\oplus1\oplus2`$、すなわち用いるブロック上で $`a_1^*a_1+a_2^*a_2=N`$ であることを確かめ、$`p(t)`$ を両方の公式と照合する。

```python
fock = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double', charges=[0, 1, 2])   # N = 0, 1, 2
r2 = math.sqrt(2)
a1 = {(0, 1): [[1.0, 0.0]], (1, 2): [[r2, 0.0, 0.0], [0.0, 1.0, 0.0]]}           # a_1 : N -> N - 1
a2 = {(0, 1): [[0.0, 1.0]], (1, 2): [[0.0, 1.0, 0.0], [0.0, 0.0, r2]]}           # a_2 : N -> N - 1
blocks = lambda a: {dc: torch.tensor(K, dtype=f64) for dc, K in a.items()}
loss = InterSectorChannel.from_blocks(fock, fock, {dc: [blocks(a1)[dc], blocks(a2)[dc]] for dc in a1})
J, U_int = 1.0, 0.8                                                # hopping and interaction
H_bh = fock.from_blocks([torch.zeros(1, 1, dtype=f64), torch.tensor([[0.0, J], [J, 0.0]], dtype=f64),
                         torch.tensor([[U_int, r2 * J, 0.0], [r2 * J, 0.0, r2 * J], [0.0, r2 * J, U_int]], dtype=f64)])

assert (loss.adjoint()(fock.identity()).matrix - fock.central([0.0, 1.0, 2.0]).matrix).abs().max() < 1e-12
R = loss.transition_matrix()[0]                                    # r_dc (columns: from sector c)
Q = R - torch.diag(R.sum(0))                                       # Proposition 8.9(2), for rate 1
print(Q)

gamma = 0.5
psi2 = torch.randn(3, dtype=c128)
psi2 = psi2 / psi2.norm()
m0 = torch.zeros(1, 3, 3, 3, dtype=c128)
m0[0, 2] = torch.outer(psi2, psi2.conj())                          # a pure state with N = 2
ts = [0.0, 0.5, 1.0, 2.0, 4.0]
traj = dynamics.lindblad_evolve(DensityMatrix(fock, matrix=m0), H_bh, [loss], ts, rates=[gamma], substeps=100)
for s, rho_s in zip(ts, traj):
    q = math.exp(-gamma * s)
    p = rho_s.sector_probabilities()[0]
    assert (p - torch.tensor([(1 - q) ** 2, 2 * q * (1 - q), q ** 2], dtype=f64)).abs().max() < 1e-8
    assert (p - torch.from_numpy(sla.expm(gamma * s * Q.numpy())[:, 2])).abs().max() < 1e-8
neg = min(w.min().item() for rho_s in traj for w in rho_s.eigenvalues())
print(f"smallest eigenvalue of rho(t) under RK4: {neg:.1e}")
assert neg > -1e-7                                                 # positive up to the RK4 error

sm = torch.tensor([[0.0, 1.0], [0.0, 0.0]], dtype=c128)           # |g><e|; H = 0, rate 3, one step h = 1
decay = lambda s, X: 3.0 * (sm @ X @ sm.mH - 0.5 * (sm.mH @ sm @ X + X @ sm.mH @ sm))
step = dynamics.rk4(decay, torch.diag(torch.tensor([0.0, 1.0], dtype=c128)), [0.0, 1.0], substeps=1)[-1]
assert (step - torch.diag(torch.tensor([-3 / 8, 11 / 8], dtype=c128))).abs().max() < 1e-14   # not positive
```

RK4 は正値性を保つことが保証されていない：$`\rho'=\mathcal L\rho`$ に対し、刻み幅 $`h`$ の一ステップは
8.1節の多項式 $`T`$ による $`\rho\mapsto T(h\mathcal L)\rho`$ である（四つの段を代入せよ）。これはトレースを保存する
（すべての $`X`$ に対し $`\mathrm{Tr}\,\mathcal L(X)=0`$ であるから）が、正値である必要はない。たとえば、基底
$`|g\rangle,|e\rangle`$ をもつ $`M_2`$ 上で、$`H=0`$、レート $`\gamma`$ の一つのジャンプ $`|g\rangle\langle e|`$ とすると、$`\mathcal L`$ は
$`\mathrm{diag}(a,b)`$ を $`\gamma b\mathrm{diag}(1,-1)`$ に写すので、$`h\gamma=3`$ の一ステップは $`\mathrm{diag}(0,1)`$ を
$`\mathrm{diag}(1-T(-3),T(-3))=\mathrm{diag}(-\frac38,\frac{11}8)`$ に写す。*数値的観察*：この実行では
$`\rho(t)`$ の最小固有値は $`\approx-2\cdot10^{-8}`$ である（出力される値。刻み幅は高々 $`0.02`$）。`lindblad_channel` は厳密であるが、
セクター内のジャンプしか扱わない。最後に、セクター間のランダムな写像について命題 8.8(2) を、
粒子損失について命題 8.8(3) を、numpy で組み立てた $`M_6`$ 上の GKSL 生成子と照合して確かめる：

```python
Phi = channels.random_inter_sector_channel(fock, fock, kraus_rank=2)
rho = states.random_density_matrix(fock)
on_M = dynamics.lindblad_rhs(None, [Phi], algebra=fock)(0.0, rho.matrix)[0].numpy()
ks, K = fock.k_factors, Phi.kraus[0].numpy()                       # K[l, d, c] = K_l^{dc}, zero padded
embed = lambda Y: sla.block_diag(*[Y[c, :k, :k] for c, k in enumerate(ks)])    # (+)_c Y_c in M_6
Ls = [np.block([[K[l, d, c, :kd, :kc] for c, kc in enumerate(ks)] for d, kd in enumerate(ks)])
      for l in range(2)]                                           # L_l = sum_{d,c} J_l^{dc}
X = embed(rho.matrix[0].numpy())
full = sum(L @ X @ L.conj().T - 0.5 * (L.conj().T @ L @ X + X @ L.conj().T @ L) for L in Ls)   # GKSL on M_6
E_full = full * embed(np.ones((3, 3, 3)))                          # E(X) = sum_d P_d X P_d
assert np.abs(E_full - embed(on_M)).max() < 1e-12                  # Proposition 8.8(2)
print(f"coherences between sectors created on M_6: {np.abs(full - E_full).max():.2f}")

zero = lambda d, c: np.zeros((ks[d], ks[c]))
L_loss = [np.block([[blocks(a)[(d, c)].numpy() if (d, c) in a else zero(d, c) for c in range(3)] for d in range(3)])
          for a in (a1, a2)]                                       # a_1, a_2 on M_6: one block per block row and column
full_loss = sum(L @ X @ L.T - 0.5 * (L.T @ L @ X + X @ L.T @ L) for L in L_loss)
loss_on_M = dynamics.lindblad_rhs(None, [loss], algebra=fock)(0.0, rho.matrix)[0].numpy()
assert np.abs(full_loss - embed(loss_on_M)).max() < 1e-12         # Proposition 8.8(3): no compression needed
```

<a id="s-8-4"></a>

## 8.4 事後選択とジャンプなしの時間発展

*講義ノート：* [8.4節](../../lecture_notes/ja/08_dynamics.md#s-8-4)；[命題 8.10](../../lecture_notes/ja/08_dynamics.md#r-8-10)（ジャンプなしの確率 $`w_0`$）、[命題 8.11](../../lecture_notes/ja/08_dynamics.md#r-8-11)（「粒子が一つも失われなかった」への事後選択）。

ライブラリの `states.lueders_update(rho, P)` は $`P\rho P/\mathrm{Tr}(P\rho P)`$ と
$`\mathrm{Tr}(\rho P)`$ を返す。今度は二つのモードが異なるレート $`\gamma_1>\gamma_2`$ で漏れ出すとする。すると
$`A=\gamma_1n_1+\gamma_2n_2`$ は中心に属さず（セクター $`N=2`$ では、基底 $`|20\rangle,|11\rangle,|02\rangle`$ において
$`\mathrm{diag}(2\gamma_1,\gamma_1+\gamma_2,2\gamma_2)`$ である）、命題 8.9(2) は適用できず、$`N=2`$ に事後選択された状態は、
反エルミート部分 $`-\frac i2A`$ がそのセクターで $`1`$ の定数倍でない $`H_{\mathrm{eff}}`$ による規格化された時間発展である。命題 8.11 を $`n=0`$ について確かめ、さらに
$`S_1`$ の積分をシンプソン則で計算して $`n=1`$ についても確かめる：

```python
g1, g2 = 0.6, 0.2                                                  # mode 1 leaks faster
loss1 = InterSectorChannel.from_blocks(fock, fock, {dc: [K] for dc, K in blocks(a1).items()})
loss2 = InterSectorChannel.from_blocks(fock, fock, {dc: [K] for dc, K in blocks(a2).items()})
m1 = torch.zeros(1, 3, 3, 3, dtype=c128)
m1[0, 2, 0, 0] = 1.0                                               # |20>: both particles in mode 1
ts = [0.0, 0.5, 1.0, 2.0]
traj = dynamics.lindblad_evolve(DensityMatrix(fock, matrix=m1), H_bh, [loss1, loss2], ts, rates=[g1, g2],
                                substeps=100)

n_1 = [np.diag(x) for x in ([0.0], [1.0, 0.0], [2.0, 1.0, 0.0])]   # n_1 in each sector; n_2 = N - n_1
H_eff = [H_bh.matrix[0, N, :N + 1, :N + 1].numpy() - 0.5j * (g1 * n + g2 * (N * np.eye(N + 1) - n))
         for N, n in enumerate(n_1)]

def S0(s, N, X):                                  # e^{-isH_eff} X e^{isH_eff^*} in sector N, s a time or an array
    E = torch.linalg.matrix_exp(torch.from_numpy(-1j * np.asarray(s)[..., None, None] * H_eff[N])).numpy()
    return E @ X @ E.conj().swapaxes(-1, -2)

rho_20 = m1[0, 2].numpy()
for s, rho_s in zip(ts, traj):
    w0 = np.trace(S0(s, 2, rho_20)).real                           # probability of no loss
    post, prob = states.lueders_update(rho_s, fock.central([0.0, 0.0, 1.0]))   # post-select N = 2
    assert abs(prob.item() - w0) < 1e-8                            # Proposition 8.11, n = 0
    assert np.abs(post.matrix[0, 2].numpy() - S0(s, 2, rho_20) / w0).max() < 1e-7
closed = dynamics.von_neumann(DensityMatrix(fock, matrix=m1), H_bh, [2.0])[0]   # no losses at all
dist = 0.5 * np.abs(np.linalg.eigvalsh(S0(2.0, 2, rho_20) / w0 - closed.matrix[0, 2].numpy())).sum()
print(f"t = 2: P(no loss) = {w0:.3f} (without hopping: {math.exp(-2 * g1 * 2):.3f}); "
      f"trace distance post-selected/closed system: {dist:.3f}")

K1, K2 = blocks(a1)[(1, 2)].numpy(), blocks(a2)[(1, 2)].numpy()   # S_1(2) rho = int_0^2 S_0(2-s) Phi S_0(s) rho ds
jump = lambda X: g1 * K1 @ X @ K1.T + g2 * K2 @ X @ K2.T
s_grid = np.linspace(0.0, 2.0, 201)
S1 = scipy.integrate.simpson(S0(2.0 - s_grid, 1, jump(S0(s_grid, 2, rho_20))), x=s_grid, axis=0)
assert np.abs(S1 - traj[-1].matrix[0, 1, :2, :2].numpy()).max() < 1e-8     # Proposition 8.11, n = 1
```

*数値的観察*：$`t=2`$ において、セクター $`N=2`$ の重み $`p_2(2)=w_0(2)`$（「損失なし」）は $`0.230`$ である（上で出力）。
ホッピングがなければ（$`J=0`$）、それはちょうど $`e^{-2\gamma_1t}=e^{-2.4}=0.091`$ となるはずである。そのとき $`|20\rangle`$ は
固有値 $`U-i\gamma_1`$ をもつ $`H_{\mathrm{eff}}`$ の固有ベクトルであり、$`w_0(t)=|e^{-it(U-i\gamma_1)}|^2`$ だからである（命題 8.10）。したがって
ここではセクターの重みがハミルトニアンに依存している。事後選択された状態は、閉じた系の状態から
トレース距離 $`0.135`$（上で出力）だけ離れている：$`N=2`$ という条件づけは、ユニタリな時間発展を超えて状態を変える。

<a id="s-8-5"></a>

## 8.5 定常状態と緩和

*講義ノート：* [8.5節](../../lecture_notes/ja/08_dynamics.md#s-8-5)；[系 8.13](../../lecture_notes/ja/08_dynamics.md#r-8-13)、[命題 8.14](../../lecture_notes/ja/08_dynamics.md#r-8-14)（自然放出）。

[命題 8.14](../../lecture_notes/ja/08_dynamics.md#r-8-14) により、$`M_2`$ 上の自然放出（$`H=\omega|e\rangle\langle e|`$、レート $`\gamma`$ の一つのジャンプ
$`|g\rangle\langle e|`$）は $`e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}}`$、スペクトル
$`\{0,-\gamma,-\gamma/2\pm i\omega\}`$、スペクトルギャップ $`\gamma/2`$ をもつ。スペクトルを
`lindblad_superoperator` で、チャネルとトレース距離を `lindblad_channel` で確かめ、さらに
$`M_2\oplus M_3`$ のセクター内のランダムなジャンプについて[系 8.13](../../lecture_notes/ja/08_dynamics.md#r-8-13)を確かめる：

```python
q2 = TypeIAlgebra([2], [2], precision='double')                    # basis |g>, |e>
omega, gamma = 1.3, 0.4
H_at = q2.from_blocks([torch.tensor([[0.0, 0.0], [0.0, omega]], dtype=f64)])
sigma_minus = q2.from_blocks([torch.tensor([[0.0, 1.0], [0.0, 0.0]], dtype=f64)])   # |g><e|
gen = dynamics.lindblad_superoperator(q2, H_at, [sigma_minus], rates=[gamma])[0][0]
spec = sorted(torch.linalg.eigvals(gen).tolist(), key=lambda z: (z.real, z.imag))
print([f"{z.real:+.2f}{z.imag:+.2f}i" for z in spec])
assert np.allclose(spec, [-gamma, -gamma / 2 - 1j * omega, -gamma / 2 + 1j * omega, 0.0], atol=1e-12)

psi_at = torch.tensor([0.6, 0.8], dtype=c128)                      # rho_ee(0) = 0.64, rho_ge(0) = 0.48
rho_at = DensityMatrix(q2, matrix=torch.outer(psi_at, psi_at.conj())[None, None])
ground = DensityMatrix(q2, matrix=torch.tensor([[[[1.0, 0.0], [0.0, 0.0]]]], dtype=c128))
for s in (0.5, 2.0, 10.0):
    ch = dynamics.lindblad_channel(q2, H_at, [sigma_minus], t=s, rates=[gamma])
    UA = channels.unitary_channel(dynamics.propagator(H_at, s)).compose(
        channels.amplitude_damping_channel(q2, 1 - math.exp(-gamma * s)))
    assert (ch.superoperator()[0] - UA.superoperator()[0]).abs().max() < 1e-12
    r = ch(rho_at).matrix[0, 0]
    assert abs(r[1, 1].real - 0.64 * math.exp(-gamma * s)) < 1e-12
    assert abs(r[0, 1] - 0.48 * np.exp((1j * omega - gamma / 2) * s)) < 1e-12
    exact = math.sqrt(0.64 ** 2 * math.exp(-2 * gamma * s) + 0.48 ** 2 * math.exp(-gamma * s))
    assert abs(states.trace_distance(ch(rho_at), ground).item() - exact) < 1e-12

two = TypeIAlgebra([2, 3], [2, 3], precision='double')
H_two = two.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
L_two = [two.operator(torch.randn(1, 2, 3, 3, dtype=c128)) for _ in range(2)]
kernel = [int((torch.linalg.eigvals(G_c[0]).abs() < 1e-10).sum())
          for G_c in dynamics.lindblad_superoperator(two, H_two, L_two)]
print("stationary states per sector:", kernel)
assert all(dim >= 1 for dim in kernel)                             # Corollary 8.13
```

*数値的観察*：これらのランダムなジャンプについて、コードはセクターごとにちょうど一つの定常状態を見つける。
セクターごとに適用した定理 D.6（[講義ノート、8.5節](../../lecture_notes/ja/08_dynamics.md#s-8-5)）により、各セクターが忠実な定常状態をもち、
$`H_c`$ およびすべての $`L_{j,c}`$、$`L_{j,c}^*`$ と可換な行列が $`1`$ の定数倍だけならば、このことが成り立つ。
コードはこの二つの仮定を（浮動小数点で、上で用いたレート $`g_j=1`$ について）確かめ、それから結論を assert する。

```python
for c, k in enumerate(two.k_factors):                              # Theorem D.6 in every sector
    G_c = dynamics.lindblad_superoperator(two, H_two, L_two)[c][0].numpy()   # on row-major vec(rho_c)
    w_c, V_c = np.linalg.eig(G_c)
    sigma_c = V_c[:, np.argmin(abs(w_c))].reshape(k, k)
    sigma_c = sigma_c / np.trace(sigma_c)
    assert np.abs(G_c @ sigma_c.reshape(-1)).max() < 1e-12            # stationary
    assert np.abs(sigma_c - sigma_c.conj().T).max() < 1e-10 and np.linalg.eigvalsh(sigma_c).min() > 1e-3   # faithful
    ops = [H_two.matrix[0, c, :k, :k].numpy()] + [L.matrix[0, c, :k, :k].numpy() for L in L_two]
    ops += [A.conj().T for A in ops[1:]]
    comm = np.concatenate([np.kron(A, np.eye(k)) - np.kron(np.eye(k), A.T) for A in ops])   # X -> [A, X]
    assert int((np.linalg.svd(comm, compute_uv=False) < 1e-9).sum()) == 1   # commutant = C 1
assert kernel == [1, 1]                                            # the conclusion of Theorem D.6
```

<a id="s-8-6"></a>

## 8.6 比較の対象

QuTiP（`sesolve`、適応的 ODE ソルバーによるリンドブラッド方程式用の `mesolve`、量子軌跡用の `mcsolve`、
`steadystate`、`liouvillian`）、dynamiqs（JAX。GPU 上でバッチ化され微分可能な
シュレーディンガー・ソルバーとリンドブラッド・ソルバー）、QuantumOptics.jl、および SciPy（`scipy.linalg.expm`、および
疎なベクトル化された生成子に対して $`e^{t\mathcal L}\rho`$ を計算する `scipy.sparse.linalg.expm_multiply`）は、
同じ方程式を一つのヒルベルト空間上で解く。上の問題に対しては、それらは同じ数値を与えなければならない。
$`\mathcal M`$ 上では、それらに $`D\times D`$ のブロック対角行列を渡すことになる。このときベクトル化された生成子は
$`D^2`$ 次元の空間に作用し（疎ではあるが、一回の求解の中でセクターどうしが結合される）、ユーザーが
問題を手作業でセクターに分割しない限りそうなる。ライブラリはブロック構造を組み込んでいる：セクターごとの
厳密なプロパゲーターとリンドブラッド・チャネル（$`k_c^2\times k_c^2`$ の超演算子）、命題 8.8 の圧縮を伴う
`InterSectorChannel` としてのセクター間のジャンプ、ハミルトニアンとジャンプ
作用素のバッチ（不規則性の平均のため）、そして自動微分（最後の二つは dynamiqs にもある）。
適応的な刻み幅制御はなく（RK4 は固定刻みで、近似的にしか正値性を保たない）、軌跡
ソルバーもなく、`lindblad_channel` はセクター間のジャンプを受け付けない。それらについては、$`M_D`$ 上の QuTiP や dynamiqs
が基準となる。

<a id="x-8"></a>

## 演習問題

講義ノートの[演習問題 2 と 4](../../lecture_notes/ja/08_dynamics.md#x-8)：8.4節のセクター $`N=1`$ における $`H_{\mathrm{eff}}`$ の固有値
$`-\frac i4(\gamma_1+\gamma_2)\pm\sqrt{J^2-\delta^2}`$、$`\delta=(\gamma_1-\gamma_2)/4`$ と、$`J=|\delta|`$ でのサイズ $`2`$ のただ一つのジョルダン細胞；
および $`M_2`$ 上の純粋位相緩和（$`L=\sigma_z`$）のスペクトル $`\{0,0,-2\gamma\pm i\omega\}`$。両方を数値的に確かめよ。

<details>
<summary>解答</summary>

```python
lam = np.linalg.eigvals(H_eff[1])
delta = (g1 - g2) / 4
root = np.sqrt(J ** 2 - delta ** 2 + 0j)
assert np.allclose(sorted(lam, key=np.real), sorted(-0.25j * (g1 + g2) + np.array([root, -root]), key=np.real))
H_ep = np.array([[-0.5j * g1, delta], [delta, -0.5j * g2]])         # J = delta: the exceptional point
Y = H_ep + 0.25j * (g1 + g2) * np.eye(2)                           # H_ep minus its only eigenvalue
assert np.abs(Y @ Y).max() < 1e-15 and np.linalg.matrix_rank(Y) == 1   # one Jordan block of size 2
```

```python
deph = dynamics.lindblad_superoperator(q2, H_at, [q2.from_blocks([torch.diag(torch.tensor([1.0, -1.0], dtype=f64))])],
                                       rates=[gamma])[0][0]
ev = torch.linalg.eigvals(deph)
assert int((ev.abs() < 1e-12).sum()) == 2
assert np.allclose(sorted(ev.tolist(), key=lambda z: (z.real, z.imag))[:2],
                   [-2 * gamma - 1j * omega, -2 * gamma + 1j * omega], atol=1e-12)
```

</details>

**次へ：** [チュートリアル 9. スピン鎖](09_spin_chains.md)。
