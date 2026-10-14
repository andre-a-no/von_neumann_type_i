[English](../05_random_matrices.md) · [Español](../es/05_random_matrices.md) · [Français](../fr/05_random_matrices.md) · [Deutsch](../de/05_random_matrices.md) · [中文](../zh/05_random_matrices.md) · **日本語** · [Русский](../ru/05_random_matrices.md)

# チュートリアル 5. 環の中のランダム行列

このチュートリアルは講義ノートの[第5章](../../lecture_notes/ja/05_random_matrices.md)に対応する。同章では $\mathcal M=\bigoplus_cM_{k_c}$ のユニタリ群上の
ハール測度を構成し（[系 5.4](../../lecture_notes/ja/05_random_matrices.md#r-5-4)）、その1次と2次のモーメントを計算し（[命題 5.7](../../lecture_notes/ja/05_random_matrices.md#r-5-7)、[定理 5.10](../../lecture_notes/ja/05_random_matrices.md#r-5-10)）、
準位統計をセクターごとに計算しなければならない理由を示す（[定理 5.15](../../lecture_notes/ja/05_random_matrices.md#r-5-15)）。ここではライブラリでハール・ユニタリと
ガウス・アンサンブルを生成し、これらの結果をモンテカルロ法で確かめる：平均は標準誤差の $5$ 倍以内で厳密な値と一致しなければならない。
定理の裏付けのない結果は数値的観察と明記する。

本章を通じて、$\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$、$D=\sum_ck_c$ とし、$P_c$ は
極小中心射影、$E$ は命題 2.8 の中心への条件付き期待値とする。$U(k)$ は $k\times k$ ユニタリ行列のなす群であり、
$U(\mathcal M)=\{U\in\mathcal M:U^*U=1\}$ は $\mathcal M$ のユニタリ群である。確率変数は期待値 $\mathbb E$ をもつ確率空間上に定義され、ランダム行列の期待値は成分ごとにとる。*標準複素ガウス* 変数とは、実部と虚部が独立で
$\sim N(0,\tfrac12)$ であり、したがって $\mathbb E|z|^2=1$ となるものである（ドキュメントによれば、複素
dtype に対して `torch.randn` が生成するのはこれである。最初のコードブロックで2次のモーメントを確かめる）。*ジニブル行列* とは、成分が独立な標準複素ガウス変数である行列である。

<a id="s-5-1"></a>

## 5.1 ハール測度

*講義ノート：* [5.1節](../../lecture_notes/ja/05_random_matrices.md#s-5-1)；[定理 5.3](../../lecture_notes/ja/05_random_matrices.md#r-5-3)（ジニブル行列の QR 分解のユニタリ因子は、位相を $\Lambda=\mathrm{diag}(R'_{ii}/|R'_{ii}|)$ で固定すればハール分布に従う）、[系 5.4](../../lecture_notes/ja/05_random_matrices.md#r-5-4)（$U(\mathcal M)$ のハール・ユニタリのブロックは独立で、それぞれハール分布に従う）、[系 5.5](../../lecture_notes/ja/05_random_matrices.md#r-5-5)（$\mathbb EU_{ij}=0$、$\mathbb E|U_{ij}|^2=1/k$）。$SU(k)$ のハール測度は $e^{-i\theta/k}U$（$\theta=\arg\det U$）の分布である。

[定理 5.3](../../lecture_notes/ja/05_random_matrices.md#r-5-3) の補正 $\Lambda$ がなければ、$Q'$ の分布は QR ルーチンの符号の規約に依存する。
`torch.linalg.qr` が呼び出すルーチン（LAPACK）について、下のコードは
$\mathbb E\mathrm{Re}\,\mathrm{Tr}\,Q'$ が $0$ から標準誤差の何倍も離れていることを見いだすが、ハール・ユニタリでは
$\mathbb EU=0$ である（系 5.5）。したがってこの $Q'$ はハール分布に従わない。これはそのルーチンについての主張であり、
コードで確かめられるもので、定理ではない。

ライブラリは `random_unitary_operator` により $U(\mathcal M)$ 上のハール・ユニタリを生成する。ブロックサイズごとに一回のバッチ化された QR
分解で、すべてのセクターとサンプルを一度に処理する（GPU 上ではグラム・シュミットを用いるが、
補題 5.2 の一意性により同じ $Q$ が得られる）。測度 `'haar'` は各ブロックにおける定理 5.3 であり、
複素環に対しては `'haar_su'` が上のように $e^{-i\theta/k}$ を掛ける。`'coe'` は演習問題 3 のアンサンブル
$S=W^TW$ である（ライブラリのその他の測度は本講義では用いない）。系 5.4 と 5.5、
およびライブラリのサンプラーがちょうど定理 5.3 の補正付き QR であることを確かめる：

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, SpinChain, channels

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def mc(samples, exact, name):
    """Monte Carlo mean of real samples with its standard error; must agree within 5 standard errors."""
    mean, err = samples.mean().item(), samples.std().item() / math.sqrt(samples.numel())
    print(f"{name}: {mean:.4f} +- {err:.4f}   (theory {exact:.4f})")
    assert abs(mean - exact) <= 5 * err + 1e-12
    return mean

alg = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')        # M_1 (+) M_2 (+) M_3, D = 6
B = 100_000
U = alg.random_unitary_operator(B)                                   # Haar on U(1) x U(2) x U(3)
assert close(U.adjoint() @ U, alg.identity())                        # unitary elements of the algebra
u = [U.matrix[:, c, :k, :k] for c, k in enumerate(alg.k_factors)]   # the blocks U_c
for c, k in enumerate(alg.k_factors):                                # Corollary 5.5
    mc(u[c][:, 0, 0].real, 0.0, f"E Re (U_{c})_11")
    mc(u[c][:, 0, 0].abs() ** 2, 1 / k, f"E |(U_{c})_11|^2")
mc(u[1][:, 0, 0].abs() ** 2 * u[2][:, 0, 0].abs() ** 2, 1 / 6, "E |(U_1)_11|^2 |(U_2)_11|^2")   # independent blocks

z = torch.randn(B, dtype=dt)                                         # torch.randn draws standard complex Gaussians:
mc(z.real ** 2, 0.5, "E (Re z)^2")                                   # Re z, Im z ~ N(0, 1/2), uncorrelated
mc(z.imag ** 2, 0.5, "E (Im z)^2")
mc(z.real * z.imag, 0.0, "E Re z Im z")

SU = alg.random_unitary_operator(4, measure='haar_su').matrix
assert torch.allclose(torch.linalg.det(SU[:, 2, :3, :3]), torch.ones(4, dtype=dt))

torch.manual_seed(1)
Z = torch.randn(B, 2, 2, dtype=dt)                                   # Ginibre matrices
torch.manual_seed(1)
V = alg.random_unitary(2, batch_size=B)                              # the library, same random numbers
Q, R = torch.linalg.qr(Z)                                            # LAPACK: diag(R) real, either sign
d = R.diagonal(dim1=-2, dim2=-1)
assert close(V, Q * (d / d.abs()).unsqueeze(-2))                     # Theorem 5.3: V = Q Lambda
raw = Q.diagonal(dim1=-2, dim2=-1).sum(-1).real                      # Re Tr Q without the correction
print(f"E Re Tr Q without the phase correction: {raw.mean():.4f} +- {raw.std() / math.sqrt(B):.4f}")
assert abs(raw.mean()) > 20 * raw.std() / math.sqrt(B)               # far from 0: not Haar
```

*数値的観察.* 補正がなければ、コードは $k=2$ に対し
$\mathbb E\mathrm{Re}\,\mathrm{Tr}\,Q'\approx-0.84$ を見いだすが、ハール・ユニタリでは $\mathbb EU=0$ である。

<a id="s-5-2"></a>

## 5.2 モーメントとワインガルテン計算

*講義ノート：* [5.2節](../../lecture_notes/ja/05_random_matrices.md#s-5-2)；[命題 5.7](../../lecture_notes/ja/05_random_matrices.md#r-5-7)（$\mathbb E[UXU^*]=E(X)$）、[補題 5.8](../../lecture_notes/ja/05_random_matrices.md#r-5-8)（すべての $V\otimes V$ の可換子環は $\mathrm{span}\{1,F\}$、$F$ はスワップ）、[定理 5.10](../../lecture_notes/ja/05_random_matrices.md#r-5-10)（2次のワインガルテン公式）、[系 5.11](../../lecture_notes/ja/05_random_matrices.md#r-5-11)。

命題 5.7 を `channels.center_expectation` と照合し、補題 5.8 を線形代数で確かめ（サンプルした三つの
$V\otimes V$ の可換子環を零空間として計算すると、それは $F$ を含み $2$ 次元であるから、
補題 5.8 によりすでに $\mathrm{span}\{1,F\}$ である）、定理 5.10 と系 5.11 を上のサンプルの
$U(3)$ ブロックを用いたモンテカルロ法で確かめる：

```python
X = alg.from_blocks([torch.randn(k, k, dtype=dt) for k in alg.k_factors])      # a fixed element
T = (U @ X @ U.adjoint()).matrix                                    # B samples of U X U*
mean, err = T.mean(0), T.std(0) / math.sqrt(B)
EX = channels.center_expectation(alg)(X).matrix[0]
random = err > 1e-12                                                 # entries that are not constant
z = ((mean - EX).abs()[random] / err[random]).max().item()
print(f"twirl over U(M) vs E(X): largest deviation {z:.2f} standard errors")
assert ((mean - EX).abs() <= 5 * err + 1e-12).all()                 # Proposition 5.7

def commutant_dim(ops, tol=1e-9):
    """Dimension of {X : XA = AX for all A in ops}, as the null space of X -> (AX - XA)_A."""
    d = ops[0].shape[0]
    eye = torch.eye(d, dtype=dt)
    ops = [A.contiguous() for A in ops]                              # torch.kron needs contiguous input
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    return int((torch.linalg.svdvals(L) < tol).sum())

k = 3
F = torch.zeros(k * k, k * k, dtype=dt)
for i in range(k):
    for j in range(k):
        F[j * k + i, i * k + j] = 1                                  # F (e_i (x) e_j) = e_j (x) e_i
VV = [torch.kron(W, W) for W in alg.random_unitary(k, batch_size=3).contiguous()]
assert all(close(F @ W, W @ F) for W in VV)
assert commutant_dim(VV) == 2                                        # Lemma 5.8: span{1, F}

u3 = u[2]                                                            # Haar on U(3)
mc(u3[:, 0, 0].abs() ** 4, 2 / (k * (k + 1)), "E|U_11|^4")
mc(u3[:, 0, 0].abs() ** 2 * u3[:, 0, 1].abs() ** 2, 1 / (k * (k + 1)), "E|U_11|^2 |U_12|^2")
mc((u3[:, 0, 0] * u3[:, 1, 1] * (u3[:, 0, 1] * u3[:, 1, 0]).conj()).real,
   -1 / ((k - 1) * k * (k + 1)), "E U_11 U_22 conj(U_12 U_21)")

n_s = 20_000                                                         # Theorem 5.10 for a random X
X2 = torch.randn(k * k, k * k, dtype=dt)
trX, trFX = X2.diagonal().sum(), (F @ X2).diagonal().sum()
alpha, beta = (trX - trFX / k) / (k ** 2 - 1), (trFX - trX / k) / (k ** 2 - 1)
UU = torch.einsum('bij,bkl->bikjl', u3[:n_s], u3[:n_s]).reshape(n_s, k * k, k * k)   # U (x) U
T2 = UU @ X2 @ UU.mH
mean, err = T2.mean(0), T2.std(0) / math.sqrt(n_s)
assert ((mean - alpha * torch.eye(k * k, dtype=dt) - beta * F).abs() <= 5 * err).all()
```

<a id="s-5-3"></a>

## 5.3 準位統計

*講義ノート：* [5.3節](../../lecture_notes/ja/05_random_matrices.md#s-5-3)：GOE $H=(G+G^T)/2$、GUE $H=(G+G^*)/\sqrt2$、および準位間隔 $s_i=\lambda_{i+1}-\lambda_i$ の間隔比 $r_i=\min(s_i,s_{i+1})/\max(s_i,s_{i+1})$；[命題 5.13](../../lecture_notes/ja/05_random_matrices.md#r-5-13)（$2\times2$ 行列に対するウィグナー推測）と[命題 5.14](../../lecture_notes/ja/05_random_matrices.md#r-5-14)（平均比は独立な準位では $2\ln2-1$、$3\times3$ の GOE と GUE ではそれぞれ $4-2\sqrt3$ と $\frac{2\sqrt3}\pi-\frac12$。ここで $\mathbb Er=I_1/I_0$、$I_j=\int_0^1x^jf_\beta(x)\,dx$）。

関数 `SpinChain.level_spacing_ratio`（`chains` 内）は準位の各行の平均間隔比を返す。
命題 5.13 と 5.14 を確かめ、$n=64$ に対する GOE の値を記録する：

```python
ratio = SpinChain.level_spacing_ratio
n2 = 200_000
alg_r = TypeIAlgebra([2], [2], complex_valued=False, precision='double')     # real symmetric 2 x 2
alg_c = TypeIAlgebra([2], [2], precision='double')
G = torch.randn(n2, 1, 2, 2, dtype=torch.float64)
Gc = torch.randn(n2, 1, 2, 2, dtype=dt)
goe2 = alg_r.operator((G + G.mT) / 2, is_self_adjoint=True)
gue2 = alg_c.operator((Gc + Gc.mH) / math.sqrt(2), is_self_adjoint=True)
cdf = {1: lambda t: 1 - math.exp(-math.pi * t ** 2 / 4),
       2: lambda t: math.erf(2 * t / math.sqrt(math.pi)) - 4 * t / math.pi * math.exp(-4 * t ** 2 / math.pi)}
for beta, H, mean_s in ((1, goe2, math.sqrt(math.pi)), (2, gue2, 4 / math.sqrt(math.pi))):
    lam = H.eigenvalues()[0]                                         # (n2, 2), ascending
    S = (lam[:, 1] - lam[:, 0]) / mean_s
    mc(S, 1.0, f"beta={beta}: E S")
    for t in (0.25, 1.0, 2.0):
        mc((S <= t).double(), cdf[beta](t), f"beta={beta}: P(S <= {t})")

levels = torch.rand(4000, 50, dtype=torch.float64)                   # independent levels
mc(ratio(levels), 2 * math.log(2) - 1, "<r>, independent uniform levels")      # Proposition 5.14(a)

n3 = 100_000                                                         # Proposition 5.14(b): 3 x 3
G = torch.randn(n3, 1, 3, 3, dtype=torch.float64)
Gc = torch.randn(n3, 1, 3, 3, dtype=dt)
goe3 = TypeIAlgebra([3], [3], complex_valued=False, precision='double').operator((G + G.mT) / 2, is_self_adjoint=True)
gue3 = TypeIAlgebra([3], [3], precision='double').operator((Gc + Gc.mH) / math.sqrt(2), is_self_adjoint=True)
h = 1e-4
xs = (torch.arange(int(1 / h), dtype=torch.float64) + 0.5) * h      # midpoint rule on [0, 1]
for beta, H, exact in ((1, goe3, 4 - 2 * math.sqrt(3)), (2, gue3, 2 * math.sqrt(3) / math.pi - 0.5)):
    f = (xs + xs ** 2) ** beta / (1 + xs + xs ** 2) ** (1 + 1.5 * beta)
    assert abs((xs * f).sum() / f.sum() - exact) < 1e-7              # E r = I_1 / I_0
    r3 = ratio(H.eigenvalues()[0])                                   # one ratio per matrix
    mc(r3, exact, f"beta={beta}: <r>, 3 x 3")
    mc((r3 <= 0.5).double(), (f * (xs <= 0.5)).sum().item() / f.sum().item(), f"beta={beta}: P(r <= 1/2), 3 x 3")

n = 64                                                               # numerical observation, no theorem
alg64 = TypeIAlgebra([n], [n], complex_valued=False, precision='double')
G = torch.randn(800, 1, n, n, dtype=torch.float64)
lam = alg64.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()[0]
r_goe = ratio(lam[:, n // 4: 3 * n // 4])                            # 800 samples of the mean ratio
print(f"<r>, GOE(64), middle half of the spectrum: {r_goe.mean():.4f} +- {r_goe.std() / math.sqrt(800):.4f}")

def same_mean(a, b, name):
    """Two independent Monte Carlo samples of the same law: means agree within 5 combined standard errors."""
    ea, eb = a.std().item() / math.sqrt(a.numel()), b.std().item() / math.sqrt(b.numel())
    print(f"{name}: {a.mean():.4f} vs {b.mean():.4f} +- {math.hypot(ea, eb):.4f}")
    assert abs(a.mean().item() - b.mean().item()) <= 5 * math.hypot(ea, eb)
```

*数値的観察.* コードは、$n=64$ の GOE のスペクトルの中央半分について $\langle r\rangle\approx0.53$ を見いだす。これは
$3\times3$ の厳密な値 $4-2\sqrt3\approx0.536$ より小さく、独立な準位の値
$2\ln2-1\approx0.386$ よりはるかに大きい。この数は計算の出力であって
定理ではない。以下ではこれを、同じサイズ $64$ の GOE ブロックに対する参照値として用いる。それらの間隔比は
まったく同じ分布をもつ。

<a id="s-5-4"></a>

## 5.4 なぜセクターを混ぜてはならないか

*講義ノート：* [5.4節](../../lecture_notes/ja/05_random_matrices.md#s-5-4)；[定理 5.15](../../lecture_notes/ja/05_random_matrices.md#r-5-15) により、異なる独立なブロックの準位は反発しない：二つの GOE ブロックでは、合併したスペクトルの最小間隔 $\delta$ について $P(\delta<\varepsilon)$ は少なくとも $\varepsilon$ のオーダーであるが、一つの GOE ブロックの内部では高々 $\varepsilon^2$ のオーダーである。

セクターが多いほど効果は強まる。多数のセクターに対する極限定理は証明しない。以下の
実験は数値的観察である。ライブラリではブロックは最初から分かれている：
`eigenvalues()` はセクターごとに一つのスペクトルを返し、それらを混ぜるには明示的な `torch.cat` が必要である。サイズ $64$ の独立な
GOE ブロック $C$ 個を一つの環に入れ、セクターごとの平均間隔比と合併したスペクトルのそれとを
比較する。各セクターはサイズ $64$ の GOE 行列であるから、その間隔比は上の `r_goe` とまったく同じ分布をもち、
コードはこれを検定する。合併したスペクトルの振る舞いは観察されるだけである：

```python
k, merged_r = 64, {}
for C in (1, 2, 4, 8, 16):
    alg_C = TypeIAlgebra([k] * C, [k] * C, complex_valued=False, precision='double')
    G = torch.randn(800 // C, C, k, k, dtype=torch.float64)
    lam = alg_C.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()    # one spectrum per sector
    per_sector = torch.stack([ratio(l[:, k // 4: 3 * k // 4]) for l in lam])
    merged = torch.cat(lam, dim=-1).sort(dim=-1).values                       # the spectrum on C^D
    r_m = ratio(merged[:, k * C // 4: 3 * k * C // 4])
    merged_r[C] = r_m.mean().item()
    print(f"C = {C:2d}: <r> per sector {per_sector.mean():.4f}, merged {merged_r[C]:.4f}"
          f" +- {r_m.std() / math.sqrt(len(r_m)):.4f}")
    same_mean(per_sector.flatten(), r_goe, f"C = {C:2d}: sectors vs GOE(64)")   # every sector is GOE(64)
assert merged_r[1] > merged_r[2] > merged_r[4] > merged_r[16]       # numerical observation
assert abs(merged_r[16] - (2 * math.log(2) - 1)) < 0.01              # numerical observation
```

*数値的観察.* どのセクターも $\langle r\rangle\approx0.53$ であり、これは 5.3節の GOE(64) の値である
（この一致はある事実の検定である：各セクターは GOE(64) 行列である）。合併したスペクトルでは、二つのセクターで
$0.42$ に下がり、$8$ と $16$ のセクターでは約 $0.39$ まで下がり、独立な準位の値 $2\ln2-1$ との差は $0.01$ 以内となる
（命題 5.14(a)）。$C$ が大きくなるにつれて合併したスペクトルが独立な準位に近づくことは、コードが観察する事柄であり、
本講義では証明しない。

同じことが物理的なモデルでも起こる。ランダム磁場中のハイゼンベルク鎖（第9章）は
磁化を保存するから、そのハミルトニアンは $\bigoplus_{N=0}^LM_{\binom LN}$ に属する（命題 9.1）。そのブロックは独立では *ない*
——すべてが同じランダム磁場に依存する——ので、定理 5.15 はそのままでは適用できない。
下の実験は、数値的観察として、同じ効果を示す。

```python
chain = SpinChain(10, complex_valued=False)                         # sectors N = 0, ..., 10
lam = chain.random_field_heisenberg(W=0.5, batch_size=20).eigenvalues()
k5 = chain.sector_dim(5)                                            # 252
r_sector = ratio(lam[5][:, k5 // 4: 3 * k5 // 4]).mean().item()
spectrum = torch.cat(lam, dim=-1).sort(dim=-1).values               # all 1024 levels mixed
r_all = ratio(spectrum[:, 256:768]).mean().item()
print(f"Heisenberg chain, L = 10: <r> in sector N = 5: {r_sector:.4f}, all sectors mixed: {r_all:.4f}")
assert r_sector > 0.49 and r_all < 0.42                             # numerical observation
```

*数値的観察.* セクター $N=5$ の内部では、コードは上の GOE(64) の値に近い $\langle r\rangle\approx0.53$ を見いだし、
すべてのセクターを混ぜたスペクトルでは $2\ln2-1$ に近い $0.39$ となる。

<a id="s-5-5"></a>

## 5.5 比較の対象

ハール・ユニタリ：`scipy.stats.unitary_group` と `ortho_group`、QuTiP のハール・ランダムなユニタリ、
そして位相補正を施した任意の `numpy`/`torch` の QR。
ワインガルテン計算：記号計算パッケージ IntU（Mathematica；Z. Puchała, J. A. Miszczak）、RTNI
（Mathematica と Python；M. Fukuda, R. König, I. Nechita, J. Phys. A 52 (2019), 425303）、および
Python パッケージ Haarpy は任意の次数のモーメントを記号的に計算する。定理 5.10 はその
2次の場合であり、定理 5.12（付録 B で証明）はそれらが実装する一般の公式である。準位統計は通常数行の `numpy` で計算され、
QuSpin のような厳密対角化コードはセクターごとにスペクトルを与える。ライブラリが
付け加えるのは、環のユニタリ群を一つの対象として扱うことである：あらゆるサイズの独立なハール・ブロックをバッチ化された呼び出しで生成し、
ツワールを条件付き期待値
`channels.center_expectation` と同一視し、スペクトルをセクターごとに返す。そのため、セクターを混ぜることは決して
既定の動作にならない。

<a id="x-5"></a>

## 演習問題

講義ノートの[演習問題](../../lecture_notes/ja/05_random_matrices.md#x-5)の数値的な部分：
1. $\mathcal M=M_1\oplus M_2\oplus M_3$ とし、$U$ を $U(\mathcal M)$ 上でハール分布に従うとして、$\mathbb E|\mathrm{Tr}\,U|^2=C=3$ と
   $\mathbb E|\mathrm{Tr}\,U|^4=17$ を確かめよ。
2. $k=3$ に対し $3$ 次のワインガルテン行列（[定理 5.12](../../lecture_notes/ja/05_random_matrices.md#r-5-12)）を数値的に計算し、
   $\mathbb E|U_{11}|^6=\frac1{10}$ を確かめよ。
3. **円形直交アンサンブル**（`measure='coe'`）は、$W$ を $U(k)$ 上でハール分布に従うとして $S=W^TW$ で与えられる。ライブラリの
   COE が $W^TW$ であること、$S$ がユニタリかつ対称であること、$m=1,2,3$ に対し $\mathbb E\mathrm{Tr}\,S^m=0$ であること、
   固有位相の平均密度が一様であることを確かめよ。$k=32$ の COE の固有位相について $\langle r\rangle$ を計算し、
   5.3節の GOE の値と比較せよ（数値的観察）。

<details>
<summary>解答</summary>

*演習問題 1.*

```python
tr = U.trace                                                        # Tr U = sum_c Tr U_c, shape (B,)
mc(tr.abs() ** 2, 3.0, "E|Tr U|^2")
mc(tr.abs() ** 4, 17.0, "E|Tr U|^4")
```

*演習問題 2.*

```python
import itertools

def n_cycles(p):                                                    # number of cycles of a permutation
    seen, count = set(), 0
    for i in range(len(p)):
        count += i not in seen
        while i not in seen:
            seen.add(i)
            i = p[i]
    return count

perms = list(itertools.permutations(range(3)))                      # S_3; s.index(j) = s^{-1}(j)
Gram = torch.tensor([[3.0 ** n_cycles([s.index(p[i]) for i in range(3)]) for p in perms] for s in perms],
                    dtype=torch.float64)                            # k^{#cycles(s^{-1} p)}, k = 3
Wg3 = torch.linalg.inv(Gram)                                        # the Weingarten matrix, n = k = 3
assert abs(Wg3.sum().item() - 6 / (3 * 4 * 5)) < 1e-12
mc(u3[:, 0, 0].abs() ** 6, 0.1, "E|U_11|^6")
```

*演習問題 3.*

   *数値的観察.* コードは $k=32$ の COE について $\langle r\rangle\approx0.53$ を見いだし、これは 5.3節の GOE(64) の値と
   $0.015$ 以内で一致する。この一致は計算の出力である。本講義の
   定理で COE の間隔統計を GOE のそれと結びつけるものはない（円形
   アンサンブルについては F. J. Dyson, J. Math. Phys. 3 (1962), 140–156 を参照）。

```python
alg_s = TypeIAlgebra([32], [32], precision='double')
torch.manual_seed(2)
S = alg_s.random_unitary_operator(400, measure='coe').matrix[:, 0]
torch.manual_seed(2)
W = alg_s.random_unitary(32, batch_size=400)                         # same random numbers
assert close(S, W.mT @ W)                                            # the library's COE is W^T W
assert close(S, S.mT) and close(S @ S.mH, torch.eye(32, dtype=dt))
for m in (1, 2, 3):                                                  # E Tr S^m = 0
    trm = torch.linalg.matrix_power(S, m).diagonal(dim1=-2, dim2=-1).sum(-1)
    mc(trm.real, 0.0, f"E Re Tr S^{m}")
    mc(trm.imag, 0.0, f"E Im Tr S^{m}")
phases = torch.linalg.eigvals(S).angle().sort(dim=-1).values         # in (-pi, pi]
mc(((phases > 0) & (phases <= 1)).sum(-1).double(), 32 / (2 * math.pi), "E #eigenphases in (0, 1]")
s = torch.diff(torch.cat([phases, phases[:, :1] + 2 * math.pi], dim=-1), dim=-1)   # around the circle
r = torch.minimum(s, s.roll(-1, -1)) / torch.maximum(s, s.roll(-1, -1))
r_coe = r.mean(-1)
print(f"<r>, COE(32) eigenphases: {r_coe.mean():.4f}; GOE(64), Section 5.3: {r_goe.mean():.4f}")
assert abs(r_coe.mean() - r_goe.mean()) < 0.015                      # numerical observation
```

</details>

**次へ：** [チュートリアル 6](06_quantum_channels.md)：ライブラリにおける量子チャネル。
