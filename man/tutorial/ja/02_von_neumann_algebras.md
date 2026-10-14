[English](../02_von_neumann_algebras.md) · [Español](../es/02_von_neumann_algebras.md) · [Français](../fr/02_von_neumann_algebras.md) · [Deutsch](../de/02_von_neumann_algebras.md) · [中文](../zh/02_von_neumann_algebras.md) · **日本語** · [Русский](../ru/02_von_neumann_algebras.md)

# チュートリアル 2. 有限次元フォン・ノイマン環

このチュートリアルは講義ノートの[第2章](../../lecture_notes/ja/02_von_neumann_algebras.md)に対応する。同章は、$`\mathbb C^n`$ 上の任意のフォン・ノイマン環が
ユニタリな基底変換を除いて
$`\bigoplus_cM_{k_c}(\mathbb C)\otimes1_{m_c}`$ であることを証明する（[定理 2.6](../../lecture_notes/ja/02_von_neumann_algebras.md#r-2-6)）。
ここでは対称性の群の例で可換子環を数値的に計算し、ライブラリがそのような環、その中心、および中心の上への条件付き期待値をどのように格納するかを見る。

<a id="s-2-4"></a>

## 2.4 例：対称性の群とその可換子環

*講義ノート：* [2.4節](../../lecture_notes/ja/02_von_neumann_algebras.md#s-2-4)；[定理 2.2](../../lecture_notes/ja/02_von_neumann_algebras.md#r-2-2)、[定理 2.6](../../lecture_notes/ja/02_von_neumann_algebras.md#r-2-6)、[命題 2.7](../../lecture_notes/ja/02_von_neumann_algebras.md#r-2-7)。

講義ノートの[2.4節](../../lecture_notes/ja/02_von_neumann_algebras.md#s-2-4)は、$`(\mathbb C^2)^{\otimes3}`$ の三つのテンソル因子を置換する
$`S_3`$ の表現 $`\pi`$ を分解する：$`\pi(S_3)'\cong M_4\oplus M_2`$ の次元は
$`20`$、$`\pi(S_3)''\cong M_1\oplus M_2`$ の次元は $`5`$ である。下のコードは二つの次元を数値的に確認し、
さらに $`\pi(S_3)'`$ の可換子環が再び $`\pi(S_3)''`$ であること
（[定理 2.2](../../lecture_notes/ja/02_von_neumann_algebras.md#r-2-2)）も確認する：集合 $`S`$ の可換子環は線形写像 $`X\mapsto(AX-XA)_{A\in S}`$ の
核であり、その次元はここでは浮動小数点演算で、許容誤差より大きい特異値を数えることによって計算される。

```python
import itertools
import torch

torch.manual_seed(0)
dt = torch.complex128

def permutation_operator(perm, n_qubits=3):
    """Unitary on (C^2)^{(x)3} that permutes the tensor factors."""
    dim = 2 ** n_qubits
    P = torch.zeros(dim, dim, dtype=dt)
    for bits in itertools.product([0, 1], repeat=n_qubits):
        src = int(''.join(map(str, bits)), 2)
        dst = int(''.join(str(bits[perm[i]]) for i in range(n_qubits)), 2)
        P[dst, src] = 1
    return P

group = [permutation_operator(p) for p in itertools.permutations(range(3))]

def commutant_basis(ops, tol=1e-10):
    """Orthonormal basis (Hilbert-Schmidt) of {X : X A = A X for all A in ops}."""
    d = ops[0].shape[0]
    eye = torch.eye(d, dtype=dt)
    # vec(AX - XA) = (A (x) 1 - 1 (x) A^T) vec(X) for row-major vec
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    _, s, Vh = torch.linalg.svd(L)
    null = Vh[(s > tol).sum():].conj().resolve_conj()   # rows span the null space
    return [v.reshape(d, d) for v in null]

def span_dim(ops, tol=1e-10):
    return int((torch.linalg.svdvals(torch.stack([A.reshape(-1) for A in ops])) > tol).sum())

comm = commutant_basis(group)
print(len(comm), span_dim(group), len(commutant_basis(comm)))
assert len(comm) == 20                                   # Step 1: pi(S3)'  = M_4 (+) M_2
assert span_dim(group) == 5                              # Step 2: pi(S3)'' = M_1 (+) M_2
assert len(commutant_basis(comm)) == 5                   # double commutant theorem
```

重複度はスペクトルに現れる（[講義ノート、2.4節](../../lecture_notes/ja/02_von_neumann_algebras.md#s-2-4)、
[定理 2.6](../../lecture_notes/ja/02_von_neumann_algebras.md#r-2-6) から）：
任意の自己共役な $`K\in\pi(S_3)'`$ は対称部分空間 $`F`$ の補空間上で偶数の重複度の固有値をもち、
$`Z`$ をジニブル行列、$`E`$ を群平均として $`K=E(Z+Z^*)`$ とおくと、確率 1 で
$`K`$ は $`F`$ 上に $`4`$ 個の単純固有値を、その補空間上に $`2`$ 個の二重縮退した固有値をもつ。コードは
縮退した対（これはすべての $`K`$ について成り立つ）を確かめ、さらにサンプルした $`K`$ について、二つの対と $`F`$ 上の四つの
固有値が相異なることを、$`F`$ への直交射影 $`P_F=\frac1{|G|}\sum_g\pi(g)`$ を用いて確かめる。

```python
def group_average(X):
    return sum(g @ X @ g.conj().T for g in group) / len(group)

H0 = torch.randn(8, 8, dtype=dt)                        # Ginibre: Re, Im independent N(0, 1/2)
Hsym = group_average(H0 + H0.conj().T)                  # K = E(Z + Z^*), a random symmetric Hamiltonian
assert all(torch.allclose(g @ Hsym, Hsym @ g) for g in group)

P_sym = sum(group) / len(group)                         # P_F, the projection onto the symmetric subspace
w, Q = torch.linalg.eigh(P_sym)
Q_rest = Q[:, w.real < 0.5]                             # orthonormal basis of the complement (dim 4)
e_rest = torch.linalg.eigvalsh(Q_rest.conj().T @ Hsym @ Q_rest)
print(e_rest)
assert Q_rest.shape[1] == 4
assert torch.allclose(e_rest[0], e_rest[1]) and torch.allclose(e_rest[2], e_rest[3])   # even multiplicities
e_sym = torch.linalg.eigvalsh(Q[:, w.real > 0.5].conj().T @ Hsym @ Q[:, w.real > 0.5])   # eigenvalues on F
assert e_sym.shape[0] == 4 and e_rest[2] - e_rest[1] > 1e-6                           # this sample: generic
assert (e_sym[1:] - e_sym[:-1]).min() > 1e-6                                           # 4 simple eigenvalues on F
assert (e_sym[:, None] - e_rest[None, :]).abs().min() > 1e-6                           # all six distinct
```

<a id="s-2-5"></a>

## 2.5 ライブラリにおける環

*講義ノート：* [2.3節](../../lecture_notes/ja/02_von_neumann_algebras.md#s-2-3)、[定理 2.6](../../lecture_notes/ja/02_von_neumann_algebras.md#r-2-6)；中心：[2.2節](../../lecture_notes/ja/02_von_neumann_algebras.md#s-2-2)、[補題 2.3](../../lecture_notes/ja/02_von_neumann_algebras.md#r-2-3)。

ライブラリは定理 2.6 の抽象的な環 $`\bigoplus_cM_{k_c}(\mathbb C)`$ を
$`\bigoplus_c\mathbb C^{k_c}`$ 上で、重複度なしに格納する：作用素はそのブロック $`A_c`$ のリストであり、共通のサイズにパディングされる。上の例の可換子環は
`TypeIAlgebra([4, 2], [4, 2])` である。

```python
from torch_vn_algebra import TypeIAlgebra, channels

alg = TypeIAlgebra([4, 2], [4, 2], precision='double')    # M_4 (+) M_2
print(alg.k_factors, alg.C)

A = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
print(A.matrix.shape)              # (1, 2, 4, 4): block 1 occupies the top-left 2 x 2 corner
assert A.matrix[0, 1, 2:, :].abs().max() == 0 and A.matrix[0, 1, :, 2:].abs().max() == 0
```

中心は元 $`\sum_cz_cP_c`$ からなる。`alg.central` がそれらを構成し、
それらはすべての元と可換である。

```python
def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

P0, P1 = alg.central([1.0, 0.0]), alg.central([0.0, 1.0])   # minimal central projections
Z = alg.central([2.0, -1.0])
assert P0.is_projection and P1.is_projection
assert close(P0 + P1, alg.identity())
assert close(Z @ A, A @ Z)
assert close(A @ P0 + A @ P1, A)                         # A = sum_c A P_c
```

<a id="s-2-6"></a>

## 2.6 条件付き期待値

*講義ノート：* [2.6節](../../lecture_notes/ja/02_von_neumann_algebras.md#s-2-6)。

[講義ノート、命題 2.8](../../lecture_notes/ja/02_von_neumann_algebras.md#r-2-8) により、$`E(A)=\sum_c\frac{\mathrm{Tr}\,A_c}{k_c}\,P_c`$
は $`\bigoplus_cM_{k_c}`$ からその中心の上への条件付き期待値である：これは正であり、中心上で恒等写像であり、
中心上の双加群写像であり、すべてのトレース $`\sum_cw_c\mathrm{Tr}\,A_c`$ を保つ。

ライブラリではこの写像はチャネル `channels.center_expectation(alg)` である。これを作用素に適用すると古典的な（セクター平均した）部分が得られる。

```python
E = channels.center_expectation(alg)
EA = E(A)
assert close(E(EA), EA)                                  # idempotent
assert close(E(Z @ A @ P1), Z @ EA @ P1)                 # bimodule property
A2 = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
assert close(EA @ A2, A2 @ EA)                           # E(A) commutes with a fresh element
for c, k in enumerate(alg.k_factors):                    # E(A)_c = (Tr A_c / k_c) 1, so Tr E(A)_c = Tr A_c
    trA = A.matrix[0, c].diagonal().sum()
    assert close(EA.matrix[0, c, :k, :k], trA / k * torch.eye(k, dtype=dt))
    assert abs(EA.matrix[0, c].diagonal().sum() - trA) < 1e-10
```

<a id="s-2-7"></a>

## 2.7 比較の対象

2.4節の例は任意の線形代数パッケージで確かめることができる。
表現論のソフトウェア（GAP, SageMath）は分解
$`\pi\cong\bigoplus\rho_\lambda\otimes1_{m_\lambda}`$ を記号的に与える。物理学の厳密対角化コード（QuSpin、QuTiP の対称性モジュール）は同じブロック構造を用いて、ハミルトニアンをセクターごとに対角化する。本ライブラリは構造を出発点とする：ひとたびブロックがわかれば、その後のすべての計算——関数、トレース、状態、チャネル、ダイナミクス——は環の中にとどまる。

<a id="x-2"></a>

## 演習問題

講義ノートの[演習問題 2](../../lecture_notes/ja/02_von_neumann_algebras.md#x-2)：$`\mathbb C^3`$ 上の巡回シフト $`S`$、
$`Se_i=e_{i+1\bmod3}`$ について、$`\{S\}'`$ と $`\{S\}''`$ はともに $`\mathrm{span}\{1,S,S^2\}`$、すなわち巡回行列全体に等しい。
次元を数値的に確かめよ。

<details>
<summary>解答</summary>

```python
S = torch.roll(torch.eye(3, dtype=dt), 1, dims=0)
assert len(commutant_basis([S])) == 3
assert len(commutant_basis(commutant_basis([S]))) == 3
assert span_dim([torch.eye(3, dtype=dt), S, S @ S]) == 3
```

</details>

**次へ：** [チュートリアル 3. トレース、状態、エントロピー](03_traces_and_states.md)
