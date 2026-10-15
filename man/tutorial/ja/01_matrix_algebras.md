[English](../01_matrix_algebras.md) · [Español](../es/01_matrix_algebras.md) · [Français](../fr/01_matrix_algebras.md) · [Deutsch](../de/01_matrix_algebras.md) · [中文](../zh/01_matrix_algebras.md) · **日本語** · [Русский](../ru/01_matrix_algebras.md)

# チュートリアル 1. 行列環と関数計算

このチュートリアルは講義ノートの[第1章](../../lecture_notes/ja/01_matrix_algebras.md)に対応する。第1章では、$`M_n(\mathbb C)`$ におけるスペクトル定理、関数計算、極分解、トレースノルムの双対性を証明する。ここでは、これらの結果の一つ一つをライブラリで検証する。

**本章で学ぶこと**：指定したスペクトルをもつ作用素の構成法、作用素への関数の適用法、作用素のバッチの格納法、そして計算を厳密な結果と照合して検証する方法。

<a id="s-1-5"></a>

## 1.5 ライブラリにおける作用素

*講義ノート：* [1.2節](../../lecture_notes/ja/01_matrix_algebras.md#s-1-2)と[1.3節](../../lecture_notes/ja/01_matrix_algebras.md#s-1-3)：[定理 1.3](../../lecture_notes/ja/01_matrix_algebras.md#r-1-3)、[系 1.4](../../lecture_notes/ja/01_matrix_algebras.md#r-1-4)、[定理 1.5](../../lecture_notes/ja/01_matrix_algebras.md#r-1-5)、[定理 1.6](../../lecture_notes/ja/01_matrix_algebras.md#r-1-6)、[定理 1.11](../../lecture_notes/ja/01_matrix_algebras.md#r-1-11)、[命題 1.12](../../lecture_notes/ja/01_matrix_algebras.md#r-1-12)、[命題 1.13](../../lecture_notes/ja/01_matrix_algebras.md#r-1-13)。

ライブラリは環 $`\bigoplus_cM_{k_c}(\mathbb C)`$（第2章）を扱う。単一のブロックは $`C=1`$ の場合であり、`TypeIAlgebra([n], [n])` が $`M_n(\mathbb C)`$ である。（二つのリストはそれぞれ名目上のブロックサイズと実際のブロックサイズであり、当面は両者は等しい。）

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra

torch.manual_seed(0)
M3 = TypeIAlgebra([3], [3], precision='double')       # the algebra M_3(C), complex128

# a positive operator with eigenvalues 1, 2, 4 and Haar-random eigenvectors
A = M3.operator_from_eigenvalues(lambda d: torch.tensor([[1.0, 2.0, 4.0]]),
                                 force_positive=True)
print(A.matrix.shape)          # (batch, C, k_max, k_max) = (1, 1, 3, 3)
print(A.eigenvalues()[0])      # one tensor per block, eigenvalues in ascending order
```

各作用素は形状 `(batch, C, k_max, k_max)` の一つのテンソルとして格納される。すなわち作用素のバッチであって、各作用素は `C` 個のブロックのリストである。ここでは `batch = C = 1` である。この構成はまさに定理 1.3 を逆向きに読んだものである：$`U`$ をハール分布に従うユニタリとして $`A=U\Lambda U^*`$（ライブラリは $`U`$ を定理 5.3 の QR 構成で生成する）。フラグ `force_positive=True` は指定された固有値に対して検査される（そのいずれかが $`-10^{-12}`$ を下回るとライブラリはエラーを送出する）。系 1.4(3) により、このとき結果は正である（ここでのように指定された固有値が $`\ge0`$ であれば厳密に、一般にはこの許容誤差の範囲で）。ライブラリはこのことを記録し、改めて検査はしない。同じ理由で、上の `A.eigenvalues()` は指定されたスペクトルをソートして返すだけで何も計算しない。これは検証ではない。一方、以下の $`\mathrm{spec}\,f(A)=f(\mathrm{spec}\,A)`$ の検証は実際に計算を行う：`apply_function` は $`A`$ の行列を対角化し、得られた行列 $`e^A`$ の固有値は一から計算される。

自己共役作用素に対する関数計算は `apply_function` である（各ブロックを `eigh` で対角化する）。関数 $`f`$ は要素ごとに作用しなければならない：同じサイズ $`k`$ のすべてのセクターの固有値を、形 `(batch, n·k)` の一つのテンソルとして受け取り（$`n`$ はそのようなセクターの個数）、同じ形のテンソルを返す。よく使う関数には専用のメソッドがある。定理 1.5 と 1.6 を検証しよう。

```python
def close(X, Y, tol=1e-10):
    """Maximal entry of |X - Y| for two operators (or tensors) is below tol."""
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

R = A.sqrt()
assert close(R @ R, A)                                  # Theorem 1.6
assert close(A.log().expm(), A)                         # exp(log A) = A, Theorem 1.5(e)
assert close(A.apply_function(lambda w: w ** 3), A @ A @ A)   # Theorem 1.5(b)
assert close(A.power(-1), A.inverse())                  # f(t) = 1/t
spec = A.apply_function(torch.exp).eigenvalues()[0]     # spec f(A) = f(spec A)
assert torch.allclose(spec, torch.exp(torch.tensor([[1.0, 2.0, 4.0]], dtype=spec.dtype)))
```

`R @ R` はこの時点ではまだ何も計算しない。積、和、関数は計算の手順（レシピ）を組み立てるだけであり、`.matrix` やトレースなどの数値が要求されたときに評価される。

**微分.** [定理 1.11](../../lecture_notes/ja/01_matrix_algebras.md#r-1-11)（ダレツキー–クレイン）により、自己共役な方向 $`E`$ への $`f(A)`$ の微分は $`U(L\circ U^*EU)U^*`$ である。ここで $`L_{ij}=f^{[1]}(\lambda_i,\lambda_j)`$ は差分商からなるレヴナー行列である（固有値が等しいときは $`f'(\lambda_i)`$）。`functional.divided_differences` は固有値と $`f`$、$`f'`$ の値から $`L`$ を計算し、autograd はそれを用いて `apply_function` を微分する。縮退したスペクトルにおいて、密な指数関数の中心差分商と公式を比較する。中心差分商の誤差は $`O(s^2)`$ に $`\varepsilon/s`$ 程度の丸めを加えたものである：

```python
from torch_vn_algebra import functional

lam = torch.tensor([1.0, 1.0, 2.0], dtype=torch.float64)          # a degenerate spectrum
Q = torch.linalg.qr(torch.randn(3, 3, dtype=torch.complex128))[0]  # a unitary
A0 = Q @ torch.diag(lam).to(torch.complex128) @ Q.mH
E = torch.randn(3, 3, dtype=torch.complex128)
E = (E + E.mH) / 2                                                 # a self-adjoint direction
L = functional.divided_differences(lam, torch.exp(lam), torch.exp(lam))
D = Q @ (L * (Q.mH @ E @ Q)) @ Q.mH                                # Theorem 1.11 for f = exp
s = 1e-6
fd = (torch.linalg.matrix_exp(A0 + s * E) - torch.linalg.matrix_exp(A0 - s * E)) / (2 * s)
print(f"{(D - fd).abs().max().item():.1e}")
assert (D - fd).abs().max() < 1e-8
```

距離 $`h`$ の近い固有値 $`a,b`$ では、商 $`\hat q`$ は桁落ちで精度を失い、中点値 $`\hat m=(f'(a)+f'(b))/2`$ の誤差は高々 $`T=h^2\max|f'''|/12`$ に丸めを加えたものである（[命題 1.12](../../lecture_notes/ja/01_matrix_algebras.md#r-1-12)）。`divided_differences` は $`|\hat q-\hat m|\le2R`$ のとき、かつそのときに限り $`\hat m`$ をとる。ここで $`R`$ は[命題 1.13](../../lecture_notes/ja/01_matrix_algebras.md#r-1-13) の丸め誤差の評価で、評価誤差には $`\varepsilon(|\hat f_a|+|a\hat f'_a|)`$ を用いる。このときその誤差は高々 $`\min(T,3R)`$ である。$`c=10^4`$ 付近の $`f(x)=e^{x-c}`$ でこれを確かめる。ここでは $`f^{[1]}(a,b)=e^{a-c}(e^h-1)/h`$ を `expm1` で正確に計算できる。最後の距離 $`h=0.06`$ は以前の版の閾値 $`\varepsilon^{1/3}\max(|a|,|b|)`$ を下回り、以前の版はそこで中点値をとっていた：

```python
eps = torch.finfo(torch.float64).eps
c = 1e4
f = lambda x: torch.exp(x - c)
for h in (1e-12, 1e-9, 1e-6, 1e-4, 1e-2, 0.06):
    w = torch.tensor([c, c + h], dtype=torch.float64)
    h = (w[1] - w[0]).item()                                       # the distance actually stored
    fw = f(w)                                                      # f' = f
    exact = math.expm1(h) / h                                      # f[a, b], as e^(a - c) = 1
    q = ((fw[1] - fw[0]) / h).item()                               # quotient
    m = (fw.sum() / 2).item()                                      # midpoint value
    y = functional.divided_differences(w, fw, fw)[0, 1].item()     # the library's choice
    R = eps * (fw.abs() + (w * fw).abs()).sum().item() / h + 4 * eps * abs(q)
    T = h ** 2 * math.exp(h) / 12 + 2 * eps                        # Proposition 1.12, plus rounding
    print(f"h = {h:.0e}: quotient {abs(q - exact):.0e}, midpoint {abs(m - exact):.0e}, library {abs(y - exact):.0e}")
    assert abs(y - exact) <= min(T, 3 * R) + 2 * eps              # Proposition 1.13(b)
```

*数値的観察*：商の誤差は $`h\le10^{-6}`$ で最大 $`5\cdot10^{-10}`$、$`h=0.06`$ で $`10^{-15}`$ である。中点値は $`h\le10^{-9}`$ では丸めを除いて正確で、$`h=0.06`$ では誤差 $`3\cdot10^{-4}`$ である。ライブラリは $`h=10^{-4}`$ を除くすべての行で二つの誤差の小さいほうを得る。$`h=10^{-4}`$ では商のほうがよい（$`10^{-12}`$）にもかかわらず中点値（$`8\cdot10^{-10}`$）をとる。$`R`$ が $`\varepsilon|a|\approx2\cdot10^{-12}`$ だけ摂動された引数を許すためで、これは `torch.exp` には悲観的だが、$`\tau x`$ を丸める $`e^{-i\tau x}`$ のような関数には必要である。保証 $`\min(T,3R)`$ はすべての行で成り立つ。

<a id="s-1-6"></a>

## 1.6 ノルムとトレースノルムの双対性

*講義ノート：* [1.4節](../../lecture_notes/ja/01_matrix_algebras.md#s-1-4)：[定理 1.8](../../lecture_notes/ja/01_matrix_algebras.md#r-1-8)、三つのノルム、[定理 1.9](../../lecture_notes/ja/01_matrix_algebras.md#r-1-9)。

一般の（正規とは限らない）行列は `alg.operator(tensor)` によって導入する。三つのノルムを PyTorch で直接計算した特異値と照合し、さらに極分解（定理 1.8）と双対性（定理 1.9）を検証する。

```python
G = torch.randn(1, 1, 3, 3, dtype=torch.complex128)    # Ginibre matrix: no symmetry at all
X = M3.operator(G)
s = torch.linalg.svdvals(G[0, 0])                       # singular values, descending

assert abs(X.operator_norm().item() - s[0].item()) < 1e-12
assert abs(X.frobenius_norm().item() - s.norm().item()) < 1e-12
assert abs(X.trace_norm().item() - s.sum().item()) < 1e-12

absX = X.abs()
W = X @ absX.inverse()                                  # |X| is invertible almost surely, see below
assert close(W.adjoint() @ W, M3.identity())            # W is unitary
assert abs((X @ W.adjoint()).trace.item() - s.sum().item()) < 1e-10    # Tr(X W*) = ||X||_1
assert abs(X.trace.item()) <= X.trace_norm().item()     # |Tr X| <= ||X||_1
```

複素数の dtype を指定した `torch.randn` は、実部と虚部が独立に $`\sim N(0,\tfrac12)`$ に従う独立な成分を生成するので、`G` はジニブル行列であり、密度をもつ（命題 A.16(b)）。
$`|X|`$ はほとんど確実に可逆である：$`(\det|X|)^2=\det(|X|^2)=\det(X^*X)=|\det X|^2`$ であり、$`|\det X|^2`$ は成分の実部と虚部の多項式で、$`X=1`$ において $`1`$ に等しいから、それが零となる確率は零である（系 A.7(b)。これは定理 5.3 の第一段階でもある）。

<a id="s-1-7"></a>

## 1.7 バッチ：多数の作用素を一度に

*講義ノート：* [1.7節](../../lecture_notes/ja/01_matrix_algebras.md#s-1-7)、[命題 1.10](../../lecture_notes/ja/01_matrix_algebras.md#r-1-10)。

バッチ軸こそがこのライブラリの存在理由である。ランダム行列理論、モンテカルロ法、マルチスタート最適化には数千個の作用素が必要であり、バッチ化された一回の呼び出しがそれらにわたる Python のループに取って代わる。例として、**ガウス型ユニタリアンサンブル**（GUE）の厳密なモーメントを検証する。

[講義ノート、命題 1.10](../../lecture_notes/ja/01_matrix_algebras.md#r-1-10)により、$`G`$ が独立な成分をもち、各成分の実部と虚部が独立に $`\sim N(0,\tfrac12)`$ に従い、$`X=(G+G^*)/\sqrt{2n}`$ であるならば、$`\mathbb E\,\tfrac1n\mathrm{Tr}\,X^4=2+\tfrac1{n^2}`$ である。

```python
n, B = 4, 200_000
G = torch.complex(torch.randn(B, 1, n, n), torch.randn(B, 1, n, n)).to(torch.complex128) / math.sqrt(2)
M4 = TypeIAlgebra([n], [n], precision='double')
X = M4.operator((G + G.conj().transpose(-2, -1)) / math.sqrt(2 * n), is_self_adjoint=True)

X2 = X @ X
m4 = (X2 @ X2).trace.real / n                           # shape (B,): one value per sample
mean, err = m4.mean().item(), m4.std().item() / math.sqrt(B)
print(f"E tr X^4 / n = {mean:.4f} +- {err:.4f}   (exact {2 + 1 / n**2:.4f})")
assert abs(mean - (2 + 1 / n**2)) < 5 * err
```

モンテカルロ法による検証は統計的検定である：推定値は厳密値と標準誤差の数倍以内で一致しなければならない。変数 $`\tfrac1n\mathrm{Tr}\,X^4`$ はガウス変数の多項式であるから、有限の分散 $`\sigma^2`$ をもつ（補題 A.23：ガウス分布のすべてのモーメントは有限である。積は $`|y_1\cdots y_r|\le\sum_i|y_i|^r`$ で扱う）。チェビシェフの不等式により、$`B`$ 個の独立な標本の平均が厳密値から $`5\sigma/\sqrt B`$ より大きくずれる確率は高々 $`1/25`$ である。コードでは $`\sigma`$ の代わりに標本標準偏差を用い、固定したシードによって検証は再現可能になる。これが本チュートリアルにおける証明済みの結果との照合すべてに共通する型である：厳密な主張、数値的推定、そして明示的な許容誤差。

<a id="s-1-8"></a>

## 1.8 比較の対象

単一の行列については、`numpy.linalg` と `scipy.linalg`（`sqrtm`, `expm`, `logm`,
`polar`）が参照実装であり、上の検証はそれらを用いても繰り返すことができる。本ライブラリは、後で重要になる三つのものを付け加える：バッチ化（1.7節）、ブロック構造を担う作用素（第2章）、そしてその演算を通じた自動微分（PyTorch の autograd）である。

<a id="x-1"></a>

## 演習問題

講義ノートの[演習問題 1 と 3](../../lecture_notes/ja/01_matrix_algebras.md#x-1)を数値的に確かめる：
- ランダムな自己共役行列 $`A`$ について、$`f=\sin`$, $`g=\exp`$ に対し $`(fg)(A)=f(A)g(A)`$ を確かめよ（演習問題 1）。
- $`A,B>0`$ の **幾何平均** $`A\#B=A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}`$ は、$`XA^{-1}X=B`$ のただ一つの正の解
  $`X`$ である（演習問題 3）。この方程式と $`A\#B=B\#A`$ を数値的に確かめよ。

<details>
<summary>解答</summary>

```python
Hs = M3.operator_from_eigenvalues(lambda d: torch.randn(1, d), force_self_adjoint=True)
assert close(Hs.apply_function(lambda w: torch.sin(w) * torch.exp(w)),
             Hs.apply_function(torch.sin) @ Hs.apply_function(torch.exp))
```

```python
def geometric_mean(A, B):
    a, ai = A.sqrt(), A.power(-0.5)
    return a @ (ai @ B @ ai).sqrt() @ a

P = M3.operator_from_eigenvalues(lambda d: 0.5 + torch.rand(1, d), force_positive=True)
Q = M3.operator_from_eigenvalues(lambda d: 0.5 + torch.rand(1, d), force_positive=True)
GM = geometric_mean(P, Q)
assert close(GM, geometric_mean(Q, P), 1e-9)
assert close(GM @ P.inverse() @ GM, Q, 1e-9)
```

</details>

**次へ：** [チュートリアル 2](02_von_neumann_algebras.md)：ライブラリにおける可換子環、中心、環 $`\bigoplus_cM_{k_c}(\mathbb C)`$。
