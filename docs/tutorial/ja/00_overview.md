[English](../00_overview.md) · [Español](../es/00_overview.md) · [Français](../fr/00_overview.md) · [Deutsch](../de/00_overview.md) · [中文](../zh/00_overview.md) · **日本語** · [Русский](../ru/00_overview.md)

# 0. 本講義の主題

フォン・ノイマン環は、ほとんどもっぱら紙と鉛筆で研究されている。有限次元では、そのどれもが全行列環の直和である（第2章）から、原理的には何でも計算できる——しかし実際には、環はたいてい一つの大きな行列へと平坦化され、その構造（セクター、中心、多数のトレース）は計算から消えてしまう。本講義とライブラリ `torch_vn_algebra` は逆の方針をとる：**環そのものがデータ型** であり、その中での計算を効率的に実行する。

## 0.1 基本的な観察

有限次元フォン・ノイマン環は、ユニタリな基底変換を除いて
$$
\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c}
$$
である（定理 2.6）。抽象的な環としては $\bigoplus_cM_{k_c}(\mathbb C)$ であり、その元は
$C$ 個の正方ブロックのリストである。これはまさに *セクター軸* をもつテンソルの形である：
$\mathcal M$ の元は形状 $(C,k,k)$ の配列であり（ブロックは $k=\max_ck_c$ までパディングする）、$B$ 個の元——モンテカルロ標本、最適化の初期点、時間ステップ、乱れの実現——は形状 $(B,C,k,k)$ の配列をなす。現代の配列ライブラリ（CPU および GPU 上の PyTorch）は、ニューラルネットワークが必要とするため、まさにこのような「小さな行列のバッチ」のために作られている。

この同一視のもとで、環の演算は、厳密に近似なしに、バッチ化された配列演算となる。

| 環 $\mathcal M$ において | 形状 $(B,C,k,k)$ の配列上で |
|---|---|
| 積 $AB$、和、随伴 | $(B,C)$ にわたるバッチ化された行列積、和、共役転置 |
| 関数計算 $f(A)$（第1章） | ブロックのバッチ化された固有値分解 |
| トレース $\sum_cw_c\operatorname{Tr}A_c$（第3章） | ブロックごとの対角和と、セクター軸にわたる重み付き和 |
| 中心、中心射影（第2章） | セクターごとに一つの数 |
| 中心の上への条件付き期待値 | ブロックのトレースを $k_c$ で割り、ブロードキャストで戻す |
| ユニタリ群 $\prod_cU(k_c)$ 上のハール測度 | ガウス型ブロックのバッチ化された QR 分解（第5章） |
| 状態、チャネル、ダイナミクス（第3, 6, 8章） | 同じ演算の合成 |

二つの帰結により、これは単なる帳簿付け以上のものとなる。

- **計算量.** $\mathcal M$ における積またはスペクトル分解の計算量は $\sum_ck_c^3$ 回の演算であり、
  $M_D(\mathbb C)$（$D=\sum_ck_c$）への埋め込みでの $D^3$ と比べられる。$C$ 個の等しいブロックの場合、その比は
  $C^2$ 倍である。積、和、トレースは $B\cdot C$ 個のブロックすべてを一回の呼び出しで処理する。分解（固有値、特異値）は相異なるブロックサイズごとに一回の呼び出しとして実行される——すべてのセクターが同じサイズなら一回の呼び出しである。これにより GPU を休ませない。
- **構成による正しさ.** 計算されるものはすべて $\mathcal M$ の中にとどまる：ランダムな
  ユニタリは $\prod_cU(k_c)$ 上でハール分布に従い（これは $U(D)$ の中の零集合であり、
  $U(D)$ からのサンプリングでは決して到達しない）、勾配ステップは環から出ず、トレース、密度、
  条件付き期待値は偶然の埋め込みではなく環そのものを参照する。

## 0.2 最初の計算

同じ計算を $\mathcal M=\bigoplus_{c=1}^{64}M_{16}(\mathbb C)$ と、その
$M_{1024}(\mathbb C)$ への埋め込みとで行って比較する：ランダムな正の元の平方根と、そのトレースである。

```python
import time
import torch
from torch_vn_algebra import TypeIAlgebra

torch.manual_seed(0)
C, k = 64, 16
alg = TypeIAlgebra([k] * C, [k] * C, precision='double')       # M = (+)_{c=1}^{64} M_16(C)

A = alg.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(1, d), force_positive=True)

t0 = time.perf_counter()
root_blocks = A.sqrt().matrix                                    # (1, 64, 16, 16): 64 small eigenproblems
t_blocks = time.perf_counter() - t0

big = torch.block_diag(*A.matrix[0])                             # the same element in M_1024(C)
t0 = time.perf_counter()
w, V = torch.linalg.eigh(big)                                    # one large eigenproblem
root_big = (V * w.sqrt()) @ V.conj().T
t_big = time.perf_counter() - t0

assert torch.allclose(torch.block_diag(*root_blocks[0]), root_big, atol=1e-10)   # same result
assert abs(A.sqrt().Tr_blunt().item() - root_big.trace().item()) < 1e-9
print(f"in the algebra: {t_blocks * 1e3:.1f} ms, as one 1024 x 1024 matrix: {t_big * 1e3:.1f} ms")
```

結果は丸め誤差の範囲で一致する。環を用いる版のほうが速い——この大きさでは典型的な CPU で約十倍である（本講義をテストした4コアのマシンで約 16 ms 対 200 ms。64 個のブロックすべてが一回の呼び出しで固有値ソルバーに渡される）——が、演算回数の比 $C^2=4096$ が示唆するよりはるかに小さい。一回の大きな固有値ソルバーの呼び出しは、多数の小さな呼び出しよりもハードウェアを有効に使うからである。この差はセクター数とともに広がり、メモリの節約は $C$ 倍である。$B$ 個の元のバッチでは、環を用いる版は依然として一回の呼び出しであるが、大きな行列を用いる版は $B$ 個の大きな固有値問題となる。

## 0.3 既存のものと新しいもの

有限次元フォン・ノイマン環の *内部で* 数値計算を行うソフトウェアを調査した。見つかったものは四つのグループに分かれ、いずれも有用であるが目的は異なる。

- **汎用の行列ソフトウェア**（NumPy/SciPy、QuTiP その他の量子ツールボックス）は一度に一つの
  行列を扱う。ブロック構造があっても、それは手作業で維持しなければならない。
- **対称テンソルネットワークのライブラリ** は、テンソルをアーベル対称性の電荷セクター上のブロック疎配列として格納し、
  中には GPU と自動微分のための PyTorch バックエンドをもつものもある
  （例えば YASTN [1]；ほかに TeNPy, ITensor, TensorKit）。技術的にはこれらが最も近い親戚である。その
  対象は縮約すべきテンソルネットワークであって作用素環ではない：それらはトレースの族、状態と密度、中心と条件付き期待値、
  $\prod_cU(k_c)$ 上のハール測度、セクター間のチャネルを提供しない。
- 行列 $*$-環の **数値的ブロック対角化** [2, 3, 4] は、与えられた環の分解
  $\mathcal M\cong\bigoplus_cM_{k_c}\otimes1_{m_c}$ を求める——典型的には対称性をもつ半正定値計画問題の大きさを縮小するためである。これは相補的である：それらのアルゴリズムはブロックを見つけ、
  本ライブラリはその中で計算する。群環の厳密な記号的分解は計算機代数に存在する
  （例えば GAP パッケージ Wedderga）。
- **非可換多項式最適化**（NCSOStools, ncpol2sdpa）は、非可換変数の多項式の固有値とトレースの上下界を、すべての行列サイズにわたって一度に、半正定値緩和によって与える [5]。これは異なる問い（すべての環で成り立つ上下界）に異なる手段（凸緩和）で答えるものである。

我々の知る限り、有限次元フォン・ノイマン環を、バッチ化され微分可能な数値計算の基本データ型とする汎用ライブラリは存在しない。したがって新規性は新しい数学ではなく **効率的な計算** にある：§0.1 の同一視を、環の構造が保たれ、すべての演算がバッチ化されたカーネル（分解では相異なるブロックサイズごとに一つ）として実行されるように実装したことである。関連するソフトウェアをご存じであれば、ぜひお知らせいただきたい。

## 0.4 本講義の読み方

| 章 | 内容 |
|---|---|
| 1–3 | 数学：行列環、構造定理、トレース、状態、エントロピー——完全な証明付き |
| 4–8 | これらの環はどこから来て、その中で何を計算するのか：対称性、ランダム行列、チャネル、テンソル積、ダイナミクス |
| 9–12 | 応用：スピン鎖、量子計算、最適化、同変機械学習 |

各章は結果を証明とともに述べ、ライブラリのコードを示し、コードを厳密な結果と照合して検証し、比較すべき標準的な道具を挙げる。

## 参考文献

1. M. M. Rams, G. Wójtowicz, A. Sinha, J. Hasik, *YASTN: Yet another symmetric tensor networks; A
   Python library for Abelian symmetric tensor network calculations*, SciPost Phys. Codebases 52
   (2025).
2. K. Murota, Y. Kanno, M. Kojima, S. Kojima, *A numerical algorithm for block-diagonal
   decomposition of matrix $*$-algebras*, Technical Report METR 2007-52, University of Tokyo (2007).
3. T. Maehara, K. Murota, *A numerical algorithm for block-diagonal decomposition of matrix
   $*$-algebras, Part II: general algorithm*, Optimization Online (2009).
4. E. de Klerk, C. Dobre, D. V. Pasechnik, *Numerical block diagonalization of matrix
   $*$-algebras with application to semidefinite programming*, Math. Program. 129 (2011) 91–111.
5. I. Klep, J. Povh, *Semidefinite programming and sums of hermitian squares of noncommutative
   polynomials*, J. Pure Appl. Algebra 214 (2010) 740–749.

**次へ：**[第1章](01_matrix_algebras.md)。
