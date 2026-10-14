[English](../00_overview.md) · [Español](../es/00_overview.md) · [Français](../fr/00_overview.md) · [Deutsch](../de/00_overview.md) · [中文](../zh/00_overview.md) · **日本語** · [Русский](../ru/00_overview.md)

# チュートリアル 0. 本講義の主題

このチュートリアルは講義ノートの[第0章](../../lecture_notes/ja/00_overview.md)に対応する。第0章では、フォン・ノイマン環 $`\bigoplus_cM_{k_c}(\mathbb C)`$ の元を形状 $`(C,k,k)`$ の配列と同一視し、その演算がバッチ化された配列演算であることを証明する。ここでは、これがライブラリ `torch_vn_algebra` において何をもたらすかを見る：その計算量、最初の計算、そして本ライブラリと既存のソフトウェアとの関係である。

フォン・ノイマン環は、主に紙と鉛筆で研究されている。有限次元では、そのどれもが全行列環の直和と $`*`$-同型である（定理 2.6）から、原理的には何でも計算できる——しかしそのような環を数値的に扱うとき、環はしばしば一つの大きな行列へと平坦化され、その構造（セクター、中心、多数のトレース）は計算から消えてしまう。本講義とライブラリ `torch_vn_algebra` は逆の方針をとる：**環そのものがデータ型** であり、その中での計算を効率的に実行する。

<a id="s-0-1"></a>

## 0.1 基本的な観察

*講義ノート：* [0.1節](../../lecture_notes/ja/00_overview.md#s-0-1)。そこでは[定理 2.6](../../lecture_notes/ja/02_von_neumann_algebras.md#r-2-6)、[定理 3.1](../../lecture_notes/ja/03_traces_and_states.md#r-3-1)、[定理 5.3](../../lecture_notes/ja/05_random_matrices.md#r-5-3)、[系 5.4](../../lecture_notes/ja/05_random_matrices.md#r-5-4)から、次の表の最初の六行が恒等式であることを証明している：

| 環 $`\mathcal M`$ において | 形状 $`(B,C,k,k)`$ の配列上で |
|---|---|
| 積 $`AB`$、和、随伴 | $`(B,C)`$ にわたるバッチ化された行列積、和、共役転置 |
| 関数計算 $`f(A)`$（第1章） | ブロックのバッチ化された固有値分解 |
| トレース $`\sum_cw_c\mathrm{Tr}\,A_c`$（第3章） | ブロックごとの対角和と、セクター軸にわたる重み付き和 |
| 中心、中心射影（第2章） | セクターごとに一つの数 |
| 中心の上への条件付き期待値 | ブロックのトレースを $`k_c`$ で割り、ブロードキャストで戻す |
| ユニタリ群 $`\prod_cU(k_c)`$ 上のハール測度 | ガウス型ブロックのバッチ化された QR 分解と、定理 5.3 の位相補正 |
| 状態、チャネル、ダイナミクス（第3, 6, 8章） | 上の演算から組み立てられる |

二つの帰結により、これは単なる帳簿付け以上のものとなる。

- **計算量.** $`\mathcal M`$ におけるスペクトル分解はブロックごとに計算される（表の 2 行目）。標準的な密行列用固有値ソルバー（PyTorch が CPU 上で呼び出す LAPACK）が用いる浮動小数点演算の回数は、行列サイズの 3 乗のように増える。これはそれらのアルゴリズムの記述であって、本講義の結果ではない。$`k\times k`$ の固有値問題の計算量が定数倍の $`k^3`$ であるというモデルでは、$`\mathcal M`$ におけるスペクトル分解の計算量は $`\sum_ck_c^3`$ であり、
  $`M_D(\mathbb C)`$（$`D=\sum_ck_c`$）への埋め込みでの $`D^3`$ と比べられる。サイズ $`k`$ の等しいブロックが $`C`$ 個の場合、その比は
  $`(Ck)^3/(Ck^3)=C^2`$ である。積はゼロで埋めたテンソル上で実行される。二つの $`k\times k`$ 行列の教科書的な積は $`k^3`$ 回の乗算を用いる（$`k^2`$ 個の成分のそれぞれが $`k`$ 個の積の和である）から、$`\mathcal M`$ における積は $`C\,k_{\max}^3`$ 回を用い、これはすべてのブロックが同じサイズのとき $`\sum_ck_c^3`$ に等しい。元は $`D^2`$ 個ではなく $`C\,k_{\max}^2`$ 個の数として格納される。等しいブロックの場合、その比は $`(Ck)^2/(Ck^2)=C`$ である。コードが実際に行うこと：積、和、トレースは $`(B,C,k_{\max},k_{\max})`$ の配列全体に対する一回のテンソル演算である。エルミート行列の固有値分解、特異値、QR 分解は、相異なるサイズ $`k`$ ごとに、そのサイズの $`B\cdot n_k`$ 個のブロック（$`n_k`$ はそのサイズのセクターの数）すべてのスタックに一回のバッチ呼び出しで適用される——すべてのセクターが同じサイズなら呼び出しは一回だけである。$`t\ge2`$ 個の PyTorch スレッドをもつ CPU では、スタックが少なくとも $`2t`$ 個のブロックをもち、かつ（ブロック数）$`\,\cdot k^3\ge2^{18}`$ であれば（モジュール `cost`、関数 `_run`）、ライブラリはそのような呼び出しを 高々 $`t`$ 個の塊に分割して並列スレッドで実行する。結果は同じ行列に対する同じ LAPACK 呼び出しである。目的は GPU、あるいは CPU のすべてのコアを休ませないことである。
- **構成による正しさ.** ゼロでパディングされた形状 $`(C,k_{\max},k_{\max})`$ の配列はどれも $`\mathcal M`$ の元であるから、ブロックから計算されるものはすべて $`\mathcal M`$ の中にとどまる。ライブラリのランダムなユニタリは、定理 5.3 の構成によって各ブロックで独立に生成される（`random_unitary_operator` はブロックごとのサンプラー `random_unitary` をブロックサイズごとに一回呼び出し、5.1 節で `random_unitary` をこの構成と照合している）。これは系 5.4 により $`\prod_cU(k_c)`$ のハール測度である（理想的なガウス乱数に対して）。$`C\ge2`$ のとき、$`U(D)`$ のハール・ユニタリがこの部分群に属する確率はゼロである（5.1 節、系 5.4 の後）。勾配ステップはブロックを変えるので、環から出ることはありえない。そしてトレースは、明示的に選んだ重みをもつ定理 3.1 の汎関数 $`\sum_cw_c\mathrm{Tr}\,A_c`$ であり（密度と条件付き期待値もそれに関してとる）、$`\mathcal M`$ がたまたま作用している何らかの $`\mathbb C^n`$ のトレース——それは $`w_c=m_c`$ と固定してしまう（2.3 節）——ではない。

<a id="s-0-2"></a>

## 0.2 最初の計算

*講義ノート：* [0.1節](../../lecture_notes/ja/00_overview.md#s-0-1)、表の 2 行目：$`f(A)=\bigoplus_cf(A_c)`$、ブロックごとに一回の固有値分解。

同じ計算を $`\mathcal M=\bigoplus_{c=1}^{64}M_{16}(\mathbb C)`$ と、その
$`M_{1024}(\mathbb C)`$ への埋め込みとで行って比較する：ランダムな正の元の平方根と、そのトレースである。
`TypeIAlgebra(n_factors, k_factors)` はブロックサイズの二つのリストを受け取る：名目上のサイズ $`n_c`$ と、実際に使われるブロックのサイズ $`k_c\le n_c`$ である（1.5 節）。本講義が `TypeIAlgebra` で構成する環はどれも、二つのリストが等しい。

```python
import time
import torch
from torch_vn_algebra import TypeIAlgebra

torch.manual_seed(0)
C, k = 64, 16
alg = TypeIAlgebra([k] * C, [k] * C, precision='double')       # M = (+)_{c=1}^{64} M_16(C)

A = alg.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(1, d), force_positive=True)
blocks = A.matrix                                                # build A before timing

t0 = time.perf_counter()
root_blocks = A.sqrt().matrix                                    # (1, 64, 16, 16): 64 small eigenproblems
t_blocks = time.perf_counter() - t0

big = torch.block_diag(*blocks[0])                               # the same element in M_1024(C)
t0 = time.perf_counter()
w, V = torch.linalg.eigh(big)                                    # one large eigenproblem
root_big = (V * w.sqrt()) @ V.conj().T
t_big = time.perf_counter() - t0

assert torch.allclose(torch.block_diag(*root_blocks[0]), root_big, atol=1e-10)   # row 2: f(A) = (+)_c f(A_c)
assert abs(A.sqrt().Tr_blunt().item() - root_big.trace().item()) < 1e-9
assert big.numel() == C * blocks.numel()                          # storage: a factor C = 64
print(f"in the algebra: {t_blocks * 1e3:.1f} ms, as one 1024 x 1024 matrix: {t_big * 1e3:.1f} ms")
```

結果は丸め誤差の範囲で一致し（表の 2 行目）、記憶容量の比は上で数えたとおり $`C`$ である。64 個のブロックは形状 $`(64,16,16)`$ のスタックに対する一回のバッチ呼び出しとして固有値ソルバーに渡される。PyTorch のスレッド数が $`2\le t\le32`$ のとき、ライブラリはそれを 高々 $`t`$ 個の塊に分割し（$`64\cdot16^3=2^{18}`$ はモジュール `cost` のしきい値にちょうど達する）、スレッドが一つのときは分割しない。*数値的観察*（時間はハードウェアに依存し、コードがそれを出力する）：二台の 4 コア CPU で、環を用いる版は 3–10 ms、大きな行列を用いる版は 400–760 ms かかり、その比は数十倍から数百倍であった（二台目のマシンでは 4 スレッドでも 1 スレッドでも同様）。測定された比は、3 乗の計算量モデルの比 $`C^2=4096`$ をはるかに下回る：演算回数は実行時間のモデルであって実行時間そのものではなく、これらのマシンでは、小さな固有値問題は大きな固有値問題よりもはるかに低い毎秒演算回数で実行された。計算量モデルでは、比 $`C^2`$ はセクター数とともに増大する。$`B`$ 個の元のバッチでは、環を用いる版は依然として一回のバッチ呼び出しであるが、大きな行列を用いる版はサイズ $`1024`$ の固有値問題を $`B`$ 個解かなければならない。

<a id="s-0-3"></a>

## 0.3 既存のものと新しいもの

有限次元フォン・ノイマン環の *内部で* 数値計算を行うソフトウェアを調査した。見つかったものは四つのグループに分かれ、いずれも有用であるが目的は異なる。

- **汎用の行列ソフトウェア**（NumPy/SciPy、QuTiP その他の量子ツールボックス）は一度に一つの
  行列を扱う。ブロック構造があっても、それは手作業で維持しなければならない。
- **対称テンソルネットワークのライブラリ** は、テンソルを対称性の電荷セクター上のブロック疎配列として格納する（YASTN [1]、TeNPy、ITensor ではアーベル対称性、TensorKit は非アーベル対称性も扱う）。YASTN は GPU と自動微分のための PyTorch バックエンドをもつ。技術的にはこれらが最も近い親戚である。その
  対象は縮約すべきテンソルネットワークであって作用素環ではない：それらのドキュメントで我々が確認した限り、トレースの族、状態と密度、中心と条件付き期待値、
  $`\prod_cU(k_c)`$ 上のハール測度、セクター間のチャネルは提供されていない。
- 行列 $`*`$-環の **数値的ブロック対角化** [2, 3, 4] は、与えられた環の分解
  $`\mathcal M\cong\bigoplus_cM_{k_c}\otimes1_{m_c}`$ を求める——典型的には対称性をもつ半正定値計画問題の大きさを縮小するためである。これは相補的である：それらのアルゴリズムはブロックを見つけ、
  本ライブラリはその中で計算する。本ライブラリの入力段階（モジュール `decompose`、4.3 節）は、可換子環のランダムな元の固有空間がブロックを見つけるという [2, 3] のアイデアを用いる。4.3 節では、厳密な算術ではこれが確率 1 で成功することを証明し（補題 4.6、命題 4.7、補題 4.8）、命題 4.9 は `Decomposition.verify()` が数値的に検査する判定基準である。群環の厳密な記号的分解は計算機代数に存在する
  （例えば GAP パッケージ Wedderga）。
- **非可換多項式最適化**（NCSOStools, ncpol2sdpa）は、非可換変数の多項式の固有値とトレースの上下界を、すべての行列サイズにわたって一度に、半正定値緩和によって与える [5]。これは異なる問い（すべての行列サイズで成り立つ上下界）に異なる手段（凸緩和）で答えるものである。

*講義ノート：* ランダム化された分解アルゴリズムとその判定基準 — [4.3節](../../lecture_notes/ja/04_symmetries.md#s-4-3)、[命題 4.7](../../lecture_notes/ja/04_symmetries.md#r-4-7)、[命題 4.9](../../lecture_notes/ja/04_symmetries.md#r-4-9)。

この概観は我々が見つけたソフトウェアの記述であり、定理ではない。我々の知る限り、有限次元フォン・ノイマン環を、バッチ化され微分可能な数値計算の基本データ型とする汎用ライブラリは存在しない。したがって新規性は新しい数学ではなく **効率的な計算** にある：§0.1 の同一視を、環の構造が保たれ、演算がバッチ化されたカーネル（分解では相異なるブロックサイズごとに一回の呼び出しで、マルチコア CPU では §0.1 のサイズのしきい値を超えるとスレッドごとの塊に分割される）として実行されるように実装したことである。関連するソフトウェアをご存じであれば、ぜひお知らせいただきたい。

<a id="s-0-4"></a>

## 0.4 本講義の読み方

| 章 | 内容 |
|---|---|
| 1–3 | 数学：行列環、構造定理、トレース、状態、エントロピー——完全な証明付き |
| 4–8 | これらの環はどこから来て、その中で何を計算するのか：対称性、ランダム行列、チャネル、テンソル積、ダイナミクス |
| 9–12 | 応用：スピン鎖、量子計算、最適化、同変機械学習 |
| 13 | 基礎：ジョルダン代数と、量子力学が *-環を用いる理由 |
| 付録 | 各章が番号で引用する補助的な結果の完全な証明（例えば定理 A.24） |

各チュートリアルは、同じ番号の[講義ノート](../../lecture_notes/ja/00_overview.md)の章に対応する：ライブラリのコードを示し、そこで証明された結果とコードを数値的に照合して検証し、比較すべき標準的な道具を挙げる。各節は、その節と用いる結果へリンクする *講義ノート：* の行で始まる。チュートリアルのコードブロックは、一つの Python セッションで順に実行することを想定している。講義ノートが証明しない数値（例えば実行時間）は、そのチュートリアル自身のコードの *数値的観察* としてのみ現れる。

## 参考文献

1. M. M. Rams, G. Wójtowicz, A. Sinha, J. Hasik, *YASTN: Yet another symmetric tensor networks; A
   Python library for Abelian symmetric tensor network calculations*, SciPost Phys. Codebases 52
   (2025).
2. K. Murota, Y. Kanno, M. Kojima, S. Kojima, *A numerical algorithm for block-diagonal
   decomposition of matrix $`*`$-algebras with application to semidefinite programming*, Japan J.
   Indust. Appl. Math. 27 (2010) 125–160（プレプリント：Technical Report METR 2007-52, University of
   Tokyo, 2007）.
3. T. Maehara, K. Murota, *A numerical algorithm for block-diagonal decomposition of matrix
   $`*`$-algebras with general irreducible components*, Japan J. Indust. Appl. Math. 27 (2010)
   263–293.
4. E. de Klerk, C. Dobre, D. V. Pasechnik, *Numerical block diagonalization of matrix
   $`*`$-algebras with application to semidefinite programming*, Math. Program. 129 (2011) 91–111.
5. I. Klep, J. Povh, *Semidefinite programming and sums of hermitian squares of noncommutative
   polynomials*, J. Pure Appl. Algebra 214 (2010) 740–749.

**次へ：** [チュートリアル 1](01_matrix_algebras.md)：ライブラリにおける作用素、作用素の関数、バッチ。
