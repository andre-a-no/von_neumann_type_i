[English](../11_optimisation.md) · [Español](../es/11_optimisation.md) · [Français](../fr/11_optimisation.md) · [Deutsch](../de/11_optimisation.md) · [中文](../zh/11_optimisation.md) · **日本語** · [Русский](../ru/11_optimisation.md)

# チュートリアル 11. 最適化とトレース不等式

このチュートリアルは講義ノートの[第11章](../../lecture_notes/ja/11_optimisation.md)に対応する。同章は三つのトレース汎関数の極値を証明し、ユニタリの
サンプリングではそれらが見つからない理由を説明する。ここでは、ライブラリの制約集合のパラメータ化
（`optimize.PositiveParam`、`UnitaryParam`、`SelfAdjointParam`）とマルチスタート最適化器（`optimize.extremize`）を、
それらの定理と照合する。$\Delta$ をマイケルソン・コントラスト、$\ell(\delta)=(1-\delta)/(1+\delta)$ とすると、三つの汎関数は次のとおりである：

| | 汎関数 | 変数 |
|---|---|---|
| 実験 1 | $z=\lvert\mathrm{Tr}(XUY)\rvert-\mathrm{Tr}(XY)$ | $X,Y\ge0$、$U\in\mathcal M$ はユニタリ |
| 実験 2 | $z=\mathrm{Tr}\lvert XY\rvert-\mathrm{Tr}(\lvert X\rvert\,\lvert Y\rvert)$ | $X=UX_0$、$Y=VY_0$、$X_0,Y_0\ge0$、$U,V$ はユニタリ |
| 実験 3 | $z=\mathrm{Tr}(Y\lvert X\rvert Y)-\mathrm{Tr}\lvert YXY\rvert$ | $X=X^*$、$Y\ge0$ |

ここで $|X|$、$|Y|$ は $\lambda_{\max}=1$ に正規化し、それらのコントラストをパラメータとする。**包絡線** とは
$\overline z(\delta)=\sup\{z:\Delta(|X|)=\delta\}$ と $\underline z(\delta)=\inf\{z:\Delta(|X|)=\delta\}$ のことであり、その他の
変数はすべて自由に動かす。

<a id="s-11-1"></a>

## 11.1 制約集合とそのパラメータ化

*講義ノート：* [11.1節](../../lecture_notes/ja/11_optimisation.md#s-11-1)；[補題 11.1](../../lecture_notes/ja/11_optimisation.md#r-11-1)、[命題 11.2](../../lecture_notes/ja/11_optimisation.md#r-11-2)。

コントラスト $\delta$ を指定した `optimize.PositiveParam` は、講義ノートの[11.1節](../../lecture_notes/ja/11_optimisation.md#s-11-1)のパラメータ化を実装する：
各スタートごとに二つの異なるスロット $s_\pm$（`spectrum.pin_max`、`spectrum.pin_min`。どの順序対も確率 $1/(D(D-1))$ をもつ）と
ハール分布に従うランダムなユニタリ $U_0$（`unitary.base`。`random_unitary_operator` で抽出する。[補題 11.10](../../lecture_notes/ja/11_optimisation.md#r-11-10)）を
抽出し、パラメータ $\theta$（`spectrum.theta`）と $G$（`unitary.G`）を次のように写す。
$$
X(\theta,G)=U\mathrm{diag}(\lambda)U^*,\qquad U=U_0\,e^{G-G^*},\qquad
\lambda_{s_+}=1,\quad\lambda_{s_-}=\ell(\delta),\quad\lambda_s=\ell+(1-\ell)\tfrac{1+\sin\theta_s}2\ \text{otherwise}.
$$
`UnitaryParam` は写像 $G\mapsto U_0e^{G-G^*}$ 単独であり、`SelfAdjointParam` は符号を加えて $X=U\mathrm{diag}(\sigma\lambda)U^*$ とする。
[命題 11.2](../../lecture_notes/ja/11_optimisation.md#r-11-2) により、像は、抽出されたスロットに $1$ と $\ell(\delta)$ をもつコントラスト殻 $\mathcal P_\delta$ の作用素全体である。
$[a,b]$ 内の自由なコントラストはパラメータが一つ増えるだけである。ライブラリの既定値は $[a,b]=[0,0.999]$ であり、以下のコードもこれを用いる。

コードは制約を検証し、さらに命題 11.2 を構成的に検証する：証明のパラメータが目標の作用素を再現する。

```python
import math
import scipy.linalg
import torch
from torch_vn_algebra import TypeIAlgebra, optimize as opt

torch.manual_seed(0)
alg = TypeIAlgebra([3, 3], [3, 3], precision='double')          # M_3 (+) M_3, D = 6
ell = lambda delta: (1 - delta) / (1 + delta)                     # lambda_min for contrast delta
d = torch.tensor([0.0, 0.25, 0.5, 0.9], dtype=torch.float64)
x = opt.PositiveParam(alg, 4, delta=d).operator()                 # four starts, one contrast each
assert torch.allclose(x.michelson_contrast, d, atol=1e-12)        # the constraint holds exactly
assert torch.allclose(x.lambda_max, torch.ones(4, dtype=torch.float64), atol=1e-12)

P = opt.PositiveParam(alg, 1, delta=torch.tensor([0.5], dtype=torch.float64))
lo = ell(0.5)
mu = lo + (1 - lo) * torch.rand(1, 2, 3, dtype=torch.float64)    # target spectrum in [lo, 1] ...
(cp, ip), (cm, im) = P.spectrum.pin_max[0].tolist(), P.spectrum.pin_min[0].tolist()
mu[0, cp, ip], mu[0, cm, im] = 1.0, lo                            # ... with 1 and lo in the drawn slots
W = alg.random_unitary_operator(1).matrix                         # target eigenvectors
target = W @ torch.diag_embed(mu.to(W.dtype)) @ W.mH
with torch.no_grad():                                             # the parameters of the proof
    P.spectrum.theta.copy_(torch.asin((2 * (mu - lo) / (1 - lo) - 1).clamp(-1, 1)))
    V = P.unitary.base.mH @ W                                     # U0^* W = exp(S), G = S / 2
    S = torch.stack([torch.from_numpy(scipy.linalg.logm(V[0, c].numpy())) for c in range(2)])
    P.unitary.G.copy_(S[None] / 2)
err = (P.operator().matrix - target).abs().max().item()
print(f"slots of 1 and lo: {(cp, ip)}, {(cm, im)};  reconstruction error {err:.1e}")
assert err < 1e-10
```

<a id="s-11-2"></a>

## 11.2 フォン・ノイマンのトレース不等式

*講義ノート：* [11.2節](../../lecture_notes/ja/11_optimisation.md#s-11-2)；[定理 11.4](../../lecture_notes/ja/11_optimisation.md#r-11-4)、[系 11.5](../../lecture_notes/ja/11_optimisation.md#r-11-5)、[命題 11.6](../../lecture_notes/ja/11_optimisation.md#r-11-6)、[系 11.7](../../lecture_notes/ja/11_optimisation.md#r-11-7)。

[講義ノート、定理 11.4](../../lecture_notes/ja/11_optimisation.md#r-11-4)（フォン・ノイマンのトレース不等式）により $|\mathrm{Tr}(AB)|\le\sum_is_i(A)s_i(B)$ であり、
[系 11.5](../../lecture_notes/ja/11_optimisation.md#r-11-5) により、固有値 $p,q$ をもつ自己共役な $P,Q$ に対し $\sum_ip^\downarrow_iq^\uparrow_i\le\mathrm{Tr}(PQ)\le\sum_ip^\downarrow_iq^\downarrow_i$、$P,Q\ge0$ に対し $\mathrm{Tr}(PQ)\le\|PQ\|_1\le\sum_ip^\downarrow_iq^\downarrow_i$ である。
スペクトルを固定した実験 2 については、[命題 11.6](../../lecture_notes/ja/11_optimisation.md#r-11-6) により $\sup z=-\inf z=\sum_c\sum_ix^\downarrow_{c,i}(y^\downarrow_{c,i}-y^\uparrow_{c,i})$ であり、これは $X_0=\mathrm{diag}(x^\downarrow)$、$V=J$（反転置換）、$Y_0=\mathrm{diag}(y^\uparrow)$
（それぞれ $\mathrm{diag}(y^\downarrow)$）で達成される。また[系 11.7](../../lecture_notes/ja/11_optimisation.md#r-11-7) は包絡線 $\sup z=-\inf z=\sum_c\lfloor k_c/2\rfloor(1-\ell_X)(1-\ell_Y)$ を与える。

ここで $\Delta(|X|)=\delta$ に対し $1-\ell_X=2\delta/(1+\delta)$ である。$k_c$ が奇数のとき中央の固有値は関係しない。コードは
$k=3$、$C=2$ について定理 11.4、系 11.5、および命題 11.6 の極値配置を検証する。

```python
M5 = TypeIAlgebra([5], [5], precision='double')
A = M5.operator(torch.randn(2000, 1, 5, 5, dtype=torch.complex128))         # Ginibre, no structure
B = M5.operator(torch.randn(2000, 1, 5, 5, dtype=torch.complex128))
sA, sB = torch.linalg.svdvals(A.matrix[:, 0]), torch.linalg.svdvals(B.matrix[:, 0])   # decreasing
assert ((A @ B).trace.abs() <= (sA * sB).sum(-1) + 1e-12).all()             # Theorem 11.4
Ph, Qh = (M5.operator_from_eigenvalues(lambda k: torch.randn(2000, k), batch_size=2000, force_self_adjoint=True)
          for _ in range(2))
p, q = Ph.eigenvalues()[0], Qh.eigenvalues()[0]                             # increasing order
t = (Ph @ Qh).trace.real
assert (t <= (p * q).sum(-1) + 1e-12).all() and (t >= (p * q.flip(-1)).sum(-1) - 1e-12).all()   # Cor. 11.5(a)
Pp, Qp = Ph @ Ph, Qh @ Qh                                                   # positive, spectra p^2, q^2
nuc = (Pp @ Qp).trace_norm()
assert ((Pp @ Qp).trace.real <= nuc + 1e-12).all()                          # Cor. 11.5(b)
assert (nuc <= ((p ** 2).sort(-1)[0] * (q ** 2).sort(-1)[0]).sum(-1) + 1e-10).all()

alg3 = TypeIAlgebra([3, 3], [3, 3], precision='double')
x_dn = torch.tensor([1.0, 0.6, ell(0.5)], dtype=torch.float64)              # Delta = 0.5, middle value 0.6
y_dn = torch.tensor([1.0, 0.7, ell(0.8)], dtype=torch.float64)              # Delta = 0.8, middle value 0.7
J = torch.eye(3, dtype=torch.float64).flip(0)                               # reversal permutation
blocks = lambda m: alg3.from_blocks([m, m])                                 # the same block in both sectors

def z_exp2(X0, V, Y0):
    return ((X0 @ V @ Y0).trace_norm() - (X0 @ Y0).trace.real).item()

X0 = blocks(torch.diag(x_dn))
z_max = z_exp2(X0, blocks(J), blocks(torch.diag(y_dn.flip(0))))
z_min = z_exp2(X0, blocks(J), blocks(torch.diag(y_dn)))
exact = 2 * (x_dn * (y_dn - y_dn.flip(0))).sum().item()                     # Proposition 11.6
envelope = 2 * (1 - ell(0.5)) * (1 - ell(0.8))                              # Corollary 11.7, floor(3/2) = 1
print(f"Exp. 2, k = 3, C = 2: z_max = {z_max:.6f}, z_min = {z_min:.6f}, formula {envelope:.6f}")
assert abs(z_max - exact) < 1e-12 and abs(z_min + exact) < 1e-12 and abs(exact - envelope) < 1e-12
```

<a id="s-11-3"></a>

## 11.3 等号条件つきのトレース不等式

*講義ノート：* [11.3節](../../lecture_notes/ja/11_optimisation.md#s-11-3)；[補題 11.8](../../lecture_notes/ja/11_optimisation.md#r-11-8)、[命題 11.9](../../lecture_notes/ja/11_optimisation.md#r-11-9)。

[講義ノート、命題 11.9](../../lecture_notes/ja/11_optimisation.md#r-11-9) により、実験 3 では $z=\mathrm{Tr}(Y|X|Y)-\mathrm{Tr}|YXY|\ge0$ であり、等号は $YX_+Y^2X_-Y=0$ のときかつそのときに限り
成り立つ。特に $XY=YX$ ならば等号が成り立つ。$X=\mathrm{diag}(1,-1,0)$ と、行が $(1,0,a)$, $(0,1,0)$, $(a,0,1)$, $|a|<1$ である $Y$ に対しては、$a\ne0$ ならば $XY\ne YX$ で
あるにもかかわらず等号が成り立つ。コードは両方を（$a=0.5$ で）検証する。`X.abs()` は $|X|$、`trace_norm()` は $\|\cdot\|_1$ である。

```python
def z_exp3(X, Y):
    return (Y @ X.abs() @ Y).trace.real - (Y @ X @ Y).trace_norm()

Xs = alg3.operator_from_eigenvalues(lambda k: torch.randn(500, k), batch_size=500, force_self_adjoint=True)
Yp = alg3.operator_from_eigenvalues(lambda k: torch.rand(500, k), batch_size=500, force_positive=True)
assert (z_exp3(Xs, Yp) >= -1e-12).all()                                     # Proposition 11.9

M3 = TypeIAlgebra([3], [3], precision='double')
Xe = M3.operator(torch.diag(torch.tensor([1.0, -1.0, 0.0]))[None, None])
Ye = M3.operator(torch.tensor([[1.0, 0, 0.5], [0, 1, 0], [0.5, 0, 1]])[None, None])
print(f"equality without commutation: z = {z_exp3(Xe, Ye).item():.1e}, "
      f"||XY - YX|| = {(Xe @ Ye - Ye @ Xe).operator_norm().item():.2f}")
assert abs(z_exp3(Xe, Ye).item()) < 1e-12 and (Xe @ Ye - Ye @ Xe).operator_norm().item() > 0.1
```

<a id="s-11-4"></a>

## 11.4 サンプリング対最適化

*講義ノート：* [11.4節](../../lecture_notes/ja/11_optimisation.md#s-11-4)；[補題 11.10](../../lecture_notes/ja/11_optimisation.md#r-11-10)、[命題 11.11](../../lecture_notes/ja/11_optimisation.md#r-11-11)、[系 11.12](../../lecture_notes/ja/11_optimisation.md#r-11-12)、[補題 11.13](../../lecture_notes/ja/11_optimisation.md#r-11-13)。

$|\mathrm{Tr}(AU)|$ は、その最大値 $\|A\|_1$ と比べて典型的にどれほどの大きさか。ライブラリは $U$ を、
**ジニブル行列** $Z$（成分は独立で密度
$\pi^{-1}e^{-|z|^2}$ をもつ。Mezzadri の処方）の QR 分解 $Z=QR$（$R_{ii}>0$）の因子 $Q$ として、各ブロックで独立にサンプリングする。

[講義ノート、命題 11.11](../../lecture_notes/ja/11_optimisation.md#r-11-11) により、ハール分布に従う $U\in U(k)$ に対し $\mathbb E|\mathrm{Tr}(AU)|^2=\|A\|_2^2/k$ であり、
[系 11.12](../../lecture_notes/ja/11_optimisation.md#r-11-12) により、スペクトルが $[a,b]$、$a>0$ に含まれる $X,Y\ge0$ に対し、$|\mathrm{Tr}(XUY)|$ の二乗平均平方根は
$U$ についての最大値の高々 $\min(1/\sqrt k,b^2/(a^2k))$ 倍であり、$\mathbb P\big(|\mathrm{Tr}(XUY)|\ge\mathrm{Tr}(XY)\big)\le b^4/(a^4k^2)$ である。

したがって、$N$ 個のハール・サンプルが $|\mathrm{Tr}(XUY)|\le\mathrm{Tr}(XY)$ の破れを見つける確率は高々
$Nb^4/(a^4k^2)$ である（$N$ 個の事象の和）。破れは可換でないすべての組に対して存在する（命題 3.9）にもかかわらず、である。
この上界が意味をもつのは $k>\sqrt N\,b^2/a^2$ のときだけである（$b/a=6$、$N=2000$ では $k>1609$）。数値的には、サンプルは
それよりはるかに小さい $k$ ですでに破れを見逃す（下記）。コードは命題 11.11 をモンテカルロ法で検証し（この
シード固定の実行では、標本平均は厳密な値から標準誤差の5倍以内にある）、$[0.2,1.2]$ 上で一様なスペクトル（$b/a=6$）に対して、
系 11.12 の比、$\sup z$、および $2000$ 個のサンプルのうち最良の $z$ を出力する。後二者はいずれも $\mathrm{Tr}(XY)$ に対する相対値である。

```python
def random_pair(alg_k):
    pos = lambda k: 0.2 + torch.rand(1, k)                                  # spectrum in [a, b] = [0.2, 1.2]
    return tuple(alg_k.operator_from_eigenvalues(pos, force_positive=True) for _ in range(2))

M6 = TypeIAlgebra([6], [6], precision='double')
Xf, Yf = random_pair(M6)
N = 20_000
t2 = (Yf @ Xf @ M6.random_unitary_operator(N)).trace.abs() ** 2
m2, err, exact2 = t2.mean().item(), t2.std().item() / math.sqrt(N), (Yf @ Xf).frobenius_norm().item() ** 2 / 6
print(f"E|Tr(AU)|^2 = {m2:.4f} +- {err:.4f}   (exact ||A||_2^2 / k = {exact2:.4f})")
assert abs(m2 - exact2) < 5 * err                                           # Proposition 11.11

print(" k   rms/max   sup z/TrXY   best of 2000 samples")
for k in (2, 4, 8, 16, 32):
    Mk = TypeIAlgebra([k], [k], precision='double')
    Xk, Yk = random_pair(Mk)
    trXY, nuc = (Xk @ Yk).trace.real.item(), (Yk @ Xk).trace_norm().item()
    ratio = (Yk @ Xk).frobenius_norm().item() / (math.sqrt(k) * nuc)
    best = (Xk @ Mk.random_unitary_operator(2000) @ Yk).trace.abs().max().item()
    print(f"{k:2d}   {ratio:.4f}    {(nuc - trXY) / trXY:.2e}     {(best - trXY) / trXY:+.3f}")
    assert ratio <= min(1 / math.sqrt(k), 36 / k) + 1e-12                   # Corollary 11.12
    assert nuc > trXY and best < nuc and (k < 4 or best < trXY)            # from k = 4 on: no violation found
```

*数値的観察*（出力される表。このシード固定の実行についてのみ assert する）：比は証明された上界
$\min(1/\sqrt k,b^2/(a^2k))$ を下回り続け、$k$ が2倍になるごとにほぼ半分になる。すなわち $1/k$ のように減少する。
また $k=4$ 以降では $2000$ 個のサンプルのうち最良のものでも $z=0$ にさえ届かない。これに対し勾配探索は、下記の実行で
両方の極値に近づく（これも数値的観察である）。組 `Xf, Yf` については、$U$ に関する両方の極値が分かっている。

[講義ノート、補題 11.13](../../lecture_notes/ja/11_optimisation.md#r-11-13) により、$k\ge2$ ならば任意の $A\in M_k(\mathbb C)$ に対し $\min_U|\mathrm{Tr}(AU)|=0$ である。

したがって $\inf_Uz=-\mathrm{Tr}(XY)$、$\sup_Uz=\|YX\|_1-\mathrm{Tr}(XY)$ である（命題 3.9）。`extremize` を一回
呼び出せば両方が近似される。$12$ 個のスタートが $z$ を最大化し、$12$ 個が $-z$ を最大化する（この実行では相対精度
$10^{-3}$ まで。下で assert する）。命題 11.11 により $|\mathrm{Tr}(XUY)|$ の二乗平均平方根は $\|YX\|_2/\sqrt k$ である。
系 11.12 の証明における評価（スペクトルは $[a,b]$ 内）により $\mathbb E|\mathrm{Tr}(XUY)|^2\le b^4$ かつ
$\mathrm{Tr}(XY)\ge ka^2$ であるから、マルコフの不等式により $\varepsilon>0$ に対し
$\mathbb P\big(z\ge-(1-\varepsilon)\mathrm{Tr}(XY)\big)\le b^4/(\varepsilon^2a^4k^2)$ が成り立つ：$k$ が大きいとき、
$z$ は高い確率でその下限に近い。
*数値的観察*（このシード固定の実行について assert する）：サンプリングは下限の $1\,\%$ 以内に達するが、
$z=0$ にさえ届かない。一方、上限は正である。

```python
trXY = (Xf @ Yf).trace.real
sup_exact, inf_exact = ((Yf @ Xf).trace_norm() - trXY).item(), -trXY.item()
U6 = opt.UnitaryParam(M6, 24)                                               # 24 starts
pm = torch.tensor([1.0, -1.0], dtype=torch.float64).repeat_interleave(12)   # 12 maximise z, 12 maximise -z
best = opt.extremize(lambda: pm * ((Xf @ U6.operator() @ Yf).trace.abs() - trXY), [U6],
                     maximize=True, rounds=6, steps=30).reshape(2, 12).max(1)[0]
best_up, best_down = best[0].item(), -best[1].item()
z_mc = (Xf @ M6.random_unitary_operator(N) @ Yf).trace.abs() - trXY
print(f"sup z: exact {sup_exact:.6e}, ascent {best_up:.6e}, best sample {z_mc.max().item():+.4f}")
print(f"inf z: exact {inf_exact:.6f}, descent {best_down:.6f}, worst sample {z_mc.min().item():+.6f}")
assert best_up <= sup_exact + 1e-12 and best_down >= inf_exact - 1e-12    # attained values: one-sided
assert sup_exact - best_up < 1e-3 * sup_exact and best_down - inf_exact < 1e-3 * abs(inf_exact)
assert z_mc.min().item() - inf_exact < 1e-2 * abs(inf_exact)               # sampling finds the typical side ...
assert z_mc.max().item() < 0                                               # ... but misses the violation
```

<a id="s-11-5"></a>

## 11.5 包絡線のためのマルチスタート探索

*講義ノート：* [11.5節](../../lecture_notes/ja/11_optimisation.md#s-11-5)；[命題 11.14](../../lecture_notes/ja/11_optimisation.md#r-11-14)、[系 11.7](../../lecture_notes/ja/11_optimisation.md#r-11-7)、[命題 11.9](../../lecture_notes/ja/11_optimisation.md#r-11-9)。

`optimize.extremize(objective, params, maximize)` は Adam を `steps` ステップずつ `rounds` ラウンド実行し（学習率は
`lr * decay**round`、損失は $\mp\sum_b\text{objective}_b$）、各ラウンドの後にユニタリ・チャートを再中心化し、
**スタートごとに**、勾配法の反復点と最終点（一度）で計算された目的関数の値のうち最大のもの（`maximize=False` のときは
最小のもの）を返す。バッチの各要素は、それぞれ固有のスロット、符号、$U_0$、コントラストをもつ独立なスタートである：
Adam は座標ごとに作用するから、目的関数の第 $b$ 成分が第 $b$ のパラメータだけに依存するならば（以下のすべての目的関数が
そうである）、スタート同士は相互作用しない。
`SelfAdjointParam` の符号は離散的である：ラウンドの間に、`discrete_step` はスタートごとに単一の符号反転を八回
（それぞれ一様にランダムなスロットで）、一つずつ順に試し、目的関数を改善する反転をそのつど保持する。

厳密な算術では、すべての反復は許容的である（命題 11.2 とその後の注意）から、返される値はすべて許容点における
$z$ の値である。浮動小数点ではそれは丸め誤差を除いてその値である（下の assert は $10^{-12}$ を許す）。したがって最大化は上限の
**下界** を、最小化は下限の上界を与え、それ以上のものではない——局所探索は局所極値で止まりうる。
逆向きの評価には、系 11.7 や命題 11.9 のような証明か、凸緩和（11.6節）が必要である。計算された包絡線は
(i) 閉じた形と照合し、(ii) 小さな解の直和のような手で構成した配置と照合し（包絡線はそれを下回ってはならない）、
(iii) 構造化されたサンプル（2準位のスペクトル、置換基底）と照合する。後者の値もまた達成される値である。

[講義ノート、命題 11.14](../../lecture_notes/ja/11_optimisation.md#r-11-14) により、すべての $k_c\ge2$ ならば、実験 1 の下側包絡線は $\underline z(\delta)=-\big(D-2\delta/(1+\delta)\big)$ である。

コードは $M_2\oplus M_2$ 上で命題 11.14 を、$M_3\oplus M_3$ 上で系 11.7 を（両方のコントラストを指定、
$k=3$ は奇数）、そして命題 11.9 を検証する。上限と下限は、上と同様に一回の呼び出しで求める。片側の assert（証明された
極値を越える値が返されないこと）はこれらの定理の帰結である。探索が極値から所定の許容誤差以内に達するという assert は、
このシード固定の実行の数値的観察である。キーワード `is_self_adjoint` は
トレースノルムのアルゴリズム（特異値、あるいは $YXY$ に対しては固有値）を、各ステップで判定することなく選択する。最後の
検証は型 (ii) である：実験 3 の $z$ は直和について加法的であり、直和のコントラストは各成分のコントラストの
最大値である（演習 3）から、コードと同様に $\Delta(Y)$ を $[0,0.999]$ 内で自由とした実験 3 の包絡線について、
$$
\overline z_{M_2\oplus M_2}(\delta)\ \ge\ \overline z_{M_2}(\delta)+\sup_{\delta'\le\delta}\overline z_{M_2}(\delta');
$$
特に、$M_2$ 上でコントラスト $\delta$ と $\delta'\le\delta$ において達成される二つの値の和は、$M_2\oplus M_2$ 上でコントラスト
$\delta$ において達成される値である。

```python
torch.manual_seed(1)
S = 6                                                                       # starts per grid point
grid = torch.tensor([0.0, 0.5, 0.9], dtype=torch.float64)
alg22 = TypeIAlgebra([2, 2], [2, 2], precision='double')                    # D = 4

# Exp. 1, lower envelope along Delta(X); Delta(Y) free in [0, 0.999]
X1, Y1, U1 = (opt.PositiveParam(alg22, 3 * S, delta=grid.repeat_interleave(S)),
              opt.PositiveParam(alg22, 3 * S), opt.UnitaryParam(alg22, 3 * S))
def z1():
    x, y = X1.operator(), Y1.operator()
    return (x @ U1.operator() @ y).trace.abs() - (x @ y).trace.real
inf1 = opt.extremize(z1, [X1, Y1, U1], maximize=False, rounds=4, steps=30).reshape(3, S).min(1)[0]
exact1 = -(4 - 2 * grid / (1 + grid))                                       # Proposition 11.14
print("Exp. 1  inf:", [f"{v:.4f}" for v in inf1.tolist()], " exact:", [f"{v:.4f}" for v in exact1.tolist()])
assert (inf1 >= exact1 - 1e-12).all() and (inf1 - exact1).max() < 4e-3

# Exp. 2 on M_3 (+) M_3, Delta(X0) = 0.5, Delta(Y0) = 0.8: sup and inf in one call
sign = torch.tensor([1.0, -1.0], dtype=torch.float64).repeat_interleave(S)
X2, Y2, V2 = (opt.PositiveParam(alg3, 2 * S, delta=torch.full((2 * S,), 0.5, dtype=torch.float64)),
              opt.PositiveParam(alg3, 2 * S, delta=torch.full((2 * S,), 0.8, dtype=torch.float64)),
              opt.UnitaryParam(alg3, 2 * S))
def z2():
    x0, y0 = X2.operator(), Y2.operator()
    M = alg3.operator((x0 @ V2.operator() @ y0).matrix, is_self_adjoint=False)
    return sign * (M.trace_norm() - (x0 @ y0).trace.real)
best2 = opt.extremize(z2, [X2, Y2, V2], maximize=True, rounds=4, steps=30).reshape(2, S).max(1)[0]
sup2, inf2 = best2[0].item(), -best2[1].item()
print(f"Exp. 2  sup {sup2:.6f}, inf {inf2:.6f}, exact +-{envelope:.6f}")
assert sup2 <= envelope + 1e-12 and inf2 >= -envelope - 1e-12               # Corollary 11.7: one-sided ...
assert envelope - sup2 < 1e-3 * envelope and inf2 + envelope < 1e-3 * envelope   # ... and reached

# Exp. 3: maximise sign * z at the contrasts `deltas` of |X|; Delta(Y) free in [0, 0.999]
def exp3(alg_, deltas, sign):
    X3, Y3 = opt.SelfAdjointParam(alg_, len(deltas), delta=deltas), opt.PositiveParam(alg_, len(deltas))
    def z3():
        y = Y3.operator()
        M = alg_.operator((y @ X3.operator() @ y).matrix, is_self_adjoint=True)
        return sign * ((y @ X3.abs_operator() @ y).trace.real - M.trace_norm())
    return opt.extremize(z3, [X3, Y3], maximize=True, rounds=4, steps=30).reshape(2, S).max(1)[0]

best3 = exp3(alg22, torch.full((2 * S,), 0.5, dtype=torch.float64), sign)
sup_M22, inf3 = best3[0].item(), -best3[1].item()
s_M2 = exp3(TypeIAlgebra([2], [2], precision='double'),
            torch.tensor([0.0, 0.5], dtype=torch.float64).repeat_interleave(S), 1.0)
direct_sum = (s_M2[1] + s_M2.max()).item()                                 # delta = 0.5, best delta' <= 0.5
conj = [(1 - dl) * (2 * 0.999 / 1.999) ** 2 for dl in (0.0, 0.5)]          # unproven closed form, d = 0.999
print(f"Exp. 3  inf {inf3:.1e};  sup for M_2 at delta = 0, 0.5: {s_M2[0]:.5f}, {s_M2[1]:.5f}"
      f"  (closed form {conj[0]:.5f}, {conj[1]:.5f})")
print(f"        sup for M_2 (+) M_2 at 0.5: search {sup_M22:.5f}, direct sum {direct_sum:.5f}")
assert inf3 >= -1e-12 and inf3 < 1e-6                                       # Proposition 11.9
assert sup_M22 >= direct_sum - 5e-3 * direct_sum                           # (ii), up to 0.5 %
```

*数値的観察*（この実行で出力される値）。$M_2\oplus M_2$ 上の探索は直和を約 $0.2\,\%$ 下回る：探索は完全には収束しておらず、
$\overline z_{M_2\oplus M_2}(0.5)$ の最良の下界は大きい方の数である——これを検出することが検証 (ii)
の目的である。$M_2$ に対する値は $(1-\delta)\big(2d/(1+d)\big)^2=(1-\delta)(1-\ell(d))^2$ と $10^{-4}$ まで一致する。ここで
$d=0.999$ は $Y$ に許される最大のコントラストである。この閉じた形は本講義の定理**ではない**：解説論文はこれを
（実環において）数値的に観察したが証明はしておらず、コードはこれを出力するだけで、何も assert しない。より大きなブロックに
ついては、実験 1 と実験 3 の上側包絡線の閉じた形を我々は知らない。その場合、計算された曲線は下界（達成される値）であって、
それ以上のものではない。

<a id="s-11-6"></a>

## 11.6 比較の対象

**Manopt**、**pymanopt**、**geoopt** は行列多様体（シュティーフェル多様体、直交群、ユニタリ群）を、
QR、極分解、ケーリーなどのレトラクションや、リーマン的ソルバー（共役勾配法、信頼領域法、リーマン
Adam）とともに提供しており、再中心化したチャートでの Adam に代わる手段となる（ここでは性能の比較は行わない）。
**cvxpy** は作用素変数の凸問題を SDP ソルバーに渡し、ソルバーは最適値を双対証明書とともに返す。どちらもソルバーの
許容誤差の範囲での話である。緩和——より大きな実行可能集合上の凸問題——は上限を上から抑える。これは局所探索では
得られない方向である。**NCTSSOS**、**ncpol2sdpa**、**NCSOStools** は、非可換
多項式の固有値とトレースをすべての次元で一度に評価する（0.3節）。ライブラリが加えるのは、$\bigoplus_cM_{k_c}$ のスペクトル制約と
コントラスト制約の厳密なパラメータ化（極値はどのセクターにもありうる）、ブロック構造を
保つユニタリ、そして多数のスタートと制約値に対する、CPU または GPU 上での一回のバッチ計算である。ライブラリは極値の
候補と上限の下界を生み出す。残りは証明か緩和が補わなければならない。

<a id="x-11"></a>

## 演習問題

講義ノートの[演習問題 1 と 2](../../lecture_notes/ja/11_optimisation.md#x-11)の数値的検証：歪エルミートな $S$ のケーリー変換 $C(S)=(1-S)(1+S)^{-1}$ はユニタリであり、
固有値 $-1$ をもたず、$S=(1-C)(1+C)^{-1}$ である（演習 1）。実環の `random_unitary_operator`（`complex_valued=False`）で
サンプリングした、$O(k)$ 上でハール分布に従う $O$ に対し、$\mathbb E\mathrm{Tr}(AO)^2=\|A\|_2^2/k$ である（演習 2）。

<details>
<summary>解答</summary>

```python
Sk = torch.randn(4, 4, dtype=torch.complex128)
Sk = Sk - Sk.mH                                                   # skew-adjoint
I4 = torch.eye(4, dtype=torch.complex128)
Cay = (I4 - Sk) @ torch.linalg.inv(I4 + Sk)
assert torch.allclose(Cay.mH @ Cay, I4, atol=1e-12)
assert torch.allclose((I4 - Cay) @ torch.linalg.inv(I4 + Cay), Sk, atol=1e-10)
assert torch.linalg.eigvals(Cay).add(1).abs().min() > 1e-6        # -1 is never an eigenvalue
```

```python
Mr = TypeIAlgebra([5], [5], complex_valued=False, precision='double')
Ar = Mr.operator(torch.randn(1, 1, 5, 5, dtype=torch.float64))
tr2 = (Ar @ Mr.random_unitary_operator(20_000)).trace ** 2
assert abs(tr2.mean().item() - Ar.frobenius_norm().item() ** 2 / 5) < 5 * tr2.std().item() / math.sqrt(20_000)
```

</details>

**次へ：** [チュートリアル 12](12_equivariant_learning.md)：可換子環の元としての同変層を、ライブラリで扱う。
