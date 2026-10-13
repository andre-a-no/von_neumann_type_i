[English](../09_spin_chains.md) · [Español](../es/09_spin_chains.md) · [Français](../fr/09_spin_chains.md) · [Deutsch](../de/09_spin_chains.md) · [中文](../zh/09_spin_chains.md) · **日本語** · [Русский](../ru/09_spin_chains.md)

# 9. スピン鎖

$L$ 個のスピン $\frac12$ からなる鎖は $(\mathbb C^2)^{\otimes L}$ 上に住み、次元 $2^L$ の密な線形代数は
$L\approx14$ で行き詰まる。対称性はこの限界を押し上げる：磁化は空間を $L+1$ 個のセクターに分け
（第4章）、並進は各セクターをさらに分け、$SU(2)$ はブロックをもう一度簡約する。本章では XXZ 鎖に対して
これらのセクターを構成し、それらが基底状態について何を語るかを証明し、自由フェルミオンとベーテ仮説によって
ライブラリを検証し、さらに密行列（クリロフ法）と乱れのない鎖（乱れ）の先へ進む。

**本章で学ぶこと**：鎖の環がなぜ $\bigoplus_NM_{\binom LN}$ であるか、これらの環の塔がなぜ
パスカルの三角形をブラッテリ図とし、指数 $4\cos^2\frac\pi{2L+4}$ をもつか、運動量セクターの数え方、
マーシャルの符号則と、それが基底状態の運動量とスピンについて導く帰結（リーブ–マティス）、XX 鎖の
ジョルダン–ウィグナーによる解、二マグノンのベーテ方程式、ランチョス法がなぜ単調に収束するか、
そして準位統計と自由フェルミオンのエンタングルメントの測り方。

本章を通じて、サイトは $0,\dots,L-1$ とし、$(\mathbb C^2)^{\otimes L}$ の基底は配置
$|b\rangle=|b_0\dots b_{L-1}\rangle$ からなり、上向きスピンに対して $b_i=1$ とする（ライブラリの規約：ラベルの
第 $i$ ビットがサイト $i$ である）。$s^a_i$（$a=x,y,z$）はサイト $i$ 上のスピン行列 $\frac12\sigma_a$、
$s^\pm_i=s^x_i\pm is^y_i$、$n_i=s^z_i+\frac12$、$S^a=\sum_is^a_i$ とし、$N=S^z+\frac L2$ は上向きスピンの
個数である。鎖は **開いた** 鎖（ボンド $(i,i+1)$、$i<L-1$）か **周期的** な鎖（ボンド $(L-1,0)$ も含む。このとき
$L\ge3$ で、添字は $L$ を法とする）である。$J>0$（反強磁性）、$\hbar=1$ とし、$\log$ は自然対数である。

## 9.1 XXZ 鎖とその $U(1)$ セクター

異方性 $\Delta\in\mathbb R$ と磁場 $h_i\in\mathbb R$ をもつ **XXZ 鎖** は
$$
H=J\sum_{(i,j)}\Big[\tfrac12\big(s^+_is^-_j+s^-_is^+_j\big)+\Delta\,s^z_is^z_j\Big]+\sum_ih_is^z_i ,
$$
である。和はボンドにわたってとる。$\frac12(s^+_is^-_j+s^-_is^+_j)=s^x_is^x_j+s^y_is^y_j$ は上向きスピンをボンドに沿って動かす。
$\Delta=1$ ではハイゼンベルク鎖 $J\sum\mathbf s_i\cdot\mathbf s_j$、$\Delta=0$ では XX 鎖である。

**命題 9.1.** $[H,S^z]=0$。したがって $H\in U(1)'\cong\bigoplus_{N=0}^LM_{\binom LN}(\mathbb C)$ であり、
ブロック $N$ は上向きスピンを $N$ 個もつ配置の張る空間に作用する。

*証明.* 対角項は対角な $S^z$ と可換である。$i\neq j$ のとき、$s^+_is^-_j$ は $l\ne i,j$ に対する
$s^z_l$ と可換であり、$[s^z,s^\pm]=\pm s^\pm$（4.5節）により
$[s^z_i+s^z_j,s^+_is^-_j]=[s^z_i,s^+_i]s^-_j+s^+_i[s^z_j,s^-_j]=s^+_is^-_j-s^+_is^-_j=0$ である。
$s^-_is^+_j$ については共役をとればよい。残りは命題 4.10(2) である。$\square$

`SpinChain(L).xxz` は各セクターの配置（`chain.labels[N]`）から直接ブロックを構成する。
これを、クロネッカー積から構成した完全な $2^L\times2^L$ 行列と比較する。

```python
import functools, itertools, math, time
import numpy as np
import torch
from scipy.optimize import brentq
from torch_vn_algebra import SpinChain, SU2Chain, DensityMatrix, inclusions, krylov

torch.manual_seed(0)
f64, c128 = torch.float64, torch.complex128
close = lambda a, b, tol=1e-10: bool((torch.as_tensor(a) - torch.as_tensor(b)).abs().max() < tol)
def op(L, ops):                             # product of one-site operators {site: 2 x 2}; site 0 = lowest bit
    return functools.reduce(torch.kron, [ops.get(j, torch.eye(2, dtype=f64)) for j in reversed(range(L))])

sp = torch.tensor([[0., 0.], [1., 0.]], dtype=f64)                    # s^+ : |0> = down -> |1> = up
sm, sz = sp.T.contiguous(), torch.diag(torch.tensor([-0.5, 0.5], dtype=f64))
L = 6
chain = SpinChain(L, boundary='periodic')
h = 0.3 * torch.randn(L, dtype=f64)                                  # random fields
H_full = sum(0.5 * (op(L, {i: sp, j: sm}) + op(L, {i: sm, j: sp})) + 0.7 * op(L, {i: sz, j: sz})
             for i, j in chain.bonds) + sum(h[i] * op(L, {i: sz}) for i in range(L))   # J = 1, Delta = 0.7
Sz = sum(op(L, {i: sz}) for i in range(L))
assert close(H_full @ Sz, Sz @ H_full)                               # Proposition 9.1
H = chain.xxz(J=1.0, Delta=0.7, h=h)                                 # an element of (+)_N M_binom(6, N)
for N, k in enumerate(chain.algebra.k_factors):
    lab = chain.labels[N]                                            # configurations of sector N
    assert close(H.matrix[0, N, :k, :k], H_full[lab][:, lab])
```

**鎖の塔.** $\mathcal A_L=\bigoplus_NM_{\binom LN}$ を $L$ 個のサイトに作用させる。サイトを一つ加える写像
$x\mapsto x\otimes1_2$ は $\mathcal A_L$ を $\mathcal A_{L+1}$ の中に写す。$x\otimes1$ は
$S^z\otimes1+1\otimes s^z$ と可換だからである。これにより塔 $\mathcal A_0\subset\mathcal A_1\subset\mathcal A_2\subset\cdots$ が得られる。

**命題 9.2.** 包含 $\mathcal A_L\subset\mathcal A_{L+1}$ の包含行列は $(L+2)\times(L+1)$
行列 $\Lambda_L$ であり、$\lambda_{N,N}=\lambda_{N,N-1}=1$、他の成分はすべて $0$ である。ブロックの
大きさはパスカルの法則 $\binom{L+1}N=\binom LN+\binom L{N-1}$ に従い、塔のブラッテリ図は
パスカルの三角形であり、$\Lambda_L$ の二部グラフは $2L+3$ 個の頂点をもつ道である。

*証明.* $L+1$ 個のサイト上で上向きスピンを $N$ 個もつ配置は、最初の $L$ 個のサイトに上向きスピンが
$N$ 個あって新しいスピンが下向きのものと、$N-1$ 個あって新しいスピンが上向きのものである。この二つの部分空間の上で
$x\otimes1$ はそれぞれ $x_N$ と $x_{N-1}$ として作用するので、$\mathcal A_{L+1}$ のブロック $N$ は
$\mathcal A_L$ のブロック $N$ と $N-1$ を一度ずつ含む。次元を足し合わせるとパスカルの法則になる。二部グラフでは
$\mathcal A_{L+1}$ のセクター $N$ は $\mathcal A_L$ の $N-1$ と $N$ に結ばれ、これにより道
$0'-0-1'-1-\cdots-L'-L-(L+1)'$ が得られる（プライムは $\mathcal A_{L+1}$ のもの）。$\square$

ブラッテリ図の頂上からレベル $L$ のセクター $N$ に至る道は一つの配置 $b$ であり、その数は
$\binom LN$、すなわちブロックの大きさである。包含行列 $\Lambda$ をもつ連結な包含に対して、
$\|\Lambda\|^2$ を **指数** とする。これはマルコフ・トレースの指数である（V. Jones, *Invent. Math.* 72 (1983); F. Goodman,
P. de la Harpe, V. Jones, *Coxeter graphs and towers of algebras*, Springer 1989, Ch. 2–3）。
$M_k\otimes1\subset M_k\otimes M_m$ に対しては、これは 7.2節の $m^2$ である。`Inclusion.index` がこれを計算する。

**補題 9.3.** 頂点 $1,\dots,n$ をもつ道の隣接行列 $A_n$ は固有値
$2\cos\frac{\pi k}{n+1}$、$k=1,\dots,n$ と、固有ベクトル $v_k(j)=\sin\frac{\pi jk}{n+1}$ をもつ。したがって
$\|A_n\|=2\cos\frac\pi{n+1}$ である。

*証明.* $v_k(0)=v_k(n+1)=0$ とおく。$\sin(a-b)+\sin(a+b)=2\sin a\cos b$ により
$(A_nv_k)(j)=v_k(j-1)+v_k(j+1)=2\cos\frac{\pi k}{n+1}\,v_k(j)$ であり、$v_k(1)=\sin\frac{\pi k}{n+1}\neq0$ であるから $v_k\ne0$ である。
$n$ 個の固有値は相異なる（$\cos$ は $(0,\pi)$ 上で狭義単調減少）ので、これらで固有値はすべてであり、
対称行列 $A_n$ のノルムは固有値の絶対値の最大値で、$k=1$ と $k=n$ で達成される。$\square$

**系 9.4.** $[\mathcal A_{L+1}:\mathcal A_L]=\|\Lambda_L\|^2=4\cos^2\frac\pi{2L+4}$ であり、これは増加して $4$ に近づく。

*証明.* 頂点を命題 9.2 の道に沿って並べると、隣接行列
$B=\begin{pmatrix}0&\Lambda_L\\\Lambda_L^T&0\end{pmatrix}$ は $A_{2L+3}$ になる（置換による相似）。補題 1.1 により
$\|B\|^2=\|B^2\|$ かつ $\|\Lambda_L\Lambda_L^T\|=\|\Lambda_L^T\Lambda_L\|=\|\Lambda_L\|^2$ であり、$B^2=\Lambda_L\Lambda_L^T\oplus\Lambda_L^T\Lambda_L$ であるから
$\|B\|=\|\Lambda_L\|$ である。$n=2L+3$ とした補題 9.3 が値を与える。$\square$

数 $4\cos^2\frac\pi n$（$n\ge3$）は、部分因子の指数がとりうる $4$ 未満の値にちょうど一致する
（Jones, 前掲論文）。ここでは偶数の $n$ が、補題 9.3 という初等的な理由によって現れる。

```python
tower = inclusions.chain_tower(8)                                    # A_0 c A_1 c ... c A_8
print(inclusions.bratteli_diagram(tower).levels[:6])                 # Pascal's triangle
for L_, inc in enumerate(tower):
    sizes = [math.comb(L_ + 1, N) for N in range(L_ + 2)]
    assert inc.M.k_factors == sizes and inc.is_connected
    assert list(inc.matrix @ np.array(inc.N.k_factors)) == sizes      # unitality = Pascal's rule
    n = 2 * L_ + 3
    A_path = np.eye(n, k=1) + np.eye(n, k=-1)
    assert abs(np.linalg.eigvalsh(A_path).max() - 2 * math.cos(math.pi / (n + 1))) < 1e-12   # Lemma 9.3
    assert abs(inc.index() - 4 * math.cos(math.pi / (2 * L_ + 4)) ** 2) < 1e-12           # Corollary 9.4
```

> **物理的意味.** 相互作用が $z$ 軸まわりの回転で不変であれば、$S^z$ は保存される。$N$ は上向きスピンの数、
> あるいは 9.4節のジョルダン–ウィグナー写像の後では粒子の数を数える。塔の中では、長い鎖の各
> 電荷セクターに短い鎖の二つのセクターが流れ込む。

## 9.2 並進、運動量、マーシャルの符号則

周期的な鎖の上で、$T$ を並進 $T|b_0b_1\dots b_{L-1}\rangle=|b_{L-1}b_0\dots b_{L-2}\rangle$ とする。
これはサイト $i$ をサイト $i+1$ に移す。$T$ は置換行列であるからユニタリであり、$T^L=1$ かつ
$Ts^a_iT^*=s^a_{i+1}$ である。したがって $TS^zT^*=S^z$ であり、一様な磁場 $h_i=h$ に対しては $THT^*=H$ である：$T$ はボンドをボンドに写す。
その固有値は $e^{ik}$ で、**運動量** は $k=2\pi m/L$、$m=0,\dots,L-1$ である。

**命題 9.5（運動量セクター）.** $N$ を固定する。上向きスピンを $N$ 個もつ配置 $a$ に対し、$R_a$ を
その軌道 $\{T^ja\}$ の大きさとする。(a) $R_a$ は $L$ を割り切る。$mR_a\equiv0\pmod L$ を満たす $k=2\pi m/L$ に対し、ベクトル
$|a,k\rangle=\sum_{j=0}^{R_a-1}e^{-ikj}T^j|a\rangle$ は零でなく $T|a,k\rangle=e^{ik}|a,k\rangle$ を満たし、これら
$R_a$ 個のベクトルは軌道の張る空間の直交基底をなす。(b) セクター $N$ 上での $T$ の固有値 $e^{2\pi im/L}$ の
重複度は $d_{N,m}=\#\{\text{軌道 }:mR_a\equiv0\pmod L\}$ であり、
$S^z$ および $T$ と可換な任意の作用素——たとえば並進不変な XXZ 鎖——はブロック対角であり、$\bigoplus_N\bigoplus_mM_{d_{N,m}}$ に属する。

*証明.* (a) 巡回群 $\{T^j\}\cong\mathbb Z_L$ における $a$ の安定化群は部分群であり、$R_a$ はその
指数であるから $L$ の約数である。ベクトル $T^j|a\rangle$（$0\le j<R_a$）は正規直交であり、$T^{R_a}|a\rangle=|a\rangle$ である。
$e^{ikR_a}=1$ ならば、$j=R_a$ の項は $j=0$ の項に等しいので
$T|a,k\rangle=e^{ik}\sum_{j=1}^{R_a}e^{-ikj}T^j|a\rangle=e^{ik}|a,k\rangle$ である。$|a\rangle$ の係数は $1$ である。許される
$m\in\{0,\dots,L-1\}$ は $L/R_a$ の倍数で、$R_a$ 個あり、ユニタリな $T$ の相異なる固有値に対応するので、これらのベクトルは直交し、
軌道の張る $R_a$ 次元の空間を張る。(b) セクター $N$ はその軌道の張る空間の直交和である。
$S^z$ および $T$ と可換な $X$ はそれらのスペクトル射影 $P_N$、$Q_m$ と可換である（定理 1.5(d)）から、
$X=\sum_{N,m}P_NQ_mXP_NQ_m$ である。$\square$

`xxz_momentum` は、各軌道から代表元 $a$ を一つずつとった正規化された基底 $|a,k\rangle$ でブロックを構成する。
軌道を独立に数え、スペクトルを密なセクターと比較する。

```python
L = 12
ring = SpinChain(L, boundary='periodic')
rot = lambda x: ((x << 1) | (x >> (L - 1))) & ((1 << L) - 1)        # T on labels: site i -> i + 1
for N in range(L + 1):
    orbits = {frozenset(itertools.accumulate(range(L - 1), lambda x, _: rot(x), initial=a))
              for a in ring.labels[N].tolist()}                       # {a, Ta, ..., T^{L-1} a}
    d = [sum((m * len(o)) % L == 0 for o in orbits) for m in range(L)]
    assert d == ring.momentum_dims(N)                                   # Proposition 9.5(b)
print(ring.momentum_dims(6))
Hk = ring.xxz_momentum(J=1.0, Delta=1.0, h=0.2, sector=6)               # 12 blocks instead of one of size 924
E_k = torch.sort(torch.cat([e[0] for e in Hk.eigenvalues()])).values
assert close(E_k, torch.linalg.eigvalsh(ring.xxz(1.0, 1.0, h=0.2, sector=6).matrix[0, 0]))
```

$924$ 個の状態からなるセクター $N=6$ は、大きさ $80,75,78,76,78,75,80,75,78,76,78,75$ のブロックに分かれる。

**補題 9.6（ペロン–フロベニウス、対称な場合）.** $M$ を、$i\ne j$ に対して $M_{ij}\ge0$ であり、
そのグラフ（$i\sim j$ であるのは $M_{ij}>0$ のとき、かつそのときに限る）が連結な実対称行列とする。このとき最大固有値 $\lambda$ は単純であり、
すべての成分が $>0$ である固有ベクトルをもつ。

*証明.* $\lambda=\max_{\|x\|=1}\langle x,Mx\rangle$ であり、最大値は
$\lambda$ 固有空間の単位ベクトルでちょうど達成される（正規直交固有基底で展開せよ）。そのような $x$ に対し、絶対値のベクトル $|x|$ は
$\langle|x|,M|x|\rangle\ge\operatorname{Re}\langle x,Mx\rangle=\lambda$ を満たす。対角項は一致し、$i\ne j$ に対して
$M_{ij}|x_i||x_j|\ge\operatorname{Re}M_{ij}\bar x_ix_j$ だからである。よって $|x|$ は固有ベクトルである。$|x|_i=0$ ならば
$0=(\lambda-M_{ii})|x|_i=\sum_{j\ne i}M_{ij}|x|_j$ から $i$ のすべての隣接点 $j$ について $|x|_j=0$ が従い、
連結性により $x=0$ となる。したがって零でない固有ベクトルは零成分をもたない。しかし次元 $\ge2$ の固有空間は
$x_1=0$ を満たす零でないベクトルを含む。よって $\lambda$ は単純であり、$|x|>0$ は固有ベクトルである。$\square$

**定理 9.7（マーシャルの符号則）.** 鎖は開いているか、$L$ が偶数の周期的な鎖であるとする。$B$ を
奇数番目のサイトの集合、$N_B(b)$ を $b$ の $B$ 上の上向きスピンの数とし、$J>0$、$\Delta,h_i\in\mathbb R$ は任意とする。
各セクター $N$ において $H_N$ の最小固有値は単純であり、その固有ベクトルは、上向きスピンを $N$ 個もつ
すべての配置 $b$ に対して $(-1)^{N_B(b)}\psi(b)>0$ となるように選べる。

*証明.* $V|b\rangle=(-1)^{N_B(b)}|b\rangle$ とおく。$V$ は対角であり、$V=V^*=V^{-1}$ である。各ボンドは偶数サイトと
奇数サイトを結ぶ（環については $L$ が偶数だから）ので、$s^+_is^-_j$ は $N_B$ を $\pm1$ だけ変え、
$Vs^+_is^-_jV=-s^+_is^-_j$ である。よって $M=-VH_NV$ は、一つのボンドの二つのスピンの交換だけ異なる配置の間で
非対角成分 $J/2>0$ をもち、それ以外では $0$ である。そのグラフは連結である：隣接サイト $(i,i+1)$（$i<L-1$）の
交換はサイトのすべての置換を生成するので、上向きスピンを $N$ 個もつ任意の二つの配置を結ぶ
（等しいスピンの交換は何もしない）。補題 9.6 により、$M$ の最大固有値——$H_N$ の最小固有値の符号を変えたもの——は
単純で、固有ベクトル $\varphi>0$ をもつ。$\psi=V\varphi$ とすればよい。$\square$

**系 9.8（基底状態の運動量）.** $L$ が偶数、$J>0$ で磁場が一様な周期的な鎖の上で、
セクター $N$ の基底状態の運動量は、$N$ が偶数なら $k=0$、$N$ が奇数なら $k=\pi$ である。

*証明.* $T$ は $H_N$ と可換であり、基底状態 $\psi=V\varphi$ は単純であるから、$T\psi=e^{ik}\psi$ である。ここで
$(T\psi)(b)=\psi(T^{-1}b)$ であり、$T^{-1}b$ はサイト $i$ にスピン $b_{i+1}$ をもつ。$L$ は偶数なので、奇数の $i$ に対して $i+1$ は偶数であり、
$N_B(T^{-1}b)=N-N_B(b)$、したがって $(VT\psi)(b)=(-1)^{N_B(b)}(-1)^{N-N_B(b)}\varphi(T^{-1}b)=(-1)^N\varphi(T^{-1}b)$ である。また
$VT\psi=e^{ik}\varphi$ である。$\varphi$ と $\varphi\circ T^{-1}$ はともに正であるから、$e^{ik}=(-1)^N$ である。$\square$

```python
for bc in ('open', 'periodic'):
    c = SpinChain(10, bc, complex_valued=False)
    for N in range(11):
        w, V = torch.linalg.eigh(c.xxz(1.0, 0.7, h=0.3 * torch.randn(10, dtype=f64), sector=N).matrix[0, 0])
        x = (-1.0) ** c.bits(N)[:, 1::2].sum(-1) * V[:, 0]                  # (-1)^{N_B(b)} psi(b)
        assert (x * x.sum().sign()).min() > 0 and (len(w) == 1 or w[1] - w[0] > 1e-8)   # Theorem 9.7
for N in range(1, L):
    lowest = [e[0, 0].item() for e in ring.xxz_momentum(1.0, 0.5, sector=N).eigenvalues()]
    m0 = int(np.argmin(lowest))
    assert m0 == (0 if N % 2 == 0 else L // 2)                              # Corollary 9.8: k = 0 or pi
    assert sorted(lowest)[1] - lowest[m0] > 1e-8
```

> **物理的意味.** 反強磁性体では、格子が基底状態の符号を決める：スピンをボンドの向こうへ動かすと
> 符号が変わる。このため二部格子の反強磁性体に対する量子モンテカルロ法には負符号問題がない。
> 結合が同じ副格子のサイトどうしを結ぶフラストレートした鎖では、これは成り立たない。

## 9.3 $\Delta=1$ における $SU(2)$

$\Delta=1$ かつ $h=0$ のとき、$H=J\sum\mathbf s_i\cdot\mathbf s_j$ は $S^x,S^y,S^z$ と可換である（命題 4.15）から、
$H\cong\bigoplus_S1_{2S+1}\otimes H_S$ であり、ブロックは $H_S\in M_{m_S}$、$m_S=\binom L{L/2-S}-\binom L{L/2-S-1}$
（命題 4.13）である。モジュール `su2` は最高ウェイトベクトルの上で $H_S$ を構成する（命題 4.15(b)）。
$E_0(S)$ を $H_S$ の最小固有値とする。

**定理 9.9（鎖に対するリーブ–マティスの順序）.** 鎖は定理 9.7 と同様とし、$\Delta=1$、$h_i=0$、
$J>0$ とする。$L$ 個のスピンのとりうる各全スピン $S$ に対し、セクター $S^z=S$ の基底状態は全スピン $S$ をもつ。したがって
$E_0(S)$ は $H_S$ の単純な固有値であり、$S_{\min}\in\{0,\frac12\}$ として
$E_0(S_{\min})<E_0(S_{\min}+1)<\dots<E_0(\frac L2)$ が成り立つ：基底状態は最小スピンの単一の多重項である。

*証明.* $\psi$ をセクター $S^z=S$ の基底状態とする。これは定理 9.7 により単純である。$\mathbf S^2$ は
$H$ と可換でセクターを保つので、$\mathbf S^2\psi$ は一次元の基底状態空間に属する：$\psi$ はある全
スピン $S'$ をもつ。$\chi$ を、対 $(0,1),(2,3),\dots,(L-2S-2,L-2S-1)$ の上の一重項
$\frac1{\sqrt2}(|{\uparrow\downarrow}\rangle-|{\downarrow\uparrow}\rangle)$ の積で、残りの $2S$ 個のサイトが上向きのものとする。$S^+$ は各因子を消し、
$S^z\chi=S\chi$ であるから、$\mathbf S^2\chi=(S^-S^++(S^z)^2+S^z)\chi=S(S+1)\chi$ である（4.5節）。各対において上向きスピンは
偶数サイト（係数 $+$）か奇数サイト（$-$）にあるので、$(-1)^{N_B(b)}\chi(b)$ はその台に含まれるすべての $b$ に対して
同じ符号をもつ。定理 9.7 の $V,\varphi$ を用いれば $\langle\psi,\chi\rangle=\langle\varphi,V\chi\rangle\neq0$ である。
$\mathbf S^2$ の相異なる固有値に属する固有ベクトルは直交するので、$S'=S$ である。命題 4.15(a) と系 4.12 により、
スピン $S'\ge S$ の各多重項はセクター $S^z=S$ とちょうど一つのベクトルで交わるので、そこでの $H$ のスペクトルは
$H_{S'}$（$S'\ge S$）のスペクトルの（重複度込みの）和集合である。その最小固有値は単純で、固有ベクトルは
スピン $S$ をもつから、それは $E_0(S)$ であり、$H_S$ において単純で、$H_{S'}$（$S'>S$）の固有値ではない：$E_0(S)<E_0(S')$。$\square$

リーブとマティスはさらに多くを証明した：副格子スピンが $S_A,S_B$ の連結な二部格子で、結合が
副格子間では反強磁性的、副格子内では強磁性的であれば、$E_0(S)$ は $S\ge|S_A-S_B|$ で増加し、
基底状態はスピン $|S_A-S_B|$ をもつ（E. Lieb, D. Mattis, *J. Math. Phys.* 3 (1962) 749。証明なしで
述べる）。**一重項–三重項ギャップ** $E_0(1)-E_0(0)$ は定理 9.9 により正である。環の上では
$L(E_0(1)-E_0(0))\to\pi^2/2$ であり、対数補正を伴う（J. des Cloizeaux, J. Pearson, *Phys. Rev.* 128
(1962) 2131; I. Affleck, D. Gepner, H. Schulz, T. Ziman, *J. Phys. A* 22 (1989) 511。証明なしで述べる）。
十二個のスピンについてブロックの大きさを比較し、定理 9.9 を $U(1)$ セクターで確かめ、ギャップを計算する。

```python
sym = SU2Chain(12, boundary='periodic', complex_valued=False)
H_su2 = sym.heisenberg(J=1.0)                                       # (+)_S M_{m_S}, S = 0, ..., 6
u1_sizes = [math.comb(12, N) for N in range(13)]
print(sym.multiplicities, max(u1_sizes), sum(k ** 3 for k in u1_sizes), sum(m ** 3 for m in sym.multiplicities))
E0_S = sym.lowest_energies(H_su2)[0]                                # E_0(S), S = 0, ..., 6
E0_N = [torch.linalg.eigvalsh(ring.xxz(1.0, 1.0, sector=6 + S).matrix[0, 0])[0] for S in range(7)]
assert close(E0_S, torch.stack(E0_N))                               # Theorem 9.9: sector S^z = S has spin S
assert (E0_S[1:] > E0_S[:-1] + 1e-8).all()                          # E_0(S) strictly increasing
gap = {}
for n in (6, 8, 10, 12):
    s = SU2Chain(n, 'periodic', complex_valued=False)
    e = [torch.linalg.eigvalsh(s.heisenberg(spin=S).matrix[0, 0])[0].item() for S in (0, 1)]
    gap[n] = n * (e[1] - e[0])
print({n: round(g, 4) for n, g in gap.items()}, round(math.pi ** 2 / 2, 4))
assert 0 < gap[6] < gap[8] < gap[10] < gap[12] < math.pi ** 2 / 2
```

$SU(2)$ ブロックの大きさは $m_S=132,297,275,154,54,11,1$ で、最大の $U(1)$ ブロックの $924$ と比べて小さく、
対角化のコスト $\sum k^3$ は $2.0\cdot10^9$ から $5.3\cdot10^7$ に減る。$L$ とギャップの積はゆっくりと増加し、
$L=6,\dots,12$ に対して $4.108,4.181,4.232,4.270$ と、$\pi^2/2\approx4.935$ に向かう。

**相関.** 各サイトのスピンは階数一のテンソル演算子をなし、**ウィグナー–エッカートの定理** により、
多重項間の行列要素はクレブシュ–ゴルダン係数と、磁気量子数によらない **換算行列要素**
$\langle S',a\|\mathbf s_i\|S,b\rangle$ との積である（A. R. Edmonds, *Angular
Momentum in Quantum Mechanics*, Princeton 1957, Ch. 5。証明なしで述べる）。`spin_correlation` はこれを用いて、
最高ウェイトベクトルだけから各 $H_S$ 上の $\mathbf s_i\cdot\mathbf s_j$ を構成する。`'direct'` 法は
代わりにボンド演算子を射影する。両者を、$U(1)$ セクターで書き下した基底状態と比較する。

```python
gs = sym.ground_state(H_su2)
S0, coef = gs.spin.item(), gs.coefficients[0, :sym.multiplicity(gs.spin.item())]
pairs = [(0, r) for r in range(12)]
corr_we = sym.expectation(sym.spin_correlation(pairs), S0, coef)
corr_dir = sym.expectation(sym.spin_correlation(pairs, method='direct'), S0, coef)
psi = sym.embed(S0, coef)                                            # the ground state in the sector N = 6
lab, z = sym.u1.labels[6], sym.u1.bits(6).to(f64) - 0.5
def s_dot_s(i, j):                                                   # <psi| s_i . s_j |psi>, from the definition
    if i == j:
        return 0.75
    flip = ((lab >> i) & 1) != ((lab >> j) & 1)
    partner = sym.u1.index_of[lab[flip] ^ ((1 << i) | (1 << j))]
    return ((psi ** 2 * z[:, i] * z[:, j]).sum() + 0.5 * (psi[flip] * psi[partner]).sum()).item()
assert S0 == 0 and close(corr_we, corr_dir)
assert close(corr_we, torch.tensor([s_dot_s(0, r) for r in range(12)], dtype=f64))
assert ((-1.0) ** torch.arange(1, 12) * corr_we[1:] > 0).all()       # signs (-1)^r
print([round(x, 4) for x in corr_we.tolist()])
```

相関 $\langle\mathbf s_0\cdot\mathbf s_r\rangle=0.75,-0.4489,0.1878,-0.1658,0.1209,-0.1224,0.1070,\dots$ は
定理 9.7 と 9.9 が強制するとおり符号が交互に変わる：一重項は回転で不変であり、回転は座標軸を
入れ替える。またセクター $S^z=0$ では $\langle s^+_0s^+_r\rangle=0$ であるから、$\langle\mathbf s_0\cdot\mathbf s_r\rangle=3\langle s^x_0s^x_r\rangle=
\frac34\langle s^+_0s^-_r+s^-_0s^+_r\rangle=\frac34\sum_{b_0\ne b_r}\psi(b)\psi(b')$ である。ここで $b'$ は $b$ からサイト $0$ と $r$ の
スピンを交換して得られる。$N_B(b')-N_B(b)$ は奇数の $r$ では $\pm1$、偶数の $r$ では $0$ であるから、各項は符号
$(-1)^r$ をもつ。減衰は遅く、$L\to\infty$ で $(-1)^r\sqrt{\log r}/r$ のようにふるまう（Affleck ら, 前掲論文。証明
なしで述べる）。

## 9.4 厳密な検証：自由フェルミオンとベーテ仮説

**ジョルダン–ウィグナー.** 開いた鎖の上で $c_j=\big(\prod_{l<j}(-1)^{n_l}\big)s^-_j$ とおく。

**定理 9.10（XX 鎖は自由である）.** (a) $\{c_i,c_j^*\}=\delta_{ij}$、$\{c_i,c_j\}=0$、$c_j^*c_j=n_j$。
(b) 開いた XX 鎖は $H=\frac J2\sum_{j<L-1}(s^+_js^-_{j+1}+s^-_js^+_{j+1})=\sum_{i,j}h_{ij}c_i^*c_j$ であり、
$h=\frac J2A_L$ である。(c) $u_k$ を $h$ の正規直交な実固有ベクトルで固有値が
$\varepsilon_k=J\cos\frac{\pi k}{L+1}$ のもの（補題 9.3）とし、$b_k=\sum_ju_k(j)c_j$ とおく。$K=\{k_1<\dots<k_N\}$ に対し、ベクトル
$b_K^*|0\rangle=b_{k_1}^*\cdots b_{k_N}^*|0\rangle$（$|0\rangle$ はすべてのスピンが下向きの状態）はセクター
$N$ の正規直交固有基底をなし、エネルギーは $\sum_{k\in K}\varepsilon_k$ である。

*証明.* (a) 異なるサイト上の作用素は可換であり、一つのサイト上では $(-1)^n$ は $s^\pm$ と反可換である。
$i<j$ のとき、$c_i$ と $c_j^*$（または $c_j$）はサイト $i$ で $s^-_i$ と $(-1)^{n_i}$ として出会い——一回の反交換——
それ以外では可換であるから、反可換である。$i=j$ のとき、ストリングの二乗は $1$ であり、
$c_jc_j^*+c_j^*c_j=s^-_js^+_j+s^+_js^-_j=1$、$c_j^*c_j=s^+_js^-_j=n_j$ である。(b) $s^+_j\neq0$ となるのは $n_j=0$ の上だけであるから
$c_j^*c_{j+1}=s^+_j(-1)^{n_j}s^-_{j+1}=s^+_js^-_{j+1}$ である。共役をとればよい。(c) $\{b_k,b_l^*\}=\sum_ju_k(j)u_l(j)=\delta_{kl}$ かつ
$\{b_k,b_l\}=0$ であり、$h$ のスペクトル分解により $H=\sum_k\varepsilon_kb_k^*b_k$ である。正準反交換関係（CAR）から
$[b_l^*b_l,b_k^*]=\delta_{lk}b_k^*$、よって $[H,b_k^*]=\varepsilon_kb_k^*$ である。$H|0\rangle=0$ と合わせて
$Hb_K^*|0\rangle=(\sum_K\varepsilon_k)b_K^*|0\rangle$ を得る。$b_k|0\rangle=0$（各 $c_j$ はスピンを下げる）であるから、CAR を用いて消滅演算子を
右へ移すと、$N$ についての帰納法により $\langle b_{K'}^*0,b_K^*0\rangle=\delta_{KK'}$ を得る。各 $c_j^*$ は
$N$ を一つ上げるので、$b_K^*|0\rangle$ はセクター $|K|$ に属する：その次元の空間の中の $\binom LN$ 個の正規直交ベクトルである。$\square$

```python
L = 10
xx = SpinChain(L)
eps = [math.cos(math.pi * k / (L + 1)) for k in range(1, L + 1)]
for N in range(L + 1):
    E = torch.linalg.eigvalsh(xx.xxz(J=1.0, Delta=0.0, sector=N).matrix[0, 0])
    E_ff = sorted(sum(K) for K in itertools.combinations(eps, N))     # sums over N-element subsets
    assert close(E, torch.tensor(E_ff, dtype=f64))                    # all 1024 levels, Theorem 9.10
    assert abs(E[0].item() - sum(sorted(eps)[:N])) < 1e-12           # ground state: the N lowest modes
```

**マグノン.** $h=0$ の環の上で、完全に偏極した状態 $|0\rangle$ のエネルギーは $E_F=J\Delta L/4$ である。$x$ にある一つの上向き
スピンは二つのボンドを壊す：$H|x\rangle=(E_F-J\Delta)|x\rangle+\frac J2(|x-1\rangle+|x+1\rangle)$（$L\ge3$）。

**命題 9.11（一マグノン）.** 環のセクター $N=1$ において、$|k\rangle=\sum_xe^{-ikx}|x\rangle$、
$k=2\pi m/L$ は $T|k\rangle=e^{ik}|k\rangle$ と $H|k\rangle=\big(E_F+J(\cos k-\Delta)\big)|k\rangle$ を満たす。

*証明.* $e^{ikL}=1$ を用いると $\sum_xe^{-ikx}(|x-1\rangle+|x+1\rangle)=\sum_y(e^{-ik(y+1)}+e^{-ik(y-1)})|y\rangle=2\cos k\,|k\rangle$ である。
同じ添字のつけかえで $T|k\rangle=e^{ik}|k\rangle$ が得られる。$\square$

二つのマグノンは散乱する。$x_1<x_2$ に対し、上向きスピンが $x_1,x_2$ にある状態を $|x_1,x_2\rangle$ と書く。隣り合う対は
二つのボンドを、離れた対は四つのボンドを壊す。

**命題 9.12（二マグノン；ベーテ方程式）.** $k_1,k_2\in\mathbb C$、$K=k_1+k_2$ とし、
$(A_{12},A_{21})\neq0$ が (i) $A_{12}(1+e^{iK}-2\Delta e^{ik_2})+A_{21}(1+e^{iK}-2\Delta e^{ik_1})=0$ と
(ii) $A_{21}=e^{ik_2L}A_{12}$、$A_{12}=e^{ik_1L}A_{21}$ を満たすとする。$\mathbb Z^2$ 上で
$\psi(x_1,x_2)=A_{12}e^{i(k_1x_1+k_2x_2)}+A_{21}e^{i(k_2x_1+k_1x_2)}$ とし、
$\Psi=\sum_{0\le x_1<x_2<L}\psi(x_1,x_2)|x_1,x_2\rangle$ とおく。$\Psi\neq0$ ならば
$H\Psi=\big(E_F+J(\cos k_1+\cos k_2-2\Delta)\big)\Psi$ である。分母が消えなければ、(i) と (ii) は次のことをいう：
$$
e^{ik_1L}=-\frac{1+e^{iK}-2\Delta e^{ik_1}}{1+e^{iK}-2\Delta e^{ik_2}},\qquad
e^{ik_2L}=-\frac{1+e^{iK}-2\Delta e^{ik_2}}{1+e^{iK}-2\Delta e^{ik_1}} .
$$

*証明.* $E=E_F+J(\cos k_1+\cos k_2-2\Delta)$ とおく。(1) すべての $(x_1,x_2)\in\mathbb Z^2$ に対し、各平面波について
$\psi(x_1\pm1,x_2)$ と $\psi(x_1,x_2\pm1)$ の和は $2(\cos k_1+\cos k_2)\psi(x_1,x_2)$ である。(2) 代入すると、
$\psi(x,x)+\psi(x+1,x+1)-2\Delta\psi(x,x+1)$ は (i) の左辺の $e^{iKx}$ 倍であるから、すべての $x$ に対して消える。
(3) (ii) により $\psi(x_2,x_1+L)=A_{12}e^{ik_2L}e^{i(k_2x_1+k_1x_2)}+A_{21}e^{ik_1L}e^{i(k_1x_1+k_2x_2)}=\psi(x_1,x_2)$ であり、
これを二回適用して $\psi(x_1+L,x_2+L)=\psi(x_1,x_2)$ を得る。配置 $\{y_1<y_2\}$ の
$D=\{x_1<x_2<x_1+L\}$ における代表元は $(y_1+nL,y_2+nL)$ と $(y_2+nL,y_1+(n+1)L)$ であるから、$\psi$ は配置の関数である。
$(x_1,x_2)\in D$ とする。$2\le x_2-x_1\le L-2$ ならば、四つのホッピングはすべて $D$ にとどまり、対角成分は $E_F-2J\Delta$ であるから、(1) により
$(H\Psi)(x_1,x_2)=(E_F-2J\Delta)\psi+\frac J2\cdot2(\cos k_1+\cos k_2)\psi=E\psi$ である。$x_2=x_1+1$ ならば、ホッピング
$(x_1+1,x_2)$、$(x_1,x_2-1)$ は阻まれ、対角成分は $E_F-J\Delta$ である。(1) と (2) により
$(H\Psi)(x_1,x_2)=(E_F-J\Delta)\psi+\frac J2\big(2(\cos k_1+\cos k_2)\psi-\psi(x_1,x_1)-\psi(x_1+1,x_1+1)\big)=E\psi$ である。
$x_2=x_1+L-1$ ならば、代表元 $(x_2,x_1+L)$ を用いる。これは前の場合に当たる（$L\ge3$）。$\square$

命題 9.12 は、解が固有ベクトルを与えることをいうだけであり、すべての固有ベクトルがこのようにして得られるとはいっていない。
$k_{1,2}=K/2\pm q$、$c=\cos\frac K2$ と書く。すると $1+e^{iK}-2\Delta e^{ik_{1,2}}=2e^{iK/2}(c-\Delta e^{\pm iq})$ であるから、
$A_{12}=c-\Delta e^{iq}$、$A_{21}=-(c-\Delta e^{-iq})$ は (i) を解く。(ii) の二つの式を掛け合わせると
$e^{iKL}=1$、$K=2\pi n/L$ となり、このとき (ii) は $F(q)=(-1)^ne^{iqL}(c-\Delta e^{-iq})+c-\Delta e^{iq}=0$ と同値である。
これは一変数の一つの方程式である。その根を、$c-\Delta e^{\pm iq}$ が互いに複素共役であるために
$e^{-iqL/2}F(q)$ が実数または純虚数となる $q\in(0,\pi)$ で探し（散乱状態）、また $F$ が実数となる直線
$q=i\kappa$、$\pi+i\kappa$（$\kappa>0$）の上で探す（束縛状態）。
$L=9$、$\Delta=0.6$ に対して各根について $\Psi$ を構成し、$36$ 個の状態からなるセクター $N=2$ と比較する。

```python
L, Delta = 9, 0.6
ring9 = SpinChain(L, boundary='periodic')
H2 = ring9.xxz(J=1.0, Delta=Delta, sector=2).matrix[0, 0]
x1, x2 = np.triu_indices(L, 1)                                       # all pairs x1 < x2
rows = ring9.index_of[torch.tensor((1 << x1) | (1 << x2))]
def roots(g, grid):                                                  # sign changes of g on the grid, refined
    v = np.array([g(t) for t in grid])
    return [brentq(g, grid[i], grid[i + 1], xtol=1e-15) for i in np.nonzero(v[:-1] * v[1:] < 0)[0]]
E_bethe, vecs, n_bound = [], [], 0
for n in range(L):
    K, c = 2 * np.pi * n / L, math.cos(np.pi * n / L)
    F = lambda q: (-1) ** n * np.exp(1j * q * L) * (c - Delta * np.exp(-1j * q)) + c - Delta * np.exp(1j * q)
    G = lambda q: (lambda z: z.real + z.imag)(np.exp(-0.5j * q * L) * F(q))    # real, or imaginary
    bound = [s + 1j * t for s in (0, np.pi)                          # q = i kappa and pi + i kappa
             for t in roots(lambda t: F(s + 1j * t).real, np.linspace(1e-3, 20, 4000))]
    n_bound += len(bound)
    for q in roots(G, (np.arange(1, 2000) + 0.5 ** 0.5) * np.pi / 2000) + bound:
        psi = np.exp(0.5j * K * (x1 + x2)) * ((c - Delta * np.exp(1j * q)) * np.exp(1j * q * (x1 - x2))
                                              - (c - Delta * np.exp(-1j * q)) * np.exp(-1j * q * (x1 - x2)))
        v = torch.zeros(36, dtype=c128)
        v[rows] = torch.tensor(psi / np.linalg.norm(psi))
        E = 9 * Delta / 4 + 2 * (c * np.cos(q) - Delta)                 # E_F + J (cos k1 + cos k2 - 2 Delta)
        assert abs(E.imag) < 1e-12 and close(H2 @ v, E.real * v, 1e-9)   # Proposition 9.12
        E_bethe.append(E.real)
        vecs.append(v)
assert torch.linalg.matrix_rank(torch.stack(vecs)) == 36             # a complete eigenbasis
assert close(torch.sort(torch.tensor(E_bethe)).values, torch.linalg.eigvalsh(H2))
print(len(E_bethe), n_bound)
ring8, Phi = SpinChain(8, 'periodic'), torch.zeros(28, dtype=c128)     # even L: Phi = sum_x (-1)^x |x, x+1>
for x in range(8):
    Phi[ring8.index_of[(1 << x) | (1 << (x + 1) % 8)]] = (-1) ** x
assert close(ring8.xxz(1.0, Delta, sector=2).matrix[0, 0] @ Phi, (8 * Delta / 4 - Delta) * Phi)   # E_F - J Delta
```

$L=9$ では $36$ 個の固有状態すべてがベーテ方程式から得られる：散乱状態が $34$ 個、束縛状態が $2$ 個である。
完全性はここで確かめられたのであって証明されたのではなく、偶数の $L$ では成り立たない。$\Delta\ne0$ ならば、$\Phi=\sum_x(-1)^x|x,x+1\rangle$ は
$T\Phi=-\Phi$ と $H\Phi=(E_F-J\Delta)\Phi$ を満たす（$|x-1,x+1\rangle$ と $|x,x+2\rangle$ へのホッピングは対ごとに打ち消し合う）。
ベーテ・ベクトルは $\psi(x_1-1,x_2-1)=e^{-iK}\psi(x_1,x_2)$ であるから $T\Psi=e^{-iK}\Psi$ を満たす。$e^{iK}=-1$ ならば
$\cos k_1+\cos k_2=2\cos\frac K2\cos\frac{k_1-k_2}2=0$ であり、そのエネルギーは $E_F-2J\Delta$ である。よって $\Phi$ はすべての
ベーテ・ベクトルと直交する（コードの最後の数行は $L=8$ に対して $H\Phi$ を確かめている）。$\Delta=0.6$、$4\le L\le13$ では、上の
探索は奇数の $L$ に対してはすべての状態を、偶数の $L$ に対しては $\Phi$ 以外のすべての状態を見つける。$N>2$ 個のマグノンに対しては、仮説は $N!$ 個の項と $N$ 個の連立したベーテ方程式をもつ（H. Bethe, *Z. Phys.* 71 (1931)
205; M. Karbach, G. Müller, *Computers in Physics* 11 (1997) 36）。

> **物理的意味.** マグノンはスピン波である：一つの反転したスピンが運動量 $k$ をもって環全体に
> 非局在化したものである。二つのマグノンは位相のずれを伴って散乱する（実数の $q$）か、束縛される（複素数の $q$）。ベーテの洞察は、
> $N$ 個のマグノンの散乱が二体の位相に因子分解されることであった。彼の方程式は、たとえば無限に長いハイゼンベルク鎖の
> サイトあたりの基底状態エネルギー $\frac14-\log2$ を与える（L. Hulthén, *Ark. Mat. Astron.
> Fys.* 26A (1938) no. 11。証明なしで述べる）。

## 9.5 ランチョス法とクリロフ法

セクター $N=L/2$ の次元は $L=16$ で $12870$、$L=20$ で $184756$ であるが、$H$ は疎であり（各行に $O(L)$ 個の非零
成分、`xxz_sparse`）、クリロフ法は積 $Hv$ だけを必要とする。$H$ をサイズ $n$ のエルミート行列とし、
$v\ne0$ とする。**クリロフ空間** は $\mathcal K_m=\operatorname{span}\{v,Hv,\dots,H^{m-1}v\}$ である。

**命題 9.13（ランチョス）.** $\dim\mathcal K_m=m$ と仮定する。$q_1=v/\|v\|$ とし、$j<m$ に対して $q_{j+1}$ を
$Hq_j$ の $\mathcal K_j$ に直交する成分を正規化したものとする。このとき $q_1,\dots,q_m$ は
$\mathcal K_m$ の正規直交基底である。$Q_m=(q_1\cdots q_m)$ とすると、$T_m=Q_m^*HQ_m$ は実対称三重対角行列で、対角成分は
$\alpha_j=\langle q_j,Hq_j\rangle$、非対角成分は $\beta_j=\langle q_{j+1},Hq_j\rangle>0$ である。また $HQ_m=Q_mT_m+r_me_m^T$ であり、ここで
$r_m\perp\mathcal K_m$ が消えるのは $\mathcal K_{m+1}=\mathcal K_m$ のとき、かつそのときに限る（そうでなければ $r_m=\beta_mq_{m+1}$、$q_{m+1}=r_m/\|r_m\|$）。

*証明.* $\operatorname{span}\{q_1,\dots,q_j\}=\mathcal K_j$ と仮定する（$j=1$ では正しい）。$Hq_i\in\mathcal K_{i+1}$ であるから
$H\mathcal K_j\subseteq\mathcal K_j+\mathbb CHq_j$ であり、したがって $\mathcal K_{j+1}=\mathcal K_j+\mathbb CH^jv\subseteq\mathcal K_j+\mathbb CHq_j\subseteq\mathcal K_{j+1}$、
すなわち $\mathcal K_{j+1}=\mathcal K_j+\mathbb CHq_j$ である。$j<m$ では次元が増えるので、$Hq_j$ の
$\mathcal K_j$ に直交する成分は零でなく、帰納法が進む。$j=m$ ではその成分 $r_m$ が消えるのは
$\mathcal K_{m+1}=\mathcal K_m$ のとき、かつそのときに限る。$i>j+1$ に対しては $q_i\perp\mathcal K_{j+1}\ni Hq_j$ であるから $(T_m)_{ij}=0$ であり、$T_m$ は
エルミートであるから $(T_m)_{ji}=0$ である。$\alpha_j$ は実数であり、$\beta_j=\langle q_{j+1},Hq_j\rangle$ はその成分のノルムで $>0$ である。よって
$j<m$ に対して $Hq_j=\beta_{j-1}q_{j-1}+\alpha_jq_j+\beta_jq_{j+1}$、また $Hq_m=\beta_{m-1}q_{m-1}+\alpha_mq_m+r_m$ である。$\square$

$T_m$ の固有値 $\theta^{(m)}_1\le\dots\le\theta^{(m)}_m$ を **リッツ値** という。

**定理 9.14（単調収束）.** $\lambda_1\le\dots\le\lambda_n$ を $H$ の固有値とし、
$\dim\mathcal K_m=m$ とする。このとき $\theta^{(m)}_j\ge\lambda_j$ であり、$\dim\mathcal K_{m+1}=m+1$ ならば $\theta^{(m)}_j\ge\theta^{(m+1)}_j$
（$j\le m$）である。$v$ が $\lambda_1$ の固有空間に零でない成分をもち、$\mathcal K_{m+1}=\mathcal K_m$ ならば
$\theta^{(m)}_1=\lambda_1$ である。

*証明.* クーラン–フィッシャー（標準的な線形代数）により、$\theta^{(m)}_j=\min_{U}\max_{x\in U,\|x\|=1}\langle x,Hx\rangle$
であり、最小は $j$ 次元の $U\subseteq\mathcal K_m$ にわたってとる。$Q_m$ は $\mathcal K_m$ の上への等長写像であり、
$\langle y,T_my\rangle=\langle Q_my,HQ_my\rangle$ だからである。同じ式を $U\subseteq\mathcal K_{m+1}$ にわたってとると $\theta^{(m+1)}_j$ が、
$U\subseteq\mathbb C^n$ にわたってとると $\lambda_j$ が得られる。より多くの部分空間にわたって最小化すれば値は小さくなる。
$\mathcal K_{m+1}=\mathcal K_m$ ならば $\mathcal K_m$ は $H$ 不変であり、任意の多項式 $p$ に対して $p(H)v$ を含み、
とくに $P_1v\ne0$ を含む（$P_1$ は $\lambda_1$ のスペクトル射影で、定理 1.5 と同様にラグランジュ補間により
$H$ の多項式である）。よって $x=P_1v/\|P_1v\|\in\mathcal K_m$ において $\theta^{(m)}_1\le\langle x,Hx\rangle=\lambda_1$ である。$\square$

$\theta^{(m)}_1\to\lambda_1$ の速さは、スペクトルの広がりに対するギャップ $\lambda_2-\lambda_1$ の比による
（チェビシェフ多項式によるカニエル–ペイジ–サードの評価。Y. Saad, *Numerical Methods for Large Eigenvalue
Problems*, SIAM 2011, Ch. 6。証明なしで述べる）。浮動小数点演算では $q_j$ は直交性を失う。`krylov.lanczos` は
完全な再直交化を行う。XX 鎖については厳密な基底状態エネルギーがわかっている（定理 9.10）。

```python
L, N = 12, 6
xx = SpinChain(L, complex_valued=False)
Hs = xx.xxz_sparse(J=1.0, Delta=0.0, sector=N)                      # sparse: diagonal + hopping
alpha, beta, Q = krylov.lanczos(Hs.matvec, torch.randn(1, xx.sector_dim(N), dtype=f64), 40)
T = torch.diag(alpha[0]) + torch.diag(beta[0], 1) + torch.diag(beta[0], -1)
assert close(Q[0] @ Q[0].T, torch.eye(40, dtype=f64)) and close(Q[0] @ Hs.to_dense()[0] @ Q[0].T, T)
ritz = torch.tensor([torch.linalg.eigvalsh(T[:m, :m])[0].item() for m in range(1, 41)])   # T_m = leading block
E_xx = sum(sorted(math.cos(math.pi * k / (L + 1)) for k in range(1, L + 1))[:N])
assert (ritz[1:] <= ritz[:-1] + 1e-12).all() and (ritz >= E_xx - 1e-12).all()             # Theorem 9.14
print([f"{(ritz[m - 1] - E_xx).item():.1e}" for m in (5, 10, 20, 40)])
```

最小のリッツ値の誤差は単調に減少し、$5$、$10$、$20$ ステップの後でそれぞれ $0.42$、$3.1\cdot10^{-3}$、$3.1\cdot10^{-6}$ となり、
$40$ ステップの後には丸め誤差の水準に達する。次元は $924$ である。`krylov.ground_state` は
残差 $\|H\psi-E\psi\|$ が $10^{-10}$ を下回るまでランチョス法を再始動する。これを、`scripts/chains.py` に九桁まで
表にされたハイゼンベルク環のエネルギー（付随する論文の参照値）と照合する。$L=12$ までは密な方法で、
$L=16$ ではランチョス法で確かめる。

```python
RING = {4: -2.0, 6: -2.802775638, 8: -3.651093409, 10: -4.515446354, 12: -5.387390917,
        14: -6.263549534, 16: -7.142296361}
for n in (4, 6, 8, 10, 12):
    e0 = torch.linalg.eigvalsh(SpinChain(n, 'periodic').heisenberg(sector=n // 2).matrix[0, 0])[0].item()
    assert abs(e0 - RING[n]) < 5e-10
    print(n, round(e0 / n, 6))
t0 = time.time()
E16, psi16, res = krylov.ground_state(SpinChain(16, 'periodic', complex_valued=False).xxz_sparse(sector=8),
                                      generator=torch.Generator().manual_seed(0))
print(f"L = 16: E0 = {E16.item():.10f}, residual {res.item():.1e}, {time.time() - t0:.2f} s; "
      f"infinite chain: {0.25 - math.log(2):.6f}")
assert abs(E16.item() - RING[16]) < 5e-10 and res.item() < 1e-9
```

$L=16$（次元 $12870$）では、ランチョス法は一つの CPU スレッドで一秒を十分に下回る時間で終わる。サイトあたりのエネルギー
$-0.5,-0.467129,-0.456387,-0.451545,-0.448949$（$L=4,\dots,12$）は、無限に長い鎖の値 $\frac14-\log2=-0.443147$ に
下から近づく。

**時間発展.** $e^{-i\tau H}v$ に対しては $\mathcal K_m$ に射影する：$w_\tau=\|v\|\,Q_me^{-i\tau T_m}e_1$。

**命題 9.15.** $\dim\mathcal K_m=m$ とする。任意の $\tau$ に対して $\|w_\tau\|=\|v\|$ かつ $\langle w_\tau,Hw_\tau\rangle=\langle v,Hv\rangle$ である。
$\mathcal K_m$ が $H$ 不変ならば $w_\tau=e^{-i\tau H}v$ である。

*証明.* $Q_m$ は等長写像であり、$e^{-i\tau T_m}$ はユニタリである（$T_m$ は実対称）から、ノルムについての主張が従う。
$Q_m^*HQ_m=T_m$ は $e^{-i\tau T_m}$ と可換であるから、
$\langle w_\tau,Hw_\tau\rangle=\|v\|^2\langle e^{-i\tau T_m}e_1,T_me^{-i\tau T_m}e_1\rangle=\|v\|^2\langle e_1,T_me_1\rangle=\langle v,Hv\rangle$ である。
$\mathcal K_m$ が不変ならば $HQ_m=Q_mT_m$（命題 9.13、$r_m=0$）であるから、$H^jQ_m=Q_mT_m^j$ かつ $e^{-i\tau H}Q_m=Q_me^{-i\tau T_m}$ である。
これを $e_1$ に作用させ、$Q_me_1=v/\|v\|$ を用いればよい。$\square$

したがって、クリロフ法の各ステップは、正確であろうとなかろうとノルムとエネルギーを保存する。精度を得るには $m$ が
$\tau\|H\|$ に比べて大きい必要があり、`krylov.evolve` は時間区間を分割し、事後的に誤差を評価する
（M. Hochbruck, C. Lubich, *SIAM J. Numer. Anal.* 34 (1997) 1911。証明なしで述べる）。ハイゼンベルク環の上で
ネール状態 $|{\uparrow\downarrow\uparrow\downarrow}\cdots\rangle$ からクエンチを行い、交替磁化
$m_s=\frac2L\sum_i(-1)^i\langle s^z_i\rangle$ を追跡し、$L=12$ については密な伝播子と比較する。

```python
times, m_s = [0.0, 0.5, 1.0, 2.0, 4.0], {}
stag = lambda c, N, out: [((p[0].abs() ** 2) @ (c.bits(N).to(f64) - 0.5) @ torch.tensor(
    [(-1.0) ** i for i in range(c.L)], dtype=f64) * 2 / c.L).item() for p in out]
for n in (12, 16):
    c = SpinChain(n, 'periodic')
    Hn = c.xxz_sparse(J=1.0, Delta=1.0, sector=n // 2)
    psi0 = c.vector_in_sector('10' * (n // 2))[None]
    out = krylov.evolve(Hn, psi0, times)                               # (times, batch, dim)
    energy = torch.stack([(p.conj() * Hn.matvec(p)).sum(-1).real for p in out])
    assert close(out.norm(dim=-1), 1.0) and close(energy, energy[0])   # Proposition 9.15
    m_s[n] = stag(c, n // 2, out)
    print(n, [round(x, 4) for x in m_s[n]])
    if n == 12:
        w, U = torch.linalg.eigh(c.xxz(1.0, 1.0, sector=6).matrix[0, 0])
        for t, p in zip(times, out):
            assert close(p[0], U @ (torch.exp(-1j * t * w) * (U.mH @ psi0[0])))
print(f"{max(abs(a - b) for a, b in zip(m_s[12][:4], m_s[16][:4])):.1e}")   # L = 12 vs 16 up to t = 2
```

$m_s$ は $1$ から $0.770$、$0.279$ を経て $t=2$ で $-0.183$ まで減衰し、二つのサイズで $1.7\cdot10^{-5}$ の範囲で一致する。
$t=4$ では両者は異なり（$0.026$ 対 $0.035$）、小さいほうの環の有限サイズ効果が見えてくる。

## 9.6 乱れ：準位統計とエンタングルメント

**ランダム磁場ハイゼンベルク鎖** は、$[-W,W]$ 上の一様分布に従う独立な磁場 $h_i$ をもつ。これは $S^z$ を保存し
（命題 9.1）、一般にはそれ以外は何も保存しない：磁場は並進、反転、スピン反転、
$SU(2)$ を破る。したがって準位統計はセクター $N$ ごとに計算する（定理 5.15 と 5.4節：セクターを混ぜると
独立な準位を偽造してしまう）。弱い乱れでは準位は GOE と同様に反発し、$\langle r\rangle\approx0.5307$ である（5.3節）。
強い乱れでは固有状態が局在し（**多体局在**）、準位は独立になり、
$\langle r\rangle=2\log2-1\approx0.3863$ となる（命題 5.14）。セクター $N=L/2$ の中央の半分で $\langle r\rangle$ を計算し、
独立な実現にわたる標準誤差を一つ添える（各実現が平均比を一つ与える）。

```python
gen = torch.Generator().manual_seed(1)
r_mean = {}
for n, R in ((8, 200), (10, 100)):
    c = SpinChain(n, 'periodic', complex_valued=False)
    k = c.sector_dim(n // 2)
    for W in (0.5, 1, 2, 4, 8, 16):
        E = torch.linalg.eigvalsh(c.random_field_heisenberg(W, R, sector=n // 2, generator=gen).matrix[:, 0])
        r = SpinChain.level_spacing_ratio(E[:, k // 4: 3 * k // 4])     # one value per realisation
        r_mean[n, W] = (r.mean().item(), (r.std() / math.sqrt(R)).item())
    print(n, '  '.join(f"W={W}: {m:.3f}+-{s:.3f}" for (n_, W), (m, s) in r_mean.items() if n_ == n))
assert abs(r_mean[10, 0.5][0] - 0.5307) < 0.015 and abs(r_mean[10, 16][0] - 0.3863) < 0.015
assert r_mean[10, 1][0] > r_mean[10, 2][0] > r_mean[10, 4][0] > r_mean[10, 8][0]
```

$L=10$ では平均比は $0.530$（$W=0.5$）から $0.484$（$W=2$）、$0.410$（$W=4$）を経て $0.384$
（$W=16$）まで変化し、標準誤差は $0.003$–$0.004$ である。$L=8$ では $W=0.5$ での値（$0.507$）が
$W=1$ での値（$0.525$）より小さい：$W=0$ では鎖は可積分で上記のすべての対称性をもち、GOE の値を得るには
乱れが準位間隔に比べて大きくなければならない。クロスオーバーは幅広く、$L$ とともにずれていく。有限サイズ
スケーリングはそれを $W\approx3.7$ 付近に置き（D. Luitz, N. Laflorencie, F. Alet, *Phys. Rev. B* 91 (2015) 081103）、
局在相が $L\to\infty$ で生き残るかどうかは議論が続いている（J. Šuntajs, J. Bonča, T. Prosen, L. Vidmar,
*Phys. Rev. E* 102 (2020) 062144）。

**自由フェルミオンのエンタングルメント.** 開いた XX 鎖の固有状態について、最初の $n$ 個のサイトの
エンタングルメント・エントロピーは $n\times n$ 行列によって決まる。

**定理 9.16（ペシェルの公式）.** $\Phi=b_{k_1}^*\cdots b_{k_N}^*|0\rangle$ を定理 9.10 のとおりとし、
$C_{ij}=\langle\Phi,c_i^*c_j\Phi\rangle$、$C_A=(C_{ij})_{i,j<n}$ とおく。このとき $0\le C_A\le1$ であり、サイト
$A=\{0,\dots,n-1\}$ のエンタングルメント・エントロピーは、$\nu$ が $C_A$ の固有値を重複度込みで動くとして、次で与えられる。
$$
S_A=\sum_\nu \eta(\nu),\qquad \eta(\nu)=-\nu\log\nu-(1-\nu)\log(1-\nu)\quad(0\log0=0).
$$

*証明.* 実ベクトル $x$ に対し $b(x)=\sum_jx_jc_j$ とおく。すると $b(x)^*b(y)^*=-b(y)^*b(x)^*$ かつ
$\{b(x),b(y)^*\}=\langle x,y\rangle$ である。$U\subseteq\mathbb R^L$ を $u_{k_a}$ の張る空間とし、$P_A$、$P_B=1-P_A$ を
$\{0,\dots,n-1\}$ とその補集合への座標射影とする。$U$ 上の形式 $\langle P_Ax,P_Ay\rangle$ を対角化すると、
$\langle P_A\phi_a,P_A\phi_b\rangle=\nu_a\delta_{ab}$ を満たす $U$ の正規直交基底 $\phi_1,\dots,\phi_N$ が得られ、したがって
$\langle P_B\phi_a,P_B\phi_b\rangle=(1-\nu_a)\delta_{ab}$ である。$(x_1,\dots,x_N)\mapsto b(x_1)^*\cdots b(x_N)^*|0\rangle$ は
多重線形かつ交代的であるから、$f_a=b(P_A\phi_a/\sqrt{\nu_a})$、$g_a=b(P_B\phi_a/\sqrt{1-\nu_a})$ として
$\Phi=\pm\prod_a\big(\sqrt{\nu_a}f_a^*+\sqrt{1-\nu_a}\,g_a^*\big)|0\rangle$ である（係数が零の項は省く）。
展開すると、部分集合 $I$ にわたる和として $\Phi=\sum_I\pm\prod_{a\in I}\sqrt{\nu_a}\prod_{a\notin I}\sqrt{1-\nu_a}\,F_I^*G_{I^c}^*|0\rangle$ となる。
ここで $F_I^*$、$G_{I^c}^*$ は順序づけられた積である。$j\ge n$ に対しては $c_j^*$ のストリングは $A$ 全体を含む：
$c_j^*=(-1)^{N_A}\otimes\tilde c_j^*$ であり、$\tilde c_j$ はサイト $\ge n$ のジョルダン–ウィグナー演算子である。一方
$i<n$ に対しては $c_i^*=\hat c_i^*\otimes1$ である。$(-1)^{N_A}$ は $|0_A\rangle$ の上で $1$ であるから、
$F_I^*G_{I^c}^*|0\rangle=\xi_I\otimes\chi_{I^c}$ であり、$\xi_I=\hat F_I^*|0_A\rangle$、$\chi_{I^c}=\tilde G_{I^c}^*|0_B\rangle$ は
定理 9.10 の CAR の計算により正規直交系である。これはシュミット分解（補題 7.7）であり、その
重み $p_I=\prod_{a\in I}\nu_a\prod_{a\notin I}(1-\nu_a)$ は独立なベルヌーイ変数の分布であるから、
$S_A=-\sum_Ip_I\log p_I=\sum_a\eta(\nu_a)$ である。最後に $C_{ij}=\sum_a\phi_a(i)\phi_a(j)$
（$c_i^*c_j$ を $b(\phi_a)$ で展開せよ）であるから、$C_A=\sum_a(P_A\phi_a)(P_A\phi_a)^T$ は
直交するベクトル $P_A\phi_a$ の上で固有値 $\nu_a$ をもち、それ以外では $0$ であり、$\eta(0)=0$ である。$\square$

`entanglement_entropy` は縮約状態をセクターごとに計算し（それはブロックの磁化と可換である、
命題 7.8）、そのエントロピーを求める。セクター $N$ の基底状態は $N$ 個の最低のモードからなる $\Phi$ であり、
$\varepsilon_k$ が相異なるので単純である。

```python
L, N = 12, 6
xx = SpinChain(L)
w, V = torch.linalg.eigh(xx.xxz(J=1.0, Delta=0.0, sector=N).matrix[0, 0])
rho = DensityMatrix(xx.sector_algebra(N), matrix=torch.outer(V[:, 0], V[:, 0].conj())[None, None])
eps_k, u = np.linalg.eigh(np.eye(L, k=1) / 2 + np.eye(L, k=-1) / 2)  # h = A_L / 2, ascending
C = u[:, :N] @ u[:, :N].T                                            # correlation matrix of the N lowest modes
bern = lambda nu: -(nu * np.log(nu) + (1 - nu) * np.log(1 - nu))
S_peschel = [bern(np.clip(np.linalg.eigvalsh(C[:n, :n]), 1e-15, 1 - 1e-15)).sum() for n in range(1, L)]
S_lib = [xx.entanglement_entropy(rho, n, sector=N).item() for n in range(1, L)]
assert np.allclose(S_lib, S_peschel, atol=1e-9)                       # Theorem 9.16
print([round(s, 4) for s in S_lib[:6]])
```

エントロピー $0.6931,0.5199,0.7456,0.6283,0.7719,0.6560$（$n=1,\dots,6$）は $n$ の偶奇とともに振動し、
$n/L$ を固定すると $\frac16\log L$ のように増大する。これは開いた鎖に対する共形場理論の法則である（I. Peschel,
*J. Phys. A* 36 (2003) L205; P. Calabrese, J. Cardy, *J. Stat. Mech.* (2004) P06002。証明なしで述べる）。

## 9.7 比較の対象

QuSpin（Python）は磁化、運動量、パリティ、スピン反転のセクターでスピン鎖のハミルトニアンを構成し、
SciPy の疎な固有値ソルバーと行列指数関数を用いる。専用の厳密対角化コード
（XDiag、SPINPACK、HΦ）は大規模な疎問題向けに最適化されており、格子の対称性と大型計算機を用いて
そのようなコードは約 40–50 スピンに達している。ITensor、TeNPy、TensorKit.jl（$SU(2)$ の融合木をもつ）は
行列積状態によって長い鎖を扱う。このライブラリが付け加えるものは控えめで具体的である：鎖を
`TypeIAlgebra` $\bigoplus_NM_{\binom LN}$ として扱い、第3–8章の状態、チャネル、エントロピー、ダイナミクスを
セクターごとに適用できること。指数をもつ包含の塔。最高ウェイトベクトルから得られる $SU(2)$ ブロックと
ウィグナー–エッカートによる相関。一度の呼び出しで乱れの実現のバッチを扱えること。そしてセクターから決して
出ない縮約状態。反転やスピン反転のセクターはなく、運動量セクターは一様な XXZ 鎖に限られ、密な
方法は $L\approx16$ まで、クリロフ法は $L\approx24$ まで（大型 GPU では約 $26$ まで）であり、テンソルネットワークはない。

## 演習問題

1. （運動量セクターを手で求める。）$L=4$、$N=2$ に対し、六つの配置の $T$ による軌道と
   次元 $d_{2,m}$（$m=0,\dots,3$）を求めよ。
2. （和則。）$\psi$ を $L$ が偶数のハイゼンベルク環の基底状態とする。
   $E_0=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle$ かつ $\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle=0$ を示せ。
3. （対称性が多すぎる。）セクター $N=L/2$ において、すべての係数が $1$ に等しいベクトルは
   $(S^-)^{L/2}|{\uparrow\cdots\uparrow}\rangle$ の定数倍であり、ハイゼンベルク環のエネルギー
   $JL/4$ の固有ベクトルであること、また任意のクリロフ空間 $\mathcal K_m(H,v)$ は $v$ を含む任意の $H$ 不変部分空間に含まれることを示せ。
   このベクトルから始めたとき、`krylov.ground_state` は何を返すか。

<details>
<summary>解答</summary>

1. $0011\to1001\to1100\to0110\to0011$ は $R=4$ の軌道であり、$0101\to1010\to0101$ は $R=2$ の軌道である。
   命題 9.5 により、前者はすべての $m$ に寄与し、後者は $2m\equiv0\pmod4$ を満たす $m$、すなわち $m=0,2$ に寄与する：
   $d_{2,m}=2,1,2,1$。

```python
assert SpinChain(4, 'periodic').momentum_dims(2) == [2, 1, 2, 1]
```

2. 基底状態は単純な一重項であり（定理 9.9）、$T$ は $H$ と可換であるから $T\psi=e^{ik}\psi$ である。
   $\mathbf s_i\cdot\mathbf s_{i+1}=T^i(\mathbf s_0\cdot\mathbf s_1)T^{*i}$ であるから、$\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle$ は $L$ 本のボンドすべてで等しく、
   $E_0=J\sum_i\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle$ である。$\mathbf S^2\psi=0$ であるから、
   再び並進不変性により $0=\langle\mathbf S^2\rangle=\sum_{i,j}\langle\mathbf s_i\cdot\mathbf s_j\rangle=L\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle$ である。

```python
assert abs(12 * corr_we[1].item() - gs.energy.item()) < 1e-10 and abs(corr_we.sum().item()) < 1e-10
```

3. $S^-$ は上向きスピンの一つを下げる。完全に偏極した状態にこれを $L/2$ 回作用させると、上向きスピンを
   $L/2$ 個もつ各配置に $(L/2)!$ 通りの順序で到達し、どれも係数は $1$ である。$|{\uparrow\cdots\uparrow}\rangle$ では
   $L$ 本のボンドすべてが平行でエネルギーは $JL/4$ であり、$S^-$ は $H$ と可換であるから、このベクトルは同じ
   エネルギーの固有ベクトルである。$W$ が不変で $v\in W$ ならば、すべての $j$ に対して $H^jv\in W$ である。ここではすべての $m$ に対して $\mathcal K_m=\mathbb Cv$ であり、
   `ground_state` は $E=JL/4$、すなわちスペクトルの最上端を返す（各 $\mathbf s_i\cdot\mathbf s_j\le\frac14$）。ランダムな初期
   ベクトルは確率一でこれを避ける。

```python
c10 = SpinChain(10, 'periodic', complex_valued=False)
E_sym, _, _ = krylov.ground_state(c10.xxz_sparse(1.0, 1.0, sector=5), v0=torch.ones(1, 252, dtype=f64))
assert abs(E_sym.item() - 2.5) < 1e-10 and E_sym.item() > RING[10] + 7     # far from E_0
```

</details>

**次へ：**[第10章](10_quantum_computing.md)ではセクターを用いて量子計算を扱う：デコヒーレンスフリー
部分空間とノイズレス部分系、作用素環による誤り訂正、そして保存されるハミング重みを尊重するミキサー——
本章の環 $\bigoplus_wM_{\binom nw}$ が新しい役割で現れる。
