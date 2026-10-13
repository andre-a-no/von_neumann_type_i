[English](../06_quantum_channels.md) · [Español](../es/06_quantum_channels.md) · [Français](../fr/06_quantum_channels.md) · [Deutsch](../de/06_quantum_channels.md) · [中文](../zh/06_quantum_channels.md) · **日本語** · [Русский](../ru/06_quantum_channels.md)

# 6. 量子チャネル

量子チャネルとは、状態が変化しうる最も一般的な仕方である：ユニタリ発展、結果が忘れられる測定、
環境との相互作用、粒子の損失などがその例である。数学的には、それは完全正値かつトレース保存の線形写像である。本章では、
有限次元フォン・ノイマン環の間のそのような写像について標準的な記述——チョイ行列、クラウス作用素、
スタインスプリング拡張——を証明し、続いて直和に特有の事柄を調べる：セクター間の写像、第3章のトレースに関する双対、
そしてチャネルの古典的部分である。第2章と第3章ではすでにライブラリのチャネルを用いた。ここではその理論を述べる。

**本章で学ぶこと**：正値性だけでは不十分である理由、チョイ行列・クラウス作用素・超演算子・等長作用素の間を
行き来する方法、写像がセクターを変えるとき「トレース保存」がトレースに依存する理由、デコヒーレンスを生き延びる
オブザーバブルは何か、そしてランダムなチャネルが純粋状態に平均として何をするか。

本章を通じて、$\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ および
$\mathcal N=\bigoplus_{d=1}^{C'}M_{k'_d}(\mathbb C)$ とし、それぞれの極小中心射影を $P_c$ と
$Q_d$ とする。$A=\bigoplus_cA_c$、$M_{k'\times k}$ は複素 $k'\times k$ 行列の空間、
$E_{ij}=|i\rangle\langle j|$ は行列単位、$\operatorname{Tr}$ は blunt トレース
$\sum_c\operatorname{Tr}A_c$、$\operatorname{tr}_w=\sum_cw_c\operatorname{Tr}A_c$ は忠実な
トレース（定理 3.1）である。ブロック行列は $\mathbb C^n\otimes\mathbb C^k$ に作用し、ブロックの添字を先に書く。

## 6.1 完全正値写像とチョイの定理

線形写像 $\Phi:\mathcal M\to\mathcal N$ は、$A\ge0$ ならば $\Phi(A)\ge0$ となるとき **正値** であるという。
成分が $\mathcal M$ に属する $n\times n$ 行列 $[A_{ij}]$ は再び行列環の直和をなす：
$[A_{ij}]\mapsto([(A_{ij})_c]_{ij})_c$ により $M_n(\mathcal M)\cong\bigoplus_cM_n(M_{k_c})$ であり、
したがって $[A_{ij}]\ge0$ とは、すべてのブロック行列 $[(A_{ij})_c]\in M_{nk_c}$ が正であることを意味する。
$\Phi_n:[A_{ij}]\mapsto[\Phi(A_{ij})]$ が正値のとき $\Phi$ は **$n$-正値** であるといい、すべての $n$ に対し
$n$-正値であるとき **完全正値**（CP）であるという。$\operatorname{Tr}\Phi(A)=\operatorname{Tr}A$ のとき
**トレース保存**（TP）、$\Phi(1)=1$ のとき **単位的** であるといい、CP かつ TP であるとき **チャネル** という。
任意の長方行列 $T$ に対し $X\ge0$ ならば $TXT^*\ge0$ であることを繰り返し用いる。実際
$\langle y,TXT^*y\rangle=\langle T^*y,XT^*y\rangle$ である。

**命題 6.1.** $k\ge2$ のとき、$M_k$ 上の転置 $\theta(A)=A^T$ は正値であるが
$2$-正値ではない。

*証明.* $A\ge0$ に対し $\langle x,A^Tx\rangle=\sum_{ij}\bar x_iA_{ji}x_j=\langle\bar x,A\bar x\rangle\ge0$ である。
$\Omega=e_1\otimes e_1+e_2\otimes e_2\in\mathbb C^2\otimes\mathbb C^k$ とすると
$|\Omega\rangle\langle\Omega|=[E_{ij}]_{i,j=1}^2\ge0$ であり、$\theta_2$ はこれを
$F=[E_{ji}]_{i,j=1}^2$ に写す。$F$ は $e_1\otimes e_2-e_2\otimes e_1$ をその $-1$ 倍に写す。$\square$

**定理 6.2（チョイ）.** 線形写像 $\Phi:M_k\to M_{k'}$ に対し、次は同値である：
(a) $\Phi$ は CP である。(b) $\Phi$ は $k$-正値である。(c) **チョイ行列**
$$
J(\Phi)=\sum_{i,j=1}^kE_{ij}\otimes\Phi(E_{ij})=[\Phi(E_{ij})]_{i,j=1}^k\in M_k(M_{k'})
$$
は正である。(d) $\Phi(A)=\sum_lK_lAK_l^*$ を満たす $K_1,\dots,K_r\in M_{k'\times k}$（**クラウス作用素**）が存在し、
$r=\operatorname{rank}J(\Phi)\le kk'$ ととれる。

*証明.* (a)$\Rightarrow$(b) は自明である。(b)$\Rightarrow$(c)：
$\Omega=\sum_ie_i\otimes e_i\in\mathbb C^k\otimes\mathbb C^k$ に対し
$|\Omega\rangle\langle\Omega|=[E_{ij}]\ge0$ であり、$J(\Phi)=\Phi_k([E_{ij}])$ である。(c)$\Rightarrow$(d)：
定理 1.3 により、$r=\operatorname{rank}J(\Phi)$ として $J(\Phi)=\sum_{l=1}^rv_lv_l^*$ と書ける（$v_l$ は
非零固有値の固有ベクトルを固有値の平方根倍したもの）。任意の $v\in\mathbb C^k\otimes\mathbb C^{k'}$
は、ただ一つの $K\in M_{k'\times k}$ により $v=\sum_ie_i\otimes Ke_i$ と書ける（$Ke_i$ は $v$ の第 $i$ ブロック）。
また $vv^*$ の $(i,j)$ ブロックは $Ke_ie_j^*K^*=KE_{ij}K^*$ である。$v_l$ をこのように
$K_l$ を用いて書き、ブロックを比較すると $\Phi(E_{ij})=\sum_lK_lE_{ij}K_l^*$ を得るので、線形性により (d) が従う。
(d)$\Rightarrow$(a)：$\Phi_n([A_{ij}])=\sum_l(1_n\otimes K_l)[A_{ij}](1_n\otimes K_l)^*\ge0$。$\square$

**系 6.3.** $\Phi(A)=\sum_{l=1}^sL_lAL_l^*$ ならば $s\ge\operatorname{rank}J(\Phi)$ である。したがって
クラウス作用素の最小個数、すなわち **クラウス階数** は $\operatorname{rank}J(\Phi)$ である。

*証明.* 上の計算を逆にたどると、
$w_l=\sum_ie_i\otimes L_le_i$ として $J(\Phi)=\sum_{l=1}^sw_lw_l^*$ であり、これは階数 $\le1$ の行列 $s$ 個の和である。$\square$

直和の間の写像は、ブロック間の写像を成分とする行列である。$\iota_c:M_{k_c}\to\mathcal M$
を行列をブロック $c$ に置く写像、$\pi_d:\mathcal N\to M_{k'_d}$ をブロック $d$ を取り出す写像とし、
$\Phi_{dc}=\pi_d\circ\Phi\circ\iota_c:M_{k_c}\to M_{k'_d}$ とおく。$A=\sum_c\iota_c(A_c)$ であるから
$\Phi(A)_d=\sum_c\Phi_{dc}(A_c)$ である。

**命題 6.4.** $\Phi:\mathcal M\to\mathcal N$ が $n$-正値であるための必要十分条件は、すべての $\Phi_{dc}$ が $n$-正値であることである。したがって
$\Phi$ が CP であるための必要十分条件は、$K^{dc}_l\in M_{k'_d\times k_c}$（$l=1,\dots,r$）が存在して
$$
\Phi(A)_d=\sum_c\sum_{l=1}^rK^{dc}_l\,A_c\,K^{dc*}_l ,
$$
となることであり、$r=\max_{d,c}\operatorname{rank}J(\Phi_{dc})$ で十分である。$\Phi$ が TP であるための必要十分条件は、
すべての $c$ に対し $\sum_{d,l}K^{dc*}_lK^{dc}_l=1_{k_c}$ となることである。

*証明.* すべての $\Phi_{dc}$ が $n$-正値で $[A_{ij}]\ge0$ ならば、各 $X^c=[(A_{ij})_c]\ge0$ であり、
$\Phi_n([A_{ij}])$ のブロック $d$ は $\sum_c(\Phi_{dc})_n(X^c)\ge0$ である。逆に、
$M_n(M_{k_c})$ の $[X_{ij}]\ge0$ に対し、$M_n(\mathcal M)$ において $[\iota_c(X_{ij})]\ge0$ である（そのブロック $c$ は $[X_{ij}]$ で、
他のブロックは零）。その像のブロック $d$ は $[\Phi_{dc}(X_{ij})]$ である。クラウス形は各組 $(d,c)$ に定理 6.2
を適用し、零作用素を補って共通の $r$ にそろえれば得られる。最後に、トレースの巡回性により
$\operatorname{Tr}\Phi(A)=\sum_c\operatorname{Tr}\big(A_c\sum_{d,l}K^{dc*}_lK^{dc}_l\big)$ であり、これがすべての $A$ に対し
$\sum_c\operatorname{Tr}A_c$ に等しいのは各和が $1_{k_c}$ であるときかつそのときに限る（行列単位をとればよい：
$\operatorname{Tr}(E_{ij}X)=X_{ji}$）。$\square$

ライブラリには二つのクラスがある。`Channel` は **セクター保存** 写像（$C'=C$、$k'=k$、
$d\ne c$ に対し $\Phi_{dc}=0$）を格納する。そのクラウス作用素 $K_l=\bigoplus_cK^{cc}_l$ は $\mathcal M$ に属し、
形状 `(batch, r, C, k_max, k_max)` のテンソルとして格納される。`InterSectorChannel` は
命題 6.4 の一般の写像を形状 `(batch, r, C', C, k'_max, k_max)` で格納する。`Channel.choi()` は各セクターについて
$J(\Phi_{cc})$ を返し、`Channel.superoperator()` は行優先のベクトル
$\operatorname{vec}(X)_{(a,b)}=X_{ab}$ に作用する $\Phi_{cc}$ の行列 $S_c$ を返す。どちらも同じ数
$\Phi(E_{ij})_{ab}=J_{(i,a),(j,b)}=S_{(a,b),(i,j)}$ を異なる順序で含んでいる。`Channel.from_superoperator`
は定理 6.2 の証明のようにチョイ行列を対角化し、CP でない写像を拒否する。

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
```

> **物理的意味.** 系は決して孤立していない：$n$ 準位の、手を触れられていない補助系と
> 相関しているかもしれない。物理的過程はこの組に $\Phi_n=\mathrm{id}_n\otimes\Phi$ として作用し、
> 状態を状態に写さなければならないので、$\Phi$ は完全正値でなければならない。転置は正値であるが、
> エンタングルした状態 $|\Omega\rangle\langle\Omega|/2$ を固有値 $-1/2$ をもつ作用素に写す。
> チョイ行列は最大エンタングル状態の像である：ただ一つの状態が完全正値性を決定する。

## 6.2 スタインスプリング拡張とクラウス作用素の自由度

$M_k$ から $M_{k'}$ への $\Phi(A)=\sum_{l=1}^rK_lAK_l^*$ に対し $\Phi^*(B)=\sum_lK_l^*BK_l$ とおく。
巡回性により $\operatorname{Tr}(\Phi(A)B)=\operatorname{Tr}(A\Phi^*(B))$ であるから、$\Phi^*$ はトレースに関する
$\Phi$ の双対であり（6.3節）、$K_l$ の選び方に依存しない。
$\mathbb C^r\otimes\mathbb C^{k'}$（**環境** $\mathbb C^r$ が先）上で
$\operatorname{Tr}_EX=\sum_l(e_l^*\otimes1)X(e_l\otimes1)$ とおく。

**定理 6.5（スタインスプリング、有限次元）.**
$V=\sum_le_l\otimes K_l:\mathbb C^k\to\mathbb C^r\otimes\mathbb C^{k'}$ とおく。このとき
$$
\Phi(A)=\operatorname{Tr}_E(VAV^*),\qquad\Phi^*(B)=V^*(1_r\otimes B)V,
$$
であり、$\Phi$ が TP であるための必要十分条件は $V$ が等長作用素であることである。逆に、任意の $V$ は第一の式により CP 写像を定める。

*証明.* $(e_l^*\otimes1)V=K_l$ であるから $\operatorname{Tr}_E(VAV^*)=\sum_lK_lAK_l^*$ である。逆に
$\operatorname{Tr}_E(VAV^*)$ は $K_l=(e_l^*\otimes1)V$ をもつクラウス形である。次に
$V^*(1\otimes B)V=\sum_{l,l'}\langle e_l,e_{l'}\rangle K_l^*BK_{l'}=\Phi^*(B)$ であり、
$V^*V=\sum_lK_l^*K_l$ は $\Phi$ が TP であるときかつそのときに限り $1$ である（命題 6.4）。$\square$

一般には（証明なしに述べるが）、単位的 C\*-環 $\mathcal A$ からヒルベルト空間 $H$ 上の有界作用素への
CP 写像 $\Psi$ は、$\mathcal A$ のヒルベルト空間 $K$ 上の表現 $\pi$ と有界な $V:H\to K$ により
$\Psi(a)=V^*\pi(a)V$ と書ける（W. F. Stinespring, Proc. AMS 6 (1955)
211; V. Paulsen, *Completely Bounded Maps and Operator Algebras*, Cambridge 2002, Thm. 4.1）。上では
$\pi(B)=1_r\otimes B$ である。直和の間の写像に対しては、
$V_c=\sum_{d,l}e_l\otimes\iota'_dK^{dc}_l:\mathbb C^{k_c}\to\mathbb C^r\otimes\bigoplus_d\mathbb C^{k'_d}$
（$\iota'_d$ は第 $d$ 直和成分の埋め込み）とおくと、$\Phi(A)_d$ は
$\sum_c\operatorname{Tr}_E(V_cA_cV_c^*)$ の $(d,d)$ ブロックとなり、$\Phi$ が TP であるのはすべての $V_c$ が
等長作用素であるときかつそのときに限る。`random_channel` と `random_inter_sector_channel` はこのようにしてチャネルを生成する（6.5節）。

**補題 6.6.** $X,Y\in M_{p\times n}$ かつ $XX^*=YY^*$ ならば、あるユニタリ $W\in M_n$ により $Y=XW$ である。

*証明.* すべての $z$ に対し $\|X^*z\|^2=\langle z,XX^*z\rangle=\|Y^*z\|^2$ であるから、$X^*z\mapsto Y^*z$ は
$X^*$ の値域から $Y^*$ の値域の上へのwell-definedな線形等長写像である。定理 1.8 の証明と同様に、
これは $UX^*=Y^*$ を満たすユニタリ $U$ に拡張される。$W=U^*$ とすればよい。$\square$

**定理 6.7（ユニタリ自由度）.** $K_1,\dots,K_n$ と $L_1,\dots,L_n$ を $M_{k'\times k}$ の元とする
（短い方の族は零で補う）。このとき、すべての $A$ に対し $\sum_lK_lAK_l^*=\sum_jL_jAL_j^*$ となるための必要十分条件は、
あるユニタリ $u\in M_n$ により $L_j=\sum_lu_{jl}K_l$ となることである。

*証明.* $L_j=\sum_lu_{jl}K_l$ ならば
$\sum_jL_jAL_j^*=\sum_{l,l'}(u^*u)_{l'l}K_lAK_{l'}^*=\sum_lK_lAK_l^*$ である。逆に、等しい写像は
等しいチョイ行列をもつ：$v_l=\sum_ie_i\otimes K_le_i$、
$w_j=\sum_ie_i\otimes L_je_i$ として $\sum_lv_lv_l^*=\sum_jw_jw_j^*$ である。$v_l$ と $w_j$ を列とする行列に補題 6.6 を適用すると
$w_j=\sum_lW_{lj}v_l$ を得る。$K\mapsto\sum_ie_i\otimes Ke_i$ は線形かつ単射であるから
$L_j=\sum_lW_{lj}K_l$ であり、$u=W^T$ はユニタリである。$\square$

コードでは、$V$ の環境の添字はクラウス・テンソルの最初の添字である：

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

> **物理的意味.** あらゆるチャネルは、系と固定された状態にある環境との等長的な相互作用の後に
> 環境を無視することとして得られる。環境を基底 $e_l$ で測定すると、確率 $\operatorname{Tr}K_l\rho K_l^*$ で
> 結果 $l$ が得られ、状態は $K_l\rho K_l^*/\operatorname{Tr}(K_l\rho K_l^*)$ となる：クラウス作用素は起こりうる「量子ジャンプ」である。
> 別の基底をとれば別のジャンプが得られるが、チャネルは同じである（定理 6.7）。最小の環境の次元は
> $\operatorname{rank}J$ である。

## 6.3 トレースに関する双対

**定義.** $\operatorname{tr}_w$、$\operatorname{tr}_v$ をそれぞれ $\mathcal M$、
$\mathcal N$ 上の忠実なトレースとする。線形写像 $\Phi:\mathcal M\to\mathcal N$ の $(w,v)$ に関する **双対** とは、
$$
\operatorname{tr}_v\big(\Phi(x)\,y\big)=\operatorname{tr}_w\big(x\,\Phi^\sharp(y)\big)\qquad(x\in\mathcal M,\ y\in\mathcal N).
$$
を満たす線形写像 $\Phi^\sharp:\mathcal N\to\mathcal M$ のことである。
これは存在して一意である：$y$ を固定すると $x\mapsto\operatorname{tr}_v(\Phi(x)y)$ は線形汎関数であるから、
定理 3.2 により、ただ一つの $\rho=:\Phi^\sharp(y)$ に対して $\operatorname{tr}_w(\rho x)=\operatorname{tr}_w(x\rho)$ に等しく、
一意性により $\Phi^\sharp$ は線形となる。blunt トレースの場合には $\Phi^*$ と書く。
$z_w=\sum_cw_cP_c$、$z_v=\sum_dv_dQ_d$ とおく（中心に属し、可逆である）。

**命題 6.8.** $\Phi$ がクラウス作用素 $K^{dc}_l$ をもつならば（命題 6.4）、
$$
\Phi^\sharp(y)_c=\sum_d\frac{v_d}{w_c}\sum_lK^{dc*}_l\,y_d\,K^{dc}_l,\qquad\text{i.e.}\qquad
\Phi^\sharp(y)=z_w^{-1}\,\Phi^*(z_v\,y).
$$
特に $\Phi^\sharp$ は CP であり、クラウス作用素 $(v_d/w_c)^{1/2}K^{dc*}_l$ をもつ。

*証明.* 各ブロックにおけるトレースの巡回性により
$$
\operatorname{tr}_v(\Phi(x)y)=\sum_dv_d\sum_{c,l}\operatorname{Tr}\big(K^{dc}_lx_cK^{dc*}_ly_d\big)
=\sum_cw_c\operatorname{Tr}\Big(x_c\sum_d\frac{v_d}{w_c}\sum_lK^{dc*}_ly_dK^{dc}_l\Big),
$$
であり、一意性からこの式が得られる。$v=w=1$ のときこれは $\Phi^*(y)_c=\sum_{d,l}K^{dc*}_ly_dK^{dc}_l$ となる。
因子 $v_d$ と $1/w_c$ は $z_v$ と $z_w^{-1}$ を掛けることにあたる。$\square$

すべての $x$ に対し $\operatorname{tr}_v\Phi(x)=\operatorname{tr}_w(x)$ のとき、$\Phi$ は **$(w,v)$ に関してトレース保存** であるという。

**系 6.9.** (a) $\Phi$ が $(w,v)$ に関して TP であるための必要十分条件は $\Phi^\sharp(1)=1$、すなわち
すべての $c$ に対し $\sum_{d,l}(v_d/w_c)K^{dc*}_lK^{dc}_l=1_{k_c}$ となることである。(b) $\Phi$ が単位的であるための必要十分条件は、$\Phi^\sharp$ が
$(v,w)$ に関して TP であることである。(c) $\Phi^\sharp$ の $(v,w)$ に関する双対は $\Phi$ である。(d)
$\Phi:\mathcal M\to\mathcal M$ がセクター保存で $v=w$ ならば、すべての $w$ に対し $\Phi^\sharp=\Phi^*$ であり、
$\Phi$ が一つのトレースに関して TP であることと、すべてのトレースに関して TP であることは同値である。

*証明.* (a) $\operatorname{tr}_v\Phi(x)=\operatorname{tr}_w(x\Phi^\sharp(1))$ であり、
すべての $x$ に対し $\operatorname{tr}_w(x(\Phi^\sharp(1)-1))=0$ となるのは $\Phi^\sharp(1)=1$ のときかつそのときに限る（
$x=(\Phi^\sharp(1)-1)^*$ ととる。$\operatorname{tr}_w$ は忠実である）。あとは命題 6.8 を用いる。(b)
$\operatorname{tr}_w\Phi^\sharp(y)=\operatorname{tr}_v(\Phi(1)y)$ がすべての
$y$ に対し $\operatorname{tr}_v(y)$ に等しいのは $\Phi(1)=1$ のときかつそのときに限る。(c) トレースの性質により
$\operatorname{tr}_w(\Phi^\sharp(y)x)=\operatorname{tr}_v(y\Phi(x))$ である。一意性を用いる。(d) 現れるのは $d=c$ の組のみであり、そこでは $v_d/w_c=1$ である。$\square$

セクターを変える写像に対しては (d) は成り立たない。$\mathcal M=\mathcal N=M_1\oplus M_2$（セクター $0$：
粒子なし、セクター $1$：二つのサイト上の一つの粒子）とし、$\Phi$ を **粒子損失**
$\Phi(x)=(x_0+\operatorname{Tr}x_1)\oplus0$ とする。そのクラウス作用素は $K^{00}=1$、$K^{01}_l=\langle l|$ である。
これは blunt トレースに関して TP である。`Tr_norm` に対しては $w=v=(1,\tfrac12)$、$\Phi^\sharp(1)=1\oplus2\cdot1_2$
であり、$\operatorname{tr}_v\Phi(x)=x_0+\operatorname{Tr}x_1\ne x_0+\tfrac12\operatorname{Tr}x_1=\operatorname{tr}_w(x)$ である。
ライブラリでは、`Channel.adjoint()` は $\Phi^*$ であり、`InterSectorChannel.adjoint(trace=...)` は
両方の環で同じ名前のトレース（`'blunt'`、`'norm'`、`'tau_vN'`）を用いたときの $\Phi^\sharp$ である。`tau_vN` の重み
$1/(Ck_c)$ と $1/(C'k'_d)$ はセクターの個数にもよる。ランダムな写像
$M_2\oplus M_3\to M_1\oplus M_2\oplus M_2$ と三つのトレースについて命題 6.8 を確かめる：

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

> **物理的意味.** トレースに依存しない対象は、オブザーバブルに作用するハイゼンベルク写像
> $\Psi:\mathcal N\to\mathcal M$ である：過程の前の状態 $\omega$ は過程の後には
> $\omega\circ\Psi$ となり、$\Psi$ は CP かつ単位的である。密度はトレースに依存する
> （定理 3.2）。$\omega$ の密度を $\omega\circ\Psi$ の密度に写す写像は、選んだトレースに関する
> $\Psi$ の双対であり、系 6.9(b) によりそれらのトレースに関して TP である。環が
> 重複度をもって作用する場合、物理的なトレースの重みは $w_c=m_c$ である（第2章）。上の粒子損失は、状態を数える
> blunt トレースに関しては確率を保存するが、セクターを等しく重みづける `Tr_norm` に関しては保存しない。

## 6.4 デコヒーレンス

ライブラリの標準的なチャネルは各ブロックの内部で作用する（$p,\gamma\in[0,1]$）：

| チャネル | $\Phi(\rho)_c$ | ブロック $c$ のクラウス作用素 |
|---|---|---|
| `dephasing_channel` | $(1-p)\rho_c+p\operatorname{diag}\rho_c$ | $\sqrt{1-p}\,1$, $\sqrt p\,E_{ii}$ |
| `depolarizing_channel` | $(1-p)\rho_c+p\operatorname{Tr}(\rho_c)\,1/k_c$ | $\sqrt{1-p}\,1$, $\sqrt{p/k_c}\,E_{ij}$ |
| `center_expectation` | $\operatorname{Tr}(\rho_c)\,1/k_c$ | $\sqrt{1/k_c}\,E_{ij}$ |
| `amplitude_damping_channel` | 準位 $i>0$ は $0$ へ減衰する | $E_{00}+\sqrt{1-\gamma}\sum_{i>0}E_{ii}$, $\sqrt\gamma\,E_{0i}$ |
| `random_mixed_unitary_channel` | $\sum_jp_jU_j\rho U_j^*$ | $\sqrt{p_j}\,U_j$ |

クラウス形は直接計算で確かめられる。例えば
$\sum_{ij}\frac1kE_{ij}\rho E_{ji}=\frac1k\sum_{ij}\rho_{jj}E_{ii}$ である。五つとも TP であり、
$\gamma>0$ の振幅減衰を除いてすべて単位的である。中心の上への条件付き期待値 $E$
（命題 2.8）は $p=1$ の脱分極である。セクターの測定
$X\mapsto\sum_cP_cXP_c$ は $\mathcal M$ 上で恒等写像である（3.3節）のと異なり、$E$ は $p_c$ だけを残す。

**定理 6.10（単位的チャネルの不動点）.** $K_l\in\mathcal M$ として $\Phi(X)=\sum_lK_lXK_l^*$
が TP かつ単位的であるとする。このとき $\Phi(X)=X$ となるための必要十分条件は、$X$ がすべての $K_l$ および $K_l^*$ と可換であることである。特に
不動点全体は $\mathcal M$ のフォン・ノイマン部分環をなす。

*証明.* $X$ がすべての $K_l$ と可換ならば $\Phi(X)=X\sum_lK_lK_l^*=X$ である。逆に $\Phi(X^*)=\Phi(X)^*$
であるから、不動点の自己共役部分 $(X+X^*)/2$ と $(X-X^*)/2i$ も不動点である。そこで
$X=X^*$ とする。$Y_l=K_lX-XK_l$ とおくと
$$
\sum_lY_lY_l^*=\sum_l\big(K_lX^2K_l^*-K_lXK_l^*X-XK_lXK_l^*+XK_lK_l^*X\big)=\Phi(X^2)-X^2,
$$
であり、$\Phi$ が TP であるからそのトレースは $0$ である。$\operatorname{Tr}Y_lY_l^*=\|Y_l\|_2^2\ge0$ であるから、すべての $Y_l=0$ である。
共役をとれば、$X$ は $K_l^*$ とも可換である。このような $X$ は自己共役な集合の可換子環と $\mathcal M$ との共通部分をなす：
これは $1$ を含む $*$-環であり、したがってフォン・ノイマン環である
（定理 2.2）。$\square$

したがって、位相緩和（$p>0$）は各ブロックの対角部分をちょうど保ち、脱分極（$p>0$）と $E$ は
ちょうど中心を保つ。単位性がなければ定理は成り立たない：振幅減衰は $E_{00}$ を固定するが、
$E_{00}$ はジャンプ $E_{01}$ と可換でない（下では、その不動点の空間と可換子環はいずれもブロックごとに
次元 $1$ であるが、それぞれ $\mathbb CE_{00}$ と $\mathbb C1$ である）。

**定理 6.11（単位的チャネルはエントロピーを増加させる）.** $\Phi:\mathcal M\to\mathcal N$ を正値、
TP かつ単位的とし、$\rho$ と $\Phi(\rho)$ を固有値 $\lambda_i$ と $\mu_j$ をもつ状態とする。
このとき、ある二重確率行列 $\Delta$ により $\mu=\Delta\lambda$ となり、したがって $[0,1]$ 上の任意の凹関数 $f$ に対し
$\sum_jf(\mu_j)\ge\sum_if(\lambda_i)$ である。特に $S(\Phi(\rho))\ge S(\rho)$ かつ
$\operatorname{Tr}\Phi(\rho)^2\le\operatorname{Tr}\rho^2$ である。

*証明.* 各ブロックを対角化する：$\rho=\sum_i\lambda_i|e_i\rangle\langle e_i|$、
$\Phi(\rho)=\sum_j\mu_j|f_j\rangle\langle f_j|$。ここで $\bigoplus_c\mathbb C^{k_c}$
と $\bigoplus_d\mathbb C^{k'_d}$ の正規直交基底は、各ベクトルがそれぞれ一つの直和成分に属するようにとる。すると
$|e_i\rangle\langle e_i|\in\mathcal M$ である。このとき $\mu_j=\langle f_j,\Phi(\rho)f_j\rangle=\sum_i\Delta_{ji}\lambda_i$
であり、$\Delta_{ji}=\langle f_j,\Phi(|e_i\rangle\langle e_i|)f_j\rangle\ge0$、
$\sum_j\Delta_{ji}=\operatorname{Tr}\Phi(|e_i\rangle\langle e_i|)=1$、
$\sum_i\Delta_{ji}=\langle f_j,\Phi(1)f_j\rangle=1$ である（したがって $\Delta$ は正方行列である：全成分の和をとればわかる）。
有限形のイェンセンの不等式 $f(\sum_i\Delta_{ji}\lambda_i)\ge\sum_i\Delta_{ji}f(\lambda_i)$（第3章の演習問題 2 の解答と同様に、
凹性の定義から項数に関する帰納法で従う）により、
$\sum_jf(\mu_j)\ge\sum_j\sum_i\Delta_{ji}f(\lambda_i)=\sum_if(\lambda_i)$ である。$f(t)=-t\log t$ と
$f(t)=-t^2$ をとればよい。$\square$

セクター間の写像の古典的部分はその **遷移行列**
$T_{dc}(\rho)=\operatorname{Tr}\Phi_{dc}(\rho_c)/\operatorname{Tr}\rho_c$（$\rho_c\ne0$）である。ライブラリの
`transition_matrix()` は、すべての入力セクターの最大混合状態に対する $T=T(1)$ である。

**命題 6.12.** $\Phi:\mathcal M\to\mathcal N$ を正値かつ TP とし、$\rho$ を
$p_c=\operatorname{Tr}\rho_c>0$ なる状態、$q_d=\operatorname{Tr}\Phi(\rho)_d$ とする。このとき $T(\rho)$ は列
確率行列であり、$q=T(\rho)p$ である。中心の上への条件付き期待値 $E_{\mathcal M}$、$E_{\mathcal N}$ を用いると
$$
E_{\mathcal N}\Phi E_{\mathcal M}(\rho)=\bigoplus_d\frac{(Tp)_d}{k'_d}\,1_{k'_d},\qquad T=T(1).
$$

*証明.* $\Phi_{dc}$ は正値であるから $T_{dc}\ge0$ であり、
$\sum_dT_{dc}=\operatorname{Tr}\Phi(\iota_c\rho_c)/\operatorname{Tr}\rho_c=1$ である。次に
$q_d=\sum_c\operatorname{Tr}\Phi_{dc}(\rho_c)=\sum_cT_{dc}(\rho)p_c$ である。最後に
$E_{\mathcal M}\rho=\bigoplus_cp_c1_{k_c}/k_c$ であるから、
$\Phi(E_{\mathcal M}\rho)_d=\sum_cp_c\Phi_{dc}(1_{k_c}/k_c)$ のトレースは $\sum_cT_{dc}p_c$ であり、$E_{\mathcal N}$
はブロック $d$ をそのトレースと $1/k'_d$ の積で置き換える。$\square$

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

完全な振幅減衰のもとで、$M_2\oplus M_3$ の最大混合状態のエントロピーは
$\log5\approx1.609$ から $H(\tfrac25,\tfrac35)\approx0.673$ へ下がる：単位的でないチャネルは状態を純粋化しうる。

> **物理的意味.** 量子ビットに対し、位相緩和はコヒーレンス $\rho_{01}$ に $1-p$ を掛け、
> 占有数は変えない（純粋位相緩和、時定数 $T_\varphi$）。
> $1-\gamma=e^{-t/T_1}$ の振幅減衰は励起準位を減衰させ、コヒーレンスに
> $\sqrt{1-\gamma}=e^{-t/2T_1}$ を掛ける。ここから $1/T_2=1/(2T_1)+1/T_\varphi$ が従う。定理 6.10 は
> デコヒーレンスのない情報の代数的な形である：すべてのジャンプ作用素と可換なものは生き残る（第10章）。
> 命題 6.12 によれば、中心を通して見ると、セクター間の写像は電荷の上の古典的なマルコフ
> 行列である。

## 6.5 ランダムなチャネル

**ハール・ランダムな等長作用素** $V:\mathbb C^k\to\mathbb C^N$ とは、ライブラリ（`random_isometry`）と同様に、
$G=QR$ の因子 $Q$ のことである。ここで $G$ は $N\times k$ のジニブル行列（成分は独立で、密度
$\pi^{-1}e^{-|z|^2}$ をもつ）、$Q^*Q=1$、$R$ は対角成分が正の上三角行列である。
`random_channel(alg, r)` はこのような $V_c:\mathbb C^{k_c}\to\mathbb C^r\otimes\mathbb C^{k_c}$ を
各ブロックで独立に生成し、定理 6.5 のクラウス作用素を用いる。

**補題 6.13.** $V$ は $U(N)$ のハール・ユニタリの最初の $k$ 列と同じ分布をもつ。
特に、任意のユニタリ $U\in M_N$ に対し、$UV$ は $V$ と同じ分布をもつ。

*証明.* $Z$ を最初の $k$ 列が $G$ である $N\times N$ のジニブル行列とし、$E\in M_{N\times k}$ を
$1_N$ の最初の $k$ 列とする。定理 5.3 により、$Z$ はほとんど確実に可逆であり、$Z=Q_ZR_Z$ で
$Q_Z$ はハール分布に従う。$R_Z$ は上三角であるから、$G=ZE=(Q_ZE)(E^*R_ZE)$ は
求める型の分解であり、しかもそれはただ一つである：$G=QR=Q'R'$ ならば $R^*R=G^*G=R'^*R'$ であるから、$R'R^{-1}$
はユニタリかつ対角成分が正の上三角行列であり、したがって $1$ である（補題 5.2 の証明）。ゆえに $V=Q_ZE$ であり、
$UQ_Z$ は $Q_Z$ と同じ分布をもつから、$UV=(UQ_Z)E$ は $V$ と同じ分布をもつ。$\square$

**補題 6.14.** $u$ を $\mathbb C^N$ のランダムな単位ベクトルとし、その分布はすべての
ユニタリのもとで不変であるとする。このとき $\mathbb E[u_a\bar u_b]=\delta_{ab}/N$ であり、
$$
\mathbb E\big[u_a\bar u_b\,u_c\bar u_d\big]=\frac{\delta_{ab}\delta_{cd}+\delta_{ad}\delta_{cb}}{N(N+1)} .
$$

*証明.* $\operatorname{diag}(e^{i\theta_1},\dots,e^{i\theta_N})$ のもとでの不変性により、
これらのモーメントにはそれぞれ $e^{i(\theta_a-\theta_b)}$ と $e^{i(\theta_a-\theta_b+\theta_c-\theta_d)}$ が掛かるので、
$a=b$、あるいはそれぞれ多重集合として $\{a,c\}=\{b,d\}$ でない限り零となる。置換行列のもとでの不変性により
$\mathbb E|u_a|^2=1/N$ であり、残るのは $\alpha=\mathbb E|u_1|^4$ と $\beta=\mathbb E|u_1|^2|u_2|^2$ だけである。
第一に、$\sum_c|u_c|^2=1$ から $\alpha+(N-1)\beta=\mathbb E|u_1|^2=1/N$ を得る。第二に、$(u_1+u_2)/\sqrt2$ は
あるユニタリ $U$ に対する $Uu$ の第一座標であるから、$\mathbb E|u_1+u_2|^4=4\alpha$ である。
$s=|u_1|^2+|u_2|^2$、$z=u_1\bar u_2$ とおくと $|u_1+u_2|^4=s^2+2s(z+\bar z)+(z+\bar z)^2$ である。
不変性 $u_1\mapsto-u_1$ と $u_1\mapsto iu_1$ から $\mathbb E[s(z+\bar z)]=0$ と
$\mathbb Ez^2=0$ が得られるので、$4\alpha=(2\alpha+2\beta)+2\beta$、すなわち $\alpha=2\beta$ である。したがって
$\beta=1/(N(N+1))$、$\alpha=2\beta$ であり、これはすべての場合に上の式を与える。（ここでは $N\ge2$ を用いた。
$N=1$ のときは $|u_1|=1$ であり、両辺とも $1$ に等しい。）$\square$

**定理 6.15（ランダムなチャネルの出力純度）.** $\psi$ をサイズ
$k$ のセクターにおける純粋状態とし、$\Phi$ を `random_channel` でクラウス階数 $r$ として生成する。このとき
$$
\mathbb E\operatorname{Tr}\Phi(\psi)^2=\frac{k+r}{kr+1},\qquad\mathbb E\,\Phi(\psi)=\frac1k1_k .
$$

*証明.* $\psi$ のセクターにおいて $\Phi(\psi)=\operatorname{Tr}_E|u\rangle\langle u|$ であり、$u=V\psi$ は
単位ベクトルで、補題 6.13 によりその分布はユニタリ不変である。$\mathbb C^r\otimes\mathbb C^k$ の添字を
組 $(e,x)$ で表す。すると $\Phi(\psi)_{xy}=\sum_eu_{ex}\bar u_{ey}$ であり、その平均は $r\delta_{xy}/(kr)$ である。また
$\operatorname{Tr}\Phi(\psi)^2=\sum_{x,y,e,f}u_{ex}\bar u_{ey}u_{fy}\bar u_{fx}$ である。補題 6.14 を
$N=kr$ として用いると、第一項は $x=y$（添字の組は $kr^2$ 個）を、第二項は $e=f$（$k^2r$ 個）を強いるので、
平均は $(kr^2+k^2r)/(kr(kr+1))$ である。$\square$

$r=1$ のときチャネルはユニタリであり、純度は $1$ である。$r\to\infty$ のとき純度は
最大混合状態の純度 $1/k$ に近づく。証明は $u$ の分布のユニタリ不変性しか用いていないので、
ランダムな純粋状態 $u=g/\|g\|$（$g$ はガウス的）の縮約状態の平均純度に関するルブキンの結果も与える。
これは第7章の定理 7.11 でガウス積分を用いて証明される。補題 6.14 の数
$\alpha,\beta$ は、ハール・ユニタリの成分に関する系 5.11 の数と同じである。

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

> **物理的意味.** ランダムなチャネルは、次元 $r$ の環境に由来する一般的で構造のない雑音のモデルである。
> すでに $r=4$ で、純粋なクトリット状態の純度は平均して $7/13\approx0.54$ まで下がる。
> 雑音がなければ $1$、最小値は $1/3$ である。超選択された系では、雑音は各セクターで
> 別々に作用し、問題となる次元は全空間の次元ではなくセクターの次元である。

## 6.6 比較の対象

QuTiP（`to_choi`, `to_kraus`, `to_super`, `to_stinespring`）、Qiskit（`qiskit.quantum_info`：
クラス `Kraus`, `Choi`, `SuperOp`, `Stinespring` と `random_quantum_channel`。Qiskit Aer の
`amplitude_damping_error` などの雑音モデル）、toqito（`toqito.channel_ops`, `toqito.channel_props`）および
QETLAB（`ChoiMatrix`, `KrausOperators`, `ApplyMap`, `RandomSuperoperator`）は、6.1–6.2節の表現の間の変換、
完全正値性の判定、そして（その多くは）ランダムなチャネルの生成を、
全行列環の間の写像 $M_n\to M_m$ について行う。
それらの規約は互いに異なる（列方向か行方向の積み重ね、チョイ行列の因子の順序）ので、
数値を比べる前に規約を比べること。$\mathcal M$ 上のセクター保存チャネルは、それらのパッケージでは
$M_D$（$D=\sum_ck_c$）上のチャネルとして扱え、その超演算子のサイズはセクターごとの $k_c^2\times k_c^2$
ではなく $D^2\times D^2$ となる。ライブラリが付け加えるのはブロック構造そのものである：セクターごとに、
バッチとして格納されたクラウス作用素、異なる環の間の写像（`InterSectorChannel`）、
第3章の各トレースに関する双対、そしてセクター間の遷移行列である。

## 演習問題

1. （ワーナー–ホレボ）$M_k$（$k\ge2$）上の $\Phi(\rho)=\big(\operatorname{Tr}(\rho)1-\rho^T\big)/(k-1)$ は、
   転置から作られているにもかかわらず、クラウス階数 $k(k-1)/2$ の単位的チャネルであることを示せ。
   `Channel.from_superoperator` で確かめよ。
2. $0<p\le1$ に対し、$M_k$ 上の位相緩和チャネルのクラウス階数は $k$ であり、ライブラリが格納する
   クラウス作用素の個数より一つ少ないことを示せ。
3. $X$ を $\mathbb C^k$ 上の巡回シフト、$Z=\operatorname{diag}(1,\omega,\dots,\omega^{k-1})$、
   $\omega=e^{2\pi i/k}$ とする。$\frac1{k^2}\sum_{a,b=0}^{k-1}X^aZ^b\rho Z^{-b}X^{-a}=\operatorname{Tr}(\rho)1/k$ を示せ：
   $M_k$ の中心の上への条件付き期待値は混合ユニタリ・チャネルである。
4. 中間の環で同じトレースを用いるとき $(\Psi\circ\Phi)^\sharp=\Phi^\sharp\circ\Psi^\sharp$ であることを示し、
   ライブラリの三つのトレースについて確かめよ。

<details>
<summary>解答</summary>

1. $F$ を入れ替え作用素とすると、$J(\theta)=\sum_{ij}E_{ij}\otimes E_{ji}=F$、$J(\rho\mapsto\operatorname{Tr}(\rho)1)=1\otimes1$
   であるから、$J(\Phi)=(1-F)/(k-1)=2P_-/(k-1)$ である。ここで $P_-=(1-F)/2$ は
   反対称ベクトルの上への射影である（$F=F^*$、$F^2=1$）。よって $J\ge0$ で階数は $k(k-1)/2$ であり（定理 6.2、
   系 6.3）、$\operatorname{Tr}\Phi(\rho)=\operatorname{Tr}\rho$、$\Phi(1)=1$ は直接確かめられる。

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

2. $J=(1-p)|\Omega\rangle\langle\Omega|+p\sum_iE_{ii}\otimes E_{ii}$ は
   $e_i\otimes e_i$ の張る空間の外では零であり、この $k$ 次元空間上では行列 $(1-p)\mathbf 1\mathbf 1^T+p1$
   （$\mathbf 1$ はすべての成分が $1$ のベクトル）であって、その固有値は $p$ と $(1-p)k+p$ で、すべて正である。よって
   $\operatorname{rank}J=k$ であるが、ライブラリは $1+k_{\max}$ 個の作用素を格納する。

```python
ranks = [int((torch.linalg.eigvalsh(J_c[0]) > 1e-10).sum()) for J_c in zoo['dephasing'].choi()]
print(ranks, zoo['dephasing'].kraus_rank)                 # [2, 3] versus 4 stored operators
assert ranks == alg.k_factors
```

3. $\operatorname{Tr}((X^aZ^b)^*X^{a'}Z^{b'})=k\,\delta_{aa'}\delta_{bb'}$ である：$a=a'$ でなければ
   $Z^{-b}X^{a'-a}Z^{b'}$ の対角成分は零であり、$\operatorname{Tr}Z^{b'-b}=\sum_j\omega^{j(b'-b)}=k\delta_{bb'}$ である。よって
   $B_{ab}=X^aZ^b/\sqrt k$ は $\operatorname{Tr}(A^*B)$ に関する $M_k$ の正規直交基底であり、任意の
   二つの正規直交基底はユニタリで結ばれる。定理 6.7 により、$\sum_mB_m\rho B_m^*$ はすべての正規直交基底について
   同じである。基底 $E_{ab}$ に対しては
   $\sum_{ab}E_{ab}\rho E_{ba}=\operatorname{Tr}(\rho)1$ である。$k$ で割ればよい。

```python
Xs = torch.roll(torch.eye(k, dtype=dt), 1, dims=0)
Zs = torch.diag(torch.exp(2j * math.pi * torch.arange(k, dtype=torch.float64) / k)).to(dt)
W = [torch.linalg.matrix_power(Xs, a) @ torch.linalg.matrix_power(Zs, b) / k for a in range(k) for b in range(k)]
twirl = channels.kraus_channel([M3.operator(Wab[None, None]) for Wab in W])
r3op = states.random_density_matrix(M3)
assert close(twirl(r3op), channels.center_expectation(M3)(r3op))
```

4. すべての $x$ と $z$ に対し $\operatorname{tr}(\Psi(\Phi(x))z)=\operatorname{tr}(\Phi(x)\Psi^\sharp(z))=\operatorname{tr}(x\,\Phi^\sharp(\Psi^\sharp(z)))$
   である。双対の一意性から主張が従う。

```python
Psi2 = channels.random_inter_sector_channel(alg_out, alg, kraus_rank=2)
for t in trace:
    assert close((Psi2 @ Psi).adjoint(trace=t)(x), Psi.adjoint(trace=t)(Psi2.adjoint(trace=t)(x)))
```

</details>

**次へ：**[第7章](07_tensor_products.md)ではテンソル積とエンタングルメントを扱う：部分トレース——
異なる環の間のチャネルであり、その双対は埋め込み $A\mapsto A\otimes1$ である——、ペイジの
公式、そして電荷の融合。
