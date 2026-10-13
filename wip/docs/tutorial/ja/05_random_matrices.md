[English](../05_random_matrices.md) · [Español](../es/05_random_matrices.md) · [Français](../fr/05_random_matrices.md) · [Deutsch](../de/05_random_matrices.md) · [中文](../zh/05_random_matrices.md) · **日本語** · [Русский](../ru/05_random_matrices.md)

# 5. 環の中のランダム行列

ランダム行列理論は「典型的な」作用素がどのような姿をしているかを問う。環 $\bigoplus_cM_{k_c}$ においては、その答えはセクターを尊重しなければならない：自然なランダム・ユニタリは大きな行列環 $M_D$ のものではなく環のものであり、自然なランダム・ハミルトニアンは互いに独立なブロックをもつ。本章では環のユニタリ群上のハール測度を構成し、その1次と2次のモーメントを厳密に計算し、そのうえでセクターを無視するとスペクトル統計で何がおかしくなるかを示す。

**本章で学ぶこと**：ガウス行列の QR 分解がなぜハール・ユニタリを与えるか（そしてなぜ位相を固定しなければならないか）、環のユニタリ群にわたる平均が中心への条件付き期待値であること、2次のワインガルテン公式、ウィグナー推測と間隔比、そしてセクターが混ざったスペクトルで計算した準位統計がなぜ独立な準位のそれのように見えるか。

本章を通じて、$\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$、$D=\sum_ck_c$ とし、$P_c$ は
極小中心射影、$E$ は命題 2.8 の中心への条件付き期待値とする。$U(k)$ は $k\times k$ ユニタリ行列のなす群であり、
$U(\mathcal M)=\{U\in\mathcal M:U^*U=1\}$ は $\mathcal M$ のユニタリ群である。確率変数は期待値 $\mathbb E$ をもつ確率空間上に定義され、ランダム行列の期待値は成分ごとにとる。*標準複素ガウス* 変数とは、実部と虚部が独立で
$\sim N(0,\tfrac12)$ であり、したがって $\mathbb E|z|^2=1$ となるものである（複素
dtype に対して `torch.randn` が生成するのはこれである）。*ジニブル行列* とは、成分が独立な標準複素ガウス変数である行列である。

## 5.1 ハール測度

$G$ を行列のコンパクト群（ある $U(n)$ の閉部分群）とする。$G$ 上のボレル確率測度
$\mu$ が **左不変** であるとは、すべての $g\in G$ とボレル集合
$B$ に対し $\mu(gB)=\mu(B)$ となること、同値なことに、すべての有界ボレル関数 $f$ に対し $\int f(gh)\,d\mu(h)=\int f\,d\mu$ となることである。
**右不変** は $Bg$ と $f(hg)$ を用いて同様に定義する。

**補題 5.1（一意性）.** $\mu$ が $G$ 上の左不変ボレル確率測度、$\nu$ が右不変ボレル確率測度ならば、$\mu=\nu$ である。したがって $G$ 上の左不変ボレル確率測度は高々一つであり、存在すればそれは右不変でもあり、$g\mapsto g^{-1}$ のもとで不変である。

*証明.* $f:G\to\mathbb R$ を連続（したがって有界）とする。フビニの定理により
$\int\!\big(\int f(gh)\,d\mu(h)\big)d\nu(g)=\int\!\big(\int f(gh)\,d\nu(g)\big)d\mu(h)$ である。左辺の内側の
積分は左不変性により $\int f\,d\mu$ であり、右辺の内側の積分は右不変性により $\int f\,d\nu$ である。
よって $\int f\,d\mu=\int f\,d\nu$ である。距離空間上の二つのボレル確率測度は、連続関数の積分がすべて等しければ一致する：連続関数
$\max(0,1-n\operatorname{dist}(\cdot,K))$ は閉集合 $K$ の指示関数に単調減少しながら収束するので、二つの測度は閉集合上で一致する（優収束定理）。閉集合はボレル集合族を生成し、共通部分について閉じている
（ディンキンの $\pi$–$\lambda$ 定理）。$\mu_1,\mu_2$ が左不変ならば、
$(Bg)^{-1}=g^{-1}B^{-1}$ であるから $\nu(B)=\mu_2(B^{-1})$ は右不変である。したがって $\mu_1=\nu=\mu_2$ であり、
$\mu_1$ は右不変かつ逆元をとる操作で不変である。$\square$

**定義.** $G$ の **ハール測度** とは、その左不変ボレル確率測度のことである。この分布に従う
ランダムな元を **ハール分布に従う** という。

任意のコンパクト群はハール測度をもつ（G. B. Folland, *A Course in Abstract Harmonic
Analysis*, 2nd ed., CRC Press, 2016, Section 2.2 を参照；ここでは証明なしに述べる）。$U(k)$ についてはこの一般的な定理は必要ない：定理 5.3 が測度を構成する。

**補題 5.2（QR 分解）.** 可逆な $Z\in M_k(\mathbb C)$ はすべて、$Q$ をユニタリ、$R$ を対角成分が正の上三角行列として $Z=QR$ と書ける。$Q$ と $R$ は一意であり、
$Z$ に連続的に依存する。

*証明.* $Z$ の列 $z_1,\dots,z_k$ は線形独立であり、これにグラム・シュミットの直交化を適用する：
$w_j=z_j-\sum_{i<j}\langle q_i,z_j\rangle q_i$ として $q_j=w_j/\|w_j\|$。ここで
$z_j\notin\operatorname{span}\{z_1,\dots,z_{j-1}\}=\operatorname{span}\{q_1,\dots,q_{j-1}\}$ であるから $w_j\neq0$ である。すると
$z_j=\sum_{i<j}\langle q_i,z_j\rangle q_i+\|w_j\|q_j$、すなわち $R_{jj}=\|w_j\|>0$ として $Z=QR$ であり、
各段階は可逆行列の上で連続である。$Q_1R_1=Q_2R_2$ ならば、
$T=Q_2^*Q_1=R_2R_1^{-1}$ はユニタリかつ対角成分が正の上三角行列である。その逆行列
$T^*$ は下三角であり、かつ上三角でもある（上三角行列の逆行列は上三角である）から、$T$ は対角行列である。対角成分が正の対角ユニタリ行列は $1$ である。$\square$

**定理 5.3（ガウス行列からのハール・ユニタリ）.** $Z$ を $k\times k$
ジニブル行列とする。このとき $Z$ は確率1で可逆であり、補題 5.2 の因子 $Q$ は
$U(k)$ 上でハール分布に従う。さらに、$Z=Q'R'$ が $Q'$ ユニタリ、$R'$
上三角である任意の QR 分解（例えば LAPACK が返すもので、その対角成分は実数で符号は正負どちらもありうる）ならば、
$\Lambda=\operatorname{diag}(R'_{ii}/|R'_{ii}|)$ として $Q=Q'\Lambda$ はハール分布に従う。

*証明.* *可逆性.* $\det Z$ は $Z$ の $2k^2$ 個の実座標についての零でない多項式であり、
零でない実多項式の零点集合はルベーグ測度零である（変数の個数に関する帰納法：
ある $p_j\neq0$ により $p=\sum_jp_j(x')x_n^j$ と書く。$x'$ が $p_j$ の測度零の零点集合の外にあれば、根 $x_n$ は有限個しかない。フビニの定理により結論を得る）。$Z$ は密度をもつから、
$P(\det Z=0)=0$ である。

*ガウス分布の不変性.* $M_k(\mathbb C)\cong\mathbb R^{2k^2}$ 上のルベーグ測度に関する $Z$ の密度は
$\pi^{-k^2}e^{-\|Z\|_2^2}$ である。$V\in U(k)$ に対し、写像
$Z\mapsto VZ$ は実線形で $\|Z\|_2$ を保つ（1.4節）から、
ヤコビアン $1$ の $\mathbb R^{2k^2}$ の直交写像であり、密度を保つ：$VZ$ は $Z$ と同じ分布をもつ。

*同変性.* $Z=QR$ ならば、$VZ=(VQ)R$ は $VZ$ の補題 5.2 の分解であるから、
$Q(VZ)=VQ(Z)$ である。したがって任意のボレル集合 $B\subseteq U(k)$ に対し
$P(VQ(Z)\in B)=P(Q(VZ)\in B)=P(Q(Z)\in B)$ となる：$Q(Z)$ の分布は左不変であり、
補題 5.1 によりそれはハール測度である（これは $U(k)$ のハール測度の存在の証明にもなっている）。
最後に $Q'R'=(Q'\Lambda)(\Lambda^*R')$ であり、$\Lambda^*R'$ は対角成分
$|R'_{ii}|>0$ の上三角行列であるから、一意性により $Q'\Lambda=Q(Z)$ である。$\square$

実ガウス行列に対して同じ証明を行えば、直交群 $O(k)$ のハール測度が得られる。
補正 $\Lambda$ がなければ、結果はハール分布に *従わない*：LAPACK の
$R'_{ii}$ の符号に関する規約が偏りを生み、それは下のコードで確認できる
（F. Mezzadri, *How to generate random matrices from the classical compact groups*, Notices AMS 54
(2007), 592–604）。

**系 5.4（環のハール測度）.** $U(\mathcal M)=\{\bigoplus_cU_c:U_c\in U(k_c)\}\cong\prod_cU(k_c)$
であり、$U=\bigoplus_cU_c$ が $U(\mathcal M)$ 上でハール分布に従うための必要十分条件は、ブロック $U_c$ が独立であり、
各 $U_c$ が $U(k_c)$ 上でハール分布に従うことである。

*証明.* $U^*U=\bigoplus_cU_c^*U_c$ であるから、$U$ がユニタリであることと各ブロックがユニタリであることは同値である。$U_c$ が
独立でハール分布に従うとする。$V=\bigoplus_cV_c\in U(\mathcal M)$ に対し、$VU$ のブロック $V_cU_c$ は
独立であり（独立な変数の関数）、$V_cU_c$ は $U_c$ と同じ分布をもつから、
$VU$ は $U$ と同じ分布をもつ。したがって $U$ の分布は左不変であり、補題 5.1 によりそれはコンパクト群 $U(\mathcal M)$ のハール
測度である。逆に、一意性により、ハール分布に従う任意の $U$ はこの分布をもつ。$\square$

群 $U(\mathcal M)$ は $U(D)$ よりはるかに小さい：$C\ge2$ のとき、$U(D)$ のハール・ユニタリが
$\mathcal M$ に属する確率は零である（定理 5.3 と補題 5.2 のグラム・シュミットによる構成により、
その第1列はガウス・ベクトル $z$ により $z/\|z\|$ と表されるから、その成分はほとんど確実にすべて非零である）。$SU(k)$ のハール測度も同じ考え方で得られる：$U$ が
$U(k)$ 上でハール分布に従い $\theta=\arg\det U$ ならば $V=e^{-i\theta/k}U\in SU(k)$ であり、$W\in SU(k)$ に対し
$\det(WU)=\det U$ から $V(WU)=WV(U)$ を得る。$WU$ は $U$ と同じ分布をもつから、$V$ の分布は
$SU(k)$ 上で左不変であり、補題 5.1 により $SU(k)$ のハール測度である。

**系 5.5（成分の低次モーメント）.** $U$ が $U(k)$ 上でハール分布に従うならば、
$\mathbb EU_{ij}=0$ かつ $\mathbb E|U_{ij}|^2=1/k$ である。

*証明.* $-1\in U(k)$ であるから、$-U$ は $U$ と同じ分布をもち、$\mathbb EU=-\mathbb EU$ である。置換
行列 $P$ に対し $PU$ は $U$ と同じ分布をもつから、$\mathbb E|U_{ij}|^2$ は $i$ によらない。
$\sum_i|U_{ij}|^2=1$ であるから、それは $1/k$ に等しい。$\square$

ライブラリは `random_unitary_operator` により $U(\mathcal M)$ 上のハール・ユニタリを生成する。ブロックサイズごとに一回のバッチ化された QR
分解で、すべてのセクターとサンプルを一度に処理する（GPU 上ではグラム・シュミットを用いるが、
補題 5.2 の一意性により同じ $Q$ が得られる）。測度 `'haar'` は各ブロックにおける定理 5.3 であり、
複素環に対しては `'haar_su'` が上のように $e^{-i\theta/k}$ を掛ける。`'coe'` と
`'cse'` はダイソンの円形アンサンブルである（前者は演習問題 3）。系 5.4 と 5.5、
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

補正がなければ $k=2$ に対し $\mathbb E\operatorname{Re}\operatorname{Tr}Q\approx-0.84$ であるが、
ハール・ユニタリでは $\mathbb EU=0$ である。

## 5.2 モーメントとワインガルテン計算

**補題 5.6.** $T\in M_k(\mathbb C)$ がすべての対角ユニタリ行列およびすべての置換
行列と可換ならば、$T\in\mathbb C1$ である。特に $U(k)$ の可換子環は $\mathbb C1$ である。

*証明.* 相異なる $e^{i\theta_j}$ をもつ $\Delta=\operatorname{diag}(e^{i\theta_1},\dots,e^{i\theta_k})$ に対し、
$\Delta T=T\Delta$ から $(e^{i\theta_i}-e^{i\theta_j})T_{ij}=0$ を得るので、$T$ は
対角行列である。置換行列 $P_\sigma e_j=e_{\sigma(j)}$ に対し、$P_\sigma TP_\sigma^{-1}=T$ から
$T_{\sigma(j)\sigma(j)}=T_{jj}$ を得るので、対角成分はすべて等しい。$\square$

**命題 5.7（ハール・ツワールは中心への条件付き期待値である）.** $U$ を
$U(\mathcal M)$ 上でハール分布に従うとする。すべての $X\in\mathcal M$ に対し
$$
\mathbb E\,[UXU^*]=E(X)=\bigoplus_c\frac{\operatorname{Tr}X_c}{k_c}\,1_{k_c}.
$$

*証明.* $(UXU^*)_c=U_cX_cU_c^*$ であり、$U_c$ は $U(k_c)$ 上でハール分布に従う（系 5.4）。
$T_c=\mathbb E[U_cX_cU_c^*]$ とおく（成分は有界である）。$V\in U(k_c)$ に対し $VU_c$ は
$U_c$ と同じ分布をもつから、$VT_cV^*=\mathbb E[(VU_c)X_c(VU_c)^*]=T_c$ である。補題 5.6 により $T_c=\lambda1$ であり、
$k_c\lambda=\operatorname{Tr}T_c=\mathbb E\operatorname{Tr}(U_cX_cU_c^*)=\operatorname{Tr}X_c$ である。$\square$

正規化は命題 2.8 および `channels.center_expectation` のものと同じである。
$X\in\mathcal M$ が $U(\mathcal M)$ と可換ならば $UXU^*=X$ であるから、$X=E(X)$ は中心に属する。ハール
平均は第2章の演習問題 3 の群平均の連続版である。代わりに $U$ が
$U(D)$ 上でハール分布に従うならば、同じ証明（$C=1$）により $\mathbb E[UXU^*]=\frac{\operatorname{Tr}X}{D}1$ となる：
セクターの重みは失われ、$UXU^*$ は $\mathcal M$ に属しさえしない。

> **物理的意味.** 環のハール・ランダムなユニタリは、超選択則を尊重する
> ランダムな時間発展である。それについて平均をとると、各セクター内部のすべてがかき混ぜられ、
> ちょうど古典的部分、すなわち各電荷の重みだけが残る。これが量子
> 情報理論における「ツワール」であり、それは第2章の中心への条件付き期待値である。

2次のモーメントには、$\mathbb C^k\otimes\mathbb C^k$ 上の $\{V\otimes V\}$ の可換子環が必要である。
$F$ を **スワップ**、$F(x\otimes y)=y\otimes x$ とする。これは自己共役かつユニタリであり、
$\operatorname{Tr}1=k^2$、$\operatorname{Tr}F=\sum_{i,j}\langle e_i\otimes e_j,e_j\otimes e_i\rangle=k$、$F^2=1$ である。

**補題 5.8.** $k\ge2$ とする。$\mathbb C^k\otimes\mathbb C^k$ 上の行列 $X$ がすべての $V\in U(k)$ に対して $V\otimes V$
と可換であるための必要十分条件は、$a,b\in\mathbb C$ により $X=a1+bF$ と書けることである。また $1,F$ は線形独立である。

*証明.* $F(V\otimes V)(x\otimes y)=Vy\otimes Vx=(V\otimes V)F(x\otimes y)$ であるから、$1$ と $F$ は
すべての $V\otimes V$ と可換である。逆に $X$ がそれらと可換であるとし、
$X_{ij,pq}=\langle e_i\otimes e_j,X(e_p\otimes e_q)\rangle$ と書く。
*対角ユニタリ:* $\Delta\otimes\Delta$ は $e_i\otimes e_j$ に $e^{i(\theta_i+\theta_j)}$ を掛けるから、
$(e^{i(\theta_i+\theta_j)}-e^{i(\theta_p+\theta_q)})X_{ij,pq}=0$ である。多重集合 $\{i,j\}$ と
$\{p,q\}$ が異なれば、ある添字 $m$ が前者に後者より $n\in\{\pm1,\pm2\}$ 回多く現れる。
$\theta_m=\pi/2$、他のすべての $\theta=0$ とすると、括弧は $e^{i(\theta_p+\theta_q)}(e^{in\pi/2}-1)\neq0$ となるから、
$X_{ij,pq}=0$ である。したがって $X(e_i\otimes e_i)=a_ie_i\otimes e_i$ であり、$i\neq j$ に対しては
$X(e_i\otimes e_j)=b_{ij}\,e_i\otimes e_j+c_{ij}\,e_j\otimes e_i$ である。
*置換:* $X(P_\sigma\otimes P_\sigma)=(P_\sigma\otimes P_\sigma)X$ を $e_i\otimes e_j$ に適用すると
$a_{\sigma(i)}=a_i$、$b_{\sigma(i)\sigma(j)}=b_{ij}$、$c_{\sigma(i)\sigma(j)}=c_{ij}$ を得る。置換は相異なる添字の順序対の上に
推移的に作用するから、$a_i=a'$、$b_{ij}=b'$、$c_{ij}=c'$ である。よって
$i\ne j$ なる $e_i\otimes e_j$ の張る空間の上で $X=b'1+c'F$ である。
*回転:* $R$ を $Re_1=(e_1+e_2)/\sqrt2$、$Re_2=(-e_1+e_2)/\sqrt2$、$m\ge3$ に対し $Re_m=e_m$ となるユニタリとする。
すると $(R\otimes R)(e_1\otimes e_1)=\frac12(e_1\otimes e_1+e_1\otimes e_2+e_2\otimes e_1+e_2\otimes e_2)$ であり、
$X(R\otimes R)(e_1\otimes e_1)=\frac12\big(a'e_1\otimes e_1+(b'+c')(e_1\otimes e_2+e_2\otimes e_1)+a'e_2\otimes e_2\big)$
は $(R\otimes R)X(e_1\otimes e_1)=a'(R\otimes R)(e_1\otimes e_1)$ に等しくなければならない。
$e_1\otimes e_2$ の係数を比較すると $a'=b'+c'$ を得るが、これは $b'1+c'F$ の $e_i\otimes e_i$ 上での値である。よって $X=b'1+c'F$ である。
最後に、$1$ と $F$ は $e_1\otimes e_2$ をそれぞれ線形独立なベクトル $e_1\otimes e_2$ と $e_2\otimes e_1$ に写す。$\square$

**系 5.9（2次のシューア・ワイル双対性）.** $k\ge2$ とし、$P_\pm=\frac12(1\pm F)$ とおく。
これらは対称部分空間と反対称部分空間
$\mathrm{Sym}=\{x:Fx=x\}$、$\mathrm{Alt}=\{x:Fx=-x\}$ の上への射影であり、次元はそれぞれ $k(k+1)/2$ と $k(k-1)/2$ である。
$U(k)$ の表現 $V\mapsto V\otimes V$ は $\mathrm{Sym}$ 上でも $\mathrm{Alt}$ 上でも既約であり、
二つの表現は同値でない。

*証明.* $F=F^*$ と $F^2=1$ により、$P_\pm$ は $P_++P_-=1$ を満たす射影であり、値域は $\mathrm{Sym}$、
$\mathrm{Alt}$ である。その次元は $\operatorname{Tr}P_\pm=(k^2\pm k)/2$ である。
$\{V\otimes V\}$ の張る空間 $\mathcal A$ は $*$-環であり（$1\otimes1$、積 $VW\otimes VW$、随伴 $V^*\otimes V^*$）、
補題 5.8 によりその可換子環は $\mathcal A'=\operatorname{span}\{1,F\}=\operatorname{span}\{P_+,P_-\}$ である。
$K\subseteq\mathrm{Sym}$ がすべての $V\otimes V$ で不変であるとする。補題 2.1 によりその射影は $\mathcal A'$ に属するから、
$x,y\in\{0,1\}$ により $xP_++yP_-$ と書ける（射影だからである）。それは $\mathrm{Alt}$ 上で零であるから $y=0$ であり、
$K=0$ または $K=\mathrm{Sym}$ である。同じ議論が $\mathrm{Alt}$ にも適用できる。
$T:\mathrm{Sym}\to\mathrm{Alt}$ が二つの表現の間の絡作用素ならば、$TP_+$（$\mathbb C^k\otimes\mathbb C^k$ 上の
作用素とみなす）は $\mathcal A$ と可換であるから、$TP_+=xP_++yP_-$ である。それは
$\mathrm{Alt}$ 上で零であるから $y=0$ であり、その値域は $\mathrm{Alt}\cap\mathrm{Sym}=0$ に含まれるから $T=0$ である。$\square$

**定理 5.10（2次のワインガルテン公式）.** $k\ge2$ とし、$U$ を $U(k)$ 上でハール分布に従うとし、
$X$ を $\mathbb C^k\otimes\mathbb C^k$ 上の行列とする。このとき
$$
\mathbb E\,[(U\otimes U)X(U\otimes U)^*]=\alpha\,1+\beta F,\qquad
\alpha=\frac{\operatorname{Tr}X-\operatorname{Tr}(FX)/k}{k^2-1},\quad
\beta=\frac{\operatorname{Tr}(FX)-\operatorname{Tr}X/k}{k^2-1}.
$$
同値なことに、群 $S_2=\{e,\tau\}$ 上で $\mathrm{Wg}(e)=\frac1{k^2-1}$、$\mathrm{Wg}(\tau)=-\frac1{k(k^2-1)}$ とすると
$$
\mathbb E\big[U_{i_1j_1}U_{i_2j_2}\overline{U_{i'_1j'_1}U_{i'_2j'_2}}\big]
=\sum_{\sigma,\pi\in S_2}\delta_{i_1i'_{\sigma(1)}}\delta_{i_2i'_{\sigma(2)}}\,
\delta_{j_1j'_{\pi(1)}}\delta_{j_2j'_{\pi(2)}}\,\mathrm{Wg}(\sigma\pi^{-1}).
$$

*証明.* $T=\mathbb E[(U\otimes U)X(U\otimes U)^*]$ とおく。命題 5.7 と同様に、$T$ は
すべての $V\otimes V$ と可換であるから、補題 5.8 により $T=\alpha1+\beta F$ である。トレースをとると
$\operatorname{Tr}T=\operatorname{Tr}X$ であり、また $F$ は $U\otimes U$ と可換であるから
$\operatorname{Tr}(FT)=\mathbb E\operatorname{Tr}((U\otimes U)^*F(U\otimes U)X)=\operatorname{Tr}(FX)$ である。
これより $k^2\alpha+k\beta=\operatorname{Tr}X$ と $k\alpha+k^2\beta=\operatorname{Tr}(FX)$ を得る。これは
行列式 $k^2(k^2-1)\neq0$ の連立一次方程式であり、その解は述べたとおりである。
成分については $X=|e_{j_1}\otimes e_{j_2}\rangle\langle e_{j'_1}\otimes e_{j'_2}|$ ととる：すると
$\langle e_{i_1}\otimes e_{i_2},T\,e_{i'_1}\otimes e_{i'_2}\rangle$ が左辺であり、$\operatorname{Tr}X=\delta_{j_1j'_1}\delta_{j_2j'_2}$
（$\pi=e$ の項）、$\operatorname{Tr}(FX)=\delta_{j_1j'_2}\delta_{j_2j'_1}$（$\pi=\tau$ の項）であり、
$1$、$F$ はそれぞれ $\delta_{i_1i'_1}\delta_{i_2i'_2}$（$\sigma=e$ の項）と
$\delta_{i_1i'_2}\delta_{i_2i'_1}$（$\sigma=\tau$）を寄与する。項をまとめると
$\alpha=\sum_\pi\mathrm{Wg}(e\pi^{-1})[\pi]$、$\beta=\sum_\pi\mathrm{Wg}(\tau\pi^{-1})[\pi]$ となる。
ここで $[\pi]$ は $\pi$ に対応する $j$ のデルタの積を表す。$\square$

行列 $\big(\mathrm{Wg}(\sigma^{-1}\pi)\big)_{\sigma,\pi}$ は、
$\langle A,B\rangle=\operatorname{Tr}(A^*B)$ に関する $1$ と $F$ のグラム行列
$\begin{pmatrix}k^2&k\\k&k^2\end{pmatrix}$ の逆行列である。一般化されるのはこの形である。

**系 5.11.** $k\ge2$ に対し
$$
\mathbb E|U_{11}|^4=\frac2{k(k+1)},\qquad
\mathbb E|U_{11}|^2|U_{12}|^2=\frac1{k(k+1)},\qquad
\mathbb E\big[U_{11}U_{22}\overline{U_{12}U_{21}}\big]=-\frac1{(k-1)k(k+1)}.
$$

*証明.* 定理 5.10 において、(i) すべての添字が $1$ のとき：デルタの積はすべて $1$ であるから、
モーメントは $2(\mathrm{Wg}(e)+\mathrm{Wg}(\tau))=2\frac{k-1}{k(k^2-1)}$ である。(ii) $i=(1,1)$、$j=(1,2)$、
$i'=(1,1)$、$j'=(1,2)$：$\pi=e$ の項だけが残り、$\mathrm{Wg}(e)+\mathrm{Wg}(\tau)$ を与える。
(iii) $i=(1,2)$、$j=(1,2)$、$i'=(1,2)$、$j'=(2,1)$：$\sigma=e$、$\pi=\tau$ の項だけが残り、
$\mathrm{Wg}(\tau)$ を与える。$\square$

最後のモーメントは負である：ハール・ユニタリの成分は独立ではなく、
相関は相対的に $1/k$ のオーダーである。環については、異なるブロックの成分は
独立であるから（系 5.4）、セクターをまたぐ混合モーメントは単一ブロックのモーメントの積に分解する。

任意の次数のモーメントについても、
$\{V^{\otimes n}\}$ の可換子環がわかれば同じ推論が適用できる。それは $n$ 個のテンソル因子の置換によって張られる
（シューア・ワイル双対性。R. Goodman, N. R. Wallach, *Symmetry, Representations, and Invariants*,
Springer, 2009, Chapter 9 を参照）。この結果は証明なしに述べる。

**定理 5.12（ワインガルテン公式；証明なしに述べる）.** $n\ge1$ とし、$U$ を $U(k)$ 上でハール分布に従うとする。このとき
$$
\mathbb E\big[U_{i_1j_1}\cdots U_{i_nj_n}\overline{U_{i'_1j'_1}\cdots U_{i'_nj'_n}}\big]
=\sum_{\sigma,\pi\in S_n}\prod_{m=1}^n\delta_{i_mi'_{\sigma(m)}}\delta_{j_mj'_{\pi(m)}}\;\mathrm{Wg}(\sigma\pi^{-1},k),
$$
ここで $\mathrm{Wg}(\cdot,k)$ は $S_n$ 上の関数で、巡回型のみに依存する。$n\le k$ のとき
行列 $\big(\mathrm{Wg}(\sigma^{-1}\pi,k)\big)_{\sigma,\pi}$ はグラム
行列 $\big(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)}\big)_{\sigma,\pi}$ の逆行列である。$n>k$ のときこの行列は特異であり、
その擬似逆行列をとる。$U$ と $\overline U$ の個数が異なるモーメントは零である。（D. Weingarten,
J. Math. Phys. 19 (1978), 999–1001; B. Collins, Int. Math. Res. Not. 2003, 953–982; B. Collins,
P. Śniady, Comm. Math. Phys. 264 (2006), 773–795.）

命題 5.7 を `channels.center_expectation` と照合し、補題 5.8 の可換子環を
厳密な線形代数で確かめ、定理 5.10 と系 5.11 を上のサンプルの $U(3)$ ブロックを用いたモンテカルロ法で確かめる：

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

## 5.3 準位統計

$H$ を固有値
$\lambda_1\le\dots\le\lambda_n$ をもつランダムな自己共役 $n\times n$ 行列とする。その **準位間隔** は $s_i=\lambda_{i+1}-\lambda_i$ である。標準的なアンサンブルが
二つある：**GOE** は $H=(G+G^T)/2$ で、$G$ は成分が独立な $N(0,1)$ の実行列である。
**GUE** は $H=(G+G^*)/\sqrt2$ で、$G$ はジニブル行列である（1.7節）。いずれにおいても、対角成分は
$N(0,1)$ に従い、対角線より上の成分の実部と虚部は $N(0,\tfrac12)$ に従い、
$H$ の密度は $e^{-\operatorname{Tr}H^2/2}$ に比例する。

**命題 5.13（ウィグナー推測）.** $n=2$、$s=\lambda_2-\lambda_1$、$S=s/\mathbb Es$ とする。
GOE（$\beta=1$）と GUE（$\beta=2$）に対し、$S$ の密度はそれぞれ
$$
p_1(S)=\frac\pi2\,S\,e^{-\pi S^2/4},\qquad p_2(S)=\frac{32}{\pi^2}\,S^2\,e^{-4S^2/\pi},
$$
であり、分布関数はそれぞれ $P(S\le t)=1-e^{-\pi t^2/4}$ と
$P(S\le t)=\operatorname{erf}(2t/\sqrt\pi)-\frac{4t}\pi e^{-4t^2/\pi}$ である。特に
$\varepsilon\to0$ のとき $P(S\le\varepsilon)$ は $\varepsilon^{\beta+1}$ のオーダーである。

*証明.* $H=\begin{pmatrix}a&b\\\bar b&c\end{pmatrix}$ に対し固有値は
$\frac{a+c}2\pm\frac12\sqrt{(a-c)^2+4|b|^2}$ であるから、$s=\sqrt{(a-c)^2+4|b|^2}$ である。*GOE:* $x=a-c$ と
$y=2b$ は独立な $N(0,2)$ 変数であり、極座標で
$P(s\le t)=\int_{x^2+y^2\le t^2}\frac{e^{-(x^2+y^2)/4}}{4\pi}\,dx\,dy=\int_0^t\frac r2e^{-r^2/4}dr=1-e^{-t^2/4}$ となる。
すると $\mathbb Es=\int_0^\infty P(s>t)\,dt=\int_0^\infty e^{-t^2/4}dt=\sqrt\pi$ であり、
$P(S\le t)=P(s\le\sqrt\pi t)=1-e^{-\pi t^2/4}$ である。これを微分すればよい。*GUE:* $x=a-c$、$y=2\operatorname{Re}b$、
$y'=2\operatorname{Im}b$ は独立な $N(0,2)$ 変数であり、球座標で
$s=\sqrt{x^2+y^2+y'^2}$ の密度は $4\pi r^2(4\pi)^{-3/2}e^{-r^2/4}=\frac{r^2e^{-r^2/4}}{2\sqrt\pi}$ である。その平均は
$\frac1{2\sqrt\pi}\int_0^\infty r^3e^{-r^2/4}dr=\frac{8}{2\sqrt\pi}=\frac4{\sqrt\pi}$ である。
$r=4S/\sqrt\pi$ を代入すると $p_2$ が得られ、$p_2$ を部分積分すると
（$a=4/\pi$ として $\int x^2e^{-ax^2}dx=-\frac{x}{2a}e^{-ax^2}+\frac1{2a}\int e^{-ax^2}dx$）
分布関数が得られる。小さい $\varepsilon$ での振る舞いは $p_\beta(S)\sim cS^\beta$ から従う。$\square$

指数は次元を数えている：実対称 $2\times2$ 行列が縮退するのは二つの実
条件（$a=c$、$b=0$）が成り立つときかつそのときに限り、複素エルミート行列では三つの条件が成り立つときかつそのときに限る。任意の $n$ に対し、固有値の同時
密度は $\prod_{i<j}|\lambda_i-\lambda_j|^\beta e^{-\sum_i\lambda_i^2/2}$ に比例し、
$n\to\infty$ でスペクトルのバルクにおける間隔分布はある極限（ゴーダン・メータ分布）に収束する。それは
推測に近いが等しくはない（M. L. Mehta, *Random Matrices*, 3rd ed., Elsevier, 2004；証明なしに述べる）。

> **物理的意味.** ウィグナーは、複雑な量子系（重い原子核）のエネルギー準位は
> 統計的に、適切な対称性をもつランダムなハミルトニアンの固有値のようにふるまうと提唱した：
> 時間反転対称性があれば実対称（$\beta=1$）、なければ複素エルミート
> （$\beta=2$）である。このとき準位は *反発* する：同じ対称性をもつ二つの準位はほとんど決して近づかない
> （フォン・ノイマンとウィグナー、1929年）。反発は量子カオスの典型的な特徴であり、独立な
> 準位（反発なし）は可積分系または局在系の特徴である（ボヒガス・ジャンノーニ・シュミット予想と
> ベリー・テイバー予想）。

間隔は局所的な準位密度に依存し、それはスペクトルに沿って変化する。間隔を
$p_\beta$ と比較するには **アンフォールディング**、すなわち計数関数の滑らかな近似 $\bar N$ による変数変換 $\lambda\mapsto\bar N(\lambda)$ が必要であり、
これにより平均間隔が $1$ になる。$\bar N$ の選び方は偏りの原因となる。
連続する間隔の比はこれを回避する。

**定義.** **間隔比** を $r_i=\min(s_i,s_{i+1})/\max(s_i,s_{i+1})\in[0,1]$ と定める
（V. Oganesyan, D. A. Huse, Phys. Rev. B 75 (2007), 155111）。これは
$\lambda\mapsto a\lambda+b$（$a>0$）で不変であるから、密度が二つの間隔にわたってほぼ一定である範囲では
アンフォールディングを必要としない。

**命題 5.14（独立な準位）.** $s,s'$ を独立で同じ指数
分布に従うとする（ポアソン過程の連続する二つの間隔のように）。このとき $r=\min(s,s')/\max(s,s')$ は
$[0,1]$ 上で $P(r\le t)=\frac{2t}{1+t}$、密度 $\frac2{(1+t)^2}$ をもち、
$$
\mathbb Er=2\ln2-1\approx0.3863 .
$$

*証明.* スケーリングにより、レートは $1$ としてよい。$x=s/s'$ に対し $P(x\le t)=\mathbb E\,P(s\le ts'\mid s')=\mathbb E(1-e^{-ts'})=1-\frac1{1+t}=\frac t{1+t}$ である。
$t<1$ に対し事象 $\{x\le t\}$ と $\{1/x\le t\}$ は互いに素であり、$x$ と $1/x$ は同じ分布をもち
（$s,s'$ は交換可能である）、$\{r\le t\}$ はそれらの和集合であるから、$P(r\le t)=\frac{2t}{1+t}$ である。よって
$\mathbb Er=\int_0^1\frac{2t}{(1+t)^2}dt=2\Big[\ln(1+t)+\frac1{1+t}\Big]_0^1=2\ln2-1$ である。$\square$

区間上の $n$ 個の独立な一様準位についても同じことが厳密に成り立つ。それらの間隔は、独立な指数変数 $e_j$ による
$(e_1,\dots,e_{n+1})/\sum_je_j$ と同じ分布をもち、したがって各
比は $s_i/s_{i+1}=e_i/e_{i+1}$ となるからである（W. Feller, *An Introduction to Probability Theory and Its
Applications*, Vol. II, Wiley, 1971, Chapter I；証明なしに述べる）。GOE と GUE については
$n\to\infty$ の極限で $\langle r\rangle=0.5307$ と $0.5996$ である（数値的な値。推測の $3\times3$ 版
は $4-2\sqrt3\approx0.5359$ と $2\sqrt3/\pi-\frac12\approx0.6027$ を与える。Y. Y. Atas,
E. Bogomolny, O. Giraud, G. Roux, Phys. Rev. Lett. 110 (2013), 084101；証明なしに述べる）。

関数 `SpinChain.level_spacing_ratio`（`chains` 内）は準位の各行の平均間隔比を返す。
命題 5.13 と 5.14、および GOE の値を確かめる：

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
mc(ratio(levels), 2 * math.log(2) - 1, "<r>, independent uniform levels")

n = 64
alg64 = TypeIAlgebra([n], [n], complex_valued=False, precision='double')
G = torch.randn(800, 1, n, n, dtype=torch.float64)
lam = alg64.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()[0]
mc(ratio(lam[:, n // 4: 3 * n // 4]), 0.5307, "<r>, GOE(64), middle half of the spectrum")
```

## 5.4 なぜセクターを混ぜてはならないか

$H=\bigoplus_cH_c\in\mathcal M$ を **独立なブロック** をもつランダムな自己共役元、
例えば GOE ブロックとする。$\mathbb C^D$ 上の作用素として、そのスペクトルは各ブロックのスペクトルの
和集合である。ブロックの内部では準位は反発する。異なるブロックの準位は反発しないことを示す。

各ブロックが有界連続な **準位密度** $\rho_c$ をもつと仮定する：
有界ボレル関数 $\varphi$ に対し $\mathbb E\sum_i\varphi(\lambda_i(H_c))=\int\varphi\,\rho_c$（したがって
$\int\rho_c=k_c$）。GOE と GUE については、これは上で引用した同時密度から従う。

**定理 5.15（セクター間の反発の不在）.** $c\neq d$ と任意の有界ボレル関数
$\psi:\mathbb R^2\to\mathbb R$ に対し
$$
\mathbb E\sum_{i,j}\psi\big(\lambda_i(H_c),\lambda_j(H_d)\big)=\iint\psi(x,y)\,\rho_c(x)\rho_d(y)\,dx\,dy .
$$
その結果、$N_\varepsilon$ を異なるブロックに属し距離が
$<\varepsilon$ である準位の対の個数とすると、$\varepsilon\to0$ のとき $\mathbb EN_\varepsilon/(2\varepsilon)\to\sum_{c<d}\int\rho_c\rho_d$ であり、
合併したスペクトルの最小間隔 $\delta$ は次を満たす。
$$
P(\delta<\varepsilon)\ \ge\ \frac{\mathbb EN_\varepsilon}{\sum_{c<d}k_ck_d}
=\frac{2\varepsilon\sum_{c<d}\int\rho_c\rho_d}{\sum_{c<d}k_ck_d}+o(\varepsilon).
$$

*証明.* 順序づけた固有値は行列の連続関数であるから、
$\Lambda_c=(\lambda_i(H_c))_i$ と $\Lambda_d$ は独立なランダム・ベクトルである。フビニの定理により
（$\Lambda_d$ について条件づけて）、$\mathbb E\sum_{i,j}\psi(\Lambda_{c,i},\Lambda_{d,j})=\mathbb E\sum_j g(\Lambda_{d,j})$ である。ここで
$g(y)=\mathbb E\sum_i\psi(\Lambda_{c,i},y)=\int\psi(x,y)\rho_c(x)\,dx$ は有界なボレル関数である。$g$ に
$\rho_d$ の定義を適用すればよい。$\psi(x,y)=1\{|x-y|<\varepsilon\}$ とすると
$$
\frac{\mathbb EN_\varepsilon}{2\varepsilon}=\sum_{c<d}\int\rho_c(x)\Big(\frac1{2\varepsilon}\int_{x-\varepsilon}^{x+\varepsilon}\rho_d(y)\,dy\Big)dx ;
$$
括弧内は連続性により $\rho_d(x)$ に収束し、$\sup\rho_d$ で抑えられ、$\rho_c$ は可積分であるから、
優収束定理により極限が得られる。最後に $N_\varepsilon\le\sum_{c<d}k_ck_d$ であるから
$\mathbb EN_\varepsilon\le\big(\sum_{c<d}k_ck_d\big)P(N_\varepsilon\ge1)$ であり、$N_\varepsilon\ge1$ ならば $\delta<\varepsilon$ である。$\square$

したがって、スペクトルが重なり合う二つのブロックに対しては $P(\delta<\varepsilon)$ は少なくとも $\varepsilon$ のオーダーであるのに対し、
$2\times2$ GOE ブロックの内部ではそれは $1-e^{-\varepsilon^2/4}\le\varepsilon^2/4$ であり
（命題 5.13、正規化していない間隔）、$n\times n$ GOE ブロックではそれは
$\varepsilon^2$ のオーダーである（同時密度から従う。ここでは証明しない）。すべてのブロックが完全な GOE であるにもかかわらず、
合併したスペクトルでは準位反発の検定が失敗する。

*多数のセクター（証明ではなく導出）.* 定常なスペクトルに対しては、長さ $s$ の区間に
準位が一つもない確率 $\mathcal G(s)$ が間隔の密度を決める：平均間隔が $1$ のとき
$p(s)=\mathcal G''(s)$ である（Mehta, *op. cit.*）。独立なスペクトルの重ね合わせでは、これらのギャップ
確率は掛け合わされる。$C$ 個のセクターがそれぞれ準位の割合 $1/C$ を担い、各々がそれ自身のアンフォールドされた単位でギャップ
確率 $\mathcal G_1$ をもつならば、
$\mathcal G_1(x)=1-x+\int_0^x(x-t)p_1(t)\,dt$ として $\mathcal G(s)=\mathcal G_1(s/C)^C$ である。$0$ で2回微分すると
$$
p(0)=1-\frac1C+\frac{p_1(0)}C=1-\frac1C ,
$$
となり、これはアンフォールドされたスペクトルに対する定理 5.15 の類似である。また $\mathcal G_1(s/C)^C\to e^{-s}$ であり、これはポアソン過程のギャップ確率である：
$C\to\infty$ のとき、合併したスペクトルは独立な準位のように見える
（N. Rosenzweig, C. E. Porter, Phys. Rev. 120 (1960), 1698–1714; M. V. Berry, M. Robnik, J. Phys. A 17
(1984), 2413–2421。ポアソン極限の厳密な形は疎な点過程の重ね合わせに関する極限定理である。
B. Grigelionis, Theory Probab. Appl. 8 (1963) を参照）。

ライブラリではブロックは最初から分かれている：`eigenvalues()` はセクターごとに一つのスペクトルを返し、
それらを混ぜるには明示的な `torch.cat` が必要である。サイズ $64$ の独立な GOE ブロック $C$ 個を
一つの環に入れ、セクターごとの平均間隔比と合併したスペクトルのそれとを比較する：

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
    assert abs(per_sector.mean() - 0.5307) < 0.01                    # every sector is GOE
assert merged_r[1] > merged_r[2] > merged_r[4] > merged_r[16]
assert abs(merged_r[16] - (2 * math.log(2) - 1)) < 0.01              # close to independent levels
```

どのセクターも $\langle r\rangle\approx0.53$ である。合併したスペクトルでは、二つのセクターで $0.42$ に下がり、
$8$ と $16$ では独立な準位の値 $0.39$ まで下がる。

同じことが物理的なモデルでも起こる。ランダム磁場中のハイゼンベルク鎖（第9章）は
磁化を保存するから、そのハミルトニアンは $\bigoplus_{N=0}^LM_{\binom LN}$ に属する。弱い乱れに対しては、
各セクター内部の準位は GOE 統計に従う（数値的に。A. Pal, D. A. Huse, Phys. Rev. B 82 (2010),
174411）。そのブロックは独立では *ない*——すべてが同じランダム磁場に依存する——ので、
定理 5.15 はそのままでは適用できない。それでも実験は同じ効果を示す。

```python
chain = SpinChain(10, complex_valued=False)                         # sectors N = 0, ..., 10
lam = chain.random_field_heisenberg(W=0.5, batch_size=20).eigenvalues()
k5 = chain.sector_dim(5)                                            # 252
r_sector = ratio(lam[5][:, k5 // 4: 3 * k5 // 4]).mean().item()
spectrum = torch.cat(lam, dim=-1).sort(dim=-1).values               # all 1024 levels mixed
r_all = ratio(spectrum[:, 256:768]).mean().item()
print(f"Heisenberg chain, L = 10: <r> in sector N = 5: {r_sector:.4f}, all sectors mixed: {r_all:.4f}")
assert r_sector > 0.49 and r_all < 0.42
```

セクター $N=5$ の内部では比は $0.52$（$0.5307$ よりわずかに小さい有限サイズの値）であり、
すべてのセクターを混ぜたスペクトルでは $0.39$ となる。

> **物理的意味.** 保存する量子数が異なる準位は相互作用しないので、
> 自由に交差する。反発するのは同じ対称性の準位だけである。したがって対称性セクターを混ぜた
> スペクトルは、各セクターのダイナミクスが完全にカオス的であってもポアソン的な統計を示し、
> 分離されていない対称性は可積分性や局在と取り違えられうる。
> 「準位統計を計算する前にすべての対称性を分離せよ」という経験則は、
> $M_D$ ではなく環 $\bigoplus_cM_{k_c}$ の中で計算すべきだという主張にほかならない。

## 5.5 比較の対象

ハール・ユニタリ：`scipy.stats.unitary_group`（定理 5.3 の補正付き QR）と `ortho_group`、
QuTiP のハール・ランダムなユニタリ、そして位相補正を施した任意の `numpy`/`torch` の QR。
ワインガルテン計算：記号計算パッケージ IntU（Mathematica；Z. Puchała, J. A. Miszczak）、RTNI
（Mathematica と Python；M. Fukuda, R. König, I. Nechita, J. Phys. A 52 (2019), 425303）、および
Python パッケージ Haarpy は任意の次数のモーメントを記号的に計算する。本章で証明した公式はその
2次の場合である。準位統計は通常数行の `numpy` で計算され、
QuSpin のような厳密対角化コードはセクターごとにスペクトルを与える。ライブラリが
付け加えるのは、環のユニタリ群を一つの対象として扱うことである：あらゆるサイズの独立なハール・ブロックをバッチ化された呼び出しで生成し、
ツワールを条件付き期待値
`channels.center_expectation` と同一視し、スペクトルをセクターごとに返す。そのため、セクターを混ぜることは決して
既定の動作にならない。

## 演習問題

1. $U$ を $U(\mathcal M)$ 上でハール分布に従うとし、$\operatorname{Tr}U=\sum_c\operatorname{Tr}U_c$ とする。
   $\mathbb E|\operatorname{Tr}U|^2=C$ かつ $\mathbb E|\operatorname{Tr}U|^4=\sum_cm(k_c)+2C(C-1)$ であることを示せ。ただし
   $m(1)=1$、$k\ge2$ に対し $m(k)=2$ である。$U(D)$（$D\ge2$）のハール・ユニタリに対する $1$ と $2$ と比較せよ。
   $\mathcal M=M_1\oplus M_2\oplus M_3$ について確かめよ。
2. 定理 5.12 を用いて、$n\le k$ に対し $\mathbb E|U_{11}|^{2n}=\frac{n!}{k(k+1)\cdots(k+n-1)}$ であることを示せ。
   $k=3$ に対し3次のワインガルテン行列を数値的に計算し、$\mathbb E|U_{11}|^6=\frac1{10}$ を確かめよ。
3. **円形直交アンサンブル**（`measure='coe'`）は、$W$ を $U(k)$ 上でハール分布に従うとして $S=W^TW$ で与えられる。
   $S$ がユニタリかつ対称であること、その分布がすべての
   $V\in U(k)$ に対して $S\mapsto V^TSV$ のもとで不変であること、そしてその固有位相の平均密度が一様であることを示せ。
   その固有位相が GOE と同様に $\langle r\rangle\approx0.53$ をもつことを数値的に確かめよ。

<details>
<summary>解答</summary>

1. $t_c=\operatorname{Tr}U_c$ と書く。$t_c$ は独立である（系 5.4）。$U_c$ に
   $e^{i\theta}$ を掛けてもその分布は変わらず、$t_c$ には $e^{i\theta}$ が掛かる。したがってモーメント
   $\mathbb E[t_at_b\bar t_c\bar t_d]$ は、各セクターが $(a,b)$ と
   $(c,d)$ に同じ回数だけ現れる場合、すなわち多重集合として $\{a,b\}=\{c,d\}$ である場合を除いて零であり、同様に $a\neq b$ に対し $\mathbb E[t_a\bar t_b]=0$ である。
   命題 5.7 により $\mathbb E[(U_c)_{ii}\overline{(U_c)_{jj}}]=\langle e_i,\mathbb E[U_ce_ie_j^*U_c^*]e_j\rangle=\delta_{ij}/k_c$ であるから、
   $\mathbb E|t_c|^2=1$、$\mathbb E|\operatorname{Tr}U|^2=C$ である。$\mathbb E|\sum_at_a|^4$ において、
   $a=b$ の項は $\sum_a\mathbb E|t_a|^4$ を与え、$a\neq b$ の項は
   $2\sum_{a\neq b}\mathbb E|t_a|^2\,\mathbb E|t_b|^2=2C(C-1)$ を与える（$(c,d)$ の選び方が二通り）。$k=1$ に対しては $|t|=1$ である。
   $k\ge2$ に対しては $\mathbb E|t|^4=\sum_{i,j,i',j'}\mathbb E[U_{ii}U_{jj}\overline{U_{i'i'}U_{j'j'}}]$ である。
   定理 5.10 において、$\sigma=\pi$ の二つの項は同一の仕方で $\{i',j'\}=\{i,j\}$ を強制し、和をとるとそれぞれ
   $k^2\mathrm{Wg}(e)$ を与える。$\sigma\neq\pi$ の二つの項は $i=j=i'=j'$ を強制し、それぞれ
   $k\,\mathrm{Wg}(\tau)$ を与える。合計は $2\big(\frac{k^2}{k^2-1}-\frac1{k^2-1}\big)=2$ である。$C=1$、
   $k=D$ の場合が $U(D)$ に対する $1$ と $2$ を与える。$M_1\oplus M_2\oplus M_3$ に対しては
   $\mathbb E|\operatorname{Tr}U|^4=1+2+2+12=17$ であり、$U(6)$ のハール・ユニタリに対する $2$ と対照的である。

```python
tr = U.trace                                                        # Tr U = sum_c Tr U_c, shape (B,)
mc(tr.abs() ** 2, 3.0, "E|Tr U|^2")
mc(tr.abs() ** 4, 17.0, "E|Tr U|^4")
```

2. すべての添字が $1$ のとき、デルタはすべて $1$ であるから
   $\mathbb E|U_{11}|^{2n}=\sum_{\sigma,\pi}\mathrm{Wg}(\sigma\pi^{-1})=\mathbf 1^TW\mathbf 1$ である。ここで
   $W=\Gamma^{-1}=(\mathrm{Wg}(\sigma^{-1}\pi))$、$\Gamma=(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)})$ である。ここでは
   $\mathrm{Wg}$ が巡回型のみに依存すること、および $\sigma\pi^{-1}=\sigma(\pi^{-1}\sigma)\sigma^{-1}$ が
   $\pi^{-1}\sigma=(\sigma^{-1}\pi)^{-1}$ と共役であり、後者は $\sigma^{-1}\pi$ と同じ巡回型をもつことを用いた。$\Gamma$ の
   各行和は $c_n(k)=\sum_{\rho\in S_n}k^{\#\mathrm{cycles}(\rho)}$ であるから、$\Gamma\mathbf 1=c_n(k)\mathbf 1$、
   $W\mathbf 1=\mathbf 1/c_n(k)$、$\mathbf 1^TW\mathbf 1=n!/c_n(k)$ である。最後に、帰納法により
   $c_n(x)=x(x+1)\cdots(x+n-1)$ である：$\{1,\dots,n\}$ の置換は
   $\{1,\dots,n-1\}$ の置換から、$n$ を不動点として加える（巡回が一つ増え、因子 $x$）か、$n$ を
   $n-1$ 個の元のいずれかの後ろにその巡回の中で挿入する（$n-1$ 通り、巡回の個数は同じ）ことにより得られる。これは全単射であるから、
   $c_n(x)=(x+n-1)c_{n-1}(x)$ である。

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

3. $W^T$ はユニタリであるから $S^*S=W^*\bar WW^TW=1$ であり、$S^T=W^TW=S$ である。$V\in U(k)$ に対し $WV$ は
   $W$ と同じ分布をもち（右不変性、補題 5.1）、$(WV)^T(WV)=V^TSV$ である。$V=e^{i\theta/2}1$ とすると $e^{i\theta}S$
   は $S$ と同じ分布をもつから、固有位相の集合の分布は回転で不変であり、平均密度も同様である。
   したがって平均密度は一様である。特に間隔はアンフォールディングを必要としない。COE の局所
   統計は $k$ が大きいとき GOE のそれと一致する（F. J. Dyson, J. Math. Phys. 3
   (1962), 140–156; Mehta, *op. cit.*；証明なしに述べる）。

```python
alg_s = TypeIAlgebra([32], [32], precision='double')
S = alg_s.random_unitary_operator(400, measure='coe').matrix[:, 0]
assert close(S, S.mT) and close(S @ S.mH, torch.eye(32, dtype=dt))
phases = torch.linalg.eigvals(S).angle().sort(dim=-1).values
s = torch.diff(torch.cat([phases, phases[:, :1] + 2 * math.pi], dim=-1), dim=-1)   # around the circle
r = torch.minimum(s, s.roll(-1, -1)) / torch.maximum(s, s.roll(-1, -1))
mc(r.mean(-1), 0.5307, "<r>, COE(32) eigenphases")
```

</details>

**次へ：**[第6章](06_quantum_channels.md) では量子チャネルを扱う：クラウス表現とスタインスプリング
表現、トレースに関する双対、そしてセクター間の写像。ここでハール平均として現れた条件付き
期待値 $E$ が最初の例である。
