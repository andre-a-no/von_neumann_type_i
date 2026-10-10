[English](../03_traces_and_states.md) · [Español](../es/03_traces_and_states.md) · [Français](../fr/03_traces_and_states.md) · [Deutsch](../de/03_traces_and_states.md) · [中文](../zh/03_traces_and_states.md) · **日本語** · [Русский](../ru/03_traces_and_states.md)

# 3. トレース、状態、エントロピー

$M_n(\mathbb C)$ 上のトレースは本質的に一つしかない。$\bigoplus_cM_{k_c}$ 上には
$C$ 個のパラメータをもつ族があり、その選択は重要である：密度、エントロピー、写像の双対はすべてそれに依存する。本章ではトレースを分類し、状態を密度で表し、基本的なエントロピー不等式を証明し、最後にもっともらしいが偽であるトレース不等式を取り上げる。

**本章で学ぶこと**：ライブラリの三つのトレースとその間の変換方法、状態の古典的部分、フォン・ノイマン・エントロピーのセクターへの分解、クラインの不等式とその帰結、そしてマイケルソン・コントラスト。

本章を通じて、$\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$（すべての $k_c\ge1$）、
$D=\sum_ck_c$ とし、$P_c$ は極小中心射影、$A=\bigoplus_cA_c$ とする。

## 3.1 トレース

$\mathcal M$ 上の **トレース** とは、すべての
$A,B$ に対し $\tau(AB)=\tau(BA)$ を満たす線形汎関数 $\tau$ のことである。すべての $A$ に対し $\tau(A^*A)\ge0$ のとき **正** であるといい、さらに
$\tau(A^*A)=0$ ならば $A=0$ となるとき **忠実** であるという。

**定理 3.1.** $\mathcal M$ 上の任意のトレースは、$w\in\mathbb C^C$ により
$\operatorname{tr}_w(A)=\sum_cw_c\operatorname{Tr}A_c$ の形をもつ。それが正であるための必要十分条件はすべての $w_c\ge0$ であり、忠実であるための必要十分条件はすべての $w_c>0$ である。

*証明.* $E^c_{ij}$ を第 $c$ ブロックの行列単位とする。トレース $\tau$ に対し
$\tau(E^c_{ij})=\tau(E^c_{i1}E^c_{1j})=\tau(E^c_{1j}E^c_{i1})=\delta_{ij}\tau(E^c_{11})$ である。
$w_c=\tau(E^c_{11})$ とおくと、線形性により $\tau(A)=\sum_c\sum_i(A_c)_{ii}w_c$ である。逆に
任意の $\operatorname{tr}_w$ はトレースである。正値性：$E^c_{11}=(E^c_{11})^*E^c_{11}$ であるから
$w_c\ge0$ は必要であり、$\operatorname{tr}_w(A^*A)=\sum_cw_c\|A_c\|_2^2$ により十分性と忠実性の判定条件が示される。$\square$

ライブラリは三つの忠実なトレースを実装している。

| 名称 | 重み $w_c$ | $1$ での値 | 意味 |
|---|---|---|---|
| `Tr_blunt` | $1$ | $D$ | $\bigoplus_c\mathbb C^{k_c}$ のトレース（重複度 $m_c=1$） |
| `Tr_norm` | $1/k_c$ | $C$ | 各ブロックを正規化したもの |
| `tau_vN` | $1/(Ck_c)$ | $1$ | セクターの重みが等しく $1/C$ であるトレース状態 |

第2章により、$\mathcal M$ が重複度 $m_c$ で作用する空間のトレースは
$w_c=m_c$ とした $\operatorname{tr}_w$ である。したがって忠実なトレースとは「各セクターをどれだけ数えるか」の選択である。

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, DensityMatrix, states, channels

torch.manual_seed(0)
alg = TypeIAlgebra([2, 3, 5], [2, 3, 5], precision='double')    # M_2 (+) M_3 (+) M_5, D = 10, C = 3
U = alg.random_unitary_operator(4)
A = alg.operator_from_eigenvalues(lambda d: torch.randn(4, d), batch_size=4, force_self_adjoint=True)
B = U @ A                                                       # a batch of four non-normal operators
one = alg.identity()

for f, value in ((lambda X: X.Tr_blunt(), 10), (lambda X: X.Tr_norm(), 3), (lambda X: X.tau_vN(), 1)):
    assert abs(f(one).item() - value) < 1e-12                   # normalisation
    assert (f(A @ B) - f(B @ A)).abs().max() < 1e-10            # trace property
assert (A.Tr_norm() - 3 * A.tau_vN()).abs().max() < 1e-12       # Tr_norm = C tau_vN
```

## 3.2 状態と密度

**状態** とは、$\mathcal M$ 上の線形汎関数 $\omega$ であって、正
（$\omega(A^*A)\ge0$）かつ正規化されている（$\omega(1)=1$）もののことである。

**定理 3.2.** $\operatorname{tr}_w$ を忠実なトレースとする。任意の線形汎関数
$\omega$ に対し、すべての $A$ について $\omega(A)=\operatorname{tr}_w(\rho_wA)$ を満たす $\rho_w\in\mathcal M$ がただ一つ存在する。汎関数 $\omega$ が正であるための必要十分条件は $\rho_w\ge0$ であり、状態であるための必要十分条件はさらに
$\operatorname{tr}_w(\rho_w)=1$ となることである。二つの忠実なトレースに対し、密度はブロックごとに次の関係にある。
$$
w_c\,(\rho_w)_c=w'_c\,(\rho_{w'})_c .
$$

*証明.* $\langle X,Y\rangle_w=\operatorname{tr}_w(X^*Y)$ は $\mathcal M$ 上の内積である
（定理 3.1）。有限次元空間
$\mathcal M$ におけるリースの表現定理により、$\omega(A)=\langle R,A\rangle_w=\operatorname{tr}_w(R^*A)$ を満たす $R$ がただ一つ存在する。
$\rho_w=R^*$ とおく。$\rho_w\ge0$ ならば
$\omega(A^*A)=\sum_cw_c\operatorname{Tr}(A_c\,(\rho_w)_c\,A_c^*)\ge0$ である。実際
$A_c(\rho_w)_cA_c^*\ge0$ である（系 1.7 を $C=(\rho_w)_c^{1/2}A_c^*$ として適用）。逆に
$\omega$ が正であるとし、$x\in\mathbb C^{k_c}$ をとり、$A$ をブロック $c$ における階数1の作用素
$|x\rangle\langle x|$ とする。すると $A=A^*A/\|x\|^2$ であり
$0\le\omega(A)=w_c\langle x,(\rho_w)_cx\rangle$ であるから $(\rho_w)_c\ge0$ である。正規化は
$\omega(1)=\operatorname{tr}_w(\rho_w)$ である。最後に、$\sum_cw_c\operatorname{Tr}((\rho_w)_cA_c)$ と
$\sum_cw'_c\operatorname{Tr}((\rho_{w'})_cA_c)$ はともにすべての $A$ に対して $\omega(A)$ に等しいから、一意性により
$w_c(\rho_w)_c=w'_c(\rho_{w'})_c$ である。$\square$

ライブラリは状態を
`Tr_blunt` に関する密度 $\rho=\rho_{\mathrm{blunt}}$ によって格納する（クラス `DensityMatrix`）。`density(trace)` と `DensityMatrix.from_density`
が変換を行う。状態の値は規約に依存しない。

```python
rho = states.random_density_matrix(alg, batch_size=4)          # Hilbert-Schmidt random states
for trace, tr in (('blunt', lambda X: X.Tr_blunt()), ('norm', lambda X: X.Tr_norm()),
                  ('tau_vN', lambda X: X.tau_vN())):
    d = alg.operator(rho.density(trace))
    assert (tr(d @ A) - (rho @ A).Tr_blunt()).abs().max() < 1e-10       # omega(A) for all traces
    assert (tr(d) - 1).abs().max() < 1e-10                               # tr_w(rho_w) = 1
    back = DensityMatrix.from_density(alg, rho.density(trace), trace=trace)
    assert (back.matrix - rho.matrix).abs().max() < 1e-12
```

> **物理的意味.** $\omega(A)$ はオブザーバブル $A$ の期待値である。
> blunt トレースに関する密度 $\rho$ は、$\bigoplus_c\mathbb C^{k_c}$ 上の量子力学における通常の密度行列である。他のトレースに関する密度は、
> 同じ環が重複度をもって作用するとき（第2章）、たとえばスピン系の
> 対称部分と非対称部分の上で現れる。

## 3.3 状態の古典的部分

$\omega$ の中心への制限はセクター上の確率分布
$p_c=\omega(P_c)=\operatorname{Tr}\rho_c$ である。$p_c>0$ ならば、セクター $c$ における **条件付き状態** は
$\mathbb C^{k_c}$ 上の密度 $\sigma_c=\rho_c/p_c$ をもち、
$$
\rho=\bigoplus_cp_c\,\sigma_c
$$
である（$p_c=0$ ならば、$\rho_c\ge0$ かつ $\operatorname{Tr}\rho_c=0$ であるから $\rho_c=0$ であり、その項は
省く）。
$\mathcal M$ 上の状態は、各セクターにおける量子状態の古典的な混合である。「どのセクターか」を測定してもそれは乱されない：測定チャネル
$\rho\mapsto\sum_cP_c\rho P_c$ は $\mathcal M$ 上で恒等写像である。なぜなら
$\mathcal M$ の任意の元はブロック対角だからである。これが超選択則の操作的な内容である。

```python
p = rho.sector_probabilities()                                   # shape (batch, C)
assert torch.allclose(p.sum(-1), torch.ones(4, dtype=p.dtype))
P = [alg.central([1.0 if c == d else 0.0 for d in range(alg.C)]) for c in range(alg.C)]
measure_sector = channels.measurement_channel(P)
assert (measure_sector(rho).matrix - rho.matrix).abs().max() < 1e-12   # no disturbance
```

## 3.4 エントロピー

（blunt）密度 $\rho$ をもつ状態の **フォン・ノイマン・エントロピー** は
$S(\rho)=-\operatorname{Tr}\rho\log\rho=-\sum_i\lambda_i\log\lambda_i$ である。ここで和は
$\rho$ の固有値にわたり、$0\log0=0$ とする。

**定理 3.3（セクターの混合のエントロピー）.** $p_c$ と $\sigma_c$ を上のとおりとする
（$p_c=0$ のセクターは省く）と、
$$
S(\rho)=H(p)+\sum_cp_c\,S(\sigma_c),\qquad H(p)=-\sum_cp_c\log p_c .
$$

*証明.* $\rho$ の非零固有値は、$\sigma_c$ の非零
固有値 $\mu$（$p_c>0$ のセクター）に対する数 $p_c\mu$ である。零固有値は両辺に $0$ を寄与する。したがって
$-\sum_c\sum_{\mu>0}p_c\mu\log(p_c\mu)=-\sum_cp_c\log p_c\sum_{\mu>0}\mu-\sum_cp_c\sum_{\mu>0}\mu\log\mu$
であり、$\sum_{\mu>0}\mu=\operatorname{Tr}\sigma_c=1$ である。$\square$

他のトレースに関するエントロピーは $S_w(\rho)=-\operatorname{tr}_w(\rho_w\log\rho_w)$ である。
定理 3.2 により $(\rho_w)_c=\rho_c/w_c$ であるから、
$$
S_w(\rho)=-\sum_cw_c\operatorname{Tr}\Big(\frac{\rho_c}{w_c}\log\frac{\rho_c}{w_c}\Big)
=S(\rho)+\sum_cp_c\log w_c .
$$
エントロピーは状態 *と* トレースの性質である：同じ状態でも
`Tr_blunt` と `tau_vN` に関しては異なるエントロピーをもつ。

```python
S = rho.von_neumann_entropy()
sigma = [rho.matrix[:, c, :k, :k] / p[:, c, None, None] for c, k in enumerate(alg.k_factors)]
S_sec = [-(w * torch.log(w)).sum(-1) for w in (torch.linalg.eigvalsh(s) for s in sigma)]
H = -(p * torch.log(p)).sum(-1)
assert torch.allclose(S, H + sum(p[:, c] * S_sec[c] for c in range(alg.C)))    # Theorem 3.3

w_tau = torch.tensor([1 / (alg.C * k) for k in alg.k_factors], dtype=torch.float64)
d_tau = alg.operator(rho.density('tau_vN'), is_positive=True)
S_tau = -(d_tau @ d_tau.log()).tau_vN().real
assert torch.allclose(S_tau, S + (p * torch.log(w_tau)).sum(-1))
```

## 3.5 相対エントロピーとクラインの不等式

密度 $\rho,\sigma$（$\sigma>0$）をもつ状態に対し、**相対エントロピー** を
$S(\rho\|\sigma)=\operatorname{Tr}\rho\log\rho-\operatorname{Tr}\rho\log\sigma$ と定める。ここで
$\operatorname{Tr}\rho\log\rho:=\sum_i\lambda_i\log\lambda_i$ は $\rho$ の固有値にわたる和であり、
$0\log0=0$ とする（各ブロックについても同様）。これはトレースに依存しない：
因子 $1/w_c$ が相殺するので $\log(\rho_w)_c-\log(\sigma_w)_c=\log\rho_c-\log\sigma_c$ であり、
また $\operatorname{tr}_w(\rho_wX)=\operatorname{Tr}(\rho X)$ である。

**補題（接線不等式）.** $f:I\to\mathbb R$ を区間 $I$ 上の凸関数とし、
$b\in I$ で微分可能とする。このときすべての $a\in I$ に対し $f(a)\ge f(b)+(a-b)f'(b)$ であり、$f$ が狭義凸ならば
$a\neq b$ に対して不等号は狭義である。

*証明.* $t\in(0,1]$ に対し、凸性から $f(b+t(a-b))\le(1-t)f(b)+tf(a)$、すなわち
$\big(f(b+t(a-b))-f(b)\big)/t\le f(a)-f(b)$ を得る。$t\to0$ とすれば不等式を得る。$f$
が狭義凸で $a\neq b$ のとき、$m=(a+b)/2$ とおくと
$f(a)-f(b)>2\big(f(m)-f(b)\big)\ge2(m-b)f'(b)=(a-b)f'(b)$ である。第一段は $m$ における狭義
凸性により、第二段は今証明した不等式を $m$ に適用したものである。$\square$

**定理 3.4（クラインの不等式）.** $f:I\to\mathbb R$ を区間 $I$ 上の凸関数とし、
$X,Y$ をスペクトルが $I$ に含まれる自己共役作用素とし、$f$ は $Y$ のすべての固有値で微分可能とする。
このとき $\operatorname{Tr}[f(X)-f(Y)-(X-Y)f'(Y)]\ge0$ である。$f$ が狭義凸ならば、等号は
$X=Y$ のときに限り成り立つ。

*証明.* 正規直交固有基底により $X=\sum_ia_i|u_i\rangle\langle u_i|$、$Y=\sum_jb_j|v_j\rangle\langle v_j|$ と書き、
$c_{ij}=|\langle u_i,v_j\rangle|^2$ とおく。すると $\sum_ic_{ij}=\sum_jc_{ij}=1$ である。
このとき $\operatorname{Tr}f(X)=\sum_{ij}c_{ij}f(a_i)$、$\operatorname{Tr}f(Y)=\sum_{ij}c_{ij}f(b_j)$、
$\operatorname{Tr}Xf'(Y)=\sum_{ij}c_{ij}a_if'(b_j)$、$\operatorname{Tr}Yf'(Y)=\sum_{ij}c_{ij}b_jf'(b_j)$ である。
したがって左辺は $\sum_{ij}c_{ij}[f(a_i)-f(b_j)-(a_i-b_j)f'(b_j)]$ に等しく、これは接線不等式により
非負の項の和である。$f$ が狭義凸で和が零ならば、
$c_{ij}>0$ のときは常に $a_i=b_j$ であるから、すべての $j$ に対し $Xv_j=\sum_ia_i\langle u_i,v_j\rangle u_i=b_jv_j$ となり、
$X=Y$ である。$\square$

**系 3.5.** $S(\rho\|\sigma)\ge0$ であり、等号は $\rho=\sigma$ のときかつそのときに限り成り立つ。特に
$S(\rho)\le\log D$ であり、等号は $\rho=1/D$ のときかつそのときに限り成り立つ。

*証明.* 定理 3.4 を $f(t)=t\log t$, $f(0)=0$ に適用する。これは
$I=[0,\infty)$ 上で狭義凸であり（$(0,\infty)$ 上で $f''(t)=1/t>0$ であり、$f$ は $0$ で連続）、
$(0,\infty)\supseteq\operatorname{spec}\sigma$ 上で微分可能で $f'(t)=\log t+1$ である。
$X=\rho$, $Y=\sigma$ とすると、左辺は
$\operatorname{Tr}\rho\log\rho-\operatorname{Tr}\sigma\log\sigma-\operatorname{Tr}(\rho-\sigma)(\log\sigma+1)
=S(\rho\|\sigma)$ である。$\operatorname{Tr}\rho=\operatorname{Tr}\sigma=1$ だからである。第二の主張については
$\sigma=1/D$ ととる：$0\le S(\rho\|1/D)=-S(\rho)+\log D$。$\square$

**定理 3.6（連鎖律；中心への制限に関する単調性）.** $\rho,\sigma$ を
$\sigma>0$ なる状態とし、$p_c=\operatorname{Tr}\rho_c$、$q_c=\operatorname{Tr}\sigma_c>0$、$\hat\sigma_c=\sigma_c/q_c$ とおき、
$p_c>0$ に対し $\hat\rho_c=\rho_c/p_c$ とおく。このとき、和を $p_c>0$ なる $c$ にわたってとると（その他の
セクターは $0\log0=0$ により $0$ を寄与する）、
$$
S(\rho\|\sigma)=\sum_cp_c\log\frac{p_c}{q_c}+\sum_cp_c\,S(\hat\rho_c\|\hat\sigma_c)
\ \ge\ \sum_cp_c\log\frac{p_c}{q_c}=S(E\rho\|E\sigma),
$$
が成り立つ。ここで $E$ は中心の上への条件付き期待値である（命題 2.7）。

*証明.* $p_c>0$ に対し $\operatorname{Tr}\rho_c\log\rho_c=p_c\log p_c+p_c\operatorname{Tr}\hat\rho_c\log\hat\rho_c$
（固有値は $p_c\mu$、定理 3.3 と同様）かつ $\log\sigma_c=(\log q_c)1+\log\hat\sigma_c$ であるから、
$\operatorname{Tr}\rho_c\log\sigma_c=p_c\log q_c+p_c\operatorname{Tr}\hat\rho_c\log\hat\sigma_c$ であり、
$\operatorname{Tr}\rho_c\log\rho_c-\operatorname{Tr}\rho_c\log\sigma_c=p_c\log(p_c/q_c)+p_cS(\hat\rho_c\|\hat\sigma_c)$ である。
$p_c=0$ ならば $\rho_c=0$ であり、そのセクターは $0$ を寄与する。
$c$ について和をとる。第二の和は各ブロックにおける系 3.5 により $\ge0$ である。最後に
$E\rho=\bigoplus_c(p_c/k_c)1_{k_c}$、$E\sigma=\bigoplus_c(q_c/k_c)1_{k_c}$ であり、
$\operatorname{Tr}E\rho(\log E\rho-\log E\sigma)=\sum_ck_c\frac{p_c}{k_c}\log\frac{p_c}{q_c}$ である。$\square$

したがって、セクター内部の量子情報を忘れることは、二つの状態の識別可能性を減らすことしかできない。残るのはセクター分布の古典的なカルバック・ライブラー情報量である。

```python
sigma_state = states.random_density_matrix(alg, batch_size=4)
q = sigma_state.sector_probabilities()
D_rel = states.relative_entropy(rho, sigma_state)
kl = (p * torch.log(p / q)).sum(-1)
E = channels.center_expectation(alg)
assert (D_rel >= kl - 1e-12).all()                                       # Theorem 3.6
assert torch.allclose(states.relative_entropy(E(rho), E(sigma_state)), kl)
assert (rho.von_neumann_entropy() <= math.log(10) + 1e-12).all()         # Corollary 3.5
assert abs(states.maximally_mixed_state(alg).von_neumann_entropy().item() - math.log(10)) < 1e-12
```

## 3.6 マイケルソン・コントラスト

$A\ge0$, $A\neq0$ に対し、**マイケルソン・コントラスト** を
$\Delta(A)=(\lambda_{\max}-\lambda_{\min})/(\lambda_{\max}+\lambda_{\min})\in[0,1]$ と定める。ここで
両端の固有値はすべてのセクターにわたってとる（また $\Delta(0)=0$ とする）。これは $A$ が
単位元のスカラー倍からどれだけ離れているかを測る：$\Delta(A)=0$ $\iff$ $A=\lambda1$。

**命題 3.7.** $\Phi:\mathcal M\to\mathcal M$ が線形、正かつ単位的
（$\Phi(1)=1$）ならば、すべての $A\ge0$ に対し $\Delta(\Phi(A))\le\Delta(A)$ である。

*証明.* $\lambda_{\min}1\le A\le\lambda_{\max}1$ であり、正値性と単位性から
$\lambda_{\min}1\le\Phi(A)\le\lambda_{\max}1$ を得る。したがって
$\lambda_{\min}(\Phi(A))\ge\lambda_{\min}$ かつ $\lambda_{\max}(\Phi(A))\le\lambda_{\max}$ である。
$\Phi(A)=0$ ならば主張は自明である。そうでなければ、$r=\lambda_{\min}/\lambda_{\max}$ とおくと
$\Delta=(1-r)/(1+r)$ は $r$ について減少であり、$r(\Phi(A))\ge r(A)$ である。$\square$

```python
A_pos = alg.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(4, d), batch_size=4, force_positive=True)
Phi = channels.random_mixed_unitary_channel(alg, 3, batch_size=4)       # unital
assert (Phi(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
assert (E(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
```

## 3.7 もっともらしいが偽である不等式

$M_k$ の正の $X,Y$ とユニタリ $U$ に対し、$|\operatorname{Tr}(XUY)|\le\operatorname{Tr}(XY)$ は成り立つだろうか。これは
$U=1$ のときは成り立ち、コーシー・シュワルツ型の主張のように見える。しかし $X$
と $Y$ が可換でないときは常に偽である。

**命題 3.8.** $\max_U|\operatorname{Tr}(XUY)|=\|YX\|_1\ge\operatorname{Tr}(XY)$ であり、等号は
$XY=YX$ のときかつそのときに限り成り立つ。

*証明.* $\operatorname{Tr}(XUY)=\operatorname{Tr}(YX\,U)$ であるから、定理 1.9 により最大値は $\|YX\|_1$ であり、
$\|YX\|_1\ge|\operatorname{Tr}(YX)|=\operatorname{Tr}(XY)$ である。ここで
$X^{1/2}YX^{1/2}=(Y^{1/2}X^{1/2})^*(Y^{1/2}X^{1/2})\ge0$（系 1.7）より
$\operatorname{Tr}(XY)=\operatorname{Tr}(X^{1/2}YX^{1/2})\ge0$ であることを用いた。$XY=YX$ ならば
$X^{1/2}$ は $Y$ と可換であるから（定理 1.5(d)）、$YX=X^{1/2}YX^{1/2}\ge0$ であり、
$\|YX\|_1=\operatorname{Tr}(YX)$ である（$P\ge0$ に対し $\|P\|_1=\operatorname{Tr}P$、1.4節）。逆に
$\|YX\|_1=\operatorname{Tr}(YX)$ と仮定する。$|YX|=0$ ならば $YX=0=(YX)^*=XY$ である。そうでなければ
$YX=W|YX|$ と書く。定理 1.9 により、$\operatorname{Tr}(YX\cdot1)=\|YX\|_1$ は $U=1$ に対する
コーシー・シュワルツの段階で等号が成り立つことを意味する：
$|YX|^{1/2}$ と $W|YX|^{1/2}$ は非負の因子で比例し、その因子は
ノルムを比較すると $1$ であるから、$W|YX|^{1/2}=|YX|^{1/2}$ かつ $YX=W|YX|=|YX|\ge0$ である。正
行列は自己共役であるから、$YX=(YX)^*=XY$ である。$\square$

ライブラリは極分解によって最大値を与える $U$ を求める。解説論文は、$U$ のモンテカルロ・サンプリングではそれに届かず、$k$ が大きくなるほどその差が広がることを数値的に示している
——これは論文の中心的な主題の一つである。以下のコードは一つの組についてこれを示す。

```python
Mk = TypeIAlgebra([6], [6], precision='double')
X = Mk.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
Y = Mk.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
YX = Y @ X
W = YX @ YX.abs().inverse()                                          # polar part of YX
gap = (X @ W.adjoint() @ Y).trace.abs() - (X @ Y).trace.real
assert abs(gap.item() - (YX.trace_norm() - (X @ Y).trace.real).item()) < 1e-10
assert gap.item() > 0                                                # the inequality fails
samples = Mk.random_unitary_operator(20_000)
best = (X @ samples @ Y).trace.abs().max() - (X @ Y).trace.real
print(f"exact sup {gap.item():.4f}, best of 20000 Haar samples {best.item():.4f}")
assert best.item() < gap.item()
```

## 3.8 比較の対象

QuTiP（`entropy_vn`, `entropy_relative`）や量子情報のパッケージ（例えば
`qutip`, `toqito`）は、単一のヒルベルト空間上の密度行列のエントロピーと相対エントロピーを計算する。$\bigoplus_c\mathbb C^{k_c}$ 上ではそれらは本章の blunt トレースの規約と一致する。それらが提供しないのは定理 3.1 のトレースの族と定理 3.2 の変換であり、これらは環が重複度をもって作用するとき、あるいはセクターに異なる重みをつけなければならないときに直ちに必要となる。

## 演習問題

1. トレース状態 $\tau_{\mathrm{vN}}$（セクターの重みが等しい）の
   $\tau_{\mathrm{vN}}$ に関する密度が $1$ であることを示し、その blunt 密度とエントロピー $S$ を計算せよ。
2. $\operatorname{tr}_w$ を忠実とする。$S_w(\rho)\le\log\operatorname{tr}_w(1)$ であり、等号が
   $\rho_w=1/\operatorname{tr}_w(1)$ のときかつそのときに限り成り立つことを示せ。
3. $M_D$ 上の（$\mathcal M$ 上ではない！）測定チャネル $\rho\mapsto\sum_cP_c\rho P_c$ は状態を乱すことを示し、二つのセクターの重ね合わせについて測定前後の $S(\rho)$ を計算せよ。

<details>
<summary>解答</summary>

1. $\tau_{\mathrm{vN}}(A)=\tau_{\mathrm{vN}}(1\cdot A)$ であるから、一意性により $\rho_{\tau}=1$ である。
   定理 3.2 により blunt 密度は $\bigoplus_c1_{k_c}/(Ck_c)$ であり、$p_c=1/C$、
   $\sigma_c=1/k_c$ であるから、定理 3.3 により $S=\log C+\frac1C\sum_c\log k_c$ である。

```python
tau = states.tracial_state(alg)
assert (tau.density('tau_vN') - one.matrix).abs().max() < 1e-12
assert abs(tau.von_neumann_entropy().item() - (math.log(3) + sum(math.log(k) for k in (2, 3, 5)) / 3)) < 1e-12
```

2. 和は $p_c>0$ なる $c$ にわたる。このとき
   $S_w(\rho)=S(\rho)+\sum_cp_c\log w_c=-\sum_cp_c\log p_c+\sum_cp_cS(\sigma_c)+\sum_cp_c\log w_c
   \le\sum_cp_c\log\frac{k_cw_c}{p_c}\le\log\sum_{p_c>0}k_cw_c\le\log\operatorname{tr}_w(1)$ である。ここで
   $S(\sigma_c)\le\log k_c$（系 3.5）と、凹関数 $\log$ に対する有限形のイェンセンの不等式
   $\sum_cp_c\log x_c\le\log\sum_cp_cx_c$（これは凹性の定義から項数に関する帰納法で従う）を
   $x_c=k_cw_c/p_c$ として用いた。$\log$ は狭義凹であるから、この不等式はすべての $x_c$ が等しい場合を除いて狭義である。すべてで等号が成り立つためには、すべての $p_c>0$、$\sigma_c=1/k_c$ かつ $k_cw_c/p_c$ が一定、すなわち
   $p_c=k_cw_c/\operatorname{tr}_w(1)$、すなわち $(\rho_w)_c=\rho_c/w_c=1_{k_c}/\operatorname{tr}_w(1)$ でなければならない。
   逆に $\rho_w=1/\operatorname{tr}_w(1)$ に対しては
   $S_w(\rho)=-\operatorname{tr}_w(1)\cdot\frac1{\operatorname{tr}_w(1)}\log\frac1{\operatorname{tr}_w(1)}=\log\operatorname{tr}_w(1)$ である。

3. $M_D$ において、二つの
   セクターの単位ベクトルの純粋な重ね合わせ $\psi=(\psi_1\oplus\psi_2)/\sqrt2$ は $S=0$ をもつ。測定後の状態は
   $\frac12|\psi_1\rangle\langle\psi_1|\oplus\frac12|\psi_2\rangle\langle\psi_2|$ であり $S=\log2$ である。しかし
   $\mathcal M$ に制限すると、$\psi$ のベクトル状態は混合と一致する：
   すべての $A\in\mathcal M$ に対し
   $\langle\psi,A\psi\rangle=\tfrac12\langle\psi_1,A_1\psi_1\rangle+\tfrac12\langle\psi_2,A_2\psi_2\rangle$
   である。$\mathcal M$ のいかなるオブザーバブルもコヒーレンスを検出しない。これが超選択則の内容である。

</details>

**次へ：**第4章（準備中）では対称性と表現を扱う：物理系への群作用がどのようにしてセクターを生み出すか、そしてライブラリがスピン系や量子ビットレジスタに対応する環をどのように構成するか。
