[English](../02_von_neumann_algebras.md) · [Español](../es/02_von_neumann_algebras.md) · [Français](../fr/02_von_neumann_algebras.md) · [Deutsch](../de/02_von_neumann_algebras.md) · [中文](../zh/02_von_neumann_algebras.md) · **日本語** · [Русский](../ru/02_von_neumann_algebras.md)

# 2. 有限次元フォン・ノイマン環

フォン・ノイマン環とは、「可換子環を二度とる」操作で閉じた作用素環である。有限次元ではこれは純粋に代数的な条件であり、本章の主定理は、そのような環はすべて、ユニタリな基底変換を除いて、全行列環の直和（各成分はある重複度で繰り返される）であることを述べる。これがライブラリの土台となる構造である。

**本章で学ぶこと**：二重可換子環定理、完全な証明付きの構造定理、対称性の群がどのようにしてそのような環を生み出すか、中心とは何か、そしてライブラリがこれらすべてをどのように表現するか。

本章を通じて、*$\mathbb C^n$ 上の $*$-環* とは、積と随伴について閉じていて単位元 $1$ を含む複素線形部分空間
$\mathcal M\subseteq M_n(\mathbb C)$ を意味する。

## 2.1 可換子環と二重可換子環定理

集合 $S\subseteq M_n(\mathbb C)$ に対し、その **可換子環** を
$S'=\{X\in M_n(\mathbb C):XA=AX\ \forall A\in S\}$ と定める。これは常に単位的な環である。
$S$ が随伴について閉じている（$S^*=S$）ならば $S'$ は $*$-環である。実際、
すべての $A\in S$ に対する $XA=AX$ から、すべての $A\in S$ に対して $X^*A^*=A^*X^*$ が従う。$S$ の任意の元は $S'$ の任意の元と可換であるから $S\subseteq S''$ である。これを $S'$ に適用すると
$S'\subseteq S'''$ を得て、$S\subseteq S''$ の可換子環をとる（包含関係は逆転する）と
$S'''\subseteq S'$ を得る。したがって $S'=S'''$ である。

環が $\mathbb C^n$ の部分空間 $E$ に作用するときは、正規直交基底によって $E$ を
$\mathbb C^{\dim E}$ と同一視し、$\mathbb C^n$ 上で述べたすべての結果がそれに適用できるようにする。

**定義.** $\mathbb C^n$ 上の **フォン・ノイマン環** とは、$\mathcal M=\mathcal M''$ を満たす $*$-環
$\mathcal M$ のことである。その **中心** は
$Z(\mathcal M)=\mathcal M\cap\mathcal M'$ である。**因子環**（factor）とは、
$Z(\mathcal M)=\mathbb C1$ を満たすフォン・ノイマン環のことである。

無限次元では $\mathcal M=\mathcal M''$ は位相的な条件（弱作用素位相での閉性）である。有限次元ではこれは自動的に成り立つ。

**補題 2.1.** $\mathcal M$ を $\mathbb C^n$ 上の $*$-環とし、$K\subseteq\mathbb C^n$ を
$\mathcal MK\subseteq K$ を満たす部分空間とする。このとき $K$ への直交射影 $P$ は
$\mathcal M'$ に属する。

*証明.* $A\in\mathcal M$ に対し、$AK\subseteq K$ かつ $A^*K\subseteq K$ である（$A^*\in\mathcal M$ だから）。
第二の包含から $AK^\perp\subseteq K^\perp$ が従う：$y\in K^\perp$, $x\in K$ に対し
$\langle x,Ay\rangle=\langle A^*x,y\rangle=0$ である。したがって $A$ は $K$ を $K$ に、$K^\perp$
を $K^\perp$ に写し、これは $AP=PAP=PA$ を意味する。$\square$

**定理 2.2（二重可換子環定理）.** $\mathbb C^n$ 上の任意の $*$-環 $\mathcal M$ は
$\mathcal M=\mathcal M''$ を満たす。

*証明.* $\mathcal M\subseteq\mathcal M''$ である。$T\in\mathcal M''$ とする。

*ステップ1：一つのベクトル.* $x\in\mathbb C^n$ に対し、部分空間 $K=\mathcal Mx$ は
$\mathcal M$ で不変であり、$x=1x$ を含む。補題 2.1 によりその射影 $P$ は $\mathcal M'$ に属するから、
$TP=PT$ かつ $Tx=TPx=PTx\in K$ である：$Tx=Ax$ を満たす $A\in\mathcal M$ が存在する。

*ステップ2：すべてのベクトルを同時に.* $\mathbb C^n\otimes\mathbb C^n=(\mathbb C^n)^{\oplus n}$ 上で
$*$-環 $\tilde{\mathcal M}=\{A\oplus\dots\oplus A:A\in\mathcal M\}$ を考える。ブロック行列
$[X_{ij}]_{i,j=1}^n$ がすべての $A\oplus\dots\oplus A$ と可換であるための必要十分条件は、すべての $i,j$ に対して $X_{ij}A=AX_{ij}$ となることであるから、$\tilde{\mathcal M}'=\{[X_{ij}]:X_{ij}\in\mathcal M'\}$ である。作用素
$\tilde T=T\oplus\dots\oplus T$ はそのような任意の $[X_{ij}]$ と可換である。なぜなら
$TX_{ij}=X_{ij}T$ だからである。よって $\tilde T\in\tilde{\mathcal M}''$ である。ステップ1を $\tilde{\mathcal M}$、
$\tilde T$ およびベクトル $\xi=(e_1,\dots,e_n)$ に適用すると、
$(Te_1,\dots,Te_n)=(Ae_1,\dots,Ae_n)$、すなわち $T=A\in\mathcal M$ となる $A\in\mathcal M$ が存在する。$\square$

したがって有限次元では「フォン・ノイマン環」と「行列からなる単位的 $*$-環」は同じものである。それでも二重可換子環定理は正しい観点を与える：物理的には、オブザーバブルの環はきわめて多くの場合、可換子環として *与えられる* ——すなわち対称性と可換な作用素全体として（2.4節）。

## 2.2 中心

**補題 2.3.** $\mathcal Z$ を $\mathbb C^n$ 上の可換な $*$-環とする。零でなく互いに直交し
$\sum_cP_c=1$ を満たす射影 $P_1,\dots,P_C\in\mathcal Z$ であって、$\mathcal Z=\operatorname{span}\{P_1,\dots,P_C\}$ となるものが存在する。それらは $\mathcal Z$ の極小な非零射影であり、順序を除いて一意である。

*証明.* 任意の $z\in\mathcal Z$ は正規である（$z^*\in\mathcal Z$ は $z$ と可換）。
$z_1,\dots,z_r$ を $\mathcal Z$ を張る元とする。$\mathbb C^n$ を同時固有空間に分解する。
$z_1$ の固有空間は互いに直交し $\mathbb C^n$ を張る（定理 1.3）。それらは
$z_2,\dots,z_r$ およびその随伴（これらは $z_1$ と可換）で不変であるから、各
$z_j$ の $z_1$ の固有空間への制限は再び正規であり、$z_1$ の各固有空間を $z_2$ の固有空間に分解でき、以下同様である。その結果、直交分解
$\mathbb C^n=E_1\oplus\dots\oplus E_C$ が得られる。ここで各 $E_c$ は零でない部分空間であり、その上で任意の
$z\in\mathcal Z$ はスカラー $\chi_c(z)$ として作用する。ただし同じ関数
$\chi_c:\mathcal Z\to\mathbb C$ をもつ部分空間は併合しておき、$c\neq d$ に対し $\chi_c\neq\chi_d$ となるようにする。
$P_c$ を $E_c$ への射影とする。すると任意の $z\in\mathcal Z$ に対し $z=\sum_c\chi_c(z)P_c$ である。

残るは $P_c\in\mathcal Z$ を示すことである。$d\neq c$ に対し
$\chi_c(z_{cd})\neq\chi_d(z_{cd})$ となる $z_{cd}\in\mathcal Z$ を選ぶ。積
$\prod_{d\neq c}\big(z_{cd}-\chi_d(z_{cd})1\big)/\big(\chi_c(z_{cd})-\chi_d(z_{cd})\big)$
は $\mathcal Z$ に属し、$E_c$ 上では $1$ として、すべての $E_d$, $d\neq c$ 上では $0$ として作用する。すなわちこれは
$P_c$ である（$C=1$ のときは $P_1=1$ であり、空積は $1$ である）。極小性と一意性：射影 $\sum_c a_cP_c$ は $a_c\in\{0,1\}$ を満たすから、
$\mathcal Z$ の非零射影は $P_c$ たちの空でない集合の和である。$\square$

これを $\mathcal Z=Z(\mathcal M)$ に適用したとき、射影 $P_c$ を $\mathcal M$ の **極小中心射影** という。任意の $A\in\mathcal M$ は
$A=\sum_cAP_c$（$AP_c=P_cAP_c$）と分解されるから、
$\mathcal M=\bigoplus_c\mathcal MP_c$ である（$E_c=P_c\mathbb C^n$ 上に作用する環の直交直和）。

**補題 2.4.** 各 $\mathcal M_c=\mathcal MP_c$ は、単位元 $P_c$ をもつ $E_c$ 上の環として、因子環である。

*証明.* $y\in\mathcal M_c$ が $\mathcal M_c$ と可換であるとする。$A\in\mathcal M$ に対し、$P_c$ が中心的であることと $y=yP_c=P_cy$ を用いて
$Ay=AP_cy=(AP_c)y=y(AP_c)=yP_cA=yA$ である。
よって $y\in Z(\mathcal M)$、$y=\sum_d\chi_d(y)P_d$ であり、$y=yP_c=\chi_c(y)P_c$ である。$\square$

## 2.3 構造定理

**補題 2.5（因子環）.** $\mathcal N$ を空間 $E$ 上の因子環とする。整数
$k,m\ge1$ とユニタリ $V:\mathbb C^k\otimes\mathbb C^m\to E$ であって
$V^*\mathcal NV=M_k(\mathbb C)\otimes1_m$ となるものが存在する。さらに $V^*\mathcal N'V=1_k\otimes M_m(\mathbb C)$ である。

*証明.* *ステップ1：極小射影.* $\mathcal N$ は零でない射影（例えば $1$）を含む。
$e\in\mathcal N$ を階数最小の非零射影とする。
$e\mathcal Ne=\mathbb Ce$ を主張する。$x=x^*\in e\mathcal Ne$ とする。$x$ の非零固有値 $\lambda$ に対するスペクトル射影 $q$ は $q=p(x)$ と書ける。ここで $p$ は
$p(\lambda)=1$ かつ $\mu\in(\operatorname{spec}x\cup\{0\})\setminus\{\lambda\}$ に対し $p(\mu)=0$ を満たすラグランジュ多項式である
（定理 1.5(a) の証明）。$p(0)=0$ かつ $x=ex=xe$ であるから $q\in\mathcal N$ かつ $q=eqe$ であり、よって
$q$ の値域は $e$ の値域に含まれる。極小性により $q=e$ または $q=0$ である。したがって $x$ は高々一つの非零固有値をもち、その固有空間は $e$ の値域であるから、
$x=\lambda e$ である。$e\mathcal Ne$ の一般の元は、
$x_{1,2}=x_{1,2}^*\in e\mathcal Ne$ により $x_1+ix_2$ と書ける。

*ステップ2：$e$ は $E$ を生成する.* $F=\operatorname{span}\mathcal Ne E$ とおく。これは
$\mathcal N$ で不変であるから、その射影 $p_F$ は $\mathcal N'$ に属する（補題 2.1）。また $F$ は
$*$-環 $\mathcal N'$ でも不変である：$b\in\mathcal N'$ に対し $bAe\xi=Aeb\xi\in F$。よって補題 2.1 を
$\mathcal N'$ に適用すると $p_F\in\mathcal N''=\mathcal N$ である（定理 2.2）。したがって
$p_F\in Z(\mathcal N)=\mathbb C1$ であり、$p_F\neq0$（$eE\subseteq F$ だから）なので $F=E$ である。

*ステップ3：行列単位.* $v^*v=e$ を満たす $v\in\mathcal N$ を *$e$-等長元* と呼ぶ。このとき $vv^*$
は $\mathcal N$ の射影である。（まず $v=ve$ である：任意の $x$ に対し
$\|v(1-e)x\|^2=\langle(1-e)x,v^*v(1-e)x\rangle=\langle(1-e)x,e(1-e)x\rangle=0$。よって
$vv^*v=ve=v$ かつ $(vv^*)^2=vev^*=vv^*$ である。以下の $v_i$ についてはこれから $e_iv_i=v_i$ および
$v_i^*e_i=v_i^*$ が従う。）値域 $e_i=v_iv_i^*$ が互いに直交する
$e$-等長元 $v_1=e,v_2,\dots,v_k$ を、$k$ が最大となるように選ぶ
（値域は同じ次元 $\operatorname{rank}e$ をもつ直交部分空間であるから、
$k\le\dim E/\operatorname{rank}e$ である）。$\sum_ie_i=1$ を主張する。そうでなければ $g=1-\sum_ie_i\neq0$ は
$\mathcal N$ の射影である。ステップ2により、ある $A\in\mathcal N$ に対して $gAe\neq0$ である（すべての $A$ について $gAe=0$ ならば、
すべての $A,\xi$ について $gAe\xi=0$ となり、$g$ は $F=E$ 上で消える）。$y=gAe$ とおく。すると $y^*y\in e\mathcal Ne$ は正かつ非零であるから、
$\lambda>0$ により $y^*y=\lambda e$ となり、$v=y/\sqrt\lambda$ は値域が
$gE$ に含まれ、すべての $e_i$ と直交する $e$-等長元である——これは最大性に矛盾する。

$i\neq j$ に対し $v_i^*v_j=v_i^*e_ie_jv_j=0$ である。$E_{ij}=v_iv_j^*$ とおく。すると
$E_{ij}E_{kl}=v_iv_j^*v_kv_l^*=\delta_{jk}v_iev_l^*=\delta_{jk}E_{il}$、
$E_{ij}^*=E_{ji}$ かつ $\sum_iE_{ii}=1$ である。$A\in\mathcal N$ に対し
$A=\sum_{ij}e_iAe_j=\sum_{ij}v_i(v_i^*Av_j)v_j^*$ かつ $v_i^*Av_j\in e\mathcal Ne=\mathbb Ce$ であるから、
$a_{ij}\in\mathbb C$ により $A=\sum_{ij}a_{ij}E_{ij}$ と書ける。よって $\mathcal N=\operatorname{span}\{E_{ij}\}$ である。

*ステップ4：ユニタリ.* $m=\operatorname{rank}e$ とし、$f_1,\dots,f_m$ を
$eE$ の正規直交基底とする。ベクトル $v_if_r$ は正規直交である：
$\langle v_if_r,v_jf_s\rangle=\langle f_r,v_i^*v_jf_s\rangle=\delta_{ij}\langle f_r,ef_s\rangle=\delta_{ij}\delta_{rs}$。
$E=\bigoplus_ie_iE$ かつ $e_iE=v_ieE$ であるから、これらは $E$ を張る。$V(|i\rangle\otimes f_r)=v_if_r$ と定義する。ここで $|1\rangle,\dots,|k\rangle$ は
$\mathbb C^k$ の標準基底であり、$f_r$ を通じて $\mathbb C^m$ を $eE$ と同一視する。
すると $E_{ij}v_lf_r=\delta_{jl}v_if_r$、すなわち $V^*E_{ij}V=|i\rangle\langle j|\otimes1_m$ であり、
$V^*\mathcal NV=M_k(\mathbb C)\otimes1_m$ である。

*ステップ5：可換子環.* $X\in M_k\otimes M_m$ をブロック行列
$X=\sum_{ij}|i\rangle\langle j|\otimes X_{ij}$ として書く。固定した $a,b$ に対し
$(|a\rangle\langle b|\otimes1)X=\sum_j|a\rangle\langle j|\otimes X_{bj}$ かつ
$X(|a\rangle\langle b|\otimes1)=\sum_i|i\rangle\langle b|\otimes X_{ia}$ である。両者が等しければ、
$j\neq b$ として位置 $(a,j)$ のブロックを比較すると $X_{bj}=0$ を得て、位置
$(a,b)$ のブロックを比較すると $X_{bb}=X_{aa}$ を得る。$a,b$ は任意であるから、$Y=X_{11}$ として $X=1_k\otimes Y$ である。
逆に任意の $1_k\otimes Y$ は $M_k\otimes1_m$ と可換である。$\square$

**定理 2.6（有限次元フォン・ノイマン環の構造）.** $\mathcal M$ を
$\mathbb C^n$ 上のフォン・ノイマン環とする。$C\ge1$、$\sum_ck_cm_c=n$ を満たす整数
$k_c,m_c\ge1$、およびユニタリ
$V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n$ が存在して
$$
V^*\mathcal MV=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},\qquad
V^*\mathcal M'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),
$$
かつ $V^*Z(\mathcal M)V=\operatorname{span}\{V^*P_1V,\dots,V^*P_CV\}$ となる。ここで $P_1,\dots,P_C$ は
$\mathcal M$ の極小中心射影（補題 2.3）であり、$V^*P_cV$ は第 $c$ 直和成分への射影である。組 $(k_c,m_c)$ は順序を除いて一意である：$P_c$ は一意であり（補題 2.3）、
$k_c^2=\dim\mathcal MP_c$ かつ $m_c=\operatorname{rank}P_c/k_c$ である。特に $*$-環として $\mathcal M\cong\bigoplus_cM_{k_c}(\mathbb C)$ であり、
$\dim\mathcal M=\sum_ck_c^2$、$\dim\mathcal M'=\sum_cm_c^2$ である。

*証明.* 補題 2.3 とその後の注意により、$E_c$ 上に作用する
$\mathcal M_c$ によって $\mathcal M=\bigoplus_c\mathcal M_c$ であり、各 $\mathcal M_c$ は因子環である（補題 2.4）。補題 2.5 により
$V_c^*\mathcal M_cV_c=M_{k_c}\otimes1_{m_c}$ となるユニタリ $V_c$ が得られる。
$V=\bigoplus_cV_c$ とおく。可換子環について：$X\in\mathcal M'$ は各 $P_c$ と可換であるから
$X=\sum_cP_cXP_c$ であり、$P_cXP_c$ は $\mathcal M_c$ と可換である。逆にそのような和はすべて
$\mathcal M$ と可換である。補題 2.5 により $E_c$ 上の $\mathcal M_c'$ が同定される。$\square$

数 $k_c$ は環そのものを記述する。**重複度** $m_c$ は環が
$M_n(\mathbb C)$ の中にどのように置かれているかを記述する：それらは代数的構造からは見えないが、
$\mathcal M$ に制限した $\mathbb C^n$ のトレースを決定する。
$$
\operatorname{Tr}_{\mathbb C^n}(A)=\sum_cm_c\operatorname{Tr}(A_c),\qquad A=V\Big(\bigoplus_cA_c\otimes1_{m_c}\Big)V^*.
$$
これが、ライブラリがトレースの族を扱う理由である（第3章）。

> **物理的意味.** $P_c$ は **超選択セクター** あるいは **対称性
> セクター** である：$\mathcal M$ のいかなるオブザーバブルも異なるセクターを結びつけないので、
> 「どのセクターか」の測定は他のすべての測定と両立する。中心は系の
> *古典的* な部分である。重複度 $m_c$ は縮退度である：ブロック $A_c$ の各
> 固有値は、$A_c$ における重複度の $m_c$ 倍の重複度で $\mathbb C^n$ に現れる。

## 2.4 例：対称性の群とその可換子環

有限群 $G$ がユニタリ表現 $\pi$ によって $\mathbb C^n$ に作用しているとする。対称性と可換な作用素全体は
$\pi(G)'$ をなし、これはフォン・ノイマン環である。その
可換子環 $\pi(G)''$ は $\pi(g)$ たちの張る空間である（その張る空間に定理 2.2 を適用する。これは
$*$-環である：$\pi(e)=1$、積 $\pi(g)\pi(h)=\pi(gh)$、随伴
$\pi(g)^*=\pi(g^{-1})$ を含む）。

**命題 2.7.** 定理 2.6 を $\mathcal M=\pi(G)''$ に適用する。すると
$V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$ であり、ここで各 $\rho_c$ は $\mathbb C^{k_c}$ 上の
$G$ の既約ユニタリ表現であり、$c\neq d$ に対し $\rho_c,\rho_d$ は同値でない。したがって
$$
\pi\cong\bigoplus_c\rho_c\otimes1_{m_c},\qquad
\pi(G)''\cong\bigoplus_cM_{\dim\rho_c},\qquad
\pi(G)'\cong\bigoplus_cM_{m_c}:
$$
対称な作用素 $\pi(G)'$ のブロックの大きさは既約表現の重複度であり、
$k_c$ と $m_c$ の役割は環とその可換子環の間で入れ替わる。

*証明.* $\pi(g)\in\mathcal M$ であるから、ある
$\rho_c(g)\in M_{k_c}$ により $V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$ と書け、$\pi$ がユニタリ表現であるから $g\mapsto\rho_c(g)$ もユニタリ表現である。
$\pi(g)$ たちは $\mathcal M$ を張るので、$\rho_c(g)$ たちは $M_{k_c}$ を張り、より強く、写像
$\sum_ga_g\pi(g)\mapsto\big(\sum_ga_g\rho_c(g)\big)_c$ は $\bigoplus_cM_{k_c}$ の上への全射である。
*既約性：* すべての $\rho_c(g)$ で不変な $\mathbb C^{k_c}$ の部分空間は、それらの張る $M_{k_c}$ で不変であるから、$0$ または $\mathbb C^{k_c}$ である。*非同値性：* ある可逆な $T$ と $c\neq d$ に対して、すべての $g$ について
$\rho_d(g)=T\rho_c(g)T^{-1}$ ならば、像の任意の元
$(x_{c'})_{c'}$ は $x_d=Tx_cT^{-1}$ を満たし、これは
$x_c=1$, $x_d=0$ なる元を除外する——全射性に矛盾する。最後の二つの同型は
定理 2.6 である。$\square$

表現の既約表現への分解は同値を除いて一意であり、したがって $m_c$ はいかなる分解においても *その* 重複度である。*証明.* $\rho,\rho'$ を既約とし、$T$ を絡作用素、すなわちすべての $g$ に対して $T\rho(g)=\rho'(g)T$ を満たすものとする。$T$ の核と像は不変であるから、$T=0$ であるか $T$ は可逆である（シューアの補題）。$\rho'=\rho$ ならば、$T$ の固有値 $\lambda$ に対して $T-\lambda1$ に同じことが当てはまるので、$T=\lambda1$ である。よって $\rho$ から $\rho'$ への絡作用素の空間の次元は、$\rho\cong\rho'$ なら $1$、そうでなければ $0$ である。既約な $\sigma_j$（重複を許す）により $\pi\cong\bigoplus_j\sigma_j$ と書けば、既約な $\rho$ から $\pi$ への絡作用素の空間の次元は $\#\{j:\sigma_j\cong\rho\}$ となり、これは $\pi$ と $\rho$ のみで決まる数である。$\rho_c$ が互いに非同値である命題 2.7 の分解では、この数は $\rho\cong\rho_c$ なら $m_c$ であり、$\rho$ がどの $\rho_c$ とも同値でなければ $0$ である。$\square$

$G=S_3$ が $(\mathbb C^2)^{\otimes3}$ の三つのテンソル因子を置換する場合、すなわち三つの
量子ビットあるいは三つのスピン $\tfrac12$ の場合を考える。分解は次のとおりである：対称部分空間（全スピン $\tfrac32$、次元 $4$）は重複度 $4$ の自明表現を担い、その補空間（次元 $4$）は重複度 $2$ の二次元標準表現を担い、符号表現は現れない。したがって
$\pi(S_3)''\cong M_1\oplus M_2$（次元 $5$）かつ $\pi(S_3)'\cong M_4\oplus M_2$
（次元 $20$）である。

*証明.* $\pi(S_3)''\cong\bigoplus_cM_{k_c}\otimes1_{m_c}$ と書く（定理 2.6）。すると $\sum_ck_cm_c=8$、$\sum_ck_c^2=\dim\pi(S_3)''=5$、$\sum_cm_c^2=\dim\pi(S_3)'=20$ である。二つの次元は下のコードが厳密な線形代数によって計算する（集合 $S$ の可換子環は線形写像 $X\mapsto(AX-XA)_{A\in S}$ の核である）。数 $20$ は手で数えることもできる：$X$ がすべての $\pi(g)$ と可換であることは、その成分 $X_{xy}$（$x,y\in\{0,1\}^3$）が $S_3$ の対 $(x,y)$ への作用、すなわち $\{0,1\}^2$ の元の三つ組への作用の軌道上で一定であることと同値であり、そのような軌道は $\binom63=20$ 個ある（$4$ 個の元から作る大きさ $3$ の多重集合）。$\sum_ck_c^2=5$ から $k=(1,2)$ または $k=(1,1,1,1,1)$ である。前者の場合、$m_1+2m_2=8$ と $m_1^2+m_2^2=20$ から $m=(4,2)$ を得る。後者の場合、和が $8$ で平方和が $20$ の五つの正の整数は $m=(4,1,1,1,1)$ でなければならない。しかしそのとき $\pi(S_3)''$ は可換になるが、互換 $(12)$ と $(23)$ は可換でない：$\pi((12))\pi((23))|001\rangle=|100\rangle$ であるが $\pi((23))\pi((12))|001\rangle=|010\rangle$ である。よって $k=(1,2)$、$m=(4,2)$ である。一次元表現 $\rho_1$ は自明表現か符号表現である。符号表現ではない：各基底ベクトル $|b_1b_2b_3\rangle$ は二つの等しいビットをもつから、ある互換 $t$ で固定され、$A\pi(t)=-A$ を満たす $A=\sum_g\operatorname{sgn}(g)\pi(g)$ について $A|b_1b_2b_3\rangle=A\pi(t)|b_1b_2b_3\rangle=-A|b_1b_2b_3\rangle$ となる。よって $A=0$ であり、符号表現は現れない。したがって $\rho_1$ は自明表現で、その $4$ 個のコピーは対称部分空間を張り、$\rho_2$ は二次元の既約表現（命題 2.7）、すなわち標準表現である。$\square$

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
assert len(comm) == 20                                   # pi(S3)'  = M_4 (+) M_2
assert span_dim(group) == 5                              # pi(S3)'' = M_1 (+) M_2
assert len(commutant_basis(comm)) == 5                   # double commutant theorem
```

重複度はスペクトルに現れる。群平均
$X\mapsto\frac1{|G|}\sum_g\pi(g)X\pi(g)^*$ は任意の作用素を可換子環に写す
（演習問題 3）。したがって $\pi(S_3)'$ のランダムな自己共役元は（ほとんど確実に）対称部分空間上に $4$ 個の単純
固有値をもち、その補空間上に $2$ 個の二重縮退した固有値をもつ。そこでの $\pi(S_3)'$ のブロックが $1_2\otimes M_2$ だからである。コードは縮退した対を確かめる。

```python
def group_average(X):
    return sum(g @ X @ g.conj().T for g in group) / len(group)

H0 = torch.randn(8, 8, dtype=dt)
Hsym = group_average(H0 + H0.conj().T)                  # a random symmetric Hamiltonian
assert all(torch.allclose(g @ Hsym, Hsym @ g) for g in group)

P_sym = sum(group) / len(group)                         # projection onto the symmetric subspace
w, Q = torch.linalg.eigh(P_sym)
Q_rest = Q[:, w.real < 0.5]                             # orthonormal basis of the complement (dim 4)
e_rest = torch.linalg.eigvalsh(Q_rest.conj().T @ Hsym @ Q_rest)
print(e_rest)
assert Q_rest.shape[1] == 4
assert torch.allclose(e_rest[0], e_rest[1]) and torch.allclose(e_rest[2], e_rest[3])
```

## 2.5 ライブラリにおける環

ライブラリは定理 2.6 の抽象的な環 $\bigoplus_cM_{k_c}(\mathbb C)$ を
$\bigoplus_c\mathbb C^{k_c}$ 上で、重複度なしに格納する：作用素はそのブロック $A_c$ のリストであり、共通のサイズにパディングされる。上の例の可換子環は
`TypeIAlgebra([4, 2], [4, 2])` である。

```python
from torch_vn_algebra import TypeIAlgebra, channels

alg = TypeIAlgebra([4, 2], [4, 2], precision='double')    # M_4 (+) M_2
print(alg.k_factors, alg.C)

A = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
print(A.matrix.shape)              # (1, 2, 4, 4): block 1 occupies the top-left 2 x 2 corner
assert A.matrix[0, 1, 2:, :].abs().max() == 0 and A.matrix[0, 1, :, 2:].abs().max() == 0
```

中心は元 $\sum_cz_cP_c$ からなる。`alg.central` がそれらを構成し、
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

## 2.6 条件付き期待値

$\mathcal M$ からフォン・ノイマン部分環
$\mathcal N\subseteq\mathcal M$ の上への **条件付き期待値** とは、線形写像 $E:\mathcal M\to\mathcal N$ であって、正
（$A\ge0\Rightarrow E(A)\ge0$）であり、$B\in\mathcal N$ に対して $E(B)=B$ を満たし、かつ
$\mathcal N$-双加群写像であるもの：$B_1,B_2\in\mathcal N$ に対し $E(B_1AB_2)=B_1E(A)B_2$。

**命題 2.8.** $\mathcal M=\bigoplus_cM_{k_c}$ および $\mathcal N=Z(\mathcal M)$ に対し、
$$
E(A)=\sum_c\frac{\operatorname{Tr}A_c}{k_c}\,P_c
$$
は条件付き期待値である。これは
$\sum_cw_c\operatorname{Tr}A_c$ の形のすべてのトレースを保つ。

*証明.* 線形性は明らかである。$A\ge0$ ならば $\operatorname{Tr}A_c\ge0$ であるから $E(A)\ge0$ である。
$B=\sum_cb_cP_c$ に対し $E(B)=\sum_cb_ck_cP_c/k_c=B$ である。$B_i=\sum_cb^{(i)}_cP_c$ に対し
$\operatorname{Tr}(B_1AB_2)_c=b^{(1)}_cb^{(2)}_c\operatorname{Tr}A_c$ であり、これから双加群性が従う。最後に、すべての $c$ に対し $\operatorname{Tr}(E(A)_c)=\operatorname{Tr}A_c$ である。$\square$

ライブラリではこの写像はチャネル `channels.center_expectation(alg)` である。これを作用素に適用すると古典的な（セクター平均した）部分が得られる。

```python
E = channels.center_expectation(alg)
EA = E(A)
assert close(E(EA), EA)                                  # idempotent
assert close(E(Z @ A @ P1), Z @ EA @ P1)                 # bimodule property
A2 = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
assert close(EA @ A2, A2 @ EA)                           # E(A) commutes with a fresh element
for c in range(alg.C):                                   # Tr E(A)_c = Tr A_c
    assert abs(EA.matrix[0, c].diagonal().sum() - A.matrix[0, c].diagonal().sum()) < 1e-10
```

2.4節の群平均は $M_n(\mathbb C)$ の中での可換子環
$\pi(G)'$ の上への条件付き期待値であり、$\operatorname{Tr}$ を保つ（演習問題 3）。
（条件付き期待値は一般には一意でない：$\mathbb C1$ の上へのものは、任意の状態
$\varphi$ から $A\mapsto\varphi(A)1$ として一つ得られる。）両者は同じ考え方の実例である：
ある構造を「見ない」作用素への射影である。

> **物理的意味.** 完全な行列環 $M_D(\mathbb C)$, $D=\sum_ck_c$ の上で、写像 $X\mapsto\sum_cP_cXP_c$ は
> $\mathcal M$ の上への条件付き期待値である：これは超選択セクター間のコヒーレンスを消す。電荷（セクター）を
> 測定して結果を記録しないときに起こるのがまさにこれである。中心の上への条件付き期待値 $E$ はさらに先へ進む：
> 古典的な情報、すなわち各セクターの重み $p_c$ だけを残し、各セクター内の状態を最大混合状態に置き換える。
> 状態 $\omega$ と状態 $\omega\circ E$ は、「どのセクターか？」だけを問うすべての観測量の上で一致する。

## 2.7 比較の対象

2.4節の例は任意の線形代数パッケージで確かめることができる。
表現論のソフトウェア（GAP, SageMath）は分解
$\pi\cong\bigoplus\rho_\lambda\otimes1_{m_\lambda}$ を記号的に与える。物理学の厳密対角化コード（QuSpin、QuTiP の対称性モジュール）は同じブロック構造を用いて、ハミルトニアンをセクターごとに対角化する。本ライブラリは構造を出発点とする：ひとたびブロックがわかれば、その後のすべての計算——関数、トレース、状態、チャネル、ダイナミクス——は環の中にとどまる。

## 演習問題

1. $M_k(\mathbb C)$ の中心が $\mathbb C1$ であること（行列単位
   $|i\rangle\langle j|$ を用いよ）、および $Z(\bigoplus_cM_{k_c})$ の次元が $C$ であることを示せ。
2. $S$ を $\mathbb C^3$ 上の巡回シフト $Se_i=e_{i+1\bmod3}$ とする。
   $\{S\}'$ と $\{S\}''$ を決定し、次元を数値的に確かめよ。
3. 有限群のユニタリ表現 $\pi$ に対し、
   $E(X)=\frac1{|G|}\sum_g\pi(g)X\pi(g)^*$ が
   $M_n(\mathbb C)$ から $\pi(G)'$ の上への、$\operatorname{Tr}$ を保つ条件付き期待値であることを示せ。

<details>
<summary>解答</summary>

1. $X$ がすべての $|i\rangle\langle j|$ と可換ならば
   $X|i\rangle\langle j|=|i\rangle\langle j|X$ である。両辺を $e_j$ に適用して比較すると、
   すべての $i,j$ に対し $Xe_i=X_{jj}e_i$ を得るので、$X=X_{11}1$ である。直和においては、元
   $\bigoplus X_c$ が中心的であるのはすべての $X_c$ が中心的であるときかつそのときに限るから、中心は
   $\{\bigoplus_cz_c1_{k_c}\}$ であり、その次元は $C$ である。

2. $S$ は固有値 $1,\omega,\omega^2$（$\omega=e^{2\pi i/3}$）をもつユニタリであり、各固有値は単純である。
   $S$ と可換な作用素はその固有空間を保つので、固有基底において
   対角である：$\{S\}'\cong\mathbb C^3$、すなわち巡回行列全体である。すると $\{S\}''$ も
   $\cong\mathbb C^3$ である。なぜなら $S$ の固有基底において、すべての
   対角行列と可換な行列は対角だからである（相異なる対角成分をもつ $D$ について $XD$ と $DX$ の成分を比較せよ）。それは $1,S,S^2$ によって張られる（ヴァンデルモンド：固有値は相異なる）。

```python
S = torch.roll(torch.eye(3, dtype=dt), 1, dims=0)
assert len(commutant_basis([S])) == 3
assert len(commutant_basis(commutant_basis([S]))) == 3
assert span_dim([torch.eye(3, dtype=dt), S, S @ S]) == 3
```

3. $E$ は線形かつ正である（写像 $X\mapsto UXU^*$ の和）。$h\in G$ に対し
   $\pi(h)E(X)\pi(h)^*=\frac1{|G|}\sum_g\pi(hg)X\pi(hg)^*=E(X)$ であるから $E(X)\in\pi(G)'$ である。
   $X\in\pi(G)'$ ならば各項は $X$ に等しいので $E(X)=X$ である。$B_1,B_2\in\pi(G)'$ に対し
   $\pi(g)B_1XB_2\pi(g)^*=B_1\pi(g)X\pi(g)^*B_2$ であり、これから双加群性が従う。
   最後に、すべての $g$ に対し $\operatorname{Tr}(\pi(g)X\pi(g)^*)=\operatorname{Tr}X$ であるから、
   $\operatorname{Tr}E(X)=\operatorname{Tr}X$ である。

</details>

**次へ：**[第3章](03_traces_and_states.md) では
$\bigoplus_cM_{k_c}$ 上のトレースを分類し、状態と密度を導入し、後で用いるエントロピーとトレースの不等式を証明する。
