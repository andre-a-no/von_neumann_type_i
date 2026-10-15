[English](../02_von_neumann_algebras.md) · [Español](../es/02_von_neumann_algebras.md) · [Français](../fr/02_von_neumann_algebras.md) · [Deutsch](../de/02_von_neumann_algebras.md) · [中文](../zh/02_von_neumann_algebras.md) · **日本語** · [Русский](../ru/02_von_neumann_algebras.md)

# 2. 有限次元フォン・ノイマン環

フォン・ノイマン環とは、「可換子環を二度とる」操作で閉じた作用素環である。有限次元ではこれは純粋に代数的な条件であり、本章の主定理は、そのような環はすべて、ユニタリな基底変換を除いて、全行列環の直和（各成分はある重複度で繰り返される）であることを述べる。これがライブラリの土台となる構造である。

*ライブラリでは：* 環 $`\bigoplus_cM_{k_c}(\mathbb C)`$ はクラス `TypeIAlgebra` である — [チュートリアル、2.5節](../../tutorial/ja/02_von_neumann_algebras.md#s-2-5)。

**本章で学ぶこと**：二重可換子環定理、完全な証明付きの構造定理、対称性の群がどのようにしてそのような環を生み出すか、中心とは何か、そして条件付き期待値とは何か。本章の数値的な検算は
[チュートリアル 2](../../tutorial/ja/02_von_neumann_algebras.md) にある。

本章を通じて、*$`\mathbb C^n`$ 上の $`*`$-環* とは、積と随伴について閉じていて単位元 $`1`$ を含む複素線形部分空間
$`\mathcal M\subseteq M_n(\mathbb C)`$ を意味する。

<a id="s-2-1"></a>

## 2.1 可換子環と二重可換子環定理

集合 $`S\subseteq M_n(\mathbb C)`$ に対し、その **可換子環** を
$`S'=\{X\in M_n(\mathbb C):XA=AX\ \forall A\in S\}`$ と定める。これは常に単位的な環である。
$`S`$ が随伴について閉じている（$`S^*=S`$）ならば $`S'`$ は $`*`$-環である。実際、
すべての $`A\in S`$ に対する $`XA=AX`$ から、すべての $`A\in S`$ に対して $`X^*A^*=A^*X^*`$ が従う。$`S`$ の任意の元は $`S'`$ の任意の元と可換であるから $`S\subseteq S''`$ である。これを $`S'`$ に適用すると
$`S'\subseteq S'''`$ を得て、$`S\subseteq S''`$ の可換子環をとる（包含関係は逆転する）と
$`S'''\subseteq S'`$ を得る。したがって $`S'=S'''`$ である。

*計算：* 零空間としての可換子環の数値計算 — [チュートリアル、2.4節](../../tutorial/ja/02_von_neumann_algebras.md#s-2-4)。*ライブラリでは：* `decompose.commutant_dimension` — [チュートリアル、4.2節](../../tutorial/ja/04_symmetries.md#s-4-2)；`decompose.commutant_basis` — [チュートリアル、4.3節](../../tutorial/ja/04_symmetries.md#s-4-3)；雑音を含む行列に対して零空間のどの許容誤差が正しいかは定理 4.25 — [チュートリアル、4.9節](../../tutorial/ja/04_symmetries.md#s-4-9)。

環が $`\mathbb C^n`$ の部分空間 $`E`$ に作用するときは、正規直交基底によって $`E`$ を
$`\mathbb C^{\dim E}`$ と同一視し、$`\mathbb C^n`$ 上で述べたすべての結果がそれに適用できるようにする。

**定義.** $`\mathbb C^n`$ 上の **フォン・ノイマン環** とは、$`\mathcal M=\mathcal M''`$ を満たす $`*`$-環
$`\mathcal M`$ のことである。その **中心** は
$`Z(\mathcal M)=\mathcal M\cap\mathcal M'`$ である。**因子環**（factor）とは、
$`Z(\mathcal M)=\mathbb C1`$ を満たすフォン・ノイマン環のことである。

この定義は無限次元ヒルベルト空間上の作用素のためのものであり、本講義では無限次元は扱わない。有限次元では条件 $`\mathcal M=\mathcal M''`$ は自動的に成り立つ。

<a id="r-2-1"></a>

**補題 2.1.** $`\mathcal M`$ を $`\mathbb C^n`$ 上の $`*`$-環とし、$`K\subseteq\mathbb C^n`$ を
$`\mathcal MK\subseteq K`$ を満たす部分空間とする。このとき $`K`$ への直交射影 $`P`$ は
$`\mathcal M'`$ に属する。

*証明.* $`A\in\mathcal M`$ に対し、$`AK\subseteq K`$ かつ $`A^*K\subseteq K`$ である（$`A^*\in\mathcal M`$ だから）。
第二の包含から $`AK^\perp\subseteq K^\perp`$ が従う：$`y\in K^\perp`$, $`x\in K`$ に対し
$`\langle x,Ay\rangle=\langle A^*x,y\rangle=0`$ である。したがって $`A`$ は $`K`$ を $`K`$ に、$`K^\perp`$
を $`K^\perp`$ に写し、これは $`AP=PAP=PA`$ を意味する。$`\square`$

<a id="r-2-2"></a>

**定理 2.2（二重可換子環定理）.** $`\mathbb C^n`$ 上の任意の $`*`$-環 $`\mathcal M`$ は
$`\mathcal M=\mathcal M''`$ を満たす。

*証明.* $`\mathcal M\subseteq\mathcal M''`$ である。$`T\in\mathcal M''`$ とする。

*ステップ1：一つのベクトル.* $`x\in\mathbb C^n`$ に対し、部分空間 $`K=\mathcal Mx`$ は
$`\mathcal M`$ で不変であり、$`x=1x`$ を含む。補題 2.1 によりその射影 $`P`$ は $`\mathcal M'`$ に属するから、
$`TP=PT`$ かつ $`Tx=TPx=PTx\in K`$ である：$`Tx=Ax`$ を満たす $`A\in\mathcal M`$ が存在する。

*ステップ2：すべてのベクトルを同時に.* $`\mathbb C^n\otimes\mathbb C^n=(\mathbb C^n)^{\oplus n}`$ 上で
$`*`$-環 $`\tilde{\mathcal M}=\{A\oplus\dots\oplus A:A\in\mathcal M\}`$ を考える。ブロック行列
$`[X_{ij}]_{i,j=1}^n`$ がすべての $`A\oplus\dots\oplus A`$ と可換であるための必要十分条件は、すべての $`i,j`$ に対して $`X_{ij}A=AX_{ij}`$ となることであるから、$`\tilde{\mathcal M}'=\{[X_{ij}]:X_{ij}\in\mathcal M'\}`$ である。作用素
$`\tilde T=T\oplus\dots\oplus T`$ はそのような任意の $`[X_{ij}]`$ と可換である。なぜなら
$`TX_{ij}=X_{ij}T`$ だからである。よって $`\tilde T\in\tilde{\mathcal M}''`$ である。ステップ1を $`\tilde{\mathcal M}`$、
$`\tilde T`$ およびベクトル $`\xi=(e_1,\dots,e_n)`$ に適用すると、
$`(Te_1,\dots,Te_n)=(Ae_1,\dots,Ae_n)`$、すなわち $`T=A\in\mathcal M`$ となる $`A\in\mathcal M`$ が存在する。$`\square`$

*計算：* 例における定理 2.2 — [チュートリアル、2.4節](../../tutorial/ja/02_von_neumann_algebras.md#s-2-4)。

したがって有限次元では「フォン・ノイマン環」と「行列からなる単位的 $`*`$-環」は同じものである。それでも二重可換子環定理は正しい観点を与える：物理的には、オブザーバブルの環はきわめて多くの場合、可換子環として *与えられる* ——すなわち対称性と可換な作用素全体として（2.4節）。

<a id="s-2-2"></a>

## 2.2 中心

<a id="r-2-3"></a>

**補題 2.3.** $`\mathcal Z`$ を $`\mathbb C^n`$ 上の可換な $`*`$-環とする。零でなく互いに直交し
$`\sum_cP_c=1`$ を満たす射影 $`P_1,\dots,P_C\in\mathcal Z`$ であって、$`\mathcal Z=\mathrm{span}\{P_1,\dots,P_C\}`$ となるものが存在する。それらは $`\mathcal Z`$ の極小な非零射影であり、順序を除いて一意である。

*証明.* 任意の $`z\in\mathcal Z`$ は正規である（$`z^*\in\mathcal Z`$ は $`z`$ と可換）。
$`z_1,\dots,z_r`$ を $`\mathcal Z`$ を張る元とする。$`\mathbb C^n`$ を同時固有空間に分解する。
$`z_1`$ の固有空間は互いに直交し $`\mathbb C^n`$ を張る（定理 1.3）。それらは
$`z_2,\dots,z_r`$ およびその随伴（これらは $`z_1`$ と可換）で不変であるから、各
$`z_j`$ の $`z_1`$ の固有空間への制限は再び正規であり、$`z_1`$ の各固有空間を $`z_2`$ の固有空間に分解でき、以下同様である。その結果、直交分解
$`\mathbb C^n=E_1\oplus\dots\oplus E_C`$ が得られる。ここで各 $`E_c`$ は零でない部分空間であり、その上で任意の
$`z\in\mathcal Z`$ はスカラー $`\chi_c(z)`$ として作用する。ただし同じ関数
$`\chi_c:\mathcal Z\to\mathbb C`$ をもつ部分空間は併合しておき、$`c\neq d`$ に対し $`\chi_c\neq\chi_d`$ となるようにする。
$`P_c`$ を $`E_c`$ への射影とする。すると任意の $`z\in\mathcal Z`$ に対し $`z=\sum_c\chi_c(z)P_c`$ である。

残るは $`P_c\in\mathcal Z`$ を示すことである。$`d\neq c`$ に対し
$`\chi_c(z_{cd})\neq\chi_d(z_{cd})`$ となる $`z_{cd}\in\mathcal Z`$ を選ぶ。積
$`\prod_{d\neq c}\big(z_{cd}-\chi_d(z_{cd})1\big)/\big(\chi_c(z_{cd})-\chi_d(z_{cd})\big)`$
は $`\mathcal Z`$ に属し、$`E_c`$ 上では $`1`$ として、すべての $`E_d`$, $`d\neq c`$ 上では $`0`$ として作用する。すなわちこれは
$`P_c`$ である（$`C=1`$ のときは $`P_1=1`$ であり、空積は $`1`$ である）。極小性と一意性：射影 $`\sum_c a_cP_c`$ は $`a_c\in\{0,1\}`$ を満たすから、
$`\mathcal Z`$ の非零射影は $`P_c`$ たちの空でない集合の和である。$`\square`$

これを $`\mathcal Z=Z(\mathcal M)`$ に適用したとき、射影 $`P_c`$ を $`\mathcal M`$ の **極小中心射影** という。任意の $`A\in\mathcal M`$ は
$`A=\sum_cAP_c`$（$`AP_c=P_cAP_c`$）と分解されるから、
$`\mathcal M=\bigoplus_c\mathcal MP_c`$ である（$`E_c=P_c\mathbb C^n`$ 上に作用する環の直交直和）。

*ライブラリでは：* 中心の元 $`\sum_cz_cP_c`$（極小中心射影を含む）は `TypeIAlgebra.central` で構成される — [チュートリアル、2.5節](../../tutorial/ja/02_von_neumann_algebras.md#s-2-5)。

<a id="r-2-4"></a>

**補題 2.4.** 各 $`\mathcal M_c=\mathcal MP_c`$ は、単位元 $`P_c`$ をもつ $`E_c`$ 上の環として、因子環である。

*証明.* $`y\in\mathcal M_c`$ が $`\mathcal M_c`$ と可換であるとする。$`A\in\mathcal M`$ に対し、$`P_c`$ が中心的であることと $`y=yP_c=P_cy`$ を用いて
$`Ay=AP_cy=(AP_c)y=y(AP_c)=yP_cA=yA`$ である。
よって $`y\in Z(\mathcal M)`$、$`y=\sum_d\chi_d(y)P_d`$ であり、$`y=yP_c=\chi_c(y)P_c`$ である。$`\square`$

<a id="s-2-3"></a>

## 2.3 構造定理

<a id="r-2-5"></a>

**補題 2.5（因子環）.** $`\mathcal N`$ を空間 $`E`$ 上の因子環とする。整数
$`k,m\ge1`$ とユニタリ $`V:\mathbb C^k\otimes\mathbb C^m\to E`$ であって
$`V^*\mathcal NV=M_k(\mathbb C)\otimes1_m`$ となるものが存在する。さらに $`V^*\mathcal N'V=1_k\otimes M_m(\mathbb C)`$ である。

*証明.* *ステップ1：極小射影.* $`\mathcal N`$ は零でない射影（例えば $`1`$）を含む。
$`e\in\mathcal N`$ を階数最小の非零射影とする。
$`e\mathcal Ne=\mathbb Ce`$ を主張する。$`x=x^*\in e\mathcal Ne`$ とする。$`x`$ の非零固有値 $`\lambda`$ に対するスペクトル射影 $`q`$ は $`q=p(x)`$ と書ける。ここで $`p`$ は
$`p(\lambda)=1`$ かつ $`\mu\in(\mathrm{spec}\,x\cup\{0\})\setminus\{\lambda\}`$ に対し $`p(\mu)=0`$ を満たすラグランジュ多項式である
（定理 1.5(a),(b)：$`p`$ は $`\mathrm{spec}\,x`$ 上で $`\{\lambda\}`$ の指示関数と一致する）。
$`p(0)=0`$ であるから、ある多項式 $`r`$ により $`p(t)=t\,r(t)`$ と書け、よって $`q=x\,r(x)=r(x)\,x\in\mathcal N`$ であり、$`x=ex=xe`$ から
$`q=eq=qe`$ を得る。よって $`q`$ の値域は $`e`$ の値域に含まれる。$`q\neq0`$ ならば、極小性により
$`\mathrm{rank}\,q\ge\mathrm{rank}\,e`$ であるから、二つの値域は一致し $`q=e`$ である。相異なる固有値に対するスペクトル
射影は互いに直交するから、$`x`$ は高々一つの非零固有値
$`\lambda`$ をもち、そのスペクトル射影は $`e`$ である。よって $`x=\lambda e`$ である（非零固有値がなければ $`x=0=0\cdot e`$）。$`e\mathcal Ne`$ の一般の元は、
$`x_{1,2}=x_{1,2}^*\in e\mathcal Ne`$ により $`x_1+ix_2`$ と書ける。

*ステップ2：$`e`$ は $`E`$ を生成する.* $`F=\mathrm{span}\,\mathcal Ne E`$ とおく。これは
$`\mathcal N`$ で不変であるから、その射影 $`p_F`$ は $`\mathcal N'`$ に属する（補題 2.1）。また $`F`$ は
$`*`$-環 $`\mathcal N'`$ でも不変である：$`b\in\mathcal N'`$ に対し $`bAe\xi=Aeb\xi\in F`$。よって補題 2.1 を
$`\mathcal N'`$ に適用すると $`p_F\in\mathcal N''=\mathcal N`$ である（定理 2.2）。したがって
$`p_F\in Z(\mathcal N)=\mathbb C1`$ であり、$`p_F\neq0`$（$`eE\subseteq F`$ だから）なので $`F=E`$ である。

*ステップ3：行列単位.* $`v^*v=e`$ を満たす $`v\in\mathcal N`$ を *$`e`$-等長元* と呼ぶ。このとき $`vv^*`$
は $`\mathcal N`$ の射影である。（まず $`v=ve`$ である：任意の $`x`$ に対し
$`\|v(1-e)x\|^2=\langle(1-e)x,v^*v(1-e)x\rangle=\langle(1-e)x,e(1-e)x\rangle=0`$。よって
$`vv^*v=ve=v`$ かつ $`(vv^*)^2=vev^*=vv^*`$ である。以下の $`v_i`$ についてはこれから $`e_iv_i=v_i`$ および
$`v_i^*e_i=v_i^*`$ が従う。）値域 $`e_i=v_iv_i^*`$ が互いに直交する
$`e`$-等長元 $`v_1=e,v_2,\dots,v_k`$ を、$`k`$ が最大となるように選ぶ
（値域は同じ次元 $`\mathrm{rank}\,e`$ をもつ直交部分空間であるから、
$`k\le\dim E/\mathrm{rank}\,e`$ である）。$`\sum_ie_i=1`$ を主張する。そうでなければ $`g=1-\sum_ie_i\neq0`$ は
$`\mathcal N`$ の射影である。ステップ2により、ある $`A\in\mathcal N`$ に対して $`gAe\neq0`$ である（すべての $`A`$ について $`gAe=0`$ ならば、
すべての $`A,\xi`$ について $`gAe\xi=0`$ となり、$`g`$ は $`F=E`$ 上で消える）。$`y=gAe`$ とおく。すると $`y^*y\in e\mathcal Ne`$ は正かつ非零であるから、
$`\lambda>0`$ により $`y^*y=\lambda e`$ となり、$`v=y/\sqrt\lambda`$ は値域が
$`gE`$ に含まれ、すべての $`e_i`$ と直交する $`e`$-等長元である——これは最大性に矛盾する。

$`i\neq j`$ に対し $`v_i^*v_j=v_i^*e_ie_jv_j=0`$ である。$`E_{ij}=v_iv_j^*`$ とおく。すると
$`E_{ij}E_{kl}=v_iv_j^*v_kv_l^*=\delta_{jk}v_iev_l^*=\delta_{jk}E_{il}`$、
$`E_{ij}^*=E_{ji}`$ かつ $`\sum_iE_{ii}=1`$ である。$`A\in\mathcal N`$ に対し
$`A=\sum_{ij}e_iAe_j=\sum_{ij}v_i(v_i^*Av_j)v_j^*`$ かつ $`v_i^*Av_j\in e\mathcal Ne=\mathbb Ce`$ であるから、
$`a_{ij}\in\mathbb C`$ により $`A=\sum_{ij}a_{ij}E_{ij}`$ と書ける。よって $`\mathcal N=\mathrm{span}\{E_{ij}\}`$ である。

*ステップ4：ユニタリ.* $`m=\mathrm{rank}\,e`$ とし、$`f_1,\dots,f_m`$ を
$`eE`$ の正規直交基底とする。ベクトル $`v_if_r`$ は正規直交である：
$`\langle v_if_r,v_jf_s\rangle=\langle f_r,v_i^*v_jf_s\rangle=\delta_{ij}\langle f_r,ef_s\rangle=\delta_{ij}\delta_{rs}`$。
$`E=\bigoplus_ie_iE`$ かつ $`e_iE=v_ieE`$ であるから、これらは $`E`$ を張る。$`V(|i\rangle\otimes f_r)=v_if_r`$ と定義する。ここで $`|1\rangle,\dots,|k\rangle`$ は
$`\mathbb C^k`$ の標準基底であり、$`f_r`$ を通じて $`\mathbb C^m`$ を $`eE`$ と同一視する。
すると $`E_{ij}v_lf_r=\delta_{jl}v_if_r`$、すなわち $`V^*E_{ij}V=|i\rangle\langle j|\otimes1_m`$ であり、
$`V^*\mathcal NV=M_k(\mathbb C)\otimes1_m`$ である。

*ステップ5：可換子環.* $`X\in M_k\otimes M_m`$ をブロック行列
$`X=\sum_{ij}|i\rangle\langle j|\otimes X_{ij}`$ として書く。固定した $`a,b`$ に対し
$`(|a\rangle\langle b|\otimes1)X=\sum_j|a\rangle\langle j|\otimes X_{bj}`$ かつ
$`X(|a\rangle\langle b|\otimes1)=\sum_i|i\rangle\langle b|\otimes X_{ia}`$ である。両者が等しければ、
$`j\neq b`$ として位置 $`(a,j)`$ のブロックを比較すると $`X_{bj}=0`$ を得て、位置
$`(a,b)`$ のブロックを比較すると $`X_{bb}=X_{aa}`$ を得る。$`a,b`$ は任意であるから、$`Y=X_{11}`$ として $`X=1_k\otimes Y`$ である。
逆に任意の $`1_k\otimes Y`$ は $`M_k\otimes1_m`$ と可換である。$`\square`$

<a id="r-2-6"></a>

**定理 2.6（有限次元フォン・ノイマン環の構造）.** $`\mathcal M`$ を
$`\mathbb C^n`$ 上のフォン・ノイマン環とする。$`C\ge1`$、$`\sum_ck_cm_c=n`$ を満たす整数
$`k_c,m_c\ge1`$、およびユニタリ
$`V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n`$ が存在して
```math
V^*\mathcal MV=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},\qquad
V^*\mathcal M'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),
```
かつ $`V^*Z(\mathcal M)V=\mathrm{span}\{V^*P_1V,\dots,V^*P_CV\}`$ となる。ここで $`P_1,\dots,P_C`$ は
$`\mathcal M`$ の極小中心射影（補題 2.3）であり、$`V^*P_cV`$ は第 $`c`$ 直和成分への射影である。組 $`(k_c,m_c)`$ は順序を除いて一意である：最初の性質をもつ *任意の* ユニタリ $`V`$ に対し、
$`V(\cdot)V^*`$ は直和成分への射影を $`P_c`$ たちに写す（証明を参照）。
よって $`k_c^2=\dim\mathcal MP_c`$ と $`m_c=\mathrm{rank}\,P_c/k_c`$ は $`\mathcal M`$ によって決まる。特に $`*`$-環として $`\mathcal M\cong\bigoplus_cM_{k_c}(\mathbb C)`$ であり、
$`\dim\mathcal M=\sum_ck_c^2`$、$`\dim\mathcal M'=\sum_cm_c^2`$ である。

*証明.* 補題 2.3 とその後の注意により、$`E_c`$ 上に作用する
$`\mathcal M_c`$ によって $`\mathcal M=\bigoplus_c\mathcal M_c`$ であり、各 $`\mathcal M_c`$ は因子環である（補題 2.4）。補題 2.5 により
$`V_c^*\mathcal M_cV_c=M_{k_c}\otimes1_{m_c}`$ となるユニタリ $`V_c`$ が得られる。
$`V=\bigoplus_cV_c`$ とおく。可換子環について：$`X\in\mathcal M'`$ は各 $`P_c`$ と可換であるから
$`X=\sum_cP_cXP_c`$ であり、$`P_cXP_c`$ は $`\mathcal M_c`$ と可換である。逆にそのような和はすべて
$`\mathcal M`$ と可換である。補題 2.5 により $`E_c`$ 上の $`\mathcal M_c'`$ が同定され、$`V`$ で共役をとると
$`V^*\mathcal M'V`$ の式を得る。

*一意性.* $`V`$ を $`V^*\mathcal MV=\mathcal B:=\bigoplus_{c=1}^{C'}M_{k_c}\otimes1_{m_c}`$ を満たす任意のユニタリとし、
$`Q_c`$ を第 $`c`$ 直和成分への射影とする。$`\mathcal B`$ の元 $`\bigoplus_cX_c`$（$`X_c\in M_{k_c}\otimes1_{m_c}`$）が
中心的であることは、各 $`X_c`$ が $`M_{k_c}\otimes1_{m_c}`$ と可換であること、すなわち（補題 2.5、ステップ5）$`X_c\in(M_{k_c}\otimes1_{m_c})\cap(1_{k_c}\otimes M_{m_c})=\mathbb C1`$
と同値である（元 $`a\otimes1=1\otimes b`$ はすべての対角ブロックが $`b`$ に等しく、すべての非対角ブロックが $`0`$ であるから、$`a`$ はスカラーである）。
よって $`Z(\mathcal B)=\mathrm{span}\{Q_1,\dots,Q_{C'}\}`$ であり、$`Z(\mathcal M)=VZ(\mathcal B)V^*`$ は、和が $`1`$ である零でない
互いに直交する射影 $`VQ_cV^*`$ によって張られる。そのそれぞれは $`Z(\mathcal M)`$ の非零射影であるから、（補題 2.3）
$`P_d`$ たちのある空でない集合 $`S_c`$ の和である。$`S_c`$ たちは互いに素であり（直交性）、すべての $`d`$ を覆う（和が $`1`$）。さらに
$`C'=\dim Z(\mathcal M)=C`$ であるから、各 $`S_c`$ は一元集合である。したがって番号を付け替えれば $`VQ_cV^*=P_c`$ であり、ゆえに
$`\dim\mathcal MP_c=\dim(M_{k_c}\otimes1_{m_c})=k_c^2`$ かつ $`\mathrm{rank}\,P_c=k_cm_c`$ である。$`\square`$

*ライブラリでは：* 環 $`\bigoplus_cM_{k_c}(\mathbb C)`$ は `TypeIAlgebra` である — [チュートリアル、2.5節](../../tutorial/ja/02_von_neumann_algebras.md#s-2-5)。ユニタリ $`V`$ と組 $`(k_c,m_c)`$ は `decompose.wedderburn`（フロントエンド `from_group`, `from_generators`）が計算する — [チュートリアル、4.3節](../../tutorial/ja/04_symmetries.md#s-4-3)。

数 $`k_c`$ は環そのものを記述する。**重複度** $`m_c`$ は環が
$`M_n(\mathbb C)`$ の中にどのように置かれているかを記述する：それらは代数的構造からは見えないが、
$`\mathcal M`$ に制限した $`\mathbb C^n`$ のトレースを決定する。
```math
\mathrm{Tr}_{\mathbb C^n}(A)=\sum_cm_c\mathrm{Tr}(A_c),\qquad A=V\Big(\bigoplus_cA_c\otimes1_{m_c}\Big)V^*.
```
これが、ライブラリがトレースの族を扱う理由である（第3章）。同様に、自己共役な $`A`$ に対し、
$`A`$ の $`V(\mathbb C^{k_c}\otimes\mathbb C^{m_c})`$ への制限はブロック $`A_c`$ の固有値をもち、その重複度はそれぞれ
$`A_c`$ における重複度の $`m_c`$ 倍である。$`A_c\otimes1_{m_c}`$ は $`A_c`$ の $`m_c`$ 個のコピーだからである。

*ライブラリでは：* トレース `Tr_blunt`, `Tr_norm`, `tau_vN` — [チュートリアル、3.1節](../../tutorial/ja/03_traces_and_states.md#s-3-1)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. $`P_c`$ は **超選択セクター** あるいは **対称性
> セクター** である：$`\mathcal M`$ のいかなるオブザーバブルも異なるセクターを結びつけないので、
> 「どのセクターか」の測定は他のすべての測定と両立する。中心は系の
> *古典的* な部分である。重複度 $`m_c`$ は $`\mathcal M`$ のすべてのオブザーバブルの縮退度である。

<a id="s-2-4"></a>

## 2.4 例：対称性の群とその可換子環

有限群 $`G`$ がユニタリ表現 $`\pi`$ によって $`\mathbb C^n`$ に作用しているとする。対称性と可換な作用素全体は
$`\pi(G)'`$ をなし、これはフォン・ノイマン環である。その
可換子環 $`\pi(G)''`$ は $`\pi(g)`$ たちの張る空間である（その張る空間に定理 2.2 を適用する。これは
$`*`$-環である：$`\pi(e)=1`$、積 $`\pi(g)\pi(h)=\pi(gh)`$、随伴
$`\pi(g)^*=\pi(g^{-1})`$ を含む）。

<a id="r-2-7"></a>

**命題 2.7.** 定理 2.6 を $`\mathcal M=\pi(G)''`$ に適用する。すると
$`V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}`$ であり、ここで各 $`\rho_c`$ は $`\mathbb C^{k_c}`$ 上の
$`G`$ の既約ユニタリ表現であり、$`c\neq d`$ に対し $`\rho_c,\rho_d`$ は同値でない。したがって
```math
\pi\cong\bigoplus_c\rho_c\otimes1_{m_c},\qquad
\pi(G)''\cong\bigoplus_cM_{\dim\rho_c},\qquad
\pi(G)'\cong\bigoplus_cM_{m_c}:
```
対称な作用素 $`\pi(G)'`$ のブロックの大きさは既約表現の重複度であり、
$`k_c`$ と $`m_c`$ の役割は環とその可換子環の間で入れ替わる。

*証明.* $`\pi(g)\in\mathcal M`$ であるから、ある
$`\rho_c(g)\in M_{k_c}`$ により $`V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}`$ と書け、$`\pi`$ がユニタリ表現であるから $`g\mapsto\rho_c(g)`$ もユニタリ表現である。
$`\pi(g)`$ たちは $`\mathcal M`$ を張るので、$`\rho_c(g)`$ たちは $`M_{k_c}`$ を張り、より強く、写像
$`\sum_ga_g\pi(g)\mapsto\big(\sum_ga_g\rho_c(g)\big)_c`$ は $`\bigoplus_cM_{k_c}`$ の上への全射である。
*既約性：* すべての $`\rho_c(g)`$ で不変な $`\mathbb C^{k_c}`$ の部分空間は、それらの張る $`M_{k_c}`$ で不変であるから、$`0`$ または $`\mathbb C^{k_c}`$ である。*非同値性：* ある可逆な $`T`$ と $`c\neq d`$ に対して、すべての $`g`$ について
$`\rho_d(g)=T\rho_c(g)T^{-1}`$ ならば、像の任意の元
$`(x_{c'})_{c'}`$ は $`x_d=Tx_cT^{-1}`$ を満たし、これは
$`x_c=1`$, $`x_d=0`$ なる元を除外する——全射性に矛盾する。最後の二つの同型は
定理 2.6 である。$`\square`$

*ライブラリでは：* 有限群の表現の分解、`decompose.from_group` — [チュートリアル、4.3節](../../tutorial/ja/04_symmetries.md#s-4-3)。

表現の既約表現への分解は同値を除いて一意であり、したがって $`m_c`$ はいかなる分解においても *その* 重複度である。*証明.* $`\rho,\rho'`$ を既約とし、$`T`$ を絡作用素、すなわちすべての $`g`$ に対して $`T\rho(g)=\rho'(g)T`$ を満たすものとする。$`T`$ の核と像は不変であるから、$`T=0`$ であるか $`T`$ は可逆である（シューアの補題）。$`\rho'=\rho`$ ならば、$`T`$ の固有値 $`\lambda`$ に対して $`T-\lambda1`$ に同じことが当てはまるので、$`T=\lambda1`$ である。よって $`\rho`$ から $`\rho'`$ への絡作用素の空間の次元は、$`\rho\cong\rho'`$ なら $`1`$（可逆な絡作用素 $`S`$ に対し、$`T\mapsto S^{-1}T`$ によってそれらは $`\rho`$ から自身への絡作用素、すなわち $`\mathbb C1`$ と同一視される）、そうでなければ $`0`$ である（零でない絡作用素は可逆、すなわち同値を与えるはずだからである）。既約な $`\sigma_j`$（重複を許す）により $`\pi\cong\bigoplus_j\sigma_j`$ と書けば、既約な $`\rho`$ から $`\bigoplus_j\sigma_j`$ への絡作用素は $`\rho`$ から $`\sigma_j`$ への絡作用素の組 $`(T_j)_j`$ であり、同値写像との合成は次元を変えない。よって $`\rho`$ から $`\pi`$ への絡作用素の空間の次元は $`\#\{j:\sigma_j\cong\rho\}`$ となり、これは $`\pi`$ と $`\rho`$ のみで決まる数である。$`\rho_c`$ が互いに非同値である命題 2.7 の分解では、この数は $`\rho\cong\rho_c`$ なら $`m_c`$ であり、$`\rho`$ がどの $`\rho_c`$ とも同値でなければ $`0`$ である。$`\square`$

$`G=S_3`$ が $`(\mathbb C^2)^{\otimes3}`$ の三つのテンソル因子を置換する場合、すなわち三つの
量子ビットあるいは三つのスピン $`\tfrac12`$ の場合を考える：$`x\in\{0,1\}^3`$ に対し $`\pi(\sigma)|x\rangle=|\sigma\cdot x\rangle`$ であり、ここで
$`(\sigma\cdot x)_{\sigma(i)}=x_i`$ である。$`F=\{v:\pi(g)v=v\ \forall g\}`$ を **対称部分空間** とする。分解は次のとおりである：$`F`$（次元 $`4`$）は重複度 $`4`$ の自明表現を担い、その直交補空間（次元 $`4`$）は重複度 $`2`$ の二次元 **標準表現**（座標の和が $`0`$ である $`\mathbb C^3`$ のベクトル上の $`S_3`$ の置換表現）を担い、符号表現は現れない。したがって
$`\pi(S_3)''\cong M_1\oplus M_2`$（次元 $`5`$）かつ $`\pi(S_3)'\cong M_4\oplus M_2`$
（次元 $`20`$）である。

*証明.* $`V^*\pi(S_3)''V=\bigoplus_cM_{k_c}\otimes1_{m_c}`$ と書く（定理 2.6）。すると
$`\sum_ck_cm_c=8`$ かつ $`\sum_cm_c^2=\dim\pi(S_3)'`$ である。

*ステップ1：$`\dim\pi(S_3)'=20`$.* $`\pi(g)^*|g\cdot x\rangle=|x\rangle`$ であるから、行列
$`\pi(g)X\pi(g)^*`$ は位置 $`(g\cdot x,g\cdot y)`$ に成分 $`X_{xy}`$ をもつ。よって $`X\in\pi(S_3)'`$、すなわちすべての $`g`$ に対し
$`\pi(g)X\pi(g)^*=X`$ であることは、$`(x,y)\mapsto X_{xy}`$ が、$`g\cdot(x,y)=(g\cdot x,g\cdot y)`$ による $`S_3`$ の
対 $`(x,y)`$ への作用の軌道上で一定であることと同値である。軌道の指示行列は
$`\pi(S_3)'`$ の基底をなすから、その次元は軌道の個数である。対 $`(x,y)`$ は四元集合 $`\{0,1\}^2`$ の元の三つ組
$`\big((x_1,y_1),(x_2,y_2),(x_3,y_3)\big)`$ であり、$`S_3`$ は三つの成分を置換するから、軌道は四つの元から作る
大きさ $`3`$ の多重集合である：$`\{a,a,a\}`$ の形のものが $`4`$ 個、$`b\neq a`$ として $`\{a,a,b\}`$ の形のものが $`4\cdot3=12`$ 個、
相異なる元からなる $`\{a,b,c\}`$ の形のものが $`4`$ 個、合計 $`20`$ 個である。

*ステップ2：$`k=(1,2)`$, $`m=(4,2)`$.* $`k_c\ge1`$ であるから $`\sum_cm_c\le\sum_ck_cm_c=8`$ である。$`M=\max_cm_c`$ とおく。
すると $`20=\sum_cm_c^2\le M\sum_cm_c\le8M`$ かつ $`M^2\le20`$ であるから、$`M\in\{3,4\}`$ である。$`M=4`$ ならば、残りの
$`m_c`$ の平方和は $`4`$ であるから、それらは $`(2)`$ または $`(1,1,1,1)`$ である。$`M=3`$ ならば、残りの $`m_c`$ の平方和は $`11`$、
和は高々 $`5`$ である。それらがすべて $`\le2`$ ならば平方和は高々
$`2\cdot5=10`$ となるから、そのうち一つは $`3`$ であり、残りの平方和は $`2`$、すなわち $`(1,1)`$ である。
したがって $`m`$ は順序を除いて $`(4,2)`$, $`(4,1,1,1,1)`$, $`(3,3,1,1)`$ のいずれかである。後の二つの場合、
$`\sum_cm_c=8=\sum_ck_cm_c`$ からすべての $`k_c=1`$ となり、$`\pi(S_3)''\cong\bigoplus_cM_1`$ は
可換になるが、互換 $`(12)`$ と $`(23)`$ は可換でない：
$`\pi((12))\pi((23))|001\rangle=|100\rangle`$ であるが $`\pi((23))\pi((12))|001\rangle=|010\rangle`$ である。よって
$`m=(4,2)`$ であり（$`m_1=4`$ となるように直和成分に番号を付ける）、$`k_1,k_2\ge1`$ のもとで $`4k_1+2k_2=8`$ から $`k=(1,2)`$ を得る。特に
$`\dim\pi(S_3)''=\sum_ck_c^2=5`$ である。

*ステップ3：自明表現と符号表現.* $`\rho_2`$ は零でないベクトル $`u`$ を固定しない：もし固定すれば、直線
$`\mathbb Cu`$ は $`\mathbb C^2`$ の不変部分空間となり、既約性に反する
（命題 2.7）。$`v_1\in\mathbb C^1\otimes\mathbb C^4`$ および
$`v_2=\sum_{r=1}^2u_r\otimes f_r\in\mathbb C^2\otimes\mathbb C^2`$（$`f_1,f_2`$ は標準基底）により $`V^*v=v_1\oplus v_2`$ と書く。すると
$`v\in F`$ であることは、すべての $`g`$ と $`r`$ に対し $`\rho_1(g)v_1=v_1`$ かつ $`\rho_2(g)u_r=u_r`$ であることと同値であり、よって任意の $`v\in F`$ は
$`v_2=0`$ を満たす。零でないベクトル $`|000\rangle`$ は $`F`$ に属するから、スカラー $`\rho_1(g)`$ は零でないベクトルを固定する：
$`\rho_1`$ は自明表現であり、$`F=V(\mathbb C^1\otimes\mathbb C^4)`$ は次元 $`4`$ をもち、その直交
補空間は $`V(\mathbb C^2\otimes\mathbb C^2)`$ である。符号表現は
$`\rho_1`$ とも（$`\mathrm{sgn}(12)=-1`$ であるから）$`\rho_2`$ とも（次元が異なる）同値でないから、上で示した
重複度の一意性により、それは $`0`$ 回現れる。

*ステップ4：$`\rho_2`$ は標準表現である.* ベクトル $`e_1=|100\rangle`$, $`e_2=|010\rangle`$,
$`e_3=|001\rangle`$ は、$`\sigma`$ が位置を置換するのと同じように $`\pi(\sigma)`$ によって置換されるから、
$`W=\{\sum_ix_ie_i:\sum_ix_i=0\}`$ は不変であり、$`\pi`$ はその上に標準表現として作用する。
$`W\perp F`$ である：$`v\in F`$ と $`w=\sum_ix_ie_i\in W`$ に対し、すべての $`g`$ について $`\langle v,w\rangle=\langle\pi(g)v,w\rangle=\langle v,\pi(g)^{-1}w\rangle`$
であるから $`\langle v,w\rangle=\langle v,\frac16\sum_g\pi(g)w\rangle`$ であり、
$`\frac16\sum_g\pi(g)w=\frac13\big(\sum_ix_i\big)(e_1+e_2+e_3)=0`$ である（六つの置換のうちちょうど二つが各 $`e_i`$ を各 $`e_j`$ に写す）。
よって $`V^*W\subseteq\mathbb C^2\otimes\mathbb C^2`$ である。写像
$`T_r:\mathbb C^2\otimes\mathbb C^2\to\mathbb C^2`$, $`u\otimes f\mapsto\langle f_r,f\rangle u`$ は
$`T_r(\rho_2(g)\otimes1)=\rho_2(g)T_r`$ を満たす。$`T_rV^*|_W\neq0`$ となる $`r`$ を選ぶ（両方とも消えるならば $`V^*W`$ は
$`0`$ となってしまう）。すると $`T_rV^*|_W`$ は $`\pi|_W`$ から $`\rho_2`$ への零でない絡作用素であり、その像は $`\mathbb C^2`$ の零でない
不変部分空間、したがって $`\mathbb C^2`$ であり、$`\dim W=2`$ であるから全単射である。よって $`\rho_2`$ は
標準表現と同値である（ついでに標準表現が既約であることもわかる）。$`\square`$

*計算：* 次元 $`20`$ と $`5`$、およびこの例における定理 2.2 — [チュートリアル、2.4節](../../tutorial/ja/02_von_neumann_algebras.md#s-2-4)。*ライブラリでは：* `decompose.permutation_representation`, `decompose.commutant_dimension` — [チュートリアル、4.2節](../../tutorial/ja/04_symmetries.md#s-4-2)。

重複度はスペクトルに現れる。任意の自己共役な $`K\in\pi(S_3)'`$ は、自己共役な $`\kappa_1\in M_4`$, $`\kappa_2\in M_2`$ により
$`V^*KV=(1_1\otimes\kappa_1)\oplus(1_2\otimes\kappa_2)`$ の形をもつ
（定理 2.6）。$`F`$ 上でのその固有値は $`\kappa_1`$ の固有値であり、$`F`$ の補空間上では
$`\kappa_2`$ の固有値がそれぞれ二回ずつ現れる。$`1_2\otimes\kappa_2`$ は $`\kappa_2`$ の二つのコピーだからである。よって $`F`$ の補空間上の
$`K`$ の *すべての* 固有値は偶数の重複度をもつ。ランダムな $`K`$ についてはさらに多くのことが成り立つ。群平均
$`E(X)=\frac1{|G|}\sum_g\pi(g)X\pi(g)^*`$ は $`M_8(\mathbb C)`$ を可換子環の上へ写し、そこでは
恒等写像である（演習問題 3）。よってジニブル行列 $`Z`$（成分が独立な標準複素ガウス変数）に対し、
$`K=E(Z+Z^*)`$ は自己共役元のなす実ベクトル空間 $`\pi(S_3)'_h`$ 上に密度をもつ
（定理 A.15 とその後の注意）。命題 4.7(a)（先の章への参照である。その証明は
定理 2.6 と付録 A を用い、この段落は用いない）により、確率 1 で、$`\kappa_1`$ と $`\kappa_2`$ の六つの
固有値は互いに相異なる。このとき $`K`$ は $`F`$ 上に $`4`$ 個の単純固有値をもち、その補空間上に
$`2`$ 個の二重縮退した固有値をもつ。最後に、$`P_F=\frac1{|G|}\sum_g\pi(g)`$ は $`F`$ への直交射影である：すべての $`h`$ に対し $`\pi(h)P_F=P_F`$ であるから、その
値域は $`F`$ に含まれる。$`v\in F`$ に対し $`P_Fv=v`$ である。また置換 $`g\to g^{-1}`$ により $`P_F^*=P_F`$ である。

*計算：* ランダムな $`K`$ の二重縮退した固有値 — [チュートリアル、2.4節](../../tutorial/ja/02_von_neumann_algebras.md#s-2-4)。

<a id="s-2-5"></a>

## 2.5 ライブラリにおける環

*ライブラリでは：* `TypeIAlgebra` が $`\bigoplus_cM_{k_c}(\mathbb C)`$、その作用素、およびその中心をどのように格納するか — [チュートリアル、2.5節](../../tutorial/ja/02_von_neumann_algebras.md#s-2-5)。

<a id="s-2-6"></a>

## 2.6 条件付き期待値

$`\mathcal M`$ からフォン・ノイマン部分環
$`\mathcal N\subseteq\mathcal M`$（$`\mathcal M`$ に含まれ、同じ単位元をもつ $`*`$-環）の上への **条件付き期待値** とは、線形写像 $`E:\mathcal M\to\mathcal N`$ であって、正
（$`A\ge0\Rightarrow E(A)\ge0`$）であり、$`B\in\mathcal N`$ に対して $`E(B)=B`$ を満たし、かつ
$`\mathcal N`$-双加群写像であるもの：$`B_1,B_2\in\mathcal N`$ に対し $`E(B_1AB_2)=B_1E(A)B_2`$。

<a id="r-2-8"></a>

**命題 2.8.** $`\bigoplus_c\mathbb C^{k_c}`$ 上の $`\mathcal M=\bigoplus_cM_{k_c}`$（$`P_c`$ は第 $`c`$ ブロックへの射影）
および $`\mathcal N=Z(\mathcal M)=\mathrm{span}\{P_1,\dots,P_C\}`$（演習問題 1）に対し、
```math
E(A)=\sum_c\frac{\mathrm{Tr}\,A_c}{k_c}\,P_c
```
は条件付き期待値である。これは
$`\sum_cw_c\mathrm{Tr}\,A_c`$ の形のすべてのトレースを保つ。

*証明.* 線形性は明らかである。$`A\ge0`$ ならば $`\mathrm{Tr}\,A_c=\sum_i\langle e_i,Ae_i\rangle\ge0`$（第 $`c`$ ブロックの標準基底ベクトル $`e_i`$ にわたる和）であるから $`E(A)\ge0`$ である。
$`B=\sum_cb_cP_c`$ に対し $`E(B)=\sum_cb_ck_cP_c/k_c=B`$ である。$`B_i=\sum_cb^{(i)}_cP_c`$ に対し
$`\mathrm{Tr}(B_1AB_2)_c=b^{(1)}_cb^{(2)}_c\mathrm{Tr}\,A_c`$ であり、これから双加群性が従う。最後に、すべての $`c`$ に対し $`\mathrm{Tr}(E(A)_c)=\mathrm{Tr}\,A_c`$ である。$`\square`$

*ライブラリでは：* $`E`$ は `channels.center_expectation` である — [チュートリアル、2.6節](../../tutorial/ja/02_von_neumann_algebras.md#s-2-6)。

2.4節の群平均は $`M_n(\mathbb C)`$ の中での可換子環
$`\pi(G)'`$ の上への条件付き期待値であり、$`\mathrm{Tr}`$ を保つ（演習問題 3）。同様に、完全な行列環
$`M_D(\mathbb C)`$, $`D=\sum_ck_c`$（これはブロック対角行列として $`\mathcal M=\bigoplus_cM_{k_c}`$ を含む）の上で、
写像 $`E_{\mathcal M}(X)=\sum_cP_cXP_c`$ は $`\mathrm{Tr}`$ を保つ $`\mathcal M`$ の上への条件付き期待値である。
*証明.* $`P_cXP_c`$ は第 $`c`$ 対角ブロックの外で消えるから、$`E_{\mathcal M}(X)\in\mathcal M`$ である。$`X\ge0`$ に対し
$`\langle v,P_cXP_cv\rangle=\langle P_cv,XP_cv\rangle\ge0`$ である。元
$`B\in\mathcal M`$ はブロック対角であるから、$`BP_c=P_cB`$ かつ $`B=\sum_cP_cBP_c`$ である。よって $`E_{\mathcal M}(B)=B`$ であり、
$`B_1,B_2\in\mathcal M`$ に対し $`E_{\mathcal M}(B_1XB_2)=\sum_cB_1P_cXP_cB_2=B_1E_{\mathcal M}(X)B_2`$ である。最後に
$`\sum_c\mathrm{Tr}(P_cXP_c)=\sum_c\mathrm{Tr}(P_cX)=\mathrm{Tr}\,X`$ である。$`\square`$

*ライブラリでは：* 群平均は `decompose.twirl` である — [チュートリアル、4.3節](../../tutorial/ja/04_symmetries.md#s-4-3)。

$`\mathcal M`$ 上の **状態** $`\omega`$（すべての $`A`$ に対し $`\omega(A^*A)\ge0`$ かつ $`\omega(1)=1`$ を満たす線形汎関数。
状態は第3章で扱う）に対し、命題 2.8 と線形性から
$`\omega(E(A))=\sum_c\omega(P_c)\mathrm{Tr}(A_c)/k_c`$ を得る。また中心の上では $`E(z)=z`$ であるから、そこで
$`\omega\circ E=\omega`$ である。（条件付き期待値は一般には一意でない：$`\mathbb C1`$ の上へのものは、任意の状態
$`\varphi`$ から $`A\mapsto\varphi(A)1`$ として一つ得られる。$`M_2(\mathbb C)`$ 上では、状態 $`A\mapsto A_{11}`$ と $`A\mapsto A_{22}`$
が二つの異なるものを与える。）これらはすべて同じ考え方の実例である：
ある構造を「見ない」作用素への射影である。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. 写像 $`E_{\mathcal M}`$ は超選択セクター間のコヒーレンスを消す。電荷（セクター）を
> 測定して結果を記録しないときに起こるのがまさにこれである。中心の上への条件付き期待値 $`E`$ はさらに先へ進む：
> 古典的な情報、すなわち各セクターの重み $`\omega(P_c)`$ だけを残し、各セクター内の状態を最大混合状態（正規化されたトレース
> $`\mathrm{Tr}(A_c)/k_c`$）に置き換える。

<a id="s-2-7"></a>

## 2.7 比較の対象

*計算：* 他のソフトウェア（GAP, SageMath, QuSpin, QuTiP）との比較 — [チュートリアル、2.7節](../../tutorial/ja/02_von_neumann_algebras.md#s-2-7)。

<a id="x-2"></a>

## 演習問題

1. $`M_k(\mathbb C)`$ の中心が $`\mathbb C1`$ であること（行列単位
   $`|i\rangle\langle j|`$ を用いよ）、および $`Z(\bigoplus_cM_{k_c})`$ の次元が $`C`$ であることを示せ。
2. $`S`$ を $`\mathbb C^3`$ 上の巡回シフト $`Se_i=e_{i+1\bmod3}`$ とする。
   $`\{S\}'`$ と $`\{S\}''`$ を決定せよ。
3. 有限群のユニタリ表現 $`\pi`$ に対し、
   $`E(X)=\frac1{|G|}\sum_g\pi(g)X\pi(g)^*`$ が
   $`M_n(\mathbb C)`$ から $`\pi(G)'`$ の上への、$`\mathrm{Tr}`$ を保つ条件付き期待値であることを示せ。

<details>
<summary>解答</summary>

1. $`X`$ がすべての $`|i\rangle\langle j|`$ と可換ならば
   $`X|i\rangle\langle j|=|i\rangle\langle j|X`$ である。両辺を $`e_j`$ に適用して比較すると、
   すべての $`i,j`$ に対し $`Xe_i=X_{jj}e_i`$ を得るので、$`X=X_{11}1`$ である。直和においては、元
   $`\bigoplus X_c`$ が中心的であるのはすべての $`X_c`$ が中心的であるときかつそのときに限るから、中心は
   $`\{\bigoplus_cz_c1_{k_c}\}`$ であり、その次元は $`C`$ である。

2. $`S`$ は固有値 $`1,\omega,\omega^2`$（$`\omega=e^{2\pi i/3}`$）をもつユニタリであり、各固有値は単純である：
   ベクトル $`u_j=\sum_i\omega^{-ij}e_i`$（$`j=0,1,2`$）は $`Su_j=\omega^ju_j`$ を満たし、互いに直交する。
   $`S`$ と可換な作用素はその固有空間を保つので、固有基底において
   対角であり、逆にすべての対角行列は $`S`$ と可換である：$`\{S\}'\cong\mathbb C^3`$。これは
   $`1,S,S^2`$ を含み、それらは固有基底において $`\lambda=(1,\omega,\omega^2)`$ として $`\mathrm{diag}(\lambda^0)`$, $`\mathrm{diag}(\lambda)`$,
   $`\mathrm{diag}(\lambda^2)`$ であり、一次独立である（ヴァンデルモンド：固有値は相異なる）。よって
   $`\{S\}'=\mathrm{span}\{1,S,S^2\}`$、すなわち巡回行列全体である。すると $`\{S\}''`$ も
   $`\{S\}'`$ である。なぜなら $`S`$ の固有基底において、すべての
   対角行列と可換な行列は対角だからである（相異なる対角成分をもつ $`D`$ について $`XD`$ と $`DX`$ の成分を比較せよ）。

3. $`E`$ は線形かつ正である（写像 $`X\mapsto UXU^*`$ の和）。$`h\in G`$ に対し
   $`\pi(h)E(X)\pi(h)^*=\frac1{|G|}\sum_g\pi(hg)X\pi(hg)^*=E(X)`$ であるから $`E(X)\in\pi(G)'`$ である。
   $`X\in\pi(G)'`$ ならば各項は $`X`$ に等しいので $`E(X)=X`$ である。$`B_1,B_2\in\pi(G)'`$ に対し
   $`\pi(g)B_1XB_2\pi(g)^*=B_1\pi(g)X\pi(g)^*B_2`$ であり、これから双加群性が従う。
   最後に、すべての $`g`$ に対し $`\mathrm{Tr}(\pi(g)X\pi(g)^*)=\mathrm{Tr}\,X`$ であるから、
   $`\mathrm{Tr}\,E(X)=\mathrm{Tr}\,X`$ である。

*計算：* 演習問題 2 の数値的検算 — [チュートリアル、演習問題](../../tutorial/ja/02_von_neumann_algebras.md#x-2)。

</details>

**次へ：**[第3章](03_traces_and_states.md) では
$`\bigoplus_cM_{k_c}`$ 上のトレースを分類し、状態と密度を導入し、後で用いるエントロピーとトレースの不等式を証明する。
