[English](../05_random_matrices.md) · [Español](../es/05_random_matrices.md) · [Français](../fr/05_random_matrices.md) · [Deutsch](../de/05_random_matrices.md) · [中文](../zh/05_random_matrices.md) · **日本語** · [Русский](../ru/05_random_matrices.md)

# 5. 環の中のランダム行列

ランダム行列理論は「典型的な」作用素がどのような姿をしているかを問う。環 $\bigoplus_cM_{k_c}$ においては、その答えはセクターを尊重しなければならない：自然なランダム・ユニタリは大きな行列環 $M_D$ のものではなく環のものであり、自然なランダム・ハミルトニアンは互いに独立なブロックをもつ。本章では環のユニタリ群上のハール測度を構成し、その1次と2次のモーメントを厳密に計算し、そのうえでセクターを無視するとスペクトル統計で何がおかしくなるかを示す。

**本章で学ぶこと**：ガウス行列の QR 分解がなぜハール・ユニタリを与えるか（そしてなぜ位相を固定しなければならないか）、環のユニタリ群にわたる平均が中心への条件付き期待値であること、2次のワインガルテン公式、ウィグナー推測、間隔比と $3\times3$ ガウス行列に対するその厳密な値、そしてセクターが混ざったスペクトルで計算した準位統計がなぜ独立な準位のそれのように見えるか。

*計算：* 本章の数値的検証 — [チュートリアル 5](../../tutorial/ja/05_random_matrices.md)。

本章を通じて、$\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$、$D=\sum_ck_c$ とし、$P_c$ は
極小中心射影、$E$ は命題 2.8 の中心への条件付き期待値とする。$U(k)$ は $k\times k$ ユニタリ行列のなす群であり、
$U(\mathcal M)=\{U\in\mathcal M:U^*U=1\}$ は $\mathcal M$ のユニタリ群である。確率変数は期待値 $\mathbb E$ をもつ確率空間上に定義され、ランダム行列の期待値は成分ごとにとる。*標準複素ガウス* 変数とは、実部と虚部が独立で
$\sim N(0,\tfrac12)$ であり、したがって $\mathbb E|z|^2=1$ となるものである。*ジニブル行列* とは、成分が独立な標準複素ガウス変数である行列である。

<a id="s-5-1"></a>

## 5.1 ハール測度

$G$ を行列のコンパクト群（ある $U(n)$ の閉部分群）とする。$G$ 上のボレル確率測度
$\mu$ が **左不変** であるとは、すべての $g\in G$ とボレル集合
$B$ に対し $\mu(gB)=\mu(B)$ となること、同値なことに、すべての有界ボレル関数 $f$ に対し $\int f(gh)\,d\mu(h)=\int f\,d\mu$ となることである。
**右不変** は $Bg$ と $f(hg)$ を用いて同様に定義する。

<a id="r-5-1"></a>

**補題 5.1（一意性）.** $\mu$ が $G$ 上の左不変ボレル確率測度、$\nu$ が右不変ボレル確率測度ならば、$\mu=\nu$ である。したがって $G$ 上の左不変ボレル確率測度は高々一つであり、存在すればそれは右不変でもあり、$g\mapsto g^{-1}$ のもとで不変である。

*証明.* $f:G\to\mathbb R$ を連続（したがって有界）とする。フビニの定理により
$\int\!\big(\int f(gh)\,d\mu(h)\big)d\nu(g)=\int\!\big(\int f(gh)\,d\nu(g)\big)d\mu(h)$ である。左辺の内側の
積分は左不変性により $\int f\,d\mu$ であり、右辺の内側の積分は右不変性により $\int f\,d\nu$ である。
よって $\int f\,d\mu=\int f\,d\nu$ である。距離空間上の二つのボレル確率測度は、連続関数の積分がすべて等しければ一致する（系 A.2：連続関数
$\max(0,1-n\mathrm{dist}(\cdot,K))$ は閉集合 $K$ の指示関数に単調減少しながら収束し、閉集合は
ディンキンの補題（補題 A.1）により測度を決定する）。$\mu_1,\mu_2$ が左不変ならば、
$(Bg)^{-1}=g^{-1}B^{-1}$ であるから $\nu(B)=\mu_2(B^{-1})$ は右不変である。したがって $\mu_1=\nu=\mu_2$ であり、
$\mu_1$ は右不変かつ逆元をとる操作で不変である。$\square$

**定義.** $G$ の **ハール測度** とは、その左不変ボレル確率測度のことである。この分布に従う
ランダムな元を **ハール分布に従う** という。

補題 5.1 は存在については何も述べていない。本講義で必要となるのは $U(k)$、$U(1)$、
$U(\mathcal M)$、$SU(k)$、$O(k)$、および補題 C.2 の対角群 $T$ のハール測度だけであり、それぞれ
明示的に構成される：$U(k)$ は定理 5.3、$U(\mathcal M)$ は系 5.4、$SU(k)$ はその直後、
$O(k)$ は定理 B.17、$T$ は補題 C.2 で構成される（B.7節を参照）。

<a id="r-5-2"></a>

**補題 5.2（QR 分解）.** 可逆な $Z\in M_k(\mathbb C)$ はすべて、$Q$ をユニタリ、$R$ を対角成分が正の上三角行列として $Z=QR$ と書ける。$Q$ と $R$ は一意であり、
$Z$ に連続的に依存する。

*証明.* $Z$ の列 $z_1,\dots,z_k$ は線形独立であり、これにグラム・シュミットの直交化を適用する：
$w_j=z_j-\sum_{i<j}\langle q_i,z_j\rangle q_i$ として $q_j=w_j/\|w_j\|$。ここで
$z_j\notin\mathrm{span}\{z_1,\dots,z_{j-1}\}=\mathrm{span}\{q_1,\dots,q_{j-1}\}$ であるから $w_j\neq0$ である。すると
$z_j=\sum_{i<j}\langle q_i,z_j\rangle q_i+\|w_j\|q_j$、すなわち $R_{jj}=\|w_j\|>0$ として $Z=QR$ であり、
各段階は可逆行列の上で連続である。$Q_1R_1=Q_2R_2$ ならば、
$T=Q_2^*Q_1=R_2R_1^{-1}$ はユニタリかつ対角成分が正の上三角行列である。その逆行列
$T^*$ は下三角であり、かつ上三角でもある（上三角行列の逆行列は上三角である）から、$T$ は対角行列である。対角成分が正の対角ユニタリ行列は $1$ である。$\square$

<a id="r-5-3"></a>

**定理 5.3（ガウス行列からのハール・ユニタリ）.** $Z$ を $k\times k$
ジニブル行列とする。このとき $Z$ は確率1で可逆であり、補題 5.2 の因子 $Q$ は
$U(k)$ 上でハール分布に従う。さらに、$Z=Q'R'$ が $Q'$ ユニタリ、$R'$
上三角である任意の QR 分解（例えば `torch.linalg.qr` が返すもの）ならば、
$\Lambda=\mathrm{diag}(R'_{ii}/|R'_{ii}|)$ として $Q=Q'\Lambda$ はハール分布に従う。

*証明.* *可逆性.* $|\det Z|^2$ は $Z$ の $2k^2$ 個の実座標についての実多項式で、恒等的に零ではなく
（$Z=1$ で $1$ に等しい）、$Z$ は密度 $\pi^{-k^2}e^{-\|Z\|_2^2}$ をもつ
（命題 A.16(b)）。よって系 A.7(b) により $P(\det Z=0)=0$ である。

*ガウス分布の不変性.* $V\in U(k)$ に対し、写像 $Z\mapsto VZ$ は実線形で
$\|Z\|_2$ を保つ（1.4節）から、$M_k(\mathbb C)\cong\mathbb R^{2k^2}$ の直交写像であり、
ルベーグ測度を保ち（補題 A.3(b)）、密度も保つ：$VZ$ は $Z$ と同じ分布をもつ
（命題 A.16(b)）。

*同変性.* $Z=QR$ ならば、$VZ=(VQ)R$ は $VZ$ の補題 5.2 の分解であるから、
$Q(VZ)=VQ(Z)$ である。$Q$ は可逆行列の上で連続である（補題 5.2）から $Q(Z)$ は $U(k)$ のランダムな
元であり、任意のボレル集合 $B\subseteq U(k)$ に対し
$P(VQ(Z)\in B)=P(Q(VZ)\in B)=P(Q(Z)\in B)$ となる：$Q(Z)$ の分布は左不変であり、
補題 5.1 によりそれはハール測度である（これは $U(k)$ のハール測度の存在の証明にもなっている）。
最後に $Q'R'=(Q'\Lambda)(\Lambda^*R')$ であり、$\Lambda^*R'$ は対角成分
$|R'_{ii}|>0$ の上三角行列であるから、一意性により $Q'\Lambda=Q(Z)$ である。$\square$

実ガウス行列に対して同じ証明を行えば、補正 $\mathrm{diag}(\mathrm{sgn}\,R'_{ii})$ のもとで
直交群 $O(k)$ のハール測度が得られる。これは定理 B.17 で詳しく述べる。
補正 $\Lambda$ がなければ、$Q'$ の分布は QR ルーチンの符号の規約に依存し、ハール測度であるとは限らない（$\Lambda$ を用いる方法は F. Mezzadri, *How to
generate random matrices from the classical compact groups*, Notices AMS 54 (2007), 592–604 による）。

*計算：* ライブラリのサンプラー `TypeIAlgebra.random_unitary` は定理 5.3 の補正付き QR であり、`torch.linalg.qr` の補正なしの $Q'$ はハール分布に従わない — [チュートリアル、5.1節](../../tutorial/ja/05_random_matrices.md#s-5-1)。

<a id="r-5-4"></a>

**系 5.4（環のハール測度）.** $U(\mathcal M)=\{\bigoplus_cU_c:U_c\in U(k_c)\}\cong\prod_cU(k_c)$
であり、$U=\bigoplus_cU_c$ が $U(\mathcal M)$ 上でハール分布に従うための必要十分条件は、ブロック $U_c$ が独立であり、
各 $U_c$ が $U(k_c)$ 上でハール分布に従うことである。

*証明.* $U^*U=\bigoplus_cU_c^*U_c$ であるから、$U$ がユニタリであることと各ブロックがユニタリであることは同値である。$U_c$ が
独立でハール分布に従うとする。$V=\bigoplus_cV_c\in U(\mathcal M)$ に対し、$VU$ のブロック $V_cU_c$ は
独立であり（独立な変数の関数）、$V_cU_c$ は $U_c$ と同じ分布をもつから、
$VU$ は $U$ と同じ分布をもつ。したがって $U$ の分布は左不変であり、補題 5.1 によりそれはコンパクト群 $U(\mathcal M)$ のハール
測度である。逆に、一意性により、ハール分布に従う任意の $U$ はこの分布をもつ。$\square$

*ライブラリでは：* 独立なブロックをもつ $U(\mathcal M)$ のハール・ユニタリは `TypeIAlgebra.random_unitary_operator` で生成する — [チュートリアル、5.1節](../../tutorial/ja/05_random_matrices.md#s-5-1)。

群 $U(\mathcal M)$ は $U(D)$ よりはるかに小さい：$C\ge2$ のとき、$U(D)$ のハール・ユニタリ $U$ が
$\mathcal M$ に属する確率は零である。実際、定理 5.3 と補題 5.1 により、$U$ は $D\times D$ ジニブル行列 $Z$ の
因子 $Q$ と同じ分布をもち、補題 5.2 のグラム・シュミットによる構成により、$Q$ の第1列は $Z$ の第1列 $z$ を用いて
$z/\|z\|$ と表される。$z$ は密度をもつ $\mathbb C^D$ の標準複素ガウス・ベクトルである（命題 A.16(a)）。各 $i$ に対し
$\{z_i=0\}$ は $\mathbb C^D$ の真部分空間であるから $P(z_i=0)=0$ であり（系 A.7(c)）、第1列の成分はほとんど確実に
すべて非零である。$\mathcal M$ の元の第1列は、最初のブロックの行以外では成分が零であり、$C\ge2$ のときそのような
行が存在する。

$SU(k)$ のハール測度も同じ考え方で得られる：$U$ が $U(k)$ 上でハール分布に従い $\theta=\arg\det U\in(-\pi,\pi]$（$U$ のボレル関数）ならば $V=e^{-i\theta/k}U\in SU(k)$ であり、$W\in SU(k)$ に対し
$\det(WU)=\det U$ から $V(WU)=WV(U)$ を得る。$WU$ は $U$ と同じ分布をもつから、$V$ の分布は
$SU(k)$ 上で左不変であり、補題 5.1 により $SU(k)$ のハール測度である。

<a id="r-5-5"></a>

**系 5.5（成分の低次モーメント）.** $U$ が $U(k)$ 上でハール分布に従うならば、
$\mathbb EU_{ij}=0$ かつ $\mathbb E|U_{ij}|^2=1/k$ である。

*証明.* $-1\in U(k)$ であるから、$-U$ は $U$ と同じ分布をもち、$\mathbb EU=-\mathbb EU$ である。置換
行列 $P$ に対し $PU$ は $U$ と同じ分布をもつから、$\mathbb E|U_{ij}|^2$ は $i$ によらない。
$\sum_i|U_{ij}|^2=1$ であるから、それは $1/k$ に等しい。$\square$

*計算：* 系 5.4 と 5.5 のモンテカルロ法による検証、および $SU(k)$ のハール測度（`measure='haar_su'`）— [チュートリアル、5.1節](../../tutorial/ja/05_random_matrices.md#s-5-1)。

<a id="s-5-2"></a>

## 5.2 モーメントとワインガルテン計算

<a id="r-5-6"></a>

**補題 5.6.** $T\in M_k(\mathbb C)$ がすべての対角ユニタリ行列およびすべての置換
行列と可換ならば、$T\in\mathbb C1$ である。特に $U(k)$ の可換子環は $\mathbb C1$ である。

*証明.* 相異なる $e^{i\theta_j}$ をもつ $\Delta=\mathrm{diag}(e^{i\theta_1},\dots,e^{i\theta_k})$ に対し、
$\Delta T=T\Delta$ から $(e^{i\theta_i}-e^{i\theta_j})T_{ij}=0$ を得るので、$T$ は
対角行列である。置換行列 $P_\sigma e_j=e_{\sigma(j)}$ に対し、$P_\sigma TP_\sigma^{-1}=T$ から
$T_{\sigma(j)\sigma(j)}=T_{jj}$ を得るので、対角成分はすべて等しい。$\square$

<a id="r-5-7"></a>

**命題 5.7（ハール・ツワールは中心への条件付き期待値である）.** $U$ を
$U(\mathcal M)$ 上でハール分布に従うとする。すべての $X\in\mathcal M$ に対し
$$
\mathbb E\,[UXU^*]=E(X)=\bigoplus_c\frac{\mathrm{Tr}\,X_c}{k_c}\,1_{k_c}.
$$

*証明.* $(UXU^*)_c=U_cX_cU_c^*$ であり、$U_c$ は $U(k_c)$ 上でハール分布に従う（系 5.4）。
$T_c=\mathbb E[U_cX_cU_c^*]$ とおく（成分は有界である）。$V\in U(k_c)$ に対し $VU_c$ は
$U_c$ と同じ分布をもつから、$VT_cV^*=\mathbb E[(VU_c)X_c(VU_c)^*]=T_c$ である。補題 5.6 により $T_c=\lambda1$ であり、
$k_c\lambda=\mathrm{Tr}\,T_c=\mathbb E\mathrm{Tr}(U_cX_cU_c^*)=\mathrm{Tr}\,X_c$ である。$\square$

正規化は命題 2.8 および `channels.center_expectation` のものと同じである。
$X\in\mathcal M$ が $U(\mathcal M)$ と可換ならば $UXU^*=X$ であるから、$X=E(X)$ は中心に属する。ハール
平均は第2章の演習問題 3 の群平均の連続版である。代わりに $U$ が
$U(D)$ 上でハール分布に従うならば、同じ証明（$C=1$）により $\mathbb E[UXU^*]=\frac{\mathrm{Tr}\,X}{D}1$ となる：
セクターの重みは失われ、$UXU^*$ は $\mathcal M$ に属するとさえ限らない（$C\ge2$ のとき、$X=P_1$ とし、$U$ を
$e_1$ と $e_{k_1+1}$ の張る平面内の角 $\pi/4$ の回転とすると、$UP_1U^*$ の $(1,k_1+1)$ 成分は $\frac12$ である）。

*計算：* 命題 5.7 のモンテカルロ法による検証と `channels.center_expectation` との照合 — [チュートリアル、5.2節](../../tutorial/ja/05_random_matrices.md#s-5-2)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. 環のハール・ランダムなユニタリは、超選択則を尊重する
> ランダムな時間発展である。それについて平均をとると、各セクター内部のすべてがかき混ぜられ、
> ちょうど古典的部分、すなわち各電荷の重みだけが残る。これが量子
> 情報理論における「ツワール」であり、それは第2章の中心への条件付き期待値である。

2次のモーメントには、$\mathbb C^k\otimes\mathbb C^k$ 上の $\{V\otimes V\}$ の可換子環が必要である。
$F$ を **スワップ**、$F(x\otimes y)=y\otimes x$ とする。これは自己共役かつユニタリであり、
$\mathrm{Tr}\,1=k^2$、$\mathrm{Tr}\,F=\sum_{i,j}\langle e_i\otimes e_j,e_j\otimes e_i\rangle=k$、$F^2=1$ である。

<a id="r-5-8"></a>

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

*計算：* サンプルした三つの $V\otimes V$ の可換子環を用いた、補題 5.8 の線形代数による検証 — [チュートリアル、5.2節](../../tutorial/ja/05_random_matrices.md#s-5-2)。

<a id="r-5-9"></a>

**系 5.9（2次のシューア・ワイル双対性）.** $k\ge2$ とし、$P_\pm=\frac12(1\pm F)$ とおく。
これらは対称部分空間と反対称部分空間
$\mathrm{Sym}=\{x:Fx=x\}$、$\mathrm{Alt}=\{x:Fx=-x\}$ の上への射影であり、次元はそれぞれ $k(k+1)/2$ と $k(k-1)/2$ である。
$U(k)$ の表現 $V\mapsto V\otimes V$ は $\mathrm{Sym}$ 上でも $\mathrm{Alt}$ 上でも既約であり、
二つの表現は同値でない。

*証明.* $F=F^*$ と $F^2=1$ により、$P_\pm$ は $P_++P_-=1$ を満たす射影であり、値域は $\mathrm{Sym}$、
$\mathrm{Alt}$ である。その次元は $\mathrm{Tr}\,P_\pm=(k^2\pm k)/2$ である。
$\{V\otimes V\}$ の張る空間 $\mathcal A$ は $*$-環であり（$1\otimes1$、積 $VW\otimes VW$、随伴 $V^*\otimes V^*$）、
補題 5.8 によりその可換子環は $\mathcal A'=\mathrm{span}\{1,F\}=\mathrm{span}\{P_+,P_-\}$ である。
$K\subseteq\mathrm{Sym}$ がすべての $V\otimes V$ で不変であるとする。補題 2.1 によりその射影は $\mathcal A'$ に属するから、
$x,y\in\{0,1\}$ により $xP_++yP_-$ と書ける（射影だからである）。それは $\mathrm{Alt}$ 上で零であるから $y=0$ であり、
$K=0$ または $K=\mathrm{Sym}$ である。同じ議論が $\mathrm{Alt}$ にも適用できる。
$T:\mathrm{Sym}\to\mathrm{Alt}$ が二つの表現の間の絡作用素ならば、$TP_+$（$\mathbb C^k\otimes\mathbb C^k$ 上の
作用素とみなす）は $\mathcal A$ と可換であるから、$TP_+=xP_++yP_-$ である。それは
$\mathrm{Alt}$ 上で零であるから $y=0$ であり、その値域は $\mathrm{Alt}\cap\mathrm{Sym}=0$ に含まれるから $T=0$ である。$\square$

<a id="r-5-10"></a>

**定理 5.10（2次のワインガルテン公式）.** $k\ge2$ とし、$U$ を $U(k)$ 上でハール分布に従うとし、
$X$ を $\mathbb C^k\otimes\mathbb C^k$ 上の行列とする。このとき
$$
\mathbb E\,[(U\otimes U)X(U\otimes U)^*]=\alpha\,1+\beta F,\qquad
\alpha=\frac{\mathrm{Tr}\,X-\mathrm{Tr}(FX)/k}{k^2-1},\quad
\beta=\frac{\mathrm{Tr}(FX)-\mathrm{Tr}\,X/k}{k^2-1}.
$$
同値なことに、群 $S_2=\{e,\tau\}$ 上で $\mathrm{Wg}(e)=\frac1{k^2-1}$、$\mathrm{Wg}(\tau)=-\frac1{k(k^2-1)}$ とすると
$$
\mathbb E\big[U_{i_1j_1}U_{i_2j_2}\overline{U_{i'_1j'_1}U_{i'_2j'_2}}\big]
=\sum_{\sigma,\pi\in S_2}\delta_{i_1i'_{\sigma(1)}}\delta_{i_2i'_{\sigma(2)}}\,
\delta_{j_1j'_{\pi(1)}}\delta_{j_2j'_{\pi(2)}}\,\mathrm{Wg}(\sigma\pi^{-1}).
$$

*証明.* $T=\mathbb E[(U\otimes U)X(U\otimes U)^*]$ とおく。命題 5.7 と同様に、$T$ は
すべての $V\otimes V$ と可換であるから、補題 5.8 により $T=\alpha1+\beta F$ である。トレースをとると
$\mathrm{Tr}\,T=\mathrm{Tr}\,X$ であり、また $F$ は $U\otimes U$ と可換であるから
$\mathrm{Tr}(FT)=\mathbb E\mathrm{Tr}((U\otimes U)^*F(U\otimes U)X)=\mathrm{Tr}(FX)$ である。
これより $k^2\alpha+k\beta=\mathrm{Tr}\,X$ と $k\alpha+k^2\beta=\mathrm{Tr}(FX)$ を得る。これは
行列式 $k^2(k^2-1)\neq0$ の連立一次方程式であり、その解は述べたとおりである。
成分については $X=|e_{j_1}\otimes e_{j_2}\rangle\langle e_{j'_1}\otimes e_{j'_2}|$ ととる：すると
$\langle e_{i_1}\otimes e_{i_2},T\,e_{i'_1}\otimes e_{i'_2}\rangle$ が左辺であり、$\mathrm{Tr}\,X=\delta_{j_1j'_1}\delta_{j_2j'_2}$
（$\pi=e$ の項）、$\mathrm{Tr}(FX)=\delta_{j_1j'_2}\delta_{j_2j'_1}$（$\pi=\tau$ の項）であり、
$1$、$F$ はそれぞれ $\delta_{i_1i'_1}\delta_{i_2i'_2}$（$\sigma=e$ の項）と
$\delta_{i_1i'_2}\delta_{i_2i'_1}$（$\sigma=\tau$）を寄与する。項をまとめると
$\alpha=\sum_\pi\mathrm{Wg}(e\pi^{-1})[\pi]$、$\beta=\sum_\pi\mathrm{Wg}(\tau\pi^{-1})[\pi]$ となる。
ここで $[\pi]$ は $\pi$ に対応する $j$ のデルタの積を表す。$\square$

行列 $\big(\mathrm{Wg}(\sigma^{-1}\pi)\big)_{\sigma,\pi}$ は、
$\langle A,B\rangle=\mathrm{Tr}(A^*B)$ に関する $1$ と $F$ のグラム行列
$\begin{pmatrix}k^2&k\\k&k^2\end{pmatrix}$ の逆行列である。一般化されるのはこの形である。

*計算：* $k=3$ に対する定理 5.10 のモンテカルロ法による検証 — [チュートリアル、5.2節](../../tutorial/ja/05_random_matrices.md#s-5-2)。

<a id="r-5-11"></a>

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

*計算：* $k=3$ に対する系 5.11 の三つのモーメントのモンテカルロ法による検証 — [チュートリアル、5.2節](../../tutorial/ja/05_random_matrices.md#s-5-2)。

最後のモーメントは負である：ハール・ユニタリの成分は独立ではなく、
相関は相対的に $1/k$ のオーダーである。環については、異なるブロックの成分は
独立であるから（系 5.4）、セクターをまたぐ混合モーメントは単一ブロックのモーメントの積に分解する。

任意の次数のモーメントについても、
$\{V^{\otimes n}\}$ の可換子環がわかれば同じ推論が適用できる。それは $n$ 個のテンソル因子の置換によって張られる
（シューア・ワイル双対性、定理 B.6(c)）。付録 B がこれを実行する。

<a id="r-5-12"></a>

**定理 5.12（任意の次数のワインガルテン公式）.** $n\ge1$ とし、$U$ を $U(k)$ 上でハール分布に従うとする。このとき
$$
\mathbb E\big[U_{i_1j_1}\cdots U_{i_nj_n}\overline{U_{i'_1j'_1}\cdots U_{i'_nj'_n}}\big]
=\sum_{\sigma,\pi\in S_n}\prod_{m=1}^n\delta_{i_mi'_{\sigma(m)}}\delta_{j_mj'_{\pi(m)}}\;\mathrm{Wg}(\sigma\pi^{-1},k),
$$
ここで $\mathrm{Wg}(\cdot,k)$ は $S_n$ 上の関数で、巡回型のみに依存する。$n\le k$ のとき
行列 $\big(\mathrm{Wg}(\sigma^{-1}\pi,k)\big)_{\sigma,\pi}$ はグラム
行列 $\big(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)}\big)_{\sigma,\pi}$ の逆行列である。$n>k$ のときこの行列は特異であり、
その擬似逆行列をとる。$U$ と $\overline U$ の個数が異なるモーメントは零である。

*証明.* これは $d=k$ とした定理 B.13 である：(a) は $\mathrm{Wg}$ を巡回型の関数として $W_{\sigma\pi}=\mathrm{Wg}(\sigma^{-1}\pi)$ を与え、
(b) はグラム行列が可逆であることと $n\le k$ が同値であることを、(d) はモーメントの公式を、(e) は
釣り合わないモーメントが零であることを与える。$\square$

（この公式は D. Weingarten, J. Math. Phys. 19 (1978), 999–1001 と B. Collins, Int. Math. Res.
Not. 2003, 953–982 による。B. Collins, P. Śniady, Comm. Math. Phys. 264 (2006), 773–795 も参照。）

*計算：* $k=3$ に対する $3$ 次のワインガルテン行列と $\mathbb E|U_{11}|^6=\frac1{10}$（演習問題 2）— [チュートリアル、演習問題](../../tutorial/ja/05_random_matrices.md#x-5)。

<a id="s-5-3"></a>

## 5.3 準位統計

$H$ を固有値
$\lambda_1\le\dots\le\lambda_n$ をもつランダムな自己共役 $n\times n$ 行列とする。その **準位間隔** は $s_i=\lambda_{i+1}-\lambda_i$ である。標準的なアンサンブルが
二つある：**GOE** は $H=(G+G^T)/2$ で、$G$ は成分が独立な $N(0,1)$ の実行列である。
**GUE** は $H=(G+G^*)/\sqrt2$ で、$G$ はジニブル行列である（1.7節）。いずれにおいても、対角成分は
$N(0,1)$ に従い、対角線より上の成分は GOE では $N(0,\tfrac12)$ に従い、GUE では実部と虚部が独立にそれぞれ $N(0,\tfrac12)$ に従い、
成分 $h_{ij}$（$i\le j$）は独立であり、
$H$ の密度は $e^{-\mathrm{Tr}\,H^2/2}$ に比例する（命題 A.17。比例関係は
$\mathrm{Tr}\,H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2$ から従う）。

<a id="r-5-13"></a>

**命題 5.13（ウィグナー推測）.** $n=2$、$s=\lambda_2-\lambda_1$、$S=s/\mathbb Es$ とする。
GOE（$\beta=1$）と GUE（$\beta=2$）に対し、$S$ の密度はそれぞれ
$$
p_1(S)=\frac\pi2\,S\,e^{-\pi S^2/4},\qquad p_2(S)=\frac{32}{\pi^2}\,S^2\,e^{-4S^2/\pi},
$$
であり、分布関数はそれぞれ $P(S\le t)=1-e^{-\pi t^2/4}$ と
$P(S\le t)=\mathrm{erf}(2t/\sqrt\pi)-\frac{4t}\pi e^{-4t^2/\pi}$ である。特に
$\varepsilon\to0$ のとき $P(S\le\varepsilon)$ は $\varepsilon^{\beta+1}$ のオーダーである。

*証明.* $H=\begin{pmatrix}a&b\\\bar b&c\end{pmatrix}$ に対し固有値は
$\frac{a+c}2\pm\frac12\sqrt{(a-c)^2+4|b|^2}$ であるから、$s=\sqrt{(a-c)^2+4|b|^2}$ である。命題 A.17 により、$a$、$c$ と
$b$ の実座標は独立で、$a,c\sim N(0,1)$ であり、$b$ の各実座標は $N(0,\frac12)$ に従う。よって
$a-c\sim N(0,2)$（命題 A.14(c)）であり、$2\mathrm{Re}\,b$、$2\mathrm{Im}\,b$ は $N(0,2)$ に従い、これらはすべて独立である。
*GOE:* $x=a-c$ と $y=2b$ は独立な $N(0,2)$ 変数であり、極座標で（系 A.19、
$n=2$、$s_2=2\pi$）
$P(s\le t)=\int_{x^2+y^2\le t^2}\frac{e^{-(x^2+y^2)/4}}{4\pi}\,dx\,dy=\int_0^t\frac r2e^{-r^2/4}dr=1-e^{-t^2/4}$ となる。
すると $\mathbb Es=\int_0^\infty P(s>t)\,dt=\int_0^\infty e^{-t^2/4}dt=\sqrt\pi$ であり、
$P(S\le t)=P(s\le\sqrt\pi t)=1-e^{-\pi t^2/4}$ である。これを微分すればよい。*GUE:* $x=a-c$、$y=2\mathrm{Re}\,b$、
$y'=2\mathrm{Im}\,b$ は独立な $N(0,2)$ 変数であり、定理 A.20（$n=3$、$s_3=4\pi$）により
$s=\sqrt{x^2+y^2+y'^2}$ の密度は $4\pi r^2(4\pi)^{-3/2}e^{-r^2/4}=\frac{r^2e^{-r^2/4}}{2\sqrt\pi}$ である。その平均は
$\frac1{2\sqrt\pi}\int_0^\infty r^3e^{-r^2/4}dr=\frac{8}{2\sqrt\pi}=\frac4{\sqrt\pi}$ である。
$r=4S/\sqrt\pi$ を代入すると $p_2$ が得られ、$p_2$ を部分積分すると
（$a=4/\pi$ として $\int x^2e^{-ax^2}dx=-\frac{x}{2a}e^{-ax^2}+\frac1{2a}\int e^{-ax^2}dx$）
分布関数が得られる。小さい $\varepsilon$ での振る舞いは $p_\beta(S)\sim cS^\beta$ から従う。$\square$

*計算：* $2\times2$ の GOE と GUE に対する命題 5.13 のモンテカルロ法による検証 — [チュートリアル、5.3節](../../tutorial/ja/05_random_matrices.md#s-5-3)。

指数は次元を数えている：実対称 $2\times2$ 行列が縮退するのは二つの実
条件（$a=c$、$b=0$）が成り立つときかつそのときに限り、複素エルミート行列では三つの条件が成り立つときかつそのときに限る。任意の $n$ に対し、GOE（$\beta=1$）と
GUE（$\beta=2$）の順序づけた固有値 $\lambda_1<\dots<\lambda_n$ は同時密度
$$
Z_{n,\beta}^{-1}\prod_{i<j}(\lambda_j-\lambda_i)^\beta\,e^{-\sum_i\lambda_i^2/2}
$$
をもち（定理 C.7。付録 C でワイル積分公式を用いて証明される）、最小間隔 $\delta$ は
$0<\varepsilon\le1$ に対し $c\,\varepsilon^{\beta+1}\le P(\delta<\varepsilon)\le C\varepsilon^{\beta+1}$ を満たす。ここで定数
$0<c\le C$ は $n$ と $\beta$ に依存する（命題 C.9）。極限 $n\to\infty$ は本講義では扱わない
（M. L. Mehta, *Random Matrices*, 3rd ed., Elsevier, 2004 を参照）。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. ウィグナーは、
> 複雑な量子系（重い原子核）のエネルギー準位は統計的に、適切な対称性をもつランダムなハミルトニアンの固有値のようにふるまうと提唱した：
> 時間反転対称性があれば実対称（$\beta=1$）、なければ複素エルミート
> （$\beta=2$）である。このとき準位は *反発* する：同じ対称性をもつ二つの準位はほとんど決して近づかない
> （フォン・ノイマンとウィグナー、1929年）。反発は量子カオスの典型的な特徴であり、独立な
> 準位（反発なし）は可積分系または局在系の特徴である。これらは物理学におけるボヒガス・ジャンノーニ・シュミット
> 予想とベリー・テイバー予想であり、本講義では証明も使用もしない。

間隔は局所的な準位密度に依存し、それはスペクトルに沿って変化する。間隔を
$p_\beta$ と比較するには **アンフォールディング**、すなわち計数関数の滑らかな近似 $\bar N$ による変数変換 $\lambda\mapsto\bar N(\lambda)$ が必要であり、
これにより平均間隔が $1$ になる。結果は $\bar N$ の選び方に依存する。
連続する間隔の比はこの選択を回避する。

**定義.** **間隔比** を $r_i=\min(s_i,s_{i+1})/\max(s_i,s_{i+1})\in[0,1]$ と定める
（V. Oganesyan, D. A. Huse, Phys. Rev. B 75 (2007), 155111）。これは
$\lambda\mapsto a\lambda+b$（$a>0$）で不変であるから、局所的な平均間隔に依存しない。これが
アンフォールディングなしで用いられる理由である。

*ライブラリでは：* 準位の各行の平均間隔比は `SpinChain.level_spacing_ratio` で得られる — [チュートリアル、5.3節](../../tutorial/ja/05_random_matrices.md#s-5-3)。

比には三つの準位が必要であるから、間隔比をもつ最小のガウス行列は $3\times3$ である。それらについては、
$2\times2$ に対する推測と同様に、平均間隔比を厳密に計算できる。

<a id="r-5-14"></a>

**命題 5.14（平均間隔比）.** (a) *独立な準位.* $s,s'$ を独立で同じ指数
分布に従うとする（独立な指数変数 $e_j$ から作った準位 $e_0+\dots+e_{k-1}$ の連続する二つの間隔のように。A.7節）。このとき
$r=\min(s,s')/\max(s,s')$ は $[0,1]$ 上で $P(r\le t)=\frac{2t}{1+t}$、密度 $\frac2{(1+t)^2}$ をもち、
$$
\mathbb Er=2\ln2-1\approx0.3863 .
$$
(b) *$3\times3$ の GOE と GUE.* $H$ を $3\times3$ の GOE（$\beta=1$）または GUE（$\beta=2$）行列とし、その固有値を
$\lambda_1<\lambda_2<\lambda_3$、間隔を $s_1=\lambda_2-\lambda_1$、$s_2=\lambda_3-\lambda_2$ とし、$r=\min(s_1,s_2)/\max(s_1,s_2)$ とする。
このとき $x=s_2/s_1$ は
$$
f_\beta(x)=\frac{(x+x^2)^\beta}{(1+x+x^2)^{1+3\beta/2}}\qquad\text{on }(0,\infty),
$$
に比例する密度をもち、$r$ は $(0,1)$ 上で密度 $f_\beta/I_0$ をもつ。ここで $I_0=\int_0^1f_\beta$ である。さらに
$$
\mathbb Er=4-2\sqrt3\approx0.5359\quad(\text{GOE}),\qquad
\mathbb Er=\frac{2\sqrt3}\pi-\frac12\approx0.6027\quad(\text{GUE}).
$$

*証明.* (a) スケーリングにより、レートは $1$ としてよい。$x=s/s'$ に対し、独立性とトネリの定理により
$P(x\le t)=\int_0^\infty e^{-y}(1-e^{-ty})\,dy=1-\frac1{1+t}=\frac t{1+t}$ である。
$t<1$ に対し事象 $\{x\le t\}$ と $\{1/x\le t\}$ は互いに素であり、$x$ と $1/x$ は同じ分布をもち
（$s,s'$ は交換可能である）、$\{r\le t\}$ はそれらの和集合であるから、$P(r\le t)=\frac{2t}{1+t}$ である。よって
$\mathbb Er=\int_0^1\frac{2t}{(1+t)^2}dt=2\Big[\ln(1+t)+\frac1{1+t}\Big]_0^1=2\ln2-1$ である。

(b) *$x$ の分布.* 定理 C.7 により、ベクトル $\lambda=(\lambda_1,\lambda_2,\lambda_3)$ は
$C=\{\lambda_1<\lambda_2<\lambda_3\}$ 上で密度 $Z^{-1}\Delta(\lambda)^\beta e^{-|\lambda|^2/2}$ をもつ。ここで $\Delta(\lambda)=\prod_{i<j}(\lambda_j-\lambda_i)$ である。線形写像
$$
(m,s_1,s_2)\mapsto\Big(m-\tfrac{2s_1+s_2}3,\ m+\tfrac{s_1-s_2}3,\ m+\tfrac{s_1+2s_2}3\Big)
$$
は行列式 $1$ をもち、$\mathbb R\times(0,\infty)^2$ を $C$ の上に写す。その逆は $m=\frac13\sum_i\lambda_i$、
$s_1=\lambda_2-\lambda_1$、$s_2=\lambda_3-\lambda_2$ である。これらの変数では $\Delta=s_1s_2(s_1+s_2)$、
$|\lambda|^2=3m^2+\frac23(s_1^2+s_1s_2+s_2^2)$ である（展開すればよい）。補題 A.3(c) とフビニの定理（$m$ について積分する）により、
$(s_1,s_2)$ は $(0,\infty)^2$ 上で、ある定数 $c>0$ に対し密度 $c\,\big(s_1s_2(s_1+s_2)\big)^\beta e^{-(s_1^2+s_1s_2+s_2^2)/3}$ をもつ。
$q(x)=1+x+x^2$ とおく。ボレル集合 $A\subseteq(0,\infty)$ に対し、内側の積分で $s_2=xs_1$ と置換し、
トネリの定理を用いると
$$
P(x\in A)=c\int_0^\infty\!\!\int_A s_1^{3\beta+1}(x+x^2)^\beta e^{-s_1^2q(x)/3}\,dx\,ds_1
=c\,\gamma_\beta\int_A(x+x^2)^\beta\Big(\frac3{q(x)}\Big)^{1+3\beta/2}dx ,
$$
となる。ここで因子は置換 $s_1=(3/q(x))^{1/2}u$ から生じ、$\gamma_\beta=\int_0^\infty u^{3\beta+1}e^{-u^2}du\in(0,\infty)$ である。
よって $x$ は密度 $f_\beta/\int_0^\infty f_\beta$ をもつ。

*$[0,1]$ への帰着.* $f_\beta(1/x)=x^2f_\beta(x)$ である（分子と分母に $x$ のべきを掛ければよい）。
$x=1/y$ とすると、ボレル集合 $A\subseteq(0,1)$ に対し $\int_{1/A}f_\beta(x)\,dx=\int_Af_\beta(1/y)y^{-2}dy=\int_Af_\beta$ である。ここで
$1/A=\{1/y:y\in A\}$ である。特に、$I_j=\int_0^1x^jf_\beta(x)\,dx$ として $\int_0^\infty f_\beta=2I_0$ である。$r=\min(x,1/x)$ かつ
$P(x=1)=0$ であるから $\{r\in A\}=\{x\in A\}\cup\{x\in1/A\}$（互いに素）であり、したがって
$$
P(r\in A)=\frac{2\int_Af_\beta}{2I_0}=\int_A\frac{f_\beta}{I_0},\qquad \mathbb Er=\frac{I_1}{I_0} .
$$
$x+x^2=q-1$ であるから $f_\beta=g_\beta(q)$ であり、ここで $g_1(u)=u^{-3/2}-u^{-5/2}$、$g_2(u)=(u-1)^2u^{-4}$ である。
$q'=2x+1$ を用いて $x=\frac12(q'(x)-1)$ と書き、$u=q(x)$ と置換する（$q(0)=1$ から $q(1)=3$ まで）：
$$
I_1=\frac12\int_1^3g_\beta(u)\,du-\frac12I_0,\qquad \mathbb Er=\frac{\int_1^3g_\beta}{2I_0}-\frac12 .
$$
ここで $\int_1^3g_1=\big[-2u^{-1/2}+\frac23u^{-3/2}\big]_1^3=\frac43-\frac{16}{9\sqrt3}$、
$\int_1^3g_2=\big[-u^{-1}+u^{-2}-\frac13u^{-3}\big]_1^3=\frac8{81}$ である。

*積分 $I_0$.* $L_a=\int_0^1q^{-a}dx$ とおく。$(2x+1)^2=4q-3$ から
$\frac d{dx}\big[(2x+1)q^{-a}\big]=2q^{-a}-a(2x+1)^2q^{-a-1}=(2-4a)q^{-a}+3aq^{-a-1}$ であり、$[0,1]$ 上で積分すると
$$
3a\,L_{a+1}=3^{1-a}-1+(4a-2)L_a .
$$
GOE：$f_1=q^{-3/2}-q^{-5/2}$。$a=\frac12$ とすると $L_{3/2}=\frac23(\sqrt3-1)$、$a=\frac32$ とすると
$L_{5/2}=\frac29\big(\frac1{\sqrt3}-1+4L_{3/2}\big)=\frac{18\sqrt3-22}{27}$。よって $I_0=L_{3/2}-L_{5/2}=\frac4{27}$ であり、
$\mathbb Er=\frac{27}8\big(\frac43-\frac{16}{9\sqrt3}\big)-\frac12=\frac92-2\sqrt3-\frac12=4-2\sqrt3$ である。
GUE：$f_2=q^{-2}-2q^{-3}+q^{-4}$。$q=(x+\frac12)^2+\frac34$ であるから
$L_1=\frac2{\sqrt3}\big[\arctan\frac{2x+1}{\sqrt3}\big]_0^1=\frac2{\sqrt3}\big(\frac\pi3-\frac\pi6\big)=\frac\pi{3\sqrt3}$ であり、
$a=1,2,3$ での漸化式から $L_2=\frac23L_1=\frac{2\pi}{9\sqrt3}$、$L_3=L_2-\frac19$、$L_4=\frac19\big(10L_3-\frac89\big)$ を得る。よって
$I_0=L_2-2L_3+L_4=\frac{2\pi}{81\sqrt3}$ であり、$\mathbb Er=\frac{81\sqrt3}{4\pi}\cdot\frac8{81}-\frac12=\frac{2\sqrt3}\pi-\frac12$ である。$\square$

(b) の値は Y. Y. Atas, E. Bogomolny, O. Giraud, G. Roux によって求められた（Phys. Rev. Lett. 110
(2013), 084101）。より大きい $n$ に対する平均比は本講義では計算しない。

(a) は区間上の $n$ 個の独立な一様準位について厳密に成り立つ：それらの $n+1$ 個の間隔（両端の
二つの隙間を含む）は、独立な指数変数 $e_j$ による $(e_0,\dots,e_n)/\sum_je_j$ と同じ分布をもつから、連続する間隔の
各比 $s_i/s_{i+1}$ は $e_i/e_{i+1}$ と同じ分布に従い、各間隔比 $r_i$ は (a) の $r$ と同じ分布に従う
（定理 A.25 と系 A.26(b),(d)）。

*計算：* 命題 5.14 のモンテカルロ法と $f_\beta$ の数値積分による検証、および数値的に観察した $n=64$ の GOE の平均比 — [チュートリアル、5.3節](../../tutorial/ja/05_random_matrices.md#s-5-3)。

<a id="s-5-4"></a>

## 5.4 なぜセクターを混ぜてはならないか

$H=\bigoplus_cH_c\in\mathcal M$ を **独立なブロック** をもつランダムな自己共役元、
例えば GOE ブロックとする。$\mathbb C^D$ 上の作用素として、そのスペクトルは各ブロックのスペクトルの
和集合である。ブロックの内部では準位は反発する。異なるブロックの準位は反発しないことを示す。

各ブロックが有界連続な **準位密度** $\rho_c$ をもつと仮定する：
有界ボレル関数 $\varphi$ に対し $\mathbb E\sum_i\varphi(\lambda_i(H_c))=\int\varphi\,\rho_c$（したがって
$\int\rho_c=k_c$）。GOE ブロックと GUE ブロックについては、これは系 C.8 により成り立つ。系 C.8 は定理 C.7 の同時
密度から $\rho_c$ を導く。

<a id="r-5-15"></a>

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

*証明.* 順序づけた固有値は行列の連続関数であるから（系 A.11）、
$\Lambda_c=(\lambda_i(H_c))_i$ と $\Lambda_d$ は独立なランダム・ベクトルである。フビニの定理により
（それらの同時分布はそれぞれの分布の直積である）、$\mathbb E\sum_{i,j}\psi(\Lambda_{c,i},\Lambda_{d,j})=\mathbb E\sum_j g(\Lambda_{d,j})$ である。ここで
$g(y)=\mathbb E\sum_i\psi(\Lambda_{c,i},y)=\int\psi(x,y)\rho_c(x)\,dx$ は有界なボレル関数である。$g$ に
$\rho_d$ の定義を適用すればよい。$\psi(x,y)=1\{|x-y|<\varepsilon\}$ とすると
$$
\frac{\mathbb EN_\varepsilon}{2\varepsilon}=\sum_{c<d}\int\rho_c(x)\Big(\frac1{2\varepsilon}\int_{x-\varepsilon}^{x+\varepsilon}\rho_d(y)\,dy\Big)dx ;
$$
括弧内は連続性により $\rho_d(x)$ に収束し、$\sup\rho_d$ で抑えられ、$\rho_c$ は可積分であるから、
優収束定理により極限が得られる。最後に $N_\varepsilon\le\sum_{c<d}k_ck_d$ であるから
$\mathbb EN_\varepsilon\le\big(\sum_{c<d}k_ck_d\big)P(N_\varepsilon\ge1)$ であり、$N_\varepsilon\ge1$ ならば $\delta<\varepsilon$ である。$\square$

したがって、二つの GOE または GUE ブロックに対しては $P(\delta<\varepsilon)$ は少なくとも $\varepsilon$ のオーダーである：それらの準位密度は
いたるところ正であり（系 C.8 において、$(x,y)$ の二つの座標が一致する零集合の外では $q(x,\cdot)>0$）、
したがって $\int\rho_c\rho_d>0$ である。$2\times2$ GOE ブロックの内部ではそれは $1-e^{-\varepsilon^2/4}\le\varepsilon^2/4$ であり
（命題 5.13、正規化していない間隔）、$n\times n$ GOE ブロックの内部ではそれは高々定数倍の
$\varepsilon^2$ である（命題 C.9）。すべてのブロックが完全な GOE であるにもかかわらず、
合併したスペクトルでは準位反発の検定が失敗する。

*計算：* 独立な GOE ブロック $C$ 個に対し、合併したスペクトルの平均間隔比は $C$ が大きくなるにつれて独立な準位の値 $2\ln2-1$ に向かって下がり、ブロックが独立でないランダム磁場中のハイゼンベルク鎖（第9章）も同じ効果を示す。どちらも数値的観察である — [チュートリアル、5.4節](../../tutorial/ja/05_random_matrices.md#s-5-4)。

> **物理的意味**（物理的解釈であり、本講義の定理ではない）. 保存する量子数が異なる準位は相互作用しないので、
> 自由に交差する。反発するのは同じ対称性の準位だけである。したがって対称性セクターを混ぜた
> スペクトルは、各セクターのダイナミクスが完全にカオス的であってもポアソン的な統計を示し、
> 分離されていない対称性は可積分性や局在と取り違えられうる。
> 「準位統計を計算する前にすべての対称性を分離せよ」という経験則は、
> $M_D$ ではなく環 $\bigoplus_cM_{k_c}$ の中で計算すべきだという主張にほかならない。

<a id="s-5-5"></a>

## 5.5 比較の対象

*ライブラリでは：* 他のソフトウェアとの比較 — [チュートリアル、5.5節](../../tutorial/ja/05_random_matrices.md#s-5-5)。

<a id="x-5"></a>

## 演習問題

1. $U$ を $U(\mathcal M)$ 上でハール分布に従うとし、$\mathrm{Tr}\,U=\sum_c\mathrm{Tr}\,U_c$ とする。
   $\mathbb E|\mathrm{Tr}\,U|^2=C$ かつ $\mathbb E|\mathrm{Tr}\,U|^4=\sum_cm(k_c)+2C(C-1)$ であることを示せ。ただし
   $m(1)=1$、$k\ge2$ に対し $m(k)=2$ である。$U(D)$（$D\ge2$）のハール・ユニタリに対する $1$ と $2$ と比較し、
   $\mathcal M=M_1\oplus M_2\oplus M_3$ について値を求めよ。
2. 定理 5.12 を用いて、$n\le k$ に対し $\mathbb E|U_{11}|^{2n}=\frac{n!}{k(k+1)\cdots(k+n-1)}$ であることを示せ
   （系 B.14 はこれを $n>k$ に拡張する）。
3. **円形直交アンサンブル**は、$W$ を $U(k)$ 上でハール分布に従うとして $S=W^TW$ で与えられる。
   $S$ がユニタリかつ対称であること、およびその分布がすべての
   $V\in U(k)$ に対して $S\mapsto V^TSV$ のもとで不変であることを示せ。$\lambda_1,\dots,\lambda_k$ を $S$ の固有値（重複度を込めて）とし、
   単位円 $\mathbb T$ のボレル集合 $B$ に対し $N_B(S)=\#\{j:\lambda_j\in B\}$ とおく。$N_B(S)$ が
   確率変数であること、$\nu(B)=\frac1k\mathbb EN_B(S)$ が正規化された弧長であること（固有位相の平均密度は
   一様である）、そしてすべての整数 $m\neq0$ に対し $\mathbb E\mathrm{Tr}\,S^m=0$ であることを示せ。

<details>
<summary>解答</summary>

1. $t_c=\mathrm{Tr}\,U_c$ と書く。$t_c$ は独立である（系 5.4）。$U_c$ に
   $e^{i\theta}$ を掛けてもその分布は変わらず、$t_c$ には $e^{i\theta}$ が掛かる。したがってモーメント
   $\mathbb E[t_at_b\bar t_c\bar t_d]$ は、各セクターが $(a,b)$ と
   $(c,d)$ に同じ回数だけ現れる場合、すなわち多重集合として $\{a,b\}=\{c,d\}$ である場合を除いて零であり、同様に $a\neq b$ に対し $\mathbb E[t_a\bar t_b]=0$ である。
   命題 5.7 により $\mathbb E[(U_c)_{ii}\overline{(U_c)_{jj}}]=\langle e_i,\mathbb E[U_ce_ie_j^*U_c^*]e_j\rangle=\delta_{ij}/k_c$ であるから、
   $\mathbb E|t_c|^2=1$、$\mathbb E|\mathrm{Tr}\,U|^2=C$ である。$\mathbb E|\sum_at_a|^4$ において、
   $a=b$ の項は $\sum_a\mathbb E|t_a|^4$ を与え、$a\neq b$ の項は
   $2\sum_{a\neq b}\mathbb E|t_a|^2\,\mathbb E|t_b|^2=2C(C-1)$ を与える（$(c,d)$ の選び方が二通り）。$k=1$ に対しては $|t|=1$ である。
   $k\ge2$ に対しては $\mathbb E|t|^4=\sum_{i,j,i',j'}\mathbb E[U_{ii}U_{jj}\overline{U_{i'i'}U_{j'j'}}]$ である。
   定理 5.10 において、$\sigma=\pi$ の二つの項は同一の仕方で $\{i',j'\}=\{i,j\}$ を強制し、和をとるとそれぞれ
   $k^2\mathrm{Wg}(e)$ を与える。$\sigma\neq\pi$ の二つの項は $i=j=i'=j'$ を強制し、それぞれ
   $k\,\mathrm{Wg}(\tau)$ を与える。合計は $2\big(\frac{k^2}{k^2-1}-\frac1{k^2-1}\big)=2$ である。$C=1$、
   $k=D$ の場合が $U(D)$ に対する $1$ と $2$ を与える。$M_1\oplus M_2\oplus M_3$ に対しては
   $\mathbb E|\mathrm{Tr}\,U|^4=1+2+2+12=17$ であり、$U(6)$ のハール・ユニタリに対する $2$ と対照的である。

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

3. $W^T$ はユニタリであるから $S^*S=W^*\bar WW^TW=1$ であり、$S^T=W^TW=S$ である。$V\in U(k)$ に対し $WV$ は
   $W$ と同じ分布をもち（右不変性、補題 5.1）、$(WV)^T(WV)=V^TSV$ である。$V=e^{i\alpha/2}1$ とすると、すべての実数 $\alpha$ に対し
   $e^{i\alpha}S$ は $S$ と同じ分布をもつ。

   *可測性.* 定理 1.3 により、ユニタリ $Y$ を用いて $S=Y\mathrm{diag}(\lambda_1,\dots,\lambda_k)Y^*$ と書ける。$v$ が
   $Sv=\mu v$ を満たす単位ベクトルならば、$v=\sum_jc_jYe_j$ と書くと $c_j(\lambda_j-\mu)=0$ となるから、$\mu\in K$ のとき $v$ は
   $E_K=\mathrm{span}\{Ye_j:\lambda_j\in K\}$ に属する。また $\dim E_K=N_K(S)$ である。したがって $K\subseteq\mathbb T$
   と $m\ge1$ に対し：$N_K(S)\ge m$ であることと、$S$ が固有値が $K$ に属する $m$ 個の正規直交固有ベクトルをもつことは同値である。$K$ が
   閉集合ならば、この性質をもつユニタリの集合 $F_{K,m}$ は $U(k)$ の中で閉じている：$S_l\to S$ であり、$S_l$ が
   固有値 $\mu_i^{(l)}\in K$（$i\le m$）をもつ正規直交固有ベクトル $v_i^{(l)}$ をもつならば、ある部分列はすべての $i\le m$ に対し
   $v_i^{(l)}\to v_i$、$\mu_i^{(l)}\to\mu_i\in K$ を満たし（コンパクト集合
   $(\text{単位球面}\times K)^m$ におけるボルツァーノ・ワイエルシュトラスの定理）、極限は正規直交で $Sv_i=\mu_iv_i$ を満たす。
   $0\le\theta_1(S)\le\dots\le\theta_k(S)<2\pi$ を固有位相とし、固有値が
   $e^{i\theta_j(S)}$ となるようにする（系 1.4：$|\lambda_j|=1$）。$0\le t<2\pi$ に対し
   弧 $K_t=\{e^{i\varphi}:0\le\varphi\le t\}$ は閉集合であり、$\{\theta_j\le t\}=\{N_{K_t}\ge j\}=F_{K_t,j}$ であるから、各
   $\theta_j$ は $U(k)$ 上のボレル関数である。したがって $N_B(S)=\sum_j1_B(e^{i\theta_j(S)})$ は確率変数であり、
   $\nu(B)=\frac1k\sum_jP(e^{i\theta_j(S)}\in B)$ は $\mathbb T$ 上のボレル確率測度、すなわち $k$ 個の固有値の分布の
   平均である。

   *一様性.* $e^{i\alpha}S$ の固有値は $e^{i\alpha}\lambda_j$ であるから $N_B(e^{i\alpha}S)=N_{e^{-i\alpha}B}(S)$ であり、
   $e^{i\alpha}S$ は $S$ と同じ分布をもつから $\nu(e^{-i\alpha}B)=\frac1k\mathbb EN_B(e^{i\alpha}S)=\nu(B)$ である。よって $\nu$ は
   回転で不変であり、命題 A.22（$N=1$）によりそれは正規化された弧長である：
   長さ $|B|$ のすべての弧に対し $\mathbb EN_B(S)=k\,|B|/2\pi$ である。同様に $\mathrm{Tr}\,S^m=\sum_j\lambda_j^m$（$S^{-1}=S^*$ であるから
   $m<0$ に対しても）は $S$ の有界連続関数であり、すべての $\alpha$ に対し
   $\mathbb E\mathrm{Tr}\,S^m=\mathbb E\mathrm{Tr}(e^{i\alpha}S)^m=e^{im\alpha}\,\mathbb E\mathrm{Tr}\,S^m$ であるから、
   $m\neq0$ に対しそれは零である：これらは $\nu$ のフーリエ係数の $k$ 倍である。特に平均
   固有位相密度は一定であり、円周に沿った間隔はアンフォールディングを必要としない。

</details>

*計算：* 演習問題 1–3 の数値的な部分 — [チュートリアル、演習問題](../../tutorial/ja/05_random_matrices.md#x-5)。

**次へ：**[第6章](06_quantum_channels.md) では量子チャネルを扱う：クラウス表現とスタインスプリング
表現、トレースに関する双対、そしてセクター間の写像。ここでハール平均として現れた条件付き
期待値 $E$ が最初の例である。
