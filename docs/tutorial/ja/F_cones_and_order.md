[English](../F_cones_and_order.md) · [Español](../es/F_cones_and_order.md) · [Français](../fr/F_cones_and_order.md) · [Deutsch](../de/F_cones_and_order.md) · [中文](../zh/F_cones_and_order.md) · **日本語** · [Русский](../ru/F_cones_and_order.md)

# 付録 F. 対称錐、順序同型とジョルダン写像

第13章では、オブザーバブルの順序とそのジョルダン積を三つの定理によって結びつけた：ケッヒャー–ヴィンベルクの
定理（定理 13.8）、カディソンの定理（定理 13.10）、そして $M_n(\mathbb C)$ のジョルダン自己同型の記述
（命題 13.11 は $n=2$ の場合にそれを証明する）である。本付録では、一つ目を完全な一般性で証明し、残りの二つを
本講義の環 $\bigoplus_cM_{k_c}(\mathbb C)$ に対して、すべてのブロックサイズについて証明する。その途中で、
13.6節が必要とするユークリッド的ジョルダン代数についての事実、すなわちスペクトル定理、平方の錐の自己双対性、
および $P(y)Q\subseteq Q$ を証明する。

**本付録で学ぶこと**：すべてのユークリッド的ジョルダン代数が冪結合性とスペクトル定理をもつ理由、その平方の錐が
自己双対かつ等質である理由、自己双対で等質な錐だけから、その特性関数の対数の3階導関数を通じてジョルダン積が
どのように復元されるか、正値性を両方向に保つ単位的な線形全単射が平方を保つ理由、そして $\bigoplus_cM_{k_c}$ の
そのような写像がすべて等しいブロックを置換し、各ブロック上で $UxU^*$ または $Ux^TU^*$ として作用する理由。

線形代数、微積分、そして番号で引用する第1章、第6章、第13章の結果のほかに、F.3節ではルベーグ積分
（ファトゥの補題、優収束定理と積分記号下での微分、フビニ–トネリの定理、線形変数変換の公式）を用いる。
それ以外は何も仮定しない。付録 E のジョルダン枠によるスペクトル定理は下の定理 F.3 を精密化するものであるが、
ここでの議論はそれに依存しない。
**ユークリッド空間** とは、内積 $\langle\cdot,\cdot\rangle$ をもつ有限次元実ベクトル空間 $V\ne\{0\}$ のことである。
$\|x\|$ はノルムであり（$V$ の元に対する一本線 $|x|$ は定理 F.3 の意味での絶対値を表し、決してノルムを
表さない）、$\|A\|$ は $V$ の線形写像 $A$ の作用素ノルム、$A^T$ はその随伴
（$\langle Ax,y\rangle=\langle x,A^Ty\rangle$）であり、$A\ge0$（$A>0$）は $A=A^T$ が半正定値
（正定値）であることを意味する。$e^A=\exp A=\sum_{k\ge0}A^k/k!$ とする。

## F.1 ユークリッド的ジョルダン代数

13.1節と同様に、**ジョルダン代数** とは双線形で可換な積 $\circ$ をもつ実ベクトル空間であって、
$x\circ(x^2\circ y)=x^2\circ(x\circ y)$、すなわち $[L(x),L(x^2)]=0$ を満たすものである。ここで $L(x)y=x\circ y$、$x^2=x\circ x$ であり、
$P(x)=2L(x)^2-L(x^2)$ とする。

**定義.** **ユークリッド的ジョルダン代数** とは、単位元 $e$（$e\circ x=x$）と、すべての $L(x)$ が対称となる内積
$\langle x\circ y,z\rangle=\langle y,x\circ z\rangle$ をもつ有限次元ジョルダン代数 $V$ のことである。
その **平方の錐** は $Q=\{x^2:x\in V\}$ である。

例：$M_{sa}$（$\bigoplus_cM_{k_c}(\mathbb C)$ の自己共役部分）と $\operatorname{Herm}_n(\mathbb R)$ に
$\langle x,y\rangle=\operatorname{Tr}(xy)$ を入れたもの。これは実数値で $\operatorname{Tr}(x^2)=\sum_{ij}|x_{ij}|^2$ を満たす（補題 13.1。
$\operatorname{Tr}((x\circ y)z)=\frac12\operatorname{Tr}(xyz+yxz)=\operatorname{Re}\operatorname{Tr}(xyz)$ は、
$\operatorname{Tr}(yxz)=\overline{\operatorname{Tr}((yxz)^*)}=\overline{\operatorname{Tr}(zxy)}=\overline{\operatorname{Tr}(xyz)}$
であるから、巡回置換で不変であり、同じ理由で $x\leftrightarrow z$ でも不変であるので、
$x,y,z$ について対称である）。スピン因子 $V_m$（補題 13.3 によりジョルダン代数）に内積
$\langle x,y\rangle=\operatorname{tr}(x\circ y)=2(ab+u\cdot v)$（$x=(a,u)$、$y=(b,v)$）を入れたもの。実際、$z=(c,w)$ に対し
$\operatorname{tr}((x\circ y)\circ z)=2(abc+c\,u\cdot v+a\,v\cdot w+b\,u\cdot w)$ は
$x,y,z$ について対称である。そして座標ごとの積をもつ $\mathbb R^N$（演習問題 1）。

**補題 F.1.** すべてのジョルダン代数において、すべての $a,b,d,x,y$ に対し：

1. $[L(a),L(b\circ d)]+[L(b),L(d\circ a)]+[L(d),L(a\circ b)]=0$；
2. $L(x^2\circ y)=L(x^2)L(y)-2L(x)L(y)L(x)+2L(x\circ y)L(x)$。

*証明.* (1) $s,t,r\in\mathbb R$ と $x=sa+tb+rd$ に対し、$[L(x),L(x^2)]=0$ は $(s,t,r)$ についての多項式恒等式であるから、
各係数は零である。$str$ の係数は (1) の左辺の2倍である。(2) (1) で $a=b=x$ とおくと
$2[L(x),L(x\circ d)]+[L(d),L(x^2)]=0$ となる。これを $y$ に適用すると
$2x\circ((x\circ d)\circ y)-2(x\circ d)\circ(x\circ y)+d\circ(x^2\circ y)-x^2\circ(d\circ y)=0$、すなわちすべての $d$ に対し
$\big(2L(x)L(y)L(x)-2L(x\circ y)L(x)+L(x^2\circ y)-L(x^2)L(y)\big)d=0$ となる。$\square$

**命題 F.2（冪結合性）.** $V$ を単位元 $e$ をもつジョルダン代数、$x\in V$ とし、$x^0=e$、
$x^{k+1}=x\circ x^k$ とする。このとき、すべての $p,q\ge0$ に対し $x^p\circ x^q=x^{p+q}$ であり、各 $L(x^p)$ は
可換な作用素 $L(x)$ と $L(x^2)$ の多項式である。したがって実多項式 $p,q$ に対し $p(x)\circ q(x)=(pq)(x)$ であり、作用素
$L(p(x))$ どうしは可換である。

*証明.* ジョルダン恒等式により $L(x)$ と $L(x^2)$ は可換である。$N\ge2$ についての帰納法で次を示す：(A$_N$)
$i+j\le N$ に対し $x^i\circ x^j=x^{i+j}$、および (B$_N$) $i\le N$ に対し $L(x^i)$ は $L(x),L(x^2)$ の多項式である。どちらも
$N=2$ では明らかである。$N$ について両者を仮定し、$i,j\ge1$ として $i+j=N+1$ とする（$i=0$ または $j=0$ のときは示すことはない）。
$i\le N$ であるから $L(x^i)$ は $L(x)$ と可換であり、(A$_N$) により
$$
x^i\circ x^j=L(x^i)L(x)x^{j-1}=L(x)L(x^i)x^{j-1}=x\circ x^N=x^{N+1}.
$$
これが (A$_{N+1}$) である。$y=x^{N-1}$ とした補題 F.1(2) と $x^2\circ x^{N-1}=x^{N+1}$ から
$L(x^{N+1})=L(x^2)L(x^{N-1})-2L(x)L(x^{N-1})L(x)+2L(x^N)L(x)$ が得られ、これが (B$_{N+1}$) である。最後の文は
双線形性から従う。$\square$

**定理 F.3（スペクトル定理）.** $V$ をユークリッド的ジョルダン代数とし、$x\in V$ とする。
$k\ge1$、実数 $\lambda_1>\dots>\lambda_k$ および零でない $c_1,\dots,c_k\in V$ で
$$
c_i\circ c_j=\delta_{ij}c_i,\qquad \textstyle\sum_ic_i=e,\qquad x=\sum_i\lambda_ic_i .
$$
を満たすものが一意的に存在する。$c_i$ は $x$ の多項式であり、すべての実多項式 $p$ に対し $p(x)=\sum_ip(\lambda_i)c_i$ である。

$(c_i)$ を **直交冪等元の完全系**、$\lambda_i$ を $x$ の **固有値** と呼び、
$f:\{\lambda_1,\dots,\lambda_k\}\to\mathbb R$ に対し $f(x)=\sum_if(\lambda_i)c_i$ とおく。例えば $|x|$、$x_\pm=\max(\pm x,0)$、
$\sqrt x$（すべての $\lambda_i\ge0$ のとき）、$\exp x$、そして $\log x$（すべての $\lambda_i>0$ のとき）である。

*証明.* 存在。$\mathbb R[x]=\{p(x)\}$ は $e$ を含む部分空間であり、命題 F.2 により
$p(x)\circ q(x)=(pq)(x)$ を満たす。よって $M=L(x)|_{\mathbb R[x]}$ は $\mathbb R[x]$ をそれ自身に写し、対称であり、
$p(M)q(x)=p(x)\circ q(x)$ を満たす。対称作用素は実固有値で対角化可能である。$M$ の固有値を
$\lambda_1>\dots>\lambda_k$ とし、$m(t)=\prod_i(t-\lambda_i)$ をその最小多項式とする。すると $m(x)=m(M)e=0$ である。ラグランジュ
多項式 $\ell_i(t)=\prod_{j\ne i}\frac{t-\lambda_j}{\lambda_i-\lambda_j}$ を用いると、多項式
$\ell_i\ell_j-\delta_{ij}\ell_i$、$\sum_i\ell_i-1$、$\sum_i\lambda_i\ell_i-t$ はすべての $\lambda_l$ で零となるので、
$m$ の倍元である。$x$ で値をとれば、$c_i=\ell_i(x)$ が三つの性質をもつことがわかる。$c_i=0$ ならば、すべての $q$ に対し
$\ell_i(M)q(x)=\ell_i(x)\circ q(x)=0$ となるので、$\deg\ell_i<\deg m$ であるにもかかわらず $\ell_i(M)=0$ となる。よって $c_i\ne0$ である。
最後に $p(x)=\sum_ip(\lambda_i)c_i$ である：これは $p=1$ と $p=t$ に対して正しく、$x\circ\sum_i\lambda_i^rc_i=\sum_i\lambda_i^{r+1}c_i$ である。

一意性。相異なる $\mu_j$ と零でない直交冪等元の完全系 $d_j$ によって $x=\sum_j\mu_jd_j$ と書けたとする。
上と同様に $p(x)=\sum_jp(\mu_j)d_j$ である。$d_j$ は一次独立である（$\sum_ja_jd_j=0$ に $d_i$ を掛ける）から、
$p(x)=0$ であることと、すべての $j$ に対し $p(\mu_j)=0$ であることは同値である。よって $m(x)=0$ から $\{\mu_j\}\subseteq\{\lambda_i\}$ が従い、
$\tilde m(t)=\prod_j(t-\mu_j)$ は $\tilde m(M)q(x)=\tilde m(x)\circ q(x)=0$ を満たすので、$m$ は $\tilde m$ を割り切り、$\{\lambda_i\}\subseteq\{\mu_j\}$ である。
したがって $d_j=\ell_j(x)=c_j$ である。$\square$

**補題 F.4（ピアス分解）.** $V$ をユークリッド的ジョルダン代数とし、$c$ を冪等元（$c^2=c$）とする。

1. $L(c)(2L(c)-I)(L(c)-I)=0$ である。したがって $V$ は固有空間 $V(c,\nu)=\{z:c\circ z=\nu z\}$
   （$\nu\in\{0,\frac12,1\}$）の直交和であり、$L(c)\ge0$ である。
2. $V(c,1)\circ V(c,0)=\{0\}$。$\nu\in\{0,1\}$ に対し $V(c,\nu)\circ V(c,\nu)\subseteq V(c,\nu)$ かつ $V(c,\nu)\circ V(c,\frac12)\subseteq V(c,\frac12)$。
   $V(c,\frac12)\circ V(c,\frac12)\subseteq V(c,0)+V(c,1)$。
3. $s=2c-e$ は $s^2=e$ を満たし、$P(s)$ は $V$ のジョルダン自己同型である：それは $V(c,0)+V(c,1)$ 上で $+1$、
   $V(c,\frac12)$ 上で $-1$ であり、$P(s)(z\circ w)=P(s)z\circ P(s)w$ を満たす。
4. $x=\sum_{i=1}^k\lambda_ic_i$ を定理 F.3 のとおりとする。$L(c_i)$ たちは可換であり、$V$ は部分空間
   $V_{ii}=\{z:c_i\circ z=z\}$ と $V_{ij}=\{z:c_i\circ z=c_j\circ z=\frac12z\}$（$i<j$）の直交和である。
   $a=\sum_la_lc_l$（$a_l\in\mathbb R$）に対し、$L(a)$ は $V_{ij}$（$i\le j$）上でスカラー $\frac12(a_i+a_j)$ として、$P(a)$ は
   スカラー $a_ia_j$ として作用する。

*証明.* (1) $x=y=c$ とした補題 F.1(2) により $L(c)=L(c)^2-2L(c)^3+2L(c)^2$、すなわち $2L(c)^3-3L(c)^2+L(c)=0$ であり、これは
述べた因数分解である。よって対称作用素 $L(c)$ の固有値は $\{0,\frac12,1\}$ に属する。

(2) $a=b=c$、$d=z\in V(c,\nu)$ とした補題 F.1(1) により $2\nu[L(c),L(z)]+[L(z),L(c)]=0$ となるので、
$\nu\in\{0,1\}$ ならば $L(z)$ は $L(c)$ と可換であり、したがってすべての $V(c,\mu)$ を保つ。これから二番目と三番目の包含関係、および
$V(c,1)\circ V(c,0)\subseteq V(c,0)\cap V(c,1)=\{0\}$ が得られる。$z,w\in V(c,\frac12)$ に対しては、$(a,b,d)=(z,w,c)$ とした補題 F.1(1)
を $c$ に適用する：$L(w\circ c)=\frac12L(w)$、$L(c\circ z)=\frac12L(z)$ であるから最初の二項は打ち消し合い、
$u=z\circ w$ に対し $[L(c),L(u)]c=L(c)^2u-L(c)u=0$ となる。(1) により $u=u_1+u_{1/2}+u_0$ と書けば、これは $-\frac14u_{1/2}=0$ である。

(3) $s^2=4c-4c+e=e$ であり、$L(s)=2L(c)-I$ は $V(c,1),V(c,\frac12),V(c,0)$ 上でそれぞれ $1,0,-1$、そして $P(s)=2L(s)^2-I$ である。
$z',w'\in V(c,0)+V(c,1)$、$z'',w''\in V(c,\frac12)$ として $z=z'+z''$、$w=w'+w''$ と書く。(2) により $z'\circ w'+z''\circ w''$ は
$V(c,0)+V(c,1)$ に、$z'\circ w''+z''\circ w'$ は $V(c,\frac12)$ に属するので、
$P(s)(z\circ w)=(z'-z'')\circ(w'-w'')=P(s)z\circ P(s)w$ である。

(4) $c_l$ は $x$ の多項式であるから、$L(c_l)$ たちは可換である（命題 F.2）。それらは対称であるから、共通の
正規直交固有基底をもつ。共通の固有ベクトル上で、$L(c_l)$ は (1) により固有値 $\varepsilon_l\in\{0,\frac12,1\}$ をもち、
$\sum_lL(c_l)=L(e)=I$ であるから $\sum_l\varepsilon_l=1$ である。よって、ある一つの $\varepsilon_i=1$ で他は零であるか、
ある一組 $i<j$ に対し $\varepsilon_i=\varepsilon_j=\frac12$ で他は零であるかのいずれかである。よって $V$ は
共通固有空間 $W_{ii}$（そこでは $L(c_i)=1$、$l\ne i$ に対し $L(c_l)=0$）と $W_{ij}$（そこでは $L(c_i)=L(c_j)=\frac12$、
それ以外では $L(c_l)=0$）の直交和であり、$W_{ij}=V_{ij}$ である：$c_i\circ z=z$（それぞれ $c_i\circ z=c_j\circ z=\frac12z$）を満たすベクトルは、
$L(c_i)$（それぞれ組 $(L(c_i),L(c_j))$）が異なる値をとる他の共通固有空間に成分をもたないからである。
$a^2=\sum_la_l^2c_l$ であるから、$V_{ij}$（$i<j$）上で $L(a)=\frac12(a_i+a_j)$、$L(a^2)=\frac12(a_i^2+a_j^2)$、
$P(a)=\frac12(a_i+a_j)^2-\frac12(a_i^2+a_j^2)=a_ia_j$ であり、$V_{ii}$ 上で $L(a)=a_i$、$P(a)=2a_i^2-a_i^2$ である。$\square$

## F.2 平方の錐は対称錐である

$\Omega$ をユークリッド空間 $V$ 内の空でない開 **凸錐**（$x\in\Omega$、$t>0$ に対し $tx\in\Omega$）とする。その
（開）**双対錐** は $\Omega^*=\{y:\langle x,y\rangle>0\text{（すべての }x\in\bar\Omega\setminus\{0\}\text{ に対し）}\}$ である。$\Omega^*=\Omega$ のとき $\Omega$ は
**自己双対** であるという。$G(\Omega)=\{g\in GL(V):g\Omega=\Omega\}$ とし、$G(\Omega)$ が $\Omega$ 上に推移的に作用するとき $\Omega$ は
**等質** であるという。**対称錐** とは、自己双対かつ等質な錐のことである。$\Omega_1\subseteq\Omega_2$ ならば
$\Omega_2^*\subseteq\Omega_1^*$ である。

**定理 F.5（自己双対性）.** $V$ をユークリッド的ジョルダン代数とし、$y=\sum_i\lambda_ic_i$（定理 F.3）とする。

1. 次は同値である：(i) $y\in Q$；(ii) すべての $\lambda_i\ge0$；(iii) $L(y)\ge0$；(iv) すべての $x\in Q$ に対し
   $\langle y,x\rangle\ge0$。特に $Q$ は閉凸錐である。
2. 次は同値である：(i) $y\in\Omega:=\operatorname{int}Q$；(ii) すべての $\lambda_i>0$；(iii) $L(y)>0$；(iv)
   すべての $x\in Q\setminus\{0\}$ に対し $\langle y,x\rangle>0$。さらに $e\in\Omega$、$\bar\Omega=Q$、$\Omega^*=\Omega$ である：
   $\Omega$ は自己双対な開凸錐である。
3. $x,y\in Q$ かつ $\langle x,y\rangle=0$ ならば $x\circ y=0$ である。
4. $x^2=0$ ならば $x=0$ であり、$Q\cap(-Q)=\{0\}$ である。

*証明.* $\langle y,x^2\rangle=\langle y\circ x,x\rangle=\langle L(y)x,x\rangle$ と
$\langle x^2,e\rangle=\langle x,x\circ e\rangle=\|x\|^2$ を用いる。後者は (1)(ii) とあわせて (4) を与える。一意性により
$-y$ の固有値は $-\lambda_i$ たちだからである。

(1) (i)$\Rightarrow$(iii)：$x=\sum_i\mu_id_i$ により $y=x^2$ ならば、$y=\sum_i\mu_i^2d_i$ であり、補題 F.4(1) により
$L(y)=\sum_i\mu_i^2L(d_i)\ge0$ である。(iii)$\Leftrightarrow$(iv) は最初の恒等式による。
(iii)$\Rightarrow$(ii)：$0\le\langle L(y)c_i,c_i\rangle=\lambda_i\|c_i\|^2$。(ii)$\Rightarrow$(i)：
$y=(\sum_i\sqrt{\lambda_i}c_i)^2$。(iv) により、$Q$ は $0$ を通る閉半空間の共通部分である。

(2) (ii)$\Rightarrow$(iii)：$L(y)=\sum_i\lambda_iL(c_i)\ge(\min_i\lambda_i)\sum_iL(c_i)=(\min_i\lambda_i)I$。
(iii)$\Rightarrow$(iv)：零でない $x\in Q$ は $z\ne0$ により $z^2$ と書け、$\langle y,x\rangle=\langle L(y)z,z\rangle>0$ である。
(iv)$\Rightarrow$(ii)：$\langle y,c_i\rangle=\lambda_i\|c_i\|^2$ かつ $c_i=c_i^2\in Q\setminus\{0\}$。(iii)$\Rightarrow$(i)：$m>0$ により
$L(y)\ge mI$ ならば $L(y')\ge(m-\|L\|\,\|y-y'\|)I$ である。ここで $\|L\|$ は線形写像
$b\mapsto L(b)$ のノルムである。よって $y$ のある近傍は
$L>0$ を満たす元からなり、それらは (1) により $Q$ に属する。(i)$\Rightarrow$(iii)：(1) により $L(y)\ge0$ である。$z\ne0$ に対し $\langle L(y)z,z\rangle=0$ ならば、
$x=z^2\ne0$ かつ $\langle y,x\rangle=0$ である。小さな $\varepsilon>0$ に対し $y-\varepsilon x\in Q$ であるが、
$\langle y-\varepsilon x,x\rangle=-\varepsilon\|x\|^2<0$ となり、(1)(iv) に矛盾する。次に $L(e)=I>0$ であり、$y\in Q$ に対し
$y+\varepsilon e\in\Omega$ は $y$ に収束するので $\bar\Omega=Q$ である。すると (iv) により $\Omega^*=\{y:\langle y,x\rangle>0\ \forall x\in Q\setminus\{0\}\}=\Omega$ である。
$\Omega=\{L>0\}$ は凸である。

(3) $\lambda_i\ge0$ として $x=\sum_i\lambda_ic_i$ とする。すると $0=\sum_i\lambda_i\langle c_i,y\rangle$ であり、(1) によりすべての項は $\ge0$ である。
$\lambda_i>0$ ならば $0=\langle y,c_i\rangle=\langle L(y)c_i,c_i\rangle$ かつ $L(y)\ge0$ であるから $L(y)c_i=0$ である。よって
$x\circ y=\sum_i\lambda_i\,c_i\circ y=0$ である。$\square$

**定理 F.6（等質性と二次表現）.** $V$ をユークリッド的ジョルダン代数とし、$Q$ と
$\Omega$ を定理 F.5 のとおりとする。

1. すべての $b\in V$ に対し：$e^{L(b)}Q=Q$、$e^{L(b)}\Omega=\Omega$、$e^{L(b)}e=\exp b$。
2. すべての $y\in\Omega$ は $\exp(\log y)$ に等しい。したがって $\Omega$ は等質であり、対称錐である。
3. $y\in\Omega$ に対し $P(y)=e^{2L(\log y)}$ である。よって $P(y)Q=Q$、$P(y)\Omega=\Omega$、$P(y)e=y^2$ である。
4. すべての $y\in V$ に対し $P(y)Q\subseteq Q$ である。

*証明.* (1) $0\le s<1/\|L(b)\|$ とする（$L(b)=0$ ならば任意の $s\ge0$）。対称作用素 $I-sL(b)$ の固有値は
$\ge1-s\|L(b)\|>0$ である。$x\in Q$ とし、$z=(I-sL(b))^{-1}x=\sum_j\nu_jd_j$（定理 F.3）とする。$z\circ d_j=\nu_jd_j$、
$\langle d_l,d_j\rangle=\langle d_l,d_j\circ d_j\rangle=\langle d_l\circ d_j,d_j\rangle=\delta_{lj}\|d_j\|^2$ および
$\langle b,d_j\rangle=\langle b,d_j^2\rangle=\langle L(b)d_j,d_j\rangle\le\|L(b)\|\,\|d_j\|^2$ を用いると
$$
0\le\langle x,d_j\rangle=\langle z-s\,b\circ z,d_j\rangle=\nu_j\|d_j\|^2-s\langle b,z\circ d_j\rangle
=\nu_j\big(\|d_j\|^2-s\langle b,d_j\rangle\big),
$$
となる。ここで左側の不等式は定理 F.5(1)(iv) である。括弧は $\ge(1-s\|L(b)\|)\|d_j\|^2>0$ であるから、すべての $\nu_j\ge0$ であり、
$z\in Q$ である。したがって $t\ge0$ と $n>t\|L(b)\|$ に対し $(I-\frac tnL(b))^{-n}Q\subseteq Q$ である。$n\to\infty$ のときこの作用素は
$e^{tL(b)}$ に収束し（$L(b)$ の正規直交固有基底において：$(1-t\beta/n)^{-n}\to e^{t\beta}$）、$Q$ は閉であるから
$e^{tL(b)}Q\subseteq Q$ である。$b$ の代わりに $-b$ を用いると
$e^{-L(b)}Q\subseteq Q$ となるので $e^{L(b)}Q=Q$ であり、同相写像 $e^{L(b)}$ は $\operatorname{int}Q$ をそれ自身の上に写す。
最後に $L(b)^re=b^r$（ジョルダン冪）であるから、定理 F.3 により $b=\sum_i\beta_id_i$ に対し
$e^{L(b)}e=\sum_rb^r/r!=\sum_ie^{\beta_i}d_i=\exp b$ である。

(2) $\lambda_i>0$ として $y=\sum_i\lambda_ic_i$ ならば、$\log y=\sum_i(\log\lambda_i)c_i$ であり $\exp\log y=y$ である。よって
$y=e^{L(\log y)}e$ であり、(1) により $e^{L(\log y)}\in G(\Omega)$ である。自己双対性は定理 F.5(2) である。

(3) $y$ の $c_i$ に対する補題 F.4(4) により、$P(y)$ と $e^{2L(\log y)}$ はどちらも $V_{ij}$ 上で
$\lambda_i\lambda_j=e^{\log\lambda_i+\log\lambda_j}$ として作用する。(1) と $P(y)e=2y\circ y-y^2=y^2$ を用いればよい。

(4) $y=\sum_i\lambda_ic_i$、$c=\sum_{\lambda_i\ge0}c_i$（冪等元）、$\varepsilon_i=\pm1$ により $s=2c-e=\sum_i\varepsilon_ic_i$、
そして $|y|=\sum_i|\lambda_i|c_i$ とする。補題 F.4(4) により、$P(s)$、$P(|y|)$、$P(y)$ は $V_{ij}$ 上でそれぞれ
$\varepsilon_i\varepsilon_j$、$|\lambda_i\lambda_j|$、$\lambda_i\lambda_j$ として作用するので、$P(y)=P(s)P(|y|)$ である。$\eta>0$ に対し
$|y|+\eta e\in\Omega$ であるから、(3) により $P(|y|+\eta e)Q\subseteq Q$ である。$P(\cdot)$ は多項式であり $Q$ は閉であるから
$P(|y|)Q\subseteq Q$ である（$\eta\to0$）。$P(s)$ はジョルダン自己同型である（補題 F.4(3)）から、$P(s)x^2=(P(s)x)^2\in Q$ である。$\square$

13.6節では、$V_m$ において三つの性質から $\operatorname{tr}|P(y)x|\le\operatorname{tr}P(y)|x|$ を証明した。それらの性質は、
線形形式 $\tau(x)=\langle x,e\rangle$ に対してすべてのユークリッド的ジョルダン代数で成り立つ。この線形形式は
$M_{sa}$ 上では $\operatorname{Tr}$、$V_m$ 上では $\operatorname{tr}$ である（上の内積に関して）：

**系 F.7（すべてのユークリッド的ジョルダン代数における命題 13.14）.** ユークリッド的ジョルダン代数の $x,y$ に対し
$\tau(|P(y)x|)\le\tau(P(y)|x|)$ であり、$x\in Q$ または $-x\in Q$ ならば等号が成り立つ。

*証明.* $x_\pm\in Q$、$x=x_+-x_-$、$|x|=x_++x_-$ である（定理 F.3）。$A=P(y)x_+$ と $B=P(y)x_-$ は $Q$ に属し（定理 F.6(4)）、
$z=P(y)x=A-B$ である。$z=\sum_i\mu_ic_i$、$c=\sum_{\mu_i>0}c_i$ とすると、$|z|=z\circ(2c-e)$ であり、
$\tau(a\circ b)=\langle a,b\rangle$ である。$\langle Q,Q\rangle\ge0$（定理 F.5(1)）かつ $c,e-c\in Q$ であるから
$$
\tau(|z|)=\langle A,c\rangle-\langle A,e-c\rangle-\langle B,c\rangle+\langle B,e-c\rangle\le\tau(A)+\tau(B)=\tau(P(y)|x|).
$$
$x\in Q$ ならば $x_-=0$ かつ $z=A\in Q$ であるから $|z|=z=P(y)|x|$ である。$-x\in Q$ の場合も同様である。$\square$

コードは $\operatorname{Herm}_4(\mathbb C)$ とローレンツ錐 $Q\subseteq V_4$ に対して定理 F.5 を（負のペアリングは
(1), $\langle y,c_i\rangle=\lambda_i\|c_i\|^2<0$ の証拠 $c_i$ である）、$V_4$ において定理 F.6 を（座標
$(a,u)$ において作用素 $L(b)$ と $P(y)$ は $5\times5$ 対称行列である）、そして
$\operatorname{Herm}_5(\mathbb R)$ において系 F.7 を検証する：

```python
import math
import torch
from torch_vn_algebra import SpinFactor

torch.manual_seed(0)
dt = torch.complex128
f64 = torch.float64
close = lambda A, B, tol=1e-10: (A - B).abs().max().item() < tol

G = torch.randn(2, 500, 4, 4, dtype=dt)
X, Y = G[0] @ G[0].mH, G[1] @ G[1].mH                             # squares in Herm_4(C)
assert ((X @ Y).diagonal(dim1=-2, dim2=-1).sum(-1).real > -1e-10).all()   # <Q, Q> >= 0
Z = torch.randn(500, 4, 4, dtype=dt)
Z = Z + Z.mH
w, U = torch.linalg.eigh(Z)
c1 = U[:, :, :1] @ U[:, :, :1].mH                                  # eigenprojection of the smallest eigenvalue
pair = (Z @ c1).diagonal(dim1=-2, dim2=-1).sum(-1).real            # <Z, c1> = lambda_min ||c1||^2
assert close(pair, w[:, 0]) and (pair < 0).float().mean() > 0.9   # witnesses of Z not in Q

V = SpinFactor(4)                                                  # <x, y> = tr(x o y) = V.inner
x, y = V.random(2000, cone=True), V.random(2000, cone=True)
assert (V.inner(x, y) > -1e-12).all()
z = V.random(2000)
lam, idem = V.spectral_decomposition(z)                            # (lambda-, lambda+), (c-, c+); ||c-||^2 = 1
assert close(V.inner(z, idem[:, 0]), lam[:, 0])
assert ((lam[:, 0] < 0) == ~V.in_cone(z, tol=0.0)).all()

I5 = torch.eye(5, dtype=f64)
Lmat = lambda b: V.product(b.expand(5, 5), I5).T                   # columns b o e_k
Pmat = lambda b: V.quadratic_representation(b.expand(5, 5), I5).T
for _ in range(20):
    b, yi, yg = V.random(), V.random(cone=True), V.random()
    E = torch.linalg.matrix_exp(Lmat(b))
    assert close(E @ V.unit(), V.exp(b))                           # Theorem F.6(1)
    assert V.in_cone((E @ x.T).T, tol=1e-9).all()
    assert close(Pmat(yi), torch.linalg.matrix_exp(2 * Lmat(V.log(yi))))      # Theorem F.6(3)
    lam_g, idem_g = V.spectral_decomposition(yg)
    s = (torch.where(lam_g >= 0, 1.0, -1.0)[:, None] * idem_g).sum(0)
    assert close(Pmat(yg), Pmat(s) @ Pmat(V.apply_function(yg, torch.abs)))  # P(y) = P(s) P(|y|)
    assert V.in_cone((Pmat(yg) @ x.T).T, tol=1e-9).all()                     # Theorem F.6(4)

xs = torch.randn(300, 5, 5, dtype=f64)
xs, ys = xs + xs.mT, torch.randn(300, 5, 5, dtype=f64)
ys = ys + ys.mT                                                    # Corollary F.7 in Herm_5(R), tau = Tr
absm = lambda A: (lambda w, U: U @ torch.diag_embed(w.abs()) @ U.mT)(*torch.linalg.eigh(A))
lhs = absm(ys @ xs @ ys).diagonal(dim1=-2, dim2=-1).sum(-1)
rhs = (ys @ absm(xs) @ ys).diagonal(dim1=-2, dim2=-1).sum(-1)
assert (rhs - lhs > -1e-9 * rhs).all()
```

## F.3 ケッヒャー–ヴィンベルクの定理

F.2節では、ユークリッド的ジョルダン代数の平方の錐の内部が対称錐であることを証明した。逆は
錐からジョルダン積を再構成する。定理 F.12 までは、$\Omega$ は次元 $N$ のユークリッド空間 $V$ 内の自己双対な開凸錐であり、
$dy$ は正規直交座標でのルベーグ測度である（座標の直交変換は $|\det|=1$ を満たす）。

**補題 F.8.** $\Omega$ を自己双対とする。

1. $0\notin\Omega$ である。$x,y\in\bar\Omega$ に対し $\langle x,y\rangle\ge0$ である。各 $x\in\Omega$ に対し、すべての
   $y\in\bar\Omega$ について $\langle x,y\rangle\ge\delta\|y\|$ となる $\delta>0$ が存在する。
2. $x\in\bar\Omega\setminus\Omega$ ならば、$\langle x,u\rangle=0$ となる $u\in\bar\Omega\setminus\{0\}$ が存在する。$y\in\Omega$、
   $u\in\bar\Omega$、$t\ge0$ ならば $y+tu\in\Omega$ である。
3. $\Omega$ は $\bar\Omega$ の内部である。
4. $G(\Omega)$ は $GL(V)$ の閉部分群であり、$g\in G(\Omega)$ ならば $g^T\in G(\Omega)$ である。

*証明.* (1) $0\in\Omega$ ならば $\Omega$ は $0$ を中心とする球を含み、したがって $\Omega=V$ となるが、$V^*=\emptyset$ である。
$x,y\in\Omega$ に対しては $y\in\Omega^*$ であるから $\langle x,y\rangle>0$ であり、連続性により $\bar\Omega$ 上で $\langle x,y\rangle\ge0$ である。
閉球 $\bar B(x,\delta)$ が $\Omega$ に含まれ、$y\in\bar\Omega\setminus\{0\}$ ならば、
$\langle x-\delta y/\|y\|,y\rangle\ge0$ である。(2) $x\notin\Omega=\Omega^*$ から
$\langle u,x\rangle\le0$ となる $u\in\bar\Omega\setminus\{0\}$ が得られ、(1) により $\langle u,x\rangle\ge0$ である。$x'\in\bar\Omega\setminus\{0\}$ に対し
$\langle x',y+tu\rangle\ge\langle x',y\rangle>0$ であるから $y+tu\in\Omega^*=\Omega$ である。(3) $\bar B(y,\varepsilon)\subseteq\bar\Omega$ ならば、
$x\in\bar\Omega\setminus\{0\}$ に対し $\langle x,y-\varepsilon x/\|x\|\rangle\ge0$ であるから $\langle x,y\rangle\ge\varepsilon\|x\|>0$ であり、
$y\in\Omega^*=\Omega$ である。(4) $g\in G(\Omega)$ に対し、$g$ は $\bar\Omega\setminus\{0\}$ をそれ自身の上に写すので、$y\in\Omega$ と
$x\in\bar\Omega\setminus\{0\}$ に対し $\langle x,g^Ty\rangle=\langle gx,y\rangle>0$ である：$g^T\Omega\subseteq\Omega$。$g^{-1}$ に対しても同様にして
$g^T\Omega=\Omega$ を得る。$g_k\in G(\Omega)$ が $g\in GL(V)$ に収束するならば、$g\Omega\subseteq\bar\Omega$ は開であるから、(3) により
$g\Omega\subseteq\Omega$ である。同様に $g^{-1}=\lim g_k^{-1}$ は $\Omega$ を $\Omega$ の中に写すので、$g\Omega=\Omega$ である。$\square$

$\Omega$ の **特性関数** は $\varphi(x)=\int_{\Omega^*}e^{-\langle x,y\rangle}dy=\int_\Omega e^{-\langle x,y\rangle}dy$ である。
$f=\log\varphi$ とおき、$H(x)$ を $\langle H(x)u,v\rangle=D^2f(x)[u,v]$ を満たす対称作用素とし、
$x^*=-\nabla f(x)$、すなわち $Df(x)[u]=-\langle x^*,u\rangle$ とする。

**補題 F.9（特性関数）.** $\Omega$ を自己双対とする。

1. $\varphi$ は $\Omega$ 上で有限、正、かつ $C^\infty$ であり、積分記号下で微分できる。
2. $t>0$ に対し $\varphi(tx)=t^{-N}\varphi(x)$ であり、$g\in G(\Omega)$ に対し $\varphi(gx)=|\det g|^{-1}\varphi(x)$ である。
3. $x\in\Omega$ に対し $H(x)>0$ かつ $x^*\in\Omega$ である。
4. $H(x)x=x^*$ かつ $D^3f(x)[x,u,v]=-2\langle H(x)u,v\rangle$ である。$g\in G(\Omega)$ に対し $(gx)^*=g^{-T}x^*$ かつ
   $H(gx)=g^{-T}H(x)g^{-1}$ である。
5. ある定数 $c>0$ に対し $\varphi(x)\ge c\|x\|^{-N}$ であり、$x_k\in\Omega$ が
   $\bar\Omega\setminus\Omega$ の点に収束するならば $\varphi(x_k)\to\infty$ である。

*証明.* (1) $x\in\Omega$ とし、$\delta$ を補題 F.8(1) のとおりとする。$\|x'-x\|<\delta/2$ と $y\in\Omega$ に対し
$\langle x',y\rangle\ge\frac\delta2\|y\|$ であるから、被積分関数の $x'$ についての $k$ 階導関数
$(-1)^k\prod_{i}\langle u_i,y\rangle e^{-\langle x',y\rangle}$ は $\prod_i\|u_i\|\,\|y\|^ke^{-\delta\|y\|/2}$ で抑えられ、これは
可積分である（$\|y\|\ge N^{-1/2}\sum_i|y_i|$ とフビニの定理）。優収束定理により積分記号下での微分が
許される。(2) $\Omega$ を $\Omega$ の上に写す置換 $y=y'/t$ と、補題
F.8(4) により $\Omega$ を $\Omega$ の上に写す置換 $y'=g^Ty$ を行う。後者では $\langle gx,y\rangle=\langle x,y'\rangle$ かつ $dy=|\det g|^{-1}dy'$ である。
(3) $\mu_x$ を $\Omega$ 上の確率測度 $e^{-\langle x,y\rangle}dy/\varphi(x)$ とする。(1) により
$Df(x)[u]=-\mathbb E\langle u,y\rangle$ であり、$D^2f(x)[u,u]=\mathbb E\langle u,y\rangle^2-(\mathbb E\langle u,y\rangle)^2$ は
$\langle u,y\rangle$ の分散である。これが零になるのは、$\langle u,y\rangle$ が $\mu_x$-ほとんどいたるところで定数 $\gamma$ に等しい、
すなわち $\Omega$ 上でルベーグ測度についてほとんどいたるところそうである（$\mu_x$ の密度は正である）場合に限る。集合
$\{y\in\Omega:\langle u,y\rangle\ne\gamma\}$ は開かつ零集合であるから空であり、したがって線形関数 $\langle u,\cdot\rangle$ はある球の上で
定数となり、これは $u=0$ を強いる。また $x^*=\mathbb E\,y$ であるから、$z\in\bar\Omega\setminus\{0\}$ に対し $\langle z,x^*\rangle=\mathbb E\langle z,y\rangle>0$ である：
$x^*\in\Omega^*=\Omega$。
(4) (2) により $f(tx)=f(x)-N\log t$ であるから、$\Omega$ 上で $Df(x)[x]=-N$ である。$x$ について $u$ 方向に微分すると
$D^2f(x)[x,u]+Df(x)[u]=0$、すなわち $\langle H(x)x,u\rangle=\langle x^*,u\rangle$ であり、さらに $v$ 方向に微分すると
$D^3f(x)[x,u,v]+2D^2f(x)[u,v]=0$ である。$f(gx)=f(x)-\log|\det g|$ を1回および2回微分すると
$\langle(gx)^*,gu\rangle=\langle x^*,u\rangle$ と $\langle H(gx)gu,gv\rangle=\langle H(x)u,v\rangle$ が得られる。
(5) $\|y\|\le1/\|x\|$ に対し $\langle x,y\rangle\le\|x\|\|y\|\le1$ であるから、$\varphi(x)\ge e^{-1}\operatorname{vol}(\Omega\cap B(0,1/\|x\|))=c\|x\|^{-N}$
である。ただし $c=e^{-1}\operatorname{vol}(\Omega\cap B(0,1))>0$ である。$x_k\to x_0\in\bar\Omega\setminus\Omega$ とする。ファトゥの補題により
$\liminf\varphi(x_k)\ge\int_\Omega e^{-\langle x_0,y\rangle}dy$ である。$u$ を補題 F.8(2) のとおりにとり、球 $B(y_1,r)\subseteq\Omega$
と互いに交わらない球 $B_l=B(y_1+2lru/\|u\|,r)$（$l\ge0$）をとる。これらは補題 F.8(2) により $\Omega$ に含まれ、$B_l$ 上で
$\langle x_0,y\rangle=\langle x_0,y-2lru/\|u\|\rangle\le\|x_0\|(\|y_1\|+r)$ である。よって積分は無限大である。$\square$

**補題 F.10（1パラメータ群）.** $G$ を $GL(V)$ の閉部分群とし、
$\mathfrak g=\{Z:e^{tZ}\in G\text{（すべての }t\in\mathbb R\text{ に対し）}\}$ とする。

1. $c:[0,s_0)\to G$ が $s\to0^+$ のとき $c(s)=I+sZ+o(s)$ を満たすならば、$Z\in\mathfrak g$ である。
2. $\mathfrak g$ はリー環である：$[Z,W]=ZW-WZ$ について閉じた線形部分空間である。

*証明.* (1) $t>0$ を固定し、$A_n=c(t/n)$、$B_n=e^{tZ/n}$ とおく。すると $\|A_n-B_n\|=o(1/n)$ であり、ある $K$ と大きな $n$ に対し
$\|A_n\|,\|B_n\|\le1+K/n$ であるから
$$
\|A_n^n-e^{tZ}\|=\Big\|\sum_{r=0}^{n-1}A_n^r(A_n-B_n)B_n^{n-1-r}\Big\|\le n\,e^K\|A_n-B_n\|\to0 .
$$
$A_n^n\in G$ であり $G$ は $GL(V)$ の中で閉であるから $e^{tZ}\in G$ である。また $e^{-tZ}=(e^{tZ})^{-1}\in G$ である。(2) スカラー倍は明らかである。
(1) を $c(s)=e^{sZ}e^{sW}=I+s(Z+W)+O(s^2)$ に、また $a=\sqrt s$ として $c(s)=e^{aZ}e^{aW}e^{-aZ}e^{-aW}$ に適用する：各
因子を2次まで展開すると、$a$ の1次の項は打ち消し合い、$a^2$ の項は足し合わせて $ZW-WZ$ となるので、
$c(s)=I+s[Z,W]+O(s^{3/2})$ である。$\square$

**命題 F.11（基点）.** $\Omega$ を対称錐、$G=G(\Omega)$ とし、$\mathfrak g$ を補題 F.10 のとおりとする。

1. $x\in\Omega$ に対し $x^{**}=x$ かつ $H(x^*)=H(x)^{-1}$ である。
2. $e^*=e$ を満たす $e\in\Omega$ が存在し、そのような $e$ はすべて $H(e)=I$ を満たす。
3. $x\in\Omega$ と、$ge=x$ を満たす任意の $g\in G$ に対し：$H(x)^{-1}=gg^T\in G$。
4. $Z\in\mathfrak g$ かつ $Ze=0$ ならば $Z^T=-Z$ である。

*証明.* (1) 写像 $x\mapsto x^*$ は $\Omega$ から $\Omega$ への滑らかな写像であり（補題 F.9(3)）、その導関数は $-H(x)$ である。
$g^{-T}\in G$（補題 F.8(4)）であるから、補題 F.9(2,4) により
$f(gx)+f((gx)^*)=f(x)-\log|\det g|+f(g^{-T}x^*)=f(x)+f(x^*)$ となり、$F(x)=f(x)+f(x^*)$ は軌道
$G x=\Omega$ 上で定数である。よって $0=DF(x)[u]=-\langle x^*,u\rangle+\langle x^{**},H(x)u\rangle$、すなわち $H(x)x^{**}=x^*=H(x)x$
（補題 F.9(4)）であり、$H(x)$ は可逆であるから $x^{**}=x$ である。$x^{**}=x$ を微分すると $H(x^*)H(x)=I$ が得られる。
(2) $F_2(x)=f(x)+\frac12\|x\|^2$ は $F_2(x)\ge\log c-N\log\|x\|+\frac12\|x\|^2$ を満たし（補題 F.9(5)）、これは
$\|x\|\to\infty$ のとき $\infty$ に発散する。よって $S=\{x\in\Omega:F_2(x)\le F_2(x_1)\}$（$x_1\in\Omega$）は有界であり、$V$ の中で閉でもある。
$\bar\Omega\setminus\Omega$ 内の極限点があれば $f$ は $S$ 上で非有界になるからである（補題 F.9(5)）。連続な $F_2$ はコンパクト集合 $S$ 上で、
開集合 $\Omega$ 内のある点 $e$ で最小値をとり、そこで $0=\nabla F_2(e)=-e^*+e$ である。すると (1) により $H(e)=H(e^*)=H(e)^{-1}$
であるから $H(e)^2=I$ であり、$H(e)>0$ であるから $H(e)=I$ である。
(3) 補題 F.9(4) と (2) により $H(x)=H(ge)=g^{-T}g^{-1}$ であり、補題 F.8(4) により $gg^T\in G$ である。
(4) $e^{tZ}\in G$ は $e$ を固定するので $I=H(e)=e^{-tZ^T}H(e)e^{-tZ}=e^{-tZ^T}e^{-tZ}$ である。$t=0$ で微分すればよい。$\square$

**定理 F.12（ケッヒャー–ヴィンベルク）.** $\Omega$ をユークリッド空間 $V$ 内の空でない開凸錐とする。$\Omega$ が
対称錐であるための必要十分条件は、$(V,\circ,\langle\cdot,\cdot\rangle)$ が
ユークリッド的ジョルダン代数となり、$\Omega$ がその平方の錐の内部となるような積 $\circ$ が $V$ 上に存在することである。対称錐に対しては、
$e^*=e$ を満たす任意の $e\in\Omega$ を単位元とし、積を
$$
\langle u\circ v,w\rangle=-\tfrac12\,D^3\log\varphi(e)[u,v,w].
$$
とすることができる。

*証明.* 「ならば」の向き（積の存在から対称性）は定理 F.5(2) と F.6(2) である。逆に $\Omega$ を対称錐とし、$e$ を命題 F.11(2) のとおりとし、
$\langle L(u)v,w\rangle=-\frac12D^3f(e)[u,v,w]$ により $L(u)$ を定め、$u\circ v=L(u)v$ とおく。

*ステップ1：可換性、対称性、単位元.* $D^3f(e)$ は三つの引数について対称であるから $u\circ v=v\circ u$ かつ
$L(u)^T=L(u)$ である。補題 F.9(4) と $H(e)=I$ により $\langle e\circ v,w\rangle=\langle H(e)v,w\rangle=\langle v,w\rangle$ である：$e$ は単位元である。

*ステップ2：$L(u)\in\mathfrak g$.* $R(x)=H(x)^{-1}$ は $G$ に属し（命題 F.11(3)）、滑らかで、$R(e)=I$ である。その
導関数は $DR(e)[u]=-H(e)^{-1}DH(e)[u]H(e)^{-1}=-DH(e)[u]=2L(u)$ である。$\langle DH(e)[u]v,w\rangle=D^3f(e)[u,v,w]$ だからである。
$c(s)=R(e+su)$（$e+su\in\Omega$ となるよう $s\ge0$ は小さくとる）と閉部分群 $G$（補題 F.8(4)）に対する補題 F.10(1) により
$2L(u)\in\mathfrak g$、したがって $L(u)\in\mathfrak g$ である。

*ステップ3：一意性.* $Y\in\mathfrak g$ が対称ならば $Y=L(Ye)$ である：対称な $Y-L(Ye)\in\mathfrak g$ は $e$ を零に写す
（$L(w)e=w\circ e=w$）ので、命題 F.11(4) により反対称でもあり、したがって $0$ である。

*ステップ4：内部導分.* $D=[L(x),L(y)]$ とする。これは $\mathfrak g$ に属し（補題 F.10(2)）、$D^T=-D$ であり、
$De=x\circ y-y\circ x=0$ である。$z\in V$ に対し、$[D,L(z)]\in\mathfrak g$ は対称であり $[D,L(z)]e=Dz-L(z)De=Dz$ であるから、
ステップ3により $[D,L(z)]=L(Dz)$ である。これを $w$ に適用すると $D(z\circ w)=Dz\circ w+z\circ Dw$ となる。

*ステップ5：ジョルダン恒等式.* $D=[L(y),L(x)]$ に対し $Dx=y\circ x^2-x\circ(x\circ y)$ および
$D(x^2)=y\circ(x\circ x^2)-x\circ(y\circ x^2)$ である。ステップ4により $D(x^2)=2x\circ Dx$、すなわちすべての $y$ に対し
$y\circ(x\circ x^2)=3x\circ(x^2\circ y)-2x\circ(x\circ(x\circ y))$ である：
$$
L(x\circ x^2)=3L(x)L(x^2)-2L(x)^3 .
$$
すべての $L$ は対称であるから、転置をとると $L(x\circ x^2)=3L(x^2)L(x)-2L(x)^3$ となり、$[L(x),L(x^2)]=0$ である。よって $V$ は
ユークリッド的ジョルダン代数である。$\Omega_J$ をその平方の錐の内部とする。

*ステップ6：$\Omega=\Omega_J$.* $b\in V$ に対し、ステップ2により $e^{L(b)}\in G$ であるから $\exp b=e^{L(b)}e\in\Omega$ である（定理 F.6(1)）。
定理 F.6(2) により $\Omega_J\subseteq\Omega$ である。両方の錐は $\langle\cdot,\cdot\rangle$ について自己双対である（定理 F.5(2)）から、
$\Omega=\Omega^*\subseteq\Omega_J^*=\Omega_J$ である。$\square$

この証明は J. Faraut, A. Korányi, *Analysis on Symmetric Cones* (Oxford 1994) の Chapters I–III の構成に従っている。
彼らが $G(\Omega)$ のリー群構造を用いるところで、補題 F.10 と命題 F.11(3) が必要な
1パラメータ群を直接に作り出す。この定理は M. Koecher (1957) と E. B. Vinberg (1960) によって証明された。

**例 F.13（半正定値計画の錐）.** $V=\operatorname{Herm}_n(\mathbb R)$、$\langle x,y\rangle=\operatorname{Tr}(xy)$、
$N=n(n+1)/2$ とし、$\Omega$ を正定値行列全体とする。これは定理 F.5 と F.6 により対称錐である。このとき
$\varphi(x)=\varphi(1)\det(x)^{-(n+1)/2}$、$\kappa=\sqrt{(n+1)/2}$ として $e=\kappa1$ であり、定理 F.12 の積は
$u\circ_\kappa v=\frac1{2\kappa}(uv+vu)$ である。よって $x\mapsto x/\kappa$ は通常のジョルダン代数の上への同型である。

*証明.* 対称行列の固有射影によるスペクトル分解 $x=\sum_i\lambda_ip_i$ は
定理 F.3 の分解である（一意性。$p_i$ たちは可換であるから $p_i\circ p_j=p_ip_j$）。よって定理 F.5(2) により、平方の錐の
内部は確かに正定値行列全体である。$a\in\Omega$ に対し、$P(a)x=axa$
（命題 13.2(1)）は $G(\Omega)$ に属する（定理 F.6(3)）。$a$ の正規直交固有基底において、
行列 $E_{ii}$ と $(E_{ij}+E_{ji})/\sqrt2$（$i<j$）は、固有値 $a_ia_j$ をもつ $P(a)$ の固有ベクトルの正規直交基底をなすので、
$\det P(a)=\prod_{i\le j}a_ia_j=\det(a)^{n+1}$ である（各 $a_i$ は $n+1$ 回現れる）。補題 F.9(2) により
$\varphi(a^2)=\varphi(P(a)1)=\det(a)^{-(n+1)}\varphi(1)$ であり、すべての $x\in\Omega$ は $a=\sqrt x$ により $a^2$ と書ける。
$D\log\det(x)[u]=\operatorname{Tr}(x^{-1}u)$（$\det(x+su)=\det x\,\det(1+sx^{-1}u)=\det x\,(1+s\operatorname{Tr}(x^{-1}u)+O(s^2))$ から）
と $D(x^{-1})[u]=-x^{-1}ux^{-1}$ を用いると、
$x^*=\kappa^2x^{-1}$ である。よって $e^*=e$ は、$e$ が正定値で $e^2=\kappa^2 1$ であること、すなわち $e=\kappa1$ を意味する。$D^2f(x)[v,w]=\kappa^2\operatorname{Tr}(x^{-1}vx^{-1}w)$ であり、
$D^3f(e)[u,v,w]=-\kappa^2\kappa^{-3}\operatorname{Tr}(uvw+vuw)=-\frac2\kappa\operatorname{Tr}((u\circ v)w)$ である。$\square$

コードは $\operatorname{Herm}_3(\mathbb R)$ の正規直交座標において $f(x)=-\frac{n+1}2\log\det x$ を自動微分することで $e$ における
$D^3f$ を計算し、命題 F.11 と定理 F.12 の主張を検証する：

```python
from torch.autograd.functional import jacobian, hessian

n = 3
pairs = [(i, j) for i in range(n) for j in range(i, n)]
Nd = len(pairs)                                                    # N = 6
Bs = torch.zeros(Nd, n, n, dtype=f64)                              # orthonormal basis for Tr(xy)
for k, (i, j) in enumerate(pairs):
    Bs[k, i, j] = Bs[k, j, i] = 1.0 if i == j else 1 / math.sqrt(2)
mat = lambda v: torch.einsum('k,kij->ij', v, Bs)
vec = lambda A: torch.einsum('kij,ij->k', Bs, A)
f = lambda v: -(n + 1) / 2 * torch.logdet(mat(v))                  # log phi up to a constant (Example F.13)
kappa = math.sqrt((n + 1) / 2)
e = vec(kappa * torch.eye(n, dtype=f64))
assert close(-jacobian(f, e), e)                                   # e* = e
assert close(hessian(f, e), torch.eye(Nd, dtype=f64))             # H(e) = I (Proposition F.11(2))
T = jacobian(lambda v: hessian(f, v, create_graph=True), e)        # T[a, b, c] = D^3 f(e)[B_a, B_b, B_c]
prod = lambda u, v: torch.einsum('a,b,abc->c', u, v, -0.5 * T)     # Theorem F.12
u, v = torch.randn(2, Nd, dtype=f64)
assert close(prod(e, v), v)                                        # unit
assert close(prod(u, v), vec(mat(u) @ mat(v) + mat(v) @ mat(u)) / (2 * kappa))
u2 = prod(u, u)
assert close(prod(u, prod(u2, v)), prod(u2, prod(u, v)))           # Jordan identity
```

## F.4 順序同型：有限次元におけるカディソンの定理

$\mathcal M=\bigoplus_cM_{k_c}(\mathbb C)$ と $\mathcal N=\bigoplus_dM_{k'_d}(\mathbb C)$ をブロック対角行列として考える。
$x\ge0$ は $x=x^*$ で固有値が非負であること（同値なことに、すべてのブロックが半正定値であること）を意味し、
$[0,1]=\{a\in\mathcal M_{sa}:0\le a\le1\}$ とする。**射影** とは $p=p^*=p^2$ のことである。線形写像 $\Phi$ は、
$x\ge0$ ならば $\Phi(x)\ge0$ となるとき **正値** であるという。**ジョルダン $*$-同型** とは、$\Phi(x^*)=\Phi(x)^*$ と
$\Phi(x\circ y)=\Phi(x)\circ\Phi(y)$ を満たす線形全単射のことである（13.5節）。凸集合 $K$ の点 $a$ は、相異なる二点 $b,b'\in K$ の
中点 $\frac12(b+b')$ でないとき **端点** であるという（同値なことに、$b\ne b'$
かつ $0<t<1$ として $tb+(1-t)b'$ の形でないこと：そのような $a$ は、小さな $\varepsilon>0$ に対し $a\pm\varepsilon(b-b')\in K$ の中点である）。

**補題 F.14.** $\mathcal M$ において：

1. 凸集合 $[0,1]$ の端点は、ちょうど射影である。
2. 射影 $p,q$ に対し、$pq=0$、$p+q\le1$、$p\circ q=0$ は同値である。
3. 射影 $p,q$ に対し、$p\le q$、$pq=p$、$p\circ q=p$ は同値である。極小な非零射影は（あるブロック内の）階数1の
   射影である。
4. 射影 $p$ が中心的であるための必要十分条件は、すべての $x\in\mathcal M$ に対し $pxp+(1-p)x(1-p)=x$ となることである。中心射影は
   ブロック単位 $P_c$ の和であり、極小な非零中心射影は $P_c$ たちである。

*証明.* (1) $p$ を射影とし、$a,b\in[0,1]$ により $p=\frac12(a+b)$ とする。$\xi\in\ker p$ に対し
$0=\langle\xi,a\xi\rangle+\langle\xi,b\xi\rangle$ で、両項とも $\ge0$ であるから $\|a^{1/2}\xi\|^2=0$、よって $a\xi=0$ である。$p$ の値域の
$\xi$ に対しては、$1-a,1-b$ について同じ議論をすれば $a\xi=\xi$ を得る。よって $a=p$ である。逆に、$a\in[0,1]$ が射影でなければ、
固有値 $\lambda\in(0,1)$ をもつ。スペクトル射影 $q_\lambda\in\mathcal M$（$a$ の多項式）と
$\varepsilon=\min(\lambda,1-\lambda)$ により $a\pm\varepsilon q_\lambda\in[0,1]$ であり、$a$ はそれらの中点である。
(2) $pq=0$ ならば $qp=(pq)^*=0$ であり、$p+q$ は射影で、$1-(p+q)=(1-p-q)^2\ge0$ である。また $p\circ q=0$ である。
$p+q\le1$ ならば、すべての $\xi$ に対し $\|qp\xi\|^2=\langle p\xi,qp\xi\rangle\le\langle p\xi,(1-p)p\xi\rangle=0$ であるから $qp=0=pq$ である。
$p\circ q=0$ ならば $pq=-qp$ であるから $pq=p(pq)=-pqp=(qp)p=qp$ であり、$pq=-pq=0$ である。
(3) $pq=p$ ならば $qp=p$、$p\circ q=p$、そして $q-p=(q-p)^2\ge0$ である。$p\circ q=p$ ならば、$pq+qp=2p$ に $p$ を
左右から掛けると $pq+pqp=2p=pqp+qp$ となるので $pq=qp=p$ である。$p\le q$ ならば $p+(1-q)\le1$ であり、(2) により
$p(1-q)=0$ である。階数 $\ge2$ の射影は、その値域内の直線であって一つのブロックに含まれるものへの射影 $\xi\xi^*$ を
真に上から抑える（$\xi$ は $p$ のあるブロックの値域内の単位ベクトル）。階数1の射影 $p$ は一つのブロックに含まれ、
$0\ne q\le p$ ならば $q=pq$ の値域は $p$ の値域に含まれるので $q=p$ である。
(4) $x=pxp+px(1-p)+(1-p)xp+(1-p)x(1-p)$ であり、$px(1-p)+(1-p)xp=0$ であることは両項がともに零であること（$p$ と $1-p$ を掛ける）、
すなわち $px=pxp=xp$ であることと同値である。中心元は各ブロックのすべての行列単位と可換であるから、ブロック
$c$ 上でスカラー $z_c$ であり、射影に対しては $z_c\in\{0,1\}$ である。$\square$

**定理 F.15（カディソン、有限次元）.** $\Phi:\mathcal M\to\mathcal N$ を
$\Phi(1)=1$ を満たす線形全単射で、$\Phi$ と $\Phi^{-1}$ がともに正値であるものとする。このとき $\Phi$ はジョルダン $*$-同型である。逆に、すべての
ジョルダン $*$-同型は単位的であり、それとその逆写像はともに正値である。

*証明.* すべての $x\in\mathcal M_{sa}$ は $x_\pm\ge0$ により $x_+-x_-$ と書けるので、$\Phi(x)$ は自己共役である。
$h,k$ を自己共役として $x=h+ik$ と書けば $\Phi(x^*)=\Phi(h)-i\Phi(k)=\Phi(x)^*$ である。$\Phi^{-1}$ についても同じことが成り立つので、$\Phi$ は
$x\ge0\iff\Phi(x)\ge0$ を満たす実線形全単射 $\mathcal M_{sa}\to\mathcal N_{sa}$ である。$\Phi(1)=1$ であるから、$\Phi$ は $[0,1]$ を $[0,1]$ の上に写す。
凸集合の間の線形全単射は端点を端点の上に写す（$\Phi(a)=\frac12(b+b')$ ならば
$a=\frac12(\Phi^{-1}b+\Phi^{-1}b')$）ので、補題 F.14(1) により $\Phi$ は射影を射影の上に写し、補題
F.14(2) により $pq=0\iff p+q\le1\iff\Phi(p)+\Phi(q)\le1\iff\Phi(p)\Phi(q)=0$ である。$x=\sum_i\lambda_ip_i\in\mathcal M_{sa}$（スペクトル
分解、$i\ne j$ に対し $p_ip_j=0$）に対し $\Phi(x)^2=\sum_i\lambda_i^2\Phi(p_i)=\Phi(x^2)$ である。偏極
$x\circ y=\frac12((x+y)^2-x^2-y^2)$ により $\mathcal M_{sa}$ 上で $\Phi(x\circ y)=\Phi(x)\circ\Phi(y)$ が得られ、複素
双線形性により $\mathcal M=\mathcal M_{sa}+i\mathcal M_{sa}$ 上でも成り立つ。

逆に、すべての $y$ に対し $\Phi(1)\circ\Phi(y)=\Phi(y)$ であり $\Phi$ は全射であるから、$\Phi(1)$ は $\circ$ の単位元であり、
$\Phi(1)=\Phi(1)\circ1=1$ である。すべての $x\ge0$ は自己共役な $h=\sqrt x$ により $h^2$ と書け（定理 1.6）、
$\Phi(h)$ は自己共役であるから $\Phi(h^2)=\Phi(h)^2\ge0$ である。$\Phi^{-1}$ もジョルダン $*$-同型である。$\square$

カディソンはこの定理をすべての単位的 C\*-環に対して証明した（R. V. Kadison, *A generalized Schwarz inequality and algebraic
invariants for operator algebras*, Ann. of Math. 56 (1952) 494–503）。その版は本講義では用いない。

```python
def blocks_rand(ks, psd=False):
    out = []
    for k in ks:
        A = torch.randn(k, k, dtype=dt)
        out.append(A @ A.mH if psd else A + A.mH)
    return out

ks = (2, 2, 3)                                                     # M = M_2 (+) M_2 (+) M_3
Us = [torch.linalg.qr(torch.randn(k, k, dtype=dt))[0] for k in ks]
def Phi(x):                                                        # swap the M_2 blocks, transpose two blocks
    x1, x2, x3 = x
    return [Us[0] @ x2.mT @ Us[0].mH, Us[1] @ x1 @ Us[1].mH, Us[2] @ x3.mT @ Us[2].mH]
def Phi_inv(y):
    y1, y2, y3 = y
    return [Us[1].mH @ y2 @ Us[1], (Us[0].mH @ y1 @ Us[0]).mT, (Us[2].mH @ y3 @ Us[2]).mT]
is_psd = lambda x: all(torch.linalg.eigvalsh(b)[0] > -1e-10 for b in x)
one = [torch.eye(k, dtype=dt) for k in ks]
assert all(close(a, b) for a, b in zip(Phi(one), one))            # unital
for _ in range(50):
    p = blocks_rand(ks, psd=True)
    assert is_psd(Phi(p)) and is_psd(Phi_inv(p))                   # positive in both directions
x, y = blocks_rand(ks), blocks_rand(ks)
assert all(close(a, b) for a, b in zip(Phi_inv(Phi(x)), x))
jordan = lambda a, b: [0.5 * (s @ t + t @ s) for s, t in zip(a, b)]
assert all(close(a, b) for a, b in zip(Phi(jordan(x, y)), jordan(Phi(x), Phi(y))))   # Theorem F.15
assert not close(Phi([s @ t for s, t in zip(x, y)])[0], (Phi(x)[0] @ Phi(y)[0]))      # not multiplicative
```

## F.5 $\bigoplus_cM_{k_c}$ のジョルダン $*$-同型

**補題 F.16.** $\Phi:\mathcal M\to\mathcal N$ をジョルダン $*$-同型とする。このとき $\Phi(1)=1$、
$\Phi(xyx)=\Phi(x)\Phi(y)\Phi(x)$ であり、$\Phi$ は $\mathcal M$ の射影を $\mathcal N$ の射影の上に全単射に写し、
$pq=0$、$p\le q$、極小性、中心性、および中心射影の中での極小性を保存し、かつ反映する。

*証明.* 定理 F.15 と同様に $\Phi(1)=1$ である。すべての結合代数において $xyx=2x\circ(x\circ y)-(x\circ x)\circ y$ であり、これが
二つ目の主張を与える。$p=p^*=p\circ p$ は保たれ、補題 F.14(2,3,4) の特徴づけも保たれる。それらは
$\circ$、三重積 $xyx$、$1$、そして環全体にわたる量化子だけを用いるからである。$\Phi^{-1}$ もジョルダン
$*$-同型である。$\square$

**定理 F.17（$M_n$ のジョルダン $*$-自己同型）.** すべてのジョルダン $*$-同型 $\Phi:M_n(\mathbb C)\to M_n(\mathbb C)$
は、ユニタリ $U$ により $x\mapsto UxU^*$ または $x\mapsto Ux^TU^*$ と表される。$n\ge2$ のとき、前者の写像は乗法的、後者は
反乗法的であるから、両方の形をもつ写像はない。

*証明.* $n=1$ のときは $\Phi=\mathrm{id}$ である。$n\ge2$ とし、$E_{ij}$ を行列単位とする。補題 F.16 により $\Phi(E_{ii})$ は
互いに直交する階数1の射影であり、正規直交な $\xi_i$ により $\Phi(E_{ii})=\xi_i\xi_i^*$ と書ける。
$We_i=\xi_i$ を満たすユニタリ $W$ により、ジョルダン $*$-同型 $\Phi_1(x)=W^*\Phi(x)W$ は $\Phi_1(E_{ii})=E_{ii}$ を満たす。

$i\ne j$ とする。空間 $S_{ij}=\operatorname{span}\{E_{ij},E_{ji}\}$ は、
$p=E_{ii}+E_{jj}$ として $\{x:pxp=x,\ E_{ii}xE_{ii}=E_{jj}xE_{jj}=0\}$ に等しいので、補題 F.16 により $\Phi_1(S_{ij})\subseteq S_{ij}$ である。$\Phi_1(E_{ij})=\alpha E_{ij}+\beta E_{ji}$ と書く。
$E_{ij}\circ E_{ij}=0$ から $0=\Phi_1(E_{ij})^2=\alpha\beta(E_{ii}+E_{jj})$ となるので、$\alpha,\beta$ のちょうど一方が零でない
（$\Phi_1$ は単射である）。$\Phi_1(E_{ij})\in\mathbb CE_{ij}$ ならば $\varepsilon(i,j)=+$、$\Phi_1(E_{ij})\in\mathbb CE_{ji}$ ならば $-$ とおく。
$\Phi_1(E_{ji})=\Phi_1(E_{ij})^*$ であるから $\varepsilon(j,i)=\varepsilon(i,j)$ である。相異なる $i,j,k$ に対し
$E_{ij}\circ E_{jk}=\frac12E_{ik}\ne0$ である。$\varepsilon(i,j)=+$ かつ $\varepsilon(j,k)=-$ ならば
$\Phi_1(E_{ij})\circ\Phi_1(E_{jk})\in\mathbb C\,E_{ij}\circ E_{kj}=\{0\}$ となり矛盾する。$(-,+)$ についても同様である。よって
$i,j,k$ が相異なるときは常に $\varepsilon(i,j)=\varepsilon(j,k)$ であり、$\varepsilon$ は定数である：二つの組は添字を共有するか、
第三の組を通じてつながっているからである。$\varepsilon=-$ ならば $\Phi_1$ を $\Phi_1\circ T$（$T$ は転置で、命題 13.9 により
ジョルダン $*$-自己同型）で置き換える。これは $\varepsilon=+$ をもち、なお $E_{ii}$ を固定する。

これですべての $i,j$ に対し $\Phi_1(E_{ij})=\alpha_{ij}E_{ij}$ であり、$\alpha_{ii}=1$、$\alpha_{ji}=\bar\alpha_{ij}$ である。
$E_{ij}\circ E_{ji}=\frac12(E_{ii}+E_{jj})$ から $|\alpha_{ij}|^2=1$ が、$E_{ij}\circ E_{jk}=\frac12E_{ik}$（$i,j,k$ は相異なる）から
$\alpha_{ij}\alpha_{jk}=\alpha_{ik}$ が得られる。$u_i=\alpha_{i1}$ とおくと $\alpha_{ij}=\alpha_{i1}\alpha_{1j}=u_i\bar u_j$ となる（
$i=j$、$i=1$ または $j=1$ の場合も成り立つ）。よって $D=\operatorname{diag}(u_i)$ として $\Phi_1(x)=DxD^*$ であり、$\Phi(x)=WDxD^*W^*$ または
$\Phi(x)=WDx^TD^*W^*$ である。両方の形が一致するならば $T$ は乗法的となるが、$(E_{12}E_{21})^T=E_{11}\ne E_{22}=E_{12}^TE_{21}^T$ である。
$\square$

これは、13.5節で引用した N. Jacobson, C. E. Rickart (Trans. Amer. Math. Soc. 69 (1950) 479–502) と I. N.
Herstein (Trans. Amer. Math. Soc. 81 (1956) 331–341) の結果を、有限次元の場合にすべての
$n$ について証明するものである。命題 13.11 は $n=2$ の場合である。$\operatorname{Herm}_n(\mathbb C)$ のジョルダン自己同型は
$\Phi(x+iy)=\Phi(x)+i\Phi(y)$ によって $M_n$ のジョルダン $*$-自己同型に拡張されるので、これも二つの形のいずれかである。

**定理 F.18（直和）.** $\Phi:\bigoplus_{c=1}^CM_{k_c}\to\bigoplus_{d=1}^{C'}M_{k'_d}$ をジョルダン
$*$-同型とする。このとき $C'=C$ であり、$k'_{\sigma(c)}=k_c$ を満たす全単射 $\sigma$、ユニタリ $U_c\in M_{k_c}$、および
ブロックの集合 $S$ が存在して
$$
\Phi(x)_{\sigma(c)}=U_cx_cU_c^*\ (c\notin S),\qquad \Phi(x)_{\sigma(c)}=U_cx_c^TU_c^*\ (c\in S).
$$
となる。逆に、そのような写像はすべてジョルダン $*$-同型であり、$S$ は $k_c=1$ のブロックを除いて一意である。

*証明.* 補題 F.16 と補題 F.14(4) により、極小中心射影の間の全単射 $\sigma$ によって $\Phi(P_c)=P'_{\sigma(c)}$ となるので、$C=C'$ である。
$x=P_cxP_c$ ならば $\Phi(x)=P'_{\sigma(c)}\Phi(x)P'_{\sigma(c)}$ であるから、$\Phi$ はブロック $c$ をブロック
$\sigma(c)$ の中に写し、$\Phi^{-1}$ はブロック $\sigma(c)$ をブロック $c$ の中に写す。よって $k_c^2=k'^2_{\sigma(c)}$ であり、制限は
ジョルダン $*$-同型 $M_{k_c}\to M_{k_c}$ であって、定理 F.17 により二つの形のいずれかである。$x=\sum_cP_cxP_c$ であるから、これで
$\Phi$ が定まる。逆はブロックごとに成り立ち（命題 13.9）、$S$ の一意性は定理 F.17 から従う。$\square$

**系 F.19（単位的順序同型）.** $\Phi$ と $\Phi^{-1}$ がともに正値であるような単位的線形全単射 $\Phi:\mathcal M\to\mathcal N$ は、
ちょうど定理 F.18 の写像である。そのような $\Phi$ について、$2$-正値であること、
完全正値であること、$S$ が $k_c\ge2$ のブロックを含まないこと、$\Phi$ が $*$-同型であることは、すべて同値である。

*証明.* 最初の文は定理 F.15 と F.18 である。$S$ が $k_c\ge2$ のブロックを含まなければ、すべてのブロック写像は
$x\mapsto U_cxU_c^*$ であり（$k_c=1$ に対しては $x^T=x$）、よって $\Phi$ は完全正値（命題 6.4）かつ乗法的である。
$k_c\ge2$ なる $c\in S$ があり $\Phi$ が $2$-正値ならば、ブロック写像 $\Phi_{\sigma(c)c}(a)=U_ca^TU_c^*$ もそうであり
（命題 6.4）、$a\mapsto U_c^*\Phi_{\sigma(c)c}(a)U_c=a^T$ もそうである。$[A_{ij}]\ge0$ ならば
$[U^*A_{ij}U]=(1\otimes U)^*[A_{ij}](1\otimes U)\ge0$ だからである。これは命題 6.1 に矛盾する。定理 F.17 により、そのような $\Phi$ は
乗法的でもない。$\square$

したがって単位的順序同型は、ユニタリを除いて、どの等しいブロックを入れ替え、どのブロックを
転置するかによって決まり、系 F.19 により、転置されたブロックこそが完全正値性を妨げるものである。コードは $M_3$ のジョルダン
自己同型 $x\mapsto Ux^TU^*$ を構成し、定理 F.17 の証明の手順によって $U$ を（位相を除いて）復元する：

```python
k = 3
E = lambda i, j: torch.zeros(k, k, dtype=dt).index_put_((torch.tensor(i), torch.tensor(j)), torch.tensor(1.0, dtype=dt))
U0 = torch.linalg.qr(torch.randn(k, k, dtype=dt))[0]
Psi = lambda A: U0 @ A.mT @ U0.mH                                  # the unknown Jordan *-automorphism

W = torch.stack([torch.linalg.eigh(Psi(E(i, i)))[1][:, -1] for i in range(k)], dim=1)   # Psi(E_ii) = xi_i xi_i^*
Psi1 = lambda A: W.mH @ Psi(A) @ W
assert all(close(Psi1(E(i, i)), E(i, i)) for i in range(k))
transposed = Psi1(E(0, 1))[1, 0].abs().item() > 0.5               # epsilon = - : Psi1(E_01) is a multiple of E_10
Psi2 = (lambda A: Psi1(A.mT)) if transposed else Psi1
alpha = torch.stack([Psi2(E(i, 0))[i, 0] for i in range(k)])       # u_i = alpha_{i1}
for i in range(k):
    for j in range(k):
        assert close(Psi2(E(i, j)), alpha[i] * alpha[j].conj() * E(i, j))
Ufound = W @ torch.diag(alpha)
X = torch.randn(k, k, dtype=dt)
assert transposed and close(Psi(X), Ufound @ X.mT @ Ufound.mH)     # Theorem F.17, second form
phase = (U0.mH @ Ufound)[0, 0]
assert close(U0.mH @ Ufound, phase * torch.eye(k, dtype=dt)) and abs(abs(phase) - 1) < 1e-12
print("U recovered up to the phase", complex(phase))
```

## 演習問題

1. 座標ごとの積と標準内積をもつ $\mathbb R^N$ が、平方の錐が $\mathbb R^N_{\ge0}$ であるような
   ユークリッド的ジョルダン代数であることを示し、系 F.19 から、
   $\mathbb C^N=\bigoplus_{c=1}^NM_1$ の単位的順序同型は座標の置換であることを導け。
2. （**正値だけでは不十分である。**）$0<t<1$ に対し、$M_n$（$n\ge2$）上で $\Psi_t(x)=(1-t)x+\frac tn\operatorname{Tr}(x)1$ とする。
   $\Psi_t$ が単位的な正値全単射であること、$\Psi_t^{-1}$ が正値でないこと、そして $\Psi_t$ が
   射影を射影に写さないことを示せ。

<details>
<summary>解答</summary>

1. 積は可換かつ結合的であるからジョルダン積であり、単位元は $(1,\dots,1)$ である。$L(x)$ は対角行列
   $\operatorname{diag}(x)$ であり、これは対称である。平方は座標が非負のベクトルである。すべての $k_c=1$ とした
   系 F.19 により $\Phi(x)_{\sigma(c)}=U_cx_cU_c^*=x_c$ である（$M_1$ のユニタリは位相因子である）から、
   $\Phi$ は座標を置換する。逆に、置換は単位的順序同型である。

2. $\Psi_t(1)=(1-t)1+t1=1$ であり、$x\ge0$ に対しては両項が $\ge0$ である。$\operatorname{Tr}\Psi_t(x)=\operatorname{Tr}x$ であるから
   $\Psi_t^{-1}(y)=\frac1{1-t}\big(y-\frac tn\operatorname{Tr}(y)1\big)$ である。階数1の射影 $p$ に対し、$\Psi_t^{-1}(p)$ は
   $\ker p$ 上で固有値 $-\frac{t}{n(1-t)}<0$ をもち、$\Psi_t(p)$ は固有値 $\frac tn\notin\{0,1\}$ をもつ。これは
   定理 F.15 と整合的である：射影を射影に写すことを強いるのは、逆写像の正値性である。

```python
nn, t = 3, 0.25
p = E(0, 0)
Psi_t = lambda A: (1 - t) * A + t / nn * torch.trace(A) * torch.eye(nn, dtype=dt)
Psi_t_inv = lambda A: (A - t / nn * torch.trace(A) * torch.eye(nn, dtype=dt)) / (1 - t)
assert close(Psi_t_inv(Psi_t(X)), X)
assert abs(torch.linalg.eigvalsh(Psi_t_inv(p))[0].item() + t / (nn * (1 - t))) < 1e-12
assert not close(Psi_t(p) @ Psi_t(p), Psi_t(p))                   # not a projection
```

</details>

**戻る：**[第13章](13_jordan_algebras.md)、13.4–13.6節：定理 F.12 は定理 13.8 を証明し、定理 F.15 は
環 $\bigoplus_cM_{k_c}(\mathbb C)$ に対して定理 13.10 を証明し、定理 F.17 は命題 13.11 をすべての
$n$ に拡張し、系 F.7 は命題 13.14 をすべてのユークリッド的ジョルダン代数に拡張する。
