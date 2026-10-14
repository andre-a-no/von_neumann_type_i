[English](../F_cones_and_order.md) · [Español](../es/F_cones_and_order.md) · [Français](../fr/F_cones_and_order.md) · [Deutsch](../de/F_cones_and_order.md) · **中文** · [日本語](../ja/F_cones_and_order.md) · [Русский](../ru/F_cones_and_order.md)

# 附录 F. 对称锥、序同构与若尔当映射

第 13 章通过三个定理把可观测量的序与它们的若尔当积联系起来：Koecher–Vinberg 定理（定理 13.8）、
Kadison 定理（定理 13.10），以及对 $M_n(\mathbb C)$ 的若尔当自同构的描述（命题 13.11 对 $n=2$ 证明了它）。本附录
在完全一般的情形下证明第一个定理，并对本课程的代数 $\bigoplus_cM_{k_c}(\mathbb C)$ 在所有块大小下证明另外两个。
在此过程中，它证明了 13.6 节所需的关于欧几里得若尔当代数的事实：谱定理、平方锥的自对偶性，以及
$P(y)Q\subseteq Q$。

**你将学到**为什么每个欧几里得若尔当代数都具有幂结合性并有谱定理；为什么它的平方锥是自对偶且齐性的；
如何仅从一个自对偶齐性锥出发，通过其特征函数的对数的三阶导数恢复若尔当积；为什么双向保持正性的
保单位线性双射保持平方；以及为什么 $\bigoplus_cM_{k_c}$ 的每个这样的映射都置换相等的块，并在每个块上
以 $UxU^*$ 或 $Ux^TU^*$ 的方式作用。

除线性代数、微积分以及按编号引用的第 1、6、13 章的结果外，F.3 节还用到 Lebesgue 积分（Fatou 引理、
控制收敛与积分号下求导、Fubini–Tonelli 定理、线性变量替换公式）；此外不假定其他知识。附录 E 中借助若尔当标架的
谱定理是下面定理 F.3 的加细；这里没有任何内容依赖于它。
**欧几里得空间**是带有内积 $\langle\cdot,\cdot\rangle$ 的有限维实向量空间 $V\ne\{0\}$；$\|x\|$ 是范数（对 $V$
的元素，单竖线 $|x|$ 将表示定理 F.3 意义下的绝对值，而绝不表示范数），$\|A\|$ 是 $V$ 的线性映射 $A$ 的算子范数，
$A^T$ 是它的伴随（$\langle Ax,y\rangle=\langle x,A^Ty\rangle$），而 $A\ge0$（$A>0$）表示 $A=A^T$ 是半正定的
（正定的）。$e^A=\exp A=\sum_{k\ge0}A^k/k!$。

## F.1 欧几里得若尔当代数

与 13.1 节一样，**若尔当代数**是一个实向量空间，带有双线性交换乘积 $\circ$，满足
$x\circ(x^2\circ y)=x^2\circ(x\circ y)$，即 $[L(x),L(x^2)]=0$，其中 $L(x)y=x\circ y$，$x^2=x\circ x$；并且
$P(x)=2L(x)^2-L(x^2)$。

**定义。** **欧几里得若尔当代数**是具有单位元 $e$（$e\circ x=x$）以及一个内积的有限维若尔当代数 $V$，
该内积使每个 $L(x)$ 都是对称的：$\langle x\circ y,z\rangle=\langle y,x\circ z\rangle$。
它的**平方锥**是 $Q=\{x^2:x\in V\}$。

例子：$M_{sa}$（$\bigoplus_cM_{k_c}(\mathbb C)$ 的自伴部分）与 $\operatorname{Herm}_n(\mathbb R)$，
取 $\langle x,y\rangle=\operatorname{Tr}(xy)$，它是实的且 $\operatorname{Tr}(x^2)=\sum_{ij}|x_{ij}|^2$（引理 13.1；
$\operatorname{Tr}((x\circ y)z)=\frac12\operatorname{Tr}(xyz+yxz)=\operatorname{Re}\operatorname{Tr}(xyz)$，因为
$\operatorname{Tr}(yxz)=\overline{\operatorname{Tr}((yxz)^*)}=\overline{\operatorname{Tr}(zxy)}=\overline{\operatorname{Tr}(xyz)}$，
它在循环置换下不变，并且出于同样的理由在 $x\leftrightarrow z$ 下不变，从而关于
$x,y,z$ 对称）；自旋因子 $V_m$（由引理 13.3 是若尔当代数），取内积
$\langle x,y\rangle=\operatorname{tr}(x\circ y)=2(ab+u\cdot v)$，其中 $x=(a,u)$，$y=(b,v)$，因为对 $z=(c,w)$，
$\operatorname{tr}((x\circ y)\circ z)=2(abc+c\,u\cdot v+a\,v\cdot w+b\,u\cdot w)$ 关于
$x,y,z$ 对称；以及带有逐坐标乘积的 $\mathbb R^N$（习题 1）。

**引理 F.1.** 在每个若尔当代数中，对所有 $a,b,d,x,y$：

1. $[L(a),L(b\circ d)]+[L(b),L(d\circ a)]+[L(d),L(a\circ b)]=0$；
2. $L(x^2\circ y)=L(x^2)L(y)-2L(x)L(y)L(x)+2L(x\circ y)L(x)$。

*证明。* (1) 对 $s,t,r\in\mathbb R$ 与 $x=sa+tb+rd$，$[L(x),L(x^2)]=0$ 是关于 $(s,t,r)$ 的多项式恒等式，所以每个
系数都为零；$str$ 的系数是 (1) 左边的两倍。(2) 在 (1) 中令 $a=b=x$：
$2[L(x),L(x\circ d)]+[L(d),L(x^2)]=0$。作用于 $y$ 得
$2x\circ((x\circ d)\circ y)-2(x\circ d)\circ(x\circ y)+d\circ(x^2\circ y)-x^2\circ(d\circ y)=0$，即对所有 $d$ 有
$\big(2L(x)L(y)L(x)-2L(x\circ y)L(x)+L(x^2\circ y)-L(x^2)L(y)\big)d=0$。$\square$

**命题 F.2（幂结合性）。** 设 $V$ 是具有单位元 $e$ 的若尔当代数，$x\in V$，$x^0=e$，
$x^{k+1}=x\circ x^k$。则对所有 $p,q\ge0$ 有 $x^p\circ x^q=x^{p+q}$，并且每个 $L(x^p)$ 都是两两可交换的算子
$L(x)$ 与 $L(x^2)$ 的多项式。因此对实多项式 $p,q$ 有 $p(x)\circ q(x)=(pq)(x)$，并且各算子
$L(p(x))$ 彼此可交换。

*证明。* 由若尔当恒等式，$L(x)$ 与 $L(x^2)$ 可交换。我们对 $N\ge2$ 归纳证明：(A$_N$) 对 $i+j\le N$ 有
$x^i\circ x^j=x^{i+j}$，以及 (B$_N$) 对 $i\le N$，$L(x^i)$ 是 $L(x),L(x^2)$ 的多项式。两者对 $N=2$ 都是
显然的。假设它们对 $N$ 成立，并设 $i+j=N+1$，$i,j\ge1$（$i=0$ 或 $j=0$ 时无需证明）。
由于 $i\le N$，$L(x^i)$ 与 $L(x)$ 可交换，所以由 (A$_N$)
$$
x^i\circ x^j=L(x^i)L(x)x^{j-1}=L(x)L(x^i)x^{j-1}=x\circ x^N=x^{N+1}.
$$
这就是 (A$_{N+1}$)。对 $y=x^{N-1}$ 应用引理 F.1(2)，并利用 $x^2\circ x^{N-1}=x^{N+1}$，得到
$L(x^{N+1})=L(x^2)L(x^{N-1})-2L(x)L(x^{N-1})L(x)+2L(x^N)L(x)$，这就是 (B$_{N+1}$)。最后一句由
双线性得出。$\square$

**定理 F.3（谱定理）。** 设 $V$ 是欧几里得若尔当代数，$x\in V$。存在唯一的
$k\ge1$、实数 $\lambda_1>\dots>\lambda_k$ 以及非零的 $c_1,\dots,c_k\in V$，使得
$$
c_i\circ c_j=\delta_{ij}c_i,\qquad \textstyle\sum_ic_i=e,\qquad x=\sum_i\lambda_ic_i .
$$
各 $c_i$ 是 $x$ 的多项式，并且对每个实多项式 $p$ 有 $p(x)=\sum_ip(\lambda_i)c_i$。

我们称 $(c_i)$ 为**正交幂等元的完全系**，称 $\lambda_i$ 为 $x$ 的**特征值**，并对
$f:\{\lambda_1,\dots,\lambda_k\}\to\mathbb R$ 令 $f(x)=\sum_if(\lambda_i)c_i$；例如 $|x|$、$x_\pm=\max(\pm x,0)$、
$\sqrt x$（若所有 $\lambda_i\ge0$）、$\exp x$ 以及 $\log x$（若所有 $\lambda_i>0$）。

*证明。* 存在性。$\mathbb R[x]=\{p(x)\}$ 是一个包含 $e$ 的子空间，并且由命题 F.2 满足
$p(x)\circ q(x)=(pq)(x)$。所以 $M=L(x)|_{\mathbb R[x]}$ 把 $\mathbb R[x]$ 映入自身，是对称的，并满足
$p(M)q(x)=p(x)\circ q(x)$。对称算子可对角化且特征值为实数；设 $\lambda_1>\dots>\lambda_k$
是 $M$ 的特征值，$m(t)=\prod_i(t-\lambda_i)$ 是它的极小多项式。则 $m(x)=m(M)e=0$。借助 Lagrange
多项式 $\ell_i(t)=\prod_{j\ne i}\frac{t-\lambda_j}{\lambda_i-\lambda_j}$，多项式
$\ell_i\ell_j-\delta_{ij}\ell_i$、$\sum_i\ell_i-1$ 与 $\sum_i\lambda_i\ell_i-t$ 在每个 $\lambda_l$ 处为零，所以它们都是
$m$ 的倍式；在 $x$ 处取值表明 $c_i=\ell_i(x)$ 具有这三条性质。若 $c_i=0$，则对所有 $q$ 有
$\ell_i(M)q(x)=\ell_i(x)\circ q(x)=0$，所以 $\ell_i(M)=0$，尽管 $\deg\ell_i<\deg m$；因此 $c_i\ne0$。
最后 $p(x)=\sum_ip(\lambda_i)c_i$：它对 $p=1$ 与 $p=t$ 成立，并且 $x\circ\sum_i\lambda_i^rc_i=\sum_i\lambda_i^{r+1}c_i$。

唯一性。设 $x=\sum_j\mu_jd_j$，其中 $\mu_j$ 互不相同，$d_j$ 是由非零正交幂等元组成的完全系。
同上，$p(x)=\sum_jp(\mu_j)d_j$。各 $d_j$ 线性无关（用 $d_i$ 乘 $\sum_ja_jd_j=0$），所以
$p(x)=0$ 当且仅当对所有 $j$ 有 $p(\mu_j)=0$。于是由 $m(x)=0$ 得 $\{\mu_j\}\subseteq\{\lambda_i\}$，而
$\tilde m(t)=\prod_j(t-\mu_j)$ 满足 $\tilde m(M)q(x)=\tilde m(x)\circ q(x)=0$，所以 $m$ 整除 $\tilde m$，从而 $\{\lambda_i\}\subseteq\{\mu_j\}$。
于是 $d_j=\ell_j(x)=c_j$。$\square$

**引理 F.4（Peirce 分解）。** 设 $V$ 是欧几里得若尔当代数，$c$ 是幂等元（$c^2=c$）。

1. $L(c)(2L(c)-I)(L(c)-I)=0$。因此 $V$ 是特征子空间 $V(c,\nu)=\{z:c\circ z=\nu z\}$
   （$\nu\in\{0,\frac12,1\}$）的正交和，并且 $L(c)\ge0$。
2. $V(c,1)\circ V(c,0)=\{0\}$；对 $\nu\in\{0,1\}$ 有 $V(c,\nu)\circ V(c,\nu)\subseteq V(c,\nu)$ 与 $V(c,\nu)\circ V(c,\frac12)\subseteq V(c,\frac12)$；
   $V(c,\frac12)\circ V(c,\frac12)\subseteq V(c,0)+V(c,1)$。
3. $s=2c-e$ 满足 $s^2=e$，并且 $P(s)$ 是 $V$ 的若尔当自同构：它在 $V(c,0)+V(c,1)$ 上为 $+1$，在
   $V(c,\frac12)$ 上为 $-1$，且 $P(s)(z\circ w)=P(s)z\circ P(s)w$。
4. 设 $x=\sum_{i=1}^k\lambda_ic_i$ 如定理 F.3 所示。各 $L(c_i)$ 两两可交换，并且 $V$ 是子空间
   $V_{ii}=\{z:c_i\circ z=z\}$ 与 $V_{ij}=\{z:c_i\circ z=c_j\circ z=\frac12z\}$（$i<j$）的正交和。对
   $a=\sum_la_lc_l$（$a_l\in\mathbb R$），$L(a)$ 在 $V_{ij}$（$i\le j$）上作用为数乘 $\frac12(a_i+a_j)$，而 $P(a)$
   作用为数乘 $a_ia_j$。

*证明。* (1) 在引理 F.1(2) 中取 $x=y=c$，得 $L(c)=L(c)^2-2L(c)^3+2L(c)^2$，即 $2L(c)^3-3L(c)^2+L(c)=0$，这就是
所述的因式分解。所以对称算子 $L(c)$ 的特征值属于 $\{0,\frac12,1\}$。

(2) 在引理 F.1(1) 中取 $a=b=c$，$d=z\in V(c,\nu)$，得 $2\nu[L(c),L(z)]+[L(z),L(c)]=0$，所以当
$\nu\in\{0,1\}$ 时 $L(z)$ 与 $L(c)$ 可交换，从而保持每个 $V(c,\mu)$。这给出第二和第三个包含关系，以及
$V(c,1)\circ V(c,0)\subseteq V(c,0)\cap V(c,1)=\{0\}$。对 $z,w\in V(c,\frac12)$，把取 $(a,b,d)=(z,w,c)$ 的引理 F.1(1)
作用于 $c$：由于 $L(w\circ c)=\frac12L(w)$ 且 $L(c\circ z)=\frac12L(z)$，前两项相消，于是对 $u=z\circ w$ 有
$[L(c),L(u)]c=L(c)^2u-L(c)u=0$。按 (1) 写 $u=u_1+u_{1/2}+u_0$，这就是 $-\frac14u_{1/2}=0$。

(3) $s^2=4c-4c+e=e$，$L(s)=2L(c)-I$ 在 $V(c,1),V(c,\frac12),V(c,0)$ 上分别为 $1,0,-1$，且 $P(s)=2L(s)^2-I$。写
$z=z'+z''$，$w=w'+w''$，其中 $z',w'\in V(c,0)+V(c,1)$，$z'',w''\in V(c,\frac12)$。由 (2)，$z'\circ w'+z''\circ w''$ 属于
$V(c,0)+V(c,1)$，而 $z'\circ w''+z''\circ w'$ 属于 $V(c,\frac12)$，所以
$P(s)(z\circ w)=(z'-z'')\circ(w'-w'')=P(s)z\circ P(s)w$。

(4) 各 $c_l$ 是 $x$ 的多项式，所以各 $L(c_l)$ 两两可交换（命题 F.2）；它们是对称的，因而有公共的
规范正交特征基。在一个公共特征向量上，由 (1)，$L(c_l)$ 有特征值 $\varepsilon_l\in\{0,\frac12,1\}$，并且
$\sum_l\varepsilon_l=1$，因为 $\sum_lL(c_l)=L(e)=I$。所以要么某一个 $\varepsilon_i=1$ 而其余为零，要么
对某一对 $i<j$ 有 $\varepsilon_i=\varepsilon_j=\frac12$ 而其余为零。所以 $V$ 是公共特征子空间
$W_{ii}$（其上 $L(c_i)=1$，对 $l\ne i$ 有 $L(c_l)=0$）与 $W_{ij}$（其上 $L(c_i)=L(c_j)=\frac12$，
其余 $L(c_l)=0$）的正交和，并且 $W_{ij}=V_{ij}$：满足 $c_i\circ z=z$（相应地 $c_i\circ z=c_j\circ z=\frac12z$）的向量在
其他公共特征子空间中没有分量，因为在那些子空间上 $L(c_i)$（相应地，数对 $(L(c_i),L(c_j))$）取不同的
值。由于
$a^2=\sum_la_l^2c_l$，在 $V_{ij}$（$i<j$）上 $L(a)=\frac12(a_i+a_j)$，$L(a^2)=\frac12(a_i^2+a_j^2)$，且
$P(a)=\frac12(a_i+a_j)^2-\frac12(a_i^2+a_j^2)=a_ia_j$；在 $V_{ii}$ 上 $L(a)=a_i$，$P(a)=2a_i^2-a_i^2$。$\square$

## F.2 平方锥是对称锥

设 $\Omega$ 是欧几里得空间 $V$ 中的非空开**凸锥**（对 $x\in\Omega$，$t>0$ 有 $tx\in\Omega$）。它的
（开）**对偶锥**是 $\Omega^*=\{y:\langle x,y\rangle>0\text{ 对所有 }x\in\bar\Omega\setminus\{0\}\text{ 成立}\}$；若 $\Omega^*=\Omega$，则称 $\Omega$ 是
**自对偶的**。令 $G(\Omega)=\{g\in GL(V):g\Omega=\Omega\}$；若 $G(\Omega)$ 在 $\Omega$ 上可迁地作用，则称 $\Omega$ 是**齐性的**。
**对称锥**是自对偶的齐性锥。若 $\Omega_1\subseteq\Omega_2$，
则 $\Omega_2^*\subseteq\Omega_1^*$。

**定理 F.5（自对偶性）。** 设 $V$ 是欧几里得若尔当代数，$y=\sum_i\lambda_ic_i$（定理 F.3）。

1. 以下各条等价：(i) $y\in Q$；(ii) 所有 $\lambda_i\ge0$；(iii) $L(y)\ge0$；(iv) 对所有 $x\in Q$ 有
   $\langle y,x\rangle\ge0$。特别地，$Q$ 是闭凸锥。
2. 以下各条等价：(i) $y\in\Omega:=\operatorname{int}Q$；(ii) 所有 $\lambda_i>0$；(iii) $L(y)>0$；(iv)
   对所有 $x\in Q\setminus\{0\}$ 有 $\langle y,x\rangle>0$。此外 $e\in\Omega$，$\bar\Omega=Q$，且 $\Omega^*=\Omega$：
   $\Omega$ 是自对偶的开凸锥。
3. 若 $x,y\in Q$ 且 $\langle x,y\rangle=0$，则 $x\circ y=0$。
4. $x^2=0$ 蕴含 $x=0$，并且 $Q\cap(-Q)=\{0\}$。

*证明。* 我们使用 $\langle y,x^2\rangle=\langle y\circ x,x\rangle=\langle L(y)x,x\rangle$ 与
$\langle x^2,e\rangle=\langle x,x\circ e\rangle=\|x\|^2$；后者连同 (1)(ii) 给出 (4)，因为由唯一性，
$-y$ 的特征值是各 $-\lambda_i$。

(1) (i)$\Rightarrow$(iii)：若 $y=x^2$，$x=\sum_i\mu_id_i$，则 $y=\sum_i\mu_i^2d_i$，并由引理 F.4(1)，
$L(y)=\sum_i\mu_i^2L(d_i)\ge0$。(iii)$\Leftrightarrow$(iv) 由第一个恒等式得出。
(iii)$\Rightarrow$(ii)：$0\le\langle L(y)c_i,c_i\rangle=\lambda_i\|c_i\|^2$。(ii)$\Rightarrow$(i)：
$y=(\sum_i\sqrt{\lambda_i}c_i)^2$。由 (iv)，$Q$ 是一族过 $0$ 的闭半空间的交。

(2) (ii)$\Rightarrow$(iii)：$L(y)=\sum_i\lambda_iL(c_i)\ge(\min_i\lambda_i)\sum_iL(c_i)=(\min_i\lambda_i)I$。
(iii)$\Rightarrow$(iv)：非零的 $x\in Q$ 是 $z^2$，其中 $z\ne0$，并且 $\langle y,x\rangle=\langle L(y)z,z\rangle>0$。
(iv)$\Rightarrow$(ii)：$\langle y,c_i\rangle=\lambda_i\|c_i\|^2$，且 $c_i=c_i^2\in Q\setminus\{0\}$。(iii)$\Rightarrow$(i)：若
$L(y)\ge mI$，$m>0$，则 $L(y')\ge(m-\|L\|\,\|y-y'\|)I$，其中 $\|L\|$ 是线性映射
$b\mapsto L(b)$ 的范数，所以 $y$ 的某个邻域由满足
$L>0$ 的元素组成，而由 (1) 它们属于 $Q$。(i)$\Rightarrow$(iii)：由 (1) 有 $L(y)\ge0$；若 $\langle L(y)z,z\rangle=0$ 且 $z\ne0$，则
$x=z^2\ne0$ 且 $\langle y,x\rangle=0$；对充分小的 $\varepsilon>0$，$y-\varepsilon x\in Q$，但
$\langle y-\varepsilon x,x\rangle=-\varepsilon\|x\|^2<0$，与 (1)(iv) 矛盾。其次，$L(e)=I>0$，并且对 $y\in Q$，
$y+\varepsilon e\in\Omega$ 趋于 $y$，所以 $\bar\Omega=Q$；于是由 (iv)，$\Omega^*=\{y:\langle y,x\rangle>0\ \forall x\in Q\setminus\{0\}\}=\Omega$。
$\Omega=\{L>0\}$ 是凸的。

(3) 设 $x=\sum_i\lambda_ic_i$，$\lambda_i\ge0$。则 $0=\sum_i\lambda_i\langle c_i,y\rangle$，由 (1) 其中每一项都 $\ge0$。
若 $\lambda_i>0$，则 $0=\langle y,c_i\rangle=\langle L(y)c_i,c_i\rangle$，而 $L(y)\ge0$，所以 $L(y)c_i=0$。因此
$x\circ y=\sum_i\lambda_i\,c_i\circ y=0$。$\square$

**定理 F.6（齐性与二次表示）。** 设 $V$ 是欧几里得若尔当代数，$Q$ 与
$\Omega$ 如定理 F.5 所示。

1. 对每个 $b\in V$：$e^{L(b)}Q=Q$，$e^{L(b)}\Omega=\Omega$，且 $e^{L(b)}e=\exp b$。
2. 每个 $y\in\Omega$ 都等于 $\exp(\log y)$。因此 $\Omega$ 是齐性的，从而是对称锥。
3. 对 $y\in\Omega$，$P(y)=e^{2L(\log y)}$；所以 $P(y)Q=Q$，$P(y)\Omega=\Omega$，且 $P(y)e=y^2$。
4. 对每个 $y\in V$ 有 $P(y)Q\subseteq Q$。

*证明。* (1) 设 $0\le s<1/\|L(b)\|$（若 $L(b)=0$ 则 $s\ge0$ 任意）。对称算子 $I-sL(b)$ 的特征值
$\ge1-s\|L(b)\|>0$。设 $x\in Q$，$z=(I-sL(b))^{-1}x=\sum_j\nu_jd_j$（定理 F.3）。利用 $z\circ d_j=\nu_jd_j$、
$\langle d_l,d_j\rangle=\langle d_l,d_j\circ d_j\rangle=\langle d_l\circ d_j,d_j\rangle=\delta_{lj}\|d_j\|^2$ 以及
$\langle b,d_j\rangle=\langle b,d_j^2\rangle=\langle L(b)d_j,d_j\rangle\le\|L(b)\|\,\|d_j\|^2$，得
$$
0\le\langle x,d_j\rangle=\langle z-s\,b\circ z,d_j\rangle=\nu_j\|d_j\|^2-s\langle b,z\circ d_j\rangle
=\nu_j\big(\|d_j\|^2-s\langle b,d_j\rangle\big),
$$
其中左边的不等式就是定理 F.5(1)(iv)。括号 $\ge(1-s\|L(b)\|)\|d_j\|^2>0$，所以所有 $\nu_j\ge0$，且
$z\in Q$。因此，对 $t\ge0$ 与 $n>t\|L(b)\|$，有 $(I-\frac tnL(b))^{-n}Q\subseteq Q$。当 $n\to\infty$ 时该算子趋于
$e^{tL(b)}$（在 $L(b)$ 的规范正交特征基下：$(1-t\beta/n)^{-n}\to e^{t\beta}$），而 $Q$ 是闭的，所以
$e^{tL(b)}Q\subseteq Q$。以 $-b$ 代替 $b$，得
$e^{-L(b)}Q\subseteq Q$，所以 $e^{L(b)}Q=Q$，并且同胚 $e^{L(b)}$ 把 $\operatorname{int}Q$ 映到自身上。
最后 $L(b)^re=b^r$（若尔当幂），所以由定理 F.3，对 $b=\sum_i\beta_id_i$ 有
$e^{L(b)}e=\sum_rb^r/r!=\sum_ie^{\beta_i}d_i=\exp b$。

(2) 若 $y=\sum_i\lambda_ic_i$，$\lambda_i>0$，则 $\log y=\sum_i(\log\lambda_i)c_i$，且 $\exp\log y=y$。所以
$y=e^{L(\log y)}e$，其中由 (1) 有 $e^{L(\log y)}\in G(\Omega)$；自对偶性就是定理 F.5(2)。

(3) 对 $y$ 的各 $c_i$ 应用引理 F.4(4)，$P(y)$ 与 $e^{2L(\log y)}$ 在 $V_{ij}$ 上都作用为
$\lambda_i\lambda_j=e^{\log\lambda_i+\log\lambda_j}$。再利用 (1) 以及 $P(y)e=2y\circ y-y^2=y^2$。

(4) 设 $y=\sum_i\lambda_ic_i$，$c=\sum_{\lambda_i\ge0}c_i$（一个幂等元），$s=2c-e=\sum_i\varepsilon_ic_i$，其中
$\varepsilon_i=\pm1$，并设 $|y|=\sum_i|\lambda_i|c_i$。由引理 F.4(4)，$P(s)$、$P(|y|)$ 与 $P(y)$ 在 $V_{ij}$ 上分别作用为
$\varepsilon_i\varepsilon_j$、$|\lambda_i\lambda_j|$ 与 $\lambda_i\lambda_j$，所以 $P(y)=P(s)P(|y|)$。对 $\eta>0$，
$|y|+\eta e\in\Omega$，所以由 (3) 有 $P(|y|+\eta e)Q\subseteq Q$；$P(\cdot)$ 是多项式且 $Q$ 是闭的，所以
$P(|y|)Q\subseteq Q$（$\eta\to0$）。$P(s)$ 是若尔当自同构（引理 F.4(3)），所以 $P(s)x^2=(P(s)x)^2\in Q$。$\square$

13.6 节在 $V_m$ 中由三条性质证明了 $\operatorname{tr}|P(y)x|\le\operatorname{tr}P(y)|x|$。对线性型
$\tau(x)=\langle x,e\rangle$，这些性质在每个欧几里得若尔当代数中都成立；它在
$M_{sa}$ 上就是 $\operatorname{Tr}$，在 $V_m$ 上就是 $\operatorname{tr}$（取上面的内积）：

**推论 F.7（每个欧几里得若尔当代数中的命题 13.14）。** 对欧几里得若尔当代数中的 $x,y$，
$\tau(|P(y)x|)\le\tau(P(y)|x|)$，且当 $x\in Q$ 或 $-x\in Q$ 时等号成立。

*证明。* $x_\pm\in Q$，$x=x_+-x_-$，$|x|=x_++x_-$（定理 F.3）。$A=P(y)x_+$ 与 $B=P(y)x_-$ 属于 $Q$（定理 F.6(4)），
且 $z=P(y)x=A-B$。设 $z=\sum_i\mu_ic_i$，$c=\sum_{\mu_i>0}c_i$；则 $|z|=z\circ(2c-e)$，并且
$\tau(a\circ b)=\langle a,b\rangle$。由于 $\langle Q,Q\rangle\ge0$（定理 F.5(1)）且 $c,e-c\in Q$，
$$
\tau(|z|)=\langle A,c\rangle-\langle A,e-c\rangle-\langle B,c\rangle+\langle B,e-c\rangle\le\tau(A)+\tau(B)=\tau(P(y)|x|).
$$
若 $x\in Q$，则 $x_-=0$ 且 $z=A\in Q$，所以 $|z|=z=P(y)|x|$；$-x\in Q$ 时同理。$\square$

代码对 $\operatorname{Herm}_4(\mathbb C)$ 以及洛伦兹锥 $Q\subseteq V_4$ 检验定理 F.5（负的配对值就是
(iv)$\Rightarrow$(ii) 中的见证 $c_i$），在 $V_4$ 中检验定理 F.6（在坐标
$(a,u)$ 下，算子 $L(b)$ 与 $P(y)$ 是 $5\times5$ 对称矩阵），并在
$\operatorname{Herm}_5(\mathbb R)$ 中检验推论 F.7：

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

## F.3 Koecher–Vinberg 定理

F.2 节证明了欧几里得若尔当代数的平方锥的内部是对称锥。其逆命题
从锥重建若尔当积。直到定理 F.12 为止，$\Omega$ 是维数为 $N$ 的欧几里得空间 $V$ 中的自对偶开凸锥，
$dy$ 是规范正交坐标下的 Lebesgue 测度（正交坐标变换满足 $|\det|=1$）。

**引理 F.8.** 设 $\Omega$ 是自对偶的。

1. $0\notin\Omega$；对 $x,y\in\bar\Omega$ 有 $\langle x,y\rangle\ge0$；对每个 $x\in\Omega$，存在 $\delta>0$，使得对所有
   $y\in\bar\Omega$ 有 $\langle x,y\rangle\ge\delta\|y\|$。
2. 若 $x\in\bar\Omega\setminus\Omega$，则存在 $u\in\bar\Omega\setminus\{0\}$ 使得 $\langle x,u\rangle=0$。若 $y\in\Omega$，
   $u\in\bar\Omega$ 且 $t\ge0$，则 $y+tu\in\Omega$。
3. $\Omega$ 是 $\bar\Omega$ 的内部。
4. $G(\Omega)$ 是 $GL(V)$ 的闭子群，并且 $g\in G(\Omega)$ 蕴含 $g^T\in G(\Omega)$。

*证明。* (1) 若 $0\in\Omega$，则 $\Omega$ 包含一个以 $0$ 为中心的球，从而 $\Omega=V$；但 $V^*=\emptyset$。对
$x,y\in\Omega$，由于 $y\in\Omega^*$，有 $\langle x,y\rangle>0$；由连续性，在 $\bar\Omega$ 上 $\langle x,y\rangle\ge0$。若
闭球 $\bar B(x,\delta)$ 含于 $\Omega$ 且 $y\in\bar\Omega\setminus\{0\}$，则
$\langle x-\delta y/\|y\|,y\rangle\ge0$。(2) $x\notin\Omega=\Omega^*$ 给出 $u\in\bar\Omega\setminus\{0\}$，使得
$\langle u,x\rangle\le0$，而由 (1) 有 $\langle u,x\rangle\ge0$。对 $x'\in\bar\Omega\setminus\{0\}$，
$\langle x',y+tu\rangle\ge\langle x',y\rangle>0$，所以 $y+tu\in\Omega^*=\Omega$。(3) 若 $\bar B(y,\varepsilon)\subseteq\bar\Omega$，
则对 $x\in\bar\Omega\setminus\{0\}$ 有 $\langle x,y-\varepsilon x/\|x\|\rangle\ge0$，所以 $\langle x,y\rangle\ge\varepsilon\|x\|>0$，从而
$y\in\Omega^*=\Omega$。(4) 对 $g\in G(\Omega)$，$g$ 把 $\bar\Omega\setminus\{0\}$ 映到自身上，所以对 $y\in\Omega$ 与
$x\in\bar\Omega\setminus\{0\}$ 有 $\langle x,g^Ty\rangle=\langle gx,y\rangle>0$：$g^T\Omega\subseteq\Omega$。对 $g^{-1}$ 同样论证
得 $g^T\Omega=\Omega$。若 $g_k\in G(\Omega)$ 收敛到 $g\in GL(V)$，则 $g\Omega\subseteq\bar\Omega$ 是开集，所以由 (3) 有
$g\Omega\subseteq\Omega$；同样 $g^{-1}=\lim g_k^{-1}$ 把 $\Omega$ 映入 $\Omega$，所以 $g\Omega=\Omega$。$\square$

$\Omega$ 的**特征函数**是 $\varphi(x)=\int_{\Omega^*}e^{-\langle x,y\rangle}dy=\int_\Omega e^{-\langle x,y\rangle}dy$。
令 $f=\log\varphi$，设 $H(x)$ 是满足 $\langle H(x)u,v\rangle=D^2f(x)[u,v]$ 的对称算子，并令
$x^*=-\nabla f(x)$，即 $Df(x)[u]=-\langle x^*,u\rangle$。

**引理 F.9（特征函数）。** 设 $\Omega$ 是自对偶的。

1. $\varphi$ 在 $\Omega$ 上有限、为正且为 $C^\infty$，并且可以在积分号下求导。
2. 对 $t>0$ 有 $\varphi(tx)=t^{-N}\varphi(x)$，对 $g\in G(\Omega)$ 有 $\varphi(gx)=|\det g|^{-1}\varphi(x)$。
3. 对 $x\in\Omega$ 有 $H(x)>0$ 与 $x^*\in\Omega$。
4. $H(x)x=x^*$ 且 $D^3f(x)[x,u,v]=-2\langle H(x)u,v\rangle$；对 $g\in G(\Omega)$ 有 $(gx)^*=g^{-T}x^*$ 与
   $H(gx)=g^{-T}H(x)g^{-1}$。
5. 对某个常数 $c>0$ 有 $\varphi(x)\ge c\|x\|^{-N}$，并且若 $x_k\in\Omega$ 收敛到
   $\bar\Omega\setminus\Omega$ 中的一点，则 $\varphi(x_k)\to\infty$。

*证明。* (1) 设 $x\in\Omega$，$\delta$ 如引理 F.8(1) 所示。对 $\|x'-x\|<\delta/2$ 与 $y\in\Omega$ 有
$\langle x',y\rangle\ge\frac\delta2\|y\|$，所以被积函数关于 $x'$ 的 $k$ 阶导数
$(-1)^k\prod_{i}\langle u_i,y\rangle e^{-\langle x',y\rangle}$ 以 $\prod_i\|u_i\|\,\|y\|^ke^{-\delta\|y\|/2}$ 为界，后者
可积（$\|y\|\ge N^{-1/2}\sum_i|y_i|$ 以及 Fubini 定理）。控制收敛定理允许在积分号下
求导。(2) 作替换 $y=y'/t$，它把 $\Omega$ 映到 $\Omega$ 上；以及替换 $y'=g^Ty$，由引理
F.8(4) 它把 $\Omega$ 映到 $\Omega$ 上，且 $\langle gx,y\rangle=\langle x,y'\rangle$，$dy=|\det g|^{-1}dy'$。
(3) 设 $\mu_x$ 是 $\Omega$ 上的概率测度 $e^{-\langle x,y\rangle}dy/\varphi(x)$。由 (1)，
$Df(x)[u]=-\mathbb E\langle u,y\rangle$，而 $D^2f(x)[u,u]=\mathbb E\langle u,y\rangle^2-(\mathbb E\langle u,y\rangle)^2$ 是
$\langle u,y\rangle$ 的方差。它为零仅当 $\langle u,y\rangle$ $\mu_x$-几乎处处等于某个常数 $\gamma$，即在 $\Omega$ 上
Lebesgue-几乎处处如此（$\mu_x$ 的密度为正）；集合
$\{y\in\Omega:\langle u,y\rangle\ne\gamma\}$ 是开的零测集，因而是空集，所以线性函数 $\langle u,\cdot\rangle$ 在一个球上
为常数，这迫使 $u=0$。并且 $x^*=\mathbb E\,y$，所以对 $z\in\bar\Omega\setminus\{0\}$ 有
$\langle z,x^*\rangle=\mathbb E\langle z,y\rangle>0$：$x^*\in\Omega^*=\Omega$。
(4) 由 (2)，$f(tx)=f(x)-N\log t$，所以在 $\Omega$ 上 $Df(x)[x]=-N$。对 $x$ 沿 $u$ 求导：
$D^2f(x)[x,u]+Df(x)[u]=0$，即 $\langle H(x)x,u\rangle=\langle x^*,u\rangle$；再沿 $v$ 求导一次：
$D^3f(x)[x,u,v]+2D^2f(x)[u,v]=0$。对 $f(gx)=f(x)-\log|\det g|$ 求一阶和二阶导数，得
$\langle(gx)^*,gu\rangle=\langle x^*,u\rangle$ 与 $\langle H(gx)gu,gv\rangle=\langle H(x)u,v\rangle$。
(5) 对 $\|y\|\le1/\|x\|$ 有 $\langle x,y\rangle\le\|x\|\|y\|\le1$，所以 $\varphi(x)\ge e^{-1}\operatorname{vol}(\Omega\cap B(0,1/\|x\|))=c\|x\|^{-N}$，
其中 $c=e^{-1}\operatorname{vol}(\Omega\cap B(0,1))>0$。设 $x_k\to x_0\in\bar\Omega\setminus\Omega$。由 Fatou 引理，
$\liminf\varphi(x_k)\ge\int_\Omega e^{-\langle x_0,y\rangle}dy$。取 $u$ 如引理 F.8(2) 所示，取球 $B(y_1,r)\subseteq\Omega$
以及互不相交的球 $B_l=B(y_1+2lru/\|u\|,r)$，$l\ge0$；由引理 F.8(2) 它们含于 $\Omega$，并且在 $B_l$ 上
$\langle x_0,y\rangle=\langle x_0,y-2lru/\|u\|\rangle\le\|x_0\|(\|y_1\|+r)$。所以该积分是无穷大。$\square$

**引理 F.10（单参数群）。** 设 $G$ 是 $GL(V)$ 的闭子群，
$\mathfrak g=\{Z:e^{tZ}\in G\text{ 对所有 }t\in\mathbb R\text{ 成立}\}$。

1. 若 $c:[0,s_0)\to G$ 当 $s\to0^+$ 时满足 $c(s)=I+sZ+o(s)$，则 $Z\in\mathfrak g$。
2. $\mathfrak g$ 是李代数：一个在 $[Z,W]=ZW-WZ$ 下封闭的线性子空间。

*证明。* (1) 固定 $t>0$，令 $A_n=c(t/n)$，$B_n=e^{tZ/n}$。则 $\|A_n-B_n\|=o(1/n)$，并且对某个 $K$ 与充分大的 $n$ 有
$\|A_n\|,\|B_n\|\le1+K/n$，所以
$$
\|A_n^n-e^{tZ}\|=\Big\|\sum_{r=0}^{n-1}A_n^r(A_n-B_n)B_n^{n-1-r}\Big\|\le n\,e^K\|A_n-B_n\|\to0 .
$$
$A_n^n\in G$ 且 $G$ 在 $GL(V)$ 中是闭的，所以 $e^{tZ}\in G$；又 $e^{-tZ}=(e^{tZ})^{-1}\in G$。(2) 数乘是显然的。
对 $c(s)=e^{sZ}e^{sW}=I+s(Z+W)+O(s^2)$ 应用 (1)，并对 $c(s)=e^{aZ}e^{aW}e^{-aZ}e^{-aW}$（$a=\sqrt s$）应用 (1)：把每个
因子展开到二阶，$a$ 阶的项相消，而 $a^2$ 阶的项加起来是 $ZW-WZ$，所以
$c(s)=I+s[Z,W]+O(s^{3/2})$。$\square$

**命题 F.11（基点）。** 设 $\Omega$ 是对称锥，$G=G(\Omega)$，$\mathfrak g$ 如引理 F.10 所示。

1. 对 $x\in\Omega$ 有 $x^{**}=x$ 与 $H(x^*)=H(x)^{-1}$。
2. 存在 $e\in\Omega$ 使得 $e^*=e$，并且每个这样的 $e$ 都满足 $H(e)=I$。
3. 对 $x\in\Omega$ 以及任意满足 $ge=x$ 的 $g\in G$：$H(x)^{-1}=gg^T\in G$。
4. 若 $Z\in\mathfrak g$ 且 $Ze=0$，则 $Z^T=-Z$。

*证明。* (1) 映射 $x\mapsto x^*$ 是从 $\Omega$ 到 $\Omega$ 的光滑映射（引理 F.9(3)），其导数为 $-H(x)$。由于
$g^{-T}\in G$（引理 F.8(4)），引理 F.9(2,4) 给出
$f(gx)+f((gx)^*)=f(x)-\log|\det g|+f(g^{-T}x^*)=f(x)+f(x^*)$，所以 $F(x)=f(x)+f(x^*)$ 在轨道
$G x=\Omega$ 上为常数。因此 $0=DF(x)[u]=-\langle x^*,u\rangle+\langle x^{**},H(x)u\rangle$，即 $H(x)x^{**}=x^*=H(x)x$
（引理 F.9(4)），又因 $H(x)$ 可逆，得 $x^{**}=x$。对 $x^{**}=x$ 求导得 $H(x^*)H(x)=I$。
(2) $F_2(x)=f(x)+\frac12\|x\|^2$ 满足 $F_2(x)\ge\log c-N\log\|x\|+\frac12\|x\|^2$（引理 F.9(5)），当
$\|x\|\to\infty$ 时它趋于 $\infty$。所以 $S=\{x\in\Omega:F_2(x)\le F_2(x_1)\}$（$x_1\in\Omega$）是有界的；它在 $V$ 中是闭的，因为
$\bar\Omega\setminus\Omega$ 中的极限点会使 $f$ 在 $S$ 上无界（引理 F.9(5)）。连续函数 $F_2$ 在紧集 $S$ 上于开集 $\Omega$
中的某点 $e$ 处达到最小值，在该点 $0=\nabla F_2(e)=-e^*+e$。于是由 (1)，$H(e)=H(e^*)=H(e)^{-1}$，
所以 $H(e)^2=I$，又因 $H(e)>0$，得 $H(e)=I$。
(3) 由引理 F.9(4) 与 (2)，$H(x)=H(ge)=g^{-T}g^{-1}$，而由引理 F.8(4)，$gg^T\in G$。
(4) $e^{tZ}\in G$ 固定 $e$，所以 $I=H(e)=e^{-tZ^T}H(e)e^{-tZ}=e^{-tZ^T}e^{-tZ}$；在 $t=0$ 处求导。$\square$

**定理 F.12（Koecher–Vinberg）。** 设 $\Omega$ 是欧几里得空间 $V$ 中的非空开凸锥。则 $\Omega$ 是
对称锥，当且仅当 $V$ 上存在乘积 $\circ$，使得 $(V,\circ,\langle\cdot,\cdot\rangle)$ 是
欧几里得若尔当代数，且 $\Omega$ 是其平方锥的内部。对于对称锥，可以取任意满足 $e^*=e$ 的
$e\in\Omega$ 作为单位元，并取乘积
$$
\langle u\circ v,w\rangle=-\tfrac12\,D^3\log\varphi(e)[u,v,w].
$$

*证明。* "当"的方向就是定理 F.5(2) 与 F.6(2)。反之，设 $\Omega$ 是对称锥，$e$ 如命题 F.11(2) 所示，并
由 $\langle L(u)v,w\rangle=-\frac12D^3f(e)[u,v,w]$ 定义 $L(u)$，令 $u\circ v=L(u)v$。

*第 1 步：交换性、对称性、单位元。* $D^3f(e)$ 关于它的三个变元对称，所以 $u\circ v=v\circ u$，且
$L(u)^T=L(u)$。由引理 F.9(4) 与 $H(e)=I$，$\langle e\circ v,w\rangle=\langle H(e)v,w\rangle=\langle v,w\rangle$：$e$ 是单位元。

*第 2 步：$L(u)\in\mathfrak g$。* $R(x)=H(x)^{-1}$ 属于 $G$（命题 F.11(3)），是光滑的，且 $R(e)=I$。它的
导数是 $DR(e)[u]=-H(e)^{-1}DH(e)[u]H(e)^{-1}=-DH(e)[u]=2L(u)$，因为 $\langle DH(e)[u]v,w\rangle=D^3f(e)[u,v,w]$。
对 $c(s)=R(e+su)$（$s\ge0$ 充分小，使得 $e+su\in\Omega$）以及闭子群 $G$（引理 F.8(4)）应用引理 F.10(1)，
得 $2L(u)\in\mathfrak g$，从而 $L(u)\in\mathfrak g$。

*第 3 步：唯一性。* 若 $Y\in\mathfrak g$ 是对称的，则 $Y=L(Ye)$：对称的 $Y-L(Ye)\in\mathfrak g$ 把 $e$ 映为零
（$L(w)e=w\circ e=w$），所以由命题 F.11(4) 它也是反对称的，因而为 $0$。

*第 4 步：内导子。* 令 $D=[L(x),L(y)]$。它属于 $\mathfrak g$（引理 F.10(2)），$D^T=-D$，且
$De=x\circ y-y\circ x=0$。对 $z\in V$，$[D,L(z)]\in\mathfrak g$ 是对称的，且 $[D,L(z)]e=Dz-L(z)De=Dz$，所以
由第 3 步，$[D,L(z)]=L(Dz)$；作用于 $w$：$D(z\circ w)=Dz\circ w+z\circ Dw$。

*第 5 步：若尔当恒等式。* 对 $D=[L(y),L(x)]$，有 $Dx=y\circ x^2-x\circ(x\circ y)$ 与
$D(x^2)=y\circ(x\circ x^2)-x\circ(y\circ x^2)$。第 4 步给出 $D(x^2)=2x\circ Dx$，即对所有 $y$ 有
$y\circ(x\circ x^2)=3x\circ(x^2\circ y)-2x\circ(x\circ(x\circ y))$：
$$
L(x\circ x^2)=3L(x)L(x^2)-2L(x)^3 .
$$
所有 $L$ 都是对称的，所以取转置得 $L(x\circ x^2)=3L(x^2)L(x)-2L(x)^3$，从而 $[L(x),L(x^2)]=0$。所以 $V$ 是
欧几里得若尔当代数；设 $\Omega_J$ 是其平方锥的内部。

*第 6 步：$\Omega=\Omega_J$。* 对 $b\in V$，由第 2 步有 $e^{L(b)}\in G$，所以 $\exp b=e^{L(b)}e\in\Omega$（定理 F.6(1)）。由
定理 F.6(2)，$\Omega_J\subseteq\Omega$。两个锥关于 $\langle\cdot,\cdot\rangle$ 都是自对偶的（定理 F.5(2)），所以
$\Omega=\Omega^*\subseteq\Omega_J^*=\Omega_J$。$\square$

这一证明沿用了 J. Faraut, A. Korányi, *Analysis on Symmetric Cones* (Oxford 1994) 第 I–III 章的结构；
在他们使用 $G(\Omega)$ 的李群结构之处，引理 F.10 与命题 F.11(3) 直接产生所需的
单参数群。该定理由 M. Koecher (1957) 与 E. B. Vinberg (1960) 证明。

**例 F.13（半定规划的锥）。** 设 $V=\operatorname{Herm}_n(\mathbb R)$，$\langle x,y\rangle=\operatorname{Tr}(xy)$，
$N=n(n+1)/2$，$\Omega$ 为正定矩阵全体，由定理 F.5 与 F.6 它是对称锥。则
$\varphi(x)=\varphi(1)\det(x)^{-(n+1)/2}$，$e=\kappa1$，其中 $\kappa=\sqrt{(n+1)/2}$，而定理 F.12 的乘积是
$u\circ_\kappa v=\frac1{2\kappa}(uv+vu)$，所以 $x\mapsto x/\kappa$ 是到通常的若尔当代数上的同构。

*证明。* 对称矩阵按其特征投影的谱分解 $x=\sum_i\lambda_ip_i$ 就是
定理 F.3 的分解（唯一性；由于各 $p_i$ 可交换，$p_i\circ p_j=p_ip_j$），所以由定理 F.5(2)，平方锥的
内部确实是正定矩阵的集合。对 $a\in\Omega$，$P(a)x=axa$
（命题 13.2(1)）属于 $G(\Omega)$（定理 F.6(3)）。在 $a$ 的规范正交特征基下，
矩阵 $E_{ii}$ 与 $(E_{ij}+E_{ji})/\sqrt2$（$i<j$）构成 $P(a)$ 的特征向量的规范正交基，对应的
特征值为 $a_ia_j$，所以 $\det P(a)=\prod_{i\le j}a_ia_j=\det(a)^{n+1}$（每个 $a_i$ 出现 $n+1$ 次）。引理 F.9(2) 给出
$\varphi(a^2)=\varphi(P(a)1)=\det(a)^{-(n+1)}\varphi(1)$，并且每个 $x\in\Omega$ 都是 $a^2$，其中 $a=\sqrt x$。利用
$D\log\det(x)[u]=\operatorname{Tr}(x^{-1}u)$（由 $\det(x+su)=\det x\,\det(1+sx^{-1}u)=\det x\,(1+s\operatorname{Tr}(x^{-1}u)+O(s^2))$ 得出）
与 $D(x^{-1})[u]=-x^{-1}ux^{-1}$：
$x^*=\kappa^2x^{-1}$；所以 $e^*=e$ 意味着 $e^2=\kappa^2 1$ 且 $e$ 正定，即 $e=\kappa1$；$D^2f(x)[v,w]=\kappa^2\operatorname{Tr}(x^{-1}vx^{-1}w)$，并且
$D^3f(e)[u,v,w]=-\kappa^2\kappa^{-3}\operatorname{Tr}(uvw+vuw)=-\frac2\kappa\operatorname{Tr}((u\circ v)w)$。$\square$

代码在 $\operatorname{Herm}_3(\mathbb R)$ 的规范正交坐标下，通过对 $f(x)=-\frac{n+1}2\log\det x$ 自动微分计算 $e$ 处的
$D^3f$，并检验命题 F.11 与定理 F.12 的结论：

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

## F.4 序同构：有限维中的 Kadison 定理

设 $\mathcal M=\bigoplus_cM_{k_c}(\mathbb C)$，$\mathcal N=\bigoplus_dM_{k'_d}(\mathbb C)$，视为分块对角矩阵；
$x\ge0$ 表示 $x=x^*$ 且特征值非负（等价地，每个块都是半正定的），并且
$[0,1]=\{a\in\mathcal M_{sa}:0\le a\le1\}$。**投影**是指 $p=p^*=p^2$。线性映射 $\Phi$ 称为**正的**，如果
$x\ge0$ 蕴含 $\Phi(x)\ge0$；**若尔当 $*$-同构**是满足 $\Phi(x^*)=\Phi(x)^*$ 与
$\Phi(x\circ y)=\Phi(x)\circ\Phi(y)$ 的线性双射（13.5 节）。凸集 $K$ 中的点 $a$ 称为**端点**，如果它不是
两个不同点 $b,b'\in K$ 的中点 $\frac12(b+b')$（等价地，它不能写成 $tb+(1-t)b'$，其中 $b\ne b'$
且 $0<t<1$：这样的 $a$ 是 $a\pm\varepsilon(b-b')\in K$ 的中点，其中 $\varepsilon>0$ 充分小）。

**引理 F.14.** 在 $\mathcal M$ 中：

1. 凸集 $[0,1]$ 的端点恰好是投影。
2. 对投影 $p,q$：$pq=0$ 当且仅当 $p+q\le1$，当且仅当 $p\circ q=0$。
3. 对投影 $p,q$：$p\le q$ 当且仅当 $pq=p$，当且仅当 $p\circ q=p$。极小非零投影是（某个块中的）秩一
   投影。
4. 投影 $p$ 是中心的，当且仅当对所有 $x\in\mathcal M$ 有 $pxp+(1-p)x(1-p)=x$。中心投影是
   块单位 $P_c$ 的和，而极小的非零中心投影就是各 $P_c$。

*证明。* (1) 设 $p$ 是投影，$p=\frac12(a+b)$，$a,b\in[0,1]$。对 $\xi\in\ker p$，
$0=\langle\xi,a\xi\rangle+\langle\xi,b\xi\rangle$，两项都 $\ge0$，所以 $\|a^{1/2}\xi\|^2=0$，从而 $a\xi=0$；对 $p$ 的值域中的
$\xi$，对 $1-a,1-b$ 作同样论证得 $a\xi=\xi$。所以 $a=p$。反之，若 $a\in[0,1]$ 不是投影，
则它有特征值 $\lambda\in(0,1)$；取谱投影 $q_\lambda\in\mathcal M$（$a$ 的多项式）以及
$\varepsilon=\min(\lambda,1-\lambda)$，则 $a\pm\varepsilon q_\lambda\in[0,1]$，而 $a$ 是它们的中点。
(2) 若 $pq=0$，则 $qp=(pq)^*=0$，$p+q$ 是投影，且 $1-(p+q)=(1-p-q)^2\ge0$；同时 $p\circ q=0$。若
$p+q\le1$，则对每个 $\xi$ 有 $\|qp\xi\|^2=\langle p\xi,qp\xi\rangle\le\langle p\xi,(1-p)p\xi\rangle=0$，所以 $qp=0=pq$。若
$p\circ q=0$，则 $pq=-qp$，所以 $pq=p(pq)=-pqp=(qp)p=qp$，从而 $pq=-pq=0$。
(3) $pq=p$ 蕴含 $qp=p$、$p\circ q=p$ 以及 $q-p=(q-p)^2\ge0$。若 $p\circ q=p$，用 $p$ 从左边和从右边乘
$pq+qp=2p$，得 $pq+pqp=2p=pqp+qp$，所以 $pq=qp=p$。若 $p\le q$，则 $p+(1-q)\le1$，由 (2) 得
$p(1-q)=0$。秩 $\ge2$ 的投影严格大于到其值域中某条直线上的投影 $\xi\xi^*$，这条直线位于
一个块内（$\xi$ 是 $p$ 的某个块的值域中的单位向量）。秩一投影 $p$ 位于一个块内，并且若
$0\ne q\le p$，则 $q=pq$ 的值域含于 $p$ 的值域中，所以 $q=p$。
(4) $x=pxp+px(1-p)+(1-p)xp+(1-p)x(1-p)$，而 $px(1-p)+(1-p)xp=0$ 当且仅当两项都为零（用 $p$ 与 $1-p$ 去乘），
当且仅当 $px=pxp=xp$。中心元与每个块的所有矩阵单位可交换，所以它在块 $c$ 上是数乘 $z_c$，
而对投影有 $z_c\in\{0,1\}$。$\square$

**定理 F.15（Kadison，有限维）。** 设 $\Phi:\mathcal M\to\mathcal N$ 是线性双射，
$\Phi(1)=1$，且 $\Phi$ 与 $\Phi^{-1}$ 都是正的。则 $\Phi$ 是若尔当 $*$-同构。反之，每个
若尔当 $*$-同构都是保单位的，并且它及其逆都是正的。

*证明。* 每个 $x\in\mathcal M_{sa}$ 都是 $x_+-x_-$，其中 $x_\pm\ge0$，所以 $\Phi(x)$ 是自伴的；写成 $x=h+ik$，
$h,k$ 自伴，则 $\Phi(x^*)=\Phi(h)-i\Phi(k)=\Phi(x)^*$。对 $\Phi^{-1}$ 同样成立，所以 $\Phi$ 是实线性
双射 $\mathcal M_{sa}\to\mathcal N_{sa}$，且 $x\ge0\iff\Phi(x)\ge0$。由于 $\Phi(1)=1$，$\Phi$ 把 $[0,1]$ 映到 $[0,1]$ 上。
凸集之间的线性双射把端点映到端点上（若 $\Phi(a)=\frac12(b+b')$，则
$a=\frac12(\Phi^{-1}b+\Phi^{-1}b')$），所以由引理 F.14(1)，$\Phi$ 把投影映到投影上，并由引理
F.14(2)，$pq=0\iff p+q\le1\iff\Phi(p)+\Phi(q)\le1\iff\Phi(p)\Phi(q)=0$。对 $x=\sum_i\lambda_ip_i\in\mathcal M_{sa}$（谱
分解，$i\ne j$ 时 $p_ip_j=0$），$\Phi(x)^2=\sum_i\lambda_i^2\Phi(p_i)=\Phi(x^2)$。极化公式
$x\circ y=\frac12((x+y)^2-x^2-y^2)$ 给出 $\mathcal M_{sa}$ 上的 $\Phi(x\circ y)=\Phi(x)\circ\Phi(y)$，再由复
双线性推广到 $\mathcal M=\mathcal M_{sa}+i\mathcal M_{sa}$ 上。

反之，对所有 $y$ 有 $\Phi(1)\circ\Phi(y)=\Phi(y)$，且 $\Phi$ 是满射，所以 $\Phi(1)$ 是 $\circ$ 的单位元，从而
$\Phi(1)=\Phi(1)\circ1=1$。每个 $x\ge0$ 都是 $h^2$，其中 $h=\sqrt x$ 自伴（定理 1.6），并且由于 $\Phi(h)$ 自伴，
$\Phi(h^2)=\Phi(h)^2\ge0$；$\Phi^{-1}$ 也是若尔当 $*$-同构。$\square$

Kadison 对所有有单位的 C\*-代数证明了这一定理（R. V. Kadison, *A generalized Schwarz inequality and algebraic
invariants for operator algebras*, Ann. of Math. 56 (1952) 494–503）；本课程不使用那个版本。

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

## F.5 $\bigoplus_cM_{k_c}$ 的若尔当 $*$-同构

**引理 F.16.** 设 $\Phi:\mathcal M\to\mathcal N$ 是若尔当 $*$-同构。则 $\Phi(1)=1$，
$\Phi(xyx)=\Phi(x)\Phi(y)\Phi(x)$，并且 $\Phi$ 把 $\mathcal M$ 的投影双射地映到 $\mathcal N$ 的投影上，
保持并反映 $pq=0$、$p\le q$、极小性、中心性以及在中心投影中的极小性。

*证明。* 与定理 F.15 中一样，$\Phi(1)=1$。在每个结合代数中 $xyx=2x\circ(x\circ y)-(x\circ x)\circ y$，这
给出第二个结论。$p=p^*=p\circ p$ 被保持，引理 F.14(2,3,4) 的刻画也被保持，因为它们
只用到 $\circ$、三重积 $xyx$、$1$ 以及对整个代数的量词；$\Phi^{-1}$ 也是若尔当
$*$-同构。$\square$

**定理 F.17（$M_n$ 的若尔当 $*$-自同构）。** 每个若尔当 $*$-同构 $\Phi:M_n(\mathbb C)\to M_n(\mathbb C)$
都是 $x\mapsto UxU^*$ 或 $x\mapsto Ux^TU^*$，其中 $U$ 为酉矩阵。对 $n\ge2$，第一种映射是乘性的，第二种是
反乘性的，所以没有映射同时具有两种形式。

*证明。* 对 $n=1$，$\Phi=\mathrm{id}$。设 $n\ge2$，$E_{ij}$ 为矩阵单位。由引理 F.16，各 $\Phi(E_{ii})$ 是
两两正交的秩一投影，$\Phi(E_{ii})=\xi_i\xi_i^*$，其中 $\xi_i$ 规范正交；借助满足
$We_i=\xi_i$ 的酉矩阵 $W$，若尔当 $*$-同构 $\Phi_1(x)=W^*\Phi(x)W$ 满足 $\Phi_1(E_{ii})=E_{ii}$。

设 $i\ne j$。空间 $S_{ij}=\operatorname{span}\{E_{ij},E_{ji}\}$ 等于 $\{x:pxp=x,\ E_{ii}xE_{ii}=E_{jj}xE_{jj}=0\}$，其中
$p=E_{ii}+E_{jj}$，所以由引理 F.16，$\Phi_1(S_{ij})\subseteq S_{ij}$。写 $\Phi_1(E_{ij})=\alpha E_{ij}+\beta E_{ji}$。由
$E_{ij}\circ E_{ij}=0$：$0=\Phi_1(E_{ij})^2=\alpha\beta(E_{ii}+E_{jj})$，所以 $\alpha,\beta$ 中恰有一个非零
（$\Phi_1$ 是单射）。若 $\Phi_1(E_{ij})\in\mathbb CE_{ij}$，令 $\varepsilon(i,j)=+$；若 $\Phi_1(E_{ij})\in\mathbb CE_{ji}$，令其为 $-$；
由于 $\Phi_1(E_{ji})=\Phi_1(E_{ij})^*$，有 $\varepsilon(j,i)=\varepsilon(i,j)$。对互不相同的 $i,j,k$，
$E_{ij}\circ E_{jk}=\frac12E_{ik}\ne0$；若 $\varepsilon(i,j)=+$ 且 $\varepsilon(j,k)=-$，则
$\Phi_1(E_{ij})\circ\Phi_1(E_{jk})\in\mathbb C\,E_{ij}\circ E_{kj}=\{0\}$，矛盾；对 $(-,+)$ 同理。所以
只要 $i,j,k$ 互不相同就有 $\varepsilon(i,j)=\varepsilon(j,k)$，从而 $\varepsilon$ 是常数：两个指标对要么有一个公共
指标，要么通过第三个指标对相连。若 $\varepsilon=-$，则把 $\Phi_1$ 换成 $\Phi_1\circ T$（$T$ 为转置，由命题 13.9
它是若尔当 $*$-自同构），它满足 $\varepsilon=+$，并且仍固定各 $E_{ii}$。

现在对所有 $i,j$ 有 $\Phi_1(E_{ij})=\alpha_{ij}E_{ij}$，其中 $\alpha_{ii}=1$，$\alpha_{ji}=\bar\alpha_{ij}$。由
$E_{ij}\circ E_{ji}=\frac12(E_{ii}+E_{jj})$：$|\alpha_{ij}|^2=1$；由 $E_{ij}\circ E_{jk}=\frac12E_{ik}$（$i,j,k$ 互不相同）：
$\alpha_{ij}\alpha_{jk}=\alpha_{ik}$。令 $u_i=\alpha_{i1}$，得 $\alpha_{ij}=\alpha_{i1}\alpha_{1j}=u_i\bar u_j$（对
$i=j$、$i=1$ 或 $j=1$ 也成立）。所以 $\Phi_1(x)=DxD^*$，其中 $D=\operatorname{diag}(u_i)$，并且 $\Phi(x)=WDxD^*W^*$ 或
$\Phi(x)=WDx^TD^*W^*$。若两种形式一致，$T$ 就会是乘性的，但 $(E_{12}E_{21})^T=E_{11}\ne E_{22}=E_{12}^TE_{21}^T$。
$\square$

这在有限维情形下、对所有 $n$ 证明了 13.5 节所引用的 N. Jacobson、C. E. Rickart（Trans. Amer. Math. Soc. 69 (1950) 479–502）
与 I. N. Herstein（Trans. Amer. Math. Soc. 81 (1956) 331–341）的结果；命题 13.11 是 $n=2$ 的情形。
$\operatorname{Herm}_n(\mathbb C)$ 的若尔当自同构可以由 $\Phi(x+iy)=\Phi(x)+i\Phi(y)$ 延拓为 $M_n$ 的
若尔当 $*$-自同构，所以它也具有这两种形式之一。

**定理 F.18（直和）。** 设 $\Phi:\bigoplus_{c=1}^CM_{k_c}\to\bigoplus_{d=1}^{C'}M_{k'_d}$ 是若尔当
$*$-同构。则 $C'=C$，并且存在满足 $k'_{\sigma(c)}=k_c$ 的双射 $\sigma$、酉矩阵 $U_c\in M_{k_c}$ 以及
一个块的集合 $S$，使得
$$
\Phi(x)_{\sigma(c)}=U_cx_cU_c^*\ (c\notin S),\qquad \Phi(x)_{\sigma(c)}=U_cx_c^TU_c^*\ (c\in S).
$$
反之，每个这样的映射都是若尔当 $*$-同构，并且 $S$ 在相差 $k_c=1$ 的块的意义下唯一。

*证明。* 由引理 F.16 与引理 F.14(4)，对极小中心投影之间的某个双射 $\sigma$ 有
$\Phi(P_c)=Q_{\sigma(c)}$，所以 $C=C'$。若 $x=P_cxP_c$，则 $\Phi(x)=Q_{\sigma(c)}\Phi(x)Q_{\sigma(c)}$，所以 $\Phi$ 把块 $c$ 映入块
$\sigma(c)$，而 $\Phi^{-1}$ 把块 $\sigma(c)$ 映入块 $c$；因此 $k_c^2=k'^2_{\sigma(c)}$，并且其限制是
若尔当 $*$-同构 $M_{k_c}\to M_{k_c}$，由定理 F.17 它具有两种形式之一。由于 $x=\sum_cP_cxP_c$，这就
决定了 $\Phi$。逆命题逐块成立（命题 13.9），而 $S$ 的唯一性由定理 F.17 得出。$\square$

**推论 F.19（保单位的序同构）。** 使 $\Phi$ 与 $\Phi^{-1}$ 都是正的保单位线性双射
$\Phi:\mathcal M\to\mathcal N$ 恰好是定理 F.18 中的映射。这样的 $\Phi$ 是 $2$-正的，当且仅当它是
完全正的，当且仅当 $S$ 不包含 $k_c\ge2$ 的块，当且仅当 $\Phi$ 是 $*$-同构。

*证明。* 第一句就是定理 F.15 与 F.18。若 $S$ 不含 $k_c\ge2$ 的块，则每个块映射都是
$x\mapsto U_cxU_c^*$（对 $k_c=1$，$x^T=x$），所以 $\Phi$ 是完全正的（命题 6.4）并且是乘性的。若
$c\in S$，$k_c\ge2$，且 $\Phi$ 是 $2$-正的，则块映射 $\Phi_{\sigma(c)c}(a)=U_ca^TU_c^*$ 也是如此
（命题 6.4），并且 $a\mapsto U_c^*\Phi_{\sigma(c)c}(a)U_c=a^T$ 也是，因为 $[A_{ij}]\ge0$ 蕴含
$[U^*A_{ij}U]=(1\otimes U)^*[A_{ij}](1\otimes U)\ge0$；这与命题 6.1 矛盾。由定理 F.17，这样的 $\Phi$ 也
不是乘性的。$\square$

所以，在相差酉矩阵的意义下，保单位的序同构由它交换哪些相等的块以及转置哪些块决定，
而由推论 F.19，被转置的块恰恰是妨碍完全正性的因素。代码构造 $M_3$ 的若尔当
自同构 $x\mapsto Ux^TU^*$，并按照定理 F.17 证明的步骤恢复 $U$（在相差一个相位的意义下）：

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

## 习题

1. 证明带有逐坐标乘积与标准内积的 $\mathbb R^N$ 是欧几里得若尔当代数，
   其平方锥是 $\mathbb R^N_{\ge0}$，并由推论 F.19 推出
   $\mathbb C^N=\bigoplus_{c=1}^NM_1$ 的保单位序同构就是坐标的置换。
2. （**仅有正性是不够的。**）对 $0<t<1$，在 $M_n$（$n\ge2$）上令 $\Psi_t(x)=(1-t)x+\frac tn\operatorname{Tr}(x)1$。证明
   $\Psi_t$ 是保单位的正双射，$\Psi_t^{-1}$ 不是正的，并且 $\Psi_t$ 不把
   投影映为投影。

<details>
<summary>解答</summary>

1. 该乘积是交换且结合的，从而是若尔当乘积，单位元为 $(1,\dots,1)$，而 $L(x)$ 是对角矩阵
   $\operatorname{diag}(x)$，它是对称的。平方元就是坐标非负的向量。由所有 $k_c=1$ 情形下的
   推论 F.19，$\Phi(x)_{\sigma(c)}=U_cx_cU_c^*=x_c$（$M_1$ 的酉元是相位因子），所以
   $\Phi$ 置换坐标；反之，置换都是保单位的序同构。

2. $\Psi_t(1)=(1-t)1+t1=1$；对 $x\ge0$，两项都 $\ge0$。由于 $\operatorname{Tr}\Psi_t(x)=\operatorname{Tr}x$，
   $\Psi_t^{-1}(y)=\frac1{1-t}\big(y-\frac tn\operatorname{Tr}(y)1\big)$。对秩一投影 $p$，$\Psi_t^{-1}(p)$
   在 $\ker p$ 上有特征值 $-\frac{t}{n(1-t)}<0$，而 $\Psi_t(p)$ 有特征值 $\frac tn\notin\{0,1\}$。这
   与定理 F.15 相容：正是逆映射的正性迫使投影被映为投影。

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

**返回：** [第 13 章](13_jordan_algebras.md)，13.4–13.6 节：定理 F.12 证明了定理 13.8，定理 F.15
对代数 $\bigoplus_cM_{k_c}(\mathbb C)$ 证明了定理 13.10，定理 F.17 把命题 13.11 推广到所有
$n$，而推论 F.7 把命题 13.14 推广到每个欧几里得若尔当代数。
