[English](../B_schur_weyl.md) · [Español](../es/B_schur_weyl.md) · [Français](../fr/B_schur_weyl.md) · [Deutsch](../de/B_schur_weyl.md) · **中文** · [日本語](../ja/B_schur_weyl.md) · [Русский](../ru/B_schur_weyl.md)

# 附录 B. Schur–Weyl 对偶与 Weingarten 演算

酉群 $`U(d)`$ 通过 $`V\mapsto V^{\otimes n}`$ 作用在 $`(\mathbb C^d)^{\otimes n}`$ 上，对称群 $`S_n`$ 则通过置换各张量因子作用在
其上。这两族算子中的每一族都是另一族的交换子（**Schur–Weyl 对偶**）。本附录完整地证明这一对偶，并由它导出
任意阶 Haar 酉矩阵的矩的 Weingarten 公式。在此过程中，它证明了第 4、5、9 和 11 章从这里引用的结论：对所有 $`L`$
成立的定理 4.14(b)，$`S_L`$ 在总自旋的重数空间上的表示的不可约性，4.5 节的 Catalan 数，定理 5.12，9.3 节
所用的 $`SU(2)`$ 的 Wigner–Eckart 定理，以及第 11 章所用的正交群的 Haar 测度。

证明路线是初等的：一个极化恒等式把置换的交换子等同于所有 $`A^{\otimes n}`$（$`A\in M_d(\mathbb C)`$）的张成；
一个 Cayley 变换论证表明酉矩阵 $`V^{\otimes n}`$ 张成同一个空间；双交换子定理（定理 2.2）再把这一点变成对偶。
随后 Weingarten 公式之所以成立，是因为 Haar 平均就是到交换子上的正交投影。

本附录的数值检验见[教程，附录 B](../../tutorial/zh/B_schur_weyl.md)。

在本附录中，$`d,n\ge1`$，$`\mathcal H_n=(\mathbb C^d)^{\otimes n}`$，其规范正交基为
$`e_{\mathbf i}=e_{i_1}\otimes\dots\otimes e_{i_n}`$，其中多重指标 $`\mathbf i=(i_1,\dots,i_n)\in[d]^n`$ 被看作映射
$`[n]\to[d]`$，$`m\mapsto i_m`$（$`[d]=\{1,\dots,d\}`$）。$`\mathrm{End}\,\mathcal H_n=M_{d^n}(\mathbb C)`$ 带有
Hilbert–Schmidt 内积 $`\langle X,Y\rangle=\mathrm{Tr}(X^*Y)`$。对 $`X\in M_d`$ 和 $`1\le m\le n`$，记
$`X^{(m)}=1\otimes\dots\otimes X\otimes\dots\otimes1`$（$`X`$ 位于第 $`m`$ 个位置），以及 $`d\Pi(X)=\sum_{m=1}^nX^{(m)}`$。记
```math
\mathcal T_n=\mathrm{span}\{A^{\otimes n}:A\in M_d(\mathbb C)\},\qquad
\mathcal U_n=\mathrm{span}\{V^{\otimes n}:V\in U(d)\}.
```
Kronecker 积的乘积按因子逐个计算，$`(A_1\otimes\dots\otimes A_n)(B_1\otimes\dots\otimes B_n)=A_1B_1\otimes\dots\otimes A_nB_n`$，
并且 $`(A_1\otimes\dots\otimes A_n)^*=A_1^*\otimes\dots\otimes A_n^*`$。置换按映射复合，
$`(\sigma\tau)(m)=\sigma(\tau(m))`$，$`\#\text{cycles}(\sigma)`$ 表示 $`\sigma`$ 的循环个数（包括不动点）。

<a id="s-B-1"></a>

## B.1 张量因子的置换

对 $`\sigma\in S_n`$，令 $`\pi(\sigma)e_{\mathbf i}=e_{\mathbf i\circ\sigma^{-1}}`$，即位于第 $`m`$ 个槽位的基向量移到第
$`\sigma(m)`$ 个槽位。这是 `decompose.permutation_representation` 以及第 2 章和第 4 章的约定。

*库中：* 算子 $`\pi(\sigma)`$ 由 `decompose.permutation_representation` 构造——[教程，B.3 节](../../tutorial/zh/B_schur_weyl.md#s-B-3)。

<a id="r-B-1"></a>

**引理 B.1.** (a) $`\pi`$ 是 $`S_n`$ 的酉表示，并且对所有 $`x_m\in\mathbb C^d`$ 有
$`\pi(\sigma)(x_1\otimes\dots\otimes x_n)=x_{\sigma^{-1}(1)}\otimes\dots\otimes x_{\sigma^{-1}(n)}`$。
(b) 对所有 $`A_m\in M_d`$ 有 $`\pi(\sigma)(A_1\otimes\dots\otimes A_n)\pi(\sigma)^*=A_{\sigma^{-1}(1)}\otimes\dots\otimes A_{\sigma^{-1}(n)}`$。
特别地，$`\pi(\sigma)`$ 与每个 $`A^{\otimes n}`$ 及每个 $`d\Pi(X)`$ 可交换，并且
$`\pi(\sigma)X^{(m)}\pi(\sigma)^*=X^{(\sigma(m))}`$。
(c) $`\mathrm{Tr}\,\pi(\sigma)=d^{\#\text{cycles}(\sigma)}`$，$`\langle\pi(\sigma),\pi(\tau)\rangle=d^{\#\text{cycles}(\sigma^{-1}\tau)}`$。
循环个数在共轭和取逆下不变，并且两个置换共轭当且仅当它们具有相同的循环型。
(d) 这 $`n!`$ 个算子 $`\pi(\sigma)`$ 线性无关当且仅当 $`n\le d`$。

*证明。* (a) $`\mathbf i\mapsto\mathbf i\circ\sigma^{-1}`$ 是 $`[d]^n`$ 的双射，所以 $`\pi(\sigma)`$ 置换一组规范正交基，
因而是酉的。$`\pi(\sigma)\pi(\tau)e_{\mathbf i}=\pi(\sigma)e_{\mathbf i\circ\tau^{-1}}=e_{\mathbf i\circ\tau^{-1}\circ\sigma^{-1}}=\pi(\sigma\tau)e_{\mathbf i}`$。
展开 $`x_m=\sum_ax_m(a)e_a`$，得 $`x_1\otimes\dots\otimes x_n=\sum_{\mathbf i}\prod_mx_m(i_m)\,e_{\mathbf i}`$；作用 $`\pi(\sigma)`$ 并
代入 $`\mathbf j=\mathbf i\circ\sigma^{-1}`$，得到 $`\sum_{\mathbf j}\prod_mx_m(j_{\sigma(m)})\,e_{\mathbf j}=\sum_{\mathbf j}\prod_{p}x_{\sigma^{-1}(p)}(j_p)\,e_{\mathbf j}`$，
这正是所断言的乘积向量。
(b) 由 (a)，在乘积向量上：$`\pi(\sigma)^*=\pi(\sigma^{-1})`$ 把 $`x_{\sigma(m)}`$ 放入第 $`m`$ 个槽位，算子在那里作用
$`A_m`$，然后 $`\pi(\sigma)`$ 把第 $`m`$ 个槽位移到第 $`\sigma(m)`$ 个槽位；所以第 $`p`$ 个槽位最终是 $`A_{\sigma^{-1}(p)}x_p`$。乘积
向量张成 $`\mathcal H_n`$。对 $`X^{(m)}`$，因子 $`X`$ 从第 $`m`$ 个槽位移到第 $`\sigma(m)`$ 个槽位，所以 $`d\Pi(X)=\sum_mX^{(m)}`$ 保持不变。
(c) $`\pi(\sigma)`$ 是置换矩阵，而 $`e_{\mathbf i}`$ 不动当且仅当对所有 $`p`$ 有 $`i_{\sigma^{-1}(p)}=i_p`$，即当且仅当
$`\mathbf i`$ 在 $`\sigma`$ 的每个循环上取常值：这样的 $`\mathbf i`$ 共有 $`d^{\#\text{cycles}(\sigma)}`$ 个。于是
$`\langle\pi(\sigma),\pi(\tau)\rangle=\mathrm{Tr}(\pi(\sigma^{-1})\pi(\tau))=\mathrm{Tr}\,\pi(\sigma^{-1}\tau)`$。
$`\kappa\rho\kappa^{-1}`$ 的循环是 $`\rho`$ 的循环在 $`\kappa`$ 下的像，而 $`\rho^{-1}`$ 沿相反方向遍历同样的循环。反之，若
$`\rho,\rho'`$ 具有相同的循环型，把它们的循环按长度排序、上下对齐写成两行，并令 $`\kappa`$ 把第一行的每个元素映到它
下方的元素；则 $`\kappa\rho\kappa^{-1}=\rho'`$。
(d) 若 $`n\le d`$，令 $`\mathbf u=(1,2,\dots,n)`$。对不同的 $`\sigma`$，向量 $`\pi(\sigma)e_{\mathbf u}=e_{\mathbf u\circ\sigma^{-1}}`$ 是
不同的基向量，所以把 $`\sum_\sigma c_\sigma\pi(\sigma)=0`$ 作用于 $`e_{\mathbf u}`$ 即得所有 $`c_\sigma=0`$。若 $`n>d`$，
令 $`\mathcal A=\sum_\sigma\mathrm{sgn}(\sigma)\pi(\sigma)`$。每个 $`\mathbf i\in[d]^n`$ 都对某对 $`a\neq b`$ 有 $`i_a=i_b`$；对
对换 $`t=(a\,b)`$，有 $`\pi(t)e_{\mathbf i}=e_{\mathbf i}`$ 且 $`\mathcal A\pi(t)=\sum_\sigma\mathrm{sgn}(\sigma)\pi(\sigma t)=-\mathcal A`$。
因此对所有 $`\mathbf i`$ 有 $`\mathcal Ae_{\mathbf i}=\mathcal A\pi(t)e_{\mathbf i}=-\mathcal Ae_{\mathbf i}=0`$：$`\mathcal A=0`$ 是一个非平凡的
线性关系。$`\square`$

*计算：* 引理 B.1(a) 与 (d)——[教程，B.3 节](../../tutorial/zh/B_schur_weyl.md#s-B-3)。

<a id="s-B-2"></a>

## B.2 极化：置换的交换子

<a id="r-B-2"></a>

**引理 B.2（极化）。** 对 $`A_1,\dots,A_n\in M_d(\mathbb C)`$，记 $`|\varepsilon|=\sum_m\varepsilon_m`$，有
```math
\sum_{\sigma\in S_n}A_{\sigma(1)}\otimes\dots\otimes A_{\sigma(n)}=
\sum_{\varepsilon\in\{0,1\}^n}(-1)^{n-|\varepsilon|}\Big(\sum_{m=1}^n\varepsilon_mA_m\Big)^{\otimes n}.
```

*证明。* 由多重线性，$`\big(\sum_m\varepsilon_mA_m\big)^{\otimes n}=\sum_f\varepsilon_{f(1)}\cdots\varepsilon_{f(n)}\,A_{f(1)}\otimes\dots\otimes A_{f(n)}`$，
求和遍及所有映射 $`f:[n]\to[n]`$。由于当 $`r\ge1`$ 时 $`\varepsilon_m^r=\varepsilon_m`$，右边中项 $`f`$ 的系数为
```math
\sum_{\varepsilon}\prod_{m=1}^n(-1)^{1-\varepsilon_m}\prod_{m\in f([n])}\varepsilon_m
=\prod_{m\notin f([n])}\Big(\sum_{\varepsilon_m=0,1}(-1)^{1-\varepsilon_m}\Big)\prod_{m\in f([n])}\Big(\sum_{\varepsilon_m=0,1}(-1)^{1-\varepsilon_m}\varepsilon_m\Big).
```
第一类因子是 $`-1+1=0`$，第二类是 $`1`$。所以当 $`f`$ 是满射（即置换）时系数为 $`1`$，否则为 $`0`$。$`\square`$

*计算：* $`n=3`$ 时的引理 B.2——[教程，B.3 节](../../tutorial/zh/B_schur_weyl.md#s-B-3)。

$`\mathrm{End}\,\mathcal H_n`$ 上的**对称化算子**是 $`\mathcal S(X)=\frac1{n!}\sum_\sigma\pi(\sigma)X\pi(\sigma)^*`$。

<a id="r-B-3"></a>

**命题 B.3.** $`\mathcal S`$ 是 $`\mathrm{End}\,\mathcal H_n`$ 到 $`\pi(S_n)'`$ 上的正交投影，并且
```math
\pi(S_n)'=\mathcal T_n,\qquad\dim\pi(S_n)'=\binom{d^2+n-1}n .
```

*证明。* $`\mathcal S`$ 是有限群 $`S_n`$ 在表示 $`\pi`$ 中的 twirl，因此是到 $`\pi(S_n)'`$ 上的正交投影（取
$`\pi_1=\pi_2=\pi`$ 的定理 12.5；另见 4.2 节）。由引理 B.1(b)，$`\mathcal T_n\subseteq\pi(S_n)'`$。反之，
$`\mathrm{End}\,\mathcal H_n`$ 由矩阵单位的乘积 $`A_1\otimes\dots\otimes A_n`$ 张成，而由引理 B.1(b) 及代换 $`\sigma\to\sigma^{-1}`$，
$`\mathcal S(A_1\otimes\dots\otimes A_n)=\frac1{n!}\sum_\sigma A_{\sigma(1)}\otimes\dots\otimes A_{\sigma(n)}`$，由引理 B.2 它属于
$`\mathcal T_n`$。所以 $`\pi(S_n)'=\mathcal S(\mathrm{End}\,\mathcal H_n)\subseteq\mathcal T_n`$。

*维数。* 令 $`P=[d]\times[d]`$，对 $`\alpha\in P^n`$ 令 $`E_\alpha=E_{\alpha_1}\otimes\dots\otimes E_{\alpha_n}`$，其中矩阵单位为
$`E_{(a,b)}=|e_a\rangle\langle e_b|`$：这是 $`\mathrm{End}\,\mathcal H_n`$ 的一组基。由引理 B.1(b)，
$`\mathcal S(E_\alpha)=\frac1{n!}\sum_\sigma E_{\alpha\circ\sigma^{-1}}`$，对轨道 $`O_\alpha=\{\alpha\circ\sigma\}`$，它等于
$`|O_\alpha|^{-1}\sum_{\beta\in O_\alpha}E_\beta`$（轨道的每个元素都恰对同样数目 $`n!/|O_\alpha|`$ 个 $`\sigma`$ 出现）。不同轨道
上的和涉及互不相交的基元素集合，所以 $`\mathcal S`$ 的值域有一组以轨道为指标的基。两个元组属于同一轨道当且仅当
$`P`$ 的每个元素在两者中出现的次数相同：轨道就是由 $`P`$ 的 $`N=d^2`$ 个元素构成的大小为 $`n`$ 的多重集。把 $`P`$ 排序为
$`1,\dots,N`$，这样的多重集就是一个序列 $`1\le c_1\le\dots\le c_n\le N`$，而 $`c_m\mapsto c_m+m-1`$ 是到 $`[N+n-1]`$ 中严格
递增序列（即 $`n`$ 元子集）的双射；它们共有 $`\binom{N+n-1}n`$ 个。$`\square`$

*计算：* 对较小的 $`d`$ 和 $`n`$ 验证 $`\dim\pi(S_n)'=\binom{d^2+n-1}n`$——[教程，B.3 节](../../tutorial/zh/B_schur_weyl.md#s-B-3)。

<a id="s-B-3"></a>

## B.3 Schur–Weyl 对偶

<a id="r-B-4"></a>

**引理 B.4（在酉群上为零的多项式）。** 设 $`p:M_d(\mathbb C)\to\mathbb C`$ 是关于 $`d^2`$ 个矩阵元的多项式，且是
$`n\ge1`$ 次齐次的。若对所有 $`V\in U(d)`$ 有 $`p(V)=0`$，则 $`p=0`$。

*证明。* *第 1 步（Cayley 变换）。* 对 $`A\in M_d`$ 令 $`q(A)=p\big((1+iA)\mathrm{adj}(1-iA)\big)`$，其中
$`\mathrm{adj}`$ 表示伴随矩阵（adjugate）；$`q`$ 是关于 $`A`$ 的矩阵元的多项式。设 $`H`$ 为厄米矩阵。它的特征值
都是实数（推论 1.4），所以 $`1-iH`$ 可逆，$`\mathrm{adj}(1-iH)=\det(1-iH)(1-iH)^{-1}`$（Cramer 法则），并由
齐次性得 $`q(H)=\det(1-iH)^n\,p(C)`$，其中 $`C=(1+iH)(1-iH)^{-1}`$。矩阵 $`1\pm iH`$ 及其逆都是 $`H`$ 的函数，因而彼此
可交换（定理 1.5），并且 $`(1+iH)^*=1-iH`$，所以
$`C^*C=(1+iH)^{-1}(1-iH)(1+iH)(1-iH)^{-1}=1`$。因此 $`q(H)=0`$。

*第 2 步（从厄米矩阵到所有矩阵）。* 取厄米矩阵构成的实向量空间的一组基 $`B_1,\dots,B_{d^2}`$。它也是
$`M_d(\mathbb C)`$ 在 $`\mathbb C`$ 上的一组基，因为每个 $`A`$ 都等于 $`\frac{A+A^*}2+i\frac{A-A^*}{2i}`$，且维数相同。
$`r(t)=q(\sum_kt_kB_k)`$ 是 $`t\in\mathbb C^{d^2}`$ 的多项式，并在 $`\mathbb R^{d^2}`$ 上为零。这样的多项式的所有系数
都是 $`0`$：对一个变量，它有无穷多个根；对 $`N`$ 个变量，写 $`r=\sum_jr_j(t')t_N^j`$——对每个实的 $`t'`$，关于 $`t_N`$ 的
多项式在 $`\mathbb R`$ 上为零，所以每个 $`r_j`$ 在 $`\mathbb R^{N-1}`$ 上为零，再用归纳法即可。因此对每个
$`A=\sum_kt_kB_k\in M_d`$ 有 $`q(A)=r(t)=0`$。

*第 3 步（反演变换）。* 设 $`B\in M_d`$ 使 $`1+B`$ 可逆，并令 $`A=-i(1+B)^{-1}(B-1)`$。则
$`1-iA=1-(1+B)^{-1}(B-1)=2(1+B)^{-1}`$，$`1+iA=(1+B)^{-1}2B`$；由于 $`B`$ 与 $`1+B`$ 可交换，
$`(1+iA)(1-iA)^{-1}=(1+B)^{-1}B(1+B)=B`$。所以 $`0=q(A)=\det(1-iA)^np(B)`$，而 $`\det(1-iA)\neq0`$：$`p(B)=0`$。

*第 4 步。* 对任意 $`B\neq0`$ 和 $`0<t<1/\|B\|`$，$`-1`$ 不是 $`tB`$ 的特征值（若 $`tBv=-v`$，则除非 $`v=0`$，否则
$`\|v\|=t\|Bv\|\le t\|B\|\|v\|<\|v\|`$），所以 $`t^np(B)=p(tB)=0`$。最后，由齐次性 $`p(0)=0`$。$`\square`$

<a id="r-B-5"></a>

**引理 B.5（Newton 恒等式）。** 设 $`x_1,\dots,x_n`$ 是某个代数中两两可交换的元素，$`x_S=\prod_{m\in S}x_m`$，
$`e_k=\sum_{|S|=k}x_S`$（$`e_0=1`$），$`p_j=\sum_mx_m^j`$。则对 $`1\le k\le n`$ 有 $`k\,e_k=\sum_{j=1}^k(-1)^{j-1}e_{k-j}\,p_j`$。

*证明。* 对 $`1\le j\le k`$，令 $`a_j=\sum x_Sx_m^j`$，求和遍及满足 $`|S|=k-j`$ 且 $`m\notin S`$ 的对 $`(S,m)`$，并令
$`a_{k+1}=0`$。在 $`e_{k-j}p_j=\sum_{|S|=k-j}\sum_mx_Sx_m^j`$ 中，$`m\notin S`$ 的项给出 $`a_j`$，而 $`m\in S`$ 的项写成
$`S'=S\setminus\{m\}`$ 的形式后为 $`x_{S'}x_m^{j+1}`$，它们给出 $`a_{j+1}`$。所以交错和裂项相消为 $`a_1`$，而
$`a_1=\sum_{|T|=k}|T|\,x_T=k\,e_k`$（每个 $`k`$ 元集 $`T`$ 由 $`m\in T`$ 的 $`k`$ 种选择产生）。$`\square`$

*计算：* $`x_m=A^{(m)}`$、$`n=4`$ 时的引理 B.5，它由各个 $`d\Pi(A^j)`$ 重建 $`A^{\otimes4}`$——[教程，B.3 节](../../tutorial/zh/B_schur_weyl.md#s-B-3)。

<a id="r-B-6"></a>

**定理 B.6（Schur–Weyl 对偶）。** 设 $`d,n\ge1`$。
(a) $`\pi(S_n)'=\mathcal T_n=\mathcal U_n=\mathrm{span}\{V^{\otimes n}:V\in SU(d)\}`$。这个空间是一个 $`*`$-代数，所以它就是
由 $`V^{\otimes n}`$（$`V\in U(d)`$）生成的代数。
(b) 对每个张成 $`M_d`$ 的集合 $`\mathcal X\subseteq M_d`$，由
$`\{d\Pi(X):X\in\mathcal X\}`$ 生成的含单位元的代数是 $`\mathcal T_n`$。
(c) $`\{V^{\otimes n}:V\in U(d)\}'=\{V^{\otimes n}:V\in SU(d)\}'=\{d\Pi(X):X\in M_d\}'=\mathrm{span}\,\pi(S_n)`$。

*证明。* (a) $`\pi(S_n)'=\mathcal T_n`$ 就是命题 B.3，并且 $`\mathcal U_n\subseteq\mathcal T_n`$。假设 $`\mathcal U_n\neq\mathcal T_n`$。则
某个 $`A_0^{\otimes n}\notin\mathcal U_n`$；令 $`Y\neq0`$ 为 $`A_0^{\otimes n}`$ 正交于 $`\mathcal U_n`$ 的分量。函数
$`p(A)=\langle Y,A^{\otimes n}\rangle`$ 是关于 $`A`$ 的矩阵元的 $`n`$ 次齐次多项式（$`A^{\otimes n}`$ 的每个矩阵元都是 $`A`$ 的
$`n`$ 个矩阵元的乘积），它在 $`U(d)`$ 上为零，而 $`p(A_0)=\|Y\|^2\neq0`$，与引理 B.4 矛盾。
对 $`V\in U(d)`$，取 $`\theta`$ 使 $`\det V=e^{i\theta}`$；则 $`W=e^{-i\theta/d}V\in SU(d)`$ 且 $`V^{\otimes n}=e^{in\theta/d}W^{\otimes n}`$，
所以两个张成相同。$`\mathcal U_n`$ 包含 $`1=1^{\otimes n}`$、$`V^{\otimes n}W^{\otimes n}=(VW)^{\otimes n}`$ 和
$`(V^{\otimes n})^*=(V^*)^{\otimes n}`$，所以它是 $`*`$-代数，而每个包含这些 $`V^{\otimes n}`$ 的代数都包含它们的张成。

(b) 设 $`\mathcal A`$ 为该代数。$`d\Pi`$ 是线性的，所以 $`\mathcal A`$ 对所有 $`X\in M_d`$ 包含 $`d\Pi(X)`$。由引理 B.1(b)，
每个 $`d\Pi(X)`$ 都属于含单位元的代数 $`\pi(S_n)'=\mathcal T_n`$，所以 $`\mathcal A\subseteq\mathcal T_n`$。反之，固定 $`A\in M_d`$，
并把引理 B.5 应用于两两可交换的 $`x_m=A^{(m)}`$：则 $`e_n=A^{\otimes n}`$，$`p_j=d\Pi(A^j)\in\mathcal A`$，并且对 $`k`$
归纳可知每个 $`e_k=\frac1k\sum_j(-1)^{j-1}e_{k-j}p_j`$ 都属于 $`\mathcal A`$。所以 $`A^{\otimes n}\in\mathcal A`$，
$`\mathcal T_n\subseteq\mathcal A`$。

(c) 一个集合的交换子等于它所生成的代数的交换子。由 (a) 和 (b)，这三个集合都生成
$`\mathcal T_n=\pi(S_n)'`$，其交换子是 $`\pi(S_n)''`$。$`\pi(S_n)`$ 的张成是一个 $`*`$-代数
（$`\pi(\sigma)\pi(\tau)=\pi(\sigma\tau)`$，$`\pi(\sigma)^*=\pi(\sigma^{-1})`$，$`\pi(e)=1`$），所以由定理 2.2，
$`\pi(S_n)''=\mathrm{span}\,\pi(S_n)`$。$`\square`$

*计算：* 对较小的 $`d`$ 和 $`n`$ 验证定理 B.6(a) 与 (c) 中的维数——[教程，B.3 节](../../tutorial/zh/B_schur_weyl.md#s-B-3)。

<a id="r-B-7"></a>

**推论 B.7（定理 4.14(b)）。** 对 $`L`$ 个自旋 $`\frac12`$ 的总自旋 $`S_a=\sum_is_a^{(i)}`$，由 $`S_x,S_y,S_z`$ 生成的代数是
$`\{S_a\}''=\pi(S_L)'`$，并且 $`\{S_a\}'=\mathrm{span}\,\pi(S_L)`$。

*证明。* $`M_2=\mathrm{span}\{1,s_x,s_y,s_z\}`$（Pauli 矩阵与 $`1`$ 构成一组基），$`d\Pi(1)=L\cdot1`$，$`d\Pi(s_a)=S_a`$。
由定理 B.6(b)，由 $`S_x,S_y,S_z`$ 生成的含单位元的代数 $`\mathcal A`$ 是 $`\mathcal T_L=\pi(S_L)'`$。由于 $`S_a`$ 是厄米的，
$`\mathcal A`$ 是 $`*`$-代数（生成元乘积的伴随就是顺序颠倒的乘积），所以由定理 2.2，
$`\mathcal A=\mathcal A''\supseteq\{S_a\}''`$；又因为 $`\{S_a\}''`$ 是包含 $`S_a`$ 的含单位元的代数，$`\{S_a\}''\supseteq\mathcal A`$。
最后 $`\{S_a\}'=\{S_a\}'''=\pi(S_L)''=\mathrm{span}\,\pi(S_L)`$（2.1 节与定理 B.6(c)）。$`\square`$

由 $`d=2`$ 时的定理 B.6(a)，同一个代数由转动 $`V^{\otimes L}`$（$`V\in SU(2)`$）张成；这不需要 4.5 节中提到的
$`SU(2)`$ 的指数参数化。

<a id="s-B-4"></a>

## B.4 $`(\mathbb C^d)^{\otimes n}`$ 的分解

<a id="r-B-8"></a>

**定理 B.8.** 存在 $`C\ge1`$、满足 $`\sum_ck_cm_c=d^n`$ 的整数 $`k_c,m_c\ge1`$，以及酉映射
$`W:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathcal H_n`$，使得
```math
W^*\pi(\sigma)W=\bigoplus_{c=1}^C\rho_c(\sigma)\otimes1_{m_c},\qquad W^*V^{\otimes n}W=\bigoplus_{c=1}^C1_{k_c}\otimes\tau_c(V)
\qquad(\sigma\in S_n,\ V\in U(d)),
```
其中 $`\rho_c`$ 是 $`S_n`$ 在 $`\mathbb C^{k_c}`$ 上两两不等价的不可约酉表示，$`\tau_c`$ 是 $`U(d)`$ 在 $`\mathbb C^{m_c}`$ 上
两两不等价的不可约酉表示。此外
$`W^*\mathrm{span}\,\pi(S_n)W=\bigoplus_cM_{k_c}\otimes1_{m_c}`$，$`W^*\mathcal U_nW=\bigoplus_c1_{k_c}\otimes M_{m_c}`$，
```math
\sum_ck_c^2=\dim\mathrm{span}\,\pi(S_n),\qquad\sum_cm_c^2=\binom{d^2+n-1}n,
```
并且当 $`n\le d`$ 时 $`\sum_ck_c^2=n!`$。

*证明。* 对有限群 $`S_n`$ 应用命题 2.7，并对 $`\mathcal M=\mathrm{span}\,\pi(S_n)=\pi(S_n)''`$ 应用定理 2.6，即得
$`W`$、具有所述性质的 $`\rho_c`$，以及 $`W^*\pi(S_n)'W=\bigoplus_c1_{k_c}\otimes M_{m_c}`$。由定理 B.6(a)，
$`\pi(S_n)'=\mathcal U_n`$。所以对某些 $`\tau_c(V)\in M_{m_c}`$ 有 $`W^*V^{\otimes n}W=\bigoplus_c1\otimes\tau_c(V)`$，并且由于
$`V\mapsto V^{\otimes n}`$ 是酉表示，$`V\mapsto\tau_c(V)`$ 也是。由于 $`V^{\otimes n}`$ 张成 $`\mathcal U_n`$，线性映射
$`\sum_la_lV_l^{\otimes n}\mapsto(\sum_la_l\tau_c(V_l))_c`$ 映满 $`\bigoplus_cM_{m_c}`$；$`\tau_c`$ 的不可约性与互不等价性的
证明与命题 2.7 的证明逐字相同。维数来自定理 2.6、命题 B.3 和引理 B.1(d)。$`\square`$

*计算：* $`(\mathbb C^3)^{\otimes3}`$ 上 $`S_3`$ 的定理 B.8，分解由 `decompose.from_group` 完成——[教程，B.4 节](../../tutorial/zh/B_schur_weyl.md#s-B-4)。

对 $`d=2`$，各扇区由总自旋来标识：

<a id="r-B-9"></a>

**推论 B.9（$`L`$ 个自旋 $`\frac12`$）。** 设 $`V`$ 是对 $`L`$ 个自旋 $`\frac12`$ 的总自旋满足
$`V^*S_aV=\bigoplus_SS^{(S)}_a\otimes1_{m_S}`$ 的酉矩阵（推论 4.12，命题 4.13）。则
$`V^*\pi(\sigma)V=\bigoplus_S1_{2S+1}\otimes\rho_S(\sigma)`$，其中 $`\rho_S`$ 是 $`S_L`$ 的两两不等价的不可约酉表示，
其维数为
```math
m_S=\binom Lj-\binom L{j-1},\qquad j=\tfrac L2-S\in\{0,1,\dots,\lfloor L/2\rfloor\}.
```
特别地，$`\dim\mathrm{span}\,\pi(S_L)=\sum_Sm_S^2`$，$`\dim\pi(S_L)'=\sum_S(2S+1)^2=\binom{L+3}3`$。

*证明。* 由定理 4.3，$`V^*\{S_a\}'V=\bigoplus_S1_{2S+1}\otimes M_{m_S}`$，而由推论 B.7，$`\{S_a\}'=\mathrm{span}\,\pi(S_L)`$。所以
$`V^*\pi(\sigma)V=\bigoplus_S1\otimes\rho_S(\sigma)`$，其中 $`\rho_S`$ 是酉表示，并且映射
$`\sum_\sigma a_\sigma\pi(\sigma)\mapsto(\sum_\sigma a_\sigma\rho_S(\sigma))_S`$ 映满 $`\bigoplus_SM_{m_S}`$；不可约性与互不等价性
的证明同命题 2.7。$`m_S`$ 的值来自命题 4.13。由定理 4.3，
$`\dim\{S_a\}'=\sum_Sm_S^2`$，$`\dim\{S_a\}''=\sum_S(2S+1)^2`$；由推论 B.7，$`\{S_a\}'=\mathrm{span}\,\pi(S_L)`$，
$`\{S_a\}''=\pi(S_L)'`$，而由 $`d=2`$ 时的命题 B.3，后者的维数为 $`\binom{4+L-1}L=\binom{L+3}3`$。$`\square`$

*计算：* $`L=5,6`$ 时的数对 $`(m_S,2S+1)`$——[教程，B.4 节](../../tutorial/zh/B_schur_weyl.md#s-B-4)。

<a id="r-B-10"></a>

**推论 B.10（Catalan 数）。** 对 $`L\ge1`$，$`\displaystyle\sum_Sm_S^2=\binom{2L}L\frac1{L+1}`$；这就是
$`(\mathbb C^2)^{\otimes L}`$ 上的 $`\dim\{S_a\}'=\dim\mathrm{span}\,\pi(S_L)`$。

*证明。* 记 $`N=2^L`$，并令 $`\mathrm{vec}:M_N\to\mathbb C^N\otimes\mathbb C^N`$ 为线性双射
$`|x\rangle\langle y|\mapsto x\otimes\bar y`$。则 $`\mathrm{vec}(AX-XA)=(A\otimes1-1\otimes A^T)\mathrm{vec}\,X`$，因为
$`A|x\rangle\langle y|B=|Ax\rangle\langle B^*y|\mapsto Ax\otimes B^T\bar y`$。所以 $`\{S_a\}'`$ 被映到
$`S_a\otimes1-1\otimes S_a^T`$（$`a=x,y,z`$）的公共核上，其中 $`S_a^T=\sum_i(s_a^T)^{(i)}`$。由 $`s_x^T=s_x`$，$`s_y^T=-s_y`$，$`s_z^T=s_z`$ 以及
$`\sigma_ys_x\sigma_y=-s_x`$，$`\sigma_ys_y\sigma_y=s_y`$，$`\sigma_ys_z\sigma_y=-s_z`$ 得 $`\sigma_y(-s_a^T)\sigma_y=s_a`$，所以酉矩阵
$`Y=1\otimes\sigma_y^{\otimes L}`$ 给出 $`Y(S_a\otimes1-1\otimes S_a^T)Y^*=S_a\otimes1+1\otimes S_a`$，即
$`(\mathbb C^2)^{\otimes2L}=\mathbb C^N\otimes\mathbb C^N`$ 上 $`2L`$ 个自旋 $`\frac12`$ 的总自旋。在 $`2L`$ 个自旋的推论 4.12 的分解中，
被三个分量全部零化的向量也被 $`\mathbf S^2`$ 零化，而它在自旋 $`S\ge\frac12`$ 的直和项上等于 $`S(S+1)\ne0`$；
反之，三个分量在自旋 $`0`$ 的直和项上为零。所以公共核就是这个直和项，其维数为
$`m_0(2L)=\binom{2L}L-\binom{2L}{L-1}=\binom{2L}L\big(1-\frac L{L+1}\big)`$（命题 4.13）。维数的第一个等式是定理 4.3，
第二个是推论 B.7。$`\square`$

*计算：* $`L=5,6`$ 时的推论 B.10——[教程，B.4 节](../../tutorial/zh/B_schur_weyl.md#s-B-4)。

<a id="s-B-5"></a>

## B.5 Weingarten 演算

<a id="r-B-11"></a>

**引理 B.11（Haar twirl 是到交换子上的正交投影）。** 设 $`U`$ 在 $`U(d)`$ 上服从 Haar 分布，并对
$`X\in\mathrm{End}\,\mathcal H_n`$ 令 $`\Phi(X)=\mathbb E[U^{\otimes n}X(U^{\otimes n})^*]`$。则 $`\Phi`$ 是
$`\mathrm{End}\,\mathcal H_n`$ 到 $`\{V^{\otimes n}:V\in U(d)\}'=\mathrm{span}\,\pi(S_n)`$ 上的正交投影。

*证明。* $`U^{\otimes n}X(U^{\otimes n})^*`$ 的矩阵元是 $`U`$ 和 $`\bar U`$ 的矩阵元的多项式，而这些矩阵元的绝对值以 $`1`$
为界，所以 $`\Phi`$ 是定义良好的线性映射。对 $`V\in U(d)`$，$`VU`$ 与 $`U`$ 同分布，所以
$`V^{\otimes n}\Phi(X)(V^{\otimes n})^*=\mathbb E[(VU)^{\otimes n}X((VU)^{\otimes n})^*]=\Phi(X)`$：值域含于交换子中。若 $`Y`$ 属于
交换子，则 $`U^{\otimes n}Y(U^{\otimes n})^*=Y`$，$`\Phi(Y)=Y`$。所以 $`\Phi`$ 是幂等的，其值域就是交换子。由
引理 5.1，$`U^*=U^{-1}`$ 与 $`U`$ 同分布，因此
$`\langle Y,\Phi(X)\rangle=\mathbb E\mathrm{Tr}\big(((U^*)^{\otimes n}YU^{\otimes n})^*X\big)=\langle\Phi(Y),X\rangle`$：$`\Phi`$ 是自伴的。
自伴的幂等映射就是到其值域上的正交投影
（$`\langle X-\Phi(X),\Phi(Z)\rangle=\langle\Phi(X-\Phi(X)),Z\rangle=0`$）。由定理 B.6(c)，交换子是 $`\mathrm{span}\,\pi(S_n)`$。$`\square`$

*计算：* 用 Haar 酉矩阵上的蒙特卡罗方法验证引理 B.11——[教程，B.5 节](../../tutorial/zh/B_schur_weyl.md#s-B-5)。

为了用张成族 $`\pi(\sigma)`$（当 $`n>d`$ 时它线性相关）来写出这个投影，我们使用厄米矩阵 $`G`$ 的
**伪逆**：$`G^+=f(G)`$，其中当 $`\lambda\neq0`$ 时 $`f(\lambda)=1/\lambda`$，且 $`f(0)=0`$（泛函演算，定理 1.5）。
它是厄米的；若 $`G`$ 可逆，则 $`G^+=G^{-1}`$；$`GG^+=\mathbf 1_{\neq0}(G)`$ 是到 $`(\ker G)^\perp`$ 上的投影；若
$`Gv=\lambda v`$ 且 $`\lambda\neq0`$，则 $`G^+v=v/\lambda`$；并且对酉矩阵 $`R`$ 有 $`(R^*GR)^+=R^*G^+R`$（若 $`G=Q\Lambda Q^*`$，
则 $`R^*GR=(R^*Q)\Lambda(R^*Q)^*`$）。

<a id="r-B-12"></a>

**引理 B.12.** 设 $`v_1,\dots,v_N`$ 是有限维内积空间 $`\mathcal K`$ 中的向量，$`G_{ab}=\langle v_a,v_b\rangle`$ 是它们的
Gram 矩阵，$`W=G^+`$。则到 $`\mathrm{span}\{v_a\}`$ 上的正交投影是
$`x\mapsto\sum_{a,b}v_a\,W_{ab}\,\langle v_b,x\rangle`$。

*证明。* 令 $`T:\mathbb C^N\to\mathcal K`$，$`Tc=\sum_ac_av_a`$。则 $`T^*x=(\langle v_b,x\rangle)_b`$，$`G=T^*T`$，所论映射是
$`Q=TG^+T^*`$。由于 $`\langle c,Gc\rangle=\|Tc\|^2`$，$`\ker G=\ker T`$，所以对到 $`\ker T`$ 上的投影 $`\Pi`$ 有
$`TG^+G=T(1-\Pi)=T`$。因此 $`QT=T`$：$`Q`$ 在 $`\mathrm{ran}\,T`$ 上是恒等映射。$`Q`$ 在 $`(\mathrm{ran}\,T)^\perp=\ker T^*`$ 上为零，
且 $`Q^*=Q`$。所以 $`Q`$ 是到 $`\mathrm{ran}\,T`$ 上的正交投影。$`\square`$

<a id="r-B-13"></a>

**定理 B.13（任意阶的 Weingarten 公式）。** 设 $`n,d\ge1`$，$`U`$ 在 $`U(d)`$ 上服从 Haar 分布，$`G`$ 是
$`n!\times n!`$ 矩阵 $`G_{\sigma\tau}=d^{\#\text{cycles}(\sigma^{-1}\tau)}`$（$`\sigma,\tau\in S_n`$），$`W=G^+`$。
(a) $`W_{\sigma\tau}=\mathrm{Wg}(\sigma^{-1}\tau)`$，其中 $`\mathrm{Wg}(\rho)=\mathrm{Wg}(\rho,d)=W_{e\rho}`$ 是 $`S_n`$ 上的函数，它只依赖于
$`\rho`$ 的循环型；特别地 $`\mathrm{Wg}(\rho^{-1})=\mathrm{Wg}(\rho)`$。
(b) $`G`$ 可逆当且仅当 $`n\le d`$，此时 $`W=G^{-1}`$。
(c) 对 $`X\in\mathrm{End}\,\mathcal H_n`$，
$`\mathbb E[U^{\otimes n}X(U^{\otimes n})^*]=\sum_{\sigma,\tau}\mathrm{Wg}(\sigma^{-1}\tau)\mathrm{Tr}(\pi(\tau)^*X)\,\pi(\sigma)`$。
(d) 对所有 $`\mathbf i,\mathbf j,\mathbf i',\mathbf j'\in[d]^n`$，
```math
\mathbb E\big[U_{i_1j_1}\cdots U_{i_nj_n}\overline{U_{i'_1j'_1}\cdots U_{i'_nj'_n}}\big]
=\sum_{\sigma,\tau\in S_n}\prod_{m=1}^n\delta_{i_mi'_{\sigma(m)}}\delta_{j_mj'_{\tau(m)}}\;\mathrm{Wg}(\sigma\tau^{-1}),
```
即求和遍及满足 $`\mathbf i=\mathbf i'\circ\sigma`$ 且 $`\mathbf j=\mathbf j'\circ\tau`$ 的对。
(e) $`U`$ 与 $`\bar U`$ 的因子个数 $`p\neq q`$ 不同的矩为零。

这就是定理 5.12（取 $`k=d`$），现在它得到了证明。

*证明。* (a) 记 $`g(\rho)=d^{\#\text{cycles}(\rho)}`$，于是 $`G_{\sigma\tau}=g(\sigma^{-1}\tau)`$，并令 $`\Lambda_\kappa,R_\kappa`$ 为
$`\mathbb C^{S_n}`$ 上满足 $`\Lambda_\kappa e_\sigma=e_{\kappa\sigma}`$，$`R_\kappa e_\sigma=e_{\sigma\kappa}`$ 的置换矩阵。则
$`(\Lambda_\kappa^*G\Lambda_\kappa)_{\sigma\tau}=g(\sigma^{-1}\kappa^{-1}\kappa\tau)=G_{\sigma\tau}`$，
$`(R_\kappa^*GR_\kappa)_{\sigma\tau}=g(\kappa^{-1}\sigma^{-1}\tau\kappa)=G_{\sigma\tau}`$（引理 B.1(c)）。由伪逆的性质，$`W`$ 具有
同样的两种不变性：$`W_{\kappa\sigma,\kappa\tau}=W_{\sigma\tau}`$，$`W_{\sigma\kappa,\tau\kappa}=W_{\sigma\tau}`$。在第一种中取
$`\kappa=\sigma^{-1}`$ 得 $`W_{\sigma\tau}=\mathrm{Wg}(\sigma^{-1}\tau)`$；第二种给出 $`\mathrm{Wg}(\kappa^{-1}\rho\kappa)=\mathrm{Wg}(\rho)`$，所以
$`\mathrm{Wg}`$ 在共轭类上为常数，而共轭类就是循环型（引理 B.1(c)）；$`\rho^{-1}`$ 与 $`\rho`$ 有相同的循环型。
(b) $`G`$ 是诸 $`\pi(\sigma)`$ 的 Gram 矩阵（引理 B.1(c)），取 $`T`$ 如引理 B.12 的证明中，则 $`G=T^*T`$ 可逆当且仅当
$`T`$ 是单射，即当且仅当诸 $`\pi(\sigma)`$ 线性无关：由引理 B.1(d)，当且仅当 $`n\le d`$。
(c) 取 $`\mathcal K=\mathrm{End}\,\mathcal H_n`$ 和 $`v_\sigma=\pi(\sigma)`$，应用引理 B.11 和引理 B.12。
(d) Kronecker 积的矩阵元为 $`\langle e_{\mathbf i},U^{\otimes n}e_{\mathbf j}\rangle=\prod_mU_{i_mj_m}`$，所以对
$`X=|e_{\mathbf j}\rangle\langle e_{\mathbf j'}|`$，该矩等于 $`\langle e_{\mathbf i},\Phi(X)e_{\mathbf i'}\rangle`$。这里
$`\mathrm{Tr}(\pi(\tau)^*X)=\langle\pi(\tau)e_{\mathbf j'},e_{\mathbf j}\rangle=\delta_{\mathbf j,\mathbf j'\circ\tau^{-1}}`$，
$`\langle e_{\mathbf i},\pi(\sigma)e_{\mathbf i'}\rangle=\delta_{\mathbf i,\mathbf i'\circ\sigma^{-1}}`$。由 (c)，该矩等于
$`\sum_{\sigma,\tau}\delta_{\mathbf i,\mathbf i'\circ\sigma^{-1}}\delta_{\mathbf j,\mathbf j'\circ\tau^{-1}}\mathrm{Wg}(\sigma^{-1}\tau)`$；代换
$`\sigma\to\sigma^{-1}`$，$`\tau\to\tau^{-1}`$ 即得公式。
(e) 对实数 $`\theta`$，$`e^{i\theta}1\in U(d)`$，所以 $`e^{i\theta}U`$ 与 $`U`$ 同分布，于是该矩等于它自身乘以
$`e^{i(p-q)\theta}`$。取 $`\theta`$ 使 $`e^{i(p-q)\theta}\neq1`$ 即可。$`\square`$

对 $`n=2\le d`$，$`G=\begin{pmatrix}d^2&d\\d&d^2\end{pmatrix}`$，由 $`G^{-1}`$ 得 $`\mathrm{Wg}(e)=\frac1{d^2-1}`$，
$`\mathrm{Wg}(\tau)=-\frac1{d(d^2-1)}`$：这就是定理 5.10。对 $`d=1`$，所有 $`G_{\sigma\tau}=1`$，$`W=\frac1{(n!)^2}G`$，而相位
$`U=U_{11}`$ 的矩为 $`\sum_{\sigma,\tau}\frac1{(n!)^2}=1=|U_{11}|^{2n}`$，理应如此。

*计算：* 定理 B.13，其中 $`\mathrm{Wg}`$ 由伪逆计算，包括 $`d=2`$ 时 $`G`$ 奇异的 $`3`$ 阶矩——[教程，B.5 节](../../tutorial/zh/B_schur_weyl.md#s-B-5)。

<a id="r-B-14"></a>

**推论 B.14.** 对所有 $`n,d\ge1`$，$`\mathbb E|U_{11}|^{2n}=\dfrac{n!}{d(d+1)\cdots(d+n-1)}`$。

*证明。* 在定理 B.13(d) 中令所有指标都为 $`1`$，该矩为 $`\sum_{\sigma,\tau}\mathrm{Wg}(\sigma\tau^{-1})=\sum_{\sigma,\tau}W_{\sigma^{-1},\tau^{-1}}=\mathbf 1^TW\mathbf 1`$（(a) 部分）。$`G`$ 的每个
行和都是 $`c_n(d)=\sum_{\rho\in S_n}d^{\#\text{cycles}(\rho)}\ge d^n>0`$，所以 $`G\mathbf 1=c_n(d)\mathbf 1`$，$`W\mathbf 1=\mathbf 1/c_n(d)`$，
该矩等于 $`n!/c_n(d)`$。最后 $`c_n(x)=x(x+1)\cdots(x+n-1)`$（第 5 章习题 2 的解答）。$`\square`$

这把第 5 章习题 2 推广到 $`G`$ 奇异的情形 $`n>d`$。

*计算：* $`n=3`$、$`d=2`$ 时的推论 B.14——[教程，B.5 节](../../tutorial/zh/B_schur_weyl.md#s-B-5)。

<a id="s-B-6"></a>

## B.6 张量算子与 $`SU(2)`$ 的 Wigner–Eckart 定理

9.3 节由 Wigner–Eckart 定理构造自旋关联函数。我们对 4.5 节的自旋表示证明它。自旋 $`k`$ 是定理 4.11 中
$`\mathbb C^{2k+1}`$ 上的表示 $`S^{(k)}_a`$，其基为 $`|k,q\rangle`$，$`q=k,\dots,-k`$；回忆
$`S^{(k)}_\pm|k,q\rangle=\sqrt{k(k+1)-q(q\pm1)}\,|k,q\pm1\rangle`$（对 $`S_-=S_+^*`$，这由 $`S_+`$ 的公式推出）。

**定义。** 设 $`(S_a)`$ 是 $`\mathbb C^N`$ 上的自旋表示，$`k\in\{0,\frac12,1,\dots\}`$。**$`k`$ 阶张量算子**是一族
$`N\times N`$ 矩阵 $`(T_q)_{q=-k,\dots,k}`$，满足
```math
[S_z,T_q]=q\,T_q,\qquad[S_\pm,T_q]=\sqrt{k(k+1)-q(q\pm1)}\;T_{q\pm1}
```
（当 $`q=\pm k`$ 时系数为零，此时 $`T_{\pm(k+1)}`$ 无定义）。**矢量算子**是满足 $`[S_a,T_b]=i\sum_c\varepsilon_{abc}T_c`$ 的三元组
$`(T_x,T_y,T_z)`$，即 $`[S_x,T_y]=iT_z`$，$`[S_x,T_z]=-iT_y`$ 以及轮换，并且
$`[S_a,T_a]=0`$。

<a id="r-B-15"></a>

**引理 B.15.** (a) 若 $`(T_a)`$ 是矢量算子，则它的**球分量** $`T_{\pm1}=\mp\frac1{\sqrt2}(T_x\pm iT_y)`$，
$`T_0=T_z`$ 构成一个 $`1`$ 阶张量算子，并且对两个矢量算子有
$`\sum_aT_aT'_a=\sum_q(-1)^qT_qT'_{-q}`$。
(b) 对 $`L`$ 个自旋 $`\frac12`$ 及总自旋 $`S_a`$，每个格点的自旋 $`(s^{(i)}_a)_a`$ 以及 $`(S_a)_a`$ 本身都是矢量算子。
(c) $`(T_q)`$ 是 $`k`$ 阶张量算子，当且仅当由 $`\Psi(v\otimes|k,q\rangle)=T_qv`$ 定义的映射
$`\Psi:\mathbb C^N\otimes\mathbb C^{2k+1}\to\mathbb C^N`$ 对 $`a=x,y,z`$ 满足 $`\Psi\,(S_a\otimes1+1\otimes S^{(k)}_a)=S_a\Psi`$。

*证明。* (a) 由这些关系得 $`[S_z,T_x\pm iT_y]=\pm(T_x\pm iT_y)`$，$`[S_\pm,T_z]=\mp(T_x\pm iT_y)`$，
$`[S_\pm,T_x\mp iT_y]=\pm2T_z`$，$`[S_\pm,T_x\pm iT_y]=0`$（例如 $`[S_x+iS_y,T_x-iT_y]=-i[S_x,T_y]+i[S_y,T_x]=T_z+T_z`$）。用
球分量表示：$`[S_z,T_q]=qT_q`$，$`[S_\pm,T_0]=\sqrt2\,T_{\pm1}`$，$`[S_\pm,T_{\mp1}]=\sqrt2\,T_0`$，$`[S_\pm,T_{\pm1}]=0`$，
这正是 $`k=1`$ 时的定义关系（$`\sqrt{2-q(q\pm1)}`$ 在 $`q=0`$ 和 $`q=\mp1`$ 时为 $`\sqrt2`$，在 $`q=\pm1`$ 时为 $`0`$）。
对乘积，$`T_{+1}T'_{-1}+T_{-1}T'_{+1}=-\frac12\big((T_x+iT_y)(T'_x-iT'_y)+(T_x-iT_y)(T'_x+iT'_y)\big)=-(T_xT'_x+T_yT'_y)`$。
(b) $`[S_a,s^{(i)}_b]=[s^{(i)}_a,s^{(i)}_b]`$（其他格点可交换），而自旋 $`\frac12`$ 的矩阵满足自旋关系；对 $`S_a`$ 则对 $`i`$
求和。(c) $`\Psi(S_a\otimes1+1\otimes S^{(k)}_a)(v\otimes|k,q\rangle)=T_qS_av+\sum_{q'}(S^{(k)}_a)_{q'q}T_{q'}v`$，所以关于 $`S_a`$ 的
缠结关系就是对所有 $`q`$ 有 $`[S_a,T_q]=\sum_{q'}(S^{(k)}_a)_{q'q}T_{q'}`$。对 $`S_z,S_\pm`$，这些就是定义关系，而关于
$`S_x,S_y,S_z`$ 的关系与关于 $`S_z,S_\pm`$ 的关系可以互相线性组合得到。$`\square`$

*计算：* 对六个自旋 $`\frac12`$ 的格点自旋验证引理 B.15(a) 与 (b)——[教程，B.6 节](../../tutorial/zh/B_schur_weyl.md#s-B-6)。

<a id="r-B-16"></a>

**定理 B.16（Wigner–Eckart）。** 设 $`(S_a)`$ 是 $`\mathbb C^N`$ 上的自旋表示，$`V`$ 是满足
$`V^*S_aV=\bigoplus_SS^{(S)}_a\otimes1_{m_S}`$ 的酉矩阵（推论 4.12），$`|S,M;b\rangle=V(|S,M\rangle\otimes f_b)`$ 是相应的规范正交
基（$`f_b`$ 是 $`\mathbb C^{m_S}`$ 的标准基），$`(T_q)`$ 是 $`k`$ 阶张量算子。设 $`S,S'`$ 出现（$`m_S,m_{S'}\ge1`$）。
(a) 若 $`S'\notin\{|S-k|,|S-k|+1,\dots,S+k\}`$，则对所有 $`M',M,q,a,b`$ 有 $`\langle S',M';a|T_q|S,M;b\rangle=0`$。
(b) 否则，从自旋 $`S\otimes`$ 自旋 $`k`$ 到自旋 $`S'`$ 的缠结算子构成一维空间；对每个这样的非零
$`C=C^{S'}_{S,k}:\mathbb C^{2S+1}\otimes\mathbb C^{2k+1}\to\mathbb C^{2S'+1}`$，存在唯一的 $`m_{S'}\times m_S`$ 矩阵
$`R=\langle S'\|T\|S\rangle`$，即**约化矩阵元**，使得
```math
\langle S',M';a|\,T_q\,|S,M;b\rangle=\langle S',M'|\,C\,(|S,M\rangle\otimes|k,q\rangle)\rangle\;R_{ab}\qquad\text{for all }M',M,q,a,b .
```
(c) 除非 $`M'=M+q`$，矩阵元为零。在情形 (b) 中，存在 $`q`$ 使 $`\langle S',S'|C(|S,S'-q\rangle\otimes|k,q\rangle)\rangle\neq0`$，
所以 $`R`$，以及随之 $`T`$ 在多重态 $`S`$ 与 $`S'`$ 之间的整个块，由 $`S'`$ 的最高权向量与 $`S`$ 的每个拷贝中的一个向量之间的
矩阵元 $`\langle S',S';a|T_q|S,S'-q;b\rangle`$ 决定。

*证明。* 固定 $`a\le m_{S'}`$，$`b\le m_S`$，令 $`J:\mathbb C^{2S+1}\to\mathbb C^N`$，$`Jx=V(x\otimes f_b)`$（在直和项 $`S`$ 中），并令
$`K:\mathbb C^N\to\mathbb C^{2S'+1}`$ 为 $`x\mapsto V(x\otimes f_a)`$（在直和项 $`S'`$ 中）的伴随。用 $`\alpha\in\{x,y,z\}`$ 表示
分量指标。由 $`V^*S_\alpha V`$ 的形式，$`S_\alpha J=JS^{(S)}_\alpha`$，并且对 $`S'`$ 的同一关系取伴随
（自旋矩阵是厄米的）得 $`KS_\alpha=S^{(S')}_\alpha K`$。令 $`t=K\Psi(J\otimes1)`$，其中 $`\Psi`$ 取自引理 B.15(c)；它的矩阵
元为 $`\langle S',M'|t(|S,M\rangle\otimes|k,q\rangle)\rangle=\langle S',M';a|T_q|S,M;b\rangle`$，并且
```math
t\,(S^{(S)}_\alpha\otimes1+1\otimes S^{(k)}_\alpha)=K\Psi(S_\alpha\otimes1+1\otimes S^{(k)}_\alpha)(J\otimes1)=KS_\alpha\Psi(J\otimes1)=S^{(S')}_\alpha\,t,
```
所以 $`t`$ 是从 $`\rho=`$ 自旋 $`S\otimes`$ 自旋 $`k`$（自旋表示
$`S^{(S)}_\alpha\otimes1+1\otimes S^{(k)}_\alpha`$，4.5 节）到 $`\sigma=`$ 自旋 $`S'`$ 的缠结算子。取伴随是
$`\mathrm{Hom}(\rho,\sigma)`$ 与 $`\mathrm{Hom}(\sigma,\rho)`$ 之间的线性双射（$`t\rho(\alpha)=\sigma(\alpha)t`$ 的伴随是
$`\rho(\alpha)t^*=t^*\sigma(\alpha)`$，因为所涉及的矩阵都是厄米的）。由推论 4.12 和第 4 章习题 2（其解答仅由推论 4.12
证明），$`\rho`$ 酉等价于 $`\bigoplus_j(\text{spin }j)`$，其中每个
$`j\in\{|S-k|,|S-k|+1,\dots,S+k\}`$ 各出现一次。这些自旋表示是不可约的且两两不等价（定理 4.11），所以这是定理 4.3
意义下的分解，而该定理的最后一个结论给出：若 $`S'`$ 是这些 $`j`$ 之一，则 $`\dim\mathrm{Hom}(\sigma,\rho)=1`$，
否则 $`\mathrm{Hom}(\sigma,\rho)=0`$（由定理 4.11，此时自旋 $`S'`$ 与它们都不等价）。
这就给出 (a)，并且在情形 (b) 中对唯一的数 $`R_{ab}`$ 有 $`t=R_{ab}C`$。
(c) 关于 $`S_z`$ 的缠结关系给出 $`M'\langle S',M'|t(|S,M\rangle\otimes|k,q\rangle)\rangle=(M+q)\langle S',M'|t(|S,M\rangle\otimes|k,q\rangle)\rangle`$。
$`C\neq0`$ 的值域是不可约的自旋 $`S'`$ 的一个非零不变子空间，因而是全空间，所以对某个基向量有
$`\langle S',S'|C(|S,M\rangle\otimes|k,q\rangle)\rangle\neq0`$，再对缠结算子 $`C`$ 应用同一关系得 $`M=S'-q`$。
除以这个数即可。$`\square`$

*计算：* 以不依赖基的形式，在 `decompose.from_generators` 计算出的分解中验证定理 B.16(a) 与 (b)——[教程，B.6 节](../../tutorial/zh/B_schur_weyl.md#s-B-6)。

由引理 B.15(a)，诸如
$`\mathbf s^{(i)}\cdot\mathbf s^{(j)}=\sum_q(-1)^q s^{(i)}_qs^{(j)}_{-q}`$ 这样的标量积在多重态 $`S`$ 与 $`S'`$ 之间的块是
$`\sum_{S''}\sum_q(-1)^q(s^{(i)}_q\text{ 的 }S''\to S'\text{ 块})(s^{(j)}_{-q}\text{ 的 }S\to S''\text{ 块})`$，而每个因子都由
定理 B.16 通过一个约化矩阵元和一个 Clebsch–Gordan 映射给出。

*库中：* 这就是 `su2.spin_correlation` 背后的原理，第 9 章把它与 `'direct'` 方法进行比较检验——[教程，9.3 节](../../tutorial/zh/09_spin_chains.md#s-9-3)。

<a id="s-B-7"></a>

## B.7 本课程中使用的 Haar 测度

本课程使用以下群的 Haar 测度：$`U(k)`$（定理 5.3），$`U(1)`$（$`k=1`$ 的情形；第 5 章习题 3），
$`U(\mathcal M)\cong\prod_cU(k_c)`$（推论 5.4），$`SU(k)`$（5.1 节，推论 5.4 之后），正交群
$`O(k)`$（5.1 节，第 11 章习题 2，以及实代数），以及引理 C.2 中的对角群 $`T`$。它们的唯一性都由引理 5.1 给出，因为
每个都是某个 $`U(n)`$ 的闭子群；存在性在每种情形下都由显式构造给出（对 $`T`$ 见引理 C.2），而 $`O(k)`$ 的构造如下。
本课程中任何地方都没有用到紧群 Haar 测度的一般存在性定理。

<a id="r-B-17"></a>

**定理 B.17（正交群的 Haar 测度）。** 设 $`Z`$ 是 $`k\times k`$ 实矩阵，其矩阵元相互独立且服从
$`N(0,1)`$。则 $`Z`$ 以概率一可逆，$`Z=OR`$，其中 $`O`$ 是正交矩阵，$`R`$ 是对角元为正的实上三角矩阵，并且 $`O`$ 在
$`O(k)`$ 上服从 Haar 分布。若 $`Z=O'R'`$ 是任意实 QR 分解（$`O'`$ 正交，$`R'`$ 为实上三角），则
$`O'\mathrm{diag}(\mathrm{sgn}\,R'_{ii})=O`$。对服从 Haar 分布的 $`O`$ 和 $`X\in M_k(\mathbb C)`$，
$`\mathbb E[O]=0`$，$`\mathbb E[OXO^T]=\frac{\mathrm{Tr}\,X}k1`$；特别地
$`\mathbb E[O_{ij}O_{lm}]=\delta_{il}\delta_{jm}/k`$。

*证明。* $`\det Z`$ 是 $`k^2`$ 个实矩阵元的多项式，且不恒为零（$`\det1=1`$），而 $`Z`$ 有密度，所以
$`\det Z\neq0`$ 几乎必然成立（推论 A.7(b)，与定理 5.3 的证明中一样）。对实的列作 Gram–Schmidt 正交化得到实的 $`q_j`$ 和实的
系数，即实分解 $`Z=OR`$，$`R_{jj}>0`$；由引理 5.2 中的唯一性，它就是引理 5.2 的分解，所以 $`O`$
连续依赖于 $`Z`$。$`Z`$ 的密度 $`(2\pi)^{-k^2/2}e^{-\|Z\|_2^2/2}`$ 在 $`Z\mapsto VZ`$（$`V\in O(k)`$）下不变，这是
$`\mathbb R^{k^2}`$ 上保持 $`\|Z\|_2^2=\mathrm{Tr}(Z^TZ)`$ 的线性映射，因而是正交的，其 Jacobi 行列式为 $`1`$；
所以 $`VZ`$ 与 $`Z`$ 同分布（命题 A.14(a)）。由于 $`VZ=(VO)R`$，唯一性给出 $`O(VZ)=VO(Z)`$，所以 $`O`$ 的分布在紧群
$`O(k)`$ 上是左不变的，从而可以应用引理 5.1。符号修正就是定理 5.3 证明的最后部分，其中
$`\mathrm{sgn}\,R'_{ii}=R'_{ii}/|R'_{ii}|`$（由于 $`Z`$ 可逆，$`R'_{ii}\neq0`$）。由于 $`-1\in O(k)`$，$`-O`$ 与 $`O`$ 同分布，所以
$`\mathbb E[O]=-\mathbb E[O]=0`$（矩阵元的绝对值以 $`1`$ 为界）。最后，由于 $`VO`$ 与 $`O`$ 同分布，$`T=\mathbb E[OXO^T]`$ 对
$`V\in O(k)`$ 满足 $`VTV^T=T`$。取 $`V`$ 为第 $`i`$ 个位置为 $`-1`$、其余为 $`1`$ 的对角矩阵，得对 $`j\neq i`$ 有 $`T_{ij}=-T_{ij}`$；
取置换矩阵，得所有 $`T_{ii}`$ 相等（与引理 5.6 中一样）。所以 $`T=\lambda1`$，且 $`k\lambda=\mathrm{Tr}\,T=\mathrm{Tr}\,X`$。对
$`X=|e_j\rangle\langle e_m|`$，有 $`(OXO^T)_{il}=O_{ij}O_{lm}`$。$`\square`$

*计算：* 用蒙特卡罗方法验证定理 B.17 的二阶矩，以及为何需要符号修正——[教程，B.7 节](../../tutorial/zh/B_schur_weyl.md#s-B-7)。

**返回：** [第 4 章](04_symmetries.md)（定理 4.14），[第 5 章](05_random_matrices.md)（定理 5.12），
[第 9 章](09_spin_chains.md)（关联函数）。
