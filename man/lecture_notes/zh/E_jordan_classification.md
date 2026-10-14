[English](../E_jordan_classification.md) · [Español](../es/E_jordan_classification.md) · [Français](../fr/E_jordan_classification.md) · [Deutsch](../de/E_jordan_classification.md) · **中文** · [日本語](../ja/E_jordan_classification.md) · [Русский](../ru/E_jordan_classification.md)

# 附录 E. 欧几里得若尔当代数的分类

第 13 章使用了有限维形式实若尔当代数的 Jordan–von Neumann–Wigner 分类。本附录仅从若尔当恒等式和
有限维线性代数出发完整地证明它。所走的是经典路线（J. Faraut, A. Korányi, *Analysis on Symmetric Cones*,
Chapters II–V；K. McCrimmon, *A Taste of Jordan Algebras*）：幂结合性与一个谱引理给出幂等元以及一个结合的、
正的迹形式；关于若尔当标架的 Peirce 分解把代数分成若干块 $`J_{ij}`$，它们的乘积服从严格的规则；在秩 $`r\ge3`$
的单代数中，这些规则把 Peirce 空间 $`J_{13}`$ 变成一个合成代数，由 Hurwitz 定理，它是
$`\mathbb R,\mathbb C,\mathbb H`$ 或 $`\mathbb O`$，而整个代数就是其上的
$`\mathrm{Herm}_r`$。E.6 节还证明了关于实结合可除代数的 Frobenius 定理，
E.8 节证明了 Albert 代数 $`\mathrm{Herm}_3(\mathbb O)`$ 是一个不能嵌入
任何结合代数的若尔当代数。

**你将学到**为什么幂结合性、迹形式的正性以及单位元的存在仅由形式实性就能推出；Peirce 规则，以及它们
如何迫使单代数的坐标像数一样相乘；用 Cayley–Dickson 倍增证明 Hurwitz 定理；以及为什么八元数仅在秩 $`3`$
时给出一个新的单代数。

*计算：* 本附录的数值检验——[教程，附录 E](../../tutorial/zh/E_jordan_classification.md)。

**约定。** 所有代数都是带有双线性乘积的有限维实向量空间。**若尔当代数**
如 13.1 节所述：$`x\circ y=y\circ x`$ 且 $`(x^2\circ y)\circ x=x^2\circ(y\circ x)`$，$`x^2=x\circ x`$；不假定
有单位元。$`L(x)y=x\circ y`$，$`[A,B]=AB-BA`$，$`\mathrm{Tr}`$ 是线性映射的迹。**子代数**是
在 $`\circ`$ 下封闭的子空间，**理想**是满足 $`J\circ I\subseteq I`$ 的子空间 $`I`$；若
$`J\ne0`$ 且它仅有的理想是 $`0`$ 与 $`J`$，则称 $`J`$ 是**单的**。**幂等元**是满足 $`c^2=c`$ 的 $`c`$；若 $`c\circ c'=0`$，则称幂等元 $`c,c'`$
**正交**。若 $`\sum_ix_i^2=0`$ 蕴含所有 $`x_i=0`$，则称 $`J`$ 是**形式实的**；若
它有单位元 $`1`$ 以及满足 $`\langle x\circ y,z\rangle=\langle y,x\circ z\rangle`$ 的内积，则称它是**欧几里得的**。
多项式恒等式按如下方式**线性化**：若 $`p(t)=\sum_kp_kt^k`$ 具有向量系数，且对所有实数 $`t`$ 有 $`p(t)=0`$，
则所有 $`p_k=0`$（在 $`\deg p+1`$ 个不同的点上取值，并求 Vandermonde 矩阵的逆）。

<a id="s-E-1"></a>

## E.1 线性化的若尔当恒等式与幂结合性

<a id="r-E-1"></a>

**引理 E.1.** 在每个若尔当代数中，对所有 $`x,y,z,w`$：

1. $`[L(x),L(x^2)]=0`$ 且 $`[L(x^2),L(z)]+2[L(x\circ z),L(x)]=0`$；
2. $`L(x^2\circ y)=2L(x\circ y)L(x)+L(x^2)L(y)-2L(x)L(y)L(x)`$；
3. $`L((x\circ w)\circ y)=L(w\circ y)L(x)+L(x\circ y)L(w)+L(x\circ w)L(y)-L(w)L(y)L(x)-L(x)L(y)L(w)`$；
4. $`L((x\circ w)\circ y)-L(x\circ(w\circ y))=[L(w),[L(x),L(y)]]`$；特别地，
   $`\mathrm{Tr}\,L((x\circ w)\circ y)=\mathrm{Tr}\,L(x\circ(w\circ y))`$。

*证明。* (1) 若尔当恒等式说的是 $`L(x)L(x^2)y=L(x^2)L(x)y`$。把 $`x`$ 换成 $`x+tz`$；$`t`$ 的系数为
$`[L(z),L(x^2)]+2[L(x),L(x\circ z)]=0`$。(2) 在若尔当恒等式本身中把 $`x`$ 换成 $`x+tz`$；$`t`$ 的系数为
```math
2((x\circ z)\circ y)\circ x+(x^2\circ y)\circ z=2(x\circ z)\circ(y\circ x)+x^2\circ(y\circ z).
```
把它们看作 $`z`$ 的线性映射，四项依次是 $`2L(x)L(y)L(x)z`$、$`L(x^2\circ y)z`$、$`2L(x\circ y)L(x)z`$ 和
$`L(x^2)L(y)z`$。(3) 在 (2) 中把 $`x`$ 换成 $`x+sw`$，取 $`s`$ 的系数并除以 $`2`$。(4) 把 (3) 照原样写一次，
再交换 $`x`$ 与 $`y`$ 写一次。两个右端的前三项相同，并且
$`(y\circ w)\circ x=x\circ(w\circ y)`$，所以差为
$`-L(w)L(y)L(x)-L(x)L(y)L(w)+L(w)L(x)L(y)+L(y)L(x)L(w)=[L(w),[L(x),L(y)]]`$，这是一个交换子，其迹为 $`0`$。
$`\square`$

*计算：* $`\mathrm{Herm}_3(\mathbb C)`$ 中的 (4)——[教程，E.3 节](../../tutorial/zh/E_jordan_classification.md#s-E-3)。

<a id="r-E-2"></a>

**引理 E.2（幂结合性）。** 令 $`x^1=x`$，$`x^{k+1}=x\circ x^k`$。则对所有
$`p,q\ge1`$ 有 $`x^p\circ x^q=x^{p+q}`$，并且算子 $`L(x^k)`$ 两两交换。因此，对没有常数项的实多项式 $`f,g`$ 以及
$`f(x):=\sum_kf_kx^k`$，有 $`(fg)(x)=f(x)\circ g(x)`$，且 $`\mathbb R[x]=\{f(x)\}`$ 是一个结合的交换子代数。

*证明。* 由引理 E.1(1)，$`L(x)`$ 与 $`L(x^2)`$ 交换；设 $`\mathcal A`$ 是由它们生成的交换算子
代数。我们对 $`n\ge2`$ 归纳证明：只要 $`p+q\le n`$ 就有 $`x^p\circ x^q=x^{p+q}`$，并且对
$`k\le n`$ 有 $`L(x^k)\in\mathcal A`$。$`n=2`$ 时这就是定义。设它对 $`n`$ 成立，并设 $`p+q=n+1`$。若 $`q=1`$，由
交换性得 $`x^p\circ x=x^{n+1}`$。若 $`q\ge2`$，则
$`x^p\circ x^q=L(x^p)L(x)x^{q-1}=L(x)L(x^p)x^{q-1}=L(x)x^{n}=x^{n+1}`$。其次，对 $`y=x^{n-1}`$ 使用引理 E.1(2)，并利用
$`x^2\circ x^{n-1}=x^{n+1}`$，得到 $`L(x^{n+1})=2L(x^n)L(x)+L(x^2)L(x^{n-1})-2L(x)^2L(x^{n-1})\in\mathcal A`$。其余
由双线性推出。$`\square`$

*库中：* 算子的若尔当积与若尔当幂，`jordan.jordan_product`、`jordan.jordan_power`——[教程，13.1 节](../../tutorial/zh/13_jordan_algebras.md#s-13-1)。

<a id="s-E-2"></a>

## E.2 幂等元与 Peirce 分解

<a id="r-E-3"></a>

**引理 E.3（谱引理）。** 设 $`J`$ 是形式实的，$`x\in J`$。则存在 $`k\ge0`$、互不相同的非零实数
$`\lambda_1,\dots,\lambda_k`$ 以及两两正交的非零幂等元 $`c_1,\dots,c_k\in\mathbb R[x]`$，使得
```math
x=\sum_{i=1}^k\lambda_ic_i,\qquad x^n=\sum_{i=1}^k\lambda_i^nc_i\quad(n\ge1),\qquad\mathbb R[x]=\mathrm{span}\{c_1,\dots,c_k\}.
```
特别地，每个非零的形式实若尔当代数都含有一个非零幂等元。

*证明。* 只有一项的形式实性说的是：$`z^2=0`$ 蕴含 $`z=0`$。令 $`K=\{q\in\mathbb R[t]:(tq)(x)=0\}`$。由
引理 E.2，$`K`$ 是 $`\mathbb R[t]`$ 的理想：若 $`q\in K`$，$`s=s_0+s_1`$ 且 $`s_1(0)=0`$，则
$`(tqs)(x)=s_0(tq)(x)+(tq)(x)\circ s_1(x)=0`$。由于 $`x,x^2,\dots,x^{\dim J+1}`$ 线性相关，$`K\ne0`$，所以对某个首一多项式 $`\mu`$ 有
$`K=\mu\,\mathbb R[t]`$。(i) $`\mu(0)\neq0`$：若 $`\mu=th`$，则 $`w=(th)(x)`$ 满足
$`w^2=(t\cdot\mu h)(x)=0`$，因为 $`\mu h\in K`$；于是 $`w=0`$，$`h\in K`$，从而 $`\mu\mid h`$，按次数这不可能。(ii) $`\mu`$ 没有
重因子：若 $`\mu=g^2h`$ 且 $`\deg g\ge1`$，则 $`w=(tgh)(x)`$ 满足 $`w^2=(t\cdot\mu\cdot th)(x)=0`$，所以 $`gh\in K`$，
按次数这不可能。(iii) $`\mu`$ 没有非实根：若 $`\mu=((t-\alpha)^2+\beta^2)h`$ 且 $`\beta\ne0`$，则
$`a=(t(t-\alpha)h)(x)`$ 与 $`b=(\beta th)(x)`$ 满足 $`a^2+b^2=(t\cdot\mu\cdot th)(x)=0`$，所以 $`b=0`$ 且 $`h\in K`$，这不可能。
因此 $`\mu=\prod_{i=1}^k(t-\lambda_i)`$，其中 $`\lambda_i`$ 是互不相同的非零实数。令
```math
p_i(t)=\frac t{\lambda_i}\prod_{j\ne i}\frac{t-\lambda_j}{\lambda_i-\lambda_j},\qquad c_i=p_i(x).
```
多项式 $`p_ip_l-\delta_{il}p_i`$ 等于 $`t`$ 乘以一个在每个 $`\lambda_j`$ 处为零的多项式，因而等于 $`t`$ 乘以
$`K`$ 的一个元素；所以 $`c_i\circ c_l=\delta_{il}c_i`$。同样，$`t-\sum_i\lambda_ip_i=t\big(1-\sum_i\lambda_ip_i/t\big)`$，而
括号中的多项式在每个 $`\lambda_j`$ 处为零，所以 $`x=\sum_i\lambda_ic_i`$。若 $`c_i=0`$，则 $`p_i/t\in K`$，其次数为 $`k-1<\deg\mu`$；
所以 $`c_i\ne0`$。$`x^n`$ 的公式由正交性推出，并且它表明 $`\mathbb R[x]\subseteq\mathrm{span}\{c_i\}`$。
$`\square`$

*计算：* 在 $`\mathrm{Herm}_3(\mathbb C)`$ 中仅用若尔当幂求出幂等元 $`p_i(x)`$——[教程，E.3 节](../../tutorial/zh/E_jordan_classification.md#s-E-3)。

<a id="r-E-4"></a>

**命题 E.4（关于一个幂等元的 Peirce 分解）。** 设 $`c`$ 是若尔当代数 $`J`$ 的幂等元。

1. $`L(c)(2L(c)-1)(L(c)-1)=0`$。因此 $`J=J_1(c)\oplus J_{1/2}(c)\oplus J_0(c)`$，其中
   $`J_\alpha(c)=\{x:c\circ x=\alpha x\}`$，且 $`c\in J_1(c)`$。
2. 对 $`z\in J_0(c)\cup J_1(c)`$ 有 $`[L(c),L(z)]=0`$，对 $`y\in J_{1/2}(c)`$ 有 $`[L(c),L(y^2)]=0`$。
3. $`J_1\circ J_1\subseteq J_1`$，$`J_0\circ J_0\subseteq J_0`$，$`J_1\circ J_0=0`$，$`(J_1+J_0)\circ J_{1/2}\subseteq J_{1/2}`$
   且 $`J_{1/2}\circ J_{1/2}\subseteq J_1+J_0`$（均相对于 $`c`$）。
4. 若 $`c'`$ 是与 $`c`$ 正交的幂等元，则 $`c+c'`$ 是幂等元，$`[L(c),L(c')]=0`$，并且
   $`J_1(c')\subseteq J_0(c)`$。

*证明。* (1) 对 $`x=y=c`$ 使用引理 E.1(2)，得 $`L(c)=3L(c)^2-2L(c)^3`$。被一个只有
单根的多项式零化的线性映射可对角化，其特征值都在这些根之中。(2) 对 $`x=c`$ 使用引理 E.1(1)，得
$`[L(c),L(z)]=2[L(c),L(c\circ z)]`$，即对 $`z\in J_\alpha(c)`$ 有 $`(1-2\alpha)[L(c),L(z)]=0`$。对 $`x=y`$
与 $`z=c`$ 使用引理 E.1(1)，得 $`[L(y^2),L(c)]+2[L(y\circ c),L(y)]=[L(y^2),L(c)]+[L(y),L(y)]=[L(y^2),L(c)]=0`$。(3) 由 (2)，若 $`z\in J_0\cup J_1`$，则 $`L(z)`$
把 $`L(c)`$ 的每个特征空间映到自身。这给出第一、第二和第四条规则，以及
$`J_1\circ J_0\subseteq J_1\cap J_0=0`$。对 $`y\in J_{1/2}`$，$`L(y^2)`$ 与 $`L(c)`$ 交换，所以 $`c\circ y^2=y^2\circ c`$
与 $`c`$ 一样是 $`L(c)`$ 属于特征值 $`1`$ 的特征向量；写 $`y^2=a_1+a_{1/2}+a_0`$，则 $`c\circ y^2=a_1+\frac12a_{1/2}`$ 属于 $`J_1`$，
所以 $`a_{1/2}=0`$。极化（$`y\to y+y'`$）给出最后一条规则。(4) $`(c+c')^2=c+2c\circ c'+c'=c+c'`$，并且
$`c'\in J_0(c)`$，所以由 (2)，$`L(c)`$ 与 $`L(c')`$ 交换，且可同时对角化；特别地，$`J_1(c')`$ 由
公共特征向量张成。在一个公共特征向量上，
若 $`L(c)`$ 的特征值为 $`\alpha`$，$`L(c')`$ 的特征值为 $`1`$，则 $`L(c+c')`$ 的特征值为 $`\alpha+1\in\{0,\frac12,1\}`$，所以
$`\alpha=0`$。$`\square`$

*计算：* $`\mathrm{Herm}_3(\mathbb C)`$ 中的 (1) 以及 $`L(c)`$ 的特征值——[教程，E.3 节](../../tutorial/zh/E_jordan_classification.md#s-E-3)。

<a id="s-E-3"></a>

## E.3 迹形式；欧几里得 = 形式实

<a id="r-E-5"></a>

**定理 E.5.** 设 $`J`$ 是形式实的，$`\tau(x,y)=\mathrm{Tr}\,L(x\circ y)`$。

1. $`\tau`$ 是对称、双线性且**结合的**：$`\tau(x\circ y,z)=\tau(x,y\circ z)`$。
2. $`\tau`$ 是正定的；对每个非零幂等元 $`c`$，$`\tau(c,c)=\mathrm{Tr}\,L(c)\ge1`$。
3. $`J`$ 有单位元。
4. $`J`$ 是欧几里得的，内积为 $`\tau`$。反之，每个欧几里得若尔当代数都是形式实的。

*证明。* (1) 引理 E.1(4)。(2) 由命题 E.4(1)，$`\mathrm{Tr}\,L(c)=\dim J_1(c)+\frac12\dim J_{1/2}(c)\ge1`$，因为
$`0\ne c\in J_1(c)`$。对引理 E.3 中的 $`x=\sum_i\lambda_ic_i`$，$`x^2=\sum_i\lambda_i^2c_i`$，且
$`\tau(x,x)=\mathrm{Tr}\,L(x^2)=\sum_i\lambda_i^2\mathrm{Tr}\,L(c_i)`$，除非 $`k=0`$（即 $`x=0`$），它 $`>0`$。
(3) 设 $`J\neq0`$，$`e`$ 是使 $`\dim J_1(e)`$ 最大的非零幂等元（引理 E.3）。若 $`J_0(e)\ne0`$，它是一个
形式实子代数（命题 E.4(3)），从而含有非零幂等元 $`f`$；$`f`$ 与 $`e`$ 正交，所以
$`e+f`$ 是幂等元，$`J_1(e)\subseteq J_1(e+f)`$（因为 $`f\circ J_1(e)\subseteq J_0\circ J_1=0`$），且
$`f\in J_1(e+f)\setminus J_1(e)`$，与最大性矛盾。所以 $`J_0(e)=0`$。于是对 $`y\in J_{1/2}(e)`$ 有
$`y^2\in J_1(e)`$，而结合性给出
$`\frac12\tau(y,y)=\tau(e\circ y,y)=\tau(e,y^2)=\mathrm{Tr}\,L(e\circ y^2)=\mathrm{Tr}\,L(y^2)=\tau(y,y)`$，所以 $`y=0`$。
因此 $`J=J_1(e)`$，$`e`$ 是单位元。(4) 由 (1)–(3)，$`\tau(x\circ y,z)=\tau(y,x\circ z)`$。反之，若 $`J`$ 是欧几里得的且
$`\sum_ix_i^2=0`$，则 $`0=\sum_i\langle x_i\circ x_i,1\rangle=\sum_i\langle x_i,x_i\rangle`$。$`\square`$

*计算：* $`\mathrm{Herm}_3(\mathbb C)`$ 中的 (1) 与 (2)——[教程，E.3 节](../../tutorial/zh/E_jordan_classification.md#s-E-3)；对每种类型检验 $`\tau`$ 的正性——[教程，E.9 节](../../tutorial/zh/E_jordan_classification.md#s-E-9)。

所以“形式实”与“欧几里得”是同一类有限维代数。从现在起 $`J`$ 是形式实的，
$`1`$ 是它的单位元，$`\tau`$ 是它的迹形式。

<a id="s-E-4"></a>

## E.4 若尔当标架与谱定理

非零幂等元若不是两个正交非零幂等元之和，则称为**本原的**。**若尔当标架**
是两两正交的本原幂等元组成的族 $`(c_1,\dots,c_r)`$，满足 $`\sum_ic_i=1`$。

<a id="r-E-6"></a>

**引理 E.6.** (1) 若 $`c`$ 是本原的，则 $`J_1(c)=\mathbb Rc`$。(2) 两两正交的非零幂等元
线性无关；每个非零幂等元都是两两正交的本原幂等元之和；并且 $`J`$ 有一个
若尔当标架。

*证明。* (1) $`J_1(c)`$ 是以 $`c`$ 为单位元的形式实子代数（命题 E.4(3)）。对 $`x\in J_1(c)`$，由引理 E.3 写
$`x=\sum_i\lambda_ic_i`$，其中 $`c_i\in\mathbb R[x]\subseteq J_1(c)`$。则 $`c-c_i`$ 是与 $`c_i`$ 正交的
幂等元（$`(c-c_i)^2=c-2c_i+c_i`$ 且 $`c_i\circ(c-c_i)=0`$），所以由本原性 $`c-c_i=0`$。不同的 $`c_i`$ 互相正交，
所以至多有一个，从而 $`x\in\mathbb Rc`$。(2) 若 $`\sum_ia_if_i=0`$，与 $`f_j`$ 相乘得 $`a_jf_j=0`$。所以
幂等元 $`e`$ 是尽可能多个两两正交的非零幂等元之和。若其中之一 $`f`$
等于 $`g+g'`$，其中 $`g,g'`$ 是正交的非零幂等元，则对其他每个成员 $`f'`$ 有 $`g=f\circ g\in J_1(f)\subseteq J_0(f')`$
（命题 E.4(4)），$`g'`$ 也一样，所以可以用 $`g,g'`$ 代替 $`f`$。取 $`e=1`$。$`\square`$

<a id="r-E-7"></a>

**定理 E.7（谱定理）。** 对每个 $`x\in J`$，存在若尔当标架 $`(c_1,\dots,c_r)`$ 和实数
$`\lambda_1,\dots,\lambda_r`$（不必互不相同），使得 $`x=\sum_i\lambda_ic_i`$。

*证明。* 由引理 E.3 写 $`x=\sum_{i=1}^k\lambda_ie_i`$，再添上 $`e_0=1-\sum_ie_i`$，$`\lambda_0=0`$；$`e_0`$ 是
与所有 $`e_i`$ 正交的幂等元。把每个非零的 $`e_i`$ 拆成两两正交的本原幂等元
（引理 E.6(2)）。$`e_i`$ 的一块 $`f`$ 落在 $`J_1(e_i)\subseteq J_0(e_l)`$（$`l\ne i`$）中，而 $`e_l`$ 的一块落在
$`J_1(e_l)`$ 中；由于 $`J_0(e_l)\circ J_1(e_l)=0`$，所有块两两正交。它们的和为 $`\sum_ie_i=1`$。$`\square`$

$`J`$ 的所有若尔当标架都有相同个数的元素，这一点在推论 E.28 中证明。

<a id="r-E-8"></a>

**定理 E.8（关于标架的 Peirce 分解）。** 设 $`(c_1,\dots,c_r)`$ 是若尔当标架，$`t_i=\tau(c_i,c_i)`$，并且
对 $`i\ne j`$ 令 $`J_{ij}=J_{ji}=J_{1/2}(c_i)\cap J_{1/2}(c_j)`$。

1. 诸 $`L(c_k)`$ 两两交换，$`J=\bigoplus_i\mathbb Rc_i\oplus\bigoplus_{i<j}J_{ij}`$，并且对 $`x\in J_{ij}`$ 有
   $`c_k\circ x=\frac12(\delta_{ki}+\delta_{kj})x`$。
2. 对互不相同的 $`i,j,k`$：$`J_{ij}\circ J_{jk}\subseteq J_{ik}`$；若 $`\{i,j\}\cap\{k,l\}=\emptyset`$，则 $`J_{ij}\circ J_{kl}=0`$。
3. $`J_{ij}`$ 上存在正定二次型 $`N`$，使得 $`x^2=N(x)(c_i+c_j)`$；记其极化形式为
   $`\langle\cdot,\cdot\rangle`$，则对 $`x,y\in J_{ij}`$ 有 $`x\circ y=\langle x,y\rangle(c_i+c_j)`$ 与 $`\tau(x,y)=(t_i+t_j)\langle x,y\rangle`$。若 $`J_{ij}\ne0`$，则 $`t_i=t_j`$。

*证明。* (1) 由命题 E.4(4)，诸 $`L(c_k)`$ 两两交换，它们可对角化，且和为 $`L(1)=\mathrm{id}`$。一个公共
特征向量的特征值 $`\alpha_k\in\{0,\frac12,1\}`$ 满足 $`\sum_k\alpha_k=1`$：要么某一个 $`\alpha_i=1`$ 而其余
为零，此时该向量属于 $`J_1(c_i)=\mathbb Rc_i`$（引理 E.6(1)）；要么对某一对有 $`\alpha_i=\alpha_j=\frac12`$ 而
其余为零。所以 $`J_{ij}`$ 恰好是第二种公共特征空间，而 $`J`$ 是公共
特征空间的直和。(2) 设 $`x\in J_{ij}`$，$`y\in J_{jk}`$。相对于 $`c_i`$，$`x\in J_{1/2}`$ 且 $`y\in J_0`$，所以
$`x\circ y\in J_{1/2}(c_i)`$；相对于 $`c_k`$，同样有 $`x\circ y\in J_{1/2}(c_k)`$。因此 $`x\circ y\in J_{ik}`$。若改为
$`y\in J_{kl}`$，则对 $`p=i,j,k,l`$，$`x\circ y`$ 属于 $`J_{1/2}(c_p)`$；四个特征值都为
$`\frac12`$ 的公共特征向量不存在，所以 $`x\circ y=0`$。(3) 对 $`x,y\in J_{ij}`$，对 $`p=i,j`$ 有 $`x\circ y\in J_0(c_p)+J_1(c_p)`$，
对 $`k\notin\{i,j\}`$ 有 $`x\circ y\in J_0(c_k)`$（命题 E.4(3)）；由 (1)，$`x\circ y\in\mathbb Rc_i+\mathbb Rc_j`$。写
$`x^2=\alpha c_i+\beta c_j`$。取 $`y=c_i`$ 的若尔当恒等式给出 $`(x^2\circ c_i)\circ x=\frac12\alpha x`$ 以及
$`x^2\circ(c_i\circ x)=\frac12x^2\circ x=\frac14(\alpha+\beta)x`$，所以当 $`x\ne0`$ 时 $`\alpha=\beta`$；令 $`N(x)=\alpha`$。由于
$`\tau(c_i,1)=\mathrm{Tr}\,L(c_i)=t_i`$，$`\tau(x,x)=\tau(x^2,1)=N(x)(t_i+t_j)`$；所以 $`N`$ 是二次型，由定理 E.5(2)
它是正定的，极化给出 $`x\circ y`$ 与 $`\tau`$ 的公式。最后
$`N(x)t_i=\tau(x^2,c_i)=\tau(x,x\circ c_i)=\frac12\tau(x,x)`$，对 $`j`$ 也一样；所以只要存在某个 $`x\ne0`$，就有 $`t_i=t_j`$。$`\square`$

*计算：* 对每种类型，把 Peirce 空间作为诸 $`L(c_i)`$ 的公共特征空间求出——[教程，E.9 节](../../tutorial/zh/E_jordan_classification.md#s-E-9)。

<a id="s-E-5"></a>

## E.5 单代数与秩二

<a id="r-E-9"></a>

**定理 E.9.** $`J`$ 是单理想的直和 $`J^{(1)}\oplus\dots\oplus J^{(s)}`$，且当
$`p\ne q`$ 时 $`J^{(p)}\circ J^{(q)}=0`$。每个 $`J^{(p)}`$ 都是形式实若尔当代数。

*证明。* 若 $`I`$ 是理想，则它关于 $`\tau`$ 的正交补 $`I^\perp`$ 也是理想：
$`\tau(J\circ I^\perp,I)=\tau(I^\perp,J\circ I)=0`$。因此 $`J=I\oplus I^\perp`$ 且 $`I\circ I^\perp\subseteq I\cap I^\perp=0`$，所以
$`J`$ 是代数 $`I`$ 与 $`I^\perp`$ 的直和，而一个直和项的理想也是 $`J`$ 的理想。设 $`I`$ 是
维数最小的非零理想；则 $`I`$ 是单的，再对维数归纳，应用于
形式实代数 $`I^\perp`$，即得结论。$`\square`$

<a id="r-E-10"></a>

**引理 E.10.** 设 $`(c_i)`$ 是若尔当标架，$`i,j,k,l`$ 互不相同，$`u\in J_{ij}`$，$`w\in J_{jk}`$，$`z\in J_{kl}`$。则
(1) $`u\circ(u\circ w)=\frac14N(u)w`$，(2) $`(u\circ w)\circ z=u\circ(w\circ z)`$。

*证明。* (1) 把取 $`x=u`$ 的引理 E.1(1) 作用于 $`c_j`$。利用 $`L(u^2)=N(u)L(c_i+c_j)`$、$`c_j\circ w=\frac12w`$、
$`c_j\circ u=\frac12u`$ 以及 $`c_j\circ(u\circ w)=0`$（因为 $`u\circ w\in J_{ik}`$）：
$`[L(u^2),L(w)]c_j=\frac14N(u)w-\frac12N(u)w`$，$`2[L(u\circ w),L(u)]c_j=(u\circ w)\circ u`$；两者之和为零。(2) 把
引理 E.1(3)（其中 $`(x,w,y)`$ 换成 $`(u,z,w)`$）作用于 $`c_j`$。这里 $`u\circ z=0`$（定理 E.8(2)），$`c_j\circ u=\frac12u`$
且 $`c_j\circ z=0`$，所以 $`0=\frac12(z\circ w)\circ u-\frac12z\circ(w\circ u)`$。$`\square`$

<a id="r-E-11"></a>

**命题 E.11.** 设 $`(c_1,\dots,c_r)`$ 是 $`J`$ 的若尔当标架，若 $`J_{ij}\ne0`$，则用一条边连接 $`i\neq j`$。

1. $`J`$ 是单的当且仅当这个图是连通的。
2. 若 $`J`$ 是单的且 $`r\ge2`$，则所有 $`J_{ij}`$（$`i\ne j`$）有相同的维数 $`d\ge1`$，所有 $`t_i`$ 都等于
   $`t=1+\frac{(r-1)d}2`$，且 $`\dim J=r+d\,\frac{r(r-1)}2=rt`$。

*证明。* (1) 若图不连通，设 $`S`$ 是一个连通分支，$`I`$ 是由 $`i,j\in S`$ 的 $`c_i`$ 与 $`J_{ij}`$
张成的空间，$`I'`$ 是补集的相应空间。则 $`J=I\oplus I'`$（所有 $`i\in S`$、$`j\notin S`$ 的 $`J_{ij}`$ 都为零），
并且由定理 E.8，$`I\circ I\subseteq I`$，$`I\circ I'=0`$；所以 $`I`$ 是理想，且 $`\ne0,J`$。反之，设图
连通，$`I\ne0`$ 是理想。$`I`$ 在两两交换的可对角化算子 $`L(c_k)`$ 下不变，所以它是它与各公共特征空间
之交的和。若 $`0\ne x\in I\cap J_{ij}`$，则 $`c_i\circ x^2=N(x)c_i\in I`$；所以某个 $`c_i`$
属于 $`I`$。若 $`c_i\in I`$ 且 $`0\ne x\in J_{ij}`$，则 $`x=2c_i\circ x\in I`$，$`c_j=c_j\circ x^2/N(x)\in I`$。由
连通性，所有 $`c_j\in I`$，所以 $`1\in I`$，$`I=J`$。(2) 设 $`u\in J_{ij}`$，$`N(u)=1`$，$`k\notin\{i,j\}`$。由
定理 E.8(2)，$`2L(u)`$ 把 $`J_{jk}`$ 映到 $`J_{ik}`$，把 $`J_{ik}`$ 映到 $`J_{jk}`$，而由引理 E.10(1)（也在交换 $`i,j`$
之后使用），它的平方在这两个空间上都是恒等映射。所以只要 $`J_{ij}\ne0`$，就有 $`\dim J_{ik}=\dim J_{jk}`$。于是边 $`ij`$、$`jk`$
迫使存在边 $`ik`$；从而连通图是完全图，每次改变一个指标即可看出所有 $`\dim J_{ij}`$
相等。最后，$`L(c_i)`$ 在 $`\mathbb Rc_i`$ 上的特征值为 $`1`$，在 $`r-1`$ 个空间 $`J_{ij}`$ 上为 $`\frac12`$，在其他地方
为 $`0`$，所以 $`t_i=\mathrm{Tr}\,L(c_i)=1+\frac{(r-1)d}2`$。$`\square`$

*计算：* 对每种类型检验 (2) 的维数公式——[教程，E.9 节](../../tutorial/zh/E_jordan_classification.md#s-E-9)。

在本节余下部分中，$`J`$ 是单的，$`(c_1,\dots,c_r)`$ 是 $`r\ge2`$ 的若尔当标架，对 $`x\in J_{ij}`$、
$`y\in J_{jk}`$（$`i,j,k`$ 互不相同）记
```math
m(x,y)=2\,x\circ y\in J_{ik}.
```
由定理 E.8(3) 与命题 E.11(2)，在每个 $`J_{ij}`$ 上 $`\tau(x,y)=2t\langle x,y\rangle`$。

<a id="r-E-12"></a>

**引理 E.12.** 对互不相同的 $`i,j,k`$ 以及 $`x\in J_{ij}`$、$`y\in J_{jk}`$、$`z\in J_{ik}`$：

1. $`N(m(x,y))=N(x)N(y)`$ 且 $`m(x,m(x,y))=N(x)y`$；
2. $`\langle m(x,y),z\rangle=\langle y,m(x,z)\rangle=\langle x,m(y,z)\rangle`$；
3. 对 $`l\notin\{i,j,k\}`$ 与 $`v\in J_{kl}`$，$`m(m(x,y),v)=m(x,m(y,v))`$。

*证明。* (2) 这三个数是 $`\tau(x\circ y,z)/t`$、$`\tau(y,x\circ z)/t`$ 与 $`\tau(x,y\circ z)/t`$，由
结合性与交换性它们相等。(1) 第二个公式就是引理 E.10(1)。于是
$`2tN(m(x,y))=4\tau(x\circ y,x\circ y)=4\tau(y,x\circ(x\circ y))=N(x)\tau(y,y)=2tN(x)N(y)`$。(3) 引理 E.10(2)。$`\square`$

<a id="r-E-13"></a>

**定理 E.13（秩一与秩二）。** 设 $`J`$ 是单的，具有 $`r`$ 个元素的若尔当标架。若 $`r=1`$，则
$`J=\mathbb Rc_1\cong\mathbb R`$。若 $`r=2`$，则 $`J\cong V_{d+1}`$，其中 $`d=\dim J_{12}\ge1`$（自旋因子如 13.2 节所述）。

*证明。* $`r=1`$：定理 E.8(1)。$`r=2`$：$`J=\mathbb R1\oplus W`$，其中 $`W=\mathbb Rf\oplus J_{12}`$，$`f=c_1-c_2`$。由
定理 E.8，对 $`x,y\in J_{12}`$ 有 $`f\circ f=1`$、$`f\circ x=\frac12x-\frac12x=0`$ 以及 $`x\circ y=\langle x,y\rangle1`$。所以
对 $`w,w'\in W`$ 有 $`w\circ w'=B(w,w')1`$，其中 $`B(\alpha f+x,\beta f+y)=\alpha\beta+\langle x,y\rangle`$ 是一个内积。
若 $`U:W\to\mathbb R^{d+1}`$ 是等距映射，则 $`a1+w\mapsto(a,Uw)`$ 把
$`(a1+w)\circ(b1+w')=(ab+B(w,w'))1+aw'+bw`$ 映为 $`V_{d+1}`$ 中的乘积。$`\square`$

*库中：* 自旋因子 $`V_m`$ 即 `SpinFactor`——[教程，13.2 节](../../tutorial/zh/13_jordan_algebras.md#s-13-2)。

<a id="s-E-6"></a>

## E.6 合成代数：Hurwitz 定理与 Frobenius 定理

**合成代数**是具有单位元 $`1`$ 以及正定二次型 $`N`$ 的代数 $`A`$，满足
$`N(xy)=N(x)N(y)`$。我们把极化形式记为 $`\langle\cdot,\cdot\rangle`$（$`N(x)=\langle x,x\rangle`$），
并令 $`\mathrm{Re}\,x=\langle x,1\rangle`$，$`\bar x=2\mathrm{Re}(x)1-x`$。由于 $`N(1)=1`$（引理 E.14(1)，其
证明不用到共轭），这是关于直线 $`\mathbb R1`$ 的反射，所以 $`x\mapsto\bar x`$ 是等距的
对合，且 $`\mathrm{Re}\,\bar x=\mathrm{Re}\,x`$。合成代数的同构是指等距的代数
同构。**结合子**是 $`(x,y,z)=(xy)z-x(yz)`$。

<a id="r-E-14"></a>

**引理 E.14.** 在合成代数中，对所有 $`x,y,z,w`$：

1. $`N(1)=1`$，$`\langle xy,xz\rangle=N(x)\langle y,z\rangle`$，$`\langle yx,zx\rangle=\langle y,z\rangle N(x)`$，以及
   $`\langle xy,wz\rangle+\langle wy,xz\rangle=2\langle x,w\rangle\langle y,z\rangle`$；
2. $`\langle xy,z\rangle=\langle y,\bar xz\rangle=\langle x,z\bar y\rangle`$；
3. $`\bar x(xy)=N(x)y=(yx)\bar x`$；因此 $`\bar xx=x\bar x=N(x)1`$，$`x^2-2\mathrm{Re}(x)x+N(x)1=0`$，
   $`xy+yx=2\mathrm{Re}(x)y+2\mathrm{Re}(y)x-2\langle x,y\rangle1`$，
   $`\bar x(wy)+\bar w(xy)=2\langle x,w\rangle y`$ 以及 $`(yx)\bar w+(yw)\bar x=2\langle x,w\rangle y`$；
4. $`\overline{xy}=\bar y\bar x`$，$`\mathrm{Re}(xy)=\mathrm{Re}(yx)`$，以及
   $`\mathrm{Re}((xy)z)=\mathrm{Re}(x(yz))`$；
5. $`x(xy)=x^2y`$ 且 $`(yx)x=yx^2`$；结合子是交错的（交换两个变元时它改变符号）。

*证明。* (1) $`N(1)=N(1)^2>0`$。把 $`N(xy)=N(x)N(y)`$ 分别对 $`y`$、对 $`x`$ 极化，再把 $`\langle xy,xz\rangle=N(x)\langle y,z\rangle`$ 对 $`x`$ 极化。
(2) (1) 的最后一个恒等式取 $`w=1`$ 即为 $`\langle xy,z\rangle+\langle y,xz\rangle=2\mathrm{Re}(x)\langle y,z\rangle`$，
即 $`\langle xy,z\rangle=\langle y,\bar xz\rangle`$；把 $`\langle yx,zx\rangle=\langle y,z\rangle N(x)`$ 对 $`x`$ 极化，
并令两个变元之一等于 $`1`$，得 $`\langle yx,z\rangle=\langle y,z\bar x\rangle`$。(3) 由 (2) 与 (1)，
对所有 $`z`$ 有 $`\langle\bar x(xy),z\rangle=\langle xy,xz\rangle=N(x)\langle y,z\rangle`$，右侧同理。令 $`y=1`$
并代入 $`\bar x=2\mathrm{Re}(x)1-x`$；然后极化。(4) 由 (2) 以及共轭的等距性，
$`\langle\bar y\bar x,z\rangle=\langle\bar x,yz\rangle=\langle\bar x\bar z,y\rangle=\langle\bar z,xy\rangle=\langle z,\overline{xy}\rangle`$。
其次 $`\mathrm{Re}(xy)=\langle y,\bar x\rangle=\langle\bar y,x\rangle=\mathrm{Re}(yx)`$，并且
$`\mathrm{Re}((xy)z)=\langle xy,\bar z\rangle=\langle y,\bar x\bar z\rangle=\langle yz,\bar x\rangle=\mathrm{Re}(x(yz))`$。
(5) 由 (3)，$`x(xy)=(2\mathrm{Re}(x)-\bar x)(xy)=(2\mathrm{Re}(x)x-N(x))y=x^2y`$，右侧同理。所以
$`(x,x,y)=(y,x,x)=0`$；极化后，交换前两个或后两个变元时结合子改变符号，
而这些对换生成所有置换。$`\square`$

*计算：* 在随机八元数上检验 $`N(xy)=N(x)N(y)`$、(2)、(4) 与 (5)——[教程，E.6 节](../../tutorial/zh/E_jordan_classification.md#s-E-6)。

**Cayley–Dickson 构造。** 对具有单位元、对合 $`b\mapsto\bar b`$ 与二次型 $`N`$ 的代数 $`B`$，令
$`\mathrm{CD}(B)=B\times B`$，其中
```math
(a,b)(c,d)=(ac-\bar db,\;da+b\bar c),\qquad N(a,b)=N(a)+N(b),\qquad\overline{(a,b)}=(\bar a,-b).
```
从 $`\mathbb R`$（$`\bar x=x`$，$`N(x)=x^2`$）出发，我们定义 $`\mathbb C=\mathrm{CD}(\mathbb R)`$、
$`\mathbb H=\mathrm{CD}(\mathbb C)`$ 与 $`\mathbb O=\mathrm{CD}(\mathbb H)`$，维数分别为 $`1,2,4,8`$，每次倍增
都使用刚刚定义的对合与二次型（由下面的引理 E.15(1)，它们就是该
合成代数的共轭与范数）。在
$`\mathbb C`$ 中，$`(a,b)(c,d)=(ac-bd,ad+bc)`$：这就是复数 $`a+bi`$。$`B`$ 是子代数 $`B\times0`$。

*计算：* $`\mathbb R^{2^k}`$ 上的 Cayley–Dickson 乘积（函数 `cd`）——[教程，E.6 节](../../tutorial/zh/E_jordan_classification.md#s-E-6)。

<a id="r-E-15"></a>

**引理 E.15.** (1) 若 $`B`$ 是结合的合成代数，则 $`\mathrm{CD}(B)`$ 是合成代数，
单位元为 $`(1,0)`$，共轭为 $`\overline{(a,b)}=(\bar a,-b)`$。(2) $`\mathbb C`$ 是交换且结合的。
$`M(a,b)=\begin{pmatrix}a&-b\\\bar b&\bar a\end{pmatrix}`$ 是一个单的代数同态 $`\mathbb H\to M_2(\mathbb C)`$，
满足 $`M(\bar q)=M(q)^*`$；所以 $`\mathbb H`$ 是结合的；它不是交换的。$`\mathbb O`$ 不是结合的。(3)
$`\mathbb R,\mathbb C,\mathbb H,\mathbb O`$ 是合成代数，而 $`\mathbb R,\mathbb C,\mathbb H`$ 是可除代数：
每个 $`x\ne0`$ 都有双边逆 $`\bar x/N(x)`$，所以 $`xy=0`$ 蕴含 $`x=0`$ 或 $`y=0`$。

*证明。* (1) $`N((a,b)(c,d))=N(a)N(c)+N(d)N(b)+N(d)N(a)+N(b)N(c)-2\langle ac,\bar db\rangle+2\langle da,b\bar c\rangle`$，并且由
引理 E.14(2)，$`\langle ac,\bar db\rangle=\langle d(ac),b\rangle`$，$`\langle da,b\bar c\rangle=\langle(da)c,b\rangle`$，当 $`B`$ 结合时
两者相等。$`(1,0)`$ 是单位元，且 $`2\langle(a,b),(1,0)\rangle(1,0)-(a,b)=(\bar a,-b)`$。(2) 展开乘积
并利用 $`\mathbb C`$ 的交换性，$`M(a,b)M(c,d)=M(ac-\bar db,da+b\bar c)`$；$`M`$ 是线性单射，且
$`M(\bar a,-b)=M(a,b)^*`$。在 $`\mathbb H`$ 中，$`(i,0)(0,1)=(0,i)`$，但 $`(0,1)(i,0)=(0,-i)`$。在 $`\mathbb O=\mathbb H\times\mathbb H`$ 中令
$`p=((i,0),0)`$，$`q=((0,1),0)`$，$`s=(0,1)`$；则 $`(pq)s=(0,(0,i))`$，而 $`p(qs)=(0,(0,-i))`$。(3) 由 (1)、(2) 与归纳法。若
$`x\ne0`$，由引理 E.14(3)，$`x(\bar x/N(x))=(\bar x/N(x))x=1`$；在结合代数中，$`xy=0`$ 于是给出
$`y=(\bar x/N(x))(xy)=0`$。$`\square`$

*计算：* $`\mathbb O`$ 的非结合性与同态 $`M`$——[教程，E.6 节](../../tutorial/zh/E_jordan_classification.md#s-E-6)。

<a id="r-E-16"></a>

**引理 E.16（倍增）。** 设 $`A`$ 是合成代数，$`B\subseteq A`$ 是满足 $`1\in B`$ 的子代数，$`i\in A`$ 满足
$`i\perp B`$ 且 $`N(i)=1`$。对 $`a,b,c,d\in B`$：$`Bi\perp B`$，$`N(a+bi)=N(a)+N(b)`$，并且
```math
(a+bi)(c+di)=(ac-\bar db)+(da+b\bar c)\,i .
```
所以 $`B\oplus Bi`$ 是维数为 $`2\dim B`$ 的子代数，并且只要 $`h:B\to D`$ 是合成代数的同构，$`a+bi\mapsto(h(a),h(b))`$ 就是同构
$`B\oplus Bi\cong\mathrm{CD}(D)`$。此外 $`B`$ 是
结合的；若 $`B\oplus Bi`$ 是结合的，则 $`B`$ 是交换的；若 $`B\oplus Bi`$ 是交换的，则 $`B=\mathbb R1`$。

*证明。* $`B`$ 含有 $`b`$ 时也含有 $`\bar b=2\mathrm{Re}(b)1-b`$。由引理 E.14，$`\langle bi,a\rangle=\langle i,\bar ba\rangle=0`$，
$`N(bi)=N(b)`$，$`\bar i=-i`$，并且对 $`c\in B`$ 有 $`ci+ic=2\mathrm{Re}(c)i`$，即 $`ci=i\bar c`$。引理 E.14(3) 的
恒等式依次给出：取 $`(x,w,y)=(\bar a,i,b)`$，得 $`a(ib)=i(\bar ab)`$；由于 $`ib=\bar bi`$ 且
$`i(\bar ab)=\overline{\bar ab}\,i=(\bar ba)i`$，把 $`b`$ 换成 $`\bar b`$ 得 $`a(bi)=(ba)i`$。取 $`(y,x,w)=(a,i,\bar c)`$，
得 $`(ai)c=(a\bar c)i`$。取 $`x=i`$，得 $`i(iy)=-y`$。取 $`(x,w,y)=(ai,i,\bar d)`$，其中 $`\overline{ai}=-ai`$（因为 $`ai\perp1`$），且
$`\langle ai,i\rangle=\mathrm{Re}\,a`$（引理 E.14(1)），得
```math
-(ai)(i\bar d)-i((ai)\bar d)=2\mathrm{Re}(a)\bar d .
```
这里 $`i\bar d=di`$，$`i((ai)\bar d)=i((ad)i)=i(i\,\overline{ad})=-\bar d\bar a`$，所以
$`(ai)(di)=\bar d\bar a-2\mathrm{Re}(a)\bar d=-\bar da`$。用这四条规则展开 $`(a+bi)(c+di)`$ 即得该公式，
并且在 $`h`$ 下 $`N`$ 与乘积都与 $`\mathrm{CD}(D)`$ 中的相符。由于 $`N`$ 在
$`B\oplus Bi\subseteq A`$ 上是乘性的，引理 E.15(1) 证明中的计算给出：对所有
$`b\in B`$ 有 $`\langle d(ac)-(da)c,b\rangle=0`$，所以 $`B`$ 是结合的。若 $`B\oplus Bi`$ 是结合的，则 $`(ab)i=a(bi)=(ba)i`$，所以 $`ab=ba`$。若它是交换的，
则 $`ai=ia=\bar ai`$，所以 $`a=\bar a\in\mathbb R1`$。$`\square`$

<a id="r-E-17"></a>

**定理 E.17（Hurwitz）。** 每个合成代数都恰好同构于 $`\mathbb R,\mathbb C,\mathbb H,\mathbb O`$ 之一。

*证明。* 令 $`B_0=\mathbb R1\cong\mathbb R`$。只要 $`B_k\ne A`$，就选取满足 $`N(i_{k+1})=1`$ 的 $`i_{k+1}\perp B_k`$，并令
$`B_{k+1}=B_k\oplus B_ki_{k+1}`$；由引理 E.16，$`B_k\cong\mathrm{CD}^k(\mathbb R)`$ 的维数为 $`2^k`$。若 $`B_4`$
存在，引理 E.16 将使 $`B_3`$ 结合，从而 $`B_2`$ 交换，从而 $`B_1=\mathbb R1`$，但后者的维数为
$`1\ne2`$。所以对某个 $`k\le3`$ 有 $`A=B_k`$。这四个代数的维数各不相同。$`\square`$

<a id="r-E-18"></a>

**引理 E.18.** 设 $`A`$ 是合成代数。(1) $`A`$ 的任意两个元素都落在某个结合子代数中。
(2)（**Moufang 恒等式**）对所有 $`u,x,y`$，$`(ux)(yu)=u((xy)u)=(u(xy))u`$。

*证明。* (1) 设 $`p,q\in A`$；不妨设 $`p\notin\mathbb R1`$（否则交换二者，若两者都是实的，就用 $`\mathbb R1`$）。
$`B=\mathrm{span}\{1,p\}`$ 是子代数（$`p^2=2\mathrm{Re}(p)p-N(p)1`$），是维数为 $`2`$ 的合成代数，
所以 $`B\cong\mathbb C`$（定理 E.17）。若 $`q\in B`$，则已证毕；否则设 $`i`$ 是 $`q`$ 正交于 $`B`$ 的分量
归一化后的向量。由引理 E.16，$`B\oplus Bi\cong\mathrm{CD}(\mathbb C)=\mathbb H`$ 是结合的，且含有 $`p,q`$。
(2) 结合子是交错的（引理 E.14(5)），所以 $`(u,xy,u)=0`$，这就是第二个等式。令 $`X=(u,x,y)`$。
由引理 E.14(4)，$`\mathrm{Re}\,X=0`$，并且 $`\langle u,X\rangle=\langle u\bar y,ux\rangle-N(u)\langle1,xy\rangle =N(u)(\langle\bar y,x\rangle-\langle\bar x,y\rangle)=0`$；所以引理 E.14(3) 中极化后的二次恒等式给出
$`uX+Xu=2\mathrm{Re}(u)X`$，即 $`\bar uX=Xu`$。展开定义可知，在任何代数中
```math
(wx,y,z)-(w,xy,z)+(w,x,yz)=w(x,y,z)+(w,x,y)z .
```
取 $`(w,x,y,z)=(u,u,x,y)`$，并利用 $`(u^2,x,y)=2\mathrm{Re}(u)X`$：$`(u,ux,y)=2\mathrm{Re}(u)X-uX=\bar uX=Xu`$。取
$`(w,x,y,z)=(u,x,y,u)`$，并利用 $`(ux,y,u)=(u,ux,y)=Xu`$ 与 $`(x,y,u)=X`$（循环置换）：$`(u,x,yu)=uX+Xu-Xu=uX`$。因此
$`(ux)(yu)-u((xy)u)=(u,x,yu)-u(x,y,u)=uX-uX=0`$。$`\square`$

*计算：* 随机八元数上的 Moufang 恒等式——[教程，E.6 节](../../tutorial/zh/E_jordan_classification.md#s-E-6)。

<a id="r-E-19"></a>

**定理 E.19（Frobenius）。** 每个具有单位元 $`1\ne0`$ 且没有零因子的有限维结合实代数 $`A`$
都同构于 $`\mathbb R`$、$`\mathbb C`$ 或 $`\mathbb H`$。

*证明。* (i) 对 $`x\in A`$，幂 $`1,x,x^2,\dots`$ 线性相关；设 $`\mu`$ 是满足 $`\mu(x)=0`$ 的次数最低的首一实多项式。
若 $`\mu=fg`$ 且 $`\deg f,\deg g\ge1`$，则 $`f(x)g(x)=0`$ 迫使 $`f(x)=0`$ 或 $`g(x)=0`$；所以 $`\mu`$
不可约，因而（代数基本定理）$`\mu=t-\alpha`$ 或 $`\mu=(t-\alpha)^2+\beta^2`$，$`\beta>0`$。令
$`V=\{v\in A:v^2\in(-\infty,0]\cdot1\}`$。则 $`A=\mathbb R1+V`$（在第二种情形中 $`x-\alpha\in V`$），并且
$`\mathbb R1\cap V=0`$。(ii) $`V`$ 是子空间：它在数乘下封闭，我们证明对 $`u,v\in V`$ 有 $`u+v\in V`$。若
$`u,v`$ 线性相关，则 $`u+v`$ 是其中一个的倍数；所以设它们线性无关。则 $`1,u,v`$ 线性无关：若
$`1=\alpha u+\gamma v`$，不妨设 $`\alpha\ne0`$，则 $`u^2=(1-2\gamma v+\gamma^2v^2)/\alpha^2\in\mathbb R1`$ 迫使 $`\gamma=0`$，
从而 $`u\in\mathbb R1\cap V=0`$。由 (i)，$`(u\pm v)^2=p_\pm(u\pm v)+q_\pm1`$，其中 $`p_\pm,q_\pm`$ 是实数；相加得
$`2u^2+2v^2=(p_++p_-)u+(p_+-p_-)v+(q_++q_-)1`$，所以 $`p_\pm=0`$。若 $`q_+>0`$，则
$`(u+v-\sqrt{q_+})(u+v+\sqrt{q_+})=0`$，从而 $`u+v\in\mathbb R1`$，这不可能；所以 $`u+v\in V`$。(iii) 对 $`u,v\in V`$，
$`uv+vu=(u+v)^2-u^2-v^2=-2\langle u,v\rangle1`$ 定义了一个对称双线性型，满足 $`\langle u,u\rangle=-u^2\ge0`$，并且
仅当 $`u^2=0`$，即 $`u=0`$ 时 $`\langle u,u\rangle=0`$：这是一个内积。正交的 $`u,v\in V`$ 反交换。(iv) 若 $`V=0`$，
则 $`A=\mathbb R`$。若 $`\dim V=1`$，则 $`V=\mathbb Ri`$，$`i^2=-1`$，且 $`A\cong\mathbb C`$。若 $`\dim V\ge2`$，取规范正交的
$`i,j\in V`$ 并令 $`k=ij`$。则 $`k^2=ijij=-i^2j^2=-1`$，所以 $`k\in V`$，$`ik+ki=-j+j=0`$，$`jk+kj=i-i=0`$：$`i,j,k`$ 是
规范正交的，且 $`ij=-ji=k`$，$`jk=-kj=i`$，$`ki=-ik=j`$。若存在与 $`i,j,k`$ 正交的单位向量 $`l\in V`$，则
$`lk=(li)j=-i(lj)=(ij)l=kl`$，而 $`lk=-kl`$，所以 $`kl=0`$，这不可能。因此 $`A=\mathrm{span}\{1,i,j,k\}`$。在 $`\mathbb H`$ 中，
元素 $`I=(i,0)`$、$`J=(0,1)`$、$`K=(0,i)`$ 有相同的乘法表（引理 E.15(2) 给出 $`IJ=K=-JI`$，而
$`I^2=J^2=-1`$、$`K=IJ`$ 借助结合性决定其余部分），所以 $`A\cong\mathbb H`$。$`\square`$

<a id="s-E-7"></a>

## E.7 坐标化

对 $`D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}`$、$`d=\dim D`$ 与 $`n\ge1`$，令
$`\mathrm{Herm}_n(D)=\{X\in M_n(D):X_{ji}=\overline{X_{ij}}\}`$，乘积为 $`X\circ Y=\frac12(XY+YX)`$，其中
$`(XY)_{ik}=\sum_jX_{ij}Y_{jk}`$（定义它不需要 $`D`$ 的结合性）。对角元满足 $`\bar x=x`$，所以
它们是实的。设 $`E_{ij}`$ 是矩阵单位，对 $`i\ne j`$ 与 $`a\in D`$ 令 $`a[ij]=aE_{ij}+\bar aE_{ji}`$；于是
$`a[ji]=\bar a[ij]`$。

<a id="r-E-20"></a>

**引理 E.20.** $`\mathrm{Herm}_n(D)`$ 在 $`\circ`$ 下封闭，有单位元 $`\sum_iE_{ii}`$，维数为
$`n+d\frac{n(n-1)}2`$，并且对 $`a,b\in D`$、前两个公式中任意的 $`i,j,k`$ 以及其余公式中互不相同的 $`i,j,k,l`$（$`a[ij]`$ 中总有 $`i\ne j`$）：
```math
E_{ii}\circ E_{jj}=\delta_{ij}E_{ii},\quad E_{kk}\circ a[ij]=\tfrac12(\delta_{ki}+\delta_{kj})a[ij],\quad
a[ij]\circ b[ij]=\langle a,b\rangle(E_{ii}+E_{jj}),\quad a[ij]\circ b[jk]=\tfrac12(ab)[ik],\quad a[ij]\circ b[kl]=0 .
```
此外 $`(X|Y)=\mathrm{Re}\,\mathrm{tr}(XY)`$ 是满足 $`(X\circ Y|Z)=(Y|X\circ Z)`$ 的内积。

*证明。* 矩阵乘法：$`a[ij]\,b[ij]=a\bar bE_{ii}+\bar abE_{jj}`$，且 $`a\bar b+b\bar a=\bar ab+\bar ba=2\langle a,b\rangle`$
（把 $`x\bar x=\bar xx=N(x)`$ 极化）；$`a[ij]\,b[jk]=abE_{ik}`$，$`b[jk]\,a[ij]=\bar b\bar aE_{ki}=\overline{ab}E_{ki}`$；
其他乘积是直接的。所有基元素的乘积都属于 $`\mathrm{Herm}_n(D)`$，因此所有乘积都属于它。其次，
对 $`X\ne0`$，$`(X|X)=\sum_{i,j}\mathrm{Re}(X_{ij}\overline{X_{ij}})=\sum_{i,j}N(X_{ij})>0`$。令
$`T(X,Y,Z)=\mathrm{Re}\,\mathrm{tr}((XY)Z)=\sum_{i,j,k}\mathrm{Re}((X_{ij}Y_{jk})Z_{ki})`$。由引理 E.14(4)，
$`\mathrm{Re}\,\mathrm{tr}(AB)=\mathrm{Re}\,\mathrm{tr}(BA)`$，$`T(X,Y,Z)=T(Z,X,Y)`$，并且取共轭
（$`\mathrm{Re}\,\bar p=\mathrm{Re}\,p`$，$`\overline{(pq)r}=\bar r(\bar q\bar p)`$），再利用 $`X,Y,Z`$ 是厄米的，
得 $`T(X,Y,Z)=\sum\mathrm{Re}(Z_{ik}(Y_{kj}X_{ji}))=T(Z,Y,X)`$。所以 $`T`$ 是对称的，从而
$`(X\circ Y|Z)=\frac12(T(X,Y,Z)+T(Y,X,Z))=T(X,Y,Z)`$ 关于 $`X,Y,Z`$ 对称。$`\square`$

*计算：* 由结构常数给出的 $`\mathrm{Herm}_n(D)`$（函数 `herm`）——[教程，E.3 节](../../tutorial/zh/E_jordan_classification.md#s-E-3)；$`\mathrm{Herm}_3(\mathbb O)`$ 中的内积 $`(X|Y)`$——[教程，E.8 节](../../tutorial/zh/E_jordan_classification.md#s-E-8)。

<a id="r-E-21"></a>

**引理 E.21（定向）。** 设 $`P_{12},P_{13},P_{23}`$ 是内积空间，$`P_{ji}=P_{ij}`$，并且对 $`\{1,2,3\}`$ 中互不相同元素的每个
有序三元组 $`(i,j,k)`$，设 $`m:P_{ij}\times P_{jk}\to P_{ik}`$ 是双线性的，使得对
$`x\in P_{ij}`$、$`y\in P_{jk}`$、$`z\in P_{ik}`$ 有 $`m(x,y)=m(y,x)`$（右边是对应于 $`(k,j,i)`$ 的映射）以及 $`\langle m(x,y),z\rangle=\langle y,m(x,z)\rangle`$。设 $`D`$ 是合成代数，$`\iota_{ij}:D\to P_{ij}`$（$`i<j`$）是线性的
等距双射，满足 $`m(\iota_{12}a,\iota_{23}b)=\iota_{13}(ab)`$；令 $`\iota_{ji}(a)=\iota_{ij}(\bar a)`$。则对所有有序三元组和所有 $`a,b\in D`$ 有
$`m(\iota_{ij}a,\iota_{jk}b)=\iota_{ik}(ab)`$。

这些假设在以下情形成立：(i) 对单的 $`J`$ 的若尔当标架中的三个指标，取 Peirce 空间与 $`m=2\circ`$
（交换性与引理 E.12(2)）；(ii) 在 $`\mathrm{Herm}_3(D)`$ 中，取 $`P_{ij}=\{a[ij]\}`$、
$`\langle a[ij],b[ij]\rangle=\langle a,b\rangle`$（定义合理，因为共轭是等距的）以及 $`m=2\circ`$，因为
$`m(a[ij],b[jk])=(ab)[ik]`$，并且对 $`c[ik]`$ 有 $`\langle ab,c\rangle=\langle b,\bar ac\rangle`$（引理 E.14(2)）。

*证明。* 若断言对 $`(i,j,k)`$ 成立，则它对 $`(k,j,i)`$ 也成立：
$`m(\iota_{kj}a,\iota_{ji}b)=m(\iota_{ij}\bar b,\iota_{jk}\bar a)=\iota_{ik}(\bar b\bar a)=\iota_{ik}(\overline{ab})=\iota_{ki}(ab)`$。
所以在 $`(1,2,3)`$ 已知的前提下，只需处理 $`(2,1,3)`$ 与 $`(1,3,2)`$。对所有 $`c\in D`$，利用伴随性、$`(1,2,3)`$ 以及
引理 E.14(2)，
```math
\langle m(\iota_{21}a,\iota_{13}b),\iota_{23}c\rangle=\langle\iota_{13}b,m(\iota_{12}\bar a,\iota_{23}c)\rangle=\langle b,\bar ac\rangle=\langle ab,c\rangle=\langle\iota_{23}(ab),\iota_{23}c\rangle,
```
```math
\langle m(\iota_{13}a,\iota_{32}b),\iota_{12}c\rangle=\langle\iota_{13}a,m(\iota_{12}c,\iota_{23}\bar b)\rangle=\langle a,c\bar b\rangle=\langle ab,c\rangle=\langle\iota_{12}(ab),\iota_{12}c\rangle,
```
而 $`\iota_{23},\iota_{12}`$ 是满射。（在第二行中，先用 $`m(x,y)=m(y,x)`$，再对处于 $`x`$ 位置的
$`y=\iota_{32}b`$ 使用伴随性。）$`\square`$

<a id="r-E-22"></a>

**定理 E.22（坐标化）。** 设 $`J`$ 是单的，具有若尔当标架 $`(c_1,\dots,c_r)`$，$`r\ge3`$。则存在
维数为 $`d`$ 的合成代数 $`D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}`$ 以及满足 $`\Phi(E_{ii})=c_i`$ 的代数同构
$`\Phi:\mathrm{Herm}_r(D)\to J`$。若 $`r\ge4`$，则 $`D`$ 是结合的，即 $`D\ne\mathbb O`$。

*证明。* *第 1 步（一个三角形）。* 选取满足 $`N=1`$ 的 $`e_{12}\in J_{12}`$、$`e_{23}\in J_{23}`$（命题 E.11）。在
$`D'=J_{13}`$ 上定义 $`p\cdot q=m(m(p,e_{23}),m(e_{12},q))`$，这里用到乘积 $`J_{13}\times J_{32}\to J_{12}`$、
$`J_{21}\times J_{13}\to J_{23}`$ 与 $`J_{12}\times J_{23}\to J_{13}`$。由引理 E.12(1)，$`1'=m(e_{12},e_{23})`$ 是单位元：
$`p\cdot1'=m(m(p,e_{23}),e_{23})=p`$，$`1'\cdot q=m(e_{12},m(e_{12},q))=q`$；并且 $`N(p\cdot q)=N(p)N(q)`$。所以 $`(D',N)`$ 是
合成代数，由定理 E.17，存在同构 $`h:D\to D'`$，其中 $`D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}`$。
令 $`\iota_{13}=h`$、$`\iota_{12}=m(h(\cdot),e_{23})`$、$`\iota_{23}=m(e_{12},h(\cdot))`$；它们是等距双射
（引理 E.12(1)），并且 $`m(\iota_{12}a,\iota_{23}b)=h(a)\cdot h(b)=\iota_{13}(ab)`$。由引理 E.21，对 $`k=3`$ 有：

> $`(\mathrm T_k)`$ 存在等距双射 $`\iota_{ij}:D\to J_{ij}`$，$`i<j\le k`$，使得在 $`\iota_{ji}(a)=\iota_{ij}(\bar a)`$ 之下，
> 对所有互不相同的 $`i,j,l\le k`$ 和所有 $`a,b\in D`$ 有 $`m(\iota_{ij}a,\iota_{jl}b)=\iota_{il}(ab)`$。

*第 2 步（添加一个指标）。* 设 $`(\mathrm T_{k-1})`$ 成立，$`4\le k\le r`$。选取满足 $`N(v)=1`$ 的 $`v\in J_{1k}`$，并令
对 $`2\le i<k`$ 有 $`\iota_{ik}(a)=m(\iota_{i1}a,v)`$，以及 $`\iota_{1k}(a)=m(\iota_{12}a,\iota_{2k}1)`$；由
引理 E.12(1)，它们是等距双射。利用引理 E.12(3)（每次都是四个互不相同的指标）与 $`(\mathrm T_{k-1})`$：对互不相同的
$`i,j\in\{2,\dots,k-1\}`$，
```math
m(\iota_{ij}a,\iota_{jk}b)=m(\iota_{ij}a,m(\iota_{j1}b,v))=m(m(\iota_{ij}a,\iota_{j1}b),v)=m(\iota_{i1}(ab),v)=\iota_{ik}(ab).
```
特别地 $`\iota_{2k}(b)=m(\iota_{23}b,\iota_{3k}1)`$，因此
```math
m(\iota_{12}a,\iota_{2k}b)=m(m(\iota_{12}a,\iota_{23}b),\iota_{3k}1)=m(\iota_{13}(ab),\iota_{3k}1)
=m(m(\iota_{12}(ab),\iota_{23}1),\iota_{3k}1)=m(\iota_{12}(ab),\iota_{2k}1)=\iota_{1k}(ab).
```
对 $`3\le j<k`$，$`\iota_{1j}a=m(\iota_{12}a,\iota_{2j}1)`$，所以
$`m(\iota_{1j}a,\iota_{jk}b)=m(\iota_{12}a,m(\iota_{2j}1,\iota_{jk}b))=m(\iota_{12}a,\iota_{2k}b)=\iota_{1k}(ab)`$。所以该关系
对所有三元组 $`i<j<k`$ 成立，再对每个三角形 $`\{i,j,k\}`$ 应用引理 E.21，即得 $`(\mathrm T_k)`$。

*第 3 步（结合性）。* 若 $`r\ge4`$，由 $`(\mathrm T_4)`$ 与引理 E.12(3)，
$`\iota_{14}((ab)c)=m(m(\iota_{12}a,\iota_{23}b),\iota_{34}c)=m(\iota_{12}a,m(\iota_{23}b,\iota_{34}c))=\iota_{14}(a(bc))`$。

*第 4 步（同构）。* 借助 $`(\mathrm T_r)`$ 定义线性映射
$`\Phi(\sum_i\alpha_iE_{ii}+\sum_{i<j}a_{ij}[ij])=\sum_i\alpha_ic_i+\sum_{i<j}\iota_{ij}(a_{ij})`$；由
定理 E.8(1) 它是双射，并且对所有 $`i\ne j`$ 有 $`\Phi(a[ij])=\iota_{ij}(a)`$。由定理 E.8 与 $`(\mathrm T_r)`$，
$`c_i\circ c_j=\delta_{ij}c_i`$，$`c_k\circ\iota_{ij}(a)=\frac12(\delta_{ki}+\delta_{kj})\iota_{ij}(a)`$，
$`\iota_{ij}(a)\circ\iota_{ij}(b)=\langle a,b\rangle(c_i+c_j)`$，$`\iota_{ij}(a)\circ\iota_{jk}(b)=\frac12\iota_{ik}(ab)`$，以及
$`\iota_{ij}(a)\circ\iota_{kl}(b)=0`$：这正是引理 E.20 的乘法表。由双线性，$`\Phi(X\circ Y)=\Phi(X)\circ\Phi(Y)`$。$`\square`$

<a id="s-E-8"></a>

## E.8 Albert 代数

<a id="r-E-23"></a>

**引理 E.23（$`\mathrm{Herm}_3(D)`$ 的自同构）。** 设 $`D`$ 是合成代数，$`f_{12},f_{23},f_{13}`$ 是
$`D`$ 的等距线性双射，且对所有 $`a,b`$ 有 $`f_{13}(ab)=f_{12}(a)f_{23}(b)`$。则
```math
\Phi\Big(\sum_i\alpha_iE_{ii}+\sum_{i<j}a_{ij}[ij]\Big)=\sum_i\alpha_iE_{ii}+\sum_{i<j}f_{ij}(a_{ij})[ij]
```
是 $`(\mathrm{Herm}_3(D),\circ)`$ 的自同构。对每个满足 $`N(u)=1`$ 的 $`u\in D`$，
$`(f_{12},f_{23},f_{13})=(L_u,R_u,B_u)`$ 就是这样的三元组，其中 $`L_ua=ua`$，$`R_ua=au`$，$`B_ua=uau`$。

*证明。* 在例子 (ii) 中对 $`\iota_{ij}(a)=f_{ij}(a)[ij]`$ 应用引理 E.21：于是对所有 $`i\ne j`$ 有 $`\Phi(a[ij])=\iota_{ij}(a)`$，
并且对所有互不相同的 $`i,j,k`$ 有 $`2\iota_{ij}(a)\circ\iota_{jk}(b)=\iota_{ik}(ab)`$。由于 $`\Phi`$ 固定诸 $`E_{ii}`$，且诸 $`f_{ij}`$
是等距的，$`\Phi`$ 保持引理 E.20 乘法表中的每个乘积，因而保持 $`\circ`$。至于该三元组，$`L_u,R_u,B_u`$
是等距的（$`N(ua)=N(a)`$）、双射的，而 $`u(ab)u=(ua)(bu)`$ 就是 Moufang 恒等式（引理 E.18(2)）。$`\square`$

*计算：* 对随机的单位八元数 $`u`$ 检验自同构 $`(L_u,R_u,B_u)`$——[教程，E.8 节](../../tutorial/zh/E_jordan_classification.md#s-E-8)。

<a id="r-E-24"></a>

**定理 E.24.** (1) 对 $`D\in\{\mathbb R,\mathbb C,\mathbb H\}`$ 与 $`n\ge1`$，$`\mathrm{Herm}_n(D)`$ 是欧几里得若尔当
代数，并且 $`\mathrm{Herm}_n(\mathbb H)`$ 同构于 $`\mathrm{Herm}_{2n}(\mathbb C)`$ 的一个若尔当子代数。
(2) **Albert 代数** $`\mathrm{Herm}_3(\mathbb O)`$（维数 $`27`$）与 $`\mathrm{Herm}_2(\mathbb O)`$ 是
欧几里得若尔当代数。(3) 在其中每一个代数中，$`(E_{11},\dots,E_{nn})`$ 是若尔当标架，其 Peirce 空间
$`J_{ij}=\{a[ij]:a\in D\}`$ 的维数为 $`d=\dim D`$；该代数是单的；对 $`D=\mathbb R,\mathbb C,\mathbb H,\mathbb O`$ 有 $`\mathrm{Herm}_1(D)\cong\mathbb R`$ 与
$`\mathrm{Herm}_2(D)\cong V_{d+1}`$。

*证明。* (1) $`M_n(D)`$ 是结合实代数（引理 E.15(2)），所以由引理 E.20 与引理 13.1，
$`\mathrm{Herm}_n(D)`$ 是若尔当代数，并由引理 E.20 是欧几里得的。把每个矩阵元 $`X_{ij}\in\mathbb H`$ 换成
$`2\times2`$ 块 $`M(X_{ij})`$，得到一个单的代数同态 $`M_n(\mathbb H)\to M_{2n}(\mathbb C)`$（引理 E.15(2)），它
把厄米矩阵映为厄米矩阵，因为 $`M(\bar q)=M(q)^*`$；所以它保持 $`\circ`$。

(2) 在 $`\mathrm{Herm}_3(\mathbb O)`$ 上令 $`J(x,y)=(x^2\circ y)\circ x-x^2\circ(y\circ x)`$；自同构 $`\Phi`$ 满足
$`J(\Phi x,\Phi y)=\Phi J(x,y)`$。设 $`x=\sum_i\alpha_iE_{ii}+p[12]+q[23]+s[13]`$。我们复合三个具有
引理 E.23 形式的自同构，使得 $`(1,3)`$ 与 $`(2,3)`$ 处的矩阵元变为实数。(a) 若 $`s\notin\mathbb R`$，写
$`s=|s|(\cos\theta+v\sin\theta)`$，其中 $`v\perp1`$，$`N(v)=1`$，于是 $`v^2=-1`$（引理 E.14(3)）；对
$`u=\cos\frac\theta2-v\sin\frac\theta2`$，在交换子代数 $`\mathrm{span}\{1,v\}\cong\mathbb C`$ 中计算，得
$`B_us=u^2s=|s|`$。(b) 现在 $`s`$ 是实的。选取满足 $`N(u)=1`$ 且 $`u\perp q`$ 的 $`u\perp1`$（因为 $`\dim\mathbb O=8`$，这是可能的）；
则 $`B_us=su^2=-s`$ 是实的，并且 $`\mathrm{Re}(qu)=\langle u,\bar q\rangle=0`$，因为 $`\bar q\in\mathrm{span}\{1,q\}`$。
(c) 现在 $`q\perp1`$；若 $`q\ne0`$，取 $`u=q/|q|`$：$`qu=q^2/|q|=-|q|`$ 与 $`B_us=-s`$ 都是实的。所以存在自同构
$`\Phi`$，使得 $`\Phi x`$ 除了 $`(1,2)`$ 处可能的 $`p'`$ 之外，矩阵元都是实的。由于 $`J(\Phi x,\cdot)`$ 是线性的，只需
证明对 $`y=E_{ii}`$ 以及 $`y=c[ij]`$（$`c\in\mathbb O`$）有 $`J(\Phi x,y)=0`$。由引理 E.18(1)，$`p'`$ 与 $`c`$ 落在一个结合
子代数 $`Q`$ 中；于是 $`\Phi x,y\in\mathrm{Herm}_3(Q)`$，它是 $`\mathrm{Herm}_3(\mathbb O)`$ 的子代数，位于
结合代数 $`M_3(Q)`$ 之中，所以由引理 13.1，$`J(\Phi x,y)=0`$。因此对所有 $`y`$ 有 $`J(x,\Phi^{-1}y)=\Phi^{-1}J(\Phi x,y)=0`$：
若尔当恒等式成立，而引理 E.20 使该代数成为欧几里得的。第三行与第三列为零的矩阵构成
一个子代数（引理 E.20），其单位元为 $`E_{11}+E_{22}`$，这就是 $`\mathrm{Herm}_2(\mathbb O)`$。

(3) 由引理 E.20，诸 $`E_{ii}`$ 是两两正交的幂等元，其和为 $`1`$；$`L(E_{kk})`$ 在 $`E_{jj}`$ 上作用为乘以
$`\delta_{kj}`$，在 $`a[ij]`$ 上作用为乘以 $`\frac12(\delta_{ki}+\delta_{kj})`$，所以 $`J_1(E_{ii})=\mathbb RE_{ii}`$，$`J_{1/2}(E_{ii})\cap J_{1/2}(E_{jj})=\{a[ij]\}`$。$`E_{ii}`$ 是本原的：若 $`E_{ii}=g+g'`$ 是两个正交的非零幂等元之和，则
$`g=E_{ii}\circ g\in\mathbb RE_{ii}`$，所以 $`g=E_{ii}`$，$`g'=0`$。命题 E.11 的图是完全图，所以该代数是
单的。对 $`n=1,2`$ 应用定理 E.13。$`\square`$

*计算：* Albert 代数的若尔当恒等式，以及在一个随机元素上执行 (2) 的证明中的步骤 (a)–(c)——[教程，E.8 节](../../tutorial/zh/E_jordan_classification.md#s-E-8)。

<a id="r-E-25"></a>

**推论 E.25.** 对 $`n\ge4`$，$`\mathrm{Herm}_n(\mathbb O)`$ 不是若尔当代数。

*证明。* 若它是若尔当代数，则引理 E.10(2) 适用：其证明只用到若尔当恒等式和引理 E.20 中的乘积
（取 $`c_k=E_{kk}`$）。于是对所有 $`a,b,c\in\mathbb O`$ 有 $`\frac14((ab)c)[14]=(a[12]\circ b[23])\circ c[34]=a[12]\circ(b[23]\circ c[34])=\frac14(a(bc))[14]`$，
与引理 E.15(2) 矛盾。$`\square`$

*计算：* 若尔当恒等式在 $`\mathrm{Herm}_4(\mathbb O)`$ 中不成立——[教程，E.8 节](../../tutorial/zh/E_jordan_classification.md#s-E-8)。

<a id="r-E-26"></a>

**定理 E.26（Albert）。** 设 $`D`$ 是合成代数。若存在结合实代数 $`A`$ 以及
满足 $`\varphi(x\circ y)=\frac12(\varphi(x)\varphi(y)+\varphi(y)\varphi(x))`$ 的线性单射 $`\varphi:\mathrm{Herm}_3(D)\to A`$，
则 $`D`$ 是结合的。特别地，Albert 代数是**例外的**：它不同构于任何结合代数（例如 $`M_N(\mathbb C)`$）的
若尔当子代数。

*证明。* 令 $`e_i=\varphi(E_{ii})`$，$`e=e_1+e_2+e_3=\varphi(1)`$。则 $`e_i^2=e_i`$，且对 $`i\ne j`$ 有 $`e_ie_j+e_je_i=0`$；
分别从左边、从右边乘以 $`e_i`$，得 $`e_ie_j+e_ie_je_i=0=e_ie_je_i+e_je_i`$，所以
$`e_ie_j=e_je_i=0`$。对 $`x=\varphi(X)`$，$`ex+xe=2x`$ 给出 $`exe=ex=xe=x`$，所以 $`x=\sum_{p,q}x_{pq}`$，其中
$`x_{pq}=e_pxe_q`$。对 $`X=a[ij]`$，关系 $`e_i\circ x=e_j\circ x=\frac12x`$ 与 $`e_k\circ x=0`$（$`k\ne i,j`$）表明
$`\frac12([p=i]+[q=i])x_{pq}=\frac12x_{pq}`$ 等等，这样只剩下 $`x_{ij}`$ 与 $`x_{ji}`$。对有序对定义
$`\alpha_{ij}(a)=e_i\varphi(a[ij])e_j`$。对互不相同的 $`i,j,k`$，$`\varphi(a[ij])\varphi(b[jk])=\alpha_{ij}(a)\alpha_{jk}(b)`$ 是
$`\varphi(a[ij])\varphi(b[jk])+\varphi(b[jk])\varphi(a[ij])=\varphi((ab)[ik])`$ 中唯一属于 $`e_iAe_k`$ 的项，所以
```math
\alpha_{ij}(a)\,\alpha_{jk}(b)=\alpha_{ik}(ab).
```
由 $`\varphi(1[13])^2=\varphi(E_{11}+E_{33})`$ 得 $`\alpha_{13}(1)\alpha_{31}(1)=e_1`$，$`\alpha_{31}(1)\alpha_{13}(1)=e_3`$；
由 $`\varphi(a[13])^2=N(a)(e_1+e_3)`$ 得 $`\alpha_{13}(a)\alpha_{31}(\bar a)=N(a)e_1`$，所以 $`\alpha_{13}`$ 是单射。现在
```math
\alpha_{13}(ay)=\alpha_{12}(a)\alpha_{23}(y)=\alpha_{13}(a)\alpha_{32}(1)\alpha_{21}(1)\alpha_{13}(y)=\alpha_{13}(a)\alpha_{31}(1)\alpha_{13}(y),
```
因此，令 $`\pi=\alpha_{31}(1)`$ 并利用 $`A`$ 中的结合性，
$`\alpha_{13}((ab)c)=\alpha_{13}(a)\pi\alpha_{13}(b)\pi\alpha_{13}(c)=\alpha_{13}(a(bc))`$。所以 $`(ab)c=a(bc)`$。由引理 E.15(2)，
$`\mathbb O`$ 不是结合的。$`\square`$

Albert 的原始论证（A. A. Albert, Ann. of Math. 35 (1934) 65–73）就是这一类型的。

<a id="s-E-9"></a>

## E.9 分类

<a id="r-E-27"></a>

**定理 E.27（Jordan–von Neumann–Wigner）。**

1. 每个有限维形式实（等价地，欧几里得）若尔当代数都是单理想的直和，
   而每个单的这样的代数都同构于某个 $`n\ge1`$ 的 $`\mathrm{Herm}_n(\mathbb R)`$、$`\mathrm{Herm}_n(\mathbb C)`$ 或
   $`\mathrm{Herm}_n(\mathbb H)`$，或同构于 $`m\ge2`$ 的自旋因子 $`V_m`$，或同构于 Albert 代数
   $`\mathrm{Herm}_3(\mathbb O)`$。
2. 代数 $`\mathbb R`$；$`V_m`$（$`m\ge2`$）；$`\mathrm{Herm}_n(D)`$（$`D=\mathbb R,\mathbb C,\mathbb H`$，$`n\ge3`$）；
   $`\mathrm{Herm}_3(\mathbb O)`$ 都是单的、形式实的，且两两不同构，所以每个单的形式
   实若尔当代数都恰好同构于其中之一。此外，对 $`D=\mathbb R,\mathbb C,\mathbb H,\mathbb O`$ 有 $`\mathrm{Herm}_1(D)\cong\mathbb R`$ 与
   $`\mathrm{Herm}_2(D)\cong V_{\dim D+1}`$。
3. 在单的形式实若尔当代数中，所有若尔当标架有相同的元素个数 $`r`$，它们的所有 Peirce
   空间 $`J_{ij}`$ 有相同的维数 $`d`$：对 $`\mathbb R`$ 有 $`r=1`$，对 $`V_m`$ 有 $`(r,d)=(2,m-1)`$，对
   $`\mathrm{Herm}_n(D)`$（$`n\ge2`$）为 $`(n,\dim D)`$。

*证明。* (1) 由定理 E.9，$`J`$ 是单的形式实理想的直和。这样的理想有某个大小 $`r`$ 的若尔当标架
（引理 E.6(2)）；若 $`r=1`$，它 $`\cong\mathbb R=\mathrm{Herm}_1(\mathbb R)`$；若 $`r=2`$，它 $`\cong V_{d+1}`$
（定理 E.13）；否则它 $`\cong\mathrm{Herm}_r(D)`$，其中若 $`r\ge4`$ 则 $`D\ne\mathbb O`$，若 $`r=3`$ 则 $`D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}`$
（定理 E.22）。

(2) 关于 $`\mathrm{Herm}_n(D)`$，见定理 E.24。由引理 13.3，$`V_m`$ 是若尔当代数，并且对
$`\langle(a,u),(b,v)\rangle=ab+u\cdot v`$ 是欧几里得的，因为 $`\langle x\circ y,z\rangle=abc+c\,u\cdot v+a\,v\cdot w+b\,u\cdot w`$ 关于
$`x=(a,u)`$、$`y=(b,v)`$、$`z=(c,w)`$ 对称。对单位向量 $`w`$，$`c_{\pm w}=\frac12(1,\pm w)`$ 是若尔当标架
（定理 13.4(1)；$`c_{\pm w}`$ 是本原的，因为两个正交非零幂等元之和 $`c_{w'}+c_{-w'}`$ 等于 $`1`$），
其 $`J_{12}=\{(0,v):v\perp w\}`$，维数为 $`m-1\ge1`$（各维数之和为 $`m+1`$）；所以由命题 E.11，$`V_m`$
是单的。为了区分这些代数，设 $`s(J)`$ 是两两正交的非零
幂等元的最大个数（由引理 E.6(2) 它是有限的）；$`s`$ 与 $`\dim`$ 在同构下不变。显然 $`s(\mathbb R)=1`$。在
$`V_m`$ 中，非零幂等元是 $`1`$ 与诸 $`c_w`$，其中 $`1\circ c_w=c_w`$，而 $`c_w\circ c_{w'}=0`$ 仅当 $`w'=-w`$
（定理 13.4(1)）；所以 $`s(V_m)=2`$。对 $`D\in\{\mathbb R,\mathbb C,\mathbb H\}`$，$`s(\mathrm{Herm}_n(D))=n`$：诸
$`E_{ii}`$ 给出 $`s\ge n`$。反之，正交幂等元 $`p_k\ne0`$ 对 $`k\ne l`$ 满足 $`p_kp_l=0`$（同
定理 E.26 的证明），所以 $`D^n`$ 的右 $`D`$-子模 $`V_k=p_kD^n`$ 非零且独立：若 $`\sum_kv_k=0`$，
$`v_k\in V_k`$，则 $`v_l=p_l\sum_kv_k=0`$。非零的 $`V_k`$ 含有某个 $`v\ne0`$，从而含有实维数为 $`d`$ 的 $`vD`$
（$`D`$ 没有零因子）；所以 $`sd\le nd`$。对 $`\mathrm{Herm}_3(\mathbb O)`$，诸 $`E_{ii}`$ 给出 $`s\ge3`$；而 $`s`$ 个
两两正交的非零幂等元，连同 $`1`$ 减去它们之和，可以如定理 E.7 的证明那样细分为一个
有 $`r'\ge s`$ 个元素的若尔当标架。由命题 E.11(2) 与定理 E.22，$`27=r'+d'\frac{r'(r'-1)}2`$，其中 $`d'\ge1`$，并且
若 $`r'\ge4`$，则 $`d'\in\{1,2,4\}`$。对 $`r'=4,5,6`$，右边分别属于 $`\{10,16,28\}`$、$`\{15,25,45\}`$、$`\{21,36,66\}`$，而
对 $`r'\ge7`$，它超过 $`27`$；所以 $`s\le r'=3`$。数对 $`(s,\dim)`$ 为 $`(1,1)`$、$`(2,m+1)`$、$`(n,\frac{n(n+1)}2)`$、
$`(n,n^2)`$、$`(n,n(2n-1))`$（$`n\ge3`$）以及 $`(3,27)`$，它们互不相同，因为 $`\frac{n(n+1)}2<n^2<n(2n-1)`$ 且 $`6,9,15\ne27`$。
秩 $`1`$ 与秩 $`2`$ 中的同构就是定理 E.24(3)。

(3) 若 $`F`$ 是单的 $`J`$ 的一个有 $`r_F`$ 个元素的若尔当标架，则用 $`F`$ 来执行 (1) 的证明，就给出 $`J`$ 与
(2) 中一个秩为 $`r_F`$ 的模型之间的同构，而该模型的 $`s`$ 等于 $`r_F`$。所以对每个标架有 $`r_F=s(J)`$，从而
对 $`r\ge2`$ 有 $`d=(\dim J-r)/\binom r2`$（命题 E.11(2)）。$`\square`$

*计算：* 对每种类型计算若尔当标架、(3) 中的 Peirce 维数 $`(r,d)`$ 以及 $`\dim J=r+d\,\frac{r(r-1)}2`$——[教程，E.9 节](../../tutorial/zh/E_jordan_classification.md#s-E-9)。

<a id="r-E-28"></a>

**推论 E.28（秩与迹）。** 在每个形式实若尔当代数中，所有若尔当标架都有相同的
元素个数，称为**秩**。若 $`J`$ 是秩为 $`r`$ 的单代数，则 $`\mathrm{tr}\,x=\frac r{\dim J}\mathrm{Tr}\,L(x)`$ 对每个本原幂等元 $`c`$ 满足
$`\mathrm{tr}\,c=1`$，因此对定理 E.7 的每个分解有 $`\mathrm{tr}(\sum_i\lambda_ic_i)=\sum_i\lambda_i`$，
并且 $`\mathrm{tr}(x\circ y)`$ 是一个结合的内积。

*证明。* 设 $`J=\bigoplus_pJ^{(p)}`$ 如定理 E.9 所述，$`c`$ 是本原的。它的分量 $`c^{(p)}\in J^{(p)}`$ 是
两两正交的幂等元（$`c^2=\sum_p(c^{(p)})^2`$），所以恰有一个非零。因此 $`J`$ 的一个若尔当标架中
落在 $`J^{(p)}`$ 里的成员构成 $`J^{(p)}`$ 的一个若尔当标架（它们之和为 $`1`$ 的分量 $`1^{(p)}`$），而 $`J`$ 的秩
是诸 $`J^{(p)}`$ 的秩之和（定理 E.27(3)）。在单的 $`J`$ 中，本原的 $`c`$ 属于某个若尔当标架
（如定理 E.7 的证明那样细分 $`1-c`$），所以由命题 E.11(2)，$`\mathrm{Tr}\,L(c)=t=\dim J/r`$（对 $`r=1`$，
$`J=\mathbb Rc`$，$`\mathrm{Tr}\,L(c)=1`$）。最后，$`\mathrm{tr}(x\circ y)=\frac r{\dim J}\tau(x,y)`$（定理 E.5）。$`\square`$

<a id="s-E-10"></a>

## E.10 本附录为第 13 章证明了什么

| 第 13 章中的陈述 | 在此处的证明 |
|---|---|
| 形式实 $`\Leftrightarrow`$ 欧几里得（13.4 节） | 定理 E.5 |
| 定理 13.7：分解为单理想、列表、“恰好其中之一” | 定理 E.9、E.27 |
| 对 $`\mathbb K=\mathbb R,\mathbb C,\mathbb H,\mathbb O`$，$`\mathrm{Herm}_2(\mathbb K)\cong V_{\dim\mathbb K+1}`$ | 定理 E.24(3) |
| $`\mathrm{Herm}_r(\mathbb H)`$ 嵌入 $`\mathrm{Herm}_{2r}(\mathbb C)`$ | 定理 E.24(1) |
| Albert 代数是若尔当代数，并且是例外的 | 定理 E.24(2)、E.26 |
| $`\mathrm{Herm}_r(\mathbb O)`$ 仅当 $`r\le3`$ 时是若尔当代数 | 定理 E.24(2)、推论 E.25 |
| 借助若尔当标架的谱定理、秩、$`\mathrm{tr}\,z=\sum_i\mu_i`$、$`\mathrm{tr}(x\circ y)`$ 的结合性（13.6 节） | 定理 E.7、推论 E.28 |

**注记。** 该分类归功于 P. Jordan、J. von Neumann 和 E. Wigner，*On an algebraic generalization of the
quantum mechanical formalism*, Ann. of Math. 35 (1934) 29–64，而 $`\mathrm{Herm}_3(\mathbb O)`$ 的例外性
归功于 A. A. Albert（同一卷，65–73）。上面的证明遵循 Faraut–Korányi, *Analysis on Symmetric
Cones*（Oxford 1994）第 II–V 章的结构；Hurwitz 定理（Nachr. Ges. Wiss. Göttingen 1898, 309–316）的倍增证明取自
J. H. Conway, D. A. Smith, *On Quaternions and Octonions*（2003）以及 T. A. Springer, F. D. Veldkamp, *Octonions, Jordan
Algebras and Exceptional Groups*（2000），Frobenius 定理（1878）的证明则是 R. S. Palais 的证明，Amer. Math. Monthly
75 (1968) 366–368。延伸阅读（本课程中任何地方都未用到）：C. M. Glennie, Pacific J. Math. 16 (1966)
47–59，它通过多项式恒等式以不同的途径处理例外性；以及 McCrimmon 的 *A Taste of Jordan
Algebras*（2004），关于历史与一般理论。
