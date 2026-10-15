[English](../01_matrix_algebras.md) · [Español](../es/01_matrix_algebras.md) · [Français](../fr/01_matrix_algebras.md) · [Deutsch](../de/01_matrix_algebras.md) · **中文** · [日本語](../ja/01_matrix_algebras.md) · [Русский](../ru/01_matrix_algebras.md)

# 1. 矩阵代数与函数演算

本章固定记号，并证明后续各章都要用到的关于 $`M_n(\mathbb C)`$ 的事实：谱定理、函数演算及其导数、
极分解以及由奇异值定义的范数。这些结果的库代码和数值检验见[教程 1](../../tutorial/zh/01_matrix_algebras.md)。

**约定。** 内积 $`\langle x,y\rangle=\sum_i\bar x_iy_i`$ 关于第二个变元是线性的。$`A^*`$ 是共轭
转置，因此 $`\langle x,Ay\rangle=\langle A^*x,y\rangle`$。$`e_1,\dots,e_n`$ 是标准基，$`1`$ 是单位
矩阵，$`\mathrm{Tr}`$ 是通常的迹。

<a id="s-1-1"></a>

## 1.1 代数 $`M_n(\mathbb C)`$

$`M_n(\mathbb C)`$ 是复 $`n\times n`$ 矩阵构成的代数，带有对合
$`A\mapsto A^*`$ 和算子范数 $`\|A\|=\sup_{\|x\|=1}\|Ax\|`$。称元素 $`A`$ 为

| 名称 | 定义 |
|---|---|
| 正规的 | $`A^*A=AA^*`$ |
| 自伴的 | $`A=A^*`$ |
| 正的，$`A\ge0`$ | 对所有 $`x`$，$`\langle x,Ax\rangle\ge0`$ |
| 酉的 | $`A^*A=AA^*=1`$ |
| 投影 | $`A=A^*=A^2`$ |

若 $`B-A\ge0`$，记 $`A\le B`$；若 $`A\ge0`$ 且 $`A`$ 可逆，记 $`A>0`$。

*库中：* 具有指定谱的自伴或正算子，`TypeIAlgebra.operator_from_eigenvalues`——[教程，1.5 节](../../tutorial/zh/01_matrix_algebras.md#s-1-5)；算子范数，`operator_norm`——[教程，1.6 节](../../tutorial/zh/01_matrix_algebras.md#s-1-6)。

<a id="r-1-1"></a>

**引理 1.1（C\*-恒等式）。** $`\|A^*\|=\|A\|`$ 且 $`\|A^*A\|=\|A\|^2`$。

*证明。* $`\|Ax\|=\sup_{\|y\|=1}|\langle y,Ax\rangle|=\sup_{\|y\|=1}|\langle A^*y,x\rangle|`$，
所以 $`\|A\|=\sup_{\|x\|=\|y\|=1}|\langle A^*y,x\rangle|=\|A^*\|`$。其次，当 $`\|x\|=1`$ 时
$`\|Ax\|^2=\langle x,A^*Ax\rangle\le\|A^*A\|`$，因此 $`\|A\|^2\le\|A^*A\|\le\|A^*\|\,\|A\|=\|A\|^2`$。$`\square`$

> **物理意义**（物理解释，而非本课程的定理）。在量子力学中，具有 $`n`$ 个可区分构型的系统由 $`\mathbb C^n`$ 描述。自伴矩阵是
> *可观测量*（能量、自旋、格点上的位置）；它们的特征值是测量的可能结果。投影对应于是/否
> 问题。酉矩阵是对称变换和时间演化，对哈密顿量 $`H`$ 而言即 $`e^{-itH}`$。迹为一的正矩阵是
> *态*（密度矩阵，第 3 章）。

<a id="s-1-2"></a>

## 1.2 谱定理

<a id="r-1-2"></a>

**引理 1.2（Schur）。** 对每个 $`A\in M_n(\mathbb C)`$，存在酉矩阵 $`U`$ 使得
$`U^*AU`$ 为上三角矩阵。

*证明。* 对 $`n`$ 归纳；$`n=1`$ 是平凡的。$`A`$ 的特征多项式有一个根 $`\lambda`$（代数基本定理），
所以对某个单位向量 $`v`$ 有 $`Av=\lambda v`$。将 $`v`$ 扩充为规范正交基 $`v,w_2,\dots,w_n`$，并令
$`V`$ 为以这些向量为列的酉矩阵。则
$`V^*AV=\begin{pmatrix}\lambda&*\\0&A'\end{pmatrix}`$，其中 $`A'\in M_{n-1}(\mathbb C)`$。
由归纳假设，存在酉矩阵 $`U'`$ 使 $`U'^*A'U'`$ 为上三角矩阵，于是
$`U=V\,(1\oplus U')`$ 满足要求。$`\square`$

<a id="r-1-3"></a>

**定理 1.3（谱定理）。** $`A`$ 是正规的当且仅当对某个酉矩阵 $`U`$ 有
$`A=U\,\mathrm{diag}(\lambda_1,\dots,\lambda_n)\,U^*`$。这些
$`\lambda_i`$ 是 $`A`$ 的计重数的特征值。

*证明。* 若 $`A=U\Lambda U^*`$，$`\Lambda`$ 为对角矩阵，则
$`A^*A=U\bar\Lambda\Lambda U^*=U\Lambda\bar\Lambda U^*=AA^*`$。反之，设 $`A`$
正规，$`T=U^*AU`$ 为上三角矩阵（引理 1.2）；$`T`$ 是正规的。比较
$`(1,1)`$ 元：$`(T^*T)_{11}=|t_{11}|^2`$，$`(TT^*)_{11}=\sum_j|t_{1j}|^2`$，所以当 $`j>1`$ 时
$`t_{1j}=0`$。于是 $`T=t_{11}\oplus T'`$，其中 $`T'`$ 正规且为上三角矩阵，
由归纳法知 $`T`$ 是对角矩阵。$`U^*AU`$ 的对角线由
$`\det(A-\lambda)=\det(T-\lambda)`$ 的根组成。$`\square`$

*库中：* 把定理 1.3 反过来读，即具有指定特征值、且 $`U`$ 服从 Haar 分布的 $`A=U\Lambda U^*`$，就是 `TypeIAlgebra.operator_from_eigenvalues`——[教程，1.5 节](../../tutorial/zh/01_matrix_algebras.md#s-1-5)。

<a id="r-1-4"></a>

**推论 1.4.** 设 $`A`$ 正规，特征值为 $`\lambda_i`$。则

1. $`\|A\|=\max_i|\lambda_i|`$；
2. $`A`$ 自伴当且仅当所有 $`\lambda_i\in\mathbb R`$；$`A`$ 为酉矩阵当且仅当所有 $`|\lambda_i|=1`$；
   $`A`$ 为投影当且仅当所有 $`\lambda_i\in\{0,1\}`$；
3. $`A\ge0`$ 当且仅当 $`A`$ 自伴且所有 $`\lambda_i\ge0`$。此外，每个正矩阵都是自伴的，因而是正规的，所以 (3)
   在*所有*矩阵中刻画了正性，并且下面的函数演算适用于每个 $`A\ge0`$。

*证明。* 令 $`y=U^*x`$，则 $`\|Ax\|^2=\sum_i|\lambda_i|^2|y_i|^2`$ 且 $`\|y\|=\|x\|`$，
由此得 (1)，当 $`x=Ue_i`$ 且 $`|\lambda_i|`$ 最大时等号成立。在 (2) 中，$`A^*=U\bar\Lambda U^*`$，
$`A^*A=U|\Lambda|^2U^*`$，$`A^2=U\Lambda^2U^*`$，这些条件成立当且仅当它们对 $`\Lambda`$
逐元成立。对于 (3)，首先，每个正矩阵 $`A`$ 都是自伴的，因为对所有 $`x`$，
$`\langle x,(A-A^*)x\rangle=\langle x,Ax\rangle-\overline{\langle x,Ax\rangle}=0`$，
而对所有 $`x`$ 满足 $`\langle x,Dx\rangle=0`$ 的矩阵 $`D`$ 必为零（极化恒等式：
$`4\langle x,Dy\rangle=\sum_{k=0}^3 i^{-k}\langle x+i^ky,D(x+i^ky)\rangle`$）。其次，
$`\langle x,Ax\rangle=\sum_i\lambda_i|y_i|^2`$，它对所有 $`x`$ 都 $`\ge0`$ 当且仅当所有 $`\lambda_i\ge0`$。$`\square`$

<a id="s-1-3"></a>

## 1.3 函数演算

设 $`A=U\Lambda U^*`$ 正规，$`f`$ 是 $`\mathrm{spec}\,A`$ 上的复值函数。定义
```math
f(A)=U\,\mathrm{diag}\big(f(\lambda_1),\dots,f(\lambda_n)\big)\,U^*.
```

<a id="r-1-5"></a>

**定理 1.5.** (a) $`f(A)`$ 与 $`U`$ 的选取无关。(b) 对多项式 $`p`$，$`p(A)`$ 就是通常意义下
$`A`$ 的多项式。(c) $`f\mapsto f(A)`$ 是线性且可乘的，$`\bar f(A)=f(A)^*`$，且 $`\mathrm{spec}\,f(A)=f(\mathrm{spec}\,A)`$。
(d) 若 $`B`$ 与 $`A`$ 可交换，则 $`B`$ 与 $`f(A)`$ 可交换。(e) 若 $`g`$ 定义在
$`f(\mathrm{spec}\,A)`$ 上，则 $`g(f(A))=(g\circ f)(A)`$。

*证明。* (b) $`A^k=U\Lambda^kU^*`$，由线性得 $`p(A)=Up(\Lambda)U^*`$。(a) 谱是有限集，
所以由 Lagrange 插值可得多项式 $`p`$，使得在 $`\mathrm{spec}\,A`$ 上 $`p=f`$：若 $`\mu_1,\dots,\mu_r`$ 是互不相同的特征值，取
$`p(t)=\sum_if(\mu_i)\prod_{j\neq i}(t-\mu_j)/(\mu_i-\mu_j)`$。于是对*每个*
对角化 $`A`$ 的酉矩阵 $`U`$ 都有 $`U\,f(\Lambda)\,U^*=U\,p(\Lambda)\,U^*=p(A)`$，而 $`p(A)`$ 的定义
不依赖于 $`U`$。(c) 对对角矩阵 $`f(\Lambda)`$ 成立，且在用 $`U`$ 共轭时保持；$`Uf(\Lambda)U^*`$
的谱就是 $`f(\Lambda)`$ 的对角线。(d) $`B`$ 与 $`A`$ 的每个多项式可交换，特别地与
$`p(A)=f(A)`$ 可交换。(e) $`f(A)=Uf(\Lambda)U^*`$ 是正规的，且由同一个 $`U`$ 对角化，所以由 (a)
得 $`g(f(A))=U\,g(f(\Lambda))\,U^*=(g\circ f)(A)`$。$`\square`$

<a id="r-1-6"></a>

**定理 1.6（平方根）。** 对 $`A\ge0`$，恰有一个 $`B\ge0`$ 满足
$`B^2=A`$，即由函数演算定义的 $`B=\sqrt A`$。

*证明。* 由定理 1.5(c) 和推论 1.4(3)，$`\sqrt A\ge0`$ 且 $`(\sqrt A)^2=A`$。
设 $`B\ge0`$，$`B^2=A`$，$`B=V\,\mathrm{diag}(\mu_j)\,V^*`$，$`\mu_j\ge0`$。则
$`A=V\,\mathrm{diag}(\mu_j^2)\,V^*`$，所以 $`\mathrm{spec}\,A=\{\mu_j^2\}`$。设 $`p`$ 是
在 $`\mathrm{spec}\,A`$ 上满足 $`p(t)=\sqrt t`$ 的多项式。则
$`p(A)=\sqrt A`$，且 $`p(A)=p(B^2)=V\,\mathrm{diag}(p(\mu_j^2))\,V^*=V\,\mathrm{diag}(\mu_j)\,V^*=B`$。$`\square`$

*计算：* 用 `apply_function`、`sqrt`、`log`、`expm` 和 `power` 检验定理 1.5 和定理 1.6——[教程，1.5 节](../../tutorial/zh/01_matrix_algebras.md#s-1-5)。

<a id="r-1-7"></a>

**推论 1.7.** $`A\ge0`$ 当且仅当对某个 $`T`$ 有 $`A=T^*T`$。（若 $`A\ge0`$，取 $`T=\sqrt A`$；
反之 $`\langle x,T^*Tx\rangle=\|Tx\|^2\ge0`$。）

**导数。** 梯度方法要对 $`f(A)`$ 关于 $`A`$ 求导。对自伴的 $`A`$，导数是与 $`f`$ 的差商矩阵的逐元素乘积。

<a id="r-1-11"></a>

**定理 1.11（Daleckii–Krein）。** 设 $`A=U\,\mathrm{diag}(\lambda_1,\dots,\lambda_n)\,U^*`$ 自伴，$`I`$ 是包含 $`\mathrm{spec}\,A`$
的开区间，$`f:I\to\mathbb C`$ 连续可微。设
```math
f^{[1]}(x,y)=\frac{f(x)-f(y)}{x-y}\quad(x\ne y),\qquad f^{[1]}(x,x)=f'(x)
```
为 $`f`$ 的**差商**（divided difference），$`L_{ij}=f^{[1]}(\lambda_i,\lambda_j)`$ 为 **Loewner 矩阵**。则对每个自伴的 $`E`$，
```math
\frac{d}{dt}\Big|_{t=0}f(A+tE)=U\big(L\circ U^*EU\big)U^*,
```
其中 $`(X\circ Y)_{ij}=X_{ij}Y_{ij}`$ 是逐元素（Schur）乘积。

*证明。* 由微积分基本定理，对所有 $`x,y\in I`$ 有 $`f^{[1]}(x,y)=\int_0^1f'(y+s(x-y))\,ds`$（$`x=y`$ 时两边都是
$`f'(x)`$），所以 $`f^{[1]}`$ 在 $`I\times I`$ 上连续。重排 $`U`$ 的列会以同样方式置换 $`L`$ 和 $`U^*EU`$ 的行与列，不改变右边，
因此可设 $`\lambda_1\le\dots\le\lambda_n`$。取 $`\delta>0`$ 使 $`[\lambda_1-\delta,\lambda_n+\delta]\subset I`$。对
$`0<|t|\,\|E\|<\delta`$，写 $`A+tE=V\,\mathrm{diag}(\mu_1,\dots,\mu_n)\,V^*`$，$`\mu_1\le\dots\le\mu_n`$；则
$`|\mu_j-\lambda_j|\le|t|\,\|E\|`$（定理 A.10），故所有 $`\mu_j\in I`$。令 $`W=U^*V`$（酉矩阵），$`\tilde E=U^*EU`$。由于
$`U^*f(A)=\mathrm{diag}(f(\lambda_i))\,U^*`$ 且 $`f(A+tE)V=V\,\mathrm{diag}(f(\mu_j))`$，
```math
\big(U^*(f(A+tE)-f(A))V\big)_{ij}=\big(f(\mu_j)-f(\lambda_i)\big)W_{ij}=f^{[1]}(\lambda_i,\mu_j)\,(\mu_j-\lambda_i)\,W_{ij},
```
对 $`f(x)=x`$ 作同样计算得 $`(\mu_j-\lambda_i)W_{ij}=t\,(\tilde EW)_{ij}`$。记 $`(L_t)_{ij}=f^{[1]}(\lambda_i,\mu_j)`$，这就是精确恒等式
```math
\frac1t\big(f(A+tE)-f(A)\big)=U\big((L_t\circ\tilde EW)\,W^*\big)U^*.
```
当 $`t\to0`$ 时，由于 $`f^{[1]}`$ 连续且 $`\mu_j\to\lambda_j`$，$`L_t\to L`$ 逐元素成立。令：若 $`\lambda_i=\lambda_j`$ 则
$`D_{ij}=W_{ij}`$，否则 $`D_{ij}=0`$；并令 $`R=W-D`$。若 $`\lambda_i\ne\lambda_j`$ 且 $`|t|\,\|E\|\le|\lambda_i-\lambda_j|/2`$，则
$`|\mu_j-\lambda_i|\ge|\lambda_i-\lambda_j|/2`$ 且 $`|(\tilde EW)_{ij}|\le\|\tilde EW\|=\|E\|`$，从而
$`|W_{ij}|\le2|t|\,\|E\|/|\lambda_i-\lambda_j|`$。因此 $`R\to0`$，且 $`DD^*=(W-R)(W-R)^*\to WW^*=1`$。由于只要 $`\lambda_k=\lambda_j`$
就有 $`L_{ij}=L_{ik}`$，而除非 $`\lambda_k=\lambda_j`$ 否则 $`D_{kj}=0`$，所以
$`(L\circ\tilde ED)_{ij}=\sum_kL_{ik}\tilde E_{ik}D_{kj}=((L\circ\tilde E)D)_{ij}`$。于是
```math
(L_t\circ\tilde EW)W^*-L\circ\tilde E=\big((L_t-L)\circ\tilde EW\big)W^*+(L\circ\tilde ER)W^*+(L\circ\tilde ED)R^*+(L\circ\tilde E)(DD^*-1),
```
而每一项都趋于 $`0`$，因为 $`W`$、$`D`$、$`\tilde E`$ 的元素分别以 $`1`$、$`1`$、$`\|E\|`$ 为界。$`\square`$

对实值 $`f`$，$`L`$ 是实矩阵，并且对每个 $`G`$ 有 $`\mathrm{Tr}\big(G^*U(L\circ U^*EU)U^*\big)=\mathrm{Tr}(\Gamma^*E)`$，其中
$`\Gamma=U(L\circ U^*GU)U^*`$（在基 $`U`$ 下把两边写成对 $`i,j`$ 的求和即可）：$`A\mapsto\mathrm{Re}\,\mathrm{Tr}(G^*f(A))`$ 沿方向
$`E`$ 的导数是 $`\mathrm{Re}\,\mathrm{Tr}(\Gamma^*E)`$。

*库中：* `apply_function` 借助定理 1.11 求导（通过 autograd）；Loewner 矩阵由 `functional.divided_differences` 计算。

**舍入。** 当 $`\lambda_i,\lambda_j`$ 相近时，商 $`(f(\lambda_i)-f(\lambda_j))/(\lambda_i-\lambda_j)`$ 因相消而损失精度，当
$`\lambda_i=\lambda_j`$ 时它没有定义。这时导数是更好的数据。

<a id="r-1-12"></a>

**命题 1.12（中点值）。** 设 $`a<b`$，$`h=b-a`$，$`f`$ 在 $`[a,b]`$ 上三次连续可微（允许取复值）。则
```math
f^{[1]}(a,b)-\frac{f'(a)+f'(b)}2=-\frac{h^2}2\int_0^1s(1-s)\,f'''(a+sh)\,ds,
\qquad
\Big|f^{[1]}(a,b)-\frac{f'(a)+f'(b)}2\Big|\le\frac{h^2}{12}\max_{[a,b]}|f'''|,
```
对 $`3`$ 次多项式等号成立。

*证明。* 令 $`g(s)=f'(a+sh)`$，使 $`f^{[1]}(a,b)=\frac1h\int_a^bf'=\int_0^1g`$；令 $`\varphi(s)=s(1-s)/2`$，使
$`\varphi(0)=\varphi(1)=0`$，$`\varphi'(0)=\frac12=-\varphi'(1)`$，$`\varphi''=-1`$。两次分部积分得
$`\int_0^1\varphi g''=-\int_0^1\varphi'g'=-[\varphi'g]_0^1+\int_0^1\varphi''g=\frac{g(0)+g(1)}2-\int_0^1g`$。由于
$`g''(s)=h^2f'''(a+sh)`$，这就是所述恒等式；$`\int_0^1\varphi=\frac1{12}`$ 给出界，而对三次多项式 $`f'''`$ 是常数。$`\square`$

<a id="r-1-13"></a>

**命题 1.13（商还是中点值）。** 设浮点运算的单位舍入误差为 $`u\le1/20`$：对浮点数 $`x,y`$ 和
$`\circ\in\{+,-,\times,/\}`$，计算结果为 $`(x\circ y)(1+\theta)`$，$`|\theta|\le u`$（无上溢和下溢）。设 $`a\ne b`$ 为实浮点数，
$`h=|a-b|`$，$`x=f^{[1]}(a,b)`$，$`\hat f_a,\hat f_b`$ 为计算值，满足 $`|\hat f_a-f(a)|\le\delta_a`$、$`|\hat f_b-f(b)|\le\delta_b`$。

(a) $`\hat f_a-\hat f_b`$ 与 $`a-b`$ 的计算商 $`\hat q`$ 满足 $`|\hat q-x|\le R:=(\delta_a+\delta_b)/h+5u|\hat q|`$。

(b) 设 $`|\hat m-x|\le T`$（例如 $`\hat m`$ 为计算出的中点值，$`T`$ 为命题 1.12 的界加上其舍入误差）。若
$`|\hat q-\hat m|\le2R`$ 则令 $`y=\hat m`$，否则令 $`y=\hat q`$。则 $`|y-x|\le\min(T,3R)`$。

(c) 若 $`\delta_a+\delta_b\le\Delta`$ 且 $`T=Mh^2+T_0`$，则对每个 $`h`$ 有 $`|y-x|\le(3\Delta)^{2/3}M^{1/3}+T_0+15u|\hat q|`$。两个误差项
$`3\Delta/h`$ 与 $`Mh^2`$ 在 $`h_*=(3\Delta/M)^{1/3}`$ 处相等。

*证明。* (a) 令 $`\tilde q=(\hat f_a-\hat f_b)/(a-b)`$。三次运算给出
$`\hat q=\tilde q\,(1+\theta_1)(1+\theta_3)/(1+\theta_2)=\tilde q\,(1+\theta)`$，$`|\theta_i|\le u`$，故
$`-3u\le\frac{(1-u)^2}{1+u}-1\le\theta\le\frac{(1+u)^2}{1-u}-1=\frac{3u+u^2}{1-u}\le4u`$。于是
$`|\hat q-\tilde q|=|\theta|\,|\tilde q|\le\frac{4u}{1-4u}|\hat q|\le5u|\hat q|`$，且
$`|\tilde q-x|=|(\hat f_a-f(a))-(\hat f_b-f(b))|/h\le(\delta_a+\delta_b)/h`$。(b) 若 $`y=\hat m`$，则 $`|y-x|\le T`$ 且
$`|y-x|\le|\hat m-\hat q|+|\hat q-x|\le3R`$。若 $`y=\hat q`$，则 $`|y-x|\le R`$，且
$`T\ge|\hat m-x|\ge|\hat m-\hat q|-|\hat q-x|>2R-R=R`$。(c) $`\min(T,3R)\le\min(Mh^2,3\Delta/h)+T_0+15u|\hat q|`$；当
$`h\ge h_*`$ 时 $`3\Delta/h\le3\Delta/h_*`$，当 $`h<h_*`$ 时 $`Mh^2<Mh_*^2`$，而 $`3\Delta/h_*=Mh_*^2=(3\Delta)^{2/3}M^{1/3}`$。$`\square`$

该证明对复值 $`f`$ 逐字成立（实部与虚部分别舍入，$`a-b`$ 为实数；此时 $`\theta`$ 为复数，界相同）。规则 (b) 需要由数据
算出的 $`R`$，但不需要 $`T`$：不必知道 $`f'''`$，其误差至多是两个界中较小者的三倍。固定阈值（当且仅当 $`h\le h_*`$ 时取
中点值）则需要 $`M=\max|f'''|/12`$。本库早期版本的阈值 $`h\le\varepsilon^{1/3}\max(|a|,|b|)`$（$`\varepsilon=2u`$）就属这一类，
相当于 $`M`$ 取 $`\Delta/\max(|a|,|b|)^3`$ 的量级：对导数随 $`|x|`$ 缩放的 $`x^p`$ 或 $`\log x`$ 是对的，对在与 $`|x|`$ 无关的尺度上
变化的函数则是错的。对 $`a\approx c=10^4`$ 处的 $`f(x)=e^{x-c}`$，它在 $`h\approx0.06`$ 以内都取中点值，相对误差可达
$`h^2/12\approx3\cdot10^{-4}`$（在那里命题 1.12 除一个因子 $`e^{h}`$ 外是等式）；对 $`e^{-i\tau x}`$ 和带有大常数平移的哈密顿量
也是如此。

*库中：* `functional.divided_differences` 对每一对本征值应用规则 (b)，取 $`\delta_a=\varepsilon\,(|\hat f_a|+|a\hat f'_a|)`$，并以
$`4\varepsilon|\hat q|`$ 代替 $`5u|\hat q|`$。$`\delta_a`$ 的选取是关于 $`f`$ 求值的假设，而非定理：在至多扰动 $`\varepsilon|a|`$ 的
自变量处，相对误差至多 $`\varepsilon`$（一阶近似；第二项涵盖诸如 $`e^{-i\tau x}`$ 这类要对 $`\tau x`$ 舍入的函数）。同一规则也
用于正规矩阵的复本征值，这不在证明范围之内。对相等的本征值以及本征值求解器无法分辨的本征值对
$`h\le n\varepsilon\max_l|\lambda_l|`$，也取中点值：由命题 1.12，其误差至多为 $`(n\varepsilon\max_l|\lambda_l|)^2\max|f'''|/12`$，而对于
不连续的 $`f`$（如示性函数），商是没有意义的。

*计算：* 在简并谱处将定理 1.11 与差商比较，并对 $`c=10^4`$ 的 $`e^{x-c}`$ 检验命题 1.12 和 1.13——[教程，1.5 节](../../tutorial/zh/01_matrix_algebras.md#s-1-5)。

<a id="s-1-4"></a>

## 1.4 极分解与奇异值

对任意 $`A`$，矩阵 $`A^*A`$ 是正的（推论 1.7）；它的平方根
$`|A|=(A^*A)^{1/2}`$ 称为 $`A`$ 的**绝对值**。

<a id="r-1-8"></a>

**定理 1.8（极分解）。** 每个 $`A\in M_n(\mathbb C)`$ 都可以写成
$`A=W|A|`$，其中 $`W`$ 为酉矩阵。

*证明。* 对每个 $`x`$，$`\||A|x\|^2=\langle x,|A|^2x\rangle=\langle x,A^*Ax\rangle=\|Ax\|^2`$。
因此 $`W_0:|A|x\mapsto Ax`$ 是从 $`|A|`$ 的值域 $`R_1`$ 到 $`A`$ 的值域 $`R_2`$ 上的良定义的线性
等距映射：若 $`|A|x=|A|x'`$，则
$`\|A(x-x')\|=\||A|(x-x')\|=0`$。特别地 $`\dim R_1=\dim R_2`$，所以二者的正交补也有相同
维数；令 $`W_1:R_1^\perp\to R_2^\perp`$ 为它们之间的任一等距映射（把一组规范正交基映为一组
规范正交基）。令
$`W=W_0\oplus W_1`$，即对 $`u_1\in R_1`$、$`u_2\in R_1^\perp`$ 有 $`W(u_1+u_2)=W_0u_1+W_1u_2`$。由于
$`W_0u_1\in R_2\perp R_2^\perp\ni W_1u_2`$，由勾股定理，对所有 $`u`$ 有 $`\|Wu\|^2=\|u_1\|^2+\|u_2\|^2=\|u\|^2`$。
因此对所有 $`u`$ 有 $`\langle u,(W^*W-1)u\rangle=0`$，由极化恒等式（推论 1.4(3) 的证明）得 $`W^*W=1`$；于是
$`W`$ 可逆且 $`W^{-1}=W^*`$，所以也有 $`WW^*=1`$：$`W`$ 是酉的，且 $`W|A|x=W_0|A|x=Ax`$。$`\square`$

*计算：* 在一个 Ginibre 矩阵上检验定理 1.8，其中 $`|A|`$ 由 `abs` 计算——[教程，1.6 节](../../tutorial/zh/01_matrix_algebras.md#s-1-6)。

$`|A|`$ 的特征值 $`s_1\ge\dots\ge s_n\ge0`$ 称为 $`A`$ 的**奇异值**。全书使用以下三种范数：
```math
\|A\|=s_1,\qquad \|A\|_2=\Big(\sum_is_i^2\Big)^{1/2}=(\mathrm{Tr}\,A^*A)^{1/2},
\qquad\|A\|_1=\sum_is_i=\mathrm{Tr}|A|.
```
第一个等式成立是因为对每个 $`x`$ 有 $`\|Ax\|=\||A|x\|`$，所以由推论 1.4(1)，
$`\|A\|=\||A|\|=s_1`$；第二个等式成立是因为
$`\mathrm{Tr}\,A^*A=\mathrm{Tr}|A|^2=\sum_is_i^2`$。

<a id="r-1-9"></a>

**定理 1.9（迹范数的对偶性）。**
$`\|A\|_1=\max_{U}|\mathrm{Tr}(AU)|`$（最大值取遍酉矩阵 $`U`$），且最大值在
$`U=W^*`$ 处取得，其中 $`A=W|A|`$。特别地 $`|\mathrm{Tr}\,A|\le\|A\|_1`$。

*证明。* $`\mathrm{Tr}(AW^*)=\mathrm{Tr}(W|A|W^*)=\mathrm{Tr}|A|`$。对任意酉矩阵 $`U`$，
由迹的循环性以及关于 $`\langle X,Y\rangle_{HS}=\mathrm{Tr}(X^*Y)`$ 的 Cauchy–Schwarz 不等式，
```math
|\mathrm{Tr}(AU)|=|\mathrm{Tr}(|A|^{1/2}\,UW\,|A|^{1/2})|
=|\langle|A|^{1/2},\,UW|A|^{1/2}\rangle_{HS}|
\le\||A|^{1/2}\|_2\,\|UW|A|^{1/2}\|_2=\mathrm{Tr}|A|,
```
这是因为对酉矩阵 $`UW`$ 有 $`\|UW X\|_2=\|X\|_2`$，且 $`\||A|^{1/2}\|_2^2=\mathrm{Tr}|A|`$。$`\square`$

特别地，上述三个量都是范数。$`\|\cdot\|`$ 是算子范数；
$`\|\cdot\|_2`$ 是内积 $`\langle X,Y\rangle_{HS}`$ 诱导的范数；对于
$`\|\cdot\|_1`$，由定理 1.9 得
$`\|A+B\|_1=\max_U|\mathrm{Tr}(AU)+\mathrm{Tr}(BU)|\le\|A\|_1+\|B\|_1`$ 且
$`\|cA\|_1=\max_U|c\mathrm{Tr}(AU)|=|c|\,\|A\|_1`$，而
$`\|A\|_1=0`$ 迫使 $`|A|=0`$，从而 $`A=W|A|=0`$。对 $`P\ge0`$，由定理 1.6 有 $`|P|=(P^2)^{1/2}=P`$，
所以 $`\|P\|_1=\mathrm{Tr}\,P`$。

*计算：* 三种范数（`operator_norm`、`frobenius_norm`、`trace_norm`）与定理 1.9——[教程，1.6 节](../../tutorial/zh/01_matrix_algebras.md#s-1-6)。

<a id="s-1-5"></a>

## 1.5 库中的算子

*库中：* 代数 $`M_n(\mathbb C)`$ 即 `TypeIAlgebra([n], [n])`，其算子如何存储，以及函数演算 `apply_function`——[教程，1.5 节](../../tutorial/zh/01_matrix_algebras.md#s-1-5)。

<a id="s-1-6"></a>

## 1.6 范数与迹范数对偶

*计算：* 在一个 Ginibre 矩阵上检验1.4 节的范数、极分解和对偶性——[教程，1.6 节](../../tutorial/zh/01_matrix_algebras.md#s-1-6)。

<a id="s-1-7"></a>

## 1.7 批量：同时处理许多算子

作为一个可以用许多随机算子同时估计的量的例子，我们计算**高斯酉系综**（GUE）的一个精确矩。

设 $`G`$ 的各元素相互独立，且实部与虚部相互独立、服从
$`N(0,\tfrac12)`$，从而 $`\mathbb E|G_{ij}|^2=1`$，并令 $`H=(G+G^*)/\sqrt2`$。

<a id="r-1-10"></a>

**命题 1.10.** $`\mathbb E\,\mathrm{Tr}\,H^4=2n^3+n`$。因此对
$`X=H/\sqrt n`$，$`\ \mathbb E\,\tfrac1n\mathrm{Tr}\,X^4=2+\tfrac1{n^2}`$。

*证明。* *Isserlis（Wick）公式*（定理 A.24，$`m=4`$ 的情形）。$`G`$ 的各元素的实部与虚部是相互独立的
中心化高斯变量，所以它们是同一个中心化高斯向量的坐标，而 $`H`$ 与 $`\bar H`$ 的元素是它们的复线性组合。
于是定理 A.24 对任意四个这样的变量 $`Z_1,\dots,Z_4`$（允许重复）给出
```math
\mathbb E[Z_1Z_2Z_3Z_4]=\mathbb E[Z_1Z_2]\,\mathbb E[Z_3Z_4]+\mathbb E[Z_1Z_3]\,\mathbb E[Z_2Z_4]+\mathbb E[Z_1Z_4]\,\mathbb E[Z_2Z_3],
```
即对 $`\{1,2,3,4\}`$ 的三种配对求和，每种配对贡献各对期望的乘积（不取复共轭）。附录 A 中的证明
把矩母函数的幂级数与高斯矩相比较，不需要在期望号下求导。

*矩的计算。* 配对期望为
$`\mathbb E[H_{ij}H_{kl}]=\tfrac12\mathbb E[(G_{ij}+\bar G_{ji})(G_{kl}+\bar G_{lk})]=\delta_{il}\delta_{jk}`$，
因为 $`\mathbb E[G_{ij}G_{kl}]=0`$ 且 $`\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}`$（命题 A.16(a)；
当 $`(i,j)\neq(k,l)`$ 时这两个元素相互独立且中心化，而对 $`G_{ij}=x+iy`$ 有
$`\mathbb E[G_{ij}^2]=\mathbb E[x^2-y^2]+2i\,\mathbb E[xy]=0`$ 和 $`\mathbb E|G_{ij}|^2=\tfrac12+\tfrac12=1`$）。
在 $`\mathrm{Tr}\,H^4=\sum_{ijkl}H_{ij}H_{jk}H_{kl}H_{li}`$ 中，各配对给出：
$`\{12\}\{34\}`$：$`\delta_{ik}\delta_{ki}`$，求和得 $`n^3`$（自由指标为 $`i=k`$、$`j`$、$`l`$）；
$`\{14\}\{23\}`$：$`\delta_{jl}\delta_{jl}`$，求和得 $`n^3`$；
$`\{13\}\{24\}`$：$`\delta_{il}\delta_{jk}\delta_{ji}\delta_{kl}`$，迫使 $`i=j=k=l`$，求和得 $`n`$。
除以 $`n\cdot n^2`$ 即得第二个结论。$`\square`$

两种*非交叉*配对 $`\{12\}\{34\}`$ 和 $`\{14\}\{23\}`$（这些对可以画成点 $`1,2,3,4`$ 上方互不相交的弧）
对每个 $`n`$ 都给出首项 $`2`$；交叉配对 $`\{13\}\{24\}`$ 给出 $`1/n^2`$ 修正项。

*计算：* 用蒙特卡罗方法在 $`200\,000`$ 个样本的批量上检验命题 1.10——[教程，1.7 节](../../tutorial/zh/01_matrix_algebras.md#s-1-7)。

<a id="s-1-8"></a>

## 1.8 可与之比较的工具

*库中：* 可供比较的参考实现——[教程，1.8 节](../../tutorial/zh/01_matrix_algebras.md#s-1-8)。

<a id="x-1"></a>

## 习题

1. 对正规的 $`A`$ 以及 $`\mathrm{spec}\,A`$ 上的函数 $`f,g`$，证明 $`(fg)(A)=f(A)g(A)`$。
2. 利用定理 1.9 证明 $`\|A\|_1=\|A^*\|_1`$ 且 $`\|AB\|_1\le\|A\|\,\|B\|_1`$。
3. $`A,B>0`$ 的**几何平均**定义为 $`A\#B=A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}`$。
   证明 $`A\#B`$ 是方程 $`XA^{-1}X=B`$ 的唯一正解 $`X`$，并由此推出
   $`A\#B=B\#A`$。

<details>
<summary>解答</summary>

1. 由定理 1.5(c) 即得。

2. $`|\mathrm{Tr}(A^*U)|=|\overline{\mathrm{Tr}(U^*A)}|=|\mathrm{Tr}(AU^*)|`$，而
   $`U\mapsto U^*`$ 是酉群到自身的双射，所以定理 1.9 中的两个最大值相等。对第二个结论，设
   $`B=V|B|`$。对酉矩阵 $`U`$，
   $`|\mathrm{Tr}(ABU)|=|\mathrm{Tr}(|B|^{1/2}UAV|B|^{1/2})|    \le\||B|^{1/2}\|_2\,\|UAV|B|^{1/2}\|_2\le\|A\|\mathrm{Tr}|B|`$，这里用到了
   $`\|YX\|_2\le\|Y\|\,\|X\|_2`$（对 $`X`$ 的每一列应用 $`\|Yx\|\le\|Y\|\,\|x\|`$）以及对酉矩阵
   $`U,V`$ 有 $`\|UAV\|=\|A\|`$（酉矩阵保持向量的长度）。

3. 记 $`T=A^{-1/2}XA^{-1/2}`$。则 $`XA^{-1}X=B`$ 当且仅当 $`T^2=A^{-1/2}BA^{-1/2}`$，并且
   $`X\ge0`$ 当且仅当 $`T\ge0`$，因为 $`\langle x,Xx\rangle=\langle A^{1/2}x,TA^{1/2}x\rangle`$ 且
   $`A^{1/2}`$ 可逆（这里 $`A^{\pm1/2}`$ 由 $`\mathrm{spec}\,A\subset(0,\infty)`$ 上的函数 $`t^{\pm1/2}`$ 通过函数演算
   给出，所以由定理 1.5(c) 有 $`A^{-1/2}A^{1/2}=1`$）。此外
   $`A^{-1/2}BA^{-1/2}\ge0`$，因为 $`\langle x,A^{-1/2}BA^{-1/2}x\rangle=\langle A^{-1/2}x,BA^{-1/2}x\rangle\ge0`$。
   由定理 1.6，$`T=(A^{-1/2}BA^{-1/2})^{1/2}`$ 是唯一的正解，即
   $`X=A\#B`$。方程 $`XA^{-1}X=B`$ 的每个解 $`X`$ 都可逆，
   因为 $`\det(X)^2\det(A)^{-1}=\det B\neq0`$。所以 $`XA^{-1}X=B`$ 等价于
   $`XB^{-1}X=A`$（由 $`A^{-1}=X^{-1}BX^{-1}`$ 两边取逆），再由唯一性得
   $`A\#B=B\#A`$。

</details>

*计算：* 习题 1 和习题 3 的数值检验——[教程，习题](../../tutorial/zh/01_matrix_algebras.md#x-1)。

**下一章：** [第 2 章](02_von_neumann_algebras.md)把 $`M_n`$ 换成矩阵代数的直和，并证明每个
有限维冯·诺依曼代数都具有这种形式。
