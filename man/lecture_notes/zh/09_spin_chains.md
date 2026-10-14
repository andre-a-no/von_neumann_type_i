[English](../09_spin_chains.md) · [Español](../es/09_spin_chains.md) · [Français](../fr/09_spin_chains.md) · [Deutsch](../de/09_spin_chains.md) · **中文** · [日本語](../ja/09_spin_chains.md) · [Русский](../ru/09_spin_chains.md)

# 9. 自旋链

由 $`L`$ 个自旋 $`\frac12`$ 组成的链位于 $`(\mathbb C^2)^{\otimes L}`$ 上，而维数为 $`2^L`$ 的稠密线性代数在
$`L\approx14`$ 以上就变得代价高昂（一个稠密矩阵有 $`4^L`$ 个元素）。对称性能推高这一极限：磁化把空间分成 $`L+1`$ 个扇区（第 4 章），平移把每个扇区再分一次，
而 $`SU(2)`$ 又一次约化这些块。本章为 XXZ 链构造这些扇区，证明它们关于基态说明了什么，用自由费米子和 Bethe
拟设检验本库，并超越稠密矩阵（Krylov 方法）和纯净的链（无序）。

**你将学到**为什么链的代数是 $`\bigoplus_NM_{\binom LN}`$，为什么这些代数构成的塔以 Pascal 三角形为 Bratteli 图、
以 $`4\cos^2\frac\pi{2L+4}`$ 为指标，如何计数动量扇区，Marshall 符号规则及其对基态动量与自旋的推论
（Lieb–Mattis），XX 链的 Jordan–Wigner 解，双磁振子 Bethe 方程，为什么 Lanczos 方法单调收敛，以及如何测量
能级统计和自由费米子的纠缠。本章的数值检验见
[教程 9](../../tutorial/zh/09_spin_chains.md)。

在本章中，格点为 $`0,\dots,L-1`$，$`(\mathbb C^2)^{\otimes L}`$ 的基由构型
$`|b\rangle=|b_0\dots b_{L-1}\rangle`$ 组成，自旋向上时 $`b_i=1`$（本库的约定：标号的第 $`i`$ 位
对应格点 $`i`$），$`s^a_i`$（$`a=x,y,z`$）是格点 $`i`$ 上的自旋矩阵 $`\frac12\sigma_a`$，
$`s^\pm_i=s^x_i\pm is^y_i`$，$`n_i=s^z_i+\frac12`$，$`S^a=\sum_is^a_i`$，而 $`N=S^z+\frac L2`$ 是向上自旋的
个数。链是**开链**（键 $`(i,i+1)`$，$`i<L-1`$）或**周期链**（还有键 $`(L-1,0)`$；此时
$`L\ge3`$，指标取模 $`L`$）。$`J>0`$（反铁磁），$`\hbar=1`$，$`\log`$ 为自然对数。

<a id="s-9-1"></a>

## 9.1 XXZ 链及其 $`U(1)`$ 扇区

各向异性为 $`\Delta\in\mathbb R`$、外场为 $`h_i\in\mathbb R`$ 的 **XXZ 链**是
```math
H=J\sum_{(i,j)}\Big[\tfrac12\big(s^+_is^-_j+s^-_is^+_j\big)+\Delta\,s^z_is^z_j\Big]+\sum_ih_is^z_i ,
```
求和遍及所有键；$`\frac12(s^+_is^-_j+s^-_is^+_j)=s^x_is^x_j+s^y_is^y_j`$ 使一个向上自旋沿键移动。
$`\Delta=1`$ 时它是 Heisenberg 链 $`J\sum\mathbf s_i\cdot\mathbf s_j`$，$`\Delta=0`$ 时是 XX 链。

<a id="r-9-1"></a>

**命题 9.1.** $`[H,S^z]=0`$。因此 $`H\in U(1)'\cong\bigoplus_{N=0}^LM_{\binom LN}(\mathbb C)`$，其中
块 $`N`$ 作用在由含 $`N`$ 个向上自旋的构型张成的空间上。

*证明。* 对角项与对角的 $`S^z`$ 可交换。对 $`i\neq j`$，$`s^+_is^-_j`$ 与 $`s^z_l`$（$`l\ne i,j`$）可交换，并且由
$`[s^z,s^\pm]=\pm s^\pm`$（4.5 节），$`[s^z_i+s^z_j,s^+_is^-_j]=[s^z_i,s^+_i]s^-_j+s^+_i[s^z_j,s^-_j]=s^+_is^-_j-s^+_is^-_j=0`$；
对 $`s^-_is^+_j`$ 取伴随即可。因此 $`H`$ 与 $`N=S^z+\frac L2`$ 可交换，后者具有整数特征值 $`0,\dots,L`$，其特征空间由
含 $`N`$ 个向上自旋的构型张成，维数为 $`\binom LN`$；其余部分是取 $`Q=N`$ 的命题 4.10(2)。$`\square`$

*计算：* 命题 9.1，以及 `SpinChain.xxz` 由各扇区的构型（`chain.labels[N]`）构造的块 $`H_N`$，与完整的 $`2^L\times2^L`$ 矩阵对照——[教程，9.1 节](../../tutorial/zh/09_spin_chains.md#s-9-1)。

**链的塔。** 令 $`\mathcal A_L=\bigoplus_NM_{\binom LN}`$ 作用在 $`L`$ 个格点上。添加一个格点，
$`x\mapsto x\otimes1_2`$，把 $`\mathcal A_L`$ 映入 $`\mathcal A_{L+1}`$，因为 $`x\otimes1`$ 与
$`S^z\otimes1+1\otimes s^z`$ 可交换。这就给出一个塔 $`\mathcal A_0\subset\mathcal A_1\subset\mathcal A_2\subset\cdots`$。

<a id="r-9-2"></a>

**命题 9.2.** 包含 $`\mathcal A_L\subset\mathcal A_{L+1}`$ 的包含矩阵是 $`(L+2)\times(L+1)`$
矩阵 $`\Lambda_L`$，其中 $`\lambda_{N,N}=\lambda_{N,N-1}=1`$，其余元素均为 $`0`$；块的
大小服从 Pascal 法则 $`\binom{L+1}N=\binom LN+\binom L{N-1}`$，塔的 Bratteli 图是
Pascal 三角形，而 $`\Lambda_L`$ 的二部图是一条有 $`2L+3`$ 个顶点的路径。

*证明。* $`L+1`$ 个格点上含 $`N`$ 个向上自旋的构型，一类是前 $`L`$ 个格点上有 $`N`$ 个向上自旋而新自旋向下，
另一类是前 $`L`$ 个格点上有 $`N-1`$ 个向上自旋而新自旋向上。在这两个子空间上
$`x\otimes1`$ 分别作用为 $`x_N`$ 和 $`x_{N-1}`$，所以 $`\mathcal A_{L+1}`$ 的块 $`N`$ 各包含一次 $`\mathcal A_L`$ 的
块 $`N`$ 与块 $`N-1`$；维数相加正是 Pascal 法则。在二部图中，$`\mathcal A_{L+1}`$ 的扇区 $`N`$
与 $`\mathcal A_L`$ 的 $`N-1`$ 和 $`N`$ 相连，这就给出路径
$`0'-0-1'-1-\cdots-L'-L-(L+1)'`$（带撇号的属于 $`\mathcal A_{L+1}`$）。$`\square`$

从 Bratteli 图顶端到第 $`L`$ 层扇区 $`N`$ 的一条路径就是一个构型 $`b`$；这样的路径共有
$`\binom LN`$ 条，恰为块的大小。对包含矩阵为 $`\Lambda`$ 的连通包含，我们把 $`\|\Lambda\|^2`$（算子范数）
**定义**为**指标**；对 $`M_k\otimes1\subset M_k\otimes M_m`$，包含矩阵是 $`1\times1`$ 矩阵 $`(m)`$，所以指标就是
7.2 节的 $`m^2`$。`Inclusion.index` 计算它。（这个名称借自 V. Jones 的指标理论，*Invent. Math.* 72 (1983) 1；
在本课程中，有限维代数之包含的指标就是这个定义，仅此而已。）

*库中：* 指标是 `Inclusion.index`——[教程，9.1 节](../../tutorial/zh/09_spin_chains.md#s-9-1)。

<a id="r-9-3"></a>

**引理 9.3.** 顶点为 $`1,\dots,n`$ 的路径的邻接矩阵 $`A_n`$ 的特征值为
$`2\cos\frac{\pi k}{n+1}`$，$`k=1,\dots,n`$，对应的特征向量为 $`v_k(j)=\sin\frac{\pi jk}{n+1}`$；因此
$`\|A_n\|=2\cos\frac\pi{n+1}`$。

*证明。* 令 $`v_k(0)=v_k(n+1)=0`$。由
$`\sin(a-b)+\sin(a+b)=2\sin a\cos b`$，$`(A_nv_k)(j)=v_k(j-1)+v_k(j+1)=2\cos\frac{\pi k}{n+1}\,v_k(j)`$，并且由于 $`v_k(1)=\sin\frac{\pi k}{n+1}\neq0`$，有 $`v_k\ne0`$。这 $`n`$ 个特征值
互不相同（$`\cos`$ 在 $`(0,\pi)`$ 上严格递减），所以它们就是全部特征值，而对称矩阵
$`A_n`$ 的范数是特征值的最大模，在 $`k=1`$ 和 $`k=n`$ 处取到。$`\square`$

<a id="r-9-4"></a>

**推论 9.4.** $`[\mathcal A_{L+1}:\mathcal A_L]=\|\Lambda_L\|^2=4\cos^2\frac\pi{2L+4}`$，它递增趋于 $`4`$。

*证明。* 沿命题 9.2 的路径给顶点排序，邻接矩阵
$`B=\begin{pmatrix}0&\Lambda_L\\\Lambda_L^T&0\end{pmatrix}`$ 就变成 $`A_{2L+3}`$（一个置换相似）。由引理 1.1，
$`\|B\|^2=\|B^2\|`$ 且 $`\|\Lambda_L\Lambda_L^T\|=\|\Lambda_L^T\Lambda_L\|=\|\Lambda_L\|^2`$，而 $`B^2=\Lambda_L\Lambda_L^T\oplus\Lambda_L^T\Lambda_L`$，所以
$`\|B\|=\|\Lambda_L\|`$；取 $`n=2L+3`$ 的引理 9.3 给出其值。$`\square`$

*计算：* 对 $`L\le8`$，用 `inclusions.chain_tower` 和 `inclusions.bratteli_diagram` 检验命题 9.2、引理 9.3 和推论 9.4——[教程，9.1 节](../../tutorial/zh/09_spin_chains.md#s-9-1)。

> **物理意义**（物理解释，而非本课程的定理）。只要相互作用在绕 $`z`$ 轴的转动下不变，$`S^z`$ 就是守恒的；
> $`N`$ 计数向上自旋，或者在 9.4 节的 Jordan–Wigner 映射之后计数粒子。在塔中，较长链的每个荷扇区
> 都由较短链的两个扇区供给。

<a id="s-9-2"></a>

## 9.2 平移、动量与 Marshall 符号规则

在周期链上，令 $`T`$ 为平移 $`T|b_0b_1\dots b_{L-1}\rangle=|b_{L-1}b_0\dots b_{L-2}\rangle`$，
它把格点 $`i`$ 移到格点 $`i+1`$。$`T`$ 是置换矩阵，所以是酉的，$`T^L=1`$ 且
$`Ts^a_iT^*=s^a_{i+1}`$。因此 $`TS^zT^*=S^z`$，并且对均匀外场 $`h_i=h`$ 有 $`THT^*=H`$：$`T`$ 把键映为键。
因此它的特征值具有 $`e^{ik}`$ 的形式，其中**动量** $`k=2\pi m/L`$，$`m=0,\dots,L-1`$。

<a id="r-9-5"></a>

**命题 9.5（动量扇区）。** 固定 $`N`$。对含 $`N`$ 个向上自旋的构型 $`a`$，令 $`R_a`$ 为
其轨道 $`\{T^ja\}`$ 的大小。(a) $`R_a`$ 整除 $`L`$；对满足 $`mR_a\equiv0\pmod L`$ 的 $`k=2\pi m/L`$，向量
$`|a,k\rangle=\sum_{j=0}^{R_a-1}e^{-ikj}T^j|a\rangle`$ 非零且 $`T|a,k\rangle=e^{ik}|a,k\rangle`$，并且这
$`R_a`$ 个向量构成轨道张成空间的一组正交基。(b) $`T`$ 在扇区 $`N`$ 上的特征值 $`e^{2\pi im/L}`$
的重数为 $`d_{N,m}=\#\{\text{轨道 }:mR_a\equiv0\pmod L\}`$，并且每个与
$`S^z`$ 和 $`T`$ 都可交换的算子——例如平移不变的 XXZ 链——都是块对角的，属于 $`\bigoplus_N\bigoplus_mM_{d_{N,m}}`$。

*证明。* (a) $`a`$ 在循环群 $`\{T^j\}\cong\mathbb Z_L`$ 中的稳定子是一个子群，而 $`R_a`$ 是它的
指数，是 $`L`$ 的因子；向量 $`T^j|a\rangle`$（$`0\le j<R_a`$）规范正交，且 $`T^{R_a}|a\rangle=|a\rangle`$。若
$`e^{ikR_a}=1`$，则 $`T|a,k\rangle=e^{ik}\sum_{j=1}^{R_a}e^{-ikj}T^j|a\rangle=e^{ik}|a,k\rangle`$，因为 $`j=R_a`$ 的项等于
$`j=0`$ 的项；$`|a\rangle`$ 的系数为 $`1`$。容许的 $`m\in\{0,\dots,L-1\}`$ 是
$`L/R_a`$ 的倍数，共 $`R_a`$ 个，对应酉算子 $`T`$ 的互不相同的特征值，所以这些向量相互正交，并张成
轨道所张成的 $`R_a`$ 维空间。(b) 扇区 $`N`$ 是其各轨道张成空间的正交和，每个都在 $`T`$ 下不变，并且由 (a)，
若 $`mR_a\equiv0\pmod L`$，特征值 $`e^{2\pi im/L}`$ 在 $`a`$ 的轨道张成的空间中恰出现一次，否则不出现。与 $`S^z`$ 和 $`T`$
都可交换的 $`X`$ 与它们的谱投影 $`P_N`$、$`Q_m`$ 可交换（定理 1.5(d)），所以
$`X=\sum_{N,m}P_NQ_mXP_NQ_m`$。$`\square`$

*计算：* $`L=12`$ 时的命题 9.5(b)，以及 `SpinChain.xxz_momentum` 构造的动量块——[教程，9.2 节](../../tutorial/zh/09_spin_chains.md#s-9-2)。

<a id="r-9-6"></a>

**引理 9.6（Perron–Frobenius，对称情形）。** 设 $`M`$ 是实对称矩阵，对
$`i\ne j`$ 有 $`M_{ij}\ge0`$，且其图（$`i\sim j`$ 当且仅当 $`M_{ij}>0`$）连通。则其最大特征值 $`\lambda`$ 是单的，
并有一个各分量全都 $`>0`$ 的特征向量。

*证明。* $`\lambda=\max_{\|x\|=1}\langle x,Mx\rangle`$，且最大值恰好在
$`\lambda`$-特征空间的单位向量上取到（在一组规范正交特征基中展开）。对这样的 $`x`$，由各分量模构成的向量 $`|x|`$
满足 $`\langle|x|,M|x|\rangle\ge\mathrm{Re}\langle x,Mx\rangle=\lambda`$，因为对角项相同，且对 $`i\ne j`$ 有
$`M_{ij}|x_i||x_j|\ge\mathrm{Re}\,M_{ij}\bar x_ix_j`$；所以 $`|x|`$ 是特征向量。若 $`|x|_i=0`$，则
$`0=(\lambda-M_{ii})|x|_i=\sum_{j\ne i}M_{ij}|x|_j`$ 迫使 $`i`$ 的所有邻点 $`j`$ 都有 $`|x|_j=0`$，再由
连通性得 $`x=0`$。所以非零特征向量没有为零的分量；但维数 $`\ge2`$ 的特征空间
含有一个满足 $`x_1=0`$ 的非零向量。因此 $`\lambda`$ 是单的，且 $`|x|>0`$ 是特征向量。$`\square`$

我们对任意二部的耦合排布证明符号规则；链是其特例。

<a id="r-9-7"></a>

**定理 9.7（Marshall 符号规则）。** 设格点 $`\{0,\dots,L-1\}`$ 被分成两个集合 $`A\sqcup B`$，并设
```math
H=\sum_{i<j}\Big[\tfrac{J_{ij}}2\big(s^+_is^-_j+s^-_is^+_j\big)+K_{ij}\,s^z_is^z_j\Big]+\sum_ih_is^z_i
```
其中 $`K_{ij},h_i`$ 为实数，$`J_{ij}`$ 为实数，且当 $`i,j`$ 位于不同集合时 $`J_{ij}\ge0`$，位于同一集合时 $`J_{ij}\le0`$，
并且格点上的图 $`G`$（$`i\sim j`$ 当且仅当 $`J_{ij}\ne0`$）是连通的。令 $`N_B(b)`$ 为 $`b`$ 在 $`B`$ 上的向上自旋数。
在每个扇区 $`N`$ 中，$`H_N`$ 的最低特征值是单的，并且其特征向量可以选得使对每个含 $`N`$ 个向上自旋的
构型 $`b`$ 都有 $`(-1)^{N_B(b)}\psi(b)>0`$。

特别地，对 $`J>0`$、$`\Delta,h_i`$ 任意的 XXZ 链，若链是开链或 $`L`$ 为偶数的周期链，取 $`A`$ 为偶数格点、
$`B`$ 为奇数格点，则上述结论成立：此时在键上 $`J_{ij}=J`$，$`K_{ij}=J\Delta`$，其余为 $`0`$；每条键都连接一个
偶数格点和一个奇数格点（对环而言是因为 $`L`$ 为偶数），并且 $`G`$ 包含路径 $`0-1-\cdots-(L-1)`$。

*证明。* 令 $`V|b\rangle=(-1)^{N_B(b)}|b\rangle`$；$`V`$ 是对角的，$`V=V^*=V^{-1}`$。算子 $`s^+_is^-_j`$ 把一个向上自旋
从 $`j`$ 移到 $`i`$；若 $`i,j`$ 位于不同集合，它使 $`N_B`$ 改变 $`\pm1`$，否则使 $`N_B`$ 不变，所以相应地
$`Vs^+_is^-_jV=\mp s^+_is^-_j`$，对 $`s^-_is^+_j`$ 也一样。两个不同的构型 $`b,b'`$ 由这样一项相连，当且仅当 $`b'`$ 由
$`b`$ 交换两个格点 $`i,j`$ 上相反的自旋得到，而此时 $`i,j`$ 由 $`b,b'`$ 确定（即它们相异的格点）。因此
$`M=-VH_NV`$ 的非对角元为 $`J_{ij}/2`$（不同集合）或 $`-J_{ij}/2`$（同一集合），即在两种情形下都是
$`|J_{ij}|/2\ge0`$，它为正当且仅当在 $`G`$ 中 $`i\sim j`$；而在以其他方式相异的构型之间为 $`0`$。$`M`$ 的图是连通的。
事实上，用 $`(i\,j)`$ 表示交换格点 $`i,j`$ 上的自旋；对互不相同的 $`a,b,c`$ 有 $`(a\,c)=(a\,b)(b\,c)(a\,b)`$，所以沿
$`G`$ 中从 $`i`$ 到 $`j`$ 的一条路径归纳可知，每个 $`(i\,j)`$ 都是沿 $`G`$ 的边的交换的复合，而每个这样的交换要么
固定一个构型（自旋相同），要么沿 $`M`$ 的图的一条边移动。若含 $`N`$ 个向上自旋的 $`b\ne b'`$，则存在格点
$`i`$（在 $`b`$ 中向上、在 $`b'`$ 中向下）和格点 $`j`$（在 $`b`$ 中向下、在 $`b'`$ 中向上）；$`(i\,j)b`$ 与 $`b'`$ 相异的格点比
$`b`$ 少两个，所以有限多步这样的操作就从 $`b`$ 到达 $`b'`$。由引理 9.6，$`M`$ 的最大特征值——即 $`H_N`$ 最低
特征值的相反数——是单的，并有特征向量 $`\varphi>0`$；取 $`\psi=V\varphi`$。$`\square`$

<a id="r-9-8"></a>

**推论 9.8（基态的动量）。** 在 $`L`$ 为偶数、$`J>0`$ 且外场均匀的周期链上，
扇区 $`N`$ 的基态在 $`N`$ 为偶数时动量为 $`k=0`$，在 $`N`$ 为奇数时动量为 $`k=\pi`$。

*证明。* $`T`$ 与 $`H_N`$ 可交换，且基态 $`\psi=V\varphi`$ 是单的，所以 $`T\psi=e^{ik}\psi`$。而
$`(T\psi)(b)=\psi(T^{-1}b)`$，且 $`T^{-1}b`$ 在格点 $`i`$ 上的自旋为 $`b_{i+1}`$；由于 $`L`$ 为偶数，对奇数 $`i`$，$`i+1`$ 是偶数，
所以 $`N_B(T^{-1}b)=N-N_B(b)`$，且 $`(VT\psi)(b)=(-1)^{N_B(b)}(-1)^{N-N_B(b)}\varphi(T^{-1}b)=(-1)^N\varphi(T^{-1}b)`$。另一方面
$`VT\psi=e^{ik}\varphi`$。$`\varphi`$ 与 $`\varphi\circ T^{-1}`$ 都是正的，所以 $`e^{ik}=(-1)^N`$。$`\square`$

*计算：* 对带随机场的开链和周期链检验定理 9.7，对 $`L=12`$ 检验推论 9.8——[教程，9.2 节](../../tutorial/zh/09_spin_chains.md#s-9-2)。

> **物理意义**（物理解释，而非本课程的定理）。在反铁磁体中，晶格决定了基态的符号：使一个自旋越过
> 一条键会改变符号。物理学家利用这一点在二部反铁磁体上进行没有符号问题的量子 Monte Carlo 计算；
> 对于阻挫链，即有反铁磁耦合连接同一集合的两个格点时，定理 9.7 的假设不成立。

<a id="s-9-3"></a>

## 9.3 $`\Delta=1`$ 时的 $`SU(2)`$

对 $`\Delta=1`$ 和 $`h=0`$，$`H=J\sum\mathbf s_i\cdot\mathbf s_j`$ 与 $`S^x,S^y,S^z`$ 可交换（命题 4.15），所以
$`H\cong\bigoplus_S1_{2S+1}\otimes H_S`$，其中块 $`H_S\in M_{m_S}`$，$`m_S=\binom L{L/2-S}-\binom L{L/2-S-1}`$
（命题 4.13），而模块 `su2` 在最高权向量上构造 $`H_S`$（命题 4.15(b)）。
对每个 $`J_{ij}`$ 为实数的 $`H=\sum_{i<j}J_{ij}\,\mathbf s_i\cdot\mathbf s_j`$ 同样成立：由命题 4.15(a)，每个
$`\mathbf s_i\cdot\mathbf s_j`$ 都属于 $`\{S_a\}'`$。对出现的自旋 $`S`$，令 $`E_0(S)`$ 为 $`H_S`$ 的最低特征值。

*库中：* 块 $`H_S`$ 由 `SU2Chain` 构造——[教程，9.3 节](../../tutorial/zh/09_spin_chains.md#s-9-3)。

<a id="r-9-9"></a>

**定理 9.9（Lieb–Mattis）。** 设 $`H=\sum_{i<j}J_{ij}\,\mathbf s_i\cdot\mathbf s_j`$，其中 $`A,B,J_{ij}`$ 如定理 9.7 中所述
（所以 $`K_{ij}=J_{ij}`$，$`h_i=0`$），$`n_A=|A|\ge n_B=|B|`$，且 $`S_*=\frac{n_A-n_B}2`$。
(a) 对 $`S_*\le S\le\frac L2`$（$`S-\frac L2`$ 为整数），扇区 $`S^z=S`$ 的基态具有总自旋 $`S`$；
对 $`0\le M\le S_*`$（$`M-\frac L2`$ 为整数），扇区 $`S^z=M`$ 的基态具有总自旋 $`S_*`$。
(b) 因此对 $`S\ge S_*`$，$`E_0(S)`$ 是 $`H_S`$ 的单特征值，$`E_0(S_*)<E_0(S_*+1)<\dots<E_0(\frac L2)`$，并且对
$`S<S_*`$ 有 $`E_0(S)>E_0(S_*)`$：$`H`$ 的基态空间是单个多重态，其自旋为 $`S_*`$。
对 $`J>0`$ 的 Heisenberg 链（开链，或 $`L`$ 为偶数的周期链），$`A`$ 和 $`B`$ 是偶数格点和奇数格点，
$`L`$ 为偶数时 $`S_*=0`$，$`L`$ 为奇数时 $`S_*=\frac12`$（开链），并且 $`E_0(S_*)<E_0(S_*+1)<\dots<E_0(\frac L2)`$。

*证明。* (a) 设 $`M\ge0`$，$`M-\frac L2\in\mathbb Z`$，$`s=\max(M,S_*)`$，并设 $`\psi`$ 为扇区 $`S^z=M`$ 的基态，
由定理 9.7 它是单的：$`\psi=V\varphi`$，其中 $`\varphi>0`$。$`\mathbf S^2`$ 与 $`H`$ 和 $`S^z`$ 都可交换，所以它保持该扇区
的一维基态空间：$`\psi`$ 具有某个总自旋 $`S'`$。我们构造该扇区中的一个向量 $`\chi`$，使
$`\mathbf S^2\chi=s(s+1)\chi`$ 且 $`\langle\psi,\chi\rangle\ne0`$；由于 $`\mathbf S^2`$ 属于不同特征值的特征向量相互正交，
这就给出 $`S'=s`$。令 $`p=\frac L2-s`$；则 $`0\le p\le n_B`$，因为 $`s\ge S_*`$ 意味着
$`p\le\frac{n_A+n_B}2-\frac{n_A-n_B}2`$。取互不相同的 $`a_1,\dots,a_p\in A`$ 和 $`b_1,\dots,b_p\in B`$，令 $`R`$ 为其余
$`2s`$ 个格点，并令
```math
\chi_0=\bigotimes_{l=1}^p\tfrac1{\sqrt2}\big(|{\uparrow}_{a_l}{\downarrow}_{b_l}\rangle-|{\downarrow}_{a_l}{\uparrow}_{b_l}\rangle\big)\otimes|{\uparrow\cdots\uparrow}\rangle_R ,
\qquad \chi=(S^-)^{s-M}\chi_0 .
```
每个单态都被 $`s^\pm_{a_l}+s^\pm_{b_l}`$ 湮灭，$`R`$ 上的向上自旋被 $`s^+_i`$ 湮灭，所以 $`S^+\chi_0=0`$；结合
$`S^z\chi_0=s\chi_0`$ 得 $`\mathbf S^2\chi_0=(S^-S^++(S^z)^2+S^z)\chi_0=s(s+1)\chi_0`$（4.5 节），又由于 $`\mathbf S^2`$ 与 $`S^-`$
可交换，$`\mathbf S^2\chi=s(s+1)\chi`$。写 $`S^-=\sum_l(s^-_{a_l}+s^-_{b_l})+S^-_R`$，各项相互可交换，前面各项湮灭
单态，所以 $`\chi`$ 是这些单态与 $`(S^-_R)^{s-M}|{\uparrow\cdots\uparrow}\rangle_R=(s-M)!\sum_D|D\rangle`$ 之积，后者的求和
遍及 $`R`$ 上恰有 $`s-M`$ 个自旋向下的构型（这样的构型存在，因为 $`s-M\le2s=|R|`$）；于是 $`\chi\ne0`$，
$`S^z\chi=M\chi`$，且 $`\chi`$ 的系数是实的。其符号：一个单态在其向上自旋位于 $`A`$ 上时贡献 $`+`$，位于 $`B`$ 上时
贡献 $`-`$，而 $`R`$ 部分是正的。若 $`M=s`$，$`R`$ 的格点全部向上；若 $`M<s`$，则 $`s=S_*`$，$`p=n_B`$，且
$`R\subseteq A`$。在两种情形下，$`N_B(b)`$ 都等于向上自旋位于 $`B`$ 上的单态个数加上一个常数 $`|R\cap B|`$
（若 $`M=s`$）或 $`0`$，所以 $`(-1)^{N_B(b)}\chi(b)`$ 对所有满足 $`\chi(b)\ne0`$ 的 $`b`$ 符号相同。因此
$`\langle\psi,\chi\rangle=\langle\varphi,V\chi\rangle=\pm\sum_b\varphi(b)|\chi(b)|\ne0`$。

(b) 由命题 4.15(a)、推论 4.12 和定理 4.11，自旋为 $`S'`$ 的多重态在 $`S'\ge|M|`$ 时与扇区 $`S^z=M`$ 恰好交于
一个向量，否则不相交，所以 $`H`$ 在扇区 $`S^z=M\ge0`$ 上的谱是各 $`H_{S'}`$（$`S'\ge M`$）的谱的并（计重数）。
由 (a)，它的最低特征值是单的，并属于自旋 $`s=\max(M,S_*)`$：所以 $`E_0(s)`$ 在 $`H_s`$ 中是单的，并且对所有
$`S'\ge M`$、$`S'\ne s`$ 有 $`E_0(s)<E_0(S')`$。对 $`M=S\ge S_*`$，这给出单性以及对 $`S'>S`$ 有 $`E_0(S)<E_0(S')`$；
对 $`M=S<S_*`$，它给出 $`E_0(S_*)<E_0(S)`$。所以 $`E_0(S_*)`$ 是 $`H`$ 的最低特征值，只在 $`H_{S_*}`$ 中取到，并且在
其中是单的。对链而言，偶数格点比奇数格点多 $`0`$ 个（$`L`$ 为偶数）或 $`1`$ 个（$`L`$ 为奇数）。$`\square`$

定理 9.7 可追溯到 W. Marshall（*Proc. Roy. Soc. A* 232 (1955) 48），定理 9.9 可追溯到 E. Lieb 和
D. Mattis（*J. Math. Phys.* 3 (1962) 749）。由定理 9.9，偶数 Heisenberg 环的**单态–三重态能隙**
$`E_0(1)-E_0(0)`$ 为正。

*计算：* 对 $`L=12`$ 用 $`U(1)`$ 扇区检验定理 9.9，单态–三重态能隙，以及在一个不是链的二部图上检验定理 9.7 和 9.9——[教程，9.3 节](../../tutorial/zh/09_spin_chains.md#s-9-3)。

**关联。** 在每个偶数 Heisenberg 环上，基态关联 $`\langle\mathbf s_0\cdot\mathbf s_r\rangle`$ 符号交替；
这是定理 9.7 和 9.9 所强制的。设 $`\psi`$ 为实的归一化基态，它是单态
（定理 9.9），因而是扇区 $`S^z=0`$ 的基态。由 $`\mathbf S^2\psi=0`$ 得
$`\|S^\pm\psi\|^2=\langle\psi,S^\mp S^\pm\psi\rangle=\langle\psi,(\mathbf S^2-(S^z)^2\mp S^z)\psi\rangle=0`$，所以 $`S^\pm\psi=0`$，
并且对每个 $`X`$ 有 $`\langle\psi,[S^+,[S^-,X]]\psi\rangle=0`$（四项中的每一项都有 $`S^\pm`$ 紧挨着 $`\psi`$）。对
$`X=s^z_0s^z_r`$（$`r\ne0`$），关系式 $`[s^\pm,s^z]=\mp s^\pm`$、$`[s^+,s^-]=2s^z`$ 给出 $`[S^-,X]=s^-_0s^z_r+s^z_0s^-_r`$ 和
$`[S^+,[S^-,X]]=4s^z_0s^z_r-s^+_0s^-_r-s^-_0s^+_r`$。因此 $`\langle s^z_0s^z_r\rangle=\frac14\langle s^+_0s^-_r+s^-_0s^+_r\rangle`$，并且
```math
\langle\mathbf s_0\cdot\mathbf s_r\rangle=\langle s^z_0s^z_r\rangle+\tfrac12\langle s^+_0s^-_r+s^-_0s^+_r\rangle
=\tfrac34\sum_{b_0\ne b_r}\psi(b)\psi(b') ,
```
其中 $`b'`$ 由 $`b`$ 交换格点 $`0`$ 与 $`r`$ 的自旋得到。由定理 9.7，$`\psi(b)\psi(b')`$ 的符号与
$`(-1)^{N_B(b)+N_B(b')}`$ 相同，而 $`N_B(b')-N_B(b)`$ 对奇数 $`r`$ 为 $`\pm1`$、对偶数 $`r`$ 为 $`0`$；该求和非空，
所以对 $`0<r<L`$，$`\langle\mathbf s_0\cdot\mathbf s_r\rangle`$ 非零且符号为 $`(-1)^r`$。

*计算：* 用 `SU2Chain.spin_correlation` 计算 $`L=12`$ 时的关联 $`\langle\mathbf s_0\cdot\mathbf s_r\rangle`$，它利用了 Wigner–Eckart 定理（定理 B.16）——[教程，9.3 节](../../tutorial/zh/09_spin_chains.md#s-9-3)。

<a id="s-9-4"></a>

## 9.4 精确检验：自由费米子与 Bethe 拟设

**Jordan–Wigner。** 在开链上令 $`c_j=\big(\prod_{l<j}(-1)^{n_l}\big)s^-_j`$。

<a id="r-9-10"></a>

**定理 9.10（XX 链是自由的）。** (a) $`\{c_i,c_j^*\}=\delta_{ij}`$，$`\{c_i,c_j\}=0`$，且 $`c_j^*c_j=n_j`$。
(b) 开放 XX 链为 $`H=\frac J2\sum_{j<L-1}(s^+_js^-_{j+1}+s^-_js^+_{j+1})=\sum_{i,j}h_{ij}c_i^*c_j`$，其中
$`h=\frac J2A_L`$。(c) 设 $`u_k`$ 是 $`h`$ 的规范正交实特征向量，特征值为
$`\varepsilon_k=J\cos\frac{\pi k}{L+1}`$（引理 9.3），并令 $`b_k=\sum_ju_k(j)c_j`$。对 $`K=\{k_1<\dots<k_N\}`$，向量
$`b_K^*|0\rangle=b_{k_1}^*\cdots b_{k_N}^*|0\rangle`$（$`|0\rangle`$ 为所有自旋向下）构成扇区
$`N`$ 的一组规范正交特征基，能量为 $`\sum_{k\in K}\varepsilon_k`$。

*证明。* (a) 不同格点上的算子可交换，而在同一格点上 $`(-1)^n`$ 与 $`s^\pm`$ 反对易。对
$`i<j`$，$`c_i`$ 与 $`c_j^*`$（或 $`c_j`$）在格点 $`i`$ 上以 $`s^-_i`$ 和 $`(-1)^{n_i}`$ 的形式相遇——一次反对易——
在其他地方可交换，所以它们反对易。对 $`i=j`$，弦算子的平方为 $`1`$，且
$`c_jc_j^*+c_j^*c_j=s^-_js^+_j+s^+_js^-_j=1`$，$`c_j^*c_j=s^+_js^-_j=n_j`$，$`c_j^2=(s^-_j)^2=0`$；$`i>j`$ 的情形
交换 $`i,j`$ 的角色（并取伴随）即得。(b) $`c_j^*c_{j+1}=s^+_j(-1)^{n_j}s^-_{j+1}=s^+_js^-_{j+1}`$，
因为 $`s^+_j`$ 仅在 $`n_j=0`$ 上非零；再取伴随。(c) $`\{b_k,b_l^*\}=\sum_ju_k(j)u_l(j)=\delta_{kl}`$ 且
$`\{b_k,b_l\}=0`$，并由 $`h`$ 的谱分解得 $`H=\sum_k\varepsilon_kb_k^*b_k`$。由正则反对易关系（CAR），
$`[b_l^*b_l,b_k^*]=\delta_{lk}b_k^*`$，所以 $`[H,b_k^*]=\varepsilon_kb_k^*`$；结合 $`H|0\rangle=0`$ 得
$`Hb_K^*|0\rangle=(\sum_K\varepsilon_k)b_K^*|0\rangle`$。由于 $`b_k|0\rangle=0`$（每个 $`c_j`$ 降低一个自旋），用 CAR 把湮灭算子
移到右边，对 $`N`$ 归纳即得 $`\langle b_{K'}^*0,b_K^*0\rangle=\delta_{KK'}`$。每个 $`c_j^*`$ 使
$`N`$ 增加一，所以 $`b_K^*|0\rangle`$ 位于扇区 $`|K|`$ 中：在维数为 $`\binom LN`$ 的空间中有 $`\binom LN`$ 个规范正交向量。$`\square`$

*计算：* $`L=10`$ 的开 XX 链的全部 $`1024`$ 个能级，与定理 9.10 对照——[教程，9.4 节](../../tutorial/zh/09_spin_chains.md#s-9-4)。

**磁振子。** 在 $`h=0`$ 的环上，完全极化态 $`|0\rangle`$ 的能量为 $`E_F=J\Delta L/4`$。位于 $`x`$ 的一个向上
自旋破坏两条键：$`H|x\rangle=(E_F-J\Delta)|x\rangle+\frac J2(|x-1\rangle+|x+1\rangle)`$（$`L\ge3`$）。

<a id="r-9-11"></a>

**命题 9.11（单磁振子）。** 在环的扇区 $`N=1`$ 中，$`|k\rangle=\sum_xe^{-ikx}|x\rangle`$，
$`k=2\pi m/L`$，满足 $`T|k\rangle=e^{ik}|k\rangle`$ 和 $`H|k\rangle=\big(E_F+J(\cos k-\Delta)\big)|k\rangle`$。

*证明。* 利用 $`e^{ikL}=1`$，$`\sum_xe^{-ikx}(|x-1\rangle+|x+1\rangle)=\sum_y(e^{-ik(y+1)}+e^{-ik(y-1)})|y\rangle=2\cos k\,|k\rangle`$；
同样的重新标号给出 $`T|k\rangle=e^{ik}|k\rangle`$。$`\square`$

两个磁振子会发生散射。对 $`x_1<x_2`$，用 $`|x_1,x_2\rangle`$ 表示向上自旋位于 $`x_1,x_2`$ 的态；一对相邻的向上自旋
破坏两条键，一对分开的向上自旋破坏四条。

<a id="r-9-12"></a>

**命题 9.12（双磁振子；Bethe 方程）。** 设 $`k_1,k_2\in\mathbb C`$，$`K=k_1+k_2`$，且
$`(A_{12},A_{21})\neq0`$ 满足 (i) $`A_{12}(1+e^{iK}-2\Delta e^{ik_2})+A_{21}(1+e^{iK}-2\Delta e^{ik_1})=0`$ 和
(ii) $`A_{21}=e^{ik_2L}A_{12}`$，$`A_{12}=e^{ik_1L}A_{21}`$。令
$`\psi(x_1,x_2)=A_{12}e^{i(k_1x_1+k_2x_2)}+A_{21}e^{i(k_2x_1+k_1x_2)}`$（定义在 $`\mathbb Z^2`$ 上），并令
$`\Psi=\sum_{0\le x_1<x_2<L}\psi(x_1,x_2)|x_1,x_2\rangle`$。若 $`\Psi\neq0`$，则
$`H\Psi=\big(E_F+J(\cos k_1+\cos k_2-2\Delta)\big)\Psi`$。若分母不为零，则 (i) 和 (ii) 等价于
```math
e^{ik_1L}=-\frac{1+e^{iK}-2\Delta e^{ik_1}}{1+e^{iK}-2\Delta e^{ik_2}},\qquad
e^{ik_2L}=-\frac{1+e^{iK}-2\Delta e^{ik_2}}{1+e^{iK}-2\Delta e^{ik_1}} .
```

*证明。* 令 $`E=E_F+J(\cos k_1+\cos k_2-2\Delta)`$。(1) 对所有 $`(x_1,x_2)\in\mathbb Z^2`$，每个平面波给出的
$`\psi(x_1\pm1,x_2)`$ 与 $`\psi(x_1,x_2\pm1)`$ 之和为 $`2(\cos k_1+\cos k_2)\psi(x_1,x_2)`$。(2) 代入可知，
$`\psi(x,x)+\psi(x+1,x+1)-2\Delta\psi(x,x+1)`$ 等于 $`e^{iKx}`$ 乘以 (i) 的左边，所以它对所有 $`x`$ 为零。
(3) 由 (ii)，$`\psi(x_2,x_1+L)=A_{12}e^{ik_2L}e^{i(k_2x_1+k_1x_2)}+A_{21}e^{ik_1L}e^{i(k_1x_1+k_2x_2)}=\psi(x_1,x_2)`$，
两次应用此式得 $`\psi(x_1+L,x_2+L)=\psi(x_1,x_2)`$。构型 $`\{y_1<y_2\}`$ 在
$`D=\{x_1<x_2<x_1+L\}`$ 中的代表元是 $`(y_1+nL,y_2+nL)`$ 与 $`(y_2+nL,y_1+(n+1)L)`$，所以 $`\psi`$ 是构型的函数。
设 $`(x_1,x_2)\in D`$。若 $`2\le x_2-x_1\le L-2`$，四种跳跃都留在 $`D`$ 中，对角元为 $`E_F-2J\Delta`$，所以由 (1)，
$`(H\Psi)(x_1,x_2)=(E_F-2J\Delta)\psi+\frac J2\cdot2(\cos k_1+\cos k_2)\psi=E\psi`$。若 $`x_2=x_1+1`$，跳跃
$`(x_1+1,x_2)`$、$`(x_1,x_2-1)`$ 被阻挡，对角元为 $`E_F-J\Delta`$；由 (1) 和 (2)，
$`(H\Psi)(x_1,x_2)=(E_F-J\Delta)\psi+\frac J2\big(2(\cos k_1+\cos k_2)\psi-\psi(x_1,x_1)-\psi(x_1+1,x_1+1)\big)=E\psi`$。
若 $`x_2=x_1+L-1`$，则用代表元 $`(x_2,x_1+L)`$，它属于前一种情形（$`L\ge3`$）。$`\square`$

命题 9.12 只说明解给出特征向量；它并没有说所有特征向量都以这种方式得到。
记 $`k_{1,2}=K/2\pm q`$ 和 $`c=\cos\frac K2`$。则 $`1+e^{iK}-2\Delta e^{ik_{1,2}}=2e^{iK/2}(c-\Delta e^{\pm iq})`$，
所以 $`A_{12}=c-\Delta e^{iq}`$，$`A_{21}=-(c-\Delta e^{-iq})`$ 是 (i) 的解。由 (ii) 和 $`(A_{12},A_{21})\ne0`$，
两个系数都非零，所以把 (ii) 的两个方程相乘得
$`e^{iKL}=1`$，$`K=2\pi n/L`$，于是 (ii) 等价于 $`F(q)=(-1)^ne^{iqL}(c-\Delta e^{-iq})+c-\Delta e^{iq}=0`$，
这是一个单变量方程。我们在
$`q\in(0,\pi)`$ 中寻找它的根，这里由于 $`c-\Delta e^{\pm iq}`$ 互为复共轭，$`e^{-iqL/2}F(q)`$ 为实数或纯虚数
（散射态）；我们也在直线 $`q=i\kappa`$、$`\pi+i\kappa`$（$`\kappa>0`$）上寻找，那里 $`F`$ 是实的（束缚态）。

对偶数 $`L`$ 和 $`\Delta\ne0`$，Bethe 向量不是完备的。若 $`\Delta\ne0`$，
$`\Phi=\sum_x(-1)^x|x,x+1\rangle`$ 满足 $`T\Phi=-\Phi`$ 和 $`H\Phi=(E_F-J\Delta)\Phi`$（到 $`|x-1,x+1\rangle`$ 和
$`|x,x+2\rangle`$ 的跳跃成对抵消）。Bethe 向量满足 $`T\Psi=e^{-iK}\Psi`$，因为 $`\psi(x_1-1,x_2-1)=e^{-iK}\psi(x_1,x_2)`$。若
$`e^{-iK}\ne-1`$，则 $`\Psi\perp\Phi`$，因为它们是酉算子 $`T`$ 属于不同特征值的特征向量；若
$`e^{iK}=-1`$，则 $`\cos k_1+\cos k_2=2\cos\frac K2\cos\frac{k_1-k_2}2=0`$，$`\Psi`$ 的能量为 $`E_F-2J\Delta\ne E_F-J\Delta`$，
从而 $`\Psi\perp\Phi`$，因为它们是 $`H`$ 属于不同特征值的特征向量。所以 $`\Phi`$ 与命题 9.12 的所有 Bethe 向量正交，
因此这些向量永远张不成偶数环的扇区 $`N=2`$。双磁振子的计算是 H. Bethe 拟设（*Z. Phys.* 71 (1931) 205）的 $`N=2`$ 情形；我们不使用其他情形。

*计算：* $`L=9`$、$`\Delta=0.6`$ 时 $`F`$ 的根，它们给出扇区 $`N=2`$ 的全部 $`36`$ 个态，以及 $`L=8`$ 时的向量 $`\Phi`$——[教程，9.4 节](../../tutorial/zh/09_spin_chains.md#s-9-4)。

> **物理意义**（物理解释，而非本课程的定理）。磁振子是自旋波：一个翻转的自旋以动量 $`k`$ 离域在
> 整个环上。两个磁振子发生带相移的散射（实的 $`q`$）或形成束缚态（复的 $`q`$）。

<a id="s-9-5"></a>

## 9.5 Lanczos 与 Krylov 方法

扇区 $`N=L/2`$ 的维数在 $`L=16`$ 时为 $`12870`$，在 $`L=20`$ 时为 $`184756`$，但 $`H`$ 是稀疏的（每行 $`O(L)`$ 个非零
元素，`xxz_sparse`），而 Krylov 方法只需要乘积 $`Hv`$。设 $`H`$ 是 $`n`$ 阶厄米矩阵，
$`v\ne0`$。**Krylov 空间**为 $`\mathcal K_m=\mathrm{span}\{v,Hv,\dots,H^{m-1}v\}`$。

<a id="r-9-13"></a>

**命题 9.13（Lanczos）。** 设 $`\dim\mathcal K_m=m`$。令 $`q_1=v/\|v\|`$，并对 $`j<m`$ 令 $`q_{j+1}`$ 为
$`Hq_j`$ 正交于 $`\mathcal K_j`$ 的分量归一化后的向量。则 $`q_1,\dots,q_m`$ 是
$`\mathcal K_m`$ 的一组规范正交基；记 $`Q_m=(q_1\cdots q_m)`$，则 $`T_m=Q_m^*HQ_m`$ 是实对称三对角矩阵，其对角元为
$`\alpha_j=\langle q_j,Hq_j\rangle`$，非对角元为 $`\beta_j=\langle q_{j+1},Hq_j\rangle>0`$；并且 $`HQ_m=Q_mT_m+r_me_m^T`$，其中
$`r_m\perp\mathcal K_m`$ 为零当且仅当 $`\mathcal K_{m+1}=\mathcal K_m`$（否则 $`r_m=\beta_mq_{m+1}`$，$`q_{m+1}=r_m/\|r_m\|`$）。

*证明。* 设 $`\mathrm{span}\{q_1,\dots,q_j\}=\mathcal K_j`$（对 $`j=1`$ 成立）。则 $`H\mathcal K_j\subseteq\mathcal K_j+\mathbb CHq_j`$，
因为 $`Hq_i\in\mathcal K_{i+1}`$；所以 $`\mathcal K_{j+1}=\mathcal K_j+\mathbb CH^jv\subseteq\mathcal K_j+\mathbb CHq_j\subseteq\mathcal K_{j+1}`$，从而
$`\mathcal K_{j+1}=\mathcal K_j+\mathbb CHq_j`$。对 $`j<m`$ 维数增长，所以 $`Hq_j`$ 正交于
$`\mathcal K_j`$ 的分量非零，归纳得以继续；对 $`j=m`$，该分量 $`r_m`$ 为零当且仅当
$`\mathcal K_{m+1}=\mathcal K_m`$。对 $`i>j+1`$，$`q_i\perp\mathcal K_{j+1}\ni Hq_j`$，所以 $`(T_m)_{ij}=0`$，又由于 $`T_m`$ 是
厄米的，$`(T_m)_{ji}=0`$。$`\alpha_j`$ 是实的，而 $`\beta_j=\langle q_{j+1},Hq_j\rangle`$ 是该分量的范数，$`>0`$。所以
对 $`j<m`$ 有 $`Hq_j=\beta_{j-1}q_{j-1}+\alpha_jq_j+\beta_jq_{j+1}`$，且 $`Hq_m=\beta_{m-1}q_{m-1}+\alpha_mq_m+r_m`$。$`\square`$

$`T_m`$ 的特征值 $`\theta^{(m)}_1\le\dots\le\theta^{(m)}_m`$ 称为 **Ritz 值**。

<a id="r-9-14"></a>

**定理 9.14（单调收敛）。** 设 $`\lambda_1\le\dots\le\lambda_n`$ 为 $`H`$ 的特征值，且
$`\dim\mathcal K_m=m`$。则 $`\theta^{(m)}_j\ge\lambda_j`$，并且若 $`\dim\mathcal K_{m+1}=m+1`$，则 $`\theta^{(m)}_j\ge\theta^{(m+1)}_j`$
（$`j\le m`$）。若 $`v`$ 在 $`\lambda_1`$ 的特征空间中有非零分量且 $`\mathcal K_{m+1}=\mathcal K_m`$，则
$`\theta^{(m)}_1=\lambda_1`$。

*证明。* 由推论 A.9（压缩的 Courant–Fischer 定理，定理 A.8）应用于 $`Q=Q_m`$，
$`\theta^{(m)}_j=\min_{U}\max_{x\in U,\|x\|=1}\langle x,Hx\rangle`$，其中 $`U`$ 取遍 $`\mathcal K_m=\mathrm{ran}\,Q_m`$ 的 $`j`$ 维子空间。
同一公式对 $`U\subseteq\mathcal K_{m+1}`$ 给出 $`\theta^{(m+1)}_j`$，对 $`U\subseteq\mathbb C^n`$ 给出 $`\lambda_j`$
（定理 A.8）；在更多子空间上取最小给出更小的值。若
$`\mathcal K_{m+1}=\mathcal K_m`$，则 $`\mathcal K_m`$ 是 $`H`$-不变的，并对每个多项式 $`p`$ 含有 $`p(H)v`$，
特别地含有 $`P_1v\ne0`$（$`P_1`$ 为 $`\lambda_1`$ 的谱投影，由 Lagrange
插值它是 $`H`$ 的多项式，如定理 1.5 中那样）；所以在 $`x=P_1v/\|P_1v\|\in\mathcal K_m`$ 处 $`\theta^{(m)}_1\le\langle x,Hx\rangle=\lambda_1`$。$`\square`$

小的残差保证附近有一个特征值：对厄米的 $`H`$、数 $`E`$ 和 $`\|\psi\|=1`$，
$`\|H\psi-E\psi\|^2=\sum_i|\lambda_i-E|^2|\langle u_i,\psi\rangle|^2\ge\min_i|\lambda_i-E|^2`$（$`u_i`$ 为一组规范正交特征基）。
它并不保证这个特征值是最低的。对 $`L=4`$，Heisenberg 环的基态能量是精确的：令 $`\mathbf S_A=\mathbf s_0+\mathbf s_2`$，
$`\mathbf S_B=\mathbf s_1+\mathbf s_3`$，则环为
$`H=J\,\mathbf S_A\cdot\mathbf S_B=\frac J2(\mathbf S^2-\mathbf S_A^2-\mathbf S_B^2)\ge\frac J2(0-2-2)=-2J`$，因为 $`\mathbf S^2\ge0`$，而
$`\mathbf S_A^2,\mathbf S_B^2`$ 的特征值为 $`0,2`$（取 $`L=2`$ 的命题 4.13）；在两个三重态之积中的一个自旋为 $`0`$ 的
向量上等号成立，由第 4 章习题 2（$`1\otimes1\cong2\oplus1\oplus0`$），这样的向量存在。

*计算：* 对开 XX 链检验定理 9.14，在 `SpinChain.xxz_sparse` 的稀疏扇区上使用 `krylov.lanczos`；用 `krylov.ground_state` 及其残差计算直到 $`L=16`$ 的 Heisenberg 环的基态——[教程，9.5 节](../../tutorial/zh/09_spin_chains.md#s-9-5)。

**时间演化。** 对 $`e^{-i\tau H}v`$，我们投影到 $`\mathcal K_m`$ 上：$`w_\tau=\|v\|\,Q_me^{-i\tau T_m}e_1`$。

<a id="r-9-15"></a>

**命题 9.15.** 设 $`\dim\mathcal K_m=m`$。对每个 $`\tau`$，$`\|w_\tau\|=\|v\|`$ 且 $`\langle w_\tau,Hw_\tau\rangle=\langle v,Hv\rangle`$。若
$`\mathcal K_m`$ 是 $`H`$-不变的，则 $`w_\tau=e^{-i\tau H}v`$。

*证明。* $`Q_m`$ 是等距，$`e^{-i\tau T_m}`$ 是酉的（$`T_m`$ 是实对称的），这给出范数的结论。
由于 $`Q_m^*HQ_m=T_m`$ 与 $`e^{-i\tau T_m}`$ 可交换，
$`\langle w_\tau,Hw_\tau\rangle=\|v\|^2\langle e^{-i\tau T_m}e_1,T_me^{-i\tau T_m}e_1\rangle=\|v\|^2\langle e_1,T_me_1\rangle=\langle v,Hv\rangle`$。若
$`\mathcal K_m`$ 是不变的，则 $`HQ_m=Q_mT_m`$（命题 9.13，$`r_m=0`$），所以 $`H^jQ_m=Q_mT_m^j`$，再对指数级数求和
（命题 8.2(1)）得 $`e^{-i\tau H}Q_m=Q_me^{-i\tau T_m}`$；将其作用于
$`e_1`$，并注意 $`Q_me_1=v/\|v\|`$。$`\square`$

*计算：* 对在 Krylov 空间中逐步推进的 `krylov.evolve` 检验命题 9.15，用的是 Néel 态的淬火，并与稠密传播子对照——[教程，9.5 节](../../tutorial/zh/09_spin_chains.md#s-9-5)。

<a id="s-9-6"></a>

## 9.6 无序：能级统计与纠缠

**随机场 Heisenberg 链**的外场 $`h_i`$ 相互独立，在 $`[-W,W]`$ 上均匀分布。它保持 $`S^z`$ 守恒
（命题 9.1），所以能级统计按扇区 $`N`$ 分别计算（定理 5.15 和 5.4 节：混合扇区会伪造出独立的能级）。
若两个外场不同，则 $`THT^*-H=\sum_i(h_{i-1}-h_i)s^z_i\ne0`$，所以平移不再是对称性。

*计算：* 随机场链的平均比值 $`\langle r\rangle`$，与 GOE 矩阵以及独立能级比较（数值观察）——[教程，9.6 节](../../tutorial/zh/09_spin_chains.md#s-9-6)。

**自由费米子的纠缠。** 对开放 XX 链的本征态，前 $`n`$ 个格点的纠缠熵
由一个 $`n\times n`$ 矩阵决定。

<a id="r-9-16"></a>

**定理 9.16（Peschel 公式）。** 设 $`\Phi=b_{k_1}^*\cdots b_{k_N}^*|0\rangle`$ 如定理 9.10 中所述，
$`C_{ij}=\langle\Phi,c_i^*c_j\Phi\rangle`$，$`C_A=(C_{ij})_{i,j<n}`$。则 $`0\le C_A\le1`$，并且格点
$`A=\{0,\dots,n-1\}`$ 的纠缠熵为（$`\nu`$ 取遍 $`C_A`$ 的特征值，计重数）
```math
S_A=\sum_\nu \eta(\nu),\qquad \eta(\nu)=-\nu\log\nu-(1-\nu)\log(1-\nu)\quad(0\log0=0).
```

*证明。* 对实向量 $`x`$ 令 $`b(x)=\sum_jx_jc_j`$；则 $`b(x)^*b(y)^*=-b(y)^*b(x)^*`$，且
$`\{b(x),b(y)^*\}=\langle x,y\rangle`$。令 $`U\subseteq\mathbb R^L`$ 为各 $`u_{k_a}`$ 张成的空间，$`P_A`$、$`P_B=1-P_A`$ 为到
$`\{0,\dots,n-1\}`$ 及其补集上的坐标投影。在 $`U`$ 上对角化形式 $`\langle P_Ax,P_Ay\rangle`$，得到 $`U`$ 的一组规范正交基
$`\phi_1,\dots,\phi_N`$，满足 $`\langle P_A\phi_a,P_A\phi_b\rangle=\nu_a\delta_{ab}`$，其中 $`\nu_a=\|P_A\phi_a\|^2\in[0,1]`$，因而
$`\langle P_B\phi_a,P_B\phi_b\rangle=(1-\nu_a)\delta_{ab}`$。由于 $`(x_1,\dots,x_N)\mapsto b(x_1)^*\cdots b(x_N)^*|0\rangle`$ 是
多重线性且交错的，$`\Phi=\pm\prod_a\big(\sqrt{\nu_a}f_a^*+\sqrt{1-\nu_a}\,g_a^*\big)|0\rangle`$，其中
$`f_a=b(P_A\phi_a/\sqrt{\nu_a})`$，$`g_a=b(P_B\phi_a/\sqrt{1-\nu_a})`$（系数为零的项略去）。
展开得 $`\Phi=\sum_I\pm\prod_{a\in I}\sqrt{\nu_a}\prod_{a\notin I}\sqrt{1-\nu_a}\,F_I^*G_{I^c}^*|0\rangle`$，求和遍及
子集 $`I`$，其中 $`F_I^*`$、$`G_{I^c}^*`$ 为有序乘积。对 $`j\ge n`$，$`c_j^*`$ 的弦包含整个
$`A`$：$`c_j^*=(-1)^{N_A}\otimes\tilde c_j^*`$，其中 $`\tilde c_j`$ 是格点 $`\ge n`$ 的 Jordan–Wigner 算子，而对
$`i<n`$，$`c_i^*=\hat c_i^*\otimes1`$。由于 $`(-1)^{N_A}`$ 在 $`|0_A\rangle`$ 上为 $`1`$，
$`F_I^*G_{I^c}^*|0\rangle=\xi_I\otimes\chi_{I^c}`$，其中 $`\xi_I=\hat F_I^*|0_A\rangle`$，$`\chi_{I^c}=\tilde G_{I^c}^*|0_B\rangle`$，由定理 9.10 的 CAR 计算，
它们都是规范正交族。这是一个 Schmidt 分解（引理 7.7），其
权重 $`p_I=\prod_{a\in I}\nu_a\prod_{a\notin I}(1-\nu_a)`$（符号归入 $`\xi_I`$），即独立 Bernoulli 变量的分布，所以，由于
$`\log p_I=\sum_{a\in I}\log\nu_a+\sum_{a\notin I}\log(1-\nu_a)`$ 且 $`\sum_{I\ni a}p_I=\nu_a`$，
$`S_A=-\sum_Ip_I\log p_I=\sum_a\eta(\nu_a)`$。最后 $`C_{ij}=\sum_a\phi_a(i)\phi_a(j)`$
（把 $`c_i^*c_j`$ 用 $`b(\phi_a)`$ 展开），所以 $`C_A=\sum_a(P_A\phi_a)(P_A\phi_a)^T`$ 在
正交向量 $`P_A\phi_a`$ 上的特征值为 $`\nu_a`$，在其他地方为 $`0`$，而 $`\eta(0)=0`$。$`\square`$

定理 9.16 归功于 I. Peschel（*J. Phys. A* 36 (2003) L205）。

*计算：* 对 $`L=12`$ 的开 XX 链的基态检验定理 9.16，与 `SpinChain.entanglement_entropy` 对照——[教程，9.6 节](../../tutorial/zh/09_spin_chains.md#s-9-6)。

<a id="s-9-7"></a>

## 9.7 可与之比较的工具

*计算：* 与其他软件（QuSpin、精确对角化程序、张量网络）的比较——[教程，9.7 节](../../tutorial/zh/09_spin_chains.md#s-9-7)。

<a id="x-9"></a>

## 习题

1. （手算动量扇区。）对 $`L=4`$ 和 $`N=2`$，求六个构型在 $`T`$ 下的轨道
   以及维数 $`d_{2,m}`$，$`m=0,\dots,3`$。
2. （求和规则。）设 $`\psi`$ 为 $`L`$ 为偶数的 Heisenberg 环的基态。证明
   $`E_0=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle`$ 且 $`\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle=0`$。
3. （对称性过多。）证明在扇区 $`N=L/2`$ 中，所有系数都等于 $`1`$ 的向量是
   $`(S^-)^{L/2}|{\uparrow\cdots\uparrow}\rangle`$ 的倍数，并且是 Heisenberg 环的能量为
   $`JL/4`$ 的特征向量；再证明每个 Krylov 空间 $`\mathcal K_m(H,v)`$ 都包含于每个含有 $`v`$ 的 $`H`$-不变子空间中。
   以这个向量为起点时，`krylov.ground_state` 返回什么？

<details>
<summary>解答</summary>

1. $`0011\to1001\to1100\to0110\to0011`$ 是一条 $`R=4`$ 的轨道，$`0101\to1010\to0101`$ 是一条 $`R=2`$ 的轨道。由
   命题 9.5，前者对每个 $`m`$ 都有贡献，后者对满足 $`2m\equiv0\pmod4`$ 即 $`m=0,2`$ 的 $`m`$ 有贡献：
   $`d_{2,m}=2,1,2,1`$。

2. 基态是单的单态（定理 9.9），且 $`T`$ 与 $`H`$ 可交换，所以 $`T\psi=e^{ik}\psi`$。由于
   $`\mathbf s_i\cdot\mathbf s_{i+1}=T^i(\mathbf s_0\cdot\mathbf s_1)T^{*i}`$，$`\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle`$ 对全部 $`L`$ 条键都相同，从而
   $`E_0=J\sum_i\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle`$。由于 $`\mathbf S^2\psi=0`$，
   再次利用平移不变性，$`0=\langle\mathbf S^2\rangle=\sum_{i,j}\langle\mathbf s_i\cdot\mathbf s_j\rangle=L\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle`$。

3. $`S^-`$ 降低其中一个向上自旋；把它作用于完全极化态 $`L/2`$ 次，会以 $`(L/2)!`$ 种次序到达每个
   含 $`L/2`$ 个向上自旋的构型，每种次序的系数都是 $`1`$。$`|{\uparrow\cdots\uparrow}\rangle`$ 的全部
   $`L`$ 条键都是平行的，能量为 $`JL/4`$，而 $`S^-`$ 与 $`H`$ 可交换，所以该向量是具有相同
   能量的特征向量。若 $`W`$ 是不变的且 $`v\in W`$，则对所有 $`j`$ 有 $`H^jv\in W`$。这里对所有 $`m`$ 有 $`\mathcal K_m=\mathbb Cv`$，
   `ground_state` 返回 $`E=JL/4`$，即谱的顶端（由命题 4.15(a)，每个
   $`\mathbf s_i\cdot\mathbf s_j=\frac12F_{ij}-\frac14\le\frac14`$，因为 $`F_{ij}`$ 是特征值为 $`\pm1`$ 的对称置换矩阵）。随机的初始向量
   可以避开这个陷阱：高斯向量 $`v`$ 有密度，所以它几乎必然在 $`\lambda_1`$ 的特征空间中有非零分量 $`P_1v`$
   （推论 A.7(c)；对实的 $`v`$，在 $`W=\mathbb R^n`$ 中把它应用于真子空间 $`\{v:P_1v=0\}`$，它之所以是真子空间，
   是因为实矩阵 $`P_1\ne0`$ 有一个非零列），于是由定理 9.14，$`H`$-不变的 Krylov 空间 $`\mathcal K_m(H,v)`$ 的
   最低 Ritz 值为 $`\lambda_1`$，而这里不变的 $`\mathcal K_1=\mathbb Cv`$ 的 Ritz 值为 $`JL/4`$。

</details>

*计算：* 习题 1–3 的数值检验——[教程，习题](../../tutorial/zh/09_spin_chains.md#x-9)。

**下一章：** [第 10 章](10_quantum_computing.md)用扇区讨论量子计算：无退相干
子空间与无噪声子系统、算子代数纠错，以及保持守恒的 Hamming 重量的混合算子——
本章的代数 $`\bigoplus_wM_{\binom nw}`$ 以新的角色出现。
