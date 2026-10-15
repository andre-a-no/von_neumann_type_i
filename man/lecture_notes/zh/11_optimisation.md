[English](../11_optimisation.md) · [Español](../es/11_optimisation.md) · [Français](../fr/11_optimisation.md) · [Deutsch](../de/11_optimisation.md) · **中文** · [日本語](../ja/11_optimisation.md) · [Русский](../ru/11_optimisation.md)

# 11. 优化与迹不等式

迹不等式是关于算子泛函极值的命题。3.7 节给出了一个看似合理却不成立的迹不等式（命题 3.9）；一次采样运行（[教程，3.7 节](../../tutorial/zh/03_traces_and_states.md#s-3-7)）在数值上没有找到反例。本章解释为什么采样在高维中会失败，并以在精确参数化的约束集上的优化
取而代之；在极值有闭式表达的情形，我们给出证明，使优化器的结果能与定理相对照。

**你将学到**如何对酉算子以及具有给定对比度的谱进行参数化，以及为什么每个容许算子都能被取到；
冯·诺依曼迹不等式及由它给出的精确极值；一个带有取等条件的迹不等式；为什么 Haar 样本在均方意义下远离
极值；以及多起点搜索究竟证明了什么。本章的数值检验以及本库的优化器见
[教程 11](../../tutorial/zh/11_optimisation.md)。

在本章中，$`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$，$`D=\sum_ck_c\ge2`$，$`\mathrm{Tr}`$ 是
blunt 迹（3.1 节），$`\|A\|_1=\sum_c\|A_c\|_1`$ 是把 $`A`$ 看作
$`M_D(\mathbb C)`$ 中块对角矩阵时的迹范数。$`A\in M_k(\mathbb C)`$ 的奇异值为 $`s_1(A)\ge\dots\ge s_k(A)`$；$`x^\downarrow`$ 与
$`x^\uparrow`$ 分别是实向量 $`x`$ 的递减重排与递增重排。$`\Delta`$ 是 Michelson
对比度（3.6 节），对 $`\delta\in[0,1)`$ 令
```math
\ell(\delta)=\frac{1-\delta}{1+\delta},\qquad 1-\ell(\delta)=\frac{2\delta}{1+\delta};
```
$`\lambda_{\max}=1`$ 的正算子的对比度为 $`\ell(\lambda_{\min})`$，因此它的对比度为 $`\delta`$ 当且仅当
$`\lambda_{\min}=\ell(\delta)`$，因为 $`\ell`$ 是 $`[0,1]`$ 上严格递减的对合。例子是配套论文中的三个泛函
（*From sampling to optimisation* 一节）：

| | 泛函 | 变量 |
|---|---|---|
| 实验 1 | $`z=\lvert\mathrm{Tr}(XUY)\rvert-\mathrm{Tr}(XY)`$ | $`X,Y\ge0`$，$`U\in\mathcal M`$ 为酉算子 |
| 实验 2 | $`z=\mathrm{Tr}\lvert XY\rvert-\mathrm{Tr}(\lvert X\rvert\,\lvert Y\rvert)`$ | $`X=UX_0`$，$`Y=VY_0`$，$`X_0,Y_0\ge0`$，$`U,V`$ 为酉算子 |
| 实验 3 | $`z=\mathrm{Tr}(Y\lvert X\rvert Y)-\mathrm{Tr}\lvert YXY\rvert`$ | $`X=X^*`$，$`Y\ge0`$ |

其中 $`|X|`$、$`|Y|`$ 归一化为 $`\lambda_{\max}=1`$，它们的对比度作为参数。**包络**为
$`\overline z(\delta)=\sup\{z:\Delta(|X|)=\delta\}`$ 与 $`\underline z(\delta)=\inf\{z:\Delta(|X|)=\delta\}`$，其余
变量均自由变化。

<a id="s-11-1"></a>

## 11.1 约束集及其参数化

**定义。** 对 $`\delta\in[0,1)`$ 与 $`\ell=\ell(\delta)`$，**对比度壳层**为
```math
\mathcal P_\delta=\{X\in\mathcal M:X\ge0,\ \lambda_{\max}(X)=1,\ \Delta(X)=\delta\}
=\{X=X^*:\mathrm{spec}\,X\subseteq[\ell,1],\ 1,\ell\in\mathrm{spec}\,X\}.
```
$`\mathcal M`$ 的**槽位**是 $`D`$ 个数对 $`(c,i)`$，$`i\le k_c`$；对 $`\lambda\in\mathbb R^D`$，
$`\mathrm{diag}(\lambda)=\bigoplus_c\mathrm{diag}(\lambda_{c,1},\dots,\lambda_{c,k_c})`$。

**参数化。** 对给定的对比度 $`\delta`$，搜索对每个起点抽取两个不同的槽位
$`s_\pm=(c_\pm,i_\pm)`$（先在 $`D`$ 个槽位中均匀地抽取 $`s_+`$，再在其余 $`D-1`$ 个槽位中均匀地抽取 $`s_-`$，因此
每个由不同槽位构成的有序对的概率都是 $`1/(D(D-1))`$）和一个 Haar 随机酉算子 $`U_0\in\mathcal M`$
（引理 11.10），并把 $`\theta\in\mathbb R^D`$ 与
$`G=\bigoplus_cG_c`$ 映为
```math
X(\theta,G)=U\mathrm{diag}(\lambda)U^*,\qquad U=U_0\,e^{G-G^*},\qquad
\lambda_{s_+}=1,\quad\lambda_{s_-}=\ell,\quad\lambda_s=\ell+(1-\ell)\tfrac{1+\sin\theta_s}2\ \text{otherwise}.
```
单独的酉算子则由 $`G\mapsto U_0e^{G-G^*}`$ 参数化。与 sigmoid 不同，$`(1+\sin\theta)/2`$ 在有限的 $`\theta`$ 处
就能取到 $`0`$ 和 $`1`$，这一点很方便，因为本章的极值谱（推论 11.7、命题 11.14）位于 $`\ell`$ 和 $`1`$。
$`[a,b]`$ 中的自由对比度只是多一个参数，$`\delta=a+(b-a)(1+\sin\eta)/2`$，例如取 $`[a,b]=[0,0.999]`$。在这一限制下取到的值在没有限制时也能取到，因此它仍然是无限制上确界的下界
（无限制下确界的上界）。

*库中：* 这一参数化就是 `optimize.PositiveParam`，而 $`G\mapsto U_0e^{G-G^*}`$ 就是 `optimize.UnitaryParam`——[教程，11.1 节](../../tutorial/zh/11_optimisation.md#s-11-1)。

<a id="r-11-1"></a>

**引理 11.1.** 对每个酉算子 $`U_0\in\mathcal M`$，$`G\mapsto U_0e^{G-G^*}`$ 映满 $`\mathcal M`$ 的
酉群。

*证明。* 设 $`U`$ 是酉算子。由定理 1.3（逐块应用）与推论 1.4(2)，
$`U_0^*U=\bigoplus_cW_c\mathrm{diag}(e^{i\varphi_{c,j}})_jW_c^*`$，其中 $`W_c`$ 为酉矩阵，$`\varphi_{c,j}\in\mathbb R`$。于是
$`S=\bigoplus_cW_c\mathrm{diag}(i\varphi_{c,j})_jW_c^*`$ 满足 $`S^*=-S`$ 与 $`e^S=U_0^*U`$（级数
$`\sum_nS^n/n!`$ 在特征基中逐项计算），而 $`G=S/2`$ 给出 $`G-G^*=S`$。$`\square`$

<a id="r-11-2"></a>

**命题 11.2（参数化是满射）。** 对固定的 $`\delta`$、$`s_\pm`$ 与 $`U_0`$，
$`(\theta,G)\mapsto X(\theta,G)`$ 的像是 $`\{X\in\mathcal P_\delta:1\in\mathrm{spec}\,X_{c_+},\ \ell\in\mathrm{spec}\,X_{c_-}\}`$，
并且 $`\mathcal P_\delta`$ 是这些像对所有不同槽位对的并。

*证明。* 像包含在该集合中：$`X`$ 的特征值 $`\lambda_s\in[\ell,1]`$，且 $`\lambda_{s_\pm}`$ 位于块
$`c_\pm`$ 中。反之，对角化 $`X_c=W_c\mathrm{diag}(\mu_c)W_c^*`$ 并对每个 $`\mu_c`$ 作置换（把 $`W_c`$ 换成
$`W_c\Pi_c`$，$`\Pi_c`$ 为置换矩阵），使 $`\mu_{s_+}=1`$ 且 $`\mu_{s_-}=\ell`$。这是可行的：若 $`c_+=c_-`$
且 $`\delta>0`$，则 $`1\ne\ell`$ 是同一块的两个特征值，占据不同的位置；若 $`\delta=0`$，
所有特征值都是 $`1`$。若 $`\ell<1`$，对其余槽位令 $`\theta_s=\arcsin\big(2(\mu_s-\ell)/(1-\ell)-1\big)`$，
使 $`\lambda=\mu`$（若 $`\ell=1`$，$`\theta`$ 任意）。引理 11.1 给出满足 $`U_0e^{G-G^*}=\bigoplus_cW_c`$ 的 $`G`$。
*并：* 对 $`X\in\mathcal P_\delta`$，选取含特征值 $`1`$、$`\ell`$ 的块 $`c_\pm`$，并在其中选取不同的槽位；
若 $`c_+=c_-`$ 且 $`\delta>0`$，该块有两个不同的特征值，所以 $`k_{c_+}\ge2`$；若 $`\delta=0`$，则 $`X=1`$，
任取两个槽位即可。$`\square`$

*计算：* 构造性地检验命题 11.2：证明中的参数重现一个目标算子——[教程，11.1 节](../../tutorial/zh/11_optimisation.md#s-11-1)。

由于每一对槽位都以正概率被抽到，壳层中的每一点都以正概率位于某个起点的像中，并且（在精确算术中）
约束在每一步迭代中都成立。
对自伴的 $`X`$（实验 3），使用 $`X=U\mathrm{diag}(\sigma\lambda)U^*`$，其中符号 $`\sigma\in\{\pm1\}^D`$ 对每个起点抽取
（在每个块中，负号的个数在 $`\{0,\dots,k_c\}`$ 中均匀分布，其位置均匀随机，因此每种符号模式都有正概率），
并且只通过 11.5 节中的符号翻转改变。于是
$`|X|=U\mathrm{diag}(\lambda)U^*`$（定理 1.5(e)，因为 $`\lambda>0`$），并且对所有符号与槽位取并得到
$`\{X=X^*:|X|\in\mathcal P_\delta\}`$：对这样的 $`X`$，对角化 $`X_c=W_c\mathrm{diag}(\nu_c)W_c^*`$，于是
$`|X|_c=W_c\mathrm{diag}(|\nu_c|)W_c^*`$，其中 $`|\nu_{c,i}|\ge\ell>0`$，再对 $`\mu=|\nu|`$ 与
$`\sigma_s=\nu_s/|\nu_s|`$ 应用命题 11.2 的证明。只有 $`S=G-G^*`$ 起作用。

*库中：* `optimize.SelfAdjointParam`——[教程，11.1 节](../../tutorial/zh/11_optimisation.md#s-11-1)。

坐标卡 $`S\mapsto U_0e^S`$ 只在 $`S=0`$ 附近性质良好。令 $`\mathfrak u=\{S\in\mathcal M:S^*=-S\}`$，
$`\mathcal M_h=\{H\in\mathcal M:H^*=H\}`$，它们是实子空间，且 $`\mathcal M=\mathfrak u\oplus\mathcal M_h`$
（$`X=\frac{X-X^*}2+\frac{X+X^*}2`$）。*在 $`0`$ 附近，$`S\mapsto U_0e^S`$ 是 $`\mathfrak u`$ 的一个开球 $`N\ni0`$ 到
酉群中 $`U_0`$ 的一个邻域上的同胚。* 事实上，$`\Psi(S,H)=U_0e^Se^H`$ 是 $`C^1`$ 映射
$`\mathfrak u\times\mathcal M_h\to\mathcal M`$（级数 $`\sum_{m,n\ge0}U_0S^mH^n/(m!\,n!)`$ 以及其各项关于实坐标的
偏导数构成的级数在有界集上一致收敛），而
$`e^{tS}e^{tH}=1+t(S+H)+O(t^2)`$ 表明它在 $`0`$ 处的微分为 $`(S,H)\mapsto U_0(S+H)`$，这是到
$`\mathcal M`$ 上的双射。由反函数定理，存在 $`\mathfrak u`$ 中的开球 $`N\ni0`$ 与 $`\mathcal M_h`$ 中的开球 $`N_h\ni0`$，
使得 $`\Psi`$ 是 $`N\times N_h`$ 到开集 $`O\subseteq\mathcal M`$ 上的同胚。若
$`\Psi(S,H)`$ 是酉的，则 $`e^H=e^{-S}U_0^*\Psi(S,H)`$ 也是酉的（在正规矩阵 $`S`$ 的特征基中——其特征值为纯虚数——
$`e^{\pm S}`$ 是对角的，对角元的模为 $`1`$，并且二者互逆）；但 $`e^H`$ 具有正特征值 $`e^{h_i}`$，$`h_i`$ 为 $`H`$ 的特征值，
所以 $`e^{h_i}=1`$（推论 1.4(2)），从而 $`H=0`$。
因此 $`O\cap U(\mathcal M)=U_0e^N`$，其中 $`U(\mathcal M)`$ 是 $`\mathcal M`$ 的酉群：这是 $`U_0`$ 在
$`U(\mathcal M)`$ 中的一个邻域，而 $`S\mapsto U_0e^S`$ 是连续双射 $`N\to U_0e^N`$，其逆映射，即 $`\Psi^{-1}`$ 的
第一个分量，是连续的。远离 $`0`$ 时该坐标卡不是单射：对每个投影 $`P\in\mathcal M`$ 有 $`e^{2\pi iP}=1=e^0`$。
因此搜索在每一轮之后都把坐标卡**重新居中**（$`U_0\leftarrow U_0e^{G-G^*}`$，$`G\leftarrow0`$）。另外两个坐标卡
在 $`0`$ 处具有相同的微分 $`S\mapsto U_0S`$。Cayley 坐标卡
$`S\mapsto U_0C(-S/2)`$（习题 1）：$`(1+tS/2)(1-tS/2)^{-1}=(1+tS/2)(1+tS/2+O(t^2))=1+tS+O(t^2)`$。极分解回缩，
即把 $`S`$ 映为 $`A=U_0(1+S)`$ 的酉因子（定理 1.8）：$`A^*A=(1-S)(1+S)=1+S^*S\ge1`$ 可逆，所以 $`|A|`$
可逆，且酉因子为 $`A|A|^{-1}`$；把 $`S`$ 换成 $`tS`$，有 $`|A|^{-1}=(1+t^2S^*S)^{-1/2}=1+O(t^2)`$
（在 $`S^*S`$ 的特征基中计算），所以该酉因子为 $`U_0(1+tS)+O(t^2)`$。

> **物理意义**（物理解释，而非本课程的定理）。Michelson 对比度就是干涉条纹的可见度 $`(I_{\max}-I_{\min})/(I_{\max}+I_{\min})`$。
> 对 $`\lambda_{\max}=1`$ 的 $`0\le X\le1`$——滤光片、衰减器、探测器效率——
> 它确定了动态范围 $`\lambda_{\min}/\lambda_{\max}`$，即条件数的倒数。

<a id="s-11-2"></a>

## 11.2 冯·诺依曼迹不等式

元素 $`c_{ij}\ge0`$ 且所有行和与列和都为 $`1`$ 的矩阵 $`(c_{ij})`$ 称为**双随机**矩阵；对酉矩阵 $`G`$，
$`(|G_{ij}|^2)`$ 与 $`(|G_{ji}|^2)`$ 都是双随机的，因为 $`G`$ 的行和列都是单位向量。

<a id="r-11-3"></a>

**引理 11.3（重排）。** 设 $`p_1\ge\dots\ge p_k`$，$`q_1\ge\dots\ge q_k`$ 为实数，$`(c_{ij})`$ 为双
随机矩阵。则 $`\sum_ip_iq_{k+1-i}\le\sum_{i,j}p_iq_jc_{ij}\le\sum_ip_iq_i`$。

*证明。* *上界。* 对 $`m\le k`$ 令 $`r_j=\sum_{i\le m}c_{ij}\in[0,1]`$，于是 $`\sum_jr_j=m`$。由于当
$`j\le m`$ 时 $`q_j\ge q_m`$，当 $`j>m`$ 时 $`q_j\le q_m`$，
```math
\sum_{j\le m}q_j-\sum_jq_jr_j=\sum_{j\le m}q_j(1-r_j)-\sum_{j>m}q_jr_j\ge q_m\Big(m-\sum_jr_j\Big)=0 .
```
记 $`x_i=\sum_jq_jc_{ij}`$，上式即 $`\sum_{i\le m}x_i\le\sum_{j\le m}q_j`$，且 $`m=k`$ 时取等号。由 Abel
求和以及 $`p_m-p_{m+1}\ge0`$，
```math
\sum_ip_ix_i=p_k\sum_{i\le k}x_i+\sum_{m<k}(p_m-p_{m+1})\sum_{i\le m}x_i
\le p_k\sum_{j\le k}q_j+\sum_{m<k}(p_m-p_{m+1})\sum_{j\le m}q_j=\sum_ip_iq_i .
```
*下界。* 对 $`p`$、递减的 $`q'_j=-q_{k+1-j}`$ 以及双随机矩阵
$`c'_{ij}=c_{i,k+1-j}`$ 应用上界：$`-\sum_{ij}p_iq_jc_{ij}=\sum_{ij}p_iq'_jc'_{ij}\le\sum_ip_iq'_i=-\sum_ip_iq_{k+1-i}`$。$`\square`$

每个 $`A\in M_k(\mathbb C)`$ 都有**奇异值分解** $`A=W\Sigma V^*`$，$`W,V`$ 为酉矩阵，
$`\Sigma=\mathrm{diag}(s_i(A))`$：写 $`A=W'|A|`$（定理 1.8），$`|A|=V\Sigma V^*`$（定理 1.3），$`W=W'V`$。

<a id="r-11-4"></a>

**定理 11.4（冯·诺依曼迹不等式）。** 对 $`A,B\in M_k(\mathbb C)`$，
$`|\mathrm{Tr}(AB)|\le\sum_is_i(A)\,s_i(B)`$。

*证明。* 设 $`A=W_A\Sigma_AV_A^*`$，$`B=W_B\Sigma_BV_B^*`$，$`a_i=s_i(A)`$，$`b_j=s_j(B)`$，$`G=V_A^*W_B`$，$`H=V_B^*W_A`$。
由循环性，$`\mathrm{Tr}(AB)=\mathrm{Tr}(\Sigma_AG\Sigma_BH)=\sum_{i,j}a_iG_{ij}b_jH_{ji}`$，因此由
$`|g||h|\le\frac12(|g|^2+|h|^2)`$ 与 $`a_ib_j\ge0`$，
```math
|\mathrm{Tr}(AB)|\le\sum_{i,j}a_ib_j\,d_{ij},\qquad d_{ij}=\tfrac12\big(|G_{ij}|^2+|H_{ji}|^2\big).
```
$`(d_{ij})`$ 是双随机的，引理 11.3 给出 $`\sum_{ij}a_ib_jd_{ij}\le\sum_ia_ib_i`$。$`\square`$

<a id="r-11-5"></a>

**推论 11.5.** (a) 对特征值为 $`p,q`$ 的自伴矩阵 $`P,Q\in M_k(\mathbb C)`$：
$`\sum_ip^\downarrow_iq^\uparrow_i\le\mathrm{Tr}(PQ)\le\sum_ip^\downarrow_iq^\downarrow_i`$。
(b) 对 $`P,Q\ge0`$：$`\mathrm{Tr}(PQ)\le\|PQ\|_1\le\sum_ip^\downarrow_iq^\downarrow_i`$。

*证明。* (a) 与定理 3.5 的证明一样，$`\mathrm{Tr}(PQ)=\sum_{ij}p_iq_jc_{ij}`$，其中对按递减顺序排列的特征基，
$`c_{ij}=|\langle u_i,v_j\rangle|^2`$ 是双随机的；应用引理 11.3（其下界为
$`\sum_ip^\downarrow_iq^\downarrow_{k+1-i}=\sum_ip^\downarrow_iq^\uparrow_i`$）。(b) 第一个不等式就是
定理 1.9。由定理 1.9 与定理 11.4，$`\|PQ\|_1=\max_U|\mathrm{Tr}(P\,QU)|\le\sum_is_i(P)s_i(QU)`$，而
$`s_i(QU)=q^\downarrow_i`$，因为 $`|QU|=(U^*Q^2U)^{1/2}=U^*QU`$（定理 1.6）。$`\square`$

*计算：* 在随机矩阵上检验定理 11.4 与推论 11.5——[教程，11.2 节](../../tutorial/zh/11_optimisation.md#s-11-2)。

因此对 $`M_k(\mathbb C)`$ 中谱为 $`x,y`$ 的 $`X,Y\ge0`$，命题 3.9 中的违背量 $`\|YX\|_1-\mathrm{Tr}(XY)`$
至多为 $`\sum_ix^\downarrow_i(y^\downarrow_i-y^\uparrow_i)`$（对 $`(Y,X)`$ 应用推论 11.5(b)，再用 (a) 的下界）。
对实验 2，这个数恰好就是精确答案：

<a id="r-11-6"></a>

**命题 11.6（固定谱时的实验 2）。** 在实验 2 中，$`|X|=X_0`$，$`|Y|=Y_0`$，且
$`z=\|X_0VY_0\|_1-\mathrm{Tr}(X_0Y_0)`$。若 $`(X_0)_c`$ 的谱 $`x_c`$ 与 $`(Y_0)_c`$ 的谱 $`y_c`$ 固定，而
特征基与 $`V`$ 变化，则
```math
\sup z=-\inf z=\sum_c\sum_ix^\downarrow_{c,i}\big(y^\downarrow_{c,i}-y^\uparrow_{c,i}\big),
```
并且二者都能取到。

*证明。* $`X^*X=X_0U^*UX_0=X_0^2`$，所以 $`|X|=X_0`$（定理 1.6）；同理 $`|Y|=Y_0`$，而 $`|UA|=|A|`$ 给出
$`\mathrm{Tr}|XY|=\|X_0VY_0\|_1`$。所有变量都按块作用，且 $`z`$ 是各块之和，故不妨设 $`C=1`$。
令 $`Y'=VY_0V^*\ge0`$（谱为 $`y`$）；由于 $`|AV|=V^*|A|V`$，$`\|X_0VY_0\|_1=\|X_0Y'V\|_1=\|X_0Y'\|_1`$。对数对
$`(X_0,Y')`$ 与 $`(X_0,Y_0)`$ 应用推论 11.5，并利用 $`\|X_0Y'\|_1\ge\mathrm{Tr}(X_0Y')`$，得到
```math
\sum_ix^\downarrow_iy^\uparrow_i-\sum_ix^\downarrow_iy^\downarrow_i\le\|X_0Y'\|_1-\mathrm{Tr}(X_0Y_0)\le\sum_ix^\downarrow_iy^\downarrow_i-\sum_ix^\downarrow_iy^\uparrow_i .
```
设 $`J`$ 为反序置换，$`Je_i=e_{k+1-i}`$；则 $`J\mathrm{diag}(w)=\mathrm{diag}(w_k,\dots,w_1)J`$，并且
$`\mathrm{diag}(w)J`$ 的迹范数为 $`\sum_i|w_i|`$。取 $`X_0=\mathrm{diag}(x^\downarrow)`$ 与 $`V=J`$：
对 $`Y_0=\mathrm{diag}(y^\uparrow)`$，$`X_0VY_0=\mathrm{diag}(x^\downarrow_iy^\downarrow_i)J`$ 且
$`\mathrm{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\uparrow_i`$（上界）；对 $`Y_0=\mathrm{diag}(y^\downarrow)`$，
$`X_0VY_0=\mathrm{diag}(x^\downarrow_iy^\uparrow_i)J`$ 且 $`\mathrm{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\downarrow_i`$
（下界）。$`\square`$

<a id="r-11-7"></a>

**推论 11.7（实验 2 的包络）。** 设 $`\ell_X,\ell_Y\in(0,1]`$，且某个 $`k_c\ge2`$。在所有满足
$`\lambda_{\max}=1`$、$`\lambda_{\min}(X_0)=\ell_X`$、$`\lambda_{\min}(Y_0)=\ell_Y`$ 的 $`X_0,Y_0\ge0`$ 以及所有酉算子上，
```math
\sup z=-\inf z=\sum_c\Big\lfloor\frac{k_c}2\Big\rfloor(1-\ell_X)(1-\ell_Y).
```

*证明。* 在一个块中记 $`x_i=x^\downarrow_{c,i}`$，$`y_i=y^\downarrow_{c,i}`$，$`k=k_c`$，于是 $`y^\uparrow_{c,i}=y_{k+1-i}`$。
把 $`i`$ 与 $`k+1-i`$ 配对，$`x_i(y_i-y_{k+1-i})+x_{k+1-i}(y_{k+1-i}-y_i)=(x_i-x_{k+1-i})(y_i-y_{k+1-i})`$，而当
$`k`$ 为奇数时，中间的指标贡献 $`x_i(y_i-y_i)=0`$。因此
```math
\sum_ix_i(y_i-y_{k+1-i})=\sum_{i\le\lfloor k/2\rfloor}(x_i-x_{k+1-i})(y_i-y_{k+1-i})\le\Big\lfloor\frac k2\Big\rfloor(1-\ell_X)(1-\ell_Y),
```
因为 $`x_i\in[\ell_X,1]`$，$`y_i\in[\ell_Y,1]`$ 都是递减的，所以两个因子分别位于 $`[0,1-\ell_X]`$ 与 $`[0,1-\ell_Y]`$ 中。
若在每个块中，最大的 $`\lfloor k_c/2\rfloor`$ 个特征值为 $`1`$，最小的 $`\lfloor k_c/2\rfloor`$ 个为 $`\ell_X`$（相应地为
$`\ell_Y`$），中间那个（$`k_c`$ 为奇数时）在 $`[\ell_X,1]`$（相应地 $`[\ell_Y,1]`$）中任意，则等号成立；由于某个 $`k_c\ge2`$，这些谱的全局极值为 $`1`$ 和 $`\ell_X`$（$`\ell_Y`$）。由命题 11.6
即得结论。$`\square`$

这里对 $`\Delta(|X|)=\delta`$ 有 $`1-\ell_X=2\delta/(1+\delta)`$；当 $`k_c`$ 为奇数时，中间的特征值无关紧要。

*计算：* $`k=3`$、$`C=2`$ 时命题 11.6 的极值构型与推论 11.7 的包络——[教程，11.2 节](../../tutorial/zh/11_optimisation.md#s-11-2)；多起点搜索求得的包络——[教程，11.5 节](../../tutorial/zh/11_optimisation.md#s-11-5)。

<a id="s-11-3"></a>

## 11.3 一个带有取等条件的迹不等式

对 $`X=X^*`$，令 $`X_\pm=f_\pm(X)`$，其中 $`f_+(t)=\max(t,0)`$，$`f_-(t)=\max(-t,0)`$。则 $`X_\pm\ge0`$，$`X=X_+-X_-`$，
$`|X|=X_++X_-`$，且 $`X_+X_-=0`$（定理 1.5(c)；对 $`|X|=(X^2)^{1/2}=h(X)`$，其中
$`h(t)=\sqrt{t^2}=|t|=f_+(t)+f_-(t)`$，用定理 1.5(e)；正性由推论 1.4(3) 得到）。

<a id="r-11-8"></a>

**引理 11.8.** 对 $`M_n(\mathbb C)`$ 中的 $`A,B\ge0`$，$`\|A-B\|_1\le\mathrm{Tr}\,A+\mathrm{Tr}\,B`$，等号成立当且仅当
$`AB=0`$。

*证明。* 令 $`H=A-B`$，$`Q=g(H)`$，其中 $`g`$ 在 $`(0,\infty)`$ 上为 $`1`$，在其他地方为 $`0`$；$`Q`$ 是投影，因为 $`g=\bar g=g^2`$
（定理 1.5(c)）。由于 $`|t|=t(2g(t)-1)`$，有 $`|H|=H(2Q-1)`$（定理 1.5(c)、(e)），再对投影 $`R=Q,1-Q`$ 利用
$`\mathrm{Tr}(AR)=\mathrm{Tr}(RAR)`$（循环性），得到
```math
\|H\|_1=\mathrm{Tr}\big(A(2Q-1)\big)-\mathrm{Tr}\big(B(2Q-1)\big)
=\mathrm{Tr}\,A+\mathrm{Tr}\,B-2\mathrm{Tr}\big((1-Q)A(1-Q)\big)-2\mathrm{Tr}(QBQ).
```
被减去的迹都 $`\ge0`$，而 $`\mathrm{Tr}(RAR)=\|A^{1/2}R\|_2^2`$ 为零当且仅当 $`AR=0`$；所以等号成立
当且仅当 $`A(1-Q)=0`$ 且 $`BQ=0`$，此时 $`AB=AQB=A(BQ)^*=0`$。反之，若 $`AB=0`$，则 $`BA=0`$，
$`(A-B)^2=A^2+B^2=(A+B)^2`$，由定理 1.6 得 $`|A-B|=A+B`$。$`\square`$

<a id="r-11-9"></a>

**命题 11.9（实验 3）。** 对 $`\mathcal M`$ 中的 $`X=X^*`$ 与 $`Y\ge0`$，
$`\mathrm{Tr}|YXY|\le\mathrm{Tr}(Y|X|Y)`$，等号成立当且仅当 $`YX_+Y^2X_-Y=0`$，特别地，当 $`XY=YX`$ 时成立。

*证明。* 在 $`M_D(\mathbb C)`$ 中，$`A=YX_+Y`$ 与 $`B=YX_-Y`$ 是正的（取 $`T=X_\pm^{1/2}Y`$ 应用推论 1.7），$`YXY=A-B`$，
且 $`\mathrm{Tr}(Y|X|Y)=\mathrm{Tr}\,A+\mathrm{Tr}\,B`$；应用引理 11.8。若 $`XY=YX`$，则 $`Y`$ 与
$`X_\pm`$ 可交换（定理 1.5(d)），从而 $`AB=Y^4X_+X_-=0`$。$`\square`$

因此在实验 3 中 $`z\ge0`$，且 $`\underline z=0`$（取 $`Y=1`$）。取等条件比可交换性更弱：对
$`X=\mathrm{diag}(1,-1,0)`$ 以及行为 $`(1,0,a)`$、$`(0,1,0)`$、$`(a,0,1)`$（$`-1<a<1`$）的 $`Y`$，
$`AB=(Ye_1)\langle Ye_1,Ye_2\rangle(Ye_2)^*=0`$，而当 $`a\ne0`$ 时 $`XY\ne YX`$。

*计算：* 在随机矩阵对上检验命题 11.9，以及不可交换而等号成立的例子——[教程，11.3 节](../../tutorial/zh/11_optimisation.md#s-11-3)。

<a id="s-11-4"></a>

## 11.4 采样与优化

与其最大值 $`\|A\|_1`$ 相比，$`|\mathrm{Tr}(AU)|`$ 的典型大小是多少？我们把 $`U`$ 取为
**Ginibre 矩阵** $`Z`$（元素相互独立，密度为
$`\pi^{-1}e^{-|z|^2}`$；Mezzadri 的方法）的 QR 分解 $`Z=QR`$（$`R_{ii}>0`$）中的因子 $`Q`$，在每个块中独立抽取。

<a id="r-11-10"></a>

**引理 11.10.** $`Q`$ 在 $`U(k)`$ 上服从 Haar 分布；特别地，对每个固定的酉矩阵 $`V\in M_k(\mathbb C)`$，$`VQ`$ 与
$`Q`$ 同分布。各块相互独立的酉算子 $`\bigoplus_cQ_c`$ 在 $`\mathcal M`$ 的酉群上服从 Haar 分布。

*证明。* 这就是定理 5.3（左不变性是其证明中的等变性步骤；该分布由唯一性引理 5.1 认定为 Haar 测度，
其测度论步骤是推论 A.2）与推论 5.4。$`\square`$

*库中：* 这一采样器就是 `TypeIAlgebra.random_unitary_operator`——[教程，11.4 节](../../tutorial/zh/11_optimisation.md#s-11-4)。

<a id="r-11-11"></a>

**命题 11.11.** 对 Haar 分布的 $`U\in U(k)`$ 与 $`A\in M_k(\mathbb C)`$：$`\mathbb E\mathrm{Tr}(AU)=0`$，且
$`\mathbb E|\mathrm{Tr}(AU)|^2=\|A\|_2^2/k`$。

*证明。* 由推论 5.5，$`\mathbb EU_{ij}=0`$。对 $`X=|e_j\rangle\langle e_m|`$ 应用命题 5.7（取 $`C=1`$），得到
一阶 Weingarten 公式
$`\mathbb E[U_{ij}\bar U_{lm}]=\mathbb E(UXU^*)_{il}=\delta_{il}\delta_{jm}/k`$；对 $`k\ge2`$，它也可由
二阶公式（定理 5.10）得到：令 $`i_2=i'_2`$，$`j_2=j'_2`$ 并对二者求和，因为
$`\sum_{i_2,j_2}|U_{i_2j_2}|^2=k`$。利用 $`\mathrm{Tr}(AU)=\sum_{ij}A_{ji}U_{ij}`$，
$`\mathbb E|\mathrm{Tr}(AU)|^2=\sum_{ijlm}A_{ji}\bar A_{ml}\mathbb E[U_{ij}\bar U_{lm}]=\frac1k\sum_{ij}|A_{ji}|^2`$。$`\square`$

<a id="r-11-12"></a>

**推论 11.12.** 设 $`M_k(\mathbb C)`$ 中的 $`X,Y\ge0`$ 的谱含于 $`[a,b]`$，$`a>0`$，$`U`$ 服从 Haar 分布。则
```math
\frac{\big(\mathbb E|\mathrm{Tr}(XUY)|^2\big)^{1/2}}{\max_V|\mathrm{Tr}(XVY)|}=\frac{\|YX\|_2}{\sqrt k\,\|YX\|_1}\le\min\Big(\frac1{\sqrt k},\frac{b^2}{a^2k}\Big),
\qquad
\mathbb P\big(|\mathrm{Tr}(XUY)|\ge\mathrm{Tr}(XY)\big)\le\frac{b^4}{a^4k^2}.
```

*证明。* 其中的等式由命题 11.11（取 $`A=YX`$，因为 $`\mathrm{Tr}(XUY)=\mathrm{Tr}(AU)`$）与命题 3.9 得到。
由于 $`\|A\|_2^2=\sum_is_i^2\le(\sum_is_i)^2=\|A\|_1^2`$（1.4 节），该比值 $`\le1/\sqrt k`$。
由于 $`\|Av\|=\||A|v\|`$（1.4 节），对角化 $`|A|`$ 给出 $`s_k(A)=\min_{\|v\|=1}\|Av\|`$ 与
$`s_1(A)=\max_{\|v\|=1}\|Av\|`$；而对单位向量 $`v`$ 有 $`\|Xv\|^2=\langle v,X^2v\rangle\in[a^2,b^2]`$（对 $`X^2-a^2`$ 与
$`b^2-X^2`$ 应用推论 1.4(3)），对 $`Y`$ 同理。因此 $`a^2\le s_k(A)\le s_1(A)\le b^2`$，$`\|A\|_2^2\le kb^4`$，且
$`\|A\|_1\ge ka^2`$。最后 $`\mathrm{Tr}(XY)=\mathrm{Tr}(X^{1/2}YX^{1/2})\ge a\mathrm{Tr}\,X\ge ka^2`$，而对
$`|\mathrm{Tr}(XUY)|^2`$ 应用 Markov 不等式，可把该概率界定为 $`(kb^4/k)/(ka^2)^2`$。$`\square`$

因此，$`N`$ 个 Haar 样本找到 $`|\mathrm{Tr}(XUY)|\le\mathrm{Tr}(XY)`$ 的反例的概率至多为
$`Nb^4/(a^4k^2)`$（$`N`$ 个事件的并），尽管对每一对不可交换的矩阵这样的反例都存在（命题 3.9）。这个界
只有在 $`k>\sqrt N\,b^2/a^2`$ 时才有信息量（对 $`b/a=6`$ 与 $`N=2000`$：$`k>1609`$）。

*计算：* 用 Monte Carlo 检验命题 11.11，并对 $`k\le32`$ 把推论 11.12 中的比值与 $`2000`$ 个样本中的最好者比较——[教程，11.4 节](../../tutorial/zh/11_optimisation.md#s-11-4)。

<a id="r-11-13"></a>

**引理 11.13.** 若 $`k\ge2`$，则对每个 $`A\in M_k(\mathbb C)`$ 有 $`\min_U|\mathrm{Tr}(AU)|=0`$。

*证明。* 设 $`A=W\Sigma V^*`$，$`P`$ 为循环移位 $`Pe_i=e_{(i\bmod k)+1}`$，当 $`k\ge2`$ 时其对角线为零。对
$`U=VPW^*`$，$`\mathrm{Tr}(AU)=\mathrm{Tr}(\Sigma P)=0`$。$`\square`$

所以 $`\inf_Uz=-\mathrm{Tr}(XY)`$，$`\sup_Uz=\|YX\|_1-\mathrm{Tr}(XY)`$（命题 3.9）。由命题 11.11，$`|\mathrm{Tr}(XUY)|`$ 的均方根为 $`\|YX\|_2/\sqrt k`$；利用推论 11.12 的证明中的
界（谱含于 $`[a,b]`$），$`\mathbb E|\mathrm{Tr}(XUY)|^2\le b^4`$ 且 $`\mathrm{Tr}(XY)\ge ka^2`$，所以 Markov 不等式给出
$`\mathbb P\big(z\ge-(1-\varepsilon)\mathrm{Tr}(XY)\big)\le b^4/(\varepsilon^2a^4k^2)`$（$`\varepsilon>0`$）：对大的 $`k`$，
$`z`$ 以高概率接近其下确界。

*计算：* 用梯度上升与梯度下降求关于 $`U`$ 的两个极值，并与 $`20\,000`$ 个 Haar 样本比较——[教程，11.4 节](../../tutorial/zh/11_optimisation.md#s-11-4)。

> **物理意义**（物理解释，而非本课程的定理）。推论 11.12 是统计力学所称的**典型性**的一个较弱但可证明的
> 实例：在高维中，人们预期几乎每个态或演化都给出宏观量几乎相同的值。一个精确（最优）的不等式由精细调节的
> 构型决定——对齐的特征基、两能级的谱——而随机酉算子（随机线路、随机散射矩阵）很少产生这样的构型。

<a id="s-11-5"></a>

## 11.5 包络的多起点搜索

**多起点搜索**从许多独立的起点出发，分轮运行梯度方法，每个起点有自己的槽位、符号、$`U_0`$ 与对比度；
每轮之后对酉坐标卡重新居中，并对每个起点返回目标函数在各迭代点处的最优值。自伴参数化的符号是离散的：
在两轮之间，通过在随机槽位上的单个符号翻转来改变它们，每次翻转若能改进目标函数则予以保留。

*库中：* `optimize.extremize`——[教程，11.5 节](../../tutorial/zh/11_optimisation.md#s-11-5)。

在精确算术中，每一步迭代都是容许的（命题 11.2 及其后的说明），所以返回的每个值都是 $`z`$ 在某个容许点处的值；
在浮点运算中，它在舍入误差范围内等于该值。因此最大化给出上确界的
**下界**，最小化给出下确界的上界，仅此而已——局部搜索可能停在局部极值处。另一个方向需要证明，
如推论 11.7 与命题 11.9，或者凸松弛（11.6 节）。计算得到的包络要接受以下检验：(i) 与闭式表达比较；(ii) 与
手工构造的构型比较，例如较小解的直和，包络不得低于它们；(iii) 与
结构化样本（两能级谱、置换基）比较，它们的值是另外一些可取到的值。

<a id="r-11-14"></a>

**命题 11.14（实验 1 的下包络）。** 设所有 $`k_c\ge2`$。在 $`X\in\mathcal P_\delta`$、满足
$`\lambda_{\max}(Y)=1`$ 的 $`Y\ge0`$ 以及酉算子 $`U\in\mathcal M`$ 上，$`\ \min z=-\big(D-2\delta/(1+\delta)\big)`$。

*证明。* $`Y\le1`$ 给出 $`X^{1/2}YX^{1/2}\le X`$，所以 $`\mathrm{Tr}(XY)\le\mathrm{Tr}\,X\le(D-1)+\ell(\delta)`$，因为
$`X`$ 的一个特征值为 $`\ell(\delta)`$，其余都 $`\le1`$；因此 $`z\ge-\mathrm{Tr}(XY)\ge-(D-2\delta/(1+\delta))`$。
当 $`X`$ 的特征值为一个 $`\ell(\delta)`$、其余都为 $`1`$，$`Y=1`$，且 $`U=\bigoplus_cU_c`$ 满足
$`\mathrm{Tr}(X_cU_c)=0`$（引理 11.13）时，等号成立。$`\square`$

类型 (ii) 的检验可以利用直和：实验 3 的 $`z`$ 对直和是可加的，而直和的对比度是各直和项对比度中
最大者（习题 3），所以对于实验 3 的包络（$`\Delta(Y)`$ 在固定区间 $`[0,b]`$ 中自由变化），
```math
\overline z_{M_2\oplus M_2}(\delta)\ \ge\ \overline z_{M_2}(\delta)+\sup_{\delta'\le\delta}\overline z_{M_2}(\delta');
```
特别地，在 $`M_2`$ 上于对比度 $`\delta`$ 与 $`\delta'\le\delta`$ 处取到的两个值之和，是在 $`M_2\oplus M_2`$ 上于对比度
$`\delta`$ 处取到的值。

*计算：* 用多起点搜索检验 $`M_2\oplus M_2`$ 上的命题 11.14、$`M_3\oplus M_3`$ 上的推论 11.7 以及命题 11.9，并做一项类型 (ii) 的检验——[教程，11.5 节](../../tutorial/zh/11_optimisation.md#s-11-5)。

对大小为 $`2`$ 的块，实验 1 与实验 3 的上包络可以精确求出。在 $`M_2(\mathbb C)`$ 中，本征基只通过一个数起作用：
对正交归一本征基 $`e_1,e_2`$ 与 $`f_1,f_2`$，双随机矩阵 $`(|\langle e_i,f_j\rangle|^2)`$ 的两行为
$`(\gamma,1-\gamma)`$ 与 $`(1-\gamma,\gamma)`$，$`\gamma=|\langle e_1,f_1\rangle|^2\in[0,1]`$，而且每个 $`\gamma`$
都可由实基取到（$`f_1=\cos t\,e_1+\sin t\,e_2`$，$`f_2=-\sin t\,e_1+\cos t\,e_2`$，$`\gamma=\cos^2t`$）。

<a id="r-11-15"></a>

**命题 11.15（$`M_2`$ 上的实验 1 与实验 3）。** 设 $`M_2(\mathbb C)`$ 中的 $`Y\ge0`$ 的特征值为
$`y_1\ge y_2\ge0`$，本征向量为 $`f_1,f_2`$，且 $`y_1>0`$。

(a) 设 $`X\ge0`$ 的特征值为 $`x_1\ge x_2\ge0`$，$`x_1>0`$，本征向量为 $`e_1,e_2`$。则对每个酉矩阵 $`U`$，
```math
\lvert\mathrm{Tr}(XUY)\rvert-\mathrm{Tr}(XY)\ \le\ \frac{(x_1-x_2)^2(y_1-y_2)^2}{4(x_1+x_2)(y_1+y_2)} .
```
若 $`x_1>x_2`$ 且 $`y_1>y_2`$，则等号成立当且仅当 $`|\mathrm{Tr}(XUY)|=\|YX\|_1`$ 且
$`\gamma=\gamma_1:=\frac14+\frac{x_1y_2+x_2y_1}{2(x_1+x_2)(y_1+y_2)}`$；否则该界为 $`0`$，且它被取到当且仅当
$`|\mathrm{Tr}(XUY)|=\|YX\|_1`$。

(b) 设 $`X=X^*`$，$`|X|`$ 的特征值为 $`x_1\ge x_2\ge0`$，$`x_1>0`$。则
```math
\mathrm{Tr}(Y|X|Y)-\mathrm{Tr}|YXY|\ \le\ \frac{2x_1x_2}{x_1+x_2}\,(y_1-y_2)^2 .
```
若 $`X`$ 有正负两种符号的特征值且 $`y_1>y_2`$，令 $`e_1`$、$`e_2`$ 为 $`X`$ 对应于模为 $`x_1`$ 与 $`x_2`$ 的特征值的
本征向量；则等号成立当且仅当 $`\gamma=\gamma_3:=\frac{x_1y_2+x_2y_1}{(x_1+x_2)(y_1+y_2)}`$。

对上述任意特征值，两个界都能取到，且可由实的 $`X`$、$`Y`$、$`U`$ 取到（在 (b) 中取 $`X=\mathrm{diag}(x_1,-x_2)`$）；
特别地，它们在 $`M_2(\mathbb R)`$ 中相同。

*证明。* 记 $`z`$ 为左边，并令 $`S=x_1y_2+x_2y_1`$，$`p=(x_1-x_2)(y_1-y_2)\ge0`$，
$`q=(x_1+x_2)(y_1+y_2)>0`$，于是 $`q-2S=p`$。

(a) 由定理 1.9，$`|\mathrm{Tr}(XUY)|=|\mathrm{Tr}(YX\,U)|\le\|YX\|_1`$，当 $`U=W^*`$、$`YX=W|YX|`$ 时取等号。
与推论 11.5 的证明一样，对 $`P=\sum_ia_i|e_i\rangle\langle e_i|`$、$`Q=\sum_jb_j|f_j\rangle\langle f_j|`$ 有
$`\mathrm{Tr}(PQ)=\sum_{i,j}a_ib_j|\langle e_i,f_j\rangle|^2`$；因此
$`\mathrm{Tr}(XY)=(x_1y_2+x_2y_1)+\gamma(x_1-x_2)(y_1-y_2)=S+\gamma p`$，且
$`\mathrm{Tr}(X^2Y^2)=(x_1^2y_2^2+x_2^2y_1^2)+\gamma(x_1^2-x_2^2)(y_1^2-y_2^2)`$。$`YX`$ 的奇异值的平方 $`s_1^2,s_2^2`$
是 $`(YX)^*YX=XY^2X`$ 的特征值，所以 $`s_1^2+s_2^2=\mathrm{Tr}(X^2Y^2)`$，且
$`s_1s_2=|\det(YX)|=x_1x_2y_1y_2`$。因此
```math
\|YX\|_1^2=(s_1+s_2)^2=\mathrm{Tr}(X^2Y^2)+2x_1x_2y_1y_2=S^2+\gamma pq .
```
对 $`r\ge0`$，$`\sqrt r\le r/q+q/4`$，等号当且仅当 $`\sqrt r=q/2`$（展开 $`(\sqrt r-q/2)^2\ge0`$）。取 $`r=S^2+\gamma pq`$，
```math
z\ \le\ \|YX\|_1-\mathrm{Tr}(XY)\ \le\ \frac{S^2+\gamma pq}q+\frac q4-S-\gamma p=\frac{(q/2-S)^2}q=\frac{p^2}{4q}.
```
等号成立当且仅当 $`|\mathrm{Tr}(XUY)|=\|YX\|_1`$ 且 $`S^2+\gamma pq=q^2/4`$，即 $`\gamma pq=(q/2-S)(q/2+S)=\frac p2(q/2+S)`$。
若 $`p>0`$，这就是 $`\gamma=(q/2+S)/(2q)=\gamma_1`$，由于 $`0\le S\le q/2`$，它落在 $`[\frac14,\frac12]`$ 中；若 $`p=0`$，
它对每个 $`\gamma`$ 都成立。

(b) 把 $`X`$ 换成 $`-X`$ 两边都不变。若 $`X\ge0`$，则 $`YXY\ge0`$，$`|YXY|=YXY`$，于是 $`z=0`$；
$`X\le0`$ 时同理。否则 $`x_2>0`$，必要时把 $`X`$ 换成 $`-X`$，可设 $`Xe_1=x_1e_1`$，
$`Xe_2=-x_2e_2`$。令 $`u=\langle e_1,Y^2e_1\rangle=\gamma y_1^2+(1-\gamma)y_2^2`$，$`v=\langle e_2,Y^2e_2\rangle=s-u`$，$`s=y_1^2+y_2^2`$。
则 $`\mathrm{Tr}(Y|X|Y)=\mathrm{Tr}(|X|Y^2)=x_1u+x_2v=:A`$。自伴矩阵 $`YXY`$ 的行列式为
$`\det X(\det Y)^2=-x_1x_2y_1^2y_2^2\le0`$，所以其特征值为 $`\mu_+\ge0\ge\mu_-`$，且
```math
\mathrm{Tr}|YXY|=\mu_+-\mu_-=\sqrt{(\mu_++\mu_-)^2-4\mu_+\mu_-}=\sqrt{L^2+K^2},\qquad L=\mathrm{Tr}(XY^2)=x_1u-x_2v,\quad K=2\sqrt{x_1x_2}\,y_1y_2 .
```
令 $`\alpha=(x_1-x_2)/(x_1+x_2)`$，$`\beta=2\sqrt{x_1x_2}/(x_1+x_2)>0`$，于是 $`\alpha^2+\beta^2=1`$。由于
$`A=(x_1-x_2)u+x_2s`$，$`L=(x_1+x_2)u-x_2s`$，有 $`A-\alpha L=x_2s(1+\alpha)=2x_1x_2s/(x_1+x_2)`$ 与
$`\beta K=4x_1x_2y_1y_2/(x_1+x_2)`$，从而 $`A-\alpha L-\beta K=\frac{2x_1x_2}{x_1+x_2}(y_1-y_2)^2=:m`$。
$`\mathbb R^2`$ 中的 Cauchy–Schwarz 不等式给出 $`\sqrt{L^2+K^2}\ge\alpha L+\beta K`$，所以
```math
z=A-\sqrt{L^2+K^2}\ \le\ A-\alpha L-\beta K=m,
```
等号成立当且仅当 $`(L,K)`$ 是 $`(\alpha,\beta)`$ 的非负倍数，即（因为 $`\beta>0`$，$`K\ge0`$）当且仅当
$`\beta L=\alpha K`$，即 $`(x_1+x_2)u-x_2s=(x_1-x_2)y_1y_2`$。这确定了 $`u`$，且
$`u-y_2^2=(y_1-y_2)(x_1y_2+x_2y_1)/(x_1+x_2)`$；当 $`y_1>y_2`$ 时，由 $`u=\gamma y_1^2+(1-\gamma)y_2^2`$ 得 $`\gamma=\gamma_3\in[0,\frac12]`$。

*实的极大点。* 取 $`e_i`$ 为标准基，$`f_1,f_2`$ 为满足 $`\cos^2t=\gamma_1`$（相应地 $`\gamma_3`$）的旋转实基，
在 (a) 中由实奇异值分解 $`YX=W\Sigma V^{\mathrm T}`$（对实对称矩阵 $`(YX)^{\mathrm T}YX`$ 作对角化）取
$`U=VW^{\mathrm T}`$，于是 $`\mathrm{Tr}(UYX)=\mathrm{Tr}\,\Sigma=\|YX\|_1`$。$`\square`$

<a id="r-11-16"></a>

**推论 11.16（块大小为 $`2`$ 时实验 1 与实验 3 的包络）。** 设所有 $`k_c=2`$（于是 $`D=2C`$），
$`\delta\in[0,1)`$，$`\delta'\in[0,1]`$，$`\ell=\ell(\delta)`$，$`\ell'=\ell(\delta')`$（约定 $`\ell(1)=0`$），并让 $`Y`$
取遍满足 $`\lambda_{\max}(Y)=1`$、$`\Delta(Y)=\delta'`$（即 $`\lambda_{\min}(Y)=\ell'`$）的正算子。

(a) 实验 1，在 $`X\in\mathcal P_\delta`$、这样的 $`Y`$ 与酉算子 $`U`$ 上：
$`\ \sup z=C\,\frac{\delta^2\delta'^2}{(1+\delta)(1+\delta')}`$，$`\inf z=-(2C-1+\ell\ell')`$。

(b) 实验 3，在满足 $`|X|\in\mathcal P_\delta`$ 的 $`X=X^*`$ 与这样的 $`Y`$ 上：
$`\ \sup z=(C-\delta)\Big(\frac{2\delta'}{1+\delta'}\Big)^2`$，$`\inf z=0`$。

四个值都能取到。因此，当另一个对比度在 $`[0,d]`$（$`d\in(0,1)`$）中自由变化时：实验 1 的上包络沿
$`\Delta(X)=\delta`$ 为 $`C\delta^2d^2/((1+\delta)(1+d))`$，（由 (a) 的对称性）沿 $`\Delta(Y)=\delta'`$ 为
$`C\delta'^2d^2/((1+\delta')(1+d))`$；实验 3 的上包络沿 $`\Delta(|X|)=\delta`$ 为 $`(C-\delta)\big(2d/(1+d)\big)^2`$，沿
$`\Delta(Y)=\delta'`$ 为 $`C\big(2\delta'/(1+\delta')\big)^2`$，与 $`d`$ 无关。当 $`d\to1`$ 时，前两者递增地趋于
$`C\delta^2/(2(1+\delta))`$ 与 $`C-\delta`$，即 $`\Delta(Y)=1`$（$`Y`$ 奇异）时的值。

*证明。* 设 $`x_{c,1}\ge x_{c,2}`$ 为 $`X_c`$（在 (b) 中为 $`|X|_c`$）的特征值，$`y_{c,1}\ge y_{c,2}`$ 为
$`Y_c`$ 的特征值；它们分别落在 $`[\ell,1]`$ 与 $`[\ell',1]`$ 中，而 $`\ell`$ 是它所在块中较小的那个。注意
$`1-\ell=2\delta/(1+\delta)`$，$`1+\ell=2/(1+\delta)`$，故 $`(1-\ell)^2/(1+\ell)=2\delta^2/(1+\delta)`$ 且
$`2\ell/(1+\ell)=1-\delta`$；$`\ell'`$ 同理。

(a) *上确界。* 由三角不等式以及在每个块上应用定理 1.9，
$`|\mathrm{Tr}(XUY)|\le\sum_c|\mathrm{Tr}(X_cU_cY_c)|\le\sum_c\|Y_cX_c\|_1`$，于是由命题 11.15(a)（$`Y_c=0`$ 的块贡献
$`-\mathrm{Tr}(X_cY_c)=0`$）
```math
z\le\sum_c\frac{x_{c,1}y_{c,1}}4\,g\Big(\frac{x_{c,2}}{x_{c,1}}\Big)\,g\Big(\frac{y_{c,2}}{y_{c,1}}\Big),\qquad g(r)=\frac{(1-r)^2}{1+r}.
```
$`g`$ 在 $`[0,1]`$ 上递减，$`x_{c,2}/x_{c,1}\ge x_{c,2}\ge\ell`$，$`y_{c,2}/y_{c,1}\ge\ell'`$，且 $`x_{c,1},y_{c,1}\le1`$，所以
每一项至多为 $`g(\ell)g(\ell')/4=\delta^2\delta'^2/((1+\delta)(1+\delta'))`$。取等：每个块中
$`X_c=\mathrm{diag}(1,\ell)`$，$`Y_c`$ 的特征值为 $`1,\ell'`$ 且 $`\gamma=\gamma_1`$，$`U_c=W_c^*`$，$`Y_cX_c=W_c|Y_cX_c|`$；
这时每个块中 $`\mathrm{Tr}(X_cU_cY_c)=\|Y_cX_c\|_1\ge0`$，模长相加。
*下确界。* $`z\ge-\mathrm{Tr}(XY)`$，当 $`U=\bigoplus_cU_c`$ 满足 $`\mathrm{Tr}(Y_cX_cU_c)=0`$（引理 11.13）时取等号。由
推论 11.5(a)，$`\mathrm{Tr}(XY)\le\sum_c(x_{c,1}y_{c,1}+x_{c,2}y_{c,2})`$。若 $`\ell`$ 与 $`\ell'`$ 在同一个块中，
该块贡献至多 $`1+\ell\ell'`$，其余每块至多 $`2`$；否则这两个块分别贡献至多 $`1+\ell`$ 与 $`1+\ell'`$，而
$`2C-2+\ell+\ell'=2C-1+\ell\ell'-(1-\ell)(1-\ell')`$。所以
$`\mathrm{Tr}(XY)\le2C-1+\ell\ell'`$，当 $`X=\mathrm{diag}(1,\ell)\oplus1\oplus\dots\oplus1`$、
$`Y=\mathrm{diag}(1,\ell')\oplus1\oplus\dots\oplus1`$ 时取等号。

(b) 迹、迹范数与 $`|\cdot|`$ 都按块作用（习题 3），所以 $`z=\sum_cz_c`$，由命题 11.15(b)，
$`z_c\le H_c(y_{c,1}-y_{c,2})^2\le H_c(1-\ell')^2`$，其中 $`H_c=2x_{c,1}x_{c,2}/(x_{c,1}+x_{c,2})\le x_{c,1}\le1`$。在 $`|X|`$
具有特征值 $`\ell`$ 的块 $`c_-`$ 中，$`x_{c_-,2}=\ell`$，且 $`H_{c_-}=2\ell/(1+\ell/x_{c_-,1})\le2\ell/(1+\ell)=1-\delta`$。
因此 $`z\le(C-1+1-\delta)(1-\ell')^2`$。取等：$`X_1=\mathrm{diag}(1,-\ell)`$，$`c\ge2`$ 时 $`X_c=\mathrm{diag}(1,-1)`$，
每个 $`Y_c`$ 的特征值为 $`1,\ell'`$ 且取其所在块的 $`\gamma=\gamma_3`$（$`c\ge2`$ 时 $`\gamma_3=\frac12`$）。下确界即
命题 11.9，$`X`$、$`Y`$ 可交换时取等号。

*自由对比度。* (a) 中的上确界随 $`\delta'`$ 递增（若 $`\delta>0`$ 则严格递增），(b) 中的上确界随 $`\delta'`$ 严格递增、
随 $`\delta`$ 递减。因此在 $`[0,d]`$ 中自由的对比度上，上确界在 $`d`$ 处取到（实验 3 沿 $`\Delta(Y)`$ 时：在
$`\Delta(|X|)=0`$ 处），而 $`d\to1`$ 时的极限只在 $`\Delta(Y)=1`$ 处达到（对实验 1 要求 $`\delta>0`$）。$`\square`$

对 $`C=1`$、$`d=0.999`$，这正是配套论文最早在数值上观察到的闭式（在实代数中）。对每个 $`C`$，上包络都是
$`M_2`$ 上 $`C`$ 个解的直和所给出的值（实验 3 中 $`C-\delta=(1-\delta)+(C-1)\cdot1`$）：在 $`M_2\oplus M_2`$ 上，
上面的直和不等式取等号。第 11.7 节把这些上确界推广到任意大小的块（定理 11.20）。

*计算：* 命题 11.15 的极大点、$`C=1,2`$ 时推论 11.16 的包络及其 $`d\to1`$ 的极限，以及与之比较的多起点搜索——[教程，11.5 节](../../tutorial/zh/11_optimisation.md#s-11-5)。

<a id="s-11-6"></a>

## 11.6 可与之比较的工具

*库中：* 与其他软件的比较，其中包括 SDP 求解器，其凸松弛从上方界定上确界——[教程，11.6 节](../../tutorial/zh/11_optimisation.md#s-11-6)。

<a id="s-11-7"></a>

## 11.7 任意大小的块的包络

推论 11.16 中的上确界可推广到任意大小的块。令
```math
N=\sum_c\Big\lfloor\frac{k_c}2\Big\rfloor,
```
即各扇区中能容纳的两两正交的二维子空间的最大个数。上包络由命题 11.15 的极大点在 $`N`$ 个大小为 $`2`$ 的块上与
$`1\times1`$ 块的直和取到。证明先把谱化为两个值，再把算子化为大小为 $`1`$ 和 $`2`$ 的块（引理 11.17），在那里可用
命题 11.15——唯一的例外是实验 3 中 $`|X|`$ 属于 $`\lambda_{\min}`$ 的特征向量分散在多个块上的情形；命题 11.19 用线性化
处理这一情形。

*顶点论证。* 立方体的一个面 $`F=\{x\in[s,t]^n:x_i=t\ (i\in I),\ x_j=s\ (j\in J)\}`$ 上的凸函数 $`f`$ 在 $`F`$ 的某个顶点
（所有坐标都属于 $`\{s,t\}`$ 的点）取到它在 $`F`$ 上的最大值：$`F`$ 的每个点都是这些顶点的凸组合（把一个自由坐标写成
$`s`$ 与 $`t`$ 的凸组合，再对自由坐标的个数归纳），而 $`f`$ 在该点的值不超过它在各顶点处的最大值。

<a id="r-11-17"></a>

**引理 11.17（两个投影）。** 对投影 $`P,Q\in M_n(\mathbb C)`$，$`\mathbb C^n`$ 是若干维数为 $`1`$ 或 $`2`$ 的子空间的正交和，
每个子空间在 $`P`$ 与 $`Q`$ 下不变，且每个二维子空间上有标准正交基 $`f,f'`$，使得在其上 $`P=|f\rangle\langle f|`$，
$`Q=|q\rangle\langle q|`$，$`q=\cos\theta\,f+e^{i\psi}\sin\theta\,f'`$，$`\theta\in[0,\pi/2]`$。对实的 $`P,Q`$，一切都可取为实的。

*证明。* $`PQP`$ 在 $`\mathrm{ran}\,P`$ 上自伴；取其标准正交特征向量基 $`v_1,\dots,v_r`$，$`PQPv_i=\mu_iv_i`$（定理 1.3），
令 $`W_i=\mathrm{span}\{v_i,Qv_i\}`$。则 $`PQv_i=PQPv_i=\mu_iv_i`$，故 $`W_i`$ 在 $`P`$ 下不变，并且（由于 $`Q^2=Q`$）在
$`Q`$ 下不变。对 $`i\ne j`$，$`\langle v_i,Qv_j\rangle=\langle Qv_i,Qv_j\rangle=\langle v_i,PQPv_j\rangle=0`$，因此各 $`W_i`$ 两两正交。
它们的和包含 $`\mathrm{ran}\,P`$；其正交补 $`W`$ 在 $`Q`$ 下不变（两个投影都是自伴的）且含于 $`\ker P`$，$`Q`$ 在 $`W`$ 上的
标准正交特征向量基给出一维直和项。维数为 $`1`$ 的 $`W_i`$ 由 $`v_i`$ 张成。若 $`\dim W_i=2`$，令 $`f=v_i`$，$`f'`$ 为 $`Qv_i`$
与 $`f`$ 正交的分量的归一化；则在 $`W_i`$ 上 $`P=|f\rangle\langle f|`$，而 $`Q`$ 把 $`W_i`$ 映到 $`\mathrm{span}\{Qv_i\}`$ 上，
故在那里对某单位向量 $`q=\alpha f+\beta f'`$ 有 $`Q=|q\rangle\langle q|`$；把 $`q`$ 乘以一个相位可使 $`\alpha\ge0`$，于是
$`\alpha=\cos\theta`$，$`\beta=e^{i\psi}\sin\theta`$，$`\theta\in[0,\pi/2]`$。在实的情形，特征向量基、$`f'`$ 与 $`q`$ 都是实的，
且 $`e^{i\psi}=\pm1`$。$`\square`$

<a id="r-11-18"></a>

**引理 11.18（实验 3：单调性与凸性）。** 设 $`X=X^*`$ 可逆，$`P`$ 为到其正部的投影（$`X_+=PX`$），$`S=2P-1`$，
$`Y\ge0`$，$`z(X,Y)=\mathrm{Tr}(Y|X|Y)-\mathrm{Tr}|YXY|`$。

(a) 若 $`|X|\le|X'|`$，其中 $`|X'|`$ 与 $`P`$ 可交换且 $`X'=S|X'|`$，则 $`z(X,Y)\le z(X',Y)`$。

(b) 固定 $`X`$ 时，$`z`$ 是 $`K=Y^2`$ 的凸函数：对 $`K_1,K_2\ge0`$ 与 $`t\in[0,1]`$，
$`z\big(X,(tK_1+(1-t)K_2)^{1/2}\big)\le t\,z(X,K_1^{1/2})+(1-t)\,z(X,K_2^{1/2})`$。

*证明。* 令 $`A=YX_+Y`$，$`B=YX_-Y`$，则 $`z=\mathrm{Tr}\,A+\mathrm{Tr}\,B-\|A-B\|_1`$（命题 11.9 的证明）。

(a) $`X_+=P|X|=|X|^{1/2}P|X|^{1/2}`$，$`X_-=(1-P)|X|`$，且 $`P|X'|-P|X|=P(|X'|-|X|)P\ge0`$，对 $`1-P`$ 亦然；故对 $`X'`$ 的
算子 $`A',B'`$ 有 $`A\le A'`$，$`B\le B'`$（推论 1.7）。对 $`R\ge0`$，$`\|A+R-B\|_1\le\|A-B\|_1+\mathrm{Tr}\,R`$（三角不等式，
$`\|R\|_1=\mathrm{Tr}\,R`$），所以当 $`A`$ 或 $`B`$ 增大时 $`\mathrm{Tr}\,A+\mathrm{Tr}\,B-\|A-B\|_1`$ 不减。

(b) 对自伴的 $`H`$ 与满足 $`-1\le T\le1`$ 的自伴 $`T`$，有 $`|\mathrm{Tr}(TH)|\le\|H\|_1`$：算子 $`T\pm i(1-T^2)^{1/2}`$ 是
酉的（定理 1.5），其平均为 $`T`$，对每一个应用定理 1.9 即可。当 $`T=2Q-1`$、$`Q`$ 为到 $`H`$ 正部的投影时等号成立
（引理 11.8 的证明）。取 $`T=1-2E`$，
```math
z=2\min_{0\le E\le1}\mathrm{Tr}\big(EA+(1-E)B\big).
```
先设 $`Y`$ 可逆。则 $`N=|X|^{1/2}Y`$ 可逆，$`A=N^*PN`$，$`B=N^*(1-P)N`$，且 $`E\mapsto F=NEN^*`$ 把 $`\{0\le E\le1\}`$ 映满
$`\{0\le F\le G\}`$，$`G=NN^*=|X|^{1/2}Y^2|X|^{1/2}`$（推论 1.7，两个方向）。由于 $`\mathrm{Tr}(EA)=\mathrm{Tr}(PF)`$ 且
$`\mathrm{Tr}((1-E)B)=\mathrm{Tr}((1-P)(G-F))`$，$`z=2\min\{\mathrm{Tr}(PF+(1-P)(G-F)):0\le F\le G\}`$。若 $`0\le F_i\le G_i`$，则
$`tF_1+(1-t)F_2`$ 对 $`tG_1+(1-t)G_2`$ 可行，且目标函数关于 $`(F,G)`$ 是线性的；所以 $`z`$ 关于 $`G`$ 是凸的，而 $`G`$ 关于
$`Y^2`$ 是仿射的。若 $`K_1`$ 或 $`K_2`$ 奇异，把上述结论用于 $`K_i+\varepsilon`$ 并令 $`\varepsilon\to0`$（$`z`$ 关于 $`Y`$ 连续，
$`K\mapsto K^{1/2}`$ 连续）。$`\square`$

<a id="r-11-19"></a>

**命题 11.19（$`m`$ 个块上的实验 3）。** 设 $`\mathbb C^m\oplus\mathbb C^m`$ 上 $`X=X_+\oplus(-1)`$，其中 $`0<X_+\le1`$ 的特征值为
$`x_1,\dots,x_m`$；设 $`C=\mathrm{diag}(\cos\theta_j)`$，$`\Sigma=\mathrm{diag}(\sin\theta_j)`$，$`Q`$ 为到等距映射
$`v=\binom C\Sigma:\mathbb C^m\to\mathbb C^{2m}`$ 的值域的投影，$`b\in[0,1)`$，$`Y=b+(1-b)Q`$。则
```math
\mathrm{Tr}(Y|X|Y)-\mathrm{Tr}|YXY|\ \le\ (1-b)^2\sum_{i=1}^m\frac{2x_i}{1+x_i},
```
当 $`X_+=\mathrm{diag}(x_1,\dots,x_m)`$ 且 $`\cos^2\theta_i=(1+bx_i)/((1+x_i)(1+b))`$ 时等号成立。

*证明。* 两个线性代数事实：对方阵 $`R,T`$，$`\det(\lambda-RT)=\det(\lambda-TR)`$（$`R`$ 可逆时 $`RT=R(TR)R^{-1}`$；一般情形
由关于 $`R`$ 的连续性得到）；对右下块 $`E`$ 可逆的分块矩阵，$`\det\begin{pmatrix}R&S\\T&E\end{pmatrix}=\det E\,\det(R-SE^{-1}T)`$
（左乘 $`\begin{pmatrix}1&-SE^{-1}\\0&1\end{pmatrix}`$）。

令 $`\gamma=1-b^2`$，$`K=1+X_+`$，$`A=C^2`$，于是 $`\Sigma^2=1-A`$；$`Q=vv^*`$ 的各块为 $`C^2,\ C\Sigma,\ C\Sigma,\ \Sigma^2`$，
且 $`Y^2=b^2+\gamma Q`$（因为 $`Q^2=Q`$）。对 $`\lambda>0`$，$`\lambda-XY^2`$ 的右下块是可逆对角矩阵
$`E=\lambda+b^2+\gamma\Sigma^2=\lambda+1-\gamma A`$；由于 $`C,\Sigma,E`$ 都是对角的，且
$`(b^2+\gamma A)E-\gamma^2A\Sigma^2=\lambda(b^2+\gamma A)+b^2`$，
```math
\det(\lambda-XY^2)=\det E\,\det\Big(\lambda-X_+\big(\lambda(b^2+\gamma A)+b^2\big)E^{-1}\Big)
=\det\big((\lambda+1)(\lambda-b^2X_+)-\lambda\gamma KA\big)=p(\lambda),
```
其中 $`p(\lambda)=\det\big((\lambda+1)(\lambda-b^2X_+)-\lambda\gamma K^{1/2}AK^{1/2}\big)`$（用 $`K^{1/2}`$ 共轭矩阵
$`K^{-1}(\lambda+1)(\lambda-b^2X_+)-\lambda\gamma A`$）。对自伴矩阵
```math
\mathcal L=\begin{pmatrix}\gamma K^{1/2}AK^{1/2}+b^2X_+-1&bX_+^{1/2}\\ bX_+^{1/2}&0\end{pmatrix}
```
用右下块为 $`\lambda`$ 的同一公式得
$`\det(\lambda-\mathcal L)=\lambda^m\det(\lambda+1-\gamma K^{1/2}AK^{1/2}-b^2X_+-b^2X_+/\lambda)=p(\lambda)`$。由于
$`\det(\lambda-YXY)=\det(\lambda-XY^2)`$，$`YXY`$ 与 $`\mathcal L`$ 的特征多项式在 $`\lambda>0`$ 上一致，从而处处一致，这两个
自伴矩阵有相同的特征值：$`\mathrm{Tr}|YXY|=\|\mathcal L\|_1`$。令 $`\Delta=(1-X_+)(1+X_+)^{-1}`$，$`\Gamma=2X_+^{1/2}(1+X_+)^{-1}`$
（$`X_+`$ 的函数，故可交换，且 $`\Delta^2+\Gamma^2=1`$），以及
```math
T=\begin{pmatrix}-\Delta&\Gamma\\ \Gamma&\Delta\end{pmatrix},\qquad T^2=1 .
```
$`T`$ 自伴，特征值为 $`\pm1`$，故 $`\|\mathcal L\|_1\ge\mathrm{Tr}(T\mathcal L)`$（引理 11.18(b) 的证明）。而
$`\mathrm{Tr}(Y|X|Y)=\mathrm{Tr}(|X|Y^2)=b^2\mathrm{Tr}\,X_++\gamma\mathrm{Tr}(X_+A)+mb^2+\gamma\mathrm{Tr}(1-A)`$，并且由于
$`K^{1/2}\Delta K^{1/2}=1-X_+`$，
$`\mathrm{Tr}(T\mathcal L)=-\gamma\mathrm{Tr}((1-X_+)A)-\mathrm{Tr}(\Delta(b^2X_+-1))+2b\,\mathrm{Tr}(\Gamma X_+^{1/2})`$。含 $`A`$ 的项
相互抵消，在 $`X_+`$ 的特征向量基下
```math
\mathrm{Tr}(Y|X|Y)-\mathrm{Tr}|YXY|\le\sum_{i=1}^m\Big(b^2x_i+1+\frac{1-x_i}{1+x_i}(b^2x_i-1)-\frac{4bx_i}{1+x_i}\Big)
=(1-b)^2\sum_{i=1}^m\frac{2x_i}{1+x_i}.
```
当 $`X_+`$ 为对角时，$`X`$ 与 $`Y`$ 是 $`\mathrm{span}\{e_i,e_{m+i}\}`$ 上 $`2\times2`$ 块的直和，在其上 $`X=\mathrm{diag}(x_i,-1)`$，
命题 11.15(b) 在 $`\gamma_3=\sin^2\theta_i`$（即所述角度）时给出值 $`\frac{2x_i}{1+x_i}(1-b)^2`$。$`\square`$

由于 $`\mathcal L`$ 关于 $`A=PQP`$ 是仿射的（$`P`$ 为到第一个直和项的投影），$`z`$ 是 $`A`$ 的凹函数，而 $`T`$ 是极大点处
$`\mathcal L`$ 的符号。然而，固定角度时，关于 $`X_+`$ 的最小特征值的特征向量的最大值一般并不在单个块内取到（一个数值
观察）；只有联合最大值如此。

<a id="r-11-20"></a>

**定理 11.20（实验 1 与实验 3 的上包络）。** 设 $`N\ge1`$，$`\delta\in[0,1)`$，$`\delta'\in[0,1]`$，$`\ell=\ell(\delta)`$，
$`\ell'=\ell(\delta')`$（约定 $`\ell(1)=0`$），并令 $`Y`$ 取遍满足 $`\lambda_{\max}(Y)=1`$ 与 $`\Delta(Y)=\delta'`$ 的正算子。

(a) 实验 1，在 $`X\in\mathcal P_\delta`$、这样的 $`Y`$ 与酉算子 $`U\in\mathcal M`$ 上：
$`\ \sup z=N\,\frac{\delta^2\delta'^2}{(1+\delta)(1+\delta')}`$。

(b) 实验 3，在满足 $`|X|\in\mathcal P_\delta`$ 的 $`X=X^*`$ 与这样的 $`Y`$ 上：若所有 $`k_c`$ 为偶数，
$`\ \sup z=(N-\delta)\Big(\frac{2\delta'}{1+\delta'}\Big)^2`$；若某个 $`k_c`$ 为奇数，
$`\ \sup z=N\Big(\frac{2\delta'}{1+\delta'}\Big)^2`$。

两个上确界都可由实算子取到，在 $`\bigoplus_cM_{k_c}(\mathbb R)`$ 中相同。

*证明。* $`\mathcal M`$ 中的上界对实算子更是成立，而下面的极大点都是实的。

(a) 对每个块用定理 1.9，$`|\mathrm{Tr}(XUY)|\le\|YX\|_1`$，当 $`U=W^*`$、$`YX=W|YX|`$（$`\mathcal M`$ 中的酉元）时取等号；所以
$`\sup_Uz=F(X,Y):=\|YX\|_1-\mathrm{Tr}(XY)`$。对 $`X`$ 的固定的（与扇区相适应的）标准正交特征向量基 $`(e_s)`$，$`F`$ 是特征值
$`x\in[\ell,1]^D`$ 的凸函数（$`x`$ 的线性函数的范数减去一个线性函数），而 $`X\in\mathcal P_\delta`$ 意味着 $`x`$ 位于某个面
$`\{x_s=1,\ x_{s'}=\ell\}`$ 上；由顶点论证，对某投影 $`P\in\mathcal M`$ 有 $`F(X,Y)\le F(X',Y)`$，$`X'=\ell+(1-\ell)P`$。同理
$`F(X',Y)\le F(X',Y')`$，$`Y'=\ell'+(1-\ell')Q`$。在每个扇区中对 $`P_c,Q_c`$ 应用引理 11.17：$`X'`$ 与 $`Y'`$ 是维数为 $`1`$ 和
$`2`$ 的子空间上的直和，$`\|Y'X'\|_1`$ 与 $`\mathrm{Tr}(X'Y')`$ 亦然。一维直和项贡献 $`xy-xy=0`$。在二维直和项上，$`X'`$ 的
特征值为 $`1,\ell`$，$`Y'`$ 的特征值为 $`1,\ell'`$，故由命题 11.15(a) 其贡献至多为
$`(1-\ell)^2(1-\ell')^2/(4(1+\ell)(1+\ell'))=\delta^2\delta'^2/((1+\delta)(1+\delta'))`$（见推论 11.16 的证明）。大小为 $`k_c`$ 的
扇区至多含 $`\lfloor k_c/2\rfloor`$ 个这样的子空间。等号：$`N`$ 个块上 $`X=\mathrm{diag}(1,\ell)`$，$`Y`$ 的特征值为 $`1,\ell'`$
且 $`\gamma=\gamma_1`$，其余一维块上 $`X=Y=1`$，以及酉元 $`W^*`$。

(b) $`X`$ 可逆（$`\ell>0`$），且 $`z(-X,Y)=z(X,Y)`$。令 $`P`$ 为到 $`X`$ 正部的投影，$`S=2P-1`$。$`|X|`$ 与 $`P`$ 可交换（二者
都是 $`X`$ 的函数），因此 $`|X|`$ 有属于特征值 $`\ell`$ 的单位特征向量 $`e`$，它位于单个扇区中，且属于 $`\mathrm{ran}\,P`$ 或
$`\ker P`$；必要时把 $`X`$ 换成 $`-X`$，可设 $`e\in\mathrm{ran}\,P`$。于是 $`|X|\le|X'|:=1-(1-\ell)|e\rangle\langle e|`$，
$`X'=S|X'|=S-(1-\ell)|e\rangle\langle e|`$ 满足 $`|X'|\in\mathcal P_\delta`$（$`D\ge2`$），由引理 11.18(a)，$`z(X,Y)\le z(X',Y)`$。
由引理 11.18(b) 及对 $`Y^2`$ 在 $`[\ell'^2,1]`$ 中的特征值的顶点论证，对某投影 $`Q\in\mathcal M`$ 有 $`z(X',Y)\le z(X',Y')`$，
$`Y'=\ell'+(1-\ell')Q`$。在每个扇区中对 $`P_c,Q_c`$ 应用引理 11.17。
在不含 $`e`$ 的扇区中，$`X'=S`$ 是各子空间上的直和；一维直和项贡献 $`0`$（在那里 $`X'`$ 与 $`Y'`$ 可交换，命题 11.9），二维
直和项由命题 11.15(b)（取 $`x_1=x_2=1`$）至多贡献 $`(1-\ell')^2`$：合计至多 $`\lfloor k_c/2\rfloor(1-\ell')^2`$。
在 $`e`$ 所在的扇区中，由引理 11.18(a) 与 $`|X'|\le1`$ 得同样的界，因此若 $`k_c`$ 为奇数，贡献至多 $`\lfloor k_c/2\rfloor(1-\ell')^2`$；
若该扇区所含二维子空间少于 $`k_c/2`$ 个，则至多 $`(k_c/2-1)(1-\ell')^2`$。否则它的 $`m=k_c/2`$ 个二维子空间填满该扇区，
其基为 $`f_j,f_j'`$，$`f_j\in\mathrm{ran}\,P`$，$`f_j'\in\ker P`$。把 $`f_j'`$ 乘以 $`e^{-i\psi_j}`$ 的酉元与 $`X'`$ 可交换
（$`X'`$ 在 $`\ker P`$ 上等于 $`-1`$），并消去引理 11.17 中的相位；在基 $`(f_1,\dots,f_m,f_1',\dots,f_m')`$ 下，该扇区正是
命题 11.19 的情形，其中 $`X_+=1-(1-\ell)|e\rangle\langle e|`$，特征值为 $`\ell,1,\dots,1`$；由于 $`2\ell/(1+\ell)=1-\delta`$，其贡献
至多为 $`(1-\ell')^2(1-\delta+m-1)=(k_c/2-\delta)(1-\ell')^2`$。
对各扇区求和：若所有 $`k_c`$ 为偶数，$`z\le(N-\delta)(1-\ell')^2`$（注意 $`N-1\le N-\delta`$），否则 $`z\le N(1-\ell')^2`$；
这里 $`1-\ell'=2\delta'/(1+\delta')`$。
等号：若所有 $`k_c`$ 为偶数，在一个块上 $`X=\mathrm{diag}(\ell,-1)`$，在其余 $`N-1`$ 个块上 $`X=\mathrm{diag}(1,-1)`$，且每个块上
$`Y`$ 的特征值为 $`1,\ell'`$、$`\gamma=\gamma_3`$（命题 11.15(b)）；若某个 $`k_c`$ 为奇数，在 $`N`$ 个块上 $`X=\mathrm{diag}(1,-1)`$，
特征值 $`\ell`$ 位于某个奇数扇区的一维块上，其余一维块上为 $`1`$，$`Y`$ 在各块上同前，在一维块上为 $`1`$。$`\square`$

对 $`k_c=2`$，这就是推论 11.16 中关于上确界的部分（$`N=C`$）。对 $`C`$ 个大小为 $`k`$ 的扇区、另一个对比度在 $`[0,d]`$
（$`d\in(0,1)`$）中自由时，上确界关于 $`\delta'`$ 递增，上包络为：实验 1 中沿任一对比度为
$`C\lfloor k/2\rfloor\delta^2d^2/((1+\delta)(1+d))`$；实验 3 中沿 $`\Delta(|X|)=\delta`$，$`k`$ 为偶数时为
$`(Ck/2-\delta)\big(2d/(1+d)\big)^2`$，$`k`$ 为奇数时为 $`C\lfloor k/2\rfloor\big(2d/(1+d)\big)^2`$；沿 $`\Delta(Y)=\delta'`$（在
$`|X|=1`$ 处取到）为 $`C\lfloor k/2\rfloor\big(2\delta'/(1+\delta')\big)^2`$。$`k`$ 为偶数时，这些正是 11.5 节中直和的值：
对两个泛函，直和不等式都取等号。$`k`$ 为奇数时，带一维直和项的直和在实验 1 中是最优的，在实验 3 中却不是：$`|X|`$ 的
特征值 $`\ell`$ 放在一维块上，在那里没有代价，包络与 $`\delta`$ 无关。

*计算：* 命题 11.19（$`YXY`$ 与 $`\mathcal L`$ 的谱、该界及其等号情形）、定理 11.20 对若干块大小的极大点，以及从结构化
起点出发达到包络的多起点搜索——[教程，11.7 节](../../tutorial/zh/11_optimisation.md#s-11-7)。

<a id="x-11"></a>

## 习题

1. （**Cayley 变换。**）对 $`S^*=-S`$，证明 $`C(S)=(1-S)(1+S)^{-1}`$ 是酉的，并且 $`S\mapsto C(S)`$ 是到
   没有特征值 $`-1`$ 的酉矩阵集合上的双射。因此，与指数映射（引理 11.1）不同，Cayley 坐标卡
   $`U_0C(S)`$ 取不到 $`-U_0`$。
2. （**实 Haar 测度。**）设 $`O`$ 在 $`O(k)`$ 上服从 Haar 分布（实 Ginibre 矩阵的 QR 分解），$`A`$ 为实矩阵，证明
   $`\mathbb E\mathrm{Tr}(AO)^2=\|A\|_2^2/k`$，并用实代数加以检验。
3. （**直和。**）证明实验 3 的 $`z`$ 对直和是可加的，并且对 $`\lambda_{\max}=1`$ 的 $`A,A'\ge0`$ 有
   $`\Delta(A\oplus A')=\max(\Delta(A),\Delta(A'))`$。由此推出 11.5 节中的直和
   不等式（其中 $`\Delta(Y)`$ 在固定区间 $`[0,b]`$ 中自由变化）。

<details>
<summary>解答</summary>

1. $`S`$ 是正规的，且 $`iS`$ 是自伴的，所以 $`S`$ 的特征值是纯虚数（推论 1.4(2)），从而 $`1\pm S`$ 可逆；
   它们可交换，且 $`(1+S)^*=1-S`$。因此
   $`C^*C=(1-S)^{-1}(1+S)(1-S)(1+S)^{-1}=1`$，而 $`1+C=2(1+S)^{-1}`$ 可逆，所以 $`-1\notin\mathrm{spec}\,C`$。
   反之，对满足 $`1+U`$ 可逆的酉矩阵 $`U`$，$`S=(1-U)(1+U)^{-1}`$ 是反自伴的：
   $`S^*=(1+U^*)^{-1}(1-U^*)=(U+1)^{-1}U\,U^*(U-1)=(U+1)^{-1}(U-1)=-S`$，因为 $`U\pm1`$ 可交换。对可交换的 $`C`$ 与
   $`S`$，$`C(1+S)=1-S`$ 等价于 $`S(1+C)=1-C`$（二者都是说 $`C+CS+S=1`$）。而 $`C(S)`$ 与 $`S`$ 可交换，
   $`S=(1-U)(1+U)^{-1}`$ 与 $`U`$ 可交换；所以 $`U=C(S)`$ 当且仅当 $`S=(1-U)(1+U)^{-1}`$：这两个映射互为逆映射。

2. 由定理 B.17（另见引理 C.2），实高斯矩阵经符号修正的 QR 因子 $`O`$——本库对实代数正是这样抽样的——
   在 $`O(k)`$ 上服从 Haar 分布，所以对 $`V\in O(k)`$，$`VO`$ 与 $`O`$ 同分布。
   左乘 $`\mathrm{diag}(\pm1)`$ 给出当 $`i\ne l`$ 时 $`\mathbb E[O_{ij}O_{lm}]=0`$，置换表明
   $`\mathbb E[O_{ij}O_{im}]`$ 与 $`i`$ 无关，而 $`\sum_iO_{ij}O_{im}=\delta_{jm}`$；所以
   $`\mathbb E[O_{ij}O_{lm}]=\delta_{il}\delta_{jm}/k`$，$`\mathbb E\mathrm{Tr}(AO)^2=\sum_{ij}A_{ji}^2/k`$。

3. 块对角算子的迹与迹范数都是各块之和，且 $`|X\oplus X'|=|X|\oplus|X'|`$；
   $`\lambda_{\min}(A\oplus A')=\min(\lambda_{\min}(A),\lambda_{\min}(A'))`$，而 $`\Delta`$ 随 $`\lambda_{\min}`$ 递减。
   所以若 $`(X,Y)`$ 在 $`\delta`$ 处容许，$`(X',Y')`$ 在 $`\delta'\le\delta`$ 处容许，则 $`(X\oplus X',Y\oplus Y')`$ 在
   $`\delta`$ 处容许，且 $`z=z(X,Y)+z(X',Y')`$；对两个数对以及 $`\delta'\le\delta`$ 取上确界即可。

</details>

*计算：* 习题 1 与 2 的数值验证——[教程，习题](../../tutorial/zh/11_optimisation.md#x-11)。

**下一章：** [第 12 章](12_equivariant_learning.md)讨论等变机器学习：与对称性可交换的层是
交换子（第 2 章）的元素，它们由同样的批量化、可微分的算子构造而成。
