[English](../11_optimisation.md) · [Español](../es/11_optimisation.md) · [Français](../fr/11_optimisation.md) · [Deutsch](../de/11_optimisation.md) · **中文** · [日本語](../ja/11_optimisation.md) · [Русский](../ru/11_optimisation.md)

# 11. 优化与迹不等式

迹不等式是关于算子泛函极值的命题。3.7 节给出了一个看似合理却不成立的迹不等式（命题 3.9）；一次采样运行（[教程，3.7 节](../../tutorial/zh/03_traces_and_states.md#s-3-7)）在数值上没有找到反例。本章解释为什么采样在高维中会失败，并以在精确参数化的约束集上的优化
取而代之；在极值有闭式表达的情形，我们给出证明，使优化器的结果能与定理相对照。

**你将学到**如何对酉算子以及具有给定对比度的谱进行参数化，以及为什么每个容许算子都能被取到；
冯·诺依曼迹不等式及由它给出的精确极值；一个带有取等条件的迹不等式；为什么 Haar 样本在均方意义下远离
极值；以及多起点搜索究竟证明了什么。本章的数值检验以及本库的优化器见
[教程 11](../../tutorial/zh/11_optimisation.md)。

在本章中，$\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$，$D=\sum_ck_c\ge2$，$\mathrm{Tr}$ 是
blunt 迹（3.1 节），$\|A\|_1=\sum_c\|A_c\|_1$ 是把 $A$ 看作
$M_D(\mathbb C)$ 中块对角矩阵时的迹范数。$A\in M_k(\mathbb C)$ 的奇异值为 $s_1(A)\ge\dots\ge s_k(A)$；$x^\downarrow$ 与
$x^\uparrow$ 分别是实向量 $x$ 的递减重排与递增重排。$\Delta$ 是 Michelson
对比度（3.6 节），对 $\delta\in[0,1)$ 令
$$
\ell(\delta)=\frac{1-\delta}{1+\delta},\qquad 1-\ell(\delta)=\frac{2\delta}{1+\delta};
$$
$\lambda_{\max}=1$ 的正算子的对比度为 $\ell(\lambda_{\min})$，因此它的对比度为 $\delta$ 当且仅当
$\lambda_{\min}=\ell(\delta)$，因为 $\ell$ 是 $[0,1]$ 上严格递减的对合。例子是配套论文中的三个泛函
（*From sampling to optimisation* 一节）：

| | 泛函 | 变量 |
|---|---|---|
| 实验 1 | $z=\lvert\mathrm{Tr}(XUY)\rvert-\mathrm{Tr}(XY)$ | $X,Y\ge0$，$U\in\mathcal M$ 为酉算子 |
| 实验 2 | $z=\mathrm{Tr}\lvert XY\rvert-\mathrm{Tr}(\lvert X\rvert\,\lvert Y\rvert)$ | $X=UX_0$，$Y=VY_0$，$X_0,Y_0\ge0$，$U,V$ 为酉算子 |
| 实验 3 | $z=\mathrm{Tr}(Y\lvert X\rvert Y)-\mathrm{Tr}\lvert YXY\rvert$ | $X=X^*$，$Y\ge0$ |

其中 $|X|$、$|Y|$ 归一化为 $\lambda_{\max}=1$，它们的对比度作为参数。**包络**为
$\overline z(\delta)=\sup\{z:\Delta(|X|)=\delta\}$ 与 $\underline z(\delta)=\inf\{z:\Delta(|X|)=\delta\}$，其余
变量均自由变化。

<a id="s-11-1"></a>

## 11.1 约束集及其参数化

**定义。** 对 $\delta\in[0,1)$ 与 $\ell=\ell(\delta)$，**对比度壳层**为
$$
\mathcal P_\delta=\{X\in\mathcal M:X\ge0,\ \lambda_{\max}(X)=1,\ \Delta(X)=\delta\}
=\{X=X^*:\mathrm{spec}\,X\subseteq[\ell,1],\ 1,\ell\in\mathrm{spec}\,X\}.
$$
$\mathcal M$ 的**槽位**是 $D$ 个数对 $(c,i)$，$i\le k_c$；对 $\lambda\in\mathbb R^D$，
$\mathrm{diag}(\lambda)=\bigoplus_c\mathrm{diag}(\lambda_{c,1},\dots,\lambda_{c,k_c})$。

**参数化。** 对给定的对比度 $\delta$，搜索对每个起点抽取两个不同的槽位
$s_\pm=(c_\pm,i_\pm)$（先在 $D$ 个槽位中均匀地抽取 $s_+$，再在其余 $D-1$ 个槽位中均匀地抽取 $s_-$，因此
每个由不同槽位构成的有序对的概率都是 $1/(D(D-1))$）和一个 Haar 随机酉算子 $U_0\in\mathcal M$
（引理 11.10），并把 $\theta\in\mathbb R^D$ 与
$G=\bigoplus_cG_c$ 映为
$$
X(\theta,G)=U\mathrm{diag}(\lambda)U^*,\qquad U=U_0\,e^{G-G^*},\qquad
\lambda_{s_+}=1,\quad\lambda_{s_-}=\ell,\quad\lambda_s=\ell+(1-\ell)\tfrac{1+\sin\theta_s}2\ \text{otherwise}.
$$
单独的酉算子则由 $G\mapsto U_0e^{G-G^*}$ 参数化。与 sigmoid 不同，$(1+\sin\theta)/2$ 在有限的 $\theta$ 处
就能取到 $0$ 和 $1$，这一点很方便，因为本章的极值谱（推论 11.7、命题 11.14）位于 $\ell$ 和 $1$。
$[a,b]$ 中的自由对比度只是多一个参数，$\delta=a+(b-a)(1+\sin\eta)/2$，例如取 $[a,b]=[0,0.999]$。在这一限制下取到的值在没有限制时也能取到，因此它仍然是无限制上确界的下界
（无限制下确界的上界）。

*库中：* 这一参数化就是 `optimize.PositiveParam`，而 $G\mapsto U_0e^{G-G^*}$ 就是 `optimize.UnitaryParam`——[教程，11.1 节](../../tutorial/zh/11_optimisation.md#s-11-1)。

<a id="r-11-1"></a>

**引理 11.1.** 对每个酉算子 $U_0\in\mathcal M$，$G\mapsto U_0e^{G-G^*}$ 映满 $\mathcal M$ 的
酉群。

*证明。* 设 $U$ 是酉算子。由定理 1.3（逐块应用）与推论 1.4(2)，
$U_0^*U=\bigoplus_cW_c\mathrm{diag}(e^{i\varphi_{c,j}})_jW_c^*$，其中 $W_c$ 为酉矩阵，$\varphi_{c,j}\in\mathbb R$。于是
$S=\bigoplus_cW_c\mathrm{diag}(i\varphi_{c,j})_jW_c^*$ 满足 $S^*=-S$ 与 $e^S=U_0^*U$（级数
$\sum_nS^n/n!$ 在特征基中逐项计算），而 $G=S/2$ 给出 $G-G^*=S$。$\square$

<a id="r-11-2"></a>

**命题 11.2（参数化是满射）。** 对固定的 $\delta$、$s_\pm$ 与 $U_0$，
$(\theta,G)\mapsto X(\theta,G)$ 的像是 $\{X\in\mathcal P_\delta:1\in\mathrm{spec}\,X_{c_+},\ \ell\in\mathrm{spec}\,X_{c_-}\}$，
并且 $\mathcal P_\delta$ 是这些像对所有不同槽位对的并。

*证明。* 像包含在该集合中：$X$ 的特征值 $\lambda_s\in[\ell,1]$，且 $\lambda_{s_\pm}$ 位于块
$c_\pm$ 中。反之，对角化 $X_c=W_c\mathrm{diag}(\mu_c)W_c^*$ 并对每个 $\mu_c$ 作置换（把 $W_c$ 换成
$W_c\Pi_c$，$\Pi_c$ 为置换矩阵），使 $\mu_{s_+}=1$ 且 $\mu_{s_-}=\ell$。这是可行的：若 $c_+=c_-$
且 $\delta>0$，则 $1\ne\ell$ 是同一块的两个特征值，占据不同的位置；若 $\delta=0$，
所有特征值都是 $1$。若 $\ell<1$，对其余槽位令 $\theta_s=\arcsin\big(2(\mu_s-\ell)/(1-\ell)-1\big)$，
使 $\lambda=\mu$（若 $\ell=1$，$\theta$ 任意）。引理 11.1 给出满足 $U_0e^{G-G^*}=\bigoplus_cW_c$ 的 $G$。
*并：* 对 $X\in\mathcal P_\delta$，选取含特征值 $1$、$\ell$ 的块 $c_\pm$，并在其中选取不同的槽位；
若 $c_+=c_-$ 且 $\delta>0$，该块有两个不同的特征值，所以 $k_{c_+}\ge2$；若 $\delta=0$，则 $X=1$，
任取两个槽位即可。$\square$

*计算：* 构造性地检验命题 11.2：证明中的参数重现一个目标算子——[教程，11.1 节](../../tutorial/zh/11_optimisation.md#s-11-1)。

由于每一对槽位都以正概率被抽到，壳层中的每一点都以正概率位于某个起点的像中，并且（在精确算术中）
约束在每一步迭代中都成立。
对自伴的 $X$（实验 3），使用 $X=U\mathrm{diag}(\sigma\lambda)U^*$，其中符号 $\sigma\in\{\pm1\}^D$ 对每个起点抽取
（在每个块中，负号的个数在 $\{0,\dots,k_c\}$ 中均匀分布，其位置均匀随机，因此每种符号模式都有正概率），
并且只通过 11.5 节中的符号翻转改变。于是
$|X|=U\mathrm{diag}(\lambda)U^*$（定理 1.5(e)，因为 $\lambda>0$），并且对所有符号与槽位取并得到
$\{X=X^*:|X|\in\mathcal P_\delta\}$：对这样的 $X$，对角化 $X_c=W_c\mathrm{diag}(\nu_c)W_c^*$，于是
$|X|_c=W_c\mathrm{diag}(|\nu_c|)W_c^*$，其中 $|\nu_{c,i}|\ge\ell>0$，再对 $\mu=|\nu|$ 与
$\sigma_s=\nu_s/|\nu_s|$ 应用命题 11.2 的证明。只有 $S=G-G^*$ 起作用。

*库中：* `optimize.SelfAdjointParam`——[教程，11.1 节](../../tutorial/zh/11_optimisation.md#s-11-1)。

坐标卡 $S\mapsto U_0e^S$ 只在 $S=0$ 附近性质良好。令 $\mathfrak u=\{S\in\mathcal M:S^*=-S\}$，
$\mathcal M_h=\{H\in\mathcal M:H^*=H\}$，它们是实子空间，且 $\mathcal M=\mathfrak u\oplus\mathcal M_h$
（$X=\frac{X-X^*}2+\frac{X+X^*}2$）。*在 $0$ 附近，$S\mapsto U_0e^S$ 是 $\mathfrak u$ 的一个开球 $N\ni0$ 到
酉群中 $U_0$ 的一个邻域上的同胚。* 事实上，$\Psi(S,H)=U_0e^Se^H$ 是 $C^1$ 映射
$\mathfrak u\times\mathcal M_h\to\mathcal M$（级数 $\sum_{m,n\ge0}U_0S^mH^n/(m!\,n!)$ 以及其各项关于实坐标的
偏导数构成的级数在有界集上一致收敛），而
$e^{tS}e^{tH}=1+t(S+H)+O(t^2)$ 表明它在 $0$ 处的微分为 $(S,H)\mapsto U_0(S+H)$，这是到
$\mathcal M$ 上的双射。由反函数定理，存在 $\mathfrak u$ 中的开球 $N\ni0$ 与 $\mathcal M_h$ 中的开球 $N_h\ni0$，
使得 $\Psi$ 是 $N\times N_h$ 到开集 $O\subseteq\mathcal M$ 上的同胚。若
$\Psi(S,H)$ 是酉的，则 $e^H=e^{-S}U_0^*\Psi(S,H)$ 也是酉的（在正规矩阵 $S$ 的特征基中——其特征值为纯虚数——
$e^{\pm S}$ 是对角的，对角元的模为 $1$，并且二者互逆）；但 $e^H$ 具有正特征值 $e^{h_i}$，$h_i$ 为 $H$ 的特征值，
所以 $e^{h_i}=1$（推论 1.4(2)），从而 $H=0$。
因此 $O\cap U(\mathcal M)=U_0e^N$，其中 $U(\mathcal M)$ 是 $\mathcal M$ 的酉群：这是 $U_0$ 在
$U(\mathcal M)$ 中的一个邻域，而 $S\mapsto U_0e^S$ 是连续双射 $N\to U_0e^N$，其逆映射，即 $\Psi^{-1}$ 的
第一个分量，是连续的。远离 $0$ 时该坐标卡不是单射：对每个投影 $P\in\mathcal M$ 有 $e^{2\pi iP}=1=e^0$。
因此搜索在每一轮之后都把坐标卡**重新居中**（$U_0\leftarrow U_0e^{G-G^*}$，$G\leftarrow0$）。另外两个坐标卡
在 $0$ 处具有相同的微分 $S\mapsto U_0S$。Cayley 坐标卡
$S\mapsto U_0C(-S/2)$（习题 1）：$(1+tS/2)(1-tS/2)^{-1}=(1+tS/2)(1+tS/2+O(t^2))=1+tS+O(t^2)$。极分解回缩，
即把 $S$ 映为 $A=U_0(1+S)$ 的酉因子（定理 1.8）：$A^*A=(1-S)(1+S)=1+S^*S\ge1$ 可逆，所以 $|A|$
可逆，且酉因子为 $A|A|^{-1}$；把 $S$ 换成 $tS$，有 $|A|^{-1}=(1+t^2S^*S)^{-1/2}=1+O(t^2)$
（在 $S^*S$ 的特征基中计算），所以该酉因子为 $U_0(1+tS)+O(t^2)$。

> **物理意义**（物理解释，而非本课程的定理）。Michelson 对比度就是干涉条纹的可见度 $(I_{\max}-I_{\min})/(I_{\max}+I_{\min})$。
> 对 $\lambda_{\max}=1$ 的 $0\le X\le1$——滤光片、衰减器、探测器效率——
> 它确定了动态范围 $\lambda_{\min}/\lambda_{\max}$，即条件数的倒数。

<a id="s-11-2"></a>

## 11.2 冯·诺依曼迹不等式

元素 $c_{ij}\ge0$ 且所有行和与列和都为 $1$ 的矩阵 $(c_{ij})$ 称为**双随机**矩阵；对酉矩阵 $G$，
$(|G_{ij}|^2)$ 与 $(|G_{ji}|^2)$ 都是双随机的，因为 $G$ 的行和列都是单位向量。

<a id="r-11-3"></a>

**引理 11.3（重排）。** 设 $p_1\ge\dots\ge p_k$，$q_1\ge\dots\ge q_k$ 为实数，$(c_{ij})$ 为双
随机矩阵。则 $\sum_ip_iq_{k+1-i}\le\sum_{i,j}p_iq_jc_{ij}\le\sum_ip_iq_i$。

*证明。* *上界。* 对 $m\le k$ 令 $r_j=\sum_{i\le m}c_{ij}\in[0,1]$，于是 $\sum_jr_j=m$。由于当
$j\le m$ 时 $q_j\ge q_m$，当 $j>m$ 时 $q_j\le q_m$，
$$
\sum_{j\le m}q_j-\sum_jq_jr_j=\sum_{j\le m}q_j(1-r_j)-\sum_{j>m}q_jr_j\ge q_m\Big(m-\sum_jr_j\Big)=0 .
$$
记 $x_i=\sum_jq_jc_{ij}$，上式即 $\sum_{i\le m}x_i\le\sum_{j\le m}q_j$，且 $m=k$ 时取等号。由 Abel
求和以及 $p_m-p_{m+1}\ge0$，
$$
\sum_ip_ix_i=p_k\sum_{i\le k}x_i+\sum_{m<k}(p_m-p_{m+1})\sum_{i\le m}x_i
\le p_k\sum_{j\le k}q_j+\sum_{m<k}(p_m-p_{m+1})\sum_{j\le m}q_j=\sum_ip_iq_i .
$$
*下界。* 对 $p$、递减的 $q'_j=-q_{k+1-j}$ 以及双随机矩阵
$c'_{ij}=c_{i,k+1-j}$ 应用上界：$-\sum_{ij}p_iq_jc_{ij}=\sum_{ij}p_iq'_jc'_{ij}\le\sum_ip_iq'_i=-\sum_ip_iq_{k+1-i}$。$\square$

每个 $A\in M_k(\mathbb C)$ 都有**奇异值分解** $A=W\Sigma V^*$，$W,V$ 为酉矩阵，
$\Sigma=\mathrm{diag}(s_i(A))$：写 $A=W'|A|$（定理 1.8），$|A|=V\Sigma V^*$（定理 1.3），$W=W'V$。

<a id="r-11-4"></a>

**定理 11.4（冯·诺依曼迹不等式）。** 对 $A,B\in M_k(\mathbb C)$，
$|\mathrm{Tr}(AB)|\le\sum_is_i(A)\,s_i(B)$。

*证明。* 设 $A=W_A\Sigma_AV_A^*$，$B=W_B\Sigma_BV_B^*$，$a_i=s_i(A)$，$b_j=s_j(B)$，$G=V_A^*W_B$，$H=V_B^*W_A$。
由循环性，$\mathrm{Tr}(AB)=\mathrm{Tr}(\Sigma_AG\Sigma_BH)=\sum_{i,j}a_iG_{ij}b_jH_{ji}$，因此由
$|g||h|\le\frac12(|g|^2+|h|^2)$ 与 $a_ib_j\ge0$，
$$
|\mathrm{Tr}(AB)|\le\sum_{i,j}a_ib_j\,d_{ij},\qquad d_{ij}=\tfrac12\big(|G_{ij}|^2+|H_{ji}|^2\big).
$$
$(d_{ij})$ 是双随机的，引理 11.3 给出 $\sum_{ij}a_ib_jd_{ij}\le\sum_ia_ib_i$。$\square$

<a id="r-11-5"></a>

**推论 11.5.** (a) 对特征值为 $p,q$ 的自伴矩阵 $P,Q\in M_k(\mathbb C)$：
$\sum_ip^\downarrow_iq^\uparrow_i\le\mathrm{Tr}(PQ)\le\sum_ip^\downarrow_iq^\downarrow_i$。
(b) 对 $P,Q\ge0$：$\mathrm{Tr}(PQ)\le\|PQ\|_1\le\sum_ip^\downarrow_iq^\downarrow_i$。

*证明。* (a) 与定理 3.5 的证明一样，$\mathrm{Tr}(PQ)=\sum_{ij}p_iq_jc_{ij}$，其中对按递减顺序排列的特征基，
$c_{ij}=|\langle u_i,v_j\rangle|^2$ 是双随机的；应用引理 11.3（其下界为
$\sum_ip^\downarrow_iq^\downarrow_{k+1-i}=\sum_ip^\downarrow_iq^\uparrow_i$）。(b) 第一个不等式就是
定理 1.9。由定理 1.9 与定理 11.4，$\|PQ\|_1=\max_U|\mathrm{Tr}(P\,QU)|\le\sum_is_i(P)s_i(QU)$，而
$s_i(QU)=q^\downarrow_i$，因为 $|QU|=(U^*Q^2U)^{1/2}=U^*QU$（定理 1.6）。$\square$

*计算：* 在随机矩阵上检验定理 11.4 与推论 11.5——[教程，11.2 节](../../tutorial/zh/11_optimisation.md#s-11-2)。

因此对 $M_k(\mathbb C)$ 中谱为 $x,y$ 的 $X,Y\ge0$，命题 3.9 中的违背量 $\|YX\|_1-\mathrm{Tr}(XY)$
至多为 $\sum_ix^\downarrow_i(y^\downarrow_i-y^\uparrow_i)$（对 $(Y,X)$ 应用推论 11.5(b)，再用 (a) 的下界）。
对实验 2，这个数恰好就是精确答案：

<a id="r-11-6"></a>

**命题 11.6（固定谱时的实验 2）。** 在实验 2 中，$|X|=X_0$，$|Y|=Y_0$，且
$z=\|X_0VY_0\|_1-\mathrm{Tr}(X_0Y_0)$。若 $(X_0)_c$ 的谱 $x_c$ 与 $(Y_0)_c$ 的谱 $y_c$ 固定，而
特征基与 $V$ 变化，则
$$
\sup z=-\inf z=\sum_c\sum_ix^\downarrow_{c,i}\big(y^\downarrow_{c,i}-y^\uparrow_{c,i}\big),
$$
并且二者都能取到。

*证明。* $X^*X=X_0U^*UX_0=X_0^2$，所以 $|X|=X_0$（定理 1.6）；同理 $|Y|=Y_0$，而 $|UA|=|A|$ 给出
$\mathrm{Tr}|XY|=\|X_0VY_0\|_1$。所有变量都按块作用，且 $z$ 是各块之和，故不妨设 $C=1$。
令 $Y'=VY_0V^*\ge0$（谱为 $y$）；由于 $|AV|=V^*|A|V$，$\|X_0VY_0\|_1=\|X_0Y'V\|_1=\|X_0Y'\|_1$。对数对
$(X_0,Y')$ 与 $(X_0,Y_0)$ 应用推论 11.5，并利用 $\|X_0Y'\|_1\ge\mathrm{Tr}(X_0Y')$，得到
$$
\sum_ix^\downarrow_iy^\uparrow_i-\sum_ix^\downarrow_iy^\downarrow_i\le\|X_0Y'\|_1-\mathrm{Tr}(X_0Y_0)\le\sum_ix^\downarrow_iy^\downarrow_i-\sum_ix^\downarrow_iy^\uparrow_i .
$$
设 $J$ 为反序置换，$Je_i=e_{k+1-i}$；则 $J\mathrm{diag}(w)=\mathrm{diag}(w_k,\dots,w_1)J$，并且
$\mathrm{diag}(w)J$ 的迹范数为 $\sum_i|w_i|$。取 $X_0=\mathrm{diag}(x^\downarrow)$ 与 $V=J$：
对 $Y_0=\mathrm{diag}(y^\uparrow)$，$X_0VY_0=\mathrm{diag}(x^\downarrow_iy^\downarrow_i)J$ 且
$\mathrm{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\uparrow_i$（上界）；对 $Y_0=\mathrm{diag}(y^\downarrow)$，
$X_0VY_0=\mathrm{diag}(x^\downarrow_iy^\uparrow_i)J$ 且 $\mathrm{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\downarrow_i$
（下界）。$\square$

<a id="r-11-7"></a>

**推论 11.7（实验 2 的包络）。** 设 $\ell_X,\ell_Y\in(0,1]$，且某个 $k_c\ge2$。在所有满足
$\lambda_{\max}=1$、$\lambda_{\min}(X_0)=\ell_X$、$\lambda_{\min}(Y_0)=\ell_Y$ 的 $X_0,Y_0\ge0$ 以及所有酉算子上，
$$
\sup z=-\inf z=\sum_c\Big\lfloor\frac{k_c}2\Big\rfloor(1-\ell_X)(1-\ell_Y).
$$

*证明。* 在一个块中记 $x_i=x^\downarrow_{c,i}$，$y_i=y^\downarrow_{c,i}$，$k=k_c$，于是 $y^\uparrow_{c,i}=y_{k+1-i}$。
把 $i$ 与 $k+1-i$ 配对，$x_i(y_i-y_{k+1-i})+x_{k+1-i}(y_{k+1-i}-y_i)=(x_i-x_{k+1-i})(y_i-y_{k+1-i})$，而当
$k$ 为奇数时，中间的指标贡献 $x_i(y_i-y_i)=0$。因此
$$
\sum_ix_i(y_i-y_{k+1-i})=\sum_{i\le\lfloor k/2\rfloor}(x_i-x_{k+1-i})(y_i-y_{k+1-i})\le\Big\lfloor\frac k2\Big\rfloor(1-\ell_X)(1-\ell_Y),
$$
因为 $x_i\in[\ell_X,1]$，$y_i\in[\ell_Y,1]$ 都是递减的，所以两个因子分别位于 $[0,1-\ell_X]$ 与 $[0,1-\ell_Y]$ 中。
若在每个块中，最大的 $\lfloor k_c/2\rfloor$ 个特征值为 $1$，最小的 $\lfloor k_c/2\rfloor$ 个为 $\ell_X$（相应地为
$\ell_Y$），中间那个（$k_c$ 为奇数时）在 $[\ell_X,1]$（相应地 $[\ell_Y,1]$）中任意，则等号成立；由于某个 $k_c\ge2$，这些谱的全局极值为 $1$ 和 $\ell_X$（$\ell_Y$）。由命题 11.6
即得结论。$\square$

这里对 $\Delta(|X|)=\delta$ 有 $1-\ell_X=2\delta/(1+\delta)$；当 $k_c$ 为奇数时，中间的特征值无关紧要。

*计算：* $k=3$、$C=2$ 时命题 11.6 的极值构型与推论 11.7 的包络——[教程，11.2 节](../../tutorial/zh/11_optimisation.md#s-11-2)；多起点搜索求得的包络——[教程，11.5 节](../../tutorial/zh/11_optimisation.md#s-11-5)。

<a id="s-11-3"></a>

## 11.3 一个带有取等条件的迹不等式

对 $X=X^*$，令 $X_\pm=f_\pm(X)$，其中 $f_+(t)=\max(t,0)$，$f_-(t)=\max(-t,0)$。则 $X_\pm\ge0$，$X=X_+-X_-$，
$|X|=X_++X_-$，且 $X_+X_-=0$（定理 1.5(c)；对 $|X|=(X^2)^{1/2}=h(X)$，其中
$h(t)=\sqrt{t^2}=|t|=f_+(t)+f_-(t)$，用定理 1.5(e)；正性由推论 1.4(3) 得到）。

<a id="r-11-8"></a>

**引理 11.8.** 对 $M_n(\mathbb C)$ 中的 $A,B\ge0$，$\|A-B\|_1\le\mathrm{Tr}\,A+\mathrm{Tr}\,B$，等号成立当且仅当
$AB=0$。

*证明。* 令 $H=A-B$，$Q=g(H)$，其中 $g$ 在 $(0,\infty)$ 上为 $1$，在其他地方为 $0$；$Q$ 是投影，因为 $g=\bar g=g^2$
（定理 1.5(c)）。由于 $|t|=t(2g(t)-1)$，有 $|H|=H(2Q-1)$（定理 1.5(c)、(e)），再对投影 $R=Q,1-Q$ 利用
$\mathrm{Tr}(AR)=\mathrm{Tr}(RAR)$（循环性），得到
$$
\|H\|_1=\mathrm{Tr}\big(A(2Q-1)\big)-\mathrm{Tr}\big(B(2Q-1)\big)
=\mathrm{Tr}\,A+\mathrm{Tr}\,B-2\mathrm{Tr}\big((1-Q)A(1-Q)\big)-2\mathrm{Tr}(QBQ).
$$
被减去的迹都 $\ge0$，而 $\mathrm{Tr}(RAR)=\|A^{1/2}R\|_2^2$ 为零当且仅当 $AR=0$；所以等号成立
当且仅当 $A(1-Q)=0$ 且 $BQ=0$，此时 $AB=AQB=A(BQ)^*=0$。反之，若 $AB=0$，则 $BA=0$，
$(A-B)^2=A^2+B^2=(A+B)^2$，由定理 1.6 得 $|A-B|=A+B$。$\square$

<a id="r-11-9"></a>

**命题 11.9（实验 3）。** 对 $\mathcal M$ 中的 $X=X^*$ 与 $Y\ge0$，
$\mathrm{Tr}|YXY|\le\mathrm{Tr}(Y|X|Y)$，等号成立当且仅当 $YX_+Y^2X_-Y=0$，特别地，当 $XY=YX$ 时成立。

*证明。* 在 $M_D(\mathbb C)$ 中，$A=YX_+Y$ 与 $B=YX_-Y$ 是正的（取 $T=X_\pm^{1/2}Y$ 应用推论 1.7），$YXY=A-B$，
且 $\mathrm{Tr}(Y|X|Y)=\mathrm{Tr}\,A+\mathrm{Tr}\,B$；应用引理 11.8。若 $XY=YX$，则 $Y$ 与
$X_\pm$ 可交换（定理 1.5(d)），从而 $AB=Y^4X_+X_-=0$。$\square$

因此在实验 3 中 $z\ge0$，且 $\underline z=0$（取 $Y=1$）。取等条件比可交换性更弱：对
$X=\mathrm{diag}(1,-1,0)$ 以及行为 $(1,0,a)$、$(0,1,0)$、$(a,0,1)$（$|a|<1$）的 $Y$，
$AB=(Ye_1)\langle Ye_1,Ye_2\rangle(Ye_2)^*=0$，而当 $a\ne0$ 时 $XY\ne YX$。

*计算：* 在随机矩阵对上检验命题 11.9，以及不可交换而等号成立的例子——[教程，11.3 节](../../tutorial/zh/11_optimisation.md#s-11-3)。

<a id="s-11-4"></a>

## 11.4 采样与优化

与其最大值 $\|A\|_1$ 相比，$|\mathrm{Tr}(AU)|$ 的典型大小是多少？我们把 $U$ 取为
**Ginibre 矩阵** $Z$（元素相互独立，密度为
$\pi^{-1}e^{-|z|^2}$；Mezzadri 的方法）的 QR 分解 $Z=QR$（$R_{ii}>0$）中的因子 $Q$，在每个块中独立抽取。

<a id="r-11-10"></a>

**引理 11.10.** $Q$ 在 $U(k)$ 上服从 Haar 分布；特别地，对每个固定的酉矩阵 $V\in M_k(\mathbb C)$，$VQ$ 与
$Q$ 同分布。各块相互独立的酉算子 $\bigoplus_cQ_c$ 在 $\mathcal M$ 的酉群上服从 Haar 分布。

*证明。* 这就是定理 5.3（左不变性是其证明中的等变性步骤；该分布由唯一性引理 5.1 认定为 Haar 测度，
其测度论步骤是推论 A.2）与推论 5.4。$\square$

*库中：* 这一采样器就是 `TypeIAlgebra.random_unitary_operator`——[教程，11.4 节](../../tutorial/zh/11_optimisation.md#s-11-4)。

<a id="r-11-11"></a>

**命题 11.11.** 对 Haar 分布的 $U\in U(k)$ 与 $A\in M_k(\mathbb C)$：$\mathbb E\mathrm{Tr}(AU)=0$，且
$\mathbb E|\mathrm{Tr}(AU)|^2=\|A\|_2^2/k$。

*证明。* 由推论 5.5，$\mathbb EU_{ij}=0$。对 $X=|e_j\rangle\langle e_m|$ 应用命题 5.7（取 $C=1$），得到
一阶 Weingarten 公式
$\mathbb E[U_{ij}\bar U_{lm}]=\mathbb E(UXU^*)_{il}=\delta_{il}\delta_{jm}/k$；对 $k\ge2$，它也可由
二阶公式（定理 5.10）得到：令 $i_2=i'_2$，$j_2=j'_2$ 并对二者求和，因为
$\sum_{i_2,j_2}|U_{i_2j_2}|^2=k$。利用 $\mathrm{Tr}(AU)=\sum_{ij}A_{ji}U_{ij}$，
$\mathbb E|\mathrm{Tr}(AU)|^2=\sum_{ijlm}A_{ji}\bar A_{ml}\mathbb E[U_{ij}\bar U_{lm}]=\frac1k\sum_{ij}|A_{ji}|^2$。$\square$

<a id="r-11-12"></a>

**推论 11.12.** 设 $M_k(\mathbb C)$ 中的 $X,Y\ge0$ 的谱含于 $[a,b]$，$a>0$，$U$ 服从 Haar 分布。则
$$
\frac{\big(\mathbb E|\mathrm{Tr}(XUY)|^2\big)^{1/2}}{\max_V|\mathrm{Tr}(XVY)|}=\frac{\|YX\|_2}{\sqrt k\,\|YX\|_1}\le\min\Big(\frac1{\sqrt k},\frac{b^2}{a^2k}\Big),
\qquad
\mathbb P\big(|\mathrm{Tr}(XUY)|\ge\mathrm{Tr}(XY)\big)\le\frac{b^4}{a^4k^2}.
$$

*证明。* 其中的等式由命题 11.11（取 $A=YX$，因为 $\mathrm{Tr}(XUY)=\mathrm{Tr}(AU)$）与命题 3.9 得到。
由于 $\|A\|_2^2=\sum_is_i^2\le(\sum_is_i)^2=\|A\|_1^2$（1.4 节），该比值 $\le1/\sqrt k$。
由于 $\|Av\|=\||A|v\|$（1.4 节），对角化 $|A|$ 给出 $s_k(A)=\min_{\|v\|=1}\|Av\|$ 与
$s_1(A)=\max_{\|v\|=1}\|Av\|$；而对单位向量 $v$ 有 $\|Xv\|^2=\langle v,X^2v\rangle\in[a^2,b^2]$（对 $X^2-a^2$ 与
$b^2-X^2$ 应用推论 1.4(3)），对 $Y$ 同理。因此 $a^2\le s_k(A)\le s_1(A)\le b^2$，$\|A\|_2^2\le kb^4$，且
$\|A\|_1\ge ka^2$。最后 $\mathrm{Tr}(XY)=\mathrm{Tr}(X^{1/2}YX^{1/2})\ge a\mathrm{Tr}\,X\ge ka^2$，而对
$|\mathrm{Tr}(XUY)|^2$ 应用 Markov 不等式，可把该概率界定为 $(kb^4/k)/(ka^2)^2$。$\square$

因此，$N$ 个 Haar 样本找到 $|\mathrm{Tr}(XUY)|\le\mathrm{Tr}(XY)$ 的反例的概率至多为
$Nb^4/(a^4k^2)$（$N$ 个事件的并），尽管对每一对不可交换的矩阵这样的反例都存在（命题 3.9）。这个界
只有在 $k>\sqrt N\,b^2/a^2$ 时才有信息量（对 $b/a=6$ 与 $N=2000$：$k>1609$）。

*计算：* 用 Monte Carlo 检验命题 11.11，并对 $k\le32$ 把推论 11.12 中的比值与 $2000$ 个样本中的最好者比较——[教程，11.4 节](../../tutorial/zh/11_optimisation.md#s-11-4)。

<a id="r-11-13"></a>

**引理 11.13.** 若 $k\ge2$，则对每个 $A\in M_k(\mathbb C)$ 有 $\min_U|\mathrm{Tr}(AU)|=0$。

*证明。* 设 $A=W\Sigma V^*$，$P$ 为循环移位 $Pe_i=e_{(i\bmod k)+1}$，当 $k\ge2$ 时其对角线为零。对
$U=VPW^*$，$\mathrm{Tr}(AU)=\mathrm{Tr}(\Sigma P)=0$。$\square$

所以 $\inf_Uz=-\mathrm{Tr}(XY)$，$\sup_Uz=\|YX\|_1-\mathrm{Tr}(XY)$（命题 3.9）。由命题 11.11，$|\mathrm{Tr}(XUY)|$ 的均方根为 $\|YX\|_2/\sqrt k$；利用推论 11.12 的证明中的
界（谱含于 $[a,b]$），$\mathbb E|\mathrm{Tr}(XUY)|^2\le b^4$ 且 $\mathrm{Tr}(XY)\ge ka^2$，所以 Markov 不等式给出
$\mathbb P\big(z\ge-(1-\varepsilon)\mathrm{Tr}(XY)\big)\le b^4/(\varepsilon^2a^4k^2)$（$\varepsilon>0$）：对大的 $k$，
$z$ 以高概率接近其下确界。

*计算：* 用梯度上升与梯度下降求关于 $U$ 的两个极值，并与 $20\,000$ 个 Haar 样本比较——[教程，11.4 节](../../tutorial/zh/11_optimisation.md#s-11-4)。

> **物理意义**（物理解释，而非本课程的定理）。推论 11.12 是统计力学所称的**典型性**的一个较弱但可证明的
> 实例：在高维中，人们预期几乎每个态或演化都给出宏观量几乎相同的值。一个精确（最优）的不等式由精细调节的
> 构型决定——对齐的特征基、两能级的谱——而随机酉算子（随机线路、随机散射矩阵）很少产生这样的构型。

<a id="s-11-5"></a>

## 11.5 包络的多起点搜索

**多起点搜索**从许多独立的起点出发，分轮运行梯度方法，每个起点有自己的槽位、符号、$U_0$ 与对比度；
每轮之后对酉坐标卡重新居中，并对每个起点返回目标函数在各迭代点处的最优值。自伴参数化的符号是离散的：
在两轮之间，通过在随机槽位上的单个符号翻转来改变它们，每次翻转若能改进目标函数则予以保留。

*库中：* `optimize.extremize`——[教程，11.5 节](../../tutorial/zh/11_optimisation.md#s-11-5)。

在精确算术中，每一步迭代都是容许的（命题 11.2 及其后的说明），所以返回的每个值都是 $z$ 在某个容许点处的值；
在浮点运算中，它在舍入误差范围内等于该值。因此最大化给出上确界的
**下界**，最小化给出下确界的上界，仅此而已——局部搜索可能停在局部极值处。另一个方向需要证明，
如推论 11.7 与命题 11.9，或者凸松弛（11.6 节）。计算得到的包络要接受以下检验：(i) 与闭式表达比较；(ii) 与
手工构造的构型比较，例如较小解的直和，包络不得低于它们；(iii) 与
结构化样本（两能级谱、置换基）比较，它们的值是另外一些可取到的值。

<a id="r-11-14"></a>

**命题 11.14（实验 1 的下包络）。** 设所有 $k_c\ge2$。在 $X\in\mathcal P_\delta$、满足
$\lambda_{\max}(Y)=1$ 的 $Y\ge0$ 以及酉算子 $U\in\mathcal M$ 上，$\ \min z=-\big(D-2\delta/(1+\delta)\big)$。

*证明。* $Y\le1$ 给出 $X^{1/2}YX^{1/2}\le X$，所以 $\mathrm{Tr}(XY)\le\mathrm{Tr}\,X\le(D-1)+\ell(\delta)$，因为
$X$ 的一个特征值为 $\ell(\delta)$，其余都 $\le1$；因此 $z\ge-\mathrm{Tr}(XY)\ge-(D-2\delta/(1+\delta))$。
当 $X$ 的特征值为一个 $\ell(\delta)$、其余都为 $1$，$Y=1$，且 $U=\bigoplus_cU_c$ 满足
$\mathrm{Tr}(X_cU_c)=0$（引理 11.13）时，等号成立。$\square$

类型 (ii) 的检验可以利用直和：实验 3 的 $z$ 对直和是可加的，而直和的对比度是各直和项对比度中
最大者（习题 3），所以对于实验 3 的包络（$\Delta(Y)$ 在固定区间 $[0,b]$ 中自由变化），
$$
\overline z_{M_2\oplus M_2}(\delta)\ \ge\ \overline z_{M_2}(\delta)+\sup_{\delta'\le\delta}\overline z_{M_2}(\delta');
$$
特别地，在 $M_2$ 上于对比度 $\delta$ 与 $\delta'\le\delta$ 处取到的两个值之和，是在 $M_2\oplus M_2$ 上于对比度
$\delta$ 处取到的值。

*计算：* 用多起点搜索检验 $M_2\oplus M_2$ 上的命题 11.14、$M_3\oplus M_3$ 上的推论 11.7 以及命题 11.9，并做一项类型 (ii) 的检验——[教程，11.5 节](../../tutorial/zh/11_optimisation.md#s-11-5)。

<a id="s-11-6"></a>

## 11.6 可与之比较的工具

*库中：* 与其他软件的比较，其中包括 SDP 求解器，其凸松弛从上方界定上确界——[教程，11.6 节](../../tutorial/zh/11_optimisation.md#s-11-6)。

<a id="x-11"></a>

## 习题

1. （**Cayley 变换。**）对 $S^*=-S$，证明 $C(S)=(1-S)(1+S)^{-1}$ 是酉的，并且 $S\mapsto C(S)$ 是到
   没有特征值 $-1$ 的酉矩阵集合上的双射。因此，与指数映射（引理 11.1）不同，Cayley 坐标卡
   $U_0C(S)$ 取不到 $-U_0$。
2. （**实 Haar 测度。**）设 $O$ 在 $O(k)$ 上服从 Haar 分布（实 Ginibre 矩阵的 QR 分解），$A$ 为实矩阵，证明
   $\mathbb E\mathrm{Tr}(AO)^2=\|A\|_2^2/k$，并用实代数加以检验。
3. （**直和。**）证明实验 3 的 $z$ 对直和是可加的，并且对 $\lambda_{\max}=1$ 的 $A,A'\ge0$ 有
   $\Delta(A\oplus A')=\max(\Delta(A),\Delta(A'))$。由此推出 11.5 节中的直和
   不等式（其中 $\Delta(Y)$ 在固定区间 $[0,b]$ 中自由变化）。

<details>
<summary>解答</summary>

1. $S$ 是正规的，且 $iS$ 是自伴的，所以 $S$ 的特征值是纯虚数（推论 1.4(2)），从而 $1\pm S$ 可逆；
   它们可交换，且 $(1+S)^*=1-S$。因此
   $C^*C=(1-S)^{-1}(1+S)(1-S)(1+S)^{-1}=1$，而 $1+C=2(1+S)^{-1}$ 可逆，所以 $-1\notin\mathrm{spec}\,C$。
   反之，对满足 $1+U$ 可逆的酉矩阵 $U$，$S=(1-U)(1+U)^{-1}$ 是反自伴的：
   $S^*=(1+U^*)^{-1}(1-U^*)=(U+1)^{-1}U\,U^*(U-1)=(U+1)^{-1}(U-1)=-S$，因为 $U\pm1$ 可交换。对可交换的 $C$ 与
   $S$，$C(1+S)=1-S$ 等价于 $S(1+C)=1-C$（二者都是说 $C+CS+S=1$）。而 $C(S)$ 与 $S$ 可交换，
   $S=(1-U)(1+U)^{-1}$ 与 $U$ 可交换；所以 $U=C(S)$ 当且仅当 $S=(1-U)(1+U)^{-1}$：这两个映射互为逆映射。

2. 由定理 B.17（另见引理 C.2），实高斯矩阵经符号修正的 QR 因子 $O$——本库对实代数正是这样抽样的——
   在 $O(k)$ 上服从 Haar 分布，所以对 $V\in O(k)$，$VO$ 与 $O$ 同分布。
   左乘 $\mathrm{diag}(\pm1)$ 给出当 $i\ne l$ 时 $\mathbb E[O_{ij}O_{lm}]=0$，置换表明
   $\mathbb E[O_{ij}O_{im}]$ 与 $i$ 无关，而 $\sum_iO_{ij}O_{im}=\delta_{jm}$；所以
   $\mathbb E[O_{ij}O_{lm}]=\delta_{il}\delta_{jm}/k$，$\mathbb E\mathrm{Tr}(AO)^2=\sum_{ij}A_{ji}^2/k$。

3. 块对角算子的迹与迹范数都是各块之和，且 $|X\oplus X'|=|X|\oplus|X'|$；
   $\lambda_{\min}(A\oplus A')=\min(\lambda_{\min}(A),\lambda_{\min}(A'))$，而 $\Delta$ 随 $\lambda_{\min}$ 递减。
   所以若 $(X,Y)$ 在 $\delta$ 处容许，$(X',Y')$ 在 $\delta'\le\delta$ 处容许，则 $(X\oplus X',Y\oplus Y')$ 在
   $\delta$ 处容许，且 $z=z(X,Y)+z(X',Y')$；对两个数对以及 $\delta'\le\delta$ 取上确界即可。

</details>

*计算：* 习题 1 与 2 的数值验证——[教程，习题](../../tutorial/zh/11_optimisation.md#x-11)。

**下一章：** [第 12 章](12_equivariant_learning.md)讨论等变机器学习：与对称性可交换的层是
交换子（第 2 章）的元素，它们由同样的批量化、可微分的算子构造而成。
