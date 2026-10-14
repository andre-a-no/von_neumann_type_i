[English](../00_overview.md) · [Español](../es/00_overview.md) · [Français](../fr/00_overview.md) · [Deutsch](../de/00_overview.md) · **中文** · [日本語](../ja/00_overview.md) · [Русский](../ru/00_overview.md)

# 0. 本课程讲什么

冯·诺依曼代数大多是用纸和笔来研究的。在有限维情形，每个冯·诺依曼代数都 $*$-同构于全矩阵代数的
直和（定理 2.6），因此原则上一切都可以计算——但当这样的代数被数值处理时，它常常被展平成一个大矩阵，
它的结构（扇区、中心、众多的迹）也就在计算中消失了。本课程和库 `torch_vn_algebra` 采取相反的做法：
**代数本身就是数据类型**，并在其中高效地进行计算。

<a id="s-0-1"></a>

## 0.1 基本观察

$\mathbb C^n$ 上的冯·诺依曼代数 $\mathcal M$ 在相差一个酉基变换的意义下为
$$
\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c}
$$
（定理 2.6）。作为抽象代数，它是 $\bigoplus_cM_{k_c}(\mathbb C)$：一个元素就是由 $C$ 个方块组成的
列表。这恰好是带有*扇区轴*的张量的形状：$\mathcal M$ 的一个元素是形状为 $(C,k,k)$ 的数组（各块
用零填充到 $k=\max_ck_c$；库正是这样存储它的，见第 2.5 节），而 $B$ 个元素——蒙特卡罗样本、优化的
初始点、时间步、无序实现——构成形状为 $(B,C,k,k)$ 的数组。现代数组库（CPU 和 GPU 上的 PyTorch）
正是为这种“小矩阵的批量”提供了批量例程——矩阵乘积、特征分解、QR 分解——神经网络也使用这些例程。

*库中：* $\mathcal M$ 的元素是 `TypeIAlgebra` 的算子，存储为一个形状为 $(B,C,k_{\max},k_{\max})$ 的张量——[教程，1.5 节](../../tutorial/zh/01_matrix_algebras.md#s-1-5) 和 [教程，2.5 节](../../tutorial/zh/02_von_neumann_algebras.md#s-2-5)。

在这种对应下，代数的运算就是批量数组运算。这种对应是精确的（表中前六行的每一行都是一个恒等式，
下面由课程中的结果加以证明）；只有执行这些运算的数值例程是在浮点运算中工作的：

| 在代数 $\mathcal M$ 中 | 在形状为 $(B,C,k,k)$ 的数组上 |
|---|---|
| 乘积 $AB$、和、伴随 | 在 $(B,C)$ 上批量进行的矩阵乘积、求和、共轭转置 |
| 函数演算 $f(A)$（第 1 章） | 各块的批量特征分解 |
| 迹 $\sum_cw_c\mathrm{Tr}\,A_c$（第 3 章） | 每块求对角和，再沿扇区轴加权求和 |
| 中心、中心投影（第 2 章） | 每个扇区一个数 |
| 到中心上的条件期望 | 块的迹除以 $k_c$，再广播回去 |
| 酉群 $\prod_cU(k_c)$ 上的 Haar 测度 | 高斯随机块的批量 QR 分解，加上定理 5.3 的相位校正 |
| 态、信道、动力学（第 3、6、8 章） | 由上述运算构建 |

*前六行的证明。* 由定理 2.6，存在酉算子 $V$ 使得
$\mathcal M=V\big(\bigoplus_cM_{k_c}\otimes1_{m_c}\big)V^*$，并且
$\varphi\big(V(\bigoplus_cA_c\otimes1_{m_c})V^*\big)=\bigoplus_cA_c$ 是 $\mathcal M$ 到
$\bigoplus_cM_{k_c}(\mathbb C)$ 上的 $*$-同构（它是良定义的，因为 $A_c\otimes1_{m_c}$ 决定 $A_c$；它是乘性的，因为
$(A_c\otimes1)(B_c\otimes1)=A_cB_c\otimes1$）。库存储的是 $\varphi(A)$（第 2.5 节）。
(1) 在 $\bigoplus_cM_{k_c}$ 中，和、乘积与伴随都是逐块计算的。
(2) 设 $A\in\mathcal M$ 是正规的，$\varphi(A)=\bigoplus_cA_c$，$f$ 是 $\mathrm{spec}\,A$ 上的函数。每个 $A_c$ 都是
正规的，因为 $\varphi(A)^*\varphi(A)=\varphi(A^*A)=\varphi(AA^*)=\varphi(A)\varphi(A)^*$ 逐块成立。若
$A_cx=\lambda x$ 且 $x\neq0$，则 $V(0\oplus\dots\oplus x\otimes e_1\oplus\dots\oplus0)$（$x\otimes e_1$ 位于第
$c$ 个直和项）是 $A$ 对应于 $\lambda$ 的特征向量，故 $\mathrm{spec}\,A_c\subseteq\mathrm{spec}\,A$。设 $p$ 是一个
在 $\mathrm{spec}\,A$ 上满足 $p=f$ 的多项式（拉格朗日插值，与定理 1.5 的证明中相同）。对任意满足
$\mathrm{spec}\,X\subseteq\mathrm{spec}\,A$ 且 $X=U\Lambda U^*$ 的正规 $X$，由定理 1.5(b) 有 $f(X)=Uf(\Lambda)U^*=Up(\Lambda)U^*=p(X)$；
因此 $f(A)=p(A)$，$f(A_c)=p(A_c)$，又因为 $\varphi$ 是保单位的同态，
$$
\varphi(f(A))=p(\varphi(A))=\bigoplus_cp(A_c)=\bigoplus_cf(A_c),\qquad f(A_c)=U_cf(\Lambda_c)U_c^*
$$
对任意对角化 $A_c=U_c\Lambda_cU_c^*$（定理 1.3）成立：每块一次特征分解。
(3) $\bigoplus_cM_{k_c}$ 上的迹恰为泛函 $\sum_cw_c\mathrm{Tr}\,A_c$（定理 3.1）。
(4) $Z(\mathcal M)$ 由极小中心投影 $P_c$ 张成，而 $\varphi(P_c)$ 在块 $c$ 中为 $1_{k_c}$，在其他块中为 $0$
（定理 2.6），故 $\varphi(Z(\mathcal M))=\{\bigoplus_cz_c1_{k_c}:z\in\mathbb C^C\}$。
(5) 命题 2.8。(6) 推论 5.4（各块相互独立，每块服从 $U(k_c)$ 上的 Haar 测度）和定理 5.3
（每块由一个高斯矩阵的 QR 分解得到）。最后一行只是一个指引：态由密度表示（定理 3.2），扇区之间的信道
由块之间的 Kraus 算子表示（命题 6.4），封闭动力学由 $e^{-itH}=\bigoplus_ce^{-itH_c}$ 表示（命题 8.2）。

*库中：* 乘积以及 `sqrt`、`apply_function` 等函数——[教程，1.5 节](../../tutorial/zh/01_matrix_algebras.md#s-1-5)；迹，例如 `Tr_blunt`——[教程，3.1 节](../../tutorial/zh/03_traces_and_states.md#s-3-1)；中心，`TypeIAlgebra.central`——[教程，2.5 节](../../tutorial/zh/02_von_neumann_algebras.md#s-2-5)；到中心上的条件期望，`channels.center_expectation`——[教程，2.6 节](../../tutorial/zh/02_von_neumann_algebras.md#s-2-6)；Haar 酉元，`random_unitary_operator`——[教程，5.1 节](../../tutorial/zh/05_random_matrices.md#s-5-1)。

*计算：* 表中第 2 行在 $\bigoplus_{c=1}^{64}M_{16}(\mathbb C)$ 上的验证——[教程，0.2 节](../../tutorial/zh/00_overview.md#s-0-2)；逐块计算的代价——[教程，0.1 节](../../tutorial/zh/00_overview.md#s-0-1)。

<a id="s-0-2"></a>

## 0.2 第一个计算

*计算：* $\bigoplus_{c=1}^{64}M_{16}(\mathbb C)$ 中一个随机正元的平方根及其迹，分别逐块计算和在 $M_{1024}(\mathbb C)$ 中计算——[教程，0.2 节](../../tutorial/zh/00_overview.md#s-0-2)。

<a id="s-0-3"></a>

## 0.3 已有的工作与新的内容

*库中：* 相关的软件，以及本库在其基础上增加了什么——[教程，0.3 节](../../tutorial/zh/00_overview.md#s-0-3)。

<a id="s-0-4"></a>

## 0.4 如何阅读本课程

| 章 | 内容 |
|---|---|
| 1–3 | 数学部分：矩阵代数、结构定理、迹、态与熵——附完整证明 |
| 4–8 | 这些代数从何而来以及在其中计算什么：对称性、随机矩阵、信道、张量积、动力学 |
| 9–12 | 应用：自旋链、量子计算、优化、等变机器学习 |
| 13 | 基础：若尔当代数，以及量子力学为何使用 *-代数 |
| 附录 | 各章按编号引用的辅助结果的完整证明（例如定理 A.24） |

第 1–13 章的每一章都给出结果及其完整证明。预备知识（线性代数、实分析、测度论、初等概率论）列在
[README](../../lecture_notes/zh/README.md#prerequisites) 中；其余一切都在讲义中证明，其他章节和附录中的结果按编号引用。标题为
*物理意义*的方框是物理解释，而非讲义中的定理。

每一章在[教程](../../tutorial/zh/00_overview.md)中都有一个对应的章节，它展示库代码，用这里证明的结果对代码进行数值检验，并指出可供比较的标准工具。结果之后的 *计算：* 和 *库中：* 行链接到教程中计算或使用该结果的位置。

**下一章：** [第 1 章](01_matrix_algebras.md)。
