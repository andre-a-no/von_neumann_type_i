[English](../00_overview.md) · [Español](../es/00_overview.md) · [Français](../fr/00_overview.md) · [Deutsch](../de/00_overview.md) · **中文** · [日本語](../ja/00_overview.md) · [Русский](../ru/00_overview.md)

# 0. 本课程讲什么

冯·诺依曼代数大多是用纸和笔来研究的。在有限维情形，每个冯·诺依曼代数都 $*$-同构于全矩阵代数的
直和（定理 2.6），因此原则上一切都可以计算——但当这样的代数被数值处理时，它常常被展平成一个大矩阵，
它的结构（扇区、中心、众多的迹）也就在计算中消失了。本课程和库 `torch_vn_algebra` 采取相反的做法：
**代数本身就是数据类型**，并在其中高效地进行计算。

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

在这种对应下，代数的运算就是批量数组运算。这种对应是精确的（表中前六行的每一行都是一个恒等式，
下面由课程中的结果加以证明）；只有执行这些运算的数值例程是在浮点运算中工作的：

| 在代数 $\mathcal M$ 中 | 在形状为 $(B,C,k,k)$ 的数组上 |
|---|---|
| 乘积 $AB$、和、伴随 | 在 $(B,C)$ 上批量进行的矩阵乘积、求和、共轭转置 |
| 函数演算 $f(A)$（第 1 章） | 各块的批量特征分解 |
| 迹 $\sum_cw_c\operatorname{Tr}A_c$（第 3 章） | 每块求对角和，再沿扇区轴加权求和 |
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
(2) 设 $A\in\mathcal M$ 是正规的，$\varphi(A)=\bigoplus_cA_c$，$f$ 是 $\operatorname{spec}A$ 上的函数。每个 $A_c$ 都是
正规的，因为 $\varphi(A)^*\varphi(A)=\varphi(A^*A)=\varphi(AA^*)=\varphi(A)\varphi(A)^*$ 逐块成立。若
$A_cx=\lambda x$ 且 $x\neq0$，则 $V(0\oplus\dots\oplus x\otimes e_1\oplus\dots\oplus0)$（$x\otimes e_1$ 位于第
$c$ 个直和项）是 $A$ 对应于 $\lambda$ 的特征向量，故 $\operatorname{spec}A_c\subseteq\operatorname{spec}A$。设 $p$ 是一个
在 $\operatorname{spec}A$ 上满足 $p=f$ 的多项式（拉格朗日插值，与定理 1.5 的证明中相同）。对任意满足
$\operatorname{spec}X\subseteq\operatorname{spec}A$ 且 $X=U\Lambda U^*$ 的正规 $X$，由定理 1.5(b) 有 $f(X)=Uf(\Lambda)U^*=Up(\Lambda)U^*=p(X)$；
因此 $f(A)=p(A)$，$f(A_c)=p(A_c)$，又因为 $\varphi$ 是保单位的同态，
$$
\varphi(f(A))=p(\varphi(A))=\bigoplus_cp(A_c)=\bigoplus_cf(A_c),\qquad f(A_c)=U_cf(\Lambda_c)U_c^*
$$
对任意对角化 $A_c=U_c\Lambda_cU_c^*$（定理 1.3）成立：每块一次特征分解。
(3) $\bigoplus_cM_{k_c}$ 上的迹恰为泛函 $\sum_cw_c\operatorname{Tr}A_c$（定理 3.1）。
(4) $Z(\mathcal M)$ 由极小中心投影 $P_c$ 张成，而 $\varphi(P_c)$ 在块 $c$ 中为 $1_{k_c}$，在其他块中为 $0$
（定理 2.6），故 $\varphi(Z(\mathcal M))=\{\bigoplus_cz_c1_{k_c}:z\in\mathbb C^C\}$。
(5) 命题 2.8。(6) 推论 5.4（各块相互独立，每块服从 $U(k_c)$ 上的 Haar 测度）和定理 5.3
（每块由一个高斯矩阵的 QR 分解得到）。最后一行只是一个指引：态由密度表示（定理 3.2），扇区之间的信道
由块之间的 Kraus 算子表示（命题 6.4），封闭动力学由 $e^{-itH}=\bigoplus_ce^{-itH_c}$ 表示（命题 8.2）。

有两个推论使这不仅仅是记账上的便利。

- **代价。** $\mathcal M$ 中的谱分解是逐块计算的（表中第 2 行）。标准的稠密特征值求解器（PyTorch 在
  CPU 上调用的 LAPACK）所用的浮点运算次数按矩阵大小的三次方增长；这是对这些算法的描述，而不是
  本课程的结果。在“$k\times k$ 特征值问题的代价为常数乘以 $k^3$”的模型中，$\mathcal M$ 中的一次谱分解
  代价为 $\sum_ck_c^3$，而在嵌入 $M_D(\mathbb C)$（$D=\sum_ck_c$）中则为 $D^3$；对于 $C$ 个大小为 $k$ 的
  相等的块，二者之比为 $(Ck)^3/(Ck^3)=C^2$。乘积在补零后的张量上进行；两个 $k\times k$ 矩阵的
  教科书式乘法需要 $k^3$ 次乘法（$k^2$ 个元素，每个是 $k$ 个乘积之和），因此 $\mathcal M$ 中的一次
  乘积需要 $C\,k_{\max}^3$ 次，当所有块大小相同时它等于 $\sum_ck_c^3$。一个元素以 $C\,k_{\max}^2$ 个数
  存储，而嵌入中需要 $D^2$ 个；对于相等的块，二者之比为 $(Ck)^2/(Ck^2)=C$。代码实际的做法是：乘积、
  求和与迹是在整个 $(B,C,k_{\max},k_{\max})$ 数组上的单次张量运算；埃尔米特特征分解、奇异值和 QR
  分解则对每种不同块大小 $k$ 的全部 $B\cdot n_k$ 个块（$n_k$ 为该大小的扇区数）组成的堆栈做一次
  批量调用——当所有扇区大小相同时就只有一次调用。在有 $t\ge2$ 个 PyTorch 线程的 CPU 上，只要堆栈
  至少有 $2t$ 个块且（块数）$\,\cdot k^3\ge2^{18}$（模块 `cost`，函数 `_run`），库就会把这样的调用拆成
  $t$ 块，在并行线程中运行；结果是对相同矩阵的相同 LAPACK 调用。其目的是让 GPU 或 CPU 的所有核心
  保持满载。
- **结构上保证正确。** 每个形状为 $(C,k_{\max},k_{\max})$ 且补零的数组都是 $\mathcal M$ 的元素，因此
  由这些块计算出的一切都留在 $\mathcal M$ 中。库中的随机酉元是用定理 5.3 的构造在每个块中独立
  抽取的（`random_unitary_operator` 对每种块大小调用一次逐块采样器 `random_unitary`，第 5.1 节
  对照这一构造检验了 `random_unitary`），由推论 5.4，这就是 $\prod_cU(k_c)$ 的 Haar 测度（对于理想的
  高斯随机数）；当 $C\ge2$ 时，$U(D)$ 的 Haar 酉元以概率零落在这个子群中（第 5.1 节，推论 5.4 之后）。
  梯度步改变的是各块，因此不会离开代数。而迹是定理 3.1 中带有显式选定权重的泛函
  $\sum_cw_c\operatorname{Tr}A_c$（密度和条件期望也是相对于它们来取的），而不是 $\mathcal M$ 恰好作用于其上的
  某个 $\mathbb C^n$ 的迹——那将固定 $w_c=m_c$（第 2.3 节）。

## 0.2 第一个计算

我们把同一个计算分别在 $\mathcal M=\bigoplus_{c=1}^{64}M_{16}(\mathbb C)$ 及其在 $M_{1024}(\mathbb C)$
中的嵌入里完成：求一个随机正元的平方根及其迹。`TypeIAlgebra(n_factors, k_factors)` 接受两个块大小
列表：名义大小 $n_c$ 和实际使用的块的大小 $k_c\le n_c$（第 1.5 节）。本课程用 `TypeIAlgebra` 构造的
每个代数，其两个列表都相同。

```python
import time
import torch
from torch_vn_algebra import TypeIAlgebra

torch.manual_seed(0)
C, k = 64, 16
alg = TypeIAlgebra([k] * C, [k] * C, precision='double')       # M = (+)_{c=1}^{64} M_16(C)

A = alg.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(1, d), force_positive=True)
blocks = A.matrix                                                # build A before timing

t0 = time.perf_counter()
root_blocks = A.sqrt().matrix                                    # (1, 64, 16, 16): 64 small eigenproblems
t_blocks = time.perf_counter() - t0

big = torch.block_diag(*blocks[0])                               # the same element in M_1024(C)
t0 = time.perf_counter()
w, V = torch.linalg.eigh(big)                                    # one large eigenproblem
root_big = (V * w.sqrt()) @ V.conj().T
t_big = time.perf_counter() - t0

assert torch.allclose(torch.block_diag(*root_blocks[0]), root_big, atol=1e-10)   # row 2: f(A) = (+)_c f(A_c)
assert abs(A.sqrt().Tr_blunt().item() - root_big.trace().item()) < 1e-9
assert big.numel() == C * blocks.numel()                          # storage: a factor C = 64
print(f"in the algebra: {t_blocks * 1e3:.1f} ms, as one 1024 x 1024 matrix: {t_big * 1e3:.1f} ms")
```

两者的结果在舍入误差范围内一致（表中第 2 行），存储之比为 $C$，与上面的计数相符。64 个块作为对
形状为 $(64,16,16)$ 的堆栈的一次批量调用交给特征值求解器；在 $2\le t\le32$ 个 PyTorch 线程下，库把它
拆成 $t$ 块（$64\cdot16^3=2^{18}$ 恰好达到模块 `cost` 的阈值），单线程时则不拆分。*数值观察*（计时
取决于硬件；代码会把它们打印出来）：在两台 4 核 CPU 上，代数版本耗时 3–10 ms，大矩阵版本耗时
400–760 ms，倍数在几十到几百之间（在第二台机器上，4 线程和单线程都是如此）。测得的倍数远低于三次
代价模型给出的比值 $C^2=4096$：运算次数是运行时间的模型，而不是运行时间本身；在这些机器上，小型
特征值问题每秒执行的运算次数远低于大型问题。在代价模型中，比值 $C^2$ 随扇区数增长。对于 $B$ 个
元素的批量，代数版本仍然只需一次批量调用；大矩阵版本则必须求解 $B$ 个大小为 $1024$ 的特征值问题。

## 0.3 已有的工作与新的内容

我们查找过在有限维冯·诺依曼代数*内部*进行数值计算的软件。找到的软件可分为四类，它们都很有用，
但目的各不相同：

- **通用矩阵软件**（NumPy/SciPy、QuTiP 以及其他量子工具箱）一次处理一个矩阵。即使存在块结构，
  也必须手工维护。
- **对称张量网络库**把张量存储为按某种对称性的荷扇区划分的块稀疏数组（在 YASTN [1]、TeNPy 和
  ITensor 中是阿贝尔对称性；TensorKit 也处理非阿贝尔对称性）；YASTN 带有支持 GPU 和自动微分的
  PyTorch 后端。从技术上看，它们是最接近的亲缘。它们处理的对象是待收缩的张量网络，而不是算子代数：
  就我们在其文档中所见，它们不提供迹、态和密度的族，不提供中心和条件期望、$\prod_cU(k_c)$ 上的 Haar
  测度，也不提供扇区之间的信道。
- 矩阵 $*$-代数的**数值块对角化** [2, 3, 4] 求出给定代数的分解
  $\mathcal M\cong\bigoplus_cM_{k_c}\otimes1_{m_c}$——通常是为了缩小具有对称性的半定规划的规模。
  这与本库是互补的：那些算法找出块，本库在块中计算。本库的输入步骤（模块 `decompose`，第 4.3 节）
  采用了 [2, 3] 的思想，即交换子中一个随机元素的特征空间能找出这些块；第 4.3 节证明了在精确算术下
  这以概率一成功（引理 4.6、命题 4.7、引理 4.8），而命题 4.9 是 `Decomposition.verify()` 数值检验的
  判据。群代数的精确符号分解在计算机代数中已有实现（例如 GAP 程序包 Wedderga）。
- **非交换多项式优化**（NCSOStools、ncpol2sdpa）通过半定松弛 [5]，对所有矩阵规模同时给出
  非交换变量多项式的特征值和迹的界。它以不同的手段（凸松弛）回答一个不同的问题（对所有矩阵规模都
  成立的界）。

这一综述是对我们所找到的软件的描述，而不是定理。据我们所知，还没有哪个通用库把有限维冯·诺依曼代数
作为基本数据类型，用于批量的、可微的数值计算。因此，新意不在于新的数学，而在于**高效的计算**：把
§0.1 的对应关系实现出来，使代数的结构得以保持，并使运算作为批量内核运行（对于分解，每种不同的块
大小一次调用，在多核 CPU 上超过 §0.1 的规模阈值时再拆成每个线程一块）。如果您知道相关的软件，我们
很乐意了解。

## 0.4 如何阅读本课程

| 章 | 内容 |
|---|---|
| 1–3 | 数学部分：矩阵代数、结构定理、迹、态与熵——附完整证明 |
| 4–8 | 这些代数从何而来以及在其中计算什么：对称性、随机矩阵、信道、张量积、动力学 |
| 9–12 | 应用：自旋链、量子计算、优化、等变机器学习 |
| 13 | 基础：若尔当代数，以及量子力学为何使用 *-代数 |
| 附录 | 各章按编号引用的辅助结果的完整证明（例如定理 A.24） |

第 1–13 章的每一章都给出结果及其证明，展示库代码，用课程中证明的结果对代码进行数值检验，并指出
可供比较的标准工具。预备知识（线性代数、实分析、测度论、初等概率论）列在
[README](README.md#预备知识) 中；其余一切都在课程中证明，其他章节和附录中的结果按编号引用。课程
不加证明的数值（例如计时）只作为该章自身代码的*数值观察*出现，而标题为*物理意义*的方框是物理解释，
而非本课程的定理。

## 参考文献

1. M. M. Rams, G. Wójtowicz, A. Sinha, J. Hasik, *YASTN: Yet another symmetric tensor networks; A
   Python library for Abelian symmetric tensor network calculations*, SciPost Phys. Codebases 52
   (2025).
2. K. Murota, Y. Kanno, M. Kojima, S. Kojima, *A numerical algorithm for block-diagonal
   decomposition of matrix $*$-algebras with application to semidefinite programming*, Japan J.
   Indust. Appl. Math. 27 (2010) 125–160（预印本：Technical Report METR 2007-52, University of
   Tokyo, 2007）.
3. T. Maehara, K. Murota, *A numerical algorithm for block-diagonal decomposition of matrix
   $*$-algebras with general irreducible components*, Japan J. Indust. Appl. Math. 27 (2010)
   263–293.
4. E. de Klerk, C. Dobre, D. V. Pasechnik, *Numerical block diagonalization of matrix
   $*$-algebras with application to semidefinite programming*, Math. Program. 129 (2011) 91–111.
5. I. Klep, J. Povh, *Semidefinite programming and sums of hermitian squares of noncommutative
   polynomials*, J. Pure Appl. Algebra 214 (2010) 740–749.

**下一章：** [第 1 章](01_matrix_algebras.md)。
