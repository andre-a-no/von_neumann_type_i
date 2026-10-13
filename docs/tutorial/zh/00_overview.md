[English](../00_overview.md) · [Español](../es/00_overview.md) · [Français](../fr/00_overview.md) · [Deutsch](../de/00_overview.md) · **中文** · [日本語](../ja/00_overview.md) · [Русский](../ru/00_overview.md)

# 0. 本课程讲什么

冯·诺依曼代数几乎完全是用纸和笔来研究的。在有限维情形，每个冯·诺依曼代数都是全矩阵代数的
直和（第 2 章），因此原则上一切都可以计算——但实际上，代数通常被展平成一个大矩阵，它的结构
（扇区、中心、众多的迹）也就在计算中消失了。本课程和库 `torch_vn_algebra` 采取相反的做法：
**代数本身就是数据类型**，并在其中高效地进行计算。

## 0.1 基本观察

有限维冯·诺依曼代数在相差一个酉基变换的意义下为
$$
\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c}
$$
（定理 2.6）。作为抽象代数，它是 $\bigoplus_cM_{k_c}(\mathbb C)$：一个元素就是由 $C$ 个方块组成的
列表。这恰好是带有*扇区轴*的张量的形状：$\mathcal M$ 的一个元素是形状为 $(C,k,k)$ 的数组（各块
填充到 $k=\max_ck_c$），而 $B$ 个元素——蒙特卡罗样本、优化的初始点、时间步、无序实现——构成
形状为 $(B,C,k,k)$ 的数组。现代数组库（CPU 和 GPU 上的 PyTorch）正是为这种“小矩阵的批量”而
设计的，因为神经网络需要它们。

在这种对应下，代数的运算就是批量数组运算，而且是精确的，没有任何近似：

| 在代数 $\mathcal M$ 中 | 在形状为 $(B,C,k,k)$ 的数组上 |
|---|---|
| 乘积 $AB$、和、伴随 | 在 $(B,C)$ 上批量进行的矩阵乘积、求和、共轭转置 |
| 函数演算 $f(A)$（第 1 章） | 各块的批量特征分解 |
| 迹 $\sum_cw_c\operatorname{Tr}A_c$（第 3 章） | 每块求对角和，再沿扇区轴加权求和 |
| 中心、中心投影（第 2 章） | 每个扇区一个数 |
| 到中心上的条件期望 | 块的迹除以 $k_c$，再广播回去 |
| 酉群 $\prod_cU(k_c)$ 上的 Haar 测度 | 高斯随机块的批量 QR 分解（第 5 章） |
| 态、信道、动力学（第 3、6、8 章） | 同样的运算，加以复合 |

有两个推论使这不仅仅是记账上的便利。

- **代价。** 在 $\mathcal M$ 中做一次谱分解需要 $\sum_ck_c^3$ 次运算（各块按大小分组），而在嵌入
  $M_D(\mathbb C)$（$D=\sum_ck_c$）中则需要 $D^3$ 次；对于 $C$ 个相等的块，二者相差 $C^2$ 倍。
  乘积与求和在补零后的张量上进行，代价为 $C\,k_{\max}^3$，当所有块大小相同时它等于 $\sum_ck_c^3$。
  乘积、求和与迹由一次调用处理全部 $B\cdot C$ 个块；分解（特征值、奇异值）则对每种不同的块大小各做
  一次批量调用——当所有扇区大小相同时就只需一次调用——而在多核 CPU 上，库会把每次这样的调用
  拆成每个线程一块。这能让 GPU 或 CPU 的所有核心保持满载。
- **结构上保证正确。** 计算出的一切都留在 $\mathcal M$ 中：随机酉元在 $\prod_cU(k_c)$ 上服从
  Haar 分布（这是 $U(D)$ 中的零测集，对 $U(D)$ 采样永远得不到它），梯度步不会离开代数，而迹、
  密度和条件期望都是相对于代数本身的，而不是相对于某个偶然的嵌入。

## 0.2 第一个计算

我们把同一个计算分别在 $\mathcal M=\bigoplus_{c=1}^{64}M_{16}(\mathbb C)$ 及其在 $M_{1024}(\mathbb C)$
中的嵌入里完成：求一个随机正元的平方根及其迹。`TypeIAlgebra(n_factors, k_factors)` 接受两个块大小
列表：名义大小 $n_c$ 和实际使用的块的大小 $k_c\le n_c$；在本课程中两个列表总是相同的（第 1.5 节）。

```python
import time
import torch
from torch_vn_algebra import TypeIAlgebra

torch.manual_seed(0)
C, k = 64, 16
alg = TypeIAlgebra([k] * C, [k] * C, precision='double')       # M = (+)_{c=1}^{64} M_16(C)

A = alg.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(1, d), force_positive=True)

t0 = time.perf_counter()
root_blocks = A.sqrt().matrix                                    # (1, 64, 16, 16): 64 small eigenproblems
t_blocks = time.perf_counter() - t0

big = torch.block_diag(*A.matrix[0])                             # the same element in M_1024(C)
t0 = time.perf_counter()
w, V = torch.linalg.eigh(big)                                    # one large eigenproblem
root_big = (V * w.sqrt()) @ V.conj().T
t_big = time.perf_counter() - t0

assert torch.allclose(torch.block_diag(*root_blocks[0]), root_big, atol=1e-10)   # same result
assert abs(A.sqrt().Tr_blunt().item() - root_big.trace().item()) < 1e-9
print(f"in the algebra: {t_blocks * 1e3:.1f} ms, as one 1024 x 1024 matrix: {t_big * 1e3:.1f} ms")
```

两者的结果在舍入误差范围内一致。在 CPU 上代数版本快数十倍（在一台 4 核机器上约为 10 ms 对
400 ms；确切数字取决于硬件——代码会把它们打印出来）。64 个块作为一次批量调用交给特征值求解器，
在多核 CPU 上库会把这次调用拆成每个线程一块。但加速仍远不及运算次数之比 $C^2=4096$ 所暗示的那么多，
因为一次大型特征值求解比许多次小型求解更能充分利用硬件。扇区数越多，差距越大，而内存节省为 $C$ 倍。
对于 $B$ 个元素的批量，代数版本仍然只需一次批量调用；大矩阵版本则变成 $B$ 个大型特征值问题。

## 0.3 已有的工作与新的内容

我们查找过在有限维冯·诺依曼代数*内部*进行数值计算的软件。找到的软件可分为四类，它们都很有用，
但目的各不相同：

- **通用矩阵软件**（NumPy/SciPy、QuTiP 以及其他量子工具箱）一次处理一个矩阵。即使存在块结构，
  也必须手工维护。
- **对称张量网络库**把张量存储为按阿贝尔对称性的荷扇区划分的块稀疏数组，其中一些带有支持 GPU
  和自动微分的 PyTorch 后端（例如 YASTN [1]；还有 TeNPy、ITensor、TensorKit）。从技术上看，它们
  是最接近的亲缘。它们处理的对象是待收缩的张量网络，而不是算子代数：它们不提供迹、态和密度的族，
  不提供中心和条件期望、$\prod_cU(k_c)$ 上的 Haar 测度，也不提供扇区之间的信道。
- 矩阵 $*$-代数的**数值块对角化** [2, 3, 4] 求出给定代数的分解
  $\mathcal M\cong\bigoplus_cM_{k_c}\otimes1_{m_c}$——通常是为了缩小具有对称性的半定规划的规模。
  这与本库是互补的：那些算法找出块，本库在块中计算，并实现了其中一种算法作为输入步骤（`decompose`，第 4 章）。群代数的精确符号分解在计算机代数中已有实现
  （例如 GAP 程序包 Wedderga）。
- **非交换多项式优化**（NCSOStools、ncpol2sdpa）通过半定松弛 [5]，对所有矩阵规模同时给出
  非交换变量多项式的特征值和迹的界。它以不同的手段（凸松弛）回答一个不同的问题（对所有代数都
  成立的界）。

据我们所知，还没有哪个通用库把有限维冯·诺依曼代数作为基本数据类型，用于批量的、可微的数值计算。
因此，新意不在于新的数学，而在于**高效的计算**：把 §0.1 的对应关系实现出来，使代数的结构得以
保持，并使每个运算都作为批量内核运行（对于分解，每种不同的块大小一个内核，在多核 CPU 上再拆成每个线程一块）。如果您知道相关的软件，我们很乐意了解。

## 0.4 如何阅读本课程

| 章 | 内容 |
|---|---|
| 1–3 | 数学部分：矩阵代数、结构定理、迹、态与熵——附完整证明 |
| 4–8 | 这些代数从何而来以及在其中计算什么：对称性、随机矩阵、信道、张量积、动力学 |
| 9–12 | 应用：自旋链、量子计算、优化、等变机器学习 |

每一章都给出结果及其证明，展示库代码，用精确结果检验代码，并指出可供比较的标准工具。

## 参考文献

1. M. M. Rams, G. Wójtowicz, A. Sinha, J. Hasik, *YASTN: Yet another symmetric tensor networks; A
   Python library for Abelian symmetric tensor network calculations*, SciPost Phys. Codebases 52
   (2025).
2. K. Murota, Y. Kanno, M. Kojima, S. Kojima, *A numerical algorithm for block-diagonal
   decomposition of matrix $*$-algebras*, Technical Report METR 2007-52, University of Tokyo (2007).
3. T. Maehara, K. Murota, *A numerical algorithm for block-diagonal decomposition of matrix
   $*$-algebras, Part II: general algorithm*, Optimization Online (2009).
4. E. de Klerk, C. Dobre, D. V. Pasechnik, *Numerical block diagonalization of matrix
   $*$-algebras with application to semidefinite programming*, Math. Program. 129 (2011) 91–111.
5. I. Klep, J. Povh, *Semidefinite programming and sums of hermitian squares of noncommutative
   polynomials*, J. Pure Appl. Algebra 214 (2010) 740–749.

**下一章：** [第 1 章](01_matrix_algebras.md)。
