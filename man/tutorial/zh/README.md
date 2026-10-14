[English](../README.md) · [Español](../es/README.md) · [Français](../fr/README.md) · [Deutsch](../de/README.md) · **中文** · [日本語](../ja/README.md) · [Русский](../ru/README.md)

# 教程：在讲义的背景下使用 torch_vn_algebra

本教程在[讲义](../../lecture_notes/zh/README.md)的背景下，展示如何在库 `torch_vn_algebra` 中对有限维冯·诺依曼代数
$\bigoplus_cM_{k_c}(\mathbb C)$ 进行计算。每个教程配合讲义中编号相同的一章：它展示那里证明的结果对应的库代码，用精确
结果（定理、闭合公式或独立计算）检验代码，并在*可与之比较的工具*一节中列出完成同一任务的标准工具。每一节都以一行
*讲义：* 开头，链接到该节及其所用的结果。只是计算出来的数值（例如大 $n$ 时的能级统计）标记为*数值观察*，从不用作
前提或测试目标。除讲义的[预备知识](../../lecture_notes/zh/README.md#prerequisites)之外，本教程还假定读者熟悉 Python。

教程中的全部代码都由测试套件（`tests/test_man.py`）执行，因此每个例子都能在当前版本的
库中运行。

## 如何运行代码

```bash
pip install -e .                # from the root of the repository
python -m pytest tests/test_man.py -q
```

一章中的所有代码块可以按顺序粘贴到同一个 Python 会话或 Jupyter 笔记本中运行。同一章的
代码块共享变量。

## 目录

| | 章 | 主要结果 |
|---|---|---|
| 0 | [本课程讲什么](00_overview.md) | 代数作为数据类型：块结构 = 批量张量；代价；已有的工作与新的内容 |
| 1 | [矩阵代数与函数演算](01_matrix_algebras.md) | 谱定理、函数演算、极分解、迹范数对偶、GUE 的一个矩 |
| 2 | [有限维冯·诺依曼代数](02_von_neumann_algebras.md) | 双交换子定理、结构定理、中心、对称群与交换子、条件期望 |
| 3 | [迹、态与熵](03_traces_and_states.md) | 迹的分类、密度、扇区的熵、Klein 不等式、在到中心的条件期望下的单调性、Michelson 对比度 |
| 4 | [对称性与表示](04_symmetries.md) | Schur 引理、荷扇区、$U(1)$ 与 $SU(2)$、自旋系统、实型、复型与四元数型 |
| 5 | [代数中的随机矩阵](05_random_matrices.md) | $\prod_cU(k_c)$ 上的 Haar 测度、Weingarten 演算、能级统计以及为何不能混合扇区 |
| 6 | [量子信道](06_quantum_channels.md) | Kraus 表示与 Stinespring 表示、关于迹的对偶、扇区之间的映射、退相干 |
| 7 | [张量积与纠缠](07_tensor_products.md) | 偏迹、Page 公式、荷的融合 |
| 8 | [动力学](08_dynamics.md) | Schrödinger 方程与 Lindblad 方程、扇区间的跃迁、后选择 |
| 9 | [自旋链](09_spin_chains.md) | 磁化与动量扇区、Lanczos 方法、Bethe 拟设检验、无序 |
| 10 | [量子计算](10_quantum_computing.md) | 无退相干子空间、算子代数纠错、带约束的混合算子 |
| 11 | [优化与迹不等式](11_optimisation.md) | 带约束的参数化、采样与优化的对比 |
| 12 | [等变机器学习](12_equivariant_learning.md) | 作为交换子的等变层 |
| 13 | [若尔当代数：量子力学为何使用 *-代数](13_jordan_algebras.md) | 若尔当积与二次表示、自旋因子与洛伦兹锥、Jordan–von Neumann–Wigner 分类、转置与部分转置、若尔当形式的迹不等式 |
| A | [分析与概率中的若干事实](A_analysis.md) | 幂级数与零测集、高斯向量、球面测度、特征值不等式、Isserlis–Wick 公式 |
| B | [Schur–Weyl 对偶与 Weingarten 演算](B_schur_weyl.md) | Schur–Weyl 对偶、Weingarten 函数、Wigner–Eckart 定理、$O(k)$ 上的 Haar 测度 |
| C | [特征值密度与 Page 公式](C_eigenvalue_densities.md) | Weyl 积分公式、GOE/GUE、Wishart 密度与诱导密度、Page 公式 |
| D | [量子动力学半群的生成元](D_gksl.md) | GKSL 定理、直和上的生成元、Frigerio 定理、Bény–Kempf–Kribs |
| E | [欧几里得若尔当代数的分类](E_jordan_classification.md) | Peirce 分解、Hurwitz 定理与 Frobenius 定理、坐标化、Albert 代数、Jordan–von Neumann–Wigner 定理 |
| F | [对称锥、序同构与若尔当映射](F_cones_and_order.md) | 平方锥、Koecher–Vinberg 定理、Kadison 定理、$\bigoplus_cM_{k_c}$ 的若尔当自同构 |

数学部分——定义、定理和完整证明——见[讲义](../../lecture_notes/zh/README.md)，其中还有预备知识、延伸阅读和本库的配套论文。
