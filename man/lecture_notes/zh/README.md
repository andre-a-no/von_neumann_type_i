[English](../README.md) · [Español](../es/README.md) · [Français](../fr/README.md) · [Deutsch](../de/README.md) · **中文** · [日本語](../ja/README.md) · [Русский](../ru/README.md)

# 有限维冯·诺依曼代数：讲义

本课程讲授有限维冯·诺依曼代数（von Neumann 代数）——即分块对角矩阵构成的代数
$`\bigoplus_cM_{k_c}(\mathbb C)`$——的理论，以及它们在量子物理、量子信息、随机矩阵理论和
机器学习中的应用：对称性扇区与超选择扇区、双交换子定理与 Wedderburn 结构定理、迹、态与熵、Haar 测度与 Weingarten
演算、量子信道（Kraus、Stinespring、Choi）、Lindblad（GKSL）方程、自旋链、量子纠错、等变神经网络以及若尔当代数。

本课程同时面向两类读者：

- **数学家**将看到一个自成体系的文本：文中陈述的每个定理、公式、常数和不等式，都在课程中
  从下面列出的预备知识出发得到证明——或在该章之内，或在按编号引用的附录 A–F 之一中。只是
  计算出来的数值（例如大 $`n`$ 时的能级统计）标记为*数值观察*，从不用作前提或测试目标。参考
  文献只用于历史溯源或延伸阅读，从不代替证明。
- **物理学家和工程师**可以在标有*物理意义*的方框中读到每个概念的物理含义。这些方框是物理解释，而非本课程的定理：其他地方证明的任何结论都不
  依赖于它们。可以直接改写使用的可运行代码见[教程](../../tutorial/zh/README.md)。

第 1–13 章的结构相同：带完整证明的定义与结果、算出的例子，以及附有解答的习题。在库 `torch_vn_algebra`
计算或检验的每个结果之后，都有一行 *计算：* 或 *库中：*，链接到[教程](../../tutorial/zh/README.md)中完成这一工作的位置。

## 目录

| | 章 | 主要结果 |
|---|---|---|
| 0 | [本课程讲什么](00_overview.md) | 代数作为数据类型：块结构 = 批量张量；代价；已有的工作与新的内容 |
| 1 | [矩阵代数与函数演算](01_matrix_algebras.md) | 谱定理、函数演算、极分解、迹范数对偶、GUE 的一个矩 |
| 2 | [有限维冯·诺依曼代数](02_von_neumann_algebras.md) | 双交换子定理、结构定理、中心、对称群与交换子、条件期望 |
| 3 | [迹、态与熵](03_traces_and_states.md) | 迹的分类、密度、扇区的熵、Klein 不等式、在到中心的条件期望下的单调性、Michelson 对比度 |
| 4 | [对称性与表示](04_symmetries.md) | Schur 引理、荷扇区、$`U(1)`$ 与 $`SU(2)`$、自旋系统、实型、复型与四元数型 |
| 5 | [代数中的随机矩阵](05_random_matrices.md) | $`\prod_cU(k_c)`$ 上的 Haar 测度、Weingarten 演算、能级统计以及为何不能混合扇区 |
| 6 | [量子信道](06_quantum_channels.md) | Kraus 表示与 Stinespring 表示、关于迹的对偶、扇区之间的映射、退相干 |
| 7 | [张量积与纠缠](07_tensor_products.md) | 偏迹、Page 公式、荷的融合 |
| 8 | [动力学](08_dynamics.md) | Schrödinger 方程与 Lindblad 方程、扇区间的跃迁、后选择 |
| 9 | [自旋链](09_spin_chains.md) | 磁化与动量扇区、Lanczos 方法、Bethe 拟设检验、无序 |
| 10 | [量子计算](10_quantum_computing.md) | 无退相干子空间、算子代数纠错、带约束的混合算子 |
| 11 | [优化与迹不等式](11_optimisation.md) | 带约束的参数化、采样与优化的对比 |
| 12 | [等变机器学习](12_equivariant_learning.md) | 作为交换子的等变层 |
| 13 | [若尔当代数：量子力学为何使用 *-代数](13_jordan_algebras.md) | 若尔当积与二次表示、自旋因子与洛伦兹锥、Jordan–von Neumann–Wigner 分类、转置与部分转置、若尔当形式的迹不等式 |
| A | [分析与概率中的若干事实](A_analysis.md) | 幂级数与零测集、高斯向量、球面测度、特征值不等式、Isserlis–Wick 公式 |
| B | [Schur–Weyl 对偶与 Weingarten 演算](B_schur_weyl.md) | Schur–Weyl 对偶、Weingarten 函数、Wigner–Eckart 定理、$`O(k)`$ 上的 Haar 测度 |
| C | [特征值密度与 Page 公式](C_eigenvalue_densities.md) | Weyl 积分公式、GOE/GUE、Wishart 密度与诱导密度、Page 公式 |
| D | [量子动力学半群的生成元](D_gksl.md) | GKSL 定理、直和上的生成元、Frigerio 定理、Bény–Kempf–Kribs |
| E | [欧几里得若尔当代数的分类](E_jordan_classification.md) | Peirce 分解、Hurwitz 定理与 Frobenius 定理、坐标化、Albert 代数、Jordan–von Neumann–Wigner 定理 |
| F | [对称锥、序同构与若尔当映射](F_cones_and_order.md) | 平方锥、Koecher–Vinberg 定理、Kadison 定理、$`\bigoplus_cM_{k_c}`$ 的若尔当自同构 |

<a id="prerequisites"></a>

## 预备知识

本课程只假定以下内容，其余一切都加以证明。

- $`\mathbb R`$ 与 $`\mathbb C`$ 上的（有限维）线性代数，以及代数基本定理。
- 标准的实分析与多元微积分：极限、连续性、$`\mathbb R^n`$ 中的紧性、微分、中值定理、Taylor
  公式、反函数定理与隐函数定理。
- 标准的测度论与积分论：Lebesgue 测度、乘积测度与 Fubini–Tonelli 定理、控制收敛与单调
  收敛、$`\mathbb R^n`$ 中的变量替换、推前测度、密度、独立性。
- 初等概率论：期望、方差、独立性、高斯分布及其矩母函数、Markov 不等式与 Chebyshev 不等式、
  大数定律。

不要求具备物理或泛函分析知识。

## 延伸阅读

M. Takesaki, *Theory of Operator Algebras I* (Springer, 2002)；V. Jones,
*von Neumann algebras*，讲义；D. Petz, *Quantum Information Theory and Quantum
Statistics* (Springer, 2008)；J. Watrous, *The Theory of Quantum Information* (Cambridge,
2018)；M. L. Mehta, *Random Matrices* (Elsevier, 2004)。

库 `torch_vn_algebra` 及其在本讲义背景下的使用是[教程](../../tutorial/zh/README.md)的主题。本库的配套论文见
[`man/article/vn_algebra_pytorch.pdf`](https://github.com/andre-a-no/von_neumann_type_i/blob/main/man/article/vn_algebra_pytorch.pdf)。
