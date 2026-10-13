[English](../README.md) · [Español](../es/README.md) · [Français](../fr/README.md) · [Deutsch](../de/README.md) · **中文** · [日本語](../ja/README.md) · [Русский](../ru/README.md)

# 有限维冯·诺依曼代数：基于 `torch_vn_algebra` 的课程

本课程讲授有限维冯·诺依曼代数——即代数
$\bigoplus_cM_{k_c}(\mathbb C)$——的理论，以及它们在量子物理、量子信息、随机矩阵理论和
机器学习中的应用，并在每一步展示如何用 `torch_vn_algebra` 进行计算。

本课程同时面向两类读者：

- **数学家**将看到所用每个结果的完整证明。除标准线性代数、代数基本定理以及分析与概率论中的
  初等事实（极限、凸性、高斯分布的矩母函数）之外，没有任何结论是不加证明而直接引用的；
  凡是某个结果先陈述后证明，或者只加以应用，文中都会明确指出。
- **物理学家和工程师**可以在标有*物理意义*的方框中读到每个概念的物理含义，并获得可以直接
  改写使用的可运行代码。

第 1–3 章的结构相同：带证明的理论；每个结果对应的库代码；用精确结果（定理、闭合公式或独立
计算）对代码进行的检验；一节*可与之比较的工具*，列出完成同一任务的标准工具；以及附有解答的
习题。

课程中的全部代码都由测试套件（`tests/test_tutorial.py`）执行，因此每个例子都能在当前版本的
库中运行。

## 目录

| | 章 | 主要结果 |
|---|---|---|
| 0 | [本课程讲什么](00_overview.md) | 代数作为数据类型：块结构 = 批量张量；代价；已有的工作与新的内容 |
| 1 | [矩阵代数与函数演算](01_matrix_algebras.md) | 谱定理、函数演算、极分解、迹范数对偶、GUE 的一个矩 |
| 2 | [有限维冯·诺依曼代数](02_von_neumann_algebras.md) | 双交换子定理、结构定理、中心、对称群与交换子、条件期望 |
| 3 | [迹、态与熵](03_traces_and_states.md) | 迹的分类、密度、扇区的熵、Klein 不等式、在到中心的条件期望下的单调性、Michelson 对比度 |
| 4 | 对称性与表示 *（编写中）* | Schur 引理、荷扇区、$U(1)$ 与 $SU(2)$、自旋系统 |
| 5 | 代数中的随机矩阵 *（编写中）* | $\prod_cU(k_c)$ 上的 Haar 测度、Weingarten 演算、能级统计以及为何不能混合扇区 |
| 6 | [量子信道](06_quantum_channels.md) | Kraus 表示与 Stinespring 表示、关于迹的对偶、扇区之间的映射、退相干 |
| 7 | [张量积与纠缠](07_tensor_products.md) | 偏迹、Page 公式、荷的融合 |
| 8 | 动力学 *（编写中）* | Schrödinger 方程与 Lindblad 方程、扇区间的跃迁、后选择 |
| 9 | 自旋链 *（编写中）* | 磁化与动量扇区、Lanczos 方法、Bethe 拟设检验、无序 |
| 10 | 量子计算 *（编写中）* | 无退相干子空间、算子代数纠错、带约束的混合算子 |
| 11 | 优化与迹不等式 *（编写中）* | 带约束的参数化、采样与优化的对比 |
| 12 | 等变机器学习 *（编写中）* | 作为交换子的等变层 |

## 预备知识

$\mathbb C$ 上的线性代数（特征值、内积、酉矩阵）、基础概率论以及 Python。不要求具备物理
或泛函分析知识。

## 如何运行代码

```bash
pip install -e .                # from the root of the repository
python -m pytest tests/test_tutorial.py -q
```

一章中的所有代码块可以按顺序粘贴到同一个 Python 会话或 Jupyter 笔记本中运行。同一章的
代码块共享变量。

## 延伸阅读

M. Takesaki, *Theory of Operator Algebras I* (Springer, 2002)；V. Jones,
*von Neumann algebras*，讲义；D. Petz, *Quantum Information Theory and Quantum
Statistics* (Springer, 2008)；J. Watrous, *The Theory of Quantum Information* (Cambridge,
2018)；M. L. Mehta, *Random Matrices* (Elsevier, 2004)。物理背景的数学化表述见
[写给数学家的物理](../../physics_for_mathematicians.zh.md)。本库的配套论文见
[`paper/main.pdf`](https://github.com/andre-a-no/von_neumann_type_i/blob/main/paper/main.pdf)。
