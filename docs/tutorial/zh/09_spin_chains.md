[English](../09_spin_chains.md) · [Español](../es/09_spin_chains.md) · [Français](../fr/09_spin_chains.md) · [Deutsch](../de/09_spin_chains.md) · **中文** · [日本語](../ja/09_spin_chains.md) · [Русский](../ru/09_spin_chains.md)

# 9. 自旋链

由 $L$ 个自旋 $\frac12$ 组成的链位于 $(\mathbb C^2)^{\otimes L}$ 上，而维数为 $2^L$ 的稠密线性代数在
$L\approx14$ 处就止步了。对称性能推高这一极限：磁化把空间分成 $L+1$ 个扇区（第 4 章），平移把每个扇区再分一次，
而 $SU(2)$ 又一次约化这些块。本章为 XXZ 链构造这些扇区，证明它们关于基态说明了什么，用自由费米子和 Bethe
拟设检验本库，并超越稠密矩阵（Krylov 方法）和纯净的链（无序）。

**你将学到**为什么链的代数是 $\bigoplus_NM_{\binom LN}$，为什么这些代数构成的塔以 Pascal 三角形为 Bratteli 图、
以 $4\cos^2\frac\pi{2L+4}$ 为指标，如何计数动量扇区，Marshall 符号规则及其对基态动量与自旋的推论
（Lieb–Mattis），XX 链的 Jordan–Wigner 解，双磁振子 Bethe 方程，为什么 Lanczos 方法单调收敛，以及如何测量
能级统计和自由费米子的纠缠。

在本章中，格点为 $0,\dots,L-1$，$(\mathbb C^2)^{\otimes L}$ 的基由构型
$|b\rangle=|b_0\dots b_{L-1}\rangle$ 组成，自旋向上时 $b_i=1$（本库的约定：标号的第 $i$ 位
对应格点 $i$），$s^a_i$（$a=x,y,z$）是格点 $i$ 上的自旋矩阵 $\frac12\sigma_a$，
$s^\pm_i=s^x_i\pm is^y_i$，$n_i=s^z_i+\frac12$，$S^a=\sum_is^a_i$，而 $N=S^z+\frac L2$ 是向上自旋的
个数。链是**开链**（键 $(i,i+1)$，$i<L-1$）或**周期链**（还有键 $(L-1,0)$；此时
$L\ge3$，指标取模 $L$）。$J>0$（反铁磁），$\hbar=1$，$\log$ 为自然对数。

## 9.1 XXZ 链及其 $U(1)$ 扇区

各向异性为 $\Delta\in\mathbb R$、外场为 $h_i\in\mathbb R$ 的 **XXZ 链**是
$$
H=J\sum_{(i,j)}\Big[\tfrac12\big(s^+_is^-_j+s^-_is^+_j\big)+\Delta\,s^z_is^z_j\Big]+\sum_ih_is^z_i ,
$$
求和遍及所有键；$\frac12(s^+_is^-_j+s^-_is^+_j)=s^x_is^x_j+s^y_is^y_j$ 使一个向上自旋沿键移动。
$\Delta=1$ 时它是 Heisenberg 链 $J\sum\mathbf s_i\cdot\mathbf s_j$，$\Delta=0$ 时是 XX 链。

**命题 9.1.** $[H,S^z]=0$。因此 $H\in U(1)'\cong\bigoplus_{N=0}^LM_{\binom LN}(\mathbb C)$，其中
块 $N$ 作用在由含 $N$ 个向上自旋的构型张成的空间上。

*证明。* 对角项与对角的 $S^z$ 可交换。对 $i\neq j$，$s^+_is^-_j$ 与 $s^z_l$（$l\ne i,j$）可交换，并且由
$[s^z,s^\pm]=\pm s^\pm$（4.5 节），$[s^z_i+s^z_j,s^+_is^-_j]=[s^z_i,s^+_i]s^-_j+s^+_i[s^z_j,s^-_j]=s^+_is^-_j-s^+_is^-_j=0$；
对 $s^-_is^+_j$ 取伴随即可。其余部分是命题 4.10(2)。$\square$

`SpinChain(L).xxz` 直接由每个扇区的构型（`chain.labels[N]`）构造这些块。
我们把它们与用 Kronecker 积构造的完整 $2^L\times2^L$ 矩阵作比较：

```python
import functools, itertools, math, time
import numpy as np
import torch
from scipy.optimize import brentq
from torch_vn_algebra import SpinChain, SU2Chain, DensityMatrix, inclusions, krylov

torch.manual_seed(0)
f64, c128 = torch.float64, torch.complex128
close = lambda a, b, tol=1e-10: bool((torch.as_tensor(a) - torch.as_tensor(b)).abs().max() < tol)
def op(L, ops):                             # product of one-site operators {site: 2 x 2}; site 0 = lowest bit
    return functools.reduce(torch.kron, [ops.get(j, torch.eye(2, dtype=f64)) for j in reversed(range(L))])

sp = torch.tensor([[0., 0.], [1., 0.]], dtype=f64)                    # s^+ : |0> = down -> |1> = up
sm, sz = sp.T.contiguous(), torch.diag(torch.tensor([-0.5, 0.5], dtype=f64))
L = 6
chain = SpinChain(L, boundary='periodic')
h = 0.3 * torch.randn(L, dtype=f64)                                  # random fields
H_full = sum(0.5 * (op(L, {i: sp, j: sm}) + op(L, {i: sm, j: sp})) + 0.7 * op(L, {i: sz, j: sz})
             for i, j in chain.bonds) + sum(h[i] * op(L, {i: sz}) for i in range(L))   # J = 1, Delta = 0.7
Sz = sum(op(L, {i: sz}) for i in range(L))
assert close(H_full @ Sz, Sz @ H_full)                               # Proposition 9.1
H = chain.xxz(J=1.0, Delta=0.7, h=h)                                 # an element of (+)_N M_binom(6, N)
for N, k in enumerate(chain.algebra.k_factors):
    lab = chain.labels[N]                                            # configurations of sector N
    assert close(H.matrix[0, N, :k, :k], H_full[lab][:, lab])
```

**链的塔。** 令 $\mathcal A_L=\bigoplus_NM_{\binom LN}$ 作用在 $L$ 个格点上。添加一个格点，
$x\mapsto x\otimes1_2$，把 $\mathcal A_L$ 映入 $\mathcal A_{L+1}$，因为 $x\otimes1$ 与
$S^z\otimes1+1\otimes s^z$ 可交换。这就给出一个塔 $\mathcal A_0\subset\mathcal A_1\subset\mathcal A_2\subset\cdots$。

**命题 9.2.** 包含 $\mathcal A_L\subset\mathcal A_{L+1}$ 的包含矩阵是 $(L+2)\times(L+1)$
矩阵 $\Lambda_L$，其中 $\lambda_{N,N}=\lambda_{N,N-1}=1$，其余元素均为 $0$；块的
大小服从 Pascal 法则 $\binom{L+1}N=\binom LN+\binom L{N-1}$，塔的 Bratteli 图是
Pascal 三角形，而 $\Lambda_L$ 的二部图是一条有 $2L+3$ 个顶点的路径。

*证明。* $L+1$ 个格点上含 $N$ 个向上自旋的构型，一类是前 $L$ 个格点上有 $N$ 个向上自旋而新自旋向下，
另一类是前 $L$ 个格点上有 $N-1$ 个向上自旋而新自旋向上。在这两个子空间上
$x\otimes1$ 分别作用为 $x_N$ 和 $x_{N-1}$，所以 $\mathcal A_{L+1}$ 的块 $N$ 各包含一次 $\mathcal A_L$ 的
块 $N$ 与块 $N-1$；维数相加正是 Pascal 法则。在二部图中，$\mathcal A_{L+1}$ 的扇区 $N$
与 $\mathcal A_L$ 的 $N-1$ 和 $N$ 相连，这就给出路径
$0'-0-1'-1-\cdots-L'-L-(L+1)'$（带撇号的属于 $\mathcal A_{L+1}$）。$\square$

从 Bratteli 图顶端到第 $L$ 层扇区 $N$ 的一条路径就是一个构型 $b$；这样的路径共有
$\binom LN$ 条，恰为块的大小。对包含矩阵为 $\Lambda$ 的连通包含，我们取
$\|\Lambda\|^2$ 作为**指标**，即 Markov 迹的指标（V. Jones, *Invent. Math.* 72 (1983)；F. Goodman,
P. de la Harpe, V. Jones, *Coxeter graphs and towers of algebras*, Springer 1989, Ch. 2–3）；对
$M_k\otimes1\subset M_k\otimes M_m$，它就是 7.2 节的 $m^2$。`Inclusion.index` 计算它。

**引理 9.3.** 顶点为 $1,\dots,n$ 的路径的邻接矩阵 $A_n$ 的特征值为
$2\cos\frac{\pi k}{n+1}$，$k=1,\dots,n$，对应的特征向量为 $v_k(j)=\sin\frac{\pi jk}{n+1}$；因此
$\|A_n\|=2\cos\frac\pi{n+1}$。

*证明。* 令 $v_k(0)=v_k(n+1)=0$。由
$\sin(a-b)+\sin(a+b)=2\sin a\cos b$，$(A_nv_k)(j)=v_k(j-1)+v_k(j+1)=2\cos\frac{\pi k}{n+1}\,v_k(j)$，并且由于 $v_k(1)=\sin\frac{\pi k}{n+1}\neq0$，有 $v_k\ne0$。这 $n$ 个特征值
互不相同（$\cos$ 在 $(0,\pi)$ 上严格递减），所以它们就是全部特征值，而对称矩阵
$A_n$ 的范数是特征值的最大模，在 $k=1$ 和 $k=n$ 处取到。$\square$

**推论 9.4.** $[\mathcal A_{L+1}:\mathcal A_L]=\|\Lambda_L\|^2=4\cos^2\frac\pi{2L+4}$，它递增趋于 $4$。

*证明。* 沿命题 9.2 的路径给顶点排序，邻接矩阵
$B=\begin{pmatrix}0&\Lambda_L\\\Lambda_L^T&0\end{pmatrix}$ 就变成 $A_{2L+3}$（一个置换相似）。由引理 1.1，
$\|B\|^2=\|B^2\|$ 且 $\|\Lambda_L\Lambda_L^T\|=\|\Lambda_L^T\Lambda_L\|=\|\Lambda_L\|^2$，而 $B^2=\Lambda_L\Lambda_L^T\oplus\Lambda_L^T\Lambda_L$，所以
$\|B\|=\|\Lambda_L\|$；取 $n=2L+3$ 的引理 9.3 给出其值。$\square$

数 $4\cos^2\frac\pi n$（$n\ge3$）恰好是子因子的指标在 $4$ 以下所能取到的值
（Jones，前引文献）；这里出现的是偶数 $n$，原因是引理 9.3 这一初等事实。

```python
tower = inclusions.chain_tower(8)                                    # A_0 c A_1 c ... c A_8
print(inclusions.bratteli_diagram(tower).levels[:6])                 # Pascal's triangle
for L_, inc in enumerate(tower):
    sizes = [math.comb(L_ + 1, N) for N in range(L_ + 2)]
    assert inc.M.k_factors == sizes and inc.is_connected
    assert list(inc.matrix @ np.array(inc.N.k_factors)) == sizes      # unitality = Pascal's rule
    n = 2 * L_ + 3
    A_path = np.eye(n, k=1) + np.eye(n, k=-1)
    assert abs(np.linalg.eigvalsh(A_path).max() - 2 * math.cos(math.pi / (n + 1))) < 1e-12   # Lemma 9.3
    assert abs(inc.index() - 4 * math.cos(math.pi / (2 * L_ + 4)) ** 2) < 1e-12           # Corollary 9.4
```

> **物理意义。** 只要相互作用在绕 $z$ 轴的转动下不变，$S^z$ 就是守恒的；$N$ 计数向上自旋，
> 或者在 9.4 节的 Jordan–Wigner 映射之后计数粒子。在塔中，较长链的每个荷扇区
> 都由较短链的两个扇区供给。

## 9.2 平移、动量与 Marshall 符号规则

在周期链上，令 $T$ 为平移 $T|b_0b_1\dots b_{L-1}\rangle=|b_{L-1}b_0\dots b_{L-2}\rangle$，
它把格点 $i$ 移到格点 $i+1$。$T$ 是置换矩阵，所以是酉的，$T^L=1$ 且
$Ts^a_iT^*=s^a_{i+1}$。因此 $TS^zT^*=S^z$，并且对均匀外场 $h_i=h$ 有 $THT^*=H$：$T$ 把键映为键。
它的特征值为 $e^{ik}$，其中**动量** $k=2\pi m/L$，$m=0,\dots,L-1$。

**命题 9.5（动量扇区）。** 固定 $N$。对含 $N$ 个向上自旋的构型 $a$，令 $R_a$ 为
其轨道 $\{T^ja\}$ 的大小。(a) $R_a$ 整除 $L$；对满足 $mR_a\equiv0\pmod L$ 的 $k=2\pi m/L$，向量
$|a,k\rangle=\sum_{j=0}^{R_a-1}e^{-ikj}T^j|a\rangle$ 非零且 $T|a,k\rangle=e^{ik}|a,k\rangle$，并且这
$R_a$ 个向量构成轨道张成空间的一组正交基。(b) $T$ 在扇区 $N$ 上的特征值 $e^{2\pi im/L}$
的重数为 $d_{N,m}=\#\{\text{轨道 }:mR_a\equiv0\pmod L\}$，并且每个与
$S^z$ 和 $T$ 都可交换的算子——例如平移不变的 XXZ 链——都是块对角的，属于 $\bigoplus_N\bigoplus_mM_{d_{N,m}}$。

*证明。* (a) $a$ 在循环群 $\{T^j\}\cong\mathbb Z_L$ 中的稳定子是一个子群，而 $R_a$ 是它的
指数，是 $L$ 的因子；向量 $T^j|a\rangle$（$0\le j<R_a$）规范正交，且 $T^{R_a}|a\rangle=|a\rangle$。若
$e^{ikR_a}=1$，则 $T|a,k\rangle=e^{ik}\sum_{j=1}^{R_a}e^{-ikj}T^j|a\rangle=e^{ik}|a,k\rangle$，因为 $j=R_a$ 的项等于
$j=0$ 的项；$|a\rangle$ 的系数为 $1$。容许的 $m\in\{0,\dots,L-1\}$ 是
$L/R_a$ 的倍数，共 $R_a$ 个，对应酉算子 $T$ 的互不相同的特征值，所以这些向量相互正交，并张成
轨道所张成的 $R_a$ 维空间。(b) 扇区 $N$ 是其各轨道张成空间的正交和。与 $S^z$ 和 $T$ 都可交换的 $X$
与它们的谱投影 $P_N$、$Q_m$ 可交换（定理 1.5(d)），所以
$X=\sum_{N,m}P_NQ_mXP_NQ_m$。$\square$

`xxz_momentum` 在归一化的基 $|a,k\rangle$ 中构造这些块，每条轨道取一个代表元 $a$。我们
独立地计数轨道，并把谱与稠密扇区作比较：

```python
L = 12
ring = SpinChain(L, boundary='periodic')
rot = lambda x: ((x << 1) | (x >> (L - 1))) & ((1 << L) - 1)        # T on labels: site i -> i + 1
for N in range(L + 1):
    orbits = {frozenset(itertools.accumulate(range(L - 1), lambda x, _: rot(x), initial=a))
              for a in ring.labels[N].tolist()}                       # {a, Ta, ..., T^{L-1} a}
    d = [sum((m * len(o)) % L == 0 for o in orbits) for m in range(L)]
    assert d == ring.momentum_dims(N)                                   # Proposition 9.5(b)
print(ring.momentum_dims(6))
Hk = ring.xxz_momentum(J=1.0, Delta=1.0, h=0.2, sector=6)               # 12 blocks instead of one of size 924
E_k = torch.sort(torch.cat([e[0] for e in Hk.eigenvalues()])).values
assert close(E_k, torch.linalg.eigvalsh(ring.xxz(1.0, 1.0, h=0.2, sector=6).matrix[0, 0]))
```

含 $924$ 个态的扇区 $N=6$ 分裂成大小为 $80,75,78,76,78,75,80,75,78,76,78,75$ 的块。

**引理 9.6（Perron–Frobenius，对称情形）。** 设 $M$ 是实对称矩阵，对
$i\ne j$ 有 $M_{ij}\ge0$，且其图（$i\sim j$ 当且仅当 $M_{ij}>0$）连通。则其最大特征值 $\lambda$ 是单的，
并有一个各分量全都 $>0$ 的特征向量。

*证明。* $\lambda=\max_{\|x\|=1}\langle x,Mx\rangle$，且最大值恰好在
$\lambda$-特征空间的单位向量上取到（在一组规范正交特征基中展开）。对这样的 $x$，由各分量模构成的向量 $|x|$
满足 $\langle|x|,M|x|\rangle\ge\operatorname{Re}\langle x,Mx\rangle=\lambda$，因为对角项相同，且对 $i\ne j$ 有
$M_{ij}|x_i||x_j|\ge\operatorname{Re}M_{ij}\bar x_ix_j$；所以 $|x|$ 是特征向量。若 $|x|_i=0$，则
$0=(\lambda-M_{ii})|x|_i=\sum_{j\ne i}M_{ij}|x|_j$ 迫使 $i$ 的所有邻点 $j$ 都有 $|x|_j=0$，再由
连通性得 $x=0$。所以非零特征向量没有为零的分量；但维数 $\ge2$ 的特征空间
含有一个满足 $x_1=0$ 的非零向量。因此 $\lambda$ 是单的，且 $|x|>0$ 是特征向量。$\square$

**定理 9.7（Marshall 符号规则）。** 设链为开链，或 $L$ 为偶数的周期链；令 $B$ 为
奇数格点的集合，$N_B(b)$ 为 $b$ 在 $B$ 上的向上自旋数；设 $J>0$，$\Delta,h_i\in\mathbb R$ 任意。
在每个扇区 $N$ 中，$H_N$ 的最低特征值是单的，并且其特征向量可以选得使
对每个含 $N$ 个向上自旋的构型 $b$ 都有 $(-1)^{N_B(b)}\psi(b)>0$。

*证明。* 令 $V|b\rangle=(-1)^{N_B(b)}|b\rangle$；$V$ 是对角的，$V=V^*=V^{-1}$。每条键都连接一个偶数格点和一个
奇数格点（对环而言是因为 $L$ 为偶数），所以 $s^+_is^-_j$ 使 $N_B$ 改变 $\pm1$，且
$Vs^+_is^-_jV=-s^+_is^-_j$。于是 $M=-VH_NV$ 在相差一条键上两个自旋交换的构型之间
的非对角元为 $J/2>0$，其余非对角元为 $0$。它的图是连通的：相邻格点 $(i,i+1)$（$i<L-1$）的交换
生成格点的所有置换，因而连接任意两个含 $N$ 个向上自旋的构型
（交换两个相同的自旋不起作用）。由引理 9.6，$M$ 的最大特征值——即 $H_N$ 最低特征值的
相反数——是单的，并有特征向量 $\varphi>0$；取 $\psi=V\varphi$。$\square$

**推论 9.8（基态的动量）。** 在 $L$ 为偶数、$J>0$ 且外场均匀的周期链上，
扇区 $N$ 的基态在 $N$ 为偶数时动量为 $k=0$，在 $N$ 为奇数时动量为 $k=\pi$。

*证明。* $T$ 与 $H_N$ 可交换，且基态 $\psi=V\varphi$ 是单的，所以 $T\psi=e^{ik}\psi$。而
$(T\psi)(b)=\psi(T^{-1}b)$，且 $T^{-1}b$ 在格点 $i$ 上的自旋为 $b_{i+1}$；由于 $L$ 为偶数，对奇数 $i$，$i+1$ 是偶数，
所以 $N_B(T^{-1}b)=N-N_B(b)$，且 $(VT\psi)(b)=(-1)^{N_B(b)}(-1)^{N-N_B(b)}\varphi(T^{-1}b)=(-1)^N\varphi(T^{-1}b)$。另一方面
$VT\psi=e^{ik}\varphi$。$\varphi$ 与 $\varphi\circ T^{-1}$ 都是正的，所以 $e^{ik}=(-1)^N$。$\square$

```python
for bc in ('open', 'periodic'):
    c = SpinChain(10, bc, complex_valued=False)
    for N in range(11):
        w, V = torch.linalg.eigh(c.xxz(1.0, 0.7, h=0.3 * torch.randn(10, dtype=f64), sector=N).matrix[0, 0])
        x = (-1.0) ** c.bits(N)[:, 1::2].sum(-1) * V[:, 0]                  # (-1)^{N_B(b)} psi(b)
        assert (x * x.sum().sign()).min() > 0 and (len(w) == 1 or w[1] - w[0] > 1e-8)   # Theorem 9.7
for N in range(1, L):
    lowest = [e[0, 0].item() for e in ring.xxz_momentum(1.0, 0.5, sector=N).eigenvalues()]
    m0 = int(np.argmin(lowest))
    assert m0 == (0 if N % 2 == 0 else L // 2)                              # Corollary 9.8: k = 0 or pi
    assert sorted(lowest)[1] - lowest[m0] > 1e-8
```

> **物理意义。** 在反铁磁体中，晶格决定了基态的符号：使一个自旋越过一条键会
> 改变符号。这使得二部反铁磁体的量子 Monte Carlo 没有符号问题；
> 对于阻挫链，即有耦合连接同一子格的格点时，这一点不再成立。

## 9.3 $\Delta=1$ 时的 $SU(2)$

对 $\Delta=1$ 和 $h=0$，$H=J\sum\mathbf s_i\cdot\mathbf s_j$ 与 $S^x,S^y,S^z$ 可交换（命题 4.15），所以
$H\cong\bigoplus_S1_{2S+1}\otimes H_S$，其中块 $H_S\in M_{m_S}$，$m_S=\binom L{L/2-S}-\binom L{L/2-S-1}$
（命题 4.13），而模块 `su2` 在最高权向量上构造 $H_S$（命题 4.15(b)）。
令 $E_0(S)$ 为 $H_S$ 的最低特征值。

**定理 9.9（链的 Lieb–Mattis 序）。** 设链如定理 9.7 中所述，$\Delta=1$，$h_i=0$，
$J>0$。对 $L$ 个自旋的每个总自旋 $S$，扇区 $S^z=S$ 的基态具有总自旋 $S$。因此
$E_0(S)$ 是 $H_S$ 的单特征值，并且 $E_0(S_{\min})<E_0(S_{\min}+1)<\dots<E_0(\frac L2)$，其中
$S_{\min}\in\{0,\frac12\}$：基态是自旋最小的单个多重态。

*证明。* 设 $\psi$ 为扇区 $S^z=S$ 的基态，由定理 9.7 它是单的。$\mathbf S^2$ 与
$H$ 可交换并保持该扇区，所以 $\mathbf S^2\psi$ 位于一维的基态空间中：$\psi$ 具有某个总
自旋 $S'$。令 $\chi$ 为在格点对
$(0,1),(2,3),\dots,(L-2S-2,L-2S-1)$ 上的单态 $\frac1{\sqrt2}(|{\uparrow\downarrow}\rangle-|{\downarrow\uparrow}\rangle)$ 之积，其余 $2S$ 个格点自旋向上。$S^+$ 湮灭每个因子，
且 $S^z\chi=S\chi$，所以 $\mathbf S^2\chi=(S^-S^++(S^z)^2+S^z)\chi=S(S+1)\chi$（4.5 节）。在每一对上，向上自旋
位于偶数格点（系数 $+$）或奇数格点（$-$），所以 $(-1)^{N_B(b)}\chi(b)$ 对其支撑中的所有 $b$
符号相同；取定理 9.7 中的 $V,\varphi$，有 $\langle\psi,\chi\rangle=\langle\varphi,V\chi\rangle\neq0$。
$\mathbf S^2$ 属于不同特征值的特征向量相互正交，所以 $S'=S$。由命题 4.15(a) 和推论 4.12，每个
自旋 $S'\ge S$ 的多重态与扇区 $S^z=S$ 恰好交于一个向量，所以 $H$ 在该扇区的谱是
各 $H_{S'}$（$S'\ge S$）的谱的并（计重数）。它的最低特征值是单的，其特征向量的
自旋为 $S$，所以它就是 $E_0(S)$，在 $H_S$ 中是单的，且不是 $H_{S'}$（$S'>S$）的特征值：$E_0(S)<E_0(S')$。$\square$

Lieb 和 Mattis 证明了更多：对子格自旋为 $S_A,S_B$ 的连通二部格子，若耦合在
子格之间为反铁磁、在子格内部为铁磁，则 $E_0(S)$ 在 $S\ge|S_A-S_B|$ 时递增，
且基态的自旋为 $|S_A-S_B|$（E. Lieb, D. Mattis, *J. Math. Phys.* 3 (1962) 749；不加
证明地陈述）。由定理 9.9，**单态–三重态能隙** $E_0(1)-E_0(0)$ 为正；在环上
$L(E_0(1)-E_0(0))\to\pi^2/2$，带有对数修正（J. des Cloizeaux, J. Pearson, *Phys. Rev.* 128
(1962) 2131；I. Affleck, D. Gepner, H. Schulz, T. Ziman, *J. Phys. A* 22 (1989) 511；不加证明地陈述）。
我们比较十二个自旋时的块大小，用 $U(1)$ 扇区检验定理 9.9，并计算能隙：

```python
sym = SU2Chain(12, boundary='periodic', complex_valued=False)
H_su2 = sym.heisenberg(J=1.0)                                       # (+)_S M_{m_S}, S = 0, ..., 6
u1_sizes = [math.comb(12, N) for N in range(13)]
print(sym.multiplicities, max(u1_sizes), sum(k ** 3 for k in u1_sizes), sum(m ** 3 for m in sym.multiplicities))
E0_S = sym.lowest_energies(H_su2)[0]                                # E_0(S), S = 0, ..., 6
E0_N = [torch.linalg.eigvalsh(ring.xxz(1.0, 1.0, sector=6 + S).matrix[0, 0])[0] for S in range(7)]
assert close(E0_S, torch.stack(E0_N))                               # Theorem 9.9: sector S^z = S has spin S
assert (E0_S[1:] > E0_S[:-1] + 1e-8).all()                          # E_0(S) strictly increasing
gap = {}
for n in (6, 8, 10, 12):
    s = SU2Chain(n, 'periodic', complex_valued=False)
    e = [torch.linalg.eigvalsh(s.heisenberg(spin=S).matrix[0, 0])[0].item() for S in (0, 1)]
    gap[n] = n * (e[1] - e[0])
print({n: round(g, 4) for n, g in gap.items()}, round(math.pi ** 2 / 2, 4))
assert 0 < gap[6] < gap[8] < gap[10] < gap[12] < math.pi ** 2 / 2
```

$SU(2)$ 块的大小为 $m_S=132,297,275,154,54,11,1$，而最大的 $U(1)$ 块为 $924$，
对角化的开销 $\sum k^3$ 从 $2.0\cdot10^9$ 降到 $5.3\cdot10^7$。$L$ 乘以能隙增长得很慢，
对 $L=6,\dots,12$ 依次为 $4.108,4.181,4.232,4.270$，趋向 $\pi^2/2\approx4.935$。

**关联。** 各格点的自旋构成一个秩一张量算子，由 **Wigner–Eckart 定理**，它们在
多重态之间的矩阵元等于一个 Clebsch–Gordan 系数乘以一个与磁量子数无关的**约化矩阵元**
$\langle S',a\|\mathbf s_i\|S,b\rangle$（A. R. Edmonds, *Angular
Momentum in Quantum Mechanics*, Princeton 1957, Ch. 5；不加证明地陈述）。`spin_correlation` 利用它
仅从最高权向量出发在每个 $H_S$ 上构造 $\mathbf s_i\cdot\mathbf s_j$；`'direct'` 方法
则改为投影键算子。我们把两者与在 $U(1)$ 扇区中显式写出的基态作比较：

```python
gs = sym.ground_state(H_su2)
S0, coef = gs.spin.item(), gs.coefficients[0, :sym.multiplicity(gs.spin.item())]
pairs = [(0, r) for r in range(12)]
corr_we = sym.expectation(sym.spin_correlation(pairs), S0, coef)
corr_dir = sym.expectation(sym.spin_correlation(pairs, method='direct'), S0, coef)
psi = sym.embed(S0, coef)                                            # the ground state in the sector N = 6
lab, z = sym.u1.labels[6], sym.u1.bits(6).to(f64) - 0.5
def s_dot_s(i, j):                                                   # <psi| s_i . s_j |psi>, from the definition
    if i == j:
        return 0.75
    flip = ((lab >> i) & 1) != ((lab >> j) & 1)
    partner = sym.u1.index_of[lab[flip] ^ ((1 << i) | (1 << j))]
    return ((psi ** 2 * z[:, i] * z[:, j]).sum() + 0.5 * (psi[flip] * psi[partner]).sum()).item()
assert S0 == 0 and close(corr_we, corr_dir)
assert close(corr_we, torch.tensor([s_dot_s(0, r) for r in range(12)], dtype=f64))
assert ((-1.0) ** torch.arange(1, 12) * corr_we[1:] > 0).all()       # signs (-1)^r
print([round(x, 4) for x in corr_we.tolist()])
```

关联 $\langle\mathbf s_0\cdot\mathbf s_r\rangle=0.75,-0.4489,0.1878,-0.1658,0.1209,-0.1224,0.1070,\dots$
符号交替，这正是定理 9.7 和 9.9 所强制的：单态在转动下不变，而转动置换
坐标轴，并且在扇区 $S^z=0$ 中 $\langle s^+_0s^+_r\rangle=0$，所以 $\langle\mathbf s_0\cdot\mathbf s_r\rangle=3\langle s^x_0s^x_r\rangle=
\frac34\langle s^+_0s^-_r+s^-_0s^+_r\rangle=\frac34\sum_{b_0\ne b_r}\psi(b)\psi(b')$，其中 $b'$ 由 $b$ 交换
格点 $0$ 与 $r$ 的自旋得到；由于 $N_B(b')-N_B(b)$ 对奇数 $r$ 为 $\pm1$、对偶数 $r$ 为 $0$，每一项的符号都是
$(-1)^r$。衰减很慢，当 $L\to\infty$ 时如 $(-1)^r\sqrt{\log r}/r$（Affleck 等，前引文献；不加
证明地陈述）。

## 9.4 精确检验：自由费米子与 Bethe 拟设

**Jordan–Wigner。** 在开链上令 $c_j=\big(\prod_{l<j}(-1)^{n_l}\big)s^-_j$。

**定理 9.10（XX 链是自由的）。** (a) $\{c_i,c_j^*\}=\delta_{ij}$，$\{c_i,c_j\}=0$，且 $c_j^*c_j=n_j$。
(b) 开放 XX 链为 $H=\frac J2\sum_{j<L-1}(s^+_js^-_{j+1}+s^-_js^+_{j+1})=\sum_{i,j}h_{ij}c_i^*c_j$，其中
$h=\frac J2A_L$。(c) 设 $u_k$ 是 $h$ 的规范正交实特征向量，特征值为
$\varepsilon_k=J\cos\frac{\pi k}{L+1}$（引理 9.3），并令 $b_k=\sum_ju_k(j)c_j$。对 $K=\{k_1<\dots<k_N\}$，向量
$b_K^*|0\rangle=b_{k_1}^*\cdots b_{k_N}^*|0\rangle$（$|0\rangle$ 为所有自旋向下）构成扇区
$N$ 的一组规范正交特征基，能量为 $\sum_{k\in K}\varepsilon_k$。

*证明。* (a) 不同格点上的算子可交换，而在同一格点上 $(-1)^n$ 与 $s^\pm$ 反对易。对
$i<j$，$c_i$ 与 $c_j^*$（或 $c_j$）在格点 $i$ 上以 $s^-_i$ 和 $(-1)^{n_i}$ 的形式相遇——一次反对易——
在其他地方可交换，所以它们反对易。对 $i=j$，弦算子的平方为 $1$，且
$c_jc_j^*+c_j^*c_j=s^-_js^+_j+s^+_js^-_j=1$，$c_j^*c_j=s^+_js^-_j=n_j$。(b) $c_j^*c_{j+1}=s^+_j(-1)^{n_j}s^-_{j+1}=s^+_js^-_{j+1}$，
因为 $s^+_j$ 仅在 $n_j=0$ 上非零；再取伴随。(c) $\{b_k,b_l^*\}=\sum_ju_k(j)u_l(j)=\delta_{kl}$ 且
$\{b_k,b_l\}=0$，并由 $h$ 的谱分解得 $H=\sum_k\varepsilon_kb_k^*b_k$。由正则反对易关系（CAR），
$[b_l^*b_l,b_k^*]=\delta_{lk}b_k^*$，所以 $[H,b_k^*]=\varepsilon_kb_k^*$；结合 $H|0\rangle=0$ 得
$Hb_K^*|0\rangle=(\sum_K\varepsilon_k)b_K^*|0\rangle$。由于 $b_k|0\rangle=0$（每个 $c_j$ 降低一个自旋），用 CAR 把湮灭算子
移到右边，对 $N$ 归纳即得 $\langle b_{K'}^*0,b_K^*0\rangle=\delta_{KK'}$。每个 $c_j^*$ 使
$N$ 增加一，所以 $b_K^*|0\rangle$ 位于扇区 $|K|$ 中：在维数为 $\binom LN$ 的空间中有 $\binom LN$ 个规范正交向量。$\square$

```python
L = 10
xx = SpinChain(L)
eps = [math.cos(math.pi * k / (L + 1)) for k in range(1, L + 1)]
for N in range(L + 1):
    E = torch.linalg.eigvalsh(xx.xxz(J=1.0, Delta=0.0, sector=N).matrix[0, 0])
    E_ff = sorted(sum(K) for K in itertools.combinations(eps, N))     # sums over N-element subsets
    assert close(E, torch.tensor(E_ff, dtype=f64))                    # all 1024 levels, Theorem 9.10
    assert abs(E[0].item() - sum(sorted(eps)[:N])) < 1e-12           # ground state: the N lowest modes
```

**磁振子。** 在 $h=0$ 的环上，完全极化态 $|0\rangle$ 的能量为 $E_F=J\Delta L/4$。位于 $x$ 的一个向上
自旋破坏两条键：$H|x\rangle=(E_F-J\Delta)|x\rangle+\frac J2(|x-1\rangle+|x+1\rangle)$（$L\ge3$）。

**命题 9.11（单磁振子）。** 在环的扇区 $N=1$ 中，$|k\rangle=\sum_xe^{-ikx}|x\rangle$，
$k=2\pi m/L$，满足 $T|k\rangle=e^{ik}|k\rangle$ 和 $H|k\rangle=\big(E_F+J(\cos k-\Delta)\big)|k\rangle$。

*证明。* 利用 $e^{ikL}=1$，$\sum_xe^{-ikx}(|x-1\rangle+|x+1\rangle)=\sum_y(e^{-ik(y+1)}+e^{-ik(y-1)})|y\rangle=2\cos k\,|k\rangle$；
同样的重新标号给出 $T|k\rangle=e^{ik}|k\rangle$。$\square$

两个磁振子会发生散射。对 $x_1<x_2$，用 $|x_1,x_2\rangle$ 表示向上自旋位于 $x_1,x_2$ 的态；一对相邻的向上自旋
破坏两条键，一对分开的向上自旋破坏四条。

**命题 9.12（双磁振子；Bethe 方程）。** 设 $k_1,k_2\in\mathbb C$，$K=k_1+k_2$，且
$(A_{12},A_{21})\neq0$ 满足 (i) $A_{12}(1+e^{iK}-2\Delta e^{ik_2})+A_{21}(1+e^{iK}-2\Delta e^{ik_1})=0$ 和
(ii) $A_{21}=e^{ik_2L}A_{12}$，$A_{12}=e^{ik_1L}A_{21}$。令
$\psi(x_1,x_2)=A_{12}e^{i(k_1x_1+k_2x_2)}+A_{21}e^{i(k_2x_1+k_1x_2)}$（定义在 $\mathbb Z^2$ 上），并令
$\Psi=\sum_{0\le x_1<x_2<L}\psi(x_1,x_2)|x_1,x_2\rangle$。若 $\Psi\neq0$，则
$H\Psi=\big(E_F+J(\cos k_1+\cos k_2-2\Delta)\big)\Psi$。若分母不为零，则 (i) 和 (ii) 等价于
$$
e^{ik_1L}=-\frac{1+e^{iK}-2\Delta e^{ik_1}}{1+e^{iK}-2\Delta e^{ik_2}},\qquad
e^{ik_2L}=-\frac{1+e^{iK}-2\Delta e^{ik_2}}{1+e^{iK}-2\Delta e^{ik_1}} .
$$

*证明。* 令 $E=E_F+J(\cos k_1+\cos k_2-2\Delta)$。(1) 对所有 $(x_1,x_2)\in\mathbb Z^2$，每个平面波给出的
$\psi(x_1\pm1,x_2)$ 与 $\psi(x_1,x_2\pm1)$ 之和为 $2(\cos k_1+\cos k_2)\psi(x_1,x_2)$。(2) 代入可知，
$\psi(x,x)+\psi(x+1,x+1)-2\Delta\psi(x,x+1)$ 等于 $e^{iKx}$ 乘以 (i) 的左边，所以它对所有 $x$ 为零。
(3) 由 (ii)，$\psi(x_2,x_1+L)=A_{12}e^{ik_2L}e^{i(k_2x_1+k_1x_2)}+A_{21}e^{ik_1L}e^{i(k_1x_1+k_2x_2)}=\psi(x_1,x_2)$，
两次应用此式得 $\psi(x_1+L,x_2+L)=\psi(x_1,x_2)$。构型 $\{y_1<y_2\}$ 在
$D=\{x_1<x_2<x_1+L\}$ 中的代表元是 $(y_1+nL,y_2+nL)$ 与 $(y_2+nL,y_1+(n+1)L)$，所以 $\psi$ 是构型的函数。
设 $(x_1,x_2)\in D$。若 $2\le x_2-x_1\le L-2$，四种跳跃都留在 $D$ 中，对角元为 $E_F-2J\Delta$，所以由 (1)，
$(H\Psi)(x_1,x_2)=(E_F-2J\Delta)\psi+\frac J2\cdot2(\cos k_1+\cos k_2)\psi=E\psi$。若 $x_2=x_1+1$，跳跃
$(x_1+1,x_2)$、$(x_1,x_2-1)$ 被阻挡，对角元为 $E_F-J\Delta$；由 (1) 和 (2)，
$(H\Psi)(x_1,x_2)=(E_F-J\Delta)\psi+\frac J2\big(2(\cos k_1+\cos k_2)\psi-\psi(x_1,x_1)-\psi(x_1+1,x_1+1)\big)=E\psi$。
若 $x_2=x_1+L-1$，则用代表元 $(x_2,x_1+L)$，它属于前一种情形（$L\ge3$）。$\square$

命题 9.12 只说明解给出特征向量；它并没有说所有特征向量都以这种方式得到。
记 $k_{1,2}=K/2\pm q$ 和 $c=\cos\frac K2$。则 $1+e^{iK}-2\Delta e^{ik_{1,2}}=2e^{iK/2}(c-\Delta e^{\pm iq})$，
所以 $A_{12}=c-\Delta e^{iq}$，$A_{21}=-(c-\Delta e^{-iq})$ 是 (i) 的解。把 (ii) 的两个方程相乘得
$e^{iKL}=1$，$K=2\pi n/L$，于是 (ii) 等价于 $F(q)=(-1)^ne^{iqL}(c-\Delta e^{-iq})+c-\Delta e^{iq}=0$，
这是一个单变量方程。我们在
$q\in(0,\pi)$ 中寻找它的根，这里由于 $c-\Delta e^{\pm iq}$ 互为复共轭，$e^{-iqL/2}F(q)$ 为实数或纯虚数
（散射态）；我们也在直线 $q=i\kappa$、$\pi+i\kappa$（$\kappa>0$）上寻找，那里 $F$ 是实的（束缚态）；
对 $L=9$、$\Delta=0.6$，我们为每个根构造 $\Psi$，并与含 $36$ 个态的扇区 $N=2$ 作比较：

```python
L, Delta = 9, 0.6
ring9 = SpinChain(L, boundary='periodic')
H2 = ring9.xxz(J=1.0, Delta=Delta, sector=2).matrix[0, 0]
x1, x2 = np.triu_indices(L, 1)                                       # all pairs x1 < x2
rows = ring9.index_of[torch.tensor((1 << x1) | (1 << x2))]
def roots(g, grid):                                                  # sign changes of g on the grid, refined
    v = np.array([g(t) for t in grid])
    return [brentq(g, grid[i], grid[i + 1], xtol=1e-15) for i in np.nonzero(v[:-1] * v[1:] < 0)[0]]
E_bethe, vecs, n_bound = [], [], 0
for n in range(L):
    K, c = 2 * np.pi * n / L, math.cos(np.pi * n / L)
    F = lambda q: (-1) ** n * np.exp(1j * q * L) * (c - Delta * np.exp(-1j * q)) + c - Delta * np.exp(1j * q)
    G = lambda q: (lambda z: z.real + z.imag)(np.exp(-0.5j * q * L) * F(q))    # real, or imaginary
    bound = [s + 1j * t for s in (0, np.pi)                          # q = i kappa and pi + i kappa
             for t in roots(lambda t: F(s + 1j * t).real, np.linspace(1e-3, 20, 4000))]
    n_bound += len(bound)
    for q in roots(G, (np.arange(1, 2000) + 0.5 ** 0.5) * np.pi / 2000) + bound:
        psi = np.exp(0.5j * K * (x1 + x2)) * ((c - Delta * np.exp(1j * q)) * np.exp(1j * q * (x1 - x2))
                                              - (c - Delta * np.exp(-1j * q)) * np.exp(-1j * q * (x1 - x2)))
        v = torch.zeros(36, dtype=c128)
        v[rows] = torch.tensor(psi / np.linalg.norm(psi))
        E = 9 * Delta / 4 + 2 * (c * np.cos(q) - Delta)                 # E_F + J (cos k1 + cos k2 - 2 Delta)
        assert abs(E.imag) < 1e-12 and close(H2 @ v, E.real * v, 1e-9)   # Proposition 9.12
        E_bethe.append(E.real)
        vecs.append(v)
assert torch.linalg.matrix_rank(torch.stack(vecs)) == 36             # a complete eigenbasis
assert close(torch.sort(torch.tensor(E_bethe)).values, torch.linalg.eigvalsh(H2))
print(len(E_bethe), n_bound)
ring8, Phi = SpinChain(8, 'periodic'), torch.zeros(28, dtype=c128)     # even L: Phi = sum_x (-1)^x |x, x+1>
for x in range(8):
    Phi[ring8.index_of[(1 << x) | (1 << (x + 1) % 8)]] = (-1) ** x
assert close(ring8.xxz(1.0, Delta, sector=2).matrix[0, 0] @ Phi, (8 * Delta / 4 - Delta) * Phi)   # E_F - J Delta
```

对 $L=9$，全部 $36$ 个本征态都来自 Bethe 方程：$34$ 个散射态和 $2$ 个束缚态。
完备性在这里是检验出来的，而不是证明出来的，并且对偶数 $L$ 它不成立。若 $\Delta\ne0$，$\Phi=\sum_x(-1)^x|x,x+1\rangle$
满足 $T\Phi=-\Phi$ 和 $H\Phi=(E_F-J\Delta)\Phi$（到 $|x-1,x+1\rangle$ 和 $|x,x+2\rangle$ 的跳跃成对抵消）。
Bethe 向量满足 $T\Psi=e^{-iK}\Psi$，因为 $\psi(x_1-1,x_2-1)=e^{-iK}\psi(x_1,x_2)$；若 $e^{iK}=-1$，则
$\cos k_1+\cos k_2=2\cos\frac K2\cos\frac{k_1-k_2}2=0$，其能量为 $E_F-2J\Delta$。所以 $\Phi$ 与所有
Bethe 向量正交（代码的最后几行对 $L=8$ 检验了 $H\Phi$）。取 $\Delta=0.6$ 和 $4\le L\le13$，上面的
搜索对奇数 $L$ 找到所有态，对偶数 $L$ 找到除 $\Phi$ 以外的所有态。对 $N>2$ 个磁振子，拟设有 $N!$ 项和 $N$ 个耦合的 Bethe 方程（H. Bethe, *Z. Phys.* 71 (1931)
205；M. Karbach, G. Müller, *Computers in Physics* 11 (1997) 36）。

> **物理意义。** 磁振子是自旋波：一个翻转的自旋以动量 $k$ 离域在整个环上。
> 两个磁振子发生带相移的散射（实的 $q$）或形成束缚态（复的 $q$）。Bethe 的洞见在于，
> 对 $N$ 个磁振子，散射分解为两体相位之积；例如，他的方程给出无限长 Heisenberg 链
> 每个格点的基态能量 $\frac14-\log2$（L. Hulthén, *Ark. Mat. Astron.
> Fys.* 26A (1938) no. 11；不加证明地陈述）。

## 9.5 Lanczos 与 Krylov 方法

扇区 $N=L/2$ 的维数在 $L=16$ 时为 $12870$，在 $L=20$ 时为 $184756$，但 $H$ 是稀疏的（每行 $O(L)$ 个非零
元素，`xxz_sparse`），而 Krylov 方法只需要乘积 $Hv$。设 $H$ 是 $n$ 阶厄米矩阵，
$v\ne0$。**Krylov 空间**为 $\mathcal K_m=\operatorname{span}\{v,Hv,\dots,H^{m-1}v\}$。

**命题 9.13（Lanczos）。** 设 $\dim\mathcal K_m=m$。令 $q_1=v/\|v\|$，并对 $j<m$ 令 $q_{j+1}$ 为
$Hq_j$ 正交于 $\mathcal K_j$ 的分量归一化后的向量。则 $q_1,\dots,q_m$ 是
$\mathcal K_m$ 的一组规范正交基；记 $Q_m=(q_1\cdots q_m)$，则 $T_m=Q_m^*HQ_m$ 是实对称三对角矩阵，其对角元为
$\alpha_j=\langle q_j,Hq_j\rangle$，非对角元为 $\beta_j=\langle q_{j+1},Hq_j\rangle>0$；并且 $HQ_m=Q_mT_m+r_me_m^T$，其中
$r_m\perp\mathcal K_m$ 为零当且仅当 $\mathcal K_{m+1}=\mathcal K_m$（否则 $r_m=\beta_mq_{m+1}$，$q_{m+1}=r_m/\|r_m\|$）。

*证明。* 设 $\operatorname{span}\{q_1,\dots,q_j\}=\mathcal K_j$（对 $j=1$ 成立）。则 $H\mathcal K_j\subseteq\mathcal K_j+\mathbb CHq_j$，
因为 $Hq_i\in\mathcal K_{i+1}$；所以 $\mathcal K_{j+1}=\mathcal K_j+\mathbb CH^jv\subseteq\mathcal K_j+\mathbb CHq_j\subseteq\mathcal K_{j+1}$，从而
$\mathcal K_{j+1}=\mathcal K_j+\mathbb CHq_j$。对 $j<m$ 维数增长，所以 $Hq_j$ 正交于
$\mathcal K_j$ 的分量非零，归纳得以继续；对 $j=m$，该分量 $r_m$ 为零当且仅当
$\mathcal K_{m+1}=\mathcal K_m$。对 $i>j+1$，$q_i\perp\mathcal K_{j+1}\ni Hq_j$，所以 $(T_m)_{ij}=0$，又由于 $T_m$ 是
厄米的，$(T_m)_{ji}=0$。$\alpha_j$ 是实的，而 $\beta_j=\langle q_{j+1},Hq_j\rangle$ 是该分量的范数，$>0$。所以
对 $j<m$ 有 $Hq_j=\beta_{j-1}q_{j-1}+\alpha_jq_j+\beta_jq_{j+1}$，且 $Hq_m=\beta_{m-1}q_{m-1}+\alpha_mq_m+r_m$。$\square$

$T_m$ 的特征值 $\theta^{(m)}_1\le\dots\le\theta^{(m)}_m$ 称为 **Ritz 值**。

**定理 9.14（单调收敛）。** 设 $\lambda_1\le\dots\le\lambda_n$ 为 $H$ 的特征值，且
$\dim\mathcal K_m=m$。则 $\theta^{(m)}_j\ge\lambda_j$，并且若 $\dim\mathcal K_{m+1}=m+1$，则 $\theta^{(m)}_j\ge\theta^{(m+1)}_j$
（$j\le m$）。若 $v$ 在 $\lambda_1$ 的特征空间中有非零分量且 $\mathcal K_{m+1}=\mathcal K_m$，则
$\theta^{(m)}_1=\lambda_1$。

*证明。* 由 Courant–Fischer（标准的线性代数），$\theta^{(m)}_j=\min_{U}\max_{x\in U,\|x\|=1}\langle x,Hx\rangle$，
其中 $U$ 取遍 $\mathcal K_m$ 的 $j$ 维子空间，因为 $Q_m$ 是到 $\mathcal K_m$ 上的等距，且
$\langle y,T_my\rangle=\langle Q_my,HQ_my\rangle$。同一公式对 $U\subseteq\mathcal K_{m+1}$ 给出 $\theta^{(m+1)}_j$，
对 $U\subseteq\mathbb C^n$ 给出 $\lambda_j$；在更多子空间上取最小给出更小的值。若
$\mathcal K_{m+1}=\mathcal K_m$，则 $\mathcal K_m$ 是 $H$-不变的，并对每个多项式 $p$ 含有 $p(H)v$，
特别地含有 $P_1v\ne0$（$P_1$ 为 $\lambda_1$ 的谱投影，由 Lagrange
插值它是 $H$ 的多项式，如定理 1.5 中那样）；所以在 $x=P_1v/\|P_1v\|\in\mathcal K_m$ 处 $\theta^{(m)}_1\le\langle x,Hx\rangle=\lambda_1$。$\square$

$\theta^{(m)}_1\to\lambda_1$ 的快慢取决于能隙 $\lambda_2-\lambda_1$ 与谱的展宽之比
（借助 Chebyshev 多项式的 Kaniel–Paige–Saad 界；Y. Saad, *Numerical Methods for Large Eigenvalue
Problems*, SIAM 2011, Ch. 6；不加证明地陈述）。在浮点运算中 $q_j$ 会失去正交性；`krylov.lanczos`
进行完全再正交化。对 XX 链，精确的基态能量是已知的（定理 9.10）：

```python
L, N = 12, 6
xx = SpinChain(L, complex_valued=False)
Hs = xx.xxz_sparse(J=1.0, Delta=0.0, sector=N)                      # sparse: diagonal + hopping
alpha, beta, Q = krylov.lanczos(Hs.matvec, torch.randn(1, xx.sector_dim(N), dtype=f64), 40)
T = torch.diag(alpha[0]) + torch.diag(beta[0], 1) + torch.diag(beta[0], -1)
assert close(Q[0] @ Q[0].T, torch.eye(40, dtype=f64)) and close(Q[0] @ Hs.to_dense()[0] @ Q[0].T, T)
ritz = torch.tensor([torch.linalg.eigvalsh(T[:m, :m])[0].item() for m in range(1, 41)])   # T_m = leading block
E_xx = sum(sorted(math.cos(math.pi * k / (L + 1)) for k in range(1, L + 1))[:N])
assert (ritz[1:] <= ritz[:-1] + 1e-12).all() and (ritz >= E_xx - 1e-12).all()             # Theorem 9.14
print([f"{(ritz[m - 1] - E_xx).item():.1e}" for m in (5, 10, 20, 40)])
```

最低 Ritz 值的误差单调下降，在 $5$、$10$ 和 $20$ 步之后分别为 $0.42$、$3.1\cdot10^{-3}$ 和 $3.1\cdot10^{-6}$，
在 $40$ 步之后达到舍入误差水平，而维数为 $924$。`krylov.ground_state`
反复重启 Lanczos，直到残差 $\|H\psi-E\psi\|$ 低于 $10^{-10}$。我们用 `scripts/chains.py` 中列出的九位有效数字的
Heisenberg 环能量（配套论文的参考值）来检验它：对 $L\le12$ 用稠密方法，
对 $L=16$ 用 Lanczos：

```python
RING = {4: -2.0, 6: -2.802775638, 8: -3.651093409, 10: -4.515446354, 12: -5.387390917,
        14: -6.263549534, 16: -7.142296361}
for n in (4, 6, 8, 10, 12):
    e0 = torch.linalg.eigvalsh(SpinChain(n, 'periodic').heisenberg(sector=n // 2).matrix[0, 0])[0].item()
    assert abs(e0 - RING[n]) < 5e-10
    print(n, round(e0 / n, 6))
t0 = time.time()
E16, psi16, res = krylov.ground_state(SpinChain(16, 'periodic', complex_valued=False).xxz_sparse(sector=8),
                                      generator=torch.Generator().manual_seed(0))
print(f"L = 16: E0 = {E16.item():.10f}, residual {res.item():.1e}, {time.time() - t0:.2f} s; "
      f"infinite chain: {0.25 - math.log(2):.6f}")
assert abs(E16.item() - RING[16]) < 5e-10 and res.item() < 1e-9
```

对 $L=16$（维数 $12870$），Lanczos 在单个 CPU 线程上用时远少于一秒。每个格点的能量
$-0.5,-0.467129,-0.456387,-0.451545,-0.448949$（$L=4,\dots,12$）从下方趋近无限长链的值
$\frac14-\log2=-0.443147$。

**时间演化。** 对 $e^{-i\tau H}v$，我们投影到 $\mathcal K_m$ 上：$w_\tau=\|v\|\,Q_me^{-i\tau T_m}e_1$。

**命题 9.15.** 设 $\dim\mathcal K_m=m$。对每个 $\tau$，$\|w_\tau\|=\|v\|$ 且 $\langle w_\tau,Hw_\tau\rangle=\langle v,Hv\rangle$。若
$\mathcal K_m$ 是 $H$-不变的，则 $w_\tau=e^{-i\tau H}v$。

*证明。* $Q_m$ 是等距，$e^{-i\tau T_m}$ 是酉的（$T_m$ 是实对称的），这给出范数的结论。
由于 $Q_m^*HQ_m=T_m$ 与 $e^{-i\tau T_m}$ 可交换，
$\langle w_\tau,Hw_\tau\rangle=\|v\|^2\langle e^{-i\tau T_m}e_1,T_me^{-i\tau T_m}e_1\rangle=\|v\|^2\langle e_1,T_me_1\rangle=\langle v,Hv\rangle$。若
$\mathcal K_m$ 是不变的，则 $HQ_m=Q_mT_m$（命题 9.13，$r_m=0$），所以 $H^jQ_m=Q_mT_m^j$ 且 $e^{-i\tau H}Q_m=Q_me^{-i\tau T_m}$；将其作用于
$e_1$，并注意 $Q_me_1=v/\|v\|$。$\square$

所以每个 Krylov 步都保持范数和能量，不论精确与否；精度要求 $m$ 相对于
$\tau\|H\|$ 足够大，而 `krylov.evolve` 细分时间区间并作后验误差估计
（M. Hochbruck, C. Lubich, *SIAM J. Numer. Anal.* 34 (1997) 1911；不加证明地陈述）。我们在 Heisenberg 环上
对 Néel 态 $|{\uparrow\downarrow\uparrow\downarrow}\cdots\rangle$ 作淬火，跟踪交错磁化
$m_s=\frac2L\sum_i(-1)^i\langle s^z_i\rangle$，并对 $L=12$ 与稠密传播子作比较：

```python
times, m_s = [0.0, 0.5, 1.0, 2.0, 4.0], {}
stag = lambda c, N, out: [((p[0].abs() ** 2) @ (c.bits(N).to(f64) - 0.5) @ torch.tensor(
    [(-1.0) ** i for i in range(c.L)], dtype=f64) * 2 / c.L).item() for p in out]
for n in (12, 16):
    c = SpinChain(n, 'periodic')
    Hn = c.xxz_sparse(J=1.0, Delta=1.0, sector=n // 2)
    psi0 = c.vector_in_sector('10' * (n // 2))[None]
    out = krylov.evolve(Hn, psi0, times)                               # (times, batch, dim)
    energy = torch.stack([(p.conj() * Hn.matvec(p)).sum(-1).real for p in out])
    assert close(out.norm(dim=-1), 1.0) and close(energy, energy[0])   # Proposition 9.15
    m_s[n] = stag(c, n // 2, out)
    print(n, [round(x, 4) for x in m_s[n]])
    if n == 12:
        w, U = torch.linalg.eigh(c.xxz(1.0, 1.0, sector=6).matrix[0, 0])
        for t, p in zip(times, out):
            assert close(p[0], U @ (torch.exp(-1j * t * w) * (U.mH @ psi0[0])))
print(f"{max(abs(a - b) for a, b in zip(m_s[12][:4], m_s[16][:4])):.1e}")   # L = 12 vs 16 up to t = 2
```

$m_s$ 从 $1$ 经过 $0.770$ 和 $0.279$ 衰减到 $t=2$ 时的 $-0.183$，两种尺寸的结果相差在 $1.7\cdot10^{-5}$ 以内；在
$t=4$ 时它们不同（$0.026$ 对 $0.035$）：较小的环的有限尺寸效应显现出来。

## 9.6 无序：能级统计与纠缠

**随机场 Heisenberg 链**的外场 $h_i$ 相互独立，在 $[-W,W]$ 上均匀分布。它保持 $S^z$ 守恒
（命题 9.1），而在一般情形下不再有其他守恒量：外场破坏了平移、反射、自旋翻转和
$SU(2)$。所以能级统计按扇区 $N$ 分别计算（定理 5.15 和 5.4 节：混合扇区会伪造出
独立的能级）。对弱无序，能级像 GOE 中那样相互排斥，$\langle r\rangle\approx0.5307$（5.3 节）；
对强无序，本征态局域化（**多体局域化**），能级变得相互独立，
$\langle r\rangle=2\log2-1\approx0.3863$（命题 5.14）。我们在扇区 $N=L/2$ 的中间一半能级中计算 $\langle r\rangle$，
并给出关于独立实现的一个标准误差（每个实现给出一个平均比值）：

```python
gen = torch.Generator().manual_seed(1)
r_mean = {}
for n, R in ((8, 200), (10, 100)):
    c = SpinChain(n, 'periodic', complex_valued=False)
    k = c.sector_dim(n // 2)
    for W in (0.5, 1, 2, 4, 8, 16):
        E = torch.linalg.eigvalsh(c.random_field_heisenberg(W, R, sector=n // 2, generator=gen).matrix[:, 0])
        r = SpinChain.level_spacing_ratio(E[:, k // 4: 3 * k // 4])     # one value per realisation
        r_mean[n, W] = (r.mean().item(), (r.std() / math.sqrt(R)).item())
    print(n, '  '.join(f"W={W}: {m:.3f}+-{s:.3f}" for (n_, W), (m, s) in r_mean.items() if n_ == n))
assert abs(r_mean[10, 0.5][0] - 0.5307) < 0.015 and abs(r_mean[10, 16][0] - 0.3863) < 0.015
assert r_mean[10, 1][0] > r_mean[10, 2][0] > r_mean[10, 4][0] > r_mean[10, 8][0]
```

对 $L=10$，平均比值从 $0.530$（$W=0.5$）经过 $0.484$（$W=2$）和 $0.410$（$W=4$）变到 $0.384$
（$W=16$），标准误差为 $0.003$–$0.004$。对 $L=8$，$W=0.5$ 处的值（$0.507$）低于
$W=1$ 处的值（$0.525$）：在 $W=0$ 时链是可积的，并具有上述所有对称性，而要得到 GOE 的值，
无序须比能级间距大得多。这一交叉过渡很宽，并随 $L$ 漂移；有限尺寸
标度把它定在 $W\approx3.7$ 附近（D. Luitz, N. Laflorencie, F. Alet, *Phys. Rev. B* 91 (2015) 081103），而
局域相在 $L\to\infty$ 时是否存活仍有争议（J. Šuntajs, J. Bonča, T. Prosen, L. Vidmar,
*Phys. Rev. E* 102 (2020) 062144）。

**自由费米子的纠缠。** 对开放 XX 链的本征态，前 $n$ 个格点的纠缠熵
由一个 $n\times n$ 矩阵决定。

**定理 9.16（Peschel 公式）。** 设 $\Phi=b_{k_1}^*\cdots b_{k_N}^*|0\rangle$ 如定理 9.10 中所述，
$C_{ij}=\langle\Phi,c_i^*c_j\Phi\rangle$，$C_A=(C_{ij})_{i,j<n}$。则 $0\le C_A\le1$，并且格点
$A=\{0,\dots,n-1\}$ 的纠缠熵为（$\nu$ 取遍 $C_A$ 的特征值，计重数）
$$
S_A=\sum_\nu \eta(\nu),\qquad \eta(\nu)=-\nu\log\nu-(1-\nu)\log(1-\nu)\quad(0\log0=0).
$$

*证明。* 对实向量 $x$ 令 $b(x)=\sum_jx_jc_j$；则 $b(x)^*b(y)^*=-b(y)^*b(x)^*$，且
$\{b(x),b(y)^*\}=\langle x,y\rangle$。令 $U\subseteq\mathbb R^L$ 为各 $u_{k_a}$ 张成的空间，$P_A$、$P_B=1-P_A$ 为到
$\{0,\dots,n-1\}$ 及其补集上的坐标投影。在 $U$ 上对角化形式 $\langle P_Ax,P_Ay\rangle$，得到 $U$ 的一组规范正交基
$\phi_1,\dots,\phi_N$，满足 $\langle P_A\phi_a,P_A\phi_b\rangle=\nu_a\delta_{ab}$，因而
$\langle P_B\phi_a,P_B\phi_b\rangle=(1-\nu_a)\delta_{ab}$。由于 $(x_1,\dots,x_N)\mapsto b(x_1)^*\cdots b(x_N)^*|0\rangle$ 是
多重线性且交错的，$\Phi=\pm\prod_a\big(\sqrt{\nu_a}f_a^*+\sqrt{1-\nu_a}\,g_a^*\big)|0\rangle$，其中
$f_a=b(P_A\phi_a/\sqrt{\nu_a})$，$g_a=b(P_B\phi_a/\sqrt{1-\nu_a})$（系数为零的项略去）。
展开得 $\Phi=\sum_I\pm\prod_{a\in I}\sqrt{\nu_a}\prod_{a\notin I}\sqrt{1-\nu_a}\,F_I^*G_{I^c}^*|0\rangle$，求和遍及
子集 $I$，其中 $F_I^*$、$G_{I^c}^*$ 为有序乘积。对 $j\ge n$，$c_j^*$ 的弦包含整个
$A$：$c_j^*=(-1)^{N_A}\otimes\tilde c_j^*$，其中 $\tilde c_j$ 是格点 $\ge n$ 的 Jordan–Wigner 算子，而对
$i<n$，$c_i^*=\hat c_i^*\otimes1$。由于 $(-1)^{N_A}$ 在 $|0_A\rangle$ 上为 $1$，
$F_I^*G_{I^c}^*|0\rangle=\xi_I\otimes\chi_{I^c}$，其中 $\xi_I=\hat F_I^*|0_A\rangle$，$\chi_{I^c}=\tilde G_{I^c}^*|0_B\rangle$，由定理 9.10 的 CAR 计算，
它们都是规范正交族。这是一个 Schmidt 分解（引理 7.7），其
权重 $p_I=\prod_{a\in I}\nu_a\prod_{a\notin I}(1-\nu_a)$ 是独立 Bernoulli 变量的分布，所以
$S_A=-\sum_Ip_I\log p_I=\sum_a\eta(\nu_a)$。最后 $C_{ij}=\sum_a\phi_a(i)\phi_a(j)$
（把 $c_i^*c_j$ 用 $b(\phi_a)$ 展开），所以 $C_A=\sum_a(P_A\phi_a)(P_A\phi_a)^T$ 在
正交向量 $P_A\phi_a$ 上的特征值为 $\nu_a$，在其他地方为 $0$，而 $\eta(0)=0$。$\square$

`entanglement_entropy` 逐个扇区计算约化态（它与块的磁化可交换，
命题 7.8）及其熵。扇区 $N$ 的基态是由 $N$ 个最低模构成的 $\Phi$，
由于各 $\varepsilon_k$ 互不相同，它是单的：

```python
L, N = 12, 6
xx = SpinChain(L)
w, V = torch.linalg.eigh(xx.xxz(J=1.0, Delta=0.0, sector=N).matrix[0, 0])
rho = DensityMatrix(xx.sector_algebra(N), matrix=torch.outer(V[:, 0], V[:, 0].conj())[None, None])
eps_k, u = np.linalg.eigh(np.eye(L, k=1) / 2 + np.eye(L, k=-1) / 2)  # h = A_L / 2, ascending
C = u[:, :N] @ u[:, :N].T                                            # correlation matrix of the N lowest modes
bern = lambda nu: -(nu * np.log(nu) + (1 - nu) * np.log(1 - nu))
S_peschel = [bern(np.clip(np.linalg.eigvalsh(C[:n, :n]), 1e-15, 1 - 1e-15)).sum() for n in range(1, L)]
S_lib = [xx.entanglement_entropy(rho, n, sector=N).item() for n in range(1, L)]
assert np.allclose(S_lib, S_peschel, atol=1e-9)                       # Theorem 9.16
print([round(s, 4) for s in S_lib[:6]])
```

熵 $0.6931,0.5199,0.7456,0.6283,0.7719,0.6560$（$n=1,\dots,6$）随 $n$ 的奇偶性振荡，
并在固定 $n/L$ 时如 $\frac16\log L$ 增长，这是开链的共形场论规律（I. Peschel,
*J. Phys. A* 36 (2003) L205；P. Calabrese, J. Cardy, *J. Stat. Mech.* (2004) P06002；不加证明地陈述）。

## 9.7 可与之比较的工具

QuSpin（Python）在磁化、动量、宇称和自旋反演的扇区中构造自旋链哈密顿量，
并使用 SciPy 的稀疏本征求解器和矩阵指数；专门的精确对角化程序
（XDiag、SPINPACK、HΦ）针对大型稀疏问题作了优化，借助格点对称性和大型计算机，
这类程序已达到约 40–50 个自旋；ITensor、TeNPy 和 TensorKit.jl（带有 $SU(2)$ 融合树）用
矩阵乘积态处理长链。本库所增加的东西不多而具体：把链作为
`TypeIAlgebra` $\bigoplus_NM_{\binom LN}$，使第 3–8 章的态、信道、熵和动力学
可以逐扇区应用；带有指标的包含塔；由最高权向量得到的 $SU(2)$ 块以及
Wigner–Eckart 关联；一次调用处理一批无序实现；以及从不离开扇区的约化态。
它没有反射或自旋翻转扇区，动量扇区仅适用于均匀 XXZ 链，稠密
方法到 $L\approx16$，Krylov 方法到 $L\approx24$（在大型 GPU 上约 $26$），并且没有张量网络。

## 习题

1. （手算动量扇区。）对 $L=4$ 和 $N=2$，求六个构型在 $T$ 下的轨道
   以及维数 $d_{2,m}$，$m=0,\dots,3$。
2. （求和规则。）设 $\psi$ 为 $L$ 为偶数的 Heisenberg 环的基态。证明
   $E_0=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle$ 且 $\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle=0$。
3. （对称性过多。）证明在扇区 $N=L/2$ 中，所有系数都等于 $1$ 的向量是
   $(S^-)^{L/2}|{\uparrow\cdots\uparrow}\rangle$ 的倍数，并且是 Heisenberg 环的能量为
   $JL/4$ 的特征向量；再证明每个 Krylov 空间 $\mathcal K_m(H,v)$ 都包含于每个含有 $v$ 的 $H$-不变子空间中。
   以这个向量为起点时，`krylov.ground_state` 返回什么？

<details>
<summary>解答</summary>

1. $0011\to1001\to1100\to0110\to0011$ 是一条 $R=4$ 的轨道，$0101\to1010\to0101$ 是一条 $R=2$ 的轨道。由
   命题 9.5，前者对每个 $m$ 都有贡献，后者对满足 $2m\equiv0\pmod4$ 即 $m=0,2$ 的 $m$ 有贡献：
   $d_{2,m}=2,1,2,1$。

```python
assert SpinChain(4, 'periodic').momentum_dims(2) == [2, 1, 2, 1]
```

2. 基态是单的单态（定理 9.9），且 $T$ 与 $H$ 可交换，所以 $T\psi=e^{ik}\psi$。由于
   $\mathbf s_i\cdot\mathbf s_{i+1}=T^i(\mathbf s_0\cdot\mathbf s_1)T^{*i}$，$\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle$ 对全部 $L$ 条键都相同，从而
   $E_0=J\sum_i\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle$。由于 $\mathbf S^2\psi=0$，
   再次利用平移不变性，$0=\langle\mathbf S^2\rangle=\sum_{i,j}\langle\mathbf s_i\cdot\mathbf s_j\rangle=L\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle$。

```python
assert abs(12 * corr_we[1].item() - gs.energy.item()) < 1e-10 and abs(corr_we.sum().item()) < 1e-10
```

3. $S^-$ 降低其中一个向上自旋；把它作用于完全极化态 $L/2$ 次，会以 $(L/2)!$ 种次序到达每个
   含 $L/2$ 个向上自旋的构型，每种次序的系数都是 $1$。$|{\uparrow\cdots\uparrow}\rangle$ 的全部
   $L$ 条键都是平行的，能量为 $JL/4$，而 $S^-$ 与 $H$ 可交换，所以该向量是具有相同
   能量的特征向量。若 $W$ 是不变的且 $v\in W$，则对所有 $j$ 有 $H^jv\in W$。这里对所有 $m$ 有 $\mathcal K_m=\mathbb Cv$，
   `ground_state` 返回 $E=JL/4$，即谱的顶端（每个 $\mathbf s_i\cdot\mathbf s_j\le\frac14$）；随机的初始
   向量以概率一避免这种情况。

```python
c10 = SpinChain(10, 'periodic', complex_valued=False)
E_sym, _, _ = krylov.ground_state(c10.xxz_sparse(1.0, 1.0, sector=5), v0=torch.ones(1, 252, dtype=f64))
assert abs(E_sym.item() - 2.5) < 1e-10 and E_sym.item() > RING[10] + 7     # far from E_0
```

</details>

**下一章：** [第 10 章](10_quantum_computing.md)用扇区讨论量子计算：无退相干
子空间与无噪声子系统、算子代数纠错，以及保持守恒的 Hamming 重量的混合算子——
本章的代数 $\bigoplus_wM_{\binom nw}$ 以新的角色出现。
