[English](../13_jordan_algebras.md) · [Español](../es/13_jordan_algebras.md) · [Français](../fr/13_jordan_algebras.md) · [Deutsch](../de/13_jordan_algebras.md) · **中文** · [日本語](../ja/13_jordan_algebras.md) · [Русский](../ru/13_jordan_algebras.md)

# 13. 若尔当代数：量子力学为何使用 *-代数

两个可观测量的乘积 $xy$ 不是可观测量，除非 $x$ 与 $y$ 对易；但对称化乘积
$x\circ y=\frac12(xy+yx)$ 总是可观测量。1933 年 P. Jordan 提议把带有这一乘积的可观测量作为量子力学的
基本结构，随后 Jordan、von Neumann 与 Wigner 对由此产生的有限维代数进行了分类。本章在证明简短之处给出
完整证明，发展有限维理论，并追问若尔当结构遗漏了什么。答案是对易子的符号：单个系统的若尔当代数只在
相差一个符号的意义下决定对易子（13.5 节），但复合系统需要这个符号，而且它决定一个态是否纠缠。

**你将学到**若尔当积与二次表示 $P(y)x=yxy$（它与 $\circ$ 不同，保持正性）；自旋因子、它们的谱理论以及
洛伦兹锥；不是 $*$-代数自伴部分的厄米矩阵若尔当子代数；Jordan–von Neumann–Wigner 分类；为什么转置是
若尔当自同构而不是 $*$-自同构，以及它在复合系统上起什么作用；以及第 11 章中的一个迹不等式，用若尔当语言
写出并在自旋因子中加以证明。

在本章中，$\operatorname{Herm}_n(\mathbb K)$ 是 $\mathbb K=\mathbb R$ 或 $\mathbb C$ 上 $n\times n$ 厄米矩阵构成的
实向量空间，$M_{sa}$ 是 $\mathcal M=\bigoplus_cM_{k_c}(\mathbb C)$ 的自伴部分，
$X,Y,Z$ 是 Pauli 矩阵，$\sigma=(X,Y,Z)$，于是对 $u\in\mathbb R^3$ 有 $u\cdot\sigma=u_1X+u_2Y+u_3Z$。
$[a,b]=ab-ba$，$|u|$ 是 $u\in\mathbb R^m$ 的欧几里得范数。

## 13.1 可观测量与若尔当积

**定义。** **若尔当代数**是一个实向量空间 $J$，带有双线性、交换的乘积 $x\circ y$，
满足**若尔当恒等式** $(x^2\circ y)\circ x=x^2\circ(y\circ x)$，其中 $x^2=x\circ x$。若它同构于某个结合代数中
对 $\frac12(xy+yx)$ 封闭的子空间，则称它是**特殊的**。在 $J$ 上，令
$L(y)x=y\circ x$，**二次表示**为 $P(y)=2L(y)^2-L(y^2)$。

**引理 13.1.** 设 $\mathcal A$ 是实或复的结合代数，$x\circ y=\frac12(xy+yx)$。则
$\circ$ 满足若尔当恒等式，且若尔当幂 $x^{\circ1}=x$，$x^{\circ(k+1)}=x\circ x^{\circ k}$ 等于
结合幂 $x^k$。特别地，$M_{sa}$、$\operatorname{Herm}_n(\mathbb C)$ 与
$\operatorname{Herm}_n(\mathbb R)$ 都是特殊若尔当代数。

*证明。* 展开得 $4(x^2\circ y)\circ x=x^2yx+yx^3+x^3y+xyx^2$ 与
$4x^2\circ(y\circ x)=x^2yx+x^3y+yx^3+xyx^2$，是同样的四项。由归纳法，
$x\circ x^k=\frac12(x^{k+1}+x^{k+1})=x^{k+1}$。若 $x,y$ 是厄米的，则 $(xy+yx)^*=yx+xy$，所以这三个空间
对 $\circ$ 封闭；对实对称的 $x,y$，乘积是实的。$\square$

**命题 13.2.** 设 $x,y,z\in M_{sa}$。

1. $P(y)x=yxy$，且**基本公式** $P(P(y)x)=P(y)P(x)P(y)$ 成立。
2. $\operatorname{Tr}P(y)x=\operatorname{Tr}y^2x$。
3. 若 $x\ge0$，则 $P(y)x\ge0$。若 $y$ 可逆，则 $P(y)$ 是正锥到自身的双射，其
   逆为 $P(y^{-1})$。
4. 正元素的若尔当积不一定是正的：对 $p=\begin{pmatrix}1&0\\0&0\end{pmatrix}$ 与
   $q=\frac12\begin{pmatrix}1&1\\1&1\end{pmatrix}$，$\det(p\circ q)=-\frac1{16}$。

*证明。* (1) $2y\circ(y\circ x)=\frac12(y^2x+2yxy+xy^2)$，$y^2\circ x=\frac12(y^2x+xy^2)$；两式相减。于是由结合律，
$P(P(y)x)z=yxy\,z\,yxy=P(y)P(x)P(y)z$。(2) 迹的循环性。(3)
$yxy=(x^{1/2}y)^*(x^{1/2}y)\ge0$（推论 1.7），且 $P(y)P(y^{-1})=P(y^{-1})P(y)=\mathrm{id}$。(4)
$p\circ q=\frac14\begin{pmatrix}2&1\\1&0\end{pmatrix}$。$\square$

将 $P(y)$ 写成 $2L(y)^2-L(y^2)$ 后，(1) 中的每个恒等式都只用 $\circ$ 表达，所以它们在每个特殊若尔当代数中成立；基本公式在每个若尔当代数中
都成立（K. McCrimmon, *A Taste of Jordan Algebras*, Springer 2004，我们不在此
重复其证明）。模块 `torch_vn_algebra.jordan` 在 `TypeIAlgebra` 的算子上实现这些映射，惰性求值，
并在所有扇区上批量进行：

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, SpinFactor
from torch_vn_algebra import jordan as jd

torch.manual_seed(0)
dt = torch.complex128
jp, P = jd.jordan_product, jd.quadratic_representation

def close(A, B, tol=1e-10):
    A = A.matrix if hasattr(A, 'matrix') else A
    B = B.matrix if hasattr(B, 'matrix') else B
    return (A - B).abs().max().item() < tol

alg = TypeIAlgebra([2, 3], [2, 3], precision='double')           # M_2 (+) M_3
def sa(B=4):
    return alg.operator_from_eigenvalues(lambda d: torch.randn(B, d), batch_size=B, force_self_adjoint=True)
x, y, z = sa(), sa(), sa()
x2 = jp(x, x)
assert close(jp(jp(x2, y), x), jp(x2, jp(y, x)))                  # Jordan identity (Lemma 13.1)
assert close(jd.jordan_power(x, 5), x @ x @ x @ x @ x)            # power associativity
assert close(P(y, x), 2 * jp(y, jp(y, x)) - jp(jp(y, y), x))      # P(y) = 2L(y)^2 - L(y^2)
assert close(P(P(y, x), z), P(y, P(x, P(y, z))))                  # fundamental formula
assert close(P(y, x).trace, (jp(y, y) @ x).trace)                 # Tr P(y)x = Tr y^2 x
pos = alg.operator_from_eigenvalues(lambda d: torch.rand(4, d), batch_size=4, force_positive=True)
assert all((e > -1e-12).all() for e in P(y, pos).eigenvalues())  # Proposition 13.2(3)

M2 = TypeIAlgebra([2], [2], precision='double')
p = M2.operator(torch.tensor([[[[1.0, 0], [0, 0]]]]), is_self_adjoint=True, is_positive=True)
q = M2.operator(torch.tensor([[[[0.5, 0.5], [0.5, 0.5]]]]), is_self_adjoint=True, is_positive=True)
print(jp(p, q).eigenvalues()[0])                                  # (1 -+ sqrt 2)/4: not positive
assert abs(torch.linalg.det(jp(p, q).matrix[0, 0]).real.item() + 1 / 16) < 1e-14
```

> **物理意义。** $x\circ y$ 是可观测量的对称化乘积：它是唯一一个对 $x,y$ 对称、且在 $y=x$ 时化为 $x^2$ 的
> 组合 $\alpha xy+\beta yx$，并且它仍是可观测量。二次表示就是 **Lüders 规则**：测量一个是/否问题 $p$ 并保留
> 结果"是"，会把密度 $\rho$ 变为 $P(p)\rho/\operatorname{Tr}P(p)\rho=p\rho p/\operatorname{Tr}(p\rho)$。命题 13.2(4)
> 表明，$\circ$ 本身对正元素没有这样的操作含义。

## 13.2 自旋因子与洛伦兹锥

**定义。** 对 $m\ge1$，**自旋因子** $V_m=\mathbb R\oplus\mathbb R^m$ 的元素为 $x=(a,u)$，
乘积为 $(a,u)\circ(b,v)=(ab+u\cdot v,\;av+bu)$，单位为 $1=(1,0)$。它的**迹**为 $\operatorname{tr}x=2a$。

这一乘积是双线性且交换的；若尔当恒等式可由一个矩阵表示得出。

**引理 13.3（反对易的对称元）。** 设 $p=\lfloor m/2\rfloor$，$D=2^p$（若 $m=1$ 则 $D=2$）。在
$(\mathbb C^2)^{\otimes p}$ 上，对 $j\le p$ 令 $s_{2j-1}=Z^{\otimes(j-1)}\otimes X\otimes1$，$s_{2j}=Z^{\otimes(j-1)}\otimes Y\otimes1$，
若 $m$ 为奇数，再令 $s_m=Z^{\otimes p}$（若 $m=1$ 则 $s_1=Z$）。则 $s_i=s_i^*$ 且
$s_is_j+s_js_i=2\delta_{ij}1$。映射 $\phi(a,u)=a1+\sum_iu_is_i$ 是单射，并满足
$\phi(x\circ y)=\phi(x)\circ\phi(y)$；因此 $V_m$ 是特殊若尔当代数。对 $m\ge2$，诸 $s_i$ 作为 $*$-代数生成
$M_D(\mathbb C)$。

*证明。* Pauli 矩阵的平方为 $1$ 且两两反对易；$Z$ 与 $X$ 和 $Y$ 也反对易。两个不同的串
$s_i,s_j$ 恰好在一个格点上以反对易的方式不同：对 $s_{2j-1},s_{2j}$，是在格点 $j$（$X$ 对
$Y$）；对指标为 $j<j'$ 的串，是在格点 $j$（$X$ 或 $Y$ 对 $Z$），而在其他格点上它们带有相同的
因子或因子 $1$；对 $s_m=Z^{\otimes p}$ 也是如此。因此对 $i\ne j$ 有 $s_is_j=-s_js_i$，且
$s_i^2=1$。于是 $(u\cdot s)\circ(v\cdot s)=\sum_{ij}u_iv_j\,s_i\circ s_j=(u\cdot v)1$，这给出
$\phi(x)\circ\phi(y)=\phi(x\circ y)$。每个 $s_i$ 都是无迹的（它是带有一个无迹因子的张量积），且
$\operatorname{Tr}(s_is_j)=\operatorname{Tr}(s_i\circ s_j)=D\delta_{ij}$，所以 $1,s_1,\dots,s_m$ 关于
$\operatorname{Tr}(ab)$ 正交，$\phi$ 是单射。若尔当恒等式在 $\operatorname{Herm}_D$ 中成立（引理 13.1），并经由 $\phi$
传回。最后 $s_{2j-1}s_{2j}=1\otimes\cdots\otimes XY\otimes\cdots\otimes1=iZ_j$（$Z$ 位于格点
$j$），所以生成的 $*$-代数包含 $Z_j$，进而包含 $X_j=Z_1\cdots Z_{j-1}s_{2j-1}$ 与 $Y_j$，从而包含 Pauli 矩阵的所有
张量积，它们构成 $M_D(\mathbb C)$ 的一组基。$\square$

在相差符号的意义下，$s_{2j-1}$ 与 $s_{2j}$ 就是 9.4 节中 Jordan–Wigner 费米子 $c_j$ 的
Majorana 算子 $c_j+c_j^*$ 与 $i(c_j-c_j^*)$。

**定理 13.4（$V_m$ 的谱理论）。** 设 $x=(a,u)\in V_m$，$\lambda_\pm=a\pm|u|$。

1. $V_m$ 的幂等元是 $0$、$1$ 以及 $|w|=1$ 的 $c_w=\frac12(1,w)$；且 $c_w\circ c_{w'}=0$ 当且仅当 $w'=-w$。
2. 若 $u\ne0$，则 $x=\lambda_+e_++\lambda_-e_-$，其中 $e_\pm=c_{\pm u/|u|}$；若 $u=0$，则 $x=a1$。若 $x=\mu_1c_w+\mu_2c_{-w}$，
   则 $\{\mu_1,\mu_2\}=\{\lambda_+,\lambda_-\}$。我们称 $\lambda_\pm$ 为 $x$ 的**特征值**。
3. 对 $f:\{\lambda_+,\lambda_-\}\to\mathbb R$，令 $f(x)=f(\lambda_+)e_++f(\lambda_-)e_-$（若 $u=0$ 则为 $f(a)1$）。则
   $x^{\circ k}=\lambda_+^ke_++\lambda_-^ke_-$，$\phi(f(x))=f(\phi(x))$，$\phi(x)$ 的特征值为 $\lambda_\pm$，
   重数各为 $D/2$，且 $\operatorname{Tr}\phi(x)=\frac D2\operatorname{tr}x$。此外
   $\operatorname{tr}x=\lambda_++\lambda_-$，$\det x:=\lambda_+\lambda_-=a^2-|u|^2$。
4. 平方元构成的锥 $Q=\{x^2:x\in V_m\}$ 就是**洛伦兹锥** $\{(a,u):a\ge|u|\}=\{x:\lambda_-\ge0\}=\{x:\phi(x)\ge0\}$。
   它是**自对偶的**：$y\in Q$ 当且仅当对所有 $x\in Q$ 有 $\operatorname{tr}(x\circ y)\ge0$。
5. $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$；因此对每个 $y\in V_m$ 有 $P(y)Q\subseteq Q$。

*证明。* (1) $(a,u)^2=(a^2+|u|^2,2au)$ 等于 $(a,u)$ 当且仅当 $u=0$、$a\in\{0,1\}$，或 $a=\frac12$、$|u|=\frac12$。
$c_w\circ c_{w'}=\frac14(1+w\cdot w',\,w+w')$ 为零当且仅当 $w'=-w$。(2)
$\lambda_+e_++\lambda_-e_-=(a,0)+|u|(0,u/|u|)=x$，且 $\mu_1c_w+\mu_2c_{-w}=(\frac{\mu_1+\mu_2}2,\frac{\mu_1-\mu_2}2w)$，所以
$\mu_1+\mu_2=2a$，$|\mu_1-\mu_2|=2|u|$。(3) 由 (1)，$(\alpha e_++\beta e_-)\circ(\gamma e_++\delta e_-)=\alpha\gamma
e_++\beta\delta e_-$，由归纳法得到幂的公式。令 $w=u/|u|$，则 $(w\cdot s)^2=1$，所以
$\phi(e_\pm)=\frac12(1\pm w\cdot s)$ 是和为 $1$、乘积为 $0$ 的正交投影，迹各为 $D/2$；因此
$\phi(x)=\lambda_+\phi(e_+)+\lambda_-\phi(e_-)$ 是一个谱分解，且
$f(\phi(x))=f(\lambda_+)\phi(e_+)+f(\lambda_-)\phi(e_-)=\phi(f(x))$（定理 1.5）。(4) $x^2$ 的特征值为
$\lambda_\pm^2\ge0$；反之，若 $\lambda_-\ge0$，则 $x=(\sqrt x)^2$。$\lambda_-\ge0$ 当且仅当 $a\ge|u|$，由 (3) 又当且仅当
$\phi(x)\ge0$。关于自对偶性，$\operatorname{tr}(x\circ y)=2(ab+u\cdot v)$。若 $a\ge|u|$ 且 $b\ge|v|$，则
$ab\ge|u||v|\ge-u\cdot v$。若 $b<|v|$，则 $x=(|v|,-v)\in Q$ 给出 $2(|v|b-|v|^2)<0$（若 $v=0$、$b<0$，取 $x=1$）。
(5) $P(y)x=2y\circ(y\circ x)-y^2\circ x$，且 $\phi$ 保持 $\circ$，所以由命题 13.2(1) 得
$\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$；再应用 (4) 与命题 13.2(3)。$\square$

**命题 13.5.** (a) 对 $m=3$，$\phi(a,u)=a1+u\cdot\sigma$ 是若尔当同构
$V_3\cong\operatorname{Herm}_2(\mathbb C)$。(b) $V_m$ 的若尔当自同构恰为映射 $(a,u)\mapsto(a,Ou)$，
$O\in O(m)$。(c) 对 $m\ge2$，$\phi(V_m)$ 的维数为 $m+1$，而它生成的 $*$-代数的自伴部分
维数为 $D^2=4^{\lfloor m/2\rfloor}$；两者仅在 $m=3$ 时相等。

*证明。* (a) 对 $m=3$，$(s_1,s_2,s_3)=(X,Y,Z)$，且 $\phi$ 是两个 $4$ 维空间之间的单射。(b) 映射
$(a,u)\mapsto(a,Ou)$ 保持 $u\cdot v$，从而保持 $\circ$。反之，设 $\alpha$ 是若尔当自同构。则
$\alpha(1)$ 是 $V_m$ 的单位，所以 $\alpha(1)=1$。由 $x^2$ 的公式，$x^2=1$ 当且仅当 $x=\pm1$ 或 $x=(0,w)$，其中
$|w|=1$。对 $|w|=1$，$\alpha(0,w)^2=1$ 且 $\alpha(0,w)\ne\pm1=\alpha(\pm1)$，所以 $\alpha(0,w)=(0,w')$，$|w'|=1$。
把 $\alpha(0,u)$ 写成 $(\ell(u),Ou)$，其中 $\ell,O$ 线性，上述结论即 $\ell=0$ 且 $|Ou|=|u|$，即 $O\in O(m)$。
(c) 引理 13.3；对 $m=2$ 与 $m\ge4$ 有 $4^{\lfloor m/2\rfloor}>m+1$。$\square$

```python
V = SpinFactor(5)
xs, ys = V.random(1000), V.random(1000)                          # (1000, 6): batches of (a, u)
x2 = V.square(xs)
assert (V.product(V.product(x2, ys), xs) - V.product(x2, V.product(ys, xs))).abs().max() < 1e-12
lam, e = V.spectral_decomposition(xs)                            # (lambda-, lambda+), (e-, e+)
assert (lam[:, 1, None] * e[:, 1] + lam[:, 0, None] * e[:, 0] - xs).abs().max() < 1e-12
Xm, Ym, Dm = V.to_matrix(xs), V.to_matrix(ys), V.representation_dim  # phi into Herm_4(C)
assert (V.to_matrix(V.product(xs, ys)) - 0.5 * (Xm @ Ym + Ym @ Xm)).abs().max() < 1e-12
assert (torch.linalg.eigvalsh(Xm) - lam.repeat_interleave(Dm // 2, dim=-1)).abs().max() < 1e-12
assert (V.to_matrix(V.exp(xs)) - torch.linalg.matrix_exp(Xm)).abs().max() < 1e-9       # phi(f(x)) = f(phi(x))
cone = V.random(1000, cone=True)                                 # squares
assert V.in_cone(cone).all() and (torch.linalg.eigvalsh(V.to_matrix(cone)) > -1e-12).all()
assert (V.square(V.sqrt(cone)) - cone).abs().max() < 1e-10
O, _ = torch.linalg.qr(torch.randn(5, 5, dtype=torch.float64))  # Proposition 13.5(b)
rot = lambda w: torch.cat([w[..., :1], w[..., 1:] @ O.T], dim=-1)
assert (rot(V.product(xs, ys)) - V.product(rot(xs), rot(ys))).abs().max() < 1e-12

for m in range(1, 7):                                            # Proposition 13.5(c)
    r = jd.compare_with_star_algebra(jd.anticommuting_symmetries(m))
    print(f"m = {m}: dim V_m = {r['jordan_dim']}, dim of the self-adjoint part of the *-algebra = {r['star_dim']}")
    assert r['jordan_dim'] == m + 1 and r['equal'] == (m in (1, 3))
```

对 $m=1$，$s_1=Z$ 生成对角矩阵，维数为 $2=m+1$，而 $V_1\cong\mathbb R\oplus\mathbb R$ 不是
单的。打印出的维数是：$m=2$ 时为 $3,4$；$m=4$ 时为 $5,16$；$m=6$ 时为 $7,64$。

> **物理意义。** 由自对偶性，$V_m$ 的态（正的、归一化的线性泛函）是
> $x\mapsto\operatorname{tr}(\rho\circ x)$，其中 $\rho=\frac12(1,r)\in Q$，$|r|\le1$：态空间是一个
> $m$ 维的球。对 $m=3$，它就是一个量子比特的 Bloch 球（命题 13.5(a)）。洛伦兹锥是 $1+m$ 维狭义相对论中的
> 未来光锥，而命题 13.5(b) 说，若尔当自同构就是固定时间轴的空间正交变换。

锥 $Q$ 是**二阶锥规划**（SOCP）的约束集：在若干洛伦兹锥的乘积上的线性优化。SOCP 的
内点法用到的恰恰是定理 13.4：特征值、$\det$，以及用于缩放的 $P(y)$
（F. Alizadeh, D. Goldfarb, *Second-order cone programming*, Math. Program. 95 (2003) 3–51）。

## 13.3 不是 *-子代数的若尔当子代数

若尔当子代数 $J\subseteq\operatorname{Herm}_n(\mathbb C)$ 生成一个复 $*$-代数 $\mathcal A$，且
$J\subseteq\mathcal A_{sa}$。命题 13.5(c) 给出了 $J\ne\mathcal A_{sa}$ 的例子；这里再给出两个。

**命题 13.6.** (a) $\operatorname{Herm}_n(\mathbb R)$ 是 $\operatorname{Herm}_n(\mathbb C)$ 的若尔当子代数，
维数为 $n(n+1)/2$；它生成的 $*$-代数是 $M_n(\mathbb C)$，其自伴部分的维数为 $n^2$。
(b) 对 $x,y\in\operatorname{Herm}_2(\mathbb C)$，存在酉矩阵 $U$ 使得 $UxU^*$ 与 $UyU^*$ 都是实
对称的。由 $1,x,y$ 生成的若尔当代数是 $\operatorname{span}\{1,x,y\}$；若 $x,y,1$ 线性
无关，则它是 $U^*\operatorname{Herm}_2(\mathbb R)U\cong V_2$，维数为 $3$，而由
$x,y$ 生成的 $*$-代数是 $M_2(\mathbb C)$，维数为 $4$。

*证明。* (a) 引理 13.1。对 $i\ne j$，$E_{ij}=E_{ii}(E_{ij}+E_{ji})$ 是实对称矩阵的乘积，所以
生成的代数包含所有矩阵单位。(b) 写 $x=a1+u\cdot\sigma$，$y=b1+v\cdot\sigma$，并取单位向量
$n\perp u,v$。矩阵 $n\cdot\sigma$ 与 $Y$ 都是厄米的、无迹的，且平方为 $1$，所以两者的
特征值都是 $1,-1$，从而对某个酉矩阵 $U$ 有 $U(n\cdot\sigma)U^*=Y$（定理 1.3）。用 $U$ 共轭保持
$\operatorname{Tr}$ 与 $\circ$，所以它把 $u\cdot\sigma$ 映为 $(Ru)\cdot\sigma$，其中 $R$ 线性，且
$Ru\cdot Rv=\frac12\operatorname{Tr}\big((Ru\cdot\sigma)(Rv\cdot\sigma)\big)=u\cdot v$：$R\in O(3)$，且 $Rn=e_2$。
因此 $Ru,Rv\perp e_2$，而 $UxU^*=a1+(Ru)_1X+(Ru)_3Z$ 是实的，因为 $X,Z$ 是实的；$UyU^*$ 同理。其次，
$x\circ y=(ab+u\cdot v)1+av\cdot\sigma+bu\cdot\sigma\in\operatorname{span}\{1,x,y\}$，$x^2,y^2$ 也类似，所以该张成空间
是封闭的。若它的维数为 $3$，则它在 $U$ 下的共轭就是整个 $\operatorname{Herm}_2(\mathbb R)=\operatorname{span}\{1,X,Z\}$，
由引理 13.3（$X,Z$ 反对易）它 $\cong V_2$，且由 (a) 它生成 $M_2(\mathbb C)$。$\square$

```python
for n in range(2, 6):                                            # Proposition 13.6(a), a random real pair
    G = torch.randn(2, n, n, dtype=torch.float64)
    r = jd.compare_with_star_algebra(G + G.mT)
    assert (r['jordan_dim'], r['star_dim']) == (n * (n + 1) // 2, n * n)

H = torch.randn(3, 2, 2, dtype=dt)
H = H + H.mH                                                     # three random Hermitian 2 x 2 matrices
r2, r3 = jd.compare_with_star_algebra(H[:2]), jd.compare_with_star_algebra(H)
print("two generators:", r2['jordan_dim'], r2['star_dim'], "  three:", r3['jordan_dim'], r3['star_dim'])
assert (r2['jordan_dim'], r2['star_dim'], r3['jordan_dim']) == (3, 4, 4)

V3 = SpinFactor(3)
u, v = V3.from_matrix(H[0])[1:], V3.from_matrix(H[1])[1:]        # Bloch vectors
nrm = torch.linalg.cross(u, v)
nrm = nrm / nrm.norm()
_, W1 = torch.linalg.eigh(V3.to_matrix(torch.cat([torch.zeros(1, dtype=torch.float64), nrm])))
_, W2 = torch.linalg.eigh(jd.anticommuting_symmetries(3)[1])     # Y
Um = W2 @ W1.mH                                                  # U (n . sigma) U* = Y
assert all((Um @ A @ Um.mH).imag.abs().max() < 1e-12 for A in H[:2])
```

三个一般的 $2\times2$ 厄米矩阵生成整个 $\operatorname{Herm}_2(\mathbb C)$：两列中的维数都是 $4$。
因此在 $\mathcal M$ 的一个 $M_2$ 块中，两个可观测量永远不会揭示复结构；第三个则
可以。

## 13.4 分类

若尔当代数称为**形式实的**，如果 $\sum_ix_i^2=0$ 蕴含所有 $x_i=0$；在有限维中，这
等价于它是**欧几里得的**，即具有满足 $\langle x\circ y,z\rangle=\langle y,x\circ z\rangle$ 的内积
[FK94]。$M_{sa}$（取 $\operatorname{Tr}(xy)$）与 $V_m$（取 $\operatorname{tr}(x\circ y)$）都是
欧几里得的。

**定理 13.7（Jordan–von Neumann–Wigner；不加证明地引用）。** 每个有限维形式实若尔当
代数都是单若尔当代数的直和，而每个单的若尔当代数恰好同构于下列之一：$\mathbb R$；
$\mathbb K=\mathbb R,\mathbb C,\mathbb H$ 且 $r\ge3$ 时的 $\operatorname{Herm}_r(\mathbb K)$；自旋因子
$V_m$，$m\ge2$；**Albert 代数** $\operatorname{Herm}_3(\mathbb O)$，即八元数上 $3\times3$ 厄米矩阵构成的代数
（维数为 $27$）。在秩 $2$ 的情形，对
$\mathbb K=\mathbb R,\mathbb C,\mathbb H,\mathbb O$ 有 $\operatorname{Herm}_2(\mathbb K)\cong V_{\dim\mathbb K+1}$。（P. Jordan, J. von Neumann, E. Wigner, *On an algebraic
generalization of the quantum mechanical formalism*, Ann. of Math. 35 (1934) 29–64；现代的证明见 [FK94,
Chapter V]。）

除 Albert 代数外，以上所有代数都是特殊的：把四元数写成 $2\times2$ 复矩阵，$\operatorname{Herm}_r(\mathbb H)$ 就嵌入
$\operatorname{Herm}_{2r}(\mathbb C)$；$V_m$ 则由引理 13.3 嵌入。
Albert 代数是**例外的**：它不同构于任何结合代数的若尔当子代数
（A. A. Albert, *On a certain algebra of quantum mechanics*, Ann. of Math. 35 (1934) 65–73）。一个可以
简短陈述的理由是：C. M. Glennie 找到了一些多项式恒等式，它们在每个特殊若尔当代数中成立，但在
$\operatorname{Herm}_3(\mathbb O)$ 中不成立（*Some identities valid in special Jordan algebras but not valid in all Jordan
algebras*, Pacific J. Math. 16 (1966) 47–59；参见 McCrimmon, *A Taste of Jordan Algebras*）。八元数不满足
结合律，而 $\operatorname{Herm}_r(\mathbb O)$ 仅在 $r\le3$ 时满足若尔当恒等式，这就是
列表止于此处的原因。所以 $\operatorname{Herm}_3(\mathbb O)$ 不能用实或复矩阵（也不能在任何结合代数中）
以乘积 $\frac12(xy+yx)$ 忠实地表示；软件必须直接处理八元数矩阵元。

**定理 13.8（Koecher–Vinberg；不加证明地引用）。** 欧几里得若尔当代数平方锥的内部
恰好就是**对称锥**：欧几里得空间中自对偶且
齐性（其线性自同构群可迁地作用）的开凸锥。[FK94, Chapter III.]

定理 13.4(4) 与习题 1 对洛伦兹锥证明了这两条性质。因此半定规划的锥
（$\operatorname{Herm}_n(\mathbb R)^+$ 或 $\operatorname{Herm}_n(\mathbb C)^+$ 的乘积）、SOCP 的锥以及线性
规划的锥（$\mathbb R_+^n$，即 $\mathbb R^n=\mathbb R\oplus\dots\oplus\mathbb R$ 的锥）都具有这种形式，而
内点法可以对所有欧几里得若尔当代数一次性地写出。

[FK94]: J. Faraut, A. Korányi, *Analysis on Symmetric Cones*, Oxford University Press, 1994.

## 13.5 若尔当映射与 *-映射

冯·诺依曼代数之间的**若尔当 $*$-同态**是满足 $\Phi(x^*)=\Phi(x)^*$ 与
$\Phi(x\circ y)=\Phi(x)\circ\Phi(y)$ 的复线性映射 $\Phi$；在自伴部分上它是若尔当同态，而自伴部分之间的每个若尔当
同态都可以由 $\Phi(x+iy)=\Phi(x)+i\Phi(y)$ 延拓为若尔当 $*$-同态。

**命题 13.9（转置）。** 在 $\mathcal M$ 上令 $T(A)=\bigoplus_cA_c^T$。(a) $T$ 是保单位的若尔当
$*$-自同构，且 $T(AB)=T(B)T(A)$；它保持谱，所以 $T$ 与 $T^{-1}=T$ 都是正的。(b) 在 $M_k$ 上，
$k\ge2$，Choi 矩阵 $\sum_{ij}E_{ij}\otimes T(E_{ij})$ 就是交换算子 $F(e_a\otimes e_b)=e_b\otimes e_a$，其
特征值为 $+1$（重数 $k(k+1)/2$）与 $-1$（重数 $k(k-1)/2$）。所以 $T$ 不是完全正的。

*证明。* (a) $(A^T)^*=(A^*)^T$，$(AB)^T=B^TA^T$，因此 $(A\circ B)^T=B^T\circ A^T=A^T\circ B^T$；
$\det(A^T-\lambda)=\det(A-\lambda)$。(b) $\sum_{ij}E_{ij}e_a\otimes E_{ji}e_b=e_b\otimes e_a$。$F^2=1$，且 $F$ 在对称张量上
为 $+1$，在反对称张量上为 $-1$，两者张成整个空间。由定理 6.2（以及命题 6.1），$T$ 不是
CP 的。$\square$

所以若尔当自同构不一定是 $*$-自同构。除这一现象之外，正性决定了若尔当
结构：

**定理 13.10（Kadison；不加证明地引用）。** 设 $\Phi$ 是有单位的 C\*-代数之间的线性双射，
$\Phi(1)=1$，且 $\Phi$ 与 $\Phi^{-1}$ 都是正的。则 $\Phi$ 是若尔当 $*$-同构。（R. V. Kadison,
*A generalized Schwarz inequality and algebraic invariants for operator algebras*, Ann. of Math. 56 (1952) 494–503,
Corollary 5。在 *Isometries of operator algebras*, Ann. of Math. 54 (1951) 325–338 中，他已证明从一个有单位的
C\*-代数到另一个的满线性等距 $\Phi$ 等于酉元 $\Phi(1)$ 乘以一个若尔当 $*$-同构。）

逆命题是显然的：若尔当 $*$-同构把自伴元的平方 $h^2=h\circ h$ 映为平方，
从而把正锥映到正锥上。对于因子 $M_n$，若尔当自同构已被完全确定：
$M_n(\mathbb C)$ 的每个若尔当 $*$-自同构都是 $A\mapsto UAU^*$ 或 $A\mapsto UA^TU^*$，其中 $U$ 为酉矩阵
（N. Jacobson, C. E. Rickart, *Jordan homomorphisms of rings*, Trans. Amer. Math. Soc. 69 (1950) 479–502；I. N.
Herstein, *Jordan homomorphisms*, Trans. Amer. Math. Soc. 81 (1956) 331–341）。我们证明 $n=2$ 的情形。

**命题 13.11.** $\operatorname{Herm}_2(\mathbb C)$ 的每个若尔当自同构都是 $A\mapsto UAU^*$ 或
$A\mapsto UA^TU^*$，其中 $U\in U(2)$。

*证明。* 由命题 13.5(a,b)，若尔当自同构的作用形如 $a1+u\cdot\sigma\mapsto a1+(Ou)\cdot\sigma$，其中
$O\in O(3)$。在这一意义下，记 $R_U\in O(3)$ 为 $A\mapsto UAU^*$ 的矩阵；$R_{UV}=R_UR_V$。(i) $\det R_U=1$：
$U=W\operatorname{diag}(e^{i\theta_1},e^{i\theta_2})W^*$ 经由路径
$U_t=W\operatorname{diag}(e^{it\theta_1},e^{it\theta_2})W^*$ 与 $1$ 相连，而 $\det R_{U_t}\in\{\pm1\}$ 连续地依赖于 $t$。
(ii) 每个 $R\in SO(3)$ 都是某个 $R_U$：与命题 13.6 中一样，存在 $U_1$ 使得
$U_1ZU_1^*=(Re_3)\cdot\sigma$，即 $R_{U_1}e_3=Re_3$。于是 $R_{U_1}^{-1}R\in SO(3)$ 固定 $e_3$，所以它是绕 $e_3$
转某个角度 $\beta$ 的旋转。直接计算表明，对 $U_2=\operatorname{diag}(e^{-i\beta/2},e^{i\beta/2})$，
$U_2XU_2^*=\cos\beta\,X+\sin\beta\,Y$，$U_2YU_2^*=-\sin\beta\,X+\cos\beta\,Y$，所以 $R=R_{U_1U_2}$。(iii)
转置固定 $X,Z$ 并把 $Y\mapsto-Y$：$O_T=\operatorname{diag}(1,-1,1)$，$\det O_T=-1$。因此若
$\det O=1$，则 $O=R_U$；若 $\det O=-1$，则 $OO_T=R_U$，即该自同构为 $A\mapsto UA^TU^*$。$\square$

**复合系统。** 一个因子的若尔当结构无法区分 $T$ 与恒等映射。对于两个系统，这种
多义性变得可见：

**命题 13.12.** 设 $x,x'\in\operatorname{Herm}_n(\mathbb C)$，$y,y'\in\operatorname{Herm}_m(\mathbb C)$。在作用于
$\mathbb C^n\otimes\mathbb C^m$ 的 $\operatorname{Herm}_{nm}(\mathbb C)$ 中，
$$
(x\otimes y)\circ(x'\otimes y')=(x\circ x')\otimes(y\circ y')+\tfrac14[x,x']\otimes[y,y'].
$$
令 $\theta=\mathrm{id}\otimes T$（**部分转置**），$A\circ_\theta B=\theta(\theta(A)\circ\theta(B))$。则
$\circ_\theta$ 是 $\operatorname{Herm}_{nm}(\mathbb C)$ 上的若尔当积，经由 $\theta$ 与 $\circ$ 同构，它
在所有局域可观测量对 $x\otimes1$、$1\otimes y$ 上与 $\circ$ 一致，且
$(x\otimes y)\circ_\theta(x'\otimes y')=(x\circ x')\otimes(y\circ y')-\frac14[x,x']\otimes[y,y']$。它的平方锥是
$\theta(\operatorname{Herm}_{nm}^+)$。对 $n=m=2$ 与 $\Omega=e_1\otimes e_1+e_2\otimes e_2$，态
$\omega=\frac12|\Omega\rangle\langle\Omega|\ge0$ 不在这个锥中：$\theta(\omega)=\frac12F$ 有特征值
$-\frac12$。

*证明。* 令 $j=x\circ x'$，$c=\frac12[x,x']$，$k=y\circ y'$，$d=\frac12[y,y']$：则 $xx'=j+c$，$x'x=j-c$，对 $y$ 也
类似，所以 $\frac12(xx'\otimes yy'+x'x\otimes y'y)=\frac12\big((j+c)\otimes(k+d)+(j-c)\otimes(k-d)\big)=j\otimes k+c\otimes d$。
$\theta$ 是 $\operatorname{Herm}_{nm}$ 的实线性对合：它把 $x\otimes y$ 映为 $x\otimes y^T$，而诸
$x\otimes y$ 在 $\mathbb R$ 上张成 $\operatorname{Herm}_{nm}$（$\operatorname{Herm}_n,\operatorname{Herm}_m$ 的实基的乘积 $h_\alpha\otimes g_\beta$
构成 $M_{nm}$ 的一组复基，而厄米矩阵在这组基下的
系数是实的）。所以 $\circ_\theta$ 是若尔当积，$\theta$ 是同构。由于 $(y^T\circ y'^T)^T=y\circ y'$ 且
$[y^T,y'^T]^T=-[y,y']$，把 $\theta$ 作用于 $(x\otimes y^T)\circ(x'\otimes y'^T)$ 的公式，就得到
$\circ_\theta$ 的公式。在局域可观测量上，对易子项为零（$[1,1]=0$，$[x,1]=0$）。平方元为
$A\circ_\theta A=\theta(\theta(A)^2)$，又因 $\theta^2=\mathrm{id}$，$\omega$ 属于 $\theta(\operatorname{Herm}_{nm}^+)$ 当且仅当 $\theta(\omega)\ge0$。
最后 $\theta(\omega)=\frac12\sum_{ij}E_{ij}\otimes E_{ji}=\frac12F$
（命题 13.9(b)；这正是命题 6.1 证明中的计算）。$\square$

这两个乘积是同样好的若尔当代数，它们在每次只对一个子系统进行的一切测量上都一致，
却在哪些联合算子是正的这一点上，也就是在复合系统的态是什么这一点上，意见不同。区分它们的是
对易子项的符号，而这个符号由 $*$-代数 $M_n\otimes M_m$ 确定，而不是由各因子的若尔当代数确定。（$\circ_\theta$ 是
结合乘积 $A\star B=\theta(\theta(A)\theta(B))$ 的若尔当积，对于该乘积有 $(x\otimes y)\star(x'\otimes y')=xx'\otimes y'y$：
即代数 $M_n\otimes M_m^{\rm op}$。）锥 $\theta(\operatorname{Herm}_{nm}^+)$
给出了 **PPT 判据**（A. Peres, Phys. Rev. Lett. 77 (1996) 1413）：可分态
$\rho=\sum_kp_ka_k\otimes b_k$（$p_k\ge0$，$a_k,b_k\ge0$）同时属于两个锥，因为
$\theta(\rho)=\sum_kp_ka_k\otimes b_k^T$，其中 $b_k^T\ge0$（命题 13.9(a)），且
$a\otimes b=(a^{1/2}\otimes b^{1/2})^2\ge0$。所以不属于 $\theta(\operatorname{Herm}_4^+)$ 的 $\omega$ 是
纠缠的。关于哪些若尔当代数的复合存在的结果，参见 H. Hanche-Olsen, *On the structure and tensor
products of JC-algebras*, Canad. J. Math. 35 (1983) 1059–1074，以及 H. Barnum, M. Graydon, A. Wilce,
*Composites and categories of Euclidean Jordan algebras*, Quantum 4 (2020) 359。

**动力学。** 人们或许希望时间演化会迫使出现结合乘积。事实并非如此：

**命题 13.13.** 设 $H\in\operatorname{Herm}_n(\mathbb C)$，$\delta(z)=i[H,z]$。(a) $\delta$ 把
$\operatorname{Herm}_n$ 映到自身，并且是 $\circ$ 的导子；$e^{t\delta}(z)=e^{itH}ze^{-itH}$ 是若尔当
自同构。(b) 对 $x,y,z\in\operatorname{Herm}_n$，$[L(x),L(y)]z=\frac14[[x,y],z]$。(c) 对某些 $x_k,y_k\in\operatorname{Herm}_n$，
$\delta=4\sum_k[L(x_k),L(y_k)]$：哈密顿动力学仅由若尔当积生成。

*证明。* (a) $(i[H,z])^*=-i[z,H]=i[H,z]$；$[H,\cdot]$ 是结合乘积的导子，从而也是 $\circ$ 的导子。
(b) $4x\circ(y\circ z)=xyz+xzy+yzx+zyx$；减去交换 $x,y$ 后的同一表达式。(c) 由 (b)，只需把
$iH-\frac{\operatorname{Tr}H}n i1$（它与 $iH$ 有相同的对易子）写成若干 $[x_k,y_k]$ 之和。对 $i\ne j$ 以及
$a=E_{ij}+E_{ji}$，$b=i(E_{ij}-E_{ji})$：$[a,E_{ii}]=E_{ji}-E_{ij}$，$[b,E_{ii}]=-i(E_{ij}+E_{ji})$，
$[a,b]=2i(E_{jj}-E_{ii})$；它们在 $\mathbb R$ 上张成无迹反厄米矩阵。$\square$

若尔当代数所不包含的，是可观测量与它所生成的动力学之间的**对应** $H\mapsto i[H,\cdot]$。
例如在 $\operatorname{Herm}_n(\mathbb R)$ 中，导子
$[L(x),L(y)]=\frac14[[x,y],\cdot]$ 来自实反对称矩阵 $[x,y]$，而它不是 $i$ 乘以一个实对称
矩阵。E. M. Alfsen 与 F. W. Shultz 证明了：JBW-代数（形式实若尔当代数的一种无穷维推广；在有限维中
两个概念一致）是某个冯·诺依曼代数的自伴部分，当且仅当
它允许这样一个"动力学对应"；他们还对 C\*-代数的态空间给出了几何刻画
（*Geometry of State Spaces of Operator Algebras*, Birkhäuser, 2003；不加证明地引用）。这样的对应不是
唯一的：在 $\operatorname{Herm}_n(\mathbb C)$ 上，$H\mapsto i[H,\cdot]$ 与 $H\mapsto-i[H,\cdot]$ 都满足要求；
转置交换这两者，因为 $T(i[H,T(z)])=-i[H^T,z]$，而它们分别属于 $M_n$ 与它的反代数——这正是命题 13.12 的
两个符号。

```python
x, y = sa(), sa()
T = jd.transpose_map
assert close(T(jp(x, y)), jp(T(x), T(y))) and close(T(x @ y), T(y) @ T(x))   # Proposition 13.9(a)
assert not close(T(x @ y), T(x) @ T(y))

k = 3
Mk = TypeIAlgebra([k], [k], precision='double')
units = torch.eye(k * k, dtype=dt).reshape(k * k, 1, k, k)                     # E_ij, batch index i*k + j
TE = T(Mk.operator(units)).matrix[:, 0]
choi = TE.reshape(k, k, k, k).permute(0, 2, 1, 3).reshape(k * k, k * k)       # block (i, j) = T(E_ij)
swap = torch.eye(k * k, dtype=dt).reshape(k, k, k, k).permute(0, 1, 3, 2).reshape(k * k, k * k)
w = torch.linalg.eigvalsh(choi)
assert torch.equal(choi, swap) and (w < 0).sum() == k * (k - 1) // 2          # Proposition 13.9(b)

sig = jd.anticommuting_symmetries(3)                                           # X, Y, Z
O_T = torch.stack([V3.from_matrix(s.T)[1:] for s in sig], dim=1)               # transpose on Bloch vectors
assert torch.allclose(O_T, torch.diag(torch.tensor([1.0, -1.0, 1.0], dtype=torch.float64)))

def ptranspose(R, n, m):                                                       # theta = id (x) T
    return R.reshape(n, m, n, m).permute(0, 3, 2, 1).reshape(n * m, n * m)
def herm(n):
    A = torch.randn(n, n, dtype=dt)
    return A + A.mH
jm = lambda A, B: 0.5 * (A @ B + B @ A)
cm = lambda A, B: A @ B - B @ A
a, a2, b, b2 = herm(2), herm(2), herm(3), herm(3)
local = torch.kron(jm(a, a2), jm(b, b2))
assert close(jm(torch.kron(a, b), torch.kron(a2, b2)), local + 0.25 * torch.kron(cm(a, a2), cm(b, b2)))
theta = lambda R: ptranspose(R, 2, 3)
assert close(theta(jm(theta(torch.kron(a, b)), theta(torch.kron(a2, b2)))),
             local - 0.25 * torch.kron(cm(a, a2), cm(b, b2)))                  # Proposition 13.12
Om = torch.tensor([1.0, 0, 0, 1.0], dtype=dt) / math.sqrt(2)
omega = torch.outer(Om, Om.conj())                                             # Bell state
print("spectrum of the partial transpose of the Bell state:", torch.linalg.eigvalsh(ptranspose(omega, 2, 2)))
assert abs(torch.linalg.eigvalsh(ptranspose(omega, 2, 2))[0] + 0.5) < 1e-12

Hm, z1, z2 = herm(3), herm(3), herm(3)
delta = lambda A: 1j * cm(Hm, A)
assert close(delta(jm(z1, z2)), jm(delta(z1), z2) + jm(z1, delta(z2)))         # Proposition 13.13(a)
assert close(jm(z1, jm(z2, Hm)) - jm(z2, jm(z1, Hm)), 0.25 * cm(cm(z1, z2), Hm))   # (b)
```

> **物理意义。** 对于无自旋粒子，转置（在实基中的复共轭）就是**时间
> 反演**。它是单个系统可观测量的一个对称性，但它是反酉的，若只作用于一个
> 子系统，则不是物理操作：命题 13.12 表明，它把 Bell 态变成一个具有负特征值的
> 算子。量子力学使用复 $*$-代数，是因为复合系统需要带符号的对易子，
> 而不仅仅是对称化乘积。

## 13.6 回到迹不等式

命题 11.9（实验 3）断言，对 $X=X^*$ 与 $Y\ge0$ 有 $\operatorname{Tr}|YXY|\le\operatorname{Tr}(Y|X|Y)$（这里
$X,Y$ 是 $\mathcal M$ 的算子，不是 Pauli 矩阵）。其证明从未用到 $Y\ge0$：取
$T=X_\pm^{1/2}Y$，则 $YX_\pm Y=T^*T$ 只需要 $Y=Y^*$，而引理 11.8 被应用于 $YX_+Y$ 与 $YX_-Y$。所以它对所有 $Y=Y^*$ 都成立，
用若尔当语言写出来就是
$$
\operatorname{tr}|P(y)x|\le\operatorname{tr}P(y)|x|,
$$
它在每个欧几里得若尔当代数中都有意义。我们在 $V_m$ 中用一个只用到三条性质的论证来
证明它。

**命题 13.14（自旋因子中的实验 3）。** 对 $x,y\in V_m$，
$\operatorname{tr}|P(y)x|\le\operatorname{tr}P(y)|x|$，且当 $x\in Q$ 或 $-x\in Q$ 时等号成立。

*证明。* 令 $x_\pm=f_\pm(x)$，其中 $f_+(t)=\max(t,0)$，$f_-(t)=\max(-t,0)$（定理 13.4(3)）：$x_\pm\in Q$，$x=x_+-x_-$，
$|x|=x_++x_-$。令 $A=P(y)x_+$，$B=P(y)x_-$；由定理 13.4(5) 它们属于 $Q$，且 $P(y)x=A-B$。设
$P(y)x=\mu_1c_1+\mu_2c_2$ 是一个谱分解（定理 13.4(2)；$c_1=c_w$，$c_2=c_{-w}$，所以 $c_1+c_2=1$，
$c_1\circ c_2=0$），$c$ 是 $\mu_i>0$ 的那些 $c_i$ 之和，并令
$s=2c-1=\sum_i\operatorname{sgn}'(\mu_i)c_i$，其中当 $t>0$ 时 $\operatorname{sgn}'(t)=1$，否则为 $-1$。则
$s\circ P(y)x=\sum_i|\mu_i|c_i=|P(y)x|$（由定理 13.4(2,3)，当
$\mu_1=\mu_2$ 时这与 $w$ 的选取无关），且 $c,1-c\in Q$（幂等元是平方元）。把自对偶性中容易的一半（定理 13.4(4)：
对 $a,b\in Q$ 有 $\operatorname{tr}(a\circ b)\ge0$）应用于 $\operatorname{tr}(A\circ(1-c))$ 与 $\operatorname{tr}(B\circ c)$，得
$$
\operatorname{tr}|P(y)x|=\operatorname{tr}(A\circ c)-\operatorname{tr}(A\circ(1-c))-\operatorname{tr}(B\circ c)+\operatorname{tr}(B\circ(1-c))
\le\operatorname{tr}A+\operatorname{tr}B=\operatorname{tr}P(y)|x|.
$$
若 $x\in Q$，则 $x_-=0$，$B=0$，$|P(y)x|=A$，所以两边都等于 $\operatorname{tr}A$；$-x\in Q$ 时同理。
$\square$

证明只用到：(i) 对每个 $y$ 有 $P(y)Q\subseteq Q$；(ii) 谱定理：每个 $z$ 都可写成
$\sum_i\mu_ic_i$，其中 $c_i$ 是和为 $1$ 的正交幂等元，且 $f(z)=\sum_if(\mu_i)c_i$（因此 $x_\pm$ 与诸
幂等元都是平方元）；(iii) $\operatorname{tr}$ 是线性的，$\operatorname{tr}(\sum_i\mu_ic_i)=\sum_i\mu_i\operatorname{tr}c_i$，
其中 $\operatorname{tr}c_i\ge0$，并且对 $a,b\in Q$ 有 $\operatorname{tr}(a\circ b)\ge0$。在秩为 $r$ 的欧几里得若尔当代数中，
(ii) 对 $r$ 个本原幂等元 $c_i$（一个若尔当标架）成立，迹为
$\operatorname{tr}z=\sum_i\mu_i$，形式 $\operatorname{tr}(x\circ y)$ 是结合的，即
$\operatorname{tr}((x\circ y)\circ z)=\operatorname{tr}(x\circ(y\circ z))$，并且 (i) 成立 [FK94, Chapters II–III]；
我们不证明这些事实。于是 (iii) 的最后一部分可由 (i) 推出：展开 $P(a^{1/2})=2L(a^{1/2})^2-L(a)$
并利用结合性，得 $\operatorname{tr}(a\circ b)=\operatorname{tr}P(a^{1/2})b$，这是 $Q$ 中一个元素的迹，
因而 $\ge0$。所以该不等式在每个欧几里得若尔当代数中都成立——包括
$\operatorname{Herm}_r(\mathbb H)$ 与 Albert 代数，它们不若尔当同构于任何复
$*$-代数的自伴部分。在 $M_{sa}$ 中，用 $\operatorname{Tr}$ 与 $P(y)x=yxy$ 进行同样的论证，就重新证明了
命题 11.9。在 $V_m$ 中，该不等式也可以经由 $\phi$ 从命题 11.9 得出，
因为 $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$ 且 $\phi(|x|)=|\phi(x)|$；代码对两者都做了检验：

```python
for m in (2, 4, 7):
    Vm = SpinFactor(m)
    xs, ys = Vm.random(20_000), Vm.random(20_000)
    lhs = Vm.trace(Vm.apply_function(Vm.quadratic_representation(ys, xs), torch.abs))
    rhs = Vm.trace(Vm.quadratic_representation(ys, Vm.apply_function(xs, torch.abs)))
    definite = Vm.determinant(xs) >= 0                                 # x or -x in the cone
    print(f"m = {m}: min gap {(rhs - lhs).min().item():.1e}, max gap {(rhs - lhs).max().item():.1f}, "
          f"definite x: {definite.double().mean().item():.3f}")
    assert (rhs - lhs >= -1e-12 * rhs).all()                            # Proposition 13.14
    assert ((rhs - lhs)[definite].abs() <= 1e-12 * rhs[definite]).all()
    Xs, Ys = Vm.to_matrix(xs[:200]), Vm.to_matrix(ys[:200])
    tn = torch.linalg.svdvals(Ys @ Xs @ Ys).sum(-1)                     # Tr |phi(y) phi(x) phi(y)|
    assert (tn - Vm.representation_dim / 2 * lhs[:200]).abs().max() < 1e-10 * tn.max()
```

差值 $\operatorname{tr}P(y)|x|-\operatorname{tr}|P(y)x|$ 对定号的 $x$ 为零，对典型的
不定 $x$ 为正（对这些 Gauss 样本，打印出的最大值的量级为 $10^1$–$10^2$）。

## 13.7 可与之比较的工具

据我们所知，目前没有主流的若尔当代数数值库；本章中的运算
通常是用矩阵手工写出的。计算机代数涵盖了代数方面：SageMath 的 `JordanAlgebra` 可以构造
特殊若尔当代数、对称双线性形式的若尔当代数（其中包括自旋因子），以及由一个
八元数代数构造的 27 维例外若尔当代数；Magma 也有若尔当代数的构造函数，
包括例外的那一个。这些是精确的代数工具，而不是批量化的数值工具。洛伦兹锥在
优化中是标准的：SOCP 求解器（ECOS、SCS、Clarabel、MOSEK，可通过 cvxpy 的 `SOC` 约束调用）使用 $V_m$ 的谱
分解（用于向锥投影）或它的二次表示（用于内点法的 Nesterov–Todd 缩放），
半定规划求解器对 $\operatorname{Herm}_n$ 也是如此。本库补充的是：若尔当积、二次
表示与转置，作为 $\bigoplus_cM_{k_c}$ 上批量化、可微的算子；带有
函数演算与矩阵表示的自旋因子；以及对给定矩阵所生成的若尔当代数的随机化计算，
并与它们所生成的 $*$-代数进行比较。

## 习题

1. （**洛伦兹锥的齐性。**）设 $y\in V_m$，$\det y\ne0$。证明 $y^{-1}=(a,-u)/\det y$ 满足
   $y\circ y^{-1}=1$，$P(y)P(y^{-1})=\mathrm{id}$，$P(y)$ 把 $Q$ 的内部映到自身上，并且
   对内部的 $x$ 有 $P(x^{1/2})1=x$。由此推出 $Q$ 的内部是齐性的。
2. （**各向同性态。**）对 $p\in[0,1]$，在 $\mathbb C^2\otimes\mathbb C^2$ 上令 $\rho_p=p\,\omega+(1-p)\frac14$
   （命题 13.12）。计算 $\theta(\rho_p)$ 的谱，并证明当 $p>\frac13$ 时 $\rho_p$ 是纠缠的。
3. （**边界点的乘积。**）设 $\mathbb R^m$ 中 $|w|=|w'|=1$。证明 $(1,w)\circ(1,w')\in Q$ 当且仅当
   $w'=w$ 或 $w'=-w$：洛伦兹锥的两条边界射线的若尔当积几乎从不落在锥中。
4. （**实量子力学。**）证明对 $x,y\in\operatorname{Herm}_n(\mathbb R)$，$\operatorname{Herm}_n(\mathbb R)$ 的导子
   $[L(x),L(y)]$ 是 $z\mapsto[K,z]$，其中 $K$ 是实反对称矩阵，并证明对
   $n\ge2$，这样的映射除非为 $0$，否则永远不是 $z\mapsto i[H,z]$（$H\in\operatorname{Herm}_n(\mathbb R)$）。

<details>
<summary>解答</summary>

1. $(a,u)\circ(a,-u)=(a^2-|u|^2,0)$，由此得 $y\circ y^{-1}=1$。$\phi(y)$ 与 $\phi(y^{-1})$ 对易（两者都属于
   $\operatorname{span}\{1,u\cdot s\}$），所以 $\phi(y)\phi(y^{-1})=\phi(y)\circ\phi(y^{-1})=1$，再由定理 13.4(5) 得
   $\phi(P(y)P(y^{-1})x)=\phi(x)$；而 $\phi$ 是单射。$P(y)$ 把 $Q$ 映入 $Q$，它的逆 $P(y^{-1})$ 也是如此，所以
   $P(y)Q=Q$；作为线性同胚，它把内部映到内部上。最后
   $P(y)1=2y\circ y-y^2=y^2$，所以 $P(x^{1/2})1=x$，而当 $x$ 是内点时 $P(x^{1/2})$ 可逆
   （$\det x^{1/2}=\sqrt{\lambda_+\lambda_-}>0$）。任意两个内点都由
   $P(x'^{1/2})P(x^{1/2})^{-1}$ 联系起来。

```python
V4 = SpinFactor(4)
xi = V4.random(100, cone=True)                                   # interior points (almost surely)
r = V4.sqrt(xi)
assert (V4.quadratic_representation(r, V4.unit((100,))) - xi).abs().max() < 1e-10
yv = V4.random(100)
back = V4.quadratic_representation(V4.inverse(yv), V4.quadratic_representation(yv, xi))
assert (back - xi).abs().max() < 1e-8 * xi.abs().max()
assert V4.in_cone(V4.quadratic_representation(yv, xi)).all()
```

2. 由命题 13.12，$\theta(\rho_p)=\frac p2F+\frac{1-p}4$。$F$ 的特征值为 $1$（三重）与 $-1$（一重），所以
   谱为 $\frac{1+p}4$（三重）与 $\frac{1-3p}4$，后者为负当且仅当 $p>\frac13$；此时由 PPT 判据，$\rho_p$ 不是
   可分的。（对 $2\times2$ 系统，逆命题也成立：M., P. and R. Horodecki, Phys.
   Lett. A 223 (1996) 1；我们不加证明。）

```python
for pp in (0.2, 1 / 3, 0.5):
    rho_p = pp * omega + (1 - pp) * torch.eye(4, dtype=dt) / 4
    lam_min = torch.linalg.eigvalsh(ptranspose(rho_p, 2, 2))[0].item()
    assert abs(lam_min - (1 - 3 * pp) / 4) < 1e-12
```

3. $(1,w)\circ(1,w')=(1+t,\,w+w')$，其中 $t=w\cdot w'\in[-1,1]$，$|w+w'|=\sqrt{2+2t}$。它属于 $Q$ 当且仅当
   $1+t\ge\sqrt{2(1+t)}$，即当且仅当 $t=-1$（乘积为 $0$）或 $1+t\ge2$，即 $t=1$。

```python
V3c = SpinFactor(3)
wv = torch.nn.functional.normalize(torch.randn(2, 1000, 3, dtype=torch.float64), dim=-1)
prod = V3c.product(V3c.element(1.0, wv[0]), V3c.element(1.0, wv[1]))
assert not V3c.in_cone(prod).any()
assert V3c.in_cone(V3c.product(V3c.element(1.0, wv[0]), V3c.element(1.0, wv[0]))).all()
```

4. 由命题 13.13(b)，$[L(x),L(y)]z=[K,z]$，其中 $K=\frac14[x,y]$，且 $K^T=\frac14[y^T,x^T]=-K$ 是实的。对实的 $H$
   与 $z$，$i[H,z]$ 是纯虚的，而 $[K,z]$ 是实的，所以两者只有在都为零时才相等。因此在实
   量子力学中，动力学的生成元 $K$ 不是可观测量。

</details>

**下一步。** 这是本课程的最后一章。有限维
理论的标准参考书是 [FK94]；K. McCrimmon, *A Taste of Jordan Algebras* (Springer 2004) 讲述了代数理论与历史，
H. Hanche-Olsen and E. Størmer, *Jordan Operator Algebras* (Pitman 1984) 以及 Alfsen–Shultz (2003) 讲述了
无穷维理论及其与 C\*-代数的关系。关于复合系统在从操作公理重建
量子理论中的作用，参见上面引用的 Barnum 与 Wilce 的论文。关于第 11 章之外的迹不等式，参见
E. Carlen, *Trace inequalities and quantum entropy* (Contemp. Math. 529, 2010)。二次
表示与洛伦兹锥在锥优化中再次出现：参见 Alizadeh–Goldfarb (2003) 以及其中综述的
若尔当代数内点法。
