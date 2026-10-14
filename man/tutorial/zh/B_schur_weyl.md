[English](../B_schur_weyl.md) · [Español](../es/B_schur_weyl.md) · [Français](../fr/B_schur_weyl.md) · [Deutsch](../de/B_schur_weyl.md) · **中文** · [日本語](../ja/B_schur_weyl.md) · [Русский](../ru/B_schur_weyl.md)

# 教程，附录 B. Schur–Weyl 对偶与 Weingarten 演算

本教程配合讲义的[附录 B](../../lecture_notes/zh/B_schur_weyl.md)，该附录证明了 Schur–Weyl 对偶（[定理 B.6](../../lecture_notes/zh/B_schur_weyl.md#r-B-6)）和任意阶的
Weingarten 公式（[定理 B.13](../../lecture_notes/zh/B_schur_weyl.md#r-B-13)）。我们检验对偶中的维数、本库计算出的 $`(\mathbb C^d)^{\otimes n}`$ 的分解、
Weingarten 函数与 Haar 酉矩阵上蒙特卡罗平均的比较、Wigner–Eckart 定理，以及正交群的 Haar 测度。

<a id="s-B-3"></a>

## B.3 Schur–Weyl 对偶

*讲义：* [B.3 节](../../lecture_notes/zh/B_schur_weyl.md#s-B-3)；[引理 B.1](../../lecture_notes/zh/B_schur_weyl.md#r-B-1)、[引理 B.2](../../lecture_notes/zh/B_schur_weyl.md#r-B-2)、[命题 B.3](../../lecture_notes/zh/B_schur_weyl.md#r-B-3)、[引理 B.5](../../lecture_notes/zh/B_schur_weyl.md#r-B-5)、[定理 B.6](../../lecture_notes/zh/B_schur_weyl.md#r-B-6)。

代码对较小的 $`d`$ 和 $`n`$ 比较四个维数：$`\mathrm{span}\,\pi(S_n)`$ 的维数、两个随机的 $`V^{\otimes n}`$ 的交换子的维数
（它包含 $`\{V^{\otimes n}\}'=\mathrm{span}\,\pi(S_n)`$，并且在所有打印出的情形中与之相等），以及许多随机的
$`A^{\otimes n}`$ 和 $`V^{\otimes n}`$ 的张成的维数，并与 $`\binom{d^2+n-1}n`$ 比较（[命题 B.3](../../lecture_notes/zh/B_schur_weyl.md#r-B-3)，[定理 B.6](../../lecture_notes/zh/B_schur_weyl.md#r-B-6)）。
算子 $`\pi(\sigma)`$ 由 `decompose.permutation_representation(n, d)` 构造；代码检验它遵循[引理 B.1(a)](../../lecture_notes/zh/B_schur_weyl.md#r-B-1) 的约定，
并检验[引理 B.2](../../lecture_notes/zh/B_schur_weyl.md#r-B-2) 以及[引理 B.5](../../lecture_notes/zh/B_schur_weyl.md#r-B-5) 的 Newton 递推，后者用于定理 B.6(b) 的证明。

```python
import functools
import itertools
import math
import torch
import torch_vn_algebra.decompose as dc

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    return (torch.as_tensor(X) - torch.as_tensor(Y)).abs().max().item() < tol

def haar(d, B):                                        # Theorem 5.3: QR of Ginibre matrices, phases fixed
    Q, R = torch.linalg.qr(torch.randn(B, d, d, dtype=dt))
    r = R.diagonal(dim1=-2, dim2=-1)
    return (Q * (r / r.abs()).unsqueeze(-2)).contiguous()

def power(A, n):                                       # A^{(x) n} for a batch A of shape (B, d, d)
    out = A
    for _ in range(n - 1):
        out = torch.einsum('bij,bkl->bikjl', out, A).reshape(A.shape[0], out.shape[1] * A.shape[1], -1)
    return out

def rank(X, tol=1e-9):                                 # dimension of the span of a batch of matrices
    return torch.linalg.matrix_rank(X.reshape(X.shape[0], -1), atol=tol).item()

def commutant_dim(ops, tol=1e-9):                      # dim {X : XA = AX for all A in ops}
    eye = torch.eye(ops.shape[-1], dtype=dt)
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    return int((torch.linalg.svdvals(L) < tol).sum())

for d, n in ((2, 2), (2, 3), (2, 4), (3, 2), (3, 3)):
    P = dc.permutation_representation(n, d)            # the n! operators pi(sigma)
    N_sym = math.comb(d * d + n - 1, n)                # Proposition B.3
    span_A = rank(power(torch.randn(N_sym + 5, d, d, dtype=dt), n))
    span_V = rank(power(haar(d, N_sym + 5), n))
    comm_V = commutant_dim(power(haar(d, 2), n))       # commutant of two random V^{(x) n}
    print(d, n, rank(P), comm_V, span_A, span_V, N_sym)
    assert span_A == span_V == N_sym == commutant_dim(P)                     # Theorem B.6(a)
    assert comm_V == rank(P)              # numerical observation: two V suffice here (B.6(c) gives >=)
    assert rank(P) == (math.factorial(n) if n <= d else math.comb(2 * n, n) // (n + 1))   # B.1(d), B.10

A = torch.randn(3, 2, 2, dtype=dt)                     # Lemma B.2 for n = 3
lhs = sum(functools.reduce(torch.kron, [A[s[0]], A[s[1]], A[s[2]]]) for s in itertools.permutations(range(3)))
rhs = sum((-1) ** (3 - sum(eps)) * functools.reduce(torch.kron, [sum(e * a for e, a in zip(eps, A))] * 3)
          for eps in itertools.product((0, 1), repeat=3))
assert close(lhs, rhs)

xs = torch.randn(3, 2, dtype=dt)                       # Lemma B.1(a): pi(s) x_1 (x) x_2 (x) x_3 = (x) x_{s^-1(p)}
for s in itertools.permutations(range(3)):
    inv = [s.index(m) for m in range(3)]
    Ps = dc.permutation_representation(3, 2, [list(s)])[0]
    assert close(Ps @ functools.reduce(torch.kron, list(xs)), functools.reduce(torch.kron, [xs[inv[p]] for p in range(3)]))

def d_pi(X, n):                                        # sum_m 1 (x) .. (x) X (x) .. (x) 1
    eye = torch.eye(X.shape[0], dtype=dt)
    return sum(functools.reduce(torch.kron, [X if l == m else eye for l in range(n)]) for m in range(n))

n, A = 4, torch.randn(2, 2, dtype=dt)
e = [torch.eye(2 ** n, dtype=dt)]
for k in range(1, n + 1):                              # Lemma B.5 with x_m = A^(m), p_j = dPi(A^j)
    e.append(sum((-1) ** (j - 1) * e[k - j] @ d_pi(torch.linalg.matrix_power(A, j), n)
                 for j in range(1, k + 1)) / k)
assert close(e[n], functools.reduce(torch.kron, [A] * n), 1e-8)
```

<a id="s-B-4"></a>

## B.4 $`(\mathbb C^d)^{\otimes n}`$ 的分解

*讲义：* [B.4 节](../../lecture_notes/zh/B_schur_weyl.md#s-B-4)；[定理 B.8](../../lecture_notes/zh/B_schur_weyl.md#r-B-8)、[推论 B.9](../../lecture_notes/zh/B_schur_weyl.md#r-B-9)、[推论 B.10](../../lecture_notes/zh/B_schur_weyl.md#r-B-10)。

代码用 `decompose.from_group` 分解 $`(\mathbb C^3)^{\otimes3}`$ 上的 $`\pi(S_3)`$，其 `k_factors` 与
`multiplicities` 就是 $`k_c`$ 与 $`m_c`$，并检验[定理 B.8](../../lecture_notes/zh/B_schur_weyl.md#r-B-8) 的维数；然后对 $`L=5,6`$ 分解
$`(\mathbb C^2)^{\otimes L}`$ 上的 $`\pi(S_L)`$，把数对 $`(k_c,m_c)`$ 与[推论 B.9](../../lecture_notes/zh/B_schur_weyl.md#r-B-9) 的 $`(m_S,2S+1)`$ 比较，并把
$`\sum_ck_c^2`$ 与[推论 B.10](../../lecture_notes/zh/B_schur_weyl.md#r-B-10) 的 Catalan 数比较。`from_group` 背后的算法在[教程 4，4.3 节](04_symmetries.md#s-4-3)
中解释。

```python
D33 = dc.from_group(dc.permutation_representation(3, 3))          # S_3 on (C^3)^{(x)3}
print(D33)
assert D33.verify()['ok']
assert sum(k * k for k in D33.k_factors) == 6 == math.factorial(3)            # Theorem B.8, n <= d
assert sum(m * m for m in D33.multiplicities) == math.comb(11, 3)              # 165

for L in (5, 6):
    DL = dc.from_group(dc.permutation_representation(L, 2))
    pairs = sorted((math.comb(L, j) - (math.comb(L, j - 1) if j else 0), L - 2 * j + 1)
                   for j in range(L // 2 + 1))                                 # (m_S, 2S + 1), S = L/2 - j
    print(L, DL.k_factors, DL.multiplicities)
    assert sorted(zip(DL.k_factors, DL.multiplicities)) == pairs               # Corollary B.9
    assert sum(k * k for k in DL.k_factors) == math.comb(2 * L, L) // (L + 1)  # Corollary B.10
    assert sum(m * m for m in DL.multiplicities) == math.comb(L + 3, 3)
```

<a id="s-B-5"></a>

## B.5 Weingarten 演算

*讲义：* [B.5 节](../../lecture_notes/zh/B_schur_weyl.md#s-B-5)；[引理 B.11](../../lecture_notes/zh/B_schur_weyl.md#r-B-11)、[定理 B.13](../../lecture_notes/zh/B_schur_weyl.md#r-B-13)、[推论 B.14](../../lecture_notes/zh/B_schur_weyl.md#r-B-14)。

代码由伪逆计算 $`\mathrm{Wg}`$，检验[定理 B.13(a)、(b)](../../lecture_notes/zh/B_schur_weyl.md#r-B-13) 和
[定理 5.10](../../lecture_notes/zh/05_random_matrices.md#r-5-10)，并把 $`d=2`$ 时 $`n=3`$ 阶的整个矩张量（一个 $`G`$ 奇异的情形）以及
随机 $`X`$ 的 twirl（[引理 B.11](../../lecture_notes/zh/B_schur_weyl.md#r-B-11)）与 Haar 酉矩阵上的蒙特卡罗平均进行比较；其中的元素
$`\mathbb E|U_{11}|^6=\frac{3!}{2\cdot3\cdot4}`$ 就是[推论 B.14](../../lecture_notes/zh/B_schur_weyl.md#r-B-14)。

```python
def n_cycles(p):                                        # number of cycles of a permutation in one-line form
    seen, count = set(), 0
    for i in range(len(p)):
        count += i not in seen
        while i not in seen:
            seen.add(i)
            i = p[i]
    return count

def compose(s, t):                                     # (s t)(m) = s(t(m))
    return tuple(s[m] for m in t)

def inverse(s):
    return tuple(sorted(range(len(s)), key=lambda m: s[m]))

def weingarten(n, d):
    perms = list(itertools.permutations(range(n)))
    G = torch.tensor([[float(d) ** n_cycles(compose(inverse(s), t)) for t in perms] for s in perms],
                     dtype=torch.float64)
    return perms, G, torch.linalg.pinv(G, hermitian=True)

for n, d in ((2, 2), (3, 2), (3, 3), (4, 2), (4, 3)):
    perms, G, W = weingarten(n, d)
    pos = {p: a for a, p in enumerate(perms)}
    Wg = {p: W[0, pos[p]].item() for p in perms}        # perms[0] is the identity
    assert all(abs(W[a, b] - Wg[compose(inverse(s), t)]) < 1e-12
               for a, s in enumerate(perms) for b, t in enumerate(perms))       # Theorem B.13(a)
    assert all(abs(Wg[p] - Wg[compose(compose(inverse(r), p), r)]) < 1e-12 for p in perms for r in perms)
    assert (torch.linalg.matrix_rank(G).item() == math.factorial(n)) == (n <= d)  # Theorem B.13(b)
for d in (2, 3, 5):                                     # Theorem 5.10
    _, _, W = weingarten(2, d)
    assert close(W[0], torch.tensor([1 / (d * d - 1), -1 / (d * (d * d - 1))], dtype=torch.float64), 1e-12)

n, d, B = 3, 2, 100_000
perms, G, W = weingarten(n, d)
P = dc.permutation_representation(n, d, [list(p) for p in perms])        # pi(sigma), same order as W
Phi = lambda X: torch.einsum('st,tij,ij,skl->kl', W.to(dt), P.conj(), X, P)   # Theorem B.13(c)
U3 = power(haar(d, B), n)                                                 # B samples of U^{(x)3}
moments = torch.einsum('bx,by->xy', U3.reshape(B, -1), U3.reshape(B, -1).conj()) / B
second = torch.einsum('bx,by->xy', U3.reshape(B, -1).abs() ** 2, U3.reshape(B, -1).abs() ** 2) / B
err = (second - moments.abs() ** 2).clamp(min=0).sqrt() / math.sqrt(B)   # standard errors
N = d ** n
exact = torch.stack([Phi(torch.outer(torch.eye(N, dtype=dt)[j], torch.eye(N, dtype=dt)[jp]))
                     for j in range(N) for jp in range(N)])               # [(j, j'), i, i']
exact = exact.reshape(N, N, N, N).permute(2, 0, 3, 1).reshape(N * N, N * N)   # [(i, j), (i', j')]
z = ((moments - exact).abs() / (err + 1e-15)).max().item()
print(f"{N ** 4} moments of order 3 for d = 2: largest deviation {z:.2f} standard errors")
assert ((moments - exact).abs() <= 5 * err + 1e-12).all()                 # Theorem B.13(d), singular G
assert abs(exact[0, 0].item() - 1 / 4) < 1e-12                            # Corollary B.14: 3!/(2*3*4)

X = torch.randn(N, N, dtype=dt)
T = Phi(X)
coef = torch.linalg.lstsq(P.reshape(len(perms), -1).T, T.reshape(-1, 1)).solution
assert close(Phi(T), T) and close(torch.einsum('s,sij->ij', coef[:, 0], P), T)   # a projection into span pi(S_n)
Vn = power(haar(d, 1), n)[0]
assert close(Vn @ T, T @ Vn)                                             # = {V^{(x)n}}'
assert close(Phi(P[3]), P[3])
Y = U3 @ X @ U3.mH                                                        # samples of U X U*
assert ((Y.mean(0) - T).abs() <= 5 * Y.std(0) / math.sqrt(B) + 1e-12).all()   # Lemma B.11
```

<a id="s-B-6"></a>

## B.6 张量算子与 $`SU(2)`$ 的 Wigner–Eckart 定理

*讲义：* [B.6 节](../../lecture_notes/zh/B_schur_weyl.md#s-B-6)；[引理 B.15](../../lecture_notes/zh/B_schur_weyl.md#r-B-15)、[定理 B.16](../../lecture_notes/zh/B_schur_weyl.md#r-B-16)。

代码对六个自旋 $`\frac12`$ 的格点自旋检验[引理 B.15](../../lecture_notes/zh/B_schur_weyl.md#r-B-15)（`decompose.total_spin`、
`decompose.spin_matrices`），并以不依赖基的形式检验[定理 B.16](../../lecture_notes/zh/B_schur_weyl.md#r-B-16)：在 `decompose.from_generators`
计算出的分解中（其扇区基，即 `V` 的列块，不必是 $`|S,M\rangle`$），把 $`(T_q)_q`$ 在两个多重态之间的块看作行为
$`(M',M,q)`$、列为 $`(a,b)`$ 的矩阵，它的秩至多为一，并且在 (a) 的范围之外为零。该定理是 `su2.spin_correlation`
背后的原理，[教程 9，9.3 节](09_spin_chains.md#s-9-3)把它与 `'direct'` 方法进行比较检验。

```python
L = 6
S6 = dc.total_spin(L)
Sp, Sm = S6[0] + 1j * S6[1], S6[0] - 1j * S6[1]
D6 = dc.from_generators(S6)
print(D6)
off = [0]
for k, m in zip(D6.k_factors, D6.multiplicities):
    off.append(off[-1] + k * m)
sx, sy, sz = dc.spin_matrices(0.5)
site = lambda x, i: functools.reduce(torch.kron, [x if l == i else torch.eye(2, dtype=dt) for l in range(L)])
spherical = lambda X, Y, Z: {1: -(X + 1j * Y) / math.sqrt(2), 0: Z, -1: (X - 1j * Y) / math.sqrt(2)}
for i in (0, 3):
    T = spherical(site(sx, i), site(sy, i), site(sz, i))
    zero = torch.zeros_like(T[0])
    for q in (-1, 0, 1):                                                   # Lemma B.15(a), (b)
        assert close(S6[2] @ T[q] - T[q] @ S6[2], q * T[q])
        assert close(Sp @ T[q] - T[q] @ Sp, math.sqrt(2 - q * (q + 1)) * T.get(q + 1, zero))
        assert close(Sm @ T[q] - T[q] @ Sm, math.sqrt(2 - q * (q - 1)) * T.get(q - 1, zero))
    for c, (kc, mc) in enumerate(zip(D6.k_factors, D6.multiplicities)):
        for e, (ke, me) in enumerate(zip(D6.k_factors, D6.multiplicities)):
            Vc, Ve = D6.V[:, off[c]:off[c + 1]], D6.V[:, off[e]:off[e + 1]]
            block = torch.stack([(Ve.mH @ T[q] @ Vc).reshape(ke, me, kc, mc) for q in (-1, 0, 1)])
            sv = torch.linalg.svdvals(block.permute(0, 1, 3, 2, 4).reshape(3 * ke * kc, me * mc))
            S, S2 = (kc - 1) / 2, (ke - 1) / 2
            if abs(S - S2) <= 1 and S + S2 >= 1:
                assert sv.numel() == 1 or sv[1] < 1e-10                        # Theorem B.16(b): rank one
            else:
                assert sv[0] < 1e-10                                           # Theorem B.16(a)
X, Y, Z = (site(s, 0) for s in (sx, sy, sz))
X2, Y2, Z2 = (site(s, 1) for s in (sx, sy, sz))
T, T2 = spherical(X, Y, Z), spherical(X2, Y2, Z2)
assert close(X @ X2 + Y @ Y2 + Z @ Z2, sum((-1) ** q * T[q] @ T2[-q] for q in (-1, 0, 1)))   # Lemma B.15(a)
```

<a id="s-B-7"></a>

## B.7 本课程中使用的 Haar 测度

*讲义：* [B.7 节](../../lecture_notes/zh/B_schur_weyl.md#s-B-7)；[定理 B.17](../../lecture_notes/zh/B_schur_weyl.md#r-B-17)。

由[定理 B.17](../../lecture_notes/zh/B_schur_weyl.md#r-B-17)，实高斯矩阵的 QR 分解 $`Z=OR`$ 在 $`R_{ii}>0`$ 的归一化下，其正交因子 $`O`$ 在 $`O(k)`$ 上服从
Haar 分布，且 $`\mathbb E[O_{ij}O_{lm}]=\delta_{il}\delta_{jm}/k`$；对任何其他 QR 分解 $`Z=O'R'`$，
$`O=O'\mathrm{diag}(\mathrm{sgn}\,R'_{ii})`$。

代码用蒙特卡罗方法检验二阶矩、$`\mathbb E\mathrm{Tr}\,O=0`$ 以及 $`\mathbb E(\mathrm{Tr}\,O)^2=\sum_{i,l}\mathbb E[O_{ii}O_{ll}]=1`$。
它还记录了关于 `torch.linalg.qr` 的 QR 例程的一个*数值观察*（这是关于该例程的陈述，而不是定理）：不作符号修正时，
其因子 $`Q'`$ 的 $`\mathbb E\mathrm{Tr}\,Q'`$ 与 $`0`$ 相差许多个标准误差，所以它不服从 Haar 分布。

```python
k, B = 3, 100_000
Zr = torch.randn(B, k, k, dtype=torch.float64)
Q, R = torch.linalg.qr(Zr)                              # LAPACK: diag(R) of either sign
O = Q * torch.sign(R.diagonal(dim1=-2, dim2=-1)).unsqueeze(-2)     # Theorem B.17
assert close(O.mT @ O, torch.eye(k, dtype=torch.float64).expand(B, k, k))
assert (R.diagonal(dim1=-2, dim2=-1) * torch.sign(R.diagonal(dim1=-2, dim2=-1)) > 0).all()
second = torch.einsum('bij,blm->ijlm', O, O) / B                   # E O_ij O_lm
expected = torch.einsum('il,jm->ijlm', torch.eye(k, dtype=torch.float64), torch.eye(k, dtype=torch.float64)) / k
err = (torch.einsum('bij,blm->ijlm', O ** 2, O ** 2) / B - second ** 2).clamp(min=0).sqrt() / math.sqrt(B)
assert ((second - expected).abs() <= 5 * err).all()
tr, tr_raw = O.diagonal(dim1=-2, dim2=-1).sum(-1), Q.diagonal(dim1=-2, dim2=-1).sum(-1)
print(f"E Tr O = {tr.mean():.4f}, E (Tr O)^2 = {(tr ** 2).mean():.4f};  without the sign correction "
      f"E Tr Q = {tr_raw.mean():.4f}")
assert abs(tr.mean()) < 5 * tr.std() / math.sqrt(B)                # E Tr O = 0
assert abs((tr ** 2).mean() - 1) < 5 * (tr ** 2).std() / math.sqrt(B)   # E (Tr O)^2 = sum_i 1/k = 1
assert abs(tr_raw.mean()) > 20 * tr_raw.std() / math.sqrt(B)       # LAPACK's signs: not Haar
```

**下一章：** [教程，附录 C. 特征值密度与 Page 公式](C_eigenvalue_densities.md)。
