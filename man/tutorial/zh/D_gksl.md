[English](../D_gksl.md) · [Español](../es/D_gksl.md) · [Français](../fr/D_gksl.md) · [Deutsch](../de/D_gksl.md) · **中文** · [日本語](../ja/D_gksl.md) · [Русский](../ru/D_gksl.md)

# 教程，附录 D. 量子动力学半群的生成元

本教程配合讲义的[附录 D](../../lecture_notes/zh/D_gksl.md)，该附录由生成元的 Choi 矩阵证明了 GKSL 定理（[定理 D.4](../../lecture_notes/zh/D_gksl.md#r-D-4)）
及其在有扇区的代数上的版本（[定理 D.7](../../lecture_notes/zh/D_gksl.md#r-D-7)），以及算子代数纠错定理（[定理 D.9](../../lecture_notes/zh/D_gksl.md#r-D-9)）。我们用本库的
Lindblad 例程和 `numpy`，把这些证明作为算法在随机例子上运行。

<a id="s-D-2"></a>

## D.2 条件完全正性

*讲义：* [D.2 节](../../lecture_notes/zh/D_gksl.md#s-D-2)；[引理 D.1](../../lecture_notes/zh/D_gksl.md#r-D-1)、[定理 D.3](../../lecture_notes/zh/D_gksl.md#r-D-3)。

我们在 $`M_3`$ 上一个带任意跳跃的 GKSL 生成元上检验定理 D.3，它由 `lindblad_superoperator` 构造
（即 $`\mathcal L`$ 作用在按行展开的 $`\mathrm{vec}\,\rho`$ 上的矩阵 $`S`$，因此与 6.1 节一样
$`J_{(i,a),(j,b)}=S_{(a,b),(i,j)}`$），并且也在 $`\theta-\mathrm{id}`$ 上检验：

```python
import math
import numpy as np
import scipy.linalg as sla
import torch
from torch_vn_algebra import TypeIAlgebra, channels, dynamics

torch.manual_seed(0)
rng = np.random.default_rng(0)
f64, c128 = torch.float64, torch.complex128

def choi_of(S, n):
    """J_{(i,a),(j,b)} = S_{(a,b),(i,j)}: the Choi matrix of the map with row-major matrix S (Section 6.1)."""
    return S.reshape(n, n, n, n).transpose(2, 0, 3, 1).reshape(n * n, n * n)

def superop(f, n):
    """Row-major matrix of a linear map f on n x n matrices."""
    E = np.eye(n * n).reshape(n * n, n, n)
    return np.stack([f(X).reshape(-1) for X in E], axis=1)

vec_of = lambda K: K.T.reshape(-1)                        # (1 (x) K) Omega = sum_i e_i (x) K e_i
op_of = lambda v, n: v.reshape(n, n).T                    # its inverse (Lemma D.1(b))
cplx = lambda *s: rng.standard_normal(s) + 1j * rng.standard_normal(s)

n = 3
M3 = TypeIAlgebra([n], [n], precision='double')
H = cplx(n, n)
H = (H + H.conj().T) / 2
Ls = [cplx(n, n) / 2 for _ in range(2)]                   # arbitrary jumps, not traceless
to_op = lambda X: M3.operator(torch.from_numpy(X).reshape(1, 1, n, n))
S = dynamics.lindblad_superoperator(M3, to_op(H), [to_op(L) for L in Ls])[0][0].numpy()
J = choi_of(S, n)
q = vec_of(np.eye(n)) / math.sqrt(n)                      # the unit vector Omega / sqrt(n)
Pperp = np.eye(n * n) - np.outer(q, q)
assert np.abs(J - J.conj().T).max() < 1e-12               # Hermiticity preserving (Lemma D.1(d))
print(f"min eigenvalue of J(L): {np.linalg.eigvalsh(J).min():.2f}; "
      f"on the complement of Omega: {np.linalg.eigvalsh(Pperp @ J @ Pperp).min():.1e}")
assert np.linalg.eigvalsh(J).min() < -0.1                 # L itself is not CP ...
assert np.linalg.eigvalsh(Pperp @ J @ Pperp).min() > -1e-12        # ... but CCP (Theorem D.3)
for t in (0.01, 0.1, 1.0):
    assert np.linalg.eigvalsh(choi_of(sla.expm(t * S), n)).min() > -1e-12

L_theta = superop(lambda X: X.T, n) - np.eye(n * n)     # L = transpose - id: Hermiticity preserving, not CCP
anti = vec_of(np.eye(n)[:, [1]] @ np.eye(n)[[0], :] - np.eye(n)[:, [0]] @ np.eye(n)[[1], :])   # e1(x)e2 - e2(x)e1
assert abs(anti @ q) < 1e-15 and abs(anti @ choi_of(L_theta, n) @ anti + 2) < 1e-12
X = cplx(n, n)
for t in (0.01, 0.1, 1.0):
    T_t = sla.expm(t * L_theta)
    b = (1 - math.exp(-2 * t)) / 2
    assert abs(anti @ choi_of(T_t, n) @ anti + (1 - math.exp(-2 * t))) < 1e-12        # <v, J(e^{tL}) v>
    assert np.abs(T_t - ((1 - b) * np.eye(n * n) + b * superop(lambda Y: Y.T, n))).max() < 1e-12
    assert abs(np.trace((T_t @ X.reshape(-1)).reshape(n, n)) - np.trace(X)) < 1e-12   # trace preserving
    assert np.linalg.eigvalsh(choi_of(T_t, n)).min() < -0.9 * b                     # but not CP
```

<a id="s-D-3"></a>

## D.3 GKSL 定理

*讲义：* [D.3 节](../../lecture_notes/zh/D_gksl.md#s-D-3)；[定理 D.4](../../lecture_notes/zh/D_gksl.md#r-D-4)、[定理 D.5](../../lecture_notes/zh/D_gksl.md#r-D-5)、[定理 D.6](../../lecture_notes/zh/D_gksl.md#r-D-6)。

定理 D.4 的证明就是一个算法：由 $`J(\mathcal L)`$ 在 $`\Omega^\perp`$ 上的压缩得到 $`\Phi`$，再由 $`\xi_1`$ 得到 $`G`$，
然后得到 $`H`$。我们对上面的生成元运行这个算法，用本库重建生成元，并与按定理 D.5(iii) 和 (i) 计算出的
原始数据的归一化形式进行比较：

```python
def gksl_from_generator(S, n, tol=1e-10):
    """(H, [L_j]) with Tr H = Tr L_j = 0 from the matrix of a GKSL generator: the proof of Theorem D.4."""
    J = choi_of(S, n)
    q = vec_of(np.eye(n)) / math.sqrt(n)
    Pperp = np.eye(n * n) - np.outer(q, q)
    w, U = np.linalg.eigh(Pperp @ J @ Pperp)              # = J(Phi) >= 0
    jumps = [math.sqrt(x) * op_of(U[:, k], n) for k, x in enumerate(w) if x > tol]
    xi = J @ q - 0.5 * (q.conj() @ J @ q) * q
    G = op_of(xi, n) / math.sqrt(n)
    return 0.5j * (G - G.conj().T), jumps

H_rec, L_rec = gksl_from_generator(S, n)
S_rec = dynamics.lindblad_superoperator(M3, to_op(H_rec), [to_op(L) for L in L_rec])[0][0].numpy()
assert np.abs(S_rec - S).max() < 1e-12                    # the same generator
assert abs(np.trace(H_rec)) < 1e-12 and all(abs(np.trace(L)) < 1e-12 for L in L_rec)

a = [-np.trace(L) / n for L in Ls]                        # Theorem D.5(iii): make the jumps traceless
L_norm = [L + aj * np.eye(n) for L, aj in zip(Ls, a)]
H_norm = H - 0.5j * sum(np.conj(aj) * L - aj * L.conj().T for L, aj in zip(Ls, a))
H_norm = H_norm - np.trace(H_norm) / n * np.eye(n)        # Theorem D.5(i)
assert np.abs(H_norm - H_norm.conj().T).max() < 1e-12
assert np.abs(H_rec - H_norm).max() < 1e-12               # H is unique (Theorem D.5(b))
jump_choi = lambda Ks: sum(np.outer(vec_of(K), vec_of(K).conj()) for K in Ks)
assert np.abs(jump_choi(L_rec) - jump_choi(L_norm)).max() < 1e-12      # so is sum_j L_j . L_j^*
Wr, Wn = np.stack([vec_of(K) for K in L_rec], 1), np.stack([vec_of(K) for K in L_norm], 1)
u = np.linalg.lstsq(Wr, Wn, rcond=None)[0]                # L_norm_j = sum_k u_kj L_rec_k
assert np.abs(Wr @ u - Wn).max() < 1e-12 and np.abs(u.conj().T @ u - np.eye(2)).max() < 1e-12   # unitary
print("eigenvalues of the normalised H:", np.round(np.linalg.eigvalsh(H_rec), 4), " jumps:", len(L_rec))
```

由[定理 D.6](../../lecture_notes/zh/D_gksl.md#r-D-6)，当与 $`H`$ 以及所有 $`L_j`$、$`L_j^*`$ 可交换的只有 $`1`$ 的倍数时，忠实稳态 $`\sigma`$ 是唯一的；
此时还有 $`\ker\mathcal L^*=\mathbb C1`$。对上面的随机生成元，代码检验了假设与结论（$`\mathcal L^*`$ 在 $`A^T`$ 的按行展开
向量上的矩阵是 $`S^T`$: $`\mathrm{vec}(\mathcal L^*(A)^T)=S^T\mathrm{vec}(A^T)`$），并检验了
[讲义，D.3 节](../../lecture_notes/zh/D_gksl.md#s-D-3)的第二个例子：在 $`M_3`$ 上，$`H=0`$，$`L_1=|e_1\rangle\langle e_3|`$，$`L_2=|e_2\rangle\langle e_3|`$，
此时与 $`L_j,L_j^*`$ 可交换的只有 $`1`$ 的倍数，但 $`\dim\ker\mathcal L=4`$：

```python
def kernel_dim(S, tol=1e-9):
    return int((np.linalg.svd(S, compute_uv=False) < tol).sum())

def commutant_dim(ops, n):                    # dim {X : XA = AX for all A in ops}
    return kernel_dim(np.concatenate([np.kron(A, np.eye(n)) - np.kron(np.eye(n), A.T) for A in ops]))

w, V = np.linalg.eig(S)
sigma = V[:, np.argmin(abs(w))].reshape(n, n)             # row-major: vec(X)_{(a,b)} = X_ab
sigma = sigma / np.trace(sigma)
assert np.abs(S @ sigma.reshape(-1)).max() < 1e-12        # stationary
assert np.abs(sigma - sigma.conj().T).max() < 1e-10 and np.linalg.eigvalsh(sigma).min() > 1e-3   # faithful
assert commutant_dim([H] + Ls + [L.conj().T for L in Ls], n) == 1
assert kernel_dim(S) == 1 and kernel_dim(S.T) == 1        # Theorem D.6
assert np.abs(S.T @ np.eye(n).reshape(-1)).max() < 1e-12  # ker L^* = C 1

e3 = np.eye(3)
L_abs = [np.outer(e3[0], e3[2]), np.outer(e3[1], e3[2])]       # L_1 = |e_1><e_3|, L_2 = |e_2><e_3|
S_abs = superop(lambda R: sum(L @ R @ L.T - 0.5 * (L.T @ L @ R + R @ L.T @ L) for L in L_abs), 3)
assert commutant_dim(L_abs + [L.T for L in L_abs], 3) == 1
print("second example: dim ker L =", kernel_dim(S_abs))
assert kernel_dim(S_abs) == 4                             # no faithful stationary state, no uniqueness
```

<a id="s-D-4"></a>

## D.4 有扇区的代数

*讲义：* [D.4 节](../../lecture_notes/zh/D_gksl.md#s-D-4)；[定理 D.7](../../lecture_notes/zh/D_gksl.md#r-D-7)、[推论 D.8](../../lecture_notes/zh/D_gksl.md#r-D-8)。

我们在 $`M_1\oplus M_2\subseteq M_3`$ 上，对 `lindblad_rhs` 的一个生成元检验[定理 D.7](../../lecture_notes/zh/D_gksl.md#r-D-7)，它有一个扇区内部的
跳跃和一个随机的扇区之间的跳跃，检验按照证明的五个步骤进行；我们还检验[推论 D.8](../../lecture_notes/zh/D_gksl.md#r-D-8) 的对偶生成元。正如证明
之后所指出的，这里 $`J(\mathcal L\circ E)`$ 在 $`\Omega^\perp`$ 上不是正的：

```python
alg = TypeIAlgebra([1, 2], [1, 2], precision='double')    # M = M_1 (+) M_2 in M_3
D, ks = 3, alg.k_factors
P = [np.diag([1.0, 0, 0]), np.diag([0, 1.0, 1.0])]
E = lambda X: sum(Pc @ X @ Pc for Pc in P)               # the conditional expectation onto M
def blocks(X):                                            # dense block-diagonal 3 x 3 -> library layout
    out = torch.zeros(1, 2, 2, 2, dtype=c128)
    out[0, 0, :1, :1] = torch.from_numpy(X[:1, :1])
    out[0, 1] = torch.from_numpy(X[1:, 1:])
    return out
dense = lambda Y: sla.block_diag(Y[0, 0, :1, :1].numpy(), Y[0, 1].numpy())
Hs = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
inside = alg.operator(torch.randn(1, 2, 2, 2, dtype=c128))       # a jump inside the sectors
between = channels.random_inter_sector_channel(alg, alg, kraus_rank=2)
rhs = dynamics.lindblad_rhs(Hs, [inside, between], rates=[0.7, 1.3], algebra=alg)
L_M = lambda X: dense(rhs(0.0, blocks(X)))               # the generator on M (X block diagonal)
J = choi_of(superop(lambda X: L_M(E(X)), D), D)          # J(L o E) on M_3
p = [vec_of(Pc) for Pc in P]
Q = np.stack([pc / np.linalg.norm(pc) for pc in p], 1)    # orthonormal q_c = p_c / sqrt(k_c)
Pw = np.eye(D * D) - Q @ Q.T                              # projection onto the complement of W
assert np.linalg.eigvalsh(Pw @ J @ Pw).min() > -1e-12     # Step 2
w, U = np.linalg.eigh(Pw @ J @ Pw)                        # Step 3 (Lemma D.2)
jumps = [math.sqrt(x) * op_of(U[:, k], D) for k, x in enumerate(w) if x > 1e-10]
xis = [J @ Q[:, c] - 0.5 * sum((Q[:, d] @ J @ Q[:, c]) * Q[:, d] for d in range(2)) for c in range(2)]
G = E(sum(op_of(xi, D) @ Pc / math.sqrt(k) for xi, Pc, k in zip(xis, P, ks)))   # Step 4
Phi = lambda X: E(sum(K @ X @ K.conj().T for K in jumps))
Phi_dual_one = E(sum(K.conj().T @ K for K in jumps))      # Phi^*(1)
H_M = 0.5j * (G - G.conj().T)
assert np.abs(G + G.conj().T + Phi_dual_one).max() < 1e-12      # Step 5
for _ in range(3):
    X = E(cplx(D, D))
    gksl = Phi(X) - 0.5 * (Phi_dual_one @ X + X @ Phi_dual_one) - 1j * (H_M @ X - X @ H_M)
    assert np.abs(L_M(X) - gksl).max() < 1e-12            # Theorem D.7
    dual = lambda A: E(sum(K.conj().T @ A @ K for K in jumps)) + G.conj().T @ A + A @ G   # Corollary D.8
    assert abs(np.trace(L_M(X) @ (Y := E(cplx(D, D)))) - np.trace(X @ dual(Y))) < 1e-12
Pomega = np.eye(D * D) - np.outer(sum(p), sum(p)) / D     # projection onto the complement of Omega
print(f"numerical observation: min eigenvalue of J(L o E) on the complement of Omega: "
      f"{np.linalg.eigvalsh(Pomega @ J @ Pomega).min():.2f}")
assert np.linalg.eigvalsh(Pomega @ J @ Pomega).min() < -0.1
```

<a id="s-D-5"></a>

## D.5 算子代数纠错

*讲义：* [D.5 节](../../lecture_notes/zh/D_gksl.md#s-D-5)；[定理 D.9](../../lecture_notes/zh/D_gksl.md#r-D-9)、[推论 D.10](../../lecture_notes/zh/D_gksl.md#r-D-10)。

具有 $`C\ge2`$ 个扇区且某个 $`k_c\ge2`$ 的可纠正 $`\mathfrak A`$ 是一个**混合**码（[讲义，D.5 节](../../lecture_notes/zh/D_gksl.md#s-D-5)）。
作为例子，我们在 $`\mathbb C^8`$ 中一个随机的五维码上取
$`\mathfrak A=(M_2\otimes1_2)\oplus M_1`$，噪声通过 $`\sum_iB_{c,i}^*B_{c,i}=1_{m_c}`$ 作用在规范因子 $`\mathbb C^{m_c}`$ 上，
把各扇区映到相互正交的子空间，打乱 $`\mathcal C^\perp`$，然后再被一个随机酉矩阵混合（定理 6.7）。
Knill–Laflamme 条件不成立，[定理 D.9](../../lecture_notes/zh/D_gksl.md#r-D-9) 的条件成立，而证明中的恢复还原了每个 $`A\in\mathfrak A`$，
却不能还原 $`\mathcal C`$ 上的一般算子：

```python
D, dims = 8, [(2, 2), (1, 1)]                             # A = (M_2 (x) 1_2) (+) (M_1 (x) 1_1), dim C = 5
U8 = np.linalg.qr(cplx(D, D))[0]
Vc = np.linalg.qr(cplx(D, 5))[0]                          # orthonormal basis |c; a, b> of the code
W = [Vc[:, :4], Vc[:, 4:]]                                # sector c: columns indexed by (a, b), row-major
Pcode = sum(Wc @ Wc.conj().T for Wc in W)
Y = [U8[:, :4], U8[:, 4:6]]                               # orthogonal images C^{k_c} (x) C^{s_c}, s_c = 2
r = 3
B = [np.linalg.qr(cplx(2 * r, m))[0].reshape(r, 2, m) for _, m in dims]   # sum_i B_i^* B_i = 1_{m_c}
E_code = [sum(Yc @ np.kron(np.eye(k), Bc[i]) @ Wc.conj().T for Yc, Bc, Wc, (k, m) in zip(Y, B, W, dims))
          for i in range(r)]
kraus = E_code + [np.linalg.qr(cplx(D, D))[0] @ (np.eye(D) - Pcode)]
mix = np.linalg.qr(cplx(r + 1, r + 1))[0]                 # Theorem 6.7: an equivalent family
kraus = [sum(mix[i, l] * K for l, K in enumerate(kraus)) for i in range(r + 1)]
assert np.abs(sum(K.conj().T @ K for K in kraus) - np.eye(D)).max() < 1e-12      # a channel

units = []                                                 # matrix units of A
for Wc, (k, m) in zip(W, dims):
    for a1 in range(k):
        for a2 in range(k):
            units.append(Wc @ np.kron(np.eye(k)[:, [a1]] @ np.eye(k)[[a2], :], np.eye(m)) @ Wc.conj().T)
comm = lambda X: max(np.abs(X @ A - A @ X).max() for A in units)
assert max(comm(Pcode @ Ki.conj().T @ Kj @ Pcode) for Ki in kraus for Kj in kraus) < 1e-12   # Theorem D.9
KL_dev = max(np.abs(Pcode @ Ki.conj().T @ Kj @ Pcode - np.trace(Pcode @ Ki.conj().T @ Kj) / 5 * Pcode).max()
             for Ki in kraus for Kj in kraus)
print(f"Knill-Laflamme condition violated by {KL_dev:.2f}")
assert KL_dev > 0.05

R_kraus, Qsum = [], np.zeros((D, D), complex)             # Steps 1-3 of the proof
for Wc, (k, m) in zip(W, dims):
    F = [Ki @ Wc[:, b::m] for Ki in kraus for b in range(m)]   # F_(i,b) e_a = E_i |c; a, b>
    Lam = np.array([[(Fa.conj().T @ Fb)[0, 0] for Fb in F] for Fa in F])
    assert max(np.abs(Fa.conj().T @ Fb - Lam[x, y] * np.eye(k)).max()
               for x, Fa in enumerate(F) for y, Fb in enumerate(F)) < 1e-12      # Step 1
    d, u = np.linalg.eigh(Lam)
    iota = Wc[:, 0::m]                                    # e_a -> |c; a, 1>
    for kk in range(len(d)):
        if d[kk] > 1e-12:
            Vk = sum(u[x, kk] * Fx for x, Fx in enumerate(F)) / math.sqrt(d[kk])
            R_kraus.append(iota @ Vk.conj().T)
            Qsum += Vk @ Vk.conj().T
R_kraus.append(np.eye(D) - Qsum)
assert np.abs(sum(R.conj().T @ R for R in R_kraus) - np.eye(D)).max() < 1e-10    # R is a channel
dual = lambda X: sum(Ki.conj().T @ R.conj().T @ X @ R @ Ki for R in R_kraus for Ki in kraus)
for _ in range(3):
    A = sum(complex(z) * Ue for z, Ue in zip(cplx(len(units)), units))
    assert np.abs(Pcode @ dual(A) @ Pcode - A).max() < 1e-10                       # A is corrected
B_gen = Pcode @ cplx(D, D) @ Pcode
print(f"a generic operator on the code: error {np.abs(Pcode @ dual(B_gen) @ Pcode - B_gen).max():.2f}")
assert np.abs(Pcode @ dual(B_gen) @ Pcode - B_gen).max() > 0.1
assert max(comm(R @ Ki @ Pcode) for R in R_kraus for Ki in kraus) < 1e-10   # proof of necessity: W_mi in A'
```

**下一章：** [教程，附录 E. 欧几里得若尔当代数的分类](E_jordan_classification.md)。
